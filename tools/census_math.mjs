// census_math.mjs
//
// The small amount of mathematics shared by the data build (build_site_data.mjs) and the website,
// which gets a copy of this file: keys, singular fibres, and the factors of P2.
// Polynomials are arrays of coefficients, constant term first.  Every coefficient of L and P2
// is far below 2^53, so ordinary numbers are exact.

export function decodeKey(key) {
  // The model of an orbit key, as FibrationFromKeyF2 reads it.
  //   Input:  key, an integer in [0, 2^32)
  //   Output: { a1, a3, a4, a6 }, each a list of [s exponent, t exponent] of its monomials
  const form = (n, d) => {
    const out = [];
    for (let i = 0; i <= d; i++) if (Math.floor(n / 2 ** (d - i)) % 2) out.push([d - i, i]);
    return out;
  };
  return { a1: form(Math.floor(key / 2 ** 29), 2), a3: form(Math.floor(key / 2 ** 22) % 128, 6),
           a4: form(Math.floor(key / 2 ** 13) % 512, 8), a6: form(key % 8192, 12) };
}

export function parseFibres(text) {
  // The singular_fibres column as records.
  //   Input:  text, "deg:f:c:type;deg:f:c:type;..."
  //   Output: a list of { deg, f, c, type } with deg, f, c integers
  return text.split(';').map(r => { const [deg, f, c, type] = r.split(':');
    return { deg: Number(deg), f: Number(f), c: Number(c), type }; });
}

export function componentCycles(type, c) {
  // Frobenius orbit sizes on the non-identity components of a fibre; ComponentCycles in
  // functions_census.m.
  //   Input:  type, the Kodaira symbol ("I4", "IV*", "I0*"); c, the Tamagawa number
  //   Output: a list of positive integers adding up to (number of components - 1)
  const ones = (k) => Array(k).fill(1);
  if (type === 'II') return [];
  if (type === 'III') return [1];
  if (type === 'IV') return c === 3 ? [1, 1] : [2];
  if (type === 'IV*') return c === 3 ? ones(6) : [1, 1, 2, 2];
  if (type === 'III*') return ones(7);
  if (type === 'II*') return ones(8);
  if (type.endsWith('*')) {
    const n = Number(type.slice(1, -1));
    if (n === 0) {
      if (c === 4) return ones(4);
      if (c === 2) return [1, 1, 2];
      if (c === 1) return [1, 3];
    } else {
      if (c === 4) return ones(n + 4);
      if (c === 2) return [...ones(n + 2), 2];
    }
  } else {
    const n = Number(type.slice(1));
    if (n <= 1) return [];
    if (c === n) return ones(n - 1);
    if (n % 2 === 1 && c === 1) return Array((n - 1) / 2).fill(2);
    if (n % 2 === 0 && c === 2) return [1, ...Array((n - 2) / 2).fill(2)];
  }
  throw new Error(`unexpected Tamagawa number ${c} for Kodaira type ${type}`);
}

export function typeOrder(type) {
  // A sort position for Kodaira types: I_n, then I_n*, then II, III, IV, IV*, III*, II*.
  //   Input:  type, the Kodaira symbol
  //   Output: a number
  const named = { 'II': 2000, 'III': 2001, 'IV': 2002, 'IV*': 2003, 'III*': 2004, 'II*': 2005 };
  if (type in named) return named[type];
  return type.endsWith('*') ? 1000 + Number(type.slice(1, -1)) : Number(type.slice(1));
}

export function configuration(records) {
  // The fibre configuration over the algebraic closure, in a canonical form.
  //   Input:  records, from parseFibres
  //   Output: a string such as "16 I1, I4, IV"
  const n = new Map();
  for (const r of records) n.set(r.type, (n.get(r.type) || 0) + r.deg);
  return [...n].sort((a, b) => typeOrder(a[0]) - typeOrder(b[0]))
    .map(([t, k]) => (k > 1 ? `${k} ${t}` : t)).join(', ');
}

export function normaliseConfiguration(text) {
  // A configuration typed by a person, in the canonical form.
  //   Input:  text, such as "iv, 16 I1,I4" or "16I1 I4 IV"
  //   Output: the canonical string, or null if it isn't a valid configuration
  const parts = text.toUpperCase().match(/(\d+\s*)?(IV\*|III\*|II\*|IV|III|II|I\d+\*?)/g);
  if (!parts || parts.join('').replace(/\s/g, '') !== text.toUpperCase().replace(/[\s,;]/g, '')) return null;
  const n = new Map();
  for (const p of parts) {
    const m = p.match(/^(\d+)?\s*(.+)$/);
    n.set(m[2], (n.get(m[2]) || 0) + (m[1] ? Number(m[1]) : 1));
  }
  return [...n].sort((a, b) => typeOrder(a[0]) - typeOrder(b[0]))
    .map(([t, k]) => (k > 1 ? `${k} ${t}` : t)).join(', ');
}

export const isReducible = (type) => type !== 'I1' && type !== 'II';

export function polyMul(a, b) {
  // The product of two polynomials.
  //   Input:  a, b, coefficient lists, constant term first
  //   Output: their product, as a coefficient list
  const out = new Array(a.length + b.length - 1).fill(0);
  a.forEach((x, i) => { if (x) b.forEach((y, j) => { out[i + j] += x * y; }); });
  return out;
}

export function fibreFactors(records) {
  // The factors 1 - (2u)^(k deg P) of P2 from the fibres, one per Frobenius orbit of components.
  //   Input:  records, from parseFibres
  //   Output: a list of exponents e = k * deg P, each meaning the factor 1 - (2u)^e
  const out = [];
  for (const r of records) for (const k of componentCycles(r.type, r.c)) out.push(k * r.deg);
  return out.sort((a, b) => a - b);
}

export function p2FromData(L, records) {
  // P2 from the L-function and the fibres; P2FromData in functions_census.m, with q = 2.
  //   Input:  L, a coefficient list; records, from parseFibres
  //   Output: P2, a coefficient list
  let p = polyMul([1, -4, 4], L);
  for (const e of fibreFactors(records)) {
    const f = new Array(e + 1).fill(0); f[0] = 1; f[e] = -(2 ** e);
    p = polyMul(p, f);
  }
  return p;
}

function divideExact(p, d) {
  // p / d when d divides p in Z[u] and d has constant term 1; otherwise null.
  if (p.length < d.length) return null;
  const q = new Array(p.length - d.length + 1).fill(0), r = p.slice();
  for (let i = 0; i < q.length; i++) {
    q[i] = r[i];
    if (q[i]) for (let j = 0; j < d.length; j++) r[i + j] -= q[i] * d[j];
  }
  return r.every(x => x === 0) ? q : null;
}

const cyclotomic = (() => {
  // The cyclotomic polynomials of degree at most 22, with constant term 1 (so 1 - u for n = 1),
  // in increasing order of n.  The largest such n is 66.
  const phi = [null], out = [];
  const N = 200;
  for (let n = 1; n <= N; n++) {
    let p = new Array(n + 1).fill(0); p[0] = 1; p[n] = -1;     // 1 - u^n = prod_{d|n} rev(Phi_d)
    for (let d = 1; d < n; d++) if (n % d === 0) p = divideExact(p, phi[d]);
    phi[n] = p;
    if (p.length - 1 <= 22) out.push({ n, poly: p });
  }
  return out;
})();
// psi[n](u) = prod over primitive n-th roots z of (1 - 2 z u)
const cyclotomicScaled = cyclotomic.map(({ n, poly }) => ({ n, poly: poly.map((c, i) => c * 2 ** i) }));

export function splitP2(p) {
  // Splits a polynomial into the factors whose inverse roots are 2 * (root of unity) and the rest.
  //   Input:  p, a coefficient list with constant term 1 (P2 or L)
  //   Output: { algebraic: [{ n, mult, poly }], rest }, where n is the order of the root of unity
  //           and rest is the remaining factor, not necessarily irreducible
  let rest = p.slice();
  const algebraic = [];
  for (const { n, poly } of cyclotomicScaled) {
    let mult = 0;
    for (;;) {
      const q = divideExact(rest, poly);
      if (!q) break;
      rest = q; mult++;
    }
    if (mult) algebraic.push({ n, mult, poly });
  }
  return { algebraic, rest };
}

const gcd = (a, b) => (b ? gcd(b, a % b) : a);

// The n with Phi_n of degree at most 22, and the degree of Phi_n.
export const CYCLOTOMIC_ORDERS = cyclotomic.map(c => c.n);
export function cyclotomicDegree(n) {
  const c = cyclotomic.find(x => x.n === n);
  if (!c) throw new Error(`Phi_${n} has degree above 22`);
  return c.poly.length - 1;
}
export const cyclotomicCoefficients = (n) => cyclotomic.find(x => x.n === n).poly.slice();

export function cyclotomicPart(p2) {
  // The cyclotomic part of P2: P2(u) = prod_n Phi_n(2u)^e_n * Q(u), where Phi_n is the n-th
  // cyclotomic polynomial with constant term 1 (so Phi_1(x) = 1 - x) and Q has no inverse root
  // of the form 2 * (root of unity).
  //   Input:  p2, a coefficient list
  //   Output: { factors: [{ n, e }], text, nsDegree, rest }, where text is canonical, such as
  //           "1^4 3^2 6"; nsDegree is the lcm of the n, the degree of the smallest extension of
  //           F_2 over which the Neron-Severi group of the surface over the closure is defined
  const { algebraic, rest } = splitP2(p2);
  const factors = algebraic.map(({ n, mult }) => ({ n, e: mult }));
  return {
    factors, rest,
    text: factors.map(({ n, e }) => (e > 1 ? `${n}^${e}` : `${n}`)).join(' '),
    nsDegree: factors.reduce((m, { n }) => m / gcd(m, n) * n, 1),
  };
}

export function parseCyclotomic(text) {
  // A cyclotomic part typed by a person, in the canonical form of cyclotomicPart.
  //   Input:  text, such as "1^4 3^2 6", "Phi1^4 Phi3^2 Phi6" or "1^4, 6, 3^2"
  //   Output: the canonical string, or null if it isn't one
  const clean = text.replace(/phi_?|Φ_?|\{|\}|_/gi, ' ').replace(/[,*·]/g, ' ').trim();
  if (!clean || !/^(\d+\s*(\^\s*\d+)?\s*)+$/.test(clean)) return null;
  const e = new Map();
  for (const m of clean.matchAll(/(\d+)\s*(?:\^\s*(\d+))?/g)) {
    const n = Number(m[1]), k = m[2] ? Number(m[2]) : 1;
    if (!cyclotomic.some(c => c.n === n) || k < 1) return null;
    e.set(n, (e.get(n) || 0) + k);
  }
  return [...e].sort((a, b) => a[0] - b[0]).map(([n, k]) => (k > 1 ? `${n}^${k}` : `${n}`)).join(' ');
}

// ---------------------------------------------------------------------------------------
// LMFDB labels of zeta functions
// ---------------------------------------------------------------------------------------
// The LMFDB's table of motives of K3 type over finite fields (github.com/roed-math/k3-fq) stores
// P(T) = 2 det(1 - Frob T / 2 | H^2) = 2 P2(T/2), whose roots lie on the unit circle.  Its label
// is 2.rho.c.t: rho is the geometric Picard number, c numbers the cyclotomic part of P among
// those of degree rho, and t encodes the lower half of the coefficients of the remaining factor.
// The functions below follow k3fq/labels.py in that repository.

export function lmfdbPoly(p2) {
  // P(T) = 2 P2(T/2).
  //   Input:  p2, the coefficient list of P2, of length 23
  //   Output: the coefficient list of P, with constant term 2 and last coefficient 2 or -2
  if (p2.length !== 23) throw new Error('P2 must have degree 22');
  const P = p2.map((c, i) => c / 2 ** (i - 1));
  if (!P.every(Number.isInteger)) throw new Error('2 P2(T/2) is not integral');
  return P;
}

const tailCounts = (() => {
  // N[i][k] = the number of products of cyclotomic polynomials of degree k that only use
  // cyclotomic[i], cyclotomic[i+1], ...
  const n = cyclotomic.length;
  const N = Array.from({ length: n + 1 }, () => new Array(23).fill(0));
  N[n][0] = 1;
  for (let i = n - 1; i >= 0; i--) {
    const ph = cyclotomic[i].poly.length - 1;
    for (let k = 0; k <= 22; k++) for (let e = 0; e * ph <= k; e++) N[i][k] += N[i + 1][k - e * ph];
  }
  return N;
})();

function letterCode(n) {
  // A non-negative integer in base 26 with the digits a to z (a = 0).
  let out = '';
  do { out = String.fromCharCode(97 + n % 26) + out; n = Math.floor(n / 26); } while (n);
  return out;
}

function balancedDigits(n) {
  // An integer in balanced base 26, digits from -13 to 12, most significant first.
  const digits = [];
  do {
    const r = ((n + 13) % 26 + 26) % 26 - 13;
    digits.unshift(r);
    n = (n - r) / 26;
  } while (n);
  return digits;
}

export function zetaLabel(p2) {
  // The LMFDB label of a zeta function.
  //   Input:  p2, the coefficient list of P2, of length 23
  //   Output: a string such as "2.20.a.o"
  let rest = lmfdbPoly(p2);
  const exps = cyclotomic.map(({ poly }) => {
    let e = 0;
    for (;;) {
      const q = divideExact(rest, poly);
      if (!q) break;
      rest = q; e++;
    }
    return e;
  });
  const rho = 23 - rest.length;
  // The position of the exponents among those of degree rho, largest first in the order of n.
  let index = 0, k = rho;
  cyclotomic.forEach(({ poly }, i) => {
    const ph = poly.length - 1;
    for (let f = exps[i] + 1; f * ph <= k; f++) index += tailCounts[i + 1][k - f * ph];
    k -= exps[i] * ph;
  });
  const label = `2.${rho}.${letterCode(index)}`;
  const code = transcendentalCode(rest);
  return code ? label + '.' + code : label;
}

function transcendentalCode(rest) {
  // The last part of an LMFDB label: the coefficients of T, ..., T^(d/2) of the transcendental
  // part 2 Q(T/2), each in balanced base 26 and padded to a common width.
  //   Input:  rest, the coefficient list of 2 Q(T/2), palindromic with constant term 2
  //   Output: a string, '' when the degree is 0
  const d = rest.length - 1;
  if (d === 0) return '';
  if (d % 2 || rest[0] !== 2 || rest.some((c, i) => c !== rest[d - i])) throw new Error('unexpected transcendental factor');
  const lower = rest.slice(1, d / 2 + 1).map(balancedDigits);
  const width = Math.max(...lower.map(x => x.length));
  return lower.map(x => 'n'.repeat(width - x.length) + x.map(v => String.fromCharCode(110 + v)).join('')).join('');
}

// ---------------------------------------------------------------------------------------
// Transcendental parts
// ---------------------------------------------------------------------------------------
// The transcendental part of P2 is the factor Q with no inverse root of the form 2 * (root of
// unity).  It is named "d.code": its degree d = 22 - rho(X over the closure), and the last part
// of the LMFDB label of any zeta function with this Q.  When d = 0 the name is "0".

export function transcendentalId(p2) {
  // The name of the transcendental part of P2.  Input: P2's coefficients.  Output: a string.
  const { rest } = splitP2(p2);
  const code = transcendentalCode(rest.map((c, i) => c * 2 / 2 ** i));
  return code ? `${rest.length - 1}.${code}` : '0';
}

export function transcendentalPoly(id) {
  // The transcendental part in the LMFDB's normalisation, 2 Q(T/2), from its name.
  //   Input:  id, such as "20.nlnonmnonm"
  //   Output: the coefficient list, or null if id isn't a valid name
  if (id === '0') return [2];
  const m = /^(\d+)\.([a-z]+)$/.exec(id);
  if (!m) return null;
  const d = Number(m[1]), k = d / 2;
  if (!d || d % 2 || d > 22 || m[2].length % k) return null;
  const w = m[2].length / k;
  const chunks = Array.from({ length: k }, (_, i) => m[2].slice(i * w, (i + 1) * w));
  if (w > 1 && chunks.every(c => c[0] === 'n')) return null;     // padded more than necessary
  const lower = chunks.map(c => [...c].reduce((n, ch) => 26 * n + ch.charCodeAt(0) - 110, 0));
  return [2, ...lower, ...lower.slice(0, -1).reverse(), 2];
}

export function parseTranscendental(text) {
  // A transcendental part typed by a person, as its name.
  //   Input:  text: a name ("20.nlnonmnonm"), an LMFDB label ("2.2.a.nlnonmnonm"), or the
  //           coefficients, constant term first, either of Q(u) (constant term 1) or of 2 Q(T/2)
  //           (constant term 2)
  //   Output: the name, or null if text is none of these
  const s = text.trim();
  const label = /^2\.(\d+)\.[a-z]+(?:\.([a-z]+))?$/.exec(s);
  if (label) {
    const id = label[2] ? `${22 - Number(label[1])}.${label[2]}` : '0';
    return transcendentalPoly(id) ? id : null;
  }
  if (/^\d+\.[a-z]+$/.test(s) || s === '0') return transcendentalPoly(s) ? s : null;
  const nums = s.match(/-?\d+/g);
  if (!nums || nums.join(' ') !== s.replace(/[\[\],]/g, ' ').trim().split(/\s+/).join(' ')) return null;
  let c = nums.map(Number);
  if (c[0] === 1) c = c.map((x, i) => x * 2 / 2 ** i);
  if (c[0] !== 2 || !c.every(Number.isInteger)) return null;
  try {
    const code = transcendentalCode(c);
    return code ? `${c.length - 1}.${code}` : '0';
  } catch { return null; }
}

export function geometricRank(rhoGeom, records) {
  // The rank of E over the function field of P^1 over the algebraic closure, by Shioda-Tate:
  // rho(X over the closure) = 2 + sum over places of deg P (components - 1) + rank.
  //   Input:  rhoGeom, the geometric Picard number; records, from parseFibres
  //   Output: the rank, an integer
  return rhoGeom - 2 - records.reduce((s, r) => s + r.deg * componentCycles(r.type, r.c).reduce((a, b) => a + b, 0), 0);
}

export function rowFromCsv(cols) {
  // One surface as an object with the column names of the site's search table.
  //   Input:  cols, the fields of a row of the site's row files
  //   Output: an object; height is null for a supersingular surface
  const [key, a1, torsion, fibres, L, L_sign, rank, zeta_sign, rho_f2, rho_geom, height, br_disc,
         pts_f2, pts_f4, pts_f8, pts_f16, pts_f32] = cols;
  const records = parseFibres(fibres);
  const lc = L.split(' ').map(Number);
  return {
    key: Number(key), a1, torsion, fibres, L: lc, records,
    config: configuration(records),
    n_places: records.length,
    n_singular: records.reduce((s, r) => s + r.deg, 0),
    n_reducible: records.reduce((s, r) => s + (isReducible(r.type) ? r.deg : 0), 0),
    L_degree: lc.length - 1,
    L_sign: Number(L_sign), analytic_rank: Number(rank), zeta_sign: Number(zeta_sign),
    rho_f2: Number(rho_f2), rho_geom: Number(rho_geom),
    rank_geom: geometricRank(Number(rho_geom), records),
    height: height === 'inf' ? null : Number(height), br_disc: Number(br_disc),
    pts_f2: Number(pts_f2), pts_f4: Number(pts_f4), pts_f8: Number(pts_f8),
    pts_f16: Number(pts_f16), pts_f32: Number(pts_f32),
  };
}
