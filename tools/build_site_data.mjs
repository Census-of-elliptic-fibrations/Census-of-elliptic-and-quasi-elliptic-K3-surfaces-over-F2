// build_site_data.mjs
//
// Turns the census CSV into the data files of the census website.  Run from the repository root:
//
//     cd tools && npm install && cd ..
//     node tools/build_site_data.mjs <census.csv> <website folder>
//
// It writes, in <website folder>:
//     data/surfaces.parquet   one row per surface, the columns the search filters on
//     data/fibres.parquet     one row per bad place, for the per-fibre filters
//     data/zeta.parquet       one row per zeta function that occurs, with its LMFDB label
//     data/cyclotomic.parquet one row per cyclotomic factor of each zeta function
//     data/rows/NNNNN.csv     1000 surfaces each, with L; what a surface page loads
//     data/rows/index.json    the first key of each of those files
//     data/missing_keys.txt   census models with no row yet
//     data/stats.json         the fixed charts of the statistics tab
//     data/manifest.json      counts and where the data came from
//     js/census_math.mjs      a copy of tools/census_math.mjs
//
// It also checks every row: the fibres add up to Euler number 24 by Ogg's formula, and P2 computed
// from L and the fibres equals the stored P2.  On a sample it checks the Picard numbers too.
// Surfaces with the same P2 must agree on everything that P2 determines.
// census.txt.gz is the list of all models, used to find what is missing.  Takes a few minutes.

import fs from 'node:fs';
import path from 'node:path';
import zlib from 'node:zlib';
import readline from 'node:readline';
import { execSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';
import { DuckDBInstance } from '@duckdb/node-api';
import { parseFibres, componentCycles, configuration, isReducible, p2FromData, splitP2, geometricRank, zetaLabel, cyclotomicPart, cyclotomicDegree, transcendentalId } from './census_math.mjs';

const [csvPath, siteDir] = process.argv.slice(2);
if (!csvPath || !siteDir) throw new Error('usage: node tools/build_site_data.mjs <census.csv> <website folder>');
const here = path.dirname(fileURLToPath(import.meta.url));
const repo = path.dirname(here);
const dataDir = path.join(siteDir, 'data');
const rowsDir = path.join(dataDir, 'rows');
const tmpDir = path.join(siteDir, '.build');
for (const d of [rowsDir, tmpDir, path.join(siteDir, 'js')]) fs.mkdirSync(d, { recursive: true });
for (const f of fs.readdirSync(rowsDir)) fs.rmSync(path.join(rowsDir, f));

// ---------------------------------------------------------------------------------------
// All models, from census.txt.gz
// ---------------------------------------------------------------------------------------

function keyOfModel(line) {
  // The orbit key of a line of census.txt.gz.
  //   Input:  line, "[a1, a2, a3, a4, a6]" with forms in s, t
  //   Output: the key, an integer
  const mask = (form, d) => {
    if (form === '0') return 0;
    let m = 0;
    for (const term of form.split(' + ')) {
      let e = 0;
      for (const f of term.split('*')) { const [v, x] = f.split('^'); if (v === 't') e = x ? Number(x) : 1; }
      m += 2 ** (d - e);
    }
    return m;
  };
  const [a1, , a3, a4, a6] = line.trim().slice(1, -1).split(', ');
  return mask(a1, 2) * 2 ** 29 + mask(a3, 6) * 2 ** 22 + mask(a4, 8) * 2 ** 13 + mask(a6, 12);
}
const A1_NAME = { 0: '0', 4: 's^2', 2: 's*t', 7: 's^2 + s*t + t^2' };
const allKeys = zlib.gunzipSync(fs.readFileSync(path.join(repo, 'census.txt.gz'))).toString().trim().split('\n').map(keyOfModel);
const expectedByA1 = {};
for (const k of allKeys) { const a = A1_NAME[Math.floor(k / 2 ** 29)]; expectedByA1[a] = (expectedByA1[a] || 0) + 1; }
const allSet = new Set(allKeys);
console.log('models in census.txt.gz:', allKeys.length);

// ---------------------------------------------------------------------------------------
// One pass over the CSV
// ---------------------------------------------------------------------------------------

const count = () => new Map();
const bump = (m, k, n = 1) => m.set(k, (m.get(k) || 0) + n);
const H = { rank: count(), rank_geom: count(), rho_f2: count(), rho_geom: count(), L_sign: count(), zeta_sign: count(),
            torsion: count(), br: count(), ns_degree: count(), config: count(), L_degree: count(), n_places: count(),
            singular: count(), reducible: count(), height: count(),
            pts_f2: count(), pts_f4: count(), pts_f8: count(), pts_f16: count(), pts_f32: count() };
const WIDTH = { pts_f2: 1, pts_f4: 2, pts_f8: 5, pts_f16: 10, pts_f32: 25 };
const byA1 = {};
const surfacesOut = fs.createWriteStream(path.join(tmpDir, 'surfaces.tsv'));
const fibresOut = fs.createWriteStream(path.join(tmpDir, 'fibres.tsv'));
const firstKeys = [];
let shard = [], shardNo = 0, n = 0, problems = 0, lastKey = -1, sampled = 0;
const have = new Set();
const zetas = new Map();      // P2 as text -> { id, label, n, first, facts }
const fail = (key, why) => { if (problems++ < 20) console.log('PROBLEM', key, why); };

const write = (s, text) => new Promise(res => { if (s.write(text)) res(); else s.once('drain', res); });
const flushShard = () => {
  if (!shard.length) return;
  firstKeys.push(Number(shard[0].slice(0, shard[0].indexOf(','))));
  fs.writeFileSync(path.join(rowsDir, String(shardNo++).padStart(5, '0') + '.csv'), shard.join('\n') + '\n');
  shard = [];
};

const rl = readline.createInterface({ input: fs.createReadStream(csvPath), crlfDelay: Infinity });
let header = null;
for await (const line of rl) {
  if (!header) {
    header = line.split(',');
    if (header.length !== 18 || header[0] !== 'key') throw new Error('unexpected header: ' + line);
    continue;
  }
  const c = line.split(',');
  if (c.length !== 18) { fail(c[0], 'has ' + c.length + ' columns'); continue; }
  const [key, a1, torsion, fibres, L, L_sign, rank, P2, zeta_sign, rho_f2, rho_geom, height, br_disc,
         pts_f2, pts_f4, pts_f8, pts_f16, pts_f32] = c;
  const k = Number(key);
  if (k <= lastKey) throw new Error('keys are not ascending at ' + key);
  lastKey = k;
  if (!allSet.has(k)) fail(key, 'is not in census.txt.gz');
  have.add(k);

  const records = parseFibres(fibres);
  const lc = L.split(' ').map(Number);
  // Checks
  const rankGeom = geometricRank(Number(rho_geom), records);
  if (rankGeom < Number(rank) || rankGeom > lc.length - 1) fail(key, `geometric rank ${rankGeom} against rank ${rank}, deg L ${lc.length - 1}`);
  let euler = 0;
  for (const r of records) euler += r.deg * (r.f + componentCycles(r.type, r.c).reduce((a, b) => a + b, 0));
  if (euler !== 24) fail(key, 'Euler number ' + euler);
  if (p2FromData(lc, records).join(' ') !== P2) fail(key, 'P2 differs from L and the fibres');
  if (n % 50 === 0) {
    const { algebraic } = splitP2(P2.split(' ').map(Number));
    const geo = algebraic.reduce((s, a) => s + a.mult * (a.poly.length - 1), 0);
    const ar = (algebraic.find(a => a.n === 1) || { mult: 0 }).mult;
    if (geo !== Number(rho_geom) || ar !== Number(rho_f2)) fail(key, `Picard numbers: ${geo} ${ar} vs ${rho_geom} ${rho_f2}`);
    const fromL = splitP2(lc).algebraic;
    if (fromL.reduce((s, a) => s + a.mult * (a.poly.length - 1), 0) !== rankGeom) fail(key, 'geometric rank differs from the inverse roots of L');
    sampled++;
  }

  // The zeta function: surfaces are grouped by P2, numbered in order of first appearance.
  const facts = [zeta_sign, rho_f2, rho_geom, height === 'inf' ? '' : height, br_disc, pts_f2, pts_f4, pts_f8, pts_f16, pts_f32].join('\t');
  // The Picard numbers are checked once per zeta function, so for every surface.
  let z = zetas.get(P2);
  if (!z) {
    const p2c = P2.split(' ').map(Number), cyc = cyclotomicPart(p2c);
    const geo = cyc.factors.reduce((s, f) => s + f.e * cyclotomicDegree(f.n), 0);
    const ar = (cyc.factors.find(f => f.n === 1) || { e: 0 }).e;
    if (geo !== Number(rho_geom) || ar !== Number(rho_f2)) fail(key, `Picard numbers: ${geo} ${ar} vs ${rho_geom} ${rho_f2}`);
    zetas.set(P2, z = { id: zetas.size, label: zetaLabel(p2c), n: 0, first: key, facts, cyc, trc: transcendentalId(p2c) });
  } else if (z.facts !== facts) fail(key, 'differs from surface ' + z.first + ', which has the same P2, in an invariant of P2');
  z.n++;

  const config = configuration(records);
  const nSing = records.reduce((s, r) => s + r.deg, 0);
  const nRed = records.reduce((s, r) => s + (isReducible(r.type) ? r.deg : 0), 0);
  const hText = height === 'inf' ? '' : height;
  await write(surfacesOut, [key, a1, torsion, config, records.length, nSing, nRed, lc.length - 1, L_sign, rank, rankGeom, zeta_sign,
    rho_f2, rho_geom, hText, br_disc, pts_f2, pts_f4, pts_f8, pts_f16, pts_f32, z.id, z.cyc.text, z.cyc.nsDegree].join('\t') + '\n');
  for (const r of records) await write(fibresOut, [key, r.deg, r.f, r.c, r.type].join('\t') + '\n');

  shard.push([key, a1, torsion, fibres, L, L_sign, rank, zeta_sign, rho_f2, rho_geom, height, br_disc,
              pts_f2, pts_f4, pts_f8, pts_f16, pts_f32].join(','));
  if (shard.length === 1000) flushShard();

  // Statistics
  const a = byA1[a1] ||= { models: 0, rank: 0, rankGeom: 0, singular: 0, maxSingular: 0, reducible: 0, maxReducible: 0, rho: 0, supersingular: 0, ordinary: 0, torsion: 0 };
  a.models++; a.rank += Number(rank); a.rankGeom += rankGeom; a.singular += nSing; a.maxSingular = Math.max(a.maxSingular, nSing);
  a.reducible += nRed; a.maxReducible = Math.max(a.maxReducible, nRed); a.rho += Number(rho_geom);
  if (height === 'inf') a.supersingular++;
  if (height === '1') a.ordinary++;
  if (torsion !== 'trivial') a.torsion++;
  const j0 = a1 === '0';
  bump(H.rank, Number(rank)); bump(H.rank_geom, rankGeom); bump(H.rho_f2, Number(rho_f2)); bump(H.rho_geom, Number(rho_geom));
  bump(H.L_sign, L_sign); bump(H.zeta_sign, zeta_sign); bump(H.torsion, torsion);
  bump(H.br, Math.floor(Math.log2(Number(br_disc)))); bump(H.config, config);
  bump(H.L_degree, lc.length - 1); bump(H.n_places, records.length); bump(H.ns_degree, z.cyc.nsDegree);
  const split = (m, x, second) => { const e = m.get(x) || [0, 0]; e[second ? 1 : 0]++; m.set(x, e); };
  split(H.singular, nSing, j0); split(H.reducible, nRed, j0);
  split(H.height, height === 'inf' ? 'inf' : Number(height), a1 === 's*t' || a1 === 's^2 + s*t + t^2');
  for (const [name, v] of [['pts_f2', pts_f2], ['pts_f4', pts_f4], ['pts_f8', pts_f8], ['pts_f16', pts_f16], ['pts_f32', pts_f32]]) {
    bump(H[name], Math.floor(Number(v) / WIDTH[name]) * WIDTH[name]);
  }
  if (++n % 200000 === 0) console.log(n, 'rows');
}
flushShard();
await new Promise(res => surfacesOut.end(res));
await new Promise(res => fibresOut.end(res));
const zetaOut = fs.createWriteStream(path.join(tmpDir, 'zeta.tsv'));
const cycOut = fs.createWriteStream(path.join(tmpDir, 'cyclotomic.tsv'));
const labels = new Set();
const parts = { cyc: new Map(), trc: new Map() };     // part -> surfaces and zeta functions with it
for (const [P2, z] of zetas) {
  if (labels.has(z.label)) fail(z.first, 'shares the label ' + z.label + ' with another zeta function');
  labels.add(z.label);
  await write(zetaOut, [z.id, z.label, P2, z.facts, z.n, z.first, z.cyc.text, z.cyc.nsDegree, z.trc].join('\t') + '\n');
  for (const [m, part] of [[parts.cyc, z.cyc.text], [parts.trc, z.trc]]) {
    const e = m.get(part) || { surfaces: 0, zetas: 0 };
    e.surfaces += z.n; e.zetas++;
    m.set(part, e);
  }
  for (const f of z.cyc.factors) await write(cycOut, [z.id, f.n, f.e].join('\t') + '\n');
}
await new Promise(res => zetaOut.end(res));
await new Promise(res => cycOut.end(res));
console.log('rows:', n, ' problems:', problems, ' Picard numbers checked for all', zetas.size, 'zeta functions; geometric ranks on a sample of', sampled);
if (problems) throw new Error('stopping: fix the problems above before publishing');

// ---------------------------------------------------------------------------------------
// Parquet files
// ---------------------------------------------------------------------------------------

const db = await DuckDBInstance.create(':memory:');
const con = await db.connect();
const P = (f) => path.join(tmpDir, f).replace(/\\/g, '/');
const O = (f) => path.join(dataDir, f).replace(/\\/g, '/');
await con.run(`COPY (SELECT * FROM read_csv('${P('surfaces.tsv')}', delim='\t', header=false, quote='', auto_detect=false, columns={
    'key':'UINTEGER','a1':'VARCHAR','torsion':'VARCHAR','config':'VARCHAR','n_places':'UTINYINT','n_singular':'UTINYINT','n_reducible':'UTINYINT',
    'L_degree':'UTINYINT','L_sign':'TINYINT','analytic_rank':'UTINYINT','rank_geom':'UTINYINT','zeta_sign':'TINYINT','rho_f2':'UTINYINT','rho_geom':'UTINYINT',
    'height':'UTINYINT','br_disc':'UINTEGER','pts_f2':'USMALLINT','pts_f4':'USMALLINT','pts_f8':'USMALLINT','pts_f16':'USMALLINT','pts_f32':'USMALLINT','zid':'UINTEGER','cyc':'VARCHAR','ns_degree':'USMALLINT'})
  ORDER BY key) TO '${O('surfaces.parquet')}' (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 50000)`);
// Sorted by type, so that a filter on one Kodaira type reads only that part of the file.
await con.run(`COPY (SELECT * FROM read_csv('${P('fibres.tsv')}', delim='\t', header=false, quote='', auto_detect=false, columns={
    'key':'UINTEGER','deg':'UTINYINT','f':'UTINYINT','c':'UTINYINT','type':'VARCHAR'})
  ORDER BY type, c, f, deg, key) TO '${O('fibres.parquet')}' (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 50000)`);
// Sorted by the number of surfaces, the order in which the zeta functions tab first shows them.
await con.run(`COPY (SELECT * FROM read_csv('${P('zeta.tsv')}', delim='\t', header=false, quote='', auto_detect=false, columns={
    'zid':'UINTEGER','label':'VARCHAR','p2':'VARCHAR','zeta_sign':'TINYINT','rho_f2':'UTINYINT','rho_geom':'UTINYINT',
    'height':'UTINYINT','br_disc':'UINTEGER','pts_f2':'USMALLINT','pts_f4':'USMALLINT','pts_f8':'USMALLINT','pts_f16':'USMALLINT','pts_f32':'USMALLINT',
    'n_surfaces':'UINTEGER','first_key':'UINTEGER','cyc':'VARCHAR','ns_degree':'USMALLINT','trc':'VARCHAR'})
  ORDER BY n_surfaces DESC, zid) TO '${O('zeta.parquet')}' (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 50000)`);
// One row per cyclotomic factor of each zeta function, for the cyclotomic factor filters.
await con.run(`COPY (SELECT * FROM read_csv('${P('cyclotomic.tsv')}', delim='\t', header=false, quote='', auto_detect=false, columns={
    'zid':'UINTEGER','n':'UTINYINT','e':'UTINYINT'})
  ORDER BY n, e, zid) TO '${O('cyclotomic.parquet')}' (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 50000)`);

// ---------------------------------------------------------------------------------------
// Everything else
// ---------------------------------------------------------------------------------------

const missing = allKeys.filter(k => !have.has(k)).sort((x, y) => x - y);
fs.writeFileSync(path.join(dataDir, 'missing_keys.txt'), missing.join('\n') + '\n');
fs.writeFileSync(path.join(rowsDir, 'index.json'), JSON.stringify(firstKeys));
const sorted = (m) => [...m].sort((x, y) => (Number(x[0]) - Number(y[0])) || (String(x[0]) < String(y[0]) ? -1 : 1));
const numeric = (m) => sorted(m).map(([x, v]) => [x, v]);
const stats = {
  rank: numeric(H.rank), rank_geom: numeric(H.rank_geom), rho_f2: numeric(H.rho_f2), rho_geom: numeric(H.rho_geom), L_sign: numeric(H.L_sign),
  zeta_sign: numeric(H.zeta_sign), L_degree: numeric(H.L_degree), n_places: numeric(H.n_places), ns_degree: numeric(H.ns_degree),
  torsion: [...H.torsion].sort((x, y) => y[1] - x[1]),
  br: numeric(H.br), singular: numeric(H.singular), reducible: numeric(H.reducible),
  height: [...H.height].sort((x, y) => (x[0] === 'inf') - (y[0] === 'inf') || x[0] - y[0]),
  pts_f2: numeric(H.pts_f2), pts_f4: numeric(H.pts_f4), pts_f8: numeric(H.pts_f8), pts_f16: numeric(H.pts_f16), pts_f32: numeric(H.pts_f32),
  widths: WIDTH, byA1,
  configurations: { distinct: H.config.size, top: [...H.config].sort((x, y) => y[1] - x[1]).slice(0, 30) },
  // The 30 most common cyclotomic and transcendental parts: [part, surfaces, zeta functions]
  ...Object.fromEntries(Object.entries(parts).map(([name, m]) => [name + 'Parts', {
    distinct: m.size,
    top: [...m].sort((x, y) => y[1].surfaces - x[1].surfaces || y[1].zetas - x[1].zetas).slice(0, 30).map(([p, e]) => [p, e.surfaces, e.zetas]),
  }])),
};
fs.writeFileSync(path.join(dataDir, 'stats.json'), JSON.stringify(stats));
// The release name, from the census.manifest.json next to the CSV when there is one.  The census
// commit is CENSUS_COMMIT when set (for data taken from a commit that is not checked out).
const releaseFile = path.join(path.dirname(csvPath), 'census.manifest.json');
const release = fs.existsSync(releaseFile) ? JSON.parse(fs.readFileSync(releaseFile, 'utf8')).release : null;
const cycOrders = new Set();
for (const z of zetas.values()) for (const f of z.cyc.factors) cycOrders.add(f.n);
const mb = (f) => Math.round(fs.statSync(path.join(dataDir, f)).size / 1e5) / 10;
fs.writeFileSync(path.join(dataDir, 'manifest.json'), JSON.stringify({
  rows: n, models: allKeys.length, missing: missing.length, zetas: zetas.size, expectedByA1,
  rowsByA1: Object.fromEntries(Object.entries(byA1).map(([a, v]) => [a, v.models])),
  torsionGroups: [...H.torsion.keys()].sort(),
  kodairaTypes: [...new Set([...H.config.keys()].flatMap(c => c.split(', ').map(x => x.replace(/^\d+ /, ''))))],
  cyclotomicOrders: [...cycOrders].sort((x, y) => x - y),
  sizesMB: { surfaces: mb('surfaces.parquet'), fibres: mb('fibres.parquet'), zeta: mb('zeta.parquet'), cyclotomic: mb('cyclotomic.parquet') },
  source: path.basename(csvPath), release,
  censusCommit: process.env.CENSUS_COMMIT || execSync('git rev-parse --short HEAD', { cwd: repo }).toString().trim(),
  built: new Date().toISOString().slice(0, 10),
}, null, 1));
fs.copyFileSync(path.join(here, 'census_math.mjs'), path.join(siteDir, 'js', 'census_math.mjs'));
fs.rmSync(tmpDir, { recursive: true });
for (const f of ['surfaces.parquet', 'fibres.parquet', 'zeta.parquet', 'cyclotomic.parquet', 'missing_keys.txt', 'stats.json', 'manifest.json']) {
  console.log(f.padEnd(20), (fs.statSync(path.join(dataDir, f)).size / 1e6).toFixed(2), 'MB');
}
console.log('row files:', shardNo, ' missing keys:', missing.length, ' zeta functions:', zetas.size);
