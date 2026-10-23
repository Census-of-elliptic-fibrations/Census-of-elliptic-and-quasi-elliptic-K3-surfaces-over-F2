# Census of elliptic and quasi-elliptic K3 surfaces over F2

This Magma code computes the 1,397,598 isomorphism classes of elliptic K3 surfaces with a section over $\mathbf F_2$ (zero section fixed, automorphisms of the base allowed) in under two minutes. It also lists the 720 Weierstrass models with $a_1 = a_3 = 0$, which give the quasi-elliptic fibrations, a phenomenon of characteristic 2 and 3.

## What is counted

Each surface has a Weierstrass model $y^2 + a_1xy + a_3y = x^3 + a_2x^2 + a_4x + a_6$ in $\mathbf P(1,1,4,6)$. Here $a_i$ is a binary form of degree $2i$ in $s,t$ over $\mathbf F_2$. Models are identified under $\mathrm{GL}(2,\mathbf F_2)$ acting on $(s,t)$ and the changes $x \mapsto x+q$, $y \mapsto y+px+r$, with $p,q,r$ forms of degree 2, 4, 6. A model with $\Delta \ne 0$ that is minimal at every place is an elliptic K3 surface.

When $a_1 = a_3 = 0$, the substitution $x \mapsto x + a_2$ removes $a_2$ and leaves $y^2 = x^3 + a_4x + a_6$. In characteristic 2 this has $\Delta = 0$ identically: the generic fibre over $\mathbf F_2(T)$ is regular but geometrically a cuspidal cubic, and every fibre is singular. Such a model is not counted among the elliptic surfaces above. When it is a K3 surface, its fibration is quasi-elliptic. The remaining changes of coordinates reduce $a_4$ and $a_6$ to a normal form without the monomials $s^8, s^4t^4, t^8$ in $a_4$ and the even monomials $s^{12-2k}t^{2k}$ in $a_6$, since these can be absorbed by $x \mapsto x + p^2$, $y \mapsto y + px + r$ (the changes that keep $a_2 = 0$).


## Files

Code:

- `elliptic_fibration_orbits_f2.m`: the orbit enumeration. `FibrationOrbitPrefixesF2` lists the 88 pairs $(a_1,a_3)$, `FibrationFiberOrbitKeysF2` returns one key per orbit over a pair (dropping $\Delta = 0$ by default), and `FibrationFromKeyF2` turns a key back into a model.
- `k3_minimality_f2.m`: `IsEverywhereMinimalF2`, which returns whether a model has $\Delta \ne 0$ and is minimal at the four candidate places, and the places where it is not minimal.
- `orbit-fib.m`: runs both over all 88 pairs and writes the everywhere minimal models to `census.txt`, then compresses it with `gzip` to `census.txt.gz`.
- `functions_census.m`: the functions computing the data collected for each surface, with notes on useful built-in Magma functions. `ZetaFunction(E)` returns $P_2(u)$ and $Z(X,u) = 1/((1-u)P_2(u)(1-q^2u))$ for the elliptic surface with generic fibre $E$, from `LFunction(E)` and the Frobenius action on the components of the reducible fibres. `PicardNumber(E)` returns $\rho(X_{\bar{\mathbf F}_2})$ and $\rho(X_{\mathbf F_2})$ from $P_2$, via the Tate conjecture (proved for elliptic K3 surfaces by Artin and Swinnerton-Dyer).
- `test_zeta_picard.m`: checks `ZetaFunction` and `PicardNumber` on the first 100 curves of `mini_census.m` (degree, functional equation, parity of $\rho$) and prints the distribution of Picard numbers.
- `timing_tests.m`: times built-in Magma functions (`LocalInformation`, `Conductor`, `TorsionSubgroup`, `LFunction`, ...) on the generic fibres of the models in `mini_census.txt`, and estimates the time for the full census. The results go to `timing_report.txt`.

Data:

- `census.txt.gz`: the 1,397,598 elliptic K3 surfaces, one model per line as the list $[a_1,a_2,a_3,a_4,a_6]$ with $a_2 = 0$, e.g. `[0, 0, t^6, s^6*t^2, s^11*t]` (4.3 MB compressed).
- `mini_census.txt`: a random sample of 1000 lines of the census (seed 2026), in the same format, for quick tests.
- `mini_census.m`: the same sample as Magma objects. After `load "mini_census.m";`, `models[i]` is the $i$-th model in `P<s,t>` and `curves[i]` its generic fibre over `K<T>` $= \mathbf F_2(T)$, $T = t/s$.
- `mini_census_stats.csv`: statistics of the 1000 models of `mini_census.txt`, one row per model: the torsion of $E(\mathbf F_2(T))$, the number of bad places and of singular fibres, the Kodaira types of the singular fibres, the degree of `LFunction(E)` and the sign of its functional equation, the analytic rank, and the Picard numbers over $\mathbf F_2$ and $\bar{\mathbf F}_2$.
- `quasi-elliptic.txt`: the 720 models with $a_1 = a_3 = 0$, one per line in the same format, e.g. `[0, 0, 0, 0, s*t^11]`. Which of them are K3 surfaces is not yet checked: the first line, $[0,0,0,0,0]$, is $y^2 = x^3$ and is not.
- `quasi_elliptic.m`: the same models as Magma objects. After `load "quasi_elliptic.m";`, `models[i]` is the $i$-th model and `curves[i]` its generic fibre $y^2 = x^3 + a_4x + a_6$ as a plane curve over $\mathbf F_2(T)$. That curve is not smooth, so `EllipticCurve` and the functions built on it do not apply.

## List of tasks that have been done
- Count the 1,397,598 elliptic K3 surfaces over $\mathbf F_2$ (with a section, up to isomorphism).
- List the 720 models with $a_1 = a_3 = 0$ (up to isomorphism), among which are the quasi-elliptic K3 surfaces.
- Program the following functions for the elliptic surfaces:
  - `ZetaFunction` 
  - `PicardNumber` 
  - `LFunctionSign` (sign of the functional equation of the L-function)

## Wish list
These are functions that would be nice to have, but are not yet implemented:
- Singular fibres of the quasi-elliptic surfaces over $\mathbf F_2$.
- Artin invariants of the quasi-elliptic surfaces over $\mathbf F_2$.
- Product of the size of the Brauer group and the discriminant of the Néron-Severi lattice of the elliptic surfaces over $\mathbf F_2$ (through the Shioda-Tate formula).
- Angle ranks?
- Brauer groups of the surfaces.
- Mordell-Weil groups of the fibrations.
- Néron-Severi lattices of the surfaces.
- Artin height of the surfaces.


## Method

In characteristic 2 the change $x \mapsto x+a_2$ kills $a_2$, leaving 32 coefficient bits. The rest is the staged quotient and stabilizer enumeration of Lemma 1.1 and Algorithm 1.2 in Auel, Kulkarni, Petok and Weinbaum, "A census of cubic fourfolds over $\mathbf F_2$", [arXiv:2306.09908](https://arxiv.org/abs/2306.09908). It first finds the 88 orbits of pairs $(a_1,a_3)$. For each pair, its stabilizer acts affinely on the quotient of the $(a_4,a_6)$ space by the shifts $y \mapsto y+r$, which has at most $2^{16}$ points. Non-minimality at a place forces $v(\Delta) \ge 12$ and $\deg \Delta = 24$, so minimality only needs checking at $s$, $t$, $s+t$ and $s^2+st+t^2$.

## Usage

Run the scripts from the repository root. They do not end with `quit;`, so they can also be loaded into a Magma session with `load "file.m";`; to run one and exit, use `magma -b file.m < /dev/null`.

Run `magma -b orbit-fib.m < /dev/null` to write the list of surfaces to `census.txt.gz`. It needs `gzip` on the path. Read the file with `zcat census.txt.gz`.

Run `magma -b timing_tests.m < /dev/null` to time the built-in functions on the sample, and `magma -b test_zeta_picard.m < /dev/null` to check `ZetaFunction` and `PicardNumber`. To try functions by hand, `load "mini_census.m";` (or `load "quasi_elliptic.m";`) and `load "functions_census.m";` in a Magma session.

## Results

| Models | Orbits |
| --- | ---: |
| All | 1,402,468 |
| $\Delta = 0$ | 736 |
| $\Delta \ne 0$ | 1,401,732 |
| Non-minimal at exactly one place | 4,129 |
| Non-minimal at two places | 5 |
| Elliptic K3 (everywhere minimal) | 1,397,598 |
| K3 with $a_1 = 0$ ($j = 0$) | 173,977 |
| K3 with $a_1 = s^2$ | 524,493 |
| K3 with $a_1 = st$ | 523,262 |
| K3 with $a_1 = s^2+st+t^2$ | 175,866 |

| Quasi-elliptic case | Models |
| --- | ---: |
| $a_1 = a_3 = 0$, up to isomorphism (`quasi-elliptic.txt`) | 720 |
| Of these, quasi-elliptic K3 surfaces | not yet determined |

## Checks

- The 88 pairs, the 1,402,468 orbits and the 736 with $\Delta = 0$ agree with an independent Burnside count over all $2^{37}$ models.
- The minimality filter agrees with Magma's Tate algorithm (`LocalInformation`) on all 97,036 census representatives with $v(\Delta) \ge 12$ at a candidate place, and on 20,000 random others.
- Not yet done: a comparison with an outside count, and the automorphism groups needed for a mass formula check.

