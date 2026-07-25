# JET_OBSTRUCTION — profiling jet-lifting obstruction depth on the J6 control

> **STATUS (2026-07-24, `jet_obstruction.py`): CONTROL REPRODUCES — and it
> surfaced a live serializer defect.** Two independent things happened.
>
> **(1) The J6 control reproduces.** The obstruction profiler recovers exactly
> the shape `J6_MSOLVE.md` reports: `X_2` non-empty and 0-dimensional, `X_3`
> empty, obstruction localised at the depth-3 step — for all 4 states, over 4
> primes, and (depth 3) over Q via msolve. Regime verdict: **A (bounded
> obstruction)**. `python jet_obstruction.py --quiet` exits 0.
>
> **(2) BLOCKING, and the reason this file exists: the polynomials J6 actually
> handed to its solvers were not the polynomials it meant to hand them.**
> `blowup_diagnosis.sing_poly_intcoeff` — and its mod-p twin
> `modular_triage.poly_to_singular_modp` — **silently truncate rational
> coefficients toward zero** under this environment (sympy 1.14.0 +
> python-flint 0.9.0). The class quartic
> `c0_0^4 + c0_0^3/4 + 5c0_0^2/32 + 15c0_0/128 + 195/2048` is serialized as
> **`c0_0^4`**. 18 of the 22 polynomials in the J6 systems are corrupted.
> Both helpers sit on the char-0 *and* mod-p paths of `j6_msolve.py`,
> `msolve_bridge.py`, `blowup_diagnosis.emit_program`, and (via
> `poly_to_singular_modp`) `alt_bridge.py`, `bridge_sweep.py`,
> `bridge_sweep_verify.py`, `full_system_bridge.py`,
> `full_system_bridge_verify.py`, `modular_triage.py`.
>
> **(3) The J6 verdicts survive re-derivation.** Re-run with an exact
> serializer, the corrected systems give the same verdicts: depth 2 satisfiable,
> depth 3 the unit ideal, all 4 states. So `J6_MSOLVE.md`'s *conclusion* stands;
> its *computation* does not, and every other result routed through those two
> helpers needs re-checking. **No existing file was modified by this lane.**

New files (uncommitted): `jet_obstruction.py`, `jet_obstruction_results.json`,
this doc. READ-ONLY on every existing artifact.

---

## 0. What is and is not established

**Established here.**

- The two serializer defects, minimally reproduced (§1). Demonstrated live on
  the J6 systems by the deliverable script: the bugcheck line and the full
  per-polynomial discrepancy list are emitted on every run and recorded under
  `serializer_bug` in the results JSON.
- Provenance: the J6 generators, class relations and saturation are re-derived
  through the original committed code path (`alt_hunt_depth2.reconstruct_general`
  + `ConvolutionDescent.master_coefficient`) and are **string-identical** to
  `j6_msolve_results.json`, 4/4 states (§2). So the *recorded systems* are
  right; only their *serialization to the solver* was wrong.
- The obstruction profile of §4, over F_p for p ∈ {10007, 10009, 32003, 100019},
  computed from exactly-serialized polynomials.
- Over Q: msolve returns `[-1]` (empty in the algebraic closure) on the
  corrected **depth-3** system for all 4 states (§6).

**NOT established here.**

- That the other results routed through the two defective helpers are wrong —
  or right. Not audited. Only the J6 systems were checked.
- The char-0 **depth-2** satisfiability for 3 of the 4 states (§6): msolve was
  run to completion on the corrected depth-2 system for one state only; host
  CPU contention (79–97% from unrelated jobs) starved the rest. Depth-2
  non-emptiness for all four is established **mod p only**.
- Any regime claim beyond depth 3. Only 3 accumulated master coefficients are
  recorded, so "bounded obstruction" is a statement about this truncation, not
  a proof about arbitrarily deep jets.
- A harder second target. Not reached — see §8.

---

## 1. The serializer defect

`blowup_diagnosis.sing_poly_intcoeff` intends to clear one global rational
denominator and emit an integer-coefficient Singular/msolve string. Its
denominator-clearing step is

```python
expr = sp.cancel(sp.sympify(expr))
num, den = sp.fraction(expr)
...
for monom, coeff in poly.terms():
    c = int(coeff)
```

In this environment `sp.cancel` does **not** put a sum over a common
denominator, so `den` is 1 and `int(coeff)` truncates. Minimal repro, actually
run:

```
sympy 1.14.0
python-flint 0.9.0
cancel(x/2+1/3) = x/2 + 1/3 | fraction den = 1
cancel(quartic) = c0_0**4 + c0_0**3/4 + 5*c0_0**2/32 + 15*c0_0/128 + 195/2048
fraction den   = 1
Poly terms     = [((4,), 1), ((3,), 1/4), ((2,), 5/32), ((1,), 15/128), ((0,), 195/2048)]
int() of each  = [1, 0, 0, 0, 0]
```

`modular_triage.poly_to_singular_modp` has the identical structure and the
identical defect. Both, called on live J6 objects:

```
expr                       : c0_0**4 + c0_0**3/4 + 5*c0_0**2/32 + 15*c0_0/128 + 195/2048
bd.sing_poly_intcoeff      : c0_0^4
mt.poly_to_singular_modp(p): c0_0^4

expr                       : 3981312*E**21/221 - 2305843009213693952*E**8/400329564123571875
bd.sing_poly_intcoeff      : 18014*E^21-5*E^8
mt.poly_to_singular_modp(p): 8007*E^21+10002*E^8
```

The correct reductions mod 10007 are `9140*E^21 + 6635*E^8`; the helper emits
`8007*E^21 + 10002*E^8` (i.e. `int(3981312/221) = 18014` reduced, and
`int(-5.76…) = -5` reduced). This is not a rounding nuisance — for the class
relations it deletes every term but the leading one, replacing a separable
quartic by a nilpotent `c0_0^4 = 0`.

`jet_obstruction.py` reports the discrepancy on every run:

```
== serializer bugcheck: 18 of 22 polynomials differ from blowup_diagnosis.sing_poly_intcoeff
```

| state | polynomial | exact serializer | `sing_poly_intcoeff` |
|---|---|---|---|
| `…a9…#state3` | class_relation[0] | `2048*c0_0^4+512*c0_0^3+320*c0_0^2+240*c0_0+195` | `c0_0^4` |
| `…a8…#state0` | class_relation[0] | `256*c0_0*c0_1+32*c0_0-128*c0_1^3-32*c0_1^2-20*c0_1-15` | `2*c0_0*c0_1-1*c0_1^3` |
| `…a8…#state0` | class_relation[1] | `2048*c0_0^2-2048*c0_0*c0_1^2-512*c0_0*c0_1-320*c0_0+195` | `c0_0^2-1*c0_0*c0_1^2` |

Only the saturation `E*X*w-1` (integral) survives unchanged, in all 4 states —
hence 18 of 22.

`jet_obstruction.py` uses its own serializer, which clears the lcm of the
coefficient denominators exactly and then **verifies the emitted string by
reparsing it** and asserting equality with `expr * den`.

---

## 2. Method

For a state, the depth-N truncation is the affine scheme built exactly as
`j6_msolve.py` builds it (accumulated master coefficients + class relations +
Rabinowitsch saturation, the latter two present at every depth):

    X_N = V( f_1, …, f_N, class relations, w·Π(scalars) − 1 ),   A_N = k[vars]/I_N

and π_N : X_{N+1} → X_N forgets the newly-introduced coefficients. The
first-order obstruction at x ∈ X_N is the cokernel of the derivative of the new
equations with respect to the new unknowns.

**Degeneracy, stated loudly.** In this construction the class relations and the
saturation are present at *every* depth, so the ambient variable set does not
grow with N: **the set of new unknowns is empty at every step**. The relative
Jacobian is therefore a 1×0 matrix — rank 0, cokernel the whole new-equation
line at every point. π_N is a closed immersion, not a fibration. That is not a
profiler bug; it is the reason the correct obstruction object here is the
multiplication operator

    M_{f_{N+1}} : A_N → A_N,     coker M_{f_{N+1}} = A_{N+1}

whose rank and cokernel are plain linear algebra over F_p on a deg(X_N)-square
matrix. `X_{N+1}` is empty ⟺ `M_{f_{N+1}}` is invertible ⟺ `f_{N+1}` is a
**unit** in `A_N` ⟺ `f_{N+1}` vanishes at no point of `X_N`.

The profiler nevertheless records the *first-appearance filtration* of the
unknowns, because that is the datum a depth-staged construction would turn into
actual columns:

| state | E | X | c0_1 | c0_0 |
|---|---|---|---|---|
| `a9…#state3` | 1 | 1 | — | 2 |
| `a8…#state0` | 1 | 1 | 2 | 3 |

So `a8` *does* introduce `c0_0` only at depth 3 — but the class relations pin it
from depth 1 onward, so it is not new in the ideal-theoretic sense.

**Pipeline.** Singular (WSL) computes `std`, `dim`, `vdim`, `kbase`, and the
normal forms `NF(f·b_j)` for `b_j` in the standard monomial basis; numpy does
Gaussian elimination mod p for the rank, cokernel and determinant. Point
structure uses the multiplication matrices of the variables (so any linear form
is a free numpy combination) plus a Krylov minimal-polynomial computation.

**Provenance.** Generators are read verbatim from `j6_msolve_results.json`;
`--replay` re-derives them through the original committed code path and asserts
string-identity. Result, actually run:

```
  replay sub2:a9_b1000_T1_sz1_dz1_gz-#state3: OK (25.1s)
  replay sub2:a8_b1100_T1_sz1_dz1_gz-#state0: OK (23.9s)
  replay sub1:a9_b1000_T1_sz1_dz1_gz-#state3: OK (23.8s)
  replay sub1:a8_b1100_T1_sz1_dz1_gz-#state0: OK (23.9s)
  replay: 4/4 verbatim
```

(polys, all three master coefficients, class relations and saturation each
compared as strings.)

---

## 3. Primes, and why they are safe

Forbidden primes are computed per state, not assumed: a prime is unsafe if it
divides a coefficient denominator (reduction undefined) or the integer content
of the denominator-cleared generator (the generator would reduce to 0). For all
four states the profiler reports

```
forbidden [2, 3, 5, 13, 17]
```

Independently, the class quartic `2048c^4+512c^3+320c^2+240c+195` is
**irreducible over Q** with discriminant

```
disc = 12837954459480883200 = 2^36 · 3^2 · 5^2 · 13^3 · 17^3
```

— the same prime set. So for every prime used the marked-root scheme has good,
separable reduction; the quartic stays squarefree and the class algebra keeps
its char-0 degree. Primes used: **10007, 10009, 32003, 100019** (the first
three continue `J6_MSOLVE.md` §5's choice; 32003 is Singular's default). The
saturation factors are `E` and `X`, and their non-vanishing is imposed by the
generator `E*X*w−1` inside the ideal, so it holds in every characteristic by
construction rather than by luck.

Raw Singular confirmation that the class relation now survives serialization
(depth-1 GB, `sub2:a9`, p = 10007):

```
@GB1|E*X*w-1,c0_0^4+2502*c0_0^3-938*c0_0^2+4300*c0_0+992,X^10*w^9+4957*E^17*w-3377*E^13*X^4-…
```

`c0_0^4+2502*c0_0^3-938*c0_0^2+4300*c0_0+992` is the monic quartic mod 10007.
Under the defective serializer this generator was `c0_0^4`.

---

## 4. The obstruction profile

Identical at all four primes. `deg X_N` is `dim_k A_N` (`vdim`); `inf` means
positive-dimensional. "rel rank" is the rank of the relative Jacobian w.r.t.
new unknowns (a 1×0 matrix — see §2); "obstruction rank/coker" are those of
`M_{f_{N+1}}` on `A_N`.

| state | depth | master deg | #eqs | #new unk | dim X_N | deg X_N | rel rank | obstruction rank | coker dim | lifts |
|---|---|---|---|---|---|---|---|---|---|---|
| `sub2:a9…#state3` | 1 | 250 | 3 | 0 | 1 | inf | 0 | n/a | n/a | yes |
| `sub2:a9…#state3` | 2 | 249 | 4 | 0 | 0 | 680 | 0 | 680 | **0** | **NO** |
| `sub2:a9…#state3` | 3 | 248 | 5 | — | −1 (empty) | 0 | — | — | — | — |
| `sub2:a8…#state0` | 1 | 250 | 4 | 0 | 1 | inf | 0 | n/a | n/a | yes |
| `sub2:a8…#state0` | 2 | 249 | 5 | 0 | 0 | 1020 | 0 | 1020 | **0** | **NO** |
| `sub2:a8…#state0` | 3 | 248 | 6 | — | −1 (empty) | 0 | — | — | — | — |
| `sub1:a9…#state3` | 1 | 250 | 3 | 0 | 1 | inf | 0 | n/a | n/a | yes |
| `sub1:a9…#state3` | 2 | 249 | 4 | 0 | 0 | 680 | 0 | 680 | **0** | **NO** |
| `sub1:a9…#state3` | 3 | 248 | 5 | — | −1 (empty) | 0 | — | — | — | — |
| `sub1:a8…#state0` | 1 | 250 | 4 | 0 | 1 | inf | 0 | n/a | n/a | yes |
| `sub1:a8…#state0` | 2 | 249 | 5 | 0 | 0 | 1020 | 0 | 1020 | **0** | **NO** |
| `sub1:a8…#state0` | 3 | 248 | 6 | — | −1 (empty) | 0 | — | — | — | — |

The two windows (`sub1`, `sub2`) give numerically identical profiles — same
degrees, same ranks, same determinants — consistent with them carrying the same
admissible split.

**Step 1→2** is a proper divisor cut: `dim` drops 1 → 0, `A_1` is
infinite-dimensional so no finite multiplication matrix exists; the profiler
reports the dimension drop instead of inventing a rank. Lifts exist.

**Step 2→3 is the obstruction.** `M_{f_3}` has full rank on the 680- (resp.
1020-) dimensional algebra `A_2`, cokernel 0, and nonzero determinant:

| state | p=10007 | p=10009 | p=32003 | p=100019 |
|---|---|---|---|---|
| `a9…#state3` | 4716 | 314 | 324 | 28165 |
| `a8…#state0` | 5088 | 476 | 6928 | 35828 |

`det M_{f_3} = Π f_3(x_i)^{m_i}` over the points of `X_2`, so a nonzero
determinant is exactly the statement **"the depth-3 master coefficient vanishes
at no point of `X_2`"**. That is the geometric replacement for "the Gröbner
basis timed out": the obstruction is not a subvariety cutting `X_2` down, it is
a nowhere-vanishing function on it.

The linear algebra agrees with an independent Gröbner computation in every
case (`coker_matches_vdim_next: true`, 4 states × 4 primes). Raw markers
(`sub2:a9`, p = 10007):

```
@DIM1|1
@VDIM1|-1
@SIZE1|5
@DIM2|0
@VDIM2|680
@SIZE2|114
@DIM3|-1
@VDIM3|0
@SIZE3|1
@GB3|1
@KN2|680
```

`@GB3|1` is Singular's reduced GB of the depth-3 ideal: the unit ideal.

Honest note on cost: at these sizes mod p, `std(I_3)` is itself cheap
(sub-second), so the linear-algebra route is **not** a speed win here. Its value
is that it returns the obstruction *object* — an invertible operator with a
determinant — rather than a yes/no.

---

## 5. Point structure of `X_2`

By pure linear algebra on the multiplication matrices (no factorization of the
big system), cross-checked against an FGLM lex eliminant:

| state | deg `X_2` (with multiplicity) | distinct points (lower bound) | `A_2` reduced? | `E`-eliminant degree |
|---|---|---|---|---|
| `a9…#state3` | 680 | ≥ 476 | **no** | 374 |
| `a8…#state0` | 1020 | ≥ 714 | **no** | 510 |

The lower bound is `deg` of the squarefree part of the minimal polynomial of a
random linear form (2 independent random forms, agreeing); it equals the true
count iff that form separates the points, so it is a **bound**, not a count.
Non-reducedness *is* proved: a finite reduced algebra has squarefree minimal
polynomials for every element, and these do not (680 vs 476, 1020 vs 714). So
`X_2` carries genuine multiplicity.

The FGLM cross-check (`stdfglm` to lex, then `factorize`) gives the univariate
eliminant in `E` of degree 374 (a9) / 510 (a8); its factor-degree multiset is
prime-dependent, as expected, while the degree is stable across all four primes.

---

## 6. Char-0 re-run with the exact serializer

`j6_msolve.py`'s msolve inputs were built by the defective serializer, so its
char-0 verdicts are computed on a different system than intended. Re-running
msolve over Q on the **corrected** systems (same script, `--msolve`):

| state | depth | verdict | wall | msolve raw output |
|---|---|---|---|---|
| `sub2:a9…#state3` | 2 | NOT_EMPTY | 233.1 s | `[0, [1,\n[]\n]]:` |
| `sub2:a9…#state3` | 3 | **EMPTY** | 0.6 s | `[-1]:` |
| `sub2:a8…#state0` | 3 | **EMPTY** | 0.7 s | `[-1]:` |
| `sub1:a9…#state3` | 3 | **EMPTY** | 0.5 s | `[-1]:` |
| `sub1:a8…#state0` | 3 | **EMPTY** | 1.0 s | `[-1]:` |

`[-1]` is msolve's "no solution in the algebraic closure of Q". **So all four
J6 kills at depth 3 survive re-derivation with exact arithmetic.**

Two honest observations:

- The corrected depth-2 system cost msolve **233.1 s**, against the **1.9 s**
  recorded in `J6_MSOLVE.md` §2 for the same state. The recorded J6 timings are
  timings of the truncated systems, not of the intended ones.
- Depth-2 msolve was completed for one state only; the other three were starved
  by unrelated host load and are **not reported**. Depth-2 non-emptiness for all
  four rests on the mod-p profile of §4 (`dim = 0`, `deg = 680`/`1020` > 0 at
  four primes).
- I did not chase msolve's exact output grammar for `[0, [1,\n[]\n]]:`; the
  empty inner list plausibly denotes zero *real* solutions, which is consistent
  with §4's 0-dimensional degree-680 scheme over `F̄_p`. What is load-bearing is
  only that it is not `[-1]`.

---

## 7. Regime verdict

**Regime A — bounded obstruction — for all four states**, at every prime tried
and (depth 3) over Q:

> `X_3` is empty; beyond depth 2 nothing lifts. The obstruction is carried by a
> single element `f_3` that is a **unit** in the 680- (resp. 1020-) dimensional
> coordinate algebra `A_2` — it vanishes nowhere on `X_2`, and the cokernel of
> multiplication by it is 0.

This is the strong form of regime A: not a divisor trimming `X_2`, but a
nowhere-zero function killing it outright. Regimes B (eventual formal
smoothness) and C (periodic/templated) are **not** in play here, and the
profiler would have reported `B?` or `UNDETERMINED` had no depth been empty.

Scope: this is a statement about the 3-term truncation that J6 records. It says
nothing about deeper jets, because deeper master coefficients were never
computed — see §8.

---

## 8. What was NOT checked

- **The blast radius of the serializer defect.** Only the J6 systems were
  checked. `poly_to_singular_modp` and `sing_poly_intcoeff` are called from
  `alt_bridge.py`, `bridge_sweep.py`, `bridge_sweep_verify.py`,
  `full_system_bridge.py`, `full_system_bridge_verify.py`, `modular_triage.py`,
  `msolve_bridge.py`, `blowup_diagnosis.py`, `j6_msolve.py`. Whether any of
  those pass polynomials with non-integral rational coefficients — and are
  therefore also affected — was **not** audited. The defect is inert on
  integral input.
- **Nothing was fixed.** Per lane scope, no existing file was modified.
  `jet_obstruction.py` carries its own verified serializer; the two defective
  helpers are untouched and still on their callers' paths.
- **Depth > 3.** Only 3 accumulated master coefficients exist in the record.
  "Bounded obstruction" is established for this truncation; the profiler cannot
  and does not claim a bound on all jets.
- **Char-0 depth-2 for 3 of 4 states** (§6).
- **mod-p ⇒ char-0 is not claimed.** Unit ideal at four primes corroborates but
  does not prove the char-0 verdict; the load-bearing char-0 statement is
  msolve's `[-1]`, which was obtained for depth 3 in all four states. This
  inherits `J6_MSOLVE.md` [J3]'s stance verbatim.
- **The distinct-point counts are lower bounds**, not counts (§5).
- **No kill certificate was extracted.** A Singular `lift` of `1` through the
  depth-3 ideal was attempted mod p and **timed out at 600 s** (`rc 124`) — the
  cofactor certificate route remains open, as `J6_MSOLVE.md` §6 [J1] already
  warned.
- **No harder target was profiled.** The intended second target (a small mixed
  phase-D cell) was not reached: no artifact in the repo other than
  `j6_msolve_results.json` records an explicit per-depth generator set
  (`phase_d_states_sub{1,2}.json` carry only degree data — `deg_d1`, `deg_d2`,
  `deg_e`, `deg_g`), so a second target requires re-deriving systems through the
  cascade engine, which additionally routes through the two defective
  serializers. That is the natural next lane, and it should be run *after* the
  serializers are repaired.
- **`sub1` vs `sub2` identity.** The two windows produce numerically identical
  profiles; this was observed, not explained.

---

## 9. Reproduction

```
python jet_obstruction.py --quiet                 # self-check, exit 0 iff the
                                                  #   J6 control reproduces
python jet_obstruction.py                         # full profile + table + JSON
python jet_obstruction.py --replay                # + re-derive gens from source
python jet_obstruction.py --msolve                # + char-0 msolve re-run
JO_MSOLVE_DEPTHS=3 JO_MSOLVE_CAP=300 python jet_obstruction.py --msolve
```

Env: `JO_PRIMES` (default `10007,10009,32003,100019`), `JO_SING_CAP` (per-call
WSL wall s, default 300), `JO_MSOLVE_CAP`, `JO_MSOLVE_DEPTHS`.
Output: `jet_obstruction_results.json`. Singular and msolve run WSL-side under
`timeout` (msolve additionally under `ulimit -v 8G`), and the relay catches
`TimeoutExpired` so a Windows-side kill cannot leave the run in an undefined
state. Actual self-check on this machine:

```
SELF-CHECK PASS: J6 control reproduces (depth-2 satisfiable, depth-3 empty, obstruction at step 2->3)
EXIT=0
```
