# Phase F, work item F2 — the divisor-reconstruction pilot (alternate regime)

**Date:** 2026-07-22 (late).
**Status:** EXPERIMENT — pilot computations, not a branch-closing claim.
**New files (uncommitted):** `phase_f2_pilot.py`, `phase_f2_pilot_verify.py`,
this document. READ-ONLY on every audited artifact.
**Engines:** `phase_f2_pilot.py` (reconstruct + saturated Gröbner);
`phase_f2_pilot_verify.py` (independent re-derivation of the Pilot A chain,
**ALL CHECKS PASS**). Both build on the audited level-0 graded coefficient
`h_0 = cascade_engine.MONOMIALS[0]` (26 monomials, weight-20 homogeneous,
weights `(d2,d1,sigma,e) = (2,3,4,5)`).

## 0. The question and the one-line answer

`ALT_RESIDUE_CONGRUENCES.md` proved that every alternate-regime survivor's
**only** depth-1 obligation is a level-0 leading cancellation of `h_0` at
infinity (the bottom close `E^21 h_0 + u r_0 = 0`), and that all 19 distinct
`h_0` initial forms carry an **all-nonzero rational point** — *because the
leading coefficients `(D,X,S,E) = (lc d2, lc d1, lc sigma, lc e)` were treated as
FREE*. Phase F asks whether that freedom survives once the finite-place divisors
are honoured.

**Answer (this pilot): it does not.** For states whose divisors force
**defect 0** (`deg p = sum_s v_s(p)`), each `p in {d1,sigma,e}` is
`p = lambda_p * prod_s (y-s)^{v_s(p)}` *exactly*, so not just its leading
coefficient but **every** coefficient is a determined multiple of the single
scalar `lambda_p`. Feeding those determined coefficients into the **full**
level-0 tie tower (the depth-1 initial form *plus* the deeper convolution
coefficients the tie depth demands) turns two previously-UNOBSTRUCTED constraints
into **exact KILLS** — the residue system has no solution with the leading
coefficients nonzero, over `Q` or over `Q[r]/(q)`.

**Divisor reconstruction changes the alternate endgame.**

## 1. Reconstruction principle (every step exact)

Places `S = { t : y=-1 ; r_1..r_4 : roots of q = 2048 y^4 - 512 y^3 + 320 y^2 -
240 y + 195 }`. The t-place is the linear place `y=-1`, so `t^m = (y+1)^m`; a
marked q-root `r` contributes `(y-r)^m`. A defect-0 polynomial is therefore

```
p = lambda_p * (y+1)^{v_t(p)} * prod_j (y - r_j)^{v_{r_j}(p)},   deg p = sum of the exponents,
```

with `lambda_p` its single free leading scalar. The finite-place data of a state
(the `ALT_REGIME_L2.md` sec.2 cones and the branch b-pattern
`b_i = v_{r_i}(e)`, plus `v_t(e)=a`) fixes the exponents. At infinity `h_0(y)`
then has degree 60 (the sub1 cap `H_0`), and the bottom close forces its top
`depth` coefficients to vanish, where `depth` is the state's L0 tie depth
(`ALT_RESIDUE_CONGRUENCES.md`). Reading those coefficients as a polynomial system
in `(lambda_p)` and the marked root, and saturating by the leading scalars, is
the Phase-F test. The verdict is read off an exact **saturated Gröbner basis**:
unit ideal after saturating the leading coefficients ⇒ **KILL**.

The top coefficients are computed by truncated convolution of the reversed
factor series (so the degree-60 `h_0` is never expanded in full); the verifier
re-derives them independently and cross-checks the initial form against
`alt_residue_congruences.json`.

## 2. Pilot A — `a11_b3000_T1`, state `(deg d2,d1,sigma,e) = (0,9,12,14)`

**Forced divisors.** `X_min = 9 = deg d1` with the *unique* tight split
`v_t(d1)=5, v_{r}(d1)=4` ⇒ `d1 = X (y+1)^5 (y-r)^4` (defect 0, **forced**).
`deg e = 14 = a + b = 11 + 3` ⇒ `e = E (y+1)^11 (y-r)^3` (defect 0, **forced**).
`d2 = D` (`deg d2 = 0`). Both-tight choice `v_t(sigma)=12` ⇒
`sigma = S (y+1)^12` (defect 0). Support 8, tie depth 14. Here `r` is a single
**marked** root (b3000 activates one q-place); arithmetic in `Q[r]/(q)`.

**Depth-1 (initial form).** Reproduces the audited support-8 form exactly:

```
j0  =  2187 S^2 (4 S^3 + X^4)  = 0      ==>   X^4 = -4 S^3,
```

which — with `(D,X,S,E)` free — has the rational torus point `(X,S)=(4,-4)`
(`4·(-4)^3 + 4^4 = 0`); this is the ALT_RESIDUE "UNOBSTRUCTED" fact.

**Deeper coefficients.** The next tie coefficients `j1, j2, j3` (coefficients of
`y^59, y^58, y^57`) involve **only** `(X,S,E,r)` — not `D` — because the
d2-carrying monomials sit at low y-degree. Writing them out (reduced mod `q(r)`):

```
j1 = -2916 ( 2 E X^5 + 12 E S^3 X - 180 S^5 + 12 S^2 X^4 r - 33 S^2 X^4 ) = 0,
j2 = 1458  ( 21 E^2 S X^2 + ... + 180 S^2 X^4 r^2 - 1056 S^2 X^4 r + 1419 S^2 X^4 ) = 0,
j3 = ...   (a determined element of Q[r]/(q)[X,S,E]).
```

**VERDICT: KILL (exact).** The saturated ideal

```
( j0, j1, j2, j3, q(r), w·X·S·E − 1 )  =  (1)      (unit ideal)
```

so **no** `(X,S,E)` with `X S E != 0` over `Q[r]/(q)` solves even the depth-3
truncation. In fact depth 3 already suffices: `(j0,j1,j2,q,sat)=(1)` too. The
obstruction is **universal in the marked point**: the same ideal with `r` a free
indeterminate (no `q(r)=0`) is still the unit ideal, so it holds for *every*
marked point, not only the roots of `q`.

*Honest scope.* `d1, e` are forced defect-0; `sigma` defect-0 is the both-tight
choice, so this kills the **both-tight (all-defect-0) stratum** of the state. The
sigma-relaxed configuration is Pilot C.

## 3. Pilot B — `a14_b0000_T2`, state `(6, d1≡0, 12, 14)`

**Forced divisors (no choices).** T2 forces `v_t(sigma) = w = 3a−30 = 12 =
deg sigma` ⇒ `sigma = S (y+1)^12` (defect 0, **forced**). `deg e = 14 = a = 14`,
`b0000` ⇒ `e = E (y+1)^14` (defect 0, **forced**). `d1 ≡ 0` (T2). `d2` is a
**free** degree-6 polynomial `sum_{i=0}^6 D_i y^i` (its divisor is unconstrained
in this branch). Support 0, tie depth 14. Everything is over `Q` (no marked
root). Coherence `sigma = 4 d0 − d2^2` merely fixes `d0 = (sigma + d2^2)/4`; `d2`
stays free.

**Depth-1 (initial form).** Reproduces the audited support-0 form:

```
j0 = 12 S^2 (4 D6^2 + 9 S)^2 (5 D6^2 + 9 S) = 0   (D6 = lc d2).
```

On the torus (`S != 0`) this gives two branches `S = −4 D6^2/9` and
`S = −5 D6^2/9`.

**VERDICT: KILL (exact) — genuine whole-state kill.** Fixing `D6 = 1` (the weight
scaling) on each branch and saturating `E != 0` (with `D0..D5` free):

```
branch 5 D6^2 + 9 S :  ideal = (1)   (unit)
branch 4 D6^2 + 9 S :  ideal = (1)   (unit)
```

so **no** completion of the depth-14 tower exists with `lc(d2), lc(sigma), lc(e)`
all nonzero — even though `d2` carries a **full free degree-6 cofactor** (7
coefficients). Because the tightness of `sigma` and `e` here is *forced* (not a
substratum choice), this closes the **entire state** `(6,∅,12,14)`.

## 4. Control — reconstruction is discriminating, not a sledgehammer

Reconstruct `d1` defect-0 as above but leave **e GENERIC** (degree 14, not
reduced to its forced divisor) and `sigma` at its forced floor
`sigma = (y+1)^3 · C`, `C` generic degree 9. The depth-14 tower is now
**triangular**: coefficient `j_k` is linear in a fresh top coefficient of `e`
(and of `C`), with a nonzero pivot at every level. Exact stepwise solution gives

```
all 13 pivots nonzero  ==>  the full depth-14 tie solves for e14..e2 with
lc(d1)=1, lc(sigma)=c9 != 0 :   UNOBSTRUCTED.
```

This is the decisive control: the engine finds solutions **when the freedom is
present**, so the Pilot A/B kills are a genuine consequence of removing that
freedom by divisor reconstruction — not an artefact. It also localises the
mechanism: **reconstructing `e` (and `sigma`) to defect-0 removes exactly the
free top coefficients that solved the deep tie.**

## 5. Pilot C — the sigma-relaxed `a11_b3000_T1 (0,9,12,14)`

Same state as Pilot A but with `sigma` at its **freest** admissible form
(`v_t(sigma)=3`, the a=11 rectangle floor; `sigma = (y+1)^3 · C`, `C` generic
degree 9). Now `d1, e` remain forced defect-0 but `sigma` carries a degree-9
cofactor (10 coefficients). The tower is triangular: `C`'s 10 coefficients and
the initial-form relation absorb `j0..j9`, leaving **4 residual conditions
`j10..j13` in only `(D,E)`** — an over-determined `4-in-2` system over
`Q(r, c9)` (with `X=1`, `4 c9^3 + 1 = 0`).

**VERDICT: obstructed (NARROWED to empty).** An exact per-root resultant/nsolve
search over all four marked roots and all three `c9` branches finds no `(D,E)`
with `E != 0` satisfying the 4 residuals [exact Gröbner certification of the
whole-state kill is the immediate follow-up; the structural cause is transparent:
with `e` **forced** defect-0 the tower-solving top-`e` coefficients of the control
are gone, and `sigma`'s cofactor can only discharge `j0..j9`]. So the entire state
`(0,9,12,14)` — not merely its both-tight stratum — appears to be killed once the
forced defect-0 structure of `d1` and `e` is imposed.

## 6. Assessment — does divisor reconstruction change the alternate endgame?

**Yes, decisively, for the defect-0 (tight-divisor) frontier.**

- The depth-1 layer is unobstructed *precisely because* it ignores the
  divisor-forced sub-leading coefficients. Pilot A exhibits a support (support 8,
  `987`+ states carry the sibling support 15, `438` carry support 13, …) whose
  depth-1 rational point survives but whose depth-3 lift under reconstruction is
  the **unit ideal**. The kill lives at depth 2–3, not depth 1 — exactly the
  layer the census flagged as "necessary but far from sufficient"
  (`ALT_RESIDUE_CONGRUENCES.md` K4).
- Pilot B is a **forced** whole-state kill: `sigma` and `e` are tight with no
  choice, and even a free degree-6 `d2` cannot rescue the tower. This is the
  first alternate-regime *state* closed by the residue congruences that
  `ALT_COMBINED.md` J1 explicitly left "out of scope."
- The control proves the method is discriminating: retained freedom ⇒
  UNOBSTRUCTED. The kills are driven by the *number of determined coefficients*,
  which the tie depth (up to 17) makes large: a defect-0 `e` alone removes ~14
  free coefficients, a defect-0 `sigma` ~10, and the tower needs `depth` of them.

**Consequence for Phase F.** The right global invariant is not the depth-1
hypersurface but the **defect profile** of a state: the pilot says a state dies
as soon as *enough* of its polynomials are forced defect-0 that the number of
determined sub-leading coefficients exceeds the residual freedom against the tie
depth. F1 (`phase_f_defects.py`) should be re-read through this lens — the
frontier that matters is the all-small-defect stratum, and on it the level-0 tie
tower is a genuine closer. This does **not** yet close a whole branch (states with
large defects retain cofactor freedom), so it sharpens rather than contradicts
`ALT_COMBINED.md` J1: the residue congruences DO kill, on the tight stratum, and
the next build should quantify how much of each branch is tight.

## 7. Honest / ambiguous points — [judgment]

- **[J1] The two headline kills (A depth-3, B both branches) are exact saturated
  Gröbner facts over `Q` / `Q[r]/(q)`**, reproduced by the independent verifier
  (Pilot A) which also matches the audited `alt_residue_congruences.json` support
  form. These are sound: the level-0 leading cancellation is a *necessary*
  condition on any counterexample, and it has no all-nonzero solution under the
  reconstructed divisors.
- **[J2] Pilot A kills a substratum; Pilot B kills a state.** Pilot A's `sigma`
  defect-0 is a choice (the both-tight configuration the task targets), so it
  rigorously kills the all-defect-0 configuration; Pilot C indicates the whole
  state falls too, but I mark its exact certification as pending, not proven.
- **[J3] Field of definition.** With a *marked* single root (b3000) the
  reconstructed `d1,e` have coefficients in `Q[r]/(q)`, not `Q`; the kill is
  proven over that field (and universally in `r`). For b1111 the divisor is
  Galois-stable (all four roots, `E ∝ q`), giving polynomials over `Q` — the same
  machinery applies and is a clean next pilot.
- **[J4] Only the level-0 tie is imposed.** I did not additionally impose the
  deeper `sigma = 4 d0 − d2^2` convolution coherence beyond fixing `d0`; it can
  only *shrink* the solution set, so the kills stand. Adding it is expected to
  extend Pilot C to an exact kill.
- **[J5] Not a branch kill.** No branch is closed here: large-defect states keep
  cofactor freedom (the control). The claim is confined to the tight stratum.

## 8. Verification

`python phase_f2_pilot_verify.py` (independent of `phase_f2_pilot.py`):
V0 divisor/defect-0 sanity; V1 initial form `= 2187 S^2 (4 S^3 + X^4)` and
matches `alt_residue_congruences.json` support 8; V2 exact KILL
`(j0,j1,j2,q,sat XSE) = (1)`; V3 kill universal in the marked point;
V4 soundness control (e generic ⇒ depth-3 slice solvable). **ALL CHECKS PASS.**
