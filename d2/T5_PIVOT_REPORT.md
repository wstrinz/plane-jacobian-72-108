# T5 pivoted lift — python-flint reimplementation report (2026-07-21)

Task: redo the pivoted t-adic lift of `t5_pivot.sing` (which timed out at
pivot slice 5 of 10 after 90 min) in python-flint with truncated-series
arithmetic, export the consistency equations, and run staged slimgb for a
per-base exact certificate. Tool: `t5_pivot_flint.py`. Base point
(unchanged, from `t5_pivot_gen.py` seed 20260721):
(a0,b0,c0,e0) = (20066, 6066, 31791, 22015) over F_32003,
grad_d0 h0 = 17048, grad_dm1 h0 = 27622 (inverses 17370 / 27284, matching the
Singular generator's header).

## 1. The lift itself: solved — 0.17 s for all 10 slices

Windows are lists of `nmod_mpoly` coefficients (index = t-power) in the 18
frees a1..a4, b1..b6, e1..e8; series products are truncated convolutions and
h_0 is evaluated Horner-in-d0 so every big multiplication has one small
factor. Per-slice truncation (slice j only needs arithmetic mod t^{j+1})
makes the pivot loop essentially free:

| slice | pivot | terms | Singular ref | time |
|---|---|---|---|---|
| 1 | c1 | 3 | 3 ✓ | 0.00s |
| 2 | c2 | 9 | 9 ✓ | 0.00s |
| 3 | c3 | 22 | 22 ✓ | 0.00s |
| 4 | c4 | 51 | 51 ✓ | 0.00s |
| 5 | c5 | 107 | 107 ✓ | 0.01s |
| 6 | c6 | 217 | — | 0.01s |
| 7 | c7 | 415 | — | 0.01s |
| 8 | c8 | 771 | — | 0.02s |
| 9 | e9 | 1379 | — | 0.04s |
| 10 | e10 | 2407 | — | 0.07s |

Every slice was re-checked to vanish after pivot insertion (asserted). The
first five term counts match Singular's `size(piv)` exactly — same
polynomials. What took Singular >90 minutes takes 0.17 s here; the Singular
bottleneck really was the repeated full `subst` in the 23-variable qring, not
the object sizes.

## 2. Structural discovery: the lift is DENSE in a weight grading

Give a_i, b_i, e_i weight i (the t-power each free sits at). The whole
construction is weight-graded: t-slice j of any window object is
weighted-homogeneous of weight j. Let P(j) = dim of the weight-j monomial
space in the 18 frees (generating function ∏ 1/(1-q^{w_i}),
weights 1..4, 1..6, 1..8). The computed pivots are EXACTLY dense: term count
of pivot j = P(j) for all j = 1..10 (3, 9, 22, 51, 107, 217, 415, 771, 1379,
2407). The consistency equations are dense up to a handful of accidental
cancellations (e.g. Q_29 has 4,812,201 terms vs P(29) = 4,812,342).

This single fact controls the whole pipeline's feasibility:

| object | weight | dense size |
|---|---|---|
| Q_11 .. Q_29 (block-1 eqs) | 11..29 | 4,086 .. 4,812,342 |
| g1 entries = H[30..40] | 30..40 | 6,535,586 .. 103,447,659 |
| block-2 eqs (E³g1 + U·H1 slices) | ≤ 59 | up to P(59) ≈ 6.84 × 10⁹ |
| terminal eqs | ≤ 240-ish | astronomically dense |

## 3. Block 1: complete (19 equations)

Full evaluation H = h_0(windows) mod t^30 (lossless for slices ≤ 29) took
226 s; slices 0..10 verified identically zero; equations are
Q_j = [t^j] h_0(windows), j = 11..29 — 19 polynomials in 18 variables,
weighted-homogeneous of weight j, total degree j, sizes 4,086 → 4,812,201
terms (17.1 M terms overall, 520 MB as Singular text).

## 4. Blocks 2–7 + terminal: stopped by the 5M term cap (by rule, and rightly)

The first full-truncation (t^41) run aborted exactly where the grading
predicts: a weight-30 power-table entry reached 6,504,802 terms
(≈ P(30) = 6,535,586 > 5,000,000 cap). g1 alone needs eleven components of
sizes P(30)..P(40) (6.5M–103M terms); block-2 equations live at weight up to
59 with ~10⁹–10¹⁰ monomials. **The swell is not an implementation artifact —
the pivoted objects genuinely fill their weight spaces**, so no re-encoding of
the same objects can compress them. Any continuation past block 1 needs a
different mathematical idea (e.g. eliminating/reducing modulo the block-1
ideal first, or working with a parametrized/implicit representation of g1),
not better arithmetic.

## 5. Gröbner attempt (Singular 4.2.1 in WSL, timeout 3600 s)

Staged driver `t5_pivot_gb.sing`:
* stage A — `std` on Q_11..Q_16 with `degBound=16` (cheap unit hunt),
* stage B — full `slimgb` on Q_11..Q_16,
* stage C — `slimgb` on all 19 (loads the 520 MB file).

A subideal containing 1 certifies the full system, so any stage can deliver
the money result.

Results:
* **Stage A**: finished in seconds; deg-capped GB of size 16; **no unit**.
* **Stage B**: stalled. slimgb's F4-style eliminations climbed the weight
  ladder (453k irreducible monomials at rank 30, 2.0M at rank 88, 2.8M at
  rank 192) and the process was **killed at ~19+ GB RSS after ~50 min**
  (WSL memory ceiling / OOM), still inside a weight-≈30 reduction.
* **Stage C**: never reached.

Per the task's contingency ("if slimgb stalls on block 1 alone, report and
stop") this is where the computation ends. But the stall is NOT the real
story — see below.

## 6. Structural verdict: 1 ∈ block-1 ideal is IMPOSSIBLE — the money result
## cannot exist before block 2

The block-1 equations Q_11..Q_29 are weighted-homogeneous of weights
11..29 > 0 (verified numerically: term counts = dense dimensions P(j) of the
weight-j spaces). For any ideal generated by homogeneous elements of positive
weight, every member has zero constant term (the weight-0 component of
p·Q_j is (weight-(−j) component of p)·Q_j = 0), hence **1 is never in the
block-1 ideal, for every base point**. Geometrically: the trivial deformation
(all 18 frees = 0, all pivots = 0, window ≡ constant base) satisfies
h_0(windows) ≡ 0 exactly, so the block-1 variety always contains the origin.

Consequences:
* The original `t5_pivot.sing` block-1 "check for 1" was structurally doomed
  from the start, independent of its performance problems.
* Stage A/B "no unit" outcomes above were forced a priori; the only
  informative block-1 questions are dim/degree of the ideal, not 1 ∈ I.
* The FIRST equation that can produce a certificate is block-2's t^0 slice:
  N2[0] = e0³·g1[0] + u(0)·h_1(base), with u(0) = −1/2 and
  h_1(base) = 8491 ≠ 0 (verified). It forces the weight-30 dense polynomial
  g1[0] = H[30] to equal a NONZERO constant — this is exactly branch (i)'s
  "β_0 = 30 exactly" condition, it breaks the grading, and it excludes the
  origin. Everything certifying lives at weight ≥ 30, i.e. beyond the term
  cap (H[30] alone is ≈ P(30) = 6.5M terms).

Recommended next steps (mathematical, not implementational):
1. Raise the cap to ~7M for the single slice H[30] (TLEN=31 run, ~6 min) and
   study the augmented system {Q_11..Q_29, H[30] − c}; note the GB will face
   the same weight-30-scale linear algebra that OOM'd stage B, so this wants
   the graded structure exploited (weight-by-weight linear algebra on
   P(w)-dimensional spaces) rather than generic slimgb.
2. Alternatively attack branch (i) through the block-1 ideal's positive-
   dimensional structure: compute dim/Hilbert data of ⟨Q_11..Q_29⟩ from a
   weight-graded elimination, and intersect with the β_0 = 30 condition
   symbolically.
3. The g-cascade in the ORIGINAL 32 window unknowns (no pivoting) keeps
   objects at degree ≤ 10 + 3(ℓ−1) (T5_NP.md bookkeeping) — the pivot
   substitution is what inflates everything to dense weight spaces. A
   Gröbner attack on the un-pivoted cascade blocks may be better conditioned
   than the pivoted one, since density here is a consequence of eliminating
   c_1..c_8, e_9, e_10.

## 7. Files

* `t5_pivot_flint.py` — the tool (lift + export + driver generation).
* `t5_pivot_eqs.txt` — all 19 block-1 equations (520 MB, chunked
  `EQ[k]=EQ[k]+(...)` lines, Singular-executable).
* `t5_pivot_eqs_small.txt` — Q_11..Q_16 only (2.5 MB) for fast GB staging.
* `t5_pivot_gb.sing` — staged slimgb driver.
* Run logs in the session scratchpad (`t5_run2.log`, `t5_gb_run.log`).
