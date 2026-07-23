# Batch convolution descent over the sub2 Phase-D residual states

**Status: PENDING AUDIT.** Same-author layer over `convolution_descent.py` and
`phase_d_states_sub2.json`; no independent audit has been performed. Every kill
below is a *candidate* kill until audited.

Runner: `batch_convolution_sub2.py` (new file; the landed driver is untouched).
Artifacts: `batch_convolution_sub2.json` (final, merged),
`batch_convolution_sub2_pass1_ungauged.json`,
`batch_convolution_sub2_gauge_raw.json`,
`batch_convolution_sub2_gauge_resume.json` (raw passes).

## 1. Worklist and dedup

`phase_d_states_sub2.json`: 220 surviving flag cases, **7888 raw degree
states**. The master identity `f31 = sum_f Phi^f e^(21-3f) h_f == 0` is
independent of the `g_zero` flag structure, so a degree-state verdict applies
to every flag case sharing the same variable zero-flags and degrees. Deduping
to unique tuples `(a_t, d1_zero(=branch==T2), sigma_zero, d2_zero, deg_d2,
deg_d1, deg_sigma, deg_e)` gives **1782 unique states** (flag/degree
consistency of the worklist checked exactly: 0 mismatches):

| triage tier | description | unique states |
|---|---|---|
| 1 | all T2 | 100 |
| 2 | T1 with deg_e = a_t (constant-E analogues) | 485 |
| 3 | T1 with sigma_zero or d2_zero | 252 |
| 4 | everything else | 945 |

## 2. Ansatz and descent window (per state)

* `e = (y+1)^a_t * (generic poly of degree m = deg_e - a_t)` — encodes
  `v_t(e) = a_t` exactly. **The q-root support conditions on e are DROPPED**:
  a sound over-approximation (strictly larger ansatz family), so a kill here
  kills the original flag case, but an UNRESOLVED here does **not** certify
  the original survives.
* `d1 = 0` for T2, else generic of degree `deg_d1`; `sigma = 0` iff
  `sigma_zero`, else generic of degree `deg_sigma`; `d2 = 0` iff `d2_zero`,
  else generic of degree `deg_d2`; `d0 = (d2^2+sigma)/4`.
* `c = -1/6630` fixed. Start degree
  `= 1 + max_f(34f + (21-3f)·deg_e + maxdeg_y h_f at the state)`
  (validated against the driver self-test shape); floor `= start - 14`.
* Forced substitution sending the **leading coefficient** of d2/d1/sigma/e's
  generic part to 0 is verdict `STATE_KILLED_BY_DEGREE_DROP` (the recorded
  degree state itself is contradicted); a leading coefficient forced to a
  nonzero value continues the chain.

## 3. Passes, budgets, and the gauge fix

| pass | mode | wall budget | attempted (triage indices) | census |
|---|---|---|---|---|
| 1 | ungauged | 45 min | 146 (0–145) | 145 UNRESOLVED, 1 SKIPPED_BUDGET |
| 2 | gauge | 45 min | 138 (0–137) | 7 CONTRADICTION, 6 DEGREE_DROP, 124 UNRESOLVED, 1 SKIPPED_BUDGET |
| 3 | gauge resume | 15 min | 56 (138–193) | 33 CONTRADICTION, 14 DEGREE_DROP, 2 FORCED, 6 UNRESOLVED, 1 SKIPPED_BUDGET |

Per-state timeout 90 s throughout (process-isolated hard kill; skips recorded,
never silent).

**Why pass 1 produced no kills.** With every coefficient of e a forceable
unknown, the overall scale of e couples all top-degree equations: the first
nonzero master coefficient always contains >= 2 unknowns, so the driver's
sound single-unknown forcing rule can never fire. All 145 completed states
stalled at their first nonzero coefficient. This is a property of the ansatz
parameterization, not evidence about the states.

**Gauge fix (passes 2–3).** The leading coefficient of e's generic part is
declared a *nonzero parameter* `gamma` instead of an unknown. Soundness:
`gamma != 0` is exactly the degree-exactness of e in the recorded state
(`deg e = deg_e`); gamma remains a fully general symbol, so no generality is
lost. Residuals whose only free symbol is gamma are necessary polynomial
constraints (every master coefficient must vanish and prior forced
substitutions are uniquely determined); the runner accumulates them and, if
their gcd — with `gamma^k` content stripped, legitimate since gamma != 0 —
becomes a nonzero constant, no admissible gamma exists: `CONTRADICTION`.

Homogeneity context (verified computationally in-session): each `h_f` is
weighted-homogeneous for `(d2,d1,d0,e) -> (l^2, l^3, l^4, l^5)` and the master
sum is homogeneous of weight 125 only if also `c -> l^17 c`. At fixed
`c = -1/6630` the scaling is **not** a symmetry — gamma is genuinely
constrained (through `gamma^17` in constant-E states), not normalizable away.
The "gauge" is therefore a parameterization choice justified by degree
exactness, not by a symmetry reduction.

## 4. Final verdict census (194 unique states attempted)

| verdict | unique states | raw states covered |
|---|---|---|
| CONTRADICTION (candidate kill, PENDING AUDIT) | 40 | — |
| STATE_KILLED_BY_DEGREE_DROP (candidate kill, PENDING AUDIT) | 20 | — |
| kills subtotal | **60** | **70** |
| UNRESOLVED | 130 | — |
| FORCED (floor reached, no kill claimed) | 2 | — |
| SKIPPED_BUDGET | 2 | — |

Raw-state coverage of all attempted states: 303 + resume ≈ 4.6% of 7888;
kill coverage 70 raw states.

## 5. Candidate kills (all PENDING AUDIT)

All 60 kills are tier-2 constant-E T1 states (deg_e = a_t). The pattern is
fully regular:

**a_t = 9, deg_e = 9, deg_d1 = 2 (cell T1 a9 b0000), 38 kills:**

| deg_d2 | deg_sigma 0–2 | deg_sigma 3–5 |
|---|---|---|
| 0, 1, 3, and d2_zero | CONTRADICTION (incompatible gamma constraints, stop 236) | DEGREE_DROP (lc sigma s3/s4/s5 forced to 0 at 237/239/241) |
| 2 | CONTRADICTION at deg_sigma 0, 2 (deg_sigma 1 = SKIPPED_BUDGET) | DEGREE_DROP s3/s4/s5 |
| 4 | DEGREE_DROP s1–s5 at deg_sigma 1–5 | (deg_sigma 0 / sigma_zero: FORCED, see section 6) |
| sigma_zero (deg_d2 0–3) | CONTRADICTION | — |

**a_t = 10, deg_e = 10, deg_d2 = 0 (cell T1 a10 b0000), 22 kills:**
CONTRADICTION for deg_d1 in {0,1,2} x deg_sigma in {0,…,6}, plus
(deg_d1, deg_sigma) = (3, 0). All by incompatible gamma constraints
(typical pair: `A1*gamma^17 + B1` and `A2*gamma^17 + B2` with
`B1/A1 != B2/A2`, gcd = 1).

Full step-by-step chains (every forced substitution, every gamma constraint,
factored residuals) are recorded per state in `batch_convolution_sub2.json`
under `kills_pending_audit[*].gauge_detail`.

## 6. Structure of the survivors (attempted but not killed)

* **All 100 T2 states: UNRESOLVED.** The first nonzero coefficient factors as
  `gamma^j * (1105*gamma^5 + 512*s_k^2)` (constant-E) or close analogues —
  a genuine two-unknown branch (`s_k^2 = -1105/512 * gamma^5`, solvable over R
  when gamma < 0). Not forceable by the driver's rules; honestly open.
* **All 24 a_t = 8 constant-E T1 states: UNRESOLVED** with residual
  `8192*b_i^2 + 9945*gamma^3*s_j^2` where `b_i`, `s_j` are the *leading*
  coefficients of d1 and sigma. Observation (not a claim): for gamma > 0 this
  is a positive-definite sum forcing both leading coefficients to 0, i.e. a
  degree drop; only the gamma < 0 branch can survive. A sign-split refinement
  of the descent would likely convert these to kills or half-kills.
* **2 FORCED states** (a_t = 9, deg_d2 = 4, deg_d1 = 2, deg_sigma in
  {sigma_zero, 0}): the descent forces exactly the known T5_90_T1 constant-E
  chain (`b2 = 3315*gamma^4/2048`, …, matching the landed self-test at
  `c = -1/6630`) down to the floor 232. The landed gate kills that chain at
  degree 226; our uniform floor (start-14) stops 6 degrees short. A deeper
  floor for these two states is expected to decide them.
* **2 SKIPPED_BUDGET:** (a9, deg_d2 2, deg_d1 2, deg_sigma 1, deg_e 9) and
  (a10, deg_d2 0, deg_d1 3, deg_sigma 1, deg_e 10) — per-state/wall timeout
  boundary artifacts, in the middle of otherwise fully-killed blocks.

## 7. Honest coverage statement

* Attempted: triage indices 0–193 of 1782 = **194 unique states** (≈ 11%):
  all of tier 1 (100 T2), and within tier 2 all a_t = 8 (24), all a_t = 9
  (41), and the first 29 of 420 a_t = 10 states.
* **Not attempted:** 391 remaining tier-2 a_t = 10 states, all 252 tier-3
  states, all 945 tier-4 states (m >= 1 generic-E states are also where the
  per-state cost grows: measured 10 s at m = 0 up to ~34 s at m = 5).
* UNRESOLVED means only that this driver's sound forcing rules stalled at the
  recorded residual (included, factored, per state); it is **not** evidence
  of survival — especially since the q-root support conditions were dropped.
* Kills are candidate kills from a same-author pipeline (worklist generator,
  engine caps, and descent runner share authorship); they require an
  independent audit (recompute the recorded coefficient chains from the
  source-linked h_f) before being counted in the proof inventory.
* Verdict reproducibility: `python batch_convolution_sub2.py` (env:
  `BATCH_GAUGE=1`, `BATCH_START_INDEX`, `BATCH_TOTAL_WALL_BUDGET`,
  `BATCH_OUT`) and `python batch_convolution_sub2.py merge`.
