# Batch convolution descent over the sub1 Phase-D residual states

**Status: PENDING AUDIT.** Same-author layer over `convolution_descent.py`,
`phase_d_states_sub1.json`, and the landed sub2 batch machinery
(`batch_convolution_sub2.py`, imported and reused verbatim — not modified).
Every kill below (transferred or fresh) is a *candidate* kill until audited.

Runner: `batch_convolution_sub1.py` (new file). Artifacts:
`batch_convolution_sub1.json` (final), `batch_convolution_sub1_gauge_raw.json`
(incremental checkpoint, written after every state; the run is resumable).

## 1. Worklist and dedup census

`phase_d_states_sub1.json`: 1145 surviving flag cases, **44117 raw degree
states** (sub1 window caps d2<=6, d1<=9, sigma<=12, e<=15). Deduping to unique
tuples `(a_t, d1_zero(=branch==T2), sigma_zero, d2_zero, deg_d2, deg_d1,
deg_sigma, deg_e)` — justified exactly as in the sub2 run by the g_zero
independence of the master identity — gives **4994 unique states**:

| triage tier | description | unique states |
|---|---|---|
| 1 | all T2 | 220 |
| 2 | T1 with deg_e = a_t (constant-E) | 708 |
| 3 | T1 with sigma_zero or d2_zero | 517 |
| 4 | everything else | 3549 |

By `a` (columns = tiers 1/2/3/4):

| a | total | T2 | T1 const-E | T1 zero-flag | other |
|---|---|---|---|---|---|
| 2 | 16 | 0 | 0 | 2 | 14 |
| 3 | 74 | 0 | 0 | 4 | 70 |
| 4 | 182 | 0 | 16 | 15 | 151 |
| 5 | 304 | 1 | 0 | 29 | 274 |
| 6 | 546 | 9 | 38 | 58 | 441 |
| 7 | 738 | 32 | 0 | 86 | 620 |
| 8 | 976 | 47 | 64 | 129 | 736 |
| 9 | 1045 | 55 | 55 | 143 | 792 |
| 10 | 1113 | 76 | 535 | 51 | 451 |

(Tier-2 gaps at a in {5,7}: no T1 case at those `a` reaches deg_e = a_t
inside the sub1 caps.)

## 2. Transfer pass (free kills from sub2 round 1)

The master identity `f31 = sum_f Phi^f e^(21-3f) h_f` is
**window-independent** (same h_f tables, same Phi); only the degree caps and
stratum data differ between sub1 and sub2. A sub2 degree-state verdict
therefore applies verbatim to any sub1 state with the identical dedup tuple.

All **194** sub2-attempted unique tuples occur among the 4994 sub1 tuples
(kills and non-kills alike). Inherited verdicts, recorded in
`batch_convolution_sub1.json` under `transferred_kills` /
`transfer_verdict_census`:

| transferred verdict | unique tuples |
|---|---|
| CONTRADICTION (kill, PENDING AUDIT) | 40 |
| STATE_KILLED_BY_DEGREE_DROP (kill, PENDING AUDIT) | 20 |
| **transferred kills subtotal** | **60** (sub1 raw coverage **70**) |
| UNRESOLVED | 130 |
| FORCED | 2 |
| SKIPPED_BUDGET | 2 |

The 60 transferred kills are the sub2 round-1 tier-2 constant-E blocks:
38 at (a=9, deg_e=9, deg_d1=2) and 22 at (a=10, deg_e=10, deg_d2=0).
They inherit sub2's pending-audit status; the identical computation (same
ansatz, same h_f, same c) is the proof object, already recorded step-by-step
in `batch_convolution_sub2.json`. Matched non-kill tuples were excluded from
the fresh schedule so budget went only to genuinely new tuples.

## 3. Fresh gauge run

Same pipeline as sub2 passes 2-3: gauge mode (lc of e's generic part = nonzero
parameter gamma, grounded in the state's degree-exactness), c = -1/6630,
q-root support conditions dropped (sound over-approximation), floor = start-14,
per-state 90 s timeout in an isolated worker process, incremental checkpoint
after **every** state. Wall budget 2520 s (used 2521 s; median 11 s/state,
max 65 s).

Schedule (703 fresh targets): (a) fresh T2 at a in {9,10} — 89 states, all
attempted; (b) fresh T1 constant-E at a=10 — first 60 of 520 a-in-{9,10}
states attempted before the wall budget expired; (c) constant-E a<9 (94
states) — not reached.

### Fresh verdict census (149 unique states attempted)

| verdict | unique states |
|---|---|
| CONTRADICTION (candidate kill, PENDING AUDIT) | 68 |
| UNRESOLVED | 80 |
| SKIPPED_BUDGET | 1 |

All 68 kills are CONTRADICTIONs by exactly two incompatible gamma
constraints (a pair `A_i*gamma^17 + B_i` with `B_1/A_1 != B_2/A_2`, gcd = 1),
stopping at degree 249 of start 251. No degree-drop kills this run.

**T2 constant-E kills at a = 10 (28 kills, new territory).** In the sub2
window every T2 state went UNRESOLVED at the two-unknown branch
`gamma^j*(1105*gamma^5 + 512*s_k^2)`. In the sub1 window the T2 constant-E
(deg_e = 10 = a) states with deg_sigma <= 6 die instead: for
deg_d2 in {d2_zero, 0, 1, 2} x deg_sigma in {0..6} (7 each, 28 total) the
top two master coefficients are gamma-only and incompatible. The boundary is
exact: at deg_sigma in {7, 8} the leading sigma coefficient enters the top
coefficient (residual `... + 2718625650243010560*gamma^6*s7 + ...`) and the
state is honestly UNRESOLVED (29 such m=0 states, plus all 32 m>=1 T2 states
at a=10 and all 13 T2 states at a=9).

**T1 constant-E kills at a = 10, d2_zero cell (40 kills).** deg_d1 in {0..4}
x deg_sigma in {sigma_zero, 0..6} (8 each): same incompatible-gamma^17
mechanism (example pair: `-6561*gamma^17 - 2305843009213693952/400329564123571875`
and `-1640250*gamma^17 - 35164105890508832768/26688637608238125`). This
complements the 22 transferred a=10 kills, which sit in the deg_d2 = 0 cell.
UNRESOLVED boundary again exact: deg_sigma in {7, 8} (s7/s8 enters the top
coefficient) and deg_d1 = 5 (b5 enters); 19 such states.

**1 SKIPPED_BUDGET:** (a=10, d2_zero, deg_d1=5, deg_sigma=8, deg_e=10) — the
per-state timeout was clipped to the last 7 s of the wall budget; an artifact
of the budget boundary, not of the state.

Full step-by-step chains (every gamma constraint, factored residuals) are in
`batch_convolution_sub1.json` under `fresh_kills_pending_audit[*].steps` /
`fresh_states[*]`.

## 4. Combined verdict census

| bucket | unique tuples | raw sub1 states |
|---|---|---|
| transferred kills (PENDING AUDIT) | 60 | 70 |
| fresh kills (PENDING AUDIT) | 68 | 70 |
| **kills total** | **128** | **140** |
| transferred non-kills (UNRESOLVED/FORCED/SKIPPED) | 134 | — |
| fresh UNRESOLVED | 80 | — |
| fresh SKIPPED_BUDGET | 1 | — |
| attempted total | 343 of 4994 (6.9%) | — |

## 5. Honest coverage boundary

* **Not attempted: 4651 unique tuples** (93.1%), by tier:
  T2 a<=8: 31; T1 constant-E: 554 (446 at a=10, 14 at a=9, all 94 at
  a in {4,6,8}; the other 154 tier-2 tuples were covered by the 94
  transferred sub2 verdicts + 60 fresh attempts); T1 zero-flag: all 517;
  other T1: all 3549. Raw-state kill coverage is 140 of 44117 (0.3%) — the
  dedup leverage runs the other way in sub1 (44117 raw / 4994 unique ≈ 8.8x),
  so per-tuple kills are worth ~1.1 raw states here versus what the flag-case
  counting suggests; flag-case elimination requires killing *every* state of
  a case.
* UNRESOLVED means the driver's sound forcing rules stalled at the recorded
  residual (included, factored, per state); it is **not** evidence of
  survival — especially since the q-root support conditions were dropped.
* All kills are candidate kills from a same-author pipeline (worklist
  generator, engine caps, sub2 machinery, and this runner share authorship);
  they require an independent audit (recompute the recorded coefficient
  chains from the source-linked h_f) before entering the proof inventory.
  Transferred kills additionally inherit the sub2 round-1 audit obligation —
  auditing the 60 sub2 chains discharges both windows at once.
* Reproducibility / resume: `python batch_convolution_sub1.py` (env:
  `BATCH_TOTAL_WALL_BUDGET`, `BATCH_PER_STATE_TIMEOUT`); the run resumes
  from `batch_convolution_sub1_gauge_raw.json` (checkpoint after every
  state) and rewrites `batch_convolution_sub1.json` at the end. The obvious
  next budget: the remaining 446 a=10 constant-E states (the killing block
  was cut mid-stride by the wall budget) and the 94 constant-E states at
  a in {8,6,4}.
