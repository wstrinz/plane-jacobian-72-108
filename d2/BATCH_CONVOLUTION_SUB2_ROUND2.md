# Batch convolution descent, round 2: the remaining tier-2 a_t=10 states

**Status: PENDING AUDIT.** Same-author layer over the round-1 pipeline
(`batch_convolution_sub2.py`, itself a layer over `convolution_descent.py`);
every kill below is a *candidate* kill until independently audited.

NEW files only: the landed round-1 artifacts (`batch_convolution_sub2.json`
and friends) are frozen under an auditor lane and untouched. Runner:
`batch_convolution_sub2_round2.py` (imports the committed round-1 module
unmodified; adds per-state incremental checkpointing with atomic swap — the
round-1 100-minute stall-with-no-checkpoint failure mode is closed). Merge:
`round2_merge.py`. Artifact: `batch_convolution_sub2_round2.json`
(part-2 raw pass preserved in `batch_convolution_sub2_round2_part2.json`).

## 1. Scope and where we stopped

Deterministic triage order identical to round 1 (1782 unique states). The
round-1 batch had attempted indices 0–193; the 391 remaining tier-2
`a_t = 10, deg_e = 10` constant-E T1 states occupy indices 194–584.

| part | index window | wall budget | attempted | stop |
|---|---|---|---|---|
| 1 | [194, 585) | 40 min | 204 (194–397) | wall budget exhausted at 398 |
| 2 | [398, 585) | 14 min | 58 (398–455) | wall budget exhausted at 456 |

Per-state timeout 90 s, process-isolated hard kill, gauge mode on
(`lc(e generic part) = gamma`, nonzero parameter grounded in degree-exactness),
gamma-gcd constraint accumulator — the exact mechanism that killed 22 of the
first 29 a10 states in round 1. `c = -1/6630`, floor = start - 14, q-root
support conditions dropped (sound over-approximation; kills transfer,
UNRESOLVED does not certify survival).

**Stopped exactly at triage index 456.** Unattempted after round 2:
129 tier-2 a10 states (indices 456–584, the `deg_d2 = 4` tail plus
`deg_d2 in {5,...}` blocks), all 252 tier-3, all 945 tier-4. **Tier 3 was
not started** (priority-2 field audit and the a10 continuation consumed the
remaining compute budget).

## 2. Verdict census (262 unique states attempted, all a_t=10 tier 2)

| verdict | unique states | raw states |
|---|---|---|
| CONTRADICTION (candidate kill, PENDING AUDIT) | **82** | **86** |
| UNRESOLVED | 178 | — |
| SKIPPED_BUDGET | 2 | — |

All 82 kills are `incompatible gamma constraints` with exactly two
accumulated constraints (the round-1 pattern: `A1*gamma^17 + B1` and
`A2*gamma^17 + B2` with `B1/A1 != B2/A2`, gcd 1), every one stopping at
degree 249. Full chains per state under `kills_pending_audit[*]`.

Kill structure is fully regular in `(deg_d2, deg_d1)` (7 = all deg_sigma
values in the cell):

| deg_d2 \ deg_d1 | 0 | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|---|
| 0 (round-1 remainder) | — | — | — | 5 | 7 | open | open |
| 1 | 7 | 7 | 7 | 7 | 7 | open | open |
| 2 | 7 | 7 | 7 | 7 | 7 | open | open |
| 3 | open | open | open | open | open* | open | open |
| 4 | open* | — | — | — | — | — | — |

(*the two SKIPPED_BUDGET states, `(3,4,7)` at index 397 and `(4,4,2)` at
index 455, are per-state/wall boundary artifacts inside otherwise-uniform
blocks.) With round 1, the a10 constant-E block now stands at
**104 kills / 291 attempted / 129 unattempted**.

Survivor structure: UNRESOLVED concentrates at `deg_d1 in {5,6}`,
`deg_sigma in {7,8}`, and the whole `deg_d2 in {3,4}` band; the stalling
residuals couple exactly one leading coefficient (`a3`, `b5`, `b6`, `s7`,
`s8`) with powers of gamma — one-unknown-but-gamma-coefficient shapes the
landed forcing rule soundly declines. Mean 12.3 s/state.

## 3. The a8 sign-split: CONDITIONAL observation, NOT a kill

All 24 round-1 `a_t = 8` constant-E T1 survivors stall on

```
8192*b_lc^2 + 9945*gamma^3*s_lc^2
```

with `b_lc`, `s_lc`, `gamma` the leading y-coefficients of `d1`, `sigma`,
and `e`'s generic part. A gamma-sign positivity split assumes these live in
an ordered (real-embeddable) field. **The repo does not establish that:**

* `FIELD_SPLIT_AUDIT.md` (executive finding): the ambient problem is "over
  an arbitrary characteristic-zero field `K` (and in particular over `C`)".
* The residue-lemma kills C08/C20 are justified by a confinement statement
  that is explicitly a **q-place** phenomenon. `RESIDUE_LEMMAS_DEPTH.md` §0
  (decisive structural fact): the t-place residues and jets "are **free
  rational Taylor coefficients** of the unknown polynomials `d2,d1,σ,e` ...
  They are *not* confined to the `q`-splitting field — that confinement is
  the `q`-place phenomenon behind the only two kills C08/C20"; §5:
  "Arithmetic (square-class) obstructions live only at the `q`-place, where
  residues are confined to the `q`-splitting field `Q(√17)`."
* Even where confinement holds it is arithmetic, not order-theoretic:
  `RESIDUE_LEMMAS.md` §4 kills C08/C20 because the square classes 105 / 170
  are missing from the unique quadratic subfield `Q(sqrt(17))`, and states
  outright that "the obstruction is arithmetic, not a sign/no-real-point
  obstruction"; §6: "No multi-term shape is sign-killed over the reals."

Our `b_lc, s_lc, gamma` are *global* leading coefficients of polynomials
over `K` — precisely the free-coefficient kind, not q-place residues — so
the C08/C20 confinement does not apply to them, and no other repo statement
embeds them in a real-embeddable field. Recorded verdict (also embedded as
`a8_sign_split_audit` in the JSON, `kill_credit_claimed: 0`):

* **Over an ordered `K`, gamma > 0:** positive-definite sum forces
  `b_lc = s_lc = 0`, a double degree drop — a kill *of that branch only*,
  conditional on the field.
* **Over an ordered `K`, gamma < 0:** open; exact witness
  `(gamma, s, b) = (-1, 1, sqrt(9945/8192))` (sympy-verified residual 0).
* **Over `C` (the actual ambient generality):** no kill; exact witness
  `(gamma, s, b) = (1, 1, i*sqrt(9945/8192))` (sympy-verified residual 0).

The 24 a8 states therefore remain UNRESOLVED. Honesty over yield.

## 4. Honest coverage after rounds 1+2

* Attempted (both rounds): 194 + 262 = 456 of 1782 unique states (25.6%).
  Candidate kills to date: 60 + 82 = **142 unique** (156 raw), all tier-1/2
  mechanisms, all PENDING AUDIT.
* Unattempted globally: 129 tier-2 a10 (indices 456–584), 252 tier-3,
  945 tier-4. The round-2 JSON's `unattempted_by_tier_global` field counts
  round-2 attempts only (self-contained round; see its `coverage_note`).
* UNRESOLVED is not evidence of survival (forcing-rule stall, q-root
  support dropped); each stalling residual is recorded, factored, per state.
* Same-author caveat as round 1: worklist, engine caps, and runner share
  authorship; kills need independent recomputation from the source-linked
  `h_f` before entering the proof inventory.

## 5. Reproducibility

```
BATCH2_START_INDEX=194 BATCH2_END_INDEX=585 BATCH2_WALL_BUDGET=2400 \
  python batch_convolution_sub2_round2.py     # part 1 (checkpoints per state)
BATCH2_START_INDEX=398 BATCH2_WALL_BUDGET=840 \
  BATCH2_OUT=batch_convolution_sub2_round2_part2.json \
  python batch_convolution_sub2_round2.py     # part 2
python round2_merge.py                        # merge -> round-2 JSON
```

Resuming the remaining 129 a10 states: `BATCH2_START_INDEX=456` (write to a
part-3 file and extend `round2_merge.py`'s part list).
