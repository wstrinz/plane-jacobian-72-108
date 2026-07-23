# Cascade engine Phase B — q-place descent to level 4

**Date:** 2026-07-22
**Code:** `cascade_engine.py` (engine), `test_cascade_engine.py` (regressions),
consuming `cascade_signature.py` monomial tables and `split_place_ledger.json`.
**Output:** `cascade_cones.json` (depth 4), `cascade_cones_depth5.json` (depth 5).

## Result

Exact valuation descent at the four split places through cascade levels
7/6/5/4, coupled by the global degree budgets, on the 420 open T1/T2
branch records of the geometric ledger:

| depth (lowest level) | killed | surviving |
|---:|---:|---:|
| terminal only | 0 (ledger baseline) | 420 |
| 5 (levels 7/6/5) | 352 | 68 |
| 4 (levels 7/6/5/4) | **390** | **30** |

**Every stratum with `a_t <= 4` is killed** — the entire low-multiplicity
mass (195 open strata) that the one-stratum-at-a-time program could not
reach. Survivors live at `a_t in [5,10]` with small q-support (all have
`b_max <= 3` and at most one place with `b_i >= 2`; no surviving vector
contains `b_i = 2`).

**Status of the kills: engine-proven and independently audited.** The
engine's semantics are exact necessary conditions (see below), its
terminal layer reproduces the independently generated
`split_place_ledger.json` on all 654 records, and two kills have been
re-derived fully by hand (records below). The independent audit
(`audit_cascade_kills.py`, authored separately against the semantics
specification only, with conservative relaxed-tie semantics and its own
parser for `f31_graded.txt`) agrees on all 654 terminal records and all
420 open-branch verdicts: 390 kills confirmed, 30 survivors confirmed,
zero disagreements (`CASCADE_KILL_AUDIT.md`). The human-readable
compression of the kills into two lemma families is
`CASCADE_CONE_LEMMAS.md`.

## Semantics (soundness direction)

Ground truth is the verified cascade (`t5_multiplace_verify.py` checks 5–7):

```text
t^v g_{l+1} = ehat^3 g_l + u^l h_l,   u = c q,  v = 30-3a,
T1 terminal: ehat^3 g_7 = -u^7 (8192 d1^2)
T2 terminal: ehat^3 g_6 = -u^6 (-3072 sigma^2)
deg g_l <= 10+3a; deg d2 <= 4, deg d1 <= 6, deg sigma <= 8.
```

At a root `p_i` of `q`: `t` is a unit, `v_i(u) = 1`, `v_i(e) = b_i` exactly.
The engine enumerates exact local states `(v_i(d2), v_i(d1), v_i(sigma))`
plus global zero flags (`sigma == 0`, `d2 == 0`, `g_l == 0` identically), and
solves each level identity ultrametrically:

- unique minimum valuation on the right side ⇒ forced valuation;
- tied minimum ⇒ the valuation may rise **only** through a recorded residue
  obligation (`monomial_tie_rise`, `term_cancellation`,
  `identical_vanishing`, `exact_identity`) — never silently;
- `h_l == 0` identically is admitted only when at least two monomials
  survive the zero flags (a single nonzero monomial cannot vanish).

Local chains at one place interact with the other places **only** through
the budget sums `sum_i v_i(.) <= deg cap`, so Pareto reduction of completed
chains and the four-place DFS join are exact for existence. The engine
therefore over-approximates the true solution set: a kill means *no
consistent valuation profile exists*; a survivor comes with an explicit
witness profile and its residue obligations (Phase C input).

## Hand-audited kills

1. `(a=0, b=(2,1,1,1), T2)`, killed at level 5. Terminal identity forces
   `r_6 = 6 + 2 v_i(sigma) - 3 b_i` per place; the budget
   `sum r_6 = 9 + 2 sum v_i(sigma) <= 10` forces `v_i(sigma) = 0` for all i,
   so `r_6 = (0,3,3,3)`. At the `b=2` place the level-5 identity needs
   valuation 0 on `ehat^3 g_5 + u^5 h_5`, but both terms have valuation
   `>= 5`; the zero-flag escapes (`g_5 == 0`, `h_5 == 0`) still leave
   valuation `>= 5`. Contradiction.
2. `(a=0, b=(2,2,1,1), T1)`, killed at level 6. Terminal budget
   `sum r_7 = 10 + 2 sum v_i(d1) <= 10` forces `v_i(d1) = 0`,
   `r_7 = (1,1,4,4)`. At a `b=2` place the level-6 identity needs valuation
   1 from terms of valuation `>= 6`. Contradiction.

Both audits confirm the dominant kill pattern: **terminal budget rigidity
pins the g-valuations at high-`b` places to small values, and the next
level's identity cannot produce them** — a two-line cone lemma, not 390
separate arguments. Extracting the explicit cone inequalities that cover
all 390 kills is Phase B follow-up work; the row-level witnesses are in
`cascade_cones.json`.

## Regression suite (`test_cascade_engine.py`, in `run_tests.sh`)

1. Terminal reproduction: engine at terminal depth agrees with
   `split_place_ledger.json` feasibility on all 327 strata × both branches.
2. Hand-worked level-6 descent example: all three ultrametric cases with
   their exact obligation depths; unique-monomial levels admit no rise.
3. Depth monotonicity on a sample: adding a level never revives a kill.

## Consistency checks against the by-hand program

- The five uniform-frontier records (`(6,(1,1,1,1)) T1`, `(8,0) T1/T2`,
  `(9,0) T1`, `(10,0) T1`) all **survive** — the engine does not
  contradict the cells the by-hand program found genuinely hard (the
  `a=9 T1` cell resisted exact descent to level 2 in `T5_90_T1.md`).
- Depth-5 survivor set (68) strictly contains the depth-4 set (30).

## What the 30 survivors need (Phase C/D)

Every surviving case carries explicit obligations; none survive
obligation-free. The heavy T1 cases (up to 25 zero-flag cases at
`(8,0),(9,0),(10,0)`) are exactly the known rigid tail. Recommended order:

1. Independent audit of the 390 kills (separately written checker, ideally
   a different author/agent than the engine).
2. Cone-lemma extraction: group the 390 row kills into explicit inequality
   lemmas (the two audited patterns likely cover most).
3. Phase C residue systems for the 30 survivors, starting from the T2
   cases (4–5 zero-flag cases each, low obligation counts) and pairing the
   T1 tail with the existing `a=9` local-descent machinery and infinity
   (Phase D).
