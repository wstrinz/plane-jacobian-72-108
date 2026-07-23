# Phase F2 follow-up: do the four d2-freedom-threshold survivors die?

**Date:** 2026-07-23. **Status: ENGINE+LEMMA-PROVEN PENDING AUDIT** (same-author
layer over `phase_f2_pilot.py` / `phase_f2_scale.py` mechanism, `cascade_engine.MONOMIALS[0]`,
`alt_combined.json`, `alt_residue_congruences.json`; no independent audit).
**New files (uncommitted):** `d2_threshold.py` (independent engine + driver),
`d2_threshold_verify.py` (independent verifier, **ALL CHECKS PASS**),
this doc. READ-ONLY on every audited artifact; nothing committed.

**Headline.** Two of the four survivors (`deg d2=5` in each branch) **die** — the
level-0 tie tower is the unit ideal at depth 8 with `d2` fully free, refuting the
"free d2 rescues the tower" reading of `PHASE_F2_SCALE.md` [J2]; they were
survivors only by the census's depth-2 cost cap. The other two (`deg d2=6`) are
PENDING for compute cost. Each of the two entirely-defect-0 T2 branches thus goes
6/8 -> **7/8 killed**; **no whole branch is fully closed** (I do not claim one).

## 0. The question

`PHASE_F2_SCALE.md` killed 6/8 states in each of the two entirely-defect-0 T2
branches `a11_b3100_T2` and `a12_b1110_T2`; the four survivors are exactly the
`deg d2 in {5,6}` states of those two branches:

| state | branch | field | census verdict | why survived |
|---|---|---|---|---|
| `a11_b3100_T2` `deg d2=5` (sup12 idx3008) | T2 | `Q[r1,r2]/(q,q), r1!=r2` | NARROWED (depth 2) | census cap `md=2` (two roots) |
| `a11_b3100_T2` `deg d2=6` (sup14 idx3009) | T2 | `Q[r1,r2]/(q,q), r1!=r2` | PENDING_HEAVY | census guard `two-root + free deg-6 d2` |
| `a12_b1110_T2` `deg d2=5` (sup12 idx3065) | T2 | `Q[r]/(q)` | NARROWED (depth 2) | census cap `md=2` (single root) |
| `a12_b1110_T2` `deg d2=6` (sup14 idx3066) | T2 | `Q[r]/(q)` | PENDING_TIMEOUT (380s @ depth 2) | census cap `md=2` |

`PHASE_F2_SCALE.md` [J2] read these as a **d2-freedom threshold**: "a free d2
cofactor of degree >= 5 rescues the tower." The task: test whether that is really
freedom, or merely the census's cost-driven depth cap.

## 1. Soundness finding: d2 has NO forced structure to impose (the crux)

The task asked to impose `d2`'s **own forced** finite-place valuations
(`d2 = (y+1)^{v_t} prod (y-r_j)^{v_j} * cofactor of degree delta_d2`) and check
whether `delta_d2 < 5` closes the survivors. **It does not, because there is
nothing forced to impose:**

- The alternate-regime finite-place cones constrain **only** `v_P(d1)` (T1) and
  `v_P(sigma)` — never `v_P(d2)`. Source: `ALT_REGIME_L2.md` sec.2, which states
  verbatim "There is no additional restriction on `k` [`= v(d2)`] at this
  order-only depth: every listed pair survives ... for all `k=0,...,6`."
- `PHASE_F_DEFECTS.md`'s **alt** per-delta table (lines 61-64) lists ONLY `d1`
  and `sigma`; `d2`'s defect is not tracked for alt at all.
- The audited witnesses in `alt_combined.json` for both branches carry only a
  `sigma` split (`Zmin=12`) and the `e` b-pattern; `d2` gets no valuation entry:
  - `a11_b3100_T2`: `sigma` witness `[t:3, q(b=3):7, q(b=1):2]`;
    `e = E (y+1)^11 (y-r1)^3 (y-r2)^1`.
  - `a12_b1110_T2`: `sigma` witness `[t:6, q(b=1):2, q(b=1):2, q(b=1):2]`;
    `e = E (y+1)^12 * (q/(2048(y-r)))^1`.

So for every survivor `delta_d2 = deg d2 in {5,6}` (full defect); imposing forced
`d2`-valuations is **impossible** (none exist) — **per the task's own step 2, "if
`delta_d2` is still >= 5 the lemma fails honestly," this route reports honestly.**
The only exact algebraic constraint on `d2`, the coherence
`sigma = 4 d0 - d2^2`, is already baked into `h_0 = MONOMIALS[0]` (it is the
substitution that defines how `sigma` enters) and, at the `h_0` tie, "adds no new
relation beyond fixing `d0`" (`PHASE_F2_SCALE.md` pilot_c; `ALT_RESIDUE_CONGRUENCES.md`
K4). **The "impose d2's forced structure" lemma therefore cannot, by itself,
close these states.**

## 2. The sound kill that DOES work: push depth with d2 fully free

The survivors survive the census **only** because of cost-driven depth caps
(`phase_f2_scale.py`: `md=2` for marked-root states, `PENDING_HEAVY` for
two-root+deg-6-d2), not because free `d2` mathematically rescues the tower. The
Pilot-B precedent is decisive: `a14_b0000_T2` with a **free degree-6 `d2`** over
`Q` is NARROWED at depth 2 but the depth-9 saturated GB is the **unit ideal**
(`phase_f2_scale.json` note; `PHASE_F2_PILOT.md` sec.3). Leaving `d2` fully free
is the conservative/sound choice (extra freedom only enlarges the solution set),
so a kill under free `d2` is valid whatever `d2`'s true divisor.

`d2_threshold.py` rebuilds the level-0 tie tower independently (only the audited
`MONOMIALS[0]` shared, as the pilot does), reconstructs each survivor over its
minimal field with `d2` a **fully free** polynomial of its degree, saturates the
leading scalars `(S, E, lc d2)` nonzero (plus `r1 != r2` for two roots), and pushes
the saturated grevlex Gröbner test past depth 2 to the true tie depth 17.
The engine reproduces every census `j0` initial form exactly (cross-check).

### Per-survivor verdicts

| state | field | initial form `j0` | verdict | kill depth |
|---|---|---|---|---|
| `a12_b1110_T2` `deg d2=5` | `Q[r]/(q)` | `-2187(3E^4 - 4S^5)` | **KILLED** | **8** (5.3s) |
| `a11_b3100_T2` `deg d2=5` | `Q[r1,r2]/(q,q), r1!=r2` | `-2187(3E^4 - 4S^5)` | **KILLED** | **8** (783s) |
| `a12_b1110_T2` `deg d2=6` | `Q[r]/(q)` | irreducible 8-term `(D6,S,E)` | **PENDING (cost)** | — |
| `a11_b3100_T2` `deg d2=6` | `Q[r1,r2]/(q,q), r1!=r2` | irreducible 8-term `(D6,S,E)` | **PENDING (cost)** | — |

**Both `deg d2=5` survivors die: the depth-8 saturated ideal is the UNIT IDEAL**
after saturating `S,E,lc(d2)` nonzero, with all of `d2`'s coefficients free —
over the single-root field AND over the two-distinct-root field. The kill depth
is identical (8) in both fields; the two-root computation is only more expensive
(219s for the depth-8 GB vs 1.5s), never mathematically different. Depth trace
(two-root): not-unit at depths 1-7, unit at depth 8. The `deg d2=5` initial form
`-2187(3E^4 - 4S^5)` is `d2`-free (exactly why depth 1 does not kill and the
census, capped at `md=2`, called it NARROWED); the kill lives at depth 8. **This
refutes the "free d2 rescues the tower" reading of `PHASE_F2_SCALE.md` [J2] for
`deg d2=5`: it was a depth-cap artifact, not a real rescue.**

**Both `deg d2=6` survivors are PENDING for cost, not mathematics.** The `deg
d2=6` initial form couples all of `(D6,S,E)` and is **irreducible over Q** (so
Pilot B's factor-branch shortcut does not apply); with a free degree-6 `d2` (7
coefficients) the saturated GB is heavy even at low depth — single-root reached
only depth 2 in 655s (matching the census `PENDING_TIMEOUT`), two-root hit a
~12 GB RAM blowup before finishing depth 1 (matching `PENDING_HEAVY`). Reaching
the required depth (~8-9) is not feasible in-session. **Strong inference they die
too:** they are the marked-root analogues of Pilot B `a14_b0000_T2` (free deg-6
`d2`, `sigma`+`e` forced defect-0), which is a proven whole-state kill at depth 9
over `Q` (`phase_f2_scale.json` note; `PHASE_F2_PILOT.md` sec.3). Only the
number-field arithmetic makes them heavier; the obstruction structure is the same.
I do NOT claim these two as kills — they remain honestly PENDING.

## 3. Independent verification (`d2_threshold_verify.py`, ALL CHECKS PASS)

Re-derives the `a12_b1110_T2` `deg d2=5` kill by a chain rebuilt from scratch:
independent literal reconstruction (the complement divisor `q/(2048(y-r))` via
`sp.div`, checked exact mod `q(r)`), an independent top-window recomputation of
`h_0`'s leading coefficients (scalars `S,E` factored out of the powered factors,
`d2^k` kept as a free multivariate factor, coefficients pulled with `sympy.Poly`
and reduced mod `q(r)` — a different code path and reduction strategy from the
engine's reversed-series arrays), and an independent saturated Gröbner test.
Checks that PASS: complement exact mod `q(r)`; `deg e=15, deg sigma=12, deg d2=5`;
top coefficient nonzero; initial form `c_0 = -6561E^4 + 8748S^5` (matches census
support-12); **depth 7 NOT the unit ideal, depth 8 IS the unit ideal**; `d2`
sub-leading coefficients enter the tower as genuinely free unknowns. **ALL CHECKS
PASS**, independently confirming the depth-8 kill with `d2` fully free.

## 4. Branch-closure statement

**Each branch advances from 6/8 to 7/8 killed; NEITHER is fully closed yet.**
Per branch (8 states, `deg d2 in {none,0,1,2,3,4,5,6}`):

| branch | census killed (`deg d2 in none..4`) | + `deg d2=5` (this work) | `deg d2=6` | status |
|---|---|---|---|---|
| `a11_b3100_T2` | 6 | **+1 KILLED (depth 8)** | PENDING (cost) | **7/8 killed** |
| `a12_b1110_T2` | 6 | **+1 KILLED (depth 8)** | PENDING (cost) | **7/8 killed** |

So the two entirely-defect-0 T2 branches are now **7/8 killed each** (up from the
scale lane's 6/8), with the sole survivor in each the `deg d2=6` state — PENDING
for computational cost only. **A WHOLE-BRANCH CLOSURE is NOT achieved** and I do
not claim one: honestly, one state per branch is uncertified. The claim I stand
behind (**ENGINE+LEMMA-PROVEN PENDING AUDIT**): the `deg d2=5` states, the last
NARROWED survivors, are exact saturated-Gröbner kills at depth 8 with `d2` free.

Whole-stratum note: neither branch's stratum loses all branches — other branches
in the same `(a,sum_b)` strata remain OPEN (`ALT_COMBINED.md` census), so no
whole-stratum closure follows either.

## 5. Honest / ambiguous points — [judgment]

- **[J1] The "impose d2 forced valuations" route fails honestly — d2 is genuinely
  free in these branches.** No audited artifact forces any `v_P(d2)`
  (`ALT_REGIME_L2.md` sec.2; `PHASE_F_DEFECTS.md` alt table; `alt_combined.json`
  witnesses). `delta_d2 = deg d2 in {5,6}`. I do NOT impose any `d2`-valuation;
  doing so would be unsound.
- **[J2] The kills that DO land use `d2` FULLY FREE** — conservative and sound
  whatever `d2`'s true divisor, exactly the Pilot-B model. They are exact
  saturated-Gröbner facts (unit ideal after saturating `S,E,lc(d2) != 0`).
- **[J3] The survivors were survivors only by the census's depth cap of 2**, not
  by a real d2-freedom rescue. `PHASE_F2_SCALE.md` [J2]'s "free d2 rescues the
  tower" is refuted for `deg d2=5`: the tower is the unit ideal at depth 8. This
  sharpens, not contradicts, the scale lane — its own conservative `md=2` cap
  hid the kill.
