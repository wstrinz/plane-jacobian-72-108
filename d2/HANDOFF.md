# HANDOFF: the (72,108) endgame — plane Jacobian Conjecture bound 108 → 125

> **STATUS (2026-07-22):** Frontier figures below are superseded. The current frontier is **26 sub2 cells + 171 sub1 branches + 27 alternate-regime branches** — see `cascade_cones_qt_inf_rl.json`, `cascade_cones_sub1_qt_inf_rl.json`, `alt_inf_sweep.json`, `STATE.md` (and the generated `FRONTIER.md`). The f37 status stated below is understated: the whole f37 branch is CLOSED by the ideal-membership result (`f31` lies in the pre-resultant ideal, both subcases) — see `F37_SATURATION_REPORT.md` — not merely a free-family exclusion.

> **Stale (2026-07-22).** This handoff predates the cascade-engine campaign. Current sources of truth: `STATE.md` (tail entries), `PROOF_INVENTORY.md` (claim graph and trust tiers), and the frontier artifacts it lists. The 235-strata/420-branch frontier described below is superseded (f37 closed; sub2 at 26 cells; sub1 measured and audited).

**Read this first, then STATE.md (full research log with proofs and caveats).**

## Priority-zero correction (2026-07-22)

Read `FIELD_SPLIT_AUDIT.md` before using any T5 survivor count. The scalar
label `a_q=v_q(e)` is not stable after base change: over an algebraic closure
the quartic `q` splits into four simple places, and partial-support vectors
were absent from the old 21-row ledger. Existing `(a_t,a_q)` reductions must
be read as conditional on geometric coprimality or uniform `q^r` divisibility
until they are ported to the split-place ledger. The repaired sigma-locus
theorem is field-stable; the new `a_t=7` result covers the geometrically
`q`-coprime branch. `SPLIT_PLACE_LEDGER.md` is now generated from all 327
geometric multiplicity vectors; terminal levels kill 81 of them outright.
The bare `f37` identity has the exact free family `d2=d1=0`, but the compact
pre-resultant system excludes every lift of that family; see
`F37_FREE_FAMILY_SYSTEM.md`. Correctly scoped proofs leave 235 strata and 420
live T1/T2 branch records, so the old description of a handful of remaining
cells is retired.

Post-ledger exact work kills the `a_t=9` T2 branch and the constant-cofactor
cell of `a_t=9` T1. The remaining nonconstant T1 cell satisfies
`s^6|sigma`, `s^3|R`, and has its quadratic `W` fixed before the level-3
descent. This cell is now the regression pilot for the cone-level engine in
`CASCADE_ENGINE_PLAN.md`, not the template for hundreds of hand proofs.

## Mission
The only open case below degree 125 for a plane Jacobian Conjecture
counterexample is the GGV–Horruitiner (arXiv:2204.14178, Prop 4.3) case
(8,28), left open in 2022 "for lack of computing power." We (Will + Claude,
chat sessions, July 2026) reduced it to: do four window-constrained
polynomials d2,d1,d0,dm1 exist with f(d2,d1,d0,dm1,Phi) ≡ 0 in K[y], for
f ∈ {f31, f37} (files here) and Phi = f1·C4^28 explicit? The d₋₁≡0 branch is
PROVEN impossible (STATE.md). Strong controlled numerics still support
infeasibility, but exact completion now runs through the split-place ledger
and lower cascade. The immediate job is to eliminate valuation cones in bulk
and isolate the few residue systems that genuinely require bespoke algebra.
If everything closes, the theorem is: **no counterexample below degree 125**
— the first movement of Moh's 1983 bound in 43 years beyond the 2022 step.

## Inventory
- STATE.md            — research log: derivation, proofs, caveats, audit list
- f31_deg31.txt       — genuine master factor (102 terms, weight 125)
- f37_deg37.txt       — understudy factor (618 terms, weight 134); must also
                        be shown infeasible (f31·f37·d₋₁²¹ ≡ 0 over a domain)
- jetlift.py          — VALIDATED harness (controls pass); run modes inside
- regenerate_system.py, run_singular.sh — rebuild everything from scratch
- t4_state.pkl, candidate_*.pkl, jetlift_f31_final.pkl — session artifacts

## Environment
pip: sympy scipy numpy.  apt: singular (for regeneration/certificates).
Long runs are the point — you are not time-boxed like the chat sessions were.

## Tasks, in order
T1. Implement Phase A of `CASCADE_ENGINE_PLAN.md`: extract `h_0,...,h_7`
    from `f31_graded.txt`, rewrite in `(d2,d1,sigma,e)`, emit monomial tables,
    and replay every existing exact branch proof as a regression test.
T2. Extend `split_place_ledger.py` through levels 6, 5, and 4. Eliminate
    unique local minima, retain the exact leading-coefficient equation for
    every tie, combine the four split roots by degree-budget dynamic
    programming, and report cone certificates rather than isolated rows.
T3. Continue the nonconstant `a_t=9` T1 cell through level 3 and below as
    the first residue-transition regression. Do not generalize a valuation
    rise unless its residue cancellation has an independently checked
    coefficient identity.
T4. Apply the f37 hybrid rule beyond its free family: use resultants to find
    components, then restore the smallest original equations needed to test
    whether each component lifts.
T5. Finish the upstream audit flagged in STATE.md, especially the subcase-1
    window bound, before treating any subcase-1 computation as coverage.
    Numerical sweeps remain reconnaissance, not the proof program.
T6. After the (72,108) case closes, compare a symbolic-j F_2 signature and
    the repeated `(8,28)` corner at degree 144. This is the decision gate
    between a family/template program and a degree-bound-only program.

## Rules of engagement
- No external claims, posts, or contact until T1–T6 are done; then the
  writeup conversation includes Guccione–Guccione–Horruitiner–Valqui
  (this is their program and their open question).
- Log every run; save every best-minimum pickle; never overwrite STATE.md
  history — append.
- If numerics and a derivation step ever disagree, the derivation is guilty
  until proven innocent: we lost hours twice to absolute-vs-relative
  tolerance bugs before suspecting anything real. All the numerical rules
  in jetlift.py's docstring are blood-bought; read them.
