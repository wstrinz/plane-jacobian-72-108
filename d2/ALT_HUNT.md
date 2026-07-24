# ALT_HUNT — depth-2 residue closure of the 17 HUNT (BM-candidate) cells

**Date:** 2026-07-23. **Status: ALL VERDICTS PENDING AUDIT** (same-author
layer over `cascade_engine.py`, `convolution_descent.py`, `phase_f2_sub2.py`,
`phase_f_defects.json`, `phase_d_states_{sub2,sub1}.json`,
`s_unit_results.json`; no independent audit yet).
**New files (uncommitted):** `alt_hunt_depth2.py` (runner),
`alt_hunt_results.json` (full census with per-kill audit trail), this doc.
READ-ONLY on every existing artifact; nothing committed.

## 0. What this closes

`DIVISOR_LEMMAS.md` §6 (the HUNT) localised the 17 s-unit BM-candidate cells
(`s_unit_results.json` `census.candidate_rows`) to one depth-2 coefficient
each: the depth-1 top-degree S-unit relation is a single satisfiable
hypersurface (`lc(e)¹⁷` forced), and the completed kill needs the next master
coefficient(s) under the defect-0 divisor reconstruction. This lane executes
exactly that program over every **fully-forced** (all-core-defect-0) state of
the 17 cells — 49 states — and decides 45 of them.

**Headline: 45/49 states KILLED (37 at depth 2, 8 at depth 3);
13 of the 17 HUNT cells now have every fully-forced state killed — including
all 11 T2 cells.** The weight-scaling `deg σ = 8` states (the degree-250
homogeneous point that §6 flagged as the hard case, with its 254-leaf subsum
tree) die at depth 3: the subsum tree never has to be walked, because the
divisor reconstruction pins every coefficient before the S-unit structure
matters.

## 1. Mechanism (per state, everything exact)

1. **Valuation split.** The state's degrees are split into per-place
   valuations using the cascade-engine Pareto place profiles — the same
   machinery `phase_f_defects.py` used, generalised over the window
   (`ce.CONFIGS[win]`, sub2 AND sub1; `phase_f2_sub2.place_val_options`
   hardwired sub2). The enumeration of admissible splits is exhaustive.
2. **Class-polynomial reconstruction (the new step).** Group the four
   `q`-roots by joint exponent profile `(v_e, v_d1, v_σ)`. Each class `C` of
   size `n` contributes the monic degree-`n` factor `ψ_C` of `q/2048` formed
   by its roots. Parameterize every class polynomial except the largest with
   unknown monic coefficients, derive the largest as the exact quotient
   `quo(q/2048, Π ψ_i)`, and take the remainder coefficients
   `rem(q/2048, Π ψ_i) = 0` as the defining relations. Because `q` is
   squarefree, classes are automatically root-disjoint (no saturation
   needed), and the solution variety is exactly the set of root partitions
   matching the profile multiset — so **a unit ideal kills every Galois
   assignment of the split at once**. Unknowns = 4 − (largest class size)
   ≤ 3; observed: 8 kills with 0 unknowns, 19 with 1, 18 with 2.
   This uniformly subsumes what previously required separate treatments:
   marked single roots (`Q[r]/(q)`), σ-valuations on e-unmarked roots
   (`phase_f2_sub2`'s SKIPPED class), two-marked `b1100`, three-marked
   `b1110`/`b3110` (complement-root: ONE unknown), and Galois-stable `b1111`
   (`ψ = q/2048`, zero unknowns, over `Q`).
3. **Exhaustive split disjunction.** If a state's admissible split is not
   unique, it dies iff EVERY split dies (the enumeration is exhaustive;
   splits equivalent under permuting same-profile roots are deduped to one
   class representative — sound by step 2). After Galois dedup every state
   in this census had exactly ONE representative split.
4. **Master-identity walk.** The reconstructed `(d1, σ, e)` (d2 ≡ 0 in all
   17 cells) feed `convolution_descent.ConvolutionDescent`; master
   coefficients are walked from the top degree (250, or 375 for the
   `deg e = 15` cells) downward, accumulated with the class relations and
   the scalar saturation `w·Π(scalars) − 1`, and the verdict read off an
   exact grevlex Gröbner basis: **unit ideal ⇒ KILLED**.

The typical depth-2 pattern is exactly the §6 prediction: `c_top` is pure in
`E` (forcing `E¹⁷` to a fixed rational), and `c_{top−1}` is linear in the
class unknown with nonzero coefficient, forcing a rational root of the
irreducible `q` — contradiction. For the 0-unknown (Galois-stable) splits
the two coefficients are already jointly inconsistent in the scalars.

## 2. Census (49 fully-forced states across the 17 cells)

| verdict | count | note |
|---|---:|---|
| **KILLED** | **45** | 37 at depth 2, 8 at depth 3; max single-GB 117 s, total GB time ≈ 7 min |
| OPEN (cost) | 4 | all `deg d1 = 6` T1 states — the `phase_f2_sub2` [J6] grevlex blowup, unresolved at 600 s budget / 1500 s wall |

**Cells with ALL fully-forced states killed: 13/17** — all 11 T2 cells
(`sub{1,2}` × `a9_b1000`, `a8_b1100`, `a7_b1110`, plus `sub1` ×
`a10_b0000`, `a10_b3110`, `a10_b5000`, `a6_b1111`, `a6_b1111_gz5`) and the
T1 cells `sub2:a10_b0000_T1_sz1_dz1`, `sub1:a10_b0000_T1_sz1_dz1`.

The 4 open states (2 windows × {`a9_b1000_T1` deg-6-d1, `a8_b1100_T1`
deg-6-d1}) keep 4 T1 cells open. They are **cost, not mathematics**: the
same reconstruction applies; the single grevlex GB call exceeds the wall.
The campaign's named cure for exactly this shape is the msolve bridge
(`BLOWUP_DIAGNOSIS.md`) — one msolve pass over the 4 recorded systems (the
generators are in `alt_hunt_results.json`) is the obvious next dent.

## 3. Relation to prior coverage (no double-counting)

- `phase_f2_sub2.json` had previously killed exactly ONE of these 49 states
  (`sub2:a9_b1000_T1_sz1_dz1#state0`, depth 2) and SKIPPED the rest of the
  overlap (ambiguous / unmarked-root). This lane re-derives that kill
  independently (same depth) — a cross-check, plus 44 new state kills.
- The 12 `b0000`/T-only kills over `Q` whose degree tuples coincide with
  batch tuples are *reconstructed-divisor* kills, a strictly finer statement
  than the batch's degree-only kills; overlap tuples were NOT excluded here
  because the object differs (divisor-pinned vs free coefficients).
- The `s_unit_layer` census marked all 17 cells "bm_candidate_pending_
  subsums" — this lane REPLACES the pending-subsum route entirely: no BM
  bound and no subsum tree is used; the kills are ideal-membership facts.

## 4. Audit trail (for the certificate/auditor round)

Every KILLED split records in `alt_hunt_results.json`: the reconstructed
polynomials (exact, with class unknowns), the class relations
(`rem(q/2048, Πψ) = 0` coefficients), the saturation polynomial, and EVERY
accumulated master coefficient as an exact string (`gens`, with its absolute
degree). A spec-only auditor can therefore re-derive each kill without this
file's machinery: rebuild the master coefficients from `f31_graded.txt` under
the recorded ansatz, check they match `gens`, and certify the unit ideal
(Singular `lift` / the `kill_certificate_tools.py` pipeline once landed).

## 5. Honest points — [judgment]

- **[J1] Same-author verdicts.** Runner and engine share authorship with the
  cascade/convolution stack; PENDING AUDIT until a spec-only auditor
  re-derives a sample (the audit trail is deliberately complete).
- **[J2] The class-polynomial variety is exactly the split set.** Soundness
  rests on: (i) admissible-split enumeration is exhaustive (same Pareto
  machinery as the audited defect computation); (ii) monic factorizations of
  squarefree `q` correspond bijectively to root partitions; (iii) only
  leading scalars are saturated. A kill is a nonexistence proof for every
  Galois assignment of every admissible split of the state.
- **[J3] Fully-forced stratum only.** The 49 states are the all-core-defect-0
  (fully forced) states of the 17 cells — the HUNT's own scope. Higher-defect
  states in these cells are NOT touched (they were never HUNT residuals; the
  s-unit census `n_states_all0` counts exactly these 49).
- **[J4] The 4 open states are a compute wall, not a gap in the method** —
  same reconstruction, GB too big for sympy's Buchberger at the wall.
  msolve/Singular next; systems recorded.
- **[J5] What "cell closed" means.** "CLOSED(all0)" = every fully-forced
  state of that HUNT cell is killed, i.e. the cell's BM-candidate residual
  is gone. It does NOT mean the enclosing branch closes (branches carry
  other, non-fully-forced states outside the HUNT scope).

## 6. Reproduction

```
python alt_hunt_depth2.py            # full sweep (resumable checkpoint)
AH_REDO_OPEN=1 AH_BUDGET=600 AH_HARD=1500 python alt_hunt_depth2.py
                                     # long-budget retry of OPEN states
```
