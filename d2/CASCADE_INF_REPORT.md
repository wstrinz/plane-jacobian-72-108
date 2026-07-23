# The max-plus infinity layer — semantics specification (Phase D)

**Date:** 2026-07-22
**Code:** `cascade_engine.py` (`--with-inf`); regressions `test_cascade_inf.py`;
tie-equation layer `cascade_inf_ties_verify.py`.
**Purpose of this document:** a complete, code-free statement of the
infinity-place semantics, sufficient for an independently authored audit
checker (the analogue of `CASCADE_ENGINE_REPORT.md` for the finite descent).
Every rule below is a NECESSARY condition on a putative solution; the layer
over-approximates the solution set, so a kill is final (pending audit) and a
survivor carries explicit obligations.

## 1. Ground truth

All levels of the graded identity (verified in `t5_multiplace_verify.py`,
checks 5–7, from source-parsed `f31_graded.txt` tables):

```text
t^v g_1 = h_0                                  (anchor, level 0)
t^v g_{l+1} = ehat^3 g_l + u^l h_l,  l = 1..6  (u = c q, v = 30-3a)
T1 terminal:  ehat^3 g_7 = -u^7 h_7,           h_7 = 8192 d1^2
T2 terminal:  ehat^3 g_6 = -u^6 h_6|_{d1=0},   h_6|_{d1=0} = -3072 sigma^2
```

with `e = t^a ehat`, `v_t(e) = a` exact, `c = -1/6630`,
`q = 2048y^4 - 512y^3 + 320y^2 - 240y + 195` (`lc(q) = 2048`,
`q(-1) = 3315`), `t = y+1` monic. Standard regime `a <= 10` (`v >= 0`).

Degrees: `deg t^v = v`, `deg u = 4`, `deg ehat = deg e - a`,
`lc(u) = c*lc(q) = -1024/3315` (equal to Phi's top coefficient, since
`Phi = c t^30 q`).

## 2. The place at infinity

`v_inf = -deg`, so the ultrametric min-plus semantics of the finite places
dualize to max-plus over degrees:

- **Products are exact:** `deg(t^v g_{l+1}) = v + deg g_{l+1}`,
  `deg(ehat^3 g_l) = 3(deg e - a) + deg g_l`, `deg(u^l h_l) = 4l + deg h_l`.
- **Monomial degrees are exact** given per-variable degree assignments
  `(deg d2, deg d1, deg sigma, deg e)`: the `e`-slot of every `h_l`
  monomial costs `deg e` (h_l is evaluated at the original variables).
- **h_l degree:** the maximum over the surviving monomials (zero-flagged
  variables drop their monomials). A UNIQUE maximum achiever forces
  `deg h_l` to that maximum, and h_l cannot vanish identically. Several
  achievers permit `deg h_l` to DROP below the maximum only through
  leading-coefficient cancellation, floored at 0 (a nonzero polynomial has
  nonnegative degree); `h_l == 0` identically is possible only via total
  cancellation. Every granted drop is recorded as an obligation
  (`degree_tie_drop` with the tied monomials, or `identical_vanishing`).
- **Sum rule:** in `t^v g_{l+1} = ehat^3 g_l + u^l h_l`, the left degree
  equals the maximum of the two right-side degrees unless they tie; a tie
  permits the left degree to drop below the tie value only through
  cancellation of the two leading coefficients, recorded as a
  `leading_cancellation` obligation of depth = (tie value − left degree),
  with both sides' exact leading forms (the `(-1024/3315)^l` constants and
  the tied h-monomials) in the record.
- **Zero flags** are the same global statements as in the finite descent:
  `g_l == 0` turns the level identity into an exact identity between the
  remaining two terms (degrees match exactly; for `g_{l+1} == 0` the two
  leading coefficients must also cancel — recorded as `exact_identity`
  with both leading forms).
- **Terminal:** degrees match exactly:
  `3(deg e - a) + deg g_T = 4T + deg h_T` with `deg h_7 = 2 deg d1` (T1)
  resp. `deg h_6 = 2 deg sigma` (T2), single monomials.
- **Chain range:** the infinity chain runs the FULL ladder — terminal down
  to level 1 — and closes with the level-0 anchor `deg h_0 = v + deg g_1`
  (or `h_0 == 0` identically when `g_1 == 0`). Levels below the finite
  descent depth carry no outer zero flag; both zero/nonzero branches for
  `g_3, g_2, g_1` are explored and the choice is recorded in the witness.

## 3. Degrees as first-class unknowns (the sandwich)

The degree assignment is enumerated within

```text
sum_places v_p(x)  <=  deg x  <=  cap(x)      for x in {d2, d1, sigma},
a + sum_i b_i      <=  deg e  <=  e_cap,
sum_places v_p(g_l) <= deg g_l <= g_cap(l),
```

where the left-hand sums range over the four q-root places (and t when
coupled), taken from the same audited witness machinery as the finite
join. Operationally the join runs the finite budget DFS with every cap
REPLACED by the chosen degree (identically-zero variables keep their
original cap slot; the constraint is vacuous there because all finite
valuations are infinite). A branch/flag case survives only if SOME degree
assignment admits BOTH a consistent infinity chain and a finite-place
selection fitting under it. Caps: sub2 `(d2,d1,sigma,e) <= (4,6,8,10)`,
`deg g_l <= 10+3a`; sub1 per `sub1_cascade_verify.py` with the per-level
min(forward, backward) g-caps.

## 4. Stage 2: exact tie equations

The depth-one equation of an infinity degree tie is coefficient extraction
at the TOP: writing `d2 = y^K(D + D_1 y^{-1} + ...)` etc., the tied-sum
leading coefficient must vanish — the SAME polynomial in `(D,X,S,E)` as
the residue-lemma initial form (IF) of `RESIDUE_LEMMAS.md` section 1,
with leading coefficients in place of local residues. Consequences,
verified in `cascade_inf_ties_verify.py`:

- the 21 CONSTRAINT lemmas (C01–C07, C09–C19, C21–C23) pin infinity ties
  to the identical hypersurfaces (backbone P6/P10/P11 checked exactly);
- the two arithmetic KILLS C08 (level 5) and C20 (level 4) forbid the
  corresponding infinity drops outright under `--residue-kills`: at
  infinity the unknown leading coefficients lie in the BASE field, and
  C08/C20 have no torus points over Q (a fortiori none over the
  q-splitting field).

## 5. Soundness notes

- The layer only ever ADDS necessary conditions to the audited q/t join:
  survivors with `--with-inf` are a subset of the corresponding q+t
  survivors (regression R4), and a branch with no finite witness is dead
  without consulting degrees.
- A depth-0 tie (left degree equals the tie value) is generic and carries
  no obligation, exactly as in the finite descent.
- A `FORCED`/narrowed outcome is a reduction, not a kill; kills require
  either an empty option set at some level, an unsatisfiable terminal or
  anchor, or an empty degree sandwich.
- Everything is derived from the source-parsed `h_f` tables; no
  coefficient is hand-copied in the engine or the checkers.

## 6. Regression ladder (all PASS, wired into run_tests.sh)

| rung | content | source of truth |
|:--|:--|:--|
| R0 | max-plus semantics units | hand-worked level-6 example |
| R1 | a_t=9 T2: deg e=9 dead all states; linear E open with obligations | T5_90_T2.md §2 |
| R2 | exactly 43/50 constant-cell kills; survivors d=2 × {σ≡0, z≤5} | T5_90_T1.md §3 |
| R3 | all seven T2-column margins dead; pattern-A/B states open | T5_T2_COLUMN.md §2–6 |
| R4 | joint q+t+inf ⊆ audited qt survivors; automatic G5 narrowing (10,8) | cascade_cones_qt.json |
| A–D | tie-equation constants, hypersurface identification, forbidden drops, label reconstruction | source q, RESIDUE_LEMMAS.md |
