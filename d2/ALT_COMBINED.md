# Combination layer: intersecting the alternate-regime degree sweep with the finite-place lemmas

Date: 2026-07-22
Engine: `alt_combined.py` -> `alt_combined.json`
Verifier: `alt_combined_verify.py` (PASSES)
Inputs (both already audited, never edited here):
- degree layer: `alt_inf_sweep.json` (4690 surviving degree states of 38360; semantics `ALT_INF_SWEEP.md` / `ALT_REGIME_INF.md`)
- finite-place layer: `ALT_REGIME.md` (terminal + first-level parity lemmas), `ALT_REGIME_L2.md` (h6/h5 order cones section 2, residual normal form section 5)

This is the step `ALT_INF_SWEEP.md` [judgment] **J5** explicitly flagged as
"the next step ... left to the finite-place engine": the two layers were derived
independently and had **never been intersected**. This document intersects them.

## 1. Mechanism (every step cited; no new lemma)

Each surviving degree state of a branch fixes `deg d1` and `deg sigma` (and
`deg d2`, `deg e`). The finite-place lemmas constrain, at the t-place and at
every active q-place `p` (`b_i>0`), the LOCAL valuations
`x = v_P(d1)`, `z = v_P(sigma)` of any surviving counterexample to lie in an
explicit allowed cone. Valuations at distinct linear places add toward the
polynomial degree (`ALT_REGIME.md` "Orders at distinct places add toward
polynomial degree"; `ALT_REGIME_L2.md` (R1)/(R2)):

```
deg d1    >= sum_P v_P(d1)    = X,        deg sigma >= sum_P v_P(sigma) = Z.
```

**Kill rule.** A degree state dies iff NO admissible per-place selection of cone
pairs `(x_P, z_P)` satisfies `X <= deg d1` AND `Z <= deg sigma` simultaneously
(T1), resp. `Z <= deg sigma` (T2). Equivalently: the finite-place lower bounds
are incompatible with the state's enumerated degrees. This is sound (an
over-approximation of the finite-place layer): a state is killed only when even
the most favourable admissible valuation assignment cannot reach its degrees, so
a kill is a genuine proof of infeasibility.

### Cones used (transcribed verbatim; source per row)

| place | T1 allowed `(v_P(d1), v_P(sigma))` | source |
|:--|:--|:--|
| `t`, `a=11` | `5<=x<=9, 3<=z<=12` | `ALT_REGIME_L2.md` sec.2 [C] |
| `t`, `a=12` | `(3,0),(4,1),..,(8,5)`, or `x=9, 6<=z<=12` | `ALT_REGIME_L2.md` sec.2 [C] |
| `t`, `a=14` | `(6,0),(7,1),(8,2),(9,3)` | `ALT_REGIME_L2.md` sec.2 [C] |
| q `b=1` | `(1,0),(2,1)`, or `3<=x<=9, 2<=z<=12` | `ALT_REGIME_L2.md` sec.2 [C] |
| q `b=3` | `(4,0),(5,1),..,(9,5)` | `ALT_REGIME_L2.md` sec.2 [C] |
| `sigma=0` col | `t,a=11:5..9`; `t,a=12:{9}`; `t,a=14:none`; `q b=1:3..9`; `q b=3:none` | `ALT_REGIME_L2.md` sec.2 [C] |

T1 t/q minima are the degree image of the `h_7 = 8192 d1^2` top-anchor parity
(`2 v_t(d1) >= w`, odd-`s` refinement) and the q-terminal law
(`ALT_REGIME.md` "Terminal plus first-level local lemmas").

| place | T2 allowed `v_P(sigma)` | source |
|:--|:--|:--|
| `t` | `z >= w = 3a-30` (=3,6,9,12 for a=11,12,13,14) | `ALT_REGIME.md` T2 first level; `ALT_REGIME_L2.md` sec.2 (O2) |
| q `b=1` | `z >= 2` | `ALT_REGIME.md` T2 q-lemma |
| q `b=2` | impossible | `ALT_REGIME.md` T2 q-lemma |
| q `b=3` | `z = 7` (residue cancellation) | `ALT_REGIME.md` T2 q-lemma |
| q `b=4` | impossible | `ALT_REGIME.md` T2 q-lemma |

The T2 t-bound is the `sigma^2` terminal structure of the flipped chain
(`T r_5 = E^3 h_6 = -3072 E^3 sigma^2`).

## 2. Result

**State reduction census: 4690 -> 3102** (1588 states killed by the
intersection, a 33.9% further reduction of the degree-layer survivors).

**Whole-branch kills: 0.** Every one of the 27 branches retains at least one
degree state after intersection, so none is promoted to `WHOLE_BRANCH_KILL`.
This is an honest negative result, not a gap in the engine: the verifier's
PART 0 independently recomputes all 27 per-branch counts and PART 4 exhibits an
explicit survivor in the tightest branch (see [judgment] J1).

| id | a | sum_b | br | before | killed | remaining | verdict |
|:--|--:|--:|:--|--:|--:|--:|:--|
| a11_b0000_T1 | 11 | 0 | T1 | 591 | 183 | 408 | OPEN |
| a11_b1000_T1 | 11 | 1 | T1 | 487 | 186 | 301 | OPEN |
| a11_b1100_T1 | 11 | 2 | T1 | 395 | 193 | 202 | OPEN |
| a11_b1110_T1 | 11 | 3 | T1 | 319 | 186 | 133 | OPEN |
| a11_b1111_T1 | 11 | 4 | T1 | 259 | 179 | 80 | OPEN |
| a11_b3000_T1 | 11 | 3 | T1 | 319 | 213 | 106 | OPEN |
| a12_b0000_T1 | 12 | 0 | T1 | 459 | 12 | 447 | OPEN |
| a12_b1000_T1 | 12 | 1 | T1 | 370 | 45 | 325 | OPEN |
| a12_b1100_T1 | 12 | 2 | T1 | 296 | 61 | 235 | OPEN |
| a12_b1110_T1 | 12 | 3 | T1 | 238 | 74 | 164 | OPEN |
| a12_b3000_T1 | 12 | 3 | T1 | 238 | 94 | 144 | OPEN |
| a14_b0000_T1 | 14 | 0 | T1 | 227 | 14 | 213 | OPEN |
| a14_b1000_T1 | 14 | 1 | T1 | 175 | 31 | 144 | OPEN |
| a11_b0000_T2 | 11 | 0 | T2 | 30 | 3 | 27 | OPEN |
| a11_b1000_T2 | 11 | 1 | T2 | 27 | 5 | 22 | OPEN |
| a11_b1100_T2 | 11 | 2 | T2 | 24 | 7 | 17 | OPEN |
| a11_b1110_T2 | 11 | 3 | T2 | 22 | 9 | 13 | OPEN |
| a11_b1111_T2 | 11 | 4 | T2 | 20 | 11 | 9 | OPEN |
| a11_b3000_T2 | 11 | 3 | T2 | 22 | 10 | 12 | OPEN |
| a11_b3100_T2 | 11 | 4 | T2 | 20 | 12 | 8 | OPEN |
| a12_b0000_T2 | 12 | 0 | T2 | 27 | 6 | 21 | OPEN |
| a12_b1000_T2 | 12 | 1 | T2 | 24 | 8 | 16 | OPEN |
| a12_b1100_T2 | 12 | 2 | T2 | 22 | 10 | 12 | OPEN |
| a12_b1110_T2 | 12 | 3 | T2 | 20 | 12 | 8 | OPEN |
| a13_b0000_T2 | 13 | 0 | T2 | 22 | 7 | 15 | OPEN |
| a13_b1000_T2 | 13 | 1 | T2 | 20 | 9 | 11 | OPEN |
| a14_b0000_T2 | 14 | 0 | T2 | 17 | 8 | 9 | OPEN |
| **total** | | | | **4690** | **1588** | **3102** | **0 whole kills** |

### The constraint that did the most work

| kills | cited constraint |
|--:|:--|
| **1283** | `ALT_REGIME_L2.md` sec.2 [C] **T1 t-row cone** (+ `ALT_REGIME.md` T1 anchor `2 v_t(d1)>=w`, odd-`s` parity, `h_7=8192 d1^2`) |
| 188 | `ALT_REGIME_L2.md` sec.2 [C] T1 `sigma=0` column |
| 93 | `ALT_REGIME.md` T2 q-lemmas (`b=1:z>=2; b=3:z=7`) + `ALT_REGIME_L2.md` sec.2 (O2) |
| 24 | `ALT_REGIME.md` T2 t-lemma `v_t(sigma)>=w` + sec.2 (O2) |

The **T1 t-place cone did the overwhelming majority of the work (81% of all
kills)**: it forces `deg d1 >= 5` (a=11) / `>= 3` (a=12) / `>= 6` (a=14), and
for a=11 additionally `deg sigma >= 3` (the rectangle floor `z >= 3`), and the
degree sweep had kept states with `deg d1` as low as 2 and `deg sigma` as low as
0. The T2 side is small only because those branches had few surviving degree
states to begin with (17-30 each).

### Merged obligations for the 3102 survivors

Each survivor in `alt_combined.json` carries: (i) its infinity-layer
leading-cancellation obligation count (`inf_n_obligations`, from the sweep), and
(ii) a `finite_place_witness` -- an explicit admissible per-place `(x_P, z_P)`
selection realising `X <= deg d1`, `Z <= deg sigma`. Both layers are now
simultaneously satisfiable for these states; each still owes the **exact residue
congruences** `(D_t)`/`(D_p)` through `h_5` (`ALT_REGIME_L2.md` section 5) and
the closing `E^21 h_0 + u r_0 = 0`, which neither layer discharges.

## 3. Honest / ambiguous points -- [judgment]

- **[judgment] J1 -- no whole-branch kill; the intersection is a census reducer,
  not a branch closer (on these two layers).** The reason is structural and
  provable: for every open branch the forced minimum total order stays within
  the sub1 window caps. The largest T1 requirement is `deg d1 >= 9` (a11 b1111,
  a11 b3000) and `deg d1 <= 9` is the cap, so a state with `deg d1 = 9` (and
  `deg sigma >= 3`) survives; the largest T2 requirement is `sum v_P(sigma) = 12`
  (a11 b3100, a12 b1110, a14 b0000) and `deg sigma <= 12` is the cap, so
  `deg sigma = 12` survives. No branch's forced minimum EXCEEDS the cap, hence no
  whole-branch kill. A whole-branch kill would need a branch whose `X_min > 9` or
  `Z_min > 12`; those strata (e.g. b-patterns with a second `b>=2` place) were
  already removed by `ALT_REGIME_L2.md` before the 27-branch frontier and are not
  in the denominator. So on the **degree x finite-order** layers alone the
  frontier is 27 OPEN branches / 3102 states; closing a branch needs the residue
  congruences (the leading-cancellation coefficient conditions), which are a
  strictly finer fact than any valuation lower bound and are out of scope here.

- **[judgment] J2 -- I used the SHARPER `ALT_REGIME_L2.md` section 2 cones, not
  the weaker `ALT_REGIME.md` first-level bounds, wherever both exist.** E.g. for
  a=11 T1, `ALT_REGIME.md` gives only `v_t(d1) >= 3` (from `2x >= w`, `w=3`),
  whereas `ALT_REGIME_L2.md` section 2's h6/h5 cone gives `v_t(d1) >= 5`. Both
  are proven; L2 is strictly stronger and cited, so I use it. The engine's
  citation strings name both docs on the relevant rows.

- **[judgment] J3 -- the `(x,z)` coupling matters and is applied jointly.** The
  a=11 T1 t-cone is a rectangle whose `z`-floor is 3, so it kills not only small
  `deg d1` but also `deg sigma in {0,1,2}` even when `deg d1` is large (kill C in
  the verifier). The a=12/a=14 t-cones are staircases containing `(x_min, 0)`, so
  they impose no independent `sigma` floor. I compute joint feasibility by a
  reachability DP over per-place `(x,z)` sums (pruned to the `<=9`, `<=12` caps),
  not by two independent 1-D bounds, so these couplings are honoured exactly.

- **[judgment] J4 -- the residual normal form (`ALT_REGIME_L2.md` sec.5, (R0)-(R2))
  adds no kill beyond the order bounds.** Its degree content is `deg d1 >= X`,
  `deg sigma >= Z`, `X <= 9`, `Z <= 12` (used above) and `deg F <= 15-a-sum b`
  with `deg_E = sum b + deg F`. The last is already enforced by the sweep's
  enumeration (`deg_E >= sum b`, `deg e <= 15`), so it produces no new state kill;
  it is retained only as an obligation label on survivors. The `gcd(F, tq)=1`
  clause is a finite-place exactness statement with no degree consequence on the
  enumerated variables.

- **[judgment] J5 -- soundness direction.** The finite-place cones are NECESSARY
  conditions on any counterexample (proven lower bounds). A state violating them
  hosts no counterexample and is validly killed. A surviving state is NOT a
  counterexample -- it merely passes both necessary layers; the sufficient
  content (residue congruences) is untouched. So "OPEN" here means "not refuted
  by degree OR finite-order," strictly weaker than "a solution exists."

- **[judgment] J6 -- T2 b=2/b=4 impossibility is implemented but never fires**
  in the 27-branch set (no open T2 branch has a `b_i in {2,4}`); it is present so
  the engine is correct if the branch list is ever widened.

## 4. Verification (`alt_combined_verify.py`, PASSES)

Independent of `alt_combined.py` (re-transcribes the cones from the docs, never
imports the engine):

- **PART 0** recomputes all 27 per-branch `(killed, remaining)` counts from the
  sweep survivors and asserts they equal `alt_combined.json` exactly, plus the
  `4690 -> 3102` totals.
- **PART 1** hand-derives kill A: `a11_b0000_T1` state `(d2,d1,sigma,e)=(5,2,10,11)`
  dies because `deg d1 = 2 < 5 = min v_t(d1)` (a=11 T1 t-cone).
- **PART 2** hand-derives kill B: `a14_b0000_T2` state `(6,-,5,15)` dies because
  `deg sigma = 5 < 12 = w = v_t(sigma)` (a=14 T2).
- **PART 3** hand-derives kill C (the coupling): `a11_b0000_T1` `(5,6,2,11)` has
  `deg d1 = 6 >= 5` yet dies because `deg sigma = 2 < 3` = the a=11 rectangle
  `z`-floor.
- **PART 4** verifies the whole-branch verdict: 0 whole kills, every branch OPEN,
  and an explicit survivor of the tightest branch `a11_b3100_T2`
  (`deg sigma = 12`, hand `Z_min = 3+7+2 = 12`) shown feasible.

All checks PASS.
