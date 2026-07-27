# PASSPORT_75_125_REPAIR — the (5,20) chart exponent, the guard, and the fallout

2026-07-26. Repo `d2_plane_72_108`, task HEAD **`170b836`**. Repairs the standing
`(5,20)`-corner reduction inputs that `passport_75_125.py` (81/81, committed
`170b836`) and `PASSPORT_75_125.md` showed to be wrong, and **guards the rule
that produced them** so the next case cannot re-break it.

`passport_75_125.py` was not modified.

---

## 0. HEADLINE

> **The root cause was a dictionary applied outside its domain, and it is now a
> guarded function that raises.** The repo derived a corner's chart data from
> GGV5's **final chain corner** `A_1 = (a\l_final, b_final)` by the dictionary
> `(t, q) = (l_final, b_final)`. That dictionary is valid **exactly** on the
> **retraction shape** `b_0 = l_chart*(a_0 - 1)` with `l_chart = ceil(b_0/a_0)`.
> It holds at `(8,28)` and `(9,24)`; it **fails** at `(7,21)` and `(5,20)`.
>
> ```
>   l  5 -> 4     kappa  3 -> 2     C  y^2(y^3+1) -> y     q = ord_y C  2 -> 1
> ```
>
> `polygon_reduction.final_corner_dictionary()` now **raises** off the
> precondition. `polygon_reduction.case_f2` derives `l` instead of reading it,
> and — new — the `(5,20)` corner now has **computed** reduced polygons:
>
> ```
>   Delta' = {(0,0),(3,0),(4,1),(0,5)}
>   N(P) = 3*Delta' = {(0,0),(9,0),(12,3),(0,15)}      N(Q) = 5*Delta' = {(0,0),(15,0),(20,5),(0,25)}
>   [P,Q] = x^2      reduced degrees (15,25)
> ```
>
> And the same engine, at `(m,n) = (2,3)`, reproduces **all three** of GGV3's
> published integers for the sibling `(50,75)` at this very corner.
>
> **`weight_lemma_75_125.py`'s verdict STANDS** — re-run from scratch on the
> rebuilt G-system, under a **corrected criterion**, and *strengthened on both
> legs*. Details in §5; this was the outcome the brief said to check hardest.

---

## 1. THE ROOT CAUSE, AND WHY A CONSTANT-FIX WOULD HAVE BEEN INSUFFICIENT

`polygon_reduction.py` §0b now carries the whole story in source. The four rows
where both sides are independently known (`PASSPORT_75_125.md` P6):

| corner | `A_1` | `l_final` | `l_chart` | `b_final` | `ord C` | dictionary | `b_0 = l_chart(a_0-1)`? |
|---|---|---|---|---|---|---|---|
| `(8,28)` | `(11\4,7)` | 4 | **4** | 7 | **7** | VALID | yes, `28 = 4*7` |
| `(9,24)` | `(11\3,8)` | 3 | **3** | 8 | **8** | VALID | yes, `24 = 3*8` |
| `(7,21)` | `(11\7,2)` | 7 | **3** | 2 | **1** | BROKEN | no, `21 != 3*6` |
| `(5,20)` | `(7\5,2)` | 5 | **4** | 2 | **1** | BROKEN | no, `20 != 4*4` |

**The trap, encoded in `has_retraction`'s docstring because it is the thing that
actually went wrong.** The precondition must be tested *for the `l` that will be
used*, **not** as "does some integer `l` satisfy `b_0 = l(a_0-1)`". At `(5,20)`
some `l` does — `l = 5` gives `5*4 = 20` — and that coincidence **is** how `l = 5`
entered this repo. With the correct `l = ceil(20/5) = 4` one has `4*4 = 16 != 20`
and the test fails, which is the right answer. A repair that fixed four constants
and left the existential reading in place would have been re-broken immediately.

### The guard (`polygon_reduction.py` §0b)

| function | contract |
|---|---|
| `chart_exponent(a0, b0)` | `l = ceil(b0/a0)` — **INFERRED**, see §7 |
| `has_retraction(a0, b0, l=None)` | `b0 == l*(a0-1)` for the `l` in use |
| `final_corner_dictionary(a0, b0, l_final, b_final)` | returns `(t,q)` on the shape; **raises `FinalCornerDictionaryError`** off it |
| `corner_chart_data(a0, b0, ...)` | routes both shapes: retracted → `(l, l-2, a0, b_final)`; not → `(l, l-2, 1, 1)` and `monomial=True` |

Checked in **both directions** (`polygon_reduction_verify.py` `RG`, 21 checks):
it **returns** `(4,7)` at `(8,28)` and `(3,8)` at `(9,24)`, and **raises** at
`(7,21)` and `(5,20)`. It also reproduces GGHV22's published `(7,21)` chart data
`(t,kappa,deg C,ord C) = (3,1,1,1)` — an independent published confirmation from a
corner this lane did not otherwise touch.

Every repaired consumer now takes its inputs **through** `corner_chart_data`, so
the values cannot drift back by hand-editing a constant.

---

## 2. FILES CHANGED, AND HOW

**Root cause (1).** `polygon_reduction.py` — added §0b (the guard, ~120 lines);
`case_f2` now derives `l`, adds two new branch-manifest branches (the chart
exponent, the retraction), records the superseded `l=5` as an **explicitly
EXCLUDED option with the retraction reason attached**, and emits **computed**
pre-inversion and reduced polygons where it previously had none.
`_f2_forcing_divisor` gained a `deg g == 0` regime plus an independent
uniqueness check (general polynomial ansatz, no `f = A y^rho g^e` shape assumed).

**Named in the brief (5).** `phi_75_125.py`, `c_series_75_125.py`,
`g_system_75_125.py` (+ regenerated `g_system_75_125.json`), `f2_tower.py`,
`weight_lemma_75_125.py`.

**Found by search, same lane (2).** `window_functions_75_125.py` (rewritten; its
central object does not exist — §6), `f2_family_verify.py` (family-level closed
form).

**Verifiers (9).** `polygon_reduction_verify.py` (+ the new `RG` section),
`phi_75_125_verify.py`, `c_series_75_125_verify.py`, `g_system_75_125_verify.py`
(+ new §G), `f2_tower_verify.py`, `window_functions_75_125_verify.py`,
`phi_corner4_verify.py`, `phi_f14_verify.py`, `case_compiler_verify.py`.

**Downstream registries / control rows (7).** `case_compiler.py` (+ 3 regenerated
`case_dossier_*.json`), `family_grammar.py`, `q_window_theorem.py`,
`chain_survey_verify.py`, `phi_corner4.py`, `phi_f14.py`, `caps_audit.py` (`J6`).

**Docs (6 + this one).** Correction banners on `PHI_75_125.md`,
`C_SERIES_75_125.md`, `G_SYSTEM_75_125.md`, `F2_TOWER.md`,
`WEIGHT_LEMMA_75_125.md`, `WINDOW_FUNCTIONS_75_125.md` (which also gained the
three named refutations of §6).

### Two corner-law generalizations the repair forced

Neither is cosmetic; `(5,20)` is the first case in the repo to exercise either,
and both are shared by `phi_corner4.py`, `phi_f14.py` and `case_compiler.py`.

- **`mult_and_cofactor` — the residual-free branch `dg = a0 - q = 0`.** When `C`
  is a monomial there is no residual `g`, hence no `(y+1)` place and no residual
  cofactor: `mult = 0`, `cof = gap`. The old rule `r = a0-q-1` returns `-1` and
  `mult = e+N = 80`, `cof = -80` — nonsense. Checked discriminating: the `dg > 0`
  branch still lands `(108,144)` on `(550,205,69,276)` exactly.
- **`gap_effective` — the resonance gap.** `gap = (q-1) - a0/t` is the resonant
  degree minus the pure-ansatz degree of `f`; it is an extra unit factor only when
  it is a **positive integer** (that is `(72,108)`, `gap = 4`). At the repaired
  `(5,20)` corner `gap = -1/4`: negative and non-integral, so no resonance sits at
  or above the ansatz degree and `gap_effective = 0`. Confirmed independently —
  `phi_75_125_verify.py` §E solves the ODE with a general 12-coefficient ansatz
  and gets the unique solution `f = (1/3) y^3`, of degree exactly `rho = 3`.
  `case_compiler.law_signature` previously refused any non-integral gap; it now
  accepts a *negative* one and flags it, and still refuses a positive one.

---

## 3. RE-VERIFICATION

Every file re-verified passing in the repo (not the sandbox):

| checker | result |
|---|---|
| `polygon_reduction_verify.py` | **ALL 94** (`R1 R2 R3 RG` + branch manifest) |
| `phi_75_125_verify.py` | **ALL 64** |
| `c_series_75_125_verify.py` | **ALL 48** |
| `g_system_75_125_verify.py` | **ALL 53** |
| `f2_tower_verify.py` | **ALL 23** |
| `f2_family_verify.py` | **ALL 20** |
| `weight_lemma_75_125.py --quiet` | **45/45** |
| `window_functions_75_125_verify.py` | **ALL 46** |
| `phi_corner4_verify.py` | **ALL 43** |
| `phi_f14_verify.py` | **ALL 37** |
| `case_compiler_verify.py` | **ALL 74** |
| `family_grammar_verify.py` | **ALL 210** |
| `chain_survey_verify.py` | **ALL 26** |
| `caps_audit.py --quiet` | **70/70** |
| `zeta_tail_verify.py` | **ALL 34** |

`run_tests.sh` was **not** run (a suite was in flight for this session's whole
duration; the analysis was completed first and the edits applied only after the
in-flight run had passed the gated files).

**Anti-vacuity.** Because the registry lane found a structural change silently
making a baseline vacuous, every new negative check here is paired with a
positive control on the same machinery:

- the guard **raises** at `(5,20)`/`(7,21)` **and returns** at `(8,28)`/`(9,24)`;
- the superseded `f = -(1/9) y^5 (y^3+1)^3` is checked to **not** solve the
  repaired ODE (so "the repaired `f` solves it" is discriminating);
- `U(w)` **raises** at `(75,125)` **and returns 238** at `(72,108)`;
- the `dg = 0` branch is checked alongside `(108,144)` landing unchanged;
- the corrected syzygy criterion **detects** `(72,108)`'s K-syzygy (§5);
- the `VACUOUS` Galois verdict is paired with F9 `(56,84)` still reporting `KILLS`;
- `R3`'s "residual gauge DISSOLVED" is a string test **plus** an independent
  numeric check that `deg g = deg C - ord C = 0`;
- the "no vertical top face" test is paired with the published `(8,28)` hull,
  where two vertices *do* attain max `x`.

Two stale-string checks that the repair would otherwise have left
**misleadingly passing** were replaced rather than re-pointed:
`polygon_reduction_verify.py`'s `"REOPENED" in joined` (the word survives inside
"was REOPENED … now DISSOLVED"), and `weight_lemma`'s `A5b`/`A5c`, which asserted
`deg_slope == 14` and `lambda` non-integral — both now assert the opposite, with
reasons.

---

## 4. THE RECOMPUTED DOWNSTREAM QUANTITIES

Family conventions: `a = j+2`, `b = 2a-1`, `t = 4`, `kappa = 2`, `C = y`.
`(75,125)` is `a = 3` (`j = 1`); `(50,75)` is `a = 2` (`j = 0`).

| quantity | recorded | recomputed | verdict |
|---|---|---|---|
| `l` = `t` | 5 | **4** | **CHANGED** |
| `kappa` | 3 | **2** | **CHANGED** |
| `C`, `deg C`, `ord C` | `y^2(y^3+1)`, 5, 2 | **`y`, 1, 1** | **CHANGED** |
| `f` | `-(1/9) y^5 (y^3+1)^3` | **`(1/3) y^3`**, i.e. `A = 1/a` | **CHANGED** |
| `N_j` | `(3j+4)(5j+9)` = 98 at `j=1` | **`(3j+4)(4j+7)` = 77** | **CHANGED** — *shape survives*: `(3a-2)(5a-1) -> (3a-2)(4a-1)`, one factor tracking `t: 5 -> 4` |
| `Phi` | `-(1/9) y^201 (y^3+1)^101` | **`(1/3) y^80`**, a monomial | **CHANGED** |
| signature | `(504,201,101,202)` | **`(80,80,0,0)`** | **CHANGED** |
| `M` (forcing slice) | 36 | **29** = `12a-7` | **CHANGED** |
| `clear = aM-b` | 103 | **82** | **CHANGED** |
| `q_window` | `5j+7` = 12 at `j=1` | **`12j+17` = 29** (`= 12a-7 = M_a`) | **CHANGED** |
| `L(w) = ceil(67w/12)`, period 12 | — | **`ceil(80w/29)`, period 29** | **CHANGED** (and see §6 (R1)) |
| `deg_slope`, affine deg cap `U(w)=14w` | 14, affine | **`80/29`, NOT affine** | **REFUTED** — no affine degree cap exists |
| stripped slope `lambda` | `101/12` | **`0`** | **CHANGED**, and now degenerate |
| slice-sum lemma (`clear = aM-b`) | — | verified on **40** slices | **SURVIVES** verbatim |
| `N` formula `a[t(a+b)-(kappa+1)]-2b` | — | agrees with the built tower | **SURVIVES** |
| judgment item 3 (non-integral per-term slice index `(b-1)/a = 4/3`) | dissolved | still dissolved | **SURVIVES** |
| block rule "+5 generators / +5 spares" | — | **+4 / +4** (step `= t`) | **CHANGED** |
| generator / spare counts | 10 / 9 | **8 / 7** (`a t-kappa-2`, `(a-1)t-1`) | **CHANGED** |
| nesting `coeff_{d0^2}(G1^{a=3}) = (10/3) G1^{a=2}` | — | **still exactly `10/3`** | **SURVIVES** |
| `Phi` recurrence | `(a/(a+1)) C^{30a+3}` | **`(a/(a+1)) C^{24a+2}`** | **CHANGED** (shape survives) |
| tower-step verdict BLOCK-OBSTRUCTION | — | unchanged | **SURVIVES** |
| period coprimality `gcd(q_a, q_{a+1}) = 1` | `gcd(7,12)` | **`gcd(17,29)`** | **SURVIVES** (different numbers) |
| frac-denominator sets `{1,7} -> {1,2,3,4,6,12}` | — | **`{1,17} -> {1,29}`** | **REFUTED** — both periods PRIME |
| `q_window` THEOREM `M/gcd(M,H)` | — | exact, inputs only changed | **SURVIVES** |
| `(50,75)` `N`, signature | 36, `(189,75,38,76)` | **28, `(30,30,0,0)`** | **CHANGED** |
| Galois-transfer input `H2 = y^2-y+1` (C2, disc −3), "kills C08/C20" | — | **no residual exists** | **REFUTED / VACUOUS** |
| `POTENTIAL_PROBE` a-only prediction `(550, 69)` | DIFFERS | still DIFFERS | **SURVIVES**, *but the argument changed* — see below |
| `phi_corner4` fit `(a,b,t,kappa,a0,q)` closed forms | — | need the `dg=0` branch | **CHANGED (generalized)** |

**Two verdicts whose *reasoning* changed even though the verdict held.** Worth
flagging because a reader could otherwise carry the old argument forward:

1. The `POTENTIAL_PROBE` refutation. Under `deg C = 5`, `550` was excluded by a
   **divisibility** argument (`(550-14)/5` is not an integer). With `deg C = 1`
   every integer degree is reachable and that argument is **gone**; `550` is now
   excluded only because the tower **forces** `N = 77`. The `mult_(y+1)`
   refutation, by contrast, got *stronger*: `Phi` has no `(y+1)` place at all, so
   `mult = 0 != 69` for **every** `N`.
2. `sigma = 0`. It was "true but tautological". It is now **false as premised** —
   `deg_slope = 80/29` is not an integer, so `CAPS_AUDIT.md` §5's `deg_slope =
   504/36 = 14` is not merely definitional, it is wrong.

**UNDETERMINED.** One item, unchanged in status by this repair: whether the
actual `(75,125)` window cone dips **below** the `80/29` ray at some non-`Phi`
weight. That still needs the deeper Newton polygon of `P` (the bridge
construction), not the `u`-grading. If anything the *premise* it rests on is now
weaker: at `(72,108)` `Phi` sat strictly inside a 2-dimensional cone
(`204 < 238`), whereas here it is a monomial with no interior, so "`Phi` realises
the extreme ray" has less content — which is exactly what §6 (R3) says.

---

## 5. `weight_lemma_75_125.py` — VERDICT STANDS, and a criterion was wrong

Its conclusion — **the `(72,108)` `Phi`-divisor relation is NOT family-level; the
mechanism does not transfer to `(75,125)`** — was re-derived from scratch on the
**rebuilt 8-generator, 7-spare** G-system. It holds, and both legs got stronger.

**But the search criterion had to be replaced, and this is a finding in its own
right.** Section B asserted the graded piece has **full column rank** (nullity 0)
and read that as "no relation `c*Phi = e*B`". Nullity 0 is **sufficient but not
necessary**: a *pure inter-generator syzygy* — a nullspace vector with no
`Phi`-carrying column in its support — also lowers the rank while giving `c = 0`,
i.e. **no `Phi`-relation at all**. On the repaired system such syzygies **exist**:
the first appears at weight **43**, supported on `G1` and `G2` only. The old test
would therefore have raised a **false alarm** on correct data.

The criterion is now the right one:

> a relation `c*Phi = e*B` with `c != 0` exists **iff** the `Phi`-carrying columns
> are **dependent modulo the span of the other columns**.

Implemented as `phi_relation_exists` via a split incremental rank. The negative is
**proved over `Q`**: independence mod `p` implies independence over `Q` (clear a
rational dependency to integers of content 1; it cannot vanish mod `p`). Results:

- **`B3b`** — at every weight `29..45`, all `Phi`-carrying columns stay
  independent modulo the rest ⇒ `c = 0` forced ⇒ **no relation**. Also none for
  `e` replaced by `d0`, `d1`, `d2` (`B4`).
- **`B0b`** — the same test **detects** `(72,108)`'s K-syzygy (the `Phi` column is
  dependent at every weight). Positive control: the test is not vacuous.
- **`B3c`** — recorded explicitly that the superseded criterion **fails** on this
  system, so the replacement was necessary rather than cosmetic.

**The ord-side obstruction became total.** `q_window = 29 = M` **exactly**, so no
split `0 < w_e < M` has carry 0 — the superseded model still left the two escapes
`w_e in {12,24}`. And with `lambda = 0` the weight lemma's own interval
`[max(0, D - lambda(W-w_e)), lambda*w_e] = [80, 0]` is **EMPTY**, so the lemma now
*forbids* the relation rather than predicting it.

**Its explanation of *why* is reinforced, as `PASSPORT_75_125.md` anticipated.**
`(72,108)` is the **sporadic** case (`(8,28)`, in GGV5's "9 other pairs" table, in
no family) and sits on the **retraction** shape; `(75,125)` is the `F_2` family
member at `(5,20)` and does **not**. Different reduction shapes, so there was never
a family relation to inherit. The distinguishing invariant is unchanged in kind —
`(72,108)` has `q_window = 1` and `(75,125)` does not — only its value moved
(12 → 29).

**Nothing in the `(72,108)` program depends on any file changed here.** Checked:
the `(72,108)` control paths (`full_system_bridge.gsystem`,
`paper_src/upstream_facts.json`, the audited `Phi = (238,204,30,4)`, the
`(8,28)` polygons, `build_gsystem(2,3,4,7,204)`) are untouched and still land
exactly. `(8,28)` satisfies the retraction shape, so the guard **returns** there
and its behaviour is bit-identical. `case_dossier_GGHV_72_108.json` regenerated
**byte-identical**.

---

## 6. `window_functions_75_125.py` — its central object does not exist

Three named refutations, all machine-checked (46/46), all recorded in
`WINDOW_FUNCTIONS_75_125.md`:

- **(R1) Period 12 is refuted.** `q_window = 12a-7` (17, 29), not `5a-3` (7, 12).
  **Both are PRIME**, so the "fractional-denominator classes `{2,3,4,6,12}`" /
  "divisor lattice of the period" reading has no counterpart — the sets are just
  `{1,17}` and `{1,29}`. They remain **coprime**, so the *qualitative*
  incommensurability conclusion of `F2_TOWER.md` survives.
- **(R2) There is no affine y-degree cap.** `deg_slope = 80/29` is not an integer.
  `window_law` used to **assert** integrality — it encoded the `(72,108)` shape as
  a universal law. It now reports `deg_affine` instead, and `U(w)` raises rather
  than returning a bogus cap (while still returning 238 at `(72,108)`).
- **(R3) The two-slope cone COLLAPSES to a ray.** `C = y` a monomial ⇒ `Phi` a
  monomial ⇒ `ord_y(Phi) = deg_y(Phi)` ⇒ both slopes are `80/29` and `lambda = 0`.
  Under the extreme-ray premise the caps *pinch* (`L(w) > U_ray(w)` unless
  `29 | w`) — a demonstration that the premise does not transfer, not a window
  system.

**Direct answer to the external reviewer's "integrality resonance" story.** It was
built on the **period-7 → period-12 carry structure**. Both numbers are wrong. The
replacement periods are **17 → 29, both prime**, so the carry structure has no
divisor lattice to resonate with, and the specific `{2,3,4,6,12}` fragmentation is
refuted. What *does* survive is coarser and, arguably, cleaner: consecutive
periods are coprime, and `q_window(a) = M_a` **exactly** at every rung — which is
precisely why the ord-side carry obstruction is total at every rung rather than
just at `a = 3`.

---

## 7. PROVED / CHECKED / INFERRED

### PROVED (mathematical argument, machine-verified, not machine-dependent)

- **P1. The retraction criterion is a geometric iff.** The edge `{(0,1),(b0,a0)}`
  collapses to a vertical face under `(a,b) -> (l b - a, b)` exactly when
  `b0 = l(a0-1)`. Off it the reduced polygon has no vertical top face, so
  `deg C = 1` and `C` is a monomial. (`PASSPORT_75_125.md` Q5; re-checked here
  against the computed `Delta'`, whose max-`x` vertex is unique, versus the
  published `(8,28)` hull, where two vertices attain max `x`.)
- **P2. `kappa = l - 2` is forced, not chosen.** The fused-chart Jacobian
  `(x^-1, x^l y + shears)` is `-x^(l-2)` for *any* shears
  (`composite_charts.py`), so once `l` is fixed `kappa` is not a branch.
- **P3. The slice-sum lemma.** Every monomial of the `u^M` slice of `S^b` carries
  `c`-exponent exactly `a*M - b`, independent of the individual coefficient
  indices. It never mentions `C`, so this repair moves `N`'s value without
  touching the argument. Verified on all 40 reachable slices.
- **P4. The corrected syzygy criterion, and its negative.** A relation with
  `c != 0` exists iff the `Phi`-carrying columns are dependent modulo the rest;
  independence mod `p` implies independence over `Q`. Hence "no relation at
  weights 29..45" is proved, not sampled.
- **P5. `q_window(a) = M_a` for the whole F2 family.** `ord_y(Phi_a) =
  2(3a-1)(2a-1)` and `M_a = 12a-7`; since `4(3a-1) = M+3`, `6(2a-1) = M+1`,
  `M` is odd and `M != 0 mod 3`, one gets `gcd = 1`. Verified symbolically for
  `a = 2..12`.
- **P6. The residual-free branch.** `deg g = deg C - ord C = 0` ⇒ `g` is a monic
  constant ⇒ no `(y+1)` place and no residual cofactor ⇒ `mult = 0`, `cof = gap`.
  The 2026-07-24 "ramified vs unramified gauge" reopening is not resolved but
  **DISSOLVED**: it presupposed a cubic `g`, which presupposed `deg C = 5`, which
  presupposed a retraction this corner does not have.

### CHECKED (exact machine computation, reproducible)

- **C1. GGV3's three published integers, reproduced by the engine.**
  `1406.0886_GGV3.tex:1723-1727` — verbatim: *"`[P_1,Q_1]=x^2`, `deg(P_1)=10` and
  `deg(Q_1)=15`."* `polygon_reduction.case_f2(0)` outputs bracket `x^2` and reduced
  polygons of total degree 10 and 15. `l = 5` predicts `x^3` and `(20,30)`,
  contradicting all three. **This is the decisive external control and it is a
  reduction GGV3 *performs*, not merely tabulates.**
- **C2. The guard, both directions**, plus GGHV22's published `(7,21)` chart data
  reproduced (`(3,1,1,1)`).
- **C3. Every repaired file passing** — §3 table (~900 checks total).
- **C4. Controls untouched.** `(72,108)` `(238,204,30,4)`, `(108,144)`
  `(550,205,69,276)`, `(56,84)` `(377,107,54,216)`, `(66,231)`, `(48,64)`, the
  published `(8,28)` polygons, and the `(72,108)` polygon-route `deg` slope 14 all
  still land exactly.
- **C5. The old criterion demonstrably mis-fires** on the repaired system
  (weight 43, nullity 1, support `{G1,G2}`, no `Phi` column).
- **C6. The independent-corroboration lane.** A concurrent lane
  (`f2_branch_manifest.py`, untracked, running during this session) states in its
  own header the *same* corrected inputs — `l = t = 4, kappa = 2, q = 1, C = y` —
  derived from scratch and explicitly **not** imported from any `*_75_125.py`, with
  the same GGV3 citation. Independent agreement on the repaired data.

### INFERRED (stated as such, not proved here)

- **N1. `l = ceil(b_0/a_0)`.** THE load-bearing inference. Validated on five
  published GGHV22 reductions and **pinned at `(5,20)` by GGV3's published
  `(50,75)` reduction**, but I did not locate a single published proposition
  stating it in this form; `CORNER_RESOLVENT.md` §5.1 correctly records that no
  general dictionary exists in the literature. **Not citable as a published
  proposition.** `chart_exponent`'s docstring says so in source.
- **N2. `gap_effective`.** That a negative or non-integral resonance gap
  contributes no extra unit factor is *checked* at `(5,20)` (the general-ansatz ODE
  solve returns `deg f = rho` exactly) and is *consistent* with all five landed
  points, but it is not a proved lemma about the resonance structure in general.
- **N3. Extreme-ray premise.** That `Phi` realises the minimal-ord ray of the
  window cone is inherited verbatim from `(72,108)` and remains a premise — and
  §6 (R3) argues it has *less* content here.
- **N4. The `en`-split exclusion at `(5,20)`** (Prop 8.2(2)) is carried over from
  `PASSPORT_75_125.md` Q3, not re-derived here.

---

## 8. OPEN, FLAGGED, DELIBERATELY NOT FIXED

- **The `(7,21)` corner carries the identical defect, and it is out of this lane's
  scope.** `l_chart = ceil(21/7) = 3`, the retraction shape fails
  (`21 != 3*6`), and GGHV22 publishes `phi_3(y) = y x^3` and `[P,Q] = x`, i.e.
  `t=3, kappa=1, C=y`. So `F9`'s standing `t=7, kappa=5, q=2` is **wrong** by the
  same mechanism. `q_window_theorem.py`'s `F9 a=2 = (56,84)` row is left
  numerically unchanged (so its census figures stay reproducible) and **flagged in
  its docstring**; `case_compiler.compile_case` now attaches a
  `[SUSPECT DICTIONARY, NOT REPAIRED]` judgment to **every** family whose corner
  fails the retraction shape, naming the two published counterexamples. Repairing
  them would move landed signatures (`F9 (377,107,54,216)`, `F3`–`F6` at `(5,20)`,
  `F10`, `F11`) for which this lane has no external control. **Corners affected:
  `(5,20)` — repaired — and `(7,21)`, `(9,27)`, plus any other census row failing
  the shape.**
- **The `(5,20)` corner's *other* families.** `F3`–`F6` share `A_0 = (5,20)` and
  `A_0' = (1,0)`, so their chart exponent is also 4 and their `C` is also `y`. They
  are flagged by the same mechanism, not repaired.
- **`t = 5` has left the corner census entirely.** `phi_f14_verify.py`'s coverage
  check now reads `t in {3,4,7}`: the only two `t=5` points were `(75,125)` and
  `(50,75)`, both at `(5,20)`. Recorded there in source.
- The `(9,24)` case-2 branch gap and the `(75,125)` Belyi-emptiness result from
  `PASSPORT_75_125.md` are untouched; the emptiness verdict does not depend on any
  value repaired here (it was already computed at both `l = 4` and `l = 5`).
