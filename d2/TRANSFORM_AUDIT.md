# TRANSFORM_AUDIT.md — lost necessary conditions, pipeline-wide (2026-07-25)

**Files:** `transform_audit.py` (99 exact sympy checks; `--quiet`, exit 0 on pass).
Read-only over every existing artifact; this lane wrote only these two files.

## Why this exists

The `positive_slice` lane found that the `d3`-killing shift is **rational in `y`**,
so a solution of the shifted G-system need not reconstruct an ORIGINAL polynomial
`P` on the Newton polygon — the inverse condition was never imposed, and that one
omission is why a cell survived everything. This file asks the same question of
**every** step that changes coordinates rationally, drops an equation, or weakens
a system.

**Scope.** The `d3`-killing shift itself and the positive slices of `P = C²` are
OUT OF SCOPE (owned by `positive_slice*`); they are cited, never redone.

## Verdict up front

> **One gap, and it is a real one: the dropped `j = 4` Q-slice is NOT two free
> definitions. It collapses to a single clean necessary condition — the
> "λ row", `G4`, at u-weight 192, the missing rung of the ladder
> 156, 168, 180, **[192]**, 204 that `full_system_bridge.G_generators()` asserts.
> It is imposed nowhere, it is not an algebraic consequence of the four used
> rows, and it is worth +32 (sub2) / +48 (sub1) coefficient equations — about a
> 26 % increase in the size of the whole system.**

Everything else audited **loses nothing**, including all five controls.

---

## 1. The table

Rows are: what the step does · invertible? · condition lost · imposed elsewhere? ·
what the missing condition says.

| # | Transform / step | Where | Invertible? | Condition lost | Imposed elsewhere? | Status |
|---|---|---|---|---|---|---|
| **F1** | **Drop the `j = 4` Q-slice `(D̃³)₋₄ + λC₄²⁸ = 0`** (and its partner `(D̃²)₋₈`) | `regenerate_system.py:22`, `verify_derivation.py` §E, `T6_SELECTION_AUDIT.md` §2 | n/a (equation dropped) | **λ ∈ K is a CONSTANT**, i.e. `G4_stripped = −λ·y⁴(y+1)²⁸` | **NO** — named as "unused ammunition" in `T6_SELECTION_AUDIT.md` §3 and `STATE.md:204`, never written down in closed form, never fed to any engine | **CONDITION LOST** |
| F2 | Use `[P,Q] = x²` at **leading order only** (the f1-ODE ⇒ Φ) | `verify_derivation.py` §A | n/a | every bracket slice `n ≤ 1` | **NO** | **OPEN** (see §4) |
| F3 | Transcribe `N(Q)` and then never consume it | `paper_src/upstream_facts.json`, `window_caps_verify.py` W0 | n/a | `supp(Q)` inside its hull; corner `(2,1)`; `deg_y Q₁₂ ∈ [21,24]` | **NO** (only transcription-checked, `envelope_bounds_verify.py`) | **CONDITION AVAILABLE**, unquantified |
| **C1** | **`dm4` elimination** `sol4 = −R(S/e + d2) − d1·e/2` | `regenerate_system.py:29`, `generators.json` | **YES** on `{d₋₁ ≠ 0}` (linear, `coeff(G1,dm4) = 3·d₋₁`) | only the branch `d₋₁ = 0` | **YES** — closed outright, both legs, `AUDIT.md` §A.3, re-derived here | **NO LOSS** (mandated control) |
| **C2** | **Deep-spare eliminations** `dm5…dm16` from `(D̃²)₋ₖ = 0` | `regenerate_system.py:25-26` | YES (each slice is `2·dm_{k+4} + rest`) | the eliminated spares' own window caps | **not needed** — they are implied **term-by-term** | **NO LOSS** |
| **C3** | **Window-floor strip** `d_{4−k} = y^{12k}·d̃_{4−k}` | `full_system_bridge.py` WEIGHT | **YES**, a bijection on `{ord ≥ 12k}` | none — the rows are weighted-homogeneous, so `G(full)=0 ⟺ G(stripped)=0` | n/a | **NO LOSS** |
| **C4** | **Saturation / Rabinowitsch** `I + ⟨1 − z·f⟩` | `bridge_sweep.py`, `saturated_cell.py`, `alt_hunt_depth2.py`, … | exact encoding of `V(I) \ V(f)` | **none in the leak sense** — the auxiliary variable admits **no** solutions the original lacks | n/a | **NO LOSS** (but see §5: it *adds a hypothesis*) |
| **C5** | Bracket-slice machinery (this lane's tool) | `transform_audit.py` C5 | n/a | n/a | n/a | **CALIBRATED** — reproduces `verify_derivation.py` §A's f1-ODE exactly |
| A1 | α-strip / WLOG `Q → Q − α₂P − α₀`, `P → P + (2/3)α₁` | `T6_SELECTION_AUDIT.md` §4 | YES (affine, `[P,Q]`-preserving) | none *for the caps* (see §5.1) | — | **NO LOSS**, on the premise `ℓ(Q−C³) = α_k(x⁴C₄)^k` |
| A2 | `C₄ = y⁷(a₀+a₁y) → y⁷(y+1)` normalization | `T6_SELECTION_AUDIT.md` §4 | YES over `K̄`; needs roots over small `K` | field-of-definition, not a polynomial condition | `FIELD_SCOPE_AUDIT.md` / `FIELD_SCOPE_REPAIR.md` own this | **out of scope here** (cited) |
| A3 | `(5,20)` / general polygon reduction | `polygon_reduction.py`, `POLYGON_REDUCTION.md` | branch manifest, every option FOLLOWED/EXCLUDED with a reason | polygon-layer judgment already **retired**; residual-gauge branch completeness **REOPENED** 2026-07-24 | tracked there | **no new loss found** |
| A4 | Resultant + `factor_list` `[0]`-pick (`Ain/Bin`) | `regenerate_system.py:31-36` | resultant is a sound weakening; the pick is a case split | the `d₋₁` factor branch | **YES** — `A = d₋₁·Ah`, `B = d₋₁·Bh`, single `dm2`-factor each, `d₋₁ = 0` closed (`AUDIT.md` §A.1–3) | **NO LOSS** (legacy path; live path is the ideal-membership certificate C11) |
| A5 | Engine-level "sound over-approximations": q-root support conditions dropped; top-N `f31` coefficients | `batch_convolution_sub{1,2}.py:21`, `saturated_cell.py:37` | n/a | genuine conditions, dropped **deliberately** | **NO**, and labelled as such in-code | **CONDITION LOST by design** (sound for kills; costs proving power) |
| A6 | Positive slices of `P = C²` after the rational shift | `positive_slice.py` | — | — | — | **OTHER LANE** — cited, not duplicated |

---

## 2. F1 — the λ row, in full

### 2.1 What actually happens at `j = 4`

The P-side slices `(C²)₋ₖ = P₋ₖ = 0` hold for **every** `k ≥ 1`, and each one
*defines* the fresh spare `dm_{k+4}` (coefficient 2, exact division). Running that
chain for `k = 1..12` and substituting into the Q-side slices gives

```
Gj := (D̃³)_-j  after the P-side substitutions,    u-weight  12*(12+j):

  G1  156   used                 (D̃³)_-1 = 0
  G2  168   used                 (D̃³)_-2 = 0
  G3  180   used                 (D̃³)_-3 = 0
  G4  192   DROPPED              (D̃³)_-4 + λ·C₄²⁸ = 0
  G5  204   used                 (D̃³)_-5 + Φ      = 0
  G6  216   dropped (see F2)
  G7  228   dropped (see F2)
```

The repo's stated reason for dropping `j = 4` is that it "shares the undetermined
unknown `dm12`" with the dropped `(D̃²)₋₈` (`T6_SELECTION_AUDIT.md` §2,
`G_SYSTEM_75_125.md` §2). **That reason does not survive contact with the
arithmetic:** once `(D̃²)₋₈ = 0` is used — it *is* a true equation, `P₋₈ = 0` —
`dm12` **cancels identically** from `(D̃³)₋₄`. Checked in
`transform_audit.py` F1. What is left is

```
G4 = -(3/2) * ( 2·d0·d₋₁·d₋₃ + d0·d₋₂² + 2·d1·d₋₂·d₋₃
                + d2·d₋₃² + d₋₁²·d₋₂ − d₋₄² )
```

a five-monomial quintic in exactly the **seven G-system window variables** — no
`dm12`, no spare, nothing new.

### 2.2 The lost condition

`λ` is the coefficient `α₋₁` of the α-strip: **a scalar of `K`**, not a function
of `y`. So `G4 + λ·C₄²⁸ = 0` with `λ` constant, and since `C₄²⁸ = y¹⁹⁶(y+1)²⁸`
while `G4` is u-homogeneous of weight 192:

> **LOST CONDITION (F1).** In the stripped coordinates the pipeline already
> works in,
> ```
> G4_stripped(d̃2, d̃1, d̃0, d̃₋₁, d̃₋₂, d̃₋₃, d̃₋₄)  =  −λ · y⁴ · (y+1)²⁸ ,   λ ∈ K.
> ```
> Equivalently: the quintic `G4` has **exactly** the divisor `4·(0) + 28·(−1)`
> and no other zero, at its maximal allowed degree.

This is a *divisibility-plus-degree* condition — precisely the species of
argument `T3_WINDOW_AUDIT.md` §4 identifies as the paper's own endgame for the
cases it closes, and precisely the species the four used rows (pure vanishing
conditions) do not contain.

### 2.3 It is not implied by the used rows

Two machine-checked witnesses (`transform_audit.py` F1):

1. **Ideal non-membership.** The point
   `(d2,d1,d0,d₋₁,d₋₂,d₋₃,d₋₄) = (2Φ−3, −1/3, 1, 1, 0, 1, 1/6)`
   satisfies `G1 = G2 = G3 = G5body+Φ = 0` identically in `Φ`, and there
   `G4 = −3Φ + 37/24 ≠ 0`. Hence `G4 ∉ ⟨G1,G2,G3,G5body+Φ⟩`.
2. **Not implied up to scale either.** The 1-parameter family
   `(2Φ−3w, −1/(3w), w², 1, 0, w, 1/(6w))` lies on all four used rows for every
   `w`, and on it `G4 = −(3/2)(−w³ + 2Φw² − 1/(36w²))` — non-constant. So the
   used rows pin neither `G4` nor its divisor.

### 2.4 How much freedom it removes

Every row is homogeneous in **both** gradings, so a u-weight-`W` row has stripped
degree cap `W/6` (sub2, `deg ≤ 14k`) resp. `W/4` (sub1, `deg ≤ 15k`) — verified
monomial-by-monomial. Counting coefficient equations:

| row | u-weight | sub2 equations | sub1 equations |
|---|---:|---:|---:|
| G1 | 156 | 27 | 40 |
| G2 | 168 | 29 | 43 |
| G3 | 180 | 31 | 46 |
| G5 | 204 | 35 | 52 |
| **used total** | | **122** | **181** |
| **G4 (λ row)** | **192** | **33 − 1 = 32** | **49 − 1 = 48** |
| | | **+26 %** | **+27 %** |

(the `−1` is the free scalar `λ`). The λ row is bigger than `G1`, `G2` or `G3`
and second only to `G5`: **adding it is like adding a whole extra generator.**

### 2.5 Why this is the same failure mode as the `d3`-shift gap

Both are "the system is *sound* — dropping equations only weakens a necessary
condition — but the dropped content is exactly what would empty the survivor."
`positive_slice` found the inverse of a rational coordinate change; this is the
inverse of an *elimination*: the pipeline eliminated `λ` by dropping its defining
equation, and never re-imposed the only thing that made `λ` special — that it is
a constant.

**Recommendation:** feed `G4_stripped + λ·y⁴(y+1)²⁸ = 0` (one new scalar unknown
`λ`, 33 y-coefficient rows) into `full_system_bridge.build_state_system` as a
fifth generator, and re-run the surviving cells. It costs one unknown and buys
32/48 equations.

---

## 3. The controls (calibration evidence)

The method is only trustworthy if it says "no loss" where no loss exists.

* **C1 (mandated) — `dm4` elimination. NO LOSS.** `G1` is degree 1 in `dm4` with
  `coeff = 3·d₋₁`; `G1 = 3·d₋₁·(dm4 − T)` with `T = −R(S/e + d2) − d1·e/2`
  exactly (re-derived, then matched to `generators.json`'s `sol4`). The map
  `dm4 ↔ sol4` is a bijection on `{d₋₁ ≠ 0}`, and `d₋₁ = 0` is closed outright:
  `G1|_{d₋₁=0} = 3d₋₂d₋₃`, and both legs force `Φ = 0`. **The method returns the
  required verdict on the control.**
* **C2 — deep-spare eliminations. NO LOSS.** For every `k = 1..12` and both
  regimes, *every monomial* of the solved `dm_{k+4}` has ord-floor exactly
  `48 + 12(k+4)` and deg-ceiling exactly `56 + 14(k+4)` (sub2) / `60 + 15(k+4)`
  (sub1) — i.e. exactly that spare's own window cap. The elimination is
  **bi-graded exact**; re-imposing the eliminated spares' caps would add nothing.
* **C3 — window-floor strip. NO LOSS.** `WEIGHT` *is* the proven order floor
  `12k`; the four rows are u-homogeneous, so `G(full) = y^W·G(stripped)` and
  vanishing is **equivalent**, not merely implied. Strip/unstrip round-trips
  exactly, and `Φ_full = y²⁰⁴·(c t³⁰ q)`.
* **C4 — Rabinowitsch. NO LOSS.** `⟨ab, 1−za⟩ ∩ K[a,b] = ⟨b⟩ = ⟨ab⟩ : a^∞`, and
  every point of `V(I)` with `a ≠ 0` lifts via `z = 1/a`. The auxiliary variable
  is an exact encoding of `V(I) \ V(f)`, never larger.
* **C5 — bracket-slice machinery.** The `n = 2` slice of `[P,F] = x²`, with
  `f1 := C₄³F₋₅`, is **exactly** `8y(y+1)f1′ − 14(8y+7)f1 − y⁸(y+1)² = 0`,
  i.e. `verify_derivation.py` §A. The tool used for F2 is calibrated against the
  repo's own audited result.

---

## 4. F2 — the sub-leading bracket slices, and one retraction

`[P,Q] = [C², C³ + λC⁻¹ + F] = [C²,F]`, so the whole bracket is `[P,F] = x²`. The
`x^n` coefficient pairs `P_i` with `F₋ⱼ` at `i − j − 1 = n`. The pipeline uses
**only `n = 2`** — the single pair `(i,j) = (8,5)`, which is the f1-ODE that
produces `Φ`. Every slice `n ≤ 1` is an unused necessary condition.

**A near miss, recorded and retracted.** The `n = 1` slice is fed by `(7,5)` and
`(8,6)`. Post-shift `P̃₇ = 2C₄c̃₃ = 0`, so *if the bracket were still exactly `x²`
in the shifted coordinates* the slice would be the homogeneous ODE
`8C₄²F̃′₋₆ + 6(C₄²)′F̃₋₆ = 0`, whose kernel is `A·C₄^{−3/2} ∉ K(y)` (odd
multiplicities 7 and 1) — forcing `F̃₋₆ = 0`, and hence, via
`(C⁻¹)₋₆ = −c₂/C₄²`, the clean λ-free, Φ-free extra row

```
G6 + d2·G4 = 0            (u-weight 216)
```

**This does not hold as stated.** The `d3`-killing shift is
`(x,y) ↦ (x − s(y), y)` with Jacobian 1, so
`[P̃,Q̃](x,y) = [P,Q](x−s,y) = (x−s)²`, **not** `x²`. Its `x¹` coefficient is
`−2s = −c₃/(2C₄) ≠ 0`, so the `n = 1` slice is an *inhomogeneous* ODE that
couples `F̃₋₆` to the **unshifted `D₃`** — which is exactly the inverse-shift
variable `h = D₃` the `positive_slice` lane reintroduces. So this row is not free
content; it buys an equation and pays an unknown, and it lands squarely in the
other lane's territory. **Recorded as OPEN, not shipped as a gap.**

(Nothing already landed depends on this either way: checked that the derivation
path uses the bracket only through leading forms. F1 is unaffected — the `j = 4`
Q-slice is a *polynomiality* condition on `Q̃`, and `Q̃ = Q∘φ` is still a
polynomial in `x` of the same degree, so all `D3(j)` rows remain valid under the
shift. Likewise `F̃₋₅ = F₋₅`, so `Φ` is shift-invariant.)

---

## 5. Two things that are *not* leaks but should be said once

### 5.1 The α-strip's modification of `P` is polygon-safe

`P → P + (2/3)α₁` changes only the `(0,0)` coefficient. All three direction
maxima that drive the caps are `≥ 0` (`v₋₁,₁ ≤ 8`, `v₋₂,₁ ≤ 0`, `v₂,₋₁ ≤ 2`) and
`(0,0)` contributes `0` in each, so adding or killing the constant term can only
preserve or improve every bound. The `(0,8)` corner — the one that actually
attains the sub1 deg maximum — is untouched. No condition is lost. (The premise
that `ℓ(Q − C³)` is `α_k(x⁴C₄)^k` for `k ∈ (−2,3)` is separate, and remains
outline-only: `T6_SELECTION_AUDIT.md` §4.)

### 5.2 Saturation adds a hypothesis; it does not lose one

A UNIT saturated ideal proves emptiness **only on `{f ≠ 0}`**. That is the
opposite of the failure mode hunted here, but it is a live obligation: every
saturation factor must be discharged on the cell. The repo already separates
these correctly — `γ ≠ 0` and `lc ≠ 0` are cell/state data, while `G(r) ≠ 0` is
flagged in `saturated_cell.py` as an inherited *ansatz* assumption with
`--no-sat-Gr` to drop it. No change recommended; recorded so the distinction
stays explicit.

---

## 6. Ranking (by freedom removed)

| rank | missing condition | freedom removed | confidence |
|---|---|---|---|
| **1** | **F1 — λ-constancy: `G4_str = −λ y⁴(y+1)²⁸`** | **+32 (sub2) / +48 (sub1) coefficient equations, ~26 %** | **PROVED necessary; PROVED not implied by the used rows** |
| 2 | (out of scope) inverse of the `d3` shift / positive slices of `P = C²` | empties `a10_b0000_T1` outright | owned by `positive_slice` |
| 3 | F2 — bracket slices `n ≤ 1` | one equation per slice, but one unknown (`D₃`) too; unbounded ladder | INFERRED, mechanism verified, net gain OPEN |
| 4 | F3 — `N(Q)` support/corner conditions | unquantified; a Q-side sibling of #2 | OPEN |
| 5 | A5 — engine over-approximations (q-root support, top-N) | per-run, deliberately traded for tractability | known, in-code |

---

## 7. PROVED / CHECKED / INFERRED / OPEN inventory

**PROVED** (exact symbolic identity or explicit witness, in `transform_audit.py`)

1. `dm12` cancels identically from `(D̃³)₋₄` once `(D̃²)₋₈ = 0` is used; the
   dropped pair leaves exactly one residual condition.
2. `G4 = −(3/2)(2d0d₋₁d₋₃ + d0d₋₂² + 2d1d₋₂d₋₃ + d2d₋₃² + d₋₁²d₋₂ − d₋₄²)`,
   u-homogeneous of weight 192.
3. `G4 ∉ ⟨G1,G2,G3,G5body+Φ⟩` (explicit witness point), and `G4` is non-constant
   on a 1-parameter subfamily of the used variety.
4. Coefficient counts 27/29/31/35 (sub2) and 40/43/46/52 (sub1) for the used
   rows; 33/49 for `G4`; every row deg-homogeneous, so the caps are attained.
5. `G1 = 3d₋₁(dm4 − T)` with `T = −R(S/e + d2) − d1·e/2`; `sol4` is the unique
   root; the `d₋₁ = 0` branch forces `Φ = 0` on both legs.
6. Bi-graded exactness of the `dm5…dm16` eliminations (both regimes).
7. u-homogeneity of `G1,G2,G3,G5` ⇒ the strip is an equivalence.
8. `⟨ab, 1−za⟩ ∩ K[a,b] = ⟨b⟩`; lifting via `z = 1/a`.
9. The `n = 2` bracket slice is exactly the repo's f1-ODE.
10. The `n = 1` slice's `(8,6)` operator has kernel `A·C₄^{−3/2} ∉ K(y)`.

**CHECKED** (verified against the repo's own artifacts, not re-proved from scratch)

11. `G1,G2,G3,G5body` rebuilt here equal `generators.json` exactly.
12. `full_system_bridge.py` asserts generator weights `{156,168,180,204}` —
    192 is structurally absent from the consumed system.
13. `window_caps_verify.py` / `caps_audit.py` read only `N(P)`; `N(Q)` is
    transcription-checked and discarded.
14. λ-constancy appears in the repo only as prose "unused ammunition"
    (`T6_SELECTION_AUDIT.md` §3.1, `STATE.md:204`), never in closed form and
    never in any generator set.
15. The legacy resultant/`factor_list` path's `[0]`-pick and its `d₋₁` debris are
    already audited (`AUDIT.md` §A.1–3).

**INFERRED** (follows from premises this lane did not re-prove)

16. `λ ∈ K` is a constant — from the α-strip / GGV1 leading-form premises
    (`T6_SELECTION_AUDIT.md` §4). **F1 stands or falls with this**, and it is the
    same premise the whole `Q = C³ + λC⁻¹ + F` normal form already rests on.
17. `[P̃,Q̃] = (x−s)²` under the shift (Jacobian-1 shear, standard).
18. The α-strip's `P → P + (2/3)α₁` is polygon-safe for the caps (§5.1).

**OPEN**

19. Whether the λ row, once imposed, empties any of the 224 surviving cells.
    Not attempted here.
20. The net content of the bracket slices `n ≤ 1` after paying for the unshifted
    `D₃` (§4) — entangled with the `positive_slice` lane's inverse shift.
21. `N(Q)`'s support/corner consequences, unquantified (F3).
22. `C₄`-normalization field-of-definition (A2) — `FIELD_SCOPE_*` owns it.
23. Residual-gauge branch completeness in `polygon_reduction` — REOPENED there
    2026-07-24, not a new finding of this lane.

---

## 8. Running it

```
python -u transform_audit.py            # full report (99 checks)
python -u transform_audit.py --quiet    # exit 0 iff every check passes
```

The `dm4` control is asserted inside the script (section C1); if the method ever
flags it, the method is wrong and the run must be discarded.
