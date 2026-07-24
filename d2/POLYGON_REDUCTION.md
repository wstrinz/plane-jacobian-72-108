# POLYGON_REDUCTION.md — the polygon-reduction compiler

**Files:** `polygon_reduction.py` (the compiler), `polygon_reduction_verify.py`
(the exact PASS/FAIL checker; `--quiet`, exit 0 on pass). Pure sympy, exact.

## What this module is

Every derived corner in this program (PHI_75_125, PHI_CORNER4, PHI_F14, PHI_F7,
C_SERIES_75_125, G_SYSTEM_75_125) carries the same standing judgment, copied
verbatim from `case_compiler.py` and the `PHI_*` docs:

> **[judgment] unreduced polygon** — the standard type-II.b root-shift + Laurent
> chart reduction is *assumed* (t = l, κ = l−2, deg C = a0); it is written out in
> full in no paper except the published (8,28) reduction of GGHV22.

An external review named building this compiler the single most important missing
module: it turns that *assumption* into a *derivation*, converting each
conditional corner model into a theorem about the original case. This module is
that front end. Given a GGV-chain case `(A0, A0', chain/(m,n) data)` it emits the
complete reduction — transform sequence, bracket exponent, both reduced Newton
polygons, the full branch manifest, and the corner signature the corner law
consumes.

It **builds on** the already-proved **fused-chart lemma** of
`composite_charts.py` (STEP 2): every chain transformation is either a root-shift
shear `y → y + λ x^(-s/l)` (Jacobian 1, bracket-preserving) or the one final
presentation map, and fusing them gives

```
(X, Y) = (x^-1,  x^l y + Σ_i λ_i x^(e_i)),      Jacobian  −x^(l−2).
```

We do **not** rederive that lemma; the compiler confirms the Jacobian for the
concrete `l` of each case (`fused_jacobian`) and derives `κ = l − 2` from it.

## The engine

The only genuinely computational step is the final Laurent inversion. Under
`φ(x)=x^-1, φ(y)=x^l y` a monomial `x^a y^b` maps to `x^(-a)(x^l y)^b =
x^(l·b − a) y^b`, i.e.

```
(a, b)  ⟼  (l·b − a,  b)        # invert_vertex
```

The shears only move lower/interior vertices already recorded in the
pre-inversion vertex set, so acting on the vertex set is exact. `compile_reduction`
(a) confirms the chart Jacobian is `−x^(l−2)`, sets `κ = l−2` and the bracket
`[P,Q]=x^κ`, and (b) pushes every retained pre-inversion polygon through
`invert_polygon`. The pre-inversion bracket is a nonzero constant (the flip
contributes `−1`, every root-shift contributes `+1`), so the final bracket is
`−(const)·x^(l−2)`, normalised to `x^(l−2)`.

## R1 — the published (72,108)/(8,28) reduction (the validation)

Source: GGHV22 `paper_src/2204.14178.tex`, Proposition **Case (8,28)**, lines
1000–1311; pinned in `paper_src/upstream_facts.json`.

**Transform sequence** (l = 4):

1. `φ1` flip `x ↔ y` — Jacobian −1 (line 1012).
2. `φ2` root shift `y → y + λ x^(-s)`, s ∈ {2,3} — clears the flipped lower edge;
   Jacobian 1 (lines 1073–1083).
3. `φ3` edge root shift `y → y + α x^(-4)` — reduces `{(28,8),(0,1)}` to
   `{(28,8),(24,7)}`, from the edge form `y(x^4 y − α)^7` (line 1132).
4. `φ` **final Laurent inversion** `x → x^-1, y → x^4 y` — Jacobian `−x^2`, giving
   `[φP,φQ] = −[P,Q]·x^2`, hence `[P,Q] = x^2` (line 1229).

**Bracket:** κ = l − 2 = **2**, so `[P,Q] = x^2` — derived from the Jacobian,
matching the pinned `bracket_case`.

**Reduced polygons** (computed by `(a,b) ↦ (4b−a, b)` from the pre-inversion
feet at lines 1137–1186):

| output | N(P) | N(Q) |
|---|---|---|
| **sub2** (cases a, b) | `{(0,0),(1,0),(8,14),(8,16)}` | `{(0,0),(2,1),(12,21),(12,24)}` |
| **sub1** (case c) | sub2 ∪ `{(0,8)}` | sub2 ∪ `{(0,12)}` |

Both match `upstream_facts.json` exactly (the extra `(0,8)`/`(0,12)` corner of
sub1 is the image of the intermediate `2(16,4)`/`3(16,4)` vertex of case c).

**Branch manifest** (4 branches, 12 options — the paper's own case split, each
tagged FOLLOW/EXCLUDE with a reason):

- **Pred_P(1,0) direction / shift depth** — `(1,-2)` FOLLOW, `(1,-3)` FOLLOW,
  deeper/shallower EXCLUDE (GGV6 Prop 2.5; else `deg_x P(x,0) ≤ 0` contradicts
  vd Essen Prop 10.2.6).
- **leading-form factor count** — case a) 1 factor via Pred=(2,-7), case b) 1
  factor via Pred=(1,-3), case c) 2 factors → intermediate corner `(16,4)`; 3
  distinct factors EXCLUDE (GGV2 Prop 3.12(2), a mult-6 factor would appear).
  a/b → sub2, c → sub1.
- **opposite vertex (a,b) & GGV1 Prop 8.2 exponent k** — `(24,7), k=1` FOLLOW;
  interior endpoints `{(17,5),(10,3),(3,1)}` EXCLUDE (parallel-edge closure);
  `k=2` EXCLUDE (edges can't be parallel).
- **leading forms proportional?** — both `en(P)~en(Q)` and `en(P)≁en(Q)` FOLLOW;
  they merge onto the same `k=1` conclusion (inessential branch).

**Corner signature:** t = 4, κ = 2, a0 = 8, q = 7 (from the `y(x^4y−α)^7` edge),
residual foot `(0,4)`.

## R2 — F2 j=0 = (50,75) (consistent with GGV3)

Corner `A0=(5,20)`, `A0'=(1,0)`, final corner `(7/5, 2)`, l = 5, `(m,n)=(2,3)`,
reduced pair `(a,b)=(2,3)`. Standard single-Laurent A0'=(1,0) chart: flip +
`y → y + λ x^(-2)` root shift + inversion `x → x^-1, y → x^5 y`.

- κ = l − 2 = **3**, `[P,Q] = x^3`.
- Corner signature **t=5, κ=3, a0=5, q=2, c(y)=y²(y³+1)** — matches the landed
  data exactly. The residual cubic `g = y³+1 = (y+1)(y²−y+1)` is recovered by
  solving the corner-144 forcing ODE (`_f2_forcing_divisor`).
- Φ divisor signature `(189, 75, 38, 76)`, tower length N = 36 — identical to
  `phi_corner4.py`'s landed `(50,75)` point.

GGV3 §5 discards this exact `(50,75)` case using this reduction, so for this
corner the chart is effectively published; the "unreduced polygon" judgment is
**retired at the polygon layer**.

## R3 — F2 j=1 = (75,125) (the target)

`A0=(5,20)`, `A0'=(1,0)`, final corner `(7/5,2)`, l = 5, `(m,n)=(3,5)`, reduced
pair `(a,b)=(3,5)`.

**The decisive observation:** the polygon reduction depends only on the corner
`A0=(5,20)` and `A0'=(1,0)` — *not* on `(m,n)`, which merely scales the polygon.
So F2 j=1 uses the **identical chart** as F2 j=0: same flip, same
`y → y + λ x^(-2)` shift, same inversion `(x^-1, x^5 y)`. The fused-chart lemma
forces κ = l − 2 = 3 unconditionally.

- Corner signature **t=5, κ=3, a0=5, q=2, c(y)=y²(y³+1)** — same corner ⇒ same
  residual divisor.
- Φ divisor signature `(504, 201, 101, 202)` — matches `phi_75_125.py`.

**Judgment verdict — DISCHARGED at the polygon layer.** Since the chart is the
standard single-Laurent A0'=(1,0) chart (identical to F2 j=0, which GGV3 treats),
and κ=l−2 is forced by the proved fused-chart lemma, the PHI_75_125 **judgment
item 2 ("unreduced polygon")** is discharged: **the (75,125) model is
unconditional at the polygon layer.**

**The honest boundary.** The one choice geometry alone does not pin is the
common-root gauge of the residual cubic g (unramified `g=y³+1` vs a ramified
double-root g).

> **CORRECTION (2026-07-24): this gauge choice is REOPENED, not resolved.** The
> previous text argued "deg g = 3 is odd, a cubic always has a real root, so the
> unramified gauge is realizable and selected." **That implication does not hold:
> realizability of the unramified branch does not establish branch
> COMPLETENESS.** A repeated-root (ramified) branch can *coexist* with the
> unramified one — the existence of a real root of `g` shows the unramified gauge
> is *available*, not that it is *forced* or *unique*. Our own `dg=3` work
> exhibits exactly such coexistence: multiplicities `mu=1,2,3` are simultaneously
> realized (FAMILY_GRAMMAR.md §3, F12; MU_RUNGS / ZETA_TAIL μ-graded law). So the
> residual-divisor / common-root gauge **branch completeness is REOPENED as a
> standing judgment** at the forcing/residual layer.

**What is discharged vs reopened at R3:**

- **Chart + `κ` — DISCHARGED** (unchanged): identical chart to F2 j=0, and
  `κ=l−2=3` forced by the fused-chart lemma. The polygon-layer "unreduced
  polygon" judgment is retired.
- **Selected multiplicity `q=2` — likely DISCHARGED by chain data**: fixed by the
  chain-table row `A0=(5,20) → (7/5,2)`, `k=1` (GGV5 line 1679). This is a
  chain-combinatorics selection, not a residual-gauge choice.
- **Residual-divisor / common-root gauge branch completeness — REOPENED**: the
  odd-degree real-root argument does not close it; a ramified double-root branch
  is not excluded. This is a forcing-layer judgment, not a polygon-layer flag.

## Branch-manifest sizes

| case | branches | options | followed | excluded | output shapes |
|---|---|---|---|---|---|
| (8,28) R1 | 4 | 12 | 8 | 4 | 2 (sub1, sub2) |
| F2 j=0 R2 | 3 | 6 | 4 | 2 | 1 |
| F2 j=1 R3 | 3 | 6 | 4 | 2 | 1 |

(2026-07-24: the residual-gauge branch now keeps **both** the unramified and the
ramified option OPEN — the ramified option is no longer marked EXCLUDED, since
branch completeness is reopened. This raises followed 3→4, excluded 3→2 for R2/R3.)

## Judgments retired vs surviving

| judgment | (8,28) | (50,75) | (75,125) |
|---|---|---|---|
| unreduced polygon (chart) | published (retired) | **retired** (GGV3-published corner) | **retired / discharged** (identical chart to j=0; κ=l−2 forced) |
| selected multiplicity q | audited | fixed by chain row | **likely discharged** (chain-table row `(5,20)→(7/5,2)`, k=1) |
| residual-gauge / branch completeness | n/a | **reopened** (see 2026-07-24 correction) | **REOPENED** (odd-degree real root does not force uniqueness; ramified branch not excluded) |
| forcing-polynomial identification (corner-144 correspondence) | audited | **surviving** (forcing layer, not a polygon flag) | **surviving** (forcing layer, not a polygon flag) |

The compiler retires the **polygon-layer** judgment for all three cases. What
survives for the F2 cases is the *separate* forcing-polynomial identification
(the corner-144 correspondence, audited only for (72,108)), which lives at the
forcing layer and is untouched by this front end.

## Running it

```
python3 polygon_reduction.py              # full human-readable dossiers
python3 polygon_reduction_verify.py       # 56 exact checks, verbose
python3 polygon_reduction_verify.py --quiet   # exit 0 on pass, silent
```

`polygon_reduction_verify.py` checks R1 against `paper_src/upstream_facts.json`
(with a transcribed fallback if the copyrighted `.tex`/json is absent in a public
clone), R2/R3 against the landed corner data (`phi_corner4.py`, `phi_75_125.py`),
and asserts branch-manifest completeness (every branch keeps ≥1 followed option
with a cited reason; the published (8,28) exclusions — k=2, three-factor, deep
shift — are all present; both (8,28) output shapes are retained).
