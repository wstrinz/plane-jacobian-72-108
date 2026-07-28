# Where GGV3's `−10` comes from: the primitivity depth is a Newton-polygon vertex

> ### SUPERSEDED IN PART, 2026-07-28 (later the same day) — read this first
>
> **§1's "one maximally degenerate corner" and §6's "`N(P₁)` is known at only two
> corners" are RETIRED.** `N(P₁)` is known at **five**, one of them in print:
> `passport_75_125.PUB["7_21"]` already held GGHV22's own published reduced
> polygon (`2204.14178.tex:1313-1320`) and was unnoticed, and `(9,36)` and the
> off-class `(7,42)` were derived by an `A_0'`-recovery lemma. See
> `SECOND_CORNER.md` and `second_corner_probe.py` (67/67).
>
> **The eight candidate formulas collapse to one.** With `j* = a₀`, the depth must
> be divisible by `a₀`; `(7,21)` alone kills four, and the three formal survivors
> are *provably identical* wherever the law is non-vacuous. So the derivation
> below stands, and its entire remaining content is the single γ-dependent
> integer **δ**.
>
> **δ is still not determined, and is deliberately not fitted.** The two published
> charts (`γ=3, δ=2` and `γ=2, δ=3`) both sit at `(5,20)`, where
> `a₀ = l+1 = ρ = 5` makes four candidate rules coincide; and the recovered γ at
> the new corners is **4** and **6**, so neither has a published δ to borrow.

**Status: 2026-07-28. Partial result — the mechanism is derived and controlled;
the depth law has one test, and the loose end is named in §5.**

---

## 0. The one integer the 125 programme is stuck on

`ENDPOINT_CONTRACT.md` §2's kill predicate is **identical** on both sides of the
(50,75) comparison. Only its input differs:

| input | `required_nonzero` | `forced_floor` | predicate |
|---|---|---|---|
| GGV3 γ=3 contract | `[(-1,3), (-2,4), (0,-10)]` | `{0:-6, -1:3, -2:4}` | **fires at `(0,-10)`** |
| our transferred class-row data | `[]` | `{0:9}` | nothing |

Two of GGV3's three required-nonzeros are already **derived in-repo**:
`MOH_CONTROL_50_75.md` §2 Steps 3–4 get `C_{-1}=ay³`, `C_{-2}=by⁴` with `ab≠0`
from `(a3)`'s unit `y⁷` and `(a5)`'s cap. Only the third — `c_{0,-10} ≠ 0`, GGV3's
corner-primitivity condition `(a6)` — has no derivation anywhere.

**That single integer is the whole obstruction.** Every condition our own §8
imposes is *closed* (equations, degree caps, order floors), and no boolean
combination of closed conditions can require a coefficient to be **nonzero**.

## 1. Why it cannot be fitted

The obvious move is to read `−10` off the corner data. It fails, and instructively.
At `(5,20)` **every** plausible formula collides:

```
-2*a0  -2*degC  -b0/2  -(a0+degC)  -degP1  -t*kappa-2  -(b0-2t-2)  -2t-2
 -10     -10     -10      -10       -10       -10         -10       -10
```

because there `a₀ = deg C = 5`, `b₀/2 = 2a₀ = 10`, `2t+2 = t·κ+2 = 10`,
`deg P₁ = 10`. And `(5,20)` is the **only** corner where GGV3 publishes
primitivity, so the calibration set is one maximally degenerate point.

`polygon_reduction.all_reductions()` returns three reductions at two corners —
`(8,28)` (ours, and the unique `t=4` corner that RETRACTS) and `(5,20)` twice.
There is no second point to fit against.

**So the number has to be derived, not fitted.**

## 2. GGV3 hands us the chart chain explicitly

This was the surprise. The map is not a research construction to be reverse
engineered — it is stated in the source (`1406.0886_GGV3.tex:1739-1742`):

> If `γ=3`, then applying to `(P₁,Q₁)` first the automorphism `x ↦ xy³`,
> `y ↦ y⁻²` of `K[x,y,y⁻¹]`, and then the automorphism `x ↦ x−G`, `y ↦ y` for
> some suitable `G ∈ K[y,y⁻¹]` …

and for `γ=2` (tex:1777-1780) the pair is `x ↦ xy²`, `y ↦ y⁻³`. Write the
substitution as

```
    x ↦ x·y^γ ,      y ↦ y^(−δ)          (γ,δ) = (3,2) or (2,3)
```

so a monomial transforms as **`x^i y^j ↦ x^i y^(γi − δj)`**.

*Consistency check.* That substitution has Jacobian `−δ·γ + 0 = −2` at `(γ,δ)=(3,2)`,
so `[P₁,Q₁] = x²` becomes `−2x²y⁶`, and after `x ↦ x−G` that is
`μ y⁶ (x−G)²` with `μ = −2` — which is exactly what `(a2)` prints. The y-exponent
`3γ−δ−1` gives `6` at `(3,2)` and `2` at `(2,3)`, matching `(a2)` and `(b2)`.

## 3. The derivation

`(a1)` gives `P = C²`, so `N(C) = ½·N(P₁)`. From the **computed** reduction
(`polygon_reduction.case_f2(0)`, not fitted):

```
N(P₁) = {(0,0), (0,10), (6,0), (8,2)}          deg P₁ = 10
N(C)  = {(0,0), (0, 5), (3,0), (4,1)}          deg C  = 5
cross-check: 3·N(C) == N(Q₁)                   ✓
```

The **x⁰ row** of `C` is the segment from `(0,0)` to `(0,5)` — i.e. `y^j` for
`j = 0…5`. Under the substitution, `x⁰ y^j ↦ x⁰ y^(−δj)`. With `δ = 2`:

```
predicted x⁰ support:  [0, −2, −4, −6, −8, −10]
GGV3 (a6) prints    :  c_{0,2}y² + c_{0,0} + c_{0,-2}y⁻² + ⋯ + c_{0,-10}y⁻¹⁰
MATCH on support    :  True
MATCH on depth      :  True
```

> **The primitivity depth is `−δ · deg_y(C|_{x=0})`**, and `c_{0,-10} ≠ 0` says
> exactly that the Newton-polygon vertex `(0,5)` of `C` is **attained**.

That is what corner primitivity *means*: the corner is a genuine vertex, its
leading coefficient nonzero. It is not an extra hypothesis bolted onto the
window system — it is the polygon's own corner, transported.

## 4. The control that makes this more than a coincidence

A fitted formula reproduces one integer. This reproduces **the whole support,
including its spacing**, and the spacing is a second, independent prediction: the
step in every window series must equal `δ`.

| series | γ | δ | observed step | predicted |
|---|---|---|---|---|
| `(a6)` `C₀` | 3 | 2 | **2** | 2 ✓ |
| `(b6)` `C₁` | 2 | 3 | **3** | 3 ✓ |
| `(b5)` `C₋₁` | 2 | 3 | **3** | 3 ✓ |

Three series, two charts, two different `δ`, all correct — including `(b5)`'s
eight-term run `1, −2, −5, −8, −11, −14, −17, −20`. The step-2 parity of `(a6)`
and the step-3 parity of `(b5)`/`(b6)` are *forced* by `y ↦ y^{−δ}` and by
nothing else in sight.

## 5. What is NOT established — read before citing

* **The depth law has one test.** `(b5)` and `(b6)` constrain `C₋₁` and `C₁`, not
  `C₀`; γ=2's `C₀` carries no published condition. They test the **step**, not
  the **depth**. So: step law — two charts; depth law — one instance.
* **The leading x-power is unreconciled.** `(a4)` has `C = x² + ⋯`, `(b4)` has
  `C = x³ + ⋯`, and `N(C)` above has x-degree 4. The depth argument uses only the
  x⁰ row, where these do not interfere, but a full **coefficient-level** bridge
  must resolve which x-graded piece of `C` becomes which `C_{-k}`. Until it does,
  this derivation gives the *slot*, not a map carrying one of our §8 witnesses to
  a value in that slot.
* **`(a1)`–`(a6)` remain GGV3's, asserted without proof** ("We do not provide
  proofs for this first part", tex:1716). This derives where `(a6)`'s number comes
  from; it does not prove `(a6)`.
* **γ is not pinned by the corner** (`gamma_from_corner.py`, 43 checks): `(5,20)`
  admits `γ ∈ {2,3,4}` and GGV3 asserts `{2,3}` without proof. Since `δ` depends
  on `γ`, **the primitivity depth is γ-dependent**, so any runnable `(a6)` test
  must be run per-γ — including the undischarged `γ = 4`.

## 6. What this changes

The compiler step is now **mechanical wherever a reduced polygon exists**:

```
get N(P₁)  →  halve it  →  read the y-axis vertex j*  →  depth = −δ · j*
```

The bottleneck is therefore **not** the chart bridge, and not the fitting: it is
that `N(P₁)` is known at only two corners. `SESSION_HANDOFF.md` §5 lists "a
reduced polygon at a class corner other than `(5,20)`" as item 4; the degeneracy
in §1 above says it is item **1**, and the other items are downstream of it.

A second corner would also discriminate sharply — the candidate families spread
to `−10 / −16 / −22` at `(8,32)`, `−10 / −18 / −26` at `(9,36)`. Note `−2a₀` and
`−b₀/2` are *identically equal on the whole class* (`b₀ = 4a₀`), so separating
those two needs an **off-class** corner; `(8,28)` is off-class (`b₀ = 28 ≠ 32`),
where they differ (`−16` vs `−14`).

## Files

| file | role |
|---|---|
| `PRIMITIVITY_DEPTH.md` | this writeup |
| `polygon_reduction.py` | `case_f2(0)` supplies the computed `N(P₁)` |
| `MOH_CONTROL_50_75.md` | §2 derives the other two required-nonzeros; §3 shows why our side cannot state this one |
| `ENDPOINT_CONTRACT.md` | §2 the kill predicate this feeds |
| `gamma_from_corner.py` | γ is not pinned by the corner, hence δ is not either |
