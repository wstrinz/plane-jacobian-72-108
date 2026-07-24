# The F2 certificate tower: the (50,75) kill does NOT extend to (75,125) — BLOCK-OBSTRUCTION

## Verdict

**BLOCK-OBSTRUCTION.** The a=2 case (50,75) kill (GGV3 §5) is reproduced
**exactly**, and its certificate kind is identified: a small set of *terminal
coefficient equations* whose elimination forces a corner window-coefficient to
vanish — a **bigraded (u-weight, y-order) window-depth contradiction**, not a
scalar syzygy on the G-system generators. The **algebraic** five-block layer of
the D-transform G-system extends cleanly to a=3 (+5 forcing generators, +5 spare
window unknowns, the Φ recurrence, and a literal generator nesting). But the
**kill layer** — the y-order / window-cap layer where the (50,75) contradiction
actually lives — does **not** extend by the fixed five-block rule: the
window-denominator invariant `q_window = 5a-3` jumps `7 → 12` with `gcd(7,12)=1`
(incommensurate window lattices), and the y-order fractional-denominator set of
the forcing slices fragments from `{1,7}` at a=2 to `{1,2,3,4,6,12}` at a=3.

**Consequence for (75,125): it is NOT killed by tower extension.** The specific
obstructing block is the window-cap / y-order layer — the 5 new forcing
generators plus the deepened Φ slice, which carry the new fractional-denominator
classes `{2,3,4,6,12}` absent at a=2. This **is** the "new geometry at 125" the
review predicted: an incommensurate period-12 window lattice. (75,125) would
require a fresh period-12 window compiler, not a fixed-block lift of the (50,75)
certificate.

Everything below is exact (sympy). Independent checker `f2_tower_verify.py`
(`--quiet`, **exit 0**, all checks pass).

---

## CORRECTION (2026-07-24): "incommensurate" → COPRIME / NONALIGNED (obstruction softened)

> The word **"incommensurate lattices"** used below (Verdict, §2b, §3) is
> **imprecise and overstates the obstruction.** Periods **7** and **12** are
> `gcd(7,12)=1` **coprime**, hence the two window lattices `(1/7)Z` and
> `(1/12)Z` are **finite-index sublattices of a common refinement** — they share
> the refinement `(1/84)Z` (`lcm(7,12)=84`). They are therefore **commensurate**
> (any two rank-1 rational lattices are); the correct descriptor is
> **coprime / nonaligned periods**, not "incommensurate."

**What the obstruction actually is, precisely.** The BLOCK-OBSTRUCTION verdict is
an obstruction to a **period-PRESERVING block map** — a fixed five-block lift that
carries the a=2 certificate's single period-7 window pattern onto the a=3 system
**while keeping the period fixed at 7**. That fails because a=3's forcing slices
populate the divisor classes `{2,3,4,6,12}` of period 12, which have **no image
mod 7**. This is a genuine obstruction to the *fixed-block, period-locked*
extension, and only to that.

**What is NOT obstructed (the softening).** The obstruction says **nothing** about
towers built after **refining to the common lattice**. Passing both rungs to the
refinement `(1/84)Z` — or, equivalently, to the full two-coordinate window
lattice `Z^2 (u-weight, y-order)` before projecting to any single period — makes
the a=2 pattern and the a=3 pattern live on **one** lattice, and the block map is
no longer required to preserve a period. **A refined-lattice tower remains OPEN**:
nothing here rules it out. This refined-lattice / period-84 (equivalently
bigraded `Z^2`) window compiler is precisely **the engine's target** — the object
the "fresh period-12 window compiler" of §3 was gesturing at, now correctly stated
as a common-refinement construction rather than an impossibility.

*(Committed STATE.md history that records the old "incommensurate" language is
left unchanged — it is a historical entry; this block is the standing
correction.)*

---

## 1. The a=2 certificate (GGV3 §5), reproduced exactly

GGV3 (arXiv:1406.0886) §5 kills (50,75) in **two reduced charts** `γ=3, γ=2`
(the automorphisms `x↦xy³,y↦y^-2` and `x↦xy²,y↦y^-3` of the (5,20) corner). In
both, an auxiliary Laurent series `Z` supplies "E-equations" — `E_k=(Z²)_{-k}`
(the **linear window**, `P=C²`) and `E=(Z³+λZ^-1)_{-k}` (the **forcing window**,
`Q=C³+λC^-1+F`). This is exactly the linear-window / forcing-window split of our
D-transform G-system: `E_k=0` ↔ the linear eliminations; the forcing `E=-F` ↔
the generators.

**γ=3 (window-depth kill).** Eliminating the deep window unknowns `C_{-3..-7}`
from `E_1..E_5=0, E_6=-F_{-1}, E_7=-F_{-2}` gives (checker §B, exact)
`F_{-1}=-3C_{-1}C_{-2}` and the forcing relation
`3C_0C_{-1}² − 3C_{-2}² − 2λ = 2F_{-2}`. Imposing the corner window caps (a5)
`C_{-1}=ay³, C_{-2}=by⁴` yields

```
C_0 = b²y²/a² + 2f_8 y²/(3a²) + 2f_6/(3a²) + 2f_4/(3a²y²) + 2f_2/(3a²y⁴) + 2λ/(3a²y⁶),
```

whose **lowest y-power is `y^-6`**. Corner primitivity (a6) demands
`c_{0,-10} ≠ 0`, i.e. a nonzero `y^-10` coefficient. The forced `C_0` has none —
**contradiction**. The kill is a mismatch between the *forced* window depth
(`y^-6`) and the *required* window depth (`y^-10`).

**γ=2 (terminal elimination + square/window-depth kill).** The 13-equation
terminal system (paper lines 2027–2042), on eliminating `{a,e_{-10},e_{-7},
e_{-4},e_{-1},λ}`, gives the elimination ideal (checker §A, exact **Gröbner**)

```
⟨ g_{-2}^5, g_{-2}^4 g_{-5}, g_{-2}^2 g_{-5}^2, g_{-2} g_{-5}^3, g_{-5}^4 ⟩,
```

so `g_{-2}^5 = g_{-5}^4 = 0` (matching the paper verbatim), hence
`g_{-2}=g_{-5}=0` and `F_{-4}=0`. Then the derived relation
`3a²C_1²y² = (4a³−8)y³ − 8λ`, whose right side must be a **perfect square in
`K((1/y))`**, has odd top y-degree 3 unless `4a³−8=0`, forcing `a³=2`; and then
`C_1` is homogeneous of y-degree `-1`, i.e. `e_{-10}=0` — contradicting
primitivity (b6) `e_{-10}≠0` (as `C_0=0` was forced). Again a **window-depth**
contradiction: the forced support (`y^-1` only) falls short of the required depth
(`y^-10`).

**Certificate kind (from the menu, in order tried).** Not a left-null covector,
Koszul syzygy, or determinant minor on the *scalar* G-system (that system has the
trivial origin solution — the scalar ideal does not die). It is **a small set of
terminal coefficient equations** whose elimination ideal forces a corner
window-coefficient to vanish — realized in the y-graded (window-cap) layer. The
kill is intrinsically **bigraded** `(u-weight, y-order)`; the y-order axis carries
it.

---

## 2. The tower step: what transfers and what does not

Build the a=2 and a=3 systems in our D-transform `t=5` chart with the landed
parametric builder `g_system_75_125.build_gsystem` (read-only import).

### 2a. The algebraic block layer — transfers EXACTLY

| quantity | a=2 = (50,75) | a=3 = (75,125) | rule |
|---|---:|---:|---|
| forcing generators | 5 (`j=1,2,3,4,6`) | 10 (`j=1..9,11`) | **+5** |
| spare window unknowns | 4 (`dm2..dm5`) | 9 (`dm2..dm10`) | **+5** |
| Φ u-slice `M` | 21 | 36 | `bt+jphi` |
| u-weights of generators | 16..21 (skip 20) | 26..36 (skip 35) | AP |

- **Φ recurrence** `Φ_{a+1} = (a/(a+1)) C^{30a+3} Φ_a` holds exactly (checker §E,
  a=2,3,4), `C=y²(y³+1)`.
- **Generator nesting** (checker §F): the a=2 generator `G1` is the leading
  `d0²` block of the a=3 generator `G1`:
  `coeff_{d0²}(G1^{a=3}) = (10/3)·G1^{a=2}`. The linear-window jump `S²→S³`
  multiplies the old block by the extra `S`-copy and appends the new deep terms.

So the **structure** — five-block growth, the AP of u-weights, the Φ recurrence,
the nested generators — is a genuine fixed-block tower.

### 2b. The kill layer — does NOT transfer

The (50,75) kill of §1 lives **entirely** in the y-order / window-cap layer. That
layer is governed by the window-denominator invariant

```
q_window := denom( W_step ),   W_step = ord_y(Φ)/M = (30a²−24a+3)/(5(2a−1)+ (5a−4)),
```

which for the F2 family is **`5a−3`** (a=2: `W_step=25/7`, `q_window=7`; a=3:
`W_step=67/12`, `q_window=12`), always with `gcd(2a−1, 5a−3)=1`. The two rungs'
window lattices are **incommensurate**: `gcd(7,12)=1`.

Concretely, the y-orders of the forcing slices (physical order `W_step × u-weight`)
have fractional denominators (checker §G):

| a | forcing-slice y-order fractional denominators |
|---|---|
| 2 | `{1, 7}` — a single quasi-period-7 class |
| 3 | `{1, 2, 3, 4, 6, 12}` — the full divisor lattice of 12 |

At a=2 every non-Φ generator sits at one uniform fractional class (denom 7); at
a=3 the generators **spread across six** fractional denominators. The
certificate's y-order bookkeeping — a single period-7 window pattern at a=2 — has
no image under a fixed five-block extension into the period-12 divisor lattice.
The 5 new generators (and Φ's deepening from slice 21 to 36) introduce the new
fractional classes `{2,3,4,6,12}` that simply do not exist mod 7.

**This is the named obstructing block:** the window-cap / y-order layer (the 5 new
forcing generators + the deepened Φ slice). The algebra nests; the arithmetic of
the window lattice does not.

---

## 3. What it means for (75,125)

- The (75,125) case is **not killed** by extending the (50,75) certificate. The
  tower's algebraic skeleton lifts, but the actual contradiction is carried by an
  incommensurate window lattice that the fixed-block rule cannot produce.
- The obstruction is *located and named*: it is the period jump `7 → 12`
  (`q_window = 5a-3`, `gcd=1`), materialised as the new fractional y-order classes
  `{2,3,4,6,12}` on the five new forcing generators — the concrete form of the
  review-predicted "new geometry at 125."
- A kill of (75,125) along these lines would require a **fresh period-12 window
  compiler** (the bigraded / two-coordinate window lattice flagged in
  `f2_family_verify.py`'s window-denominator correction and
  `G_SYSTEM_75_125.md`'s a=3-boundary result), not a lift of the a=2 object. This
  is consistent with (75,125) remaining open.

### Honest scope / judgment notes

1. **[a=2 presentation nuance]** The reproduced a=2 kill is GGV3's own, in its
   two **reduced** charts (γ=2, γ=3). It is *not* verbatim a certificate on our
   D-transform `t=5` G-system: that scalar system does not die (trivial origin
   solution). The GGV3 charts supply the y-order/window layer the scalar
   G-system lacks. The linear/forcing-window *split* matches ours exactly; the
   decisive step is y-graded. This is the "A2-reconstruction nuance."
2. **[inherited]** The (5,20)→(7/5,2) reduction, `C=y²(y³+1)`,
   `Φ_a=−(1/(3a))y^{...}(y³+1)^{...}`, and the F2 family theorem are the landed
   phase-1 objects (`C_SERIES_75_125.md`, `f2_family_verify.py`), themselves
   conditional on the standard unreduced-polygon chart (judgment 2 there).
   Unchanged here.
3. **[derived, this lane]** The five-block growth, Φ recurrence, generator
   nesting, `q_window=5a-3` law, and the `{1,7} → {1,2,3,4,6,12}` fractional
   fragmentation are computed exactly from the built systems and checked in
   `f2_tower_verify.py`. The verdict BLOCK-OBSTRUCTION answers the well-posed
   tower question negatively for the *fixed-block extension of the a=2
   certificate*; it does **not** claim (75,125) is unkillable by other means.

---

## Files

- `F2_TOWER.md` — this writeup.
- `f2_tower.py` — the construction: reproduces the GGV3 §5 a=2 kill (both
  charts), builds the a=2/a=3 G-systems, exhibits the five-block growth, Φ
  recurrence, generator nesting, and the window-lattice obstruction. Exact
  sympy; run end to end.
- `f2_tower_verify.py` — independent exact PASS/FAIL checker (`--quiet`, exit 0);
  §A GGV3 γ=2 terminal elimination `→ g_{-2}^5=g_{-5}^4=0`; §B γ=3 forced-`C_0`
  window-depth contradiction; §C the square / e_{-10}=0 window-depth obstruction;
  §D five-block growth; §E Φ recurrence; §F generator nesting; §G the
  `q_window=5a-3` jump and fractional-denominator fragmentation.
