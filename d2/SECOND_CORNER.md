# The second corner

**Status: 2026-07-28. Delivered — one reduced Newton polygon that was already in
print and unnoticed, plus two newly derived ones, plus a proof that the residue
of the degeneracy cannot be removed by any corner.**

**Checker:** `second_corner_probe.py` — **67/67**, `--quiet`, exit 0, < 2 s,
read-only. Mutation-tested: perturbing the published `(7,21)` vertex list makes
it exit 1 on three checks.

---

## 0. The answer in one table

| corner `A_0` | `A_0'` | `l = t` | `κ` | `Δ' = N(C)` | `j*` | provenance |
|---|---|---|---|---|---|---|
| `(4,12)` | `(1,0)` | 3 | 1 | `{(0,0),(2,0),(3,1),(0,4)}` | 4 | derived |
| `(5,20)` | `(1,0)` | 4 | 2 | `{(0,0),(3,0),(4,1),(0,5)}` | 5 | derived; **GGV3's 3 integers** |
| **`(7,21)`** | `(1,0)` | 3 | 1 | `{(0,0),(2,0),(3,1),(0,7)}` | **7** | **PUBLISHED, GGHV22** |
| **`(7,42)`** | `(1,0)` | 6 | 4 | `{(0,0),(5,0),(6,1),(0,7)}` | **7** | **NEW** |
| **`(9,36)`** | `(1,0)` | 4 | 2 | `{(0,0),(3,0),(4,1),(0,9)}` | **9** | **NEW** |

`N(P) = m·Δ'`, `N(Q) = n·Δ'`, `N(C) = Δ'`, `[P₁,Q₁] = x^κ`, `deg P₁ = m·a₀`.

`j*` is `deg_y(C|_{x=0})` — the input `PRIMITIVITY_DEPTH.md` §6 needs.

Two of the three corners the task named are **blocked**, for two different and
precisely located reasons (§5). `(8,32)` is not merely hard: it is outside the
chart class entirely.

---

## 1. The correction that matters most: `(7,21)` was already in the repo

`PRIMITIVITY_DEPTH.md` §1 says

> `polygon_reduction.all_reductions()` returns three reductions at two corners …
> There is no second point to fit against.

That is true of `polygon_reduction.py` and **false of the repo**.
`passport_75_125.PUB` carries six rows transcribed from GGHV22 at four corners,
and one of them is a usable second point:

> **GGHV22, `paper_src/2204.14178.tex:1313-1320` (Proposition, Case (7,21)):**
> `[P,Q] = x`, `N(P) = {(0,0), 2(2,0), 2(3,1), 2(0,7)}`,
> `N(Q) = {(0,0), 3(2,0), 3(3,1), 3(0,7)}`.
> Figure caption, `:1384`: *"The transformation of ½N(P) = ⅓N(Q)"*.
> Proof, `:1388-1395`: the starting polygon `{(0,0),(1,0),(7,21),(0,7)}` times
> `(m,n)=(2,3)`; flip; `φ₂(y) = y + λx⁻²`; `φ₃(x)=x⁻¹, φ₃(y)=yx³`; result
> `{(0,0),(2,0),(3,1),(0,7)}`.

The paper writes the polygons **as `m` and `n` times one polygon** and names that
polygon in its own caption. So `N(C) = {(0,0),(2,0),(3,1),(0,7)}` and
`j* = 7 ≠ 5` — with no derivation of ours involved at all.

Why it was missed: of the six `PUB` rows, five are **en-split**
(GGV1 Prop 8.2(2)) — `(8,28)` ×2, `(9,24)` ×2, `(9,27)` — where `N(P)` is *not*
`m·Δ'` (`(8,28)`'s `N(P)` contains `(1,0)`, which is not divisible by `m=2`), so
no `N(C)` exists there and no `j*` can be read. `(7,21)` is the **only**
proportional row. Checks `A1`–`A6`, with the non-divisibility asserted row by row
as a mutation control.

**Consequence, immediately (§4): `(7,21)` alone kills four of the eight candidate
depth formulas.** The discrimination the 125 programme was waiting for did not
need a new derivation — it needed the repo to read its own control.

---

## 2. The missing input at a sporadic corner, and how to recover it

GGV5's **family** tables (`tex:1678-1694`, `:1709-1715`) print `A_0'`. Its
**sporadic** tables — the "9 other possible pairs with a complete chain of
length 1" (`tex:1828-1836`) and the length-2 table (`:1848-1858`) — print `A_0`,
`A_1`, `(m,n)`, `max deg` and **no
`A_0'` column**. `A_0'` is exactly the datum `passport_75_125.Reduction`'s rule
(r1) consumes, so its absence is what has kept every sporadic corner unreduced.

### The `A_0'`-recovery lemma

Given `A_0 = (u,v)` and the **printed** final corner `A_1 = (a\l, b)`:

1. `ρ := l`. GGV1 Prop `final` (7) writes `A^(1) = (a'/ρ, b')`, so the printed
   denominator **is** `ρ`.
2. Select the unique GGV1 branch `(f₁,f₂) = μ·(u,v)` whose
   `(ρ,σ) = dir(f₁−1, f₂−1)` has that `ρ`.
3. Solve `v_{ρ,σ}(A_0') = v_{ρ,σ}(A_0)` for `0 ≤ s' < r' < u`.
4. `γ := b` (the second coordinate of `A^(1)` is always `γ`).
5. **Test:** GGV1 (7), `A^(1) = A_0' + (γ−s')·(−σ/ρ, 1)`, must reproduce the
   printed **first** coordinate `a/l`.

Step 5 is a genuine prediction — one non-trivial rational per row.

### Calibration — 20 rows, then a refusal

| | |
|---|---|
| rows where GGV5 prints **both** `A_0'` and a fractional `A_1` | **20** (F₁–F₁₇, F₂₂–F₂₄) |
| recovered `A_0'` equals the printed one, and (7) reproduces `a/l` | **20 / 20** |
| rows the recovery **refuses** | **exactly** F₁₈–F₂₁ |

The refusal is the sharpest part. F₁₈–F₂₁ are exactly the four families **GGV5
itself proves cannot come from a standard `(m,n)`-pair** (`tex:1726-1786`, with
`(ρ₀,σ₀) = dir(A_0−A_0') = (1,0)` at `:1728`): their
`A_1` is an *integer* corner equal to `A_0'`, so `r' < u` fails by construction.
A recovery that reproduces 20 published values and refuses precisely the four
published impossibilities was not fitted to either. Checks `B1`, `B2`, `B2b`.

**Mutation controls.** Shifting the printed numerator `a` of `A_1` by ±1 makes
the recovery return nothing on all 20 rows (`B3`); every branch with the wrong
`ρ` predicts a different `A_1` (`B4`, 5+ wrong branches). Internal control: the
recovery reproduces `A_0' = (1,0)` at `(8,28)`, which the repo already knows by
an independent route (`B6`).

### Applied to the nine sporadic length-1 rows

| `A_0` | `A_1` | recovered `A_0'` | `γ` |
|---|---|---|---|
| `(7,35)` | `(19\7,5)` | `(2,0)` | 5 |
| **`(7,42)`** | `(13\7,6)` | **`(1,0)`** | 6 |
| `(8,28)` | `(7\4,3)` / `(11\4,7)` | `(1,0)` | 3 / 7 |
| **`(9,36)`** | `(17\9,4)` | **`(1,0)`** | 4 |
| `(11,33)` | `(19\4,8)` | `(3,1)` | 8 |
| `(12,33)` | `(11\3,8)` | `(1,0)` | 8 |

`(7,42)` and `(9,36)` land on `A_0' = (1,0)` — the one shape rule (r1) is
anchored on. `(12,33)` does too, but retracts (§5.4).

---

## 3. The class closed form

For `A_0' = (1,0)` **and** `a₀ | b₀`, everything is forced:

```
l = b0/a0,  mu = l-1,  c = b0 - mu*a0 = a0,  q = gcd(a0,b0) = a0,
zdeg = a0/q = 1      => rule (r5) admits NO two-factor split branch
b0 = l*(a0-1)        => impossible (it would need l = 0), so NO retraction,
                        hence C is a MONOMIAL
Prop 8.2(2) en-split => ILLEGAL for every k and both assignments
=>  Delta' = {(0,0), (l-1,0), (l,1), (0,a0)},   kappa = l-2,
    N(P) = m*Delta',  N(Q) = n*Delta',  N(C) = Delta',  j* = a0.
```

`deg C = max(l−1, l+1, a₀) = a₀`, because GGV1 Prop `u(u-1)`
(`1401.1784_GGV1.tex:3631-3632`, *"Then `v ≤ u(u−1)` and `u ≥ 4`"*) forces
`l ≤ a₀−1`. Checks `C1`–`C3`.

**Controls on the closed form.**

- `C0` — **rule (r1) is not a house convention.** GGHV22 *prints* the starting
  polygon `Δ = {(0,0),(1,0),A_0,(0,c)}` at four corners, and the rule's
  `c = b₀ − μa₀` comes out right at each:
  `(9,27)` `c=9` (`tex:471`), `(9,24)` `c=6` (`:682`), `(8,28)` `c=4` (`:1010`),
  `(7,21)` `c=7` (`:1388`).
- `C0b` — the one printed `Δ` carrying an **extra** vertex is `(9,27)`'s, and
  that vertex is `(9,24)` — which is that chain's `A_1` **and** its `A_0'`. That
  is the single published data point bearing on (r1) off `A_0' = (1,0)`, and it
  is what §5.2 turns on.
- `C4` — it reproduces GGHV22's published `(7,21)` polygon **exactly**. The form
  was not fitted to it: `(7,21)` was reached by the same three rules used at
  `(5,20)`.
- `C5` — **three** code paths agree on all five corners: the closed form,
  `passport_75_125.Reduction`, and a from-scratch flip/shift/invert written in
  the checker.
- `C6` — `κ = l−2` (the FUSED-CHART LEMMA of `composite_charts.py`: Jacobian
  `−x^(l−2)` for *any* shears) reproduces **both** published brackets: `[P,Q]=x`
  at `(7,21)` (`l=3`) and `[P₁,Q₁]=x²` at `(5,20)` (`l=4`, GGV3 `tex:1725`).
- `C7` — GGV3's other two integers at `(5,20)`: `deg P₁ = 10`, `deg Q₁ = 15`.
- `C8`/`C8b` — the en-split branch is illegal at all 40 (corner, k, swap) combos
  in the class, and the same test **passes** at `(8,28)` and `(12,33)`, so it
  discriminates.
- `C9` — the root-shift depth branch (`s = μ` vs `s = μ−1`, GGV6 Prop 2.5) moves
  **only** the foot vertex; `j* = a₀` is `s`-invariant, so the depth prediction
  does not depend on that unresolved branch.
- `C10` — `l = 2` or `4` at `(7,21)` fails to reproduce the published polygon, so
  the `l` rule is pinned there by print (a second pin, independent of GGV3's).

### The two new corners

```
(9,36)   Delta' = {(0,0),(3,0),(4,1),(0,9)}   kappa = 2
   (m,n) = (2,3):  N(P) = {(0,0),(6,0),(8,2),(0,18)}   deg 18
                   N(Q) = {(0,0),(9,0),(12,3),(0,27)}  deg 27
   (m,n) = (3,2):  the same two polygons, swapped
   GGV5 max{deg P,deg Q} = 3*45 = 135  -- matches the printed row

(7,42)   Delta' = {(0,0),(5,0),(6,1),(0,7)}   kappa = 4
   (m,n) = (2,3):  N(P) = {(0,0),(10,0),(12,2),(0,14)}  deg 14
                   N(Q) = {(0,0),(15,0),(18,3),(0,21)}  deg 21
   (m,n) = (3,2):  the same two polygons, swapped
   GGV5 max{deg P,deg Q} = 3*49 = 147  -- matches the printed row
```

These are **not** `(5,20)` replayed: the class now spans three `a₀` (5, 7, 9) and
three `l` (3, 4, 6), hence three `κ` (1, 2, 4). `(7,42)` in particular is the
**off-class** corner the task asked for — `b₀ = 6a₀`, so `−2a₀` and `−b₀/2` are no
longer identically equal there.

---

## 4. What this settles about the primitivity depth

The derived law is `depth = −δ · j*`, and `j* = a₀` on the whole class. So **`a₀`
must divide the depth**. Running `PRIMITIVITY_DEPTH.md` §1's eight candidates:

| candidate | `(5,20)` | `(4,12)` | `(7,21)` | `(7,42)` | `(9,36)` | verdict |
|---|---|---|---|---|---|---|
| `−2a₀` | −10 | −8 | −14 | −14 | −18 | **survives** (`δ=2`) |
| `−2·degC` | −10 | −8 | −14 | −14 | −18 | **survives** (identical) |
| `−(a₀+degC)` | −10 | −8 | −14 | −14 | −18 | **survives** (identical) |
| `−degP₁` | −10 | −12 | −14 | −14 | −18 | survives; killed by §4.2 |
| `−b₀/2` | −10 | **−6** | **−21/2** | −21 | −18 | **REFUTED** |
| `−t·κ−2` | −10 | **−5** | **−5** | **−26** | **−10** | **REFUTED** |
| `−(b₀−2t−2)` | −10 | **−4** | **−13** | −28 | **−26** | **REFUTED** |
| `−2t−2` | −10 | **−8** | **−8** | −14 | **−10** | **REFUTED** |

(Bold = not a multiple of `j*`; `−b₀/2` at `(7,21)` is not even an integer.)

Per-corner kill counts: `(5,20)`: 0 · `(4,12)`: 2 · **`(7,21)`: 4** · `(7,42)`: 1
· `(9,36)`: 3. Checks `E1`–`E3`.

### 4.1 The residual degeneracy is a **theorem**, not a gap

`−2a₀`, `−2·degC` and `−(a₀+degC)` are **identically equal wherever the law is
non-vacuous**, and no corner can ever separate them:

* `j* > 0` requires the top vertex `(l·a₀ − b₀, a₀)` to sit on the `y`-axis, i.e.
  `a₀ | b₀`;
* `a₀ | b₀` forces `c = a₀` (§3);
* GGV1 Prop `u(u-1)` forces `deg C = a₀` (§3).

So the separation would need `a₀ ∤ b₀`. Exactly **six** census corners have that:
`(6,15)`, `(8,28)`, `(9,21)`, `(9,24)`, `(12,30)`, `(12,33)`. **Four of them
retract**, and at a retracting corner `Δ'` has a *vertical* top face and no
`y`-axis vertex above the origin — `j* = 0`, the law is vacuous (verified). The
remaining two, `(9,21)` and `(12,30)`, both have `A_0' = (2,0)` and are blocked by
§5.2. Checks `E4`, `F3`–`F3c`.

### 4.2 `−degP₁` falls to an `(m,n)` argument, not a corner

At `(5,20)` the rows `(m,n) = (2,3)` and `(3,5)` share one corner, one `A_0'`, one
chart — hence one `γ` and one `δ` — yet give `−degP₁ = −10` and `−15`. A depth
that is a function of `(corner, γ)` cannot be both. (`E5`. This is
**CONDITIONAL** on the depth being `(m,n)`-independent; it is the one kill here
that is not pure arithmetic.)

### 4.3 Net

**Eight corner-data formulas collapse to one: `−2a₀ = −2j*`, whose entire content
is `δ = 2`.** The depth law's residue is a single `γ`-dependent integer — exactly
what `PRIMITIVITY_DEPTH.md` §5 flags, now with everything else stripped away.

**`δ` is NOT determined here, and this file does not claim it.** Declared
non-discriminating (`E7`): the class corners cannot test whether `δ` depends on
`γ`, because GGV3's two published charts (`γ=3, δ=2` and `γ=2, δ=3`) both sit at
`(5,20)` and only the `γ=3` one carries a printed depth. The recovered `γ` at the
new corners is 4 (`9,36`) and 6 (`7,42`) — neither is 2 or 3, so **neither new
corner has a published `δ` to borrow**. Fitting one from `a₀−γ` (which happens to
give `δ=2` at `(5,20)`, `γ=3`) would be a 2-anchor guess on a maximally
degenerate point — `a₀ = l+1 = ρ = 5` there — and is deliberately **not** done.

---

## 5. The blockers, as checked facts

### 5.1 `(8,32)` — outside the chart class, not merely hard

**No `A_0'` exists on any branch**: there is no `(r',s')` with `0 ≤ s' < r' < 8`
of equal `(ρ,σ)`-valuation on any surviving `(f₁,f₂)`. The corner cannot be given
a first chart at all (`F1`).

The reason is structural. `(8,32)`'s chain is length 2 with `A_1 = (8,28)` an
**integer** corner, and `v_{1,0}(8,32) = v_{1,0}(8,28) = 8`, i.e.
`(ρ₀,σ₀) = (1,0)` — the F₁₈–F₂₁ shape (GGV5 `tex:1731`), which is the `ℓ_{10}`
layer, **not** the type-II.b root-shift + Laurent chart this machinery
implements (`F1b`).

Exhaustively (`F1c`): across all 21 census corners the "no `A_0'`" verdict lands
on exactly `(6,18)`, `(6,24)`, `(8,32)`, `(8,40)` — and every one of them is a
corner whose chain's first link ends at an integer `A_1`. It is that shape being
detected, not a search failure.

This also **corrects** the task's framing: `(8,32)`'s unreachability is not that
"no branch survives at `u=8`" in the shared-`f=(4,16)` sense; it is that *every*
branch fails, because the corner's first link is a different kind of step.

### 5.2 `(10,40)` — one undetermined polygon vertex

> **RESOLVED 2026-07-28, same day — and the question was a false dichotomy.**
> `r1_vertex_reading.py` (12/12) scores the two readings below against all four
> printed polygons: **A scores 3/4 and B scores 3/4, both failing at the same
> corner, `(9,27)`.** Its printed set has FIVE vertices,
> `{(0,0),(1,0),(9,24),(9,27),(0,9)}`, containing `(1,0)` **and** `A_0' = (9,24)`
> at once — A cannot produce `(9,24)`, B cannot produce `(1,0)`. The **union**
> `Δ = hull{(0,0), (1,0), A_0', A_0, (0,c)}` scores **4/4** and is the only one of
> the three consistent with print.
>
> So the fifth published `Δ` this section asks for is **not needed**: the fourth
> already discriminates, against both stated readings. Where `A_0'` lies on the
> x-axis — exactly the blocked `(2,0)` case — `(1,0)` is absorbed into the segment
> to `(2,0)`, so U and B coincide and only A differs; that is why the ambiguity
> looked undecidable from the three on-axis corners alone.
>
> At `(10,40)`, U predicts `hull{(0,0),(2,0),(10,40),(0,10)}`, consuming no
> unpublished datum since `c = 10` is pinned independently (below).
> **Scope:** one discriminating instance; U is the unique survivor of the three
> readings tested, not proved to be the unique rule fitting four points. The other
> five blocked corners still need `c` there. The text below is retained as the
> statement of the blocker as it stood.


`A_0' = (2,0)` **is** recovered, and anchored **twice**: both printed final
corners `(16\5,6)` and `(18\5,8)` are reproduced by the same `A_0'` with `γ = 6`
and `8` (`F2`). So the datum the sporadic table omits is available.

What blocks it is one vertex of rule (r1). (r1) writes
`Δ = hull{(0,0), (1,0), A_0, (0,c)}`. Two readings of that `(1,0)`:

* **A** — the vertex is always `(1,0)`, a normalisation;
* **B** — the vertex is `A_0'`.

(r1) is printed **only** at the four corners of `C0`: `(7,21)`, `(8,28)`,
`(9,24)`, `(9,27)`. At the first three `A_0' = (1,0)` and the two readings
**coincide** — nothing in print distinguishes them (`F2b`). The fourth,
`(9,27)`, is the single data point off `(1,0)`: its printed `Δ` carries the extra
vertex `(9,24)`, which is exactly that chain's `A_0' = A_1` — evidence for **B**,
but one instance, at a corner whose `A_0'` GGV5 does not print (`F2b2`).

At `(10,40)` the readings **disagree**:
`hull{(0,0),(1,0),(10,40),(0,10)}` vs `hull{(0,0),(2,0),(10,40),(0,10)}` (`F2c`).
Same blocker at every `A_0' ≠ (1,0)` corner: `(7,35)`, `(8,24)`, `(9,21)`,
`(10,40)`, `(11,33)`, `(12,30)`.

> **This is the precise missing input.** It is one polygon vertex at one corner
> class, and it is a missing **published** datum — a *fifth* printed `Δ`, at any
> corner with `A_0' = (2,0)`, in the form of `C0`'s four, would settle it — not missing
> mathematics. It is *not* the "no general chain-to-chart dictionary" problem of
> `CURRENT_STATUS.md` correction 3; that obstruction is about chains, and this
> one is a single vertex.

Partial consolation (`F2d`): both readings put the top vertex at
`(l·a₀ − b₀, a₀) = (0,10)`, so **`j* = 10` at `(10,40)` either way**. What is
undetermined is the rest of `Δ'`, not the depth input. If a future step needs only
`j*`, `(10,40)` is usable now.

### 5.3 The length-2 corners generally

`(8,32)`, `(8,40)`, `(9,27)`, `(9,36)`'s *second* row, `(12,36)` — every length-2
chain whose first link ends at an integer `A_1` is out of scope by §5.1. The
length-2 chains with a **fractional** `A_1` — `(10,40)` ×2, `(12,30)`, `(12,36)`
×2 — all recover an `A_0'` cleanly (`(2,0)`, `(2,0)`, `(3,0)`) and are blocked
only by §5.2.

### 5.4 `(12,33)` — a third computed polygon that does not help

`A_0' = (1,0)` recovered, chain length 1, so the engine applies and gives
`Δ' = {(0,0),(2,0),(3,11),(3,12)}`, `κ = 1`. But `(12,33)` **retracts**
(`33 = 3·11`), so `Δ'` has a vertical top face, `j* = 0`, and it does not feed the
depth law; and its en-split branch is **legal**, so the proportional branch is not
forced and `N(C)` is not even well defined there. Recorded, not used (`F4`).

---

## 6. Evidence boundary

| level | what |
|---|---|
| **PROVED** | `κ = l−2` (FUSED-CHART LEMMA, `composite_charts.py`). The class forcing of §3 given (r1)+(r4)–(r7): `a₀\|b₀ ⟹ l=b₀/a₀, c=a₀, zdeg=1, no retraction`. `deg C = a₀` from GGV1 Prop `u(u-1)`. §4.1's identity of the three survivors, and the theorem that no corner can separate them. |
| **EXACT-CHECKED** | The 20/20 `A_0'` calibration and the 4/4 refusal. The four printed `Δ`'s of `C0`. The three-code-path agreement on five corners. The en-split illegality (40 branches). The `s`-invariance of `j*`. The eight-candidate divisibility table. Every mutation control. 67/67 in `second_corner_probe.py`. |
| **CITATION-LEVEL** | GGHV22's `(7,21)` polygon and `[P,Q]=x` (`2204.14178.tex:1313-1320`), and its four printed `Δ`'s (`:471, 682, 1010, 1388`); GGV3's `[P₁,Q₁]=x², deg 10, deg 15` (`1406.0886_GGV3.tex:1723-1727`); GGV1 Prop `final` (5)–(9) and Prop `u(u-1)` (`1401.1784_GGV1.tex:3631-3632`); GGV5's tables and its F₁₈–F₂₁ impossibility. None re-derived here. |
| **INFERRED** | Rule (r1) is **printed at four corners** (`C0`) — its *generalisation* to a fifth is what is inferred, and the undetermined part is exactly the one vertex of §5.2; `CORNER_RESOLVENT.md` §5.1's "no general dictionary exists" still stands for chains. `l = ceil(b₀/a₀)` inherits `CORNER_ATLAS.md` §4's INFERRED status, now with two published pins (`(5,20)` via GGV3, `(7,21)` via GGHV22's own `φ₃(y) = yx³`). The recovery lemma's step 1 — "the printed denominator of `A_1` **is** `ρ`" — is INFERRED, with **20 anchors**. |
| **NOT ESTABLISHED** | `δ` at any corner other than `(5,20)`. `δ = a₀ − γ` is a **2-anchor guess on a degenerate point** and is explicitly not adopted. `(a1)`–`(a6)` remain GGV3's, asserted without proof (`tex:1716`). Whether `P₁ = C^m` (as opposed to `P₁ = C²` only) holds off `m = 2`. |

---

## 7. What to do next

1. **`δ`, and nothing else.** §4.3 reduces the whole depth question to one
   `γ`-dependent integer. The cheapest attack is GGV3 `tex:1739-1742` /
   `:1777-1780` read as a *rule* rather than two instances: what determines the
   pair `(γ,δ)` from the chart? At `(5,20)` four candidate rules
   (`a₀−γ`, `l+1−γ`, `ρ−γ`, …) coincide because `a₀ = l+1 = ρ = 5`. At `(9,36)`
   they spread to `δ ∈ {1,5}` — so **`(9,36)` discriminates `δ`-rules even though
   it does not discriminate depth formulas.** That is the highest-leverage use of
   the new corner.
2. **One published `A_0' = (2,0)` reduction** unblocks six corners at once
   (§5.2), including both `(10,40)` rows and both corners that would otherwise
   separate `−2a₀` from `−2c`.
3. `PRIMITIVITY_DEPTH.md` §1 and §6, and `SESSION_HANDOFF.md` §5 item 4, should be
   updated: the bottleneck as stated ("`N(P₁)` is known at only two corners") is
   retired. It is known at five, one of them in print.

## Files

| file | role |
|---|---|
| `SECOND_CORNER.md` | this writeup |
| `second_corner_probe.py` | the checker, 67/67 |
| `passport_75_125.py` | `Reduction` (the engine) and `PUB` (the six published rows, incl. `(7,21)`) |
| `gamma_from_corner.py` | GGV1 Prop `final` (5)–(9), on which the recovery lemma is built |
| `composite_charts.py` | FUSED-CHART LEMMA, `κ = l−2` |
| `PRIMITIVITY_DEPTH.md` | the depth law this feeds; §1 and §6 superseded |
| `CORNER_ATLAS.md` | the 34-row census and the `t`-rule evidence base |
