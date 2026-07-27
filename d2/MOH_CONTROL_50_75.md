# The (50,75) tension resolves: our §8 witnesses survive there, and what kills them is a condition our machinery cannot state

**The tension was not a contradiction. The transfer is SOUND and INCOMPLETE, in a
nameable way — and the name is the one `YPLACE_TRANSFER.md` §8 already wrote down
as its untouched lead #2.**

Two committed results looked incompatible:

* `yplace_transfer.py` (57/57, `9afcb79`) — at a class row's `y`-place the whole
  of `PROOF_72_108` §§3–7 transfers; the row lands at `k = 0`, Cor 8.5; and **four
  explicit witnesses satisfy the entire §8.1 `k = 0` system with residual 0**
  inside every transferred cap.
* `moh_discards.py` (21/21, `a54a63e`) — `F_2(2,3)/75` **is** `(50,75)`, which
  GGV3 §5 kills outright in two γ-charts; and `(5,20)` is a class corner
  (`b₀ = 4a₀`).

**Answer.** Every condition our transfer imposes is **closed** — five polynomial
equations, six degree caps, five order floors. GGV3's kill is **open**: corner
primitivity `(a6)` requires one *named* window coefficient to be **nonzero** at a
depth the window equations force to vanish. A degree cap and an order floor are
satisfied by *more* vanishing, never less, so no combination of the conditions we
transfer can require a coefficient to be nonzero. Our witnesses are points of a
strictly weaker condition set than "is a germ", and a point of it need not lift.

Three results, in order of value:

| | statement | status |
|---|---|---|
| **1** | **The witnesses survive at `(50,75)` — and `(50,75)` is where they are MOST constrained.** The caps they are tested against are corner `(5,20)`'s own, the tightest of the four class corners; the other three are `15/4, 4, 17/4`. Nothing in the transferred chain reads `(deg P, deg Q)`, so the same four points serve `F_2(2,3)/75` and its `P↔Q` swap `F_3(3,2)/75` **identically** — no asymmetry exists to find. | **EXACT-CHECKED** |
| **2** | **The killer is corner primitivity `(a6)`, and the whole γ=3 kill is DERIVED here from `(a1)`–`(a6)`, not replayed.** `E_1..E_8` (γ=3) and `E_1..E_13` (γ=2) are computed from the `Z`-series definition and reproduce GGV3 **term for term**; the two displayed conclusions come out of the elimination. Zero margin: `deg_y C_{-1} + deg_y C_{-2} ≤ 3 + 4 = 7 = deg_y F_{-1}`. | **PROVED** from `(a1)`–`(a6)`; those premises **CITATION-LEVEL** |
| **3** | **It does NOT transfer to `(75,125)` — but the γ layer does, verbatim.** Same `A_0 = (5,20)`, same unique branch, same `A_0' = (1,0)`, same `γ ∈ {2,3,4}`, same standing γ=4 obligation. The kill layer moves: the computed reduced pair is `(15,25)`, not `(10,15)`, and the E-depth law `K − L = c(n−m)` goes from `c` to `2c`. | **EXACT-CHECKED** / **PROVED** (the depth law) |

**Checker:** `moh_control_50_75.py` — **48/48**, `--quiet` exit 0, ~12 s.
**NOT GATED**, stated so it is not inferred from absence: this lane did not touch
`tools/suite_manifest.py` (it lives outside `d2_plane_72_108/`, and
`EXPECTED_TOTAL` is a shared count). **Nor is `yplace_transfer.py` gated** —
grepping the manifest for `yplace` returns only `slice_phi_yplace.py`, a
different file. So **half one of the tension currently rests on an ungated
checker**; it was re-run green (57/57) while writing this. Both should be added
by whoever next edits the manifest.
**Sources of truth, nothing retyped where a module has it:** corner data from
`polygon_reduction.corner_chart_data` (guarded) and `.case_f2(0)/.case_f2(1)`
(computed reduced polygons); γ-admissibility from `gamma_from_corner.analyse`;
row data from `corner_atlas.json`; `(a1)`–`(a6)`/`(b1)`–`(b6)` and the published
`E`-lists quoted from `paper_src/1406.0886_GGV3.tex` at their point of use;
§8.1's generators retyped from `PROOF_72_108.md` (8.1.1) as printed.

---

## 0. Identity — the numbers describe the case the label names

This is the check the repo has been burned on, so it is first and it is
redundant on purpose.

* GGV3 §5 says, in its own words, that it verifies `deg(P₀)=50`, `deg(Q₀)=75`;
  that `A_0 = (5,20)`; and that the reduced pair satisfies `[P₁,Q₁] = x²`,
  `deg P₁ = 10`, `deg Q₁ = 15` [A1–A3].
* `corner_atlas.json`'s `F_2(2,3)/75` carries `A_0 = (5,20)`, `max_deg = 75`,
  `(m,n) = (2,3)` — and `50:75 = 2:3` [A4].
* `(5,20)` satisfies `b₀ = 4a₀`, the **class** shape, not `(8,28)`'s retraction
  shape `b₀ = 4a₀ − 4` [A5]; `polygon_reduction.corner_chart_data(5,20)` returns
  `t = 4, κ = 2, C = y` monomial, `retraction = False` — exactly
  `yplace_transfer`'s class-row chart `(a,b,t,κ) = (2,3,4,2)` [A6].
* `polygon_reduction.case_f2(0)` is itself **labelled** `(50,75)`
  (`tag = 'F2_j0_50_75'`, `signature['degs'] = (50,75)`), with `A_0' = (1,0)`,
  reduced pair `(2,3)`, `l = 4`, bracket `x²` — so the in-repo reduction and
  GGV3 §5 agree on the **chart**, not merely on the corner [A7].
* **Cross-check between two objects that must agree.** The *computed* reduced
  polygons give `max_y N(P₁) = 10`, `max_y N(Q₁) = 15` — which **are** GGV3's
  published `deg P₁ = 10`, `deg Q₁ = 15` [A8]. Two independent closed forms agree
  too: `deg P₁ = m(t+1) = 10 = deg P / a₀`, `deg Q₁ = n(t+1) = 15 = deg Q / a₀`
  [A9].
* `F_3(3,2)/75` is the `P↔Q` swap: same corner, same `t, κ`, and the atlas's own
  sorted reduced pair `(D_P,D_Q) = (2,3)` is **identical** [A10].
* The class of nine splits **8 + 1** exactly as `YPLACE_TRANSFER.md` scopes it:
  eight rows with reduced pair `(2,3)` — of which the two `/75` rows are two —
  and one with `(3,5)`, namely `F_2(3,5)/125 = (75,125)` [A12].

> **So `(50,75)` is INSIDE `yplace_transfer`'s stated scope, and `(75,125)` is
> outside it by that document's own words.** The tension is entirely inside the
> eight, and two of the eight are settled in the literature.

---

## 1. The witnesses survive at `(50,75)` — at the tightest caps in the class

The four §8 witnesses were rebuilt here from `PROOF_72_108` (8.1.1) as printed
(`Π = 1`, `B = Πv = v`, `u := γd₂`, `w := ½γ²d₁Π`), not imported. For
`(A,z,ζ,γ) = (y,2,1,1), (y,3,1,1), (y²,4,1/3,1), (0,5,2,1)` the full `k=0`
system `g₁ = g₂ = g₃ = □ = 0` **and** `(*) FZ = (1/6)γ⁵y⁹` holds with residual
**exactly 0**, inside all six degree caps and all five order floors [B4].

And the caps are not generic. The slope reader used here is a deliberate
re-implementation, calibrated first on the **published** control: fed `(8,28)`'s
computed reduced hulls it returns `(σ,τ) = (2,1)` and `(2,2)`, giving
`ord_y h_k ≥ 12k`, `deg_y h_k ≤ 15k / 14k` — `PROOF` §2.6(iii)'s published
`ord D_{j_x} ≥ 48 − 12j_x` and `λ = 3 / 2` [B0]. Fed the **computed** reduced
polygon `N(P₁) = {(0,0),(0,10),(6,0),(8,2)}` at `(5,20)` it returns
`(σ,τ) = (1,−1)`, so §2.6's affine induction gives `ord_y h_k ≥ k`,
`deg_y h_k ≤ 3k`, cap-`λ = 2` [B1] — which is what the shipped atlas records for
this row (`G3.lam = 2`, PASS) [B2]. The six caps `(A,u,v,w,d₀,C) ≤
(9,6,12,9,12,15)` are that slope-3 ledger `3w = (6,9,12,15,18,21,24)` stripped by
`y⁹` [B3]. `YPLACE_TRANSFER` §D5 records the other three class corners at deg
slopes `15/4, 4, 17/4` — all looser [B5].

> **`(50,75)` is the corner where the transferred machinery is most constrained,
> and the witnesses survive there.** So the answer to "does `(50,75)` differ from
> the generic class row in a way that already kills them?" is **no** — and it is
> a strong no, because the "generic class row" caps `yplace_transfer` used were
> `(5,20)`'s all along.

Nothing in the transferred chain reads `(deg P, deg Q)` — the witnesses depend
only on `(a,b,t,κ,C) = (2,3,4,2,y)`, `ord_y Φ = 30` and the `(5,20)` caps. So
they are the **same four points** at `F_2(2,3)/75` and at `F_3(3,2)/75` [B6],
and the `P↔Q` asymmetry test returns nothing: there is no asymmetric object to
differ.

### The inventory that settles it

`PROOF` §8 at `k = 0` states **exactly five** nonvanishing conditions — `γ ≠ 0`
(§8.4), `μ ≠ 0` (§8.5), `ζ ≠ 0` (8.4.1), `F ≠ 0` and `Z ≠ 0` (§8.3 `(*)`) — and
the sixth candidate, `A(r) ≠ 0` at roots of `Π` (§8.2), is **vacuous** because
`Π = 1`. All four witnesses satisfy every one [B7].

**There is no sixth. §8 never requires a named coefficient at a named depth to be
nonzero.** That is the whole story.

---

## 2. The killer, derived from `(a1)`–`(a6)`

`f2_tower.a2_certificate()` writes GGV3's γ=2 system down as **13 literals** and
supplies `a³ = 2` as a **given**; it references `T`, `KAPPA`, `QC`, `C`,
`ordPhi`, `Nof`, `build_gsystem` zero times each (SESSION_HANDOFF's "REPLAY
TRAP"). This section does not do that.

**Step 1 — the `E`-system, from the definition.** From `(a4)`
`C = x² + C_0 + C_{-1}x^{-1} + ⋯` and `(a1)` `P = C²`, `Q = C³ + λC^{-1} + F`,
form `E_k := (Z²)_{-k}` for `k = 1..5` and `E_{5+k} := (Z³ + λZ^{-1})_{-k}` for
`k = 1..3`. This reproduces GGV3's published `E_1..E_8` **term for term** [C1],
including the paper's own remark that `Z_{-7}` is the deepest coefficient
occurring and occurs only in `E_5` and `E_8` [C2]. That validates the series
conventions before anything is eliminated.

**Step 2 — the elimination.** Imposing `E_1 = ⋯ = E_5 = 0` collapses `E_6` to
exactly `3C_{-1}C_{-2}`, so `(a3)`'s `E_6 = −F_{-1}` **is** the paper's displayed
`F_{-1} = −3C_{-1}C_{-2}` [C3]. Then `E_7 = −F_{-2}` reads

```
3 C_0 C_{-1}^2  -  3 C_{-2}^2  -  2*lambda  -  2*F_{-2}  =  0
```

[C4]. **This is stronger than what GGV3 prints.** The paper displays
`C_0(3C_0C_{-1}² − 3C_{-2}² − 2λ) = 2C_0F_{-2}` and therefore has to branch on
"either `C_0 = 0` or …". The spurious `C_0` factor is an artefact of their
elimination, not of the system: multiplying our relation by `C_0` recovers the
printed form exactly [C4]. (Both branches die anyway — `C_0 = 0` also has
`c_{0,-10} = 0`.)

**Step 3 — where the kill is generated.** `(a3)` supplies `F_{-1} = y^7`, a
**unit** of `K[y,y^{-1}]`. The unit group of `K[y,y^{-1}]` is `{c y^n}`, so
`C_{-1}C_{-2} = −y^7/3` being a unit forces **both factors to be monomials**
[C5]. This step needs no cap at all.

**Step 4 — the zero margin.** `(a5)` is `deg_y(C_{-k}) ≤ k+2`, so
`deg_y C_{-1} ≤ 3` and `deg_y C_{-2} ≤ 4`. Two monomials whose exponents sum to
`7` with `e₁ ≤ 3`, `e₂ ≤ 4` force `(e₁,e₂) = (3,4)` — **the unique solution, one
integer of slack in neither place** [C6]. Hence `C_{-1} = a y³`, `C_{-2} = b y⁴`,
`ab ≠ 0`. Mutation: had `(a5)` read `k+3` the caps would be `4, 5`, `4+5 = 9 > 7`,
and three solutions survive — no monomial forcing, no floor, no kill [C9].

**Step 5 — the forced floor and the contradiction.**

```
C_0 = (3 C_{-2}^2 + 2 F_{-2} + 2 lambda) / (3 C_{-1}^2)
    = b^2 y^2/a^2 + (2/3a^2)( f_8 y^2 + f_6 + f_4 y^-2 + f_2 y^-4 + lambda y^-6 )
```

y-support `{2, 0, −2, −4, −6}`: **forced floor `j_min(0) = −6`**, reached only by
the `λ` term [C7]. `(a6)` declares `C_0 = c_{0,2}y² + ⋯ + c_{0,-10}y^{-10}` with
`c_{0,-10} ≠ 0` — a **required-nonzero at depth −10**. `−10 < −6`, so
`c_{0,-10} = 0`. **KILL**, margin 4 = two chart-steps [C8].

**The schema is chart-general, the instance is not.** The same construction from
`(b4)` `Z = x³ + Z_1x + Z_0 + ⋯` reproduces GGV3's γ=2 list `E_1..E_13` term for
term, including the `−λZ_1` of `E_13` and `Z_{-11}` occurring only in `E_8` and
`E_13` [C11]; its required-nonzero sits at a **different** slot — `(b6)`
`e_{-10} ≠ 0`, `(b5)` `c_{-1,1} ≠ 0` [C12]. Three slots across two charts,
one predicate. `ENDPOINT_CONTRACT.md` §3's cross-check claim is confirmed.

> **Evidence boundary.** `(a1)`–`(a6)` are **asserted** in GGV3 — *"We do not
> provide proofs for this first part"* (tex:1716). What is proved here is: given
> `(a1)`–`(a6)`, the kill follows, and follows from a derivation rather than a
> transcription. This is a **reproduction of a published kill from its stated
> premises**, not an independent proof that `(50,75)` is dead. That `(50,75)` is
> dead remains **CITATION-LEVEL** (Moh, via GGV5, via `moh_discards.py`).

---

## 3. Why our transferred machinery cannot see it

Every condition the transferred §8 imposes is **closed**: five polynomial
equations (`g₁, g₂, g₃, □, (*)`) and eleven inequalities (six degree caps, five
order floors) [D1]. A degree cap and an order floor are both satisfied by *more*
vanishing, never less. No boolean combination of them can require a coefficient
to be nonzero.

**The decisive comparison — one predicate, two inputs.** Run
`ENDPOINT_CONTRACT.md` §2's kill predicate

```
KILL  <=>  exists (s,j) in required_nonzero  with  j < forced_floor[s]
```

on both sides [D2]:

| input | `required_nonzero` | `forced_floor` | predicate returns |
|---|---|---|---|
| GGV3 γ=3 contract | `[(-1,3), (-2,4), (0,-10)]` | `{0: -6, -1: 3, -2: 4}` | **`[(0,-10)]`** |
| our transferred class-row data | **`[]`** | `{0: 9}` | `[]` |

**The predicate is identical. The missing datum is the primitivity list.**

And it cannot fire for a structural reason, not an oversight. Grant our side its
two strongest coefficient facts — `yplace_transfer` §E6's `e = γy⁹` and
`B = y²¹/γ`, both **forced monomials**. A forced monomial marks every slot above
the bottom as `forbidden` and the bottom one as nonzero, so its nonzero slot sits
**exactly at** the forced floor (`9 = 9`, `21 = 21`), never below it, and the
predicate still returns `[]` [D4]. **A required-nonzero can only fire when it is
DEEPER than the floor — and that is a fact no forced-monomial statement can
produce.**

**Mechanical confirmation**, in the repo's "read the module, not its output"
discipline: `yplace_transfer.py` references `(a6)`, `c_{0,-10}`, `e_{-10}`,
`required-nonzero` and `forced floor` **zero times each**. Its one use of the
word *primitive* is check A3 — *"`R = x^t C` is primitive … it is no `d`-th power
for `d ≥ 2`"* — a statement about `R` not being a proper power, **not** a
statement that a named window coefficient at a named depth is nonzero [D5].
Different notion, and the only one we carry.

> **Hence no contradiction.** The witnesses are points of a strictly weaker
> condition set than "is a germ". `yplace_transfer` said exactly this ("the §8
> witnesses are points of the §8.1 reduced system with the caps, **not germs**").
> What `(50,75)` adds is **which** missing ingredient does the work — and it is
> the same one `PROOF` §7.4(c) names and `YPLACE_TRANSFER.md` §8 lead #2 opens
> and leaves untouched: **leading-coefficient non-vanishing** [D4].

That is the value of the external control. Before this, "the missing ingredient
is leading coefficients" was our own inference from the shape of §7.4(c). Now
there is a case with a **published answer** where that is exactly and only what
killed it.

---

## 4. Mutation controls

**#1, the mandatory one: the killer must not kill `(72,108)`.** `(72,108)` dies
at Cor 8.5 by a pure **degree** count: `deg(μt³q) = 3 + deg q = 7 > 6`, attained
uniquely, for every `(deg A, z)` with `z ∈ [2,6]` — 200 pairs, zero survivors
[E1]. That kill consumes **no primitivity input whatsoever**. Removing `(a6)` —
which is what our transfer does — leaves `(72,108)` dead. The diagnosis does not
damage the closed case.

**#2, non-vacuity.** The diagnosis does not predict that everything survives:
feeding a class row's `Q_Π = 1` (`deg q = 0`) into the same count gives `3 ≤ 6`
and **25** surviving pairs, against `(72,108)`'s **0** [E2]. The same closed
machinery already separates the two cases.

**#3, the γ-chart route is not even runnable at `(72,108)`.**
`gamma_from_corner` (43 checks, calibrated on 28 published GGV1 data points)
gives at `(5,20)` **one** surviving branch `f = (4,16)`, `A_0' = (1,0)`,
`γ ∈ {2,3,4}`; at `(8,28)` one branch `f = (6,21)` with **six** admissible γ
[E3]. GGV3's two-chart argument has no `(8,28)` counterpart to be spuriously
fired.

**#4, the class corners are NOT interchangeable** — so the `(50,75)` kill must
**not** be assumed to spread across the eight [E4]:

| corner | surviving branches | shape |
|---|---|---|
| `(5,20)` | **1** — `f=(4,16)`, `A_0'=(1,0)`, `d=3` | `γ ∈ {2,3,4}` — the GGV3 §5 shape |
| `(8,32)` | **0** — no branch survives the corner conditions at all | — |
| `(9,36)` | **3** — `A_0' = (4,1), (2,1), (1,0)` | γ-sets `∅`, `{3..7}`, `{2,3,4}` |
| `(10,40)` | **1** — `A_0' = (2,0)` | `γ ∈ {2,4,5,6,7,8}` — six charts |

Only `(5,20)` has the unique-branch / `A_0'=(1,0)` / three-γ shape GGV3 §5
analyses. At `(9,36)` one of the three branches (`f=(7,28)`, `A_0'=(1,0)`,
`γ ∈ {2,3,4}`) is chart-identical in *shape* to `(5,20)`'s — that is a lead, not
a result, because the other surviving branch would also have to be handled.

**#5, the honest negative, stated so it is not inferred from absence.** Nothing
here shows the `(50,75)` witnesses fail a test we can **run**. The witnesses live
in `PROOF`'s D-transform chart (`t=4, κ=2, C=y`); `(a1)`–`(a6)` live in GGV3's
γ-reduced chart, three automorphisms downstream. **No in-repo map carries a point
of one chart to a point of the other**, so "these four points violate `(a6)`" is
**INFERRED** from "`(50,75)` is dead", not computed [E5]. Building that map is
the concrete deliverable this file identifies.

---

## 5. Does it transfer to `(75,125)`?

### The γ layer transfers verbatim

`F_2(2,3)/75`, `F_3(3,2)/75` and `F_2(3,5)/125` carry the **same** `A_0 = (5,20)`
in the atlas, and the γ layer is a function of `A_0` **alone** — GGV1's
conditions (5)–(9) read only `(u,v)`. So all three get the same unique branch
`f=(4,16)`, the same `A_0' = (1,0)`, the same `d = gcd(f₁−1,f₂−1) = 3`, the same
bound `γ ≤ 4` and the same admissible set `{2,3,4}` — **including the same
standing γ = 4 obligation** [F1]. Whatever closes γ at `(50,75)` closes it at
`(75,125)`, and vice versa. That is a real, if narrow, transfer.

### The kill layer does not

`polygon_reduction.case_f2(1)` is labelled `(75,125)` and its **computed**
reduced polygons give `deg P₁ = 15`, `deg Q₁ = 25` — **not** `(10,15)` [F2].
Every one of `(a1)`–`(a6)` is a statement about a reduced pair of degrees
`(10,15)`.

**The E-depth law (PROVED).** In `(Z^m)_{-k}` with `Z = x^c + ⋯` the only term
linear in a deep coefficient is `m·x^{c(m-1)}·Z_{-j}`, so the deepest index
reached is `K + c(m−1)` on the `P`-side and `L + c(n−1)` on the `Q`-side.
Equating them (the paper's own "`Z_{-7}` / `Z_{-11}` is the lowest coefficient
which appears") gives

```
K - L  =  c * (n - m).
```

Both published charts confirm: γ=3 (`c=2`): `5−3 = 2 = 2·1`, deepest `7 = 5+2`;
γ=2 (`c=3`): `8−5 = 3 = 3·1`, deepest `11 = 8+3` [F3].

At `(m,n) = (3,5)` the law forces `K − L = 2c`, and the `P`-side family becomes
`(Z³)_{-k}`, not `(Z²)_{-k}` [F4]. **So the `(75,125)` E-system is a different
system, not the `(50,75)` one reindexed** — and in particular there is no derived
counterpart of the collapse `E_6 ↦ 3C_{-1}C_{-2}`, which is what generates the
entire kill.

*(Scope: the `(a1)` shape at `(3,5)` — `P = C^m`, `Q = C^n + λC^{m-n} + F` — is
**INFERRED** from one anchor, `(2,3)`, verified there in both charts. The
consequence is robust regardless: `m = 3` alone changes the P-side family.)*

### Verdict, and what it would take

**The killer does NOT transfer as it stands.** What transfers is (i) the γ layer
verbatim, and (ii) the **schema** "required-nonzero below the forced floor".
What is missing is exactly three derivations at `(m,n) = (3,5)` [F5]:

1. **`(a3)`** — the leading forcing term `F_{-1} = y^7`. This is the **unit** that
   makes `C_{-1}, C_{-2}` monomials at all; without it Step 3 above does not
   start. *(Note the derivable-looking pattern across the two charts: `(a2)`'s
   `μ y^{2γ+γ-δ-1}` gives `y^6` at `(γ,δ)=(3,2)` and `y^2` at `(2,3)`, matching
   both; and the `F` exponent is that `+1` in both. Two data points — **not**
   asserted here as a law.)*
2. **`(a5)`** — the cap law `deg_y(C_{-k}) ≤ k+2`. This is precisely **Step 2 of
   the γ-window compiler**, which `SESSION_HANDOFF.md` records as NOT STARTED.
3. **`(a6)`** — the primitivity depth `−10`. There is **no in-repo derivation of
   the required-nonzero slot at any corner**; `bigrade_annotator.build_R1` takes
   it as contract *input* and computes only the floor against it.

Two of the three are not derived at `(2,3)` either — GGV3 asserts them. **So the
gap is not a transcription gap; it is the compiler.** Building it at `(2,3)`,
where the answer is published and the margins (`3+4=7`, `−10 < −6`) are known,
is now the cheapest available calibration for it — and that is the first time
this campaign has had one in the monomial regime.

### The one integer that already discriminates the rows

The shipped atlas separates them by the same pair: `G3`'s gate `λ ≥ m` **PASSES**
at `F_2(2,3)/75` (`λ = 2`, `m = 2`) and **FAILS** at `F_2(3,5)/125` (`λ = 2`,
`m = 3`) [F6]. Same corner, same `λ`; only `(m,n)` differs. `YPLACE_TRANSFER`'s
scope note (`"(75,125) = (3,5,4) has a = 3 … none of this applies"`),
`weight_free_transfer`'s rank 18 vs 19, `f2_tower`'s period 17 → 29, the reduced
degrees `(10,15)` vs `(15,25)`, and the E-depth law `c` vs `2c` are **five
readings of one fact**: `(m,n)` is `(2,3)` there and `(3,5)` here.

---

## 6. Status, honestly

| item | status |
|---|---|
| `(50,75)` = `F_2(2,3)/75` = `polygon_reduction.case_f2(0)`, corner `(5,20)`, `b₀=4a₀`, one of the eight `(2,3)` class rows | **EXACT-CHECKED**, three independent sources agreeing (paper text, atlas, computed reduction) |
| `deg P₁ = 10`, `deg Q₁ = 15` from the **computed** reduced polygon, matching GGV3's published values | **EXACT-CHECKED** (cross-check, not a fit) |
| the four §8 witnesses satisfy the §8.1 `k=0` system at `(50,75)`'s own caps, residual 0 | **EXACT-CHECKED**, rebuilt from `PROOF` (8.1.1) as printed |
| `(5,20)`'s caps are the tightest of the four class corners | **EXACT-CHECKED** at `(5,20)`; the other three slopes are `YPLACE_TRANSFER` §D5's, **INFERRED** there (no computed polygon) |
| `F_3(3,2)/75` is the `P↔Q` swap and the transferred objects are bit-identical, so no asymmetry exists | **EXACT-CHECKED** |
| `E_1..E_8` and `E_1..E_13` derived from the `Z`-series reproduce GGV3 term for term | **PROVED** (ring identities) |
| `F_{-1} = −3C_{-1}C_{-2}` and `3C_0C_{-1}² = 3C_{-2}² + 2F_{-2} + 2λ` derived by elimination; the latter **without** GGV3's spurious `C_0` factor | **PROVED** from `(a1)`–`(a5)` |
| monomial forcing by the unit `y^7`, then `(3,4)` uniquely by `(a5)`; zero margin | **PROVED**, with the `k+3` mutation control |
| forced floor `−6`; `(a6)`'s required `−10`; KILL with margin 4 | **PROVED** given `(a1)`–`(a6)` |
| `(a1)`–`(a6)` themselves; "Moh ruled `(50,75)` out" | **CITATION-LEVEL.** GGV3 declines to prove them; we have not read [M] |
| our transferred §8 imposes only closed conditions; its `required_nonzero` list is empty; the contract predicate cannot fire | **EXACT-CHECKED** (the inventory is complete against `PROOF` §8.1–8.5) |
| **therefore the tension is not a contradiction** | **established** |
| "these four witnesses violate `(a6)`" | **INFERRED** from `(50,75)` being dead. No in-repo map carries a point of `PROOF`'s chart to a point of GGV3's γ-chart. **NOT COMPUTED.** |
| the γ layer is identical at `(50,75)`, `(75,50)` and `(75,125)` | **EXACT-CHECKED** via `gamma_from_corner` |
| reduced pair `(15,25)` at `(75,125)`; E-depth law `K−L = c(n−m)` | **EXACT-CHECKED** / **PROVED** |
| the `(a1)` shape at `(3,5)` | **INFERRED**, one anchor |
| the killer transfers to `(75,125)` | **NO** — and the three missing derivations are named |

### Negatives, plainly

* **No case is closed by this file, and the frontier does not move.** `(50,75)`
  and `(75,50)` were already settled (`moh_discards.py`); the other six `(2,3)`
  class rows are **not** thereby settled — mutation control #4 shows their corner
  data are genuinely different, and `(8,32)` has no surviving branch at all.
* **The resolution is a statement about our machinery, not about the class.** It
  explains why the §8 witnesses exist; it does not produce a kill for anything.
* **The one thing that would have been a bigger finding is absent.** If `(50,75)`
  had differed from the generic class row in a way that already killed the
  witnesses, we would have had a discriminator to test at `(75,125)`. It does
  not: `(5,20)` is where `yplace_transfer` instantiated the caps in the first
  place.

### What would unblock it, precisely

1. **A chart bridge.** A map carrying a point of `PROOF`'s D-transform chart
   (`t=4, κ=2, C=y`) to a point of GGV3's γ-reduced chart. With it, "witness ⟹
   violates `(a6)`" becomes computable, and the `(a6)` obstruction becomes a
   *runnable* test at every class row rather than a citation. This is the single
   highest-value item this file opens.
2. **`(a5)` at `(2,3)`, derived.** Step 2 of the γ-window compiler, now with a
   published answer to calibrate against and a known zero margin (`3+4=7`) that
   any correct derivation must reproduce exactly.
3. **The required-nonzero slot, derived.** Nothing in the repo derives `−10` (or
   `(b5)`'s `−20`, or `(b6)`'s `−10`) from corner data. Three published slots
   across two charts is a real, if small, calibration set.

---

## Files

| file | role |
|---|---|
| `MOH_CONTROL_50_75.md` | this writeup |
| `moh_control_50_75.py` | the checker — 48/48, `--quiet` exit 0, ~12 s |
| `yplace_transfer.py` / `YPLACE_TRANSFER.md` | half one of the tension: the §8 witnesses |
| `moh_discards.py` | half two: `F_2(2,3)/75` **is** `(50,75)`, settled in the literature |
| `ENDPOINT_CONTRACT.md` | the three-status contract and the kill predicate used in §3 |
| `f2_tower.py` / `F2_TOWER.md` | the *replay* this file replaces with a derivation |
| `gamma_from_corner.py` | the γ layer, and mutation controls #3 and #4 |
| `polygon_reduction.py` | the computed reduced polygons at `(50,75)` and `(75,125)` |
| `paper_src/1406.0886_GGV3.tex` | §5, tex:1708–2070 — `(a1)`–`(a6)`, `(b1)`–`(b6)`, both `E`-lists |
