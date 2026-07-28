# Makar-Limanov 2025, Theorem 2 and the Lemma on divisibility, vs. the 34-row census

**Lane:** literature-import / filter evaluation. **Date:** 2026-07-28.
**New files (2, both intended):** `makar_limanov_filter.py`, `MAKAR_LIMANOV_FILTER.md`.
Read-only over every existing artifact. No concurrent-lane file touched.

```
python makar_limanov_filter.py           # full narrative + census table
python makar_limanov_filter.py --quiet   # SUMMARY + exit code only   -> 49/49, exit 0
```

---

## 0. Verdict, up front

**Neither obstruction kills anything.**

| | |
|---|---|
| Theorem 2, rows refuted | **0 / 34** |
| Theorem 2, rows refuted among the 10 already-settled calibration rows | **0 / 10** |
| `F_2(3,5)/125` (the unique open row **at** the 125 bound) | **SURVIVES** |
| `(8,28)/(3,2)/108` (the one row with a full reduction *and* a proved death) | **SURVIVES T2** — i.e. T2 is blind to the one death we can check it against |
| Lemma on divisibility, applicability to our corners | **DOES NOT TRANSFER** — its conclusion is exhibited **false** under `[P,Q] = x²` |
| Rows where observed admissible-edge count exceeds the predicted bound | **0** (no contradiction, no bug, no kill) |

This is a **clean negative**, and it is calibrated: §4 shows the filter is live
(it kills synthetic corners through the identical code path), and §3 shows that
the one reading of T2 which *would* have killed 14 rows — including the flagship
— is refuted by GGV5's own normal form on **all 34 rows**.

**The 125 bound does not move.** The one thing T2 leaves behind is a
falsifiable *shape* constraint at the `(5,20)` corner (§5.3), not a kill.

---

## 1. The source, and the two results verbatim

**L. Makar-Limanov, "On the shape of a counterexample to the two-dimensional
Jacobian conjecture", _Serdica Math. J._ **51** (2025) 299–314,**
doi `10.55630/serdica.2025.51.299-314`.
Open galley: `https://serdica.math.bas.bg/index.php/serdica/article/download/300/153`.

**Full text obtained and read.** This closes the single open item recorded in
`ML_RESTRICTION.md` §6 ("Obtain the Serdica 51(3–4) PDF and reconcile the
verbatim lemma statement"). See §6.1 for the reconciliation.

### (T2) Theorem 2 — p. 311, verbatim

> **Theorem 2.** If `f` is reduced and `v₀ = ∏_{i=1}^{k} p_i^{δ_i} (a,b)` where
> all `p_i` are prime numbers and `gcd(a,b) = 1` then `N(f)` has at most
> `Σ_i δ_i − 1` admissible edges.

### (DIV) Lemma on divisibility — p. 305, verbatim

> **Lemma on divisibility.** If `J(σ,τ) = σ` where `σ, τ` are `w` homogeneous,
> `w(σ) > 0`, `w(τ) > 0` and `σ` is not a monomial then `w(σ)` doesn't divide
> `w(τ)`.

### The four definitions T2 depends on, verbatim

| term | ML's text | page |
|---|---|---|
| **reduced** | "From now on we will be looking at counterexamples where `supp(f)` has such a monomial [`x^{d_x} y^{d_y}`, `d_x = deg_x(f)`, `d_y = deg_y(f)`]. Since we can make an automorphism `x → y, y → −x` we will also assume that `d_x = m < d_y = n`. Additionally, if `w_{1,0}(f) = x^m p(y)` we can make a substitution `y → y − c` such that the order of `p(y−c)` is larger than `m` … Let us call such an `f` **reduced**." | 304–305 |
| **the rectangle** | "Since `f` is reduced … points of the Newton polygon `N(f)` of `f` belong to the rectangle with vertices `(0,0), (m,0), (m,n), (0,n)`." | 307 |
| **right edges** | "There are two edges `e₀` and `e₁` with the vertex `(m,n)` … extension of `e₀` intersects the `y` axis and extension of `e₁` intersects the `x` axis. We will call them left and right leading edges accordingly." / "The boundary of the Newton polygon `N(f)` of `f` consists of right edges `e₁, e₂, …, e_k` and left edges `e₀, e₋₁, …, e₋ₗ`." | 307 |
| **admissible / the bisectrix** | "We will be interested in the right edges which have at least one vertex **above the bisectrix of the first quadrant**. We will call this edges **admissible**." | 308 |
| **`v₀`** | "For a monomial `μ` denote by `\|μ\|` its degree vector and by `v₀` the vertex of `N(f)` corresponding to `f̂`." (and p. 302: "If `a ∈ C[x,y]` denote by `â` its leading monomial.") | 310 |

For reduced `f` the `(1,1)`-leading monomial is `x^m y^n`, so

> **`v₀ = (deg_x f, deg_y f) = (m,n)` with `m < n`.**

`Σ_i δ_i = Ω(gcd(v₀))`, the number of prime factors **with multiplicity**.

---

## 2. Two readings declared (the paper is ambiguous, and it matters)

### 2.1 The quantifier on `s`

T2 as *literally* quantified would read `0 ≤ −1` whenever `gcd(v₀) = 1` and
`N(f)` happens to have no admissible edge — i.e. it would by itself refute every
counterexample with coprime corner coordinates. That is **not** what the proof
establishes. The proof opens

> "Now we can finish the proof. For **the last admissible edge** `d_s` is at
> least the denominator of `γ = deg_y(g)/deg_y(f)` …"

i.e. it assumes `s ≥ 1` admissible edges exist and bounds `s`. **The reading
used here is the proof's:**

> **(T2′)** if `N(f)` has `s ≥ 1` admissible edges then `s ≤ Ω(gcd(v₀)) − 1`.

**This costs nothing**, because `s ≥ 1` is *forced*: `e₁` is a right edge
carrying the vertex `v₀ = (m,n)` (p. 307), and reducedness gives `m < n`, so
`v₀` lies strictly above the bisectrix and `e₁` is admissible. Hence

> **(T2-KILL)** a reduced `f` is impossible unless **`Ω(gcd(v₀)) ≥ 2`**.

That single inequality is the *entire* content of T2 as a filter
(`makar_limanov_filter.py` check `A5`). It is checked, with the floor `s ≥ 1`
implemented separately (`t2_floor`) and a **mutation control** (`A3m`) showing
that dropping the `m < n` half of "reduced" collapses the floor to 0 and makes
T2 unable to kill anything at all.

### 2.2 "above the bisectrix"

Read as **strictly** above (`y > x`). Immaterial for `v₀`, where `m < n` is
strict.

---

## 3. The load-bearing step: what `v₀` **is** on our rows

A filter of this kind is won or lost on the identification of `v₀`. Ours is
**not** a guess — it is forced by a verbatim sentence of GGV5.

### 3.1 GGV5's normal form *is* ML's reduced form

`paper_src/1708.07936_GGV5.tex:250`, verbatim:

> "If this conjecture is false, then there exist `P,Q ∈ L` such that
> `[P,Q] = K^×`, and there exist `m,n,a,b ∈ ℕ`, such that `m,n > 1` are coprime,
> `a < b`, **the support of `P` is contained in the rectangle with vertices
> `{(0,0), m(a,0), m(a,b), m(0,b)}`**, the support of `Q` is contained in the
> rectangle with vertices `{(0,0), n(a,0), n(a,b), n(0,b)}`, **the point
> `m(a,b)` is in the support of `P`** and the point `n(a,b)` is in the support
> of `Q`. Note that `deg(P) = m(a+b)` and `deg(Q) = n(a+b)`."

and `:398-400`: "`A₀ = (1/m) en₁₀(P)` … This point `A₀` corresponds to `(a,b)`
in the introduction."

Pure arithmetic on that sentence gives, with `A₀ = (a₀,b₀)`:

1. `supp(P) ⊆` rectangle ⟹ `deg_x P ≤ m a₀`, `deg_y P ≤ m b₀`; and
   `m(a₀,b₀) ∈ supp(P)` ⟹ the reverse. Hence
   **`(deg_x P, deg_y P) = m·A₀`** and **`(deg_x Q, deg_y Q) = n·A₀`**.
2. `x^{deg_x P} y^{deg_y P} ∈ supp(P)` — **exactly ML's defining property of
   "reduced"**. GGV5's normal form *is* ML's reduced form (ML's extra `y → y−c`
   shift moves no vertex).
3. `a < b ⟹ m a₀ < m b₀`: ML's `m < n` normalisation holds for **both** `P` and `Q`.
4. `m,n > 1` coprime ⟹ `γ = n/m` is neither an integer nor the reciprocal of
   one — **ML's own normalisation of `γ`**, satisfied verbatim.

> ### **`v₀(P) = m·A₀`  and  `v₀(Q) = n·A₀`.  NOT `A₀`.**

**Grade: PROVED** (arithmetic on a pinned published sentence), cross-checked
three ways in the checker: `B1` reproduces `max_deg` on all 34 rows from
`deg P = m(a+b)`; `B2`/`B2b` verify ML's hypotheses row by row; `B3` verifies
**ML's own Lemma on similarity** — `γ·v₀(P) = v₀(Q)`, integral — on all 34 rows.

### 3.2 The naive reading, and why it is a trap

The obvious thing to do is set `v₀ := A₀`. **It is wrong on all 34 rows**, and
it is wrong in the dangerous direction.

| refutation | result |
|---|---|
| Reducedness forces `deg f = (v₀)_x + (v₀)_y`. With `v₀ := A₀` that gives `deg f = a₀+b₀`, but GGV5 gives `deg P = m(a₀+b₀)` with `m ≥ 2`. | fails on **34/34** (`B4`) |
| ML's Lemma on similarity needs `denom(γ) = min(m,n)` to divide `gcd(v₀)`. With `v₀ := A₀`, it doesn't. | fails on **21/34**, including the flagship (`B4b`) |

And what the naive reading would have *claimed*:

> **MISFIRE CONTROL (`B5`, `B5b`).** `v₀ := A₀` "kills" **14** rows —
> `F_2(2,3)/75`, **`F_2(3,5)/125`**, `F_3(3,2)/75`, `F_7(2,7)/147`,
> `F_8(3,7)/147`, `F_9(2,3)/84`, `F_9(3,5)/140`, `F_11(2,5)/140`,
> `F_17(2,3)/99`, `(7,35)/(2,3)/126`, `(7,42)/(3,2)/147`, `(7,42)/(2,3)/147`,
> `(11,33)/(2,3)/132`, `(12,33)/(2,3)/135` — i.e. it kills the **flagship open
> row** while **sparing `(8,28)/(3,2)/108`**, the one row this campaign actually
> closed. **Anti-correlated with truth on the two rows whose answer we know.**
> All 14 are artefacts. **Do not cite them.**

A second way to manufacture the same fake kills: use `ω` (distinct primes)
instead of `Ω` (with multiplicity). Mutating that one token in the checker
produces **15** kills including both `(8,28)/(3,2)/108` and `F_2(3,5)/125`, and
trips 16 checks (verified out-of-tree). `A1m` guards it.

---

## 4. The census

`v₀(P) = m·A₀`, `v₀(Q) = n·A₀`; `Ω = Ω(gcd(v₀))`; `bnd = Ω − 1` is T2's
predicted maximum admissible-edge count; the forced floor is `1` everywhere.

| row | `A₀` | `v₀(P)` | `Ω` | bnd | `v₀(Q)` | `Ω` | bnd | T2 |
|---|---|---|---|---|---|---|---|---|
| `F_1(3,4)/64` | (4,12) | (12,36) | 3 | 2 | (16,48) | 4 | 3 | alive |
| `F_1(5,7)/112` | (4,12) | (20,60) | 3 | 2 | (28,84) | 3 | 2 | alive |
| `F_2(2,3)/75` | (5,20) | (10,40) | 2 | 1 | (15,60) | 2 | 1 | alive |
| **`F_2(3,5)/125`** | **(5,20)** | **(15,60)** | **2** | **1** | **(25,100)** | **2** | **1** | **alive** |
| `F_3(3,2)/75` | (5,20) | (15,60) | 2 | 1 | (10,40) | 2 | 1 | alive |
| `F_7(2,7)/147` | (6,15) | (12,30) | 2 | 1 | (42,105) | 2 | 1 | alive |
| `F_8(3,7)/147` | (6,15) | (18,45) | 2 | 1 | (42,105) | 2 | 1 | alive |
| `F_9(2,3)/84` | (7,21) | (14,42) | 2 | 1 | (21,63) | 2 | 1 | alive |
| `F_9(3,5)/140` | (7,21) | (21,63) | 2 | 1 | (35,105) | 2 | 1 | alive |
| `F_11(2,5)/140` | (7,21) | (14,42) | 2 | 1 | (35,105) | 2 | 1 | alive |
| `F_17(2,3)/99` | (9,24) | (18,48) | 2 | 1 | (27,72) | 2 | 1 | alive |
| `F_22(2,3)/96` | (8,24) | (16,48) | 4 | 3 | (24,72) | 4 | 3 | alive |
| `F_24(3,4)/128` | (8,24) | (24,72) | 4 | 3 | (32,96) | 5 | 4 | alive |
| `(7,35)/(2,3)/126` | (7,35) | (14,70) | 2 | 1 | (21,105) | 2 | 1 | alive |
| `(7,42)/(3,2)/147` | (7,42) | (21,126) | 2 | 1 | (14,84) | 2 | 1 | alive |
| `(7,42)/(2,3)/147` | (7,42) | (14,84) | 2 | 1 | (21,126) | 2 | 1 | alive |
| `(8,28)/(3,4)/144` | (8,28) | (24,84) | 3 | 2 | (32,112) | 4 | 3 | alive |
| **`(8,28)/(3,2)/108`** | **(8,28)** | **(24,84)** | **3** | **2** | **(16,56)** | **3** | **2** | **alive** |
| `(9,36)/(3,2)/135` | (9,36) | (27,108) | 3 | 2 | (18,72) | 3 | 2 | alive |
| `(9,36)/(2,3)/135` | (9,36) | (18,72) | 3 | 2 | (27,108) | 3 | 2 | alive |
| `(11,33)/(2,3)/132` | (11,33) | (22,66) | 2 | 1 | (33,99) | 2 | 1 | alive |
| `(12,33)/(2,3)/135` | (12,33) | (24,66) | 2 | 1 | (36,99) | 2 | 1 | alive |
| `(8,32)/(3,2)/120` | (8,32) | (24,96) | 4 | 3 | (16,64) | 4 | 3 | alive |
| `(8,40)/(3,2)/144` | (8,40) | (24,120) | 4 | 3 | (16,80) | 4 | 3 | alive |
| `(9,27)/(2,3)/108` | (9,27) | (18,54) | 3 | 2 | (27,81) | 3 | 2 | alive |
| `(9,36)/(2,3)/135` | (9,36) | (18,72) | 3 | 2 | (27,108) | 3 | 2 | alive |
| `(10,40)/(3,2)/150` ×2 | (10,40) | (30,120) | 3 | 2 | (20,80) | 3 | 2 | alive |
| `(12,30)/(3,2)/126` | (12,30) | (36,90) | 3 | 2 | (24,60) | 3 | 2 | alive |
| `(12,36)/(2,3)/144` ×4 | (12,36) | (24,72) | 4 | 3 | (36,108) | 4 | 3 | alive |
| `(12,36)/(3,2)/144` | (12,36) | (36,108) | 4 | 3 | (24,72) | 4 | 3 | alive |

**0 kills.** The reason is structural, not accidental (`C3`, `C3b`):

```
Ω(gcd(v₀(P))) = Ω(m) + Ω(gcd(A₀)),    m > 1 (GGV5),    gcd(A₀) ∈ {3,…,12}
```

so `Ω ≥ 2` **automatically** on every row. T2's kill condition requires
`gcd(a,b) = 1` together with a prime multiplier, and **no published corner in
the `max_deg ≤ 150` range has `gcd(a,b) = 1`.**

> **MUTATION CONTROL (`C4`, `C4b`).** The *same* code path **does** kill
> synthetic rows with `gcd(A₀) = 1` and a prime multiplier — `A₀ ∈ {(1,4),
> (3,8), (2,5), (5,12)}` are all refuted — and correctly spares the boundary
> case `gcd(A₀) = 1` with both multipliers composite. **The 0/34 is a fact
> about the census, not a dead filter.**

### 4.1 A general corollary worth keeping (`C5`) — **PROVED**

> T2 + GGV5's normal form ⟹ **for any counterexample, if `gcd(a,b) = 1` then
> both `m` and `n` must be composite.**

Non-vacuous in general; unused here, because every corner in the table already
has `gcd(a,b) ≥ 3`. It is the *only* place T2 can ever bite in this framework,
so it is where to point T2 if the census is ever extended past `max_deg = 150`.

---

## 5. Calibration

The task's standing rule: *a filter that misfires on a known answer is not
usable.* Three calibrations were run.

### 5.1 Against the 10 rows settled below the 125 bound — `D1`, `D2`

The calibration set is `gghv_sub125.py`'s computed partition (nine discarded
upstream, one closed by this campaign).

> **T2 kills 0 of the 10 known-dead rows.** It reproduces **no** death that is
> already established. **It therefore has no demonstrated discriminating power
> on this census.**

This is the honest headline. A filter that cannot recover a single known death
is not evidence about the open rows either way.

### 5.2 Against `(8,28)/(3,2)/108` — the sharpest single test — `D3`

This is the one row with a **published full reduction** *and* a **proved
death**. T2 sees:

```
v₀(P) = (24,84)   Ω(gcd) = Ω(12) = 3   bound 2
v₀(Q) = (16,56)   Ω(gcd) = Ω(8)  = 3   bound 2      ⟹ ALIVE
```

**The one death we can check T2 against is invisible to it.**

*By-product (`D3b`), worth recording separately:* this **derives
original-coordinate corner data the repo never wrote down** — for the `(72,108)`
case, `(deg_x P, deg_y P) = (24,84)` (deg 108) and `(deg_x Q, deg_y Q) = (16,56)`
(deg 72). Previously the repo carried only `A₀ = (8,28)` and the post-chart
reduced polygons; the original-coordinate corner existed only as a rule.
**Grade: PROVED** from GGV5:250.

### 5.3 Observed vs. predicted admissible-edge counts — `D4`, `D4b`, `D4m`

Where the repo has an original-coordinate polygon at all, it is the unit polygon
of `passport_75_125.py` rule **r1**:
`Δ = {(0,0), (1,0), (a₀,b₀), (0,c)}`, `c = b₀ − ⌊(b₀−1)/a₀⌋·a₀`, with
`N(P) = m·Δ`, `N(Q) = n·Δ`. Written out in print in this repo only at `(5,20)`:
`Δ = {(0,0),(1,0),(5,20),(0,5)}`.

> **Grade: INFERRED.** The repo itself flags r1's chart exponent as unpublished,
> and the "+ published extra corners" clause means the vertex list may be
> **incomplete**. So the observed counts below are **lower bounds** on the true
> admissible-edge count. That direction is exactly the one needed to test for a
> *contradiction* with T2's upper bound, so the test is sound even though the
> number is not authoritative.

| row | poly | `N(·)` hull | observed | predicted |
|---|---|---|---|---|
| `F_2(2,3)/75` | P | (0,0),(2,0),(10,40),(0,10) | 1 | 1 |
| `F_2(2,3)/75` | Q | (0,0),(3,0),(15,60),(0,15) | 1 | 1 |
| **`F_2(3,5)/125`** | **P** | **(0,0),(3,0),(15,60),(0,15)** | **1** | **1** |
| **`F_2(3,5)/125`** | **Q** | **(0,0),(5,0),(25,100),(0,25)** | **1** | **1** |
| `F_3(3,2)/75` | P | (0,0),(3,0),(15,60),(0,15) | 1 | 1 |
| `F_3(3,2)/75` | Q | (0,0),(2,0),(10,40),(0,10) | 1 | 1 |
| `(8,28)/(3,4)/144` | P | (0,0),(3,0),(24,84),(0,12) | 1 | 2 |
| `(8,28)/(3,4)/144` | Q | (0,0),(4,0),(32,112),(0,16) | 1 | 3 |
| `(8,28)/(3,2)/108` | P | (0,0),(3,0),(24,84),(0,12) | 1 | 2 |
| `(8,28)/(3,2)/108` | Q | (0,0),(2,0),(16,56),(0,8) | 1 | 2 |

**No row has observed > predicted.** Nothing is contradicted; nothing is killed;
no bug is exposed.

> **MUTATION CONTROL (`D4m`).** Inserting one genuine extra hull vertex above
> the bisectrix into the `(5,20)` polygon raises the observed count to **2**,
> which **exceeds** the predicted **1** and is flagged. The `D4` pass is a real
> test, not a detector that always returns 1.

### 5.4 The one residue T2 leaves at the 125 corner — `D5b`, `D5c`

At `A₀ = (5,20)` the predicted bound equals the forced floor:

> **`bound = 1 = floor` ⟹ `N(P)` and `N(Q)` have EXACTLY ONE admissible right
> edge.**

Concretely, for `F_2(3,5)/125`: **no right-boundary vertex of `N(P)` other than
`(15,60)` may lie strictly above the bisectrix `y = x`** (and likewise `(25,100)`
for `N(Q)`). The r1 polygon satisfies this with zero slack. This is a genuine,
falsifiable shape constraint on the flagship case — the only usable output T2
produced — but it is **not** a kill, and it constrains a polygon detail nothing
in the current program is trying to vary.

---

## 6. (DIV): the Lemma on divisibility

### 6.1 Reconciliation with `ML_RESTRICTION.md` (2026-07-24) — **now closed**

That lane tested (DIV) from an external *paraphrase*, could not obtain the PDF,
and left one open item: obtain the text and reconcile. **Done.** Two findings:

1. **The paraphrase was faithful.** The verbatim statement (§1) matches
   `ML_RESTRICTION.md` §2's rendering, including the degree-bookkeeping
   consequence `w(τ) = α+β` — which the relation `J(σ,τ) = σ` *forces*, since
   `w(J(σ,τ)) = w(σ)+w(τ)−w(xy)`. So the lemma's conclusion is literally
   "`w(σ) ∤ (α+β)`". (`E1`)
2. **The missing-hypothesis finding was correct, and is not an artefact of the
   paraphrase.** ML's own weight convention is "`α, β ∈ ℤ, gcd(α,β) = 1`", which
   **admits the diagonal `(1,1)`**, and there the statement is false:

   ```
   w = (1,1),   ρ = x + y   (not a monomial),   τ = xy + y²
   J(ρ,τ) = (1)(x+2y) − (1)(y) = x + y = ρ        ✓ exactly
   w(ρ) = 1 > 0,  w(τ) = 2 > 0,  and  1 | 2       ✗ contradicts the conclusion
   ```
   **Grade: EXACT-CHECKED** (`E2`, `E2b`). The proof's first step ("`σ̂` and `τ̂`
   are algebraically dependent") is inherited from the surrounding Dixmier
   context, not from the lemma's stated hypotheses; `x` and `xy` are
   independent. The operative hypotheses are **positive + primitive +
   non-diagonal**.
3. Over 14 positive primitive **non-diagonal** weights, no counterexample exists
   (`E3`). **Grade: RECONNAISSANCE ONLY** — a finite coefficient grid, not a
   proof. Consistent with the lemma holding in the regime a Newton-polygon shape
   lemma is about.

### 6.2 Does it transfer to our setting? **No — and the reason is not a technicality**

Our corners carry **`[P,Q] = x^κ` with `κ = 2`**, not `J = 1`
(`polygon_reduction.py` R1/R2/R3; GGHV22 `2204.14178.tex:1001`; GGV3
`1406.0886_GGV3.tex:1726`). The relation available at a corner is therefore
`J(ρ,τ) = x²ρ`, and degree bookkeeping forces `w(τ) = 3α+β` rather than `α+β`.

The natural hope is that a rescaling repairs the normalisation. **It does not,
and the failure is constructive:**

> **EXACT COUNTEREXAMPLE (`E4b`).** At the **positive, primitive, non-diagonal**
> weight `w = (1,3)`:
> ```
> ρ = x³ + y        (w-homogeneous, w(ρ) = 3, NOT a monomial)
> τ = (x³y + y²)/3  (w-homogeneous, w(τ) = 6)
> J(ρ,τ) = 3x²·(x³+2y)/3 − 1·(3x²y/3) = x⁵ + x²y = x²(x³+y) = x²ρ   ✓ exactly
> w(ρ) = 3   DIVIDES   w(τ) = 6                                     ✗
> ```
> **The divisibility conclusion is FALSE under our bracket normalisation.**
> A second, independent instance at `w = (1,2)`: `ρ = x(x²+y)²`,
> `τ = (x³y + xy²)/5`, `w(ρ) = w(τ) = 5` (`E4d`).
> Dividing `τ` by `x²` does not restore `J(ρ,τ′) = ρ` (`E4c`).
>
> **MUTATION CONTROL (`E4m`):** the identical search at `κ = 0` (ML's own
> normalisation) at `w = (1,3)` and `w = (1,2)` finds **nothing**. The
> counterexamples are produced by the `x²` twist, not by the search.

**Grade: EXACT-CHECKED**, hand-verifiable in one line.

### 6.3 An independent structural reason, twist aside — `E5`, `E5b`

Even granting the untwisted lemma, its hypothesis is unreachable on our
polygons. Strictly-positive-weight **edges** of the computed reduced polygons:

| polygon | positive-weight edges |
|---|---|
| `(8,28)` sub1/sub2, `N(P)` and `N(Q)` (4 polygons) | **NONE** — the origin is a vertex, so every positive weight is maximised at a single monomial |
| `(5,20)` `(2,3)` and `(3,5)`, `N(P)` and `N(Q)` (4 polygons) | **exactly one each, at the DIAGONAL weight `(1,1)`** |

So at `(8,28)` there is no positive-weight edge for (DIV) to speak about at all
(reproducing `ML_RESTRICTION.md` §3 independently), and at `(5,20)` the only one
sits at **precisely the weight class where (DIV)'s conclusion is false even at
`κ = 0`** (§6.1). **(DIV) has no bite at the 125 corner either.**

> **MUTATION CONTROL (`E5m`).** On a probe polygon with a genuine non-diagonal
> positive-weight edge, the same detector finds it.

---

## 7. Evidence ledger

| claim | grade |
|---|---|
| The verbatim text of T2, (DIV), and the four definitions | **CITATION-LEVEL** — Serdica 51 (2025) 299–314, full text read |
| `Σ_i δ_i = Ω(gcd(v₀))`; T2-as-a-filter reduces to `Ω(gcd(v₀)) ≥ 2` | **PROVED** |
| `s ≥ 1` is forced (`e₁` carries `v₀`, `m < n`) | **PROVED** from ML pp. 307–308 |
| The `s ≥ 1` reading of T2's quantifier | **DECLARED** (§2.1) — the proof's reading; the literal one is stronger than what is proved |
| `v₀(P) = m·A₀`, `v₀(Q) = n·A₀`; GGV5's normal form *is* ML's reduced form | **PROVED** from `1708.07936_GGV5.tex:250` + `:398` |
| `v₀ := A₀` is refuted on 34/34 rows; its 14 "kills" are artefacts | **PROVED** / **EXACT-CHECKED** |
| The census (`Ω`, bounds, 0 kills over 34 rows, 0 over the 10 calibration rows) | **EXACT-CHECKED** (49/49, exit 0) |
| The filter is live — synthetic corners *are* killed by the same code | **EXACT-CHECKED** |
| `(deg_x P, deg_y P) = (24,84)`, `(deg_x Q, deg_y Q) = (16,56)` for `(72,108)` | **PROVED** from GGV5:250 (new to the repo) |
| Observed admissible-edge counts (all = 1, none exceeding predicted) | **INFERRED** polygon (rule r1, repo-flagged unpublished, possibly incomplete) → the counts are **lower bounds**, which is the safe direction for a contradiction test |
| (DIV) as literally stated is false at `w = (1,1)` | **EXACT-CHECKED** |
| (DIV) holds over the non-diagonal weights scanned | **RECONNAISSANCE ONLY** — finite grid, not a proof |
| (DIV) does not transfer: its conclusion is false under `[P,Q] = x²` | **EXACT-CHECKED** (two explicit instances + mutation control) |
| No positive-weight edge at `(8,28)`; only the diagonal one at `(5,20)` | **EXACT-CHECKED** on the computed reduced polygons |

---

## 8. Blunt verdict

**Neither obstruction kills anything, and the negative is calibrated.**

* **Theorem 2 kills 0 of the 34 rows and 0 of the 10 rows whose death is already
  established.** It cannot see the one death this campaign proved. On this
  census it is structurally toothless: `Ω(gcd(v₀)) = Ω(m) + Ω(gcd(A₀)) ≥ 2`
  holds automatically because `m > 1` and every published corner in range has
  `gcd(a,b) ≥ 3`. Its only possible bite here is the corollary of §4.1, which
  needs `gcd(a,b) = 1` — a corner shape that does not occur below
  `max_deg = 150`.
* **`F_2(3,5)/125` survives.** T2 does leave a sharp shape residue there
  (exactly one admissible right edge, §5.3), but that constrains a polygon
  detail, not the coefficient system, and **the 125 bound does not move.**
* **The Lemma on divisibility does not transfer.** Under `[P,Q] = x²` its
  conclusion is demonstrably false; and independently, our reduced polygons
  carry no non-diagonal positive-weight edge for it to act on. The lane opened
  by `ML_RESTRICTION.md` on 2026-07-24 is now closed against the verbatim text,
  with the same answer it reached from the paraphrase.
* **The most valuable thing this lane produced is a trap, disarmed.** The
  obvious identification `v₀ := A₀` "kills" 14 rows including the flagship,
  spares the row we proved dead, and is refuted by GGV5's own normal form on
  every row. Anyone importing Makar-Limanov's Theorem 2 into this program
  without §3 will get a spectacular and entirely false result.

**Frontier impact: none.** No row is removed, promoted, or contradicted. The
open frontier remains **24**.
