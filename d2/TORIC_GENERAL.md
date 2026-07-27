# Is `6*W*Z = e^5` family-level? — the general power window, tested at `(3,5)`

**Answer: NO on all counts. `(72,108)`'s toric syzygy is a `(m,n,t) = (2,3,4)`
resonance, the exponent tracks `t` and not `m+n`, `(75,125)` admits no identity of
that shape at any exponent, and the `(75,125)` slice route closes too (`lam = 0`).**

**Checker:** `toric_general.py` — **38/38**, `--quiet` exit 0, ~80 s.
**Prerequisites:** `TORIC_SYZYGY.md` / `toric_syzygy.py` (the `(72,108)` identity),
`CONTACT_LEMMA.md` (the `lam >= m` gate), `PASSPORT_75_125_REPAIR.md` (the repaired
`(75,125)` inputs).
**Sources of truth:** generators from `g_system_75_125.build_gsystem` (whose `(2,3,4)`
output is checked against `g_system_75_125.published_72108`); `lam` arithmetic from
`window_functions_75_125.window_law` / `.family`. Nothing is retyped.

Notation: `(m,n) = (a,b)` are the reduced powers (linear window `S^a`, forcing window
`S^b`); `t` is the chart exponent (`l`); `D = (a-1)*t` is the spare inventory
`dm1..dmD`; `e = dm1`, `R = dm2`, `S = dm3`, `T = dm4`;
`W = e*S - R^2`, `Z = e*T - R*S`.

---

## 0. The three answers, up front

### Q1. Is the exponent `m + n`? — **NO.**

The exponent is forced by the **chart exponent `t` alone**, through a weight-divisibility
condition, and `m`/`n` do not enter it. In the `u`-grading `w(s_i) = t - i` (which is
where the whole system lives), `w(W) = 2t+4`, `w(Z) = 2t+5`, `w(e) = t+1`. So an identity
of the shape `const * W * Z = e^k` **requires**

```
4t + 9  =  k*(t+1) ,     i.e.   (t+1) | (4t+9)  =  4(t+1) + 5 ,   i.e.  (t+1) | 5
```

whose only solution with `t >= 2` is **`t = 4`, and then `k = 5` necessarily**.

At `(72,108)` the number 5 is over-determined: `m+n = 5`, `t+1 = 5`, `(4t+9)/(t+1) = 5`,
`2n-1 = 5`, `m*n-1 = 5` — **five different formulas, one data point**. They separate at
`(3,5)` (`8`, `5`, `5`, `9`, `14`) and at `t = 5,6`. Holding `(m,n) = (2,3)` fixed and
moving `t` to `3, 5, 6` **destroys the identity entirely** — same `m`, same `n`, no
identity at any exponent. That is the direct refutation of the `m+n` reading.

### Q2. Does `(3,5)` admit the identity? — **NO. Completely, not "below a bound".**

At `(75,125)` we have `t = 4` and `D = 8`, so the window-Hankel minors have `u`-weights
`12..22`, products of two of them have weights `24..44`, and `e^k` has weight `5k`. Hence
the **complete** list of weight-admissible exponents is `k ∈ {5,6,7,8}` — a *closed*
list: no other `k` is even expressible. `m+n = 8` is in that list, so the prediction got a
fair test. Exact rational linear algebra in each of the four graded pieces kills all four,
testing **all** weight-matched minor pairs simultaneously (15 pairs at `k=8`, 67 at `k=7`,
36 at `k=6`, 1 at `k=5`):

| `k` | pairs tested | `dim`-source `|I_W|` | verdict |
|---|---|---|---|
| 5 (`= t+1`, the direct `W*Z` analogue) | 1 | 5 | **no** |
| 6 | 36 | 28 | **no** |
| 7 | 67 | 118 | **no** |
| **8 (`= m+n`)** | **15** | **383** | **no** |

So **`X*Y = e^8` is false**, and so is every other product-of-two-minors relation, at every
exponent that weight even permits. The negative survives widening the Hankel to admit the
state coefficients `d_0` (then `d_0,d_1`) as entries — 253 and 501 candidate pairs, still
nothing — and survives relaxing "product of two minors" to "any element of the toric ideal
with coefficients in the window variables" (`k <= 7`; that relaxation subsumes 3x3
catalecticant determinants, which lie in the same ideal by Laplace expansion).

### Q3. What is `lam` at `(75,125)`? — **`lam = 0`. The `(H-cap)` gate `lam >= m` FAILS.**

`lam := deg_slope - W_step = (deg_y(Phi) - ord_y(Phi))/M`.

* `(72,108)`: `(238 - 204)/17 = 2 = m`. Exactly at equality, zero margin — as flagged.
* `(75,125)` (repaired): `Phi = (1/3)*y^80`, `M = 29`, so `(80 - 80)/29 = **0**`, and
  `0 < 3 = m`.

The reason is one line and it is a **family** fact, not an `a = 3` accident: the repaired
`C` is the monomial `y`, so `deg_C - ord_C = 0`, so
`deg_y(Phi) - ord_y(Phi) = N*(deg_C - ord_C) = 0` for every rung `a = 2..8`.

**Consequence:** `CONTACT_LEMMA.md` §4.3's open flag closes **NEGATIVE**, and §4.2 (which
asserts all three gates hold at `(3,5)`) is **void as written**. The slice-cascade route at
`(75,125)` is closed. This was the cheap gating number and it came out against us.

---

## 1. The general derivation

### 1.1 What the G-system is, parametrically

`g_system_75_125.build_gsystem(a,b,t,...)` is a convolution recipe. Its closed form — which
is what makes a *general* statement possible — is this. Put `u = 1/x` and work with Laurent
series at `x = infinity`.

> **Presentation.** Let `p` be **monic of degree `a*t`** with `[x^(a t - 1)] p = 0`. Set
> `s := p^(1/a)` (the branch with `s = x^t + ...`). Then the G-system's variables are
> exactly the Laurent coefficients of `s`,
> ```
> d_i = [x^i] s   (i = t-2 .. 0),        dm_k = [x^-k] s   (k = 1 .. D),   D = (a-1)*t
> ```
> and the generators are exactly the deep Laurent coefficients of `p^(b/a)`:
> ```
> G_j  =  [x^-j] p^(b/a)          j = 1 .. j_Phi,   j_Phi = D + 1
> ```
> with `j = j_Phi - 1` skipped and `+Phi` added at `j = j_Phi`.

Machine-verified with residual exactly `0`: symbolically at `(2,3,4)` and `(2,3,5)`
(all generators, ratio exactly `1`, no fudge constant), and at `(3,5,4)` — where the
12-symbol expansion is minutes-scale — at four independent exact rational `p`-vectors, all
eight generators. Mutation: `p^((b±1)/a)` and `p^((b+a)/a)` all leave a residual.

*Why the linear window disappears.* The linear window imposes `[x^-k] s^a = 0` for
`k = 1..K_lin` except the one skipped slice; taking `s := p^(1/a)` exactly imposes all of
them including the skipped one. That is a specialisation, but a free one: the variable it
pins, `dm_{(a-1)t + K_lin - 1}` (`dm12` at `(72,108)`, `dm24` at `(75,125)`), does not occur
in any generator. So the generators are unchanged.

Three structural consequences, all checked (`A4`):

* the spare inventory is `dm1..dm_{(a-1)t}` and there are exactly `D - 1` **`Phi`-free**
  generators — a count depending on `(a,t)` only;
* `n = b` enters **only** through the exponent `b/a` and through the weight shift
  `w(G_j) = b*t + j`;
* the `u`-grading is uniform: `w(s_i) = t - i` for *every* coefficient, state and window
  alike (so `w(d_2,d_1,d_0,e,R,S,T) = (2,3,4,5,6,7,8)` at `t=4`, the familiar list).

### 1.2 The invariant form of the question

The minors `W`, `Z` are 2x2 minors of the Hankel/catalecticant matrix on `(dm1,...,dmD)`.
Let `J` be the ideal they generate — the ideal of the **cone over the rational normal
curve**, whose variety is `dm_k = lam*c^(k-1)`, i.e. *the tail of `s` is a single simple
pole* `lam/(x-c)`. Let `I = (G_1,...,G_{D-1})` be the `Phi`-free generators.

> A toric syzygy `sum_j c_j*G_j = theta - const*e^k` with `theta` in the toric ideal exists
> **iff** `e^k ∈ I + J`, **iff** `e` vanishes on `V(I) ∩ V(J)`, **iff** the generators force
> `lam = 0` on the "the window tail is one simple pole" locus.

This is the honest generalisation of `TORIC_SYZYGY.md`'s Corollary 3.6 (`W = 0` forces
`e = 0`). Note it subsumes 3x3 catalecticant determinants, which lie in `J` by Laplace
expansion — so the search below is broader than "products of two 2x2 minors".

Everything is `u`-weight homogeneous, so membership in a fixed weight is **pure linear
algebra over `Q`**: `I_W` is spanned by `{monomial * G_j}` with `w(monomial) = W - w(G_j)`.
That makes every negative below **complete within its graded piece**, not a heuristic
search.

### 1.3 The `(2,3)` positive control — constant and cofactors recovered, not typed

At `(2,3,4)` the minors have weights `12, 13, 14`; products span `24..28`; the only multiple
of `w(e) = 5` is `25`. So `k = 5` is forced **before any algebra**, and the unique
weight-25 pair is `W*Z`.

Because `e^5 ∉ I_25` (checked) and there is exactly one candidate product, the constant is
**uniquely determined**. Exact normal-form comparison modulo `I_25` returns

```
c = 6
```

and solving (not typing) for cofactors restricted to monomials in `(e,R)` returns

```
2*e^2*G3  -  4*e*R*G2  +  2*R^2*G1   =   6*W*Z  -  e^5
```

— i.e. `TORIC_SYZYGY.md` R2, reconstructed from the parametric builder with no constant
hardcoded. Mutations: dropping `G3` kills the certificate; `6*W^2` and `6*Z^2` in place of
`6*W*Z` both fail.

### 1.4 Is the constant a function of `(m,n)`?

`6 = m*n*(n-m)` at `(2,3)`, and `6 = 2*c_2` with the contact lemma's
`c_2 = m*n*(n-m)/2 = 3`. Both fit.

**Both are unfalsifiable.** In the entire swept family, `(2,3,4)` is the **only** case that
carries such an identity at all (§2). One data point cannot discriminate between two
formulas, or between them and the constant function `6`. Reported as an observation, **not**
as a law. `toric_general.py` D2 asserts precisely this — that the number of available data
points is one.

---

## 2. The sweeps

### 2.1 Fixed `(m,n) = (2,3)`, varying `t` — the identity is a `t = 4` fact

| `t` | `D` | weight-admissible `k` | candidate pairs tested | identity? |
|---|---|---|---|---|
| 3 | 3 | `{5}` (only `W^2`; `W*Z` needs `dm4`) | 1 | **none** |
| **4** | 4 | `{5}` (uniquely `W*Z`) | 1 | **HIT, `k = 5`, `c = 6`** |
| 5 | 5 | `{5,6}` | 5 | **none** |
| 6 | 6 | `{5,6}` | 10 | **none** |

Same `m`, same `n`. The identity is present only at `t = 4`, consistent with the weight law
of §0/Q1: at `t != 4` the specific product `W*Z` is **weight-forbidden** (its weight
`4t+9` is not a multiple of `w(e) = t+1`), and no *other* minor pair rescues it either —
the pairs that *are* weight-admissible were all tested and all fail. Every rung had
candidates to test (mutation control `C6`), so the empty rows are not vacuous.

### 2.2 Fixed `t = 4`, varying `(m,n)` — `(2,3)` is alone

| `(m,n)` | `D` | admissible `k` | identity? | in the checker? |
|---|---|---|---|---|
| **(2,3)** | 4 | `{5}` | **HIT at `k=5`, `c=6`** | yes (`D1`) |
| (2,5) | 4 | `{5}` | none | yes (`D1`) |
| (3,4) | 8 | `{5,6,7,8}` | none | yes (`D1`) |
| **(3,5)** — *this is `(75,125)`* | 8 | `{5,6,7,8}` | **none** | yes (`C2`, `D1`) |
| (4,5) | 12 | `{5,...,12}` | none for `k <= 9` (**partial**) | run separately, `k<=9` only |
| (2,7) | 4 | `{5}` | none, and *provably*: `I_25 = 0` | covered by P3 below |

There is also a **proved** necessary condition that does involve `n`: the identity's weight
`4t+9` must be at least the lowest generator weight `n*t + 1`, so

```
n  <=  4 + 8/t          (at t = 4:  n <= 6)
```

otherwise `I_{4t+9} = 0` and no relation of that weight can exist for any reason. That is
what kills `(2,7)` trivially. It does **not** kill `(2,5)`, `(3,4)`, `(3,5)`, `(4,5)` — those
are genuine algebraic failures, not weight bookkeeping.

### 2.3 The looser invariant: minimal `k` with `e^k ∈ I + J`

Relaxing the shape all the way (toric cofactors allowed to involve the state variables
`d_i`) gives a coarser invariant, and it *still* refuses to be `m+n`:

| case | minimal `k`, window-only cofactors | minimal `k`, any cofactors | `m+n` |
|---|---|---|---|
| `(2,3,4)` = `(72,108)` | **5** | 4 | 5 |
| `(2,3,5)` | none `<= 6` | 4 | 5 |
| `(2,3,6)` | none `<= 6` | 4 | 5 |
| **`(3,5,4)`** = `(75,125)` | **none `<= 7`** | **6** (and `k=3,4,5` fail, so minimal) | **8** |
| `(2,5,4)` | none `<= 8` | none `<= 8` | 7 |
| `(2,3,3)` | none | none | 5 |

(The `(2,3,5)`, `(2,3,6)`, `(2,5,4)`, `(2,3,3)` rows were produced by the same routine
outside the checker, which carries the `(2,3,4)` and `(3,5,4)` rows — `C5`.)

The `(72,108)` column reads `5`; the `(75,125)` column reads `6`, not `8`. Neither column is
`m+n`, and even the *existence* of some toric relation fails at `(2,3,3)` and `(2,5,4)` —
so it is not a family-level phenomenon in any form.

**What `(75,125)` does satisfy:** `e^6 ∈ I + J`, with toric cofactors that must involve
`d0,d1,d2`. It is a real relation, minimal at exponent 6. It is *not* of the `(72,108)`
shape — the thing that made `6*W*Z = e^5` usable was that its cofactor `W` is *itself a
minor*, which is exactly the property that fails here.

---

## 3. What the `(2,3)` identity actually depends on

This is the finding that justifies treating the next case as an independent endgame.

1. **`t = 4`** — through `(t+1) | (4t+9)`. Necessary, and it has no other solution.
   *Not* a fact about `(m,n)`.
2. **`(m,n) = (2,3)` specifically** — at `t = 4` the weight slot `k = 5` exists for every
   `(m,n)` with `n <= 6` and `D >= 4`, and *only* `(2,3)` satisfies it. So there is a second,
   purely algebraic condition on top of the weight arithmetic, and `(3,5)` fails it.
   The sharpest form of what fails: at `(3,5)` a toric relation **does** exist
   (`e^6 ∈ I + J`), but its cofactors **must** involve the state variables `d_0,d_1,d_2`.
   At `(2,3,4)` the cofactor of `Z` is `W` — *itself a minor*. That "the cofactor is a
   minor" property is the entire reason `6*W*Z = e^5` is usable (it is what turns the
   relation into a statement about a *product*, hence into the `v(W) + v(Z) = 5 v(e)`
   divisor law and the `Pi^4` contact order), and it is exactly the property `(3,5)` lacks.
3. **Nothing else.** The parts of `TORIC_SYZYGY.md` that said the identity is `Phi`-free,
   cap-free, slice-free and branch-free all stand — the identity is exactly as cheap as
   claimed at `(72,108)`. It is simply not *transferable*.

The right comparison is the `Phi`-divisor relation, which `weight_lemma_75_125.py` (45/45)
already showed is not family-level. **The toric syzygy joins it.** Both are `(2,3)`-at-`t=4`
resonances. The structural things that *do* transfer are the ones `g_system_75_125.py`
already established: the slice equations, the generator/spare counts, the `u`-grading. The
*algebra on top of them* does not.

Combined with **Q3** (`lam = 0`, `(H-cap)` fails), two of the three bridges from `(72,108)`
to `(75,125)` are now explicitly cut. `(75,125)` should be planned as an independent
endgame.

---

## 4. PROVED / CHECKED / INFERRED

### PROVED

* **P1.** Weight law. In the `u`-grading, `const*W*Z = e^k` requires `(t+1) | (4t+9)`, whose
  unique solution for `t >= 2` is `t = 4`, forcing `k = 5`. Elementary; `B4`.
  (Rests on the weight-homogeneity of the G-system under `w(s_i) = t-i`, which `A4`
  **verifies monomial by monomial** across 21 generators in four cases, with a mutation
  showing that grading is forced rather than one choice among many.)
* **P2.** The admissible-exponent window is finite and computable from `(t, D)`:
  minor weights lie in `[2t+4, 2t+2D-2]`, so `k` ranges over the integers in
  `[(4t+8)/(t+1), (4t+4D-4)/(t+1)]` and **no other `k` is expressible**. At `(3,5,4)` this is
  exactly `{5,6,7,8}`. This is what makes §0/Q2 a complete negative rather than a bounded
  search. `C1`.
* **P3.** Necessary condition `n*t + 1 <= 4t + 9` for the `W*Z` shape, i.e. `n <= 4 + 8/t`.
* **P4.** `deg_y(Phi) - ord_y(Phi) = N*(deg_C - ord_C)`, so `C` a monomial forces
  `lam = 0`. `E5` checks it for rungs `a = 2..8`.

### CHECKED (exact rational arithmetic; complete within the stated graded pieces)

* **C1.** The parametric presentation `G_j = [x^-j] p^(b/a)`: symbolic at `(2,3,4)`,
  `(2,3,5)`; at four exact rational points at `(3,5,4)`. `A2`, `A3`.
  *(The `(3,5,4)` leg is point-verified, not symbolic — flagged in the checker.)*
* **C2.** `(2,3,4)`: `k = 5`, `c = 6`, cofactors `(2e^2, -4eR, 2R^2)`, all recovered by
  solving. `B2`, `B3`.
* **C3.** `(3,5,4)`: **no** product-of-two-minors identity at any of `k = 5,6,7,8`; survives
  widening the Hankel by 1 and 2 state coefficients; no window-coefficient toric relation for
  `k <= 7`. `C2`, `C4`, `C5`.
* **C4.** The two sweeps of §2.1 / §2.2.
* **C5.** Minimal `e^k ∈ I + J` exponents of §2.3.
* **C6.** `lam = 2` at `(72,108)`, `lam = 0` at `(75,125)`; `(H-cap)` fails. `E1`–`E5`.

### INFERRED

*(I1 and I2 are load-bearing nowhere. **I3 IS load-bearing: it is the premise under Q3.**
Read Q3 as "`lam = 0` under the same premise that `(72,108)` proves and `(75,125)` inherits";
the alternative reading of that premise closes the route even harder — see I3.)*

* **I1.** That `(2,3,4)` is the unique member of the whole two-parameter family. Only finitely
  many `(m,n,t)` were swept (`t = 3..6` at `(2,3)`; `n <= 7`, `m <= 4` at `t = 4`), and the
  `(4,5)` row is only partial (`k <= 9` of `k <= 12`). A larger sweep could in principle find
  another; nothing in the checked data suggests one. **This is the only claim in the document
  that a wider search could overturn — and it would not touch Q1, Q2 or Q3**, each of which is
  settled *within* `(2,3,t)` / `(3,5,4)` by complete graded arguments.
* **I2.** Any formula for the constant `6`. One data point; unfalsifiable. See §1.4.
* **I3.** `lam = 0` at `(75,125)` inherits one premise from `(72,108)` — that `Phi` realises
  the extreme (minimal ord/weight) ray of the window cone. At `(72,108)` that is proved
  (`window_caps_verify.py` W5); at `(75,125)` it is inherited. The missing input that would
  make it PROVED is the **unreduced Newton polygon of `P` at the corner `A0 = (5,20)`** — the
  `(75,125)` analogue of premise `[P1]` — which no paper carries. Note the premise is
  *self-inconsistent* there (`WINDOW_FUNCTIONS_75_125.md` R3 shows the caps pinch unless
  `29 | w`), which argues the `(72,108)` window architecture has no `(75,125)` counterpart at
  all, rather than one with `lam = 0`. Either reading closes the route.

### Doc corrections this run forces

* `CONTACT_LEMMA.md` §4.2 ("gates `gcd(m,n)=1`, `lam>=m`, `N_Q>=D_P+D_Q` — all satisfied at
  `(3,5)`") is **false post-repair** and is self-contradicted by its own §4.3 flag.
  §4.3 should be closed **NEGATIVE**.
* `contact_lemma.py:1115` hardcodes `(m,n,ell,lam) = (3,5,5,3)`; `ell = 5` and `lam = 3` are
  both pre-repair. With repaired `ell = 4, lam = 0` its `capok` side-flag is FALSE.
* `CONTACT_LEMMA.md:549` `[C75]`, `:469`, `CAPS_AUDIT.md:210-215`, `WEIGHT_LEMMA_75_125.md:93`
  and `G_SYSTEM_75_125.md:167-174` all carry pre-repair `(75,125)` numbers.
  `PASSPORT_75_125_REPAIR.md:208` and `WINDOW_FUNCTIONS_75_125.md` R3 are the authoritative,
  post-repair statements (`lam: 101/12 -> 0`).

---

## 5. Reproduce

```
python toric_general.py            # verbose
python toric_general.py --quiet    # exit 0 iff all pass
python toric_general.py --fast     # skips group D and the t=6 rung
```
