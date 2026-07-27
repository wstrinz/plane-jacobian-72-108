# Is there machinery that works BECAUSE `C` is a monomial?

**Answer: partly yes, and the framing that motivated the question is wrong.
Monomiality is not ONE systematic bias — it is TWO logically independent chart
defects, and only one of them is what killed the four documented mechanisms.
The other one, which nobody had isolated, is the sole reason the `(72,108)`
Φ-divisor mechanism cannot transfer; it makes that obstruction TOTAL rather than
partial, and it MAXIMALLY STRENGTHENS the one piece of machinery that survived
(the depth ledger's forced floor).**

**Checker:** `monomial_window_law.py` — **56/56**, `--quiet` exit 0, ~4 min
(the two `build_gsystem(2,3,4,·,·)` builds dominate).
**Prerequisites:** `q_window_theorem.py` (the `M/gcd(M,H)` invariant and its Bezout
lemma), `window_functions_75_125.py` (`(R1)`–`(R3)`, `(S1)`–`(S3)`),
`weight_lemma_75_125.py` §(4) (the ord-side carry obstruction),
`MINIMAL_CORE.md` §4 (the `q_window | w(e)` criterion), `ENDPOINT_CONTRACT.md`
(the depth ledger), `corner_atlas.json` (the 34 published rows).
**Sources of truth:** the `(M,H)` and `q_window` arithmetic from
`q_window_theorem`; `ord_y(Φ)`/`M`/`N`/`ρ` from `window_functions_75_125.family`;
the 34 rows' `(t, κ, ord C, m, n)` re-read out of `corner_atlas.json` and
**re-derived**, never trusted; the G-systems from
`g_system_75_125.build_gsystem`; `ord_y(Φ) = 204` at `(72,108)` from
`AT_LE9_AUDIT.md` B7 / `AUDIT.md` §A5 (the `f1` ODE solution times `C4^28`).

Notation: `t` = chart exponent, `κ = t-2` (the standard class, and what
`polygon_reduction` returns for **all 34** rows), `(a,b) = sorted(m,n)`,
`q := ord_y C`, `M := t(a+b) - (κ+1)`, `H := q(a+b) - 1`,
`q_window := M/gcd(M,H)`, `α := ord_y(Φ)`.

---

## 0. Headline, in five lines

| # | statement | status |
|---|---|---|
| **A** | `ord_y(Φ) = a·q·M − H` — an **exact identity**, so `denom(ord_y Φ / M) = q_window`. The analytic window denominator IS the combinatorial corner invariant. | **PROVED** from the repo's `ρ`, `N` formulas (whose generality is **INFERRED**, §6) |
| **B** | Monomial `C` (`q = 1`) with `κ = t−2` forces the Bezout corner integer to `−1`, hence `gcd(M,H) = 1` and **`q_window = M`, maximal, at every corner and every family member**. | **PROVED** |
| **C** | `q_window = 1` requires `ord_y C ≥ t` (given `a+b > t−2`, true on all 34 rows). A monomial has `ord_y C = 1`. So the `(72,108)` regime is **arithmetically impossible** at a monomial corner, not merely unobserved. | **PROVED** |
| **D** | Therefore the ord-side **carry obstruction is TOTAL** at every monomial corner: *every* split of `M`, every `k`, carry `∈ [1, k−1]`, never `0` — decided with **no split enumeration and no `w(e)` datum**. | **PROVED** (carry arithmetic) + **EXACT-CHECKED** on 9337 splits |
| **E** | Monomiality = **thinness** (`deg C − ord C = 0`) **⊕ shallowness** (`ord C = 1`), logically independent. The four documented deaths all consume thinness; the `q_window` death consumes shallowness and has **zero `deg C` dependence**. | **PROVED** + all four quadrants realised |

And the **positive**: `ceil(αw/q_window)` is the depth ledger's forced floor, and
its gain over the affine ray is `((−αw) mod q_window)/q_window`, identically `0`
**iff** `q_window = 1`. So at a monomial corner the floor is strictly higher than
the affine prediction at **every** admissible weight, by exactly `(q_window−1)/2`
units over one period — and `q_window = M` is maximal there. The
`ENDPOINT_CONTRACT` kill predicate fires when a required-nonzero coefficient sits
strictly **below** the forced floor, so raising the floor can only create kills.
**Monomiality collapses the two-slope cone (which needs two slopes) and at the
same time maximally strengthens the one-slope floor.** That is the precise
mechanism behind the lane's hypothesis "the depth ledger is monomial-compatible
and the cone is not."

---

## 1. The bridge identity (A) — and the flagged negative it dissolves

`window_functions_75_125.family` records two formulas:

```
ρ := ord_y(f) = q(b-a) + 1          N := a*M - 2b          ord_y(Φ) = ρ + N*q
```

Substituting is a one-liner:

```
ord_y(Φ) + H  =  [q(b-a)+1] + [aM-2b]q + [q(a+b)-1]
              =  q[(b-a) + aM - 2b + (a+b)]
              =  a·q·M .
```

> **BRIDGE IDENTITY.  `ord_y(Φ) = a·q·M − H`, exactly** (not a congruence), in
> `(t, κ, q, a, b)`.  [checker A1, symbolic; A2, 24 948 numeric points]

Two consequences.

**(1) `gcd(M, ord_y Φ) = gcd(M, H)` identically**, so

```
denom( ord_y(Φ) / M )   =   M / gcd(M, ord_y Φ)   =   M / gcd(M, H)   =   q_window .
```

The **analytic** window denominator of `window_functions_75_125` — the object the
carry obstruction actually consumes — **is** the **combinatorial** corner
invariant of `q_window_theorem`. The repo previously had these as two objects
that happened to agree on four tabulated `KNOWN_CASES`. They are now the same
object. [A5, swept; E5, on all 34 published rows]

**(2) `MINIMAL_CORE.md`'s "Negative I could not overcome" is dissolved.** That
note says a `q_window` sweep across the 34 candidates is impossible because
`upstream_facts.json` carries Newton polygons for the `(8,28)` corner only, so
`ord_y(Φ)` is unavailable elsewhere. The identity supplies `ord_y(Φ)` from corner
data alone. In particular:

> **`ord_y(Φ) = 204` at `(72,108)` is predicted from `(t,κ,a,b,q) = (4,2,2,3,7)`
> with no Newton polygon:** `2·7·17 − 34 = 204`. The disjoint route
> (`f1 = −y^8(y+1)^2 q(y)/6630` solved from the ODE, times `C4^28 = y^196 t^28`)
> gives `8 + 196 = 204`.  [A3]

That is a cross-check between two derivations that share no inputs, and it is
tight: `ρ → ρ+1` or `N → N+1` breaks it at all four probes [MUT A]. It also
pins `ρ = ord_y(f1) = 8` and `N = 28`, matching `g4_row.py:251`'s `C4^28` [A4].
The last column of §5's table is `ord_y(Φ)` for all 34 rows — the first time the
repo has had it off the `(8,28)` corner.

---

## 2. Monomial rigidity (B) — one integer explains it

`q_window_theorem`'s Bezout lemma: `t·H − q·M = q(κ+1) − t`, so
`gcd(M,H)` divides the **fixed corner integer** `q(κ+1) − t`. That integer is the
**only** place `ord_y C` enters the window arithmetic. In the standard class
`κ = t−2`:

```
q = 1   ==>   q(κ+1) - t  =  (t-1) - t  =  -1      (independent of t, a, b)
```

so `gcd(M,H) | 1`, i.e. **`gcd(M,H) = 1` and `q_window = M` exactly.** Maximal.
[B1, B2, B3 over `t = 2..20`, `a+b = 3..60`]

Two mutation controls show this is a statement about `q = 1` and `κ = t−2`, not
about the code: with `q = 2` the corner integer becomes `t−2` and 71 swept points
have `gcd(M,H) > 1`; with `κ = t−1` the corner integer is `0`, the divisibility
lemma is vacuous, and 528 swept points have `q_window < M` [MUT B, B6]. `κ = t−2`
holds on **all 34** published rows [E2].

Combined with (A): `gcd(α, M) = gcd(M,H) = 1` at every monomial corner [B4] —
which is exactly the hypothesis (D) needs.

---

## 3. The integral regime is inaccessible (C)

`M − H = (t−q)(a+b) − κ` identically [C1], and `q_window = 1 ⟺ M | H`.

* **Monomial:** `H = a+b−1`, `M = t(a+b)−(t−1)`, so
  `M − H = (t−1)(a+b) − (t−2) ≥ t > 0`. Hence `0 < H < M` and `M | H` is
  impossible. **`q_window = 1` cannot occur at a monomial corner, for any `t ≥ 2`,
  any `(a,b)`.** [C2]
* **General necessary condition:** `M | H ⟹ M ≤ H ⟹ (t−q)(a+b) ≤ κ = t−2`. If
  `a+b > t−2` this forces `t − q ≤ 0`, i.e. **`ord_y C ≥ t`**. [C3, 1138 integral
  sweep points, 0 violations] The hypothesis is load-bearing, not decoration: 53
  integral sweep points *do* have `ord_y C < t`, and **every one** of them has
  `a+b ≤ t−2` (there `M = H` exactly) [MUT C]. On the published population
  `min(a+b) = 5` and `max t = 6`, so `t−2 ≤ 4 < 5` and the condition applies to
  every row [E5b]; the three `q_window = 1` rows have `(q,t) = (8,3), (7,4), (8,3)`
  [E6].

**Quantitatively — and this is the sharpest way to see the class of nine.** At the
corner shape `(t, κ, a+b) = (4, 2, 5)`, shared by `(72,108)` and **eight of the
nine** class rows, `M = 17` and

```
q_window = 1   ⟺   5q ≡ 1  (mod 17)   ⟺   ord_y C ≡ 7  (mod 17) ,
```

whose minimal solution is **`7` — exactly `(72,108)`'s `ord_y C`.** [C5] At the
`(75,125)` shape `(4,2,8)`, `M = 29` and the condition is `ord_y C ≡ 11 (mod 29)`
[C6]. So `(72,108)` is the **minimal integral point of its own corner shape**, and
the class of nine differs from it in exactly one integer: `ord_y C = 1` instead of
`7`. Not "a different regime" — one integer, off by six.

---

## 4. The carry obstruction is TOTAL (D), and the discriminating pair

> **TOTAL CARRY LEMMA.** Let `gcd(α, M) = 1` and `M = w_1 + … + w_k` with all
> `w_i ≥ 1`, `k ≥ 2`. Write `r_i := (−α w_i) mod M`. Each `r_i ≠ 0` (because
> `M ∤ w_i` and `gcd(α,M)=1`), `Σ r_i ≡ 0 (mod M)`, and `r_i ≤ M−1`, so
> ```
> Σ_i ceil(α w_i / M)  −  α   =   (Σ_i r_i)/M   ∈   [1, k-1] .
> ```
> **Never zero.** By (B), `gcd(α,M) = 1` at every monomial corner, so the
> obstruction is total there — for every split, every `k`, with **no split
> enumeration and no `w(e)` datum**.

Exact-checked on 9337 splits (`k = 2,3,4`) at six monomial corners including
`M = 17` `(50,75)`, `M = 29` `(75,125)`, `M = 13, 21, 25`: **zero** have carry `0`
[D1]. The contrast: at `(72,108)`, all 696 splits have carry **exactly `0`** [D2].

`MINIMAL_CORE.md` §4's criterion is "the obstruction vanishes iff `q_window | w(e)`".
Under either reading of the admissible `w(e)` range — `1..M−1` (the split
enumeration) or `{2,…,t+1}` (`window_functions.state_uweights`) — **no class-of-nine
row has an escape** [D4], because `q_window = M`.

### The discriminating pair: `(50,75)` vs `(72,108)`

`(50,75) = F_2(2,3)` has `(m,n) = (2,3)`; `(72,108) = (8,28)/(3,2)` has
`(m,n) = (3,2)`. Both give `(a,b,t) = (2,3,4)`, `κ = 2`, `M = 17`. And
`g_system_75_125.build_gsystem` forms its generators from `(a,b,t)` alone.
Verified, not assumed:

> `build_gsystem(2,3,4,1,30)` and `build_gsystem(2,3,4,7,204)` are **identical on
> every structural field** — `Gs` (the generators), `Klin`, `M`, the state/spare
> inventory, `sub`, `homog`, and the `u`-weight function symbol by symbol. The
> **only** fields that differ are `q` (1 vs 7), `ordPhi` (30 vs 204) and `W_step`
> (`30/17` vs `12`) — i.e. exactly the weight normalisation.  [D3b′, D3b″; MUT D
> confirms the comparison is sensitive: changing `t` to 5 *does* change `Gs`]

So **the K-syzygy `c·Φ = e·B` exists as an algebraic relation at `(50,75)` too.**
On the published split `(w_e, w_B) = (5,12)`:

```
(72,108)  α=204 :  ceil(204·5/17) + ceil(204·12/17)  =  60 + 144  =  204  =  α ,  carry 0
(50,75)   α= 30 :  ceil( 30·5/17) + ceil( 30·12/17)  =   9 +  22  =   31  >   30 ,  carry 1
```

[D3b] **At `(50,75)` the algebra permits the syzygy and the arithmetic forbids it,
and the single differing input is `ord_y C = 7` vs `1`.** This is the cleanest
isolation of the mechanism's failure mode the repo has: same ideal, same `M`, same
split, one integer moved.

Mutation control that this is a monomial fact and not a property of `carry()`: at
`F_7` (`t=3, κ=1, q=4, (a,b)=(2,7)`, `M=25`, `α=165`, `q_window=5`) the 2-splits
with carry `0` are **exactly** `w ∈ {5,10,15,20}` — the multiples of `q_window`
[MUT D].

### What this resolves in the atlas

Atlas gate `G5` records **31 UNKNOWN**, each citing the missing `w(e)` split
enumeration. The total-carry lemma decides them without it:

* `w(e) ∈ 1..M−1`: an escape exists at **exactly the 6 non-monomial rows**; the
  obstruction is TOTAL at **28/34**. [E7a]
* `w(e) ∈ {2,…,t+1}`: an escape exists at **exactly 5 rows** — `F_7`
  (`q_window = 5 > t+1 = 4`) also loses its escape. [E7b]

Every row with an escape is non-monomial, in both readings.

---

## 5. All 34 rows, with `ord_y(Φ)` for the first time

`mono` = retraction shape fails ⇒ `C = y`. `q = ord_y C`. `escape` = the `w(e)`
with `q_window | w(e)`, `0 < w(e) < M`. Every number re-derived from the atlas's
own transcribed `(t, κ, ord_C, m, n)`, not read from the stored `G5` block [E1].

```
 id                       mono   t  q   M   H  q_win  ord_y(Phi)  escape w(e)
 F_1(3,4)/64              y      3  1  19   6    19         51   NONE
 F_1(5,7)/112             y      3  1  34  11    34        159   NONE
 F_2(2,3)/75              y      4  1  17   4    17         30   NONE      <- (50,75)
 F_2(3,5)/125             y      4  1  29   7    29         80   NONE      <- (75,125) flagship
 F_3(3,2)/75              y      4  1  17   4    17         30   NONE
 F_7(2,7)/147             -      3  4  25  35     5        165   [5,10,15,20]
 F_8(3,7)/147             -      3  5  28  49     4        371   [4,8,12,16,20]
 F_9(2,3)/84              y      3  1  13   4    13         22   NONE
 F_9(3,5)/140             y      3  1  22   7    22         59   NONE
 F_11(2,5)/140            y      3  1  19   6    19         32   NONE
 F_17(2,3)/99             -      3  8  13  39     1        169   ALL
 F_22(2,3)/96             y      3  1  13   4    13         22   NONE
 F_24(3,4)/128            y      3  1  19   6    19         51   NONE
 (7,35)/(2,3)/126         y      5  1  21   4    21         38   NONE
 (7,42)/(3,2)/147         y      6  1  25   4    25         46   NONE
 (7,42)/(2,3)/147         y      6  1  25   4    25         46   NONE
 (8,28)/(3,4)/144         -      4  3  25  20     5        205   [5,10,15,20]
 (8,28)/(3,2)/108         -      4  7  17  34     1        204   ALL       <- CLOSED
 (9,36)/(3,2)/135         y      4  1  17   4    17         30   NONE
 (9,36)/(2,3)/135         y      4  1  17   4    17         30   NONE
 (11,33)/(2,3)/132        y      3  1  13   4    13         22   NONE
 (12,33)/(2,3)/135        -      3  8  13  39     1        169   ALL
 (8,32)/(3,2)/120         y      4  1  17   4    17         30   NONE
 (8,40)/(3,2)/144         y      5  1  21   4    21         38   NONE
 (9,27)/(2,3)/108         y      3  1  13   4    13         22   NONE
 (9,36)/(2,3)/135         y      4  1  17   4    17         30   NONE
 (10,40)/(3,2)/150        y      4  1  17   4    17         30   NONE
 (10,40)/(3,2)/150        y      4  1  17   4    17         30   NONE
 (12,30)/(3,2)/126        y      3  1  13   4    13         22   NONE
 (12,36)/(2,3)/144        y      3  1  13   4    13         22   NONE   (x5 rows)
 (12,36)/(3,2)/144        y      3  1  13   4    13         22   NONE
```

> **`C` monomial ⟺ `q_window = M` on 34/34 rows** [E4]. **Only the forward
> direction is PROVED.** The converse is a **population fact, not a theorem**:
> 2443 abstract points with `q ≥ 2` also have `q_window = M` [MUT E]. Reported as
> **EXACT-CHECKED on the published population**, never as proved.

The **class of nine** (`t = 4`, retraction fails): eight rows at `M = 17`, and
`(75,125)` at `M = 29`. All nine have `q_window = M` and `gcd(α,M) = 1` [B5].

---

## 6. Status, honestly

| item | status |
|---|---|
| Bridge identity `ord_y(Φ) = a·q·M − H` as an algebraic consequence of `ρ = q(b−a)+1`, `N = a·M − 2b`, `ord_y Φ = ρ + N·q` | **PROVED** (symbolic) |
| **that those three formulas hold at a general corner** | **INFERRED.** `window_functions_75_125` states them inside the F2-family block (`ORD_C = 1`). They are confirmed at `(72,108)` (`q = 7`, `b−a = 1`, reproducing the published `204`) and across F2 rungs `a = 2..8` (`q = 1`, `b−a = a−1` varying) — i.e. on a 1-dimensional slice in each variable **separately**. The joint `(q, b−a)` dependence is untested. **This is the load-bearing soft spot of the whole writeup.** |
| `q = 1, κ = t−2 ⟹ gcd(M,H) = 1 ⟹ q_window = M` | **PROVED** |
| `q_window = 1` impossible at a monomial corner; `q_window = 1 ∧ a+b > t−2 ⟹ q ≥ t` | **PROVED** |
| Total-carry lemma (carry `∈ [1,k−1]`, never 0, when `gcd(α,M)=1`) | **PROVED**; exact-checked on 9337 splits |
| **"carry ≥ 1 ⟹ no Φ-divisor relation"** | **PREMISE-DEPENDENT.** It inherits the extreme-ray premise of `window_functions_75_125` ("`Φ` realises the minimal-ord ray of the window cone"), unchanged in status. Independently corroborated at `(75,125)` by `weight_lemma_75_125` §B's direct graded search, which finds no relation under the corrected criterion. The arithmetic prediction and the direct search agree at the one corner where both are computable. |
| `(50,75)` and `(72,108)` have identical G-systems | **EXACT-CHECKED** (field-by-field object comparison, with an `(a,b,t)`-sensitivity control) |
| `q_window` computed for all 34 rows, `(M,H,q_window)` re-derived | **EXACT-CHECKED** |
| `C` monomial ⟺ `q_window = M` | **⟹ PROVED; ⟸ EXACT-CHECKED on the 34 rows only, with an explicit abstract counterexample** |
| Thinness ⊕ shallowness independence, four quadrants realised | **PROVED / EXACT-CHECKED** |
| Floor-gain law `Σ_{w=1}^{q_window−1} gain(w) = (q_window−1)/2` | **PROVED** (bijection of residues), checked at six corners |
| **"a raised floor produces kills at the class of nine"** | **CLAIMED, NOT CHECKED.** §7. |

---

## 7. The defect decomposition (E), and what it means for repairs

The four documented deaths at a monomial corner:

| mechanism | dies because |
|---|---|
| slice cascade | `λ = N(deg C − ord C)/M = 0` |
| window cone | `ord_y Φ = deg_y Φ`, cone → ray `(R3)` |
| toric identity | needs the retracting corner |
| F2 closed form | `dg = deg C − ord C = 0` |

**All four consume `deg C − ord C = 0` — THINNESS.** None of them reads `ord C`
except through that difference [G2].

The `q_window` death consumes something else entirely:

```
q_window = M / gcd(M, q(a+b)-1) ,    M = t(a+b)-(κ+1)
```

**No formula in that chain mentions `deg_y C`** [G1]. It reads `ord_y C = 1`
alone — **SHALLOWNESS**.

All four quadrants are realised by explicit chart data [G3]:

| `deg C` | thin? | shallow? | `q_window` | `λ` | example |
|---|---|---|---|---|---|
| `C = y` | yes | yes | `M` (17) | 0 | **the class of nine** |
| `C = y^8`, `t=3, κ=1, (a,b)=(2,3)` | yes | no | **1** | 0 | a **monomial** `C` with an **integral window** [G4] |
| `C = y(y+1)`, `t=4, κ=2, (a,b)=(2,3)` | no | yes | `M` (17) | ≠0 | residual restored, `q_window` unmoved [G5] |
| `C = y^7(y+1)` | no | no | **1** | 2 | `(72,108)` |

**Consequences for repair strategy — this is the decision-relevant part.**

1. **Monomiality per se does NOT kill the divisor mechanism.** Row 2 is a monomial
   `C` with `q_window = 1`. What kills it is `ord_y C = 1`.
2. **Restoring a residual to `C` is not a repair for `q_window`.** Row 3 has
   `λ ≠ 0` (cascade and cone alive) and still `q_window = M` (divisor mechanism
   still dead). Any attempt to revive the F2 closed form by finding a residual
   leaves the Φ-divisor lane exactly as dead as before.
3. **Deepening `C` is not a repair for `λ`.** Row 2 has `q_window = 1` and still
   `λ = 0`.
4. Therefore **the "one systematic bias" reading in the lane brief is refuted.**
   There are two orthogonal defects requiring two different repairs, and neither
   repair helps the other [MUT G].

---

## 8. The positive: monomiality maximally strengthens the depth-ledger floor

`(S1)` of `window_functions_75_125` gives the lower `y`-order cap as
`L(w) = ceil(α w / q_window)` — and `ENDPOINT_CONTRACT.md` §2's kill predicate is

```
KILL  ⇔  a required-nonzero coefficient's y-order lies strictly BELOW the forced floor.
```

`L` **is** that floor. Its gain over the affine ray is

```
gain(w) = ceil(αw/q_window) - αw/q_window = ((-α w) mod q_window)/q_window ,
```

which is `> 0` for every `w` with `q_window ∤ w` [F1] and is **identically 0 iff
`q_window = 1`** [F3]. Because `(−αw) mod q_window` is a bijection of
`{1,…,q_window−1}` when `gcd(α,q_window) = 1`,

```
Σ_{w=1}^{q_window-1} gain(w)  =  (q_window - 1)/2      exactly.
```

[F4 at six corners; F2 the monomial specialisation, where `q_window = M`, giving
`(M−1)/2` = 8 at `M = 17` and 14 at `M = 29`; MUT F shows quoting it with `M` at a
cancelling corner would be wrong]

At `(72,108)` — where every mechanism in this repo was calibrated — the sum is
**0**: the floor is exactly affine and there is no gain at all. **`(72,108)` is
the case with the weakest possible floor, and it is the only `t = 4` row that
retracts.** The class of nine has the maximal one.

This explains *why* the depth ledger survived the `2adb92a` chart repair untouched
while the cone did not, in structural rather than historical terms:

* the **cone** needs **two** slopes (`ord`-lower and `deg`-upper). Thinness makes
  them coincide. Cone → ray. Dead, and irreparably so at fixed `deg C = ord C`.
* the **floor** needs **one** slope, and shallowness makes that one slope
  maximally quasi-affine. The ledger is not merely monomial-compatible — it is
  **stronger** at a monomial corner.

Consistently with this: the one live depth-ledger kill in the repo,
`(50,75)`/`ENDPOINT_CONTRACT.md` §3 (`f2_tower_verify` §B, exact PASS), is
**itself at a monomial corner** — `C = y`, `q = 1`.

**What is NOT established.** That the raised floor actually *fires* the kill
predicate at any of `(8,32)`, `(9,36)`, `(10,40)`. That needs each corner's
reduced γ-chart caps and its `required-nonzero` primitivity slot — GGV3 §3/§5
data that exists in-repo for `(50,75)` only. §8 establishes that the floor is
higher and that a higher floor can only create kills; it does **not** claim a
kill. Marked **CLAIMED, NOT CHECKED** in §6 deliberately.

---

## 9. Negatives, stated plainly

* **No new closure.** Nothing here kills a case. The strengthenings are on
  obstructions to *our* mechanisms, plus one quantified improvement to the depth
  ledger's floor that is not yet cashed at any new corner.
* **The Φ-divisor lane is now provably dead class-wide, not just at `(75,125)`.**
  Under the extreme-ray premise, all 28 monomial rows admit no Φ-divisor relation
  at any split. Anyone planning to port the `(72,108)` K-syzygy to the class of
  nine should stop. **What would revive it:** a corner with `ord_y C ≥ t` — and by
  §3 the retraction shape is what supplies large `ord_y C` (`ord C = b_final`
  there, `1` off it). So the lane is coextensive with the retracting corners, of
  which GGV5 publishes six, one of them already closed.
* **The `(a,b,t)` coincidence cuts both ways.** `(50,75)` sharing `(72,108)`'s
  entire G-system is an opportunity as well as a diagnosis: any *algebraic*
  consequence of that ideal proved for `(72,108)` transfers verbatim to `(50,75)`,
  `F_3(3,2)/75`, and every other `(a,b) = (2,3)`, `t = 4` row — **eight of the nine**.
  Only weight-normalised consequences fail to transfer. **This is the strongest
  surviving lead and it was not visible before.**
* **The `ρ`/`N` generality gap (§6) is unresolved** and is the one place this
  writeup could be wrong in a way that matters. A single further corner with an
  independently derived `ord_y(Φ)` — `(8,28)/(3,4)/144`, where the bridge predicts
  `205` — would settle it. `PHI_KNOWN` in `corner_atlas.py` has one entry; a
  second would turn the bridge from INFERRED to CHECKED-in-two-directions.
* **`(8,32)` remains the outlier.** `gamma_from_corner.py` finds no surviving
  branch there at all; nothing here changes that, and its window arithmetic
  (`M = 17`, `q_window = 17`, `α = 30`) is *identical* to `(50,75)`'s, so the
  window layer cannot be what distinguishes them.

---

## Files

| file | role |
|---|---|
| `MONOMIAL_WINDOW_LAW.md` | this writeup |
| `monomial_window_law.py` | the checker — 56/56, `--quiet` exit 0 |
| `q_window_theorem.py` | the `M/gcd(M,H)` invariant and Bezout lemma this builds on |
| `window_functions_75_125.py` | `ρ`, `N`, `ord_y Φ`, `(R1)`–`(R3)`, `(S1)`–`(S3)` |
| `weight_lemma_75_125.py` | §(4), the carry obstruction; §B, the independent graded search |
| `corner_atlas.json` / `corner_atlas.py` | the 34 rows and gate `G5`, whose 31 UNKNOWNs §4 resolves |
| `ENDPOINT_CONTRACT.md` | the depth ledger and kill predicate that §8 strengthens |
| `MINIMAL_CORE.md` §4 | the `q_window \| w(e)` criterion and the flagged negative §1 dissolves |
