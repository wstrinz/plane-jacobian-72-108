# The class of nine at its own `y`-place: the three routes are DEAD, and the obstruction is one integer

**Answer to the assigned question: no route closes. But the reason is not the one
the target was built on.**

The eight `(a,b,t) = (2,3,4)` class rows are not "`(72,108)` minus the cascade".
They are literally **`(72,108)` with `C_4 = y^7(y+1)` replaced by `C = y`**, and at
the place where `C` is thin — the `y`-place — the *entire* spine §§2–7 transfers,
**cascade included**. So the class rows reach exactly the same place `(72,108)`
reaches: `ord_y(e) = 9`, `Π = 1`, `k = 0`, Corollary 8.5. And there the transfer
stops, at a single integer: Corollary 8.5's engine is `deg(μ t³ q) = 7 > 6`, and
at a class row `q` is replaced by `1`, so that degree is `3 ≤ 6`.

Three results, in order of value:

| | statement | status |
|---|---|---|
| **1** | **All three routes `{5}`, `{1,2,7}`, `{2,6,7}` are REFUTED.** An explicit family satisfies *every* slice condition and *both* polygon caps — everything the class row's `y`-place proves about valuations — while having `ord_y(h₁) = 1`, `ord_y(h₂) = 3`, `ord_y(h₅) = 9`. | **EXACT-CHECKED**, scope stated in §6 |
| **2** | **`WEIGHT_FREE_TRANSFER.md` rows #5, #6 and G1/G2 are WRONG.** The slice families and the cascade are *not* lost: they hold at the class rows' `y`-place with the same exponents, and `λ = 2`, not `0`. Recomputed here from scratch, levels 2–12. | **PROVED** (same imports as §6.1) / **EXACT-CHECKED** |
| **3** | **§8 cannot kill the class rows either.** `k = 0` forced, §8.5 vacuous, Cor 8.5 fails; four explicit witnesses satisfy the whole §8.1 system *and every transferred cap*. The calibration control shows the same machinery **does** kill `(72,108)`. | **EXACT-CHECKED**, scope stated in §7 |

**Checker:** `yplace_transfer.py` — **57/57**, `--quiet` exit 0, ~3 min
(`--fast` skips cascade levels 10, 12: 55/55 in ~8 s).
**Prerequisites:** `PROOF_72_108.md` §§2.1–2.6, 5, 6.1–6.3, 7.1–7.5, 8.1–8.7;
`WEIGHT_FREE_TRANSFER.md` (the target being audited); `polygon_reduction.py`
(the computed reduced polygons and the chart guard).
**Sources of truth, nothing retyped:** reduced Newton polygons from
`polygon_reduction.all_reductions()` (with `case_8_28` as the *published*
control); chart data from `polygon_reduction.corner_chart_data`; window
arithmetic from `window_functions_75_125`; generators from
`g_system_75_125.build_gsystem` / `.published_72108`; every `(72,108)` number
quoted at its point of use with the `PROOF_72_108.md` section.

**Scope note.** Everything below is about the **eight** `(a,b,t) = (2,3,4)` rows.
`(75,125) = (3,5,4)` has `a = 3`; the slice identity of §2 has different exponents
there and none of this applies. `WEIGHT_FREE_TRANSFER.md` §2.1 already showed it
inherits neither syzygy.

---

## 1. The chart identity: one ODE, two `C`'s

`PROOF_72_108.md` §2.1–2.2 is a function of `(a,b,t,κ)` and `C` alone. At
`(a,b,t,κ) = (2,3,4,2)` — shared by `(72,108)` **and** all four class-row corners
`(5,20)`, `(8,32)`, `(9,36)`, `(10,40)`, each of which
`polygon_reduction.corner_chart_data` returns with `deg C = ord C = 1`, i.e.
**`C = y`** [A1] — the following are *identical*, not merely analogous [A2]:

```
v_{1,0}(P) = a*t = 8      v_{1,0}(Q) = b*t = 12      3*8 = 2*12  (alignment (2,3))
[P,Q] = x^kappa = x^2     v_{1,0}(F) = kappa+1-a*t = -5  exactly, and 4 does not divide -5
M = t(a+b)-(kappa+1) = 17          N = a*M - 2b = 28
```

and `R = x^t C` is primitive in both cases [A3]. The **only** differing input is
`C`. Feed both into the one forcing ODE

```
a { t*C*f' - [t(b-a)+kappa+1]*C'*f }  =  C^(b-a+1)
```

solved with a *general* polynomial ansatz (uniqueness asserted, no shape
assumed):

| `C` | the ODE | its unique polynomial solution | `Φ = f·C^28` |
|---|---|---|---|
| `y^7(y+1)` | `8y(y+1)f' − 14(8y+7)f = y^8(y+1)^2` | `f₁ = −y^8(y+1)^2 q(y)/6630` | `ord_y = 204`, `deg = 238`, `mult_t = 30`, `[t^30]Φ = −1/2` |
| `y` | `8y f' − 14 f = y^2` | `f = y^2/2` | `Φ = (1/2)y^30`, `ord_y = deg_y = 30` |

The first row **is** `PROOF` §2.2 and Lemma 2.2 verbatim, quartic and `6630` and
all [A4, A5]; the second is the class row [A6], cross-checked against the bridge
identity `ord_y Φ = a q M − H`, which also shows `ord_y Φ = 30` is independent of
the corner `a₀` and therefore holds at all four class-row corners [A7]. And
`build_gsystem(2,3,4,1,30)` reproduces the published `(72,108)` generators term
by term [A8].

> **`Φ`'s residual factor is the single difference.** `q(y)` of degree 4 becomes
> `1` of degree 0. Remember that number.

---

## 2. `(2.5.1)` is PLACE-BLIND — this is the pivot

`PROOF` §2.5 imports the two slice families `t^{2n−2} | p_n`, `t^{2n−3} | r_n`
as `[I1]`, `[I2]`, and `WEIGHT_FREE_TRANSFER.md` row #5 reads their exponents as
"the `t`-place ord/deg data of `C_4 = y^7 t`". They are not that. With
`c_j = D_j / C^{2(t−j)−1}` and `h_k := D_{t−k}` (the definitions of §2.3/2.5),

```
      P_M  =  [u^n] H^2 / C^(2n-2) ,   n = 2t - M
(script-C^3)_M = [u^n] H^3 / C^(2n-3) ,   n = 3t - M
```

**identically, with `C` and every `D_j` a free symbol** [B1, B1b], and the
exponents are exact — shifting either by `±1` breaks it at every slice [B1c].
So polynomiality of the slices *is* `C`-divisibility, and at a place `β` with
`mult_β(C) = μ` the `β`-primary content is `β^{μ(2n−2)}`, `β^{μ(2n−3)}`.

* At `(72,108)`, `mult_t(C_4) = 1`, so the `t`-part is `t^{2n−2}`, `t^{2n−3}`:
  **that is `(2.5.1)` verbatim** [B2].
* At a class row, `mult_y(C) = 1`, so the `y`-part is `y^{2n−2}`, `y^{2n−3}`:
  **the same condition set, with the same exponents** [B3].
* The `n`-ranges match too: `min_i N(P) = 0` in both computed reduced polygons, so
  `p_n = 0` for `n ≥ a t + 1 = 9`; and the `Q`-slices run `n = 2..b t + 3 = 15`
  in both [B4].

**What does the work is `mult_β(C) = 1`, not monomiality.** Mutation: at `C = y²`
the exponents *double* to `y^{4n−4}`, `y^{4n−6}` [B3b, H1], so a deeper monomial
corner would inherit nothing.

---

## 3. Therefore the cascade transfers — recomputed from scratch

`WEIGHT_FREE_TRANSFER.md` classifies #5 (slice cokernels) and #6 (the cascade) as
**"WEIGHT-DEPENDENT, LOST"**, on the ground that `CONTACT_LEMMA`'s gate `λ ≥ m`
fails because `λ = 0`. Both premises are wrong: the gate is asked of the wrong
`λ` (§4), and the cascade needs no gate at all — it needs the two divisibility
families of §2, which the class row has.

Run at the class row's `y`-place from scratch — absorption
`h_n = −½q_n + y^{2n−2}g_n` for `n ≤ 8` (which makes every `P`-condition hold
identically), `h_n = −½q_n` for `n ≥ 9` (the `(P<)` vanishing), then impose
`y^{2n−3} | r_n` level by level:

| level `n = 2m` | obstructed coefficient | the jet | forces |
|---|---|---|---|
| 2 | `y^0` only | `(3/2)([y^0]h₁)²` | `ord_y h₁ ≥ 1` |
| 4 | `y^4` only | `(3/2)([y^2]h₂)²` | `ord_y h₂ ≥ 3` |
| 6 | `y^8` only | `(3/2)([y^4]h₃)²` | `ord_y h₃ ≥ 5` |
| 8 | `y^12` only | `(3/2)([y^6]h₄)²` | `ord_y h₄ ≥ 7` |
| **10** | `y^16` only | `(3/2)([y^8]h₅)²` | **`ord_y h₅ ≥ 9`** |
| **12** | `y^20` only | `(3/2)([y^10]h₆)²` | **`ord_y h₆ ≥ 11`** |

[C2–C7]. Every firing jet is a **perfect square** with a *constant* leading
coefficient, linear in the one fresh coefficient `g_{m,0}` — `PROOF` §6.1's own
shape, reproduced at a different place; the level-12 jet is a perfect square only
*after* the lower levels are imposed, which is §6.1's "the ordering is
load-bearing", also reproduced. Lemma 6.1 is pure `min`-arithmetic over the
absorption and transfers as such: `ord_y h₇ ≥ 12`, `ord_y h₈ ≥ 13` [C8].

> **The class rows' `y`-place carries `PROOF (6.2.1)` entire:
> `(1,3,5,7,9,11,12,13)`.** **PROVED**, consuming exactly what §6.1 consumes at
> `(72,108)` — Prop 2.1 (`[QQ1]`, which is what gives the `Q`-column at all) and
> the identification `h₅ = e` (`[I3]`) — and no more.

---

## 4. The polygon caps: `λ = 2`, not `0`

`window_functions_75_125` `(R3)` and `WEIGHT_FREE_TRANSFER` G1/G2 report `λ = 0`
and "the affine degree cap does not exist" at a class row. That `λ` is
`deg_y Φ/M − ord_y Φ/M`, i.e. **both slopes read off `Φ`** under the extreme-ray
premise — and `Φ` is a monomial, so of course they coincide. The **polygon** says
something else.

`PROOF` §2.6's cap lemma for `a = 2` is, in closed form, an affine induction that
closes identically:

```
ord_y P_M >= sigma*M - m   ==>   ord_y h_k >= k*(a*q   - sigma)
deg_y P_M <= tau  *M + c   ==>   deg_y h_k <= k*(a*degC - tau)
```

with `σ`, `τ` the terminal slopes of the lower / upper hull of `N(P)`. Fed
`polygon_reduction`'s **computed** reduced polygons:

| case | `N(P)` (computed) | `σ` | `τ` | `ord_y h_k ≥` | `deg_y h_k ≤` | `λ` |
|---|---|---|---|---|---|---|
| `(8,28)` sub1 | `{(0,0),(1,0),(0,8),(8,14),(8,16)}` | 2 | 1 | `12k` | `15k` | 3 |
| `(8,28)` sub2 | `{(0,0),(1,0),(8,14),(8,16)}` | 2 | 2 | `12k` | `14k` | 2 |
| class row `(5,20)` | `2·{(0,0),(3,0),(4,1),(0,5)}` | 1 | −1 | **`k`** | **`3k`** | **2** |

The first two rows **are** `PROOF` §2.6(i)'s three direction functionals
(`max(2i−j) = 2`, `max(j−i) = 8`, `max(j−2i) = 0`) and §2.6(iii)'s
`ord D_{j_x} ≥ 48−12j_x`, `deg ≤ 15k / 14k`, `λ = 3 / 2` — recovered, not quoted
[D1, D2]. **This is the calibration control for the whole section.** The third
row is the class row [D3]: **two distinct integral slopes, `λ = 2`, and `Φ`
strictly inside the cone** (`17 = 1·M ≤ 30 ≤ 3·M = 51`) [D3b]. The other three
class-row corners run through the same construction with the **same** ord slope 1
and a *weaker* deg slope (`15/4`, `4`, `17/4`), so `(5,20)`'s deg caps are the
tightest of the four [D5].

Two consequences.

* **Step #11 has a counterpart after all.** Unstripped at a class row the deg
  ledger is `deg(d₂,d₁,d₀,e,R,S,T) ≤ 3w = (6,9,12,15,18,21,24)` — *numerically
  the config-(1) table of §2.6* [E3]. The shift preserves it term by term because
  the per-step slope **is** `cap(w=1) = 3`.
* **The floor `L` is not what the polygon gives.** `L(w) = ceil(30w/17)` is
  `(2,4,…,9,11,13)`; the polygon ray is `w`. It is the *cascade* of §3, not the
  polygon, that supplies `(6.2.1)` — and `(6.2.1)` is still one unit under `L` at
  `w = 1,2,7`, exactly as `WEIGHT_FREE_TRANSFER` §3 measured [D4].

---

## 5. §7 transfers entire: `ord_y(e) = 9` exactly

The dictionary `(7.1.1)` is generalized-binomial algebra in `t = 4` alone,
re-derived here from `tilde D_j = Σ_m binom(m,m−j) D_m θ^{m−j}` [E1]; fed the
`y`-place profile, Lemma 7.4's four `min`'s give the paper's ledger **number for
number**, `ord_y(d₂,d₁,e,R,S,T) ≥ (2,3,9,10,11,12)` [E2]; the bracket collapse
is weight-free [E4]. Hence, with `ord_y Φ = 30` *exactly* and the K-syzygy
`2Φ = eB` (`WEIGHT_FREE_TRANSFER` B1, proved at the class row's own generators),

```
a + ord_y(B) = 30   (an EQUALITY),   the four terms of B on (21,21,21,22)
a = 9 survives;  a >= 10 refuted for every a in [10,59]
```

[E5] — so **`ord_y(e) = 9` exactly at every class row**, by Theorems 6.2 and 7.2
transferred. Then Theorem 3.4 pins `e` completely: `e | 2Φ = y^30` gives
`e = γ y^n`, and `n = 9`, so

```
e = gamma * y^9   (a monomial; Pi = 1, k = deg Pi = 0 FORCED)
B = 2Phi/e = y^21/gamma   (a monomial, ord = deg = 21)
```

[E6]. This is *more* than `(72,108)` has at its `t`-place, where `B = 2Φ/e` is
`t^{21}` times a nonconstant unit.

---

## 6. The three routes are REFUTED

Re-enumerated independently: over `(6.2.1)`, the minimal weight sets whose `+1`
upgrade makes `ord(e) + min_B > 30` are exactly `{5}`, `{1,2,7}`, `{2,6,7}`, and
`B` sees only weights `1,2,5,6,7`, so the list is complete [F1, F1b] —
`WEIGHT_FREE_TRANSFER` §4's enumeration is **confirmed**.

Now the witness. Put `h_k = c_k y^{2k−1}` (`c_k ≠ 0`) for `k = 1..8` and let
`h₉..h₁₆` be **forced** by the `(P<)` vanishing `[u^n]H^2 = 0`. Then:

* `y^{2n−2} | p_n` for `n = 2..8`, `p_n = 0` for `n = 9..16`, and
  `y^{2n−3} | r_n` for **all** `n = 2..15` — every slice condition of §2 [F2];
* both polygon caps of §4: `k ≤ ord_y h_k = 2k−1` and `deg_y h_k = 2k−1 ≤ 3k`
  [F2b].

So it is a point of **everything the class row's `y`-place proves about
valuations**. And on it

```
ord_y(h1) = 1 < 2 ,   ord_y(h2) = 3 < 4 ,   ord_y(h5) = 9 < 10 .
```

**Every one of the three routes is false on it** [F3]: `{5}` needs `ord e ≥ 10`,
`{1,2,7}` needs `ord h₁ ≥ 2`, `{2,6,7}` needs `ord h₂ ≥ 4`. The reason the
witness satisfies everything is a one-line identity that transfers verbatim from
§6.1's own sharpness argument: with `u = v/y²`,

```
[u^n]H^2 = 2 y^(2n-1) c_n + y^(2n-2) [v^n]Chat^2
[u^n]H^3 = 3 y^(2n-1) c_n + 3 y^(2n-2) [v^n]Chat^2 + y^(2n-3) [v^n]Chat^3
```

[F3b].

> **`{1,2,7}` — reported as the cheapest route, needing "exactly the floor at
> three weights" — is not merely unproved. It is FALSE on a point of its own
> hypothesis set.** So is `{2,6,7}`, and so is `{5}`.

**Scope, stated exactly.** What is refuted is: *the routes are consequences of
the two slice families together with the polygon caps* — i.e. of the entire
valuation-supplying input at a class row's `y`-place. A route could still in
principle be proved by *also* using the `G`-system equations (a proof by
contradiction may use `G₅ = 0`); the witness is not shown to be a `G`-point and
this document does **not** claim the class rows are nonempty. But `PROOF` §7.4(c)
now transfers too [F4]: at `a = 9` two terms of `B` sit on 21, `h₈` does not occur
in `B`, deeper cascade levels supply **valuations only**, and killing `a = 9`
requires a *non-vanishing* statement about `[y^{21}]B` — leading coefficients,
which is exactly what the three routes are not.

**The calibration control on this section.** The refutation does *not* fire at
`(72,108)`'s `t`-place spuriously, because it is the *same* refutation: the same
witness family is `PROOF` §6.1's own sharpness family, and §6.1/§7.3 already
record `v_t(h_k) = 2k−1` attained for `k = 1..4` and `v_t(h₂) ≥ 3` sharp. The two
places are isomorphic (§2), so the sharpness that `(72,108)` publishes *is* the
class rows' sharpness. Nothing here proves `(72,108)` dies at §3 — the floor `L`
is not proved at either place, and §4 says explicitly what the polygon does and
does not give.

---

## 7. …and §8 cannot kill them either

Because §§3–7 all transfer, the class rows arrive exactly where `(72,108)`
arrives, at the `k = 0` case:

* `rad Φ = y` alone, so `Π = 1`, `k = 0` **forced**; `k = 1,2,3,4` do not arise
  and §8.5's `Π²`-support test is **vacuous** [G1] (this part
  `WEIGHT_FREE_TRANSFER` D4 had right).
* The §8.4 window transfers unchanged: `ord_y Z ≥ 2`, `ord_y F ≥ 3`, hence
  `z = ord_y(Z) ∈ [2,6]` — the paper's own window [G2].
* Corollary 8.5's boxed row transfers: `γu = μ π³ Q_Π − 6A² + 3ζ π^z`,
  `deg u ≤ 6`.

**The calibration control, which is the most important check in this file.** At
`(72,108)`, `deg(μ t³ q) = 7` exactly. Exhaustively over `deg A ∈ {A=0} ∪ [0,39]`
and `z ∈ [2,6]`, *some* contributor attains a degree `> 6` **uniquely** in every
single case, so `deg u ≤ 6` is contradicted every time: **Corollary 8.5 fires,
and this machinery gives `(72,108)` the right answer** [G3]. It also reproduces
§8.6's two published zero-margin sensitivities — the kill switches off at `z = 7`
and at cap `7` [G3b].

**At a class row it has no engine.** `Q_Π = 1`, so `deg(μ y³ Q_Π) = 3 ≤ 6`, and
**25** of the `(deg A, z)` pairs — every `deg A ≤ 3` (and `A = 0`) against every
`z ∈ [2,6]` — pass outright [G4]. The failure is located at one integer:
`deg(π³Q_Π) = 3 + deg q`, the kill needs `> 6`, i.e. `deg q ≥ 4`; `(72,108)` has
exactly 4 (zero margin, as §0.4 says) and a class row has 0 [G4b].

And it is not merely that the *argument* fails — the *system* is satisfiable.
For `(A,z,ζ,γ) = (y,2,1,1)`, `(y,3,1,1)`, `(y²,4,1/3,1)` and `(0,5,2,1)` the full
§8.1 system at `k = 0` — `g₁ = g₂ = g₃ = □ = 0` **and** `(*) FZ = (1/6)γ⁵y⁹` —
holds with residual **exactly 0**, respecting **all six** transferred degree caps
(`deg A,u,v,w,d₀,C ≤ 9,6,12,9,12,15`) and **all five** transferred orders
(`ord A,u,v,w,C ≥ 1,2,2,3,3`) [G5]. `g₃ = 0` comes out automatically, which is a
cross-check rather than luck: Theorem 8.1 forces it once `g₂ = 0` and `(*)` hold,
and the witnesses were built from the boxed row and `(*)` only [G5b]. `Z` and `F`
are the required nonzero monomials with `deg F + deg Z = 9 = 9+4k` [G5c].

For example (`A = y`, `z = 3`, `ζ = γ = 1`):

```
A  = y            u = 2y^2(2y-3)          v = -y^2(y-1)
w  = y^3(y^3-12y+24)/6                    d0 = -y^4(y^3-3y^2+9)/3
CT = -y^3(y^3+6y-6)/6                     Z = y^3      F = y^6/6
```

> **The class of nine is not merely un-killed at §3. It is un-killable by the
> whole of `PROOF` §§4–8 as it stands** [G6] — and the discriminant is `deg q`,
> not the place, not `λ`, not the caps [H3, H4].

---

## 8. Status, honestly

| item | status |
|---|---|
| a class row is `(72,108)` with `C = y`: `(a,b,t,κ)`, `v(P)`, `v(Q)`, `[P,Q]`, `v(F) = −5`, `M`, `N` and the forcing ODE all shared | **PROVED** (each is a formula in `(a,b,t,κ)`); **EXACT-CHECKED**, with the `(72,108)` row reproducing `PROOF` §2.2 including `q(y)`, `6630`, `204`, `238`, `30`, `−1/2` |
| `P_M = [u^n]H²/C^{2n−2}`, `(𝒞³)_M = [u^n]H³/C^{2n−3}`, `C` and the `D_j` free | **PROVED** (ring identity); mutation-controlled at `±1` |
| hence `(2.5.1)`'s exponents attach to any place with `mult_β(C) = 1`, so the class rows' `y`-place inherits the condition set verbatim | **PROVED** |
| the cascade at the class row's `y`-place, levels 2–12, `ord_y h_k ≥ 2k−1`; Lemma 6.1's 12, 13 | **PROVED**, consuming exactly `[QQ1]` + `[I3]` as §6.1 does; **EXACT-CHECKED** by recomputation from scratch, with the perfect-square / one-fresh-coefficient shape reproduced |
| **correction:** `WEIGHT_FREE_TRANSFER` rows #5, #6 ("GENUINELY LOST") are wrong | **established** |
| polygon caps `ord_y h_k ≥ k`, `deg_y h_k ≤ 3k`, `λ = 2` at corner `(5,20)` | **EXACT-CHECKED** from `polygon_reduction`'s **computed** reduced polygon, with the `(8,28)` sub1/sub2 slopes `(12,15)`/`(12,14)` recovered as the control. Inherits `polygon_reduction`'s standing judgments (the rule `l = ceil(b0/a0)` is INFERRED there, though externally controlled by GGV3 §5 at this corner) |
| the same ord slope 1 at `(8,32)`, `(9,36)`, `(10,40)`; deg slopes `15/4, 4, 17/4` | **INFERRED** — same construction, but no in-repo *computed* polygon at those corners. The ord side is what the kill needs and it is uniform |
| **correction:** `λ = 0` / "no affine degree cap" is an artefact of reading both slopes off a monomial `Φ` | **established** |
| `ord_y(e) = 9` exactly, `e = γy⁹`, `B = y^21/γ` at every class row | **PROVED** given the transferred profile and `WEIGHT_FREE_TRANSFER` B1's K-syzygy |
| the three-route enumeration is complete | **EXACT-CHECKED** (independent re-derivation; confirms `WEIGHT_FREE_TRANSFER` H3) |
| **all three routes are FALSE on an explicit point of the slice families + polygon caps** | **EXACT-CHECKED.** Scope: they are not consequences of that input. A proof using `G₅ = 0` itself is not excluded — but §7.4(c), which now transfers, says the missing ingredient is leading-coefficient non-vanishing, not valuations |
| Cor 8.5 fires at `(72,108)` for every `(deg A, z)`, `z ∈ [2,6]` | **EXACT-CHECKED** — the calibration control, and it also reproduces §8.6's own `z ≤ 7` / cap-7 sensitivities |
| Cor 8.5 fails at a class row; the §8.1 `k = 0` system is satisfiable within every transferred cap | **EXACT-CHECKED**, four witnesses, residual 0 |
| "the class rows are nonempty" | **NOT CLAIMED and NOT SUPPORTED.** The §8 witnesses are points of the §8.1 reduced system with the caps, not germs; §§4's transfer to germs is not attempted here |

### Negatives, plainly

* **No case is closed. Eight rows stay open, and this document makes them harder,
  not easier**: it removes the target the previous lane left (all three routes)
  and shows the natural continuation (§8) is also blocked.
* **The good news and the bad news are again the same fact.** `mult_y(C) = 1` is
  what makes the cascade transfer *and* what makes `Φ`'s residual factor trivial
  (the forcing ODE's right side is `C^{b−a+1} = y²`). One cannot have the cascade
  without losing `q`.
* **`WEIGHT_FREE_TRANSFER`'s numbers were right; its classification was not.**
  Its §3 measurement (one unit at `w = 1,2,7`) is confirmed [D4]; what is wrong is
  #5/#6 "lost", G1/G2 `λ = 0`, and — as a consequence — the framing of §4's routes
  as the cheap way in. They were the cheap way in only because the cascade was
  believed absent; with the cascade present they are one unit past a bound that is
  *attained*.
* **The `λ ≥ m` gate of `CONTACT_LEMMA` is asked of the wrong object.** The `λ`
  that gate needs is the polygon's stripped slope, which is 2 here, not the
  `deg Φ/M − ord Φ/M` of `window_functions`, which is 0. Whether the gate then
  *passes* is not settled here — §3 bypasses it by deriving the slice families
  directly — but the recorded reason for its failure is not a reason.

### What would unblock it, precisely

Two targets, both now sharply stated.

1. **A degree-4 substitute for `q` at `k = 0`.** Cor 8.5 needs one contributor of
   degree `> deg u cap` in `γu = μπ³Q_Π − 6A² + 3ζπ^z`. At a class row the only
   candidates are `−6A²` (even degree, cancellable against `3ζπ^z` exactly when
   `2 deg A = z ≤ 6`) and the caps themselves. So: **either** lower the class
   row's `deg d₂` cap below 3 — impossible if the polygon caps are attained, and
   §2.6 says attainment is what forbids lowering a cap — **or** find a second
   equation that forbids `2 deg A = z`. The witnesses of §7 show which pairs must
   be excluded: `(deg A, z) ∈ {(0..3)} × {2..6}` and `A = 0`.
2. **Leading-coefficient information at `[y^{21}]B`**, per §7.4(c) — and here the
   class rows have something `(72,108)` does not: `e = γy⁹` and `B = y^{21}/γ` are
   *forced monomials*, so **every** coefficient of `B` above `y^{21}` vanishes,
   not just the bottom one. That is 15 equations on leading coefficients that
   `(72,108)`'s `t`-place cannot state. This is the one lead this document opens
   rather than closes, and it is untouched.

---

## Files

| file | role |
|---|---|
| `YPLACE_TRANSFER.md` | this writeup |
| `yplace_transfer.py` | the checker — 57/57, `--quiet` exit 0, ~3 min (`--fast`: 55/55, ~8 s) |
| `WEIGHT_FREE_TRANSFER.md` | the target audited here: §4's routes confirmed-then-refuted, rows #5/#6 and G1/G2 corrected |
| `PROOF_72_108.md` | §§2.1–2.6, 5, 6.1–6.3, 7.1–7.5, 8.1–8.7 — every number quoted at its point of use |
| `polygon_reduction.py` | the computed reduced polygons (`case_8_28` = the published control, `case_f2(0)` = the class-row corner) and the chart guard |
| `window_functions_75_125.py` | `window_law` / `family`; §4 corrects the `λ = 0` reading of its `(R3)` |
