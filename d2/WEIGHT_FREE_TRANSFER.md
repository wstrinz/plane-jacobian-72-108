# Which of the eleven machine steps are WEIGHT-FREE — and what that buys at the class of nine

**Answer: 3 transfer verbatim, 5 are mixed with a named transferable core, 3 are
genuinely lost. But the interesting number is not 3/5/3. It is ONE.**

The eight class rows with `(a,b,t) = (2,3,4)` inherit the entire *algebraic* spine
of `(72,108)` — the K-syzygy, the shift dictionary, the bracket collapse, the toric
syzygy — and, because `Φ` is a **monomial** there, they inherit a *strengthened*
Theorem 3.4 which collapses the five-case analysis to the single case `k = 0`.
What they do not inherit is the slice cascade. And the whole of the remaining gap
is **one unit of `y`-order at three specific weights**: supply it and eight of the
nine die at §3 in one line, with §§4–11 never invoked.

**Checker:** `weight_free_transfer.py` — **58/58**, `--quiet` exit 0, ~35 s
(the `(3,5,4)` build in §C dominates).
**Prerequisites:** `PROOF_72_108.md` §13.2 (the eleven steps), `MONOMIAL_WINDOW_LAW.md`
(the bridge identity, the total-carry lemma, the thinness ⊕ shallowness split),
`TORIC_GENERAL.md` Q1 (the `t = 4` resonance), `CONTACT_LEMMA.md` (the `λ ≥ m` gate).
**Sources of truth, nothing retyped:** generators from
`g_system_75_125.build_gsystem` (with `.published_72108` as the label control);
window arithmetic from `window_functions_75_125.window_law` / `.family`; every
`(72,108)` number quoted at its point of use with the `PROOF_72_108.md` section.

Notation as in `MONOMIAL_WINDOW_LAW.md`: `q := ord_y C`, `M := t(a+b) − (κ+1)`,
`α := ord_y Φ`, `q_window := M/gcd(M, α)`, and the window floor
`L(w) := ceil(α w / q_window)` (`window_functions_75_125` `(S1)`). The `u`-weights
are `w(h_k) = k`, so `w(d_2,d_1,d_0,e,R,S,T) = (2,3,4,5,6,7,8)` and `w(Φ) = M`.

---

## 0. The operational definition of "weight-free"

A step is **weight-free** iff its statement survives with **`Φ` a free symbol**.
That is not a slogan: `build_gsystem` forms the generators from `(a,b,t)` alone and
`Φ` enters additively in one slice, so an identity that holds with `Φ` free cannot
have consumed a valuation, a cap, a window, or `λ` — there is nothing for those to
attach to. Every WEIGHT-FREE row below is verified in that form, and verified at
**`(50,75)`'s own generators**, not at `(72,108)`'s.

Three guards make the label honest [A1–A5]:

* `build_gsystem(2,3,4,1,30)` reproduces `g_system_75_125.published_72108()` term
  by term — so the object called "(50,75)" really is the published `(72,108)`
  system;
* the `(50,75)` and `(72,108)` builds differ in **exactly** `{q, ordPhi, W_step}`
  — only the weight normalisation;
* moving `t` to 5 **does** change the generators, so the agreement is not an
  artefact of the comparison.

---

## 1. The eleven rows

| # | step | § | checker | class | what the weight-dependence consumes |
|---|---|---|---|---|---|
| 1 | the K-syzygy expansion | 3.1 | `divisor_syzygy.py` 7/7 | **WEIGHT-FREE** | nothing. A polynomial identity in `Z[d2,d1,d0,e,R,S,T,Φ]`. Re-derived here from the class row's own generators, `Φ` free, residual exactly 0 [B1], with two mutation controls [B1b, B1c]. Corollary 3.2 (ideal *equality*) transfers with it [B2]. |
| 2 | forcing-ODE solution, `q` irreducible | 2.2 | `at_le9_audit.py` B1–B9 | **MIXED** | **Core:** `Φ ≠ 0` — the only thing Lemma 3.3 consumes; at a class row `Φ = (1/2)y^30`, so the `e ≡ 0` branch is empty there too, by the same one line, char ≠ 2 [D1, D2]. **Wrapper:** the specific ODE and the quartic `q` (`q(−1) = 3315`, `disc ≠ 0`, `deg = 4`). The wrapper does not transfer and does not need to: the class row's replacement is the empty word [D3]. |
| 3 | shift triangularity at index `−1` | 2.3 | `window_caps_verify.py` W2/W3/W5 | **WEIGHT-FREE** | nothing. Generalized-binomial algebra in the chart exponent `t` alone: `C(m,m−j) = 0` for `m ≥ 0 > j`, hence `D̃₋₁ = D₋₁` with no `θ` [B5]. `t = 4` is shared by all nine. All five dictionary relations (7.1.1) are re-derived from the binomials [B6.*], and triangularity is confirmed to **fail** at index `−2` [B5b]. **Correction to the lane prior**, which did not flag this row: it is weight-free. |
| 4 | **the cap lemma** | 2.6 | `caps_audit.py` 70/70 | **MIXED — and this is the pivot** | Lemma 2.5 is *two* lemmas with opposite fates. **Ord half (core, transfers in form and gets STRONGER):** `ord D_{j_x} ≥ 48 − 12 j_x`. New cross-check: under `w = 4 − j_x` this is **identically** the window floor `L(w)` at `α/q_window = 12` — two repo objects that had never been compared [E2]. At a class row the same floor reads `L(w) = ceil(30w/17)`, strictly above the affine ray. **Deg half (wrapper, LOST):** `deg d_j ≤ λw` with `λ = 3/2`. At a class row `λ = 0` and `deg_slope = 30/17 ∉ Z`, so the affine degree cap **does not exist** [G1, G2]. *(Classified only; not deep-dived, per instruction. The one substantive remark: the paper states the two halves as one lemma, and the transfer statement is only expressible once they are separated.)* |
| 5 | slice-family cokernel ranks `2n−3` | 5.1–5.2 | `slice_obstruction_basis.py` S4 | **WEIGHT-DEPENDENT, LOST** | the two slice families (2.5.1) `t^{2n−2} | p_n`, `t^{2n−3} | r_n`, whose exponents are the `t`-place ord/deg data of `C_4 = y^7 t`, plus `CONTACT_LEMMA`'s gate `λ ≥ m`. At every class row `λ = 0 < 2 = m`: gate fails. The `2·3 − 3·2 = 0` cancellation itself is weight-free (`contact_lemma.py` A1–A5, symbolic in `(m,n)`) but has nothing to hang on. |
| 6 | cascade lowest jets to level 12 | 6.1 | `slice_obstruction_audit.py` 56/56 | **WEIGHT-DEPENDENT, LOST** | same gate — it is the level-by-level execution of #5. Its *output*, the profile (6.2.1), is exactly what the class rows need and do not have. See §3: this row is where the entire remaining gap lives, and its size is now measured. |
| 7 | the bracket collapse | 7.2 | `syzygy_collision.py` X6.1 | **WEIGHT-FREE** | nothing. A coefficient identity given #3's dictionary: `−3/8 + 3/16 + 3/16 = 0`. Verified here from the generalized-binomial dictionary, residual 0 [B7, B7b], with the `−3h₁h₅h₆` mutation control [B7c]. |
| 8 | **the valuation ledger** `v_t(d2,d1,R,S,T) ≥ (2,3,10,11,12)` | 7.5 | `sub1_spine9.py` P16 | **MIXED — and it gets STRONGER, not lost** | **Core:** the four `min`'s of Lemma 7.4 over the dictionary (7.1.1) — pure weight-free bookkeeping. **Wrapper:** the profile it is fed. Fed (6.2.1) the min's give the paper's `(2,3,10,11)` and `v_t(T) ≥ 12`; fed the class row's floor they give `(4,6,11,13)` and `≥ 15` — uniformly higher [F6, F6b]. **Correction to the lane prior**, which listed this as a weight-dependent suspect: it is weight-dependent, but the dependence is on an *input*, and the class row's input is a better one. |
| 9 | the `a = 9` divisions and `FZ = (1/6)γ⁵t⁹Π⁴` | 8.1, 8.3 | `sub1_spine9.py` P7 | **MIXED** | **Core:** Theorem 3.5, `6WZ = e⁵` — `Φ`-free, cap-free; §8.3 says outright that the cofactor identity *is* Theorem 3.5 in normalised coordinates and survives `t⁹ → t^a`. Verified at the class row's own generators, residual 0, with `Φ` absent from the combination [B3], plus §3.5's companion `W²` syzygy [B4]. The exponent 5 also transfers: `TORIC_GENERAL` Q1's condition `(t+1) | (4t+9)` has the unique solution `t = 4`, shared by all nine [D6]. **Wrapper:** the `t⁹/t¹⁸/t²⁷` normalisation (consumes `a_t = 9`) and `Π`. |
| 10 | the marked-support feasibility test (40 pairs) | 8.5 | `sub1_spine9.py` P12, `spine9_audit.py` E1–E7 | **MIXED — and VACUOUS at a class row, not lost** | **Core:** the rank-1 / `Π²`-divisibility criterion, pure arithmetic of `q`. **Wrapper:** the existence of `Π` and the `z`-window. At a class row `rad Φ = y` alone, so `e = γ y^n` and `Π = 1`: `k = 0` is **forced** and the four cases `k = 1,2,3,4` that this step kills **do not arise** [D4]. This is the strengthening of Theorem 3.4 that monomiality hands over free. |
| 11 | the degree ledger | 8.6 | `sub1_spine9.py` P13–P15 | **WEIGHT-DEPENDENT, LOST — and it is the one that matters** | every entry of the §8.6 table is a `λ = 3` cap (`deg A ≤ 9`, `deg u ≤ 6`, `deg v ≤ 12−k`, `deg w ≤ 9+k`), plus `(*deg)`. At `λ = 0` there is no table [G3]. Worse: the *only* case surviving #10's collapse is `k = 0`, Corollary 8.5 — a pure degree dichotomy with **three independent zero-margin dependencies**, all of them weight data [G4]. So the one case that arises at the class of nine is killed at `(72,108)` by exactly the machinery that has no counterpart. |

**Split: N = 3 weight-free (1, 3, 7); M = 5 mixed (2, 4, 8, 9, 10); K = 3 lost
(5, 6, 11).** Of the five mixed, #2 and #10 need no replacement at all (their
wrappers become trivial resp. vacuous), and #8's wrapper improves. The real
replacement list is two items long: #4's ord half, and — only if that fails —
#5/#6 plus #11.

---

## 2. Two upgrades from classification to result

### 2.1 The K-syzygy re-derived at `(50,75)`, and it is a `(2,3,4)` RESONANCE

`(50,75)`'s generators, built by `build_gsystem(2,3,4,1,30)`, satisfy

```
2*(G5 + d2*G3 + d1*G2 + d0*G1)  =  2*Phi  -  e*(d2*e^2 + 3*e*S + 3*R^2)
```

with residual exactly 0 and `Φ` a free symbol [B1]. So the whole of §3 —
Theorem 3.1, Corollary 3.2, Lemma 3.3, Theorem 3.4 — holds at all eight
`(a,b,t) = (2,3,4)` class rows. **PROVED** (it is a ring identity; the checker is
a re-derivation, not evidence).

But it is **not family-level**. An exhaustive `u`-homogeneous search for a
cofactor system realising `c·Φ = e·B` — every monomial of the forced weight
`M − (bt+j)` on every generator, `Φ`-coefficient normalised to 1 — returns:

* at `(2,3,4)`: **consistent**, unique solution, and the solution *is* the
  published `(d2, d1, d0)` [C1 — the positive control];
* at `(3,5,4)` = `(75,125)`: **inconsistent**, `rank A = 18 < rank[A|b] = 19` over
  21 unknowns and 241 equations. **No identity of that shape exists, at any
  cofactors** [C2].

So the K-syzygy joins `6WZ = e⁵` as a `(2,3,4)`-only resonance
(`TORIC_GENERAL.md` Q1 established the latter; nobody had tested the former).
**The class of nine is 8 + 1, not 9.** `(75,125)` inherits neither syzygy and
needs its own mechanism — though §E10 covers it by the total-carry lemma alone,
without any syzygy: every 2-split of `M = 29` has carry `≥ 1`.

### 2.2 Theorem 3.4 STRENGTHENS: the five cases collapse to one

At a class row `Φ = (1/2)y^30`, a monomial [D1]. Hence `rad Φ = y`, one linear
factor, against `rad Φ = t·q` of degree 5 at `(72,108)` [D3]. Theorem 3.4's
`e | 2Φ` therefore gives

```
(72,108) :  e = gamma * t^a * Pi ,  Pi | q squarefree ,  k = deg Pi in {0,1,2,3,4}
class row:  e = gamma * y^n       ,  Pi = 1            ,  k = 0  FORCED
```

[D4]. The four cases `k = 1,2,3,4` — the ones §8.5 and §8.6 kill by weight-free
arithmetic, with wide margins and no cascade input — **do not arise**. The one
case that does is `k = 0`, Corollary 8.5, which `PROOF_72_108.md` §8.7's own
per-case table marks as the *only* one that consumes the cascade, with zero margin
in three independent places [D5, G4].

**That is the bad news, stated precisely.** Monomiality deletes exactly the four
cases whose kills transfer and keeps exactly the one whose kill does not.

---

## 3. The measurement: `(72,108)`'s `t`-place is the class rows' `y`-place

This is the part that converts "lost" into a number.

`MONOMIAL_WINDOW_LAW`'s bridge identity `ord Φ = a·q·M − H`, evaluated at
`(a,b,t,κ) = (2,3,4,2)` — i.e. at `(72,108)` — at the **two places of `C_4 = y^7 t`**:

```
q = 7  (the y-place, ord_y C_4 = 7)  ->  204  =  ord_y(Phi)      [PROOF sec.2.2]
q = 1  (the t-place, mult_t C_4 = 1) ->   30  =  v_t(Phi)        [PROOF Lemma 2.3]
```

Both are the paper's own numbers [F1]. So **`(72,108)`'s `t`-place is a
monomial-type place**: thin (`mult_t C_4 = deg_t C_4 = 1`) and shallow, with

```
(alpha, M, q_window)  =  (30, 17, 17)   at (72,108)'s t-place
(alpha, M, q_window)  =  (30, 17, 17)   at the eight class rows' y-place
```

**identically** [F2]. The `t`-place at `(72,108)` is a faithful model of the class
rows' `y`-place. And `(72,108)` is the one corner in the repo where a
monomial-type place is backed not by the extreme-ray *premise* but by a *proof* —
the cascade. So we can measure how strong that proof is against the premise.

| `k` | cascade law `2k−1` | profile (6.2.1) | floor `L(k) = ceil(30k/17)` |
|---|---|---|---|
| 1 | 1 | 1 | **2** |
| 2 | 3 | 3 | **4** |
| 3 | 5 | 5 | **6** |
| 4 | 7 | 7 | **8** |
| 5 | 9 | 9 | 9 |
| 6 | 11 | 11 | 11 |
| 7 | 13 | 12 | 13 |
| 8 | 15 | 13 | 15 |

The cascade **law** `2k−1` and the floor **agree exactly at `k = 5,6,7,8`**, the
floor is higher by exactly 1 at `k = 1,2,3,4`, and the cascade is higher for
`k ≥ 9` (they cross at `4k = 17`) [F4]. The profile the paper actually records
uses Lemma 6.1's weaker 12, 13 at `k = 7,8`, so against the floor it is one short
at `k = 1,2,3,4,7`, equal at `5,6`, and two short at `8` — where `h₈` does not
occur in the bracket at all [F4b].

**And that one unit is the whole of §§4–11.** Under the bracket collapse (#7),
every monomial of `B` has `u`-weight 12 [E6], and:

```
profile (6.2.1) :  min_B = 21 ,  a_t + v_t(B) = 9 + 21 = 30 = v_t(Phi)   NO contradiction
floor           :  min_B = 22 ,       9 + 22 = 31 > 30                   CONTRADICTION
```

[F5, F5b — the profile row reproduces `PROOF_72_108.md` §7.3's table
`(≥21, ≥21, ≥21, ≥22)` exactly]. So:

> If the window floor held at `(72,108)`'s `t`-place, `(72,108)` would die at
> **§3, in one line**: `2Φ = eB` is impossible on valuation grounds, and §§4–11
> — the cascade, the collision, the five cases, the caps, the support test, the
> degree ledger — would never be written.

The authors needed thirty pages instead. That is the honest measure of how much
stronger the extreme-ray premise is than the strongest thing this repository has
ever *proved* at a monomial-type place: **one unit, at four or five weights.**

Same arithmetic at the class rows' `y`-place, where it is the actual kill
[E7, E8, E9]:

```
(72,108) y-place, q_window = 1 :  L(5) + L(12) = 60 + 144 = 204 = ord_y Phi   carry 0
class row y-place, q_window = 17:  L(5) + L(12) =  9 +  22 =  31 > 30         carry 1
```

and all four bracket monomials land on `22 = L(12)` exactly [E8].

---

## 4. The replacement ledger — what has to be proved, exactly

Enumerating over the profile (6.2.1) every set of weights raised by exactly `+1`,
and asking which sets make `ord(e) + min_B > 30`, gives **three minimal sets and
no others** [H1–H4]:

| route | what must be proved at the class row's `y`-place | relative to the floor |
|---|---|---|
| `{5}` | `ord(e) ≥ 10` | one unit **past** the floor (`L(5) = 9`) |
| `{1,2,7}` | `ord(h₁) ≥ 2`, `ord(h₂) ≥ 4`, `ord(h₇) ≥ 13` | **exactly** the floor at three weights, profile-strength at the other two |
| `{2,6,7}` | `ord(h₂) ≥ 4`, `ord(h₆) ≥ 12`, `ord(h₇) ≥ 13` | floor at two, one past it at `w = 6` |

`B` only ever sees weights `1, 2, 5, 6, 7`, so this list is complete. Any one of
the three routes empties **eight of the nine** at §3.

The `{1,2,7}` route is the cheapest and is **strictly weaker than "prove the whole
floor"**: it needs the floor at 3 of the 5 relevant weights and nothing at all
beyond the profile at the other 2 — and at `w = 5, 6` the profile *is* the floor
already [H4].

If none of the three is available and only a cascade analogue (profile strength)
can be built, the class rows land on Corollary 8.5, the `k = 0` degree dichotomy —
which needs `deg d₂ ≤ 6`, i.e. an affine degree cap, which at `λ = 0` **does not
exist** [G2, G4]. That branch is a dead end, not a longer road.

---

## 5. Status, honestly

| item | status |
|---|---|
| K-syzygy holds at the eight `(2,3,4)` class rows, from their own generators, `Φ` free | **PROVED** (ring identity); **EXACT-CHECKED** with 2 mutation controls |
| Ideal equality (Cor 3.2), toric syzygy (Thm 3.5), §3.5's `W²` companion, all at the class row | **PROVED / EXACT-CHECKED** |
| Shift triangularity at `−1` and all five dictionary relations (7.1.1) from generalized binomials in `t` | **PROVED** (symbolic), **EXACT-CHECKED** |
| Bracket collapse under that dictionary | **PROVED**, with mutation control |
| **No K-syzygy at `(3,5,4)`** — the search is inconsistent, with a working positive control at `(2,3,4)` | **EXACT-CHECKED**. Scope: identities of the shape `c·Φ = e·B` that are `u`-homogeneous with polynomial cofactors. A relation of some *other* shape is not excluded. |
| `Φ = (1/2)y^30` at `(50,75)`, hence `Π = 1`, `k = 0` forced, five cases → one | **EXACT-CHECKED** — but it inherits `window_functions_75_125`/`c_series_75_125`'s F2-family derivation of `Φ`. At `(8,32)`, `(9,36)`, `(10,40)` the value `ord_y Φ = 30` comes from the bridge, whose `ρ`/`N` generality is **INFERRED** (`MONOMIAL_WINDOW_LAW` §6, the flagged soft spot). So "five cases collapse to one" is EXACT-CHECKED at `(50,75)` and **INFERRED** at the other seven. |
| `ord D_{j_x} ≥ 48 − 12j_x` **is** the window floor `L(w)` at `α/q = 12` | **EXACT-CHECKED** (new cross-check between `caps_audit` §2.6(iii) and `window_functions` `(S1)`) |
| `(72,108)`'s `t`-place has window data `(30,17,17)`, identical to the class rows' `y`-place | **PROVED** from the bridge, using the bridge at `q = 1`; both endpoints (`204`, `30`) are the paper's own numbers |
| cascade law `2k−1` = floor at `k = 5..8`, one below at `k = 1..4` | **EXACT-CHECKED** |
| "the floor would kill `(72,108)` at §3" | **EXACT-CHECKED arithmetic**, and it is a *conditional*: it does not refute the floor premise, because `(72,108)` is in fact empty. What it establishes is that the premise is strictly stronger than the cascade. |
| the three minimal upgrade routes are complete | **EXACT-CHECKED** over the `+1` model on the weights `B` sees. It is a model: a real upgrade might raise several weights at once, which only helps. |
| classification of steps 4, 5, 6, 11 as consuming `λ` / the caps | **INFERRED** from reading §§2.6, 5, 6, 8.6 plus the `λ = 0` / non-integral `deg_slope` facts, which are **EXACT-CHECKED**. No `caps_audit.py` code was modified or deep-dived. |

### Negatives, plainly

* **No case is closed here.** Every transferred step is a *hypothesis-side* result:
  it says the class rows carry the same ideal, not that the ideal is empty.
* **`(75,125)` is worse off than it looked.** It shares neither syzygy. The
  "class of nine" is really 8 + 1 for every purpose in this document.
* **The good news and the bad news are the same fact.** `Φ` monomial deletes the
  four transferable kills (`k = 1,2,3,4`) and keeps the one non-transferable kill
  (`k = 0`). Monomiality is not neutral about *which* case survives.
* **The `ρ`/`N` generality gap propagates here.** Seven of the eight rows'
  `ord_y Φ = 30` rests on it. `(50,75)` does not.

---

## 6. Files

| file | role |
|---|---|
| `WEIGHT_FREE_TRANSFER.md` | this writeup |
| `weight_free_transfer.py` | the checker — 58/58, `--quiet` exit 0 |
| `PROOF_72_108.md` §13.2 | the eleven steps being classified |
| `MONOMIAL_WINDOW_LAW.md` | the bridge identity, total-carry lemma, thinness ⊕ shallowness |
| `TORIC_GENERAL.md` | the `t = 4` resonance for `6WZ = e⁵`; §2.1 here is its K-syzygy counterpart |
| `g_system_75_125.py` | `build_gsystem`, `published_72108` |
| `window_functions_75_125.py` | `window_law`, `family`, `(R1)`–`(R3)`, `(S1)` |
