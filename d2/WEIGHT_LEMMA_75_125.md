> ## !! SUPERSEDED IN PART -- 2026-07-26 REPAIR !!
>
> The `(5,20)` corner data this document consumes is **wrong**, and every number
> derived from it below has moved.  The chart exponent is `l = ceil(b0/a0) = 4`,
> **not** the denominator `5` of GGV5's final chain corner `(7\5,2)`: the
> dictionary `(t, q) = (l_final, b_final)` is valid only on the **retraction
> shape** `b0 = l*(a0-1)`, and `(5,20)` fails it (`20 != 4*4 = 16`).
>
> | quantity | this document | CORRECT |
> |---|---|---|
> | `l` = `t` | 5 | **4** |
> | `kappa` | 3 | **2** |
> | `C` | `y^2(y^3+1)`, `deg C = 5` | **`y`**, a monomial, `deg C = 1` |
> | `q = ord_y C` | 2 | **1** |
> | `N` (75,125) | 98 | **77** |
> | `Phi` (75,125) | `-(1/9) y^201 (y^3+1)^101` | **`(1/3) y^80`** |
> | signature | `(504,201,101,202)` | **`(80,80,0,0)`** |
> | `N` (50,75) | 36 | **28** |
> | signature (50,75) | `(189,75,38,76)` | **`(30,30,0,0)`** |
> | `q_window` | `5a-3` (7, 12) | **`12a-7`** (17, 29) |
>
> Decisive external evidence: GGV3 `1406.0886` sec.5 (`paper_src/
> 1406.0886_GGV3.tex:1723-1727`) performs this very reduction on the sibling
> `(50,75)` and publishes `[P_1,Q_1] = x^2`, `deg(P_1) = 10`, `deg(Q_1) = 15`.
> `l = 4` reproduces all three; `l = 5` contradicts all three.
>
> The root cause is now **guarded**: `polygon_reduction.final_corner_dictionary()`
> raises off the retraction shape.  Read **`PASSPORT_75_125_REPAIR.md`** for what
> survived, what changed, and what is undetermined, and treat the corresponding
> `*.py` (which has been repaired) as authoritative over this prose.

# WEIGHT_LEMMA_75_125 — does the Φ-divisor mechanism transfer to (75,125)?

2026-07-25. Checker: `weight_lemma_75_125.py` (**38/38**, `--quiet`, exit 0; ~48 s;
`--fast` truncates the deep sweeps). Read-only — no existing file touched.

## Verdict up front

| question | answer |
|---|---|
| Does a Φ-divisor relation `c·Φ = e·B` exist at (75,125)? | **NO.** Complete graded search of the built G-system ideal: nullity **0** at every weight `W = 36..50`, for `e = dm1` and for each of `d0,d1,d2,d3`. |
| Does the forcing fire? | **VACUOUS** — there is no relation, so the weight lemma has no hypothesis to discharge. `σ = 0` predicts nothing. |
| Is `σ = 0` at (75,125) corroboration? | **No — it is TAUTOLOGICAL.** `window_functions_75_125.window_law` *defines* `deg_slope := deg_y(Φ)/M`, so `σ = deg_slope·M − deg_y(Φ) ≡ 0` for every case by construction. |
| Is there a structural reason? | **Yes, and it is on the ORD side, not the deg side.** The quasi-affine lower cap has quasi-period `q_window = 12`; `e = dm1` has u-weight `6`, and `L(6)+L(30) = 202 > 201 = ord_y(Φ)`. The mechanism needs `q_window \| w_e`. At (72,108) `q_window = 1`, so that is free. |

**Net: the divisor mechanism is NOT family-level. (72,108) is special, and the
thing that makes it special is `q_window = 1`** — the integrality of
`W_step = ord_y(Φ)/M`, which `G_SYSTEM_75_125.md`'s own correction block already
called "a friendly coincidence of that corner, not generic `a=2` behavior."

Two independent routes agree: the algebraic search (§3) finds no relation, and
the cap arithmetic (§5) says a relation is exactly what the carry forbids.

---

## 1. The (75,125) weight data, re-derived from primitives

`weight_lemma_75_125.py` §A. Nothing below is imported as a stored value except
the generator strings themselves (which are then re-checked, A4).

**Case parameters** `(a,b,t,κ,q_mult) = (3,5,5,3,2)`, corner `(5,20)→(7/5,2)`,
`C = y²(y³+1)` (A0).

**Φ, re-solved from the forcing ODE** (A1–A2). The operator identity
`a{t·c·f′ − [t(b−a)+κ+1]·c′·f} = c^(b−a+1)` instantiates to `15·c·f′ − 42·c′·f = c³`
(A1a). A 15-variable linear solve over `deg ≤ 14` returns a **unique** polynomial
solution (A1b), `f = −(1/9)·y⁵(y³+1)³` (A1c). With
`N = a[t(a+b) − (κ+1)] − 2b = 98` (A2a),

```
Φ = f·C^98 = −(1/9)·y^201·(y³+1)^101
signature (deg_y, ord_y, mult_(y+1), cofactor) = (504, 201, 101, 202)   (A2b,c)
```

matching `PHI_75_125.md` exactly, re-derived here rather than read.

**The u-grading, from `S = Σ_m d_m u^(t−m)`** — so `w(d_m) = t − m` (A3):

| symbol | `d3` | `d2` | `d1` | `d0` | `e = dm1` | `dm2..dm10` | `Φ` |
|---|--:|--:|--:|--:|--:|--:|--:|
| u-weight | 2 | 3 | 4 | 5 | **6** | 7 … 15 | **36** |

`w(Φ) = M = b·t + jφ = 36` (A2d). Every one of the ten generators is
**re-checked u-homogeneous** from its stored string, weight `= b·t + j` (A4a–c),
giving the AP `26,27,28,29,30,31,32,33,34,36`; and **Φ occurs only in `G11`, with
coefficient exactly 1** (A4d — the standing stale-`2Φ` guard, transplanted).

**The two slopes** (A5):

```
W_step    = ord_y(Φ)/M = 201/36 = 67/12     q_window = 12 = 5a−3   NON-integral
deg_slope = deg_y(Φ)/M = 504/36 = 14        = 5a−1                 integral
lambda    = deg_slope − W_step = 101/12                            NOT an integer
```

> **The (72,108) stripping factor has no (75,125) counterpart.** At (72,108) the
> stripped slope is `λ = 14 − 12 = 2` (sub2), an integer, which is what makes
> `deg_stripped ≤ λ·w` additive over monomials. Here it is `101/12` (A5c). The
> lemma therefore **cannot be run in stripped coordinates at all** at (75,125).
> It *can* be run in **full y-degree coordinates**, because `deg_slope = 14` is
> an integer and `U(w) = 14w` is affine — which is the version §4 uses.

**Contrast with (72,108), including the disjoint-route point** (A6). From the
pinned Prop 4.3 sub2 polygon, `max(j − 2i)` over the hull `= 0` (A6a), whence
`deg C_{4−k} ≤ 8 − 2k` and, after the D-transform,
`deg D_{4−k} ≤ (8−2k) + 8(7−2(4−k)) = 14k` — **slope 14 from the polygon alone,
with Φ never consulted** (A6b). Independently `deg_y(Φ)/M = 238/17 = 14` (A6c).
Two disjoint computations agreeing. At (75,125) there is **no such polygon**:
`upstream_facts.json` carries the `(8,28)` corner only (A6d), so the (75,125)
`14` has *only* the Φ-defined provenance. This is §4's point.

---

## 2. Does an analogue of the K-syzygy exist? — the method

This question comes **first**, and it is decided by a complete linear-algebra
search of the graded ideal, not by numerology.

**Setup.** A relation `c·Φ = e·B` that is an exact ideal identity means
`c·Φ − e·B = Σ h_i G_i`. Everything in sight is weight-homogeneous, so each
weight component is separately such an identity; fix weight `W`. Reduce mod `e`:
writing `h_i = h_i⁰ + e·(…)` with `h_i⁰` free of `e`, and `c = c⁰ + e·c′`,

```
c⁰·Φ  ≡  Σ h_i⁰ G_i   (mod e),      h_i⁰  weight-homogeneous of weight W − w(G_i).
```

Since `Φ` occurs only in `G11` and with coefficient 1, `c⁰ = h_11⁰`, and the
Φ-free part must vanish mod `e`. So the search is:

> **columns** = `(Φ-free part of m·G_i)` reduced mod `e`, one column per
> generator `G_i` and per monomial `m` of weight `W − w(G_i)` **not divisible by
> `e`** (monomials divisible by `e` are automatically `≡ 0 mod e`, contribute
> only to `B`, and constrain nothing);
> **nullity 0 ⟺ every solution has `c⁰ = 0` ⟺ `c ∈ (e)` ⟺ the relation is vacuous.**

**Rigour of the rank test.** Ranks are computed mod the prime `p = 2⁶¹−1`.
`rank_p ≤ rank_ℚ` always, so **`rank_p = #columns` PROVES nullity 0 over ℚ**. The
weight-36 case is additionally confirmed by an exact rational nullspace (B3a).

**The control, and its sensitivity.** The identical engine on (72,108) at
weight 17 returns **exact nullity 1** (B1a) with the recovered relation

```
G5 + d2·G3 + d1·G2 + d0·G1        (B1b — the published K-syzygy, uniquely)
2·(that)  ==  2Φ − e·(d2·e² + 3·e·S + 3·R²)     (B1c — residual exactly 0)
```

and at weight `17+k` the kernel dimension equals **exactly** the number of
weight-`k` monomials, `k = 0..9` (B2) — i.e. the kernel is precisely the
`K`-multiples, nothing more, nothing less. The engine is exactly sensitive.

Four mutation tests (run out-of-band, not in the checker):

| mutation | expected | observed |
|---|---|---|
| (72,108) `G3 += dm1³` (an `e`-multiple) | syzygy survives | nullity 1 ✅ |
| (72,108) `G3 += d0·d1·dm4` (`e`-free) | syzygy destroyed | nullity 0 ✅ |
| (75,125) `G11body := dm1·dm10·dm5·d0` (`e`-multiple) | relation appears | nullity 1 ✅ |
| (75,125) `G11body := dm2·dm10·dm5·d1` (`e`-free) | none | nullity 0 ✅ |

So the (75,125) zero below is a real zero, not an engine that cannot find things.

---

## 3. The result: NO relation exists

`weight_lemma_75_125.py` §B, on the canonical `g_system_75_125.json` (whose own
verifier `g_system_75_125_verify.py` passes, exit 0, re-run here).

**Divisor `e = dm1`** (B3):

| `k = w(c)` | 0 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| `W = 36+k` | 36 | 38 | 39 | 40 | 41 | 42 | 43 | 44 | 45 | 46 | 47 | 48 | 49 | 50 |
| columns | 37 | 64 | 82 | 108 | 138 | 177 | 223 | 283 | 352 | 441 | 543 | 673 | 822 | 1008 |
| **nullity** | **0** | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | **0** |

(`k = 1` is skipped: there is **no** monomial of u-weight 1, so `c = 0` is forced
outright there.)

**Other divisor candidates** (B4): the same search with `e` replaced by `d0`,
`d1`, `d2` or `d3`, weights 36..46 — **nullity 0 throughout**. No state variable
divides `c·Φ` on this ideal.

### 3b. Why — the structural reason, and a concrete witness

The (72,108) K-syzygy is not an accident of that ideal; it comes from a slice
identity. With `S = Σ d_m u^(t−m)`, `d_t = 1`, `d_{t−1} = 0`, the linear window
`S^a` kills `[u^j]S^a` for `j = a·t+1 … a·t+((b−a)t+jφ)`, skipping one. At
(72,108) (`a=2, b=3, t=4`) that makes `[u^17]S⁴ = [u^17](S²·S²) = 0`, and
expanding `S⁴ = S·S³` gives

```
Σ_k d_{t−k}·[u^(17−k)]S³ = 0 .
```

The `k = 1` term dies because `d_3 = 0` (which is also why the *skipped*
generator never appears). The **head** `k = 0,2,3,4` is `G5, d2·G3, d1·G2, d0·G1`;
the **tail** `k ≥ 5` is `e·[u^12]S³ + R·[u^11]S³ + …`, and — this is the crux —
**those tail slices lie BELOW the generator range**, which starts at `u^13`. So
the identity survives modulo the ideal as a genuine statement about `Φ`, and the
tail turns out to be `e`-divisible.

At (75,125) the same identity exists — `b+1 = 6 = 2a`, and `[u^36]S⁶ = 0` for the
same reason (the only surviving cross-term pairs the skipped `[u^35]S³` with
`[u^1]S³ = 3d_t²d_{t−1} = 0`). But the head/tail split lands differently:

```
(G11 − Φ) + d3·G9 + d2·G8 + d1·G7 + d0·G6 + e·G5 + dm2·G4 + dm3·G3 + dm4·G2 + dm5·G1
          + [ tail below the generator range ]  =  0
```

The tail slices `[u^26..u^30]S⁵` **are** the generators `G1..G5`. The identity
therefore collapses into the ideal and isolates nothing: modulo the ideal it says
`Φ ≡ Σ_{k≥11} d_{t−k}[u^(36−k)]S⁵`, with no `e` in sight.

**Concrete witness** (B5): the `e`-free head combination
`G11 + d3·G9 + d2·G8 + d1·G7 + d0·G6 + dm2·G4 + dm3·G3 + dm4·G2 + dm5·G1`, with
`Φ` removed and `e → 0`, leaves a **nonzero 149-term residue**. At (72,108) the
same construction's residue is exactly `0` — that *is* the K-syzygy.

---

## 4. What `σ = 0` at (75,125) actually says: nothing

`weight_lemma_75_125.py` §C1–C2. Mechanically the number is right:

```
σ = deg_slope·W − D = 14·36 − 504 = 0 .
```

But `window_functions_75_125.window_law(ordPhi, M, degPhi)` sets
`deg_slope := degPhi/M`. So `σ = (degPhi/M)·M − degPhi ≡ 0` **identically, for
every case, by construction**. It is not a measurement. The CAPS_AUDIT entry
that logged this as `[inference]` was right to; the sharper statement is that
`σ = 0` here is not merely premise-driven but *definitionally empty*.

The content the number *would* have is: **is `14` the true upper window-degree
cap slope at (75,125)?** At (72,108) that question has an independent answer
(A6a–b: the Prop 4.3 polygon, Φ never consulted), and the agreement `14·17 = 238`
is real. At (75,125) the corresponding polygon does not exist in any source
(`PHI_75_125.md` judgment 2, A6d) — so there is nothing for `deg_y(Φ)/M` to agree
*with*.

**Conditional, and moot** (C2): *if* a relation existed *and* `14` were the true
cap, the forcing would fire and pin `deg_y(e) = 14·6 = 84` (stripped `50`). §3
says there is no relation, so this pin is never reached.

---

## 5. The ord-side obstruction — a carry corollary to the weight lemma

The weight lemma as stated in `CAPS_AUDIT.md` §5 uses only the **upper** (degree)
cap, and that cap is affine at (75,125). The thing that breaks is the **lower**
(y-order) cap, which is quasi-affine — and it breaks the relation's *arithmetic*,
not its degree bookkeeping.

Under the inherited extreme-ray premise, `L(w) = ceil(αw/q)` with `α = 67`,
`q = 12`, `gcd(α,q) = 1`, and `L(36) = 201 = ord_y(Φ)` exactly (C3a). Then:

> **Carry corollary.** Let `c·Φ = e·B` with `c` a nonzero constant and
> `w_e + w_B = W`, with `q \| W` and `ord_y(Φ) = L(W)`. On any lift
> `ord(e) + ord(B) = L(W)`, while the caps give
> `ord(e) + ord(B) ≥ L(w_e) + L(W − w_e) = L(W) + carry(w_e, W−w_e)`. Hence
>
> ```
> carry(w_e, W − w_e) = 0 ,   and since gcd(α,q)=1 this holds iff  q | w_e .
> ```
>
> *(Elementary: write `αw_e/q = A + r/q`; `r = 0` gives carry 0, `r ≠ 0` gives
> carry 1, and `r = 0 ⟺ q | w_e`.)*

At (75,125), `W = M = 29`, `w_e = t+1 = 5` (**REPAIRED 2026-07-26**; the values
below read `W = 36`, `w_e = 6`, `α/q = 67/12` until then — the superseded `(5,20)`
chart, `2adb92a`):

```
L(5) + L(24) = 14 + 67 = 81  >  80 = L(29) = ord_y(Φ)           carry = 1   (C3b)
carry(w_e, 29−w_e) = 0  for  NO  w_e in 1..28                               (C3c)
```

**The repaired reading is strictly stronger.** `q_window = 29 = M` exactly, so no
split is a multiple of the period and the carry is 1 on *every* admissible split —
whereas the superseded `67/12` against `M = 36` would have let it vanish at
`w_e ∈ {12, 24}`.

`e = dm1` sits at `w_e = 5`. **The mechanism is ord-obstructed for the actual
`e`.** At (72,108) `q_window = 1`, so carry `= 0` for *every* split (C4) — the
ord side is exactly balanced there, which is the same fact recorded in
`DIVISOR_SYZYGY.md` §3b as `t^a | R,S,T`.

**Family consequence** (C5). The F2-family law `q_window = 5a − 3` gives
`7, 12, 17, 22, 27, …` — **never 1**. So the carry obstruction is *generic in the
family*, and (72,108)'s `q_window = 1` is exceptional to its corner. The
mechanism's real precondition is not "`σ = 0`" but "**`q_window` divides
`w_e`**", which at `q_window = 1` is invisible.

**Logical direction, stated carefully.** This corollary does *not* by itself
prove the relation absent; it proves that *if* the relation existed, (75,125)
would have no lifts at all. §3 is the load-bearing negative. The carry is (i) an
explanation of why one should not have expected the relation, (ii) an independent
consistency check that agrees with §3, and (iii) a genuine sharpening of the
weight lemma to quasi-affine regimes.

---

## 6. PROVED / CHECKED / INFERRED / CONDITIONAL

**PROVED (exact, unconditional given the stated inputs)**

* The graded-search reduction of §2: a weight-homogeneous ideal identity
  `c·Φ − e·B = Σ h_i G_i` reduces mod `e` to exactly the linear system solved.
  Nullity 0 ⟹ `c ∈ (e)` ⟹ the relation is vacuous.
* `rank_p ≤ rank_ℚ`, so full rank mod `2⁶¹−1` proves nullity 0 over ℚ. The
  weight-36 case is also proved by an exact rational nullspace.
* **No relation `c·Φ = e·B` (nor with `e` replaced by `d0,d1,d2,d3`) exists as an
  exact ideal identity in the built (75,125) G-system, for any multiplier `c` of
  u-weight ≤ 14.**
* The carry corollary of §5 (elementary number theory, given `gcd(α,q)=1`).
* The u-weight assignments, generator homogeneity, `w(Φ) = 36`, and
  `coeff(G11, Φ) = 1`, all recomputed from the stored generator strings.
* Φ's signature `(504, 201, 101, 202)` from the uniquely-solvable forcing ODE and
  `N = 98`.
* At (72,108): the K-syzygy is the **unique** weight-17 relation, and the kernel
  at weight `17+k` is exactly the `K`-multiples (`k = 0..9`).

**CHECKED (computed, not proved in general)**

* Sensitivity of the search engine: four mutation tests, §2 table.
* `deg_slope = 14` from the (72,108) Newton polygon by brute hull enumeration
  (A6a–b) — a re-derivation of an already-proved statement, not new.

**INFERRED**

* The *reason* for the negative (§3b: head/tail split of the `[u^{b+1}·t·…}]`
  slice identity). The slice identity `[u^36]S⁶ = 0` and the index arithmetic are
  exact, but the claim that this is the *only* source of a Φ-relation is an
  interpretation — §3's search is what actually rules relations out.

**CONDITIONAL ON A PREMISE**

* **Everything is conditional on the built (75,125) objects.** The corner
  `(5,20)→(7/5,2)`, `C = y²(y³+1)`, `Φ = −(1/9)y^201(y³+1)^101`, and hence the
  whole G-system, rest on the **standard unreduced-polygon reduction** of the
  `(5,20)` corner, which no paper performs (`PHI_75_125.md` judgment 2,
  `C_SERIES_75_125.md`). If that reduction is not the right one, the negative is
  a negative about the wrong object.
* **§5 is conditional on the extreme-ray premise** (`ord ≥ (α/q)w` for every
  window object, with equality on Φ) — `WINDOW_FUNCTIONS_75_125.md` §5. Note the
  direction: under the premise `L = ceil` is the *smallest* admissible integer
  cap, so a true cap can only be **higher**, and the carry obstruction only gets
  worse. The premise is therefore not load-bearing *against* the obstruction.
* `deg_slope = 14` at (75,125) is **Φ-defined**, not independently derived (§4).
  §3's negative does not use it; §4's conditional pin does.

---

## 7. What would falsify this

1. **A relation at higher weight.** The sweep covers `w(c) ≤ 14`. A relation with
   a multiplier `c` of u-weight ≥ 15 would falsify the headline. (The (72,108)
   relation has `w(c) = 0`; nothing suggests a high-weight one, but the search is
   finite.)
2. **A relation that holds on the variety but is not an ideal identity.** The
   search is complete for exact identities — the standard the K-syzygy meets — but
   not for membership in the radical or in a saturation `(I : e^∞)`. Settling that
   needs a Gröbner computation on 15 variables with 973-term generators, not
   attempted here. **This is the one gap in the negative.**
3. **A different `g_system_75_125.json`.** If the (5,20) reduction premise is
   wrong, the built ideal is wrong and everything here is about a different ring.
4. **A `w_e` with `q_window | w_e`.** If the relevant divisor at (75,125) were a
   window symbol of u-weight 12 or 24 rather than `e = dm1` at weight 6, §5's
   obstruction would lift. No such object plays `e`'s role in the state.
5. **An independent (75,125) window-degree cap ≠ 14.** That would kill §4's
   conditional pin (already moot) and would also disturb §5's `α/q`.

---

## 8. Reproduce

```
python -u weight_lemma_75_125.py            # full report, 38 checks, ~48 s
python -u weight_lemma_75_125.py --quiet    # exit 0 iff all pass
python -u weight_lemma_75_125.py --fast     # truncated sweeps (~25 s)
```

Files created by this lane: `weight_lemma_75_125.py`, `WEIGHT_LEMMA_75_125.md`.
Consumed read-only: `g_system_75_125.json`, `paper_src/upstream_facts.json`,
`full_system_bridge.py` (the (72,108) control generators). Nothing modified.
