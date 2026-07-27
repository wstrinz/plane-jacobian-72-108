# Is the bridge identity general, or only true on a slice?

**Answer: GENERAL. `205` is confirmed independently at `(8,28)/(3,4)/144`, and the
weaker of the two load-bearing formulas is now PROVED rather than inferred — by a
two-line local argument at `y=0` that nobody had made. The soft spot flagged in
`MONOMIAL_WINDOW_LAW.md` §6 as "the one place this writeup could be wrong in a
way that matters" is closed. v0.4.0 needs no erratum.**

**Checker:** `bridge_generality.py` — **59/59**, `--quiet` exit 0, ~55 s.
**Prerequisites:** `polygon_reduction.corner_chart_data` (the retraction guard —
the *only* source of chart data here), `corner_atlas.json` (the 34 rows, re-read),
`corner144_verify.py` (the `(550,205,69,276)` cross-check),
`c_series_75_125.py` §2 (the slice-sum lemma this generalises),
`monomial_window_law.py` / `MONOMIAL_WINDOW_LAW.md` (the identity under test).
**Sources of truth:** the forcing ODE re-derived here by direct bracket
differentiation over 144 `(a,b,t)` triples; `f` by a *fully generic* linear solve;
`N` from the *built* D-transform tower; chart data always through the guard.

Notation as in `MONOMIAL_WINDOW_LAW.md`: `t` chart exponent, `κ = t−2`,
`(a,b) = sorted(m,n)`, `q := ord_y C`, `a0 := deg_y C`, `dg := a0 − q`,
`e := b−a+1`, `coef := t(b−a)+κ+1`, `M := t(a+b) − (κ+1)`, `H := q(a+b) − 1`,
`ρ := ord_y(f)`, `α := ord_y(Φ)`.

---

## 0. Headline

| # | statement | status |
|---|---|---|
| **A** | `ord_y(Φ) = 205` at `(8,28)/(3,4)/144`, derived with **neither** `ρ = q(b−a)+1` **nor** `N = a·M−2b`. Equals the bridge prediction `a·q·M − H = 3·3·25 − 20`. | **EXACT-CHECKED** [F1–F6] |
| **B** | `ρ = q(b−a)+1` holds at **every** corner with `κ ≥ 0`, `q ≥ 1`, `t ≠ q(κ+1)`, and for **every** residual `g` with `g(0) ≠ 0`. | **PROVED** [C1–C7] |
| **C** | The unique excluded locus is `t = q(κ+1)`; on the standard class `κ = t−2` its only integer point is `(t,κ,q) = (2,0,2)`. No published row is there (`min t = 3`). | **PROVED** [C4, C5, G6b, G6c] |
| **D** | `ρ` is **branch-independent**. `phi_f7.py`'s ramified-vs-complex-pair ambiguity really does move `mult_(y+1)` and the cofactor — it **cannot** move `ord_y`, hence cannot move the bridge or anything downstream. | **PROVED** + 4 branches at 2 corners [C7, D6] |
| **E** | `N = a·M − 2b` follows from the D-transform exponent `a·w − 1`, which is here **derived** from the `a`-th-root expansion of `P` (and is *exactly* attained, so `clear` is tight, not merely sufficient) plus additivity of an affine exponent. | **PROVED** [E1, E2] + bounded sweep [G7] |
| **F** | All **34** atlas rows (**15** distinct chart signatures), including all **six** non-monomial rows, now have an independently derived `ord_y(Φ)`, every one equal to `a·q·M − H`. `PHI_KNOWN` had **one** entry. | **EXACT-CHECKED** [F7–F10] |

---

## 1. Why the previous confirmation was a cross, not a region

`MONOMIAL_WINDOW_LAW.md` §6 states the problem exactly: `ρ` and `N` were
confirmed at `(72,108)` (`q=7`, `b−a=1`) and along the F2 rungs (`q=1`, `b−a`
varying). Plot those in the `(q, b−a)` plane and you get a horizontal line plus a
vertical line — **a cross through one point**. Any formula of the shape
`ρ = q(b−a) + 1 + λ·(q−1)(b−a−1)` fits *every* previously confirmed point, for
*any* `λ`. That is the precise content of the flag, and it is a real risk: `λ ≠ 0`
would move `ord_y(Φ)` at exactly the six non-monomial rows and nowhere else — i.e.
at exactly the rows where the downstream consequences live.

So the test has to move both coordinates at once, and it has to do so without
using either formula.

---

## 2. The decisive test, and what makes it independent

`(8,28)/(3,4)/144` is the right target for the reason the brief gives: it is a
different `(m,n) = (3,4)` at the **same** corner `A0 = (8,28)` as the closed
`(8,28)/(3,2)/108`, so the corner geometry is held fixed while `(a,b,q)` moves
(`q: 7 → 3`, `a+b: 5 → 7`).

Chart data comes **only** from `polygon_reduction.corner_chart_data(8,28,4,3)`:
`t=4, κ=2, deg C=8, ord C=3`, retraction **holds**. Re-read from
`corner_atlas.json`, not assumed [A3].

**Step 1 — the ODE is a construction fact, not a formula.** Differentiating the
bracket directly with symbolic `c, f`:

```
[ x^(at) c^a , x^s f/c^b ]  =  x^kappa * a { t c f' - [t(b-a)+kappa+1] c' f }
```

verified over **144** `(a,b,t)` triples, `s = κ+1−a·t` [B1]. At the target this is
`12 c f' − 21 c' f = c^2` with `c = y^3(y^5+1)` [B3].

**Step 2 — `f` by a fully generic linear solve.** Eighteen unknowns `f_0..f_17`
— no `y^ρ` ansatz, no `g^e` shape, no monomial assumption. `sympy.linsolve`
returns a **unique** solution (zero free coefficients) [D1]:

```
f  =  -(1/15) y^4 (y^5+1)^2 ,        ord_y(f) = 4   READ OFF the polynomial.
```

This reproduces `corner144_verify.py`'s `f` exactly [D4], and the same machinery
reproduces the published `(72,108)` `f = −y^8(y+1)^2 q(y)/6630` [D3] and
`phi_f7.py`'s F7 `f` [D5] — so the solver is calibrated against three
independently-derived polynomials before being trusted anywhere new.

**Step 3 — `N` from the built tower.** `d_w := c_w · c^(a·w−1)`, and the `u`-slice
`M` of `S^b` is `c`-homogeneous with exponent `clear = a·M − b`, so
`N = clear − b`. Read off the built tower: `M = 25`, `clear = 71`, **`N = 67`**
[E2, E4]. The tower's own input `a·w − 1` is not quoted from the paper's `a=2`
recurrence — it is derived in §4 below.

**Step 4 — combine, on actual polynomials.**

```
ord_y(Phi)  =  ord_y(f) + N * ord_y(C)  =  4 + 67*3  =  205 .
Phi = f * C^67 = -(1/15) y^205 (y^5+1)^69 ,  deg 550,  mult_(y+1) 69.
```

and the bridge identity predicts `a·q·M − H = 3·3·25 − 20 = **205**` [F2–F6].

> **`205` CONFIRMED.** Independent derivation `==` bridge prediction, and the
> expanded `Φ` reproduces `corner144_verify.py`'s signature `(550,205,69,276)`.

Nothing in steps 1–4 reads `ρ = q(b−a)+1` or `N = a·M−2b`. The mutation controls
confirm the test has teeth: `ρ → ρ±1` and `N → N±1` each break the identity at
**all nine** derived corners [MUT A, MUT B].

---

## 3. `ρ` is PROVED — the local argument

This is the real result. `ord_y(f)` is a purely **local** invariant at `y=0`, and
the ODE localises cleanly. Write `c = y^q g` with `g(0) ≠ 0`. Then
`c' = q y^(q−1) g + y^q g'`, so dividing the forcing equation by `y^(q−1)` and
reading the `y^K` coefficients gives a **triangular** recursion,

```
a * sum_{i=0..dg} g_i [ t(k-i) - coef(q+i) ] f_{k-i}  =  [g^e]_{k+q-1-q*e} ,
```

whose pivot (the `i=0` term) is `a·g_0·(t·k − coef·q)`. The right-hand side first
becomes nonzero at `k + q − 1 − q·e = 0`, i.e. at

```
k = q(e-1) + 1 = q(b-a) + 1 =: rho_0 .
```

**The two-line proof.** By induction `f_k = 0` for every `k < ρ₀`, *unless* some
`k < ρ₀` kills the pivot. But

```
t*k = coef*q   and   k < rho_0
    <==>   t | q(kappa+1)   and   q(kappa+1) < t                       [C1, C2]
```

and `q(κ+1) ≥ 1` for `κ ≥ 0, q ≥ 1`, so that asks for a **positive multiple of
`t` smaller than `t`** — impossible [C3, 1728 swept points; G6, 4032 swept
points]. Hence `f_k = 0` below `ρ₀`, and at `k = ρ₀`:

```
a * g_0 * (t*rho_0 - coef*q) * f_{rho_0}  =  [g^e]_0  =  g_0^e ,
t*rho_0 - coef*q  =  t - q(kappa+1)   identically ,                     [C1]
==>   f_{rho_0}  =  g_0^(e-1) / [ a (t - q(kappa+1)) ]   !=  0 .
```

So `ord_y(f) = q(b−a)+1` **exactly**, at every corner, for every `g` with
`g(0) ≠ 0` — checked symbolically with a free `g` at all nine derived corners
[C7] and on 48 sampled abstract points [G8b].

### Three corollaries that matter

1. **The excluded locus.** `t = q(κ+1)` kills the pivot at `k = ρ₀` itself, and
   then the recursion has **no** solution. On `κ = t−2` this reads `q(t−1) = t`,
   whose only integer point is `(t,κ,q) = (2,0,2)` [C4]. `min t` on the 34 rows
   is `3`, so no published row is affected [C5]. The abstract sweep **locates**
   this locus rather than asserting it is empty — 17 `(t,κ,q)` triples off the
   standard class, exactly one on it [G6b, G6c].

2. **`ρ` reads `q` and `b−a`, and nothing else.** Not `t`, not `κ`, not `deg C`.
   Pinned from both sides: `q → q+1` moves `ord_y(f)` to `(q+1)(b−a)+1` at 5/5
   corners [MUT G]; `dg → dg+2` (so `deg C` moves, `q` fixed) leaves it
   **unchanged** at 5/5 [MUT H]. A check that moved under MUT H would be secretly
   reading `deg C`.

3. **The excluded locus is minus the Bezout corner integer.** `t − q(κ+1)` is the
   negative of `q(κ+1) − t`, the integer that `MONOMIAL_WINDOW_LAW.md` §2 shows
   controls `gcd(M,H)` [C6]. The *same* integer governs the window arithmetic and
   the solvability of the forcing ODE. At a monomial corner it is `+1`, so `(R)`
   is safest exactly where the class of nine lives.

---

## 4. `N` is derived, not extrapolated

`c_series_75_125.py` §2 already derives `N` from the built tower via the slice-sum
invariant, and its `(108,144)` control row already lands `205`. What was still
quoted rather than derived was the tower's own input, the D-exponent `a·w − 1`
(the paper's clearing recurrence, generalised in `a`). It follows from the
`a`-th root:

```
P = x^(at) c^a ( 1 + sum_{w>=1} pi_w z^w ) ,   pi_w = p_w/c^a ,  z = x^-1 ,
C = P^(1/a) = x^t * c * ( 1 + sum pi_w z^w )^(1/a) .
```

The `z^w` coefficient of the root is a sum of products of **at most `w`** factors
`pi_{w_i}` with `Σ w_i = w`, so its `c`-denominator is `c^(a·w)` at worst; times
the overall `c` that is `c^(a·w−1)`. Checked for `a = 2..5`, `w = 1..5`: the
exponent is **exactly** `a·w − 1`, not merely at most [E1]. Exactness matters —
it makes `clear` tight, so `N` is the true exponent and not an over-clearing.

Then the slice-sum is one line: `d_{w_1}···d_{w_b}` at `Σ w_i = M` carries
`Σ (a·w_i − 1) = a·M − b`, **independent of the individual `w_i`**. So the slice
is `c`-homogeneous, `clear = a·M − b`, `N = clear − b = a·M − 2b` [E2, E3].

Both inputs are shown load-bearing: perturbing the D-exponent to `a·w−1±1` moves
`N` at 18/18 mutants and breaks the bridge every time [MUT C]; making it
**non-affine** destroys `c`-homogeneity of the slice at 9/9 corners [MUT D] — so
the invariant is exactly "affine exponent + additivity", not an artefact of how
the tower is coded.

**What stays CLAIMED.** That `Φ` *is* the cleared slice-`M` object of `S^b` at
`M = b·t + j`. That is a construction identification, unchanged in status by this
file, and it is the honest residual premise on the `N` side. This file settles the
**arithmetic** of the bridge, not the geometry of `Φ`.

---

## 5. The joint direction — tested, with an honest qualification

Five corners move `q` **and** `b−a` together, all independently derived [G1, G4]:

```
 corner            q   b-a   dg   ord_y(Phi)   a*q*M - H
 F_7 (42,147)      4    5     2       165         165
 F_8 (63,147)      5    4     1       371         371
 F14 (66,231)      4    5     5       165         165
 F15 (99,231)      5    4     4       371         371
 F16 (99,165)      7    2     2       407         407
```

Both `dg` parities are represented — odd `dg` (forced `g = y^dg+1`, simple root at
`−1`) and even `dg` (the ramified branch) — so the branch structure is exercised
[G5].

> **NEGATIVE, stated plainly.** Those five are **COLLINEAR**: all have
> `q + (b−a) = 9`. This is *forced*, not accidental. All five sit at `t=3`,
> `k=1`, `p = q+l`, where GGV5's Diophantine `(m+n)qk − n(ql−p) = k` reduces to
> `q(n−m) = 3n−1`, hence `b−a = (3n−1)/q`; the rows with `q+(b−a) = 9` are exactly
> `n = (9q−q²+1)/3` [G3, G3b]. **GGV5's `v11 ≤ 35` tables contain no off-line
> joint corner.** So "the strongest available joint test" is genuinely
> constrained by the published population, and I could not do better *with real
> corners*.

Two things close the gap anyway:

* **The nine-corner tested set is affinely 2-dimensional.** `(3,1), (7,1), (8,1),
  (4,5), (5,4), (7,2)` span the plane — no line contains them [G3c]. The cross of
  §1 has become a region.
* **The abstract sweep fills the rectangle.** `q ∈ 1..8` × `b−a ∈ 1..6`, all 48
  combinations, across `t ∈ 2..8` and four `κ` regimes — **4032** points, and the
  pivot condition that §3 proves equivalent to `(R)` holds at every one [G6, G8].
  These points are not all GGV5 corners; but the object in doubt was a
  **formula**, and a formula is exactly what an abstract sweep can settle.

And the last step of the bridge is a **symbolic identity**, not a fit:
`(q(b−a)+1) + (a·M−2b)·q ≡ a·q·M − H` in all five symbols [G9]. So the chain is
two derived inputs plus one exact algebraic step — there is no fitting anywhere.

**The strongest joint test that remains undone** is an off-line joint corner,
which needs a corner outside GGV5's `v11 ≤ 35` length-1 `A0'=(1,0)` class:
either an `A0'=(2,0)` family (F12, F13) or a length-2 chain (F18–F24). Both need
a composite reduction chart that exists in no paper (`PHI_CORNER4.md` §2). Given
§3's proof, that test would now be corroboration rather than evidence — the proof
covers the whole plane, so no further point can be decisive.

---

## 6. All 34 atlas rows, independently derived

The 34 rows carry only **15** distinct chart signatures [A2]. Every one is derived
here — generic ODE solve + built tower — and every one agrees with `a·q·M − H`
[F7]:

```
  t  kap  a0   q  (a,b)   dg    M     H     N  ordPhi  bridge  rows
  3    1   1   1  (2,3)    0   13     4    20      22      22    10
  3    1   1   1  (2,5)    0   19     6    28      32      32     1
  3    1   1   1  (3,4)    0   19     6    49      51      51     2
  3    1   1   1  (3,5)    0   22     7    56      59      59     1
  3    1   1   1  (5,7)    0   34    11   156     159     159     1
  3    1   6   4  (2,7)    2   25    35    36     165     165     1   <- F_7  joint
  3    1   6   5  (3,7)    1   28    49    70     371     371     1   <- F_8  joint
  3    1   9   8  (2,3)    1   13    39    20     169     169     1   <- F_17
  3    1  12   8  (2,3)    4   13    39    20     169     169     1   <- (12,33)
  4    2   1   1  (2,3)    0   17     4    28      30      30     8
  4    2   1   1  (3,5)    0   29     7    77      80      80     1   <- (75,125)
  4    2   8   3  (3,4)    5   25    20    67     205     205     1   <- THE TEST
  4    2   8   7  (2,3)    1   17    34    28     204     204     1   <- CLOSED
  5    3   1   1  (2,3)    0   21     4    36      38      38     2
  6    4   1   1  (2,3)    0   25     4    44      46      46     2
```

All **six** non-monomial rows are covered [F8]; `204` and `205` are both
reproduced [F9], and `204` has the further disjoint corroboration of the `f1`-ODE
route (`AT_LE9_AUDIT.md` B7). The 28 monomial rows are the `dg = 0` signatures,
where `C = y` and `ord_y(Φ) = (b−a+1) + N` [F10].

This replaces `PHI_KNOWN`'s single entry with 15 signatures. Two label-integrity
guards are worth noting because this repo has been burned by exactly that: the
degree recipe `max(m,n)·(a0+b0) = max_deg` is re-checked on all 34 rows [A4], and
**MUT F** shows the numbers are *sensitive* to which chart dictionary is used —
stale data moves `ord_y(Φ)` from `51 → 205` at F1, `30 → 112` at F3, `22 → 107`
at F9. These are not label-agnostic arithmetic that would pass on anything.

---

## 7. Blast radius: none. What this means for v0.4.0

Because `205` is **confirmed**, the three shipped consequences stand, and their
status is now *stronger* than when they shipped:

| consequence | was | now |
|---|---|---|
| (i) Φ-divisor carry obstruction TOTAL at all 28 monomial rows | rested on INFERRED `ρ`, `N` | `ρ` **PROVED**; `N` proved modulo the unchanged Φ-identification. Unaffected. |
| (ii) 31 UNKNOWN verdicts resolved at atlas gate `G5` | same | **all 31 survive.** `G5` consumes `q_window = M/gcd(M,H)`, which by the identity equals `denom(α/M)`; the identity is now checked at all 15 signatures. |
| (iii) class-wide closure of the K-syzygy route | same | **survives.** It consumes `gcd(α,M) = gcd(M,H) = 1` at monomial corners, and monomial corners are `q=1`, `dg=0` — the *safest* region for `(R)` (Bezout integer `+1`, §3 corollary 3). |

**No erratum is required for v0.4.0.** The public release's claims are unchanged;
what changed is that the load-bearing lemma behind them moved from INFERRED to
PROVED. `MONOMIAL_WINDOW_LAW.md` §6's flagged row and §9's "unresolved" bullet are
now discharged and are annotated in place.

One additional strengthening falls out for free: because `ρ` is
**branch-independent** [D6], the branch ambiguity that `phi_f7.py` flags as "a
judgment item" cannot reach `ord_y(Φ)`. Anything downstream that consumes only
`ord_y(Φ)` — which is all of (i), (ii), (iii) — is therefore *insensitive* to that
judgment item. Only `mult_(y+1)` and the cofactor are exposed to it.

---

## 8. Negatives and open items, stated plainly

* **No new closure.** Nothing here kills a case. This is a hardening of a lemma
  that many things stand on, plus one genuinely new proof.
* **The joint corners GGV5 publishes are collinear** (§5). Real limitation of the
  population, not of the method; dissolved by the proof, but worth knowing if
  someone later wants a purely empirical joint check.
* **`Φ` = the cleared slice-`M` object of `S^b` stays CLAIMED** (§4). This is the
  honest residual premise on the `N` side and this file does not touch it. So the
  status of `N` is "proved *given* that identification", whereas `ρ` is
  unconditional.
* **The extreme-ray premise of `window_functions_75_125.py` is untouched** and
  remains the premise for "carry ≥ 1 ⟹ no Φ-divisor relation". Nothing here
  strengthens or weakens it.
* **LIVE INCONSISTENCY FOUND, unrelated to the bridge.** `phi_corner4.py` and
  `phi_f7.py` use **pre-repair** chart data at F1, F2, F3, F5, F9, F10: they take
  `t = l` from GGV5's table and apply the final-corner dictionary at corners that
  do **not** retract, where `polygon_reduction.corner_chart_data` returns the
  MONOMIAL data instead [A5b]. This is not cosmetic — `phi_corner4.py`'s VERDICT
  claims five reproduced points and **three** of them (`(50,75)`, `(75,125)`,
  `(56,84)`) sit at refused corners. At `(50,75)` it implies `ord_y(Φ) = 75`
  where the repaired route and `corner_atlas.json` both give `30`.
  `PASSPORT_75_125_REPAIR.md` lists `phi_corner4.py` as swept; the survey table
  and the two `derive(...)` calls were not. **Nothing in this file rests on
  those rows** — every number here goes through the guard — and the two files'
  F7/F14/F16 points are fine because `(6,15)` and `(9,24)` do retract. But the
  files should not be cited for the six stale rows, and someone should decide
  whether to repair or banner them. Flagged, not fixed: it is another lane's file.

---

## Files

| file | role |
|---|---|
| `BRIDGE_GENERALITY.md` | this writeup |
| `bridge_generality.py` | the checker — 59/59, `--quiet` exit 0, ~55 s |
| `MONOMIAL_WINDOW_LAW.md` §6 | the flag this discharges |
| `corner144_verify.py` | the `(550,205,69,276)` cross-check for `205` |
| `c_series_75_125.py` §2 | the slice-sum lemma generalised in §4 |
| `polygon_reduction.py` | `corner_chart_data`, the retraction guard used throughout |
| `phi_f7.py` | the branch ambiguity §3 corollary 2 neutralises for `ord_y` |
| `phi_corner4.py` | **stale chart data at six rows** — see §8 |
