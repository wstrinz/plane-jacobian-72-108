# POSITIVE_SLICE.md — the positive-slice obstruction, and the emptiness of the last standard-sub2 cell

2026-07-25. Checkers: `positive_slice.py` (63/63, `--quiet`, exit 0),
`positive_slice_verify.py` (79/79, `--quiet`, exit 0 — separately authored, and
it imports nothing from `positive_slice.py`). Read-only on every existing
artifact: this lane wrote only `POSITIVE_SLICE.md`, `positive_slice.py`,
`positive_slice_verify.py`, `positive_slice_stage.json`. Pure sympy — no
Singular, no msolve, no WSL, no subprocess. Nothing was written to
`state_kill_ledger.json`, `proof_dag.json`, `phase_d_states*.json`, or
`frontier_rebuild.py`.

> **HEADLINE.** `SPINE.md` reduces standard sub2 to ONE cell, `a10_b0000_T1`
> (19 flag cases / 828 states `rl`; 21 / 850 `norl`). That cell survives because
> the surviving family satisfies the four canonical G rows **identically** — no
> Gröbner engine can close it, because the ideal is genuinely non-empty.
>
> The missing equations are not in the G-system at all. The `d3`-killing shift
> that *defines* the G-system's coordinates is **rational in `y`**, so a solution
> of the shifted system need not reconstruct a polynomial `P = C^2` on the Newton
> polygon. Undoing the shift and demanding that the positive-`x` slices of `P`
> come out polynomial gives three conditions. On `a10_b0000` they have **no
> common solution in any characteristic-zero field**, by a resultant over `Q`.
>
> **`a10_b0000_T1` is EMPTY. Standard sub2 is EMPTY.** No large elimination, and
> no field-scope caveat.
>
> **Two things a reader should attack first.** (i) The whole result rests on
> `t^a | dm2,dm3,dm4` (§7 [Q6]), which `SPINE.md` §8 itself flags for
> second-party adjudication — it is not re-proved here. (ii) The identification
> of the G-system indeterminates with the *shifted* window variables (§7 [Q8]);
> §3.3 gives four corroborations, but it is a convention, not a theorem proved
> in this file.

---

## 0. What the reviewer proposed, and what survived contact

The external review proposed the obstruction and the three equations (A),(B),(C).
This lane was asked to **derive**, not to accept. The outcome:

| claim under test | verdict |
|---|---|
| `P_i = y^(2i-2) [u^(8-i)] H(u)^2 / t^(14-2i)` | **DERIVED** from the repo's own D-transform (§2). Reproduces all nine slices `P_0..P_8` exactly on genuine data. |
| `D2* = d2 + (3/8)h^2`, `D1* = d1 + (1/2)h d2 + (1/16)h^3`, `D0* = d0 + (1/4)h d1 + (1/16)h^2 d2 + (1/256)h^4` | **DERIVED** from the *same* transformation `window_caps_verify.py` W3 uses (§3). All three correct as stated. |
| shift → unshift recovers the original positive slices | **CONTROL PASSES** (§4), 4 + 3 independent random polygon-supported instances across the two checkers. |
| the three constant terms are (A),(B),(C) | **CONFIRMED, with a correction of emphasis** (§5.3): (A) is exact; (B) and (C) are *not* the raw conditions — they are the raw conditions reduced modulo (A) and then **divided by `X` and `X^2`**. That division is precisely why `X = 0` must be a separate horn. |
| the system is empty over every char-0 field, `res(p,q) = 561971200` | **CONFIRMED**, twice, by two disjoint mechanisms (§6). |
| any two of the three conditions suffice | **NO.** All three are load-bearing; explicit witnesses for every pair (§8). |

Nothing in the reviewer's algebra was found to be wrong. The one substantive
correction is the (B)/(C) clearing factors, which the horn analysis already
covered but which must be stated, because they are what makes `X = 0` a case.

---

## 1. The problem the G-system cannot see

The window variables are the coefficients of the Laurent square root
`C = sum_{j<=4} c_j x^j` of `P = C^2`, in the D-transform of
`window_caps_verify.py` W2:

```
c_j = D_j * C4^(2j-7) ,     D_4 = 1 ,     C4 = y^7*(y+1) = y^7*t
```

and the bridge's *stripped* coordinate is `d_j = D_j / y^(48-12j)` (weight `12k`,
`k = 4-j`; `full_system_bridge.WEIGHT`).

The G-system is not written in those coordinates. It is written **after** the
`d3`-killing shift

```
x -> x - s ,      s = c_3/(4*C4) = D_3/(4*C4^2)                 (W3)
```

which is what makes `d3` disappear (`generators.json`'s `variable_order` is
`d2, d1, d0, dm1..dm4, Phi` — there is no `d3` row, and that is not an accident).

**`s` is rational in `y`.** The shift is therefore invertible on coefficient
*sequences* but is not a change of variable that preserves "comes from a
polygon-supported polynomial". A point of `V(G1,G2,G3,G5)` is a shifted
coefficient sequence; to be a counterexample germ it must un-shift to an honest
`C` whose square has polynomial slices. That is an extra, non-empty demand, and
it is invisible to any Gröbner basis of the G-system ideal.

---

## 2. The slice formula, derived

`P = C^2` gives, with no shift anywhere,

```
P_M = sum_{i+j=M} c_i c_j = C4^(2M-14) * sum_{i+j=M} D_i D_j
```

Put `H(u) := sum_{j<=4} d_j u^(4-j)`. Since `D_j u^(4-j) = d_j (y^12 u)^(4-j)`,

```
[u^(8-M)] H(u)^2 = y^(-12(8-M)) * sum_{i+j=M} D_i D_j
```

and with `C4 = y^7*t` the two exponents collapse:

```
y :   7*(2M-14) + 12*(8-M)  =  2M-2
t :   2M-14
```

> ```
> P_M  =  y^(2M-2) * [u^(8-M)] H(u)^2 / t^(14-2M)
> ```

`P_M` is a polynomial, and `gcd(y, t) = 1`, so for `M <= 6`

> ```
> t^(14-2M)  |  [u^(8-M)] H(u)^2
> ```

At `M = 7, 8` the exponent is `<= 0` — no condition. At `M = 6, 5, 4` the
conditions are `t^2 | [u^2]H^2`, `t^4 | [u^3]H^2`, `t^6 | [u^4]H^2`, and those
slices involve **only** `d_j` with `j >= 0`: no window spare `dm1..dm4` enters
the obstruction at all (`positive_slice.py` P1.5).

*This is the whole reason the obstruction is cheap.* It touches `d2, d1, d0` and
one new quantity, and nothing else.

`M <= 3` gives further conditions (and at `M = 0` a `y^2` condition as well).
**They are not used.** Only the three above are needed.

---

## 3. The inverse shift, derived

### 3.1 The general map

`window_caps_verify.py` W3's transformation is

```
X_j  =  sum_{m=j..4} binom(m, m-j) * src_m * theta^(m-j)
```

which is exactly the coefficient map of `x -> x + theta` (checked against a
literal substitution on a generic degree-4 polynomial: `positive_slice.py` P2a,
`positive_slice_verify.py` V2.3). The forward, `d3`-killing shift is this map at
`theta = -D_3/4`; it sends `D_3 -> D_3 + 4*(-D_3/4) = 0`.

The inverse is therefore the **same map at `theta = +D_3/4`**, and nothing needs
to be hand-coded. `positive_slice.py` proves it by applying the map twice
(P2b.2); `positive_slice_verify.py` proves it structurally, from the group law of
the 9×9 shift matrices

```
M(a) * M(b) = M(a+b) ,       M(0) = I                            (V2.1, V2.2)
```

### 3.2 The literal formulas

Writing `h := D_3` (the ORIGINAL, pre-shift stripped `D_3`) and
`(D~_4, D~_3, D~_2, D~_1, D~_0) = (1, 0, d2, d1, d0)`:

> ```
> D3* = h
> D2* = d2 + (3/8)*h^2
> D1* = d1 + (1/2)*h*d2 + (1/16)*h^3
> D0* = d0 + (1/4)*h*d1 + (1/16)*h^2*d2 + (1/256)*h^4
> ```

exactly as the review states. Hence

> ```
> [u^2]H^2 = 2*d2 + (7/4)*h^2
> [u^3]H^2 = 2*d1 + 3*h*d2 + (7/8)*h^3
> [u^4]H^2 = 2*d0 + (5/2)*h*d1 + d2^2 + (15/8)*h^2*d2 + (35/128)*h^4
> ```

The formulas are scale-invariant under stripping (`y^(48-12j)` cancels on both
sides), so they hold verbatim in the bridge's stripped coordinates.

### 3.3 `h` is a genuinely new unknown — the one convention premise

`h` is the quantity the shift **destroys**. It cannot be recovered from the
G-system, and the G-system imposes nothing on it. What the window does impose is
its degree: `deg d_3 <= 8 - 2*3 = 2`. In the controls `h(-1)` takes the values
`4, -11/2, -5, 11/2` — it is not forced to vanish.

The obstruction below uses **only** `eta := h(-1)`, one free scalar. So "no `eta`
works" is the same as "no `h` works", and the argument does not depend on `h`'s
degree at all.

That the G-system's `(d2,d1,d0,dm1..dm4)` are the **shifted** variables `D~_j` is
premise **[Q8]**. It is a convention, corroborated four ways:

1. `full_system_bridge.py` states it outright: *"G1,G2,G3,G5body =
   (D~^3)_{-1,-2,-3,-5} after the (D~^2) linear substitutions"* — tilde, explicitly.
2. `generators.json`'s `variable_order` has **no `d3`** (P0.3, V0.2). Only the
   shifted system is missing that variable.
3. `window_caps_verify.py` W3 exists precisely to show the caps survive the
   shift, and W4 re-checks every cap "on the unshifted AND shifted variables" —
   the bridge's caps are the shifted ones.
4. The stripped caps match on the nose: `deg d_j <= 8-2j` gives
   `4, 6, 8, 10, 12, 14, 16` for `j = 2..-4`, exactly
   `full_system_bridge.STRIP_DEGCAP` plus the `d2/d1/d0` rows (P3.2, V1.5).

---

## 4. THE POSITIVE CONTROL (non-negotiable — and it passes)

Both checkers independently generate polygon-supported `P` over `Q` — corners
loaded from `paper_src/upstream_facts.json`, never transcribed; `P_8 = C4^2`
forced; every other slice random inside the hull — and then run the full
round trip.

| control | `positive_slice.py` | `positive_slice_verify.py` |
|---|---|---|
| Laurent square root | quadratic D-recursion (`verify_derivation` §C) | **global binomial series** `sum_n binom(1/2,n) w^n` |
| slice identity | stripped `y/t` form | **unstripped** `P_M*C4^(14-2M) = sum D_i D_j` |
| shift inverse | apply the map twice | **9×9 shift-matrix group law** |
| instances | 4 seeds | 3 seeds |

Results (P3.1–P3.9, V1.2–V1.7, V2.6–V2.9), all exact, all instances:

* the D-recursion / binomial series agree, and both attain the certified window
  caps `ord >= 12k`, `deg <= 14k` on every row `k = 0..8`;
* stripping by `y^(48-12j)` is legal and lands exactly on the bridge's caps;
* **`P_i = y^(2i-2)[u^(8-i)]H^2/t^(14-2i)` reproduces all nine slices
  `P_0..P_8` exactly**;
* the forward shift kills `D~_3`;
* **shift then unshift returns the original stripped `D`'s exactly**;
* **the original positive slices are recovered exactly from the unshifted-back
  data** — in the verifier, through the unstripped convolution, so the stripped
  form is not reused;
* the three divisibilities **hold** on genuine data. They are necessary
  conditions, not vacuous ones.

Both checkers **abort with exit 1** if any control fails, before reaching the
obstruction. Neither does.

---

## 5. The `n = 0` family, and the three equations

### 5.1 Re-derived from `generators.json`

`spine.py` is **not** imported. `G5 = G5body + Phi` is rebuilt and
`coeff(G5, Phi) == 1` is asserted (the stale-`2*Phi` guard). With `n = 0`
(`Rm = 1`, `a = 10`, `Q = q`) and `dm1 = gamma*t^10`, `dm2 = t^10*A`,
`dm3 = t^10*B`, `dm4 = t^10*C`, `Phi = c*t^30*q`, unassisted `sp.div` gives
`t^20 | G1,G2,G3` and `t^30 | K := 2*(G5 + d2*G3 + d1*G2 + d0*G1)` exactly, and
the quotients are exactly SPINE's `g1, g2, g3, kbox` (residual 0; P4.2).

`d1` is **never set to zero**, so everything below covers T1 and T2 alike.

Then, all with `gamma != 0`:

```
g1 = 0   =>   C  = -A*(d2 + B/gamma) - (1/2)*gamma*d1
g2 = 0   =>   d0 = (d2*A^2 + 2*A*C + B^2)/gamma^2        (deg d0 <= 8 never used)
```

and the elimination certificate, with the cofactor **produced**, not quoted:

> ```
> F*Z - (1/6)*gamma^5*t^10  =  gamma^2 * g3hat ,
> F := A*(gamma*d2 + 2B) + (1/2)*gamma^2*d1 ,   Z := A^2 - gamma*B
> ```

(`positive_slice.py` P4.4. `positive_slice_verify.py` V3.4 instead **reduces**
`F*Z - (1/6)gamma^5 t^10` to 0 modulo a Gröbner basis of `<g1,g2,g3>`, supplying
no cofactor at all.)

### 5.2 Degree exactness, then `y = -1`

From the certified caps `deg A <= 2`, `deg B <= 4`, `deg d2 <= 4`, `deg d1 <= 6`:

```
deg F <= 6 ,   deg Z <= 4 ,   6 + 4 = 10 = deg(t^10)     EXACTLY
```

Zero slack. The right-hand side is nonzero, so both factors are nonzero, both
caps are attained, and — being polynomials whose product is a constant times
`t^10` — each is a constant times a power of `t`:

```
Z = zeta*t^4 ,   F = phi*t^6 ,   phi*zeta = gamma^5/6 != 0
```

so **`Z(-1) = 0` and `F(-1) = 0`**. Together with `kbox(-1) = 0` and the `g2`
solve, the values at `y = -1` are **uniquely** forced, with `alpha := A(-1)` and
`gamma` free (P4.7; V3.8 solves the 5×5 system and finds a unique solution, and
V3.9 confirms the remaining row `g3` is then automatically satisfied):

```
beta   = alpha^2/gamma
delta2 = -(6*alpha^2*gamma + 1)/gamma^3
delta1 =  2*alpha*(4*alpha^2*gamma + 1)/gamma^4
delta0 = -alpha^2*(3*alpha^2*gamma + 1)/gamma^5
```

(`q(-1) = 3315`, `mu*q(-1) = -1/gamma`.)

### 5.3 The three equations

Substituting into §3.2 and taking the constant term (only the *weakest*
consequence of each divisibility is used — `t | ·`, not `t^2 | ·`, `t^4 | ·`,
`t^6 | ·`) gives, with `eta := h(-1)`:

```
N(6) =  7*eta^2*gamma^3 - 48*alpha^2*gamma - 8
N(5) =  128*alpha^3*gamma - 144*alpha^2*eta*gamma^2 + 32*alpha
        + 7*eta^3*gamma^4 - 24*eta*gamma
N(4) =  3840*alpha^4*gamma^2 + 2560*alpha^3*eta*gamma^3
        - 1440*alpha^2*eta^2*gamma^4 + 1280*alpha^2*gamma + 640*alpha*eta*gamma^2
        + 35*eta^4*gamma^6 - 240*eta^2*gamma^3 + 128
```

With `X = alpha^2*gamma`, `Y = alpha*eta*gamma^2` (the rewrite is verified by
substituting back, monomial by monomial):

```
N(6) * alpha^2*gamma^1 = 7Y^2 - 48X^2 - 8X                                = (A)
N(5) * alpha^3*gamma^2 = 128X^3 - 144X^2Y + 32X^2 - 24XY + 7Y^3
N(4) * alpha^4*gamma^2 = 3840X^4 + 2560X^3Y + 1280X^3 - 1440X^2Y^2
                         + 640X^2Y + 128X^2 - 240XY^2 + 35Y^4
```

> **(A) is reproduced exactly.**
> **(B) and (C) are the reductions modulo (A), after dividing out `X` and `X^2`:**
> ```
> E(5) mod (A)  =    16   * X   * (8X^2 - 6XY + 2X - Y)              = 16 X (B)
> E(4) mod (A)  =  (-64/7)* X^2 * (480X^2 - 280XY + 160X - 70Y + 11) = -(64/7) X^2 (C)
> ```

**This is the one place where the brief's presentation needs sharpening.** (B)
and (C) are consequences of the raw conditions *only away from `X = 0`*. The
reviewer's horn analysis handles `X = 0` and is therefore correct — but the
divisions are what make the horn necessary, and they should be stated. Both
checkers avoid the issue entirely in their primary route (§6.1, §6.2), which
never divides by anything.

---

## 6. Emptiness — two disjoint proofs

### 6.1 Raw, saturated (`positive_slice.py` P6.1)

The ideal `<N(6), N(5), N(4), w*gamma - 1>` in `Q[alpha, eta, gamma, w]` has
Gröbner basis `{1}`. No clearing, no horns, no coordinate change. A unit ideal
over `Q` means no solution over any field containing `Q` — in particular none
over `C`.

### 6.2 Resultants only, no Gröbner basis anywhere (`positive_slice_verify.py` §V5)

* **Horn `alpha = 0`.** Then `X = 0` and the `(X,Y)` picture is blind, so the raw
  equations are used. `N(6)` reads `7W - 8 = 0` with `W := eta^2*gamma^3`, so
  `W = 8/7`; `N(4)` is `35W^2 - 240W + 128` in the *same* `W`, and at `W = 8/7`
  that is **`-704/7 != 0`**. Contradiction. (Notably this horn needs only `P_6`
  and `P_4`.)
* **Main branch `alpha != 0`.** Then `X != 0` and the multipliers
  `alpha^p*gamma^r` are units, so the raw system is equivalent to the `(X,Y)`
  system. `lc_Y(E(6)) = 7`, a nonzero constant, so the resultant criterion is
  exact. Eliminating `Y`:

  ```
  res_Y(E6, E5) = -50176  * X^3 * (320X^3 + 160X^2 + 29X + 2)
  res_Y(E6, E4) = -200704 * X^4 * (307200X^4 + 204800X^3 + 42240X^2 + 2080X - 121)
  gcd            =  50176 * X^3                      -- a MONOMIAL in X
  ```

  A common solution therefore needs `X = 0`, i.e. `alpha = 0`, already refuted.

### 6.3 Why C08/C20 cannot reopen this

The horn-free route of the brief is also reproduced: off `X = 0` and
`6X + 1 = 0`, (B) gives `Y = 2X(4X+1)/(6X+1)`, and then

```
(A) -> -4X * p ,   p = 320X^3 + 160X^2 + 29X + 2
(C) ->        q ,   q = 640X^3 + 320X^2 + 86X + 11
resultant(p, q) = 561971200 != 0        (gcd(p,q) = 1 over Q)
```

`6X + 1 = 0` is refuted by (B) alone: `(B)|_{X=-1/6} = -1/9`, independent of `Y`.

**Every step is a polynomial identity, a gcd, or a resultant over `Q`.** No
square class is taken, no splitting field is entered, no "the discriminant is a
square" argument appears. A gcd of 1 over `Q` means no common root in **any**
extension. The `FIELD_SCOPE_AUDIT` / C08 / C20 downgrade — which converts KILLs
into CONSTRAINTs over a general characteristic-zero base field — has no purchase
here. (It can still reopen branches *elsewhere*; the net frontier count is not
this file's to quote. See §8.)

---

## 7. PROVED / CHECKED / INFERRED

**PROVED here** — an exact polynomial identity, a finite exhaustive enumeration,
or an exact resultant, machine-checked in char 0 by *both* checkers over the
premises below:

* the slice formula `P_M = y^(2M-2)[u^(8-M)]H^2/t^(14-2M)`, and hence the three
  divisibility conditions (§2);
* the inverse-shift formulas `D2*, D1*, D0*`, and the three slice polynomials
  (§3.2);
* the round trip: shift ∘ unshift = identity, and the recovery of the original
  positive slices from genuine polygon-supported data (§4);
* the `n = 0` row factorisations, the `C` and `d0` eliminations, the certificate
  `F*Z = (1/6)gamma^5 t^10` with its cofactor, and the degree exactness (§5.1–5.2);
* the unique `y = -1` values, and the three constant-term equations (§5.2–5.3);
* their identification with (A),(B),(C), including the clearing factors (§5.3);
* emptiness over every characteristic-zero field, twice (§6);
* the ablation, with explicit witnesses (§8).

**CHECKED** — reproduced from an existing artifact without re-proving it:

| tag | statement | source |
|---|---|---|
| [Q1] | the canonical `G1,G2,G3,G5body` | `generators.json` (loaded, never transcribed) |
| [Q2] | `Phi = c*t^30*q`, `c = -1/6630`, `q(-1) = 3315` | `verify_derivation.py` §A |
| [Q3] | window caps `ord >= 12k`, `deg <= 14k`; `D_j = C_j*C4^(7-2j)` | `window_caps_verify.py` W2 |
| [Q4] | the `d3`-killing shift and its D-coordinate form | `window_caps_verify.py` W3 |
| [Q5] | `e = gamma*t^a*Rm`, `a + deg Rm = 10`; `n = 0 => Rm = 1, a = 10` | `DIVISOR_SYZYGY.md` |
| [Q6] | **`t^a \| dm2, dm3, dm4` on both branches** | `SPINE.md` §8 |
| [Q7] | the Prop-4.3 sub2 corner set | `paper_src/upstream_facts.json` |
| [Q8] | the G-system indeterminates are the **shifted** stripped `D~_j` | convention; four corroborations in §3.3 |

**INFERRED** — nothing. Where the argument stops it stops; see §9.

**The two load-bearing imports, stated plainly.**

* **[Q6] is the hinge.** Without `t^10 | dm2,dm3,dm4` the `n = 0` rows do not
  factor and there is no forced family to evaluate. `SPINE.md` §8 proves it (hand
  proof + three machine checks + an ablation control) and *itself* asks for
  second-party adjudication, because it upgrades a recorded T2-only status. This
  file inherits that request in full. **If [Q6] falls, this file falls.**
* **[Q8] is a convention, not a theorem proved here.** §3.3's corroborations are
  strong — `full_system_bridge.py` says "D~" in so many words, and there is no
  `d3` variable — but a reader who wants to attack the result should attack here
  or at [Q6], not at the algebra.

---

## 8. Controls and ablation

**Ablation — no TWO of the three conditions suffice.** Both checkers agree, by
different mechanisms (P7.1/P7.3: non-unit saturated ideals; V6.1: explicit
witnesses).

| pair kept | verdict | witness |
|---|---|---|
| `{P_6, P_5}` | **satisfiable** | `res_Y = X^3 * (320X^3+160X^2+29X+2)`; tail degree 3 |
| `{P_6, P_4}` | **satisfiable** | `res_Y = X^4 * (307200X^4+204800X^3+42240X^2+2080X-121)`; tail degree 4 |
| `{P_5, P_4}` | **satisfiable** | `res_Y = X^6 * (921600X^4+678400X^3+185880X^2+22260X+961)`; tail degree 4 |
| all three | **UNIT ideal** | — |

Each tail has positive degree, so it has a root `X0 != 0` in `Qbar`; since
`lc_Y(E(6)) = 7 != 0`, a matching `Y0` exists; and `(alpha, eta, gamma) =
(1, Y0/X0^2, X0)` is an honest witness with `gamma != 0`. So all three conditions
are load-bearing and the obstruction is **exactly** determined, not
over-determined by accident.

**Admissibility.** With `d2, d1, d0` unconstrained the three conditions are
satisfiable (explicitly, `d2 = d1 = d0 = h = 0`). They are not self-contradictory;
the contradiction genuinely comes from the SPINE forcing. (P7.4, V6.3.)

**Cross-corroboration against an independently established kill.** Nothing here
sets `d1 = 0`, so the argument must also empty `a10_b0000_T2` — which
`SPINE.md` §6.6 already proved, by an unrelated T2-only route (`A | t^a`, hence
`A = lambda*t^2`, hence `kbox` at `y = -1`). It does:

```
d1 = 0  =>  delta1 = 2*alpha*(4X+1)/gamma^4 = 0  =>  alpha = 0  or  X = -1/4
alpha = 0   : refuted by the horn (-704/7 != 0)
X = -1/4    : (A) gives 7Y^2 = 1, and (C) collapses to the CONSTANT 1 != 0
              (every Y-term cancels)
```

Two entirely different derivations reaching the same verdict on the same cell is
the strongest control available, because *"no survivors" is exactly the shape a
bug takes*. (P7.5, V6.4–V6.5.)

**Conservatism.** Only the **constant term** of each divisibility is used — `t | ·`
rather than the full `t^2 | ·`, `t^4 | ·`, `t^6 | ·`. The conditions at
`M <= 3` are not used at all. The obstruction has slack it never spends.

---

## 9. Frontier impact, the compiler stage, and what is NOT claimed

### 9.1 Read-only census

`positive_slice.py` P8 reads `phase_d_states_sub2_divfilter.json` and
`phase_d_states_sub2_norl_divfilter.json` — read-only, no writes:

| window | C08/C20 | cell | flag cases | states |
|---|---|---|---:|---:|
| sub2 | ON (`rl`) | `a10_b0000_T1` | 19 | 828 |
| sub2 | OFF (`norl`) | `a10_b0000_T1` | 21 | 850 |

With `stage3_spine` closing the other nine columns, this stage closes the last
one. **Standard sub2 goes from ONE cell to EMPTY, in both censuses.**

### 9.2 The compiler stage — deliberately NOT wired in

The brief asks for a stage after SPINE and a frontier regeneration; it also says
not to modify existing files, and another lane owns the alternate regime.
`frontier_rebuild.py`'s `STAGES` is a hard-coded list covering **both** windows,
so editing it would collide with that lane.

Resolution: the stage is emitted as a **drop-in record**,
`positive_slice_stage.json`, in `frontier_rebuild.STAGES`' exact schema:

```json
{"id": "stage4_positive_slice",
 "title": "Positive-slice obstruction (inverse d3-shift polynomiality), n = 0",
 "checker": "python positive_slice.py --quiet && python positive_slice_verify.py --quiet",
 "dead": {"sub2": ["a10_b0000_T1", "a10_b0000_T2"], "sub1": []},
 "applies_after": "stage3_spine"}
```

**`FRONTIER_REBUILD.md` and `frontier_rebuild.py` were not touched, and the
frontier artifacts were NOT regenerated.** Appending this record to `STAGES` and
re-running `python frontier_rebuild.py` is a one-line change that the owning lane
should make. Until then the frontier artifacts do not reflect this result. That
is a deliberate, flagged omission, not an oversight.

### 9.3 What this does not do

1. **sub1 is untouched, and the argument does not transfer as stated.** The
   degree exactness `deg F + deg Z = 10` is a sub2 coincidence (`SPINE.md` §9.4
   makes the same point about its own hinge); under the sub1 caps it is gone. The
   slice formula itself is regime-independent, but the forced `y = -1` values are
   not.
2. **Nothing is entered into the ledger or the DAG.** No `state_kill_ledger.json`,
   no `proof_dag.json`, no `phase_d_states*.json` was written. §9.1 is a
   read-only census and has not been audited.
3. **The net frontier count is not this file's to quote.** C08/C20 can reopen
   branches elsewhere; this result and that downgrade must not be netted.
4. **This is one arithmetic fact about one specific quartic**, exactly as
   `SPINE.md` §9.2 says of its five kills. `res(p, q) = 561971200 != 0` is a
   property of `q = 2048y^4 - 512y^3 + 320y^2 - 240y + 195`. Nothing here says a
   different `(d,e)` pair would behave the same way.
5. **[Q6] is not re-proved here** (§7). It is the single point of failure.
6. **No Gröbner basis over the G-system, no modular arithmetic, no solver.**
   Nothing external was run, so there are no aborts, timeouts, or exit codes to
   report.

---

## 10. Reproduce

```
cd d2_plane_72_108
python -u positive_slice.py                 # 63/63, full derivation + report
python -u positive_slice.py --quiet         # exit 0 iff every check passes
python -u positive_slice_verify.py          # 79/79, independent re-derivation
python -u positive_slice_verify.py --quiet
```

Both are read-only and pure sympy; together they run in well under a minute.
Neither imports the other. Both abort with exit 1 — before reaching the
obstruction — if the positive control of §4 fails.
