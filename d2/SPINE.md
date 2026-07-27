# SPINE.md — the forcing-divisor spine, and the emptiness of nine of the ten sub2 columns

2026-07-25. Checkers: `spine.py` (47/47, `--quiet`), `spine_verify.py`
(21/21, independent re-derivation, `--quiet`). Read-only on every existing
artifact — this lane wrote only `SPINE.md`, `spine.py`, `spine_verify.py`.

> **HEADLINE.** The five-family reduction **holds**, exactly, for every
> `n = 0..4`, derived from `generators.json` and verified with residual 0. It
> also does far more than reduce: it **closes nine of the ten surviving sub2
> columns**, on both branches, with no Gröbner basis, in exact rational
> arithmetic.
>
> `a9_b1000_T2` — the target of this brief — is **EMPTY**, and it is empty *at
> the top stratum*: the reduced system forces `deg d2 = 4` and `deg d0 = 8` (the
> top stratum's own coordinates), pins the whole spare space down to a
> 3-parameter family on which `G1 = G2 = G3 = 0` hold **identically**, and then
> kills that family by evaluating the boxed identity at `y = -1`.
>
> **The one survivor is `a10_b0000_T1`** (19 flag cases / 828 states — the
> largest column). `d1`, free of degree 6, is precisely what saves it; §9.1 says
> exactly how.

---

## 0. Verdict table

| `n` | family | `a` | T2 (`d1 = 0`) | T1 (`d1 != 0`) | killed by |
|---:|---|---:|---|---|---|
| 0 | `a10_b0000` | 10 | **EMPTY** | **OPEN** ← the residue | §6.6 (T2 only) |
| 1 | `a9_b1000` | 9 | **EMPTY** | **EMPTY** | §6.5 (+§6.6, §7 on T2) |
| 2 | `a8_b1100` | 8 | **EMPTY** | **EMPTY** | §6.5 |
| 3 | `a7_b1110` | 7 | **EMPTY** | **EMPTY** | §6.4 |
| 4 | `a6_b1111` | 6 | **EMPTY** | **EMPTY** | §6.4 |

Read-only census of `phase_d_states_sub2.json` (`spine_verify.py` V7): the ten
columns surviving the divisor filter are

```
a6_b1111_T1  10c/354s     a8_b1100_T1  15c/592s     a9_b1000_T2   4c/ 55s
a7_b1110_T1  11c/481s     a8_b1100_T2   3c/ 39s     a10_b0000_T1 19c/828s
a7_b1110_T2   3c/ 25s     a9_b1000_T1  15c/692s
```

Nine of them (61 flag cases / 2238 states) are closed here. **`a10_b0000_T1`
(19 / 828) is the entire residue.**

Status words, in the sense the handoff asks for:

* **PROVED** — an exact polynomial identity or a finite exhaustive enumeration,
  machine-checked here in char 0, over premises listed in §8.
* **CHECKED** — reproduced from an existing artifact without re-proving it: the
  K-syzygy, `e | 2Phi`, `deg e = 10`, the canonical `G5 = G5body + Phi`, the
  window caps.
* **INFERRED** — nothing in this file. Where an argument stops it stops; see §9.

**This file contradicts one recorded status.** `GSYSTEM_CELL.md` §7.9 and
`DIVISOR_SYZYGY.md` §4 record `t^a | dm2,dm3,dm4` as **T2-only, open on T1**.
§8 proves it on **both** branches, by a strictly stronger argument, with a hand
proof and two independent machine checks plus an ablation control. That upgrade
is what removes all conditionality from the T1 column of the table above, so it
should be adjudicated before the table is quoted.

---

## 1. Where the five families come from

`DIVISOR_SYZYGY.md` proves the universal K-syzygy

```
2*(G5 + d2*G3 + d1*G2 + d0*G1)  ==  2*Phi - e*(d2*e^2 + 3*e*S + 3*R^2)
```

(`e = dm1`, `R = dm2`, `S = dm3`, `T = dm4`), hence `e | 2*Phi` on any G-system
lift; hence with `Phi = c*t^30*q`, `t = y+1`, `c = -1/6630`, `q` the fixed
squarefree irreducible quartic, `q(-1) = 3315 != 0`,

```
e = gamma * t^a * Rm ,     Rm | q  squarefree ,     a + deg Rm = deg e = 10.
```

That filter kills 140 of 220 sub2 flag cases. Its survivors are exactly the
cells it cannot see, indexed by one integer

```
n := deg Rm = 10 - a   in  {0,1,2,3,4}
```

`Q := q/Rm`, `deg Q = 4-n`, `Q(r_i) != 0` at each marked root (`q` squarefree).

---

## 2. The parametrisation, derived

With `t^a | dm2, dm3, dm4` (§8 — **proved here on both branches**),

```
e = gamma * t^a * Rm ,   dm2 = t^a * A ,   dm3 = t^a * B ,   dm4 = t^a * C
```

The certified sub2 stripped window caps (`full_system_bridge.WEIGHT` /
`STRIP_DEGCAP`: cap `2k` at weight `12k`) give

| object | weight | cap | after stripping `t^a` |
|---|---:|---:|---|
| `d2` | 24 | 4 | — |
| `d1` | 36 | 6 | — |
| `d0` | 48 | 8 | — (**eliminated; the cap is never used**) |
| `e = dm1` | 60 | 10 | `deg Rm = n` |
| `dm2` | 72 | 12 | `deg A <= n+2` |
| `dm3` | 84 | 14 | `deg B <= n+4` |
| `dm4` | 96 | 16 | `deg C <= n+6` |
| `Phi` | 204 | 34 | — |

Substituting into the canonical generators (loaded from `generators.json`, the
de-pickled term list — never a hand transcription) the four rows factor
**exactly**:

```
G1 =   3    * t^(2a) * g1
G2 =  (3/2) * t^(2a) * g2
G3 =   3    * t^(2a) * g3
K  = -gamma * t^(3a) * Rm * kbox
```

```
g1   = (1/2)*gamma^2*d1*Rm^2 + gamma*Rm*(d2*A + C) + A*B
g2   = d2*A^2 + 2*A*C + B^2 - gamma^2*d0*Rm^2
g3   = -gamma*d0*Rm*A - (1/2)*d1*A^2 + B*C - (1/6)*gamma^3*t^a*Rm^3
kbox = 3*A^2 + gamma^2*d2*Rm^2 + 3*gamma*Rm*B - (2c/gamma)*t^(3n)*Q
```

Residual **exactly 0** for `n = 0..4` on **both** branches (`spine.py` S2, S5;
`spine_verify.py` V2, V2b — the latter by unassisted `sp.div` with the quotient
not supplied). `t^(3a) = t^(30-3n)` is what cancels `Phi`'s `t^30` and leaves
`t^(3n)` on the boxed row.

---

## 3. The boxed identity

`kbox = 0` is the brief's boxed identity, and it comes out on the nose:

> ```
> 3*A^2 + gamma^2*d2*Rm^2 + 3*gamma*Rm*B  ==  (2c/gamma) * t^(3n) * Q
> ```
> `2c/gamma = -1/(3315*gamma)`.

**Every term has degree exactly `4 + 2n`** — not "at most"; it is perfectly
graded (`spine.py` S3, S4; `spine_verify.py` V3):

| `n` | `2*deg A` | `deg d2 + 2n` | `n + deg B` | `deg t^(3n)Q = 3n + (4-n)` |
|---:|---:|---:|---:|---:|
| 0 | 4 | 4 | 4 | 4 |
| 1 | 6 | 6 | 6 | 6 |
| 2 | 8 | 8 | 8 | 8 |
| 3 | 10 | 10 | 10 | 10 |
| 4 | 12 | 12 | 12 | 12 |

`dm4` does **not** appear. That is the first sign the boxed identity alone
cannot close a cell (§7.1).

### 3b. In the `(A, Sbar)` coordinates

§6.1 forces `Rm | B`; writing `B = Rm*v` and `Sbar := v/gamma` (so
`dm3 = e*Sbar` — this **is** the coordinator's `e | S`, reached here by a
different route, §8.3) the boxed row collapses to two blocks:

> ```
> 3*A^2 + gamma^2*Rm^2*(d2 + 3*Sbar)  ==  (2c/gamma) * t^(3n) * Q
> ```

residual 0 for all `n` (`spine.py` S21). Likewise `g1 = 0` solves for `C`:

```
C = -A*(d2 + Sbar) - (1/2)*gamma*d1*Rm            (i.e. dm4 is not a spare)
```

which is **literally** `T = -R*(S/e + d2) - d1*e/2`, verified identical
(`spine.py` S19).

---

## 4. The spare collapse — confirmed, then superseded

**As the brief asked.** Against the actual ansatz
(`full_system_bridge.STRIP_DEGCAP["sub2"]`, `{dm2:12, dm3:14, dm4:16}` →
`13+15+17 = 45`), the `t^a` reduction alone gives (`spine.py` S6):

| `n` | family | `a` | full | after `t^a` | `45 - 3a` |
|---:|---|---:|---:|---|---:|
| 0 | `a10_b0000` | 10 | 45 | `{3,5,7}` = **15** | 15 |
| 1 | `a9_b1000` | 9 | 45 | `{4,6,8}` = **18** | 18 |
| 2 | `a8_b1100` | 8 | 45 | `{5,7,9}` = **21** | 21 |
| 3 | `a7_b1110` | 7 | 45 | `{6,8,10}` = **24** | 24 |
| 4 | `a6_b1111` | 6 | 45 | `{7,9,11}` = **27** | 27 |

Matches the brief exactly, and reproduces `GSYSTEM_CELL.md` §4.5's `45 → 18`.

**But that count is superseded.** With `dm4` determined (§3b) and `dm3 = e*Sbar`
(so `deg Sbar = 14 - 10 = 4`), the honest free-spare count is (`spine.py` S22):

| `n` | family | `A` | `Sbar` | `dm4` | total | `= n+8` | old `45-3a` |
|---:|---|---:|---:|---:|---:|---:|---:|
| 0 | `a10_b0000` | 3 | 5 | 0 | **8** | 8 | 15 |
| 1 | `a9_b1000` | 4 | 5 | 0 | **9** | 9 | 18 |
| 2 | `a8_b1100` | 5 | 5 | 0 | **10** | 10 | 21 |
| 3 | `a7_b1110` | 6 | 5 | 0 | **11** | 11 | 24 |
| 4 | `a6_b1111` | 7 | 5 | 0 | **12** | 12 | 27 |

**`45 -> n+8`.** For `a9_b1000` that is `45 -> 9`, not `45 -> 18`.

---

## 5. The elimination certificate

The three non-boxed rows collapse to one. Both branches; `spine.py` S7 /
`spine_verify.py` V4, residual exactly 0. With `u := gamma*d2`,
`w := (1/2)*gamma^2*d1*Rm` (`w = 0` on T2), `B = Rm*v`, and

```
F := A*(u + 2*v) + w                 Z := A^2 - gamma*Rm^2*v
```

there is an explicit **polynomial** cofactor identity (`g1r := g1/Rm`):

```
F*Z - (1/6)*gamma^5*t^a*Rm^4  =  (A^2 + Z)*g1r - gamma*A*g2 + gamma^2*Rm*g3
```

so on the variety

> ```
> F * Z  ==  (1/6) * gamma^5 * t^a * Rm^4                                  (*)
> ```

`d0` and `C` are both **eliminated**: the cap `deg d0 <= 8` is never used and no
`sigma` ansatz enters — ansatz-free in the sense `GSYSTEM_CELL.md` §7.5 asks
for. `spine_verify.py` V4b re-confirms (\*) at 40 random exact-rational points
of `{g1r = g2 = g3 = 0}`.

---

## 6. The emptiness argument

Throughout: `gamma != 0` (else `e = 0`, `2*Phi = 0`), `c != 0`,
`mu := 2c/gamma != 0`, `q` squarefree with `q(-1) != 0`, `Rm` a squarefree
divisor of `q` with roots `r_1..r_n != -1`, `Q = q/Rm` with `Q(r_i) != 0`.

### 6.1 The marked roots force `A(r_i) != 0` and `Rm | B`

`kbox = 0` at `y = r_i` (where `Rm` vanishes):

```
3*A(r_i)^2 = mu * (r_i+1)^(3n) * Q(r_i)
```

For `n >= 1` the right side is a product of nonzero factors, so **`A(r_i) != 0`**
and `gcd(A, Rm) = 1`. Then `g1 = 0` at `y = r_i` reads `A(r_i)*B(r_i) = 0`, so
**`B(r_i) = 0`**, i.e. `Rm | B`. Write `B = Rm*v`, `deg v <= 4`. (At `n = 0`,
`Rm = 1`; both statements are vacuous and `v = B`.) Consequently
`Z(r_i) = A(r_i)^2 != 0`, so **`gcd(Z, Rm) = 1`**.

*This step is also the local proof of `e | dm3`* — see §8.3.

### 6.2 Degree exactness — the caps have zero slack

```
deg A <= n+2     deg F <= max(deg A + 4, deg d1 + n) = n+6     deg Z <= 2n+4
```

and `deg( t^a * Rm^4 ) = a + 4n = 10 + 3n`. Now

```
(n+6) + (2n+4) = 3n + 10 = deg RHS
```

**exactly**. The right side of (\*) is nonzero, so both factors are nonzero and
their degrees add: **every cap is attained**, `deg F = n+6`, `deg Z = 2n+4`
(`spine.py` S8). This is where the certified cap `deg d2 <= 4` becomes
load-bearing, and it is tight precisely at the **top** of the `d2`
stratification.

### 6.3 `Z | t^a`

From (\*), `Z | t^a*Rm^4`; by §6.1 `gcd(Z, Rm) = 1`; so **`Z | t^a`**, hence
`Z = zeta*t^(2n+4)`, `zeta != 0`, and `2n+4 = deg Z <= a = 10-n`.

### 6.4 `n <= 2` — kills `n = 3, 4` on **both** branches

| `n` | `deg Z = 2n+4` | `a = 10-n` | |
|---:|---:|---:|---|
| 0 | 4 | 10 | ok |
| 1 | 6 | 9 | ok |
| 2 | 8 | 8 | ok |
| 3 | **10** | 7 | **impossible → EMPTY** |
| 4 | **12** | 6 | **impossible → EMPTY** |

`a7_b1110` and `a6_b1111` die from nothing but a degree count (`spine.py` S10).

### 6.5 Feeding `Z` back into the boxed row — kills `n = 2` and `n = 1`

Substituting `A^2 = zeta*t^(2n+4) + gamma*Rm^2*v` and `B = Rm*v` into
`kbox = 0` gives, exactly (`spine.py` S11),

```
gamma * Rm^2 * (u + 6v)  =  t^(3n) * V_n ,     V_n := mu*Q - 3*zeta*t^(4-n)
```

`gcd(Rm, t) = 1`, so **`Rm^2 | V_n`**, with `deg V_n <= 4-n`.

* **`n = 2`:** `deg V_2 <= 2 < 4 = deg Rm^2`, so `V_2 = 0`; at `y = -1`,
  `mu*Q(-1) = 0`. But `mu != 0` and `Q(-1) != 0` (`q(-1) = 3315`).
  **Contradiction — `a8_b1100` is EMPTY.**
* **`n = 1`:** `deg V_1 <= 3`, `deg Rm^2 = 2`, so `V_1 = Rm^2*(linear)`:

  ```
  mu*Q(r)  - 3*zeta*(r+1)^3 = 0
  mu*Q'(r) - 9*zeta*(r+1)^2 = 0
  ```

  a 2×2 homogeneous system in `(mu, zeta)` with `mu != 0`, so its determinant
  `3*(r+1)^2*((r+1)*Q'(r) - 3*Q(r))` vanishes; `r != -1`, so
  `(r+1)*Q'(r) = 3*Q(r)`. With `q = (y-r)*Q`: `Q(r) = q'(r)`,
  `Q'(r) = q''(r)/2`, hence

  ```
  (r+1)*q''(r) = 6*q'(r)   and   q(r) = 0.
  ```

  But `W3 := (y+1)*q'' - 6*q' = -24576*y^3 + 30720*y^2 - 6272*y + 2080`, and
  **`gcd(q, W3) = 1`** — equivalently
  `res(q, W3) = 28619860707246607140126720 != 0`, so they share no root in
  **any** field extension. **Contradiction — `a9_b1000` is EMPTY.**
* **`n = 0`:** `Rm = 1`, so §6.5 gives nothing; see §6.6.

Nothing in §6.3–6.5 uses `d1 = 0`: these are branch-independent.

### 6.6 The T2-only step — `A | t^a`, and `n = 0`

On T2, `w = 0`, so `F = A*(u+2v)` and (\*) shows **`A` divides `t^a*Rm^4`**;
with `gcd(A, Rm) = 1` and `deg A = n+2` attained (§6.2),

```
A = lambda * t^(n+2) ,   lambda != 0.
```

**Corollary (`spine.py` S24):** `dm2 = t^a*A = lambda*(y+1)^(a+n+2) =
lambda*(y+1)^12` for **every** `n`. The flagged inference `R = c*(y+1)^rho` on
T2 is therefore a **theorem**, with `rho = 12` exactly, and `deg R` is at its
cap.

Cancelling `A` from (\*): `(u+2v)*Z = (gamma^5/(6*lambda))*t^(8-2n)*Rm^4`, so
`Z | t^(8-2n)` and `2n+4 <= 8-2n`, i.e. **`n <= 1`** (`spine.py` S15) — a
second, independent kill of `n = 2,3,4` on T2.

At **`n = 0`** (`Rm = 1`, `a = 10`): `A = lambda*t^2`, `Z = zeta*t^4` forces
`B = ((lambda^2-zeta)/gamma)*t^4`, and `u + 2B = (gamma^5/(6*lambda*zeta))*t^4`
forces `u = tau*t^4`. Every term of `kbox = 3*A^2 + gamma*u + 3*gamma*B - mu*q`
except `mu*q` is then divisible by `t^4`, so at `y = -1`: `mu*q(-1) = 0`.
**Contradiction — `a10_b0000_T2` is EMPTY** (`spine.py` S16).

---

## 7. The top stratum of `a9_b1000` — what the reduction forces there

### 7.1 The boxed identity **alone** does not bite (measured, not guessed)

At `n = 1`, top stratum (`deg d2 = 4`, `deg sigma = 8`), `kbox = 0` is 7 scalar
equations in degree 0..6. But `gamma^2*d2*Rm^2` carries `d2`'s five free
coefficients, so one can simply **solve for `d2`**:

```
d2 = [ mu*t^3*Q - 3*A^2 - 3*gamma*Rm*B ] / (gamma^2 * Rm^2)
```

whenever the numerator is divisible by `(y-r)^2`; the derivative condition then
solves for `B(r)`. The entire content of the boxed row at the top stratum is the
**single** scalar condition `3*A(r)^2 = mu*(r+1)^3*q'(r)`, satisfiable over `C`
for any `A(r) != 0`.

**So the divisor syzygy by itself does not close the cell at the top stratum.**
It is the reduction, not the mechanism. Saying otherwise would repeat exactly
the overreach the handoff warns about.

### 7.2 What does bite: the elimination certificate, and it bites at the top

At `n = 1`, `a = 9`: `deg A = 3`, `deg F = 7`, `deg Z = 6`,
`deg(t^9*Rm^4) = 13 = 7 + 6`. **Zero slack.** Running §6.6 then §6.3 at `n = 1`
on T2 forces, with no freedom left:

| object | forced value | degree | cap |
|---|---|---:|---:|
| `A` | `lambda * t^3` — so `dm2 = lambda*(y+1)^12` | 3 | 3 |
| `v` (`B = Rm*v`) | `0` — so **`dm3 = 0` identically**, `Sbar = 0` | — | — |
| `u = gamma*d2` | `kappa*Rm^4`, `kappa = gamma^5/(6*lambda^3)` | 4 | 4 |
| `d2` | `(kappa/gamma)*(y-r)^4` | **4** | 4 |
| `C` | `-(lambda*kappa/gamma)*t^3*Rm^4` | 7 | 7 |
| `d0` | `-(lambda^2*kappa/gamma^3)*t^6*Rm^2` | **8** | 8 |

Two things to read off.

1. **`deg d2 = 4` and `deg d0 = 8` are forced.** The surviving family lives
   *only* at the **top stratum**; every lower `d2`-stratum of `a9_b1000_T2` is
   already empty at the step "`v = 0`, hence `deg u = 4`". This is the shape
   `SESSION_HANDOFF.md`'s SPEC demands and the opposite of Φ-depth's bottom-up
   behaviour.
2. **The family is genuine for three of the four rows.** Substituting back gives
   `g1 = 0`, `g2 = 0` identically and
   `g3 = -Rm^3*t^9*(gamma^5 - 6*kappa*lambda^3)/(6*gamma^2)`, which vanishes
   exactly at `kappa = gamma^5/(6*lambda^3)`. So `<G1,G2,G3>` is **not** empty
   here — it is a 3-parameter family `(gamma, lambda, r)`. The cell is closed by
   the **fourth** row alone:

   ```
   kbox |_{y=-1}  =  gamma^6 * Rm(-1)^6 / (6*lambda^3)  !=  0
   ```

   `Rm(-1) = -1-r != 0`, `lambda != 0`, `gamma != 0`.
   **`a9_b1000_T2` is EMPTY, at the top stratum and everywhere below it.**
   (`spine.py` S18; re-checked concretely over `Q[r]/(q)` — `G1 = G2 = G3 = 0`
   exactly, `G5 = K/2 != 0`.)

The 9 free spares collapse to `{lambda}`: `dm2 = lambda*(y+1)^12`, `dm3 = 0`,
`dm4 = -(lambda*kappa/gamma)*(y+1)^12*(y-r)^4`. The `deg sigma = 8` half of the
pilot cell — which cost the G-system engine `> 3600 s` at one extra ring
variable and was never reached at `deg d2 = 4` — is settled here by evaluating
one polynomial at `y = -1`.

### 7.3 On T1

`a9_b1000_T1` dies too, but by §6.5 rather than §6.6: `Z = zeta*t^6` and the
boxed row give `Rm^2 | (mu*Q - 3*zeta*t^3)`, hence `(r+1)*q''(r) = 6*q'(r)`,
hence the same `gcd(q, W3) = 1` contradiction. §6.5 never uses `d1 = 0`, so the
T1 kill needs nothing beyond §8.

---

## 8. `t^a | dm2, dm3, dm4` — proved on **both** branches

This is the one load-bearing import, and it is the one place where this file
**upgrades a recorded status**. The record (`GSYSTEM_CELL.md` §7.9,
`DIVISOR_SYZYGY.md` §4) is: proved on T2 (`d1 = 0`) for `a in {7,8,9}` by a
two-relation valuation enumeration (the identity plus the `dm4`-eliminated
`H3`); **open on T1**. That enumeration is `tpower_divisibility.py`.

### 8.1 The hand proof

Let `rho = v_t(R)`, `s = v_t(S)`, `tau = v_t(T)` at `t = y+1`, with
`v_t(e) = a <= 10` and `v_t(2*Phi) = 30`. Valuations add on products, and a sum
of terms whose minimal valuation is attained **once** cannot vanish.

**Step 1 — `rho >= a`.** Suppose `rho < a`. The row `K = 0` has terms of
valuation `30`, `v(d2)+3a`, `2a+s`, `a+2*rho`. Since `rho <= a-1`,
`a+2*rho <= 3a-2 < 3a <= v(d2)+3a`, and `3a-2 <= 28 < 30`; so the minimum is
`min(2a+s, a+2*rho) < 3a` and must be attained twice, forcing

```
s = 2*rho - a       (and s >= 0, so rho >= a/2).
```

Now `G1 = 0` has terms `v(d1)+2a`, `v(d2)+a+rho`, `a+tau`, `rho+s = 3*rho-a`.
Because `rho < a`, `3*rho-a < a+rho <= v(d2)+a+rho` and `3*rho-a < 2a <=
v(d1)+2a`. So the tie must be with `a+tau`, forcing

```
tau = 3*rho - 2a    (and tau >= 0, so rho >= 2a/3).
```

Finally `G3 = 0` has terms `v(d0)+a+rho`, `v(d1)+2*rho`, `3a`, and
`s+tau = 5*rho-3a`. Because `rho < a`,

```
5*rho-3a < 3a ,   5*rho-3a < a+rho ,   5*rho-3a < 2*rho ,   5*rho-3a >= 0.
```

So the `3*S*T` term is a **strict unique minimum** and `G3 != 0`.
**Contradiction.**

**Step 2 — `s >= a`.** Given `rho >= a` and `s < a`: in `K = 0`,
`2a+s <= 3a-1 < 3a <= min(v(d2)+3a, a+2*rho)` and `3a-1 <= 29 < 30`, so `2a+s`
is a unique minimum. **Contradiction.**

**Step 3 — `tau >= a`.** Given `rho, s >= a` and `tau < a`: in `G1 = 0`,
`a+tau <= 2a-1 < 2a <= min(v(d1)+2a, v(d2)+a+rho, rho+s)`, a unique minimum.
**Contradiction.**

`d1` appears only inside terms that are bounded **below**, so **no step uses
`d1 = 0`** — the proof is branch-independent. It uses no degree caps at all,
only `v_t(e) = a <= 10` and `v_t(2*Phi) = 30`. (The `a <= 10` is needed twice,
to keep `3a-2` and `3a-1` below 30; it is exactly `deg e = 10`.)

### 8.2 Three machine checks and a control

* `spine.py` S9 — the **repo's** two-relation enumeration, re-implemented from
  scratch, extended to `a = 6..10` on T2. Clean.
* `spine.py` S23 — the four-row enumeration `G1,G2,G3,K` with valuation
  ceilings from the certified caps, **both branches**, `a = 6..10`: zero
  counterexample configurations.
* `spine_verify.py` V5 — the same, but with the term-valuation vectors
  **extracted automatically from `generators.json`** rather than hand-written,
  so a transcription slip cannot be shared between the two files. Same result.
* `spine_verify.py` V5b — the §8.1 hand proof re-checked step by step as
  arithmetic assertions, all `a = 6..10`, all `rho < a`.

**Controls, because "no survivors" is exactly the shape a bug takes.**
(i) *Admissibility*: for every `a` and both branches the scan finds hundreds of
**non**-counterexample configurations that pass all four rows (e.g.
`rho=s=tau=a`, `v(d0)=v(d2)=0`), so the row tests are not vacuously restrictive.
(ii) *Ablation*: dropping `G1`, `G3` or `K` makes counterexamples reappear at
`a = 9`, T1 — each is load-bearing. Dropping `G2` does **not**: `G2` is
redundant, which matches §8.1 using only `K`, `G1`, `G3`.

Why the repo's version stops at T2 and this one does not: the record's
enumeration uses `H3 = dm1*G3|_{dm4}` , whose extra `d1*R^2` and `d1*e^2*S`
terms wreck the case analysis when `d1 != 0`. Keeping `dm4` and using `G1` and
`G3` **separately** avoids that entirely — `d1` then appears only in terms that
are bounded below and can never be the unique minimum.

**Adjudication note.** This is a strictly stronger claim than the record, and it
is the hinge of the whole T1 column of §0. It should be checked by a second
party before the T1 verdicts are entered anywhere. Everything it rests on
(`v_t(e) = a`, `v_t(2*Phi) = 30`, `a <= 10`, the canonical generators) is
already established elsewhere in the repo.

### 8.3 Cross-corroboration with the mid-lane results

Two results arrived from the coordinator while this lane was running. Both are
**already inside this reduction**, reached by a different mechanism — verified
identical, residual 0:

| coordinator's form | this file's form | check |
|---|---|---|
| `e \| S` (Sylvester-resultant integral dependence) | `Rm \| B` from `kbox` + `g1` at the marked roots (§6.1), giving `dm3 = e*(v/gamma)` | `spine.py` S20 |
| `T = -R*(S/e + d2) - d1*e/2` annihilates `G1` | `C = -(A*(u+v) + w)/gamma`, the `g1`-elimination of §5 | `spine.py` S19, residual 0 |

So `Sbar = v/gamma`, and the reduction was already carrying both. What did need
correcting is the **spare count** (§4): `45 -> n+8`, not `45 - 3a`.

`e*R | Phi` and `e | d1*R^2` were **not used** and are not adjudicated here. The
flagged T2 inference `R = c*(y+1)^rho` **is** adjudicated: §6.6 proves it, with
`rho = 12`.

---

## 9. What the five-family picture does **not** resolve

1. **`a10_b0000_T1`.** `n = 0` makes `Rm = 1`, so §6.1 and §6.3–6.5 are all
   vacuous, and the T2 rescue in §6.6 needs `d1 = 0`. The certificate survives
   as `F*Z = (1/6)*gamma^5*t^10` with `Z = zeta*t^4` and
   `F = A*(u+2B) + (1/2)*gamma^2*d1` forced to `phi*t^6` — and `d1`, free of
   degree 6, has exactly the seven coefficients needed to absorb that. The
   reduced boxed row `3*zeta*t^4 + gamma*(u + 6B) = mu*q` likewise has five free
   coefficients against five equations. **`d1` is what saves `n = 0` on T1**,
   which is why §6.6 has no T1 analogue. 19 flag cases / 828 states, the largest
   column, and now the entire residue of sub2.
2. **The kills are five separate arithmetic facts, not one uniform statement.**
   `n = 3,4` die on a degree count; `n = 2` on `q(-1) != 0`; `n = 1` on
   `gcd(q, (y+1)q''-6q') = 1`; `n = 0` (T2) on `q(-1) != 0` again. Each is a
   different property of the specific quartic. Nothing here says a *different*
   `(d,e)` pair would behave the same way.
3. **The families are low-dimensional, not points, before the kills.** Even
   after `A`, `Sbar`, `C`, `d0` are pinned, `(gamma, lambda, kappa, r)` and — at
   `n >= 2` — the coefficients of `v` and `u` remain.
4. **`sub1` is untouched, and the argument does not transfer.**
   `DIVISOR_SYZYGY.md` §4 already records that sub1's caps make the degree
   forcing vacuous. Worse, §6.2's zero-slack count `(n+6)+(2n+4) = 3n+10` is a
   **sub2 coincidence**: under the sub1 caps `{dm2:18, dm3:21, dm4:24}` the cap
   sum exceeds the RHS degree and the exactness — the hinge of everything — is
   gone.
5. **Nothing here is entered into the ledger or the DAG.** No
   `state_kill_ledger.json`, no `proof_dag.json`, no `phase_d_states*.json` was
   written. The §0 impact figures are a read-only census
   (`spine_verify.py` V7), not a ledger entry, and have not been audited.
6. **The C08/C20 field-scope objection is untouched.** It does not threaten this
   argument — every step is a polynomial identity or a gcd/resultant over `Q`,
   and a gcd of 1 over `Q` means no common root in any extension, so the
   conclusions are genuinely about `C`-points. But it can reopen branches
   elsewhere, so the net frontier count is not this file's to quote.
7. **No Gröbner basis, no modular arithmetic, no solver.** Nothing here is
   mod-`p` reconnaissance; there are no aborts or timeouts to report because
   nothing external was run.

---

## 10. Reproduce

```
cd d2_plane_72_108
python -u spine.py              # 47/47, full derivation and report
python -u spine.py --quiet
python -u spine_verify.py       # independent re-derivation
python -u spine_verify.py --quiet
```

Both are read-only and pure sympy; no Singular, no msolve, no WSL.
`spine_verify.py` takes several minutes (symbolic-coefficient polynomials of
degree ~30 at `n = 4`).
