# CONTACT_LEMMA.md — the general contact identity, the cascade as one lemma, and its reach

2026-07-26. Checker: `contact_lemma.py` (`--quiet`, exit 0, **64 checks**, ~90 s,
pure sympy).

This lane wrote only `CONTACT_LEMMA.md` and `contact_lemma.py`. It **modified no
existing file**: not `PROOF_72_108.md`, not `SLICE_OBSTRUCTION.md`, not
`slice_obstruction_basis.py` / `slice_obstruction_audit.py`, not `run_tests.sh`,
not `tools/`, not any ledger / DAG / frontier artefact. No Singular, no msolve,
no WSL, no subprocess, no solver — so there are no exit codes to misread.

> **HEADLINE.**
>
> 1. **The identity is true, and in closed form.** For all `m, n`
>    ```
>    m*H^n - n*H^m + (n-m)  =  sum_{j>=2} c_j*(H-1)^j ,
>        c_j = m*binom(n,j) - n*binom(m,j),      c_2 = m*n*(n-m)/2 .
>    ```
>    The contact order at `H=1` is **exactly 2** for every pair except
>    `(m,n) = (1,2)` (where `G = (H-1)^2` on the nose). At `(2,3)`: `c_2 = 3`,
>    `c_3 = 2`, i.e. `2H^3-3H^2+1 = 3K^2+2K^3 = (H-1)^2(2H+1)`. **PROVED
>    symbolically in `(m,n)`**, plus all 91 pairs `1<=m<n<=14` exactly, plus 25
>    mutants that all lose the double root.
>
> 2. **The general bound is `v_t(h_k) >= m*k - 1`, NOT `2k-1`.** `2k-1` is the
>    `m = 2` specialisation. At `(3,5)` the correct profile is `3k-1`. The
>    exponents themselves generalise to
>    ```
>    t^(m*k - r)  |  [u^k] H^r          for r = m (the P family) and r = n (the Q family),
>    ```
>    so `P_k = m*k-m`, `Q_k = m*k-n` — **intercept/power = 1 for both**, and that
>    common `1` is the whole mechanism. The general cokernel is `m*k-n` (`2k-3`
>    at `(2,3)`).
>
> 3. **The cascade is two lemmas, not one.** LEMMA A (the fixed point) is
>    general and easy. The *forcing* splits: the audited "one clean perfect
>    square per level" engine completes the profile **iff `(m,n) = (2,3)`**
>    (LEMMA C + gate E9). For every other pair the last step per index
>    degenerates, and what closes it is **LEMMA B, a leading-order rigidity
>    statement gated on `gcd(m,n) = 1`**.
>
> 4. **`gcd(m,n) = 1` is the arithmetic gate.** Exactly the shape
>    `MINIMAL_CORE.md` said to look for. For every non-coprime pair the rigidity
>    conclusion is **false**, with an explicit counterexample `Psi = U^(1/g)`.
>    Both live corners are coprime: `gcd(2,3) = gcd(3,5) = 1`.
>
> 5. **The answer to task 4: YES, the slice obstruction does apply at `(3,5)`,**
>    even though `passport_75_125.py`'s Belyi gate `u*kappa = m+n-1` fails there.
>    They are different mechanisms and they do not share a gate. But the
>    `(2,3)` *engine* does **not** transfer verbatim: at `(3,5)` index 1 needs
>    **two** steps (a square then a **cube**), and index 2's last step **stalls**
>    for the clean-jet engine and needs LEMMA B.
>
> 6. **A new sharp negative.** The profile is **not unconditional**. An explicit
>    witness `H = (1 + a*t^(m(d-1))*u^d)^(1/m)` satisfies `(P)`, `(P0)` and
>    *every* `(Q)` condition in range and yet has `v_t(h_d) = m*d-m < m*d-1`.
>    At `(2,3)` it lands at `d = 8` with `v_t(h_8) = 14 < 15`, **inside the
>    degree caps in both windows**. So `slice_obstruction_audit.py` G9's
>    "this cascade can advance at most `h_1..h_7`" is **SHARP, not
>    conservative** — a fact the audit recorded as a scope note and could not
>    prove. `a_t >= 9` uses `k = 5`, two indices inside the ceiling.

---

## 0. What was asked, and what survived contact

| claim under test | verdict |
|---|---|
| `m*H^n - n*H^m + (n-m)` has vanishing constant and linear terms, quadratic coefficient `m*n*(n-m)/2` | **PROVED**, symbolically in `(m,n)` (§1) |
| "because the two power maps share a tangent direction after weighting" | **CONFIRMED, and that is literally the proof**: `G'(1) = m*n - n*m` (§1.1) |
| the cascade is one filtered-power-series lemma | **PARTLY.** It is *three* statements with different scopes (§2). A single lemma covering all `(m,n)` would be false — see §4.3 |
| the bound is `v_t(h_m) >= 2m-1` in general | **NO.** It is `v_t(h_k) >= m*k-1`. `2k-1` is the `m=2` case (§2.1) |
| the cokernel is `2n-3` in general | **NO.** It is `m*k-n`; `2k-3` at `(2,3)` (§2.4) |
| "reusable engine: any corner with `P = C^m`, `Q = C^n + ...` admits the same obstruction" | **QUALIFIED YES**, under three gates: `gcd(m,n)=1`, `lam >= m`, `N_Q >= D_P+D_Q` (§3) |
| does the obstruction apply at `(3,5)`? | **YES** (§4), by LEMMA B; the clean-jet engine alone does **not** get there |
| the `(3/4)` in the audited table is a second constant | **NO.** It is `c_2 = 3` with the square's `1/2` cleared: `3*(g-A^2/2)^2 = (3/4)*(A^2-2g)^2` (§5, E3) |

Nothing in the external review's algebra was found wrong. Two corrections of
substance: the bound and the cokernel are `m`-dependent, and the *forcing*
mechanism is `(2,3)`-special even though the *identity* is universal.

---

## 1. The contact identity

### 1.1 Statement and proof

Let `G(H) := m*H^n - n*H^m + (n-m)` over any commutative ring, `m` and `n` any
integers (no coprimality, no ordering needed here).

**Proposition 1.** `G(1) = 0`, `G'(1) = 0`, and for every `j >= 1`
```
G^(j)(1)/j!  =  m*binom(n,j) - n*binom(m,j) .
```
Hence `(H-1)^2 | G`, with quadratic coefficient `c_2 = m*n*(n-m)/2`.

*Proof.* `G(1) = m - n + (n-m) = 0`. Differentiating,
`G'(H) = m*n*H^(n-1) - n*m*H^(m-1)`, so `G'(1) = m*n - n*m = 0` — **this is the
"shared tangent direction": the weights `m` and `n` are chosen precisely so the
two first derivatives `n*H^(n-1)` and `m*H^(m-1)`, weighted by `m` and `n`,
become the same number `m*n` at `H=1`.** For `j >= 1`,
`d^j/dH^j (m*H^n) = m*(n)_j*H^(n-j)` with `(n)_j` the falling factorial, so
`G^(j)(1)/j! = m*(n)_j/j! - n*(m)_j/j! = m*binom(n,j) - n*binom(m,j)`. At `j=1`
that is `m*n - n*m = 0`; at `j=2` it is `m*n(n-1)/2 - n*m(m-1)/2 = m*n*(n-m)/2`.
At `j=0` the closed form returns `m-n`, and the additive constant `(n-m)` cancels
it — **that is the only role the constant plays.** □

**Corollary (contact order).** For `m < n`: `c_j > 0` for all `2 <= j <= n`,
`c_n = m`, and `c_j = 0` for `j > n`. In particular the contact order is
**exactly 2** unless `c_3 = (m*n/6)[(n-1)(n-2) - (m-1)(m-2)] = 0`, which for
`m<n` happens iff `m+n = 3`, i.e. only at `(m,n) = (1,2)`, where
`G = (H-1)^2` identically.

**At `(m,n) = (2,3)`:** `c_2 = 2*3*1/2 = 3`, `c_3 = 2`, so
`2H^3 - 3H^2 + 1 = 3K^2 + 2K^3 = (H-1)^2*(2H+1)` — exactly
`SLICE_OBSTRUCTION.md` §2. **The `3` in `3*([t^(2m-2)]h_m)^2` is `c_2`.**

### 1.2 How it is checked, and why the check is not vacuous

`contact_lemma.py` A1–A14:

* **symbolic in `(m,n)`** (A1–A6): `sympy` differentiates `m*H^n - n*H^m` with
  `m,n` free symbols; `G(1)`, `G'(1)`, `G''(1)/2`, `G'''(1)/6` and the closed
  form for `j = 1..8` are verified as identities in `(m,n)`, not sampled;
* **exact, all 91 pairs** `1 <= m < n <= 14` (A7–A9): the full `K`-expansion,
  `c_2`, and the multiplicity-exactly-2 statement with the `(1,2)` exception;
* **MUTATION (A11), 25 mutants.** *This is the check that matters*, because
  "some polynomial vanishes to order 2 at 1" would be a weak claim if it held
  for everything. It does not. For each of `(2,3),(3,5),(2,5),(3,4),(4,7)`:
  the unweighted `H^n - H^m`, the mis-weighted `(m+1)H^n - nH^m + ...`,
  `mH^n - (n+1)H^m + ...`, the weight-swapped `nH^n - mH^m + ...`, and the
  wrong-constant `... + (n-m)+1` **all lose the double root**. The cancellation
  is a property of the `m/n` weighting, not of differences of powers.
* **A12** kills the plausible-but-wrong `c_2 = m*n*(n-m)` (no `/2`): it agrees
  with the truth on **no** pair, so A8 is a genuine equality test.
* **A13**: `m = n` degenerates `G` to `0` identically — `m != n` is load-bearing.

This is the specific trap the brief warned about ("a Riemann–Hurwitz test that
was an identity for all inputs"). §1 asserts an identity **and** exhibits inputs
for which the corresponding statement is false.

---

## 2. The cascade, as lemmas with explicit hypotheses

### 2.0 Setting, and where the exponents come from (§B of the checker)

`C_SERIES_75_125.md` §4 fixes the stripping convention at every corner:
`d_k := c_k * c^(m*(ell-k) - 1)` with `c := c_ell` the leading coefficient, `m`
the P-power (`a`), `ell` the chart parameter (`t = l`). Setting
`H(u) := sum_{j>=0} d_{ell-j} u^j` (so `h_0 = 1`) and `u := 1/(c^m x)` gives
`C = c*x^ell*H(u)` and, **by a one-line convolution** (verified symbolically,
B1):
```
(C^r)_M  =  c^( r - m*k ) * [u^k] H^r ,        k := r*ell - M .
```
With `v_t(c) = 1` (true at both corners: `C4 = y^7(y+1)` at `(72,108)`,
`C = y^2(y^3+1)` at `(75,125)`, both with a simple zero at `y=-1`),
polynomiality of the slice is
```
       ***   t^(m*k - r)  |  [u^k] H^r   ***
```
so, writing `p_k := [u^k]H^m`, `r_k := [u^k]H^n`:
```
(P)   t^(m*k - m) | p_k     k = 1..N_P ,      N_P = m*ell
(P0)  p_k = 0               k > N_P                            [P has no x^(<0)]
(Q)   t^(m*k - n) | r_k     k = 1..N_Q ,      N_Q = (n+1)*ell - 1
```
At `(m,n,ell) = (2,3,4)` this is `t^(2k-2) | p_k` for `k<=8`, `p_k = 0` for
`k>=9`, `t^(2k-3) | r_k` for `k<=15` — **`SLICE_OBSTRUCTION.md` §1 and
`slice_obstruction_audit.py` G9 verbatim** (B2, B6). At `(3,5,5)`:
`t^(3k-3) | p_k` for `k<=15`, `t^(3k-5) | r_k` for `k<=29` (B3, B6).

**The key structural fact (B4).** Both intercepts are the power itself:
`P_k = m*k - m*1`, `Q_k = m*k - n*1`. The ratio intercept/power is `1` for both
families, automatically, for every `(m,n)`. §2.1 is that observation.

### 2.1 LEMMA A — the fixed point (general, PROVED)

**LEMMA A.** Let `T` be a valuation ring of characteristic `0`, `v = v_t`. If
`v_t(h_k) >= m*k - 1` for all `k`, then for every power `r >= 1` and every `k`
```
v_t([u^k] H^r)  >=  m*k - r ,     with EQUALITY for all k >= r.
```
In particular `(P)` and `(Q)` hold **identically**, and both are **sharp**.

*Proof.* Substitute `u = v/t^m` and write `h_k = t^(m*k-1)*A_k`. Then
`K := H-1 = sum_k t^(m*k-1)A_k v^k/t^(m*k) = Ahat(v)/t` is homogeneous of
`t`-weight `-1`, so
`[u^k]H^r = t^(m*k) * sum_{j=0}^{r} binom(r,j) t^(-j) [v^k]Ahat^j`, whose lowest
term is `j = r`, giving `m*k - r`; it is nonzero exactly when
`[v^k]Ahat^r != 0`, i.e. for `k >= r` (`Ahat` starts at `v^1`). □

Three consequences, all checked (C1–C4):

* **the cascade constrains, it never contradicts** — this is
  `SLICE_OBSTRUCTION.md` §3.2 generalised, and it is why the slice system is
  never emptied by the cascade;
* **`m*k-1` is the unique scaling-compatible profile.** MUTATION C3: tighten the
  Q intercept by one (`n -> n-1`) and the profile fails the family at *every*
  level `k >= n`. The intercept pair `(m,n)` is forced, not merely allowed;
* **the profile is attained** (C4): an explicit rational instance with
  `v_t(h_k) = m*k-1` exactly, at `(2,3)`, `(3,5)` and `(3,4)`.

**The general bound is therefore `m*k - 1`.** At `(2,3)`: `2k-1`. At `(3,5)`:
`3k-1`. It is *not* `2k-1` in general, and it is not `2m-1` — the brief's `2m-1`
is the `(2,3)` value with the index renamed.

### 2.2 The elimination (general, PROVED)

`h_0 = 1` gives `p_k = m*h_k + q_k^P` and `r_k = n*h_k + q_k^Q`, so the fresh
coefficient cancels in `m*r_k - n*p_k` because `m*n - n*m = 0` — which is
exactly `G'(1) = 0` read off the `u`-grading. And (D2)
```
m*r_k - n*p_k  =  [u^k] ( sum_{j>=2} c_j K^j )        exactly,
```
with the `c_j` of §1. **The contact identity *is* the stacked family.** Verified
for `k = 2..10` at `(2,3),(3,5),(2,5),(3,4),(4,7)`. Since `P_k > Q_k` whenever
`m<n` (B5), the reduction is lossless: given `(P)`, `t^(Q_k) | r_k` is
*equivalent* to `t^(Q_k) | [u^k](sum c_j K^j)`.

### 2.3 The cokernel (general, PROVED) and the cap gate

**Proposition 2.** Let `h_k` be a polynomial of `y`-degree `<= lam*k`
(`lam*k+1` free coefficients). If `lam >= m` then, by exact rank over `Q`,
```
coker(P only) = 0 ,    coker(Q only) = 0 ,    coker(STACKED) = Q_k = m*k - n .
```
At `(2,3)`: `0 / 0 / 2k-3` for `k = 2..8` in **both** windows (`lam = 3` sub1,
`lam = 2` sub2) — reproducing `slice_obstruction_audit.py` E4 and
`SLICE_OBSTRUCTION.md` §2 **on the nose** (D3, printed row by row).

**GATE (H-cap): `lam >= m`.** MUTATION D5: with `lam = m-1` the P-only cokernel
becomes nonzero from level `4` (`(2,3)`), `5` (`(3,5)`), `6` (`(4,7)`) on, and
the entire "each side alone absorbs everything" picture collapses. `(2,3)`
satisfies the gate in both windows, **sub2 exactly at equality (`lam = 2 = m`)**
— a hypothesis the audited case meets invisibly, and the first of the three
hidden hypotheses the brief asked for.

### 2.4 LEMMA C — the per-index clean-jet cascade

This is the audited engine, restated generally. State at index `r`: the profile
`v_t(h_i) >= m*i-1` established for `i < r`, current bound `v_t(h_r) >= V`.
Parametrise
```
h_i = t^(m*i-1)*A_i            i < r
h_r = t^V * X                                        (X free)
h_i = -q_i/m + t^(m*i-m)*g_i   r < i <= N_P          (g_i free)
h_i = -q_i/m                   i > N_P               [(P0)]
```
The level-`r` P condition is **automatic** here (`v_t(q_r) >= m*r-2 >= m*r-m` for
`m>=2`), so `h_r` ranges over every series with `v_t >= V`; all other P
conditions hold identically; the residue is `t^(m*L-n) | E_L` with
`E_L := [u^L](sum_j c_j K^j)`.

**Availability (E8).** The pure term `c_j*X_V^j` sits at level `j*r` and
`t`-order `j*V`; it is a *required* jet iff `j*V < Q_{j*r}`, i.e.
```
V  <  m*r - n/j .
```
Taking `j = n` this holds at every `V <= m*r-2`, so a candidate condition
reaching the full profile exists for **every** `(m,n)`. **Availability is never
the obstacle. Cleanliness is.**

**A step is CLEAN** when some required-nonzero jet equals `(unit)*X_V^p`: then
`X_V = 0` is forced, with no factorisation into coprime pieces, hence **no case
split and no silently chosen branch**.

**Gate E9 — where `(2,3)` is special.** The square step (`j=2`) alone reaches
only `V = m*r - floor(n/2)`. For `m >= 2` that equals the target `m*r-1` **iff
`n = 3`, i.e. iff `(m,n) = (2,3)`.** For every other pair the last step per
index needs a strictly higher power — and higher powers collide (§4.2).

### 2.5 LEMMA B — leading-order rigidity, and the `gcd` gate (general, PROVED)

This is the general replacement for the degenerate step, and it is the real
content of this file.

**LEMMA B.** Let `char = 0`, `m < n`, and suppose the current knowledge is a
**linear** weight `v_t(h_i) >= q'*i` for every `i >= 1`. Put `g := m - q' >= 1`,
`B_i := [t^(q'*i)]h_i`, and `Psi(v) := 1 + sum_{i>=1} B_i v^i`. Assume

* **(B1) `gcd(m,n) = 1`;**
* **(B2) `N_Q >= D_P + D_Q`**, where `D_P := floor(m/g)`, `D_Q := floor(n/g)`.

Then `Psi` is a **polynomial** of degree `<= min(floor(D_P/m), floor(D_Q/n))`.
In particular at the critical slope `q' = m-1` (`g = 1`, so `D_P = m`,
`D_Q = n`) the degree bound is `1`:
```
   ***  Psi = 1 + B_1*v ,  i.e.  B_i = 0  and  v_t(h_i) >= (m-1)*i + 1
        for every i >= 2.  ***
```

*Proof.* Because the weight is linear, `[t^(q'*L)]([u^L]H^r) = [v^L]Psi^r`.
`(P)` forces that to vanish whenever `q'*L < m*L - m`, i.e. `L > m/g`; so
`R := Psi^m` is a **polynomial of degree `<= D_P`**. `(Q)` forces
`[v^L]Psi^n = 0` for `n/g < L <= N_Q`.

*Step 1 (propagation).* `y := Psi^n = R^(n/m)` satisfies `m*R*y' = n*R'*y`, i.e.
with `R = 1 + sum_{j=1}^{D_P} r_j v^j`,
```
m*(k+1)*y_{k+1}  =  sum_{j=1}^{D_P} r_j * ( n*j - m*(k+1-j) ) * y_{k+1-j} .
```
This is a `(D_P+1)`-term recurrence with the `y_{k+1}` coefficient `m*(k+1)`
invertible in characteristic `0`. So **`D_P` consecutive vanishing coefficients
propagate forever.** `(Q)` supplies the run `L = D_Q+1, ..., D_Q+D_P` precisely
when `N_Q >= D_P + D_Q`, hence `y = Psi^n` is a polynomial of degree `<= D_Q`.

*Step 2 (rigidity).* `R^n = (Psi^m)^n = (Psi^n)^m = y^m` with `R, y` polynomials.
For every irreducible `p`, `n*val_p(R) = m*val_p(y)`, so `m | n*val_p(R)`, and
**by `gcd(m,n)=1`**, `m | val_p(R)`. Hence `R = U^m` with `U` a polynomial,
normalised by `U(0)=1`; then `y^m = U^(m*n)` and `y(0)=1` give `y = U^n`, and
`Psi^m = U^m` with `Psi(0)=U(0)=1` gives `Psi = U`. Degrees:
`m*deg U <= D_P` and `n*deg U <= D_Q`. □

**GATE (B1) `gcd(m,n) = 1` is load-bearing — MUTATION F7.** If `g := gcd(m,n) > 1`,
take any polynomial `U` of degree `g` and `Psi := U^(1/g)`. Then
`Psi^m = U^(m/g)` and `Psi^n = U^(n/g)` are **polynomials** of degree `<= D_P`
and `<= D_Q`, so both families are satisfied at the critical slope — **and
`B_2 != 0`.** Verified at `(2,4), (3,6), (4,6), (2,6)`. The rigidity conclusion
is false for every non-coprime pair. This is the arithmetic gate
`MINIMAL_CORE.md` predicted, in exactly that shape: it is `gcd(m,n) = 1`, not
`m = 2`, and it is not about sporadicity.

Structurally it is obvious in hindsight: if `g > 1` then `C^m` and `C^n` are both
powers of `C^g`, so the two families are not independent and stacking them buys
nothing new.

**GATE (B2) `N_Q >= D_P + D_Q`** (F3). At the critical slope: `(2,3)` needs
`N_Q >= 5` and has `15`; `(3,5)` needs `N_Q >= 8` and has `29`. Both comfortably
satisfied — the second hidden hypothesis, and the one that would bite first at
small `ell`.

**LEMMA B also re-proves an audited step independently (F10).** At `(2,3)` the
state `v_t(h_1) >= 1` induces the linear weight `w(i) = i = (m-1)*i`, and
LEMMA B then gives `Psi = 1+B_1 v`, i.e. `v_t(h_2) >= 3 = 2*2-1` — the audit's
level-4 step, by a completely different route (a `(D_P+1)`-term recurrence plus
unique factorisation, versus a perfect-square jet).

---

## 3. The lemma, assembled

**MAIN LEMMA (the cascade).** Let `char = 0`, `m < n`, `gcd(m,n) = 1`, and let
`H = 1 + sum_{k>=1} h_k u^k` over `Frac(T)` satisfy `(P)`, `(P0)`, `(Q)` of §2.0
with ranges `N_P = m*ell`, `N_Q = (n+1)*ell-1`, and with `h_k` of `y`-degree
`<= lam*k`, `lam >= m`. Then

```
v_t(h_k)  >=  m*k - 1        for  k <= K(m,n,ell) ,
```

and this is the exact fixed point of the pair of families (LEMMA A), attained
(C4), with stacked cokernel `m*k-n` per level (Prop. 2). The forcing is by
LEMMA C where the step is clean and LEMMA B at the critical slope `q' = m-1`;
`K` is bounded above by the reach ceiling of §4.1, which is **sharp at `(2,3)`**.

**The three honest hypotheses**, as the brief asked for them:

| gate | statement | why the audited case satisfies it invisibly |
|---|---|---|
| **(H-cap)** | `lam >= m` | `lam in {2,3}`, `m = 2`; sub2 sits exactly at equality |
| **(H-gcd)** | `gcd(m,n) = 1` | `gcd(2,3)=1`; false for `(2,4),(3,6),(4,6),...` |
| **(H-range)** | `N_Q >= D_P + D_Q` | needs `15 >= 5`; only binds at small `ell` |

Plus one **ceiling**, which is not a hypothesis but a hard limit: §4.1.

---

## 4. Reach

### 4.1 The ceiling — where the profile is genuinely FALSE (new, PROVED)

The profile is not unconditional. For a unit `a` and an integer `d` put
```
H(u) = ( 1 + beta*u^d )^(1/m) ,        beta = a*t^(m*(d-1)) .
```
Then `H^m = 1 + beta*u^d` **exactly**, so `[u^k]H^m = 0` for all `k != 0,d` — so
`(P0)` holds for `d <= N_P`, and at `k = d` the P condition is met with
**equality** (`v_t(beta) = m*d-m = P_d`). And
`[u^(j*d)]H^n = binom(n/m,j)*beta^j` has `v_t = j*m*(d-1)`, which meets
`Q_{j*d} = m*j*d - n` iff `j*m <= n`. So with
```
j0 := floor(n/m) + 1 ,
```
**every** `(Q)` condition in range holds as soon as `j0*d > N_Q`. But
`v_t(h_d) = m*(d-1) = m*d - m < m*d - 1`.

**Therefore the profile can hold at best for `k <= floor(N_Q / j0)`.**

| case | `m,n,ell` | `gcd` | `N_P` | `N_Q` | `j0` | **ceiling** | first violated `k` |
|---|---|---|---|---|---|---|---|
| **(72,108)** | 2,3,4 | 1 | 8 | 15 | 2 | **7** | `d = 8`, `v_t(h_8) = 14 < 15` |
| **(75,125)** | 3,5,5 | 1 | 15 | 29 | 2 | **14** | `d = 15`, `v_t(h_15) = 42 < 44` |
| — | 3,4,4 | 1 | 12 | 19 | 2 | 9 | `d = 10`, `27 < 29` |
| — | 2,5,4 | 1 | 8 | 23 | 3 | 7 | `d = 8`, `14 < 15` |
| — | 4,5,4 | 1 | 16 | 23 | 2 | 11 | `d = 12`, `44 < 47` |
| — | 2,4,4 | **2** | 8 | 19 | 3 | 6 | `d = 7`, `12 < 13` (moot: (H-gcd) already fails) |

**This makes `slice_obstruction_audit.py` G9 sharp.** G9 records, as a scope
note, "this cascade can advance at most `h_1..h_7`". G1/G2 here show the
witness at `k = 8` **exists**, satisfies every imposed condition, and sits
**inside the degree caps in both windows** (`deg_y h_8 = 14 <= 2*8 = 16` sub2,
`<= 3*8 = 24` sub1). So `h_1..h_7` is not conservatism — it is the truth. The
audited `a_t >= 9` uses `k = 5`, **two indices inside the ceiling**, and is
untouched.

MUTATION G4: for **every** `d` at or below the ceiling the same construction
*violates* a `(Q)` condition at level `j0*d <= N_Q`. So G1 is not an accident of
one `d`, and `floor(N_Q/j0)` is exactly the crossover.

### 4.2 `(3,5)` — the answer to task 4

**Does the slice obstruction apply at `(3,5)` even though the Belyi layer does
not?** `passport_75_125.py` (81/81) shows the top-band gate `u*kappa = m+n-1`
FAILS on every primitive functional at `(75,125)`, so `D^3/A^5` gives no Belyi
map there. **That is a different mechanism with a different gate, and it does not
transfer.** The slice obstruction's gates are `gcd(m,n)=1`, `lam >= m`,
`N_Q >= D_P+D_Q` — all satisfied at `(3,5)`. So:

**YES, the slice obstruction does apply at `(3,5)`.** But the *engine* changes,
in three ways, all computed:

1. **Index 1 needs TWO steps, and one of them is a CUBE** (E10). The profile
   target at `k=1` is `m-1 = 2`, while `(P)` only gives `0`; the square step
   (level 2, `15*X_0^2`, `c_2 = 15`) advances `0 -> 1`, and the **cube** step
   (level 3, `15*X_1^3`) advances `1 -> 2`. Higher-power jets are a real
   `m >= 3` phenomenon. They are still single-component, so **still no
   branching** — that part of the `(2,3)` picture survives.
2. **Index 2's last step STALLS for the clean-jet engine** (E11/E12). The square
   step (level 4, `15*X_3^2`) advances `V = 3 -> 4`; then `V = 4 -> 5` has no
   clean jet at levels 6 or 7. The level-6 jet is an irreducible cubic in
   `(A1_2, X_4, g3_6)`:
   ```
   A1_2^4*X_4 - 42*A1_2^3*g3_6 - 9*A1_2^2*X_4^2 - 162*A1_2*X_4*g3_6
        - 27*X_4^3 - 27*g3_6^2  =  0
   ```
   because the pure term `c_3*X_4^3` at `t^12` is **tied** by `g3_6^2`
   (`v_t(h_3) = 6` and `2*6 = 3*4`). Collision, not absence.
3. **LEMMA B closes it** (F9/F4/F5). The stall state
   (`v_t(h_1)>=2`, `v_t(h_2)>=4`) induces the weight `w(i) = 2*i = (m-1)*i`
   for `i = 1..10` — **exactly the critical slope**. LEMMA B then says: `Psi^3`
   is a cubic `R`, `Psi^5` a quintic, `R^5 = (Psi^5)^3`, `gcd(3,5)=1`, so
   `R = (1+c*v)^3`. Computed independently (F4): the three conditions
   `[v^L]R^(5/3) = 0` for `L = 6,7,8` have solution set **exactly** the cube
   locus — verified both ways (the cubes satisfy them; pinning `r1` the *only*
   solution is `(r2,r3) = (r1^2/3, r1^3/27)`), with MUTATION F6 confirming four
   explicit non-cubes each violate one. Then `Psi = 1+c*v`, so `B_2 = 0` and
   `v_t(h_2) >= 5 = 3*2-1`. **Step taken.**

So at `(3,5)` the obstruction bites, but the proof needs a *variety* argument
(unique factorisation on the leading-order system) where `(2,3)` needed only a
perfect square. A general lemma stated as "the lowest jet is always a perfect
square of the first unresolved coefficient" would be **false at `(3,5)`** — that
is the hypothesis `(2,3)` satisfies invisibly, and gate E9 names it.

### 4.3 Where the engine does NOT reach

* **`gcd(m,n) > 1`: dead** (F7). `(2,4)`, `(3,6)`, `(4,6)`, `(2,6)`, ... The
  rigidity step has explicit counterexamples. No amount of extra levels helps.
* **`lam < m`: dead** (D5). The P family stops being absorbable and the whole
  "0 / 0 / stacked" architecture is gone.
* **beyond the ceiling `floor(N_Q/j0)`: dead** (§4.1), with a witness.
* **`(3,5)` indices `r >= 3`: NOT COMPUTED.** LEMMA C's square step gives
  `V = 3r-2` and LEMMA B closes `r <= 2`. For `r >= 3` the intermediate state is
  no longer a linear weight, so neither lemma applies as stated, and this lane
  did not resolve it. **Flagged, not claimed.**
* **the `(3,5)` degree cap `lam` is an unknown input.** `(75,125)` has selected
  multiplicity `q = 2` but this lane did not locate the analogue of `[Q3]`'s
  `lam in {2,3}` for that corner. `lam >= 3` is *required* by (H-cap) since
  `m = 3`. **If `lam < 3` at `(75,125)`, the obstruction does not apply there
  and §4.2 is void.** That is the single load-bearing unknown in the reach
  analysis and it should be checked before anyone builds on §4.2.

---

## 5. The `(2,3)` positive control — not optional, and complete

Every audited `(2,3)` fact, reproduced by the general machinery specialised at
`(m,n) = (2,3)` (checker H1 gates on all 64 checks passing):

| audited fact | source | reproduced by | value |
|---|---|---|---|
| cokernels `0 / 0 / 2n-3`, both windows | audit E4, §2 | **D3** | `k=2..8: 0/0/1,3,5,7,9,11,13` |
| `2H^3-3H^2+1 = 3K^2+2K^3` | §2 | **A10** | `c_2=3`, `c_3=2`, `=(H-1)^2(2H+1)` |
| fresh coefficient cancels in `2r_n-3p_n` | audit E2, §2 | **D1** | `k=2..10` |
| stacked family `= [u^n](3K^2+2K^3)` | audit E3 | **D2** | exact |
| odd levels `3,5,7,9` empty | audit F1, §3 | **E4** | all `None` |
| even levels fire at `t^(2L-4)` | audit F2, §3 | **E5** | `t^0,t^4,t^8,t^12,t^16` |
| jet `= 3*([t^(2m-2)]h_m)^2`, linear in the fresh `g` with unit coefficient | audit F2/F3 | **E1** | `3*X_(2r-2)^2` at level `2r`, `r=1..5` |
| the `(3/4)*(...)^2` rows | §3 table | **E3** | `3*(g-A^2/2)^2 = (3/4)*(A^2-2g)^2` — **same `3`** |
| `v_t(h_k) >= 2k-1`, `k=1..5` | audit F4/F5, §3 | **E1/E2** | `1,3,5,7,9` |
| ordering load-bearing; would branch without it | audit F7 | **E6** | `t^14` jet, **2** coprime factors |
| no reliance on `(P0)` | audit F6 | **E7** | identical jet with `g_9,g_10` free |
| profile satisfiable and sharp | §3.2, audit F8/F9 | **C1/C2/C4** | equality at every `k>=r` |
| "advance at most `h_1..h_7`" | audit G9 | **G1/G2** | **and it is SHARP** |

**The `3/4` reconciliation, explicitly.** In this file's parametrisation
`h_r = t^V*X` and the jet is literally `c_2*X^2 = 3*X^2`. In the audit's
parametrisation `h_2 = -q_2/2 + t^2 g_2`, so
`X = [t^2]h_2 = g2_0 - g1_1^2/2` and
`3*X^2 = (3/4)*(g1_1^2 - 2*g2_0)^2`. **The coefficient is `c_2 = 3` at every
even level; the `3/4` is denominator-clearing inside the square, not a second
constant.** (E3, an exact polynomial identity.)

**And `a_t >= 9` is untouched.** `v_t(h_5) >= 9` is reproduced (E1, `r=5`,
level 10, `t^16`, `3*X_8^2`); with `[Q8]`/S3.4's `h_5 = dm1 = e` and
`a_t = v_t(e)` `[QC1]`, `a_t >= 9`. Nothing in this file weakens either import —
and `k=5` is two indices inside the newly-proved ceiling `k <= 7`.

---

## 6. PROVED / CHECKED / INFERRED

**PROVED here** — symbolic identity, exact linear algebra over `Q`, exact
polynomial factorisation, or a written proof whose every computational step is
machine-checked:

* the contact identity, its closed form, `c_2 = m*n*(n-m)/2`, and contact order
  exactly `2` with the `(1,2)` exception — **symbolically in `(m,n)`** (§1);
* the slice-exponent family `t^(m*k-r) | [u^k]H^r`, hence `P_k = m*k-m`,
  `Q_k = m*k-n`, `N_P = m*ell`, `N_Q = (n+1)*ell-1` (§2.0);
* **LEMMA A**: the profile `v_t(h_k) >= m*k-1` is the fixed point of both
  families, sharply, for every `(m,n)`, and is attained (§2.1);
* the elimination `m*r_k - n*p_k = [u^k](sum c_j K^j)` (§2.2);
* **Prop. 2**: cokernels `0 / 0 / (m*k-n)` by exact rank, and the gate
  `lam >= m` with its break levels (§2.3);
* **LEMMA C** and gate E9: the clean square step completes the profile iff
  `(m,n) = (2,3)` (§2.4);
* **LEMMA B**: leading-order rigidity, with the `(D_P+1)`-term recurrence
  (verified against honest series expansions) and the unique-factorisation step;
  and the gate `gcd(m,n) = 1` with explicit counterexamples for `g>1` (§2.5);
* the reach ceiling `k <= floor(N_Q/j0)`, with a witness at the first violated
  index and a mutation showing every smaller `d` fails (§4.1);
* the `(3,5)` computations: two-step index 1 with a cube, the index-2 stall and
  its explicit cubic, the induced critical slope, and the cube-locus rigidity
  closing the step (§4.2);
* the full `(2,3)` positive control of §5, including the `3/4` reconciliation
  and the negative/robustness controls E6/E7.

**CHECKED** — reproduced from an existing artefact without re-proving it:

| tag | statement | source |
|---|---|---|
| [Q3] | `lam in {2,3}`, `D_j = C_j*C4^(7-2j)`, `C4 = y^7*(y+1)` both windows | `window_caps_verify.py` W2/W5 |
| [QQ1] | `Q = C^3 + lambda*C^-1 + F`, `v_{1,0}(F) = -5` | `PROOF_INVENTORY.md` C3 (**2/4**) |
| [Q8] | `h_5 = dm1 = e` (the shift is triangular across zero) | `SLICE_OBSTRUCTION.md` §3.1 |
| [QC1] | `a_t = v_t(e)` | `divisor_filter.py` |
| [C75] | `(75,125)`: `(a,b,t,kappa) = (3,5,5,3)`, `C = y^2(y^3+1)`, `s = -11` | `C_SERIES_75_125.md` §2–3 (**[judgment 1/2]**) |
| [P75] | the Belyi top-band gate `u*kappa = m+n-1` fails on every primitive functional at `(75,125)` | `passport_75_125.py` (81/81) |

The two `(72,108)` imports this file *depends on for its interpretation*, and
does not re-prove, are `[QQ1]` and `[Q8]` — the same two
`SLICE_OBSTRUCTION.md` flags. This file does not improve their status. The
`(3,5)` reach analysis inherits `[C75]`'s **judgment 1/2** conditionality: if the
`(5,20)` reduction is not the standard type-II.b chart, `ell` and hence `N_P`,
`N_Q` change and §4.2's numbers move.

**INFERRED** — nothing here is asserted as a result. The forward-looking items
are labelled as such: `(3,5)` indices `r >= 3` (§4.3), and the value of `lam` at
`(75,125)` (§4.3), which is the one gate this lane could not evaluate.

---

## 7. What this does NOT do

1. **Nothing is entered into the ledger, the DAG, or the frontier.** No existing
   file was written. No stage record is emitted. `a_t >= 9` and the `(72,108)`
   exclusion are exactly where `SLICE_OBSTRUCTION.md` left them; this file
   *explains* the mechanism and *bounds* it, it does not extend the kill.
2. **It does not upgrade `[QQ1]` or `[Q8]`.** The Q column still rests on the
   `alpha`-strip WLOG at confidence 2/4, and the census delta still rests on
   `h_5 = dm1`.
3. **It does not prove the `(3,5)` slice obstruction as a whole** — only indices
   `1` and `2`, and only conditionally on `lam >= 3` at `(75,125)`, which is
   **unverified**. §4.2 is a mechanism result, not a corner kill.
4. **It says nothing about the Belyi layer.** The `(75,125)` `u*kappa` failure
   `[P75]` is cited, not re-derived, and the two mechanisms are shown to have
   *different* gates — which is an argument that they are independent, not that
   either is right.
5. **The alternate regime is untouched.** `a_t >= 9` is the wrong shape for
   `a_t >= 11` (`SLICE_OBSTRUCTION.md` §6.2) and nothing here changes that.
6. **The unused slack is still unused.** The `y`-order conditions, the exact
   equations `r_n = 0` for `n = 13,14,15` (audit G8's correction), and the
   `Phi` relation at `n = 17` are not exploited. In particular the ceiling of
   §4.1 is a ceiling **for the imposed set only**; the witness there has
   `r_13 != 0`, so it is not a genuine `(P,Q)` pair — the same honest caveat the
   audit records at G7. It bounds the strength of the conditions actually used.
7. **No Gröbner basis is load-bearing.** F4's variety statement is confirmed by
   direct solve with `r1` pinned at four values, and the *proof* (§2.5 step 2) is
   unique factorisation, not a computation.

---

## 8. Reproduce

```
cd d2_plane_72_108
python contact_lemma.py                # full report, 64 checks, ~90 s
python contact_lemma.py --quiet        # exit 0 iff every check passes
```

Read-only, pure sympy, no external tools. The expensive parts are the `(2,3)`
`r=5` level-10 jet (~8 s) and the `(3,5)` level-7 stall scan (~60 s).
