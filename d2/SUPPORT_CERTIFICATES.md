# SUPPORT_CERTIFICATES.md — the last irreducible machine step, as Bézout identities

Checker: `support_certificates.py` (**55/55**, `--quiet` exits 0 iff all pass;
`--no-xcheck` skips the `spine9_audit.py` cross-check). Read-only: writes no
file, mutates no repo artifact, changes no verdict.

---

## 0. Headline

`MINIMAL_CORE.md` §3 records nine machine steps in the (72,108) proof, of which
**exactly one** was not hand-checkable in principle: step 8, the marked-support
feasibility test over the 40 `(k,z)` pairs, decided in `spine9_audit.py` §E by a
saturated Gröbner basis returning the unit ideal.

**That step is now gone.** The whole test reduces to:

| `k` | cell | what decides it now | machine? |
|---:|---|---|---|
| 0 | `a9_b0000_T1` | degree dichotomy: an odd degree 7 cannot be cancelled | **no** |
| 1 | `a9_b1000_T1` | 5 Bézout identities in `Z[y]`, degree ≤ 4 | **no** |
| 2 | `a9_b1100_T1` | 5 Bézout identities in `Z[p]`, degree ≤ 6 | **no** |
| 3 | `a9_b1110_T1` | degree + `t`-valuation bookkeeping | **no** |
| 4 | `a9_b1111_T1` | support pins `z = 3` by monomial comparison, then `22 > 17` | **no** |

The external review's reduction of `k = 2` to short univariate Bézout identities
**holds**, and generalises: the same shape works for `k = 1`, and both are
cleaner than the review claimed, because a Bézout identity needs **no
irreducibility and no field theory at all** (§5).

One discrepancy with the review, reported not reconciled: its stated
`m_2(p) = -288p^3 + 504p^2 - 73p - 8` is **not a rank minor of this problem**
(§6). Its `m_3` is exact. **No verdict is affected** — `z = 2` still dies, by a
certificate computed here.

---

## 1. The condition being decided

Downstream of the cofactor identity `F·Z = (1/6)γ⁵t⁹Π⁴` and `Z = ζ·t^z`
(`spine9_audit` C3/D6–D10), the surviving constraint on a marked factorisation
is

```
    q = Pi * Q ,   deg Pi = k ,   Pi^2 | ( mu*t^3*Q - 3*zeta*t^z ) ,   mu != 0 != zeta
```

with `t = y + 1`, `q = 2048y^4 - 512y^3 + 320y^2 - 240y + 195`, and `mu`, `zeta`
**scalars** of the coefficient field. The condition is homogeneous of degree 1
in `(mu, zeta)` jointly, so every scalar rescaling of `Pi` or `Q` is free
(checker A3, E13).

The valuation ledger (`spine9_audit` G3/G11, from the audited cascade rows
`h1..h5 = 1,3,5,7,9` alone) gives **`2 <= z <= 6`**. All tables below sweep
`z = 0..9` anyway — the extra five columns are free corroboration.

**The key observation.** `Pi^2 | X` with `mu, zeta` both nonzero is a
**rank-one condition** on the pair of columns of `X`'s coefficient vector in
`K[y]/(Pi^2)`. Rank can only *drop* under the further reduction
`K[y]/(Pi^2) -> K[y]/(Pi)`, so a nonvanishing minor of the **`Pi`-level**
`2 x 2` block already kills the case. That is what makes the certificates short:
one divides by `Pi`, not by `Pi^2`.

---

## 2. `k = 1` — five elementary gcds, each with a Bézout certificate

Take `Pi = y - r` with `q(r) = 0`, `Q = q/(y-r)`, `N := t^3 Q`. Then
`Pi^2 | (mu*N - 3*zeta*t^z)` is the `2 x 2` system

```
    mu*N(r)  = 3*zeta*t(r)^z
    mu*N'(r) = 3*zeta*z*t(r)^(z-1)
```

whose determinant vanishing is `z*N(r) - t(r)*N'(r) = 0`.

**Lemma (checker B1, PROVED here with `z` a free symbol, so for all `z` at once).**
Modulo `q(r) = 0`,

```
    z*N(r) - t(r)*N'(r)  =  -(t(r)^3 / 2) * n_z(r) ,
    n_z(y) := (y+1)*q''(y) - 2*(z-3)*q'(y) .
```

Since `q(-1) = 3315 != 0`, `t(r) != 0`, so **`k = 1` is infeasible at `z` iff
`gcd(q, n_z) = 1`**. All ten gcds are `1` (checker B3).

### The certificates: `a_z(y)·q(y) + c_z(y)·n_z(y) = N_z`

`q(y) = 2048y⁴ − 512y³ + 320y² − 240y + 195`.

| `z` | `n_z(y)` |
|---:|---|
| **2** | `40960y³ + 18432y² − 1152y + 160` |
| **3** | `24576y³ + 21504y² − 2432y + 640` |
| **4** | `8192y³ + 24576y² − 3712y + 1120` |
| **5** | `−8192y³ + 27648y² − 4992y + 1600` |
| **6** | `−24576y³ + 30720y² − 6272y + 2080` |

```
z = 2   a_2 = -33280*y^2 + 244224*y + 197176
        c_2 =   1664*y^3 -  13376*y^2 -   480*y +  195
        N_2 = 38480520

z = 3   a_3 =  39936*y^2 + 181632*y + 150512
        c_3 =  -3328*y^3 -  11392*y^2 +   360*y -  195
        N_3 = 29225040

z = 4   a_4 =     64*y^2 +    272*y +    221
        c_4 =    -16*y^3 -     16*y^2
        N_4 = 43095

z = 5   a_5 = -93184*y^2 + 246656*y + 166416
        c_5 = -23296*y^3 -  11136*y^2 -   840*y + 1365
        N_5 = 34635120

z = 6   a_6 = -44544*y^2 +  34560*y +  16712
        c_6 =  -3712*y^3 -    832*y^2 +           435
        N_6 = 4163640
```

Outside the window (free corroboration; same shape):

```
z = 0   n_0 =  73728*y^3 + 12288*y^2 +  1408*y -  800
        a_0 = -569088*y^2 + 1371072*y + 1181332
        c_0 =   15808*y^3 -   44672*y^2 -  13680*y + 3705      N_0 = 227395740
z = 1   n_1 =  57344*y^3 + 15360*y^2 +   128*y -  320
        a_1 = -1025024*y^2 + 2798720*y + 2303792
        c_1 =    36608*y^3 -  118912*y^2 -  19800*y + 6435     N_1 = 447180240
z = 7   n_7 = -40960*y^3 + 33792*y^2 -  7552*y + 2560
        a_7 = -2462720*y^2 + 1372544*y + 576176
        c_7 =  -123136*y^3 -    2176*y^2 +  13320*y + 21645    N_7 = 167765520
z = 8   n_8 = -57344*y^3 + 36864*y^2 -  8832*y + 3040
        a_8 =  -69888*y^2 +  41792*y +  21548
        c_8 =   -2496*y^3 +    512*y^2 +    720*y +   585      N_8 = 5980260
z = 9   n_9 = -73728*y^3 + 39936*y^2 - 10112*y + 3520
        a_9 = -6349824*y^2 + 4833408*y + 3149104
        c_9 =  -176384*y^3 +   82816*y^2 +  95400*y + 51675    N_9 = 795971280
```

`z = 4` is worth staring at: `(64y²+272y+221)·q + (−16y³−16y²)·n_4 = 43095`.
That is the entire kill of `a9_b1000_T1` at `z = 4`.

---

## 3. `k = 2` — the resolvent sextic and five rank minors

Normalise `q/2048 = Pi * Q_Pi` with both factors monic:

```
    Pi   = y^2 + p*y + B ,           Q_Pi = y^2 + (-1/4 - p)*y + D ,
    B = (128p^3 + 32p^2 + 20p + 15) / (32*(8p+1)) ,
    D = ( 64p^3 + 32p^2 + 14p -  5) / (16*(8p+1)) .
```

**Checker C1 (PROVED).** Matching coefficients, the `y^1..y^4` coefficients of
`2048(8p+1)^2 (Pi·Q_Pi − q/2048)` vanish identically and the `y^0` coefficient is
**exactly**

```
    r(p) = 32768p^6 + 24576p^5 + 16384p^4 + 5632p^3 - 10080p^2 - 2680p - 495 .
```

So the factorisation exists **iff `r(p) = 0`** — one residual, as the review said.
`deg r = 6` is the number of 2-subsets of the four roots of `q`, and `8p+1` is
invertible (`r(-1/8) = -2601/8 != 0`). `B`, `D` coincide with `spine9_audit`'s
own audited `s_of` / `ss_of` (checker C6), so this is the *same* `Pi`, not a
re-parametrisation.

**The minor.** With `u := rem((y+1)^3 Q_Pi, Pi)` and `w := rem(-3(y+1)^z, Pi)`
(each degree ≤ 1, hence a 2-vector over `Q[p]/(r)`), set

```
    m_z(p) := primitive integer form of  u_0*w_1 - u_1*w_0 .
```

`Pi | (mu t^3 Q_Pi − 3 zeta t^z)` says `mu·u = 3·zeta·w`; with `mu, zeta` both
nonzero that forces `m_z(p) = 0`. So **`m_z(p) != 0` kills the case**.

### The five `m_z`

```
m_2 = -1024000*p^5 -  274432*p^4 -  856320*p^3 +  201952*p^2 +  191224*p +   45771
m_3 =    -4096*p^5 +    3072*p^4 -    5312*p^3 +    2944*p^2 -     210*p -     313
m_4 =    14336*p^5 +   66048*p^4 -   18176*p^3 +   32560*p^2 -   14454*p -    9529
m_5 =   253952*p^5 +  297984*p^4 +   50560*p^3 +   99616*p^2 -   79488*p -   46789
m_6 = 32768000*p^5 + 22812672*p^4 + 11051776*p^3 + 4124320*p^2 - 7420416*p - 3468037
```

All five have degree exactly 5 and `gcd(m_z, r) = 1` (checker C9). *This corrects
the review's expectation that `z = 2` and `z = 3` give cubics: in the natural
minor family only one object anywhere in the problem is a cubic, and it is the
`Pi^2`-level minor at `z = 3` — see §6.*

### The certificates: `a_z(p)·r(p) + c_z(p)·m_z(p) = N_z`

```
z = 2   a_2 =    525568000*p^4 -    12971776*p^3 +   392825408*p^2 -  235872392*p -  82702785
        c_2 =     16818176*p^5 +     7691264*p^4 +     4542720*p^3 +     230560*p^2 -   6470240*p -  220950
        N_2 = 30824776125

z = 3   a_3 =       -32768*p^4 +       34816*p^3 -       53312*p^2 +      38448*p -     13665
        c_3 =      -262144*p^5 -      114688*p^4 -       94720*p^3 -      28800*p^2 +     68880*p -   18900
        N_3 = 12679875

z = 4   a_4 =       630784*p^4 +     2839808*p^3 -     1028160*p^2 +    1901336*p -    702455
        c_4 =     -1441792*p^5 -      929792*p^4 -      783360*p^3 -     371840*p^2 +    227120*p -   84600
        N_4 = 1153868625

z = 5   a_5 =      9586688*p^4 +    12423424*p^3 +     5153760*p^2 +    6981568*p -    335677
        c_5 =     -1236992*p^5 -     1079296*p^4 -      973056*p^3 -     571840*p^2 +     41080*p -   56340
        N_5 = 2802252375

z = 6   a_6 = 57876635648000*p^4 + 68505833800192*p^3 + 55317624278016*p^2 + 38185333328656*p + 8225371821707
        c_6 =   -57876635648*p^5 -    71620320256*p^4 -    66253893888*p^3 -    46308651040*p^2 -    8000500280*p - 1778836320
        N_6 = 2097511122958875
```

`z = 3` is the cheapest and `z = 6` the dearest; the largest integer anywhere is
16 digits. Outside the window `z = 0,1,7,8,9` also certify (cofactors up to 23
digits) — printed in full by `python support_certificates.py`.

**Robustness (checker C18).** The *swapped* assignment — taking `Pi` to be the
other quadratic factor, i.e. `p <-> -1/4 - p` — is also infeasible for every `z`,
so no choice of which factor is called `Pi` escapes.

---

## 4. `k = 3`, `k = 4`, `k = 0` — no certificate needed

### `k = 3`: structural (checker D1–D4)

`deg Pi = 3`, so `deg Q = 1` and `deg Pi^2 = 6`. Inputs: `deg q = 4`, `q`
squarefree, `q(-1) != 0`. Nothing else.

* **`z <= 5`.** `deg(mu t^3 Q − 3 zeta t^z) <= max(4, 5) = 5 < 6`, so the
  numerator must vanish *identically*: `mu t^3 Q = 3 zeta t^z`, i.e.
  `Q = (3 zeta/mu) t^(z-3)`. `deg Q = 1` forces `z = 4` and `Q = const·t`, hence
  `t | q`, i.e. `q(-1) = 0`. **Contradiction.** (`z = 3` would make `Q`
  constant, `z < 3` would make it non-polynomial; both absurd.)
* **`z = 6`.** The numerator has degree exactly 6 (the `t^3 Q` term has degree
  4 < 6, so the leading coefficient is `-3 zeta != 0`). With `deg Pi^2 = 6` the
  divisibility forces `numerator = const·Pi^2`. But
  `v_t(numerator) = min(3, 6) = 3` exactly, because `v_t(t^3 Q) = 3`
  (`Q(-1) != 0`, else `q(-1) = 0`). Hence `t^3 | Pi^2`, so `t | Pi`, so
  `q(-1) = 0`. **Contradiction.**

Since the whole admissible window `2 <= z <= 6` is covered, `k = 3` uses **no
arithmetic of `q` beyond `q(-1) != 0`** — so no certificate is possible or
needed. (`spine9_audit` E7 independently confirms this by finding the same kill
on three unrelated squarefree quartics.)

### `k = 4`: monomial comparison, then the degree ledger (checker D5–D8)

`Pi = q/2048`, `Q = 2048`, `deg Pi^2 = 8`. For every `z <= 6`,
`deg(2048 mu t^3 − 3 zeta t^z) = max(3, z) <= 6 < 8`, so the numerator must
vanish identically: `2048 mu t^3 = 3 zeta t^z`. That happens **iff `z = 3`**
(with `3 zeta = 2048 mu`, both nonzero — consistent). So the support test *pins*
`z = 3` by comparing two monomials, no Gröbner and no certificate.

Then `(*deg)` gives `deg F = 9 + 4k − z = 9 + 16 − 3 = 22`, while the certified
sub1 caps give `deg F <= 17`. **`22 > 17`.** Independent second route: `z <= 6`
alone gives `deg F >= 19 > 17`.

### `k = 0`: the degree dichotomy (checker D9–D12)

`Pi = 1`, so the support condition is vacuous. Eliminating `v` from
`Z = A^2 − gamma·v` in the boxed row gives (residual 0)

```
    gamma*u = mu*t^3*q - 6*A^2 + 3*zeta*t^z ,       deg(gamma*u) = deg d2 <= 6 .
```

The right side has three contributors with nonzero leading coefficients:
`mu t^3 q` of degree **7** (`lc = 2048 mu != 0`), `−6A^2` of **even** degree
`2 deg A`, and `3 zeta t^z` of degree `z <= 6`. A degree-7 term cannot be
cancelled by an even-degree square (7 is odd) nor by a term of degree ≤ 6.
So `deg(RHS) = 7 > 6`. **Contradiction.**

Both inputs are load-bearing with **zero margin**: the kill switches off at
`deg d2 <= 7` and at `z <= 7` (checker D12), so this is not a blanket argument.

---

## 5. Why a Bézout identity is stronger than "the ideal is the unit ideal"

Beyond checkability, the certificate form removes two hypotheses the
field-theoretic phrasing appears to need.

**Checker E10.** From `a_z(x) f(x) + c_z(x) m_z(x) = N_z` with `N_z` a nonzero
integer, specialise at **any** root `x_0` of `f`, in any commutative ring where
`N_z` is not a zero divisor:

```
    c_z(x_0) * m_z(x_0) = N_z != 0     =>     m_z(x_0) != 0 .
```

So the kill needs **neither** the irreducibility of `q` **nor** the
irreducibility of `r(p)` **nor** `K = Q[p]/(r)` being a field. Those facts
(themselves machine computations) only *explain* why one certificate suffices for
all six values of `p`; they are not load-bearing. What remains is one polynomial
expansion and the observation `N_z != 0`.

That is exactly the `ring_nf`-plus-`norm_num` shape a Lean proof wants, and it is
why "a nonzero constant on the right" is a genuine characteristic-0 fact rather
than reconnaissance.

---

## 6. Adjudication of the external review

| review claim | verdict |
|---|---|
| `B`, `D` are those two rational functions of `p` | **CONFIRMED** exactly (C1) |
| the factorisation leaves exactly one residual, the sextic `r(p)` | **CONFIRMED** exactly — the residual *is* `r(p)`, over denominator `2048(8p+1)^2` (C1) |
| `r(p)` is irreducible over `Q` | **CONFIRMED** (C2) — and shown **not load-bearing** (E10) |
| `m_3(p) = 320p^3 - 1104p^2 + 1460p - 401` is a rank minor coprime to `r` | **CONFIRMED exactly**: it is `-64/3` times the `Pi^2`-level minor on rows `(y^2, y^3)` at `z = 3`, and the only sub-quintic minor anywhere in the family (C15) |
| `m_2(p) = -288p^3 + 504p^2 - 73p - 8` is a rank minor coprime to `r` | **HALF FALSE — see below** |
| `z = 4,5,6` reduce to degree-5 polynomials, coprime | **CONFIRMED** (C9/C11) — and so do `z = 2,3` in the `Pi`-level family |
| `k = 2` reduces to five short univariate Bézout identities | **CONFIRMED and delivered** (§3) |
| `k = 1` needs only five elementary gcd checks | **CONFIRMED and delivered** (§2), with the reduction lemma proved for all `z` at once |

### The `m_2` discrepancy (checker C16)

`-288p^3 + 504p^2 - 73p - 8` **is** coprime to `r` (`gcd = 1`), so a Bézout
identity for it exists. But it is **not a rational multiple of any of the 56
nonzero `Pi^2`-level rank minors, at any `z` in `[0,9]`, in either
`Pi <-> Q_Pi` assignment** — and it is not the `Pi`-level minor at any `z`
either. So it is not a rank minor of this problem, and a certificate for it
would certify nothing.

Given that `m_3` matches *exactly* (to the factor `-3/64`), the review's
framework is demonstrably the same as the one here, which makes a computational
slip in the `z = 2` entry the likely explanation. **The verdict at `z = 2` is
unaffected:** it is killed by `m_2` as computed in §3, with a verified
certificate.

---

## 7. Agreement with the standing verdicts

`support_certificates.py` §E executes `spine9_audit.py` read-only in a private
namespace and compares tables pair by pair.

| source | result |
|---|---|
| `sub1_spine9.py` | 37/37 pass; records `k = 1,2,3` infeasible for every `z in [0,9]`, `k = 4` feasible only at `z = 3` |
| `spine9_audit.py` | 81/81 pass; same table, by **two** independent routes (saturated ideal over `Q`; rank test in an explicit splitting subfield) |
| **certificates (this file)** | `k=1`: none feasible. `k=2`: none feasible. `k=3`: none feasible for `z <= 6`. `k=4`: exactly `{3}` for `z <= 6`. |

**Checker E6/E7: zero disagreements** on every `(k,z)` pair both sides decide —
`k = 1, 2` for all ten `z`, and `k = 3, 4` for `z <= 6`. (For `k = 3, 4` the
structural arguments here are stated only on the admissible window `z <= 6`; the
audited machine sweep covers `z = 7,8,9` too and also finds them infeasible, so
nothing is lost.)

### Controls (a trivially-true check would be worse than none)

| control | what it shows |
|---|---|
| B7 | 30 single-term corruptions of the `k=1` certificates (`a_z+1`, `N_z+1`, `n_z+y`) all **fail** |
| C12 | 40 single-term corruptions of the `k=2` certificates (`a_z+1`, `N_z+1`, `m_z+p`, `r+p`) all **fail** |
| B8 | NON-VACUITY, `k=1`: on the synthetic squarefree quartic `y^4+y^3-9y^2+7` (`q(-1) = -2 != 0`) the same gcd test returns **FEASIBLE at exactly `z = 3`** |
| C14/C14b | NON-VACUITY, `k=2` — *the control the audited lane lacks for `k=2`*: on `3y(y+2)(2y^2+4y+1)` with `Pi = y^2+2y` the same `Pi`-level determinant **does vanish**, at `z = 7` (where `Pi^2` genuinely divides `t^3Q - 3t^7`) and at the other odd `z`, for the understood reason that `t = y+1` takes values `±1` at the roots `0, -2` |
| C8 | every `zeta` column is nonzero and every `m_z != 0` in `Q[p]`, so the minor is not nonzero-by-degeneracy |
| E11/E12 | numerical corroboration at 40 digits: `n_z` nonvanishing at all 4 roots of `q`, `m_z` at all 6 (distinct) roots of `r`, for every `z` |
| external | mutating `r`'s constant, or `q`'s constant, or stubbing `verify_bezout` to `return True`, each drives the checker to a **nonzero exit** — the last one is caught by B7/C12 |

---

## 8. PROVED / CHECKED / INFERRED

**PROVED here, from first principles, importing nothing:**

* the `k=1` reduction lemma `z N(r) − t(r)N'(r) = −(t(r)^3/2) n_z(r) mod q(r)`,
  as a polynomial identity in `r` with `z` a **free symbol** — hence for all `z`
  simultaneously (B0–B2);
* `q/2048 = Pi·Q_Pi` with the stated `B`, `D` leaves exactly one residual and
  that residual **is** `r(p)` (C1), together with the free monic normalisation
  (E13);
* `B`, `D` are the same objects as `spine9_audit`'s audited `s_of`/`ss_of` (C6);
* **all twenty Bézout identities** `a_z f + c_z m_z = N_z != 0` over `Z`
  (B4, C10) — each one expansion;
* the specialisation principle that makes irreducibility unnecessary (E10);
* the three structural arguments for `k = 3`, `k = 4`, `k = 0` on the admissible
  window, with their load-bearing controls (D1–D12);
* the arithmetic facts about `q` and `r`: irreducibility, squarefreeness,
  `q(-1) = 3315`, `r(-1/8) != 0`, `r` separable (A1, A2, C2, C4, E12).

**CHECKED here but resting on premises proved elsewhere in the repo:**

* the standing premise `a_t = 9` (`SUB1_SPINE9.md`, `AT_LE9_AUDIT.md`);
* that `Pi^2 | (mu t^3 Q − 3 zeta t^z)` **is** the surviving condition — i.e. the
  cofactor identity `F·Z = (1/6)γ⁵t⁹Π⁴`, `Π | B`, `gcd(Z,Π) = 1`, `Z = ζ t^z`
  (`spine9_audit` C3, D1–D10);
* the valuation window `2 <= z <= 6` (`spine9_audit` G3/G11, itself derived from
  the audited cascade rows `h1..h5 = 1,3,5,7,9`);
* the degree caps used by `k = 4` (`deg F <= 17`) and `k = 0` (`deg d2 <= 6`),
  read at runtime from `cascade_engine.SUB1` / `full_system_bridge.STRIP_DEGCAP`;
* the two published verdict tables, re-executed and compared (E5–E8).

**INFERRED (judgement, not proof):**

* that the review's `m_2` is a computational slip rather than a different
  (unstated) normalisation. What is *proved* is only that it is not a rank minor
  in any assignment tried (C16); the diagnosis of *why* is inference. It costs
  nothing, since `z = 2` is certified independently.
* that this closes `MINIMAL_CORE.md` §3 step 8 as "not a machine step". The
  certificates are still *produced* by machine; the claim is that *verifying*
  them needs none. Reproducing the `m_z` from scratch by hand is a page of
  polynomial remainder arithmetic — tedious, but no engine.

---

## 9. What a writeup now needs for step 8

An appendix of:

* the definition `n_z = (y+1)q'' − 2(z−3)q'` plus the one-line lemma of §2;
* five integer Bézout triples in `Z[y]`, degree ≤ 4 (§2);
* the two rational functions `B`, `D`, the sextic `r`, and the recipe
  `m_z = u_0w_1 − u_1w_0` (§3);
* five integer Bézout triples in `Z[p]`, degree ≤ 6 (§3);
* three short paragraphs for `k = 3`, `k = 4`, `k = 0` (§4).

No Gröbner basis, no unit-ideal claim, no splitting field, and — by §5 — not even
an irreducibility test.
