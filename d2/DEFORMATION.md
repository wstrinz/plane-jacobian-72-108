# DEFORMATION — the family is IDENTIFIED, and the discriminant identity is Szegő 6.71.5

> **STATUS (2026-07-24): THE FAMILY IS NAMED — AND THE NAMING IS BAD NEWS FOR THE
> NOVELTY CLAIM.** `CORNER_RESOLVENT.md` §2.5.2 asked for the family to be
> identified in closed form, flagging that "if line 3 lets you IDENTIFY the family
> against a known classical one, say so prominently — and then immediately check
> whether that family's discriminant is already published". It does, and it is.
>
> **The monic marked polynomial is a Jacobi polynomial:**
> ```
> g/L  =  P_n^(alpha,beta)(1 + 2y),      alpha = -n - 1/c,   beta = 1/c
> ```
> so `alpha + beta = -n` identically — a one-dimensional slice of the Jacobi
> family. **Szegő's classical discriminant formula (6.71.5), restricted to that
> slice, reproduces our identity exactly** — sign, powers, everything. Verified
> symbolically in `c` for `n = 2..8` and numerically at 49 `(n,c)` points, zero
> residual.
>
> **What this does and does not kill.** It does NOT touch the correctness of
> anything in `CORNER_RESOLVENT.md` §2.5 — the identity is true and the ODE
> derivation is valid. It DOES mean the identity is a **specialisation of a
> published formula**, not a new discriminant. The prior-art table in
> §2.5.1 lists "Classical Jacobi (Szegő 6.71.5) | quadruple product over `j` |
> **different**". That verdict was **wrong**, and for an understandable reason:
> the fourth product `prod (n+j+alpha+beta)^(n-j)` collapses to `prod j^(n-j)`
> precisely when `alpha+beta = -n`, which is the only line we ever sit on. The
> "quadruple product" is a triple product here, and it is ours.
>
> **The one thing that survives as prospectively new** is the *route*: the
> inhomogeneous first-order ODE `c y(y+1) g' - (cn y + cn + 1) g = 1/2` as a
> generator, and the observation that `Delta = cn+1` is the sole
> non-leading-coefficient content. No ODE of that shape was found in the
> orthogonal-polynomial literature. But this must now be written as *"a
> first-order ODE characterisation of a known Jacobi slice"*, not as a new
> discriminant formula. See §8.

Script: `deformation_probe.py` (NEW). `python deformation_probe.py --quiet`
runs **230 checks**, exit 0 iff all pass. Nothing existing was modified.

---

## 1. The closed form (S1)

The ODE is nothing but a two-term recursion. Writing `g = sum_m a_m y^m`,
comparing coefficients of `y^m` in
`c y(y+1) g' - (cn y + cn + 1) g = 1/2` gives, with **no residue**:

```
a_0 = -1/(2(cn+1)) = -1/(2*Delta)
a_m = -c(n-m+1)/(c(n-m)+1) * a_{m-1}          m = 1..n
```

Verified against a brute-force undetermined-coefficients solve of the ODE,
**symbolically in `c`**, for `n = 1..6`, together with
`L = (-1)^(n-1)(n!/2) c^n / prod_{k=1..n}(kc+1)`.

Existence is exactly the statement that no denominator vanishes:
**`g` of degree `n` exists iff `c != -1/k` for `k = 0..n`**, i.e. iff
`prod_{k=1..n}(kc+1) != 0` and `Delta != 0`.

## 2. THE IDENTIFICATION — three equivalent names (S2)

All three verified **symbolically in `c`**, `n = 1..7`, zero residual.

| # | statement |
|---|---|
| (a) | `y^n (g/L)(1/y) = sum_{j=0}^n binom(-1/c, j) y^j` |
| (b) | `g(y) = g(0) * 2F1(-n, 1; 1-n-1/c; -y)` |
| (c) | `g/L = P_n^(-n-1/c, 1/c)(1 + 2y)` |

**(a) is the one to hold in your head.** *The marked polynomial is the reversal
of the degree-`n` truncation of the binomial series `(1+y)^(-1/c)`.*

That single sentence explains everything the file was chasing:

* **why `c = 1` is cyclotomic** — `binom(-1,j) = (-1)^j`, so the truncation is
  `1 - y + y^2 - ...`, whose reversal is the alternating sum. Not a coincidence,
  an evaluation.
* **why `Delta = cn+1` deforms `n+1`** — `Delta/c = n + 1/c`, and `1/c` is the
  (negated) binomial exponent. `Delta` is `c` times "degree plus exponent".
* **why exactly one linear factor survives** — see §3.

**(c) is the one that decides novelty.** `alpha + beta = -n` identically. This is
a genuinely distinguished line in the Jacobi parameter plane: the leading
coefficient is `binom(2n+alpha+beta, n)/2^n = binom(n,n)/2^n = 2^(-n)`, so
`alpha+beta = -n` is the **last** line on which `P_n^(alpha,beta)` still has
degree `n` — one step further (`alpha+beta = -n-1, ..., -2n`) and the degree
drops. It is also where the standard three-term recurrence degenerates: sympy's
`jacobi()` raises `ZeroDivisionError` on these parameters, because the recurrence
coefficient `(a+b+2i-1)` vanishes. That degeneracy is why this slice is easy to
overlook, and is the most charitable explanation for §2.5.1's mis-verdict.

## 3. THE DISCRIMINANT IDENTITY IS SZEGŐ 6.71.5 ON THAT LINE (S3)

Szegő, *Orthogonal Polynomials*, (6.71.5):

```
disc_x P_n^(a,b) = 2^(-n(n-1)) prod_{j=1}^n j^(j-2n+2) (j+a)^(j-1) (j+b)^(j-1) (n+j+a+b)^(n-j)
```

Set `a = -n-1/c`, `b = 1/c`. Then `n+j+a+b = j`, so the fourth product merges
into the first:

```
disc_y(g/L) = 2^(n(n-1)) * disc_x P_n^(a,b)
            = prod_{j=1}^n j^(2-n) (j - n - 1/c)^(j-1) (j + 1/c)^(j-1)
```

and collecting `(cj+1)` powers — exponent `(j-1) + (n-j-1) = n-2` for `j < n`,
and `n-1` for `j = n` — gives exactly

```
disc(g/L) = (-1)^(n(n-1)/2) * Delta * (prod_{k=1}^n (kc+1))^(n-2) / (n!^(n-2) c^(n(n-1)))
```

which is `CORNER_RESOLVENT.md` §2.5's identity divided by `L^(2n-2)`. The sign
`(-1)^(n(n-1)/2)` falls straight out of `prod_{k=0}^{n-1} (-1)^(n-k-1)`.

**Machine verification (all zero-residual):**

| check | range | result |
|---|---|---|
| `2^(n(n-1)) * Szegő == our monic identity`, symbolic in `c` | `n = 2..8` | IDENTICAL |
| `our monic identity * L^(2n-2) == the published identity`, symbolic in `c` | `n = 2..8` | IDENTICAL |
| full chain against `sympy.discriminant` | 49 `(n,c)` points | OK |

*Methodological note, worth recording.* The first run of this check reported
FAIL everywhere for `n >= 5`. The cause was `2**(-n*(n-1))` in Python evaluating
to a **float**; `sp.Integer(2)**(-k)` fixes it. A float prefactor times an exact
integer of 40+ digits silently loses precision and produces a confident wrong
answer. Any Szegő-formula check in this repo must use exact rationals.

## 4. GALOIS GROUPS ACROSS THE FAMILY (S4)

`sympy.galois_group` on the primitive integer model.

**`n = 4`, integer `c` in `[-12, 32]`:**

| group | `c` values |
|---|---|
| **A4** | 2, 6, 12, 20, 30 |
| **C4** | 1 |
| **S4** | everything else (37 values) |

**`n = 4`, rational `c = p/q`, `\|p\| <= 60`, `q <= 12` (894 members, `--full`):**

| group | count | `c` values (all of them, where few) |
|---|---:|---|
| S4 | 875 | generic |
| **A4** | 15 | 2, 6, 12, 20, 30, 42, 56, 3/4, 15/4, 35/4, -2/9, 4/9, 10/9, 28/9, 40/9 |
| **C4** | 2 | **1**, **-1/5** |
| **D4** | 2 | **1/3**, **-1/7** |
| V4 | **0** | never observed |
| reducible | 0 | — |

**Findings.**

1. **A4 occurs exactly on the square-`Delta` locus.** The forward direction is
   forced by the theorem (`n` even ⇒ `[disc] = [(-1)^(n/2) Delta]`, so square
   `Delta` ⇒ square disc ⇒ `Gal <= A4`). The *converse* is the content: of the 15
   swept members with `4c+1` a rational square, **all 15 are A4 and none is V4**.
   Over the integers `4c+1 = square` means `c = j(j+1)`: `c = 2, 6, 12, 20, 30,
   42, 56` give `Delta = 9, 25, 49, 81, 121, 169, 225`. So `(a,c) = (12,2)`,
   the GGHV-new A4 case, is the **first member of an infinite A4 sub-family**,
   not an isolated accident.
2. **The other transitive subgroups DO occur — the question in the brief is
   answered YES.** `C4` at `c = 1` (that is `Phi_10`, Galois `(Z/10)^x = C4`) and
   at `c = -1/5`; `D4` at `c = 1/3` and `c = -1/7`.
3. **Every non-generic, non-A4 case has `1/c` an integer** (`1, 1/3, -1/5, -1/7`
   ⇒ `-1/c = -1, -3, 5, 7`). Under §2(a) these are exactly the members that are
   genuine **integer-exponent** truncated binomials. That is where the
   truncated-binomial literature (Filaseta–Moy and successors) lives, and it is
   where the arithmetic is special.
4. **Other degrees, integer `c` in `[-10,12]`:** `n=3` is S3 everywhere except
   `c=1` (reducible); `n=5` is S5 everywhere except `c=1` (reducible); `n=6` is
   S6 everywhere except `c=1`, which is `Phi_14` with group **C6**. No cyclic
   cubic (`C3`) was found: for every integer `c` tested at `n=3` the discriminant
   is **negative**, so a square discriminant — hence `C3` — is impossible there.

## 5. IRREDUCIBILITY (S5)

| sweep | members | reducible |
|---|---:|---:|
| integer `c` in `[-15,15]`, `n = 2..10` | 261 | 7 |
| rational `c = p/q`, `q <= 7`, `\|p\| <= 20`, `n = 2..8` | 1023 | 7 |

**Answer to the brief's question: reducibility never occurs at `c >= 2`.** Every
reducible member found sits at `c = 1` or at `c < 0`.

| `n` | `c` | `Delta` | factorisation | reason |
|---|---|---:|---|---|
| 3,5,7,9,11 | 1 | `n+1` | cyclotomic | `n+1` composite — the classical case |
| 8 | 1 | 9 | `(y^2-y+1)(y^6-y^3+1)` | ditto |
| 2 | -5, -13, -25, -5/2, -17/2 | -9, -25, -49, -4, -16 | two linear | `n=2`: `Delta = -square` ⟺ square disc |
| 3,5,7,… | `-2/n` (odd `n`) | **-1** | `(2y+1) * (…)` | see below |
| 3 | -1/4 | 1/4 | `(y+2)(y^2+2y+2)` | integer-exponent truncated binomial |
| 5 | -1/6 | 1/6 | `(y+2)(y^2+y+1)(y^2+3y+3)` | ditto |

**An infinite reducible family, apparently not previously noted here.** On the
line `Delta = -1`, i.e. `c = -2/n`:

```
n ODD  :  g(-1/2) = 0        -- rational root, g REDUCIBLE   (verified n = 3..21)
n EVEN :  g(-1/2) = (-1)^(n/2) / 2^n  != 0, and c = -2/n = -1/(n/2) is itself
          a degeneracy point of the ODE  (verified n = 2..14)
```

So the `Delta = -1` line is reducible for every odd `n` and does not exist for
even `n`. This is a clean parity dichotomy that mirrors the even/odd split in the
discriminant law itself.

## 6. ROOT GEOMETRY (S6)

The brief asked whether the roots stay on or near a circle. **They stay in a thin
annulus, and the geometric mean of their moduli has an exact closed form:**

```
(prod |root_i|)^(1/n) = ( prod_{k=1}^{n-1} (k + 1/c) / (n! * c) )^(1/n)
```

verified to 25 digits at 32 `(n,c)` points. At `c = 1` the product telescopes to
`n!/n! = 1` — **the unit circle, exactly**, recovering the roots of unity.

| `n` | `c` | min \|root\| | max \|root\| | max/min | geom. mean |
|---:|---:|---:|---:|---:|---:|
| 4 | 1 | 1.0000000 | 1.0000000 | 1.000 | 1.0000000 |
| 4 | 2 | 0.7030070 | 0.7438230 | 1.058 | 0.7231269 |
| 4 | 10 | 0.3929250 | 0.4396140 | 1.119 | 0.4156144 |
| 4 | 1/2 | 1.4450458 | 1.5474028 | 1.071 | 1.4953488 |
| 4 | -2 | 0.3375354 | 0.8177525 | **2.423** | 0.4445699 |
| 12 | 2 | 0.8422605 | 0.8939163 | 1.061 | 0.8589001 |
| 12 | 10 | 0.6629199 | 0.7414379 | 1.118 | 0.6876949 |
| 12 | -2 | 0.6079182 | 0.9326482 | 1.534 | 0.6614017 |

**Description.** For `c > 0` the roots sit on a near-circle of radius `< 1`
(`c > 1`) or `> 1` (`c < 1`), with modulus spread never exceeding ~12% across the
whole positive-`c` sweep. The circle inflates back toward radius 1 as `n` grows at
fixed `c`. For `c < 0` the annulus is genuinely wider (up to 2.4x at `n=4`,
narrowing with `n`). This is the expected Jentzsch/Szegő picture for partial sums
of a power series with radius of convergence 1 (`(1+y)^(-1/c)` is singular at
`y = -1`): the zeros of the truncations cluster on the circle of convergence, so
the reversals cluster on the unit circle too. **The `c = 1` unit circle is the
only exact case; everywhere else it is approximate but tight.**

## 7. OTHER SPECIAL VALUES OF `c` (S7)

| locus | what happens | verified |
|---|---|---|
| `c = 1` | monic `g` = alternating sum `(y^(n+1) ± 1)/(y ± 1)`; `\|disc\| = (n+1)^(n-1)` | `n = 2..8` |
| `c -> infinity` | `c * (reversal - 1) -> -sum_{j=1}^n (-y)^j / j` = **truncated logarithm**; the family itself collapses (monic `g -> y^n`) | `n = 3..6` |
| `c = -1/n` (`Delta = 0`) | monic `g = (y+1)^n`, `disc = 0` | `n = 2..8` |
| `c = -1/m`, `m > n` | reversal is **exactly** `sum_{j<=n} binom(m,j) y^j` — the Filaseta–Moy truncated binomial expansion of `(1+y)^m` | `n = 3..6`, `m = n+1..n+4` |
| `c = -1/m`, `m < n` | the ODE has **no** degree-`n` polynomial solution (the truncation of `(1+y)^m` is `(1+y)^m` itself, degree `m`) | `n = 5,6,7` |

Three of these are worth calling out.

**`c = -1/n` explains `Delta = 0` completely.** The brief asked "`Delta = 0` when
`c = -1/n` — what happens there?". Answer: the ODE loses its polynomial solution
(`a_0 = -1/(2*Delta)` blows up), but the **monic** normalisation is regular there
— the `(cn+1)` denominator never appears in `a_m/a_n` — and its limit is
`(y+1)^n`. The `n`-fold root at `y = -1` is where the discriminant's `Delta`
factor vanishes. `Delta` is not an abstract linear form; it is the equation of the
wall where the deformed cyclotomic collapses onto a single point.

**`c = -1/m` is the truncated-binomial literature.** `CORNER_RESOLVENT.md`
§2.5.1 lists Filaseta–Moy's truncated binomial as "does NOT subsume ours (one
factor vs many)". The relation is now exact rather than analogical: their
`P_{n,k}` **is** our member at `c = -1/n_theirs`. It does not subsume ours (we
have a full rational parameter, they have a subfamily indexed by integers), but
it is not a distant neighbour either — it is a sub-line. And it is a sub-line
whose Galois groups have received serious attention (Filaseta–Moy; and
`arXiv:2304.12658`, which uses Faltings' theorem to prove the `k=6` truncation
has group `S6` for all but finitely many `n`). None of that literature appears to
use the Jacobi identification, which is why the discriminant question there is
treated as hard.

**`c -> infinity` is the truncated logarithm.** §2.5.1's fifth "neighbour"
(arXiv:2401.14138, truncated logarithmic, "no closed form; explicitly
intractable") is literally the boundary point of our line. That is not evidence
that we found something new near it; it is evidence the whole neighbourhood is
one object.

## 8. THE ARITHMETIC OF `Delta` (S8)

**Primality.** `Delta = cn+1` is a linear form in `c` with `gcd(1,n) = 1`, so
Dirichlet gives infinitely many prime `Delta` for every `n`; nothing special.
Counts of prime `Delta` for `c = 1..39`: `n=3`: 13; `n=4`: 17; `n=5`: 10;
`n=6`: 24. The four GGHV-shaped cases `Delta = 9, 13, 17, 29` are three primes and
one square, which is unremarkable at these densities.

**Ramification — one conjecture refuted, a sharper one survives.**

The natural guess is *"every prime dividing `Delta` ramifies in `K = Q[y]/(g)`"*.
**REFUTED**, 5 counterexamples in 124 computed maximal orders
(`sympy round_two`):

| `n` | `c` | `Delta` | unramified `p` | `v_p(Delta)` | excused because |
|---:|---:|---:|---:|---:|---|
| 2 | 4 | 9 | 3 | 2 | `v_p` **even** |
| 2 | 12 | 25 | 5 | 2 | `v_p` **even** |
| 3 | -3 | -8 | 2 | 3 | `p \| n!` |
| 3 | -9 | -26 | 2 | 1 | `p \| n!` |
| 3 | 13 | 40 | 2 | 3 | `p \| n!` |

Every counterexample has either even `v_p(Delta)` (where the square-class law says
nothing) or `p | n! c` (where `L`'s own primes interfere). The surviving statement,
**63 confirmations and 0 counterexamples**:

> **CONJECTURE (well-supported, not proved).** If `p` exactly divides `Delta`
> (`v_p(Delta) = 1`) and `p` divides neither `n!` nor `c`, then `p` is
> **totally ramified** in `K`: `v_p(d_K) = n - 1`.

This is stronger than the squarefree law and does not follow from it: the law
constrains the square class of `disc(g)`, whereas `disc(g) = [O_K:Z[alpha]]^2 d_K`
leaves open whether a given `p` sits in the index rather than in `d_K`.

For `n = 4` the table is very clean — `Delta` prime always gives `p^3 \|\| d_K`:

| `c` | `Delta` | `d_K` factorisation | `v_p(d_K)` for `p \| Delta` |
|---:|---:|---|---|
| 3 | 13 | `3^4 5^2 7^2 13^3` | `13: 3` |
| 4 | 17 | `2^6 3^2 5^2 13^2 17^3` | `17: 3` |
| 7 | 29 | `2^2 5^2 11^2 29^3` | `29: 3` |
| 5 | 21 = 3·7 | `2^2 3^3 7^3 11^2` | `3: 3, 7: 3` |
| 2 | 9 = 3² | `2^6 3^4 5^2 7^2` | `3: 4` |
| 6 | 25 = 5² | `2^6 3^4 5^2 7^2 13^2 19^2` | `5: 2` |
| 12 | 49 = 7² | `2^6 3^4 7^2 13^2 37^2` | `7: 2` |

So the answer to the brief's question *"does the quartic field's ramification see
the factorisation of `Delta`?"* is **yes, and finely**: not just the square class,
but the exact ramification index at each simple prime factor.

---

## 9. HONESTY SECTION

### Computed and verified (symbolically or exactly; `--quiet` = 230 checks, exit 0)

* The coefficient recursion and `L`'s closed form, symbolic in `c`, `n = 1..6`.
* All three identifications of §2 (truncated binomial reversal, `2F1`, Jacobi),
  symbolic in `c`, `n = 1..7`.
* **Szegő 6.71.5 at `alpha+beta = -n` equals our discriminant identity**,
  symbolic in `c` for `n = 2..8`, plus 49 numeric points. Zero residual.
* All Galois-group and irreducibility tallies in §4–§5 (finite explicit sweeps).
* The geometric-mean root formula, to 25 digits, 32 points; exact `= 1` at `c=1`.
* Every §7 special-value statement, at the stated ranges.
* The 5 ramification counterexamples and the 63 confirmations of §8.

### Conjectural / unproved

* **§8's surviving ramification conjecture.** 63 confirmations, no proof. A proof
  would likely come from Newton polygons of `g` at `p`, not from the discriminant.
* **"A4 ⟺ square `Delta`" at `n = 4`.** The `⇒` direction is proved (the
  theorem forces `Gal <= A4`); the `⇐` direction (never `V4`) is 15/15 empirical.
  `V4` would need the resolvent cubic to split completely; nothing rules it out.
* **"`g` is irreducible for all `c >= 2`."** 0 counterexamples in ~1280 members,
  but this is exactly the kind of statement that is hard even for the
  integer-exponent sub-line, where it is an open research problem
  (Filaseta–Kumchev–Pasechnik). Do not state it as known.
* **"Every non-generic Galois group has `1/c` integral."** 4/4 data points. This
  is a pattern, not a finding; the sweep was `q <= 12`.

### Explicitly NOT established

* That the Szegő specialisation is *stated anywhere* in the literature. What is
  verified is that our identity **follows from** a published formula by
  substitution. Whether anyone has written down the `alpha+beta = -n` slice, or
  the resulting single-linear-factor form, is a literature question that a search
  did not settle: searches for "discriminant Jacobi alpha+beta = -n" and for
  discriminants of partial sums of binomial series with non-integer exponent
  returned nothing on point, and `arXiv:2304.12658` and
  `math/0409523` (both on the integer-exponent sub-line) contain **no** Jacobi
  identification and **no** discriminant formula. That is weak evidence the slice
  is unwritten, and it is not enough to claim it.
* Any statement about the **geometric** (GGHV/GGV) side. This file is entirely
  about the ODE family. The scope caveats of `CORNER_RESOLVENT.md` §5 are
  unchanged and unaddressed here.
* `V4` cannot occur — only that it was not observed.

### What this changes upstream

`CORNER_RESOLVENT.md` §2.5.1's neighbour table needs two corrections and
§2.5.2's framing needs rewriting. **Neither is done here** (this file creates no
edits to existing files):

1. The row *"Classical Jacobi (Szegő 6.71.5) | quadruple product over `j` |
   different"* is **wrong**. On our line the fourth product collapses and the
   formula is ours.
2. The row *"Truncated binomial (Filaseta-Moy) | does NOT subsume ours"* is
   half-wrong: it does not subsume, but it is a sub-line (`c = -1/m`), not a
   different family.
3. §2.5.2's "the deformation itself appears unwritten" should become "the
   deformation is the `alpha+beta = -n` Jacobi line; what may be unwritten is the
   first-order-ODE characterisation of it and the resulting one-linear-factor
   normal form."

The **single highest-value next check** is no longer a computation: it is to open
Szegő §6.71 and the Jacobi-parameter-degeneracy literature and find out whether
the `alpha+beta = -n` slice is named. If it is, the discriminant content of
`CORNER_RESOLVENT.md` §2.5 is a corollary of a textbook and only the ODE survives.

Checker: `deformation_probe.py` (`--quiet` for the gate, `--full` for wide sweeps).
