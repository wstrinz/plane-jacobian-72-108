# CORNER_RESOLVENT — the marked-polynomial generator and an arithmetic corner law

> **PRIOR-ART CORRECTION (2026-07-24, verified against primary sources).** Two of
> the three ingredients here are **NOT new** and must be presented as
> reinterpretation, not discovery:
>
> 1. **The corner integer `C = q(kappa+1) - t` is GGV's published `bl - a`.**
>    GGHV22 tex line 373 reads verbatim *"Here `A_0=(8,32)`, `A_1=(8,28)` and
>    `A_2=(11/4,7)`"*; with `(a,l,b) = (11,4,7)` that is `bl - a = 7*4 - 11 = 17`
>    — exactly our `C` for (72,108). The (66,99) corner `(11/3,8)` gives
>    `8*3 - 11 = 13`. GGV5 uses `bl - a` throughout its `(m,n)`-family machinery
>    (gcd / Diophantine family counts). **Both checked by hand against the
>    sources.** What GGV never do is attach it to a discriminant, a field, or a
>    Galois group.
> 2. **The forcing ODE is published in general form** — GGV1
>    (`arXiv:1401.1784`) line 1036, equation `eqq1`:
>    `z^h p = cp f + az p' f - bzf'p`, with `a := (1/rho) v_{rho,sigma}(F)` and
>    `b := (1/rho) v_{rho,sigma}(P)`, i.e. already our `(c,s)` parameterization.
>    Verified in the cached source. Valqui co-authors GGV1, GGV5 and GGHV22.
>
> **What remains genuinely ours:** (a) the DISCRIMINANT / RESOLVENT LAW of §3 —
> no occurrence of Galois, discriminant, resolvent or splitting field anywhere in
> the GGV corpus; GGHV reaches separability and stops; and (b) the CLASSIFICATION
> of §2 — that `deg f = 6` forces `deg g = 4`, that a polynomial solution exists
> iff `c | a` and is unique, and that exactly four `(a,c)` satisfy it, two of them
> not in the literature. The ODE is theirs; solving it as a *generator* and
> classifying its solutions is ours.
>
> Any write-up must cite GGV1, GGV5 and GGHV22 for (1) and (2), and lead with the
> resolvent law as the sole novel claim. See also Ramirez-Valqui
> `arXiv:2506.05697`, which carries a parity-indexed square-root ambiguity over
> marked roots — close enough that the boundary needs drawing explicitly — and
> the classical function-field Galois-case theorems (Campbell 1973, Razar 1979,
> Wright 1981, Bass-Connell-Wright 1982), which are a DIFFERENT Galois group and
> must be distinguished in one sentence or specialists will read this as a
> rediscovery.

> **STATUS (2026-07-24): ORIGINAL CONJECTURE REFUTED; A SHARPENED FORM SUPPORTED
> (27 data points, 9 delimiting counterexamples).** The instrument matters more
> than the law: the GGHV forcing ODE is reconstructed and reproduces BOTH
> published marked polynomials byte-exact, so marked polynomials can now be
> generated rather than read out of papers. The law's *interpretation* rests on a
> fitted dictionary — see §5, the named weak link.

## 1. Where this came from

An unpromising direction. A competing program (`xdLawless2/horseflow`) computes
the same GGHV case over a degree-5 field with Galois group S5 and quadratic
resolvent `Q(sqrt 663)`; ours is a degree-4 field with group S4 and resolvent
`Q(sqrt 17)`. The two are provably non-isomorphic (S5 has no subgroup of order
30, so no degree-4 subfield; S4 has no index-5 subgroup) — yet they share the
ramification set `{2,3,13,17}`, and `663 = 3*13*17`. Chasing *why* two
incompatible normalisations of the same geometry leave the same fingerprint led
here.

## 2. The generator (the real result)

GGHV22 prints **exactly one** explicit marked polynomial (tex line 1594, the
`(66,99)` / `A_0=(9,24)` case); the `(7,21)` passage at line 2056 has none
(`f1 = y^2/2`, no quartic factor). Our home case has a second, in `AUDIT.md`
item 4. Both come from one forcing ODE, with `C(y) = y^(a-1)(y+1)` and
`f1 := C^3 F_{-v}`:

```
2c * C * f1'  -  2s * C' * f1  =  C^2 ,      s = c + e + 1
```

* `(a,c,e) = (9,3,1)` → `6*C*f' - 10*C'*f = C^2` — GGHV tex 1588 exactly
* `(a,c,e) = (8,4,2)` → `8*C*f' - 14*C'*f = C^2` — `AUDIT.md` item 4 exactly

`marked_polynomial.py` solves it for arbitrary `(a,c)` and **reproduces both
anchors byte-exact** — which is what licenses using it as a generator.

Structure recovered from it:

* a polynomial solution exists **iff `c | a`**, and is then unique (dim 0);
* the solution always factors as `f1 = y^a (y+1)^2 g`, so `f = f1/C = y(y+1)g`;
* therefore GGHV's "deg f = 6, separable, `y(y+1) | f`" **forces `deg g = 4`**;
* with `n := deg g = (c-1)(a/c - 1) + (c-3)`, the equation `n = 4` has **exactly
  four** integer solutions:

| `(a,c)` | marked quartic `g` | `disc` factorisation | sqfree | Galois | `C` |
|---|---|---|---:|---|---:|
| (12,2) | `128y^4 - 64y^3 + 48y^2 - 40y + 35` | `2^24*3^6*5^2*7^2` | **1** | **A4** | 9 |
| (9,3) | `243y^4 - 81y^3 + 54y^2 - 42y + 35` | `3^16*5^2*7^2*13^3` | **13** | S4 | 13 |
| (8,4) | `2048y^4 - 512y^3 + 320y^2 - 240y + 195` | `2^36*3^2*5^2*13^2*17^3` | **17** | S4 | 17 |
| (7,7) | `2401y^4 - 343y^3 + 196y^2 - 140y + 110` | `2^2*5^2*7^12*11^2*29^3` | **29** | S4 | 29 |

Two are GGHV's published cases; **(12,2) and (7,7) are new**. All four have zero
real roots.

## 2.5 THE LAW IS A THEOREM (2026-07-24) — exact discriminant identity

The empirical law of §3 is not empirical. It follows from the ODE in closed form.
Every step below was verified symbolically on all generated cases.

**Reduction.** With `a = cm` and `f1 = y^a (y+1)^2 g`, the forcing ODE
`2c*C*f1' - 2(2c-1)*C'*f1 = C^2` (with `C = y^(a-1)(y+1)`) reduces to

```
c * y(y+1) * g'  -  (c*n*y + c*n + 1) * g  =  1/2 ,        n = deg g
```

**Boundary values** (exact, verified): `g(0) = -1/(2(cn+1))` and `g(-1) = -1/2`.
Evaluating the reduced ODE at a root `alpha_i` of `g` gives
`c*alpha_i*(alpha_i+1)*g'(alpha_i) = 1/2`. Multiplying over the roots and using
the standard `disc = (-1)^(n(n-1)/2) * L^(-1) * prod g'(alpha_i)` formula yields:

> **PROPOSITION.** Let `g` be the ODE-normalised marked polynomial with
> `n = deg g` and leading coefficient `L`, and set
> `Delta := c*n + 1 = (a-1)(c-1) - c`. Then
>
> ```
> disc(g)  =  (-1)^(n(n-1)/2) * 2^(2-n) * Delta * (L/c)^n
> ```

**Verified exactly** on `(a,c) = (12,2), (9,3), (8,4), (7,7), (9,9)` [even `n`]
and `(10,2), (6,6)` [odd `n`] — the reduced ODE, both boundary values, and the
identity itself, all with zero residual.

**Both empirical observations are corollaries:**

* **`n` EVEN:** `2^(2-n)` and `(L/c)^n` are rational squares and
  `(-1)^(n(n-1)/2) = (-1)^(n/2)`, so in `Q^x / (Q^x)^2`
  `[disc(g)] = [(-1)^(n/2) * Delta]` — the law of §3, now proved rather than
  fitted, for the whole modelled ODE family rather than 27 sample points.
* **`n` ODD:** the square class retains the normalised leading coefficient,
  `[disc(g)] = [(-1)^(n(n-1)/2) * 2 * Delta * L/c]`. So **no `Delta`-only law can
  exist at odd `n`** — our 9/9 "counterexamples" were not failures but the
  theorem's other branch.

Note `Delta = cn+1` makes the corner integer a trivial function of the ODE
parameters, which is consistent with §5's finding that it is GGV's published
`bl - a` rather than a deep invariant.

**Scope, unchanged:** this is a theorem about the marked-polynomial ODE family.
Only geometrically grounded cases inherit it. The ODE and `Delta` are prior art
(§ prior-art header); the discriminant consequence is the prospective new
content — and it must still be checked against classical discriminant formulas
for hypergeometric / Jacobi-type families, since the ODE's shape strongly
suggests an equivalent identity may exist there. **Filaseta-Moy's truncated
binomial discriminant is the first place to look, and may well subsume it.**

## 3. The law — refuted, then sharpened

**Original conjecture (REFUTED):** "squarefree part of `disc` = `C`". It fails
whenever `C` is not squarefree — decisively at `C = 9`, where the squarefree part
is **1** — and again at `C = 45` (sqfree part 5, not 45, not 3, not 15).

**Surviving form.** With `C := (a-1)(c-1) - c` and `n = deg g`:

> **If `n` is EVEN:** `disc(g) = (-1)^(n/2) * sqfree(C) * (perfect square)`,
> i.e. `Q(sqrt(disc g)) = Q(sqrt((-1)^(n/2) * C))`.
>
> **If `n` is ODD the statement is false**, and no substitute in terms of `C` is
> visible.

For the four GGHV-shaped cases `n = 4`, the sign is `+1`, so the quadratic
resolvent of the marked quartic is exactly `Q(sqrt C)`.

**Evidence:** 27/27 even-`n` cases (`c` in 2..8, `a <= 32`), including **7 blind
sign-and-value predictions made before computing** — `(16,2)→-13`, `(20,2)→+17`,
`(18,3)→-31`, `(21,3)→+37`, `(16,4)→-41`, `(20,5)→-71`, `(14,7)→-71`, all exact.
The `C = 9` case is the strongest single datum: the law *demands* a square
discriminant there (Galois inside A4), and the quartic really is A4 — a
prediction that could trivially have failed.

**Counterexamples that delimit it:** 9/9 odd-`n` cases break, and not narrowly.
`(12,3)` and `(6,6)` share `C = 19`; the even-`n` one gives `-19`, the odd-`n`
one `-546`.

**It is `C`, not `a_0`, `c`, or degree-4-ness.** Collision tests with the same
`C` but different geometry agree up to the parity sign: `C=17` at both `(8,4)`
n=4 and `(20,2)` n=8 → both `+17`; `C=71` at `(20,5)` n=14 and `(14,7)` n=10 →
both `-71`.

**Polynomial vs field discriminant is a non-issue** — provably, not empirically:
for primitive `f` of degree `n`, the monic model satisfies
`disc(F) = a_n^((n-1)(n-2)) disc(f)` with `(n-1)(n-2)` always even, and
`disc(F) = [O_K : Z[alpha]]^2 d_K`, so `sqfree(disc f) = sqfree(d_K)` identically.

## 4. The tie to the q_window theorem

`deg g = 4` ⟺ `M = C` at base pair `(2,3)` ⟺ **`q_window = 1`**. Both conditions
reduce to the same four `(a,c)`. So the GGHV-shaped cases are *exactly* the
integral ones — which connects this directly to the 51-integral-families result.

The law is **not** confined to the integral locus, though: `(a,c) = (6,3)` is
`(t,kappa,q) = (3,1,5)` with `q_window = 13` (non-integral), marked polynomial
`9y^2 - 3y + 2`, `disc = -63`, squarefree part `-7`, and `|-7| = 7 = C`.

## 5. THE WEAK LINK — RESOLVED, AND THE DICTIONARY WAS PARTLY WRONG

**Update (2026-07-24, same day).** The dictionary was checked against GGHV's
actual case table and the three worked reductions. Verdict: **two of three
components confirmed, one refuted.** A third instantiation of the forcing ODE was
found — the `(7,21)` / `(56,84)` case — and it decides the question.

| component | verdict | evidence |
|---|---|---|
| `t = c` | **CONFIRMED** 3/3 | `t` is the exponent `l` of the final Laurent chart `phi(x)=x^-1, phi(y)=x^l y`; GGHV states it per case (tex 655, 1229, 1395) |
| `kappa = t - 2` | **CONFIRMED — a theorem** | `[x^-1, x^l y] = -x^(l-2)`, so the bracket exponent is forced; GGHV's three reductions state brackets `x`, `x^2`, `x` |
| `q = a_0 - 1` | **REFUTED** | the true content is `q = ord_y C`, `a = deg C`. `q = a-1` holds iff `dg := deg C - ord_y C = 1` |
| `a = A_0[0]` | **REFUTED** | a two-anchor coincidence. For `(7,21)`, `A_0[0] = 7` but `deg C = 1` |

**The deciding case.** `(7,21)` reduces with `C = y` (a monomial — `A_0` lands on
the y-axis and `ell_{1,0}` degenerates), `c = 3`, `s = 5`, `f1 = y^2/2`. Verified:
the forcing ODE residual is **exactly 0**, so it is a genuine third instantiation.
But `q = ord_y C = 1`, not `a_0 - 1 = 6`; and `deg C = 1`, not 7.

**Corrected dictionary:**

```
t      = l, exponent of the final Laurent chart phi(x)=x^-1, phi(y)=x^l y
         (= c when m = 2; in general 2c = m*t)
kappa  = t - 2                       [PROVEN, not fitted]
q      = ord_y C                     [NOT a_0 - 1]
a      = deg C                       [NOT A_0's first coordinate]
s      = t(n-m) + kappa + 1
C_int  = q(kappa+1) - t
```

**Necessary precondition** before applying any of it (derived from the three
worked reductions, not read from the paper): **`b_0 = t*(a_0 - 1)`**, i.e.
`t = b_0/(a_0-1)` must be an integer. Verified against GGHV's table: it gives
`t = 3` for `(9,24)` and `t = 4` for `(8,28)` — both matching the read values —
and correctly flags `(7,21)` as `7/2`, non-integral, standard chart inapplicable.
It also resolves the census ambiguity of §5-old: for `a_0 = 8`, `b_0 = 28` forces
`t = 4`, `q = 7`, hence `(4,2,7)` and `C_int = 17`; the competing rows
`(4,2,2)`/`(8,6,2)` come from `(8,24)`, where `24/7` is not an integer.

**Consequence for §3.** The arithmetic law is untouched — it is an identity about
`(a-1)(c-1) - c` within the class the generator hard-codes (`C = y^(a-1)(y+1)`,
i.e. `dg = 1`, and `m = 2`). What is now correct is its SCOPE: it applies where
`dg = 1` and `m = 2`, which covers both anchors and, in the census, F1, F5, F17
and the `(8,28)` row. It must **not** be evaluated from `A_0` in general.

## 5.1 A separate finding: the census `t` may not be the chart exponent

`chain_survey.py:447` sets `t = final[1]` (`l_final`, the denominator of GGV5's
final corner `A_final = (p\l, b)`). GGHV line 433 describes that `l` as a
ramification/conjugacy index from GGV1 Cor 7.4 — *not* as the Laurent-chart
exponent. For `(7,21)` GGV5 gives `l = 7` while GGHV's chart exponent is 3.

The `b_0 = l(a_0-1)` criterion passes for F1–F8 and F14–F17 and the `(8,28)` row,
and **fails for F9, F10, F11, F12, F13** (F9/F10/F11: `21/6 = 7/2`; F12: `24/7`;
F13: `21/8`) plus the `(8,24)`-rooted F22/F23/F24 that `phi_corner4.py` already
flags. **F9 is a LANDED point** (`family_grammar.py:89`,
`("F9",0) -> "(56,84)"`), and F9/F10/F11 are not currently flagged.

### ESTABLISHED (2026-07-24, from the GGV1 arXiv source)

The two quantities are **definitively different objects**, and GGV1 proves it
inside a single section:

* **The chart exponent is an EDGE SLOPE.** GGV1 Remark 8.6 gives
  `psi_l_bar(rho,sigma) = (-rho, l*rho + sigma)`, so `psi_l_bar(-1,l) = (1,0)`:
  `l` is the slope `dx/dy` of the Newton-polygon edge being flattened to vertical.
  Verified on all three GGHV cases: edge `{(24,9),(21,8)}` -> 3; `{(28,8),(24,7)}`
  -> 4; `{(21,7),(0,0)}` -> 3. All match GGHV's stated `l`.
* **The bracket is forced:** `[phi(x),phi(y)] = [x^-1, x^l y] = -x^(l-2)`
  (verified symbolically), hence `kappa = l - 2` — which is why `kappa = t-2` was
  a theorem rather than a fit.
* **GGV5's final-corner `l` is the Laurent / ramification index** of
  `L^(l) = K[x^(1/l), x^(-1/l), y]` — a Puiseux denominator, the number of
  conjugate branches. GGV1 §1 defines it that way, and every §7 corner is written
  `(a/l, b)` for that reason.
* **The same-paper proof that they differ:** GGV1 §8 works entirely in `L^(1)`
  (ramification index **1**) while charting with exponent **3**.

So `chain_survey.py:447`'s `t := l_final` is measuring a ramification index, and
it is **the wrong quantity for chart purposes**. The `(7,21)` discrepancy (GGV5
`l = 7` vs GGHV chart `3`) is exactly this.

Also worth recording as a naming hazard: **four distinct objects share the
letters `l` and `q`** across these papers — GGV1's Laurent index `l`, GGV1's power
denominator `q_j` (Thm 7.6(3): the exponent making the leading form a `q*m`-th
power), GGHV's chart exponent, and our `ord_y C`. Our code uses `t` and `q` for
two of them.

**`q` and `a`, correctly characterised:** `a = deg_y C` and `q = ord_y C` are the
`y`-coordinates of the two corners bounding the `(-1,l)`-edge, which sit at
`(l(a-1), a)` and `(l(q-1), q)`. Generically `q = a-1` (edge form
`y(x^l y - alpha)^(a-1)`, cut to `C = y^(a-1)(y+alpha)`); the degenerate
`q = a = 1` gives `C = y`, our `(7,21)` case.

**And the precondition is in the literature.** `b_0 = l(a_0 - 1)` is the
hypothesis of **GGV1 Proposition 7.8(3)** in transposed coordinates
(`A_0 = (1,0) + r(1,rho_0)` with `A_0' = (1,0)`). The general form is the
edge-slope identity `b_0 - b_prev = l(a_0 - a_prev)`; our version is the branch
where the edge runs from `(0,1)`, and `b_0 = l*a_0` is the branch where it runs
from `(0,0)` — which is precisely `(7,21)` (`21 = 3*7`). Verified numerically.

**NO GENERAL DICTIONARY EXISTS IN THE LITERATURE.** GGV1 §8 is hard-coded to one
polygon, and its closing **Remark 8.13** says so explicitly — *"But in cases
different from `16m,16n` this cannot always be obtained"* — and exhibits a
counterexample where the required squeezing is unachievable. GGHV's ad-hoc
per-case treatment is not laziness; it reflects a real gap. So a general
`(A_0, A_1, m, n) -> (t, kappa, q, deg C)` rule is **not** available to be looked
up, and anything we want of that kind we would have to prove.

### What remains open

Our `t` is *defined* as `l_final` and the
q_window identity `t*H - q*M = q(kappa+1) - t` is proved symbolically, so it is
true whatever `(t,kappa,q)` denote. The open question is whether `M`, `H` computed
from `l_final` are the correct window numerator/denominator for the F9–F13
families, or whether those rows need the chart exponent instead. That is a
question about the *geometric meaning* of `M` and `H`, and it is not settled here.
`family_grammar_verify.py` (210 checks) and the rest of the suite pass either way,
since they do not depend on this identification.

**Next step to settle it:** GGV1 (arXiv:1401.1784) §8 — GGHV line 460 says its
reductions are "reminiscent of the procedure in section 8 of the ArXiv version of
[GGV1]" and cites GGV1 Prop 8.2. That is the only place a general
`(A_0, A_1, m, n) -> (t, kappa, q, deg C)` rule could be stated.

## 5.2 The original weak-link note (superseded by §5 above)

The dictionary `(t, kappa, q) = (c, c-2, a_0 - 1)` is **fitted from two anchors**,
not read off GGHV's Newton-polygon data. And census matching on `a_0` + base pair
does **not** pin it:

* `a_0=8`, base {2,3}: `(t,kappa,q)` ∈ {(4,2,2), (4,2,7), (8,6,2)} → `C` ∈ {2, **17**, 6}
* `a_0=9`, base {2,3}: {(3,1,8), (9,7,4)} → `C` ∈ {**13**, 23}
* `a_0=7`, base {2,3}: {(7,5,2), (7,5,5), (7,5,6)} → `C` ∈ {5, 23, **29**}
* the `C=9` case `(t,kappa,q) = (2,0,11)` has **zero** census rows

So: the arithmetic law about `(a-1)(c-1) - c` is solid and machine-checked. The
claim that this quantity **is the census corner integer** rests on the fitted
dictionary. **Pinning `(t,kappa,q)` directly from GGHV's Newton-polygon data for
a third case is the highest-value next check** — a paper-reading job, not a
computation. If the dictionary is wrong, the arithmetic law survives and the
"corner" *interpretation* collapses.

**Second scope caveat:** the ODE for arbitrary `(a,c,e)` is an *interpolation*
between GGHV's two instantiated ODEs. It reproduces both anchors exactly, but the
intermediate `(a,c)` pairs are not GGHV-certified Newton-polygon geometries. The
27 confirmations confirm an arithmetic law about that ODE family; only 2 are
certified to correspond to real Prop-4.3 cases.

## 6. What would falsify it

1. Any even-`n` case with `|sqfree(disc g)| != sqfree(|C|)` — none in 27 tries.
2. A sign violation at even `n` — none found.
3. A demonstration that the dictionary is wrong (§5). Highest value.
4. A GGHV-realised case with odd `deg g` — none exists in this family, since
   `deg f = 6` forces `n = 4`.

Checker: `marked_polynomial.py` (anchors byte-exact + the four deg-4 cases).
