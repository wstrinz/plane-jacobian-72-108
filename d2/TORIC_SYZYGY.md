# The toric syzygy `6*W*Z = e^5`, and the structure of the bare G-variety

**Checker:** `toric_syzygy.py` — **67/67**, `--quiet` exit 0.
**Paper edits landed:** `PROOF_72_108.md` §3.5 (new), §4.2a (new), §4.3 ledger, §8.3, §8.4, §10, §13.3, §14.11.
**Generators:** loaded from `bigrade_annotator._G_generators()`. Nothing retyped.

Notation throughout: `W := e*S - R^2`, `Z := e*T - R*S` (the first two 2x2 minors of the
Hankel matrix `[[e,R,S],[R,S,T]]`); `t = y+1`; `Phi` is the fixed degree-34 polynomial of
§2.2, `2*Phi = 2*c*t^30*q` with `q` squarefree of degree 4 and `q(-1) != 0`.

Where §8's normalised coordinates appear, the paper's `Z = A^2 - gamma*Pi^2*v` is written
`Z_paper` and the minor is written `Zcal`, to keep them apart.

---

## 1. The headline results

### R1. The bare variety is nonempty — and completely describable **[PROVED]**

On the locus `e != 0`, `W != 0`, the equations `G1 = G2 = G3 = 0` determine `(d0,d1,T)`
uniquely from the free data `(d2,e,R,S)`:

```
T  =  e^4/(6*W)  +  R*S/e
d1 = -2*d2*R/e  -  4*R*S/e^2  -  e^3/(3*W)
d0 =  d2*R^2/e^2  +  S*(e*S + 2*R^2)/e^3  +  e^2*R/(3*W)
```

and the remaining equation `G5 = 0` becomes **exactly the K-syzygy of Theorem 3.1**:

```
2*Phi  =  e*(d2*e^2 + 3*e*S + 3*R^2)
```

By Corollary 3.6 (below) that locus is not a chart — it is the whole variety. So
`V(G1,G2,G3,G5)` is, in its entirety, the graph over one hypersurface of three explicit
rational functions. (`toric_syzygy.py` E1, E2.)

### R2. The universal toric identity **[PROVED]**

```
2*e^2*G3  -  4*e*R*G2  +  2*R^2*G1   ==   6*(e*S - R^2)*(e*T - R*S)  -  e^5
```

residual exactly 0 in `Z[d0,d1,d2,e,R,S,T]`. Hence **every** G-point satisfies

```
6*W*Z = e^5
```

using no cap, no slice divisibility, no branch condition, no support enumeration, no degree
stratum, and **not even `Phi`** — `G5` is the only generator carrying `Phi` and it is not
among the three rows used. (B1.)

`T`'s formula in R1 says `e*T - R*S = e^5/(6*W)`. **The toric syzygy is the `T`-formula of
the chart.** That is why it costs nothing: solving `G1` for `T` *is* the syzygy. (E3.)

### R3. `W` never vanishes, so the chart is total **[PROVED]**

**Corollary 3.6.** On every G-point, `e != 0` and `W != 0`.
Proof: `e != 0` is Lemma 3.3; if `W == 0` then `e^5 = 6*W*Z = 0`. (B4.)

### R4. Theorem 8.1 is R2 in normalised coordinates **[PROVED]**

Under the §8 ansatz `e = gamma*t^9*Pi`, `R = t^9*A`, `S = t^9*Pi*v`, `T = t^9*C`, with `C`
eliminated by `g1 = 0`:

```
W    = -t^18 * Z_paper
Zcal = -t^18 * Pi * F
```

so `6*W*Zcal = e^5` reads `6*t^36*Pi*F*Z_paper = gamma^5*t^45*Pi^5`, i.e. after cancelling
`6*t^36*Pi`,

```
F * Z_paper = (1/6) * gamma^5 * t^9 * Pi^4          <- Theorem 8.1, exactly
```

The identity holds with the exponent `9` replaced by a **free symbol** `a`, so `Pi^4` is
independent of `a_t = 9`; only the `t`-power moves. **The cofactor identity was never a
lucky discovery downstream of `a_t = 9`** — it is the normalised shadow of a relation
present on the whole variety. (C1–C4.)

### R5. The `Pi^4` is a forced local contact order **[PROVED, modulo Lemma 11.5 — re-proved here]**

`e | S` (Lemma 11.5) gives `Z = e*Z1` with `Z1 = T - R*Sbar`, so `6*W*Z = e^5` becomes

```
6*W*Z1 = e^4      hence      v_beta(W) + v_beta(Z1) = 4*v_beta(e)   at every place
```

At a root `beta` of `Pi`:

| input | value | source |
|---|---|---|
| `v_beta(e)` | `1` | `e = gamma*t^9*Pi`, `Pi \| q` squarefree, `t(beta) != 0` since `q(-1) = 3315 != 0` |
| `v_beta(W)` | `0` | `e(beta) = 0` so `W(beta) = -R(beta)^2`, and `R(beta) != 0` is §8.2 |
| `v_beta(Z1)` | **`4`** | the divisor law |

Since `Z1 = -t^9*F/gamma`, this is `v_beta(F) = 4`, i.e. `Pi^4 | F`. With
`gcd(Z_paper, Pi) = 1` (from `Z_paper(beta) = A(beta)^2 != 0`), **all** of the `Pi^4` in
`(*)` sits in `F` and none in `Z_paper`. The exponent `4` traces to the `e^5` of R2, shed of
one factor of `e` by `e | S`. A non-squarefree `Pi` would force order `8`; `W(beta) = 0`
would force `3`. (D0–D6.)

---

## 2. The two explicit families of bare G-points

Both are genuine G-points in the sense of **Definition 4.1**, i.e. with the *fixed* `Phi`.

### Family (a) — `R = 0`

```
d0 = 1,  d1 = -E/3,  d2 = 2*Phi/E^3 - 3,  e = E,  R = 0,  S = E,  T = E^2/6
```

for any nonzero `E` with `E^3 | 2*Phi`. Verified with `E` **symbolic** (all four residuals
exactly 0), and the divisibility hypothesis is satisfied **exactly** by

```
E = lam * t^m,   lam != 0,   0 <= m <= 10
```

because any irreducible `p | E` with `p` not `t` would divide `q`, giving `v_p(2*Phi) = 1 < 3`.
(A1, A2.)

### Family (b) — `W` a unit, `R` arbitrary

For any `lam, w != 0` and **any** polynomial `A`:

```
e = lam,   R = A,   S = (A^2 + w)/lam,   d2 = ( 2*Phi/lam - 6*A^2 - 3*w ) / lam^2
```

with `(d0,d1,T)` from R1. Then `W = w` and `e = lam` are both units, so every denominator in
the chart is a unit and all seven coordinates are polynomials. This family contains a **free
polynomial parameter of unbounded degree**, and `R != 0` generically. Instances verified:
`A = y`, `1`, `y^2-3`, `y^3+y`, `y/2+1`. (E4.)

Family (b) is the important one. It shows there is no clever strengthening of the four
equations to be found, because on the (total) chart they are *equivalent to one equation*.

### Which admissibility hypothesis excludes them

**(A4), the cap profile of Lemma 2.5, alone — but in three cases, and the naive answer is
wrong.**

| | config (1), `lambda = 3` | config (2), `lambda = 2` |
|---|---|---|
| family (a), `m <= 9` | `deg d2 = 34 - 3m >= 7 > 6` | `>= 7 > 4` |
| family (a), `m = 10` | `deg d1 = 10 > 9` | `deg d1 = 10 > 6`; also `deg T = 20 > 16` |
| family (b), any `A` | `deg d2 >= 16 > 6` | `deg d2 >= 16 > 4` |

⚠️ **The review's framing that `d2` alone does it is FALSE at `m = 10`.** There
`deg d2 = 4`, comfortably inside both cap profiles; the kill is `deg d1 = 10 > 9`.
(A3 checks exactly this, and asserts the *negation* of the `d2`-only claim.)

The `deg d2 >= 16` bound in family (b) is sharp and is not the obvious `34`: if
`deg A != 17` then `deg(2*Phi/lam - 6*A^2) >= 34`, while `deg A = 17` permits cancellation —
matching the top eighteen coefficients is a triangular system with two sign branches, and
both leave degree exactly `16`. So (A4) would still exclude family (b) with a cap loosened
to `15`, but not to `16`. (E5.)

**Non-vacuity:** raising `lambda` to 4 lets family (a)'s `m = 9, 10` through. The exclusion
is a real consequence of the *values* `lambda = 3, 2` that Lemma 2.5 derives, not of the
shape of the test. (MUT A3.)

---

## 3. A defect found and repaired in the paper

`PROOF_72_108.md` §4.2 (old line 349) asserted the emptiness of `V(G1,G2,G3,G5)` was **open**
while, in the same sentence, pointing at "an explicit point of it" recorded in §10 — a
self-contradiction either way round. §14 item 11 repeated the "open" reading.

Investigating the cited point turned up a second, subtler issue. The witness is
`t1_branch.py` `C12_negatives`:

```
e = 1, R = 0, S = 1, T = 1/6, d0 = 1, d1 = -1/3, d2 = 0, Phi = 3/2
```

That `Phi = 3/2` is **solved for** — the checker takes `Phi` as a free indeterminate and
reads off the value that makes `G5` vanish. In the ambient ring
`Q[d0,d1,d2,e,R,S,T,Phi]` this is correct and sufficient for §10's purpose, which is an
ideal **non**-membership. But it is **not a G-point in the sense of Definition 4.1**, which
fixes `Phi` to the degree-34 polynomial of §2.2: substituting that tuple with the fixed
`Phi` leaves `G5 = Phi - 3/2 != 0`. So §4.2's parenthetical was citing, as a point of the
bare variety, an object that is not one.

Family (a) at `(lam,m) = (1,0)` is the repair. It is the same shape but with
`d2 = 2*Phi - 3` instead of `d2 = 0`. (A4 checks all of this, including that the old tuple
fails with the fixed `Phi`.)

**Landed edits.** §4.2's parenthetical replaced; new §4.2a with Proposition 4.4,
Corollary 4.5 and three remarks; §14 item 11 struck and replaced with the settled statement;
§10's bullet carries a warning about which ambient ring its witness lives in; §4.3's ledger
extended; §3.5 added for Theorem 3.5 / Corollary 3.6; §8.3 and §8.4 reframed per R4 and R5.

---

## 4. PROVED / CHECKED / INFERRED

**PROVED** — exact symbolic identity or finite exhaustive enumeration, hand-checkable in
principle, mutation-controlled:

* Theorem 3.5, `2e^2 G3 - 4eR G2 + 2R^2 G1 = 6WZ - e^5`. Hand-checkable (~40 terms).
* Corollary 3.6, `W != 0` on every G-point. One line, given Lemma 3.3.
* §10's companion `W^2 = R^4 + d2 e^2 R^2 + d1 e^3 R + d0 e^4 + (2/3)(e^2 G2 - eR G1)`.
* Proposition 4.4, the chart, and that `G5` on it is exactly the K-syzygy.
* Corollary 4.5(a), including the classification `E = lam*t^m`, `0 <= m <= 10`.
* Corollary 4.5(b), for the five instances exhibited; the *shape* is proved for all `A`
  since every denominator is a unit by construction.
* The cap exclusion table, including `deg d2 >= 16` on family (b).
* The specialisation `W = -t^18 Z_paper`, `Zcal = -t^18 Pi F`, and the collapse of
  `6WZ = e^5` onto Theorem 8.1; and the free-exponent form.
* Lemma 11.5 (`e | S`), re-proved here by the Sylvester/adjugate route rather than imported:
  `Res_R(A2,A3) = -2e[S^7 + sum e^i alpha_i S^(7-i)]`, every coefficient divisibility exact.
* The divisor law `6*W*Z1 = e^4` and `v_p(W) + v_p(Z1) = 4*v_p(e)`.
* The two UFD steps, by exhaustive enumeration over 780 and 150 factorisations respectively,
  each with its hypothesis-dropping control.

**CHECKED** — verified computationally here, but resting on a result imported from elsewhere
in the paper:

* `R(beta) != 0` at marked roots (§8.2). Re-derived here as `box mod Pi = 3A^2 - mu t^3 Q_Pi`,
  but the conclusion that all three factors are nonzero is §8.2's.
* `v_beta(e) = 1` at marked roots. Rests on Theorem 3.4's `Pi | q` squarefree.
* The §8.1 reduction (rebuilt here, residual 0, with the `t^28` falsifiability control), which
  the contact-order argument uses to identify `Z1 = -t^9 F/gamma`.
* Lemma 3.3, used in Corollary 3.6. Imported from §3.3; not re-proved here.

**INFERRED** — stated in the paper on the strength of the above, but not itself a machine
check:

* "The `Pi^4` is a forced contact order **rather than** an artifact of the `t^9`
  normalisation." The mathematical content (R5) is proved; the word "rather" is an
  interpretive claim about where the exponent comes from. What is *proved* is that the
  identity survives a free exponent `a`, which is the strongest form of the claim available.
* "The `G`-point / admissible-germ distinction is mathematically essential." Follows from
  Corollary 4.5(b) plus the observation that no (A1)-only argument can exclude a family with
  a free unbounded parameter — but "essential" is a judgement, not a theorem.
* The minimal cap profile that still excludes both families is pinned between 15 and 16 for
  family (b) only; the corresponding question for family (a) was not swept.

---

## 5. Mutation controls

Every positive check in `toric_syzygy.py` is paired with a control that **must** fail; the
ledger reports a vacuous check as a failure. Notable ones:

| control | what it kills |
|---|---|
| nine corruptions of family (a)'s seven coordinates | all leave a nonzero generator |
| seven coefficient corruptions of `(2,-4,2,6,1)` | all leave a nonzero residual |
| exponents `e^3, e^4, e^6, e^7` | all fail; the `5` is forced |
| `W <-> Z` swap; `6W^2 - e^5`, `6Z^2 - e^5` | both wrong; the two minors are not interchangeable |
| dividing `K` by `t^28*Pi` instead of `t^27*Pi` | not a polynomial (the paper's own falsifiability remark) |
| dropping either minus sign in `W = -t^18 Z_paper`, `Zcal = -t^18 Pi F` | fails |
| `t^35` in the collapse factor; dropping `Pi` from it | both fail |
| `1/5` in place of `1/6` in Theorem 8.1 and in the `T`-formula | fails |
| dropping `gcd(Z,Pi) = 1` | 776 of 780 factorisations violate `Pi^4 \| F` |
| dropping `W(beta) != 0` | 120 of 150 factorisations violate `v_beta(Z1) = 4` |
| raising `lambda` to 4, 5, 6 in the cap test | 8 `(lambda,m)` pairs survive — the test is not vacuous |
| `e^(i+1)` in the `e \| S` integrality pattern | fails; the pattern is sharp |
| `R = 0` specialisation of `6WZ - e^5` | does *not* give `-e^5`, so B4 discriminates |

**One control that legitimately does not fire, recorded rather than hidden.** On the `R = 0`
family the syzygy degenerates to `2*e^2*G3 = 6WZ - e^5`, because the `G1` and `G2` cofactors
both carry factors of `R`. Since `G3|_{R=0} = -e^3/2 + 3*S*T` contains neither `d0` nor `d1`,
a `d0`- or `d1`-perturbation of that family provably **cannot** break the law. This is a
structural fact about the slice, not a missing control; the checker asserts it explicitly
(D2) and uses `S`-, `T`-, `e`-perturbations plus four generic off-variety tuples instead.

---

## 6. What did not reproduce

Nothing failed. Two claims in the incoming review were **sharpened rather than confirmed**:

1. *"Most visibly `d2 = 2*Phi/E^3 - 3` catastrophically violates the degree cap of
   Lemma 2.5."* True for `m <= 9`, **false for `m = 10`**, where `deg d2 = 4` is inside both
   cap profiles and the kill is `deg d1 = 10 > 9`. The review flagged this as a claim to
   check rather than assume; it needed the check.
2. *"§10 records an explicit point of it."* The recorded point has `Phi` free and solved to
   `3/2`; it is not a Definition-4.1 G-point. See §3 above.

And one result came out **stronger** than the review proposed: the family in the review has
`R = 0`, but the chart (R1) yields family (b), with `R` a free polynomial of unbounded
degree, and Corollary 3.6 shows the chart covers the entire variety — so the bare variety is
not just nonempty, it is a completely described rational hypersurface.
