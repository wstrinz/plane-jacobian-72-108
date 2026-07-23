# T5 coefficient-field audit and split-place repair

**Date:** 2026-07-21
**Repository head inspected:** `75f8c50dac4b6c18b75a9ff53ef7cf1fb231c334`
**Scope:** `d2_plane_72_108`, especially the exact `f31`/`f37` T5 campaign.

> **Post-audit progress (2026-07-22).** `split_place_ledger.py` now generates
> the full geometric ledger and proves 81 strata terminally impossible.
> `F37_FREE_FAMILY_SYSTEM.md` restores three small pre-resultant equations and
> proves that no point of the bare `f37` free family lifts to the original
> system. The recommendation in item 5 below has therefore been completed.

## Executive finding

The recent exact campaign contains a genuinely powerful idea — the joint
`(y+1)`/quartic-place cascade — but its current ledger treats

```text
q = 2048 y^4 - 512 y^3 + 320 y^2 - 240 y + 195
```

as one irreducible prime of degree four because it is irreducible over
`Q`.  The ambient Jacobian problem, however, is over an arbitrary
characteristic-zero field `K` (and in particular over `C`), while the numerical
harness explicitly searches complex coefficient windows.  After base change to
an algebraic closure, `q` is a product of four distinct linear factors.  A
solution may contain some of those factors in `e=d_{-1}` without containing the
whole polynomial `q`.

Consequently:

* the scalar stratum label `a_q=v_q(e)` is not field-stable;
* the survivor counts `21 -> 15 -> 13 -> 12 -> ...` are not yet counts for the
  geometric problem;
* any step using “`q` is prime”, “one `q` costs four degrees”, or the `S_4`
  quartic-field argument needs either a descent theorem or a split-place repair.

This does **not** erase the program.  The `t=y+1` cascade, graded decompositions,
window audits, exact identities, numerical evidence, the `f37` free-family
finding, and many coprime/divisibility arguments remain useful.  It changes the
correct case space and the order in which it should be attacked.

The right replacement is to base-change once and write

```text
q = p1 p2 p3 p4,
bi = v_pi(e / t^a),
```

with the multiplicity vector `(b1,b2,b3,b4)` considered up to the `S4` symmetry
of the four roots.  The degree cap becomes

```text
a + b1 + b2 + b3 + b4 + (mass away from t and q) <= 10.
```

The existing `a_q=0,1,2` rows are only the uniform subfamilies
`(0,0,0,0)`, `(1,1,1,1)`, `(2,2,2,2)` with a cofactor geometrically coprime to
`q`; partial-support vectors are additional branches.

## What remains valid without change

1. **The `t`-cascade.**  Lemma A comes from powers of `t=y+1`, the degree windows,
   and polynomial division by `t^v`; it is field-stable.
2. **The graded identities.**  The decompositions of `f31` and `f37`, the small
   cofactors `h_f`, and their symbolic verification are exact over `Z`.
3. **The upstream audit.**  Reproduction, equation selection, window
   transcription, and the master dichotomy `f31=0 or f37=0` are not affected by
   factorization of `q`.
4. **The `f37` free family.**  The identity
   `f37|_{d1=d2=0}=0` is an integral polynomial identity.  It remains the decisive
   reason that the bare `f37` resultant cannot close its branch.
5. **The `a=10` `d1=0` result.**  Here `e=C t^10`, so `gcd(e,q)=1` geometrically;
   division by the whole polynomial `q` is legitimate even when `q` splits.
   The unit-ideal certificates are over `Q`, hence remain certificates after
   field extension.
6. **Any argument explicitly assuming `gcd(e/t^a,q)=1`.**  If every geometric
   factor of `q` is coprime to the cofactor of `e`, divisibility by the whole
   squarefree polynomial `q` follows factor by factor.  Such proofs should be
   relabeled “geometrically `q`-coprime”, not `a_q=0`.

## Results requiring re-audit

The following mechanisms are not automatically geometric:

* the 21-row `(a,a_q)` ledger and every survivor count derived from it;
* degree-starvation inequalities of the form `4(7-3a_q)>10+3a` unless the
  relevant full power of `q` has already been proved to divide the object;
* cancellations using `q | XY => q | X or q | Y` without geometric coprimality;
* the `S4` root-field/no-subfield step in the `(6,1)` T1 reduction;
* the original `sigma`-locus proof as written, because it uses a single
  `q`-valuation.

The last item is repairable, as follows.

---

## The split-place sigma-locus theorem

The argument applies to the `f31` master equation

```text
A^4 B = const * e^17,
5 A + B = const * Phi,
4 A - B = const * d2 * e^3,
```

and, with changed nonzero constants, to the `f37` sigma-locus master equation
when `d2 != 0`.

Base-change to an algebraic closure.  The divisor of `Phi` consists of one root
`t=0` of multiplicity `30` and four simple roots of `q`.

### 1. Constant cases

If `A` or `B` is constant, differentiation gives `e^16 | Phi'`.  Since

```text
Phi' = const * t^29 * (30 q + t q')
```

and the second factor has degree four, the same short argument already in the
repository forces `e` to be a scalar multiple of `t`, then evaluation at `t=0`
gives a contradiction.  This part never needed `q` to be irreducible.

### 2. Nonconstant cases force full degrees

From

```text
4 deg A + deg B = 17 deg e,
deg A, deg B <= 34,
deg(const*A + const*B) = deg Phi = 34,
deg(d2 e^3) <= 4 + 3 deg e,
```

one obtains, field-independently,

```text
deg A = deg B = 34,
deg e = 10.
```

### 3. Local valuation options

At a root of `Phi` of order `phi`, put

```text
m = v(e),  alpha = v(A),  beta = v(B).
```

Then

```text
4 alpha + beta = 17 m.
```

If `alpha != beta`, the nonzero linear combination equal to `Phi` has order
`min(alpha,beta)=phi`; if `alpha=beta`, cancellation can only raise the order,
so `alpha<=phi`.

The complete possibilities under the global degree bounds are:

```text
at t (phi=30):
    (m,alpha,beta) = (0,0,0), (5,17,17), (9,30,33)

at a simple q-root (phi=1):
    (m,alpha,beta) = (0,0,0), (1,1,13), (1,4,1),
                     (2,1,30), (5,21,1).
```

Away from the five roots of `Phi`, a root is either

* type A: `(v(e),v(A),v(B))=(4k,17k,0)`, or
* type B: `(v(e),v(A),v(B))=(m,0,17m)`.

Degree balance leaves six labeled patterns: two actual patterns plus four
symmetric copies of one exceptional pattern.

### 4. Mason--Stothers kills the two generic patterns

After dividing `G=gcd(A,B)` from the relation between `A`, `B`, and `Phi`, the
three terms are pairwise coprime.

* If no special root divides `e`, the outside masses are `8` of type A and `2`
  of type B.  The reduced degree is `34`, while the radical has at most
  `2+2+5=9` roots.  Mason gives `34 <= 8`, contradiction.
* If `v_t(e)=5`, then `v_t(A)=v_t(B)=17`; after division the reduced degree is
  `17`, and the radical has at most `1+1+5=7` roots.  Mason gives `17 <= 6`,
  contradiction.

### 5. The sole split exceptional pattern is impossible

Up to permutation of the four roots, the only pattern not killed by the crude
Mason bound is

```text
e = const * t^9 * p,
(v_t(A),v_t(B)) = (30,33),
(v_p(A),v_p(B)) = (4,1),
```

where `p` is one simple factor of `q`.

But the second linear relation has left-hand order `1` at `p`, while
`d2*e^3` has order at least `3`.  Contradiction.

Therefore the `f31` sigma-locus is empty over every characteristic-zero field.
The same proof repairs the `f37` sigma-locus theorem on `d2 != 0`; its
`d1=d2=0` free family remains, as it must.

The finite enumeration is checked by `t5_split_place_verify.py`.

---

## New result: the `a=7` geometrically q-coprime stratum is dead

Consider `f31`, subcase 2, with

```text
e = t^7 E,
deg E <= 3,
gcd(E,tq)=1,
v = 30-3*7 = 9,
deg g_l <= 31.
```

The repaired sigma-locus theorem handles `d1=0, sigma=0`.  The two remaining
branches are impossible.

### T2: `d1=0`, `sigma!=0`

The terminal pair gives

```text
g7 = 0,
E^3 g6 = 3072 c^6 q^6 sigma^2.
```

Geometric coprimality gives `g6=q^6 G`, `deg G<=7`, and

```text
E^3 G = 3072 c^6 sigma^2.                       (1)
```

Using

```text
h5 = -9216 d2 sigma^2 + 2048 t^14 E^2
```

in level 5 and eliminating `sigma^2` with (1) yields exactly

```text
E^3 (g5 + 19890 q^5 d2 G)
    = t^9 q^5 (qG - 2048 c^5 t^5 E^2).          (2)
```

The bracket on the left has degree at most `31`.  Since `gcd(E,tq)=1`, write

```text
qG - 2048 c^5 t^5 E^2 = E^3 S,
```

and (2) forces `deg S<=2`.  Hence

```text
qG = E^2(2048 c^5 t^5 + E S),
2048 c^5 t^5 + E S = q L,
deg L<=1,
G=E^2 L.
```

Substitution in (1) gives `E^5 L = const * sigma^2`, so `E L` is a square up to
a nonzero scalar.

* If `L` is constant, `E` is a scalar square and `deg E<=2`; but
  `qL-2048c^5t^5=ES` has degree `5` on the left and at most `4` on the right.
* If `L` is linear, square parity forces `L|E`; reducing the preceding equation
  modulo `L` gives `L|t^5`, hence `L` is proportional to `t`, contradicting
  `gcd(E,t)=1`.

Thus T2 is dead.

### T1: `d1!=0`

The terminal equation gives

```text
g7=q^7 H,  deg H<=3,
E^3 H = const * d1^2.
```

UFD parity gives

```text
E  = gamma * s * u^2,
H  = eta   * s * v^2,
d1 = delta * s^2 * u^3 * v,
```

where `s` is squarefree and

```text
deg s + 2 deg u <= 3,
deg s + 2 deg v <= 3.
```

After writing `g6=q^6 G6`, level 6 is

```text
E^3 G6 + c^6 h6 = t^9 q H.                     (3)
```

At a root of `s`, parity and (3) rule out every local exponent pair except

```text
v(E)=1, v(H)=3, v(d1)=3.
```

In that sole escape, (3) forces `v(sigma)>=2` and `v(G6)=0`.  Level 5 then has

```text
v(E^3 g5) >= 3,
v(h5) = 2             (the unique lowest term is 2048 e^2),
v(t^9 g6) = 0,
```

which is impossible.  Hence `s` is constant.

Now `deg u,deg v` are each `0` or `1`.  The exact degree table for the eight
terms

```text
Phi^f e^(21-3f) h_f,  f=0,...,7,
```

has a unique largest term in every one of the `2*2*9=36` cases determined by
`deg u`, `deg v`, and `deg sigma` in `0,...,8`.  The top term is always one of
`f=6` or `f=7`; the only coarse-bound tie is resolved using the explicit
`h5` formula, where `2048 e^2` has degree `18`.  A polynomial identity cannot
have a unique top-degree term.  Contradiction.

Thus T1 is dead, and so is the complete geometrically `q`-coprime `a=7`
stratum.

Again, the finite algebra and degree table are checked by
`t5_split_place_verify.py`.

---

## Correct next program

1. **Replace the scalar `a_q` ledger.**  Generate the sorted multiplicity vectors
   `(b1,b2,b3,b4)` with `a+sum(bi)<=10`.  There are 327 vectors up to permutation
   before terminal pruning.  Store them in a machine-readable ledger rather
   than hand-maintained prose.
2. **Port every proof with explicit hypotheses.**  Mark a result as one of:
   geometrically coprime, uniformly `q^r`-divisible with coprime cofactor, or
   genuinely split-place.  Do not count partial-support vectors as covered.
3. **Use integer valuation pruning first.**  The last two cascade levels give
   linear constraints on the four `q`-root valuations of `e`, `d1`/`sigma`, and
   `g6,g7`; an exact enumerator should discard impossible vectors before any
   symbolic expansion.
4. **Attack the high-`a` coprime rows next.**  The new `a=7` proof shows the
   pattern: UFD parity plus levels 6 and 5 plus infinity.  The same method should
   be tested on `a=8` and `a=9`, where `deg E<=2` and `<=1` respectively.
5. **For `f37`, reintroduce an original equation immediately.**  The resultant
   alone has the exact family `(d2,d1)=(0,0)`, so more computation on the bare
   `f37` identity cannot prove infeasibility.  Restrict one pre-resultant system
   equation to that family and determine the smallest extra generator that cuts
   it.
6. **Stop launching monolithic Gröbner jobs.**  The repository already records a
   ~20 GB msolve kill and a structural NulLA no-go.  Computation should verify
   small factorizations, valuation ledgers, and low-variable terminal systems,
   not materialize the global ideal.
7. **Reconcile `STATE.md`, `HANDOFF.md`, and the default branch.**  The current
   handoff predates the exact campaign, and several survivor maps disagree.
   Generate the human summary from the machine-readable ledger to prevent this
   recurring.

## Bottom line

The program is promising, but the strongest recent survivor-count claims are
currently conditional on a coefficient-field assumption that the plane problem
does not provide.  The correction is finite and structured rather than fatal:
work over the four geometric roots, retain the `t` cascade, and make field-stable
local proofs.  The repaired sigma theorem and the complete coprime `a=7` kill
show that the exact strategy continues to produce results after that correction.
