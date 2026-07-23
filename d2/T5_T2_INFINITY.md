# T5 T2 infinity convolution: exact partial closure of the eight open cells

**Date:** 2026-07-22  
**Branch:** `f31`, subcase (2), `d1=0`,
`sigma=4d0-d2^2 != 0`.

This note continues `T5_T2_COLUMN.md`.  It does not close an entire T2 cell,
but it strictly narrows all eight.  The exact ledger is

```text
cells killed / narrowed / unchanged:       0 / 8 / 0
flag cases killed:                         a9 b1000 G5
additional degree-state removals:          11
  pattern A (all flags in their cells):     2
  d2=0 fixed-F states:                      9
```

All computations are checked by `t5_t2_infinity_verify.py`.  Its `H_j` are
read by the established regex parser in `t5_t2_column_verify.py` from
`f31_graded.txt`; no `h_j` coefficient is copied into the new checker.

## 1. Exact residual parametrization

Work over a splitting field and normalize the four quartic factors to
`p_i=y-r_i`.  Thus `q` is a nonzero constant times `p_1...p_4`.  For a cell
with profiles `b=(b_i)` and `s=(s_i)`, put

```text
R = product p_i^b_i,                 P = product p_i^s_i,
e = epsilon t^a R F,                 sigma = zeta P Z,
F = y^f + u_(f-1)y^(f-1)+...+u_0,
Z = y^z + v_(z-1)y^(z-1)+...+v_0.
```

Here `epsilon,zeta` are nonzero, `gcd(F,tq)=1`, `gcd(Z,q)=1`, and the degrees
`(f,z,g;D,Sigma)` are exactly those in the residual lists below.  The
terminal equation and the divisibility already proved in the column become

```text
F^3 G = 3072 c^6 Z^2,                G=F^2 W,
F^5 W = 3072 c^6 Z^2.                                      (1)
```

This is also an explicit parametrization of the terminal free quotient:
`W` must be a polynomial satisfying the second equation.  With
`Q=product p_i^m_i` and `Qbar=Q/q`, the next quotients are forced, not new
independent parameters:

```text
g5 = [t^v q Qbar G - 2048 c^5 q^5 t^(2a) R^2 F^2]/(R^3 F^3)
     - 19890 d2 Qbar G,                                      (2)

g4 = [t^v g5 - c^4 q^4 h4]/(R^3 F^3).                       (3)
```

Polynomiality of (2) and (3), together with `g4=0` in case `G4`, is part of
the exact residual.  Check I5 verifies the two solved quotient formulas.

The remaining `d2` parameters are

```text
N or G4: d2=d_4 y^4+...+d_0,  d2(-1) != 0;
D:       d2=0;
G5(a9):  d2=t^3(lambda_1 y+lambda_0), lambda_0-lambda_1 != 0.
```

The last line includes `lambda_1=0`, i.e. degree three.

## 2. One exact convolution for every pattern-B state

Substitute `d0=(d2^2+sigma)/4` in the source-parsed master identity and set

```text
M(P,E,D,S) = sum_(j=0)^7 P^j E^(21-3j)
             h_j((D^2+S)/4,0,D,E).                           (4)
```

Every monomial in (4) has weighted degree 250 for weights

```text
wt(P,E,D,S)=(34,10,4,8).                                    (5)
```

This is the promised common computation.  If `x=1/y` and

```text
Phi = y^34(P_0+P_1 x+...),       e     = y^10(E_0+E_1 x+...),
d2  = y^4 (D_0+D_1 x+...),       sigma = y^8 (S_0+S_1 x+...),
```

then the master coefficients are exactly

```text
C_250 = M(P_0,E_0,D_0,S_0),                                (6)
C_249 = P_1 M_P + E_1 M_E + D_1 M_D + S_1 M_S,             (7)
```

where the derivatives in (7) are evaluated at the four leading
coefficients.  Check I1 verifies (5)--(7) term by term from all parsed
`h_j`.  In particular, (4), rather than a truncated `h5+h6` expression, is
used: lower `h_j` terms enter automatically as soon as their reversed-series
order is reached.

Equations (6) and (7) are imposed on every surviving pattern-B state.  The
precise obstruction to a uniform further descent is already visible in
(6): `D_0=d_4` is free in cases `N,G4`, and `S_0` is also free when
`Sigma=8`.  The polynomial `M(P_0,E_0,d_4,S_0)` is nonzero but has a
positive-dimensional zero set.  Consequently its cancellation does not
force a unique coefficient substitution.  Splitting all of its algebraic
components before applying (7) is the missing finite elimination; treating
the tie as noncancellation would be invalid.

## 3. Terminating descents

### 3.1 Pattern A: `(D,Sigma)=(8,3)`

Write `e=gamma E` and `sigma=beta S` with `E,S` monic of degrees 8 and 3.
Source parsing shows that every master term other than the `e^2` part of
`h5` and the `sigma^2` part of `h6` has degree at most 228.  Above that
degree the identity is therefore the nonzero common factor times

```text
2 gamma^5 E^5 - 3 Phi beta^2 S^2.                            (8)
```

Let `phi_34=lc(Phi)`.  Degree 234 first forces

```text
beta^2 = 2 gamma^5/(3 phi_34).                               (9)
```

For `a8 b0000`, put

```text
E=t^8,             S=y^3+s_2 y^2+s_1 y+s_0.
```

After (9), degrees 233, 232, and 231 successively force

```text
s_2=41/8,          s_1=1353/128,       s_0=11275/1024.
```

The normalized degree-230 coefficient is then
`191675/16384 != 0`.  Thus `(0,3,6;8,3)` is killed in all three flag cases.

For `a7 b1000`, choose the supported root `p_1=y-r` and write

```text
E=t^7(y-r),        S=(y-r)^2(y-u).
```

Degree 233 forces `u=r/2-21/8`.  The normalized degree-232 coefficient is

```text
-(16r^2+168r-273)/64.                                      (10)
```

The numerator in (10) is coprime to the source quartic `q(r)`, so it cannot
vanish at a quartic root.  Hence `(0,1,2;8,3)` is killed in all three flag
cases.  Check I2 performs every substitution and both nonzero tests.

### 3.2 Pattern B with `d2=0` and `Sigma<=6`

Here `S_0=S_1=0` in (6)--(7).  Only the source `h5` `e^2` monomial and the
source `h0` `e^4` monomial can occur at degrees 250 and 249.  The first
coefficient fixes the seventeenth power of `E_0`.  Substituting that relation
in the second coefficient forces the normalized next coefficient

```text
E_1/E_0 = (5/17)(Phi_33/Phi_34) = 35/4.                     (11)
```

For each fixed-`F` support this contradicts the quartic roots:

* `a9 b1000`: `E_1/E_0=9-r`; (11) requires `r=1/4`, but `q(1/4)!=0`.
* `a8 b1100`: `E_1/E_0=8-(r_1+r_2)`; (11) requires a pair sum `-3/4`.
  The checker takes the gcd of `q(r)` and `q(-3/4-r)` and obtains 1.
* `a7 b1110`: writing the omitted fourth root as `r`, the coefficient is
  `27/4+r`; (11) requires `r=2`, but `q(2)!=0`.

Exactly these nine `D` degree states have `Sigma<=6` and are killed at
degree 249 [I3]:

```text
a9 b1000:  z=0,1,2,3,4;       states Sigma=2,...,6;
a8 b1100:  z=0,1,2;           states Sigma=4,...,6;
a7 b1110:  z=0;               state  Sigma=6.
```

At `Sigma=7`, a sigma-linear term first appears in degree 249, so (11) is no
longer forced.  At `Sigma=8`, sigma terms already occur in degree 250.  These
are exact stopping points, not omitted cases.

### 3.3 The rigid `a9 b1000 G5` flag

The column had already forced `(f,z,g;D,Sigma)=(0,6,12;10,8)`.  Put

```text
p=p_1,       e=gamma t^9 p,       sigma=zeta p^2 Z,
d2=t^3 L,    deg Z=6,              deg L<=1.
```

The terminal equation gives `G=unit*Z^2`; the q-profile gives
`Q=unit*q^6p`.  Consequently the exact flag equation
`t^3g6=c^5q^5h5` says

```text
h5 = mu t^3 q p Z^2,              mu != 0.
```

Using the source-parsed collapse of `h5`, dividing by `t^3p`, and absorbing
nonzero normalization units gives [I4]

```text
2048 gamma^2 t^15 p
 + Z^2(-9216 zeta^2 L p^3-mu q) = 0.                         (12)
```

Thus `Z^2 | t^15p`.  But the t-place witness gives `t` not dividing `Z`,
and `gcd(Z,q)=1` gives `p` not dividing `Z`.  Since `deg Z=6`, (12) is
impossible in the splitting-field UFD.  The complete `G5` flag case is
therefore **killed**.  This is stronger than another unresolved leading
coefficient tie; it uses the exact level-5 equation after the infinity
narrowing.

## 4. Per-cell ledger

In the table, `B250/249` means that the residual has the explicit equations
(6) and (7), with its supported `e,sigma` coefficients substituted.  This is
a genuine narrowing, but not a kill: the degree-250 equation has several
algebraic leading-coefficient branches and no unique forced substitution.

| cell | flag | exact residual after this attack | verdict |
|:--|:--|:--|:--|
| `a9 b1000` | `N` | `z=0,...,6`, plus `B250/249` and (1)--(3) | narrowed |
|  | `G4` | `z=0,...,6`, `B250/249`, (1)--(3), `g4=0` | narrowed |
|  | `G5` | none, by (12) | **killed** |
|  | `D` | only `z=5,6` (`Sigma=7,8`), plus `B250/249` | narrowed |
| `a8 b0000` | `N,G4,D` | `(2,z,2z-6;10,z)`, `z=5,6,7,8`, plus `B250/249`; the `(0,3,6;8,3)` state is killed | narrowed |
| `a8 b1000` | `N,G4,D` | the four original `(1,z,2z-3;10,2+z)`, `z=3,4,5,6`, plus `B250/249` | narrowed |
| `a8 b1100` | `N,G4` | `z=0,...,4`, plus `B250/249` | narrowed |
|  | `D` | only `z=3,4` (`Sigma=7,8`), plus `B250/249` | narrowed |
| `a7 b1000` | `N,G4,D` | `(2,5,4;10,7)` and `(2,6,6;10,8)`, plus `B250/249`; the `(0,1,2;8,3)` state is killed | narrowed |
| `a7 b1100` | `N,G4,D` | `(1,3,3;10,7)`, `(1,4,5;10,8)`, plus `B250/249` | narrowed |
| `a7 b1110` | `N,G4` | `z=0,1,2`, plus `B250/249` | narrowed |
|  | `D` | only `z=1,2` (`Sigma=7,8`), plus `B250/249` | narrowed |
| `a7 b3000` | `N,G4` | `z=0,1`, plus `B250/249` | narrowed |

For `G4`, equation (3) is additionally set equal to zero.  For `D`, the
source identities `h5=2048e^2` and `h4=-12096sigma^3` remain imposed.  The
table does not discard these flag equations merely because the common
infinity notation suppresses them.

## 5. Deepest convolution degree and honest stopping point

| cell | deepest master degree checked | reason for termination |
|:--|--:|:--|
| `a9 b1000` | 249 | five `D` states close there; `G5` then closes by (12); other flags branch at (6) |
| `a8 b0000` | 230 | pattern-A nonzero coefficient; remaining pattern B branches at (6) |
| `a8 b1000` | 249 | non-unique degree-250 leading variety |
| `a8 b1100` | 249 | three `D` states close; other branches are non-unique |
| `a7 b1000` | 232 | pattern-A nonzero coefficient; remaining pattern B branches at (6) |
| `a7 b1100` | 249 | non-unique degree-250 leading variety |
| `a7 b1110` | 249 | one `D` state closes; other branches are non-unique |
| `a7 b3000` | 249 | non-unique degree-250 leading variety |

The next exact task for every still-open pattern-B entry is a component-wise
elimination of

```text
M(Phi_34,E_0,d_4,S_0)=0
```

followed by (7) on each component.  When `d_4=0` or `S_0=0`, the newly
activated `d_3` or `S_1` term in (7) must be retained.  This is precisely why
a single forced-substitution chain analogous to the T1 constant cell does
not yet exist on these branches.

## 6. Verification

Run from this directory:

```text
python t5_t2_infinity_verify.py
```

The passing checks are:

| check | content |
|:--|:--|
| I1 | parses `f31_graded.txt`; proves weighted homogeneity; checks the common degree-250/249 convolution and all 29 pattern-B states |
| I2 | performs every pattern-A forced substitution and verifies the degree-230 and degree-232 nonzero obstructions |
| I3 | derives `E_1/E_0=35/4` from parsed coefficients and verifies all nine `D` state kills against the source quartic |
| I4 | derives (12) from the parsed `h5` collapse and checks the `G5` support data used by the UFD kill |
| I5 | checks the forced `g5,g4` quotient parametrization |

The script exits nonzero on any failed source collapse, convolution identity,
state inventory, forced substitution, quartic gcd, or final nonzero test.
