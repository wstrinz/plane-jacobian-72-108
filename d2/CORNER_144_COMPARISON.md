# Corner (8,28) at degree 144: forcing comparison

## Verdict

**Partial, with a negative answer at cascade-signature level.**  The `(108,144)` case is governed by the same *parametric C-series/forcing operator* as the current `(72,108)` case, after replacing the reduced shared-power pair `(a,b)=(2,3)` by `(3,4)`.  It is not the same forcing divisor or the same split-place ledger.  The other final corner changes the selected root multiplicity from `7` to `3`, leaving a residual quartic in the leading polynomial.  Consequently the target `Phi` signature is `(deg,ord,mult_{y+1},cofactor degree)=(550,205,69,276)`, not `(238,204,30,4)`.  These exact values are checked by `corner144_verify.py` checks **“target Phi signature (550,205,69,276)”** and **“current Phi signature (238,204,30,4)”**.

Thus the experiment supports a template parameterized by **(shared-power exponents, final-root multiplicity)**, not a theorem parameterized by the corner `A0=(8,28)` alone.  This is evidence against reusing the current cascade signature unchanged, but evidence for a broader two-parameter forcing theorem.  The target result below is conditional only at the polygon-reduction boundary stated in “Scope and stopping point”; none of the displayed algebra is numerical or heuristic.

## Sources and convention

GGV5 lists the two length-one chains separately:

- `(8,28) -> (7/4,3)`, table pair `(m,n)=(3,4)`, maximum degree `144`;
- `(8,28) -> (11/4,7)`, table pair `(m,n)=(3,2)`, maximum degree `108`.

This is exactly `paper_src/1708.07936_GGV5.tex` lines 1821-1837.  The degree recipe `deg(P)=m(a+b)`, `deg(Q)=n(a+b)` is stated at lines 248-252.  The two Diophantine checks are re-run exactly by `corner144_verify.py` checks **“(7/4,3),(m,n)=(3,4) satisfies k=1”** and **“(11/4,7),(m,n)=(3,2) satisfies k=1”**.

There is a convention clash in the current case.  GGHV22 calls its original corner pair `(m,n)=(3,2)` at `paper_src/2204.14178.tex` lines 1008-1016, but after reduction its leading forms are powers `P~R^2`, `Q~R^3`; the displayed reduced polygons have the `2`-multiple for `P` and the `3`-multiple for `Q` at lines 1137-1141 and 1181-1186.  This document therefore uses

```text
(a,b) := (power of C in P, power of C in Q).
```

It calls the current reduced pair `(a,b)=(2,3)` and the target pair `(a,b)=(3,4)`.  The degree readings `(72,108)` and `(108,144)` are checked by `corner144_verify.py` checks **“reduced current exponent pair gives (72,108)”** and **“degree recipe gives target (108,144)”**.  Carrying the paper-table label instead changes names only; it does not change any formula below.

## 1. What the final corner changes

For a type-II.b edge, GGV5 writes the new corner as

```text
A1 = (k_i/(m l_i),0) + gamma*(-sigma_i/rho_i,1),
gamma = multiplicity(selected root)/m.
```

This is `paper_src/1708.07936_GGV5.tex` lines 971-982, with the root-multiplicity definition also at lines 589-598.  On the `(8,28)` edge the two listed final corners are exactly

```text
(11/4,7) = (1,0) + 7*(1/4,1),
( 7/4,3) = (1,0) + 3*(1/4,1).
```

Both identities are checked by the two `gamma` checks in `corner144_verify.py` section A.

GGHV22 explicitly identifies the multiplicity-seven edge as `y(x^4 y-alpha)^7` and shifts its root at `paper_src/2204.14178.tex` lines 1132-1136.  For multiplicity three the same length-one-chain transformation leaves a quartic residual.  Before the final Laurent map its shifted base edge has the form

```text
(Y + alpha X^-4) (X^4 Y)^3 H_4(X^4 Y + alpha).
```

The paper's final map is `X -> x^-1`, `Y -> x^4 y` (`paper_src/2204.14178.tex` lines 1228-1234).  Exact substitution gives

```text
x^4 C4(y),   C4(y)=y^3(y+alpha)H_4(y+alpha).
```

The multiplicity-seven specialization gives `x^4 y^7(y+alpha)`.  These are `corner144_verify.py` checks **“multiplicity-three edge gives C4=...”** and **“multiplicity-seven edge gives current C4=...”**.  After scaling `alpha=1`, write

```text
C4 = y^3 (y+1) h4(y),       deg h4=4,
h4(0) h4(-1) != 0.
```

The nonvanishing follows from selecting a nonzero root and removing exactly its multiplicity; the edge polynomial convention `p(0) != 0` is at `paper_src/1708.07936_GGV5.tex` lines 505-520.  The final corner is type I by lines 971-984, and type-I residual roots are simple by lines 576-587.  This is the source of the squarefree residual used in solving the target ODE.

The commutator exponent does **not** change.  The Laurent map has Jacobian `-x^2`, checked by `corner144_verify.py` check **“Jacobian of Laurent map is -x^2”**; this is also the chain-rule computation stated at `paper_src/2204.14178.tex` lines 1228-1230.  Multiplying one member by `-1` therefore normalizes the target reduced pair to `[P,Q]=x^2`.  No `x^2` assumption was imported from the current case.

## 2. C-series normalization and the common formula

The paper constructs the current C-series coefficient by coefficient and then removes commuting powers of `C`; see `paper_src/2204.14178.tex` lines 1416-1468 and 1502-1525.  Applying the same formal root recursion with a cube instead of a square gives the target comparison model

```text
P = C^3,
Q = C^4 + alpha_3 C^3 + alpha_2 C^2 + alpha_1 C
        + alpha_0 + alpha_-1 C^-1 + alpha_-2 C^-2 + F,
ell_10(C)=x^4 C4,             v_10(F)=-9.
```

This is a raw central-series form; some `alpha_j` can be removed by triangular changes, but they commute with `P` and do not enter the forcing equation.  The stopping index is exact because `v(C^-2)=-8>-9` whereas `v(C^-3)=-12<-9`; see `corner144_verify.py` check **“target commuting correction range ends at C^-2”**.

More generally, let

```text
ell(C)=x^t c(y),  P=C^a,  Q=C^b+(commuting powers)+F,
[P,Q]=x^kappa,    f=c^b F_s.
```

Direct differentiation gives

```text
s = kappa+1-a t,
a { t c f' - [t(b-a)+kappa+1] c' f } = c^(b-a+1).      (F)
```

Both specializations are independently differentiated in `corner144_verify.py` checks **“direct bracket gives general formula at (2,3)”** and **“... at (3,4)”**.  For this corner `t=4`, `kappa=2`, and both cases have `b=a+1`, so (F) becomes the single family

```text
4a c f' - 7a c' f = c^2.                                (F_a)
```

It gives

```text
current a=2:  8 c f' - 14 c' f = c^2,  s=-5;
target  a=3: 12 c f' - 21 c' f = c^2,  s=-9.
```

The obstruction degrees and the two exact equations are checked in `corner144_verify.py` sections C-E.  In particular, substituting `c=y^7(y+1)` and dividing by `y^6` recovers the current displayed ODE

```text
8y(y+1)f' - 14(8y+7)f = y^8(y+1)^2.
```

The current solution and its separable quartic are rechecked by `corner144_verify.py` checks **“current f solves its ODE”** and **“current quartic separable and avoids 0,-1”**.

## 3. Solving the target ODE

Put `g=(y+1)h4`, so `c=y^3 g`, `deg g=5`, `g(0)!=0`, and `g(-1)=0`.  Since the final residual is squarefree (GGV5 lines 576-587 and 971-984), local orders in `(F_3)` force a polynomial solution of degree at most `14` to have

```text
f = A y^4 g^2.
```

Substitution reduces the entire ODE to

```text
3 A (y g' - 5g) = 1.
```

The local orders and the degree-14 resonance are checked by `corner144_verify.py` checks **“local target orders are 4 at y and 2 at every simple g-root”** and **“infinity leading coefficient is resonant exactly at degree 14”**.  The substitution reduction is checked by **“squarefree-order ansatz reduces to 3A(yg'-5g)=1”**.  Its coefficients force `g_1=...=g_4=0`; the degree-five term is resonant; and `g(-1)=0` forces `g_0=g_5`.  Those three coefficient statements are separate exact checks in script section E.

After monic normalization,

```text
g = y^5+1 = (y+1)(y^4-y^3+y^2-y+1),
h4 = y^4-y^3+y^2-y+1,
C4 = y^3(y^5+1),
f = -(1/15)y^4(y^5+1)^2.
```

The quartic is separable and avoids `0,-1`; the displayed `f` solves the ODE; and an independent 15-variable linear solve finds it as the unique polynomial solution of degree at most `14`.  These are the four checks beginning **“normalized g...”** through **“linear solve returns displayed target f”** in `corner144_verify.py` section E.

This is a separable-quartic phenomenon, but not the current one.  In the current case the quartic `q` is a new unit cofactor of `f`; in the target case `h4` is already part of `C4`, and `f` introduces no new quartic place.

## 4. D-transformation and Phi

For `P=C^a` with `ell(C)=x^4 c`, the denominator-clearing transform is

```text
D_k := C_k c^[4a-1-a k] = C_k c^[a(4-k)-1].              (D)
```

It specializes to

```text
current: D_k=C_k C4^(7-2k),
target:  D_k=C_k C4^(11-3k).
```

All integer specializations in a range containing the used slices are checked by `corner144_verify.py` check **“D exponents are 7-2k and 11-3k”**.  Formula (D) is the cube-root version of the coefficient-clearing recurrence used in the paper's C construction at `paper_src/2204.14178.tex` lines 1427-1468.

At the target forcing slice `j=9`, the `C^4` coefficient is cleared by `C4^71`.  Since `f=C4^4 F_-9`, the target forcing divisor is

```text
Phi_144 := F_-9 C4^71 = f C4^67
         = -(1/15) y^205 (y^5+1)^69
         = -(1/15) y^205 (y+1)^69 h4^69.
```

The exponents `71` and `67` are checked by `corner144_verify.py` checks **“target Q-slice clearing exponent is 71”** and **“target obstruction degree and Phi exponent are -9,67”**.  The final factor exponents and all divisor invariants are checked without expansion by the four target-Phi checks in script section E.

Consequently:

| Signature component | current `(a,b)=(2,3)` | target `(a,b)=(3,4)` | one formula? |
|---|---:|---:|---|
| reduced bracket | `x^2` | `x^2` | yes, same Laurent Jacobian |
| leading `R` | `x^4 y^7(y+1)` | `x^4 y^3(y+1)h4` | same `x^4`, no for divisor |
| shared powers | `P~C^2,Q~C^3` | `P~C^3,Q~C^4` | yes: `(a,a+1)` |
| `v(F)` | `-5` | `-9` | yes: `3-4a` |
| D exponent | `7-2k` | `11-3k` | yes: `a(4-k)-1` |
| ODE | `8cf'-14c'f=c^2` | `12cf'-21c'f=c^2` | yes: `(4a,7a)` |
| `f` degree/order | `14 / 8` | `14 / 4` | degree yes; order depends on root multiplicity |
| Phi power after `f` | `C4^28` | `C4^67` | yes: `a[4(a+b)-3]-2b` |
| `deg Phi` | `238` | `550` | formula also needs divisor profile |
| `ord_y Phi` | `204` | `205` | not corner-only |
| `mult_(y+1) Phi` | `30` | `69` | not corner-only |
| unit cofactor | quartic `q` | `h4^69` | **no** |
| cofactor degree | `4` | `276` | **no** |

Every numerical or polynomial entry in this table is checked by `corner144_verify.py` sections C-E.  The displayed Phi-exponent formula specializes to `28` and `67` in the same checks.

## 5. Window/envelope candidates

This is the only deliberately conditional row of the comparison.  Extending the GGHV22 polygon move to selected multiplicity three gives lower base point `(8,3)` before the last Laurent map.  GGV5's final type-I endpoints are the two points in `paper_src/1708.07936_GGV5.tex` lines 1445-1455.  Parallelism swaps their assignment relative to the current case: `P` takes `(2,1)` and `Q` takes `(-1,0)`.  The resulting candidate lower edges are parallel; this is checked by `corner144_verify.py` checks **“target lower candidate edges are parallel”** and **“pre-map lower slope is 4/11”**.

With weight `w=4-k`, `deg C4=8`, `ord C4=3`, and target D exponent `3w-1`, the candidate caps are

```text
upper, no y-axis extra corner: deg D_(4-w) <= 22w,
upper, analogous extra corner: deg D_(4-w) <= 23w,
lower: ord D_(4-w) >= 8w + ceil(w/5).
```

The arithmetic is checked for `w=1,...,12` by the twelve **“envelope arithmetic”** checks in `corner144_verify.py`.  Unlike the current `14w/15w` upper and `12w` lower caps, the target lower candidate is quasipolynomial rather than affine.  It is therefore another warning against copying the existing cascade weights.  It is not promoted to a polygon theorem here.

## 6. Scope and stopping point

GGHV22 proves the complete reduction to `[P,Q]=x^2` and two explicit polygons only for the `(3,2)` table case; see `paper_src/2204.14178.tex` lines 1000-1007 and the proof at lines 1008-1311.  GGV5 lists the `(3,4)` chain but does not carry out the analogous explicit Newton-polygon reduction (`paper_src/1708.07936_GGV5.tex` lines 1821-1837).

Accordingly this experiment stops before claiming a complete `(108,144)` polygon theorem.  It carries the standard length-one-chain root shift and Laurent map far enough to determine the leading `C4`, derive and solve the forcing ODE, compute `Phi`, and propose envelope weights.  Completing all opposite-vertex and subcase enumeration would reproduce the multi-page argument at GGHV22 lines 1008-1311 and is outside this comparison timebox.

The exact conditional boundary is:

> If the GGV5 length-one chain `(8,28)->(7/4,3)` is realized by the standard type-II.b root shift and the same final Laurent chart as the published `(11/4,7)` reduction, then all target formulas and invariants above follow exactly.

This is the reading consistent with `regenerate_system.py` and `verify_derivation.py`: reduced powers are ordered by `deg P < deg Q`, and the commutator is normalized after the Laurent chart.  The alternative paper-table reading merely swaps the labels of the original pair; if it also swaps `P,Q`, the sign of the bracket changes and is removed by `Q -> -Q`.

## 7. Next falsification test

The shortest decisive test is to derive only the two target edges adjacent to the transformed `(4,3)` base point, without enumerating the full polygon.  It should confirm or refute:

```text
exceptional endpoints P:(2,1), Q:(1,0) after the Laurent map,
lower direction (4,-5),
upper alternatives yielding 22w or 23w.
```

These candidates are exactly the geometry used by script section F.  If any endpoint differs, the ODE and Phi computation survives—their inputs are the leading `C4` and `[P,Q]=x^2`—but the window candidates must be discarded.  If the residual quartic is not the final type-I polynomial to which GGV5 lines 576-587 apply, rerun the ODE without the squarefree assumption; that is the cleanest falsification of the forced cyclotomic `h4`.

## Proposed STATE.md entry

> **Corner-144 comparison (conditional forcing result).**  For the GGV5 length-one chain `(8,28)->(7/4,3)`, table pair `(3,4)`, the standard root-shift/Laurent-chart continuation has `[P,Q]=x^2`, `P=C^3`, `Q=C^4+(commuting C-powers)+F`, `v(F)=-9`, and `ell(C)=x^4 C4` with `C4=y^3(y+1)h4`.  The forcing family is `12 C4 f' -21 C4' f=C4^2`, `f=C4^4 F_-9`.  Type-I squarefreeness forces `h4=y^4-y^3+y^2-y+1`, `f=-y^4(y^5+1)^2/15`, and `Phi_144=f C4^67=-y^205(y^5+1)^69/15`, hence `(deg,ord,mult_{y+1},cofactor degree)=(550,205,69,276)`.  The operator and D-exponents fit the parametric formulas (`D_k=C_k C4^(11-3k)`), but the forcing-divisor/cascade signature does not match the `(72,108)` case.  Candidate target windows are upper `22w/23w`, lower `8w+ceil(w/5)` pending a full polygon reduction.  Exact checks: `corner144_verify.py` (51/51 pass).  Verdict: partial support for a `(power pair, final multiplicity)` template; against an `A0`-only cascade template.

