# F10 corner: exact enumeration of the real mu-rungs

## Result

For the F10 forcing ODE, the real branch support is

```
mu = 1   EMPTY
mu = 2   REALIZED (two irreducible components, four real points total)
mu = 3   EMPTY
mu = 4   REALIZED (the known fully ramified rational point)
```

Thus the even-`dg` parity claim is **proved at `dg=4`**: the `mu=1` rung is
empty.  The open intermediate rungs split: `mu=2` exists, while `mu=3` does
not over the reals.

Derivation: `mu_rungs_f10.py`.  Independent checker:
`mu_rungs_f10_verify.py` (**39/39**, including the F12 `mu=2` quartic control).
Both scripts use exact SymPy arithmetic only.

## 1. Reconstructed F10 system

The F10 row gives

```
(a,b,t,kappa,q,e,dg,rho,N,res,gap,r) =
(4,7,7,5,3,4,4,10,270,27,1,3).
```

Writing `c = y^3 g`, the full forcing ODE is

```
4 (7 c f' - 27 c' f) = c^4.                         (1)
```

The residual normalization has `ord_y(c)=q=3`, so `g(0) != 0`; after
`g=(y+1)^mu h`, this is `h(0) != 0`.  The rung condition itself is
`h(-1) != 0`.  At a root of `g` having multiplicity `m`, (1) forces the
local order of `f` to be `3m+1`.  Enumerating the root partitions of `h`
therefore loses no branch, including repeated roots away from `-1`.

For `Phi = f c^N`, the branch-independent entries are

```
deg Phi = 27 + 270*7 = 1917,
ord_y Phi = 10 + 270*3 = 820.
```

Here `e+N=274` and `e+N-1=273`.

## 2. Per-rung derivation

### mu=1: EMPTY (parity claim proved at dg=4)

All three partitions of the cubic `h` were treated.

**Squarefree cubic.**  Put

```
h = y^3 + p y^2 + q y + r,
f = y^10 (y+1)^4 h^4 (u0+u1 y).
```

After solving for `u0,u1`, the Groebner basis of the three remaining
conditions contains

```
r (2250 r^4 - 4776 r^3 + 1180 r^2 + 75 r + 2250).
```

The factor `r=0` is inadmissible (`h(0)=0`).  The quartic

```
M_111(r) = 2250 r^4 - 4776 r^3 + 1180 r^2 + 75 r + 2250
```

has exactly **zero real roots** by an exact Sturm count.

**One double and one simple root.**  For
`h=(y-z)^2(y-w)`, local orders force

```
f = y^10 (y+1)^4 (y-z)^7 (y-w)^4 u,   deg u=2.
```

The two residual conditions have resultant (eliminating `z`)

```
6144 w^6 M_21(w),
```

where

```
M_21(w) = 1250w^12 - 1500w^11 + 4400w^10 + 3850w^9
          + 6303w^8 - 8272w^7 - 11792w^6 - 8272w^5
          + 6303w^4 + 3850w^3 + 4400w^2 - 1500w + 1250.
```

`w=0` is inadmissible, and `M_21` has exactly **zero real roots**.

**Triple root.**  For `h=(y-z)^3`, the sole obstruction (up to a nonzero
rational scalar) is

```
M_3(z) = 5z^4 - 2z^3 + 8z^2 + 20z + 10,
```

which also has exactly **zero real roots**.

These are all partitions `1+1+1`, `2+1`, and `3`; hence `mu=1` is proved
empty over the real F10 branch system.

### mu=2: REALIZED

There are two distinct irreducible components.

#### A. Squarefree quadratic component

Let

```
h = y^2 + p y + r,
g = (y+1)^2 h,
f = y^10 (y+1)^7 h^4 u,
```

with

```
u = -1/(44r)
    + (5p+8r)/(88r^2) y
    + (5p^2+13pr+24r^2-6r)/(88r^3) y^2.
```

The two consistency equations are

```
10p^3 + p^2r + 3pr^2 - 17pr - 24r^3 + 12r^2 = 0,
-5p^3 - 13p^2r + 25p^2 - 24pr^2 + 36pr + 64r^2 - 30r = 0.
```

Their resultant in `r` is

```
-4000 p^3 M_sf(p),
M_sf(p) = 192p^6 - 720p^5 + 586p^4 - 668p^3
          + 1728p^2 - 1510p + 45.
```

The `p=0` factor makes the two equations compatible only at the inadmissible
`r=0`.  The polynomial `M_sf` is irreducible over `Q` and has exactly two
real roots, isolated by `(0,1)` and `(2,3)`.  On this component,

```
r = (-768p^5 + 2400p^4 - 844p^3 + 3492p^2 - 6077p + 1905)/3773.
```

Exact gcd checks modulo `M_sf` prove `r != 0`, `1-p+r != 0`, and
`p^2-4r != 0`, so both real points are admissible and `h` is genuinely
squarefree.  Substitution into the **full un-divided ODE (1)** is zero
coefficient-by-coefficient in `Q[p]/(M_sf)`.

#### B. Double-root quadratic component

Let

```
h = (y-z)^2,
g = (y+1)^2 (y-z)^2,
f = y^10 (y+1)^7 (y-z)^7 u,
```

where

```
u = 1/(44z)
    + (1-z)/(11z^2) y
    + (-12z^2+17z-12)/(44z^3) y^2
    + (-12z^3+57z^2-57z+12)/(220z^4) y^3.
```

The sole condition is the irreducible quartic

```
M_d(z) = 12z^4 + 15z^3 - 5z^2 + 15z + 12 = 0.
```

It has exactly two real roots, one in `(-2,-1)` and one in `(-1,0)`;
neither is `0` or `-1`.  The full ODE is identically zero in
`Q[z]/(M_d)`.

Both `mu=2` components have

```
Phi signature = (1917, 820, 547, 550).
```

The mu-graded law predicts

```
mult = 2*274 - 1 = 547,
cof  = 1 + 3*274 - 273 = 550,
```

so the verdict is **MATCH**.

### mu=3: EMPTY

Write `h=y+s`.  Solving the four linear coefficients of

```
f = y^10 (y+1)^10 (y+s)^4 u,   deg u=3,
```

leaves, up to a nonzero rational scalar,

```
M_mu3(s) = 10s^4 - 20s^3 + 8s^2 + 2s + 5.
```

This quartic is irreducible over `Q` and has exactly **zero real roots**.
(It has four nonreal algebraic roots, hence complex-coefficient formal
solutions, but no real F10 branch.)  Therefore `mu=3` is EMPTY.

### mu=4: REALIZED

Here `h=1`, `g=(y+1)^4`, and

```
f = y^10 (y+1)^13
    (2401y^4 + 5831y^3 + 4165y^2 + 595y - 85)/3740.
```

Direct substitution into (1) gives zero exactly.  If `F27` denotes the
resonant leading coefficient, its minimal polynomial is

```
3740 F27 - 2401.
```

The signature is

```
Phi signature = (1917, 820, 1093, 4),
```

matching the previously derived ramified point.  The mu-graded prediction is

```
mult = 4*274 - 3 = 1093,
cof  = 1 + 3*274 - 3*273 = 4,
```

so the verdict is **MATCH**.

## 3. Summary table

| mu | classification | exact parameter data | Phi signature | mu-law |
|---:|---|---|---|---|
| 1 | **EMPTY** | three partition obstructions above; each has 0 real roots | -- | -- |
| 2 | **REALIZED** | squarefree `M_sf(p)` (2 real) and double-root `M_d(z)` (2 real) | `(1917,820,547,550)` | **MATCH** |
| 3 | **EMPTY** | `10s^4-20s^3+8s^2+2s+5` has 0 real roots | -- | -- |
| 4 | **REALIZED** | `3740F27-2401`; `g=(y+1)^4` | `(1917,820,1093,4)` | **MATCH** |

## 4. Independent verification and control

`mu_rungs_f10_verify.py` does not import the derivation script.  It rebuilds
the root-partition systems, repeats every resultant/Groebner obstruction and
exact Sturm count, and substitutes each realized solution into the full ODE.
For algebraic parameters, equality is checked coefficientwise modulo the
stated irreducible minimal polynomial.

As a control, it independently reconstructs F12's standard-family `mu=2`
condition

```
195 beta^4 + 120 beta^3 - 40 beta^2 + 32 beta - 80,
```

checks its two real roots and exclusions `beta != 0,-1`, and verifies the
full F12 ODE exactly modulo that quartic.

## `[judgment]` items

1. **[real branch field]** Classification is over real coefficients, as in
   the corner/root-shift problem.  Nonreal parameter roots are recorded as
   formal complex solutions but do not realize a real rung.
2. **[residual normalization]** `h(0)!=0` is not an added branch-selection
   assumption: it is required by the defining normalization `c=y^q g` with
   `q=ord_y(c)=3`.  Eliminant powers coming only from `h(0)=0` are boundary
   artifacts.
3. **[partition completeness]** The local ODE order `ord_root(f)=3m+1` was
   used for each multiplicity `m`.  All partitions of `deg h` are treated:
   three at `mu=1`, two at `mu=2`, and one each at `mu=3,4`.
4. **[ODE versus tower selection]** This is a complete enumeration of the
   real forcing-ODE rungs under the F10 corner ansatz.  It does not build the
   unreduced polygon's tower C-series; as in PHI_F7.md, that upstream model
   selection remains a separate issue.
5. **[cofactor dependence]** `cof=deg-ord-mult`, so each realized rung adds
   one independent new signature datum (`mult`); the reported cofactor check
   is nevertheless performed explicitly against the requested mu-graded law.
