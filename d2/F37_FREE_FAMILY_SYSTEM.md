# The f37 free family is excluded by the pre-resultant system

**Date:** 2026-07-22
**Scope:** subcase (2), on the resultant locus `d2=d1=0`.

## Result

The polynomial `f37` vanishes identically when `d2=d1=0`, with `d0,e`
arbitrary in their windows.  This is a genuine free family for the bare
resultant, but **none of it lifts to the original pre-resultant system**.

The proof is field-stable: it is carried out after base change to an algebraic
closure and uses only DVR valuations, the degree windows, and characteristic
different from two and three.

## 1. Restore the smallest original system

Put

```text
e = d_{-1},  r = d_{-2},  s = d_{-3}.
```

Restrict the regenerated equations to `d2=d1=0`, solve `G1` for
`d_{-4}`, and clear the same harmless power of `e` used in the resultant
derivation.  Directly from `t4_state.pkl` one obtains

```text
H2 = -3 (d0 e^3 - e s^2 + 2 r^2 s),
H3 = -(6 d0 e^2 r + e^4 + 6 r s^2),
H5 = e (2 Phi - 3 e^2 s - 3 e r^2).
```

The branch `e=0` is already globally impossible, so cancel `e` in `H5`.
Eliminating `d0` from `H2,H3` gives the compact necessary system

```text
12 r s (r^2 - e s) = e^5,                 (P)
3 e (r^2 + e s) = 2 Phi.                  (S)
```

The subcase-(2) windows are

```text
deg e <= 10,  deg r <= 12,  deg s <= 14,
Phi = c t^30 q,  deg Phi = 34,
t=y+1,  q=p1 p2 p3 p4 squarefree.
```

All displayed identities and the elimination are checked by
`f37_free_family_verify.py` from the regenerated state rather than embedded
resultant data.

## 2. First consequence: e divides Phi

Equation (S) immediately gives

```text
e | Phi.
```

Thus `e` has no roots away from `t,p1,p2,p3,p4`.  At a simple root of
`q`, its multiplicity in `e` is either zero or one.

This already repairs the coefficient-field issue on this branch: partial
support is allowed, but it is now a subset of the four simple roots rather
than an arbitrary multiplicity vector.

## 3. Complete local valuation table

At a root of `Phi` of order `phi`, write

```text
m=v(e),  x=v(r),  y=v(s),  z=v(r^2-es).
```

Equation (P) gives

```text
x+y+z = 5m,
```

and (S) gives

```text
v(r^2+es) = phi-m.
```

For `X=r^2`, `Y=es`, characteristic zero gives

```text
min(v(X+Y),v(X-Y)) = min(v(X),v(Y)).
```

The finite local enumeration is therefore exhaustive.

At an unselected simple `q)-root (`m=0`):

```text
(x,y,z) = (0,0,0).
```

At a selected simple `q)-root (`m=1`):

```text
(x,y,z) = (0,5,0).
```

At `t`, where `phi=30` and `0<=m<=10`, the only possibilities are

```text
(m,x,y,z) =
    (0,  0,  0,  0),
    (5,  6,  7, 12),
    (9, 12, 12, 21),
    (10,10, 10, 30).
```

## 4. Infinity forces e=C t^10

Let `E=deg e`.  Since `e|Phi), equation (S) can be divided by `e):

```text
r^2 + e s = 2 Phi/(3e).
```

The left side has degree at most

```text
max(2 deg r, deg e + deg s) <= max(24,E+14) <= 24,
```

whereas the right side has degree `34-E`.  Hence `E=10`.

If `k` simple `q)-roots are selected, then `E=m+k), and the local
`q)-table forces `deg s >= y+5k`.  Combining this with the four
`t)-options and `deg s<=14` leaves only

```text
m=10, k=0,
e=C t^10, v_t(r)=v_t(s)=10, v_t(r^2-es)=30.
```

## 5. Final contradiction

Write

```text
r=t^10 R,  deg R<=2,
s=t^10 S,  deg S<=4.
```

Then

```text
r^2-es = t^20 (R^2-C S).
```

The required order `v_t(r^2-es)=30` says

```text
t^10 | R^2-C S.
```

But `deg(R^2-C S)<=4), so `R^2-C S=0).  Equation (P) then has zero left
side and nonzero right side `e^5), a contradiction.

Therefore no point of the exact `f37` free family lifts to the original
pre-resultant equations.

## Consequence for the program

The bare resultant still cannot prove its own free family empty, but restoring
only the three small equations above removes it completely.  Future `f37`
work may therefore focus on the complementary locus `(d2,d1)!=(0,0)`; it
does not need to carry a 20-coefficient exceptional family.
