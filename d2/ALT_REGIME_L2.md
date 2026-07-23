# Alternate regime: flipped levels 2 and 3

Date: 2026-07-22  
Verifier: `alt_regime_l2_verify.py`  
Scope: the 33 branches left open by `ALT_REGIME.md`.

## Result

The exact `h_6` and `h_5` order cones kill **6 of the 33** input branches.
All six are T1 branches. The remaining frontier is **27 branches**: 13 T1
and 14 T2.

This round treats every tied leading form as capable of cancelling to
arbitrary depth. That enlarges the local solution set. Thus an empty cone is
a proof of infeasibility, while a nonempty cone is recorded as an open
residue problem and is not promoted to a solution.

## 1. Exact flipped transitions

Write

```text
e=t^a E,  u=cq,  w=3a-30,  T=t^w,
G'=sum_{f=0}^7 T^(7-f) u^f E^(21-3f) h_f.
```

At `t`, `E` and `u` are units. The polynomial descending chain is

```text
h_7                         = T r_6,
E^3 h_6 + u r_6             = T r_5,
E^6 h_5 + u r_5             = T r_4,
...
E^21 h_0 + u r_0            = 0.                         (D_t)
```

Substitution telescopes `(D_t)` exactly to `G'=0` [B]. This is the
polynomial replacement for the Laurent upward transition at the flipped
place.

For comparison, define bottom-up polynomial auxiliaries

```text
g_1=T h_0,
g_(l+1)=T(E^3 g_l+u^l h_l),
E^3g_7+u^7h_7=G'.                                      (U)
```

At every q-root, `T` is a unit, so `(U)` is equivalently

```text
t^v g_(l+1)=E^3g_l+u^l h_l,  v=-w,
```

in the q-local ring. This is why the old q-local transitions and terminal
laws survive verbatim.

There is also a descending normalization at a q-root `p` with `b=v_p(E)>0`.
Put

```text
E=p^b E_0,  u=p U_0,  s=3b-1,  A=T p^s.
```

Since `f+b(21-3f)=7+(7-f)(3b-1)`, division by the common `p^7` changes `G'`
into

```text
sum A^(7-f) U_0^f E_0^(21-3f) h_f.
```

Consequently its first three exact lines are

```text
h_7                              = A r_6,
E_0^3 h_6 + U_0 r_6              = A r_5,
E_0^6 h_5 + U_0 r_5              = A r_4.               (D_p)
```

The factor `T` in `A` is a p-unit but is retained because this is an exact
identity, not merely an order equality [B]. Roots with `b=0` use the upward
q-unit form `(U)`.

## 2. Source formulas and local cone

After `sigma=4d0-d2^2`, `f31_graded.txt` gives [A]

```text
h_7 = 8192 d1^2,
h_6 = 14336 d1^2 d2 + 8192 d1 e - 3072 sigma^2,
h_5 = -12288 d1^2 d2^2 + 32256 d1^2 sigma
      +18432 d1 d2 e - 9216 d2 sigma^2 + 2048 e^2.       (H)
```

At an active place put `s=w` at `t` or `s=3b-1` at `p`, `m=a` at `t` or
`m=b` at `p`, and `x=v(d1), z=v(sigma), k=v(d2)`. For T1:

```text
v(r_6)=2x-s,
orders(h_6)={2x+k, x+m, 2z},
orders(h_5)={2x+2k, 2x+z, x+k+m, k+2z, 2m}.             (O1)
```

The verifier parses these monomials from the source, enumerates every
`x<=9`, `z<=12`, `k<=6` (also `sigma=0` and `d2=0`), and permits every tied
minimum to rise arbitrarily. The complete projected T1 possibilities through
`h_5` are [C]:

| place | finite `sigma`: allowed `(x,z)` | `sigma=0`: allowed `x` |
|:--|:--|:--|
| `t`, `a=11` | `5<=x<=9`, `3<=z<=12` | `5<=x<=9` |
| `t`, `a=12` | `(3,0),(4,1),...,(8,5)`, or `x=9, 6<=z<=12` | `x=9` |
| `t`, `a=13` | none | none |
| `t`, `a=14` | `(6,0),(7,1),(8,2),(9,3)` | none |
| q-place `b=1` | `(1,0),(2,1)`, or `3<=x<=9, 2<=z<=12` | `3<=x<=9` |
| q-place `b=2` | `x=7, 5<=z<=12` | `x=7` |
| q-place `b=3` | `(4,0),(5,1),...,(9,5)` | none |
| q-place `b=4` | none | none |

There is no additional restriction on `k` at this order-only depth: every
listed pair survives the enlarged cone for all `k=0,...,6` and for `d2=0`.

For T2, `d1=0`, so

```text
v(r_5)=2z-s,   orders(h_5)={k+2z,2m}.                   (O2)
```

The complete possibilities are exactly the first-round T2 conditions [C]:

```text
t: a=11,12,13,14 => z>=3,6,9,12 respectively;
q: b=1 => z>=2; b=2 => none; b=3 => z=7; b=4 => none.
```

For `b=3`, survival includes a genuine depth-two residue cancellation between
`U_0 r_5` and `2048 E_0^6 e^2`; it is not automatic divisibility.

## 3. Six new kills

Orders at distinct linear places add toward polynomial degree. Combining the
table above with `deg d1<=9` and `deg sigma<=12` gives these empty cones [D].

| `a` | `b` | branch | exact obstruction |
|---:|:---:|:---:|:---|
| 11 | `(2,0,0,0)` | T1 | t requires `x>=5`; the `b=2` root requires `x=7`, so `deg d1>=12` |
| 11 | `(2,1,0,0)` | T1 | the same two places already force `deg d1>=12` |
| 11 | `(3,1,0,0)` | T1 | t, `b=3`, `b=1` force `deg d1>=5+4+1=10` |
| 12 | `(2,0,0,0)` | T1 | t requires `x>=3`; the `b=2` root requires `x=7`, so `deg d1>=10` |
| 12 | `(2,1,0,0)` | T1 | the same two places already force `deg d1>=10` |
| 13 | `(0,0,0,0)` | T1 | the t-local `h_6/h_5` cone is empty for every `x<=9,z<=12,k<=6` |

No coefficient cancellation is assumed impossible here. In particular, the
`a=13` kill remains valid after all tied minima in both `h_6` and `h_5` are
allowed arbitrary rise.

## 4. Per-input-branch verdict

`K` is killed this round. `O` is open with the residual shape in section 5.
Only the 33 branches supplied by `ALT_REGIME.md` appear.

| `a` | sorted `b` | T1 | T2 |
|---:|:---:|:---:|:---:|
| 11 | `(0,0,0,0)` | O | O |
| 11 | `(1,0,0,0)` | O | O |
| 11 | `(1,1,0,0)` | O | O |
| 11 | `(1,1,1,0)` | O | O |
| 11 | `(1,1,1,1)` | O | O |
| 11 | `(2,0,0,0)` | **K** | -- |
| 11 | `(2,1,0,0)` | **K** | -- |
| 11 | `(3,0,0,0)` | O | O |
| 11 | `(3,1,0,0)` | **K** | O |
| 12 | `(0,0,0,0)` | O | O |
| 12 | `(1,0,0,0)` | O | O |
| 12 | `(1,1,0,0)` | O | O |
| 12 | `(1,1,1,0)` | O | O |
| 12 | `(2,0,0,0)` | **K** | -- |
| 12 | `(2,1,0,0)` | **K** | -- |
| 12 | `(3,0,0,0)` | O | -- |
| 13 | `(0,0,0,0)` | **K** | O |
| 13 | `(1,0,0,0)` | -- | O |
| 14 | `(0,0,0,0)` | O | O |
| 14 | `(1,0,0,0)` | O | -- |

Dashes are branches killed before this round, not part of the 33-branch
denominator.

## 5. Exact residual shapes of the 27 survivors

Work over a splitting field and write `q` (up to a unit) as `p1 p2 p3 p4`.
For every row below set

```text
R_b=product p_i^b_i,
E=R_b F,   e=t^a R_b F,
deg F <= r=15-a-sum(b),   gcd(F,tq)=1.                  (R0)
```

### T1: 13 residual branches

Choose one allowed local pair `(x_0,z_0)` from the `t,a` row of section 2 and
one `(x_i,z_i)` from the `b_i` row for every `b_i>0`, subject to

```text
X=x_0+sum x_i <=9,   Z=z_0+sum z_i <=12.                (R1)
```

Alternatively `sigma=0` is allowed only when every involved place permits
it in the last column of section 2. Then

```text
d1=t^x_0 product p_i^x_i D,       deg D<=9-X,
sigma=t^z_0 product p_i^z_i S,    deg S<=12-Z,           (R2)
```

with displayed orders exact and `D,S` coprime to the displayed factors. The
exact congruences `(D_t)` and `(D_p)` through `h_5` remain imposed.
