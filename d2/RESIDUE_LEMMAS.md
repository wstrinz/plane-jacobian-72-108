# Shared residue lemmas for P6/P10/P11

This file proves the level-6, level-5, and level-4 initial-form lemmas once
for both `cascade_cones_qt.json` (subcase (2), abbreviated **sub2**) and
`cascade_cones_sub1_qt.json` (subcase (1), abbreviated **sub1**).  The
coefficients below are obtained by rewriting the parsed `f31_graded.txt`
tables with `sigma = 4*d0-d2^2`; `residue_lemmas_verify.py` repeats that
derivation without copying any coefficient.

The conclusion is deliberately local.  A **KILL** says that the indicated
residue obligation cannot occur over the splitting field of the fixed
quartic `q`.  A **CONSTRAINT** says exactly that the displayed hypersurface
is necessary; it does not claim that the rest of the cascade or an
`identical_vanishing` polynomial identity is satisfied.

## 1. Initial-form lemma, all depths and both windows

Fix a root-place `p` and write

```
(k,x,z,b) = (v_p(d2),v_p(d1),v_p(sigma),v_p(e)),
(D,X,S,E) = (c_d2,c_d1,c_sigma,c_e).
```

Every displayed leading coefficient is nonzero.  If `T` is the tied set,
all its monomials have one common order `m`.  Writing, for example,
`d2=pi^k(D+D_1 pi+...)`, the condition
`v_p(sum(T)) >= m+delta` is equivalent to the first `delta` convolution
coefficients being zero.  At depth one this is simply

```
             sum_{A_alpha M_alpha in T} A_alpha D^alpha_d2
                 X^alpha_d1 S^alpha_sigma E^alpha_e = 0.       (IF)
```

For depth `delta >= 2`, (IF) is the coefficient numbered zero and, for
`1 <= j < delta`, the same sum is taken over all choices of jet indices in
each factor whose total is `j`.  Thus every deeper obligation contains the
depth-one equation proved here.  The proof is coefficient extraction after
factoring `pi^m`; it is independent of the window, of `a`, and of whether
the place is labelled `q` or `t`.  An `identical_vanishing` record implies
the same necessary initial equation, but generally imposes more than it.

For the three full ties, the source gives the promised shared equations:

```
P6:   14336 X^2 D + 8192 X E - 3072 S^2 = 0,
      2x+k = x+b = 2z.

P10: -12288 X^2 D^2 + 32256 X^2 S + 18432 X D E
      -9216 D S^2 + 2048 E^2 = 0,
      2x+2k = 2x+z = x+k+b = k+2z = 2b.

P11: -220752 X^4 -31232 X^2 D^3 -23616 X^2 D S
      -3072 X D^2 E +34560 X E S -5184 D^2 S^2
      +5632 D E^2 -12096 S^3 = 0,
      4x = 2x+3k = 2x+k+z = x+2k+b = x+b+z
         = 2k+2z = k+2b = 3z.
```

The P10/P11 equalities solve as `(k,x,z,b)=(2r,3r,4r,5r)`.  At a `b=0`
q-place, `r=0`.  For P6, `b=x+k`; nonnegative valuations and `b=0` again
give `k=x=z=0`.  None of P6/P10/P11 is a kill: rational torus points are,
respectively,

```
(D,X,S,E)=(-2,-1,-2,-5), (-6,-2,-4,6), (-54,-28,100,3730).
```

Hence the shared backbone lemmas are **CONSTRAINT lemmas**, with the exact
hypersurfaces P6, P10, and P11.

## 2. Coherence with sigma = 4*d0-d2^2

Put `w=v_p(d0)` and let `C0` be the leading coefficient of `d0`.

- If `w<2k`, then `z=w` and `S=4C0`.
- If `w>2k`, then `z=2k` and `S=-D^2`.
- If `w=2k` and `4C0-D^2 != 0`, then `z=2k` and
  `S=4C0-D^2`.
- If `z>2k`, necessarily `w=2k` and `4C0=D^2`; the first nonzero later
  convolution is `S`.  In particular `z<min(w,2k)` is impossible.
- If `sigma` is identically zero, `d0=d2^2/4`.  If `d2` is identically
  zero, `sigma=4d0` and the first bullet applies at every finite order.

These conditions are additional to every row below whenever the recorded
orders interact in the stated way.  They never weaken a KILL.

## 3. Global zero-flag cuts

These are all global flag combinations that actually occur in both survivor
files.  `1-term KILL` means that a requested rise would make a nonzero
leading coefficient vanish.  The T2 level-6 entry is the important pure
square: `-3072 S^2=0` forces `S=0`, contradicting its being leading.

| flags | h6 | h5 | h4 |
|:--|:--|:--|:--|
| T1, none | P6 = C04, CONSTRAINT | P10 = C10, CONSTRAINT | P11 = C23, CONSTRAINT |
| T1, `d2=0` | C01, CONSTRAINT | C05, CONSTRAINT | C18, CONSTRAINT |
| T1, `sigma=0` | C03, CONSTRAINT | C08, **KILL over the q-splitting field** | C21, CONSTRAINT |
| T1, `d2=sigma=0` | `8192 X E`, 1-term KILL | `2048 E^2`, 1-term KILL | `-220752 X^4`, 1-term KILL |
| T2 (`d1=0`) | `-3072 S^2`, 1-term KILL | C06, CONSTRAINT | C19, CONSTRAINT |
| T2 and `d2=0` | `-3072 S^2`, 1-term KILL | `2048 E^2`, 1-term KILL | `-12096 S^3`, 1-term KILL |

The one-term rows carry no survivor `monomial_tie_rise`: the engine already
forbids a rise of a unique minimum.  They therefore kill the hypothetical
flag/valuation branch rather than an additional listed survivor obligation.

## 4. Every tied-set shape in the two worklists

The inventories contain 23 distinct nonempty tied supports at levels 4--6.
The following notation records their valuation constraints without hiding
any coincidence:

```
level 6: a=2x+k, b=x+b, c=2z
level 5: a=2x+2k, b=2x+z, c=x+k+b, d=k+2z, e=2b
level 4: a=4x, b=2x+3k, c=2x+k+z, d=x+2k+b,
         e=x+b+z, f=2k+2z, g=k+2b, h=3z.
```

In the table, `[a,c]` means precisely that those displayed weights are
equal.  A common nonzero rational factor and any common nonzero leading
monomial have been divided out; consequently each equation is equivalent
to the source initial equation on the coefficient torus.

| ID | L | tied weights | primitive depth-one relation | result |
|:--|--:|:--|:--|:--|
| C01 | 6 | `[b,c]` | `8XE-3S^2=0` | CONSTRAINT |
| C02 | 6 | `[a,c]` | `14X^2D-3S^2=0` | CONSTRAINT |
| C03 | 6 | `[a,b]` | `7XD+4E=0` | CONSTRAINT |
| C04 | 6 | `[a,b,c]` | `14X^2D+8XE-3S^2=0` (P6/1024) | CONSTRAINT |
| C05 | 5 | `[b,e]` | `63X^2S+4E^2=0` | CONSTRAINT |
| C06 | 5 | `[d,e]` | `-9DS^2+2E^2=0` | CONSTRAINT |
| C07 | 5 | `[b,d,e]` | `63X^2S-18DS^2+4E^2=0` | CONSTRAINT |
| C08 | 5 | `[a,c,e]` | `6X^2D^2-9XDE-E^2=0` | **KILL over Q and the q-splitting field** |
| C09 | 5 | `[a,b,d]` | `8X^2D^2-21X^2S+6DS^2=0` | CONSTRAINT |
| C10 | 5 | `[a,b,c,d,e]` | `24X^2D^2-63X^2S-36XDE+18DS^2-4E^2=0` (P10/-512) | CONSTRAINT |
| C11 | 4 | `[e,h]` | `20XE-7S^2=0` | CONSTRAINT |
| C12 | 4 | `[a,h]` | `73X^4+4S^3=0` | CONSTRAINT |
| C13 | 4 | `[g,h]` | `88DE^2-189S^3=0` | CONSTRAINT |
| C14 | 4 | `[f,h]` | `3D^2+7S=0` | CONSTRAINT |
| C15 | 4 | `[a,e]` | `511X^3-80ES=0` | CONSTRAINT |
| C16 | 4 | `[a,b]` | `13797X^2+1952D^3=0` | CONSTRAINT |
| C17 | 4 | `[f,g]` | `81DS^2-88E^2=0` | CONSTRAINT |
| C18 | 4 | `[a,e,h]` | `511X^4-80XES+28S^3=0` | CONSTRAINT |
| C19 | 4 | `[f,g,h]` | `81D^2S^2-88DE^2+189S^3=0` | CONSTRAINT |
| C20 | 4 | `[b,d,g]` | `61X^2D^2+6XDE-11E^2=0` | **KILL over Q and the q-splitting field** |
| C21 | 4 | `[a,b,d,g]` | `13797X^4+1952X^2D^3+192XD^2E-352DE^2=0` | CONSTRAINT |
| C22 | 4 | `[a,b,c,f,h]` | `13797X^4+1952X^2D^3+1476X^2DS+324D^2S^2+756S^3=0` | CONSTRAINT |
| C23 | 4 | `[a,b,c,d,e,f,g,h]` | `13797X^4+1952X^2D^3+1476X^2DS+192XD^2E-2160XES+324D^2S^2-352DE^2+756S^3=0` (P11/-16) | CONSTRAINT |

### Proof of the classifications

For C08, divide by `(XD)^2` and put `r=E/(XD)`.  The relation becomes
`r^2+9r-6=0`, whose discriminant has square class `105`.  For C20 the same
substitution gives `11r^2-6r-61=0`, whose discriminant has square class
`170`.  The fixed quartic `q` is irreducible and has Galois group `S4`; its
discriminant has square class `17`.  Since `S4` has abelianization `C2`, its
splitting field has the unique quadratic subfield `Q(sqrt(17))`.  It
contains neither `sqrt(105)` nor `sqrt(170)`.  Thus C08 and C20 have no
torus point over the splitting field.  They do have real torus points,
using respectively

```
r=(-9+sqrt(105))/2,        r=(3+2sqrt(170))/11.
```

So the obstruction is arithmetic, not a sign/no-real-point obstruction.
Every other multi-term row has a rational torus point.  Direct witnesses
for `(D,X,S,E)`, with unused slots omitted, are:

| rows | rational witnesses |
|:--|:--|
| C01, C02, C03, C04 | `(-,-6,-4,-1)`; `(3/14,1,1,-)`; `(1,1,-,-7/4)`; `(-2,-1,-2,-5)` |
| C05, C06, C07 | `(-,1,-4/63,1)`; `(2,-,-2,-6)`; `(-6,-2,-1,-6)` |
| C09, C10 | `(13/6,13/6,169/36,-)`; `(-6,-2,-4,6)` |
| C11, C12 | `(-,1,10,35)`; `(-,(73/4)^2,-(73/4)^3,-)` |
| C13, C14 | `(189/88,-,1,1)`; `(7,-,-21,-)` |
| C15, C16 | `(-,1,1,511/80)`; `(-13797/1952,13797/1952,-,-)` |
| C17, C18, C19 | `(88/81,-,1,1)`; `(-,1,1,539/80)`; `(1,-,1/27,1/27)` |
| C21 | with `u=-13797/1792`: `(u,u,-,u^2)` |
| C22 | with `u=152/511`: `(1/u,1/u,-4/(3u^2),-)` |
| C23 | `(-54,-28,100,3730)` |

Substitution proves that each remaining zero set is a nonempty, nontrivial
hypersurface.  This also verifies that none of the multi-term cuts is killed
over the reals.  The checker performs all substitutions into equations
reconstructed from the source, not into hand-entered versions of this table.

## 5. Per-cell applicability

Cell lists below are complete but compressed. `1000T1@6-9` means the four
cells `a6 b=1000 T1` through `a9 b=1000 T1`; commas denote separated values.
For each ID, occurrences and cells are the union of its rise and
identical-vanishing uses at both recorded place kinds. A CONSTRAINT row adds
that exact equation to each carrying case. A KILL row removes every carrying
case (not automatically every cell, since a cell can have other cases).

| ID | sub2 incidence | sub1 incidence |
|:--|:--|:--|
| C01 | 36 occ/11 cells: `0000T1@8,10; 1000T1@6-8; 1100T1@6-8; 1110T1@6-7; 1111T1@6` | 219 occ/42 cells: `0000T1@2-10; 1000T1@2-10; 1100T1@2-10; 1110T1@2-10; 1111T1@2-7` |
| C02 | 48 occ/14 cells: `0000T1@8-9; 1000T1@6-9; 1100T1@6-8; 1110T1@5-7; 1111T1@6; 3000T1@7` | 272 occ/80 cells: `0000T1@2-9; 1000T1@2-9; 1100T1@2-9; 1110T1@2-9; 1111T1@2-9; 3000T1@2-9; 3100T1@2-9; 3110T1@2-9; 3111T1@2-9; 5000T1@2-9` |
| C03 | 15 occ/3 cells: `0000T1@9-10; 1000T1@9` | 212 occ/35 cells: `0000T1@2-10; 1000T1@2-10; 1100T1@2-10; 1110T1@2-6,8-10` |
| C04/P6 | 132 occ/14 cells: `0000T1@8-10; 1000T1@6-9; 1100T1@6-8; 1110T1@5-7; 3000T1@7` | 621 occ/72 cells: `0000T1@2-10; 1000T1@2-10; 1100T1@2-10; 1110T1@2-10; 3000T1@2-10; 3100T1@2-10; 3110T1@2-10; 5000T1@2-10` |
| C05 | 24 occ/6 cells: `0000T1@8-10; 1000T1@8-9; 1100T1@8` | 389 occ/109 cells: `0000T1@0-10; 1000T1@0-10; 1100T1@0-10; 1110T1@0-10; 1111T1@0-6,8; 2000T1@0-5,9-10; 2100T1@0-4,9-10; 2110T1@0-3,9-10; 3000T1@0-10; 3100T1@0-7,9-10; 3110T1@0-6,9-10; 3111T1@0-5` |
| C06 | 24 occ/10 cells: `1000T1@7-9; 1000T2@9; 1100T1@7-8; 1110T1@6-7; 1111T1@6; 3000T1@7` | 614 occ/166 cells: `0000T2@0-10; 1000T1@2-10; 1000T2@0-10; 1100T1@2-10; 1100T2@0-10; 1110T1@2-10; 1110T2@0-10; 1111T1@2-10; 1111T2@0-7; 3000T1@2-10; 3000T2@0-10; 3100T1@2-10; 3100T2@0-6,8-10; 3110T1@2-10; 3110T2@0-4,8-10; 3111T1@2-9; 5000T1@2-10; 5000T2@0-2,9-10` |
| C07 | 16 occ/10 cells: `1000T1@6-9; 1100T1@6-8; 1110T1@6-7; 1111T1@6` | 108 occ/54 cells: `1000T1@2-10; 1100T1@2-10; 1110T1@2-10; 1111T1@2-10; 3100T1@2-10; 3110T1@2-10` |
| **C08 KILL** | 15 occ/3 cells: `0000T1@9-10; 1000T1@9` | 304 occ/54 cells: `0000T1@0-10; 1000T1@0-10; 1100T1@0-10; 1110T1@0-6,8-10; 2000T1@0-10` |
| C09 | 57 occ/18 cells: `0000T1@8-10; 1000T1@6-9; 1100T1@6-8; 1110T1@5-7; 1111T1@6; 3000T1@6-7; 3100T1@6; 3110T1@5` | 410 occ/159 cells: `0000T1@1-10; 1000T1@1-10; 1100T1@1-10; 1110T1@1-10; 1111T1@1-10; 2000T1@1-9; 2100T1@1-9; 2110T1@1-9; 3000T1@1-10; 3100T1@1-10; 3110T1@1-10; 3111T1@1-9; 3300T1@1-9; 3310T1@1-8; 5000T1@1-9; 5100T1@1-9; 5110T1@1-8` |
| C10/P10 | 145 occ/17 cells: `0000T1@8-10; 1000T1@6-9; 1100T1@6-8; 1110T1@5-7; 3000T1@6-7; 3100T1@6; 3110T1@5` | 910 occ/161 cells: `0000T1@0-10; 1000T1@0-10; 1100T1@0-10; 1110T1@0-10; 1111T1@0; 2000T1@0-10; 2100T1@0-10; 2110T1@0-10; 3000T1@0-10; 3100T1@0-10; 3110T1@0-10; 3111T1@0; 3300T1@0-9; 3310T1@0-8; 5000T1@0-10; 5100T1@0-9; 5110T1@0-8` |
| C11 | none | 30 occ/23 cells: `0000T1@6-8; 1000T1@0,5-10; 1100T1@0,5-8,10; 1110T1@0,6-8,10; 1111T1@6-7` |
| C12 | 14 occ/14 cells: `0000T1@8-10; 1000T1@6-9; 1100T1@6-8; 1110T1@6-7; 1111T1@6; 3000T1@7` | 155 occ/155 cells: `0000T1@1-10; 1000T1@1-10; 1100T1@1-10; 1110T1@1-10; 1111T1@1-9; 2000T1@1-9; 2100T1@1-9; 2110T1@1-9; 3000T1@1-9; 3100T1@1-9; 3110T1@1-9; 3111T1@1-9; 3300T1@1-9; 3310T1@1-8; 5000T1@1-9; 5100T1@1-9; 5110T1@1-8` |
| C13 | none | 11 occ/11 cells: `1000T2@0-10` |
| C14 | 46 occ/27 cells: `0000T2@8; 1000T1@6-9; 1000T2@5-9; 1100T1@6-8; 1100T2@6-8; 1110T1@5-7; 1110T2@6-7; 1111T1@6; 3000T1@6-7; 3000T2@7; 3100T1@6; 3110T1@5` | 516 occ/219 cells: `0000T2@1-10; 1000T1@0-10; 1000T2@1-10; 1100T1@0-10; 1100T2@1-9; 1110T1@0-10; 1110T2@1-9; 1111T1@0-10; 1111T2@1-9; 3000T1@0-10; 3000T2@1-9; 3100T1@0-10; 3100T2@1-9; 3110T1@0-10; 3110T2@1-9; 3111T1@0-9; 3300T1@0-9; 3310T1@0-8; 5000T1@0-10; 5000T2@1-9; 5100T1@0-9; 5110T1@0-8` |
| C15 | none | 2 occ/2 cells: `0000T1@9; 1000T1@9` |
| C16 | 7 occ/7 cells: `0000T1@8-10; 1000T1@7-9; 1100T1@8` | 46 occ/46 cells: `0000T1@1-10; 1000T1@1-9; 1100T1@1-9; 1110T1@1-9; 2000T1@1-9` |
| C17 | 10 occ/6 cells: `1000T1@7-9; 1100T1@7-8; 1110T1@7` | 132 occ/69 cells: `0000T2@1-4,6-8; 1000T1@2-10; 1000T2@1-4,6-9; 1100T1@2-10; 1110T1@2-10; 1111T1@2-10; 3000T1@2-10; 3100T1@2-10` |
| C18 | 43 occ/13 cells: `0000T1@8-10; 1000T1@6-9; 1100T1@6-8; 1110T1@6-7; 3000T1@7` | 453 occ/162 cells: `0000T1@0-10; 1000T1@0-10; 1100T1@0-10; 1110T1@0-10; 1111T1@0,5; 2000T1@0-10; 2100T1@0-10; 2110T1@0-10; 3000T1@0-10; 3100T1@0-10; 3110T1@0-10; 3111T1@0; 3300T1@0-9; 3310T1@0-8; 5000T1@0-10; 5100T1@0-9; 5110T1@0-8` |
| C19 | 30 occ/12 cells: `0000T2@8; 1000T2@5-9; 1100T2@6-8; 1110T2@6-7; 3000T2@7` | 299 occ/89 cells: `0000T2@0-10; 1000T2@0-10; 1100T2@0-10; 1110T2@0-10; 1111T2@0; 3000T2@0-10; 3100T2@0-10; 3110T2@0-10; 5000T2@0-10` |
| **C20 KILL** | none | 17 occ/8 cells: `0000T1@1-4,6-9` |
| C21 | 31 occ/7 cells: `0000T1@8-10; 1000T1@7-9; 1100T1@8` | 248 occ/55 cells: `0000T1@0-10; 1000T1@0-10; 1100T1@0-10; 1110T1@0-10; 2000T1@0-10` |
| C22 | 57 occ/18 cells: `0000T1@8-10; 1000T1@6-9; 1100T1@6-8; 1110T1@5-7; 1111T1@6; 3000T1@6-7; 3100T1@6; 3110T1@5` | 419 occ/160 cells: `0000T1@1-10; 1000T1@1-10; 1100T1@1-10; 1110T1@1-10; 1111T1@1-10; 2000T1@1-9; 2100T1@1-9; 2110T1@1-9; 3000T1@1-10; 3100T1@1-10; 3110T1@1-10; 3111T1@1-9; 3300T1@1-9; 3310T1@1-8; 5000T1@1-10; 5100T1@1-9; 5110T1@1-8` |
| C23/P11 | 145 occ/17 cells: `0000T1@8-10; 1000T1@6-9; 1100T1@6-8; 1110T1@5-7; 3000T1@6-7; 3100T1@6; 3110T1@5` | 937 occ/161 cells: `0000T1@0-10; 1000T1@0-10; 1100T1@0-10; 1110T1@0-10; 1111T1@0; 2000T1@0-10; 2100T1@0-10; 2110T1@0-10; 3000T1@0-10; 3100T1@0-10; 3110T1@0-10; 3111T1@0; 3300T1@0-9; 3310T1@0-8; 5000T1@0-10; 5100T1@0-9; 5110T1@0-8` |

### 5.1 The six genuinely-new sub1 keys

The worklist's “six new residue polynomials” are six new pattern keys but
only four distinct tied supports.  Their exact instantiations are:

| new key | lemma | exact cells | effect |
|:--|:--|:--|:--|
| q L4 identical, `{-12096 S^3,34560 XSE}` (16 occ) | C11 | `1000T1@0,5,9-10; 1100T1@0,5,10; 1110T1@0,10` | CONSTRAINT |
| t L4 rise, same support (14 occ) | C11 | `0000T1@6-8; 1000T1@6-8; 1100T1@6-8; 1110T1@6-8; 1111T1@6-7` | CONSTRAINT |
| q L4 identical, `{-12096 S^3,5632 D E^2}` (11 occ) | C13 | `1000T2@0-10` | CONSTRAINT |
| t L4 identical, `{-220752 X^4,34560 XSE}` (2 occ) | C15 | `0000T1@9; 1000T1@9` | CONSTRAINT |
| t L4 identical, `{-31232 D^3X^2,-3072 D^2XE,5632 DE^2}` (13 occ) | C20 | `0000T1@1-4,6-9` | **KILL** those cases |
| t L4 rise, same C20 support (4 occ) | C20 | `0000T1@4,6-8` | **KILL** those cases |

Thus the six keys contribute three exact Phase-D hypersurface constraints
and one arithmetic kill lemma (used by two keys).  The q-place equations in
sub1 are independent of `a`; the displayed `@` ranges show that the same
lemma is being instantiated, not reproved, as `a` varies.

## 6. Impact summary

- **KILL:** C08 kills every carrying flag case in 3 sub2 cells (15
  obligations) and 54 sub1 cells (304 obligations).  C20 kills its carrying
  cases in 8 sub1 cells (17 obligations).  Singleton initial forms kill any
  hypothetical requested rise, notably T2 level 6 with `-3072 S^2`, but no
  survivor tie obligation records those already-rejected rises.
- **CONSTRAINT:** C01--C07, C09--C19, and C21--C23.  In particular P6, P10,
  and P11 are constraints, not local kills.  Their q-rise cell impacts are
  respectively 14/17/17 cells in sub2 and 72/159/159 cells in sub1.
- No multi-term shape is sign-killed over the reals.  C08 and C20 have real
  points but no points over Q or over the q-splitting field; all other
  multi-term shapes already have rational points.

Run `python residue_lemmas_verify.py` from this directory to reparse the
source, re-extract all 41/67 patterns and 23 tied supports, verify every
rational or real point and both splitting-field obstructions, check every
occurring global zero-flag cut, and audit the headline incidence counts.
