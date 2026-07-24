# The period-12 window functions for (75,125), derived exactly

## Verdict

**DERIVED.** The arithmetic layer a (75,125) window compiler needs — the
period-12 window-cap functions, their residue-class interaction table, and the
consistency with the built `Phi` — is derived exactly and closed-form. The
quasipolynomial caps that `G_SYSTEM_75_125.md` flagged as *obstructed* are pinned
to a single, clearly-stated structural premise inherited verbatim from `(72,108)`
("`Phi` realises the extreme ray of the window cone"), under which:

```
lower y-order cap :  L(w) = floor((alpha*w + beta_m)/q) = ceil(alpha*w/q)     m = w mod q
upper y-degree cap:  U(w) = deg_slope * w                       (affine)

    alpha     = 10a^2 - 8a + 1   =  67        (a=3 ;  a=2: 25)
    q         = 5a - 3           =  12        (= q_window ; a=2: 7)
    deg_slope = 5a - 1           =  14        (a=3 ;  a=2: 9)
    beta_m    = (-alpha*m) mod q :  [0,5,10,3,8,1,6,11,4,9,2,7]   (m = 0..11)
```

Only the **lower** (y-order) cap is quasipolynomial — because `W_step = 67/12`
is non-integral; the **upper** (y-degree) cap is affine because `deg_slope =
5a-1` is an integer for every `a`. This is the exact content of "quasipolynomial
window cap". `window_functions_75_125.py` is the derivation;
`window_functions_75_125_verify.py` is the independent exact checker (**37 checks,
`--quiet`, exit 0**), whose controls reproduce the a=2 period-7 analogue (against
`f2_tower.py`'s window table) and the `(72,108)` affine limit (against
`WINDOW_CAPS.md`).

---

## 1. The floor/ceiling window-cap functions

### 1.1 Where the caps come from

The `(72,108)` window variable `d_{4-k} = C_{4-k}·C4^(2k-1)` obeyed the affine
caps `ord >= 12k`, `deg <= 15k` (sub1) / `14k` (sub2), `k` = u-slice weight
(`WINDOW_CAPS.md`). Both caps have **integer** slope there because
`W_step = ord_y(Phi)/M = 204/17 = 12` is an integer: the physical y-order grading
is the u-grading scaled by 12.

At `(75,125)`, `W_step = ord_y(Phi)/M = 201/36 = 67/12` is **non-integral**
(`G_SYSTEM_75_125.md` §3). A window object of intrinsic u-slice weight `w` has
physical y-order that tracks the ray `ord = W_step * w = (67/12) w`, an integer
only when `12 | w`. The tight integer lower cap is therefore the ceiling of that
ray, which is exactly a quasipolynomial of period `q = 12`:

```
L(w) = ceil( (67/12) w ) = floor( (67 w + beta_m)/12 ),   beta_m = (-67 m) mod 12.
```

Writing `w = 12n + m` gives `L(w) = 67n + L(m)` with the twelve base values
`L(0..11) = 0,6,12,17,23,28,34,40,45,51,56,62`, i.e. the offsets
`beta_m = 12·L(m) - 67m = (-67m) mod 12 = (5m) mod 12`:

| m | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |10 |11 |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| beta_m | 0 | 5 |10 | 3 | 8 | 1 | 6 |11 | 4 | 9 | 2 | 7 |

### 1.2 The upper cap is affine

`deg_slope = deg_y(Phi)/M = 504/36 = 14`, an integer, so `U(w) = 14w` with no
residue offset — the direct analogue of `(72,108)`'s clean affine `deg <= 14k`
(the sub2 cap that `Phi` sits on). The general-`a` slope is `5a-1`, integral for
every `a`; the deg cap is **never** quasipolynomial. The whole quasipolynomial
phenomenon lives on the **lower** (y-order) cap alone.

### 1.3 The family formulas (F2 family, a = j+2)

From `Phi_a = -(1/(3a)) y^(30a^2-24a+3) (y^3+1)^(15a^2-12a+2)` at forcing slice
`M = bt + jphi = 15a-9 = 3(5a-3)` (`f2_family_verify.py`, `C_SERIES_75_125.md`):

```
ord_y(Phi_a) = 30a^2 - 24a + 3          W_step  = ord/M = (10a^2-8a+1)/(5a-3)
deg_y(Phi_a) = 3(5a-1)(5a-3)            deg/M   = 5a - 1        (integral)
```

so `alpha = 10a^2-8a+1`, `q = 5a-3`, `deg_slope = 5a-1`. These give `(67,12,14)`
at `a=3` and `(25,7,9)` at `a=2`.

---

## 2. The class-interaction table

Each window symbol carries a residue class `m = w mod 12` of its u-weight
(`w(d_m) = t-m`, `w(Phi) = M`). **Multiplication adds u-weights**, so classes
compose additively:

```
class(x · y) = (class x + class y) mod 12.
```

The 12×12 composition structure is the cyclic group `Z/12`. Occupancy:

| symbols | u-weights | classes mod 12 |
|---|---|---|
| state `d3,d2,d1,d0,e=dm1` | 2,3,4,5,6 | {2,3,4,5,6} |
| spares `dm2..dm10` = `d_-2..d_-10` | 7..15 | {7,8,9,10,11,0,1,2,3} |
| generators `G1..G9,G11` | 26..34, 36 | **{0,2,3,4,5,6,7,8,9,10}** (all but {1,11}) |
| skipped `G10` | 35 | 11 |
| `Phi` | 36 | **0** |

The forcing window is `S^5`, so every generator monomial is a product of 5
window factors whose classes sum mod 12 to the generator's class. The `beta`
offsets are a **group 1-cocycle** for this composition:

```
beta_{m1} + beta_{m2} = beta_{(m1+m2) mod 12} + 12 · carry(m1,m2),   carry in {0,1},
```

and the carry is *exactly* the ceiling superadditivity defect of the lower cap,
`carry = L(w1) + L(w2) - L(w1+w2) in {0,1}`. So the class table and the cap
arithmetic are one object: composing two window slices either preserves the ray
value (carry 0) or gains one unit of y-order (carry 1) precisely when their
fractional parts overflow.

---

## 3. Consistency checks

**(a) The known `Phi` point sits AT both caps, at equality.** `M = 36 = 3·12` is
a multiple of the period, so the quasipolynomial floor is exact there:
`L(36) = 67·36/12 = 201 = ord_y(Phi)` and `U(36) = 14·36 = 504 = deg_y(Phi)`.
The stripped `Phi` degree `U(36) - L(36) = 303` is the degree of the
`(y^3+1)^101` cofactor block — `Phi` sits on the lower cap in y-order and the
upper cap in y-degree, exactly as `(72,108)`'s tight rows do.

**(b) Reduction to a=2 (the control against `f2_tower.py`).** Under the family
substitution `a=2`, the caps become `alpha=25, q=7, deg_slope=9`,
`beta_m = (-25m) mod 7 = [0,3,6,2,5,1,4]`; `Phi_2` (signature `deg 189, ord 75`
at `M=21=3·7`) sits at `L(21)=75, U(21)=189`. The naive physical-order
fractional denominators of the forcing slices are `{1,7}` at `a=2` and
`{1,2,3,4,6,12}` at `a=3` — reproducing `f2_tower.py` §G / `F2_TOWER.md` §2b
verbatim (these are the divisors of the period that occur). A second control:
the `(72,108)` **integral limit** `W_step=12` gives `q_window=1`, `beta=[0]`,
and the caps degenerate to the affine `ord >= 12w`, `deg <= 14w` recited in
`WINDOW_CAPS.md` — the framework contains the landed affine case as its `q=1`
specialisation.

**(c) The `q_window = 5a-3` law.** `q = denom(W_step) = 5a-3` with
`gcd(alpha, q) = gcd(10a^2-8a+1, 5a-3) = 1` for all `a` (checked `a=2..6`), so
`W_step` is already in lowest terms and `gcd(q(2),q(3)) = gcd(7,12) = 1` — the
incommensurate period jump of `F2_TOWER.md`.

All of the above are exact and independently re-derived in
`window_functions_75_125_verify.py` (37 checks, exit 0), including the landed
anchors (`M`, `ord_y(Phi)`, `deg_y(Phi)` at `a=2,3` and `(72,108)`) so the
derivation is cross-checked against published signatures, not self-referential.

---

## 4. Files

- `WINDOW_FUNCTIONS_75_125.md` — this writeup.
- `window_functions_75_125.py` — the derivation: family constants, the
  `window_law` deriving `(alpha, q, beta, deg_slope)` from any case's `Phi`
  signature, the cap functions `L,U`, and the class/composition machinery.
  Run end-to-end for the (75,125) instance.
- `window_functions_75_125_verify.py` — independent exact checker (`--quiet`,
  exit 0): §1 cap functions, §2 class-interaction table, §3 consistency incl.
  the a=2 period-7 control and the (72,108) affine control.

---

## 5. The honest boundary (what is derived vs. what needs the bridge)

**Derived exactly, from the u-graded G-system + the built `Phi` + the family:**
the period `q = 12`, the numerator `alpha = 67`, the affine deg cap `14w`, the
quasi-linearity `L(w+12) = L(w)+67`, and the class-0 line (`Phi` and every
`w ≡ 0 mod 12` slice at `ord = 67w/12` exactly). The `q_window = 5a-3` law and
`gcd(alpha,q)=1`.

**Derived under one inherited premise:** the twelve `beta_m`. They equal
`(-alpha m) mod q` **iff** `Phi` realises the *extreme* (minimal ord/weight) ray
of the window cone — i.e. `ord >= (alpha/q) w` holds for every window object,
with equality on `Phi`. This is exactly the sense in which `(72,108)`'s `Phi`
"sits at its caps" (`ord 204 = 12·17`), so the premise is structural, not new;
under it `L(w) = ceil(alpha w/q)` is forced by integrality and the `beta_m` are
pinned.

**The residual — the characterised obstruction.** The only thing *not* fixable
from the u-grading alone is whether the actual `(75,125)` window cone dips
*below* the `67/12` ray at some non-`Phi` weight `w`, which would raise the
corresponding `beta_m` (equivalently, deepen the tight lower cap). Settling that
needs the deeper Newton polygon of `P` — the **unreduced-polygon data**
(`C_SERIES_75_125.md` judgment 2) that no paper carries and that only the actual
bridge construction supplies (the `(72,108)` caps were proven from Prop 4.3's two
explicit `P`-polygons; the `(75,125)` analogue of that polygon is exactly the
missing input). So the boundary is precise: **the cap SHAPE — period, slopes,
deg cap, Phi line, and the canonical extreme-ray `beta_m` — is derived; whether
the true `beta_m` exceed the canonical values at some class is the one datum the
window compiler must read off the bridge.**

## 6. `[judgment]` list

1. **[inherited]** The corner `(5,20)→(7/5,2)`, `C = y^2(y^3+1)`,
   `Phi = -(1/9) y^201 (y^3+1)^101`, and the F2-family closed form are the landed
   phase-1 objects (`C_SERIES_75_125.md`, `f2_family_verify.py`), conditional on
   the standard-chart reduction (judgment 2 there). Unchanged here.
2. **[derived]** `alpha`, `q`, `deg_slope`, the period law, the class-interaction
   table, and the `q_window = 5a-3` law are computed exactly from the built `Phi`
   slice and the u-grading, and checked against the a=2 and (72,108) controls.
3. **[derived under the inherited extreme-ray premise]** The twelve `beta_m =
   (-alpha m) mod q`. The premise (`Phi` on the extreme ray) is the `(75,125)`
   image of `(72,108)`'s "`Phi` sits at its caps"; the residual (a true `beta_m`
   above canonical) is §5's characterised boundary, needing `P`'s unreduced
   polygon.
