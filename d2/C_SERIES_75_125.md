# The (75,125) C-series, built from the polygon data: N = 98 is DERIVED

## Verdict

**N = 98 CONFIRMED.**  Building the actual C-series and its D-transform tower
for the `(75,125)` case — from the corner/chain data, independently of the
N-formula — the tower length comes out `N = 98`, matching the formula value.
The forcing divisor that emerges is

```
Phi_75_125 = f * C^98 = -(1/9) y^201 (y^3+1)^101
           = -(1/9) y^201 (y+1)^101 (y^2-y+1)^101,
signature (deg, ord_y, mult_(y+1), cofactor) = (504, 201, 101, 202).
```

This is **THE TRANSFER TEST, phase 1**: the first from-scratch instantiation of
the `(72,108)`-style series/tower pipeline on a second case.  It discharges the
pending item in `PHI_75_125.md` judgment 3 (`N = 98` was FORMULA-based pending an
actual C-series build) and **upgrades that item to DERIVED**.  All numbers are
exact (sympy) and independently re-derived by `c_series_75_125_verify.py`
(**36/36 pass, exit 0**); the construction is `c_series_75_125.py`.

The build reproduces both landed checkpoints with the identical machinery:
`(72,108)` → `N=28`, `(238,204,30,4)` (STATE.md audited ground truth) and
`(108,144)` → `N=67`, `(550,205,69,276)` (corner-144).

---

## 1. What "build the C-series" means here, and what was actually built

The case compiler (`case_compiler.py`) emits the corner signature and reads the
`Phi` prediction off the unified corner law — a **dossier**, not a construction.
This lane builds the underlying objects:

1. **Newton-polygon reduction** of the corner `(5,20) → (7/5,2)` to
   `ell(C) = x^t c`, `[P,Q] = x^kappa`  (§2).
2. **The C-series leading polynomial `C`**, as the exact solution of the forcing
   ODE that `[P,Q]=x^kappa` imposes  (§3).
3. **The D-transform tower `C^k`** up to the slice where `Phi = f·C^N` lives, and
   the derivation of `N` from that tower's structure  (§4 — the decisive step).
4. **`Phi`** read off the top of the tower  (§5).

Only step 2's chart is a standing assumption (judgment 2); steps 3–4 are exact
given the built `C`.

---

## 2. Corner geometry → reduction parameters

| quantity | value | source |
|---|---|---|
| corner `A_0` | `(5,20)`, `v11=25` | GGV5 line 1679 (F_2 j=1) — **[judgment 1]** |
| final corner | `(7/5,2) = (1,0)+2·(1/5,1)`, type-II.b | GGV5 line 1679 |
| Laurent chart | `X→x^-1, Y→x^5 y`; Jacobian `-x^3` | denominator `l=5` |
| `t` in `ell(C)=x^t c` | `t = l = 5` | `(Y-rX^-5) → x^5(y-r)` |
| `kappa` in `[P,Q]=x^kappa` | `kappa = l-2 = 3` | Laurent Jacobian |
| selected mult `q` | `q = 2` | 2nd coord of final corner |
| `deg C = a_0` | `5` | `C = y^2(y+α)H_2`, `deg H_2 = a_0-q-1 = 2` |
| C-powers `(a,b)` | `(3,5)` = sorted `(m,n)` | `P=C^3`, `Q=C^5+…+F` |

These are exactly the parameters of `PHI_75_125.md`; the point of this lane is
not to re-fix them but to **construct `C` and derive `N`** from them.

**[judgment 2: unreduced polygon]** GGV5/GGHV22 carry out the explicit reduction
only for the `(8,28)` corner; the `(5,20)` reduction is in no paper.  This build
assumes the same standard type-II.b root shift + final Laurent chart, giving
`t=l=5`, `kappa=l-2=3`, `deg C=a_0=5`, `q=2`.  This is the shared conditional
boundary of `CORNER_144_COMPARISON.md` / `PHI_75_125.md` / `PHI_CORNER4.md`, not
new to this lane.  Everything downstream is exact **given** this chart.

---

## 3. The C-series, built by solving the forcing ODE

The operator identity `[P, x^s f/c^b] = x^kappa` (`f = c^b F_s`, `s = v(F) =
kappa+1-at = -11`) is, at `(a,b,t,kappa) = (3,5,5,3)`,

```
15 c f' - 42 c' f = c^3 ,     c = C = y^2 g ,   deg g = a_0-q = 3 .
```

Local orders force `f = A y^5 g^3`; substituting collapses the entire ODE to the
single equation

```
3A (y g' - 3 g) = 1 .
```

Its coefficients force `g_1 = g_2 = 0`, make the `y^3` term resonant (free), and
`g(-1)=0` forces `g_0 = g_3`.  Monic normalisation gives the **forced C-series**

```
g   = y^3 + 1 = (y+1)(y^2 - y + 1) ,     H_2 = y^2 - y + 1  (separable, ≠0 at 0,-1),
C   = y^2 (y^3 + 1) ,                     deg C = 5 = a_0 ,
f   = -(1/9) y^5 (y^3 + 1)^3 ,            deg f = 14 = resonant degree,   A = -1/9 .
```

An independent 15-variable linear solve confirms `f` is the **unique** polynomial
solution of degree `≤ 14`.  The `∞`-leading coefficient `15d - 42·5` is resonant
exactly at `d = 14`.  (Verifier §C.)

---

## 4. THE DECISIVE STEP — the D-transform tower and where N comes from

Write the normalised C-series in the tower parameter `u`:

```
S = Σ_k d_k u^(t-k),     d_t = 1,  d_{t-1} = 0 (x-shift),   d_k := c_k · c^(a(t-k)-1),
```

where `c := C` (the leading polynomial) and `c_k` are the C-series coefficients.
`P = C^a` occupies the linear window `S^a`; the forcing term `F`, hence `Phi`,
occupies the **`C^b` tower `S^b`**, at the slice `(D~^b)_{-j}` = the coefficient
of `u^M` with

```
M = b·t + j ,     j = -s = a·t - (kappa+1) = 11 .
```

**Slice-sum invariant (the crux).**  For any monomial `d_{k_1}···d_{k_b}` in the
`u^M` coefficient of `S^b`, `Σ_i (t-k_i) = M`, so its total `c`-exponent is

```
Σ_i ( a(t-k_i) - 1 ) = a·( Σ_i (t-k_i) ) - b = a·M - b  =:  clear ,
```

which depends **only on the integer slice `u`-power `M`** — never on the
individual `k_i`.  Hence every term of the forcing slice carries the *same*
factor `c^clear`, so

```
Phi = (D~^b)_{-j} = c^clear · F_s = c^(clear-b) · (c^b F_s) = f · c^N ,
      N := clear - b = a·M - 2b .
```

At `(a,b,t) = (3,5,5)`, `M = 5·5 + 11 = 36`, `clear = 3·36 - 5 = 103`, and

```
N = clear - b = 103 - 5 = 98 .
```

The verifier builds `S^5` and checks that **every reachable `u`-slice is
`c`-homogeneous with `c`-exponent `a·M - b`, including the forcing slice `M = 36`
directly** — so `clear = 103` and `N = 98` are read off the built tower, not
assumed.

### Why this upgrades judgment 3 (the non-integral slice index)

`PHI_75_125.md` judgment 3 flagged that the corner-144 clearing derivation used a
**per-term** "forcing slice index"

```
k = t - b·t + s + (b-1)/a ,
```

which is an integer only when `a | (b-1)`.  For `(72,108)` and `(108,144)`,
`(b-1)/a = 1` (integer); for `(75,125)`, `(b-1)/a = 4/3` and `k = -89/3` is
**not** an integer — so "the `C^b` monomial sits at integer tower level `k`,
read `clear = a(t-k)-1` off it" does not transfer verbatim.

The tower build shows this concern **does not reach `N`**: `Phi` is not a single
`d_k`, it is the whole `u^M` slice, and `clear = a·M - b` is a **slice-SUM**
invariant governed by the slice `u`-power `M = b·t + j`, which is **always an
integer**.  The non-integral per-term index cannot move `Phi` to a different
slice or change `clear`.  So `N = 98` is DERIVED, and the general formula

```
N = a·M - 2b = a(bt + j) - 2b = a[ t(a+b) - (kappa+1) ] - 2b
```

is re-derived structurally (not extrapolated), coinciding with the corner-144
formula at all three points.  **PHI_75_125 judgment item 3 → derived.**

---

## 5. Φ emergence

```
Phi = f · C^98
    = -(1/9) y^201 (y^3+1)^101
    = -(1/9) y^201 (y+1)^101 (y^2-y+1)^101 .
```

Exponent bookkeeping from the built objects: `ord_y = rho + qN = 5 + 2·98 = 201`;
`mult_(y+1) = (e+N) = 3 + 98 = 101` (the `(y+1)` factor rides inside `C` with
multiplicity 1 per power); cofactor `= -(1/9) H_2^101`, degree `2·101 = 202`.
Signature `(504, 201, 101, 202)` — matching the landed target exactly.  The
residual `H_2` rides **inside** `C` (as `h_4` did at `(108,144)`), not as a new
unit place (contrast the `(72,108)` quartic).  (Verifier §E.)

---

## 6. Controls — same machinery, known checkpoints

| case | corner | `(a,b,t,kappa)` | `M=bt+j` | `clear` | derived `N` | Phi signature | status |
|---|---|---|---|---|---|---|---|
| `(72,108)` | `(8,28)` | `(2,3,4,2)` | `17` | `31` | `28` | `(238,204,30,4)` | STATE.md audited |
| `(108,144)` | `(8,28)` | `(3,4,4,2)` | `25` | `71` | `67` | `(550,205,69,276)` | corner-144 |
| **`(75,125)`** | **`(5,20)`** | **`(3,5,5,3)`** | **`36`** | **`103`** | **`98`** | **`(504,201,101,202)`** | **derived here** |

For `(72,108)` (the `r=0` resonance-gap case) the ODE solution `f` carries the
audited quartic unit cofactor, but the tower still yields `N=28` — the tower
derivation of `N` is orthogonal to the resonance-gap distinction.  For
`(108,144)` the full `C, f, N, Phi` chain is rebuilt from scratch and reproduces
`(550,205,69,276)`.  (Verifier §F.)

---

## 7. `[judgment]` list — where this is conditional

1. **[judgment 1: chain data]** Corner `(5,20)→(7/5,2)`, `(m,n)=(3,5)`, `k=1`
   from GGV5 line 1679, re-checked by the Diophantine `(m+n)qk - n(ql-p) = k`.
   Primary-source.
2. **[judgment 2: unreduced polygon]** The `(5,20)` reduction is in no paper; the
   standard type-II.b root shift + Laurent chart is assumed (`t=l=5`,
   `kappa=l-2=3`, `deg C=a_0=5`, `q=2`).  Shared conditional boundary with
   CORNER_144 / PHI_75_125 / PHI_CORNER4.  Unchanged by this lane.
3. **[judgment 3 → DERIVED]** `N = 98` is now derived from the built D-transform
   tower: `clear = a·M - b` is a slice-sum invariant at the integer forcing slice
   `M = b·t + j`, so the non-integral per-term slice index `(b-1)/a = 4/3` does
   not affect it.  The N-formula is thereby re-derived, not extrapolated.
   *(This is the item this lane discharges.)*

**Scope note.** This build derives `N` and reads off `Phi`.  Later phases of the
transfer test (the full D-transform G-system, the pre-resultant / window layer)
remain the pending items flagged in `case_compiler.py`'s transfer inventory;
they are out of scope here.  The C-series, its tower, and `N` — the pending
item 3 of `PHI_75_125.md` — now exist as a construction.

---

## Files

- `C_SERIES_75_125.md` — this writeup.
- `c_series_75_125.py` — the construction (C-series ODE solve, D-transform tower,
  N derivation, Phi, controls).  Exact sympy; run end-to-end.
- `c_series_75_125_verify.py` — independent PASS/FAIL checker (`--quiet`, exit 0);
  36 checks incl. the tower slice-homogeneity lemma at the forcing slice and both
  control checkpoints.
