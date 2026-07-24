# Sixth + seventh points: F14 confirms t=3, and the resonance-gap story generalizes

## Verdict

Two new out-of-sample derivations by the corner-144 template, both **matching the
corner law with zero fitting freedom**:

```
(66,231)  F14 j=0, corner (9,24) -> (7\3,4), (a,b)=(2,7), t=3:
    Phi = -(1/10) y^165 (y^5+1)^42
    signature (deg, ord_y, mult_(y+1), cofactor) = (375, 165, 42, 168)   MATCHES

(48,64)   F1 j=0, corner (4,12) -> (7\4,3),  (a,b)=(3,4), t=4, r=0, gap=1:
    Phi = (1/15) y^205 (y+1)^69 (4y-1)
    signature = (275, 205, 69, 1)                                        MATCHES
```

F14 is the first point at **t=3**, the first with **e=6** (previous coverage e in
{2,3}), and the first with q=4.  F1 is the first **fresh derivation in the
resonance-gap regime** (r=0, gap>0) — the class that until now contained only the
audited (72,108) exception.

The law now has **seven exact points** and one unified statement:

```
signature = ( (e*a0 - q + 1) + gap + N*a0,   (e-1)q + 1 + N*q,
              e + N,                          gap + r*(e+N) )
N   = a[t(a+b-1) + 1] - 2b            (kappa = t-2 eliminated, PHI_CORNER4)
gap = (q-1) - a0/t                    (mini-lemma, see below)
```

| case | corner | `(a,b)` | `t` | `(a0,q,e,r,gap,N)` | signature |
|---|---|---|---|---|---|
| `(72,108)` | `(8,28)` | `(2,3)` | `4` | `(8,7,2,0,4,28)` | `(238,204,30,4)` audited |
| `(108,144)` | `(8,28)` | `(3,4)` | `4` | `(8,3,2,4,0,67)` | `(550,205,69,276)` |
| `(75,125)` | `(5,20)` | `(3,5)` | `5` | `(5,2,3,2,0,98)` | `(504,201,101,202)` |
| `(56,84)` | `(7,21)` | `(2,3)` | `7` | `(7,2,2,4,0,52)` | `(377,107,54,216)` |
| `(50,75)` | `(5,20)` | `(2,3)` | `5` | `(5,2,2,2,0,36)` | `(189,75,38,76)` |
| **`(66,231)`** | **`(9,24)`** | **`(2,7)`** | **`3`** | **`(9,4,6,4,0,36)`** | **`(375,165,42,168)`** |
| **`(48,64)`** | **`(4,12)`** | **`(3,4)`** | **`4`** | **`(4,3,2,0,1,67)`** | **`(275,205,69,1)`** |

Coverage census: `t in {3,4,5,7}`, `e in {2,3,6}`, `a0 in {4,5,7,8,9}`,
`q in {2,3,4,7}`.  Checked by `phi_f14_verify.py` — **37/37 pass** — with
independent routines (trial-division signature extraction, full linear solves for
uniqueness, all four previous corners re-derived as controls).  Derivation:
`phi_f14.py`.

---

## F14: one-line collapse at t=3

ODE `6 c f' − 34 c' f = c^6`, `c = y^4 g`, `deg g = 5`.  The bracket constants are
`(t·rho − coef·q, te − coef) = (−5, 1)`, so the coefficient multipliers
`−5 + i` kill `g_1..g_4`, the top is resonant, and `g(−1)=0` + monic force
`g = y^5+1`, `A = −1/10`, `f = −(1/10) y^21 (y^5+1)^6` (deg 51, the unique
polynomial solution of degree ≤ 51 — full 52-unknown linear solve).  Then
`N = 36` and `Phi = −(1/10) y^165 (y^5+1)^42`.

The residual `H = y^4−y^3+y^2−y+1` (10th cyclotomic) is the same as
`(108,144)` and F9, because all three corners have `dg = a0−q = 5`: **the
residual is indexed by `dg` alone** (`g = y^dg + 1` throughout).

## F1: the resonance-gap regime generalizes — unit cofactor of degree exactly `gap`

For F1, `r = 0` and `gap = 1`: the pure ansatz `A y^4 (y+1)^2` (degree 6) fails at
the top — resonance is broken exactly as in `(72,108)`.  Solving the ODE
`12 c f' − 21 c' f = c^2` with **fully generic** `f` (degree allowed 2 past
resonant) gives a *unique* polynomial solution, and it factors as

```
f = (1/15) y^4 (y+1)^2 (4y−1)   —  pure-ansatz shape  ×  unit cofactor u,
u = (4y−1)/15:  deg u = gap = 1,  u(0) ≠ 0,  u(−1) ≠ 0.
```

This is the degree-1 analogue of the `(72,108)` quartic
`2048y^4−512y^3+320y^2−240y+195` (whose degree, 4, is exactly that case's gap).
The amended `r=0` law — `deg` offset by `gap`, `cofactor = gap` — now has **two
exact points, one of them derived from scratch**, and with `gap = 0` it reduces
to the ordinary law, giving the single unified statement above.

Amusing check: F1 and `(108,144)` share `(a,b,t,q) = (3,4,4,3)` and hence
`N = 67`, `ord = 205`, `mult = 69` — they differ only through `a0` (4 vs 8),
and the signatures differ exactly as the law says they must.

## Mini-lemmas (checked over all 15 standard-chart families)

```
gap = (q − 1) − a0/t          so   gap = 0  <=>  a0 = t(q−1)
dg  = a0 − q                  (residual index; g = y^dg + 1)
te − coef = t − kappa − 1 = 1 (the kappa = t−2 identity in ODE form)
```

The first gives a *corner-data criterion* for when the resonance-gap regime
occurs — no ODE needed.  The `r=0` rows of the survey are precisely F1, F5, F17
(+ the GGHV `(72,108)` corner).

## `[judgment]` list — where this is conditional

1. **`[judgment: chain data]`** Corner rows from the GGV5 `v11 <= 35` tables,
   same transcription as `phi_corner4.py`; Diophantine identity re-checked here
   for every family used.  Primary-source.
2. **`[judgment: unreduced polygon]`** As with all non-(8,28) corners: the
   `(9,24)` and `(4,12)` reductions are performed in no paper; the standard
   type-II.b root shift + Laurent chart is assumed (`t = l`, `kappa = l−2`,
   `deg C = a0`, `q` from the table).  Same conditional boundary as
   CORNER_144 / PHI_75_125 / PHI_CORNER4.
3. **`[judgment: c normalization]`** For F1 (`dg = 1`) the derivation takes
   `c = y^q (y+1)` — the analogue of the audited `(72,108)` premise
   `C4 = y^7(y+1)` and of the forced `g(−1) = 0` at every `dg > 1` corner, but
   at `dg = 1` the ODE itself does not force the root's position; the root
   shift places it at `−1` by construction.  Flagged, not proved here.
4. **`[judgment: N formula]`** Extrapolated to `t = 3`.  `(b−1)/a` is integral
   for both new points (3 for F14, 1 for F1), so the corner-144 slice-index
   picture transfers verbatim — both are in the *less*-conditional class (like
   F9, unlike `(75,125)`).  The matches are evidence for the formula, not
   independent of it; a from-scratch C-series build at one new corner remains
   the way to discharge this item entirely.
5. **`[judgment: untested regime]`** `gap > 0` with `r > 0` (F3, F7, F10, F15,
   F16) has **no derived point**; there the unified cofactor formula
   `gap + r(e+N)` is a conjecture.  This is now the *only* regime of the
   length-1 survey with no out-of-sample test, and the named next experiment
   (F7 at `(42,147)` is the smallest).
