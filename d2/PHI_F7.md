# The gap>0, r>0 regime: old conjecture REFUTED, ramified regime law found

## Verdict

The last unprobed regime of the corner law -- `gap > 0` with `r > 0` (F3, F7,
F10, F15, F16), where PHI_F14.md §5 flagged the unified cofactor formula
`gap + r(e+N)` as a conjecture -- is now tested at **four corners**.  The
verdict is **DIFFERS**, in a precisely structured way:

* the **deg and ord_y components hold on every branch** --
  `deg = res + N·a0`, `ord = rho + N·q` (res = pure + gap);
* the **mult and cofactor components fail, provably** -- and are replaced in
  this regime by the amended **ramified law** (four fresh exact points):

```
mult  = dg·(e+N) − (dg−1)          cofactor = gap + r
```

```
F7  (42,147) : f = (1/10)  y^21 (y+1)^11 (9y^2+3y−1)
               Phi sig (250, 165, 83, 2)      [old law said (250,165,42,43)]
F3  (75,50)  : f = (1/42)  y^4  (y+1)^3  (25y^2+15y−3)
               Phi sig (189, 112, 75, 2)      [old law said (189,112,38,39)]
F10 (196,112): f = (1/3740) y^10 (y+1)^13 (2401y^4+5831y^3+4165y^2+595y−85)
               Phi sig (1917, 820, 1093, 4)   [old law said (1917,820,274,823)]
F16 (99,165) : f = (1/330) y^15 (y+1)^5  (243y^4+81y^3−27y^2+15y−10)
               Phi sig (528, 407, 117, 4)     [old law said (528,407,59,62)]
```

F16 (`gap=3`) is the separating experiment for the unit degree: `gap+r = 4`
beats `dg = 2`.  And `cofactor = gap + r` **retro-explains the audited
(72,108) quartic**: its degree 4 = gap + r = 4 + 0 -- one formula now covers
the unit cofactor at all twelve known points (`gap=0` corners have `u = 1`
trivially since `gap + r·0-contribution` -- see unification note below).

Checked by `phi_f7_verify.py` -- **62/62 pass** -- with independent routes
throughout (fresh corner arithmetic; every claimed `f` checked against the
ODE by direct differentiation; uniqueness by full generic linear solve, not
the recurrence; the F3 signature cross-checked by full expansion of Phi
itself; the impossibility theorem re-proven by quadratic-discriminant
arguments; off-branch refutation spot checks).  Derivation: `phi_f7.py`.

## Why the old conjecture HAD to fail: the dg=2 obstruction theorem

With fully generic monic residual `g = y^2 + g1·y + g0`, the triangular
solve of the ODE `a·t·c·f' − a·coef·c'·f = c^e` (c = y^q·g) leaves one
resonant free coefficient (`f_res`) and exactly two consistency conditions,
both linear in it.  The eliminant factors **completely**:

```
F7 :  E ~ g0^27 · (3g0 − 2g1^2)   · (4g0 − g1^2)^6
F3 :  E ~ g0^6  · (5g0 − 3g1^2)   · (4g0 − g1^2)^2
F16:  E ~ g0^18 · (4g0 − g1^2)^3  · (54g0^2 − 126g0·g1^2 + 35g1^4)
```

(weighted-homogeneous in (g1, g0) with weights (1,2) -- root *ratios* are
forced; the overall scale is the same gauge freedom as at dg=1).  So a
polynomial solution exists only when the residual is:

1. **RAMIFIED** -- the discriminant factor `4g0 = g1^2`; the root-shift gauge
   `g(−1)=0` then forces `g = (y+1)^2` (the double root sits AT −1); or
2. a **COMPLEX-CONJUGATE PAIR** -- every ratio-factor root has `w = g0/g1^2 >
   1/4`, so `g(−1)=0` (i.e. `1 − g1 + w·g1^2 = 0`) has **no real solution**.

A simple real root at −1 -- what the old `mult = e+N` needs -- is therefore
**impossible** in this regime.  (Independently: `mult = e+N > res ≥ deg f`
rules it out by degree count alone.)  The `gap=0` corners never face this:
there the ODE forces `g = y^dg + 1` with `dg` odd -- simple root at −1.  A
checked survey mini-lemma sharpens the split: **dg is even exactly on the
`gap>0, r>0` rows** of all 15 standard-chart families.

On the complex-pair branch the second place disappears entirely:
`Phi sig = (deg, ord, 0, deg−ord)` (recorded for F7 and F3 with rational
scale representatives).  Deg and ord are branch-independent.

## Unification note

Over the twelve exact points now known, the corner-law components behave as:

| component | formula | status |
|---|---|---|
| `deg` | `res + N·a0 = pure + gap + N·a0` | **all twelve points, every branch** |
| `ord_y` | `rho + N·q` | **all twelve points, every branch** |
| `mult` | `mult_g·(e+N) − (mult_g−1)`, `mult_g = mult_{y+1}(g)` | unifies: simple (`mult_g=1`) gives `e+N`; ramified (`mult_g=dg`) gives the new formula |
| `cofactor` | `(dg − mult_g)·(e+N) + gap + (mult_g − 1)` | reduces to `gap + r(e+N)` when `mult_g=1` (since `r = dg−1`) and to `gap + r` when `mult_g=dg` |

so the *shape* `Phi = A · y^(rho+Nq) · (y+1-part) · units` persists; what the
regime changes is **how much of `g` ramifies onto the second place**.  The
`mult_g` unification rows are observations consistent with all twelve points,
not independently derived laws.

## `[judgment]` list -- where this is conditional

1. **`[judgment: chain data]`** Corner rows from the GGV5 `v11 ≤ 35` tables,
   same transcription as `phi_corner4.py`; Diophantine identity re-checked
   independently for every family used.  Primary-source.
2. **`[judgment: unreduced polygon]`** As with all non-(8,28) corners: the
   `(6,15)`, `(5,20)`, `(7,21)`, `(9,24)` reductions are performed in no
   paper; the standard type-II.b root shift + Laurent chart is assumed
   (`t = l`, `kappa = l−2`, `deg C = a0`, `q` from the table).  Same
   conditional boundary as CORNER_144 / PHI_75_125 / PHI_CORNER4 / PHI_F14.
3. **`[judgment: branch selection]`** The ODE model itself admits BOTH real
   branches.  The ramified branch `g = (y+1)^dg` is selected here by
   continuity with the audited pattern (`(y+1) | C` at every one of the
   seven previous points and in the audited (72,108) premise `C4 = y^7(y+1)`),
   and the root-shift gauge picks −1 as the double root's position.  But the
   actual tower C-series for these families is built in no paper; the
   complex-pair branch cannot be excluded without it.  Both branches are
   derived and recorded; only the ramified one is called "the" signature.
4. **`[judgment: N formula]`** `(b−1)/a` is integral for F7 (3) and F3 (1) --
   the less-conditional class -- but NOT for F10 (3/2) or F16 (4/3), which
   sit in the same more-conditional class as `(75,125)`.  All four matches of
   deg/ord are evidence for the N-formula, not independent of it.
5. **`[judgment: dg=4 branch completeness]`** The obstruction is fully
   factored (hence branch-complete) only at `dg = 2` (F7, F3, F16).  At
   `dg = 4` (F10) the ramified-branch solution is EXHIBITED (exact, unique
   for that `g`) but the full branch variety of the 3-parameter residual
   space was not enumerated.  F15 (`dg = 4`) was not attempted.
6. **`[judgment: mult_g unification]`** The last table above interpolates
   between the regimes with `mult_g`; at present `mult_g` takes only the
   values {1, dg} at derived points, so intermediate ramification
   (`1 < mult_g < dg`, possible at `dg ≥ 3` only) is unobserved and the
   interpolation is untested there.
