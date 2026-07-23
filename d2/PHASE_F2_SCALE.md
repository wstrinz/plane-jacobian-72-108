# Phase F, work item F2 -- SCALING the divisor-reconstruction kill test

**Date:** 2026-07-23. **Status: PENDING AUDIT** (same-author layer over
`phase_f2_pilot.py`, `alt_residue_congruences.json`, `convolution_descent.py`;
no independent audit). Every kill below is a *candidate* kill.
**New files (uncommitted):** `phase_f2_scale.py` (engine + alt driver),
`phase_f2_scale.json` (census), this doc, `phase_f2_scale_verify.py`
(independent verifier, **ALL CHECKS PASS**). READ-ONLY on every audited artifact.

## 0. What this scales

`PHASE_F2_PILOT.md` killed three hand-worked alternate-regime configurations by
reconstructing the finite-place divisors of the forced defect-0 polynomials and
feeding the *determined* coefficients into the level-0 `h_0` tie tower. This item
runs that test **mechanically over the whole forced-defect-0 frontier**, plus the
Galois-stable `a11_b1111` flagship, the sub2 tie-tower analogue, and the Pilot C
exact certification.

## 1. Alt front -- the reconstructible frontier

`phase_f_defects.json` marks 90 alt survivors with `delta_d1 = delta_sigma = 0`.
Adding the third relevant obligation `e` forced defect-0 (`deg e = a + sum b`,
so e's q-cofactor has no free part) leaves **40 states** in which d1, sigma AND
e are *all* forced defect-0. For every one of the 40 the ALT_REGIME_L2 sec.2 cone
admits **exactly one** valuation split (verified: `nsplit == 1` for all 40) -- so
the reconstruction `p = lambda_p (y+1)^{v_t} prod_i (y-r_i)^{v_r_i}` is
well-defined, not a choice. `d2` is left a **free** polynomial of the state's
`deg d2` (conservative: extra d2 freedom only enlarges the solution set, so a
kill under free d2 is sound -- exactly the Pilot B model).

Field of definition is minimised per state: a Galois-stable full-`q` or pure
`t`-place divisor computes **over Q**; a 3-equal-multiplicity q-factor is the
complement `q/(2048(y-r))` of a **single** marked root; a lone active root is
`Q[r]/(q)`; a genuinely-two-distinct-root divisor uses `Q[r1,r2]/(q,q)` with
`r1 != r2` saturated.

### Census (40 states, verdicts PENDING AUDIT)

| verdict | count |
|---|---:|
| **KILLED** (saturated ideal = unit after saturating the leading scalars) | **20** |
| NARROWED / UNOBSTRUCTED (tie tower not the unit ideal at the tracked depth) | 18 |
| PENDING (two-root + free deg-6 d2 GB blowup / depth-cap timeout) | 2 |

### The four killing families

| family | field | states | KILLED | kill depth | note |
|---|---|---:|---:|---:|---|
| `a11_b1111_T1` (**flagship**) | Q | 8 | 7 | 2 (`d2=5`: 8) | Galois-stable: d1 ~ (y+1)^5 q, e ~ (y+1)^11 q, sigma = S(y+1)^3 |
| `a11_b3100_T2` | Q[r1,r2], r1!=r2 | 8 | 6 | 2 | two distinct marked roots |
| `a12_b1110_T2` | Q[r]/(q) | 8 | 6 | 2 | single marked root (complement of 3-equal) |
| `a14_b0000_T2` | Q | 1 | 1 | 9 | whole-state kill = landed Pilot B |

Support `15` = `-729 E (9 E^3 + 8 X^5)` (flagship) and support `12` (the two T2
families) are the level-0 initial forms carrying the kill; each reproduces its
`alt_residue_congruences.json` catalog form exactly.

### The d2-freedom boundary (the honest limit)

Within **every** killing family the states with `deg d2 in {none,0,1,2,3,4}` die
at depth 2, but `deg d2 in {5,6}` do **not** die at any accessible depth
(NARROWED / PENDING). A free d2 cofactor of degree >= 5 rescues the tower -- the
same phenomenon as the `PHASE_F2_PILOT` control (retained freedom => UNOBSTRUCTED).
The kills are therefore driven, exactly as the pilot argued, by the number of
*determined* sub-leading coefficients exceeding the residual freedom against the
tie depth; a large free d2 restores that freedom.

### Whole cell / branch closure

**No whole cell or branch is fully closed.** But two whole **branches**,
`a11_b3100_T2` and `a12_b1110_T2`, are *entirely* defect-0 (all 8 states each)
and **6 of 8 are killed** in each -- the residual 2 are the `deg d2 in {5,6}`
survivors above. So these two branches are 3/4 closed by reconstruction alone,
the closest the layer comes to a branch kill; closing them needs the deeper
d2-vs-tie-depth argument or d2's own divisor (not tracked for alt).

## 2. Sub2 fully-forced 243 -- tie-tower analogue

The 243 sub2 states with every `core+gT` defect 0 were located. Demonstrator
(`convolution_descent.py`, the requested master-identity machinery):
`sub2:a10_b0000_T1`, degrees `(deg d2,d1,sigma,e) = (0,0,0,10)`, reconstructs to
`d2,d1,sigma` **constants** and `e = g0 (y+1)^10` -- fully determined up to four
leading scalars. The master coefficient at the top degree (241) is a **single
nonzero relation** among the four scalars -- a satisfiable hypersurface -- so the
saturated GB is not the unit ideal at depth 1 (**NARROWED**), the exact sub2
analogue of the alt free-scalar CONSTRAINT. **Coverage boundary:** the 40 b0000
states coincide with the existing generic `batch_convolution_sub2` run (no
q-support to restore); the 203 `b != 0000` states need reconstruction of the
q-support that batch *dropped*, which requires the sub2 (geometric q-coprime)
per-state finite-place split -- not recorded in `phase_f_defects.json`, DEFERRED.

## 3. Pilot C -- exact certification attempt

State `a11_b3000_T1 (0,9,12,14)` with `sigma = (y+1)^3 C` (C generic degree 9),
`d1,e` forced defect-0, `d2 = D`. With `X=1`, `j0 = 2187 c9^2(4 c9^3 + 1)` (the
`4 c9^3 + 1 = 0` torus branch). **Triangular elimination** solves `c8..c0` from
`j1..j9` (all pivots nonzero) and isolates **4 residuals `j10..j13` in `(D,E)`
over `Q(r,c9)`** -- confirming the `PHASE_F2_PILOT` sec.5 structure *exactly*. The
coherence `sigma = 4 d0 - d2^2` is already baked into `h_0` (`MONOMIALS[0]` is
written in `sigma` via that substitution), so at the infinity tie it adds no new
relation beyond fixing `d0`. **The exact unit-ideal certification STALLS**: the
saturated grevlex GB over `(c0..c9,D,E,r,w)` exceeds 6 min (through depth ~6), and
the resultant elimination over the degree-12 number field
`Q[r,c9]/(q, 4c9^3+1)` exceeds 5 min on `Res_E` of a deg-15 x deg-20 pair.
**Verdict: PENDING** -- residual structure confirmed, exact certification
unresolved, matching the pilot's "stalled GB" note.

## 4. Verification (`phase_f2_scale_verify.py`, ALL CHECKS PASS)

Independent of `phase_f2_scale.py`. Re-derives **two new** flagship kills by a
simpler chain (Galois-stable `prod_i (y-r_i) = q/2048`, everything over Q, no
`q(r)` reduction): **K1** `a11_b1111_T1`, `d2 = 0` -- initial form
`-729 E (9 E^3 + 8 X^5)` (= catalog support 15), `(j0,j1,sat X S E)` = unit ideal
over Q (depth 2); **K2** same family with a free constant `d2 = D` -- still
`(j0,j1,sat D X S E)` = unit ideal (d2-freedom at that level does not rescue it).

## 5. Honest / ambiguous points -- [judgment]

- **[J1] The 20 kills are exact saturated-Groebner facts** over the stated field
  (Q, Q[r]/(q), or Q[r1,r2]/(q,q) with r1!=r2), with **d2 left fully free** at its
  state degree -- so they are sound whatever d2's true divisor. The flagship
  `a11_b1111_T1` (7 kills over Q) is independently reproduced by the verifier.
- **[J2] The kill is a d2-freedom threshold, not universal on a support.** Every
  killing support (15, 12) kills `deg d2 <= 4` at depth 2 but leaves `deg d2 in
  {5,6}` NARROWED/PENDING. This is not a bug: it is the pilot's own mechanism (a
  large free cofactor restores the freedom the reconstruction removed). Reporting
  those as kills would be unsound; I do not.
- **[J3] No whole branch closes, but two are 6/8 closed.** `a11_b3100_T2` and
  `a12_b1110_T2` are entirely defect-0; reconstruction kills 6 of 8 states in
  each. The survivors are exactly the large-free-d2 states of [J2].
- **[J4] Two states are PENDING for cost, not mathematics.** `a11_b3100_T2`
  `deg d2=6` (two marked roots + free deg-6 d2 + distinctness saturation) and one
  `a12_b1110_T2` `deg d2=6` blew up the GB; guarded out, not resolved.
- **[J5] Pilot C is PENDING, sharpened not solved.** The triangular reduction is
  a genuine advance (it pins the 4-in-2 residual system the pilot described), but
  the exact number-field certification still stalls -- I mark it PENDING.
- **[J6] The sub2 front is a demonstration + boundary, not a census.** Only the
  mechanism is shown (one state) and the b0000 coincidence noted; the 203
  q-support-bearing states need a sub2 split reconstruction that is deferred. I do
  not claim a sub2 census.
