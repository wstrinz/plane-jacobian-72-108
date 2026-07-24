# Galois descent across the full residue library -- C20 proved, C01..C23 classified

**Date:** 2026-07-23. **Status:** the C08 pilot (`GALOIS_DESCENT_PILOT.md`)
extended to the entire 23-shape library. Checker:
`galois_library_verify.py` (69 checks; re-derives every equation from
`f31_graded.txt` through the audited `residue_lemmas_verify.py` machinery --
no coefficient is hand-copied).

## 1. Headline

1. **C20 is proved through the full S4/C4/A4/D4-V4 branching**, exactly as
   C08 was: residue discriminant class `170`; for the fixed forcing quartic
   `q` (S4, disc class `17`) the obstruction is
   `Obs = squarefree(170*17) = 10 != 1`, so **C20 KILLS**, with sharpness
   and branching-necessity witnesses (sec. 3).
2. **The kill/constraint split of the whole library is now predicted by
   support geometry, not case-by-case arithmetic.** A shape can carry an
   arithmetic obstruction on the coefficient torus iff its exponent support
   is *collinear* with ratio-degree 2 -- and exactly two shapes qualify:
   C08 and C20, the two audited kills. Every other shape is solvable over
   **every** field, so it is a CONSTRAINT for every conceivable forcing
   quartic, not just ours.
3. **No shape in the library needs a higher resolvent.** The C12-type worry
   flagged in the pilot's sec. 4 dissolves: on the free torus, C12
   (`73X^4 + 4S^3`) is *linear* in the primitive character `S^3/X^4`
   (gcd(3,4)=1 makes the character surjective over every field), hence
   solvable everywhere. Higher resolvents would matter only on constrained
   sub-loci where deeper cascade obligations fix some variables -- scoped,
   not needed for the library itself.

## 2. The classification machine

Write a shape as `R(D,X,S,E) = sum_i c_i M_i` with all source coefficients
`c_i != 0`. Let `V = {exponent(M_i) - exponent(M_0)}` in `Z^4`.

- **Collinear (`rank V = 1`).** Let `g` be the primitive direction and
  `t_i` the integer multiples, shifted so `min t_i = 0`. Then on the torus
  `R = M_0 * chi_g^shift * P(chi_g)` with `P(rho) = sum c_i rho^{t_i}` --
  the checker verifies this reduction *exactly*. Since `g` is primitive
  there is an integer vector `u` with `<g,u> = 1`; the substitution
  `x_j = tau^{u_j}` maps `chi_g` to `tau` (checked exactly), so the
  character is **surjective on F\*-points over every field F** and torus
  solvability of `R` is precisely "`P` has a nonzero root in `F`".
  - `deg P = 1`: root `-c_0/c_1 in Q*`; the checker *constructs* the torus
    point `x_j = lambda^{u_j}` and substitutes it back into the source
    equation. Solvable over every field.
  - `deg P = 2`: the discriminant of `P` is a rational constant; its square
    class `Delta` is the shape's **residue obstruction class**. Both roots
    are nonzero (constant term `!= 0`). Solvable over `F` iff
    `sqrt(Delta) in F`.
  - `deg P >= 3`: would need the Galois theory of a higher-degree
    resolvent. **Zero shapes of this kind occur.**
- **Non-collinear (`rank V >= 2`).** The audited rational witness from
  `RESIDUE_LEMMAS.md` is re-substituted exactly: a rational torus point
  exists, so the shape is solvable over every field containing `Q`.

### Census (all 23 shapes)

| class | shapes | count | verdict |
|:--|:--|--:|:--|
| QUADRATIC-OBSTRUCTION | C08 (`Delta=105`), C20 (`Delta=170`) | 2 | **KILL over fields not containing `sqrt(Delta)`**; both kill for the fixed `q` |
| LINEAR in a character | C01, C02, C03, C05, C06, C11, C12, C13, C14, C15, C16, C17 | 12 | solvable over EVERY field; CONSTRAINT always |
| MULTIDIM (rank >= 2) | C04, C07, C09, C10, C18, C19, C21, C22, C23 | 9 | rational witness; CONSTRAINT always |
| HIGHER (deg >= 3 collinear) | -- | 0 | none occur |

The obstruction subgroup of `Q*/(Q*)^2` (pilot sec. 4's "reusable machine")
is therefore: trivial for 21 shapes, `<105>` for C08, `<170>` for C20. The
**joint** subgroup `<105,170> = {1, 105, 170, 714}` meets the fixed `q`'s
split-field square-class subgroup `{1, 17}` only in `1` -- C08 and C20 kill
*simultaneously* in the same splitting field.

## 3. The C20 proof (mirror of the C08 pilot)

C20 (source-derived, primitive form `61X^2D^2 + 6XDE - 11E^2 = 0`, all of
`X, D, E` nonzero) reduces via `r = E/(XD)` to `11r^2 - 6r - 61 = 0`,
discriminant `2720 = 4^2 * 170`, square class **170**; both roots
`r = (3 +- 2 sqrt(170))/11` are nonzero (product `-61/11`). So C20 is
solvable over `F` iff `sqrt(170) in F`. Branching on the Galois type of the
forcing quartic:

- **S4 / C4** (unique quadratic subfield `Q(sqrt(disc))`): kill iff
  `Obs = squarefree(170 * disc-class) != 1`. Fixed `q`: disc class `17`,
  `170*17 = 2890 = 17^2 * 10`, class **10** -> **KILL**. (Note `170 = 2*5*17`
  shares the factor 17 with the disc class; the obstruction survives as the
  cofactor 10. A forcing quartic with disc class exactly `170` is the unique
  S4-escape.)
- **Sharpness witness.** `qs = y^4 + 7y^2 - 8y + 10` is irreducible, S4,
  with `disc(qs) = 459680 = 170 * 52^2` *exactly*. Then `Obs(qs) = 170*170`
  is a square, the obstruction vanishes, and indeed
  `sqrt(170) = Vandermonde(qs)/52 in L`, so `E = (3+2 sqrt(170))/11`,
  `X = D = 1` solves the source C20 over `split(qs)` (verified exactly).
  `Obs` is sharp in the cyclic regime, exactly as for C08.
- **A4**: no quadratic subfield, `170` is not a square -> C20 **always
  kills**. Live test: `qa = y^4 + 8y + 12` (irreducible, order 12, square
  disc) -> kill.
- **D4 / V4** (three quadratic subfields): `Obs` alone is insufficient.
  Witness: `qv = y^4 - 344y^2 + 28224`, the minimal polynomial of
  `sqrt(170) + sqrt(2)`, is irreducible V4 with *square* discriminant; the
  naive cyclic-regime obstruction would be `squarefree(170*1) = 170`,
  predicting a kill -- but the exact identity
  `(512x - x^3)^2 = 170 * 336^2  (mod qv)`
  puts `sqrt(170)` in `L = Q(x)/(qv)`, so C20 is solvable there. The Galois
  branching is necessary, not decorative.

Master criterion in all cases: **C20 kills over `L` iff `sqrt(170) not in
L`** -- same theorem shape as C08, second data point for the
configuration-space descent of `INDUCTIVE_PROGRAM.md` NEW IDEA 2.

## 4. What this changes

- `RESIDUE_LEMMAS.md`'s classifications (2 KILLs, 21 CONSTRAINTs) were
  proved there by exhibiting witnesses case by case. This file upgrades
  that to a **structural theorem**: the split is forced by support
  geometry, and the two kills are the *only possible* kills on the free
  torus no matter what forcing quartic any other case of the family
  presents. For the inductive program this is the desired shape: the
  obstruction data (`Delta`, support direction `g`) are *corner-independent
  constants of the residue relation*, while the forcing side enters only
  through `disc(q)` and the Galois type -- cleanly separated inputs for a
  case compiler.
- For any future family member: compute its forcing quartic's Galois label
  and disc class; the C08/C20 kills transfer iff the label is S4/C4/A4 and
  the disc class avoids `{105, 170}` (S4/C4) -- a two-line check.

## 5. Judgment list  [judgment]

- **[judgment] G1 -- torus scope.** All classifications are on the free
  coefficient torus (all displayed variables nonzero), which is exactly the
  arena of the residue lemmas' leading coefficients. Deeper cascade
  obligations that *constrain* variables to sub-loci could resurrect
  obstructions for the 21 solvable shapes; nothing here claims otherwise.
- **[judgment] G2 -- D4 stays witness-decided.** As in pilot J4, the D4
  branch is decided only by explicit membership witnesses; no closed
  coefficient polynomial is claimed, and no D4 instance is used anywhere.
- **[judgment] G3 -- the witness quartics are synthetic test data.** `qs`,
  `qv`, `qa` were found by search/construction; every property used
  (irreducibility, Galois label, disc identity, membership identity) is
  verified from scratch in the checker, so nothing depends on how they were
  found.
- **[judgment] G4 -- surjectivity is constructive.** The checker does not
  cite the character-surjectivity fact abstractly; it exhibits the Bezout
  vector `u` and checks `chi_g(tau^u) = tau` exactly, then builds each
  linear shape's torus point from it.

Run `python galois_library_verify.py` (add `--quiet` for the summary line)
to re-derive the 23 equations from source and re-check all 69 claims.
