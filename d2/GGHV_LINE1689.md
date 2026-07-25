# GGHV_LINE1689 — independent reproduction of GGHV22's unverifiable CAS step

> **STATUS (2026-07-24): REPRODUCED — exactly, with no extra hypothesis, no
> scalar factor, and no hidden saturation.** Verified three independent ways
> (Singular elimination, ideal-membership, saturation cross-check) plus a
> CAS-free certificate re-verified in sympy. This is a result *about someone
> else's paper*, offered as an independent verification.

## 1. What was unverifiable

GGHV22 (`arXiv:2204.14178`, Guccione–Guccione–Horruitiner–Valqui) is the source
of the premise our whole campaign rests on: any JC(2) counterexample has
`max(deg P, deg Q) >= 125` except possibly `(72,108)`. At tex line **1689** the
paper writes, verbatim:

> "We consider this as a system of 9 equations and using a CAS (for example
> Mathematica) we eliminate the variables `d_{-10}, d_{-8}, d_{-7}, d_{-6},
> d_{-5}, d_{-4}, d_{-3}, d_{-2}`, obtaining"

followed by `\eqref{ecuacion principal}`:

```
18 C_3^23 d_1 (d_-1)^6 F_-4  +  8 C_3^69 F_-4^3  +  27 d_0 (d_-1)^9  =  0
```

The reader is asked to trust an unshown Mathematica run. The step is
**load-bearing** for GGHV's exclusion of `(66,99)` and for the `(9,27)` instance
of `(72,108)`, with an analogous step at line 2056 for `(7,21)`.

Context for why this matters: GGHV is an unpublished preprint, **v1 only** since
2022-04-29 — no v2, no journal reference, no DOI, no erratum — with a Semantic
Scholar citation count of **1** (the sole citer being Nguyen). Theorem 2.1 has
had essentially no external scrutiny.

## 2. The observation that collapses the problem

`C_3` and `F_-4` occur in the nine equations **only** through the single product

```
E := F_-4 * C_3^23
```

and the target's `8 C_3^69 F_-4^3` is exactly `8 E^3`. So the elimination is not
a large computation at all: it lives in `Q[d_1, d_0, d_-1, E]` with no `y`, no
`lambda`, and no `C_3` arithmetic. The paper's framing overstates the difficulty.

After the triangular part (the `(D~^2)_{-1..-5}` and `(D~^2)_{-7}` equations
solve for `d_-4, d_-5, d_-6, d_-7, d_-8, d_-10` with **constant** leading
coefficient 2 — an isomorphism, so no saturation can enter there), the three
`Q~` equations collapse to a 3-equation / 2-unknown residual system:

```
g1 = d1*dm1^2 + 2*dm1*dm3 + dm2^2
g2 = -d0*dm1^2 + 2*dm2*dm3
g3 = 2E - 6*d0*dm1*dm3 - 3*d0*dm2^2 - 6*d1*dm2*dm3 - 3*dm1^2*dm2
```

## 3. A CAS-free certificate (two lines, printable)

```
g3 + 3*d0*g1 + 3*d1*g2  =  2E - 3*dm1^2*dm2          [exact, no remainder]
```

so the system **forces** `2E = 3*d_-1^2*d_-2`. Then `g2 = 0` gives
`d_-3 = 3*d0*d_-1^4 / (4E)`. Substituting both into `g1` and clearing
denominators yields exactly **one half** of the paper's equation:

```
(8E^3 + 18*E*d1*dm1^6 + 27*d0*dm1^9) / 2
```

— the factor `1/2` being why the printed form is primitive. Both steps
re-verified symbolically (see §6).

**This replaces "using a CAS (for example Mathematica)" with a two-line
argument a reader can check by hand.**

## 4. Machine confirmation

Singular `eliminate()` on all nine original equations, 12-variable ring, `dp`:

```
=== ELIMINATION IDEAL J = I cap Q[d1,d0,dm1,E] ===
J[1] = 27*d0*dm1^9 + 18*d1*dm1^6*E + 8*E^3
number of generators: 1
timer (1/1000 sec): 0
```

Under one second. Cross-checks:

* `reduce(T, J) = 0` **and** `reduce(J, std(ideal(T))) = 0` — the elimination
  ideal is **exactly** `(T)`, principal.
* `reduce(T, std(I)) = 0` in the FULL ring — `T` is literally a polynomial
  combination of the nine equations, so it is an **unconditional** consequence;
  no genericity assumption is needed.
* `sat(I, d_-1)` is strictly larger than `I`, but has the **same** elimination
  ideal — the result is not hiding a saturation.
* `factorize(T)` — irreducible over Q.
* An end-to-end numeric check: random exact-rational points satisfying the
  residual system, back-substituted through `d_-4 .. d_-10`, satisfy all nine
  original equations, and `T` evaluates to 0 on every trial.

## 5. Three by-products worth reporting to the authors

1. **A LaTeX typo.** Tex lines **1680** and **2043** render `2 d_{-1}0`, which
   must be `2 d_{-10}`. Confirmed by regenerating the series coefficients
   independently.
2. **The `(7,21)` case at line 2056 is the same computation specialized**, not a
   second independent CAS run: with `E := F_-4 * y^10` and `F_-4 = (1/2) y^-1`,
   i.e. `E = y^9/2`, the master equation becomes
   `y^27 + 9 y^9 d_1 d_-1^6 + 27 d_0 d_-1^9`, identical to
   `\eqref{ecuacion principal para 7}`.
3. **The derivation path matches theirs.** At tex line 2066 the authors report,
   for `(7,21)`, `3 (d_-1)^2 d_-2 = y^9` and `2 d_-3 y^9 = 3 d_0 (d_-1)^4`. With
   `E = y^9/2` these are *precisely* the two intermediate relations of the
   certificate in §3. So this is a reconstruction of what they actually did, not
   a lucky alternative route.

Which of the nine equations are used is also **forced, not arbitrary**:
`(D~^2)_{-6}` and `(D~^2)_{-8}` introduce the new unknowns `d_-9, d_-11`, and
`(Q~)_{-3}` is the unique equation carrying `lambda`. The parallel passage at
line 2055 states the exclusion explicitly, confirming the same selection.

## 6. Scope — what is and is not verified

* **[verified]** That the printed equation is a correct and unconditional
  consequence of the nine printed equations, by elimination, by ideal
  membership, by a hand certificate, and numerically.
* **[verified]** The transcription of the nine equations: they were regenerated
  from the series definitions of `D~^2`, `D~^3`, `D~^-1` and diffed against the
  printed text; all nine match.
* **[NOT audited]** The derivation of those nine equations from the
  Jacobian-Conjecture setup, and the argument *after* `ecuacion principal` that
  produces the contradiction. Those remain unchecked here.
* The result is a **necessary condition** (an ideal contraction), which is the
  direction GGHV use it in. Correct as used.

## 7. Why this was worth doing

It converts a "trust the CAS" step in an unscrutinised preprint into a checkable
two-line argument, in the one paper our entire case-elimination campaign depends
on. It does **not** touch our own target: `A_0 = (8,28)` is the instance GGHV
explicitly leave open (line 268: *"we couldn't solve the corresponding system of
polynomial equations, thus it is left open"*), and Section 5 never touches it.
