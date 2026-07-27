# Adjudication: Helali's exact exclusion of the (72,108) frontier

**Date:** 2026-07-25. **Scope:** read-only audit of an external result. New files
only (`HELALI_ADJUDICATION.md`, `helali_adjudication_check.py`). No existing file
in this repo was edited; nothing was pushed anywhere.

**Subject.** Billel Helali, *Exact Computer-Assisted Exclusion of the (72,108)
Frontier in the Two-Dimensional Jacobian Problem*, repo commit
**`c530fe44e5f53b17840110931803e7c7c5a24cde`** (2026-07-21),
**doi:10.5281/zenodo.21479814**, release bundle
`jc2_72_108_exact_replay_v1.0.1.zip`
(`sha256 232204bdb598cc2ea0368e154c8573e18bbfdc69fa631c8878de4b884b38bb18`).

**Our side.** `d2_plane_72_108` at commit
**`ccf26d59c2f663bc20a06a36ba42e2fc9c89545c`** (the task named `d2ceb1f`; the
tree had advanced by one commit, `ccf26d5`, when this audit ran).

---

## 0. VERDICT — **SUBSUMES**

**His Case 1 covers our `sub1`, and his Case 2 covers our `sub2`. His reduction
is sound, it is complete in the sense our ledger needs, and we should adopt it.**

Every step of his Case-1 chain is one of: a linear solve over a *field*, a
division by a nonzero rational *integer*, an exhaustive factorization split, an
invertible affine change of variable, or an explicitly cased `h = 0` / `h != 0`
split. There is no point in the chain at which a stratum, a component, or a
degeneration can be silently dropped. I checked each of those steps
individually and re-verified all four of his exact certificates with arithmetic
I wrote myself.

**The 13-vs-7275 tension dissolves and is not evidence against him.** His 13
equations are what survives *complete, lossless* elimination of the bracket
layers `z⁴` down to `z^{-3}` — every band pair at every one of those layers,
nothing dropped — leaving 6 parameters over a **zero-dimensional degree-35**
first block. (He does not impose layers `z^{-4}` and below at all; that only
makes his claim harder to prove, never weaker.) Our 7275 are *valuation/degree strata* of a different parametrization,
produced by a filter that is explicitly an over-approximation of the solution
set (`CASCADE_ENGINE_REPORT.md:66-69`) and that never resolves the continuous
coefficients inside a stratum (each surviving state still carries 66 spare
scalars, `G4_ROW.md:301`). Alternative **(A)** of the task framing is correct:
7275 is not a competing count of the same thing.

I found **no uncovered state**, and I now believe none exists, for a structural
reason given in §6 rather than for want of looking.

**Two things to hand him, both presentational, neither affecting the result:**
a display bug in the paper's equations (J1) and (J0) (§7.1) and a stale section
in the bundled internal report (§7.2). His *code* is right in both cases.

**One thing we can give him:** his own open issue #1 — "verify the Proposition
4.3 transcription bridge", which he calls his main remaining risk — is
**discharged on the transcription half**, machine-checked against the arXiv
source in `helali_adjudication_check.py` (T1). See §8.

---

## 1. What I actually ran

Everything below is reproducible; see §9.

| # | check | result |
|---|---|---|
| R1 | His full replay `verify_all.sh`, on Windows/native CPython **3.10.6** with a freshly installed `gmpy2 2.3.1` — a different OS, Python (his reference used 3.14), and library build from his Apple-Silicon run | **`JC2_72_108_EXACT_REPLAY_PASS`** |
| R2 | Deleted `case1_checkpoint.pkl` and **regenerated the whole Case-1 band descent from scratch** (`case1_descent_checkpoint.py`) | checkpoint **byte-for-byte identical** (`sha256 2dcf13d9…0240ba`); all 13 residual equations identical line-for-line |
| R3 | **Independently re-derived** the 6 first-block residual generators from `J4` in my own script (`fb_gen.py`), verifying every elimination pivot | identical to his `firstblock_Q_exact.sing` ideal, all 11 pivots are the integers 1,3,5,…,21 |
| R4 | **Independently recomputed** the first-block Gröbner basis over **Q** (`modStd(I,1)` — Singular's *exact* modular standard basis — then `fglm`), which he never reported dimension data for | **`DIM=0`, `VDIM=35`**, lex basis **byte-identical** to his `L[1]…L[6]` |
| R5 | Checkpoint audit: 13 equations, `s`-split exhaustiveness, `S != 0`, `E6 − λE2 = S·L²`, the discarded equations | all pass (`helali_adjudication_check.py` T6–T8) |
| R6 | **Independent replay of the hard 89 MB identity** `h = ΣT_iE_i`, written on `python-flint` doing straight polynomial arithmetic in `L[h,u1,u2]` — *not* his `gmpy2` 2010×1925 scalar-row path | 385 multiplier terms parsed, **13,410 `L`-products**, residual is exactly the monomial `h` |
| R7 | Independent replay of the three unit certificates (Case 2, and `h=0` on both branches), again on my own `flint` arithmetic | all three give exactly `1` |
| R8 | Transcription bridge against `paper_src/2204.14178.tex` | all four corner sets match character-for-character; `L^{(1)} = K[x,x^{-1},y]` confirmed |

`helali_adjudication_check.py` bundles R3, R5, R6, R7, R8 plus the band and
layer-completeness census: **46 checks, 0 failures.**

---

## 2. PROVED / CHECKED / INFERRED

### PROVED (mathematical argument, verified here, not machine-dependent)

- **P1. The normalization `a_1 = a_8 = c_8 = 1` is a genuine WLOG.** Under
  `P̃ = ρP(λx,μy)`, `Q̃ = σQ(λx,μy)` the bracket scales by `ρσλ³μ`, so `σ`
  restores `[P̃,Q̃] = x²`; the three vertex coefficients transform by the
  monomial matrix with rows `(ρ,λ,μ)`-exponents `(1,1,0)`, `(1,8,14)`,
  `(1,8,16)`, whose determinant is **14 ≠ 0**. Over an algebraically closed
  field a monomial map with nonzero exponent determinant is surjective on
  `(K*)³`, so `(a_1,a_8,c_8)` can be moved to `(1,1,1)`. The scaling preserves
  the Newton polygons exactly. `helali_adjudication_check.py` T4.
- **P2. The `J4` elimination is division-free.** Solving the eleven `D`
  coefficients triangularly divides only by the *integers* `2m−3` for
  `m = 2..12`, i.e. `1,3,5,…,21`. Hence
  `J4 ⟺ {d_m = f_m(a)} ∧ {res = 0}` with **no component lost**, `res` being 6
  honest polynomials in `a_2..a_7`. T5.
- **P3. Every layer he uses is complete.** In `t = xy²`, `z = y^{-1}` a
  monomial `x^i y^j` is `t^i z^{2i−j}`, and layer `k` of `[P,Q]` collects every
  band pair `(i,j)` with `i+j−1 = k`. At layer `k ≤ 1` exactly two bands are
  new — `P_{k−2}` and `Q_{k−1}` — and both appear as *unknown columns* in
  `case1_descent_checkpoint.py:38-57`; every other pair is already solved. T3
  enumerates all pairs for `k = 4,3,2,1,0,−1,−2,−3` and finds none missing.
- **P4. Using only layers `4 … −3` is sound.** Layers `−4 … −21` are simply
  not imposed. Dropping necessary conditions enlarges the solution set, so an
  exclusion proved on the subset is *a fortiori* valid. Likewise he never
  imposes the vertex-nonvanishing conditions at `(0,8)`, `(0,12)`, `(12,21)`,
  `(12,24)`, and he discards the (bracket-irrelevant) constant terms of `C` and
  `G`: all relaxations, all in the safe direction.
- **P5. The `s = ±c` split is exhaustive.** The first compatibility equation
  `q0` depends on `s` **alone**, has `s`-degree 4, and equals
  `κ(s²−c²)²` exactly with `κ ≠ 0` and `c² ≠ 0`. So the zero set of `q0` in `s`
  is exactly `{c, −c}` — **there is no `c = 0` component to drop**, and the
  `c ≠ 0` claim in his paper is a consequence, not an assumption. T7.
- **P6. `S ≠ 0` is established, not assumed.** `S` is read off as the `r²`
  coefficient of `E6 − λE2` and the identity `E6 − λE2 = S·L²` with
  `L = r + αh + βu_2 + γ` is *checked exactly* — it is an `assert` on the
  execution path of `case1_cascade_machine.py`, which `verify_all.sh` runs. `K`
  is a field, so `S·L² = 0 ⟺ L = 0`, eliminating `r` with **no division by a
  variable expression**. T8c/T8d.
- **P7. Nothing is lost when he trims 13 equations to 7.** After the `r`
  substitution, four of the thirteen vanish identically and the two live ones
  he drops (indices 6 and 12) are exact `K`-scalar multiples of a kept one
  (`13/9` resp. `4/9` on branch 1; `−13/9` resp. `4/9` on branch 2). T8e.
- **P8. The only division by a variable in the whole pipeline is `u_3 = N/h`,
  and `h = 0` is covered separately.** `derive_hne0.py:41-48` homogenizes by
  `h^d` and strips a common `h^{m}`, which is valid on `h ≠ 0`; the `h = 0`
  locus is excluded by its own pre-division unit certificate built from the
  seven equations *before* that division. This is exactly the place where a
  reduction of this shape usually leaks, and he plugged it.
- **P9. The branch involution is a legitimate transport.** `φ: h↦h,
  u_1↦−u_1, u_2↦−u_2` is a ring automorphism; the branch-2 system is checked to
  be `s_i·φ(E_i^{(1)})` term-by-term, so `T'_i = s_iφ(T_i)` gives
  `ΣT'_iE_i^{(2)} = φ(h) = h`. `verify_hne0_branch_symmetry.py:330-334`.

### CHECKED (machine-verified here, exactly, in characteristic zero)

- **C1.** `JC2_72_108_EXACT_REPLAY_PASS` on an independent platform (R1).
- **C2.** Byte-identical regeneration of the Case-1 descent checkpoint (R2).
- **C3.** First block is **zero-dimensional of degree exactly 35** over Q, lex
  basis byte-identical to his (R4). This is the check that closes coverage of
  the first block: combined with his own verification that the 35 conjugate
  points satisfy all six residuals, `V(res) = V(H, five linear relations)`, so
  `K = Q[u]/(H)` is the **entire** first-block solution field.
- **C4.** `1 = ΣT_iR_i` for Case 2 over `K`, replayed on my own arithmetic.
- **C5.** `1 = ΣA_iE_i|_{h=0}` on both branches over `K`, ditto, and the seven
  generators regenerate from the descent.
- **C6.** `h = ΣT_iE_i` over `L = Q[w]/(w⁵−w⁴+3w³+3w²+26)`, replayed by direct
  polynomial multiplication (13,410 `L`-products), residual exactly `h`.
- **C7.** GGHV Prop 4.3 transcription: exact, all four corner sets.

### INFERRED (believed, resting on something not re-derived here)

- **I1.** GGHV Proposition 4.3 is *exhaustive* — i.e. its own proof is correct.
  Neither program has audited that. It is the shared conditionality and he
  states it correctly.
- **I2.** Singular's `modStd(I,1)` / `fglm` are correct. I ran them on input I
  generated myself and got byte-identical output to his, and `modStd`'s
  `exactness = 1` default is documented as computing a standard basis "for
  sure", but this is still one CAS and one algorithm family. A plain Buchberger
  `std` over Q on the same ideal was launched as a further hedge and **timed
  out** (`EXIT=124` at 50 minutes — a timeout, never a verdict); see §9.1. This
  is the single largest remaining exposure in the whole audit, and it is small.
- **I3.** `python-flint`/FLINT rational arithmetic is correct. My replays and
  his `gmpy2` replay agree, and they use different libraries for the big
  integers, which is decent independent evidence.
- **I4.** The `nfmodStd` "UNIT" outputs quoted in his bundled
  `JC2_72_108_EXACT_REPORT.md` §6.4 — I did **not** rerun them, and **nothing
  in the result depends on them** (see §7.2).

---

## 3. The dictionary — his coordinates versus ours

This is the single easiest way to reach a false verdict, so it is stated
explicitly and *not* used to compare anything numerically.

| | his | ours |
|---|---|---|
| distinguished symbol `t` | `t = x·y²` — a **coordinate on the torus**, paired with `z = y^{-1}`, `[t,z] = −1`, `x² = t²z⁴` | `t = y + 1` — a **place** in `K[y]`, the linear factor of `C₄ = y⁷(y+1)` at which valuations are taken (`divisor_filter.py:23`) |
| ambient ring | `K[x,x^{-1},y]` directly | `K[y,C₄^{-1}]((x^{-1}))` (`STATE.md:27`) |
| unknowns | Newton-support coefficients grouped into `z`-bands `P = Σ_k P_k(t)z^k`, `Q = Σ_k Q_k(t)z^k` | `D`-transform window symbols `d₂,d₁,σ,e = D̃₂,D̃₁,D̃₀,D̃_{-1}` and the cascade chain `g₁..g₇` (`STATE.md:32-40`) |
| what is recorded | actual coefficients, in `K` or in a number field | `y`-**degrees**, zero flags, and valuations — an upper envelope, never individual coefficients (`ALOK_CROSSCHECK.md:121-123`) |
| granularity | points of a variety | strata over 66 continuous spare scalars (`G4_ROW.md:301`) |

**The two `t`'s are different objects and must never be identified.** Our repo
has never used band coordinates at all; the only prior mention of them is our
own summary of *his* method (`LANDSCAPE_2026_07.md:269`).

What *is* shared, and is the whole basis of the comparison, is the ambient
object: both programs study

> `S := { (P,Q) ∈ (K[x,x^{-1},y])² : [P,Q] = x², N(P) = conv{(0,0),(1,0),(8,14),(8,16),(0,8)}, N(Q) = conv{(0,0),(2,1),(12,21),(12,24),(0,12)} }`

— his **Case 1**, our **`sub1`**, GGHV **Prop 4.3 configuration (1)**. Vertex
for vertex, same bracket normalization: `ALOK_CROSSCHECK.md:62-80`,
`STATE.md:18-21`, and now machine-checked against the arXiv source in
`helali_adjudication_check.py` T1.

---

## 4. Case 2 / `sub2` — we agree, and the agreement is real

Confirmed. In Case 2 the polygons have **no negative `z`-bands at all**
(`helali_adjudication_check.py` T2e: the minimum band index is 0 for both `P`
and `Q`), so `P = Az²+Bz+C`, `Q = Dz³+Ez²+Fz+G` is *exact* and `J4…J0` is the
**complete** coefficient system, not a truncation. The elimination is then
pure linear algebra over the field `K`: `J3`,`J2` solve `B,E,C,F` with three
free parameters `r,s,h`; `J1` solves `G_1..G_12` with no freedom left and seven
compatibility conditions; `J0` contributes 18 more; four of the 25 already
generate the unit ideal, and `1 = T₁R₁+T₂R₂+T₃R₃+T₄R₄` replays exactly.

We reached the same conclusion — `sub2` EMPTY — by positive-slice liftability
(`e | Phi` → T2 divisor → SPINE → positive slice, `CURRENT_STATUS.md:140`;
graded "exact-checked + independent verifier `positive_slice_verify.py`" at
`CURRENT_STATUS.md:179`, `positive_slice.py` 63/63 plus `_verify` 79/79 at
`SESSION_HANDOFF.md:90`), a completely different mechanism.
**Two independent methods, same answer.**
That is mutual corroboration and it is meaningful evidence that his
transcription and normalization are faithful, because a mis-transcribed polygon
would be overwhelmingly unlikely to produce an empty system by *our* route as
well.

---

## 5. The six specific checks, in order of load

1. **Does his variable change relate to ours?** — §3. Established as a
   dictionary; deliberately not used for numerical comparison. His `t` ≠ our
   `t`.
2. **Is "normalized nonzero vertex coefficients" a genuine WLOG?** — **Yes.**
   P1. Vertices are nonzero by definition of the Newton polygon; the *scaling
   to 1* is the real question and it is legitimate because the exponent matrix
   has determinant 14, so the torus acts transitively on the three vertex
   coefficients over an algebraically closed field. His explicit formulas
   (`μ² = a₈/c₈`, `λ⁷ = (a₁/a₈)μ^{-14}`, `ρ = (a₁λ)^{-1}`, `σ = (ρλ³μ)^{-1}`)
   are a constructive witness and they check out. Nothing our cells allow to
   degenerate is quietly assumed nonzero: the *other* vertices are never
   assumed nonzero at all.
3. **Is the `s = c` / `s = −c` split exhaustive?** — **Yes.** P5. The
   factorization is not merely quoted from Singular: I recomputed `q0` from the
   regenerated checkpoint, confirmed it is univariate in `s` of degree 4 with
   only even powers, and verified `q0 = κ(s²−c²)²` as an exact identity in `K`,
   with `κ ≠ 0` and `c² ≠ 0`. **No third component, and in particular `c = 0`
   is not a dropped case — it is refuted.** `c ≠ 0` is derived, not assumed.
4. **`E₆ − λE₂ = S·Λ²` with `S ≠ 0`.** — **Established.** P6. It is an
   executed assertion, and I re-derived `λ`, `S`, `α`, `β`, `γ` and re-checked
   the square identity independently on both branches. (Aside: his report
   assigns `λ = 13/9` and `−13/9` to the branches; my recomputation gets them
   the other way round. That is a labelling convention, not a discrepancy —
   both values occur, one per branch.)
5. **The hard certificate.** — **Verified independently.** The mod-71 minor
   selection is *sound as he uses it*: a square rational matrix whose
   determinant is nonzero mod 71 is nonsingular over Q, full stop; and in any
   case the minor is only the construction mechanism. Soundness comes from the
   recovered solution being substituted back into **all 2010** scalar rows,
   which his `gmpy2` verifier does. I did not rely on his verifier: I wrote a
   different one that never forms the scalar system at all, multiplying the four
   multipliers by the four generators directly in `L[h,u₁,u₂]` on
   `python-flint`, and the product is **exactly the monomial `h`** — 13,410
   number-field products, one nonzero result monomial, coefficient 1. The three
   unit certificates replay the same way, all giving exactly `1`.
6. **Coverage.** — §6.

---

## 6. THE MAIN EVENT — coverage

**Question.** Does every one of our 7275 `sub1` states correspond to a point
his thirteen equations constrain?

**Answer.** The right framing is not a state-to-point map — there is provably no
faithful one, and `ALOK_CROSSCHECK.md:113-128` already recorded why (different
lattices, a nonlinear convolution between them, degrees rather than supports).
The right framing is that **both programs parametrize subsets of the same set
`S`, and his parametrization is exhaustive.**

His Case-1 chain, read as a statement about `S` itself:

```
S ≠ ∅
 ⟹ (P1, torus action)  S contains a point with a₁ = a₈ = c₈ = 1
 ⟹ (P3)                its bands satisfy layers J4, J3, J2 and layers 1,0,−1,−2,−3
 ⟹ (P2)                its (a₂..a₇) satisfy the 6 first-block residuals
 ⟹ (C3, DIM=0 VDIM=35) (a₂..a₇) is one of exactly 35 Galois-conjugate points,
                        i.e. lies over K = Q[u]/(H)
 ⟹ (field linear algebra, no branching)
                        B,E,C,F,G and the five negative bands are determined by
                        6 parameters r,s,h,u₁,u₂,u₃ over K, subject to 13 equations
 ⟹ (P5)                s = c or s = −c, exhaustively
 ⟹ (P6)                r is eliminated; (P7) 7 equations remain in (h,u₁,u₂,u₃)
 ⟹ h = 0 or h ≠ 0
      h = 0 : unit certificate  ⟹ contradiction        (C5)
      h ≠ 0 : u₃ = N/h, descend to L, hard identity forces h = 0 ⟹ contradiction (C6)
 ⟹ S = ∅.
```

Every arrow is either a proved implication (P1–P9) or an exactly checked
identity (C1–C6). **No arrow conditions on genericity, and no arrow divides by
a quantity that has not been shown nonzero.** That is the coverage proof: it
does not need to know anything about our coordinates, because it quantifies
over all of `S`.

Consequently `S = ∅`, and therefore **every one of our 34 cells / 314 flag
cases / 7275 states is empty**, as is the whole 171-cell / 44117-state raw
Phase-D universe above it, as are the 5 cells / 1664 (ON) or 2124 (OFF) states
left after `a_t = 9`.

**Why 13 and 7275 are not in competition.** Three reasons, each sufficient:

1. **Different kind of object.** A state of ours is "a residual degree
   assignment `(deg d₂, deg d₁, deg σ, deg e)`" inside a cell `(a_t, b, branch)`
   and flag case (`FRONTIER_REBUILD.md:34`) — a *stratum*, carrying 66
   continuous spare scalars before collapse and 18 after (`G4_ROW.md:301`).
   His 13 equations cut out a *variety*. A stratum count and an equation count
   are not commensurable quantities.
2. **Ours over-approximates by design.** "A kill means no consistent valuation
   profile exists; a survivor comes with an explicit witness profile"
   (`CASCADE_ENGINE_REPORT.md:66-69`). Survivors are candidates the valuation
   filter cannot refute — not solutions. Emptiness of `S` is entirely
   consistent with 7275 survivors of a necessary-condition filter.
3. **He resolves a rigidity we never computed.** His first block pins
   `(a₂..a₇)` to a **finite** set of 35 points before any of the rest happens.
   Our repo has no analogue: it contains **no statement bounding the dimension
   or degree of the `sub1` solution variety**, and `G4_ROW.md:301` records
   explicitly that a full Gröbner sweep over the 314 flag cases "was not
   attempted". His degree-35 rigidity is genuinely new information relative to
   our ledger, and it is exactly the thing that collapses 7275 strata worth of
   apparent freedom.

**What I looked for and did not find.** A gap would have to be a place where a
solution of `S` fails to reach one of his 13 equations. There are exactly five
candidate leak types in a reduction of this shape, and each is closed: torus
normalization (P1), division by a possibly-vanishing pivot (P2 — integers only;
and everything after the first block is linear algebra over a *field*, where
pivots are zero or invertible with no third case), a non-exhaustive case split
(P5), a lost component under saturation (P8 — the one division by a variable,
with its `h = 0` companion certificate), and an incomplete first-block ideal
(C3 — `DIM=0`, `VDIM=35`, so the 35 known conjugate solutions are *all* of
`V(res)`).

---

## 7. Two things to hand him

Both are presentational. **Neither changes any result, and in both cases his
code is correct.** They are worth sending because a referee will hit them.

### 7.1 The displayed equations (J1) and (J0) are incomplete *for Case 1*

In `paper/jc2_72_108_exact_exclusion.tex` §"Exact Laurent reduction" (and in
`JC2_72_108_EXACT_REPORT.md` §2), the five identities J4…J0 are presented as
"the coefficient identities of `[P,Q] = x²`" for *both* cases. For **Case 2**
that is exactly right, because there `P` and `Q` genuinely have no bands below
`z⁰` (checked: `helali_adjudication_check.py` T2e). For **Case 1** the extra
vertices `(0,8)` and `(0,12)` create bands `P_{-1},…,P_{-8}` and
`Q_{-1},…,Q_{-12}`, and the layers `z¹` and `z⁰` pick up further pairs.

Using his own band formula `[P_i z^i, Q_j z^j] = z^{i+j-1}(iP_iQ_j' − jP_i'Q_j)`,
the correct Case-1 identities are

```
(J1)  2AG' + (BF' − B'F) − 2C'E  −  P_{-1}D' − 3P_{-1}'D  = 0
(J0)  BG' − C'F  +  2A Q_{-1}' + A' Q_{-1}
                 −  P_{-1}E' − 2P_{-1}'E
                 −  2P_{-2}D' − 3P_{-2}'D               = 0
```

where `P_{-1}(t) = Σ_{i=0}^{7} p_{i,2i+1} t^i` is the `z^{-1}` band of `P`
(lattice points `(i, 2i+1) ∈ N(P)`, `i = 0..7`; e.g. `(0,1)`, which sits on the
edge from `(0,0)` to `(0,8)`), and similarly `P_{-2}` (degree ≤ 6) and `Q_{-1}`
(degree ≤ 11).

**His implementation already does this correctly.** In
`case1_descent_checkpoint.py:36-68` the loop over `k ∈ [1,0,−1,−2,−3]` sets
`r = k−2`, `s = k−1` and introduces `P_r` and `Q_s` as **new unknown columns**
at each layer — for `k = 1` those are precisely `P_{-1}` and `G = Q_0`, and the
`P_{-1}` columns are filled by `rows[m+j-1][col] += D_j·(r·j − 3m)` with
`r = −1`, which is exactly the missing `−P_{-1}D' − 3P_{-1}'D`. So the
thirteen compatibility equations are derived from the *complete* layer
equations. Only the write-up truncates.

Suggested fix: one sentence saying that J4…J0 as displayed are the complete
system for Case 2, and that in Case 1 layers `z¹` and below additionally
involve the negative bands, which the descent introduces two at a time. His
internal report already says something close to this in §6 ("The extra negative
`z`-bands are introduced in complete slices"); the paper does not.

### 7.2 `JC2_72_108_EXACT_REPORT.md` §6.4 is stale relative to the paper

The bundled report closes Case 1 by quoting `nfmodStd` returning `{1}` for the
seven-equation systems. His own `JC2_HNE0_PROGRESS.md` is appropriately
sceptical about exactly that kind of output, and the **paper and the replay do
not use it** — they use the `h = 0` unit certificates plus the hard membership
identity, which is the stronger and cleaner argument. A reader who takes §6.4
as the load-bearing step will misjudge the evidential status of the result in
*his own disfavour*. Worth a note that §6.4 is superseded.

### 7.3 Minor

- Report §6.2 says "Three remaining equations are exact constant multiples of
  another equation". After the `r`-substitution I count four of the thirteen
  vanishing identically (indices 0,1,3,7) and **two** live ones proportional to
  a kept equation (indices 6 and 12). 13 − 4 − 2 = 7, which is the stated
  seven. Cosmetic.
- `SHA256SUMS.txt` lists `jc2_72_108_exact_exclusion.tex` as
  `449f2cce…88f0d`, but the copy at `paper/` in commit `c530fe4` hashes to
  `bae3284e…425322`. Presumably the manifest refers to the release asset and
  the repo copy has since been edited. Harmless, but a reviewer verifying
  hashes will trip on it.

---

## 8. What our work adds that his does not

Stated without inflation. If his result stands — and on this audit it does —
most of our `sub1` machinery becomes *vacuously* true, and that is fine.

1. **We discharge the transcription half of his open issue #1.** This is the
   most substantive thing we can give him. `T1` in
   `helali_adjudication_check.py` parses the `[Case (8,28)]` proposition out of
   `paper_src/2204.14178.tex` and compares all four corner sets to his
   transcription: they match exactly, as does the bracket normalization
   `[P,Q] = x²` and the ring `L^{(1)} = K[x,x^{-1},y]` (GGHV line 656). Our
   repo had already recorded this as CONFIRMED VERBATIM (`STATE.md:182-186`,
   pinned in `paper_src/upstream_facts.json`) with an independent
   vertex-for-vertex census (`ALOK_CROSSCHECK.md:62-100`,
   `caps_audit.py:109-113`, `field_scope_audit.py`). What remains open in his
   issue #1 is only the *exhaustiveness* of GGHV's own proof, which is a
   question about a published paper and is open for both of us.
2. **Independent corroboration of Case 2 by a different mechanism.** §4. Two
   unrelated methods agreeing on `sub2 = ∅` is the strongest kind of evidence
   that the shared transcription is faithful.
3. **A branch-complete ledger as a hedge.** Our `sub1` universe (171 cells →
   34 cells / 314 flag cases / 7275 states, `FRONTIER_REBUILD.md:74,88`,
   re-derived at `G4_ROW.md:203,371`) is an independently constructed
   over-approximation of the same set. If his reduction is ever found to have a
   gap, the ledger is where the search resumes, and it is branch-complete by
   construction rather than by elimination.
4. **The `a_t = 9` two-sided sandwich, honestly graded.** `a_t := v_{y+1}(e)`.
   `a_t ≥ 9` is INDEPENDENTLY AUDITED (`slice_obstruction_audit.py`, 56/56,
   separately authored); `a_t ≤ 9` is exact-checked by two different mechanisms
   but by the same author, i.e. corroborated and **not** audited
   (`CURRENT_STATUS.md:164-181`). It is deliberately not wired into the
   frontier. Post-adoption it becomes a fact about our reduction rather than
   about the conjecture; it is not something to claim priority on here.
5. **Audit discipline that caught things.** The layered PROVED/CHECKED/INFERRED
   grading, the "no unaudited kill in the frontier" rule, and the habit of
   re-deriving from scratch are what produced R2–R4 above; that methodology is
   transferable and is what let this adjudication reach a verdict in one pass.

**What we should keep:** items 1–3 and 5. **What we should stop spending
effort on:** further `sub1` kill work aimed at closing (72,108). The right move
now is to say publicly and without hedging that Helali got there, cite
`doi:10.5281/zenodo.21479814`, and offer him §7 and §8.1 directly.

---

## 9. Reproduction

```bash
# tier 1 only (sympy; runs against our repo's paper_src/)
python helali_adjudication_check.py

# tier 1 + tier 2 (needs python-flint; point at his extracted bundle)
python helali_adjudication_check.py /path/to/exact_replay
```

Expected tail: `CHECKS_RUN=46  FAILURES=0` / `HELALI_ADJUDICATION_PASS`.
Tier 2 takes about 2 minutes, dominated by parsing the 89 MB certificate.

To reproduce R1/R2/R4 (not part of the checker, they mutate his tree):

```bash
# R1 his own replay
python -m pip install gmpy2 numpy python-flint sympy
cd exact_replay && PYTHON=$(which python) bash ./verify_all.sh     # -> JC2_72_108_EXACT_REPLAY_PASS

# R2 regenerate the Case-1 descent from scratch
cp -r exact_replay rebuild && rm rebuild/case1_checkpoint.pkl
cd rebuild && python case1_descent_checkpoint.py
sha256sum case1_checkpoint.pkl        # 2dcf13d924530cdc9a8728e943efdc73d003ce1c187d5cec273f6f701e0240ba

# R4 independent first-block Groebner basis over Q
#   emit the 6 residual generators, then in Singular:
#     ideal G=modStd(I,1); dim(G); vdim(G); fglm to lp
#   -> DIM=0  VDIM=35  and L[1..6] byte-identical to firstblock_Q_exact.out
```

### 9.1 One hedge that did not finish — and why it does not matter

A plain Buchberger `std` over **Q** on the same six first-block generators — a
deterministic, non-modular confirmation of R4 — was launched under a 50-minute
cap and was **killed by the timeout (`EXIT=124`, no output, RSS 1.3 GB and
still climbing)**. Per local discipline, **124 is a TIMEOUT and is never a
verdict**; nothing about the first block is called into question by it, and the
run is best read as "plain Buchberger is the wrong algorithm here", which is
presumably why he reached for `modStd` in the first place.

**It is a hedge, not a dependency.** R4 already establishes `DIM=0`,
`VDIM=35` and a byte-identical lex basis via `modStd(I,1)` — Singular's exact
modular standard basis, default `exactness = 1`, documented as computing a
standard basis "for sure" rather than "with high probability" — followed by
`fglm`, which is deterministic linear algebra on a zero-dimensional quotient.
Crucially I ran that on **input I generated myself** (`fb_gen.py`, R3) and it
reproduced his `L[1]…L[6]` byte for byte, which is a strong consistency signal
independent of who typed the ideal.

The residual exposure is therefore exactly `I2`: one CAS, one algorithm family.
Anyone who wants to close it can rerun with a larger budget, or on a second
system:

```singular
ring rd=0,(a2,a3,a4,a5,a6,a7),dp;
ideal I = <the six generators emitted by helali_adjudication_check.py T5>;
option(redSB); ideal G = std(I);          // budget hours, not minutes
dim(G); vdim(G);                          // expect 0 and 35
```

If that ever disagrees with `modStd(I,1)`, it is a Singular bug that affects
his computation and my reproduction of it equally, and should be reported
upstream as such.

---

## 10. Citations

- **His work.** B. Helali, *Exact Computer-Assisted Exclusion of the (72,108)
  Frontier in the Two-Dimensional Jacobian Problem*, 21 July 2026. Repo
  `bilLkarkariy/jc2-72-108-exact-certificates`, commit `c530fe4`.
  **doi:10.5281/zenodo.21479814.** Replay bundle
  `jc2_72_108_exact_replay_v1.0.1.zip`,
  `sha256 232204bdb598cc2ea0368e154c8573e18bbfdc69fa631c8878de4b884b38bb18`.
  Hard certificate `hard/h_certificate_exact.txt`,
  `sha256 0e48ffab32469ef8405a6945b16cf1521ddeb3c592ae4e5051968110a4dc656a`.
  Case-1 descent checkpoint,
  `sha256 2dcf13d924530cdc9a8728e943efdc73d003ce1c187d5cec273f6f701e0240ba`
  (regenerated here from scratch).
- **Upstream.** J. A. Guccione, J. J. Guccione, R. Horruitiner, C. Valqui,
  *Increasing the degree of a possible counterexample to the Jacobian
  Conjecture from 100 to 108*, arXiv:2204.14178. Proposition 4.3
  `[Case (8,28)]`, `paper_src/2204.14178.tex:1001-1006`; `L^{(1)} = K[x,x^{-1},y]`
  at `:656`.
- **Ours.** `d2_plane_72_108` at commit `ccf26d5`. Load-bearing files cited
  inline: `ALOK_CROSSCHECK.md`, `STATE.md`, `FRONTIER_REBUILD.md`,
  `CASCADE_ENGINE_REPORT.md`, `G4_ROW.md`, `CURRENT_STATUS.md`,
  `divisor_filter.py`, `caps_audit.py`.
- **This audit.** `HELALI_ADJUDICATION.md` + `helali_adjudication_check.py`,
  both new; no existing file touched.
