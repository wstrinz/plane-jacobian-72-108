# MONOMIAL_LENS — automatic detection of latent quotient coordinates

> **STATUS (2026-07-24): TOOL LANDED, FOUR REGRESSIONS PASS, AND THE R9 ANSWER
> IS A PROVED "NO".** `monomial_lens.py --quiet` exits 0 on all four required
> regressions. The GGHV22 tex-1689 collapse (`E = F_-4·C_3^23`) is rediscovered
> mechanically from the nine printed equations, and the lens additionally proves
> it is the *unique* such coordinate. The R9 z=1 dm4-eliminated H-system — and
> the pre-elimination full G-system, and the alt bridge, which is byte-identical
> to R9 — has **exponent-lattice torus dimension 0**, which is not "we searched
> and found nothing" but a theorem: *no* monomial change of coordinates, not
> even a Laurent one, can lose a single variable there.
>
> The same holds at the **scalar** level (`--deep`): 28 spare unknowns, and 37
> with the state parameters declared too — torus dimension 0 in both.
>
> **So the line-1689 collapse was a feature of that particular system, not a
> recurring mechanism, and the R9 wall is not a missed quotient coordinate.**
>
> New files (nothing existing modified): `monomial_lens.py`,
> `monomial_lens_results.json`, `monomial_lens_run.log`, `scalar_r9_lens.log`,
> this doc.

---

## 0. What was built and why

`GGHV_LINE1689.md` records that GGHV22's "using a CAS (for example Mathematica)"
step collapses to a hand computation once you notice that `C_3` and `F_-4` occur
in the nine printed equations **only** through the single product
`E := F_-4·C_3^23`. That observation was made by eye. The obvious question is
whether it is an instance of a *mechanism* we could be missing elsewhere — in
particular on our own wall, the R9 z=1 state.

`monomial_lens.py` answers that question mechanically. Given a list of sympy
polynomials and a declared variable list it:

* builds the exponent matrix over the declared variables;
* **LEVEL A** finds every maximal set `S` of variables whose exponent *columns*
  are pairwise proportional, and emits `z = ∏_{v∈S} v^{α_v}`;
* **LEVEL B** finds degenerate variable subsets (`rank(M_S) = r < |S|`) via the
  minimal supports of low-support lattice relations, and tries to realise each
  with `r` new coordinates;
* **LEVEL C** solves for the weight vectors making each input polynomial
  weighted-homogeneous;
* reports the **exponent-lattice rank** and hence the *torus dimension*
  `#vars − rank`, an upper bound on how many variables any monomial change of
  coordinates could ever remove;
* prints before/after statistics and an explicit **soundness obligation** for
  each proposal, and *runs* the machine-checkable parts of it;
* **refuses**, naming a witness monomial, whenever a candidate variable also
  occurs outside the proposed combination.

Two structural facts make the search complete rather than heuristic:

1. `rank(M_S) = 1` **iff** the exponent columns of the variables in `S` are
   pairwise proportional. Proportionality is an equivalence relation, so
   bucketing columns by their primitive direction enumerates every maximal
   rank-1 merge and misses none. LEVEL A is exhaustive by construction, not by
   exhaustion.
2. If all monomials lie in the image of a monomial map from `r < n` variables,
   the exponent rows lie in a rank-`r` sublattice, so some nonzero weight vector
   annihilates every monomial. Contrapositive: **exponent-lattice rank = n ⇒ no
   monomial coordinate change of any kind loses a variable.** This is what makes
   the R9 negative a proof rather than a failed search.
3. The same argument localises: if `rank(M_S) < |S|` for a subset `S`, extending
   the witnessing weight vector by zero gives a global null vector *supported
   inside* `S`. So every degenerate subset contains the support of a nonzero
   null-space vector — which is why LEVEL B enumerates minimal supports of
   lattice relations instead of scanning `C(n,k)` subsets, and why torus
   dimension 0 lets it return "empty" without scanning at all.

## 1. The soundness lemma (and it is a two-sided one)

Let `V` be the declared variables, `S ⊆ V`, `α ∈ Z_{>0}^S` primitive,
`W = V \ S`, and

    φ : A := Q[W][z] → B := Q[V],   z ↦ x^α = ∏_{v∈S} v^{α_v},   φ|_{Q[W]} = id.

* **(H1) EXCLUSIVITY.** Every monomial of every `f_i` restricts on `S` to
  `k_m·α` with `k_m ∈ Z_{≥0}`.
* **(H2) NONDEGENERACY.** `α_v ≥ 1` for every `v ∈ S`.

**LEMMA L.** Under (H1),(H2), with `I = (f_1..f_m) ⊆ B`, `g_i := φ^{-1}(f_i)`,
`J = (g_1..g_m) ⊆ A`, for every `U ⊆ W`:

    J ∩ Q[U]  =  I ∩ Q[U].

*Proof.* `φ` is injective on monomials (`α ≠ 0`, `S ∩ W = ∅`), so `g_i` is well
defined and `φ(J) ⊆ I` — the "⊆" direction, which is all a necessary-condition
elimination needs. For "⊇": by (H2) every `m ∈ N^S` factors uniquely as
`m = k·α + m'` with `k = min_v ⌊m_v/α_v⌋` and `m' ≱ α`, so `B` is a **free**
`φ(A)`-module on the basis `{x^{m'} : m' ∈ N^S, m' ≱ α}`, which contains `1`.
Let `π : B → φ(A)` be the `φ(A)`-linear projection onto the basis element `1`.
If `h ∈ I ∩ Q[U]` then `h = Σ c_i f_i`, and `h` lies in `φ(A)·1`, so
`h = π(h) = Σ π(c_i) f_i ∈ φ(J)`. As `φ` is injective and `h ∈ Q[U] ⊆ A`,
`h ∈ J ∩ Q[U]`. ∎

So a LEVEL A merge is **lossless in both directions**: elimination downstairs is
exactly elimination upstairs. For LEVEL B (`r ≥ 2` new coordinates) only "⊆" is
claimed — freeness of `B` over `φ(A)` can fail for a rank-`r` monoid — and the
tool prints the converse as an **undischarged obligation** rather than asserting
it. That distinction is enforced in the emitted obligation text, not just here.

The tool never applies a merge it has not checked: `apply_merge()` raises on an
(H1) failure, and REG3 exercises that path deliberately.

## 2. Results table

Run: `python monomial_lens.py` (full report), `--quiet` (exit code only),
`--deep` (adds the scalar R9 pass). Machine: Windows, Python 3.10, sympy 1.14.

| system | vars | distinct monomials | max total deg | LEVEL A merge | exponent-lattice rank / torus dim | verdict |
|---|---:|---:|---:|---|---:|---|
| REG1 GGHV22 tex-1689, **nine printed** eqs | 13 | 45 | 24 | **`E := C_3^23·F_-4`**, α=(23,1) | 12 / **1** | PASS |
| REG1b same **+ the three eqs the paper drops** | 16 | 66 | 24 | *none* (correctly refused) | 15 / 1 | PASS |
| REG2 the seven `G5body` terms | 7 | 7 | 3 | n/a (LEVEL C) | — | PASS |
| REG3 negative control (`b` inside **and** outside `a·b²`) | 3 | 5 | 7 | *none* — refused, witness `b` in `f3` | 3 / **0** | PASS |
| REG3 positive control (same, minus `f3`) | 3 | 4 | 7 | `z := a·b²`, α=(1,2) | 2 / 1 | PASS |
| REG4 R9 z=1 dm4-eliminated H-system | 7 | 23 | 5 | **none** | 7 / **0** | PASS |
| REG4b pre-elimination full G-system | 8 | 20 | 3 | **none** | 8 / **0** | PASS |
| R9 **scalar**, 28 spare unknowns (`--deep`) | 28 | 3 360 | 3 | **none** | 28 / **0** | — |
| R9 **scalar**, unknowns + 9 state params (`--deep`) | 37 | 22 817 | 9 | **none** | 37 / **0** | — |

The alt bridge is not a separate row: `alt_eliminated_system.json`'s `H` block is
verified **byte-identical** to `r9_eliminated_system.json`'s, so the REG4 row
covers it.

### REG1 — the GGHV rediscovery

The lens is given the nine equations with `C_3`, `F_-4` and `λ` **unlumped**
(the transcription is re-verified inside the tool against an independent
regeneration from the series definitions of `D̃²`, `D̃³`, `D̃^{-1}`; a mismatch
raises). It proposes

```
E := C3**23*Fm4        merges C3, Fm4   alpha = (23, 1)
  occurring vars    13 ->  12
  distinct monoms   45 ->  45
  max total degree  24 ->   3
  graph edges       22 ->  21
  (H1) exclusivity  CHECKED on 45 distinct monomials: PASS
  (H2) nondegeneracy PASS   (L1) lift exactness PASS   (L2) numeric PASS
  ==> LEMMA L applies: J ∩ Q[U] = I ∩ Q[U].  LOSSLESS in both directions.
```

and the paper's printed answer transforms as

```
18*C3^23*d1*dm1^6*Fm4 + 8*C3^69*Fm4^3 + 27*d0*dm1^9
   ->  8*E**3 + 18*E*d1*dm1**6 + 27*d0*dm1**9
```

which is exactly the form in `GGHV_LINE1689.md` §3 (checked by expansion, not by
eye). **The headline number is the degree: 24 → 3.**

Two things the lens adds that the by-eye observation did not:

* **Uniqueness.** The exponent lattice has rank 12 over 13 occurring variables,
  torus dimension **1**, with null grading `{C3: −1, Fm4: 23}`. So at most one
  variable can ever be lost, and LEVEL A realises exactly that one
  polynomially. `E = F_-4·C_3^23` is not *a* quotient coordinate for this
  system; up to scaling the null grading it is **the** one.
* **Why the paper's selection of nine equations matters** — see REG1b.

### REG1b — the control that shows the refusal is discriminating

Adding back the three equations GGHV explicitly do **not** use —
`(D̃²)_{-6}`, `(D̃²)_{-8}`, and `(Q̃)_{-3}`, the last being the unique one
carrying `λ·C_3^20` — destroys the merge, and the lens says so:

```
LEVEL A (rank-1 merges): NONE.
exponent lattice: rank 15 over 16 occurring variables, torus dim 1
   null grading: {'C3': -1, 'Fm4': 23, 'lam': 20}
LEVEL B NEW degeneracies: 1
   S=['C3','Fm4','lam'] rank=2 realisable=True gens=['C3^20*lam','C3^23*Fm4']
```

This is a genuinely new small fact about GGHV's setup: the twelve-equation
system still depends on `C_3, F_-4, λ` only through the **two** coordinates
`λ·C_3^20` and `F_-4·C_3^23` (3 variables → 2, LEVEL B, sound direction only) —
but the *single*-coordinate collapse is available only after the paper's
selection of nine. `GGHV_LINE1689.md` §5 already argued that selection is forced
rather than arbitrary; this is an independent, mechanical corroboration from a
completely different direction.

### REG2 — the u-weight structure of the seven G5body terms

The lens is given the seven monomials and **no weights**, and solves
`M w = c·1`:

```
the 7 input monomials match FACE_KILL_SWEEP sec.1's list: True
LEVEL C grading space: 7 constraints, solution space dimension 2  (unique up to scale: False)
documented weight vector lies in the space: True
per-term u-weight under the documented vector: all seven = 17
```

**What it found, and the honest caveat it also found.** The documented vector
`d0=4, d1=3, d2=2, dm1=5, dm2=6, dm3=7, dm4=8` with common weight `17` is
confirmed — every one of the seven terms has weight 17, so `FACE_KILL_SWEEP.md`
§1's table is exactly right. But the seven monomials **by themselves do not pin
those weights**: the grading space is 2-dimensional, i.e. there is a
1-dimensional space of *null* gradings, generated by

```
{d0: -5, d1: -8, d2: -11, dm1: -2, dm2: 1, dm3: 4, dm4: 7}
```

so `w + t·null` is an equally valid weighting of `G5body` for every `t`. The
generator has mixed signs, so it is a Laurent symmetry only — no polynomial
reduction follows — but the point stands: *G5body alone does not determine the
u-weights.*

Feeding the lens the **full** generator set `G1, G2, G3, G5 = Φ + G5body` closes
this:

```
solution space dimension 1; documented (w, Φ=17; 13,14,15,17) in it: True
primitive generator: {d0:4, d1:3, d2:2, dm1:5, dm2:6, dm3:7, dm4:8, Phi:17,
                      c(G1):13, c(G2):14, c(G3):15, c(G5):17}
```

Dimension 1 means unique up to overall scale, and the primitive integer
generator is *exactly* the documented weighting including `Φ = 17` and the
generator weights `13/14/15/17`. So the campaign's u-grading is recovered from
the monomials alone, with no input from the derivation — but it needs all four
generators to be determined, and any writeup that derives it from `G5` alone is
under-determined. That is a scope correction, not a refutation: the numbers in
`FACE_KILL_SWEEP.md` §1 are unchanged.

### REG3 — the negative test

```
f1 = a*b^2 + c        f2 = a^2*b^4*c - 3 = (a b^2)^2 c - 3        f3 = b + c
```

`b` occurs inside `a·b²` and, in `f3`, on its own. The lens proposes nothing and
names the blocker:

```
LEVEL A (rank-1 merges): NONE.
blocked pairs (why no merge): a / b   blocked by  b  in f3
exponent lattice: rank 3 over 3 -- torus dimension 0
```

The refusal is checked three ways, all required to pass: (i) no proposal is
emitted; (ii) forcing `α = (1,2)` makes the (H1) check report a violation with
the witness `f3, monomial b`; (iii) `apply_merge()` **raises** rather than
silently proceeding. Positive control: delete `f3` and the merge `z := a·b²`
reappears with all obligations passing — so the refusal is discriminating, not a
tool that just never proposes anything.

**Lemma L is also *executed* here, not only proved.** On the positive control
the tool eliminates `{a,b}` upstairs and `{z}` downstairs with lex Gröbner bases
and compares the contractions into `Q[c]`:

```
I cap Q[c] (eliminate a,b upstairs) = ['c**3 - 3']
J cap Q[c] (eliminate z downstairs) = ['c**3 - 3']
equal up to scalars: True
```

That equality is part of the REG3 pass condition, so a future break in the lemma
or in `apply_merge` surfaces as a test failure rather than as a wrong answer.

### REG4 — THE REAL QUESTION: does R9 contain a comparable latent coordinate?

**No, and this is a proof, not a search report.**

The R9 z=1 dm4-eliminated H-system is loaded exactly the way
`bigrade_annotator._H_generators` does it (the module is imported; nothing is
re-transcribed), i.e. verbatim from `r9_eliminated_system.json`.

```
REG4  R9 z=1 dm4-eliminated H-system
  declared vars 7 (occurring 7) | distinct monomials 23 | max total deg 5
  interaction graph: 15 edges, 1 component, max degree 6
  single-site variables: Phi  only in  dm1*Phi
  LEVEL A (rank-1 merges): NONE.
  exponent lattice: rank 7 over 7 occurring variables
     torus dimension 0 -- NO weight vector makes every monomial weight 0
  LEVEL B (rank r < |S|): no NEW degeneracies (0 subsumed by LEVEL A)
  blocked pairs: d0/dm1 blocked by dm1*dm3^2 in H2
                 d0/dm2 blocked by dm2^2*dm3 in H2  ...
```

(The LEVEL B line costs nothing here: by fact 3 of §0, torus dimension 0 makes
every subset non-degenerate, so the scan returns empty without running.)

and the same for the pre-elimination full G-system (rank 8 of 8, torus
dimension 0, no merges).

**The alt bridge.** `alt_eliminated_system.json` records
`"shares_elimination_with": "r9_eliminated_system.json"`; the tool *checks* this
rather than trusting it, and confirms `alt["H"] == r9["H"]` exactly. The alt
bridge is the same symbol-level object, so one run covers both. (The alt regime
differs only in the spare **caps**, sub1 18/21/24 vs sub2 12/14/16 — a scalar
window question, not a monomial-support one, and the lens reads only support.)

**Why torus dimension 0 settles it.** If any monomial map from fewer variables
had all 23 monomials in its image, the exponent rows would lie in a proper
sublattice and some nonzero weight vector would annihilate every monomial. Rank
7 of 7 says no such vector exists. So there is no quotient coordinate for the
R9 H-system — not a rank-1 one, not a rank-`r` one, not a Laurent one, not one
with rational exponents. The search space is empty, not merely unsearched.

**The near miss, and why it is a miss.** `Φ` is a *single-site* variable: it
occurs in exactly one monomial, `dm1·Φ`, exactly as `C_3` and `F_-4` occur only
in `F_-4·C_3^23`. That is the same local configuration that made GGHV collapse.
It fails here for one reason: `dm1` also occurs *everywhere else* (`d0·dm1³`,
`d1·dm1²·dm2`, `dm1⁴`, …), so `{Φ, dm1}` violates (H1) and the lens refuses. The
GGHV case worked because **both** partners were single-site. The lens's blocked-
pair diagnosis names this directly.

The honest reading: `E = F_-4·C_3^23` was available because two variables in
GGHV's system appeared in exactly one place in the whole system — a degenerate,
almost bookkeeping-level feature. It is not a mechanism that recurs.

## 3. Scalar level — the same answer, at the level where it would matter

The symbol-level "no" could in principle be an artifact of looking at the wrong
level: what a solve actually faces is the **scalarised** system, the
y-coefficients of `H2, H3, H5` after the R9 z=1 state substitution, in the 28
spare unknowns `R0..R12` (`dm2`, cap 12) and `S0..S14` (`dm3`, cap 14). So the
lens was pointed there too. The block below is the verbatim output of
`monomial_lens.scalar_R9()` — the function `--deep` invokes — timed at
**218.6 s** and captured in `scalar_r9_lens.log`:

```
R9 SCALAR  variables = the 28 spare unknowns R0..R12,S0..S14
  declared vars 28 (occurring 28) | distinct monomials 3360 | max total deg 3
  interaction graph: 378 edges, 1 component, max degree 27
  LEVEL A (rank-1 merges): NONE.
  exponent lattice: rank 28 over 28  ->  torus dimension 0
  blocked pairs: R0 / R1  blocked by  R0^2*S0  in H2, ...

R9 SCALAR  variables = 28 spare unknowns + 9 state parameters
  declared vars 37 (occurring 37) | distinct monomials 22817 | max total deg 9
  interaction graph: 635 edges, 1 component, max degree 36
  LEVEL A (rank-1 merges): NONE.
  exponent lattice: rank 37 over 37  ->  torus dimension 0
```

**Torus dimension 0 at both scalar levels too** — with the parameters
`gamma, r, a0..a4, g0, g1` declared as variables as well as with them held as
coefficients. So the negative is not a level artifact.

**Why this does not need the expensive y-coefficient split.** Splitting a
substituted generator into its y-coefficients *partitions* terms; it neither
creates nor destroys monomials in the unknowns. Formally: view the substituted
generator as a polynomial in the unknowns with coefficients in `Q[y, params]`;
the coefficient of an unknown-monomial is a polynomial in `y` whose `j`-th
coefficient is exactly that monomial's coefficient in the `j`-th scalar
equation, so it vanishes iff the monomial is absent from *every* scalar
equation. The union of the supports of the scalar equations therefore equals the
support of the unsplit expression, and LEVEL A / LEVEL B read only that union.
This is what turns the documented multi-hour `build_R3()` sympy trap into a
218-second run. (LEVEL C would differ, since it wants one constant per equation;
it is not used at this level.)

The tool implements the argument, does not merely assert it: `scalar_R9()`
carries it in its docstring and takes the unsplit route deliberately.

## 4. Honesty section — what is and is not established

**Established.**

* The tool runs; `python monomial_lens.py --quiet` exits **0** with all five
  reported checks PASS (the four required plus the REG1b control). Everything
  quoted in §2 is copied from an actual run, not paraphrased.
* LEMMA L is proved above and its hypotheses are machine-checked per proposal,
  including an exact lift identity `φ(g_i) − f_i ≡ 0` and an independent
  exact-rational spot check.
* LEVEL A is **complete** for rank-1 merges (column-proportionality argument,
  §0), so a LEVEL A "none" is a statement about the system, not about the
  search.
* The exponent-lattice/torus-dimension test is a **complete** obstruction for
  monomial coordinate changes of any rank, and it returns 0 for the R9 H-system,
  the full G-system, both scalar levels of the R9 z=1 state, and (by
  byte-identity) the alt bridge.
* LEMMA L is not only proved but **executed** on the REG3 positive control:
  `I ∩ Q[c] = J ∩ Q[c] = (c³ − 3)`, computed by two independent lex Gröbner
  eliminations, and that equality is part of the pass condition.
* The LEVEL B search is complete-by-construction too, not a sampled scan: every
  degenerate subset contains the support of a nonzero null-space vector, so the
  scan runs over minimal supports of lattice relations rather than over
  `C(n,k)` subsets.
* The GGHV transcription used here is re-verified inside the tool against an
  independent regeneration from the series definitions; a mismatch raises rather
  than warns.

**NOT established / explicit limits.**

* **The lens sees monomial support only.** It cannot detect a non-monomial
  change of coordinates (`u = d0 + d2²`, a resultant, a Galois descent, a
  Möbius substitution in `y`). A "no" from this tool is a "no latent *monomial*
  quotient coordinate" and nothing more. R9 may well still have exploitable
  structure of a kind the lens is blind to by construction.
* **The minimal-support enumeration is bounded.** LEVEL B enumerates integer
  combinations of the null-space basis with coefficients in `[−3, 3]`. For null
  dimension 1 (every case that arose here) that is exact; for higher dimension
  it is a bounded search and could in principle miss a smaller support. It
  cannot produce a false positive, and it is irrelevant whenever the torus
  dimension is 0, which is every negative reported here.
* **LEVEL B's converse is undischarged.** For `r ≥ 2` the tool claims only
  `J ∩ Q[U] ⊆ I ∩ Q[U]`; the equality needs freeness of `B` over `φ(A)`, which
  can fail. The REG1b `{C_3, F_-4, λ} → {λC_3^20, F_-4C_3^23}` finding is
  reported at that weaker strength. It has **not** been used for anything.
* **REG4's PASS is a pin, not evidence.** The *underlying* negative (torus
  dimension 0) is a proof, but the self-check's role is only to assert the
  *recorded* answer (0 merges, 0 new LEVEL B, torus dim 0, alt identical) so
  that a future change becomes loud. Do not cite the PASS as corroboration;
  cite the rank computation.
* **The `--deep` scalar pass is not part of `--quiet`'s exit condition.** It is
  reported in §3 from a direct `scalar_R9()` run (218.6 s, `scalar_r9_lens.log`);
  a regression there would not fail the self-check. The `--deep` *flag wiring*
  was verified separately with a stub, but a full `--deep --quiet` CLI run was
  not carried to completion in this session (the box was running several other
  lanes and the pass exceeded 25 min under contention). Nothing in §3 depends on
  that: `--deep` calls exactly the function whose output §3 quotes.
* **Nothing here is a kill, a certificate, or a new mathematical result about
  (72,108).** The only new facts are: (a) the mechanised confirmation and
  uniqueness of GGHV's `E`; (b) the REG1b two-coordinate observation; (c) the
  REG2 scope correction (G5body alone under-determines the u-weights); (d) the
  proved absence of monomial quotient coordinates in R9 / G / alt.
* **No existing artifact was modified.** The tool is read-only on
  `r9_eliminated_system.json`, `alt_eliminated_system.json` and
  `bigrade_annotator.py` (imported, not edited).
* The REG2 caveat is a correction to how the u-weights are *justified*, not to
  their values. `FACE_KILL_SWEEP.md` §1's table is confirmed term by term.

**The bottom line asked for.** The line-1689 collapse was a lucky feature of
that particular system — two variables that each occurred in exactly one
monomial of the whole system — rather than a recurring mechanism. R9 contains no
comparable latent quotient coordinate, and the lens proves it cannot.
