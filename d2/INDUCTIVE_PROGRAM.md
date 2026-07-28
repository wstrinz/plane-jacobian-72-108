# INDUCTIVE_PROGRAM.md — the grand-goal north star

**Date:** 2026-07-23. **Status:** PROGRAM SKETCH, not claims. The question
this document tracks: is there an inductive/parametric structure over the
GGV case chain strong enough to attack JC(2) itself, rather than one case?

## Why this is the only endgame that matters

The GGV case list grows with degree (24 candidates in [125,150] alone;
corners recur). Case-by-case closure raises bounds forever and proves JC(2)
never. Every technique investment should be scored against: does it
transfer parametrically?

## The conjectured three-layer parametric shape (evidence-based)

1. TROPICAL/VALUATION LAYER — parametric already in practice: the engines
   run window-generically; kill laws are affine in a (t-depth = 30-3a);
   the sub1 26-family was a-independent until infinity broke it at low a
   in an understood way. STATUS: strong evidence FOR transfer.
2. RESIDUE/ARITHMETIC LAYER — uniform in SHAPE, per-case in FIELD: the
   kill arithmetic (square classes vs the splitting field; 17-divides-5
   divisibility; r=1/4 irrationality) depends on the case's forcing
   polynomial (disc-17 quartic here; 10th cyclotomic at (108,144)).
   NEEDED: a Galois-uniform obstruction lemma. STATUS: open; corner-144
   is one data point against naive transfer.
3. TERMINAL/GLOBAL LAYER — the divisor-lemma engine is the candidate
   transfer vehicle: rank/dimension statements over confluent Vandermonde
   matrices with resultant-product determinants are POLYNOMIAL IN THE
   CASE DATA, hence parametrizable if the nonvanishing has a structural
   reason. Same for the pre-resultant question: is one resultant factor
   ALWAYS the elimination generator (f37-style excess universal)? A
   Bezout/Koszul-structured proof of the membership certificate would
   answer this in family.

## The test instances (in order of cost)

- T1: run the full pipeline on (75,125) generated FROM POLYGON DATA (no
  hand rederivation) — the reviewer's scalability test; also now
  cross-checkable against Helali's independent (75,125) systems if the
  exchange yields his variable map.
- T2: the (108,144) corner recurrence — does the same 4-generator
  pre-resultant ideal appear with changed weights? (corner-144 partial
  data exists: skeleton yes, numerics no.)
- T3: the D3-shadow experiment — push the Alpoge-Fable mechanism (affine
  in the third variable) through the plane cascade and isolate WHICH
  plane lemma kills it (Tao's scaling-exponent observation, n-2 = 0 at
  n=2, is the analytic form; we can get the algebraic form).

## Cross-program structure as evidence

Two-plus independent formulations (ours: eliminate-then-places; Helali:
direct-depth-then-triangularize; Santibanez: forced-edge covectors) that
must describe the same geometry are themselves an inductive resource:
recurring structures visible in ALL formulations (the F^2|G squeeze <->
R2|M^2; small-branch collapse to Q + one quadratic field) are the
candidates for case-independent theorems. The variable-map exchange is
therefore not diplomacy — it is data collection for layer 2.

## Anti-goals (keep honest)

- No claim that the three layers suffice; layer 2's field-dependence may
  be irreducible (in which case JC(2) needs genuinely new mathematics and
  this program's ceiling is bounds + machinery).
- Bound-pushing remains the falsifiable near-term output; the north star
  guides INVESTMENT, not claims.


## ADOPTED FROM THE THIRD EXTERNAL REVIEW (2026-07-23)

- INDUCTION VARIABLE: not the degree bound (that is enumeration order) â€”
  the real variables are chain combinatorial type, reduced power pair,
  selected multiplicities/residual divisor, and the family parameter j.
  Minimal transferable signature Sigma = (power pair, chain type,
  multiplicities, div(c), div(Phi), window functions, G-system type).
- KEY SENTENCE (verbatim): "After 108, do not attack 125 as one more
  case. Attack the F_2 family, using 125 as the first nontrivial
  specialization." Experiment A = F_2(j) pipeline with j symbolic,
  j=0 as regression control; Experiment B = the 144 corner recurrence
  (different axis: corner fixed, pair/multiplicity varying).
- NEW IDEA 1 â€” MONOTONE POTENTIAL: promote the defect vector to a
  potential on complete chains (does each transition consume section
  freedom predictably? -> long chains die uniformly). CHEAP FIRST TEST
  with existing data: defect profiles across the 108 vs 144 instances.
- NEW IDEA 2 â€” CONFIGURATION-SPACE GALOIS DESCENT: universal marked
  roots, discriminant inverted, residues as root-difference products,
  norms/elimination -> invariant obstruction polynomial factoring into
  discriminants/resultants. Uniformity of METHOD and invariant shape,
  not one fixed arithmetic contradiction. (Resolves this document's
  layer-2 open question into a program.)
- ARCHITECTURE: the CASE COMPILER (input: polygon/chain data; output:
  normalized powers, forcing ODE + c(y), Phi + divisor, D_k + windows,
  the pre-resultant G-system, place-split cascade signature, defect/
  reconstruction matrices, proof-producing obligations). AMENDMENT
  (ours): emit BOTH presentations â€” eliminated (f31-style) for the
  tropical/discovery layers, pre-resultant for terminal decisions;
  the repo's history shows each layer wants its own presentation
  ("no early elimination" overcorrects).
- ALLOCATION until 108 closes: ~70% audited full-system sweep +
  roll-up; ~25% compiler + Experiments A/B; ~5% D3-shadow.
- CALIBRATION (review's, endorsed): reusable compiler + family-level
  architecture HIGH; killing one infinite GGV family MODERATE-HIGH;
  a broader (8,28) theorem MODERATE; the three layers sufficing for
  all of JC(2) without a new global invariant LOW-MODERATE. The
  likely endpoint: a finite grammar of universal lemmas, not one
  formula.
