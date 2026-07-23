# T3 window-bound audit — RESOLVED: sub1 deg ≤ 15w is a genuine upper bound (2026-07-20)

Audit of AUDIT.md Risks 1 and 2 (the envelope/window bounds), done against the
actual LaTeX source of arXiv:2204.14178 (fetched via `arxiv.org/e-print`, single
file `Increasing_Lower_Bound_20_04_2022.tex`, 2262 lines — no WebFetch
truncation). Verdict up front:

- **Prop 4.3 transcription: CONFIRMED verbatim** (both subcases, corners,
  [P,Q] = x² normalization). AUDIT.md item C.1 retired.
- **Envelope bounds PROVEN for both subcases** by redoing the paper's own
  induction template on our polygons: sub2 deg ≤ 14k / sub1 deg ≤ 15k /
  ord ≥ 12k per variable d_{4−k}. AUDIT.md Risks 1 and 2 retired.
- Subcase-1 numerics are therefore answering the right window question.
  (The reduction itself — equation selection, normalization chain — remains
  the open T6 debt; that is Risk 3, untouched here.)

## 1. Prop 4.3 statement, from the source (lines 1000–1007)

> Proposition [Case (8,28)]. If there is a counterexample to the Jacobian
> Conjecture in the case (8,28), then there exist P,Q ∈ L⁽¹⁾ with [P,Q] = x²
> and one of the following cases holds:
> (1) N(P) = {(0,0),(1,0),(8,14),(8,16),(0,8)},
>     N(Q) = {(0,0),(2,1),(12,21),(12,24),(0,12)}.
> (2) N(P) = {(0,0),(1,0),(8,14),(8,16)},
>     N(Q) = {(0,0),(2,1),(12,21),(12,24)}.

Identical to our transcription, including the subcase numbering and the sub1
extra corners (0,8)/(0,12). The proof (lines 1008–1311) derives these from
corners {(0,0),(1,0),(8,28),(0,4)}·(m,n), (m,n)=(3,2), via the inversion
morphism φ(x)=x⁻¹, φ(y)=x⁴y with [φP,φQ] = −[P,Q]x²; its cases a)/b) land in
subcase (2) and its case c) in subcase (1).

Framing also confirmed: the intro table (lines 304–320) lists two (72,108)
cases — A₀=(9,27),(m,n)=(2,3), closed by the paper's "teorema impossible"
section, and A₀=(8,28),(m,n)=(3,2), marked `*`, the one left open. Both
Newton subcases of (8,28) are open. Our mission statement is correct.

## 2. The paper's template (their closed (9,27) case, source lines 1416–1526)

Prop "calculo de C": C = x³C₃ + x²C₂ + …, C² = P, C₃ = y⁸(y+1), with the
recursion C_{3−k} = −(1/2C₃)(P_{6−k} + Σ_{j=1}^{k−1} C_{3−j}C_{3−k+j}), and two
valuation inductions driven purely by the Newton polygon of P:
v_{−1,1}(C) = 6 (deg side, in K((y⁻¹))) and v_{3,−1}(C) = 1 (ord side, in
K((y))). Then D_k := C_k·C₃^{5−2k} ∈ K[y], and the "magic" directions
v_{17,1}(D) = 51 and v_{−13,−1}(D) = −39 — the weights 17/13 chosen so the
bound on v(D_k x^k) is k-independent — give deg(d_{3−k}) ≤ 17k and
ord(d_{3−k}) ≥ 13k after the d₂-killing shift φ(x) = x − D₂ (which preserves
both valuations because D₂ obeys its own k=1 bounds).

## 3. The same induction on our polygons (the actual audit)

Setup: C = x⁴C₄ + x³C₃′ + …, C² = P, C₄ = y⁷(y+1) (ord 7 / deg 8 forced by
corners (8,14)/(8,16); coefficients nonzero by corner-ness; normalized by the
same linear-change freedom the paper uses). Recursion
C_{4−k} = −(1/2C₄)(P_{8−k} + Σ_{j=1}^{k−1} C_{4−j}C_{4−k+j}).

Polygon inputs (max of v over N(P) corners; whole support obeys them):

| direction | sub1 corners max | sub2 corners max | role |
|---|---|---|---|
| v_{−1,1} | 8 (at (8,16) AND (0,8)) | 8 | sub1 deg side |
| v_{−2,1} | 8 (at (0,8)) — useless | **0** (at (0,0),(8,16)) | sub2 deg side |
| v_{2,−1} | 2 (edge (1,0)–(8,14)) | 2 | ord side, both |

Inductions (all three close; verified mechanically with sympy, and the
products bound always equals the P-slice bound — the closure is exact, no
slack lost):

- **sub1 deg** (K((y⁻¹)), v_{−1,1}(C₄⁻¹) = −8): hypothesis
  v_{−1,1}(C_{4−j}) ≤ 8−j ⇒ products ≤ 16−k; support v_{−1,1}(P) ≤ 8 ⇒
  v(P_{8−k}) ≤ 16−k; step gives −8 + (16−k) = 8−k. Closes.
  ⇒ deg(C_k) ≤ k+4 (k = x-exponent).
- **sub2 deg** (v_{−2,1}, same C₄⁻¹ contribution −8): hypothesis ≤ 8−2j,
  support v_{−2,1}(P) ≤ 0 ⇒ v(P_{8−k}) ≤ 16−2k; step −8+16−2k = 8−2k. Closes.
  ⇒ deg(C_k) ≤ 2k.
- **ord, both** (K((y)), v_{2,−1}(C₄⁻¹) = +7): hypothesis ≤ 2j−7, support
  v_{2,−1}(P) ≤ 2 ⇒ v(P_{8−k}) ≤ 2k−14; step 7+2k−14 = 2k−7. Closes.
  ⇒ ord(C_k) ≥ 2k−1.

D-transform: D_k := C_k·C₄^{7−2k} (D₄ = 1). Then, k-independently:

- sub1: v_{15,1}(D_k x^k) ≤ 15k + (k+4) + 8(7−2k) = **60** ⇒ deg(d_{4−k}) ≤ 15k
- sub2: v_{14,1}(D_k x^k) ≤ 14k + 2k + 8(7−2k) = **56** ⇒ deg(d_{4−k}) ≤ 14k
- ord:  v_{−12,−1}(D_k x^k) ≤ −12k − (2k−1) − 7(7−2k) = **−48** ⇒ ord(d_{4−k}) ≥ 12k

The d₃-killing shift x → x − D₃ preserves all three (D₃ obeys its own k=1
bounds), exactly as in the paper's φ step.

**Conclusion: the jetlift windows are correct.** Per variable d_{4−k},
k = 2..5: sub2 window [12k,14k] (sizes 2k+1 = 5,7,9,11), sub1 window [12k,15k]
(sizes 3k+1 = 7,10,13,16); identity slices hi = (degmul−12)·W+1 = 251/269/
376/403 for f31/f37 × sub2/sub1. All match CONFIGS in jetlift.py exactly.
Tightness witness: Φ (weight 17) has deg 238 = 14·17 and ord 204 = 12·17 —
the sub2 bounds are attained, consistent with "empirically tight."

## 4. Bonus corroborations found in the source

- The paper's "ecuacion principal" (line 1692) for their closed case is
  18·C₃²³·d₁·d₋₁⁶·F₋₄ + 8·C₃⁶⁹·F₋₄³ + 27·d₀·d₋₁⁹ = 0 — the exact structural
  analogue of our t=3 relation 18Φd₁d₋₁⁶ + 8Φ³ + 27d₀d₋₁⁹. STATE.md called
  ours "unstated in their papers"; it IS stated for their case. Our generator
  reproduces their machinery.
- Their f₁-ODE derivation (y⁹(y+1)² = 6y(y+1)f₁′ − 10(9y+8)f₁, lines
  1571–1596) is the same commutator route as our verified
  8y(y+1)f₁′ − 14(8y+7)f₁ = y⁸(y+1)².
- Their endgame (tilde-d divisibility + degree count, lines 1764–1786) is the
  T5(a) template: exactly the multiplicity/degree interplay STATE.md item T5
  proposes for converting our floors into a proof.

## 5. What this does NOT close (remaining premises, all under T6 / Risk 3)

- Existence of R with ℓ_{1,0}(P) = R², ℓ_{1,0}(Q) = R³ for [P,Q] = x²
  (GGV1 Props 1.13 + 2.1 applied in the (1,0) direction; the inequality
  v_{1,0}([P,Q]) = 2 < 8+12−1 holds, so the standard argument applies — but
  this chain is verified only in outline here, not line-by-line).
- The equation-selection argument (k=8 vacuity, λ-isolation, (D̃⁻¹)₋₅ = 0)
  and which (D̃²)/(D̃³) slices feed the master system.
- The elimination itself is already verified exactly (AUDIT.md section A).

Sources: scratchpad copies of 2204.14178 (PDF + .tex) fetched 2026-07-20;
line numbers refer to the .tex.
