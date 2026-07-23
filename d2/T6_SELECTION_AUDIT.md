# T6 equation-selection audit — VERIFIED (2026-07-20)

Audit of AUDIT.md Risk 3, the last derivation-level risk: is the 12-equation
system that regenerate_system.py feeds to Singular the *right* necessary
condition for a (72,108) counterexample? Companion artifacts:

- `verify_derivation.py` — 48 symbolic checks, all passing (run it).
- `T3_WINDOW_AUDIT.md` — polygon transcription + window bounds (Risks 1–2).

Verdict: **the equation selection, λ-isolation, clearing exponents, and the
exact bridge to the code are now verified.** Two premises remain outline-only
(section 4) — they are the paper's standard machinery (GGV1 Props 1.13/2.1)
applied to our case, not our own invention, but a line-by-line check would
require auditing GGV1 itself.

## 1. The derivation chain, as verified

Setup (STATE.md items 1–2, now independently re-derived): C ∈
K[y,C₄⁻¹]((x⁻¹)), C² = P, ℓ_{1,0}(C) = x⁴C₄, C₄ = y⁷(y+1);
Q = C³ + λC⁻¹ + F.

- **v_{1,0}(F) = −5.** Since powers of C commute under the Poisson bracket,
  [P,Q] = [C², C³ + λC⁻¹ + F] = [C², F]. So x² = [P,F], and the generic
  valuation identity v([P,F]) = v(P) + v(F) − 1 gives 2 = 8 + v(F) − 1,
  v(F) = −5. (Equality case: if [ℓP, ℓF] = 0 then ℓF would be a C-power in
  the stripped range — excluded by the α-strip, see §4.)
- **The forcing ODE** (verify_derivation.py §A): Q² − P³ − 2λP = 2C³F +
  (lower), and [P, Q²] = 2Q[P,Q] = 2Qx². Taking ℓ_{1,0} of both sides with
  f₁ := C₄³F₋₅ yields the bracket identity that reduces *exactly* to
  8y(y+1)f₁′ − 14(8y+7)f₁ = y⁸(y+1)², whose unique polynomial solution is
  STATE.md's f₁ (unique against a general degree-15 ansatz — full linear
  solve, no divisibility hand-waving). The quartic factor is separable and
  prime to y(y+1). Φ = f₁C₄²⁸ has deg 238 / ord 204 / (y+1)-mult 30 /
  trailing −1/34 / leading −1024/3315, and jetlift's hard-coded Psi is
  exactly Φ/y²⁰⁴. This is the same commutator route as the paper's f₁-ODE
  for their closed case (source lines 1555–1596), adapted t=3 → t=4.
- **λ-isolation** (§B): (C⁻¹)₋₄ = C₄⁻¹ (a unit) and (C⁻¹)₋₅ = −c₃/C₄²,
  which vanishes after the d₃-killing shift x → x − D₃/4. Hence λ enters
  the Q-side polynomiality slices ONLY at j = 4 — the dropped equation —
  and the j = 5 slice is purely (C³)₋₅ + F₋₅ = 0.
- **Polynomiality of the d's** (§C): D_k := c_k·C₄^(7−2k) obeys
  D_k = ½P_{k+4}C₄^{6−2k} − ½Σ_{i+j=k+4} D_iD_j with all C₄-exponents
  cancelling (verified for k = 3 down to −13), so D_k ∈ K[y] by downward
  induction — the mirror of the paper's D_k-polynomiality proposition.
- **The slice bridge** (§D): under d_k = c_kC₄^{7−2k}, the cleared slices
  (C²)₋ₖ·C₄^{14+2k} (k = 1..9) and (C³)₋ⱼ·C₄^{21+2j} (j = 1..5) equal
  regenerate_system.py's D2(k) and D3(j) *exactly* — the abstract derivation
  and the code compute the same system. The F-term clearing at j=5 is
  F₋₅C₄³¹ = f₁C₄²⁸ = Φ, so the used equation (D̃³)₋₅ + Φ = 0 is exact.

## 2. Why the dropped/truncated equations are sound to drop

The infeasibility argument needs only that a genuine counterexample IMPLIES
the used system has a K[y]-solution in the window. Dropping equations can
only weaken a necessary condition, so the only audit question is whether the
12 *used* equations are genuine consequences — verified above. For
completeness, what is dropped and why it costs nothing generically:

- **(D̃²)₋₈ = 0**: contains dm12 linearly (coeff 2); dm12 appears in NO used
  equation (asserted in code and re-verified). The equation just defines
  dm12.
- **(D̃³)₋₄ + λC₄²⁸ = 0**: the only equation containing λ (by λ-isolation),
  and also contains dm12 (coeff 3). Given everything else it defines λ as a
  rational expression.
- **Slices deeper than the used range**: (D̃²)₋ₖ for k ≥ 10 and (D̃³)₋ⱼ for
  j ≥ 6 each introduce a FRESH unknown (dm14, dm15, …) linearly with unit
  coefficient (2 resp. 3) — always satisfiable definitions, never
  constraints on the window variables.

## 3. Unused ammunition (relevant to T5 exact-certificate work)

Two genuine constraints of a real counterexample are NOT encoded in the used
system and are available if the certificate needs more:

1. **λ-constancy.** (D̃³)₋₄ + λC₄²⁸ = 0 defines λ, but a real solution has
   λ ∈ K a *constant*: i.e. C₄²⁸ must divide (D̃³)₋₄ in K[y] with constant
   quotient. That is a strong divisibility condition on the d's, unused.
2. **dm12 polynomiality** is automatic, but the pair of dropped equations
   overdetermines dm12 (linear in both): compatibility of (D̃²)₋₈ and
   (D̃³)₋₄ is a further unused relation among the window variables and λ.

(The paper's own endgame for their closed case uses exactly this kind of
extra structure — multiplicity of (y+1) and degree counting — see
T3_WINDOW_AUDIT.md §4.)

## 4. Remaining outline-only premises (the residue of Risk 3)

- **ℓ_{1,0}(P) = R², ℓ_{1,0}(Q) = R³** with R = x⁴C₄: from GGV1 Props 1.13 +
  2.1 in the (1,0)-direction; the required inequality
  v_{1,0}([P,Q]) = 2 < v(P) + v(Q) − 1 = 19 holds, so the argument applies —
  but GGV1's propositions themselves are trusted, not re-proven.
- **The α-strip and WLOG**: ℓ(Q − C³) = α_k(x⁴C₄)^k for k ∈ (−2,3) and the
  Remark-style replacement (Q → Q − α₂P − α₀, P → P + ⅔α₁) preserving
  [P,Q] = x² and the polygons, leaving Q = C³ + λC⁻¹ + F. Same status:
  paper-template argument, verified in outline (their source lines
  1508–1546), not line-by-line.
- **C₄ = y⁷(y+1) normalization**: corner-ness of (8,14)/(8,16) forces
  C₄ = y⁷(a₀ + a₁y), a₀a₁ ≠ 0; the linear-change normalization to y⁷(y+1)
  is the same freedom the paper invokes for C₃ = y⁸(y+1).

These sit one level below everything this repo computes: they concern
whether the C-series ansatz itself is forced. Everything from the ansatz
down to the Singular input is now verified.
