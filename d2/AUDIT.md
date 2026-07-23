# D2 audit — is this a stable base? (2026-07-20)

Pre-merge audit of the (72,108) work. Goal: separate what is **rigorously
verified** from what is **still our own unchecked derivation**, so the next
session knows exactly where the ground is solid and where to tread carefully.

Verdict up front (UPDATED 2026-07-21, after the T3/T6 audits closed):
**the reproduction, the algebraic elimination core, the polygon transcription,
the envelope/window bounds, AND the equation-selection derivation are now all
verified.** The only residue is a few deepest premises (GGV1 Props 1.13/2.1,
the α-strip WLOG, the C₄ normalization) that are the paper's own trusted
machinery — see T6_SELECTION_AUDIT.md §4. The subcase-2 infeasibility evidence
now rests on proven windows; subcase-1 windows are proven valid too, so those
numerics are meaningful. What remains for a THEOREM is not audit but T5: an
exact certificate turning the numeric positive floors into a proof.

See T3_WINDOW_AUDIT.md (windows), T6_SELECTION_AUDIT.md + verify_derivation.py
(derivation, 48 checks), and PAPER_NOTES.md (independent paper cross-check).
The per-risk status below is kept for history; resolutions are marked inline.

## A. Verified rigorously this session (exact arithmetic)

1. **Reproduction (T1).** From a clean container: t=3 generator validation
   passes, and Singular re-derives f31 (102 terms) and f37 (618 terms)
   **byte-identical** to the enclosed factors (coeff ratio exactly 1). The
   `[0]`-selection of the dm2-factor in `regenerate_system.py` is safe: A and B
   each have **exactly one** factor containing dm2 (degrees 10, 12).

2. **Master-identity logic is sound.** A = dm1 · Ah, B = dm1 · Bh (the debris
   is *exactly* dm1, verified). At any solution both dm3-resultants vanish;
   since the only surviving branch has dm1 ≠ 0 (see 3), Ah = Bh = 0, hence
   R = Res_dm2(Ah,Bh) = 0, and R = unit · f31 · f37 · dm1^21 (Singular). So a
   solution forces f31 ≡ 0 or f37 ≡ 0. The resultant chain is a genuine
   necessary condition — no gap.

3. **d₋₁ ≡ 0 is impossible (re-verified symbolically).** On the reduced
   equations: G1|_{dm1=0} = 3·dm2·dm3 ⇒ dm2≡0 or dm3≡0.
   - dm2≡0: G2 → (3/2)dm3² ⇒ dm3≡0; then (G5body+Φ) → Φ ⇒ Φ≡0.
   - dm3≡0 (dm2≢0): G3 → −(3/2)d1·dm2² ⇒ d1≡0; then (G5body+Φ) → Φ − 3d1·dm2·dm4
     ⇒ Φ≡0.
   Both contradict Φ = f₁·C₄²⁸ ≠ 0. Clean, denominator-free.

4. **ODE forcing term.** f₁ = −y⁸(y+1)²(2048y⁴−512y³+320y²−240y+195)/6630 satisfies
   8y(y+1)f₁′ − 14(8y+7)f₁ = y⁸(y+1)² exactly; the quartic is squarefree and
   is not divisible by y or (y+1). Matches STATE.md item 2.

5. **Φ and factor structure.** Φ = f₁·C₄²⁸ has deg 238, ord 204, trailing coeff
   −1/34, leading −1024/3315, (y+1)-mult 30 — all exact. f31 is
   weighted-homogeneous of weight 125, f37 of weight 134 (single weight each,
   no stray terms), under w(d2,d1,d0,dm1,Φ)=(2,3,4,5,17).

6. **Harness validity.** Positive controls (fewer conditions than dof) reach
   ~1e-7 for both f31 and f37 ⇒ the optimizer finds solutions when they exist.

## B. Corroborated against the published paper (arXiv:2204.14178)

The paper (Guccione–Guccione–Horruitiner–Valqui, "…from 100 to 108") confirms
the **framing**, nearly verbatim:
- "the only exception is the case (deg P, deg Q) = (72,108) … if one manages to
  discard this case, it would increase the lower bound from 108 up to 125."
- "for the other case with (deg(P),deg(Q))=(72,108) we couldn't solve the
  corresponding system of polynomial equations, thus it is left open."
- Worked cases in the paper are (9,27) [Prop 4.1], (9,24) [Prop 4.2], and
  (7,21) [§6] — the last is exactly what our generator is validated against.

So the mission is real and correctly stated. The paper's LaTeX source was later
obtained (arxiv.org/e-print/2204.14178) and every item in C below checked
directly against it — see PAPER_NOTES.md and T3_WINDOW_AUDIT.md.

## C. Derivation audit — RESOLVED (was "the real audit debt")

Originally these were checked only against our own transcription; the T3/T6
audits (2026-07-20/21) closed them against the actual source.

- ~~Newton-polygon subcases / Prop transcription~~ **VERIFIED**: exact,
  corner-for-corner, both subcases, N(P) and N(Q); "(8,28)" is the paper's own
  label (T3_WINDOW_AUDIT.md §1, PAPER_NOTES.md).
- ~~Envelope bounds deg ≤ 14w/15w, ord ≥ 12w~~ **PROVEN** by redoing the paper's
  valuation induction on our polygons (sub2 via v_{-2,1}, sub1 via v_{-1,1}, ord
  via v_{2,-1}); Φ attains deg=14·17, ord=12·17 as the tightness witness
  (T3_WINDOW_AUDIT.md §3).
- ~~Equation selection (k=8 vacuity, λ-isolation, (D̃⁻¹)₋₅=0)~~ **VERIFIED**:
  the 12 used equations are genuine consequences and the dropped ones only
  define spare unknowns/λ (T6_SELECTION_AUDIT.md, verify_derivation.py §B,E).
- ~~[P,Q]=x² normalization chain~~ **VERIFIED down to GGV1 premises**: the
  F-normalization, forcing ODE, and slice bridge are re-derived exactly
  (verify_derivation.py §A,C,D). Only GGV1 Props 1.13/2.1, the α-strip WLOG,
  and C₄=y⁷(y+1) remain outline-level (the paper's trusted machinery;
  T6_SELECTION_AUDIT.md §4).

## D. Internal-consistency checks (necessary, not sufficient)

- jetlift CONFIGS windows match the stated bounds exactly: sub2 sizes = 2w+1
  (deg≤14w), sub1 sizes = 3w+1 (deg≤15w), hi = (degmul−12)·W+1. f31 y-window
  is [1500,1750] (sub2) / [1500,1875] (sub1). Consistent — but this only checks
  that the code matches the claimed bound, not that the bound is correct.

## Risks (ranked) — all three RESOLVED

1. ~~Subcase-1 window bound (deg ≤ 15w) may be too small~~ **RESOLVED**: 15w
   proven a genuine upper bound (T3_WINDOW_AUDIT.md §3). Subcase-1 numerics are
   now trustworthy.
2. ~~Envelope bounds unproven for subcase-2~~ **RESOLVED**: deg≤14w, ord≥12w
   proven (same induction); Φ attains both, so they're tight.
3. ~~Whole reduction unaudited~~ **RESOLVED down to GGV1 premises**: the (D̃)
   system equals the code exactly and the equation selection is sound
   (verify_derivation.py, T6_SELECTION_AUDIT.md). Residue = GGV1's own
   propositions, trusted not re-proven.

## Recommendation (UPDATED 2026-07-21)

The audit is effectively complete: the base is a **verified reproduction + a
sound elimination engine + numerics on proven windows**. It is NOT yet "the
case is closed" — that still needs the exact certificate (T5). Priority order
from here: (1) run subcase-1 numerics for f31 and f37 (windows now valid) to
complete the evidence across all four factor×subcase branches; (2) T5 exact
certificate — the slice-anchored / (y+1)-multiplicity + degree-count route the
paper uses for its own closed cases (T3_WINDOW_AUDIT.md §4, T6 §3 ammunition);
(3) T4 base-region completeness to backstop the numerics. Only a T5 certificate
converts the positive floors into the theorem.
