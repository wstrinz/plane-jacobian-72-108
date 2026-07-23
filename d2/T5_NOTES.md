# T5 — exact certificate: exploration log and the real lead (2026-07-21)

> **STATUS (2026-07-22):** SUPERSEDED — this exploration log predates both the field-split repair (`FIELD_SPLIT_AUDIT.md`) and the f37 closure (`F37_SATURATION_REPORT.md`). It still treats f37 as a live alternative and proposes a by-hand (y+1)-adic count that has since been replaced by the split-place ledger + cascade engine. Authoritative sources: `STATE.md`, `PROOF_INVENTORY.md`, and the cascade artifacts (`cascade_cones_qt_inf_rl.json`).

Goal: turn the numeric positive floors into a PROOF that f31 ≡ 0 and f37 ≡ 0 are
both impossible in the proven windows (with Φ the fixed forcing polynomial). This
is the genuine open endgame — the step the 2022 authors could not do — so this
file records what was tried, what is ruled out, and where the tractable hope is.
It does NOT claim a proof.

## Setup recap (all now proven, see T3/T6 audits)

Master identity: f31·f37·d₋₁²¹ ≡ 0, d₋₁ ≢ 0 ⇒ **f31 ≡ 0 or f37 ≡ 0**. Windows
(per variable d_{4−k}, weight w): ord ≥ 12w, deg ≤ 14w (sub2) / 15w (sub1). f31
is weighted-homogeneous of weight 125 in (d2,d1,d0,d₋₁,Φ)=(2,3,4,5,17); every
monomial therefore has y-min-degree exactly 12·125 = 1500 and y-max-degree
≤ 14·125 = 1750 (sub2). So f31(d(y),Φ(y)) ≡ 0 is 251 scalar slice-equations
(y^1500..y^1750) in the 32 window coefficients. Massively overdetermined; the
numerics (T2, and subcase-1 this session) say no nonzero solution.

## Ruled out / dead ends

1. **The "universal cubic" 18Φd₁d₋₁⁶ + 8Φ³ + 27d₀d₋₁⁹ is NOT a consequence of the
   t=4 system.** Singular Gröbner reduction: it does not lie in the ideal
   ⟨G1,G2,G3,G5body+Φ⟩ (remainder = itself). STATE.md's "t=3 universal relation,
   y-free" does not transfer to our t=4 (8,28) system as an identity. The paper's
   actual "ecuación principal" analogue for our case IS the f31·f37·d₋₁²¹
   factorization, not this cubic — do not build a certificate on the cubic.
   (The cubic may still hold in some graded/leading sense; unverified.)

2. **Single-slice contradiction: no.** The bottom slice is
   f31(τ2,τ1,τ0,τ₋₁, −1/34) = 0 (τ = trailing coeffs, y^{12w}); the top slice is
   f31(λ2,λ1,λ0,λ₋₁, −1024/3315) = 0 (λ = leading coeffs, y^{14w}). Each is an
   *irreducible* degree-31 hypersurface (f37: degree-37), so each alone has a
   3-dim solution variety — no contradiction from one slice. (f31|_{d₋₁=0}
   = 8192 d₁²Φ⁷, consistent with the d₋₁≡0 branch being separate.)

3. **Raw degree/dimension counting: no** (as STATE.md item 6 already warned).
   251 eqns ≫ 32 unknowns kills a *generic* f, but f31's specific coefficients
   are what must be used — the paper's whole point is that counting is not enough
   for (8,28). This is exactly why it is the last open case.

## The real lead: strip the ord bound, then a (y+1)-adic order-30 argument

This is the paper's own closing technique (T3_WINDOW_AUDIT.md §4, endgame lines
1748–1786): use a proven divisibility to strip, then count (y+1)-multiplicity
against degree. We already HAVE the divisibility: the proven ord bound
y^{12w} | d_k. Concretely:

**Step 1 — strip (exact, uses only the proven ord bound).** Write
d_k = y^{12 w_k} d̃_k. Since every f31 monomial has y-order exactly 12·125=1500,
    f31(d(y), Φ(y)) = y^{1500} · f31(d̃(y), Φ̃(y)),   Φ̃ := Φ/y^{204},
so f31 ≡ 0 ⟺ f31(d̃, Φ̃) ≡ 0, where now:
- d̃_k has ord ≥ 0 and **degree ≤ 2w_k (sub2) / 3w_k (sub1)** — small: sub2 gives
  deg d̃2≤4, d̃1≤6, d̃0≤8, d̃₋₁≤10.
- Φ̃ = −(y+1)^{30}(2048y⁴−512y³+320y²−240y+195)/6630 — degree 34, and crucially
  **v_{(y+1)}(Φ̃) = 30** with (y+1)-unit cofactor (quartic(−1)=3315≠0).

**Step 2 — the pure-d factor (verified).** Setting Φ=0,
    f31|_{Φ=0} = d₋₁^{21} · h31,    f37|_{Φ=0} = d₋₁^{18} · h37,
with h31 irreducible, weighted-homogeneous of weight 20 (28 terms, degree 10 in
the d's), h37 weight 44 (145 terms). Both checked by factorization.

**Step 3 — the (y+1)-adic order count (the endgame).** Set t = y+1. Since
Φ̃ = O(t^{30}), the identity f31(d̃,Φ̃) ≡ 0 read modulo t^{30} sees ONLY the pure-d
part: d̃₋₁^{21} · h31(d̃) ≡ 0 (mod t^{30}). Let a = v_t(d̃₋₁) ≥ 0. Then
21a + v_t(h31(d̃)) ≥ 30, and at order t^{30} the Φ̃-term first enters, pinning the
top coefficient. The clean sub-case a = 0 (i.e. (y+1) ∤ d₋₁ beyond the strip)
forces **(y+1)^{30} | h31(d̃(y))** — but h31(d̃) is a polynomial of degree ≤ 40
(sub2: weight-20 in d̃'s of degree ≤ 2·weight), so this is 30 linear conditions
on a degree-≤40 object: extremely tight, plausibly empty once the remaining
slice-equations (orders t^{31}..) and the top-slice anchor are added. The a > 0
sub-cases push (y+1) into d₋₁ and are bounded the same way.

This is a *finite, concrete* elimination: h31 is explicit (weight 20, 28 terms),
the d̃_k have single-digit degrees, and the whole thing lives in K[t]/(t^N). It is
the exact analogue of the paper's y-adic strip + degree count, transported to the
(y+1)-adic place. It has NOT been carried out to a contradiction here — that is
the recommended next work item, and it is now fully set up (no missing valuation
input; the ord bound that powers it is already proven).

## Computational fallback: exact rational jet-lifting

If the valuation argument stalls, the numeric jetlift can be redone in exact
Q-arithmetic: fix an exact bottom-slice base, lift slices upward exactly, and
locate the first slice whose consistency equation is forced nonzero over the
remaining free dof. CAVEAT: the numeric obstruction is a *global* least-squares
floor (energy spread over slices y^11–y^40 relative to the base), NOT a single
inconsistent slice — so exact lifting must carry the free-dof symbolically and
test an overdetermined (nonlinear, degree-25) system for inconsistency, not just
evaluate one slice. Tractable in principle for a fixed base stratum; heavy.

## Attempted: direct finite-field Gröbner on the window system (INTRACTABLE)

Formulated the exact feasibility of f31·sub2 as a polynomial system: strip to
d̃_k (32 unknown coeffs, degrees 4/6/8/10), form f31(d̃(y),Φ̃(y)) over F_32003,
take the low y-slices y^0..y^41 (42 equations, overdetermined vs 32 unknowns),
and Gröbner for 1 ∈ ideal (⇒ no solution ⇒ proof). Singular did NOT terminate
in 500s and was killed. So brute-force Gröbner on the window system does not
close it — unsurprising for the last open case, and a useful negative: the
certificate needs structure, not raw elimination. The barrier is the 32 vars ×
degree-25 equations; the obvious reductions to try next:
- **Exact jet-lift first**: most top window coefficients enter the high slices
  *linearly* (the jetlift gradient structure). Solve/eliminate those linear
  variables exactly, leaving a much smaller nonlinear core for Gröbner.
- **Dehomogenize** (weighted-homogeneous ⇒ fix one leading coeff = 1) to cut a
  dimension, and use an elimination order.
- **Work the (y+1)-adic reduction** (above) to replace f31 (102 terms, deg 25)
  by h31 (28 terms, deg 10) in the low-order conditions.

## Attempted: exact jet-lift pilot (exact_jetlift.py) — hits expression swell

Built the exact jet-lift (exact_jetlift.py): fix a generic base
(a0,b0,c0,e0) over F_p on the bottom-slice variety, keep the 28 higher window
coeffs symbolic, impose the y-slice equations, and Gröbner-test inconsistency
(= no lift completes over that base = partial certificate). Findings:
- Base over F_p found fine; residual symbolic dof = 28 (base-fixed), i.e. the
  analytic count 21 total dof = 3 base + 18 free, here with base pinned.
- The current formulation forms the FULL substituted polynomial
  f31(d~(y),Phi~(y)) (via Singular `coeffs(F,y)`) before extracting slices, and
  that expansion does NOT complete even for a 12-slice run: the degree-25 anchor
  d~_{-1}^25 with symbolic coeffs is the blow-up. Direct Gröbner on the
  resulting (underdetermined, high-degree, 28-var) system is also intractable
  (>600s at 45 slices), matching the earlier 32-var timeout.
- LESSON: the certificate must be computed INCREMENTALLY — never form the full
  product. Compute slice j from the degree-<=j truncations, solve the (linear!)
  slice-j equation for one pivot coeff, substitute, advance. That is what
  jetlift.py does numerically; the exact version needs the same discipline over
  Q/F_p, plus aggressive reduction to keep the free-param expressions from
  swelling past the degree-25 anchor. Only the low consistency slices (j~11-20,
  degree <=j in the ~18 free params) are plausibly Gröbner-tractable, and that
  is where the numeric obstruction concentrates — so an incremental lift that
  reaches slices ~11-25 and Gröbner-tests just those is the concrete next build.

## Attempted: incremental exact jet-lift (incremental_lift.py)

Built the correct incremental formulation (never forms the full product; lifts
slice by slice from degree-≤j truncations; pivots the linear slice equations;
Gröbner-tests the consistency slices in the ~18 free params). This is the right
algorithm, but implemented in pure sympy it is too slow: multivariate
polynomial arithmetic mod p over 28 symbols, with the degree-25 anchor
d̃₋₁²⁵ feeding coefficient reductions, stalls before producing consistency
slices. It needs a low-level finite-field polynomial backend (FLINT / Singular
at the C level / a specialized incremental Gröbner), not sympy.

Net of the three computational attempts (raw 32-var Gröbner; base-fixed full-F
Gröbner; incremental sympy lift): the exact certificate is a well-defined FINITE
computation that resists all off-the-shelf tooling, for two compounding reasons
— the degree-25 anchor causes expression swell, and the ~18–28 variable count
with high degree defeats Gröbner. The two scripts (exact_jetlift.py,
incremental_lift.py) correctly encode the problem and are scaffolding for a
serious optimized attempt.

## PROGRESS (2026-07-21): graded structure + the w = Phi/d_-1^3 reformulation

A clean structural fact, verified exactly, that shrinks the whole problem:

  **f31 = sum_{f=0}^{7} Phi^f * d_-1^{21-3f} * h_f(d2,d1,d0,d_-1)**,

with h_f weighted-homogeneous of weight 20-2f and SMALL:
  f:        0   1   2   3   4   5   6   7
  d_-1 pow: 21  18  15  12   9   6   3   0     (drops by exactly 3)
  h_f terms:28  22  17  12  10   7   5   1     (h_7 = 8192 d1^2)
  h_f deg<=:40  36  32  28  24  20  16  12  (in y, on the stripped windows)

The uniform (21-3f) power is the key: set **w := Phi~ / d~_{-1}^3**. Then
  f31(d~,Phi~) = d~_{-1}^{21} * H(d~, w),   H(d~,W) := sum_f h_f(d~) W^f,
so (since d_-1 != 0) the certificate becomes

  **H(d~(y), w(y)) == 0,  i.e. w = Phi~/d~_{-1}^3 is a K(y)-rational ROOT of the
    degree-7 polynomial H(d~(y), W) whose 8 coefficients h_f are small.**

This replaces f31 (102 terms, deg 25, anchor d_-1^25) by 8 small h_f and a
degree-7 root condition — a much better object for BOTH remaining paths:
- **By-hand:** w has (y+1)-adic valuation v(w) = 30 - 3a, a := v_{y+1}(d~_{-1}).
  So 30-3a must be a slope of the (y+1)-Newton polygon of H, i.e. there are
  f1<f2 with beta_{f1}-beta_{f2} = (30-3a)(f2-f1), beta_f := v_{y+1}(h_f(d~)).
  For the generic case a=0 this forces v(h31(d~)) in [30,40] AND
  v(h_1(d~)) = v(h31(d~)) - 30 — a sharp coupled vanishing condition on two
  SMALL polynomials, the natural next thing to push (and possibly refute) by
  hand or by a small Groebner. Not yet carried to a contradiction.
- **Optimized jet-lift (for the eng session):** lift on H(d~,w) instead of f31.
  Compute w = Phi~/d~_{-1}^3 as a truncated series (d~_{-1}(0)=e0 != 0, so it
  inverts), then the slice equations use the small h_f (deg <=10) and W-powers
  up to 7 — NO degree-25 anchor. This is the reduction that was missing; it is
  the recommended starting point for the optimized implementation.

f37 has an ANALOGOUS but non-uniform grading (d_-1 powers 18,15,12,9,6,4,2,0 —
the last three break the -3 pattern), so the clean w-substitution is specific to
f31 (the genuine generic branch; numeric solutions all sit on f31=0). f37 needs
its own handling.

## Honest status of T5

Not closed. The exact certificate is the genuine open core — the step the 2022
authors could not do — and this session establishes concretely that it does not
fall to off-the-shelf computation (three independent formulations, all
intractable/too-slow). Real progress: the exact problem is set up precisely;
the cubic-relation and single-slice routes are ruled out; the (y+1)-adic
strip-and-count is set up with verified structural facts (f31|_{Φ=0}=d₋₁²¹h31,
f37|_{Φ=0}=d₋₁¹⁸h37); and the computational difficulty is now mapped (swell +
variable count). Two realistic paths remain, and the *by-hand one now looks
more promising than the computational one*:
  1. **(y+1)-adic order/degree count by hand** on the explicit small objects
     h31 (28 terms, deg 10) / h37 — the paper's own proven technique, working
     with objects orders of magnitude smaller than f31. Recommended.
  2. **Optimized incremental jet-lift** with a C-level finite-field backend —
     the scaffolding is here, but it is real engineering with no guarantee the
     final elimination is tractable.
Numeric evidence (all four branches, positive floors, proven windows) remains
strong but is not a proof.
