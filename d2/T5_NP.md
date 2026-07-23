# T5 — the (y+1)-adic Newton-polygon case tree for f31 (2026-07-21)

Continues T5_NOTES.md ("the real lead"). All structural inputs verified exactly
by `verify_graded.py` (run it; 5 check groups). Everything below is for the
**f31 branch, subcase (2)** stripped windows unless stated:
deg d̃2 ≤ 4, d̃1 ≤ 6, d̃0 ≤ 8, d̃₋₁ ≤ 10;  Φ̃ = t³⁰·u where t := y+1 and
u := −(2048y⁴−512y³+320y²−240y+195)/6630, u(−1) = −3315/6630 = −1/2 ≠ 0.
Graded pieces (f31_graded.txt): f31 = Σ_{f=0}^{7} Φ^f d₋₁^{21−3f} h_f,
h_f weighted-homogeneous of weight 20−2f, deg_y h_f(d̃) ≤ 40−4f on the window.
Write e := d̃₋₁, β_f := v_t(h_f(d̃)) ∈ [0, 40−4f] ∪ {∞}, a := v_t(e) ∈ [0,10].

Standing hypothesis for contradiction: a window tuple d̃ with e ≢ 0 satisfies
f31(d̃, Φ̃) ≡ 0.

## Lemma 1 (no total degeneration; NEW this session)

**If e ≢ 0 then H(d̃,W) := Σ_f h_f(d̃) W^f ≢ 0.** In fact h_5, h_6, h_7 cannot
all vanish when e ≢ 0:
h_7 = 8192 d1² ≡ 0 forces d̃1 ≡ 0; then h_6|_{d1=0} = −3072(4d0−d2²)² ≡ 0
forces d̃0 = d̃2²/4 (K[y] is a domain, char 0); then
h_5|_{d1=0, d0=d2²/4} = 2048 dm1², so h_5 ≡ 0 forces e ≡ 0. ∎
(All three identities are verify_graded.py check 4.)

Consequence: w := Φ̃/e³ ∈ K(y)^× is a **root of the nonzero degree-≤7
polynomial H(d̃,·)**, with v_t(w) = 30 − 3a.

## Lemma 2 (Newton polygon)

If H(θ) = 0 for θ ∈ K(y)^× with v := v_t(θ), then min_f (β_f + v·f) is
attained for at least two indices f. (Standard: the t-adic Newton polygon of H
must have a side of slope −v; equivalently, if the min were attained once, that
single term would dominate t-adically and H(θ) ≠ 0.)

## Theorem (a = 0 case tree: exactly three branches)

Suppose a = 0, so v = v_t(w) = 30. Pairs f₁ < f₂ attaining the min require
β_{f₁} − β_{f₂} = 30(f₂ − f₁), with 0 ≤ β_f ≤ 40 − 4f whenever h_f(d̃) ≢ 0.

* **Non-consecutive pairs are impossible:** f₂ − f₁ ≥ 2 needs β_{f₁} ≥ 60 > 40.
* **Pairs (f, f+1) with f ≥ 3 are impossible:** β_f ≥ 30 > 28 ≥ deg h_3(d̃)
  (and lower for f > 3).
* **The all-degenerate escape is closed by Lemma 1:** if h_0 ≡ h_1 ≡ h_2 ≡ 0
  then any surviving pair would need f ≥ 3 (impossible), so H(d̃,·) would have
  no valuation-30 root unless h_3 ≡ … ≡ h_6 ≡ 0 as well and H = 8192 d̃1² W⁷;
  then w ≠ 0 root forces d̃1 ≡ 0, and Lemma 1's cascade forces e ≡ 0. ⨯

So exactly one of:

**(i) generic:** h_0(d̃) ≢ 0, pair (0,1): β_0 = β_1 + 30 ∈ [30, 40].
    In particular **t³⁰ | h_0(d̃)** (a degree-≤40 object vanishing to order 30)
    and h_1(d̃) ≢ 0 with β_1 = β_0 − 30 ≤ 10.
**(ii)** h_0(d̃) ≡ 0 (identically!), h_1 ≢ 0, pair (1,2): β_1 = β_2 + 30 ∈
    [30, 36], h_2 ≢ 0, β_2 ≤ 6.
**(iii)** h_0 ≡ h_1 ≡ 0, h_2 ≢ 0, pair (2,3): β_2 = β_3 + 30 ∈ [30, 32],
    h_3 ≢ 0, β_3 ≤ 2.

Branches (ii)/(iii) demand the curve y ↦ d̃(y) lie identically on the
irreducible weight-20 hypersurface h31 = 0 (resp. also on h_1 = 0) — very
rigid; candidates for a direct Gröbner kill.

## The g-cascade: branch (i) upgraded to an EXACT reformulation (NEW)

For a = 0 the full identity f31(d̃,Φ̃) ≡ 0 (all 251 slices) is **equivalent** to
the following finite cascade (divide out t³⁰ blockwise; v_t(e) = 0 lets e²¹
be cancelled mod t-powers):

    g_1 := h_0(d̃) / t³⁰                                (t³⁰ | h_0(d̃);  deg g_1 ≤ 10)
    g_{ℓ+1} := (e³ g_ℓ + u^ℓ h_ℓ(d̃)) / t³⁰    ℓ = 1…6  (t³⁰ | numerator; deg ≤ 10)
    e³ g_7 + u⁷ h_7(d̃) = 0                              (exactly; a deg-≤40 object)

Derivation: 0 = Σ_f t³⁰f u^f e^{21−3f} h_f. The f = 0 term forces t³⁰ | e²¹h_0,
and a = 0 gives t³⁰ | h_0. Dividing by t³⁰ and regrouping:
0 = e^{18}(e³g_1 + u h_1) + Σ_{f≥2} t^{30(f−1)} u^f e^{21−3f} h_f, forcing
t³⁰ | e³g_1 + u h_1 (again a = 0), and so on; after the seventh division the
sum terminates and the last line must vanish identically.
Equation count: 7 blocks × 30 + 41 = **251** ✓ (matches the slice count — the
cascade is the whole system, not just a necessary condition).
Degree bookkeeping: every cascade object has y-degree ≤ 40; the coefficients
of g_ℓ are polynomials of degree 10 + 3(ℓ−1) in the 32 window unknowns.

This is the tractable exact formulation T5_NOTES.md §"incremental" was missing:
no degree-25 anchor ever appears; the deepest object is degree 31 in the
unknowns, and the first block (t³⁰ | h_0(d̃)) is 30 equations of degree 10.

## Branch (ii) has its own (tighter) cascade

With h_0(d̃) ≡ 0 and a = 0 the identity Σ_{f≥1} t³⁰f u^f e^{21−3f} h_f ≡ 0
divides down the same way: t³⁰ | h_1(d̃) (note deg h_1(d̃) ≤ 36, so the
quotient g'_2 := h_1/t³⁰ has **deg ≤ 6**), then t³⁰ | e³g'_ℓ + u^{ℓ−1} h_ℓ for
ℓ = 2…6 (numerators deg ≤ 36, quotients deg ≤ 6), then the exact terminal
identity. So branch (ii) = {h_0(d̃) ≡ 0 (41 eqs)} + {t³⁰ | h_1 (30 eqs)} +
five more 30-blocks + terminal — even more overdetermined than branch (i),
with smaller objects. Same Singular machinery applies. Branch (iii) likewise
(h_0 ≡ h_1 ≡ 0, cascade starts at h_2, quotient degrees ≤ 2).

## The a ≥ 1 strata: pair caps

For a ≥ 1 (v = 30 − 3a), pairs (f₁ < f₂) need β_{f₁} − β_{f₂} = (30−3a)(f₂−f₁)
with β_f ≤ 40 − 4f. Non-consecutive pairs need 2(30−3a) with β_{f₁} ≤ 40 − 4f₁:
impossible for a ≤ 3 (gap ≥ 42 > 40); for a ≥ 4 they start to open up.
Consecutive pairs (f, f+1) need β_f ≥ 30 − 3a, allowed only while
40 − 4f ≥ 30 − 3a, i.e. **f ≤ (10 + 3a)/4**:

    a:                1  2  3  4  5  6  7  8  9  10
    v = 30−3a:       27 24 21 18 15 12  9  6  3   0
    max f (consec):   3  4  4  5  6  7  7  7  7   7

Each stratum is a finite list of pairing conditions on the SAME small objects;
a = 10 is degenerate (d̃₋₁ = c·t¹⁰ exactly — substitute directly). All strata
share the window unknowns, so one Singular framework covers them. Not started.

## What remains open

1. **Branch (i):** show the g-cascade has no window solution. First target:
   block 1 + block 2 over F_p (30 + 30 equations, degrees 10 and 13, gauge
   fixed) — dimension should die early; then certify in char 0.
2. **Branches (ii)/(iii):** show h_0(d̃) ≡ 0 (plus h_1(d̃) ≡ 0 for (iii), plus
   the β-pairings) is impossible on the window with e ≢ 0.
3. **a ≥ 1 strata** (v = 30 − 3a ∈ {0,3,…,27}): same polygon game with e = t^a ê;
   the pair analysis must be redone per a (the caps β_f ≤ 40 − 4f still apply,
   plus β_f inherits ≥ (dm1-degree of each h_f monomial)·a lower bounds — the
   min-degree structure of h_f in dm1 will matter). Not started.
4. **f37:** grading is non-uniform in d₋₁ (powers 18,15,12,9,6,4,2,0); needs its
   own reformulation. Not started.
5. **Subcase (1) windows** (deg ≤ 3w): same statements with caps 60 − 6f; the
   pair-(f,f+1) exclusion then fails for f ≥ 3 only at f ≥ 5 — the case tree
   is larger but analogous. Not started.
