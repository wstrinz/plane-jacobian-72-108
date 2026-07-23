# T5_F37_GRADED — the bigraded structure of f37 and the multi-place transport (2026-07-21)

Continues T5_NP.md (f31 t-adic tree) and T5_MULTIPLACE.md (f31 joint cascade).
Every computational statement below is verified exactly by `t5_f37_verify.py`
(9 check groups, all pass); the group used is cited as [Vn]. Unless stated,
everything is for the **f37 branch, subcase (2)** stripped windows:
deg d̃2 ≤ 4, d̃1 ≤ 6, d̃0 ≤ 8, e := d̃₋₁ ≤ 10, e ≢ 0;  t := y+1,
q := 2048y⁴−512y³+320y²−240y+195,  Φ̃ = c·t³⁰·q with c = −1/6630,
v_t(Φ̃) = 30, v_q(Φ̃) = 1, deg Φ̃ = 34 [V1]. The identity under attack is
f37(d̃, Φ̃) ≡ 0.

**Headline results.**
1. f37 carries a TWO-EDGE bigraded structure (Newton-polygon slopes −3 and −2),
   with an exact two-substitution normal form and a working cascade (§1–§3).
2. The q-place Lemma B and the σ-locus kill (Theorem 2) transport verbatim
   (§4, §7); the σ-locus of f37 is dead in every stratum with d̃2 ≢ 0.
3. **f37 vanishes IDENTICALLY on {d̃1 ≡ 0, d̃2 ≡ 0}** (d̃0, e free) — a genuine
   2-function solution family of the f37 identity inside the windows [V6].
   f31's Lemma 1 (no total degeneration) is FALSE for f37, and *no argument
   using only the f37 identity plus the window caps can ever close the f37
   branch*: an extra constraint from the original system is REQUIRED (§6).
4. The f31 degree-starvation kills (T5_MULTIPLACE Props 1–3) do NOT fire for
   f37: zero of the 21 joint strata die that way — the degree caps are ~2×
   f31's while the terminal q-injection is weaker (ê² vs ê³) (§8).

## 1. The bigraded decomposition [V2] (PROVEN)

    f37 = Σ_{f=0}^{7} Φ^f · dm1^{p_f} · h_f(d2, d1, d0, dm1),

with p_f maximal (dm1 ∤ h_f) and h_f weighted-homogeneous of weight
wt_f = 134 − 17f − 5p_f under w(d2,d1,d0,dm1) = (2,3,4,5). Data
(`f37_graded.txt`, same format as `f31_graded.txt`):

    f            0    1    2    3    4    5    6    7
    p_f         18   15   12    9    6    4    2    0
    wt_f        44   42   40   38   36   29   22   15
    terms      145  124  106   88   78   51   25    1
    cap sub2    88   84   80   76   72   58   44   30   (deg_y h_f(d̃) ≤ 2·wt_f)
    cap sub1   132  126  120  114  108   87   66   45   (deg_y h_f(d̃) ≤ 3·wt_f)

The claimed dm1-powers 18,15,12,9,6,4,2,0 are confirmed: **step −3 for
f = 0…4, step −2 for f = 4…7**. The full Φ-coefficients k_f = dm1^{p_f}h_f
have weight 134 − 17f, so deg_y k_f(d̃) ≤ 268 − 34f (sub2).
h_0 = f37|_{Φ=0}/dm1^{18} is exactly the known pure-d factor h37
(weight 44, 145 terms) [V2].

## 2. The Newton polygon: two edges, and the exact normal form (PROVEN)

The lower convex hull of {(f, p_f)} has vertices **(0,18), (4,6), (7,0)**;
edge slopes **−3** (f ∈ [0,4]) and **−2** (f ∈ [4,7]); all eight points lie
on the hull boundary — no interior points [V3]. Consequently the single
substitution w = Φ/dm1³ no longer produces dm1-free cofactors: the points
(f, p_f) for f = 5, 6, 7 lie strictly ABOVE the slope-3 line through (0,18)
(p_f = 4, 2, 0 vs the line's 3, 0, −3), so the best single-w form is
f37 = dm1¹⁸·H̃(d̃, dm1, w) with H̃ = Σ_{f≤4} h_f W^f + dm1·h_5 W⁵ +
dm1²·h_6 W⁶ + dm1³·h_7 W⁷ — "w is a root of H̃", but the tail coefficients
carry e-powers, which shift every tail valuation by multiples of a := v_t(e)
(this is precisely the hull correction in §9). The clean split is the
**two-substitution normal form** [V4]:

    f37 = dm1¹⁸ · A(d̃, w)  +  Φ⁵ dm1⁴ · B(d̃, z),
    A(d̃,W) := Σ_{f=0}^{4} h_f W^f  (head, slope-3 edge),
    B(d̃,Z) := h_5 + h_6 Z + h_7 Z²  (tail, slope-2 edge),
    w := Φ/dm1³,   z := Φ/dm1².

Equivalently f37 = dm1¹⁸·[A(w) + dm1·w⁵·B(dm1·w)] (note z = dm1·w): the
deviation from an f31-style single-w form is the tail correction
dm1·w⁵·B(z). **Quantified**: v_t(Φ̃⁵e⁴B(d̃,z̃)) ≥ 150 + 4a (a := v_t(e)),
so t-adically the identity below order ≈150 sees ONLY the head A — an
f31-style uniform slope-3 object in the five polynomials h_0…h_4. The
valuation attack at the places (t, q, ∞) works directly on the 8 terms
Φ̃^f e^{p_f} h_f; the hull governs which terms can balance (§9).

## 3. Lemma A′ — the cascade with non-uniform steps (PROVEN [V5][V8])

Fix a := v_t(e) ≤ 9 and write e = t^a ê, t ∤ ê, deg ê ≤ 10 − a. Put
Δ := (Δ_1,…,Δ_7) = (3,3,3,3,2,2,2) (Δ_ℓ = p_{ℓ−1} − p_ℓ),
**δ_ℓ := 30 − a·Δ_ℓ** ≥ 3 > 0, and ε_f := Σ_{ℓ≤f} δ_ℓ. Substituting
e = t^a ê and dividing by t^{18a} (using 30f + a·p_f = 18a + ε_f, all f [V9]):

    0 = Σ_f t^{ε_f} u^f ê^{p_f} h_f(d̃),      u := Φ̃/t³⁰ = c·q.

**Lemma A′.** For a ≤ 9 the identity is equivalent to the finite cascade

    g_1 := h_0(d̃) / t^{δ_1}                                (t^{δ_1} | h_0 forced)
    g_{ℓ+1} := (ê^{Δ_ℓ} g_ℓ + u^ℓ h_ℓ(d̃)) / t^{δ_{ℓ+1}}    ℓ = 1…6   (forced)
    ê² g_7 + u⁷ h_7(d̃) = 0                                  (exact terminal)

with degree caps [V8]

    deg g_1, g_2, g_3, g_4 ≤ 58 + 3a,   deg g_5 ≤ 58 + 2a,
    deg g_6 ≤ 48 + 2a,   deg g_7 ≤ 38 + 2a,
    terminal balance: deg(ê²g_7) ≤ 2(10−a) + 38 + 2a = 58 = deg(u⁷h_7) cap.

*Proof.* Identical to T5_MULTIPLACE Lemma A: mod t^{δ_1} only f = 0
survives and t ∤ ê forces t^{δ_1} | h_0; substitute, divide by t^{δ_1},
iterate. The only change is the bookkeeping of the non-uniform multipliers
ê^{Δ_ℓ} and divisors t^{δ_ℓ}; the telescoping is verified symbolically for
every a = 0…9 [V5], the degree recursion in [V8]. ∎

For a = 0 this is the direct f37 analogue of the f31 g-cascade
(δ_ℓ = 30 throughout, deg g_ℓ ≤ 58, 58, 58, 58, 58, 48, 38 — vs f31's ≤ 10).
Equation count at a = 0: 4·30 + 3·30 + 59 = 269 ✓ (= slice count
y⁰…y²⁶⁸ of the stripped identity, weight 2·134).

**a = 10 (degenerate stratum): the slope-3 edge ACTIVATES.** deg e ≤ 10
forces e = C·t¹⁰ exactly (C ∈ K^×). Then 10·p_f + 30f = 180 for ALL
f ≤ 4 (this is precisely "−30/a = −3 = head slope"), and the identity
divided by t¹⁸⁰ becomes

    0 = E(d̃) + t¹⁰·u⁵C⁴h_5 + t²⁰·u⁶C²h_6 + t³⁰·u⁷h_7,
    E(d̃) := Σ_{f=0}^{4} u^f C^{18−3f} h_f(d̃) = C¹⁸·A(d̃, u/C³),

i.e. the whole head collapses to the edge polynomial A evaluated at the
t-unit u/C³, and a three-block t¹⁰-cascade remains:
G_1 := E/t¹⁰ (deg ≤ 78), G_2 := (G_1 + u⁵C⁴h_5)/t¹⁰ (deg ≤ 68),
G_3 := (G_2 + u⁶C²h_6)/t¹⁰ (deg ≤ 58), terminal G_3 + u⁷h_7 = 0.
(Unlike f31's (10,0) stratum, this is NOT t-free: the two hull slopes
differ, so only the head becomes t-degenerate.)

## 4. Lemma B′ — the q-place rains down unchanged (PROVEN [V9])

a_q := v_q(e) ∈ {0,1,2} and a + 4a_q ≤ 10: the **same 21 joint strata**
(a, a_q) as f31. Reducing 0 = Σ_f t^{ε_f}u^f ê^{p_f}h_f mod q kills every
f ≥ 1 term (u^f = c^f q^f), so q | ê¹⁸h_0(d̃):

**Lemma B′.** If a_q = 0 (any a ≤ 9) then q | h_0(d̃) and q | g_ℓ for every
ℓ = 1…7. (Same two-line induction as f31's Lemma B: q | g_ℓ and q | u^ℓ
give q | t^{δ_{ℓ+1}}g_{ℓ+1}, and q ∤ t.)

[V9] confirms the seed is genuinely forcing: 6630⁷·f37(d̃,Φ̃) ≡
6630⁷·e¹⁸h_0(d̃) (mod q) as an exact polynomial congruence, and on a random
window instance q ∤ h_0(d̃) — so f37(d̃,Φ̃) ≡ 0 really imposes the
nontrivial divisibility q | h_0(d̃) when q ∤ e.
q-side caps: δ_f := v_q(h_f(d̃)) ≤ ⌊2wt_f/4⌋ = (22,21,20,19,18,14,11,7) [V8].

## 5. The terminal collapse chain [V6] (PROVEN)

    h_7          = 221184 · d̃1⁵                       (fifth power; f31: d1²)
    h_6|_{d̃1=0} = −82944 · d̃2 · σ⁵,   σ := 4d̃0 − d̃2²  (NEW factor d̃2; σ⁵ not σ²)
    h_5|_{d̃1=0, σ=0} = 131072 · d̃2² · e⁵              (f31: 2048·e², no d̃2)
    and  **d̃2 | h_f|_{d̃1=0}  for EVERY f = 0…7**.

Terminal case tree (a ≤ 9; T-branch names as in T5_MULTIPLACE §5):

* **T1 (d̃1 ≢ 0):** exact terminal  ê² g_7 = −221184 c⁷ q⁷ d̃1⁵,  with forced
  v_q(g_7) = 7 − 2a_q + 5·v_q(d̃1) and, for p ∉ {t, q}:
  2v_p(ê) + v_p(g_7) = 5v_p(d̃1) — a *fifth-power* multiplicity identity.
* **T2 (d̃1 ≡ 0):** g_7 = 0 and the level-6 line is exact:
  ê² g_6 = 82944 c⁶ q⁶ d̃2 σ⁵. Sub-split on the FACTORED right side:
  - **σ ≡ 0, d̃2 ≢ 0** (the σ-locus): dead by Theorem 2′ below (§7).
  - **d̃2 ≡ 0**: forces d̃0 arbitrary — this is the FREE FAMILY (§6);
    the f37 identity is satisfied identically. Not killable here.
  - **d̃2σ ≢ 0**: open; v_q(g_6) = 6 − 2a_q + v_q(d̃2) + 5v_q(σ),
    deg g_6 ≤ 48 + 2a.

## 6. The free family — f37's fundamental difference (PROVEN [V6])

**Fact.** f37|_{d̃1≡0, d̃2≡0} = 0 identically (as a polynomial in d̃0, e, Φ).
Equivalently: every monomial of f37 is divisible by d1 or by d2.
(Two independent confirmations: symbolic substitution of the 618-term
polynomial, and d̃2 | h_f|_{d̃1=0} for all f [V6]. A weight check makes the
"pure (dm1,Φ)" part impossible a priori — 5m + 17p = 134 admits only
(m,p) = (20,2), (3,7), and both monomials are absent — but the full
d1/d2-divisibility statement is stronger and is what [V6] verifies.)

**Consequences.**
1. Any window tuple (d̃2, d̃1, d̃0, e) = (0, 0, φ, ψ) with arbitrary
   φ (deg ≤ 8) and ψ (deg ≤ 10) satisfies f37(d̃, Φ̃) ≡ 0 **exactly**.
   The f37 identity has a 20-dimensional (coefficient count) solution
   space inside the sub2 windows.
2. f31's Lemma 1 ("h_5, h_6, h_7 cannot all vanish when e ≢ 0") is FALSE
   for f37: on d̃1 = d̃2 = 0 ALL EIGHT h_f(d̃) vanish identically.
3. **No infeasibility proof for the f37 branch can be extracted from the
   f37 identity + window caps alone** — the campaign must import at least
   one more fact from the original (8,28) system to exclude (or handle)
   the locus {d̃1 ≡ 0, d̃2 ≡ 0}. Candidates: a proven nonvanishing
   (d̃2 ≢ 0 or d̃1 ≢ 0) from the derivation of the d_k, or a second
   independent identity binding d̃0, e on that locus.
4. **Numerics reconciliation (AUDIT ITEM).** The T3 positive floors for the
   f37 cells do not contradict this: `jetlift.py::make_base` samples the
   base (a0, b0, c0) from a continuous complex normal, so the degenerate
   bases a0 = b0 = 0 that seed this family are hit with probability zero —
   the family was simply outside the explored domain (the gradient filter
   further biases toward generic bases). But it
   means the T3 verdict "f37 infeasible" is, as a statement about the bare
   f37 window identity, **provably false on the degenerate locus** — the T3
   conclusion must be re-scoped to "infeasible off {d̃1 = d̃2 = 0}", and the
   master-identity endgame (f31 ≡ 0 or f37 ≡ 0) must treat the free family
   as a live branch to be excluded by other constraints of the system.

## 7. The σ-locus master identity and Theorem 2′ (PROVEN, uses Mason–Stothers)

On the σ-locus {d̃1 ≡ 0, d̃0 = d̃2²/4} the h_f collapse into a factored
quintic along the slope-3 edge [V7]:

    f37|_σ = 64 · d̃2² · e⁹ · ( 32·A′⁴·B′ − 27·e¹⁷ ),
    A′ := 2Φ̃ + 3 d̃2 e³,    B′ := 4Φ̃ + 3 d̃2 e³    (deg A′, B′ ≤ 34),
    B′ − A′ = 2Φ̃,           2A′ − B′ = 3 d̃2 e³,

where the underlying edge-quintic factorization is
2048X⁵ + 13824X⁴ + 36864X³ + 48384X² + 31104X + 7776 = 32(2X+3)⁴(4X+3)
(the f37 analogue of f31's 512(X+1)⁴(4X−5)) [V7]. Since e ≢ 0, the σ-locus
identity with **d̃2 ≢ 0** is equivalent to

    **32 A′⁴ B′ = 27 e¹⁷** .

(If d̃2 ≡ 0 the prefactor kills everything — that is the free family again,
since then d̃0 = d̃2²/4 = 0.)

**Theorem 2′ (σ-locus empty for d̃2 ≢ 0).** No pair (d̃2, e) with d̃2 ≢ 0,
e ≢ 0, deg d̃2 ≤ 4, deg e ≤ 10 satisfies 32A′⁴B′ = 27e¹⁷. Hence the branch
{d̃1 ≡ 0, σ ≡ 0, d̃2 ≢ 0} is infeasible **in every joint stratum**,
including a = 10.

*Proof — transport of T5_MULTIPLACE Theorem 2.* The proof there used only:
(α) the exponent pattern A⁴B = const·e¹⁷; (β) two linear combinations of
A, B equal to const·Φ̃ and const·d̃2e³; (γ) v_t(Φ̃) = 30, v_q(Φ̃) = 1,
deg Φ̃ = 34; (δ) the caps deg d̃2 ≤ 4, deg e ≤ 10; (ε) Mason–Stothers.
All five inputs are unchanged here — only the constants in (β) differ
(B′ − A′ = 2Φ̃ and 2A′ − B′ = 3d̃2e³ replace 5A + B = 9Φ̃, 4A − B = 9d̃2e³),
and no step of the f31 proof used those constants beyond nonvanishing.
Explicitly:
(ii) A′ = A₀ constant ≠ 0: 2Φ̃ = B′ − A₀ = (27/32A₀⁴)e¹⁷ − A₀; e is
non-constant; differentiating, e¹⁶ | Φ̃′ = c t²⁹(30q + tq′) with
deg(30q + tq′) = 4, forcing e = γt (16a ≤ 29 ⇒ a = 1, unit cofactor);
evaluating at t = 0 (v_t(Φ̃) = 30) gives A₀ = 0. ⨯
(iii) B′ = B₀ constant ≠ 0: A′⁴ = (27/32B₀)e¹⁷ forces 4 | v_p(e) ∀p, so
e = γE⁴, A′ = c₂E¹⁷; differentiating 2Φ̃ = B₀ − A′ gives E¹⁶ | Φ̃′, so
E = γ₂t, A′ = c₂′t¹⁷; the t⁰-coefficient of 2Φ̃ = B₀ − c₂′t¹⁷ gives
B₀ = 0. ⨯
(iv) A′, B′ non-constant: all prime factors of A′B′ divide e; a common
prime of A′, B′ would divide B′ − A′ = 2Φ̃, so gcd concentrates in {t,q};
4·deg A′ + deg B′ = 17·deg e and deg(B′ − A′) = 34, deg(2A′ − B′) ≤ 34
force (same arithmetic as f31, using 17·deg e ≤ 170) deg A′ = deg B′ = 34,
deg e = 10. t-place: 4α_t + β_t = 17a with v_t(B′−A′) = 30 forces
α_t = β_t = τ, 5τ = 17a, a ∈ {0, 5} (the α_t ≠ β_t escapes die exactly as
in f31: a = 9 gives s_A = 16/17 ∉ Z, a = 10 gives α_t = 35 > 34 or
τ = 34 > 30). q-place: 4α_q + β_q = 17a_q with v_q(B′−A′) = 1: a_q = 0
gives α_q = β_q = 0; a_q ∈ {1,2} has no integral/degree-feasible option
(identical arithmetic). Remaining cases (a, a_q) = (0,0): s_A = 8, s_B = 2,
and Mason–Stothers on the coprime triple A′ + 2c t³⁰q = B′ gives
max deg = 34 > N₀ − 1 ≤ (8/4) + 2 + 1 + 4 − 1 = 8. ⨯
(a, a_q) = (5,0): A′ = t¹⁷A₁, B′ = t¹⁷B₁, s_A = 4, s_B = 1,
A₁ + 2c t¹³q = B₁, max deg = 17 > (4/4) + 1 + 1 + 4 − 1 = 6. ⨯  ∎

This closes the f37 "T3-terminal" σ-branch globally, exactly as for f31 —
but note the d̃2 ≢ 0 hypothesis, absent in f31's Theorem 2, is essential
(f31's case (i) d̃2 ≡ 0 was a 17∤150 contradiction; f37's d̃2 ≡ 0 is the
free family, not a contradiction).

## 8. The no-kill audit: degree starvation does not fire (PROVEN [V8])

f31's Props 1–3 (T5_MULTIPLACE §3) killed six strata because the cascade
caps (10 + 3a) were smaller than the q-degree the terminal forces
(4·(7 − 3a_q) etc.). For f37 the corresponding inequalities NEVER hold:

    level 7 (kill d̃1):  needs 4(7 − 2a_q) > 38 + 2a  — max LHS 28 < 38.  never
    level 6 (kill σ):    needs 4(6 − 2a_q) > 48 + 2a  — max LHS 24 < 48.  never

[V8] checks all 21 strata: the starvation kill set is **empty**. Two
compounding reasons: f37's window caps are ≈ 2× f31's (2wt_f with
wt_0 = 44 vs 20), and the terminal multiplier is ê² (q-injection
2a_q) instead of ê³. So the f31 §3 mechanism does not transport, and the
survivor map is: **all 21 joint strata keep T1 (and T2 with d̃2σ ≢ 0)
open**; only the σ-sublocus (§7) and nothing else is dead; the free family
(§6) is structurally alive.

## 9. Newton-polygon pair conditions at the three places (PROVEN, bookkeeping)

The t-adic valuation of term f is ε_f + v_t(u^f ê^{p_f} h_f) —
increments δ_{f} between consecutive terms: 30 − 3a for f = 1…4,
30 − 2a for f = 5…7. The min over f must be attained at least twice
(else the single minimal term survives). Balances for f₁ < f₂:

    β_{f₁} − β_{f₂} = 30(f₂ − f₁) − a·(p_{f₁} − p_{f₂}),   β_f := v_t(h_f(d̃)),

the 2D-hull-corrected form of T5_NP Lemma 2 (for f31, p-drop uniform 3
gave (30−3a)(f₂−f₁); here the drop depends on which edge the pair spans).

**a = 0 case tree** (caps β_f ≤ 2wt_f): pairs need β_{f₁} = 30(f₂−f₁) + β_{f₂}:
* (f, f+1), f = 0…6: all a priori allowed (caps ≥ 44 ≥ 30).
* (f, f+2): needs β_{f₁} ≥ 60: allowed for f₁ = 0…4 (caps 88…72), dead for
  f₁ = 5 (58 < 60).
* gaps ≥ 3: need β_{f₁} ≥ 90 > 88: impossible.
* minimal-index bound: if h_0(d̃) ≢ 0 the minimal pair index f₁ ≤ 2
  (else β_0 ≥ 30f₁ > 88); each identically-vanishing h_f below f₁ costs a
  41–89-equation window condition, as in f31 branches (ii)/(iii).
So the a = 0 tree has 12 candidate pairs (vs f31's 3 branches) — coarser,
but the cascade of §3 subsumes it entirely (as for f31).

**q-place** (increments 1 − a_q·Δ_f per step): for a_q = 0 the min forces
δ_0 ≥ 1 — Lemma B′ re-derived. For a_q = 1 the increments are −2 (head
edge) / −1 (tail edge); for a_q = 2 they are −5 / −3; in both cases the
term f = 7 is the unconstrained minimum (f + a_q·p_f is minimal at f = 7
uniquely), so a balance forces the tail δ_f's upward — e.g. a_q = 2, pair
(6,7): δ_6 = δ_7 + 3 with δ_6 ≤ 11.

**∞-place**: v_∞ = −deg; with full-degree windows every term caps at
deg = 268 = 2·134 (weighted homogeneity), and the y²⁶⁸ anchor slice is
f37(λ2, λ1, λ0, λ₋₁, −1024/3315) = 0 on the leading coefficients — same
bookkeeping as T5_MULTIPLACE §7, unexploited.

## 10. Transport map — what carries over from the f31 campaign

| f31 ingredient | f37 status |
|---|---|
| graded decomposition, weights, caps | ✓ transports, now bigraded/two-edge [V2–V4] |
| single w = Φ/dm1³ substitution, H(d̃,W) root condition | ✗ — replaced by two-substitution normal form + hull balances (§2, §9) |
| Lemma A cascade | ✓ Lemma A′, non-uniform steps Δ = (3,3,3,3,2,2,2), caps 58+3a … 38+2a (§3) |
| a = 10 degenerate stratum (t-free q-cascade) | partially — slope-3 edge activates; mixed t¹⁰-block cascade, NOT t-free (§3) |
| Lemma 1 (no total degeneration) | ✗ **FALSE** — free family {d̃1 = d̃2 = 0} (§6) |
| Lemma B (q | h_0, q rains down) | ✓ verbatim, Lemma B′ (§4) |
| Props 1–3 terminal starvation (6 strata killed) | ✗ kill set EMPTY — caps too large (§8) [V8] |
| Lemma C + Theorem 2 (σ-locus master identity, Mason–Stothers) | ✓ transports verbatim: 32A′⁴B′ = 27e¹⁷, Theorem 2′, requires d̃2 ≢ 0 (§7) |
| survivor strata map | all 21 strata survive for f37 (minus σ-sublocus) (§8) |
| T5_NP a = 0 three-branch tree | 12-pair tree, subsumed by cascade (§9) |

## 11. Status: PROVEN vs CONJECTURED

**PROVEN (verified [V1]–[V9] + arguments above):** the decomposition and
all data of §1; the hull and normal form of §2; Lemma A′ and the a = 10
edge-activation cascade (§3); Lemma B′ (§4); the collapse identities and
case tree of §5; the free family and consequence 2 of §6; Theorem 2′ (§7,
external ingredient: Mason–Stothers, as in f31); the no-kill audit (§8);
the pair-condition bookkeeping (§9).

**CONJECTURED / OPEN:**
1. How to exclude (or absorb) the free family {d̃1 ≡ d̃2 ≡ 0} — REQUIRES an
   input beyond the f37 identity (§6.3). First step: audit the derivation
   of the window system for a provable nonvanishing (is d̃2 ≡ 0 ∧ d̃1 ≡ 0
   consistent with the original G-system + d₋₁ ≢ 0?). This is now the
   gating question for the whole f37 branch.
2. Killing T1/T2 in any stratum. Most structured targets: the fifth-power
   terminal ê²g_7 = −221184c⁷q⁷d̃1⁵ (prime-multiplicity parities
   2v_p(ê) + v_p(g_7) = 5v_p(d̃1) + 7·[p=q] — Mason–Stothers candidates in
   the style of Theorem 2′(iv)), and the level-6 line ê²g_6 = 82944c⁶q⁶d̃2σ⁵.
3. Exploiting the head-only window: below t-order ~150 the identity is the
   pure slope-3 system in h_0…h_4 (§2) — the first four cascade blocks
   (30 equations each at a = 0, objects of degree ≤ 88 in y, degree 44–36
   weight in d̃) are an f31-shaped Gröbner target.
4. Subcase (1) windows: caps 3wt_f = (132,…,45), deg e ≤ 15, a_q ≤ 3;
   all bookkeeping (cascade caps, starvation margins) must be redone; the
   free family persists verbatim (it is cap-independent).
5. The ∞-place and the Q̄-splitting of q (4 conjugate places): untouched.

## 12. Verification map

`t5_f37_verify.py` (all exact, sympy over Q; run from this directory;
~5–10 min, the V4/V6/V7 expansions of the 618-term f37 dominate):

| check | verifies | used in |
|---|---|---|
| V1 | q irreducible, q(−1)=3315, v_t(Φ̃)=30, v_q(Φ̃)=1, deg Φ̃=34 | §0 |
| V2 | decomposition, p_f, weights, term counts, dm1∤h_f, h_0=h37 | §1 |
| V3 | hull vertices (0,18),(4,6),(7,0); slopes −3,−2; all points on boundary | §2 |
| V4 | two-substitution normal form; caps 2wt_f / 3wt_f | §2 |
| V5 | cascade telescoping for every a = 0…9 | §3 |
| V6 | h_7, h_6|, h_5| collapses; d̃2 | h_f|_{d1=0} ∀f; **free family** | §5, §6 |
| V7 | σ-locus master identity, A′/B′ combos, edge-quintic = 32(2X+3)⁴(4X+3) | §7 |
| V8 | cascade degree caps; terminal balance; starvation kill set EMPTY; q-caps | §3, §8 |
| V9 | 30f + a·p_f = 18a + ε_f (symbolic a); mod-q congruence exact + forcing nontrivial on a random instance | §3, §4 |
