# T5_MULTIPLACE — the joint (t, q)-adic cascade for f31, subcase (2) (2026-07-21)

> **STATUS (2026-07-22):** CONDITIONAL / SCOPE-RESTRICTED. The `(a_t, a_q)` survivor maps below are valid ONLY as geometrically-q-coprime (uniform-qʳ) statements: after base change q splits into four geometric places and the scalar `a_q` ledger is not field-stable (`FIELD_SPLIT_AUDIT.md`; `STATE.md` L5–16). The authoritative replacements are the split-place ledger (`split_place_ledger.py` / `split_place_ledger.json`) and the cascade-engine artifacts (`cascade_cones_qt_inf_rl.json`).

Continues T5_NP.md (single-place t-adic tree). All computational inputs are
verified exactly by `t5_multiplace_verify.py` (8 check groups, all pass); the
check group used by each statement is cited as [Vn]. Everything below is for
the **f31 branch, subcase (2)** stripped windows:
deg d̃2 ≤ 4, d̃1 ≤ 6, d̃0 ≤ 8, e := d̃₋₁ ≤ 10, e ≢ 0, and the identity under
attack is

    0  =  Σ_{f=0}^{7} Φ̃^f e^{21−3f} h_f(d̃),      deg_y h_f(d̃) ≤ 40 − 4f,   [V2]

with h_f from `f31_graded.txt`. Everything in §§1–7 is **PROVEN** (the only
external ingredient is the Mason–Stothers theorem, cited where used); §8 lists
what is open.

## 0. The two finite places, and one new fact

Write t := y+1 and  **q := 2048y⁴ − 512y³ + 320y² − 240y + 195**.  Then

    Φ̃ = c · t³⁰ · q,     c := −1/6630,     u := Φ̃/t³⁰ = c·q .

**Fact ([V1]).** q is irreducible over Q, q(−1) = 3315 ≠ 0, disc(q) ≠ 0. Hence
t and q generate distinct finite places of Q(y), and the divisor of Φ̃ is
concentrated in exactly these two places:

    v_t(Φ̃) = 30 (exactly),      v_q(Φ̃) = 1 (exactly),      deg Φ̃ = 34.

Notation: a := a_t := v_t(e) ∈ [0,10], a_q := v_q(e) ∈ {0,1,2} (deg q = 4 and
deg e ≤ 10 force a_q ≤ 2), and since t^a q^{a_q} | e,

    a + 4·a_q ≤ 10 .                                                     [V1]

The q-side degree caps: since deg_y h_f(d̃) ≤ 40 − 4f,

    δ_f := v_q(h_f(d̃)) ≤ (40−4f)/4 = 10 − f   whenever h_f(d̃) ≢ 0.      [V2]

Joint strata: 21 pairs (a, a_q) with a + 4a_q ≤ 10. The t-side data of
T5_NP.md (v_t(w) = 30 − 3a for w := Φ̃/e³, pair caps β_f ≤ 40 − 4f) is used
unchanged; the new leverage is that **u^f = c^f q^f injects q-multiplicity f
into the f-th term**, so q-multiplicities accumulate down the cascade.

## 1. Lemma A — the cascade exists in every stratum a ≤ 9

Fix a ≤ 9 and write e = t^a ê, t ∤ ê, deg ê ≤ 10 − a, and **v := 30 − 3a ≥ 3**.

**Lemma A.** The identity is equivalent to the finite cascade

    g_1 := h_0(d̃) / t^v                                  (t^v | h_0(d̃) forced)
    g_{ℓ+1} := (ê³ g_ℓ + u^ℓ h_ℓ(d̃)) / t^v     ℓ = 1…6   (t^v | numerator forced)
    ê³ g_7 + u⁷ h_7(d̃) = 0                               (exact terminal)

with every g_ℓ ∈ K[y] of degree ≤ 40 − v = **10 + 3a**.

*Proof.* Substituting e = t^a ê into the identity and pulling out
t^{a(21−3f)}·t^{30f} = t^{21a}·t^{(30−3a)f} from term f gives (dividing by
t^{21a}, K[y] a domain)

    0 = Σ_f t^{vf} u^f ê^{21−3f} h_f(d̃) .                                [V6]

Mod t^v only f = 0 survives: t^v | ê²¹ h_0(d̃), and t ∤ ê gives t^v | h_0(d̃)
— defining g_1. Substituting h_0 = t^v g_1 and dividing by t^v,
0 = ê^{18}(ê³g_1 + u h_1) + Σ_{f≥2} t^{v(f−1)} u^f ê^{21−3f} h_f, and mod t^v
again forces t^v | ê³g_1 + u h_1 — defining g_2; iterating six times exhausts
the sum and leaves the terminal line exactly. Each division is forced, and
conversely the recursion telescopes back to the identity ([V5], symbolic).
Degrees: deg u^ℓ h_ℓ(d̃) ≤ 4ℓ + (40−4ℓ) = 40 and deg ê³g_ℓ ≤ 3(10−a) +
(10+3a) = 40, so every numerator has degree ≤ 40 and every quotient degree
≤ 40 − v = 10 + 3a ([V7]). ∎

For a = 0 this is exactly the g-cascade of T5_NP.md (v = 30, deg g_ℓ ≤ 10).
Note the cascade subsumes all three branches (i)/(ii)/(iii) of the T5_NP a=0
tree (branch (ii) is g_1 ≡ 0, branch (iii) is g_1 ≡ g_2 ≡ 0): no case split
is needed below.

## 2. Lemma B — q-divisibility rains down the cascade (a_q = 0)

**Lemma B.** If a_q = 0 (any a ≤ 9), then q | h_0(d̃) and **q | g_ℓ for every
ℓ = 1…7**. In the joint generic stratum (a, a_q) = (0, 0) in particular:
t³⁰·q | h_0(d̃), i.e.

    h_0(d̃) = t³⁰ · q · G,    deg G ≤ 40 − 30 − 4 = 6 .

*Proof.* Reduce 0 = Σ_f t^{vf} u^f ê^{21−3f} h_f(d̃) mod q: every term f ≥ 1
carries u^f = c^f q^f ≡ 0, so q | ê²¹ h_0(d̃) ([V6]); q is prime and q ∤ ê
(a_q = v_q(ê) = 0), so q | h_0(d̃). Then q | g_1 = h_0/t^v since q ∤ t.
Inductively, t^v g_{ℓ+1} = ê³ g_ℓ + u^ℓ h_ℓ with q | g_ℓ and q | u^ℓ (ℓ ≥ 1)
gives q | t^v g_{ℓ+1}, hence q | g_{ℓ+1}. The (0,0) refinement is the two
coprime divisibilities t³⁰ | h_0 (Lemma A) and q | h_0 against deg h_0(d̃)
≤ 40. ∎

## 3. The terminal collapse chain — six strata die

All three propositions work on the **terminal end** of the cascade, where the
accumulated q-power (q⁷ in u⁷) meets the smallest degree caps. Throughout,
a ≤ 9, and recall v_q(ê) = a_q, deg g_ℓ ≤ 10 + 3a, and the collapse
identities ([V3]):

    h_7(d̃) = 8192 d̃1²,
    h_6(d̃)|_{d̃1≡0} = −3072 σ²,        σ := 4d̃0 − d̃2²   (deg σ ≤ 8),
    h_5(d̃)|_{d̃1≡0, σ≡0} = 2048 e² = 2048 t^{2a} ê² .

**Proposition 1 (level 7).** If  28 − 12a_q > 10 + 3a  (condition **C1**),
then d̃1 ≡ 0 and g_7 = 0.

*Proof.* Terminal: ê³ g_7 = −c⁷ q⁷ · 8192 d̃1². If d̃1 ≢ 0 both sides are
nonzero; taking v_q: 3a_q + v_q(g_7) = 7 + 2·v_q(d̃1), so
v_q(g_7) = 7 − 3a_q + 2v_q(d̃1) ≥ 7 − 3a_q. But g_7 ≠ 0 of degree ≤ 10 + 3a
forces 4·v_q(g_7) ≤ 10 + 3a, contradicting C1. So d̃1 ≡ 0, h_7(d̃) = 0, and
the terminal collapses to ê³g_7 = 0, i.e. g_7 = 0 (domain, ê ≢ 0). ∎

**Proposition 2 (level 6).** If additionally  24 − 12a_q > 10 + 3a
(condition **C2**; note C2 ⇒ C1), then σ ≡ 0 (i.e. d̃0 = d̃2²/4) and g_6 = 0.

*Proof.* The ℓ = 6 cascade line with g_7 = 0 reads 0 = ê³g_6 + c⁶q⁶ h_6(d̃),
and h_6(d̃) = −3072σ² by Prop. 1. If σ ≢ 0: v_q gives
v_q(g_6) = 6 − 3a_q + 2v_q(σ) ≥ 6 − 3a_q, and g_6 ≠ 0 forces
4(6 − 3a_q) ≤ 10 + 3a, contradicting C2. So σ² ≡ 0, hence σ ≡ 0 (char 0
domain), and g_6 = 0. ∎

**Proposition 3 (level 5 kill).** Under C2 the stratum is **infeasible**.

*Proof.* By Props. 1–2, d̃1 ≡ 0, d̃0 = d̃2²/4, g_6 = 0. The ℓ = 5 line reads
0 = ê³g_5 + c⁵q⁵·2048 t^{2a} ê², and cancelling ê² (domain):

    ê · g_5 = −2048 c⁵ · q⁵ t^{2a}  ≠ 0 .

So every irreducible factor of ê divides q⁵t^{2a}; t ∤ ê and v_q(ê) = a_q
force **ê = C·q^{a_q}** with C ∈ Q^×. Hence
deg g_5 = (20 + 2a) − 4a_q. The cap deg g_5 ≤ 10 + 3a then demands
a ≥ 10 − 4a_q. But C2 says 3a < 14 − 12a_q, i.e. a ≤ 4 − 4a_q < 10 − 4a_q.
Contradiction. ∎

**Theorem 1 (main kill).** Every joint stratum with 3a < 14 − 12a_q is
infeasible. Explicitly ([V8]):

    KILLED:  (a, a_q) ∈ { (0,0), (1,0), (2,0), (3,0), (4,0), (0,1) } .

In particular **the joint generic stratum (a_t, a_q) = (0,0) is dead** — and
with it all three branches (i)/(ii)/(iii) of the T5_NP.md a = 0 case tree,
for every window tuple with q ∤ e. The a_t = 0 stratum of T5_NP.md survives
only inside a_q = 2, i.e. only if **q² | e** (so deg e ∈ {8,9,10} and
e = q²·ê₀ with deg ê₀ ≤ 2, t ∤ e).

*Remark.* The mechanism is exactly the "joint degree bookkeeping" hoped for:
the t-cascade compresses everything into objects g_ℓ of degree ≤ 10 + 3a,
while the u-powers pump q-multiplicity 7, 6, 5 into the last three levels;
q has degree 4, so those multiplicities cost 28, 24, 20 degrees — more than
the objects can carry when a is small.

## 4. The σ-locus dies globally (all strata)

Call **σ-locus** the window locus  { d̃1 ≡ 0 and σ = 4d̃0 − d̃2² ≡ 0, e ≢ 0 }
— precisely where the collapse chain of §3 lands, and precisely the locus of
T5_NP.md Lemma 1's degenerate cascade.

**Lemma C (master identity, [V3][V4]).** On the σ-locus the h_f collapse to

    h_f(d̃) = c_f · d̃2^{5−f} e²  (f = 1…5),    h_0(d̃) = c_0 d̃2⁵e² − 6561 e⁴,
    h_6(d̃) = h_7(d̃) = 0,        (c_0,…,c_5) = (−2560, −8192, −7168, 2048, 5632, 2048),

and the quintic Σ_f c_f X^f = 512 (X+1)⁴ (4X−5) factors. Consequently the
full identity restricted to the σ-locus is equivalent (e ≢ 0) to

    **512 · A⁴ · B  =  6561 · e^{17}**,
    A := Φ̃ + d̃2 e³,   B := 4Φ̃ − 5 d̃2 e³     (deg A, deg B ≤ 34),

with the linear-algebra inverses  5A + B = 9Φ̃  and  4A − B = 9 d̃2 e³.

*Proof.* The collapse and the factorization are [V3]; the master identity is
the verified symbolic identity Σ_f Φ^f dm1^{21−3f} h_f|_σ = dm1⁸(512A⁴B −
6561·dm1¹⁷) [V4], and dm1⁸ = e⁸ cancels since e ≢ 0. ∎

**Theorem 2 (σ-locus empty).** No pair (d̃2, e), deg d̃2 ≤ 4, deg e ≤ 10,
e ≢ 0, satisfies 512A⁴B = 6561e^{17}. Hence the σ-locus is infeasible **in
every joint stratum**, including the degenerate (10, 0).

*Proof.* A, B ≠ 0 since A⁴B is a nonzero multiple of e^{17}.

(i) *d̃2 ≡ 0.* Then 2048 Φ̃⁵ = 6561 e^{17}, so 150 = v_t(LHS) = 17a — but
17 ∤ 150 [V8]. Dead. Assume d̃2 ≢ 0 from now on.

(ii) *A constant, A = A₀ ≠ 0.* Then B = (6561/512A₀⁴) e^{17} and
9Φ̃ = 5A₀ + B. e is non-constant (else Φ̃ would be constant). Differentiate:
9Φ̃′ = const·e^{16}e′ with Φ̃′ = c t²⁹(30q + t q′) ≠ 0 and deg(30q + tq′) = 4
[V4]. So e′ ≢ 0 and e^{16} | t²⁹(30q + tq′). Writing e = t^{a}e₁ (t ∤ e₁):
e₁^{16} divides the degree-4 polynomial 30q + tq′, so e₁ is constant; and
16a ≤ 29 gives a ≤ 1, so a = 1 and e = γt. Then B = const·t^{17} and
evaluating 9Φ̃ = 5A₀ + const·t^{17} at t = 0 (where Φ̃ vanishes to order 30)
gives 5A₀ = 0. Contradiction.

(iii) *B constant, B = B₀ ≠ 0, A non-constant.* A⁴ = (6561/512B₀)e^{17}
forces 4·v_p(A) = 17·v_p(e) for every irreducible p, so 4 | v_p(e) for all p:
e = γE⁴ and A/E^{17} is a constant (a rational function whose 4th power is
constant). E is non-constant (else A is). Differentiating 9Φ̃ = 5A + B₀ =
5c₁E^{17} + B₀ gives E^{16} | Φ̃′ = ct²⁹(30q+tq′), so as in (ii) E = γ₂t and
A = c₂t^{17}. Comparing the t⁰-coefficients of 9Φ̃ = 5c₂t^{17} + B₀ (LHS has
v_t = 30) gives B₀ = 0. Contradiction.

(iv) *A, B non-constant.* All prime factors of A and of B divide e. No prime
p ∉ {t, q} divides both A and B (it would divide 5A + B = 9c t³⁰ q). Degrees:
4·deg A + deg B = 17·deg e, and deg(5A+B) = deg 9Φ̃ = 34.
  * If deg A ≠ deg B, the larger equals 34; then deg(4A − B) = 34 =
    deg(9 d̃2 e³) ≤ 4 + 3·deg e ≤ 34 forces deg e = 10, and then the degree
    equation forces the other of deg A, deg B to be 34 too — contradiction.
  * So deg A = deg B =: D ≥ 34 and 5D = 17 deg e ≤ 170, i.e. **D = 34,
    deg e = 10**.
For p | e, p ∉ {t,q}, write μ_p := v_p(e) ≥ 1; then 4v_p(A) + v_p(B) = 17μ_p
with one of v_p(A), v_p(B) zero, so p is *type A* (v_p(A) = 17μ_p/4, 4 | μ_p,
hence μ_p ≥ 4) or *type B* (v_p(B) = 17μ_p). Put s_A := Σ_{type A} μ_p deg p,
s_B := Σ_{type B} μ_p deg p, so a + 4a_q + s_A + s_B = 10.

*The t-place.* 4α_t + β_t = 17a (α_t := v_t(A), β_t := v_t(B)) and
v_t(5A + B) = 30.
  - α_t ≠ β_t forces min(α_t, β_t) = 30: if α_t = 30 < β_t then a ≥ 9; a = 9
    leaves non-(t,q) mass s_A + s_B = 1, incompatible with deg A = 34 =
    30 + (17/4)s_A (type-A primes need μ_p ≥ 4, and s_A = 16/17 ∉ Z); a = 10
    means e = γt¹⁰, so A = ct^{α_t} with α_t = deg A = 34, whence β_t = 34 =
    α_t — contradiction. If β_t = 30 < α_t then α_t = (17a−30)/4 > 30 needs
    a = 10 and α_t = 35 ≠ 34 = deg A. Dead.
  - So α_t = β_t =: τ, 5τ = 17a, hence 5 | a: a ∈ {0, 5} with τ ∈ {0, 17}
    (a = 10 gives τ = 34 > 30 = v_t(5A+B), impossible).

*The q-place.* 4α_q + β_q = 17a_q and v_q(5A + B) = 1. If a_q = 0:
α_q = β_q = 0. If a_q ≥ 1: α_q = β_q would force 5 | 17a_q (impossible for
a_q ≤ 2), so min(α_q, β_q) = 1: (α_q, β_q) = (1, 17a_q − 4) or
((17a_q − 1)/4, 1).

*Case a = 0 (τ = 0).*
  - a_q = 0: 34 = (17/4)s_A and 34 = 17 s_B give s_A = 8, s_B = 2 (and indeed
    8 + 2 = 10). Now gcd(A, B) = 1 and 5A + B = 9c t³⁰ q with A, B
    non-constant: **Mason–Stothers** (char 0: for coprime X + Y = Z, not all
    constant, max deg ≤ N₀(XYZ) − 1, N₀ = number of distinct roots) applies
    with N₀ ≤ (s_A/4) + s_B + 1 + 4 = 2 + 2 + 5 = 9 (type-A primes carry
    μ_p ≥ 4, so A has at most s_A/4 = 2 distinct roots; B at most s_B = 2;
    t³⁰ has 1; q has 4). But max deg = 34 > 9 − 1 [V8]. Dead.
  - a_q = 1: deg A = 34 = 4α_q + (17/4)s_A has no integer solution for
    α_q ∈ {1, 4} (s_A = 120/17 or 72/17). Dead.
  - a_q = 2: (α_q, β_q) = (1, 30) makes deg B ≥ 120 > 34; the other option is
    non-integral. Dead.

*Case a = 5 (τ = 17).* Then 4a_q + s_A + s_B = 5, so a_q ≤ 1.
  - a_q = 0: 34 = 17 + (17/4)s_A and 34 = 17 + 17s_B give s_A = 4, s_B = 1
    (and 5 + 4 + 1 = 10). Divide the t-part out: A = t¹⁷A₁, B = t¹⁷B₁,
    gcd(A₁, B₁) = 1, and 5A₁ + B₁ = 9c t¹³ q with deg A₁ = deg B₁ = 17.
    Mason–Stothers: N₀ ≤ (s_A/4) + s_B + 1 + 4 = 1 + 1 + 5 = 7, but
    max deg = 17 > 7 − 1 [V8]. Dead.
  - a_q = 1: 34 = 17 + 4α_q + (17/4)s_A has no integral solution for
    α_q ∈ {1, 4}. Dead. ∎

**Corollary (global dichotomy).** Any window solution of the f31 subcase-(2)
identity with e ≢ 0 must have  **d̃1 ≢ 0  or  4d̃0 − d̃2² ≢ 0**. In
particular the "T3" terminal branch (d̃1 ≡ 0 and σ ≡ 0) is void in every
surviving stratum, and T5_NP.md Lemma 1's degenerate chain is now closed one
step earlier: h_7 ≡ h_6 ≡ 0 already forces h_5(d̃) = 2048e² ≠ 0 *and* the
whole remaining identity is the (infeasible) master equation — not just a
non-vanishing certificate.

## 5. Map of the surviving strata

15 of 21 joint strata survive §§3–4 ([V8]):

    a_q = 0 : a ∈ {5, 6, 7, 8, 9}  and the degenerate a = 10
    a_q = 1 : a ∈ {1, 2, 3, 4, 5, 6}
    a_q = 2 : a ∈ {0, 1, 2}

For a ≤ 9 the terminal trichotomy is (T3 dead by Thm. 2 in all of them):

**T1** (d̃1 ≢ 0): open iff 28 − 12a_q ≤ 10 + 3a. Exact terminal
    ê³ g_7 = −8192 c⁷ q⁷ d̃1²  with forced valuation
    v_q(g_7) = 7 − 3a_q + 2v_q(d̃1) and 4·v_q(g_7) ≤ 10 + 3a.
**T2** (d̃1 ≡ 0, σ ≢ 0): open iff 24 − 12a_q ≤ 10 + 3a. Then g_7 = 0 and the
    exact level-6 identity  ê³ g_6 = 3072 c⁶ q⁶ σ²  with
    v_q(g_6) = 6 − 3a_q + 2v_q(σ) and 4·v_q(g_6) ≤ 10 + 3a.

Per-stratum consequences (writing ê = q^{a_q}·ē, q ∤ ē, t ∤ ē):

| stratum | open branches | forced structure |
|---|---|---|
| (5,0) | **T2 only** | d̃1 ≡ 0 forced; v_q(σ) = 0; g_6 = q⁶ĝ, deg ĝ ≤ 1, and **ê³ĝ = 3072c⁶σ²** (deg ê ≤ 5, deg σ ≤ 8; ĝ of degree ≤ 1 must make ê³ĝ a perfect-square multiple) |
| (6,0) | T1, T2 | T1: v_q(d̃1) = 0, g_7 = C·q⁷ (cofactor constant), so **C·ê³ = −8192c⁷d̃1²**: every v_p(ê) even, ê = unit·r², d̃1 = unit·r³, deg r ≤ 2. T2: v_q(σ) = 0, g_6 = q⁶ĝ, deg ĝ ≤ 4, ê³ĝ = 3072c⁶σ² |
| (7,0)–(9,0) | T1, T2 | same shapes with slack: g_7 = q⁷ĝ, deg ĝ ≤ 3a − 18; v_q(d̃1) = 0 for a ≤ 8, ≤ 1 for a = 9; ê³ĝ = −8192c⁷d̃1² resp. ê³ĝ′ = 3072c⁶σ² |
| (10,0) | degenerate | e = C·t¹⁰ exactly; the identity is **t-free**: 0 = Σ_f c^f C^{21−3f} q^f h_f(d̃), equivalent to the pure q-cascade G_1 := h_0(d̃)/q, G_{ℓ+1} := (C³G_ℓ + c^ℓ h_ℓ(d̃))/q (deg G_ℓ ≤ 40 − 4ℓ), terminal C³G_7 = −8192c⁷ d̃1²; σ-locus part dead |
| (1,1) | **T2 only** | d̃1 ≡ 0 forced (C1 holds: 16 > 13); v_q(σ) = 0; g_6 = q³ĝ, deg ĝ ≤ 1, and **ē³ĝ = 3072c⁶σ²** (deg ē ≤ 5) |
| (2,1)–(6,1) | T1, T2 | T1: ē³g_7 = −8192c⁷q⁴d̃1², v_q(g_7) = 4 + 2v_q(d̃1), 4v_q(g_7) ≤ 10+3a. T2: ē³g_6 = 3072c⁶q³σ² |
| (0,2)–(2,2) | T1, T2 | e = t^a q² ē with deg ē ≤ 2 − a. T1: **ē³g_7 = −8192c⁷·q·d̃1²**, v_q(g_7) = 1 + 2v_q(d̃1). T2: **ē³g_6 = 3072c⁶σ²** (the q-powers cancel exactly), v_q(g_6) = 2v_q(σ) |

All a_q = 0 survivors additionally carry Lemma B: q | g_ℓ for all ℓ, on top of
the upstream cascade blocks (levels 1–4), which are so far unexploited.

The most rigid targets are **(5,0)** and **(1,1)**: a single open branch,
d̃1 ≡ 0 forced, and the whole terminal end compressed into ê³ĝ = const·σ²
with deg ĝ ≤ 1 — a near-perfect-power condition (for each prime p:
3v_p(ê) + v_p(ĝ) = 2v_p(σ) with Σv_p(ĝ)·deg p ≤ 1) plus five untouched
cascade blocks upstream.

## 6. Newton-polygon pair conditions at the q-place

The analog of T5_NP.md Lemma 2 at q (same proof: a nonzero root w = Φ̃/e³ of
H(d̃, W) = Σ_f h_f(d̃)W^f with v_q(w) = 1 − 3a_q requires the min of
δ_f + (1 − 3a_q)f to be attained at least twice, δ_f := v_q(h_f(d̃)) ≤ 10 − f):

* **a_q = 0** (v_q(w) = 1): pairs (f₁ < f₂) need δ_{f₁} − δ_{f₂} = f₂ − f₁.
  Taking the minimal term: δ_0 ≥ 1, i.e. q | h_0(d̃) (Lemma B re-derived).
  Caps are loose (δ_{f₁} ≤ 10 − f₁), so most pairs are a priori allowed.
* **a_q = 1** (v_q(w) = −2): δ_{f₂} − δ_{f₁} = 2(f₂ − f₁); the cap gives
  f₂ + 2(f₂ − f₁) ≤ 10.
* **a_q = 2** (v_q(w) = −5): δ_{f₂} − δ_{f₁} = 5(f₂ − f₁). Non-consecutive
  pairs need δ_{f₂} ≥ 10 with δ_{f₂} ≤ 10 − f₂: impossible. So only
  **consecutive pairs (f, f+1) with f ≤ 4** and δ_{f+1} = δ_f + 5 survive —
  as tight as the t-side table of T5_NP.md.

Joint statement: the single support S := {f : h_f(d̃) ≢ 0} must
simultaneously admit a slope pair for v_t(w) = 30 − 3a (T5_NP.md caps
β_f ≤ 40 − 4f) *and* a slope pair for v_q(w) = 1 − 3a_q (caps δ_f ≤ 10 − f).
The cascade of §3 is strictly stronger than the pair conditions at the last
three levels, but for the *upstream* levels (f ≤ 4) the pair conditions are
the sharpest per-place information currently available.

## 7. The place at infinity

v_∞ := −deg on K(y). Then v_∞(Φ̃) = −34, v_∞(w) = 3·deg e − 34, and the same
Newton-polygon lemma at ∞ says: **max_f ( deg h_f(d̃) + f·(34 − 3·deg e) )
must be attained at least twice** on S. With deg e = 10 the increment is 4
per f and the caps deg h_f(d̃) ≤ 40 − 4f make deg h_f(d̃) + 4f ≤ 40 for all
f: the top anchor. When the max 40 is attained (all windows at full degree)
the y²⁵⁰-slice of the identity is the anchor equation
f31(λ2, λ1, λ0, λ₋₁, −1024/3315) = 0 on the leading coefficients (leading
coefficient of Φ̃ is −1024/3315 [V1]; cf. T5_NOTES.md item 2 — an irreducible
degree-31 hypersurface, so the anchor alone kills nothing). Note the forced
structures of §5 feed into this: e.g. in (5,0)/(1,1), d̃1 ≡ 0 kills λ1, and
in the T3-free world every survivor keeps at least one of h_6(d̃), h_7(d̃)
nonzero. The ∞-place has **not** yet been exploited beyond this bookkeeping
(open; see §8).

## 8. What remains open

1. **Survivor strata.** Kill (5,0) and (1,1) first: one branch each, terminal
   compressed to ê³ĝ = 3072c⁶σ² with deg ĝ ≤ 1, plus five upstream cascade
   blocks and q | g_ℓ. The perfect-power structure (v_p parities) plus a
   Mason–Stothers count in the style of Theorem 2(iv) looks promising —
   *conjectured killable, not done*.
2. **(10,0)**: the t-free q-cascade (§5) is a fresh small system; combine
   with the t-adic valuations of h_f(d̃)|_{dm1 = Ct¹⁰}. Not started.
3. **T1 branches** ((6,0)–(9,0), (2,1)–(6,1), (0,2)–(2,2)): the
   near-perfect-power identities (d̃1² ∝ ê³-type) should submit to the same
   prime-multiplicity parity analysis. Not started.
4. **The place at ∞** beyond bookkeeping; and the middle place structure of q
   over Q̄ (4 conjugate places — unused; everything above needed only
   irreducibility over Q).
5. **f37** and **subcase-(1) windows**: the whole multi-place scheme
   transports (subcase 1: caps deg h_f(d̃) ≤ 60 − 6f, deg e ≤ 15, so
   a_q ≤ 3 and all the degree bookkeeping — cascade quotient caps, C1/C2 —
   must be redone from scratch; f37 needs its own grading first, cf.
   T5_NP.md item 4). Not started.

## 9. Verification map

`t5_multiplace_verify.py` (all exact, sympy over Q; run from this directory):

| check | verifies | used in |
|---|---|---|
| V1 | q irreducible, q(−1)=3315, disc ≠ 0, v_t(Φ̃)=30, v_q(Φ̃)=1, deg/lc of Φ̃, a_q ≤ 2 | §0, §7 |
| V2 | graded decomposition, weights, caps deg h_f(d̃) ≤ 40−4f ⇒ δ_f ≤ 10−f | §0, §6 |
| V3 | h_7, h_6|, h_5| collapse; σ-locus collapse h_f| = c_f d̃2^{5−f}e² (+ −6561e⁴); quintic = 512(X+1)⁴(4X−5) | §3, §4 |
| V4 | master identity Σ = e⁸(512A⁴B − 6561e¹⁷); 5A+B = 9Φ̃; 4A−B = 9d̃2e³; Φ̃′ = ct²⁹(30q+tq′), deg(30q+tq′) = 4 | §4 |
| V5 | cascade telescoping (sufficiency direction of Lemma A), symbolic | §1 |
| V6 | on a random instance (a=2): t-strip identity; t^v and q both divide G − ê²¹h_0(d̃) | §1, §2 |
| V7 | degree bookkeeping deg g_ℓ ≤ 10+3a for all a, ℓ | §1 |
| V8 | kill-set arithmetic = exactly the six strata; Mason margins 34 > 8, 17 > 6; 17 ∤ 150; survivor list | §3, §4, §5 |

The non-computational ingredients are standard: valuation arithmetic in the
UFD K[y], char-0 domain facts, and the Mason–Stothers theorem (§4(iv) only).
Nothing in this file is conjectural except where marked in §8.
