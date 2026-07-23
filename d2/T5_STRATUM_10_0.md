# T5_STRATUM_10_0 — the degenerate stratum (a_t, a_q) = (10, 0) (2026-07-21)

> **STATUS (2026-07-22):** CONDITIONAL / SCOPE-RESTRICTED. The `(a_t, a_q)` survivor maps below are valid ONLY as geometrically-q-coprime (uniform-qʳ) statements: after base change q splits into four geometric places and the scalar `a_q` ledger is not field-stable (`FIELD_SPLIT_AUDIT.md`; `STATE.md` L5–16). The authoritative replacements are the split-place ledger (`split_place_ledger.py` / `split_place_ledger.json`) and the cascade-engine artifacts (`cascade_cones_qt_inf_rl.json`).

Continues T5_MULTIPLACE.md (§5 row "(10,0)", §8 open item 2). Every
computational input below is verified exactly by `t5_stratum100_verify.py`
(check groups S1–S10, all pass, pure sympy over Q); citations [Sn] refer to
it. Setting: f31 branch, **subcase (2)** stripped windows,
deg d̃2 ≤ 4, d̃1 ≤ 6, d̃0 ≤ 8, e := d̃₋₁, t := y+1,
q := 2048y⁴−512y³+320y²−240y+195 (irreducible/Q, q(−1) = 3315 ≠ 0),
Φ̃ = c·t³⁰·q, c = −1/6630, and the identity under attack is

    0 = Σ_{f=0}^{7} Φ̃^f e^{21−3f} h_f(d̃),     h_f from f31_graded.txt.

**Stratum hypothesis.** a_t = v_t(e) = 10 with deg e ≤ 10 forces
**e = C·t¹⁰ exactly**, C ∈ K^× — one scalar unknown instead of 11
coefficients; q ∤ e is automatic (a_q = 0 ✓). Unknowns: 21 window
coefficients (d̃2: 5, d̃1: 7, d̃0: 9) plus C ≠ 0. Everything below is
**PROVEN** unless explicitly marked otherwise; the field is any field K of
characteristic 0 (all arguments are valuation-theoretic or unit-ideal
certificates over Q, hence valid over K̄ ⊇ Q̄).

## MAIN RESULTS

1. **The identity is t-free** (§1): on this stratum it collapses to a pure
   q-place statement F ≡ 0 of y-degree ≤ 40 — 41 slice equations in 22
   unknowns — with two equivalent finite cascades (q-adic and t¹⁰-adic).
2. **A new product structure** (§2): H₀ := Σ h_f|_{dm1=0} X^f factors as
   U₄·U₃, and the full identity becomes the **product-master**
   Ũ₄·Ũ₃ = −C²t²⁰·J (this is new structure available in every stratum's
   t¹⁰-blocks, but only (10,0) makes it global).
3. **The T2 branch (d̃1 ≡ 0) of stratum (10,0) is INFEASIBLE, char 0** (§3)
   — proven by a three-term t-adic Newton-polygon case tree on a completely
   factored master identity, plus six 3–4-variable unit-ideal certificates
   (sympy Gröbner over Q). This is the first full branch kill inside a
   surviving stratum of the T5_MULTIPLACE map.
4. **The T1 branch (d̃1 ≢ 0) is reduced to a minimal explicit system** (§4):
   the product-master + two cascades + three anchor/polygon conditions.
   NOT closed; two Singular F_p Gröbner probes (raw `std` and `slimgb`,
   30 min each) did not terminate — documented as evidence of difficulty,
   not of feasibility.

## 1. The t-free reduction and the two cascades

**Lemma 1.1 (t-power cancellation, [S1]).** For every f = 0..7,

    Φ̃^f (Ct¹⁰)^{21−3f} = c^f q^f C^{21−3f} · t²¹⁰

(the exponent bookkeeping 30f + 10(21−3f) = 210 is uniform in f). Hence

    Σ_f Φ̃^f e^{21−3f} h_f(d̃)  =  t²¹⁰ · F,
    F := Σ_f c^f C^{21−3f} q^f h_f(d̃2, d̃1, d̃0, Ct¹⁰),

and the window identity holds iff **F ≡ 0** — a t-free identity: deg F ≤ 40,
41 scalar slice equations. Equivalently, W = Φ̃/e³ = (c/C³)·q is a
**polynomial** root (degree 4, v_t = 0, v_q = 1) of H(d̃, W) = Σ h_f W^f.

**Lemma 1.2 (q-cascade, [S2]).** F ≡ 0 is equivalent to the finite cascade

    G_1 := h_0(d̃)/q                            (q | h_0(d̃) forced)
    G_{ℓ+1} := (C³G_ℓ + c^ℓ h_ℓ(d̃))/q,  ℓ = 1..6   (q | numerator forced)
    C³G_7 + c⁷·8192·d̃1² = 0                    (exact terminal)

with deg G_ℓ ≤ 40−4ℓ. (Forcing: reduce mod q as in T5_MULTIPLACE Lemma B;
sufficiency: telescoping, verified symbolically. Equation count
7·4 + 13 = 41 ✓.)

**Lemma 1.3 (t¹⁰-cascade, [S3]).** Split each h_f by dm1-layers:
h_f = Σ_k dm1^k h_{f,k}, and set H_k(d̃2,d̃1,d̃0; X) := Σ_f h_{f,k} X^f
(weight 20−5k with X of weight 2; window degree cap deg_y H_k(d̃, Bq)
≤ 40−10k, where B := c/C³). Then H₄ = −6561 (a nonzero constant!) and

    F = C²¹ [ H₀ + Ct¹⁰H₁ + C²t²⁰H₂ + C³t³⁰H₃ − 6561C⁴t⁴⁰ ]   (args: d̃, Bq),

equivalent to the cascade Q₁ := H₀/t¹⁰, Q₂ := (Q₁+CH₁)/t¹⁰,
Q₃ := (Q₂+C²H₂)/t¹⁰ (each division forced, deg Q_j ≤ 40−10j), with exact
terminal Q₃ + C³H₃ = 6561C⁴t¹⁰. Equation count 3·10 + 11 = 41 ✓.
(H₃ = 11664·d̃1·(4d̃2 − 5Bq) is tiny; see Lemma 2.1.)

*Anchors (remarks).* The y⁰ and y⁴⁰ slices of F ≡ 0 are the two anchor
hypersurface equations f31(τ2,τ1,τ0,C,−1/34) = 0 and
f31(λ2,λ1,λ0,C,−1024/3315) = 0 with the e-slot **pinned to C ≠ 0** in both
(cf. T5_MULTIPLACE §7); each alone kills nothing (irreducible hypersurfaces).

## 2. The product structure

**Lemma 2.1 (factorization of the 0-layer and the cofactor ladder, [S4]).**

    H₀ = U₄ · U₃,

    U₄ = 16X⁴ + 64X³d2 + 24X²d2² + 288X²d0 − 80Xd2³ + 576Xd0d2 − 432Xd1²
         + 25d2⁴ − 360d0d2² + 1296d0² + 216d1²d2                (weight 8)
    U₃ = 512X³d1² − 3072X²d0² + 1536X²d0d2² − 1152X²d1²d2 − 192X²d2⁴
         + 3072Xd0²d2 − 1152Xd0d1² − 1536Xd0d2³ + 1056Xd1²d2² + 192Xd2⁵
         + 6912d0³ − 4224d0²d2² + 1008d0d1²d2 + 816d0d2⁴ + 27d1⁴
         − 412d1²d2³ − 48d2⁶                                    (weight 12),

    H₁ = A₁·U₃ + B₁·U₄,   A₁ := −216·d1,
                          B₁ := 16·d1·(32X² − 56Xd2 − 36d0 + 29d2²),
    H₃ = 11664·d1·(4d2 − 5X),

and with D₂ := H₂ − A₁B₁ the **exact identity in all five variables**

    H = (U₄ + A₁·dm1)(U₃ + B₁·dm1) + dm1²·( D₂ + dm1·H₃ − 6561·dm1² ).

(H itself is irreducible — the product ansatz with dm1⁴-correction alone
does not close; D₂ ∉ ⟨U₄, U₃⟩ is the obstruction. All verified [S4].)

**Corollary 2.2 (product-master for stratum (10,0)).** With
Ũ₄ := U₄(d̃, Bq) − 216·Ct¹⁰·d̃1, Ũ₃ := U₃(d̃, Bq) + Ct¹⁰·B₁(d̃, Bq),
J := D₂(d̃, Bq) + Ct¹⁰·H₃(d̃, Bq) − 6561·C²t²⁰:

    **Ũ₄ · Ũ₃  =  − C² t²⁰ · J**,
    deg Ũ₄ ≤ 16,  deg Ũ₃ ≤ 24,  deg J ≤ 20.

First t-block: either v_t(Ũ₄) + v_t(Ũ₃) ≥ 20, or one factor ≡ 0 forcing
also J ≡ 0. This holds on ALL of (10,0) (both branches).

## 3. THEOREM — the T2 branch (d̃1 ≡ 0) of (10,0) is infeasible (char 0)

Throughout §3: d̃1 ≡ 0, σ := 4d̃0 − d̃2² (deg σ ≤ 8), and the normalized
variables  z := c·q (deg 4 exactly; t ∤ z, v_q(z) = 1),  v := C³d̃2,
S := C⁶σ,  K := C¹⁷ ≠ 0, and the six linear forms in (z,v)

    Â := z+v,  B̂' := 4z+v,  B̂ := 4z−5v,  L₄ := 10z+v,  L₅ := 2z−v,
    M := 29v−16z.

**Master identity ([S5], fully symbolic).** On T2 the stratum identity is
equivalent (C ≠ 0) to

    T₁ + T₂ + T₃ = 0,
    T₁ := 12·S²·N̂·P̂²,        N̂ := 9S − B̂B̂',   P̂ := 9S + 4Â²,
    T₂ := K·t²⁰·R̃,            R̃ := 512Â⁴B̂ − 432Â²·M·S + 2916B̂S²,
    T₃ := −6561·K²·t⁴⁰,

with deg S, N̂, P̂ ≤ 8, deg R̃ ≤ 20. (The identity comes from the exact
h_f-collapse h_f|_{d1=0, d0=(σ+d2²)/4} = σ²s_f + dm1²r_f + [f=0](−6561dm1⁴)
and the closed factorizations Σs_f X^f = 12(9σ−(4X−5d2)(4X+d2))(4(X+d2)²+9σ)²
and R̃ as displayed — all verified [S5][S6].)

Write m := v_t(S), n := v_t(N̂), p := v_t(P̂) (≤ 8 whenever nonzero),
α := v_t(Â) when Â ≢ 0. Two standing tools:

* **Newton lemma.** T₁+T₂+T₃ = 0 with T₃ ≠ 0 forces the minimum of the
  t-valuations of the (nonzero) terms to be attained at least twice;
  v(T₁) = 2m+n+2p (exact), v(T₂) = 20 + v_t(R̃), v(T₃) = 40 (exact).
  If T₂ = 0 (R̃ ≡ 0), then v(T₁) = 40 is forced.
* **Coprimality ([S7]).** The six forms are pairwise non-proportional, and
  any two combine Z-linearly to a nonzero multiple of z with t ∤ z: **at
  most one of the six forms can have v_t ≥ 1**, and no two can vanish
  identically simultaneously.

**Case tree** (exhaustive; σ ≡ 0 first, then by which form carries t):

* **σ-locus (S ≡ 0).** Master → 512Â⁴B̂ = 6561Kt²⁰. Valuation = 20 and
  degree ≤ 20 force Â = a·t⁴, B̂ = b·t⁴ ([S9]: 4d₁+d₂ = 20, d_i ≤ 4 only at
  (4,4)); then 9z = 5Â + B̂ = (5a+b)t⁴ makes q a multiple of t⁴ (q(−1)=3315≠0,
  ⨯) or z ≡ 0 (⨯). *Independent, elementary re-proof of T5_MULTIPLACE
  Theorem 2 restricted to this stratum.*
* **T2a: P̂ ≡ 0** (9S = −4Â²). Master → R̃ = 6561Kt²⁰ and
  R̃ = 128Â⁴L₄ [S6]. Same forced shape Â = at⁴, L₄ = bt⁴; L₄ − Â = 9z ⨯.
* **T2b: N̂ ≡ 0** (9S = B̂B̂'). R̃ = 4B̂²(8z−v)L₄² [S6] = 6561Kt²⁰ forces
  B̂, 8z−v, L₄ all = const·t⁴; L₄ + (8z−v) = 18z ⨯.
* **I₀: Â ≡ 0** (v = −z). Master/2187 = 4S⁵ − 12z²S⁴ + 12zKt²⁰S² − 3K²t⁴⁰
  [S8]. m = 0 forces S = 3z² exactly, then 36z⁵ = Kt²⁰ (v_t: 0 ≠ 20 ⨯);
  m ≥ 1 forces 4m ≥ 20, and {4m, 20+2m, 40} are pairwise distinct for
  m ≤ 8 ⨯ [S8].
* **B̂ ≡ 0 or B̂' ≡ 0** (v = 4z/5 or −4z). Then N̂ = 9S, n = m, Â is a
  t-unit so p = 0 for m ≥ 1: v(T₁) = 3m with v(T₂) ∈ {20+m, 20} exact —
  3m ∉ {20, 20+m, 40} and 20+m ≠ 40 for m ≤ 8 ⨯ [S9]; m = 0 gives
  v(T₁) ≤ 16 < 20 ⨯ [S7].
* **L₄ ≡ 0 / L₅ ≡ 0** (v = −10z / v = 2z). Then P̂ = N̂; m = 0 forced and
  3p ≥ 20 → p ∈ {7,8}: closed cells W = −9z (resp. 3z),
  9S + 4W² = t⁷(v₀+v₁y). **Unit-ideal certificates** [S10] ⨯.
* **Case D: all six forms t-units — or only M carries t, or M ≡ 0.** The
  form M enters only the middle R̃-group, so these three configurations
  share one argument: n ≥ 1 and p ≥ 1 impossible (v_t(P̂−N̂) = v_t(L₄L₅)
  = 0; and P̂ = 9S + 4Â² has p = 0 whenever m ≥ 1 since Â is a t-unit);
  every profile leaves v(T₁) ≤ 16 < 20, minimum unique ⨯ [S7].
* **Case B: v_t(B̂) = β ∈ [1,4].** p = 0 for m ≥ 1 (Â unit) and the
  R̃-groups have exact valuations (β, m+0, β+2m) → for m ≠ β:
  v(T₁) = 2m+min(m,β), v(T₂) = 20+min(m,β): equal only at m = 10 ⨯;
  m = β ≤ 4: v(T₁) ≤ 16 < 20 ⨯; m = 0: v(T₁) ≤ 16 ⨯ [S9].
* **Case B': v_t(B̂') = β' ∈ [1,4].** Now the group 512Â⁴B̂ is a t-unit:
  v(T₂) = 20 exactly (m ≥ 1). The polygon forces 2m+n = 20 with
  n = min(m,β') ≤ 4: **only (m, β') = (8,4)** survives [S9] — the closed
  cell S = γt⁸, B̂' = δt⁴. **Unit-ideal certificate** [S10] ⨯.
* **Case C / C': v_t(L₄) (or L₅) = λ ∈ [1,4].** m = 0 forced; the profile
  n+2p ≥ 20, n,p ≤ 8, min(n,p) ≤ λ ≤ 4 has the **unique solution
  (n,p) = (4,8)** [S9]: closed cells L₄ = δt⁴ (resp. L₅ = δt⁴) and
  P̂ = εt⁸. **Unit-ideal certificates** [S10] ⨯.
* **Case A: v_t(Â) = α ∈ [1,4]** (the fat case). Here n = 0 for m ≥ 1;
  m = 0 dies at the first block; so m ≥ 1 and the Newton lemma runs on
  exact valuations ([S9] throughout):
  - **m < 2α:** p = m and R̃-groups (4α, 2α+m, 2m) are distinct, so
    v(T₂) = 20+2m; v(T₁) = 4m < min(20+2m, 40) — minimum unique ⨯.
  - **m > 2α:** p = 2α, v(T₂) = 20+4α; v(T₁) = 2m+4α < min(20+4α, 40) ⨯.
  - **m = 2α:** write W = t^αW₁, S = t^{2α}S₁ (t ∤ W₁S₁). The regrouping

        R̃ = t⁴ᵅ·( 9z·Θ − t^α·W₁·Θ' ),
        Θ := 512W₁⁴ + 2160W₁²S₁ + 2916S₁²,
        Θ' := 2560W₁⁴ + 12528W₁²S₁ + 14580S₁²,

    and the **incompatibility lemma** Θ|_{9S₁=−4W₁²} = 128W₁⁴ ≠ 0 show
    that the two possible cancellations κ₁ := v_t(9S₁+4W₁²) ≥ 1 and
    κ₂ := v_t(9zΘ − t^αW₁Θ') ≥ 1 are mutually exclusive. Then:
    κ₁ = 0 → v(T₁) = 8α < min(20+4α, 40) ⨯; κ₁ ≥ 1 → v(T₂) = 20+4α
    exactly, and with the cap p = 2α+κ₁ ≤ 8 (P̂ ≢ 0):
    v(T₁) = 4α+2p < 20+4α = v(T₂) and v(T₁) ≤ 32 < 40 — minimum unique ⨯.

All cases are exhausted: after σ ≡ 0, T2a, T2b, any remaining tuple has
S, N̂, P̂ ≢ 0; the coprimality lemma says at most one of the six forms
vanishes identically or carries positive t-valuation (two such forms force
z ≡ 0), which is precisely the label list {I₀, B̂≡0, B̂'≡0, L₄≡0, L₅≡0,
M≡0, D, B, B', C, C', A, only-M-carries-t (in D)}. ∎

**Status of Theorem 3: PROVEN (char 0).** Machine ingredients: only the
symbolic identities [S1–S9] and six unit-ideal Gröbner certificates over Q
in ≤ 4 variables [S10], all verified inside `t5_stratum100_verify.py` by
sympy alone (no Singular, no F_p anywhere in the T2 proof).

**Corollary.** Any window solution in stratum (10,0) has **d̃1 ≢ 0** (the
T1 branch), and then the q-cascade terminal of Lemma 1.2 gives
G_7 = −8192(c⁷/C³)·d̃1² ≢ 0: the whole q-cascade is live.

## 4. The T1 branch (d̃1 ≢ 0): minimal explicit system — OPEN

What any T1 solution must satisfy (all PROVEN reductions):

1. **Product-master (Corollary 2.2):** Ũ₄Ũ₃ = −C²t²⁰J. Either
   v_t(Ũ₄) + v_t(Ũ₃) = 20 + v_t(J) with v_t ≤ deg caps (16, 24, 20), or a
   degenerate sub-branch:
   - Ũ₄ ≡ 0 ∧ J ≡ 0:  U₄(d̃,Bq) = 216Ct¹⁰d̃1  and
     D₂(d̃,Bq) + Ct¹⁰H₃(d̃,Bq) = 6561C²t²⁰;
   - Ũ₃ ≡ 0 ∧ J ≡ 0:  U₃(d̃,Bq) = −Ct¹⁰B₁(d̃,Bq)  and the same J-identity.
2. **q-cascade** (Lemma 1.2) with live terminal C³G₇ = −8192c⁷d̃1².
3. **t¹⁰-cascade** (Lemma 1.3) with terminal Q₃ + C³H₃(d̃,Bq) = 6561C⁴t¹⁰,
   where H₃(d̃,Bq) = 11664·d̃1·(4d̃2 − 5Bq) — note both terminal objects
   are degree ≤ 10 and involve d̃1 linearly/quadratically.
4. **Three-place polygon conditions** (T5_NP Lemma 2 and T5_MULTIPLACE §6–7
   specialized to w = Bq): at t (v_t(w) = 0): min_f v_t(h_f(d̃)) attained
   twice; at q (v_q(w) = 1): min_f (v_q(h_f(d̃)) + f) attained twice
   (f = 0 term: q | h_0(d̃), the cascade start); at ∞ (v_∞(w) = −4):
   max_f (deg h_f(d̃) + 4f) attained twice.
5. **Anchors**: f31(τ2,τ1,τ0,C,−1/34) = 0 and f31(λ2,λ1,λ0,C,−1024/3315) = 0
   with shared C ≠ 0.

**Computational attempts (F_32003 probes — NOT proofs, and none
terminated):** (i) the raw 41-equation system in 23 variables (window
coefficients, C, B, with 6630·B·C³+1 = 0), Singular 4.2.1 `std` and
`slimgb`, dp order: killed at the 1800 s timeout, for the full system and
for the (now superseded) T2 subsystem; (ii) the degenerate sub-branch
Ũ₄ ≡ 0 ∧ J ≡ 0 (38 lower-degree equations, same 23 variables): `std`
killed at 1700 s with memory past 3.7 GB; the Ũ₃ ≡ 0 ∧ J ≡ 0 probe showed
the same divergence pattern and was cancelled. Consistent with the
T5_NOTES.md experience that raw elimination does not close these systems;
the T2 kill above (structure first, tiny certificates last) is the
template to replicate on T1.

**Why T1 is genuinely harder:** the T2 kill lived off the complete
factorization of the master into the six linear forms in (z,v) — one
polynomial unknown v (plus S). On T1 the analogous objects U₄, U₃, D₂ are
irreducible in X with three polynomial unknowns, so the t-valuation case
tree does not reduce to linear-form bookkeeping. The most promising next
moves, in order:
1. run the three-term polygon on the 5-block form
   Ũ₄Ũ₃ + C²t²⁰D₂ + C³t³⁰H₃ − 6561C⁴t⁴⁰ = 0 (four terms, valuations
   V₁ := v(Ũ₄)+v(Ũ₃), 20+v(D₂(d̃,Bq)), 30+v_t(d̃1)+v_t(4d̃2−5Bq), 40 — the
   H₃ term is a product of LINEAR objects, and v(T₃-term) interacts with
   v_t(d̃1) which also caps h_7 = 8192d̃1² in the t-polygon of item 4);
2. the q-side: J ≡ D₂(d̃, 0) mod q and Lemma-B-style q-multiplicity
   pumping down the q-cascade combined with the product-master;
3. only then a structured Gröbner (cells, not raw).

## 5. Verification map

`t5_stratum100_verify.py` (sympy over Q, all exact; runs in a few minutes):

| check | verifies | used in |
|---|---|---|
| S1 | graded decomposition; t-power bookkeeping; t-free F | §1 |
| S2 | q-cascade telescoping, deg G_ℓ ≤ 40−4ℓ | §1 |
| S3 | dm1-layers H_k, weights/caps, H₄ = −6561; t¹⁰-regrouping + cascade | §1 |
| S4 | H₀ = U₄U₃; H₁ = A₁U₃+B₁U₄; H₃; the 5-variable product identity | §2 |
| S5 | T2 h_f-collapse; S-factorization; R̃ closed form; symbolic master | §3 |
| S6 | R̃ collapses on P̂≡0 / N̂≡0 / S≡0; q(−1) ≠ 0 | §3 |
| S7 | six forms pairwise t-coprime; case D / m=0 arithmetic | §3 |
| S8 | I₀ identity + three-term kill arithmetic | §3 |
| S9 | case-A regrouping, Θ-lemma, all polygon pairing arithmetic; forced rigid profiles for B'/C/C'; σ-locus (4,4) forcing | §3 |
| S10 | six unit-ideal certificates (B, B', C, C', L₄≡0, L₅≡0) over Q | §3 |

External inputs: f31_deg31.txt, f31_graded.txt (both pre-existing, re-verified
in S1), and T5_MULTIPLACE.md only for context (its Theorem 2 is *re-proved*
here for this stratum, in the σ-locus bullet, so §3 is self-contained).

## 6. Honest status

* **(10,0) ∩ {d̃1 ≡ 0}: infeasible — PROVEN, char 0** (§3). No F_p input,
  no unverified machine step; the only computer algebra beyond expansion is
  six ≤ 4-variable Gröbner unit-ideal certificates over Q inside sympy.
* **(10,0) ∩ {d̃1 ≢ 0}: OPEN.** Reduced to the explicit system of §4
  (product-master + two cascades + polygon/anchor conditions). Raw F_p
  Gröbner probes did not terminate (so not even probe-level evidence from
  that route); the numeric jetlift evidence for f31·sub2 as a whole (T2/T3
  campaigns) covers this stratum but is not a proof.
* Nothing in this file is conjectural; the two labels above are exact.
