# T5_STRATA_50_11 — strata (5,0) and (1,1) are INFEASIBLE (2026-07-21)

> **STATUS (2026-07-22):** CONDITIONAL / SCOPE-RESTRICTED. The `(a_t, a_q)` survivor maps below are valid ONLY as geometrically-q-coprime (uniform-qʳ) statements: after base change q splits into four geometric places and the scalar `a_q` ledger is not field-stable (`FIELD_SPLIT_AUDIT.md`; `STATE.md` L5–16). The authoritative replacements are the split-place ledger (`split_place_ledger.py` / `split_place_ledger.json`) and the cascade-engine artifacts (`cascade_cones_qt_inf_rl.json`).

Resolves open item 1 of T5_MULTIPLACE.md §8. **Both target strata are dead**:

    KILLED THIS SESSION:  (a_t, a_q) ∈ { (5,0), (1,1) } .

Everything below is **PROVEN** (nothing conjectural remains in this file). All
computational inputs are verified exactly by `t5_strata50_11_verify.py`
(8 check groups W1–W8, all pass); citations [Wn] refer to its checks, [Vn] to
`t5_multiplace_verify.py`. Setting: f31 branch, subcase (2), stripped windows
deg d̃2 ≤ 4, d̃1 ≤ 6, d̃0 ≤ 8, e := d̃₋₁ (deg ≤ 10, e ≢ 0), t := y+1,
q := 2048y⁴−512y³+320y²−240y+195, c := −1/6630, Φ̃ = c t³⁰ q, u := c q,
σ := 4d̃0 − d̃2², stratum (a, a_q) := (v_t(e), v_q(e)), and the identity under
attack is

    0 = Σ_{f=0}^{7} Φ̃^f e^{21−3f} h_f(d̃) .                              [W1]

## 0. Inherited inputs (all proven in T5_NP.md / T5_MULTIPLACE.md)

* **Lemma A** (T5_MULTIPLACE §1): for a ≤ 9, writing e = t^a ê (t ∤ ê,
  deg ê ≤ 10−a), v := 30−3a, the identity is equivalent to the cascade
  g_1 := h_0(d̃)/t^v, g_{ℓ+1} := (ê³g_ℓ + u^ℓ h_ℓ(d̃))/t^v (ℓ = 1…6),
  ê³g_7 + u⁷h_7(d̃) = 0, with all g_ℓ ∈ K[y] of degree ≤ 10+3a.
  (Telescoping re-verified [W5].)
* **Proposition 1** (T5_MULTIPLACE §3): if 4(7−3a_q) > 10+3a then d̃1 ≡ 0
  and g_7 = 0. Holds in both target strata: 28 > 25 and 16 > 13 [W4].
* **Theorem 2** (T5_MULTIPLACE §4): the σ-locus {d̃1 ≡ 0, σ ≡ 0, e ≢ 0} is
  infeasible in **every** joint stratum. (Used twice below, as a black box;
  its own computational inputs are [V3][V4][V8].)
* Collapse identities [W2]: h_7 = 8192 d̃1², h_6|_{d̃1≡0} = −3072 σ²,
  **h_5|_{d̃1≡0} = −9216 d̃2 σ² + 2048 e²**.

## 1. Common reduction: the terminal pair (ê, ĝ)

Both strata have a single open branch (T2) with d̃1 ≡ 0 and g_7 = 0 forced
(Prop. 1). The ℓ = 6 cascade line then reads 0 = ê³g_6 + u⁶h_6(d̃) =
ê³g_6 − 3072c⁶q⁶σ², i.e.

    ê³ g_6 = 3072 c⁶ q⁶ σ² ,

and σ ≢ 0 (else the tuple lies on the σ-locus, dead by Theorem 2), so
g_6 ≠ 0.

**Stratum (5,0)** (v = 15, deg g_ℓ ≤ 25, q ∤ ê, deg ê ≤ 5): taking v_q,
v_q(g_6) = 6 + 2v_q(σ). If v_q(σ) ≥ 1 then 4·v_q(g_6) ≥ 32 > 25 ≥ deg g_6,
impossible; so **v_q(σ) = 0**, g_6 = q⁶ĝ with q ∤ ĝ and
deg ĝ ≤ 25 − 24 = 1 [W4]. Dividing by q⁶:

    (5,0):   ê³ ĝ = 3072 c⁶ σ² ,    deg ĝ ≤ 1 .

**Stratum (1,1)** (v = 27, deg g_ℓ ≤ 13, ê = q ē with q ∤ ē, t ∤ ē,
deg ē ≤ 5): the line gives q³ē³g_6 = 3072c⁶q⁶σ², so ē³g_6 = 3072c⁶q³σ²,
v_q(g_6) = 3 + 2v_q(σ); v_q(σ) ≥ 1 would force 4·v_q(g_6) ≥ 20 > 13,
impossible; so **v_q(σ) = 0**, g_6 = q³ĝ, deg ĝ ≤ 13 − 12 = 1 [W4], and

    (1,1):   ē³ ĝ = 3072 c⁶ σ² ,    deg ĝ ≤ 1 .

In both strata ĝ ≠ 0: ĝ = 0 would give σ² = 0, hence σ ≡ 0 (char-0 domain),
dead by Theorem 2.

## 2. Lemma S (the level-5 squeeze): ê resp. ē is a nonzero constant

This is the new mechanism. The ℓ = 5 cascade line is
t^v g_6 = ê³g_5 + u⁵h_5(d̃), and with d̃1 ≡ 0,
h_5(d̃) = −9216 d̃2 σ² + 2048 e² [W2]. The point: the terminal relation of §1
lets the σ² inside h_5 be **absorbed into the ê³-divisible part**, leaving a
two-term remainder whose divisibility by ê³ is fatal. The needed constant is
−9216c⁵/(3072c⁶) = −3/c = **19890** [W5].

**Lemma S1.** In stratum (5,0), ê ∈ K^×.

*Proof.* e² = t^{10}ê² and g_6 = q⁶ĝ, so the ℓ = 5 line reads

    t¹⁵ q⁶ ĝ = ê³ g_5 + c⁵q⁵( −9216 d̃2 σ² + 2048 t¹⁰ ê² ) .

Substituting σ² = ê³ĝ/(3072c⁶) (terminal relation, §1) and rearranging — the
exact polynomial identity is verified symbolically [W5] —

    q⁵ t¹⁰ ( t⁵ q ĝ − 2048 c⁵ ê² )  =  ê³ · G₅ ,
    G₅ := g_5 + 19890 q⁵ d̃2 ĝ  ∈ K[y] .

Since gcd(ê, t) = gcd(ê, q) = 1, the UFD K[y] gives
ê³ | t⁵qĝ − 2048c⁵ê²; in particular ê² divides that difference, and
ê² | 2048c⁵ê², so ê² | t⁵qĝ, and coprimality again gives **ê² | ĝ**. If ê
were non-constant then deg ê² ≥ 2 > 1 ≥ deg ĝ would force ĝ = 0,
contradicting §1. So ê is a nonzero constant. ∎

**Lemma S2.** In stratum (1,1), ē ∈ K^×.

*Proof.* Identical squeeze: e² = t²q²ē², g_6 = q³ĝ, and the ℓ = 5 line

    t²⁷ q³ ĝ = q³ē³ g_5 + c⁵q⁵( −9216 d̃2 σ² + 2048 t² q² ē² )

with σ² = ē³ĝ/(3072c⁶) rearranges exactly [W5] to

    q³ [ t² ( t²⁵ ĝ − 2048 c⁵ q⁴ ē² )  −  ē³ ( g_5 + 19890 q² d̃2 ĝ ) ] = 0 .

Cancelling q³ (domain): t²(t²⁵ĝ − 2048c⁵q⁴ē²) = ē³G₅′. Since gcd(ē, t) = 1,
ē³ | t²⁵ĝ − 2048c⁵q⁴ē², hence ē² | t²⁵ĝ, hence ē² | ĝ (gcd(ē, t) = 1), and
as before deg ĝ ≤ 1 forces ē ∈ K^×. ∎

**Corollary S3 (rigidification).** In stratum (5,0): e = C t⁵; in stratum
(1,1): e = C t q — with C ∈ K^× — and in both, the terminal relation becomes
ĝ = (3072c⁶/C³)σ² with deg ĝ ≤ 1; since deg σ² is even, **σ is a nonzero
constant s**, i.e.

    d̃1 ≡ 0 ,    d̃0 = (d̃2² + s)/4 ,    s ∈ K^× ,    deg e = 5 .

Only d̃2 (deg ≤ 4), C, s remain free.

## 3. The endgame: the place at infinity kills the rigid pencil

Write H_f := h_f(d̃2, 0, (d̃2²+s)/4, e) for the h_f with d̃1 ≡ 0 and σ = s
substituted. Structure of H_f as a polynomial in (d̃2, s, e) [W3]: every
monomial d̃2^k s^m e^j has weight 2k + 4m + 5j = 20 − 2f, j is even (j = 4
only at f = 0), and — crucially, by the σ-locus collapse [V3]/[W3] — **j = 0
forces m ≥ 1** (no pure-d̃2 monomials). Hence, with deg d̃2 ≤ 4, deg e = 5,
s constant:

    deg_y H_f ≤ 32 − 4f  (f ≤ 5),    H_6 = −3072 s²  (degree 0!),   H_7 = 0 .

**Proposition E (degree domination).** No tuple with d̃1 ≡ 0, σ = s ∈ K^×,
and e of degree 5 of the forced shapes above satisfies the identity.

*Proof.* Under the substitutions the identity reads
0 = P := Σ_{f=0}^{6} T_f, T_f := Φ̃^f e^{21−3f} H_f. Degrees:

    deg T_f = 34f + 5(21−3f) + deg H_f ≤ 105 + 19f + (32−4f) = 137 + 15f
            ≤ 212   for f ≤ 5,

(verified with symbolic coefficients d̃2 = a₄y⁴+…+a₀ — the formal generic
degree bounds every specialization; the actual values are 135 + 15f ≤ 210
[W6][W7]), while the f = 6 term is **exactly**

    (5,0):  T_6 = Φ̃⁶ (Ct⁵)³ (−3072s²) = −3072 s² c⁶ C³ t¹⁹⁵ q⁶ ,
    (1,1):  T_6 = Φ̃⁶ (Ctq)³ (−3072s²) = −3072 s² c⁶ C³ t¹⁸³ q⁹ ,

of degree **219** in both cases, with leading coefficient
−3072·(1024/3315)⁶·s²C³ resp. −3072·(2048⁹/6630⁶)·s²C³ [W6][W7]. Hence the
y²¹⁹-coefficient of P equals that leading coefficient, a nonzero constant
multiple of s²C³ ≠ 0 (s, C ∈ K^×). So P ≠ 0 — contradiction. ∎

(Conceptually: this is the ∞-place Newton polygon of T5_MULTIPLACE §7. Once
deg e drops to 5, v_∞(w) = 15 − 34 = −19, the per-f increment 19 beats the
cap decay 4, and the max of deg H_f + 19f is attained **once**, at f = 6 —
where H_6 = −3072s² cannot vanish. Numeric instance checks: [W8].)

## 4. Main results

**Theorem 3.** The joint stratum **(a_t, a_q) = (5,0) is infeasible**: no
window tuple with d̃ in the subcase-(2) windows, e ≢ 0, v_t(e) = 5,
v_q(e) = 0 satisfies the f31 identity.

*Proof.* Prop. 1 forces d̃1 ≡ 0, g_7 = 0; §1 forces v_q(σ) = 0 and the
terminal pair ê³ĝ = 3072c⁶σ², deg ĝ ≤ 1, ĝ ≠ 0 (Theorem 2); Lemma S1 forces
ê = C ∈ K^×; Corollary S3 forces σ = s ∈ K^×, e = Ct⁵, d̃0 = (d̃2²+s)/4;
Proposition E shows the resulting identity has a nonvanishing y²¹⁹
coefficient. ∎

**Theorem 4.** The joint stratum **(a_t, a_q) = (1,1) is infeasible**.

*Proof.* Same chain with Lemma S2: d̃1 ≡ 0, v_q(σ) = 0, ē³ĝ = 3072c⁶σ²,
ē = C ∈ K^×, e = Ctq (degree 5), σ = s ∈ K^×, and Proposition E. ∎

**Updated survivor map** (was 15 in T5_MULTIPLACE §5; now **13**):

    a_q = 0 : a ∈ {6, 7, 8, 9}  and the degenerate a = 10
    a_q = 1 : a ∈ {2, 3, 4, 5, 6}
    a_q = 2 : a ∈ {0, 1, 2}

Every remaining stratum with a ≤ 9 has the T1 branch (d̃1 ≢ 0) open — the
"single-branch" rigidity that made (5,0)/(1,1) tractable is gone, but the
two new tools transport:

* **The level-5 squeeze** (Lemma S mechanism) applies verbatim to every T2
  branch (d̃1 ≡ 0, σ ≢ 0) of every surviving stratum with a ≤ 9: it always
  yields ê² | ĝ-type divisibilities; only the cap deg ĝ ≤ (10+3a) − 4(6−3a_q)
  grows with a (e.g. (6,0): deg ĝ ≤ 4, so it forces deg ê ≤ 2 instead of
  ê constant — still a strong reduction).
* **Proposition E generalizes**: for any stratum reduced to d̃1 ≡ 0,
  σ = s ∈ K^× and deg e = D, the f = 6 term has degree 204 + 3D and term f
  is capped by 34f + (21−3f)D + max(2D+20−4f, 32−4f, 4D·[f=0]); the f = 6
  term dominates strictly whenever (6−f)(34−3D) > deg-cap gap — in
  particular for every D ≤ 9 one checks 34 − 3D ≥ 7 beats the caps (not
  verified in-script beyond D = 5; the D = 5 instances are what Theorems 3–4
  use).

## 5. Status ledger

**PROVEN** (this file): §1 reductions, Lemmas S1–S2, Corollary S3,
Proposition E, Theorems 3–4, the 13-stratum survivor map.

**Dependencies** (proven elsewhere): Lemma A, Proposition 1, Theorem 2 of
T5_MULTIPLACE.md (Theorem 2 uses Mason–Stothers); the graded decomposition
[V2]/[W1]. Non-computational ingredients here: UFD/valuation arithmetic in
K[y] and degree bookkeeping only — no new Mason–Stothers use was needed.

**CONJECTURED** (not used anywhere above): the D ≤ 9 generalization of
Proposition E stated in §4's second bullet.

## 6. Verification map

`t5_strata50_11_verify.py` (all exact, sympy over Q; run from this
directory; ~1 min):

| check | verifies | used in |
|---|---|---|
| W1 | setup: q irreducible, Φ̃ = ct³⁰q, v_t = 30, v_q = 1, deg 34, lc −1024/3315; graded decomposition + weights | §0 |
| W2 | h_7 = 8192d̃1²; h_6\|_{d̃1=0} = −3072σ²; h_5\|_{d̃1=0} = −9216d̃2σ² + 2048e² | §0, §1, §2 |
| W3 | H_f(d̃2, s, e) monomial structure: weights, even e-powers, e⁴ only f=0, no pure-d̃2 terms; deg H_f ≤ 32−4f; H_6 = −3072s²; H_7 = 0; s⁰-parts = σ-locus constants | §3 |
| W4 | stratum arithmetic: C1 margins 28 > 25, 16 > 13; v_q(σ)=0 squeezes 32 > 25, 20 > 13; deg ĝ ≤ 1 both | §1 |
| W5 | Lemma A telescoping; −9216c⁵/(3072c⁶) = 19890; the two exact level-5 rearrangement identities | §2 |
| W6 | (5,0) endgame: T_6 = −3072s²c⁶C³t¹⁹⁵q⁶ (deg 219, lc −3072(1024/3315)⁶s²C³); symbolic per-term degrees 135+15f ≤ 210 | §3, Thm 3 |
| W7 | (1,1) endgame: T_6 = −3072s²c⁶C³t¹⁸³q⁹ (deg 219, lc −3072·2048⁹/6630⁶·s²C³); per-term degrees ≤ 210 | §3, Thm 4 |
| W8 | numeric instances (both strata): P has degree exactly 219, predicted lc, ≠ 0; agrees with direct graded evaluation | §3 |
