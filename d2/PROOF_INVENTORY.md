# PROOF INVENTORY — the (72,108) program, single source of truth

> **HISTORICAL CLAIM GRAPH.** This file is the historical, judgment-tagged claim
> graph (every claim the writeup might rest on, including superseded, conditional,
> and evidence-only nodes). For the **current** state — proven and audited nodes,
> trusted inputs, and the live 224-cell frontier — see **`CURRENT_STATUS.md`**,
> which supersedes this file for current truth.

**Generated:** 2026-07-22, from git HEAD `51ac3f0` ("Sub1 Phase C worklist"), plus the
uncommitted working-tree docs `ALT_REGIME_L2.md`, `T5_T2_INFINITY.md` and their checkers
(both already logged in `STATE.md` and wired into `run_tests.sh`).
**Scope:** every claim in `d2_plane_72_108` that the eventual writeup will rest on: what is
claimed, where it is proven, which script checks it, whether it has been independently
audited, and what depends on it. Anything inferred rather than read is marked `[inferred]`.
**Convention:** checkers marked ♦ are wired into `run_tests.sh` (the authoritative suite);
checkers marked ◇ exist in the repo but are NOT in `run_tests.sh` (see Inconsistency I8).

---

## 0. Trust-tier legend

| Tier | Meaning |
|---|---|
| **1** | Independently audited: a separately authored checker/derivation with no code sharing with the artifact it audits (the Codex-authored `audit_cascade_kills.py`, `audit_cascade_kills_sub1.py`, `audit_tplace_cases.py`; the externally authored `FIELD_SPLIT_AUDIT.md` derivations cross-checked by a local exact checker). |
| **2** | Exact checker in `run_tests.sh` (♦). Same-author verification, exact arithmetic, re-run before commit. |
| **2\*** | Exact checker exists in the repo (◇) but is NOT wired into `run_tests.sh`. Substantively tier-2 strength, procedurally weaker. |
| **3** | Document-only, judgment-tagged, conditional, or evidence-only (numerics). No exact checker enforces the claim, or the claim is explicitly non-probative. |
| **4** | Published result used as stated (not re-proven here). Exact list in §4.1. |

A row may carry a compound tier (e.g. `2/4`) when a checker verifies the finite arithmetic
but the claim also leans on a published proposition used as stated.

---

## 1. The claim graph, top down

Dependencies point upward: each node depends on the nodes listed in "Deps".
The target theorem is C0; everything else supports it.

### 1.1 Target and setup

| ID | Claim (one line) | Proven in | Checker | Audit status | Tier | Deps |
|---|---|---|---|---|---|---|
| **C0** | Target theorem (OPEN): no P,Q ∈ K[x,y] with [P,Q]=x² and the Prop-4.3 case-(8,28) Newton polygons (subcases (1),(2)); if closed, the plane JC counterexample bound rises 108→125. | STATE.md L18–23 (statement); not yet proven — see §2 frontier | — | open | — | C1–C33 |
| **C1** | The (72,108)/(8,28) case is the paper's own last open case below 125; both subcase polygons transcribed corner-for-corner verbatim from the published Prop 4.3. | T3_WINDOW_AUDIT.md §1 L17–36; AUDIT.md §B L55–68 | none (LaTeX-source comparison, in-session) | independently re-checked against arXiv source; document-only | 3/4 | GGHV22 Prop 4.3 |
| **C2** | Normalization: ℓ₁,₀(P)=R², ℓ₁,₀(Q)=R³, R=x⁴C₄, C₄=y⁷(y+1) (Premise 1). | T6_PREMISES.md §1 L48–128 (READY-WITH-CITATION) | t6_premises_verify.py ♦ (checks P1a–P1e) | checker + published props as stated | 2/4 | GGV1 Props 1.13, 2.1 |
| **C3** | α-strip WLOG: Q=C³+λC⁻¹+F with v₁,₀(F)=−5 (Premise 2), transported from the paper's (9,27) template t=3→t=4. | T6_PREMISES.md §2 L140–202 | t6_premises_verify.py ♦ (P2a–P2e) | checker + published argument as stated | 2/4 | C2; GGHV22 §4 L1508–1546 |
| **C4** | Forcing ODE: [P,Q]=x² ⇒ 8y(y+1)f₁′−14(8y+7)f₁=y⁸(y+1)², unique solution f₁=−y⁸(y+1)²·q/6630 with q the separable quartic. | STATE.md item 2 L28–31; AUDIT.md A.4 L43–45; T6_SELECTION_AUDIT.md §1 L28–36 | verify_derivation.py ◇ (§C of its 48 checks) | checker-only (not in suite) | 2\* | C2, C3 |
| **C5** | D-transformation: D_k=C_k·C₄^(7−2k) ∈ K[y]; recursion closes with exact C₄-exponent cancellation. | STATE.md item 3 L32–35; T6_SELECTION_AUDIT.md §1 L41–44 | verify_derivation.py ◇ | checker-only (not in suite) | 2\* | C2 |
| **C6** | Envelope/window bounds: deg d_{4−k} ≤ 14w (sub2) / 15w (sub1), ord ≥ 12w; proven by redoing the paper's valuation induction (v₂,₋₁, v₋₂,₁, v₋₁,₁) on our polygons; Φ attains both bounds (tight). | T3_WINDOW_AUDIT.md §3 L69–94; AUDIT.md Risks 1–2 resolved L100–104 | none committed (in-session sympy; jetlift CONFIGS consistency only) | document-only — **no committed checker; writeup gap G2** | 3 | C1, C2; template GGHV22 Props 5.2(3)(4)/5.6 |
| **C7** | Equation selection: the 12-equation window system is the correct necessary condition; k=8 vacuous, λ isolated in (D̃³)₋₄, dropped equations define only spare unknowns/λ; slice bridge equals regenerate_system.py exactly. | T6_SELECTION_AUDIT.md L10–68; AUDIT.md C L82–89 | verify_derivation.py ◇ (48 symbolic checks, all pass) | checker-only (not in suite) | 2\* | C4, C5 |
| **C8** | System regeneration is reproducible and validated against the paper's published (7,21) case verbatim; f31 (102 terms) and f37 (618 terms) regenerate byte-identical. | STATE.md T1 session L134–148; AUDIT.md A.1 L22–27 | regenerate_system.py + run_singular.sh ◇ (self-validating) | fresh-container reproduction (same program) | 2\* | C7 |

### 1.2 Elimination and the f31-only reduction

| ID | Claim | Proven in | Checker | Audit status | Tier | Deps |
|---|---|---|---|---|---|---|
| **C9** | Master identity: f31·f37·d₋₁²¹ ≡ 0 in K[y] is a sound necessary condition (A=dm1·Ah, B=dm1·Bh with unique dm2-factors; resultant chain has no gap). | AUDIT.md A.2 L29–35; STATE.md item 5 L43–48 | Singular via run_singular.sh ◇ | checker-only (Singular-trusting); largely superseded in load-bearing role by C11 | 3 | C7, C8 |
| **C10** | d₋₁ ≡ 0 is impossible (both legs force Φ=0, contradiction). | STATE.md item 5 L49–55; AUDIT.md A.3 L36–41 | re-verified symbolically in the T6 session `[inferred: within verify_derivation.py ◇]` | checker-only | 2\* | C7 |
| **C11** | **f31-only reduction (ideal membership):** f31 ∈ ⟨G1,G2,G3,G5body+Φ⟩ over Q with Φ free; the elimination ideal in (d̃2,d̃1,d̃0,d₋₁,Φ) is EXACTLY ⟨f31⟩; f37 and d₋₁²¹ are classical resultant excess — the entire f37 branch is an artifact, uniformly, for BOTH subcases. | F37_SATURATION_REPORT.md L8–23, §1 table L42–49, §5 scope L127–142 | f37_sat_verify.py ♦ (sympy re-expansion of the lift certificate, "not trusting the Gröbner engine"); independent H-system route L57–62; f37_sat_confirm.sing ◇ | exact certificate, cross-derived two ways, same author | 2 | C7, C8, C9 |
| **C12** | f37 free family d̃2=d̃1=0 does not lift to the original system (field-stable proof; now a subsumed special case of C11). | F37_FREE_FAMILY_SYSTEM.md L6–169 | f37_free_family_verify.py ♦ (checks against regenerated state, not embedded resultant data) | checker; subsumed by C11 | 2 | C7 (subsumed by C11) |
| **C13** | Graded decomposition: f31 = Σ_{f=0..7} Φ̃^f d̃₋₁^{21−3f} h_f with explicit small h_f; the unit cofactor of Φ̃ is −q/6630, so the cascade runs at t and at the q-places simultaneously. | STATE.md L275–292; T5_NP.md; T5_MULTIPLACE.md §1–2 | verify_graded.py ◇; t5_multiplace_verify.py ◇ (8 exact check groups) | checker-only (not in suite) | 2\* | C11 |

### 1.3 The field-split framework (subcase 2 and shared machinery)

| ID | Claim | Proven in | Checker | Audit status | Tier | Deps |
|---|---|---|---|---|---|---|
| **C14** | Field-split repair: q splits into four geometric places after base change; the scalar a_q ledger is not field-stable; correct case space = multiplicity vectors (b₁..b₄) up to S₄ with a+Σbᵢ ≤ cap. All pre-repair (a_t,a_q) survivor maps are conditional until reclassified. | FIELD_SPLIT_AUDIT.md (whole doc, esp. L17–58, L84–94); STATE.md status correction L5–16 | test_split_place_proofs.py ♦, test_split_place_ledger.py ♦ | externally authored audit + local exact checkers | 1 | C13 |
| **C15** | Split-place sigma-locus theorem: the T3 branch (d1=0, σ=0) is empty over every characteristic-zero field, for f31 and for f37 on d2≠0 — field-stable, both subcases (sub1 transfer with re-verified Mason–Stothers margins). | FIELD_SPLIT_AUDIT.md §sigma L98–213; sub1 transfer STATE.md L499–502, SPLIT_PLACE_LEDGER_SUB1.md L24–25 | t5_split_place_verify.py ♦; sub1 margins in sub1_cascade_verify.py ♦ | externally authored derivation + local exact checkers | 1 | C13, C14; Mason–Stothers |
| **C16** | The geometrically q-coprime a_t=7 stratum (T1 and T2) is dead. | FIELD_SPLIT_AUDIT.md L217–349 | t5_split_place_verify.py ♦ (finite algebra + 36-case degree table) | externally authored derivation + local checker | 1 | C13, C14 |
| **C17** | Sub2 terminal ledger: 327 geometric vectors (21 uniform + 306 partial); terminal valuation/degree pruning kills 81 strata; honest frontier after scoped proofs = 235 strata / 420 live T1/T2 branches; T3 dead globally. | SPLIT_PLACE_LEDGER.md L7–25 (counts L15–21); split_place_ledger.json (data) | split_place_ledger.py (generator; "do not hand-edit") + test_split_place_ledger.py ♦; terminal layer independently confirmed 654/654 by audit_cascade_kills.py ♦ | independently audited | 1 | C13, C14, C15, C16 |
| **C18** | Sub2 cascade depth-4 kills: of the 420 open branches, 390 killed, **30 survive** (a∈[5,10], bᵢ∈{0,1,3}), each with explicit residue obligations. | CASCADE_ENGINE_REPORT.md L14–36; CASCADE_KILL_AUDIT.md L3, L35–47 (654/654, 420/420, 390/390, 30/30, zero disagreements); cascade_cones.json (data) | test_cascade_engine.py ♦; **audit_cascade_kills.py ♦ (Codex-authored, spec-only, no engine code access)**; two kills re-derived by hand (CASCADE_ENGINE_REPORT.md L69–81) | **engine-proven AND independently audited** | 1 | C13, C17 |
| **C19** | Sub2 cone-lemma compression: all 390 kills fall into two certificate families — L (single-place: β=2 or β≥4 locally dead; live β∈{0,1,3}) and B (five budget patterns); verdicts match the engine on all 420 branches. | CASCADE_CONE_LEMMAS.md L3–40; cascade_cone_certificates.json (data) | cone_lemmas.py + test_cone_lemmas.py ♦ (explicitly compression, NOT independent — doc L3–5) | checker (engine-derived) | 2 | C18 |
| **C20** | Sub2 t-place coupling: t=y+1 as fifth place adds no tropical kills but kills 88/320 flag-cases → frontier 30 branches / 232 flag-cases with explicit finite-depth t-obligations (max 18). | STATE.md L472–486; cascade_cones_qt.json (data); CASCADE_TPLACE_AUDIT.md L9–23 (7,872/7,872 agreements, 0 disagreements) | engine regressions ♦; **audit_tplace_cases.py ◇ (independent, spec-only — NOT in run_tests.sh, see I8)** | independently audited (wiring gap) | 1 | C18 |

### 1.4 Subcase-2 per-cell exact kills (post-cascade)

| ID | Claim | Proven in | Checker | Audit status | Tier | Deps |
|---|---|---|---|---|---|---|
| **C21** | a_t=9 T2 (geometrically q-coprime) is infeasible (constant E + infinity domination). | T5_90_T2.md L1–85 | t5_90t2_verify.py ♦ | checker-only | 2 | C13, C16-method |
| **C22** | a_t=9 T1 constant-E cell is infeasible (exact convolution 238→226; final constant 29570349989420274657771126784·c⁵γ⁸ ≠ 0); 45/50 infinity cells excluded uniquely, leading coefficients of d1, then σ=0, v₀=525/32, d2 fully determined before the kill. | T5_90_T1.md §3 L137–246; STATE.md L329–353 | t5_90t1_verify.py ♦ + t5_90t1_constant_verify.py ♦ (sourced from f31_graded.txt) | checker-only | 2 | C13 |
| **C23** | a_t=9 T1 nonconstant cell: reduction (not closure) — s⁶\|σ, s³\|R, W fully determined; five compact local coefficient identities verified. | T5_90_T1.md §2; STATE.md L354–365 | t5_90t1_local_verify.py ♦ | checker-only; OPEN cell (in the 26-cell frontier) | 2 | C22 |
| **C24** | T2-column squeeze: 4 of the 12 T2 survivor cells proven infeasible — (5,(1,0,0,0)), (6,(1,0,0,0)), (6,(1,1,0,0)), (6,(1,1,1,0)) — 8/32 flag cases killed; under partial q-support the squeeze conclusion weakens to F²\|G (F the q-coprime remainder), CORRECTING PHASE_C_WORKLIST.md §4. | T5_T2_COLUMN.md L1–15, §1 L96–103, ledger L326–328 | t5_t2_column_verify.py ♦ (C1–C5; parses f31_graded.txt, no hand-copied h_l coefficients) | checker-only | 2 | C18, C20 |
| **C25** | T2 infinity convolution: 0 further cells killed, but all 8 open T2 cells strictly narrowed; flag case a9 b1000 g5=0 CLOSED; 11 residual degree-states removed (2 pattern-A, 9 d2=0 fixed-F). | T5_T2_INFINITY.md L7–16, §3, §4 ledger L234–248 | t5_t2_infinity_verify.py ♦ (I1–I5) | checker-only | 2 | C24 |
| **C26** | Sub2 Phase C worklist: the 232-case frontier carries 41 distinct obligation patterns (all 3621 tied monomials exact h_l members); targeting order and lemma candidates. | PHASE_C_WORKLIST.md L4–80 | phase_c_inventory.py ◇ (read-only re-derivation) | document + [judgment] tags; §4 T2 claim CORRECTED by C24; headline counts STALE (see I1) | 3 | C20 |

### 1.5 Subcase 1

| ID | Claim | Proven in | Checker | Audit status | Tier | Deps |
|---|---|---|---|---|---|---|
| **C27** | Sub1 cascade parameters: stripped caps (6,9,12,15), σ≤12, deg h_f ≤ 60−6f, budget a+Σb ≤ 15, standard regime a≤10; the proposed uniform cap deg g_l ≤ 15+3a is UNSUPPORTED for sub1 (corrected); safe terminal caps deg g7≤46 (T1) / g6≤48 (T2); sub2's 10+3a cap re-derived sound. | STATE.md L487–508; SPLIT_PLACE_LEDGER_SUB1.md L9–14 | sub1_cascade_verify.py ♦ (correction found by independent Codex derivation) | checker + independent derivation of the correction | 1 | C6, C13 |
| **C28** | Sub1 terminal ledger: 1333 strata; 136 terminal kills (all partial-support); 1197 open standard-regime (2178 branches); 26 alternate-regime strata (a∈[11,15], v<0). | SPLIT_PLACE_LEDGER_SUB1.md L16–26; split_place_ledger_sub1.json (data) | split_place_ledger_sub1.py (generator); terminal independently confirmed 2,614/2,614 by audit_cascade_kills_sub1.py ♦ | independently audited | 1 | C27, C15 |
| **C29** | Sub1 depth-4 cascade: 2178 branches → 1899 killed, **279 survive** — one a-independent family (26 (b,branch) pairs per a≤8; 24 at a=9, 21 at a=10; split 180 T1 / 99 T2). | STATE.md L569–578; CASCADE_KILL_AUDIT_SUB1.md L3–7, L87–108 (2178/2178, zero discrepancies); cascade_cones_sub1_depth4.json (data) | **audit_cascade_kills_sub1.py ♦ (Codex-authored port, no engine code access, own cap derivation)** | **engine-proven AND independently audited** | 1 | C27, C28 |
| **C30** | Sub1 cone lemmas: all 1899 kills compress to the same two families; live local values {0,1,2,3,5} (T1) / {0,1,3,5} (T2); only two budget patterns (T1:d1, T2:σ). | CASCADE_CONE_LEMMAS_SUB1.md L3–37; cascade_cone_certificates_sub1.json (data) | cone_lemmas.py --window sub1 + test_cone_lemmas.py ♦ (engine-derived) | checker (engine-derived) | 2 | C29 |
| **C31** | Sub1 t-place coupling: same 279 branches; 266/2,519 flag-cases die → **2253 flag-cases**, all with explicit t-obligations; t-depth law affine (depth = v = 30−3a, verified a=0..9). | STATE.md L579–584; cascade_cones_sub1_qt.json (data); CASCADE_TPLACE_AUDIT.md L16–23 (41,592/41,592 agreements) | audit_tplace_cases.py ◇ (independent; NOT in run_tests.sh, see I8) | independently audited (wiring gap) | 1 | C29 |
| **C32** | Sub1 Phase C worklist: 67 patterns / 21337 obligations over 2253 cases; sub2's 41 patterns ⊆ sub1's (29 identical, 12 deeper, only 6 genuinely new residue polynomials, all sub-cuts of h₄) — 91% of sub1 needs no algebra beyond the sub2 library; level-5 squeeze covers all 99 sub1 T2 cells. | PHASE_C_WORKLIST_SUB1.md L4–124, L343 | phase_c_inventory_sub1.py ◇ (read-only; 26965 tied monomials, 0 mismatches) | document + [judgment] tags | 3 | C31 |

### 1.6 Alternate regime (sub1, a ∈ [11,15])

| ID | Claim | Proven in | Checker | Audit status | Tier | Deps |
|---|---|---|---|---|---|---|
| **C33** | Flipped cascade: for v=30−3a<0 the t-order telescopes to the constant 210 (F = t²¹⁰G′), descending cascade anchored at h₇=8192·d1²; q-place transitions and terminal caps survive verbatim; first-level parity/degree kills 19/52 branches (T1 7, T2 12; 6 strata fully dead) → 33 branches / 20 strata. | ALT_REGIME.md L28–33, L120–154 | alt_regime_verify.py ♦ (7 groups; random a=12 window instance) | checker-only — **no independent audit of the alternate regime (gap G3)** | 2 | C27, C28 |
| **C34** | Flipped-cascade levels 3/2 (h6/h5): 6 more kills (all T1) → **27 residual branches (13 T1 + 14 T2)**, with terminal UFD residual normal form and survivor caps. | ALT_REGIME_L2.md L8–17, §3 L142–155, table L162–184 | alt_regime_l2_verify.py ♦ (checks A–E) | checker-only (deliverable verified locally after the authoring runtime lost its registry — STATE.md L613–620) | 2 | C33 |

### 1.7 Superseded / conditional / supporting results

| ID | Claim | Proven in | Checker | Audit status | Tier | Deps |
|---|---|---|---|---|---|---|
| **C35** | Pre-repair (a_t,a_q) strata kills — (0,0),(1,0),(2,0),(3,0),(4,0),(0,1) (T5_MULTIPLACE Thm 1), (5,0),(1,1) (T5_STRATA_50_11), the (6,x)/(a_q=1,2) column work (T5_60_T1/T2, T5_T1_AQ12, T5_STRATUM_10_0, T5_70) — valid ONLY as geometrically-q-coprime / uniform-q^r statements per the field-split repair; the split-place ledger and cascade (C17–C18) are the authoritative replacements. | T5_MULTIPLACE.md, T5_STRATA_50_11.md, T5_60_T1.md, T5_60_T2.md, T5_T1_AQ12.md, T5_STRATUM_10_0.md, T5_70.md; scoping: STATE.md L5–16, FIELD_SPLIT_AUDIT.md L84–94 | t5_multiplace_verify.py ◇, t5_strata50_11_verify.py ◇, t5_60t1_verify.py ◇, t5_60t2_verify.py ◇, t5_t1_aq12_verify.py ◇, t5_stratum100_verify.py ◇, t5_70_verify.py ◇ | CONDITIONAL (scoped); the docs themselves carry NO conditional banner (see I5) | 3 | C13; superseded by C17–C18 |
| **C36** | f37 graded/frontier analysis (grading non-uniform, h37 irreducible, ledger-starvation kills zero f37 strata, free family permanently live at resultant level). | T5_F37_GRADED.md; F37_FRONTIER.md | t5_f37_verify.py ◇ | MOOT — superseded wholesale by C11 (see I3) | 3 | superseded by C11 |
| **C37** | Numerical infeasibility evidence: all four factor×subcase window systems show genuine positive floors under central-difference polish (f37_sub1 5e-7 < f37_sub2 1.2e-6 < f31_sub1 8e-6 < f31_sub2 2.3e-5); independent two-machine replication of sub1. | STATE.md sessions L119–274 | jetlift.py control ♦ (positive control only); campaign logs/pickles (data) | EVIDENCE ONLY, explicitly not proof; f37 halves moot after C11 | 3 | C6 |
| **C38** | NulLA/linear-certificate no-go: no homogeneous block-1-only Nullstellensatz certificate exists for any N (dim ≥ 3 witness family); monolithic Gröbner/msolve routes are dead. | LINALG_CERT_REPORT.md L14–16, L106–127, L167–184; MSOLVE_REPORT.md (stub — see I6) | none committed (scratchpad F_32003 rank computations) | document-only (negative/method result, not load-bearing for C0) | 3 | C13 |
| **C39** | Corner-144 comparison (conditional): the (108,144) reappearance of corner (8,28) shares the parametric operator/D-exponent skeleton and the one-t-place-plus-quartic geometry (q-multiplicity 69), but NOT the numeric forcing signature; supports a (power-pair, multiplicity) template. | CORNER_144_COMPARISON.md L3–8, table L183–188; STATE.md L399–418 | corner144_verify.py ◇ (51/51 pass) | checker-only (not in suite); conditional on the unreduced polygon boundary; not load-bearing for C0 | 2\* | C1; GGV5 |
| **C40** | Next-cases enumeration: (75,125) is the single next case; 24 candidates in [125,150]; corner (8,28) reappears at (108,144); reproduces the paper's 10-case table exactly. | NEXT_CASES.md L143–182 | paper_src/next_cases.py ◇ | checker-only; context, not load-bearing for C0 | 2\*/4 | GGHV22, GGV5 |

### 1.8-pre Phase D infinity layer (2026-07-22, post-inventory additions)

| ID | Claim | Proven in | Checker | Audit status | Tier | Deps |
|---|---|---|---|---|---|---|
| **C41** | Max-plus infinity layer: infinity as a sixth place (v_inf = −deg); unique max forces, ties drop only as recorded obligations; full chain terminal..1 + level-0 anchor; degrees sandwiched between finite valuation sums and window caps. Regression ladder R0–R5 replays T5_90_T2 §2, T5_90_T1 §3 (43/50), all seven T5_T2_COLUMN margins, and the T2-column cell verdicts (12/12). | CASCADE_INF_REPORT.md §1–3, §6 | test_cascade_inf.py ♦ | same-author checker; audited artifacts byte-identical | 2 | C13, C18, C20, C21, C22, C24 |
| **C42** | Infinity tie equations: depth-1 initial form of an infinity degree tie = the residue-lemma initial form with leading coefficients as unknowns (backbone P6/P10/P11 verified); lc(u)=lc(Φ)=−1024/3315 from source; C08/C20 act as forbidden drops at infinity under --residue-kills (no Q-points). | CASCADE_INF_REPORT.md §4 | cascade_inf_ties_verify.py ♦ | same-author checker; leans on tier-2 residue_lemmas_verify | 2 | C41, RESIDUE_LEMMAS (C08/C20) |
| **C43** | Infinity sweeps, BOTH windows: sub2 **26 branches / 220 flag cases** (four T2 column cells killed endogenously via the C24 squeeze as a g6 budget cap); sub1 **171 branches / 1145 flag cases — 108 NEW branch kills** (a=0/1 closed entirely, a=2 keeps 2; the 26-family's a-independence breaks at low a). No unsound new survivors vs the audited rl artifacts. | STATE.md sub2/sub1-infinity entries | test_cascade_inf.py ♦ (R5); **audit_inf_cases.py ♦ (Codex-authored, spec-only from CASCADE_INF_REPORT.md, no engine access: sub2 420/420, sub1 2178/2178, zero disagreements)** | **engine-proven AND independently audited** | 1 | C41, C42, C24, C18, C20, C29 |
| **C44** | Alternate-regime degree layer: flipped chain T r_{f−1} = E^{3(7−f)} h_f + u r_f (top anchor T r_6 = h_7, bottom close E²¹h_0 + u r_0 = 0) with max-plus degree identities; degree sweep kills 33670/38360 states, 27 branches open. | ALT_REGIME_INF.md; ALT_INF_SWEEP.md | alt_regime_inf_verify.py ♦; alt_inf_sweep_verify.py ♦; **audit_alt_regime.py ♦ (Codex spec-only: 25 kills re-derived, 27 OPEN with full state partitions, exit 0 — also covers C33/C34, closing gap G3)** | **independently audited** | 1 | C33, C34 |
| **C46** | No-jet-kill theorem: all ten t-place tied supports occurring on the 1145 surviving sub1 cases have smooth rational points, so the a-quantified jet tower (affine law depth = 30−3a) is solvable to every order — t-place obligations are CONSTRAINT for all a ∈ [2,10], never local kills; Phase D closures must be global/infinity. | RESIDUE_LEMMAS_DEPTH.md | residue_lemmas_depth_verify.py ♦ (V1–V7) | same-author exact checker | 2 | C31, C43 |
| **C45** | Phase D residual worklist: complete degree-state lists per surviving sub2 flag case (220 cases / 7888 states); the 8 open T2 cells' (deg e, deg σ) sets equal T5_T2_COLUMN's R-tables exactly; a9 T1 constant-E states are exactly the 7 tied cells of T5_90_T1 §3. | phase_d_states_sub2.json | phase_d_states.py ◇ (generator; R-table cross-check executed at generation) | engine-derived worklist (not a proof artifact) | 3 | C43 |

**Claim-graph row count: 47 (C0–C46).**

### 1.8 Tier distribution

| Tier | Rows | IDs |
|---|---|---|
| 1 | 10 | C14, C15, C16, C17, C18, C20, C27, C28, C29, C31 |
| 2 | 11 | C11, C12, C19, C21, C22, C23, C24, C25, C30, C33, C34 |
| 2/4 | 2 | C2, C3 |
| 2\* | 8 | C4, C5, C7, C8, C10, C13, C39, C40 (C40 also cites tier 4) |
| 3 | 9 | C1 (3/4), C6, C9, C26, C32, C35, C36, C37, C38 |
| open (no tier) | 1 | C0 (the target itself) |

Summary: **10 tier-1, 13 tier-2 (incl. 2/4), 8 tier-2\*, 9 tier-3, C0 open (41 rows total).**
Every load-bearing kill on the current frontier path (C14–C31, C33–C34) is tier 1 or 2;
the tier-2\* and tier-3 rows are either upstream derivation audits with unwired checkers
(C4–C10, C13 — see I8) or non-load-bearing/superseded material.

---

## 2. Frontier table — the live fronts, exactly, as of 2026-07-22

> **UPDATE (2026-07-22, post-infinity-layer):** the figures in this section
> predate the max-plus infinity sweeps. The current frontier is **26 sub2
> cells (220 flag cases) + 171 sub1 branches (1145 flag cases) + 27
> alternate branches** — authoritative: the machine-generated `FRONTIER.md`
> (from `cascade_cones_qt_inf_rl.json`, `cascade_cones_sub1_qt_inf_rl.json`,
> `alt_inf_sweep.json`; closes G6). The sub1 figure (279 → 171, 108 new
> kills, a=0/1 closed) is engine-proven PENDING the spec-only audit (C43+).

Everything not listed here is closed: the f37 branch (C11), the d₋₁≡0 branch (C10), the
T3 sigma-locus (C15), all sub2 strata outside the 30 audited survivors (C18), 4 of the 12
sub2 T2 survivor cells (C24), all sub1 standard-regime branches outside the 279 (C29), and
25 of the 52 alternate-regime branches (C33+C34).

| Front | Size | Composition | Authoritative artifact(s) | Next machinery |
|---|---|---|---|---|
| **Sub2 cells** | **26 cells** | 18 T1 + 8 T2. = the 30 independently audited survivors (CASCADE_KILL_AUDIT.md L71–81; a∈[5,10]) minus the 4 T2 column kills (T5_T2_COLUMN.md L326–328). Open T2 cells: a7 {1000,1100,1110,3000}, a8 {0000,1000,1100}, a9 {1000} — all narrowed by T5_T2_INFINITY.md §4; flag case a9 b1000 g5=0 closed. The T1 tail includes the a_t=9 nonconstant cell (C23). | cascade_cones_qt.json (30-branch/232-case level) + T5_T2_COLUMN.md + T5_T2_INFINITY.md ledgers; obligations: PHASE_C_WORKLIST.md (headline stale, see I1). Post-column case count ≈ 232−8−1 = 223 `[inferred — no artifact records it]` | Phase C residue lemmas; Phase D infinity layer (needed by all a≥7 T2 ties and the T1 tail) |
| **Sub1 standard regime** | **279 branches** | One a-independent family: 26 (b,branch) pairs per a for a≤8 (17 T1 + 9 T2), 24 at a=9, 21 at a=10; overall 180 T1 / 99 T2; 2253 flag-cases with a-independent q-obligations and affine t-depth law. | cascade_cones_sub1_qt.json; CASCADE_KILL_AUDIT_SUB1.md (audit); PHASE_C_WORKLIST_SUB1.md (obligations); CASCADE_CONE_LEMMAS_SUB1.md (compression) | Reuse sub2 residue library (91% coverage, C32); level-5 squeeze for the 99 T2 cells; a-quantified lemmas |
| **Sub1 alternate regime** | **27 branches** | 13 T1 + 14 T2, in ≤20 strata, a∈[11,15]; terminal UFD residual normal form + survivor caps established. | ALT_REGIME_L2.md (verdict table L162–184); upstream ALT_REGIME.md K/O table L120–147 | Lower flipped-cascade levels; residue congruences; NEEDS an independent audit (gap G3) |

Caution: sub1's "26-family" (26 (b,branch) pairs per a) and sub2's "26 cells" are
**unrelated quantities that happen to share the number 26** — do not conflate (see I9).

---

## 3. Independent-audit map (what earns tier 1)

| Independent artifact | Author independence | Covers | In run_tests.sh? |
|---|---|---|---|
| audit_cascade_kills.py | Codex-authored from CASCADE_ENGINE_REPORT.md semantics only; own f31_graded.txt parser; no engine code access | Sub2: 654/654 terminal, 420/420 branches (390 kills, 30 survivors) | ♦ yes (--quiet) |
| audit_cascade_kills_sub1.py | Codex-authored port, zero references to cascade_engine.py, own cap derivation | Sub1: 2,614/2,614 terminal, 2,178/2,178 depth-4 (1,899 kills, 279 survivors), depth-5 ⊂ depth-4 | ♦ yes (--quiet) |
| audit_tplace_cases.py | Spec-only, no engine imports | Sub2 7,872 + sub1 41,592 flag-case verdicts, zero disagreements | ◇ **NO** (see I8) |
| FIELD_SPLIT_AUDIT.md | External proof-program audit (separate author) | Sigma-locus theorem, a_t=7 kill, field-split framework | via t5_split_place_verify.py ♦ |
| sub1 cap correction | Found by independent Codex derivation (STATE.md L493) | C27 | via sub1_cascade_verify.py ♦ |
| Two-machine numeric replication (2026-07-21 merge) | Independent campaigns, different seeds/machines | C37 (evidence only) | — |

NOT independent (explicitly): cone_lemmas.py output (engine-table compression, both
CASCADE_CONE_LEMMAS*.md L3–5); every t5_*_verify.py / alt_regime*_verify.py (same-author
exact checkers); split_place_ledger*.py (self-generating ledgers).

---

## 4. Published results used as stated (tier 4)

### 4.1 Load-bearing (in the proof chain of C0)

| Result | Used for | Invoked at |
|---|---|---|
| **GGHV22** (arXiv:2204.14178) **Prop 4.3**, case (8,28) | The case statement, subcase polygons, [P,Q]=x² normalization — the problem itself (C1) | T3_WINDOW_AUDIT.md §1 L17–26 (verbatim quote, source lines 1000–1007); T6_PREMISES.md L53–69; HANDOFF.md L30–33 |
| **GGV1** (arXiv:1401.1784) **Prop 1.13** (bracket-valuation inequality) | Premise 1 and 2 (C2, C3): triggers alignment via strict gap 2<19 | T6_PREMISES.md L31–37 (verbatim), L74, L159–179; T6_SELECTION_AUDIT.md L88–91; T3_WINDOW_AUDIT.md L112–115 |
| **GGV1** **Prop 2.1** (aligned elements are common powers) | Premise 1 (C2): existence of R with ℓ(P)=R², ℓ(Q)=R³ | T6_PREMISES.md L39–44 (verbatim), L79–115; T6_SELECTION_AUDIT.md L88–91 |
| **GGHV22 §4** (the closed (9,27) twin case, source lines 1399–1596, esp. 1508–1546) | The complete setup template transported t=3→t=4 (C2–C4) | T6_PREMISES.md L23, L155, L204; T3_WINDOW_AUDIT.md L38–48 |
| **GGHV22 §6** (the (7,21) worked case) | Generator validation target (C8) | AUDIT.md §B L63–64; STATE.md L138–139 |
| **GGHV22 Props 5.2(3)(4), 5.6** (valuation induction) | Template for the envelope-bound induction (C6) — the induction itself was redone here | STATE.md L35; T3_WINDOW_AUDIT.md §3 |
| **Mason–Stothers theorem** (char 0) | Sigma-locus kills (C15) and pre-repair Theorem 2 machinery | FIELD_SPLIT_AUDIT.md §4 L181–196; T5_MULTIPLACE.md §4 L231–245; T5_F37_GRADED.md §7; SPLIT_PLACE_LEDGER_SUB1.md L25 |

### 4.2 Cited but not load-bearing

| Result | Status | Where |
|---|---|---|
| van den Essen Prop 10.2.6 | Mentioned as unused ammunition only | STATE.md L80 |
| GGV5 (arXiv:1708.07936) chain/corner recipe | Forward-looking only (C39, C40) | CORNER_144_COMPARISON.md, NEXT_CASES.md |
| GGV2 Prop 3.29/Rmk 3.31, GGV3, GGV6 Thm 7.3, Moh 1983, Heitmann 1990, Orevkov 2001 | Historical context table | NEXT_CASES.md §4 L191–206 |
| **GGV3 §1** (cited by GGHV22 for the strip descent) | **Cited in prose but source not in paper_src/**; the argument is filled here from GGV1 Props 1.13/2.1 instead (see gap G7) | T6_PREMISES.md L170–173, L209 |

---

## 5. Inconsistency sweep

Reported only — nothing has been fixed. Ordered by writeup risk.

- **I1. PHASE_C_WORKLIST.md is stale AND contains a corrected claim, with no banner.**
  (a) §4 L296–309: "all 12 T2 survivor cells … hypotheses are met by construction … One
  lemma … reduces every T2 cell" and L300–303 (squeeze applies "verbatim") — explicitly
  corrected by T5_T2_COLUMN.md L6–8 and L100–103 (under partial q-support only F²|G holds,
  F the q-coprime remainder; 4 cells killed, 8 remain open with genuine infinity ties).
  (b) Headline counts L4, L74 ("30-branch / 232-case", "12 T2 cells" L305) predate the
  column results; the current frontier is 26 cells / 8 open T2 cells. The worklist carries
  no correction notice pointing at T5_T2_COLUMN.md.
- **I2. HANDOFF.md presents a superseded frontier as current.** L18–20: "leave 235 strata
  and 420 live T1/T2 branch records" — that was the pre-cascade frontier; current is 26
  sub2 cells + 279 sub1 branches + 27 alternate branches. L16–18 states the f37 status as
  free-family-only exclusion, understating the full closure (C11). Task list T1–T6
  (L56–76) is largely executed (T1/T2 done as the cascade engine; T5's window-bound
  concern resolved by T3_WINDOW_AUDIT.md). "Read this first" makes this the entry
  document, so the staleness is high-risk. (FIELD_SPLIT_AUDIT.md item 7 L378–381 already
  demanded this reconciliation; it has not happened for HANDOFF.md.)
- **I3. F37_FRONTIER.md reads as live but is wholly superseded.** L278: the decisive
  saturation "has not been tried"; L254–255 tags the artifact hypothesis [I] (unproven);
  §5 L294–328 recommends it as future work. F37_SATURATION_REPORT.md executes the
  stronger membership result and explicitly supersedes the stratum-by-stratum f37 program
  (L21–23, L99–105). F37_FRONTIER.md carries no supersession banner. Same (milder) issue:
  F37_FREE_FAMILY_SYSTEM.md's forward scope L171–176 ("future f37 work may focus on
  (d2,d1)≠(0,0)") is closed by C11; T5_F37_GRADED.md is moot.
- **I4. T5_NOTES.md predates both the field-split repair and the f37 closure.** L11 and
  L186–190 treat f37 as a live alternative; L192–209 ("Not closed … recommended next
  work" = the (y+1)-adic by-hand count) describes a strategy since replaced by the
  split-place ledger + cascade engine. No banner.
- **I5. The seven pre-repair strata docs carry no conditional banner.** T5_MULTIPLACE.md,
  T5_STRATA_50_11.md, T5_60_T1.md, T5_60_T2.md, T5_70.md, T5_STRATUM_10_0.md,
  T5_T1_AQ12.md state (a_t,a_q) survivor maps unconditionally; STATE.md L5–16 declares
  every such map "conditional until reclassified by the split-place audit," but only
  STATE.md says so (grep: no "conditional"/"FIELD_SPLIT" match inside T5_MULTIPLACE.md or
  T5_STRATA_50_11.md). A reader entering via those docs sees retracted-in-scope claims
  presented as unconditional. This is exactly the FIELD_SPLIT-episode drift mode.
- **I6. MSOLVE_REPORT.md is an empty placeholder.** L1–3: "(Results being filled in —
  placeholder, see final version below.)" — there is no final version below. The msolve
  results it promises live only in STATE.md prose (L306–308).
- **I7. The 30-branch → 26-cell transition is recorded nowhere as a single statement.**
  CASCADE_ENGINE_REPORT.md L18/L26 and CASCADE_KILL_AUDIT.md L45 state the audited
  frontier as 30; cascade_cones_qt.json still lists 30 branches / 232 cases; the 26-cell
  figure exists only by combining those with T5_T2_COLUMN.md L326–328, and appears as a
  number only in STATE.md L553 and CASCADE_ENGINE_PLAN.md L184. No JSON artifact reflects
  the post-column frontier; the post-column flag-case count (≈223) is nowhere recorded.
- **I8. run_tests.sh omissions.** Not wired in: **audit_tplace_cases.py (an INDEPENDENT
  audit — the only tier-1 artifact outside the suite)**, verify_derivation.py (the 48-check
  derivation audit backing C4/C5/C7/C10), verify_graded.py (C13), corner144_verify.py,
  t5_multiplace_verify.py, t5_strata50_11_verify.py, t5_60t1_verify.py, t5_60t2_verify.py,
  t5_70_verify.py, t5_stratum100_verify.py, t5_t1_aq12_verify.py, t5_f37_verify.py,
  paper_src/next_cases.py. The suite therefore does not currently re-verify the upstream
  derivation chain on which every frontier kill depends.
- **I9. "26" collision.** Sub2's frontier (26 cells) and sub1's per-a survivor family
  (26 (b,branch) pairs) are unrelated; STATE.md L575 and PHASE_C_WORKLIST_SUB1.md §3 use
  the same number for the latter. Also "26" appears a third time as the count of
  alternate-regime strata (SPLIT_PLACE_LEDGER_SUB1.md L22). High confusion potential in
  any summary.
- **I10. Minor log corruption in STATE.md.** The "A_T=9 T1 REDUCTION" entry's final
  sentence is split across the subsequent entry: L339 ends "and the five" and the
  continuation "residual cells have an explicit low-parameter form." appears at L353,
  after the intervening constant-cell entry.
- **I11. ALT_REGIME_L2.md provenance note.** STATE.md L613–620 records that the authoring
  (codex) runtime lost its job registry before reporting; deliverables were verified
  locally and landed on that basis. Not an error, but the alternate-regime chain has the
  weakest provenance of any live front (see G3).

---

## 6. Writeup gap list — what a referee will ask for

- **G1. The theorem is not yet provable: the frontier is nonempty.** 26 sub2 cells + 279
  sub1 branches + 27 alternate branches remain (all with explicit finite obligation
  lists, but obligations are not proofs — both worklists say so). This is the primary
  gap; everything else is presentation.
- **G2. Envelope-bound formal writeup and checker.** C6 (deg ≤ 14w/15w, ord ≥ 12w) is
  the foundation of every window and cap, yet exists only as T3_WINDOW_AUDIT.md prose
  with in-session sympy checks: no committed checker, no referee-grade induction
  writeup. Same for the Prop 4.3 transcription (C1, document-only).
- **G3. Independent-audit coverage for the per-cell and alternate-regime kills.** The
  cascade layer is tier 1, but every post-cascade kill (C21–C25) and the whole
  alternate regime (C33–C34) rest on same-author checkers; the alternate regime is not
  covered by any of the three independent auditors (both audit_cascade_kills* scripts
  are standard-regime only), and its authoring session had degraded provenance (I11).
  A spec-only auditor for the flipped cascade and the T2-column/infinity arguments is
  the missing piece to bring the whole live chain to tier 1.
- **G4. The Phase D infinity layer does not exist.** Every a≥7 sub2 T2 tie, the sub2 T1
  tail, and the zero-obligation sub1 tier (a10,(1,1,1,1)) are declared to need it
  (T5_T2_COLUMN.md L328-context, T5_T2_INFINITY.md L7–16, PHASE_C_WORKLIST_SUB1.md
  L370); only per-cell convolutions exist today.
- **G5. Narrative restructure around the ideal-membership result.** Most documents (and
  STATE.md items 5–6) still narrate resultant-first with f37 as a live branch; the
  writeup should lead with C11 (elimination ideal = ⟨f31⟩), which retires C9's
  Singular-trust question and the entire f37 literature in this repo. The supersession
  banners demanded by I2–I5 are the working-level version of this gap.
- **G6. Machine-generated frontier summary.** FIELD_SPLIT_AUDIT.md L378–381 prescribed
  generating the human summary from the machine ledger to stop drift; the ledgers do
  this, but the cross-front frontier (26/279/27) is still hand-maintained prose in
  STATE.md, and I1/I2/I7 show the drift recurring. A single generated FRONTIER.md (or
  this inventory regenerated) from the JSONs + kill ledgers would close it.
- **G7. Citation debt: GGV3.** T6_PREMISES.md fills GGHV22's GGV3-§1 citation from GGV1
  Props 1.13/2.1 (L170–173, L209); the writeup must either verify against GGV3 itself
  or state the substitution explicitly. Relatedly, the plan to contact the
  GGHV/GGV authors before writeup (STATE.md L88–91, HANDOFF.md L78–80) is a stated
  obligation, not yet an action.
- **G8. Status of the numerics and T4.** The four positive floors (C37) are evidence
  only; base-region completeness (T4) was never done. The writeup must either scope the
  numerics as motivation or complete T4 — a referee will ask which.
