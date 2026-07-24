# Attack on the (72,108) case — state of the derivation
### Target: the last open case below degree 125 for the plane Jacobian Conjecture
### (GGV-Horruitiner 2022, arXiv:2204.14178, Prop 4.3 case (8,28) — left open "for lack of computing power")

> **Status correction, 2026-07-22.** The later T5 entries below used the
> quartic `q` as one irreducible prime over `Q`. This is not a geometric
> stratification: after base change, `q=p1*p2*p3*p4`, and `e` may contain only
> some `pi`. Therefore every old `(a_t,a_q)` survivor map is conditional until
> reclassified by the split-place audit. See `FIELD_SPLIT_AUDIT.md`. The
> upstream reduction and the `t`-cascade remain valid. The repaired global
> sigma-locus theorem is field-stable, and the geometrically `q`-coprime
> `a_t=7` branch is newly closed. The `f37` resultant alone also has the exact
> free family `d2=d1=0`, so its generic numerical floors are not an
> infeasibility proof even at the resultant level. Restoring the three smallest
> pre-resultant equations now proves that this family has no lift to the
> original system (`F37_FREE_FAMILY_SYSTEM.md`).

## The open problem, exactly
Do there exist P, Q ∈ K[x,y], [P,Q] = x², with Newton polygons
  subcase (2): N(P) = {(0,0),(1,0),(8,14),(8,16)},        N(Q) = {(0,0),(2,1),(12,21),(12,24)}
  subcase (1): same plus corners (0,8) in N(P), (0,12) in N(Q)
If NO: the JC(2) counterexample bound rises 108 → 125 (new theorem).
If YES: explicit constraints on a candidate (72,108) counterexample.

## Derived so far (all verified computationally; audit list at bottom)
1. Normalization: l_{1,0}(P) = R², l_{1,0}(Q) = R³, R = x⁴C₄, C₄ = y⁷(y+1).
   Series C ∈ K[y,C₄⁻¹]((x⁻¹)), P = C², Q = C³ + λC⁻¹ + F, v_{1,0}(F) = −5.
2. Forcing term (Prop-5.4 analogue): [P,Q] = x² ⇒ ODE
   8y(y+1)f₁′ − 14(8y+7)f₁ = y⁸(y+1)²,  f₁ := C₄³F₋₅.
   Unique solution: f₁ = −y⁸(y+1)²(2048y⁴ − 512y³ + 320y² − 240y + 195)/6630.
   Quartic factor separable; y, (y+1) do not divide it. Same structure as their solved cases.
3. D-transformation: D_k := C_k·C₄^(7−2k) ∈ K[y] (verified on random polygon
   instances). Envelope bounds (weight w = 4−k): deg_y ≤ 14w, ord_y ≥ 12w
   (empirically tight; to be proven via v_{2,−1}(C)=1, v_{−2,1}(C)=0 induction,
   template = their Prop 5.2(3)(4) + 5.6). Generic (y+1)-multiplicity: 0.
4. Shift x ↦ x − D₃/4 gives d₃ = 0; system in d₂,d₁,d₀,d₋₁,…,d₋₁₃ (d₋₁₂ absent):
   (D̃²)₋ₖ = 0 for k = 1..7, 9;  (D̃³)₋ⱼ = 0 for j = 1,2,3;  (D̃³)₋₅ + Φ = 0,
   where Φ := F₋₅C₄³¹ = f₁C₄²⁸ ∈ K[y] explicit:
     deg Φ = 238, ord_y Φ = 204, mult_{y+1} Φ = 30,
     trailing coeff (y²⁰⁴) = −1/34, leading coeff (y²³⁸) = −1024/3315.
   (Equation (D̃³)₋₄ excluded: it only defines λ. Generator validated by
   regenerating the published (7,21) system verbatim.)
5. ELIMINATION DONE (the step that stopped the 2022 authors): linear
   substitutions (all Q-side equations are linear in d₋₄) + one Singular
   resultant. Master identity in the ideal:
     f31 · f37 · d₋₁²¹ ≡ 0 in K[y]
   with f31 = 102-term factor, f37 = 618-term factor (files enclosed).
   Since K[y] is a domain: one factor vanishes identically. Case tree:
   - d₋₁ ≡ 0: **PROVEN IMPOSSIBLE** (this session, direct substitution, no
     denominator artifacts). G1|_{d₋₁=0} = 3d₋₂d₋₃ ≡ 0, so over the domain
     K[y] either d₋₂ ≡ 0 or d₋₃ ≡ 0. Leg d₋₂ ≡ 0: G2 = (3/2)d₋₃² forces
     d₋₃ ≡ 0, then G5body ≡ 0 forces Φ = 0. Leg d₋₃ ≡ 0, d₋₂ ≢ 0:
     G3 = −(3/2)d₁d₋₂² forces d₁ ≡ 0, then G5body = −3d₁d₋₂d₋₄ ≡ 0 forces
     Φ = 0. Both contradict Φ = f₁C₄²⁸ ≠ 0. Case tree is now exactly
     {f31 ≡ 0} ∪ {f37 ≡ 0}.
   - f31 ≡ 0: THE MAIN BRANCH. Numeric solutions of the system at 3 random
     parameter sets (21 verified solutions, 60-digit precision) all satisfy
     f31 = 0 and never f37 — f31 is the genuine generic component.
   - f37 ≡ 0: must still be processed (the actual K(y)-point could a priori
     sit on this component even though generic C-points don't).
6. Structure of f31: weighted-homogeneous of weight 125 under
   w(d₂,d₁,d₀,d₋₁,Φ) = (2,3,4,5,17). 102 monomials, max Φ-power 7,
   pure anchor d₋₁²⁵. Because the envelope bounds are exactly (14w, 12w),
   every term of f31 lives in the y-window [1500, 1750] — raw degree
   counting CANNOT close the case (as in their proofs, the kill must come
   from the finer slice/multiplicity structure).

## Endgame plan (next sessions)
A. Coefficient-window system, subcase (2): envelope forces
   d₂, d₁, d₀, d₋₁ to have at most 5, 7, 9, 11 nonzero y-coefficients
   (32 unknowns total). f31(d₂(y),…,Φ(y)) ≡ 0 gives ~251 scalar equations
   (y-slices 1500..1750). Massively overdetermined. Solve slice-by-slice
   (jet lifting): bottom slice is f31(τ₂,τ₁,τ₀,τ₋₁, −1/34) = 0 in the
   trailing coefficients; lift upward; find the obstruction height or a
   solution. Repeat for f37. Same machinery, wider window, for subcase (1).
B. Symbolic confirmation that d₋₁ ≡ 0 forces Φ = 0 (small elimination).
C. Prove the envelope bounds (v_{2,−1}, v_{−2,1} induction per template);
   extra ammunition if needed: P(x,0) has x-degree exactly 1 (corner (1,0)),
   giving trailing-coefficient constraints beyond the envelope
   (their tool: van den Essen Prop 10.2.6).

## Audit obligations before any claim
- Independent re-derivation of: the F-normalization (Remark 5.3 analogue for
  [P,Q] = x²), the equation-selection argument (k = 8 vacuity, λ-isolation
  in (D̃³)₋₄, (D̃⁻¹)₋₅ = 0), and the envelope-bound induction.
- Verify subcase polygons against the paper's Prop 4.3 statement (the two
  cases and their [P,Q] = x² normalization) — transcription from the fetched
  text should be checked against the published PDF.
- The 2022 authors should be contacted before/at writeup: this is their
  program and their open question; the right venue for a result here is
  with or alongside them.

## Session log
- ODE + harness + generator validation + t=3 universal relation
  (18Φd₁d₋₁⁶ + 8Φ³ + 27d₀d₋₁⁹, y-free — unstated in their papers).
- t=4 elimination via Singular 4.3.2; factor export; genuine-factor test;
  homogeneity + anchor analysis; d₋₁ = 0 probe.
- ENDGAME SESSION 1: d₋₁ ≡ 0 branch proven impossible (see above). Window
  feasibility machinery built for f31 ≡ 0, subcase (2): FFT evaluator over
  the 32-coefficient window (A:5, B:7, C:9, E:11), weighted-scaling gauge
  to condition Ψ, penalty on ‖E‖. First LM multi-start (8 runs): all local
  minima strictly positive (best ‖G‖ ≈ 0.107), optimizer attracted toward
  the dead E→0 branch. Evidence leans infeasible; INCONCLUSIVE as proof.
  Next: jet-lifting parametrization (solve slices 0–10 exactly, optimize
  the ~21 residual dof against the 240 consistency slices), then f37, then
  subcase (1), then exact certificates for whatever the numerics indicate.
- ENDGAME SESSION 2 (jet lifting): implemented and validated the Hensel-style
  lifter — slices 0–10 solved exactly (verified to machine precision relative
  to term scale; the constant-gradient structure works as derived). Key
  numerical lessons: all tolerances must be relative (term scales reach 1e24
  at |e₀|≈7 via the 25th power); balanced gauge for Ψ (γ = √(max·|Ψ(0)|));
  null bases via SVD, not QR. Data so far (n=2 completed optimizations —
  throughput bug limits statistics, ~200s/run vs expected ~10s): relative
  consistency residuals stuck at 2.3e-1 and 1.1e-2, never approaching zero;
  >90% of residual energy in slices y^11–y^40. CAVEAT: with n=2 the
  concentration pattern may reflect optimizer dynamics, not true obstruction
  location. Next: profile/fix throughput, 50+ run statistics, then exact
  analysis of the early consistency slices if concentration persists.
- ENDGAME SESSION 3 (decisive numerics, f31 subcase (2)): fixed throughput
  (lambdified base generation, ~7s/run); custom damped Gauss–Newton replaces
  scipy LM; per-slice magnitude normalization (higher slices carry ~1e6
  multinomial factors — global scaling was wrong). POSITIVE CONTROL
  VALIDATED: with consistency truncated to 20 slices (fewer conditions than
  the 21 dof) the optimizer converges to its finite-difference floor (~2e-6)
  — solutions are found when they exist. FULL RANGE: 128 independent
  lifting runs; minima span 2.3e-5 … 2.5, median ~2.5e-3. Central-difference
  polish (noise floor ~1e-10) of the best three: trajectories DEAD FLAT
  (2.323e-5 -> 2.323e-5). Genuine positive floors, not partial convergence.
  CONCLUSION (strong evidence, not proof): the f31 window system for
  subcase (2) is infeasible — the (72,108) case is dying. Remaining for
  proof: f37 branch (same harness), subcase (1) windows, base-region
  completeness (sampled |e₀| in [0.02,4], moderate parameters), exact
  certificate, full derivation audit. Best minima saved (candidate_*.pkl).
- FRESH-REPO SESSION (T1 REPRODUCE, 2026-07-20, new workspace repo): rebuilt
  everything from scratch in a clean container (numpy 2.4.6, scipy 1.17.1,
  sympy 1.14.0, Singular 4.3.2). `regenerate_system.py`: t=3 generator
  validation OK (reproduces the published (7,21) equation verbatim); t=4
  system + Ain/Bin regenerated. `run_singular.sh`: resultant = 3228 terms;
  factorize gives factor_2 (102 terms, deg 31), factor_3 (618 terms, deg 37),
  factor_4 = m1 (mult 21). POLYNOMIAL-LEVEL COMPARISON to the enclosed
  f31/f37 (as sympy Poly, not strings): IDENTICAL — all 102 and 618 monomials
  match with coefficient ratio exactly 1 (no sign/scalar difference).
  Harness positive control (f31_sub2) reaches ~1e-7 (<1e-5) in this env.
  One harness fix: cmd_control resamples bases instead of crashing when the
  first deterministic base's lift diverges on every random start. NOTE: this
  confirms the artifacts are reproducible and self-consistent; it does NOT
  audit the derivation itself (system-selection, envelope bounds, Prop 4.3
  transcription) — those remain the T6 obligations, still open.
- FRESH-REPO SESSION (T2 f37 SUBCASE-2, 2026-07-20): harness validated on f37
  (control reaches ~1.7-2.1e-7 on 3/4 bases; one local-min at 5.7e-3 from a
  bad start — expected). Full-system stats: 50 lifting runs, min 1.233e-6,
  median 1.079e-3, NO run below the 1e-8 STOP threshold. Central-difference
  polish (noise floor ~1e-10) of the best four: DEAD FLAT plateaus
  (1.233e-6->1.168e-6; 1.487e-6->1.487e-6; 3.020e-6->3.020e-6;
  4.008e-6->3.972e-6) — genuine positive minima, ~4 orders of magnitude above
  the noise floor, NOT partial convergence. CONCLUSION (strong evidence, not
  proof): the f37 window system for subcase (2) is ALSO infeasible, matching
  f31. Both master-factor branches of subcase (2) now show genuine positive
  floors. NUANCE: f37's floor (~1e-6) sits ~20x below f31's (~2.3e-5); still
  solidly a true floor, but the smaller margin is worth remembering if an
  exact certificate is attempted (f37 is "closer" to feasible than f31).
  Best minima saved (f37_sub2_bestmin.pkl). Remaining for the case: subcase
  (1) windows for both factors (T3, audit the window bound first), base-region
  completeness (T4), exact certificate (T5), derivation audit (T6).
- FRESH-REPO SESSION (T6 PARTIAL AUDIT, 2026-07-20): pre-merge stability audit,
  full writeup in AUDIT.md. VERIFIED exactly: (i) sol4 solves G1; (ii) the
  d_-1=0 impossibility, both legs, symbolically; (iii) A=dm1*Ah, B=dm1*Bh with
  A,B each having a SINGLE dm2-factor (so regenerate's [0]-pick drops nothing),
  making the master-identity resultant chain a sound necessary condition;
  (iv) the f1 ODE solution and quartic separability; (v) Phi coeffs and f31/f37
  homogeneity (already in T1). CORROBORATED against arXiv:2204.14178: the
  (72,108)-is-last-under-125 framing and the "couldn't solve the system, left
  open" statement are in the paper nearly verbatim; the generator is validated
  against the paper's (7,21) case (its section 6). STILL UNAUDITED (could not
  extract Prop 4.3 / the (72,108) polygon from WebFetch — needs the PDF): the
  Newton-polygon subcase corners, the envelope bounds (esp. the hasty sub1
  deg<=15w), the equation-selection argument, and the [P,Q]=x^2 normalization.
  KEY RISK: if sub1's 15w upper bound is too small the sub1 window is invalid —
  do NOT trust sub1 numerics until the bound is proven. Base judged stable
  enough to merge as "verified reproduction + sound elimination + unproven-
  scaffold numerics," NOT as a closed case.
- LOCAL SESSION (T3+T6 AUDITS, 2026-07-20): obtained the FULL LaTeX source of
  arXiv:2204.14178 (e-print endpoint; no truncation). T3: Prop 4.3 ("Case
  (8,28)") transcription CONFIRMED VERBATIM — both subcases, all corners,
  [P,Q]=x^2; (8,28) confirmed as the open case (their sec. "seccion 4" closes
  the other (72,108) case, A0=(9,27)). Envelope bounds PROVEN for both
  subcases by redoing their calculo-de-C induction on our polygons:
  deg(d_{4-k}) <= 14k (sub2, via v_{-2,1}(C)<=0) / 15k (sub1, via
  v_{-1,1}(C)<=4 from the (0,8) corner), ord >= 12k (via v_{2,-1}(C)<=1);
  magic weights 14/15/12 all k-free; jetlift CONFIGS match exactly. AUDIT.md
  Risks 1 AND 2 retired — subcase-1 numerics now stand on a proven window.
  See T3_WINDOW_AUDIT.md. T6: equation-selection audit — verify_derivation.py
  (48 symbolic checks, all pass): commutator route reproduces the f1 ODE
  exactly (and [P,Q]=[C^2,F] gives v(F)=-5); lambda-isolation proven
  ((C^-1)_{-5} = -c3/C4^2 = 0 post-shift); D_k-polynomiality recursion closes
  with exact C4-exponent cancellation; cleared C-side slices == regenerate's
  D2(k)/D3(j) EXACTLY; dropped equations (k=8, j=4) only define dm12/lambda;
  deeper slices only define fresh dm14+. AUDIT.md Risk 3 reduced to two
  outline-only premises (GGV1 Prop 1.13/2.1 chain, alpha-strip WLOG) — see
  T6_SELECTION_AUDIT.md sec. 4. BONUS: the t=3 "universal relation" IS the
  paper's ecuacion principal (their closed case); their endgame
  (divisibility + degree count on tilde-d's) is the ready-made T5(a)
  template, worked example in hand. Unused ammunition for T5 noted:
  lambda-constancy (C4^28 | (D~^3)_{-4} with constant quotient) and the
  dm12 compatibility relation between the two dropped equations.
- LOCAL SESSION (T3 SUBCASE-1 NUMERICS, f31, 2026-07-21): with the 15w window
  bound now PROVEN (T3_WINDOW_AUDIT.md), ran the subcase-1 campaigns. Harness
  change (committed): JETLIFT_SEED env var so independent chunked stats runs
  sample distinct bases (the fixed seed 2026 would repeat runs across fresh
  processes), per-run checkpointing of best_*.pkl, and a `merge` mode to
  combine per-seed chunks for polish. CONTROLS PASS on the wider windows for
  BOTH factors (f31_sub1: 4/4 at ~1.2-2.7e-7; f37_sub1: 4/4 at ~1.9e-7-2.2e-6;
  threshold 1e-5) — the optimizer finds solutions on the 46-unknown sub1
  windows when they exist. f31_sub1 FULL CAMPAIGN: 57 independent lifting runs
  (8 seeds), min 8.62e-6, median 1.4e-3, NO run below the 1e-8 STOP threshold.
  Central-difference polish (noise floor ~1e-10) of the best four: DEAD FLAT
  (8.622e-6->8.334e-6; 1.803e-5 flat; 2.795e-5->2.792e-5; 3.361e-5->3.342e-5)
  — genuine positive floors, same signature as both sub2 branches.
  CONCLUSION (strong evidence, not proof): f31 subcase-1 window system is
  INFEASIBLE. Floor (~8e-6) sits below f31_sub2's (~2.3e-5), consistent with
  the wider window being "closer" to feasible. Best minima: f31_sub1_bestmin.pkl.
  f37_sub1 stats in progress (early sample: 6 runs, min 3.6e-6, all positive).
- LOCAL SESSION (T3 SUBCASE-1 NUMERICS, f37, 2026-07-21, same day): f37_sub1
  FULL CAMPAIGN: 29 independent lifting runs (10 seeds), min 4.742e-7, median
  7.9e-4, NO run below the 1e-8 STOP threshold. Central-difference polish
  (noise floor ~1e-10) of the best four: DEAD FLAT (4.742e-7->4.704e-7 with a
  rigid 3-iterate tail; 2.691e-6, 3.599e-6, 3.627e-6 all exactly flat).
  CONCLUSION (strong evidence, not proof): f37 subcase-1 window system is
  ALSO INFEASIBLE. ALL FOUR factor x subcase cells now show genuine positive
  floors — the numerical phase of T3 is COMPLETE. CAVEAT for T5: 4.7e-7 is
  the tightest floor in the project (~3 orders above noise, vs ~4 elsewhere;
  ~50x above the STOP threshold); the floor ordering
  f37_sub1 (5e-7) < f37_sub2 (1.2e-6) < f31_sub1 (8e-6) < f31_sub2 (2.3e-5)
  says wider windows + the 618-term factor sit "closest" to feasible — the
  exact certificate should attack f31_sub2 first (largest margin) and treat
  f37_sub1 as the hard case. Best minima: f37_sub1_bestmin.pkl. Remaining:
  base-region completeness (T4), exact certificate (T5), two GGV1-premise
  outlines (T6, see T6_SELECTION_AUDIT.md sec. 4).
- FRESH-REPO SESSION (SUBCASE-1 NUMERICS + T5 SETUP, 2026-07-21): with the
  window bounds now proven (Fable T3), ran subcase-1 for both factors on the
  proven [12w,15w] windows. Controls pass (~1-3e-7). f31_sub1: 57 runs, min
  1.70e-5, median 1.84e-3; central-diff polish DEAD FLAT (1.70e-5, 1.71e-5,
  2.37e-5, 3.28e-5) -> genuine positive floor. f37_sub1: 7 runs only (slow,
  ~170s/run; SMALL SAMPLE), min 1.90e-4, median 6.3e-3; polish plateaus
  positive (8.5e-5, 8.4e-4, 1.24e-3, 4.2e-3). So ALL FOUR factor x subcase
  branches now show genuine positive floors under central-diff polish:
  f31_sub2 ~2.3e-5, f37_sub2 ~1.2e-6, f31_sub1 ~1.7e-5, f37_sub1 ~8.5e-5.
  The (72,108) case is numerically infeasible everywhere on proven windows
  (strong evidence, NOT proof). Caveat: f37_sub1 has thin statistics (7 runs);
  a longer run would firm it up. Best minima saved (*_sub1_bestmin.pkl).
  T5 EXACT CERTIFICATE (exploration, T5_NOTES.md): ruled out the "universal
  cubic" as an ideal member (Singular: not in <G1,G2,G3,G5+Phi>); anchor
  slices are irreducible deg-31/37 hypersurfaces (no single-slice kill). Set
  up the (y+1)-adic endgame concretely using the PROVEN ord bound: strip
  d_k = y^{12w} d_k~ so f31(d,Phi)=y^1500 f31(d~,Phi~) with d~ of degree
  <=2w (sub2), Phi~ = -(y+1)^30(quartic)/6630; verified f31|_{Phi=0} =
  dm1^21 * h31 (h31 weight-20 irreducible) and f37|_{Phi=0} = dm1^18 * h37.
  Mod (y+1)^30 the identity collapses to dm1^21 h31(d~), launching a
  paper-style (y+1)-adic order/degree count (a=v_{y+1}(d~_{-1})=0 case forces
  (y+1)^30 | h31(d~), a deg<=40 object). NOT carried to a contradiction; this
  is the recommended next work. Also reconciled AUDIT.md with the T3/T6
  resolutions and independently re-ran verify_derivation.py (48 checks pass).
- MERGE RECONCILIATION (2026-07-21): the two 2026-07-21 entries above are
  INDEPENDENT campaigns (different machines, seeds, and run counts) that ran
  the same T3 subcase-1 task in parallel and agree on every verdict — a
  genuine replication of the sub1 infeasibility for both factors. The local
  29-run f37_sub1 campaign (min 4.742e-7, polish flat) supersedes the
  fresh-repo session's 7-run sample and resolves its "thin statistics"
  caveat. Artifact naming after merge: fresh-repo files keep canonical names
  (f31_sub1_control.log, *_sub1_bestmin.pkl, *_sub1_polish.log,
  *_sub1_stats.log); the local replication's copies carry a _local suffix;
  the local per-seed logs are f31_sub1_stats_s*.log / f37_sub1_stats_s*.log.
  Tightest known floor is now f37_sub1 at 4.7e-7 (local campaign) — the
  T5 priority ordering in the local entry stands.
- LOCAL SESSION (T5 BY-HAND CAMPAIGN, 2026-07-21, same day): the (y+1)-adic
  Newton-polygon program (T5_NP.md) plus a parallel-agent sweep produced the
  first PROVEN infeasibility strata for f31 subcase-2. Chain: (1) graded
  decomposition verified (verify_graded.py); (2) T5_NP.md: a_t=0 case tree
  (three branches; degenerate escape killed by the new h7/h6/h5 cascade
  Lemma 1), g-cascade exact reformulation (251 eqs as 7x30-blocks + terminal,
  all objects deg<=40); (3) T5_MULTIPLACE.md (agent-derived, INDEPENDENTLY
  RE-VERIFIED here: t5_multiplace_verify.py, 8 exact check groups pass): the
  unit cofactor of Phi~ IS the quartic (u = -q/6630), so the cascade runs at
  the q-adic place simultaneously; the joint kill PROVES strata
  (a_t,a_q) in {(0,0),(1,0),(2,0),(3,0),(4,0),(0,1)} infeasible — including
  every branch of the a_t=0 tree with q^2 not dividing d~_-1 — and the
  sigma-locus {d1~=0, 4d0~=d2~^2} dies in ALL strata via an explicit unit
  equation 512(Phi~+d2~e^3)^4(4Phi~-5d2~e^3) = 6561e^17 (Mason-Stothers
  margins 34>8, 17>6; M-S application LINE-AUDITED this session: case split
  (i)-(iv) exhaustive, both M-S triples pairwise coprime and nonconstant,
  type-A/B prime classification exhaustive, all integrality kills recomputed
  — Theorem 2 is SOUND).
- LOCAL SESSION (T5 STRATUM KILLS, 2026-07-21, continued): strata (5,0) and
  (1,1) PROVEN infeasible (T5_STRATA_50_11.md, independently re-verified:
  t5_strata50_11_verify.py, 8 exact check groups pass; level-5 squeeze
  h5|_{d1=0} = -9216 d2 sigma^2 + 2048 e^2 re-derived by hand here). New
  transportable tools: the level-5 squeeze (terminal relation absorbs the
  sigma^2 term, coprimality forces ê^2 | ĝ, so ê constant and sigma
  constant) and infinity-place domination (the f=6 term Phi~^6 e^3(-3072 s^2)
  has degree 219 > 210 >= every other term when deg e = 5). Survivor map now
  13 strata: a_q=0: {6,7,8,9,10}; a_q=1: {2,3,4,5,6}; a_q=2: {0,1,2}.
  15 strata survive, precisely mapped (most rigid: (5,0),(1,1)).
  Computational flank: naive/dehomogenized/weighted Grobner on blocks all
  exceed 90-min caps (Singular 4.2.1); NulLA linear-certificate route is a
  PROVEN no-go for block 1 (LINALG_CERT_REPORT.md: constant-window witness
  family, dim>=3, so block-2 coupling is mandatory; bonus: all block-1 eqs
  are BIhomogeneous under (weight, subscript-index)); msolve 0.10.1 built in
  WSL (12-thread F4) and flint-based per-base pivot lift in progress.
  Strategy docs: NEXT_CASES.md (paper-verified: (75,125) is the next case,
  24 candidates in [125,150]; corner (8,28) reappears at (108,144) so a
  corner-level kill covers both). Remaining for f31_sub2: the 15 strata;
  then f37 (non-uniform grading), sub1 (caps 60-6f), char-0 certification.
- FIELD-SPLIT REPAIR + NEW EXACT RESULTS (2026-07-22): accepted the
  coefficient-field audit: the old scalar `a_q` ledger is conditional because
  `q` splits geometrically. Added `FIELD_SPLIT_AUDIT.md` and source-linked
  proof checks. Built `split_place_ledger.py`: 327 geometric multiplicity
  vectors = 21 uniform + 306 partial-support; exact terminal valuation/degree
  pruning kills 81 complete strata (75 genuinely partial-support). After
  correctly scoped old proofs and the new `a_t=9` T2 kill, the honest f31
  subcase-(2) frontier is 235 strata / 420 live T1/T2 branches, with T3 dead
  globally. The geometrically q-coprime `a_t=7` proof is checked directly
  against `f31_graded.txt`, including the sigma=0 parity edge and all 36
  noncancelling infinity cases. NEW f37 RESULT: restricting the regenerated
  pre-resultant equations to the bare-resultant free family `d2=d1=0` gives
  `12rs(r^2-es)=e^5` and `3e(r^2+es)=2Phi`; hence `e|Phi`, local
  valuations plus infinity force `e=Ct^10`, and the final required t-order
  10 exceeds a residual degree cap 4. Thus the entire f37 free family is
  PROVEN not to lift to the original system.
- A_T=9 T1 REDUCTION (2026-07-22): `T5_90_T1.md` and the source-linked
  `t5_90t1_verify.py` reduce the remaining geometrically q-coprime uniform
  T1 branch without claiming closure. If the residual `E` is linear, exact
  level-6/5 local orders force `H=s^5 W^2`, `d1=s^4 W`, `sigma=s^3 Sigma`,
  `G6=s^2 K`, `deg W<=2`, followed by the perfect-square condition
  `s | gamma^4 t^6 + 2c delta q W`. If `E=gamma` is constant, 43 of 50
  infinity degree cells die uniquely; the tail identity kills the
  `deg sigma=5` tie, degree 239 kills `deg sigma=4`, and the two intervening
  square constraints fix the leading two coefficients in every survivor:
  `d1=-gamma^4/(4096c)*(y^2+25y/4+v0)`, with only `sigma=0` or
  `deg sigma<=3` left. Thus 45/50 constant cells are excluded and the five
- A_T=9 T1 CONSTANT CELL CLOSED (2026-07-22): exact coefficient convolution
  in `t5_90t1_constant_verify.py`, sourced from `f31_graded.txt`, continues
  the infinity reduction from degree 238 down to 226. Alternating square
  coefficients force `v0=525/32`, then `sigma=0`, and determine all of `d2`:
  `d2=(c/gamma^3)(-95200y^4+255850y^3-513451y^2
  -10656467y/8+132899897/8)`. Degrees 229 and 227 cancel exactly after the
  forced substitutions, so the checker includes the newly active f2/f1/f0
  terms rather than truncating them. Degree 226 is the nonzero constant
  `29570349989420274657771126784*c^5*gamma^8`, proving the complete
  constant-E branch infeasible. The only geometrically q-coprime a_t=9 T1
  possibility left is the nonconstant local shape `E=gamma*s`,
  `H=eta*s^5W^2`, `d1=delta*s^4W`, `sigma=s^3 Sigma`, with
  `s | gamma^4 t^6+2c delta qW`.
  residual cells have an explicit low-parameter form.
- A_T=9 T1 NONCONSTANT LOCAL DESCENT (2026-07-22): the only remaining
  geometrically q-coprime cell now has a compact field-independent level-5/4
  certificate in `t5_90t1_local_verify.py`. With `x=s` and
  `R=gamma^4*t^6+2c*delta*qW`, successive necessary coefficients are nonzero
  unit multiples of `S0^2`, `R1^2`, `S1^2`, `R2^2`, `S2^2`, where
  `sigma=x^3(S0+S1*x+S2*x^2+...)` and `R=R0+R1*x+R2*x^2+...` with `R0=0`
  already. Thus every survivor satisfies `s^6|sigma` and `s^3|R`. Since
  `deg W<=2`, the latter fixes all coefficients of `W`. This is a reduction,
  not yet a closure: the next exact target is the level-3 order-two balance,
  which determines `d2 mod s`; lower levels must then be continued or paired
  with infinity. The checker derives the h4 identity from `f31_graded.txt`
  and verifies all five compact local coefficient identities symbolically.
- LOWER-CASCADE ENGINE PIVOT (2026-07-22): reconciled the post-repair
  frontier with an external proof-program audit. The 235 surviving strata
  and 420 live T1/T2 records are now treated as lattice points in one
  valuation problem, not as hundreds of bespoke proof obligations. The
  active architecture is documented in `CASCADE_ENGINE_PLAN.md`: extract a
  source-linked cascade signature, propagate exact tropical transitions
  through levels 6/5/4, attach residue equations to every tied minimum, and
  couple the four split places with global degree budgets and infinity. The
  nonconstant `a_t=9` T1 cell remains the first regression case; a bounded
  continuation through levels 3 and 2 showed further alternating square
  conditions that determine successive coefficients of `d2`, but did not
  produce an immediate contradiction. Therefore bulk cone elimination is
  now primary, while cell-level descent continues only where it validates or
  closes an exceptional residue transition.
- CASCADE ENGINE PHASE B (2026-07-22): `cascade_engine.py` implements the
  exact q-place descent (levels 7/6/5/4, ultrametric semantics with
  recorded residue obligations, global zero flags, budget-coupled
  four-place join; soundness discussion in `CASCADE_ENGINE_REPORT.md`).
  Regressions (`test_cascade_engine.py`, in `run_tests.sh`): terminal layer
  reproduces `split_place_ledger.json` on all 654 records; hand-worked
  level-6 descent example; depth monotonicity. RESULT: of the 420 open
  branches, depth 5 kills 352 and depth 4 kills **390, leaving 30
  survivors** — every stratum with `a_t<=4` dies; survivors sit at
  `a_t in [5,10]` with `b_max<=3`, no `b_i=2`, and all carry explicit
  residue obligations (`cascade_cones.json`). Two kills re-derived fully by
  hand ((0,(2,1,1,1),T2) at level 5; (0,(2,2,1,1),T1) at level 6): the
  dominant pattern is terminal-budget rigidity pinning g-valuations at
  high-b places below what the next level can produce. STATUS: kills are
  ENGINE-PROVEN, PENDING INDEPENDENT AUDIT — the required next steps are an
  independent checker not sharing the engine's pruning code, then
  cone-lemma extraction over the 390 row kills, then Phase C residue
  systems for the 30 survivors (T2 cases first; the T1 tail pairs with the
  `a=9` local-descent machinery and infinity).
- CORNER-144 COMPARISON (2026-07-22, conditional forcing result): for the
  GGV5 length-one chain `(8,28)->(7/4,3)`, table pair `(3,4)`, the standard
  root-shift/Laurent-chart continuation has `[P,Q]=x^2`, `P=C^3`,
  `Q=C^4+(commuting C-powers)+F`, `v(F)=-9`, `ell(C)=x^4 C4` with
  `C4=y^3(y+1)h4`. The forcing family is `12 C4 f'-21 C4' f=C4^2`,
  `f=C4^4 F_-9`. Type-I squarefreeness forces `h4=y^4-y^3+y^2-y+1` (the
  10th cyclotomic), `f=-y^4(y^5+1)^2/15`, and
  `Phi_144=f C4^67=-y^205(y^5+1)^69/15`, hence
  `(deg,ord,mult_{y+1},cofactor degree)=(550,205,69,276)` vs the current
  `(238,204,30,4)`. The operator and D-exponents fit the parametric
  formulas (`D_k=C_k C4^(11-3k)`), but the forcing-divisor/cascade
  signature does not transfer numerically. NOTE the geometric skeleton DOES
  persist: the cofactor is `h4^69` — again one t-place plus one quartic of
  split places, with q-multiplicity 69 instead of 1. Candidate target
  windows are upper `22w/23w`, lower `8w+ceil(w/5)` (quasipolynomial, not
  affine) pending a full polygon reduction. Exact checks:
  `corner144_verify.py` (51/51 pass). Verdict: partial support for a
  `(power pair, final multiplicity)` template; against an `A0`-only
  cascade template. Details and the conditional boundary:
  `CORNER_144_COMPARISON.md`.
- CONE-LEMMA EXTRACTION (2026-07-22): `cone_lemmas.py` compresses the 390
  depth-4 engine kills into two certificate families and verifies the
  certificate verdict against the engine verdict on all 420 open branches
  (`test_cone_lemmas.py`, in `run_tests.sh`). Family L (single-place
  kills, 6717 of 7200 flag cases): a place with `v_p(e)=beta` admits no
  consistent local chain — unconditionally (all zero-flag cases) exactly
  when `beta=2` or `beta>=4`, so the live local values are
  `beta in {0,1,3}` (single exception: T2 a=0 also kills beta=3). This
  single arithmetic law explains why no surviving vector contains 2 and
  survivors have `b_max<=3`. Family B (483 cases): one budget dimension
  overflows — only five patterns occur (T1:d1 181, T1:g7 150, T2:sigma 96,
  T2:g6 49, T2:d2 7). Row-level certificates:
  `cascade_cone_certificates.json`; human statement:
  `CASCADE_CONE_LEMMAS.md`. NOTE: certificates are generated from the
  engine's own tables (compression, not independence); the independent
  audit is a separately authored checker in progress.
- INDEPENDENT KILL AUDIT PASSED (2026-07-22): `audit_cascade_kills.py`,
  authored independently (Codex) from the semantics specification in
  `CASCADE_ENGINE_REPORT.md` only — no access to `cascade_engine.py`, own
  parser for `f31_graded.txt`, own homogeneity self-checks, conservative
  relaxed-tie semantics — agrees with the engine on ALL records: terminal
  654/654, open branches 420/420 (390 kills, 30 survivors), zero
  disagreements, ~2 min runtime (re-run locally before commit; only
  edit made here is exit-code plumbing for the test suite). Wired into
  `run_tests.sh` (`--quiet`). The 390 depth-4 kills are therefore now
  ENGINE-PROVEN AND INDEPENDENTLY AUDITED; the honest f31 subcase-(2)
  frontier is 30 branches at `a in [5,10]`, `b_i in {0,1,3}`, each with
  explicit residue obligations (Phase C targets). Method notes:
  `CASCADE_KILL_AUDIT.md`.
- F37 BRANCH CLOSED — RESULTANT ARTIFACT (2026-07-22): the saturation
  program of `F37_FRONTIER.md` terminated in one stroke with a STRONGER
  statement: `f31 IS IN the pre-resultant ideal`
  `<G1,G2,G3,G5body+Phi>` over `Q[d~2,d~1,d~0,e,r,s,dm4,Phi]` with `Phi`
  a free indeterminate. Hence the elimination ideal in
  `(d~2,d~1,d~0,dm1,Phi)` is EXACTLY `<f31>` (Singular: principal, deg 31,
  102 terms, equal to f31 both ways), and the resultant's extra factors
  `f37` and `dm1^21` are classical excess — pure artifacts, no geometry.
  Every solution of the original system satisfies `f31=0`; the locus
  `{f37=0}\{f31=0}` contains no solutions, uniformly (subsumes the
  free-family kill and every f37 ledger stratum; the f37 numerics are
  moot). Triple verification: G-system Groebner facts; explicit `lift()`
  cofactor certificate re-expanded EXACTLY in sympy with no Groebner
  trust (`f37_sat_verify.py`, PASS, re-run before commit); independent
  H-system route `f31 in <H2,H3,H5>:dm1^inf`. Report:
  `F37_SATURATION_REPORT.md`; reproducible script `f37_sat_confirm.sing`
  (WSL note: /mnt/c broken, pipe via stdin). SCOPE: the membership is a
  polynomial identity with `Phi` free, so it specializes to ANY window
  instance — it applies to subcase (1) as well as (2) (the subcases only
  restrict which `d~k(y)` are admissible, which is downstream of the
  ideal identity; flag this reading for the writeup audit). REMAINING for
  (72,108): the f31 branch alone — 30 audited subcase-(2) survivors
  (Phase C/D) and the subcase-(1) windows (re-target the cascade engine
  with sub1 caps; the graded h_f identities are window-independent).
- T-PLACE COUPLING (2026-07-22): `cascade_engine.py` now couples the place
  `t=y+1` as a fifth place (`--with-t`; `t_place_profiles`). At t the
  identity reads `v + s_{l+1} = ultrametric(s_l, w_l)` with `v=30-3a`,
  `v_t(ehat)=v_t(u)=0` (q(-1)=3315), and the e-slot of every h_l monomial
  costing `a`; terminals are `s_7=2v_t(d1)` (T1) / `s_6=2v_t(sigma)` (T2).
  The descent generalizes the q-place case by two shift parameters; the
  audited q-only path is byte-identical (regressions re-run). RESULT
  (`cascade_cones_qt.json`): no additional tropical kills (deep t-place
  cancellations are always available to zero-budget states — the t
  constraint is a residue-level constraint, as the by-hand a=9 proofs
  already showed), but 88 of the 320 surviving flag-cases die and every
  remaining witness now carries explicit t-cancellation obligations with
  finite depths (max 18). Frontier: 30 branches / 232 flag-cases, each
  with a concrete obligation list = the Phase C worklist. New regression:
  t-terminal law, forced cancellations, place monotonicity.
- SUBCASE-(1) LEDGER + PARAMETER CORRECTION (2026-07-22): sub1 cascade
  parameters derived, CHECKED, and one corrected (`sub1_cascade_verify.py`,
  in `run_tests.sh`): stripped caps (d2,d1,d0,e)=(6,9,12,15), sigma<=12,
  deg h_f(d~)<=60-6f, budget a+sum(b)<=15, standard regime a<=10 (v>=0).
  CORRECTION (found by the independent Codex derivation): the proposed
  uniform cap `deg g_l <= 15+3a` is UNSUPPORTED for sub1 — the sub2
  induction closes only because 3(10-a)+(10+3a)=40 equals the window
  bound exactly; for sub1, 3(15-a)=45-3a exceeds v and the bound grows
  ~15/level, while the top anchor needs an unavailable LOWER bound on
  deg(ehat). The sub2 cap 10+3a remains SOUND (bottom-anchored at
  g_1=h_0/t^v, closes exactly — re-derived here to confirm the audited
  kills are unaffected). Sub1 ledger uses safe terminal caps
  deg g7<=46 (T1) / deg g6<=48 (T2). T3/SIGMA-LOCUS TRANSFERS to sub1:
  nonconstant degree triples (10,34,34) and (15,51,51), all
  Mason-Stothers margins re-verified (51>10, 34>8, 20>6/5/4, 17>5,
  17!|150). LEDGER (`split_place_ledger_sub1.py/.json/.md`): 1333 strata,
  136 terminal kills (all partial-support), 1197 open standard-regime,
  26 alternate-regime strata (a in [11,15], v<0, need a separate
  reduction). NEXT: parameterize the engine budgets per level (sub1 has
  no uniform g-cap) and run the level descent + t-place on the sub1
  ledger; expect the q-place descent to do the heavy lifting since sub1
  terminal budgets are weak.
- ALTERNATE REGIME FLIPPED CASCADE (2026-07-22): for the 26 sub1 strata
  with a in [11,15] (v=30-3a<0), the minimum t-order of the graded
  identity sits at f=7 and 21a+7|v| telescopes to the CONSTANT 210:
  `F = t^210 * G'` with `G' = sum_f t^((7-f)|v|) u^f ehat^(21-3f) h_f`,
  verified exactly on a random a=12 window instance. The descending
  cascade anchored at h_7=8192*d1^2 replaces the upward one; the q-place
  transitions, terminal identities, and terminal degree caps survive
  verbatim (t-powers are units at q-roots regardless of sign(v)).
  First-level parity/degree arguments kill 19/52 branches (T1 7/26,
  T2 12/26; 6 strata fully dead). Remaining: 33 branches / 20 strata
  needing the lower descending cascade or residue analysis. Checks:
  `alt_regime_verify.py` (7 groups, in `run_tests.sh`); derivation and
  per-stratum verdicts: `ALT_REGIME.md`.
- T6 PREMISES DISCHARGED (2026-07-22): the two outline-only premises of
  `T6_SELECTION_AUDIT.md` sec. 4 (AUDIT.md Risk 3) are now full
  arguments, both READY-WITH-CITATION (`T6_PREMISES.md`;
  `t6_premises_verify.py`, 14 exact checks, in `run_tests.sh`). Key
  finding: GGHV22 sec. 4 runs the identical setup machinery in complete
  published detail for the twin CLOSED case (9,27) — our premises are
  that argument transported t=3->t=4 (numerology v(P):6->8, v(Q):9->12,
  v(F):-4->-5); the paper leaves (72,108) open only at the final
  elimination. GGV1 source fetched (`paper_src/1401.1784_GGV1.tex`);
  Prop 1.13 (bracket-valuation inequality) and Prop 2.1 (aligned
  elements are common powers) located verbatim and used as stated.
  Premise 1 (l(P)=R^2, l(Q)=R^3, R=x^4C4, C4=y^7(y+1)): strict gap
  2<19 triggers Prop 1.13, Prop 2.1 gives R, corners force the C4 shape,
  primitivity gcd(4,7,1)=1, scaling normalizes; the lambda_Q=1 gauge-fix
  matches the paper's own (9,27) usage. Premise 2 (alpha-strip WLOG,
  v(F)=-5): GGHV22 lines 1508-1546 shifted t=3->t=4; strip terminates
  since -5 is not 0 mod 4; the Remark shift absorbs alpha1*C preserving
  [P,Q]=x^2. Derivation-audit debt for the writeup is now ZERO modulo
  published GGV1 propositions used as stated.
- T2 COLUMN SQUEEZE (2026-07-22): the level-5 squeeze applied to the
  12-cell T2 survivor column (`T5_T2_COLUMN.md`;
  `t5_t2_column_verify.py`, in `run_tests.sh`). PROVEN infeasible:
  (5,(1,0,0,0)), (6,(1,0,0,0)), (6,(1,1,0,0)), (6,(1,1,1,0)) — 4 cells
  / 8 flag cases. IMPORTANT CORRECTION to PHASE_C_WORKLIST.md's
  [judgment] that all 12 cells meet the squeeze hypotheses: under
  partial q-support the coprimality conclusion weakens to F^2 | G with
  F only the q-COPRIME remainder of e, and the a>=7 cells hit genuine
  infinity-degree ties (Phase D). Open with exact residuals: the four
  a=7 cells, three a=8 cells, and (9,(1,0,0,0)) (its g5=0 case narrowed
  to (deg F,deg Z,deg G; deg e,deg sigma)=(0,6,12;10,8)). The sub2
  frontier is now 26 cells; the remaining T2 work is
  infinity-domination (Phase D), matching the T1 tail's needs.
- SUB1 DEPTH-5 SWEEP (2026-07-22): first exact bulk result on subcase
  (1). After the engine memory hardening (staged Pareto compaction,
  restored profile cache, RSS guard, per-a-chunk checkpoints — the
  unguarded first attempt saturated RAM when orphaned by a terminal
  crash), the guarded sweep over the 1089 standard-regime open strata
  finished cleanly (no resource skips): **2178 branches -> 1606 killed,
  572 survivors** (`cascade_cones_sub1_depth5.json`). Survivor
  structure: a stable ~56-branch family per a for a<=7 (thinning to
  51/42/31 at a=8/9/10), T2-heavy (381 vs 191), and — unlike sub2 —
  b=2 and b in {4,5} appear among survivors: sub1's looser caps weaken
  the single-place kills, as predicted after the cap correction.
  STATUS: engine-proven, PENDING independent audit (the existing audit
  checker is sub2-specific; a sub1-capable port is required before
  promotion). NEXT: depth-4 sweep, then t-place, then cone-lemma
  compression and the audit port.
- SUB1 DEPTH-4 SWEEP (2026-07-22): levels 7/6/5/4 over the same 2178
  open branches: **1899 killed, 279 survivors**
  (`cascade_cones_sub1_depth4.json`, guarded run, no resource skips).
  Survivor anatomy is strikingly uniform: EXACTLY 26 branches per a for
  every a<=8 (24 at a=9, 21 at a=10) — a single a-independent family,
  ideal for cone-lemma compression. Level 4 hit T2 hardest (split now
  T1 180 / T2 99); b=4 no longer appears among survivors
  (b in {0,1,2,3,5} remain; b=2 down from 221 to 33 records). Pending
  the same sub1 audit port. NEXT: t-place coupling on the 279, then
  cone lemmas + audit.
- SUB1 T-PLACE COUPLING (2026-07-22): q+t depth-4 sweep
  (`cascade_cones_sub1_qt.json`, guarded, clean): same 279 surviving
  branches (t adds no tropical kills, matching sub2 — the t-constraint
  is residue-level), but 266 of 2519 flag-cases die (-> 2253) and every
  remaining witness carries explicit t-cancellation obligations. Sub1's
  endgame is now the same typed residue worklist shape as sub2's.
- SUB1 CONE LEMMAS (2026-07-22): `cone_lemmas.py --window sub1`
  compresses all 1899 sub1 depth-4 kills into the SAME two certificate
  families, verdicts matching the engine on all 2178 branches
  (`cascade_cone_certificates_sub1.json`, `CASCADE_CONE_LEMMAS_SUB1.md`).
  Family L: unconditional dead beta = {4} u [6, 15-a] for T1, plus {2}
  for T2 — live local values {0,1,2,3,5} (T1) / {0,1,3,5} (T2), the
  sub1 analog of sub2's {0,1,3} law (143 unconditional pairs). Family
  B: only TWO budget patterns (T1:d1, T2:sigma). No coupled
  certificates needed in either subcase. The md template is now
  window-aware; sub2 doc regenerated identically (regression PASS).
- SUB1 KILL AUDIT PASSED (2026-07-22): `audit_cascade_kills_sub1.py`
  (Codex-authored port of the spec-only checker, zero references to
  cascade_engine.py, own cap derivation) agrees on ALL records: terminal
  2614/2614 standard-regime branches, depth-4 2178/2178 (1899 kills +
  279 survivors confirmed, none undecided), depth-5 kills a subset of
  depth-4, zero cap discrepancies, ~2.5 min (re-run locally before
  commit). Wired into `run_tests.sh` (--quiet). The sub1 standard-regime
  frontier is therefore ENGINE-PROVEN AND INDEPENDENTLY AUDITED at 279
  branches (one a-independent 26-family), with cone-lemma compression
  and t-obligations already in place. Method notes:
  `CASCADE_KILL_AUDIT_SUB1.md`.
- T2 INFINITY CONVOLUTION (2026-07-22): exact top-coefficient
  convolution on the 8 open T2 cells (`T5_T2_INFINITY.md`;
  `t5_t2_infinity_verify.py`, in `run_tests.sh`): 0 cells fully killed
  but ALL 8 strictly narrowed — the a9 b1000 g5=0 flag case is CLOSED
  and 11 residual degree-states are removed (2 by pattern A across all
  flags, 9 d2=0 fixed-F states). Confirms the a>=7 T2 ties need the
  systematic Phase D infinity layer rather than per-cell convolution.
- ALTERNATE REGIME LEVELS 2/3 (2026-07-22): the descending flipped
  cascade pushed through its h6/h5 steps on the 33 open branches
  (`ALT_REGIME_L2.md`; `alt_regime_l2_verify.py`, checks A-E, in
  `run_tests.sh`): 6 new kills -> **27 residual branches**, with a
  terminal UFD residual normal form and survivor caps established for
  what remains. NOTE: the codex runtime lost its job registry before
  posting the final report; the deliverables were verified locally
  (all checks pass) and landed on that basis.
- SUB1 PHASE C WORKLIST + LIBRARY OVERLAP (2026-07-22): the sub1
  obligation inventory (`PHASE_C_WORKLIST_SUB1.md`,
  `phase_c_inventory_sub1.py`): 67 distinct patterns over 2253 cases /
  21337 obligations (all 26965 tied strings exact h_l members). THE KEY
  RESULT: sub2's 41 patterns are a SUBSET of sub1's — 29 identical, 12
  same-key deeper, and of the 26 new keys, 20 reuse sub2 tied-sets;
  only 6 carry genuinely new residue polynomials, ALL sub-cuts of h_4.
  91% of sub1's patterns need no algebra beyond the sub2 library. The
  26-family is a-parametric with a-INDEPENDENT q-obligations and an
  AFFINE t-depth law (term_cancellation depth = v = 30-3a exactly,
  verified a=0..9): residue lemmas can be proven once and quantified
  over a. The level-5 squeeze covers all 99 sub1 T2 cells. Cheapest
  meaningful targets mirror sub2's: (a9/a10,(1,0,0,0),T2),
  ((10),(0,0,0,0),T1), (a9,(1,0,0,0),T1). Zero-obligation
  (a10,(1,1,1,1)) tier is pure Phase D.
- T-PLACE AUDIT PASSED + PROOF INVENTORY (2026-07-22): the flag-case
  claims are now tier-1: `audit_tplace_cases.py` (spec-only, no engine
  access, re-run locally) agrees on ALL cases — sub2 7872/7872, sub1
  41592/41592, branch-consistency q vs q+t exact in both windows; wired
  into `run_tests.sh`. `PROOF_INVENTORY.md` added: 41-row claim graph
  (C0-C40) with trust tiers (10 tier-1, 13 tier-2, 8 tier-2* unwired
  checkers, 9 tier-3); every load-bearing kill on the frontier path is
  tier 1-2. Its inconsistency sweep triggered hygiene fixes committed
  here: correction banners on PHASE_C_WORKLIST.md (squeeze judgment
  falsified by execution), F37_FRONTIER.md (superseded by saturation),
  HANDOFF.md (stale frontier). REMAINING WRITEUP GAPS (from the
  inventory): (1) the frontier itself; (2) envelope bounds need a
  committed referee-grade checker (currently T3_WINDOW_AUDIT prose);
  (3) per-cell kills (a=9, T2 column/infinity) and the alternate regime
  are same-author checker-only — queue them for the next audit round
  alongside the infinity layer.
- RESIDUE LEMMA LIBRARY (2026-07-22): the shared tie-rise residue
  systems are now proven lemmas (`RESIDUE_LEMMAS.md`;
  `residue_lemmas_verify.py`, in `run_tests.sh`): of 23 source-derived
  tied-set shapes, 21 are CONSTRAINT lemmas (rational torus solutions
  exist — they pin leading coefficients to exact hypersurfaces for
  Phase D) and TWO are KILLS over the splitting field: **C08**
  (h_5|sigma=0 three-term cut; hits 3 sub2 + 54 sub1 q-rise cells) and
  **C20** (a new sub1 h_4 three-term cut; 8 sub1 cells). Singleton
  cuts (pure squares, e.g. h_6|d1=0) forbid rises outright. P6/P10/P11
  are constraints, incidence mapped per window. NEXT: feed C08/C20 and
  the singleton prohibitions back into the engine as forbidden-rise
  rules (w_options currently over-allows those ties) and re-sweep —
  every flag case whose survival routed through a killed rise dies.
- LEMMA-ASSISTED RE-SWEEP (2026-07-22): the engine now consumes the
  proven kill lemmas as forbidden-rise rules (`FORBIDDEN_RISES` +
  `--residue-kills`, default OFF so audited artifacts stay
  reproducible; caches key on the flag). Result
  (`cascade_cones_qt_rl.json`, `cascade_cones_sub1_qt_rl.json`): no
  whole-branch kills; flag-cases sub2 232->228, sub1 2253->2170 (87
  eliminated where survival routed through a C08/C20 tie rise).
  HONEST READING: incremental — most survival routes use
  term_cancellation (free-g) obligations, which are Phase D territory,
  exactly as the worklists filed them. The decisive next layer remains
  the infinity/max-plus build. Status of the 87 eliminations:
  engine+lemma-proven (residue_lemmas_verify is tier 2); t-place audit
  extension for the rl artifacts queued with the infinity round.
- ALTERNATE-REGIME AUDIT + ENVELOPE CHECKER (2026-07-22, session
  close): (1) cross-author independent audit of the alternate-regime
  chain (`ALT_REGIME_AUDIT.md`, `alt_regime_audit_verify.py`, 9 check
  groups, in `run_tests.sh`): ALL CONFIRMED — flipped reduction
  re-derived at a=13/14, closed telescope identity
  `G' = T^7(E^21 h_0 + u r_0)`, q-invariance proven regime-independent,
  kill list 25/25 reproduced, residual list honest (probes add zero
  kills). Two notes fixed here: the L2 checker's projection quantifier
  all->any (cones verified identical, conceptual fix), and the stale
  6->10 dead-strata figure banner. The alternate regime is now TIER 1.
  (2) `envelope_bounds_verify.py` (+`ENVELOPE_BOUNDS.md`): the
  envelope-bound induction is now a committed referee-grade checker,
  106/106 checks, w=0..5, named assumptions A1-A3 isolated; found a
  SIGN TYPO in GGHV22 lines 1462-1466 (and mirrored in
  T3_WINDOW_AUDIT.md line 55): the C-recursion is (P-sum)/(2C4), not
  -(P+sum)/(2C4) — valuation proof unaffected (sign-insensitive) but
  flag for the writeup and for the GGV authors. Writeup gaps 2 and 3
  are now CLOSED; the only remaining gap is the frontier itself.
  NEXT SESSION: the max-plus infinity layer (design in
  CASCADE_ENGINE_PLAN.md Phase D), stages: (1) pure degree layer +
  regression ladder, (2) tie equations with explicit leading constants
  + the 21 constraint hypersurfaces, (3) automated convolution descent.
- MAX-PLUS INFINITY LAYER, STAGE 1 (2026-07-22): `cascade_engine.py` now
  couples infinity as a sixth place (`--with-inf`): the level identity
  read in degrees (v_inf = -deg), unique max forces, ties drop only via
  recorded leading-cancellation obligations — the exact max-plus dual of
  the audited min-plus descent (`descend_options_inf`, `deg_h_options`).
  Degrees (deg d2, d1, sigma, e; deg g_l) are first-class unknowns
  sandwiched between the finite-place valuation sums (join runs the
  budget DFS with caps tightened to the chosen degrees) and the window
  caps; deg e >= a + sum(b). Unlike the finite places, the infinity
  chain runs terminal..1 PLUS the level-0 anchor t^v g_1 = h_0
  (t5_multiplace_verify check 5), with internal g-zero branches for
  levels 3..1 folded into the chain enumeration — this is what makes
  the master-sum degree kills (which need the f=0 term) reproducible
  per-level. REGRESSION LADDER (`test_cascade_inf.py`, in
  `run_tests.sh`) all PASS: R0 semantics units; R1 the a=9 T2 kill
  (deg e=9 dead for all states, linear E open with obligations,
  matching T5_90_T2.md's split); R2 exactly 43/50 constant-cell kills
  with survivors d=2 x {sigma=0, z<=5} (T5_90_T1.md); R3 all seven
  T2-column margins dead with the pattern-A/B states open
  (T5_T2_COLUMN.md); R4 joint q+t+inf smoke on a9 b1000 T2 — no new
  survivors vs cascade_cones_qt.json and the g5=0 case narrowed to
  (deg e, deg sigma)=(10,8) automatically, reproducing the by-hand G5
  narrowing. Audited min-plus paths byte-identical (test_cascade_engine
  re-run). Notably the per-level chain is SHARPER than the master-sum
  arguments (a6 b1100 dies at level 5 via the g-cap, not just via the
  T6 tie). Default OFF; no swept artifacts changed yet.
- MAX-PLUS INFINITY LAYER, STAGE 2 (2026-07-22): infinity tie
  obligations now carry the EXPLICIT leading-coefficient equations
  (`cascade_inf_ties_verify.py`, in `run_tests.sh`): (A) lc(u) = lc(Phi)
  = -1024/3315 verified from the source q (lc(q)=2048, c=-1/6630);
  (B) the depth-1 initial form of an infinity degree tie is the SAME
  polynomial in leading coefficients as the residue-lemma initial form —
  verified for the backbone hypersurfaces P6/P10/P11, whose full ties
  are all realized at (deg d2, deg d1, deg sigma, deg e) = (2,3,4,5);
  the 21 CONSTRAINT lemmas of RESIDUE_LEMMAS.md therefore pin infinity
  ties to the same hypersurfaces; (C) the two arithmetic KILLS C08/C20
  act as forbidden DROPS at infinity under --residue-kills (the
  unknowns at infinity are leading coefficients in the base field, and
  C08/C20 have no Q-points); (D) `leading_cancellation` and
  `exact_identity` obligations at infinity record both sides as parsed
  labels with the exact (-1024/3315)^l constants (a9-T2 witness
  reconstructs lc(ehat)^3*lc(g5) + (-1024/3315)^5*2048*lc(e)^2 = 0).
- SUB2 INFINITY SWEEP + T2 SQUEEZE (2026-07-22): full q+t+inf sweep
  over the 420 open sub2 branches (`cascade_cones_qt_inf_rl.json`,
  guarded, clean). Pure tropical infinity adds NO new kills (every
  survivor admits some consistent degree state — checked honestly:
  the T5_T2_COLUMN kills need the F^2|G squeeze, not degrees alone).
  With `--t2-squeeze` (C24's level-5 squeeze as a tightened g6 budget
  cap, deg g6 >= sum_p v_p(g6) + 2 deg F; hypothesis m_i >= 1 follows
  from the terminal law 6+2s_i-3b_i whenever no b_i = 2 — auto-skipped
  otherwise), the engine reproduces the four column cell kills
  endogenously: **26 surviving branches (18 T1 + 8 T2), 220 flag
  cases** — the first machine artifact recording the 26-cell frontier
  (closes inventory gap I7); no unsound new survivors vs the audited
  rl artifact. Every survivor witness now carries an explicit degree
  state and exact leading-coefficient tie obligations. Regression R5
  (T2-column verdict replay, 12/12) added to `test_cascade_inf.py`.
  Spec for the independent audit: `CASCADE_INF_REPORT.md` (code-free
  semantics, the analogue of CASCADE_ENGINE_REPORT.md). Status: the
  4 branch kills re-derive tier-2 C24; the artifact itself is
  engine-proven PENDING the spec-only audit extension.
- PHASE D STATE WORKLIST, SUB2 (2026-07-22): `phase_d_states.py`
  enumerates the COMPLETE residual degree-state list per surviving flag
  case (`phase_d_states_sub2.json`): 220 cases, 7888 states, each with
  its minimal-obligation chain. CROSS-VALIDATION: the engine's
  (deg e, deg sigma) sets for all 8 open T2 cells equal the
  T5_T2_COLUMN.md R-tables EXACTLY (R9, R80-R82, R71-R74) — the hand
  analysis and the machine agree on the precise residual frontier.
  This is the direct input for the Stage 3 convolution descent.
- STAGE 3, CONVOLUTION DESCENT MECHANIZED (2026-07-22):
  `convolution_descent.py` (Codex-authored to spec, verified and landed
  here; in `run_tests.sh`) generalizes t5_90t1_constant_verify.py into
  a driver: sparse exact top-coefficient convolution of the master sum
  (ALL f-terms at every target — newly-activated lower-f terms are
  never truncated), forced-pattern recognition
  (unit*(linear)^2 / unit*unknown^2 -> substitution), honest verdicts
  (FORCED chain / UNRESOLVED residual / CONTRADICTION). GATE PASSED:
  the driver DISCOVERS the full T5_90_T1 section-3 chain 238..226
  (all nine forced values) and the final nonzero degree-226 constant
  autonomously from the source tables. First new-target probe
  (a9 b1000 T2, pattern-B states z=0 and z=6): the degree-250 leading
  equations are now EXPLICIT polynomials in (lc(d2), lc(e), lc(sigma))
  — verdict honestly UNRESOLVED at depth 1, exactly as T5_T2_COLUMN
  predicted; killing them needs elimination modulo the hypersurface
  (next convolution layer), not deeper single-coefficient extraction.
- SUB1 INFINITY SWEEP — THE 26-FAMILY BREAKS AT LOW a (2026-07-22):
  full q+t+inf sweep with residue kills and the T2 squeeze over the
  2178 sub1 standard-regime branches
  (`cascade_cones_sub1_qt_inf_rl.json`, guarded, clean):
  **108 NEW branch kills** (62 T1 + 46 T2) — a=0 and a=1 are CLOSED
  ENTIRELY (26 each), a=2 keeps only 2 branches, kills taper 24/15/12/
  4/1 at a=2..6. The infinity place bites hardest exactly where
  v=30-3a is largest, breaking the a-independence of the 26-family.
  Frontier: **279 -> 171 branches** (118 T1 + 53 T2; per-a
  2,11,14,22,25,26,26,24,21 for a=2..10), flag cases
  **2170 -> 1145**. Zero unsound survivors vs the audited rl
  baseline. STATUS: engine-proven, PENDING the spec-only independent
  audit (Codex audit_inf_cases.py in progress; must be re-run against
  this FINAL artifact since the job may have started against a
  partial checkpoint).
- ALTERNATE-REGIME DEGREE SWEEP (2026-07-22): `alt_inf_sweep.py`
  (Opus-agent-authored to the ALT_REGIME_INF.md spec, verified and
  landed here) runs the flipped max-plus chain over the 27 open
  branches (`alt_inf_sweep.json`; spot-checker
  `alt_inf_sweep_verify.py` in `run_tests.sh`). HONEST RESULT: no
  whole-branch kills on the degree layer alone — the doc's
  "generic instance dies" is the no-cancellation sub-scenario; under
  the drop-with-obligation policy every branch retains some state —
  but **33670/38360 degree states (87.8%) are killed** by three
  rigorous mechanisms (bottom-close unique-max 23870; T1 top-anchor
  2 deg d1 < w 9500; dead intermediate level 280), leaving 4690
  states with explicit drop/cancellation obligations. Agent
  [judgment] notes recorded in ALT_INF_SWEEP.md (T2 top anchor
  r6=0 seeded from ALT_REGIME.md's prescription; no r_f caps,
  conservative). The alternate front's endgame is now the same
  shape as the standard fronts: finite obligations at explicit
  degree states.
- INFINITY AUDIT PASSED (2026-07-22): `audit_inf_cases.py`
  (Codex-authored from CASCADE_INF_REPORT.md ONLY — no engine or
  test-file access, own f31_graded.txt parser with homogeneity
  self-checks, own finite-descent/join/max-plus enumeration, own sub1
  cap derivation) agrees on ALL records of BOTH final artifacts:
  sub2 420/420 (26 survive, 220 cases, all 8 infinity-removed cases
  incl. the four T2 column branch kills confirmed), sub1 2178/2178
  (171 survive, 1145 cases, 1025 removed cases confirmed), exit 0,
  ~4.5 min (re-run locally before commit). Wired into `run_tests.sh`
  (--quiet). The 108 sub1 infinity kills and the sub2 26-cell
  artifact are therefore ENGINE-PROVEN AND INDEPENDENTLY AUDITED
  (C43 promoted to tier 1). Two documented spec interpretations by
  the auditor (sub1 cap recurrences derived from the graded
  identities; post-drop leading-form convention) and one metadata gap
  fixed here: sweeps now record "t2_squeeze" in the payload (current
  artifacts rely on the documented places/residue_kills convention).
- ELIMINATION DESCENT (2026-07-22): `convolution_elim.py`
  (Codex-authored to spec, verified and landed; gates wired into
  `run_tests.sh` via --gates-only, ~2 min) upgrades the convolution
  driver to Groebner elimination modulo the accumulated coefficient
  ideal, with every leading-coefficient nonzero condition enforced by
  its own Rabinowitsch variable and CONTRADICTION only at basis {1}.
  GATES: (i) plain constant-E chain reproduced; (ii) NEW POWER — the
  T5_90_T1 tied cells deg sigma = 5 and 4 are killed MECHANICALLY
  from the master coefficients alone (degrees 242-241 resp. 242-239
  + saturation), with none of the hand identities (6)/(7)/(8)
  supplied. DESIGN NOTE (Codex, correct): with c fixed the weighted
  scaling is NOT a symmetry, so no scale gauge is assumed; leading
  coefficients stay variables with nonzero constraints. R9
  exploration (a9 b1000 T2): all seven z-states remain REDUCED —
  proper 4-34-polynomial bases within 30-60s Groebner budgets, no
  false kills, first unconsumed degree recorded per state. The
  elimination tool is validated; the R9 states need either longer
  budgets or the finite-place q-support conditions to close.
- ALTERNATE REGIME FULLY TIER 1 — G3 CLOSED (2026-07-22):
  `audit_alt_regime.py` (Codex-authored spec-only from the derivation
  docs, no checker/engine access, own f31 parser) confirms the ENTIRE
  alternate-regime chain: 26 strata, all 19+6=25 branch kills
  independently re-derived under relaxed semantics, all 27 OPEN
  verdicts with COMPLETE state partitions matched (not just totals),
  633 open states re-witnessed, 1000 sampled killed states confirmed
  dead, headline 38360=4690+33670. ~30s, re-run locally, wired into
  `run_tests.sh` (--quiet). C33/C34/C44 and the alt sweep are now
  INDEPENDENTLY AUDITED. Writeup note from the auditor:
  ALT_REGIME_L2.md's "exact residual shapes" section lacks the
  14-branch T2 subsection (recoverable from its verdict table).
- A-QUANTIFIED DEPTH LEMMAS — NO-JET-KILL THEOREM (2026-07-22):
  `RESIDUE_LEMMAS_DEPTH.md` + `residue_lemmas_depth_verify.py` (in
  `run_tests.sh`, V1-V7): census of the 1145 surviving sub1 flag
  cases confirms the affine law depth = 30-3a as the dominant t-place
  cancellation depth at every a; the top-3 a-growing tie patterns
  (C09/C02/C22) have exact depth-2/3 jet systems derived from source
  and are CONSTRAINT for ALL a in [2,10] — and more strongly, ALL TEN
  occurring t-place supports have smooth rational points, so the jet
  tower C1=C2=...=0 is solvable to EVERY order: t-place obligations
  cannot kill via local jets (t is a finite rational place; contrast
  the q-place, whose splitting-field confinement yields the only two
  kills C08/C20 — both a-independent and already consumed upstream).
  CONSEQUENCE for Phase D: the remaining closures must come from
  global/infinity arguments (elimination descent, q-support), not
  deeper t-jets — the worklists' term_cancellation obligations are
  locally unobstructed. [judgment] tags in the doc.
- GGV3 CITATION DEBT CLOSED (G7) + WRITEUP SKELETON (G5)
  (2026-07-22): GGV3 identified as arXiv:1406.0886 (confirmed via
  NEXT_CASES, GGHV22's own bibliography, and arXiv), source fetched
  to `paper_src/1406.0886_GGV3.tex`; `GGV3_CITATION_CHECK.md` verdict
  SUBSTITUTION VERIFIED — GGV3 §1's strip (inside Thm 1.8, via van
  den Essen Lemma 10.2.11, (1,1)-grading) is equivalent to the GGV1
  Props 1.13/2.1 reconstruction in T6_PREMISES.md; the citation
  delegates only the iteration of a step GGHV22 itself sets up with
  those props. T6_PREMISES.md note updated (checker re-run, 14/14).
  Writeup: `WRITEUP_OUTLINE.md` — referee-facing skeleton leading
  with the ideal-membership theorem (G5), full C0-C46 -> section ->
  tier map, blocking list B1-B4 (B1 = the frontier itself; B2 = the
  I8 unwired-checker debt; B3 = Prop 4.3 transcription still
  document-only; B4 closed by this GGV3 verification). One draft
  error corrected on landing: the outline's claim that run_tests.sh
  is missing was a path confusion (the suite lives at the REPO ROOT
  and is tracked).
- ALT COMBINATION LAYER (2026-07-22): `alt_combined.py` intersects
  the audited degree sweep with the audited finite-place cones (the
  step ALT_INF_SWEEP.md J5 flagged): joint (x,z)-coupled reachability
  over the per-place cone selections vs the state's degrees. RESULT
  (`alt_combined.json`, `ALT_COMBINED.md`; `alt_combined_verify.py`
  in `run_tests.sh`): **4690 -> 3102 residual states** (1588 killed,
  33.9% further; dominant constraint: the T1 t-place cone + h7
  anchor parity, 81% of kills). Whole-branch kills: 0 — and this is
  now a PROVABLE STRUCTURAL FACT, not an engine gap: the largest
  forced minima (deg d1 >= 9; sum v(sigma) = 12) sit exactly at the
  sub1 caps, so a top-degree state always survives valuation bounds
  alone. CONSEQUENCE: closing alternate branches requires the residue
  congruences (leading-cancellation coefficient conditions) — same
  endgame class as the standard fronts. [judgment] J1-J6 in the doc.
- Q-SUPPORT ELIMINATION, FIRST PASS (2026-07-22):
  `convolution_elim_qsupport.py` (Codex-authored, re-run locally,
  results byte-matching; NOT in run_tests.sh — 17-min exploration
  tool, results recorded here) adjoins the q-root r with q(r)=0 in
  the ideal and the exact cell structure e = gamma(y+1)^9(y-r),
  (y-r)^2 | sigma for the a9 b1000 T2 priority cell. VALIDATED
  REFINEMENT: the z=0 degree-250 basis collapses 34 -> 5 polynomials
  vs the plain run, valuations v_r(e)=1, v_r(sigma)=2 checked by
  construction, and the root-support relation has nonzero remainder
  against the plain basis (genuinely new information). VERDICTS: all
  seven R9 z-states remain REDUCED — every state consumes 251
  (identity) and 250 (added), then times out testing 249 under the
  120s Groebner budget; no false kills, no contradictions. NEXT
  ATTACK on this cell: longer budgets/better ordering (e.g.
  eliminate r first via resultants) on the 5-polynomial z=0 ideal.
- BATCH CONVOLUTION, SUB2 FIRST HARVEST (2026-07-22, session close):
  `batch_convolution_sub2.py` + `batch_convolution_sub2.json` +
  `BATCH_CONVOLUTION_SUB2.md` (raw pass artifacts committed for
  provenance). Dedup 7888 raw -> 1782 unique degree states; 194
  attempted under triage+90s budgets: **60 candidate state kills
  PENDING AUDIT** (40 CONTRADICTION + 20 forced-lc degree drops),
  130 honestly UNRESOLVED, 2 FORCED, 2 skipped. THE CONTENT: 38
  kills re-derive the known a9 T1 constant-E block (C22 validation);
  **22 kills at (10,(0,0,0,0),T1) constant-E are NEW** — the
  zero-obligation "pure Phase D" cell. SOUNDNESS GROUNDING (reviewed
  on landing): lc(e-part) = gamma != 0 is the state's
  degree-exactness, NOT a scaling gauge (the identity is weight-125
  homogeneous only with c -> l^17 c, so at fixed c scaling is no
  symmetry — independently confirmed by the codex elim lane);
  gamma-only residual systems kill via gcd (no common nonzero root).
  Pass 1 (free e-scale) produced zero kills — recorded as the
  negative control. COVERAGE BOUNDARY (honest): unattempted = 391
  tier-2 a10 + 252 tier-3 + 945 tier-4 unique states; UNRESOLVED !=
  survival (q-support dropped). Noted survivor structure: all 24 a8
  constant-E states stall on 8192 b_lc^2 + 9945 gamma^3 s_lc^2 — a
  sign-split argument would likely kill half. NEXT: the spec-only
  convolution-kill auditor (promote the 60 to tier 1), then the sub1
  batch and the a10 tail. Housekeeping: 7 orphaned python workers
  killed on wrap-up.
- CONVOLUTION KILLS AUDITED 60/60 (2026-07-22, overnight round):
  `audit_convolution_kills.py` (Codex-authored spec-only from
  BATCH_CONVOLUTION_SUB2.md + the master identity, own parser with
  homogeneity self-checks, no driver/engine access, read-site
  inspection confirms it reads only f31_graded.txt and the kills
  JSON) independently re-derives ALL 60 candidate kills: CONFIRMED
  60/60, UNDECIDED 0, DISAGREEMENT 0, exit 0, ~8 min (re-run
  locally). Wired into `run_tests.sh` (--quiet). The 38 a9 + 22 a10
  constant-E state kills are now ENGINE-PROVEN AND INDEPENDENTLY
  AUDITED — the sixth spec-only auditor; the convolution harvest
  pipeline is now a tier-1 production line.
- R9 ORDERING BREAKTHROUGH, STILL OPEN (2026-07-22):
  `convolution_elim_r9.py` (Codex, landed unwired — exploration
  tool): the z=0 stall was the ORDERING, not the ideal — lex order
  (g0, gamma, a0..a4, r, rabinowitsch) completes degree 249 in 1.8s
  (six-polynomial proper basis) where grevlex timed out at 120s;
  r-elimination gives a four-polynomial r-free basis at 250.
  z=0 remains honestly unresolved (no {1}, no solution point).
  NEXT: continue the winning lex ordering to deeper coefficients
  with long budgets, and/or pipe the r-free system to Singular.
- SUB1 BATCH CONVOLUTION — FIRST T2 KILLS (2026-07-22, overnight):
  `batch_convolution_sub1.py` + artifact + md (Opus lane, verified,
  landed): 44117 raw sub1 states dedupe to 4994 unique tuples.
  TRANSFER PASS: all 194 sub2-attempted tuples recur verbatim in
  sub1 (window-independence of the master identity) — 60 kills
  transfer and inherit the CONFIRMED 60/60 audit (identical
  computations). FRESH RUN (149 states in budget): **68 new
  CONTRADICTION kills**, all incompatible-gamma^17 pairs — including
  **28 T2 constant-E kills at a=10, the first T2 state kills in
  either window** (exact block boundaries recorded: dies for
  deg_sigma <= 6, resumes UNRESOLVED at deg_sigma >= 7), plus 40 T1
  constant-E a=10 d2_zero kills. Combined: 128 candidate kills (60
  audited-by-transfer + 68 PENDING AUDIT — extend
  audit_convolution_kills.py to the sub1 file next round). Honest
  coverage: 343/4994 unique attempted; the a=10 constant-E block was
  cut MID-KILLING-STREAK by the wall budget — the overnight marshal
  job continues exactly there.
- SUB2 HARVEST ROUND 2 + A8 VERDICT (2026-07-22, overnight):
  `batch_convolution_sub2_round2.json`/md (Opus lane, verified,
  landed): **82 more candidate kills** on the a10 tail (262 more
  states attempted; kill structure fully regular — all deg_d1 <= 4
  cells across deg_d2 in {0,1,2}; survivors concentrate at deg_d1 in
  {5,6}, deg_sigma in {7,8}, deg_d2 in {3,4}). Sub2 program totals:
  **142 unique candidate kills, 456/1782 attempted (25.6%)**, all
  PENDING AUDIT (extend the sixth auditor next round). A8 SIGN-SPLIT
  RESOLVED HONESTLY — NO KILL: the splitting-field confinement
  justifying C08/C20 is a q-place phenomenon; the batch's leading
  coefficients live over arbitrary char-0 K, and the stall residual
  8192 b^2 + 9945 g^3 s^2 has an explicit sympy-verified complex
  witness (1, 1, i*sqrt(9945/8192)). Over an ordered field the
  gamma>0 branch dies (recorded as CONDITIONAL); the 24 a8 states
  remain UNRESOLVED, needing deeper coefficients or q-support.
- PHASE F ADOPTED (2026-07-22, night): external review (GPT) assessed
  and adopted as `PHASE_F_PLAN.md` — the global residue obstruction
  algebra. Key inversion of our own negative results: no-jet-kill +
  caps-attained together mean survivors are LOCALLY flexible but
  GLOBALLY divisor-rigid; defect-0 polynomials are lambda*prod(y-s)^v
  with all local leading residues DETERMINED algebraic numbers
  (root-difference products in the S4 splitting field — the C08/C20
  kill shape, manufactured at scale). Pipeline: defect compression ->
  Hermite/CRT jet gluing -> Galois-aware coordinates -> toric
  quotients (fixed-c fiber kept as an equation; resolves the gauge
  question) -> saturated initial-ideal tests (tropical VARIETY not
  prevariety). Work items F1-F6 in the plan; F1 (defect histograms
  over all surviving states, the go/no-go compute) dispatched.
  Metric shift: canonical global residue ideals closed, not raw
  state-kill counts.
- ALT DEPTH-1 CONGRUENCES: THE ENDGAME IS 19 HYPERSURFACES
  (2026-07-22, night): `alt_residue_congruences.py` + artifact + md +
  verifier (Opus lane, verified, landed; verifier in run_tests.sh)
  writes the exact depth-1 leading-coefficient systems for all 3102
  alternate states. RESULT: 0 kills, 3102 CONSTRAINT — and a sharp
  structural collapse: every survivor's obligations sit at LEVEL 0
  (the bottom close E^21 h_0 + u r_0 = 0); levels 1-6 are always
  dominated (h_6..h_1 forced at maxima, no residue equations). The
  3102 h_0 initial forms collapse to **19 distinct hypersurfaces**
  in (D,X,S,E), each with an exhibited rational point; required
  cancellation depths run to 17 (2062 states). SOUNDNESS CRUX
  (K2): 523 states carry C08/C20-shaped ties but at the MAX,
  non-obligatory — refusing drops changes nothing
  (kills_on == kills_off verified); claiming those would be unsound.
  CONSEQUENCE: the alternate front = 19 canonical level-0 ideals +
  deep convolution — precisely the Phase F target shape (F2 pilots
  should attack the tightest of the 19 with divisor reconstruction).
- F1 DEFECT HISTOGRAMS — PHASE F IS VIABLE (2026-07-22, night):
  `phase_f_defects.py/.json`, `PHASE_F_DEFECTS.md` (Opus lane,
  verified, landed): exact tightest-join defects for all 55,107
  surviving states (sub2+sub1 EXACT via complete join enumeration —
  Pareto sets are tiny; alt = sound upper bounds from forced sums);
  383,351 delta values, ZERO negatives (independent confirmation of
  the audited coupling). HEADLINE (reconstruction-relevant divisors
  d1,d2,sigma,e,terminal-g; the auxiliary g4-g6 chain excluded as a
  non-target): sub2 18.1% / alt 14.1% / sub1 8.1% of states have ALL
  defects <= 2; 243/540/438 states are FULLY FORCED (every divisor
  determined up to one scalar). Tightness rises with a_t and is 2.6x
  higher on T2 than T1. Pilot rankings recorded: fully-forced
  sub2/sub1 a9_b1000_T1 (~38 states each) and the alt a11_b3000_T1 /
  a11_b1111_T1 both-min families — exactly the F2 targets. VERDICT:
  the global-reconstruction route is live; F2 decides its power.
- F2 PILOT — FIRST RECONSTRUCTION KILLS, INCLUDING A WHOLE ALT STATE
  (2026-07-22, night): `phase_f2_pilot.py` + `PHASE_F2_PILOT.md` +
  verifier (Opus lane; verifier re-run locally, ALL PASS). THE
  MECHANISM WORKS: for defect-0 polynomials every sub-leading
  coefficient is a determined multiple of one scalar, and feeding
  them into the FULL level-0 tie tower over-determines the leading
  data — kills live at depth 2-3, invisible to the depth-1 layer
  (unobstructed precisely because it treated (D,X,S,E) as free).
  VERDICTS: Pilot A (a11_b3000_T1 both-tight stratum) KILL — unit
  ideal at depth 3, UNIVERSAL in the marked root; Pilot B
  (a14_b0000_T2 (6,-,12,14)) **GENUINE WHOLE-STATE KILL** — the
  reconstructions are FORCED (v_t(sigma)=12=deg sigma, deg e=a=14)
  and even a fully FREE degree-6 d2 cofactor cannot rescue: the
  FIRST alternate-regime state closed, by exactly the residue
  congruences J1 put out of valuation reach. CONTROL (e generic):
  UNOBSTRUCTED with explicit triangular solution — the engine
  discriminates; the kills are mathematical. Pilot C: 4 residual
  conditions in 2 unknowns, no torus point found, exact GB
  certification PENDING (needs the sigma=4d0-d2^2 coherence).
  CONTROLLING INVARIANT: the defect profile. All kills PENDING
  AUDIT. NEXT: scale F2 over F1's 1221 fully-forced states and the
  tight alt families; Galois-stable a11_b1111 pilot; wire the
  verifier into run_tests.sh.
- F2 SCALED — 20 RECONSTRUCTION KILLS, TWO BRANCHES 6/8 CLOSED
  (2026-07-23, night): `phase_f2_scale.py/.json`, `PHASE_F2_SCALE.md`,
  verifier (re-run locally, ALL PASS; wire next hygiene pass). The
  reconstructible alt frontier = 40 states (d1, sigma AND e all
  forced defect-0, each with a UNIQUE cone split — reconstruction is
  well-defined, not a choice; d2 left fully free, conservative).
  CENSUS: **20 KILLED** (saturated unit ideals over Q / Q[r]/(q) /
  Q[r1,r2]), 18 narrowed/unobstructed, 2 pending on GB cost. Killing
  families: a11_b1111_T1 **7/8 over Q** (flagship, depth 2, initial
  form -729E(9E^3+8X^5)); a11_b3100_T2 6/8 (TWO marked roots);
  a12_b1110_T2 6/8; a14_b0000_T2 (= Pilot B, reproduced by the
  general engine). STRUCTURAL BOUNDARY (J2, honest): the kill has a
  d2-FREEDOM THRESHOLD — deg d2 <= 4 dies at depth 2, >= 5 is
  rescued (the control mechanism); the 6/8 branches' survivors are
  exactly the deg d2 in {5,6} states. Pilot C: residual structure
  pinned exactly (4 residuals in (D,E) over Q(r,c9)); certification
  still stalls. Sub2: mechanism demonstrated on a10_b0000_T1
  (narrowed); the 203 b!=0000 forced states need geometric-regime
  reconstruction — deferred. All kills PENDING AUDIT. The alternate
  frontier's reconstructible core is now half dead in one pass.
- OVERNIGHT COMPUTE RESULTS (2026-07-23 morning): JOB A (long batch,
  full 5h wall, clean exit): 1162 more unique states attempted across
  the remaining pools -> **30 more kills** (15 CONTRADICTION + 15
  degree-drop; census in `batch_convolution_overnight.json`), 1131
  honestly UNRESOLVED (the tier-3/4 pools are mostly not
  constant-E-shaped — the forcing mechanism thins out there, as
  expected). JOB B (R9 lex continuation): z=0 extended 249->246 ALL
  PROPER (basis 6->9), degree 245 exceeded the 1500s per-degree
  budget — the ideal is genuinely deep even under the winning
  ordering. F3 (Singular bridge): codex delivered
  `convolution_elim_r9_singular.py` with generated+audited Singular
  programs but its sandbox was DENIED WSL access (the overnight
  internet outage also left the audit-r2 codex job queued 8.5h —
  re-dispatched this morning); the Singular attack now runs locally.
  Convolution-kill total pending audit: 68 sub1 + 82 sub2-r2 + 30
  overnight = 180.
- F2 -> SUB2 MARKED ROOTS (2026-07-23): `phase_f2_sub2.py/.json`,
  `PHASE_F2_SUB2.md`, verifier (re-run locally, ALL PASS): the
  deferred sub2 extension — reconstruction against the FULL master
  identity with per-place valuations from the exact engine profiles,
  d2 imposed where forced. CENSUS over 201 targeted states: **23
  KILLED (16 NEW)** — 15 over Q[r]/(q), exactly the deferred b1000
  deg_e=10 marked-root states absent from every batch list, + 1 over
  Q. Flagship hand-verifiable chain: the tower forces r = 1/4 in Q,
  contradicting q's irreducibility (the arithmetic-kill pattern,
  manufactured by reconstruction as Phase F predicted). KEY
  STRUCTURAL FINDING [J3]: the full master identity is STRICTLY
  stronger than the level-0 tie (a10_b0000_T1 NARROWED there, KILLED
  at depth 2 here) — the alt lane should be re-run against the full
  flipped identity. Honest coverage: 152/201 skipped (no unique
  split) + 26 pending on GB cost; NO cell or branch closes. Kills
  PENDING AUDIT.
- PUBLICATION PACKAGE (2026-07-23): publication-prep lane landed:
  `PUBLICATION_AUDIT.md` (repo root; key findings: the four tracked
  arXiv .tex sources are third-party copyrighted — REMOVE AT RELEASE
  TIME via a public branch with a links-only paper_src README, NOT
  deleted now since t6_premises_verify.py reads them; two personal
  path leaks to redact; 33 logs + 9 pickles to drop;
  math_stuff_field_audit/ is benign and already git-ignored; zero
  credentials), refreshed top-level README (current truth, honest
  "what is NOT claimed" paragraph), `VERIFICATION_QUICKSTART.md`
  (15-minute skeptic ladder with real captured outputs), and
  `CONTACT_DRAFT.md` (short + long GGHV-author notes, DRAFT — owner
  review only). Release checklist = execute the audit's
  REMOVE/REDACT items on a public branch at flip time.
- LEAN CERTIFICATE LANDED (2026-07-23): the headline theorem is now
  KERNEL-ACCEPTED LEAN 4 (`lean_certificates/`, toolchain v4.32.1,
  lake build ~3m20s cold): `Cert.f37_certificate` proves the
  integer-cleared identity Dmul*f31 = c1 G1 + c2 G2 + c3 G3 + c4 G4
  (D = 46875) in Z[d~2..Phi] via a SELF-CONTAINED verified
  sparse-poly library (packed-exponent sorted lists, structurally
  recursive ops) — no mathlib. Data generated from the SAME sources
  as f37_sat_verify.py (import, not copy) and independently
  re-verified in sympy before emission. TRUST STORY (honest, in
  LEAN_CERTIFICATE.md): full identity via native_decide (trusted
  base = Lean compiler; ~27s); pure-kernel decide cross-checks on
  real generator instances are AXIOM-FREE; pure decide on the full
  identity is scale-blocked (~6.7GB whnf), removable via Kronecker
  bignum encoding (FEASIBILITY.md). The program's headline claim now
  carries the week's trust currency; remaining Lean debt for the
  full program sketched in LEAN_CERTIFICATE.md.
- FULL KILL LEDGER AT TIER 1 (2026-07-23): `audit_convolution_kills_r2.py`
  (240/240: 68 sub1 fresh + 60 transfers tuple-verified + 82 sub2-r2
  + 30 overnight) and `audit_reconstruction_kills.py` (22/22 F2
  reconstruction kills, spec-only from the F2 documents with its own
  parser and saturation) both CONFIRMED by codex AND re-run locally,
  zero disagreements, zero undecided; both wired into `run_tests.sh`.
  THE ENTIRE KILL LEDGER — 262 state kills across three mechanisms
  (forcing convolution, gamma-gcd, divisor reconstruction) — is now
  ENGINE-PROVEN AND INDEPENDENTLY AUDITED. The audit family stands at
  EIGHT spec-only auditors, all green in the suite.
- MODULAR TRIAGE — THE DEEP END IS MOSTLY EMPTY (2026-07-23):
  `modular_triage.py/.json`, `MODULAR_TRIAGE.md` (landed): 59
  resistant subsystems x 3 good primes (10007/10009/100019; bad
  primes {2,3,5,13,17} identified from disc(q) and avoided), ZERO
  mixed verdicts. HEADLINE: 24 LIKELY-EMPTY / 31 LIKELY-SOLVABLE /
  4 timeout. R9 z=0..3: UNIT at every prime — and MORE master
  coefficients collapse the GB FASTER (8 coeffs resolve in seconds
  where 6 time out): over-determination is our friend. All 10
  sampled sub2 T2 pattern-B (a7/a8) systems: UNIT — the whole T2 tie
  frontier looks killable. 10 of 18 alt NARROWED states flip to UNIT
  at full tower depth (the rational NARROWED labels were a depth-2
  artifact). CAVEATS (honest): mod-p UNIT is a PREDICTION, not a
  proof; and the 24 a8 LIKELY-SOLVABLE verdicts used only the top
  2-3 coefficients — System 1's own lesson says run the FULL
  coefficient sets before concluding a8 needs a different argument
  (real-sign routes are unavailable over C anyway per the a8 audit).
  NEXT (dispatched): exact full-depth kills for the 10 flipped alt
  states + the a7/a8 T2 cells + full-coefficient a8 probe.
- D2-THRESHOLD RESOLVED — TWO MORE KILLS, [J2] REFUTED (2026-07-23):
  `d2_threshold.py`, `D2_THRESHOLD.md`, verifier (re-run locally,
  ALL PASS). SOUNDNESS FINDING FIRST: d2 has NO forced valuations in
  the alt branches (ALT_REGIME_L2 sec 2: all k=0..6 allowed) — the
  impose-d2 route fails honestly and nothing was imposed. THE REAL
  MECHANISM: the four scale-lane survivors were DEPTH-CAP artifacts.
  With d2 FULLY FREE (conservative), deeper towers kill:
  a12_b1110_T2 deg-d2=5 and a11_b3100_T2 deg-d2=5 both die at depth
  8 (exact saturated unit ideals over Q[r]/(q) and Q[r1,r2]
  respectively) — refuting PHASE_F2_SCALE [J2]'s d2-freedom
  threshold for deg 5. Both branches now **7/8 killed**; the sole
  survivor in each is the deg-d2=6 state, PENDING on compute only
  (irreducible j0 + marked-root field -> 12GB GB blowups; strongly
  inferred dead by the Pilot-B depth-9 analogy, NOT claimed). Kills
  ENGINE+LEMMA-PROVEN PENDING AUDIT. Queued: mod-p precheck + long
  detached rational runs for the two deg-6 states.
- TRIAGE HARVEST — 19 EXACT KILLS, ZERO REFUSALS, A8 REVERSED
  (2026-07-23): `triage_harvest.py/.json`, `TRIAGE_HARVEST.md`,
  verifier (re-run locally, PASS; independent sympy-Groebner
  cross-CAS re-derivations). (1) a7/a8 T2 pattern-B sample: **10/10
  EXACT KILLS over Q** (<2s each, 8 coefficients + q(r)) — the T2
  tie frontier is falling exactly as the modular map predicted.
  (2) Flipped alt NARROWED states: **8/10 EXACT KILLS** at full tie
  depth (incl. the over-Q flagship a11_b1111_T2 and the two-root
  a11_b3100_T2); the 2 misses are 400s-budget cost, not refusal.
  (3) R9 z=3: cost-bound over Q (astronomical integer swell), mod-p
  EMPTY stands — needs the multi-modular route. (4) **A8 REVERSAL**:
  the LIKELY-SOLVABLE triage verdict was a 3-coefficient truncation
  artifact — at 10-16 coefficients the deg_sigma 5/6 states flip
  PROPER->UNIT on all primes, and a8_dd2-inf_dd10_dsig5 is UNIT
  OVER Q in 0.92s (exact kill): a8 needed MORE EQUATIONS, not
  real/sign mathematics; deg_sigma 7/8 remain cost-bound.
  HEADLINE DISCIPLINE FACT: zero states anywhere returned
  PROPER-over-Q where mod-p said UNIT — the triage map is
  trustworthy. All 19 kills PENDING AUDIT (queue a GB-kill auditor
  extension next round: these are saturated-unit-ideal certificates,
  auditable by cofactor extraction like the NulLA style).
- PUBLIC REPO LIVE (2026-07-23): github.com/wstrinz/plane-jacobian-72-108
  (public, commit 548237e) — the clean release tree (layout: README /
  VERIFICATION / d2 / d3 / lean / docs; PUBLICATION_AUDIT exclusions
  applied: no third-party tex, no drafts, no leaks/logs). Publish
  gate: d3 verify + f37_sat_verify re-run green IN the release tree;
  the builder lane's full-suite run continues as belt-and-braces
  (any late failure = fix-forward commit). Contact drafts can now
  carry the real URL. The private workspace repo remains the working
  home; releases flow one-way into the public tree.
- FULL-SYSTEM BRIDGE — THE CONCEPTUAL GAP IS RETIRED (2026-07-23):
  `full_system_bridge.py`, `FULL_SYSTEM_BRIDGE.md`, verifier (13/13,
  re-run locally; V2 proves bridge >= cascade via the audited
  ideal-membership cofactors; V4 re-derives the pilot on fresh prime
  32003). The bridge: state -> augmented G-system with stripped spare
  unknowns dm2/dm3/dm4 (weighted-homogeneity strip proven: G-weights
  156/168/180/204), ansatz cost only 45 (sub2) / 66 (sub1) spare
  coefficients — ~122 quadratic equations per state. PILOT VERDICT:
  the a8 deg_sigma=7 state (MODULAR_TRIAGE "LIKELY-SOLVABLE") is
  **UNIT OVER EXACT Q IN 8.75s** under the pure G-system (zero f31
  coefficients needed) — the a8 germ does not lift; no surviving
  germ found. CONSEQUENCE: everything that resists f31-alone has a
  cheap second gate; the "genuinely-solvable residue" risk is now
  carried only by [judgment] items (k=6..8 window-cap transcription,
  low risk) and whatever the bridge itself fails to kill — none
  found yet. NEXT (dispatched): bridge sweep across the R9 column,
  the structural blowup cases, and every remaining UNRESOLVED class.
- GPT-PRO REVIEW IN (2026-07-23): external review of the public repo —
  verdict: "a legitimate flag ... one mathematically substantive and
  unusually clean result" (the f37 elimination as centerpiece; the
  membership identity ALONE suffices — the Groebner exactness claim
  is not load-bearing, a better referee framing). REAL DEFECTS
  CAUGHT (fix lane running): (1) clean-checkout break —
  envelope_bounds_verify reads the excluded arXiv tex, so the public
  suite cannot pass from a fresh clone (the builder had gated
  against a tree still containing it); fix = original
  upstream_facts.json + optional provenance mode + a clean-clone CI
  gate; (2) pickle in the verification path — replace with canonical
  generators.json; (3) "over every field" overclaim — the Lean
  multiplier D = 46875 = 3*5^6 blocks char 3,5; scope to char 0;
  (4) f37_sat_verify Part B is tautological (reframe, not delete);
  (5) Lean claim wording narrowed; (6) CITATION TODOs + env pinning;
  plus D3 fiber -> explicit lex-Groebner triangular basis, and one
  CURRENT_STATUS.md superseding the inventory's log-vs-truth
  ambiguity. WORLD UPDATES from the review: Tao's 2026-07-21
  digestion post analyzes D3 and identifies the dim-2 obstruction
  (resultant scaling exponent n-2 = 0 — the plane's exceptionality
  in one line); the competing Helali program's OWN audit retracted
  its (72,108) closure claim (Zenodo now: "strong machine evidence
  ... not yet a complete proof") — NO competing completed proof
  exists; fresh preprint arXiv:2607.20210: Gm-equivariant plane
  Keller maps are automorphisms. STRATEGIC (recorded for the next
  campaign): the reviewer ranks the scalable outputs as (a) a
  PARAMETRIC pre-resultant elimination theorem (Bezout/Koszul
  structure of the certificate; same ideal at (108,144)?), (b) a
  global divisor/infinity obstruction engine (rewrite ledgers as
  divisor inequalities; Riemann-Roch/Mason arguments), (c)
  proof-producing tropical search with Farkas certificates
  replacing bespoke branch scripts. Suggested GGHV opener quoted
  verbatim in the review and adopted for the contact draft.
- BLOWUP DIAGNOSIS FINAL — FIRST WHOLE BRANCH AT 8/8 + THE RAM
  ANSWER (2026-07-23): `BLOWUP_DIAGNOSIS.md` + tool bridges landed
  (msolve 0.10.1 built from source in WSL — the decisive engine;
  nulla_unit_cert.py live with an honest negative on deep systems).
  FIVE new exact char-0 kills (msolve, PENDING AUDIT; .ms systems
  archived for independent re-derivation): **a12_b1110_T2 deg-d2=6 —
  taking that branch to 8/8, the program's FIRST WHOLE-BRANCH
  CLOSURE (pending audit of this final kill + the roll-up)**;
  a11_b1111_T1 #17 (last triage indeterminate); sub2 s14/s38/s94
  (3 of the 26 pending-hard class — the rest should fall the same
  way). DIAGNOSIS TABLE (11 cases): most "blowups" were
  ENGINE-SWELL (sympy/Singular Buchberger) — msolve F4 kills them
  in seconds-to-minutes at <=1GB; the genuinely STRUCTURAL set is
  R9 z=4-6 + sub2 s268 + the a11 two-root rational certificate
  (fail mod-p under 3+ engine/order combos — reformulation via
  resultant pre-elimination + lex + wall-clock, NOT hardware).
  OWNER BOTTOM LINE: renting RAM clears nothing; the 12GB blowup
  was killed mod-p in 242MB by the right engine. Housekeeping: two
  permission-blocked WSL Singular orphans killed at session level.
- S-UNIT LAYER (F5) LANDED — A PRECISION TOOL, HONESTLY CALIBRATED
  (2026-07-23): `s_unit_layer.py`, `S_UNIT_LAYER.md`, verifier
  (re-run locally, ALL PASS; Mason + Brownawell-Masser stated
  verbatim and cited, tier-4). GATE: the sigma-locus kill re-derives
  EXACTLY through the fresh engine — and is newly framed as
  a-UNIFORM (deg A = deg B = 34 independent of the b-split: one
  inequality erases the family across both windows). NEW COMPLETED
  KILL (a-uniform, both windows): the d1=d2=sigma=0 defect-0 corner
  collapses the master identity to e^17 = const * Phi^5, forcing
  17*v(e) = 5 — no integer solution (17 does not divide 5, nor 150
  at t): the family is EMPTY by pure divisibility. CALIBRATION
  (the important honest finding): raw n-term BM on the full 88-term
  identity is useless (binom(87,2)*4 ~ 15000 >> 250) — the layer's
  power concentrates at n <= 3 AFTER structural collapse; 17 sparse
  cells have favourable BM bounds but owe 254-1022-leaf subsum trees
  (labeled bm_candidate_pending_subsums, NOT kills; bounded witness-
  DB work). CONSEQUENCE: S-unit = the precision family-killer for
  collapsed/sparse cells; the bulk stays bridge territory. Zero
  coefficient-field hazards (integer valuation arithmetic only).
- SECOND GPT-PRO PASS ON THE HELALI CLAIM (2026-07-23): correction —
  the earlier comparison had matched a DIFFERENT public program (the
  Santibanez project); there are at least TWO external (72,108)/
  (75,125) efforts. Current read on Helali: SERIOUS BUT UNVERIFIED
  announcement (Zenodo DOI mentioned, artifacts not yet retrievable);
  all technical specificity concerns (75,125) — internally
  disciplined (mod-101 discovery vs exact pivot certification
  correctly separated; (75,125) honestly left open at 578 vars /
  2209 equations) — while the (72,108) exclusion sentence has NO
  audit surface. Key epistemics adopted: exact computation on a
  chosen system != proof the system exhausts the geometry; the
  danger lives in the normalization->supports->branches chain, and
  our field-split episode is exactly the warning to contribute.
  ACTIONS: draft C superseded by C-v2 (GPT-Pro's message + the
  seven-question artifact checklist + four-outcome table, in
  CONTACT_DRAFT.md); dated external-claims note added to the public
  repo; OUR STATUS UNCHANGED — (72,108) remains open in this
  program until an independently replayable certificate exists.
- UNIFIED KILL LEDGER + ROLL-UP — THREE BRANCHES CLOSED AS COMPUTED
  FACT (2026-07-23): `state_kill_ledger.py/.json` + `frontier_rollup.py`
  + `FRONTIER_V2.md` (verified by regeneration; cross-check asserts
  per-cell totals == source state totals exactly, sub2 7888/sub1
  44117, zero drift). LEDGER: 388 distinct killed states across 16
  ingests (350 audited, 38 pending incl. all alt/corner); ZERO
  verdict conflicts anywhere; 9 multi-source dedups. CLOSURES
  (CLOSED-PENDING-AUDIT): **a12_b1110_T2 8/8** (three lanes
  recombined: scale d2<=4 + threshold d2=5 + msolve d2=6 — the
  scattered-kill recombination is the ledger's raison d'etre),
  **a11_b1111_T1 8/8** (bonus finding), **a14_b0000_T2 1/1**;
  a11_b3100_T2 at 7/8, ONE state (msolve-timeout deg_d2=6) from
  closure. HONEST SCALE: the phase-D universes hold 52,005 states,
  0.70% killed — the closures live in the alt defect-0 layer; the
  bulk remains bridge/family-lemma territory. PROCESS FINDING: 27
  redundant re-attempts exposed (overnight lane re-ground states
  already dead elsewhere) — all future worklists must consult the
  ledger. 27 a8 kills held as PENDING-AMBIGUOUS-MAP (deg_e
  unpinned) rather than counted — resolve in the audit round.
- DIVISOR-LEMMA ENGINE LANDED (2026-07-23): `divisor_lemmas.py`,
  `DIVISOR_LEMMAS.md`, verifier (re-run locally, ALL PASS). **L1 —
  the new consolidation**: defect-0 determinacy is ONE lemma — the
  confluent-Vandermonde determinant = +/- prod m! * prod
  (s_i-s_j)^{k_i k_j} != 0 by separability of t*q, over ANY field —
  replacing every per-state determinacy fact; the most (75,125)-
  transferable object yet. L2 (flagship a11_b1111 as rank+torus
  lemma), L3 (17-divides-5 corner), L4 (sigma-locus degree overflow),
  L5 (C08/C20 square classes) = exact re-proofs in divisor language
  (the consolidation goal; no new coverage claimed). JUDGMENT [J3]:
  the kills are NOT all rank statements — L1 isolates the closing
  obstruction into few variables; the obstruction itself varies
  (torus/divisibility/height/square-class). HUNT: the 17 S-unit
  candidates explained — deg sigma <= 8 sits exactly at the blocking
  cap; top datum only NARROWS (2048 Phi^5 = 6561 e^17 hypersurface);
  each needs ONE depth-2 coefficient with L1-determined sub-leadings
  (the r=1/4 mechanism) — bounded queue, PENDING.
- THIRD EXTERNAL REVIEW FOLDED IN (2026-07-23): strategy-grade pass;
  verdict "mechanism-close, not proof-close" for 108 (its near-term
  list is largely our running board — roll-up already landed beyond
  its snapshot). Adopted into INDUCTIVE_PROGRAM.md: the induction-
  variable reframe (chain type/power pair/multiplicities/j, NOT
  degree), the Sigma signature, the F_2-family framing of (75,125),
  two new ideas (monotone defect potential on chains — cheap first
  test vs the 144 data; configuration-space Galois descent for the
  arithmetic layer), the case-compiler spec with our both-
  presentations amendment, and the 70/25/5 allocation.
- REVIEW FIXES LANDED + v0.1.0 TAGGED (2026-07-23): all 8 external-
  review items implemented and gated in BOTH trees (workspace cd20305,
  public 5c9d96f, pushed; public tagged v0.1.0): upstream_facts.json
  replaces the tex dependency (12 mandatory checks offline; provenance
  optional); generators.json de-pickles the verification path
  (verified passing WITHOUT the pickle in the public tree);
  char-0 scope corrected everywhere (D = 3*5^6 reason recorded);
  Part-B reframed honestly; Lean wording narrowed; CITATION +
  requirements-lock + ENVIRONMENT pinned; D3 fiber now lex-Groebner
  triangular ([z+1/4-27/4x^2, y+3/2x, x^3-x]); CURRENT_STATUS.md is
  the single current-truth document (inventory headed as historical).
  tools/clean_clone_check.py (audit-hook file tracing, 72 files, all
  tracked) is now the FIRST line of both suites. Zenodo deposit
  remains owner-side.
- POTENTIAL PROBE LANDED (2026-07-23): `potential_probe.py`,
  `POTENTIAL_PROBE.md` (verified by re-run). VERDICT: two genuine
  MONOTONE-CANDIDATE potentials across the (72,108)->(108,144)
  corner transition — deg Phi = 64a^2-8a-2 and mult_(y+1) Phi =
  8a^2-a, BOTH cases exactly on the closed forms, no free
  parameters. The section-slack candidate is monotone but jumps by
  PLACE MIGRATION (4 -> 276 = 4*69: the residual quartic moves
  inside C4) — the "fixed quantum per transition" reading is NOT
  supported; the potential idea survives in the refined form
  "closed-form-in-a monotone quantities exist; freedom moves
  discontinuously between places at corner steps". Level-1
  (defect-vector) potentials are NEEDS-DATA (no 144 cascade).
  DECISIVE NEXT DATA POINT (queued): the light Phi = f*C4^N
  derivation for (75,125) — tests power-a invariance vs corner
  coincidence; no cascade required.
- BRIDGE SWEEP FINAL — 35 EXACT KILLS, ZERO TWO-GATE SURVIVORS
  (2026-07-23): `bridge_sweep.py/.json`, `BRIDGE_SWEEP.md`, verifier
  (re-run locally, 4/4 incl. cross-CAS legs). LOUD NEGATIVE: no
  state anywhere returned PROPER under the bridge — the first true
  two-gate survivor DOES NOT EXIST in anything affordable. **A8
  CHAPTER CLOSED** (all 24 states exact-Q, incl. the 12 deg_sigma
  7/8 the harvest could not certify — 0.6-24s each; no real/sign
  mathematics needed anywhere). a7/a8 T2 deg_e=8 classes complete
  (10/10). **R9 z=0 KILLED exact-Q — first ever in the column**
  (3.4s); z=1..6 are COST not survival (five formulations exhausted;
  z<=3 logically empty over F_p-bar by f31-in-G-ideal + triage
  UNIT). STRUCTURAL FINDING: the bridge's cost wall tracks deg_e
  hitting the sub2 cap (engine-independent) — R9 z>=1 / a10 exact
  certificates need symbolic spare-unknown elimination, not budget.
  Alt-regime honestly SKIPPED (bridge is a standard-window object;
  applying it would be unsound). Kills PENDING AUDIT.
- GALOIS DESCENT PILOT LANDED — LAYER 2 HAS ITS TEMPLATE
  (2026-07-23): `galois_descent_pilot.py`, md, verifier (re-run,
  PASS). THE LEMMA: C08 kills iff sqrt(105) is not in the splitting
  field L; for S4/C4 quartics this is the single obstruction
  polynomial Obs = 105*disc(q) not a square (product of residue-side
  and forcing-side discriminants); A4 always kills; D4/V4 need the
  resolvent subfields — PROVEN necessary by an explicit V4 witness
  where the naive Obs errs. Verified three ways: specialization to
  our disc-17 quartic (105*17 = 1785 non-square -> kill, matching
  the audit); a SHARP counter-instance (disc class 105 -> obstruction
  vanishes and an explicit all-nonzero splitting-field solution
  exists); the V4 branch witness. C20 fits verbatim (170*17: class
  10 -> kills). The reusable machine: residue relation -> obstruction
  subgroup of Q*/(Q*)^2 -> membership in the forcing divisor's
  square-class subgroup. Non-quadratic rows (C12-type) need higher
  resolvents — scoped, not claimed.


## 2026-07-23 — (75,125) Phi derived: the a-only potential is refuted, corner signature confirmed as the real index
Opus lane derived Phi for Helali's (75,125) case (family F_2, j=1, corner
(5,20)->(7/5,2), reduced C-power pair (a,b)=(3,5), non-adjacent):
  Phi_75_125 = f*C^98 = -(1/9) y^201 (y^3+1)^101,
  signature (deg, ord_y, mult_{y+1}, cofactor deg) = (504, 201, 101, 202).
The ODE 15 c f' - 42 c' f = c^3 (c = y^2(y^3+1)) collapses under the forced
ansatz f = A y^5 g^3 to g = y^3+1, f = -(1/9) y^5 (y^3+1)^3 — the unique
polynomial solution of degree <= 14.
VERDICT on POTENTIAL_PROBE's a-only prediction (deg 550 / mult 69 at a=3):
DIFFERS, and unconditionally — deg Phi = 14 + 5N (deg C = a0 = 5 here vs 8),
so 550 is unreachable for ANY integer tower length N. The monotone potential
is NOT power-only; it must be indexed by the full corner signature Sigma,
confirming the third-review caution. A 3-point fit
  deg = (e*a0 - q + 1) + N*a0,  mult = e + N,  cofactor = r(e+N),
  N = a[t(a+b)-(kappa+1)] - 2b
reproduces all three known cases exactly; t and kappa remain correlated
(need a fourth corner with different t to separate).
Caveats [judgment]: the (5,20) polygon reduction is assumed type-II.b
(done in no paper); N=98 is formula-based pending an actual C-series build;
headline verdict robust to both. (72,108) is the resonance-gap exception
(r=0: pure-ansatz deg 10 < resonant 14).
Checker: phi_75_125_verify.py (50/50, from-scratch re-derivation + both
audited (8,28) corners as controls) — wired into run_tests.sh.
Files: PHI_75_125.md, phi_75_125.py, phi_75_125_verify.py.


## 2026-07-23 — Lane D: fourth AND fifth corners MATCH the fit; kappa = t-2 is structural
Two new out-of-sample Phi derivations by the PHI_75_125 method:
  F9 j=0 (56,84), corner (7,21), (a,b)=(2,3), t=7 (first t outside {4,5}):
    Phi = -(1/10) y^107 (y^5+1)^54, signature (377,107,54,216) — MATCHES the
    six-parameter fit'''s parameter-free prediction exactly.
  F2 j=0 (50,75), corner (5,20), (a,b)=(2,3):
    Phi = -(1/6) y^75 (y^3+1)^38, signature (189,75,38,76) — MATCHES.
Five points now separate corner-dependence from pair-dependence.
kappa = t-2 is PROVEN structural on the whole standard-chart class: the
Laurent chart (X,Y)->(x^-1, x^l y) has Jacobian -x^(l-2) (verified
symbolically for generic l) and t=l, so kappa=t-2 identically for all 15
length-1 A0'''=(1,0) families. The fit reduces to five parameters with
  N = a[t(a+b-1)+1] - 2b   (reproduces all five N values),
now overdetermined and validated at t in {4,5,7}. Possible kappa != t-2
escapes: A0'''=(2,0) families (F12,F13) and length-2 chains (F18-F24) —
composite charts exist in no paper; named follow-up (composition heuristic
kappa=l2-l1 is conjecture only). F14 (t=3) is the natural sixth point.
Judgment items mirror PHI_75_125.md; the t=7 point is LESS conditional
((b-1)/a=1 so the slice-index picture transfers verbatim).
Checker: phi_corner4_verify.py (40/40, controls included) — wired into
run_tests.sh. Files: PHI_CORNER4.md, phi_corner4.py, phi_corner4_verify.py.

## 2026-07-23 — Lane A: the k=6,7,8 window caps are PROVEN; bridge judgment flags retirable
The caps load-bearing under every bridge kill (ord >= 12k; deg <= 15k sub1 /
14k sub2; stripped 3k/2k; k=6/7/8 = dm2/dm3/dm4) are now recited and proved:
the three T3 section-3 valuation inductions close as symbolic identities in k
(k=6,7,8 consume genuine slices P2,P1,P0); D-transform arithmetic gives the
caps identically in k; the d3-killing shift is the polynomial identity
  D~_j = sum_m binom(m,m-j) D_m (-D_3/4)^(m-j)
(C4 exponents cancel, D~_3 = 0, caps preserved term-by-term). End-to-end
generic corroboration on random exact-Q slices: caps hold k=2..8 both regimes,
degree caps attained 7/7 rows (tight). No new conditionality beyond the
program'''s inherited premises [P1][P2][P3] — upstream facts loaded from
sha-pinned paper_src/upstream_facts.json, not hand-copied. The k=6,7,8
[judgment] flags in FULL_SYSTEM_BRIDGE.md section 4 and BRIDGE_SWEEP.md are
retirable to the audited k=2..5 trust tier (WINDOW_CAPS.md section 5).
Checker: window_caps_verify.py (81/81) — wired into run_tests.sh.
Files: WINDOW_CAPS.md, window_caps_verify.py.


## 2026-07-23 — Lane G: F14 (t=3) MATCHES; gap>0 probe confirms the resonance-gap story; SEVEN exact points
F14 j=0 ((66,231), corner (9,24), (a,b)=(2,7), t=3, kappa=1): ODE
6cf' - 34c'f = c^6 collapses to g = y^5+1, A = -1/10:
  Phi = -(1/10) y^165 (y^5+1)^42, signature (375,165,42,168) — the law's
parameter-free prediction, exactly. First t=3 point (coverage t in
{3,4,5,7}), first e=6, first q=4. Uniqueness by full 52-unknown linear
solve. The residual is again the 10th cyclotomic: residual is indexed by
dg = a0-q alone (F14, F9, (108,144) all have dg=5).
gap>0 probe (F1, r=0, gap=1): the (72,108) resonance-gap story GENERALIZES.
Generic-f solve gives unique f = (1/15) y^4 (y+1)^2 (4y-1) — pure-ansatz
shape times a unit cofactor of degree exactly gap (degree-1 analogue of the
(72,108) quartic). Phi = (1/15) y^205 (y+1)^69 (4y-1), signature
(275,205,69,1) — MATCHES the r=0-amended law. Unified law
  (pure+gap+N*a0, rho+Nq, e+N, gap+r(e+N))
now reproduces SEVEN exact points, two in the gap regime.
Two new checked mini-lemmas: gap = (q-1) - a0/t (gap=0 iff a0 = t(q-1) —
pure corner-data criterion for the resonance-gap regime); te-coef = 1 is
the kappa=t-2 identity in ODE form.
Honest remaining hole: gap>0 WITH r>0 (F3,F7,F10,F15,F16) has no derived
point — cofactor formula is conjecture there; F7 (42,147) smallest next.
Judgment: unreduced-polygon as before; F1's c = y^q(y+1) normalization is
by-construction, not ODE-forced (flagged); both new points in the
less-conditional (b-1)/a-integral class.
Checker: phi_f14_verify.py (37/37; all four previous corners as controls
plus audited (72,108) arithmetic) — wired into run_tests.sh.
Files: PHI_F14.md, phi_f14.py, phi_f14_verify.py.


## 2026-07-23 — Lane F: C20 KILLED via Galois descent; full 23-shape library census is a structural theorem
C20 (61X^2D^2+6XDE-11E^2, source-derived) reduces to 11r^2-6r-61=0,
residue discriminant class 170; solvable over F iff sqrt(170) in F. Fixed
q (S4, disc class 17): Obs = squarefree(170*17) = 10 != 1 -> KILLS (the
shared factor 17 cancels; obstruction survives as cofactor 10; a forcing
quartic with disc class exactly 170 is the unique S4-escape). All branches
carry verified witnesses: sharpness qs = y^4+7y^2-8y+10 (S4, disc =
170*52^2, obstruction vanishes, C20 solvable over split(qs)); branching
necessity qv = y^4-344y^2+28224 (V4, square disc, naive obstruction wrong,
(512x-x^3)^2 = 170*336^2 mod qv); A4 live test qa = y^4+8y+12 -> kills.
Library census (all 23 shapes, predicted by support geometry):
  QUADRATIC-OBSTRUCTION: 2 — exactly C08 (105) and C20 (170), both killing
  simultaneously for fixed q (joint subgroup {1,105,170,714} meets {1,17}
  trivially).
  LINEAR-in-primitive-character: 12 (C01-C03,C05,C06,C11-C17) — character
  surjectivity proven constructively (Bezout vector + explicit rational
  torus point per shape) -> solvable over every field; CONSTRAINT for any
  forcing quartic.
  MULTIDIM rank>=2: 9 (C04,C07,C09,C10,C18,C19,C21-C23) — rational
  witnesses re-verified.
  HIGHER-RESOLVENT-NEEDED: 0 — the pilot'''s C12 worry dissolves (linear in
  S^3/X^4, gcd(3,4)=1 -> surjective); higher resolvents matter only on
  constrained sub-loci (scoped, not claimed).
Case-compiler transfer: C08/C20 kills transfer to any family member by a
two-line check — Galois label + disc class of its forcing quartic (kill
iff S4/C4/A4 and disc class avoids {105,170}).
Checker: galois_library_verify.py (69/69; all 23 equations re-derived from
f31_graded.txt via the audited residue_lemmas machinery; witnesses''' every
used property verified from scratch) — wired into run_tests.sh.
Files: GALOIS_LIBRARY.md, galois_library_verify.py.


## 2026-07-23 — Lane I: case compiler pilot lands; three grounded cases reproduce exactly
case_compiler.py compiles a family case into a case-dossier-v1: corner
signature (kappa=t-2 substituted), unified-law Phi prediction with
data-driven REGIME_STATUS + conjectural flags, BOTH presentations
(eliminated-f31 / pre-resultant G-system) with instantiated-vs-schematic
markers per field, the GALOIS_LIBRARY two-line transfer check (Galois
label + disc class vs {105,170}), per-mechanism transfer inventory, and a
judgment section. Canonical JSON is deterministic (byte-identical across
recompiles).
Validation (case_compiler_verify.py, 71/71): (72,108) reproduces the
audited facts exactly — signature (238,204,30,4), regime resonance_gap_r0,
forcing quartic S4/class 17, both kills AUDITED, both presentations
INSTANTIATED ((D,F,E,s)=(31,7,21,3), G-weights 156/168/180/204).
(75,125): (504,201,101,202), forcing H=y^2-y+1 (C2, class -3), transfer
CONDITIONAL, presentations schematic. F9: (377,107,54,216), forcing 10th
cyclotomic labeled C4 by the generic router, cross-checked constructively.
Conjectural flags fire exactly on F7 (gap>0&r>0), F12 (A0'!=(1,0)), F4
(k=2); all 17 family Diophantine identities check.
Transfer inventory: AS-IS — residue library + Galois descent, corner law
(standard-chart class), modular triage. PARAMETRIC — cascade engine, T2
squeeze, bridge, window caps, divisor lemmas. METHOD-ONLY — S-unit corner.
New automatic flag: dg-even cases (F7 residual y^2+1) break the g(-1)=0
mechanism every derived point uses — marked UNKNOWN, not predicted (watch
Lane H). REGIME_STATUS is the table to flip when F7 lands.
Checker wired into run_tests.sh. Files: case_compiler.py,
case_compiler_verify.py, CASE_COMPILER.md, case_dossier_GGHV_72_108.json,
case_dossier_F2_j1_75_125.json, case_dossier_F9_j0_56_84.json.


## 2026-07-23 — Lane H: gap>0,r>0 regime DIFFERS — old conjecture refuted, RAMIFIED law derived at four corners
F7 (42,147) derived signature (250,165,83,2) vs old-law (250,165,42,43):
deg/ord match on every branch; mult/cofactor conjecture REFUTED. Same split
at F3 (189,112,75,2), F16 (528,407,117,4), F10 (1917,820,1093,4).
Why it had to fail (PROVEN at dg=2): the polynomial-solvability eliminant
factors completely (F7: E ~ g0^27 (3g0-2g1^2)(4g0-g1^2)^6), leaving only a
RAMIFIED branch (g=(y+1)^2, double root forced to -1) or a complex-pair
branch (no real root at -1); a simple root at -1 — what mult=e+N needs —
is impossible; independently mult=e+N > res >= deg f by degree count.
Amended ramified law (four exact points): deg = res + N*a0 and
ord = rho + N*q survive; mult = dg(e+N) - (dg-1); cofactor = gap + r —
retro-explaining the audited (72,108) quartic (4 = gap + r). F16 (gap=3)
was the separating experiment for the unit degree. Mini-lemma: dg is even
exactly on the gap>0,r>0 rows — the regime split is corner-data-visible.
Judgment: branch selection by continuity with (y+1)|C (complex-pair branch
recorded, not excluded — no C-series exists in any paper); (b-1)/a
non-integral for F10/F16; branch completeness proven only at dg=2; F15 not
attempted. Twelve derived/audited points total; a mult_g-indexed
unification is recorded as observation only.
CASE COMPILER UPDATED (first live regime flip): REGIME_STATUS
gap_pos_r_pos -> grounded; law_signature branches to the ramified formulas;
dg-even H_note now states the proven ramified shape. Compiler now predicts
F7 at exactly the derived signature, branch-annotated; verifier check F
updated accordingly; dossiers regenerated (71/71).
Checkers: phi_f7_verify.py (62/62), case_compiler_verify.py (71/71) — both
in run_tests.sh. Files: PHI_F7.md, phi_f7.py, phi_f7_verify.py + compiler
updates.


## 2026-07-23 — Prior-art deep audit: corner-law claims are NEW; literature contains a zero-freedom postdiction
Full-text read of the GGV corpus (GGV1/GGV3/GGV5/GGHV22 TeX, 1310.8249
complete, 2506.05697 + five more at theorem level; Moh 1983 scans
unreadable — caveat recorded; Horruitiner thesis unfetchable, siblings
read). Verdicts: closed-form Phi + signatures NOT FOUND (new); signature
law + N-formula NOT FOUND even partially (new); parametric ODE family /
kappa=t-2 theorem / ramified dichotomy PARTIALLY (instance layer only);
resonance gap PARTIALLY (one instance, no law).
DECISIVE FIND (GGHV22 sec.4, L1571-1597): killing the (66,99)/(9,24)
corner they solve 6*C3*f1' - 10*C3'*f1 = C3^2 (C3=y^8(y+1)) by CAS,
print f1 = -(1/910) y^9 (y+1)^2 (243y^4-81y^3+54y^2-42y+35) — exactly the
parametric ODE at (a,b,t,kappa)=(2,3,3,1), kappa=t-2 holds, and the
f-level corner law postdicts its FULL signature with zero fitting freedom
(ord 9, mult 2, deg 15, unit cofactor degree gap+r=4, separable, coprime
to y(y+1)) — the published record confirms the law out of sample. Their
D_k transformation (L1612-1632) is the ancestor of ours. (72,108) itself
is left open in print (L268), no ODE or quartic published.
Required citations recorded in PRIOR_ART.md (GGHV22 sec.4 is the critical
one; GGV1 bracket ODE + multiplicity lemma; 1310.8249 Abel-ODE genre).
New checker: prior_art_postdiction_verify.py (12/12: family membership,
uniqueness, coefficient-exact match to the printed f1, law postdiction,
cofactor properties) — wired into run_tests.sh.
Files: PRIOR_ART.md, prior_art_postdiction_verify.py.


## 2026-07-23 — Lane J: kappa = t-2 EXTENDS to every escape family; the survey boundary moves from chart to model
Composite-chart construction for the named escapes (A0'=(2,0): F12, F13;
length-2: F18-F24): every chain transformation is a root-shift shear
(Jacobian 1, bracket-preserving), the chain performs exactly ONE final
inversion, and the composite fuses to (X,Y) = (x^-1, x^l y + sum lambda_i
x^{e_i}) with Jacobian -x^(l-2) FOR ANY shear terms (verified
symbolically; fused exponents integral >= 0 on all five escapes). The
kappa = l2 - l1 heuristic is REFUTED — it is the Jacobian of a
two-inversion chart the chain never makes. Escape kappas: F12: 2, F13: 1,
F22: 2, F23: 2, F24: 6 — all t-2.
Bonus: F18-F21 are DEAD (GGV5 itself proves no standard (m,n)-pair
realizes them; core re-verified) — the only real length-2 escapes are
F22-F24; F22/F23 use only one transformation (same chart as F12);
chain-step formula A_next = lower + q*(s/l,1) reproduces every table
corner (10 families).
The real boundary: pure-power defect zeta > 0 breaks the commuting-tail
structure by exactly [x^zeta C^a, C^b] = zeta*b x^(zeta-1) C^(a+b-1) C_y
!= 0. kappa settled; zeta-corrected tail theory is the named follow-up.
F12 Phi (CONDITIONAL-zeta): first dg-odd point in gap>0,r>0 — branch
variety is two points (g1 = (-1 +- sqrt(22))/8), both squarefree: NO
ramified branch exists at dg=3, mirror of the dg=2 obstruction. So
dg-parity has a MECHANISM: ramification is root-availability (even dg
forces double/complex roots; odd dg keeps a simple real root).
Conditional signature (814, 506, 102, 206) — OLD-law shape. F22 lands at
gap=-1 (no derived point in any negative-gap regime — flagged, nothing
claimed).
Compiler updated: A0'!=(1,0) conjectural reason now states chart settled /
CONDITIONAL-zeta (still fires, honest reason); dossiers regenerated
(71/71).
Checker: composite_charts_verify.py (41/41; resultant-route eliminant,
exact Q(sqrt(22)) arithmetic, length-1 Jacobian control) — wired into
run_tests.sh. Files: COMPOSITE_CHARTS.md, composite_charts.py,
composite_charts_verify.py.


## 2026-07-23 — Lane C: 45/49 HUNT-cell states KILLED (pending audit); 13/17 cells fully-forced dead
The 17 s-unit BM-candidate cells (11 T2 d1=d2==0, 6 T1 d2=sigma==0, sub2
AND sub1) carry 49 fully-forced states. Census: 45 KILLED (37 at depth 2,
8 at depth 3 — exactly the DIVISOR_LEMMAS section-6 prediction), 4 OPEN.
All 11 T2 cells fully dead; the 4 survivors are deg-d1=6 T1 states (both
windows x a9_b1000_T1, a8_b1100_T1) — the known [J6] grevlex blowup at
600s/1500s budgets; exact systems recorded; msolve pass is the next dent.
Enabling idea: CLASS-POLYNOMIAL root reconstruction — group q-roots by
joint exponent profile; each class is an unknown monic factor of q/2048,
largest class = exact quotient, division-remainder coefficients = the
relations. q squarefree => one unit ideal kills every Galois assignment of
every admissible split at once over Q, <=2 unknowns in practice. Dissolved
every heavy pattern incl. the deg-sigma=8 weight-scaling states (the
degree-250 homogeneous point died at depth 3, no subsum tree needed).
HONEST SCOPE: these cells live in the standard-regime windows — closes
queue item 4's depth-2 residue work / s-unit residual layer; does NOT
directly close any of the 24 alt branches (need alt-window bridge analog).
Cross-check: one prior phase_f2_sub2 kill independently reproduced; 44 new.
Audit trail: every kill records reconstructed polynomials, class
relations, saturation, and all master coefficients as exact strings in
alt_hunt_results.json — spec-only re-derivable from f31_graded.txt.
Files: alt_hunt_depth2.py (resumable), alt_hunt_results.json, ALT_HUNT.md.


## 2026-07-23 — Lane B: R9 symbolic elimination CERTIFIED; sweep verdict an honest structural NEGATIVE — wall localized
Elimination artifact: every G-system generator is linear in dm4 with pivot
3*dm1 = 3e != 0; dm4-free combinations H2 = dm1*G2 - dm2*G1, H3 = dm1*G3 -
dm3*G1, H5 = dm1*G5 + (d0*dm1+d1*dm2+d2*dm3)*G1 (5/6/12 terms, weights
228/240/264) live in the G-ideal by explicit cofactor certificate,
re-verified by exact expansion; divisibility lemma monic(e) | dm2*dm3 from
G1. Spare unknowns 45 -> 28 per state. Instantiation validated
BYTE-IDENTICAL against the independent fsb.augment route on R9 z=1; every
state build asserts exact cancellation of all 17 dm4 coefficients. Cached:
r9_eliminated_system.json (rebuild+recheck seconds: r9_symbolic_elim.py).
Sweep census (11 attempted; 45s triage + 300s exact caps;
r9_symbolic_sweep.json): R9 z=1,2,3 COST; z=4..6 not attempted (recorded
truncation); a9 T2 deg_e=10 8 cheapest states COST x8; 82/90 not attempted
(recorded truncation). ZERO kills, ZERO PROPER — no survival signal.
Decisive control: first batch state re-triaged at the bridge's own 300s
budget on p=10007: TIMEOUT — the reduction did not cross the feasibility
boundary. The deg_e=10 wall survives dm4 elimination on both families; it
is now localized to the dm2/dm3 ansatze (28 scalars) coupled to window-cap
e. Named continuation (R9_SYMBOLIC.md 3c): divisibility-forced valuation
split — v_{y+1}(dm2)+v_{y+1}(dm3) >= 9, v_r(dm2)+v_r(dm3) >= 1 — 20
exhaustive cases of 18 spare unknowns each, composing with the certified
H-system unchanged.
Ops: orphan-proof runner (WSL-side timeout + ulimit -v 8G), final WSL ps
clean; mid-lane defect (slow Poly builder) found and fixed honestly with
equation-set-identical validation before any verdict relied on it.
Files: r9_symbolic_elim.py, r9_eliminated_system.json,
r9_symbolic_sweep.py, r9_symbolic_sweep.json, R9_SYMBOLIC.md.


## 2026-07-23 — Lane M: the 4 [J6] states KILLED by msolve — 49/49 HUNT states dead, 17/17 cells CLOSED(all0)
All four deg-d1=6 T1 states (sub{2,1} x a9_b1000_T1_sz1_dz1#state3,
a8_b1100_T1_sz1_dz1#state0) KILLED at depth 3 (PENDING AUDIT), msolve
walls 0.5-4.5s per call, ~20s per state end to end. The s-unit
BM-candidate residual layer is fully gone. Per ALT_HUNT [J5],
CLOSED(all0) does not close the enclosing branches.
Replay guarantee: reconstructions replayed string-identical to Lane C's
recorded polys before anything reached msolve; systems attacked identical
to Lane C's generator sets.
STRUCTURAL FINDING (why [J6] blew up): the depth-2 truncations are
GENUINELY SATISFIABLE — msolve returns isolated solutions at depth 2 for
all four; GB engines were not slow at seeing a unit ideal, they were
computing a full GB of a satisfiable system. The depth-3 master
coefficient is the killer. On the depth-3 killing system, sympy Buchberger
and Singular std over Q still time out at 240s while msolve decides in
~1s. Corroboration: Singular std over F_p (10007/10009/100019) unit ideal
12/12, <=1s each — independent second implementation agrees. Full
end-to-end rerun reproduces KILLED x4. WSL sweep clean.
Audit note: char-0 std timing out means Singular-lift certificate
extraction for these 4 may need longer budget or an msolve-assisted route.
Files: j6_msolve.py, j6_msolve_results.json, J6_MSOLVE.md.


## 2026-07-23 — Lane K: zeta-corrected tail theory lands; the mu-LADDER unifies the laws and corrects a landed claim
One-slice rigidity theorem: the zeta-defect cannot stay in the tail (the
top commutator slice sits strictly above x^kappa, uncancellable); unique
repair puts it in the element D = x^eta C, eta = zeta/a, giving the same
parametric family with T = t+eta (kappa = t-2 chart-fixed) — the family
moves off the diagonal by exactly eta. Universal lemma: eta=-1 is ALWAYS
dead (integrates to (f/c^e)' = 1/(aTc), polynomial iff all residues of 1/c
vanish — impossible).
F12: every motivated defect 0<|eta|<=1 is NO-POLYNOMIAL-SOLUTION; the two
surviving models are eta=0 (conditional) and an eta=+2 canonical collapse
(T=6, forces g=y^3+1, N=157, signature (1292,806,162,324)); model
selection needs the actual polygon reduction. F13: only eta=0 survives.
THE MU-LADDER (corrects COMPOSITE_CHARTS section 5): branches are graded
by mu = mult_{y+1}(g). F12's standard family realizes ALL rungs: mu=1
(Lane J's points, (814,506,102,206)), mu=2 NEW (g=(y+1)^2(y-beta), beta a
root of 195b^4+120b^3-40b^2+32b-80, (814,506,203,105)), mu=3 NEW
(g=(y+1)^3, (814,506,304,4) — exactly the formulas COMPOSITE_CHARTS
declared NOT realized; that claim was a g^e-uniform-ansatz artifact, now
corrected in place). mu-graded law: mult = mu(e+N)-(mu-1),
cof = gap+r(e+N)-(mu-1)(e+N-1) — specializes identically to the
unramified law at mu=1 and PHI_F7's ramified law at mu=dg (via r = dg-1).
PHI_F7 judgment-6 observation is now the actual branch structure. Parity
refinement: even dg kills mu=1 (that theorem stands); odd dg does not
exclude higher rungs. Twelve corner-law points untouched — the mu-graded
law is a strict generalization all of them obey.
Compiler F12 flag updated to cite the enumerated surviving models;
dossiers regenerated (71/71). Named follow-up: F10 (dg=4) likely carries
intermediate rungs — one small solve each.
Checker: zeta_tail_verify.py (34/34; fresh mu=2 rebuild, exact
mod-quartic residuals) — wired into run_tests.sh.
Files: ZETA_TAIL.md, zeta_tail.py, zeta_tail_verify.py.


## 2026-07-23 — Lane R: adversarial review of zeta-tail/mu-ladder — NO REAL HOLE; five patches applied
Independent skeptic re-derived the load-bearing mathematics, including the
F12 mu=2 quartic by a DIFFERENT elimination route (augmented-matrix minors
+ gcd): same quartic up to scalar, signature (814,506,203,105) exact; mu=3
also rebuilt. Verdicts: rigidity SOUND (verified at three (a,b) x two t
beyond the author'''s instance; wording fixed tail-tail -> head-head);
eta=-1 lemma SOUND (all residues reproduced via sympy.residue;
squarefreeness not needed); mu-ladder SOUND with findings (r = dg-1 is
DEFINITIONAL so the specialization is exact algebra not discovery;
deg-ord = (e+N)dg+gap identically so cof-law == mult-law — one datum per
rung); cross-doc consistency: correction block accurate, compiler flag
completed (eta=+2 mu-support {1,3}, eta in {-2,-3} OPEN rows added).
THE ONE GAP: "dg even kills mu=1" was overreach — branch-completeness is
proven only at dg=2 (PHI_F7 judgment 5); now SCOPED (proven dg=2,
conjectured even dg>=4 pending F10 enumeration — Lane Q running) in
ZETA_TAIL.md and COMPOSITE_CHARTS.md. Defensive guard added to
zeta_tail.py'''s sub-resonant-evasion branch (verified to fire on no
current row; now refuses to conclude silently if a future row trips it).
All checkers re-verified post-patch (zeta 34/34, compiler 71/71).
Review artifacts: REVIEW_ZETA_MU.md, review_zeta_mu.py (31/31).

## 2026-07-23 — Lane L: R9 valuation split — WALL SURVIVES; obstruction is formulation-level
Exhaustiveness MACHINE-CHECKED: the 20-case R9 split covers all 5145
achievable valuation profiles (batch shapes a=7/8/9: 8/9/10 cases,
6958/6126/5295 profiles); degenerate dm2/dm3=0 representable everywhere.
Implementation finding: structured-product ansatz densified convolutions
(>25 min/case, measured, killed); replaced by the equivalent LINEAR
derivative-vanishing formulation (v_{y+1}(p)>=m iff p(-1)=...=p^(m-1)(-1)
=0) on Lane B'''s cached dm4-eliminated builds — spare dimension 28->18.
Census: z=1 4 full cases COST x4 (45s x 3-prime triage TIMEOUT uniform);
2 partial (outage); DECISIVE CONTROL: case (9,1) — the structurally
easiest, 10/13 dm2 coefficients forced — TIMEOUT at the bridge'''s own
300s budget on both completed primes. Remaining cases/columns NOT
ATTEMPTED (honest truncation). Zero kills, zero PROPER.
VERDICT: two successive sound reductions each strictly shrank the system
without crossing the mod-p feasibility boundary — the swell is driven by
the generic dm3 tail coupled to window-cap e, not spare-dimension count.
Named next (R9_VALSPLIT.md section 4): msolve/F4 on a split case (the
[J6] lesson), e'''s unknown-tail divisibility content, or linear-algebra
on the H-system coefficient matrix in forced-valuation coordinates.
Files: r9_valsplit.py, r9_valsplit_results.json, R9_VALSPLIT.md.


## 2026-07-23/24 — Codex trust round: 49/49 HUNT kills independently audited (0 disagreements); f37 replay scripts; public CI live
Independent-author audit of all 49 HUNT/J6 kills (audit_alt_hunt_kills.py,
codex-authored, imports no project code, own f31_graded.txt parser +
homogeneity checks, own class-factor reconstruction, subprocess GB with
120s hard kill): 41 FULLY-VERIFIED (independent sympy GB reproduces the
unit ideal), 8 VERIFIED-DATA-ONLY (honest GB timeout — includes the 4 J6
states known to need msolve), 0 DISAGREEMENTS, 0 UNPARSEABLE; exit 0.
The zero-disagreement streak survives a full different-author replay.
Census: audit_alt_hunt_census.json.
f37 INDEPENDENT REPLAY (codex-authored, per GPT-Pro review 4 requirement):
f37_replay_m2.m2 (Macaulay2) + f37_replay_sage.py regenerate the
pre-resultant generators FROM SCRATCH (native convolution, no
generators.json consumption), rebuild f31 via the resultant chain, and
check principality + equality-up-to-unit + membership. Sympy self-test
(f37_replay_selftest.py): 8/8 — the independent construction equals
generators.json exactly. Execution of M2/Sage awaits install (F37_REPLAY.md
has instructions); after one green run, the f37 standalone note (short
Version A) is unblocked.
PUBLIC REPO: v0.2 content live (48b85e1) + CI added (834bbfc, pushed):
fast-checks.yml on push/PR (clean-clone guard, D3, f37 certs, field-split
ledger, audit_inf_cases; measured 12m10s, all PASS locally) + weekly/manual
full-suite.yml. GPT-Pro review 4 folded: contact-draft scope fix,
CASE_COMPILER.md refresh, landscape correction (CAOS/alok/Naskrecki/
eumemic) in memory. Spend-limit casualties: Lane Q (F10 rungs, no artifacts
— re-dispatched to codex), Lane N (alt bridge, artifacts on disk
uncommitted, resume pending), local auditor + overnight batch (relaunched).


## 2026-07-24 — F10 mu-rungs (codex lane): parity PROVEN at dg=4; new mu=2 point MATCHES; the ladder has selection rules
Exact enumeration of all F10 rungs (dg=4): mu=1 EMPTY — every cubic-h
partition eliminated by exact polynomials with zero real roots, so the
even-dg parity claim is PROVEN at dg=4 (the Lane R scoping gap closes:
now proven at dg=2 AND dg=4; dg>=6 has no gap>0,r>0 survey corner).
mu=2 REALIZED — signature (1917,820,547,550), mu-graded law MATCH
(thirteenth exact point); squarefree component minimal polynomial
192p^6-720p^5+586p^4-668p^3+1728p^2-1510p+45 (two real roots), plus a
double-root component 12z^4+15z^3-5z^2+15z+12 (two real roots).
mu=3 EMPTY — obstruction 10s^4-20s^3+8s^2+2s+5 irreducible with zero real
roots. mu=4 reproduces PHI_F7's (1917,820,1093,4) with resonant
coefficient 3740*F27 = 2401.
STRUCTURAL FINDING: the mu-ladder has SELECTION RULES — F12 realizes all
three rungs, F10 realizes only {2,4} — rung support varies by corner,
strengthening the finite-regime-grammar picture over any universal law.
Docs upgraded: ZETA_TAIL.md + COMPOSITE_CHARTS.md parity scope now cites
the proof. Checker: mu_rungs_f10_verify.py (39/39, exact ODE
substitutions, F12 control) — wired into run_tests.sh.
Files: MU_RUNGS_F10.md, mu_rungs_f10.py, mu_rungs_f10_verify.py.


## 2026-07-24 — AUDIT ROUND landing: 20 certificates CERTIFIED end-to-end; ledger + roll-up regenerated; reconciliation clean
The certificate pipeline's first full end-to-end pass: 20 kills now carry
cofactor certificates (1 = sum c_i f_i) verified by the certificate-
consuming ninth auditor (pure sympy expansion, no GB engine):
AUDIT CENSUS CERTIFIED=20 EXPANSION-FAILED=0 UNCERTIFICATED=29. Zero
disagreements. Reconciliation: producer claims (11 first batch + 9
parse-fixed rerun) = 20 = disk = auditor exactly (an earlier '29 found'
count was a grep artifact matching census summary lines; no certificates
were lost). The 29 uncertificated split: 15 genuine lift timeouts (0/15
cracked even at 900s — structurally hard lifts, msolve-assisted route
dispatched) + 14 fallback-route expansion failures (std-transform lift
composition bug, distinct from the fixed GEN-shadowing bug; on the fix
list). Ledger + roll-up regenerated: unmapped 0, real conflicts 0,
ambiguous 31, redundant re-attempts 27; branch closures a11_b1111_T1 8/8,
a12_b1110_T2 8/8, a14_b0000_T2 1/1 all CLOSED-pending-audit.
Overnight cert batch stamped DONE 23:21 (honest zero at 900s).
kill_certificates/ (49 per-kill JSONs incl. the 20 certified) committed.


## 2026-07-24 — Lane N (collected by Opus): alt-window bridge EXISTS as construction; computational wall — formulation-level again
Construction verdict: EXISTS (derived + mechanically spot-checked). The
alt analog is the landed G-system augment(state, regime="sub1") with sub1
stripped caps (dm2/dm3/dm4 <= 18/21/24, 66 spare scalars) + marked-root
adjunction q(r)=0. Load-bearing novelty [J1]: BRIDGE_SWEEP section-3's
"not applicable" skip used the SUB2 cap (deg e <= 10) but alt states are
SUB1 states (deg e = 15 attained); WINDOW_CAPS' k=6,7,8 caps are
t-regime-free. Soundness chain re-verified independently by the collector
(homogeneity 156/168/180/204, phi_identical, caps, complement exact mod q;
fast builder equation-set IDENTICAL to fsb.augment).
Pilot census: 0 KILL, 0 PROPER, 2 COST — both deg-d2=6 alt stragglers
time out on EVERY engine (Singular mod-p x3, exact Q, numroot x3, msolve
char-0 ~300s each). CRITICAL HONEST FINDING: the control a12_b1110_T2 is
a KNOWN kill on the direct/f31 route, and the bridge did NOT reproduce it
— the bridge system is strictly harder in the alt regime (the +66 spare
unknowns + marked-root nonlinearity put it in the known swell family).
Inverse of the standard-window pilot (1-9s kills).
Named next (ALT_BRIDGE.md section 4a): symbolically eliminate the 66
spare unknowns BEFORE solving (the Lane B pattern); secondary numroot
recon. Alt frontier UNCHANGED by this lane; zero kills (nothing to audit).
Files: ALT_BRIDGE.md (finalized), alt_bridge.py, alt_bridge_results.json,
alt_bridge_results_symbolic.json.


## 2026-07-24 — TRANSFER TEST phase 1 (Opus lane): the (75,125) C-series is BUILT; N = 98 DERIVED, judgment retired
The actual C-series + D-transform tower for (75,125) constructed from
corner data. N read off the built tower: Phi is the whole integer slice
M = b*t + j = 36 of S^b; the SLICE-SUM INVARIANT — every monomial of a
slice carries the uniform factor c^(a*M-b), because the D-weight
a(t-k)-1 summed over a slice collapses to a*M-b independent of the
individual k_i — gives clear = 103, N = clear - b = 98. This dissolves
the exact concern behind PHI_75_125 judgment item 3 (the non-integral
per-term index (b-1)/a = 4/3 is per-term; Phi is the slice): the
corner-144 N-formula is re-derived structurally, not extrapolated.
Judgment item 3 RETIRED in PHI_75_125.md.
Phi emergence CONFIRMED: the tower yields -(1/9) y^201 (y^3+1)^101,
signature (504,201,101,202) exactly. Controls: same machinery reproduces
(72,108) (N=28 -> (238,204,30,4), audited ground truth) and (108,144)
(N=67 -> (550,205,69,276)).
Remaining boundary: the unreduced-polygon assumption (shared with all
derived corners). Later phases (D-transform G-system, window/pre-resultant
layer on (75,125)) are the named next steps of the transfer test.
Checker: c_series_75_125_verify.py (36/36) — wired into run_tests.sh.
Files: C_SERIES_75_125.md, c_series_75_125.py, c_series_75_125_verify.py.


## 2026-07-24 — Alok cross-check (Opus lane): exact setup corroboration; regime-disjointness quantified; 0 findings
The independent public program alok/jacobian-two (commit ded8e67, docs/
newton-72-108-sparse.md + scripts/newton_72_108.py L118-127) proves for
the GGHV (72,108) polygons: Case 1 >= 3 / Case 2 >= 4 total nonzero
strict-interior coefficients (exact exhaustion, 7504/3683 patterns, 5
GB=[1] certificates with recorded sha256s).
Correspondence EXACT at the setup level, verified vertex-for-vertex:
their CASE_1 = our sub1 (corners (0,8),(0,12) added), CASE_2 = our sub2;
identical [P,Q]=x^2. Our checker independently re-derived their interior
lattice census from the vertices and reproduced their table exactly
(Case 1: P=35, Q=87; Case 2: P=7, Q=21) — two independent programs are
provably attacking the identical configuration.
Per-state support mapping: NOT soundly recoverable (their (x,y)-monomial
counts vs our y-degree/zero-flag D-transform coordinates are related by a
nonlinear convolution; no per-state formula invented — honesty over
coverage).
Cross-check verdict: CORROBORATION, 0 findings. Sound guard over ALL live
states (sub1 171 branches / 44117 states; sub2 26 / 7888):
live-below-threshold = 0. No open branch of ours sits in their
proven-dead sparse region. Regime-disjointness quantified: their theorem
closes the sparse floor (<=2/<=3 interior terms), our cascade attacks the
dense body (every tracked state carries forced Phi != 0, above threshold
by construction) — complementary, not overlapping; they contradict us at
no live branch. First concrete act of the exchange-variable-maps
collaboration mode.
Checker: alok_crosscheck.py (exit 0 iff no live-below-threshold state) —
wired into run_tests.sh. Files: ALOK_CROSSCHECK.md, alok_crosscheck.py.


## 2026-07-24 — TRANSFER TEST phase 2 (Opus lane): the (75,125) G-system BUILT; the a>=3 window-cap boundary characterized
The pre-resultant G-system for (75,125) constructed by the (72,108)
recipe with derived parameters: linear window S^3 (20 eliminations),
forcing window S^5 -> 10 quintic generators G1..G9,G11 (slice j=10
skipped, its window unknown pairing with the skipped linear slice —
system closed). Spare inventory: 9 spares dm2..dm10 (formula (a-1)t-1;
the jump 3->9 is purely a: 2->3). Weights: intrinsic-u arithmetic
progression 26..34,36 (Phi at 36; 35 = skipped G10), every generator
verified weight-homogeneous. Builder control: reproduces the published
(72,108) G1,G2,G3,G5body+Phi BIT-FOR-BIT with weights 156/168/180/204.
THE OBSTRUCTION (boundary as result): the physical y-order normalization
W_step = ord_y(Phi)/M is 12 (integer) at (72,108) but 201/36 = 67/12
NON-integral at (75,125) — the exact integer y^W stripping behind the
(72,108) bridge substitution does NOT transfer; the window-cap layer is
QUASIPOLYNOMIAL at a>=3 (quasi-period 12 here), exactly the boundary
CORNER_144_COMPARISON section 5 flagged (denominator 5 there). The
abstract ideal transfers; the bridge layer needs a quasipolynomial cap
theory — the named phase-3 obstacle.
Phi-consistency: CONSISTENT in the intrinsic grading (Phi enters G11 at
slice 36 with forced weight; clear=103, N=98 matching phase 1;
ord_y(Phi) = W_step*M as a rational identity).
g_system_75_125.json (g-system-v1) instantiates exactly the fields the
case compiler marks schematic -> structure-instantiated with
window_caps OBSTRUCTED(a=3 boundary).
Checker: g_system_75_125_verify.py (47/47; (72,108) control +
anti-fabrication re-derivation of G1,G2) — wired into run_tests.sh.
Files: G_SYSTEM_75_125.md, g_system_75_125.py, g_system_75_125.json,
g_system_75_125_verify.py.


## 2026-07-24 — Makar-Limanov restriction check (Opus lane): INAPPLICABLE across the frontier; paraphrase corrected
Paper located: L. Makar-Limanov, "On the shape of a counterexample to the
two-dimensional Jacobian conjecture," Serdica Math. J. 51 (2025) 299-314,
DOI 10.55630/serdica.2025.51.299-314. CLOSED ACCESS — no arXiv/galley
reachable; the tested lemma is the review's paraphrase made precise and
SELF-CERTIFIED computationally (stated loudly in ML_RESTRICTION.md §1;
obtaining the PDF to upgrade to verbatim is the recorded open item).
Finding on the lemma itself: the paraphrase AS LITERALLY STATED is false
(counterexample J(x+y, xy+y^2) = x+y on the diagonal w(x)=w(y); also
fails non-primitive weights). Verified correct hypotheses: positive +
primitive gcd(w(x),w(y))=1 + non-diagonal. Under those: holds with no
exception in tested range.
Per-edge verdicts: all 27 objects INAPPLICABLE — the decisive geometry:
with the origin a vertex, NEITHER (8,28)-reduced polygon has any strictly
positive-weight edge (every strictly positive weight maxes at the lone
monomial vertex); the principal common-root edge has weight (1,0) (not
positive) and J(R^2,R^3) = 0 exactly — our leading pairs live in the
J=0 (algebraically dependent) regime ML handles separately, never the
J(rho,tau)=rho regime the lemma restricts. Same for all 17 survey
families (chart Jacobian -x^(l-2) => bracket monomial => J in
{0, monomial}).
Frontier impact: NONE — no branch removed, promoted, or contradicted; a
clean negative literature check that locates our frontier relative to the
newest published shape restriction.
Checker: ml_restriction_check.py (exit 0 iff no inconsistency either
direction; re-checks polygon data against upstream_facts.json) — wired
into run_tests.sh. Files: ML_RESTRICTION.md, ml_restriction_check.py.


## 2026-07-24 — Proof-gate hardening (Opus lane): 532 asserts converted across 29 checkers; artifacts byte-deterministic
The python -O vulnerability was REAL and is now closed: control
demonstration showed an original bare-assert checker returning rc=0 under
-O with a CORRUPTED f31_graded.txt (false pass), while the hardened
version correctly exits 1. 532 proof-critical assert statements across 29
checker files converted to explicit _require(cond,msg) gates (AST-based
transform with byte-offset splicing; BOM and inline/multiline asserts
handled; 0 residual asserts; _require count == original count per file;
check-count report lines byte-unchanged). audit_*.py already used raise
AssertionError (-O-safe, untouched); 16 checkers had no asserts. All 29
edited checkers re-run rc=0. Three -O gate demonstrations on corrupted
scratch copies: all fire correctly.
Determinism: state_kill_ledger.py, frontier_rollup.py, frontier_gen.py
stripped of wall-clock/mtimes -> git commit + SHA-256 of every source
artifact + schema version; SOURCE_DATE_EPOCH honored. Regenerated twice:
state_kill_ledger.json, FRONTIER_V2.md, FRONTIER.md all BYTE-IDENTICAL.
Full inventory: HARDENING_NOTES.md.


## 2026-07-24 — Alt spare-unknown elimination (Opus lane): dm4 eliminated (66->41), control STILL COSTs — wall is deeper than spare count
The certified R9 dm4-elimination transfers VERBATIM to the alt bridge
(generators are regime-independent; sub1/sub2 differ only in caps). Greedy
pivot probe: dm4 is the UNIQUE spare with a guaranteed-nonzero linear
pivot (3*dm1 = 3e); all 25 dm4 scalars eliminated symbolically at zero
soundness cost (cofactor certificates re-verified by exact expansion;
homogeneity 228/240/264; divisibility identity checked). No further
guaranteed-pivot elimination exists.
CONTROL VERDICT: COST — the known a12_b1110_T2 kill is NOT reproduced on
the reduced system (187 eq / 41 spares): numroot mod-p TIMEOUT x3, char-0
Singular std TIMEOUT 300s. Identical in shape to R9_SYMBOLIC's negative:
shrinking the system does not cross the feasibility boundary. PRIZE
(a11_b3100_T2): COST. Sweep not reached (gated on control). 0 KILL,
0 PROPER, 2 COST. Alt frontier unchanged.
Wall characterization (now shared verbatim with R9): the swell lives in
the dm2/dm3 CUBIC ansatz (H cubic in 41 spares, e to the 4th) coupled to
window-cap-attaining e (deg 15) + marked-root nonlinearity. Named next:
the divisibility valuation split (R9 sec 3c, ~18 unknowns/case), set up
here by the certified monic(e) | dm2*dm3 lemma, not yet implemented for
alt.
Engineering finding: the reduced system is intractable to BUILD in sympy
(one state > 30 min CPU); pushing all expansion into Singular via
coeffs(H_i, y) made it tractable — and re-hit the documented gm^8/N
parser trap (fixed by rational-coefficient-first emission). msolve
unusable here (needs pre-expanded input). Exact engine of record:
Singular-native std over Q.
Files: alt_elim.py, alt_eliminated_system.json, alt_elim_results.json,
ALT_ELIM.md. Orphan discipline verified clean.


## 2026-07-24 — GPT-Pro review 6 folded: TWO SCOPE CORRECTIONS + the F2 FAMILY THEOREM (first symbolic family closed form)
CORRECTION 1 (real vs complex): the F10 dg=4 parity result is a
REAL-LOCUS theorem, not a complex one — the eliminants have no real
roots, but the recorded nonreal roots give formal complex mu=1 branches.
Complex parity remains proved only at dg=2 (PHI_F7's complete
factorization). Rescoped in MU_RUNGS_F10.md (correction block),
ZETA_TAIL.md, COMPOSITE_CHARTS.md. Standing rule adopted: classify ODE
branch schemes over Qbar first; real-locus info is branch-selection
annotation; Sturm counts are NEVER complex kills.
CORRECTION 2 (the boundary is not a>=3): W_step = 75/21 = 25/7 is already
non-integral at a=2 (F2 j=0 = (50,75)) — verified. The invariant is
q_window := denom(ord_y(Phi)/M) = 5a-3 exactly for F2 (gcd(2a-1,5a-3)=1).
(72,108)'s integral 12 is that corner's coincidence. Corrected in
G_SYSTEM_75_125.md. Consequence: generalize to RATIONALLY RELATED
GRADINGS (bigraded window lattice), not an a>=3 engine; tactical route
for (75,125) is a period-12 window compiler.
THE F2 FAMILY THEOREM (review-proposed, verified here symbolically in a):
with a=j+2, b=2a-1, t=5, kappa=3, c=y^2(y^3+1):
  f_a = -(1/(3a)) y^(2a-1) (y^3+1)^a solves the forcing ODE for EVERY a
  (collapse identity y*g'-3g = -3);
  N_a = 15a^2-13a+2;
  Phi_a = -(1/(3a)) y^(30a^2-24a+3) (y^3+1)^(15a^2-12a+2);
  block recurrence Phi_{a+1} = (a/(a+1)) C^(30a+3) Phi_a.
Reproduces both landed points ((50,75): N=36, (189,75,38,76); (75,125):
N=98, (504,201,101,202)). First FAMILY-level closed form in the program —
the forcing layer of all F2 in one formula. System growth is fixed
5-blocks (5 spares + 5 generators per a-step) — the certificate-tower
question (does the a=2 kill extend by a block rule?) is now well-posed.
Checker: f2_family_verify.py (14/14, symbolic-in-a) — wired into
run_tests.sh. Review's strategic program recorded: polygon-reduction
compiler (the missing front end), F2 j=0 end-to-end reproduction as
positive control, bigraded bridge, five orthogonal-axis next cases.


## 2026-07-24 — POLYGON-REDUCTION COMPILER (Opus lane): the missing front end built; (75,125) unconditional at the polygon layer
polygon_reduction.py: input a GGV chain case, output the COMPLETE
reduction — transform sequence (root shifts + fused Laurent inversion,
building on composite_charts' proven Jacobian -x^(l-2)), the DERIVED
bracket exponent kappa = l-2, both reduced Newton polygons, an explicit
BRANCH MANIFEST (every opposite-edge/endpoint/denominator-zero
alternative with its follow-or-exclude reason — no silent selection),
and the corner signature (t,kappa,a0,q,c).
Regression contract, all three targets PASS (56/56 checks):
R1: reproduces the PUBLISHED (72,108)/(8,28) reduction exactly (both
polygons match upstream_facts.json; branch manifest reproduces the
paper's own split: k=2 excluded, three-factor excluded, deep shift
excluded via GGV6 Prop 2.5).
R2: F2 j=0 = (50,75): t=5,kappa=3,a0=5,q=2,c=y^2(y^3+1); Phi
(189,75,38,76), N=36 — matches landed data and GGV3 sec 5.
R3: F2 j=1 = (75,125) derived WITHOUT a judgment flag — the chart is
IDENTICAL to j=0 (same corner, same inversion; only (m,n) scales the
polygon, not the chart), kappa forced by the fused-chart lemma; Phi
(504,201,101,202) matches. THE "UNREDUCED POLYGON" JUDGMENT (PHI_75_125
item 2) IS DISCHARGED — the (75,125) model is unconditional at the
polygon layer. Surviving judgment correctly re-scoped to the forcing
layer (corner-144 correspondence, audited only for (72,108)).
Branch manifests: (8,28) 4 branches/12 options; F2 j=0,j=1 each 3/6.
Checker: polygon_reduction_verify.py (56/56) — wired into run_tests.sh.
Files: polygon_reduction.py, polygon_reduction_verify.py,
POLYGON_REDUCTION.md.


## 2026-07-24 — CERTIFICATE TOWER (Opus lane): BLOCK-OBSTRUCTION — the new geometry at 125 is the incommensurate window lattice
First direct attack on (75,125). Verdict: the (50,75) a=2 kill does NOT
extend by the fixed 5-block rule.
THE a=2 CERTIFICATE (major artifact in itself): GGV3 section-5's (50,75)
kill reproduced BIT-EXACT by machine in both reduced charts (exact sympy
Groebner): gamma=2 — 13-equation terminal system eliminates to
<g-2^5, ..., g-5^4> so g-2^5 = g-5^4 = 0 verbatim as the paper states,
then the square-in-K((1/y)) step forces e-10 = 0 contradicting corner
primitivity (b6); gamma=3 — window elimination bottoms C0's y-support at
y^-6, c0,-10 = 0, contradicting primitivity (a6). Certificate kind: a
small set of TERMINAL COEFFICIENT EQUATIONS forcing a corner
window-coefficient to vanish — intrinsically BIGRADED (u-weight, y-order),
a window-depth contradiction, NOT a scalar covector/syzygy on the
G-system (the scalar G-ideal has the trivial solution and does not die).
THE TOWER STEP: the algebraic block layer transfers EXACTLY (generators
5->10, spares 4->9, Phi recurrence, literal nesting
coeff_{d0^2}(G1^{a=3}) = (10/3) G1^{a=2}). The KILL layer does not: it
lives in the y-order/window-cap layer governed by q_window = 5a-3, which
jumps 7 -> 12 with gcd(7,12) = 1 — INCOMMENSURATE lattices; the
forcing-slice fractional denominators fragment {1,7} -> {1,2,3,4,6,12};
the 5 new forcing generators + deepened Phi slice carry the new classes
with no period-7 analog.
NAMED OBSTRUCTING BLOCK: the window-cap/y-order layer of the 5 new
forcing generators. Exactly the review-predicted "new geometry at 125."
CONSEQUENCE: (75,125) not killed; the required tool is now UNIQUELY
DETERMINED — the period-12 window compiler / bigraded lattice engine
(already independently flagged by G_SYSTEM_75_125's correction). Three
independent threads (transfer test, family theorem, tower experiment)
now converge on the same missing tool.
Checker: f2_tower_verify.py (19/19) — wired into run_tests.sh.
Files: F2_TOWER.md, f2_tower.py, f2_tower_verify.py.


## 2026-07-24 — FAMILY GRAMMAR SWEEP (Opus lane): all 17 families classified, ZERO irregular — the finite grammar is real
The F2 closed form generalizes to a complete grammar governed by two
j-invariant corner numbers (gap, r), a 2x2 dichotomy:
THE ONE THEOREM (proved symbolically in j): the pure ansatz
f = A y^rho (y^dg+1)^e solves the forcing ODE for ALL j iff
t-(kappa+1)q+dg = 0 iff a0 = t(q-1) iff gap = 0, with UNIVERSAL constant
A = -1/(a*dg) (F2's -1/(3a) is the dg=3 instance; the collapse mechanism
is two j-constant identities + y g'-dg g = -dg).
Census: CLOSED-FORM 8 (PURE gap=0: F2, F9, F14; COFACTOR r=0 with unit
cofactor of degree gap and universal u(-1) = -1/a: F1, F5, F6, F8, F17);
RUNG-STRUCTURED 9 (F3, F4, F7, F10, F11, F12, F13, F15, F16 — the mu=dg
fully-ramified rung f = y^rho (y+1)^(dg*e-(dg-1)) u, deg u = gap+r, is a
UNIFORM closed form in j); IRREGULAR 0. Composite escapes: only F22
(gap = -1) outside the grammar.
Cross-checks 210/210: every landed point reproduced exactly (both F2
rungs, F9, F14, F1, all four PHI_F7 ramified polynomials substituted
directly into the ODE, F12 mu in {1,2,3} and F10 mu in {2,4} via the
mu-graded law). (72,108) honestly excepted (the (8,28) resonance corner —
no family formula claimed).
Surprises recorded: dg even iff RUNG on standard rows (the parity
mini-lemma falls out of dg = a0-q + the gauge); complex-scope discipline
maintained (real-emptiness = annotation only).
The grammar table (with per-family q_window laws) is direct design input
for the bigraded window engine — the tool all three threads converge on.
Checker: family_grammar_verify.py (210/210; independent full-ODE
substitution j=0..3 + symbolic-in-j structure) — wired into run_tests.sh.
Files: FAMILY_GRAMMAR.md, family_grammar.py, family_grammar_verify.py.


## 2026-07-24 — COVERAGE PROOF-DAG v1 (Opus lane): closure is now COMPUTED; three real inconsistencies found on day one
proof_dag.py builds the machine-enforced graph certificate -> state ->
cell -> branch -> subcase -> C0 (4455 nodes / 4445 edges; deterministic,
byte-identical regeneration). Evidence levels {claimed < exact-checked <
independently-audited < certified} propagate as minima; no hand roll-up
can assert closure.
Census: states 18 certified / 41 independently-audited / 334
exact-checked / 29 claimed; branches 2404 closed@claimed, 209 open
(26 sub2 + 171 sub1 + 12 defect-0) — matching the docs; C0 OPEN; 3
defect-0 family closures confirmed. Unmapped bucket: 43, surfaced loudly
(31 ledger-ambiguous, 8 unresolved certificates, 4 j6 kills without
joinable degrees). Audit priority queue (weakest edges): (1) one cascade
spec-audit join promotes 2401 branches from claimed; (2) 334 same-author
exact-checked states await independent audit; (3) the four
judgment-referenced subcase->C0 completeness edges cap C0.
THREE REAL INCONSISTENCIES (the DAG doing its job):
1. FRONTIER's "Killed (audited)" counts same-author exact-checked kills
   as audited (overstates 292 states) — labeling fix owed.
2. CURRENT_STATUS S1a asserts independent audit of cascade branch kills;
   the audit artifacts (audit_cascade_kills*.py) EXIST and run in the
   suite but are NOT machine-joined in DAG v1 — either the join (v2) or
   the claim wording must move; until joined, the DAG honestly rates
   those branches claimed.
3. ORPHAN CERTIFICATE: harvest:a8_dd2-inf_dd10_dsig5 is CERTIFICATE-FOUND
   but maps to no ledger state — the certificate corpus is stronger than
   the ledger records; mapping resolution owed.
proof_dag_report.py exits nonzero while inconsistencies stand — NOT wired
into run_tests.sh until the three fixes land (next audit round).
Files: proof_dag.py, proof_dag.json, proof_dag_report.py, PROOF_DAG.md.


## 2026-07-24 — CHAIN NATURAL-HISTORY SURVEY (Opus lane): GGV5 reproduced exactly (+ a paper typo found); the grammar picture refined
Full GGV5 pipeline re-implemented in pure Python; at M=35 reproduces the
published tables to the letter (24 chains exact; 23/24 (m,n) families
verbatim). THE EXCEPTION IS A PROVABLE TYPO IN GGV5: F6's printed base
pair (4,10) has gcd 2, violating the paper's own coprimality definition;
the algorithm returns the coprime correction (6j+7, 16j+18) satisfying
the same Diophantine identity — an erratum worth including in any GGHV
contact. Also: GGV5's prose "2 length-2 chains" contradicts its own
7-row table (algorithm sides with the table).
Extension to v11 <= 100: 3995 family rows, 3403 canonical admissible
chains (~170s full sweep).
DEMOGRAPHY: fine censuses NEVER flatten (A0-motifs 9->147, signatures
19->982, families 24->3995 — anti-signal at the labelled level), BUT
coarse regime clusters (gap-sign, r-sign, dg-parity, length) PLATEAU at
~20 across M=65-85 (26 when length-4 chains open at M=100). Verdict: a
BOUNDED SET OF REDUCTION SHAPES decorated by unbounded numeric labels —
the finite grammar is real at the shape level, exactly matching the
FAMILY_GRAMMAR 2x2 dichotomy picture. Max chain length creeps 2->3->4
(capped by bigOmega(gcd)+1). The window-denominator law (F2's 5a-3
analogue) is NON-TRIVIAL for all 3995 families — (72,108)'s integral
W_step confirmed as a coincidence at population scale.
CORRECTION APPLIED: composite_charts' "negative gap unobserved anywhere"
was wrong at scale — 588/3995 families populate gap<0 (no derived Phi
point yet); both files updated, verifier re-run green.
Checker: chain_survey_verify.py (18/18) — wired into run_tests.sh.
Files: CHAIN_SURVEY.md, chain_survey.py, chain_survey_verify.py,
chain_survey_data.json (3.3MB).


## 2026-07-24 — Eumemic study (Opus lane): deep method parallel found; variable map built; import candidates recorded
Read of eumemic/jacobian-keller-research (pinned commit 201c2f6): the
band method works in the localized Weyl algebra A_1[x^-1] = (+)_k x^k C[E];
[D,X]=1 decomposes into degree-free ladder rung-equations on root
multisets (necklaces) probed by trace-form covectors; the program forces
a counterexample pair to band <= 2 (proven tame, their tier: internally
audited) and stalls at the band-3 W2 hatch (a3=E(E+2)(E+4), b2=E(E+3)),
bounded-closed at d=3, arbitrary-d OPEN.
Map verdict: (a) slope laws ANALOGOUS not importable (their
shifted-power wall vs our kappa=t-2 — same role, different object);
(b) W2 vs our polygons DISJOINT (JC2=>DC1 gives no transfer either way);
(c) methodologically complementary, object-disjoint — one level up from
alok (who at least shared our polygons).
KEY FINDING — the same five-step schema reached independently by both
programs: graded commutator ladder -> per-rung elimination -> proven
degree-free slope law -> bounded exact UNIT-ideal kill -> open
arbitrary-degree residual needing a degree-free covector. Their W2 unit
ideal at d=3 mirrors our bridge UNIT pilot almost exactly. The tool BOTH
programs lack is the same tool.
Import candidates (eumemic_import_candidates.json, 5 + 3 non-candidates):
(1) their Nullstellensatz-multiplier certificate re-check format — hardens
our bridge kills; (2) their msolve ^-vs-** parser caveat (a wrong glyph
fabricates a FALSE SURVIVOR — directly relevant to our J6/msolve usage);
(3) trace-form covector calculus as the template for lifting bounded
per-state kills to degree-free closure (shared direction, their recipe
also OPEN); (4) joint glossary (Q_m ~ G-system, band ~ window, necklace ~
Newton corner, Q_0=1 ~ Phi != 0); (5) export offer: our proven kappa=t-2 +
arbitrary-k window caps are closed instances of slope laws their Gap 2
leaves open — offered as proof strategy.
Their audit discipline (self-flagged doc-ahead-of-proof rows) is itself
an import. Files: EUMEMIC_MAP.md, eumemic_import_candidates.json.
