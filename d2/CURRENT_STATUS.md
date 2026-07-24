# CURRENT STATUS — the (72,108) program (current truth)

**As of:** 2026-07-24. This file states only what is *currently true*: proven and
checker-enforced nodes, trusted published inputs, independently audited results,
and the live open frontier. Engine kills that are **pending independent audit are
not counted as progress** below. For the full historical claim graph (including
superseded, conditional, and evidence-only nodes) see `PROOF_INVENTORY.md`, which
is a historical record — this file supersedes it for current state.

---

## 0. Target

**C0 (OPEN).** There is no `P,Q ∈ K[x,y]` with `[P,Q] = x²` and the Prop-4.3
case-(8,28) Newton polygons (subcases (1) and (2)). If closed, the plane-Jacobian
counterexample degree bound rises 108 → 125. **Not proven.** The remaining
obligation is the `f31` window/cascade frontier in §3.

---

## 1. Proven / checker-enforced nodes (current)

Scope note: all algebraic membership/closure statements below are
**characteristic zero** (over every `ℚ`-algebra) unless stated otherwise. The one
place this bites is the f37 closure (C11): its integer certificate carries a
denominator-clearing multiplier `D = 46875 = 3·5⁶`, so the conclusion is asserted
over fields of characteristic `≠ 3, 5` — see the field-scope note in
`F37_SATURATION_REPORT.md`.

### 1a. Independently audited (tier 1 — separately authored, no code sharing)

| Node | Statement | Checker(s) in suite |
|---|---|---|
| C14–C16 | Field-split framework: `q` splits into four geometric places after base change; case space is multiplicity vectors up to `S₄`; the `σ`-locus T3 branch is empty in char 0; the geometrically q-coprime `a_t=7` stratum is dead. | `t5_split_place_verify.py`, `test_split_place_proofs.py`, `test_split_place_ledger.py` |
| C17 | Sub2 terminal ledger (327 geometric vectors; terminal pruning), terminal layer confirmed 654/654. | `test_split_place_ledger.py`, `audit_cascade_kills.py` |
| C18 | Sub2 depth-4 cascade: 420 open → 390 killed, **30 survive**; engine-proven AND spec-only audited (654/654, 420/420, 390/390). | `test_cascade_engine.py`, `audit_cascade_kills.py` |
| C20 | Sub2 t-place coupling audited (7872/7872, 0 disagreements). | `audit_tplace_cases.py` |
| C27–C29 | Sub1 cascade caps (independent correction of the `15+3a` cap), terminal ledger (2614/2614), depth-4 cascade 2178→1899 killed, **279 survive**. | `sub1_cascade_verify.py`, `audit_cascade_kills_sub1.py` |
| C31 | Sub1 t-place coupling audited (41592/41592). | `audit_tplace_cases.py` |
| C43 | Infinity layer, both windows, spec-only audited (sub2 420/420, sub1 2178/2178, zero disagreements). | `test_cascade_inf.py`, `audit_inf_cases.py` |
| C44 | Alternate-regime degree layer audited (25 kills re-derived, 27 OPEN with full state partitions). | `alt_regime_inf_verify.py`, `alt_inf_sweep_verify.py`, `audit_alt_regime.py` |

### 1b. Exact same-author checkers in the suite (tier 2)

| Node | Statement | Checker |
|---|---|---|
| **C11** | **f37-branch closure:** `f31 ∈ ⟨G1,G2,G3,G5body+Φ⟩` over `ℚ` (Φ free); the elimination ideal in `(d̃2,d̃1,d̃0,d₋₁,Φ)` is exactly `⟨f31⟩`; f37 and `d₋₁²¹` are classical resultant excess. The whole `{f37=0}` branch is an artifact, both subcases (char 0). | `f37_sat_verify.py` (sympy re-expansion of the `lift()` certificate; parses `generators.json`, no pickle) |
| C12 | f37 free family `d̃2=d̃1=0` does not lift (subsumed by C11). | `f37_free_family_verify.py` |
| C6 | Envelope/window bounds `deg d_{4−k} ≤ 14w`(sub2)/`15w`(sub1), `ord ≥ 12w`; Φ attains both (tight). Now checker-enforced. | `envelope_bounds_verify.py` (reads `paper_src/upstream_facts.json`) |
| C21–C25 | Sub2 per-cell kills: `a_t=9` T2 dead; `a_t=9` T1 constant-E dead; nonconstant reduction; T2-column squeeze (4 of 12 cells); T2 infinity narrowing. | `t5_90t2_verify.py`, `t5_90t1_verify.py`, `t5_90t1_constant_verify.py`, `t5_90t1_local_verify.py`, `t5_t2_column_verify.py`, `t5_t2_infinity_verify.py` |
| C33–C34 | Alternate-regime flipped cascade levels (→ 27 residual branches). | `alt_regime_verify.py`, `alt_regime_l2_verify.py` |
| C41–C42 | Max-plus infinity layer + tie equations. | `test_cascade_inf.py`, `cascade_inf_ties_verify.py` |
| C46 | No-jet-kill theorem: all ten sub1 t-place tied supports have smooth rational points, so t-place obligations are constraints (never local kills) for `a ∈ [2,10]`. | `residue_lemmas_depth_verify.py` |
| — | Residue-lemma libraries (shared, feeding the kills). | `residue_lemmas_verify.py`, `alt_residue_congruences_verify.py` |

The full suite is `../run_tests.sh` (public: `run_tests.sh`); its first line is
`tools/clean_clone_check.py`, which fails loudly if any suite checker reads a file
not tracked in git.

### 1c. Trusted published inputs (tier 4 — used as stated, not re-proven)

- **GGHV22** (`arXiv:2204.14178`): Prop 4.3 case-(8,28) Newton polygons; the
  normalization and valuation-induction templates (Props 5.2/5.6). The finite
  facts the suite consumes are transcribed with line citations in
  `paper_src/upstream_facts.json`.
- **GGV1** (`arXiv:1401.1784`): Propositions 1.13 and 2.1 (common-root
  normalization `ℓ(P)=R²`, `ℓ(Q)=R³`).
- **GGV5** (`arXiv:1708.07936`): next-case / corner-144 context (not load-bearing).
- **Mason–Stothers** (abc for polynomials): used in the split-place margins.

---

## 2. Formal (Lean) certificate — current scope

`lean_certificates/` (public: `lean/`) **kernel-checks the emitted integer
certificate** `D·f31 = Σ (D·cᵢ)·Gᵢ` for node C11 **under a small custom
sparse-polynomial implementation** (`Cert/Poly.lean`); `D = 46875 = 3·5⁶`. It is
NOT an abstract-algebra formalization of ideal membership, and it is not part of
`run_tests.sh`. Trust base and the honest `add`/`mul`-faithfulness caveat are in
`LEAN_CERTIFICATE.md`.

---

## 3. Open fronts (the live frontier)

Machine-generated by `frontier_gen.py` (see `FRONTIER.md`). **224 surviving
cells/branches** remain open, with a degree-state overlay:

| Front | Surviving branches/cells | Flag cases | Degree-state overlay |
|---|---|---|---|
| Sub2 cells (q+t+inf, depth 4) | **26** (a ∈ [5,10]) | 220 | 7888 residual states (`phase_d_states_sub2.json`) |
| Sub1 standard regime (q+t+inf, depth 4) | **171** (a ∈ [2,10]) | 1145 | 44117 residual states (`phase_d_states_sub1.json`) |
| Sub1 alternate regime (a ∈ [11,15], v<0) | **27** | — (tracked as states) | 38360 states → 4690 surviving, 33670 killed |
| **Total** | **224** | 1365 | — |

State-kill overlay: the infinity/degree sweeps remove degree-states within each
surviving cell (e.g. alternate regime 33670/38360 states killed). These reduce
each cell's residual obligation but **do not yet close any of the 224 cells**; the
per-cell closures are the remaining work.

**Pending-audit note.** The engine sweeps that produced the 224 frontier also
report large "killed (pending audit)" branch counts; those kills that have NOT
been re-derived by a spec-only auditor are **excluded** from the proven/audited
tables in §1 and are not counted as closed. The independently audited layers
(C18/C20/C29/C31/C43/C44) are the ones reflected in §1a.

---

## 3b. Structure layer — the corner law (added 2026-07-23; tier 2 unless noted)

A new cross-family layer, developed on the GGV/GGHV corner data and checked by
ten new suite verifiers (all exact sympy, wired into `run_tests.sh`):

- **Corner law for the tower's last element Φ** (`PHI_75_125.md`,
  `PHI_CORNER4.md`, `PHI_F14.md`, `PHI_F7.md`): twelve exact
  derived/audited points across four regimes (t ∈ {3,4,5,7}, both gap
  regimes), unified as a **μ-graded law** (`ZETA_TAIL.md`); the ramified and
  unramified branch laws are its μ = dg and μ = 1 specializations.
- **κ = t−2 chart theorem** (`PHI_CORNER4.md`, extended in
  `COMPOSITE_CHARTS.md`): proven for the standard chart class and every
  escape family (composite charts fuse to a single inversion); the true
  boundary for A0′=(2,0) families is the ζ-defect model correction
  (`ZETA_TAIL.md`), with surviving models enumerated.
- **Prior-art audit** (`PRIOR_ART.md`): the corner-law claims checked against
  the published GGV corpus (read, not searched). Novelty verdicts recorded
  per claim, with required citations. Flagship check:
  `prior_art_postdiction_verify.py` re-derives GGHV22 §4's printed
  closed-form `f₁` and confirms the corner law **postdicts it with zero
  fitting freedom** (12 checks).
- **Galois-descent library** (`GALOIS_LIBRARY.md`): the 23-shape residue
  library classified structurally; C08 and C20 are the exactly-two
  quadratic-obstruction kills, with a two-line transfer criterion for any
  family member.
- **Case compiler** (`CASE_COMPILER.md`, three dossiers): corner data →
  compiled case dossier, validated to reproduce the audited (72,108) facts
  exactly; conjectural flags are data-driven per regime.
- **Window caps k=6,7,8** (`WINDOW_CAPS.md`): the formerly flagged extension
  is now recited and proven (81 checks); the corresponding [judgment] flags
  in `FULL_SYSTEM_BRIDGE.md`/`BRIDGE_SWEEP.md` are retired.

**Trust architecture (first end-to-end pass — see §3c for the 2026-07-24
state).** A cofactor-certificate pipeline (`kill_certificate_tools.py` + the
certificate-consuming spec auditor `audit_gb_kills.py`) extracts
engine-free-checkable certificates `1 = Σ cᵢfᵢ` for saturated-Gröbner kills.
As of 2026-07-24 the first full pass has landed: **20 of the 49 targeted kills
carry a certificate (`CERTIFICATE-FOUND`); 29 remain `NOT-YET-CERTIFICATED`**
(lift timeouts / open failures, all recorded honestly in
`kill_certificates/status_log.json`). Separately, the s-unit BM-candidate
residual layer is fully killed at engine level (`ALT_HUNT.md`, `J6_MSOLVE.md`:
49/49 states, mod-p corroboration) and has now been re-derived by an
independent spec-only auditor (`audit_alt_hunt_kills.py`, census
`audit_alt_hunt_census.json`: **49/49, zero disagreements**). These layers
are still labelled **PENDING AUDIT** in their own documents and per the rule
above are NOT counted in §1 or in the frontier table yet.

---

## 3c. Update 2026-07-24 (trust layer + transfer test; tier as noted)

Same-author checker-enforced unless stated; published here for external review.

- **Certificate architecture — first end-to-end pass.** The kill-certificate
  pipeline now produces per-kill JSON certificates under `kill_certificates/`
  (49 records), consumed by two independent readers: `audit_gb_kills.py`
  (re-expands `1 = Σ cᵢfᵢ` from the certificate, no engine trust) and
  `kill_certificate_tools.py` (producer/manifest, `kill_manifest.json`).
  Census: **20 CERTIFIED, 29 not-yet-certificated**, failures logged, nothing
  overclaimed.
- **Independent 49-kill audit (0 disagreements).** `audit_alt_hunt_kills.py`
  — a from-scratch, spec-only re-derivation with no access to the producer
  code (`alt_hunt_depth2.py`, `j6_msolve.py` are neither imported nor read) —
  re-checks all 49 forced HUNT/J6 state kills via isolated per-kill SymPy
  Gröbner subprocesses. Result: **41 FULLY-VERIFIED + 8 VERIFIED-DATA-ONLY,
  0 DISAGREEMENT, 0 UNPARSEABLE** (`audit_alt_hunt_census.json`, exit 0).
- **f37-branch closure — independent CAS replay scripts.** `F37_REPLAY.md`
  ships standalone Macaulay2 (`f37_replay_m2.m2`) and Sage
  (`f37_replay_sage.py`) re-derivations of the C11 f37-artifact theorem, plus
  a pure-Python construction self-test (`f37_replay_selftest.py`, 8/8) that
  needs no external CAS.
- **μ-ladder + parity theorem (dg = 2, 4).** `MU_RUNGS_F10.md` +
  `mu_rungs_f10_verify.py` extend the μ-graded law (`ZETA_TAIL.md`) up the
  ramified ladder and **prove the even-`dg` parity claim at `dg = 4` (F10)**:
  the μ=1 rung is empty, matching dg=2. `REVIEW_ZETA_MU.md` +
  `review_zeta_mu.py` are an adversarial re-derivation of the ζ/μ layer
  (**31/31 confirmations hold**, F12 quartic re-obtained by a different
  elimination route).
- **Transfer test — phases 1–2 on (75,125).** The corner-law machinery is
  exercised on a *second* family. Phase 1: the (75,125) C-series is built
  from the polygon data and its tower length **N = 98 is DERIVED**, matching
  the formula (`C_SERIES_75_125.md` + `c_series_75_125_verify.py`). Phase 2:
  the associated G-system is built and its weight-integrality boundary
  characterised — the quasipolynomial window-cap obstruction is located at
  **a ≥ 3** (`G_SYSTEM_75_125.md` + `g_system_75_125_verify.py`).
- **Cross-program corroboration (Alok).** `ALOK_CROSSCHECK.md` +
  `alok_crosscheck.py` reconcile the (72,108) setup and the cascade/infinity
  state artifacts against an independent parallel program: exact setup
  corroboration, regime-disjointness quantified, **0 findings**.
- **Alternate-bridge construction + honest wall.** `ALT_BRIDGE.md` +
  `alt_bridge.py` and the R9 symbolic/valuation-split status docs
  (`R9_SYMBOLIC.md`, `R9_VALSPLIT.md`, `r9_symbolic_elim.py`,
  `r9_valsplit.py`) record the attempted state-level bridge kills and their
  **honest negative outcome** (Gröbner cost wall survives; every completed
  verdict is COST, no survival signal) — **PENDING AUDIT**, not counted in §1.

---

## 4. One-line summary

The f37 branch is closed (C11, char 0; Lean-kernel-checked certificate). The
`f31` branch is an open, explicitly enumerated **224-cell frontier** with finite
per-cell residue/degree obligations, built entirely from independently audited
kills. C0 remains open.
