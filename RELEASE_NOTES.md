# 2026-07-25 — the family wave

- **Certificate-tower experiment (F2_TOWER.md)**: first machine reproduction of
  GGV3 §5's historical (50,75) kill, bit-exact in both reduced charts. Verdict
  on extending it to (75,125) by the family block rule: **BLOCK-OBSTRUCTION** —
  the algebraic layer transfers exactly, but the kill lives in the window
  lattice with period q_window = 5a−3, and gcd(7,12)=1: incommensurate. The
  required tool (a bigraded/period-12 window engine) is now uniquely determined.
- **Family grammar (FAMILY_GRAMMAR.md)**: all 17 length-1 families classified —
  8 closed-form, 9 rung-structured, 0 irregular — governed by one theorem
  (pure ansatz ⟺ gap = 0, universal constant A = −1/(a·dg)). 210-check verifier.
- **Chain natural-history survey (CHAIN_SURVEY.md)**: GGV5's enumeration
  reproduced exactly (with a documented erratum: the printed F6 base pair
  violates the paper's own coprimality; corrected (6j+7,16j+18)), extended to
  v11 ≤ 100 (3995 families). Fine censuses never stabilize; coarse regime
  clusters plateau at ~20 — a bounded grammar of reduction shapes with
  unbounded numeric labels.
- **Polygon-reduction compiler (POLYGON_REDUCTION.md)**: input chain data,
  output the complete reduction with an explicit branch manifest. Reproduces
  the published (8,28) reduction exactly; derives (50,75) and (75,125) —
  discharging the standing "unreduced polygon" judgment for (75,125).
- **Coverage proof-DAG v1 (PROOF_DAG.md)**: closure is now a computed fact
  (certificate → state → cell → branch → subcase → C0; 4455 nodes). Its
  report currently finds 3 real doc-vs-data inconsistencies (fixes in
  progress) — shipped as-is: this is the machine-honesty layer working.
- **Proof-gate hardening (HARDENING_NOTES.md)**: 532 proof-critical asserts
  across 29 checkers converted to explicit exit-1 gates (immune to python -O);
  generated artifacts now byte-deterministic (SHA-256 + git-commit provenance).
- Cross-program/literature: EUMEMIC_MAP.md (shared five-step schema with the
  Weyl-algebra program; import candidates), ML_RESTRICTION.md (Makar-Limanov
  2025 restriction: inapplicable to these polygons — machine-checked),
  ALT_ELIM.md (alt-bridge spare elimination: sound, wall persists —
  formulation-level), kill_certificate_msolve.py (mod-p reconstruction lift
  route for the hard certificates; batch execution pending).
- chain_survey_data.json is 3.3 MB (full enumeration; JSON-artifact precedent).

# Release notes — public release tree

## What this tree is

This is the clean, public-release version of the **(72,108) plane Jacobian
program** plus the **dimension-3 Jacobian counterexample** re-verification. It is
a reorganized, redacted subset of a larger private working repository, prepared
for public review by mathematicians and for auditing AI-assisted mathematics.

Layout:

| dir / file | contents |
|---|---|
| `d2/` | the (72,108) plane program: checkers, engines, JSON artifacts, reports, worklists (internal structure preserved so sibling-file relative paths still resolve) |
| `d3/` | the dimension-3 counterexample and its exact verifier |
| `lean/` | Lean 4 feasibility certificate (`lake build` clean) |
| `docs/` | a pointer index mapping a reader's path into `d2/` (canonical copies stay in `d2/`) |
| `README.md`, `VERIFICATION.md` | the claim, the headline theorem, the frontier, and the 15-minute exact-arithmetic path |
| `run_tests.sh`, `setup.sh`, `requirements.txt` | the suite and environment setup |
| `CITATION.cff`, `LICENSE` | citation metadata; MIT for code, CC-BY-4.0 for the math documents, third-party paper sources not included |

## Update — 2026-07-24 (post-v0.1.0; no tag yet, published for review)

A trust-layer + transfer-test drop on top of the 2026-07-23 corner-law layer.
New content (all git-tracked; the exact files map into `d2/`):

1. **Certificate architecture — first end-to-end pass.** Per-kill cofactor
   certificates land under `d2/kill_certificates/` (49 records, ~19 MB;
   largest single JSON ~4.7 MB — shipped, they are the trust artifact), with
   `kill_manifest.json`, the producer `kill_certificate_tools.py`, and the
   engine-free consumer `audit_gb_kills.py`. Census: **20 CERTIFIED, 29
   not-yet-certificated**, all failures logged honestly in
   `kill_certificates/status_log.json`.
2. **Independent 49-kill audit (0 disagreements).** `audit_alt_hunt_kills.py`
   + `audit_alt_hunt_census.json`: a from-scratch spec-only re-derivation
   (producer code neither imported nor read) of all 49 forced HUNT/J6 state
   kills — **41 FULLY-VERIFIED + 8 VERIFIED-DATA-ONLY, 0 DISAGREEMENT**.
3. **Independent CAS replay of the f37 theorem.** `F37_REPLAY.md` +
   `f37_replay_m2.m2` (Macaulay2) + `f37_replay_sage.py` (Sage) +
   `f37_replay_selftest.py` (pure-Python construction self-test, 8/8).
4. **μ-ladder + parity theorem at dg = 4.** `MU_RUNGS_F10.md` +
   `mu_rungs_f10*.py` prove the even-`dg` parity claim at dg=4 (F10);
   `REVIEW_ZETA_MU.md` + `review_zeta_mu.py` are an adversarial re-derivation
   (31/31 confirmations). Corrections folded into `ZETA_TAIL.md`,
   `COMPOSITE_CHARTS.md`, `PHI_75_125.md`, `CASE_COMPILER.md`.
5. **Transfer test — phases 1–2 on (75,125).** `C_SERIES_75_125.md`
   (+ verifier) derives the tower length **N = 98**; `G_SYSTEM_75_125.md`
   (+ `g_system_75_125.py/.json/_verify.py`) builds the G-system and locates
   its window-cap obstruction at **a ≥ 3**.
6. **Cross-program corroboration.** `ALOK_CROSSCHECK.md` + `alok_crosscheck.py`:
   exact setup corroboration against an independent parallel program, regime
   disjointness quantified, **0 findings**.
7. **Alt-bridge construction + honest wall (PENDING AUDIT).** `ALT_BRIDGE.md`
   + `alt_bridge.py`, `J6_MSOLVE.md` + `j6_msolve.py`, `R9_SYMBOLIC.md` /
   `R9_VALSPLIT.md` + `r9_symbolic_elim.py` / `r9_valsplit.py`: attempted
   state-level bridge kills with their honest negative outcome (Gröbner cost
   wall survives). Their PENDING-AUDIT labels are kept intact and these are
   **not** counted in the frontier accounting (`d2/CURRENT_STATUS.md` §3c).

Four new exact verifiers wired into `run_tests.sh`
(`mu_rungs_f10_verify.py`, `c_series_75_125_verify.py`, `alok_crosscheck.py`,
`g_system_75_125_verify.py`); the full suite is green from this tree, and
`tools/clean_clone_check.py` confirms every file the suite reads is tracked.
The certificate/audit tools and the CAS-replay scripts that need Macaulay2 /
Sage / msolve / Singular are shipped as artifacts but not wired into the
pure-Python suite (matching the source repo's own suite).

## Update — 2026-07-23 (post-v0.1.0; no tag yet, published for review)

Thirty-seven files added/updated from the source repository. New content:

1. **The corner-law structure layer** — closed-form Φ derivations and the
   μ-graded signature law at twelve exact points (`d2/PHI_75_125.md`,
   `d2/PHI_CORNER4.md`, `d2/PHI_F14.md`, `d2/PHI_F7.md`,
   `d2/COMPOSITE_CHARTS.md`, `d2/ZETA_TAIL.md` + derivation `.py` files),
   each with its own exact verifier in `run_tests.sh`, plus an independent
   skeptical review pass of the ζ/μ layer (`d2/REVIEW_ZETA_MU.md`).
2. **Prior-art audit** (`d2/PRIOR_ART.md`) + the zero-freedom postdiction of
   GGHV22 §4's printed `f₁` (`d2/prior_art_postdiction_verify.py`).
3. **Window caps k=6,7,8 proven** (`d2/WINDOW_CAPS.md`); the corresponding
   judgment flags in `d2/FULL_SYSTEM_BRIDGE.md`/`d2/BRIDGE_SWEEP.md` (both
   newly shipped, with `d2/full_system_bridge.py`) are retired.
4. **Galois-descent library census** (`d2/GALOIS_LIBRARY.md`) and the pilot
   **case compiler** (`d2/CASE_COMPILER.md` + three dossier JSONs).
5. **Kill-side status documents** — `d2/ALT_HUNT.md`, `d2/J6_MSOLVE.md`,
   `d2/R9_SYMBOLIC.md`: the s-unit BM-candidate residual layer fully killed
   at engine level (49/49 states) and the dm4-elimination negative result —
   **all PENDING AUDIT** and not counted in the frontier accounting
   (`d2/CURRENT_STATUS.md` §3b).

Ten new verifiers appended to `run_tests.sh`; the full suite is green from
this tree. Machine kill-records and certificate JSONs for the pending-audit
layer stay in the source repository until the certificate audit round lands.

## Provenance

Generated from the source working repository. The 2026-07-23 content was cut at
source commit `05e6609650322e1a646861584011e398ec3db338`; the 2026-07-24 trust
layer + transfer test at source commit
`8ba76adcd9f3f977ad7b763c3f37e4b74eaec501`. Only git-committed source revisions
were shipped.

Two files required by the 2026-07-23 `run_tests.sh` that were untracked in the
source working tree at that commit are included here (they are load-bearing
spec-only auditors): `d2/audit_convolution_kills_r2.py` and
`d2/audit_reconstruction_kills.py`.

## Exclusions and redactions applied

Per the source repo's internal publication audit:

1. **Third-party paper sources removed.** The four arXiv LaTeX manuscripts under
   `d2/paper_src/*.tex` (GGV1 1401.1784, GGV3 1406.0886, GGV5 1708.07936,
   GGHV22 2204.14178) are other authors' copyrighted sources and are replaced by
   a links-only `d2/paper_src/README.md`. `paper_src/next_cases.py` (original
   code) is kept.
2. **Path leaks redacted.** The one personal absolute path in a shipped file
   (`d2/F37_FRONTIER.md`) was replaced with a repo-relative path. The other leak
   lived in `q1_msolve.log`, which is dropped (see #3).
3. **Computation logs dropped.** All `*.log` run traces (33 tracked + a few
   untracked) removed; none is load-bearing. `*.log` added to `.gitignore`.
4. **Numeric pickles dropped; verification path de-pickled.** The nine
   regenerable numeric-search `*.pkl` blobs are removed. The verification path no
   longer unpickles anything: the pre-resultant generators are shipped as exact
   term lists in `d2/generators.json` (emitted once from `t4_state.pkl`), and
   `f37_sat_verify.py`, `f37_free_family_verify.py`, and the Lean exporter parse
   that JSON. `d2/t4_state.pkl` is KEPT only for optional provenance (the
   checkers confirm `generators.json` reproduces it when present); see
   `d2/T4_STATE_PROVENANCE.txt` (pickle-trust: `regenerate_system.py` rebuilds it
   from scratch). Nothing mandatory unpickles.
5. **Batch raw pass intermediates dropped, final JSONs kept.** Removed
   `batch_convolution_sub1_gauge_raw.json`,
   `batch_convolution_sub2_gauge_raw.json`,
   `batch_convolution_sub2_gauge_resume.json`,
   `batch_convolution_sub2_pass1_ungauged.json`,
   `batch_convolution_sub2_round2_part2.json`. The final batch JSONs read by the
   checkers are kept.
6. **Overnight/scratch/checkpoint files dropped.** `overnight_batch.py`,
   `overnight_r9.py`, `launch_overnight.py`, and untracked scratch
   (`triage_harvest.*`, `overnight_batch_final.json`) removed; none is imported
   by any checker.
7. **Internal/private docs excluded.** `CONTACT_DRAFT.md` (private),
   `PUBLICATION_AUDIT.md` (internal), the git-ignored `math_stuff_field_audit/`
   staging copy (its two substantive files already ship tracked as
   `d2/FIELD_SPLIT_AUDIT.md` and `d2/t5_split_place_verify.py`), and `.claude/`
   tool config are all excluded.
8. **`__pycache__/` bytecode** stripped.

No credentials, API keys, or tokens were present in the source (confirmed by the
source audit's regex sweep).

## Verification gates at release

- `SKIP_SLOW=1 bash run_tests.sh` — D3 exact verification + all D2 exact proof
  checkers pass in this tree. (`SKIP_SLOW=1` skips only the multi-minute numeric
  positive control `jetlift.py control f31_sub2`; the default run adds it.)
- `lake build` — green in `lean/`.
