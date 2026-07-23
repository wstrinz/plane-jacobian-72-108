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

## Provenance

Generated from the source working repository at commit
`05e6609650322e1a646861584011e398ec3db338`
(branch `claude/d2-jacobian-counterexamples-sgp9in`).

Two files required by `run_tests.sh` that were untracked in the source working
tree at that commit are included here (they are load-bearing spec-only auditors):
`d2/audit_convolution_kills_r2.py` and `d2/audit_reconstruction_kills.py`.

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
