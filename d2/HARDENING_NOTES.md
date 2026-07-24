# HARDENING_NOTES.md

Mechanical safety hardening from an external review. Two issues addressed. No
`git commit` performed. Concurrent-lane files (`alt_elim*`,
`kill_certificate_msolve*`, `proof_dag*`, `ml_restriction*`/`ML_RESTRICTION*`,
`f37_replay_m2.m2`) were **not** touched.

---

## ISSUE 1 — proof-critical `assert` statements stripped by `python -O`

Under `python -O`, `assert` statements are removed, so a proof checker could
report success (exit 0) without ever evaluating its decisive mathematical
condition. Every gate checker driven by `../run_tests.sh` (D2 section) whose
exit code is a proof gate was swept.

### Method

Each proof-critical `assert COND[, MSG]` was rewritten to an explicit call

```python
_require(COND, MSG-or-source-repr)
```

where a small helper injected just below each file's imports does the gating:

```python
def _require(_cond, _msg):
    """Proof-critical check: fails loudly and exits nonzero, unaffected by python -O."""
    if not _cond:
        import sys as _sys
        print("FAIL: " + str(_msg))
        _sys.exit(1)
```

This matches each file's existing fail style (they print a `PASS` summary on
success; `_require` prints `FAIL: ...` and exits nonzero on failure). The
transform was AST-based (byte/char-offset-precise), so it correctly handles
inline asserts (`B=sum(b); assert B<=15-a<=4`), multi-line asserts, files with a
UTF-8 BOM, and files with non-ASCII comments. All `print(...)` check-count
reporting lines were left byte-for-byte unchanged. Every transformed file was
re-parsed to confirm **zero** remaining `assert` statements and that the
`_require` count equals the original `assert` count.

Note: the `audit_*.py` gate checkers already used `raise AssertionError(...)`
(not `assert`), which `-O` does **not** strip — they were already safe and were
left unchanged. Files listed in `run_tests.sh` that contain no `assert`
(`t6_premises_verify.py`, `alt_combined_verify.py`, `phase_f2_pilot_verify.py`,
`envelope_bounds_verify.py`, `phi_75_125_verify.py`, `phi_corner4_verify.py`,
`phi_f7_verify.py`, `prior_art_postdiction_verify.py`, `composite_charts_verify.py`,
`zeta_tail_verify.py`, `mu_rungs_f10_verify.py`, `alok_crosscheck.py`,
`g_system_75_125_verify.py`, `case_compiler_verify.py`, `window_caps_verify.py`,
`convolution_elim.py`) needed no change. (`g_system_75_125_verify.py` has the
word "assert" only in prose comments — 0 real assert statements.)

### Files changed and assert statements replaced (ISSUE 1)

| File | asserts replaced |
|---|---:|
| t5_split_place_verify.py | 21 |
| test_split_place_proofs.py | 16 |
| test_split_place_ledger.py | 8 |
| test_cascade_signature.py | 6 |
| test_cascade_engine.py | 16 |
| test_cascade_inf.py | 25 |
| cascade_inf_ties_verify.py | 21 |
| test_cone_lemmas.py | 13 |
| residue_lemmas_depth_verify.py | 19 |
| t5_90t2_verify.py | 13 |
| t5_90t1_verify.py | 29 |
| t5_90t1_constant_verify.py | 2 |
| convolution_descent.py | 3 |
| t5_90t1_local_verify.py | 6 |
| f37_free_family_verify.py | 12 |
| f37_sat_verify.py | 3 |
| sub1_cascade_verify.py | 49 |
| alt_regime_verify.py | 31 |
| t5_t2_column_verify.py | 42 |
| t5_t2_infinity_verify.py | 42 |
| alt_regime_l2_verify.py | 20 |
| residue_lemmas_verify.py | 30 |
| alt_regime_audit_verify.py | 32 |
| alt_regime_inf_verify.py | 28 |
| alt_inf_sweep_verify.py | 2 |
| alt_residue_congruences_verify.py | 35 |
| phi_f14_verify.py | 1 |
| galois_library_verify.py | 6 |
| c_series_75_125_verify.py | 1 |
| **29 files** | **532** |

Plus one proof-critical assert converted in a generator while editing it for
ISSUE 2:

| File | assert replaced |
|---|---:|
| frontier_rollup.py | 1 (the "DRIFT" per-cell census cross-check → explicit `sys.stderr.write` + `sys.exit(1)`) |

All conversions replace the decisive condition; none were left as dev-only
invariants ("explicit is never wrong"). The transformed asserts are a mix of
certificate/identity residual checks (`sp.expand(...) == 0`, `sp.cancel(...) == 1`),
census/count checks (`len(...) == N`, enumerated-list equalities), and
degree/valuation bound checks — all mathematical claims.

---

## ISSUE 2 — nondeterministic generated artifacts

Generators embedded wall-clock times and source **file mtimes**, breaking
byte-identical regeneration from a clean clone. Fixed in the two named
generators plus `frontier_gen.py` (which embedded per-artifact mtimes — the
exact antipattern).

Replaced each embedded timestamp/mtime with a deterministic provenance block:
generating **git commit** (`git rev-parse HEAD` via subprocess, falling back to
`$GIT_COMMIT` or `"unknown"`), **SHA-256** over every source artifact actually
read (content-addressed, mtime-independent), and **schema version**.
`SOURCE_DATE_EPOCH` is honored: when set, an optional reproducible UTC build
date is emitted (formatted from the epoch); when unset, **no** wall-clock value
is embedded at all.

| File | output | change |
|---|---|---|
| state_kill_ledger.py | state_kill_ledger.json | `"generated": datetime.now().isoformat()` → deterministic `"generated"` + `"provenance"` block (git commit, `source_sha256`, `source_files`, `source_date_epoch`, `generated_utc`). `load()` now records every file read. |
| frontier_rollup.py | FRONTIER_V2.md | `now = datetime.now().strftime(...)` and the "Generated: … (local clock)" line → a "Provenance: git … \| sources sha256 … \| schema 1" line (plus "Built <date>" only when SOURCE_DATE_EPOCH is set). Also the DRIFT `assert` → explicit exit (see ISSUE 1). |
| frontier_gen.py | FRONTIER.md | `os.path.getmtime` per artifact and `datetime.now()` → per-artifact SHA-256 (first 16 hex) table + a "Provenance: git …, schema 1" line honoring SOURCE_DATE_EPOCH. |

### Determinism result (regenerate twice, diff — actual bytes)

Each generator was run twice consecutively and the two outputs compared with
`cmp`:

- `state_kill_ledger.json` — **BYTE-IDENTICAL**
- `FRONTIER_V2.md` — **BYTE-IDENTICAL**
- `FRONTIER.md` — **BYTE-IDENTICAL**

With `SOURCE_DATE_EPOCH=1700000000` set, the emitted date is the deterministic
`2023-11-14T22:13:20Z` (UTC) and two runs at that epoch are again
**BYTE-IDENTICAL**. Outputs were regenerated once more with the variable unset
to leave the repository in the non-dated deterministic form.

Scope note: many heavy-compute derivation scripts (`batch_convolution_*.py`,
`phase_f2_*.py`, `d2_threshold.py`, `msolve_bridge.py`, `overnight_*.py`,
`nulla_unit_cert.py`, `jetlift.py`, etc.) embed `elapsed_seconds`/`ts` timing in
their result JSONs. These are compute logs, not documentation artifacts; they
are **not** invoked by the test suite or the frontier regeneration story, and
several are owned by concurrent lanes — left unchanged and out of scope.

---

## Verification

- **All 29 edited checkers (ISSUE 1) plus the 3 generators (ISSUE 2) run at their
  intended exit code.** Each edited checker was executed once with its
  `run_tests.sh` arguments; every one exited **0** (see the exit-code table in
  the session report).
- **`-O` gate demonstration** on three representative edited checkers, using
  scratch copies under the system temp dir with a corrupted input (never the
  repo):
  - `test_split_place_proofs.py` (corrupted `f31_graded.txt`, `8192`→`8193`):
    edited under `python -O` → **rc=1** (gate fires). Control: the original
    bare-`assert` version under `python -O` → **rc=0** (FALSE PASS — the exact
    vulnerability), and under plain `python` → rc=1.
  - `test_cascade_engine.py` (corrupted `split_place_ledger.json`, one flipped
    terminal-feasibility verdict): edited under `python -O` → **rc=1**.
  - `alt_regime_verify.py` (corrupted `f31_graded.txt`): edited under
    `python -O` → **rc=1**.

## Files modified (summary)

29 gate checkers (ISSUE 1) + `frontier_rollup.py`, `state_kill_ledger.py`,
`frontier_gen.py` (ISSUE 2; `frontier_rollup.py` also had its one assert
converted). Regenerated artifacts left in deterministic form:
`state_kill_ledger.json`, `FRONTIER_V2.md`, `FRONTIER.md`.
