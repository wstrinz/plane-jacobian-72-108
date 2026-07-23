# Verification quickstart — the 15-minute path

For a skeptical mathematician who wants to confirm the load-bearing claims
themselves, fast, in exact arithmetic. Everything below is exact (SymPy /
Singular), exits nonzero on any failure, and needs no cluster.

## Setup (~1 min)

```bash
pip install sympy          # the two headline checks (a, b) need only sympy
```

Python 3.10+. The full suite (c) also uses numpy/scipy/mpmath — `pip install -r
requirements.txt` — and one optional step wants Singular; neither is needed for
the two checks that matter most.

## The ladder

Climb as far as your skepticism demands. Each rung is self-contained.

### (a) 2 min — D3: the dimension-3 counterexample is real

```bash
cd d3 && python3 verify.py
```
Exact re-verification (independent of the announcement) that the published
`F: C³ → C³` has constant Jacobian `det = -2` yet is non-injective, refuting
JC(3). Expected output:

```
D3 counterexample — exact verification
  [PASS] Jacobian determinant == -2   (got -2)
  [PASS] three source points are pairwise distinct
  [PASS] F(0, 0, -1/4) == (-1/4, 0, 0)
  [PASS] F(1, -3/2, 13/2) == (-1/4, 0, 0)
  [PASS] F(-1, 3/2, 13/2) == (-1/4, 0, 0)
  [PASS] fiber over (-1/4, 0, 0) has exactly the 3 collision points (found 3)

RESULT: all checks passed — JC(3) is refuted by this map.
```
Runtime ~15–20 s wall (mostly the Gröbner fiber computation). Background:
`d3/NOTES.md`.

### (b) 5 min — D2: the f37 branch is a resultant artifact

**This is the single highest surprise-per-minute check in the repo.** It is the
D2 headline result: the entire `f37` branch of the (72,108) reduction is a
classical resultant excess factor, so the case collapses from three branches to
the `f31` branch alone.

```bash
cd d2 && python3 f37_sat_verify.py
```
It reads the four cofactors that Singular's `lift()` produced and verifies,
purely in SymPy (not trusting the Gröbner engine), the exact polynomial identity
`f31 = c1·G1 + c2·G2 + c3·G3 + c4·(G5body+Φ)` — i.e. `f31` lies in the
pre-resultant ideal, hence vanishes on the entire variety, over every field.
Expected output:

```
(A) membership certificate PASS:
    f31 = c1*G1 + c2*G2 + c3*G3 + c4*(G5body+Phi)  [exact over Q]
    => f31 vanishes on the entire pre-resultant variety.
(B) master identity consistency PASS:
    master identity = f31 * (f37*dm1^21); f31 in <G-system> by (A),
    so f37 and dm1^21 are resultant excess (add no ideal content).

CONCLUSION: f31 in <G1,G2,G3,G5body+Phi>.  Every solution of the
pre-resultant system has f31 = 0; the f37 branch off {f31=0} does not
lift.  The whole f37 component is a resultant artifact.
```
Runtime ~2–2.5 min (the certificate has degree-27 cofactors with thousands of
terms; the SymPy re-expansion is the cost). Full writeup:
`d2/F37_SATURATION_REPORT.md`. To reproduce the Gröbner facts
themselves (optional, needs Singular 4.2.1): `d2/f37_sat_confirm.sing`.

Note: this checker unpickles `t4_state.pkl` to regenerate the generators
`G1..G5`. If you prefer not to trust a shipped pickle, `regenerate_system.py`
rebuilds that state from scratch first (see `d2/T4_STATE_PROVENANCE.txt`).

### (c) ~30+ min — the full suite

```bash
./run_tests.sh                 # from the repo root
# SKIP_SLOW=1 ./run_tests.sh   # skips only the ~few-min numeric positive control
```
This runs (a) and (b) plus ~30 more exact D2 checkers, then a numerical
positive control. What it covers:
- the field-split repair and split-place ledgers (`t5_split_place_verify.py`,
  `test_split_place_*.py`);
- the cascade engine and its **spec-only independent auditors**
  (`audit_cascade_kills.py`, `audit_cascade_kills_sub1.py`,
  `audit_tplace_cases.py`, `audit_inf_cases.py`, `audit_alt_regime.py`,
  `audit_convolution_kills.py`, `audit_convolution_kills_r2.py`,
  `audit_reconstruction_kills.py` — each separately authored with no access to
  the engine code it checks);
- the infinity (max-plus) layer and its tie equations (`test_cascade_inf.py`,
  `cascade_inf_ties_verify.py`);
- the per-cell exact kills (`t5_90t1_*`, `t5_90t2_verify.py`,
  `t5_t2_column_verify.py`, `t5_t2_infinity_verify.py`);
- the alternate regime (`alt_regime*_verify.py`, `alt_inf_sweep_verify.py`);
- the residue-lemma library, envelope bounds, and the f37 free-family and
  saturation results.
It ends with `ALL TESTS PASSED` (rc 0) or `TEST FAILURES`. The final
`jetlift.py control f31_sub2` step is a numeric positive control (must reach
≤ 1e-5) — evidence only, not a proof; skip it with `SKIP_SLOW=1` if you only
want the exact checks.

## Where to read next

In order:
1. `d2/PROOF_INVENTORY.md` — the single source of truth: the C0–C46 claim graph,
   each claim's checker, its independent-audit status, and a trust tier
   (1 = independently audited … 4 = published result used as stated). Start here;
   it tells you exactly which claims are load-bearing and how strong each one is.
2. `d2/STATE.md` — the full chronological derivation log behind those claims.
3. `d2/FRONTIER.md` — the machine-generated live frontier (26 + 171 + 27
   surviving cells/branches), regenerated from the cascade JSON artifacts, not
   hand-maintained.
4. `d2/WRITEUP_OUTLINE.md` — the referee-facing paper skeleton, led by the
   f37-branch closure.

A one-page reader's map is in [docs/README.md](docs/README.md).

**Honest scope.** These checks confirm what is *proven*: the D3 counterexample,
the f37-branch closure, and the individual frontier-shrinking kills. They do
**not** prove the (72,108) case is closed — the frontier above is nonempty and
JC(2) remains open. See the "What is NOT claimed" section in the top-level
`README.md`.
