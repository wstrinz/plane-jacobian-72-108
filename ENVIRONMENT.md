# Reference environment (pinned)

Exact toolchain used to run `./run_tests.sh` and rebuild the certificates.
Frozen 2026-07-23. Pinned Python packages are in `requirements-lock.txt`
(`python -m pip install -r requirements-lock.txt`); `requirements.txt` keeps the
looser floors for convenience.

## Operating system / host

- Windows 11 Pro (10.0.26200). Python and the SymPy checkers run natively on
  Windows; Singular runs inside WSL2 (Ubuntu).

## Python

| Component | Version |
|-----------|---------|
| CPython   | 3.10.6  |
| sympy     | 1.14.0  |
| mpmath    | 1.3.0   |
| numpy     | 1.26.4  |
| scipy     | 1.10.1  |
| python-flint | 0.9.0 |

The exact proof checkers (`run_tests.sh`) require only `sympy` (+ its `mpmath`
dependency). `numpy`/`scipy`/`python-flint` are used by the optimization/search
harness (`jetlift.py` positive control and the sweep tooling).

## Singular (computer algebra; optional, provenance only)

- Singular **4.2.1 (4212, 64-bit)**, built 2021-12-17, running under WSL Ubuntu.
- Linked libraries: GMP 6.2.1, NTL 11.5.1, FLINT 2.8.0.
- Invoked via `d2_plane_72_108/run_singular.sh` (pipes through `wsl.exe -d Ubuntu`).
- Singular is NOT required to run the test suite: the mandatory checkers verify
  shipped certificates in pure SymPy. Singular is only needed to *regenerate*
  the certificates from scratch (`f37_sat_confirm.sing`, `regenerate_system.py`).

## Lean (formal certificate; optional)

- Toolchain: `leanprover/lean4:v4.32.1` (pinned in `lean_certificates/lean-toolchain`).
- No mathlib dependency; the certificate library `Cert/` is self-contained.
- Build with `lake build` in `lean_certificates/`. Not part of `run_tests.sh`.

## Reproducing the environment

```
python -m pip install -r requirements-lock.txt
./run_tests.sh                      # exact suite (SymPy only)
# optional, for regeneration / formal certificate:
#   Singular 4.2.1 in WSL Ubuntu
#   elan / lean4 v4.32.1 for the Lean build
```
