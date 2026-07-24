# The plane (72,108) Jacobian case — verified reduction

Work on the **Jacobian Conjecture** (JC): *a polynomial map `F: C^n -> C^n` with
nonzero-constant Jacobian determinant is invertible.* This repository holds two
related but separate lines of work.

| dir | dim | what | status |
|-----|-----|------|--------|
| `d3/` | 3 | explicit map refuting JC(3) | **done & verified** |
| `d2/` | 2 | closing the last open plane case below degree 125 | **in progress** |

The D3 counterexample settles dimension 3 in the negative. It does **not**
resolve the plane case: JC(2) is a separate open problem, which is what the
D2 effort is about.

**Who is this for?** Mathematicians who want to *verify* the load-bearing claims
in exact arithmetic; contributors who want to attack the live frontier; and
anyone auditing AI-assisted mathematics — every claim is catalogued with its
checker and a trust tier, and the load-bearing kills are re-derived by
independent spec-only auditors.

## The 15-minute path (start here)

A skeptical mathematician can confirm the two headline results, fast, in exact
arithmetic, with just `pip install sympy`:

```bash
# (a) ~2 min — D3: the dimension-3 counterexample is real
cd d3 && python3 verify.py

# (b) ~2 min — D2 headline: the entire f37 branch is a resultant artifact
cd d2 && python3 f37_sat_verify.py
```

Both exit nonzero on any failure and need no cluster. Full ladder and expected
output: **[VERIFICATION.md](VERIFICATION.md)**.

## D3 — dimension-3 counterexample (settled)

An explicit `F: C^3 -> C^3` with `det(JF) = -2` that is non-injective (three
distinct points map to `(-1/4,0,0)`). Announced by L. Alpoge 2026-07-19; the
files here are an independent exact re-verification.

```bash
cd d3 && python3 verify.py
```
See `d3/NOTES.md`.

## D2 — the plane (72,108) case (open)

**The case.** The only open case below degree 125 for a *plane* JC
counterexample is the GGV–Horruitiner case (8,28) (arXiv:2204.14178, Prop 4.3),
left open in 2022 "for lack of computing power." The work reduces it — via the
paper's own normalization — to a window-constrained polynomial-identity
feasibility question over `K[y]` for two master factors `f31`, `f37` across two
Newton-polygon subcases (1) and (2). If infeasible everywhere, the lowest
possible plane-counterexample degree rises 108 -> 125. The problem is geometric
(it must hold over `C`), so the quartic `q` is treated as its four linear places
after base change, not as one prime; the field-split repair this required is
`d2/FIELD_SPLIT_AUDIT.md`.

**Headline theorem so far — the f37 branch is closed.** The resultant reduction
factors the master identity as `f31 · f37 · d₋₁²¹ ≡ 0`, apparently leaving three
branches. In fact the whole `f37` (and `d₋₁²¹`) factor is a **resultant
artifact**: the elimination ideal of the pre-resultant system in the
`(d̃2,d̃1,d̃0,d₋₁,Φ)` variables is *exactly* `⟨f31⟩`, and `f31` lies in the
pre-resultant ideal itself. So every solution of the original system has
`f31 = 0`, and nothing on `{f37 = 0}\{f31 = 0}` lifts — in characteristic zero
(over every `ℚ`-algebra), for both subcases. This is established by an explicit
ideal-membership certificate (Singular `lift()`), re-checked independently in
SymPy without trusting the Gröbner engine. (The scope is characteristic-zero,
not "every field": the integer certificate carries a denominator-clearing
multiplier `D = 46875 = 3·5⁶`, so dividing back to `f31 = 0` needs `3, 5`
invertible — see the field-scope note in `d2/F37_SATURATION_REPORT.md`.) See
`d2/F37_SATURATION_REPORT.md`; the
afternoon-checkable verifier is `d2/f37_sat_verify.py` (~2 min, exact). The
(72,108) case therefore reduces to the `f31` window/cascade program alone.

**The frontier.** The remaining `f31` obligation is an open, explicitly-enumerated
frontier of **26 subcase-2 cells + 171 subcase-1 branches + 27 alternate-regime
branches** (machine-generated `d2/FRONTIER.md`), each carrying finite
residue/degree obligations and all produced by independently-audited kills; a
state-level reduction (convolution + reconstruction) is in progress, with roughly
300 further kills currently in audit.

**Verification culture.** The suite is `./run_tests.sh` (D3 exact check + the D2
proof checkers). Independence is deliberate: the load-bearing cascade, infinity,
and alternate-regime kills are re-derived by spec-only auditors
(`audit_cascade_kills*.py`, `audit_tplace_cases.py`, `audit_inf_cases.py`,
`audit_alt_regime.py`, `audit_convolution_kills*.py`, `audit_reconstruction_kills.py`
— separately authored, no access to the engine code they check), wired into
`run_tests.sh`. Every claim is catalogued with its checker, audit status, and a
trust tier in `d2/PROOF_INVENTORY.md` (the single source of truth, C0–C46), and
the cross-front frontier counts are regenerated from the JSON artifacts into
`d2/FRONTIER.md` rather than hand-maintained. The referee-facing skeleton is
`d2/WRITEUP_OUTLINE.md`. A reader's map through these documents is in
[docs/README.md](docs/README.md).

Key exact entry points:
```bash
cd d2
python3 f37_sat_verify.py                    # the f37-branch closure (ideal membership)
python3 split_place_ledger.py && python3 test_split_place_ledger.py
python3 t5_90t1_constant_verify.py
```
For full regeneration and numerical reconnaissance:
```bash
python3 regenerate_system.py && bash run_singular.sh   # rebuild f31/f37 from scratch
python3 jetlift.py control f31_sub2                     # validate the harness (~1e-7)
python3 jetlift.py stats   f37_sub2 1800                # map the residual floor
```

## New (2026-07-23): the corner law — a structure layer across cases

Beyond the (72,108) kill program, this tree now carries a cross-family
**structure layer** on the GGV/GGHV corner data — a candidate seed for an
inductive treatment:

- A **corner law** for the tower's last algebraic element Φ: closed forms and a
  μ-graded signature law with **twelve exact derived points** across four
  regimes (`d2/PHI_75_125.md`, `d2/PHI_CORNER4.md`, `d2/PHI_F14.md`,
  `d2/PHI_F7.md`, `d2/ZETA_TAIL.md`), plus a proven chart theorem κ = t−2
  (`d2/COMPOSITE_CHARTS.md`).
- A **prior-art audit** against the published literature (`d2/PRIOR_ART.md`),
  including the strongest available check: GGHV22 §4 prints one closed-form
  last element (their CAS-solved `f₁` at the (9,24) corner) — and the corner
  law **postdicts their printed polynomial exactly, with zero fitting
  freedom**. That postdiction is machine-verified by
  `d2/prior_art_postdiction_verify.py`, one of ten new exact-arithmetic
  verifiers wired into `run_tests.sh`.
- A structural classification of the 23-shape residue library via Galois
  descent (`d2/GALOIS_LIBRARY.md`) and a pilot **case compiler**
  (`d2/CASE_COMPILER.md`) that reproduces the audited (72,108) facts from
  corner data alone.

These are same-author checker-enforced results (tier 2), published here for
external review; see `d2/CURRENT_STATUS.md` §3b for the honest tier/status
accounting.

## What is NOT claimed

The (72,108) case is **not closed** and nothing here proves or refutes JC(2).
The f37-branch closure and the individual cascade/infinity/alternate kills are
proven and (for the load-bearing ones) independently audited, but the target
theorem `C0` remains **open**: the frontier above is nonempty, and its per-cell
obligations are finite constraint lists, not proofs. Kills still in audit —
including the ~300 state-level convolution/reconstruction kills — are labelled as
pending and are not counted as closed. The generic numerical floors are
*evidence, not proof*. And parts of the upstream derivation (envelope/window
bounds, the Prop 4.3 transcription) still rest on same-author checks rather than
independent audit; the exact verified/unverified split and ranked risks are in
`d2/AUDIT.md`, the full log in `d2/STATE.md`.

## Setup & tests

```bash
./setup.sh        # pip deps (numpy/scipy/sympy/mpmath) + Singular (D2 regen only)
./run_tests.sh    # D3 exact verification + D2 proof checkers + harness control
# SKIP_SLOW=1 ./run_tests.sh   # skips only the ~few-min numeric positive control
```
Python 3.10+. Singular is only needed for the D2 system regeneration; the D3
verification and the D2 exact checkers run without it.

## Provenance

Will Strinz + Claude, chat/agent sessions, July 2026. The underlying research
program and the plane open question are due to Guccione–Guccione–Horruitiner–
Valqui (arXiv:2204.14178); any writeup of a D2 result belongs with or
alongside them.
