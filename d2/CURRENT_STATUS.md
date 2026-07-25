# CURRENT STATUS — the (72,108) program (current truth)

> **CORRECTIONS (2026-07-24) — read before the body below, which predates them.**
>
> 1. **"51 integral families" must be read as 51 FORMAL CHAIN-CENSUS ROWS**, not
>    51 realized integral window systems. The `q_window` identity
>    `t*H - q*M = q(kappa+1) - t` is valid and unaffected. But the census supplies
>    `t` from GGV5's final-corner `l`, which is a **Laurent/ramification index**,
>    not the Laurent-**chart** exponent (an edge slope). GGV1 §8 works in `L^(1)`
>    while charting with exponent 3, proving inside one section that the two need
>    not agree. Only **3 of 39** distinct integral corner shapes satisfy the chart
>    precondition `b_0 = t(a_0-1)`. See `CORNER_RESOLVENT.md` §5.1.
> 2. **F9–F13 lack a verified geometric chart dictionary.** F9 = (56,84) is
>    GGHV's `(7,21)` row; its true reduction has `t=3, kappa=1, q=1, deg C=1`
>    (`C = y`, a monomial) against the census's `t=7, q=2`. The landed
>    `("F9",0)` grammar point is a statement about the FORMAL model at the
>    supplied parameters. See `FAMILY_GRAMMAR.md`, "F9 RESOLVED".
> 3. **No general chain-to-chart dictionary exists in the literature.** GGV1
>    §8 is hard-coded to one polygon and its Remark 8.13 states the reduction
>    cannot always be obtained. A general family compiler would require new
>    geometry, not just implementation.
> 4. **The home (72,108) case and `F2` are unaffected** — their chart and divisor
>    data come from explicit published reductions, not from the census.
> 5. **Frontier counts:** the 224-part figure below and the 209 open nodes in the
>    operational handoff are *different representations*, not competing totals —
>    209 uses the defect-0 alt-layer representation (12 surviving families after
>    3 pending closures) while the 27-branch alt sweep remains separately
>    unjoined. Neither is a correction of the other.
> 6. **New since:** the Phi-window-depth criterion (133 new state-level kills,
>    0 of 220 flag-cases closed — see `FACE_KILL_SWEEP.md`), the face-detector
>    completeness result, and the marked-polynomial discriminant **theorem**
>    (`CORNER_RESOLVENT.md` §2.5).

**As of:** 2026-07-23 (body; see corrections above). This file states only what is *currently true*: proven and
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

## 4. One-line summary

The f37 branch is closed (C11, char 0; Lean-kernel-checked certificate). The
`f31` branch is an open, explicitly enumerated **224-cell frontier** with finite
per-cell residue/degree obligations, built entirely from independently audited
kills. C0 remains open.
