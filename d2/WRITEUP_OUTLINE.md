# WRITEUP OUTLINE — the (72,108) case, referee-facing paper skeleton

**Purpose.** A section-by-section skeleton for the eventual paper proving: *no plane
Jacobian counterexample of the (72,108)/(8,28) type (Prop 4.3 subcases (1),(2)) exists*,
raising the published JC(2) lower bound 108 → 125. Organized per gap **G5**: the paper
leads with the **f31-only reduction** (ideal membership, C11) as the organizing theorem;
the resultant chain and the entire f37 narrative (C9, C12, C36) are retired to a single
remark. This is a **drafting scaffold only** — it edits nothing and asserts no new math.

**Provenance.** Built from `PROOF_INVENTORY.md` (the C0–C46 claim graph and its tier
legend), `STATE.md` (through the 2026-07-22 infinity-layer campaign, L662–879), and the
machine-generated `FRONTIER.md`. Honest-status labels use three values:
**PROVEN** (closed, with a tier-1/2 checker or audit),
**PROVEN-PENDING-FRONTIER** (the machinery is proven but its purpose is to shrink an
open frontier that is not yet empty), **OPEN** (not yet a theorem).

**Headline status.** The target theorem **C0 is OPEN**: the frontier is nonempty. As of
the machine-generated `FRONTIER.md`, the live fronts are

| Front | Size | Flag cases | Tier of the kills that produced it |
|---|---|---|---|
| Sub2 cells | **26** (18 T1 + 8 T2), a∈[5,10] | 220 | 1 (cascade + infinity, independently audited) |
| Sub1 standard regime | **171** branches (118 T1 + 53 T2), a∈[2,10] | 1145 | 1 (cascade + infinity, independently audited) |
| Sub1 alternate regime | **27** branches (13 T1 + 14 T2) | (38360 degree-states) | 1 (flipped cascade, independently audited) |
| **Total** | **224 surviving cells/branches** | **1365** | — |

Every closure that produced this frontier is tier 1 or tier 2. What remains (G1) is to
discharge these 224 obligation-carrying survivors — the obligations are finite and
explicit, but *obligations are not proofs* (both Phase-C worklists say so).

---

## Front-matter table — every C-id → section → tier

Legend: **CP** = on the critical path to C0 via the live frontier. **⚠** = tier-2\*/3
item sitting on the critical path — a **referee risk** to be addressed before submission.
Tiers are as recorded in `PROOF_INVENTORY.md §1`; where the inventory tier is **stale**
because later work upgraded it, the current strength is noted in the last column.

| C-id | Claim (short) | Section | Inventory tier | CP? | Note / referee risk |
|---|---|---|---|---|---|
| C0 | Target theorem | 1, 8 | open | — | the theorem itself; OPEN (G1) |
| C1 | Prop 4.3 transcription | 1 | 3/4 | CP ⚠ | document-only; **not** covered by the envelope checker; needs referee-grade LaTeX-vs-PDF audit (part of G2 residue) |
| C2 | Normalization ℓ(P)=R², ℓ(Q)=R³ | 1 | 2/4 | CP | `t6_premises_verify.py` + GGV1 Props 1.13/2.1 |
| C3 | α-strip WLOG, v(F)=−5 | 1 | 2/4 | CP | `t6_premises_verify.py` + GGHV22 §4 |
| C4 | Forcing ODE, unique f₁ | 2 | 2\* | CP ⚠ | `verify_derivation.py` exists but **unwired** (I8) |
| C5 | D-transformation recursion | 2 | 2\* | CP ⚠ | `verify_derivation.py` unwired (I8) |
| C6 | Envelope/window bounds (14w/15w, 12w) | 2 | 3 | CP ⚠→**2** | **UPGRADED**: `envelope_bounds_verify.py` (+`ENVELOPE_BOUNDS.md`), 106/106, referee-grade — G2 closed. Flag the GGHV22 sign typo (L1462–1466). |
| C7 | 12-equation window selection | 2 | 2\* | CP ⚠ | `verify_derivation.py` (48 checks) unwired (I8) |
| C8 | System regeneration reproducible | 2 | 2\* | CP ⚠ | `regenerate_system.py`; validated vs published (7,21) |
| C9 | Master resultant identity | 3 (remark) | 3 | — | retired to a remark; superseded in load-bearing role by C11 |
| C10 | d₋₁≡0 impossible | 3 (remark) | 2\* | — | historically on path; subsumed by C11 |
| **C11** | **f31-only reduction (ideal membership)** | **3** | **2** | **CP** | **organizing theorem**; `f37_sat_verify.py` (♦), cross-derived two ways |
| C12 | f37 free family no-lift | 3 (remark) | 2 | — | subsumed by C11 |
| C13 | Graded cascade decomposition | 4 | 2\* | CP ⚠ | `verify_graded.py`/`t5_multiplace_verify.py` **unwired** — foundation of every cascade |
| C14 | Field-split repair | 4 | 1 | CP | externally authored audit + local checkers |
| C15 | Split-place sigma-locus theorem | 4 | 1 | CP | field-stable; Mason–Stothers |
| C16 | a_t=7 stratum dead | 4 | 1 | CP | externally authored derivation + checker |
| C17 | Sub2 terminal ledger | 5 | 1 | CP | independently audited (654/654) |
| C18 | Sub2 depth-4 cascade (30 survive) | 5 | 1 | CP | engine + `audit_cascade_kills.py` |
| C19 | Sub2 cone-lemma compression | 5 | 2 | CP | compression, explicitly NOT independent |
| C20 | Sub2 t-place coupling | 5 | 1 | CP | `audit_tplace_cases.py` (**unwired**, I8) |
| C21 | a_t=9 T2 infeasible | 6 | 2 | CP | same-author checker |
| C22 | a_t=9 T1 constant-E infeasible | 6 | 2 | CP | same-author checker (constant kill) |
| C23 | a_t=9 T1 nonconstant reduction | 6 | 2 | CP | reduction only; **open cell** in frontier |
| C24 | T2-column squeeze (4 cells) | 6 | 2 | CP | same-author checker; corrects PHASE_C_WORKLIST §4 |
| C25 | T2 infinity convolution narrowing | 6 | 2 | CP | same-author checker |
| C26 | Sub2 Phase C worklist | 8 | 3 | CP ⚠ | obligation list, judgment-tagged; **not a proof**; headline counts stale (I1) |
| C27 | Sub1 cascade parameters (+correction) | 5 | 1 | CP | independent Codex cap correction |
| C28 | Sub1 terminal ledger | 5 | 1 | CP | independently audited (2614/2614) |
| C29 | Sub1 depth-4 cascade (279 survive) | 5 | 1 | CP | engine + `audit_cascade_kills_sub1.py` |
| C30 | Sub1 cone lemmas | 5 | 2 | CP | compression, engine-derived |
| C31 | Sub1 t-place coupling | 5 | 1 | CP | `audit_tplace_cases.py` (**unwired**, I8) |
| C32 | Sub1 Phase C worklist (91% reuse) | 8 | 3 | CP ⚠ | obligation list, judgment-tagged; **not a proof** |
| C33 | Alternate regime flipped cascade | 7 | 2→**1** | CP | **UPGRADED**: `audit_alt_regime.py` → tier 1 (G3 closed) |
| C34 | Flipped cascade levels 3/2 (27 branches) | 7 | 2→**1** | CP | **UPGRADED** via `audit_alt_regime.py` |
| C35 | Pre-repair strata kills (scoped) | 6 (remark) | 3 | — | superseded by C17–C18; conditional (I5) |
| C36 | f37 graded/frontier analysis | 3 (remark) | 3 | — | MOOT — superseded by C11 (I3) |
| C37 | Numerical infeasibility floors | App. C | 3 | — | evidence only (G8); motivation, not proof |
| C38 | NulLA/linear-certificate no-go | App. B | 3 | — | negative method result; not load-bearing |
| C39 | Corner-144 comparison | 7 (remark) | 2\* | — | forward-looking; not load-bearing |
| C40 | Next-cases enumeration | 1 (remark) | 2\*/4 | — | context; not load-bearing |
| C41 | Max-plus infinity layer | 5 | 2 | CP | same-author checker; regression ladder R0–R5 |
| C42 | Infinity tie equations | 5 | 2 | CP | leans on `residue_lemmas_verify` (tier 2) |
| C43 | Infinity sweeps, both windows (26/171) | 5 | 1 | CP | engine + `audit_inf_cases.py` — **produced the current frontier** |
| C44 | Alternate-regime degree layer | 7 | 1 | CP | `audit_alt_regime.py` |
| C45 | Phase D residual worklist (sub2) | 8 | 3 | CP ⚠ | engine-derived worklist; **not a proof artifact** |
| C46 | No-jet-kill theorem | 5 | 2 | CP | `residue_lemmas_depth_verify.py`; steers Phase D to global/infinity |

**Referee-risk summary (the ⚠ rows on the critical path):**
- **Upstream derivation chain (C4, C5, C7, C8, C13):** tier-2\* — exact checkers exist
  (`verify_derivation.py`, `verify_graded.py`, `t5_multiplace_verify.py`) but are **not in
  the authoritative suite** (I8). Every frontier kill sits on this chain. **Fix: wire them.**
- **C1 (Prop 4.3 transcription):** document-only (3/4); the envelope checker does *not*
  cover it. Needs a referee-grade source audit + the GGV3 substitution disclosure (G7).
- **C6 (envelope bounds):** the inventory still lists tier 3, but `envelope_bounds_verify.py`
  now closes **G2** at referee grade (106/106). Update the tier when the inventory is regenerated.
- **C26, C32, C45 (Phase C/D worklists):** tier 3 by nature — they *are* the frontier's
  obligation lists, not closures. Their presence on the critical path is exactly gap G1.

---

## Section 1 — Setup and the published case

**(a) Content.** State the plane Jacobian problem and GGV-Horruitiner (arXiv:2204.14178)
Prop 4.3, case (8,28): the last (72,108) subcase left open below 125 "for lack of computing
power." Transcribe both Newton polygons (subcases (1),(2)) corner-for-corner and fix the
[P,Q]=x² normalization ℓ(P)=R², ℓ(Q)=R³, R=x⁴C₄, C₄=y⁷(y+1), plus the α-strip WLOG
Q=C³+λC⁻¹+F with v(F)=−5. Close with the payoff statement: if the case dies, the bound
rises 108→125.
**(b) Consumes.** C1 (case statement + polygons), C2 (normalization), C3 (α-strip WLOG);
context remark C40 (next-cases table). Tier-4 citations: GGHV22 Prop 4.3, GGV1 Props
1.13/2.1, GGHV22 §4 (the closed (9,27) twin used as the setup template).
**(c) Tiers / what a referee wants.** C2, C3 are tier 2/4 — checker `t6_premises_verify.py`
(14 exact checks) plus published propositions used as stated. C1 is tier **3/4**,
document-only: a referee will want the LaTeX-vs-published-PDF transcription audit and the
**G7 GGV3 note** — GGHV22 cites GGV3 §1 for the strip descent, whose source is not in
`paper_src/`; here it is filled from GGV1 Props 1.13/2.1, and that substitution **must be
stated explicitly** (or verified against GGV3).
**(d) Status.** **PROVEN** as setup (modulo the two tier-4 citations used as stated and the
G7 disclosure). The *case* is OPEN — that is the whole paper.

## Section 2 — The window system and its selection

**(a) Content.** Derive the forcing ODE 8y(y+1)f₁′−14(8y+7)f₁=y⁸(y+1)² with its unique
separable-quartic solution; the D-transformation D_k=C_k·C₄^(7−2k)∈K[y]; the **envelope
bounds** deg d_{4−k}≤14w (sub2)/15w (sub1), ord≥12w; and the selection of the 12-equation
window system (k=8 vacuous, λ isolated in (D̃³)₋₄, dropped equations define only spare
unknowns). End at Φ:=f₁C₄³¹∈K[y] with deg 238, ord 204, mult_{y+1} 30.
**(b) Consumes.** C4 (ODE), C5 (D-recursion), C6 (envelope bounds), C7 (equation selection),
C8 (regeneration). Tier-4: GGHV22 Props 5.2(3)(4)/5.6 (induction template — redone here).
**(c) Tiers / what a referee wants.** C4/C5/C7/C8 are **tier 2\*** — `verify_derivation.py`
(48 symbolic checks) and `regenerate_system.py` exist and pass but are **not wired into the
suite (I8)**; a referee wants them wired and re-run. C6 was the old G2 gap; it is now
**closed** by `envelope_bounds_verify.py` (referee-grade, 106/106 checks, w=0..5, named
assumptions A1–A3, +`ENVELOPE_BOUNDS.md`). **Disclose** the sign typo it found in GGHV22
L1462–1466 (mirrored in `T3_WINDOW_AUDIT.md` L55): the C-recursion is (P−sum)/(2C₄), not
−(P+sum)/(2C₄); valuation proof unaffected.
**(d) Status.** **PROVEN** (window bounds now referee-grade; note the I8 wiring debt for the
selection checkers).

## Section 3 — The f31-only reduction as the organizing theorem

**(a) Content.** This is the paper's pivot (per **G5**). State and prove the ideal-membership
theorem: over Q with Φ a free indeterminate, f31 ∈ ⟨G1,G2,G3,G5body+Φ⟩, and the elimination
ideal in (d̃2,d̃1,d̃0,d₋₁,Φ) is **exactly** ⟨f31⟩. Consequently every solution of the original
system satisfies f31=0; the classical resultant's extra factors f37 and d₋₁²¹ are **excess
artifacts** carrying no geometry — uniformly, for **both** subcases. Reduce the entire case to:
does f31=0 admit an admissible window solution?
**(b) Consumes.** C11 (the theorem). **Remark** (retiring the old narrative): C9 (resultant
identity), C10 (d₋₁≡0 impossible), C12 (f37 free-family no-lift), C36 (f37 graded analysis) —
all now subsumed or moot; present the resultant route as historical motivation only.
**(c) Tiers / what a referee wants.** C11 is **tier 2**: `f37_sat_verify.py` (♦) re-expands
the `lift()` cofactor certificate in sympy **without trusting the Gröbner engine**, and an
independent H-system saturation route confirms it. A referee wants the cofactor certificate
in an appendix and the two-route agreement noted. This section *dissolves* the referee's
"do you trust Singular's resultant?" question (C9's tier-3 weakness).
**(d) Status.** **PROVEN.**

## Section 4 — The graded cascade identity and the field-split framework

**(a) Content.** The graded decomposition f31 = Σ_{f=0..7} Φ̃^f d̃₋₁^{21−3f} h_f with explicit
small h_f and unit cofactor −q/6630, so the cascade runs simultaneously at the t-place and at
the q-places. Then the **field-split repair**: q splits into four geometric places after base
change, so the scalar (a_t,a_q) ledger is not field-stable; the correct case space is
multiplicity vectors (b₁..b₄) up to S₄. Include the field-stable **sigma-locus theorem** (T3
branch empty over every char-0 field) and the death of the geometrically q-coprime a_t=7
stratum.
**(b) Consumes.** C13 (graded decomposition), C14 (field-split repair), C15 (sigma-locus),
C16 (a_t=7 dead). Tier-4: Mason–Stothers.
**(c) Tiers / what a referee wants.** C13 is **tier 2\*** (`verify_graded.py`,
`t5_multiplace_verify.py` — **unwired**, I8) and is the load-bearing foundation of every
cascade; **flag for wiring**. C14–C16 are **tier 1** (externally authored `FIELD_SPLIT_AUDIT.md`
+ local exact checkers `t5_split_place_verify.py`). A referee wants the field-split repair
narrated *before* any ledger counts, since it invalidated the earlier (a_t,a_q) maps (I5).
**(d) Status.** **PROVEN.**

## Section 5 — The descent engines (finite places, infinity, tie/residue layer)

**(a) Content.** The mechanized q-place cascade (levels 7/6/5/4, ultrametric semantics with
recorded residue obligations, budget-coupled four-place join) plus the t=y+1 fifth place and
the **max-plus infinity sixth place** (v_inf=−deg). Present the terminal ledgers, the depth-4
sweeps, cone-lemma compression, and the infinity sweeps that produced the **current
frontier** (sub2 26 cells/220 flag cases; sub1 279→**171** branches, 108 new infinity kills,
a=0/1 closed). Include the tie/residue layer (residue-lemma library: 21 constraint lemmas +
2 kills C08/C20) and the **no-jet-kill theorem**: t-place obligations are constraints for all
a, never local kills — so remaining closures must be global/infinity.
**(b) Consumes.** Finite places: C17–C20 (sub2 ledger, cascade, cone lemmas, t-coupling),
C27–C31 (sub1 analogues). Infinity: C41 (layer), C42 (tie equations), C43 (sweeps — produced
the frontier). Tie/residue: C42 + RESIDUE_LEMMAS (C08/C20), C46 (no-jet-kill).
**(c) Tiers / what a referee wants.** The engine kills are **tier 1** across the board —
independently audited by the Codex spec-only checkers `audit_cascade_kills.py` (654/654,
420/420), `audit_cascade_kills_sub1.py` (2614/2614, 2178/2178), and `audit_inf_cases.py`
(sub2 420/420, sub1 2178/2178) — the four (now five) spec-only auditors are the paper's
tier-1 backbone. **Referee risk:** C20/C31 rest on `audit_tplace_cases.py`, which is **the
only tier-1 artifact outside the suite (I8)** — wire it. C19/C30 (cone lemmas) and C41/C42
(infinity semantics) are tier 2 (compression / same-author); C46 is tier 2. A referee wants
the code-free semantics specs (`CASCADE_ENGINE_REPORT.md`, `CASCADE_INF_REPORT.md`) that the
auditors were written against, to confirm auditor independence.
**(d) Status.** **PROVEN-PENDING-FRONTIER** — the engines are sound and audited; their
*output* is the nonempty 26+171 frontier (G1).

## Section 6 — Per-cell closures and the mechanized convolution/elimination layer

**(a) Content.** The exact post-cascade cell kills: a_t=9 T2 (constant E + infinity
domination), a_t=9 T1 constant-E (exact convolution 238→226 to a nonzero constant), the
a_t=9 T1 nonconstant reduction, the T2-column squeeze (4 of 12 cells), and the T2 infinity
convolution narrowing. Then the **mechanized layer**: `convolution_descent.py` (rediscovers
the a9 constant chain autonomously) and `convolution_elim.py` (Gröbner elimination modulo the
accumulated coefficient ideal with Rabinowitsch nonzero constraints — mechanically kills the
deg σ=5,4 tied cells).
**(b) Consumes.** C21 (a9 T2), C22 (a9 T1 constant), C23 (a9 T1 nonconstant reduction), C24
(T2-column), C25 (T2 infinity). **Remark:** C35 (pre-repair strata, scoped/conditional).
**(c) Tiers / what a referee wants.** C21–C25 are **tier 2**, same-author exact checkers
(`t5_90t2_verify.py`, `t5_90t1_verify.py`, `t5_90t1_constant_verify.py`,
`t5_t2_column_verify.py`, `t5_t2_infinity_verify.py`, all parsing `f31_graded.txt`, no
hand-copied coefficients). This is the **G3-flavoured** referee risk for the *cell* layer:
no spec-only auditor covers C21–C25 (the three standard-regime auditors are cascade-level).
Note C23 is a *reduction*, not a closure — its cell is still open in the frontier.
**(d) Status.** C21, C22, C24 closures **PROVEN** (tier 2); C23, C25 **PROVEN-PENDING-FRONTIER**
(narrowing, not closure).

## Section 7 — The alternate regime (sub1, a∈[11,15])

**(a) Content.** For v=30−3a<0 the t-order telescopes to the constant 210 (F=t²¹⁰G′), giving
a **descending** flipped cascade anchored at h₇=8192·d1²; q-place transitions and terminal
caps survive verbatim. First-level parity/degree kills 19/52 branches; levels 3/2 kill 6 more
→ **27 residual branches** (13 T1 + 14 T2). The alternate-regime degree (infinity) layer kills
33670/38360 states, leaving 27 branches open with explicit obligations.
**(b) Consumes.** C33 (flipped cascade), C34 (levels 3/2), C44 (degree layer). **Remark:** C39
(corner-144 comparison, forward-looking template evidence).
**(c) Tiers / what a referee wants.** **UPGRADED to tier 1**: `audit_alt_regime.py`
(Codex spec-only, no checker/engine access) re-derived all 25 branch kills and matched all 27
OPEN verdicts with complete state partitions — **G3 is closed**. A referee should be told the
provenance caveat (I11: the L2 authoring runtime lost its job registry; deliverables verified
locally) is now *superseded* by the independent audit. Writeup note: `ALT_REGIME_L2.md` lacks
the 14-branch T2 subsection (recoverable from its verdict table).
**(d) Status.** **PROVEN-PENDING-FRONTIER** — machinery audited (tier 1); 27 branches open.

## Section 8 — The frontier and what remains (honest)

**(a) Content.** State the frontier exactly from the machine-generated `FRONTIER.md`: **26 sub2
cells + 171 sub1 branches + 27 alternate branches = 224 survivors, 1365 flag cases**, each
carrying a finite, explicit obligation list. Explain honestly that these obligations are
*constraints and degree-states*, not yet contradictions; the no-jet-kill theorem (C46) says
the remaining kills must come from global/infinity/elimination arguments, not deeper t-jets.
State what a completed proof requires (see blocking list).
**(b) Consumes.** G1 (the open frontier); the worklists C26, C32, C45 (Phase C/D obligation
inventories). Cross-references the Section 5/6 engines that will discharge them.
**(c) Tiers / what a referee wants.** C26, C32, C45 are **tier 3** by nature — obligation
lists, judgment-tagged, explicitly not proofs. A referee wants the honest statement that the
theorem is not yet proven and the precise, machine-generated survivor list (drift-free per G6:
`FRONTIER.md` is regenerated from the JSON artifacts, not hand-maintained).
**(d) Status.** **OPEN** (this is gap G1 — the primary gap; everything else is presentation).

## Appendices

- **A. Reproducibility (the suite map).** Map each claim to its checker and the authoritative
  runner. The suite is `run_tests.sh` at the REPOSITORY ROOT (`math-stuff/run_tests.sh`,
  git-tracked; it cds into `d2_plane_72_108/` — a draft of this outline mistakenly looked for
  it inside the subdirectory). The appendix should state the root-relative invocation
  explicitly. Fold in the I8 wiring list (`verify_derivation.py`,
  `verify_graded.py`, `t5_multiplace_verify.py`; `audit_tplace_cases.py` has since been wired ♦).
- **B. Audit methodology (the spec-only auditors).** Document the independence discipline: the
  Codex-authored auditors — `audit_cascade_kills.py`, `audit_cascade_kills_sub1.py`,
  `audit_tplace_cases.py`, `audit_inf_cases.py`, `audit_alt_regime.py` — are written from the
  code-free semantics reports only, with their own `f31_graded.txt` parsers and homogeneity
  self-checks, no access to `cascade_engine.py`. Plus the externally authored
  `FIELD_SPLIT_AUDIT.md`. Note C38 (NulLA no-go) here as the negative method result.
- **C. Numerics as motivation (per G8).** The four factor×subcase window systems each show a
  genuine positive floor under central-difference polish (ordering f37_sub1 5e-7 < f37_sub2
  1.2e-6 < f31_sub1 8e-6 < f31_sub2 2.3e-5), replicated on two machines (C37). Scope this
  **explicitly as motivation, not proof**; state that base-region completeness (T4) was not
  attempted — a referee will ask which, and the answer is "motivation."

---

## Blocking list — exactly what must change for the paper to exist

Ordered by severity. The paper cannot claim C0 until **B1** is resolved; the rest are
referee-grade hardening that the paper needs to *withstand review* even of a conditional result.

1. **B1 — G1: close the frontier (the theorem itself).** Discharge all **224** survivors:
   26 sub2 cells + 171 sub1 branches + 27 alternate branches (1365 flag cases). Their
   obligations are finite and explicit, but none is yet a contradiction. Per C46 the closures
   must come from the global/infinity/elimination layer (Section 5/6 engines: the Phase D
   infinity descent + `convolution_elim.py`), not deeper t-jets. **Until B1, C0 is OPEN and
   the paper is at best a conditional/partial result.**

2. **B2 — I8: wire the remaining derivation checkers into the suite.** The authoritative
   suite `run_tests.sh` EXISTS at the repository root (git-tracked; a draft of this outline
   looked for it in the subdirectory). The real debt is the unwired ◇ checkers: include
   `verify_derivation.py` (C4/C5/C7/C10), `verify_graded.py` (C13), `t5_multiplace_verify.py`
   (C13), and above all **`audit_tplace_cases.py`** — the *only tier-1 artifact outside the
   suite*, backing C20/C31 on the critical path.

3. **B3 — G2 residue: the Prop 4.3 transcription (C1).** The envelope-bound half of G2 is
   **closed** (`envelope_bounds_verify.py`, referee-grade — confirm STATE.md L688–695). What
   remains is C1: the Newton-polygon transcription is still document-only (tier 3/4) with no
   checker. Provide a referee-grade source audit, and **disclose the GGHV22 sign typo**
   (L1462–1466) the envelope checker surfaced.

4. **B4 — G7: the GGV3 citation debt.** GGHV22 cites GGV3 §1 for the strip descent; the
   source is not in `paper_src/` and the argument is filled from GGV1 Props 1.13/2.1. The
   paper must either verify against GGV3 itself or state the substitution explicitly. Related
   standing obligation: contact the GGHV/GGV authors before writeup (their program, their open
   question).

**Already cleared (do not re-litigate):**
- **G2 (envelope bounds)** — closed by `envelope_bounds_verify.py` / `ENVELOPE_BOUNDS.md`
  (106/106; STATE.md L688–695).
- **G3 (per-cell / alternate-regime independent audit)** — the alternate-regime chain is now
  tier 1 via `audit_alt_regime.py` (STATE.md L851–862); the infinity sweeps are tier 1 via
  `audit_inf_cases.py`. *Residual G3:* the per-cell kills C21–C25 (Section 6) still rest on
  same-author checkers only — a spec-only cell auditor would be the last piece, but it is
  hardening, not a B1-class blocker.
- **G5 (narrative restructure)** — this outline enacts it: lead with C11, retire f37 to a
  remark.
- **G6 (drift-free frontier)** — `FRONTIER.md` is machine-generated from the JSON artifacts.
