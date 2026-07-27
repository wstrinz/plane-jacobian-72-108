# v0.4.0 — RETRACTION: the period-12 window functions for (75,125) are withdrawn

**Read this section before any other part of this tree.** This release exists
primarily to correct published mathematics, not to add features. If you cited or
built on the (75,125) window material from v0.3.0, it was wrong, and the corrected
statement is below.

## The retraction

v0.3.0 shipped a section headed *"Period-12 window functions for (75,125)"*, and
`d2/WINDOW_FUNCTIONS_75_125.md` was titled *"The period-12 window functions for
(75,125), derived exactly."* It claimed:

| v0.3.0 claimed | actually |
|---|---|
| `q_window = 12`, period 12 | **`q_window = 29`** |
| `W_step = ord_y(Phi)/M = 201/36 = 67/12` | **`80/29`** |
| `alpha = 67 = 10a^2-8a+1`, `beta_m = 5m mod 12` | `alpha = 80`; the `beta` law is void as stated |
| upper cap `U(w) = 14w` is **affine** | **there is no affine upper cap** — `deg_slope = 80/29` is not an integer |
| family law `q_window = 5a-3` in `{7,12,17,...}` | **`12a-7` in `{17,29,41,...}`**, and both `17` and `29` are prime |
| `N` at (75,125) `= 98`; at (50,75) `= 36` | **`77`** and **`28`** |
| `Phi = -(1/9) y^201 (y^3+1)^101` | **`(1/3) y^80`** — a monomial |

**Root cause.** The `(5,20)` corner does not satisfy the retraction shape
`b0 = t(a0-1)`: `20 != 4*4`. GGV5's final-corner dictionary
`(t,q) = (l_final, b_final)`, which we used to read chart data off the published
final chain corner, is valid **only** on that shape. The correct chart exponent is
`t = ceil(b0/a0) = 4`, not `5`, and `C = y` is a **monomial**, so
`deg C = ord C = 1` rather than `5` and `2`.

**Decisive external check.** GGV3 (`arXiv:1406.0886`) section 5 performs this very
reduction on the sibling `(50,75)` and publishes `[P_1,Q_1] = x^2`,
`deg P_1 = 10`, `deg Q_1 = 15`. `t = 4` reproduces all three; `t = 5` contradicts
all three.

**The corrected reading is stronger, not weaker.** Because `q_window = 29 = M`
*exactly*, no split with `1 <= w < M` is a multiple of the period, so the carry
obstruction is `1` on **every** admissible split and there is no weight at which it
could vanish. Under the superseded `67/12` against `M = 36` it would have vanished
at `w` in `{12,24}`. So (75,125) is *more* firmly obstructed by this mechanism than
v0.3.0 claimed.

The root cause is now **guarded in code**:
`polygon_reduction.final_corner_dictionary()` raises off the retraction shape, so
this class of error cannot recur silently.

## Second retraction: the F2 family closed form

`FAMILY_GRAMMAR.md` classified the F2 family as `CLOSED-FORM (pure)` with
`f = -1/(a*dg) * y^rho * (y^dg+1)^e`. Under the corrected chart,
`dg = deg C - ord C = 0`, so that normalising constant is `-1/(a*0)` —
**undefined** — and the collapse identity `y g' - dg*g = -dg` degenerates to
`0 = 0`, constraining nothing.

**F2 has no closed form of this shape.** It is reclassified `CHART-DEGENERATE`, and
the census of pure families drops from `{F2,F9,F14}` to `{F9,F14}`. The landed
derived points still reproduce exactly and symbolically in `j`
(`N = (3j+4)(4j+7)`, `deg = ord = 2(2j+3)(3j+5)`), and independently agree with
`window_functions_75_125.family()`. What is withdrawn is the *mechanism*, not the
data. Deriving what replaces it is open.

**`F3` is repaired too — and the repair carries its own two-sided proof.** `F3`
sits at the same `(5,20)` corner and carried the same unrepaired `t = 5`. An
earlier draft of these notes said no corrected target existed for it, so no fix
would be verifiable. That was too conservative:

- **Chart data is a property of the corner, not the family.** `chart_exponent`,
  `kappa`, and the test deciding whether `C` is a monomial are functions of `A_0`
  alone. `corner_chart_data(5,20, b_final=2)` and `(..., b_final=3)` return
  bit-identical results.
- **The old table's own inconsistency is the proof.** `ord_y(C)` is a corner
  invariant, yet the table asserted `q = 2` for F2 and `q = 3` for F3 at the *same*
  corner. That contradiction is visible without consulting any paper, and it is
  exactly the fingerprint of reading `q` off GGV5's per-row `b_final`.
- **Two-sided cross-check.** `F3 j=0` is `(75,50)`, which is `(50,75)` with
  `P` and `Q` exchanged at the same corner — the same reduction. So F2 `j=0` and
  F3 `j=0` must agree, and they do: both give `N = 28` and signature
  `(30,30,0,0)`. They are reached from *different* `(m,n)` laws whose `N`
  polynomials in `j` are genuinely different — `12j²+37j+28` versus
  `84j²+99j+28` — coinciding only at `j = 0`. That is independent confirmation,
  not a tautology.
- **Published anchor.** GGV3 §5's `[P₁,Q₁] = x²`, `deg P₁ = 10`, `deg Q₁ = 15` are
  obtained *before* the paper's `γ ∈ {2,3}` branch, so they are not the property of
  one row's `b_final`; and `F_3(3,2)/75` is that case with `P,Q` exchanged. Either
  reading forces `kappa = 2`, `l = 4`.

`F3` is therefore also `CHART-DEGENERATE`, and its closed form is withdrawn on the
same grounds as F2's. The `A10` tripwire was **replaced, not deleted**: `A10a`–`A10f`
now assert the new state and its evidence, and `A10d` *discriminates* the two charts
rather than relabelling — F3's `PHI_F7` polynomial does solve the superseded ODE
(that computation was never wrong) and does **not** solve the repaired one, so the
repair stays falsifiable. New check `A9c` cross-checks the two independent copies.

**Known-suspect, and shipping in this tree.** Eight modules still transcribe `F3`'s
`(5,20)` data with `l_final = 5` / `b_final = 3` used **as chart data**:
`phi_f7.py`, `phi_f7_verify.py`, `phi_corner4.py`, `phi_corner4_verify.py`,
`phi_f14.py`, `phi_f14_verify.py`, `case_compiler.py`, `ml_restriction_check.py`.
Each keeps its own transcription, so each is internally consistent and all stay
green — which is precisely why they are dangerous. Check `A10f` names them so the
front cannot be mistaken for closed. Treat any `F3` chart quantity from those
eight as suspect until repaired.

A latent bug surfaced during this work and is fixed: `full_ode_residual` built
`g = y**dg + 1 = 2` at `dg = 0` instead of the forced monic constant `1`. It had
never been exercised, because no `dg = 0` row had ever reached that line.

## How these survived — and what changed structurally

Both retractions were found by adversarial review, and in both cases the test
suite was green throughout. Three distinct mechanisms allowed that:

1. **A module printed `MISMATCH` and exited 0.** `family_grammar.py` detected the
   contradiction between its repaired targets and its unrepaired derivation, said
   so on stdout, and returned success. It is now **fatal on mismatch**
   (mutation-tested).
2. **A checker kept its own private copy of the answer.** `family_grammar_verify.py`
   held a second `LANDED` table, still stale, and so verified a stale copy against a
   stale derivation — passing 210/210. New check `A9` **cross-checks the two
   independent tables against each other.** Independence between a module and its
   checker is worth nothing if both keep private copies and never compare them.
3. **Gated checks described the wrong case self-consistently.** `minimal_core.py`
   F0/F3 computed with `(201, 36, 504)` and `carry(6,30,67,12)` under a label
   reading "(75,125)". The arithmetic was internally consistent, so nothing failed.
   Now pinned to the corrected `(80,29)`.

The general lesson, stated because it cost three separate defects: **a retrodiction
test confirms that a tool still finds what you already knew, and is structurally
incapable of noticing capability you have lost.** Every fix that stuck was a
*cross-check between two things that should agree*, not another assertion about one
thing.

## New in this release

- **`PROOF_72_108.md`** — the (72,108) proof written out (~40 pages including
  setup; the mathematical core is 15-20). All five `C0` subcases are
  `exact-checked`, and **no step is irreducibly machine work**:
  `support_certificates.py` retired the last one, reducing the marked-support test
  over all 40 `(k,z)` pairs to five Bezout identities of degree at most 6 plus
  degree and valuation bookkeeping — no Groebner engine, no irreducibility, no
  field theory.
- **`CORNER_ATLAS.md` + `corner_atlas.py`** — GGV5's 34 published candidate
  counterexamples run through five gates. They collapse to **six signatures**, and
  **28 of the 34 are driven by one integer test**: retraction failure at `A_0` makes
  `C` a monomial, which simultaneously refuses the chart dictionary, forces
  `lambda = 0` (killing the slice cascade), and empties the Belyi sweep. As of this
  release it is *four* mechanisms — the F2 closed form is the fourth casualty of the
  same integer fact.
  **No row in the atlas is eliminated as a counterexample by the atlas.** The gates
  are on *our mechanisms*, not on the cases.
- **`gamma_from_corner.py`** — derives the chart exponent gamma from the corner
  `A_0` instead of reading it off a paper, calibrated on **28 published data points**
  across GGV1 Tables 1-3. **Result: the corner does not pin gamma.** At `(5,20)`,
  gamma is in `{2,3,4}`; GGV3 section 5 asserts `{2,3}` and says outright *"We do
  not provide proofs for this first part."* Both halves of GGV1's dichotomy fail to
  separate them: Prop 'case II' is **inapplicable** (`gcd(a,b) = 1` for all three)
  and Prop 'extremosfinales' is **non-exclusive** (an admissible `k` exists for all
  three). This is a **located obligation, not a refutation** — GGV1 carries
  machinery not mechanised here, and check `E1` keeps a standing witness of that
  incompleteness. The gap is **class-wide** across the nine atlas rows with
  `b0 = 4*a0`.
- **`CONTACT_LEMMA.md`** — the cascade is a theorem with three previously hidden
  hypotheses (`gcd(m,n)=1`, `lambda >= m`, `N_Q >= D_P+D_Q`); the general profile
  bound is `v_t(h_k) >= m*k-1`, not `2k-1`, and the cokernel is `m*k-n`, not `2n-3`.
- **`TORIC_SYZYGY.md` / `TORIC_GENERAL.md`** — `6WZ = e^5` on every bare G-point,
  using no `Phi`, no caps and no slices; and the toric exponent condition
  `(t+1) | (4t+9)` if and only if `t = 4`, uniquely.
- **`PASSPORT_75_125_REPAIR.md`** — the full 34-file account of the chart repair.

## NEW RESULT — the bridge identity, and monomiality is TWO defects

`MONOMIAL_WINDOW_LAW.md` + `monomial_window_law.py` (56/56).

**An exact identity linking the analytic and combinatorial window invariants.**
From `rho = q(b-a)+1` and `N = a*M - 2b`:

```
ord_y(Phi)  =  a*q*M  -  H,        H := q(a+b) - 1
```

Exact, not a congruence. Hence `gcd(M, ord_y Phi) = gcd(M, H)`, so the **analytic**
window denominator `denom(ord_y Phi / M)` — the quantity the carry obstruction
actually consumes — **is** the **combinatorial** corner invariant
`q_window = M/gcd(M,H)`. These were previously two objects that happened to agree
on four tabulated cases.

This dissolves a flagged negative. `MINIMAL_CORE.md` §4 recorded a "negative I
could not overcome": no family-wide `q_window` sweep was possible because
`ord_y Phi` is published only at `(8,28)`. The identity supplies it from corner
data alone — at `(72,108)`, `2*7*17 - 34 = 204`, which the disjoint `f1`-ODE route
independently confirms. The atlas now carries `ord_y Phi` for all 34 rows.

**Monomiality is two orthogonal defects, not one bias.** This corrects the framing
that motivated the lane. All four previously documented deaths (`lambda`, the
window cone, the toric identity, the F2 closed form) consume **thinness**,
`deg C - ord C = 0`. The `q_window` behaviour consumes **shallowness**,
`ord C = 1`, and has *zero* `deg C` dependence. All four quadrants are realised:
`C = y^8` at `t=3` is a **monomial with an integral window** (`q_window = 1`), and
`C = y(y+1)` has `lambda != 0` with `q_window = M`. So restoring a residual is not
a repair for `q_window`, deepening `C` is not a repair for `lambda`, and two
independent repairs are needed rather than one.

**Monomial rigidity (PROVED).** At `q = ord_y C = 1` with `kappa = t-2`, the Bezout
corner integer is exactly `-1`, so `gcd(M,H) = 1` and `q_window = M` — maximal, at
every corner and every rung. The total-carry lemma then applies: for every split of
`M` and every `k`, the carry lies in `[1, k-1]` and is **never 0**. So the
Phi-divisor carry obstruction is **total at all 28 monomial rows**, decided with no
split enumeration and no `w(e)` datum. This resolves **31 UNKNOWN verdicts at atlas
gate `G5`**, leaving an escape at exactly the 6 non-monomial rows. Corollary:
`q_window = 1` is *impossible* at a monomial corner, since it requires
`ord_y C >= t`.

**Consequence — a lane closes.** The Phi-divisor / K-syzygy route is now provably
dead **class-wide** under the unchanged extreme-ray premise, not merely at
`(75,125)`. The K-syzygy should not be ported to the class of nine. Reviving that
route requires `ord_y C >= t`, which only the retraction shape supplies — so the
route is coextensive with the six retracting rows, one of which is already closed.

**The discriminating pair.** `(50,75)` and `(72,108)` share `(a,b,t) = (2,3,4)` and
`M = 17`, and `build_gsystem` is **field-by-field identical** on both — generators,
`Klin`, state and spare inventory, u-weights — with only `q`, `ord_y Phi` and
`W_step` differing. So the K-syzygy exists as an *algebraic* relation at `(50,75)`
too; the carry on the published split `(5,12)` is `0` versus `1`. **At `(50,75)` the
algebra permits the syzygy and the arithmetic forbids it, and exactly one integer
moved.** Quantitatively, at the `(4,2,5)` shape shared by `(72,108)` and **eight of
the nine** class rows, integrality requires `ord_y C = 7 mod 17`, whose minimal
solution is `7` — precisely `(72,108)`'s value.

**What the positive is, and what it is not.** `ceil(alpha*w/q_window)` *is* the
`ENDPOINT_CONTRACT` depth-ledger floor. Its gain over the affine ray sums to
exactly `(q_window-1)/2` per period, identically `0` **iff** `q_window = 1`. So
`(72,108)`, where the entire toolkit was calibrated, has gain `0` everywhere, while
the class of nine has the *maximal* gain: **the cone needs two slopes and dies; the
floor needs one and is strengthened.** But that the raised floor actually *fires* a
kill at `(8,32)`, `(9,36)`, `(10,40)` is **CLAIMED, NOT CHECKED** — it needs each
corner's reduced gamma-chart caps, which are in-repo for `(50,75)` only.

**Soft spot, flagged rather than buried.** The generality of `rho = q(b-a)+1` and
`N = a*M - 2b` *off* the F2 family is **INFERRED**. They are confirmed at
`(72,108)` (`q=7`, `b-a=1`) and across F2 rungs `a=2..8` (`q=1`, `b-a` varying) — a
one-dimensional slice in each variable separately, with joint dependence untested.
This is the one place the result could be wrong in a way that matters. The cheapest
decisive test: derive `ord_y Phi` independently at `(8,28)/(3,4)/144`, where the
identity predicts **205**.

## Status of the main claim, stated precisely

**(72,108) is closed**, the enumerated `f31` frontier is empty, and every kill in
the chain is independently audited. `C0` itself is recorded at level **`claimed`**,
and that is **correct by construction, not a backlog item**: exhaustiveness of the
case partition rests on GGHV22 Prop 4.3 and the field-split framework — published
mathematics that no finite bookkeeping checker can re-derive. `prop43_audit.py`
discharges the *citation*, not the *partition*. The routes above `claimed` are a
machine-checkable reformulation of the partition, or a formal proof; not a regrade.

**Prior work.** M. Helali published an independent treatment of this case on
**2026-07-21**, ahead of this work. Our adjudication of the two results is
`SUBSUMES`.

## Open problems

1. **Settle the bridge identity's generality** — derive `ord_y Phi` independently at
   `(8,28)/(3,4)/144`, where the identity predicts `205`. Cheapest high-value item
   left; it turns the identity from INFERRED to checked in two directions.
2. **Cash the eight-row G-system coincidence.** Any *weight-free algebraic*
   consequence of `(72,108)`'s ideal transfers verbatim to eight of the nine class
   rows; only weight-normalised steps fail. This wasn't visible before and is the
   strongest surviving lead. It calls for re-auditing the `(72,108)` closure to
   separate its weight-free steps from its weight-normalised ones.
3. **Test the floor-raising claim** by deriving the gamma-chart caps at `(8,32)` —
   noting that `(8,32)`'s window arithmetic is *identical* to `(50,75)`'s, so the
   window layer cannot be what distinguishes them, and whatever makes
   `gamma_from_corner.py` find no branch at `(8,32)` must come from elsewhere.
4. Exclude `gamma = 4` at `(5,20)` — class-wide over nine rows. Needs the reduction
   layer; the corner layer is proven insufficient.
5. Derive GGV3 section 5's conditions (a1)-(a6) from corner data — the gamma-window
   compiler's next step. Must be built for `gamma` in `{2,3,4}` until (4) closes.
6. What replaces the withdrawn closed form when `dg = 0`, for F2 and F3 alike.
7. Repair the eight downstream modules that still transcribe `F3`'s superseded chart.

---

# v0.3.0 — the audit wave: a clean trust graph, a new kill mechanism, two errata

**Headline:** the machine-honesty layer that shipped "working but finding
problems" in the last release now exits clean, and a new *cheap* kill mechanism
lands 133 previously-unrecorded state kills. Both errata below were found by
adversarial review of our own work, not by a checker passing.

- **Proof-DAG audit round v2 — 0 inconsistencies (`PROOF_DAG.md`).** The three
  real doc-vs-data inconsistencies the DAG found on day one are closed, and
  `proof_dag_report.py --quiet` is now a suite gate (exit 0 iff no doc claim
  exceeds what the DAG supports). The cascade branch audits are machine-joined:
  **2289 of 2401 engine-killed branches are now `independently-audited`**
  (sub2 390, sub1 1899); the residual 112 are t/inf-layer kills outside the
  q-cascade auditors' scope and honestly stay `claimed`. One orphan certificate
  turned out to be a real mapping bug (an `a8` resolver hard-coded branch T2 for
  what is a T1 signature) — fixed, and that state is now `certified`.

- **The q_window theorem (`q_window_theorem.py`) — and a surprise.** The identity
  `t·H − q·M = q(κ+1) − t` is proved symbolically, the divisibility lemma
  `gcd(M,H) | C` verified on all 3995 census rows, and the census re-run with the
  exact formula rather than the old proxy. **(72,108) is NOT the unique integral
  case:** there are **51 integral (`q_window = 1`) families across 23 corner
  shapes** — an arithmetic lattice of integral windows that the proxy census
  could not see.

- **Period-12 window functions for (75,125) (`WINDOW_FUNCTIONS_75_125.md`).**
  Lower cap `L(w) = ⌈67w/12⌉` is quasipolynomial (`α = 67 = 10a²−8a+1`,
  `β_m = 5m mod 12`); the upper cap `U(w) = 14w` is **affine** — only the lower
  cap is quasipolynomial, which is the exact content of "quasipolynomial window
  cap". Controls: the `a=2` substitution reproduces the (50,75) window table
  exactly, and the (72,108) limit degenerates to the known affine caps.

- **NEW — the Φ-window-depth kill criterion (`FACE_KILL_SWEEP.md`).** `G5` is
  u-weight-homogeneous of weight 17, so every term strips by exactly `y²⁰⁴` and
  the stripped `G5` lives at y-degree ≤ 34; the stripped `Φ` has degree *exactly*
  34 with `lc = −1024/3315 ≠ 0`. So whenever every `G5body` term stays below
  degree 34 — even with all spares at their caps — the degree-34 coefficient is
  `lc(Φ)` and the equation reads `−1024/3315 = 0`: **the state admits no spares.**
  Being a pure degree argument it needs no symbolic expansion, giving an
  O(1)-per-state test that sweeps all 52,005 states in seconds:
  **sub2 195/7888, sub1 0/44117**. 62 were already in the ledger; **133 are new**,
  and all 133 lie in still-open branches (branch-join discharged: every phase-D
  case maps to a cascade branch marked `survives`). The limits are structural:
  `a_t = 10` and *all* of sub1 are provably immune. Verifier
  `phi_depth_criterion_verify.py` (36 checks) re-derives everything from
  primitives and ties the closed form to the exact symbolic computation.
  These 133 are same-author `claimed`, pending independent audit.

- **ERRATUM 1 — `G5` normalisation.** The canonical generator is
  `G5 = G5body + Φ` (`full_system_bridge.py`; the C11 membership certificate in
  `f37_sat_verify.py` verifies `f31 == … + c4·(G5body + Φ)`).
  `FULL_SYSTEM_BRIDGE.md` line 62 stated `2·Φ + G5body`, contradicting line 114
  of the same file, and `bigrade_annotator.py` transcribed it. The two differ by
  `Φ`, **not** by a nonzero scalar, so conclusions do not transfer automatically.
  Both corrected with dated notes; no landed kill changes (the depth argument
  uses only `deg Φ = 34` and `lc ≠ 0`), but emitted certificate *values* were
  wrong. A fail-loud guard now prevents reintroduction.

- **ERRATUM 2 — the "bandwidth-2 staircase" is refuted, and explained.**
  `BIGRADED_PROBE.md` §2 claimed the walled system's spare coupling is a
  bandwidth-2 staircase. Measured on the computable sub2 home G-system
  (`window_band_probe.py`), the per-slice increment is **+3**, and it is not a
  bandwidth at all — it is **the number of spare series** (`dm2,dm3,dm4` →
  `R,S,T`; one new coefficient each per slice). The R9 system has `dm4`
  eliminated, leaving two series, hence its `+2`. The probe measured a real
  number and promoted it to a structural law it is not. What survives: the
  decomposition direction, the exact top-corner equations, the extracted `a4`
  constraint. What does not: every quantitative lattice claim in §2.

- **Endpoint contract (`ENDPOINT_CONTRACT.md`)**: the per-coefficient
  {required-nonzero | forbidden | optional} contract, with (50,75) fully
  instantiated — the historical kill fires at `c_{0,−10}`, and that rediscovery
  is now an automated regression.

- **Frontier v2 (`FRONTIER_V2.md`)**: the misleading "Killed (audited)" column is
  renamed "Killed (exact-checked, same-author)", with a separate proof-DAG
  evidence-grade table. Same-author exact checks are no longer presentable as
  independent audits.

- Suite: **63 checkers**, `clean_clone_check` first (fails loudly if any checker
  reads an untracked file), `proof_dag_report --quiet` as a consistency gate.

**Scope, stated plainly.** C0 remains **OPEN**. The 133 new kills close **0 of
220** sub2 flag-cases — closing a branch requires killing states at *all* degrees
in it, and every mechanism we have is degree-limited. That gap is the honest
description of where this program is stuck.

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
