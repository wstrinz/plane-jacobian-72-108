# Cross-program study & variable map: `eumemic/jacobian-keller-research` (DC1 band method) vs our (72,108) cascade

**Date:** 2026-07-24. **Status: READ-ONLY study, no audited artifact touched.**
New files only: this doc + `eumemic_import_candidates.json`. Nothing committed;
no existing file edited; no concurrent-lane file touched.

The external program is an independent, adversarially-audited public campaign. This
note imports it *only* as a study and a method map — exactly the posture of
`ALOK_CROSSCHECK.md`, but where alok attacked the **same** (72,108) polygons from a
disjoint regime, eumemic attacks a **different problem** (DC1 / Weyl algebra `A_1`)
linked to ours only by `JC2 ⇒ DC1`. We respect their tier labels verbatim.

- **Repo:** <https://github.com/eumemic/jacobian-keller-research>
- **Commit pinned:** `201c2f6530a774bfed10939a90ae97b5cfcf80cf`
  ("Repair DC1 evidence scope and verifier reporting", 2026-07-24), current `main` HEAD.
- **Re-fetch:** `gh api repos/eumemic/jacobian-keller-research/contents/<path> --jq .content | base64 -d`

---

## 1. The band method in three sentences

The **band method** works in the Weyl algebra `A_1 = C⟨x,∂⟩/([∂,x]-1)` localized at
`x`, `A_1[x^{-1}] = ⊕_k x^k C[E]` with `E = x∂`, where an operator's **band** is the
width of its `ad(E)`-weight support (its window of `x`-powers, since
`[E, x^k f(E)] = k·x^k f(E)`). The commutator law `[D,X]=1` decomposes into
**degree-free ladder rung-equations** `Q_m = Σ_{k+l=m}[b_l^{[k]} a_k − a_k^{[l]} b_l] =
δ_{m0}` in the shift algebra (`f^{[n]}(E)=f(E+n)`, `σ=T^{-1}`), which on
root-multisets become **necklace** identities and are probed by linear **covector**
functionals (node-evaluations `ev_ρ` and moving-sum adjoints `S_n^*`, computed
degree-free as **trace forms** `Tr_{F[E]/(p)}`). The whole DC1 program is a
*structural* reduction: force a hypothetical non-generating (counterexample) pair down
to band `≤ 2` (proven tame), where the descent stalls at band-3 **walls** — the
shifted-power wall (closed by a difference-operator Abhyankar–Moh descent) and the
minimal **singular hatch `W2`** (`a_3=E(E+2)(E+4)`, `b_2=E(E+3)`), whose
arbitrary-degree realizability is the current frontier.

**Covectors / necklaces / W2, concretely.**
- **Covector** = a linear functional on the space of ladder coefficients used to
  *annihilate* a block of the `Q_m` equations and thereby extract a scalar
  obstruction. Two kinds: point functionals `ev_ρ : f ↦ f(ρ)` and moving-sum
  adjoints. Galois-symmetric ones over the roots of a datum polynomial `p` are
  **trace forms** — rational in the coefficients of `p`, "computed with no root named"
  (`algebraic-covector.md` §1). This is their degree-free engine.
- **Necklace** = the root-multiset representation `δ(p)` of a coefficient polynomial,
  with `σ = T^{-1}` acting as a shift; the top wall becomes the necklace identity
  `(σ^k−1)·δ(b_q) = (σ^q−1)·δ(a_k)`, solved by a cofactor `g ∈ Z[σ]`
  (`band-reduction.md` §5).
- **W2 obstruction** = the minimal band-3 **singular hatch**: `a_3=E(E+2)(E+4)`,
  `b_2=E(E+3)`, the smallest Dixmier-stuck leading form `(a,b)=(3,2)`, `p=x^2ξ`
  (`band-reduction.md` §4, `w2-decisive.md`). "W2 obstruction localization" = the
  analysis of whether this hatch carries a genuine Weyl pair; the negative tail forces
  a residual factor `W=0`, giving a **bounded** kill at `d=3` (`w2-verdict.md`).

## 2. What they have PROVEN vs conjectured (their own tier labels)

Reported using their labels; where their docs run ahead of their proofs I flag it
neutrally, as they themselves do.

**Peer-reviewed / classical inputs they build on (proved elsewhere):**
- Classical band-2 **M4 non-square** theorem — proved, alg-closed char 0
  (`band2-classical-proved/`).
- Dixmier's 1968 leading-symbol lemma (they use only this; take DC1 as **open**,
  noting Zheglov arXiv:2410.06959 *claims* DC1 but a Jan-2026 survey still lists it
  open — `band-reduction.md` §3 literature note).

**Independently assembled + internally audited, NOT peer reviewed (their headline tier):**
- **Full classical band-2 theorem** (`band2-classical-full/`): every classical Keller
  pair with both ladder supports in `[-2,2]` is a polynomial automorphism.
- **Full fixed quantum band-2 theorem** (`band2-square-sector/`): every genuine `A_1`
  pair with both `E`-ladder supports in `[-2,2]` is a tame automorphism image.
  *Their frontier "furthest result." Explicitly: fixed-band restriction essential,
  "does not claim DC1."*
- **Audited J3 band-2 obstruction** (`band2-j3-provisional/`): core Laurent/localized
  theorem, independently audited & repaired.

**Independently derived, machine-checked identities, mixed proved/bounded/open (the DC1 frontier):**
- **Band-reduction framework** (`dc1-program/band-reduction.md`): PROVED components —
  floors (non-generating ⟹ band `k≥3`), Dixmier stratum (`a∤b, b∤a`, minimal `(2,3)`),
  unconditional top wall `W(k,q)`, shifted-power/singular-hatch effectivity dichotomy.
  **CONDITIONAL** conclusion — two named open gaps (§9): composite-move escape (= the
  classical DC1 core) and shifted-power descent. Honest self-label: "classifies the
  fixed points of single-transvection reduction, not proven to be all minima."
- **Shifted-power descent** (`dc1-program/shifted-power-descent.md`): Gap 2 resolved at
  band 3 for *cube-separated* `h`, arbitrary degree — but the key propagation step
  `Q_0=1 ⇒ h` constant is **audit-flagged CONDITIONAL** on the full sub-leading ansatz
  including `h|a_1` (top layer machine-forced, `h|a_1` imported not derived; a later
  commit "h|a_1 discharged" claims to close this — `d5e94aca3`). The reduction *move*
  is audit-demoted to a single-instance check.
- **Degree-free slope forcing** (`dc1-program/slope-forcing-degree-free.md`): the
  factorization `R(1) = a_2(0)·W` is **PROVED arbitrary degree**; the residual
  tail-forcing identity `a_2(0)·W ∈ √(cascade+Q_{-1..-3})` is **machine-confirmed at
  d=3 only** (d=4 GB not tractable), OPEN at arbitrary degree.
- **W2 bounded verdict** (`band3/w2-verdict.md`): the combined slope+tail system is the
  **unit ideal at raw cap `d=3`, both branches**, via a committed **exact `QQ`
  Nullstellensatz multiplier certificate** (`Σ h_i f_i = 1`, re-checked without
  trusting the solver). `d=4` unit is documentary/optional-msolve; `d=5 mod 65003`
  documentary. Explicitly **not** an arbitrary-degree theorem.
- **Algebraic-node covector calculus** (`dc1-program/algebraic-covector.md`,
  `joint-covector.md`): the trace-form tool is PROVED node-free; the W2 residual is
  **localized** to a "varying-tops coupling" obstruction — a fixed finite trace-form
  recipe producing `W` at *every* `d` is **NOT obtained** (OPEN).

**Explicitly conjectured / disclaimed:**
- **DC1 and JC2 themselves** — "No Weyl pair and no counterexample is constructed"
  appears in every ledger. Their **direct JC2 attack memo** (`jc2-attack-memo.md`) is a
  provisional reductions/reconnaissance memo with *no claimed result*; its milestone
  (affine-quotient obstruction for the 3D telescope family) is posed, "no outcome
  presently established."

**Current frontier:** W2 singular-hatch realizability at arbitrary degree — proving the
residual slope-forcing identity via a **degree-free covector for the `(a_2,b_1)`
necklace block**. Their bounded machinery kills W2 at `d=3` (and `d=3–4` jointly);
the arbitrary-degree closure is the open endgame. This is the exact structural analogue
of *our* f31-survivor endgame (§3).

## 3. The variable / method map

**Framing.** Their objects and ours live in different rings and solve different
conjectures. What is genuinely shared is the **method schema**, and it is shared
deeply. Both programs: (i) take a commutator law, (ii) grade it by a
valuation/weight, (iii) extract per-rung elimination equations, (iv) prove a
**degree-free "slope law,"** and (v) bottom out in a **bounded-degree exact UNIT-ideal
kill** plus an **open arbitrary-degree residual needing a degree-free covector.**

| role in the schema | eumemic (DC1 / `A_1`) | ours ((72,108) cascade) |
|---|---|---|
| ambient graded ring | `A_1[x^{-1}] = ⊕_k x^k C[E]`, `E=x∂` | localized `K[y, C4^{-1}]`, `C4=y^7(y+1)` |
| the grading / "band" ⟷ "window" | band `k` = `ad(E)`-weight width | window index `k`; `d_k := c_k·C4^{7-2k}` |
| commutator / defining law | `[D,X]=1` | `[P,Q]=x^2`; reduced `(D̃^3)_{-5}+Φ=0` |
| per-rung elimination equations | ladder `Q_m = Σ[b_l^{[k]}a_k − a_k^{[l]}b_l]` | G-system `G1..G5 = (D̃^3)_{-1,-2,-3,-5}` |
| shift / spare unknowns | `f^{[n]}=f(E+n)`, `σ=T^{-1}` | spare window vars `d_{-2},d_{-3},d_{-4}` |
| forcing nonzero term | moment unit `Q_0=1` (`R(1)=1`) | `Φ = c·t^30·q ≠ 0` (`c=−1/6630`) |
| the "slope law" (degree-free, PROVEN) | shifted-power wall `b_{k-1}=κ∏h^{[j]}`; `R(1)=a_2(0)·W` | corner law `κ = t−2` (chart Jacobian `−x^{l-2}`); `deg Φ = (e·a0−q+1)+N·a0` |
| bounded exact UNIT-ideal kill | W2 combined system UNIT at `d=3` (`QQ` Nullstellensatz certificate) | FULL_SYSTEM_BRIDGE: full G-system exact UNIT over `QQ` (a8 pilot) |
| open arbitrary-degree residual | degree-free covector for `(a_2,b_1)` necklace | degree-free residue lemmas / cone lemmas for the D-lattice |
| solver | Singular / `msolve` (`^` graded) | Singular via WSL / `msolve` (J6_MSOLVE) |

### 3(a) Do their slope-forcing lemmas constrain the same objects as our corner signature `(t, κ, a0, q)`? — ANALOGOUS, NOT THE SAME OBJECT

Our corner signature (`STATE.md` 2026-07-23 entries; `PHI_CORNER4.md`): `t = l` the
Laurent-chart slope (`(X,Y)=(x^{-1}, x^l y)`), `κ` with **`κ = t−2` PROVEN structural**
from the chart Jacobian `−x^{l-2}`, `a0 = deg C`, and `(e,q,a,b,r,N)` the Phi
degree/mult data.

Their slope-forcing (`slope-forcing-degree-free.md`, `shifted-power-descent.md`): `κ` is
the **shifted-power wall constant** (ratio of the two top-band leading coefficients,
`b_{k-1}=κ·∏h^{[j]}`); the "slope" `R` is the **moment slope** (`E−R | D`). Their
proven degree-free content is `R(1)=a_2(0)·W`.

**Verdict: COMPLEMENTARY / analogous, not a variable identification.** Both `κ`'s are
"slope data of a graded object satisfying a *proven degree-free* structural law," and
both descents are difference/Newton Abhyankar–Moh steps — but their `κ` is a Weyl-ladder
leading-coefficient ratio in `C[E]`, ours a Newton-polygon chart-slope in the plane. The
corner law **cannot import into their framework as a constraint on their `κ`** (different
ring), and vice versa. What *is* transferable is the **proof pattern** — "derive the
slope from a change-of-chart Jacobian rather than from ladder rungs" (our closed
`κ=t−2`) is a fully-proven instance of exactly the kind of arbitrary-degree slope law
their `shifted-power-descent` §3 and `slope-forcing` §6 leave open. See exchange #3/#4.

### 3(b) Does their W2 localization touch our (72,108) polygons? — DISJOINT

Their W2 is a specific band-3 Weyl datum; our (72,108) live branches are dense-regime
D-lattice degree-states carrying `Φ≠0`. No dictionary maps `a_3=E(E+2)(E+4)` to a
(72,108) Newton corner — they are objects in different problems. Their W2 kill (`d=3`
unit ideal) touches **none** of our branches; our 390 sub2 + 2007 sub1 engine kills
imply **nothing** about W2 realizability, and W2 implies nothing about ours.

Logic check (the useful direction): `JC2 ⇒ DC1` gives `¬DC1 ⇒ ¬JC2`, but `¬JC2` does
**not** imply `¬DC1`. So even a genuine (72,108) counterexample would not yield a DC1
counterexample, and their DC1 band-fragment proofs do not constrain (72,108). The only
shared logical node is the "no counterexample anywhere" world; **no theorem transfers in
either direction.** Verdict: **DISJOINT** (a stronger disjointness than alok's, which at
least shared the (72,108) polygons).

### 3(c) Degree-free vs our fixed-degree ledger — complementary like alok, but one level up

Alok was regime-disjoint *within* (72,108) (sparse floor vs dense body). Eumemic is
**problem-disjoint** (DC1 vs plane (72,108)), linked only by `JC2⇒DC1`. But their
degree-free method is **complementary to our fixed-degree ledger as a methodology** — it
is the *same* internal tension we already have: our WINDOW_CAPS / RESIDUE_LEMMAS / cone
lemmas are our nascent degree-free layer, our cascade/bridge engine is the fixed-degree
layer. Their covector/trace-form calculus is a more-developed instance of the degree-free
layer we are building. **Verdict: methodologically complementary and mutually
instructive; object-disjoint.**

## 4. Top exchange opportunities (precise, with translation cost)

**#1 — IMPORT as a checker (alok-style, lowest cost): committed Nullstellensatz-multiplier
certificate format + the `msolve ^`-vs-`**` discipline.**
Their `w2-verdict.md` §2 commits sparse multipliers `h_i` with `Σ h_i f_i = 1` checked by
exact multivariate coefficient collection, so the verifier certifies UNIT **without
trusting the solver** (`w2_d3_qq_certificates.txt`). Our `FULL_SYSTEM_BRIDGE.md` §5
reports "exact UNIT over `QQ` in 8.75 s" but leans on Singular's verdict. Adopting their
committed-multiplier format would raise our bridge kills to their audit tier. Bundled:
their explicit finding that `msolve`'s `**` parser fabricates a spurious "surviving
sub-locus" that `^` dispels (`joint-covector.md` §4, `hatch-census.md` §5) — a concrete
hygiene fix for our `msolve` usage (`J6_MSOLVE.md`, `MSOLVE_REPORT.md`).
*Translation cost: minimal — certificate-format + tool-flag discipline, directly portable.*

**#2 — IMPORT as a method template: the degree-free covector / trace-form calculus for
lifting our bounded bridge kills to arbitrary degree.**
Their algebraic-node Thm A′ + trace-form descent (`Tr_{F[E]/(p)}` computes
Galois-symmetric annihilators over algebraic node-sets *with no root named*,
`algebraic-covector.md` §1) is precisely the tool our program lacks to convert the
FULL_SYSTEM_BRIDGE's per-state `d`-bounded UNIT kills into a degree-free closure. Our
RESIDUE_LEMMAS / cone lemmas are our version; theirs is further along.
*Translation cost: substantial — our G-system is weighted-homogeneous in `K[y,C4^{-1}]`,
theirs shift-graded in `C[E]`; the covector-annihilation idea transfers but the ring and
the "necklace" (Newton-corner vs root-multiset) must be re-set-up. Honest boundary: their
own arbitrary-`d` recipe is OPEN, so this is a shared *research direction*, not a
finished import.*

**#3 — SHARED OBJECT worth a joint definition: the "graded commutator ladder →
slope-forcing → bounded UNIT ideal → degree-free residual" schema.**
Both programs independently arrived at the same five-step schema (§3 table) and even the
same endgame shape (bounded exact UNIT kill; open degree-free residual). A joint glossary
— their `Q_m` ↔ our `(D̃^3)_{-l}` / G-system; their band `k` ↔ our window `k`; their
necklace ↔ our Newton corner; their `κ` shifted-power wall ↔ our `κ=t−2` chart law; their
moment-unit `Q_0=1` ↔ our `Φ≠0` — would let each program cite the other's degree-free
closure lemma once proven, and sanity-check the shared "unit-ideal-at-`d`" tooling.
*Translation cost: definitional, not mathematical — it formalizes an analogy both sides
already feel.*

**#4 — EXPORT from us (offer, honest boundary): our proven `κ = t−2` slope law and the
arbitrary-`k` WINDOW_CAPS induction as a template for their open slope/propagation gaps.**
They have an OPEN "shifted-power descent propagation" (`band-reduction.md` Gap 2,
conditional on `h|a_1`) and OPEN arbitrary-degree slope forcing. Our `κ=t−2` — proved
degree-free from a change-of-chart Jacobian for generic `l` (`STATE.md` 2026-07-23 Lane
D) — and our now-PROVEN `k=6,7,8` window caps (`WINDOW_CAPS.md`) are fully-closed
instances of exactly the "slope forced degree-free" statements they are chasing.
*Translation cost: this is a technique/analogy export (different ring), **not** a theorem
that plugs into their verifier — offer it as a proof strategy, not a result.*

## 5. Honest-boundary summary

- We respect their tiers: their **headline** proved result is the fixed-band quantum
  band-2 theorem (internally audited, not peer reviewed); DC1 and JC2 are **open and
  disclaimed**; the band-reduction framework is explicitly **conditional** on the DC1
  core; W2 closure is **bounded to `d=3`** with `d≥4` documentary. Their audit discipline
  (audit-demotions, `SKIP`-not-folded-into-PASS, solver-parser caveats) is exemplary and
  is itself part of what we import (#1).
- Their docs are, in a few places, slightly ahead of the committed proofs — the
  `shifted-power-descent` `h|a_1` conditionality, the `d=4/d=5` W2 documentary rows, the
  "varying-tops coupling" open residual — and in every case *they say so themselves*. We
  report their frontier as they label it.
- **Map verdict:** 3(a) **complementary/analogous** (not importable as a direct
  constraint); 3(b) **disjoint** (no theorem transfers to/from (72,108)); 3(c)
  **methodologically complementary, object-disjoint**. The value of this program to ours
  is **method and audit-discipline transfer**, not a branch kill.

## 6. Files created (new only)

- `EUMEMIC_MAP.md` — this study + map + exchange opportunities.
- `eumemic_import_candidates.json` — machine-readable importable-claim list with their
  tier labels and our translation status.

No existing file edited; no commit; no concurrent-lane file touched.
