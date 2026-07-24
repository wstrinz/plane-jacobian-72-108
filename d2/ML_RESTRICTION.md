# Makar-Limanov quasi-homogeneous edge-form restriction vs. GGHV (72,108)

**Lane:** literature-comparison. **Date:** 2026-07-24. **Status:** landed, exit-0.
**New files (uncommitted):** `ML_RESTRICTION.md`, `ml_restriction_check.py`.
Read-only over every existing artifact; no concurrent-lane files touched.

**One-line verdict.** The Makar-Limanov restriction is **INAPPLICABLE to every
edge pair of both (8,28)-reduced polygons and to all 17 survey families** — its
defining hypothesis `J(rho,tau)=rho` with `rho` nonmonomial and *positive*
weights is never realized in our data, which sits entirely in the complementary
`J(edge forms) in {0, monomial}` regime. **Nothing in the frontier is affected.**

---

## 1. The paper

**L. Makar-Limanov, "On the shape of a counterexample to the two-dimensional
Jacobian conjecture," _Serdica Mathematical Journal_ **51** (2025), no. 3–4,
pp. 299–314. DOI `10.55630/serdica.2025.51.299-314`.** (Published 2025-12-23;
Wayne State University.)

Abstract (verbatim from the journal record): *"Polynomials f, g ∈ ℂ[x,y] is a
Jacobian pair if the Jacobian J(f,g) = 1. The Jacobian conjecture (JC) formulated
by O. H. Keller states that then ℂ[f,g] = ℂ[x,y]. In this paper further
information on the shape of the Newton polygon of f if the pair f, g is a
counterexample to JC is obtained."*

### Sourcing honesty contract (read this before trusting the lemma text below)

- **The full text is closed-access.** No open galley, no arXiv preprint, no
  author-hosted copy was reachable (checked: Serdica OJS galley URLs → 404; DOI
  resolves only to the abstract page; arXiv author listing / API unreachable in
  this environment; the paper carries no arXiv id). **I did not read the paper
  body, so I cannot quote the lemma verbatim, cannot give its lemma/theorem
  number, and cannot confirm its exact in-paper hypotheses.**
- The lemma tested here is the **external review's paraphrase**, which I have (a)
  made mathematically precise, and (b) independently **verified and corrected**
  (§2). Per the task's instruction to *report what it actually says and not force
  the match*: the honest status is **"review-paraphrased lemma, self-certified as
  a true statement under the corrected hypotheses, applied conservatively."** If
  the program later obtains the Serdica PDF, the one open item is to reconcile the
  verbatim hypotheses and lemma number with §2.
- **Provenance that IS solid:** the paper, author, year, volume, page range, and
  DOI above are from the journal's own metadata; the surrounding ML machinery
  (algebraically-dependent homogeneous forms are proportional to a rational power;
  the `rho+sigma = w(x)+w(y)` degree balance forcing `J=1`) is confirmed verbatim
  in the freely-available companion **L. Makar-Limanov, "On the Newton polytope of
  a Jacobian pair," arXiv:2106.06869 (2021)**, lines 272–278 — the same author's
  framework the 2025 lemma extends.

---

## 2. The lemma, made precise and verified (`ml_restriction_check.py` PART A)

**Review paraphrase (as delivered):** *for positive weights and homogeneous ρ,τ
satisfying J(ρ,τ)=ρ, if ρ is nonmonomial then w(ρ) does not divide w(τ).*

**Made precise.** Fix weights `w(x)=α`, `w(y)=β`. For `w`-homogeneous `ρ,τ` with
`J(ρ,τ) := ρ_x τ_y − ρ_y τ_x = ρ`, weighted-degree bookkeeping forces
`w(τ) = α+β` (both sides must be `w`-homogeneous of degree `w(ρ)`, and
`deg_w J(ρ,τ) = w(ρ)+w(τ)−α−β`). The claimed conclusion: **`ρ` nonmonomial ⟹
`w(ρ) ∤ w(τ)` (equivalently `w(ρ) ∤ (α+β)`).**

**Verification + a required hypothesis correction.** I brute-forced the existence
of a nonmonomial `ρ` with `J(ρ,τ)=ρ` over all `w(ρ) | (α+β)` for weights up to 6:

- For **primitive, non-diagonal** weights (`gcd(α,β)=1` and `α≠β`) the lemma
  **holds** with no exception in range (12 nontrivial cases; each has *no*
  nonmonomial `J=ρ` solution).
- The paraphrase **as literally stated is false** on two excluded fringes, and the
  checker records them so the gap is on the record:
  - **the diagonal `α=β`** (e.g. `α=β=1`): `ρ=x+y`, `τ=xy+y²` give
    `J(x+y, xy+y²)=x+y` with `w(ρ)=1 | w(τ)=2` — a genuine counterexample;
  - **non-primitive weights** (`gcd(α,β)>1`), which are just a rescaling of the
    diagonal/degenerate cases.

So the operative, verified hypotheses are **POSITIVE + PRIMITIVE (`gcd=1`) +
NON-DIAGONAL (`α≠β`)** — i.e. a genuine quasi-homogeneous *edge* direction, which
is exactly the regime a Newton-polygon shape lemma is about. This matches the
standard quasi-homogeneous convention ("weights with no common factor",
arXiv:2106.06869 line 135). The checker's PART A self-certifies this every run; a
regression that broke the lemma under its own hypotheses would trip the exit code.

---

## 3. Translation: our edge data in the ML weight convention

Ground truth (`paper_src/upstream_facts.json`, re-checked live by the checker):

- Reduced case **(8,28)**, bracket **`[P,Q] = x²`** (so `J(P,Q)` is the *monomial*
  `x²`).
- **sub1:** `P = {(0,0),(1,0),(8,14),(8,16),(0,8)}`,
  `Q = {(0,0),(2,1),(12,21),(12,24),(0,12)}`.
- **sub2:** `P = {(0,0),(1,0),(8,14),(8,16)}`,
  `Q = {(0,0),(2,1),(12,21),(12,24)}`.
- Principal common-root edge (weight `(1,0)`): `ℓ(P)=R²`, `ℓ(Q)=R³` with
  `R = x⁴ y⁷ (y+1)` (`C4 = y⁷(y+1)`, `FULL_SYSTEM_BRIDGE.md §1`).

**The decisive geometric fact (checker PART B).** Every boundary edge of both
reduced `P`-polygons has an outward normal with a **non-positive component**:

| subcase | edge (P) | primitive normal `w` | positive? | `ℓ_w(P)` |
|---|---|---|---|---|
| sub1/2 | `(0,0)-(1,0)` | `(0,-1)` | no (axis) | nonmono |
| sub1/2 | `(1,0)-(8,14)` | `(2,-1)` | no | nonmono |
| sub1/2 | `(8,14)-(8,16)` | `(1,0)` | no (axis) | nonmono (= R²) |
| sub1 | `(8,16)-(0,8)` | `(-1,1)` | no | nonmono |
| sub1 | `(0,8)-(0,0)` | `(-1,0)` | no (axis) | nonmono |
| sub2 | `(8,16)-(0,0)` | `(-2,1)` | no | nonmono |

Because the origin `(0,0)` is a vertex, **any strictly-positive weight `(α,β)`
attains its maximum at the single vertex `(8,16)`** (resp. `(12,24)` for `Q`) —
a *monomial*, not an edge. Hence:

- **there is no strictly-positive-weight edge at all** in the reduced polygons;
- the only nonmonomial edge forms live on **non-positive (axis / negative-slope)**
  weights, where the ML positivity hypothesis fails outright.

---

## 4. Per-edge verdicts

All 27 tested objects (9 edge pairs + principal edge + 17 families):
**INAPPLICABLE.** Reasons, per class:

**(a) Reduced-polygon edge pairs — 9 pairs, INAPPLICABLE (positivity fails).**
Each sloped/axis edge carries a weight with a non-positive component (§3 table),
so the lemma's "positive weights" hypothesis fails before `ρ`/`J` are even
consulted. There is no positive-weight edge for the lemma to bite.

**(b) Principal common-root edge — INAPPLICABLE (two independent failures).**
Weight `(1,0)` has `w(y)=0` (not positive), *and* the forms are
algebraically dependent: with `R = x⁴y⁷(y+1)`,
```
J(R², R³) = 6 R³ (R_x R_y − R_y R_x) = 0   (computed exactly in the checker).
```
So `J(ρ,τ)=0 ≠ ρ`. Our leading pair sits in the **`J=0`** (proportional-powers)
regime — the case ML handles *separately* (arXiv:2106.06869 line 272: `J=0 ⟹
g_k(e)=cF^λ`) — not the `J=ρ` regime the 2025 lemma governs.

**(c) The 17 length-1 survey families — all INAPPLICABLE.** The reduction chart
`(X,Y) ↦ (x⁻¹, xˡ y)` has Jacobian `−x^{l−2}` (verified symbolically), so in
reduced coordinates `[p,q] = −x^{l−2}` is a **monomial**, and the corner leading
forms of `p,q` are proportional powers of the common root `C = xˡ c` (giving
`J(p_lead,q_lead)=0`). Either way `J(edge forms) ∈ {0, monomial}`, never a
nonmonomial `ρ`. **The lemma does not constrain the family table.**

No edge produced **AUTOMATIC**, **REDUNDANT**, or **REMOVES-A-BRANCH**: the
hypotheses are never even jointly satisfiable on our geometry, so the conclusion
(a restriction on nonmonomial positive-weight edge forms) is vacuous here.

---

## 5. Why there is no inconsistency in either direction

The check is symmetric: exit-1 fires on *either* (i) the lemma being violated by
a live configuration in our ledger (which would mean our data or the lemma is
wrong — a `REMOVES-A-BRANCH`/`INCONSISTENT` finding), or (ii) the lemma failing
its own self-cert. Neither occurs:

- Our entire program lives in the `J(P,Q) = const` (bracket `= x²`, a monomial)
  regime. For any weight `w`, `J(P_w, Q_w)` is either `0` (edge forms share a
  common root, e.g. `R²`/`R³`) or the `w`-leading form of `x²` — always `0` or a
  **monomial**. It can never equal a **nonmonomial** `ρ = P_w`. The ML hypothesis
  `J(ρ,τ)=ρ` with `ρ` nonmonomial is therefore structurally impossible on our
  polygons — the lemma restricts a *different* edge regime than the one the
  (8,28) reduction realizes.
- The lemma is a restriction on the *shape* of a hypothetical minimal
  counterexample's polygon in a positive-weight edge direction; the (8,28)
  reduced polygon has **no positive-weight edge**, so the restriction is silent
  about it. It neither confirms nor kills any branch our ledger tracks.

**Frontier impact: none.** No live case, flag-case, or survey family is removed,
promoted, or contradicted. This is a clean *negative* literature check: the
review's high-leverage candidate does not, in fact, touch the (72,108) program —
useful to record precisely because the review flagged the hypothesis-match as
unchecked.

---

## 6. Files & reproduction

- `ML_RESTRICTION.md` — this note.
- `ml_restriction_check.py` — `--quiet`, **exit 0 iff no inconsistency** in either
  direction. PART A self-certifies the (corrected) lemma; PART B/C/D issue the
  per-edge / per-family verdicts. Self-contained (embeds the polygon data; also
  re-checks it against `paper_src/upstream_facts.json` when present).

```
python ml_restriction_check.py           # full narrative
python ml_restriction_check.py --quiet   # SUMMARY + exit code only
```

Result: `27 verdicts, all INAPPLICABLE; PART A self-cert PASS; exit 0`.

### Open item (single, honest)

Obtain the Serdica 51(3–4) PDF and reconcile the **verbatim** lemma statement +
number + exact hypotheses against §2. The *conclusion* of this lane is robust to
that reconciliation — it rests on the structural fact that our reduced polygons
have no positive-weight edge and sit in the `J∈{0,monomial}` regime, independent
of the lemma's precise wording — but the citation should be upgraded from
"review-paraphrased, self-certified" to "verbatim" once the text is in hand.
