# The plane (72,108) Jacobian case — an independent structural proof

[![Fast checks](https://github.com/wstrinz/plane-jacobian-72-108/actions/workflows/fast-checks.yml/badge.svg)](https://github.com/wstrinz/plane-jacobian-72-108/actions/workflows/fast-checks.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21534895.svg)](https://doi.org/10.5281/zenodo.21534895)

Work on the **Jacobian Conjecture** (JC): *a polynomial map `F: C^n -> C^n` with
nonzero-constant Jacobian determinant is invertible.* This repository holds two
related but separate lines of work.

| dir | dim | what | status |
|-----|-----|------|--------|
| `d3/` | 3 | explicit map refuting JC(3) | **done & verified** |
| `d2/` | 2 | the last open plane case below degree 125 | **excluded — see the conditions below** |

The D3 counterexample settles dimension 3 in the negative. It does **not**
resolve the plane case: JC(2) is a separate open problem, which is what the
D2 effort is about.

---

## Priority: Helali was first

**Billel Helali published an independent exclusion of the (72,108) case on
21 July 2026** — *Exact Computer-Assisted Exclusion of the (72,108) Frontier in
the Two-Dimensional Jacobian Problem*,
[doi:10.5281/zenodo.21479814](https://doi.org/10.5281/zenodo.21479814), with
source and certificates at `bilLkarkariy/jc2-72-108-exact-certificates` commit
`c530fe4`. **He got there first, and his work is the public record for the
exclusion.**

We audited his reduction independently and adjudicate it **SUBSUMES** ours: his
Case 1 covers configuration (1), his Case 2 covers configuration (2), his
parametrisation quantifies over the whole ambient set, and his chain is complete
in the sense our own ledger needs. We enumerated five candidate leak types
(torus normalisation, division by a possibly-vanishing pivot, non-exhaustive
case split, component loss under saturation, incomplete first block) and each is
closed by an explicit check. See `d2/HELALI_ADJUDICATION.md`.

Three further independent claims to the same exclusion appeared within five days
of Helali's (two MathOverflow answers on 07-23, one pull request on 07-25).
**We make no priority claim of any kind.**

What this repository contributes is therefore **not** the exclusion. It is an
independent structural proof of the same exclusion by a different mechanism,
whose intermediate content is finer than a yes/no answer — in particular the
exact pin `a_t = 9` on a t-adic valuation of the germ, and a five-case
classification of what an admissible germ would have to look like.

---

## D2 — the result

**The case.** The only open case below degree 125 for a *plane* JC counterexample
was the GGV–Horruitiner case (8,28) (arXiv:2204.14178, Prop 4.3), left open in
2022 for lack of computing power. GGHV22 state that discharging it *"would
increase the lower bound from 108 up to 125."*

> **Theorem A.** Let `K` be a field of characteristic 0. There is no pair
> `(P,Q)` in `K[x,x^-1,y]^2` with `[P,Q] = x^2` and
> `N(P) = conv{(0,0),(1,0),(8,14),(8,16)}`, `N(Q) = conv{(0,0),(2,1),(12,21),(12,24)}`
> (configuration (2)), nor one with the same corner sets together with
> `(0,8) in N(P)` and `(0,12) in N(Q)` (configuration (1)).
>
> **Corollary B.** Conditional on the two inputs below, any counterexample to the
> plane Jacobian conjecture has `max{deg P, deg Q} >= 125`.

The proof runs on remarkably little: **four polynomials, one polynomial
identity, two slice families, one valuation cascade, and nine integers.** It
factors through one intermediate object, an *admissible germ*, and its skeleton
is: the K-syzygy splits into `e = 0` (dead in one line) and `e | Phi`, which
forces `e = c t^a Pi` with `Pi | q` squarefree; then `a >= 9` and `a <= 9`, so
`a = 9`; then the five cases `deg Pi = 0,1,2,3,4` all die.

**The proof is [`d2/PROOF_72_108.md`](d2/PROOF_72_108.md).** Read that first; it
is self-contained and states its own conditionality up front.

**Of the eleven machine steps, none is irreducibly machine work.** Exactly one
was — the marked-support feasibility test — and `support_certificates.py`
(55/55) reduced it to five Bézout identities in `Z[y]`/`Z[p]` of degree <= 6 plus
degree and valuation bookkeeping, hand-checkable at every one of the 40 `(k,z)`
pairs.

### The two conditions, stated plainly

Corollary B is **conditional**, on exactly two things, and neither is hidden in a
footnote:

1. **GGHV22 Proposition 4.3 must be exhaustive.** Everything begins after that
   reduction. `prop43_audit.py` (20/20) finds the case tree COMPLETE and
   reproduces both Newton configurations vertex-for-vertex against a
   transcription of the source. Its one formerly uncited step is now discharged
   by citation (GGV5 Thm 2.20(8) + Prop 2.5(3), GGV2 Rems 3.8/3.9). This does
   **not** remove the conditionality: the chain consumes GGV5's algorithmically
   generated `A_1` column, and a reader who declines computational input is left
   with `gamma` in `{3,4,5,6,7}`.
2. **The alpha-strip WLOG (`[QQ1]`, Proposition 2.1) is used at the boundary of
   its stated strength.** The lower bound `a_t >= 9` consumes it five slices
   clear of the boundary; the upper bound `a_t <= 9` consumes it at `M = -5`, the
   exact slice at which the `F` column first becomes nonempty. It is used at its
   stated strength and no further — but if `[QQ1]` failed anywhere, it would fail
   here first.

Two further honesty points the proof makes at its top, repeated here because
they belong on the front page:

* **The upper bound `a_t <= 9` has one proof, not two.** The repository at one
  time recorded two; it does not have two. There is one valuation cascade, read
  twice. See §0.3 of the proof.
* **Five inputs have zero margin** — `v_t(S) >= 11`, `v_t(h_6) >= 11`,
  `v_t(Phi) = 30` on the low side, `deg d_2 <= 6` and `z <= 6`, and
  `deg R <= 12`. See §0.4.

### What the registry records

`d2/proof_dag.json` records the target node **`C0` as `closed: true` with
`subcases_closed: 5`** — five leaves, no declared gaps. Its recorded **evidence
level is `claimed`**, because the `subcase -> C0` exhaustiveness edge is
judgment-referenced and a judgment edge caps the parent. That cap is **correct by
construction, not a backlog item**: exhaustiveness of the case partition rests on
GGHV22 Prop 4.3 and the field-split framework, which is published mathematics no
finite bookkeeping checker can re-derive. The only routes above `claimed` are a
machine-checkable reformulation of the partition, or a formal proof.

Note also that the **registry's route is not the proof's route**. The registry
reaches `C0` by an enumerative branch/cell/state ledger; the proof in
`PROOF_72_108.md` is a different and shorter route, deliberately not wired into
that ledger. The enumerative ledger is additionally **field-scoped** (built with
residue kills active); the proof is not — every step is a t-adic valuation, a
degree count, or an exact identity over `Q`.

---

## The 15-minute path (start here)

A skeptical mathematician can confirm the headline results in exact arithmetic
with just `pip install sympy`. Every one of these exits nonzero on failure and
needs no cluster.

```bash
# (a) ~2 min — D3: the dimension-3 counterexample is real
cd d3 && python3 verify.py

# (b) ~20 s — the step that used to be irreducibly machine work, now Bezout
#             certificates hand-checkable at all 40 (k,z) pairs
cd d2 && python3 support_certificates.py

# (c) ~2 min — the f37 branch is a resultant artifact (ideal membership)
cd d2 && python3 f37_sat_verify.py
```

Full ladder and expected output: **[VERIFICATION.md](VERIFICATION.md)**.

## D3 — dimension-3 counterexample (settled)

An explicit `F: C^3 -> C^3` with `det(JF) = -2` that is non-injective (three
distinct points map to `(-1/4,0,0)`). Announced by L. Alpoge 2026-07-19; the
files here are an independent exact re-verification.

```bash
cd d3 && python3 verify.py
```
See `d3/NOTES.md`.

## The wider frontier — 27 cases still open

GGV5 lists **34** candidate counterexample degree pairs with
`max(deg P, deg Q) <= 150`. The accounting:

| | |
|---|---|
| settled in the literature before this campaign | **6** |
| closed by this campaign | **1** — `(8,28)/(3,2)/108` |
| **open** | **27** |

The six are `F_1(3,4)/64`, `F_2(2,3)/75`, `F_3(3,2)/75`, `F_9(2,3)/84`,
`F_17(2,3)/99` (Moh, "discarded by hand", GGV5 tex:1794) and `F_22(2,3)/96`
(GGV5's own Proposition). Checker: `d2/moh_discards.py`.

**Evidence boundary:** that exactly six rows are red in GGV5's table is
exact-checked; that Moh's five are *ruled out* is **citation-level** — Moh's 1983
paper has not been read here, only GGV5's report of it.

`d2/CORNER_ATLAS.md` maps our machinery across all 34 rows. **The atlas
eliminates no case**: a `FAIL` verdict there means *our dictionary is unusable at
that row*, not that the row is safe.

## Verification culture

The suite is `./run_tests.sh`. Independence
is deliberate: the load-bearing cascade, infinity and alternate-regime kills are
re-derived by spec-only auditors (`audit_cascade_kills*.py`,
`audit_tplace_cases.py`, `audit_inf_cases.py`, `audit_alt_regime.py`) written
without access to the engine code they check.

Every claim is catalogued with its checker, audit status and trust tier in
`d2/PROOF_INVENTORY.md`; the proof's own step-by-step provenance table is §13.3
of `d2/PROOF_72_108.md`. A reader's map is in [docs/README.md](docs/README.md).

**Upstream sources are not redistributed.** The GGV/GGHV arXiv `.tex` files are
other authors' copyrighted manuscripts (see `d2/paper_src/README.md`). Checkers
that assert "the source says X" answer from `d2/paper_src/upstream_facts.json`
and `upstream_quotes.json` — one sha256-pinned transcription each, recording the
arXiv id and line number — and re-derive every probe from the `.tex` when it is
present locally. In a clean public clone those re-checks report that they did
**not** re-derive, rather than passing silently.

## What is NOT claimed

* **Nothing here proves or refutes JC(2).** Corollary B raises a lower bound; it
  does not settle the conjecture.
* **Priority for the exclusion is Helali's**, not ours (see the top of this file).
* **Corollary B is conditional** on Prop 4.3's exhaustiveness and on `[QQ1]`,
  both stated above and analysed in §0.2 of the proof.
* **The registry's recorded level for `C0` is `claimed`, not `proved`**, and the
  gap to the assessed ceiling `exact-checked` is a reviewed hold.
* **The upper bound is single-legged** and we do not present it as corroborated.
* The enumerative ledger route is field-scoped and, as recorded, is **not valid
  over R or C**; the proof route is immune, but the two must not be conflated.
* The generic numerical floors are *evidence, not proof*.
* Parts of the upstream boundary — the GGV common-root normalisation and the
  alpha-strip WLOG — are outline-verified, not line-by-line re-proved.

The exact verified/unverified split and ranked risks are in `d2/AUDIT.md`; open
items and negatives are §14 of the proof; the full log is `d2/STATE.md`.

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
Valqui (arXiv:2204.14178); any writeup of a D2 result belongs with or alongside
them. Priority for the (72,108) exclusion is **Billel Helali's**.
