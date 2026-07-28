# docs — reader's map

This is a pointer index, not a copy. The canonical documents live in `d2/`
(single source of truth); this page maps a reader's path through them so nothing
here can diverge from the originals.

## The path: proof → priority → conditions → provenance → the wider frontier

1. **Start — the proof**
   [`d2/PROOF_72_108.md`](../d2/PROOF_72_108.md) — the independent structural
   proof that both Prop-4.3 configurations are empty, hence (conditionally) that
   any plane counterexample has `max(deg P, deg Q) >= 125`. Self-contained,
   states its own conditionality at the top, and carries its own provenance table
   (§13.3) mapping every step to the checker that verifies it. **Read this
   first.** Everything else in this list exists to support or bound it.

2. **Priority — Helali was first**
   [`d2/HELALI_ADJUDICATION.md`](../d2/HELALI_ADJUDICATION.md) — B. Helali's
   independent exclusion (21 July 2026,
   [doi:10.5281/zenodo.21479814](https://doi.org/10.5281/zenodo.21479814)) is the
   public record for this result. This document is our independent audit of his
   reduction and the adjudication **SUBSUMES**, including the five candidate leak
   types we enumerated and closed. Checker: `d2/helali_adjudication_check.py`.

3. **The conditions the result rests on**
   [`d2/PROOF_INVENTORY.md`](../d2/PROOF_INVENTORY.md) — the C0–C46 claim graph,
   every claim's checker, its independent-audit status, and a trust tier
   (1 = independently audited … 4 = published result used as stated). The two
   load-bearing imports are GGHV22 Prop 4.3 (checker `d2/prop43_audit.py`) and
   the alpha-strip WLOG `[QQ1]`, which is Proposition 2.1 of the proof.

4. **The registry, and what its level means**
   [`d2/PROOF_DAG.md`](../d2/PROOF_DAG.md) and `d2/proof_dag.json` — the machine
   claim registry. `C0` is recorded `closed: true, subcases_closed: 5` at
   evidence level `claimed`; the cap is the judgment-referenced exhaustiveness
   edge and is **correct by construction**, not a backlog item. Note the registry
   reaches `C0` by an *enumerative* route which is **not** the proof's route, and
   which is field-scoped in a way the proof is not — see §0.5 of the proof.

5. **A component result — the f37 branch is a resultant artifact**
   [`d2/F37_SATURATION_REPORT.md`](../d2/F37_SATURATION_REPORT.md) — the
   ideal-membership closure of the whole f37 branch, which removed one of three
   apparent branches from the older elimination route. Checker:
   `d2/f37_sat_verify.py` (~2 min, exact). This was the headline before the proof
   existed; it is now one component of the enumerative route.

6. **The wider frontier — 24 cases still open**
   [`d2/CORNER_ATLAS.md`](../d2/CORNER_ATLAS.md) — our machinery mapped across
   GGV5's 34 candidate degree pairs. **Ten** sit below the 125 bound and all ten
   are settled (nine discarded upstream in GGHV22's own table, the tenth being
   our `(8,28)`), leaving **24 open**, every one at `max_deg >= 125`. The figure
   is computed, not asserted: `d2/gghv_sub125.py` (14/14) joins the atlas to
   GGHV22's table. `F_2(3,5)/125` sits exactly *at* the bound and is the unique
   row there. **The atlas eliminates no case** — a `FAIL` there means our
   dictionary is unusable at that row, not that the row is safe.

7. **The full derivation log**
   [`d2/STATE.md`](../d2/STATE.md) — the chronological record behind every claim.

## Supporting documents

- **Referee-facing paper skeleton:** [`d2/WRITEUP_OUTLINE.md`](../d2/WRITEUP_OUTLINE.md)
- **Verified/unverified split and ranked risks:** [`d2/AUDIT.md`](../d2/AUDIT.md)
- **The field-split repair (quartic q as four places):** [`d2/FIELD_SPLIT_AUDIT.md`](../d2/FIELD_SPLIT_AUDIT.md)
- **Ablations — which inputs are load-bearing:** [`d2/MINIMAL_CORE.md`](../d2/MINIMAL_CORE.md)
- **The (50,75) external control:** [`d2/MOH_CONTROL_50_75.md`](../d2/MOH_CONTROL_50_75.md)
- **Forward plan (global residue obstruction algebra):** [`d2/PHASE_F_PLAN.md`](../d2/PHASE_F_PLAN.md)
- **Source papers (links only, not redistributed):** [`d2/paper_src/README.md`](../d2/paper_src/README.md)

## Two stale artifacts, named rather than left to be discovered

- [`d2/FRONTIER.md`](../d2/FRONTIER.md) is machine-generated from the
  **enumerative** route's JSON artifacts and is pinned to an older commit. Its
  cell counts describe that route's intermediate state, **not** the current
  result, and it is not regenerated for this release. Do not read it as "what is
  left open" for (72,108).
- [`d2/CURRENT_STATUS.md`](../d2/CURRENT_STATUS.md) carries a correction banner
  above a body that predates it. Where the two disagree, `PROOF_72_108.md` and
  `proof_dag.json` win.

**Naming note.** The historical documents inside `d2/` refer to their own
directory by its original working name `d2_plane_72_108/`. In this release tree
that directory is `d2/`; internal relative references between sibling files are
unchanged.

## How to verify

See the top-level [`VERIFICATION.md`](../VERIFICATION.md) for the
exact-arithmetic ladder, and run `./run_tests.sh` for the full suite.
