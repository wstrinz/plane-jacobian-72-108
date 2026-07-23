# docs — reader's map

This is a pointer index, not a copy. The canonical documents live in `d2/`
(single source of truth); this page maps a reader's path through them so nothing
here can diverge from the originals.

## The path: theorem → inventory → frontier → state log

1. **Start — what is claimed and how strong each claim is**
   [`d2/PROOF_INVENTORY.md`](../d2/PROOF_INVENTORY.md) — the C0–C46 claim graph,
   every claim's checker, its independent-audit status, and a trust tier
   (1 = independently audited … 4 = published result used as stated). The single
   source of truth.

2. **The headline theorem — f37 is a resultant artifact**
   [`d2/F37_SATURATION_REPORT.md`](../d2/F37_SATURATION_REPORT.md) — the
   ideal-membership closure of the whole f37 branch. Checker: `d2/f37_sat_verify.py`.

3. **The live frontier — what is left open**
   [`d2/FRONTIER.md`](../d2/FRONTIER.md) — machine-generated (do not hand-edit)
   from the cascade / alternate-regime JSON artifacts by `d2/frontier_gen.py`:
   26 subcase-2 cells + 171 subcase-1 branches + 27 alternate-regime branches.

4. **The full derivation log**
   [`d2/STATE.md`](../d2/STATE.md) — the chronological record behind every claim.

## Supporting documents

- **Referee-facing paper skeleton:** [`d2/WRITEUP_OUTLINE.md`](../d2/WRITEUP_OUTLINE.md)
- **Verified/unverified split and ranked risks:** [`d2/AUDIT.md`](../d2/AUDIT.md)
- **The field-split repair (quartic q as four places):** [`d2/FIELD_SPLIT_AUDIT.md`](../d2/FIELD_SPLIT_AUDIT.md)
- **Forward plan (global residue obstruction algebra):** [`d2/PHASE_F_PLAN.md`](../d2/PHASE_F_PLAN.md)
- **Source papers (links only):** [`d2/paper_src/README.md`](../d2/paper_src/README.md)

**Naming note.** The historical documents inside `d2/` refer to their own
directory by its original working name `d2_plane_72_108/`. In this release tree
that directory is `d2/`; internal relative references between sibling files are
unchanged.

## How to verify

See the top-level [`VERIFICATION.md`](../VERIFICATION.md) for the 15-minute
exact-arithmetic ladder, and run `./run_tests.sh` for the full suite.
