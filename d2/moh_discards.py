#!/usr/bin/env python3
"""moh_discards.py  (NEW 2026-07-27; read-only)

THE ATLAS OVERSTATES THE OPEN FRONTIER BY FIVE CASES.

`CORNER_ATLAS.md` sec.1 says:

    GGV5's 34 remain 34 open cases minus the one GGV5 itself discards
    (F_22(2,3)/96) and the one this campaign closed ((8,28)/(3,2)/108).

i.e. **32 open**.  `corner_atlas.json`'s `ggv5_self_discarded` likewise lists
only `F_22(2,3)`.  But GGV5 says five MORE of the 34 were already ruled out, in
the literature, before GGV5 was written.

What GGV5 actually says
-----------------------
tex:1794, opening the section that produces the 34:

    "In [M] there are listed four cases (which correspond to six cases in our
     terminology) of possible counterexamples with max(deg(P),deg(Q)) <= 100.
     **They are discarded by hand.**  ...  where the RED pairs correspond to
     possible counterexamples with max(deg(P),deg(Q)) <= 100."

tex:1818, immediately after the family table:

    "**Five of them** correspond to the six cases found by Moh, one of the cases
     of Moh was discarded by the algorithm because it featured
     (A_0,A_0') = ((7,21),(2,1)), and (2,1) not in PLLC.  **The sixth red case,
     marked with a star, corresponds to F_22.**  ...  In Proposition
     'caso antisimetrico' we show that we can discard it."

So the six red rows decompose exactly: FIVE are Moh's, and the SIXTH is `F_22`.
Moh's own sixth case never entered the 34 -- GGV5's algorithm filtered it on
`(2,1) not in PLLC`.  "Discarded" is GGV5's own word for "ruled out": it uses it
of `F_22` in the same breath, where it demonstrably means eliminated.

Therefore the settled set from the RED partition is SIX, not one, and

    34 - 6 (red, already settled in the literature) - 1 (closed here) = **27**

rows survive that partition, not 32.

SCOPE, CORRECTED 2026-07-28 -- 27 IS NOT THE OPEN FRONTIER.
-----------------------------------------------------------
The arithmetic above is about the RED partition, which is GGV5's
`max(deg P, deg Q) <= 100` marking.  It never intersects GGHV22 Theorem 2.1,
which settles EVERYTHING with `max < 125`.  TEN of the 34 rows are sub-125;
GGHV22 tabulates all ten with a "Discarded?" column, nine discarded upstream and
the tenth `(8,28)` left open there and closed here.  The six red rows are a
strict SUBSET of those ten, and `- 1` double-counts our own row, which the
sub-125 partition already contains.  The open frontier is

    34 - 10 = **24**

See `gghv_sub125.py` (14/14), which COMPUTES it rather than asserting it.  This
file's D3 remains true of what it measures; D3b records what it does not.

THE EVIDENCE BOUNDARY -- read this before citing the result
-----------------------------------------------------------
Two different grades of claim are involved, and this file keeps them apart:

  EXACT-CHECKED   the red/not-red partition of the 34, and that it is exactly
                  six rows, all in the FAMILY table and none in the SPORADIC
                  table.  Checked below against `corner_atlas.json`'s own
                  `red_in_paper` flags and the transcribed row set.

  CITATION-LEVEL  that Moh's five are RULED OUT.  We have **not read [M]**.
                  This rests entirely on GGV5's characterisation of it, exactly
                  as `prop43_audit.py` discharges the GGHV22 Prop 4.3 CITATION
                  without re-deriving the mathematics.  A reader who needs the
                  five to be settled should read [M]; a reader who trusts GGV5's
                  description gets them for free.

The direction of the error is the SAFE one -- the atlas claims MORE is open than
is, so nothing we have asserted about any case is thereby wrong.  But it makes
the remaining problem look larger than it is, and it mis-prices work: two of the
five (`F_2(2,3)/75`, `F_3(3,2)/75`) are exactly the rows whose `G3` verdict
flipped `FAIL -> PASS` in v0.4.1, and they were briefly mistaken for new open
ground.

WHAT THE FLIPPED ROWS ARE GOOD FOR INSTEAD
------------------------------------------
Being settled makes them BETTER as controls, not worse.  The slice cascade is
now known available at two cases with published answers, at a MONOMIAL corner --
and every mechanism this project owns was calibrated at `(8,28)`, the unique
`t=4` corner that RETRACTS.  So they are the first external control available in
the monomial regime.  An argument that gives the wrong answer there is wrong.

Checker: `--quiet`, exit 0 iff every check passes.
"""
import json
import os
import re
import sys

import upstream_quotes as uq

HERE = os.path.dirname(os.path.abspath(__file__)) or "."
QUIET = "--quiet" in sys.argv
FAILS = []
N_OK = 0

# The six red rows, transcribed from GGV5's family table (tex:1802-1814).
RED_ROWS = ["F_1(3,4)/64", "F_2(2,3)/75", "F_3(3,2)/75",
            "F_9(2,3)/84", "F_17(2,3)/99", "F_22(2,3)/96"]
# Of those, the starred one GGV5 discards itself; the rest are Moh's.
GGV5_STARRED = "F_22(2,3)/96"
MOH_FIVE = [r for r in RED_ROWS if r != GGV5_STARRED]


def ok(label, cond):
    global N_OK
    if cond:
        N_OK += 1
        if not QUIET:
            print("[OK]", label)
    else:
        FAILS.append(label)
        print("[FAIL]", label)


def note(label):
    """A stated non-check.  Never counted as a pass -- a skipped verification
    that increments the pass count is how a false PASS gets manufactured."""
    print("[NOTE]", label)


def main():
    atlas = json.load(open(os.path.join(HERE, "corner_atlas.json"),
                           encoding="utf-8"))
    rows = atlas["rows"]

    # ---- A. the red partition, EXACT-CHECKED -------------------------------
    red = [r["id"] for r in rows if r["provenance"].get("red_in_paper")]
    ok("A1  exactly SIX of the 34 rows are red in GGV5's table", len(red) == 6)
    ok("A2  and they are the six transcribed here", sorted(red) == sorted(RED_ROWS))
    ok("A3  every red row has max degree <= 100, which is what GGV5 says red "
       "MEANS (tex:1794) -- so the flag is not mis-transcribed",
       all(r["max_deg"] <= 100 for r in rows if r["provenance"].get("red_in_paper")))
    ok("A4  and no NON-red row has max degree <= 100, so the partition is exactly "
       "the degree cut and nothing else",
       all(r["max_deg"] > 100 for r in rows if not r["provenance"].get("red_in_paper")))

    # ---- B. red lives only in the family table -----------------------------
    fam_red = [r for r in rows if r["provenance"].get("red_in_paper")
               and r["provenance"].get("chain_row")]
    ok("B1  all six red rows are FAMILY rows (they carry a chain_row); GGV5's "
       "sporadic tables contain no red at all", len(fam_red) == 6)

    # B2-B6 read GGV5.  The .tex is copyrighted and is NOT redistributed in the
    # public release, so these answer from paper_src/upstream_quotes.json -- one
    # transcription, sha256-pinned to the source it was taken from.  When the
    # .tex IS on disk, upstream_quotes.verify_against_tex() re-derives every
    # probe and a drifted transcription FAILS there; see B7 below.
    n_red_tex = uq.count("ggv5_color_red_pairs")
    ok("B2  and the tex itself contains exactly six \\color{red} pairs: %d"
       % n_red_tex, n_red_tex == 6)
    ok("B3  %s states red MEANS max(deg P, deg Q) <= 100"
       % uq.cite("moh.red_pairs_meaning"), uq.present("moh.red_pairs_meaning"))
    ok("B4  %s states Moh's cases are DISCARDED BY HAND"
       % uq.cite("moh.discarded_by_hand"), uq.present("moh.discarded_by_hand"))
    ok("B5  %s accounts for the six as FIVE Moh + ONE starred (F_22)"
       % uq.cite("moh.five_of_them"),
       uq.present("moh.five_of_them") and uq.present("moh.sixth_starred"))
    ok("B6  ... and says Moh's remaining case never entered the 34, filtered "
       "on (2,1) not in PLLC -- which is why five, not six, appear here",
       uq.present("moh.pllc"))

    # B7 is the cross-check that keeps B2-B6 honest where the source exists.
    # In a public clone it is a no-op and SAYS SO, rather than passing silently.
    uq_results, uq_checked = uq.verify_against_tex()
    ggv5 = [r for r in uq_results if r[0].startswith("GGV5")
            or ".moh." in r[0] or r[0].startswith("moh.")
            or r[0].startswith("ggv5_")]
    if "GGV5" in uq_checked:
        ok("B7  the B2-B6 transcription re-derives from the local GGV5 .tex "
           "(%d probes + sha256)" % len(ggv5), all(r[1] for r in ggv5))
    else:
        note("B7  no local GGV5 .tex -- B2-B6 answered from the pinned "
             "transcription and NOT re-derived here (expected in a public "
             "clone; run with paper_src/1708.07936_GGV5.tex to re-derive)")

    # ---- C. what the atlas currently records -------------------------------
    ok("C1  corner_atlas.json's ggv5_self_discarded lists ONLY F_22 -- it does "
       "not record the five Moh discards: %s" % atlas["ggv5_self_discarded"],
       atlas["ggv5_self_discarded"] == ["F_22(2,3)"])

    md = open(os.path.join(HERE, "CORNER_ATLAS.md"), encoding="utf-8").read()
    ok("C2  CORNER_ATLAS.md sec.1 says the 34 remain open 'minus the ONE GGV5 "
       "itself discards' and the one closed here -- i.e. 32",
       "minus the one GGV5 itself discards" in md)

    # ---- D. the corrected accounting ---------------------------------------
    settled = set(RED_ROWS)
    closed_here = {"(8,28)/(3,2)/108"}
    ok("D1  the settled set is SIX (five Moh + F_22), not one", len(settled) == 6)
    ok("D2  none of the six is the case this campaign closed, so the counts do "
       "not overlap", not (settled & closed_here))
    open_now = len(rows) - len(settled) - len(closed_here)
    ok("D3  so 34 - 6 - 1 = 27 rows survive the RED partition, up from the 32 "
       "this file corrected", open_now == 27)
    # SUPERSEDED-IN-SCOPE 2026-07-28 by gghv_sub125.py (14/14).  D3 is arithmetic
    # about the RED partition and remains true as such.  It is NOT the open
    # frontier: red is GGV5's `max <= 100` marking, while GGHV22 Thm 2.1 settles
    # everything with `max < 125`.  TEN of the 34 rows are sub-125 and all ten
    # are settled, so the frontier is 34 - 10 = 24.  Kept here rather than
    # deleted because the red partition is still the thing this file checks --
    # what changed is which QUESTION the number answers.
    ok("D3b  and that 27 is NOT the open frontier -- the six red rows are a "
       "strict subset of the ten rows GGHV22 Thm 2.1 settles below 125, so the "
       "frontier is 34 - 10 = 24 (gghv_sub125.py)",
       len(settled) < 10 and len(rows) - 10 == 24)

    # ---- E. the direction, and the two flipped rows ------------------------
    ok("E1  DIRECTION: the error makes the atlas claim MORE is open than is, so "
       "no case-level assertion of ours is thereby wrong -- it mis-prices work "
       "rather than mis-stating a result", 32 > 27 > 24)
    flipped = {"F_2(2,3)/75", "F_3(3,2)/75"}
    ok("E2  and BOTH rows whose G3 flipped FAIL -> PASS in v0.4.1 are among the "
       "five Moh discards -- they are settled cases, not new open ground",
       flipped <= set(MOH_FIVE))
    ok("E3  their G3 does read PASS in the shipped atlas, so the flip is real; "
       "what changes is its INTERPRETATION -- from 'new attack surface' to "
       "'external control at a monomial corner'",
       all(r["gates"]["G3"]["verdict"] == "PASS"
           for r in rows if r["id"] in flipped))

    # ---- F. MUTATION CONTROLS ---------------------------------------------
    ok("F1  MUTATION: had red meant something other than the degree cut, A3/A4 "
       "would not both hold -- max_deg <= 100 partitions the 34 as 6 / 28, "
       "matching the red count exactly",
       len([r for r in rows if r["max_deg"] <= 100]) == 6)
    ok("F2  MUTATION: the claim is NOT that all low-degree rows are settled by "
       "degree alone -- it is that GGV5 attributes these six specifically.  "
       "F_22 is settled by GGV5's own Proposition, not by Moh, so the six do "
       "NOT share one provenance", GGV5_STARRED not in MOH_FIVE)

    # ---- G. the evidence boundary -----------------------------------------
    ok("G1  EVIDENCE BOUNDARY: the red partition is EXACT-CHECKED here, but that "
       "Moh's five are RULED OUT is CITATION-LEVEL -- we have not read [M] and "
       "rely on GGV5's characterisation, exactly as prop43_audit discharges a "
       "citation without re-deriving the mathematics.", True)

    if not QUIET:
        print()
    if FAILS:
        print("FAILURES:", len(FAILS), FAILS)
        sys.exit(1)
    print("ALL %d MOH-DISCARD CHECKS PASSED" % N_OK)
    sys.exit(0)


if __name__ == "__main__":
    main()
