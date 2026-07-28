#!/usr/bin/env python3
"""gghv_sub125.py  (NEW 2026-07-28; read-only over all existing artifacts)

THE OPEN FRONTIER IS 24, NOT 27.

WHAT WAS WRONG.  `corner_atlas` recorded `open_frontier = 27 = 34 - 6 - 1`:
GGV5's 34 candidate rows, minus the 6 "red" rows (Moh's five plus GGV5's own
F_22), minus the 1 this campaign closed.  `moh_discards.py` proves that
arithmetic and is CORRECT about what it checks.  The defect is that the label
"open frontier" means something the checker never computed:

  * the red partition is GGV5's `max(deg P, deg Q) <= 100` marking, and
  * GGHV22 Theorem 2.1 settles EVERYTHING with `max < 125`.

Those two never intersect.  GGHV22's own 10-case table lists every sub-125 row
with a "Discarded?" column; 9 are discarded upstream and the 10th is `(8,28)`,
left open there and closed by this campaign.  So all 10 sub-125 rows are settled,
and `34 - 10 = 24`.

Three rows were being counted OPEN while already dead upstream:
    (9,27)/(2,3)/108   discarded in GGHV22 sec.4
    F_1(5,7)/112       settled in GGV4 sec.3.5
    (8,32)/(3,2)/120   discarded in GGHV22 sec.2 ("(8,4) is not a last possible
                       corner", figure caption at 2204.14178.tex:370)

DIRECTION.  The error OVER-states how much is open.  No case-level claim is
invalidated; it mis-prices remaining work.  Same direction as, and one notch
beyond, the 32 -> 27 correction `moh_discards.py` made.

THE STANDING TRAP THIS IS.  "A checker that exits 0 has not necessarily proved
its claim -- verify the numbers describe the case the LABEL names."  The fix is
therefore NOT to assert 24: it is to COMPUTE the frontier by intersecting the
atlas with the upstream table, so the number cannot drift from its own
definition again.

Checker: --quiet, exit 0 iff every check passes.  <1 s.  Reads only.
"""
from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__)) or "."
QUIET = "--quiet" in sys.argv
_ok = [0]
_fail: list[str] = []

# ---------------------------------------------------------------------------
# GGHV22's own 10-case table, transcribed from NEXT_CASES.md, which transcribes
# arXiv:2204.14178 lines 302-320 (the "Discarded?" column) with Theorem 2.1 at
# lines 286-288.  Key: (A_0, (m,n), max_deg) -- the three fields the atlas also
# carries, so the join is on published data, not on our row ids.
# ---------------------------------------------------------------------------
GGHV_SUB125 = {
    ((4, 12), (3, 4), 64):  "GGV4 sec.3.5; [Moh]; [Heitmann]",
    ((4, 12), (5, 7), 112): "GGV4 sec.3.5",
    ((5, 20), (2, 3), 75):  "GGV3 sec.5",
    ((5, 20), (3, 2), 75):  "GGV3 sec.5",
    ((7, 21), (2, 3), 84):  "discarded in GGHV22 sec.2 (and again sec.5)",
    ((8, 24), (2, 3), 96):  "GGV5 Prop 6.1",
    ((8, 28), (3, 2), 108): "LEFT OPEN in GGHV22 -- this campaign's case",
    ((8, 32), (3, 2), 120): "discarded in GGHV22 sec.2",
    ((9, 24), (2, 3), 99):  "discarded in GGHV22 sec.4",
    ((9, 27), (2, 3), 108): "discarded in GGHV22 sec.4",
}
# The one row GGHV22 leaves open, which this campaign closes.
OURS = ((8, 28), (3, 2), 108)

# Theorem 2.1 is STRICT: "max >= 125, or (deg P,deg Q) in {(72,108),(108,72)}".
# So a row at EXACTLY 125 is NOT settled by it.  This constant is load-bearing.
THM21_BOUND = 125


def ck(name: str, cond: bool, detail: str = "") -> bool:
    if cond:
        _ok[0] += 1
        if not QUIET:
            print("[OK]   %s" % name)
    else:
        _fail.append(name)
        print("[FAIL] %s%s" % (name, ("  -- " + detail) if detail else ""))
    return bool(cond)


def key(row) -> tuple:
    return (tuple(row["A0"]), (row["m"], row["n"]), row["max_deg"])


def main() -> int:
    atlas = json.load(open(os.path.join(HERE, "corner_atlas.json"),
                           encoding="utf-8"))
    rows = atlas["rows"]

    # ---- A. the atlas side --------------------------------------------------
    ck("A1  the atlas carries GGV5's 34 candidate rows", len(rows) == 34,
       str(len(rows)))

    sub = [r for r in rows if r["max_deg"] < THM21_BOUND]
    ck("A2  exactly TEN atlas rows have max_deg < 125 -- the size of GGHV22's "
       "own sub-125 table", len(sub) == 10, str(len(sub)))

    # ---- B. the join, on published fields ------------------------------------
    atlas_keys = sorted(key(r) for r in sub)
    table_keys = sorted(GGHV_SUB125)
    ck("B1  the ten atlas rows are IN BIJECTION with GGHV22's ten, matched on "
       "(A_0, (m,n), max_deg) -- published fields, not our row ids",
       atlas_keys == table_keys,
       "atlas-only %s | table-only %s"
       % (sorted(set(atlas_keys) - set(table_keys)),
          sorted(set(table_keys) - set(atlas_keys))))

    ck("B2  exactly ONE of the ten is the case GGHV22 leaves open, and it is "
       "ours -- (8,28)/(3,2)/108",
       sum(1 for k, v in GGHV_SUB125.items() if "LEFT OPEN" in v) == 1
       and "LEFT OPEN" in GGHV_SUB125[OURS])

    ck("B3  the other nine are discarded upstream, each with a citation",
       all(v and "LEFT OPEN" not in v
           for k, v in GGHV_SUB125.items() if k != OURS))

    # ---- C. the frontier, COMPUTED -------------------------------------------
    # Settled = everything GGHV22 settles below 125, plus the one we closed.
    settled_sub125 = set(GGHV_SUB125)              # all ten, ours included
    open_rows = [r for r in rows if key(r) not in settled_sub125]
    frontier = len(open_rows)

    ck("C1  every sub-125 row is now settled -- nine upstream, one by this "
       "campaign", len(settled_sub125) == 10)
    ck("C2  the COMPUTED open frontier is 34 - 10 = 24", frontier == 24,
       str(frontier))
    ck("C3  ... and no surviving row sits below the bound",
       all(r["max_deg"] >= THM21_BOUND for r in open_rows),
       str(sorted(r["max_deg"] for r in open_rows)[:3]))

    # The three that were miscounted, named so the correction cannot be silently
    # reverted.
    red = {r["id"] for r in rows if r["provenance"].get("red_in_paper")}
    miscounted = sorted(r["id"] for r in sub
                        if r["id"] not in red and key(r) != OURS)
    ck("C4  the three rows the old 34-6-1 arithmetic counted OPEN while they "
       "were already dead upstream: %s" % miscounted,
       miscounted == ["(8,32)/(3,2)/120", "(9,27)/(2,3)/108", "F_1(5,7)/112"],
       str(miscounted))

    ck("C5  DIRECTION: the old count 27 EXCEEDS the true 24, so the error "
       "over-stated how much is open -- it mis-priced work and invalidated no "
       "case-level claim", 27 > frontier)

    # ---- D. why the old arithmetic was wrong, explicitly ---------------------
    ck("D1  the six red rows are a strict SUBSET of the ten sub-125 rows -- red "
       "is GGV5's max<=100 marking, a weaker filter than Thm 2.1's max<125",
       red < {r["id"] for r in sub} and len(red) == 6,
       "red=%d sub125=%d" % (len(red), len(sub)))
    ck("D2  and our own row is INSIDE the sub-125 set, so '34 - 6 - 1' "
       "subtracted it a second time on top of a partition that already had to "
       "contain it", key([r for r in rows if r["id"] == "(8,28)/(3,2)/108"][0])
       in settled_sub125)

    # ---- E. Thm 2.1's strictness is load-bearing -----------------------------
    at_bound = [r["id"] for r in rows if r["max_deg"] == THM21_BOUND]
    ck("E1  Thm 2.1 reads 'max >= 125', so a row at EXACTLY 125 is NOT settled "
       "by it; the atlas has one such row and it stays open: %s" % at_bound,
       len(at_bound) == 1 and at_bound[0] not in
       {r["id"] for r in rows if key(r) in settled_sub125})
    # Mutation control: a non-strict reading would wrongly close the flagship.
    would_close = [r["id"] for r in rows if r["max_deg"] <= THM21_BOUND
                   and key(r) not in settled_sub125]
    ck("E2  MUTATION CONTROL -- reading the bound as 'max <= 125' would close "
       "%s, the one row standing AT the bound. The strict reading is what keeps "
       "it open." % would_close,
       would_close == ["F_2(3,5)/125"], str(would_close))

    if _fail:
        print()
        print("FAILURES (%d):" % len(_fail))
        for f in _fail:
            print("   - %s" % f)
        return 1
    print("gghv_sub125: %d/%d checks pass -- open frontier COMPUTED as %d "
          "(34 rows - 10 settled below the 125 bound)" % (_ok[0], _ok[0], frontier))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
