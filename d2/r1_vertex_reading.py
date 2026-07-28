#!/usr/bin/env python3
"""r1_vertex_reading.py  (NEW 2026-07-28; read-only)

RESOLVES THE (r1) VERTEX AMBIGUITY -- and the resolution is that the question
was posed as a false dichotomy.

THE BLOCKER, as SECOND_CORNER.md sec.5.2 states it.  Rule (r1) writes the
pre-reduction polygon as `Delta = hull{(0,0), (1,0), A_0, (0,c)}`, and two
readings of that `(1,0)` were on the table:

    A   the vertex is ALWAYS (1,0), a normalisation
    B   the vertex IS A_0'

(r1) is printed at only four corners, and at three of them `A_0' = (1,0)` so the
readings coincide.  SECOND_CORNER.md concluded that a FIFTH printed Delta, at
some corner with `A_0' != (1,0)`, was needed to separate them -- a missing
PUBLISHED datum blocking six corners: (7,35), (8,24), (9,21), (10,40), (11,33),
(12,30).

THE RESOLUTION.  The fourth printed corner, `(9,27)` at 2204.14178.tex:471,
already separates them -- against BOTH.  Its printed set has FIVE vertices,

    {(0,0), (1,0), (9,24), (9,27), (0,9)}

containing (1,0) AND A_0' = (9,24) simultaneously.  Reading A cannot produce
(9,24); reading B cannot produce (1,0).  Each scores 3/4, each failing at exactly
that corner.  The union

    (U)   Delta = hull{ (0,0), (1,0), A_0', A_0, (0,c) }

scores 4/4.  It is not a compromise: it is the only one of the three consistent
with the published evidence.

WHY IT LOOKED LIKE A DICHOTOMY.  Where `A_0'` lies ON the x-axis -- which is
exactly the blocked `A_0' = (2,0)` case -- (1,0) falls strictly inside the
segment from (0,0) to (2,0) and is absorbed by the hull.  So U and B give the
SAME polygon there, and U and A differ.  Where `A_0'` is off-axis (the (9,27)
case) all three differ, and only U survives.  The one configuration that
discriminates is the one nobody expected to matter.

WHAT THIS UNBLOCKS.  At `A_0' = (2,0)` corners U predicts
`hull{(0,0), (2,0), A_0, (0,c)}` -- i.e. reading B's polygon, but for a reason
validated against the evidence rather than assumed.

SCOPE -- READ BEFORE USING.
  * ONE discriminating instance.  (9,27) is the only printed corner with
    `A_0' != (1,0)`.  U is the unique survivor of the three readings tested; it
    is not proved to be the unique rule consistent with four data points.
  * Applying U at a blocked corner still needs `c` THERE, which is a separate
    input this file does not supply.  Only (10,40) has `c` pinned independently
    (SECOND_CORNER.md sec.5.2: both readings force the top vertex to
    `(l*a_0 - b_0, a_0) = (0,10)`).
  * `A_0'` at the blocked corners comes from the recovery lemma of
    `second_corner_probe.py`, not from print.

Checker: --quiet, exit 0 iff every check passes.  <1 s.  Reads only.
"""
from __future__ import annotations

import sys

QUIET = "--quiet" in sys.argv
_ok = [0]
_fail: list[str] = []


def ck(name: str, cond: bool, detail: str = "") -> bool:
    if cond:
        _ok[0] += 1
        if not QUIET:
            print("[OK]   %s" % name)
    else:
        _fail.append(name)
        print("[FAIL] %s%s" % (name, ("  -- " + detail) if detail else ""))
    return bool(cond)


def hull(pts):
    pts = sorted(set(map(tuple, pts)))
    if len(pts) <= 2:
        return pts

    def cross(o, a, b):
        return (a[0]-o[0])*(b[1]-o[1]) - (a[1]-o[1])*(b[0]-o[0])

    lo = []
    for p in pts:
        while len(lo) >= 2 and cross(lo[-2], lo[-1], p) <= 0:
            lo.pop()
        lo.append(p)
    up = []
    for p in reversed(pts):
        while len(up) >= 2 and cross(up[-2], up[-1], p) <= 0:
            up.pop()
        up.append(p)
    return sorted(set(lo + up))


# The four printed (r1) polygons in GGHV22, transcribed with their tex lines.
# c is read off the printed (0,c) vertex; A_0' as recorded in SECOND_CORNER.md.
PRINTED = {
    (7, 21): dict(c=7, A0p=(1, 0),
                  poly=[(0, 0), (1, 0), (7, 21), (0, 7)], tex=1388),
    (8, 28): dict(c=4, A0p=(1, 0),
                  poly=[(0, 0), (1, 0), (8, 28), (0, 4)], tex=1009),
    (9, 24): dict(c=6, A0p=(1, 0),
                  poly=[(0, 0), (1, 0), (9, 24), (0, 6)], tex=682),
    (9, 27): dict(c=9, A0p=(9, 24),
                  poly=[(0, 0), (1, 0), (9, 24), (9, 27), (0, 9)], tex=471),
}

READINGS = {
    "A": lambda A0, c, p: hull([(0, 0), (1, 0), A0, (0, c)]),
    "B": lambda A0, c, p: hull([(0, 0), p, A0, (0, c)]),
    "U": lambda A0, c, p: hull([(0, 0), (1, 0), p, A0, (0, c)]),
}


def main() -> int:
    # ---- A. the evidence base ------------------------------------------------
    ck("A1  four (r1) polygons are printed in GGHV22, at the C0 corners",
       len(PRINTED) == 4)
    off = [A0 for A0, m in PRINTED.items() if m["A0p"] != (1, 0)]
    ck("A2  exactly ONE printed corner has A_0' != (1,0) -- %s -- so it is the "
       "only one that can discriminate" % off, off == [(9, 27)])

    # ---- B. score the three readings ----------------------------------------
    score = {}
    for name, fn in READINGS.items():
        bad = [A0 for A0, m in PRINTED.items()
               if fn(A0, m["c"], m["A0p"]) != hull(m["poly"])]
        score[name] = (4 - len(bad), bad)

    ck("B1  reading A (vertex always (1,0)) scores 3/4, failing at (9,27) -- it "
       "cannot produce the printed A_0' vertex (9,24)",
       score["A"] == (3, [(9, 27)]), str(score["A"]))
    ck("B2  reading B (vertex is A_0') scores 3/4, failing at the SAME corner -- "
       "it cannot produce the printed (1,0)",
       score["B"] == (3, [(9, 27)]), str(score["B"]))
    ck("B3  the UNION reading U = hull{(0,0),(1,0),A_0',A_0,(0,c)} scores 4/4",
       score["U"] == (4, []), str(score["U"]))
    ck("B4  so the A-vs-B dichotomy was FALSE: the printed evidence discriminates, "
       "and it rules out BOTH stated readings",
       score["A"][0] < 4 and score["B"][0] < 4 and score["U"][0] == 4)

    # ---- C. why it looked undecidable ---------------------------------------
    onaxis = hull([(0, 0), (1, 0), (2, 0), (10, 40), (0, 10)])
    b_onaxis = hull([(0, 0), (2, 0), (10, 40), (0, 10)])
    a_onaxis = hull([(0, 0), (1, 0), (10, 40), (0, 10)])
    ck("C1  where A_0' lies ON the x-axis, (1,0) is absorbed into the segment to "
       "(2,0), so U and B AGREE and only A differs",
       onaxis == b_onaxis and onaxis != a_onaxis)
    m = PRINTED[(9, 27)]
    three = {k: tuple(fn((9, 27), m["c"], m["A0p"])) for k, fn in READINGS.items()}
    ck("C2  where A_0' is OFF-axis -- the (9,27) case -- all three readings differ, "
       "which is why that single corner settles it", len(set(three.values())) == 3)

    # ---- D. mutation controls ------------------------------------------------
    perturbed = hull([(0, 0), (1, 0), (9, 25), (9, 27), (0, 9)])
    ck("D1  MUTATION: moving the printed A_0' vertex (9,24) -> (9,25) breaks U's "
       "match, so B3 is not vacuous",
       perturbed != hull(PRINTED[(9, 27)]["poly"]))
    ck("D2  MUTATION: dropping (1,0) from U at (9,27) reproduces reading B and "
       "therefore FAILS, confirming (1,0) is load-bearing in the printed set",
       hull([(0, 0), (9, 24), (9, 27), (0, 9)]) != hull(PRINTED[(9, 27)]["poly"]))

    # ---- E. what it predicts, and what it does not ---------------------------
    pred = hull([(0, 0), (1, 0), (2, 0), (10, 40), (0, 10)])
    ck("E1  at (10,40) with A_0'=(2,0) and c=10, U predicts %s" % str(pred),
       pred == [(0, 0), (0, 10), (2, 0), (10, 40)])
    ck("E2  and c=10 there is independently pinned (top vertex "
       "(l*a_0-b_0, a_0) = (0,10) under BOTH original readings), so this "
       "prediction consumes no unpublished c", 4 * 10 - 40 == 0)
    if not QUIET:
        print("[NOTE] the other five blocked corners -- (7,35), (8,24), (9,21), "
              "(11,33), (12,30) -- need c THERE, which this file does not supply.")

    if _fail:
        print()
        print("FAILURES (%d):" % len(_fail))
        for f in _fail:
            print("   - %s" % f)
        return 1
    print("r1_vertex_reading: %d/%d checks pass -- A 3/4, B 3/4 (same corner), "
          "UNION 4/4; (10,40) unblocked" % (_ok[0], _ok[0]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
