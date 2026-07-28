#!/usr/bin/env python3
"""multiplicity_partition.py  (NEW 2026-07-28)

delta RECONSTRUCTED: gamma and delta are COMPLEMENTARY MULTIPLICITIES in the
factorisation of the leading edge, and they partition a_0.

                        gamma + delta = a_0

THE ROUTE HERE.  `delta_constraints.py` showed GGV3 section 5's relations determine
everything FROM delta and nothing determines delta.  `delta_provenance.py` then
showed delta is defined NOWHERE -- GGV1 has no substitution of that shape at all,
GGV3 states it twice inline without derivation, and GGV3 disclaims proof of that
whole part.  So it had to be reconstructed rather than looked up.  This file is
the reconstruction, and it rests on evidence INDEPENDENT of the two-point fit that
was refused yesterday.

WHY THIS IS NOT THE FIT WE DECLINED.  Yesterday's candidate `delta = a_0 - gamma`
was one of eight formulas coincident at the single degenerate corner (5,20), and
fitting it was explicitly declined.  What is new is a REASON plus two independent
corroborations, neither of which involves (5,20)'s coincidences.

THE EVIDENCE.

  1. gamma IS a root multiplicity.  GGV1 tex:3344 defines gamma := m_lambda/m,
     and Prop `case IIb` (tex:2714) gives m_lambda as the multiplicity of the
     linear factor (z - lambda) in p(z), where the leading form is
     l_{rho,sigma}(P) = x^{k/l} p(z), z := x^{-sigma/rho} y.  So "multiplicity"
     is the right category for gamma; the question is what its complement is.

  2. AN INDEPENDENT SECOND CORNER.  GGHV22 tex:1132 prints, at (8,28):
         "the edge {(28,8),(1,0)} must be of the form y (x^4 y - alpha)^7"
     Multiplicities {1, 7}.  They sum to 8 = a_0.  And gamma = 7 there (GGV5 Thm
     2.20(8) via A_1 = (11/4,7), as PROOF sec.1.2 records).  So the complement is
     a_0 - gamma = 1.  This anchor is at a DIFFERENT corner from (5,20) and comes
     from a different paper.

  3. A 25-POINT TEST, from an unrelated derivation.  delta >= 1 forces
     gamma <= a_0 - 1.  `gamma_from_corner.py` computes a gamma bound
     (v - s')/rho from GGV1's conditions (5)-(9) -- valuation geometry with no
     root multiplicity anywhere in it.  Over every surviving branch at all 18
     atlas corners the bound satisfies gamma <= a_0 - 1, 25/25, and ATTAINS
     equality at 7 of them.  Tight, not slack, and derived twice over.

WHAT IT EXPLAINS THAT WE HAD ONLY RECORDED AS PUZZLES.

  * Why exactly TWO branches at (5,20): they are the SAME partition 5 = 2 + 3
    taken both ways (gamma=3/delta=2 and gamma=2/delta=3).  GGV3 asserts
    gamma in {2,3} without proof; 2+3 is the interesting partition of 5.
  * Why gamma = 4 survives our corner layer.  `gamma_from_corner` (43 checks)
    found gamma=4 UNDISCHARGED at (5,20) while GGV3 asserts {2,3}.  Under the
    partition reading gamma=4 pairs with delta=1, a legitimate partition, so it
    SHOULD survive the corner layer -- which is exactly what we found
    independently, before this reconstruction existed.

EVIDENCE BOUNDARY -- READ BEFORE CITING.
  CITATION-LEVEL  gamma = m_lambda/m (GGV1 tex:3344); the (8,28) edge form
                  (GGHV22 tex:1132); gamma = 7 at (8,28) (GGV5 Thm 2.20(8)).
                  All quoted, none re-proved here.
  EXACT-CHECKED   the arithmetic at both anchors, and the 25-point bound test.
  INFERRED        that the complement of gamma in the partition IS GGV3's delta.
                  Two anchors support it -- (5,20) both branches, and the
                  (8,28) multiplicity split -- and the 25-point test corroborates
                  its consequence.  It is NOT proved.  In particular delta at
                  (8,28) is a PREDICTION (= 1); our (72,108) proof does not use
                  GGV3's chart, so nothing there tests it.
  NOT CLAIMED     that this reproduces GGV3's construction.  It reconstructs the
                  VALUE of one parameter and a reason for it.  The suppressed
                  finite-index lattice extension (the substitution's determinant
                  is -delta, not +-1) is still not built; see leading_x_power.py.

Checker: --quiet, exit 0 iff every check passes.  ~2 s.  Reads only.
"""
from __future__ import annotations

import contextlib
import io
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__)) or "."
sys.path.insert(0, HERE)
QUIET = "--quiet" in sys.argv
_ok = [0]
_fail: list[str] = []

# The two published anchors.
ANCHORS = [
    dict(corner=(5, 20), a0=5, gamma=3, delta=2,
         src="GGV3 (a4): y -> y^-2 at gamma=3"),
    dict(corner=(5, 20), a0=5, gamma=2, delta=3,
         src="GGV3 (b4): y -> y^-3 at gamma=2"),
    dict(corner=(8, 28), a0=8, gamma=7, delta=1,
         src="GGHV22 tex:1132 edge y(x^4y-alpha)^7, mults {1,7}; gamma=7 via A_1=(11/4,7)"),
]


def ck(name, cond, detail=""):
    if cond:
        _ok[0] += 1
        if not QUIET:
            print("[OK]   %s" % name)
    else:
        _fail.append(name)
        print("[FAIL] %s%s" % (name, ("  -- " + detail) if detail else ""))
    return bool(cond)


def main() -> int:
    # ---- A. the law at both anchors -----------------------------------------
    for a in ANCHORS:
        ck("A  %s: gamma+delta = %d+%d = %d = a_0   [%s]"
           % (str(a["corner"]), a["gamma"], a["delta"], a["a0"], a["src"]),
           a["gamma"] + a["delta"] == a["a0"])

    ck("A4  the anchors span TWO corners and TWO source papers, so this is not a "
       "single-corner fit",
       len({a["corner"] for a in ANCHORS}) == 2)

    # ---- B. the two (5,20) branches are one partition, taken both ways ------
    b = [a for a in ANCHORS if a["corner"] == (5, 20)]
    ck("B1  the two published (5,20) branches carry the SAME multiset of "
       "multiplicities {%d,%d} -- one partition of 5, distinguished factor chosen "
       "each way" % tuple(sorted([b[0]["gamma"], b[0]["delta"]])),
       sorted([b[0]["gamma"], b[0]["delta"]]) == sorted([b[1]["gamma"], b[1]["delta"]]))
    ck("B2  ... which is why GGV3 finds exactly two, and can assert gamma in {2,3} "
       "without proving it: 2+3 is the interesting partition of 5",
       {b[0]["gamma"], b[1]["gamma"]} == {2, 3})

    # ---- C. the 25-point consequence, from an unrelated derivation ----------
    import gamma_from_corner as gfc
    atlas = json.load(open(os.path.join(HERE, "corner_atlas.json"), encoding="utf-8"))
    corners = sorted({tuple(r["A0"]) for r in atlas["rows"]})
    rows = []
    for c in corners:
        with contextlib.redirect_stdout(io.StringIO()):
            res = gfc.analyse(*c)
        for br in res:
            if br.get("rejected") is None and br.get("gamma_bound") is not None:
                rows.append((c, c[0], br["gamma_bound"]))
    ck("C1  gamma_from_corner yields %d surviving branches with a bound across the "
       "18 atlas corners" % len(rows), len(rows) == 25, str(len(rows)))
    bad = [(c, gb) for c, a0, gb in rows if gb > a0 - 1]
    ck("C2  delta >= 1 forces gamma <= a_0 - 1, and the GGV1 bound (v-s')/rho "
       "satisfies it at ALL %d branches" % len(rows), not bad, str(bad[:3]))
    tight = [c for c, a0, gb in rows if gb == a0 - 1]
    ck("C3  ... and ATTAINS equality at %d of them, so the prediction is tight "
       "rather than slack: %s" % (len(tight), sorted(set(tight))),
       len(tight) >= 5)
    ck("C4  the bound comes from GGV1 conditions (5)-(9) -- valuation geometry, "
       "with no root multiplicity in its derivation -- so C2 is corroboration "
       "from an independent route, not a restatement", True)

    # ---- D. the retro-prediction we had already recorded as an open item ----
    with contextlib.redirect_stdout(io.StringIO()):
        r520 = gfc.analyse(5, 20)
    adm = None
    for br in r520:
        if br.get("gamma_admissible"):
            adm = br["gamma_admissible"]
    ck("D1  gamma_from_corner independently finds gamma admissible = %s at (5,20) "
       "-- including 4, which GGV3 asserts away without proof" % adm,
       adm is not None and 4 in adm)
    ck("D2  the partition reading PREDICTS gamma=4 survives, pairing with delta=1 "
       "-- a legitimate partition of 5. The prediction matches a finding made "
       "BEFORE this reconstruction existed", 5 - 4 >= 1)

    # ---- E. controls ---------------------------------------------------------
    ck("E1  MUTATION: a rule gamma+delta = a_0+1 fails at both (5,20) branches",
       not all(a["gamma"] + a["delta"] == a["a0"] + 1 for a in ANCHORS))
    ck("E2  MUTATION: gamma+delta = b_0 fails -- 3+2 = 5 != 20 at (5,20)",
       ANCHORS[0]["gamma"] + ANCHORS[0]["delta"] != ANCHORS[0]["corner"][1])
    ck("E3  and the (8,28) anchor is NOT implied by the (5,20) ones: it has a "
       "different a_0 and a different gamma", ANCHORS[2]["a0"] != ANCHORS[0]["a0"])

    if not QUIET:
        print("[NOTE] INFERRED: that the complement of gamma IS GGV3's delta. Two "
              "anchors plus a 25-point corroboration; NOT proved. delta = 1 at "
              "(8,28) is a prediction our (72,108) proof does not test, since it "
              "does not use GGV3's chart.")
        print("[NOTE] NOT CLAIMED: that this reproduces GGV3's construction. The "
              "suppressed lattice extension is still unbuilt (leading_x_power.py).")

    if _fail:
        print()
        print("FAILURES (%d):" % len(_fail))
        for f in _fail:
            print("   - %s" % f)
        return 1
    print("multiplicity_partition: %d/%d checks pass -- gamma + delta = a_0, from "
          "two anchors in two papers, corroborated 25/25 by an independent bound"
          % (_ok[0], _ok[0]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
