#!/usr/bin/env python3
"""primitivity_depth.py  (NEW 2026-07-28; read-only)

Gates the claims of PRIMITIVITY_DEPTH.md: GGV3's required-nonzero depth `-10` is
the image of a NEWTON-POLYGON VERTEX under the chart substitution the paper
states, not a fitted constant.

The chain, from 1406.0886_GGV3.tex:1739-1742 (gamma=3) and :1777-1780 (gamma=2):

    x -> x*y^gamma ,   y -> y^(-delta)        (gamma,delta) = (3,2) or (2,3)

so  x^i y^j  ->  x^i y^(gamma*i - delta*j).  The x^0 row of C therefore maps to
y^(-delta*j), and its deepest slot is -delta * deg_y(C at x=0).

WHAT IS GATED HERE
  A  the degeneracy that makes fitting impossible (8 formulas, one value)
  B  N(C) = (1/2)N(P1), cross-checked by 3*N(C) == N(Q1)
  C  the x^0 row's image support == GGV3 (a6)'s printed support, EXACTLY
  D  the step law: step == delta in all THREE published window series
  E  mutation controls -- a wrong delta must break the support match

WHAT IS NOT GATED, because it is not established (PRIMITIVITY_DEPTH.md sec.5):
the depth law has ONE instance ((b5)/(b6) constrain C_-1 and C_1, not C_0);
the leading x-power is unreconciled; (a1)-(a6) remain GGV3's, asserted.

Checker: --quiet, exit 0 iff every check passes.  <1 s.  Reads only.
"""
from __future__ import annotations

import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__)) or "."
sys.path.insert(0, HERE)
QUIET = "--quiet" in sys.argv
_ok = [0]
_fail: list[str] = []

# GGV3's three published window series, transcribed from the .tex.
#   name, gamma, delta, printed y-support
PUBLISHED = [
    ("gamma=3 (a6)  C_0",  3, 2, [2, 0, -2, -4, -6, -8, -10]),
    ("gamma=2 (b6)  C_1",  2, 3, [-1, -4, -7, -10]),
    ("gamma=2 (b5)  C_-1", 2, 3, [1, -2, -5, -8, -11, -14, -17, -20]),
]


def ck(name: str, cond: bool, detail: str = "") -> bool:
    if cond:
        _ok[0] += 1
        if not QUIET:
            print("[OK]   %s" % name)
    else:
        _fail.append(name)
        print("[FAIL] %s%s" % (name, ("  -- " + detail) if detail else ""))
    return bool(cond)


def main() -> int:
    import polygon_reduction as pr

    red = pr.case_f2(0)
    key = "standard (proportional, Prop 8.2(1))"
    NP = red.reduced[key]["P"]
    NQ = red.reduced[key]["Q"]

    # ---- A. the fit is impossible ------------------------------------------
    a0, b0 = 5, 20
    cd = pr.corner_chart_data(a0, b0)
    t, kappa = cd["t"], cd["kappa"]
    degP1 = max(i + j for i, j in NP)
    degC = degP1 // 2
    cands = {
        "-2*a0": -2 * a0, "-2*degC": -2 * degC, "-b0/2": -b0 // 2,
        "-(a0+degC)": -(a0 + degC), "-degP1": -degP1,
        "-t*kappa-2": -(t * kappa + 2), "-(b0-2t-2)": -(b0 - 2 * t - 2),
        "-2t-2": -(2 * t + 2),
    }
    ck("A1  at (5,20) EIGHT distinct corner-data formulas all yield -10, so the "
       "depth cannot be identified by fitting: %s" % sorted(set(cands.values())),
       set(cands.values()) == {-10}, str(cands))
    ck("A2  ... because a0 == degC and b0/2 == 2*a0 == 2*t+2 == t*kappa+2 there",
       a0 == degC and b0 // 2 == 2 * a0 == 2 * t + 2 == t * kappa + 2)
    ck("A3  and only TWO corners have a reduced polygon in-repo, so there is no "
       "second point to fit against",
       len({r.A0 for r in pr.all_reductions()}) == 2,
       str(sorted({r.A0 for r in pr.all_reductions()})))

    # ---- B. N(C) from (a1) P = C^2, cross-checked on Q ----------------------
    NC = [(F(i, 2), F(j, 2)) for i, j in NP]
    ck("B1  N(C) = (1/2)N(P1) is integral -- P = C^2 is consistent with the "
       "computed polygon", all(a.denominator == 1 and b.denominator == 1
                               for a, b in NC), str(NC))
    ck("B2  CROSS-CHECK on the other generator: 3*N(C) == N(Q1) exactly, so the "
       "halving is not an assumption imposed on one polygon",
       sorted((3 * a, 3 * b) for a, b in NC) == sorted((F(i), F(j)) for i, j in NQ))
    jstar = max(b for a, b in NC if a == 0)
    ck("B3  the x^0 row of C runs to the vertex (0,%s), and deg C = %s"
       % (jstar, degC), jstar == 5 and degC == 5)

    # ---- C. the image support IS (a6)'s printed support ---------------------
    delta3 = 2
    predicted = [int(-delta3 * j) for j in range(0, int(jstar) + 1)]
    printed = [e for e in PUBLISHED[0][3] if e <= 0]
    ck("C1  substituting y -> y^-%d sends the x^0 row to %s"
       % (delta3, predicted), predicted == [0, -2, -4, -6, -8, -10])
    ck("C2  ... which is EXACTLY GGV3 (a6)'s printed support at and below y^0 "
       "(%s) -- support and depth both, not one fitted number" % printed,
       predicted == printed)
    ck("C3  so the deepest slot is -delta * deg_y(C|x=0) = %d, and "
       "c_{0,-10} != 0 says the polygon vertex (0,%s) is ATTAINED"
       % (-delta3 * int(jstar), jstar), -delta3 * int(jstar) == -10)

    # ---- D. the step law, an INDEPENDENT prediction -------------------------
    for name, g, d, sup in PUBLISHED:
        steps = sorted({sup[i] - sup[i + 1] for i in range(len(sup) - 1)})
        ck("D  %-20s (gamma=%d): printed step is %s, and delta = %d"
           % (name, g, steps, d), steps == [d], str(sup))
    ck("D4  the three series span TWO charts and TWO distinct deltas, so the "
       "step law is not a restatement of one observation",
       len({d for _, _, d, _ in PUBLISHED}) == 2
       and len({g for _, g, _, _ in PUBLISHED}) == 2)

    # ---- E. mutation controls ------------------------------------------------
    for bad in (1, 3, 4):
        got = [int(-bad * j) for j in range(0, int(jstar) + 1)]
        ck("E  MUTATION delta=%d must NOT reproduce (a6)'s support" % bad,
           got != printed, str(got))
    bad_j = [int(-delta3 * j) for j in range(0, int(jstar))]
    ck("E4  MUTATION -- truncating the x^0 row one vertex short gives %s, which "
       "misses the -10 slot entirely" % bad_j, -10 not in bad_j)

    if _fail:
        print()
        print("FAILURES (%d):" % len(_fail))
        for f in _fail:
            print("   - %s" % f)
        return 1
    print("primitivity_depth: %d/%d checks pass -- GGV3's -10 is the image of the "
          "Newton vertex (0,%s) under y -> y^-%d" % (_ok[0], _ok[0], jstar, delta3))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
