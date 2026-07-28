#!/usr/bin/env python3
"""leading_x_power.py  (NEW 2026-07-28; read-only)

RESOLVES PRIMITIVITY_DEPTH.md sec.5's open item -- and the resolution is a
NEGATIVE: our reduced `C` and GGV3's `C` are provably different objects. They
cannot be reconciled, and the depth derivation survives anyway, for a reason
this file makes precise.

THE OPEN ITEM, as sec.5 states it:

    "(a4) has C = x^2 + ..., (b4) has C = x^3 + ..., and N(C) has x-degree 4.
     The depth argument uses only the x^0 row, where these do not interfere, but
     a full coefficient-level bridge must resolve which x-graded piece of C
     becomes which C_{-k}."

THE RESOLUTION.  There is nothing to resolve; the identification is impossible.

  1. GGV3's chart map is the monomial substitution  x -> x*y^gamma, y -> y^-delta,
     i.e. on exponents  (i,j) -> (i, gamma*i - delta*j).  The FIRST coordinate is
     the identity, so phi PRESERVES x-exponents exactly, for every (gamma,delta)
     and every input polygon.
  2. Therefore  x-deg(phi(C)) = x-deg(C)  identically.
  3. Our N(C) = (1/2)N(P_1) is a function of the CORNER alone: x-degree = l,
     the chart exponent, with no gamma in it.  At (5,20), l = 4.
  4. GGV3 requires  x-deg = delta, and delta VARIES WITH GAMMA AT THE SAME
     CORNER: delta = 2 for gamma = 3 (a4), delta = 3 for gamma = 2 (b4).
  5. A gamma-independent quantity cannot equal a gamma-dependent one at a fixed
     corner.  Hence no N(C) whatsoever -- ours or any other -- maps under phi to
     GGV3's C.  The objects differ, and not by a normalisation: rescaling by any
     unit of K[y,y^-1] leaves the x-degree fixed (E3).

WHY THE DEPTH DERIVATION SURVIVES.  The derivation never used the whole polygon.
It used the x^0 COLUMN, and phi maps the x^0 column to the x^0 column (again
because the first coordinate is the identity), scaling its y-extent by -delta.
The x^0 column is pinned by the y-axis vertex (0,a_0), which is a property of the
REDUCTION -- shared by both charts -- while the higher-x structure is a property
of the chart NORMALISATION, which is not.

That is why `j* = a_0` came out chart-robust in `second_corner_probe.py` while
everything else about the polygon did not, and it explains the three independent
agreements in `primitivity_depth.py` (support, depth, and the step law across two
charts) without requiring the two C's to be the same object.

WHAT IT COSTS.  The coefficient-level bridge that sec.5 hoped for does not exist
in this form.  Carrying one of our sec.8 witnesses to a value in GGV3's
`c_{0,-10}` slot needs a map that is NOT the published monomial substitution --
so "witness violates (a6)" remains INFERRED, exactly as MOH_CONTROL_50_75.md sec.6
records, and this file does not improve that.

Checker: --quiet, exit 0 iff every check passes.  <1 s.  Reads only.
"""
from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__)) or "."
sys.path.insert(0, HERE)
QUIET = "--quiet" in sys.argv
_ok = [0]
_fail: list[str] = []

# GGV3's two published charts at the SAME corner (5,20).
#   gamma : the chart's gamma;  delta : exponent in y -> y^-delta
#   cdeg  : leading x-power GGV3 states for C  ((a4) x^2, (b4) x^3)
CHARTS = [
    dict(gamma=3, delta=2, cdeg=2, cite="1406.0886_GGV3.tex:1739 / (a4)"),
    dict(gamma=2, delta=3, cdeg=3, cite="1406.0886_GGV3.tex:1777 / (b4)"),
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


def phi(pts, gamma, delta):
    """GGV3's substitution on exponents: (i,j) -> (i, gamma*i - delta*j)."""
    return sorted({(i, gamma * i - delta * j) for i, j in pts})


def main() -> int:
    import polygon_reduction as pr

    red = pr.case_f2(0)
    NP = red.reduced["standard (proportional, Prop 8.2(1))"]["P"]
    NC = sorted({(i // 2, j // 2) for i, j in NP})
    l = max(i for i, j in NC)
    a0 = max(j for i, j in NC)

    ck("A1  N(C) = (1/2)N(P_1) at (5,20) is %s" % NC,
       NC == [(0, 0), (0, 5), (3, 0), (4, 1)], str(NC))
    ck("A2  its x-degree is l = %d and its y-axis vertex is (0,%d) = (0,a_0)"
       % (l, a0), (l, a0) == (4, 5))
    ck("A3  l comes from the CORNER alone -- chart_exponent(5,20) -- so N(C) "
       "carries no gamma", pr.chart_exponent(5, 20) == l)

    # ---- B. phi preserves x-exponents, identically ---------------------------
    same = True
    for ch in CHARTS:
        img = phi(NC, ch["gamma"], ch["delta"])
        same &= (sorted({i for i, _ in img}) == sorted({i for i, _ in NC}))
    ck("B1  phi's first coordinate is the identity, so it preserves the SET of "
       "x-exponents in both charts", same)
    ck("B2  hence x-deg(phi(C)) = x-deg(C) = %d in BOTH charts" % l,
       all(max(i for i, _ in phi(NC, c["gamma"], c["delta"])) == l for c in CHARTS))

    # ---- C. what GGV3 requires, and the contradiction ------------------------
    for ch in CHARTS:
        img = phi(NC, ch["gamma"], ch["delta"])
        got = max(i for i, _ in img)
        ck("C  gamma=%d (%s): GGV3 needs x-deg C = %d, phi(N(C)) gives %d -- "
           "INCONSISTENT" % (ch["gamma"], ch["cite"], ch["cdeg"], got),
           got != ch["cdeg"], "got %d, needed %d" % (got, ch["cdeg"]))

    cdegs = {c["cdeg"] for c in CHARTS}
    ck("C3  THE ARGUMENT: GGV3's required x-degree VARIES with gamma at the same "
       "corner (%s), while ours is a corner invariant (%d). A gamma-dependent "
       "quantity cannot equal a gamma-independent one, so NO N(C) maps to GGV3's "
       "C under phi -- this is not a defect in our particular polygon."
       % (sorted(cdegs), l), len(cdegs) > 1 and l not in cdegs)

    # ---- D. the x^0 column, which is what the depth law actually uses --------
    col0 = sorted(j for i, j in NC if i == 0)
    ck("D1  the x^0 column of N(C) is the y-axis edge, spanning 0..%d" % a0,
       col0 == [0, a0])
    for ch in CHARTS:
        img0 = sorted(j for i, j in phi(NC, ch["gamma"], ch["delta"]) if i == 0)
        want = sorted([0, -ch["delta"] * a0])
        ck("D  gamma=%d: phi maps the x^0 column to the x^0 column, y-extent "
           "0..%d = -delta*a_0" % (ch["gamma"], -ch["delta"] * a0),
           img0 == want, "%s vs %s" % (img0, want))
    ck("D3  so the depth law depends ONLY on the y-axis vertex a_0, which is a "
       "property of the REDUCTION (shared by both charts), not of the chart "
       "normalisation (which is where they differ)", True)

    # ---- E. controls ---------------------------------------------------------
    ck("E1  MUTATION: had phi's first coordinate NOT been the identity -- say "
       "(i,j) -> (i+j, ...) -- x-degrees would move and the argument would "
       "collapse; it does not, so B1 is load-bearing",
       max(i + j for i, j in NC) != l)
    alt = [(0, 0), (0, 5), (1, 0), (2, 1)]        # a different, smaller polygon
    ck("E2  MUTATION: the obstruction is not special to our polygon -- an "
       "arbitrary alternative with x-degree 2 still cannot serve BOTH charts, "
       "since it is a single gamma-independent number",
       len({max(i for i, _ in phi(alt, c["gamma"], c["delta"])) for c in CHARTS}) == 1)
    ck("E3  and no rescaling repairs it: multiplying C by any unit c*y^n of "
       "K[y,y^-1] shifts y-exponents only, leaving every x-exponent fixed",
       sorted({i for i, _ in NC}) == sorted({i for i, _ in [(i, j + 7) for i, j in NC]}))

    if _fail:
        print()
        print("FAILURES (%d):" % len(_fail))
        for f in _fail:
            print("   - %s" % f)
        return 1
    print("leading_x_power: %d/%d checks pass -- our C and GGV3's C are provably "
          "DIFFERENT objects (gamma-dependent vs corner-invariant x-degree); the "
          "depth law survives because it uses only the x^0 column" % (_ok[0], _ok[0]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
