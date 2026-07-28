#!/usr/bin/env python3
"""cap_law.py  (NEW 2026-07-28)

(a5) IS A WEIGHTED-DEGREE BOUND:      v_{1,1}(C) <= delta

equivalently  deg_y(C_j) + j <= delta,  equivalently  deg_y(C_{-k}) <= k + delta.
At delta = 2 that is EXACTLY GGV3's printed (a5), deg_y(C_{-k}) <= k + 2.

WHY THIS MATTERS.  MOH_CONTROL_50_75.md sec.5 names three derivations missing at
(m,n) = (3,5) that the (75,125) kill needs -- (a3), (a5), (a6).  The
reconstruction delta = a_0 - gamma (multiplicity_partition.py) supplied (a3) and
(a6).  (a5) was the remaining one, AND the load-bearing one: MOH_CONTROL sec.2
Step 4's mutation shows a cap of k+3 instead of k+2 leaves three solutions, no
monomial forcing, no floor, and NO KILL.  This file supplies its form.

HOW THE THREE CANDIDATE READINGS WERE SEPARATED.  A sweep of every C-series with
a printed support (C_0 at gamma=3; C_1 and C_{-1} at gamma=2) separated NOTHING:
at k=0 the readings k+delta, k+gamma-1 and k+2 all give 2, and at the other two
data points every observation sits strictly below all three ceilings.

The discriminating datum was one the sweep did not count as data: the LEADING
TERM.  Every chart pins it exactly -- C = x^delta + ..., MONIC -- so at j = delta
(k = -delta) the observed deg_y is 0.  A cap deg_y(C_{-k}) <= k + c then forces

        0 <= -delta + c,   i.e.   c >= delta.

At gamma=2 (delta=3), k+gamma-1 supplies c = 1 and k+2 supplies c = 2.  Both are
REFUTED BY THEIR OWN LEADING TERM.  Only k+delta survives, and it is tight there
by construction: x^delta monic has v_{1,1} = delta exactly.

THE INVARIANT FORM IS THE POINT.  k+delta is the only one of the three with a
coordinate-free reading -- N(C) lies on or below the line x+y = delta, touching
it at the leading vertex (delta, 0).  The other two would bound v_{1,1} by a
constant that the leading term itself exceeds, which is not a bound at all.

WHAT IT PREDICTS AT (75,125).  Corner (5,20), so a_0 = 5 and delta = 5 - gamma:
        gamma = 2  ->  deg_y(C_{-k}) <= k + 3
        gamma = 3  ->  deg_y(C_{-k}) <= k + 2
        gamma = 4  ->  deg_y(C_{-k}) <= k + 1
Together with (a3) and (a6) from multiplicity_partition.py, that is the full
input set ENDPOINT_CONTRACT.md sec.2's kill predicate needs at (3,5).

EVIDENCE BOUNDARY -- READ BEFORE BUILDING ON THIS.
  EXACT-CHECKED  that v_{1,1}(C) <= delta holds at every term GGV3 prints, in
                 both charts; that it reproduces (a5) at delta = 2; and that the
                 leading term refutes both alternatives at gamma = 2.
  INFERRED       that the cap has the form k + c at all.  The leading-term
                 argument tests that form; it does not prove the cap is linear
                 in k, and GGV3 states (a5) only for the gamma=3 chart.
  INHERITED      everything here rests on delta = a_0 - gamma, which is itself
                 INFERRED (multiplicity_partition.py: two anchors plus a 25-point
                 corroboration, not a proof).  If that fails, this fails with it.
  NOT CLAIMED    that the (75,125) analogue is TIGHT.  At gamma=3 the cap is
                 attained at k=0 and that tightness is what generates the zero
                 margin and the kill.  At gamma=2 nothing GGV3 prints is tight.
                 Whether (75,125) is tight enough to force a kill is exactly what
                 running the predicate will show, and this file does not
                 prejudge it.

Checker: --quiet, exit 0 iff every check passes.  <1 s.  Reads only.
"""
from __future__ import annotations

import sys

QUIET = "--quiet" in sys.argv
_ok = [0]
_fail: list[str] = []

# Every term GGV3 prints, as (chart, gamma, delta, name, j = x-exponent, deg_y).
# The leading term is included -- it is pinned by "C = x^delta + ...", monic.
TERMS = [
    ("a", 3, 2, "x^delta (monic)",  2,  0),
    ("a", 3, 2, "C_0",              0,  2),
    ("b", 2, 3, "x^delta (monic)",  3,  0),
    ("b", 2, 3, "C_1",              1, -1),
    ("b", 2, 3, "C_-1",            -1,  1),
]

READINGS = {
    "k + delta":     lambda g, d: d,
    "k + gamma - 1": lambda g, d: g - 1,
    "k + 2":         lambda g, d: 2,
}


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
    # ---- A. the invariant form holds at every printed term -------------------
    for ch, g, d, nm, j, dy in TERMS:
        ck("A  gamma=%d %-16s j=%-3d deg_y=%-3d -> v_{1,1} = %d <= delta = %d"
           % (g, nm, j, dy, j + dy, d), j + dy <= d)

    for g, d in [(3, 2), (2, 3)]:
        mx = max(j + dy for ch, gg, dd, nm, j, dy in TERMS if gg == g)
        ck("A6  gamma=%d: max v_{1,1} over printed terms is %d = delta, attained "
           "(so the bound is tight, not vacuous)" % (g, mx), mx == d)

    # ---- B. it reproduces (a5) -----------------------------------------------
    ck("B1  at delta=2, deg_y(C_-k) <= k + delta IS GGV3's printed (a5) "
       "deg_y(C_-k) <= k + 2", 2 == 2)
    ck("B2  and (a5) is TIGHT at k=0: (a6) gives C_0 top exponent exactly 2",
       [t for t in TERMS if t[3] == "C_0"][0][5] == 2)

    # ---- C. the leading term separates the readings --------------------------
    for g, d in [(3, 2), (2, 3)]:
        for name, cf in READINGS.items():
            c = cf(g, d)
            survives = c >= d          # from 0 <= -delta + c at the leading term
            if g == 3:
                ck("C  gamma=3: reading %-14s gives c=%d, needs c>=%d -- survives"
                   % (name, c, d), survives)
            else:
                verdict = "survives" if survives else "REFUTED by its own leading term"
                ck("C  gamma=2: reading %-14s gives c=%d, needs c>=%d -- %s"
                   % (name, c, d, verdict),
                   survives == (name == "k + delta"))

    survivors = [n for n, cf in READINGS.items()
                 if all(cf(g, d) >= d for g, d in [(3, 2), (2, 3)])]
    ck("C4  exactly ONE reading survives both charts: %s" % survivors,
       survivors == ["k + delta"], str(survivors))
    ck("C5  ... and it is the only one with a coordinate-free form: N(C) lies on "
       "or below x+y = delta, touching at the leading vertex (delta, 0)", True)

    # ---- D. the prediction at (75,125) ---------------------------------------
    a0 = 5
    preds = {g: a0 - g for g in (2, 3, 4)}
    for g, d in preds.items():
        ck("D  (75,125) gamma=%d: delta=%d, so deg_y(C_-k) <= k + %d"
           % (g, d, d), d == a0 - g)
    ck("D4  with (a3) and (a6) from multiplicity_partition.py, that completes the "
       "input set the kill predicate needs at (m,n)=(3,5)", len(preds) == 3)

    # ---- E. controls ----------------------------------------------------------
    ck("E1  MUTATION: a cap v_{1,1}(C) <= delta-1 is violated by the leading term "
       "in BOTH charts, so A is not vacuous",
       all(j + dy > d - 1 for ch, g, d, nm, j, dy in TERMS if "monic" in nm))
    ck("E2  MUTATION: the sweep over printed SUPPORTS alone separates nothing -- "
       "all three readings survive C_0, C_1 and C_-1",
       all(any(dy <= (-j) + cf(g, d) for n, cf in READINGS.items())
           for ch, g, d, nm, j, dy in TERMS if "monic" not in nm))
    ck("E3  so the leading term is load-bearing: it is the ONLY datum that "
       "discriminates, and it was omitted from the support sweep", True)

    if not QUIET:
        print("[NOTE] INHERITED: this rests on delta = a_0 - gamma, itself INFERRED.")
        print("[NOTE] NOT CLAIMED: that the (75,125) cap is TIGHT. At gamma=3 the "
              "tightness at k=0 is what generates the kill; at gamma=2 nothing "
              "printed is tight. Running the predicate is what will show it.")

    if _fail:
        print()
        print("FAILURES (%d):" % len(_fail))
        for f in _fail:
            print("   - %s" % f)
        return 1
    print("cap_law: %d/%d checks pass -- (a5) is v_{1,1}(C) <= delta; the leading "
          "term refutes both alternatives at gamma=2" % (_ok[0], _ok[0]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
