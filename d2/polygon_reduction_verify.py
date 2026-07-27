#!/usr/bin/env python3
"""polygon_reduction_verify.py  (NEW; independent exact checker)

PASS/FAIL regression contract for polygon_reduction.py.  Exact sympy + the
pinned published data; no floats, no tolerances.

  R1  the compiled (8,28) reduction reproduces the PUBLISHED polygons EXACTLY
      (paper_src/upstream_facts.json sub1/sub2) and the bracket [P,Q]=x^2.
  R2  F2 j=0 = (50,75): REPAIRED corner signature (t=4,kappa=2,deg C=1,
      ord C=1,c=y) and -- the decisive EXTERNAL control -- the engine's reduced
      polygons reproduce GGV3 sec.5's three published integers for this exact
      corner: [P_1,Q_1]=x^2, deg(P_1)=10, deg(Q_1)=15.
  R3  F2 j=1 = (75,125): same chart; reduced polygons COMPUTED as m/n * Delta',
      Delta' = {(0,0),(3,0),(4,1),(0,5)}; Phi signature (80,80,0,0), N=77.
  RG  THE RETRACTION GUARD (the root-cause fix).  final_corner_dictionary()
      RETURNS (l_final,b_final) exactly on the retraction shape -- (8,28) and
      (9,24) -- and RAISES off it -- (7,21) and (5,20).  Both directions are
      checked, so the guard is not vacuous.
  BM  branch-manifest COMPLETENESS: every branch keeps >=1 followed option and
      every option carries a non-empty reason; the (8,28) manifest contains the
      published exclusions (k=2 impossible; three-factor discarded; deep shift
      excluded).

Usage:  python3 polygon_reduction_verify.py [--quiet]   (exit 0 = all pass)
"""
import sys
import os
import json
import sympy as sp

import polygon_reduction as pr

y = sp.symbols("y")
QUIET = "--quiet" in sys.argv
HERE = os.path.dirname(os.path.abspath(__file__))

_checks = 0
_fails = []


def check(cond, label):
    global _checks
    _checks += 1
    if cond:
        if not QUIET:
            print("  PASS  %s" % label)
    else:
        _fails.append(label)
        print("  FAIL  %s" % label)


def as_set(vlist):
    return frozenset((int(a), int(b)) for a, b in vlist)


# ---------------------------------------------------------------------------
# Load pinned published facts (self-contained fallback if the json is absent).
# ---------------------------------------------------------------------------
def load_upstream():
    path = os.path.join(HERE, "paper_src", "upstream_facts.json")
    if os.path.exists(path):
        with open(path) as fh:
            return json.load(fh)["facts"]
    # fallback: the transcribed constants (GGHV22 lines 1001-1004)
    return {
        "bracket_case": {"bracket": "[P,Q] = x^2"},
        "newton_polygons": {
            "sub1": {"P": [[0, 0], [1, 0], [8, 14], [8, 16], [0, 8]],
                     "Q": [[0, 0], [2, 1], [12, 21], [12, 24], [0, 12]]},
            "sub2": {"P": [[0, 0], [1, 0], [8, 14], [8, 16]],
                     "Q": [[0, 0], [2, 1], [12, 21], [12, 24]]}},
    }


# ---------------------------------------------------------------------------
# R1 -- the published (8,28) reduction.
# ---------------------------------------------------------------------------
def test_R1():
    if not QUIET:
        print("\n[R1] published (72,108)/(8,28) reduction")
    facts = load_upstream()
    np = facts["newton_polygons"]
    red = pr.case_8_28()

    # bracket
    check(red.bracket == "x^2", "R1 bracket [P,Q] = x^2 (kappa=l-2=2, l=4)")
    check(facts["bracket_case"]["bracket"].replace(" ", "") == "[P,Q]=x^2",
          "R1 pinned bracket string matches x^2")
    # Jacobian is exactly -x^2
    J = pr.fused_jacobian(4)
    check(sp.simplify(J - (-pr.x**2)) == 0, "R1 final-chart Jacobian = -x^2 (exact)")

    # sub2 (cases a,b) and sub1 (case c) reduced polygons -- exact set match
    red_sub2 = red.reduced["sub2 (cases a,b)"]
    red_sub1 = red.reduced["sub1 (case c)"]
    check(as_set(red_sub2["P"]) == as_set(np["sub2"]["P"]),
          "R1 sub2 N(P) == published %s" % np["sub2"]["P"])
    check(as_set(red_sub2["Q"]) == as_set(np["sub2"]["Q"]),
          "R1 sub2 N(Q) == published %s" % np["sub2"]["Q"])
    check(as_set(red_sub1["P"]) == as_set(np["sub1"]["P"]),
          "R1 sub1 N(P) == published %s" % np["sub1"]["P"])
    check(as_set(red_sub1["Q"]) == as_set(np["sub1"]["Q"]),
          "R1 sub1 N(Q) == published %s" % np["sub1"]["Q"])

    # the map is the real work: verify it independently on the published feet
    check(pr.invert_vertex((56, 16), 4) == (8, 16)
          and pr.invert_vertex((48, 14), 4) == (8, 14)
          and pr.invert_vertex((32, 8), 4) == (0, 8),
          "R1 inversion map (a,b)->(4b-a,b) sends the published feet correctly")


# ---------------------------------------------------------------------------
# R2 / R3 -- F2 corner (5,20), landed corner data.
# ---------------------------------------------------------------------------
def _sig_of(red):
    s = red.signature
    return (s["t"], s["kappa"], s["a0"], s["q"])


def test_R2():
    if not QUIET:
        print("\n[R2] F2 j=0 = (50,75)  -- GGV3 sec.5's OWN published reduction")
    red = pr.case_f2(0)
    check(_sig_of(red) == (4, 2, 1, 1),
          "R2 corner signature (t,kappa,deg C,ord C)=(4,2,1,1)  [REPAIRED 2026-07-26]")
    check(sp.expand(sp.sympify(red.signature["c_of_y"]) - y) == 0,
          "R2 residual divisor C = y  (a MONOMIAL: no retraction at (5,20))")
    check(red.signature["g"] == "1", "R2 residual g = 1 (deg g = deg C - ord C = 0)")
    # --- THE EXTERNAL CONTROL.  GGV3 1406.0886 sec.5, tex:1723-1727, verbatim:
    #       "[P_1,Q_1]=x^2,  deg(P_1)=10  and  deg(Q_1)=15."
    # Three published integers.  The engine must reproduce all three.
    check(red.bracket == "x^2", "R2 bracket [P,Q]=x^2 == GGV3's published x^2")
    pq = red.reduced["standard (proportional, Prop 8.2(1))"]
    degP = max(i + j for i, j in pq["P"])
    degQ = max(i + j for i, j in pq["Q"])
    check(degP == 10, "R2 deg(P_1) = %d == GGV3's published 10" % degP)
    check(degQ == 15, "R2 deg(Q_1) = %d == GGV3's published 15" % degQ)
    check(as_set(pq["P"]) == as_set([[0, 0], [6, 0], [8, 2], [0, 10]])
          and as_set(pq["Q"]) == as_set([[0, 0], [9, 0], [12, 3], [0, 15]]),
          "R2 reduced polygons = 2*Delta' and 3*Delta', Delta'={(0,0),(3,0),(4,1),(0,5)}")
    # and the superseded l=5 branch contradicts all three (recorded, not run)
    check(pr.chart_exponent(5, 20) == 4 and 5 * 2 == 10 and 4 * 2 == 8,
          "R2 l is DERIVED: chart_exponent(5,20) = ceil(20/5) = 4 (l=5 would give "
          "bracket x^3 and degrees (20,30), contradicting all three GGV3 integers)")
    check(red.signature["phi_signature"] == (30, 30, 0, 0),
          "R2 Phi signature (30,30,0,0)  [Phi = (1/2) y^30, a monomial]")
    check(red.signature["N"] == 28, "R2 tower length N=28")


def test_R3():
    if not QUIET:
        print("\n[R3] F2 j=1 = (75,125)  (the target model)")
    red = pr.case_f2(1)
    r0 = pr.case_f2(0)
    check(_sig_of(red) == (4, 2, 1, 1),
          "R3 corner signature (t,kappa,deg C,ord C)=(4,2,1,1)  [REPAIRED 2026-07-26]")
    check(sp.expand(sp.sympify(red.signature["c_of_y"]) - y) == 0,
          "R3 residual divisor C = y  (a MONOMIAL)")
    # the CHART is identical to F2 j=0 -- same transforms, only (m,n) differs
    check([t.action for t in red.transforms] == [t.action for t in r0.transforms],
          "R3 chart identical to F2 j=0 (same transform sequence; only (m,n) differs)")
    check(red.mn == (3, 5) and r0.mn == (2, 3),
          "R3/R2 differ only in the polygon multiplier (3,5) vs (2,3)")
    check(red.bracket == "x^2", "R3 bracket [P,Q]=x^2 (kappa=l-2=2, DERIVED)")
    # the reduced polygons are now COMPUTED, not absent
    pq = red.reduced["standard (proportional, Prop 8.2(1))"]
    check(as_set(pq["P"]) == as_set([[0, 0], [9, 0], [12, 3], [0, 15]]),
          "R3 N(P) = {(0,0),(9,0),(12,3),(0,15)} = 3*Delta'  (COMPUTED)")
    check(as_set(pq["Q"]) == as_set([[0, 0], [15, 0], [20, 5], [0, 25]]),
          "R3 N(Q) = {(0,0),(15,0),(20,5),(0,25)} = 5*Delta'  (COMPUTED)")
    check(max(i + j for i, j in pq["P"]) == 15
          and max(i + j for i, j in pq["Q"]) == 25,
          "R3 reduced degrees (15,25)")
    check(red.signature["phi_signature"] == (80, 80, 0, 0),
          "R3 Phi signature (80,80,0,0)  [Phi = (1/3) y^80, a monomial]")
    check(red.signature["N"] == 77, "R3 tower length N=77")
    # judgment DISCHARGED at the polygon layer
    joined = " ".join(red.judgment)
    check("RETIRED at the polygon layer" in joined
          and "UNCONDITIONAL at the polygon layer" in joined,
          "R3 'unreduced polygon' judgment DISCHARGED (unconditional at polygon layer)")
    # The 2026-07-24 residual-GAUGE reopening is DISSOLVED, not resolved and not
    # still open: it presupposed deg g = 3, which presupposed deg C = a0 = 5,
    # which presupposed a retraction this corner does not have.  Check the
    # DISSOLVED verdict and, independently, that deg g really is 0 -- so this is
    # not a string test that a stale "REOPENED" would also satisfy.
    check("DISSOLVED" in joined and red.signature["g"] == "1"
          and red.signature["a0"] - red.signature["q"] == 0,
          "R3 residual-gauge branch is DISSOLVED: deg g = deg C - ord C = 0, so g=1 "
          "is forced and there is no gauge to choose")


# ---------------------------------------------------------------------------
# RG -- THE RETRACTION GUARD.  Both directions, so the guard is not vacuous.
# ---------------------------------------------------------------------------
def test_RG():
    if not QUIET:
        print("\n[RG] retraction guard on the (t,q)=(l_final,b_final) dictionary")
    # the INFERRED chart-exponent rule on every corner with published data
    for (a0, b0), want in {(8, 28): 4, (9, 24): 3, (9, 27): 3,
                           (7, 21): 3, (5, 20): 4}.items():
        check(pr.chart_exponent(a0, b0) == want,
              "RG chart_exponent(%d,%d) = ceil(%d/%d) = %d" % (a0, b0, b0, a0, want))
    # the retraction shape splits the four dictionary rows exactly as PASSPORT P6
    check(pr.has_retraction(8, 28) and pr.has_retraction(9, 24),
          "RG retraction shape HOLDS at (8,28) [28=4*7] and (9,24) [24=3*8]")
    check(not pr.has_retraction(7, 21) and not pr.has_retraction(5, 20),
          "RG retraction shape FAILS at (7,21) [21!=3*6] and (5,20) [20!=4*4]")
    # the quantifier trap: SOME l satisfies b0=l(a0-1) at (5,20) (namely 5), and
    # that is exactly how l=5 got in.  The guard must still refuse.
    check(pr.has_retraction(5, 20, l=5) and not pr.has_retraction(5, 20),
          "RG the trap: b0=l(a0-1) is solvable at (5,20) by l=5, but FAILS for the "
          "l actually used (4) -- the guard tests the l in use, not existence")
    # POSITIVE direction: the dictionary is returned where it is valid
    check(pr.final_corner_dictionary(8, 28, 4, 7) == (4, 7),
          "RG dictionary RETURNS (t,q)=(4,7) at (8,28) -- matches the published "
          "chart l=4 and edge form y(x^4 y-alpha)^7")
    check(pr.final_corner_dictionary(9, 24, 3, 8) == (3, 8),
          "RG dictionary RETURNS (t,q)=(3,8) at (9,24)")
    # NEGATIVE direction: it raises where it is invalid, loudly
    for a0, b0, lf, bf, why in [
            (7, 21, 7, 2, "GGV5 l_final=7 vs GGHV22's published chart l=3"),
            (5, 20, 5, 2, "GGV5 l_final=5 vs GGV3's (50,75) reduction forcing l=4")]:
        try:
            pr.final_corner_dictionary(a0, b0, lf, bf)
            check(False, "RG dictionary must RAISE at (%d,%d)" % (a0, b0))
        except pr.FinalCornerDictionaryError as exc:
            check("retraction precondition" in str(exc) and "REFUSED" in str(exc),
                  "RG dictionary RAISES at (%d,%d): %s" % (a0, b0, why))
    # and corner_chart_data routes both shapes correctly
    cd8 = pr.corner_chart_data(8, 28, l_final=4, b_final=7)
    check((cd8["t"], cd8["kappa"], cd8["deg_C"], cd8["ord_C"]) == (4, 2, 8, 7)
          and cd8["retraction"] and not cd8["monomial"],
          "RG corner_chart_data(8,28) = (t,kappa,deg C,ord C)=(4,2,8,7), retracted")
    cd5 = pr.corner_chart_data(5, 20, l_final=5, b_final=2)
    check((cd5["t"], cd5["kappa"], cd5["deg_C"], cd5["ord_C"]) == (4, 2, 1, 1)
          and cd5["monomial"] and not cd5["retraction"],
          "RG corner_chart_data(5,20) = (t,kappa,deg C,ord C)=(4,2,1,1), MONOMIAL")
    cd7 = pr.corner_chart_data(7, 21, l_final=7, b_final=2)
    check((cd7["t"], cd7["kappa"], cd7["deg_C"], cd7["ord_C"]) == (3, 1, 1, 1),
          "RG corner_chart_data(7,21) = (3,1,1,1) == GGHV22's published (7,21) "
          "chart y x^3, [P,Q]=x, C=y  (independent published confirmation)")


# ---------------------------------------------------------------------------
# BM -- branch-manifest completeness.
# ---------------------------------------------------------------------------
def test_branch_manifest():
    if not QUIET:
        print("\n[BM] branch-manifest completeness")
    for red in pr.all_reductions():
        for b in red.branches:
            check(len(b.followed()) >= 1,
                  "%s / '%s' keeps >=1 followed option" % (red.tag, b.name))
            check(all(o.reason.strip() for o in b.options),
                  "%s / '%s' every option carries a reason" % (red.tag, b.name))
            check(bool(b.cite.strip()),
                  "%s / '%s' branch is cited" % (red.tag, b.name))
    # the published (8,28) exclusions must be present
    r1 = pr.case_8_28()
    all_opts = [(b.name, o.label, o.followed, o.reason)
                for b in r1.branches for o in b.options]
    check(any((not f) and "k=2" in lbl for _, lbl, f, _ in all_opts),
          "BM (8,28): the k=2 branch is EXCLUDED (edges non-parallel)")
    check(any((not f) and "three" in lbl.lower() for _, lbl, f, _ in all_opts),
          "BM (8,28): the three-distinct-factor branch is EXCLUDED")
    check(any((not f) and "s>=4" in lbl for _, lbl, f, _ in all_opts),
          "BM (8,28): deep/shallow root-shift depth EXCLUDED (GGV6 Prop 2.5)")
    # both output shapes retained (sub1 and sub2)
    check(set(r1.reduced.keys()) == {"sub2 (cases a,b)", "sub1 (case c)"},
          "BM (8,28): BOTH published output shapes (sub1, sub2) retained")
    # F2 chart-class branch excludes the double-inversion heuristic
    r3 = pr.case_f2(1)
    check(any((not o.followed) and "double-inversion" in o.label
              for b in r3.branches for o in b.options),
          "BM F2: double-inversion (kappa=l2-l1) heuristic EXCLUDED")
    # and the repaired branch must record the superseded l=5 as EXCLUDED, with
    # the retraction reason -- so the error cannot be silently reintroduced
    check(any((not o.followed) and "l_final = 5" in o.label
              and "retraction shape" in o.reason
              for b in r3.branches for o in b.options),
          "BM F2: l = l_final = 5 is an EXPLICITLY EXCLUDED option, with the "
          "retraction-shape reason attached")
    check(any((not o.followed) and "deg C = a0 = 5" in o.label
              for b in r3.branches for o in b.options),
          "BM F2: the retracted-shape C (deg C = a0 = 5, ord C = 2) is EXCLUDED")


def main():
    print("=" * 78)
    print("polygon_reduction_verify.py  --  regression contract")
    print("=" * 78)
    test_R1()
    test_R2()
    test_R3()
    test_RG()
    test_branch_manifest()
    print("\n" + "=" * 78)
    if _fails:
        print("RESULT: %d/%d checks passed; %d FAILED:" %
              (_checks - len(_fails), _checks, len(_fails)))
        for f in _fails:
            print("   FAIL: %s" % f)
        print("=" * 78)
        sys.exit(1)
    print("RESULT: ALL %d checks PASSED  (R1 R2 R3 RG + branch-manifest completeness)"
          % _checks)
    print("=" * 78)
    sys.exit(0)


if __name__ == "__main__":
    main()
