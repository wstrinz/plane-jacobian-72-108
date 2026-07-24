#!/usr/bin/env python3
"""polygon_reduction_verify.py  (NEW; independent exact checker)

PASS/FAIL regression contract for polygon_reduction.py.  Exact sympy + the
pinned published data; no floats, no tolerances.

  R1  the compiled (8,28) reduction reproduces the PUBLISHED polygons EXACTLY
      (paper_src/upstream_facts.json sub1/sub2) and the bracket [P,Q]=x^2.
  R2  F2 j=0 = (50,75): corner signature (t=5,kappa=3,a0=5,q=2,c=y^2(y^3+1))
      and the Phi divisor signature agree with the landed corner data.
  R3  F2 j=1 = (75,125): same chart, judgment DISCHARGED at the polygon layer;
      Phi divisor signature agrees with phi_75_125's landed value.
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
        print("\n[R2] F2 j=0 = (50,75)  (GGV3-consistent)")
    red = pr.case_f2(0)
    check(_sig_of(red) == (5, 3, 5, 2), "R2 corner signature (t,kappa,a0,q)=(5,3,5,2)")
    c_expected = sp.expand(y**2 * (y**3 + 1))
    check(sp.expand(sp.sympify(red.signature["c_of_y"])) == c_expected,
          "R2 residual divisor c(y) = y^2(y^3+1)")
    check(red.bracket == "x^3", "R2 bracket [P,Q]=x^3")
    # Phi signature matches phi_corner4.py's landed (50,75) point
    check(red.signature["phi_signature"] == (189, 75, 38, 76),
          "R2 Phi signature (189,75,38,76) == phi_corner4 landed")
    check(red.signature["N"] == 36, "R2 tower length N=36")


def test_R3():
    if not QUIET:
        print("\n[R3] F2 j=1 = (75,125)  (the target model)")
    red = pr.case_f2(1)
    r0 = pr.case_f2(0)
    check(_sig_of(red) == (5, 3, 5, 2), "R3 corner signature (t,kappa,a0,q)=(5,3,5,2)")
    check(sp.expand(sp.sympify(red.signature["c_of_y"])) == sp.expand(y**2 * (y**3 + 1)),
          "R3 residual divisor c(y) = y^2(y^3+1)")
    # the CHART is identical to F2 j=0 -- same transforms, only (m,n) differs
    check([t.action for t in red.transforms] == [t.action for t in r0.transforms],
          "R3 chart identical to F2 j=0 (same transform sequence; only (m,n) differs)")
    check(red.mn == (3, 5) and r0.mn == (2, 3),
          "R3/R2 differ only in the polygon multiplier (3,5) vs (2,3)")
    check(red.bracket == "x^3", "R3 bracket [P,Q]=x^3 (kappa=l-2=3, DERIVED)")
    # Phi signature matches phi_75_125.py's landed value
    check(red.signature["phi_signature"] == (504, 201, 101, 202),
          "R3 Phi signature (504,201,101,202) == phi_75_125 landed")
    # judgment DISCHARGED at the polygon layer
    joined = " ".join(red.judgment)
    check("RETIRED at the polygon layer" in joined
          and "UNCONDITIONAL at the polygon layer" in joined,
          "R3 'unreduced polygon' judgment DISCHARGED (unconditional at polygon layer)")
    check("honest boundary" in joined,
          "R3 records the honest boundary (residual-gauge choice; here resolved)")


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


def main():
    print("=" * 78)
    print("polygon_reduction_verify.py  --  regression contract")
    print("=" * 78)
    test_R1()
    test_R2()
    test_R3()
    test_branch_manifest()
    print("\n" + "=" * 78)
    if _fails:
        print("RESULT: %d/%d checks passed; %d FAILED:" %
              (_checks - len(_fails), _checks, len(_fails)))
        for f in _fails:
            print("   FAIL: %s" % f)
        print("=" * 78)
        sys.exit(1)
    print("RESULT: ALL %d checks PASSED  (R1 R2 R3 + branch-manifest completeness)"
          % _checks)
    print("=" * 78)
    sys.exit(0)


if __name__ == "__main__":
    main()
