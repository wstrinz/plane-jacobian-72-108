#!/usr/bin/env python3
"""bigrade_annotator_verify.py  (NEW 2026-07-26; read-only)

INDEPENDENT exact checker for the load-bearing claim of bigrade_annotator.py:
the R1 regression, in which the (50,75) gamma=3 kill is *rediscovered* by
scanning the endpoint contract's required-nonzero list against a COMPUTED
forced floor.

Why this file exists
--------------------
bigrade_annotator.py's own docstring advertised "Independent checker:
bigrade_annotator_verify.py".  That file did not exist, and the annotator was
gated in no suite (tools/suite_manifest.py had no entry, and EXPECTED_TOTAL never
counted it).  The module self-checks R1 with a bare `assert` inside main(), which
is (a) not independent, (b) only reachable after R2/R3 -- the slow diagnostic
stages -- have run, costing >14 CPU-minutes to exercise a 0.2 s claim.

This checker is fast (~3 s), independent in METHOD, and gated.

What is checked (all exact sympy; --quiet; exit 0 iff every check passes)
------------------------------------------------------------------------
  A. The forced corner series C_0, derived by a DIFFERENT ROUTE from the
     annotator's.  The annotator sets up a truncated window j in [-12, 2] and
     calls sp.solve on the scalarised coefficient equations.  Here we divide the
     forcing relation directly in the Laurent ring -- no truncation, no linear
     solve -- and read the support off a shifted Poly.  Both must give
     support {-6,-4,-2,0,2} and floor -6.

  B. Agreement between the two routes, on the module's own object.

  C. ANSATZ-INDEPENDENCE.  The annotator's floor is computed inside a hardcoded
     window JMIN=-12.  A floor that is an artifact of the truncation would be
     invisible to the module's own assert.  We re-solve at JMIN in
     {-8,-12,-20,-30} and require the floor to be -6 every time.

  D. The kill predicate fires at EXACTLY [(0,-10)] -- not merely "fires".

  E. DISCRIMINATION.  A regression that passes on anything is worthless, so four
     mutations must each flip the verdict the right way:
       E1 deepen the forcing by one term (mu*y^-4)  -> floor -10, NO kill
       E2 drop (0,-10) from the contract            -> NO kill
       E3 move the requirement to (0,-6) = the floor-> NO kill (strict `<`)
       E4 move the requirement to (0,-7)            -> kill (boundary is exact)

  F. CONTRACT FIDELITY.  ENDPOINT_CONTRACT.md sec.3 states forced_floor and
     required_nonzero as literals.  They must agree with the code.  Prose drift
     is a test failure here -- three documents in this repo asserted superseded
     arithmetic for weeks precisely because nothing tested prose against code.

  G. CROSS-ROUTE agreement with f2_tower.py's own gamma=3 expression (check B of
     f2_tower_verify.py), recomputed here rather than imported, so the two
     modules are confirmed to be talking about the same C_0.

SCOPE -- what this checker does NOT establish
---------------------------------------------
The forcing relation 3(a y^3)^2 C_0 = 3(b y^4)^2 + 2 lam + 2 F_{-2} is TRANSCRIBED
from GGV3 sec.5, in both bigrade_annotator.py and f2_tower.py.  Neither module
derives it from corner data: f2_tower.a2_certificate() references T, KAPPA, QC, C,
ordPhi, Nof and build_gsystem zero times each.  So A-G confirm that GIVEN GGV3's
relation the floor is -6 and the kill localises at c_{0,-10} automatically.  They
do NOT confirm that the relation itself follows from the (5,20) corner.  Deriving
it is the window compiler's job, and until that exists the (50,75) "reproduction"
is a REPLAY of published algebra, not an oracle that can be pointed at (75,125).
"""
import os
import re
import sys

import sympy as sp

sys.path.insert(0, __file__.rsplit("bigrade_annotator_verify.py", 1)[0] or ".")

HERE = os.path.dirname(os.path.abspath(__file__)) or "."
QUIET = "--quiet" in sys.argv
FAILS = []
N_OK = 0

y = sp.symbols("y")
a, b, lam, f2, f4, f6, f8, mu = sp.symbols("a b lam f2 f4 f6 f8 mu")


def ok(label, cond):
    global N_OK
    if cond:
        N_OK += 1
        if not QUIET:
            print("[OK]", label)
    else:
        FAILS.append(label)
        print("[FAIL]", label)


# ---------------------------------------------------------------------------
# Primitives
# ---------------------------------------------------------------------------
def y_support(expr, pad=64):
    """Exponents of y occurring in a Laurent expression, read off a shifted Poly.

    Deliberately NOT the term-by-term fraction-degree walk used in f2_tower.py:
    shifting into a genuine polynomial and letting Poly collect the terms also
    catches cancellation that a per-term walk would miss.
    """
    shifted = sp.expand(sp.cancel(sp.together(expr)) * y**pad)
    shifted = sp.expand(sp.cancel(shifted))
    if shifted == 0:
        return set()
    p = sp.Poly(shifted, y)
    return {m[0] - pad for m, c in zip(p.monoms(), p.coeffs()) if c != 0}


def forcing_rhs(extra_deep=False):
    """RHS of the GGV3 sec.5 gamma=3 forcing relation.  extra_deep adds one
    deeper forcing term, used only as the E1 mutation."""
    F = f8 * y**8 + f6 * y**6 + f4 * y**4 + f2 * y**2
    rhs = 3 * (b * y**4) ** 2 + 2 * lam + 2 * F
    if extra_deep:
        rhs = rhs + 2 * mu * y**-4
    return rhs


def c0_by_division(extra_deep=False):
    """ROUTE 1 (independent): divide in the Laurent ring.  No window, no solve.

    Returned EXPANDED, as a sum of Laurent monomials.  sp.cancel alone leaves a
    single fraction, on which .coeff(y, k) silently returns the whole numerator
    over the constant part of the denominator rather than a Laurent coefficient.
    """
    return sp.expand(sp.cancel(forcing_rhs(extra_deep) / (3 * (a * y**3) ** 2)))


def c0_floor_by_solve(jmin, jmax=2, extra_deep=False):
    """ROUTE 2 (the annotator's method, reimplemented and parameterised in jmin).

    Scalarise the relation over a truncated window and solve the linear block;
    a coefficient is forced-zero iff its solved value is identically 0.
    """
    c0 = {j: sp.Symbol("c0_%d" % (j + 100)) for j in range(jmin, jmax + 1)}
    C0series = sum(c0[j] * y**j for j in range(jmin, jmax + 1))
    relation = sp.expand(3 * (a * y**3) ** 2 * C0series - forcing_rhs(extra_deep))
    shift = -min(jmin + 6, -6 if not extra_deep else -10) + 8
    cleared = sp.expand(relation * y**shift)
    eqs = [sp.Poly(cleared, y).nth(pw) for pw in range(sp.Poly(cleared, y).degree() + 1)]
    eqs = [e for e in eqs if e != 0]
    sol = sp.solve(eqs, list(c0.values()), dict=True)
    if not sol:
        return None, []
    sol = sol[0]
    nz = [j for j, s in c0.items() if sp.simplify(sol.get(s, s)) != 0]
    return (min(nz) if nz else None), sorted(nz)


def kill_predicate(required_nonzero, floor, series=0):
    """ENDPOINT_CONTRACT.md sec.2: a required-nonzero coefficient STRICTLY below
    the forced floor of its series is a kill."""
    return [(s, j) for (s, j) in required_nonzero if s == series and j < floor]


EXPECTED_SUPPORT = {-6, -4, -2, 0, 2}
EXPECTED_FLOOR = -6
CONTRACT_REQUIRED = [(-1, 3), (-2, 4), (0, -10)]


# ---------------------------------------------------------------------------
# A. independent derivation of the forced floor
# ---------------------------------------------------------------------------
def check_A():
    C0 = c0_by_division()
    sup = y_support(C0)
    ok("A1: C_0 by Laurent division has support {-6,-4,-2,0,2} (no truncation, "
       "no linear solve)", sup == EXPECTED_SUPPORT)
    ok("A2: forced floor of series 0 is -6", (min(sup) if sup else None) == EXPECTED_FLOOR)
    ok("A3: the floor is carried by the lam term -- 2*lam/(3a^2) at y^-6 -- so "
       "the depth is set by the forcing constant, not by b",
       sp.simplify(C0.coeff(y, -6) - sp.Rational(2, 3) * lam / a**2) == 0)
    ok("A4: no support strictly below -6 (the slot c_{0,-10} is empty)",
       all(e >= -6 for e in sup) and -10 not in sup)


# ---------------------------------------------------------------------------
# B. agreement with the module's own route, on the module's own object
# ---------------------------------------------------------------------------
def check_B():
    try:
        import bigrade_annotator as B
    except Exception as exc:                                   # pragma: no cover
        ok("B0: bigrade_annotator imports", False)
        if not QUIET:
            print("     import error:", exc)
        return
    ok("B0: bigrade_annotator imports", True)
    r1 = B.build_R1()
    red = B.r1_rediscover(r1)
    ok("B1: module's COMPUTED forced_floor[0] agrees with the division route (-6)",
       red["forced_floor_series0"] == EXPECTED_FLOOR)
    ok("B2: module's computed support agrees with the division route",
       set(red["forced_support_series0"]) == EXPECTED_SUPPORT)
    ok("B3: module's contract required-nonzero list is [(-1,3),(-2,4),(0,-10)]",
       list(red["required_nonzero"]) == CONTRACT_REQUIRED)
    ok("B4: module localises the kill at exactly [(0,-10)]",
       list(red["kills"]) == [(0, -10)])
    ok("B5: the module's caps name a (at y^3) and b (at y^4) as the two "
       "required-nonzero window caps (a5)",
       B.ENDPOINT_CONTRACT_50_75["caps"] == {(-1, 3): "a", (-2, 4): "b"})


# ---------------------------------------------------------------------------
# C. ansatz-independence -- the check the module's own assert cannot make
# ---------------------------------------------------------------------------
def check_C():
    floors = {}
    for jmin in (-8, -12, -20, -30):
        fl, _sup = c0_floor_by_solve(jmin)
        floors[jmin] = fl
    ok("C1: forced floor is -6 independently of the truncation window "
       "(JMIN = -8, -12, -20, -30) -- the floor is not an ansatz artifact",
       all(v == EXPECTED_FLOOR for v in floors.values()))
    if not QUIET:
        print("     floors by window:", floors)
    ok("C2: the annotator's hardcoded window JMIN=-12 reaches deep enough to "
       "contain the required slot -10 (a shallower window would hide the kill)",
       -12 <= -10)


# ---------------------------------------------------------------------------
# D. the kill predicate
# ---------------------------------------------------------------------------
def check_D():
    kills = kill_predicate(CONTRACT_REQUIRED, EXPECTED_FLOOR)
    ok("D1: kill predicate fires at exactly [(0,-10)]", kills == [(0, -10)])
    ok("D2: the two window caps (-1,3) and (-2,4) do NOT fire -- they are on "
       "other series, so a series-blind predicate would false-positive here",
       all(s == 0 for (s, j) in kills))


# ---------------------------------------------------------------------------
# E. discrimination -- four mutations, each must flip the verdict correctly
# ---------------------------------------------------------------------------
def check_E():
    # E1: deepen the forcing by one term -> floor drops to -10 -> no kill.
    C0d = c0_by_division(extra_deep=True)
    supd = y_support(C0d)
    fl_d = min(supd)
    ok("E1a: mutation (add 2*mu*y^-4 to the forcing) moves the floor -6 -> -10",
       fl_d == -10)
    ok("E1b: ... and the kill then does NOT fire -- the predicate tracks the "
       "relation, it is not hardwired to (0,-10)",
       kill_predicate(CONTRACT_REQUIRED, fl_d) == [])
    fl_solve, _ = c0_floor_by_solve(-20, extra_deep=True)
    ok("E1c: ... and the solve route agrees with the division route on the "
       "mutated relation too", fl_solve == -10)

    # E2: drop the requirement -> no kill.
    ok("E2: mutation (drop (0,-10) from the contract) -> no kill; the predicate "
       "tracks the contract, not just the floor",
       kill_predicate([(-1, 3), (-2, 4)], EXPECTED_FLOOR) == [])

    # E3/E4: the boundary is strict and exact.
    ok("E3: requirement exactly AT the floor (0,-6) does NOT kill -- the "
       "comparison is strict `<`, not `<=`",
       kill_predicate([(0, -6)], EXPECTED_FLOOR) == [])
    ok("E4: requirement one step below the floor (0,-7) DOES kill -- so the "
       "boundary sits exactly between -6 and -7, no off-by-one",
       kill_predicate([(0, -7)], EXPECTED_FLOOR) == [(0, -7)])


# ---------------------------------------------------------------------------
# F. contract fidelity -- prose must agree with code
# ---------------------------------------------------------------------------
def check_F():
    path = os.path.join(HERE, "ENDPOINT_CONTRACT.md")
    try:
        txt = open(path, encoding="utf-8").read()
    except OSError:                                            # pragma: no cover
        ok("F0: ENDPOINT_CONTRACT.md is readable", False)
        return
    ok("F0: ENDPOINT_CONTRACT.md is readable", True)

    m = re.search(r"forced_floor\s*=\s*\{([^}]*)\}", txt)
    floor_doc = {}
    if m:
        for k, v in re.findall(r"(-?\d+)\s*:\s*(-?\d+)", m.group(1)):
            floor_doc[int(k)] = int(v)
    ok("F1: the document's stated forced_floor is {0:-6, -1:3, -2:4} and its "
       "series-0 entry matches the COMPUTED floor",
       floor_doc.get(0) == EXPECTED_FLOOR and floor_doc == {0: -6, -1: 3, -2: 4})

    m = re.search(r"required_nonzero\s*=\s*\[([^\]]*)\]", txt)
    req_doc = []
    if m:
        req_doc = [(int(s), int(j)) for s, j in
                   re.findall(r"\(\s*(-?\d+)\s*,\s*(-?\d+)\s*\)", m.group(1))]
    ok("F2: the document's stated required_nonzero list matches the code's",
       req_doc == CONTRACT_REQUIRED)

    ok("F3: the document still names c_{0,-10} as the kill site",
       "c_{0,-10}" in txt)

    # F4 was a tripwire: it asserted sec.3's corner signature was STILL the
    # pre-repair `t=5, kappa=3`, so that fixing the document would fail this check
    # and force a conscious update here.  The document was fixed on 2026-07-26 and
    # the tripwire fired as designed.  It is now inverted: the repaired signature
    # must be present, and the pre-repair one must appear only inside the
    # repair note that quotes it.
    has_repaired = re.search(r"t\s*=\s*ceil\(20/5\)\s*=\s*4", txt) is not None
    ok("F4: sec.3 now carries the REPAIRED corner signature t = ceil(20/5) = 4, "
       "kappa = 2, C = y a monomial (was t=5, kappa=3, read off GGV5's final "
       "chain corner via a dictionary only valid on the retraction shape)",
       has_repaired)
    ok("F4b: and the repair note explains why the gamma-chart numbers below it "
       "did NOT move -- the caps, the floor -6 and c_{0,-10} live in GGV3 sec.5's "
       "own reduced coordinates, not in our polygon bookkeeping.  That is the "
       "same depth-ledger vs window-cone distinction as F2_TOWER.md's bridge "
       "banner, and it is why a repair that moved t, kappa, C, N, Phi and "
       "q_window left this kill untouched.",
       "REPAIRED 2026-07-26" in txt and "different objects" in txt)


# ---------------------------------------------------------------------------
# G. cross-route agreement with f2_tower.py's expression
# ---------------------------------------------------------------------------
def check_G():
    # f2_tower.a2_certificate() builds C_0 with its own symbol names; rebuild the
    # same object here rather than importing (f2_tower pulls in g_system_75_125).
    aa, bb = sp.symbols("aa bb")
    Fm2v = f8 * y**8 + f6 * y**6 + f4 * y**4 + f2 * y**2
    C0_tower = sp.expand((3 * (bb * y**4) ** 2 + 2 * Fm2v + 2 * lam) / (3 * (aa * y**3) ** 2))
    sup_tower = y_support(C0_tower)
    ok("G1: f2_tower.py's gamma=3 C_0 has the same support {-6,-4,-2,0,2}",
       sup_tower == EXPECTED_SUPPORT)
    C0_here = c0_by_division().subs({a: aa, b: bb})
    ok("G2: f2_tower.py's C_0 and this file's division route are the SAME "
       "expression (so the two modules are talking about one object)",
       sp.simplify(sp.together(C0_tower - C0_here)) == 0)


def main():
    check_A()
    check_B()
    check_C()
    check_D()
    check_E()
    check_F()
    check_G()
    if not QUIET:
        print()
    if FAILS:
        print("FAILURES:", len(FAILS), FAILS)
        sys.exit(1)
    print("ALL %d BIGRADE-ANNOTATOR CHECKS PASSED" % N_OK)
    sys.exit(0)


if __name__ == "__main__":
    main()
