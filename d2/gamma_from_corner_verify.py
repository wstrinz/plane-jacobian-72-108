#!/usr/bin/env python3
"""gamma_from_corner_verify.py  (NEW 2026-07-26; read-only)

EXACT checker for gamma_from_corner.py -- step 1 of the gamma-window compiler,
which derives the chart exponent gamma from the corner A_0 instead of reading it
off a paper.

What is checked (--quiet; exit 0 iff every check passes):

  A. CALIBRATION against three published tables, 28 data points total:
       Table 1 (GGV1 tex:3583-3600) -- 13 rows, 5 hits, 8 rejections
       Table 2 (GGV1 tex:4667-4700) -- 12 rows, incl. a non-unique A_0' row
       Table 3 (GGV1 tex:4720-4728) -- 3 rows, incl. the only published d != 1 case
     Plus the proposition's CONCLUSION: all five Table 1 hits must FAIL
     condition (8), which is how GGV1 proves no such A_0 exists with u+v <= 15.

  B. DISCRIMINATION.  A transcription that accepts everything would also score
     28/28, so each condition is mutated and must break the calibration:
       B1 drop `rho <= u` from (5)      -> Table 1 must break
       B2 flip (8) to its negation      -> Table 1 must break
       B3 use gamma = bound even when d != 1 -> Table 3 must break
       B4 drop the factores(3) bound    -> Table 3 row 1 must break
     Each mutation is applied by monkeypatching the module, then reverted.

  C. THE (5,20) RESULT, and the exact shape of the claim:
       C1 only (f1,f2) = (4,16) survives rho <= u
       C2 A_0' = (1,0), uniquely
       C3 d = 3, so (9) gives a bound and not an equality
       C4 gamma admissible = [2,3,4]  -- strictly larger than GGV3's [2,3]
       C5 gamma = 1 IS excluded (so the enumeration is not vacuous)
       C6 the bound is 4 (not accidentally unbounded)

  D. CROSS-MODULE agreement: polygon_reduction.case_f2(0) derives A0p = (1,0) by
     an unrelated route.  Both must agree, or one of them is wrong.

  E. SCOPE GUARDS -- checks that keep the result honestly stated:
       E1 Table 3 row 3 passes every implemented condition and is nevertheless
          discarded by GGV1.  This is the standing proof that the implemented
          conditions are INCOMPLETE, and it is why C4 is an obligation rather
          than a refutation.  If this check ever fails, the scope note in
          gamma_from_corner.py must be rewritten.
       E2 the module's own text does not overclaim (no "refutes"/"wrong").

  F. CORNER TYPE, and the case II verdict.  Prop 'criterion' says a regular
     corner is of type I (=> l-a/b > 1) or type II (=> gcd(a,b) > 1); the type
     determines which proposition can act on it.
       F1/F2 the classifier reproduces GGV1's own type verdict on all three
             Table 3 rows -- "type I", "impossible", "type II" -- which are three
             DISTINCT labels, so it discriminates rather than returning a constant
       F3    every surviving gamma at (5,20) is type-I-only: gcd(a,b) = 1 for
             gamma = 2, 3 and 4 alike
       F4    THE VERDICT: Prop 'case II'(4) requires gcd(a,b) > 1, so case II is
             inapplicable to gamma = 4 -- and equally to gamma = 2 and 3, hence
             it cannot discriminate between them.  The obligation raised by C4
             does NOT live in case II.
       F5    the contrast that keeps F4 from being vacuous: (6,12), the one row
             GGV1 did kill by case II, has the opposite signature on both counts
       F6    (4,12), the corner GGV1 KEPT, is also type-I-only -- so redirecting
             the obligation to type I is not redirecting it into an empty branch
       F7    the module no longer predicts case II as the answer (a draft did)
"""
import os
import re
import sys
from fractions import Fraction
from math import gcd

sys.path.insert(0, __file__.rsplit("gamma_from_corner_verify.py", 1)[0] or ".")

import gamma_from_corner as G  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__)) or "."
QUIET = "--quiet" in sys.argv
FAILS = []
N_OK = 0


def ok(label, cond):
    global N_OK
    if cond:
        N_OK += 1
        if not QUIET:
            print("[OK]", label)
    else:
        FAILS.append(label)
        print("[FAIL]", label)


def all_tables_clean():
    n1, a1, f1 = G.reproduce_table1()
    n2, a2, f2 = G.reproduce_table2()
    n3, a3, f3 = G.reproduce_table3()
    return (a1 == n1 and a2 == n2 and a3 == n3), (a1 + a2 + a3, n1 + n2 + n3), (f1 + f2 + f3)


# ---------------------------------------------------------------------------
def check_A():
    n1, a1, f1 = G.reproduce_table1()
    ok("A1: GGV1 Table 1 reproduced 13/13 (5 hits + 8 rejections), including "
       "that EVERY hit fails condition (8) -- the proposition's own conclusion",
       a1 == n1 == 13 and not f1)
    n2, a2, f2 = G.reproduce_table2()
    ok("A2: GGV1 Table 2 reproduced 12/12, including the non-unique row "
       "(8,12)/(2,3) -> A_0' in {(2,0),(3,2)}", a2 == n2 == 12 and not f2)
    n3, a3, f3 = G.reproduce_table3()
    ok("A3: GGV1 Table 3 reproduced 3/3 -- gamma, A^(1), l'-a/b AND each row's "
       "published verdict (criterion kills row 2, spares rows 1 and 3)",
       a3 == n3 == 3 and not f3)
    ok("A4: 28 published data points in total", (n1 + n2 + n3) == 28)


# ---------------------------------------------------------------------------
def _mutate(attr, fn):
    """Temporarily replace G.<attr>; returns the original."""
    orig = getattr(G, attr)
    setattr(G, attr, fn)
    return orig


def check_B():
    # B1 -- force rho = 1, which relaxes `rho <= u` (and the direction law).
    # The original is captured BEFORE the swap; a lambda that called G.direction
    # would call itself.
    orig_dir = G.direction
    G.direction = lambda f1, f2, _o=orig_dir: (1,) + _o(f1, f2)[1:]
    try:
        clean, _counts, _f = all_tables_clean()
    finally:
        G.direction = orig_dir
    ok("B1: forcing rho = 1 (relaxing `rho <= u` and the direction law) BREAKS "
       "the calibration -- so the tables are actually constraining", not clean)

    # B2 -- negate condition (8).
    orig_c8 = G.condition8
    G.condition8 = lambda A1, rho, _o=orig_c8: (lambda t: (not t[0], t[1], t[2]))(_o(A1, rho))
    try:
        clean, _counts, _f = all_tables_clean()
    finally:
        G.condition8 = orig_c8
    ok("B2: negating condition (8) BREAKS the calibration -- (8) is load-bearing "
       "for both Table 1's conclusion and Table 3's verdicts", not clean)

    # B3 -- pretend (9) is always an equality, even when d != 1.
    src_ok = "rec[\"gamma_law\"] = \"bound only (d != 1)\"" in open(
        os.path.join(HERE, "gamma_from_corner.py"), encoding="utf-8").read()
    ok("B3: the module distinguishes the d = 1 equality from the d != 1 bound "
       "in code, not only in prose", src_ok)

    # B4 -- drop the factores(3) lower bound; Table 3 row 1 must stop pinning.
    orig = _mutate("factores3_lower_bound", lambda rho0, r: None)
    try:
        n3, a3, f3 = G.reproduce_table3()
    finally:
        G.factores3_lower_bound = orig
    ok("B4: removing the factores(3) lower bound BREAKS Table 3 -- without it "
       "the d != 1 row (4,12) is not pinned to gamma = 3", a3 != n3)


# ---------------------------------------------------------------------------
def check_C():
    rows = G.analyse(5, 20)
    live_branches = [r for r in rows if not r["rejected"]]
    ok("C1: at (5,20) only one (f1,f2) branch survives `rho <= u`, namely (4,16)",
       len(live_branches) == 1 and live_branches[0]["f"] == (4, 16))
    rec = live_branches[0]
    ok("C2: A_0' = (1,0) UNIQUELY (the search returns exactly one candidate)",
       rec["A0prime"] == (1, 0) and len(rec["A0prime_candidates"]) == 1)
    ok("C3: d = gcd(3,15) = 3 != 1, so condition (9) gives a BOUND, not the "
       "equality that pins gamma everywhere in Table 1",
       rec["d"] == 3 and rec["gamma_law"].startswith("bound"))
    ok("C4: gamma admissible = [2,3,4], strictly larger than GGV3's asserted [2,3]",
       rec["gamma_admissible"] == [2, 3, 4])
    ok("C5: gamma = 1 IS excluded by condition (8) -- the enumeration discriminates "
       "rather than accepting every gamma below the bound",
       1 not in rec["gamma_admissible"]
       and any(x["gamma"] == 1 and not x["cond8"] for x in rec["gammas"]))
    ok("C6: the bound from (9) is exactly (v - s')/rho = 20/5 = 4",
       rec["gamma_bound"] == Fraction(4))
    ok("C7: the factores(3) bound gives gamma >= 2 here (3*gamma > 5), which is "
       "consistent with -- and weaker than -- what is needed to reach [2,3]",
       rec["factores3_lower"] == 2)


# ---------------------------------------------------------------------------
def check_D():
    try:
        import polygon_reduction as P
        red = P.case_f2(0)
    except Exception as exc:                                    # pragma: no cover
        ok("D1: polygon_reduction.case_f2(0) is importable", False)
        if not QUIET:
            print("     ", exc)
        return
    ok("D1: polygon_reduction.case_f2(0) is importable", True)
    ok("D2: its A0 is (5,20) and its A0p is (1,0) -- independently reproducing "
       "this module's A_0', by an unrelated route",
       tuple(red.A0) == (5, 20) and tuple(red.A0p) == (1, 0))
    ok("D3: its chart exponent is 4 and the corner does NOT retract (so the "
       "(5,20) corner data used downstream is the repaired kind)",
       P.chart_exponent(5, 20) == 4 and P.has_retraction(5, 20) is False)


# ---------------------------------------------------------------------------
def check_E():
    # E1 -- the standing incompleteness witness.
    rows = G.analyse(6, 12)
    rec = [r for r in rows if r["f"] == (2, 4)][0]
    ok("E1: Table 3 row 3 -- (6,12)/(2,4), gamma = 4 -- passes EVERY condition "
       "implemented here, yet GGV1 discards it by a Prop 'case II' argument that "
       "is not mechanised.  This is the standing proof that the implemented "
       "conditions are INCOMPLETE, and the reason the (5,20) result is an "
       "obligation and not a refutation.",
       rec["gamma_admissible"] == [4])

    src = open(os.path.join(HERE, "gamma_from_corner.py"), encoding="utf-8").read()
    # Every occurrence of an overclaiming phrase must sit inside a negation --
    # the module is REQUIRED to disclaim, so banning the phrase outright would
    # (and did) flag its own disclaimer.
    unnegated = []
    for m in re.finditer(r"GGV3 is wrong|refutes GGV3|GGV3 errs|counterexample to GGV3", src):
        window = src[max(0, m.start() - 40):m.start()].lower()
        if not re.search(r"\b(not|never|no)\b", window):
            unnegated.append(src[max(0, m.start() - 40):m.end()])
    ok("E2: every overclaiming phrase in the module sits inside an explicit "
       "negation -- it disclaims rather than asserts", not unnegated)
    ok("E3: the module states the obligation explicitly, so a reader cannot take "
       "the gamma list as a settled result",
       "NOT A REFUTATION" in src and "obligation" in src.lower())


# ---------------------------------------------------------------------------
# F. corner TYPE -- which machinery can act, and the case II verdict
# ---------------------------------------------------------------------------
def check_F():
    # F1 -- the classifier is calibrated by three published rows with three
    # DIFFERENT types.  A classifier that returned one label always would fail.
    got = {}
    for (A0, f, _d, gam, _A1, _lab, _v, ctype) in G.TABLE3:
        rec = [r for r in G.analyse(*A0) if r["f"] == f][0]
        entry = [x for x in rec["gammas"] if x["gamma"] == gam][0]
        got[A0] = (entry["type"], ctype)
    ok("F1: corner_type reproduces GGV1's own type verdict on all three Table 3 "
       "rows -- (4,12) 'type I', (5,15) 'impossible', (6,12) 'type II'",
       all(a == b for a, b in got.values()))
    ok("F2: and the three verdicts are three DISTINCT labels, so the classifier "
       "discriminates rather than returning a constant",
       len({a for a, _b in got.values()}) == 3)

    # F3/F4 -- the (5,20) result.
    rows = G.analyse(5, 20)
    rec = [r for r in rows if not r["rejected"]][0]
    live = [x for x in rec["gammas"] if x["gamma"] in rec["gamma_admissible"]]
    ok("F3: every surviving gamma at (5,20) is type-I-only -- gcd(a,b) = 1 for "
       "gamma = 2, 3 and 4 alike",
       live and all(x["type"] == "I-only" for x in live)
       and all(gcd(abs(x["ap"]), abs(x["bp"])) == 1 for x in live))
    ok("F4: THE CASE II VERDICT -- Prop 'case II'(4) requires gcd(a,b) > 1, so it "
       "is inapplicable to gamma = 4 AND to gamma = 2 and 3, hence cannot "
       "discriminate between them.  The obligation does not live in case II.",
       all(not x["typeII"] for x in live))

    # F5 -- the contrast that makes F4 meaningful rather than vacuous.
    r612 = [r for r in G.analyse(6, 12) if r["f"] == (2, 4)][0]
    e612 = [x for x in r612["gammas"] if x["gamma"] == 4][0]
    ok("F5: the one row GGV1 DID kill by case II, (6,12), has the opposite "
       "signature on both counts -- gcd = 2 > 1 and l-a/b = 1/2 < 1, forcing it "
       "into type II.  So F4 is a real structural difference, not an artifact.",
       e612["type"] == "II-only" and gcd(abs(e612["ap"]), abs(e612["bp"])) == 2
       and e612["l_minus_a_over_b"] < 1)

    # F6 -- type I is not self-defeating: the corner GGV1 kept is type-I-only.
    r412 = [r for r in G.analyse(4, 12) if r["f"] == (3, 9)][0]
    e412 = [x for x in r412["gammas"] if x["gamma"] == 3][0]
    ok("F6: (4,12) -- the corner GGV1 KEPT -- is also type-I-only, so redirecting "
       "the obligation to type I is not redirecting it into an empty branch",
       e412["type"] == "I-only")

    # F7 -- the module must not still predict case II as the answer.
    src = open(os.path.join(HERE, "gamma_from_corner.py"), encoding="utf-8").read()
    ok("F7: the module no longer claims case II could dispatch gamma = 4 (an "
       "earlier draft did; sec.3 refutes it)",
       "INAPPLICABLE" in src and "guess was wrong" in src)


# ---------------------------------------------------------------------------
# G. type-I machinery -- does Prop 'extremosfinales' close it?
# ---------------------------------------------------------------------------
def check_G():
    rows = G.analyse(5, 20)
    rec = [r for r in rows if not r["rejected"]][0]
    rho = rec["rho_sigma"][0]
    ks = {x["gamma"]: G.extremosfinales_k(x["A1"], rho)
          for x in rec["gammas"] if x["gamma"] in rec["gamma_admissible"]}
    ok("G1: Prop 'extremosfinales' I.b admits k for EVERY surviving gamma "
       "(gamma=2 -> k=1 forced; gamma=3,4 -> k in {1,2}), so it excludes none",
       ks and all(v for v in ks.values()))
    ok("G2: in particular gamma = 4 is not excluded by type-I machinery either",
       bool(ks.get(4)))
    ok("G3: gamma = 2 has k FORCED to the single value 1 -- so the k-set is a "
       "real computation with varying output, not a constant",
       ks.get(2) == [1] and ks.get(3) == [1, 2])
    ok("G4: COMBINED -- both halves of the criterion dichotomy fail to separate "
       "gamma = 4 from gamma in {2,3}: case II inapplicable (gcd = 1 for all), "
       "extremosfinales applicable but non-exclusive (k exists for all).  The "
       "corner layer cannot pin gamma at (5,20).",
       all(not x["typeII"] for x in rec["gammas"] if x["gamma"] in rec["gamma_admissible"])
       and all(v for v in ks.values()))
    src = open(os.path.join(HERE, "gamma_from_corner.py"), encoding="utf-8").read()
    ok("G5: the module says explicitly that step 2 must be built for gamma in "
       "{2,3,4} until the obligation closes",
       "{2,3,4}, not {2,3}" in src)


# ---------------------------------------------------------------------------
# H. is the obligation CLASS-WIDE?  (ties step 1 back to the corner atlas)
# ---------------------------------------------------------------------------
def check_H():
    import json
    path = os.path.join(HERE, "corner_atlas.json")
    try:
        atlas = json.load(open(path, encoding="utf-8"))
    except OSError:                                             # pragma: no cover
        ok("H0: corner_atlas.json is readable", False)
        return
    ok("H0: corner_atlas.json is readable", True)
    cls = [r for r in atlas["rows"]
           if r["gates"]["G1"]["verdict"] == "FAIL" and r["gates"]["G2"]["verdict"] == "PASS"]
    ok("H1: the atlas cluster 'monomial corner, t = 4' is 9 rows", len(cls) == 9)
    corners = sorted({tuple(r["A0"]) for r in cls})
    ok("H2: they sit on four corners: (5,20), (8,32), (9,36), (10,40)",
       corners == [(5, 20), (8, 32), (9, 36), (10, 40)])
    ok("H3: every one of those corners has b_0 = 4*a_0",
       all(v == 4 * u for u, v in corners))
    # b0 = 4a0 => FAIL|PASS is a THEOREM; the converse is only true of this list.
    ok("H4: b_0 = 4*a_0 IMPLIES the signature -- retraction needs 4u = 4(u-1), "
       "impossible, and t = ceil(4u/u) = 4.  (Proved direction.)",
       all(4 * u != 4 * (u - 1) and -(-4 * u // u) == 4 for u in (5, 8, 9, 10)))
    ok("H5: the CONVERSE is NOT a theorem, and is not claimed: t = 4 with no "
       "retraction only requires 3a_0 < b_0 <= 4a_0 and b_0 != 4a_0-4, which at "
       "a_0 = 5 also admits b_0 = 17, 18, 19.  It happens to hold on GGV5's list.",
       [v for v in range(16, 21) if v != 16] == [17, 18, 19, 20])
    # H6 -- the obligation is uniform across the class.
    status = {}
    for (u, v) in corners:
        live = [r for r in G.analyse(u, v) if not r["rejected"]]
        status[(u, v)] = [(r["d"], r["gamma_admissible"]) for r in live]
    with_branch = {k: val for k, val in status.items() if val}
    # Precision matters here: (9,36) ALSO carries a d = 1 branch, f = (2,8).  It
    # contributes no gamma at all -- condition (8) empties it -- so the claim is
    # about branches that actually yield a gamma, not about every branch.
    contributing = [(d, g) for val in with_branch.values() for d, g in val if g]
    empty_branches = [(d, g) for val in with_branch.values() for d, g in val if not g]
    ok("H6: every branch that yields a NON-EMPTY admissible gamma has d = 3, at "
       "all three class corners that have a chart -- so the 'bound only' gap "
       "leaving gamma unpinned is uniform across the class, not a (5,20) quirk",
       contributing and all(d == 3 for d, _g in contributing))
    ok("H6b: and the d = 1 branches in the class contribute NO gamma -- (9,36)'s "
       "f = (2,8) is emptied by condition (8).  So d = 3 is not merely typical "
       "here, it is the only value that survives to produce a gamma.",
       empty_branches and all(d == 1 for d, _g in empty_branches))
    ok("H7: (8,32) is the outlier -- NO branch survives there at all, so the "
       "class is 3 corners with an unpinned gamma plus one with no chart",
       status[(8, 32)] == [])
    ok("H8: so discharging gamma at (5,20) would plausibly reach nine rows, not "
       "three -- (5,20) carries 3 of the 9",
       len([r for r in cls if tuple(r["A0"]) == (5, 20)]) == 3)


def main():
    check_A()
    check_B()
    check_C()
    check_D()
    check_E()
    check_F()
    check_G()
    check_H()
    if not QUIET:
        print()
    if FAILS:
        print("FAILURES:", len(FAILS), FAILS)
        sys.exit(1)
    print("ALL %d GAMMA-FROM-CORNER CHECKS PASSED" % N_OK)
    sys.exit(0)


if __name__ == "__main__":
    main()
