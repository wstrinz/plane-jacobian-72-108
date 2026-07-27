#!/usr/bin/env python3
"""gamma_from_corner.py  (NEW 2026-07-26)

STEP 1 OF THE GAMMA-WINDOW COMPILER: derive the chart exponent `gamma` from the
corner `A_0` alone, instead of reading it off a paper.

Why this is step 1
------------------
GGV3 sec.5 analyses the (50,75) case in two charts, gamma=3 and gamma=2, and
every downstream object -- the window caps (a5), the primitivity requirement (a6)
c_{0,-10} != 0, the forced C_0, the kill -- hangs off that choice.  GGV3 states
the choice like this (tex:1722-1723):

    "using similar computations as in [GGV1, Proposition 8.3], one can check that
     necessarily gamma=3 or gamma=2"

and, for the whole construction, says outright (tex:1716-1718):

    "We do not provide proofs for this first part, since it serves only to verify
     a known case and to show the usefulness of systems like (4.4)."

So `gamma in {2,3}` is ASSERTED in the source, not proved there.  f2_tower.py
inherits the assertion: its a2_certificate() writes GGV3's equations down as
literals and references our own corner data (T, KAPPA, QC, C, ordPhi, Nof,
build_gsystem) ZERO times.  Deriving gamma is therefore the cheapest possible
falsification of the whole chain -- if the corner does not pin gamma down, every
later stage inherits the ambiguity.

What is derived here
--------------------
GGV1's Proposition "final" (tex:3400-3446) gives conditions on the tuple
(A_0, (f_1,f_2), (rho,sigma), A_0', gamma, A^(1)).  Transcribed exactly:

  (5)  u*f_2 = v*f_1  and  rho <= u
  (6)  (rho,sigma) = dir(f_1-1, f_2-1) = ((f_2-1)/d, (1-f_1)/d),  d = gcd(f_1-1,f_2-1)
  (7)  A^(1) = A_0' + (gamma - s')*(-sigma/rho, 1)
  (8)  if A^(1) = (a'/rho, b')  then  rho - a'/b' > 1  or  gcd(a',b') > 1
  (9)  gamma <= (v - s')/rho.   Moreover if d = 1 then gamma = (v - s')/rho.

with A_0 = (u,v), A_0' = (r',s') required to satisfy 0 <= s' < r' < u and
v_{rho,sigma}(u,v) = v_{rho,sigma}(r',s'), where v_{rho,sigma}(a,b) = rho*a + sigma*b.

THE CALIBRATION (this is what makes the module trustworthy).  GGV1 Table 1
(tex:3583-3600) publishes the full enumeration for every A_0=(u,v) with
v > u > 2, gcd(u,v) > 1 and u+v <= 15: thirteen (A_0, (f_1,f_2)) rows, of which
exactly five admit an A_0', with the resulting (rho,sigma), A_0', d, gamma and
A^(1) printed.  This module reproduces all thirteen rows -- both the five hits
AND the eight rejections -- from the conditions alone.  A derivation that
reproduces a published table it was not fitted to is worth more than one that
merely runs.

THE RESULT AT (5,20), stated honestly
-------------------------------------
The corner does NOT pin gamma to {2,3}.  It gives gamma in {2,3,4}:

  * Of the three candidate (f_1,f_2) = mu*(5,20), only (4,16) survives rho <= u
    (the others give rho = 7 and rho = 11, both > 5).
  * (rho,sigma) = (5,-1) and A_0' = (1,0) UNIQUELY -- which independently
    reproduces polygon_reduction.case_f2(0)'s A0p = (1,0), derived there by a
    different route.
  * d = gcd(3,15) = 3 != 1, so condition (9) gives only the BOUND gamma <= 4,
    not the equality.  This is exactly the case the table never exercises: all
    five of its hits have d = 1.
  * Condition (8) kills gamma = 1 (A^(1) = (6/5,1): rho - a'/b' = -1, and
    gcd(6,1) = 1).  It does NOT kill gamma = 4 (A^(1) = (9/5,4):
    rho - a'/b' = 11/4 > 1).

So gamma = 4 satisfies every condition of GGV1's Proposition "final".  Excluding
it needs the "similar computations as in GGV1 Prop 8.3" that GGV3 cites but does
not carry out.  Until that argument is reconstructed, "gamma in {2,3}" is
UNDERIVED here, and a gamma=4 chart -- a third chart GGV3 does not analyse --
is not ruled out by the corner data.

SCOPE -- read this before quoting the result
--------------------------------------------
**This module does NOT claim GGV3 is wrong, and the gamma=4 branch is NOT
established as real.**  What is established is narrower: the conditions
transcribed here -- which reproduce 28/28 published data points across three
tables, including the d != 1 branch and the proposition's own conclusion -- do
not by themselves exclude gamma = 4 at (5,20).

GGV1 carries substantially more machinery that is NOT mechanised here:
Propositions 'criterion' (beyond the simple form), 'case II', 'case IIb',
'impossibles', 'u(u-1)', 'encima de la diagonal', and 'factores' in full.  That
this matters is not speculation -- it is visible in the calibration: Table 3's
THIRD row, (6,12), passes every condition implemented here and is nevertheless
discarded by GGV1, via a page of Prop 'case II' analysis (mu in {1,3,5}, a
derived (rho_1,sigma_1) = dir(...), then a further regular-corner argument at
st_{3,-2}(P)).

So the honest reading is: **gamma in {2,3} is UNDERIVED at this level of the
machinery, not refuted.**  The useful output is a located obligation --
"exclude gamma = 4 at (5,20), or exhibit the third chart" -- with the precise
reason it does not follow from what is transcribed: d != 1 kills the equality in
condition (9).

WHERE THE MISSING ARGUMENT CANNOT LIVE (sec.3 of the output)
------------------------------------------------------------
A first draft of this note guessed that Prop 'case II' was where the missing
argument would be found.  **That guess was wrong, and the module now refutes it.**

Prop 'criterion' (tex:2974) says a regular corner is of type I or type II, with
l - a/b > 1 in the first case and gcd(a,b) > 1 in the second.  At (5,20) the
surviving gammas give A^(1) = (7/5,2), (8/5,3), (9/5,4), so
gcd(a,b) = gcd(7,2) = gcd(8,3) = gcd(9,4) = 1 in every case.  Prop 'case II'(4)
requires gcd(a,b) > 1.  **Prop 'case II' is therefore INAPPLICABLE to gamma = 4
-- and equally inapplicable to gamma = 2 and 3, so it cannot discriminate between
them at all.**

The contrast with the one row GGV1 did kill by case II is exact.  At (6,12),
gcd(10,4) = 2 > 1 AND l - a/b = 1/2 < 1: type I is excluded and the corner is
FORCED into type II, which is what made the case II analysis available.  (5,20)
has the opposite signature on both counts.

So the obligation REDIRECTS to type-I machinery (Prop 'extremosfinales', the
source of the l - a/b > 1 condition).  This is not a dead end: (4,12), the corner
GGV1 KEPT, is also type-I-only -- type I is the branch on which a real case
survived, not a branch that self-destructs.

Checker: gamma_from_corner_verify.py (--quiet, exit 0).
"""
from __future__ import annotations

import sys
from fractions import Fraction
from math import gcd

# ---------------------------------------------------------------------------
# GGV1 Proposition "final", conditions (5)-(9).  Transcribed, not paraphrased.
# ---------------------------------------------------------------------------


def valuation(rho, sigma, a, b):
    """v_{rho,sigma}(a,b) = rho*a + sigma*b."""
    return rho * a + sigma * b


def direction(f1, f2):
    """(6): (rho,sigma) = dir(f_1-1, f_2-1) = ((f_2-1)/d, (1-f_1)/d)."""
    d = gcd(f1 - 1, f2 - 1)
    if d == 0:
        return None, None, 0
    return (f2 - 1) // d, (1 - f1) // d, d


def candidate_f(u, v):
    """All (f_1,f_2) = mu*(u,v) with f_1 >= 2 and 0 < mu < 1, integral.

    mu = k/g with g = gcd(u,v); k = 1..g-1 enumerates every integral scaling
    strictly between 0 and 1.
    """
    g = gcd(u, v)
    out = []
    for k in range(1, g):
        f1, f2 = u * k // g, v * k // g
        if u * k % g or v * k % g:
            continue
        if f1 >= 2:
            out.append((f1, f2))
    return out


def find_A0prime(u, v, rho, sigma):
    """A_0' = (r',s') with 0 <= s' < r' < u and v_{rho,sigma} equal to A_0's.

    Returns every solution; the published table always has at most one.
    """
    target = valuation(rho, sigma, u, v)
    hits = []
    for rp in range(1, u):
        for sp in range(0, rp):
            if valuation(rho, sigma, rp, sp) == target:
                hits.append((rp, sp))
    return hits


def factores3_lower_bound(rho0, r):
    """GGV1 Prop 'factores'(3), as USED at tex:4711-4713:

        "Since gcd(rho_0, r) = gcd(4,3) = 1, we have gamma(rho_0 - 2) > rho_0,
         which implies gamma >= 3."

    Returns the least integer gamma with gamma*(rho0-2) > rho0, or None when the
    hypothesis gcd(rho0,r) = 1 fails or rho0 <= 2 (no bound).

    r is read as f_1: in the cited instance (f_1,f_2) = (3,9) and the text says
    r = 3.  That identification is INFERRED from a single worked case, so it is
    reported alongside the bound rather than silently applied -- but it is
    corroborated: it pins gamma = 3 at (4,12) exactly as published, and is
    consistent with (non-contradicting at) both other Table 3 rows.
    """
    if gcd(rho0, r) != 1 or rho0 <= 2:
        return None
    g = 1
    while g * (rho0 - 2) <= rho0:
        g += 1
    return g


def corner_type(A1, rho):
    """Prop 'criterion' (GGV1 tex:2974-2979), stated in full:

        "If (a/l,b) is the first entry of a regular corner of an (m,n)-pair in
         L^(l), then it is the first entry of a regular corner of a (possibly
         different) (m,n)-pair in L^(l) OF TYPE I OR TYPE II.  Moreover, in the
         first case l - a/b > 1, while in the second one gcd(a,b) > 1."

    So the criterion is a DISJUNCTION over corner types, and each disjunct comes
    with its own necessary condition.  Condition (8) of Prop 'final' is exactly
    this disjunction -- which is why the two are one test, not two.

    Reading it type-by-type is what makes it useful here: it says WHICH kind of
    corner a surviving gamma would have to be, and therefore WHICH machinery
    could possibly discharge it.

      'I-only'   -- gcd(a,b) = 1 kills type II, so Prop 'case II' is INAPPLICABLE
      'II-only'  -- l - a/b <= 1 kills type I; only here can case II do work
      'both'     -- neither necessary condition is violated
      'neither'  -- condition (8) fails; the corner is excluded outright

    Returns (label, typeI_possible, typeII_possible, a, b).
    """
    _okay, ap, bp = condition8(A1, rho)
    if ap is None or bp in (None, 0):
        return "neither", False, False, ap, bp
    t1 = Fraction(rho) - Fraction(ap, bp) > 1
    t2 = gcd(abs(ap), abs(bp)) > 1
    label = ("both" if (t1 and t2) else
             "I-only" if t1 else
             "II-only" if t2 else "neither")
    return label, t1, t2, ap, bp


def extremosfinales_k(A1, rho):
    """Prop 'extremosfinales' (Cases I.a/I.b, GGV1 tex:2447).

    Hypothesis [l_{rho,sigma}(P), l_{rho,sigma}(Q)] != 0 -- i.e. TYPE I, the
    complement of case II's bracket-zero hypothesis.  Then l - a/b > 1 and:

      I.a) st(P) ~ st(Q)  =>  (1/m)st(P) in (1/l)Z x N_0 and st(P) ~ (1,0)
      I.b) st(P) !~ st(Q) =>  there is k in N with k < l - a/b and
                              {st(P),st(Q)} = {(k/l, 0), (1 - k/l, 1)}

    Returns the admissible k for branch I.b.  NOTE WHAT THIS PROPOSITION IS: it
    says what a type-I corner LOOKS LIKE, not that one cannot exist.  An empty
    k-set would be an exclusion; a non-empty one is not evidence of anything.
    """
    _okay, ap, bp = condition8(A1, rho)
    if ap is None or bp in (None, 0):
        return None
    lab = Fraction(rho) - Fraction(ap, bp)
    return [k for k in range(1, int(lab) + 2) if Fraction(k) < lab]


def A1_of(A0p, gamma, rho, sigma):
    """(7): A^(1) = A_0' + (gamma - s')*(-sigma/rho, 1)."""
    rp, sp = A0p
    first = Fraction(rp) + Fraction(gamma - sp) * Fraction(-sigma, rho)
    second = sp + (gamma - sp)
    return (first, second)


def condition8(A1, rho):
    """(8): writing A^(1) = (a'/rho, b'), require rho - a'/b' > 1 or gcd(a',b')>1.

    Returns (satisfied, a', b').  a' is recovered as rho * (first coordinate),
    which must be an integer for the condition to be meaningful.
    """
    first, second = A1
    ap_frac = Fraction(rho) * first
    if ap_frac.denominator != 1:
        return False, None, None
    ap, bp = int(ap_frac), int(second)
    if bp == 0:
        return False, ap, bp
    return (Fraction(rho) - Fraction(ap, bp) > 1) or (gcd(abs(ap), abs(bp)) > 1), ap, bp


def analyse(u, v, gamma_max_scan=40):
    """Full enumeration for one corner A_0 = (u,v).  Returns one record per
    (f_1,f_2) branch, each either rejected (with a reason) or carrying its
    admissible gamma values."""
    rows = []
    for (f1, f2) in candidate_f(u, v):
        rec = {"A0": (u, v), "f": (f1, f2), "rejected": None}
        # (5) u*f_2 = v*f_1 holds by construction; assert it rather than assume.
        rec["cond5_proportional"] = (u * f2 == v * f1)
        rho, sigma, d = direction(f1, f2)
        rec["rho_sigma"], rec["d"] = (rho, sigma), d
        if not rec["cond5_proportional"]:
            rec["rejected"] = "(5) u*f_2 != v*f_1"
            rows.append(rec)
            continue
        if rho is None or rho > u:
            rec["rejected"] = "(5) rho > u"
            rows.append(rec)
            continue
        hits = find_A0prime(u, v, rho, sigma)
        rec["A0prime_candidates"] = hits
        if not hits:
            rec["rejected"] = "no A_0' with 0 <= s' < r' < u and equal valuation"
            rows.append(rec)
            continue
        rec["A0prime"] = hits[0]
        rp, sp = hits[0]
        bound = Fraction(v - sp, rho)
        rec["gamma_bound"] = bound
        lo = factores3_lower_bound(rho, f1)
        rec["factores3_lower"] = lo
        if d == 1:
            # (9) equality
            gammas = [int(bound)] if bound.denominator == 1 else []
            rec["gamma_law"] = "equality (d = 1)"
        else:
            gammas = [g for g in range(1, min(int(bound), gamma_max_scan) + 1)]
            rec["gamma_law"] = "bound only (d != 1)"
        admissible = []
        for g in gammas:
            A1 = A1_of((rp, sp), g, rho, sigma)
            okay, ap, bp = condition8(A1, rho)
            passes_lo = (lo is None) or (g >= lo)
            ctype, t1, t2, _a, _b = corner_type(A1, rho)
            admissible.append({"gamma": g, "A1": A1, "cond8": okay, "ap": ap, "bp": bp,
                               "factores3": passes_lo, "type": ctype,
                               "typeI": t1, "typeII": t2,
                               "l_minus_a_over_b": (Fraction(rho) - Fraction(ap, bp))
                               if (ap is not None and bp not in (None, 0)) else None})
        rec["gammas"] = admissible
        # Condition (8) IS Prop 'criterion' in the form GGV1 applies at tex:4740
        # ("impossible, since gcd(13,3)=1 and 1-1/3<1" is exactly (8) failing:
        # neither rho - a'/b' > 1 nor gcd(a',b') > 1).  They are one test, not
        # two, so it is applied once -- separating the gamma the tables TABULATE
        # from the gamma that SURVIVES.
        rec["gamma_pre_criterion"] = [x["gamma"] for x in admissible if x["factores3"]]
        rec["gamma_admissible"] = [x["gamma"] for x in admissible
                                   if x["cond8"] and x["factores3"]]
        rows.append(rec)
    return rows


# ---------------------------------------------------------------------------
# GGV1 Table 1 (tex:3583-3600) -- the published oracle, transcribed verbatim.
#   (u,v), (f1,f2), (rho,sigma), A_0' or None, d or None, gamma or None, A^(1)
# ---------------------------------------------------------------------------
TABLE1 = [
    ((3, 6),  (2, 4),  (3, -1), (1, 0), 1, 2, (Fraction(5, 3), 2)),
    ((3, 9),  (2, 6),  (5, -1), None, None, None, None),
    ((3, 12), (2, 8),  (7, -1), None, None, None, None),
    ((4, 6),  (2, 3),  (2, -1), (1, 0), 1, 3, (Fraction(5, 2), 3)),
    ((4, 8),  (2, 4),  (3, -1), None, None, None, None),
    ((4, 8),  (3, 6),  (5, -2), None, None, None, None),
    ((4, 10), (2, 5),  (4, -1), None, None, None, None),
    ((5, 10), (2, 4),  (3, -1), (2, 1), 1, 3, (Fraction(8, 3), 3)),
    ((5, 10), (3, 6),  (5, -2), (1, 0), 1, 2, (Fraction(9, 5), 2)),
    ((5, 10), (4, 8),  (7, -3), None, None, None, None),
    ((6, 8),  (3, 4),  (3, -2), None, None, None, None),
    ((6, 9),  (2, 3),  (2, -1), (2, 1), 1, 4, (Fraction(7, 2), 4)),
    ((6, 9),  (4, 6),  (5, -3), None, None, None, None),
]


# ---------------------------------------------------------------------------
# GGV1 Table 2 (tex:4667-4700): every (u,v) with 3 < u < v <= u(u-1),
# gcd(u,v) > 2 and 16 <= u+v <= 20.  Columns: A_0, (f1,f2), (rho,sigma), A_0'.
# "(2,0) or (3,2)" is recorded as BOTH solutions -- a row where the A_0' search
# is genuinely non-unique, which nothing in Table 1 exercises.
# ---------------------------------------------------------------------------
TABLE2 = [
    ((4, 12), (2, 6),  (5, -1),  []),
    ((4, 12), (3, 9),  (4, -1),  [(1, 0)]),
    ((5, 15), (2, 6),  (5, -1),  [(2, 0)]),
    ((5, 15), (3, 9),  (4, -1),  []),
    ((5, 15), (4, 12), (11, -3), []),
    ((6, 12), (2, 4),  (3, -1),  [(2, 0)]),
    ((6, 12), (3, 6),  (5, -2),  []),
    ((6, 12), (4, 8),  (7, -3),  []),
    ((6, 12), (5, 10), (9, -4),  []),
    ((8, 12), (2, 3),  (2, -1),  [(2, 0), (3, 2)]),
    ((8, 12), (4, 6),  (5, -3),  []),
    ((8, 12), (6, 9),  (8, -5),  []),
]

# ---------------------------------------------------------------------------
# GGV1 Table 3 (tex:4720-4728): the three surviving cases, with gamma, A^(1), d
# and l' - a/b.  THE IMPORTANT ROW IS THE FIRST: (4,12) has d = 2 != 1, so its
# gamma is NOT pinned by condition (9)'s equality -- it is pinned by (9)'s bound
# TOGETHER WITH the factores(3) lower bound.  That is exactly the branch (5,20)
# sits in, and it is the only published instance of it.
#   A_0, (f1,f2), d, gamma, A^(1), l'-a/b
# ---------------------------------------------------------------------------
# The final element records GGV1's own VERDICT on the row, from the prose that
# FOLLOWS the table (tex:4733-4772):
#   'survives'      -- "we may have a regular corner of type I at A^(1)"
#   'criterion'     -- "the second case is impossible, since gcd(13,3)=1 and
#                       1-1/3<1"  (Prop 'criterion', the simple form)
#   'caseII'        -- discarded only after a long Prop 'case II' analysis
#                      (mu in {1,3,5}, (rho_1,sigma_1) = dir(...), then a
#                      further regular-corner argument).  NOT reachable by the
#                      simple criterion, and NOT mechanised here.
# The last element is the CORNER TYPE GGV1's prose assigns, quoted:
#   row 1  "we may have a regular corner of type~I at A^(1)"       -> I-only
#   row 2  "the second case is impossible"                         -> neither
#   row 3  "there could be a regular corner of type~II"            -> II-only
# Three rows, three DIFFERENT classifications, each matching GGV1's own words.
# That is what calibrates corner_type() -- and it is why the machinery available
# to discharge a row differs from row to row.
TABLE3 = [
    ((4, 12), (3, 9), 2, 3, (Fraction(7, 4), 3), Fraction(5, 3), "survives", "I-only"),
    ((5, 15), (2, 6), 1, 3, (Fraction(13, 5), 3), Fraction(2, 3), "criterion", "neither"),
    ((6, 12), (2, 4), 1, 4, (Fraction(10, 3), 4), Fraction(1, 2), "caseII", "II-only"),
]


def reproduce_table2():
    n = agree = 0
    fails = []
    for (A0, f, rs, A0p_list) in TABLE2:
        u, v = A0
        n += 1
        rows = [r for r in analyse(u, v) if r["f"] == f]
        if not rows:
            fails.append((A0, f, "branch not enumerated"))
            continue
        rec = rows[0]
        if rec["rho_sigma"] != rs:
            fails.append((A0, f, "rho,sigma %s != %s" % (rec["rho_sigma"], rs)))
            continue
        got = rec.get("A0prime_candidates", [])
        if sorted(got) != sorted(A0p_list):
            fails.append((A0, f, "A_0' set %s != published %s" % (sorted(got), sorted(A0p_list))))
            continue
        agree += 1
    return n, agree, fails


def reproduce_table3():
    """Two things are checked per row, and they are different in kind.

    (i)  The TABLE's own numbers -- d, gamma, A^(1), l'-a/b -- which are computed
         BEFORE Prop 'criterion' is applied.
    (ii) The VERDICT that GGV1's prose then reaches on that row, which is where
         the criterion does its work.  Checking (ii) is what makes the criterion
         a discriminating test rather than a filter that could accept anything:
         it must kill row 2 and spare rows 1 and 3.
    """
    n = agree = 0
    fails = []
    for (A0, f, d, gam, A1, lab, verdict, ctype) in TABLE3:
        u, v = A0
        n += 1
        rows = [r for r in analyse(u, v) if r["f"] == f]
        rec = rows[0] if rows else None
        if rec is None or rec.get("rejected"):
            fails.append((A0, f, "rejected or missing"))
            continue
        if rec["d"] != d:
            fails.append((A0, f, "d %s != %s" % (rec["d"], d)))
            continue
        if rec["gamma_pre_criterion"] != [gam]:
            fails.append((A0, f, "gamma (pre-criterion) %s != published [%s]"
                          % (rec["gamma_pre_criterion"], gam)))
            continue
        entry = [x for x in rec["gammas"] if x["gamma"] == gam][0]
        if entry["A1"] != A1:
            fails.append((A0, f, "A^(1) %s != %s" % (entry["A1"], A1)))
            continue
        if entry["l_minus_a_over_b"] != lab:
            fails.append((A0, f, "l'-a/b %s != %s" % (entry["l_minus_a_over_b"], lab)))
            continue
        # (ii) the verdict.  A row is killed exactly when condition (8) fails.
        killed = not entry["cond8"]
        if verdict == "criterion" and not killed:
            fails.append((A0, f, "GGV1 discards this row by Prop 'criterion'; we kept it"))
            continue
        if verdict in ("survives", "caseII") and killed:
            fails.append((A0, f, "we killed a row GGV1's simple criterion does NOT kill"))
            continue
        # (iii) the corner TYPE GGV1's prose assigns to the row.
        if entry["type"] != ctype:
            fails.append((A0, f, "corner type %s != GGV1's %s" % (entry["type"], ctype)))
            continue
        agree += 1
    return n, agree, fails


def table1_rows_for(u, v):
    return [r for r in TABLE1 if r[0] == (u, v)]


def reproduce_table1():
    """Recompute every Table 1 row from the conditions.  Returns
    (n_rows, n_agree, failures)."""
    n_rows = n_agree = 0
    failures = []
    for (A0, f, rs, A0p, d, gam, A1) in TABLE1:
        u, v = A0
        rows = [r for r in analyse(u, v) if r["f"] == f]
        n_rows += 1
        if not rows:
            failures.append((A0, f, "branch not enumerated at all"))
            continue
        rec = rows[0]
        if rec["rho_sigma"] != rs:
            failures.append((A0, f, "rho,sigma %s != published %s" % (rec["rho_sigma"], rs)))
            continue
        if A0p is None:
            if rec["rejected"] is None:
                failures.append((A0, f, "published X (no A_0') but we accepted it"))
                continue
        else:
            if rec["rejected"] is not None:
                failures.append((A0, f, "published a hit but we rejected: %s" % rec["rejected"]))
                continue
            if rec["A0prime"] != A0p:
                failures.append((A0, f, "A_0' %s != published %s" % (rec["A0prime"], A0p)))
                continue
            if rec["d"] != d:
                failures.append((A0, f, "d %s != published %s" % (rec["d"], d)))
                continue
            got_gam = [x["gamma"] for x in rec["gammas"]]
            if got_gam != [gam]:
                failures.append((A0, f, "gamma %s != published [%s]" % (got_gam, gam)))
                continue
            if rec["gammas"][0]["A1"] != A1:
                failures.append((A0, f, "A^(1) %s != published %s" % (rec["gammas"][0]["A1"], A1)))
                continue
            # GGV1 tex:3616-3618: "Finally we verify that in none of this cases
            # condition (8) of Proposition final is satisfied, concluding the
            # proof."  So EVERY hit in Table 1 must fail (8) -- that failure is
            # the proposition's whole conclusion, and reproducing it is a much
            # stronger check than matching the tabulated numbers.
            if rec["gammas"][0]["cond8"]:
                failures.append((A0, f, "condition (8) SATISFIED, but GGV1 states "
                                        "it fails for every Table 1 hit"))
                continue
        n_agree += 1
    return n_rows, n_agree, failures


# ---------------------------------------------------------------------------
def _fmt(rec):
    f = rec["f"]
    if rec["rejected"]:
        return "  (f1,f2)=%-8s (rho,sigma)=%-9s d=%-2s  REJECTED: %s" % (
            f, rec["rho_sigma"], rec["d"], rec["rejected"])
    return ("  (f1,f2)=%-8s (rho,sigma)=%-9s d=%-2s  A_0'=%s  gamma<=%s  [%s]\n"
            "        admissible gamma (cond 8): %s" % (
                f, rec["rho_sigma"], rec["d"], rec["A0prime"], rec["gamma_bound"],
                rec["gamma_law"], rec["gamma_admissible"]))


def main():
    print("=" * 78)
    print("1. CALIBRATION -- reproduce GGV1 Table 1 from the conditions alone")
    print("=" * 78)
    total_n = total_a = 0
    allfails = []
    for name, fn, note in (
            ("Table 1 (u+v <= 15; 13 rows, 5 hits)", reproduce_table1,
             "exact (rho,sigma), A_0', d, gamma, A^(1); 8 rejections"),
            ("Table 2 (16 <= u+v <= 20; 12 rows)", reproduce_table2,
             "A_0' sets, including the non-unique row (8,12) -> {(2,0),(3,2)}"),
            ("Table 3 (the 3 survivors)", reproduce_table3,
             "gamma, A^(1) and l'-a/b -- row 1 is the d != 1 branch")):
        n, agree, fails = fn()
        total_n += n
        total_a += agree
        allfails += fails
        print("  %-38s %2d/%2d   %s" % (name, agree, n, note))
        for f in fails:
            print("      MISMATCH:", f)
    print("  " + "-" * 70)
    print("  TOTAL published data points reproduced: %d / %d" % (total_a, total_n))
    if not allfails:
        print("  Nothing here was fitted: the conditions were transcribed from")
        print("  Prop 'final' and the tables recomputed from them.")

    print()
    print("=" * 78)
    print("2. THE TARGET -- A_0 = (5,20), the (50,75) / (75,125) corner")
    print("=" * 78)
    rows = analyse(5, 20)
    for rec in rows:
        print(_fmt(rec))
    live = sorted({g for rec in rows if not rec["rejected"] for g in rec["gamma_admissible"]})
    print()
    print("  gamma admissible from the corner alone: %s" % live)
    print("  GGV3 sec.5 asserts (tex:1723):          [2, 3]")
    if set(live) != {2, 3}:
        extra = sorted(set(live) - {2, 3})
        print("  => NOT DERIVED.  %s survives every condition of GGV1 Prop 'final'." % extra)
        print("     The gap is condition (9): d = gcd(f1-1,f2-1) = 3 != 1, so (9) gives")
        print("     only the BOUND gamma <= 4, never the equality.  All five hits in")
        print("     Table 1 have d = 1, so the published calibration never exercises")
        print("     this branch.  Excluding gamma = 4 needs the 'similar computations")
        print("     as in GGV1 Prop 8.3' that GGV3 cites but does not carry out.")
        print()
        print("     NOT A REFUTATION.  Table 3's third row, (6,12), passes every")
        print("     condition implemented here and is STILL discarded by GGV1 -- via")
        print("     a page of Prop 'case II' analysis this module does not mechanise.")
        print("     (Whether that machinery can reach gamma = 4 is settled in sec.3")
        print("     below -- it cannot.)  The output is a located OBLIGATION, not a")
        print("     counterexample:")
        print("       exclude gamma = 4 at (5,20), or a third chart exists that")
        print("       GGV3 sec.5 does not analyse.")

    # --- which machinery could possibly discharge it? -----------------------
    print()
    print("-" * 78)
    print("3. WHICH MACHINERY APPLIES -- corner type per surviving gamma")
    print("-" * 78)
    print("  Prop 'criterion': a regular corner is of type I (=> l-a/b > 1) or")
    print("  type II (=> gcd(a,b) > 1).  The type says which proposition can act.")
    print()
    print("   gamma   A^(1)      a   b   gcd(a,b)  l-a/b   type")
    for rec in rows:
        if rec["rejected"]:
            continue
        for x in rec["gammas"]:
            if x["gamma"] not in rec["gamma_admissible"]:
                continue
            print("     %d    (%s, %d)  %2d  %2d      %d      %-6s  %s" % (
                x["gamma"], x["A1"][0], x["A1"][1], x["ap"], x["bp"],
                gcd(abs(x["ap"]), abs(x["bp"])), x["l_minus_a_over_b"], x["type"]))
    types = {x["type"] for rec in rows if not rec["rejected"]
             for x in rec["gammas"] if x["gamma"] in rec["gamma_admissible"]}
    print()
    if types == {"I-only"}:
        print("  EVERY surviving gamma is TYPE-I-ONLY: gcd(a,b) = 1 in all of them,")
        print("  and Prop 'case II'(4) requires gcd(a,b) > 1.  So **Prop 'case II'")
        print("  is INAPPLICABLE at (5,20)** -- it cannot discharge gamma = 4,")
        print("  because it cannot act on gamma = 2 or 3 either.")
        print()
        print("  Contrast Table 3 row 3, (6,12): there gcd(10,4) = 2 > 1 AND")
        print("  l-a/b = 1/2 < 1, so type I was excluded and the corner was FORCED")
        print("  into type II -- which is exactly why GGV1 could run the case II")
        print("  analysis and kill it.  (5,20) has the opposite signature.")
        print()
        print("  The obligation therefore REDIRECTS: discharging gamma = 4 needs")
        print("  TYPE-I machinery (Prop 'extremosfinales', which is what supplies")
        print("  l - a/b > 1), not Prop 'case II'.  Note (4,12) -- the corner GGV1")
        print("  kept -- is also type-I-only, so type I is not self-defeating: it is")
        print("  the branch on which a real case survived.")
    print()
    print("  Cross-check: polygon_reduction.case_f2(0) records A0p = (1,0), derived")
    print("  by a different route.  This module gets A_0' = (1,0) uniquely.")

    # --- section 4: does type-I machinery close it? -------------------------
    print()
    print("-" * 78)
    print("4. TYPE-I MACHINERY -- Prop 'extremosfinales', the other half of the")
    print("   dichotomy.  Does IT discharge gamma = 4?")
    print("-" * 78)
    print("   gamma   l-a/b   admissible k (I.b: k in N, k < l-a/b)")
    live = [(x, rec) for rec in rows if not rec["rejected"]
            for x in rec["gammas"] if x["gamma"] in rec["gamma_admissible"]]
    ks = {}
    for x, rec in live:
        kk = extremosfinales_k(x["A1"], rec["rho_sigma"][0])
        ks[x["gamma"]] = kk
        print("     %d     %-6s  %s" % (x["gamma"], x["l_minus_a_over_b"], kk))
    print()
    if all(v for v in ks.values()):
        print("  NO.  Every surviving gamma admits at least one k, so Prop")
        print("  'extremosfinales' excludes none of them.  And note WHAT THE")
        print("  PROPOSITION IS: it says what a type-I corner LOOKS LIKE (its k, its")
        print("  starting pair), not that one cannot exist.  Only an EMPTY k-set would")
        print("  be an exclusion.  GGV1 uses it the same way -- of (4,12) it concludes")
        print('  "we MAY have a regular corner of type I at A^(1)", and that is exactly')
        print("  why (4,12) SURVIVED as the answer rather than being discarded.")
        print()
        print("  COMBINED VERDICT.  Both halves of Prop 'criterion's dichotomy fail to")
        print("  separate gamma = 4 from the accepted gamma in {2,3} at (5,20):")
        print("    type II  -- INAPPLICABLE (gcd(a,b) = 1 for all three)   [sec.3]")
        print("    type I   -- APPLICABLE but non-exclusive (k exists for all three)")
        print("  So the corner-level machinery of GGV1 -- Props 'final', 'criterion',")
        print("  'case II', 'extremosfinales' -- CANNOT pin gamma at (5,20).  Whatever")
        print("  GGV3's 'similar computations as in [GGV1, Prop 8.3]' does, it must act")
        print("  BELOW the corner layer: in the reduction to (P_1,Q_1) and the degree /")
        print("  polygon bookkeeping of the reduced chart.")
        print()
        print("  That is the compiler's step 2, and this narrows its brief: the reduced")
        print("  chart is not merely where the window caps (a5) live, it is the layer")
        print("  where gamma itself gets pinned.  Until it does, step 2 must be built")
        print("  for gamma in {2,3,4}, not {2,3}.")
    print("=" * 78)


if __name__ == "__main__":
    sys.exit(main())
