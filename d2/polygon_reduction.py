#!/usr/bin/env python3
"""polygon_reduction.py  (NEW; read-only over all existing artifacts)

THE POLYGON-REDUCTION COMPILER -- the missing front end of the corner program.

Every derived corner in this program (PHI_75_125, PHI_CORNER4, PHI_F14, PHI_F7,
C_SERIES_75_125, G_SYSTEM_75_125) carries the same standing judgment:

    [judgment] unreduced polygon -- the standard type-II.b root-shift + Laurent
    chart reduction is ASSUMED (t=l, kappa=l-2, deg C=a0); it is written out in
    no paper except the published (8,28) reduction of GGHV22.

2026-07-26 REPAIR.  The clause "deg C = a0" above is only true on the
RETRACTION SHAPE b0 = l*(a0-1), and the l fed to it must be the CHART exponent
ceil(b0/a0) -- NOT the denominator of GGV5's final chain corner.  Conflating the
two put (t,kappa,C,q) = (5,3,y^2(y^3+1),2) into this repo for the (5,20) corner
where the truth is (4,2,y,1).  Section 0b below derives l and GUARDS the
dictionary; see PASSPORT_75_125_REPAIR.md.

This module converts that ASSUMPTION into a DERIVATION.  Given a GGV-chain case
(A0, A0', chain/(m,n) data) it emits the COMPLETE reduction:

  1. the exact sequence of coordinate transformations -- root-shift shears
     (Jacobian 1, bracket-preserving) followed by the ONE final Laurent
     inversion.  It builds on the FUSED-CHART LEMMA already proved in
     composite_charts.py:
         (X, Y) = (x^-1,  x^l y + sum_i lambda_i x^(e_i)),   Jacobian -x^(l-2).
     (We do NOT rederive that lemma; we confirm the Jacobian for the specific l
     of each case and quote composite_charts.py for the general statement.)
  2. the transformed bracket exponent  [P,Q] = x^kappa,  kappa = l-2 DERIVED
     from that Jacobian (never assumed).
  3. the reduced Newton polygons of BOTH P and Q -- vertex lists, COMPUTED by
     pushing the pre-inversion vertices through the inversion map
         (a, b)  |-->  (l*b - a, b).
  4. EVERY alternative branch -- opposite-edge choices, endpoint alternatives,
     factor-count / denominator branches -- enumerated explicitly, each tagged
     FOLLOWED or EXCLUDED with a reason.  This branch manifest is the point:
     silent selection is exactly what the judgment flags.
  5. the corner signature (t, kappa, a0, q, c(y)) that the corner law consumes.

REGRESSION CONTRACT (checked by polygon_reduction_verify.py):
  R1  reproduce the PUBLISHED (72,108)/(8,28) reduction EXACTLY
      (matches paper_src/upstream_facts.json sub1/sub2 and [P,Q]=x^2).
  R2  reproduce F2 j=0 = (50,75) with the REPAIRED signature
      (t=4, kappa=2, deg C=1, ord C=1, c=y) and reproduce GGV3 sec.5's THREE
      published integers for it: [P_1,Q_1]=x^2, deg(P_1)=10, deg(Q_1)=15.
  R3  derive F2 j=1 = (75,125) with the same chart and COMPUTED reduced
      polygons N(P)=3*Delta', N(Q)=5*Delta', Delta'={(0,0),(3,0),(4,1),(0,5)}.
  RG  the RETRACTION GUARD: final_corner_dictionary() must RAISE at (5,20) and
      (7,21) and must return (l_final,b_final) at (8,28) and (9,24).

Sources (line numbers = local paper_src copies, sha pinned in upstream_facts.json):
  GGHV22 paper_src/2204.14178.tex : Case (8,28) Prop lines 1000-1311; final
    Laurent map line 1229 (phi(x)=x^-1, phi(y)=x^4 y, [phi P,phi Q]=-[P,Q]x^2);
    branch enumeration lines 1020,1073-1090,1132-1186; edge form line 1132
    (y (x^4 y - alpha)^7).  Case (9,27) Prop lines 463-666 (the shared method,
    incl. the (a',b') endpoint table lines 539-562).  Case (7,21) lines
    1313-1396.
  GGV1 1401.1784 Cor 7.4 / Prop 8.2 (root shift + edge form);
  GGV6 Prop 2.5 (Pred_P(1,0) restriction);  vd Essen Prop 10.2.6.
  composite_charts.py : FUSED-CHART LEMMA (Jacobian -x^(l-2), any shears).
  phi_corner4.py / phi_75_125.py / c_series_75_125.py : landed F2 corner data.

Exact sympy only.  Run `python3 polygon_reduction.py` for the full report;
`python3 polygon_reduction_verify.py --quiet` for the PASS/FAIL contract.
"""
from __future__ import annotations

import sympy as sp
from fractions import Fraction
from dataclasses import dataclass, field
from typing import Optional

x, y = sp.symbols("x y")
BAR = "=" * 96


# ---------------------------------------------------------------------------
# 0.  The fused chart -- built on composite_charts.py's proven lemma.
# ---------------------------------------------------------------------------
def fused_jacobian(l):
    """Jacobian of the final Laurent chart (X,Y) = (x^-1, x^l y + shears).

    composite_charts.py STEP 2 proves this equals -x^(l-2) for ANY shear terms
    lambda_i x^(e_i); we CONFIRM it here for the concrete l of the case (a cheap
    exact self-check) rather than rederive the general lemma.
    """
    lam1, e1 = sp.symbols("lambda1 e1", positive=True)
    ls = sp.Integer(l)
    X = x**-1
    Y = x**ls * y + lam1 * x**e1          # a representative shear term
    J = sp.simplify(sp.diff(X, x) * sp.diff(Y, y) - sp.diff(X, y) * sp.diff(Y, x))
    assert sp.simplify(J - (-x**(ls - 2))) == 0, (l, J)
    return J


def invert_vertex(v, l):
    """Push a Newton-polygon vertex through the final Laurent inversion.

    phi(x)=x^-1, phi(y)=x^l y  sends the monomial  x^a y^b  to
    x^(-a) (x^l y)^b = x^(l b - a) y^b, i.e. the vertex (a,b) -> (l b - a, b).
    Shear terms only move interior/lower vertices already accounted for in the
    pre-inversion vertex set, so acting on the vertex set is exact.
    """
    a, b = v
    return (l * b - a, b)


def invert_polygon(vertices, l):
    return sorted({invert_vertex(v, l) for v in vertices})


# ---------------------------------------------------------------------------
# 0b.  THE CHART EXPONENT, AND THE GUARD ON GGV5's FINAL-CORNER DICTIONARY.
#
#      ROOT CAUSE of the 2026-07-26 (5,20) repair (PASSPORT_75_125_REPAIR.md).
#      This repo derived chart data from GGV5's FINAL CHAIN CORNER
#      A_1 = (a\l_final, b_final) by the dictionary  (t, q) = (l_final, b_final).
#      That dictionary is valid EXACTLY on the RETRACTION SHAPE
#
#            b0 == l_chart * (a0 - 1)
#
#      i.e. exactly when A_0 sits on the integer-slope ray through A_0' = (1,0),
#      so the edge {(0,1),(b0,a0)} of the flipped polygon collapses to a
#      VERTICAL face under the inversion.  Checked on all four rows where both
#      sides are independently known (PASSPORT_75_125.md P6):
#
#        corner   A_1        l_final  l_chart   b_final  ord C   dictionary
#        (8,28)   (11\4,7)      4       4         7       7      VALID   (shape holds)
#        (9,24)   (11\3,8)      3       3         8       8      VALID   (shape holds)
#        (7,21)   (11\7,2)      7       3         2       1      BROKEN  (shape fails)
#        (5,20)   (7\5,2)       5       4         2       1      BROKEN  (shape fails)
#
#      (7,21) is the PUBLISHED counterexample: GGV5 gives l_final = 7, while
#      GGHV22 publishes the chart phi_3(y) = y x^3 (2204.14178.tex:1394) and
#      [P,Q] = x, i.e. l_chart = 3, kappa = 1.  (5,20) is the same family shape
#      (F_2 vs F_9: both k=1, both b_final=2, both (m,n)=(j+2,2j+3)), and GGV3's
#      own published reduction of the sibling (50,75) AT THIS CORNER
#      (paper_src/1406.0886_GGV3.tex:1723-1727 -- "[P_1,Q_1]=x^2, deg(P_1)=10
#      and deg(Q_1)=15") forces l_chart = 4, kappa = 2.  Three published
#      integers, all three reproduced by l=4 and all three contradicted by l=5.
#
#      So the dictionary is not retired -- it is GUARDED.  Any consumer that
#      wants (t,q) out of a chain row must come through
#      final_corner_dictionary(), which RAISES off the precondition instead of
#      silently returning a wrong pair.
# ---------------------------------------------------------------------------
class FinalCornerDictionaryError(AssertionError):
    """(t,q) = (l_final,b_final) was applied off its retraction precondition."""


def chart_exponent(a0, b0):
    """l_chart -- the Laurent denominator of the final chart (x^-1, x^l y + shears).

    RULE  l = ceil(b0/a0), the minimal integer with l*a0 >= b0 (the
    first-quadrant condition on the inversion).

    [INFERRED]  Validated on all five published GGHV22 reductions and pinned at
    (5,20) by GGV3's published (50,75) reduction, but located in NO published
    proposition in this form -- CORNER_RESOLVENT.md sec.5.1 correctly records
    that no general dictionary exists in the literature.  Do not cite as
    published.  See PASSPORT_75_125.md N1.
    """
    a0, b0 = int(a0), int(b0)
    assert a0 > 0 and b0 > 0, (a0, b0)
    return -(-b0 // a0)                      # ceil for positive ints


def has_retraction(a0, b0, l=None):
    """The RETRACTION SHAPE test  b0 == l*(a0-1),  l = chart_exponent by default.

    Geometry: the edge {(0,1),(b0,a0)} retracts to a VERTICAL face under the
    inversion (a,b) -> (l*b - a, b) exactly on this equality; off it the reduced
    polygon has no vertical top face.

    *** NOTE THE QUANTIFIER -- this is the whole trap. ***  The test is the
    equality for the l that will ACTUALLY BE USED, not "does SOME integer l
    satisfy it".  At (5,20) some l does -- l = 5 gives 5*(5-1) = 20 -- and that
    coincidence is precisely how l = 5 entered this repo.  With the correct
    l = ceil(20/5) = 4 one has 4*(5-1) = 16 != 20 and the test FAILS, which is
    the right answer.  Never re-derive l from this equality.
    """
    if l is None:
        l = chart_exponent(a0, b0)
    return int(b0) == int(l) * (int(a0) - 1)


def final_corner_dictionary(a0, b0, l_final, b_final, who=""):
    """GUARDED  (t, q) = (l_final, b_final).  Raises off the retraction shape.

    The ONLY sanctioned way to turn a GGV5 chain row's final corner
    (a\\l_final, b_final) into chart data (t, q).  See the block comment above
    for the four-row validity table and the two published counterexamples.
    """
    l = chart_exponent(a0, b0)
    if not has_retraction(a0, b0, l):
        raise FinalCornerDictionaryError(
            "(t,q) = (l_final,b_final) = (%s,%s) REFUSED at A_0=(%s,%s)%s.  The "
            "retraction precondition b0 == l_chart*(a0-1) FAILS: %s != %s*(%s-1) "
            "= %s.  GGV5's final chain corner therefore does NOT carry this "
            "corner's chart data.  Correct chart data: l_chart = ceil(b0/a0) = "
            "%s, kappa = l-2 = %s, and (no vertical top face) C is a MONOMIAL.  "
            "Published instances of exactly this failure: (7,21) [l_final = 7 vs "
            "GGHV22's published chart l = 3] and (5,20) [l_final = 5 vs GGV3's "
            "(50,75) reduction forcing l = 4]."
            % (l_final, b_final, a0, b0, (" for " + who) if who else "",
               b0, l, a0, l * (a0 - 1), l, l - 2))
    return int(l_final), int(b_final)


def corner_chart_data(a0, b0, l_final=None, b_final=None, who=""):
    """The corner signature primitives (t, kappa, deg_C, ord_C) for one corner.

    On the retraction shape the vertical top face survives, deg C = a0, and the
    selected root multiplicity ord C is GGV5's b_final -- taken through the
    guard.  Off it there is no vertical top face, so deg C = 1 and C is a
    MONOMIAL (PASSPORT_75_125.md Q5/S2; the same shape GGHV22 publishes at
    (7,21), where FAMILY_GRAMMAR.md sec.1 already records C = y).
    """
    l = chart_exponent(a0, b0)
    if has_retraction(a0, b0, l):
        _, q = final_corner_dictionary(a0, b0, l_final, b_final, who=who)
        return dict(t=l, kappa=l - 2, deg_C=int(a0), ord_C=q, monomial=False,
                    retraction=True)
    return dict(t=l, kappa=l - 2, deg_C=1, ord_C=1, monomial=True,
                retraction=False)


# ---------------------------------------------------------------------------
# 1.  Data structures for the reduction dossier.
# ---------------------------------------------------------------------------
@dataclass
class Transform:
    name: str
    action: str          # human-readable coordinate action
    jac: str             # Jacobian (as string) and its effect on the bracket
    cite: str


@dataclass
class BranchOption:
    label: str
    followed: bool
    reason: str


@dataclass
class Branch:
    name: str
    question: str
    options: list           # list[BranchOption]
    cite: str = ""

    def followed(self):
        return [o for o in self.options if o.followed]

    def excluded(self):
        return [o for o in self.options if not o.followed]


@dataclass
class Reduction:
    tag: str
    title: str
    A0: tuple
    A0p: tuple
    mn: tuple                       # (m,n) polygon multipliers (P,Q)
    l: int                          # Laurent denominator of the final chart
    transforms: list = field(default_factory=list)
    branches: list = field(default_factory=list)
    # pre-inversion polygons, keyed by output branch label:
    pre_inversion: dict = field(default_factory=dict)
    reduced: dict = field(default_factory=dict)      # filled by compile
    bracket: Optional[str] = None
    kappa: Optional[int] = None
    signature: dict = field(default_factory=dict)
    judgment: list = field(default_factory=list)
    notes: list = field(default_factory=list)


# ---------------------------------------------------------------------------
# 2.  The compiler engine.
# ---------------------------------------------------------------------------
def compile_reduction(red: Reduction) -> Reduction:
    """Fill in the DERIVED quantities: kappa from the Jacobian, and the reduced
    polygons from the pre-inversion vertices via the inversion map."""
    J = fused_jacobian(red.l)                       # = -x^(l-2), confirmed
    red.kappa = red.l - 2
    red.bracket = "x^%d" % red.kappa
    # reduced polygons of every retained output branch:
    red.reduced = {label: {"P": invert_polygon(v["P"], red.l),
                           "Q": invert_polygon(v["Q"], red.l)}
                   for label, v in red.pre_inversion.items()}
    red.notes.insert(0,
        "final chart Jacobian = %s (composite_charts.py lemma, l=%d) "
        "=> [phiP,phiQ] = -[P,Q]*x^(l-2); the pre-inversion bracket is a nonzero "
        "constant (flip contributes -1, every root-shift contributes +1), so "
        "kappa = l-2 = %d and [P,Q] = %s." % (J, red.l, red.kappa, red.bracket))
    return red


# ---------------------------------------------------------------------------
# 3.  CASE (8,28)  --  the published reduction (R1).
#     Transcribed vertex-by-vertex from GGHV22 Prop lines 1000-1311, then the
#     final Laurent inversion is APPLIED by the engine.
# ---------------------------------------------------------------------------
def case_8_28() -> Reduction:
    red = Reduction(
        tag="GGHV_8_28",
        title="Case (8,28)  [R1: the published (72,108)/(8,28) reduction]",
        A0=(8, 28), A0p=(1, 0), mn=(2, 3), l=4)

    red.transforms = [
        Transform("phi1  (flip)", "x <-> y   (phi1(x)=y, phi1(y)=x)",
                  "Jacobian -1 : [P,Q] -> -[P,Q] (constant preserved up to sign)",
                  "GGHV22 line 1012"),
        Transform("phi2  (root shift, depth s in {2,3})",
                  "y -> y + lambda x^(-s)  (clears the lower edge {(1,0),(8,28)}"
                  " after the flip; leading form y(x^4 y - alpha)^7-type)",
                  "Jacobian 1 : bracket-preserving (GGV1 Prop 3.10-class shear)",
                  "GGHV22 lines 1073-1083"),
        Transform("phi3  (edge root shift)",
                  "y -> y + alpha x^(-4)  reduces the edge {(28,8),(0,1)} to "
                  "{(28,8),(24,7)}  (edge form y (x^4 y - alpha)^7)",
                  "Jacobian 1 : bracket-preserving",
                  "GGHV22 line 1132"),
        Transform("phi  (FINAL Laurent inversion)",
                  "x -> x^-1,  y -> x^4 y   (the ONE inversion; l=4)",
                  "Jacobian -x^2 : [phiP,phiQ] = -[P,Q] x^2  => [P,Q]=x^2",
                  "GGHV22 line 1229"),
    ]

    # ---- branch manifest (GGHV22's own case split, every option with reason) --
    red.branches = [
        Branch(
            "Pred_P(1,0) direction (root-shift depth)",
            "which lower-edge slope precedes (1,0)?",
            [BranchOption("(1,-2)", True,
                "FOLLOWED: admissible by GGV6 Prop 2.5; leads to the sub-split "
                "below after the y->y+lambda x^-2 shift."),
             BranchOption("(1,-3)", True,
                "FOLLOWED: admissible by GGV6 Prop 2.5; merges into the same "
                "three final shapes via a y->y+lambda x^-3 shift."),
             BranchOption("(1,-s), s>=4 or s=1", False,
                "EXCLUDED: GGV6 Prop 2.5 restricts Pred_P(1,0) to {(1,-2),(1,-3)}; "
                "deeper/shallower shifts give deg_x P(x,0)<=0, contradicting "
                "vd Essen Prop 10.2.6.")],
            "GGHV22 line 1020"),
        Branch(
            "leading-form factor count after the phi2 shift",
            "how many distinct linear factors z-lambda_i (z=x^2 y) does the "
            "leading form carry? -> determines the residual lower corner",
            [BranchOption("case a) single factor via Pred=(2,-7)", True,
                "FOLLOWED -> lower corners m{(-2,0),(0,0),(28,8),(0,1)}; no "
                "intermediate corner."),
             BranchOption("case b) single factor via Pred=(1,-3)", True,
                "FOLLOWED -> lower corners m{(-3,0),(0,0),(28,8),(0,1)}; the "
                "(-3,0) endpoint differs from a) but the reduced polygon "
                "coincides with a) after inversion (both = sub2)."),
             BranchOption("case c) two factors (x^3y-alpha1)(x^3y-alpha2)", True,
                "FOLLOWED -> extra intermediate corner m(16,4); yields sub1 "
                "after inversion."),
             BranchOption("three distinct factors", False,
                "EXCLUDED for the (9,27) analog by GGV2 Prop 3.12(2) "
                "(a factor of multiplicity s=6 would appear, contradicting "
                "distinctness); for (8,28) the algorithm leaves only the 1- and "
                "2-factor shapes.")],
            "GGHV22 lines 1073-1090, 1085-1089"),
        Branch(
            "opposite vertex (a,b) on the far edge / GGV1 Prop 8.2 exponent k",
            "which endpoint (a,b) and multiplicity k realise the far edge?",
            [BranchOption("(a,b)=(24,7), k=1", True,
                "FOLLOWED: k=1 gives {en(P),en(Q)} = {(-1,0),(2,1)} with a common "
                "direction; the selected far corner is (24,7)."),
             BranchOption("(a,b) in {(17,5),(10,3),(3,1)}, k=1", False,
                "EXCLUDED: these interior endpoints fail the parallel-edge "
                "closure that (24,7) satisfies (same mechanism as the (9,27) "
                "(a',b') divisibility table)."),
             BranchOption("k=2", False,
                "EXCLUDED: with k=2 the edges of P and Q cannot be made parallel "
                "-- geometrically impossible.")],
            "GGHV22 lines 1132-1136"),
        Branch(
            "leading forms proportional?  en_{rho,sigma}(P) ~ en_{rho,sigma}(Q)",
            "does the far-edge analysis need the proportional or "
            "non-proportional sub-case?",
            [BranchOption("en(P) ~ en(Q)", True,
                "FOLLOWED: GGV1 Prop 8.2 applies directly."),
             BranchOption("en(P) not~ en(Q)", True,
                "FOLLOWED: GGHV22 note that BOTH sub-cases land on the same "
                "k=1 conclusion, so the branch is inessential (merged).")],
            "GGHV22 lines 1132-1136 (cf. (9,27) lines 564-566)"),
    ]

    # ---- pre-inversion polygons (GGHV22 lines 1137-1186), scaled 2*/3* -------
    # cases a) and b) share one shape (sub2); case c) adds the (16,4) corner (sub1).
    two = lambda pts: [(2 * a, 2 * b) for a, b in pts]
    three = lambda pts: [(3 * a, 3 * b) for a, b in pts]
    core = [(28, 8), (24, 7)]                      # (a,b): 2*/3* multiplied
    red.pre_inversion = {
        "sub2 (cases a,b)": {
            "P": [(-1, 0), (0, 0)] + two(core),
            "Q": [(2, 1), (0, 0)] + three(core)},
        "sub1 (case c)": {
            "P": [(-1, 0), (0, 0)] + two([(16, 4)] + core),
            "Q": [(2, 1), (0, 0)] + three([(16, 4)] + core)},
    }

    # ---- corner signature consumed by the law -------------------------------
    # final chart l=4 => t=l=4, kappa=t-2=2; the reduced edge form y(x^4 y-a)^7
    # fixes the selected root multiplicity q=7; deg C = a0 = 8 (the (28,8) foot).
    red.signature = dict(t=4, kappa=2, a0=8, q=7,
                         c_of_y="y*(x^4 y - alpha)^7 edge form; residual foot (0,4)",
                         note="R1 target is the polygon match + bracket; the L^(1) "
                              "C-series (C_3 = y^8(y+1)) is the section-4/5 object.")

    red.judgment = [
        "[RETIRED here -- this is the ONE published case] the (8,28) reduction "
        "is written out in full in GGHV22; the compiler reproduces it exactly, "
        "validating the engine."]

    return compile_reduction(red)


# ---------------------------------------------------------------------------
# 4.  F2 corner (5,20)  --  the standard A0'=(1,0) single-Laurent chart.
#     Used by BOTH j=0 (50,75) [R2] and j=1 (75,125) [R3]; only (m,n) differs.
# ---------------------------------------------------------------------------
def _f2_forcing_divisor(a, b, t, kappa, a0, q):
    """Solve the corner-144 forcing ODE to recover the residual divisor c(y).

    a { t c f' - [t(b-a)+kappa+1] c' f } = c^(b-a+1),  c = y^q g, deg g = a0-q.

    TWO REGIMES (2026-07-26 repair).  deg g = a0 - q:
      dg > 0  -- residual g present; the ODE forces g_1..g_{dg-1} = 0, leaves the
                 top coefficient resonant (fixed monic), and g(-1) = 0 selects
                 the (y+1) common root:  g = y^dg + 1.
      dg == 0 -- NO residual at all.  c = y^q is a MONOMIAL, g = 1 is FORCED
                 (monic constant), and the whole common-root-gauge branch is
                 VACUOUS -- there is no free coefficient to gauge and no root to
                 place.  This is the (5,20) regime: deg C = 1, q = 1, C = y.
    Returns (c, g, f, A, N, sig).
    """
    e = b - a + 1
    coef = t * (b - a) + kappa + 1
    rho = (e - 1) * q + 1
    dg = a0 - q
    assert dg >= 0, ("deg g = a0 - q must be >= 0", a0, q)
    A = sp.symbols("A")
    if dg == 0:
        g_sol = sp.Integer(1)                 # FORCED: monic constant, no gauge
    else:
        g_sol = y**dg + 1
    c_sol = y**q * g_sol
    f_gen = A * y**rho * g_sol**e
    resid = sp.expand(a * t * c_sol * sp.diff(f_gen, y)
                      - a * coef * sp.diff(c_sol, y) * f_gen - c_sol**e)
    A_sols = sp.solve(sp.Poly(resid, y).all_coeffs(), A)
    A_sol = A_sols[A] if isinstance(A_sols, dict) else A_sols[0]
    f_sol = sp.expand(f_gen.subs(A, A_sol))
    assert sp.expand(a * t * c_sol * sp.diff(f_sol, y)
                     - a * coef * sp.diff(c_sol, y) * f_sol - c_sol**e) == 0
    # INDEPENDENT uniqueness check: solve the ODE with a general polynomial
    # ansatz (no f = A y^rho g^e shape assumed) and confirm the same f.
    D = e * a0 + rho + 4
    ai = sp.symbols("aa0:%d" % (D + 1))
    fgen2 = sum(ai[i] * y**i for i in range(D + 1))
    r2 = sp.Poly(sp.expand(a * t * c_sol * sp.diff(fgen2, y)
                           - a * coef * sp.diff(c_sol, y) * fgen2 - c_sol**e), y)
    sols2 = sp.solve(r2.all_coeffs(), list(ai), dict=True)
    assert len(sols2) == 1, ("polynomial solution not unique", len(sols2))
    assert sp.expand(fgen2.subs(sols2[0]) - f_sol) == 0, "ansatz != general solve"
    # Phi = f * C^N signature
    N = a * (t * (a + b) - (kappa + 1)) - 2 * b
    Phi = sp.expand(f_sol * c_sol**N)
    deg = sp.degree(Phi, y)
    ordy = min(m[0] for m in sp.Poly(Phi, y).monoms())
    m1, qq, d1 = 0, sp.Poly(Phi, y), sp.Poly(y + 1, y)
    while True:
        qq2, rem = sp.div(qq, d1)
        if not rem.is_zero:
            break
        qq, m1 = qq2, m1 + 1
    sig = (deg, ordy, m1, deg - ordy - m1)
    return c_sol, g_sol, f_sol, A_sol, N, sig


def case_f2(j) -> Reduction:
    """F2 = family F_2 of GGV5 (corner A0=(5,20), A0'=(1,0), GGV5 final (7\\5,2)).
    j=0 -> (m,n)=(2,3) -> (50,75) [R2];  j=1 -> (m,n)=(3,5) -> (75,125) [R3].

    2026-07-26 REPAIR.  l is now DERIVED (chart_exponent(5,20) = 4), not read off
    GGV5's final corner denominator (which is 5 and is NOT the chart exponent --
    see the block comment at 0b).  Consequently kappa = 2, C = y is a monomial,
    deg C = 1 and ord C = 1.  The reduced Newton polygons are now COMPUTED (the
    corner is no longer a "no vertex list" case).
    """
    A0, A0p = (5, 20), (1, 0)
    GGV5_FINAL = (7, 5, 2)                      # A_1 = (p\l_final, b_final), k=1
    l = chart_exponent(*A0)                     # = 4  DERIVED, never literal
    cd = corner_chart_data(*A0, l_final=GGV5_FINAL[1], b_final=GGV5_FINAL[2],
                           who="F2 corner (5,20)")
    assert not cd["retraction"] and cd["monomial"], cd
    t, kappa = cd["t"], cd["kappa"]
    a0, q = cd["deg_C"], cd["ord_C"]            # deg C = 1, ord C = 1  => C = y
    assert (l, t, kappa, a0, q) == (4, 4, 2, 1, 1), (l, t, kappa, a0, q)
    mu = (A0[1] - 1) // A0[0]                   # = 3 : root-shift depth, Pred (1,-mu)
    assert mu + 1 == l, (mu, l)

    if j == 0:
        mn, degs, tag = (2, 3), (50, 75), "F2_j0_50_75"
        title = "F2 j=0 = (50,75)  [R2: GGV3 sec.5's own reduction of this corner]"
    else:
        mn, degs, tag = (3, 5), (75, 125), "F2_j1_75_125"
        title = "F2 j=1 = (75,125)  [R3: the target model]"
    a, b = sorted(mn)                          # reduced pair

    red = Reduction(tag=tag, title=title, A0=A0, A0p=A0p, mn=mn, l=l)
    red.transforms = [
        Transform("phi1  (flip)", "x <-> y",
                  "Jacobian -1 : bracket constant preserved up to sign",
                  "GGV1 Cor 7.4 setup; same flip as GGHV22 line 1012"),
        Transform("phi2  (root shift, depth mu=%d)" % mu,
                  "y -> y + lambda x^(-%d)   (Pred_P(1,0) = (1,-%d), read off the "
                  "flipped lower edge (5,0)--(20,5) of direction 5*(3,1); foot "
                  "(-%d,0))" % (mu, mu, mu),
                  "Jacobian 1 : bracket-preserving",
                  "GGV1 Cor 7.4 / Prop 8.2; PASSPORT_75_125.md r4"),
        Transform("phi  (FINAL Laurent inversion)",
                  "x -> x^-1,  y -> x^%d y   (the ONE inversion; l = "
                  "chart_exponent(5,20) = ceil(20/5) = %d)" % (l, l),
                  "Jacobian -x^%d : [phiP,phiQ] = -[P,Q] x^%d  => [P,Q]=x^%d"
                  % (kappa, kappa, kappa),
                  "fused-chart lemma, composite_charts.py; cf. GGHV22 line 1229"),
    ]

    red.branches = [
        Branch(
            "chart exponent l  --  THE REPAIRED BRANCH (2026-07-26)",
            "where does the Laurent denominator l come from?",
            [BranchOption("l = chart_exponent(5,20) = ceil(20/5) = 4", True,
                "FOLLOWED.  l is the minimal integer with l*a0 >= b0 (INFERRED "
                "rule; validated on all five published GGHV22 reductions).  It "
                "reproduces GGV3's THREE published integers for the sibling "
                "(50,75) at this same corner ([P_1,Q_1]=x^2, deg 10, deg 15; "
                "1406.0886_GGV3.tex:1723-1727)."),
             BranchOption("l = l_final = 5, read off GGV5's final corner (7\\5,2)",
                False,
                "EXCLUDED, and this is the error this file used to make.  The "
                "dictionary (t,q) = (l_final,b_final) holds only on the "
                "retraction shape b0 == l*(a0-1), which FAILS here (20 != 4*4); "
                "final_corner_dictionary() now raises on it.  l = 5 predicts "
                "[P,Q] = x^3 and reduced degrees (20,30) for (50,75), "
                "contradicting all three of GGV3's published integers.  The same "
                "dictionary is refuted in print at (7,21): l_final = 7 vs "
                "GGHV22's published chart l = 3.")],
            "polygon_reduction.py sec.0b; PASSPORT_75_125.md S7-S10, P6"),
        Branch(
            "retraction / vertical top face",
            "does the edge {(0,1),(20,5)} retract to a vertical face?",
            [BranchOption("NO retraction: 20 != 4*(5-1) = 16", True,
                "FOLLOWED (a computation, not a choice): has_retraction(5,20) is "
                "False for l = 4.  Hence no vertical top face, hence deg C = 1, "
                "hence C = y is a MONOMIAL with ord C = 1.  Same shape as "
                "(7,21), where GGHV22 publishes exactly this."),
             BranchOption("retraction (deg C = a0 = 5, ord C = b_final = 2)", False,
                "EXCLUDED: that is the (8,28)/(9,24) shape.  It was assumed here "
                "before 2026-07-26 and is what produced C = y^2(y^3+1).")],
            "PASSPORT_75_125.md Q5 / S2"),
        Branch(
            "chart class",
            "is the reduction chart determined?",
            [BranchOption("standard single-Laurent chart (X,Y)=(x^-1, x^4 y + shears)",
                True,
                "FOLLOWED and FORCED: A0'=(1,0), length-1 chain => exactly one "
                "inversion; composite_charts.py proves Jacobian -x^(l-2) for ANY "
                "shears, so kappa = l-2 = 2 is not a choice once l is fixed.  "
                "Same chart class GGHV22 uses for (8,28)/(9,27)/(7,21)."),
             BranchOption("double-inversion chart (kappa=l2-l1)", False,
                "EXCLUDED: would require a second inversion the length-1 chain "
                "never performs (composite_charts.py STEP 2, heuristic killed).")],
            "composite_charts.py STEP 2; phi_corner4.py STEP 2"),
        Branch(
            "two-factor split corner",
            "can the leading form carry two distinct linear factors?",
            [BranchOption("single factor: gcd(5,20)=5, a0/gcd = 1", True,
                "FOLLOWED: the GGV1 Cor-7.4 multiplicity is gcd(a0,b0) = 5 and "
                "en(R) = (b0,a0)/5 = (4,1) is primitive, so R has z-degree "
                "a0/5 = 1 -- R = x(x^3 y - alpha) has ONE linear factor."),
             BranchOption("two distinct factors (the (8,28) case-c branch)", False,
                "EXCLUDED: z-degree 1 admits no second factor, so the "
                "(8,28)-style extra intermediate corner does not exist here.")],
            "PASSPORT_75_125.md Q4 / S3 / R-rule"),
        Branch(
            "common-root gauge for the residual g  --  now VACUOUS",
            "how is the free resonant coefficient of g fixed?",
            [BranchOption("deg g = a0 - q = 1 - 1 = 0  =>  g = 1 FORCED", True,
                "FOLLOWED.  There is no residual polynomial: g is a monic "
                "constant, so there is NO free coefficient to gauge and NO root "
                "to place.  The gauge branch that was REOPENED on 2026-07-24 "
                "(unramified g = y^3+1 vs a ramified double-root g) is not "
                "resolved -- it is DISSOLVED: it presupposed deg g = 3, which "
                "came from deg C = a0 = 5, which came from the retracted-shape "
                "assumption that this corner does not satisfy."),
             BranchOption("ramified / unramified gauge on a cubic g", False,
                "EXCLUDED: presupposes deg g = 3.  With C = y a monomial there "
                "is no g to ramify.")],
            "polygon_reduction.py sec.0b; PASSPORT_75_125.md Q5"),
    ]

    # forcing-layer residual divisor (recovered exactly)
    c_sol, g_sol, f_sol, A_sol, N, phi_sig = _f2_forcing_divisor(a, b, t, kappa, a0, q)
    assert sp.expand(c_sol - y) == 0, ("C must be the monomial y", c_sol)
    red.signature = dict(t=t, kappa=kappa, a0=a0, q=q,
                         c_of_y=str(sp.factor(c_sol)),
                         g=str(sp.factor(g_sol)),
                         reduced_pair=(a, b), degs=degs,
                         N=int(N), phi_signature=tuple(int(s) for s in phi_sig),
                         ggv5_final_corner="(%d\\%d,%d) k=1  [CHAIN data; NOT the "
                                           "chart -- see sec.0b]" % GGV5_FINAL,
                         chart_exponent_rule="l = ceil(b0/a0) = %d  [INFERRED]" % l)

    # ---- pre-inversion polygon (now EXPLICIT; PASSPORT_75_125.md sec.2) -------
    #   Delta  = {(0,0),(1,0),(5,20),(0,5)}                            [r1]
    #   flip   = {(0,0),(0,1),(20,5),(5,0)}                            [r3]
    #   lower edge (5,0)--(20,5) has direction 5*(3,1) => Pred (1,-3), foot (-3,0)
    #   pre-inversion = {(-3,0),(0,0),(0,1),(20,5)}                     [r7]
    #   inversion (i,j) -> (4j - i, j)  =>  Delta' = {(0,0),(3,0),(4,1),(0,5)}
    #   Prop 8.2(2) en-split EXCLUDED (Q3) => PROPORTIONAL: N(P)=m*Delta',
    #   N(Q)=n*Delta'.
    core = [(-mu, 0), (0, 0), (0, 1), (A0[1], A0[0])]
    red.pre_inversion = {
        "standard (proportional, Prop 8.2(1))": {
            "P": [(mn[0] * i, mn[0] * jj) for i, jj in core],
            "Q": [(mn[1] * i, mn[1] * jj) for i, jj in core]},
    }
    red = compile_reduction(red)

    # the reduced polygons must be the (m,n)-scaled Delta'
    Dp = [(0, 0), (3, 0), (4, 1), (0, 5)]
    got = red.reduced["standard (proportional, Prop 8.2(1))"]
    assert set(got["P"]) == {(mn[0] * i, mn[0] * jj) for i, jj in Dp}, got["P"]
    assert set(got["Q"]) == {(mn[1] * i, mn[1] * jj) for i, jj in Dp}, got["Q"]

    # judgment resolution
    common = [
        "[RETIRED at the polygon layer] the chart is the standard "
        "single-Laurent A0'=(1,0) chart and l is DERIVED (l = ceil(b0/a0) = 4), "
        "not read off GGV5's final corner.  kappa = l-2 = 2 follows from the "
        "fused-chart Jacobian, so it is derived, not assumed.  The reduced "
        "Newton polygons are COMPUTED: N(P) = m*{(0,0),(3,0),(4,1),(0,5)}, "
        "N(Q) = n*{(0,0),(3,0),(4,1),(0,5)}.",
        "[DISSOLVED, was REOPENED 2026-07-24] the residual common-root gauge is "
        "no longer a branch: deg g = a0 - q = 0, so g = 1 is forced and C = y is "
        "a monomial.  The reopened unramified-vs-ramified question presupposed a "
        "cubic g, which presupposed deg C = 5, which presupposed the retracted "
        "shape this corner does not have.",
        "[surviving, forcing layer only] forcing-polynomial identification "
        "(corner-144 correspondence, audited only for (72,108)) is a SEPARATE, "
        "non-polygon judgment and is untouched by this compiler.",
    ]
    if j == 0:
        red.judgment = common + [
            "[EXTERNAL CONTROL] GGV3 sec.5 performs this exact reduction and "
            "publishes [P_1,Q_1] = x^2, deg(P_1) = 10, deg(Q_1) = 15 "
            "(1406.0886_GGV3.tex:1723-1727).  All three are reproduced here; all "
            "three are contradicted by the superseded l = 5."]
    else:
        red.judgment = common + [
            "[UNCONDITIONAL at the polygon layer] F2 j=1 uses the IDENTICAL chart "
            "as F2 j=0: same corner, same A0'=(1,0), same length-1 chain, same "
            "final map (x^-1, x^4 y).  Only the (m,n) multiplier changes "
            "(2,3)->(3,5), which scales the polygon but does not touch the chart. "
            "So the j=0 external control (GGV3) transfers to j=1."]
    return red


# ---------------------------------------------------------------------------
# 5.  Registry + report.
# ---------------------------------------------------------------------------
def all_reductions():
    return [case_8_28(), case_f2(0), case_f2(1)]


def _fmt_polys(d):
    out = []
    for label, pq in d.items():
        out.append("    [%s]" % label)
        out.append("      N(P) = %s" % (pq["P"],))
        out.append("      N(Q) = %s" % (pq["Q"],))
    return "\n".join(out)


def report(red: Reduction) -> str:
    L = [BAR, red.title, BAR,
         "  input:  A0 = %s   A0' = %s   (m,n) = %s   Laurent l = %d"
         % (red.A0, red.A0p, red.mn, red.l),
         "",
         "  (1) TRANSFORMATION SEQUENCE"]
    for i, tr in enumerate(red.transforms, 1):
        L.append("    %d. %-28s %s" % (i, tr.name, tr.action))
        L.append("       %s   [%s]" % (tr.jac, tr.cite))
    L += ["",
          "  (2) TRANSFORMED BRACKET",
          "      kappa = l - 2 = %d   =>   [P,Q] = %s   (DERIVED from Jacobian "
          "-x^(l-2))" % (red.kappa, red.bracket)]
    if red.reduced:
        L += ["", "  (3) REDUCED NEWTON POLYGONS (computed via (a,b)->(%d b - a, b))"
              % red.l, _fmt_polys(red.reduced)]
    else:
        L += ["", "  (3) REDUCED NEWTON POLYGONS: corner-signature case "
              "(no published vertex list for this corner; see (5))"]
    L += ["", "  (4) BRANCH MANIFEST (%d branches, %d options total)"
          % (len(red.branches),
             sum(len(b.options) for b in red.branches))]
    for b in red.branches:
        L.append("    * %s -- %s   [%s]" % (b.name, b.question, b.cite))
        for o in b.options:
            L.append("        [%s] %s" % ("FOLLOW " if o.followed else "EXCLUDE",
                                          o.label))
            L.append("               %s" % o.reason)
    L += ["", "  (5) CORNER SIGNATURE (consumed by the corner law)"]
    for k, v in red.signature.items():
        L.append("      %-14s = %s" % (k, v))
    L += ["", "  JUDGMENT RESOLUTION"]
    for jt in red.judgment:
        L.append("      %s" % jt)
    for n in red.notes:
        L += ["", "  note: " + n]
    return "\n".join(L)


def main():
    for red in all_reductions():
        print(report(red))
        print()


if __name__ == "__main__":
    main()
