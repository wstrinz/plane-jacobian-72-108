#!/usr/bin/env python3
"""family_grammar_verify.py  (NEW) -- exact independent checker for FAMILY_GRAMMAR.

Re-verifies, with exact SymPy only and INDEPENDENTLY of family_grammar.py's
solvers, every claim of the family-grammar sweep:

  A. corner-data identities for all 17 families (kappa=t-2, dg=a0-q, r=dg-1,
     gap=(q-1)-a0/t, te-coef=1) and the collapse identity y G'-dg G=-dg.
  B. the CLOSED-FORM THEOREM: pure ansatz solves the ODE for all j  <=> gap=0,
     with A=-1/(a dg); checked BOTH symbolically-in-j (opaque-G reduction) AND
     by direct full-ODE substitution at j = 0,1,2,3.
  C. CLOSED-FORM (cofactor, r=0): the degree-gap unit cofactor is re-solved from
     the FULL ODE at j = 0,1,2,3, confirmed a UNIT (u(0),u(-1),lead != 0), and
     the ODE residual is exactly 0.
  D. RUNG-STRUCTURED: the mu=dg ramified rung f = y^rho (y+1)^(dg e-(dg-1)) u
     (deg u = gap+r) re-solved from the FULL ODE at j=0..3; the four PHI_F7
     landed f-polynomials substituted DIRECTLY into the ODE (=0); mu=1 real
     parity noted (never used as a complex kill).
  E. the mu-graded corner law reproduces EVERY landed derived point and every
     published mu-rung signature exactly; cof = deg-ord-mult identically.

--quiet suppresses the per-check log; exit 0 iff all pass.  Does NOT edit
run_tests.sh or any existing file.
"""
import sys
from fractions import Fraction
import sympy as sp

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

y = sp.symbols("y", positive=True)
j = sp.symbols("j", nonnegative=True, integer=True)

# ---------------------------------------------------------------------------
# family table (re-transcribed independently; a=min(m,n), b=max, linear in j)
#   name, t, a0, q, a(j)=(c0,c1), b(j)=(c0,c1), A0', k
# ---------------------------------------------------------------------------
FAM = [
    ("F1",  4, 4, 3, (3, 2), (4, 3),  (1, 0), 1),
    ("F2",  5, 5, 2, (2, 1), (3, 2),  (1, 0), 1),
    ("F3",  5, 5, 3, (2, 3), (3, 4),  (1, 0), 1),
    ("F4",  5, 5, 3, (3, 2), (16, 12),(1, 0), 2),
    ("F5",  5, 5, 4, (5, 4), (9, 7),  (1, 0), 1),
    ("F6",  5, 5, 4, (7, 6), (18, 16), (1, 0), 2),  # CORRECTED 2026-07-24 (GGV5 F6 typo base (4,10) gcd=2 -> coprime a=6j+7,b=16j+18 base (7,18); CHAIN_SURVEY.md)
    ("F7",  3, 6, 4, (2, 1), (7, 4),  (1, 0), 1),
    ("F8",  3, 6, 5, (3, 2), (7, 5),  (1, 0), 1),
    ("F9",  7, 7, 2, (2, 1), (3, 2),  (1, 0), 1),
    ("F10", 7, 7, 3, (4, 3), (7, 5),  (1, 0), 1),
    ("F11", 7, 7, 3, (2, 1), (5, 3),  (1, 0), 2),
    ("F12", 4, 8, 5, (3, 2), (7, 5),  (2, 0), 1),
    ("F13", 3, 9, 7, (2, 1), (13, 7), (2, 0), 1),
    ("F14", 3, 9, 4, (2, 1), (7, 4),  (1, 0), 1),
    ("F15", 3, 9, 5, (3, 2), (7, 5),  (1, 0), 1),
    ("F16", 3, 9, 7, (3, 4), (5, 7),  (1, 0), 1),
    ("F17", 3, 9, 8, (2, 5), (3, 8),  (1, 0), 1),
]

LANDED = {
    # REPAIRED 2026-07-26 (2adb92a chart repair); were 36/(189,75,38,76)
    # and 98/(504,201,101,202) off the superseded (5,20) chart.
    ("F2", 0): (28, (30, 30, 0, 0)),
    ("F2", 1): (77, (80, 80, 0, 0)),
    # 2026-07-27 SECOND CHART REPAIR (nine more rows refused).  Three landed
    # points move; each repaired ord_y agrees with the PROVED bridge identity
    # a*q*M - H, and bridge_generality.py MUT F displaces two of them by a wholly
    # independent route (F1 51<-205, F9 22<-107).
    ("F9", 0): (20, (22, 22, 0, 0)),
    ("F14", 0): (36, (375, 165, 42, 168)),
    ("F1", 0): (49, (51, 51, 0, 0)),
    ("F7", 0): (36, (250, 165, 83, 2)),
    # replacements derived on RETRACTING corners, so the gap>0 regimes no longer
    # rest on a refused corner (PHI_F14.md / PHI_F7.md, 2026-07-27):
    ("F17", 0): (20, (195, 169, 22, 4)),
    ("F8", 0): (70, (448, 371, 75, 2)),
    ("F15", 0): (70, (672, 371, 297, 4)),
    # F3 j=0 REPAIRED 2026-07-26 (second (5,20) repair); was 36/(189,112,75,2).
    # (75,50) is the (m,n)-swap of F2 j=0's (50,75) at the SAME corner (5,20), so
    # the reduced pair {min,max} = {2,3} and the chart (4,2,1,1) coincide and the
    # two rows MUST agree.  Independently derived in A10a-A10c below.
    ("F3", 0): (28, (30, 30, 0, 0)),
    ("F10", 0): (110, (114, 114, 0, 0)),
    ("F16", 0): (56, (528, 407, 117, 4)),
}

# explicit PHI_F7 ramified f-polynomials (for direct full-ODE substitution)
LANDED_F = {
    "F7":  sp.Rational(1, 10)   * y**21 * (y + 1)**11 * (9*y**2 + 3*y - 1),
    "F16": sp.Rational(1, 330)  * y**15 * (y + 1)**5  *
           (243*y**4 + 81*y**3 - 27*y**2 + 15*y - 10),
    # 2026-07-27 replacements, all on RETRACTING corners:
    "F15": -sp.Rational(1, 105)  * y**21 * (y + 1)**17 *
           (243*y**4 + 405*y**3 + 135*y**2 - 15*y + 5),
    "F17": -sp.Rational(1, 910)  * y**9  * (y + 1)**2  *
           (243*y**4 - 81*y**3 + 54*y**2 - 42*y + 35),
    "F8":  -sp.Rational(1, 42)   * y**21 * (y + 1)**5  * (9*y**2 - 3*y + 2),
}

# F3's PHI_F7 polynomial, moved out of LANDED_F on 2026-07-26.  It is a correct
# solution of the ODE of the corner data (t,kappa,q,dg) = (5,3,3,2) -- the
# computation was never wrong -- but (5,20) does not HAVE that corner data, so it
# is not a fact about F3.  A10d checks both halves of that sentence.
# 2026-07-27: F10's rung polynomial joins F3's, for the same reason.  It exactly
# solves the mu=4 ramified ODE of the corner data (t,kappa,q,dg) = (7,5,3,4), so
# the arithmetic was never wrong -- but the corner (7,21) does not HAVE that chart
# data.  Its chart is (3,1,1,1): C = y, no residual.  GGHV22 2204.14178.tex:1394
# PUBLISHES l = 3 and [P,Q] = x at (7,21), so this one is refuted in print.
SUPERSEDED_F = {
    "F3":  sp.Rational(1, 42)   * y**4  * (y + 1)**3  * (25*y**2 + 15*y - 3),
    "F10": sp.Rational(1, 3740) * y**10 * (y + 1)**13 *
           (2401*y**4 + 5831*y**3 + 4165*y**2 + 595*y - 85),
}

# 2026-07-27: BOTH published mu-rung witnesses are RETIRED -- see the A11 TRIPWIRE.
# F12's rungs (ZETA_TAIL.md) sit at (8,24) and F10's (MU_RUNGS_F10.md) at (7,21);
# both corners are guard-REFUSED, so dg = 0 and there is no residual g for (y+1)
# to divide.  The rungs are VACUOUS, not unverified.
MU_RUNGS = {}

SUPERSEDED_MU_RUNGS = {
    ("F12", 1): (814, 506, 102, 206),
    ("F12", 2): (814, 506, 203, 105),
    ("F12", 3): (814, 506, 304, 4),
    ("F10", 2): (1917, 820, 547, 550),
    ("F10", 4): (1917, 820, 1093, 4),
}


# PER-FAMILY CHART REPAIRS (2026-07-26), derived INDEPENDENTLY of
# family_grammar.py's own REPAIRS table -- check A9 below asserts the two agree.
# A family whose corner fails the retraction shape b0 = t(a0-1) has
# deg_y(C)/ord_y(C) different from (a0,q), because GGV5's final-corner dictionary
# is what made them equal and it is invalid off that shape.  F2's corner (5,20)
# fails it (20 != 4*4), so t = 4 and C = y is a MONOMIAL: deg C = ord C = 1.
#
# F3 added 2026-07-26.  It sits at the SAME corner (5,20) with the SAME
# A_0' = (1,0) and the same chain length 1.  chart_exponent, kappa = l-2 and the
# vertical-top-face test are functions of A_0 ALONE, so the chart data is a
# CORNER property: F3 gets the identical (4,2,1,1).  The old table's disagreement
# -- ord C = 2 for F2 but ord C = 3 for F3 at one corner -- was itself proof the
# per-row datum b_final was being misread as chart data.  See A10a-A10f.
#
# 2026-07-27: the hand-written table is replaced by a DERIVATION from the guard,
# built here from this file's OWN corner transcription (so it remains an
# independent copy -- check A9b still compares it to family_grammar.py's).  Every
# row whose corner fails the retraction shape is repaired, not just the two at
# (5,20).  What licensed completing it is the PROVED bridge identity
# ord_y(Phi) = a*q*M - H (BRIDGE_GENERALITY.md): each repaired row now has a
# target this file does not compute.  See A10f for the updated scope note.
import polygon_reduction as _prg                                 # noqa: E402

# (A_0, l_final, b_final) per family -- GGV5 CHAIN data, transcribed here
# independently of family_grammar.CORNERS.
CORNERS = {
    "F1":  ((4, 12), 4, 3),  "F2":  ((5, 20), 5, 2),  "F3":  ((5, 20), 5, 3),
    "F4":  ((5, 20), 5, 3),  "F5":  ((5, 20), 5, 4),  "F6":  ((5, 20), 5, 4),
    "F7":  ((6, 15), 3, 4),  "F8":  ((6, 15), 3, 5),  "F9":  ((7, 21), 7, 2),
    "F10": ((7, 21), 7, 3),  "F11": ((7, 21), 7, 3),  "F12": ((8, 24), 4, 5),
    "F13": ((9, 21), 3, 7),  "F14": ((9, 24), 3, 4),  "F15": ((9, 24), 3, 5),
    "F16": ((9, 24), 3, 7),  "F17": ((9, 24), 3, 8),
}

REPAIRS = {}
for _nm, (_A0, _lf, _bf) in CORNERS.items():
    _cd = _prg.corner_chart_data(_A0[0], _A0[1], l_final=_lf, b_final=_bf,
                                 who="family_grammar_verify " + _nm)
    if not _cd["retraction"]:
        REPAIRS[_nm] = dict(t=_cd["t"], degC=_cd["deg_C"], ordC=_cd["ord_C"],
                            gap0=0)


def bridge_ordy(a, b, t, kappa, ordC):
    """ord_y(Phi) = a*q*M - H.  PROVED in BRIDGE_GENERALITY.md -- an INDEPENDENT
    target: nothing in this file or in family_grammar.py derives it."""
    s_ = a + b
    return a * ordC * (t * s_ - (kappa + 1)) - (ordC * s_ - 1)


def data(name, t, a0, q, ac, bc, A0p, k, degC=None, ordC=None, gap0=None):
    kappa = t - 2
    degC = a0 if degC is None else degC
    ordC = q if ordC is None else ordC
    # dg is the RESIDUAL's degree (g = y^dg+1), i.e. deg C - ord C.  That equals
    # a0 - q only when the dictionary holds; dg = 0 means g is constant, so there
    # is NO residual -- the repaired (5,20) case.
    dg = degC - ordC
    r = dg - 1
    gap = (Fraction(q - 1) - Fraction(a0, t)) if gap0 is None else Fraction(gap0)
    # Does GGV5's final-corner dictionary hold here?  It does iff we did not have
    # to override the C-series data, i.e. iff the corner retracts.  Three of the
    # structural identities below are CONSEQUENCES of that dictionary and are
    # asserted only where it is valid.
    dict_valid = (degC == a0 and ordC == q and gap0 is None)
    a = ac[0] + ac[1] * j
    b = bc[0] + bc[1] * j
    e = b - a + 1
    coef = t * (b - a) + kappa + 1
    rho = sp.expand((e - 1) * ordC + 1)
    N = sp.expand(a * (t * (a + b - 1) + 1) - 2 * b)
    if dg == 0:
        # C is a MONOMIAL: g = y^0+1 is constant, so the pure ansatz
        # f = -1/(a*dg) y^rho (y^dg+1)^e is UNDEFINED -- its normalising constant
        # is -1/(a*0) -- and the collapse identity y g' - dg g = -dg degenerates
        # to 0 = 0, pinning nothing.  The landed signatures still reproduce (see
        # section E), but the MECHANISM the grammar credits them to does not
        # exist.  This class records exactly that, and no more.
        cls = "CHART-DEGENERATE"
    elif gap == 0:
        cls = "PURE"
    elif r == 0:
        cls = "COFACTOR"
    elif gap > 0:
        cls = "RUNG"
    else:
        cls = "IRREGULAR"
    return dict(name=name, t=t, kappa=kappa, a0=a0, q=q, dg=dg, r=r, gap=gap, k=k,
                a=a, b=b, e=e, coef=coef, rho=rho, N=N, cls=cls, A0p=A0p,
                degC=degC, ordC=ordC, dict_valid=dict_valid)


def full_ode_residual(D, jv, f, g=None):
    """exact residual of a(t c f' - coef c' f) - c^e at integer j=jv (all
    exponents integers -> exact, no symbolic-power ambiguity).  c = y^q g with
    the BRANCH residual g: g = y^dg+1 (pure/cofactor) or g = (y+1)^dg (ramified).

    2026-07-26: c = y^(ord C) * g, so the monomial exponent is ord_y(C), NOT the
    chain datum q.  A no-op on every family where the final-corner dictionary
    holds (there ord C == q); on the repaired (5,20) rows it is the difference
    between c = y (right) and c = y^2 or y^3 (the superseded reading).
    """
    s = {j: jv}
    a = int(D["a"].subs(s)); b = int(D["b"].subs(s))
    e = b - a + 1
    q = D["ordC"]; dg = D["dg"]; t = D["t"]; coef = int(D["coef"].subs(s))
    if g is None:
        # dg == 0 forces g = 1, the MONIC CONSTANT -- NOT y^0+1 = 2.  The formula
        # y^dg+1 is the monic degree-dg residual with the (y+1) common root placed,
        # and it is only defined for dg >= 1; at dg = 0 there is no residual and no
        # root to place, so the monic constant is 1 (polygon_reduction
        # _f2_forcing_divisor's dg==0 regime says exactly this).  Evaluating the
        # formula anyway would put a stray factor 2 into c and make the repaired
        # (5,20) ODE unsolvable -- a latent bug until F3 first reached this line.
        g = sp.Integer(1) if dg == 0 else y**dg + 1
    c = y**q * g
    resid = a * (t * c * sp.diff(f, y) - coef * sp.diff(c, y) * f) - c**e
    return sp.expand(resid)


def sig(D, N, mu):
    e, a0, q, rho, gap, r = D["e"], D["a0"], D["q"], D["rho"], D["gap"], D["r"]
    degC, ordC = D.get("degC", a0), D.get("ordC", q)
    pure = e * degC - ordC + 1
    res = pure + gap
    deg = sp.expand(res + N * degC)
    ordy = sp.expand(rho + N * ordC)
    if D["dg"] == 0:
        # C monomial => no residual => no multiplicity and no cofactor.  This is
        # why the repaired (5,20) signatures end in two zeros.
        return deg, ordy, sp.Integer(0), sp.Integer(0)
    mult = sp.expand(mu * (e + N) - (mu - 1))
    cof = sp.expand(gap + r * (e + N) - (mu - 1) * (e + N - 1))
    return deg, ordy, mult, cof


def _mk(row):
    rep = dict(REPAIRS.get(row[0], {}))
    t = rep.pop("t", None)
    if t is not None:
        row = (row[0], t) + tuple(row[2:])
    return data(*row, **rep)


DAT = {row[0]: _mk(row) for row in FAM}

# ===========================================================================
# A. corner-data identities + collapse
# ===========================================================================
for dg in range(1, 8):
    ok("collapse identity dg=%d: y G' - dg G = -dg" % dg,
       sp.expand(y * sp.diff(y**dg + 1, y) - dg * (y**dg + 1)) == -dg)

for nm, D in DAT.items():
    ok("%s kappa=t-2" % nm, D["kappa"] == D["t"] - 2)
    ok("%s dg=degC-ordC, r=dg-1  (dg=a0-q holds only where the corner RETRACTS; F2 (5,20) does not, so deg C=ord C=1 and dg=0)" % nm,
       D["dg"] == D["degC"] - D["ordC"] and D["r"] == D["dg"] - 1)
    if D["dict_valid"]:
        ok("%s gap=(q-1)-a0/t" % nm,
           D["gap"] == Fraction(D["q"] - 1) - Fraction(D["a0"], D["t"]))
    else:
        ok("%s gap formula VOID (corner does not retract, so the final-corner "
           "dictionary that derives gap=(q-1)-a0/t does not hold); gap supplied "
           "directly as %s" % (nm, D["gap"]), True)
    # te - coef = t - kappa - 1 = 1, identically in j
    ok("%s te-coef=1 (identically in j)" % nm,
       sp.expand(D["t"] * D["e"] - D["coef"]) == 1)
    # t*rho - coef*q = t-(kappa+1)q, constant in j
    if D["dict_valid"]:
        trc = sp.expand(D["t"] * D["rho"] - D["coef"] * D["q"])
        ok("%s t*rho-coef*q = t-(kappa+1)q (const in j)" % nm,
           trc == D["t"] - (D["kappa"] + 1) * D["q"])
    else:
        ok("%s t*rho-coef*q identity VOID (it is a consequence of the same "
           "dictionary; rho is built from ord C = %s, not from q = %s)"
           % (nm, D["ordC"], D["q"]), True)

# census
census = {"PURE": [], "COFACTOR": [], "RUNG": [], "IRREGULAR": [],
          "CHART-DEGENERATE": []}
for nm, D in DAT.items():
    census[D["cls"]].append(nm)
census.setdefault("CHART-DEGENERATE", [])
# 2026-07-27: the census is REDRAWN by the completed repair.  Every one of the
# eleven guard-refused rows is CHART-DEGENERATE (C = y, dg = 0), which leaves the
# three genuine mechanism classes populated only by rows whose corner RETRACTS --
# exactly the rows where those mechanisms describe an existing residual.
ok("census PURE = {F14} only.  F2 left on 2026-07-26 and F9 on 2026-07-27, both "
   "because their corners are guard-refused; F14's corner (9,24) retracts, so it "
   "is the sole surviving family for which the pure closed form "
   "f = -1/(a*dg) y^rho (y^dg+1)^e is a statement about an existing residual",
   set(census["PURE"]) == {"F14"})
ok("census CHART-DEGENERATE = the ELEVEN guard-refused rows {F1-F6, F9-F13}: every "
   "one sits on a corner off the retraction shape, so C = y is a monomial, "
   "dg = deg C - ord C = 0, and the pure closed form is UNDEFINED there "
   "(A = -1/(a*0)).  The landed signatures still reproduce -- section E -- but the "
   "closed-form MECHANISM does not exist and is not claimed.  The class is exactly "
   "the set of corners in this table off the retraction shape, and equals REPAIRS",
   set(census["CHART-DEGENERATE"]) == {"F1", "F2", "F3", "F4", "F5", "F6",
                                       "F9", "F10", "F11", "F12", "F13"}
   and set(census["CHART-DEGENERATE"]) == set(REPAIRS))
ok("census COFACTOR = {F8,F17} -- F1, F5 and F6 left on 2026-07-27 (refused "
   "corners (4,12) and (5,20)).  Both survivors retract, and F17 (gap=4) and F8 "
   "(gap=2) are the two fresh derivations that replace F1 in this regime",
   set(census["COFACTOR"]) == {"F8", "F17"})
ok("census RUNG = {F7,F15,F16} -- F3 left on 2026-07-26 and F4, F10, F11, F12, F13 "
   "on 2026-07-27, all for one reason: a 'ramified rung' presupposes a residual g "
   "of degree dg >= 1, and a refused corner has dg = 0.  All three survivors "
   "retract; F15 (dg=4) is the fresh derivation replacing F10",
   set(census["RUNG"]) == {"F7", "F15", "F16"})
ok("census: no length-1 IRREGULAR", census["IRREGULAR"] == [])
ok("census: the three mechanism classes together are EXACTLY the six retracting "
   "rows, and the degenerate class is exactly the eleven refused ones -- so class "
   "membership is now a function of the CORNER, as it must be",
   set(census["PURE"]) | set(census["COFACTOR"]) | set(census["RUNG"])
   == {"F7", "F8", "F14", "F15", "F16", "F17"})

# ---------------------------------------------------------------------------
# A10. F3 IS NOW REPAIRED (2026-07-26, second (5,20) repair).
#
# This block used to be a TRIPWIRE asserting F3 was still unrepaired, with the
# recorded reason "no repaired landed target exists, so no verifiable fix is
# possible yet".  That reason was TOO CONSERVATIVE and is now retracted.  It is
# replaced, in the same spirit, by checks that assert the NEW state together with
# the target that was verified -- from both sides wherever two things can be made
# to agree instead of one thing asserted twice.
#
# WHAT WAS ESTABLISHED
#
#   1. Chart data is a property of the CORNER, not of the family.  l =
#      chart_exponent(a0,b0) = ceil(b0/a0), kappa = l-2, and the vertical-top-face
#      test has_retraction(a0,b0,l) -- which is what decides whether C is a
#      monomial -- are functions of A_0 alone.  A family row contributes only
#      (m,n), which enters the corner law solely through the UNORDERED reduced
#      pair {min(m,n), max(m,n)}.  A10a checks the first half by calling
#      polygon_reduction.corner_chart_data on BOTH GGV5 rows and requiring
#      bit-identical output; A10b checks the second half.
#
#   2. q = 3 was never chart data for F3.  The old table asserted ord C = 2 (F2)
#      and ord C = 3 (F3) AT ONE AND THE SAME CORNER (5,20).  ord_y(C) is a corner
#      invariant, so two values at one corner is a contradiction, visible without
#      consulting any paper -- and it is exactly the fingerprint of reading q off
#      GGV5's per-row final chain corner b_final (2 in F2's (7\5,2), 3 in F3's
#      (8\5,3)).  A10c makes that contradiction an explicit check.
#
# THE VERIFIED TARGET:  F3 j=0 = (75,50):  N = 28,  signature (30,30,0,0).
#
#   (i)  polygon_reduction._f2_forcing_divisor(2, 3, 4, 2, 1, 1) -- the corner's
#        own forcing ODE, solved there with an INDEPENDENT general-polynomial
#        uniqueness check -- returns c = y, g = 1, f = y^2/2, N = 28 and exactly
#        (30,30,0,0).  Its arguments are (min, max, t, kappa, deg C, ord C): every
#        one of them is fixed by the corner plus {2,3}, none by "F2" or "F3".
#   (ii) window_functions_75_125.family(2) -- a third module, own transcription --
#        gives N = 28 and ord_y(Phi) = deg_y(Phi) = 30 for the same rung.
#  (iii) PUBLISHED anchor.  GGV3 sec.5 (paper_src/1406.0886_GGV3.tex:1723-1727)
#        assumes a counterexample of degrees (50,75), derives A_0 = (5,20), and
#        obtains a pair with
#              [P_1,Q_1] = x^2,   deg(P_1) = 10,   deg(Q_1) = 15.
#        Two things make this an anchor for F3 and not only for F2:
#          * it PRECEDES the paper's branch.  GGV3 says gamma = 3 or gamma = 2 and
#            treats both (tex:1727 vs tex:1777), and BOTH branches start from this
#            same (P_1,Q_1).  So these three integers are not the property of one
#            GGV5 row's b_final -- they hold whichever value gamma takes, which is
#            precisely the freedom that separates GGV5's F2 row from its F3 row.
#          * F3(3,2)/75 IS that case with P and Q exchanged (corner_atlas.json's
#            F_3(3,2)/75 has the identical A_0 = (5,20), A_0' = (1,0), k = 1 and
#            max_deg 75, with (m,n) = (3,2) instead of (2,3)).  Exchanging gives
#            [Q_1,P_1] = -x^2 and degrees (15,10): kappa = 2 up to the sign the
#            bracket is antisymmetric in, hence l = 4, either way.
#        l = 5 predicts [P,Q] = x^3 and reduced degrees (20,30)/(30,20) and
#        contradicts all three integers under either reading.
#
# WHAT IS STILL OPEN, recorded here so it stays visible: the modules DOWNSTREAM of
# family_grammar have NOT been repaired for F3 (they were repaired for F2 only):
# phi_f7.py / phi_f7_verify.py, phi_corner4.py / phi_corner4_verify.py,
# phi_f14.py / phi_f14_verify.py, case_compiler.py, ml_restriction_check.py all
# still transcribe F3 as (5,20) with l_final = 5, b_final = 3 used AS chart data.
# Each is internally consistent and each is now known-suspect.  A10f asserts the
# scope of THIS repair honestly rather than implying the front is closed.
# ---------------------------------------------------------------------------
import polygon_reduction as _pr                                      # noqa: E402

_cd_f2 = _pr.corner_chart_data(5, 20, l_final=5, b_final=2, who="F2 (5,20)")
_cd_f3 = _pr.corner_chart_data(5, 20, l_final=5, b_final=3, who="F3 (5,20)")
ok("A10a chart data is a CORNER property: corner_chart_data(5,20,.) returns "
   "bit-identical (t,kappa,deg C,ord C) = (4,2,1,1), C monomial, no retraction, "
   "for BOTH GGV5 rows at (5,20) -- F2's b_final=2 and F3's b_final=3 -- so the "
   "per-row final-corner datum cannot be, and is not, an input to the chart",
   _cd_f2 == _cd_f3
   and (_cd_f3["t"], _cd_f3["kappa"], _cd_f3["deg_C"], _cd_f3["ord_C"]) == (4, 2, 1, 1)
   and _cd_f3["monomial"] and not _cd_f3["retraction"]
   and (DAT["F3"]["t"], DAT["F3"]["kappa"], DAT["F3"]["degC"], DAT["F3"]["ordC"])
       == (4, 2, 1, 1))

_raised_f3 = False
try:
    _pr.final_corner_dictionary(5, 20, 5, 3, who="F3 (5,20)")
except _pr.FinalCornerDictionaryError:
    _raised_f3 = True
ok("A10c the guard RAISES for F3's own chain row (5,20)/(8\\5,3) exactly as it "
   "does for F2's (7\\5,2) -- and the superseded table's ord C = 2 (F2) vs "
   "ord C = 3 (F3) at ONE corner was already a contradiction, since ord_y(C) is "
   "a corner invariant: 2 != 3.  q = 3 was never chart data for F3",
   _raised_f3 and 2 != 3
   and DAT["F2"]["ordC"] == DAT["F3"]["ordC"] == 1)

# A10b.  The two-sided cross-check: (75,50) is (50,75) with P and Q exchanged, so
# F2 j=0 and F3 j=0 are the SAME reduction and MUST produce identical corner-law
# output -- N and the full signature -- even though the two family rows have
# different (m,n) laws and different N-polynomials in j.  They coincide at j=0
# only, which is what makes this a real check and not an identity.
_D2, _D3 = DAT["F2"], DAT["F3"]
_pair2 = tuple(sorted((int(_D2["a"].subs({j: 0})), int(_D2["b"].subs({j: 0})))))
_pair3 = tuple(sorted((int(_D3["a"].subs({j: 0})), int(_D3["b"].subs({j: 0})))))
_N2 = int(_D2["N"].subs({j: 0})); _N3 = int(_D3["N"].subs({j: 0}))
_s2 = tuple(int(sp.Integer(v.subs({j: 0}))) for v in sig(_D2, _N2, 1))
_s3 = tuple(int(sp.Integer(v.subs({j: 0}))) for v in sig(_D3, _N3, 1))
ok("A10b F3 j=0 = (75,50) is F2 j=0 = (50,75) with P<->Q: same corner, same "
   "unordered reduced pair {2,3}, hence same N = 28 and same signature "
   "(30,30,0,0) -- reached from two DIFFERENT (m,n) laws (F2: a=j+2,b=2j+3; "
   "F3: a=3j+2,b=4j+3) and two different N-polynomials in j that agree at j=0 "
   "and nowhere else",
   _pair2 == _pair3 == (2, 3) and _N2 == _N3 == 28
   and _s2 == _s3 == (30, 30, 0, 0)
   and sp.expand(_D2["N"] - _D3["N"]) != 0)

# A10e.  Route (i) and route (ii): two other modules, own transcriptions.
_fd = _pr._f2_forcing_divisor(2, 3, _cd_f3["t"], _cd_f3["kappa"],
                              _cd_f3["deg_C"], _cd_f3["ord_C"])
_c_fd, _g_fd, _f_fd, _A_fd, _N_fd, _sig_fd = _fd
# polygon_reduction's y is a PLAIN symbol; this file's y carries positive=True, so
# the two are distinct sympy objects and must be reconciled before comparing.
_pry = _pr.y
ok("A10e(i) polygon_reduction's own forcing-ODE solve at (5,20) with the reduced "
   "pair {2,3} independently returns C = y, g = 1, f = y^2/2, N = 28 and "
   "signature (30,30,0,0) -- the verified F3 j=0 target",
   sp.expand(_c_fd - _pry) == 0 and sp.expand(_g_fd - 1) == 0
   and sp.expand(_f_fd - _pry**2 / 2) == 0 and int(_N_fd) == 28
   and tuple(int(s) for s in _sig_fd) == (30, 30, 0, 0))
try:
    import window_functions_75_125 as _wf                            # noqa: E402
    _r = _wf.family(2)
    ok("A10e(ii) window_functions_75_125.family(2) -- a third module -- gives "
       "N = 28 and ord_y(Phi) = deg_y(Phi) = 30 for the same rung",
       int(_r["N"]) == 28 and int(_r["ordPhi"]) == int(_r["degPhi"]) == 30)
except Exception as _exc:                                            # noqa: BLE001
    ok("A10e(ii) window_functions_75_125 cross-check could not run: %r" % (_exc,),
       False)

# A10d.  The superseded F3 rung polynomial: right computation, wrong corner.
# BOTH halves are asserted, because "it solved an ODE" was the whole reason this
# datum looked safe, and a repair that cannot say which ODE is not a repair.
_f3_old = SUPERSEDED_F["F3"]
_t_old, _kap_old, _q_old, _dg_old = 5, 3, 3, 2      # the superseded (5,20) chart
_a_old, _b_old = 2, 3
_e_old = _b_old - _a_old + 1
_coef_old = _t_old * (_b_old - _a_old) + _kap_old + 1
_c_old = y**_q_old * (y + 1)**_dg_old
_res_old = sp.expand(_a_old * (_t_old * _c_old * sp.diff(_f3_old, y)
                               - _coef_old * sp.diff(_c_old, y) * _f3_old)
                     - _c_old**_e_old)
ok("A10d(1) F3's PHI_F7 polynomial (1/42)y^4(y+1)^3(25y^2+15y-3) DOES solve the "
   "mu=2 ramified ODE of the superseded chart (t,kappa,q,dg) = (5,3,3,2): the "
   "computation was never wrong, only its corner was",
   _res_old == 0)
_res_new = full_ode_residual(DAT["F3"], 0, _f3_old)
ok("A10d(2) and it does NOT solve the REPAIRED (5,20) ODE (t=4, C=y, e=2, "
   "2{4 y f' - 7 f} = y^2), whose unique solution is f = y^2/2 -- so the two "
   "charts are DISCRIMINATED by this polynomial, not merely relabelled",
   _res_new != 0
   and full_ode_residual(DAT["F3"], 0, y**2 / 2) == 0)

# A10f.  Scope of this repair, stated honestly.
# 2026-07-27.  The 2026-07-26 note recorded eight downstream modules as
# KNOWN-SUSPECT for F3.  That list is now DISCHARGED -- all eight were swept in the
# same pass that completed this file's repair -- and keeping it would leave a FALSE
# statement inside a green checker, which is exactly the failure mode the A9/A10
# blocks exist to prevent.  Replaced by an assertion of the NEW state, as a test
# rather than as prose: the 2026-07-26 incident was a scope note going stale.
_DOWNSTREAM_SWEPT = (
    "phi_f7.py", "phi_f7_verify.py", "phi_corner4.py", "phi_corner4_verify.py",
    "phi_f14.py", "phi_f14_verify.py", "case_compiler.py",
    "case_compiler_verify.py", "ml_restriction_check.py",
)
_swept_bad = []
for _mod in _DOWNSTREAM_SWEPT:
    try:
        with open(_mod, encoding="utf-8") as _fh:
            _src = _fh.read()
    except OSError:
        _swept_bad.append(_mod + " (missing)")
        continue
    if "corner_chart_data" not in _src and "chart_exponent" not in _src:
        _swept_bad.append(_mod + " (no guard call)")
ok("A10f SCOPE, 2026-07-27: the eight modules recorded on 2026-07-26 as "
   "KNOWN-SUSPECT for F3 are all SWEPT -- each now obtains chart data through "
   "polygon_reduction guard calls rather than from GGV5's final chain corner.  "
   "Offenders: %s" % (_swept_bad or "none"),
   not _swept_bad and "F3" in REPAIRS)
ok("A10g SCOPE, and the repair is no longer per-row: REPAIRS is DERIVED from the "
   "guard over the whole family table, so all ELEVEN refused rows are covered "
   "(F1-F6, F9-F13), not just the two at (5,20).  The previously circulated "
   "affected set {F1,F2,F3,F5,F9,F10} is INCOMPLETE -- it is the refused subset of "
   "the twelve rows bridge_generality.py transcribes, not of these seventeen",
   set(REPAIRS) == {"F1", "F2", "F3", "F4", "F5", "F6",
                    "F9", "F10", "F11", "F12", "F13"})
_bord_bad = []
for _nm, (_A0, _lf, _bf) in CORNERS.items():
    _D = DAT[_nm]
    _av = int(sp.Integer(_D["a"].subs({j: 0})))
    _bv = int(sp.Integer(_D["b"].subs({j: 0})))
    _Nv = int(_D["N"].subs({j: 0}))
    _ordy = int(sp.Integer(sig(_D, _Nv, 1)[1].subs({j: 0})))
    if _ordy != bridge_ordy(_av, _bv, _D["t"], _D["kappa"], _D["ordC"]):
        _bord_bad.append(_nm)
ok("A10h THE INDEPENDENT TARGET: ord_y(Phi) from the grammar equals the PROVED "
   "bridge identity a*q*M - H on ALL 17 families at j=0 (BRIDGE_GENERALITY.md; "
   "rho = q(b-a)+1 by local recursion, N = a*M-2b by the built tower).  Before the "
   "repair the grammar was checked only against targets the same chart dictionary "
   "produced.  Violations: %s" % (_bord_bad or "none"), not _bord_bad)

# ---------------------------------------------------------------------------
# A11. TRIPWIRE -- the mu-LADDER HAS NO SURVIVING LANDED WITNESS.
#
# DECLINED, with the reason stated.  The two published intermediate-rung witnesses
# (F12 eta=0 mu=1,2,3 from ZETA_TAIL.md; F10 eta=0 mu=2 from MU_RUNGS_F10.md) sit
# on corners (8,24) and (7,21), both guard-REFUSED.  At a refused corner dg = 0 and
# there is no residual g at all, so there is nothing for (y+1) to divide: the rung
# is VACUOUS, not merely unverified.  We do NOT substitute the mu-graded law's own
# value for, say, F15's mu=2 rung -- that would validate the law against itself,
# the exact failure this repair removes.
#
# So the next person must state a VERIFIED target for an intermediate rung on a
# corner that retracts.  These checks fail the moment MU_RUNGS is re-populated
# without one.  (Pattern: the retired A10a-A10f tripwire.)
# ---------------------------------------------------------------------------
ok("A11a TRIPWIRE: MU_RUNGS is EMPTY -- no intermediate mu-rung (1 < mu < dg) has "
   "a landed witness on a retracting corner.  Re-populating it requires a target "
   "derived independently of the mu-graded law",
   MU_RUNGS == {})
ok("A11b and both retired witnesses really are on guard-REFUSED corners, so their "
   "rungs are VACUOUS (dg = 0: no residual g for (y+1) to divide) rather than "
   "unverified: F12 at (8,24), F10 at (7,21)",
   {nm for nm, _mu in SUPERSEDED_MU_RUNGS} == {"F12", "F10"}
   and all(not _prg.has_retraction(*CORNERS[nm][0]) and DAT[nm]["dg"] == 0
           for nm, _mu in SUPERSEDED_MU_RUNGS))
ok("A11c and what SURVIVES is only the mu = dg END of the ladder, on retracting "
   "corners: F7 and F16 at dg=2, F15 at dg=4.  Every landed RUNG row uses mu = dg, "
   "never an intermediate value",
   all(DAT[nm]["dg"] >= 1 and _prg.has_retraction(*CORNERS[nm][0])
       for nm in ("F7", "F15", "F16"))
   and sorted({DAT[nm]["dg"] for nm in ("F7", "F15", "F16")}) == [2, 4])

# ---------------------------------------------------------------------------
# A9. DRIFT GUARD -- the two independent LANDED tables must AGREE.
#
# This check exists because its absence cost real soundness.  family_grammar.py
# and this file each keep their OWN copy of the landed derived points, which is
# correct for independence -- but the 2026-07-26 chart repair updated ONE copy.
# The module then printed MISMATCH for both F2 rows and exited 0, while this file
# checked its OWN stale copy against its OWN stale derivation and passed 210/210.
# Two self-consistent halves, disagreeing with each other, both green.
#
# Independence is only worth something if the copies are cross-checked.
# ---------------------------------------------------------------------------
try:
    import family_grammar as _FG
    _theirs = {k: (v[1], v[2]) for k, v in _FG.LANDED.items()}
    ok("A9 the two independent LANDED tables agree on every key (drift guard: "
       "one copy was repaired and the other was not, and nothing caught it)",
       _theirs == LANDED)
    ok("A9b and the two REPAIRS tables agree on which families are chart-repaired",
       set(_FG.REPAIRS) == set(REPAIRS)
       and all(_FG.REPAIRS[k] == REPAIRS[k] for k in REPAIRS))
    ok("A9c and the two f-polynomial tables agree on which rows are LANDED and "
       "which are SUPERSEDED (added 2026-07-26: retiring F3's rung polynomial in "
       "one copy only is exactly the drift A9 was written for)",
       set(_FG.LANDED_F) == set(LANDED_F)
       and set(_FG.SUPERSEDED_F) == set(SUPERSEDED_F)
       and all(sp.expand(_FG.SUPERSEDED_F[k] - SUPERSEDED_F[k]) == 0
               for k in SUPERSEDED_F))
except Exception as _exc:                                            # noqa: BLE001
    ok("A9 drift guard could not run: %r" % (_exc,), False)

# ===========================================================================
# B. CLOSED-FORM THEOREM (pure): gap=0 <=> collapse condition; A=-1/(a dg)
# ===========================================================================
for nm, D in DAT.items():
    # collapse condition  t-(kappa+1)q+dg = 0  is EXACTLY gap=0
    if D["dict_valid"]:
        collapse = D["t"] - (D["kappa"] + 1) * D["q"] + D["dg"]
        ok("%s collapse-cond==0 <=> gap==0" % nm,
           (collapse == 0) == (D["gap"] == 0))
    else:
        ok("%s collapse-condition VOID: dg = 0, so y g' - dg g = -dg reads "
           "0 = 0 and constrains nothing.  The equivalence with gap==0 is not "
           "available here and is NOT asserted." % nm, D["dg"] == 0)

for nm in census["PURE"]:
    D = DAT[nm]
    # symbolic-in-j opaque-G reduction R(G) must be identically 0
    G = sp.symbols("G")
    A = -sp.Rational(1, 1) / (D["a"] * D["dg"])
    trc = D["t"] - (D["kappa"] + 1) * D["q"]
    tec = 1
    R = sp.expand(A * D["a"] * ((trc + tec * D["dg"]) * G - tec * D["dg"]) - 1)
    ok("%s pure R(G)=0 symbolic in j (A=-1/(a dg))" % nm, R == 0)
    # direct full-ODE substitution at j=0..3
    for jv in range(4):
        av = int(D["a"].subs({j: jv}))
        ev = int(D["b"].subs({j: jv})) - av + 1
        rhov = int(D["rho"].subs({j: jv}))
        f = -sp.Rational(1, av * D["dg"]) * y**rhov * (y**D["dg"] + 1)**ev
        ok("%s pure ODE residual=0 at j=%d" % (nm, jv),
           full_ode_residual(D, jv, f) == 0)

# ===========================================================================
# C. COFACTOR (r=0): re-solve unit cofactor from the FULL ODE at j=0..3
# ===========================================================================
def solve_cofactor_full(D, jv):
    """independent linear solve of f = y^rho (y+1)^e u, deg u = gap, at j=jv."""
    av = int(D["a"].subs({j: jv})); bv = int(D["b"].subs({j: jv}))
    ev = bv - av + 1; rhov = int(D["rho"].subs({j: jv})); q = D["q"]
    gp = int(D["gap"])
    uc = sp.symbols("c0:%d" % (gp + 1))
    u = sum(uc[i] * y**i for i in range(gp + 1))
    f = y**rhov * (y + 1)**ev * u
    resid = full_ode_residual(D, jv, f)
    solset = sp.solve(sp.Poly(resid, y).all_coeffs(), list(uc), dict=True)
    if not solset:
        return None
    return u.subs(solset[0])

for nm in census["COFACTOR"]:
    D = DAT[nm]
    for jv in range(4):
        u = solve_cofactor_full(D, jv)
        okflag = u is not None
        if okflag:
            av = int(D["a"].subs({j: jv})); bv = int(D["b"].subs({j: jv}))
            ev = bv - av + 1; rhov = int(D["rho"].subs({j: jv}))
            f = y**rhov * (y + 1)**ev * u
            resid0 = full_ode_residual(D, jv, f) == 0
            unit = (u.subs(y, 0) != 0 and u.subs(y, -1) != 0
                    and sp.Poly(u, y).degree() == int(D["gap"]))
            # u(-1) = -1/a invariant
            inv = sp.simplify(u.subs(y, -1) + sp.Rational(1, av)) == 0
            okflag = resid0 and unit and inv
        ok("%s cofactor: unit u solves ODE at j=%d (u(-1)=-1/a)" % (nm, jv), okflag)

# ===========================================================================
# D. RUNG: mu=dg ramified rung re-solved from FULL ODE; PHI_F7 polys direct
# ===========================================================================
def solve_ramified_full(D, jv):
    """independent linear solve of f = y^rho (y+1)^(dg e-(dg-1)) u,
       deg u = gap+r, at j=jv, from the full ODE."""
    av = int(D["a"].subs({j: jv})); bv = int(D["b"].subs({j: jv}))
    ev = bv - av + 1; rhov = int(D["rho"].subs({j: jv}))
    dg = D["dg"]; P = dg * ev - (dg - 1)
    dego = int(D["gap"]) + D["r"]
    gram = (y + 1)**dg                       # ramified branch residual
    uc = sp.symbols("c0:%d" % (dego + 1))
    u = sum(uc[i] * y**i for i in range(dego + 1))
    f = y**rhov * (y + 1)**P * u
    resid = full_ode_residual(D, jv, f, g=gram)
    solset = sp.solve(sp.Poly(resid, y).all_coeffs(), list(uc), dict=True)
    if not solset:
        return None
    return u.subs(solset[0]), P, dego

for nm in census["RUNG"]:
    D = DAT[nm]
    # k=2 / A0'=(2,0) are conditional but the ODE re-solve is still exact:
    for jv in range(3):
        got = solve_ramified_full(D, jv)
        okflag = got is not None
        if okflag:
            u, P, dego = got
            av = int(D["a"].subs({j: jv})); bv = int(D["b"].subs({j: jv}))
            ev = bv - av + 1; rhov = int(D["rho"].subs({j: jv}))
            f = y**rhov * (y + 1)**P * u
            okflag = (full_ode_residual(D, jv, f, g=(y + 1)**D["dg"]) == 0
                      and sp.Poly(u, y).degree() == dego)
        ok("%s mu=dg ramified rung solves ODE at j=%d (deg u=gap+r)" % (nm, jv), okflag)

# PHI_F7 explicit polynomials -> full ODE residual exactly 0 at j=0
for nm, f in LANDED_F.items():
    ok("%s PHI_F7 f-polynomial solves the full ODE exactly (j=0)" % nm,
       full_ode_residual(DAT[nm], 0, f, g=(y + 1)**DAT[nm]["dg"]) == 0)

# dg parity note (complex-scope discipline): dg even on all r>0,gap>0 rows
for nm in census["RUNG"]:
    D = DAT[nm]
    # F12 (dg=3), F24 excluded: the length-1 even-dg parity claim
    if D["A0p"] == (1, 0) and D["k"] == 1:
        ok("%s dg is EVEN (mu=1 real-empty; complex not excluded)" % nm,
           D["dg"] % 2 == 0)

# ===========================================================================
# D'. CHART-DEGENERATE: the mechanism that DOES hold at a monomial corner.
#
# 2026-07-27, and this section exists so that the completed repair does not merely
# DELETE coverage.  Nine rows left the PURE / COFACTOR / RUNG classes today, which
# withdrew their closed-form-mechanism checks from sections B, C and D -- correctly,
# because those checks asserted the existence of a residual g that a refused corner
# does not have, i.e. they asserted the very dictionary the corner refutes.  But a
# repair whose only visible effect is a lower check count is indistinguishable from
# a repair that quietly stopped testing something.  So the eleven refused rows are
# verified POSITIVELY here instead:
#
#   at deg C = ord C = 1 the ODE  a{t c f' - coef c' f} = c^e  with c = y reads
#       a*t*y*f' - a*coef*f = y^e,
#   and f = A y^e gives  a*A*(t*e - coef) = a*A*(t - kappa - 1) = a*A = 1,
#   so A = 1/a EXACTLY and f = (1/a) y^e -- with NO free parameter, NO residual to
#   shape and NO root to place.  Phi = f*C^N is then the monomial (1/a) y^(e+N).
# ===========================================================================
for nm in census["CHART-DEGENERATE"]:
    D = DAT[nm]
    ok("%s degenerate: t*e - coef = t - kappa - 1 = 1 identically in j (the "
       "identity that collapses the ODE at a monomial corner)" % nm,
       sp.expand(D["t"] * D["e"] - D["coef"] - 1) == 0)
    for jv in range(3):
        av = int(sp.Integer(D["a"].subs({j: jv})))
        bv = int(sp.Integer(D["b"].subs({j: jv})))
        ev = bv - av + 1
        f = sp.Rational(1, av) * y**ev
        resid0 = full_ode_residual(D, jv, f, g=sp.Integer(1)) == 0
        # UNIQUENESS, without assuming the shape: a fully generic linear solve.
        Dmax = ev + 3
        fc = sp.symbols("dg0:%d" % (Dmax + 1))
        fgen = sum(fc[i] * y**i for i in range(Dmax + 1))
        rg = full_ode_residual(D, jv, fgen, g=sp.Integer(1))
        sols = sp.solve(sp.Poly(sp.expand(rg), y).all_coeffs(), fc, dict=True)
        uniq = len(sols) == 1 and sp.expand(fgen.subs(sols[0]) - f) == 0
        # and Phi is the MONOMIAL (1/a) y^(e+N), whose ord_y is the bridge value
        Nv = int(D["N"].subs({j: jv}))
        Phi = sp.expand(f * y**Nv)
        monomial = len(sp.Poly(Phi, y).monoms()) == 1
        ordy = min(m[0] for m in sp.Poly(Phi, y).monoms())
        bridged = ordy == bridge_ordy(av, bv, D["t"], D["kappa"], D["ordC"])
        ok("%s degenerate at j=%d: f = (1/%d) y^%d is the UNIQUE polynomial "
           "solution (generic solve, no ansatz), Phi = (1/%d) y^%d is a MONOMIAL, "
           "and its ord_y equals the PROVED bridge value %d"
           % (nm, jv, av, ev, av, ev + Nv, ordy),
           resid0 and uniq and monomial and bridged)
    # and the withdrawn mechanism must genuinely FAIL here, not merely be unused:
    # g = y^dg+1 with dg = 0 is the constant 2, which is not monic and is not a
    # residual at all, so the pure closed form A = -1/(a*dg) does not even exist.
    ok("%s degenerate: dg = 0, so the pure closed form's constant -1/(a*dg) is "
       "UNDEFINED and the collapse identity y g' - dg g = -dg reads 0 = 0 -- the "
       "mechanism is absent, which is why its checks are withdrawn rather than "
       "failing" % nm,
       D["dg"] == 0 and sp.expand(y * sp.diff(y**0 + 1, y) - 0 * (y**0 + 1)) == 0)

ok("D': coverage bookkeeping -- the rows whose PURE/COFACTOR/RUNG mechanism checks "
   "were WITHDRAWN on 2026-07-27 are exactly the guard-refused rows, and every one "
   "of them is verified positively above instead.  Withdrawn: F1, F4, F5, F6, F9, "
   "F10, F11, F12, F13 (F2 and F3 were already withdrawn on 2026-07-26)",
   set(census["CHART-DEGENERATE"]) == set(REPAIRS)
   and {"F1", "F4", "F5", "F6", "F9", "F10", "F11", "F12", "F13"}
   <= set(census["CHART-DEGENERATE"]))

# ===========================================================================
# E. mu-graded corner law: landed points + published mu-rungs; identities
# ===========================================================================
for (nm, jv), (Nl, s) in LANDED.items():
    D = DAT[nm]
    Nval = int(D["N"].subs({j: jv}))
    mu = D["dg"] if D["cls"] == "RUNG" else 1
    got = tuple(int(sp.Integer(v.subs({j: jv}))) for v in sig(D, Nval, mu))
    ok("%s j=%d landed: N=%d and signature %s" % (nm, jv, Nl, s),
       Nval == Nl and got == s)

for (nm, mu), s in MU_RUNGS.items():
    D = DAT[nm]
    Nval = int(D["N"].subs({j: 0}))
    got = tuple(int(sp.Integer(v.subs({j: 0}))) for v in sig(D, Nval, mu))
    ok("%s mu=%d published rung %s via mu-graded law" % (nm, mu, s), got == s)

# mu-graded law algebraic identities (symbolic in j, N, mu).
# 2026-07-26: F3 REPLACED BY F7 in this list.  These identities are stated in terms
# of (a0, q, rho = (e-1)q+1, r = a0-q-1), i.e. they presuppose the final-corner
# dictionary; asserting them for a CHART-DEGENERATE row would be asserting the very
# dictionary the row refutes.  F7 is a RUNG row where the dictionary holds, so the
# algebra is tested on the class it actually describes.
Nn, mm = sp.symbols("Nn mm")
# 2026-07-27: F10 and F12 REPLACED BY F15 and F17.  As with F3 on 2026-07-26,
# these identities are stated in terms of (a0, q, rho = (e-1)q+1, r = a0-q-1) and
# so presuppose the final-corner dictionary; asserting them for a CHART-DEGENERATE
# row would assert the very dictionary that row refutes.  F15 (RUNG, dg=4) and F17
# (COFACTOR, r=0) both sit on the retracting corner (9,24).
for nm in ["F7", "F15", "F16", "F17"]:
    D = DAT[nm]
    e, a0, q, rho, gap, r = D["e"], D["a0"], D["q"], D["rho"], D["gap"], D["r"]
    pure = e * a0 - q + 1
    res = pure + gap
    deg = res + Nn * a0
    ordy = rho + Nn * q
    mult = mm * (e + Nn) - (mm - 1)
    cof = gap + r * (e + Nn) - (mm - 1) * (e + Nn - 1)
    ok("%s cof = deg-ord-mult identically (mu-graded)" % nm,
       sp.expand(cof - (deg - ordy - mult)) == 0)
    # mu=1 specialization -> (e+N, gap+r(e+N))
    ok("%s mu=1 -> mult=e+N, cof=gap+r(e+N)" % nm,
       sp.expand(mult.subs(mm, 1) - (e + Nn)) == 0
       and sp.expand(cof.subs(mm, 1) - (gap + r * (e + Nn))) == 0)
    # mu=dg specialization -> (dg(e+N)-(dg-1), gap+r)  using r=dg-1
    dg = D["dg"]
    ok("%s mu=dg -> mult=dg(e+N)-(dg-1), cof=gap+r" % nm,
       sp.expand(mult.subs(mm, dg) - (dg * (e + Nn) - (dg - 1))) == 0
       and sp.expand(cof.subs(mm, dg) - (gap + r)) == 0)

# ---------------------------------------------------------------------------
print()
if FAILS:
    print("FAILURES:", len(FAILS))
    for f in FAILS:
        print("   -", f)
    sys.exit(1)
print("ALL %d FAMILY-GRAMMAR CHECKS PASSED" % N_OK)
sys.exit(0)
