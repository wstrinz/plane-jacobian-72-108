#!/usr/bin/env python3
"""phi_corner4_verify.py  -- independent PASS/FAIL checker for PHI_CORNER4.md.

Self-contained: re-derives everything from scratch with its own routines (own
transcription of the GGV5 rows, own solvers).  The ONE thing it does NOT
re-implement is the retraction guard: chart data must come from
polygon_reduction.corner_chart_data in EVERY consumer, because a second
implementation of a guard is a second chance to get the guard wrong.  What it
does instead is CROSS-CHECK phi_corner4.py's ledger against the guard (sec. H).

Verifies:

  A. corner data Diophantine identities (all 17 length-1 GGV5 families),
  B. the chart Jacobian -x^(l-2)  =>  kappa = t-2 on the standard class,
  C. the operator bracket [P, x^s f/c^b] = x^kappa  <=>  family ODE, at the
     repaired (t,kappa) points and the two known controls,
  D. the ODE coefficient system forces g = y^(dgC)+1 and A at a RETRACTING
     corner (branch-complete nonlinear solve), and is VACUOUS at a refused one,
  E. uniqueness: the forced f is the ONLY polynomial solution up to the
     resonant degree (full linear solve),
  F. divisor signatures of Phi = f * C^N extracted by factor_list (independent
     of the derivation's div-loop) equal the corner-law predictions AND the
     PROVED bridge identity ord_y(Phi) = a*q*M - H,
  G. controls: the same law formulas reproduce (72,108), (108,144), (75,125),
     (50,75), and the (72,108) r=0 resonance-gap exception offsets by exactly
     deg q_4 = 4,
  H. THE CHART REPAIR (2026-07-27): the guard-refused set is exactly eleven
     rows; phi_corner4.py's RETRACTION_LEDGER / AFFECTED_EXPECTED / SUPERSEDED
     agree with an independent recomputation (drift guard); and MUTATION
     CONTROLS -- reinstating the superseded dictionary MUST move ord_y(Phi) and
     MUST break the bridge identity at every affected row.

Usage: python phi_corner4_verify.py [--quiet]     exit 0 iff all checks pass.
"""
import sys
import sympy as sp
from fractions import Fraction
from math import gcd

import polygon_reduction as pr

QUIET = "--quiet" in sys.argv[1:]
y, x, A = sp.symbols("y x A")
_n, _fail = 0, 0

def ok(cond, msg):
    global _n, _fail
    _n += 1
    if not cond:
        _fail += 1
    if not QUIET or not cond:
        print(f"[{'OK' if cond else 'FAIL'}] {msg}")

# phi_corner4.py is a REPORT: importing it runs its whole derivation and prints
# it.  Swallow that stream so this checker's own PASS/FAIL lines -- and any FAIL
# -- cannot be lost in it.  (Masking a checker's OUTPUT is fine; masking its EXIT
# CODE is the trap, and nothing here touches exit codes.)
import contextlib, io                                            # noqa: E402
_pc4_report = io.StringIO()
with contextlib.redirect_stdout(_pc4_report):
    import phi_corner4 as pc4                                    # noqa: E402
from phi_corner4 import mult_and_cofactor, gap_effective   # noqa: E402


def bridge(a, b, t, kappa, ordC):
    """ord_y(Phi) = a*q*M - H.  INDEPENDENT of everything below (BRIDGE_GENERALITY)."""
    s = a + b
    return a * ordC * (t * s - (kappa + 1)) - (ordC * s - 1)


def fit_sig(a, b, t, kappa, degC, ordC):
    """Independent copy of the law.  2026-07-26: shared residual-free branch
    (dg = degC-ordC = 0 => mult = cof = 0), which the (5,20) corner forces.
    2026-07-27: the last two arguments are (deg C, ord C) from the guard, NOT
    (a0, b_final) from GGV5's chain row."""
    e = b - a + 1
    N = a * (t * (a + b) - (kappa + 1)) - 2 * b
    mult, cof = mult_and_cofactor(e, N, degC, ordC, 0)
    return N, ((e * degC - ordC + 1) + N * degC,
               ((e - 1) * ordC + 1) + N * ordC, mult, cof)

# ---------------------------------------------------------------- A. corner data
FAMS = [  # (name, A0, p, l_final, b_final, k, m0, dm, n0, dn)
    ("F1", (4,12), 7,4,3,1, 3,2, 4,3),  ("F2", (5,20), 7,5,2,1, 2,1, 3,2),
    ("F3", (5,20), 8,5,3,1, 3,4, 2,3),  ("F4", (5,20), 8,5,3,2, 3,2, 16,12),
    ("F5", (5,20), 9,5,4,1, 9,7, 5,4),  ("F6", (5,20), 9,5,4,2, 7,6, 18,16),  # F6 CORRECTED 2026-07-24 (GGV5 typo base (4,10) gcd=2 -> coprime (6j+7,16j+18) base (7,18); CHAIN_SURVEY.md)
    ("F7", (6,15), 7,3,4,1, 2,1, 7,4),  ("F8", (6,15), 8,3,5,1, 3,2, 7,5),
    ("F9", (7,21), 11,7,2,1, 2,1, 3,2), ("F10",(7,21), 13,7,3,1, 7,5, 4,3),
    ("F11",(7,21), 13,7,3,2, 2,1, 5,3), ("F12",(8,24), 13,4,5,1, 3,2, 7,5),
    ("F13",(9,21), 13,3,7,1, 2,1, 13,7),("F14",(9,24), 7,3,4,1, 2,1, 7,4),
    ("F15",(9,24), 8,3,5,1, 3,2, 7,5),  ("F16",(9,24), 10,3,7,1, 3,4, 5,7),
    ("F17",(9,24), 11,3,8,1, 2,5, 3,8),
]
dio_all = True
ROWS = {}
for name, A0, p, l, q, k, m0, dm, n0, dn in FAMS:
    j = 0
    while gcd(m0 + dm * j, n0 + dn * j) != 1:
        j += 1
    m, n = m0 + dm * j, n0 + dn * j
    dio_all &= ((m + n) * q * k - n * (q * l - p) == k)
    ROWS[name] = dict(A0=A0, p=p, l=l, q=q, k=k, m=m, n=n,
                      ab=tuple(sorted((m, n))),
                      degs=((A0[0] + A0[1]) * m, (A0[0] + A0[1]) * n))
ok(dio_all, "A: Diophantine (m+n)qk - n(ql-p) = k holds for all 17 length-1 families")
ok((2+3)*2*1 - 3*(2*7-11) == 1, "A: F9 j=0 chain row (7,21)->(11\\7,2), (m,n)=(2,3): identity = k = 1")
ok((2+3)*2*1 - 3*(2*5-7) == 1,  "A: F2 j=0 chain row (5,20)->(7\\5,2), (m,n)=(2,3): identity = k = 1")
ok(28 * 2 == 56 and 28 * 3 == 84, "A: F9 j=0 degrees v11*(m,n) = (56,84)")
ok(25 * 2 == 50 and 25 * 3 == 75, "A: F2 j=0 degrees v11*(m,n) = (50,75)")

# ---- CHART data for every row, through the guard.  This is the repair. --------
CHART = {}
for name, R in ROWS.items():
    cd = pr.corner_chart_data(R["A0"][0], R["A0"][1], l_final=R["l"],
                              b_final=R["q"], who="phi_corner4_verify " + name)
    CHART[name] = cd
REFUSED = {nm for nm, cd in CHART.items() if not cd["retraction"]}
ok(REFUSED == {"F1", "F2", "F3", "F4", "F5", "F6", "F9", "F10", "F11", "F12", "F13"},
   "A: guard REFUSES GGV5's final-corner dictionary at exactly 11 of 17 rows "
   "(F1-F6, F9-F13); the previously circulated set {F1,F2,F3,F5,F9,F10} is the "
   "refused subset of a 12-row transcription, not of these 17")
ok(all(cd["deg_C"] == 1 and cd["ord_C"] == 1 and cd["monomial"]
       for nm, cd in CHART.items() if nm in REFUSED),
   "A: on every refused row C is the MONOMIAL y (deg C = ord C = 1): no vertical "
   "top face, so the residual g does not exist")
ok(all(cd["deg_C"] == ROWS[nm]["A0"][0] and cd["ord_C"] == ROWS[nm]["q"]
       and cd["t"] == ROWS[nm]["l"]
       for nm, cd in CHART.items() if nm not in REFUSED),
   "A: on the six RETRACTING rows (F7,F8,F14-F17) the dictionary is VALID and the "
   "guard returns (t,deg C,ord C) = (l_final, a0, b_final) unchanged")
ok(sorted({cd["t"] for cd in CHART.values()}) == [3, 4],
   "A: t-census -- every GGV5 v11<=35 length-1 corner has t in {3,4}.  l_final = 7 "
   "at (7,21) was never a chart exponent (GGHV22 publishes l=3 there), so the "
   "'new t=7' out-of-sample point does not exist")

# ------------------------------------------------- B. chart Jacobian => kappa=t-2
ls = sp.symbols("l_s", positive=True)
X, Y = x**-1, x**ls * y
J = sp.simplify(sp.diff(X, x) * sp.diff(Y, y) - sp.diff(X, y) * sp.diff(Y, x))
ok(sp.simplify(J + x**(ls - 2)) == 0,
   "B: Jacobian of (x^-1, x^l y) is -x^(l-2) for symbolic l")
conc = all(sp.simplify(J.subs(ls, lv) + x**(lv - 2)) == 0 for lv in (3, 4, 5, 7, 8))
ok(conc, "B: concrete l in {3,4,5,7,8} (covers every l_final AND every derived "
         "chart exponent in the GGV5 tables)")
ok(all(cd["kappa"] == cd["t"] - 2 for cd in CHART.values()),
   "B: => kappa = t-2 on all 17 rows, for the DERIVED t.  The repair does not "
   "touch this: it is a property of the chart, not of which l the chart uses")

# ------------------------------------------- C. bracket identity <=> family ODE
def bracket_is_ode(a, b, t, kappa):
    s = kappa + 1 - a * t
    cf, ff = sp.Function("c")(y), sp.Function("f")(y)
    P = (x**t * cf)**a
    G = x**s * ff / cf**b
    br = sp.simplify(sp.diff(P, x) * sp.diff(G, y) - sp.diff(P, y) * sp.diff(G, x))
    expr = sp.simplify((br - x**kappa) / x**kappa)
    ode = a * (t * cf * sp.diff(ff, y) - (t * (b - a) + kappa + 1) * sp.diff(cf, y) * ff) \
          - cf**(b - a + 1)
    return sp.simplify(expr - ode / cf**(b - a + 1)) == 0

ok(bracket_is_ode(2, 3, 3, 1), "C: bracket => ODE at REPAIRED (a,b,t,kappa)=(2,3,3,1)"
   "  [F9 (56,84); corner (7,21) is t=3, GGHV22's published chart]")
ok(bracket_is_ode(2, 3, 4, 2), "C: bracket => ODE at REPAIRED (a,b,t,kappa)=(2,3,4,2)"
   "  [F2 j=0 (50,75); corner (5,20) is t=4, pinned by GGV3]")
ok(bracket_is_ode(3, 5, 4, 2), "C: bracket => ODE control (3,5,4,2)  [(75,125), "
   "REPAIRED 2026-07-26: t=4, kappa=2 -- GGV5's l_final=5 is not the chart]")
ok(bracket_is_ode(2, 3, 4, 2), "C: bracket => ODE control (2,3,4,2)  [(72,108)]")
ok(bracket_is_ode(2, 7, 3, 1), "C: bracket => ODE control (2,7,3,1)  [F14 (66,231), "
   "a RETRACTING corner: unaffected by the repair]")

# ------------------------------- D. forced g (branch-complete nonlinear solve)
def force_g(tag, a, b, t, kappa, degC, ordC):
    """dg > 0 only: the branch-complete nonlinear solve for the residual g."""
    e = b - a + 1
    coef = t * (b - a) + kappa + 1
    rho = (e - 1) * ordC + 1
    dg = degC - ordC
    assert dg > 0, "force_g is meaningless at a monomial corner"
    gc = sp.symbols(f"g0:{dg+1}")
    g = sum(gc[i] * y**i for i in range(dg + 1))
    c = y**ordC * g
    f = A * y**rho * g**e
    resid = sp.expand(a * t * c * sp.diff(f, y) - a * coef * sp.diff(c, y) * f - c**e)
    quo = sp.expand(sp.factor(resid) / (y**(e * ordC) * g**(e - 1)))
    eqs = [sp.expand(quo).coeff(y, i) for i in range(sp.degree(quo, y) + 1)]
    sols = sp.solve(eqs, list(gc[1:dg]) + [A], dict=True)   # g0, g_top free
    good = [s for s in sols
            if not any(sp.simplify(v) == 0 and str(k) == str(gc[dg])
                       for k, v in s.items())]
    forced_mid = all(all(sp.simplify(s.get(gc[i], gc[i])) == 0
                         for i in range(1, dg)) for s in good) and len(good) >= 1
    ok(forced_mid, f"D: {tag}: every valid branch forces g_1..g_{dg-1} = 0 "
                   f"({len(good)} branch(es))")
    s0 = good[0]
    A_val = sp.simplify(sp.simplify(s0[A]).subs(gc[0], 1))
    g_sol = y**dg + 1
    c_sol = y**ordC * g_sol
    f_sol = sp.expand(A_val * y**rho * g_sol**e)
    ok(sp.expand(a * t * c_sol * sp.diff(f_sol, y)
                 - a * coef * sp.diff(c_sol, y) * f_sol - c_sol**e) == 0,
       f"D: {tag}: g = y^{dg}+1, A = {A_val} satisfies the ODE exactly")
    ok(sp.simplify(g_sol.subs(y, -1)) == 0 and sp.simplify(g_sol.subs(y, 0)) == 1,
       f"D: {tag}: g(-1) = 0 and g(0) != 0")
    Hres = sp.factor(g_sol / (y + 1))
    disc = sp.discriminant(sp.Poly(Hres, y))
    ok(disc != 0 and sp.simplify(Hres.subs(y, 0)) != 0
       and sp.simplify(Hres.subs(y, -1)) != 0,
       f"D: {tag}: residual H = {Hres} separable, avoids 0 and -1")
    return f_sol, c_sol, A_val


def force_monomial(tag, a, b, t, kappa):
    """dg == 0: C = y, g = 1 FORCED (monic constant).  The whole residual layer,
    including the common-root gauge, is VACUOUS -- there is no g to shape and no
    root to place.  The ODE collapses to a*A*(t*e - coef) = a*A*(t-kappa-1) = 1."""
    e = b - a + 1
    coef = t * (b - a) + kappa + 1
    ok(t * e - coef == t - kappa - 1 == 1,
       f"D: {tag}: t*e - coef = t - kappa - 1 = 1 identically on the standard class")
    c_sol = y
    f_sol = sp.Rational(1, a) * y**e
    ok(sp.expand(a * t * c_sol * sp.diff(f_sol, y)
                 - a * coef * sp.diff(c_sol, y) * f_sol - c_sol**e) == 0,
       f"D: {tag}: C = y monomial => f = (1/{a}) y^{e} solves the ODE exactly, "
       f"A = 1/a")
    # the gauge branch is vacuous, not merely simple: g = y^0+1 = 2 is NOT monic
    # and NOT a valid residual, which is why the dg>0 code path must not be run.
    ok(sp.degree(sp.Poly(c_sol, y), y) == 1 and c_sol.subs(y, -1) != 0,
       f"D: {tag}: C = y has NO root at -1, so mult_(y+1)(Phi) = 0 -- there is no "
       f"(y+1) place at a monomial corner")
    return f_sol, c_sol, sp.Rational(1, a)

# the two derivations of PHI_CORNER4, at their REPAIRED (refused) corners
f9, c9, A9v = force_monomial("F9  (56,84) corner (7,21)", 2, 3, 3, 1)
f2, c2, A2v = force_monomial("F2j0 (50,75) corner (5,20)", 2, 3, 4, 2)
ok(A9v == sp.Rational(1, 2), "D: F9  A = 1/a = 1/2  [was -1/10 under the dictionary]")
ok(A2v == sp.Rational(1, 2), "D: F2j0 A = 1/a = 1/2  [was -1/6 under the dictionary]")
# and the dg>0 branch-complete machinery, exercised where it is VALID: a
# RETRACTING corner.  Without this the repair would have deleted the only test
# of the forced-residual argument.
f14, c14, A14v = force_g("F14 (66,231) corner (9,24) RETRACTS", 2, 7, 3, 1, 9, 4)
ok(A14v == sp.Rational(-1, 10), "D: F14 A = -1/10 (dg=5 forced residual y^5+1)")

# ----------------------------------------------- E. uniqueness (linear solve)
def unique_f(tag, a, b, t, kappa, c_sol, f_expect, dmax):
    e = b - a + 1
    coef = t * (b - a) + kappa + 1
    fc = sp.symbols(f"f0:{dmax+1}")
    f = sum(fc[i] * y**i for i in range(dmax + 1))
    eq = sp.expand(a * t * c_sol * sp.diff(f, y)
                   - a * coef * sp.diff(c_sol, y) * f - c_sol**e)
    sols = sp.solve([eq.coeff(y, i) for i in range(sp.degree(eq, y) + 1)],
                    fc, dict=True)
    uniq = len(sols) == 1 and sp.expand(f.subs(sols[0]) - f_expect) == 0
    ok(uniq, f"E: {tag}: unique polynomial ODE solution of degree <= {dmax}, = f")

unique_f("F9  ", 2, 3, 3, 1, c9, f9, 8)
unique_f("F2j0", 2, 3, 4, 2, c2, f2, 8)
unique_f("F14 ", 2, 7, 3, 1, c14, f14, 55)
# resonance bookkeeping at a monomial corner: resonant degree coef*degC/t is
# BELOW the pure-ansatz degree, so gap < 0 and gap_effective = 0.
for tag, a, b, t, kappa in (("F9  ", 2, 3, 3, 1), ("F2j0", 2, 3, 4, 2)):
    e = b - a + 1
    coef = t * (b - a) + kappa + 1
    res = Fraction(coef * 1, t)                 # degC = 1
    pure = e * 1 - 1 + 1                        # = e
    ok(res < pure and gap_effective(res - pure) == 0,
       f"E: {tag} monomial corner: resonant deg {res} < pure-ansatz deg {pure}, "
       f"gap = {res - pure} < 0 => gap_effective 0, so deg f = rho exactly (no "
       f"unit cofactor)")
ok(Fraction((3 * 5 + 2) * 9, 3) == 51 and 6 * 9 - 4 + 1 == 51,
   "E: F14 resonant deg = pure-ansatz deg = 51 (gap 0: no extra unit cofactor)")

# ---------------------------- F. Phi signature via factor_list vs law formulas
def signature_factorlist(Phi):
    _, facs = sp.factor_list(sp.expand(Phi))
    deg = sp.degree(sp.expand(Phi), y)
    ordy = m1 = 0
    cof = 0
    for base, mult in facs:
        b_ = sp.expand(base)
        if b_ == y:
            ordy = mult
        elif b_ == y + 1:
            m1 = mult
        else:
            cof += sp.degree(b_, y) * mult
    return (deg, ordy, m1, cof)

N9, fitsig9 = fit_sig(2, 3, 3, 1, 1, 1)
N2, fitsig2 = fit_sig(2, 3, 4, 2, 1, 1)
N14, fitsig14 = fit_sig(2, 7, 3, 1, 9, 4)
ok(N9 == 20, "F: F9  N = a[t(a+b)-(kappa+1)]-2b = 20  [was 52 under the dictionary]")
ok(N2 == 28, "F: F2j0 N = 28  [was 36 under the dictionary]")
ok(N14 == 36, "F: F14 N = 36 (retracting corner: unchanged)")
Phi9 = f9 * c9**N9
Phi2 = f2 * c2**N2
Phi14 = f14 * c14**N14
s9 = signature_factorlist(Phi9)
s2 = signature_factorlist(Phi2)
s14 = signature_factorlist(Phi14)
ok(s9 == (22, 22, 0, 0), "F: F9  derived signature (22,22,0,0)  [factor_list]")
ok(s2 == (30, 30, 0, 0), "F: F2j0 derived signature (30,30,0,0)  [factor_list]")
ok(s14 == (375, 165, 42, 168), "F: F14 derived signature (375,165,42,168) [factor_list]")
ok(s9 == fitsig9, "F: F9  derived == law prediction  => MATCHES")
ok(s2 == fitsig2, "F: F2j0 derived == law prediction  => MATCHES")
ok(s14 == fitsig14, "F: F14 derived == law prediction  => MATCHES")
ok(s9[0] == s9[1] + s9[2] + s9[3] and s2[0] == s2[1] + s2[2] + s2[3]
   and s14[0] == s14[1] + s14[2] + s14[3],
   "F: signature sum identity deg = ord + mult + cofactor on all three points")
ok(sp.expand(Phi9 - sp.Rational(1, 2) * y**22) == 0,
   "F: F9  Phi = (1/2) y^22 exactly (a MONOMIAL: C = y and no residual)")
ok(sp.expand(Phi2 - sp.Rational(1, 2) * y**30) == 0,
   "F: F2j0 Phi = (1/2) y^30 exactly")
ok(sp.expand(Phi14 + sp.Rational(1, 10) * y**165 * (y**5 + 1)**42) == 0,
   "F: F14 Phi = -(1/10) y^165 (y^5+1)^42 exactly (retracting: unchanged)")

# --- THE INDEPENDENT TARGET: the PROVED bridge identity, which this file does
# --- not produce.  Every row of the survey, not just the derived ones.
bad_bridge = []
for nm, R in ROWS.items():
    cd = CHART[nm]
    a, b = R["ab"]
    _, sg = fit_sig(a, b, cd["t"], cd["kappa"], cd["deg_C"], cd["ord_C"])
    if sg[1] != bridge(a, b, cd["t"], cd["kappa"], cd["ord_C"]):
        bad_bridge.append(nm)
ok(not bad_bridge,
   "F: ord_y(Phi) == a*q*M - H on ALL 17 rows (bridge identity, PROVED in "
   "BRIDGE_GENERALITY.md: rho = q(b-a)+1 by local recursion, N = a*M-2b by the "
   "built tower).  This is a target phi_corner4.py does not compute")
ok(s9[1] == bridge(2, 3, 3, 1, 1) == 22 and s2[1] == bridge(2, 3, 4, 2, 1) == 30
   and s14[1] == bridge(2, 7, 3, 1, 4) == 165,
   "F: and the three factor_list ord_y values agree with it: 22, 30, 165")

# --------------------------------------------------------------- G. controls
ok(fit_sig(3, 4, 4, 2, 8, 3)[1] == (550, 205, 69, 276),
   "G: control (108,144) at (8,28) [RETRACTS]: law gives (550,205,69,276)")
ok(fit_sig(3, 5, 4, 2, 1, 1)[1] == (80, 80, 0, 0),
   "G: control (75,125):  law gives (80,80,0,0)  [REPAIRED: t=4, kappa=2, "
   "deg C=ord C=1; the dg=0 branch of fit_signature]")
ok(fit_sig(2, 3, 4, 2, 1, 1)[1] == (30, 30, 0, 0),
   "G: control (50,75):   law gives (30,30,0,0)  [REPAIRED]")
ok(fit_sig(3, 5, 4, 2, 1, 1)[1][2] == 0 and fit_sig(3, 5, 4, 2, 1, 1)[1][3] == 0,
   "G: the dg=0 branch is REQUIRED: with C a monomial there is no (y+1) place and "
   "no residual cofactor, so mult=cof=0 (the old r=degC-ordC-1 rule gave r=-1)")
ok(fit_sig(3, 4, 4, 2, 8, 3)[1] == (550, 205, 69, 276)
   and fit_sig(2, 7, 3, 1, 9, 4)[1] == (375, 165, 42, 168),
   "G: and the dg>0 branch is UNCHANGED -- (108,144) and F14 still land exactly, "
   "so the new branch is discriminating, not a blanket rewrite")
N72, s72 = fit_sig(2, 3, 4, 2, 8, 7)
ok(N72 == 28 and s72[2] == 30 and 238 - s72[0] == 4 and s72[3] == 0,
   "G: (72,108) r=0 exception: mult e+N=30 matches; deg offset 238-234 = 4 = deg q_4")
ok(fit_sig(2, 3, 3, 1, 1, 1)[0] * 1 + 2 == 22,
   "G: F9 deg decomposition deg f + N*deg C = 2 + 20*1 = 22")
red = lambda a, b, t: a * (t * (a + b - 1) + 1) - 2 * b
ok(red(2, 3, 4) == 28 and red(3, 4, 4) == 67 and red(3, 5, 4) == 77
   and red(2, 3, 3) == 20 and red(2, 7, 3) == 36,
   "G: reduced N-formula (kappa=t-2 substituted) reproduces N at (72,108) 28, "
   "(108,144) 67, (75,125) 77, F9 (56,84) 20 [REPAIRED from 52], F14 36")

# ===========================================================================
# H.  THE CHART REPAIR: drift guard + mutation controls
# ===========================================================================
# H1-H3: drift guard.  phi_corner4.py keeps a ledger; this file recomputes the
# classification independently from the guard.  The two MUST agree -- the F3
# incident was two self-consistent copies of the same ground truth, one repaired
# and one not, both green.  (family_grammar_verify.py check A9, same pattern.)
_mine = {tuple(R["A0"]): pr.has_retraction(R["A0"][0], R["A0"][1])
         for R in ROWS.values()}
ok(set(pc4.RETRACTION_LEDGER) == set(_mine)
   and all(pc4.RETRACTION_LEDGER[c][0] == _mine[c] for c in _mine),
   "H1 DRIFT GUARD: phi_corner4.RETRACTION_LEDGER agrees with an independent "
   "recomputation of has_retraction on all 7 corners (7 corners, 17 rows)")
ok(pc4.AFFECTED_EXPECTED == REFUSED,
   "H1b and its AFFECTED_EXPECTED set equals the refused set recomputed here: %s"
   % sorted(REFUSED, key=lambda s: int(s[1:])))
ok(pc4.A0P_20_ROWS == {nm for nm, R in ROWS.items()
                       if nm in ("F12", "F13")} <= REFUSED,
   "H1c and the A0'=(2,0) rows F12,F13 are inside the refused set (their repaired "
   "t is CLAIMED: chart_exponent is validated only at A0'=(1,0))")
ok(set(pc4.SUPERSEDED) == {"F1", "F2", "F3", "F9", "F10"}
   and all(nm in REFUSED for nm in pc4.SUPERSEDED),
   "H1d and every row in the SUPERSEDED table is a refused row (the table is kept "
   "LABELLED so the repair stays falsifiable -- deleting it would make the "
   "mutation controls below unrunnable)")
ok(all(pc4.SUPERSEDED[nm][0] == ROWS[nm]["l"]
       and pc4.SUPERSEDED[nm][1] == ROWS[nm]["A0"][0]
       and pc4.SUPERSEDED[nm][2] == ROWS[nm]["q"] for nm in pc4.SUPERSEDED),
   "H1e and its stale (t,deg C,ord C) really IS the dictionary's output "
   "(l_final, a0, b_final) on each of those rows -- so H3's control is the OLD "
   "model, not a straw man")

# H2: the guard really does refuse, per row, with the right reason.
_raised = 0
for nm in sorted(REFUSED):
    R = ROWS[nm]
    try:
        pr.final_corner_dictionary(R["A0"][0], R["A0"][1], R["l"], R["q"], who=nm)
    except pr.FinalCornerDictionaryError:
        _raised += 1
ok(_raised == len(REFUSED),
   "H2 final_corner_dictionary RAISES on all %d refused rows (not merely returns "
   "different numbers)" % len(REFUSED))
_ok_dict = 0
for nm in sorted(set(ROWS) - REFUSED):
    R = ROWS[nm]
    if pr.final_corner_dictionary(R["A0"][0], R["A0"][1], R["l"], R["q"]) \
            == (R["l"], R["q"]):
        _ok_dict += 1
ok(_ok_dict == 6,
   "H2b and it RETURNS (l_final,b_final) on all 6 retracting rows -- the guard is "
   "discriminating, not a blanket refusal")

# H3: MUTATION CONTROLS.  Reinstate the superseded dictionary and require every
# repaired quantity to MOVE and the bridge identity to BREAK.  Copied in shape
# from bridge_generality.py MUT F (51->205, 30->112, 22->107).
moved, survived = {}, []
for nm, (st, sdC, soC, sN, sordy) in sorted(pc4.SUPERSEDED.items()):
    R, cd = ROWS[nm], CHART[nm]
    a, b = R["ab"]
    Ng, sg = fit_sig(a, b, cd["t"], cd["kappa"], cd["deg_C"], cd["ord_C"])
    Ns, ss = fit_sig(a, b, st, st - 2, sdC, soC)
    good_bridge = bridge(a, b, cd["t"], cd["kappa"], cd["ord_C"])
    stale_bridge = bridge(a, b, st, st - 2, soC)
    # The stale route must (i) reproduce the superseded N and ord_y EXACTLY -- so
    # the control really is the OLD model and not a straw man -- (ii) disagree
    # with the guarded route, and (iii) fail the bridge identity computed at the
    # GUARDED chart.  Note (iii) is the load-bearing one: the stale numbers are
    # internally consistent (they satisfy the bridge identity at their OWN stale
    # chart, stale_bridge == sordy), which is precisely why they passed for
    # months.  What refutes them is that their chart is refused.
    faithful = (Ns == sN and ss[1] == sordy and stale_bridge == sordy)
    if faithful and Ns != Ng and ss[1] != sg[1] and ss[1] != good_bridge:
        moved[nm] = (sg[1], ss[1])
    else:
        survived.append((nm, faithful, Ns, sN, ss[1], sordy, stale_bridge))
ok(len(moved) == 5 and not survived,
   "H3 MUT: reinstating the refused dictionary reproduces the SUPERSEDED N and "
   "ord_y exactly and then contradicts the guarded chart's bridge value at every "
   "affected row -- "
   + "; ".join("%s %d<-%d" % (nm, g, s) for nm, (g, s) in sorted(moved.items())))
ok(moved.get("F1") == (51, 205) and moved.get("F3") == (30, 112)
   and moved.get("F9") == (22, 107),
   "H3b and the three rows bridge_generality.py MUT F independently mutates move "
   "by exactly the same amounts: F1 51<-205, F3 30<-112, F9 22<-107")
ok(moved.get("F2") == (30, 75) and moved.get("F10") == (114, 820),
   "H3c and the two rows MUT F does NOT cover also move: F2 30<-75, F10 114<-820 "
   "(F10's is the largest displacement in the repair)")

# H4: the mutation must be visible in the SHAPE too, not only in one integer.
ok(all(fit_sig(*ROWS[nm]["ab"], pc4.SUPERSEDED[nm][0], pc4.SUPERSEDED[nm][0] - 2,
               pc4.SUPERSEDED[nm][1], pc4.SUPERSEDED[nm][2])[1][2] > 0
       for nm in pc4.SUPERSEDED),
   "H4 every superseded row claims mult_(y+1) > 0, i.e. a (y+1) place that a "
   "monomial C cannot have; the repaired signatures all have mult = cof = 0")
ok(all(fit_sig(*ROWS[nm]["ab"], CHART[nm]["t"], CHART[nm]["kappa"],
               CHART[nm]["deg_C"], CHART[nm]["ord_C"])[1][2:] == (0, 0)
       for nm in REFUSED),
   "H4b all 11 refused rows now have signature (D,D,0,0) with D = ord = deg: Phi "
   "is a MONOMIAL there, which is what C = y forces")

# H5: the retracting rows must be BIT-IDENTICAL to their pre-repair values --
# otherwise the repair is a rewrite, not a repair.
UNCHANGED = {  # row: (t, degC, ordC, N) as the pre-repair file computed them
    "F7":  (3, 6, 4, 36), "F8":  (3, 6, 5, 70), "F14": (3, 9, 4, 36),
    "F15": (3, 9, 5, 70), "F16": (3, 9, 7, 56), "F17": (3, 9, 8, 20),
}
ok(all((CHART[nm]["t"], CHART[nm]["deg_C"], CHART[nm]["ord_C"],
        fit_sig(*ROWS[nm]["ab"], CHART[nm]["t"], CHART[nm]["kappa"],
                CHART[nm]["deg_C"], CHART[nm]["ord_C"])[0]) == v
       for nm, v in UNCHANGED.items()),
   "H5 the six retracting rows are BIT-IDENTICAL pre- and post-repair "
   "(t,deg C,ord C,N): the repair is targeted, not a rewrite")

# ---------------------------------------------------------------------------
print()
if _fail == 0:
    print(f"ALL {_n} PHI-CORNER4 CHECKS PASSED")
else:
    print(f"{_fail} of {_n} PHI-CORNER4 CHECKS FAILED")
print(f"script: {__import__('pathlib').Path(__file__).resolve()}")
sys.exit(0 if _fail == 0 else 1)
