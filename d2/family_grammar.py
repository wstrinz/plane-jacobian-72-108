#!/usr/bin/env python3
"""family_grammar.py  (NEW; read-only over all existing artifacts)

THE FAMILY GRAMMAR SWEEP -- generalize the F2 symbolic closed form
(f2_family_verify.py) to EVERY family in the GGV5 v11<=35 survey.

Central theorem proved here (symbolic in the family parameter j):

  For a corner (t, kappa=t-2, a0, q) the forcing ODE
      a { t c f' - [t(b-a)+kappa+1] c' f } = c^(b-a+1),   c = y^q g
  admits the F2-style UNIFORM ansatz  f = A y^rho g^e  (rho=(e-1)q+1,
  g = y^dg+1, dg=a0-q, e=b-a+1)  with a SINGLE closed form valid for all j
  IF AND ONLY IF  gap := (q-1) - a0/t = 0  (equivalently a0 = t(q-1)).
  Then A = -1/(a*dg), driven entirely by the collapse identity
      y g' - dg g = -dg      (g = y^dg + 1, always).

The trichotomy (derived, not assumed):

  * gap = 0            -> CLOSED-FORM (pure):     f = -1/(a dg) y^rho (y^dg+1)^e
  * r = a0-q-1 = 0     -> CLOSED-FORM (cofactor): dg=1 forced, f = y^rho (y+1)^e u,
    with gap > 0          u a UNIT cofactor of degree gap, coeffs rational in j
  * r > 0 AND gap > 0  -> RUNG-STRUCTURED: the unramified ansatz fails (g=y^dg+1
                          has NO simple root at -1 when dg is even, or the wrong
                          residual when dg is odd), but the mu = dg fully-ramified
                          rung DOES give a uniform closed form
                          f = y^rho (y+1)^(dg e-(dg-1)) u, deg u = gap+r.
  * gap < 0            -> IRREGULAR: res < pure, regime unobserved (F22 escape).

Every claimed closed form is re-verified symbolically in j and every landed
derived point is reproduced exactly by family_grammar_verify.py.

Complex-scope discipline (MU_RUNGS_F10 correction, adopted): branch schemes are
classified over Qbar.  dg even => the mu=1 REAL locus is empty (PHI_F7 at dg=2
is a complete factorization over C; MU_RUNGS_F10 at dg=4 is REAL-only); this is
branch-selection annotation, never a complex kill.  The mu=dg ramified rung is
rational and is the uniform representative we derive.

Sources (transcribed, not re-derived): GGV5 v11<=35 tables via phi_corner4.py;
F2 template f2_family_verify.py; ramified law PHI_F7.md; mu-ladder ZETA_TAIL.md
/ MU_RUNGS_F10.md; escapes composite_charts.py.  Exact sympy throughout.
"""
import sys

import sympy as sp
from fractions import Fraction

y = sp.symbols("y", positive=True)
j = sp.symbols("j", nonnegative=True, integer=True)

BAR = "=" * 100

# ---------------------------------------------------------------------------
# 0. Family table.  (name, t=l, a0, q, a(j)=(a0c,dac), b(j)=(b0c,dbc), A0prime, k)
#    a = min(m,n), b = max(m,n) as functions of j (derived from the survey's
#    (m,n)(j); order verified for j>=0).  kappa = t-2 (fused-chart lemma).
# ---------------------------------------------------------------------------
FAMILIES = [
    # name    t  a0  q     a(j)      b(j)     A0'     k
    ("F1",   4,  4, 3,  (3, 2),  (4, 3),  (1, 0), 1),
    ("F2",   5,  5, 2,  (2, 1),  (3, 2),  (1, 0), 1),
    ("F3",   5,  5, 3,  (2, 3),  (3, 4),  (1, 0), 1),
    ("F4",   5,  5, 3,  (3, 2),  (16, 12),(1, 0), 2),
    ("F5",   5,  5, 4,  (5, 4),  (9, 7),  (1, 0), 1),
    ("F6",   5,  5, 4,  (7, 6),  (18, 16),(1, 0), 2),  # CORRECTED 2026-07-24: GGV5 F6 typo (a,b) base (4,10) gcd=2; coprime family a=6j+7, b=16j+18 (CHAIN_SURVEY.md)
    ("F7",   3,  6, 4,  (2, 1),  (7, 4),  (1, 0), 1),
    ("F8",   3,  6, 5,  (3, 2),  (7, 5),  (1, 0), 1),
    ("F9",   7,  7, 2,  (2, 1),  (3, 2),  (1, 0), 1),
    ("F10",  7,  7, 3,  (4, 3),  (7, 5),  (1, 0), 1),
    ("F11",  7,  7, 3,  (2, 1),  (5, 3),  (1, 0), 2),
    ("F12",  4,  8, 5,  (3, 2),  (7, 5),  (2, 0), 1),
    ("F13",  3,  9, 7,  (2, 1),  (13, 7), (2, 0), 1),
    ("F14",  3,  9, 4,  (2, 1),  (7, 4),  (1, 0), 1),
    ("F15",  3,  9, 5,  (3, 2),  (7, 5),  (1, 0), 1),
    ("F16",  3,  9, 7,  (3, 4),  (5, 7),  (1, 0), 1),
    ("F17",  3,  9, 8,  (2, 5),  (3, 8),  (1, 0), 1),
]

# length-2 composite escapes (chart underived: CONDITIONAL, treated lightly)
ESCAPES = [
    # name    t  a0  q     (a,b) at j=0 only    note
    ("F22",  4,  8, 2,  (2, 3),  "gap<0 IRREGULAR"),
    ("F23",  4,  8, 4,  (2, 7),  "gap>0,r>0 RUNG (conditional)"),
    ("F24",  8,  8, 3,  (3, 4),  "gap>0,r>0 RUNG dg odd (conditional)"),
]

# landed / published derived points to reproduce (case, family, j, signature)
LANDED = {
    # PURE (gap=0)
    # 2026-07-26 REPAIR (PASSPORT_75_125_REPAIR.md): the (5,20) corner has t=4, kappa=2, C=y (deg C=ord C=1), NOT t=5, kappa=3, C=y^2(y^3+1).  GGV5's final chain corner (7\\5,2) is chart data only on the retraction shape b0=l(a0-1), which (5,20) fails; l = ceil(20/5) = 4.  Both F2 rows below move: (75,125) N 98->77, sig (504,201,101,202)->(80,80,0,0), lc -1/9->1/3; (50,75) N 36->28, sig (189,75,38,76)->(30,30,0,0), lc -1/6->1/2.  The guard lives in polygon_reduction.py sec.0b.
    ("F2", 0): ("(50,75)",  28, (30, 30, 0, 0)),      # REPAIRED 2026-07-26
    ("F2", 1): ("(75,125)", 77, (80, 80, 0, 0)),      # REPAIRED 2026-07-26
    ("F9", 0): ("(56,84)",  52, (377, 107, 54, 216)),
    ("F14", 0): ("(66,231)", 36, (375, 165, 42, 168)),
    # CHART-DEGENERATE (dg=0), same corner (5,20) as F2 -- see REPAIRS["F3"]
    # 2026-07-26 (second repair): F3 j=0 is (75,50), the (m,n)-SWAP of F2 j=0's
    # (50,75) at the SAME corner, with the SAME unordered reduced pair {2,3}.
    # The corner law depends on (t,kappa,deg C,ord C) and on {min(m,n),max(m,n)}
    # only, so F3 j=0 MUST land on F2 j=0's numbers exactly: N 36 -> 28,
    # sig (189,112,75,2) -> (30,30,0,0).  Verified three independent ways in
    # family_grammar_verify.py A10.
    ("F3", 0): ("(75,50)",  28, (30, 30, 0, 0)),      # REPAIRED 2026-07-26
    # COFACTOR (r=0)
    ("F1", 0): ("(48,64)",  67, (275, 205, 69, 1)),
    # RUNG mu=dg ramified (PHI_F7)
    ("F7", 0): ("(42,147)", 36, (250, 165, 83, 2)),
    ("F10", 0): ("(196,112)", 270, (1917, 820, 1093, 4)),
    ("F16", 0): ("(99,165)", 56, (528, 407, 117, 4)),
}

# explicit landed f-polynomials (PHI_F7 / MU_RUNGS) for direct ODE re-check
LANDED_F = {
    "F7":  sp.Rational(1, 10)   * y**21 * (y + 1)**11 * (9*y**2 + 3*y - 1),
    "F10": sp.Rational(1, 3740) * y**10 * (y + 1)**13 *
           (2401*y**4 + 5831*y**3 + 4165*y**2 + 595*y - 85),
    "F16": sp.Rational(1, 330)  * y**15 * (y + 1)**5  *
           (243*y**4 + 81*y**3 - 27*y**2 + 15*y - 10),
}

# 2026-07-26.  F3's PHI_F7 rung polynomial (phi_f7.py:32) is NOT retired as a
# computation -- it exactly solves the mu=2 ramified ODE of the corner data
# (t,kappa,q,dg) = (5,3,3,2).  It is retired as a statement ABOUT (5,20):
# that corner's chart is (t,kappa,deg C,ord C) = (4,2,1,1), and this polynomial
# does not solve the resulting ODE.  Kept here, labelled, so the discriminating
# check in family_grammar_verify.py A10d has something to run against; deleting
# it would make the repair unfalsifiable.
SUPERSEDED_F = {
    "F3":  sp.Rational(1, 42)   * y**4  * (y + 1)**3  * (25*y**2 + 15*y - 3),
}

# F12 mu-rungs at eta=0 (ZETA_TAIL) and F10 all real rungs (MU_RUNGS)
MU_RUNGS = {
    # (family, eta, mu) : signature
    ("F12", 0, 1): (814, 506, 102, 206),
    ("F12", 0, 2): (814, 506, 203, 105),
    ("F12", 0, 3): (814, 506, 304, 4),
    ("F10", 0, 2): (1917, 820, 547, 550),
    ("F10", 0, 4): (1917, 820, 1093, 4),
}


# ---------------------------------------------------------------------------
# derived per-family quantities
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# PER-FAMILY CHART REPAIRS (2026-07-26).
#
# A family whose corner does NOT satisfy the retraction shape b0 = t(a0-1) has
# deg_y(C)/ord_y(C) DIFFERENT from (a0, q): GGV5's final-corner dictionary, which
# is what made them equal, is invalid off that shape.  Only rows listed here are
# affected; every other family keeps the historical defaults exactly.
#
# F2 -- corner (5,20).  20 != 4*4, so it does not retract.  Repaired chart:
# t = 4 (was 5, whence kappa = t-2 = 2, matching PASSPORT_75_125_REPAIR.md), and
# C = y a MONOMIAL so deg C = ord C = 1 (was 5 and 2).  gap stays 0 -- the family
# is still PURE -- but it can no longer be computed as (q-1) - a0/t, because that
# formula also assumed the dictionary.  Verified: this reproduces BOTH landed
# points exactly and symbolically in j, N_j = (3j+4)(4j+7) and
# deg = ord = 2(2j+3)(3j+5) = 12j^2+38j+30, which independently agrees with
# window_functions_75_125.family()'s ord_y(Phi) = 12a^2-10a+2 at a = j+2.
# Populated by the landed-point cross-check; a non-empty list is FATAL (see the
# end of main).  Until 2026-07-26 a MISMATCH was printed and swallowed, which is
# how both F2 rows stayed contradictory through a chart repair AND a green suite.
#
# F3 -- SAME corner (5,20), repaired the same way later on 2026-07-26.  The
# operative principle, now established rather than assumed:
#
#     CHART DATA IS A PROPERTY OF THE CORNER (A_0, A_0'), NOT OF THE FAMILY.
#
# l = chart_exponent(a0,b0) = ceil(b0/a0), kappa = l-2, and "is the top face
# vertical" (which decides whether C is a monomial) are all functions of A_0
# alone; the family row contributes only (m,n).  F2 and F3 both have
# A_0 = (5,20), A_0' = (1,0) and chain length 1, so they MUST share
# (t,kappa,deg C,ord C) = (4,2,1,1) -- and polygon_reduction.corner_chart_data
# returns bit-identical dicts for the two rows (checked, A10 below).
#
# That is also why the OLD data was self-refuting before any paper was consulted:
# it gave ord C = q = 2 for F2 and ord C = q = 3 for F3 AT ONE AND THE SAME
# CORNER.  Two different values of a corner invariant at one corner is a
# contradiction, and it is the fingerprint of reading q off GGV5's per-row final
# chain corner b_final (2 for F2's (7\5,2), 3 for F3's (8\5,3)) instead of off
# the corner's own reduced polygon.  So q = 3 was never chart data for F3.
#
# Verified target for F3 (three routes, family_grammar_verify.py A10a-A10f):
#   (i)   polygon_reduction._f2_forcing_divisor(2,3,4,2,1,1) -> N=28, (30,30,0,0);
#   (ii)  window_functions_75_125.family(2)   -> N=28, ord=deg=30;
#   (iii) GGV3 sec.5 (1406.0886_GGV3.tex:1723-1727) publishes, for A_0=(5,20),
#         [P_1,Q_1] = x^2, deg P_1 = 10, deg Q_1 = 15.  Those three integers are
#         obtained BEFORE the paper's gamma in {2,3} branch (both branches start
#         from the same (P_1,Q_1)), so they are not specific to one GGV5 row, and
#         F3(3,2)/75 is that same case with P and Q exchanged: [Q_1,P_1] = -x^2,
#         degrees (15,10).  Either reading forces kappa = 2, hence l = 4.
MISMATCHES = []

REPAIRS = {
    "F2": dict(t=4, degC=1, ordC=1, gap0=0),
    "F3": dict(t=4, degC=1, ordC=1, gap0=0),
}


def lin(coef, jj):
    return coef[0] + coef[1] * jj


def family_data(name, t, a0, q, ac, bc, A0p, k, degC=None, ordC=None, gap0=None):
    """degC/ordC/gap0 REPAIR HOOK, added 2026-07-26.

    This function used `a0` for TWO different things -- the corner's first
    coordinate AND deg_y(C) -- and `q` for both a chart datum and ord_y(C).
    Those coincide only under GGV5's final-corner dictionary, which is valid
    ONLY on the retraction shape b0 = t(a0-1).  Off that shape they diverge:
    at (5,20), deg C = ord C = 1 (C is the MONOMIAL y) while a0 = 5.

    The three optional arguments default to the historical behaviour
    (degC = a0, ordC = q, gap computed), so EVERY family row not passing them
    is bit-for-bit unchanged.  F2 passes them because its corner does not
    retract.  See PASSPORT_75_125_REPAIR.md and polygon_reduction.py sec.0b.
    """
    kappa = t - 2
    degC0 = a0 if degC is None else degC
    ordC0 = q if ordC is None else ordC
    # dg is the RESIDUAL's degree (g = y^dg + 1), i.e. deg C - ord C.  That equals
    # a0 - q only when the final-corner dictionary holds.  dg = 0 means g is a
    # constant: there is NO residual, which is the repaired (5,20) case (C = y).
    dg = degC0 - ordC0
    r = dg - 1
    gap = Fraction(q - 1) - Fraction(a0, t) if gap0 is None else Fraction(gap0)
    degC = a0 if degC is None else degC
    ordC = q if ordC is None else ordC
    a = lin(ac, j)
    b = lin(bc, j)
    e = b - a + 1
    coef = t * (b - a) + kappa + 1
    rho = (e - 1) * ordC + 1
    N = a * (t * (a + b - 1) + 1) - 2 * b
    if dg == 0:
        # C is a MONOMIAL: g = y^0+1 is constant, so the pure closed form
        # f = -1/(a*dg) y^rho (y^dg+1)^e is UNDEFINED (A = -1/(a*0)) and the
        # collapse identity y g' - dg g = -dg reads 0 = 0, pinning nothing.  The
        # landed signatures still reproduce; the MECHANISM does not exist.
        cls = "CHART-DEGENERATE"
    elif gap == 0:
        cls = "CLOSED-FORM (pure)"
    elif r == 0:
        cls = "CLOSED-FORM (cofactor)"
    elif gap > 0:
        cls = "RUNG-STRUCTURED"
    else:
        cls = "IRREGULAR"
    return dict(name=name, t=t, kappa=kappa, a0=a0, q=q, dg=dg, r=r, gap=gap,
                a=a, b=b, e=e, coef=coef, rho=rho, N=sp.expand(N), k=k,
                A0p=A0p, cls=cls, ac=ac, bc=bc, degC=degC, ordC=ordC)


# ---------------------------------------------------------------------------
# solvers (all exact, symbolic in j)
# ---------------------------------------------------------------------------
def pure_residual(D):
    """ANALYTIC reduced ODE residual of f = -1/(a dg) y^rho (y^dg+1)^e, symbolic in j.

    Using  a(t c f' - coef c' f) - c^e = y^(qe) G^e * R(G)  with G = y^dg+1 and the
    collapse identity  y G' = dg G - dg,  the whole residual reduces to the scalar
        R(G) = A a[(t*rho - coef*q) + (te-coef)*dg] G  -  A a (te-coef) dg  -  1
    (a degree<=1 polynomial in the OPAQUE symbol G, coeffs rational in j).  Symbolic-
    exponent cancellation is thus avoided entirely; R(G) is identically 0 iff the
    closed form solves the ODE for all j.  Returns R(G).
    """
    dg, e, rho, a, t, coef, q = D["dg"], D["e"], D["rho"], D["a"], D["t"], D["coef"], D["q"]
    G = sp.symbols("G")
    A = -sp.Rational(1, 1) / (a * dg)
    trc = sp.expand(t * rho - coef * q)     # = t - (kappa+1) q  (constant in j)
    tec = sp.expand(t * e - coef)           # = t - kappa - 1 = 1
    R = A * a * ((trc + tec * dg) * G - tec * dg) - 1
    return sp.expand(R)


def cofactor_solve(D):
    """solve unit cofactor u (deg gap) for r=0 (dg=1); returns (u, sol)."""
    dg, e, rho, a, t, coef, q, gap = (D["dg"], D["e"], D["rho"], D["a"],
                                      D["t"], D["coef"], D["q"], int(D["gap"]))
    g = y**dg + 1
    gp = sp.diff(g, y)
    uc = sp.symbols("u0:%d" % (gap + 1))
    u = sum(uc[i] * y**i for i in range(gap + 1))
    f_over = y**rho * u
    tcf = a * t * y**q * g * (sp.diff(f_over, y) * g + e * f_over * gp)
    ccf = a * coef * (q * y**(q - 1) * g + y**q * gp) * f_over * g
    src = y**(q * e) * g
    red = sp.expand((tcf - ccf - src) / y**(q * e))
    P = sp.Poly(red, y)
    sol = sp.solve([sp.expand(co) for co in P.all_coeffs()], list(uc), dict=True)
    return (u.subs(sol[0]) if sol else None), (sol[0] if sol else None), uc


def ramified_solve(D, mu=None):
    """solve mu=dg fully-ramified rung f = y^rho (y+1)^(dg e-(dg-1)) u,
       deg u = gap+r ; returns u (symbolic in j)."""
    dg, e, rho, a, t, coef, q, gap, r = (D["dg"], D["e"], D["rho"], D["a"], D["t"],
                                         D["coef"], D["q"], int(D["gap"]), D["r"])
    m = dg if mu is None else mu
    P = m * e - (m - 1)
    dego = gap + r
    uc = sp.symbols("w0:%d" % (dego + 1))
    u = sum(uc[i] * y**i for i in range(dego + 1))
    up = sp.diff(u, y)
    # reduced residual after dividing by y^(qe)(y+1)^(dg e) with g=(y+1)^dg
    red = (a * t * (y * (y + 1) * up + rho * (y + 1) * u + P * y * u)
           - a * coef * (q * (y + 1) + dg * y) * u - 1)
    red = sp.expand(red)
    Pp = sp.Poly(red, y)
    sol = sp.solve([sp.expand(co) for co in Pp.all_coeffs()], list(uc), dict=True)
    return (u.subs(sol[0]) if sol else None), P, dego


def signature(D, N, mu):
    """mu-graded corner law signature (deg, ord, mult, cof)."""
    e, a0, q, rho, gap, r, dg = D["e"], D["a0"], D["q"], D["rho"], D["gap"], D["r"], D["dg"]
    # deg_y/ord_y of Phi are driven by deg_y(C)/ord_y(C), NOT by the corner
    # coordinates.  They default to (a0, q) so every retracting family is
    # unchanged; see family_data's repair hook.
    degC, ordC = D.get("degC", a0), D.get("ordC", q)
    pure = e * degC - ordC + 1
    res = pure + gap
    deg = sp.expand(res + N * degC)
    ordy = sp.expand(rho + N * ordC)
    if dg == 0:
        # C is a MONOMIAL: g = y^0 + 1 is a constant, so there is no residual to
        # carry a multiplicity or a cofactor.  Both are identically 0 -- which is
        # exactly the shape of the repaired (5,20) signatures (30,30,0,0) and
        # (80,80,0,0), and is why they have two trailing zeros at all.
        return deg, ordy, sp.Integer(0), sp.Integer(0)
    mult = sp.expand(mu * (e + N) - (mu - 1))
    cof = sp.expand(gap + r * (e + N) - (mu - 1) * (e + N - 1))
    return deg, ordy, mult, cof


def wstep_denom(D, N):
    """W_step = ord_y(Phi)/M, M = t(a+b)-(kappa+1); reduced denominator in j.

    2026-07-26: ord_y(Phi) = rho + N*ord_y(C), so the multiplier is ord C, not
    the chain datum q.  Identical for every family where the final-corner
    dictionary holds (there ord C == q); it differs exactly on the repaired
    corners, where it now agrees with f2_family_verify.py check E
    (a=2: 30/17, a=3: 80/29).
    """
    ordy = D["rho"] + N * D["ordC"]
    M = D["t"] * (D["a"] + D["b"]) - (D["kappa"] + 1)
    W = sp.cancel(ordy / M)
    return sp.denom(W)


# ===========================================================================
# main derivation report
# ===========================================================================
if __name__ == "__main__":
    print(BAR)
    print("THE FAMILY GRAMMAR SWEEP -- 17 length-1 families + 3 length-2 escapes")
    print(BAR)
    print("Collapse identity (all families):  y*(y^dg)' - dg*(y^dg+1) = -dg  ==>")
    for dg in (1, 2, 3, 4, 5):
        val = sp.expand(y * sp.diff(y**dg + 1, y) - dg * (y**dg + 1))
        print(f"   dg={dg}:  y g' - {dg} g = {val}")
    print("\nCollapse THEOREM: pure ansatz f=A y^rho (y^dg+1)^e solves the ODE for")
    print("all e  <=>  t - (kappa+1)q + dg = 0  <=>  a0 = t(q-1)  <=>  gap = 0,")
    print("and then A = -1/(a*dg).  (te - coef = t-kappa-1 = 1 for the whole class.)")

    data = {}
    for row in FAMILIES:
        name = row[0]
        rep = REPAIRS.get(name, {})
        if "t" in rep:
            row = (row[0], rep["t"]) + tuple(row[2:])
        D = family_data(*row, **{k: v for k, v in rep.items() if k != "t"})
        data[D["name"]] = D

    # ---- census + grammar table ----
    print("\n" + BAR)
    print("GRAMMAR TABLE (per-family local obstruction type)")
    print(BAR)
    hdr = (f"{'fam':4}{'t':>2}{'kap':>4}{'a0':>3}{'q':>2}{'dg':>3}{'r':>2}"
           f"{'gap':>5}{'k':>2}  {'class':22} {'N_j':>26}")
    print(hdr)
    for nm in [r[0] for r in FAMILIES]:
        D = data[nm]
        print(f"{D['name']:4}{D['t']:>2}{D['kappa']:>4}{D['a0']:>3}{D['q']:>2}"
              f"{D['dg']:>3}{D['r']:>2}{str(D['gap']):>5}{D['k']:>2}  "
              f"{D['cls']:22} {str(sp.factor(D['N'])):>26}")

    census = {"CLOSED-FORM (pure)": [], "CLOSED-FORM (cofactor)": [],
              "RUNG-STRUCTURED": [], "IRREGULAR": [], "CHART-DEGENERATE": []}
    for nm in [r[0] for r in FAMILIES]:
        census[data[nm]["cls"]].append(nm)
    print("\nCENSUS (length-1):")
    for c, lst in census.items():
        print(f"   {c:24}: {len(lst):2}  {', '.join(lst)}")

    # ---- per-family derivations ----
    for nm in [r[0] for r in FAMILIES]:
        D = data[nm]
        print("\n" + BAR)
        print(f"{nm}  --  {D['cls']}   (t={D['t']} kappa={D['kappa']} a0={D['a0']} "
              f"q={D['q']} dg={D['dg']} r={D['r']} gap={D['gap']})")
        print(BAR)
        print(f"   a(j) = {D['a']},  b(j) = {D['b']},  e(j) = {sp.expand(D['e'])},  "
              f"rho(j) = {sp.expand(D['rho'])}")
        print(f"   N_j = {sp.factor(D['N'])}   (block step Delta N = "
              f"{sp.expand(D['N'].subs(j, j+1) - D['N'])}, block size t = {D['t']})")
        Nsym = D["N"]
        if D["cls"] == "CHART-DEGENERATE":
            # Printing NOTHING here is what let F2 sit silently in the report
            # after 86d8fb0.  State the withdrawal explicitly instead.
            deg, ordy, mult, cof = signature(D, Nsym, 1)
            print(f"   corner does NOT retract: deg C = ord C = {D['degC']}, so "
                  f"C is a MONOMIAL and dg = 0 (a0 = {D['a0']} is the corner's "
                  f"first coordinate, NOT deg C).")
            print("   NO closed form is claimed: f = -1/(a*dg) y^rho (y^dg+1)^e "
                  "has constant -1/(a*0) and the")
            print("   collapse identity y g' - dg g = -dg degenerates to 0 = 0, "
                  "so it pins nothing.")
            print(f"   Phi_j is a MONOMIAL: deg = ord = {sp.factor(ordy)}, "
                  f"mult = cof = 0.")
            print(f"   W_step denominator: {sp.factor(wstep_denom(D, Nsym))}")
        elif D["cls"] == "CLOSED-FORM (pure)":
            res = pure_residual(D)
            print(f"   f_j = -1/(({D['a']})*{D['dg']}) y^rho (y^{D['dg']}+1)^e   "
                  f"[A = -1/(a*dg) = {sp.simplify(-sp.Rational(1,1)/(D['a']*D['dg']))}]")
            print(f"   ODE residual (symbolic in j): {res}   -> "
                  f"{'VANISHES' if res == 0 else 'NONZERO'}")
            deg, ordy, mult, cof = signature(D, Nsym, 1)
            print(f"   Phi_j signature (mu=1): deg={sp.factor(deg)}, ord={sp.factor(ordy)},")
            print(f"       mult=e+N={sp.factor(mult)}, cof=r(e+N)={sp.factor(cof)}")
            print(f"   W_step denominator (analogue of F2's 5a-3): "
                  f"{sp.factor(wstep_denom(D, Nsym))}")
        elif D["cls"] == "CLOSED-FORM (cofactor)":
            u, sol, uc = cofactor_solve(D)
            um1 = sp.simplify(u.subs(y, -1))
            u0 = sp.simplify(u.subs(y, 0))
            print(f"   f_j = y^rho (y+1)^e * u(y),  u a UNIT cofactor of degree "
                  f"gap={int(D['gap'])}:")
            print(f"       u = {sp.simplify(u)}")
            print(f"       u(0) = {u0} (!=0),  u(-1) = {um1} (= -1/a),  UNIT ok")
            deg, ordy, mult, cof = signature(D, Nsym, 1)
            print(f"   Phi_j signature (mu=1): deg={sp.factor(deg)}, ord={sp.factor(ordy)}, "
                  f"mult=e+N={sp.factor(mult)}, cof=gap={cof}")
            print(f"   W_step denominator: {sp.factor(wstep_denom(D, Nsym))}")
        elif D["cls"] == "RUNG-STRUCTURED":
            paritynote = ("dg EVEN: mu=1 real locus empty (complex not excluded); "
                          "ramified mu=dg is the rational rung"
                          if D["dg"] % 2 == 0 else
                          "dg ODD: mu=1 (squarefree) branch available; "
                          "ramified mu=dg also derived")
            print(f"   Unramified ansatz g=y^dg+1 FAILS (gap={D['gap']}!=0). {paritynote}.")
            u, P, dego = ramified_solve(D)
            if u is not None:
                print(f"   mu=dg RAMIFIED closed form (symbolic in j):")
                print(f"       f_j = y^rho (y+1)^({sp.expand(P)}) u(y),  deg u = gap+r = {dego}")
                print(f"       u(j=0) = {sp.nsimplify(sp.expand(u.subs(j, 0)))}")
            for mu in (1, D["dg"]):
                deg, ordy, mult, cof = signature(D, Nsym, mu)
                tag = "mu=1 (unram)" if mu == 1 else f"mu=dg={D['dg']} (ramified)"
                print(f"   mu-graded sig {tag:18}: "
                      f"({sp.factor(deg)}, {sp.factor(ordy)}, "
                      f"{sp.expand(mult)}, {sp.expand(cof)})")
            if D["A0p"] != (1, 0):
                print("   [CONDITIONAL: A0'=(2,0) zeta-tail model; N-formula off-diagonal]")
            if D["k"] == 2:
                print("   [CONDITIONAL: k=2, N-formula unverified upstream]")

    # ---- cross-check landed points ----
    print("\n" + BAR)
    print("CROSS-CHECK: landed derived points reproduced by the family formulas")
    print(BAR)
    for (nm, jj), (label, Nl, sig) in LANDED.items():
        D = data[nm]
        Nval = int(D["N"].subs(j, jj))
        mu = D["dg"] if D["cls"] == "RUNG-STRUCTURED" else 1
        got = tuple(int(sp.Integer(v.subs(j, jj))) for v in signature(D, Nval, mu))
        ok = (Nval == Nl and got == sig)
        if not ok:
            MISMATCHES.append((nm, jj, label, Nval, got, Nl, sig))
        print(f"   {nm} j={jj} {label:11} N={Nval:4} sig={got}  "
              f"{'OK' if ok else 'MISMATCH vs %s' % (sig,)}")

    print("\n" + BAR)
    print("mu-RUNG signatures (ZETA_TAIL F12 eta=0; MU_RUNGS_F10) via mu-graded law")
    print(BAR)
    for (nm, eta, mu), sig in MU_RUNGS.items():
        D = data[nm]
        Nval = int(D["N"].subs(j, 0))
        got = tuple(int(sp.Integer(v.subs(j, 0))) for v in signature(D, Nval, mu))
        print(f"   {nm} eta={eta} mu={mu}: law={got}  published={sig}  "
              f"{'OK' if got == sig else 'MISMATCH'}")

    print("\n" + BAR)
    print("LENGTH-2 ESCAPES (composite charts underived -- CONDITIONAL)")
    print(BAR)
    for nm, t, a0, q, (a, b), note in ESCAPES:
        dg = a0 - q
        r = a0 - q - 1
        gap = Fraction(q - 1) - Fraction(a0, t)
        cls = ("IRREGULAR (gap<0: res<pure, unobserved)" if gap < 0 else
               "RUNG-STRUCTURED (conditional)")
        print(f"   {nm}: (a,b)=({a},{b}) t={t} a0={a0} q={q} dg={dg} r={r} "
              f"gap={gap} -> {cls}")
    print("\n(72,108) is the (8,28) resonance EXCEPTION -- r=0,gap=4 audited quartic,")
    print("NOT a generic member of these 17 families; reproduced by no family formula")
    print("(handled honestly: it is its own corner, cof=gap=4 in the r=0-amended law).")
    print("\nDERIVATION COMPLETE -- see FAMILY_GRAMMAR.md; checker family_grammar_verify.py")


    # A landed point the family formulas do NOT reproduce means the grammar and the
    # landed data disagree, and that must be FATAL.  Before 2026-07-26 the mismatch
    # was printed and the module exited 0, so both F2 rows contradicted their own
    # targets through a chart repair AND a green full-tier suite.
    if MISMATCHES:
        print(chr(10) + BAR)
        print("FATAL: %d landed point(s) NOT reproduced by the family formulas"
              % len(MISMATCHES))
        for nm, jj, label, Nval, got, Nl, sig in MISMATCHES:
            print("   %s j=%d %s: derived N=%s sig=%s  vs landed N=%s sig=%s"
                  % (nm, jj, label, Nval, got, Nl, sig))
        print(BAR)
        sys.exit(1)
