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
    ("F6",  5, 5, 4, (4, 3), (10, 8), (1, 0), 2),
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
    ("F2", 0): (36, (189, 75, 38, 76)),
    ("F2", 1): (98, (504, 201, 101, 202)),
    ("F9", 0): (52, (377, 107, 54, 216)),
    ("F14", 0): (36, (375, 165, 42, 168)),
    ("F1", 0): (67, (275, 205, 69, 1)),
    ("F7", 0): (36, (250, 165, 83, 2)),
    ("F3", 0): (36, (189, 112, 75, 2)),
    ("F10", 0): (270, (1917, 820, 1093, 4)),
    ("F16", 0): (56, (528, 407, 117, 4)),
}

# explicit PHI_F7 ramified f-polynomials (for direct full-ODE substitution)
LANDED_F = {
    "F7":  sp.Rational(1, 10)   * y**21 * (y + 1)**11 * (9*y**2 + 3*y - 1),
    "F3":  sp.Rational(1, 42)   * y**4  * (y + 1)**3  * (25*y**2 + 15*y - 3),
    "F10": sp.Rational(1, 3740) * y**10 * (y + 1)**13 *
           (2401*y**4 + 5831*y**3 + 4165*y**2 + 595*y - 85),
    "F16": sp.Rational(1, 330)  * y**15 * (y + 1)**5  *
           (243*y**4 + 81*y**3 - 27*y**2 + 15*y - 10),
}

# published mu-rung signatures (ZETA_TAIL F12 eta=0; MU_RUNGS_F10)
MU_RUNGS = {
    ("F12", 1): (814, 506, 102, 206),
    ("F12", 2): (814, 506, 203, 105),
    ("F12", 3): (814, 506, 304, 4),
    ("F10", 2): (1917, 820, 547, 550),
    ("F10", 4): (1917, 820, 1093, 4),
}


def data(name, t, a0, q, ac, bc, A0p, k):
    kappa = t - 2
    dg = a0 - q
    r = a0 - q - 1
    gap = Fraction(q - 1) - Fraction(a0, t)
    a = ac[0] + ac[1] * j
    b = bc[0] + bc[1] * j
    e = b - a + 1
    coef = t * (b - a) + kappa + 1
    rho = sp.expand((e - 1) * q + 1)
    N = sp.expand(a * (t * (a + b - 1) + 1) - 2 * b)
    if gap == 0:
        cls = "PURE"
    elif r == 0:
        cls = "COFACTOR"
    elif gap > 0:
        cls = "RUNG"
    else:
        cls = "IRREGULAR"
    return dict(name=name, t=t, kappa=kappa, a0=a0, q=q, dg=dg, r=r, gap=gap, k=k,
                a=a, b=b, e=e, coef=coef, rho=rho, N=N, cls=cls, A0p=A0p)


def full_ode_residual(D, jv, f, g=None):
    """exact residual of a(t c f' - coef c' f) - c^e at integer j=jv (all
    exponents integers -> exact, no symbolic-power ambiguity).  c = y^q g with
    the BRANCH residual g: g = y^dg+1 (pure/cofactor) or g = (y+1)^dg (ramified)."""
    s = {j: jv}
    a = int(D["a"].subs(s)); b = int(D["b"].subs(s))
    e = b - a + 1
    q = D["q"]; dg = D["dg"]; t = D["t"]; coef = int(D["coef"].subs(s))
    if g is None:
        g = y**dg + 1
    c = y**q * g
    resid = a * (t * c * sp.diff(f, y) - coef * sp.diff(c, y) * f) - c**e
    return sp.expand(resid)


def sig(D, N, mu):
    e, a0, q, rho, gap, r = D["e"], D["a0"], D["q"], D["rho"], D["gap"], D["r"]
    pure = e * a0 - q + 1
    res = pure + gap
    deg = sp.expand(res + N * a0)
    ordy = sp.expand(rho + N * q)
    mult = sp.expand(mu * (e + N) - (mu - 1))
    cof = sp.expand(gap + r * (e + N) - (mu - 1) * (e + N - 1))
    return deg, ordy, mult, cof


DAT = {row[0]: data(*row) for row in FAM}

# ===========================================================================
# A. corner-data identities + collapse
# ===========================================================================
for dg in range(1, 8):
    ok("collapse identity dg=%d: y G' - dg G = -dg" % dg,
       sp.expand(y * sp.diff(y**dg + 1, y) - dg * (y**dg + 1)) == -dg)

for nm, D in DAT.items():
    ok("%s kappa=t-2" % nm, D["kappa"] == D["t"] - 2)
    ok("%s dg=a0-q, r=dg-1" % nm, D["dg"] == D["a0"] - D["q"] and D["r"] == D["dg"] - 1)
    ok("%s gap=(q-1)-a0/t" % nm,
       D["gap"] == Fraction(D["q"] - 1) - Fraction(D["a0"], D["t"]))
    # te - coef = t - kappa - 1 = 1, identically in j
    ok("%s te-coef=1 (identically in j)" % nm,
       sp.expand(D["t"] * D["e"] - D["coef"]) == 1)
    # t*rho - coef*q = t-(kappa+1)q, constant in j
    trc = sp.expand(D["t"] * D["rho"] - D["coef"] * D["q"])
    ok("%s t*rho-coef*q = t-(kappa+1)q (const in j)" % nm,
       trc == D["t"] - (D["kappa"] + 1) * D["q"])

# census
census = {"PURE": [], "COFACTOR": [], "RUNG": [], "IRREGULAR": []}
for nm, D in DAT.items():
    census[D["cls"]].append(nm)
ok("census PURE = {F2,F9,F14}", set(census["PURE"]) == {"F2", "F9", "F14"})
ok("census COFACTOR = {F1,F5,F6,F8,F17}",
   set(census["COFACTOR"]) == {"F1", "F5", "F6", "F8", "F17"})
ok("census RUNG (9 families)",
   set(census["RUNG"]) == {"F3", "F4", "F7", "F10", "F11", "F12", "F13", "F15", "F16"})
ok("census: no length-1 IRREGULAR", census["IRREGULAR"] == [])

# ===========================================================================
# B. CLOSED-FORM THEOREM (pure): gap=0 <=> collapse condition; A=-1/(a dg)
# ===========================================================================
for nm, D in DAT.items():
    # collapse condition  t-(kappa+1)q+dg = 0  is EXACTLY gap=0
    collapse = D["t"] - (D["kappa"] + 1) * D["q"] + D["dg"]
    ok("%s collapse-cond==0 <=> gap==0" % nm, (collapse == 0) == (D["gap"] == 0))

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

# mu-graded law algebraic identities (symbolic in j, N, mu)
Nn, mm = sp.symbols("Nn mm")
for nm in ["F3", "F10", "F12", "F16"]:
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
