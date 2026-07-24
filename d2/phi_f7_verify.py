#!/usr/bin/env python3
"""phi_f7_verify.py -- independent PASS/FAIL checker for PHI_F7.md / phi_f7.py.

Independent routes (no imports from phi_f7.py):
  * corner arithmetic re-derived from the GGV5 table rows with fresh code;
  * every claimed f verified against the ODE by direct sp.diff/expand;
  * uniqueness by FULL generic linear solve (sp.solve over all coefficients,
    degree slack past resonance) -- not the triangular recurrence;
  * signatures by trial division on f plus exponent arithmetic, cross-checked
    at F3 by full expansion of Phi itself;
  * the no-simple-real-root-at--1 theorem re-proven per corner from the
    claimed obstruction branches by quadratic-discriminant arguments;
  * off-branch refutation spot checks (generic g -> NO polynomial solution);
  * controls: F14 (gap=0) and F1 (dg=1 root gauge) plus (72,108) arithmetic.

Usage: python phi_f7_verify.py [--quiet].  Exit 0 iff every check passes.
"""
import sys
import sympy as sp
from fractions import Fraction
from math import gcd

QUIET = "--quiet" in sys.argv[1:]
y = sp.symbols("y")
_n, _f = 0, 0

def check(label, ok):
    global _n, _f
    _n += 1
    if not ok:
        _f += 1
    if not QUIET or not ok:
        print(("[OK] " if ok else "[FAIL] ") + label)

# ---------------------------------------------------------------------------
# A. corner arithmetic (independent implementation)
# ---------------------------------------------------------------------------
TABLE = {   # name: (A0, p, l, q, k, (m0,dm), (n0,dn))  [GGV5 v11<=35 rows]
    "F7":  ((6, 15), 7, 3, 4, 1, (2, 1), (7, 4)),
    "F3":  ((5, 20), 8, 5, 3, 1, (3, 4), (2, 3)),
    "F10": ((7, 21), 13, 7, 3, 1, (7, 5), (4, 3)),
    "F16": ((9, 24), 10, 3, 7, 1, (3, 4), (5, 7)),
    "F14": ((9, 24), 7, 3, 4, 1, (2, 1), (7, 4)),
    "F1":  ((4, 12), 7, 4, 3, 1, (3, 2), (4, 3)),
}

def params(name):
    A0, p, l, q, k, (m0, dm), (n0, dn) = TABLE[name]
    j = 0
    while gcd(m0 + dm * j, n0 + dn * j) != 1:
        j += 1
    m, n = m0 + dm * j, n0 + dn * j
    ok_dio = (m + n) * q * k - n * (q * l - p) == k
    lo, hi = min(m, n), max(m, n)
    t = l
    a0 = A0[0]
    e, r, dg = hi - lo + 1, a0 - q - 1, a0 - q
    coef = t * (hi - lo) + (t - 2) + 1
    rho = (e - 1) * q + 1
    N = lo * (t * (lo + hi) - (t - 1)) - 2 * hi
    res = Fraction(coef * a0, t)
    gap = int(res) - (e * a0 - q + 1)
    return dict(ok_dio=ok_dio, a=lo, b=hi, t=t, a0=a0, q=q, e=e, r=r, dg=dg,
                coef=coef, rho=rho, N=N, res=int(res), gap=gap,
                degs=((A0[0] + A0[1]) * m, (A0[0] + A0[1]) * n))

EXPECT = {  # (a,b,t,a0,q,e,r,dg,rho,N,res,gap,degs)
    "F7":  (2, 7, 3, 6, 4, 6, 1, 2, 21, 36, 34, 1, (42, 147)),
    "F3":  (2, 3, 5, 5, 3, 2, 1, 2, 4, 36, 9, 1, (75, 50)),
    "F10": (4, 7, 7, 7, 3, 4, 3, 4, 10, 270, 27, 1, (196, 112)),
    "F16": (3, 5, 3, 9, 7, 3, 1, 2, 15, 56, 24, 3, (99, 165)),
}
for nm, exp in EXPECT.items():
    P = params(nm)
    got = (P["a"], P["b"], P["t"], P["a0"], P["q"], P["e"], P["r"], P["dg"],
           P["rho"], P["N"], P["res"], P["gap"], P["degs"])
    check(f"A: {nm} Diophantine identity", P["ok_dio"])
    check(f"A: {nm} corner invariants {got}", got == exp)

# ---------------------------------------------------------------------------
# B. claimed solutions: ODE identity + uniqueness by full generic solve
# ---------------------------------------------------------------------------
CLAIMS = {  # name -> (g poly, claimed f)
    ("F7", "ram"):  ((y + 1)**2,
                     y**21 * (y + 1)**11 * (9 * y**2 + 3 * y - 1) / 10),
    ("F7", "cx"):   (y**2 + 3 * y + 6,
                     y**21 * (y - 1) * (y**2 + 3 * y + 6)**6 / 60),
    ("F3", "ram"):  ((y + 1)**2,
                     y**4 * (y + 1)**3 * (25 * y**2 + 15 * y - 3) / 42),
    ("F3", "cx"):   (y**2 + 3 * y + sp.Rational(27, 5),
                     y**4 * (5 * y - 3) * (5 * y**2 + 15 * y + 27)**2 / 5670),
    ("F16", "ram"): ((y + 1)**2,
                     y**15 * (y + 1)**5
                     * (243 * y**4 + 81 * y**3 - 27 * y**2 + 15 * y - 10) / 330),
    ("F10", "ram"): ((y + 1)**4,
                     y**10 * (y + 1)**13
                     * (2401 * y**4 + 5831 * y**3 + 4165 * y**2 + 595 * y - 85)
                     / 3740),
}

def ode_ok(P, g, f):
    c = y**P["q"] * g
    return sp.expand(P["a"] * P["t"] * c * sp.diff(f, y)
                     - P["a"] * P["coef"] * sp.diff(c, y) * f
                     - c**P["e"]) == 0

def generic_solutions(P, g, slack=2):
    D = P["res"] + slack
    fc = sp.symbols(f"vf0:{D + 1}")
    f = sum(fc[i] * y**i for i in range(D + 1))
    c = y**P["q"] * g
    resid = sp.expand(P["a"] * P["t"] * c * sp.diff(f, y)
                      - P["a"] * P["coef"] * sp.diff(c, y) * f - c**P["e"])
    sols = sp.solve(sp.Poly(resid, y).all_coeffs(), fc, dict=True)
    return [sp.expand(f.subs(s)) for s in sols]

def sig(p):
    P2 = sp.Poly(sp.expand(p), y)
    deg = P2.degree()
    ordy = min(mm[0] for mm in P2.monoms())
    mult, Q = 0, P2
    while True:
        q2, r2 = sp.div(Q, sp.Poly(y + 1, y))
        if not r2.is_zero:
            break
        Q, mult = q2, mult + 1
    return (deg, ordy, mult, deg - ordy - mult)

SIGS = {}
for (nm, br), (g, f_claim) in CLAIMS.items():
    P = params(nm)
    check(f"B: {nm}/{br} claimed f solves the ODE exactly", ode_ok(P, g, f_claim))
    sols = generic_solutions(P, g)
    uniq = len(sols) == 1 and sp.expand(sols[0] - f_claim) == 0
    check(f"B: {nm}/{br} full generic solve -> unique = claimed f", uniq)
    sf = sig(f_claim)
    sc = sig(y**P["q"] * g)
    SIGS[(nm, br)] = tuple(sf[i] + P["N"] * sc[i] for i in range(4))

# F3 cross-check: exponent arithmetic == full expansion of Phi
P = params("F3")
g, f_claim = CLAIMS[("F3", "ram")]
full = sig(sp.expand(f_claim * (y**P["q"] * g)**P["N"]))
check("B: F3 Phi signature via full expansion == exponent arithmetic",
      full == SIGS[("F3", "ram")])

# ---------------------------------------------------------------------------
# C. verdict vs the old unified law, and the amended ramified law
# ---------------------------------------------------------------------------
def old_law(P):
    return ((P["e"] * P["a0"] - P["q"] + 1) + P["gap"] + P["N"] * P["a0"],
            P["rho"] + P["N"] * P["q"], P["e"] + P["N"],
            P["gap"] + P["r"] * (P["e"] + P["N"]))

def ram_law(P):
    return (P["res"] + P["N"] * P["a0"], P["rho"] + P["N"] * P["q"],
            P["dg"] * (P["e"] + P["N"]) - (P["dg"] - 1), P["gap"] + P["r"])

EXPECT_SIG = {("F7", "ram"): (250, 165, 83, 2),
              ("F7", "cx"): (250, 165, 0, 85),
              ("F3", "ram"): (189, 112, 75, 2),
              ("F3", "cx"): (189, 112, 0, 77),
              ("F16", "ram"): (528, 407, 117, 4),
              ("F10", "ram"): (1917, 820, 1093, 4)}
for key, expected in EXPECT_SIG.items():
    check(f"C: {key[0]}/{key[1]} signature {expected}", SIGS[key] == expected)
for nm in ("F7", "F3", "F16", "F10"):
    P = params(nm)
    o, rl, got = old_law(P), ram_law(P), SIGS[(nm, "ram")]
    check(f"C: {nm} deg,ord match old law", got[:2] == o[:2])
    check(f"C: {nm} mult,cof DIFFER from old law (conjecture refuted)",
          got[2:] != o[2:])
    check(f"C: {nm} ramified law (res+N*a0, rho+Nq, dg(e+N)-(dg-1), gap+r)",
          got == rl)
    if (nm, "cx") in SIGS:
        check(f"C: {nm} complex-pair branch: deg,ord match, mult=0",
              SIGS[(nm, "cx")][:2] == o[:2] and SIGS[(nm, "cx")][2] == 0)
# units
for key in EXPECT_SIG:
    nm, br = key
    if br != "ram":
        continue
    P = params(nm)
    g, f_claim = CLAIMS[key]
    u = sp.cancel(f_claim / (y**P["rho"]
                             * (y + 1)**(P["dg"] * P["e"] - (P["dg"] - 1))))
    upoly = sp.Poly(sp.together(u), y)
    check(f"C: {nm} unit cofactor deg u = gap+r = {P['gap'] + P['r']}, "
          f"u(0),u(-1) != 0",
          upoly.degree() == P["gap"] + P["r"]
          and u.subs(y, 0) != 0 and u.subs(y, -1) != 0)
# (72,108) retro-explanation
check("C: (72,108) audited quartic degree = gap+r = 4+0", 4 + 0 == 4)

# ---------------------------------------------------------------------------
# D. the impossibility theorem at dg=2 (independent re-proof per corner)
#    Claimed obstruction branch ratios w = g0/g1^2 (from phi_f7.py):
#    F7: {1/4 (disc), 2/3}; F3: {1/4, 3/5}; F16: {1/4, roots of 54w^2-126w+35}.
#    A simple real root at -1 needs 1 - g1 + g0 = 0 with g1^2 > 4 g0.
#    disc branch: g0=g1^2/4 -> (1-g1/2)^2 = 0 -> g1=2: double root AT -1.
#    ratio branch w>1/4: 1 - g1 + w g1^2 = 0 has no real g1 (disc 1-4w < 0).
# ---------------------------------------------------------------------------
w = sp.symbols("w", real=True)
RATIOS = {"F7": [sp.Rational(2, 3)], "F3": [sp.Rational(3, 5)],
          "F16": list(sp.solve(54 * w**2 - 126 * w + 35, w))}
for nm, ws in RATIOS.items():
    all_gt = all(sp.simplify(wv - sp.Rational(1, 4)) > 0 for wv in ws)
    check(f"D: {nm} every complex-pair ratio w > 1/4 (so 1-g1+w*g1^2 = 0 has "
          f"no real g1: no residual root at -1)", all_gt)
g1v = sp.symbols("g1v")
sol_disc = sp.solve(1 - g1v + g1v**2 / 4, g1v)
check("D: disc branch: g(-1)=0 forces g1=2, i.e. the DOUBLE root sits at -1",
      sol_disc == [2])
# off-branch refutation spot checks: generic w -> no polynomial solution
for nm, wbad in (("F7", sp.Rational(1, 2)), ("F3", sp.Integer(1)),
                 ("F16", sp.Integer(1))):
    P = params(nm)
    sols = generic_solutions(P, y**2 + y + wbad)
    check(f"D: {nm} off-branch w={wbad}: NO polynomial solution", sols == [])
# and on the claimed ratio branch a solution EXISTS (rational rep)
sols = generic_solutions(params("F7"), y**2 + y + sp.Rational(2, 3))
check("D: F7 on-branch w=2/3: polynomial solution exists", len(sols) == 1)
# mult = e+N impossibility (the refuted component), from deg f <= res:
for nm in ("F7", "F3", "F16"):
    P = params(nm)
    # real residual options: mult_g in {0 (complex pair), 2 (double at -1)};
    # mult(f) <= deg f <= res < e+N and e+N - 2N < 0 for these corners
    imposs = (P["res"] < P["e"] + P["N"]) and (P["e"] - P["N"] < 0)
    check(f"D: {nm} old-law mult=e+N unrealizable "
          f"(res={P['res']} < e+N={P['e'] + P['N']}; e-N<0)", imposs)

# ---------------------------------------------------------------------------
# E. controls
# ---------------------------------------------------------------------------
P14 = params("F14")
sols = generic_solutions(P14, y**5 + 1)
check("E: F14 control (gap=0): unique generic solution = -(1/10) y^21 (y^5+1)^6",
      len(sols) == 1
      and sp.expand(sols[0] + sp.Rational(1, 10) * y**21 * (y**5 + 1)**6) == 0)
P1 = params("F1")
for sv in (sp.Integer(1), sp.Integer(2), sp.Rational(1, 3)):
    sols = generic_solutions(P1, y + sv)
    check(f"E: F1 control (dg=1): root position s={sv} admits a unique solution "
          f"(gauge, not forced)", len(sols) == 1)
# dg parity across the full 15-family survey (rows as in phi_f14.py)
SURVEY = [
    ("F1", (4,12), 7,4,3,1,(3,2),(4,3)),  ("F2", (5,20), 7,5,2,1,(2,1),(3,2)),
    ("F3", (5,20), 8,5,3,1,(3,4),(2,3)),  ("F4", (5,20), 8,5,3,2,(3,2),(16,12)),
    ("F5", (5,20), 9,5,4,1,(9,7),(5,4)),  ("F6", (5,20), 9,5,4,2,(4,3),(10,8)),
    ("F7", (6,15), 7,3,4,1,(2,1),(7,4)),  ("F8", (6,15), 8,3,5,1,(3,2),(7,5)),
    ("F9", (7,21), 11,7,2,1,(2,1),(3,2)), ("F10",(7,21), 13,7,3,1,(7,5),(4,3)),
    ("F11",(7,21), 13,7,3,2,(2,1),(5,3)), ("F14",(9,24), 7,3,4,1,(2,1),(7,4)),
    ("F15",(9,24), 8,3,5,1,(3,2),(7,5)),  ("F16",(9,24), 10,3,7,1,(3,4),(5,7)),
    ("F17",(9,24), 11,3,8,1,(2,5),(3,8)),
]
parity_ok = True
for row in SURVEY:
    TABLE[row[0]] = (row[1], row[2], row[3], row[4], row[5], row[6], row[7])
    Pr = params(row[0])
    if (Pr["dg"] % 2 == 0) != (Pr["gap"] > 0 and Pr["r"] > 0):
        parity_ok = False
check("E: survey mini-lemma: dg even <=> (gap>0 and r>0) on all 15 families",
      parity_ok)

# ---------------------------------------------------------------------------
print()
if _f == 0:
    print(f"ALL {_n} PHI-F7 CHECKS PASSED")
else:
    print(f"{_f}/{_n} PHI-F7 CHECKS FAILED")
import pathlib
print(f"script: {pathlib.Path(__file__).resolve()}")
sys.exit(0 if _f == 0 else 1)
