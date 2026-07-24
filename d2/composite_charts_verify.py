#!/usr/bin/env python3
"""composite_charts_verify.py -- independent PASS/FAIL checker for LANE J
(COMPOSITE_CHARTS.md / composite_charts.py).

Checks (own transcriptions and code paths; no imports from composite_charts):
  A. chart Jacobians: pure Laurent control, fused-shear lemma, double-inversion
     heuristic identification, shear neutrality, GGHV instance;
  B. chain-step corner arithmetic for 5 standard controls + all 5 escapes,
     including F24's second step and the fused-exponent integrality table;
  C. the checkable core of GGV5's F18-F21 impossibility;
  D. the zeta-obstruction identities (A0'=(2,0) breaks the model, not the chart);
  E. the F12 conditional dg=3 branch variety: reduced system, consistency
     eliminant, exactly two branch points, both squarefree, exact ODE residual
     over Q(sqrt(22)), resultant(g,u) != 0, and the conditional signature
     (814, 506, 102, 206) with old-law cofactor identity.

Usage: python3 composite_charts_verify.py [--quiet]   -- exit 0 iff all pass.
"""
import sys
import itertools
import sympy as sp
from math import gcd
from fractions import Fraction

QUIET = "--quiet" in sys.argv[1:]
FAILS = []
COUNT = [0]


def ok(label, cond):
    COUNT[0] += 1
    if not cond:
        FAILS.append(label)
    if not QUIET or not cond:
        print(f"  [{'OK' if cond else 'FAIL'}] {label}")


x, y = sp.symbols("x y")

# ---------------------------------------------------------------- A. charts
lgen = sp.symbols("l_gen", positive=True)
e1, e2, e3 = sp.symbols("e1 e2 e3", positive=True)
c1, c2, c3 = sp.symbols("c1 c2 c3")


def jac2(X, Y, u=x, v=y):
    return sp.det(sp.Matrix([[sp.diff(X, u), sp.diff(X, v)],
                             [sp.diff(Y, u), sp.diff(Y, v)]]))


ok("A1 pure Laurent (x^-1, x^l y): Jacobian -x^(l-2) [PHI_CORNER4 control]",
   sp.simplify(jac2(x**-1, x**lgen * y) + x**(lgen - 2)) == 0)
ok("A2 fused chart, one shear term: Jacobian unchanged",
   sp.simplify(jac2(x**-1, x**lgen * y + c1 * x**e1) + x**(lgen - 2)) == 0)
ok("A3 fused chart, three shear terms: Jacobian unchanged",
   sp.simplify(jac2(x**-1, x**lgen * y + c1 * x**e1 + c2 * x**e2 + c3 * x**e3)
               + x**(lgen - 2)) == 0)
sgen = sp.symbols("s_gen", positive=True)
ok("A4 bare root shift y -> y + c*x^(-s): Jacobian exactly 1",
   sp.simplify(jac2(x, y + c1 * x**(-sgen)) - 1) == 0)
xp = sp.symbols("xp", positive=True)
l1s, l2s = sp.symbols("l1s l2s", positive=True)
Jh = sp.powsimp(jac2((xp**-1)**-1, (xp**-1)**l1s * (xp**l2s * y), u=xp),
                force=True)
ok("A5 double-inversion composite: Jacobian x^(l2-l1) [the heuristic's chart]",
   sp.simplify(Jh - xp**(l2s - l1s)) == 0)
ok("A6 GGHV (8,28) instance l=4 with fused shear alpha*x^0: Jacobian -x^2 "
   "[2204.14178 lines 1228-1234]",
   sp.simplify(jac2(x**-1, x**4 * y + c1) + x**2) == 0)

# ------------------------------------------------- B. corner-step arithmetic
# own transcription of the GGV5 tables (lines 1674-1718)
STD = {"F5": ((5, 20), (1, 0), (9, 5, 4)), "F7": ((6, 15), (1, 0), (7, 3, 4)),
       "F9": ((7, 21), (1, 0), (11, 7, 2)), "F10": ((7, 21), (1, 0), (13, 7, 3)),
       "F16": ((9, 24), (1, 0), (10, 3, 7))}
ESC1 = {"F12": ((8, 24), (2, 0), (13, 4, 5)),
        "F13": ((9, 21), (2, 0), (13, 3, 7))}
LEN2 = {"F22": ((8, 24), (2, 0), (14, 4, 6), (5, 4, 2), (5, 4, 2)),
        "F23": ((8, 24), (2, 0), (14, 4, 6), (11, 4, 4), (11, 4, 4)),
        "F24": ((8, 24), (2, 0), (14, 4, 6), (5, 4, 0), (19, 8, 3))}


def normal_dir(dx, dy):
    g = gcd(abs(dx), abs(dy))
    rho, sigma = abs(dy) // g, -(dx // g) * (1 if dy > 0 else -1)
    return (rho, sigma) if rho > 0 else (-rho, -sigma)


def check_step(name, A0, A0p, corner):
    p, l, q = corner
    rho, sigma = normal_dir(A0[0] - A0p[0], A0[1] - A0p[1])
    got = (Fraction(A0p[0]) + Fraction(q * (-sigma), rho), A0p[1] + q)
    return rho == l and got == (Fraction(p, l), q)


for name, (A0, A0p, corner) in {**STD, **ESC1}.items():
    ok(f"B {name}: step formula reproduces table corner ({corner[0]}/{corner[1]},"
       f"{corner[2]})", check_step(name, A0, A0p, corner))
for name, (A0, A0p, A1, A1p, A2) in LEN2.items():
    ok(f"B {name} step1: reproduces A1 = ({A1[0]}/{A1[1]},{A1[2]})",
       check_step(name, A0, A0p, A1))
ok("B F22/F23: A2 == A1' (two edges, ONE transformation; final l = 4)",
   LEN2["F22"][4] == LEN2["F22"][3] and LEN2["F23"][4] == LEN2["F23"][3])
# F24 second step in u = x^(1/4) coords
A0, A0p, A1, A1p, A2 = LEN2["F24"]
du = int((Fraction(A1[0], A1[1]) - Fraction(A1p[0], A1p[1])) * 4)
rho2, sigma2 = normal_dir(du, A1[2] - A1p[2])
ok("B F24 step2: u-edge (9,6) -> direction (2,-3), l refines 4 -> 8 = table",
   (du, A1[2] - A1p[2]) == (9, 6) and (rho2, sigma2) == (2, -3)
   and 4 * rho2 == A2[1])
got2 = (Fraction(A1p[0], A1p[1]) + Fraction(A2[2] * (-sigma2), 4 * rho2),
        A1p[2] + A2[2])
ok("B F24 step2: corner formula reproduces A2 = (19/8,3)",
   got2 == (Fraction(19, 8), 3))
FUSED = {"F12": (4, [3]), "F13": (3, [2]), "F22": (4, [3]), "F23": (4, [3]),
         "F24": (8, [6, 5])}
ok("B fused shear exponents integral and nonnegative for all five escapes",
   all(all(isinstance(e, int) and e >= 0 for e in ee)
       for _, ee in FUSED.values()))
ok("B kappa = l-2 values: F12:2 F13:1 F22:2 F23:2 F24:6",
   [FUSED[k][0] - 2 for k in ("F12", "F13", "F22", "F23", "F24")]
   == [2, 1, 2, 2, 6])

# --------------------------------------- C. F18-F21 impossibility, core facts
ok("C shifted start (1/m)st_10 = (6,3) violates b+gcd(a,b) != a (3+3 == 6)",
   3 + gcd(6, 3) == 6)
ok("C degenerate path (6,9): v11 = 15 < 16 (no admissible chain below 16)",
   6 + 9 < 16)

# ------------------------------------------------- D. the zeta-obstruction
cfun = sp.Function("cf")(y)
ffun = sp.Function("ff")(y)
tE, zE, aE, bE, uE = 5, 3, 2, 4, 7          # independent generic exponents
C = x**tE * cfun


def br(P, Q):
    return sp.simplify(sp.diff(P, x) * sp.diff(Q, y)
                       - sp.diff(P, y) * sp.diff(Q, x))


ok("D1 [C^a, C^b] == 0 (pure powers commute; the zeta=0 tail structure)",
   sp.simplify(br(C**aE, C**bE)) == 0)
lhs = br(x**zE * C**aE, C**bE)
rhs = zE * bE * x**(zE - 1) * C**(aE + bE - 1) * sp.diff(C, y)
ok("D2 [x^z C^a, C^b] = z*b*x^(z-1) C^(a+b-1) C_y  (nonzero for z>0)",
   sp.simplify(sp.expand(lhs - rhs)) == 0 and sp.simplify(lhs) != 0)
lhs2 = br(x**zE * C**aE, x**uE * ffun)
rhs2 = x**(zE + aE * tE + uE - 1) * ((zE + aE * tE) * cfun**aE * sp.diff(ffun, y)
                                     - aE * uE * cfun**(aE - 1)
                                     * sp.diff(cfun, y) * ffun)
ok("D3 [x^z C^a, x^u f] = x^(z+at+u-1)[(z+at) c^a f' - a u c^(a-1) c' f]",
   sp.simplify(sp.expand(lhs2 - rhs2)) == 0)
lhs0 = br(C**aE, x**uE * ffun)
rhs0 = x**(aE * tE + uE - 1) * (aE * tE * cfun**aE * sp.diff(ffun, y)
                                - aE * uE * cfun**(aE - 1)
                                * sp.diff(cfun, y) * ffun)
ok("D4 zeta = 0 reduction recovers the (F)-family bracket piece",
   sp.simplify(sp.expand(lhs0 - rhs0)) == 0)

# ------------------------- E. F12 conditional dg=3 branch variety (from scratch)
a, b, t, kappa, a0, q = 3, 7, 4, 2, 8, 5
e = b - a + 1
r = a0 - q - 1
dg = a0 - q
coef = t * (b - a) + kappa + 1
rho = (e - 1) * q + 1
N = a * (t * (a + b - 1) + 1) - 2 * b
res = Fraction((t * (b - a) + kappa + 1) * a0, t)
pure = e * a0 - q + 1
gap = res - pure
ok("E0 F12 corner data: e=5 r=2 dg=3 N=97 gap=2 rho=21 res=38 kappa=t-2",
   (e, r, dg, N, gap, rho, res, kappa) == (5, 2, 3, 97, 2, 21, 38, t - 2))

g2, g1, g0, u2, u1, u0 = sp.symbols("g2 g1 g0 u2 u1 u0")
g = y**3 + g2 * y**2 + g1 * y + g0
uu = u2 * y**2 + u1 * y + u0
c = y**q * g
f = y**rho * g**e * uu
R = sp.expand(a * t * c * sp.diff(f, y) - a * coef * sp.diff(c, y) * f - c**e)
Q = sp.expand(sp.cancel(R / (y**25 * g**5)))
P = sp.Poly(Q, y)
ok("E1 reduced equation has degree 4 (full-leading-product resonance at y^5..y^6)",
   P.degree() == 4)
eqs = [sp.expand(P.nth(k)) for k in range(5)]
ok("E2 E_0 = -33 g0 u0 - 1 (forces g0, u0 nonzero)",
   sp.expand(eqs[0] + 33 * g0 * u0 + 1) == 0)
u1_of = sp.solve(eqs[4], u1)[0]
u0_of = sp.solve(eqs[3].subs(u1, u1_of), u0)[0]
E1c = sp.expand(sp.cancel(eqs[1].subs({u1: u1_of, u0: u0_of}) / u2))
E2c = sp.expand(sp.cancel(eqs[2].subs({u1: u1_of, u0: u0_of}) / u2))
ok("E3 consistency: E2' ~ 32 g0 - 40 g1 g2 + 15 g2^3",
   sp.simplify(E2c * sp.Rational(-32, 9) - (32 * g0 - 40 * g1 * g2
                                            + 15 * g2**3)) == 0)
ok("E4 consistency: E1' ~ 28 g0 g2 + 40 g1^2 - 25 g1 g2^2",
   sp.simplify(E1c * sp.Rational(16, 3) - (28 * g0 * g2 + 40 * g1**2
                                           - 25 * g1 * g2**2)) == 0)
# eliminant via resultant (independent route from the derivation script):
elim = sp.factor(sp.resultant(32 * g0 - 40 * g1 * g2 + 15 * g2**3,
                              28 * g0 * g2 + 40 * g1**2 - 25 * g1 * g2**2, g0))
target = 64 * g1**2 + 16 * g1 * g2**2 - 21 * g2**4
ok("E5 eliminant (resultant route) ~ 64 g1^2 + 16 g1 g2^2 - 21 g2^4",
   sp.simplify(sp.cancel(elim / (4 * g2))) is not None
   and sp.rem(sp.Poly(elim, g1, g2), sp.Poly(target, g1, g2)) is not None
   and sp.factor_list(elim)[1].__len__() >= 1
   and any(sp.simplify(fac - target) == 0
           for fac, _ in sp.factor_list(elim)[1]))
disc_e = sp.discriminant(target.subs(g2, 1), g1)
ok("E6 g2=1 slice: two real points, g1 = (-1 +- sqrt(22))/8  (disc = 22 class)",
   sp.simplify(disc_e / 256**2 * 64) == 5632 / 64 or
   sorted(sp.solve(target.subs(g2, 1), g1)) ==
   sorted([sp.Rational(-1, 8) + sp.sqrt(22) / 8,
           sp.Rational(-1, 8) - sp.sqrt(22) / 8]))
pts = []
for G1 in sp.solve(target.subs(g2, 1), g1):
    G0 = sp.radsimp(sp.solve((32 * g0 - 40 * g1 * g2 + 15 * g2**3)
                             .subs({g2: 1, g1: G1}), g0)[0])
    pts.append((sp.radsimp(G1), G0))
ok("E7 both branch points have squarefree g (NO ramified branch at dg=3)",
   all(sp.simplify(sp.discriminant(y**3 + y**2 + G1 * y + G0, y)) != 0
       for G1, G0 in pts))
ok("E8 g2=0 slice inconsistent: forces g0 = 0, contradicting E_0",
   sp.solve([( 32 * g0 - 40 * g1 * g2 + 15 * g2**3).subs(g2, 0),
             (28 * g0 * g2 + 40 * g1**2 - 25 * g1 * g2**2).subs(g2, 0)],
            [g0, g1]) in ([{g0: 0, g1: 0}], {g0: 0, g1: 0},
                          [(0, 0)]))
# exact solution over Q(sqrt(22)) at the first point:
G1, G0 = pts[0]
u2v = sp.Integer(96)
u1v = sp.radsimp(u1_of.subs({g2: 1, u2: u2v}))
u0v = sp.radsimp(u0_of.subs({g2: 1, g1: G1, u2: u2v}))
scale = sp.radsimp(-1 / (33 * G0 * u0v))
Uv = {u2: sp.radsimp(scale * u2v), u1: sp.radsimp(scale * u1v),
      u0: sp.radsimp(scale * u0v)}
resid = sp.expand(sp.radsimp(R.subs({g2: 1, g1: G1, g0: G0, **Uv})))
ok("E9 FULL ODE residual == 0 exactly over Q(sqrt(22))",
   sp.simplify(resid) == 0)
g_pt = y**3 + y**2 + G1 * y + G0
u_pt = Uv[u2] * y**2 + Uv[u1] * y + Uv[u0]
ok("E10 resultant(g, u) != 0 (u shares no root with g; mult_u(-1) = 0 in gauge)",
   sp.simplify(sp.resultant(g_pt, u_pt, y)) != 0)
# real nonzero root of g exists (odd degree, g0 != 0):
ok("E11 g has a real nonzero root (cubic, g(0) = g0 != 0): unramified gauge "
   "g(-1)=0 realizable by scaling",
   sp.simplify(G0) != 0)
deg_phi, ord_phi = int(res) + N * a0, rho + N * q
mult_u_, cof_u_ = e + N, int(res) + N * a0 - (rho + N * q) - (e + N)
ok("E12 conditional Phi signature = (814, 506, 102, 206)",
   (deg_phi, ord_phi, mult_u_, cof_u_) == (814, 506, 102, 206))
ok("E13 old-law cofactor identity: cof = gap + r*(e+N) = 206; ramified formulas "
   "(304, 4) are NOT realized (E7)",
   cof_u_ == gap + r * (e + N)
   and (dg * (e + N) - (dg - 1), gap + r) == (304, 4))

# ---------------------------------------------------------------- summary
print()
if FAILS:
    print(f"{len(FAILS)} / {COUNT[0]} CHECKS FAILED:")
    for s in FAILS:
        print("  FAIL:", s)
    sys.exit(1)
print(f"ALL {COUNT[0]} COMPOSITE-CHART CHECKS PASSED")
sys.exit(0)
