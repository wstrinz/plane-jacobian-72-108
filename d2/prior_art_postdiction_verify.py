#!/usr/bin/env python3
"""Postdiction check: the corner law reproduces the ONE closed-form last
element in the published literature with zero fitting freedom.

GGHV22 (arXiv:2204.14178) sec.4, killing the (66,99)/(9,24)-corner cases,
derives (their L1571-1597, notation adapted):

    6*C3*f1' - 10*C3'*f1 = C3^2,   C3 = y^8*(y+1),

states the solution is unique, and prints

    f1 = -(1/910) * y^9 * (y+1)^2 * (243y^4 - 81y^3 + 54y^2 - 42y + 35).

Checks here (all exact, sympy):
  A. their ODE is the parametric family a{t c f' - [t(b-a)+kappa+1] c' f}
     = c^(b-a+1) at (a,b,t,kappa) = (2,3,3,1) -- i.e. kappa = t-2 holds on
     the published instance;
  B. the ODE has a UNIQUE polynomial solution of degree <= 15 and it equals
     their printed f1 coefficient-for-coefficient;
  C. the f-level corner law postdicts its full signature at the (66,99)/F17
     corner (a0,q,t)=(9,8,3), (a,b)=(2,3): e=2, r=0, gap=4, dg=1 ->
     ord_y = (e-1)q+1 = 9, mult_{y+1} = e = 2 (unramified, dg odd),
     deg = (e*a0-q+1)+gap = 15, unit cofactor of degree gap+r = 4;
  D. the quartic cofactor is separable and coprime to y(y+1) (the property
     GGHV22's kill actually uses), mirroring the audited (72,108) quartic.

Zero fitting freedom: the law's parameters come from corner data alone.
Source: PRIOR_ART.md; the law: PHI_CORNER4.md / PHI_F14.md / PHI_F7.md.
"""
import sys
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

y = sp.symbols("y")

# --- A. parametric-family membership, kappa = t-2 ---------------------------
a, b, t, kappa = 2, 3, 3, 1
ok("A: kappa = t-2 at the published instance", kappa == t - 2)
c = sp.expand(y**8 * (y + 1))          # C3
f = sp.Function("f")(y)
family_lhs = a * (t * c * sp.diff(f, y) - (t * (b - a) + kappa + 1) * sp.diff(c, y) * f)
family_rhs = c ** (b - a + 1)
published_lhs = 6 * c * sp.diff(f, y) - 10 * sp.diff(c, y) * f
ok("A: family ODE at (2,3,3,1) is exactly 6*C3*f' - 10*C3'*f = C3^2",
   sp.simplify(family_lhs - published_lhs) == 0
   and sp.simplify(family_rhs - c**2) == 0)

# --- B. unique polynomial solution of degree <= 15 == printed f1 ------------
coeffs = sp.symbols("a0:16")
fp = sum(coeffs[i] * y**i for i in range(16))
eq = sp.expand(6 * c * sp.diff(fp, y) - 10 * sp.diff(c, y) * fp - c**2)
sols = sp.solve([eq.coeff(y, k) for k in range(sp.degree(eq, y) + 1)],
                list(coeffs), dict=True)
ok("B: linear system has exactly one solution (uniqueness, deg <= 15)",
   len(sols) == 1)
f1_derived = sp.expand(fp.subs(sols[0]))
f1_published = sp.expand(sp.Rational(-1, 910) * y**9 * (y + 1)**2
                         * (243*y**4 - 81*y**3 + 54*y**2 - 42*y + 35))
ok("B: derived solution equals GGHV22's printed f1 exactly",
   sp.expand(f1_derived - f1_published) == 0)

# --- C. f-level corner-law postdiction --------------------------------------
a0, q = 9, 8
e = b - a + 1
r = a0 - q - 1
gap = (q - 1) - sp.Rational(a0, t)
dg = a0 - q
ok("C: corner data gives e=2, r=0, gap=4, dg=1",
   (e, r, gap, dg) == (2, 0, 4, 1))
p = sp.Poly(f1_published, y)
ord_y = min(m[0] for m in p.monoms())
mult = sp.Poly(f1_published, y).as_expr()
mult_t = 0
rem = f1_published
while sp.simplify(rem.subs(y, -1)) == 0:
    rem = sp.cancel(rem / (y + 1))
    rem = sp.expand(rem)
    mult_t += 1
ok("C: ord_y = (e-1)q+1 = 9", ord_y == (e - 1) * q + 1 == 9)
ok("C: mult_{y+1} = e = 2 (unramified: dg odd)", mult_t == e == 2)
ok("C: deg = (e*a0-q+1)+gap = 15", p.degree() == (e * a0 - q + 1) + gap == 15)

# --- D. the unit cofactor: degree gap+r, separable, coprime to y(y+1) -------
cof = sp.expand(sp.cancel(f1_published / (y**9 * (y + 1)**2)))
cq = sp.Poly(cof, y)
ok("D: unit cofactor degree = gap + r = 4", cq.degree() == gap + r == 4)
ok("D: cofactor separable (disc != 0)", sp.discriminant(cof, y) != 0)
ok("D: cofactor coprime to y(y+1) (unit at both places)",
   cof.subs(y, 0) != 0 and cof.subs(y, -1) != 0)
ok("D: cofactor is a rational multiple of 243y^4-81y^3+54y^2-42y+35",
   sp.simplify(cof / (243*y**4 - 81*y**3 + 54*y**2 - 42*y + 35)).is_Rational)

print()
if FAILS:
    print("FAILURES:", len(FAILS))
    sys.exit(1)
print("ALL %d PRIOR-ART POSTDICTION CHECKS PASSED" % N_OK)
