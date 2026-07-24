#!/usr/bin/env python3
"""F2-family symbolic forcing theorem (first family-level closed form).

For the F2 family (a = j+2, b = 2j+3 = 2a-1; corner t=5, kappa=3,
c = y^2(y^3+1), g = y^3+1) the forcing ODE

    a { 5 c f' - (5a-1) c' f } = c^a

is solved for EVERY a by the closed form

    f_a = -(1/(3a)) y^(2a-1) g^a,

the tower length is N_a = 15a^2 - 13a + 2, hence

    Phi_a = -(1/(3a)) y^(30a^2-24a+3) (y^3+1)^(15a^2-12a+2),

with the block recurrence  Phi_{a+1} = (a/(a+1)) C^(30a+3) Phi_a,
C = y^2(y^3+1).  (Conditional on the standard-chart reduction, as for
every derived corner; source: GPT-Pro review 6, verified here.)

Also checks the WINDOW-DENOMINATOR law that corrects G_SYSTEM_75_125's
original "a>=3 boundary" claim: W_step = ord_y(Phi_a)/M_a is non-integral
for ALL a >= 2 with reduced denominator exactly q_window = 5a-3
(a=2: 25/7; a=3: 67/12) — the integral W_step = 12 at (72,108) is a
friendly coincidence of that corner, not generic a=2 behavior.

Checks A-E, exact sympy; --quiet; exit 0 iff all pass.
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
a = sp.symbols("a", positive=True, integer=True)
g = y**3 + 1
c = y**2 * g

# A. the closed form solves the ODE symbolically in a
f = -sp.Rational(1, 3) / a * y**(2 * a - 1) * g**a
lhs = a * (5 * c * sp.diff(f, y) - (5 * a - 1) * sp.diff(c, y) * f)
ok("A: ODE residual vanishes symbolically in a",
   sp.simplify(lhs - c**a) == 0)
ok("A: collapse identity y*g' - 3g = -3", sp.expand(y * sp.diff(g, y) - 3 * g) == -3)

# B. Phi_a signature formulas against both landed points
def phi(aa):
    return -sp.Rational(1, 3 * aa) * y**(2 * aa - 1) * g**aa * c**(15 * aa**2 - 13 * aa + 2)

landed = {2: (36, 189, 75, 38, 76), 3: (98, 504, 201, 101, 202)}
for aa, (N_l, deg_l, ord_l, mult_l, cof_l) in landed.items():
    Na = 15 * aa**2 - 13 * aa + 2
    ok("B: N_%d = %d" % (aa, N_l), Na == N_l)
    p = sp.Poly(sp.expand(phi(aa)), y)
    ordy = min(m[0] for m in p.monoms())
    mult, rem = 0, p.as_expr()
    while sp.simplify(rem.subs(y, -1)) == 0:
        rem = sp.expand(sp.cancel(rem / (y + 1)))
        mult += 1
    ok("B: Phi_%d signature (%d,%d,%d,%d)" % (aa, deg_l, ord_l, mult_l, cof_l),
       (p.degree(), ordy, mult, p.degree() - ordy - mult) == (deg_l, ord_l, mult_l, cof_l))

# C. closed-form exponents match the general formulas symbolically
deg_form = 30 * a**2 - 24 * a + 3
ord_form = 2 * a - 1 + 2 * (15 * a**2 - 13 * a + 2)
ok("C: ord_y(Phi_a) = 30a^2-24a+3 symbolically",
   sp.simplify(ord_form - deg_form) == 0 or sp.expand(ord_form) == sp.expand(30*a**2 - 24*a + 3))

# D. block recurrence Phi_{a+1} = (a/(a+1)) C^(30a+3) Phi_a at a = 2..4
for aa in (2, 3, 4):
    ok("D: recurrence at a=%d" % aa,
       sp.simplify(phi(aa + 1) - sp.Rational(aa, aa + 1) * c**(30 * aa + 3) * phi(aa)) == 0)

# E. window-denominator law (corrects the "a>=3 boundary")
for aa in (2, 3, 4, 5):
    bb = 2 * aa - 1
    M = bb * 5 + (aa * 5 - 4)
    W = sp.Rational(30 * aa**2 - 24 * aa + 3, M)
    ok("E: a=%d W_step=%s non-integral, denom = 5a-3 = %d" % (aa, W, 5 * aa - 3),
       (not W.is_integer) and sp.denom(W) == 5 * aa - 3 and sp.gcd(2 * aa - 1, 5 * aa - 3) == 1)
# ((72,108)'s integral W_step = 12 is that corner's own coincidence — recorded
# in G_SYSTEM_75_125.md; it is not an instance of the F2 law above.)

print()
if FAILS:
    print("FAILURES:", len(FAILS))
    sys.exit(1)
print("ALL %d F2-FAMILY CHECKS PASSED" % N_OK)
