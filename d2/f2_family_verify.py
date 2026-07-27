#!/usr/bin/env python3
"""F2-family symbolic forcing theorem (first family-level closed form).

*** REPAIRED 2026-07-26 (PASSPORT_75_125_REPAIR.md).  The (5,20) corner has
t = ceil(b0/a0) = 4, kappa = 2, C = y (a MONOMIAL, deg C = ord C = 1) -- not
t = 5, kappa = 3, C = y^2(y^3+1).  GGV5's final chain corner (7\\5,2) is chart
data only on the retraction shape b0 = l(a0-1), which (5,20) fails.  The
family-level THEOREM survives in shape (a closed form for every a, a block
recurrence, a linear window-denominator law) and every constant in it changes. ***

For the F2 family (a = j+2, b = 2j+3 = 2a-1; corner t=4, kappa=2, C = y,
so there is no residual g at all) the forcing ODE

    a { 4 c f' - (4a-1) c' f } = c^a

is solved for EVERY a by the closed form

    f_a = (1/a) y^a,

the tower length is N_a = (3a-2)(4a-1) = 12a^2 - 11a + 2, hence

    Phi_a = (1/a) y^(12a^2-10a+2),   a MONOMIAL,

with the block recurrence  Phi_{a+1} = (a/(a+1)) C^(24a+2) Phi_a, C = y.
(Conditional on the standard-chart reduction and on the INFERRED rule
l = ceil(b0/a0), as for every derived corner.)

Also checks the WINDOW-DENOMINATOR law: W_step = ord_y(Phi_a)/M_a is
non-integral for ALL a >= 2 with reduced denominator exactly
q_window = 12a-7 = M_a (a=2: 30/17; a=3: 80/29) — the integral W_step = 12 at
(72,108) is a friendly coincidence of that corner, not generic a=2 behavior.
Note q_window == M_a exactly, which is what makes the ord-side carry
obstruction of weight_lemma_75_125.py total at every rung.

In the superseded model these read N_a = 15a^2-13a+2 and q_window = 5a-3
(7 at a=2, 12 at a=3).  Both are wrong; the periods are 17 and 29, and both are
PRIME, so any argument resting on the divisor structure of 12 is void.

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
T, KAPPA = 4, 2
c = y                                      # C = y, a MONOMIAL  [REPAIRED]

# A0. the corner inputs come from the retraction guard, both directions
import polygon_reduction as pr             # noqa: E402
_cd = pr.corner_chart_data(5, 20, l_final=5, b_final=2, who="f2_family_verify")
ok("A0: guard gives (t,kappa,deg C,ord C) = (4,2,1,1), C a MONOMIAL, no retraction",
   (_cd["t"], _cd["kappa"], _cd["deg_C"], _cd["ord_C"]) == (4, 2, 1, 1)
   and _cd["monomial"] and not _cd["retraction"])
try:
    pr.final_corner_dictionary(5, 20, 5, 2)
    _raised = False
except pr.FinalCornerDictionaryError:
    _raised = True
ok("A0: the superseded (t,q)=(l_final,b_final)=(5,2) RAISES, and the dictionary "
   "still RETURNS at (8,28) -- guard checked in both directions",
   _raised and pr.final_corner_dictionary(8, 28, 4, 7) == (4, 7))

# A. the closed form solves the ODE symbolically in a
f = sp.Rational(1, 1) / a * y**a
lhs = a * (T * c * sp.diff(f, y) - (T * a - 1) * sp.diff(c, y) * f)
ok("A: ODE residual vanishes symbolically in a  (a{4 c f' - (4a-1) c' f} = c^a)",
   sp.simplify(lhs - c**a) == 0)
ok("A: collapse identity a*t - coef = 4a - (4a-1) = 1, so A = 1/a",
   sp.expand(T * a - (T * a - 1)) == 1)
# and the superseded closed form must FAIL the repaired ODE -- discriminating
_fold = -sp.Rational(1, 3) / a * y**(2 * a - 1) * (y**3 + 1)**a
ok("A: the superseded f_a = -(1/(3a)) y^(2a-1)(y^3+1)^a does NOT solve it",
   sp.simplify(a * (T * c * sp.diff(_fold, y)
                    - (T * a - 1) * sp.diff(c, y) * _fold) - c**a) != 0)

# B. Phi_a signature formulas against both landed points
def phi(aa):
    return sp.Rational(1, aa) * y**(12 * aa**2 - 10 * aa + 2)

landed = {2: (28, 30, 30, 0, 0), 3: (77, 80, 80, 0, 0)}
for aa, (N_l, deg_l, ord_l, mult_l, cof_l) in landed.items():
    Na = 12 * aa**2 - 11 * aa + 2
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
deg_form = 12 * a**2 - 10 * a + 2
ord_form = a + (12 * a**2 - 11 * a + 2)          # deg f + N * ord C  (ord C = 1)
ok("C: ord_y(Phi_a) = deg_y(Phi_a) = 12a^2-10a+2 symbolically (Phi_a is a monomial)",
   sp.expand(ord_form - deg_form) == 0)
ok("C: N_a factors as (3a-2)(4a-1) = 12a^2-11a+2  [was (3a-2)(5a-1); the 5a-1 "
   "became 4a-1 because t went 5 -> 4]",
   sp.expand((3 * a - 2) * (4 * a - 1) - (12 * a**2 - 11 * a + 2)) == 0)

# D. block recurrence Phi_{a+1} = (a/(a+1)) C^(24a+2) Phi_a at a = 2..4
for aa in (2, 3, 4):
    ok("D: recurrence at a=%d (C^(24a+2) = C^%d)" % (aa, 24 * aa + 2),
       sp.simplify(phi(aa + 1) - sp.Rational(aa, aa + 1) * c**(24 * aa + 2) * phi(aa)) == 0)

# E. window-denominator law (corrects the "a>=3 boundary")
for aa in (2, 3, 4, 5):
    bb = 2 * aa - 1
    M = bb * T + (aa * T - KAPPA - 1)
    W = sp.Rational(12 * aa**2 - 10 * aa + 2, M)
    ok("E: a=%d W_step=%s non-integral, denom = 12a-7 = %d = M" % (aa, W, 12 * aa - 7),
       (not W.is_integer) and sp.denom(W) == 12 * aa - 7 == M)
ok("E: gcd(ord_y(Phi_a), M_a) = 1 for a=2..12, so q_window = M_a exactly and the "
   "window fraction never reduces",
   all(sp.gcd(12 * aa**2 - 10 * aa + 2, 12 * aa - 7) == 1 for aa in range(2, 13)))
ok("E: consecutive periods are coprime (gcd(12a-7, 12(a+1)-7) = gcd(m, m+12) | 12 "
   "and 12a-7 is coprime to 12), a=2..8",
   all(sp.gcd(12 * aa - 7, 12 * (aa + 1) - 7) == 1 for aa in range(2, 9)))
# ((72,108)'s integral W_step = 12 is that corner's own coincidence — recorded
# in G_SYSTEM_75_125.md; it is not an instance of the F2 law above.)

print()
if FAILS:
    print("FAILURES:", len(FAILS))
    sys.exit(1)
print("ALL %d F2-FAMILY CHECKS PASSED" % N_OK)
