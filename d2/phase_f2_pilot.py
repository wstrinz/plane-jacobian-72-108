#!/usr/bin/env python3
"""Phase F, work item F2 -- the divisor-reconstruction pilot (alternate regime).

The experiment (PHASE_F_PLAN.md F2).  ALT_RESIDUE_CONGRUENCES.md established that
every surviving alternate-regime state's only depth-1 obligation is a
level-0 leading cancellation of  h_0  at infinity, and that all 19 h_0 initial
forms carry an all-nonzero rational point WHEN THE LEADING COEFFICIENTS (D,X,S,E)
= (lc d2, lc d1, lc sigma, lc e) ARE TREATED AS FREE.  This script asks the
Phase-F question: for the states whose finite-place divisors force
DEFECT 0 (deg p == sum_s v_s(p)), the polynomials p in {d1,sigma,e,d2} are
p = lambda * prod_s (y - s)^{v_s(p)} EXACTLY, so BOTH the leading coefficient
(lambda) AND every sub-leading coefficient (lambda * elementary symmetric
function of the marked roots) are DETERMINED.  Substituting those determined
coefficients into the FULL level-0 tie tower (the depth-1 initial form PLUS the
deeper convolution coefficients the tie depth demands) tests whether the residue
system still has an all-nonzero solution.

Places:  S = { t : y = -1 ,  r_1..r_4 : roots of q = 2048 y^4 - 512 y^3
             + 320 y^2 - 240 y + 195 }.
Regime facts used (all audited):
  * h_0(d2,d1,sigma,e) is cascade_engine.MONOMIALS[0] (the level-0 graded
    coefficient; 26 monomials, weight-20 homogeneous with weights
    (d2,d1,sigma,e) = (2,3,4,5)).
  * At the sub1 caps deg h_0 = 60; the bottom close E^21 h_0 + u r_0 = 0
    (ALT_REGIME_INF.md (I0)) forces the top `depth` coefficients of h_0(y) to
    vanish, where `depth` is the state's L0 tie depth
    (ALT_RESIDUE_CONGRUENCES.md).
  * e = t^a * E ; a q-root r_i carries v_{r_i}(e) = b_i (branch b-pattern),
    v_t(e) = a ; a T1 state carries v_t(d1), v_{r_i}(d1) in the
    ALT_REGIME_L2.md sec.2 cones; v_t(sigma) likewise.

Verdicts are read off exact saturated Groebner bases over Q (marked root via
Q[r]/(q)): unit ideal after saturating the leading coefficients  ==>  KILL.
Nothing here is committed; READ-ONLY on every audited artifact.
"""
import sys
sys.path.insert(0, '.')
import sympy as sp
from cascade_engine import MONOMIALS

y, s, r = sp.symbols('y s r')
Q = 2048 * r**4 - 512 * r**3 + 320 * r**2 - 240 * r + 195   # the fixed quartic
QP = sp.Poly(Q, r)
TOTAL_DEG = 60          # deg h_0 in y at the sub1 caps (H_0 = 60 - 6*level, level 0)
DEPTH = 14              # top coefficients tracked; >= every pilot's tie depth

# ---------------------------------------------------------------------------
# reversed-truncated-series engine: top DEPTH coefficients of h_0(y) about
# y = infinity, computed by truncated convolution (h_0 is never expanded whole).
# ---------------------------------------------------------------------------

def rev_of_poly(poly, deg):
    """Top-DEPTH coefficients of `poly` (y-degree `deg`), hi->lo (index i = y^{deg-i})."""
    P = sp.Poly(sp.expand(poly), y)
    return [P.coeff_monomial(y**(deg - i)) if deg - i >= 0 else sp.Integer(0)
            for i in range(DEPTH)]


def ser_mul(a, b):
    out = [sp.Integer(0)] * DEPTH
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            if i + j >= DEPTH or bj == 0:
                continue
            out[i + j] += ai * bj
    return [sp.expand(c) for c in out]


def ser_pow(a, n):
    res = [sp.Integer(1)] + [sp.Integer(0)] * (DEPTH - 1)
    base = list(a)
    while n:
        if n & 1:
            res = ser_mul(res, base)
        n >>= 1
        if n:
            base = ser_mul(base, base)
    return res


def h0_top(factors, degs, drop_d1=False):
    """List C with C[j] = coefficient of y^{TOTAL_DEG-j} in h_0, j = 0..DEPTH-1.

    factors: name -> (poly, y-degree) for names in {d2,d1,sig,e}.
    degs:    (deg d2, deg d1, deg sigma, deg e) (the tropical weights).
    drop_d1: T2 states (d1 == 0): those monomials vanish.
    """
    revs = {nm: rev_of_poly(f, dg) for nm, (f, dg) in factors.items()}
    cache = {}

    def fp(nm, n):
        if n == 0:
            return [sp.Integer(1)] + [sp.Integer(0)] * (DEPTH - 1)
        if (nm, n) not in cache:
            cache[(nm, n)] = ser_pow(revs[nm], n)
        return cache[(nm, n)]

    total = [sp.Integer(0)] * DEPTH
    for (k, x, z, b), coef in MONOMIALS[0]:
        if drop_d1 and x != 0:
            continue
        dM = degs[0] * k + degs[1] * x + degs[2] * z + degs[3] * b
        shift = TOTAL_DEG - dM
        if shift < 0 or shift >= DEPTH:
            continue
        term = [sp.Integer(1)] + [sp.Integer(0)] * (DEPTH - 1)
        for cnt, nm in ((k, 'd2'), (x, 'd1'), (z, 'sig'), (b, 'e')):
            if cnt:
                term = ser_mul(term, fp(nm, cnt))
        term = [sp.Integer(int(coef)) * c for c in term]
        for i in range(DEPTH - shift):
            total[shift + i] += term[i]
    return [sp.expand(c) for c in total]


def redq(expr):
    """Reduce modulo q(r) = 0 (marked-root arithmetic in Q[r]/(q))."""
    return sp.rem(sp.Poly(sp.expand(expr), r), QP).as_expr()


def is_unit_ideal(gens, order_vars, order='grevlex'):
    G = sp.groebner(gens, *order_vars, order=order)
    return sp.Integer(1) in G.exprs, G


# ---------------------------------------------------------------------------
# reconstruction helpers.  t-place is y = -1, so t^m = (y+1)^m; a marked q-root
# r contributes (y-r)^m.  lambda is the single leading scalar of a defect-0 poly.
# ---------------------------------------------------------------------------
def poly_from_divisor(scalar, t_mult, marked_mult):
    """scalar * (y+1)^t_mult * (y-r)^marked_mult  (a defect-0 factor)."""
    return sp.expand(scalar * (y + 1)**t_mult * (y - r)**marked_mult)


# ---------------------------------------------------------------------------
# PILOT A -- a11_b3000_T1, state (deg d2,d1,sigma,e) = (0,9,12,14).
#   d1  forced defect 0: v_t=5, v_{r}=4  (X_min = 5+4 = 9 = deg d1, unique split)
#   e   forced defect 0: v_t=a=11, v_{r}=b=3, deg e = 11+3 = 14
#   sigma both-tight (v_t = 12) defect 0: sigma = S (y+1)^12
#   d2  = D (deg 0 -> constant); support 8 = {d1^4 sigma^2, sigma^5}, tie depth 14.
# ---------------------------------------------------------------------------
def pilot_A():
    D, X, S, E, w = sp.symbols('D X S E w')
    facs = {
        'd2':  (D, 0),
        'd1':  (poly_from_divisor(X, 5, 4), 9),
        'sig': (poly_from_divisor(S, 12, 0), 12),
        'e':   (poly_from_divisor(E, 11, 3), 14),
    }
    C = [redq(c) for c in h0_top(facs, (0, 9, 12, 14))]
    # depth-1 initial form (support 8):
    j0 = sp.factor(C[0])
    # saturated test: does the tie tower admit X,S,E,r with X S E != 0 ?
    killed4, _ = is_unit_ideal(C[:4] + [Q, w * X * S * E - 1], (X, S, E, r, w))
    killed3, _ = is_unit_ideal(C[:3] + [Q, w * X * S * E - 1], (X, S, E, r, w))
    return {
        'name': 'A  a11_b3000_T1 (0,9,12,14) all-defect-0 (both-tight)',
        'initial_form': j0,
        'depth_needed': 3 if killed3 else (4 if killed4 else None),
        'verdict': 'KILL' if killed4 else 'not killed by j0..3',
        'note': 'd1,e forced defect-0; sigma defect-0 is the both-tight choice.',
    }


# ---------------------------------------------------------------------------
# PILOT B -- a14_b0000_T2, state (6, d1==0, 12, 14).
#   sigma FORCED defect 0: v_t(sigma) = w = 3a-30 = 12 = deg sigma
#   e     FORCED defect 0: v_t(e) = a = 14 = deg e  (b0000, deg_E = 0)
#   d1 == 0 (T2);  d2 FREE degree 6.  support 0, tie depth 14.
# The two torus branches come from j0 = 12 S^2 (4 D6^2 + 9 S)^2 (5 D6^2 + 9 S).
# ---------------------------------------------------------------------------
def pilot_B():
    S, E, w = sp.symbols('S E w')
    Dc = sp.symbols('D0:7')                     # d2 coefficients, D6 = lc
    d2 = sum(Dc[i] * y**i for i in range(7))
    facs = {
        'd2':  (d2, 6),
        'sig': (sp.expand(S * (y + 1)**12), 12),
        'e':   (sp.expand(E * (y + 1)**14), 14),
    }
    C = h0_top(facs, (6, 0, 12, 14), drop_d1=True)
    j0 = sp.factor(C[0])
    branches = {'5 D6^2 + 9 S': sp.Rational(-5, 9), '4 D6^2 + 9 S': sp.Rational(-4, 9)}
    verdicts = {}
    for label, Sval in branches.items():
        sub = {Dc[6]: 1, S: Sval}               # D6 = 1 fixes the weight scaling
        eqs = [sp.expand(c.subs(sub)) for c in C[1:]]
        eqs = [e for e in eqs if e != 0]
        killed, _ = is_unit_ideal(eqs + [w * E - 1],
                                  list(Dc[:6]) + [E, w])
        verdicts[label] = 'KILL' if killed else 'survives'
    return {
        'name': 'B  a14_b0000_T2 (6,-,12,14) FORCED all-defect-0, d2 free deg 6',
        'initial_form': j0,
        'branch_verdicts': verdicts,
        'verdict': 'KILL' if all(v == 'KILL' for v in verdicts.values()) else 'partial',
        'note': 'sigma,e FORCED defect-0 (not a choice); whole-state kill.',
    }


# ---------------------------------------------------------------------------
# CONTROL -- reconstruct d1 defect-0 only; leave e GENERIC (defect free).
# Demonstrates the engine does NOT over-kill: the deep tie is solved by e's free
# top coefficients (triangular, one per level).  UNOBSTRUCTED.
# ---------------------------------------------------------------------------
def control_unobstructed():
    D, X, c9 = sp.symbols('D X c9')
    Cc = list(sp.symbols('c0:9')) + [c9]
    Ee = sp.symbols('e0:15')
    sig = (y + 1)**3 * sum(Cc[i] * y**i for i in range(10))   # v_t(sigma)=3 (forced floor)
    egen = sum(Ee[i] * y**i for i in range(15))               # e GENERIC deg 14
    facs = {'d2': (D, 0), 'd1': (poly_from_divisor(X, 5, 4), 9),
            'sig': (sp.expand(sig), 12), 'e': (egen, 14)}
    C = [redq(c) for c in h0_top(facs, (0, 9, 12, 14))]
    c9p = sp.Poly(4 * c9**3 + 1, c9)             # j0 with X=1 forces 4 c9^3 + 1 = 0

    def rc(e):
        e = sp.rem(sp.Poly(sp.expand(e), r), QP).as_expr()
        return sp.rem(sp.Poly(sp.expand(e), c9), c9p).as_expr()

    val = {X: sp.Integer(1), Cc[9]: c9}
    import random
    random.seed(7)
    for c in Cc[:9]:
        val[c] = sp.Rational(random.randint(-3, 3) or 1, random.randint(1, 3))
    val[D] = sp.Integer(2)
    val[Ee[0]] = sp.Integer(1)
    newest = {j: Ee[15 - j] for j in range(1, 14)}   # j1->e14 ... j13->e2
    carry, pivots = {}, []
    for j in range(1, 14):
        eq = rc(C[j].subs(val).subs(carry))
        v = newest[j]
        P = sp.Poly(eq, v)
        a = rc(P.coeff_monomial(v))
        pivots.append(a != 0)
        inv = sp.invert(sp.Poly(a, c9), c9p).as_expr() if a.has(c9) else 1 / a
        carry[v] = rc(-P.coeff_monomial(sp.Integer(1)) * inv)
    return {
        'name': 'CONTROL  d1 defect-0, e GENERIC (not reconstructed)',
        'verdict': 'UNOBSTRUCTED' if all(pivots) else 'inconclusive',
        'note': ('all 13 tower pivots nonzero: the depth-14 tie solves for '
                 'e14..e2 with lc(d1)=1, lc(sigma)=c9!=0.  Reconstructing e '
                 'removes exactly these free coefficients.'),
    }


def main():
    for fn in (pilot_A, pilot_B, control_unobstructed):
        res = fn()
        print('=' * 72)
        print('PILOT', res['name'])
        for k, v in res.items():
            if k == 'name':
                continue
            print(f'  {k}: {v}')
    print('=' * 72)


if __name__ == '__main__':
    main()
