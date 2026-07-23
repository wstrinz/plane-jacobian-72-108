#!/usr/bin/env python3
"""Independent verifier for the Phase-F2 SUB2 divisor-reconstruction kills.

Does NOT import phase_f2_sub2 or convolution_descent.  It re-derives TWO NEW
marked-root sub2 KILLS from the a9_b1000_T1 cell by an independent route:

  * it reads the AUDITED source h_f formulas (t5_90t1_verify.load_h) and the
    master modulus Phi = c (y+1)^30 q, c = -1/6630;
  * it reconstructs the state's forced defect-0 factors directly
        e     = E (y+1)^9 (y - r)      (v_t = a = 9, v_marked = b_0 = 1)
        d1    = X (y - r)              (v_marked = 1)
        sigma = S                      (constant)
        d2    = 0  (K1)  or  d2 = D constant (K2);
  * it extracts the TWO top coefficients of the master identity
        f31(y) = sum_{f=0}^7 Phi^f e^(21-3f) h_f
    by an INDEPENDENT reversed-series short convolution (each factor is reversed
    about y=infinity and the per-f products are aligned by absolute degree),
    then reduces modulo q(r) with the single marked root r adjoined;
  * it certifies the saturated Groebner KILL over Q[r]/(q), and additionally
    walks the explicit HAND-STYLE elimination chain that closes the kill:
        top coeff  ==>  E^17 is forced to a fixed nonzero RATIONAL,
        next coeff ==>  r is then forced to a RATIONAL value,
        but q is irreducible over Q  ==>  r not in Q  ==>  CONTRADICTION.

Checks (K1 d2==0, K2 d2==D free constant):
  V0  reconstructed factors realise the claimed defect-0 divisors;
  Vq  q is irreducible over Q (so a marked root is never rational);
  V1  top coefficient depends only on E and forces E^17 rational (E != 0);
  V2  saturated ideal (c_top, c_next, q(r), sat) = unit ideal  ==>  KILL;
  V3  hand chain: after E^17 is fixed, c_next forces r rational (nonzero r
      coefficient), contradicting Vq.

Exits nonzero on any failure.
"""
import sys
import json
sys.path.insert(0, '.')
import sympy as sp
import t5_90t1_verify as base

y = base.y
r, X, S, E, D, w = sp.symbols('r X S E D w')
C = sp.Rational(-1, 6630)
PHI = sp.expand(C * (y + 1)**30 * base.q)
QR = sp.Poly(2048 * r**4 - 512 * r**3 + 320 * r**2 - 240 * r + 195, r)
H = base.load_h()
NTOP = 2

results = []
def check(tag, cond, detail=''):
    ok = bool(cond)
    results.append(ok)
    print(f'[{"PASS" if ok else "FAIL"}] {tag}' + (f'  -- {detail}' if detail else ''),
          flush=True)


# ------- independent reversed-series top coefficients of the master identity --
def rev(poly):
    P = sp.Poly(sp.expand(poly), y)
    d = P.degree()
    return d, [P.coeff_monomial(y**(d - i)) for i in range(NTOP)]


def cmul(a, b):
    out = [sp.Integer(0)] * NTOP
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            if i + j < NTOP and bj != 0:
                out[i + j] += ai * bj
    return [sp.expand(c) for c in out]


def cpow(a, n):
    res = [sp.Integer(1)] + [sp.Integer(0)] * (NTOP - 1)
    for _ in range(n):
        res = cmul(res, a)
    return res


def master_top(polys):
    """Top NTOP coefficients of f31(y); returns (TD, [c_TD, c_{TD-1}, ...])."""
    d0 = sp.expand((polys['d2']**2 + polys['sigma']) / 4)
    subs = {base.d0: d0, base.d1: polys['d1'], base.d2: polys['d2'],
            base.e: polys['e']}
    dPhi, revPhi = rev(PHI)
    dE, revE = rev(polys['e'])
    terms = []
    for f in range(8):
        hf = sp.expand(H[f].subs(subs))
        if hf == 0:
            continue
        dh, revH = rev(hf)
        prod = cpow(revPhi, f)
        prod = cmul(prod, cpow(revE, 21 - 3 * f))
        prod = cmul(prod, revH)
        TDf = dPhi * f + dE * (21 - 3 * f) + dh
        terms.append((TDf, prod))
    TD = max(t[0] for t in terms)
    out = [sp.Integer(0)] * NTOP
    for TDf, prod in terms:
        shift = TD - TDf
        for i in range(NTOP - shift):
            if shift + i < NTOP:
                out[shift + i] += prod[i]
    return TD, [sp.expand(o) for o in out]


def redq(expr):
    return sp.rem(sp.Poly(sp.expand(expr), r), QR).as_expr()


# ------- Vq : q irreducible over Q -------------------------------------------
qfac = sp.factor_list(QR.as_expr())
check('Vq q(r) irreducible over Q (marked root never rational)',
      len(qfac[1]) == 1 and qfac[1][0][1] == 1,
      f'factors={[str(f) for f, _ in qfac[1]]}')


def run_kill(tag, d2_expr, sat_scalars, order):
    polys = {
        'e': sp.expand(E * (y + 1)**9 * (y - r)),
        'd1': sp.expand(X * (y - r)),
        'sigma': S,
        'd2': d2_expr,
    }
    # V0 defect-0 divisor degrees
    check(f'{tag}.V0 e defect0 (deg 10 = 9 + 1)',
          sp.Poly(polys['e'], y).degree() == 10)
    check(f'{tag}.V0 d1 defect0 (deg 1 = 0 + 1)',
          sp.Poly(polys['d1'], y).degree() == 1)

    TD, Craw = master_top(polys)
    c_top = redq(Craw[0])
    c_next = redq(Craw[1])
    check(f'{tag}.V1 top degree == 250', TD == 250)

    # V1 : top coefficient is pure in E and forces E^17 to a fixed rational
    ctop_E = sp.factor(c_top)
    poly_in = sp.Poly(sp.expand(c_top / E**8), E)
    # expect A*E^17 + B, A,B nonzero rationals
    A = poly_in.coeff_monomial(E**17)
    B = poly_in.coeff_monomial(1)
    is_pure = (sp.expand(c_top - E**8 * (A * E**17 + B)) == 0
               and A != 0 and B != 0 and not c_top.has(r))
    check(f'{tag}.V1 top coeff = E^8 (A E^17 + B), A,B rational nonzero, r-free',
          is_pure, f'A={A}, B={B}')

    # V2 : saturated Groebner kill over Q[r]/(q)
    sat = w * sp.prod(sat_scalars) - 1
    G = sp.groebner([c_top, c_next, QR.as_expr(), sat], *order, order='grevlex')
    check(f'{tag}.V2 KILL over Q[r]/(q): (c_top,c_next,q(r),sat) = unit ideal',
          sp.Integer(1) in G.exprs)

    # V3 : explicit hand chain.  Substitute the forced E^17 = -B/A into c_next.
    E17 = sp.Rational(-B, A)
    # c_next = E^8 [ (a r + b) E^17 + (c r + d) ]; reduce using E^17 = E17.
    cn = sp.expand(c_next / E**8)
    cn = sp.Poly(cn, E)
    # collect the E^17 part and the E^0 part (cn has only E^17 and E^0 terms)
    part17 = cn.coeff_monomial(E**17)
    part0 = cn.coeff_monomial(1)
    check(f'{tag}.V3a c_next / E^8 involves only E^17 and E^0 powers',
          sp.expand(cn.as_expr() - (part17 * E**17 + part0)) == 0)
    lin = sp.expand(part17 * E17 + part0)          # linear in r after E^17 fixed
    linP = sp.Poly(lin, r)
    coeff_r = linP.coeff_monomial(r)
    const_r = linP.coeff_monomial(1)
    check(f'{tag}.V3b after E^17 fixed, c_next is linear in r with NONZERO '
          f'r-coefficient (forces r rational)', coeff_r != 0 and linP.degree() <= 1,
          f'r-coeff={coeff_r}')
    r_forced = sp.Rational(-const_r, coeff_r) if coeff_r != 0 else None
    check(f'{tag}.V3c forced r = {r_forced} is rational, but q(r)=0 has no '
          f'rational root  ==>  CONTRADICTION (KILL)',
          r_forced is not None and QR.as_expr().subs(r, r_forced) != 0)


print('=' * 64)
print('K1 : a9_b1000_T1, d2 identically zero, state (deg d1,sigma,e)=(1,0,10)')
run_kill('K1', sp.Integer(0), [E, X, S], [X, S, E, r, w])
print('-' * 64)
print('K2 : a9_b1000_T1, d2 = D free nonzero constant (imposed reconstruction '
      'still kills)')
run_kill('K2', D, [E, X, S, D], [X, S, E, D, r, w])

print('=' * 64)
if all(results):
    print('ALL CHECKS PASS')
    sys.exit(0)
print('SOME CHECKS FAILED')
sys.exit(1)
