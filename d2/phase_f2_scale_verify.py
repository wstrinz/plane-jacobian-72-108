#!/usr/bin/env python3
"""Independent verifier for the Phase-F2 SCALE alt front (phase_f2_scale.py).

Does NOT import phase_f2_scale.  It re-derives TWO NEW divisor-reconstruction
KILLS from the flagship Galois-stable a11_b1111_T1 family by an independent
route: it reads the audited level-0 graded coefficient h_0 = MONOMIALS[0],
reconstructs the state's forced defect-0 factors directly (d1 propto (y+1)^5 q,
e propto (y+1)^11 q, sigma = S (y+1)^3 -- all over Q, no marked root), extracts
the top coefficients of h_0(y) by reversing each factor and short-convolving,
and certifies the saturated Groebner kill.  Everything is over Q (the b1111
divisor is Galois-stable: prod_i (y-r_i) = q/2048), so no q(r) reduction is
needed -- a genuinely different, simpler chain than phase_f2_pilot_verify.

Checks:
  K1  a11_b1111_T1, d2 identically zero, state (deg d1,sigma,e)=(9,3,15):
      * V0 factors realise the claimed defect-0 divisors (deg = valuation sum);
      * V1 level-0 initial form == -729 E (9 E^3 + 8 X^5), matching
        alt_residue_congruences.json support 15;
      * V2 KILL: (j0, j1, sat X S E) = unit ideal over Q  (depth 2).
  K2  same family, d2 = D a free (degree-0) constant:  even granting a nonzero
      free d2 leading scalar the tie tower dies -- (j0, j1, sat D X S E) = unit
      ideal over Q (depth 2).  A distinct kill: d2-freedom does not rescue it.

Exits nonzero on any failure.
"""
import sys, json
sys.path.insert(0, '.')
import sympy as sp
from cascade_engine import MONOMIALS

y, X, S, E, D, w = sp.symbols('y X S E D w')
Q = 2048 * y**4 - 512 * y**3 + 320 * y**2 - 240 * y + 195   # 2048 prod_i (y-r_i)
NTOP = 2


def rev(poly, deg):
    P = sp.Poly(sp.expand(poly), y)
    return [P.coeff_monomial(y**(deg - i)) if deg - i >= 0 else sp.Integer(0)
            for i in range(NTOP)]


def cmul(a, b):
    out = [sp.Integer(0)] * NTOP
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            if i + j < NTOP:
                out[i + j] += ai * bj
    return [sp.expand(c) for c in out]


def cpow(a, n):
    res = [sp.Integer(1)] + [sp.Integer(0)] * (NTOP - 1)
    for _ in range(n):
        res = cmul(res, a)
    return res


def top_h0(fac_deg, TD, drop_d2=False):
    R = {nm: rev(f, dg) for nm, (f, dg) in fac_deg.items()}
    degs = {nm: dg for nm, (f, dg) in fac_deg.items()}
    out = [sp.Integer(0)] * NTOP
    for (k, x, z, b), coef in MONOMIALS[0]:
        if drop_d2 and k != 0:
            continue
        dM = degs['d2'] * k + degs['d1'] * x + degs['sig'] * z + degs['e'] * b
        shift = TD - dM
        if shift < 0 or shift >= NTOP:
            continue
        term = [sp.Integer(1)] + [sp.Integer(0)] * (NTOP - 1)
        for cnt, nm in ((k, 'd2'), (x, 'd1'), (z, 'sig'), (b, 'e')):
            if cnt:
                term = cmul(term, cpow(R[nm], cnt))
        for i in range(NTOP - shift):
            out[shift + i] += sp.Integer(int(coef)) * term[i]
    return [sp.expand(c) for c in out]


results = []
def check(tag, cond, detail=''):
    results.append(bool(cond))
    print(f'[{"PASS" if cond else "FAIL"}] {tag}' + (f'  -- {detail}' if detail else ''),
          flush=True)


# Galois-stable q-part:  prod_i (y - r_i) = q / 2048.
qfac = sp.expand(Q / 2048)
d1 = sp.expand(X * (y + 1)**5 * qfac)       # v_t=5, each root mult 1  (deg 9)
sig = sp.expand(S * (y + 1)**3)             # v_t=3                    (deg 3)
e = sp.expand(E * (y + 1)**11 * qfac)       # v_t=11, each root mult 1 (deg 15)

check('V0 d1 defect0 (deg 9 = 5 + 4)', sp.Poly(d1, y).degree() == 9)
check('V0 sigma defect0 (deg 3)', sp.Poly(sig, y).degree() == 3)
check('V0 e defect0 (deg 15 = 11 + 4)', sp.Poly(e, y).degree() == 15)

# leading degree TD = max dM over monomials (d2 dropped for K1)
def TDof(degs, drop_d2):
    m = None
    for (k, x, z, b), coef in MONOMIALS[0]:
        if drop_d2 and k != 0:
            continue
        dM = degs[0]*k + degs[1]*x + degs[2]*z + degs[3]*b
        m = dM if m is None else max(m, dM)
    return m

# ---- K1 : d2 identically zero -------------------------------------------
TD1 = TDof((0, 9, 3, 15), True)
fac1 = {'d2': (sp.Integer(0), 0), 'd1': (d1, 9), 'sig': (sig, 3), 'e': (e, 15)}
C1 = top_h0(fac1, TD1, drop_d2=True)
j0 = sp.factor(C1[0])
expected = sp.factor(-729 * E * (9 * E**3 + 8 * X**5))
check('K1.V1 initial form == -729 E (9 E^3 + 8 X^5)',
      sp.simplify(j0 - expected) == 0, f'got {j0}')
try:
    cat = json.load(open('alt_residue_congruences.json'))['support_catalog']
    s15 = next(x for x in cat if x['support_id'] == 15)
    fac15 = sp.factor(sp.sympify(s15['factored'], locals={'D': D, 'X': X, 'S': S, 'E': E}))
    check('K1.V1 matches alt_residue_congruences.json support 15',
          sp.simplify(fac15 - expected) == 0, f'json: {s15["factored"]}')
except Exception as ex:
    check('K1.V1 json cross-check', False, str(ex))
G1 = sp.groebner([c for c in C1 if c != 0] + [w * X * S * E - 1],
                 X, S, E, w, order='grevlex')
check('K1.V2 KILL over Q (j0, j1, sat X S E) = unit ideal', sp.Integer(1) in G1.exprs)

# ---- K2 : d2 = D a free nonzero constant --------------------------------
TD2 = TDof((0, 9, 3, 15), False)
fac2 = {'d2': (D, 0), 'd1': (d1, 9), 'sig': (sig, 3), 'e': (e, 15)}
C2 = top_h0(fac2, TD2, drop_d2=False)
check('K2.V1 same initial form (d2 does not enter leading tie)',
      sp.simplify(sp.factor(C2[0]) - expected) == 0)
G2 = sp.groebner([c for c in C2 if c != 0] + [w * D * X * S * E - 1],
                 D, X, S, E, w, order='grevlex')
check('K2.V2 KILL over Q with free d2 (j0, j1, sat D X S E) = unit ideal',
      sp.Integer(1) in G2.exprs)

print('-' * 60)
if all(results):
    print('ALL CHECKS PASS')
    sys.exit(0)
print('SOME CHECKS FAILED')
sys.exit(1)
