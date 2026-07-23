#!/usr/bin/env python3
"""Independent verifier for the Phase-F2 pilot (phase_f2_pilot.py).

It does NOT import phase_f2_pilot.  It re-derives the Pilot A chain by an
independent route: it reads the audited level-0 graded coefficient
h_0 = cascade_engine.MONOMIALS[0], reconstructs the state's defect-0 factors,
and extracts the top coefficients of h_0(y) by reversing each SMALL factor
polynomial (deg <= 14) and short-convolving -- a different implementation from
the pilot's engine and cheap enough to run under load.  It checks:

  (V0) the reconstructed factors realise the claimed defect-0 divisors;
  (V1) the level-0 initial form equals the AUDITED support-8 form
       2187 S^2 (4 S^3 + X^4) recorded in alt_residue_congruences.json;
  (V2) the exact KILL: (j0,j1,j2, q(r), w*X*S*E - 1) is the unit ideal;
  (V3) the kill is universal in the marked point (q(r)=0 not imposed);
  (V4) soundness control: with e GENERIC the depth-3 slice stays solvable
       (nonzero triangular pivots) -- the kill is a consequence of divisor
       reconstruction, not of the engine.

Exits nonzero on any failure.
"""
import sys
import json
sys.path.insert(0, '.')
import sympy as sp
from cascade_engine import MONOMIALS

y, r, D, X, S, E, w = sp.symbols('y r D X S E w')
Q = 2048 * r**4 - 512 * r**3 + 320 * r**2 - 240 * r + 195
QP = sp.Poly(Q, r)
NTOP = 3                       # depth-3 truncation suffices for the kill


def rev(poly, deg):
    """Top-NTOP coefficients of `poly` (y-degree `deg`), index i -> y^{deg-i}."""
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


def top_h0(fac_deg, drop_d1=False):
    """Top-NTOP coefficients of h_0(y); fac_deg: name -> (poly, deg)."""
    R = {nm: rev(f, dg) for nm, (f, dg) in fac_deg.items()}
    degs = {nm: dg for nm, (f, dg) in fac_deg.items()}
    total_deg = 60
    out = [sp.Integer(0)] * NTOP
    for (k, x, z, b), coef in MONOMIALS[0]:
        if drop_d1 and x != 0:
            continue
        dM = degs['d2'] * k + degs['d1'] * x + degs['sig'] * z + degs['e'] * b
        shift = total_deg - dM
        if shift < 0 or shift >= NTOP:
            continue
        term = [sp.Integer(1)] + [sp.Integer(0)] * (NTOP - 1)
        for cnt, nm in ((k, 'd2'), (x, 'd1'), (z, 'sig'), (b, 'e')):
            if cnt:
                term = cmul(term, cpow(R[nm], cnt))
        for i in range(NTOP - shift):
            out[shift + i] += sp.Integer(int(coef)) * term[i]
    return [sp.expand(c) for c in out]


def redq(e):
    return sp.rem(sp.Poly(sp.expand(e), r), QP).as_expr()


def finv(a):
    return sp.invert(sp.Poly(a, r), QP).as_expr()


results = []


def check(tag, cond, detail=''):
    results.append(bool(cond))
    print(f'[{"PASS" if cond else "FAIL"}] {tag}' + (f'  -- {detail}' if detail else ''),
          flush=True)


# --- Pilot A: a11_b3000_T1 (0,9,12,14), all-defect-0 -------------------------
d1 = sp.expand(X * (y + 1)**5 * (y - r)**4)
sig = sp.expand(S * (y + 1)**12)
e = sp.expand(E * (y + 1)**11 * (y - r)**3)
fac = {'d2': (D, 0), 'd1': (d1, 9), 'sig': (sig, 12), 'e': (e, 14)}

check('V0 d1 defect0 (deg 9 = 5+4)', sp.Poly(d1, y).degree() == 9 == 5 + 4)
check('V0 sigma defect0 (deg 12 = 12)', sp.Poly(sig, y).degree() == 12)
check('V0 e defect0 (deg 14 = 11+3)', sp.Poly(e, y).degree() == 14 == 11 + 3)

C = [redq(c) for c in top_h0(fac)]

j0 = sp.factor(C[0])
expected = sp.factor(2187 * S**2 * (4 * S**3 + X**4))
check('V1 initial form == 2187 S^2 (4 S^3 + X^4)', sp.simplify(j0 - expected) == 0,
      f'got {j0}')
try:
    cat = json.load(open('alt_residue_congruences.json'))['support_catalog']
    s8 = next(x for x in cat if x['support_id'] == 8)
    fac8 = sp.factor(sp.sympify(s8['factored'], locals={'D': D, 'X': X, 'S': S, 'E': E}))
    check('V1 matches alt_residue_congruences.json support 8',
          sp.simplify(fac8 - expected) == 0, f'json: {s8["factored"]}')
except Exception as ex:                                     # pragma: no cover
    check('V1 json cross-check', False, str(ex))

G = sp.groebner(C + [Q, w * X * S * E - 1], X, S, E, r, w, order='grevlex')
check('V2 KILL: (j0,j1,j2,q, sat XSE) = unit ideal', sp.Integer(1) in G.exprs,
      f'GB={list(G.exprs)[:3]}')

G2 = sp.groebner(C + [w * X * S * E - 1], X, S, E, r, w, order='grevlex')
check('V3 KILL universal in marked point r (no q imposed)',
      sp.Integer(1) in G2.exprs)

# (V4) soundness control: e generic (top 3 coeffs symbolic) => depth-3 slice
# solvable at the torus point (X,S)=(4,-4) via triangular pivots.
Ee = sp.symbols('e12 e13 e14')
egen = Ee[2] * y**14 + Ee[1] * y**13 + Ee[0] * y**12
facg = {'d2': (D, 0), 'd1': (d1, 9), 'sig': (sig, 12), 'e': (egen, 14)}
Cg = [redq(c) for c in top_h0(facg)]
base = {X: sp.Integer(4), S: sp.Integer(-4)}
ok_init = (redq(Cg[0].subs(base)) == 0)
j1 = sp.Poly(redq(Cg[1].subs(base)), Ee[2])
piv1 = redq(j1.coeff_monomial(Ee[2]))
e14 = redq(-j1.coeff_monomial(sp.Integer(1)) * finv(piv1)) if piv1 != 0 else None
piv2 = 0
if e14 is not None:
    j2 = sp.Poly(redq(Cg[2].subs(base).subs({Ee[2]: e14})), Ee[1])
    piv2 = redq(j2.coeff_monomial(Ee[1]))
check('V4 control: e generic => depth-3 slice solvable, X,S,lc(e)!=0 '
      '(nonzero triangular pivots)', ok_init and piv1 != 0 and piv2 != 0)

print('-' * 60)
if all(results):
    print('ALL CHECKS PASS')
    sys.exit(0)
print('SOME CHECKS FAILED')
sys.exit(1)
