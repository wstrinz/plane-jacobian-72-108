#!/usr/bin/env python3
"""d2-freedom threshold audit: do the four PHASE_F2_SCALE survivors die when the
level-0 tie tower is pushed past the census's cost-driven depth cap of 2, with
d2 left FULLY FREE (the sound / conservative choice)?

The four survivors are the deg d2 in {5,6} states of the two entirely-defect-0
T2 branches
    a11_b3100_T2   (two distinct marked roots r1 != r2, field Q[r1,r2]/(q,q))
    a12_b1110_T2   (single marked complement root r,   field Q[r]/(q))
Both branches force (audited):
  * e  defect-0 by the b-pattern:   v_t(e)=a, v_{r_i}(e)=b_i           (ALT_INF_SWEEP/b-pattern)
  * sigma defect-0 by the unique finite-place witness (nsplit==1):
      a11_b3100_T2: sigma = S (y+1)^3 (y-r1)^7 (y-r2)^2               (ALT_COMBINED fw)
      a12_b1110_T2: sigma = S (y+1)^6 prod_{3 roots}(y-r_i)^2         (ALT_COMBINED fw)
  * d1 == 0 (T2).
d2 is NOT constrained by any audited artifact (the finite-place cones of
ALT_REGIME_L2 sec 2 bound only v_P(d1) [T1] and v_P(sigma); PHASE_F_DEFECTS.md's
alt per-delta table lists ONLY d1 and sigma).  So d2's defect delta_d2 == deg d2
in {5,6}: there is NOTHING forced to impose on d2.  We therefore leave d2 fully
free (conservative: extra d2 freedom only enlarges the solution set, so a kill
under free d2 is sound) and push the depth.

Engine is independent of phase_f2_scale.py (only cascade_engine.MONOMIALS[0], an
audited artifact, is shared -- same as phase_f2_pilot.py).  READ-ONLY; uncommitted.
"""
import sys, time
sys.path.insert(0, '.')
import sympy as sp
from cascade_engine import MONOMIALS

y, r, r1, r2, S, E, w, wd = sp.symbols('y r r1 r2 S E w wd')
Q_COEFFS = [2048, -512, 320, -240, 195]

def qpoly(v):
    return sum(c * v**(4 - i) for i, c in enumerate(Q_COEFFS))

QY = qpoly(y)

# T2 monomials of h_0 (drop every d1-carrying monomial, x==0):
T2_MONS = [((k, z, b), int(c)) for (k, x, z, b), c in MONOMIALS[0] if x == 0]

TD = 60   # every T2 monomial has y-degree 6k+12z+15b == 60 at the deg d2=6 caps;
          # for deg d2 < 6 the max over monomials is still 60 (sig^5, e^4).

# ---- truncated reversed-series machinery (top ntop coeffs, hi->lo) ----------
def rev(poly, deg, ntop):
    P = sp.Poly(sp.expand(poly), y)
    return [P.coeff_monomial(y**(deg - i)) if deg - i >= 0 else sp.Integer(0)
            for i in range(ntop)]

def smul(a, b, ntop):
    out = [sp.Integer(0)] * ntop
    for i, ai in enumerate(a):
        if ai == 0: continue
        for j, bj in enumerate(b):
            if i + j >= ntop or bj == 0: continue
            out[i + j] += ai * bj
    return [sp.expand(c) for c in out]

def spow(a, n, ntop):
    res = [sp.Integer(1)] + [sp.Integer(0)] * (ntop - 1)
    base = list(a)
    while n:
        if n & 1: res = smul(res, base, ntop)
        n >>= 1
        if n: base = smul(base, base, ntop)
    return res

def h0_top(d2_poly, deg_d2, sig_poly, e_poly, ntop):
    """Top ntop coefficients of h_0 (T2) about y=infinity, index i = y^{60-i}."""
    revs = {'d2': rev(d2_poly, deg_d2, ntop) if deg_d2 is not None else None,
            'sig': rev(sig_poly, 12, ntop),
            'e': rev(e_poly, 15, ntop)}
    cache = {}
    def fp(nm, n):
        if n == 0: return [sp.Integer(1)] + [sp.Integer(0)] * (ntop - 1)
        if (nm, n) not in cache: cache[(nm, n)] = spow(revs[nm], n, ntop)
        return cache[(nm, n)]
    total = [sp.Integer(0)] * ntop
    for (k, z, b), coef in T2_MONS:
        if k and deg_d2 is None: continue
        dM = (deg_d2 or 0) * k + 12 * z + 15 * b
        shift = TD - dM
        if shift < 0 or shift >= ntop: continue
        term = [sp.Integer(1)] + [sp.Integer(0)] * (ntop - 1)
        for cnt, nm in ((k, 'd2'), (z, 'sig'), (b, 'e')):
            if cnt: term = smul(term, fp(nm, cnt), ntop)
        term = [coef * c for c in term]
        for i in range(ntop - shift):
            total[shift + i] += term[i]
    return [sp.expand(c) for c in total]

# ---- field reducers ---------------------------------------------------------
def reducer(root_vars):
    if not root_vars:
        return lambda e: sp.expand(e)
    if root_vars == [r]:
        QR = sp.Poly(qpoly(r), r)
        return lambda e: sp.rem(sp.Poly(sp.expand(e), r), QR).as_expr()
    QR1 = sp.Poly(qpoly(r1), r1); QR2 = sp.Poly(qpoly(r2), r2)
    def red2(e):
        e = sp.rem(sp.Poly(sp.expand(e), r1), QR1).as_expr()
        e = sp.rem(sp.Poly(sp.expand(e), r2), QR2).as_expr()
        return e
    return red2

# ---- reconstruction of the two survivor branches ---------------------------
def build_state(bid, deg_d2):
    """Return (d2_poly, sig_poly, e_poly, root_vars, Dc)."""
    if bid == 'a11_b3100_T2':
        # two distinct marked roots: r1 (b=3), r2 (b=1)
        e_poly = E * (y + 1)**11 * (y - r1)**3 * (y - r2)**1
        sig_poly = S * (y + 1)**3 * (y - r1)**7 * (y - r2)**2
        root_vars = [r1, r2]
    elif bid == 'a12_b1110_T2':
        # three equal-mult roots -> complement of the single unmarked root r:
        #   prod_{3 active}(y-r_i) = q/(2048 (y-r))
        comp = sp.div(sp.Poly(QY, y), sp.Poly(y - r, y))[0].as_expr() / 2048
        e_poly = E * (y + 1)**12 * comp**1
        sig_poly = S * (y + 1)**6 * comp**2
        root_vars = [r]
    else:
        raise ValueError(bid)
    if deg_d2 is None:
        d2_poly, Dc = None, None
    elif deg_d2 == 0:
        D = sp.Symbol('D'); d2_poly, Dc = D, [D]
    else:
        Dc = list(sp.symbols(f'D0:{deg_d2 + 1}'))
        d2_poly = sum(Dc[i] * y**i for i in range(deg_d2 + 1))
    return sp.expand(d2_poly) if d2_poly is not None else None, \
           sp.expand(sig_poly), sp.expand(e_poly), root_vars, Dc

# ---- depth-pushing saturated-GB kill test ----------------------------------
def kill_test(bid, deg_d2, max_depth, budget=600.0, verbose=True):
    d2_poly, sig_poly, e_poly, root_vars, Dc = build_state(bid, deg_d2)
    red = reducer(root_vars)
    C = [red(c) for c in h0_top(d2_poly, deg_d2, sig_poly, e_poly, max_depth)]

    # which d2 coefficients can appear in the top-`max_depth` window?
    d2_ring = []
    if Dc is not None:
        # d2^k monomial with k>=1 sits at y-degree 6k+... ; its sub-leading
        # coeffs shift it down.  Keep the top `max_depth` d2 coeffs to be safe.
        d2_ring = Dc[max(0, len(Dc) - max_depth):]
    sat_scalars = [S, E] + ([Dc[-1]] if Dc is not None else [])
    order_vars = [S, E] + list(d2_ring) + list(root_vars)
    root_gens = [qpoly(v) for v in root_vars]
    distinct = []
    if root_vars == [r1, r2]:
        distinct = [wd * (r1 - r2) - 1]
        order_vars = order_vars + [wd]
    order_vars = order_vars + [w]
    sat = w * sp.prod(sat_scalars) - 1

    results = []
    gens = []
    t0 = time.time()
    for d in range(max_depth):
        if C[d] != 0:
            gens.append(C[d])
        if not gens:
            results.append((d + 1, None, 0.0)); continue
        ts = time.time()
        G = sp.groebner(gens + root_gens + distinct + [sat],
                        *order_vars, order='grevlex')
        el = time.time() - ts
        unit = sp.Integer(1) in G.exprs
        results.append((d + 1, unit, round(el, 2)))
        if verbose:
            print(f'    depth {d+1:2d}: unit_ideal={unit}  ({el:.1f}s, '
                  f'ngens={len(gens)})', flush=True)
        if unit:
            return {'verdict': 'KILLED', 'kill_depth': d + 1,
                    'elapsed': round(time.time() - t0, 1), 'trace': results}
        if time.time() - t0 > budget:
            return {'verdict': 'PENDING_TIMEOUT', 'depth_reached': d + 1,
                    'elapsed': round(time.time() - t0, 1), 'trace': results}
    return {'verdict': 'NARROWED_or_UNOBSTRUCTED', 'depth_reached': max_depth,
            'elapsed': round(time.time() - t0, 1), 'trace': results}

def j0_form(bid, deg_d2):
    d2_poly, sig_poly, e_poly, root_vars, Dc = build_state(bid, deg_d2)
    C = h0_top(d2_poly, deg_d2, sig_poly, e_poly, 1)
    red = reducer(root_vars)
    return red(C[0])

if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--bid', required=True)
    ap.add_argument('--degd2', type=int, required=True)
    ap.add_argument('--maxdepth', type=int, default=17)
    ap.add_argument('--budget', type=float, default=600.0)
    args = ap.parse_args()
    print(f'{args.bid}  deg d2 = {args.degd2}  (tie depth 17)')
    print('  j0 =', sp.factor(j0_form(args.bid, args.degd2)))
    res = kill_test(args.bid, args.degd2, args.maxdepth, budget=args.budget)
    print('  VERDICT:', res['verdict'], res.get('kill_depth', ''),
          f"({res['elapsed']}s)")
