#!/usr/bin/env python3
"""Phase F, work item F2 -- SCALING the divisor-reconstruction kill test.

Extends phase_f2_pilot.py from the 3 hand-worked pilots to the full frontier:

  (1) ALT FRONT: every alternate-regime survivor whose relevant divisors
      (d1, sigma, e) are ALL forced defect-0 (deg p == sum_s v_s(p)).  Because
      defect-0 forces the finite-place split UNIQUELY (verified: every such
      state has exactly one admissible ALT_REGIME_L2 sec.2 cone split), each of
      d1, sigma, e is  p = lambda_p * prod_s (y-s)^{v_s(p)}  EXACTLY.  We
      reconstruct them, feed the determined coefficients into the state's
      level-0 tie tower (support_id + L0_tie_depth from
      alt_residue_congruences.json), leave d2 a FREE polynomial of its state
      degree (conservative: enlarging d2 only shrinks the kill set, so a kill
      under free d2 is sound), saturate the leading scalars nonzero, and read
      the verdict off an exact saturated Groebner basis:
          unit ideal after saturation  ==>  KILL.
      Field of definition is minimised per state: pure t-place / Galois-stable
      full-q divisors compute over Q; a 3-equal-multiplicity q-factor is the
      complement q/(2048(y-r)) of a single marked root; a lone active root is a
      single marked root in Q[r]/(q); genuinely-two-distinct-root divisors use
      two marked roots r1,r2 with q(r1)=q(r2)=0 and r1!=r2 saturated.

  Flagship sub-run: the Galois-stable a11_b1111 family (d1 propto (y+1)^5 q,
  e propto (y+1)^11 q, sigma = S(y+1)^3) -- everything over Q, no marked root.

Verdicts are labelled PENDING AUDIT.  READ-ONLY on every audited artifact;
nothing here is committed.
"""
import sys
import json
import time
sys.path.insert(0, '.')
import sympy as sp
from cascade_engine import MONOMIALS

ROOT = '.'
y, r, r1, r2 = sp.symbols('y r r1 r2')
Q_COEFFS = [2048, -512, 320, -240, 195]           # 2048 y^4 -512 y^3 +320 y^2 -240 y +195


def qpoly(var):
    return sum(c * var**(4 - i) for i, c in enumerate(Q_COEFFS))


QY = qpoly(y)
QR = sp.Poly(qpoly(r), r)
QR1 = sp.Poly(qpoly(r1), r1)
QR2 = sp.Poly(qpoly(r2), r2)


# --------------------------------------------------------------------------
#  h_0 top-`ntop` coefficients about y = infinity by truncated convolution.
#  TD (leading degree) and the set of dropped monomials are per-state.
# --------------------------------------------------------------------------
def total_deg(degs, drop_d1, drop_sig):
    drop_d2 = degs[0] is None
    m = None
    for (k, x, z, b), coef in MONOMIALS[0]:
        if drop_d1 and x != 0:
            continue
        if drop_sig and z != 0:
            continue
        if drop_d2 and k != 0:
            continue
        dM = ((degs[0] or 0) * k + (degs[1] or 0) * x
              + (degs[2] or 0) * z + (degs[3] or 0) * b)
        m = dM if m is None else max(m, dM)
    return m


def rev_of_poly(poly, deg, ntop):
    P = sp.Poly(sp.expand(poly), y)
    return [P.coeff_monomial(y**(deg - i)) if deg - i >= 0 else sp.Integer(0)
            for i in range(ntop)]


def ser_mul(a, b, ntop):
    out = [sp.Integer(0)] * ntop
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            if i + j >= ntop or bj == 0:
                continue
            out[i + j] += ai * bj
    return [sp.expand(c) for c in out]


def ser_pow(a, n, ntop):
    res = [sp.Integer(1)] + [sp.Integer(0)] * (ntop - 1)
    base = list(a)
    while n:
        if n & 1:
            res = ser_mul(res, base, ntop)
        n >>= 1
        if n:
            base = ser_mul(base, base, ntop)
    return res


def h0_top(factors, degs, TD, ntop, drop_d1=False, drop_sig=False):
    revs = {nm: rev_of_poly(f, dg, ntop) for nm, (f, dg) in factors.items()}
    cache = {}

    def fp(nm, n):
        if n == 0:
            return [sp.Integer(1)] + [sp.Integer(0)] * (ntop - 1)
        if (nm, n) not in cache:
            cache[(nm, n)] = ser_pow(revs[nm], n, ntop)
        return cache[(nm, n)]

    drop_d2 = 'd2' not in factors
    total = [sp.Integer(0)] * ntop
    for (k, x, z, b), coef in MONOMIALS[0]:
        if drop_d1 and x != 0:
            continue
        if drop_sig and z != 0:
            continue
        if drop_d2 and k != 0:
            continue
        dM = degs[0] * k + degs[1] * x + degs[2] * z + degs[3] * b
        shift = TD - dM
        if shift < 0 or shift >= ntop:
            continue
        term = [sp.Integer(1)] + [sp.Integer(0)] * (ntop - 1)
        for cnt, nm in ((k, 'd2'), (x, 'd1'), (z, 'sig'), (b, 'e')):
            if cnt:
                term = ser_mul(term, fp(nm, cnt), ntop)
        term = [sp.Integer(int(coef)) * c for c in term]
        for i in range(ntop - shift):
            total[shift + i] += term[i]
    return [sp.expand(c) for c in total]


# --------------------------------------------------------------------------
#  root-ideal reduction (single r / two r1,r2)
# --------------------------------------------------------------------------
def make_reducer(root_vars):
    if not root_vars:
        return lambda e: sp.expand(e)
    if root_vars == [r]:
        return lambda e: sp.rem(sp.Poly(sp.expand(e), r), QR).as_expr()

    def red2(e):
        e = sp.rem(sp.Poly(sp.expand(e), r1), QR1).as_expr()
        e = sp.rem(sp.Poly(sp.expand(e), r2), QR2).as_expr()
        return e
    return red2


# --------------------------------------------------------------------------
#  divisor reconstruction.  q-part per polynomial from its per-root mults.
#  field is minimised: returns (poly, root_vars).
# --------------------------------------------------------------------------
def qpart(mults, rootvar_for_complement):
    """Build prod_i (y-r_i)^{mults[i]} using the minimal field.

    mults: dict {root_index -> multiplicity} (only positive entries).
    returns (poly_in_y, list_of_root_vars_used).
    Strategy per group of active indices sharing the polynomial:
      - all 4 active & equal mult m       -> (q/2048)^m over Q
      - 3 active & equal mult m           -> (q/(2048 (y-rc)))^m , rc single var
      - 1 active                          -> (y-r)^m , single var r
      - 2 active                          -> (y-r1)^m1 (y-r2)^m2 , vars r1,r2
    """
    active = sorted(mults)
    if not active:
        return sp.Integer(1), []
    if len(active) == 4 and len(set(mults.values())) == 1:
        m = next(iter(mults.values()))
        return sp.expand((QY / 2048)**m), []
    if len(active) == 3 and len(set(mults.values())) == 1:
        m = next(iter(mults.values()))
        rc = rootvar_for_complement
        base = sp.div(sp.Poly(QY, y), sp.Poly((y - rc), y))[0].as_expr() / 2048
        return sp.expand(base**m), [rc]
    if len(active) == 1:
        m = mults[active[0]]
        return sp.expand((y - r)**m), [r]
    if len(active) == 2:
        vs = [r1, r2]
        poly = sp.Integer(1)
        used = []
        for k, idx in enumerate(active):
            poly *= (y - vs[k])**mults[idx]
            used.append(vs[k])
        return sp.expand(poly), used
    raise ValueError(f'unsupported active-root pattern {mults}')


def reconstruct(a, b, split, branch, degs, drop_d1, drop_sig):
    """Return factors dict {d2,d1,sig,e -> (poly,deg)} plus root_vars, scalars.

    split for T1: tuple of (x,z) per place [t, r0..r3]; for T2: tuple of z per
    place (d1==0).  b: branch b-pattern (v_root(e)).  a: v_t(e).
    """
    D, X, S, E = sp.symbols('D X S E')
    # per-root multiplicities for e, d1, sigma
    e_m = {i: b[i] for i in range(4) if b[i] > 0}
    if branch == 'T1':
        x_t = split[0][0]
        z_t = split[0][1] or 0
        d1_m = {i: split[i + 1][0] for i in range(4) if split[i + 1][0] > 0}
        sig_m = {i: (split[i + 1][1] or 0) for i in range(4)
                 if (split[i + 1][1] or 0) > 0}
    else:  # T2 : d1 == 0 ; split is z per place
        z_t = split[0]
        x_t = 0
        d1_m = {}
        sig_m = {i: split[i + 1] for i in range(4) if split[i + 1] > 0}

    # decide a single shared complement variable if any poly needs the
    # 3-equal-mult complement pattern (they share the same excluded root).
    def complement_root(mults):
        active = set(mults)
        if len(active) == 3 and len(set(mults.values())) == 1:
            return (set(range(4)) - active).pop()
        return None
    comp_idx = None
    for mm in (e_m, d1_m, sig_m):
        ci = complement_root(mm)
        if ci is not None:
            comp_idx = ci
    # the complement variable name r is reused across polys (same excluded root)
    root_vars = set()

    def build(scalar, t_mult, mults):
        qp, rv = qpart(mults, r)      # complement uses r
        root_vars.update(rv)
        return sp.expand(scalar * (y + 1)**t_mult * qp)

    e_poly = build(E, a, e_m)
    factors = {'e': (e_poly, degs[3])}
    if not drop_d1:
        factors['d1'] = (build(X, x_t, d1_m), degs[1])
    if not drop_sig:
        factors['sig'] = (build(S, z_t, sig_m), degs[2])
    # d2 : free polynomial of degree deg_d2 (conservative)
    Dc = None
    if degs[0] is not None:
        if degs[0] == 0:
            factors['d2'] = (D, 0)
            Dc = [D]
        else:
            Dc = list(sp.symbols(f'D0:{degs[0] + 1}'))
            factors['d2'] = (sum(Dc[i] * y**i for i in range(degs[0] + 1)), degs[0])
    scalars = [s for s, present in
               ((X, not drop_d1), (S, not drop_sig), (E, True)) if present]
    # d2 leading scalar for saturation only if d2 actually enters the window
    return factors, sorted(root_vars, key=str), scalars, Dc


# --------------------------------------------------------------------------
#  per-state kill test (incremental depth, saturated GB)
# --------------------------------------------------------------------------
def d2_in_window(degs, TD, depth, drop_d1, drop_sig):
    if degs[0] is None:
        return False
    for (k, x, z, b), coef in MONOMIALS[0]:
        if k == 0:
            continue
        if drop_d1 and x != 0:
            continue
        if drop_sig and z != 0:
            continue
        dM = (degs[0] * k + (degs[1] or 0) * x
              + (degs[2] or 0) * z + (degs[3] or 0) * b)
        if TD - depth + 1 <= dM <= TD:
            return True
    return False


def kill_test(rec, TD, depth, degs, drop_d1, drop_sig, budget=90.0,
              max_depth=None):
    factors, root_vars, scalars, Dc = rec
    ntop = depth if max_depth is None else min(depth, max_depth)
    degs_num = tuple(d if d is not None else 0 for d in degs)
    C = h0_top(factors, degs_num,
               TD, ntop, drop_d1=drop_d1, drop_sig=drop_sig)
    red = make_reducer(root_vars)
    C = [red(c) for c in C]

    # saturation variables
    w = sp.Symbol('w')
    sat_scalars = list(scalars)
    d2_needed = d2_in_window(degs, TD, ntop, drop_d1, drop_sig)
    # only the top-ntop d2 coefficients can enter the truncated window
    d2_ring = []
    if Dc is not None and d2_needed:
        d2_ring = Dc[max(0, len(Dc) - ntop):]
        sat_scalars = sat_scalars + [Dc[-1]]     # lc(d2)
    elif Dc is not None and degs[0] == 0:
        d2_ring = list(Dc)
    order_vars = list(scalars) + list(d2_ring) + list(root_vars) + [w]
    sat = w * sp.prod(sat_scalars) - 1
    root_gens = []
    if r in root_vars:
        root_gens.append(qpoly(r))
    if r1 in root_vars:
        root_gens.append(qpoly(r1))
    if r2 in root_vars:
        root_gens.append(qpoly(r2))
    distinct = []
    if r1 in root_vars and r2 in root_vars:
        wd = sp.Symbol('wd')
        distinct = [wd * (r1 - r2) - 1]
        order_vars = order_vars[:-1] + [wd, w]

    j0_factored = sp.factor(C[0]) if not root_vars else C[0]

    t0 = time.time()
    kill_depth = None
    gens = []
    for kdx in range(len(C)):
        if C[kdx] != 0:
            gens.append(C[kdx])
        if not gens:
            continue
        allgens = gens + root_gens + distinct + [sat]
        try:
            G = sp.groebner(allgens, *order_vars, order='grevlex')
        except Exception as ex:
            return {'verdict': 'ERROR', 'error': str(ex),
                    'depth_reached': kdx + 1}
        if sp.Integer(1) in G.exprs:
            kill_depth = kdx + 1
            break
        if time.time() - t0 > budget:
            return {'verdict': 'PENDING_TIMEOUT', 'depth_reached': kdx + 1,
                    'elapsed': round(time.time() - t0, 1),
                    'j0': str(j0_factored)}
    if kill_depth is not None:
        return {'verdict': 'KILLED', 'kill_depth': kill_depth,
                'tie_depth': depth, 'elapsed': round(time.time() - t0, 1),
                'j0': str(j0_factored), 'field': field_label(root_vars)}
    # not killed within tracked coefficients
    return {'verdict': 'NARROWED_or_UNOBSTRUCTED', 'depth_reached': len(C),
            'tie_depth': depth, 'elapsed': round(time.time() - t0, 1),
            'j0': str(j0_factored), 'field': field_label(root_vars)}


def field_label(root_vars):
    if not root_vars:
        return 'Q'
    if root_vars == [r]:
        return 'Q[r]/(q)'
    return 'Q[r1,r2]/(q,q), r1!=r2'


# ==========================================================================
#  DRIVER : alt defect-0 census.  Cone splits (ALT_REGIME_L2 sec.2).
# ==========================================================================
def _t_pairs_T1(a, sz):
    if a == 11:
        return [(x, None) for x in range(5, 10)] if sz else [(x, z) for x in range(5, 10) for z in range(3, 13)]
    if a == 12:
        if sz: return [(9, None)]
        return [(3 + i, i) for i in range(6)] + [(9, z) for z in range(6, 13)]
    if a == 13: return []
    if a == 14:
        return [] if sz else [(6, 0), (7, 1), (8, 2), (9, 3)]
    return []


def _q_pairs_T1(b, sz):
    if b == 0: return [(0, 0)]
    if b == 1:
        return [(x, None) for x in range(3, 10)] if sz else [(1, 0), (2, 1)] + [(x, z) for x in range(3, 10) for z in range(2, 13)]
    if b == 2:
        return [(7, None)] if sz else [(7, z) for z in range(5, 13)]
    if b == 3: return [(4 + i, i) for i in range(6)]
    return []


def _enum_T1(a, b, X, Z, sz):
    places = [_t_pairs_T1(a, sz)] + [_q_pairs_T1(bi, sz) for bi in b]
    sols = []
    def rec(i, cx, cz, acc):
        if cx > X or (not sz and cz > Z): return
        if i == len(places):
            if cx == X and (sz or cz == Z): sols.append(tuple(acc))
            return
        for (x, z) in places[i]:
            acc.append((x, z)); rec(i + 1, cx + x, cz + (z or 0), acc); acc.pop()
    rec(0, 0, 0, []); return sols


def _t_z_T2(a): return {11: 3, 12: 6, 13: 9, 14: 12}.get(a)


def _enum_T2(a, b, Z):
    tz = _t_z_T2(a)
    if tz is None: return []
    per = [list(range(tz, Z + 1))]
    for bi in b:
        if bi == 0: per.append([0])
        elif bi == 1: per.append(list(range(2, Z + 1)))
        elif bi == 3: per.append([7])
        else: return []
    sols = []
    def rec(i, c, acc):
        if c > Z: return
        if i == len(per):
            if c == Z: sols.append(tuple(acc))
            return
        for z in per[i]: acc.append(z); rec(i + 1, c + z, acc); acc.pop()
    rec(0, 0, []); return sols


def _dv(x): return None if x in ('-inf', None) else int(x)


def load_targets():
    rc = json.load(open('alt_residue_congruences.json'))['states']
    ac = json.load(open('alt_combined.json'))['branches']
    flat = []
    for br in ac:
        for rs in br['remaining_states']:
            flat.append((br['id'], br['a'], tuple(br['b']), br['branch'],
                         br.get('sum_b'), rs['state'],
                         rs.get('finite_place_witness', {})))
    tgts = []
    for i, st in enumerate(rc):
        bid, a, b, branch, sum_b, state, fw = flat[i]
        dd1 = _dv(state['deg_d1']); dsig = _dv(state['deg_sigma'])
        dege = int(state['deg_e']); dd2 = _dv(state['deg_d2'])
        X = fw.get('X', fw.get('Xmin')); Z = fw.get('Z', fw.get('Zmin'))
        sz = (dsig is None)
        d1_0 = (branch == 'T2' and dd1 is None) or \
               (dd1 is not None and X is not None and dd1 - X == 0)
        sig_0 = sz or (Z is not None and dsig is not None and dsig - Z == 0)
        e_0 = (dege == a + sum_b)
        if not (d1_0 and sig_0 and e_0):
            continue
        splits = (_enum_T1(a, b, X, Z if not sz else 0, sz) if branch == 'T1'
                  else _enum_T2(a, b, Z or 0))
        assert len(splits) == 1, (bid, splits)
        tgts.append(dict(bid=bid, a=a, b=b, branch=branch, sum_b=sum_b,
                         degs=(dd2, dd1, dsig, dege), split=splits[0], sz=sz,
                         support=st['L0_tie_support_id'],
                         depth=st['L0_tie_depth'], idx=i))
    return tgts


def run_state(t, budget=90.0, max_depth=None):
    degs = t['degs']
    drop_d1 = (t['branch'] == 'T2')
    drop_sig = t['sz']
    TD = total_deg(degs, drop_d1, drop_sig)
    rec = reconstruct(t['a'], t['b'], t['split'], t['branch'], degs,
                      drop_d1, drop_sig)
    root_vars = rec[1]
    # cost guard: two marked roots + a large free d2 is a GB blowup even at
    # depth 2 (distinctness saturation x 2 root vars x free-d2 top coeffs).
    if len(root_vars) >= 2 and (degs[0] or 0) >= 6:
        return {'verdict': 'PENDING_HEAVY', 'reason': 'two-root + free deg-6 d2',
                'field': field_label(root_vars), 'TD': TD}
    # per-field depth policy: over-Q systems are cheap -> full tie depth;
    # single marked root -> moderate; two marked roots blow up -> depth 2.
    if len(root_vars) == 0:
        # over Q: allow deep tower, but large free-d2 GBs blow up -> cap at 9
        md = min(t['depth'], 2 if (degs[0] or 0) >= 6 else 12)
    elif len(root_vars) == 1:
        md = min(t['depth'], 2)
    else:
        md = 2
    if max_depth is not None:
        md = min(md, max_depth)
    res = kill_test(rec, TD, t['depth'], degs, drop_d1, drop_sig,
                    budget=budget, max_depth=md)
    res['TD'] = TD
    res['field'] = field_label(root_vars)
    res['capped_depth'] = md
    return res


def main():
    import os
    tgts = load_targets()
    print(f'{len(tgts)} alt defect-0 target states', flush=True)
    ckpt = 'phase_f2_scale.json'
    done = {}
    if os.path.exists(ckpt):
        try:
            done = {r['key']: r for r in json.load(open(ckpt)).get('alt_states', [])}
        except Exception:
            done = {}
    results = []
    budget = float(os.environ.get('F2_BUDGET', '70'))
    md = os.environ.get('F2_MAXDEPTH')
    md = int(md) if md else None
    for n, t in enumerate(tgts):
        key = f"{t['bid']}#sup{t['support']}#idx{t['idx']}"
        if key in done:
            results.append(done[key]); continue
        t0 = time.time()
        res = run_state(t, budget=budget, max_depth=md)
        res.update(key=key, bid=t['bid'], support=t['support'],
                   tie_depth=t['depth'], degs=list(t['degs']), branch=t['branch'])
        results.append(res)
        print(f"[{n+1}/{len(tgts)}] {key} deg{t['degs']} dep{t['depth']} "
              f"-> {res['verdict']} "
              f"(kd={res.get('kill_depth')}, {round(time.time()-t0,1)}s)",
              flush=True)
        json.dump({'alt_states': results}, open(ckpt, 'w'), indent=1)
    json.dump({'alt_states': results}, open(ckpt, 'w'), indent=1)
    from collections import Counter
    print('CENSUS', dict(Counter(r['verdict'] for r in results)))


if __name__ == '__main__':
    main()
