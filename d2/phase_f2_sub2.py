#!/usr/bin/env python3
"""Phase F, work item F2 -- SUB2 divisor-reconstruction kill test (extension).

PHASE_F2_SCALE.md killed 20 ALT states by divisor reconstruction but DEFERRED
the sub2 front's  b != 0000  forced states: "need a sub2 geometric-regime split
reconstruction that isn't recorded".  This file is that extension.

MECHANISM (the requested master-identity machinery, convolution_descent.py).
A sub2 state's polynomials live in the STANDARD regime with the master identity

    f31 = sum_{f=0}^{7} Phi^f * e^(21-3f) * h_f  ==  0 ,   Phi = c (y+1)^30 q,
    c = -1/6630 ,   q = 2048 y^4 - 512 y^3 + 320 y^2 - 240 y + 195 .

For a FULLY-FORCED state at cell (a, b, T1) whose core divisors are all defect 0
(deg p == sum_places v_place(p)) AND admit a UNIQUE simultaneous place split
supported on {t = y+1 ; the marked q-roots where b_j > 0}, every polynomial is
determined up to its single leading scalar:

    e     = E (y+1)^a       prod_j (y - r_j)^{b_j}        (defect-0 e; always)
    d1    = X (y+1)^{vt(d1)} prod_j (y - r_j)^{vj(d1)}
    sigma = S (y+1)^{vt(sig)} prod_j (y - r_j)^{vj(sig)}
    d2    = D (y+1)^{vt(d2)} prod_j (y - r_j)^{vj(d2)}    (IMPOSED when forced;
            record the choice; for a d2_zero cell d2 == 0).

Unlike the ALT lane (which kept d2 FREE, conservative), sub2 IMPOSES d2's
reconstruction when its defect is 0 -- recorded per state as "d2_mode".

For b = (1,0,0,0) exactly ONE root is marked; arithmetic is in Q[r]/(q(r)) with
q(r) = 0 adjoined to the ideal.  b = 0000 is over Q (no marked root; those states
coincide with the generic batch_convolution_sub2 degree tuples -- overlap noted).

THE TOWER.  With every polynomial determined up to scalars, feed them into the
convolution_descent master-identity engine and WALK the coefficients from the top
degree downward, accumulating them as an ideal, reducing mod q(r), saturating the
leading scalars nonzero (w * prod(scalars) - 1).  Verdict off an exact saturated
grevlex Groebner basis:  unit ideal  ==>  KILL.  Walk until unit ideal, a stall
(not the unit ideal at the tracked depth), or the time budget.

Verdicts are PENDING AUDIT.  READ-ONLY on every audited artifact; nothing here is
committed.
"""
import sys
import os
import json
import time
import multiprocessing as mp
from itertools import product

sys.path.insert(0, '.')
import sympy as sp
import cascade_engine as ce
import convolution_descent as cd

y = cd.y
r = sp.Symbol('r')
D, X, S, E, w = sp.symbols('D X S E w')
Q_COEFFS = [2048, -512, 320, -240, 195]
QR = sp.Poly(sum(c * r**(4 - i) for i, c in enumerate(Q_COEFFS)), r)
QR_EXPR = QR.as_expr()
DEPTH = 4                          # cascade descent depth (mirrors phase_f_defects)
C_VAL = sp.Rational(-1, 6630)


def deg_val(x):
    return None if x in ('-inf', None) else int(x)


# --------------------------------------------------------------------------
#  per-place valuation options (individual v(d1), v(sigma), v(d2)) from the
#  cascade engine Pareto profiles -- the exact machinery phase_f_defects used.
# --------------------------------------------------------------------------
def place_val_options(case):
    a = case['a_t']
    b_vec = tuple(case['b'])
    branch = case['branch']
    sz, dz = case['sigma_zero'], case['d2_zero']
    gz = tuple(case['g_zero_levels'])
    config = ce.CONFIGS['sub2']
    terminal = ce.T1_TERMINAL if branch == 'T1' else ce.T2_TERMINAL
    g_zero = {lv: (lv in gz) for lv in range(terminal - 1, DEPTH - 1, -1)}
    g_zero[terminal] = False
    r_cap = 10 + 3 * a

    def trip(p):
        return (0 if p.x == ce.INF else int(p.x),
                0 if p.z == ce.INF else int(p.z),
                0 if p.k == ce.INF else int(p.k))
    opts = []
    for b in b_vec:
        profs = ce.place_profiles(b, branch, r_cap, DEPTH, sz, dz, g_zero,
                                  config, a)
        opts.append(sorted(set(trip(p) for p in profs)))
    tp = ce.t_place_profiles(a, branch, r_cap, DEPTH, sz, dz, g_zero, config)
    opts.append(sorted(set(trip(p) for p in tp)))     # last entry = t place
    return opts


def unique_split(case, st, deltas=None):
    """Return (vt_d1,vr_d1..., etc) unique simultaneous split, or None.

    Returns dict place->(vd1,vsig,vd2) for places [root0..root3, t] iff exactly
    one simultaneous split reproduces the state degrees AND every polynomial's
    per-place divisor is unique.  Otherwise None (ambiguous / no split).
    """
    opts = place_val_options(case)
    branch = case['branch']
    sz, dz = case['sigma_zero'], case['d2_zero']
    # target the FORCED (defect-0) portion of each degree: deg - defect.  A
    # defect-d polynomial's divisor is forced only on deg-d of its degree; the
    # remaining d is a free cofactor added at reconstruction.
    dd = deltas or {}
    tx = None if branch == 'T2' else deg_val(st['deg_d1']) - dd.get('d1', 0)
    tz = None if sz else deg_val(st['deg_sigma']) - dd.get('sigma', 0)
    tk = None if dz else deg_val(st['deg_d2']) - dd.get('d2', 0)
    sols = []
    for combo in product(*opts):
        if tx is not None and sum(c[0] for c in combo) != tx:
            continue
        if tz is not None and sum(c[1] for c in combo) != tz:
            continue
        if tk is not None and sum(c[2] for c in combo) != tk:
            continue
        sols.append(combo)
    if not sols:
        return None, 'no_simultaneous_split', len(sols)
    d1div = set(tuple(c[0] for c in s) for s in sols)
    sigdiv = set(tuple(c[1] for c in s) for s in sols)
    d2div = set(tuple(c[2] for c in s) for s in sols)
    if len(d1div) > 1 or len(sigdiv) > 1 or len(d2div) > 1:
        return None, 'ambiguous_split', len(sols)
    combo = sols[0]
    return combo, 'unique', len(sols)


# --------------------------------------------------------------------------
#  reconstruction:  divisor -> determined polynomial (single leading scalar).
#  places are [root0,root1,root2,root3, t].  A marked root (b_j>0) uses the
#  single algebraic variable r (b=1000 => one marked root).
# --------------------------------------------------------------------------
U1, U2, U3 = sp.symbols('u1 u2 u3')      # linear-cofactor extra unknowns


def reconstruct(case, st, combo, deltas=None):
    a = case['a_t']
    b_vec = tuple(case['b'])
    branch = case['branch']
    sz, dz = case['sigma_zero'], case['d2_zero']
    dd = deltas or {}
    marked = [j for j in range(4) if b_vec[j] > 0]        # e-marked roots
    # per-place valuations
    vt_d1, vt_sig, vt_d2 = combo[4]
    root_val = {j: combo[j] for j in range(4)}            # (vd1,vsig,vd2)

    # every nonzero divisor valuation must sit on t or a marked root, else the
    # reconstruction would need an extra (unmarked, Galois) marked root.
    for j in range(4):
        if j in marked:
            continue
        if any(root_val[j]):
            return None, 'valuation_on_unmarked_root'
    if len(marked) > 1:
        return None, 'multi_marked_root_heavy'
    if dd.get('e', 0) != 0:
        return None, 'e_not_defect0'
    if any(dd.get(k, 0) >= 2 for k in ('d1', 'sigma', 'd2')):
        return None, 'defect_ge2_cofactor'

    rv = marked[0] if marked else None
    cofactors = []                # extra (unsaturated) unknowns

    def build(scalar, vt, vr, delta, u):
        rootpart = (y - r)**vr if rv is not None else sp.Integer(1)
        base_poly = scalar * (y + 1)**vt * rootpart
        if delta == 1:            # free linear cofactor: one extra unknown
            cofactors.append(u)
            base_poly = base_poly * (y - u)
        return sp.expand(base_poly)

    # e (defect-0, always): v_t(e)=a, v_{marked}(e)=b
    e_expr = build(E, a, b_vec[rv] if rv is not None else 0, 0, None)
    # d1
    if branch == 'T2':
        d1_expr = sp.Integer(0)
    else:
        d1_expr = build(X, vt_d1, root_val[rv][0] if rv is not None else 0,
                        dd.get('d1', 0), U1)
    # sigma
    if sz:
        sig_expr = sp.Integer(0)
    else:
        sig_expr = build(S, vt_sig, root_val[rv][1] if rv is not None else 0,
                         dd.get('sigma', 0), U2)
    # d2 : IMPOSE reconstruction when forced defect-0 (record); zero when flagged
    if dz:
        d2_expr = sp.Integer(0)
        d2_mode = 'zero_flag'
    else:
        d2_expr = build(D, vt_d2, root_val[rv][2] if rv is not None else 0,
                        dd.get('d2', 0), U3)
        d2_mode = ('reconstructed_defect0' if dd.get('d2', 0) == 0
                   else 'reconstructed_defect1_lincofactor')

    scalars = [E]
    if branch != 'T2':
        scalars.append(X)
    if not sz:
        scalars.append(S)
    if not dz:
        scalars.append(D)
    return (dict(d2=d2_expr, d1=d1_expr, sigma=sig_expr, e=e_expr),
            scalars, rv, d2_mode, cofactors)


# --------------------------------------------------------------------------
#  master-identity kill test
# --------------------------------------------------------------------------
def redq(expr, marked):
    if marked is None:
        return sp.expand(expr)
    return sp.rem(sp.Poly(sp.expand(expr), r), QR).as_expr()


def engine_top(eng):
    """Exact top degree of the master identity (max achievable target)."""
    top = None
    for f in range(8):
        eng.term_coefficient(f, 0)              # populate caches
        ph = eng._cached_power(eng._phi_powers, eng.phi, f)
        ep = eng._cached_power(eng._e_powers, eng.e_poly, 21 - 3 * f)
        hf = eng._h[f]
        if not ph or not ep or not hf:
            continue
        m = max(ph) + max(ep) + max(hf)
        top = m if top is None else max(top, m)
    return top


def kill_test(polys, scalars, marked, budget, max_coeffs, cofactors=()):
    params = (r,) if marked is not None else ()
    params = params + tuple(cofactors)          # cofactors are free, not forced
    ans = cd.build_ansatz(d2=polys['d2'], d1=polys['d1'], e=polys['e'],
                          sigma=polys['sigma'], parameters=params)
    eng = cd.ConvolutionDescent(ans, c=C_VAL)
    top = engine_top(eng)
    order_vars = (list(scalars) + list(cofactors)
                  + ([r] if marked is not None else []) + [w])
    sat = w * sp.prod(scalars) - 1              # only leading scalars saturated
    extra = [QR_EXPR] if marked is not None else []
    gens = []
    j_top = None
    t0 = time.time()
    for n in range(max_coeffs):
        target = top - n
        mc = redq(eng.master_coefficient(target), marked)
        if j_top is None and mc != 0:
            j_top = str(sp.factor(mc)) if marked is None else str(mc)
        if mc != 0:
            gens.append(mc)
        if not gens:
            continue
        try:
            G = sp.groebner(gens + extra + [sat], *order_vars, order='grevlex')
        except Exception as ex:
            return {'verdict': 'ERROR', 'error': str(ex)[:200],
                    'depth': n + 1, 'top_degree': top}
        if sp.Integer(1) in G.exprs:
            return {'verdict': 'KILLED', 'kill_depth': n + 1, 'top_degree': top,
                    'ngens': len(gens), 'elapsed': round(time.time() - t0, 1),
                    'j_top': j_top}
        if time.time() - t0 > budget:
            return {'verdict': 'PENDING_TIMEOUT', 'depth': n + 1,
                    'top_degree': top, 'elapsed': round(time.time() - t0, 1),
                    'j_top': j_top}
    return {'verdict': 'NARROWED', 'depth': max_coeffs, 'top_degree': top,
            'ngens': len(gens), 'elapsed': round(time.time() - t0, 1),
            'j_top': j_top}


def _worker(q, polys, scalars, marked, budget, max_coeffs, cofactors):
    try:
        q.put(kill_test(polys, scalars, marked, budget, max_coeffs, cofactors))
    except Exception as ex:                    # pragma: no cover
        q.put({'verdict': 'ERROR', 'error': str(ex)[:200]})


def kill_test_guarded(polys, scalars, marked, budget, max_coeffs, hard,
                      cofactors=()):
    """Run kill_test in a spawned process; terminate if it exceeds `hard` s.

    A single grevlex Groebner call can outrun the internal per-coefficient
    budget (a high-degree marked-root GB blowup); this hard wall keeps such
    states as PENDING_HARD_TIMEOUT instead of stalling the census.
    """
    ctx = mp.get_context('spawn')
    q = ctx.Queue()
    p = ctx.Process(target=_worker,
                    args=(q, polys, scalars, marked, budget, max_coeffs,
                          tuple(cofactors)))
    p.start()
    p.join(hard)
    if p.is_alive():
        p.terminate()
        p.join()
        return {'verdict': 'PENDING_HARD_TIMEOUT', 'hard_timeout': hard}
    try:
        return q.get_nowait()
    except Exception:
        return {'verdict': 'ERROR', 'error': 'worker produced no result'}


# --------------------------------------------------------------------------
#  target enumeration
# --------------------------------------------------------------------------
def batch_kill_sigs():
    sigs = {}
    for fn in ('batch_convolution_sub2.json', 'batch_convolution_sub2_round2.json',
               'batch_convolution_overnight.json'):
        if not os.path.exists(fn):
            continue
        d = json.load(open(fn))
        for k in d.get('kills_pending_audit', []):
            s = (k['a_t'], k['branch'], str(k['deg_d1']), str(k['deg_d2']),
                 str(k['deg_sigma']), str(k['deg_e']), bool(k.get('d2_zero')),
                 bool(k.get('sigma_zero')), bool(k.get('d1_zero')))
            sigs.setdefault(s, fn)
    return sigs


def state_sig(case, st):
    return (case['a_t'], case['branch'], str(st['deg_d1']), str(st['deg_d2']),
            str(st['deg_sigma']), str(st['deg_e']), bool(case['d2_zero']),
            bool(case['sigma_zero']), case['branch'] == 'T2')


def core_deltas(delta):
    core = ['d1', 'sigma', 'd2', 'e', 'g7']
    return {k: delta[k] for k in core if k in delta}


def cell_matches(case, targets):
    a = case['a_t']
    b = ''.join(map(str, case['b']))
    br = case['branch']
    dz = case['d2_zero']
    for (ta, tb, tbr, tdz) in targets:
        if a == ta and b == tb and br == tbr and (tdz is None or dz == tdz):
            return True
    return False


def load_targets(target_cells, max_defect=0):
    data = json.load(open('phase_d_states_sub2.json'))
    defects = json.load(open('phase_f_defects.json'))['windows']['sub2']['cases']
    dmap = {c['cellid']: c for c in defects}
    out = []
    for case in data['cases']:
        if not cell_matches(case, target_cells):
            continue
        cellid = (f"sub2:a{case['a_t']}_b{''.join(map(str, case['b']))}"
                  f"_{case['branch']}_sz{int(case['sigma_zero'])}"
                  f"_dz{int(case['d2_zero'])}"
                  f"_gz{'.'.join(map(str, case['g_zero_levels'])) or '-'}")
        deltas = dmap[cellid]['state_deltas']
        for i, st in enumerate(case['states']):
            d = core_deltas(deltas[i])
            if '_error' in deltas[i]:
                continue
            mx = max(d.values()) if d else 0
            if mx > max_defect:
                continue
            # reconstruction-relevant deltas: only d1,sigma,d2 (the master
            # identity ignores the g-chain).  e must be defect 0 (checked later).
            pdelta = {k: deltas[i].get(k, 0) for k in ('d1', 'sigma', 'd2', 'e')
                      if k in deltas[i]}
            out.append((cellid, case, st, i, mx, pdelta))
    return out


TARGET_CELLS = [
    (9, '1000', 'T1', False),   # a9_b1000_T1 dz0  (marked single root)
    (9, '1000', 'T1', True),    # a9_b1000_T1 dz1  (d2 == 0)
    (10, '0000', 'T1', False),  # a10_b0000_T1 dz0 (over Q; batch overlap)
    (9, '0000', 'T1', False),   # a9_b0000_T1  dz0 (over Q; batch overlap)
]


def main():
    budget = float(os.environ.get('F2SUB2_BUDGET', '30'))
    max_coeffs = int(os.environ.get('F2SUB2_MAXCOEFFS', '10'))
    max_defect = int(os.environ.get('F2SUB2_MAXDEFECT', '0'))
    ckpt = 'phase_f2_sub2.json'
    bsigs = batch_kill_sigs()
    targets = load_targets(TARGET_CELLS, max_defect=max_defect)
    print(f'{len(targets)} candidate target states (max_defect={max_defect})',
          flush=True)
    done = {}
    if os.path.exists(ckpt):
        try:
            done = {r['key']: r for r in json.load(open(ckpt)).get('states', [])}
        except Exception:
            done = {}
    results = []
    for n, (cellid, case, st, idx, mx, pdelta) in enumerate(targets):
        key = f'{cellid}#state{idx}'
        if key in done:
            results.append(done[key])
            continue
        rec = {'key': key, 'cellid': cellid, 'state_idx': idx,
               'degs': [st['deg_d1'], st['deg_sigma'], st['deg_d2'],
                        st['deg_e']], 'max_core_defect': mx,
               'poly_deltas': pdelta}
        sig = state_sig(case, st)
        rec['batch_overlap'] = bsigs.get(sig)
        combo, why, nsol = unique_split(case, st, pdelta)
        rec['split_status'] = why
        rec['n_simultaneous_splits'] = nsol
        if combo is None:
            rec['verdict'] = 'SKIPPED_' + why.upper()
            results.append(rec)
            _log(n, targets, rec)
            json.dump({'states': results}, open(ckpt, 'w'), indent=1)
            continue
        recon = reconstruct(case, st, combo, pdelta)
        if recon[0] is None:
            rec['verdict'] = 'SKIPPED_' + recon[1].upper()
            results.append(rec)
            _log(n, targets, rec)
            json.dump({'states': results}, open(ckpt, 'w'), indent=1)
            continue
        polys, scalars, marked, d2_mode, cofactors = recon
        rec['d2_mode'] = d2_mode
        rec['n_cofactors'] = len(cofactors)
        rec['field'] = 'Q' if marked is None else 'Q[r]/(q)'
        rec['scalars'] = [str(s) for s in scalars]
        hard = float(os.environ.get('F2SUB2_HARD', '40'))
        try:
            res = kill_test_guarded(polys, scalars, marked, budget,
                                    max_coeffs, hard, cofactors=cofactors)
        except Exception as ex:
            res = {'verdict': 'ERROR', 'error': str(ex)[:200]}
        rec.update(res)
        results.append(rec)
        _log(n, targets, rec)
        json.dump({'states': results}, open(ckpt, 'w'), indent=1)
    json.dump({'states': results}, open(ckpt, 'w'), indent=1)
    from collections import Counter
    print('CENSUS', dict(Counter(r['verdict'] for r in results)), flush=True)


def _log(n, targets, rec):
    print(f"[{n+1}/{len(targets)}] {rec['key']} deg{rec['degs']} "
          f"def{rec['max_core_defect']} -> {rec['verdict']} "
          f"(kd={rec.get('kill_depth')}, field={rec.get('field')}, "
          f"top={rec.get('top_degree')}, {rec.get('elapsed')}s)", flush=True)


if __name__ == '__main__':
    main()
