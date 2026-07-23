#!/usr/bin/env python3
"""MODULAR TRIAGE (reconnaissance, NOT proof).

For each resistant ideal system we build the exact generators in sympy, reduce
coefficients modulo several good primes p, and compute a Groebner basis over
F_p in Singular (WSL).  Per prime the verdict is UNIT (basis = 1, i.e. the
variety is empty over \bar F_p) or PROPER (basis != 1, dim recorded).

PREDICTIONS (labelled as such, never proofs):
  LIKELY-EMPTY     : UNIT for every prime.
  LIKELY-SOLVABLE  : PROPER for every prime with consistent dimension.
  MIXED            : disagreement across primes (bad prime or genuine subtlety).

IMPORTANT CAVEAT baked into the interpretation: a mod-p Groebner basis sees
emptiness over the ALGEBRAIC CLOSURE of F_p.  A system that is empty only over
R (a real/positive-definiteness obstruction) will still be PROPER mod p.  Such
systems are flagged; "LIKELY-SOLVABLE" here means "has a solution over some
field", not necessarily over R.

New file, uncommitted.  READ-ONLY on every imported module/artifact.
"""
from __future__ import annotations

import json
import re
import subprocess
import time
from typing import Iterable

import sympy as sp

import convolution_descent as cd
import convolution_elim as ce
import convolution_elim_qsupport as qs

Y = cd.y
# primes chosen with >= 2 distinct roots of q so two-marked-root systems can be
# specialized; all avoid the only bad primes {2,3,5,13,17} (disc/lc/6630).
PRIMES = [10007, 10009, 100019]
Q_COEFFS = [2048, -512, 320, -240, 195]
WSL = ('wsl.exe', '-d', 'Ubuntu', '--', 'bash', '-lc', 'cd $HOME && Singular -q')
SING_TIMEOUT = 60.0


# --------------------------------------------------------------------------
#  mod-p reduction of an exact (rational-coefficient) sympy polynomial to a
#  Singular integer-coefficient string.  Denominators must be invertible mod p.
# --------------------------------------------------------------------------
def poly_to_singular_modp(expr: sp.Expr, gens: list[sp.Symbol], p: int) -> str:
    expr = sp.cancel(sp.sympify(expr))
    if expr == 0:
        return '0'
    num, den = sp.fraction(expr)
    poly = sp.Poly(sp.expand(num), *gens)
    # a single global denominator (rational constant) is typical after cancel
    den_poly = sp.Poly(sp.expand(den), *gens)
    if den_poly.total_degree() != 0:
        raise ValueError(f'non-constant denominator: {den}')
    den_val = int(den_poly.coeff_monomial(1))
    if den_val % p == 0:
        raise ValueError(f'denominator {den_val} divisible by prime {p}')
    den_inv = pow(den_val % p, -1, p)
    terms = []
    for monom, coeff in poly.terms():
        c = int(coeff)
        if c % p == 0:
            continue
        cm = (c * den_inv) % p
        factors = [str(cm)] if (cm != 1 or all(e == 0 for e in monom)) else []
        for g, e in zip(gens, monom):
            if e == 1:
                factors.append(g.name)
            elif e > 1:
                factors.append(f'{g.name}^{e}')
        terms.append('*'.join(factors) if factors else str(cm))
    return ('+'.join(terms)).replace('+-', '-') or '0'


def q_roots_mod_p(p: int, n: int) -> list[int]:
    """Return up to n distinct roots of q modulo p."""
    yv = sp.Symbol('y')
    qp = sp.Poly(sum(c * yv**(4 - i) for i, c in enumerate(Q_COEFFS)),
                 yv, modulus=p)
    return [int(r) % p for r in qp.ground_roots()][:n]


def build_singular_program(gens: list[sp.Expr], ring_vars: list[sp.Symbol],
                           p: int, want_dim: bool = True,
                           sat_factors: list[sp.Expr] | None = None) -> str:
    var_txt = ','.join(v.name for v in ring_vars)
    lines = ['LIB "elim.lib";', f'ring R = {p},({var_txt}),dp;']
    members = []
    for i, g in enumerate(gens):
        s = poly_to_singular_modp(g, ring_vars, p)
        if s == '0':
            continue
        lines.append(f'poly g{i} = {s};')
        members.append(f'g{i}')
    if not members:
        members = ['0']
    lines.append(f'ideal I = {",".join(members)};')
    if sat_factors:
        nz = '*'.join('(' + poly_to_singular_modp(f, ring_vars, p) + ')'
                      for f in sat_factors if poly_to_singular_modp(f, ring_vars, p) != '0')
        lines.append(f'poly nz = {nz if nz else "1"};')
        lines.append('ideal Is = sat(I,nz)[1];')
        lines.append('ideal G = std(Is);')
    else:
        lines.append('ideal G = std(I);')
    lines.append('int u = (reduce(1,G)==0);')
    lines.append('"@@UNIT";')
    lines.append('u;')
    if want_dim:
        lines.append('"@@DIM";')
        lines.append('dim(G);')
    lines.append('quit;')
    return '\n'.join(lines) + '\n'


def run_singular(program: str, timeout: float = SING_TIMEOUT) -> dict:
    t0 = time.monotonic()
    try:
        cp = subprocess.run(WSL, input=program, text=True, encoding='utf-8',
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            timeout=timeout, check=False)
        out = (cp.stdout or '').replace('\x00', '')
        err = (cp.stderr or '').replace('\x00', '')
        status = 'ok' if cp.returncode == 0 else 'proc_error'
    except subprocess.TimeoutExpired:
        return {'status': 'timeout', 'verdict': 'TIMEOUT',
                'unit': None, 'dim': None, 'wall': round(time.monotonic() - t0, 2)}
    combined = out + '\n' + err
    um = re.search(r'@@UNIT\s*\r?\n\s*(-?\d+)', combined)
    dm = re.search(r'@@DIM\s*\r?\n\s*(-?\d+)', combined)
    unit = None if um is None else bool(int(um.group(1)))
    dim = None if dm is None else int(dm.group(1))
    verdict = 'UNIT' if unit else ('PROPER' if unit is False else 'PARSE_FAIL')
    return {'status': status, 'verdict': verdict, 'unit': unit, 'dim': dim,
            'wall': round(time.monotonic() - t0, 2),
            'stderr': err.strip()[:400] if status != 'ok' else ''}


def ring_vars_of(gens: Iterable[sp.Expr], extra: Iterable[sp.Symbol] = ()) -> list[sp.Symbol]:
    syms: set[sp.Symbol] = set(extra)
    for g in gens:
        syms |= sp.sympify(g).free_symbols
    syms.discard(Y)
    return sorted(syms, key=lambda s: s.name)


# ==========================================================================
#  SYSTEM 1 : R9 z=0..6 q-supported ideals
# ==========================================================================
def build_system1() -> list[dict]:
    subsystems = []
    ncoef = 8  # more generators over-determine the unit ideal -> collapses faster
    for z in range(7):
        state = qs.build_qsupport_ansatz(z)
        engine = ce.HighCoefficientEngine(
            state.ansatz, start_degree=qs.START_DEGREE, target_count=ncoef,
            c=ce.DEFAULT_C)
        degrees = list(range(qs.START_DEGREE, qs.START_DEGREE - ncoef, -1))
        coeffs = [qs.quotient_reduce(engine.master_coefficient(d)) for d in degrees]
        coeffs = [c for c in coeffs if c != 0]
        Gr = qs.quotient_reduce(state.G.subs(qs.y, qs.r))
        sat_factors = [qs.gamma, state.g_coefficients[-1], Gr]
        gens = list(coeffs)  # marked root r specialized numerically per prime
        rv = ring_vars_of(gens + sat_factors)
        subsystems.append({'name': f'R9_z{z}', 'gens': gens, 'ring_vars': rv,
                           'sat_factors': sat_factors, 'root_syms': [qs.r],
                           'note': f'z={z}, {len(coeffs)} master coeffs, sat(gamma,g_last,G(r)), r->root'})
    return subsystems


# ==========================================================================
#  SYSTEM 2 : a8 constant-E gauge stall states (24)
# ==========================================================================
def build_system2() -> list[dict]:
    d = json.load(open('batch_convolution_sub2.json'))
    states = [s for s in d['states']
              if s['a_t'] == 8 and s['branch'] == 'T1'
              and s['deg_e'] == 8 and s['final_verdict'] == 'UNRESOLVED']
    gamma = sp.Symbol('gamma')
    out = []
    for s in states:
        dd1, dsig = int(s['deg_d1']), int(s['deg_sigma'])
        gd = s['gauge_detail']
        start = int(gd['start'])
        stop = int(gd['stopping_degree'])
        e = gamma * (Y + 1)**8
        degrees = {'d1': dd1, 'sigma': dsig}
        if s['d2_zero']:
            d2arg = {'d2': sp.Integer(0)}
        else:
            degrees['d2'] = int(s['deg_d2'])
            d2arg = {}
        ansatz = cd.build_ansatz(e=e, degrees=degrees, parameters=(gamma,),
                                 **d2arg)
        eng = cd.ConvolutionDescent(ansatz, c=ce.DEFAULT_C)
        gens = []
        # accumulate nonzero master coefficients across the recorded window
        for dgr in range(start, stop - 3, -1):
            c = sp.expand(eng.master_coefficient(dgr))
            if c != 0:
                gens.append(c)
            if len(gens) >= 3:
                break
        w = sp.Symbol('w')
        gens_full = list(gens) + [w * gamma - 1]
        rv = ring_vars_of(gens_full, extra=[w, gamma])
        name = f"a8_dd2{s['deg_d2']}_dd1{dd1}_dsig{dsig}"
        out.append({'name': name, 'gens': gens_full, 'ring_vars': rv,
                    'note': f'{len(gens)} master coeffs (gauge, gamma sat)'})
    return out


# ==========================================================================
#  SYSTEM 3 : alt NARROWED/UNOBSTRUCTED reconstruction tie-towers (18)
# ==========================================================================
def build_system3() -> list[dict]:
    import phase_f2_scale as f2
    narrowed_keys = set()
    js = json.load(open('phase_f2_scale.json'))['alt_states']
    key_meta = {}
    for rrec in js:
        if rrec['verdict'].startswith('NARROWED'):
            narrowed_keys.add(rrec['key'])
            key_meta[rrec['key']] = rrec
    tgts = f2.load_targets()
    out = []
    for t in tgts:
        key = f"{t['bid']}#sup{t['support']}#idx{t['idx']}"
        if key not in narrowed_keys:
            continue
        degs = t['degs']
        drop_d1 = (t['branch'] == 'T2')
        drop_sig = t['sz']
        TD = f2.total_deg(degs, drop_d1, drop_sig)
        rec = f2.reconstruct(t['a'], t['b'], t['split'], t['branch'], degs,
                             drop_d1, drop_sig)
        factors, root_vars, scalars, Dc = rec
        # push depth deeper than the rational cap: mod-p GB is cheap
        depth = min(int(t['depth']), 12)
        degs_num = tuple(x if x is not None else 0 for x in degs)
        C = f2.h0_top(factors, degs_num, TD, depth,
                      drop_d1=drop_d1, drop_sig=drop_sig)
        # assemble the Rabinowitsch-saturated generator list (over F_p keep r
        # as a ring variable and adjoin q(r); do NOT pre-reduce)
        w = sp.Symbol('w')
        sat_scalars = list(scalars)
        d2_needed = f2.d2_in_window(degs, TD, depth, drop_d1, drop_sig)
        d2_ring = []
        if Dc is not None and d2_needed:
            d2_ring = Dc[max(0, len(Dc) - depth):]
            sat_scalars = sat_scalars + [Dc[-1]]
        elif Dc is not None and degs[0] == 0:
            d2_ring = list(Dc)
        sat = w * sp.prod(sat_scalars) - 1
        gens = [c for c in C if c != 0]
        root_gens = []
        for rv_ in (f2.r, f2.r1, f2.r2):
            if rv_ in root_vars:
                root_gens.append(f2.qpoly(rv_))
        distinct = []
        extra_vars = [w]
        if f2.r1 in root_vars and f2.r2 in root_vars:
            wd = sp.Symbol('wd')
            distinct = [wd * (f2.r1 - f2.r2) - 1]
            extra_vars.append(wd)
        all_gens = gens + [sat]  # root_gens/distinct auto-satisfied when r->root
        rvars = ring_vars_of(all_gens, extra=list(scalars) + list(d2_ring)
                             + extra_vars)
        out.append({'name': key.split('#')[0] + '_' + key.split('#')[1],
                    'gens': all_gens, 'ring_vars': rvars,
                    'root_syms': list(root_vars),
                    'note': f"field={f2.field_label(root_vars)}, "
                            f"depth={depth}, {len(gens)} tie coeffs"})
    return out


# ==========================================================================
#  SYSTEM 4 : sub2 T2 pattern-B tie states (sample of 10, a7/a8 cells)
#  q-support ansatz analogous to R9 (single / no marked root only).
# ==========================================================================
def _quotient_reduce_r(expr, rsym):
    expr = sp.cancel(sp.sympify(expr))
    if expr == 0 or not expr.has(rsym):
        return expr
    num, den = sp.fraction(expr)
    if den.has(rsym):
        raise ValueError('r in denominator')
    QR = sum(c * rsym**(4 - i) for i, c in enumerate(Q_COEFFS))
    rem = sp.rem(sp.Poly(sp.expand(num), rsym), sp.Poly(QR, rsym))
    return sp.cancel(rem.as_expr() / den)


def build_system4() -> list[dict]:
    d = json.load(open('phase_d_states_sub2.json'))
    r = sp.Symbol('r')
    gamma = sp.Symbol('gamma')
    QR = sum(c * r**(4 - i) for i, c in enumerate(Q_COEFFS))
    # single/no-marked-root T2 cells: (a_t, b) -> we sample states within.
    wanted = [(8, [0, 0, 0, 0]), (8, [1, 0, 0, 0]),
              (7, [1, 0, 0, 0]), (7, [3, 0, 0, 0])]
    # pick per-cell how many states to sample to total 10
    quota = {(8, (0, 0, 0, 0)): 2, (8, (1, 0, 0, 0)): 3,
             (7, (1, 0, 0, 0)): 3, (7, (3, 0, 0, 0)): 2}
    out = []
    seen_global = set()
    taken_by_key = {k: 0 for k in quota}
    for case in d['cases']:
        if case['branch'] != 'T2':
            continue
        key = (case['a_t'], tuple(case['b']))
        if key not in quota or case['g_zero_levels']:
            continue
        a_t = case['a_t']
        b = case['b']
        sumb = sum(b)
        active = [i for i in range(4) if b[i] > 0]
        n = quota[key]
        seen = seen_global
        for st in case['states']:
            if taken_by_key[key] >= n:
                break
            dd2 = st['deg_d2']
            dsig = st['deg_sigma']
            if dsig in ('-inf', None):
                continue
            dsig = int(dsig)
            deg_G = dsig - 2 * sumb
            if deg_G < 0:
                continue
            sig = (a_t, tuple(b), dd2, dsig)
            if sig in seen:
                continue
            seen.add(sig)
            # build e and sigma with single marked root r (or none)
            root_factor = (Y - r)**sumb if active else sp.Integer(1)
            e = gamma * (Y + 1)**a_t * root_factor
            Gcoeffs = tuple(sp.symbols(f'g0:{deg_G + 1}'))
            G = sum(cc * Y**i for i, cc in enumerate(Gcoeffs))
            sig_root = (Y - r)**(2 * sumb) if active else sp.Integer(1)
            sigma = sp.expand(sig_root * G)
            dd2i = None if dd2 in ('-inf', None) else int(dd2)
            if dd2i is None:
                ansatz = cd.build_ansatz(d2=sp.Integer(0), d1=sp.Integer(0),
                                         e=e, sigma=sigma)
            else:
                ansatz = cd.build_ansatz(degrees={'d2': dd2i}, d1=sp.Integer(0),
                                         e=e, sigma=sigma, prefixes={'d2': 'a'})
            eng = cd.ConvolutionDescent(ansatz, c=ce.DEFAULT_C)
            # adaptive downward scan: zero coefficients above the top master
            # degree are cheap (empty convolution); collect first 6 nonzero.
            coeffs = []
            for dgr in range(260, 200, -1):
                c = _quotient_reduce_r(eng.master_coefficient(dgr), r)
                if c != 0:
                    coeffs.append(c)
                if len(coeffs) >= 6:
                    break
            w = sp.Symbol('w')
            sat_factors = [gamma]
            if Gcoeffs:
                sat_factors.append(Gcoeffs[-1])
            if active:
                Gr = _quotient_reduce_r(G.subs(Y, r), r)
                sat_factors.append(Gr)
            sat_prod = sp.Integer(1)
            for f in sat_factors:
                sat_prod *= f
            sat_gen = _quotient_reduce_r(w * sat_prod - 1, r)
            gens = list(coeffs) + [sat_gen]  # r specialized numerically per prime
            rv = ring_vars_of(gens, extra=[w, gamma])
            bstr = ''.join(map(str, b))
            name = f"sub2T2_a{a_t}_b{bstr}_dd2{dd2}_dsig{dsig}"
            out.append({'name': name, 'gens': gens, 'ring_vars': rv,
                        'root_syms': [r] if active else [],
                        'note': f'{len(coeffs)} master coeffs, '
                                f'{"1 root" if active else "no root"}, '
                                f'deg_G={deg_G}'})
            taken_by_key[key] += 1
    return out


# ==========================================================================
#  RUNNER
# ==========================================================================
def classify(prime_verdicts: list[dict]) -> str:
    vs = [pv['verdict'] for pv in prime_verdicts if pv.get('verdict') in ('UNIT', 'PROPER')]
    if not vs:
        return 'INDETERMINATE'
    if all(v == 'UNIT' for v in vs):
        return 'LIKELY-EMPTY'
    if all(v == 'PROPER' for v in vs):
        dims = {pv.get('dim') for pv in prime_verdicts if pv['verdict'] == 'PROPER'}
        return 'LIKELY-SOLVABLE' if len(dims) == 1 else 'LIKELY-SOLVABLE(dim-vary)'
    return 'MIXED'


def run_all(system_filter=None):
    builders = {'1': build_system1, '2': build_system2,
                '3': build_system3, '4': build_system4}
    result = {}
    for sysid, builder in builders.items():
        if system_filter and sysid not in system_filter:
            continue
        print(f'=== building system {sysid} ===', flush=True)
        t0 = time.monotonic()
        subs = builder()
        print(f'  {len(subs)} subsystems built in {time.monotonic()-t0:.1f}s', flush=True)
        sysres = []
        for ss in subs:
            pvs = []
            root_syms = ss.get('root_syms', [])
            base_sat = ss.get('sat_factors')
            timed_out = False
            for p in PRIMES:
                if timed_out:
                    pvs.append({'status': 'skip', 'verdict': 'TIMEOUT',
                                'reason': 'skipped after earlier timeout',
                                'unit': None, 'dim': None, 'prime': p})
                    continue
                try:
                    roots = q_roots_mod_p(p, len(root_syms))
                    if len(roots) < len(root_syms):
                        rr = {'status': 'skip', 'verdict': 'SKIP',
                              'reason': f'need {len(root_syms)} distinct roots, '
                                        f'q has {len(roots)} mod {p}',
                              'unit': None, 'dim': None}
                    else:
                        subst = {s: v for s, v in zip(root_syms, roots)}
                        gens = [sp.sympify(g).subs(subst) for g in ss['gens']]
                        satf = ([sp.sympify(f).subs(subst) for f in base_sat]
                                if base_sat else None)
                        rvars = ring_vars_of(gens + (satf or []))
                        prog = build_singular_program(gens, rvars, p,
                                                      sat_factors=satf)
                        rr = run_singular(prog)
                except Exception as ex:
                    rr = {'status': 'build_error', 'verdict': 'ERROR',
                          'error': str(ex)[:200], 'unit': None, 'dim': None}
                rr['prime'] = p
                if rr.get('verdict') == 'TIMEOUT':
                    timed_out = True
                pvs.append(rr)
            pred = classify(pvs)
            rec = {'name': ss['name'], 'note': ss['note'],
                   'nvars': len(ss['ring_vars']),
                   'primes': pvs, 'prediction': pred}
            sysres.append(rec)
            vstr = ','.join(f"{pv['prime']}:{pv['verdict']}" for pv in pvs)
            print(f"  {ss['name']}: {vstr} -> {pred}", flush=True)
        result[f'system{sysid}'] = sysres
        json.dump(result, open('modular_triage.json', 'w'), indent=1)
    json.dump(result, open('modular_triage.json', 'w'), indent=1)
    return result


if __name__ == '__main__':
    import sys
    which = sys.argv[1] if len(sys.argv) > 1 else 'all'
    if which == 'test':
        ss = build_system1()[0]
        prog = build_singular_program(ss['gens'], ss['ring_vars'], PRIMES[0])
        print('vars:', [v.name for v in ss['ring_vars']])
        print(run_singular(prog))
    else:
        filt = None if which == 'all' else set(which)
        run_all(filt)
