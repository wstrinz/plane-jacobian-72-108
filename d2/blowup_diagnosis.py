#!/usr/bin/env python3
"""BLOWUP DIAGNOSIS: swell vs. structural characterization of resistant GB cases.

For each known Groebner-blowup case we answer the owner's question --
"would N more GB of RAM finish this, and what is N -- or is RAM not the fix?" --
via the protocol:

  (1) mod-p Groebner basis (p=10007) with a GENEROUS budget, recording wall time,
      Singular peak heap (memory(2)) and process peak RSS (/usr/bin/time -v).
  (2) if mod-p finishes fast but the rational computation blows up  -> SWELL
        (RAM / multi-modular lifting is the fix; estimate N).
      if mod-p ALSO struggles                                        -> STRUCTURAL
        (RAM is not the fix; needs reformulation: resultant pre-elim,
         weight/lex orders, variable reduction).
  (3) for the rational side, run Singular over Q under increasing `ulimit -v`
      caps to get an honest peak-RSS / "N GB" estimate (or "unbounded-looking").

This driver REUSES the audited-adjacent builders read-only:
  modular_triage.build_system1  (R9 z=0..6)
  modular_triage.build_system3  (alt reconstruction tie towers, incl a11_b1111_T1)
  d2_threshold                  (deg-d2 in {5,6} alt tie towers, free d2)
  phase_f2_sub2                 (sub2 divisor-reconstruction master-identity towers)

NEW file, uncommitted, READ-ONLY on every imported module/artifact.
Usage:  python blowup_diagnosis.py <case>          # run one case's mod-p probe
        python blowup_diagnosis.py ratcurve <bid> <degd2> [caps_gb...]
        python blowup_diagnosis.py list
"""
from __future__ import annotations
import json, os, re, subprocess, sys, time
import sympy as sp

import modular_triage as mt

Q_COEFFS = [2048, -512, 320, -240, 195]
BIGP = 10007
HOME_WSL = ('wsl.exe', '-d', 'Ubuntu', '--', 'bash', '-lc')


# --------------------------------------------------------------------------
#  Singular runner with time + memory capture (mod-p or over-Q)
# --------------------------------------------------------------------------
def run_singular_measured(program: str, timeout: float, vcap_gb: float | None = None) -> dict:
    """Run a Singular program via WSL under /usr/bin/time -v, optional ulimit -v.
    The program should itself print @@UNIT / @@DIM / @@MEM markers."""
    ulimit = f'ulimit -v {int(vcap_gb * 1024 * 1024)}; ' if vcap_gb else ''
    inner = f'cd $HOME && {ulimit}/usr/bin/time -v Singular -q'
    cmd = HOME_WSL + (inner,)
    t0 = time.monotonic()
    try:
        cp = subprocess.run(cmd, input=program, text=True, encoding='utf-8',
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            timeout=timeout, check=False)
    except subprocess.TimeoutExpired:
        return {'verdict': 'TIMEOUT', 'wall': round(time.monotonic() - t0, 1),
                'unit': None, 'dim': None, 'rss_kb': None, 'sing_mem_mb': None}
    out = (cp.stdout or '').replace('\x00', '')
    err = (cp.stderr or '').replace('\x00', '')
    comb = out + '\n' + err
    um = re.search(r'@@UNIT\s*\r?\n\s*(-?\d+)', comb)
    dm = re.search(r'@@DIM\s*\r?\n\s*(-?\d+)', comb)
    mm = re.search(r'@@MEM\s*\r?\n\s*(-?\d+)', comb)
    rss = re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)', err)
    unit = None if um is None else bool(int(um.group(1)))
    # detect OOM / ulimit kill
    killed = ('std::bad_alloc' in comb or 'Cannot allocate' in comb
              or 'out of memory' in comb.lower() or cp.returncode == 137
              or (rss is None and um is None and dm is None))
    verdict = ('UNIT' if unit else ('PROPER' if unit is False else
               ('OOM' if killed else 'PARSE_FAIL')))
    return {'verdict': verdict, 'wall': round(time.monotonic() - t0, 1),
            'unit': unit, 'dim': (int(dm.group(1)) if dm else None),
            'rss_kb': (int(rss.group(1)) if rss else None),
            'sing_mem_mb': (round(int(mm.group(1)) / 1048576, 1) if mm else None),
            'rc': cp.returncode, 'tail': comb.strip()[-300:] if verdict in
            ('OOM', 'PARSE_FAIL') else ''}


def q_root_mod_p(p, n=1):
    yv = sp.Symbol('y')
    qp = sp.Poly(sum(c * yv**(4 - i) for i, c in enumerate(Q_COEFFS)), yv, modulus=p)
    return [int(x) % p for x in qp.ground_roots()][:n]


def emit_program(gens, ring_vars, char, sat_factors=None, minpoly_gens=None,
                 want_dim=True):
    """Build a Singular program (char=p for mod-p, char=0 for over-Q).
    minpoly_gens: extra ideal generators (e.g. q(r)) kept in the ideal for
    number-field arithmetic over Q; over F_p roots are specialized numerically."""
    var_txt = ','.join(v.name for v in ring_vars)
    if char == 0:
        conv = lambda g: mt.poly_to_singular_modp(g, ring_vars, 2147483647)  # placeholder unused
    lines = ['LIB "elim.lib";', f'ring R = {char},({var_txt}),dp;']
    members = []
    def poly_str(g):
        if char == 0:
            return sing_poly_intcoeff(g, ring_vars)
        return mt.poly_to_singular_modp(g, ring_vars, char)
    idx = 0
    for g in list(gens) + list(minpoly_gens or []):
        s = poly_str(g)
        if s == '0':
            continue
        lines.append(f'poly g{idx} = {s};')
        members.append(f'g{idx}'); idx += 1
    if not members:
        members = ['0']
    lines.append(f'ideal I = {",".join(members)};')
    if sat_factors:
        parts = [f'({poly_str(f)})' for f in sat_factors if poly_str(f) != '0']
        nz = '*'.join(parts) if parts else '1'
        lines.append(f'poly nz = {nz};')
        lines.append('ideal Is = sat(I,nz)[1];')
        lines.append('ideal G = std(Is);')
    else:
        lines.append('ideal G = std(I);')
    lines += ['int u = (reduce(1,G)==0);', '"@@UNIT";', 'u;']
    if want_dim:
        lines += ['"@@DIM";', 'dim(G);']
    lines += ['"@@MEM";', 'memory(2);', 'quit;']
    return '\n'.join(lines) + '\n'


def sing_poly_intcoeff(expr, gens):
    """Exact integer-coefficient Singular string over Q (clear a global rational
    denominator; number-field vars r/r1/r2 kept as ring variables)."""
    expr = sp.cancel(sp.sympify(expr))
    if expr == 0:
        return '0'
    num, den = sp.fraction(expr)
    poly = sp.Poly(sp.expand(num), *gens)
    denp = sp.Poly(sp.expand(den), *gens)
    if denp.total_degree() != 0:
        raise ValueError(f'non-constant denominator {den}')
    terms = []
    for monom, coeff in poly.terms():
        c = int(coeff)
        if c == 0:
            continue
        factors = [str(c)] if (c != 1 or all(e == 0 for e in monom)) else []
        for g, e in zip(gens, monom):
            if e == 1:
                factors.append(g.name)
            elif e > 1:
                factors.append(f'{g.name}^{e}')
        terms.append('*'.join(factors) if factors else str(c))
    return ('+'.join(terms)).replace('+-', '-') or '0'


# ==========================================================================
#  CASE BUILDERS -> dict(gens, ring_vars, sat_factors, root_syms, minpoly, note)
# ==========================================================================
def case_R9(z, ncoef=8):
    """Build only the requested z (build_system1 builds all z=0..6, too slow)."""
    import convolution_elim_qsupport as qs
    import convolution_elim as ce
    state = qs.build_qsupport_ansatz(z)
    engine = ce.HighCoefficientEngine(state.ansatz, start_degree=qs.START_DEGREE,
                                      target_count=ncoef, c=ce.DEFAULT_C)
    degrees = list(range(qs.START_DEGREE, qs.START_DEGREE - ncoef, -1))
    coeffs = [qs.quotient_reduce(engine.master_coefficient(d)) for d in degrees]
    coeffs = [c for c in coeffs if c != 0]
    Gr = qs.quotient_reduce(state.G.subs(qs.y, qs.r))
    sat_factors = [qs.gamma, state.g_coefficients[-1], Gr]
    rv = mt.ring_vars_of(coeffs + sat_factors)
    return {'gens': coeffs, 'ring_vars': rv, 'sat_factors': sat_factors,
            'root_syms': [qs.r],
            'minpoly': [sum(c * qs.r**(4 - i) for i, c in enumerate(Q_COEFFS))],
            'note': f'z={z}, {len(coeffs)} master coeffs, sat(gamma,g_last,G(r)), r->root'}


def case_sys3(name_sub):
    for ss in mt.build_system3():
        if name_sub in ss['name']:
            return {'gens': ss['gens'], 'ring_vars': ss['ring_vars'],
                    'sat_factors': None, 'root_syms': ss.get('root_syms', []),
                    'note': ss['name'] + ' | ' + ss['note']}
    raise KeyError(name_sub)


def case_deg6_alt(bid, deg_d2=6, depth=12):
    """deg-d2 alt tie tower (free d2) built from d2_threshold, as generator list."""
    import d2_threshold as dt
    d2_poly, sig_poly, e_poly, root_vars, Dc = dt.build_state(bid, deg_d2)
    red = dt.reducer(root_vars)
    C = [red(c) for c in dt.h0_top(d2_poly, deg_d2, sig_poly, e_poly, depth)]
    gens = [c for c in C if c != 0]
    d2_ring = Dc[max(0, len(Dc) - depth):] if Dc is not None else []
    sat_scalars = [dt.S, dt.E] + ([Dc[-1]] if Dc is not None else [])
    root_syms = list(root_vars)
    return {'gens': gens, 'sat_factors': sat_scalars,
            'root_syms': root_syms,
            'ring_vars': mt.ring_vars_of(gens, extra=[dt.S, dt.E] + list(d2_ring)
                                         + list(root_vars)),
            'minpoly': [dt.qpoly(v) for v in root_vars],
            'distinct': (dt.r1, dt.r2) if root_vars == [dt.r1, dt.r2] else None,
            'note': f'{bid} deg_d2={deg_d2} depth={depth} free-d2 tie tower, '
                    f'{len(gens)} tie coeffs'}


def case_sub2(state_key, max_coeffs=6):
    """Reconstruct one sub2 PENDING state and emit its top master coefficients."""
    import phase_f2_sub2 as f2
    import convolution_descent as cd
    tgts = f2.load_targets(f2.TARGET_CELLS, max_defect=1)
    want_idx = int(state_key.split('#state')[1])
    want_cell = state_key.split('#')[0]
    match = None
    for (cellid, case, st, idx, mx, pdelta) in tgts:
        if cellid == want_cell and idx == want_idx:
            match = (case, st, pdelta); break
    if match is None:
        raise KeyError(state_key)
    case, st, pdelta = match
    combo, why, nsol = f2.unique_split(case, st, pdelta)
    if combo is None:
        raise RuntimeError(f'no unique split: {why}')
    polys, scalars, marked, d2_mode, cofactors = f2.reconstruct(case, st, combo, pdelta)
    params = ((f2.r,) if marked is not None else ()) + tuple(cofactors)
    ans = cd.build_ansatz(d2=polys['d2'], d1=polys['d1'], e=polys['e'],
                          sigma=polys['sigma'], parameters=params)
    eng = cd.ConvolutionDescent(ans, c=f2.C_VAL)
    top = f2.engine_top(eng)
    gens = []
    for n in range(max_coeffs):
        mc = f2.redq(eng.master_coefficient(top - n), marked)
        if mc != 0:
            gens.append(mc)
    order_scalars = list(scalars) + list(cofactors)
    root_syms = [f2.r] if marked is not None else []
    return {'gens': gens, 'sat_factors': list(scalars),
            'root_syms': root_syms,
            'ring_vars': mt.ring_vars_of(gens, extra=order_scalars + root_syms),
            'minpoly': [f2.QR_EXPR] if marked is not None else [],
            'note': f'{state_key} degs={[st["deg_d1"],st["deg_sigma"],st["deg_d2"],st["deg_e"]]}'
                    f' d2_mode={d2_mode} {len(gens)} master coeffs, field='
                    + ('Q[r]/(q)' if marked is not None else 'Q')}


CASES = {
    'R9_z4': lambda: case_R9(4), 'R9_z5': lambda: case_R9(5), 'R9_z6': lambda: case_R9(6),
    'a11_b1111_T1_17': lambda: case_sys3('a11_b1111_T1_sup17'),
    'a11_b3100_T2_d6': lambda: case_deg6_alt('a11_b3100_T2', 6),
    'a12_b1110_T2_d6': lambda: case_deg6_alt('a12_b1110_T2', 6),
    'sub2_s14': lambda: case_sub2('sub2:a9_b1000_T1_sz0_dz0_gz-#state14'),
    'sub2_s38': lambda: case_sub2('sub2:a9_b1000_T1_sz0_dz0_gz-#state38'),
    'sub2_s94': lambda: case_sub2('sub2:a9_b1000_T1_sz0_dz0_gz-#state94'),
    'sub2_s263': lambda: case_sub2('sub2:a9_b1000_T1_sz0_dz0_gz-#state263'),
    'sub2_s268': lambda: case_sub2('sub2:a9_b1000_T1_sz0_dz0_gz-#state268'),
}


def run_modp(case_name, timeout=600.0, p=BIGP):
    spec = CASES[case_name]()
    root_syms = spec.get('root_syms', [])
    roots = q_root_mod_p(p, len(root_syms))
    subst = {s: v for s, v in zip(root_syms, roots)}
    gens = [sp.sympify(g).subs(subst) for g in spec['gens']]
    satf = ([sp.sympify(f).subs(subst) for f in spec['sat_factors']]
            if spec.get('sat_factors') else None)
    rvars = mt.ring_vars_of(gens + (satf or []))
    prog = emit_program(gens, rvars, p, sat_factors=satf)
    res = run_singular_measured(prog, timeout)
    res.update({'case': case_name, 'char': p, 'nvars': len(rvars),
                'ngens': len(gens), 'note': spec['note']})
    return res


def run_ratcurve(case_name, caps_gb, timeout=600.0):
    """Over-Q Singular under increasing ulimit -v caps -> honest N-GB estimate."""
    spec = CASES[case_name]()
    gens = list(spec['gens'])
    minp = list(spec.get('minpoly', []))
    satf = spec.get('sat_factors')
    distinct = spec.get('distinct')
    rvars = spec['ring_vars']
    # add a distinctness saturation var for the two-root field
    sat_extra = []
    if distinct is not None:
        wd = sp.Symbol('wd'); r1, r2 = distinct
        gens = gens + [wd * (r1 - r2) - 1]
        rvars = mt.ring_vars_of(gens + minp + (satf or []), extra=[wd])
    out = []
    for cap in caps_gb:
        prog = emit_program(gens, rvars, 0, sat_factors=satf, minpoly_gens=minp)
        r = run_singular_measured(prog, timeout, vcap_gb=cap)
        r['cap_gb'] = cap
        out.append(r)
        print(f'  cap={cap}GB -> {r["verdict"]} wall={r["wall"]}s '
              f'rss={r.get("rss_kb")}kb', flush=True)
        if r['verdict'] in ('UNIT', 'PROPER'):
            break  # finished within this cap; N found
    return {'case': case_name, 'note': spec['note'], 'curve': out}


if __name__ == '__main__':
    if len(sys.argv) < 2 or sys.argv[1] == 'list':
        print('cases:', ' '.join(CASES)); sys.exit(0)
    cmd = sys.argv[1]
    if cmd == 'ratcurve':
        bidcase = sys.argv[2]
        caps = [float(x) for x in sys.argv[3:]] or [2.0, 4.0, 8.0]
        res = run_ratcurve(bidcase, caps)
    else:
        to = float(os.environ.get('DIAG_TIMEOUT', '600'))
        res = run_modp(cmd, timeout=to)
        print(json.dumps(res, indent=1))
    # append to results file
    path = 'blowup_diagnosis_results.json'
    allres = json.load(open(path)) if os.path.exists(path) else []
    allres.append({'ts': time.strftime('%H:%M:%S'), 'cmd': ' '.join(sys.argv[1:]),
                   'result': res})
    json.dump(allres, open(path, 'w'), indent=1)
