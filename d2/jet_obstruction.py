#!/usr/bin/env python3
"""jet_obstruction.py -- JET-LIFTING OBSTRUCTION DEPTH profiler.

Replaces "the Groebner basis timed out" with a geometric diagnosis.

SETUP.  A cascade state supplies, per depth n, one accumulated master
coefficient f_n (top degree, top-1, top-2, ...).  Following the committed
J6 construction (`j6_msolve.py`), the depth-N truncation is the affine scheme

    X_N = V( f_1, ..., f_N, class relations, w*prod(scalars) - 1 )

and pi_N : X_{N+1} -> X_N is the forgetful map.  X_{N+1} = X_N cut by the
one new equation f_{N+1}.

OBSTRUCTION.  At a point x of X_N the derivative of the NEW equations with
respect to the NEW unknowns is the relative tangent map; its cokernel is the
first-order obstruction to lifting x.  This profiler computes, per depth:

  * dim X_N and (when dim = 0) deg X_N = dim_k k[vars]/I_N, plus a
    linear-algebra lower bound on the number of DISTINCT points;
  * the set of NEW unknowns introduced at the step (see [D1] below);
  * the obstruction operator.  When no new unknowns are introduced, pi_N is
    a closed immersion, the relative tangent map has zero columns, and the
    first-order obstruction degenerates to the full multiplication operator

        M_{f_{N+1}} :  A_N --> A_N ,    A_N = k[vars]/I_N ,

    whose cokernel is exactly A_{N+1}.  rank M and coker M are then plain
    linear algebra over F_p on a deg(X_N) x deg(X_N) matrix -- NO second
    Groebner basis is needed to decide whether depth N+1 is empty.
  * lifts exist at the step  <=>  coker != 0  <=>  X_{N+1} nonempty.

REGIMES.  A = bounded obstruction (some X_{N0} empty; beyond N0 nothing
lifts).  B = eventual formal smoothness (pi_N generically smooth+surjective,
stable fibre dimension; the tail is formally UNOBSTRUCTED).  C = periodic /
templated (ranks and components repeat).

[D1] DEGENERACY, stated loudly.  In the J6 construction the class relations
and the Rabinowitsch saturation are present at EVERY depth, so the ambient
variable set is constant in N and the set of new unknowns is EMPTY at every
step.  The relative Jacobian is therefore an (#new eqs) x 0 matrix: rank 0,
cokernel = the whole new-equation space at every point, i.e. "nothing lifts
unless f_{N+1} happens to vanish there".  That is not a bug in the profiler
and it is not a vacuous statement -- it is the reason the multiplication
operator M_{f_{N+1}} is the correct obstruction object here.  The profiler
still reports the FIRST-APPEARANCE filtration of the unknowns (which class
coefficient first shows up in which f_n) as data, because that is the datum
a construction with depth-staged unknowns would turn into new columns.

PROVENANCE.  Generators, class relations and saturation are read VERBATIM as
strings from `j6_msolve_results.json` (the committed record of the J6 pass).
With --replay the whole chain is re-derived through the original committed
code path (alt_hunt_depth2.reconstruct_general + ConvolutionDescent) and
asserted string-identical; this costs ~16 s/state.

New file; READ-ONLY on all existing artifacts.  Usage:
    python jet_obstruction.py                 # full profile + JSON + table
    python jet_obstruction.py --replay        # + re-derive gens from source
    python jet_obstruction.py --quiet         # self-check, exit 0 iff the
                                              #   J6 control reproduces
Env: JO_PRIMES (comma list), JO_SING_CAP (per-call WSL wall s, default 300).
Output: jet_obstruction_results.json
"""
import json
import os
import re
import subprocess
import sys
import time

import numpy as np
import sympy as sp

import blowup_diagnosis as bd

HERE = os.path.dirname(os.path.abspath(__file__))
QUIET = '--quiet' in sys.argv
REPLAY = '--replay' in sys.argv
MSOLVE = '--msolve' in sys.argv
SING_CAP = int(os.environ.get('JO_SING_CAP', '300'))
PRIMES = [int(x) for x in os.environ.get(
    'JO_PRIMES', '10007,10009,32003,100019').split(',')]
MSOLVE_CAP = int(os.environ.get('JO_MSOLVE_CAP', '600'))
MSOLVE_DEPTHS = [int(x) for x in
                 os.environ.get('JO_MSOLVE_DEPTHS', '2,3').split(',')]
J6_JSON = os.path.join(HERE, 'j6_msolve_results.json')
RESULTS = os.path.join(HERE, 'jet_obstruction_results.json')
MAXDEPTH = 3            # J6 records 3 accumulated master coefficients


def log(*a):
    if not QUIET:
        print(*a, flush=True)


# --------------------------------------------------------------------------
# Singular driver (WSL; /mnt/c is broken there, so the script is piped to
# $HOME via stdin and run from $HOME)
# --------------------------------------------------------------------------
def singular(prog, cap=None, tag='jo'):
    cap = cap or SING_CAP
    cmd = (f'cat > $HOME/_{tag}.sing; cd $HOME && '
           f'timeout {cap} Singular -q $HOME/_{tag}.sing')
    t0 = time.monotonic()
    cp = subprocess.run(('wsl.exe', '-d', 'Ubuntu', '--', 'bash', '-lc', cmd),
                        input=prog, stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE, text=True, encoding='utf-8',
                        timeout=cap + 90, check=False)
    kv, raw = {}, []
    for ln in cp.stdout.splitlines():
        raw.append(ln)
        if ln.startswith('@') and '|' in ln:
            k, v = ln.split('|', 1)
            kv[k] = v
    return {'kv': kv, 'raw': raw, 'rc': cp.returncode,
            'stderr': cp.stderr.strip()[:2000],
            'wall': round(time.monotonic() - t0, 2)}


# --------------------------------------------------------------------------
# polynomial string <-> exponent-vector parsing (F_p)
# --------------------------------------------------------------------------
_TERM = re.compile(r'([+-]?)([^+-]+)')


def parse_poly(s, idx, p):
    """Singular expanded-poly string -> {exponent tuple: coeff mod p}."""
    s = s.replace(' ', '')
    if s in ('0', ''):
        return {}
    d = {}
    for sign, body in _TERM.findall(s):
        c, mon = 1, [0] * len(idx)
        for fac in body.split('*'):
            if not fac:
                continue
            if '^' in fac:
                b, e = fac.split('^')
                e = int(e)
            else:
                b, e = fac, 1
            if b.isdigit():
                c = c * pow(int(b), e, p) % p
            elif b in idx:
                mon[idx[b]] += e
            else:
                raise ValueError(f'unparsable factor {b!r} in {s[:80]!r}')
        t = tuple(mon)
        d[t] = (d.get(t, 0) + (-c % p if sign == '-' else c)) % p
    return {k: v for k, v in d.items() if v}


# --------------------------------------------------------------------------
# dense linear algebra mod p (numpy int64; p < 2^31 so products fit)
# --------------------------------------------------------------------------
def rref_rank_det(A, p):
    """Gaussian elimination mod p -> (rank, det or None if not square/singular)."""
    A = np.array(A, dtype=np.int64) % p
    rows, cols = A.shape
    r, det, swaps = 0, 1, 0
    for c in range(cols):
        nz = np.nonzero(A[r:, c])[0]
        if nz.size == 0:
            det = 0
            continue
        piv = r + int(nz[0])
        if piv != r:
            A[[r, piv]] = A[[piv, r]]
            swaps += 1
        pv = int(A[r, c])
        det = det * pv % p
        inv = pow(pv, p - 2, p)
        A[r] = A[r] * inv % p
        col = A[r + 1:, c].copy()
        nzr = np.nonzero(col)[0]
        if nzr.size:
            A[r + 1 + nzr] = (A[r + 1 + nzr] - np.outer(col[nzr], A[r])) % p
        r += 1
        if r == rows:
            break
    if rows != cols or r < rows:
        det = 0 if r < min(rows, cols) else None
    else:
        det = det * (-1 if swaps % 2 else 1) % p
    return r, det


def cyclic_min_poly(M, v, p):
    """Coefficients (low degree first) of the monic minimal polynomial of M on
    the cyclic module generated by v.  Incremental elimination on the Krylov
    sequence; `coords[i]` records basis[i] as a polynomial in M applied to v,
    kept as a numpy array so the bookkeeping is a matvec, not a python loop."""
    n = M.shape[0]
    basis = np.zeros((n + 1, n), dtype=np.int64)
    coords = np.zeros((n + 1, n + 1), dtype=np.int64)
    pivots, m = [], 0                       # m = number of basis vectors
    cur = np.array(v, dtype=np.int64) % p
    for k in range(n + 1):
        red = cur.copy()
        rep = np.zeros(m, dtype=np.int64)
        for i in range(m):                  # forward substitution (sequential)
            c = int(red[pivots[i]])
            if c:
                rep[i] = c
                red = (red - c * basis[i]) % p
        prev = (-(rep @ coords[:m])) % p if m else np.zeros(n + 1, dtype=np.int64)
        prev[k] = (prev[k] + 1) % p
        nz = np.nonzero(red)[0]
        if nz.size == 0:                    # dependency -> minimal polynomial
            co = [int(x) for x in prev[:k + 1]]
            lead = co[k]
            inv = pow(lead, p - 2, p)
            return [x * inv % p for x in co]
        pc = int(nz[0])
        inv = pow(int(red[pc]), p - 2, p)
        basis[m] = red * inv % p
        coords[m] = prev * inv % p
        pivots.append(pc)
        m += 1
        cur = M.dot(cur) % p
    raise RuntimeError('cyclic_min_poly: no dependency within n+1 steps')


def squarefree_deg(coeffs, p):
    """deg of the squarefree part of a univariate poly given low-degree-first
    coefficients over F_p -- i.e. its number of distinct roots in the closure."""
    t = sp.Symbol('t')
    f = sp.Poly(list(reversed(coeffs)), t, modulus=p)
    g = sp.gcd(f, f.diff(t))
    return f.degree() - (g.degree() if g.degree() > 0 else 0)


# --------------------------------------------------------------------------
# system assembly (verbatim strings out of the J6 record)
# --------------------------------------------------------------------------
def load_state(rec):
    V = [sp.Symbol(n) for n in rec['ring_vars']]
    loc = {v.name: v for v in V}
    ex = [sp.sympify(g['coefficient'], locals=loc) for g in rec['gens']]
    rels = [sp.sympify(x, locals=loc) for x in rec['class_relations']]
    sat = sp.sympify(rec['saturation'], locals=loc)
    return {'key': rec['key'], 'vars': V, 'gen_deg': [g['degree'] for g in rec['gens']],
            'gens': ex, 'rels': rels, 'sat': sat,
            'kill_depth': rec.get('kill_depth'),
            'msolve': [(a['depth'], a['msolve'], a.get('wall'),
                        a.get('out_head', '')[:70]) for a in rec['msolve_attempts']]}


def prime_safety(st):
    """Primes that MUST be avoided, and the first-appearance filtration.

    A prime is unsafe if it divides a coefficient denominator (the reduction
    of that generator is then undefined) or if it divides the integer content
    of the denominator-cleared generator (the generator would reduce to 0 and
    the mod-p ideal would be strictly smaller than the true reduction)."""
    bad, detail = set(), []
    for e in list(st['gens']) + list(st['rels']) + [st['sat']]:
        pol = sp.Poly(sp.expand(e), *st['vars'])
        den = sp.Integer(1)
        for c in pol.coeffs():
            den = sp.ilcm(den, sp.Rational(c).q)
        cont = sp.Integer(0)
        for c in pol.coeffs():
            cont = sp.gcd(cont, sp.Integer(sp.Rational(c) * den))
        dp_ = sorted(int(q) for q in sp.factorint(den))
        cp_ = sorted(int(q) for q in sp.factorint(cont)) if cont != 1 else []
        bad |= set(dp_) | set(cp_)
        detail.append({'den_primes': dp_, 'content_primes': cp_})
    st['prime_detail'] = detail
    # first-appearance depth of each variable across f_1, f_2, ...
    first, seen = {}, set()
    for n, e in enumerate(st['gens'], start=1):
        for s in sorted(e.free_symbols, key=str):
            if s.name not in seen:
                seen.add(s.name)
                first[s.name] = n
    for v in st['vars']:
        first.setdefault(v.name, None)
    return sorted(bad), first


def new_unknowns(st, n):
    """Unknowns introduced by the step X_n -> X_{n+1} under the J6 ideal
    construction (class relations + saturation present at every depth)."""
    old = set()
    for e in st['gens'][:n] + st['rels'] + [st['sat']]:
        old |= {s.name for s in e.free_symbols}
    return sorted({s.name for s in st['gens'][n].free_symbols} - old)


# --------------------------------------------------------------------------
# EXACT denominator-clearing Singular serializer.
#
# [B1] We do NOT use blowup_diagnosis.sing_poly_intcoeff here: under this
# environment (sympy 1.14.0 + python-flint 0.9.0) it silently truncates
# rational coefficients to integers, because its denominator-clearing step
#     num, den = sp.fraction(sp.cancel(expr))
# is a no-op -- sp.cancel(x/2 + 1/3) returns `x/2 + 1/3`, not `(3x+2)/6`, so
# `den` is 1, `sp.Poly(num)` keeps Rational coefficients, and the subsequent
# `int(coeff)` truncates toward zero (195/2048 -> 0, 3981312/221 -> 18014).
# `bugcheck()` below demonstrates this on the live J6 generators on every run.
# --------------------------------------------------------------------------
def sing(e, V):
    """Exact integer-coefficient Singular polynomial string.

    Clears the lcm of the coefficient denominators (a unit of Q, so the
    ideal is unchanged) and VERIFIES the emitted string by reparsing it."""
    pol = sp.Poly(sp.expand(e), *V)
    if pol.is_zero:
        return '0'
    den = sp.Integer(1)
    for c in pol.coeffs():
        den = sp.ilcm(den, sp.Rational(c).q)
    terms = []
    for monom, coeff in pol.terms():
        c = sp.Rational(coeff) * den
        assert c.q == 1, (coeff, den)
        c = int(c)
        if c == 0:
            continue
        facs = [] if (c == 1 and any(monom)) else [str(c)]
        for g, k in zip(V, monom):
            if k == 1:
                facs.append(g.name)
            elif k > 1:
                facs.append(f'{g.name}^{k}')
        terms.append('*'.join(facs))
    s = ('+'.join(terms)).replace('+-', '-') or '0'
    back = sp.sympify(s.replace('^', '**'), locals={v.name: v for v in V})
    assert sp.expand(back - sp.expand(e) * den) == 0, f'serializer mismatch: {s[:80]}'
    return s


def bugcheck(states):
    """Quote the discrepancy between the exact serializer and the committed
    blowup_diagnosis.sing_poly_intcoeff on the live J6 systems."""
    rows = []
    for st in states:
        V = st['vars']
        for label, e in ([(f'class_relation[{i}]', x)
                          for i, x in enumerate(st['rels'])]
                         + [(f'master_coefficient[deg {d}]', g)
                            for d, g in zip(st['gen_deg'], st['gens'])]
                         + [('saturation', st['sat'])]):
            good = sing(e, V)
            try:
                bad = bd.sing_poly_intcoeff(e, V)
            except Exception as exc:                        # noqa: BLE001
                bad = f'<raised {type(exc).__name__}: {exc}>'
            if good != bad:
                rows.append({'state': st['key'], 'poly': label,
                             'exact_terms': good.count('+') + good.count('-') + 1,
                             'bd_terms': bad.count('+') + bad.count('-') + 1,
                             'exact_head': good[:110], 'bd_head': bad[:110]})
    return rows


# --------------------------------------------------------------------------
# the profiler
# --------------------------------------------------------------------------
def profile(st, p, points=True):
    V, vn = st['vars'], [v.name for v in st['vars']]
    idx = {n: i for i, n in enumerate(vn)}
    G = [sing(e, V) for e in st['gens']]
    base = [sing(e, V) for e in st['rels']] + [sing(st['sat'], V)]
    L = [f'ring R = {p},({",".join(vn)}),dp;']
    for N in range(1, MAXDEPTH + 1):
        L += [f'ideal I{N} = {",".join(base + G[:N])};',
              f'ideal Gb{N} = std(I{N});',
              f'"@DIM{N}|"+string(dim(Gb{N}));',
              f'"@VDIM{N}|"+string(vdim(Gb{N}));',
              f'"@SIZE{N}|"+string(size(Gb{N}));',
              f'"@GB{N}|"+string(Gb{N});']
    # obstruction operator at every step whose source algebra is 0-dimensional
    for N in range(1, MAXDEPTH):
        L += [f'if (dim(Gb{N})==0 && vdim(Gb{N})>0) {{',
              f'  ideal B{N} = kbase(Gb{N}); int n{N} = size(B{N});',
              f'  "@KN{N}|"+string(n{N}); int j{N};',
              f'  for(j{N}=1;j{N}<=n{N};j{N}++)'
              f'{{ "@KB{N}_"+string(j{N})+"|"+string(B{N}[j{N}]); }}',
              f'  poly f{N} = reduce({G[N]}, Gb{N});',
              f'  for(j{N}=1;j{N}<=n{N};j{N}++)'
              f'{{ "@MF{N}_"+string(j{N})+"|"+string(reduce(f{N}*B{N}[j{N}], Gb{N})); }}']
        # multiplication matrix of every variable -> any linear form for free
        if points:
            for vi, v in enumerate(vn):
                L += [f'  for(j{N}=1;j{N}<=n{N};j{N}++)'
                      f'{{ "@MV{N}_{vi}_"+string(j{N})+"|"'
                      f'+string(reduce({v}*B{N}[j{N}], Gb{N})); }}']
        L += ['}']
    # independent cross-check on the point structure: lex GB by FGLM, factor
    # the univariate eliminant in the last variable
    lexvars = [n for n in vn if n != vn[0]] + [vn[0]]
    L += [f'ring RL = {p},({",".join(lexvars)}),lp;',
          'if (dim(std(fetch(R,I2)))==0) {',
          '  ideal GL = stdfglm(fetch(R,I2)); int k; poly u = 0;',
          # the univariate eliminant is the generator in ONE variable; in a
          # lex GB that is the one with the smallest leading monomial, which
          # Singular lists FIRST -- locate it by variable count, not position
          '  for(k=1;k<=size(GL);k++)'
          '{ if(size(variables(GL[k]))==1){ u = GL[k]; break; } }',
          '  if (u != 0) {',
          '    "@ELIMVAR|"+string(variables(u));',
          '    "@ELIMDEG|"+string(deg(u));',
          '    list fa = factorize(u); string s = "";',
          '    for(k=1;k<=size(fa[1]);k++)'
          '{ s = s + string(deg(fa[1][k])) + "^" + string(fa[2][k]) + " "; }',
          '    "@ELIMFACT|"+s;',
          '  }',
          '}', 'setring R;']
    L.append('quit;')
    r = singular('\n'.join(L) + '\n', tag='jo_' + re.sub(r'\W', '_', st['key'])[:40])
    kv = r['kv']
    if r['rc'] != 0 or '@DIM1' not in kv:
        return {'prime': p, 'error': 'singular failed', 'rc': r['rc'],
                'stderr': r['stderr'], 'wall': r['wall']}

    depths = []
    for N in range(1, MAXDEPTH + 1):
        d, vd = int(kv[f'@DIM{N}']), int(kv[f'@VDIM{N}'])
        gb = kv[f'@GB{N}']
        depths.append({
            'depth': N, 'master_degree': st['gen_deg'][N - 1],
            'n_equations': N + len(st['rels']) + 1,
            'dim': d, 'vdim': vd, 'gb_size': int(kv[f'@SIZE{N}']),
            'empty': (d == -1), 'gb_is_unit_ideal': (gb.strip() == '1'),
            'raw_dim': f'@DIM{N}|{kv[f"@DIM{N}"]}',
            'raw_vdim': f'@VDIM{N}|{kv[f"@VDIM{N}"]}',
            'raw_gb_head': f'@GB{N}|{gb[:120]}'})

    steps = []
    for N in range(1, MAXDEPTH):
        src, dst = depths[N - 1], depths[N]
        nu = new_unknowns(st, N)
        step = {'step': f'{N}->{N + 1}', 'new_unknowns': nu,
                'n_new_unknowns': len(nu), 'n_new_equations': 1,
                'rel_jac_shape': [1, len(nu)], 'rel_jac_rank': 0 if not nu else None,
                'dim_src': src['dim'], 'dim_dst': dst['dim']}
        if nu:
            step['note'] = ('new unknowns present -- relative tangent map has '
                            'columns; NOT the J6 shape, see [D1]')
        else:
            step['note'] = ('no new unknowns: pi is a closed immersion, relative '
                            'tangent map is 1x0 (rank 0, cokernel = the whole '
                            'new-equation line at every point) -- obstruction '
                            'carried by the multiplication operator [D1]')
        if f'@KN{N}' in kv:
            n = int(kv[f'@KN{N}'])
            B = [parse_poly(kv[f'@KB{N}_{j}'], idx, p) for j in range(1, n + 1)]
            pos = {}
            for j, b in enumerate(B):
                assert len(b) == 1 and list(b.values()) == [1], (N, j, b)
                pos[next(iter(b))] = j
            def mat(pref):
                M = np.zeros((n, n), dtype=np.int64)
                for j in range(1, n + 1):
                    for m, c in parse_poly(kv[f'{pref}{N}_{j}'], idx, p).items():
                        M[pos[m], j - 1] = c
                return M
            Mf = mat('@MF')
            rk, det = rref_rank_det(Mf, p)
            step.update({
                'obstruction_operator': f'M_f{N + 1} on A_{N}',
                'source_dim_k': n, 'op_rank': rk, 'op_cokernel_dim': n - rk,
                'op_det_mod_p': int(det) if det is not None else None,
                'op_is_unit': (rk == n),
                'lifts_exist': (n - rk) > 0,
                'coker_matches_vdim_next': (n - rk) == dst['vdim']})
            # ---- point structure of X_N by pure linear algebra --------------
            # M_v for each variable v; any linear form's operator is then a
            # numpy linear combination.  deg minpoly(l) and its squarefree
            # part bound the number of DISTINCT points from below (with
            # equality iff l separates them); a non-squarefree minpoly PROVES
            # A_N is non-reduced (a reduced finite algebra has squarefree
            # minimal polynomials for every element).
            try:
                if not points or f'@MV{N}_0_1' not in kv:
                    raise RuntimeError('point structure not requested for this '
                                       'prime (computed on the first prime only)')
                Mvs = []
                for vi in range(len(vn)):
                    M = np.zeros((n, n), dtype=np.int64)
                    for j in range(1, n + 1):
                        for m, c in parse_poly(kv[f'@MV{N}_{vi}_{j}'], idx, p).items():
                            M[pos[m], j - 1] = c
                    Mvs.append(M)
                one = pos.get(tuple([0] * len(vn)))
                if one is None:
                    raise RuntimeError('1 is not in the standard monomial basis')
                e1 = np.zeros(n, dtype=np.int64)
                e1[one] = 1
                rng = np.random.default_rng(20260724 + p)
                trials, best, nonred = [], 0, False
                for _ in range(2):
                    co = rng.integers(1, p, size=len(vn))
                    Ml = np.zeros((n, n), dtype=np.int64)
                    for c, M in zip(co, Mvs):
                        Ml = (Ml + int(c) * M) % p
                    mp = cyclic_min_poly(Ml, e1, p)
                    dmp, dsf = len(mp) - 1, squarefree_deg(mp, p)
                    trials.append({'coeffs': [int(x) for x in co],
                                   'minpoly_deg': dmp, 'squarefree_deg': dsf})
                    best = max(best, dsf)
                    nonred |= (dmp != dsf)
                step['point_structure'] = {
                    'deg_with_multiplicity': n,
                    'distinct_points_lower_bound': best,
                    'non_reduced_proved': bool(nonred),
                    'random_linear_form_trials': trials}
            except Exception as exc:                       # noqa: BLE001
                step['point_structure_error'] = f'{type(exc).__name__}: {exc}'
        else:
            step.update({'obstruction_operator': None,
                         'lifts_exist': not dst['empty'],
                         'reason': f'A_{N} is not finite dimensional '
                                   f'(dim X_{N} = {src["dim"]}); the finite '
                                   f'multiplication matrix does not exist. '
                                   f'Reporting the dimension drop instead.',
                         'dim_drop': (src['dim'] - dst['dim'])
                                     if dst['dim'] >= 0 else None})
        steps.append(step)
    elim = None
    if '@ELIMDEG' in kv:
        elim = {'variable': kv.get('@ELIMVAR'), 'degree': int(kv['@ELIMDEG']),
                'factor_degrees_with_multiplicity': kv.get('@ELIMFACT', '').strip(),
                'engine': 'Singular stdfglm (FGLM to lex) + factorize'}
    return {'prime': p, 'wall': r['wall'], 'depths': depths, 'steps': steps,
            'x2_lex_eliminant': elim,
            'raw_singular_markers': [ln[:160] for ln in r['raw']
                                     if ln.startswith(('@DIM', '@VDIM', '@SIZE',
                                                       '@GB', '@KN', '@ELIM'))]}


def msolve_depth(st, n, budget=None):
    """Re-run msolve over Q on the depth-n system built with the EXACT
    serializer (j6_msolve.py built the same system with the truncating one).
    Returns the verdict plus msolve's raw stdout head."""
    V = st['vars']
    polys = [s for s in (sing(e, V) for e in
                         st['gens'][:n] + st['rels'] + [st['sat']]) if s != '0']
    prog = ','.join(v.name for v in V) + '\n0\n' + ',\n'.join(polys) + '\n'
    budget = budget or MSOLVE_CAP
    tag = re.sub(r'\W', '_', st['key'])[:40] + f'_d{n}'
    fn, on = f'jo_{tag}.ms', f'jo_{tag}.out'
    t0 = time.monotonic()
    try:
        cp = subprocess.run(
            ('wsl.exe', '-d', 'Ubuntu', '--', 'bash', '-lc',
             f'cat > $HOME/{fn}; cd $HOME && ulimit -v 8388608; '
             f'timeout {budget} $HOME/msolve/msolve -f $HOME/{fn} -o $HOME/{on}; '
             f'echo "@@EXIT:$?"; cat $HOME/{on} 2>/dev/null'),
            input=prog, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, encoding='utf-8', timeout=budget + 120, check=False)
    except subprocess.TimeoutExpired:
        # orphan-proofing: the WSL-side `timeout` still owns the process
        return {'depth': n, 'verdict': 'COST', 'exit': None,
                'wall': round(time.monotonic() - t0, 1), 'ms_file': fn,
                'n_polys': len(polys),
                'raw_out_head': '<relay timeout; WSL-side timeout pending>'}
    out = cp.stdout or ''
    m = re.search(r'@@EXIT:(\d+)', out)
    exitc = int(m.group(1)) if m else None
    body = out.split('@@EXIT:', 1)[-1]
    body = body.split('\n', 1)[1].strip() if '\n' in body else ''
    if exitc == 124:
        verdict = 'COST'
    elif body.startswith('[-1]'):
        verdict = 'EMPTY'
    elif body:
        verdict = 'NOT_EMPTY'
    else:
        verdict = 'ERROR'
    return {'depth': n, 'verdict': verdict, 'exit': exitc,
            'wall': round(time.monotonic() - t0, 1), 'ms_file': fn,
            'n_polys': len(polys), 'raw_out_head': body[:200]}


def regime(runs):
    """A / B / C from the per-prime profiles (all primes must agree)."""
    for run in runs:
        if 'error' in run:
            return 'UNDETERMINED', 'a prime run errored'
    empt = [min([d['depth'] for d in r['depths'] if d['empty']], default=None)
            for r in runs]
    if len(set(empt)) != 1:
        return 'UNDETERMINED', f'primes disagree on first empty depth: {empt}'
    N0 = empt[0]
    if N0 is not None:
        return 'A', (f'X_{N0} is empty over every prime tried; beyond depth '
                     f'{N0 - 1} nothing lifts (bounded obstruction)')
    dims = [tuple(d['dim'] for d in r['depths']) for r in runs]
    if len(set(dims)) == 1 and len(set(dims[0])) == 1:
        return 'B?', (f'no depth is empty and dim X_N is constant = {dims[0][0]}; '
                      f'consistent with eventual formal smoothness, NOT '
                      f'established (only {MAXDEPTH} depths seen)')
    return 'UNDETERMINED', f'no depth empty, dims {dims[0]} -- more depths needed'


# --------------------------------------------------------------------------
# provenance replay through the original committed code path
# --------------------------------------------------------------------------
def replay_check(j6):
    import alt_hunt_depth2 as ah
    res = json.load(open(os.path.join(HERE, 'alt_hunt_results.json'),
                         encoding='utf-8'))
    open_states = [s for s in res['states'] if s['verdict'] == 'OPEN']
    targets = {cid: (win, case, idxs)
               for win, cid, case, idxs in ah.load_hunt_targets()}
    byk = {r['key']: r for r in j6['results']}
    out = []
    for stt in open_states:
        t0 = time.monotonic()
        _, case, _ = targets[stt['cellid']]
        combo = tuple(tuple(c) for c in stt['splits'][0]['combo'])
        recon, why = ah.reconstruct_general(case, combo)
        if recon is None:
            out.append({'key': stt['key'], 'ok': False, 'why': why})
            continue
        polys, scalars, unknowns, relations = recon
        ok_polys = ({k: str(v) for k, v in polys.items()}
                    == stt['splits'][0]['polys'])
        ans = ah.cd.build_ansatz(d2=polys['d2'], d1=polys['d1'], e=polys['e'],
                                 sigma=polys['sigma'],
                                 parameters=tuple(unknowns))
        eng = ah.cd.ConvolutionDescent(ans, c=ah.f2.C_VAL)
        top = ah.f2.engine_top(eng)
        rec, ok_gens = byk[stt['key']], []
        for n in range(MAXDEPTH):
            mc = sp.expand(eng.master_coefficient(top - n))
            want = [g for g in rec['gens'] if g['degree'] == top - n]
            ok_gens.append(bool(want) and str(mc) == want[0]['coefficient'])
        ok_rel = [str(x) for x in relations] == rec['class_relations']
        ok_sat = (str(sp.expand(ah.w * sp.Mul(*scalars) - 1))
                  == rec['saturation'])
        out.append({'key': stt['key'], 'ok': bool(ok_polys and all(ok_gens)
                                                  and ok_rel and ok_sat),
                    'polys_verbatim': ok_polys, 'gens_verbatim': ok_gens,
                    'relations_verbatim': ok_rel, 'saturation_verbatim': ok_sat,
                    'top_degree': top, 'secs': round(time.monotonic() - t0, 1)})
        log(f'  replay {stt["key"]}: '
            f'{"OK" if out[-1]["ok"] else "MISMATCH"} ({out[-1]["secs"]}s)')
    return out


# --------------------------------------------------------------------------
def main():
    j6 = json.load(open(J6_JSON, encoding='utf-8'))
    states = [load_state(r) for r in j6['results']]
    out = {'schema': 1,
           'item': 'jet-lifting obstruction depth profile -- J6 control',
           'source': 'j6_msolve_results.json (generators read verbatim)',
           'primes': PRIMES, 'maxdepth': MAXDEPTH, 'states': []}

    out['serializer_bug'] = {
        'helper': 'blowup_diagnosis.sing_poly_intcoeff',
        'env': f'sympy {sp.__version__}',
        'symptom': 'rational coefficients truncated toward zero',
        'root_cause': 'sp.cancel(x/2+1/3) returns x/2+1/3 (fraction den = 1) '
                      'in this environment, so the denominator-clearing step '
                      'is a no-op and the following int(coeff) truncates',
        'on_critical_path_of': 'j6_msolve.py (builds the msolve .ms input)',
        'discrepancies': bugcheck(states)}
    log(f'== serializer bugcheck: '
        f'{len(out["serializer_bug"]["discrepancies"])} of '
        f'{sum(len(s["rels"]) + len(s["gens"]) + 1 for s in states)} polynomials '
        f'differ from blowup_diagnosis.sing_poly_intcoeff')

    if REPLAY:
        log('== provenance replay through alt_hunt_depth2 + ConvolutionDescent')
        out['replay'] = replay_check(j6)
        if not all(r['ok'] for r in out['replay']):
            print('REPLAY MISMATCH -- aborting', file=sys.stderr)
            json.dump(out, open(RESULTS, 'w'), indent=1)
            return 2
        log('  replay: 4/4 verbatim')

    ok_all = True
    for st in states:
        log(f'== {st["key"]}  vars {[v.name for v in st["vars"]]}')
        bad, first = prime_safety(st)
        used = [p for p in PRIMES if p not in bad]
        # the point structure is descriptive (not part of the control criteria);
        # compute it on the first prime only, and never in --quiet self-check
        runs = [profile(st, p, points=(i == 0 and not QUIET))
                for i, p in enumerate(used)]
        for r in runs:
            if 'error' in r:
                log(f'   p={r["prime"]}: ERROR {r.get("stderr","")[:200]}')
                continue
            d = r['depths']
            log(f'   p={r["prime"]} ({r["wall"]}s): '
                + ' '.join(f'X_{x["depth"]}(dim {x["dim"]},deg {x["vdim"]})'
                           for x in d))
            for s in r['steps']:
                log(f'      {s["step"]}: new_unknowns={s["n_new_unknowns"]} '
                    f'rank={s.get("op_rank")} coker={s.get("op_cokernel_dim")} '
                    f'lifts={s["lifts_exist"]}')
        reg, why = regime(runs)
        # ---- J6 control criteria -----------------------------------------
        crit = {}
        for r in runs:
            if 'error' in r:
                crit[r['prime']] = {'ok': False, 'why': 'singular error'}
                continue
            d = {x['depth']: x for x in r['depths']}
            s2 = [x for x in r['steps'] if x['step'] == '2->3'][0]
            c = {'depth2_satisfiable': d[2]['dim'] == 0 and d[2]['vdim'] > 0,
                 'depth3_empty': d[3]['empty'] and d[3]['gb_is_unit_ideal']
                                 and d[3]['vdim'] == 0,
                 'obstruction_at_step_2to3': s2.get('op_cokernel_dim') == 0
                                             and s2.get('op_is_unit') is True,
                 'linalg_matches_gb': s2.get('coker_matches_vdim_next') is True}
            c['ok'] = all(c.values())
            crit[r['prime']] = c
            ok_all &= c['ok']
        srec = {'key': st['key'], 'vars': [v.name for v in st['vars']],
                'master_degrees': st['gen_deg'],
                'n_class_relations': len(st['rels']),
                'forbidden_primes': bad, 'primes_used': used,
                'first_appearance_depth': first,
                'j6_recorded_kill_depth': st['kill_depth'],
                'j6_recorded_msolve': st['msolve'],
                'regime': reg, 'regime_reason': why,
                'control_criteria': crit, 'runs': runs}
        if MSOLVE:
            srec['msolve_rerun_exact_serializer'] = [
                msolve_depth(st, n) for n in MSOLVE_DEPTHS]
            for a in srec['msolve_rerun_exact_serializer']:
                log(f'   msolve(exact) depth {a["depth"]}: {a["verdict"]} '
                    f'({a["wall"]}s) raw={a["raw_out_head"][:40]!r}')
        out['states'].append(srec)
        log(f'   -> regime {reg}: {why}')
        json.dump(out, open(RESULTS, 'w'), indent=1)

    out['control_reproduces'] = bool(ok_all)
    json.dump(out, open(RESULTS, 'w'), indent=1)
    if not QUIET:
        print('\n' + table(out))
    print(('SELF-CHECK PASS: J6 control reproduces '
           '(depth-2 satisfiable, depth-3 empty, obstruction at step 2->3)')
          if ok_all else 'SELF-CHECK FAIL: J6 control did NOT reproduce')
    return 0 if ok_all else 1


def table(out):
    rows = ['| state | depth | #eqs | #new unk | dim X_N | deg X_N | rel rank '
            '| obstruction rank | coker dim | lifts |',
            '|---|---|---|---|---|---|---|---|---|---|']
    for s in out['states']:
        run = next((r for r in s['runs'] if 'error' not in r), None)
        if run is None:
            continue
        steps = {x['step']: x for x in run['steps']}
        for d in run['depths']:
            N = d['depth']
            st = steps.get(f'{N}->{N + 1}')
            rows.append(
                f'| `{s["key"]}` | {N} | {d["n_equations"]} | '
                + (f'{st["n_new_unknowns"]} | ' if st else '- | ')
                + f'{d["dim"]} | {d["vdim"] if d["vdim"] >= 0 else "inf"} | '
                + (f'{st["rel_jac_rank"]} | ' if st else '- | ')
                + (f'{st.get("op_rank", "n/a")} | '
                   f'{st.get("op_cokernel_dim", "n/a")} | '
                   f'{"yes" if st["lifts_exist"] else "NO"} |' if st
                   else '- | - | - |'))
    return '\n'.join(rows)


if __name__ == '__main__':
    sys.exit(main())
