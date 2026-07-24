#!/usr/bin/env python3
"""j6_msolve.py -- msolve pass on the 4 OPEN [J6] states of ALT_HUNT.md.

Lane C (alt_hunt_depth2.py) killed 45/49 fully-forced HUNT states; the 4
survivors (2 windows x {a9_b1000_T1 deg-6-d1, a8_b1100_T1 deg-6-d1}) are the
known [J6] grevlex blowup: sympy's Buchberger exceeds the wall on a single GB
call.  BLOWUP_DIAGNOSIS.md's named cure for exactly this shape is msolve
(F4 + multi-modular + rational reconstruction over Q).  This runner attacks
EXACTLY the systems Lane C recorded:

  - the split combo and reconstructed polys are read VERBATIM from
    alt_hunt_results.json (never re-derived);
  - alt_hunt_depth2.reconstruct_general(case, combo) is REPLAYED and its
    polys asserted string-identical to the recorded ones (replay guarantee --
    same code path, same reconstruction, or loud abort);
  - where Lane C recorded accumulated master coefficients ('gens', the
    sub1:a9 state at depth 2), the replayed walk is asserted to reproduce
    those exact strings before anything is handed to msolve;
  - the msolve system per depth n is gens[..n] + class relations +
    Rabinowitsch saturation w*prod(scalars)-1 -- identical generator set to
    alt_hunt_depth2.kill_test_record, only the GB engine differs.

msolve output semantics (msolve_bridge.py): '[-1]' => empty over the
algebraic closure => KILLED (same soundness as the unit-ideal grevlex kill:
the saturated system has no solution, so no Galois assignment of the split
survives).  Solutions/positive-dim at depth n just mean "not yet dead at
depth n" -- we deepen (more master coefficients), NOT a survival signal.
Verdicts: KILLED(depth), CONSTRAINED_NOT_KILLED_DEPTH<k> (honest: depth cap
reached, deeper coefficients exist beyond), COST (msolve/walk timeout).

Orphan-proofing: msolve runs WSL-side under `timeout` + `ulimit -v 8G`, so
a Windows-relay kill cannot leave a live WSL process; a final WSL ps sweep
is printed.

New file; READ-ONLY on all existing artifacts.  Usage:
    python j6_msolve.py [--quiet]
Env: J6_STATE_CAP (per-state wall s, default 1200), J6_MAXCOEFFS (default 8).
Output: j6_msolve_results.json, plus j6_*.ms/.out staged in WSL $HOME.
"""
import json
import os
import re
import subprocess
import sys
import time

import sympy as sp

import alt_hunt_depth2 as ah
import blowup_diagnosis as bd

QUIET = '--quiet' in sys.argv
STATE_CAP = float(os.environ.get('J6_STATE_CAP', '1200'))
MAXCOEFFS = int(os.environ.get('J6_MAXCOEFFS', '8'))
MSOLVE = '$HOME/msolve/msolve'
RESULTS = 'j6_msolve_results.json'
w = ah.w


def log(*a):
    if not QUIET:
        print(*a, flush=True)


def wsl(cmd, timeout=None, inp=None):
    return subprocess.run(('wsl.exe', '-d', 'Ubuntu', '--', 'bash', '-lc', cmd),
                          input=inp, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, text=True, encoding='utf-8',
                          timeout=timeout, check=False)


def run_msolve(tag, gens, ring_vars, budget):
    """One msolve run, WSL-side timeout + 8G vcap.  Returns (verdict, info)."""
    var_line = ','.join(v.name for v in ring_vars)
    polys = [s for s in (bd.sing_poly_intcoeff(g, ring_vars) for g in gens)
             if s != '0']
    prog = var_line + '\n0\n' + ',\n'.join(polys) + '\n'
    fn, on = f'j6_{tag}.ms', f'j6_{tag}.out'
    wsl(f'cat > $HOME/{fn}', inp=prog)
    tsec = max(10, int(budget))
    t0 = time.monotonic()
    try:
        cp = wsl(f'cd $HOME && ulimit -v 8388608; timeout {tsec} {MSOLVE} '
                 f'-f $HOME/{fn} -o $HOME/{on}; echo "@@EXIT:$?"; '
                 f'cat $HOME/{on} 2>/dev/null',
                 timeout=tsec + 60)
    except subprocess.TimeoutExpired:
        return 'COST', {'note': 'relay timeout (WSL-side timeout pending)',
                        'wall': round(time.monotonic() - t0, 1)}
    wall = round(time.monotonic() - t0, 1)
    out = cp.stdout or ''
    mex = re.search(r'@@EXIT:(\d+)', out)
    exitc = int(mex.group(1)) if mex else None
    body = out.split('@@EXIT:', 1)[-1]
    body = body.split('\n', 1)[1].strip() if '\n' in body else ''
    info = {'wall': wall, 'exit': exitc, 'ms_file': fn,
            'out_head': body[:160], 'nvars': len(ring_vars),
            'ngens': len(polys)}
    if exitc == 124:
        return 'COST', info
    if body.startswith('[-1]'):
        return 'EMPTY', info
    if re.match(r'\[1,\s*\d+,\s*-1', body):
        return 'NOT_EMPTY', dict(info, kind='positive_dim')
    if body:
        return 'NOT_EMPTY', dict(info, kind='parametrization')
    return 'ERROR', info


def attack_state(st, case):
    key = st['key']
    split = st['splits'][0]
    assert st['n_splits'] == 1, (key, st['n_splits'])
    combo = tuple(tuple(c) for c in split['combo'])
    t_state = time.monotonic()

    # --- replay guarantee: same code path, string-identical polys ----------
    recon, why = ah.reconstruct_general(case, combo)
    if recon is None:
        return {'key': key, 'verdict': 'ERROR',
                'error': f'reconstruct_general refused: {why}'}
    polys, scalars, unknowns, relations = recon
    replayed = {k: str(v) for k, v in polys.items()}
    if replayed != split['polys']:
        return {'key': key, 'verdict': 'ERROR',
                'error': 'REPLAY MISMATCH: reconstructed polys differ from '
                         'recorded strings',
                'diff': {k: (split['polys'].get(k, '<absent>')[:80],
                             replayed.get(k, '<absent>')[:80])
                         for k in set(replayed) | set(split['polys'])
                         if replayed.get(k) != split['polys'].get(k)}}

    # --- engine, exactly as kill_test_record --------------------------------
    ans = ah.cd.build_ansatz(d2=polys['d2'], d1=polys['d1'], e=polys['e'],
                             sigma=polys['sigma'],
                             parameters=tuple(unknowns))
    eng = ah.cd.ConvolutionDescent(ans, c=ah.f2.C_VAL)
    top = ah.f2.engine_top(eng)
    sat = sp.expand(w * sp.Mul(*scalars) - 1)
    ring_vars = sorted(set(scalars) | set(unknowns) | {w},
                       key=lambda v: v.name)
    recorded = {g['degree']: g['coefficient'] for g in split.get('gens', [])}

    gens, gen_strs, attempts = [], [], []
    verdict, kill_depth = None, None
    for n in range(MAXCOEFFS):
        target = top - n
        mc = sp.expand(eng.master_coefficient(target))
        if mc != 0:
            if target in recorded:            # replay vs Lane C's record
                if str(mc) != recorded[target] and \
                        sp.expand(mc - sp.sympify(recorded[target])) != 0:
                    return {'key': key, 'verdict': 'ERROR',
                            'error': f'REPLAY MISMATCH: master coefficient '
                                     f'at degree {target} differs from '
                                     f'recorded string'}
            gens.append(mc)
            gen_strs.append({'degree': target, 'coefficient': str(mc)})
        if not gens or n == 0:
            continue
        remaining = STATE_CAP - (time.monotonic() - t_state)
        if remaining < 15:
            verdict = 'COST'
            attempts.append({'depth': n + 1, 'skipped': 'state cap reached'})
            break
        tag = key.replace(':', '_').replace('#', '_') + f'_d{n + 1}'
        v, info = run_msolve(tag, gens + list(relations) + [sat], ring_vars,
                             min(remaining - 10, STATE_CAP))
        attempts.append(dict(info, depth=n + 1, msolve=v))
        log(f'    depth {n + 1}: msolve {v} ({info.get("wall")}s)')
        if v == 'EMPTY':
            verdict, kill_depth = 'KILLED', n + 1
            break
        if v == 'COST':
            verdict = 'COST'
            break
        if v == 'ERROR':
            verdict = 'ERROR'
            break
        # NOT_EMPTY at truncated depth: deepen (not a survival signal)
    else:
        verdict = f'CONSTRAINED_NOT_KILLED_DEPTH{MAXCOEFFS}'

    out = {'key': key, 'window': st['window'], 'cellid': st['cellid'],
           'state_idx': st['state_idx'], 'combo': split['combo'],
           'n_class_unknowns': split['n_class_unknowns'],
           'verdict': verdict, 'top_degree': top,
           'ring_vars': [v.name for v in ring_vars],
           'saturation': str(sat),
           'class_relations': [str(x) for x in relations],
           'gens': gen_strs, 'msolve_attempts': attempts,
           'elapsed': round(time.monotonic() - t_state, 1)}
    if kill_depth:
        out['kill_depth'] = kill_depth
    return out


def main():
    res = json.load(open('alt_hunt_results.json', encoding='utf-8'))
    open_states = [s for s in res['states'] if s['verdict'] == 'OPEN']
    log(f'{len(open_states)} OPEN states to attack')
    assert len(open_states) == 4, [s['key'] for s in open_states]

    targets = {cellid: (win, case, idxs)
               for win, cellid, case, idxs in ah.load_hunt_targets()}
    out = {'schema': 1, 'item': 'j6_msolve pass on ALT_HUNT OPEN states',
           'source': 'alt_hunt_results.json (systems replayed verbatim)',
           'state_cap_s': STATE_CAP, 'maxcoeffs': MAXCOEFFS, 'results': []}
    for st in open_states:
        log(f'== {st["key"]}')
        _, case, _ = targets[st['cellid']]
        rec = attack_state(st, case)
        out['results'].append(rec)
        log(f'   -> {rec["verdict"]}'
            + (f' (depth {rec["kill_depth"]})' if rec.get('kill_depth') else ''))
        json.dump(out, open(RESULTS, 'w'), indent=1)

    # orphan sweep
    ps = wsl("ps -eo pid,etime,comm | grep -E 'msolve|Singular' | grep -v grep"
             " ; true", timeout=60)
    out['wsl_ps_after'] = (ps.stdout or '').strip() or 'clean'
    json.dump(out, open(RESULTS, 'w'), indent=1)
    log('WSL ps after:', out['wsl_ps_after'])

    census = {}
    for r_ in out['results']:
        census[r_['verdict']] = census.get(r_['verdict'], 0) + 1
    log('CENSUS:', census)
    return 0


if __name__ == '__main__':
    sys.exit(main())
