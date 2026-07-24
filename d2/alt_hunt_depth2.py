"""alt_hunt_depth2.py -- depth-2 residue attack on the 17 HUNT (BM-candidate) cells.

DIVISOR_LEMMAS.md sec.6 localised the 17 s_unit BM-candidate cells
(s_unit_results.json census.candidate_rows) to one depth-2 coefficient each:
the depth-1 (top-degree) S-unit relation is a satisfiable hypersurface, and the
completed kill needs the next master coefficient(s) under the defect-0 divisor
reconstruction.  This runner executes exactly that program, per fully-forced
state, using the audited machinery:

  - valuation splits: cascade_engine Pareto place profiles (the same machinery
    phase_f_defects.py used), generalised over the window (sub2 AND sub1 --
    phase_f2_sub2.place_val_options hardwires sub2);
  - master-identity walk + saturated-Groebner kill test: the phase_f2_sub2
    pattern (convolution_descent engine, top-down coefficients, unit ideal =>
    kill), extended with:
      (a) full audit trail: EVERY accumulated master coefficient recorded as
          an exact string so a spec-only auditor can re-derive the kill;
      (b) CLASS-POLYNOMIAL reconstruction: the four q-roots are grouped by
          joint exponent profile (v_e, v_d1, v_sig); each class of size n
          contributes the monic degree-n factor psi_C of q/2048 formed by
          its roots.  Every class polynomial except the largest is
          parameterized with unknown monic coefficients, the largest is the
          exact quotient, and the remainder coefficients of
          rem(q/2048, prod psi_i) = 0 are the defining relations.  The
          solution variety is exactly the set of root partitions matching
          the profile multiset (q squarefree => classes automatically
          disjoint), so a unit ideal kills EVERY Galois assignment at once.
          Unknowns = 4 - largest class size <= 3; every previously "heavy"
          pattern (unmarked-root valuations, b1100/b1110/b3110 multi-marked
          support, b1111 stable) reduces to this uniformly, over Q.
      (c) EXHAUSTIVE SPLIT DISJUNCTION: when a state's admissible split is
          not unique, the Pareto enumeration of splits is exhaustive, so the
          state dies iff EVERY split dies.  Splits equivalent under Galois
          permutation of same-profile roots are deduped to one class
          representative (sound by (b)).

Verdict vocabulary (per state, ALL PENDING AUDIT):
  KILLED       every admissible split reaches a unit saturated ideal
               (contradiction; generators recorded per split)
  CONSTRAINED  exact nonzero constraints accumulated, no unit ideal in budget
  OPEN         not decidable here (reason recorded): >= 3 distinct roots
               needed (b1110/b3110 e-support), timeouts.

Soundness notes (mirror PHASE_F2_SUB2.md):
  - the master identity f31 = sum_f Phi^f e^{21-3f} h_f == 0 is a necessary
    condition on any state, T1 or T2 (branch only chose the descent path);
  - only leading scalars (and (r-s) in two-root mode) are saturated; a unit
    ideal is a sound kill for the reconstructed defect-0 family;
  - a state is KILLED only if every admissible divisor split is killed.

Usage:  python alt_hunt_depth2.py [--quiet]
Env:    AH_BUDGET (GB budget s/split, default 60), AH_HARD (hard wall
        s/split, default 300), AH_MAXCOEFFS (depth cap, default 8),
        AH_REDO_OPEN=1 (reprocess prior OPEN states after an upgrade).
Output: alt_hunt_results.json (checkpointed per state, resumable).
"""
import json
import os
import sys
import time
import multiprocessing as mp
from itertools import product

import sympy as sp
import cascade_engine as ce
import convolution_descent as cd
import phase_f2_sub2 as f2

DEPTH = 4
y = cd.y
r = f2.r
s = sp.Symbol('s')
E, S, X, w = f2.E, f2.S, f2.X, f2.w
QY = sum(c * y ** (4 - i) for i, c in enumerate(f2.Q_COEFFS))  # q(y), lc 2048
QS = sum(c * s ** (4 - i) for i, c in enumerate(f2.Q_COEFFS))  # q(s)

RESULTS = 'alt_hunt_results.json'


# --------------------------------------------------------------------------
# window-generalised valuation splits (phase_f2_sub2.place_val_options is
# hardwired to config sub2; identical otherwise, r_cap = 10+3a both windows
# exactly as phase_f_defects.case_join_vectors)
# --------------------------------------------------------------------------
def place_val_options(case, win):
    a = case['a_t']
    b_vec = tuple(case['b'])
    branch = case['branch']
    sz, dz = case['sigma_zero'], case['d2_zero']
    gz = tuple(case['g_zero_levels'])
    config = ce.CONFIGS[win]
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
    opts.append(sorted(set(trip(p) for p in tp)))
    return opts


def relevant_dims(case):
    dims = []
    if case['branch'] != 'T2':
        dims.append(0)
    if not case['sigma_zero']:
        dims.append(1)
    if not case['d2_zero']:
        dims.append(2)
    return dims


def all_splits(case, st, win):
    """Every admissible simultaneous split of a defect-0 state's degrees,
    deduped (i) on the relevant (nonzero-poly) components and (ii) under
    Galois permutation of the e-UNMARKED roots (sound: unmarked roots are
    treated generically, only q(rho)=0 imposed).

    The Pareto profile enumeration is exhaustive over admissible per-place
    valuations, so this list covers every divisor the state could carry."""
    opts = place_val_options(case, win)
    st_d1 = f2.deg_val(st['deg_d1'])
    st_sig = f2.deg_val(st['deg_sigma'])
    st_d2 = f2.deg_val(st['deg_d2'])
    branch, sz, dz = case['branch'], case['sigma_zero'], case['d2_zero']
    tx = None if branch == 'T2' else st_d1
    tz = None if sz else st_sig
    tk = None if dz else st_d2
    sols = []
    for combo in product(*opts):
        if tx is not None and sum(c[0] for c in combo) != tx:
            continue
        if tz is not None and sum(c[1] for c in combo) != tz:
            continue
        if tk is not None and sum(c[2] for c in combo) != tk:
            continue
        sols.append(combo)
    dims = relevant_dims(case)
    b_vec = tuple(case['b'])
    seen, reduced = set(), []
    for combo in sols:
        rel = lambda c: tuple(c[i] for i in dims)
        # Galois signature: multiset of joint per-root profiles (e-valuation
        # + relevant divisor valuations) + the t place.  Roots with equal
        # profiles are interchangeable; the class-polynomial reconstruction
        # covers every conjugate assignment of a signature at once.
        prof = tuple(sorted((b_vec[j],) + rel(combo[j]) for j in range(4)))
        sig = (rel(combo[4]), prof)
        if sig not in seen:
            seen.add(sig)
            reduced.append(combo)
    return reduced


# --------------------------------------------------------------------------
# reconstruction via CLASS POLYNOMIALS.
#
# Group the four q-roots by their joint exponent profile (v_e, v_d1, v_sig).
# Each class C of size n contributes the monic degree-n factor psi_C of
# q/2048 formed by its roots.  Parameterize every class polynomial except the
# largest with unknown monic coefficients and derive the largest as the exact
# quotient  psi_big = quo(q/2048, prod psi_i);  the remainder coefficients
# rem(q/2048, prod psi_i) = 0 are the defining relations.  The solution
# variety is exactly the set of root partitions matching the profile multiset
# (q squarefree => classes automatically root-disjoint, no saturation
# needed), so a unit saturated ideal kills EVERY Galois assignment of the
# split at once.  Unknown count = 4 - (largest class size) <= 3.
# --------------------------------------------------------------------------
QM = sp.expand(QY / 2048)     # monic q, exact rational coefficients


def reconstruct_general(case, combo):
    """(polys, scalars, unknowns, relations) or (None, reason)."""
    a = case['a_t']
    b_vec = tuple(case['b'])
    branch = case['branch']
    sz, dz = case['sigma_zero'], case['d2_zero']
    if not dz:
        return None, 'd2_nonzero_unsupported'    # no HUNT cell carries d2
    vt_d1, vt_sig, vt_d2 = combo[4]

    # per-root joint profile (v_e, v_d1, v_sig); zero components for polys
    # that are identically zero in this cell
    prof = []
    for j in range(4):
        ve = b_vec[j]
        vd1 = 0 if branch == 'T2' else combo[j][0]
        vsig = 0 if sz else combo[j][1]
        prof.append((ve, vd1, vsig))
    classes = {}
    for p in prof:
        classes[p] = classes.get(p, 0) + 1

    # largest class is derived by quotient; the rest are parameterized
    items = sorted(classes.items(), key=lambda kv: -kv[1])
    big_prof, big_n = items[0]
    param = items[1:]
    unknowns = []
    psis = {}
    P = sp.Integer(1)
    for ci, (pr, n) in enumerate(param):
        cs = [sp.Symbol(f'c{ci}_{i}') for i in range(n)]
        unknowns.extend(cs)
        psi = y ** n + sum(cs[i] * y ** i for i in range(n))
        psis[pr] = sp.expand(psi)
        P = sp.expand(P * psi)
    quoP, remP = sp.div(QM, P, y)
    psis[big_prof] = sp.expand(quoP)
    relations = [c for c in sp.Poly(remP, y).all_coeffs() if c != 0]

    def contrib(dim):
        p = sp.Integer(1)
        for pr, psi in psis.items():
            if pr[dim]:
                p = p * psi ** pr[dim]
        return p

    e_expr = sp.expand(E * (y + 1) ** a * contrib(0))
    d1_expr = sp.Integer(0) if branch == 'T2' else \
        sp.expand(X * (y + 1) ** vt_d1 * contrib(1))
    sig_expr = sp.Integer(0) if sz else \
        sp.expand(S * (y + 1) ** vt_sig * contrib(2))
    scalars = [E]
    if branch != 'T2':
        scalars.append(X)
    if not sz:
        scalars.append(S)
    return (dict(d2=sp.Integer(0), d1=d1_expr, sigma=sig_expr, e=e_expr),
            scalars, unknowns, relations), 'ok'


# --------------------------------------------------------------------------
# kill test -- master-identity walk with full generator recording.
# --------------------------------------------------------------------------
def kill_test_record(polys, scalars, unknowns, relations, budget,
                     max_coeffs):
    ans = cd.build_ansatz(d2=polys['d2'], d1=polys['d1'], e=polys['e'],
                          sigma=polys['sigma'], parameters=tuple(unknowns))
    eng = cd.ConvolutionDescent(ans, c=f2.C_VAL)
    top = f2.engine_top(eng)
    order_vars = list(scalars) + list(unknowns) + [w]
    sat = sp.expand(w * sp.Mul(*scalars) - 1)
    gens, gen_strs = [], []
    t0 = time.time()
    for n in range(max_coeffs):
        target = top - n
        mc = sp.expand(eng.master_coefficient(target))
        if mc != 0:
            gens.append(mc)
            gen_strs.append({'degree': target, 'coefficient': str(mc)})
        if not gens:
            continue
        try:
            G = sp.groebner(gens + list(relations) + [sat], *order_vars,
                            order='grevlex')
        except Exception as ex:
            return {'verdict': 'ERROR', 'error': str(ex)[:200],
                    'depth': n + 1, 'top_degree': top, 'gens': gen_strs}
        if sp.Integer(1) in G.exprs:
            return {'verdict': 'KILLED', 'kill_depth': n + 1,
                    'top_degree': top, 'ngens': len(gens),
                    'elapsed': round(time.time() - t0, 1), 'gens': gen_strs,
                    'saturation': str(sat),
                    'class_relations': [str(x) for x in relations]}
        if time.time() - t0 > budget:
            return {'verdict': 'PENDING_TIMEOUT', 'depth': n + 1,
                    'top_degree': top,
                    'elapsed': round(time.time() - t0, 1), 'gens': gen_strs}
    return {'verdict': 'CONSTRAINED', 'depth': max_coeffs, 'top_degree': top,
            'ngens': len(gens), 'elapsed': round(time.time() - t0, 1),
            'gens': gen_strs}


def _worker_entry(q, polys_s, scalars_s, unknowns_s, relations_s, budget,
                  max_coeffs):
    try:
        polys = {k: sp.sympify(v) for k, v in polys_s.items()}
        scalars = [sp.sympify(x) for x in scalars_s]
        unknowns = [sp.sympify(x) for x in unknowns_s]
        relations = [sp.sympify(x) for x in relations_s]
        q.put(kill_test_record(polys, scalars, unknowns, relations, budget,
                               max_coeffs))
    except Exception as ex:                       # pragma: no cover
        q.put({'verdict': 'ERROR', 'error': str(ex)[:200]})


def kill_test_guarded(polys, scalars, unknowns, relations, budget,
                      max_coeffs, hard):
    ctx = mp.get_context('spawn')
    q = ctx.Queue()
    p = ctx.Process(target=_worker_entry,
                    args=(q, {k: sp.srepr(v) for k, v in polys.items()},
                          [sp.srepr(x) for x in scalars],
                          [sp.srepr(x) for x in unknowns],
                          [sp.srepr(x) for x in relations],
                          budget, max_coeffs))
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
# per-state driver: exhaustive disjunction over admissible splits
# --------------------------------------------------------------------------
def decide_state(case, st, win, budget, max_coeffs, hard):
    sols = all_splits(case, st, win)
    if not sols:
        return {'verdict': 'OPEN', 'reason': 'no_simultaneous_split',
                'n_splits': 0}
    per_split = []
    all_killed = True
    for combo in sols:
        entry = {'combo': [list(c) for c in combo]}
        recon, why = reconstruct_general(case, combo)
        if recon is None:
            entry['verdict'] = 'OPEN'
            entry['reason'] = why
            per_split.append(entry)
            all_killed = False
            continue
        polys, scalars, unknowns, relations = recon
        entry['n_class_unknowns'] = len(unknowns)
        entry['polys'] = {kk: str(vv) for kk, vv in polys.items()}
        res = kill_test_guarded(polys, scalars, unknowns, relations, budget,
                                max_coeffs, hard)
        entry.update(res)
        per_split.append(entry)
        if res['verdict'] != 'KILLED':
            all_killed = False
    if all_killed:
        verdict = 'KILLED'
    elif any(e['verdict'] in ('KILLED', 'CONSTRAINED') for e in per_split):
        verdict = 'CONSTRAINED'
    else:
        verdict = 'OPEN'
    out = {'verdict': verdict, 'n_splits': len(sols), 'splits': per_split,
           'mechanism': 'exhaustive_split_disjunction'}
    if verdict == 'KILLED':
        out['kill_depth'] = max(e.get('kill_depth') or 0 for e in per_split)
    else:
        open_reasons = sorted(set(e.get('reason') or e['verdict']
                                  for e in per_split
                                  if e['verdict'] != 'KILLED'))
        out['reason'] = ';'.join(open_reasons)[:200]
    return out


# --------------------------------------------------------------------------
# target enumeration: the 17 census cells, fully-forced (all-0) states only
# --------------------------------------------------------------------------
def all0_indices(defcase, branch):
    gkey = 'g7' if branch == 'T1' else 'g6'
    core = ('d1', 'sigma', 'd2', 'e')
    idxs = []
    for i, d in enumerate(defcase['state_deltas']):
        if '_error' in d:
            continue
        if all(abs(d.get(k, 0)) < 1e-9 for k in core) \
                and abs(d.get(gkey, 0)) < 1e-9:
            idxs.append(i)
    return idxs


def load_hunt_targets():
    rows = json.load(open('s_unit_results.json'))['census']['candidate_rows']
    defects = json.load(open('phase_f_defects.json'))['windows']
    dmap = {c['cellid']: c for wd in defects.values() for c in wd['cases']}
    states = {win: json.load(open(f'phase_d_states_{win}.json'))
              for win in ('sub2', 'sub1')}
    targets = []
    for row in rows:
        win = row['window']
        cellid = row['cellid']
        defcase = dmap[cellid]
        case = None
        for c in states[win]['cases']:
            cid = (f"{win}:a{c['a_t']}_b{''.join(map(str, c['b']))}"
                   f"_{c['branch']}_sz{int(c['sigma_zero'])}"
                   f"_dz{int(c['d2_zero'])}"
                   f"_gz{'.'.join(map(str, c['g_zero_levels'])) or '-'}")
            if cid == cellid:
                case = c
                break
        assert case is not None, cellid
        idxs = all0_indices(defcase, case['branch'])
        assert len(idxs) == row['n_states_all0'], \
            (cellid, len(idxs), row['n_states_all0'])
        targets.append((win, cellid, case, idxs))
    return targets


def main():
    quiet = '--quiet' in sys.argv
    budget = float(os.environ.get('AH_BUDGET', '60'))
    hard = float(os.environ.get('AH_HARD', '300'))
    max_coeffs = int(os.environ.get('AH_MAXCOEFFS', '8'))
    redo_open = os.environ.get('AH_REDO_OPEN') == '1'
    targets = load_hunt_targets()

    done = {}
    if os.path.exists(RESULTS):
        try:
            done = {rec['key']: rec
                    for rec in json.load(open(RESULTS)).get('states', [])}
        except Exception:
            done = {}

    results = []
    nstates = sum(len(idxs) for _, _, _, idxs in targets)
    k = 0
    for win, cellid, case, idxs in targets:
        for i in idxs:
            k += 1
            key = f'{cellid}#state{i}'
            if key in done and not (redo_open
                                    and done[key]['verdict'] != 'KILLED'):
                results.append(done[key])
                continue
            st = case['states'][i]
            rec = {'key': key, 'window': win, 'cellid': cellid,
                   'state_idx': i,
                   'degs': [st['deg_d1'], st['deg_sigma'], st['deg_d2'],
                            st['deg_e']]}
            rec.update(decide_state(case, st, win, budget, max_coeffs, hard))
            results.append(rec)
            _save(results)
            _log(quiet, k, nstates, rec)

    from collections import Counter
    cen = Counter(rec['verdict'] for rec in results)
    bycell = {}
    for rec in results:
        bycell.setdefault(rec['cellid'], []).append(rec['verdict'])
    cells_closed = sorted(c for c, v in bycell.items()
                          if all(x == 'KILLED' for x in v))
    summary = {'state_census': dict(cen),
               'n_states': len(results),
               'cells_closed_all0': cells_closed,
               'n_cells_closed_all0': len(cells_closed),
               'n_cells': len(bycell)}
    _save(results, summary)
    print('STATE CENSUS', dict(cen))
    print(f"HUNT CELLS with all fully-forced states KILLED: "
          f"{len(cells_closed)}/{len(bycell)}")
    for c in cells_closed:
        print('  CLOSED(all0):', c)
    return 0


def _save(results, summary=None):
    payload = {'schema': 2,
               'item': 'alt_hunt_depth2 -- 17 HUNT cells, depth-2 residue',
               'status': 'ALL VERDICTS PENDING AUDIT',
               'states': results}
    if summary:
        payload['summary'] = summary
    json.dump(payload, open(RESULTS, 'w'), indent=1)


def _log(quiet, k, n, rec):
    if quiet:
        return
    print(f"[{k}/{n}] {rec['key']} deg{rec['degs']} -> {rec['verdict']}"
          f" (kd={rec.get('kill_depth')}, nsplits={rec.get('n_splits')},"
          f" {rec.get('reason', '')})", flush=True)


if __name__ == '__main__':
    sys.exit(main())
