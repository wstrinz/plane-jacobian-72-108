#!/usr/bin/env python3
"""bridge_sweep.py -- sweep the FULL-SYSTEM BRIDGE across the resistant core.

This is the endgame-tool sweep asked for in the mission brief. It reuses the
landed, verified bridge (``full_system_bridge.augment`` -- the G-system
augmentation validated in FULL_SYSTEM_BRIDGE.md) and EXTENDS it in exactly one
place: the marked-root sub2 T2 states (the R9 column and the a7/a8 T2 cells)
need the exact quotient ring Q[r]/(q(r)) and a saturation by the state's genuine
nonzero scalars (gamma, lc(G), G(r)). The core bridge has no notion of a marked
root, so this file builds that wrapper, documented below, mirroring the pilot's
stripping discipline exactly (regime=sub2 stripped ansatz, phi = c t^30 q).

CONSTRUCTION (marked-root R9 / T2 pattern-B bridge)
---------------------------------------------------
For an R9 state (a9 b1000 T2 pattern-B, z=0..6):

    e     = gamma (y+1)^9 (y-r),   d1 = 0,
    sigma = (y-r)^2 G(y),  deg G = z,   deg d2 <= 4,   d0 = (d2^2+sigma)/4

built by ``convolution_elim_qsupport.build_qsupport_ansatz(z)`` -- the LANDED,
audited R9 construction (its check_valuations gate proves v_r(e)=1, v_r(sigma)=2).
This is already a cascade *stripped* ansatz that sits exactly at the sub2 window
caps (deg e=10 = 2*5, deg sigma<=8 = 2*4, deg d2=4 = 2*2), so ``augment`` applies
verbatim: it substitutes (d2,d1,d0,e) into G1,G2,G3,(G5body+Phi), introduces the
spare window unknowns dm2,dm3,dm4 as bounded stripped polynomials, and returns
every y-coefficient as an equation. r rides along symbolically inside those
coefficients.

The two marked-root additions (the ONLY departure from the pilot):
  1. adjoin q(r) = 0  (the exact Q[r]/(q) representation -- r a ring variable);
  2. saturate by (gamma, lc(G), G(r)) via a single Rabinowitsch w*prod-1 (the
     genuine nonzero scalars: gamma!=0 the gauge, lc(G)!=0 so deg G is exact,
     G(r)!=0 so the marked root has multiplicity exactly 2 in sigma).
Both are taken verbatim from the landed R9 saturation
(``build_qsupport_ansatz.saturation_factors`` and Q_R), and from the exact-Q
harvest recipe (TRIAGE_HARVEST.md). No coefficient is emitted as a fraction.

Verdicts: mod-p triage over 3 good primes (10007,10009,100019 -- the triage
primes, all avoiding the bad set {2,3,5,13,17}); exact UNIT over Q => KILL
(PENDING AUDIT). A UNIT G-system means the germ does NOT lift through the full
window system -- a strictly stronger statement than the f31-alone / master-
coefficient kill (the bridge is >= the cascade, F37_SATURATION_REPORT fact [5]).

New file, uncommitted. READ-ONLY on every imported module/artifact.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import sympy as sp

import full_system_bridge as fsb
import convolution_elim_qsupport as qs
import convolution_descent as cd
import modular_triage as mt

ROOT = Path(__file__).resolve().parent
y = fsb.y
PRIMES = (10007, 10009, 100019)


# --------------------------------------------------------------------------
#  ring-var collection and Singular emission (marked-root aware)
# --------------------------------------------------------------------------
def _ring_vars(eqs, sat_factors, extra=()):
    syms = set(extra)
    for e in list(eqs) + list(sat_factors):
        syms |= sp.sympify(e).free_symbols
    syms.discard(y)
    return sorted(syms, key=sp.default_sort_key)


def _emit(eqs, sat_factors, ring_vars, *, char):
    """Integer-cleared (Q) or mod-p Singular program that saturates by the
    PRODUCT of sat_factors (Rabinowitsch w) and unit-tests the ideal."""
    w = sp.Symbol("w")
    rv = list(ring_vars) + [w]
    var_txt = ",".join(v.name for v in rv)
    lines = ['LIB "elim.lib";', f"ring R = {char},({var_txt}),dp;"]
    members = []
    for i, g in enumerate(eqs):
        if char > 0:
            s = mt.poly_to_singular_modp(g, rv, char)
        else:
            s = fsb._to_singular(g, rv)
        if s in ("0", ""):
            continue
        lines.append(f"poly g{i} = {s};")
        members.append(f"g{i}")
    prod = sp.Integer(1)
    for f in sat_factors:
        prod = prod * f
    nz = sp.expand(w * prod - 1)
    nzs = (mt.poly_to_singular_modp(nz, rv, char) if char > 0
           else fsb._to_singular(nz, rv))
    lines.append(f"poly nz = {nzs};")
    if not members:
        members = ["0"]
    lines.append(f"ideal I = {','.join(members)},nz;")
    lines.append("ideal G = std(I);")
    lines.append("int u = (reduce(1,G)==0);")
    lines.append('"@@UNIT";')
    lines.append("u;")
    lines.append('"@@DIM";')
    lines.append("dim(G);")
    lines.append("quit;")
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------
#  R9 marked-root bridge state
# --------------------------------------------------------------------------
import pickle

CACHE = ROOT / "__pycache__"


def build_r9_bridge(z, *, regime="sub2", use_cache=True):
    # cache is this script's OWN self-generated build output (trusted local
    # data), only to skip the ~200s symbolic G5+Phi expansion on re-runs.
    cache_f = CACHE / f"bridge_r9_z{z}_{regime}.pkl"
    if use_cache and cache_f.exists():
        d = pickle.loads(cache_f.read_bytes())
        return d
    st = qs.build_qsupport_ansatz(z)
    aug = fsb.augment(st.ansatz, regime=regime)
    eqs = [e for e in aug["equations"] if e != 0]
    eqs.append(sp.expand(qs.Q_R))                 # adjoin q(r)=0
    sat_factors = list(st.saturation_factors)     # (gamma, lc(G), G(r))
    ring_vars = _ring_vars(eqs, sat_factors)
    d = {
        "label": f"R9_z{z}", "z": z, "regime": regime,
        "equations": eqs, "sat_factors": sat_factors, "ring_vars": ring_vars,
        "n_equations": len(eqs), "n_unknowns": len(ring_vars),
        "spare_unknowns": len(aug["spare_unknowns"]),
    }
    try:
        cache_f.write_bytes(pickle.dumps(d))
    except Exception:
        pass
    return d


def presolve_bridge(bs):
    """Sound linear pre-elimination (landed fsb.linear_presolve): solve every
    equation linear in one unknown with a nonzero CONSTANT pivot and substitute
    everywhere -- an equivalent rewriting of the ideal that removes most of the
    45 spare scalar unknowns before Singular, killing the coefficient swell.
    The substitution map is applied to the marked-root saturation factors too."""
    eqs, remaining, subs = fsb.linear_presolve(bs["equations"], bs["ring_vars"])
    sat = [sp.expand(sp.sympify(f).subs(subs)) for f in bs["sat_factors"]]
    rv = _ring_vars(eqs, sat)
    out = dict(bs)
    out.update({"equations": eqs, "sat_factors": sat, "ring_vars": rv,
                "n_equations": len(eqs), "n_unknowns": len(rv),
                "eliminated": len(subs)})
    return out


def _run(prog, timeout):
    return mt.run_singular(prog, timeout=timeout)


# --------------------------------------------------------------------------
#  msolve fallback for coefficient swell over Q (F4 + multi-modular + rational
#  reconstruction).  Emits MY OWN augmented system (q(r) adjoined, leading
#  scalars Rabinowitsch-saturated) -- does not touch the blowup lane's CASES.
#  msolve output [-1] == empty variety over Qbar == KILL.
# --------------------------------------------------------------------------
_MS = "$HOME/msolve/msolve"


def _emit_msolve(eqs, sat_factors, ring_vars, char=0):
    w = sp.Symbol("w")
    rv = list(ring_vars) + [w]
    prod = sp.Integer(1)
    for f in sat_factors:
        prod *= f
    gens = list(eqs) + [sp.expand(w * prod - 1)]
    polys = []
    for g in gens:
        s = fsb._to_singular(g, rv)   # exact integer-cleared string, ^ and *
        if s not in ("0", ""):
            polys.append(s)
    if not polys:
        polys = ["0"]
    var_line = ",".join(v.name for v in rv)
    return var_line + "\n" + str(char) + "\n" + ",\n".join(polys) + "\n"


def exact_bridge_msolve(bs, *, timeout=300.0, char=0):
    import subprocess
    prog = _emit_msolve(bs["equations"], bs["sat_factors"], bs["ring_vars"], char)
    tag = bs["label"].replace(" ", "")
    fname = f"bridge_{tag}.ms"
    outname = f"bridge_{tag}.out"
    subprocess.run(("wsl.exe", "-d", "Ubuntu", "--", "bash", "-lc",
                    f"cat > $HOME/{fname}"), input=prog, text=True,
                   encoding="utf-8", check=True)
    run = ("wsl.exe", "-d", "Ubuntu", "--", "bash", "-lc",
           f"cd $HOME && {_MS} -f $HOME/{fname} -o $HOME/{outname}; "
           f"cat $HOME/{outname}")
    t0 = time.monotonic()
    try:
        cp = subprocess.run(run, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            text=True, encoding="utf-8", timeout=timeout,
                            check=False)
    except subprocess.TimeoutExpired:
        return {"verdict": "TIMEOUT", "wall": round(time.monotonic() - t0, 1),
                "engine": "msolve"}
    body = (cp.stdout or "").strip()
    if body.startswith("[-1]"):
        verdict = "UNIT"        # empty == KILL (align with Singular UNIT label)
    elif __import__("re").match(r"\[1,\s*\d+,\s*-1", body):
        verdict = "INFINITE_SOL"
    elif body:
        verdict = "HAS_SOL"
    else:
        verdict = "NO_OUTPUT"
    return {"verdict": verdict, "wall": round(time.monotonic() - t0, 1),
            "engine": "msolve", "out_head": body[:120]}


def triage_bridge(bs, *, primes=PRIMES, timeout=35.0):
    out = []
    for p in primes:
        prog = _emit(bs["equations"], bs["sat_factors"], bs["ring_vars"], char=p)
        rr = _run(prog, timeout)
        rr["prime"] = p
        out.append(rr)
        print(f"    p={p}: {rr['verdict']} dim={rr.get('dim')} "
              f"({rr.get('wall')}s)", flush=True)
    return {"primes": out, "prediction": mt.classify(out)}


def triage_bridge_numroot(bs, *, primes=PRIMES, timeout=90.0):
    """Fast mod-p reconnaissance: specialize the marked root r to a numeric root
    of q mod p (the modular_triage.build_system1 route) instead of carrying r
    symbolically.  This removes the r-nonlinearity that swells the symbolic std,
    so the pure G-system verdict is reachable in seconds.  Evidence over F_p-bar,
    not a Q certificate."""
    r = sp.Symbol("r")
    eqs0 = [e for e in bs["equations"] if e != sp.expand(qs.Q_R) and e != 0]
    sat0 = bs["sat_factors"]
    out = []
    for p in primes:
        roots = mt.q_roots_mod_p(p, 1)
        if not roots:
            out.append({"verdict": "SKIP", "prime": p}); continue
        subst = {r: roots[0]}
        eqs = [sp.sympify(e).subs(subst) for e in eqs0]
        sat = [sp.sympify(f).subs(subst) for f in sat0]
        rv = _ring_vars(eqs, sat)
        prog = _emit(eqs, sat, rv, char=p)
        rr = _run(prog, timeout)
        rr["prime"] = p
        out.append(rr)
        print(f"    (numroot) p={p}: {rr['verdict']} dim={rr.get('dim')} "
              f"({rr.get('wall')}s)", flush=True)
    return {"primes": out, "prediction": mt.classify(out)}


def exact_bridge_minpoly(bs, *, timeout=300.0):
    """Exact run over the NUMBER FIELD Q(r) = Q[r]/(q) via Singular's native
    minpoly: r becomes a coefficient-field parameter (exact algebraic
    arithmetic), q(r) leaves the ideal, the ring shrinks by one variable.
    Mathematically identical to adjoining q(r) (q is irreducible over Q), but
    std now does field arithmetic in r -- the anti-swell formulation."""
    r = sp.Symbol("r")
    qr = sp.expand(qs.Q_R)
    eqs = [e for e in bs["equations"] if e != 0 and sp.expand(e - qr) != 0]
    sat = bs["sat_factors"]
    w = sp.Symbol("w")
    rv = [v for v in bs["ring_vars"] if v != r] + [w]
    var_txt = ",".join(v.name for v in rv)
    lines = ['LIB "elim.lib";',
             f"ring R = (0,r),({var_txt}),dp;",
             "minpoly = 2048*r^4-512*r^3+320*r^2-240*r+195;"]
    members = []
    for i, g in enumerate(eqs):
        s = fsb._to_singular(g, rv + [r])
        if s in ("0", ""):
            continue
        lines.append(f"poly g{i} = {s};")
        members.append(f"g{i}")
    prod = sp.Integer(1)
    for f in sat:
        prod *= f
    nz = fsb._to_singular(sp.expand(w * prod - 1), rv + [r])
    lines.append(f"poly nz = {nz};")
    if not members:
        members = ["0"]
    lines.append(f"ideal I = {','.join(members)},nz;")
    lines.append("ideal G = std(I);")
    lines.append("int u = (reduce(1,G)==0);")
    lines.append('"@@UNIT";')
    lines.append("u;")
    lines.append('"@@DIM";')
    lines.append("dim(G);")
    lines.append("quit;")
    prog = "\n".join(lines) + "\n"
    rr = _run(prog, timeout)
    rr["engine"] = "singular-minpoly"
    print(f"    exact-Q(r) minpoly: {rr['verdict']} dim={rr.get('dim')} "
          f"({rr.get('wall')}s)", flush=True)
    return rr


def exact_bridge(bs, *, timeout=300.0):
    prog = _emit(bs["equations"], bs["sat_factors"], bs["ring_vars"], char=0)
    rr = _run(prog, timeout)
    print(f"    exact-Q: {rr['verdict']} dim={rr.get('dim')} ({rr.get('wall')}s)",
          flush=True)
    return rr


# --------------------------------------------------------------------------
#  runner
# --------------------------------------------------------------------------
def run_r9_column(zs=range(7), *, exact_timeout=300.0):
    records = []
    for z in zs:
        t0 = time.monotonic()
        bs = build_r9_bridge(z)
        print(f"\nR9 z={z}: {bs['n_equations']} eqs, {bs['n_unknowns']} vars "
              f"(spare {bs['spare_unknowns']}); build {time.monotonic()-t0:.1f}s",
              flush=True)
        tri = triage_bridge(bs)
        print(f"  mod-p: {tri['prediction']}", flush=True)
        # The exact-Q Rabinowitsch std is the real certificate and (z=0) is FAST
        # even where mod-p std times out, so attempt it unless mod-p positively
        # says the ideal is PROPER on every prime (LIKELY-SOLVABLE).
        exact = None
        if tri["prediction"] != "LIKELY-SOLVABLE":
            exact = exact_bridge(bs, timeout=exact_timeout)
        verdict = "COST"
        if exact is not None and exact.get("verdict") == "UNIT":
            verdict = "KILLED"
        elif exact is not None and exact.get("verdict") == "PROPER":
            verdict = "PROPER"
        elif tri["prediction"] == "LIKELY-SOLVABLE":
            verdict = "PROPER(modp)"
        rec = {"label": bs["label"], "z": z, "n_equations": bs["n_equations"],
               "n_unknowns": bs["n_unknowns"], "spare_unknowns": bs["spare_unknowns"],
               "modp": tri, "exact": exact, "verdict": verdict}
        records.append(rec)
        print(f"  ==> R9 z={z}: {verdict}", flush=True)
        p = ROOT / "bridge_sweep.json"
        # re-read before each write: concurrent lanes checkpoint the same file
        prior = json.load(open(p)) if p.exists() else {}
        prior = prior if isinstance(prior, dict) else {}
        prior["R9_column"] = records
        json.dump(prior, open(p, "w"), indent=1, default=str)
    return records


# --------------------------------------------------------------------------
#  Target 2: the a8 constant-E constant states (no marked root; gamma-sat only).
#  These use the CORE bridge verbatim -- fsb.augment / fsb.exact_kill -- since
#  there is no q(r).  The pilot did index 0 (dsig5) and index 2 (dsig7); this
#  sweeps the whole 24, especially the deg_sigma 7/8 class that f31-deep could
#  not reach an exact Q certificate for (TRIAGE_HARVEST.md Target 4).
# --------------------------------------------------------------------------
def _a8_states():
    d = json.load(open(ROOT / "batch_convolution_sub2.json"))
    return [s for s in d["states"]
            if s["a_t"] == 8 and s["branch"] == "T1"
            and s["deg_e"] == 8 and s["final_verdict"] == "UNRESOLVED"]


def run_a8_column(indices=None, *, regime="sub2", exact_timeout=300.0):
    states = _a8_states()
    n = len(states)
    idxs = list(indices) if indices is not None else list(range(n))
    records = []
    prior = []
    p = ROOT / "bridge_sweep.json"
    if p.exists():
        prior = json.load(open(p))
    prior = prior if isinstance(prior, dict) else {}
    for index in idxs:
        s = states[index]
        t0 = time.monotonic()
        ans, meta, label = fsb.pilot_state(index)
        aug = fsb.augment(ans, regime=regime)
        dsig = int(meta["deg_sigma"])
        print(f"\na8 idx={index} {label} dsig={dsig}: {aug['n_equations']} eqs, "
              f"{aug['n_unknowns']} vars; build {time.monotonic()-t0:.1f}s",
              flush=True)
        tri = fsb.triage(aug)
        print(f"  mod-p: {tri['prediction']}", flush=True)
        exact = None
        if tri["prediction"] == "LIKELY-EMPTY":
            exact = fsb.exact_kill(aug, timeout=exact_timeout)
            print(f"    exact-Q: {exact['verdict']} ({exact.get('wall')}s)",
                  flush=True)
        verdict = ("KILLED" if exact and exact.get("verdict") == "UNIT"
                   else "PROPER" if exact and exact.get("verdict") == "PROPER"
                   else "PROPER(modp)" if tri["prediction"] == "LIKELY-SOLVABLE"
                   else "COST")
        rec = {"label": label, "index": index, "deg_sigma": dsig,
               "deg_d1": int(meta["deg_d1"]), "d2_zero": meta["d2_zero"],
               "n_equations": aug["n_equations"], "n_unknowns": aug["n_unknowns"],
               "modp": tri, "exact": exact, "verdict": verdict}
        records.append(rec)
        print(f"  ==> a8 idx={index}: {verdict}", flush=True)
        # re-read before each write: concurrent lanes checkpoint the same file
        prior = json.load(open(p)) if p.exists() else {}
        prior = prior if isinstance(prior, dict) else {}
        prior["a8_column"] = records
        json.dump(prior, open(p, "w"), indent=1, default=str)
    return records


# --------------------------------------------------------------------------
#  Target 4 (in-scope part): remaining UNRESOLVED classes of the sub2 batch
#  census, cheapest-first sample.  These are STANDARD-window sub2 states
#  (e = (y+1)^a_t * generic tail, gauge lc = gamma) -- bridge-applicable,
#  unlike the alt-regime System-3 states (see BRIDGE_SWEEP.md section 3).
#  Excludes the a8 T1 deg_e=8 class (Target 2, swept in full).
# --------------------------------------------------------------------------
def _batch_unresolved_sample(limit=20):
    d = json.load(open(ROOT / "batch_convolution_sub2.json"))
    out = []
    for s in d["states"]:
        if s["final_verdict"] != "UNRESOLVED":
            continue
        if s["a_t"] == 8 and s["branch"] == "T1" and s["deg_e"] == 8:
            continue  # Target 2
        m = int(s["deg_e"]) - int(s["a_t"])
        def _d(v):
            return 0 if v in ("-inf", None) else int(v)
        cost = (m, _d(s["deg_sigma"]) + _d(s["deg_d1"]) + _d(s["deg_d2"]),
                int(s["deg_e"]))
        out.append((cost, s))
    out.sort(key=lambda t: t[0])
    return [s for _c, s in out[:limit]]


def _batch_ansatz(s):
    """Rebuild the batch census gauged ansatz: e = (y+1)^a_t * tail,
    lc(tail) frozen to the nonzero parameter gamma (the census gauge)."""
    a_t, m = int(s["a_t"]), int(s["deg_e"]) - int(s["a_t"])
    gamma = sp.Symbol("gamma", nonzero=True)
    gs = tuple(sp.symbols(f"g0:{m}")) + (gamma,)
    e_expr = sp.expand((y + 1) ** a_t * sum(g * y**i for i, g in enumerate(gs)))
    zero, degrees = [], {}
    if s["d2_zero"]:
        zero.append("d2")
    else:
        degrees["d2"] = int(s["deg_d2"])
    if s["d1_zero"]:
        zero.append("d1")
    else:
        degrees["d1"] = int(s["deg_d1"])
    if s.get("sigma_zero"):
        zero.append("sigma")
        sig_arg = sp.Integer(0)
    else:
        degrees["sigma"] = int(s["deg_sigma"])
        sig_arg = None
    return cd.build_ansatz(e=e_expr, degrees=degrees, zero=zero,
                           sigma=sig_arg, parameters=(gamma,))


def run_batch_sample(limit=20, *, regime="sub2", exact_timeout=300.0):
    sample = _batch_unresolved_sample(limit)
    records = []
    p = ROOT / "bridge_sweep.json"
    prior = json.load(open(p)) if p.exists() else {}
    prior = prior if isinstance(prior, dict) else {}
    for s in sample:
        label = (f"batch_a{s['a_t']}{s['branch']}_e{s['deg_e']}"
                 f"_d2{s['deg_d2']}_d1{s['deg_d1']}_sig{s['deg_sigma']}")
        t0 = time.monotonic()
        try:
            ans = _batch_ansatz(s)
            aug = fsb.augment(ans, regime=regime)
        except Exception as ex:
            records.append({"label": label, "verdict": "BUILD_ERROR",
                            "error": str(ex)[:200]})
            print(f"\n{label}: BUILD_ERROR {ex}", flush=True)
            continue
        print(f"\n{label}: {aug['n_equations']} eqs, {aug['n_unknowns']} vars; "
              f"build {time.monotonic()-t0:.1f}s", flush=True)
        tri = fsb.triage(aug)
        print(f"  mod-p: {tri['prediction']}", flush=True)
        exact = None
        if tri["prediction"] != "LIKELY-SOLVABLE":
            exact = fsb.exact_kill(aug, timeout=exact_timeout)
            print(f"    exact-Q: {exact['verdict']} ({exact.get('wall')}s)",
                  flush=True)
        verdict = ("KILLED" if exact and exact.get("verdict") == "UNIT"
                   else "PROPER" if exact and exact.get("verdict") == "PROPER"
                   else "PROPER(modp)" if tri["prediction"] == "LIKELY-SOLVABLE"
                   else "COST")
        rec = {"label": label, "state": {k: s[k] for k in
               ("a_t", "branch", "deg_e", "deg_d2", "deg_d1", "deg_sigma",
                "d2_zero", "d1_zero", "sigma_zero", "tier")},
               "n_equations": aug["n_equations"], "n_unknowns": aug["n_unknowns"],
               "modp": tri, "exact": exact, "verdict": verdict}
        records.append(rec)
        print(f"  ==> {label}: {verdict}", flush=True)
        # re-read before each write: concurrent lanes checkpoint the same file
        prior = json.load(open(p)) if p.exists() else {}
        prior = prior if isinstance(prior, dict) else {}
        prior["batch_sample"] = records
        json.dump(prior, open(p, "w"), indent=1, default=str)
    return records


if __name__ == "__main__":
    args = sys.argv[1:]
    if args and args[0] == "r9":
        zs = [int(x) for x in args[1:]] if len(args) > 1 else range(7)
        run_r9_column(zs)
    elif args and args[0] == "a8":
        idxs = [int(x) for x in args[1:]] if len(args) > 1 else None
        run_a8_column(idxs)
    elif args and args[0] == "batch":
        lim = int(args[1]) if len(args) > 1 else 20
        run_batch_sample(lim)
    else:
        run_r9_column()
