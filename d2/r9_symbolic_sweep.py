#!/usr/bin/env python3
"""r9_symbolic_sweep.py -- per-state driver for the dm4-ELIMINATED reduced
bridge system (r9_symbolic_elim.py), attacking the deg_e = 10 cost wall:
the R9 column z >= 1 and the 90 unswept deg_e = 10 T2 batch states.

WHAT IS TESTED PER STATE (all sound necessary conditions; see
r9_symbolic_elim.py for the certificates):

    [ every y-coefficient of H2, H3, H5 on the state's stripped ansatz,
      with dm2, dm3 as bounded stripped polynomials (sub2 caps 12 / 14) ]
  + [ rem(dm2*dm3, monic(e), y) == 0  coefficient-wise   (G1 divisibility) ]
  + [ q(r) = 0 and the state's nonzero saturations, marked-root states ]

dm4's 17 spare scalar unknowns NEVER appear: they were eliminated once,
symbolically, before grounding.  VERDICT SEMANTICS (asymmetric, honest):
  UNIT over Q      -> state KILLED (candidate kill, PENDING AUDIT)
  PROPER anywhere  -> INCONCLUSIVE (the reduced system is weaker than the
                      full bridge; never a survival signal)
  TIMEOUT          -> COST (pure Groebner cost, no information)

Reuses bridge_sweep's emission / triage / exact-Q machinery verbatim (same
integer-cleared Singular, same Rabinowitsch saturation, same minpoly and
numroot formulations) by presenting the same state-dict shape.  The ONLY
runtime substitution: this module's Singular runner wraps the WSL call in
`timeout` + `ulimit -v` INSIDE WSL, so a Windows-side timeout can never
orphan a WSL Singular process (see cascade-engine-long-runs memory).

New file, uncommitted.  READ-ONLY on every landed module/artifact; does not
touch kill_certificate_tools.py / audit_gb_kills.py / kill_manifest.json.
"""
from __future__ import annotations

import json
import pickle
import re
import subprocess
import sys
import time
from pathlib import Path

import sympy as sp

import full_system_bridge as fsb
import bridge_sweep as bsw
import convolution_elim_qsupport as qs
import convolution_descent as cd
import modular_triage as mt

ROOT = Path(__file__).resolve().parent
CACHE = ROOT / "__pycache__"
OUT = ROOT / "r9_symbolic_sweep.json"
ELIM = ROOT / "r9_eliminated_system.json"

y = fsb.y
R = sp.Symbol("r")
PRIMES = (10007, 10009, 100019)
ULIMIT_KB = 8_000_000          # 8 GB virtual-memory cap per Singular process
CAPS = fsb.STRIP_DEGCAP["sub2"]  # dm2 <= 12, dm3 <= 14 (all targets are sub2)


# --------------------------------------------------------------------------
#  orphan-proof Singular runner (WSL-side timeout + ulimit); same parse as
#  modular_triage.run_singular.  Installed over bridge_sweep._run so every
#  reused triage/exact routine goes through it.
# --------------------------------------------------------------------------
def run_singular_guarded(program: str, timeout: float = 60.0) -> dict:
    t0 = time.monotonic()
    inner = (f"cd $HOME && ulimit -v {ULIMIT_KB}; "
             f"timeout {int(timeout)} Singular -q")
    cmd = ("wsl.exe", "-d", "Ubuntu", "--", "bash", "-lc", inner)
    try:
        cp = subprocess.run(cmd, input=program, text=True, encoding="utf-8",
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            timeout=timeout + 30, check=False)
        out = (cp.stdout or "").replace("\x00", "")
        err = (cp.stderr or "").replace("\x00", "")
        status = "ok" if cp.returncode == 0 else f"rc{cp.returncode}"
    except subprocess.TimeoutExpired:
        return {"status": "relay_timeout", "verdict": "TIMEOUT", "unit": None,
                "dim": None, "wall": round(time.monotonic() - t0, 2)}
    combined = out + "\n" + err
    um = re.search(r"@@UNIT\s*\r?\n\s*(-?\d+)", combined)
    dm = re.search(r"@@DIM\s*\r?\n\s*(-?\d+)", combined)
    unit = None if um is None else bool(int(um.group(1)))
    dim = None if dm is None else int(dm.group(1))
    if unit is None and cp.returncode == 124:      # WSL-side `timeout` fired
        verdict = "TIMEOUT"
    else:
        verdict = "UNIT" if unit else ("PROPER" if unit is False
                                       else "PARSE_FAIL")
    return {"status": status, "verdict": verdict, "unit": unit, "dim": dim,
            "wall": round(time.monotonic() - t0, 2),
            "stderr": err.strip()[:300] if verdict == "PARSE_FAIL" else ""}


bsw._run = run_singular_guarded          # runtime install, no file edited


# --------------------------------------------------------------------------
#  the cached symbolic H-system
# --------------------------------------------------------------------------
def load_H() -> dict[str, sp.Expr]:
    names = {"d2": fsb.D2, "d1": fsb.D1, "d0": fsb.D0, "dm1": fsb.DM1,
             "dm2": fsb.DM2, "dm3": fsb.DM3, "Phi": fsb.PHI}
    data = json.load(open(ELIM))
    return {k: sp.sympify(v, locals=names) for k, v in data["H"].items()}


def _spare_polys():
    rs = sp.symbols(f"R0:{CAPS['dm2'] + 1}")
    ss = sp.symbols(f"S0:{CAPS['dm3'] + 1}")
    dm2 = sum(c * y ** i for i, c in enumerate(rs))
    dm3 = sum(c * y ** i for i, c in enumerate(ss))
    return dm2, dm3, list(rs) + list(ss)


def _coeff_list(expr: sp.Expr) -> list[sp.Expr]:
    """Ascending y-coefficient list of an (already expanded) polynomial."""
    expr = sp.expand(expr)
    if expr == 0:
        return [sp.Integer(0)]
    p = sp.Poly(expr, y)
    return list(reversed(p.all_coeffs()))


def _convolve(a: list[sp.Expr], b: list[sp.Expr]) -> list[sp.Expr]:
    """Coefficient convolution of two y-polynomials; a is the small factor."""
    out = [sp.Integer(0)] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            if bj == 0:
                continue
            out[i + j] += sp.expand(ai * bj)
    return out


def _h_coeffs(cof_a, G_a, cof_b, G_b) -> list[sp.Expr]:
    """y-coefficients of cof_a*G_a + cof_b*G_b (all coefficient lists).
    The cofactors are the SMALL known polynomials (state e, spare ansatz,
    window products); the H build never forms a large symbolic product."""
    ta, tb = _convolve(cof_a, G_a), _convolve(cof_b, G_b)
    n = max(len(ta), len(tb))
    ta += [sp.Integer(0)] * (n - len(ta))
    tb += [sp.Integer(0)] * (n - len(tb))
    return [sp.expand(x + z) for x, z in zip(ta, tb)]


def _g_coeff_lists(ansatz, dm2, dm3):
    """Coefficient lists of the substituted G1,G2,G3,G5, built directly by
    per-term convolution of the ABSTRACT generators (5-12 terms each) --
    bypassing fsb.augment's monolithic sp.expand, which took ~15 min/state
    on deg_e=10 ansaetze and dominated the sweep.  Same generators
    (fsb.gsystem(), never hand-copied), same substitution values."""
    msyms = list(sp.symbols(f"M0:{CAPS['dm4'] + 1}"))
    dm4 = sum(c * y ** i for i, c in enumerate(msyms))
    vals = {fsb.D2: _coeff_list(ansatz.d2), fsb.D1: _coeff_list(ansatz.d1),
            fsb.D0: _coeff_list(ansatz.d0), fsb.DM1: _coeff_list(ansatz.e),
            fsb.DM2: _coeff_list(dm2), fsb.DM3: _coeff_list(dm3),
            fsb.DM4: _coeff_list(dm4),
            fsb.PHI: _coeff_list(fsb.phi_stripped())}
    out = {}
    for name, g in fsb.gsystem().items():
        acc = [sp.Integer(0)]
        for term in sp.Add.make_args(sp.expand(g)):
            coeff, factors = term.as_coeff_mul()
            lst = [sp.sympify(coeff)]
            for f in factors:
                base, ex = f.as_base_exp()
                if base.is_number:
                    lst = [c * base ** ex for c in lst]
                    continue
                for _ in range(int(ex)):
                    lst = _convolve(lst, vals[base])
            n = max(len(acc), len(lst))
            acc += [sp.Integer(0)] * (n - len(acc))
            lst += [sp.Integer(0)] * (n - len(lst))
            acc = [a + b for a, b in zip(acc, lst)]
        out[name] = acc
    return out, set(msyms)


def build_reduced(ansatz, *, label: str, monic_e: sp.Expr,
                  sat_factors, extra_eqs=(), use_cache=True) -> dict:
    """Instantiate the dm4-eliminated system for one state.

    Route: build the substituted G1,G2,G3,G5 coefficient lists by per-term
    convolution (_g_coeff_lists -- validated against the fsb.augment route on
    R9 z=1), then form the certified combinations H2 = dm1*G2 - dm2*G1,
    H3 = dm1*G3 - dm3*G1, H5 = dm1*G5 + (d0*dm1+d1*dm2+d2*dm3)*G1 in
    y-COEFFICIENT space with the small cofactor lists.  Every M-symbol (dm4
    coefficient) must cancel exactly -- asserted per coefficient, which
    doubles as an independent check of the symbolic elimination identity."""
    # pickle cache holds this script's OWN self-generated build output only
    # (trusted local data, same pattern as bridge_sweep.build_r9_bridge);
    # it merely skips the symbolic expansion on re-runs.
    cache_f = CACHE / f"r9red_{label}.pkl"
    if use_cache and cache_f.exists():
        return pickle.loads(cache_f.read_bytes())
    dm2, dm3, spare_unk = _spare_polys()
    G, msyms = _g_coeff_lists(ansatz, dm2, dm3)
    e_c = _coeff_list(ansatz.e)
    neg_dm2_c = [-c for c in _coeff_list(dm2)]
    neg_dm3_c = [-c for c in _coeff_list(dm3)]
    cof5_c = _coeff_list(ansatz.d0 * ansatz.e + ansatz.d1 * dm2
                         + ansatz.d2 * dm3)
    Hlists = {
        "H2": _h_coeffs(e_c, G["G2"], neg_dm2_c, G["G1"]),
        "H3": _h_coeffs(e_c, G["G3"], neg_dm3_c, G["G1"]),
        "H5": _h_coeffs(e_c, G["G5"], cof5_c, G["G1"]),
    }
    eqs: list[sp.Expr] = []
    sizes = {}
    for hname, lst in Hlists.items():
        n0 = len(eqs)
        for c in lst:
            if c == 0:
                continue
            assert not (c.free_symbols & msyms), \
                f"{hname}: dm4 coefficients failed to cancel"
            eqs.append(c)
        sizes[hname] = len(eqs) - n0
    # G1 divisibility: monic(e) | dm2*dm3  =>  rem coefficients all vanish
    prod = sp.expand(dm2 * dm3)
    _, remainder = sp.div(prod, monic_e, y)
    n0 = len(eqs)
    for c in sp.Poly(sp.expand(remainder), y).all_coeffs():
        c = sp.expand(c)
        if c != 0:
            eqs.append(c)
    sizes["divisibility"] = len(eqs) - n0
    eqs.extend(sp.expand(e) for e in extra_eqs)
    sat = list(sat_factors)
    rv = bsw._ring_vars(eqs, sat)
    d = {"label": label, "equations": eqs, "sat_factors": sat,
         "ring_vars": rv, "n_equations": len(eqs), "n_unknowns": len(rv),
         "sizes": sizes, "n_spare": len(spare_unk)}
    try:
        cache_f.write_bytes(pickle.dumps(d))
    except Exception:
        pass
    return d


# --------------------------------------------------------------------------
#  state builders
# --------------------------------------------------------------------------
def build_r9_reduced(z: int, *, use_cache=True) -> dict:
    st = qs.build_qsupport_ansatz(z)
    return build_reduced(
        st.ansatz, label=f"R9_z{z}", monic_e=(y + 1) ** 9 * (y - R),
        sat_factors=list(st.saturation_factors),
        extra_eqs=[sp.expand(qs.Q_R)], use_cache=use_cache)


def dege10_t2_states():
    d = json.load(open(ROOT / "batch_convolution_sub2.json"))
    out = []
    for s in d["states"]:
        if (s["final_verdict"] == "UNRESOLVED" and s["branch"] == "T2"
                and int(s["deg_e"]) == 10 and int(s["a_t"]) in (7, 8, 9)):
            out.append(s)

    def _d(v):
        return 0 if v in ("-inf", None) else int(v)
    # cheapest-first: small tail first (m = deg_e - a_t), then total degree
    out.sort(key=lambda s: (int(s["deg_e"]) - int(s["a_t"]),
                            _d(s["deg_sigma"]) + _d(s["deg_d1"])
                            + _d(s["deg_d2"])))
    return out


def batch_label(s):
    return (f"batch_a{s['a_t']}{s['branch']}_e{s['deg_e']}"
            f"_d2{s['deg_d2']}_d1{s['deg_d1']}_sig{s['deg_sigma']}")


def build_batch_reduced(s, *, use_cache=True) -> dict:
    ans = bsw._batch_ansatz(s)
    gamma = sp.Symbol("gamma", nonzero=True)
    return build_reduced(
        ans, label=batch_label(s), monic_e=(y + 1) ** int(s["a_t"]),
        sat_factors=[gamma], use_cache=use_cache)


# --------------------------------------------------------------------------
#  per-state protocol
# --------------------------------------------------------------------------
def _record(rec):
    prior = json.load(open(OUT)) if OUT.exists() else {"states": []}
    prior["states"] = [x for x in prior["states"]
                       if x["label"] != rec["label"]] + [rec]
    json.dump(prior, open(OUT, "w"), indent=1, default=str)


def _kill_payload(bs):
    """Integer-cleared generator strings + saturation, enough to re-emit the
    exact system later for a cofactor certificate."""
    rv = list(bs["ring_vars"]) + [sp.Symbol("w")]
    gens = []
    for g in bs["equations"]:
        s_ = fsb._to_singular(g, rv)
        if s_ not in ("0", ""):
            gens.append(s_)
    return {"ring_vars": [v.name for v in rv],
            "generators_integer_cleared": gens,
            "sat_factors": [str(f) for f in bs["sat_factors"]]}


def attack(bs, *, marked_root: bool, exact_timeout=300.0, triage_timeout=45.0):
    print(f"\n{bs['label']}: {bs['n_equations']} eqs, {bs['n_unknowns']} vars "
          f"(spare={bs['n_spare']}, dm4 eliminated); sizes={bs['sizes']}",
          flush=True)
    t0 = time.monotonic()
    if marked_root:
        tri = bsw.triage_bridge_numroot(bs, timeout=triage_timeout)
    else:
        tri = bsw.triage_bridge(bs, timeout=triage_timeout)
    print(f"  mod-p: {tri['prediction']}", flush=True)
    exact = None
    verdicts = [p.get("verdict") for p in tri["primes"]]
    # exact-Q is attempted only on positive mod-p evidence (some prime went
    # UNIT, none PROPER): empirically exact GB on the same formulation is
    # never cheaper than mod-p, so an all-TIMEOUT triage means the exact
    # attempt would burn its full budget for nothing -- record COST instead.
    if "PROPER" not in verdicts and any(v == "UNIT" for v in verdicts):
        if marked_root:
            exact = bsw.exact_bridge_minpoly(bs, timeout=exact_timeout)
            if exact.get("verdict") == "PARSE_FAIL":
                # formulation problem, not cost -- one fallback attempt
                exact2 = bsw.exact_bridge(bs, timeout=exact_timeout)
                if exact2.get("verdict") in ("UNIT", "PROPER"):
                    exact = exact2
        else:
            exact = bsw.exact_bridge(bs, timeout=exact_timeout)
    if exact and exact.get("verdict") == "UNIT":
        verdict = "KILLED-REDUCED (PENDING AUDIT)"
    elif (exact and exact.get("verdict") == "PROPER") or "PROPER" in verdicts:
        verdict = "REDUCED-PROPER (INCONCLUSIVE -- weaker than full bridge)"
    else:
        verdict = "COST"
    rec = {"label": bs["label"], "n_equations": bs["n_equations"],
           "n_unknowns": bs["n_unknowns"], "sizes": bs["sizes"],
           "modp": tri, "exact": exact, "verdict": verdict,
           "wall": round(time.monotonic() - t0, 1)}
    if verdict.startswith("KILLED"):
        rec["kill_system"] = _kill_payload(bs)
    _record(rec)
    print(f"  ==> {bs['label']}: {verdict} ({rec['wall']}s)", flush=True)
    return rec


# --------------------------------------------------------------------------
#  drivers
# --------------------------------------------------------------------------
def run_r9(zs=range(1, 7), **kw):
    recs = []
    for z in zs:
        t0 = time.monotonic()
        try:
            bs = build_r9_reduced(z)
        except Exception as ex:
            recs.append({"label": f"R9_z{z}", "verdict": "BUILD_ERROR",
                         "error": str(ex)[:200]})
            _record(recs[-1])
            print(f"  R9_z{z}: BUILD_ERROR {ex}", flush=True)
            continue
        print(f"  [build z={z}: {time.monotonic()-t0:.1f}s]", flush=True)
        recs.append(attack(bs, marked_root=True, **kw))
    return recs


def run_batch(limit=None, budget_min=None, **kw):
    states = dege10_t2_states()
    print(f"deg_e=10 T2 UNRESOLVED census: {len(states)} states", flush=True)
    t_start = time.monotonic()
    recs = []
    for s in (states if limit is None else states[:limit]):
        if budget_min and (time.monotonic() - t_start) / 60 > budget_min:
            print(f"  budget ({budget_min} min) exhausted; "
                  f"{len(states)-len(recs)} states not attempted", flush=True)
            break
        t0 = time.monotonic()
        try:
            bs = build_batch_reduced(s)
        except Exception as ex:
            recs.append({"label": batch_label(s), "verdict": "BUILD_ERROR",
                         "error": str(ex)[:200]})
            _record(recs[-1])
            continue
        print(f"  [build {bs['label']}: {time.monotonic()-t0:.1f}s]",
              flush=True)
        recs.append(attack(bs, marked_root=False, **kw))
    return recs


def census():
    if not OUT.exists():
        print("no sweep record yet")
        return
    d = json.load(open(OUT))
    from collections import Counter
    c = Counter(x["verdict"].split(" ")[0] for x in d["states"])
    print(dict(c))
    for x in d["states"]:
        print(f"  {x['label']}: {x['verdict']} ({x.get('wall','?')}s)")


if __name__ == "__main__":
    args = sys.argv[1:]
    if args and args[0] == "pilot":
        z = int(args[1]) if len(args) > 1 else 1
        bs = build_r9_reduced(z)
        attack(bs, marked_root=True)
    elif args and args[0] == "r9":
        zs = [int(x) for x in args[1:]] if len(args) > 1 else range(1, 7)
        run_r9(zs)
    elif args and args[0] == "batch":
        lim = int(args[1]) if len(args) > 1 else None
        bud = float(args[2]) if len(args) > 2 else None
        run_batch(lim, bud)
    elif args and args[0] == "census":
        census()
    elif args and args[0] == "all":
        bud = float(args[1]) if len(args) > 1 else 60.0
        run_r9()
        run_batch(budget_min=bud)
    else:
        run_r9()
        run_batch(budget_min=60)
