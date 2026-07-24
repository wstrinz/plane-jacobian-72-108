#!/usr/bin/env python3
"""alt_elim.py -- the NAMED CURE for the alt-bridge wall (ALT_BRIDGE.md sec.4a):
symbolically ELIMINATE the spare window unknowns FIRST, before grounding the
state, exactly the pattern r9_symbolic_elim.py / R9_SYMBOLIC.md certified for
the standard window.

WHY (ALT_BRIDGE.md sec.4)
-------------------------
The alt-window bridge is a SOUND construction (the sub1-caps G-system, 66 spare
scalars, marked-root adjunction) but BOTH pilots COST on every engine -- the
control (a known kill) was NOT computationally reproduced.  The cost wall is
engine-independent, so the cure is a STRUCTURALLY REDUCED system, not a bigger
budget.

THE ELIMINATION (regime-independent; certified in r9_symbolic_elim.py)
----------------------------------------------------------------------
The four pre-resultant generators G1,G2,G3,G5(=G5body+Phi) are LOADED once from
generators.json and are t-regime-INDEPENDENT: the alt (sub1) bridge and the
standard (sub2) bridge share the SAME generators; only the spare-poly degree
CAPS differ (sub1 18/21/24 vs sub2 12/14/16).  Hence the certified dm4
elimination transfers verbatim.  Greedy pivot census (alt_elim probe, this
file):

    spare  linear-in generator      pivot                 guaranteed nonzero?
    dm4    G1                        3*dm1 = 3*e           YES (e is known !=0)
    dm3    G1 / G3                   3*dm2 / 3*dm4         NO (spare-valued)
    dm2    G1                        3*d2*dm1 + 3*dm3      NO (spare-valued)

dm2,dm3 also appear only QUADRATICALLY in G2,G5 (constant pivots, but degree 2).
=> dm4 is the UNIQUE spare admitting a guaranteed-nonzero-pivot LINEAR
elimination.  Using dm1 = e (guaranteed nonzero) as the multiplier,

    H2 := dm1*G2 - dm2*G1            (dm4-free, in <G1,G2,G3,G5>)
    H3 := dm1*G3 - dm3*G1            (dm4-free)
    H5 := dm1*G5 + (d0*dm1 + d1*dm2 + d2*dm3)*G1     (dm4-free)

are dm4-FREE ideal members (cofactor certificate re-verified by exact
expansion).  Plus the G1 divisibility lemma  monic(e) | dm2*dm3  contributes
rem(dm2*dm3, monic(e)) == 0 coefficient-wise, recovering most of G1 without dm4.

SUB1 CENSUS: of the 66 spare scalars (dm2:19, dm3:22, dm4:25), the 25 dm4
scalars are eliminated SYMBOLICALLY at zero soundness cost (66 -> 41).  A
per-state scalar linear presolve (fsb.linear_presolve, constant pivots only)
mops up further.  VERDICT SEMANTICS (asymmetric, honest, from R9_SYMBOLIC):
  UNIT over Q       -> state KILLED (candidate, PENDING AUDIT)
  PROPER anywhere   -> INCONCLUSIVE (reduced is weaker than full bridge; drops
                       G1's "dm4 is a capped polynomial" -- never survival)
  TIMEOUT           -> COST (pure Groebner cost)

Orphan discipline: WSL-side `timeout` + `ulimit -v 8G` inside WSL, so a
Windows-relay death can never orphan a WSL Singular/msolve process.

New file, uncommitted.  READ-ONLY on every imported module/artifact; does NOT
touch alt_bridge.py, r9_* files, kill_certificate*, proof_dag*,
g_system_75_125*, alok_*.  Reuses alt_bridge.build_state (the audited alt-state
reconstruction) and r9_symbolic_elim.eliminate (the certified elimination).
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
import time
from pathlib import Path

import sympy as sp

import full_system_bridge as fsb
import modular_triage as mt
import alt_bridge as ab
import r9_symbolic_elim as r9e

ROOT = Path(__file__).resolve().parent
ELIM_OUT = ROOT / "alt_eliminated_system.json"
RES_OUT = ROOT / "alt_elim_results.json"

y = fsb.y
QY = ab.QY
PRIMES = (10007, 10009, 100019)
ULIMIT_KB = 8_000_000          # 8 GB virtual-memory cap per WSL process
CAPS1 = fsb.STRIP_DEGCAP["sub1"]   # dm2:18, dm3:21, dm4:24  (66 spare scalars)


# --------------------------------------------------------------------------
#  STEP 1: the symbolic elimination (regime-independent) + sub1 census
# --------------------------------------------------------------------------
def eliminate_and_cache() -> dict:
    """Re-derive + exact-expansion-verify the dm4 elimination (reusing the
    certified r9_symbolic_elim.eliminate) and cache it, framed for the alt
    (sub1) bridge with the sub1 spare-scalar census."""
    res = r9e.eliminate()          # runs every exact-expansion certificate check
    H = res.pop("_H_sympy")
    # independent greedy pivot census over the loaded generators
    g = fsb.gsystem()
    pivot_census = {}
    for spare in (fsb.DM4, fsb.DM3, fsb.DM2):
        rows = {}
        for n, e in g.items():
            p = sp.Poly(e, spare)
            deg = p.degree()
            piv = sp.sstr(sp.expand(p.nth(deg))) if deg > 0 else None
            rows[n] = {"deg": deg, "pivot": piv}
        pivot_census[str(spare)] = rows
    scalars = {"dm2": CAPS1["dm2"] + 1, "dm3": CAPS1["dm3"] + 1,
               "dm4": CAPS1["dm4"] + 1}
    census = {
        "regime": "sub1",
        "caps": CAPS1,
        "spare_scalars_total": sum(scalars.values()),      # 66
        "spare_scalars_per": scalars,                      # 19 / 22 / 25
        "eliminated_symbolically": {"dm4": scalars["dm4"]},  # 25
        "remaining_after_dm4": scalars["dm2"] + scalars["dm3"],  # 41
        "unique_guaranteed_pivot_spare": "dm4",
        "pivot_of_G1_in_dm4": "3*dm1 (= 3*e, guaranteed nonzero)",
        "greedy_note": ("dm4 is the UNIQUE spare with a guaranteed-nonzero "
                        "linear pivot; dm3/dm2 pivots are spare-valued or the "
                        "term is quadratic (see pivot_census) -- no further "
                        "guaranteed-nonzero-pivot linear elimination exists. "
                        "Per-state scalar presolve (constant pivots) removes "
                        "more at instantiation."),
    }
    out = {
        "schema": "alt-eliminated-system-v1",
        "shares_elimination_with": "r9_eliminated_system.json",
        "H": res["H"], "cofactors": res["cofactors"], "pivots": res["pivots"],
        "weights": res["weights"], "checks": res["checks"],
        "pivot_census": pivot_census,
        "sub1_census": census,
        "note": res["note"],
    }
    ELIM_OUT.write_text(json.dumps(out, indent=1), encoding="utf-8")
    return {"H": H, "meta": out}


# --------------------------------------------------------------------------
#  H-system loader (the certified dm4-FREE combinations, symbol-mapped)
# --------------------------------------------------------------------------
def _load_H():
    names = {"d2": fsb.D2, "d1": fsb.D1, "d0": fsb.D0, "dm1": fsb.DM1,
             "dm2": fsb.DM2, "dm3": fsb.DM3, "Phi": fsb.PHI}
    data = json.load(open(ELIM_OUT))
    return {k: sp.sympify(v, locals=names) for k, v in data["H"].items()}


# --------------------------------------------------------------------------
#  instantiation via coefficient-list convolution (the r9_symbolic_sweep
#  route, sub1 caps): build the substituted G1,G2,G3,G5 coefficient lists
#  ONCE with dm4 := 0 (valid -- the dm4 terms cancel identically in every H,
#  so H = dm1*G_i|_{dm4=0} - cof*G1|_{dm4=0}; no M-symbols, no blow-up), then
#  form the certified combinations H2,H3,H5 by a SINGLE convolution with the
#  small cofactor lists.  Generators LOADED from fsb.gsystem() (never copied).
# --------------------------------------------------------------------------
def _coeff_list(expr):
    expr = sp.expand(expr)
    if expr == 0:
        return [sp.Integer(0)]
    return list(reversed(sp.Poly(expr, y).all_coeffs()))


def _convolve(a, b, red=None):
    out = [sp.Integer(0)] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            if bj == 0:
                continue
            out[i + j] += sp.expand(ai * bj)
    # Reduce mod q(r_i) after each convolution: keeps the r-degree of every
    # coefficient <= 3, bounding the symbolic blow-up of e^3,e^4 over the two
    # marked roots.  SOUND: q(r_i)=0 is adjoined to the ideal, so reducing any
    # coefficient mod q(r_i) does not change the variety.
    if red is not None:
        out = [red(c) if c != 0 else c for c in out]
    return out


def _h_coeffs(cof_a, G_a, cof_b, G_b, red=None):
    ta, tb = _convolve(cof_a, G_a, red), _convolve(cof_b, G_b, red)
    n = max(len(ta), len(tb))
    ta += [sp.Integer(0)] * (n - len(ta))
    tb += [sp.Integer(0)] * (n - len(tb))
    out = [sp.expand(x + z) for x, z in zip(ta, tb)]
    if red is not None:
        out = [red(c) if c != 0 else c for c in out]
    return out


def _g_coeff_lists_sub1(d2p, d1p, d0p, ep, dm2, dm3, red=None):
    """Substituted G1,G2,G3,G5 coefficient lists with dm4 := 0."""
    vals = {fsb.D2: _coeff_list(d2p), fsb.D1: _coeff_list(d1p),
            fsb.D0: _coeff_list(d0p), fsb.DM1: _coeff_list(ep),
            fsb.DM2: _coeff_list(dm2), fsb.DM3: _coeff_list(dm3),
            fsb.DM4: [sp.Integer(0)],
            fsb.PHI: _coeff_list(fsb.phi_stripped())}
    out = {}
    for name, g in fsb.gsystem().items():
        acc = [sp.Integer(0)]
        for term in sp.Add.make_args(sp.expand(g)):
            coeff, factors = term.as_coeff_mul()
            lst = [sp.sympify(coeff)]
            for f in factors:
                b_, ex = f.as_base_exp()
                if b_.is_number:
                    lst = [c * b_ ** ex for c in lst]
                    continue
                for _ in range(int(ex)):
                    lst = _convolve(lst, vals[b_], red)
            n = max(len(acc), len(lst))
            acc += [sp.Integer(0)] * (n - len(acc))
            lst += [sp.Integer(0)] * (n - len(lst))
            acc = [a + b for a, b in zip(acc, lst)]
        out[name] = [red(c) if red and c != 0 else c for c in acc]
    return out


def _spare_polys_sub1():
    rs = sp.symbols(f"R0:{CAPS1['dm2'] + 1}")   # 19 scalars
    ss = sp.symbols(f"S0:{CAPS1['dm3'] + 1}")   # 22 scalars
    dm2 = sum(c * y ** i for i, c in enumerate(rs))
    dm3 = sum(c * y ** i for i, c in enumerate(ss))
    return dm2, dm3, list(rs) + list(ss)


def _monic_e(bid, red):
    """The monic (E/S-stripped) e polynomial for the divisibility lemma,
    reduced mod q(r_i)."""
    if bid == "a11_b3100_T2":
        return red(sp.expand((y + 1) ** 11 * (y - ab.r1) ** 3 * (y - ab.r2)))
    if bid == "a12_b1110_T2":
        comp = sp.div(sp.Poly(QY, y), sp.Poly(y - ab.r, y))[0].as_expr() / 2048
        return red(sp.expand((y + 1) ** 12 * comp))
    raise ValueError(bid)


import os
import pickle
# The pickle cache holds ONLY this script's own self-generated build output
# (trusted local data, same pattern as r9_symbolic_sweep.build_reduced); it
# merely skips the slow symbolic build on re-runs.  Never loads external data.
_CACHE = Path(os.environ.get("ALT_ELIM_CACHE", ROOT / "__pycache__"))


def build_reduced(bid: str, deg_d2: int, *, use_cache=True) -> dict:
    """Instantiate the dm4-eliminated (sub1) reduced system for an alt state."""
    cache_f = _CACHE / f"alt_elim_{bid}_dd{deg_d2}.pkl"
    if use_cache and cache_f.exists():
        return pickle.loads(cache_f.read_bytes())
    d2p, sigp, ep, root_vars, Dc = ab.build_state(bid, deg_d2)
    red = ab.reducer(root_vars)
    d1p = sp.Integer(0)
    d0p = red(sp.expand((d2p ** 2 + sigp) / 4))
    dm2, dm3, spare_unk = _spare_polys_sub1()
    rred = red if root_vars else None
    G = _g_coeff_lists_sub1(d2p, d1p, d0p, ep, dm2, dm3, rred)
    e_c = _coeff_list(ep)
    neg_dm2_c = [-c for c in _coeff_list(dm2)]
    neg_dm3_c = [-c for c in _coeff_list(dm3)]
    cof5_c = _coeff_list(red(d0p * ep) + d1p * dm2 + d2p * dm3)
    Hlists = {
        "H2": _h_coeffs(e_c, G["G2"], neg_dm2_c, G["G1"], rred),
        "H3": _h_coeffs(e_c, G["G3"], neg_dm3_c, G["G1"], rred),
        "H5": _h_coeffs(e_c, G["G5"], cof5_c, G["G1"], rred),
    }
    eqs, sizes = [], {}
    for hname, lst in Hlists.items():
        n0 = len(eqs)
        for c in lst:
            if c != 0:
                eqs.append(c)
        sizes[hname] = len(eqs) - n0
    # G1 divisibility lemma: monic(e) | dm2*dm3  =>  rem coefficients vanish
    monic_e = _monic_e(bid, red)
    _, remainder = sp.div(sp.expand(dm2 * dm3), monic_e, y)
    n0 = len(eqs)
    for c in sp.Poly(sp.expand(remainder), y).all_coeffs():
        c = red(sp.expand(c))
        if c != 0:
            eqs.append(c)
    sizes["divisibility"] = len(eqs) - n0
    # marked-root minimal polynomials
    n0 = len(eqs)
    for rv in root_vars:
        eqs.append(ab.qpoly(rv))
    sizes["q_roots"] = len(eqs) - n0
    # dedup
    uniq, seen = [], set()
    for e in eqs:
        k = sp.sstr(e)
        if k not in seen:
            seen.add(k)
            uniq.append(e)
    eqs = uniq
    sat = [ab.E, ab.S, Dc[-1]]
    if len(root_vars) == 2:
        sat.append(root_vars[0] - root_vars[1])
    ring_vars = mt.ring_vars_of(eqs, extra=[s for f in sat
                                            for s in sp.sympify(f).free_symbols])
    d = {
        "label": f"{bid}_degd2_{deg_d2}", "bid": bid, "deg_d2": deg_d2,
        "regime": "sub1", "equations": eqs, "sat_factors": sat,
        "root_vars": root_vars, "ring_vars": ring_vars,
        "n_equations": len(eqs), "n_unknowns": len(ring_vars),
        "n_spare": len(spare_unk), "sizes": sizes,
        "state_scalars": [v.name for v in list(Dc)],
    }
    try:
        cache_f.write_bytes(pickle.dumps(d))
    except Exception:
        pass
    return d


def presolve(bs: dict) -> dict:
    """Scalar linear presolve (constant pivots only) over the spare + d2
    coefficients -- sound equivalent rewriting; further shrinks the system."""
    unk = [v for v in bs["ring_vars"]
           if re.fullmatch(r"[RS]\d+|D\d+", v.name)]
    eqs2, remaining, submap = fsb.linear_presolve(bs["equations"], unk)
    # sat factors may reference eliminated D-coeffs (Dc[-1]); rewrite them
    sat2 = [sp.expand(sp.sympify(f).subs(submap)) for f in bs["sat_factors"]]
    sat2 = [f for f in sat2 if sp.sympify(f).free_symbols]
    rv2 = mt.ring_vars_of(eqs2, extra=[s for f in sat2
                                       for s in sp.sympify(f).free_symbols])
    out = dict(bs)
    out.update({"equations": eqs2, "sat_factors": sat2, "ring_vars": rv2,
                "n_equations": len(eqs2), "n_unknowns": len(rv2),
                "eliminated_scalar": len(submap)})
    return out


# --------------------------------------------------------------------------
#  NATIVE Singular emission of the reduced (dm4-eliminated) H-system.
#  The cubic-in-41-spare-symbol expansion is intractable in sympy, so it is
#  pushed entirely into Singular (C-speed): the abstract generators G1..G5
#  (loaded from fsb.gsystem(), dm4:=0) and the small state polys are DEFINED as
#  Singular polynomials, the certified combinations H2,H3,H5 are formed, and
#  their y-coefficients are extracted with coeffs(H_i, y) as ideal generators.
#  This is the exact reduced system -- the only difference from the sympy path
#  is WHERE the polynomial arithmetic happens.  The H-coefficients alone are a
#  SOUND necessary system (UNIT => state killed); the G1 divisibility lemma is
#  an optional strengthening, omitted here because forming it needs a y-only
#  division that is itself costly -- soundness of the kill test is unaffected.
# --------------------------------------------------------------------------
def _sing(expr) -> str:
    """Singular string with the rational coefficient (and its '/') written
    FIRST in every term, so '^' is never followed by '/'.  This dodges the
    Singular `y^N/M` parser trap (memory: 'gm^8/N parser trap' -- Singular
    misreads y^2/128 as y^(2/128))."""
    expr = sp.expand(expr)
    if expr == 0:
        return "0"
    parts = []
    for t in sp.Add.make_args(expr):
        c, rest = t.as_coeff_Mul()
        c = sp.Rational(c)
        cs = f"{c.p}" if c.q == 1 else f"{c.p}/{c.q}"
        if rest == 1:
            parts.append(f"({cs})")
        else:
            rs = sp.sstr(rest).replace("**", "^").replace(" ", "")
            parts.append(f"{cs}*{rs}")
    return "+".join(parts).replace("+-", "-")


_GSTR = None


def _gen_strings():
    """Abstract G1..G5 as Singular strings in the window-var poly names
    (d2,d1,d0,dm1,dm2,dm3,dm4,Phi); loaded once, never hand-copied."""
    global _GSTR
    if _GSTR is None:
        _GSTR = {n: _sing(g) for n, g in fsb.gsystem().items()}
    return _GSTR


def build_state_native(bid: str, deg_d2: int) -> dict:
    """Cheap: only the small state polys (e, sigma, d2) + metadata; the heavy
    H expansion is deferred to Singular."""
    d2p, sigp, ep, root_vars, Dc = ab.build_state(bid, deg_d2)
    phi = fsb.phi_stripped()
    return {"label": f"{bid}_degd2_{deg_d2}", "bid": bid, "deg_d2": deg_d2,
            "regime": "sub1", "d2p": d2p, "sigp": sigp, "ep": ep, "phi": phi,
            "root_vars": root_vars, "Dc": list(Dc), "n_spare": 41}


def emit_native(st: dict, *, char: int, roots=None) -> str:
    """Singular program for the reduced H-system.  roots: optional list of
    numeric values for the marked roots (numroot reconnaissance mod p); when
    given, the marked roots + q(r) + (r1-r2) saturation are specialized away."""
    root_vars = list(st["root_vars"])
    numeric = roots is not None
    subst = dict(zip(root_vars, roots)) if numeric else {}
    ep = sp.expand(st["ep"].subs(subst)) if numeric else st["ep"]
    sigp = sp.expand(st["sigp"].subs(subst)) if numeric else st["sigp"]
    d2p = st["d2p"]
    Dc = st["Dc"]
    gen = _gen_strings()
    spare_r = [f"R{i}" for i in range(CAPS1["dm2"] + 1)]   # dm2 coeffs
    spare_s = [f"S{i}" for i in range(CAPS1["dm3"] + 1)]   # dm3 coeffs
    dcoef = [v.name for v in Dc]
    ring_vars = (["y"] + ([v.name for v in root_vars] if not numeric else [])
                 + spare_r + spare_s + dcoef + ["E", "Sg", "w"])
    var_txt = ",".join(ring_vars)
    L = [f"ring ELR = {char},({var_txt}),dp;"]
    # state polys (Sg is the sigma scale symbol 'S'; renamed to avoid clash)
    def _rs(expr):
        return _sing(sp.sympify(expr).subs(ab.S, sp.Symbol("Sg")))
    L.append(f"poly d2 = {_sing(d2p)};")
    L.append(f"poly dm1 = {_rs(ep)};")
    L.append(f"poly sig = {_rs(sigp)};")
    L.append("poly d0 = (d2^2+sig)/4;")
    L.append("poly d1 = 0;")
    L.append("poly dm4 = 0;")
    L.append("poly dm2 = " + "+".join(f"{c}*y^{i}" for i, c in enumerate(spare_r)) + ";")
    L.append("poly dm3 = " + "+".join(f"{c}*y^{i}" for i, c in enumerate(spare_s)) + ";")
    L.append(f"poly Phi = {_sing(st['phi'])};")
    for n in ("G1", "G2", "G3", "G5"):
        L.append(f"poly {n} = {gen[n]};")
    L.append("poly H2 = dm1*G2 - dm2*G1;")
    L.append("poly H3 = dm1*G3 - dm3*G1;")
    L.append("poly H5 = dm1*G5 + (d0*dm1 + d1*dm2 + d2*dm3)*G1;")
    L.append("ideal I = ideal(coeffs(H2,y)) + ideal(coeffs(H3,y)) "
             "+ ideal(coeffs(H5,y));")
    if not numeric:
        for rv in root_vars:
            L.append(f"I = I + ideal({_sing(ab.qpoly(rv))});")
    # saturation: E, Sg(=sigma scale), leading d2 coeff, (r1-r2) if two roots
    sat = ["E", "Sg", dcoef[-1]]
    if not numeric and len(root_vars) == 2:
        sat.append(f"({root_vars[0].name}-{root_vars[1].name})")
    L.append("poly nz = w*" + "*".join(sat) + "-1;")
    L.append("I = I + ideal(nz);")
    L.append('"@@NEQ"; size(I);')
    L.append("ideal G = std(I);")
    L.append("int u = (reduce(1,G)==0);")
    L.append('"@@UNIT"; u; "@@DIM"; dim(G);')
    L.append("quit;")
    return "\n".join(L) + "\n"


# --------------------------------------------------------------------------
#  orphan-proof runners (WSL-side timeout + ulimit -v 8G)
# --------------------------------------------------------------------------
def run_singular_guarded(program: str, timeout: float) -> dict:
    t0 = time.monotonic()
    inner = (f"cd $HOME && ulimit -v {ULIMIT_KB}; "
             f"timeout {int(max(5, timeout))}s Singular -q")
    cmd = ("wsl.exe", "-d", "Ubuntu", "--", "bash", "-lc", inner)
    try:
        cp = subprocess.run(cmd, input=program, text=True, encoding="utf-8",
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            timeout=timeout + 30, check=False)
    except subprocess.TimeoutExpired:
        return {"status": "relay_timeout", "verdict": "TIMEOUT", "unit": None,
                "dim": None, "wall": round(time.monotonic() - t0, 1)}
    wall = round(time.monotonic() - t0, 1)
    combined = ((cp.stdout or "") + "\n" + (cp.stderr or "")).replace("\x00", "")
    um = re.search(r"@@UNIT\s*\r?\n\s*(-?\d+)", combined)
    dm = re.search(r"@@DIM\s*\r?\n\s*(-?\d+)", combined)
    nq = re.search(r"@@NEQ\s*\r?\n\s*(-?\d+)", combined)
    neq = int(nq.group(1)) if nq else None
    if um is None:
        return {"status": "timeout", "verdict": "TIMEOUT", "unit": None,
                "dim": None, "neq": neq, "wall": wall}
    unit = bool(int(um.group(1)))
    return {"status": "ok", "verdict": "UNIT" if unit else "PROPER",
            "unit": unit, "dim": None if dm is None else int(dm.group(1)),
            "neq": neq, "wall": wall}


_MS = "$HOME/msolve/msolve"


def exact_msolve_guarded(bs: dict, *, timeout=300.0, char=0) -> dict:
    w = sp.Symbol("w")
    rv = list(bs["ring_vars"]) + [w]
    prod = sp.Integer(1)
    for f in bs["sat_factors"]:
        prod = prod * f
    gens = list(bs["equations"]) + [sp.expand(w * prod - 1)]
    polys = [s for s in (fsb._to_singular(g, rv) for g in gens)
             if s not in ("0", "")]
    prog = (",".join(v.name for v in rv) + "\n" + str(char) + "\n"
            + ",\n".join(polys) + "\n")
    tag = bs["label"].replace(" ", "")
    fname, outname = f"altelim_{tag}.ms", f"altelim_{tag}.out"
    subprocess.run(("wsl.exe", "-d", "Ubuntu", "--", "bash", "-lc",
                    f"cat > $HOME/{fname}"), input=prog, text=True,
                   encoding="utf-8", check=True)
    tsec = int(max(5, timeout))
    run = ("wsl.exe", "-d", "Ubuntu", "--", "bash", "-lc",
           f"cd $HOME && ulimit -v {ULIMIT_KB}; "
           f"timeout {tsec}s {_MS} -f $HOME/{fname} -o $HOME/{outname}; "
           f"cat $HOME/{outname} 2>/dev/null")
    t0 = time.monotonic()
    try:
        cp = subprocess.run(run, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            text=True, encoding="utf-8", timeout=timeout + 30,
                            check=False)
    except subprocess.TimeoutExpired:
        return {"verdict": "TIMEOUT", "wall": round(time.monotonic() - t0, 1),
                "engine": "msolve"}
    body = (cp.stdout or "").strip()
    wall = round(time.monotonic() - t0, 1)
    if body.startswith("[-1]"):
        verdict = "UNIT"
    elif re.match(r"\[1,\s*\d+,\s*-1", body):
        verdict = "INFINITE_SOL"
    elif body:
        verdict = "HAS_SOL"
    else:
        verdict = "TIMEOUT"
    return {"verdict": verdict, "wall": wall, "engine": "msolve",
            "out_head": body[:120]}


def triage_numroot(bs, *, primes=PRIMES, timeout=45.0):
    """mod-p reconnaissance: marked roots specialized to numeric q-roots mod p
    (strips r-nonlinearity).  Evidence over F_p-bar, not a Q certificate."""
    root_syms = list(bs["root_vars"])
    nroots = len(root_syms)
    qeqs = {sp.sstr(ab.qpoly(rv)) for rv in root_syms}
    eqs0 = [e for e in bs["equations"] if sp.sstr(e) not in qeqs]
    out = []
    for p in primes:
        roots = mt.q_roots_mod_p(p, nroots)
        if len(roots) < nroots:
            out.append({"verdict": "SKIP", "prime": p})
            print(f"    (numroot) p={p}: SKIP (q<{nroots} roots)", flush=True)
            continue
        subst = dict(zip(root_syms, roots[:nroots]))
        eqs = [sp.expand(sp.sympify(e).subs(subst)) for e in eqs0]
        sat = [sp.expand(sp.sympify(f).subs(subst)) for f in bs["sat_factors"]]
        sat = [f for f in sat if sp.sympify(f).free_symbols]
        rv = mt.ring_vars_of(eqs, extra=[s for f in sat
                                         for s in sp.sympify(f).free_symbols])
        rr = run_singular_guarded(ab.emit(eqs, sat, rv, char=p), timeout)
        rr["prime"] = p
        rr["roots"] = [int(x) for x in roots[:nroots]]
        out.append(rr)
        print(f"    (numroot) p={p}: {rr['verdict']} dim={rr.get('dim')} "
              f"({rr['wall']}s)", flush=True)
    return out


def exact_singular(bs, *, timeout=300.0):
    rr = run_singular_guarded(
        ab.emit(bs["equations"], bs["sat_factors"], bs["ring_vars"], char=0),
        timeout)
    print(f"    exact-Q Singular: {rr['verdict']} dim={rr.get('dim')} "
          f"({rr['wall']}s)", flush=True)
    return rr


def _kill_payload(bs):
    rv = list(bs["ring_vars"]) + [sp.Symbol("w")]
    gens = [s for s in (fsb._to_singular(g, rv) for g in bs["equations"])
            if s not in ("0", "")]
    return {"ring_vars": [v.name for v in rv],
            "generators_integer_cleared": gens,
            "sat_factors": [str(f) for f in bs["sat_factors"]]}


# --------------------------------------------------------------------------
#  per-state protocol: numroot mod-p (45s x3) -> msolve char-0 (300s)
#                      -> Singular exact-Q (300s)
# --------------------------------------------------------------------------
def _kill_payload_native(st):
    return {"program_char0": emit_native(st, char=0),
            "note": "re-run in Singular; @@UNIT 1 == empty variety == kill"}


def triage_numroot_native(st, *, primes=PRIMES, timeout=45.0):
    nroots = len(st["root_vars"])
    out = []
    for p in primes:
        roots = mt.q_roots_mod_p(p, nroots)
        if len(roots) < nroots:
            out.append({"verdict": "SKIP", "prime": p})
            print(f"    (numroot) p={p}: SKIP (q<{nroots} roots)", flush=True)
            continue
        prog = emit_native(st, char=p, roots=list(roots[:nroots]))
        rr = run_singular_guarded(prog, timeout)
        rr["prime"] = p
        rr["roots"] = [int(x) for x in roots[:nroots]]
        out.append(rr)
        print(f"    (numroot) p={p}: {rr['verdict']} dim={rr.get('dim')} "
              f"neq={rr.get('neq')} ({rr['wall']}s)", flush=True)
    return out


def attack_native(st, *, triage_timeout=45.0, exact_timeout=300.0):
    print(f"\n{st['label']}: reduced dm4-eliminated H-system (sub1, "
          f"{st['n_spare']} spare after elim), native Singular coeffs route",
          flush=True)
    t0 = time.monotonic()
    rec = {k: st[k] for k in ("label", "bid", "deg_d2", "regime", "n_spare")}
    print("  numroot mod-p triage (45s x3):", flush=True)
    rec["triage_numroot"] = triage_numroot_native(st, timeout=triage_timeout)
    nrv = {t["verdict"] for t in rec["triage_numroot"] if t["verdict"] != "SKIP"}
    rec["n_equations"] = next((t.get("neq") for t in rec["triage_numroot"]
                               if t.get("neq")), None)
    # exact char-0 verdict via Singular-native std over Q (msolve needs
    # pre-expanded polys -- exactly the intractable sympy step -- so the exact
    # engine here is Singular's own F4/std over Q on the coeffs()-built ideal).
    print("  exact char-0 Singular-native std over Q (300s):", flush=True)
    rec["exact_singular"] = run_singular_guarded(
        emit_native(st, char=0), exact_timeout)
    sing = rec["exact_singular"]
    print(f"    exact-Q: {sing['verdict']} dim={sing.get('dim')} "
          f"({sing['wall']}s)", flush=True)
    sv = sing["verdict"]
    if sv == "UNIT" or nrv == {"UNIT"}:
        if sv == "UNIT":
            rec["verdict"] = "KILLED-REDUCED (PENDING AUDIT)"
            rec["kill_system"] = _kill_payload_native(st)
        else:
            rec["verdict"] = "MODP-UNIT (numroot; char-0 PENDING)"
    elif "PROPER" in nrv or sv == "PROPER":
        rec["verdict"] = "REDUCED-PROPER (INCONCLUSIVE -- weaker than full bridge)"
    else:
        rec["verdict"] = "COST"
    rec["wall"] = round(time.monotonic() - t0, 1)
    print(f"  ==> {st['label']}: {rec['verdict']} ({rec['wall']}s)", flush=True)
    return rec


# --------------------------------------------------------------------------
#  drivers
# --------------------------------------------------------------------------
def _save(results):
    RES_OUT.write_text(json.dumps(results, indent=1, default=str),
                       encoding="utf-8")


def run_pilots(results):
    for bid, dd in (("a12_b1110_T2", 6), ("a11_b3100_T2", 6)):
        print(f"\n[state {bid} deg_d2={dd}] building state polys ...", flush=True)
        st = build_state_native(bid, dd)
        results["cases"].append(attack_native(st))
        results["role"][st["label"]] = ("control" if bid.startswith("a12")
                                        else "prize")
        _save(results)
    return results


def sweep_frontier(results, *, budget_min=100):
    """If the control KILLED, sweep additional alt frontier states within a
    wall-clock budget."""
    t_start = time.monotonic()
    for bid in ("a12_b1110_T2", "a11_b3100_T2"):
        for dd in (5, 4, 3, 2, 1, 0):
            if (time.monotonic() - t_start) / 60 > budget_min:
                print(f"  sweep budget ({budget_min} min) exhausted", flush=True)
                _save(results)
                return results
            lbl = f"{bid}_degd2_{dd}"
            if any(c["label"] == lbl for c in results["cases"]):
                continue
            print(f"\n[sweep {lbl}] ...", flush=True)
            try:
                st = build_state_native(bid, dd)
            except Exception as ex:
                results["cases"].append({"label": lbl, "verdict": "BUILD_ERROR",
                                         "error": str(ex)[:200]})
                _save(results)
                continue
            results["cases"].append(attack_native(st))
            _save(results)
    return results


def main(argv):
    mode = argv[1] if len(argv) > 1 else "all"
    print("== STEP 1: symbolic elimination + sub1 census ==", flush=True)
    elim = eliminate_and_cache()
    for c in elim["meta"]["checks"]:
        print(f"  [OK] {c}", flush=True)
    cen = elim["meta"]["sub1_census"]
    print(f"  sub1: {cen['spare_scalars_total']} spare scalars "
          f"({cen['spare_scalars_per']}); eliminated dm4="
          f"{cen['eliminated_symbolically']['dm4']} -> "
          f"{cen['remaining_after_dm4']} remain; unique-pivot spare="
          f"{cen['unique_guaranteed_pivot_spare']}", flush=True)
    print(f"  wrote {ELIM_OUT.name}", flush=True)
    if mode == "elim":
        return
    results = {"schema": "alt-elim-v1", "role": {},
               "elimination_census": cen, "cases": []}
    _save(results)
    print("\n== STEP 2: re-run the two pilots on the REDUCED system ==",
          flush=True)
    run_pilots(results)
    # STEP 3: sweep frontier IFF the control reproduced its known kill
    control = next((c for c in results["cases"]
                    if c.get("bid") == "a12_b1110_T2"), None)
    ctrl_kill = control and control["verdict"].startswith("KILLED")
    results["control_verdict"] = control["verdict"] if control else None
    if mode in ("all", "sweep") and ctrl_kill:
        print("\n== STEP 3: CONTROL KILLED -> sweeping alt frontier ==",
              flush=True)
        sweep_frontier(results, budget_min=110)
    elif mode in ("all", "sweep"):
        print("\n== STEP 3: control did NOT reproduce the kill -> "
              "NO sweep (wall is deeper than spare count) ==", flush=True)
    print("\n== census ==", flush=True)
    for c in results["cases"]:
        print(f"  {c['label']}: {c['verdict']}", flush=True)
    _save(results)


if __name__ == "__main__":
    main(sys.argv)
