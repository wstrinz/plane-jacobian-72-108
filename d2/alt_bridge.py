#!/usr/bin/env python3
"""ALT_BRIDGE -- the full-system bridge analog for the ALTERNATE regime.

CONSTRUCTION (the result; full derivation in ALT_BRIDGE.md):
the analog IS the landed G-system with the SUB1 window caps.  The alternate
regime (ALT_REGIME.md: subcase (1), a = v_t(e) in 11..15, v = 30-3a < 0) is a
t-adic sub-locus of the SUB1 window; every ingredient of the bridge is
t-regime-independent:

  1. The generators G1,G2,G3,G5body+Phi are (D~^3)_{-1,-2,-3,-5} identities --
     x-level consequences of C^2 = P and the d3-killing shift.  No step of
     their derivation (T6_SELECTION_AUDIT.md / regenerate_system.py) consumes
     the t-adic profile of e or the standard reduction F = t^(21a)G.
  2. f31 lies in the G-ideal (exact cofactor certificate,
     full_system_bridge_verify.py V2), and ALT_REGIME.md's own survival table
     row 1 states f31 "survives verbatim; it is window-independent."
  3. The stripped Phi is the SAME object in both regimes: ALT_REGIME.md writes
     Phi~ = t^30*u with u = c*q -- literally full_system_bridge.phi_stripped().
  4. The window caps (ord >= 12k, deg <= 15k sub1) are proven in
     WINDOW_CAPS.md from [P1][P2][P3] alone -- bidegree valuation inductions
     on the C-recursion, D-transform arithmetic, and the shift identity.  No
     t-adic input anywhere.  ALT_REGIME.md itself consumes the k=2..5 sub1
     caps (deg d1<=9, deg sigma<=12, deg e<=15, deg d2<=6).
  5. BRIDGE_SWEEP.md sec.3's inapplicability argument, examined: its cap
     objection is against the SUB2 caps (deg e <= 10) -- alt states are SUB1
     states (deg e <= 15, satisfied with equality); its Phi/reduction
     objection concerns the CASCADE bookkeeping (F = t^210 G'), which the
     bridge never uses; the ALT_REGIME "do not transfer" row is about the
     cascade's per-level g_l caps / t-coupling, also never used here.

So: augment(alt-state, regime="sub1"), marked roots adjoined via q(r)=0,
saturation by the state's genuine nonzero scalars.  Soundness is inherited
from FULL_SYSTEM_BRIDGE.md sec.4 with "sub2" -> "sub1" and the five points
above (mechanically spot-checked in check_soundness()).

PILOT STATES (exact audited ansaetze; d2_threshold.py reconstruction, itself
audited by d2_threshold_verify.py):
  a12_b1110_T2 deg_d2=6 -- CONTROL: already killed exact char-0 by msolve
     (msolve_bridge_results.json); a fresh mechanism must reproduce it.
  a11_b3100_T2 deg_d2=6 -- THE PRIZE: the last open state of a 7/8 branch
     (tie-tower: 12GB RAM blowup; msolve: TIMEOUT).  A kill closes the first
     alt branch end-to-end under this mechanism.

All kills PENDING AUDIT.  New files only; read-only on all landed modules.
"""
from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

import sympy as sp

import convolution_descent as cd
import full_system_bridge as fsb
import modular_triage as mt
import t5_90t1_verify as base

ROOT = Path(__file__).resolve().parent
y = sp.Symbol("y")
E, S, r, r1, r2 = sp.symbols("E S r r1 r2")
QY = base.q  # 2048y^4 - 512y^3 + 320y^2 - 240y + 195

PRIMES = (10007, 10009, 100019)


def qpoly(v):
    return sp.expand(QY.subs(y, v))


# --------------------------------------------------------------------------
#  audited alt-state reconstruction (d2_threshold.py construction, verbatim;
#  that reconstruction is independently audited by d2_threshold_verify.py)
# --------------------------------------------------------------------------
def build_state(bid: str, deg_d2: int):
    if bid == "a11_b3100_T2":
        e_poly = E * (y + 1) ** 11 * (y - r1) ** 3 * (y - r2)
        sig_poly = S * (y + 1) ** 3 * (y - r1) ** 7 * (y - r2) ** 2
        root_vars = [r1, r2]
    elif bid == "a12_b1110_T2":
        comp = sp.div(sp.Poly(QY, y), sp.Poly(y - r, y))[0].as_expr() / 2048
        e_poly = E * (y + 1) ** 12 * comp
        sig_poly = S * (y + 1) ** 6 * comp ** 2
        root_vars = [r]
    else:
        raise ValueError(bid)
    Dc = list(sp.symbols(f"D0:{deg_d2 + 1}"))
    d2_poly = sum(Dc[i] * y ** i for i in range(deg_d2 + 1))
    # Reduce e and sigma mod q(r_i) UP FRONT (cheap on these small polys, and
    # it keeps every downstream product at r-degree <= 3; sound because
    # q(r_i)=0 is adjoined to the ideal).  The per-equation reduction of the
    # ~180 expanded G-coefficients is prohibitively slow -- residual r-powers
    # in the emitted system are left for Singular to reduce.
    red = reducer(root_vars)
    return (sp.expand(d2_poly), red(sig_poly), red(e_poly), root_vars, Dc)


def reducer(root_vars):
    """Reduce mod q(r_i) for every adjoined root (keeps r-degree <= 3)."""
    polys = [(rv, sp.Poly(qpoly(rv), rv)) for rv in root_vars]

    def red(e):
        e = sp.expand(e)
        for rv, QR in polys:
            e = sp.rem(sp.Poly(e, rv), QR).as_expr()
        return sp.expand(e)

    return red


# --------------------------------------------------------------------------
#  soundness spot-checks (the mechanical parts of the construction argument)
# --------------------------------------------------------------------------
def check_soundness():
    out = {}
    # (1) generators weighted-homogeneous with the landed weights
    out["homogeneity"] = fsb.check_homogeneity()
    # (3) alt's Phi~ = t^30 * (c*q) is the bridge's stripped Phi, identically
    phi_alt = sp.expand((y + 1) ** 30 * (fsb.C_GENUINE * QY))
    assert sp.expand(phi_alt - fsb.phi_stripped()) == 0
    out["phi_identical"] = True
    # (4) the caps consumed are the proven sub1 column of WINDOW_CAPS.md
    assert fsb.STRIP_DEGCAP["sub1"] == {"dm2": 18, "dm3": 21, "dm4": 24}
    out["sub1_caps"] = fsb.STRIP_DEGCAP["sub1"]
    # (5) the alt pilot states live INSIDE the sub1 window (deg caps 3k)
    for bid in ("a11_b3100_T2", "a12_b1110_T2"):
        d2p, sigp, ep, rvs, _ = build_state(bid, 6)
        assert sp.degree(ep, y) == 15      # = 3*5, sub1 cap, attained
        assert sp.degree(sigp, y) == 12    # = 3*4, sub1 cap, attained
        assert sp.degree(d2p, y) == 6      # = 3*2, sub1 cap, attained
    out["alt_states_in_sub1_window"] = True
    # complement divisor exactness for the a12 state: comp*(y-r)*2048 == q
    comp = sp.div(sp.Poly(QY, y), sp.Poly(y - r, y))[0].as_expr() / 2048
    red = reducer([r])
    assert red(sp.expand(comp * (y - r) * 2048 - QY)) == 0
    out["complement_exact_mod_q"] = True
    return out


# --------------------------------------------------------------------------
#  FAST builder: y-convolution with lazily-factored coefficients.
#  fsb.augment's sp.expand on whole G-polynomials is O(huge); here each
#  generator monomial is convolved per y-degree and only the final
#  per-coefficient expressions are expanded.  Generators still LOADED from
#  fsb.gsystem() and decomposed programmatically -- never hand-copied.
#  Validated against fsb.augment on a small synthetic state (check E below).
# --------------------------------------------------------------------------
def _ydict(expr):
    d = {}
    for mono, coeff in sp.Poly(sp.expand(expr), y).terms():
        d[mono[0]] = coeff
    return d


def _yconv(a, b):
    out = {}
    for i, ca in a.items():
        for j, cb in b.items():
            out[i + j] = out.get(i + j, 0) + ca * cb
    return out


def _ypow(a, n):
    out = {0: sp.Integer(1)}
    for _ in range(n):
        out = _yconv(out, a)
    return out


def fast_gsystem_equations(subs_dicts):
    """Every y-coefficient of each generator under the substitution
    var -> ydict.  subs_dicts: {fsb-symbol: ydict}."""
    eqs = []
    for _name, g in fsb.gsystem().items():
        acc = {}
        for term in sp.Add.make_args(sp.expand(g)):
            coeff = sp.Integer(1)
            cur = {0: sp.Integer(1)}
            for b_, ex in term.as_powers_dict().items():
                if b_.is_number:
                    coeff *= b_ ** ex
                    continue
                cur = _yconv(cur, _ypow(subs_dicts[b_], int(ex)))
            for i, cv in cur.items():
                acc[i] = acc.get(i, 0) + coeff * cv
        for i in sorted(acc):
            e = sp.expand(acc[i])
            if e != 0:
                eqs.append(e)
    return eqs


def build_bridge_fast(bid: str, deg_d2: int, *, regime: str = "sub1"):
    d2p, sigp, ep, root_vars, Dc = build_state(bid, deg_d2)
    red = reducer(root_vars)
    d0p = red(sp.expand((d2p ** 2 + sigp) / 4))
    spare, spare_unk = fsb.build_spare(regime)
    subs_dicts = {
        fsb.D2: _ydict(d2p), fsb.D1: {},
        fsb.D0: _ydict(d0p), fsb.DM1: _ydict(ep),
        fsb.DM2: _ydict(spare[fsb.DM2]), fsb.DM3: _ydict(spare[fsb.DM3]),
        fsb.DM4: _ydict(spare[fsb.DM4]), fsb.PHI: _ydict(fsb.phi_stripped()),
    }
    eqs_raw = fast_gsystem_equations(subs_dicts)
    eqs, seen = [], set()
    for eq in eqs_raw:
        key = sp.sstr(eq)
        if key not in seen:
            seen.add(key)
            eqs.append(eq)
    for rv in root_vars:
        eqs.append(qpoly(rv))
    sat = [E, S, Dc[-1]]
    if len(root_vars) == 2:
        sat.append(root_vars[0] - root_vars[1])
    ring_vars = mt.ring_vars_of(eqs, extra=[s for f in sat
                                            for s in sp.sympify(f).free_symbols])
    return {
        "label": f"{bid}_degd2_{deg_d2}", "bid": bid, "deg_d2": deg_d2,
        "regime": regime, "equations": eqs, "sat_factors": sat,
        "ring_vars": ring_vars, "n_equations": len(eqs),
        "n_unknowns": len(ring_vars), "n_spare": len(spare_unk),
        "builder": "fast-yconv",
    }


def check_fast_builder():
    """Cross-validate the fast builder against fsb.augment on a small
    synthetic sub2 state (cheap for augment).  Equation SETS must agree."""
    Dc = sp.symbols("D0")
    d2p = Dc
    sigp = sp.expand(S * (y + 1) ** 2)
    ep = sp.expand(E * (y + 1) ** 3)
    ansatz = cd.build_ansatz(d2=d2p, d1=sp.Integer(0), sigma=sigp, e=ep,
                             unknowns=(Dc, E, S), parameters=())
    aug = fsb.augment(ansatz, regime="sub2")
    slow = {sp.sstr(sp.expand(e)) for e in aug["equations"] if e != 0}
    spare, _ = fsb.build_spare("sub2")
    d0p = sp.expand((d2p ** 2 + sigp) / 4)
    subs_dicts = {
        fsb.D2: _ydict(d2p), fsb.D1: {}, fsb.D0: _ydict(d0p),
        fsb.DM1: _ydict(ep),
        fsb.DM2: _ydict(spare[fsb.DM2]), fsb.DM3: _ydict(spare[fsb.DM3]),
        fsb.DM4: _ydict(spare[fsb.DM4]), fsb.PHI: _ydict(fsb.phi_stripped()),
    }
    fast = {sp.sstr(e) for e in fast_gsystem_equations(subs_dicts)}
    assert fast == slow, (len(fast), len(slow),
                          list(fast - slow)[:2], list(slow - fast)[:2])
    return True


def build_bridge(bid: str, deg_d2: int, *, regime: str = "sub1"):
    d2p, sigp, ep, root_vars, Dc = build_state(bid, deg_d2)
    ansatz = cd.build_ansatz(
        d2=d2p, d1=sp.Integer(0), sigma=sigp, e=ep,
        unknowns=tuple(Dc) + (E, S), parameters=tuple(root_vars))
    aug = fsb.augment(ansatz, regime=regime)
    eqs, seen = [], set()
    for eq in aug["equations"]:
        if eq == 0:
            continue
        key = sp.sstr(eq)
        if key in seen:
            continue
        seen.add(key)
        eqs.append(eq)
    for rv in root_vars:
        eqs.append(qpoly(rv))
    sat = [E, S, Dc[-1]]
    if len(root_vars) == 2:
        sat.append(root_vars[0] - root_vars[1])
    ring_vars = mt.ring_vars_of(eqs, extra=[s for f in sat
                                            for s in sp.sympify(f).free_symbols])
    return {
        "label": f"{bid}_degd2_{deg_d2}", "bid": bid, "deg_d2": deg_d2,
        "regime": regime, "equations": eqs, "sat_factors": sat,
        "ring_vars": ring_vars, "n_equations": len(eqs),
        "n_unknowns": len(ring_vars),
        "n_spare": len(aug["spare_unknowns"]),
    }


# --------------------------------------------------------------------------
#  emission + orphan-proof runner (WSL-side timeout; Windows relay death can
#  otherwise leave the WSL Singular alive)
# --------------------------------------------------------------------------
def emit(eqs, sat_factors, ring_vars, *, char: int):
    w = sp.Symbol("w")
    rv = list(ring_vars) + [w]
    var_txt = ",".join(v.name for v in rv)
    lines = [f"ring R = {char},({var_txt}),dp;"]
    members = []
    for i, g in enumerate(eqs):
        s = (mt.poly_to_singular_modp(g, rv, char) if char > 0
             else fsb._to_singular(g, rv))
        if s in ("0", ""):
            continue
        lines.append(f"poly g{i} = {s};")
        members.append(f"g{i}")
    prod = sp.Integer(1)
    for f in sat_factors:
        prod = prod * f
    s = (mt.poly_to_singular_modp(sp.expand(w * prod - 1), rv, char)
         if char > 0 else fsb._to_singular(sp.expand(w * prod - 1), rv))
    lines.append(f"poly nz = {s};")
    members.append("nz")
    lines.append(f"ideal I = {','.join(members)};")
    lines.append("ideal G = std(I);")
    lines.append("int u = (reduce(1,G)==0);")
    lines.append('"@@UNIT";')
    lines.append("u;")
    lines.append('"@@DIM";')
    lines.append("dim(G);")
    lines.append("quit;")
    return "\n".join(lines) + "\n"


def run_singular(program: str, timeout: float) -> dict:
    """WSL-side `timeout` guards the Singular process itself."""
    tsec = max(5, int(timeout))
    cmd = ("wsl.exe", "-d", "Ubuntu", "--", "bash", "-lc",
           f"cd $HOME && timeout {tsec}s Singular -q")
    t0 = time.monotonic()
    try:
        cp = subprocess.run(cmd, input=program, text=True, encoding="utf-8",
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            timeout=timeout + 20, check=False)
    except subprocess.TimeoutExpired:
        return {"status": "timeout", "verdict": "TIMEOUT", "unit": None,
                "dim": None, "wall": round(time.monotonic() - t0, 1)}
    wall = round(time.monotonic() - t0, 1)
    combined = (cp.stdout or "") + "\n" + (cp.stderr or "")
    combined = combined.replace("\x00", "")
    import re
    um = re.search(r"@@UNIT\s*\r?\n\s*(-?\d+)", combined)
    dm = re.search(r"@@DIM\s*\r?\n\s*(-?\d+)", combined)
    if um is None:
        # WSL-side timeout kills Singular before it prints -> honest TIMEOUT
        return {"status": "timeout", "verdict": "TIMEOUT", "unit": None,
                "dim": None, "wall": wall}
    unit = bool(int(um.group(1)))
    return {"status": "ok", "verdict": "UNIT" if unit else "PROPER",
            "unit": unit, "dim": None if dm is None else int(dm.group(1)),
            "wall": wall}


def triage_numroot(bs, *, primes=PRIMES, timeout=90.0):
    """Mod-p reconnaissance with the marked roots specialized to numeric
    roots of q mod p (removes the r-nonlinearity that swells symbolic std).
    Evidence over F_p-bar, not a Q certificate (bridge_sweep.py pattern)."""
    root_syms = [v for v in bs["ring_vars"] if v.name in ("r", "r1", "r2")]
    nroots = len(root_syms)
    qeqs = {sp.sstr(qpoly(rv)) for rv in root_syms}
    eqs0 = [e for e in bs["equations"] if sp.sstr(e) not in qeqs]
    out = []
    for p in primes:
        roots = mt.q_roots_mod_p(p, nroots)
        if len(roots) < nroots:
            out.append({"verdict": "SKIP", "prime": p})
            print(f"    (numroot) p={p}: SKIP (q has <{nroots} roots)",
                  flush=True)
            continue
        subst = dict(zip(root_syms, roots[:nroots]))
        eqs = [sp.expand(sp.sympify(e).subs(subst)) for e in eqs0]
        sat = [sp.expand(sp.sympify(f).subs(subst)) for f in bs["sat_factors"]]
        sat = [f for f in sat if f.free_symbols]
        rv = mt.ring_vars_of(eqs, extra=[s for f in sat
                                         for s in f.free_symbols])
        rr = run_singular(emit(eqs, sat, rv, char=p), timeout)
        rr["prime"] = p
        rr["roots"] = [int(x) for x in roots[:nroots]]
        out.append(rr)
        print(f"    (numroot) p={p}: {rr['verdict']} dim={rr.get('dim')} "
              f"({rr['wall']}s)", flush=True)
    return out


_MS = "$HOME/msolve/msolve"


def exact_msolve(bs, *, timeout=300.0, char=0):
    """Exact char-0 verdict via msolve (F4 + multi-modular).  Output [-1]
    means empty variety over Qbar == KILL (bridge_sweep.py pattern)."""
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
    fname, outname = f"altbridge_{tag}.ms", f"altbridge_{tag}.out"
    subprocess.run(("wsl.exe", "-d", "Ubuntu", "--", "bash", "-lc",
                    f"cat > $HOME/{fname}"), input=prog, text=True,
                   encoding="utf-8", check=True)
    tsec = max(5, int(timeout))
    run = ("wsl.exe", "-d", "Ubuntu", "--", "bash", "-lc",
           f"cd $HOME && timeout {tsec}s {_MS} -f $HOME/{fname} "
           f"-o $HOME/{outname}; cat $HOME/{outname} 2>/dev/null")
    t0 = time.monotonic()
    try:
        cp = subprocess.run(run, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, text=True,
                            encoding="utf-8", timeout=timeout + 30,
                            check=False)
    except subprocess.TimeoutExpired:
        return {"verdict": "TIMEOUT", "wall": round(time.monotonic() - t0, 1),
                "engine": "msolve"}
    body = (cp.stdout or "").strip()
    wall = round(time.monotonic() - t0, 1)
    import re as _re
    if body.startswith("[-1]"):
        verdict = "UNIT"
    elif _re.match(r"\[1,\s*\d+,\s*-1", body):
        verdict = "INFINITE_SOL"
    elif body:
        verdict = "HAS_SOL"
    else:
        verdict = "TIMEOUT"
    return {"verdict": verdict, "wall": wall, "engine": "msolve",
            "out_head": body[:120]}


def triage(bs, *, primes=PRIMES, timeout=45.0):
    out = []
    for p in primes:
        rr = run_singular(emit(bs["equations"], bs["sat_factors"],
                               bs["ring_vars"], char=p), timeout)
        rr["prime"] = p
        out.append(rr)
        print(f"    p={p}: {rr['verdict']} dim={rr.get('dim')} "
              f"({rr['wall']}s)", flush=True)
    return out


def exact_kill(bs, *, timeout=300.0):
    rr = run_singular(emit(bs["equations"], bs["sat_factors"],
                           bs["ring_vars"], char=0), timeout)
    print(f"    exact Q: {rr['verdict']} dim={rr.get('dim')} "
          f"({rr['wall']}s)", flush=True)
    return rr


# --------------------------------------------------------------------------
#  driver
# --------------------------------------------------------------------------
def run_case(bid, deg_d2, *, triage_timeout=45.0, exact_timeout=300.0):
    print(f"[{bid} deg_d2={deg_d2}] building bridge system...", flush=True)
    t0 = time.monotonic()
    bs = build_bridge_fast(bid, deg_d2)
    print(f"  built: {bs['n_equations']} equations, {bs['n_unknowns']} ring "
          f"vars ({bs['n_spare']} spare) in "
          f"{round(time.monotonic() - t0, 1)}s", flush=True)
    rec = {k: bs[k] for k in ("label", "bid", "deg_d2", "regime",
                              "n_equations", "n_unknowns", "n_spare")}
    print("  mod-p triage:", flush=True)
    rec["triage"] = triage(bs, timeout=triage_timeout)
    verdicts = {t["verdict"] for t in rec["triage"]}
    if verdicts == {"UNIT"}:
        print("  all primes UNIT -> exact Q attempt:", flush=True)
        rec["exact"] = exact_kill(bs, timeout=exact_timeout)
        rec["verdict"] = ("KILLED_EXACT_Q" if rec["exact"]["verdict"] == "UNIT"
                          else ("MODP_UNIT_EXACT_" + rec["exact"]["verdict"]))
    elif "PROPER" in verdicts:
        rec["verdict"] = "PROPER_MODP"   # loud: would be a survival signal
        print("  !! PROPER at some prime -- NOT attempting exact; "
              "survival-signal candidate", flush=True)
    else:
        print("  triage inconclusive (timeouts) -> exact Q attempt anyway:",
              flush=True)
        rec["exact"] = exact_kill(bs, timeout=exact_timeout)
        rec["verdict"] = ("KILLED_EXACT_Q" if rec["exact"]["verdict"] == "UNIT"
                          else "COST")
    return rec


def main(argv):
    which = argv[1] if len(argv) > 1 else "pilot"
    results = {"schema": "alt-bridge-v1", "soundness": None, "cases": []}
    print("== soundness spot-checks ==", flush=True)
    sc = check_soundness()
    results["soundness"] = {k: (v if not isinstance(v, dict) else dict(v))
                            for k, v in sc.items()}
    print(f"  {sc}", flush=True)
    print("== fast-builder cross-validation vs fsb.augment ==", flush=True)
    t0 = time.monotonic()
    assert check_fast_builder()
    results["fast_builder_validated"] = True
    print(f"  equation sets IDENTICAL ({round(time.monotonic()-t0,1)}s)",
          flush=True)
    if which == "pilot":
        # control first (known exact kill via msolve), then the prize
        for bid, dd in (("a12_b1110_T2", 6), ("a11_b3100_T2", 6)):
            results["cases"].append(run_case(bid, dd))
            Path(ROOT / "alt_bridge_results.json").write_text(
                json.dumps(results, indent=1, default=str), encoding="utf-8")
    elif which == "retry":
        # swell fallback for the pilot cases: numroot mod-p triage + msolve
        # char-0 exact (the engines that beat the marked-root swell before)
        for bid, dd in (("a12_b1110_T2", 6), ("a11_b3100_T2", 6)):
            print(f"[retry {bid} deg_d2={dd}] building...", flush=True)
            bs = build_bridge_fast(bid, dd)
            rec = {k: bs[k] for k in ("label", "bid", "deg_d2", "regime",
                                      "n_equations", "n_unknowns", "n_spare")}
            print("  numroot mod-p triage:", flush=True)
            rec["triage_numroot"] = triage_numroot(bs)
            print("  exact char-0 msolve:", flush=True)
            rec["exact_msolve"] = exact_msolve(bs, timeout=300.0)
            print(f"    msolve: {rec['exact_msolve']['verdict']} "
                  f"({rec['exact_msolve']['wall']}s)", flush=True)
            nr = {t["verdict"] for t in rec["triage_numroot"]
                  if t["verdict"] != "SKIP"}
            ms = rec["exact_msolve"]["verdict"]
            if ms == "UNIT":
                rec["verdict"] = "KILLED_EXACT_Q_MSOLVE"
            elif nr == {"UNIT"}:
                rec["verdict"] = "MODP_UNIT_NUMROOT_EXACT_PENDING"
            elif "PROPER" in nr or ms in ("HAS_SOL", "INFINITE_SOL"):
                rec["verdict"] = "SURVIVAL_SIGNAL"
            else:
                rec["verdict"] = "COST"
            results["cases"].append(rec)
            Path(ROOT / "alt_bridge_results.json").write_text(
                json.dumps(results, indent=1, default=str), encoding="utf-8")
    elif which == "sweep":
        # extend across the two branches' remaining deg_d2 states as
        # cross-checks (census-killed; a second mechanism's agreement)
        for bid in ("a12_b1110_T2", "a11_b3100_T2"):
            for dd in (5, 4, 3):
                results["cases"].append(run_case(bid, dd))
                Path(ROOT / "alt_bridge_results.json").write_text(
                    json.dumps(results, indent=1, default=str),
                    encoding="utf-8")
    print("\n== census ==", flush=True)
    for c in results["cases"]:
        print(f"  {c['label']}: {c['verdict']}", flush=True)
    Path(ROOT / "alt_bridge_results.json").write_text(
        json.dumps(results, indent=1, default=str), encoding="utf-8")


if __name__ == "__main__":
    main(sys.argv)
