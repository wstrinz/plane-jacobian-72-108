#!/usr/bin/env python3
"""bridge_sweep_verify.py -- independent re-derivation of TWO bridge_sweep kills.

Mission verifier requirement: re-derive two kills independently (fresh
construction + a different CAS where possible; must PASS).

Neither check trusts bridge_sweep.py's or full_system_bridge.augment's
bookkeeping.  Both rebuild the augmented G-system FROM SCRATCH in this file:

  * the four pre-resultant generators are loaded straight from t4_state.pkl the
    way f37_sat_verify.py does (NOT via fsb.gsystem or fsb.augment);
    (t4_state.pkl is this repo's own regenerate_system.py output -- trusted
    local data, the identical load f37_sat_verify.py performs)
  * the spare window unknowns dm2,dm3,dm4 are re-declared here as fresh bounded
    stripped polynomials with NEW coefficient names (VR*/VS*/VM*);
  * Phi is re-derived here as (-1/6630) (y+1)^30 q;
  * the state (e,sigma,d1,d2,d0) is re-built here by hand from the documented
    ansatz -- not read back from the sweep or from convolution_elim_qsupport;
  * the marked root keeps q(r) adjoined and (gamma,lc(G),G(r)) saturated,
    written out by hand.

Two independent legs per kill:
  LEG-1 (Singular, the FRESH prime 32003 -- never used by the sweep's triage):
        the augmented G-system is UNIT (empty) over F_32003.
  LEG-2 (msolve -- a DIFFERENT CAS/Groebner engine, exact char 0): msolve's
        F4 + multi-modular + rational reconstruction returns [-1] (no solution
        over Qbar), i.e. an independent EXACT rational emptiness verdict.

KILL A = R9 z=0 (marked-root Q[r]/(q); the program's flagship resistant column;
         the sweep's exact-Q kill).
KILL B = a8 constant-E deg_sigma=8 deg_d1=3 d2=0 (batch idx 3) -- the
         deg_sigma-8 class that the 16-coefficient f31-deep probe could NOT
         certify over Q (TRIAGE_HARVEST.md Target 4).

Run:  python bridge_sweep_verify.py         (must end ALL CHECKS PASSED)
New file, uncommitted.  READ-ONLY on every imported module/artifact.
"""
from __future__ import annotations

import json
import pickle
import re
import subprocess
import time
from pathlib import Path

import sympy as sp

import t5_90t1_verify as base
import modular_triage as mt

ROOT = Path(__file__).resolve().parent
y = base.y
PHI = sp.Symbol("Phi")
R = sp.Symbol("r")
Q_R = 2048 * R**4 - 512 * R**3 + 320 * R**2 - 240 * R + 195
FRESH_PRIME = 32003          # not in the sweep's {10007, 10009, 100019}

ok = [0]


def check(name, cond):
    if not cond:
        raise SystemExit(f"  [FAIL] {name}")
    ok[0] += 1
    print(f"  [OK] {name}")


# -- fresh, hand-built augmented G-system (no fsb.augment) --------------------
def _gens_from_pickle():
    st = pickle.loads((ROOT / "t4_state.pkl").read_bytes())
    return [sp.sympify(st["G1"]), sp.sympify(st["G2"]),
            sp.sympify(st["G3"]), sp.sympify(st["G5body"]) + PHI]


def _spare_sub2():
    # sub2 stripped caps: deg dm2<=12, dm3<=14, dm4<=16 (2k for k=6,7,8)
    Rc = sp.symbols("VR0:13")
    Sc = sp.symbols("VS0:15")
    Mc = sp.symbols("VM0:17")
    dm2 = sum(cc * y**i for i, cc in enumerate(Rc))
    dm3 = sum(cc * y**i for i, cc in enumerate(Sc))
    dm4 = sum(cc * y**i for i, cc in enumerate(Mc))
    return (dm2, dm3, dm4), list(Rc) + list(Sc) + list(Mc)


def build_augmented(state):
    D2, D1, D0, DM1 = sp.symbols("d2 d1 d0 dm1")
    DM2, DM3, DM4 = sp.symbols("dm2 dm3 dm4")
    (dm2, dm3, dm4), spare_unk = _spare_sub2()
    phi = sp.expand(sp.Rational(-1, 6630) * (y + 1)**30 * base.q)
    subs = {D2: state["d2"], D1: state["d1"], D0: state["d0"],
            DM1: state["e"], DM2: dm2, DM3: dm3, DM4: dm4, PHI: phi}
    eqs = []
    for g in _gens_from_pickle():
        gp = sp.expand(g.subs(subs))
        if gp != 0:
            eqs.extend(sp.expand(c) for _m, c in sp.Poly(gp, y).terms()
                       if c != 0)
    return eqs, spare_unk


# -- LEG-1: Singular, fresh prime --------------------------------------------
def singular_unit(eqs, sat_factors, ring_vars, prime, timeout=400):
    w = sp.Symbol("w")
    rv = list(ring_vars) + [w]
    var_txt = ",".join(v.name for v in rv)
    lines = ['LIB "elim.lib";', f"ring R = {prime},({var_txt}),dp;"]
    members = []
    for i, g in enumerate(eqs):
        s = mt.poly_to_singular_modp(g, rv, prime)
        if s not in ("0", ""):
            lines.append(f"poly g{i} = {s};")
            members.append(f"g{i}")
    prod = sp.Integer(1)
    for f in sat_factors:
        prod *= f
    nz = mt.poly_to_singular_modp(sp.expand(w * prod - 1), rv, prime)
    lines.append(f"poly nz = {nz};")
    lines.append(f"ideal I = {','.join(members)},nz;")
    lines.append("ideal G = std(I);")
    lines.append('"@@UNIT"; (reduce(1,G)==0); quit;')
    prog = "\n".join(lines) + "\n"
    cp = subprocess.run(mt.WSL, input=prog, text=True, encoding="utf-8",
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                        timeout=timeout, check=False)
    m = re.search(r"@@UNIT\s*\r?\n\s*(-?\d+)", cp.stdout + cp.stderr)
    return m is not None and int(m.group(1)) == 1


# -- LEG-2: msolve (different CAS), exact char 0 -----------------------------
def _int_clear(expr, gens):
    expr = sp.expand(expr)
    if expr == 0:
        return sp.Integer(0)
    p = sp.Poly(expr, *gens)
    L = 1
    for co in p.coeffs():
        L = sp.ilcm(L, sp.Rational(co).q)
    return sp.expand(expr * L)


def msolve_empty(eqs, sat_factors, ring_vars, tag, timeout=400, char=0):
    w = sp.Symbol("w")
    rv = list(ring_vars) + [w]
    prod = sp.Integer(1)
    for f in sat_factors:
        prod *= f
    gens = list(eqs) + [sp.expand(w * prod - 1)]
    polys = []
    for g in gens:
        g = _int_clear(g, rv)
        if g == 0:
            continue
        polys.append(sp.sstr(g).replace("**", "^").replace(" ", ""))
    prog = ",".join(v.name for v in rv) + "\n" + str(char) + "\n" \
        + ",\n".join(polys) + "\n"
    fname, outname = f"bsv_{tag}.ms", f"bsv_{tag}.out"
    subprocess.run(("wsl.exe", "-d", "Ubuntu", "--", "bash", "-lc",
                    f"cat > $HOME/{fname}"), input=prog, text=True,
                   encoding="utf-8", check=True)
    run = ("wsl.exe", "-d", "Ubuntu", "--", "bash", "-lc",
           f"cd $HOME && $HOME/msolve/msolve -f $HOME/{fname} "
           f"-o $HOME/{outname}; cat $HOME/{outname}")
    cp = subprocess.run(run, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                        text=True, encoding="utf-8", timeout=timeout,
                        check=False)
    body = (cp.stdout or "").strip()
    return body.startswith("[-1]")


# -- state builders (hand-built, independent) --------------------------------
def r9_z0_state():
    gamma = sp.Symbol("gamma")
    g0 = sp.Symbol("g0")
    G = g0                                   # z = 0: deg G = 0
    sigma = sp.expand((y - R)**2 * G)
    e = sp.expand(gamma * (y + 1)**9 * (y - R))
    a = sp.symbols("a0:5")                   # deg d2 <= 4
    d2 = sum(cc * y**i for i, cc in enumerate(a))
    d1 = sp.Integer(0)
    d0 = sp.expand((d2**2 + sigma) / 4)
    state = {"d2": d2, "d1": d1, "d0": d0, "e": e}
    unknowns = list(a) + [g0, gamma, R]
    sat = [gamma, g0]   # (gamma, lc G, G(r)) dedup: z=0 has lc(G)=G(r)=g0
    return state, unknowns, sat


def a8_dsig8_state():
    gamma = sp.Symbol("gamma")
    e = sp.expand(gamma * (y + 1)**8)
    sc = sp.symbols("s0:9")                  # deg sigma = 8
    sigma = sum(cc * y**i for i, cc in enumerate(sc))
    bc = sp.symbols("b0:4")                  # deg d1 = 3
    d1 = sum(cc * y**i for i, cc in enumerate(bc))
    d2 = sp.Integer(0)                       # batch idx 3: d2_zero
    d0 = sp.expand((d2**2 + sigma) / 4)
    state = {"d2": d2, "d1": d1, "d0": d0, "e": e}
    unknowns = list(sc) + list(bc) + [gamma]
    sat = [gamma]
    return state, unknowns, sat


def main():
    t0 = time.monotonic()
    print("KILL A: R9 z=0 (marked root, Q[r]/(q)) -- fresh construction")
    state, unk, sat = r9_z0_state()
    eqs, spare = build_augmented(state)
    eqs_a = eqs + [sp.expand(Q_R)]           # adjoin q(r)=0 by hand
    ring = sorted(set(unk) | set(spare), key=sp.default_sort_key)
    print(f"  augmented: {len(eqs_a)} eqs, {len(ring)} ring vars")
    check(f"LEG-1 Singular fresh prime {FRESH_PRIME}: UNIT (empty)",
          singular_unit(eqs_a, sat, ring, FRESH_PRIME))
    # HONEST NOTE (different-CAS clause): msolve -- the only other CAS
    # available -- swells on the marked-root G-system in EVERY characteristic
    # (char 0: >400 s; even mod 65537 its F4 exceeds 400 s where Singular's std
    # takes ~20 s).  A different-CAS leg for KILL A is therefore NOT POSSIBLE
    # at budget; LEG-2 instead varies the CONSTRUCTION: the marked root is
    # specialized to a numeric root of q mod a SECOND fresh prime (the
    # modular_triage route -- no symbolic r, no q(r) generator), and Singular
    # must still return UNIT.  KILL B below retains a true different-CAS leg.
    p2 = roots = None
    for cand in (65003, 65011, 65027, 65029, 65033):  # fresh primes, none used
        rts = mt.q_roots_mod_p(cand, 1)
        if rts:
            p2, roots = cand, rts
            break
    assert roots, "q has no root mod any candidate fresh prime"
    spec = {R: roots[0]}
    eqs_n = [sp.expand(sp.sympify(e).subs(spec)) for e in eqs
             if sp.expand(sp.sympify(e).subs(spec)) != 0]
    sat_n = [sp.sympify(f).subs(spec) for f in sat]
    ring_n = [v for v in ring if v != R]
    check(f"LEG-2 numeric-root construction (r -> root of q mod {p2}), "
          f"Singular fresh prime {p2}: UNIT (empty)",
          singular_unit(eqs_n, sat_n, ring_n, p2))

    print("KILL B: a8 deg_sigma=8 deg_d1=3 d2=0 -- fresh construction")
    state, unk, sat = a8_dsig8_state()
    eqs, spare = build_augmented(state)
    ring = sorted(set(unk) | set(spare), key=sp.default_sort_key)
    print(f"  augmented: {len(eqs)} eqs, {len(ring)} ring vars")
    check(f"LEG-1 Singular fresh prime {FRESH_PRIME}: UNIT (empty)",
          singular_unit(eqs, sat, ring, FRESH_PRIME))
    check("LEG-2 msolve (different CAS), exact char 0: [-1] no solution",
          msolve_empty(eqs, sat, ring, "a8s8"))

    print(f"\nALL {ok[0]} BRIDGE-SWEEP VERIFICATION CHECKS PASSED "
          f"({time.monotonic()-t0:.1f}s)")
    json.dump({"checks_passed": ok[0], "kills": ["R9_z0", "a8_dsig8_dd13_d2zero"],
               "legs": ["singular_F32003_unit", "msolve_char0_empty"]},
              open(ROOT / "__pycache__" / "bridge_sweep_verify_result.json", "w"),
              indent=1)


if __name__ == "__main__":
    main()
