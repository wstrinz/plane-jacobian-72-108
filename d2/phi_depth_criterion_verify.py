#!/usr/bin/env python3
"""phi_depth_criterion_verify.py  --  independent checker for the Phi-depth kill

Re-derives, from primitive data and WITHOUT importing the criterion's own
arithmetic, every fact the Phi-window-depth kill rests on, then links the
closed-form criterion to the exact symbolic computation on sampled states.

Checks:
  1. STRIP CONSISTENCY.  Every term of G5body has u-weight 17, hence window floor
     12*17 = 204, hence comparable with Phi on one y-degree axis.  Recomputed here
     from the symbol weights, independently of `face_kill_sweep`.
  2. PHI.  The stripped Phi = c*t^30*q (c = -1/6630) has y-degree EXACTLY 34 and
     nonzero leading coefficient; 34 equals the sub2 stripped cap 2*17, i.e. Phi
     sits at its cap (the C6 tightness statement).
  3. CAPS.  `full_system_bridge.STRIP_DEGCAP` equals (14-12)k for sub2 and
     (15-12)k for sub1 at k = 6,7,8.
  4. CANONICAL G5.  The generator used is `G5body + Phi` (NOT `2*Phi + G5body`),
     matching `full_system_bridge.py` and the C11 certificate in
     `f37_sat_verify.py`.
  5. KILL SAMPLES (load-bearing).  For sampled states the criterion calls KILLED,
     expand the canonical G5 symbolically and assert its y-degree-34 coefficient
     is EXACTLY lc(Phi~) -- i.e. no G5body term contributes there.
  6. NON-KILL SAMPLES.  For sampled states the criterion does NOT call killed,
     assert some G5body term genuinely reaches degree >= 34, so the criterion is
     not merely failing to fire.
  7. ROLL-UP.  Per-window state counts reconstruct the declared `state_total`.

Exit 0 iff every check passes.  `--quiet` prints only failures + the summary.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__)) or "."
sys.path.insert(0, HERE)

import full_system_bridge as fsb                     # noqa: E402
import phi_depth_criterion as P                      # noqa: E402
import face_kill_sweep as F                          # noqa: E402
import bigrade_annotator as ba                       # noqa: E402
from bigrade_annotator import y                      # noqa: E402

FAILURES = []
CHECKS = [0]


def check(label, ok, detail=""):
    CHECKS[0] += 1
    if not ok:
        FAILURES.append("%s%s" % (label, (" -- " + detail) if detail else ""))
    return ok


# independent copy of the window-symbol u-weights (NOT imported from the modules
# under test): w(d_{4-k}) = k, and Phi carries weight 17.
U_WEIGHT = {"d0": 4, "d1": 3, "d2": 2, "dm1": 5, "dm2": 6, "dm3": 7, "dm4": 8}
G5BODY = (("d0", "dm1", "dm4"), ("d0", "dm2", "dm3"), ("d1", "dm2", "dm4"),
          ("d1", "dm3", "dm3"), ("d2", "dm3", "dm4"), ("dm1", "dm1", "dm3"),
          ("dm1", "dm2", "dm2"))
PHI_WEIGHT = 17
FLOOR_SLOPE = 12


def check_strip_consistency(verbose):
    ok = True
    for factors in G5BODY:
        w = sum(U_WEIGHT[f] for f in factors)
        ok &= check("strip: %s u-weight" % "*".join(factors), w == PHI_WEIGHT,
                    "got %d, expected %d" % (w, PHI_WEIGHT))
        ok &= check("strip: %s floor" % "*".join(factors),
                    FLOOR_SLOPE * w == FLOOR_SLOPE * PHI_WEIGHT)
    if verbose:
        print("  1. strip consistency: all %d G5body terms at u-weight %d, floor y^%d  %s"
              % (len(G5BODY), PHI_WEIGHT, FLOOR_SLOPE * PHI_WEIGHT, "OK" if ok else "FAIL"))
    return ok


def check_phi(verbose):
    # rebuild Phi~ from primitives rather than calling the modules under test
    c = sp.Rational(-1, 6630)
    q = 2048 * y**4 - 512 * y**3 + 320 * y**2 - 240 * y + 195
    phi = sp.expand(c * (y + 1)**30 * q)
    deg = sp.Poly(phi, y).degree()
    lc = sp.Poly(phi, y).LC()
    ok = check("phi degree == 34", deg == P.PHI_STRIPPED_DEGREE,
               "got %s" % deg)
    ok &= check("phi leading coeff nonzero", lc != 0)
    ok &= check("phi sits at sub2 cap 2*17", deg == 2 * PHI_WEIGHT,
                "deg %s vs cap %d" % (deg, 2 * PHI_WEIGHT))
    # and it must agree with the annotator's own stripped Phi
    ok &= check("phi matches bigrade_annotator._phi_stripped",
                sp.expand(phi - ba._phi_stripped()) == 0)
    if verbose:
        print("  2. Phi~: degree %s (= sub2 cap 2*17), lc = %s, matches annotator  %s"
              % (deg, lc, "OK" if ok else "FAIL"))
    return ok, lc


def check_caps(verbose):
    ok = True
    for window, slope in (("sub2", 14 - 12), ("sub1", 15 - 12)):
        for name, k in (("dm2", 6), ("dm3", 7), ("dm4", 8)):
            got = fsb.STRIP_DEGCAP[window][name]
            ok &= check("cap %s/%s" % (window, name), got == slope * k,
                        "got %d, expected %d" % (got, slope * k))
    if verbose:
        print("  3. caps: sub2 %s / sub1 %s  %s"
              % (fsb.STRIP_DEGCAP["sub2"], fsb.STRIP_DEGCAP["sub1"],
                 "OK" if ok else "FAIL"))
    return ok


def check_canonical_g5(verbose):
    Phi = sp.Symbol("Phi")
    gens = F.canonical_G_generators()
    coeff = sp.expand(gens["G5"][0]).coeff(Phi)
    ok = check("canonical G5 uses Phi coefficient 1", coeff == 1,
               "got %s (the 2*Phi form is FULL_SYSTEM_BRIDGE.md:62's error)" % coeff)
    if verbose:
        print("  4. canonical G5 = G5body + %s*Phi  %s" % (coeff, "OK" if ok else "FAIL"))
    return ok


def _generic_state(state, case):
    """Rebuild generic state polynomials from the recorded degrees."""
    def poly(prefix, deg):
        if deg is None or deg < 0:
            return sp.Integer(0)
        return sum(sp.Symbol("%s%d" % (prefix, i)) * y**i for i in range(deg + 1))

    d2 = sp.Integer(0) if case.get("d2_zero") else poly("a", state.get("deg_d2"))
    d1 = poly("b", state.get("deg_d1"))
    sigma = sp.Integer(0) if case.get("sigma_zero") else poly("s", state.get("deg_sigma"))
    deg_e = state.get("deg_e")
    # e is degree-exact with a nonzero leading coefficient; any concrete
    # representative of that degree exercises the SAME degree bookkeeping.
    e = sp.expand(sp.Symbol("gamma") * (y + 1)**int(deg_e))
    return {"d2": d2, "d1": d1, "sigma": sigma, "e": e}


def _g5_top_coefficient(polys, window):
    """Exact y-degree-34 coefficient of the canonical stripped G5 on this state."""
    caps = fsb.STRIP_DEGCAP[window]
    d0s, d1s, d2s, dm1, dm2, dm3, dm4, Phi = ba._gsystem_symbols()
    spare = {}
    for sym, name, pre in ((dm2, "dm2", "R"), (dm3, "dm3", "S"), (dm4, "dm4", "T")):
        cap = caps[name]
        spare[sym] = sum(sp.Symbol("%s%d" % (pre, i)) * y**i for i in range(cap + 1))
    subs = {d2s: polys["d2"], d1s: polys["d1"], dm1: polys["e"],
            d0s: sp.expand((polys["d2"]**2 + polys["sigma"]) / 4),
            dm2: spare[dm2], dm3: spare[dm3], dm4: spare[dm4],
            Phi: ba._phi_stripped()}
    g5 = F.canonical_G_generators()["G5"][0]
    inst = sp.expand(g5.xreplace(subs))
    return sp.expand(sp.Poly(inst, y).nth(P.PHI_STRIPPED_DEGREE))


def check_samples(verbose, lc_phi, n_kill=5, n_alive=4):
    """Link the closed-form criterion to the exact symbolic computation."""
    data = json.load(open(os.path.join(HERE, P.UNIVERSE["sub2"]), encoding="utf-8"))
    kills, alives = [], []
    for case in data["cases"]:
        for state in case["states"]:
            res = P.classify_state(state, case, "sub2")
            bucket = kills if res["killed"] else alives
            if len(bucket) < 200:
                bucket.append((state, case, res))
    # spread the samples across the list rather than taking a prefix
    def sample(rows, n):
        if not rows:
            return []
        step = max(1, len(rows) // n)
        return rows[::step][:n]

    ok = True
    for state, case, res in sample(kills, n_kill):
        top = _g5_top_coefficient(_generic_state(state, case), "sub2")
        good = sp.expand(top - lc_phi) == 0
        ok &= check("KILL sample a_t=%s degs=%s" % (case["a_t"], res["degrees"]), good,
                    "deg-34 coeff = %s, expected lc(Phi) = %s" % (str(top)[:60], lc_phi))
        if verbose:
            print("     [kill ] a_t=%-3s d0=%-4s d1=%-4s d2=%-4s e=%-4s  max=%-3s -> deg34 coeff %s  %s"
                  % (case["a_t"], res["degrees"]["d0"], res["degrees"]["d1"],
                     res["degrees"]["d2"], res["degrees"]["dm1"],
                     res["max_g5body_degree"], top, "OK" if good else "FAIL"))

    for state, case, res in sample(alives, n_alive):
        top = _g5_top_coefficient(_generic_state(state, case), "sub2")
        # not killed => G5body reaches degree 34, so the coefficient must carry
        # unknowns/parameters beyond the bare Phi constant
        good = sp.expand(top - lc_phi) != 0
        ok &= check("NON-KILL sample a_t=%s" % case["a_t"], good,
                    "deg-34 coeff collapsed to lc(Phi) though criterion did not fire")
        if verbose:
            print("     [alive] a_t=%-3s max=%-3s -> deg34 coeff has %d free symbols  %s"
                  % (case["a_t"], res["max_g5body_degree"],
                     len(sp.expand(top).free_symbols), "OK" if good else "FAIL"))
    if verbose:
        print("  5/6. sample linkage (%d kill + %d non-kill)  %s"
              % (n_kill, n_alive, "OK" if ok else "FAIL"))
    return ok


def check_rollup(verbose):
    ok = True
    for window, fn in P.UNIVERSE.items():
        data = json.load(open(os.path.join(HERE, fn), encoding="utf-8"))
        total = sum(len(c["states"]) for c in data["cases"])
        declared = data.get("state_total")
        ok &= check("rollup %s" % window, declared is None or total == declared,
                    "counted %d, declared %s" % (total, declared))
        if verbose:
            print("  7. rollup %-5s counted %-6d declared %-6s  %s"
                  % (window, total, declared, "OK" if declared == total else "FAIL"))
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()
    verbose = not args.quiet

    if verbose:
        print("=" * 78)
        print("PHI-DEPTH CRITERION -- INDEPENDENT VERIFIER")
        print("=" * 78)

    check_strip_consistency(verbose)
    _, lc_phi = check_phi(verbose)
    check_caps(verbose)
    check_canonical_g5(verbose)
    check_samples(verbose, lc_phi)
    check_rollup(verbose)

    if FAILURES:
        print("\nFAILURES (%d of %d checks):" % (len(FAILURES), CHECKS[0]))
        for f in FAILURES:
            print("  -", f)
        return 1
    print("\nALL %d PHI-DEPTH CHECKS PASSED" % CHECKS[0])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
