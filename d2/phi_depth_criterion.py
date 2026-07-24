#!/usr/bin/env python3
"""phi_depth_criterion.py  --  MILESTONE 1.5b: the Phi-window-depth kill criterion

A CLOSED-FORM, O(1)-per-state specialisation of the bigraded face detector
(`face_kill_sweep.py`), obtained by noticing that the detector's kill at the
`delta(max)` / `nu(max)` face is a pure DEGREE argument and therefore needs no
symbolic expansion at all.

THE CRITERION.  `G5 = 2*Phi + G5body` is weighted-homogeneous of u-weight 17, so
every one of its terms strips by exactly y^(12*17) = y^204 and the stripped G5
lives at y-degree <= 2*17 = 34 (sub2 envelope).  The stripped `Phi` is the FIXED
forced polynomial

    Phi~ = -(y+1)^30 * (2048y^4 - 512y^3 + 320y^2 - 240y + 195) / 6630

of y-degree EXACTLY 34, with leading coefficient -1024/3315 != 0
(`FULL_SYSTEM_BRIDGE.md:50,166`; `c = -1/6630` is the forced ODE constant, not a
free parameter -- and `ALT_BRIDGE.md:50` records that the alt regime strips to the
identical object).

NORMALISATION NOTE.  The canonical generator is `G5 = G5body + Phi`
(`full_system_bridge.py:107`; the C11 membership certificate in
`f37_sat_verify.py` verifies `f31 == ... + c4*(G5body + Phi)`).
`FULL_SYSTEM_BRIDGE.md:62` states `G5 = 2*Phi + G5body`, contradicting line 114
of the same file; line 62 is erroneous and `bigrade_annotator.py:675` inherited
it.  This criterion is INSENSITIVE to that discrepancy -- it uses only
`deg Phi~ = 34` and `lc(Phi~) != 0` -- but the emitted certificate value is
-1024/3315 canonically, not -2048/3315.

Hence: if for a given state EVERY term of

    G5body = -3 d0 dm1 dm4 - 3 d0 dm2 dm3 - 3 d1 dm2 dm4 - (3/2) d1 dm3^2
             - 3 d2 dm3 dm4 - (3/2) dm1^2 dm3 - (3/2) dm1 dm2^2

has stripped y-degree STRICTLY BELOW 34 -- even with all three spare series taken
at their maximal admissible degree -- then the degree-34 coefficient of G5 is
`2*lc(Phi~) = -2048/3315`, and the equation `-2048/3315 = 0` is unsatisfiable.
The state admits NO spares: KILL.

SOUNDNESS.  The test is a SUFFICIENT condition, deliberately one-sided:

  * spares are taken at their CAP degrees (the most generous case), so a state
    that fails the test is simply not decided here -- no survivorship is claimed;
  * cancellation among G5body terms is irrelevant, because a kill requires that
    NO term reaches degree 34 in the first place;
  * state degrees (`deg_d1`, `deg_d2`, `deg_e`, `deg_sigma`) are degree-EXACT by
    the ansatz semantics (leading coefficients declared nonzero), which the suite
    enforces, so `deg d0 = max(2 deg_d2, deg_sigma)` is exact too.

Caps: stripped `deg d_{4-k} <= (14-12)k = 2k` (sub2) and `<= (15-12)k = 3k`
(sub1), from envelope node C6 (`deg <= 14w`/`15w`, `ord >= 12w`), giving spare
caps (dm2,dm3,dm4) = (12,14,16) for sub2 and (18,21,24) for sub1.

Read-only.  Writes `phi_depth_criterion.json` only.

Usage:
    python phi_depth_criterion.py            # sweep both state universes
    python phi_depth_criterion.py --quiet    # census only
"""
from __future__ import annotations

import argparse
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__)) or "."
OUT_PATH = os.path.join(HERE, "phi_depth_criterion.json")

NEG_INF = float("-inf")

PHI_STRIPPED_DEGREE = 34          # deg Phi~ ; equals the sub2 cap 2*17 exactly
WINDOW_CAP_SLOPE = {"sub2": 2, "sub1": 3}
UNIVERSE = {"sub2": "phase_d_states_sub2.json",
            "sub1": "phase_d_states_sub1.json"}

# (name, k-multiset) for the seven G5body terms; each has sum(k) == 17.
G5BODY_TERMS = (
    ("d0*dm1*dm4", ("d0", "dm1", "dm4")),
    ("d0*dm2*dm3", ("d0", "dm2", "dm3")),
    ("d1*dm2*dm4", ("d1", "dm2", "dm4")),
    ("d1*dm3^2", ("d1", "dm3", "dm3")),
    ("d2*dm3*dm4", ("d2", "dm3", "dm4")),
    ("dm1^2*dm3", ("dm1", "dm1", "dm3")),
    ("dm1*dm2^2", ("dm1", "dm2", "dm2")),
)
SYMBOL_K = {"d0": 4, "d1": 3, "d2": 2, "dm1": 5, "dm2": 6, "dm3": 7, "dm4": 8}
SPARES = ("dm2", "dm3", "dm4")


def _deg(value):
    """Normalise a stored degree; absent / '-inf' / None all mean 'zero poly'."""
    if value is None:
        return NEG_INF
    if isinstance(value, str):
        return NEG_INF if value.strip().lstrip("-").lower() in ("inf", "") else int(value)
    return int(value)


def state_degrees(state, case):
    """Stripped degrees of d2, d1, e, sigma, d0 for one phase-D state."""
    deg_d2 = NEG_INF if case.get("d2_zero") else _deg(state.get("deg_d2"))
    deg_d1 = _deg(state.get("deg_d1"))
    deg_e = _deg(state.get("deg_e"))
    deg_sigma = NEG_INF if case.get("sigma_zero") else _deg(state.get("deg_sigma"))
    # d0 = (d2^2 + sigma)/4 -- degree-exact, no cancellation between the two
    # (they are independent coordinates of the state).
    deg_d0 = max(2 * deg_d2 if deg_d2 != NEG_INF else NEG_INF, deg_sigma)
    return {"d0": deg_d0, "d1": deg_d1, "d2": deg_d2, "dm1": deg_e}


def max_g5body_degree(degs, slope):
    """Largest stripped y-degree any G5body term can reach, spares at their caps."""
    caps = {s: slope * SYMBOL_K[s] for s in SPARES}
    best, witness = NEG_INF, None
    for name, factors in G5BODY_TERMS:
        total = 0
        for f in factors:
            d = caps[f] if f in caps else degs[f]
            if d == NEG_INF:
                total = NEG_INF
                break
            total += d
        if total > best:
            best, witness = total, name
    return best, witness


def classify_state(state, case, window):
    slope = WINDOW_CAP_SLOPE[window]
    degs = state_degrees(state, case)
    best, witness = max_g5body_degree(degs, slope)
    killed = best < PHI_STRIPPED_DEGREE
    return {
        "killed": killed,
        "max_g5body_degree": None if best == NEG_INF else best,
        "binding_term": witness,
        "margin": None if best == NEG_INF else PHI_STRIPPED_DEGREE - best,
        "degrees": {k: (None if v == NEG_INF else v) for k, v in degs.items()},
    }


def sweep(window, verbose=True):
    path = os.path.join(HERE, UNIVERSE[window])
    data = json.load(open(path, encoding="utf-8"))
    total = killed = 0
    kill_rows = []
    per_at = {}
    for case in data["cases"]:
        for state in case["states"]:
            total += 1
            res = classify_state(state, case, window)
            at = case.get("a_t")
            bucket = per_at.setdefault(at, [0, 0])
            bucket[0] += 1
            if res["killed"]:
                killed += 1
                bucket[1] += 1
                if len(kill_rows) < 400:      # sample, not the whole list
                    kill_rows.append({
                        "a_t": at, "b": case.get("b"), "branch": case.get("branch"),
                        "d2_zero": case.get("d2_zero"), "sigma_zero": case.get("sigma_zero"),
                        **{k: v for k, v in res.items() if k != "killed"},
                    })
    # fail-loud: the rolled-up total must reconstruct the declared universe size
    declared = data.get("state_total")
    if declared is not None and declared != total:
        raise SystemExit("FAIL: %s rolled-up states %d != declared state_total %d"
                         % (window, total, declared))
    if verbose:
        print("  %-5s  states=%-6d  Phi-depth kills=%-5d  (%.2f%%)   [declared total %s OK]"
              % (window, total, killed, 100.0 * killed / total if total else 0.0, declared))
        for at in sorted(per_at, key=lambda x: (x is None, x)):
            n, k = per_at[at]
            if k:
                print("        a_t=%-3s  %5d states  %5d killed" % (at, n, k))
    return {"window": window, "states": total, "killed": killed,
            "declared_state_total": declared,
            "per_a_t": {str(k): {"states": v[0], "killed": v[1]}
                        for k, v in per_at.items()},
            "kill_sample": kill_rows}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    if not args.quiet:
        print("=" * 78)
        print("PHI-WINDOW-DEPTH KILL CRITERION  (Milestone 1.5b)")
        print("  kill iff max deg(G5body) < deg(Phi~) = %d, spares at caps" % PHI_STRIPPED_DEGREE)
        print("=" * 78)
    out = {"schema": "d2-phi-depth-criterion-v1",
           "phi_stripped_degree": PHI_STRIPPED_DEGREE,
           "window_cap_slope": WINDOW_CAP_SLOPE,
           "windows": {}}
    for window in ("sub2", "sub1"):
        out["windows"][window] = sweep(window, verbose=not args.quiet)

    with open(OUT_PATH, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2, sort_keys=True)
    tot = sum(w["states"] for w in out["windows"].values())
    kil = sum(w["killed"] for w in out["windows"].values())
    print("\nTOTAL: %d / %d states killed by Phi-window-depth (%.2f%%)"
          % (kil, tot, 100.0 * kil / tot if tot else 0.0))
    print("wrote", os.path.basename(OUT_PATH))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
