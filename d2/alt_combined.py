#!/usr/bin/env python3
"""Combination layer: intersect the alternate-regime max-plus DEGREE sweep
(alt_inf_sweep.json) with the FINITE-PLACE valuation constraints proved in
ALT_REGIME.md and ALT_REGIME_L2.md, for the 27 open flipped branches.

Companion: ALT_COMBINED.md.  Verifier: alt_combined_verify.py.

WHY this exists
---------------
alt_inf_sweep.json is the DEGREE (infinity) layer only.  ALT_INF_SWEEP.md
[judgment] J5 states verbatim: "Combining layers (intersecting surviving states
with the residue/valuation constraints) is the next step and is left to the
finite-place engine."  This is that engine.  The two layers were derived
independently and had never been intersected.

WHAT it does
------------
Each surviving degree state of a branch fixes  deg d1  and  deg sigma  (and
deg d2, deg e).  The finite-place lemmas of ALT_REGIME.md / ALT_REGIME_L2.md
constrain, at the t-place and at every active q-place p (b_i>0), the LOCAL
valuations  x = v_P(d1),  z = v_P(sigma)  of any surviving counterexample to lie
in an explicit allowed cone.  Because valuations at distinct linear places add
toward the polynomial degree (ALT_REGIME.md "Orders at distinct places add
toward polynomial degree"; ALT_REGIME_L2.md (R1)/(R2)):

        deg d1    >=  sum_P v_P(d1)  = X,        deg sigma >= sum_P v_P(sigma) = Z,

a degree state is KILLED by the intersection when NO admissible per-place
selection of local cone pairs (x_P, z_P) has  X <= deg d1  AND  Z <= deg sigma
simultaneously (T1), resp.  Z <= deg sigma  (T2).  That is: the finite-place
lower bounds are incompatible with the state's enumerated degrees.

Every cone used is cited to the doc+section that proves it (see CITE below).
No new unproven lemma is introduced -- the cones are transcribed verbatim from
ALT_REGIME_L2.md section 2 (the h6/h5 order tables [C]) and ALT_REGIME.md
"Terminal plus first-level local lemmas".

Sound over-approximation: we kill a state ONLY when even the most favourable
admissible valuation assignment cannot reach the state's degrees, so a kill is a
genuine proof of infeasibility.  A surviving state keeps its infinity-layer
leading-cancellation obligations AND still owes the exact residue congruences
(D_t)/(D_p) through h_5 (ALT_REGIME_L2.md section 5) -- merged below.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent

# ---------------------------------------------------------------------------
# CITATIONS (doc + section proving each cone).  Attached to every kill record.
# ---------------------------------------------------------------------------
CITE = {
    "T1_t":   "ALT_REGIME_L2.md sec.2 [C] T1 t-row (h6/h5 order cone) + "
              "ALT_REGIME.md 'Terminal plus first-level local lemmas' "
              "(T1 anchor 2 v_t(d1)>=w, odd-s parity, h7=8192 d1^2)",
    "T1_q":   "ALT_REGIME_L2.md sec.2 [C] T1 q-row (h6/h5 order cone) + "
              "ALT_REGIME.md T1 q-terminal (3b_i+v_i(g_7)=7+2v_i(d1))",
    "T1_s0":  "ALT_REGIME_L2.md sec.2 [C] T1 'sigma=0' column",
    "T2_t":   "ALT_REGIME.md 'Terminal plus first-level local lemmas' T2 "
              "(v_t(sigma)>=w) + ALT_REGIME_L2.md sec.2 (O2), sigma^2 terminal",
    "T2_q":   "ALT_REGIME.md T2 q-lemmas (b=1:z>=2; b=2 impossible; b=3:z=7; "
              "b=4 impossible) + ALT_REGIME_L2.md sec.2 (O2)",
    "L2_R":   "ALT_REGIME_L2.md sec.5 (R0)/(R1)/(R2) residual normal form "
              "(deg d1>=X, deg sigma>=Z, X<=9, Z<=12, deg F<=15-a-sum b)",
}

# ---------------------------------------------------------------------------
# FINITE-PLACE CONES  (transcribed verbatim from the cited sections).
#
# T1 finite-sigma allowed (x,z) = (v_P(d1), v_P(sigma)) pairs, per place.
#   Source: ALT_REGIME_L2.md section 2, "complete projected T1 possibilities
#   through h_5" table [C].
# ---------------------------------------------------------------------------
def t1_t_cone(a):
    """Allowed (x,z) at the t-place for T1, stratum a.  ALT_REGIME_L2 sec.2."""
    if a == 11:                       # 5<=x<=9, 3<=z<=12   (rectangle)
        return [(x, z) for x in range(5, 10) for z in range(3, 13)]
    if a == 12:                       # (3,0)..(8,5) staircase, or x=9,6<=z<=12
        pairs = [(3, 0), (4, 1), (5, 2), (6, 3), (7, 4), (8, 5)]
        pairs += [(9, z) for z in range(6, 13)]
        return pairs
    if a == 13:                       # none (t-local h6/h5 cone empty)
        return []
    if a == 14:                       # (6,0),(7,1),(8,2),(9,3)
        return [(6, 0), (7, 1), (8, 2), (9, 3)]
    raise ValueError(f"unexpected a={a}")


def t1_q_cone(b):
    """Allowed (x,z) at a q-root with v_p(E)=b>0 for T1.  ALT_REGIME_L2 sec.2."""
    if b == 1:                        # (1,0),(2,1), or 3<=x<=9, 2<=z<=12
        return [(1, 0), (2, 1)] + [(x, z) for x in range(3, 10)
                                   for z in range(2, 13)]
    if b == 2:                        # x=7, 5<=z<=12
        return [(7, z) for z in range(5, 13)]
    if b == 3:                        # (4,0),(5,1),(6,2),(7,3),(8,4),(9,5)
        return [(4, 0), (5, 1), (6, 2), (7, 3), (8, 4), (9, 5)]
    if b == 4:                        # none
        return []
    raise ValueError(f"unexpected b={b}")


def t1_t_sigma0(a):
    """Allowed x when sigma==0 identically, t-place, T1.  ALT_REGIME_L2 sec.2."""
    if a == 11:
        return list(range(5, 10))     # 5<=x<=9
    if a == 12:
        return [9]                    # x=9
    if a == 13:
        return []
    if a == 14:
        return []                     # none
    raise ValueError(f"unexpected a={a}")


def t1_q_sigma0(b):
    """Allowed x when sigma==0 identically, q-root, T1.  ALT_REGIME_L2 sec.2."""
    if b == 1:
        return list(range(3, 10))     # 3<=x<=9
    if b == 2:
        return [7]
    if b == 3:
        return []                     # none
    if b == 4:
        return []
    raise ValueError(f"unexpected b={b}")


# T2 (d1==0): only sigma matters.  z = v_P(sigma).
#   t: v_t(sigma) >= w = 3a-30 (i.e. >=3,6,9,12 for a=11,12,13,14).
#   q: b=1 -> z>=2 ; b=2 -> impossible ; b=3 -> z=7 ; b=4 -> impossible.
#   Source: ALT_REGIME.md "Terminal plus first-level local lemmas" T2 +
#           ALT_REGIME_L2.md section 2 (O2).
def t2_t_zmin(a):
    return 3 * a - 30                 # = w


def t2_q_zset(b):
    """Allowed z at a q-root for T2 (returns None if the place is IMPOSSIBLE)."""
    if b == 1:
        return ("min", 2)            # z >= 2
    if b == 2:
        return None                  # impossible: kills the whole branch
    if b == 3:
        return ("exact", 7)          # z = 7 (with residue cancellation)
    if b == 4:
        return None                  # impossible
    raise ValueError(f"unexpected b={b}")


# ---------------------------------------------------------------------------
# Feasibility of a state under the finite-place cones.
# ---------------------------------------------------------------------------
CAP_X, CAP_Z = 9, 12                  # deg d1<=9, deg sigma<=12 (sub1 window)


def _reach_pairs(place_cones):
    """Reachable (X,Z) partial sums over the places, pruned to X<=9, Z<=12.

    place_cones: list of lists of (x,z) allowed pairs, one per active place.
    Returns the frozenset of reachable (sumX, sumZ)."""
    reach = {(0, 0)}
    for cone in place_cones:
        nxt = set()
        for (X, Z) in reach:
            for (x, z) in cone:
                nX, nZ = X + x, Z + z
                if nX <= CAP_X and nZ <= CAP_Z:
                    nxt.add((nX, nZ))
        reach = nxt
        if not reach:
            return frozenset()
    return frozenset(reach)


def t1_feasible(a, bvec, deg_d1, deg_sigma):
    """Return (feasible: bool, witness or reason-tag, cite_key).

    deg_d1 is an int (T1 => d1 != 0).  deg_sigma is int or None (None=sigma==0).
    Feasible iff some admissible per-place (x,z) selection has sum_x<=deg_d1 and
    sum_z<=deg_sigma (finite sigma), resp. sum_x<=deg_d1 over the sigma=0
    columns."""
    active_b = [b for b in bvec if b > 0]

    if deg_sigma is None:                                   # sigma == 0 branch
        # every involved place must permit sigma==0, then need sum_x<=deg_d1
        xs_t = t1_t_sigma0(a)
        if not xs_t:
            return (False, "sigma0_forbidden_at_t", "T1_s0")
        min_x = min(xs_t)
        for b in active_b:
            xs_q = t1_q_sigma0(b)
            if not xs_q:
                return (False, f"sigma0_forbidden_at_q(b={b})", "T1_s0")
            min_x += min(xs_q)
        if min_x <= deg_d1:
            return (True, {"kind": "sigma0", "Xmin": min_x}, "T1_s0")
        return (False, f"sigma0 needs deg d1>={min_x} but deg d1={deg_d1}",
                "T1_s0")

    # finite sigma: joint (X,Z) feasibility
    cones = [t1_t_cone(a)]
    cite = "T1_t"
    for b in active_b:
        c = t1_q_cone(b)
        if not c:
            return (False, f"q-place b={b} cone empty", "T1_q")
        cones.append(c)
    reach = _reach_pairs(cones)
    if not reach:
        return (False, "no admissible (X,Z) within caps", "T1_t")
    ok = [(X, Z) for (X, Z) in reach if X <= deg_d1 and Z <= deg_sigma]
    if ok:
        Xmin = min(X for (X, Z) in ok)
        wit = min(ok, key=lambda p: (p[0], p[1]))
        return (True, {"kind": "finite", "X": wit[0], "Z": wit[1]}, cite)
    # infeasible: report the binding bound
    Xmin_all = min(X for (X, Z) in reach)
    Zmin_all = min(Z for (X, Z) in reach)
    # minimal Z achievable given X<=deg_d1 (if any), else X is the binding one
    xok = [(X, Z) for (X, Z) in reach if X <= deg_d1]
    if not xok:
        return (False, f"deg d1={deg_d1} < min sum v_P(d1)={Xmin_all}", "T1_t")
    Zneed = min(Z for (X, Z) in xok)
    return (False,
            f"with deg d1={deg_d1} the min sum v_P(sigma)={Zneed} > "
            f"deg sigma={deg_sigma}", "T1_t")


def t2_feasible(a, bvec, deg_sigma):
    """Return (feasible, witness/reason, cite_key).  deg_sigma is an int."""
    active_b = [b for b in bvec if b > 0]
    zmin = t2_t_zmin(a)
    detail = [("t", zmin)]
    for b in active_b:
        zs = t2_q_zset(b)
        if zs is None:
            return (False, f"q-place b={b} impossible (T2)", "T2_q")
        zmin += zs[1]
        detail.append((f"q(b={b})", zs[1]))
    if zmin > CAP_Z:
        # even the minimum forced Z exceeds the window cap deg sigma<=12
        return (False, f"forced sum v_P(sigma)={zmin} > cap 12", "T2_t")
    if zmin <= deg_sigma:
        return (True, {"kind": "T2", "Zmin": zmin, "detail": detail}, "T2_t")
    return (False, f"deg sigma={deg_sigma} < min sum v_P(sigma)={zmin}",
            "T2_t" if active_b == [] else "T2_q")


# ---------------------------------------------------------------------------
# Main intersection over all branches.
# ---------------------------------------------------------------------------
def intersect_branch(br, inf_samples):
    a = br["a"]
    bvec = tuple(br["b"])
    branch = br["branch"]
    rows = br["surviving_states_compact"]     # [d2,d1,sig,e,n_obl]
    n_before = len(rows)

    killed = []            # {state, reason, cite, cite_key}
    remaining = []         # {state, inf_n_obl, fp_witness, cite_key}
    kill_hist = {}

    for row in rows:
        d2, d1, sig, e, nobl = row
        if branch == "T1":
            feas, info, ckey = t1_feasible(a, bvec, d1, sig)
        else:
            feas, info, ckey = t2_feasible(a, bvec, sig)
        state = {"deg_d2": d2, "deg_d1": d1, "deg_sigma": sig, "deg_e": e,
                 "deg_E": e - a}
        if feas:
            remaining.append({"state": state, "inf_n_obligations": nobl,
                              "finite_place_witness": info,
                              "cite": CITE[ckey]})
        else:
            killed.append({"state": state, "reason": info,
                           "constraint": CITE[ckey], "cite_key": ckey})
            kill_hist[ckey] = kill_hist.get(ckey, 0) + 1

    verdict = "WHOLE_BRANCH_KILL" if remaining == [] else "OPEN"
    label = ("ENGINE+LEMMA-PROVEN, PENDING AUDIT"
             if verdict == "WHOLE_BRANCH_KILL" else "OPEN")
    return {
        "id": br["id"], "a": a, "b": list(bvec), "sum_b": sum(bvec),
        "branch": branch, "w": br["w"], "deg_E_range": br["deg_E_range"],
        "states_before": n_before,
        "states_killed": len(killed),
        "states_remaining": len(remaining),
        "verdict": verdict, "label": label,
        "kill_constraint_histogram": {CITE[k]: v for k, v in kill_hist.items()},
        "killed_by_cite_key": kill_hist,
        "killed_states": killed,
        "remaining_states": remaining,
    }


def main():
    sweep = json.loads((ROOT / "alt_inf_sweep.json").read_text(encoding="utf-8"))
    branches = sweep["branches"]

    results = []
    for br in branches:
        results.append(intersect_branch(br, br.get("survive_samples", [])))

    tot_before = sum(r["states_before"] for r in results)
    tot_killed = sum(r["states_killed"] for r in results)
    tot_remain = sum(r["states_remaining"] for r in results)
    whole_kills = [r["id"] for r in results
                   if r["verdict"] == "WHOLE_BRANCH_KILL"]

    # aggregate which constraint did the most work
    agg = {}
    for r in results:
        for k, v in r["killed_by_cite_key"].items():
            agg[k] = agg.get(k, 0) + v

    out = {
        "schema": {
            "version": 1,
            "description": "Combination layer: intersection of the alternate-"
                           "regime degree sweep (alt_inf_sweep.json) with the "
                           "finite-place valuation lemmas of ALT_REGIME.md and "
                           "ALT_REGIME_L2.md, over the 27 open flipped branches.",
            "inputs": ["alt_inf_sweep.json (degree layer)",
                       "ALT_REGIME.md (terminal + first-level parity lemmas)",
                       "ALT_REGIME_L2.md (h6/h5 order cones sec.2, residual "
                       "normal form sec.5)"],
            "mechanism": "deg d1 >= sum_P v_P(d1) = X, deg sigma >= sum_P "
                         "v_P(sigma) = Z (orders at distinct places add); a "
                         "state dies iff no admissible per-place cone selection "
                         "meets X<=deg d1 and Z<=deg sigma (T1), resp. Z<=deg "
                         "sigma (T2).",
            "citations": CITE,
        },
        "summary": {
            "n_branches": len(results),
            "branches_whole_killed": len(whole_kills),
            "whole_branch_kills": whole_kills,
            "states_before": tot_before,
            "states_killed_by_intersection": tot_killed,
            "states_remaining": tot_remain,
            "kill_constraint_totals": {CITE[k]: v for k, v in
                                       sorted(agg.items(),
                                              key=lambda kv: -kv[1])},
            "kill_constraint_totals_by_key": dict(sorted(agg.items(),
                                                  key=lambda kv: -kv[1])),
        },
        "branches": results,
    }
    (ROOT / "alt_combined.json").write_text(json.dumps(out, indent=1),
                                            encoding="utf-8")

    # ---- report ----
    print(f"Combination layer -- 27 branches\n")
    hdr = f"{'branch id':16} {'br':>3} {'before':>6} {'killed':>6} " \
          f"{'remain':>6}  verdict"
    print(hdr); print("-" * len(hdr))
    for r in results:
        print(f"{r['id']:16} {r['branch']:>3} {r['states_before']:>6} "
              f"{r['states_killed']:>6} {r['states_remaining']:>6}  "
              f"{r['verdict']}")
    print("-" * len(hdr))
    print(f"states: {tot_before} -> {tot_remain}  "
          f"(killed {tot_killed} by intersection)")
    print(f"whole-branch kills: {len(whole_kills)}  {whole_kills}")
    print("\nconstraint work census (kills per cited constraint):")
    for k, v in sorted(agg.items(), key=lambda kv: -kv[1]):
        print(f"  {v:5}  {CITE[k]}")
    print("\nWrote alt_combined.json")


if __name__ == "__main__":
    main()
