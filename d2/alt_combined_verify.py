#!/usr/bin/env python3
"""Independent spot-checker for the combination layer (alt_combined.py /
alt_combined.json).

This file NEVER imports alt_combined.  It re-transcribes the finite-place cones
directly from ALT_REGIME.md / ALT_REGIME_L2.md and:

  PART 0  independently recomputes the per-branch kill/remain counts from the
          sweep survivors and asserts they equal alt_combined.json exactly
          (a full independent re-derivation of the engine, not just spot checks);
  PART 1-3  re-derive THREE explicit state kills by hand-style valuation chains
          (two required by the task; a third exhibits the sigma-coupling);
  PART 4  verify the whole-branch verdict: NO branch is whole-killed, shown by
          exhibiting an explicit surviving state in the tightest branch and
          recomputing its finite-place order sum by hand.

Every assertion prints PASS; any failure raises.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
NEG = None  # sweep/JSON encode identically-zero degree as null -> None


def banner(t):
    print("\n" + "=" * 70 + f"\n{t}\n" + "=" * 70)


# ---------------------------------------------------------------------------
# Independent transcription of the finite-place cones (source-cited inline).
# These are re-typed from the docs, NOT imported from alt_combined.py.
# ---------------------------------------------------------------------------
def T1_t(a):
    # ALT_REGIME_L2.md sec.2, T1 t-row (h6/h5 order cone).
    if a == 11:
        return [(x, z) for x in range(5, 10) for z in range(3, 13)]
    if a == 12:
        return [(3, 0), (4, 1), (5, 2), (6, 3), (7, 4), (8, 5)] + \
               [(9, z) for z in range(6, 13)]
    if a == 13:
        return []
    if a == 14:
        return [(6, 0), (7, 1), (8, 2), (9, 3)]
    raise ValueError(a)


def T1_q(b):
    # ALT_REGIME_L2.md sec.2, T1 q-row.
    if b == 1:
        return [(1, 0), (2, 1)] + [(x, z) for x in range(3, 10)
                                   for z in range(2, 13)]
    if b == 2:
        return [(7, z) for z in range(5, 13)]
    if b == 3:
        return [(4, 0), (5, 1), (6, 2), (7, 3), (8, 4), (9, 5)]
    if b == 4:
        return []
    raise ValueError(b)


def T1_t_s0(a):
    return {11: list(range(5, 10)), 12: [9], 13: [], 14: []}[a]


def T1_q_s0(b):
    return {1: list(range(3, 10)), 2: [7], 3: [], 4: []}[b]


def T2_zmin_t(a):
    return 3 * a - 30                # v_t(sigma) >= w (ALT_REGIME.md T2)


def T2_q(b):
    # ALT_REGIME.md T2 q-lemmas: b=1 z>=2; b=2 impossible; b=3 z=7; b=4 imposs.
    return {1: 2, 2: None, 3: 7, 4: None}[b]


def t1_feasible(a, bvec, d1, sig):
    active = [b for b in bvec if b > 0]
    if sig is NEG:                                   # sigma == 0
        xs = T1_t_s0(a)
        if not xs:
            return False
        m = min(xs)
        for b in active:
            q = T1_q_s0(b)
            if not q:
                return False
            m += min(q)
        return m <= d1
    reach = {(0, 0)}
    for cone in [T1_t(a)] + [T1_q(b) for b in active]:
        if not cone:
            return False
        nxt = set()
        for (X, Z) in reach:
            for (x, z) in cone:
                if X + x <= 9 and Z + z <= 12:
                    nxt.add((X + x, Z + z))
        reach = nxt
        if not reach:
            return False
    return any(X <= d1 and Z <= sig for (X, Z) in reach)


def t2_feasible(a, bvec, sig):
    active = [b for b in bvec if b > 0]
    z = T2_zmin_t(a)
    for b in active:
        q = T2_q(b)
        if q is None:
            return False
        z += q
    if z > 12:
        return False
    return z <= sig


def feasible(branch, a, bvec, d1, sig):
    return t1_feasible(a, bvec, d1, sig) if branch == "T1" \
        else t2_feasible(a, bvec, sig)


# ---------------------------------------------------------------------------
sweep = json.loads((ROOT / "alt_inf_sweep.json").read_text(encoding="utf-8"))
comb = json.loads((ROOT / "alt_combined.json").read_text(encoding="utf-8"))
sweep_by = {b["id"]: b for b in sweep["branches"]}
comb_by = {b["id"]: b for b in comb["branches"]}

fails = 0


def check(cond, msg):
    global fails
    if cond:
        print(f"  PASS  {msg}")
    else:
        fails += 1
        print(f"  FAIL  {msg}")


# ---------------------------------------------------------------------------
banner("PART 0  independent recomputation of ALL per-branch counts")
tot_k = tot_r = 0
for bid, sb in sweep_by.items():
    a, bvec, branch = sb["a"], tuple(sb["b"]), sb["branch"]
    killed = remain = 0
    for (d2, d1, sig, e, nobl) in sb["surviving_states_compact"]:
        if feasible(branch, a, bvec, d1, sig):
            remain += 1
        else:
            killed += 1
    cb = comb_by[bid]
    tot_k += killed
    tot_r += remain
    ok = (killed == cb["states_killed"] and remain == cb["states_remaining"])
    if not ok:
        print(f"  MISMATCH {bid}: indep=({killed},{remain}) "
              f"json=({cb['states_killed']},{cb['states_remaining']})")
    check(ok, f"{bid}: killed={killed} remain={remain} matches json")
check(tot_k == comb["summary"]["states_killed_by_intersection"],
      f"total killed {tot_k} == summary "
      f"{comb['summary']['states_killed_by_intersection']}")
check(tot_r == comb["summary"]["states_remaining"],
      f"total remaining {tot_r} == summary "
      f"{comb['summary']['states_remaining']}")


def assert_survivor(bid, d2, d1, sig, e):
    """The state must be a SWEEP survivor (input to the intersection)."""
    rows = sweep_by[bid]["surviving_states_compact"]
    return any(r[0] == d2 and r[1] == d1 and r[2] == sig and r[3] == e
               for r in rows)


def assert_killed(bid, d2, d1, sig, e):
    ks = comb_by[bid]["killed_states"]
    for k in ks:
        s = k["state"]
        if (s["deg_d2"], s["deg_d1"], s["deg_sigma"], s["deg_e"]) == \
           (d2, d1, sig, e):
            return k
    return None


# ---------------------------------------------------------------------------
banner("PART 1  hand-derive kill A (T1 deg d1 too small)")
# a11_b0000_T1, state (deg d2,d1,sigma,e)=(5,2,10,11).
# ALT_REGIME_L2.md sec.2 T1 t-row, a=11:  every admissible (v_t(d1),v_t(sigma))
# has v_t(d1) in [5,9].  b=0000 => t is the only active place, so
#   deg d1 >= v_t(d1) >= 5   (ALT_REGIME.md: orders add toward degree).
# The state has deg d1 = 2 < 5  =>  no admissible d1 exists  =>  KILLED.
bid = "a11_b0000_T1"
check(assert_survivor(bid, 5, 2, 10, 11),
      f"{bid} (5,2,10,11) is a degree-sweep survivor (intersection input)")
minx = min(x for (x, z) in T1_t(11))
check(minx == 5, f"a=11 T1 t-cone min v_t(d1) = {minx} (=5, ALT_REGIME_L2 sec.2)")
check(2 < minx, f"state deg d1 = 2 < {minx} = forced min sum v_P(d1)")
k = assert_killed(bid, 5, 2, 10, 11)
check(k is not None and "min sum v_P(d1)=5" in k["reason"],
      f"alt_combined.json kills it with the d1-order reason: "
      f"{k['reason'] if k else None}")

# ---------------------------------------------------------------------------
banner("PART 2  hand-derive kill B (T2 deg sigma too small)")
# a14_b0000_T2, state deg sigma = 5 (e.g. (6,-,5,15)).
# T2 at t (ALT_REGIME.md 'Terminal plus first-level local lemmas'):
#   v_t(sigma) >= w = 3a-30 = 3*14-30 = 12.   b=0000 => t only, so
#   deg sigma >= v_t(sigma) >= 12.  State deg sigma = 5 < 12  =>  KILLED.
bid = "a14_b0000_T2"
check(assert_survivor(bid, 6, NEG, 5, 15),
      f"{bid} (6,-,5,15) is a degree-sweep survivor")
w = 3 * 14 - 30
check(w == 12 and T2_zmin_t(14) == 12,
      f"a=14 T2: v_t(sigma) >= w = {w} (=12)")
check(5 < 12, "state deg sigma = 5 < 12 = forced min sum v_P(sigma)")
k = assert_killed(bid, 6, NEG, 5, 15)
check(k is not None and "min sum v_P(sigma)=12" in k["reason"],
      f"alt_combined.json kills it with the sigma-order reason: "
      f"{k['reason'] if k else None}")

# ---------------------------------------------------------------------------
banner("PART 3  hand-derive kill C (T1 sigma-coupling in the a=11 rectangle)")
# a11_b0000_T1, state (5,6,2,11): deg d1 = 6 (>=5, so d1 alone is fine) but
# the a=11 T1 t-cone is the RECTANGLE 5<=x<=9, 3<=z<=12 -- it has NO pair with
# z<3, i.e. every admissible counterexample has v_t(sigma) >= 3, hence
#   deg sigma >= v_t(sigma) >= 3.   State deg sigma = 2 < 3  =>  KILLED.
bid = "a11_b0000_T1"
check(assert_survivor(bid, 5, 6, 2, 11),
      f"{bid} (5,6,2,11) is a degree-sweep survivor")
minz = min(z for (x, z) in T1_t(11))
check(minz == 3, f"a=11 T1 t-cone min v_t(sigma) = {minz} (=3, rectangle floor)")
check(2 < minz, "state deg sigma = 2 < 3 = forced min sum v_P(sigma)")
k = assert_killed(bid, 5, 6, 2, 11)
check(k is not None and "sum v_P(sigma)=3" in k["reason"],
      f"alt_combined.json kills it with the sigma-coupling reason: "
      f"{k['reason'] if k else None}")

# ---------------------------------------------------------------------------
banner("PART 4  whole-branch verdict: NO branch is whole-killed")
# The task asks to re-derive a whole-branch kill IF ONE EXISTS.  None does.
# We verify this honestly: (i) the summary reports 0 whole kills, and (ii) we
# exhibit an explicit SURVIVOR of the tightest branch and recompute its order
# sum by hand -- the branch therefore cannot be whole-killed.
check(comb["summary"]["branches_whole_killed"] == 0,
      "summary: branches_whole_killed == 0 (honest negative result)")
check(all(b["verdict"] == "OPEN" for b in comb["branches"]),
      "every branch verdict is OPEN (>=1 state survives the intersection)")

# tightest branch a11_b3100_T2: Zmin = w(t) + z(b=3) + z(b=1) = 3 + 7 + 2 = 12,
# and deg sigma <= 12, so exactly deg sigma = 12 survives.
bid = "a11_b3100_T2"
zmin = T2_zmin_t(11) + T2_q(3) + T2_q(1)
check(zmin == 12, f"{bid} hand Zmin = 3+7+2 = {zmin} (=12, at the cap)")
check(assert_survivor(bid, NEG, NEG, 12, 15),
      f"{bid} (-,-,12,15) is a degree-sweep survivor")
check(t2_feasible(11, (3, 1, 0, 0), 12),
      f"{bid} deg sigma=12 is finite-place FEASIBLE (12 <= Zmin=12) -> survives")
check(comb_by[bid]["states_remaining"] >= 1,
      f"{bid} retains {comb_by[bid]['states_remaining']} states -> not whole-killed")

# ---------------------------------------------------------------------------
banner("RESULT")
if fails == 0:
    print("ALL CHECKS PASS")
else:
    raise SystemExit(f"{fails} CHECK(S) FAILED")
