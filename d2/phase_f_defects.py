#!/usr/bin/env python3
"""Phase F work item F1 -- divisor-defect histograms.

For every surviving degree state (sub2, sub1, alternate regime) compute the
divisor defect vector

    delta_p = deg p - ( max total finite-place valuation sum of p
                        over admissible joins )

for p in {d1, sigma, d2, e, g_level}.  Finite places S = the four q-roots plus
t = y+1.  The maximum finite-place valuation sum is taken over admissible joins
(one Pareto profile per place, every dimension's summed valuation <= the
state's degree in that dimension), mirroring cascade_engine.join_places.

sub2 + sub1: EXACT.  The per-place Pareto profile sets are tiny (<=7 profiles
per place, <=96 complete joins per flag case), so the full join set is
enumerated and, per state, filtered to the joins that fit the state's degree
caps; the per-dimension maximum over that admissible set is exact.

alternate regime: the alt artifact carries, per remaining state, the forced
finite-place valuation sums X = sum_P v_P(d1) and Z = sum_P v_P(sigma) (the
tight cone minima from ALT_REGIME / ALT_REGIME_L2).  There delta_d1 = deg d1 - X
and delta_sigma = deg sigma - Z.  Because X, Z are the forced *minimum*
admissible sums, these alt deltas are UPPER bounds on the true divisor defect
(true defect uses the max sum >= X); for the terminal-pinned dimensions the
forced value is exact, so a delta of 0 is exact ("both attain forced minima").

New file, read-only on all inputs.  Not committed.
"""

from __future__ import annotations

import json
from collections import Counter, defaultdict
from itertools import product
from pathlib import Path

import cascade_engine as ce

ROOT = Path(__file__).resolve().parent
DEPTH = 4
ce.APPLY_RESIDUE_KILLS = True  # match the rl-sweep semantics of phase_d_states


def deg_val(x):
    """Decode a state degree ('-inf'/None -> None meaning identically zero)."""
    return None if x == "-inf" or x is None else int(x)


# ---------------------------------------------------------------------------
# Exact windows (sub2, sub1)
# ---------------------------------------------------------------------------

def case_join_vectors(a, b_vector, branch, sigma_zero, d2_zero, g_zero_levels,
                      config):
    """All complete admissible-place join valuation-sum vectors for a case.

    Returns (levels, join_vectors) where join_vectors is a list of tuples
    (sum_d1, sum_sigma, sum_d2, *sum_g_by_level).  These are the raw
    per-dimension summed valuations of every combination of one Pareto profile
    per finite place (4 q-roots + t); a state later keeps the ones that fit its
    degree caps.
    """
    terminal = ce.T1_TERMINAL if branch == "T1" else ce.T2_TERMINAL
    below = range(terminal - 1, DEPTH - 1, -1)
    g_zero = {lv: (lv in g_zero_levels) for lv in below}
    g_zero[terminal] = False
    levels = tuple(range(terminal, DEPTH - 1, -1))
    r_cap = 10 + 3 * a

    per_place = [
        ce.place_profiles(b, branch, r_cap, DEPTH, sigma_zero, d2_zero,
                          g_zero, config, a)
        for b in b_vector
    ]
    per_place.append(
        ce.t_place_profiles(a, branch, r_cap, DEPTH, sigma_zero, d2_zero,
                            g_zero, config)
    )
    dims = 3 + len(levels)
    vecs_per_place = [
        [p.budget_vector(levels) for p in opts] for opts in per_place
    ]
    joins = []
    for combo in product(*vecs_per_place):
        s = tuple(sum(combo[pl][i] for pl in range(len(combo)))
                  for i in range(dims))
        joins.append(s)
    return levels, joins


def state_defects(state, levels, joins, a, sum_b, branch, sigma_zero, d2_zero,
                  g_zero_levels):
    """Exact defect vector for one state.  Returns dict {name: delta}."""
    deg_d1 = deg_val(state["deg_d1"])
    deg_sigma = deg_val(state["deg_sigma"])
    deg_d2 = deg_val(state["deg_d2"])
    deg_e = int(state["deg_e"])
    deg_g = {int(k): deg_val(v) for k, v in state["deg_g"].items()}

    dims = 3 + len(levels)
    # Degree cap per join dimension; None -> unconstrained (identically-zero
    # polynomial, valuation component is 0 anyway).
    caps = [deg_d1 if branch != "T2" else None,
            None if sigma_zero else deg_sigma,
            None if d2_zero else deg_d2]
    for lv in levels:
        caps.append(None if lv in g_zero_levels else deg_g.get(lv))

    # Admissible joins: fit every constrained dimension.
    admissible = [
        v for v in joins
        if all(caps[i] is None or v[i] <= caps[i] for i in range(dims))
    ]
    if not admissible:
        return None  # should not happen for a surviving state
    maxsum = [max(v[i] for v in admissible) for i in range(dims)]

    out = {}
    names = ["d1", "sigma", "d2"] + [f"g{lv}" for lv in levels]
    for i, name in enumerate(names):
        if caps[i] is None:
            continue
        out[name] = caps[i] - maxsum[i]
    # e: finite-place valuation sum is fixed at a + sum(b_i).
    out["e"] = deg_e - (a + sum_b)
    return out


def run_window(window):
    config = ce.CONFIGS[window]
    data = json.loads((ROOT / f"phase_d_states_{window}.json").read_text())
    cases_out = []
    per_state_records = []  # (window, cellid, deltas)
    for case in data["cases"]:
        a = case["a_t"]
        b_vector = tuple(case["b"])
        branch = case["branch"]
        sigma_zero = case["sigma_zero"]
        d2_zero = case["d2_zero"]
        g_zero_levels = tuple(case["g_zero_levels"])
        levels, joins = case_join_vectors(
            a, b_vector, branch, sigma_zero, d2_zero, g_zero_levels, config)
        cellid = (f"{window}:a{a}_b{''.join(map(str, b_vector))}_{branch}"
                  f"_sz{int(sigma_zero)}_dz{int(d2_zero)}"
                  f"_gz{'.'.join(map(str, g_zero_levels)) or '-'}")
        state_deltas = []
        for st in case["states"]:
            d = state_defects(st, levels, joins, a, sum(b_vector), branch,
                              sigma_zero, d2_zero, g_zero_levels)
            if d is None:
                d = {"_error": "no_admissible_join"}
            state_deltas.append(d)
            per_state_records.append((cellid, a, branch, d))
        cases_out.append({
            "cellid": cellid, "a": a, "b": list(b_vector), "branch": branch,
            "sigma_zero": sigma_zero, "d2_zero": d2_zero,
            "g_zero_levels": list(g_zero_levels),
            "state_count": len(state_deltas),
            "state_deltas": state_deltas,
        })
    return cases_out, per_state_records


# ---------------------------------------------------------------------------
# Alternate regime
# ---------------------------------------------------------------------------

def run_alt():
    data = json.loads((ROOT / "alt_combined.json").read_text())
    cases_out = []
    per_state_records = []
    for br in data["branches"]:
        a = br["a"]
        b_vector = br["b"]
        branch = br["branch"]
        cellid = f"alt:a{a}_b{''.join(map(str, b_vector))}_{branch}"
        state_deltas = []
        for rs in br["remaining_states"]:
            st = rs["state"]
            fw = rs.get("finite_place_witness", {})
            kind = fw.get("kind")
            d = {}
            deg_d1 = deg_val(st.get("deg_d1"))
            deg_sigma = deg_val(st.get("deg_sigma"))
            # kind 'finite': X=sum v_P(d1), Z=sum v_P(sigma) (both live).
            # kind 'sigma0': sigma==0 flag, Xmin=min sum v_P(d1).
            # kind 'T2'    : d1==0 flag, Zmin=min sum v_P(sigma).
            X = fw.get("X", fw.get("Xmin"))
            Z = fw.get("Z", fw.get("Zmin"))
            if X is not None and deg_d1 is not None:
                d["d1"] = deg_d1 - X
            if Z is not None and deg_sigma is not None:
                d["sigma"] = deg_sigma - Z
            state_deltas.append(d)
            per_state_records.append((cellid, a, branch, d))
        cases_out.append({
            "cellid": cellid, "a": a, "b": b_vector, "branch": branch,
            "state_count": len(state_deltas),
            "state_deltas": state_deltas,
        })
    return cases_out, per_state_records


# ---------------------------------------------------------------------------
# Aggregation / histograms
# ---------------------------------------------------------------------------

CORE_DIMS = ("d1", "d2", "sigma", "e")


def _group_vals(d, branch, group):
    """Values of the deltas in a named dimension group for one state."""
    terminal_g = "g7" if branch == "T1" else "g6"
    if group == "core":
        keys = CORE_DIMS
    elif group == "core_gt":          # core + terminal g
        keys = CORE_DIMS + (terminal_g,)
    else:                              # "full": every computed delta
        return [v for k, v in d.items() if not k.startswith("_")]
    return [d[k] for k in keys if k in d]


def summarize(records, exact):
    """records: list of (cellid, a, branch, deltas).  Returns summary dict."""
    per_delta_hist = defaultdict(Counter)      # name -> Counter(delta)
    groups = ("core", "core_gt", "full")
    maxdelta_hist = {g: Counter() for g in groups}
    n = 0
    n_le2 = {g: 0 for g in groups}
    n_all0 = {g: 0 for g in groups}
    # per-cell / per-branch / per-a tallies keyed on the core_gt group
    by_branch = defaultdict(lambda: [0, 0])    # branch -> [total, le2(core_gt)]
    by_a = defaultdict(lambda: [0, 0])
    by_cell = defaultdict(lambda: [0, 0, 0])   # cell -> [total, le2, n_def0]
    for cellid, a, branch, d in records:
        vals_full = [v for k, v in d.items() if not k.startswith("_")]
        if not vals_full:
            continue
        n += 1
        for k, v in d.items():
            if not k.startswith("_"):
                per_delta_hist[k][v] += 1
        for g in groups:
            gv = _group_vals(d, branch, g)
            if not gv:
                continue
            maxdelta_hist[g][max(gv)] += 1
            if all(v <= 2 for v in gv):
                n_le2[g] += 1
            if all(v == 0 for v in gv):
                n_all0[g] += 1
        gv = _group_vals(d, branch, "core_gt")
        le2 = bool(gv) and all(v <= 2 for v in gv)
        by_branch[branch][0] += 1
        by_branch[branch][1] += le2
        by_a[a][0] += 1
        by_a[a][1] += le2
        c = by_cell[cellid]
        c[0] += 1
        c[1] += le2
        c[2] += sum(1 for v in vals_full if v == 0)
    return {
        "n_states": n,
        "frac_all_le2": {g: (n_le2[g] / n if n else None) for g in groups},
        "n_all_le2": n_le2,
        "n_all_zero": n_all0,
        "per_delta_hist": {k: dict(sorted(v.items(), key=lambda x: float(x[0])))
                           for k, v in per_delta_hist.items()},
        "maxdelta_hist": {g: dict(sorted(maxdelta_hist[g].items()))
                          for g in groups},
        "by_branch": {k: {"total": v[0], "all_le2_core_gt": v[1],
                          "frac": v[1] / v[0] if v[0] else None}
                      for k, v in sorted(by_branch.items())},
        "by_a": {str(k): {"total": v[0], "all_le2_core_gt": v[1],
                          "frac": v[1] / v[0] if v[0] else None}
                 for k, v in sorted(by_a.items())},
        "by_cell": {k: {"total": v[0], "all_le2_core_gt": v[1],
                        "n_defect0": v[2]}
                    for k, v in by_cell.items()},
        "exact": exact,
    }


def top_pilot_cells(records, alt_records, k=20):
    """Cell-level pilot rankings for F2.

    Returns three lists:
      by_def0_count   -- cells with the most defect-0 divisor polynomials
                         (core+gT dims; the literal "most defect-0" ranking).
      by_forced_state -- cells with the most FULLY forced states (every
                         core+gT delta = 0): the tightest reconstruction cells.
      alt_both_min    -- alternate-regime cells ranked by the number of states
                         where d1 AND sigma both attain their forced minima
                         (delta_d1 = delta_sigma = 0).
    """
    cell_info = defaultdict(lambda: {"n_def0": 0, "n_states": 0,
                                     "n_states_all0": 0, "branch": None,
                                     "a": None, "alt": False, "alt_both_min": 0})

    def add(cellid, a, branch, d, alt):
        terminal_g = "g7" if branch == "T1" else "g6"
        keys = CORE_DIMS + (terminal_g,)
        gv = [d[kk] for kk in keys if kk in d]
        if not gv:
            return
        ci = cell_info[cellid]
        ci["branch"] = branch
        ci["a"] = a
        ci["alt"] = ci["alt"] or alt
        ci["n_states"] += 1
        ci["n_def0"] += sum(1 for v in gv if v == 0)
        if all(v == 0 for v in gv):
            ci["n_states_all0"] += 1
        if alt and d.get("d1") == 0 and d.get("sigma") == 0:
            ci["alt_both_min"] += 1

    for cellid, a, branch, d in records:
        add(cellid, a, branch, d, alt=False)
    for cellid, a, branch, d in alt_records:
        add(cellid, a, branch, d, alt=True)

    def take(key, k, filt=None):
        items = [(c, i) for c, i in cell_info.items()
                 if filt is None or filt(i)]
        items.sort(key=lambda kv: key(kv[1]), reverse=True)
        return [{"cellid": c, **i} for c, i in items[:k]]

    return {
        "by_def0_count": take(
            lambda i: (i["n_def0"], i["n_states_all0"]), k),
        "by_forced_state": take(
            lambda i: (i["n_states_all0"],
                       i["n_states_all0"] / i["n_states"] if i["n_states"] else 0),
            k),
        "alt_both_min": take(
            lambda i: (i["alt_both_min"], i["n_states_all0"]), k,
            filt=lambda i: i["alt"] and i["alt_both_min"] > 0),
    }


def main():
    sub2_cases, sub2_rec = run_window("sub2")
    sub1_cases, sub1_rec = run_window("sub1")
    alt_cases, alt_rec = run_alt()

    summaries = {
        "sub2": summarize(sub2_rec, exact=True),
        "sub1": summarize(sub1_rec, exact=True),
        "alt": summarize(alt_rec, exact="upper_bound_on_defect"),
        "sub2+sub1": summarize(sub2_rec + sub1_rec, exact=True),
        "all": summarize(sub2_rec + sub1_rec + alt_rec, exact="mixed"),
    }
    pilots = top_pilot_cells(sub2_rec + sub1_rec, alt_rec, k=20)

    payload = {
        "schema": 1,
        "description": "Phase F1 divisor-defect histograms (sub2/sub1 exact, "
                       "alt upper-bound). delta_p = deg p - max finite-place "
                       "valuation sum over admissible joins; S = 4 q-roots + t.",
        "depth": DEPTH,
        "windows": {
            "sub2": {"cases": sub2_cases},
            "sub1": {"cases": sub1_cases},
            "alt": {"cases": alt_cases},
        },
        "summaries": summaries,
        "top_pilot_cells": pilots,
    }
    out = ROOT / "phase_f_defects.json"
    out.write_text(json.dumps(payload, separators=(",", ":")) + "\n",
                   encoding="utf-8")
    size = out.stat().st_size
    print(f"wrote {out.name} ({size/1e6:.1f} MB)")
    for w in ("sub2", "sub1", "alt", "sub2+sub1", "all"):
        s = summaries[w]
        f = s["frac_all_le2"]
        print(f"{w:10s} states={s['n_states']:6d} "
              f"le2[core]={f['core']:.4f} le2[core+gT]={f['core_gt']:.4f} "
              f"le2[full]={f['full']:.4f}")
    write_md(summaries, pilots)
    return payload


def _hist_row(hist, lo=0, hi=8):
    cells = []
    for k in range(lo, hi + 1):
        cells.append(str(hist.get(k, 0)))
    tail = sum(v for kk, v in hist.items() if kk > hi)
    cells.append(str(tail))
    return cells


def write_md(summaries, pilots):
    L = []
    L.append("# Phase F1 -- divisor-defect histograms\n")
    L.append("Generated by `phase_f_defects.py` (new file, not committed). "
             "Read-only over the committed Phase D / cascade / alt artifacts.\n")
    L.append("## Definition\n")
    L.append("For a polynomial p under a surviving state, "
             "`delta_p = deg p - (max total finite-place valuation sum of p "
             "over admissible joins)`, finite places S = {four q-roots, t=y+1}. "
             "The max is over admissible joins (one Pareto profile per place, "
             "every dimension's summed valuation <= the state degree), "
             "mirroring `cascade_engine.join_places`. delta = 0 means p is "
             "forced to `lambda * prod (y-s)^{v_s}` -- determined up to one "
             "scalar.\n")
    L.append("Per-place Pareto sets are tiny (<=7 profiles/place, <=96 joins/"
             "case), so **sub2 and sub1 are EXACT**. The **alt** regime uses "
             "the forced finite-place sums X=sum v_P(d1), Z=sum v_P(sigma) "
             "recorded per remaining state; delta_d1=deg d1 - X, "
             "delta_sigma=deg sigma - Z are UPPER bounds on the true defect "
             "(X,Z are forced minima), exact where the terminal law pins them.\n")
    L.append("**Dimension groups.** `core` = {d1, d2, sigma, e}; "
             "`core+gT` = core + terminal g (g7 for T1, g6 for T2); "
             "`full` = every computed delta including the intermediate g4..g6 "
             "chain. The intermediate g-chain are high-degree (~40) auxiliary "
             "polynomials whose divisors are small, so their defect is "
             "structurally large (10-30) and dominates `full`; they are NOT "
             "the residue-reconstruction targets. The Phase-F-relevant headline "
             "is `core+gT`.\n")

    L.append("## Headline fractions (states with ALL deltas in {0,1,2})\n")
    L.append("| window | states | core | core+gT | full |")
    L.append("|---|---|---|---|---|")
    for w in ("sub2", "sub1", "alt", "sub2+sub1", "all"):
        s = summaries[w]
        f = s["frac_all_le2"]
        L.append(f"| {w} | {s['n_states']} | {f['core']:.3f} | "
                 f"{f['core_gt']:.3f} | {f['full']:.4f} |")
    L.append("")
    L.append("All-zero (every delta = 0, i.e. every divisor fully forced):")
    L.append("")
    L.append("| window | core | core+gT | full |")
    L.append("|---|---|---|---|")
    for w in ("sub2", "sub1", "alt"):
        z = summaries[w]["n_all_zero"]
        L.append(f"| {w} | {z['core']} | {z['core_gt']} | {z['full']} |")
    L.append("")

    L.append("## Per-delta distributions (count of states at each defect)\n")
    for w in ("sub2", "sub1", "alt"):
        s = summaries[w]
        L.append(f"### {w}  (n={s['n_states']})\n")
        L.append("| delta | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | >8 |")
        L.append("|---|" + "---|" * 10)
        order = ["d1", "d2", "sigma", "e", "g7", "g6", "g5", "g4"]
        hists = s["per_delta_hist"]
        for name in order:
            if name not in hists:
                continue
            h = {int(float(k)): v for k, v in hists[name].items()}
            L.append(f"| {name} | " + " | ".join(_hist_row(h)) + " |")
        L.append("")

    L.append("## Max-delta-per-state distribution (core+gT group)\n")
    L.append("| window | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | >8 |")
    L.append("|---|" + "---|" * 10)
    for w in ("sub2", "sub1", "alt"):
        h = {int(k): v for k, v in
             summaries[w]["maxdelta_hist"]["core_gt"].items()}
        L.append(f"| {w} | " + " | ".join(_hist_row(h)) + " |")
    L.append("")

    L.append("## Breakdown by branch (core+gT all<=2 fraction)\n")
    L.append("| window | branch | total | all<=2 | frac |")
    L.append("|---|---|---|---|---|")
    for w in ("sub2", "sub1", "alt"):
        for br, v in summaries[w]["by_branch"].items():
            L.append(f"| {w} | {br} | {v['total']} | "
                     f"{v['all_le2_core_gt']} | {v['frac']:.3f} |")
    L.append("")

    L.append("## Breakdown by a_t (core+gT all<=2 fraction)\n")
    for w in ("sub2", "sub1", "alt"):
        L.append(f"### {w}\n")
        L.append("| a | total | all<=2 | frac |")
        L.append("|---|---|---|---|")
        for a, v in summaries[w]["by_a"].items():
            L.append(f"| {a} | {v['total']} | {v['all_le2_core_gt']} | "
                     f"{v['frac']:.3f} |")
        L.append("")

    def pilot_table(title, note, rows, cols):
        L.append(f"## {title}\n")
        L.append(note + "\n")
        L.append("| rank | " + " | ".join(cols) + " |")
        L.append("|---|" + "---|" * len(cols))
        for i, p in enumerate(rows, 1):
            L.append(f"| {i} | " + " | ".join(str(p[c]) for c in cols) + " |")
        L.append("")

    pilot_table(
        "Top 20 pilot cells (ranked by count of defect-0 divisor polynomials)",
        "Tightest cells = most forced (delta=0) core+gT divisor components; "
        "pilot candidates for F2. `n_def0` counts defect-0 dims across the "
        "cell's states.",
        pilots["by_def0_count"],
        ["cellid", "branch", "a", "n_states", "n_def0", "n_states_all0"])

    pilot_table(
        "Top 20 fully-forced cells (most states with ALL core+gT deltas = 0)",
        "Every reconstruction-relevant divisor (d1,d2,sigma,e,terminal g) "
        "forced up to a scalar -- the strongest F2 starting points.",
        pilots["by_forced_state"],
        ["cellid", "branch", "a", "n_states", "n_states_all0"])

    if pilots["alt_both_min"]:
        pilot_table(
            "Alternate-regime pilot cells (d1 AND sigma both at forced minima)",
            "`alt_both_min` = number of states in the cell where "
            "delta_d1 = delta_sigma = 0 (J1-type: both forced). These are the "
            "F2 alt-regime targets.",
            pilots["alt_both_min"],
            ["cellid", "a", "n_states", "alt_both_min", "n_states_all0"])

    (ROOT / "PHASE_F_DEFECTS.md").write_text("\n".join(L), encoding="utf-8")
    print("wrote PHASE_F_DEFECTS.md")


if __name__ == "__main__":
    main()
