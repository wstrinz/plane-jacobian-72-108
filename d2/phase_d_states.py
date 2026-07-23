#!/usr/bin/env python3
"""Phase D worklist: complete residual degree states per open flag case.

The infinity-layer sweeps (cascade_cones_qt_inf_rl.json and its sub1
analogue) record ONE witness degree state per surviving flag case.  This
generator enumerates, for every surviving (branch, flag-case), the COMPLETE
set of degree assignments (deg d2, deg d1, deg sigma, deg e) admitting both
a consistent infinity chain and a finite-place join — the exact residual
state list that the Stage 3 convolution descent must refute cell by cell.

For the eight open sub2 T2 cells this reproduces the R-tables of
T5_T2_COLUMN.md sections 4-6 (projected to (deg e, deg sigma)); the
cross-check is executed and reported when --window sub2 is used.

Output: phase_d_states_<window>.json.  Same-author layer over the engine
(not an independent audit); the honest content is the state lists.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import cascade_engine as ce

NEG_INF = ce.NEG_INF
ROOT = Path(__file__).resolve().parent


def case_states(
    a: int,
    b_vector: tuple[int, ...],
    branch: str,
    sigma_zero: bool,
    d2_zero: bool,
    g_zero_levels: tuple[int, ...],
    config,
    t2_squeeze: bool,
    depth: int = 4,
) -> list[dict]:
    """All degree states admitting a chain + join for one flag case."""

    r_cap = 10 + 3 * a
    terminal = ce.T1_TERMINAL if branch == "T1" else ce.T2_TERMINAL
    below = range(terminal - 1, depth - 1, -1)
    g_zero = {level: level in g_zero_levels for level in below}
    g_zero[terminal] = False
    aux, g_caps, _ = ce.resolve_caps(branch, r_cap, config, a)
    levels = tuple(range(terminal, depth - 1, -1))

    per_place = [
        ce.place_profiles(
            b, branch, r_cap, depth, sigma_zero, d2_zero, g_zero, config, a
        )
        for b in b_vector
    ]
    per_place.append(
        ce.t_place_profiles(
            a, branch, r_cap, depth, sigma_zero, d2_zero, g_zero, config
        )
    )
    if any(not options for options in per_place):
        return []
    dims = 3 + len(levels)
    vectors = [
        [profile.budget_vector(levels) for profile in options]
        for options in per_place
    ]
    minima = [
        tuple(min(vector[i] for vector in place) for i in range(dims))
        for place in vectors
    ]
    suffix = [tuple(0 for _ in range(dims))]
    for place_min in reversed(minima):
        last = suffix[0]
        suffix.insert(0, tuple(place_min[i] + last[i] for i in range(dims)))
    min_sums = suffix[0]

    e_cap = (ce.SUB2 if config is None else config).e_cap
    e_low = a + sum(b_vector)
    squeeze_on = (
        t2_squeeze and branch == "T2" and all(b != 2 for b in b_vector)
    )
    x_domain = (
        (NEG_INF,) if branch == "T2" else tuple(range(aux["d1"] + 1))
    )
    z_domain = (NEG_INF,) if sigma_zero else tuple(range(aux["sigma"] + 1))
    k_domain = (NEG_INF,) if d2_zero else tuple(range(aux["d2"] + 1))

    states = []
    for deg_e in range(e_low, e_cap + 1):
        squeeze_slack = 2 * (deg_e - e_low) if squeeze_on else 0
        for x_deg in x_domain:
            if x_deg != NEG_INF and x_deg < min_sums[0]:
                continue
            for z_deg in z_domain:
                if z_deg != NEG_INF and z_deg < min_sums[1]:
                    continue
                for k_deg in k_domain:
                    if k_deg != NEG_INF and k_deg < min_sums[2]:
                        continue
                    degstate = (
                        float(k_deg) if k_deg != NEG_INF else NEG_INF,
                        float(x_deg) if x_deg != NEG_INF else NEG_INF,
                        float(z_deg) if z_deg != NEG_INF else NEG_INF,
                        float(deg_e),
                    )
                    best = None
                    for profile in ce.inf_place_profiles(
                        a, branch, r_cap, depth, sigma_zero, d2_zero,
                        g_zero, degstate, config,
                    ):
                        chain = dict(profile.chain)
                        caps = (
                            aux["d1"] if x_deg == NEG_INF else x_deg,
                            aux["sigma"] if z_deg == NEG_INF else z_deg,
                            aux["d2"] if k_deg == NEG_INF else k_deg,
                            *(
                                g_caps[level]
                                if chain[level] == NEG_INF
                                else chain[level]
                                - (squeeze_slack if level == 6 else 0)
                                for level in levels
                            ),
                        )
                        if any(min_sums[i] > caps[i] for i in range(dims)):
                            continue
                        if ce._dfs_budget_witness(
                            per_place, vectors, suffix, caps
                        ) is None:
                            continue
                        count = len(profile.obligations)
                        if best is None or count < best["obligations"]:
                            best = {
                                "obligations": count,
                                "deg_g": {
                                    str(level): ce.encode_deg(value)
                                    for level, value in profile.chain
                                },
                            }
                    if best is not None:
                        states.append(
                            {
                                "deg_d2": ce.encode_deg(degstate[0]),
                                "deg_d1": ce.encode_deg(degstate[1]),
                                "deg_sigma": ce.encode_deg(degstate[2]),
                                "deg_e": deg_e,
                                "min_obligations": best["obligations"],
                                "deg_g": best["deg_g"],
                            }
                        )
    return states


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--window", choices=("sub2", "sub1"), default="sub2")
    parser.add_argument("--t2-squeeze", action="store_true", default=True)
    args = parser.parse_args()

    config = ce.CONFIGS[args.window]
    artifact = (
        "cascade_cones_qt_inf_rl.json"
        if args.window == "sub2"
        else "cascade_cones_sub1_qt_inf_rl.json"
    )
    ce.APPLY_RESIDUE_KILLS = True  # match the rl sweep semantics
    sweep = json.loads((ROOT / artifact).read_text(encoding="utf-8"))
    assert not sweep.get("partial_checkpoint"), "sweep artifact is partial"

    records = []
    total_states = 0
    for row in sweep["branches"]:
        if row["status"] != "survives":
            continue
        for case in row["survivor_cases"]:
            states = case_states(
                row["a_t"],
                tuple(row["b"]),
                row["branch"],
                case["sigma_zero"],
                case["d2_zero"],
                tuple(case["g_zero_levels"]),
                config,
                args.t2_squeeze,
            )
            assert states, (row["a_t"], row["b"], row["branch"])
            total_states += len(states)
            records.append(
                {
                    "a_t": row["a_t"],
                    "b": row["b"],
                    "branch": row["branch"],
                    "sigma_zero": case["sigma_zero"],
                    "d2_zero": case["d2_zero"],
                    "g_zero_levels": case["g_zero_levels"],
                    "state_count": len(states),
                    "states": states,
                }
            )

    payload = {
        "schema": 1,
        "description": (
            "Complete residual degree states per surviving flag case "
            f"({args.window}, q+t+inf, residue kills, T2 squeeze)"
        ),
        "window": args.window,
        "source_artifact": artifact,
        "case_count": len(records),
        "state_total": total_states,
        "cases": records,
    }
    out = ROOT / f"phase_d_states_{args.window}.json"
    out.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(f"{len(records)} cases, {total_states} states -> {out.name}")

    if args.window == "sub2":
        # Cross-check the open T2 cells against T5_T2_COLUMN.md's R-tables
        # projected to (deg e, deg sigma) with the witness (B, S) offsets.
        r_tables = {
            (9, (1, 0, 0, 0)): {(10, 2 + z) for z in range(7)},
            (8, (0, 0, 0, 0)): {(8, 3), (10, 5), (10, 6), (10, 7), (10, 8)},
            (8, (1, 0, 0, 0)): {(10, 5), (10, 6), (10, 7), (10, 8)},
            (8, (1, 1, 0, 0)): {(10, 4 + z) for z in range(5)},
            (7, (1, 0, 0, 0)): {(8, 3), (10, 7), (10, 8)},
            (7, (1, 1, 0, 0)): {(10, 7), (10, 8)},
            (7, (1, 1, 1, 0)): {(10, 6), (10, 7), (10, 8)},
            (7, (3, 0, 0, 0)): {(10, 7), (10, 8)},
        }
        for (a, b), expected in sorted(r_tables.items()):
            actual = set()
            for record in records:
                if (
                    record["branch"] == "T2"
                    and record["a_t"] == a
                    and tuple(record["b"]) == b
                ):
                    for state in record["states"]:
                        if state["deg_sigma"] != "-inf":
                            actual.add((state["deg_e"], state["deg_sigma"]))
            extra = actual - expected
            print(
                f"  T2 a{a} b{''.join(map(str, b))}: engine {sorted(actual)}"
                f" vs R-table {sorted(expected)}"
                + ("  EXTRA: " + str(sorted(extra)) if extra else "  (subset ok)"
                   if actual <= expected else "")
            )


if __name__ == "__main__":
    main()
