#!/usr/bin/env python3
"""Round-2 batch convolution-descent runner (NEW file; landed batch is FROZEN).

Wraps the committed ``batch_convolution_sub2`` machinery (imports its
``load_unique_states``, ``triage_sort``, ``tier``, ``Runner`` and the gauge
worker via ``_build_and_descend``) WITHOUT modifying it, and adds the one thing
round 1 lacked: an incremental checkpoint written after EVERY state, so a stall
never loses completed work.

Scope (env-configurable, defaults target task 1):
  BATCH2_START_INDEX (default 194)  first triage index to attempt
  BATCH2_END_INDEX   (default 585)  one past last (194..584 = the 391 remaining
                                    tier-2 a_t=10 unique states)
  BATCH2_WALL_BUDGET (default 2400) total wall seconds
  BATCH2_OUT         (default batch_convolution_sub2_round2.json)

Gauge mode is forced on (the mechanism that killed 22/29 of the first a10
block). Per-state timeout, process isolation, and the ansatz are inherited
unchanged from the landed runner.

Every kill is a CANDIDATE kill, PENDING AUDIT.
"""
from __future__ import annotations

import json
import os
import time
from pathlib import Path

import batch_convolution_sub2 as B

ROOT = Path(__file__).resolve().parent
START_INDEX = int(os.environ.get("BATCH2_START_INDEX", "194"))
END_INDEX = int(os.environ.get("BATCH2_END_INDEX", "585"))
WALL_BUDGET = float(os.environ.get("BATCH2_WALL_BUDGET", str(40 * 60)))
OUT = ROOT / os.environ.get("BATCH2_OUT", "batch_convolution_sub2_round2.json")
PER_STATE_TIMEOUT = B.PER_STATE_TIMEOUT

KILL_VERDICTS = ("CONTRADICTION", "STATE_KILLED_BY_DEGREE_DROP")


def build_payload(records_all, ordered, attempted, census, stop_reason,
                  stopped_index, done):
    """Round-1-schema-compatible payload (+ 'round' and honest coverage)."""
    kills = [s for s in attempted if s.get("verdict") in KILL_VERDICTS]
    attempted_keys = {
        (s["a_t"], s["branch"], s["sigma_zero"], s["d2_zero"],
         str(s["deg_d2"]), str(s["deg_d1"]), str(s["deg_sigma"]), s["deg_e"])
        for s in attempted
    }
    # unattempted across the WHOLE worklist (so coverage is honest, not just
    # relative to this round's index window)
    unattempted = [
        r for r in ordered
        if (r["a_t"], r["branch"], r["sigma_zero"], r["d2_zero"],
            str(r["deg_d2"]), str(r["deg_d1"]), str(r["deg_sigma"]),
            r["deg_e"]) not in attempted_keys
    ]
    tier_counts_all = {t: sum(1 for r in ordered if B.tier(r) == t)
                       for t in (1, 2, 3, 4)}
    tier_counts_open = {t: sum(1 for r in unattempted if B.tier(r) == t)
                        for t in (1, 2, 3, 4)}
    return {
        "schema": 2,
        "round": 2,
        "complete": done,
        "description": (
            "Round-2 gauge-mode convolution-descent verdicts over the "
            "remaining tier-2 a_t=10 unique Phase-D degree states (and beyond "
            "if budget allowed). Master identity f31; c=-1/6630; lc(e generic "
            "part)=gamma nonzero parameter; q-root support conditions dropped "
            "(sound over-approximation). Kills are CANDIDATE kills PENDING "
            "AUDIT. Round-1 landed artifacts are frozen and untouched."),
        "source_worklist": B.WORKLIST.name,
        "gauge_normalized": True,
        "c": str(B.C_VALUE),
        "floor_budget": B.FLOOR_BUDGET,
        "per_state_timeout_s": PER_STATE_TIMEOUT,
        "wall_budget_s": WALL_BUDGET,
        "index_window": [START_INDEX, END_INDEX],
        "raw_state_total": records_all,
        "unique_state_total": len(ordered),
        "tier_counts": tier_counts_all,
        "attempted_unique_this_round": len(attempted),
        "stopped_at_index": stopped_index,
        "stop_reason": stop_reason,
        "verdict_census": census,
        "kill_count": len(kills),
        "kill_raw_state_coverage": sum(r["raw_count"] for r in kills),
        "unattempted_unique_global": len(unattempted),
        "unattempted_by_tier_global": tier_counts_open,
        "coverage_note": (
            "attempted_unique_this_round counts ONLY states run in round 2; "
            "unattempted_by_tier_global is over the full 1782-state worklist "
            "(round-1 attempts are NOT reflected here as attempted, so this "
            "file is self-contained for round 2). UNRESOLVED never certifies "
            "survival (q-root support dropped)."),
        "kills_pending_audit": kills,
        "states": attempted,
    }


def main():
    records, raw_total = B.load_unique_states()
    ordered = B.triage_sort(records)
    for rec in ordered:
        rec["gauge"] = True

    print(f"raw {raw_total}; unique {len(ordered)}; "
          f"window [{START_INDEX},{END_INDEX}); wall {WALL_BUDGET:.0f}s")

    runner = B.Runner()
    attempted = []
    census = {}
    t_start = time.time()
    stop_reason = "completed index window"
    stopped_index = min(END_INDEX, len(ordered))

    def checkpoint(done, sreason, sidx):
        payload = build_payload(raw_total, ordered, attempted, census,
                                sreason, sidx, done)
        tmp = OUT.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
                       encoding="utf-8")
        tmp.replace(OUT)  # atomic swap so a crash mid-write can't corrupt

    try:
        for idx in range(START_INDEX, min(END_INDEX, len(ordered))):
            elapsed = time.time() - t_start
            if elapsed >= WALL_BUDGET:
                stop_reason = "wall budget exhausted"
                stopped_index = idx
                break
            rec = ordered[idx]
            remaining = WALL_BUDGET - elapsed
            timeout = min(PER_STATE_TIMEOUT, remaining)

            t0 = time.time()
            status, payload = runner.run(rec, timeout)
            dt = round(time.time() - t0, 1)

            entry = {k: rec[k] for k in (
                "a_t", "branch", "d1_zero", "sigma_zero", "d2_zero",
                "deg_d2", "deg_d1", "deg_sigma", "deg_e",
                "raw_count", "cell_count", "cells")}
            entry["tier"] = B.tier(rec)
            entry["m"] = rec["deg_e"] - rec["a_t"]
            entry["triage_index"] = idx
            entry["gauge"] = True
            entry["seconds"] = dt
            if status == "ok":
                entry.update(payload)
            elif status == "timeout":
                entry["verdict"] = "SKIPPED_BUDGET"
                entry["reason"] = f"exceeded per-state timeout {timeout:.0f}s"
            else:
                entry["verdict"] = "ERROR"
                entry["error"] = payload

            census[entry["verdict"]] = census.get(entry["verdict"], 0) + 1
            attempted.append(entry)
            checkpoint(False, "in progress", idx + 1)  # incremental: nothing lost

            tag = entry["verdict"]
            print(f"[{idx}] a{rec['a_t']} d2={rec['deg_d2']} d1={rec['deg_d1']} "
                  f"s={rec['deg_sigma']} -> {tag} (stop "
                  f"{entry.get('stopping_degree')}, {dt}s)")
    finally:
        runner.close()

    done = stopped_index >= min(END_INDEX, len(ordered))
    checkpoint(done, stop_reason, stopped_index)
    kills = sum(1 for s in attempted if s["verdict"] in KILL_VERDICTS)
    print(f"\ncensus {census}")
    print(f"attempted {len(attempted)}; kills {kills}; "
          f"stopped at index {stopped_index}: {stop_reason}")
    print(f"-> {OUT.name}")


if __name__ == "__main__":
    main()
