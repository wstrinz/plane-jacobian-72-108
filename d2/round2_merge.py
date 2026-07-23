#!/usr/bin/env python3
"""Merge the round-2 part files into the final batch_convolution_sub2_round2.json.

Part 1: indices 194-397 (main 40-min run, checkpointed in
batch_convolution_sub2_round2.json), part 2: indices 398+ (continuation,
batch_convolution_sub2_round2_part2.json). Disjoint index windows of the same
deterministic triage order; the merge concatenates states and recomputes the
census/coverage via batch_convolution_sub2_round2.build_payload.
"""
import json
from pathlib import Path

import batch_convolution_sub2 as B
import batch_convolution_sub2_round2 as R2

ROOT = Path(__file__).resolve().parent
P1 = ROOT / "batch_convolution_sub2_round2.json"
P2 = ROOT / "batch_convolution_sub2_round2_part2.json"

d1 = json.loads(P1.read_text(encoding="utf-8"))
d2 = json.loads(P2.read_text(encoding="utf-8")) if P2.exists() else None

states = list(d1["states"]) + (list(d2["states"]) if d2 else [])
idx_seen = [s["triage_index"] for s in states]
assert len(idx_seen) == len(set(idx_seen)), "overlapping part windows"
states.sort(key=lambda s: s["triage_index"])

census = {}
for s in states:
    census[s["verdict"]] = census.get(s["verdict"], 0) + 1

records, raw_total = B.load_unique_states()
ordered = B.triage_sort(records)

stop_parts = [f"part1 [{d1['index_window'][0]},{d1['index_window'][1]}): "
              f"stopped {d1['stopped_at_index']} ({d1['stop_reason']})"]
stopped_index = d1["stopped_at_index"]
if d2:
    stop_parts.append(
        f"part2 [{d2['index_window'][0]},{d2['index_window'][1]}): "
        f"stopped {d2['stopped_at_index']} ({d2['stop_reason']})")
    stopped_index = max(stopped_index, d2["stopped_at_index"])

payload = R2.build_payload(raw_total, ordered, states, census,
                           "; ".join(stop_parts), stopped_index,
                           stopped_index >= 585)
payload["index_window"] = [194, 585]
payload["wall_budget_s"] = (d1["wall_budget_s"]
                            + (d2["wall_budget_s"] if d2 else 0))
payload["parts"] = [
    {"file": P1.name, "window": d1["index_window"],
     "attempted": d1["attempted_unique_this_round"],
     "census": d1["verdict_census"], "wall_budget_s": d1["wall_budget_s"]},
]
if d2:
    payload["parts"].append(
        {"file": P2.name, "window": d2["index_window"],
         "attempted": d2["attempted_unique_this_round"],
         "census": d2["verdict_census"], "wall_budget_s": d2["wall_budget_s"]})

P1.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
              encoding="utf-8")
print("merged census", census)
print("kills", payload["kill_count"], "raw coverage",
      payload["kill_raw_state_coverage"])
print("attempted", payload["attempted_unique_this_round"],
      "stopped_at", payload["stopped_at_index"])
print("->", P1.name)
