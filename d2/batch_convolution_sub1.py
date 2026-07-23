#!/usr/bin/env python3
"""Batch convolution-descent runner over the sub1 Phase-D residual worklist.

NEW file. Reuses the landed sub2 machinery verbatim (imports
``batch_convolution_sub2`` for the per-state ansatz build + gauge descent +
process-isolated Runner + tier/triage helpers); it does NOT modify it. The
only sub1-specific logic here is (a) loading/deduping ``phase_d_states_sub1``,
(b) a free TRANSFER pass against the sub2 round-1 verdicts, and (c) a priority
schedule for the fresh gauge run with incremental checkpointing.

Window fact justifying the transfer pass: the master identity
``f31 = sum_f Phi^f e^(21-3f) h_f == 0`` is window-INDEPENDENT (same h_f, same
Phi); only the degree caps and stratum data differ between windows. So a sub2
degree-state verdict applies verbatim to any sub1 state with the identical
dedup tuple ``(a_t, d1_zero, sigma_zero, d2_zero, deg_d2, deg_d1, deg_sigma,
deg_e)``. Matching kills are inherited (same pending-audit status); matching
non-kills (UNRESOLVED/FORCED/SKIPPED) are recorded and skipped from the fresh
run so budget lands only on genuinely new tuples.

Kills are all PENDING AUDIT (same-author pipeline).
"""

from __future__ import annotations

import json
import os
import time
from pathlib import Path

import batch_convolution_sub2 as b2  # landed machinery, imported read-only

ROOT = Path(__file__).resolve().parent
WORKLIST = ROOT / "phase_d_states_sub1.json"
SUB2_JSON = ROOT / "batch_convolution_sub2.json"
OUT_JSON = ROOT / "batch_convolution_sub1.json"
CKPT_JSON = ROOT / "batch_convolution_sub1_gauge_raw.json"  # incremental

PER_STATE_TIMEOUT = float(os.environ.get("BATCH_PER_STATE_TIMEOUT", "90"))
TOTAL_WALL_BUDGET = float(os.environ.get("BATCH_TOTAL_WALL_BUDGET", str(45 * 60)))
C_VALUE = b2.C_VALUE
FLOOR_BUDGET = b2.FLOOR_BUDGET


def canon(rec) -> tuple:
    """Window-independent dedup/transfer key (deg fields stringified so the
    '-inf' sentinel for zeroed variables matches across the two artifacts)."""
    return (rec["a_t"], rec["d1_zero"], rec["sigma_zero"], rec["d2_zero"],
            str(rec["deg_d2"]), str(rec["deg_d1"]), str(rec["deg_sigma"]),
            rec["deg_e"])


def load_unique_states():
    data = json.loads(WORKLIST.read_text(encoding="utf-8"))
    uniq: dict[tuple, dict] = {}
    for case in data["cases"]:
        is_t2 = case["branch"] == "T2"
        for st in case["states"]:
            rec_seed = {
                "a_t": case["a_t"], "d1_zero": is_t2, "branch": case["branch"],
                "sigma_zero": case["sigma_zero"], "d2_zero": case["d2_zero"],
                "deg_d2": st["deg_d2"], "deg_d1": st["deg_d1"],
                "deg_sigma": st["deg_sigma"], "deg_e": st["deg_e"],
            }
            key = canon(rec_seed)
            rec = uniq.get(key)
            if rec is None:
                rec_seed["raw_count"] = 0
                rec_seed["cells"] = set()
                rec = uniq[key] = rec_seed
            rec["raw_count"] += 1
            rec["cells"].add((case["a_t"], tuple(case["b"]), case["branch"]))
    for rec in uniq.values():
        rec["cell_count"] = len(rec["cells"])
        rec["cells"] = sorted(
            {(a, "".join(map(str, b)), br) for a, b, br in rec["cells"]})
    return uniq, data["state_total"]


def load_sub2_verdicts():
    data = json.loads(SUB2_JSON.read_text(encoding="utf-8"))
    att = {}
    for s in data["states"]:
        att[canon(s)] = s
    return att


def fresh_priority(rec):
    """Schedule key for the fresh gauge run. Groups (ascending group index):
      0: T2 (tier 1) at a in {9,10}   -- cheapest first, a=10 before a=9
      1: T1 constant-E (tier 2) at a in {9,10}, a=10 first (round-1 killer)
      2: T1 constant-E (tier 2) at a < 9, descending a
    Everything else returns group 9 = not scheduled this run.
    """
    t = b2.tier(rec)
    a = rec["a_t"]
    m = rec["deg_e"] - a
    tie = (str(rec["deg_d2"]), str(rec["deg_d1"]), str(rec["deg_sigma"]))
    if t == 1 and a in (9, 10):
        return (0, m, -a, tie)
    if t == 2 and a in (9, 10):
        return (1, -a, m, tie)
    if t == 2:
        return (2, -a, m, tie)
    return (9,)


def build_final(uniq, raw_total, att, attempted):
    """Assemble batch_convolution_sub1.json (schema mirrors sub2 + transfers)."""
    records = list(uniq.values())

    def tc(recs):
        return {t: sum(1 for r in recs if b2.tier(r) == t) for t in (1, 2, 3, 4)}

    # transfers: sub1 tuples whose canon matches a sub2 attempted state
    transferred = []
    for key, rec in uniq.items():
        if key in att:
            s = att[key]
            transferred.append({
                **{k: rec[k] for k in (
                    "a_t", "branch", "d1_zero", "sigma_zero", "d2_zero",
                    "deg_d2", "deg_d1", "deg_sigma", "deg_e",
                    "raw_count", "cell_count", "cells")},
                "tier": b2.tier(rec),
                "m": rec["deg_e"] - rec["a_t"],
                "transferred_verdict": s["final_verdict"],
                "sub2_decided_by": s.get("decided_by"),
            })
    transfer_kills = [t for t in transferred if t["transferred_verdict"] in
                      ("CONTRADICTION", "STATE_KILLED_BY_DEGREE_DROP")]

    # fresh attempted verdicts
    fresh_census = {}
    for e in attempted:
        fresh_census[e["verdict"]] = fresh_census.get(e["verdict"], 0) + 1
    fresh_kills = [e for e in attempted if e["verdict"] in
                   ("CONTRADICTION", "STATE_KILLED_BY_DEGREE_DROP")]

    # combined kill census (transfer + fresh), all PENDING AUDIT
    all_kills = (
        [{"origin": "transfer", **t} for t in transfer_kills]
        + [{"origin": "fresh", **e} for e in fresh_kills])
    kill_raw = (sum(t["raw_count"] for t in transfer_kills)
                + sum(e["raw_count"] for e in fresh_kills))

    attempted_keys = set(att) | {canon(e) for e in attempted}
    unattempted = [r for k, r in uniq.items() if k not in attempted_keys]

    payload = {
        "schema": 2,
        "description": (
            "Convolution-descent verdicts over the deduped sub1 Phase-D "
            "residual degree states. Free TRANSFER pass from sub2 round-1 "
            "(master identity f31 is window-independent) + fresh gauge run on "
            "new tuples. c=-1/6630; q-root support conditions dropped (sound "
            "over-approximation); ALL kills PENDING AUDIT."),
        "source_worklist": WORKLIST.name,
        "sub2_source": SUB2_JSON.name,
        "c": str(C_VALUE),
        "floor_budget": FLOOR_BUDGET,
        "per_state_timeout_s": PER_STATE_TIMEOUT,
        "raw_state_total": raw_total,
        "unique_state_total": len(records),
        "tier_counts": tc(records),
        "tier_legend": {"1": "T2", "2": "T1 constant-E (deg_e==a_t)",
                        "3": "T1 sigma_zero or d2_zero", "4": "other T1"},
        "unique_by_a_tier": {
            str(a): {str(t): sum(1 for r in records
                                 if r["a_t"] == a and b2.tier(r) == t)
                     for t in (1, 2, 3, 4)}
            for a in sorted({r["a_t"] for r in records})},
        # transfer pass
        "transfer_count": len(transferred),
        "transfer_verdict_census": _census([t["transferred_verdict"]
                                            for t in transferred]),
        "transferred_kills": transfer_kills,
        "transferred_kill_raw_coverage": sum(t["raw_count"]
                                             for t in transfer_kills),
        # fresh pass
        "fresh_attempted": len(attempted),
        "fresh_verdict_census": fresh_census,
        "fresh_kills_pending_audit": fresh_kills,
        "fresh_kill_raw_coverage": sum(e["raw_count"] for e in fresh_kills),
        # combined
        "kill_count": len(all_kills),
        "kill_raw_state_coverage": kill_raw,
        "kills_pending_audit": all_kills,
        "attempted_unique": len(transferred) + len(attempted),
        "unattempted_unique": len(unattempted),
        "unattempted_by_tier": tc(unattempted),
        "fresh_states": attempted,
    }
    return payload


def _census(vals):
    out = {}
    for v in vals:
        out[v] = out.get(v, 0) + 1
    return out


def checkpoint(uniq, raw_total, att, attempted, meta):
    payload = build_final(uniq, raw_total, att, attempted)
    payload["_run_meta"] = meta
    CKPT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
                         encoding="utf-8")


def main():
    uniq, raw_total = load_unique_states()
    att = load_sub2_verdicts()
    records = list(uniq.values())
    print(f"raw {raw_total}; unique {len(records)}; "
          f"sub2 attempted overlap {sum(1 for k in uniq if k in att)}; "
          f"sub2 kills present {sum(1 for k in uniq if k in att and att[k]['final_verdict'] in ('CONTRADICTION','STATE_KILLED_BY_DEGREE_DROP'))}")

    # resume: reload already-finished fresh states from checkpoint
    done = {}
    if CKPT_JSON.exists():
        try:
            prev = json.loads(CKPT_JSON.read_text(encoding="utf-8"))
            for e in prev.get("fresh_states", []):
                done[canon(e)] = e
            print(f"resume: {len(done)} fresh states already checkpointed")
        except Exception:
            done = {}

    # schedule fresh targets (not in sub2, in a scheduled group)
    scheduled = []
    for key, rec in uniq.items():
        if key in att:
            continue
        pr = fresh_priority(rec)
        if pr[0] == 9:
            continue
        scheduled.append((pr, key, rec))
    scheduled.sort(key=lambda x: x[0])
    print(f"fresh scheduled targets: {len(scheduled)} "
          f"(groups: T2[9,10], constE[9,10], constE[<9])")

    runner = b2.Runner()
    attempted = list(done.values())
    t_start = time.time()
    stop_reason = "completed all scheduled targets"
    n_new = 0
    for pr, key, rec in scheduled:
        if key in done:
            continue
        elapsed = time.time() - t_start
        if elapsed >= TOTAL_WALL_BUDGET:
            stop_reason = "total wall budget exhausted"
            break
        timeout = min(PER_STATE_TIMEOUT, TOTAL_WALL_BUDGET - elapsed)
        rec["gauge"] = True
        t0 = time.time()
        status, payload = runner.run(rec, timeout)
        dt = round(time.time() - t0, 1)
        entry = {k: rec[k] for k in (
            "a_t", "branch", "d1_zero", "sigma_zero", "d2_zero",
            "deg_d2", "deg_d1", "deg_sigma", "deg_e",
            "raw_count", "cell_count", "cells")}
        entry["tier"] = b2.tier(rec)
        entry["m"] = rec["deg_e"] - rec["a_t"]
        entry["seconds"] = dt
        if status == "ok":
            entry.update(payload)
        elif status == "timeout":
            entry["verdict"] = "SKIPPED_BUDGET"
            entry["reason"] = f"exceeded per-state timeout {timeout:.0f}s"
        else:
            entry["verdict"] = "ERROR"
            entry["error"] = payload
        attempted.append(entry)
        n_new += 1
        marker = ("***KILL***" if entry["verdict"] in
                  ("CONTRADICTION", "STATE_KILLED_BY_DEGREE_DROP") else "")
        print(f"g{pr[0]} tier{entry['tier']} a={rec['a_t']} "
              f"d2={rec['deg_d2']} d1={rec['deg_d1']} sig={rec['deg_sigma']} "
              f"e={rec['deg_e']} -> {entry['verdict']} ({dt}s) {marker}")
        # INCREMENTAL CHECKPOINT after every state
        meta = {"n_new_this_run": n_new, "n_total_fresh": len(attempted),
                "elapsed_s": round(time.time() - t_start, 1),
                "stop_reason": "running"}
        checkpoint(uniq, raw_total, att, attempted, meta)
    runner.close()

    meta = {"n_new_this_run": n_new, "n_total_fresh": len(attempted),
            "elapsed_s": round(time.time() - t_start, 1),
            "stop_reason": stop_reason,
            "scheduled_total": len(scheduled)}
    payload = build_final(uniq, raw_total, att, attempted)
    payload["_run_meta"] = meta
    CKPT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
                         encoding="utf-8")
    OUT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8")
    fc = payload["fresh_verdict_census"]
    print(f"\nfresh census {fc}")
    print(f"fresh kills {len(payload['fresh_kills_pending_audit'])}; "
          f"transfer kills {len(payload['transferred_kills'])}; "
          f"combined {payload['kill_count']} "
          f"(raw coverage {payload['kill_raw_state_coverage']})")
    print(f"stop: {stop_reason} (new {n_new}, total fresh {len(attempted)})")
    print(f"-> {OUT_JSON.name}")


if __name__ == "__main__":
    main()
