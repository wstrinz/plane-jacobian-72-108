#!/usr/bin/env python3
"""frontier_gen.py — machine-generated cross-front frontier summary (writeup gap G6).

Reads the authoritative cascade / alternate-regime JSON artifacts and regenerates
FRONTIER.md: one table per front, a cross-front totals line, a generation-timestamp
line sourced from each artifact's own file metadata, and a "how to regenerate" line.

EVERY figure in FRONTIER.md is computed from the JSONs here; nothing is hand-typed
into the output. The script fails loudly (nonzero exit) if a required artifact is
missing or if any artifact is a partial checkpoint (partial_checkpoint == true).

Usage:
    python frontier_gen.py            # writes FRONTIER.md next to the artifacts
"""
import json
import os
import sys
import datetime

HERE = os.path.dirname(os.path.abspath(__file__))

# Required artifacts (absence or partial_checkpoint == fatal).
SUB2_ART = "cascade_cones_qt_inf_rl.json"
SUB1_ART = "cascade_cones_sub1_qt_inf_rl.json"
ALT_ART = "alt_inf_sweep.json"
REQUIRED = [SUB2_ART, SUB1_ART, ALT_ART]

# Optional per-flag-case state worklists (used only if present).
PHASE_D_SUB2 = "phase_d_states_sub2.json"
PHASE_D_SUB1 = "phase_d_states_sub1.json"


def die(msg):
    sys.stderr.write("frontier_gen.py: FATAL: " + msg + "\n")
    sys.exit(1)


def load(name, required=True):
    path = os.path.join(HERE, name)
    if not os.path.exists(path):
        if required:
            die("required artifact missing: %s" % name)
        return None, None
    try:
        with open(path, "r", encoding="utf-8") as fh:
            data = json.load(fh)
    except Exception as exc:  # noqa: BLE001
        die("could not parse %s: %s" % (name, exc))
    if isinstance(data, dict) and data.get("partial_checkpoint") is True:
        die("artifact %s is a partial checkpoint (partial_checkpoint=true); "
            "refusing to generate a frontier from incomplete data." % name)
    mtime = datetime.datetime.fromtimestamp(os.path.getmtime(path))
    return data, mtime


def cascade_front(data, name):
    """Summarize a cascade_cones_*_inf_rl.json front from its own branch list."""
    branches = data.get("branches")
    if not isinstance(branches, list):
        die("%s has no 'branches' list" % name)
    surv = [b for b in branches if b.get("status") == "survives"]
    t1 = [b for b in surv if b.get("branch") == "T1"]
    t2 = [b for b in surv if b.get("branch") == "T2"]
    fc = sum(int(b.get("survivor_case_count", 0)) for b in surv)
    fc1 = sum(int(b.get("survivor_case_count", 0)) for b in t1)
    fc2 = sum(int(b.get("survivor_case_count", 0)) for b in t2)
    a_vals = sorted({b.get("a_t") for b in surv if b.get("a_t") is not None})
    summary = data.get("summary", {})
    # Cross-check the engine's own survivor tally against our recount.
    reported = summary.get("surviving_branches")
    if reported is not None and int(reported) != len(surv):
        die("%s: summary.surviving_branches=%s disagrees with recount %d"
            % (name, reported, len(surv)))
    return {
        "window": data.get("window"),
        "places": data.get("places"),
        "depth": data.get("depth"),
        "surv": len(surv), "t1": len(t1), "t2": len(t2),
        "fc": fc, "fc1": fc1, "fc2": fc2,
        "a_lo": a_vals[0] if a_vals else None,
        "a_hi": a_vals[-1] if a_vals else None,
        "open_processed": summary.get("open_branches_processed"),
        "killed": summary.get("engine_killed_pending_audit"),
    }


def alt_front(data, name):
    summary = data.get("summary", {})
    branches = data.get("branches", [])
    t1 = [b for b in branches if b.get("branch") == "T1"]
    t2 = [b for b in branches if b.get("branch") == "T2"]
    n = summary.get("n_branches")
    if n is not None and int(n) != len(branches):
        die("%s: summary.n_branches=%s disagrees with recount %d"
            % (name, n, len(branches)))
    return {
        "n": len(branches), "t1": len(t1), "t2": len(t2),
        "open": summary.get("branches_OPEN"),
        "killed": summary.get("branches_KILLED"),
        "total_states": summary.get("total_degree_states"),
        "surv_states": summary.get("surviving_states"),
        "killed_states": summary.get("killed_states"),
    }


def phase_d_states(data):
    if data is None:
        return None
    return {
        "case_count": data.get("case_count"),
        "state_total": data.get("state_total"),
        "source": data.get("source_artifact"),
    }


def main():
    sub2_d, sub2_t = load(SUB2_ART)
    sub1_d, sub1_t = load(SUB1_ART)
    alt_d, alt_t = load(ALT_ART)
    pd2_d, pd2_t = load(PHASE_D_SUB2, required=False)
    pd1_d, pd1_t = load(PHASE_D_SUB1, required=False)

    sub2 = cascade_front(sub2_d, SUB2_ART)
    sub1 = cascade_front(sub1_d, SUB1_ART)
    alt = alt_front(alt_d, ALT_ART)
    pd2 = phase_d_states(pd2_d)
    pd1 = phase_d_states(pd1_d)

    # State counts where available (per-flag-case worklists).
    sub2_states = pd2["state_total"] if pd2 else None
    sub1_states = pd1["state_total"] if pd1 else None
    # Cross-check the sub2 phase-D case count against the cascade flag-case total.
    if pd2 and pd2["case_count"] is not None and int(pd2["case_count"]) != sub2["fc"]:
        die("phase_d_states_sub2 case_count=%s disagrees with sub2 flag-case total %d"
            % (pd2["case_count"], sub2["fc"]))

    total_units = sub2["surv"] + sub1["surv"] + alt["n"]
    total_flag = sub2["fc"] + sub1["fc"]  # alt tracks degree-states, not flag cases

    def dash(x):
        return "—" if x is None else str(x)

    lines = []
    lines.append("# FRONTIER — live fronts of the (72,108) program (machine-generated)")
    lines.append("")
    lines.append("> **DO NOT HAND-EDIT.** This file is regenerated by `frontier_gen.py` "
                 "from the cascade / alternate-regime JSON artifacts. Every figure below "
                 "is read from those JSONs (branch lists + their own summary blocks); "
                 "none is typed by hand. Closes writeup gap **G6**.")
    lines.append("")
    lines.append("Generated: %s (local clock)." % datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    lines.append("")
    lines.append("Source-artifact timestamps (file mtime — the artifacts carry no embedded date):")
    lines.append("")
    lines.append("| Artifact | Modified | partial_checkpoint |")
    lines.append("|---|---|---|")
    lines.append("| `%s` | %s | %s |" % (SUB2_ART, sub2_t.strftime("%Y-%m-%d %H:%M:%S"),
                                         sub2_d.get("partial_checkpoint")))
    lines.append("| `%s` | %s | %s |" % (SUB1_ART, sub1_t.strftime("%Y-%m-%d %H:%M:%S"),
                                         sub1_d.get("partial_checkpoint")))
    lines.append("| `%s` | %s | %s |" % (ALT_ART, alt_t.strftime("%Y-%m-%d %H:%M:%S"),
                                         alt_d.get("partial_checkpoint")))
    if pd2:
        lines.append("| `%s` | %s | (worklist) |" % (PHASE_D_SUB2, pd2_t.strftime("%Y-%m-%d %H:%M:%S")))
    if pd1:
        lines.append("| `%s` | %s | (worklist) |" % (PHASE_D_SUB1, pd1_t.strftime("%Y-%m-%d %H:%M:%S")))
    lines.append("")

    # ---- Sub2 front ----
    lines.append("## Sub2 cells — `%s` (%s, depth %s)" % (SUB2_ART, sub2["places"], dash(sub2["depth"])))
    lines.append("")
    lines.append("| Metric | T1 | T2 | Total |")
    lines.append("|---|---|---|---|")
    lines.append("| Surviving cells | %d | %d | **%d** |" % (sub2["t1"], sub2["t2"], sub2["surv"]))
    lines.append("| Flag cases | %d | %d | **%d** |" % (sub2["fc1"], sub2["fc2"], sub2["fc"]))
    lines.append("")
    lines.append("- a-range of survivors: a ∈ [%s, %s]." % (dash(sub2["a_lo"]), dash(sub2["a_hi"])))
    lines.append("- Engine: %s open branches processed, %s killed (pending audit)."
                 % (dash(sub2["open_processed"]), dash(sub2["killed"])))
    if pd2:
        lines.append("- Residual degree-states (`%s`): %s states over %s flag cases."
                     % (PHASE_D_SUB2, dash(sub2_states), dash(pd2["case_count"])))
    lines.append("")

    # ---- Sub1 front ----
    lines.append("## Sub1 standard regime — `%s` (%s, depth %s)" % (SUB1_ART, sub1["places"], dash(sub1["depth"])))
    lines.append("")
    lines.append("| Metric | T1 | T2 | Total |")
    lines.append("|---|---|---|---|")
    lines.append("| Surviving branches | %d | %d | **%d** |" % (sub1["t1"], sub1["t2"], sub1["surv"]))
    lines.append("| Flag cases | %d | %d | **%d** |" % (sub1["fc1"], sub1["fc2"], sub1["fc"]))
    lines.append("")
    lines.append("- a-range of survivors: a ∈ [%s, %s]." % (dash(sub1["a_lo"]), dash(sub1["a_hi"])))
    lines.append("- Engine: %s open branches processed, %s killed (pending audit)."
                 % (dash(sub1["open_processed"]), dash(sub1["killed"])))
    if pd1:
        lines.append("- Residual degree-states (`%s`): %s states." % (PHASE_D_SUB1, dash(sub1_states)))
    else:
        lines.append("- Residual degree-states: `%s` not present — no state count available." % PHASE_D_SUB1)
    lines.append("")

    # ---- Alternate regime front ----
    lines.append("## Sub1 alternate regime — `%s`" % ALT_ART)
    lines.append("")
    lines.append("| Metric | T1 | T2 | Total |")
    lines.append("|---|---|---|---|")
    lines.append("| Branches (OPEN) | %d | %d | **%d** |" % (alt["t1"], alt["t2"], alt["n"]))
    lines.append("")
    lines.append("- Verdicts: %s OPEN, %s killed." % (dash(alt["open"]), dash(alt["killed"])))
    lines.append("- Degree-states: %s total, %s surviving, %s killed."
                 % (dash(alt["total_states"]), dash(alt["surv_states"]), dash(alt["killed_states"])))
    lines.append("")

    # ---- Cross-front totals ----
    lines.append("## Cross-front totals")
    lines.append("")
    lines.append("**%d surviving cells/branches** = %d sub2 cells + %d sub1 branches + %d alternate-regime branches. "
                 "**%d flag cases** across the two flag-case fronts (%d sub2 + %d sub1); the alternate regime is "
                 "tracked as %s degree-states, not flag cases."
                 % (total_units, sub2["surv"], sub1["surv"], alt["n"],
                    total_flag, sub2["fc"], sub1["fc"], dash(alt["total_states"])))
    lines.append("")
    lines.append("## How to regenerate")
    lines.append("")
    lines.append("```")
    lines.append("python frontier_gen.py")
    lines.append("```")
    lines.append("")
    lines.append("Reads `%s`, `%s`, `%s` (required; fails loudly if any is missing or a "
                 "partial checkpoint) plus `%s` / `%s` if present. Rerun after any cascade "
                 "or alternate-regime sweep to refresh these figures."
                 % (SUB2_ART, SUB1_ART, ALT_ART, PHASE_D_SUB2, PHASE_D_SUB1))
    lines.append("")

    out = os.path.join(HERE, "FRONTIER.md")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    sys.stderr.write("frontier_gen.py: wrote %s (%d cells/branches, %d flag cases)\n"
                     % (out, total_units, total_flag))


if __name__ == "__main__":
    main()
