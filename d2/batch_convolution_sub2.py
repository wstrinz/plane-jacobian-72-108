#!/usr/bin/env python3
"""Batch convolution-descent runner over the sub2 Phase-D residual worklist.

Reads ``phase_d_states_sub2.json`` (schema in ``phase_d_states.py``), dedupes
its degree states to unique tuples

    (a_t, d1_zero(=branch==T2), sigma_zero, d2_zero,
     deg_d2, deg_d1, deg_sigma, deg_e)

- justified by the master identity f31 = sum_f Phi^f e^(21-3f) h_f == 0 being
INDEPENDENT of the g_zero flag structure, so one degree-state verdict applies
to every flag case sharing those variable zero-flags and degrees - then runs
the landed exact ``convolution_descent`` driver on each unique state.

Ansatz per state (NEW code; does not touch the landed driver):
  * e     = (y+1)**a_t * (generic poly of degree deg_e - a_t)  [encodes v_t(e)=a_t]
  * d1    = 0 for T2 (d1_zero), else generic of degree deg_d1
  * sigma = 0 if sigma_zero, else generic of degree deg_sigma
  * d2    = 0 if d2_zero, else generic of degree deg_d2
  * d0    = (d2**2 + sigma)/4   (via build_ansatz sigma path)
  The q-root support conditions on e are DROPPED: this is a sound
  over-approximation (a larger ansatz family), so a kill here kills the
  original; an UNRESOLVED here does not certify survival of the original.

Descent window: start = 1 + max_f( 34 f + (21-3 f) deg_e + maxdeg_y h_f|state ),
floor = start - FLOOR_BUDGET (default 14).  c = -1/6630 is fixed.

Interpretation of a forced leading coefficient set to zero: per the task, a
forced substitution that sends the LEADING coefficient of d2/d1/sigma/e's
generic part to 0 does NOT continue the state - it contradicts the recorded
degree, so the verdict is STATE_KILLED_BY_DEGREE_DROP and the state stops.
(A leading coeff forced to a nonzero VALUE is fine and continues.)

Verdicts recorded per attempted state:
  CONTRADICTION | STATE_KILLED_BY_DEGREE_DROP | FORCED | UNRESOLVED | SKIPPED_BUDGET

Triage order (stop when the total wall budget is spent; where we stopped is
recorded):
  1. all T2 states
  2. T1 with deg_e == a_t (constant-E analogues)
  3. T1 with sigma_zero or d2_zero
  4. everything else
Within a tier states are ordered by ascending m = deg_e - a_t (cheapest first)
so the fixed wall budget covers as many states as possible.
"""

from __future__ import annotations

import json
import multiprocessing as mp
import os
import time
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent
WORKLIST = ROOT / "phase_d_states_sub2.json"
OUT_JSON = ROOT / "batch_convolution_sub2.json"

FLOOR_BUDGET = int(os.environ.get("BATCH_FLOOR_BUDGET", "14"))
C_VALUE = sp.Rational(-1, 6630)
PER_STATE_TIMEOUT = float(os.environ.get("BATCH_PER_STATE_TIMEOUT", "90"))
TOTAL_WALL_BUDGET = float(os.environ.get("BATCH_TOTAL_WALL_BUDGET",
                                         str(45 * 60)))
MAX_STATES = int(os.environ.get("BATCH_MAX_STATES", "0"))  # 0 = no cap
# Gauge mode (BATCH_GAUGE=1): declare the leading coefficient of e's generic
# part a NONZERO PARAMETER gamma instead of a forceable unknown.
# Soundness: gamma != 0 is EXACTLY the degree-exactness of e in the recorded
# state (deg e = deg_e), so nothing is assumed beyond the state itself, and
# gamma stays a fully general symbol (no generality lost).  Residuals whose
# only free symbol is gamma are accumulated as necessary polynomial
# constraints; if they admit no common nonzero root the state is contradicted.
# Homogeneity context (verified computationally): each h_f is weighted-
# homogeneous for (d2, d1, d0, e) -> (l**2, l**3, l**4, l**5) and the master
# sum is homogeneous of weight 125 only when c -> l**17 c as well; at the
# fixed c = -1/6630 the scaling is NOT a symmetry, which is why gamma is
# genuinely constrained (through gamma**17 in constant-E states) rather than
# normalizable away.
GAUGE = os.environ.get("BATCH_GAUGE", "0") == "1"
OUT_JSON_GAUGE = ROOT / "batch_convolution_sub2_gauge_raw.json"
# Resume support: skip the first N states of the (deterministic) triage order.
START_INDEX = int(os.environ.get("BATCH_START_INDEX", "0"))
OUT_OVERRIDE = os.environ.get("BATCH_OUT", "")


# ----------------------------------------------------------------------------
# dedup + triage
# ----------------------------------------------------------------------------
def load_unique_states():
    """Return list of dict records, one per unique degree-state tuple.

    Each record carries the dedup key fields plus the list of source cells
    (a_t, b, branch) that collapsed onto it and how many raw states collapsed.
    """
    data = json.loads(WORKLIST.read_text(encoding="utf-8"))
    uniq: dict[tuple, dict] = {}
    for case in data["cases"]:
        is_t2 = case["branch"] == "T2"
        for st in case["states"]:
            key = (
                case["a_t"], is_t2, case["sigma_zero"], case["d2_zero"],
                st["deg_d2"], st["deg_d1"], st["deg_sigma"], st["deg_e"],
            )
            rec = uniq.get(key)
            if rec is None:
                rec = uniq[key] = {
                    "a_t": case["a_t"],
                    "d1_zero": is_t2,
                    "branch": case["branch"],
                    "sigma_zero": case["sigma_zero"],
                    "d2_zero": case["d2_zero"],
                    "deg_d2": st["deg_d2"],
                    "deg_d1": st["deg_d1"],
                    "deg_sigma": st["deg_sigma"],
                    "deg_e": st["deg_e"],
                    "raw_count": 0,
                    "cells": set(),
                }
            rec["raw_count"] += 1
            rec["cells"].add((case["a_t"], tuple(case["b"]), case["branch"]))
    records = list(uniq.values())
    for rec in records:
        rec["cell_count"] = len(rec["cells"])
        rec["cells"] = sorted(
            {(a, "".join(map(str, b)), br) for a, b, br in rec["cells"]}
        )
    return records, data["state_total"]


def tier(rec) -> int:
    if rec["branch"] == "T2":
        return 1
    if rec["deg_e"] == rec["a_t"]:
        return 2
    if rec["sigma_zero"] or rec["d2_zero"]:
        return 3
    return 4


def triage_sort(records):
    return sorted(records, key=lambda r: (tier(r), r["deg_e"] - r["a_t"],
                                          r["a_t"], str(r)))


# ----------------------------------------------------------------------------
# worker: one state, run in a persistent child process (hard timeout by kill)
# ----------------------------------------------------------------------------
def _build_and_descend(rec):
    """Heavy imports live here so the parent stays light and killable."""
    import t5_90t1_verify as base
    import convolution_descent as cd

    y = base.y
    a_t = rec["a_t"]
    deg_e = rec["deg_e"]
    m = deg_e - a_t
    gauge = rec.get("gauge", False)
    parameters = ()
    if gauge:
        # normalize lc of e's generic part to the nonzero parameter gamma
        gamma = sp.Symbol("gamma", nonzero=True)
        gs = tuple(sp.symbols(f"g0:{m}")) + (gamma,)
        parameters = (gamma,)
    else:
        gs = sp.symbols(f"g0:{m + 1}")
    e_expr = (y + 1) ** a_t * sum(g * y ** i for i, g in enumerate(gs))

    zero, degrees = [], {}
    # lc of e's generic part: forceable unknown unless gauge-frozen
    lead = set() if gauge else {gs[-1]}
    if rec["d2_zero"]:
        zero.append("d2")
    else:
        degrees["d2"] = rec["deg_d2"]
        lead.add(sp.Symbol(f"a{rec['deg_d2']}"))
    if rec["d1_zero"]:
        zero.append("d1")
    else:
        degrees["d1"] = rec["deg_d1"]
        lead.add(sp.Symbol(f"b{rec['deg_d1']}"))
    if rec["sigma_zero"]:
        zero.append("sigma")
    else:
        degrees["sigma"] = rec["deg_sigma"]
        lead.add(sp.Symbol(f"s{rec['deg_sigma']}"))

    ansatz = cd.build_ansatz(
        e=e_expr, degrees=degrees, zero=zero,
        sigma=(sp.Integer(0) if rec["sigma_zero"] else None),
        parameters=parameters,
    )
    eng = cd.ConvolutionDescent(ansatz, c=C_VALUE, h=base.load_h())

    # start = 1 + max_f (34 f + (21-3 f) deg_e + maxdeg_y h_f at the state)
    mx = 0
    for f in range(8):
        eng.term_coefficient(f, 0)  # forces h_f expansion into eng._h[f]
        maxh = max(eng._h[f]) if eng._h[f] else 0
        maxphi = max(eng._cached_power(eng._phi_powers, eng.phi, f))
        maxe = max(eng._cached_power(eng._e_powers, eng.e_poly, 21 - 3 * f))
        mx = max(mx, maxphi + maxe + maxh)
    start = mx + 1
    floor = start - FLOOR_BUDGET

    if gauge:
        return _gauge_descend(eng, ansatz, parameters[0], start, floor, lead)

    result = eng.descend(start, floor)

    # Post-process: a FORCED step that sets a leading coeff to 0 == degree drop.
    steps_out = []
    verdict = result.verdict
    stopping = result.stopping_degree
    drop_symbol = None
    for step in result.steps:
        sub = None
        if step.substitution is not None:
            sym, val = step.substitution
            sub = [str(sym), str(val)]
            if sym in lead and sp.simplify(val) == 0:
                # leading coeff forced to zero -> the degree state is contradicted
                steps_out.append({"degree": step.degree,
                                  "verdict": "STATE_KILLED_BY_DEGREE_DROP",
                                  "substitution": sub,
                                  "factored": str(step.factored)})
                verdict = "STATE_KILLED_BY_DEGREE_DROP"
                stopping = step.degree
                drop_symbol = str(sym)
                break
        entry = {"degree": step.degree, "verdict": step.verdict}
        if sub is not None:
            entry["substitution"] = sub
        if step.verdict in ("UNRESOLVED", "CONTRADICTION"):
            entry["factored"] = str(step.factored)
        steps_out.append(entry)

    out = {
        "verdict": verdict,
        "start": start,
        "floor": floor,
        "stopping_degree": stopping,
        "n_steps": len(steps_out),
        "steps": steps_out,
    }
    if drop_symbol is not None:
        out["degree_drop_symbol"] = drop_symbol
    # attach the terminal residual for the interesting / open verdicts
    if result.steps:
        last = result.steps[-1]
        if verdict in ("UNRESOLVED", "CONTRADICTION"):
            out["residual_factored"] = str(last.factored)
    out["substitutions"] = {str(k): str(v)
                            for k, v in result.substitutions.items()}
    return out


def _strip_gamma(poly, gamma):
    """Remove the gamma**k content factor (gamma != 0 is part of the state)."""
    expr = poly.as_expr() if hasattr(poly, "as_expr") else poly
    p = sp.Poly(sp.expand(expr), gamma)
    coeffs = p.all_coeffs()
    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs = coeffs[:-1]
    return sp.Poly(sp.Poly(coeffs, gamma).as_expr(), gamma)


def _gauge_descend(eng, ansatz, gamma, start, floor, lead):
    """Descend with lc(e-generic)=gamma as a nonzero parameter.

    Extension of the landed driver's loop (which stops UNRESOLVED at the
    first equation it cannot force): a residual whose only free symbol is
    gamma is a NECESSARY polynomial constraint on gamma - every master
    coefficient must vanish, and forced substitutions are uniquely
    determined - so such residuals are collected and the descent continues.
    If the running gcd of the collected constraints (with gamma**k stripped;
    gamma != 0 is exactly the degree-exactness of e in the state) becomes a
    nonzero constant, no admissible gamma exists: CONTRADICTION.
    """
    substitutions = {}
    steps_out = []
    constraints = []
    running_gcd = None
    for target in range(start, floor - 1, -1):
        residual = sp.factor(eng.master_coefficient(target).subs(substitutions))
        if residual == 0:
            steps_out.append({"degree": target, "verdict": "IDENTITY"})
            continue
        active = tuple(u for u in eng.ansatz.unknowns
                       if u not in substitutions and residual.has(u))
        if not active:
            if not residual.has(gamma):
                # unknown-free, gamma-free, nonzero exact rational
                steps_out.append({"degree": target, "verdict": "CONTRADICTION",
                                  "factored": str(residual)})
                return _gauge_out("CONTRADICTION", start, floor, target,
                                  steps_out, substitutions, constraints,
                                  running_gcd,
                                  note="nonzero constant coefficient")
            constraint = _strip_gamma(residual, gamma)
            constraints.append(constraint)
            steps_out.append({"degree": target, "verdict": "PARAM_CONSTRAINT",
                              "factored": str(residual),
                              "gamma_constraint": str(constraint.as_expr())})
            running_gcd = (constraint if running_gcd is None
                           else sp.gcd(running_gcd, constraint))
            running_gcd = _strip_gamma(running_gcd, gamma)
            if running_gcd.degree() == 0:
                # no nonzero gamma satisfies all necessary constraints
                steps_out.append({
                    "degree": target, "verdict": "CONTRADICTION",
                    "reason": "gamma constraints have no common nonzero root",
                    "constraints": [str(c.as_expr()) for c in constraints]})
                return _gauge_out("CONTRADICTION", start, floor, target,
                                  steps_out, substitutions, constraints,
                                  running_gcd,
                                  note="incompatible gamma constraints")
            continue
        forced = eng._forced_square(residual, active)
        if forced:
            sym, val = forced
            sub = [str(sym), str(val)]
            if sym in lead and sp.simplify(val) == 0:
                steps_out.append({"degree": target,
                                  "verdict": "STATE_KILLED_BY_DEGREE_DROP",
                                  "substitution": sub,
                                  "factored": str(residual)})
                out = _gauge_out("STATE_KILLED_BY_DEGREE_DROP", start, floor,
                                 target, steps_out, substitutions, constraints,
                                 running_gcd)
                out["degree_drop_symbol"] = str(sym)
                return out
            substitutions[sym] = val
            steps_out.append({"degree": target, "verdict": "FORCED",
                              "substitution": sub})
            continue
        steps_out.append({"degree": target, "verdict": "UNRESOLVED",
                          "factored": str(residual)})
        out = _gauge_out("UNRESOLVED", start, floor, target, steps_out,
                         substitutions, constraints, running_gcd)
        out["residual_factored"] = str(residual)
        return out
    return _gauge_out("FORCED", start, floor, floor, steps_out,
                      substitutions, constraints, running_gcd)


def _gauge_out(verdict, start, floor, stopping, steps_out, substitutions,
               constraints, running_gcd, note=None):
    out = {
        "verdict": verdict,
        "start": start,
        "floor": floor,
        "stopping_degree": stopping,
        "n_steps": len(steps_out),
        "steps": steps_out,
        "substitutions": {str(k): str(v) for k, v in substitutions.items()},
        "gamma_constraint_count": len(constraints),
    }
    if note:
        out["note"] = note
    if running_gcd is not None:
        out["gamma_constraint_gcd"] = str(running_gcd.as_expr())
    return out


def _worker(inbox, outbox):
    while True:
        item = inbox.get()
        if item is None:
            return
        try:
            outbox.put(("ok", _build_and_descend(item)))
        except Exception as exc:  # noqa: BLE001 - report, never hang the batch
            import traceback
            outbox.put(("err", f"{type(exc).__name__}: {exc}\n"
                               f"{traceback.format_exc()}"))


class Runner:
    """Persistent child process; killed and respawned only on timeout."""

    def __init__(self):
        self.ctx = mp.get_context("spawn")
        self._spawn()

    def _spawn(self):
        self.inbox = self.ctx.Queue()
        self.outbox = self.ctx.Queue()
        self.proc = self.ctx.Process(target=_worker,
                                     args=(self.inbox, self.outbox),
                                     daemon=True)
        self.proc.start()

    def run(self, rec, timeout):
        self.inbox.put(rec)
        try:
            status, payload = self.outbox.get(timeout=timeout)
        except Exception:  # queue.Empty -> timed out
            self.proc.terminate()
            self.proc.join()
            self._spawn()
            return ("timeout", None)
        return (status, payload)

    def close(self):
        try:
            self.inbox.put(None)
            self.proc.join(timeout=5)
        except Exception:
            pass
        if self.proc.is_alive():
            self.proc.terminate()


# ----------------------------------------------------------------------------
def main():
    records, raw_total = load_unique_states()
    ordered = triage_sort(records)
    if GAUGE:
        for rec in ordered:
            rec["gauge"] = True
    tier_counts = {t: sum(1 for r in ordered if tier(r) == t)
                   for t in (1, 2, 3, 4)}
    print(f"raw states {raw_total}; unique {len(ordered)}; "
          f"tiers {tier_counts}; gauge={GAUGE}")

    runner = Runner()
    attempted = []
    census = {}
    t_start = time.time()
    stop_reason = "completed all states"
    stopped_index = len(ordered)

    for idx, rec in enumerate(ordered):
        if idx < START_INDEX:
            continue
        if MAX_STATES and idx >= MAX_STATES:
            stop_reason = f"MAX_STATES cap ({MAX_STATES})"
            stopped_index = idx
            break
        elapsed = time.time() - t_start
        if elapsed >= TOTAL_WALL_BUDGET:
            stop_reason = "total wall budget exhausted"
            stopped_index = idx
            break
        remaining = TOTAL_WALL_BUDGET - elapsed
        timeout = min(PER_STATE_TIMEOUT, remaining)

        t0 = time.time()
        status, payload = runner.run(rec, timeout)
        dt = round(time.time() - t0, 1)

        entry = {k: rec[k] for k in (
            "a_t", "branch", "d1_zero", "sigma_zero", "d2_zero",
            "deg_d2", "deg_d1", "deg_sigma", "deg_e",
            "raw_count", "cell_count", "cells")}
        entry["tier"] = tier(rec)
        entry["m"] = rec["deg_e"] - rec["a_t"]
        entry["gauge"] = bool(rec.get("gauge", False))
        entry["seconds"] = dt

        if status == "ok":
            entry.update(payload)
        elif status == "timeout":
            entry["verdict"] = "SKIPPED_BUDGET"
            entry["reason"] = f"exceeded per-state timeout {timeout:.0f}s"
        else:  # err
            entry["verdict"] = "ERROR"
            entry["error"] = payload

        census[entry["verdict"]] = census.get(entry["verdict"], 0) + 1
        attempted.append(entry)

        if entry["verdict"] in ("CONTRADICTION", "STATE_KILLED_BY_DEGREE_DROP",
                                "FORCED"):
            print(f"[{idx}] tier{entry['tier']} a_t={rec['a_t']} "
                  f"deg_e={rec['deg_e']} -> {entry['verdict']} "
                  f"(stop {entry.get('stopping_degree')}, {dt}s)")
        else:
            print(f"[{idx}] tier{entry['tier']} a_t={rec['a_t']} "
                  f"deg_e={rec['deg_e']} m={entry['m']} -> "
                  f"{entry['verdict']} ({dt}s)")

    runner.close()

    payload = {
        "schema": 1,
        "source_worklist": WORKLIST.name,
        "gauge_normalized": GAUGE,
        "c": str(C_VALUE),
        "floor_budget": FLOOR_BUDGET,
        "per_state_timeout_s": PER_STATE_TIMEOUT,
        "total_wall_budget_s": TOTAL_WALL_BUDGET,
        "raw_state_total": raw_total,
        "unique_state_total": len(ordered),
        "tier_counts": tier_counts,
        "attempted_count": len(attempted),
        "start_index": START_INDEX,
        "stopped_at_index": stopped_index,
        "stop_reason": stop_reason,
        "verdict_census": census,
        "ansatz_note": (
            "q-root support conditions on e are DROPPED (sound "
            "over-approximation): a kill here kills the original flag case; an "
            "UNRESOLVED here does not certify the original survives."),
        "states": attempted,
    }
    out_path = (ROOT / OUT_OVERRIDE if OUT_OVERRIDE
                else (OUT_JSON_GAUGE if GAUGE else OUT_JSON))
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8")
    print(f"\ncensus {census}")
    print(f"attempted {len(attempted)} / {len(ordered)} "
          f"(stopped at index {stopped_index}: {stop_reason})")
    print(f"-> {out_path.name}")


def merge():
    """Merge the raw pass artifacts into the final batch_convolution_sub2.json.

    Passes (in precedence order for the final per-state verdict):
      1. gauge + gauge_resume (decisive: lc(e) nonzero-parameter + constraint
         accumulation) - disjoint index ranges of the same triage order;
      2. ungauged first pass (kept as a secondary record; its blanket
         UNRESOLVED verdicts reflect the free e-scale coupling every
         top-degree equation, not evidence of survival).
    """
    def state_key(s):
        return (s["a_t"], s["branch"], s["sigma_zero"], s["d2_zero"],
                str(s["deg_d2"]), str(s["deg_d1"]), str(s["deg_sigma"]),
                s["deg_e"])

    passes = []
    for name, fname in (
        ("ungauged", "batch_convolution_sub2_pass1_ungauged.json"),
        ("gauge", "batch_convolution_sub2_gauge_raw.json"),
        ("gauge_resume", "batch_convolution_sub2_gauge_resume.json"),
    ):
        path = ROOT / fname
        if path.exists():
            passes.append((name, json.loads(path.read_text(encoding="utf-8"))))

    merged: dict[tuple, dict] = {}
    for name, data in passes:
        for s in data["states"]:
            key = state_key(s)
            rec = merged.setdefault(key, {
                k: s[k] for k in (
                    "a_t", "branch", "d1_zero", "sigma_zero", "d2_zero",
                    "deg_d2", "deg_d1", "deg_sigma", "deg_e",
                    "raw_count", "cell_count", "cells", "tier", "m")})
            rec[f"{name}_verdict"] = s["verdict"]
            if name.startswith("gauge"):
                rec["gauge_detail"] = {
                    k: s[k] for k in (
                        "verdict", "start", "floor", "stopping_degree",
                        "steps", "substitutions", "seconds")
                    if k in s}
                for extra in ("residual_factored", "gamma_constraint_count",
                              "gamma_constraint_gcd", "degree_drop_symbol",
                              "note", "reason", "error"):
                    if extra in s:
                        rec["gauge_detail"][extra] = s[extra]

    for rec in merged.values():
        gauge_v = rec.get("gauge_verdict") or rec.get("gauge_resume_verdict")
        rec["final_verdict"] = gauge_v or rec.get("ungauged_verdict")
        rec["decided_by"] = ("gauge" if gauge_v else "ungauged_only")

    states = sorted(merged.values(),
                    key=lambda r: (r["tier"], r["m"], r["a_t"],
                                   str(r["deg_d2"]), str(r["deg_d1"]),
                                   str(r["deg_sigma"]), r["deg_e"]))
    census = {}
    for rec in states:
        census[rec["final_verdict"]] = census.get(rec["final_verdict"], 0) + 1
    kills = [r for r in states if r["final_verdict"] in
             ("CONTRADICTION", "STATE_KILLED_BY_DEGREE_DROP")]

    records, raw_total = load_unique_states()
    ordered = triage_sort(records)
    attempted_keys = set(merged)
    unattempted = [r for r in ordered
                   if (r["a_t"], r["branch"], r["sigma_zero"], r["d2_zero"],
                       str(r["deg_d2"]), str(r["deg_d1"]), str(r["deg_sigma"]),
                       r["deg_e"]) not in attempted_keys]
    tier_counts_all = {t: sum(1 for r in ordered if tier(r) == t)
                       for t in (1, 2, 3, 4)}
    tier_counts_open = {t: sum(1 for r in unattempted if tier(r) == t)
                        for t in (1, 2, 3, 4)}

    payload = {
        "schema": 2,
        "description": (
            "Convolution-descent verdicts over the deduped sub2 Phase-D "
            "residual degree states (master identity f31; c=-1/6630; "
            "q-root support conditions dropped - sound over-approximation; "
            "kills PENDING AUDIT)"),
        "source_worklist": WORKLIST.name,
        "c": str(C_VALUE),
        "floor_budget": FLOOR_BUDGET,
        "raw_state_total": raw_total,
        "unique_state_total": len(ordered),
        "tier_counts": tier_counts_all,
        "passes": [
            {"name": name,
             "gauge_normalized": bool(data.get("gauge_normalized")),
             "attempted": data["attempted_count"],
             "start_index": data.get("start_index", 0),
             "stopped_at_index": data["stopped_at_index"],
             "stop_reason": data["stop_reason"],
             "census": data["verdict_census"],
             "wall_budget_s": data["total_wall_budget_s"]}
            for name, data in passes],
        "attempted_unique": len(states),
        "unattempted_unique": len(unattempted),
        "unattempted_by_tier": tier_counts_open,
        "final_verdict_census": census,
        "kill_count": len(kills),
        "kill_raw_state_coverage": sum(r["raw_count"] for r in kills),
        "kills_pending_audit": kills,
        "states": states,
    }
    OUT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8")
    print(f"final census {census}")
    print(f"kills {len(kills)} (raw coverage "
          f"{sum(r['raw_count'] for r in kills)})")
    print(f"attempted {len(states)} unique; unattempted {len(unattempted)} "
          f"by tier {tier_counts_open}")
    print(f"-> {OUT_JSON.name}")


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "merge":
        merge()
    else:
        main()
