#!/usr/bin/env python3
"""Alternate-regime max-plus DEGREE sweep over the 27 open flipped-cascade
branches (companion: ALT_INF_SWEEP.md; derivation: ALT_REGIME_INF.md).

For each of the 27 branches left OPEN by ALT_REGIME_L2.md section 4 (13 T1 +
14 T2, strata a in [11,14]) this enumerates every degree assignment
(deg d2, deg d1, deg sigma, deg e) inside the sub1 sandwich, runs the FLIPPED
descending max-plus chain (D_t) of ALT_REGIME_INF.md (b), and decides whether
the bottom close  E^21 h_0 + u r_0 = 0  can be satisfied (a TIE at close) or is
refuted (a strict/unique maximum => nonzero leading term => contradiction).

Semantics mirror the STANDARD-regime infinity layer of cascade_engine.py
(deg_h_options / descend_options_inf / inf_place_profiles); that module is
IMPORTED (never edited) so the deg_h drop rules, obligation objects, monomial
tables (cascade_signature.load_levels), exponent order (d2,d1,sigma,e) and
DEG_U = 4 are the exact same code the audited standard layer uses.

Flipped chain identities (ALT_REGIME_INF.md (I7)/(If)/(I0)):

    (I7)  top anchor  f=7 :  w + R_6      = H_7                       (unique)
    (If)  levels    f=6..1 :  w + R_{f-1} = max( 3(7-f) deg_E + H_f, 4 + R_f )
    (I0)  bottom close f=0 :  max( 21 deg_E + H_0, 4 + R_0 ) must TIE

with w = 3a-30 = deg T, deg_E = deg e - a, H_f = deg h_f, R_f = deg r_f.
A drop below a max is allowed ONLY on a tie, recorded as an Obligation
(leading_cancellation / degree_tie_drop / identical_vanishing) exactly as the
standard layer records it.  Exact integer arithmetic throughout.

Zero-flag treatment (ALT_REGIME_L2.md section 2 case analysis):
  * branch T2  <=>  d1 == 0 identically  =>  h_7 == 0  =>  r_6 == 0
    (ALT_REGIME.md "Descending cascade", T2 paragraph): the T2 chain top is the
    level-6 identity  T r_5 = E^3 h_6, which the general (If) recursion realises
    automatically once R_6 = deg r_6 = NEG_INF is seeded (4 + R_6 = NEG_INF
    drops out, leaving the unique h-term).  We follow that, not a guess.
  * sigma == 0 and d2 == 0 flags are ENUMERATED as separate degree states
    (deg = NEG_INF).  Over-approximation is sound: we never silently restrict.
  * T3 (d1 == 0 AND sigma == 0) is excluded globally (split_place_ledger_sub1
    "T3 ... excluded globally"); those states are skipped, not counted.
  * r_f == 0 flags: no independent cap on deg r_f exists (ALT_REGIME_INF.md
    [judgment]); r_f == 0 (R_f = NEG_INF) is reachable only through a recorded
    identical-vanishing obligation on the level above, never assumed.
"""
from __future__ import annotations

import json
import time
from dataclasses import asdict
from pathlib import Path

import cascade_engine as ce

ROOT = Path(__file__).resolve().parent
NEG_INF = ce.NEG_INF            # float('-inf'); encodes "identically zero"
DEG_U = ce.DEG_U                # == 4
Obligation = ce.Obligation
deg_h_options = ce.deg_h_options

# sub1 stripped-window caps (ALT_REGIME_INF.md (c); sub1_cascade_verify.py).
CAP_D2, CAP_D1, CAP_SIGMA, CAP_E = 6, 9, 12, 15

# ---------------------------------------------------------------------------
# The 27 open branches: ALT_REGIME_L2.md section 4 "Per-input-branch verdict"
# (the O rows after the six h6/h5 kills).  13 T1 + 14 T2.
# ---------------------------------------------------------------------------
OPEN_BRANCHES = [
    # a, sorted b-vector, branch      (source: ALT_REGIME_L2.md sec.4)
    (11, (0, 0, 0, 0), "T1"), (11, (1, 0, 0, 0), "T1"),
    (11, (1, 1, 0, 0), "T1"), (11, (1, 1, 1, 0), "T1"),
    (11, (1, 1, 1, 1), "T1"), (11, (3, 0, 0, 0), "T1"),
    (12, (0, 0, 0, 0), "T1"), (12, (1, 0, 0, 0), "T1"),
    (12, (1, 1, 0, 0), "T1"), (12, (1, 1, 1, 0), "T1"),
    (12, (3, 0, 0, 0), "T1"),
    (14, (0, 0, 0, 0), "T1"), (14, (1, 0, 0, 0), "T1"),
    (11, (0, 0, 0, 0), "T2"), (11, (1, 0, 0, 0), "T2"),
    (11, (1, 1, 0, 0), "T2"), (11, (1, 1, 1, 0), "T2"),
    (11, (1, 1, 1, 1), "T2"), (11, (3, 0, 0, 0), "T2"),
    (11, (3, 1, 0, 0), "T2"),
    (12, (0, 0, 0, 0), "T2"), (12, (1, 0, 0, 0), "T2"),
    (12, (1, 1, 0, 0), "T2"), (12, (1, 1, 1, 0), "T2"),
    (13, (0, 0, 0, 0), "T2"), (13, (1, 0, 0, 0), "T2"),
    (14, (0, 0, 0, 0), "T2"),
]
assert len([b for b in OPEN_BRANCHES if b[2] == "T1"]) == 13
assert len([b for b in OPEN_BRANCHES if b[2] == "T2"]) == 14


def branch_id(a: int, b: tuple, branch: str) -> str:
    return f"a{a}_b{''.join(map(str, b))}_{branch}"


# ---------------------------------------------------------------------------
# max-plus helpers.  Degrees are ints or NEG_INF (identically-zero poly).
# ---------------------------------------------------------------------------
def _tie_label(f: int, term_kind: str) -> tuple:
    return (f"E^{3*(7-f)} h_{f}  (g-side)", f"u r_{f}  (h-side)")


def transition(f: int, term1: float, term2: float, w: int):
    """One (If) step.  Returns list of (R_prev, extra_obligations).

    R_prev is an int (deg r_{f-1} >= 0) or NEG_INF (r_{f-1} == 0).  A drop
    below the max is emitted ONLY when the two terms tie; each drop carries a
    leading_cancellation obligation (identical_vanishing for full collapse),
    matching descend_options_inf.  A forced negative degree is impossible
    (no candidate) -- a dead path.
    """
    tie = term1 == term2
    mx = term1 if term1 >= term2 else term2  # NEG_INF-safe (both may be -inf)
    out = []
    if mx == NEG_INF:
        # both terms identically zero -> r_{f-1} == 0, forced, no obligation.
        return [(NEG_INF, ())]
    base = mx - w
    if not tie:
        # unique maximum: forced, no leading cancellation available.
        if base >= 0:
            out.append((int(base), ()))
        # base < 0 with a unique nonzero leading term: deg r_{f-1} < 0 for a
        # nonzero polynomial -- impossible.  No candidate (dead path).
        return out
    # tie: the two equal-degree leading forms may add (no drop) or cancel.
    tied = _tie_label(f, "sum")
    if base >= 0:
        out.append((int(base), ()))                       # generic: no cancel
        for k in range(0, int(base)):                     # cancel to depth base-k
            out.append((k, (Obligation(f, "leading_cancellation",
                                       int(base) - k, tied),)))
    # full collapse r_{f-1} == 0 (needs the whole sum to vanish).
    out.append((NEG_INF, (Obligation(f, "identical_vanishing", 0, tied),)))
    return out


def _merge(reach: dict, key, obl: tuple) -> None:
    """Keep the fewest-obligation witness per reachable R value."""
    cur = reach.get(key)
    if cur is None or len(obl) < len(cur):
        reach[key] = obl


_CHAIN_CACHE: dict = {}


def run_chain(a: int, degstate: tuple):
    """Run (D_t) for one degree state.  Returns a result dict.

    degstate = (deg_d2, deg_d1, deg_sigma, deg_e), each int or NEG_INF.
    Result keys: verdict ('survive'/'killed'), reason (if killed),
    obligations (minimal witness list, if survive), close (tie value + terms).
    """
    key = (a, degstate)
    cached = _CHAIN_CACHE.get(key)
    if cached is not None:
        return cached

    w = 3 * a - 30
    de = degstate[3]
    deg_E = de - a                      # e never zero: de finite, deg_E in [0..4]
    flags = (degstate[2] == NEG_INF, degstate[0] == NEG_INF,
             degstate[1] == NEG_INF)   # (sigma_zero, d2_zero, d1_zero)

    # (I7) top anchor  T r_6 = h_7.
    h7opts = deg_h_options(7, degstate, flags)
    h7val, h7ob = h7opts[0]             # h_7 = 8192 d1^2: single monomial, forced
    reach: dict = {}
    if h7val == NEG_INF:
        # T2: h_7 == 0 => r_6 == 0.  Seed R_6 = NEG_INF (4+R_6 drops out at f=6).
        reach[NEG_INF] = ()
    else:
        base = h7val - w               # deg r_6 = H_7 - w, forced
        if base < 0:
            res = {"verdict": "killed",
                   "reason": f"top anchor forces deg r_6 = H_7 - w = "
                             f"{h7val} - {w} = {base} < 0 (r_6 must be a nonzero "
                             f"polynomial): 2*deg d1 = {h7val} < w = {w}",
                   "R6": base}
            _CHAIN_CACHE[key] = res
            return res
        reach[int(base)] = tuple(h7ob)

    # (If) levels f = 6 .. 1.
    for f in range(6, 0, -1):
        hopts = deg_h_options(f, degstate, flags)
        gshift = 3 * (7 - f) * deg_E
        newreach: dict = {}
        for Rf, ob in reach.items():
            term2 = NEG_INF if Rf == NEG_INF else DEG_U + Rf
            for hval, hob in hopts:
                term1 = NEG_INF if hval == NEG_INF else gshift + hval
                for Rp, extra in transition(f, term1, term2, w):
                    _merge(newreach, Rp, tuple(ob) + tuple(hob) + tuple(extra))
        if not newreach:
            res = {"verdict": "killed",
                   "reason": f"level f={f}: every branch forces deg r_{f-1} < 0 "
                             f"(unique maximum, nonzero leading term) -- no "
                             f"consistent deg r_{f-1}",
                   "dead_level": f}
            _CHAIN_CACHE[key] = res
            return res
        reach = newreach

    # (I0) bottom close  E^21 h_0 + u r_0 = 0  must be a TIE.
    h0opts = deg_h_options(0, degstate, flags)
    best = None
    close_terms = None
    for R0, ob in reach.items():
        term2 = NEG_INF if R0 == NEG_INF else DEG_U + R0
        for hval, hob in h0opts:
            term1 = NEG_INF if hval == NEG_INF else 21 * deg_E + hval
            if term1 == term2:                       # TIE -> close can vanish
                if term1 == NEG_INF:
                    extra = ()                        # 0 = 0 (h_0 == 0, r_0 == 0)
                else:
                    extra = (Obligation(0, "leading_cancellation", 0,
                                        (f"E^21 h_0  (g-side)", "u r_0 (h-side)")),)
                total = tuple(ob) + tuple(hob) + tuple(extra)
                if best is None or len(total) < len(best):
                    best = total
                    close_terms = {"term1_21degE+H0": _num(term1),
                                   "term2_4+R0": _num(term2),
                                   "R0": _num(R0),
                                   "H0_used": _num(hval)}
    if best is None:
        # No reachable (R_0, H_0) makes the close a tie: unique maximum at the
        # closing anchor for EVERY consistent chain -> nonzero leading term.
        r0vals = sorted((_num(r) for r in reach),
                        key=lambda v: (v is None, v))
        res = {"verdict": "killed",
               "reason": "bottom close E^21 h_0 + u r_0 has a unique maximum for "
                         "every reachable degree chain (21 deg_E + H_0 never "
                         "equals 4 + R_0), so its leading term is nonzero and the "
                         "sum cannot be 0",
               "reachable_R0": r0vals}
        _CHAIN_CACHE[key] = res
        return res

    res = {"verdict": "survive",
           "obligations": [asdict(o) for o in best],
           "n_obligations": len(best),
           "close": close_terms}
    _CHAIN_CACHE[key] = res
    return res


def _num(x):
    return None if x == NEG_INF else int(x)


# ---------------------------------------------------------------------------
# Degree-state enumeration inside the sub1 sandwich for one branch.
# ---------------------------------------------------------------------------
def enumerate_states(a: int, b: tuple, branch: str):
    """Yield degstate tuples (deg_d2, deg_d1, deg_sigma, deg_e).

    deg d2 in {NEG_INF(0-flag), 0..6}; deg sigma in {NEG_INF, 0..12};
    deg d1 = NEG_INF for T2 else 0..9; deg e in [a+sum b .. 15]
    (deg_E = deg e - a >= sum b since E = prod p_i^{b_i} F, deg F >= 0).
    T3 (d1 == 0 and sigma == 0) is skipped (excluded globally).
    """
    sumb = sum(b)
    d1_zero = branch == "T2"
    d1_opts = [NEG_INF] if d1_zero else list(range(0, CAP_D1 + 1))
    d2_opts = [NEG_INF] + list(range(0, CAP_D2 + 1))
    sigma_opts = [NEG_INF] + list(range(0, CAP_SIGMA + 1))
    e_lo = a + sumb
    e_opts = list(range(e_lo, CAP_E + 1))       # deg e; deg_E = deg e - a
    for de in e_opts:
        for dd1 in d1_opts:
            for dd2 in d2_opts:
                for dsig in sigma_opts:
                    if d1_zero and dsig == NEG_INF:
                        continue                # T3: excluded globally
                    yield (dd2, dd1, dsig, de)


def sweep_branch(a: int, b: tuple, branch: str, sample_cap: int = 25):
    total = surviving = killed = 0
    survive_states = []          # compact: [dd2,dd1,dsig,de,n_obl]
    killed_states = []           # compact: [dd2,dd1,dsig,de]
    survive_samples = []         # full obligations for a few witnesses
    killed_samples = []          # full reason for a few kills
    kill_reasons: dict = {}
    for st in enumerate_states(a, b, branch):
        total += 1
        res = run_chain(a, st)
        row = [_num(st[0]), _num(st[1]), _num(st[2]), _num(st[3])]
        if res["verdict"] == "survive":
            surviving += 1
            survive_states.append(row + [res["n_obligations"]])
            if len(survive_samples) < sample_cap:
                survive_samples.append(
                    {"degstate": {"deg_d2": row[0], "deg_d1": row[1],
                                  "deg_sigma": row[2], "deg_e": row[3],
                                  "deg_E": row[3] - a},
                     "close": res["close"],
                     "obligations": res["obligations"]})
        else:
            killed += 1
            killed_states.append(row)
            tag = res["reason"].split(":")[0]
            kill_reasons[tag] = kill_reasons.get(tag, 0) + 1
            if len(killed_samples) < sample_cap:
                killed_samples.append(
                    {"degstate": {"deg_d2": row[0], "deg_d1": row[1],
                                  "deg_sigma": row[2], "deg_e": row[3],
                                  "deg_E": row[3] - a},
                     "reason": res["reason"],
                     "detail": {k: v for k, v in res.items()
                                if k not in ("verdict", "reason")}})
    verdict = "OPEN" if surviving > 0 else "KILLED"
    return {
        "id": branch_id(a, b, branch),
        "a": a, "b": list(b), "sum_b": sum(b), "branch": branch,
        "w": 3 * a - 30, "deg_E_range": [sum(b), 15 - a],
        "verdict": verdict,
        "counts": {"total_degree_states": total,
                   "surviving": surviving, "killed": killed},
        "kill_reason_histogram": kill_reasons,
        "survive_samples": survive_samples,
        "killed_samples": killed_samples,
        "surviving_states_compact": survive_states,
        "killed_states_compact": killed_states,
    }


def main() -> None:
    t0 = time.time()
    results = []
    for a, b, branch in OPEN_BRANCHES:
        results.append(sweep_branch(a, b, branch))
    elapsed = time.time() - t0

    n_open = sum(r["verdict"] == "OPEN" for r in results)
    n_killed = sum(r["verdict"] == "KILLED" for r in results)
    tot_states = sum(r["counts"]["total_degree_states"] for r in results)
    tot_surv = sum(r["counts"]["surviving"] for r in results)
    tot_kill = sum(r["counts"]["killed"] for r in results)

    out = {
        "schema": {
            "version": 1,
            "description": "Alternate-regime max-plus degree sweep over the 27 "
                           "open flipped-cascade branches.",
            "derivation": "ALT_REGIME_INF.md",
            "branch_source": "ALT_REGIME_L2.md section 4 (27 O-rows: 13 T1 + 14 T2)",
            "semantics": "cascade_engine.py deg_h_options / descend_options_inf "
                         "(imported, not edited); exponent order (d2,d1,sigma,e); "
                         "DEG_U=4; NEG_INF encodes identically-zero.",
            "identities": {
                "I7": "w + deg r_6 = deg h_7 (unique, forced)",
                "If": "w + deg r_{f-1} = max(3(7-f) deg_E + deg h_f, 4 + deg r_f)",
                "I0": "max(21 deg_E + deg h_0, 4 + deg r_0) must TIE else contradiction",
            },
            "sandwich": {"deg_d2<=": CAP_D2, "deg_d1<=": CAP_D1,
                         "deg_sigma<=": CAP_SIGMA, "deg_e<=": CAP_E,
                         "deg_e>=": "a + sum(b)"},
            "state_fields": "surviving_states_compact rows = "
                            "[deg_d2, deg_d1, deg_sigma, deg_e, n_obligations]; "
                            "killed rows = [deg_d2, deg_d1, deg_sigma, deg_e]; "
                            "null = NEG_INF (identically zero).",
        },
        "summary": {
            "n_branches": len(results),
            "branches_OPEN": n_open,
            "branches_KILLED": n_killed,
            "total_degree_states": tot_states,
            "surviving_states": tot_surv,
            "killed_states": tot_kill,
            "elapsed_seconds": round(elapsed, 2),
        },
        "branches": results,
    }
    (ROOT / "alt_inf_sweep.json").write_text(
        json.dumps(out, indent=1), encoding="utf-8")

    # ---- summary table -----------------------------------------------------
    print(f"Alternate-regime max-plus degree sweep -- {len(results)} branches "
          f"({elapsed:.1f}s)\n")
    hdr = f"{'branch id':22} {'a':>2} {'sum_b':>5} {'br':>3} {'verdict':>7} " \
          f"{'states':>7} {'surv':>7} {'killed':>7}"
    print(hdr)
    print("-" * len(hdr))
    for r in results:
        c = r["counts"]
        print(f"{r['id']:22} {r['a']:>2} {r['sum_b']:>5} {r['branch']:>3} "
              f"{r['verdict']:>7} {c['total_degree_states']:>7} "
              f"{c['surviving']:>7} {c['killed']:>7}")
    print("-" * len(hdr))
    print(f"OPEN={n_open}  KILLED={n_killed}  | degree states: total={tot_states} "
          f"surviving={tot_surv} killed={tot_kill}")
    print("\nWrote alt_inf_sweep.json")


if __name__ == "__main__":
    main()
