#!/usr/bin/env python3
"""r9_valsplit.py -- the divisibility-forced VALUATION SPLIT for the
deg_e = 10 cost wall: the named continuation after the dm4 elimination
(R9_SYMBOLIC.md section 3c).

THE SPLIT (sound covering of the variety; exhaustiveness CHECKED below)
-----------------------------------------------------------------------
On the G-variety the certified divisibility lemma gives
monic(e) | dm2*dm3 in Q[y] (r9_symbolic_elim.py, identity exact).  With

    R9 column:     monic(e) = (y+1)^9 * (y - r)        (e fully known)
    batch a_t=a:   monic(e) = (y+1)^a                  (known factor only)

prime-factor valuation additivity forces, for R9,

    v_{y+1}(dm2) + v_{y+1}(dm3) >= 9   and   v_{y-r}(dm2) + v_{y-r}(dm3) >= 1,

so every point of the variety (including dm2 == 0 or dm3 == 0, which any
case represents with A = 0 / B = 0) satisfies, for SOME i in 0..9,
j in 0..1,

    dm2 = (y+1)^i     * (y-r)^j     * A,   deg A <= 12 - i - j
    dm3 = (y+1)^(9-i) * (y-r)^(1-j) * B,   deg B <= 14 - (9-i) - (1-j)

-- 20 cases, 18 spare unknowns each (13-i-j + 5+i+j), down from 28.  The
batch analogue with (y+1)^a is i = 0..a: a+1 cases of 28-a unknowns.
Implementation: each case is imposed as LINEAR valuation equations on the
cached generic dm4-eliminated build (see case_equations) -- equivalent to
the structured product ansatz, with zero re-expansion cost.  The base
build's divisibility remainder rows stay in the system (redundant under
the split, still sound); H2,H3,H5, q(r) = 0 and saturations compose
unchanged.

EXHAUSTIVENESS is a finite claim over valuation profiles and is CHECKED
exhaustively at import (check_exhaustive): for every achievable profile
(alpha, beta, rho2, rho3) with alpha+beta >= n1, rho2+rho3 >= n2 and the
degree caps, some case (i, j) covers it (i <= alpha, n1-i <= beta,
j <= rho2, n2-j <= rho3).  A state is KILLED only if ALL its cases are
exact-Q UNIT (each case verdict is a sound necessary condition under its
case hypothesis; the union of case varieties contains the state variety).

VERDICT SEMANTICS per case: UNIT => case dead; PROPER => inconclusive
(case system still weaker than the full bridge: dm4 polynomiality-of-
capped-degree is not imposed); TIMEOUT => COST.  Per state: KILLED iff
all cases dead; DENTED if some die; COST/WALL otherwise.

Machinery reuse: r9_symbolic_sweep's per-term convolution builder
(_g_coeff_lists, validated byte-identical vs fsb.augment on R9 z=1), its
orphan-proof WSL runner (timeout + ulimit -v 8G inside WSL, installed
over bridge_sweep._run at import), bridge_sweep triage/exact routines
verbatim.  New file; READ-ONLY on every landed module and on
r9_symbolic_sweep.json (this lane records to r9_valsplit_results.json).

Kills are PENDING AUDIT; each killed case stores the integer-cleared
generators + saturation factors (certificate-extractable).
"""
from __future__ import annotations

import json
import pickle
import sys
import time
from pathlib import Path

import sympy as sp

import full_system_bridge as fsb
import bridge_sweep as bsw
import convolution_elim_qsupport as qs
import r9_symbolic_sweep as rss          # installs the guarded runner

ROOT = Path(__file__).resolve().parent
CACHE = ROOT / "__pycache__"
OUT = ROOT / "r9_valsplit_results.json"

y = fsb.y
R = sp.Symbol("r")
CAPS = fsb.STRIP_DEGCAP["sub2"]          # dm2 <= 12, dm3 <= 14


# --------------------------------------------------------------------------
#  case family + checked exhaustiveness
# --------------------------------------------------------------------------
def cases_for(n1: int, n2: int):
    """All (i, j) with 0 <= i <= n1, 0 <= j <= n2 (j absent when n2 = 0)."""
    return [(i, j) for i in range(n1 + 1) for j in range(n2 + 1)]


def check_exhaustive(n1: int, n2: int, cap2: int, cap3: int) -> int:
    """Finite check: every achievable valuation profile is covered by some
    case.  Profiles: alpha = v_{y+1}(dm2) in 0..cap2, beta = v_{y+1}(dm3)
    in 0..cap3 with alpha + beta >= n1; rho2, rho3 likewise for the (y-r)
    factor with rho2 + rho3 >= n2 and alpha + rho2 <= cap2 (degree room),
    beta + rho3 <= cap3.  Returns the number of profiles checked."""
    cases = cases_for(n1, n2)
    n_checked = 0
    for alpha in range(cap2 + 1):
        for beta in range(cap3 + 1):
            if alpha + beta < n1:
                continue
            for rho2 in range(cap2 - alpha + 1):
                for rho3 in range(cap3 - beta + 1):
                    if rho2 + rho3 < n2:
                        continue
                    n_checked += 1
                    if not any(i <= alpha and n1 - i <= beta
                               and j <= rho2 and n2 - j <= rho3
                               for i, j in cases):
                        raise AssertionError(
                            f"profile (a={alpha},b={beta},r2={rho2},"
                            f"r3={rho3}) uncovered for n1={n1},n2={n2}")
    return n_checked


# --------------------------------------------------------------------------
#  case instantiation: linear valuation equations on the CACHED generic build
# --------------------------------------------------------------------------
def case_equations(i: int, j: int, n1: int, n2: int):
    """The split case as LINEAR equations on the GENERIC spare coefficients
    (the exact symbols of r9_symbolic_sweep._spare_polys).

    Over char 0, v_{y+1}(p) >= m  <=>  p(-1) = p'(-1) = ... = p^(m-1)(-1)=0,
    and v_{y-r}(p) >= 1  <=>  p(r) = 0.  So case (i, j) is

        dm2^(k)(-1) = 0 for k < i,      dm3^(k)(-1) = 0 for k < n1 - i,
        dm2(r) = 0 if j = 1  /  dm3(r) = 0 if n2 = 1 and j = 0.

    Same case variety as the structured-product ansatz dm2 = (y+1)^i(y-r)^j A
    (A, B are the generic coefficients modulo this linear system), but
    instantiation is FREE: it composes with the CACHED dm4-eliminated build
    (r9red_*.pkl) instead of re-running the symbolic expansion -- the product
    form's binomial-spread coefficients densified every convolution and took
    >25 min/case (measured on z=1 i=5, killed).  Degree-cap bookkeeping is
    unchanged (generic caps 12/14); the linear rows cut the spare dimension
    28 -> 18 (R9) / 28-a (batch), and Singular eliminates them instantly."""
    rs = sp.symbols(f"R0:{CAPS['dm2'] + 1}")
    ss = sp.symbols(f"S0:{CAPS['dm3'] + 1}")
    dm2 = sum(c * y ** k for k, c in enumerate(rs))
    dm3 = sum(c * y ** k for k, c in enumerate(ss))
    eqs = []
    for k in range(i):
        eqs.append(sp.expand(sp.diff(dm2, y, k).subs(y, -1)))
    for k in range(n1 - i):
        eqs.append(sp.expand(sp.diff(dm3, y, k).subs(y, -1)))
    if n2 == 1:
        eqs.append(sp.expand((dm2 if j == 1 else dm3).subs(y, R)))
    n_spare_case = (CAPS["dm2"] + CAPS["dm3"] + 2) - len(eqs)
    return eqs, n_spare_case


def build_case(base_bs: dict, *, label, i, j, n1, n2) -> dict:
    """Extend a CACHED dm4-eliminated state build (r9_symbolic_sweep
    build_reduced output) with the case's linear valuation equations.
    No new symbolic expansion happens here."""
    extra, n_spare = case_equations(i, j, n1, n2)
    eqs = list(base_bs["equations"]) + [e for e in extra if e != 0]
    sat = list(base_bs["sat_factors"])
    rv = bsw._ring_vars(eqs, sat)
    return {"label": label, "equations": eqs, "sat_factors": sat,
            "ring_vars": rv, "n_equations": len(eqs), "n_unknowns": len(rv),
            "sizes": dict(base_bs["sizes"], valsplit=len(extra)),
            "n_spare": n_spare, "case": [i, j]}


# --------------------------------------------------------------------------
#  recording + per-case attack (mirrors rss.attack; records to OWN file)
# --------------------------------------------------------------------------
def _record(rec):
    prior = json.load(open(OUT)) if OUT.exists() else {"states": {},
                                                       "cases": []}
    prior["cases"] = [x for x in prior["cases"]
                      if x["label"] != rec["label"]] + [rec]
    json.dump(prior, open(OUT, "w"), indent=1, default=str)


def _record_state(label, summary):
    prior = json.load(open(OUT)) if OUT.exists() else {"states": {},
                                                       "cases": []}
    prior.setdefault("states", {})[label] = summary
    json.dump(prior, open(OUT, "w"), indent=1, default=str)


def attack_case(bs, *, marked_root, exact_timeout=300.0, triage_timeout=45.0):
    print(f"\n{bs['label']}: {bs['n_equations']} eqs, {bs['n_unknowns']} vars "
          f"(spare={bs['n_spare']}, case={bs['case']})", flush=True)
    t0 = time.monotonic()
    if marked_root:
        tri = bsw.triage_bridge_numroot(bs, timeout=triage_timeout)
    else:
        tri = bsw.triage_bridge(bs, timeout=triage_timeout)
    print(f"  mod-p: {tri['prediction']}", flush=True)
    exact = None
    verdicts = [p.get("verdict") for p in tri["primes"]]
    if "PROPER" not in verdicts and any(v == "UNIT" for v in verdicts):
        if marked_root:
            exact = bsw.exact_bridge_minpoly(bs, timeout=exact_timeout)
            if exact.get("verdict") == "PARSE_FAIL":
                exact2 = bsw.exact_bridge(bs, timeout=exact_timeout)
                if exact2.get("verdict") in ("UNIT", "PROPER"):
                    exact = exact2
        else:
            exact = bsw.exact_bridge(bs, timeout=exact_timeout)
    if exact and exact.get("verdict") == "UNIT":
        verdict = "CASE-KILLED (PENDING AUDIT)"
    elif (exact and exact.get("verdict") == "PROPER") or "PROPER" in verdicts:
        verdict = "CASE-PROPER (INCONCLUSIVE)"
    else:
        verdict = "COST"
    rec = {"label": bs["label"], "case": bs["case"],
           "n_equations": bs["n_equations"], "n_unknowns": bs["n_unknowns"],
           "modp": tri, "exact": exact, "verdict": verdict,
           "wall": round(time.monotonic() - t0, 1)}
    if verdict.startswith("CASE-KILLED"):
        rec["kill_system"] = rss._kill_payload(bs)
    _record(rec)
    print(f"  ==> {bs['label']}: {verdict} ({rec['wall']}s)", flush=True)
    return rec


def state_verdict(case_recs, n_cases, attempted_all):
    kills = sum(1 for r in case_recs
                if r["verdict"].startswith("CASE-KILLED"))
    if kills == n_cases:
        return f"KILLED-VALSPLIT (PENDING AUDIT; all {n_cases} cases UNIT)"
    if kills:
        return (f"DENTED ({kills}/{n_cases} cases killed"
                f"{'' if attempted_all else '; truncated'})")
    return ("WALL SURVIVES (0 case kills"
            + ("" if attempted_all else "; truncated") + ")")


# --------------------------------------------------------------------------
#  drivers
# --------------------------------------------------------------------------
def run_r9_column(z: int, deadline, **kw):
    n1, n2 = 9, 1
    base = rss.build_r9_reduced(z)          # cached dm4-eliminated build
    recs, attempted_all = [], True
    for i, j in cases_for(n1, n2):
        if time.monotonic() > deadline:
            attempted_all = False
            print(f"  R9_z{z}: budget exhausted at case ({i},{j})",
                  flush=True)
            break
        label = f"R9_z{z}_i{i}_j{j}"
        try:
            bs = build_case(base, label=label, i=i, j=j, n1=n1, n2=n2)
        except Exception as ex:
            recs.append({"label": label, "verdict": f"BUILD_ERROR: {ex}"[:200],
                         "case": [i, j]})
            _record(recs[-1])
            continue
        recs.append(attack_case(bs, marked_root=True, **kw))
    sv = state_verdict(recs, len(cases_for(n1, n2)), attempted_all)
    _record_state(f"R9_z{z}", {"verdict": sv, "n_cases": len(recs),
                               "cases_killed": sum(
                                   1 for r in recs
                                   if r["verdict"].startswith("CASE-KILLED"))})
    print(f"\n== R9_z{z}: {sv}", flush=True)
    return recs


def run_batch_state(s, deadline, **kw):
    a = int(s["a_t"])
    n1, n2 = a, 0
    label0 = rss.batch_label(s)
    base = rss.build_batch_reduced(s)       # cached dm4-eliminated build
    recs, attempted_all = [], True
    for i, j in cases_for(n1, n2):
        if time.monotonic() > deadline:
            attempted_all = False
            print(f"  {label0}: budget exhausted at case ({i},{j})",
                  flush=True)
            break
        label = f"{label0}_i{i}"
        try:
            bs = build_case(base, label=label, i=i, j=j, n1=n1, n2=n2)
        except Exception as ex:
            recs.append({"label": label, "verdict": f"BUILD_ERROR: {ex}"[:200],
                         "case": [i, j]})
            _record(recs[-1])
            continue
        recs.append(attack_case(bs, marked_root=False, **kw))
    sv = state_verdict(recs, len(cases_for(n1, n2)), attempted_all)
    _record_state(label0, {"verdict": sv, "n_cases": len(recs),
                           "cases_killed": sum(
                               1 for r in recs
                               if r["verdict"].startswith("CASE-KILLED"))})
    print(f"\n== {label0}: {sv}", flush=True)
    return recs


def census():
    if not OUT.exists():
        print("no valsplit record yet")
        return
    d = json.load(open(OUT))
    from collections import Counter
    c = Counter(x["verdict"].split(" ")[0] for x in d["cases"])
    print("case verdicts:", dict(c))
    for lbl, s in d.get("states", {}).items():
        print(f"  {lbl}: {s['verdict']}")


def main(budget_min=150.0):
    # exhaustiveness first -- a checked claim, not an assumption
    n = check_exhaustive(9, 1, CAPS["dm2"], CAPS["dm3"])
    print(f"[OK] R9 split exhaustive: 20 cases cover {n} valuation profiles",
          flush=True)
    for a in (7, 8, 9):
        n = check_exhaustive(a, 0, CAPS["dm2"], CAPS["dm3"])
        print(f"[OK] batch a={a} split exhaustive: {a + 1} cases cover {n} "
              f"profiles", flush=True)
    deadline = time.monotonic() + budget_min * 60
    run_r9_column(1, deadline)
    for z in (2, 3):
        if time.monotonic() > deadline:
            print(f"z={z}: not attempted (budget)", flush=True)
            _record_state(f"R9_z{z}", {"verdict": "NOT ATTEMPTED (budget)"})
            continue
        run_r9_column(z, deadline)
    states = rss.dege10_t2_states()[:8]
    for s in states:
        if time.monotonic() > deadline:
            _record_state(rss.batch_label(s),
                          {"verdict": "NOT ATTEMPTED (budget)"})
            continue
        run_batch_state(s, deadline)
    census()


if __name__ == "__main__":
    args = sys.argv[1:]
    if args and args[0] == "census":
        census()
    elif args and args[0] == "exhaustive":
        n = check_exhaustive(9, 1, CAPS["dm2"], CAPS["dm3"])
        print(f"[OK] R9: 20 cases cover {n} profiles")
        for a in (7, 8, 9):
            n = check_exhaustive(a, 0, CAPS["dm2"], CAPS["dm3"])
            print(f"[OK] batch a={a}: {a + 1} cases cover {n} profiles")
        print("EXHAUSTIVENESS CHECKS PASSED")
    elif args and args[0] == "case":        # single-case debug entry
        z, i, j = int(args[1]), int(args[2]), int(args[3])
        base = rss.build_r9_reduced(z)
        bs = build_case(base, label=f"R9_z{z}_i{i}_j{j}", i=i, j=j,
                        n1=9, n2=1)
        attack_case(bs, marked_root=True)
    else:
        main(float(args[0]) if args else 150.0)
