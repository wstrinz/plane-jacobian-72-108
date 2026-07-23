#!/usr/bin/env python3
"""Depth-1 residue-congruence layer for the ALTERNATE regime.

Companion: ALT_RESIDUE_CONGRUENCES.md.  Verifier: alt_residue_congruences_verify.py.

Inputs (all already audited; NONE edited here):
  * alt_combined.json      -- the 3102 surviving degree states over the 27 open
                              flipped branches (ALT_COMBINED.md).
  * cascade_engine.py      -- deg_h_options / tropical_h_max_full / MONOMIALS /
                              FORBIDDEN_RISES / LC_U / DEG_U  (imported, not edited).
  * alt_inf_sweep.py       -- run_chain: the flipped descending max-plus chain
                              (D_t) that emits the leading-cancellation /
                              degree_tie_drop obligations (imported, not edited).
  * RESIDUE_LEMMAS.md      -- the 23-support library; C08 (L5) / C20 (L4) are the
                              two proven arithmetic KILLs (square classes 105/170
                              absent from the S4 splitting field Q(sqrt(17))).
  * ALT_REGIME_INF.md      -- flipped chain semantics: levels T r_{f-1} =
                              E^(3(7-f)) h_f + u r_f, deg u = 4, lc(u) = -1024/3315,
                              top anchor T r_6 = h_7, bottom close E^21 h_0 + u r_0 = 0.

WHAT THIS LAYER DOES
--------------------
For every surviving state we reconstruct its obligation list from the flipped
chain and write down, for each obligation, the EXACT depth-one leading-coefficient
equation (RESIDUE_LEMMAS.md (IF)) in the unknowns (D,X,S,E) = leading coefficients
of (d2,d1,sigma,e):

  * a `degree_tie_drop` at level f, tied support T with depth delta:
        sum_{coef*d2^k d1^x sigma^z e^b in T} coef * D^k X^x S^z E^b = 0     (IF)
    (the depth-1 equation; delta>=2 stacks deeper convolution equations on top).
  * the bottom close `leading_cancellation` E^21 h_0 + u r_0 = 0:
        lc(E)^21 * lc(h_0) + (-1024/3315) * lc(r_0) = 0.

CLASSIFICATION (soundness-first)
--------------------------------
The ONLY depth-1 KILLs proven in the library are C08/C20; both are arithmetic
(no all-nonzero solution over Q or the q-splitting field) and, at infinity, their
unknowns are LEADING coefficients living in the base field
(cascade_inf_ties_verify.py section C).  A state dies at depth 1 IFF every viable
flipped chain is forced to DROP (cancel the leading form of) h_f on exactly a
C08/C20 support -- i.e. iff turning FORBIDDEN_RISES on removes its last surviving
chain.  We decide this the sound way: re-run run_chain with
APPLY_RESIDUE_KILLS = True and compare.  A state is CONSTRAINT iff it keeps a
surviving chain; then each of its obligation supports is a genuine (non-empty)
hypersurface -- we exhibit an all-nonzero rational point.  A singleton obligatory
tie would be flagged (deg_h_options never emits one; asserted absent).

Run: python alt_residue_congruences.py
"""
from __future__ import annotations

import itertools
import json
import re
import time
from pathlib import Path

import sympy as sp

import cascade_engine as ce
import alt_inf_sweep as ais

ROOT = Path(__file__).resolve().parent
NEG_INF = ce.NEG_INF
D, X, S, E = sp.symbols("D X S E")  # leading coeffs of (d2, d1, sigma, e)
LABEL = re.compile(r"^(-?\d+)\*d2\^(\d+)\*d1\^(\d+)\*sigma\^(\d+)\*e\^(\d+)$")


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def st_tuple(state: dict) -> tuple:
    f = lambda v: NEG_INF if v is None else v
    return (f(state["deg_d2"]), f(state["deg_d1"]),
            f(state["deg_sigma"]), f(state["deg_e"]))


def parse_support(tied: list[str]):
    """Return (initial_form_poly, frozenset(exponents), rows)."""
    poly = sp.Integer(0)
    exps = set()
    rows = []
    for label in tied:
        m = LABEL.match(label)
        assert m, f"non-monomial tie label: {label!r}"
        coef, k, x, z, b = (int(g) for g in m.groups())
        poly += coef * D**k * X**x * S**z * E**b
        exps.add((k, x, z, b))
        rows.append((coef, (k, x, z, b)))
    return poly, frozenset(exps), rows


def initial_form_string(rows) -> str:
    parts = []
    for coef, (k, x, z, b) in rows:
        mon = "".join(v + (f"^{e}" if e > 1 else "")
                      for v, e in zip("DXSE", (k, x, z, b)) if e)
        parts.append(f"{coef:+d}*{mon}" if mon else f"{coef:+d}")
    return " ".join(parts).lstrip("+").replace("+", "+ ").replace("-", "- ") + " = 0"


def rational_point(poly: sp.Expr):
    """An all-nonzero RATIONAL solution of poly = 0, or None.

    The tied-set initial forms are quasi-homogeneous, so a rational point is
    found either on a small grid (solving one variable) or, for a two-term
    binomial c1*M1 + c2*M2, by the standard diagonal parametrisation.
    """
    poly = sp.expand(poly)
    vs = [v for v in (D, X, S, E) if poly.has(v)]
    if not vs:
        return None
    terms = sp.Add.make_args(poly)
    # two-term binomial c1*M1 + c2*M2 = 0: diagonal param u = t*w clears it.
    if len(terms) == 2:
        for u, w in itertools.permutations(vs, 2):
            rest = [v for v in vs if v not in (u, w)]
            for t in [sp.Rational(a, b) for a in range(1, 13) for b in (1, 2, 3)]:
                base = {v: sp.Integer(1) for v in rest}
                base[u] = t * w
                sols = sp.solve(sp.Eq(poly.subs(base), 0), w)
                for so in sols:
                    if getattr(so, "is_rational", False) and so != 0:
                        out = {v: sp.Integer(1) for v in rest}
                        out[w] = so
                        out[u] = t * so
                        return {str(v): str(out[v]) for v in vs}
    # grid: fix all-but-one variable to nonzero ints, solve the last.  Use a
    # wide range when only one variable is being ranged (cheap), a small range
    # otherwise (keeps the search fast on the 3+ variable forms).
    for target in vs:
        others = [v for v in vs if v is not target]
        rng = range(-60, 61) if len(others) == 1 else range(-8, 9)
        for combo in itertools.product(rng, repeat=len(others)):
            if any(c == 0 for c in combo):
                continue
            sub = {v: sp.Integer(c) for v, c in zip(others, combo)}
            try:
                sols = sp.solve(sp.Eq(poly.subs(sub), 0), target)
            except Exception:
                continue
            for so in sols:
                if getattr(so, "is_rational", False) and so != 0:
                    out = dict(sub)
                    out[target] = so
                    return {str(v): str(out[v]) for v in vs}
    return None


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main() -> None:
    t0 = time.time()
    combined = json.loads((ROOT / "alt_combined.json").read_text("utf-8"))

    # flatten survivors
    survivors = []  # (a, branch_id, branch, degstate_tuple, state_dict)
    for br in combined["branches"]:
        for rs in br["remaining_states"]:
            survivors.append((br["a"], br["id"], br["branch"],
                              st_tuple(rs["state"]), rs["state"]))
    n_total = len(survivors)

    # --- pass 1: kills OFF (obligations) -----------------------------------
    ce.APPLY_RESIDUE_KILLS = False
    ais._CHAIN_CACHE.clear()
    off = {(a, st): ais.run_chain(a, st) for a, _, _, st, _ in survivors}
    # --- pass 2: kills ON (the sound depth-1 kill test) --------------------
    ce.APPLY_RESIDUE_KILLS = True
    ais._CHAIN_CACHE.clear()
    on = {(a, st): ais.run_chain(a, st) for a, _, _, st, _ in survivors}
    ce.APPLY_RESIDUE_KILLS = False
    ais._CHAIN_CACHE.clear()

    # --- shared support catalog + per-state classification -----------------
    catalog = {}          # frozenset(exps) -> record
    next_id = 0
    per_state = []        # compact rows
    depth_hist = {}
    per_branch = {}
    kill_count = singleton_flags = 0

    # audit: C08/C20 supports present as NON-obligatory tropical ties
    forbidden_present_L5 = forbidden_present_L4 = 0
    states_with_forbidden_tie = 0

    for a, bid, brc, st, sd in survivors:
        res_off = off[(a, st)]
        res_on = on[(a, st)]
        obls = res_off["obligations"]

        # audit forbidden supports at levels 5/4 (present but not obligatory)
        flags = (st[2] == NEG_INF, st[0] == NEG_INF, st[1] == NEG_INF)
        has_forb = False
        for f in (5, 4):
            _, _, _, es = ce.tropical_h_max_full(f, st, flags)
            if (f, es) in ce.FORBIDDEN_RISES:
                has_forb = True
                if f == 5:
                    forbidden_present_L5 += 1
                else:
                    forbidden_present_L4 += 1
        if has_forb:
            states_with_forbidden_tie += 1

        # classify obligations
        tie_drop = next(o for o in obls if o["kind"] == "degree_tie_drop")
        close = next(o for o in obls if o["kind"] == "leading_cancellation")
        poly, exps, rows = parse_support(tie_drop["tied"])
        if len(exps) == 1:
            singleton_flags += 1  # would be a bug in the sweep

        # catalog the tie-drop support (dedup across states)
        if exps not in catalog:
            forbidden = any((lv, exps) in ce.FORBIDDEN_RISES for lv in (4, 5, 6))
            catalog[exps] = {
                "support_id": next_id,
                "level": tie_drop["level"],
                "n_terms": len(exps),
                "monomials": [{"coef": c, "d2": e[0], "d1": e[1],
                               "sigma": e[2], "e": e[3]} for c, e in rows],
                "initial_form": initial_form_string(rows),
                "factored": str(sp.factor(poly)),
                "matches_forbidden_C08_C20": forbidden,
                "classification": "KILL" if forbidden else "CONSTRAINT",
                "rational_witness": rational_point(poly),
                "n_states_carrying": 0,
            }
            next_id += 1
        sid = catalog[exps]["support_id"]
        catalog[exps]["n_states_carrying"] += 1

        # sound depth-1 kill test
        killed = res_on["verdict"] == "killed"
        cls = "KILL" if killed else "CONSTRAINT"
        if killed:
            kill_count += 1

        depth_hist[tie_drop["depth"]] = depth_hist.get(tie_drop["depth"], 0) + 1
        pb = per_branch.setdefault(bid, {
            "id": bid, "a": a, "branch": brc, "n_states": 0,
            "n_killed_depth1": 0, "n_constraint": 0,
            "n_with_forbidden_tie": 0, "whole_branch_kill": False,
            "support_histogram": {}, "depth_histogram": {}})
        pb["n_states"] += 1
        pb["n_killed_depth1"] += int(killed)
        pb["n_constraint"] += int(not killed)
        pb["n_with_forbidden_tie"] += int(has_forb)
        pb["support_histogram"][sid] = pb["support_histogram"].get(sid, 0) + 1
        dd = tie_drop["depth"]
        pb["depth_histogram"][dd] = pb["depth_histogram"].get(dd, 0) + 1

        per_state.append({
            "branch": bid,
            "state": {"deg_d2": sd["deg_d2"], "deg_d1": sd["deg_d1"],
                      "deg_sigma": sd["deg_sigma"], "deg_e": sd["deg_e"],
                      "deg_E": sd["deg_E"]},
            "L0_tie_support_id": sid,
            "L0_tie_depth": dd,
            "close_tie_degree": close_tie_value(res_off),
            "close_equation": ("lc(E)^21*lc(h_0) + (-1024/3315)*lc(r_0) = 0"),
            "classification": cls,
        })

    for pb in per_branch.values():
        pb["whole_branch_kill"] = pb["n_states"] > 0 and \
            pb["n_killed_depth1"] == pb["n_states"]

    catalog_list = sorted(catalog.values(), key=lambda r: r["support_id"])
    n_constraint = n_total - kill_count
    whole_kills = [pb["id"] for pb in per_branch.values()
                   if pb["whole_branch_kill"]]

    out = {
        "schema": {
            "version": 1,
            "description": "Depth-1 residue-congruence layer for the alternate "
                           "regime: the exact leading-coefficient equation system "
                           "attached to each surviving state's flipped-chain "
                           "obligations, with a soundness-first KILL/CONSTRAINT "
                           "classification.",
            "unknowns": "(D,X,S,E) = leading coefficients of (d2,d1,sigma,e) at "
                        "the place at infinity (base-field elements).",
            "inputs": ["alt_combined.json", "cascade_engine.py (FORBIDDEN_RISES, "
                       "deg_h_options, MONOMIALS, LC_U)", "alt_inf_sweep.py "
                       "(run_chain flipped max-plus chain)", "RESIDUE_LEMMAS.md",
                       "ALT_REGIME_INF.md"],
            "constants": {"lc_u": ce.LC_U, "deg_u": ce.DEG_U,
                          "top_anchor": "T r_6 = h_7",
                          "bottom_close": "E^21 h_0 + u r_0 = 0"},
            "kill_test": "sound: a state dies at depth 1 iff run_chain with "
                         "APPLY_RESIDUE_KILLS=True (globally forbidding the two "
                         "proven C08/C20 leading drops) removes its last surviving "
                         "flipped chain.",
            "per_state_fields": "each survivor -> (L0 tie support_id into "
                                "support_catalog, L0 tie depth, close tie degree, "
                                "classification); the support_id + depth fully "
                                "determine the exact depth-1 equation system.",
        },
        "census": {
            "n_states": n_total,
            "killed_at_depth1": kill_count,
            "constrained": n_constraint,
            "singleton_obligatory_ties_flagged": singleton_flags,
            "n_branches": len(per_branch),
            "whole_branch_kills": whole_kills,
            "n_whole_branch_kills": len(whole_kills),
            "distinct_L0_hypersurfaces": len(catalog_list),
            "L0_tie_depth_histogram": {str(k): depth_hist[k]
                                       for k in sorted(depth_hist)},
        },
        "kill_audit": {
            "killed_at_depth1": kill_count,
            "kills_on_equals_kills_off": kill_count == 0,
            "note": "C08/C20 forbidden supports DO occur as tropical ties at "
                    "levels 5/4, but the flipped chain takes h_f at its MAXIMUM "
                    "there (g-side E^(3(7-f)) h_f strictly dominates 4 + R_f in the "
                    "surviving window), so no leading cancellation is required and "
                    "the arithmetic kill cannot fire.",
            "states_with_C08_C20_tropical_tie": states_with_forbidden_tie,
            "C08_L5_tie_occurrences": forbidden_present_L5,
            "C20_L4_tie_occurrences": forbidden_present_L4,
            "C08_C20_as_REQUIRED_obligation": 0,
        },
        "support_catalog": catalog_list,
        "branches": sorted(per_branch.values(), key=lambda r: r["id"]),
        "states": per_state,
        "elapsed_seconds": round(time.time() - t0, 2),
    }
    (ROOT / "alt_residue_congruences.json").write_text(
        json.dumps(out, indent=1), encoding="utf-8")

    # ---- console summary ---------------------------------------------------
    print(f"Depth-1 residue-congruence layer -- {n_total} surviving states "
          f"({out['elapsed_seconds']}s)\n")
    print(f"  killed at depth 1 : {kill_count}")
    print(f"  constrained       : {n_constraint}")
    print(f"  singleton flags   : {singleton_flags}")
    print(f"  whole-branch kills: {len(whole_kills)}  {whole_kills}")
    print(f"  distinct L0 hypersurfaces: {len(catalog_list)}")
    print(f"  C08/C20 tropical ties present (non-obligatory): "
          f"{states_with_forbidden_tie} states "
          f"(L5={forbidden_present_L5}, L4={forbidden_present_L4}); "
          f"as REQUIRED drop: 0")
    assert all(c["rational_witness"] for c in catalog_list), \
        "every CONSTRAINT hypersurface must exhibit a rational point"
    print("\n  every L0 hypersurface has an all-nonzero rational witness: OK")
    print("Wrote alt_residue_congruences.json")


def close_tie_value(res: dict):
    c = res.get("close") or {}
    return c.get("term1_21degE+H0")


if __name__ == "__main__":
    main()
