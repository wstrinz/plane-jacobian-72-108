#!/usr/bin/env python3
"""divisor_filter.py -- the e | Phi divisor filter as a FIRST-CLASS PIPELINE STAGE.

WHAT THIS IS
------------
`divisor_syzygy.py` proves the universal K-syzygy of the canonical G-system

    2*(G5 + d2*G3 + d1*G2 + d0*G1)  ==  2*Phi - e*(d2*e^2 + 3*e*S + 3*R^2)

(residual exactly 0), hence on every genuine lift  2*Phi = e*(...)  and so

    e | Phi ,     Phi = -(1/6630)*(y+1)^30*q,   q SQUAREFREE, deg Phi = 34.

THIS file turns that lemma into a stage the authoritative compiler can run over a
Phase-D state universe.  It exports one predicate pair -- `cell_verdict` /
`state_verdict` -- and a `filter_universe` driver; `phase_d_states.py --divisor-filter`
calls them, so the filtered universes are produced BY the universe generator, not
by a hand count.

THE THREE CONSEQUENCES, and exactly where each bites
----------------------------------------------------
Write  e = gamma * t^a * prod_i (y - r_i)^{b_i} * (off-support factor),
t = y+1, r_1..r_4 the four (simple) roots of q.  The Phase-D universe indexes a
CELL by (a_t, b, branch) and a STATE inside it by (deg_d2, deg_d1, deg_sigma, deg_e).

  (D1) rad(e) | (y+1)*q  -- e has NO root off {-1, r_1..r_4}.
       => the off-support factor is a unit => deg e = a_t + sum(b_i) EXACTLY.
       This is a STATE-level filter.  It is NOT automatic from the universe:
       the universe enumerates deg_e over the whole window range
       [a_t + sum(b), e_cap], i.e. it admits DEFECT-d `e` carrying d free extra
       roots (PHASE_F2_SUB2.md).  The divisibility is what forces defect 0.

  (D2) b_i in {0,1}  -- each SIMPLE root of q divides e to order <= 1.
       This is a CELL-level filter: any cell whose b-vector has an entry >= 2
       dies whole, top stratum included.

  (D3) The degree count.  deg(RHS) = deg e + max(deg d2 + 2 deg e,
       deg e + deg S, 2 deg R) must be able to reach deg Phi = 34.  Since the
       bound is monotone in deg e this yields a window-dependent LOWER bound
       E_min(window); combined with (D1) it is a CELL-level filter
       a_t + sum(b_i) >= E_min.

       *** E_min IS DERIVED PER WINDOW HERE, NOT ASSUMED. ***
       sub2: caps (d2,R,S,e) = (4,12,14,10) give E_min = 10 = e_cap, so
             deg e = 10 EXACTLY and a_t + sum(b_i) = 10.
       sub1: caps (d2,R,S,e) = (6,18,21,15) give 2*deg R = 36 > 34 already at
             deg e = 0, so E_min = 0 and (D3) is VACUOUS in sub1.  The sub2
             degree forcing DOES NOT TRANSFER; check C5 asserts this so the
             transfer can never be assumed by accident.

CAP PROVENANCE (nothing typed in by hand)
-----------------------------------------
    deg d2, deg e   <- cascade_engine.CONFIGS[window].aux_caps / .e_cap
    deg R = deg dm2 <- full_system_bridge.STRIP_DEGCAP[window]["dm2"]
    deg S = deg dm3 <- full_system_bridge.STRIP_DEGCAP[window]["dm3"]
    deg Phi = 34    <- bigrade_annotator._phi_stripped(), degree recomputed here
    the syzygy       <- divisor_syzygy.syzygy_residual(), residual re-asserted here

C08/C20 INDEPENDENCE
--------------------
This filter reads only (a_t, b, deg_e).  It is ORTHOGONAL to the residue-kill
toggle (`cascade_engine.APPLY_RESIDUE_KILLS`, the C08/C20 forbidden-rise lemmas),
which changes WHICH cells/states exist, not which of them the divisor lemma
kills.  The same predicate is therefore applied to both universes; see
`frontier_rebuild.py` for the 2x2 census.

USAGE
-----
    python divisor_filter.py --quiet            # self-check, exit 0 iff all pass
    python divisor_filter.py                    # full report + census
    python divisor_filter.py --universe phase_d_states_sub2.json --out FILTERED.json
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from typing import Sequence

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

DEG_PHI_EXPECTED = 34


# ---------------------------------------------------------------------------
# caps, read from the authoritative sources
# ---------------------------------------------------------------------------
def window_caps(window: str) -> dict[str, int]:
    """(d2, R, S, e) stripped-degree caps for `window`, from the cap tables."""
    import cascade_engine as ce
    import full_system_bridge as fsb

    cfg = ce.CONFIGS[window]
    strip = fsb.STRIP_DEGCAP[window]
    return {
        "d2": cfg.aux_caps[2],
        "R": strip["dm2"],
        "S": strip["dm3"],
        "e": cfg.e_cap,
    }


def deg_phi() -> int:
    """deg Phi, recomputed from the concrete stripped Phi (window-independent)."""
    import sympy as sp
    import bigrade_annotator as ba

    return sp.Poly(ba._phi_stripped(), ba.y).degree()


def rhs_degree_table(caps: dict[str, int]) -> list[tuple[int, int]]:
    """[(E, max possible deg of e*(d2 e^2 + 3 e S + 3 R^2))] for E = 0..e_cap+1."""
    return [
        (E, E + max(caps["d2"] + 2 * E, E + caps["S"], 2 * caps["R"]))
        for E in range(0, caps["e"] + 2)
    ]


def forced_deg_e_min(caps: dict[str, int], deg_phi_value: int = DEG_PHI_EXPECTED) -> int | None:
    """Least deg e whose RHS degree bound can still reach deg Phi.

    The bound is monotone non-decreasing in E, so the least feasible E is a true
    lower bound on deg e.  Returns None if NO E <= e_cap works (window empty).
    """
    for E, rhs in rhs_degree_table(caps):
        if E > caps["e"]:
            break
        if rhs >= deg_phi_value:
            return E
    return None


# ---------------------------------------------------------------------------
# the predicates
# ---------------------------------------------------------------------------
DEATH_B_GE_2 = "b_i>=2 at a simple q-root (e | Phi forbids it)"
DEATH_DEG = "a+sum(b) below the forced deg e"
DEATH_DEFECT = "deg_e != a+sum(b) (off-support root; e | Phi forbids it)"


class DivisorFilter:
    """The e | Phi filter, specialised to one window's caps."""

    def __init__(self, window: str, deg_phi_value: int | None = None):
        self.window = window
        self.caps = window_caps(window)
        self.deg_phi = DEG_PHI_EXPECTED if deg_phi_value is None else deg_phi_value
        self.e_min = forced_deg_e_min(self.caps, self.deg_phi)
        # A window in which the degree count alone is unsatisfiable would be a
        # far stronger statement than anything claimed; refuse to proceed
        # silently if that ever happens.
        if self.e_min is None:
            raise SystemExit(
                "FATAL: window %s admits no deg e <= %d reaching deg Phi = %d"
                % (window, self.caps["e"], self.deg_phi)
            )
        self.degree_forcing_active = self.e_min > 0

    # -- cell level ---------------------------------------------------------
    def cell_verdict(self, a_t: int, b: Sequence[int]) -> tuple[bool, str]:
        """(alive, reason). Applies (D2) then (D3); the two are independent tests."""
        if any(x >= 2 for x in b):
            return False, DEATH_B_GE_2
        support = a_t + sum(b)
        if support < self.e_min:
            return False, "%s: a+sum(b)=%d < E_min=%d" % (DEATH_DEG, support, self.e_min)
        return True, "survives"

    # -- state level --------------------------------------------------------
    def state_verdict(self, a_t: int, b: Sequence[int], deg_e: int) -> tuple[bool, str]:
        alive, why = self.cell_verdict(a_t, b)
        if not alive:
            return False, why
        support = a_t + sum(b)
        if deg_e != support:
            return False, "%s: deg_e=%d != a+sum(b)=%d" % (DEATH_DEFECT, deg_e, support)
        return True, "survives"


# ---------------------------------------------------------------------------
# the stage: filter a phase-D universe
# ---------------------------------------------------------------------------
def filter_universe(universe: dict, window: str | None = None) -> dict:
    """Apply the filter to a phase_d_states_*.json payload; returns a new payload.

    Surviving flag cases keep only surviving states; flag cases with no surviving
    state are dropped, which drops their cell when every flag case in it dies.
    A `divisor_filter` block records the census and the death attribution.
    """
    win = window or universe["window"]
    filt = DivisorFilter(win)

    kept_cases = []
    cell_alive: dict[tuple, bool] = {}
    census = {
        "cells_before": 0, "cells_after": 0,
        "flagcases_before": len(universe["cases"]), "flagcases_after": 0,
        "states_before": 0, "states_after": 0,
        "cells_killed_b_ge_2": 0, "cells_killed_degree": 0,
        "cells_killed_defect": 0,
        "flagcases_killed_b_ge_2": 0, "flagcases_killed_degree": 0,
        "flagcases_killed_defect": 0,
        "states_killed_b_ge_2": 0, "states_killed_degree": 0,
        "states_killed_defect": 0,
    }
    cells_seen: dict[tuple, dict] = {}

    for case in universe["cases"]:
        cell = (case["a_t"], tuple(case["b"]), case["branch"])
        rec = cells_seen.setdefault(cell, {"states": 0, "kept": 0})
        alive_cell, why_cell = filt.cell_verdict(case["a_t"], case["b"])
        cell_alive[cell] = alive_cell
        kept_states = []
        for st in case["states"]:
            census["states_before"] += 1
            rec["states"] += 1
            alive, why = filt.state_verdict(case["a_t"], case["b"], st["deg_e"])
            if alive:
                kept_states.append(st)
                rec["kept"] += 1
            elif why.startswith(DEATH_B_GE_2):
                census["states_killed_b_ge_2"] += 1
            elif why.startswith(DEATH_DEG):
                census["states_killed_degree"] += 1
            else:
                census["states_killed_defect"] += 1
        if kept_states:
            nc = dict(case)
            nc["states"] = kept_states
            nc["state_count"] = len(kept_states)
            kept_cases.append(nc)
            census["states_after"] += len(kept_states)
            census["flagcases_after"] += 1
        elif not alive_cell:
            if why_cell.startswith(DEATH_B_GE_2):
                census["flagcases_killed_b_ge_2"] += 1
            else:
                census["flagcases_killed_degree"] += 1
        else:
            # cell alive, but every state in this flag case carried defect > 0
            census["flagcases_killed_defect"] += 1

    census["cells_before"] = len(cells_seen)
    census["cells_after"] = sum(1 for c, r in cells_seen.items() if r["kept"])
    for cell, rec in cells_seen.items():
        if rec["kept"]:
            continue
        if cell_alive[cell]:
            # cell-level tests pass, but every state in it carried defect > 0
            census["cells_killed_defect"] += 1
            continue
        alive, why = filt.cell_verdict(cell[0], cell[1])
        if why.startswith(DEATH_B_GE_2):
            census["cells_killed_b_ge_2"] += 1
        else:
            census["cells_killed_degree"] += 1

    out = dict(universe)
    out["cases"] = kept_cases
    out["case_count"] = len(kept_cases)
    out["state_total"] = census["states_after"]
    out["description"] = universe.get("description", "") + " + divisor filter (e | Phi)"
    out["divisor_filter"] = {
        "applied": True,
        "window": win,
        "caps": filt.caps,
        "deg_phi": filt.deg_phi,
        "forced_deg_e_min": filt.e_min,
        "degree_forcing_active": filt.degree_forcing_active,
        "source_lemma": "divisor_syzygy.py / DIVISOR_SYZYGY.md (2*Phi = e*(d2 e^2 + 3 e S + 3 R^2))",
        "census": census,
        "surviving_cells": sorted(
            "a%d_b%s_%s" % (c[0], "".join(map(str, c[1])), c[2])
            for c, r in cells_seen.items() if r["kept"]
        ),
    }
    return out


# ---------------------------------------------------------------------------
# self-check
# ---------------------------------------------------------------------------
def _load(fn):
    with open(os.path.join(HERE, fn), encoding="utf-8") as fh:
        return json.load(fh)


def run(verbose: bool = True) -> int:
    out, npass, nfail = [], 0, 0

    def check(name, cond, detail=""):
        nonlocal npass, nfail
        if cond:
            npass += 1
            out.append("  [OK] %s%s" % (name, (" -- " + detail) if detail and verbose else ""))
        else:
            nfail += 1
            out.append("  [FAIL] %s%s" % (name, (" -- " + detail) if detail else ""))

    # C1 -- the premise: the K-syzygy residual is exactly 0.
    import divisor_syzygy as ds
    resid, _G5 = ds.syzygy_residual()
    check("C1 K-syzygy residual is exactly 0 (the premise of e | Phi)", resid == 0,
          "residual=%s" % resid)

    # C2 -- deg Phi recomputed from the concrete stripped Phi.
    dphi = deg_phi()
    check("C2 deg Phi recomputed = %d" % DEG_PHI_EXPECTED, dphi == DEG_PHI_EXPECTED,
          "got %d" % dphi)

    # C3 -- sub2 caps read from the authoritative tables match DIVISOR_SYZYGY.md.
    c2 = window_caps("sub2")
    check("C3 sub2 caps (d2,R,S,e) = (4,12,14,10) from cascade_engine + full_system_bridge",
          c2 == {"d2": 4, "R": 12, "S": 14, "e": 10}, str(c2))

    # C4 -- sub2 degree forcing: E_min = 10 = e_cap, i.e. deg e = 10 EXACTLY.
    f2 = DivisorFilter("sub2")
    check("C4 sub2 forced deg e: E_min = 10 = e_cap (deg e = 10 exactly)",
          f2.e_min == 10 and f2.caps["e"] == 10,
          "E_min=%s  table=%s" % (f2.e_min, rhs_degree_table(c2)[6:12]))

    # C5 -- sub1 does NOT inherit it.  2*deg R = 36 > 34 already at deg e = 0.
    c1 = window_caps("sub1")
    f1 = DivisorFilter("sub1")
    check("C5 sub1 caps (d2,R,S,e) = (6,18,21,15)",
          c1 == {"d2": 6, "R": 18, "S": 21, "e": 15}, str(c1))
    check("C5b sub1 degree forcing is VACUOUS (E_min = 0): 2*deg R = %d > deg Phi = %d, "
          "so the sub2 'deg e = 10' does NOT transfer"
          % (2 * c1["R"], DEG_PHI_EXPECTED),
          f1.e_min == 0 and not f1.degree_forcing_active,
          "E_min=%s" % f1.e_min)

    # C6 -- the two sub2 death modes are independent tests, not one condition.
    dead_b, _ = f2.cell_verdict(7, (3, 0, 0, 0))
    dead_d, _ = f2.cell_verdict(8, (1, 0, 0, 0))
    live, _ = f2.cell_verdict(7, (1, 1, 1, 0))
    check("C6 sub2 cell verdicts: a7_b3000 dead (b>=2), a8_b1000 dead (a+sum b=9), "
          "a7_b1110 alive", (not dead_b) and (not dead_d) and live)

    # C7 -- in sub1 the SAME cells behave differently: a8_b1000 survives (no degree
    #       forcing) while a7_b3000 still dies (b>=2 is window-independent).
    s1_deg, _ = f1.cell_verdict(8, (1, 0, 0, 0))
    s1_b, _ = f1.cell_verdict(7, (3, 0, 0, 0))
    check("C7 sub1: a8_b1000 SURVIVES (degree forcing vacuous) but a7_b3000 still "
          "dies (b>=2 is window-independent)", s1_deg and (not s1_b))

    # C8 -- the defect-0 step is NOT automatic: the sub2 universe really does carry
    #       states with deg_e = 10 and a+sum(b) < 10.
    try:
        U2 = _load("phase_d_states_sub2.json")
    except OSError:
        U2 = None
    if U2 is None:
        out.append("  [SKIP] C8/C9 (phase_d_states_sub2.json absent)")
    else:
        # NOTE (2026-07-25, this run): the two counts DIFFER, and that resolves the
        # 3844-vs-3790 discrepancy recorded in DIVISOR_SYZYGY.md sec.3b.
        #   deg_e != a+sum(b)                       -> 3844   (all defect states)
        #   deg_e == 10 and a+sum(b) < 10           -> 3790   (defect AT the cap)
        # The 54-state gap is defect-1 e at deg_e = 8 and 9 (cells a6_b1000_T1 and
        # a8_b0000_T1).  The external review's 3844 is the correct count of states
        # the universe carries with a free extra root; 3790 is a strict subset.
        defect = sum(1 for c in U2["cases"] for s in c["states"]
                     if s["deg_e"] != c["a_t"] + sum(c["b"]))
        defect10 = sum(1 for c in U2["cases"] for s in c["states"]
                       if s["deg_e"] == 10 and c["a_t"] + sum(c["b"]) < 10)
        check("C8 defect-0 is a REAL extra step: %d sub2 states carry a free extra "
              "root (deg_e != a+sum(b)); %d of them sit at the cap deg_e = 10, so "
              "the two readings differ by %d" % (defect, defect10, defect - defect10),
              defect == 3844 and defect10 == 3790,
              "defect=%d defect10=%d" % (defect, defect10))

        # C9 -- filter is idempotent and never resurrects a state.
        F = filter_universe(U2)
        FF = filter_universe(F)
        check("C9 filter is idempotent (F(F(U)) == F(U))",
              FF["state_total"] == F["state_total"]
              and FF["case_count"] == F["case_count"],
              "%d/%d vs %d/%d" % (FF["case_count"], FF["state_total"],
                                  F["case_count"], F["state_total"]))
        cen = F["divisor_filter"]["census"]
        check("C10 sub2 census is internally consistent "
              "(before = after + killed, cells and states)",
              cen["states_before"] == cen["states_after"] + cen["states_killed_b_ge_2"]
              + cen["states_killed_degree"] + cen["states_killed_defect"]
              and cen["cells_before"] == cen["cells_after"]
              + cen["cells_killed_b_ge_2"] + cen["cells_killed_degree"]
              + cen["cells_killed_defect"],
              json.dumps(cen, sort_keys=True))

    if verbose:
        print("divisor_filter.py -- e | Phi as a pipeline stage")
        print("  sub2 caps %s  E_min=%s" % (c2, f2.e_min))
        print("  sub1 caps %s  E_min=%s" % (c1, f1.e_min))
        print()
        for line in out:
            print(line)
        print("\n%d passed, %d failed" % (npass, nfail))
    elif nfail:
        for line in out:
            if "[FAIL]" in line:
                print(line)
    return nfail


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quiet", action="store_true", help="self-check only; exit 0 iff all pass")
    ap.add_argument("--universe", help="phase_d_states_*.json to filter")
    ap.add_argument("--out", help="write the filtered universe here")
    args = ap.parse_args()

    if args.universe:
        U = _load(args.universe)
        F = filter_universe(U)
        cen = F["divisor_filter"]["census"]
        if args.out:
            with open(os.path.join(HERE, args.out), "w", encoding="utf-8") as fh:
                json.dump(F, fh, indent=2, sort_keys=True)
                fh.write("\n")
            print("wrote %s" % args.out)
        print("%s: cells %d -> %d, flag cases %d -> %d, states %d -> %d"
              % (U["window"], cen["cells_before"], cen["cells_after"],
                 cen["flagcases_before"], cen["flagcases_after"],
                 cen["states_before"], cen["states_after"]))
        raise SystemExit(0)

    nfail = run(verbose=not args.quiet)
    raise SystemExit(1 if nfail else 0)


if __name__ == "__main__":
    main()
