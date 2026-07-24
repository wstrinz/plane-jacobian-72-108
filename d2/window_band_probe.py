#!/usr/bin/env python3
"""window_band_probe.py  --  adjudicate the staircase claim on a COMPUTABLE system

BIGRADED_PROBE.md claims the walled R9 z=1 H-system has a bandwidth-2 staircase:
sweeping y-order inward from an extreme slice introduces "exactly two new spare
unknowns per slice".  Its own verifier FAILED at slice 1 (4 new spares, not 2),
and the Milestone-1 annotator's R3 lane -- built to adjudicate this -- does not
finish: expanding H2,H3,H5 with cubic spare terms is the documented sympy
intractability trap (>1h with no output).

This module sidesteps that.  The staircase claim is a statement about the
BIGRADED LATTICE (which spare coefficients touch which y-order slice), and that
lattice has the same shape on the sub2 home G-system -- which expands in seconds.
So we measure the increment there, on states we can actually compute, and get a
decisive answer to the structural question without waiting on R9.

It reuses `bigrade_annotator.r3_band_structure` verbatim -- the same measurement
R3 would have made -- so the numbers are directly comparable when R9 is
eventually reached by the Singular route.

Also reports, per CAOS's transplantable lesson (their EXP-037 hit our exact
failure mode): whether the MINIMAL or MAXIMAL slice is the right sweep anchor,
i.e. which end has the small, tightly-coupled face.

Read-only.  Writes `window_band_probe.json` only.

Usage:
    python window_band_probe.py                    # default representative state
    python window_band_probe.py --a-t 10 --deg-d1 2 --deg-d2 4 --deg-e 10 --deg-sigma 8
"""
from __future__ import annotations

import argparse
import json
import os
from collections import defaultdict

import sympy as sp

import bigrade_annotator as ba
import face_kill_sweep as F
from bigrade_annotator import y

HERE = os.path.dirname(os.path.abspath(__file__)) or "."
OUT_PATH = os.path.join(HERE, "window_band_probe.json")


def build_generic_state(deg_d1, deg_d2, deg_e, deg_sigma):
    """Generic sub2 state polynomials at the requested (exact) degrees."""
    def poly(prefix, deg):
        if deg is None or deg < 0:
            return sp.Integer(0)
        return sum(sp.Symbol("%s%d" % (prefix, i)) * y**i for i in range(deg + 1))

    return {
        "d2": poly("a", deg_d2),
        "d1": poly("b", deg_d1),
        "sigma": poly("s", deg_sigma),
        "e": sp.expand(sp.Symbol("gamma") * (y + 1)**int(deg_e)),
    }


def per_generator_staircase(system):
    """New-spare increment sweeping each generator's own y-order range upward."""
    by_gen = defaultdict(lambda: defaultdict(set))
    for e in system.eqs:
        by_gen[e.label[0]][e.nu] |= {str(s) for s in e.vs}

    out = {}
    for gname, slices in by_gen.items():
        seen, rows = set(), []
        for nu in sorted(slices):
            new = slices[nu] - seen
            seen |= slices[nu]
            rows.append({"nu": nu, "present": len(slices[nu]), "new": len(new)})
        increments = sorted({r["new"] for r in rows})
        out[gname] = {
            "n_slices": len(rows),
            "nu_range": [rows[0]["nu"], rows[-1]["nu"]] if rows else None,
            "increments_observed": increments,
            "head": rows[:6],
            "tail": rows[-3:],
        }
    return out


def anchor_comparison(system):
    """Which END of the y-order axis carries the small, tightly-coupled face?

    CAOS's EXP-037 found the naive extreme was the WRONG anchor and the correct
    diagonal block sits at the MINIMAL class.  We check both ends here.
    """
    by_nu = defaultdict(set)
    eq_count = defaultdict(int)
    for e in system.eqs:
        by_nu[e.nu] |= {str(s) for s in e.vs}
        eq_count[e.nu] += 1
    nus = sorted(by_nu)
    ends = {}
    for tag, nu in (("min", nus[0]), ("max", nus[-1])):
        ends[tag] = {"nu": nu, "n_eqs": eq_count[nu], "n_spares": len(by_nu[nu]),
                     "spares": sorted(by_nu[nu])}
    return ends


def probe(deg_d1, deg_d2, deg_e, deg_sigma, verbose=True):
    polys = build_generic_state(deg_d1, deg_d2, deg_e, deg_sigma)
    label = "sub2 generic  d1=%s d2=%s e=%s sigma=%s" % (deg_d1, deg_d2, deg_e, deg_sigma)
    system = F.build_state_system(label, polys, "sub2")

    viol, checked = ba.check_bigrade_consistency(system)
    if viol:
        raise SystemExit("FAIL: bigrade consistency violated (%d)" % len(viol))

    band = ba.r3_band_structure(system)          # the SAME function R3 would run
    per_gen = per_generator_staircase(system)
    ends = anchor_comparison(system)

    if verbose:
        print("=" * 78)
        print("WINDOW BAND PROBE  --  %s" % label)
        print("=" * 78)
        print("  bigrade consistency: %d monomials, 0 violations" % checked)
        print("  vars=%d  eqs=%d  spare caps=%s"
              % (len(system.vars), len(system.eqs), system.extra["spare_caps"]))
        print()
        print("  GLOBAL sweep (all generators on one y-order axis):")
        print("    u-axis values (generator weights):", band["u_axis_values"])
        print("    per-slice NEW-spare increments observed:", band["increments_seen"])
        print("    staircase head:")
        for row in band["staircase_head"]:
            print("       nu=%-4d present=%-3d new=%-3d %s" % row)
        print()
        print("  PER-GENERATOR sweep (each generator's own range):")
        for g in sorted(per_gen):
            d = per_gen[g]
            print("    %-3s slices=%-3d nu=%s  increments observed=%s"
                  % (g, d["n_slices"], d["nu_range"], d["increments_observed"]))
            print("        head:", ", ".join("nu=%d(+%d)" % (r["nu"], r["new"])
                                             for r in d["head"]))
        print()
        print("  ANCHOR comparison (which end is the small face?):")
        for tag in ("min", "max"):
            e = ends[tag]
            print("    %-3s slice nu=%-4d  %2d eqs  %2d spares  %s"
                  % (tag, e["nu"], e["n_eqs"], e["n_spares"],
                     ",".join(e["spares"][:6]) or "(none)"))
        print()
        inc = band["increments_seen"]
        nontrivial = [i for i in inc if i > 0]
        print("  VERDICT on the BIGRADED_PROBE bandwidth-2 staircase claim:")
        if nontrivial == [2]:
            print("    increments are uniformly +2 -- claim SUPPORTED on this system")
        else:
            print("    increments are %s, NOT uniformly +2 -- claim NOT SUPPORTED"
                  % nontrivial)

    return {"state": {"deg_d1": deg_d1, "deg_d2": deg_d2,
                      "deg_e": deg_e, "deg_sigma": deg_sigma},
            "n_vars": len(system.vars), "n_eqs": len(system.eqs),
            "spare_caps": system.extra["spare_caps"],
            "global_band": band, "per_generator": per_gen, "anchors": ends}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--deg-d1", type=int, default=2)
    ap.add_argument("--deg-d2", type=int, default=4)
    ap.add_argument("--deg-e", type=int, default=10)
    ap.add_argument("--deg-sigma", type=int, default=8)
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    out = probe(args.deg_d1, args.deg_d2, args.deg_e, args.deg_sigma,
                verbose=not args.quiet)
    with open(OUT_PATH, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2, sort_keys=True, default=str)
    print("wrote", os.path.basename(OUT_PATH))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
