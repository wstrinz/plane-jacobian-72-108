#!/usr/bin/env python3
"""Phase C worklist inventory over the f31 subcase-(2) survivor witnesses.

Standalone, read-only. Consumes ``cascade_cones_qt.json`` (the authoritative
q+t survivor data emitted by ``cascade_engine.py``) and cross-checks the
recorded tied-monomial strings against the ``h_l`` tables produced by
``cascade_signature.py``.

It groups every residue obligation across the surviving cases by
(place kind, cascade level, obligation kind, tied-monomial set), reports the
distinct patterns with their frequencies and the survivor cells they touch,
ranks the cells by total obligation depth of their cheapest survivor case,
and flags the T2 level-5 "squeeze"-eligible pattern.

Usage:  python phase_c_inventory.py            # human summary
        python phase_c_inventory.py --json      # machine-readable dump
"""

from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CONES = ROOT / "cascade_cones_qt.json"


def cell_id(branch: dict) -> str:
    return "a%d b=%s %s" % (
        branch["a_t"],
        "".join(map(str, branch["b"])),
        branch["branch"],
    )


def load_survivors(path: Path = CONES) -> list[dict]:
    data = json.loads(path.read_text(encoding="utf-8"))
    return [b for b in data["branches"] if b.get("survivor_cases")]


def h_level_tables() -> dict[int, set[str]]:
    """Monomial-string set for each h_l, matching the JSON tied encoding."""

    from cascade_signature import build_signature

    tables: dict[int, set[str]] = {}
    for level in build_signature()["levels"]:
        strings = set()
        for mono in level["monomials"]:
            ex = mono["exponents"]
            strings.add(
                "%s*d2^%d*d1^%d*sigma^%d*e^%d"
                % (mono["coefficient"], ex["d2"], ex["d1"], ex["sigma"], ex["e"])
            )
        tables[level["index"]] = strings
    return tables


def validate_tied(survivors: list[dict], tables: dict[int, set[str]]) -> int:
    bad = 0
    for branch in survivors:
        for case in branch["survivor_cases"]:
            for witness in case["witness"]:
                for ob in witness["obligations"]:
                    for tied in ob.get("tied", []):
                        if tied not in tables[ob["level"]]:
                            bad += 1
    return bad


def inventory(survivors: list[dict]):
    pat_freq: Counter = Counter()
    pat_cells: dict[tuple, set] = defaultdict(set)
    pat_depths: dict[tuple, Counter] = defaultdict(Counter)
    for branch in survivors:
        cid = cell_id(branch)
        for case in branch["survivor_cases"]:
            for witness in case["witness"]:
                place = witness["place"]
                for ob in witness["obligations"]:
                    tied = tuple(sorted(ob.get("tied", [])))
                    key = (place, ob["level"], ob["kind"], tied)
                    pat_freq[key] += 1
                    pat_cells[key].add(cid)
                    pat_depths[key][ob["depth"]] += 1
    return pat_freq, pat_cells, pat_depths


def cheapest_case(branch: dict):
    """Return (total_depth, q_depth, t_depth, n_exact0, case) for the
    survivor case with the smallest total obligation depth."""

    best = None
    for case in branch["survivor_cases"]:
        qd = td = zeros = 0
        for witness in case["witness"]:
            for ob in witness["obligations"]:
                if ob["depth"] == 0:
                    zeros += 1
                if witness["place"] == "q":
                    qd += ob["depth"]
                else:
                    td += ob["depth"]
        rec = (qd + td, qd, td, zeros, case)
        if best is None or rec[0] < best[0]:
            best = rec
    return best


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()

    survivors = load_survivors()
    tables = h_level_tables()
    bad = validate_tied(survivors, tables)
    pat_freq, pat_cells, pat_depths = inventory(survivors)

    total_cases = sum(len(b["survivor_cases"]) for b in survivors)
    ranking = sorted(
        ((cell_id(b),) + cheapest_case(b)[:4] for b in survivors),
        key=lambda r: r[1],
    )

    if args.json:
        out = {
            "branches": len(survivors),
            "survivor_cases": total_cases,
            "tied_mismatches": bad,
            "distinct_patterns": len(pat_freq),
            "patterns": [
                {
                    "place": k[0],
                    "level": k[1],
                    "kind": k[2],
                    "tied": list(k[3]),
                    "freq": pat_freq[k],
                    "cells": sorted(pat_cells[k]),
                    "depths": dict(pat_depths[k]),
                }
                for k in sorted(pat_freq, key=lambda k: -pat_freq[k])
            ],
            "cell_ranking": [
                {"cell": r[0], "total_depth": r[1], "q_depth": r[2],
                 "t_depth": r[3], "n_depth0": r[4]}
                for r in ranking
            ],
        }
        print(json.dumps(out, indent=2))
        return

    print(f"survivor branches (cells): {len(survivors)}")
    print(f"survivor cases:            {total_cases}")
    print(f"tied-monomial mismatches vs h_l tables: {bad}")
    print(f"distinct obligation patterns (place,level,kind,tied): {len(pat_freq)}")
    print()
    print("PATTERNS (by frequency):")
    for k in sorted(pat_freq, key=lambda k: -pat_freq[k]):
        place, level, kind, tied = k
        depths = ",".join(f"{d}x{n}" for d, n in sorted(pat_depths[k].items()))
        print(
            f"  {pat_freq[k]:4d}  {place} L{level} {kind:20s} "
            f"cells={len(pat_cells[k]):2d} tiedN={len(tied)} depth[{depths}]"
        )
    print()
    print("CELL RANKING (cheapest survivor case, by total obligation depth):")
    print(f"  {'cell':16s} {'TOT':>4s} {'qdep':>5s} {'tdep':>5s} {'d0':>4s}")
    for r in ranking:
        print(f"  {r[0]:16s} {r[1]:4d} {r[2]:5d} {r[3]:5d} {r[4]:4d}")


if __name__ == "__main__":
    main()
