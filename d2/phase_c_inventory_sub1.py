#!/usr/bin/env python3
"""Phase C worklist inventory over the f31 subcase-(1) survivor witnesses.

Standalone, read-only. Consumes ``cascade_cones_sub1_qt.json`` (the
authoritative q+t survivor data for subcase (1) emitted by
``cascade_engine.py --window sub1``) and cross-checks the recorded
tied-monomial strings against the ``h_l`` tables produced by
``cascade_signature.py`` (window-independent).

It groups every residue obligation across the surviving cases by
(place kind, cascade level, obligation kind, tied-monomial set), reports the
distinct patterns with frequencies and touched cells, and -- the point of the
sub1 port -- measures pattern OVERLAP against the subcase-(2) frontier
(``cascade_cones_qt.json``, the 41 patterns of ``PHASE_C_WORKLIST.md``):

  identical    : same (place,level,kind,tied) key AND same depth-set;
  same_key     : same key, but sub1 needs additional (usually deeper) depths;
  new          : key absent from sub2 -- split into "reuses a sub2 tied-set"
                 (only place/kind differs) vs "genuinely new residue poly".

It also verifies the a-invariance of the 26-branch survivor family (a<=8) and
ranks cells by total obligation depth of their cheapest survivor case.

Usage:  python phase_c_inventory_sub1.py            # human summary
        python phase_c_inventory_sub1.py --json      # machine-readable dump
"""

from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CONES_SUB1 = ROOT / "cascade_cones_sub1_qt.json"
CONES_SUB2 = ROOT / "cascade_cones_qt.json"


def cell_id(branch: dict) -> str:
    return "a%d b=%s %s" % (
        branch["a_t"],
        "".join(map(str, branch["b"])),
        branch["branch"],
    )


def load_survivors(path: Path) -> list[dict]:
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
    """Return freq / cells / depth-multiset / depth-set keyed by pattern."""

    pat_freq: Counter = Counter()
    pat_cells: dict[tuple, set] = defaultdict(set)
    pat_depths: dict[tuple, Counter] = defaultdict(Counter)
    pat_depthset: dict[tuple, set] = defaultdict(set)
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
                    pat_depthset[key].add(ob["depth"])
    return pat_freq, pat_cells, pat_depths, pat_depthset


def cheapest_case(branch: dict):
    """(total_depth, q_depth, t_depth, n_depth0, case) minimising total depth."""

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


def family_a_invariance(survivors: list[dict]):
    """Return {a: sorted set of (b,branch) pairs} and the a<=8 invariance flag."""

    fam: dict[int, set] = defaultdict(set)
    for b in survivors:
        fam[b["a_t"]].add((tuple(b["b"]), b["branch"]))
    base = fam.get(0, set())
    invariant = all(fam[a] == base for a in range(0, 9) if a in fam)
    return fam, base, invariant


def overlap(sub1_keys, sub1_depthset, sub2_keys, sub2_depthset, tables):
    k1, k2 = set(sub1_keys), set(sub2_keys)
    common = k1 & k2
    new = k1 - k2
    identical = {k for k in common if sub1_depthset[k] == sub2_depthset[k]}
    same_key = common - identical
    tied2 = {k[3] for k in k2 if k[3]}
    new_reuse = {k for k in new if k[3] and k[3] in tied2}
    new_poly = {k for k in new if not (k[3] and k[3] in tied2)}
    return {
        "sub2_keys": len(k2),
        "sub2_only": len(k2 - k1),
        "common": len(common),
        "identical": len(identical),
        "same_key_diff_depth": len(same_key),
        "new": len(new),
        "new_reusing_sub2_tied_set": len(new_reuse),
        "new_genuine_residue_poly": len(new_poly),
        "_new_poly_keys": sorted(new_poly),
        "_same_key_keys": sorted(same_key),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()

    sub1 = load_survivors(CONES_SUB1)
    sub2 = load_survivors(CONES_SUB2)
    tables = h_level_tables()
    bad = validate_tied(sub1, tables)

    f1, c1, d1, ds1 = inventory(sub1)
    f2, c2, d2, ds2 = inventory(sub2)
    ov = overlap(f1, ds1, f2, ds2, tables)

    fam, base, invariant = family_a_invariance(sub1)
    total_cases = sum(len(b["survivor_cases"]) for b in sub1)
    ranking = sorted(
        ((cell_id(b),) + cheapest_case(b)[:4] for b in sub1),
        key=lambda r: (r[1], -r[4]),
    )

    by_kind = Counter()
    by_place = Counter()
    by_level = Counter()
    for k, n in f1.items():
        by_place[k[0]] += n
        by_level[k[1]] += n
        by_kind[k[2]] += n

    if args.json:
        out = {
            "branches": len(sub1),
            "survivor_cases": total_cases,
            "tied_mismatches": bad,
            "distinct_patterns": len(f1),
            "by_kind": dict(by_kind),
            "by_place": dict(by_place),
            "by_level": dict(by_level),
            "overlap_with_sub2": {k: v for k, v in ov.items()
                                  if not k.startswith("_")},
            "new_genuine_residue_polys": [
                {"place": k[0], "level": k[1], "kind": k[2], "tied": list(k[3])}
                for k in ov["_new_poly_keys"]
            ],
            "family_26_a_invariant_for_a_le_8": invariant,
            "family_sizes_by_a": {a: len(fam[a]) for a in sorted(fam)},
            "patterns": [
                {
                    "place": k[0], "level": k[1], "kind": k[2],
                    "tied": list(k[3]), "freq": f1[k],
                    "cells": len(c1[k]), "depths": dict(d1[k]),
                }
                for k in sorted(f1, key=lambda k: -f1[k])
            ],
            "cell_ranking": [
                {"cell": r[0], "total_depth": r[1], "q_depth": r[2],
                 "t_depth": r[3], "n_depth0": r[4]}
                for r in ranking
            ],
        }
        print(json.dumps(out, indent=2))
        return

    print(f"survivor branches (cells):        {len(sub1)}")
    print(f"survivor cases:                   {total_cases}")
    print(f"tied-monomial mismatches vs h_l:  {bad}")
    print(f"distinct obligation patterns:     {len(f1)}")
    print(f"  by kind:  {dict(by_kind)}")
    print(f"  by place: {dict(by_place)}")
    print(f"  by level: {dict(by_level)}")
    print()
    print("OVERLAP WITH SUB2 (41 patterns of PHASE_C_WORKLIST.md):")
    for key in ("sub2_keys", "sub2_only", "common", "identical",
                "same_key_diff_depth", "new",
                "new_reusing_sub2_tied_set", "new_genuine_residue_poly"):
        print(f"  {key:30s} {ov[key]}")
    print()
    print(f"26-family a-invariant for a<=8:   {invariant}")
    print(f"  family sizes by a: {{{', '.join('%d:%d' % (a, len(fam[a])) for a in sorted(fam))}}}")
    print()
    print("TOP 15 PATTERNS (by frequency):")
    for k in sorted(f1, key=lambda k: -f1[k])[:15]:
        place, level, kind, tied = k
        depths = ",".join(f"{d}x{n}" for d, n in sorted(d1[k].items()))
        print(
            f"  {f1[k]:5d}  {place} L{level} {kind:20s} "
            f"cells={len(c1[k]):3d} tiedN={len(tied)} depth[{depths[:34]}]"
        )
    print()
    print("CELL RANKING (cheapest survivor case, first 20):")
    print(f"  {'cell':16s} {'TOT':>4s} {'qdep':>5s} {'tdep':>5s} {'d0':>4s}")
    for r in ranking[:20]:
        print(f"  {r[0]:16s} {r[1]:4d} {r[2]:5d} {r[3]:5d} {r[4]:4d}")


if __name__ == "__main__":
    main()
