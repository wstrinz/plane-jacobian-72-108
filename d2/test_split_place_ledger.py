#!/usr/bin/env python3
"""Independent finite checks for split_place_ledger.py."""

from __future__ import annotations

from itertools import product

import split_place_ledger as ledger

def _require(_cond, _msg):
    """Proof-critical check: fails loudly and exits nonzero, unaffected by python -O."""
    if not _cond:
        import sys as _sys
        print("FAIL: " + str(_msg))
        _sys.exit(1)



def brute_terminal(a: int, b: tuple[int, int, int, int], branch: str) -> bool:
    """Search every admissible local auxiliary valuation independently."""

    if branch == "T1":
        level, degree_cap = 7, 6
    else:
        level, degree_cap = 6, 8
    g_cap = 10 + 3 * a
    for orders in product(range(degree_cap + 1), repeat=4):
        if sum(orders) > degree_cap:
            continue
        g_orders = tuple(
            level + 2 * order - 3 * bi for order, bi in zip(orders, b)
        )
        if min(g_orders) >= 0 and sum(g_orders) <= g_cap:
            return True
    return False


def main() -> None:
    rows = ledger.build_ledger()
    stats = ledger.summarize(rows)

    _require(len(rows) == 327, "len(rows) == 327")
    _require(sum(len(ledger.sorted_q_vectors(10 - a)) for a in range(11)) == 327, "sum(len(ledger.sorted_q_vectors(10 - a)) for a in range(11)) == 327")

    for row in rows:
        a = row["a_t"]
        b = tuple(row["b"])
        _require(tuple(sorted(b, reverse=True)) == b, "tuple(sorted(b, reverse=True)) == b")
        _require(a + sum(b) + row["residual_degree_budget"] == 10, "a + sum(b) + row[\"residual_degree_budget\"] == 10")
        for branch in ("T1", "T2"):
            exact = ledger.terminal_test(a, b, branch=branch)
            brute = brute_terminal(a, b, branch)
            _require(exact.feasible == brute, (a, b, branch, exact, brute))
            _require(all(
                3 * bi + gi == exact.source_level + 2 * xi
                for bi, gi, xi in zip(
                    b,
                    exact.resulting_g_orders,
                    exact.minimum_auxiliary_orders,
                )
            ), "all( 3 * bi + gi == exact.source_level + 2 * xi for bi, gi, xi in zip( b, exact.resulting_g_orders, exact.minimum_auxiliary_orders, ) )")

    _require(stats == {
        "raw_strata": 327,
        "old_uniform_strata": 21,
        "partial_support_strata": 306,
        "T1_terminal_feasible": 197,
        "T2_terminal_feasible": 246,
        "strata_killed_by_terminal_both_branches": 81,
        "partial_support_strata_killed_by_terminal_both_branches": 75,
        "open_strata_after_terminal_and_existing_proofs": 235,
        "open_partial_support_strata": 231,
        "open_branches_after_terminal_and_existing_proofs": 420,
    }, "stats == { \"raw_strata\": 327, \"old_uniform_strata\": 21, \"partial_support_strata\": 306, \"T1_terminal_feasible\": 197, \"T2_terminal_feasible\": 246, \"strata_killed_by_terminal_both_branches\": 81, \"partial_support_strata_killed_by_terminal_both_branches\": 75, \"open_strata_after_terminal_and_existing_proofs\": 235, \"open_partial_support_strata\": 231, \"open_branches_after_terminal_and_existing_proofs\": 420, }")

    uniform_open = [
        (row["a_t"], tuple(row["b"]), tuple(row["open_branches"]))
        for row in rows
        if len(set(row["b"])) == 1 and row["open_branches"]
    ]
    _require(uniform_open == [
        (6, (1, 1, 1, 1), ("T1",)),
        (8, (0, 0, 0, 0), ("T1", "T2")),
        (9, (0, 0, 0, 0), ("T1",)),
        (10, (0, 0, 0, 0), ("T1",)),
    ], "uniform_open == [ (6, (1, 1, 1, 1), (\"T1\",)), (8, (0, 0, 0, 0), (\"T1\", \"T2\")), (9, (0, 0, 0, 0), (\"T1\",)), (10, (0, 0, 0, 0), (\"T1\",)), ]")

    print("split-place ledger: PASS")
    print("  327 raw = 21 old uniform + 306 partial-support strata")
    print("  terminal levels kill 81 strata, including 75 partial-support")
    print("  235 strata / 420 live branches remain after scoped old proofs")


if __name__ == "__main__":
    main()
