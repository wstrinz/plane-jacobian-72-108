#!/usr/bin/env python3
"""Generate the geometric q-place ledger for f31 subcase (1).

All cascade parameters are imported from sub1_cascade_verify.py, which executes
its SymPy assertions on import. No subcase-(2) proof table is used.
"""
from __future__ import annotations
from dataclasses import asdict, dataclass
from itertools import product
import json
from pathlib import Path
from typing import Any
import sympy as sp

from sub1_cascade_verify import (
    ALTERNATE_MIN_A, D1_DEG_CAP, E_DEG_CAP, MAX_A, Q_ROOT_COUNT,
    SCHEMA_VERSION, SIGMA_DEG_CAP, STANDARD_MAX_A, T1_G_DEG_CAP, T1_LEVEL,
    T2_G_DEG_CAP, T2_LEVEL, T3_REFERENCE, T3_STATUS, TERMINAL_AUX_POWER,
    TERMINAL_E_POWER,
)

ROOT = Path(__file__).resolve().parent
JSON_PATH = ROOT / "split_place_ledger_sub1.json"
MD_PATH = ROOT / "SPLIT_PLACE_LEDGER_SUB1.md"
BRANCHES = ("T1", "T2")


@dataclass(frozen=True)
class TerminalTest:
    branch: str
    source_level: int
    auxiliary: str
    auxiliary_degree_cap: int
    g_degree_cap: int
    minimum_auxiliary_orders: tuple[int, int, int, int]
    resulting_g_orders: tuple[int, int, int, int]
    minimum_auxiliary_degree: int
    minimum_g_degree: int
    feasible: bool


def sorted_q_vectors(total_cap: int) -> list[tuple[int, int, int, int]]:
    return [
        b for b in product(range(total_cap + 1), repeat=Q_ROOT_COUNT)
        if sum(b) <= total_cap and tuple(sorted(b, reverse=True)) == b
    ]


def minimum_order(rhs_q_order: int, b: int) -> int:
    """Least x>=0 with 3b+v(g)=level+2x and v(g)>=0."""
    value = max(0, (TERMINAL_E_POWER*b - rhs_q_order + TERMINAL_AUX_POWER - 1)
                // TERMINAL_AUX_POWER)
    assert sp.Integer(rhs_q_order + TERMINAL_AUX_POWER*value - TERMINAL_E_POWER*b) >= 0
    if value:
        assert sp.Integer(rhs_q_order + TERMINAL_AUX_POWER*(value-1) - TERMINAL_E_POWER*b) < 0
    return value


def terminal_test(a: int, b: tuple[int, int, int, int], *, branch: str) -> TerminalTest:
    """Apply the corrected, exact terminal valuation/degree conditions."""
    assert sp.Integer(a) <= STANDARD_MAX_A
    if branch == "T1":
        level, auxiliary = T1_LEVEL, "d1"
        auxiliary_cap, g_cap = D1_DEG_CAP, T1_G_DEG_CAP
    elif branch == "T2":
        level, auxiliary = T2_LEVEL, "sigma"
        auxiliary_cap, g_cap = SIGMA_DEG_CAP, T2_G_DEG_CAP
    else:
        raise ValueError(f"unknown terminal branch: {branch}")
    aux_orders = tuple(minimum_order(level, bi) for bi in b)
    g_orders = tuple(level + TERMINAL_AUX_POWER*xi - TERMINAL_E_POWER*bi
                     for xi, bi in zip(aux_orders, b))
    assert all(sp.Integer(order) >= 0 for order in g_orders)
    aux_degree, g_degree = sum(aux_orders), sum(g_orders)
    feasible = aux_degree <= auxiliary_cap and g_degree <= g_cap
    assert sp.Equivalent(
        sp.sympify(feasible),
        sp.And(sp.Integer(aux_degree) <= auxiliary_cap,
               sp.Integer(g_degree) <= g_cap),
    ) is sp.true
    return TerminalTest(
        branch, level, auxiliary, auxiliary_cap, g_cap, aux_orders, g_orders,
        aux_degree, g_degree, feasible,
    )


def branch_record(a: int, b: tuple[int, int, int, int], branch: str) -> dict[str, Any]:
    if a >= ALTERNATE_MIN_A:
        assert sp.Integer(a) > STANDARD_MAX_A
        return {
            "status": "alternate_regime_open",
            "reference": "v=30-3a<0; standard polynomial cascade reduction is unavailable",
            "terminal": None,
        }
    terminal = terminal_test(a, b, branch=branch)
    if terminal.feasible:
        status = "open_after_terminal"
        reference = "requires lower cascade levels or another exact argument"
    else:
        status = "proven_infeasible_terminal"
        reference = (
            f"split-place level {terminal.source_level}: "
            f"min deg {terminal.auxiliary}={terminal.minimum_auxiliary_degree} "
            f"(cap {terminal.auxiliary_degree_cap}), "
            f"min deg g{terminal.source_level}={terminal.minimum_g_degree} "
            f"(cap {terminal.g_degree_cap})"
        )
    return {"status": status, "reference": reference, "terminal": asdict(terminal)}


def support_type(b: tuple[int, int, int, int]) -> str:
    if b == (0,)*Q_ROOT_COUNT:
        return "geometrically_q_coprime"
    return "uniform_q_power" if len(set(b)) == 1 else "partial_q_support"


def build_ledger() -> list[dict[str, Any]]:
    ledger = []
    for a in range(MAX_A + 1):
        for b in sorted_q_vectors(E_DEG_CAP - a):
            branches = {branch: branch_record(a, b, branch) for branch in BRANCHES}
            open_branches = [
                branch for branch, data in branches.items()
                if data["status"] in ("open_after_terminal", "alternate_regime_open")
            ]
            if a >= ALTERNATE_MIN_A:
                stratum_status = "alternate_regime_open"
            else:
                stratum_status = "open_after_terminal" if open_branches else "proven_infeasible"
            row = {
                "a_t": a,
                "b": list(b),
                "q_multiplicity_sum": sum(b),
                "residual_degree_budget": E_DEG_CAP - a - sum(b),
                "support_type": support_type(b),
                "T3": {"status": T3_STATUS, "reference": T3_REFERENCE},
                "branches": branches,
                "open_branches": open_branches,
                "stratum_status": stratum_status,
            }
            assert sp.Integer(a + row["q_multiplicity_sum"] + row["residual_degree_budget"]) == E_DEG_CAP
            assert sp.Integer(row["residual_degree_budget"]) >= 0
            ledger.append(row)
    assert sp.Integer(len(ledger)) == 1333
    return ledger


def summarize(ledger: list[dict[str, Any]]) -> dict[str, int]:
    standard = [row for row in ledger if row["a_t"] <= STANDARD_MAX_A]
    alternate = [row for row in ledger if row["a_t"] >= ALTERNATE_MIN_A]
    terminal_killed = [
        row for row in standard
        if all(data["status"] == "proven_infeasible_terminal"
               for data in row["branches"].values())
    ]
    open_rows = [row for row in ledger if row["open_branches"]]
    stats = {
        "raw_strata": len(ledger),
        "old_uniform_strata": sum(len(set(row["b"])) == 1 for row in ledger),
        "partial_support_strata": sum(len(set(row["b"])) > 1 for row in ledger),
        "T1_terminal_feasible": sum(row["branches"]["T1"]["terminal"]["feasible"] for row in standard),
        "T2_terminal_feasible": sum(row["branches"]["T2"]["terminal"]["feasible"] for row in standard),
        "strata_killed_by_terminal_both_branches": len(terminal_killed),
        "partial_support_strata_killed_by_terminal_both_branches": sum(row["support_type"] == "partial_q_support" for row in terminal_killed),
        "open_strata_after_terminal_and_existing_proofs": len(open_rows),
        "open_partial_support_strata": sum(row["support_type"] == "partial_q_support" for row in open_rows),
        "open_branches_after_terminal_and_existing_proofs": sum(data["status"] == "open_after_terminal" for row in standard for data in row["branches"].values()),
        "alternate_regime_strata": len(alternate),
        "alternate_regime_open_branches": sum(data["status"] == "alternate_regime_open" for row in alternate for data in row["branches"].values()),
    }
    expected = {
        "raw_strata": 1333, "old_uniform_strata": 40,
        "partial_support_strata": 1293, "T1_terminal_feasible": 1007,
        "T2_terminal_feasible": 1171,
        "strata_killed_by_terminal_both_branches": 136,
        "partial_support_strata_killed_by_terminal_both_branches": 136,
        "open_strata_after_terminal_and_existing_proofs": 1197,
        "open_partial_support_strata": 1157,
        "open_branches_after_terminal_and_existing_proofs": 2178,
        "alternate_regime_strata": 26, "alternate_regime_open_branches": 52,
    }
    for key, value in expected.items():
        assert sp.Integer(stats[key]) == value
    return stats


def render_markdown(ledger: list[dict[str, Any]], stats: dict[str, int]) -> str:
    by_a = []
    for a in range(MAX_A + 1):
        rows = [row for row in ledger if row["a_t"] == a]
        killed = sum(row["stratum_status"] == "proven_infeasible" for row in rows)
        opened = sum(bool(row["open_branches"]) for row in rows)
        alternate = sum(row["stratum_status"] == "alternate_regime_open" for row in rows)
        assert sp.Integer(len(rows)) == killed + opened
        by_a.append((a, len(rows), killed, opened, alternate))
    expected_by_a = [
        (0,295,77,218,0),(1,241,38,203,0),(2,194,15,179,0),
        (3,155,5,150,0),(4,121,1,120,0),(5,94,0,94,0),
        (6,71,0,71,0),(7,53,0,53,0),(8,38,0,38,0),
        (9,27,0,27,0),(10,18,0,18,0),(11,12,0,12,12),
        (12,7,0,7,7),(13,4,0,4,4),(14,2,0,2,2),(15,1,0,1,1),
    ]
    for actual, expected in zip(by_a, expected_by_a):
        assert all(sp.Integer(x) == y for x, y in zip(actual, expected))
    support_rows = []
    for name in ("geometrically_q_coprime", "uniform_q_power", "partial_q_support"):
        rows = [row for row in ledger if row["support_type"] == name]
        killed = sum(row["stratum_status"] == "proven_infeasible" for row in rows)
        opened = sum(bool(row["open_branches"]) for row in rows)
        alternate = sum(row["stratum_status"] == "alternate_regime_open" for row in rows)
        assert sp.Integer(len(rows)) == killed + opened
        support_rows.append((name, len(rows), killed, opened, alternate))
    expected_support = [
        ('geometrically_q_coprime',16,0,16,5),
        ('uniform_q_power',24,0,24,1),
        ('partial_q_support',1293,136,1157,20),
    ]
    for actual, expected in zip(support_rows, expected_support):
        assert actual[0] == expected[0]
        assert all(sp.Integer(x) == y for x, y in zip(actual[1:], expected[1:]))
    uniform_open = [row for row in ledger if len(set(row["b"])) == 1 and row["open_branches"]]
    lines = [
        "# Geometric split-place ledger — f31 subcase (1)", "",
        "**Generated by `split_place_ledger_sub1.py`; do not hand-edit counts.**", "",
        "After base change write `q=p1*p2*p3*p4`, `a=v_t(e)`, and",
        "`b_i=v_{p_i}(e/t^a)` in decreasing order. The cap is `a+sum(b_i)<=15`.", "",
        "## Exact terminal pruning", "",
        "T1: `3b_i+v_i(g7)=7+2v_i(d1)`, with `deg d1<=9` and `deg g7<=46`.",
        "T2: `3b_i+v_i(g6)=6+2v_i(sigma)`, with `deg sigma<=12` and `deg g6<=48`.",
        "These corrected terminal caps replace the unsupported proposed bound `deg g_l<=15+3a`.",
        "Terminal tests apply only for `a<=10`; rows `a>=11` are `alternate_regime_open`.", "",
        "## Counts", "",
        f"- Raw geometric strata: **{stats['raw_strata']}**.",
        f"- Old uniform strata: **{stats['old_uniform_strata']}**; partial-support strata: **{stats['partial_support_strata']}**.",
        f"- T1 terminal-feasible: **{stats['T1_terminal_feasible']}**; T2 terminal-feasible: **{stats['T2_terminal_feasible']}**.",
        f"- Killed at terminal levels in both live branches: **{stats['strata_killed_by_terminal_both_branches']}**, including **{stats['partial_support_strata_killed_by_terminal_both_branches']}** partial-support strata.",
        f"- Open after terminal pruning (no existing subcase-(1) proofs imported): **{stats['open_strata_after_terminal_and_existing_proofs']}** strata.",
        f"- Standard-regime open branches: **{stats['open_branches_after_terminal_and_existing_proofs']}**; alternate-regime strata/branches: **{stats['alternate_regime_strata']}**/**{stats['alternate_regime_open_branches']}**.", "",
        "T3 (`d1=0`, `sigma=0`) is excluded globally by the subcase-(1)",
        "split-place degree enumeration and Mason–Stothers margins checked in `sub1_cascade_verify.py`.", "",
        "## By t-multiplicity", "",
        "| `a` | raw | terminal-killed | open | alternate-regime |",
        "|---:|---:|---:|---:|---:|",
    ]
    lines.extend(f"| {a} | {raw} | {killed} | {opened} | {alternate} |" for a, raw, killed, opened, alternate in by_a)
    lines.extend(["", "## By support type", "",
                  "| support type | raw | terminal-killed | open | alternate-regime |",
                  "|:---|---:|---:|---:|---:|"])
    lines.extend(f"| `{name}` | {raw} | {killed} | {opened} | {alternate} |" for name, raw, killed, opened, alternate in support_rows)
    lines.extend(["", "## Uniform frontier", "",
                  "| `a` | `(b1,b2,b3,b4)` | open branches |",
                  "|---:|:---:|:---|"])
    lines.extend(f"| {row['a_t']} | `{tuple(row['b'])}` | {', '.join(row['open_branches'])} |" for row in uniform_open)
    lines.extend(["", "The complete valuation witnesses and references are in",
                  "`split_place_ledger_sub1.json`.", ""])
    return "\n".join(lines)


def main() -> None:
    ledger = build_ledger()
    stats = summarize(ledger)
    payload = {"schema": SCHEMA_VERSION,
               "description": "Geometric q-root multiplicity ledger for f31 subcase (1)",
               "summary": stats, "strata": ledger}
    JSON_PATH.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    MD_PATH.write_text(render_markdown(ledger, stats), encoding="utf-8")
    print(json.dumps(stats, indent=2, sort_keys=True))
    print(f"wrote {JSON_PATH.name} and {MD_PATH.name}")


if __name__ == "__main__":
    main()
