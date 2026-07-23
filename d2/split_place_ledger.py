#!/usr/bin/env python3
"""Generate the geometric q-place ledger for the f31 subcase-(2) campaign.

After base change, q has four simple roots p_i. A stratum is recorded by

    a = v_t(e),  b_i = v_{p_i}(e/t^a),  b_1 >= ... >= b_4,

with a + sum(b_i) <= 10. The script performs the field-stable integer pruning
from the last two cascade levels and records which uniform strata are already
covered by repository proofs.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from itertools import product
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
JSON_PATH = ROOT / "split_place_ledger.json"
MD_PATH = ROOT / "SPLIT_PLACE_LEDGER.md"


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
    """Partitions of q-root multiplicity into at most four places."""

    return [
        b
        for b in product(range(total_cap + 1), repeat=4)
        if sum(b) <= total_cap and tuple(sorted(b, reverse=True)) == b
    ]


def minimum_order(rhs_q_order: int, b: int) -> int:
    """Smallest x >= 0 for which rhs_q_order + 2*x - 3*b >= 0."""

    return max(0, (3 * b - rhs_q_order + 1) // 2)


def terminal_test(
    a: int, b: tuple[int, int, int, int], *, branch: str
) -> TerminalTest:
    """Apply exact necessary valuation and degree conditions.

    T1 (d1 != 0), level 7:
        3*b_i + v_i(g7) = 7 + 2*v_i(d1),  deg d1 <= 6.

    T2 (d1 = 0, sigma != 0), level 6:
        3*b_i + v_i(g6) = 6 + 2*v_i(sigma),  deg sigma <= 8.

    The sums of the local orders lower-bound the global polynomial degrees.
    Failure is therefore a rigorous geometric stratum kill.
    """

    if branch == "T1":
        level, auxiliary, auxiliary_cap = 7, "d1", 6
    elif branch == "T2":
        level, auxiliary, auxiliary_cap = 6, "sigma", 8
    else:
        raise ValueError(f"unknown terminal branch: {branch}")

    g_cap = 10 + 3 * a
    aux_orders = tuple(minimum_order(level, bi) for bi in b)
    g_orders = tuple(level + 2 * xi - 3 * bi for xi, bi in zip(aux_orders, b))
    assert all(order >= 0 for order in g_orders)
    aux_degree = sum(aux_orders)
    g_degree = sum(g_orders)
    return TerminalTest(
        branch=branch,
        source_level=level,
        auxiliary=auxiliary,
        auxiliary_degree_cap=auxiliary_cap,
        g_degree_cap=g_cap,
        minimum_auxiliary_orders=aux_orders,
        resulting_g_orders=g_orders,
        minimum_auxiliary_degree=aux_degree,
        minimum_g_degree=g_degree,
        feasible=aux_degree <= auxiliary_cap and g_degree <= g_cap,
    )


def existing_proof(
    a: int, b: tuple[int, int, int, int], branch: str
) -> dict[str, str] | None:
    """Return an old proof only when its geometric hypotheses are exact."""

    if b == (0, 0, 0, 0):
        if 0 <= a <= 4:
            return {
                "status": "proven_infeasible",
                "reference": "T5_MULTIPLACE.md, Theorem 1 (geometric q-coprime reading)",
            }
        if a == 5:
            return {
                "status": "proven_infeasible",
                "reference": "T5_STRATA_50_11.md, Theorem 3",
            }
        if a == 6:
            return {
                "status": "proven_infeasible",
                "reference": f"T5_60_{branch}.md",
            }
        if a == 7:
            return {
                "status": "proven_infeasible",
                "reference": "FIELD_SPLIT_AUDIT.md, geometric q-coprime a=7 theorem",
            }
        if a == 9 and branch == "T2":
            return {
                "status": "proven_infeasible",
                "reference": "T5_90_T2.md",
            }
        if a == 10 and branch == "T2":
            return {
                "status": "proven_infeasible",
                "reference": "T5_STRATUM_10_0.md, d1=0 branch",
            }

    if b == (1, 1, 1, 1):
        if a == 0:
            return {
                "status": "proven_infeasible",
                "reference": "T5_MULTIPLACE.md, Theorem 1 (uniform q case)",
            }
        if a == 1:
            return {
                "status": "proven_infeasible",
                "reference": "T5_STRATA_50_11.md, Theorem 4",
            }
        if 2 <= a <= 5:
            return {
                "status": "proven_infeasible",
                "reference": "T5_T1_AQ12.md",
            }
        if a == 6 and branch == "T2":
            return {
                "status": "proven_infeasible",
                "reference": "T5_T1_AQ12.md, (6,1) T2",
            }

    if b == (2, 2, 2, 2) and 0 <= a <= 2:
        return {
            "status": "proven_infeasible",
            "reference": "T5_T1_AQ12.md, a_q=2 branches",
        }

    return None


def branch_record(
    a: int, b: tuple[int, int, int, int], branch: str
) -> dict[str, Any]:
    terminal = terminal_test(a, b, branch=branch)
    proof = existing_proof(a, b, branch)
    if not terminal.feasible:
        status = "proven_infeasible_terminal"
        reference = (
            f"split-place level {terminal.source_level}: "
            f"min deg {terminal.auxiliary}={terminal.minimum_auxiliary_degree} "
            f"(cap {terminal.auxiliary_degree_cap}), "
            f"min deg g{terminal.source_level}={terminal.minimum_g_degree} "
            f"(cap {terminal.g_degree_cap})"
        )
    elif proof:
        status = proof["status"]
        reference = proof["reference"]
    else:
        status = "open_after_terminal"
        reference = "requires lower cascade levels or another exact argument"
    return {
        "status": status,
        "reference": reference,
        "terminal": asdict(terminal),
    }


def build_ledger() -> list[dict[str, Any]]:
    ledger: list[dict[str, Any]] = []
    for a in range(11):
        for b in sorted_q_vectors(10 - a):
            uniform = len(set(b)) == 1
            support_type = (
                "geometrically_q_coprime"
                if b == (0, 0, 0, 0)
                else "uniform_q_power"
                if uniform
                else "partial_q_support"
            )
            branches = {
                branch: branch_record(a, b, branch) for branch in ("T1", "T2")
            }
            open_branches = [
                branch
                for branch, data in branches.items()
                if data["status"] == "open_after_terminal"
            ]
            ledger.append(
                {
                    "a_t": a,
                    "b": list(b),
                    "q_multiplicity_sum": sum(b),
                    "residual_degree_budget": 10 - a - sum(b),
                    "support_type": support_type,
                    "T3": {
                        "status": "proven_infeasible",
                        "reference": "FIELD_SPLIT_AUDIT.md, split-place sigma theorem",
                    },
                    "branches": branches,
                    "open_branches": open_branches,
                    "stratum_status": (
                        "open_after_terminal" if open_branches else "proven_infeasible"
                    ),
                }
            )
    assert len(ledger) == 327
    return ledger


def summarize(ledger: list[dict[str, Any]]) -> dict[str, int]:
    branches = [data for row in ledger for data in row["branches"].values()]
    terminal_killed = [
        row
        for row in ledger
        if all(
            data["status"] == "proven_infeasible_terminal"
            for data in row["branches"].values()
        )
    ]
    open_rows = [row for row in ledger if row["open_branches"]]
    return {
        "raw_strata": len(ledger),
        "old_uniform_strata": sum(len(set(row["b"])) == 1 for row in ledger),
        "partial_support_strata": sum(len(set(row["b"])) > 1 for row in ledger),
        "T1_terminal_feasible": sum(
            row["branches"]["T1"]["terminal"]["feasible"] for row in ledger
        ),
        "T2_terminal_feasible": sum(
            row["branches"]["T2"]["terminal"]["feasible"] for row in ledger
        ),
        "strata_killed_by_terminal_both_branches": len(terminal_killed),
        "partial_support_strata_killed_by_terminal_both_branches": sum(
            row["support_type"] == "partial_q_support" for row in terminal_killed
        ),
        "open_strata_after_terminal_and_existing_proofs": len(open_rows),
        "open_partial_support_strata": sum(
            row["support_type"] == "partial_q_support" for row in open_rows
        ),
        "open_branches_after_terminal_and_existing_proofs": sum(
            data["status"] == "open_after_terminal" for data in branches
        ),
    }


def render_markdown(ledger: list[dict[str, Any]], stats: dict[str, int]) -> str:
    by_a = []
    for a in range(11):
        rows = [row for row in ledger if row["a_t"] == a]
        by_a.append((a, len(rows), sum(bool(row["open_branches"]) for row in rows)))
    uniform_open = [
        row for row in ledger if len(set(row["b"])) == 1 and row["open_branches"]
    ]
    lines = [
        "# Geometric split-place ledger — f31 subcase (2)",
        "",
        "**Generated by §split_place_ledger.py§; do not hand-edit counts.**",
        "",
        "After base change write §q=p1*p2*p3*p4§, §a=v_t(e)§, and",
        "§b_i=v_{p_i}(e/t^a)§ in decreasing order. The cap is §a+sum(b_i)<=10§.",
        "",
        "## Exact terminal pruning",
        "",
        "T1: §3b_i+v_i(g7)=7+2v_i(d1)§, with §deg d1<=6§.",
        "T2: §3b_i+v_i(g6)=6+2v_i(sigma)§, with §deg sigma<=8§.",
        "In both cases §deg g_l<=10+3a§. Summing minimum local orders gives",
        "rigorous necessary degree inequalities.",
        "",
        "## Counts",
        "",
        f"- Raw geometric strata: **{stats['raw_strata']}**.",
        f"- Old uniform strata: **{stats['old_uniform_strata']}**; partial-support strata: **{stats['partial_support_strata']}**.",
        f"- T1 terminal-feasible: **{stats['T1_terminal_feasible']}**; T2 terminal-feasible: **{stats['T2_terminal_feasible']}**.",
        f"- Killed at terminal levels in both live branches: **{stats['strata_killed_by_terminal_both_branches']}**, including **{stats['partial_support_strata_killed_by_terminal_both_branches']}** partial-support strata.",
        f"- Open after terminal pruning plus correctly scoped existing proofs: **{stats['open_strata_after_terminal_and_existing_proofs']}** strata / **{stats['open_branches_after_terminal_and_existing_proofs']}** branches.",
        "",
        "T3 (§d1=0§, §sigma=0§) is excluded globally by the field-stable",
        "split-place sigma-locus theorem.",
        "",
        "## By t-multiplicity",
        "",
        "| §a§ | raw | open after exact kills |",
        "|---:|---:|---:|",
    ]
    lines.extend(f"| {a} | {raw} | {opened} |" for a, raw, opened in by_a)
    lines.extend(
        [
            "",
            "## Uniform frontier",
            "",
            "| §a§ | §(b1,b2,b3,b4)§ | open branches |",
            "|---:|:---:|:---|",
        ]
    )
    lines.extend(
        f"| {row['a_t']} | §{tuple(row['b'])}§ | {', '.join(row['open_branches'])} |"
        for row in uniform_open
    )
    lines.extend(
        [
            "",
            "The complete valuation witnesses and proof references are in",
            "§split_place_ledger.json§.",
            "",
        ]
    )
    return "\n".join(lines).replace("§", chr(96))


def main() -> None:
    ledger = build_ledger()
    stats = summarize(ledger)
    payload = {
        "schema": 1,
        "description": "Geometric q-root multiplicity ledger for f31 subcase (2)",
        "summary": stats,
        "strata": ledger,
    }
    JSON_PATH.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    MD_PATH.write_text(render_markdown(ledger, stats), encoding="utf-8")
    print(json.dumps(stats, indent=2, sort_keys=True))
    print(f"wrote {JSON_PATH.name} and {MD_PATH.name}")


if __name__ == "__main__":
    main()
