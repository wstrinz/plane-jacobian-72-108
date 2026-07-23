#!/usr/bin/env python3
"""Independent exhaustive flag-case audit of the coupled q+t cascade claims."""

from __future__ import annotations

import argparse
import itertools
import json
import re
import time
from dataclasses import dataclass
from pathlib import Path

import sympy as sp

HERE = Path(__file__).resolve().parent
INF = None
STANDARD_MAX_A = 10
_NO_VALUE = object()


@dataclass(frozen=True)
class Window:
    name: str
    q_file: str
    qt_file: str
    var_caps: tuple[int, int, int]
    e_cap: int
    h_intercept: int
    h_slope: int
    expected_branches: int
    expected_q_cases: int
    expected_qt_cases: int

    def h_cap(self, level: int) -> int:
        return self.h_intercept - self.h_slope * level


WINDOWS = (
    Window("sub2", "cascade_cones.json", "cascade_cones_qt.json",
           (4, 6, 8), 10, 40, 4, 420, 320, 232),
    Window("sub1", "cascade_cones_sub1_depth4.json", "cascade_cones_sub1_qt.json",
           (6, 9, 12), 15, 60, 6, 2178, 2519, 2253),
)
TERMINALS = {"T1": 7, "T2": 6}
SUB1_TERMINAL_CAPS = {"T1": 46, "T2": 48}


@dataclass(frozen=True)
class Flags:
    d2_zero: bool
    sigma_zero: bool
    g_zero: frozenset[int]

    def label(self) -> str:
        levels = ",".join(map(str, sorted(self.g_zero))) or "-"
        return (f"d2_zero={self.d2_zero},sigma_zero={self.sigma_zero},"
                f"g_zero={levels}")


@dataclass(frozen=True)
class HDomain:
    low: int | None
    high: int | None
    allow_infinity: bool


@dataclass(frozen=True)
class LocalProfile:
    resources: tuple[int, ...]
    vd2: int | None
    vd1: int | None
    vsigma: int | None
    gvals: tuple[tuple[int, int | None], ...]
    hvals: tuple[tuple[int, int | None], ...]


@dataclass
class WindowResult:
    window: str
    flag_total: int
    flag_agreements: int
    claimed_survivors: int
    audit_survivors: int
    branch_total: int
    branch_agreements: int
    disagreements: list[dict]
    consistency_errors: list[str]
    branch_audit_errors: list[str]
    structural_errors: list[str]
    elapsed: float


def load_h_monomials():
    """Parse h_0..h_7 independently, verify both windows, and eliminate d0."""
    text = (HERE / "f31_graded.txt").read_text(encoding="utf-8")
    d2, d1, d0, e, sigma = sp.symbols("d2 d1 d0 dm1 sigma")
    source_vars = (d2, d1, d0, e)
    weights = (2, 3, 4, 5)
    pattern = re.compile(
        r"^h_(\d+) \(weight (\d+), dm1-power (\d+)\) = (.+)$", re.MULTILINE
    )
    parsed, metadata = {}, {}
    for match in pattern.finditer(text):
        level, weight, outer_power = map(int, match.group(1, 2, 3))
        parsed[level] = sp.sympify(
            match.group(4),
            locals={"d2": d2, "d1": d1, "d0": d0, "dm1": e},
        )
        metadata[level] = (weight, outer_power)
    if sorted(parsed) != list(range(8)):
        raise AssertionError(f"expected h_0..h_7, found {sorted(parsed)}")

    result = {}
    for level in range(8):
        expected_weight = 20 - 2 * level
        if metadata[level] != (expected_weight, 21 - 3 * level):
            raise AssertionError(f"h_{level}: bad graded metadata {metadata[level]}")
        source = sp.Poly(parsed[level], *source_vars)
        for exponents, _coefficient in source.terms():
            actual_weight = sum(w * n for w, n in zip(weights, exponents))
            if actual_weight != expected_weight:
                raise AssertionError(
                    f"h_{level}: exponent {exponents} has weight {actual_weight}, "
                    f"expected {expected_weight}"
                )
        rewritten = sp.Poly(
            sp.expand(parsed[level].subs(d0, (sigma + d2**2) / 4)),
            d2, d1, sigma, e,
        )
        terms = tuple(
            (coefficient, tuple(map(int, exponents)))
            for exponents, coefficient in rewritten.terms() if coefficient
        )
        for window in WINDOWS:
            caps = (*window.var_caps, window.e_cap)
            maximum = max(
                sum(cap * exponent for cap, exponent in zip(caps, exponents))
                for _coefficient, exponents in terms
            )
            if maximum > window.h_cap(level):
                raise AssertionError(
                    f"{window.name} h_{level}: rewritten degree {maximum} exceeds "
                    f"{window.h_cap(level)}"
                )
        result[level] = terms
    return result


def g_caps(window: Window, branch: str, a_t: int) -> dict[int, int]:
    """Derive the stated safe per-level g caps without engine data."""
    if not 0 <= a_t <= STANDARD_MAX_A:
        raise ValueError("this audit is restricted to the standard regime a<=10")
    terminal = TERMINALS[branch]
    if window.name == "sub2":
        return {level: 10 + 3 * a_t for level in range(1, terminal + 1)}
    v = 30 - 3 * a_t
    ehat_cap = window.e_cap - a_t
    forward = {1: 60 - v}
    for level in range(1, terminal):
        forward[level + 1] = (
            max(3 * ehat_cap + forward[level], 60 - 2 * level) - v
        )
    backward = {terminal: SUB1_TERMINAL_CAPS[branch]}
    for level in range(terminal - 1, 0, -1):
        backward[level] = max(v + backward[level + 1], 60 - 2 * level)
    return {level: min(forward[level], backward[level])
            for level in range(1, terminal + 1)}


def h_domain(monomials, values, degree_cap: int) -> HDomain | None:
    """Conservative h valuation domain: unique minima force; ties may rise."""
    valuations = []
    for _coefficient, exponents in monomials:
        total = 0
        for exponent, value in zip(exponents, values):
            if exponent and value is INF:
                break
            if exponent:
                total += exponent * int(value)
        else:
            valuations.append(total)
    if not valuations:
        return HDomain(None, None, True)
    minimum = min(valuations)
    if valuations.count(minimum) == 1:
        return HDomain(minimum, minimum, False) if minimum <= degree_cap else None
    if minimum <= degree_cap:
        return HDomain(minimum, degree_cap, True)
    return HDomain(None, None, True)


def equation_holds(result, left, right) -> bool:
    """Ultrametric feasibility for result = left + right."""
    if result is INF:
        return left == right
    if left is INF:
        return right == result
    if right is INF:
        return left == result
    return result >= left if left == right else result == min(left, right)


def choose_h(domain: HDomain, result, other, shift: int):
    """Return one h valuation in a domain satisfying an ultrametric equation."""
    candidates = set()
    if domain.low is not None:
        low, high = shift + domain.low, shift + int(domain.high)
        probes = [low, high]
        for value in (result, other):
            if value is not INF:
                probes.extend((value, value + 1))
        candidates.update(value for value in probes if low <= value <= high)
    ordered = sorted(candidates)
    if domain.allow_infinity:
        ordered.append(INF)
    for shifted_h in ordered:
        if equation_holds(result, other, shifted_h):
            return INF if shifted_h is INF else shifted_h - shift
    return _NO_VALUE


def resource_vector(vd2, vd1, vsigma, gvals, levels):
    return (
        0 if vd2 is INF else vd2,
        0 if vd1 is INF else vd1,
        0 if vsigma is INF else vsigma,
        *(0 if gvals[level] is INF else gvals[level] for level in levels),
    )


def pareto_reduce(profiles):
    """Retain one witness per componentwise-minimal shared-budget vector."""
    by_resources = {}
    for profile in profiles:
        by_resources.setdefault(profile.resources, profile)
    frontier = []
    for resources in sorted(by_resources, key=lambda item: (sum(item), item)):
        if any(all(old <= new for old, new in zip(saved.resources, resources))
               for saved in frontier):
            continue
        frontier.append(by_resources[resources])
    return tuple(frontier)


def compact_states(states):
    """Pareto-compact partial descents having the same current g order."""
    groups = {}
    for gvals, hvals in states:
        current = gvals[-1][1]
        resources = tuple(0 if value is INF else value for _level, value in gvals)
        groups.setdefault(current, {}).setdefault(resources, (gvals, hvals))
    result = []
    for by_resources in groups.values():
        frontier = []
        for resources in sorted(by_resources, key=lambda item: (sum(item), item)):
            if any(all(old <= new for old, new in zip(saved, resources))
                   for saved in frontier):
                continue
            frontier.append(resources)
            result.append(by_resources[resources])
    return result


def variable_ranges(window: Window, branch: str, flags: Flags):
    d2_values = (INF,) if flags.d2_zero else range(window.var_caps[0] + 1)
    d1_values = range(window.var_caps[1] + 1) if branch == "T1" else (INF,)
    sigma_values = ((INF,) if flags.sigma_zero
                    else range(window.var_caps[2] + 1))
    return d2_values, d1_values, sigma_values


def local_profiles(window: Window, branch: str, a_t: int, place_parameter: int,
                   flags: Flags, monomial_data, minimum_level: int, place: str):
    """Exhaust one q or t local chain using only the stated place semantics."""
    terminal = TERMINALS[branch]
    caps = g_caps(window, branch, a_t)
    descent_levels = tuple(range(terminal - 1, minimum_level - 1, -1))
    all_levels = tuple(range(minimum_level, terminal + 1))
    profiles = []

    for vd2, vd1, vsigma in itertools.product(
            *variable_ranges(window, branch, flags)):
        auxiliary = vd1 if branch == "T1" else vsigma
        if auxiliary is INF:
            continue
        terminal_order = (terminal + 2 * auxiliary - 3 * place_parameter
                          if place == "q" else 2 * auxiliary)
        if not 0 <= terminal_order <= caps[terminal]:
            continue

        states = [(((terminal, terminal_order),), ())]
        h_values = (vd2, vd1, vsigma,
                    place_parameter if place == "q" else a_t)
        for level in descent_levels:
            domain = h_domain(monomial_data[level], h_values, window.h_cap(level))
            if domain is None:
                states = []
                break
            expanded = []
            candidates_g = ((INF,) if level in flags.g_zero
                            else range(caps[level] + 1))
            for gvals, hvals in states:
                next_order = gvals[-1][1]
                if place == "q":
                    result_order, shift = next_order, level
                else:
                    v = 30 - 3 * a_t
                    result_order = INF if next_order is INF else v + next_order
                    shift = 0
                for here_order in candidates_g:
                    if place == "q":
                        other = (INF if here_order is INF
                                 else 3 * place_parameter + here_order)
                    else:
                        other = here_order
                    h_value = choose_h(domain, result_order, other, shift)
                    if h_value is not _NO_VALUE:
                        expanded.append((gvals + ((level, here_order),),
                                         hvals + ((level, h_value),)))
            states = compact_states(expanded)
            if not states:
                break

        for gdesc, hdesc in states:
            gdict = dict(gdesc)
            profiles.append(LocalProfile(
                resource_vector(vd2, vd1, vsigma, gdict, all_levels),
                vd2, vd1, vsigma, tuple(sorted(gdesc)), tuple(sorted(hdesc)),
            ))
    return pareto_reduce(profiles)


def join_places(candidates, caps):
    """Join any number of places against the one shared global budget vector."""
    count = len(candidates)
    order = sorted(range(count), key=lambda index: len(candidates[index]))
    failed = set()

    def dfs(position, remaining):
        key = (position, remaining)
        if key in failed:
            return None
        if position == count:
            return {}
        place = order[position]
        for profile in candidates[place]:
            if all(used <= left for used, left in zip(profile.resources, remaining)):
                suffix = dfs(position + 1, tuple(
                    left - used for used, left in zip(profile.resources, remaining)
                ))
                if suffix is not None:
                    suffix[place] = profile
                    return suffix
        failed.add(key)
        return None

    found = dfs(0, caps)
    return None if found is None else [found[index] for index in range(count)]


def flag_cases(branch: str, minimum_level: int):
    """Enumerate all allowed global flags; the terminal g is never zero."""
    terminal = TERMINALS[branch]
    lower = tuple(range(minimum_level, terminal))
    sigma_choices = (False, True) if branch == "T1" else (False,)
    for d2_zero, sigma_zero in itertools.product((False, True), sigma_choices):
        for mask in range(1 << len(lower)):
            yield Flags(d2_zero, sigma_zero, frozenset(
                level for bit, level in enumerate(lower) if mask & (1 << bit)
            ))


def branch_key(row: dict):
    return row["a_t"], tuple(row["b"]), row["branch"]


def branch_label(key) -> str:
    a_t, b_values, branch = key
    return f"a={a_t};b={','.join(map(str, b_values))};{branch}"


def case_key(case: dict) -> Flags:
    return Flags(bool(case["d2_zero"]), bool(case["sigma_zero"]),
                 frozenset(map(int, case["g_zero_levels"])))


def status(row: dict) -> str:
    return "killed" if row["status"] == "engine_killed_pending_audit" else "survives"


def fmt_value(value) -> str:
    return "inf" if value is INF else str(value)


def witness_text(profiles: list[LocalProfile], b_values: tuple[int, ...]) -> str:
    chunks = []
    names = [f"q{i + 1}(b={b})" for i, b in enumerate(b_values)] + ["t"]
    for name, profile in zip(names, profiles):
        gs = ",".join(f"g{level}={fmt_value(value)}"
                      for level, value in profile.gvals)
        hs = ",".join(f"h{level}={fmt_value(value)}"
                      for level, value in profile.hvals)
        chunks.append(
            f"{name}(d2={fmt_value(profile.vd2)},d1={fmt_value(profile.vd1)},"
            f"sigma={fmt_value(profile.vsigma)};{gs};{hs})"
        )
    return " ".join(chunks)


def compare_branch_files(q_claims, qt_claims):
    """Explicitly compare audited q-only and q+t branch verdicts."""
    q_rows = {branch_key(row): row for row in q_claims["branches"]}
    qt_rows = {branch_key(row): row for row in qt_claims["branches"]}
    errors = []
    for key in sorted(set(q_rows) - set(qt_rows)):
        errors.append(f"{branch_label(key)}: missing from q+t file")
    for key in sorted(set(qt_rows) - set(q_rows)):
        errors.append(f"{branch_label(key)}: absent from q-only file")
    for key in sorted(set(q_rows) & set(qt_rows)):
        if status(q_rows[key]) != status(qt_rows[key]):
            errors.append(f"{branch_label(key)}: q-only={status(q_rows[key])}, "
                          f"q+t={status(qt_rows[key])}")
    return q_rows, qt_rows, errors


def audit_window(window: Window, monomial_data, quiet: bool) -> WindowResult:
    started = time.perf_counter()
    q_claims = json.loads((HERE / window.q_file).read_text(encoding="utf-8"))
    qt_claims = json.loads((HERE / window.qt_file).read_text(encoding="utf-8"))
    q_rows, qt_rows, consistency_errors = compare_branch_files(q_claims, qt_claims)
    structural_errors = []
    if q_claims.get("depth") != 4 or qt_claims.get("depth") != 4:
        structural_errors.append("both claim files must be complete depth-4 outputs")
    if qt_claims.get("places") != "q+t":
        structural_errors.append("q+t claim file does not declare places=q+t")
    if len(q_rows) != window.expected_branches or len(qt_rows) != window.expected_branches:
        structural_errors.append(
            f"branch count q={len(q_rows)}, q+t={len(qt_rows)}, "
            f"expected={window.expected_branches}"
        )
    q_case_count = sum(len(row.get("survivor_cases", [])) for row in q_rows.values())
    qt_case_count = sum(len(row.get("survivor_cases", [])) for row in qt_rows.values())
    if q_case_count != window.expected_q_cases:
        structural_errors.append(
            f"q-only survivor case count {q_case_count}, expected {window.expected_q_cases}"
        )
    if qt_case_count != window.expected_qt_cases:
        structural_errors.append(
            f"q+t survivor case count {qt_case_count}, expected {window.expected_qt_cases}"
        )

    q_cache, t_cache = {}, {}
    flag_total = flag_agreements = audit_survivors = branch_agreements = 0
    disagreements, branch_audit_errors = [], []
    minimum_level = 4

    for number, key in enumerate(sorted(qt_rows), 1):
        row = qt_rows[key]
        a_t, b_values, branch = key
        legal_flags = tuple(flag_cases(branch, minimum_level))
        legal_set = set(legal_flags)
        q_list = q_rows[key].get("survivor_cases", [])
        qt_list = row.get("survivor_cases", [])
        q_expected = [case_key(case) for case in q_list]
        expected = [case_key(case) for case in qt_list]
        if len(set(q_expected)) != len(q_expected):
            structural_errors.append(f"{branch_label(key)}: duplicate q-only flag cases")
        if len(set(expected)) != len(expected):
            structural_errors.append(f"{branch_label(key)}: duplicate q+t flag cases")
        for flags in set(q_expected) - legal_set:
            structural_errors.append(
                f"{branch_label(key)}: illegal q-only flag case {flags.label()}"
            )
        for flags in set(expected) - legal_set:
            structural_errors.append(
                f"{branch_label(key)}: illegal q+t flag case {flags.label()}"
            )
        for flags in set(expected) - set(q_expected):
            structural_errors.append(
                f"{branch_label(key)}: q+t case not in q-only cases: {flags.label()}"
            )
        if row.get("survivor_case_count") != len(qt_list):
            structural_errors.append(
                f"{branch_label(key)}: survivor_case_count={row.get('survivor_case_count')} "
                f"but list has {len(qt_list)}"
            )

        expected_set = set(expected)
        caps_by_level = g_caps(window, branch, a_t)
        levels = tuple(range(minimum_level, TERMINALS[branch] + 1))
        budgets = (*window.var_caps, *(caps_by_level[level] for level in levels))
        branch_has_feasible = False


        for flags in legal_flags:
            flag_total += 1
            q_candidates = []
            for b_value in b_values:
                cache_key = (branch, a_t, b_value, flags, minimum_level)
                if cache_key not in q_cache:
                    q_cache[cache_key] = local_profiles(
                        window, branch, a_t, b_value, flags,
                        monomial_data, minimum_level, "q",
                    )
                if not q_cache[cache_key]:
                    q_candidates = []
                    break
                q_candidates.append(q_cache[cache_key])

            q_found = (join_places(q_candidates, budgets)
                       if len(q_candidates) == 4 else None)
            found = None
            if q_found is not None:
                t_key = (branch, a_t, flags, minimum_level)
                if t_key not in t_cache:
                    t_cache[t_key] = local_profiles(
                        window, branch, a_t, a_t, flags,
                        monomial_data, minimum_level, "t",
                    )
                if t_cache[t_key]:
                    found = join_places(q_candidates + [t_cache[t_key]], budgets)

            actual = found is not None
            claimed = flags in expected_set
            audit_survivors += actual
            flag_agreements += actual == claimed
            branch_has_feasible |= actual
            if actual != claimed:
                disagreements.append({
                    "branch": branch_label(key),
                    "flags": flags.label(),
                    "expected": "survives" if claimed else "killed",
                    "actual": "survives" if actual else "killed",
                    "q_feasible": q_found is not None,
                    "witness": witness_text(found, b_values) if found is not None else None,
                })

        audit_branch_status = "survives" if branch_has_feasible else "killed"
        claimed_branch_status = status(row)
        if audit_branch_status == claimed_branch_status:
            branch_agreements += 1
        else:
            branch_audit_errors.append(
                f"{branch_label(key)}: q+t file={claimed_branch_status}, "
                f"flag audit={audit_branch_status}"
            )
        if not quiet:
            verdict = "AGREE" if audit_branch_status == claimed_branch_status else "DISAGREE"
            print(f"{window.name} {number:04d} {branch_label(key):<38} "
                  f"claim={claimed_branch_status:<8} audit={audit_branch_status:<8} {verdict}")

    elapsed = time.perf_counter() - started
    return WindowResult(
        window.name, flag_total, flag_agreements, qt_case_count, audit_survivors,
        len(qt_rows), branch_agreements, disagreements, consistency_errors,
        branch_audit_errors, structural_errors, elapsed,
    )


def print_result(result: WindowResult):
    print(f"{result.window}_flag_cases: {result.flag_agreements}/{result.flag_total} agree; "
          f"disagreements={len(result.disagreements)}; "
          f"claimed_survivors={result.claimed_survivors}; "
          f"audit_survivors={result.audit_survivors}")
    print(f"{result.window}_branch_audit: {result.branch_agreements}/{result.branch_total} "
          f"agree with q+t branch status; "
          f"disagreements={len(result.branch_audit_errors)}")
    print(f"{result.window}_q_vs_qt_branch_consistency: "
          f"{'PASS' if not result.consistency_errors else 'FAIL'}; "
          f"mismatches={len(result.consistency_errors)}")
    for item in result.structural_errors:
        print(f"STRUCTURAL_ERROR {result.window} {item}")
    for item in result.consistency_errors:
        print(f"BRANCH_CONSISTENCY_ERROR {result.window} {item}")
    for item in result.branch_audit_errors:
        print(f"BRANCH_AUDIT_ERROR {result.window} {item}")
    for item in result.disagreements:
        print(f"FLAG_DISAGREEMENT {result.window} {item['branch']} {item['flags']} "
              f"expected={item['expected']} actual={item['actual']} "
              f"q_feasible={item['q_feasible']}")
        if item["witness"]:
            print("  CONSERVATIVE_WITNESS", item["witness"])
    print(f"{result.window}_runtime_seconds: {result.elapsed:.3f}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quiet", action="store_true",
                        help="omit the per-branch verdict tables")
    parser.add_argument(
        "--consistency-only", action="store_true",
        help="only compare q-only and q+t branch keys and verdicts",
    )
    args = parser.parse_args()
    started = time.perf_counter()

    if args.consistency_only:
        all_errors = []
        for window in WINDOWS:
            q_claims = json.loads((HERE / window.q_file).read_text(encoding="utf-8"))
            qt_claims = json.loads((HERE / window.qt_file).read_text(encoding="utf-8"))
            q_rows, qt_rows, errors = compare_branch_files(q_claims, qt_claims)
            all_errors.extend(f"{window.name} {item}" for item in errors)
            common = set(q_rows) & set(qt_rows)
            agreements = sum(status(q_rows[key]) == status(qt_rows[key])
                             for key in common)
            print(f"{window.name}_q_vs_qt_branch_consistency: "
                  f"{'PASS' if not errors else 'FAIL'}; "
                  f"agreements={agreements}/"
                  f"{len(set(q_rows) | set(qt_rows))}; mismatches={len(errors)}")
        for item in all_errors:
            print("BRANCH_CONSISTENCY_ERROR", item)
        print("branch_consistency_overall: "
              f"{'PASS' if not all_errors else 'FAIL'}; mismatches={len(all_errors)}")
        return 1 if all_errors else 0

    monomial_data = load_h_monomials()
    print("weighted_homogeneity: PASS (h_0..h_7 parsed independently)")
    print("degree_caps: PASS (both window rewrites checked)")
    print("d0_elimination: PASS (substituted d0=(sigma+d2^2)/4)")
    for window in WINDOWS:
        for branch in TERMINALS:
            for a_t in range(STANDARD_MAX_A + 1):
                g_caps(window, branch, a_t)
    print("g_cap_derivation: PASS (both windows, T1/T2, a=0..10)")
    print("terminal_g_zero_policy: PASS (terminal excluded from every flag mask)")

    results = []
    for window in WINDOWS:
        result = audit_window(window, monomial_data, args.quiet)
        results.append(result)
        print_result(result)

    consistency_errors = [item for result in results
                          for item in result.consistency_errors]
    print("branch_consistency_overall: "
          f"{'PASS' if not consistency_errors else 'FAIL'}; "
          f"mismatches={len(consistency_errors)}")
    elapsed = time.perf_counter() - started
    print(f"runtime_seconds: {elapsed:.3f}")
    failed = any(result.disagreements or result.consistency_errors
                 or result.branch_audit_errors or result.structural_errors
                 for result in results)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
