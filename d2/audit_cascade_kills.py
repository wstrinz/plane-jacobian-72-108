#!/usr/bin/env python3
"""Independent exhaustive valuation audit for the q-place cascade claims."""

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
VAR_CAPS = (4, 6, 8)  # d2, d1, sigma


@dataclass(frozen=True)
class Flags:
    d2_zero: bool
    sigma_zero: bool
    g_zero: frozenset[int]

    def label(self) -> str:
        gs = ",".join(map(str, sorted(self.g_zero))) or "-"
        return f"d2_zero={self.d2_zero},sigma_zero={self.sigma_zero},g_zero={gs}"


@dataclass(frozen=True)
class LocalProfile:
    resources: tuple[int, ...]
    vd2: int | None
    vd1: int | None
    vsigma: int | None
    gvals: tuple[tuple[int, int | None], ...]
    hvals: tuple[tuple[int, int | None], ...]


def load_h_monomials():
    """Parse h_0..h_7, verify weights, then substitute d0=(sigma+d2^2)/4."""
    text = (HERE / "f31_graded.txt").read_text(encoding="utf-8")
    d2, d1, d0, e, sigma = sp.symbols("d2 d1 d0 dm1 sigma")
    old_vars = (d2, d1, d0, e)
    weights = (2, 3, 4, 5)
    pat = re.compile(
        r"^h_(\d+) \(weight (\d+), dm1-power (\d+)\) = (.+)$", re.MULTILINE
    )
    parsed, metadata = {}, {}
    for match in pat.finditer(text):
        level, weight, outer_power = map(int, match.group(1, 2, 3))
        parsed[level] = sp.sympify(match.group(4), locals={
            "d2": d2, "d1": d1, "d0": d0, "dm1": e
        })
        metadata[level] = (weight, outer_power)
    if sorted(parsed) != list(range(8)):
        raise AssertionError(f"expected h_0..h_7, found {sorted(parsed)}")

    result = {}
    for level in range(8):
        expected_weight = 20 - 2 * level
        if metadata[level][0] != expected_weight:
            raise AssertionError(f"h_{level}: header weight is wrong")
        if metadata[level][1] != 21 - 3 * level:
            raise AssertionError(f"h_{level}: header dm1-power is wrong")
        poly = sp.Poly(parsed[level], *old_vars)
        if not poly.terms():
            raise AssertionError(f"h_{level}: unexpectedly zero")
        for exponents, _coefficient in poly.terms():
            actual = sum(w * exponent for w, exponent in zip(weights, exponents))
            if actual != expected_weight:
                raise AssertionError(
                    f"h_{level}: exponent {exponents} has weight {actual}, "
                    f"expected {expected_weight}"
                )
        rewritten = sp.Poly(
            sp.expand(parsed[level].subs(d0, (sigma + d2**2) / 4)),
            d2, d1, sigma, e,
        )
        result[level] = tuple(
            (coefficient, tuple(map(int, exponents)))
            for exponents, coefficient in rewritten.terms() if coefficient
        )
    return result


def h_options(monomials, values, degree_cap: int) -> tuple[int | None, ...]:
    """All h-valuations under the deliberately permissive tie semantics."""
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
    if not valuations:  # all monomials were killed by global zero flags
        return (INF,)
    minimum = min(valuations)
    if valuations.count(minimum) == 1:
        return (minimum,) if minimum <= degree_cap else ()
    finite = tuple(range(minimum, degree_cap + 1)) if minimum <= degree_cap else ()
    return finite + (INF,)


def equation_holds(r_next, a_val, b_val) -> bool:
    """Ultrametric feasibility for g_next = A + B at one place."""
    if r_next is INF:
        return a_val == b_val  # exact cancellation, including 0 + 0
    if a_val is INF:
        return b_val == r_next
    if b_val is INF:
        return a_val == r_next
    return r_next >= a_val if a_val == b_val else r_next == min(a_val, b_val)


def resource_vector(vd2, vd1, vsigma, gvals, levels):
    return (
        0 if vd2 is INF else vd2,
        0 if vd1 is INF else vd1,
        0 if vsigma is INF else vsigma,
        *(0 if gvals[level] is INF else gvals[level] for level in levels),
    )


def pareto_reduce(profiles):
    """Retain one witness for every componentwise-minimal resource vector."""
    by_resources = {}
    for profile in profiles:
        by_resources.setdefault(profile.resources, profile)
    frontier = []
    for resources in sorted(by_resources, key=lambda x: (sum(x), x)):
        if any(all(a <= b for a, b in zip(old.resources, resources)) for old in frontier):
            continue
        frontier.append(by_resources[resources])
    return frontier


def local_profiles(branch: str, a_t: int, b: int, flags: Flags, monomial_data):
    """Exhaust every local integer valuation chain, then Pareto-reduce it."""
    terminal = 7 if branch == "T1" else 6
    descent_levels = tuple(range(terminal - 1, 3, -1))
    all_levels = tuple(range(4, terminal + 1))
    g_cap = 10 + 3 * a_t
    d2_values = (INF,) if flags.d2_zero else range(5)
    d1_values = range(7) if branch == "T1" else (INF,)
    sigma_values = (INF,) if flags.sigma_zero else range(9)
    profiles = []

    for vd2, vd1, vsigma in itertools.product(d2_values, d1_values, sigma_values):
        auxiliary = vd1 if branch == "T1" else vsigma
        if auxiliary is INF:
            continue
        r_terminal = terminal + 2 * auxiliary - 3 * b
        if not 0 <= r_terminal <= g_cap:
            continue

        def descend(index, gvals, hvals):
            if index == len(descent_levels):
                profiles.append(LocalProfile(
                    resource_vector(vd2, vd1, vsigma, gvals, all_levels),
                    vd2, vd1, vsigma,
                    tuple(sorted(gvals.items())), tuple(sorted(hvals.items())),
                ))
                return
            level = descent_levels[index]
            r_next = gvals[level + 1]
            candidates_r = (INF,) if level in flags.g_zero else range(g_cap + 1)
            hs = h_options(monomial_data[level], (vd2, vd1, vsigma, b), 40 - 4 * level)
            for r_here in candidates_r:
                a_val = INF if r_here is INF else 3 * b + r_here
                for vh in hs:
                    b_val = INF if vh is INF else level + vh
                    if equation_holds(r_next, a_val, b_val):
                        next_g, next_h = dict(gvals), dict(hvals)
                        next_g[level], next_h[level] = r_here, vh
                        descend(index + 1, next_g, next_h)

        descend(0, {terminal: r_terminal}, {})
    return pareto_reduce(profiles)


def join_places(candidates, caps):
    """Complete four-place DFS with memoized failed remaining-budget states."""
    order = sorted(range(4), key=lambda i: len(candidates[i]))
    failed = set()

    def dfs(position, remaining):
        key = (position, remaining)
        if key in failed:
            return None
        if position == 4:
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
    return None if found is None else [found[i] for i in range(4)]


def flag_cases(branch: str):
    terminal = 7 if branch == "T1" else 6
    lower = tuple(range(4, terminal))
    result = []
    sigma_choices = (False, True) if branch == "T1" else (False,)
    for d2_zero, sigma_zero in itertools.product((False, True), sigma_choices):
        for mask in range(1 << len(lower)):
            zeros = frozenset(level for bit, level in enumerate(lower) if mask & (1 << bit))
            result.append(Flags(d2_zero, sigma_zero, zeros))
    return result


def terminal_feasible(branch: str, a_t: int, b_values: tuple[int, ...]):
    """Exhaust terminal auxiliary orders without ledger-derived minima."""
    terminal = 7 if branch == "T1" else 6
    auxiliary_cap = 6 if branch == "T1" else 8
    g_cap = 10 + 3 * a_t
    for orders in itertools.product(range(auxiliary_cap + 1), repeat=4):
        if sum(orders) > auxiliary_cap:
            continue
        gorders = [terminal + 2 * order - 3 * b for order, b in zip(orders, b_values)]
        if min(gorders) >= 0 and sum(gorders) <= g_cap:
            return True, list(orders)
    return False, None


def fmt_value(value) -> str:
    return "inf" if value is INF else str(value)


def branch_id(row: dict) -> str:
    return f"a={row['a_t']};b={','.join(map(str, row['b']))};{row['branch']}"


def witness_text(flags: Flags, profiles: list[LocalProfile]) -> str:
    chunks = []
    for i, profile in enumerate(profiles, 1):
        gs = ",".join(f"g{level}={fmt_value(value)}" for level, value in profile.gvals)
        hs = ",".join(f"h{level}={fmt_value(value)}" for level, value in profile.hvals)
        chunks.append(
            f"p{i}(d2={fmt_value(profile.vd2)},d1={fmt_value(profile.vd1)},"
            f"sigma={fmt_value(profile.vsigma)};{gs};{hs})"
        )
    return flags.label() + " | " + " ".join(chunks)


def emit_artifact(path: str, window: str, results: list, summary: dict) -> None:
    """Write the per-branch audit verdicts so the coverage proof-DAG can machine-
    join this independent audit (promoting confirmed engine-kills to
    'independently-audited').  Deterministic (sorted keys, no timestamps)."""
    import hashlib
    src = Path(__file__).read_bytes()
    id_re = re.compile(r"a=(\d+);b=([\d,]+);(T\d)")
    recs = []
    for r in results:
        m = id_re.match(r["id"])
        recs.append({
            "a_t": int(m.group(1)),
            "b": [int(x) for x in m.group(2).split(",")],
            "branch": m.group(3),
            "claim": r["claim"], "audit": r["audit"],
            "agreement": r["agreement"],
        })
    recs.sort(key=lambda r: (r["a_t"], r["b"], r["branch"]))
    out = {
        "schema": 1,
        "generator": Path(__file__).name,
        "generator_sha256": hashlib.sha256(src).hexdigest(),
        "window": window,
        "summary": summary,
        "branches": recs,
    }
    target = path if Path(path).is_absolute() else str(HERE / path)
    with open(target, "w") as f:
        json.dump(out, f, indent=1, sort_keys=True)
    print(f"emitted audit artifact: {target} ({len(recs)} branches)")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quiet", action="store_true", help="omit the 420-row verdict table")
    parser.add_argument("--emit-artifact", nargs="?", const="audit_cascade_kills.json",
                        default=None, metavar="PATH",
                        help="write per-branch audit verdicts as JSON (for the proof-DAG join)")
    args = parser.parse_args()
    started = time.perf_counter()
    monomial_data = load_h_monomials()
    print("weighted_homogeneity: PASS (h_0..h_7 parsed independently)")
    print("d0_elimination: PASS (substituted d0=(sigma+d2^2)/4)")

    ledger = json.loads((HERE / "split_place_ledger.json").read_text(encoding="utf-8"))
    claims = json.loads((HERE / "cascade_cones.json").read_text(encoding="utf-8"))
    terminal_total = terminal_agree = 0
    terminal_disagreements = []
    for stratum in ledger["strata"]:
        b_values = tuple(stratum["b"])
        for branch in ("T1", "T2"):
            actual, _orders = terminal_feasible(branch, stratum["a_t"], b_values)
            claimed = bool(stratum["branches"][branch]["terminal"]["feasible"])
            terminal_total += 1
            terminal_agree += actual == claimed
            if actual != claimed:
                terminal_disagreements.append(
                    f"a={stratum['a_t']};b={','.join(map(str, b_values))};{branch}: "
                    f"ledger={claimed},audit={actual}"
                )
    print(f"terminal_crosscheck: {terminal_agree}/{terminal_total} agree; "
          f"disagreements={len(terminal_disagreements)}")
    if terminal_disagreements:
        for item in terminal_disagreements:
            print("TERMINAL_DISAGREEMENT", item)
        print("Refusing cascade results because the terminal sanity check failed.")
        return 2

    expected_open = {
        (s["a_t"], tuple(s["b"]), branch)
        for s in ledger["strata"] for branch in s["open_branches"]
    }
    claim_keys = {(r["a_t"], tuple(r["b"]), r["branch"]) for r in claims["branches"]}
    if expected_open != claim_keys or len(claims["branches"]) != 420:
        raise AssertionError("claim records do not equal the ledger's open branches")

    cache, results = {}, []
    for number, row in enumerate(claims["branches"], 1):
        branch = row["branch"]
        terminal = 7 if branch == "T1" else 6
        levels = tuple(range(4, terminal + 1))
        caps = (*VAR_CAPS, *((10 + 3 * row["a_t"],) * len(levels)))
        found = found_flags = None
        for flags in flag_cases(branch):
            candidates = []
            for b in row["b"]:
                key = (branch, row["a_t"], b, flags)
                if key not in cache:
                    cache[key] = local_profiles(branch, row["a_t"], b, flags, monomial_data)
                candidates.append(cache[key])
                if not cache[key]:
                    break
            if len(candidates) != 4 or not candidates[-1]:
                continue
            found = join_places(candidates, caps)
            if found is not None:
                found_flags = flags
                break

        audit_status = "survives" if found is not None else "killed"
        claim_status = "killed" if row["status"] == "engine_killed_pending_audit" else "survives"
        result = {
            "id": branch_id(row), "claim": claim_status, "audit": audit_status,
            "agreement": audit_status == claim_status,
            "witness": None if found is None else witness_text(found_flags, found),
        }
        results.append(result)
        if not args.quiet:
            print(f"{number:03d} {result['id']:<34} claim={claim_status:<8} "
                  f"audit={audit_status:<8} {'AGREE' if result['agreement'] else 'DISAGREE'}")

    agreements = sum(r["agreement"] for r in results)
    disagreements = [r for r in results if not r["agreement"]]
    audit_kills = sum(r["audit"] == "killed" for r in results)
    elapsed = time.perf_counter() - started
    print(f"open_summary: {agreements}/420 agree; disagreements={len(disagreements)}; "
          f"audit_killed={audit_kills}; audit_survives={420-audit_kills}")
    for result in disagreements:
        print(f"OPEN_DISAGREEMENT {result['id']} claim={result['claim']} audit={result['audit']}")
        if result["witness"]:
            print("  CONSERVATIVE_WITNESS", result["witness"])
    print(f"runtime_seconds: {elapsed:.3f}")
    if args.emit_artifact:
        emit_artifact(args.emit_artifact, "sub2", results, {
            "total": len(results), "agree": agreements,
            "disagreements": len(disagreements),
            "audit_killed": audit_kills, "audit_survives": 420 - audit_kills,
            "audits": "cascade_cones.json (depth-4 q-place cascade)",
        })
    return 1 if disagreements else 0


if __name__ == "__main__":
    raise SystemExit(main())
