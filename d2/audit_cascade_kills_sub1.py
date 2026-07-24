#!/usr/bin/env python3
"""Independent exhaustive valuation audit for subcase (1)."""
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
VAR_CAPS = (6, 9, 12)  # d2, d1, sigma
E_CAP = 15
STANDARD_MAX_A = 10
TERMINALS = {"T1": (7, 46), "T2": (6, 48)}
_NO_VALUE = object()

@dataclass(frozen=True)
class Flags:
    d2_zero: bool
    sigma_zero: bool
    g_zero: frozenset[int]

    def label(self) -> str:
        gs = ",".join(map(str, sorted(self.g_zero))) or "-"
        return f"d2_zero={self.d2_zero},sigma_zero={self.sigma_zero},g_zero={gs}"

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

def load_h_monomials():
    """Parse h_0..h_7, verify weights and sub1 degree bounds, eliminate d0."""
    text = (HERE / "f31_graded.txt").read_text(encoding="utf-8")
    d2, d1, d0, e, sigma = sp.symbols("d2 d1 d0 dm1 sigma")
    old_vars = (d2, d1, d0, e)
    weights = (2, 3, 4, 5)
    degree_caps = (6, 9, 12, 15)
    pattern = re.compile(
        r"^h_(\d+) \(weight (\d+), dm1-power (\d+)\) = (.+)$", re.MULTILINE
    )
    parsed, metadata = {}, {}
    for match in pattern.finditer(text):
        level, weight, outer_power = map(int, match.group(1, 2, 3))
        parsed[level] = sp.sympify(
            match.group(4), locals={"d2": d2, "d1": d1, "d0": d0, "dm1": e}
        )
        metadata[level] = (weight, outer_power)
    if sorted(parsed) != list(range(8)):
        raise AssertionError(f"expected h_0..h_7, found {sorted(parsed)}")

    result = {}
    for level in range(8):
        expected_weight = 20 - 2 * level
        expected_degree = 60 - 6 * level
        if metadata[level] != (expected_weight, 21 - 3 * level):
            raise AssertionError(f"h_{level}: bad graded metadata {metadata[level]}")
        source = sp.Poly(parsed[level], *old_vars)
        source_degrees = []
        for exponents, _coefficient in source.terms():
            actual_weight = sum(w * exponent for w, exponent in zip(weights, exponents))
            if actual_weight != expected_weight:
                raise AssertionError(
                    f"h_{level}: exponent {exponents} has weight {actual_weight}, "
                    f"expected {expected_weight}"
                )
            source_degrees.append(
                sum(cap * exponent for cap, exponent in zip(degree_caps, exponents))
            )
        if max(source_degrees) != expected_degree:
            raise AssertionError(
                f"h_{level}: derived source degree {max(source_degrees)}, expected {expected_degree}"
            )
        rewritten = sp.Poly(
            sp.expand(parsed[level].subs(d0, (sigma + d2**2) / 4)), d2, d1, sigma, e
        )
        terms = tuple(
            (coefficient, tuple(map(int, exponents)))
            for exponents, coefficient in rewritten.terms() if coefficient
        )
        if max(
            sum(cap * exponent for cap, exponent in zip(degree_caps, exponents))
            for _coefficient, exponents in terms
        ) > expected_degree:
            raise AssertionError(f"h_{level}: d0 elimination exceeded degree cap")
        result[level] = terms
    return result

def derive_g_caps(branch: str, a_t: int):
    """Independently apply the stated forward/backward recurrences."""
    terminal, terminal_cap = TERMINALS[branch]
    v = 30 - 3 * a_t
    ehat_cap = E_CAP - a_t
    if not (0 <= a_t <= STANDARD_MAX_A and v >= 0 and ehat_cap >= 0):
        raise ValueError("g caps are only defined here in the standard regime")
    forward = {1: 60 - v}
    for level in range(1, terminal):
        forward[level + 1] = (
            max(3 * ehat_cap + forward[level], 60 - 2 * level) - v
        )
    backward = {terminal: terminal_cap}
    for level in range(terminal - 1, 0, -1):
        backward[level] = max(v + backward[level + 1], 60 - 2 * level)
    safe = {
        level: min(forward[level], backward[level])
        for level in range(1, terminal + 1)
    }
    return forward, backward, safe

def h_domain(monomials, values, degree_cap: int) -> HDomain | None:
    """All h-valuations under the conservative relaxed-tie policy."""
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
        return HDomain(None, None, True)  # structural zero
    minimum = min(valuations)
    if valuations.count(minimum) == 1:
        return HDomain(minimum, minimum, False) if minimum <= degree_cap else None
    if minimum <= degree_cap:
        return HDomain(minimum, degree_cap, True)
    return HDomain(None, None, True)

def equation_holds(r_next, a_val, b_val) -> bool:
    """Ultrametric feasibility for g_next=A+B at one q-place."""
    if r_next is INF:
        return a_val == b_val
    if a_val is INF:
        return b_val == r_next
    if b_val is INF:
        return a_val == r_next
    return r_next >= a_val if a_val == b_val else r_next == min(a_val, b_val)

def choose_h(domain: HDomain, level: int, r_next, a_val):
    """Return one h in the domain satisfying the equation, or a sentinel."""
    candidates = set()
    if domain.low is not None:
        low_b, high_b = level + domain.low, level + int(domain.high)
        probes = [low_b, high_b]
        for value in (r_next, a_val):
            if value is not INF:
                probes.extend((value, value + 1))
        for b_val in probes:
            if low_b <= b_val <= high_b:
                candidates.add(b_val)
    ordered = sorted(candidates)
    if domain.allow_infinity:
        ordered.append(INF)
    for b_val in ordered:
        if equation_holds(r_next, a_val, b_val):
            return INF if b_val is INF else b_val - level
    return _NO_VALUE

def resource_vector(vd2, vd1, vsigma, gvals, levels):
    return (
        0 if vd2 is INF else vd2,
        0 if vd1 is INF else vd1,
        0 if vsigma is INF else vsigma,
        *(0 if gvals[level] is INF else gvals[level] for level in levels),
    )

def pareto_reduce(profiles):
    """Retain a witness for each componentwise-minimal resource vector."""
    by_resources = {}
    for profile in profiles:
        by_resources.setdefault(profile.resources, profile)
    frontier = []
    for resources in sorted(by_resources, key=lambda item: (sum(item), item)):
        if any(all(a <= b for a, b in zip(old.resources, resources)) for old in frontier):
            continue
        frontier.append(by_resources[resources])
    return frontier

def compact_states(states):
    """Pareto-compact partial paths having the same current g-order."""
    groups = {}
    for gvals, hvals in states:
        current = gvals[-1][1]
        resources = tuple(0 if value is INF else value for _level, value in gvals)
        groups.setdefault(current, {}).setdefault(resources, (gvals, hvals))
    result = []
    for by_resources in groups.values():
        frontier = []
        for resources in sorted(by_resources, key=lambda item: (sum(item), item)):
            if any(all(a <= b for a, b in zip(old, resources)) for old in frontier):
                continue
            frontier.append(resources)
            result.append(by_resources[resources])
    return result

def local_profiles(branch, a_t, b, flags, monomial_data, minimum_level):
    """Exhaust a local chain, using only exact dominance reductions."""
    terminal, _terminal_cap = TERMINALS[branch]
    _forward, _backward, g_caps = derive_g_caps(branch, a_t)
    descent_levels = tuple(range(terminal - 1, minimum_level - 1, -1))
    all_levels = tuple(range(minimum_level, terminal + 1))
    d2_values = (INF,) if flags.d2_zero else range(VAR_CAPS[0] + 1)
    d1_values = range(VAR_CAPS[1] + 1) if branch == "T1" else (INF,)
    sigma_values = (INF,) if flags.sigma_zero else range(VAR_CAPS[2] + 1)
    profiles = []

    for vd2, vd1, vsigma in itertools.product(d2_values, d1_values, sigma_values):
        auxiliary = vd1 if branch == "T1" else vsigma
        if auxiliary is INF:
            continue
        r_terminal = terminal + 2 * auxiliary - 3 * b
        if not 0 <= r_terminal <= g_caps[terminal]:
            continue
        states = [(((terminal, r_terminal),), ())]
        values = (vd2, vd1, vsigma, b)
        for level in descent_levels:
            domain = h_domain(monomial_data[level], values, 60 - 6 * level)
            if domain is None:
                states = []
                break
            expanded = []
            candidates_r = (INF,) if level in flags.g_zero else range(g_caps[level] + 1)
            for gvals, hvals in states:
                r_next = gvals[-1][1]
                for r_here in candidates_r:
                    a_val = INF if r_here is INF else 3 * b + r_here
                    vh = choose_h(domain, level, r_next, a_val)
                    if vh is not _NO_VALUE:
                        expanded.append(
                            (gvals + ((level, r_here),), hvals + ((level, vh),))
                        )
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
    """Complete four-place DFS with memoized failed budget states."""
    order = sorted(range(4), key=lambda index: len(candidates[index]))
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
                suffix = dfs(
                    position + 1,
                    tuple(left - used for used, left in zip(profile.resources, remaining)),
                )
                if suffix is not None:
                    suffix[place] = profile
                    return suffix
        failed.add(key)
        return None

    found = dfs(0, caps)
    return None if found is None else [found[index] for index in range(4)]

def flag_cases(branch: str, minimum_level: int):
    terminal, _terminal_cap = TERMINALS[branch]
    lower = tuple(range(minimum_level, terminal))
    sigma_choices = (False, True) if branch == "T1" else (False,)
    for d2_zero, sigma_zero in itertools.product((False, True), sigma_choices):
        for mask in range(1 << len(lower)):
            zeros = frozenset(
                level for bit, level in enumerate(lower) if mask & (1 << bit)
            )
            yield Flags(d2_zero, sigma_zero, zeros)

def terminal_feasible(branch: str, b_values: tuple[int, ...]):
    """Exhaust all four terminal auxiliary orders independently of the ledger."""
    terminal, g_cap = TERMINALS[branch]
    auxiliary_cap = VAR_CAPS[1] if branch == "T1" else VAR_CAPS[2]
    for orders in itertools.product(range(auxiliary_cap + 1), repeat=4):
        if sum(orders) > auxiliary_cap:
            continue
        gorders = tuple(
            terminal + 2 * order - 3 * b for order, b in zip(orders, b_values)
        )
        if min(gorders) >= 0 and sum(gorders) <= g_cap:
            return True, orders, gorders
    return False, None, None

def fmt_value(value) -> str:
    return "inf" if value is INF else str(value)

def branch_id(row: dict) -> str:
    return f"a={row['a_t']};b={','.join(map(str, row['b']))};{row['branch']}"

def witness_text(flags: Flags, profiles: list[LocalProfile]) -> str:
    chunks = []
    for index, profile in enumerate(profiles, 1):
        gs = ",".join(f"g{level}={fmt_value(value)}" for level, value in profile.gvals)
        hs = ",".join(f"h{level}={fmt_value(value)}" for level, value in profile.hvals)
        chunks.append(
            f"p{index}(d2={fmt_value(profile.vd2)},d1={fmt_value(profile.vd1)},"
            f"sigma={fmt_value(profile.vsigma)};{gs};{hs})"
        )
    return flags.label() + " | " + " ".join(chunks)

def claim_cap_checks(ledger, claims):
    """Check explicit ledger caps and stored finite witnesses against derivation."""
    discrepancies = []
    for stratum in ledger["strata"]:
        if stratum["a_t"] > STANDARD_MAX_A:
            continue
        for branch, (terminal, terminal_cap) in TERMINALS.items():
            item = stratum["branches"][branch]["terminal"]
            if item["g_degree_cap"] != terminal_cap:
                discrepancies.append(
                    f"{branch} level {terminal}: claimed {item['g_degree_cap']}, derived {terminal_cap}"
                )
            expected_aux = VAR_CAPS[1] if branch == "T1" else VAR_CAPS[2]
            if item["auxiliary_degree_cap"] != expected_aux:
                discrepancies.append(
                    f"{branch} auxiliary: claimed {item['auxiliary_degree_cap']}, derived {expected_aux}"
                )
    for row in claims["branches"]:
        _forward, _backward, caps = derive_g_caps(row["branch"], row["a_t"])
        for case in row.get("survivor_cases", []):
            witnesses = case.get("witness", [])
            for level in caps:
                values = [
                    witness.get("v_g", {}).get(str(level), "inf") for witness in witnesses
                ]
                finite = [value for value in values if value != "inf"]
                if any(value > caps[level] for value in finite) or sum(finite) > caps[level]:
                    discrepancies.append(
                        f"{branch_id(row)} level {level}: witness uses {values}, derived cap {caps[level]}"
                    )
    return discrepancies

def audit_claims(claims, monomial_data, quiet):
    cache, results = {}, []
    minimum_level = claims["depth"]
    for number, row in enumerate(claims["branches"], 1):
        branch = row["branch"]
        terminal, _terminal_cap = TERMINALS[branch]
        _forward, _backward, g_caps = derive_g_caps(branch, row["a_t"])
        levels = tuple(range(minimum_level, terminal + 1))
        caps = (*VAR_CAPS, *(g_caps[level] for level in levels))
        found = found_flags = None
        for flags in flag_cases(branch, minimum_level):
            candidates = []
            for b in row["b"]:
                key = (branch, row["a_t"], b, flags, minimum_level)
                if key not in cache:
                    cache[key] = local_profiles(
                        branch, row["a_t"], b, flags, monomial_data, minimum_level
                    )
                if not cache[key]:
                    break
                candidates.append(cache[key])
            if len(candidates) != 4:
                continue
            found = join_places(candidates, caps)
            if found is not None:
                found_flags = flags
                break
        audit_status = "survives" if found is not None else "killed"
        claim_status = (
            "killed" if row["status"] == "engine_killed_pending_audit" else "survives"
        )
        result = {
            "id": branch_id(row), "claim": claim_status, "audit": audit_status,
            "agreement": audit_status == claim_status,
            "witness": None if found is None else witness_text(found_flags, found),
        }
        results.append(result)
        if not quiet:
            verdict = "AGREE" if result["agreement"] else "DISAGREE"
            print(
                f"{number:04d} {result['id']:<38} claim={claim_status:<8} "
                f"audit={audit_status:<8} {verdict}"
            )
    return results

def emit_artifact(path: str, window: str, results: list, summary: dict) -> None:
    """Write per-branch audit verdicts so the coverage proof-DAG can machine-join
    this independent audit.  Deterministic (sorted keys, no timestamps)."""
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
    parser.add_argument("--quiet", action="store_true", help="omit the 2178-row table")
    parser.add_argument("--emit-artifact", nargs="?",
                        const="audit_cascade_kills_sub1.json", default=None,
                        metavar="PATH",
                        help="write per-branch audit verdicts as JSON (for the proof-DAG join)")
    args = parser.parse_args()
    started = time.perf_counter()

    monomial_data = load_h_monomials()
    print("weighted_homogeneity: PASS (h_0..h_7 parsed independently)")
    print("sub1_h_caps: PASS (deg h_f <= 60-6f, checked on source and rewrite)")
    print("d0_elimination: PASS (substituted d0=(sigma+d2^2)/4)")
    for branch in TERMINALS:
        for a_t in range(STANDARD_MAX_A + 1):
            derive_g_caps(branch, a_t)
    print("g_cap_derivation: PASS (forward/backward recurrences checked for a=0..10)")

    ledger = json.loads((HERE / "split_place_ledger_sub1.json").read_text(encoding="utf-8"))
    claims = json.loads(
        (HERE / "cascade_cones_sub1_depth4.json").read_text(encoding="utf-8")
    )
    depth5 = json.loads(
        (HERE / "cascade_cones_sub1_depth5.json").read_text(encoding="utf-8")
    )
    if claims.get("depth") != 4 or claims.get("partial_checkpoint"):
        raise AssertionError("depth-4 claims are not a complete depth-4 file")
    if depth5.get("depth") != 5 or depth5.get("partial_checkpoint"):
        raise AssertionError("depth-5 claims are not a complete depth-5 file")

    terminal_total = terminal_agree = 0
    terminal_disagreements = []
    terminal_cache = {}
    for stratum in ledger["strata"]:
        if stratum["a_t"] > STANDARD_MAX_A:
            continue
        b_values = tuple(stratum["b"])
        for branch in TERMINALS:
            key = (branch, b_values)
            if key not in terminal_cache:
                terminal_cache[key] = terminal_feasible(branch, b_values)
            actual, _orders, _gorders = terminal_cache[key]
            claimed = bool(stratum["branches"][branch]["terminal"]["feasible"])
            terminal_total += 1
            terminal_agree += actual == claimed
            if actual != claimed:
                terminal_disagreements.append(
                    f"a={stratum['a_t']};b={','.join(map(str, b_values))};{branch}: "
                    f"claimed={'feasible' if claimed else 'infeasible'}, "
                    f"audit={'feasible' if actual else 'infeasible'}"
                )
    print(
        f"terminal_crosscheck: {terminal_agree}/{terminal_total} agree; "
        f"disagreements={len(terminal_disagreements)}"
    )
    for item in terminal_disagreements:
        print("TERMINAL_DISAGREEMENT", item)

    expected_open = {
        (stratum["a_t"], tuple(stratum["b"]), branch)
        for stratum in ledger["strata"] if stratum["a_t"] <= STANDARD_MAX_A
        for branch in stratum["open_branches"]
    }
    claim_keys = {
        (row["a_t"], tuple(row["b"]), row["branch"]) for row in claims["branches"]
    }
    if expected_open != claim_keys or len(claims["branches"]) != 2178:
        raise AssertionError("depth-4 records do not equal standard-regime ledger branches")

    cap_discrepancies = claim_cap_checks(ledger, claims)
    print(f"claim_cap_crosscheck: discrepancies={len(cap_discrepancies)}")
    for item in cap_discrepancies:
        print("CAP_DISCREPANCY", item)

    results = audit_claims(claims, monomial_data, args.quiet)
    agreements = sum(result["agreement"] for result in results)
    disagreements = [result for result in results if not result["agreement"]]
    audit_kills = sum(result["audit"] == "killed" for result in results)
    print(
        f"depth4_summary: {agreements}/2178 agree; disagreements={len(disagreements)}; "
        f"audit_killed={audit_kills}; audit_survives={2178-audit_kills}"
    )
    for result in disagreements:
        print(
            f"DEPTH4_DISAGREEMENT {result['id']} claim={result['claim']} "
            f"audit={result['audit']}"
        )
        if result["witness"]:
            print("  CONSERVATIVE_WITNESS", result["witness"])

    depth5_keys = {
        (row["a_t"], tuple(row["b"]), row["branch"]): row
        for row in depth5["branches"]
    }
    if set(depth5_keys) != claim_keys:
        raise AssertionError("depth-5 and depth-4 branch sets differ")
    depth4_kills = {
        (row["a_t"], tuple(row["b"]), row["branch"])
        for row in claims["branches"]
        if row["status"] == "engine_killed_pending_audit"
    }
    depth5_kills = {
        key for key, row in depth5_keys.items()
        if row["status"] == "engine_killed_pending_audit"
    }
    subset_violations = sorted(depth5_kills - depth4_kills)
    print(
        f"depth5_subset_check: {len(depth5_kills)-len(subset_violations)}/{len(depth5_kills)} "
        f"depth-5 kills occur among depth-4 kills; violations={len(subset_violations)}"
    )
    for a_t, b_values, branch in subset_violations:
        print(
            "DEPTH5_SUBSET_VIOLATION "
            f"a={a_t};b={','.join(map(str, b_values))};{branch}"
        )

    elapsed = time.perf_counter() - started
    print(f"runtime_seconds: {elapsed:.3f}")
    failed = terminal_disagreements or disagreements or cap_discrepancies or subset_violations
    if args.emit_artifact:
        emit_artifact(args.emit_artifact, "sub1", results, {
            "total": len(results), "agree": agreements,
            "disagreements": len(disagreements),
            "audit_killed": audit_kills, "audit_survives": 2178 - audit_kills,
            "audits": "cascade_cones_sub1_depth4.json (depth-4 q-place cascade, standard regime)",
        })
    return 1 if failed else 0

if __name__ == "__main__":
    raise SystemExit(main())
