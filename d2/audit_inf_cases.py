#!/usr/bin/env python3
"""Independent specification-only audit of max-plus infinity artifacts.

Semantics come only from CASCADE_INF_REPORT.md, cited cap facts,
T5_T2_COLUMN.md section 1, and f31_graded.txt. This file neither imports nor
reads cascade_engine.py, test_cascade_inf.py, or cascade_signature.py.

Finite kills proved by the independent q+t lane are cited via its artifacts.
Every q+t survivor removed at infinity is rechecked with an independently
enumerated finite descent and five-place budget join.  The infinity side is
conservatively relaxed: any maximum tie may drop arbitrarily or vanish, and
residue-tie kills are ignored. This can reject a sharp kill but cannot falsely
confirm one by using stronger cancellation semantics.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import re
import time
from dataclasses import dataclass
from pathlib import Path

import sympy as sp

HERE = Path(__file__).resolve().parent
NEG_INF = None
TERMINAL = {"T1": 7, "T2": 6}

# Documented expectation for the --emit-artifact join (C43 / PROOF_INVENTORY.md:
# "sub1 ... 108 NEW branch kills", plus the four sub2 T2 column cells).  These
# are the branches whose emptiness is established AT the infinity layer and
# nowhere earlier -- i.e. exactly the branches the emitted artifact is able to
# support.  Pinned so that an artifact which supports NOTHING (or a different
# set) cannot be emitted silently and read as a successful join.  Only checked
# under --emit-artifact; the audit itself is unaffected.
EXPECTED_INF_ONLY_KILLS = {"sub2": 4, "sub1": 108}


@dataclass(frozen=True)
class Window:
    name: str
    artifact: str
    finite_artifact: str
    var_caps: tuple[int, int, int]
    e_cap: int
    scale: int
    terminal_caps: tuple[int, int]

    def h_cap(self, level: int) -> int:
        return self.scale * (20 - 2 * level)


WINDOWS = (
    Window("sub2", "cascade_cones_qt_inf_rl.json", "cascade_cones_qt_rl.json",
           (4, 6, 8), 10, 2, (40, 40)),
    # The cited sub1 module docstring corrects these terminal caps to 46/48;
    # the other caps are the weight-derived sub1 window.
    Window("sub1", "cascade_cones_sub1_qt_inf_rl.json",
           "cascade_cones_sub1_qt_rl.json", (6, 9, 12), 15, 3, (46, 48)),
)


@dataclass(frozen=True)
class Flags:
    d2_zero: bool
    sigma_zero: bool
    g_zero: frozenset[int]

    def label(self) -> str:
        gs = ",".join(map(str, sorted(self.g_zero))) or "-"
        return (f"d2_zero={self.d2_zero},sigma_zero={self.sigma_zero},"
                f"g_zero={gs}")


def branch_key(row: dict) -> tuple[int, tuple[int, ...], str]:
    return int(row["a_t"]), tuple(map(int, row["b"])), str(row["branch"])


def branch_label(key) -> str:
    a_t, bs, branch = key
    return f"a={a_t};b={','.join(map(str, bs))};{branch}"


def case_flags(case: dict) -> Flags:
    return Flags(bool(case["d2_zero"]), bool(case["sigma_zero"]),
                 frozenset(map(int, case["g_zero_levels"])))


def ext_int(value, context: str) -> int | None:
    if value in ("inf", "-inf"):
        return NEG_INF
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{context}: expected integer/infinity, got {value!r}")
    return value


def load_h_monomials():
    """Own h-line parser, source/rewrite homogeneity checks, and d0 rewrite."""
    text = (HERE / "f31_graded.txt").read_text(encoding="utf-8-sig")
    d2, d1, d0, e, sigma = sp.symbols("d2 d1 d0 dm1 sigma")
    old_vars, new_vars = (d2, d1, d0, e), (d2, d1, sigma, e)
    weights = (2, 3, 4, 5)
    pattern = re.compile(
        r"^h_(\d+) \(weight (\d+), dm1-power (\d+)\) = (.+)$", re.MULTILINE)
    expressions, headers = {}, {}
    for match in pattern.finditer(text):
        level, weight, power = map(int, match.group(1, 2, 3))
        if level in expressions:
            raise ValueError(f"duplicate h_{level}")
        expressions[level] = sp.sympify(
            match.group(4), locals={"d2": d2, "d1": d1, "d0": d0, "dm1": e})
        headers[level] = (weight, power)
    if sorted(expressions) != list(range(8)):
        raise ValueError(f"expected h_0..h_7, got {sorted(expressions)}")
    result = {}
    for level in range(8):
        expected_weight = 20 - 2 * level
        if headers[level] != (expected_weight, 21 - 3 * level):
            raise ValueError(f"h_{level}: bad header {headers[level]}")
        source = sp.Poly(expressions[level], *old_vars)
        for exponents, _coefficient in source.terms():
            actual = sum(w * n for w, n in zip(weights, exponents))
            if actual != expected_weight:
                raise ValueError(f"h_{level}: source weight {actual} != {expected_weight}")
        rewritten = sp.Poly(
            sp.expand(expressions[level].subs(d0, (sigma + d2**2) / 4)), *new_vars)
        terms = tuple((coefficient, tuple(map(int, exponents)))
                      for exponents, coefficient in rewritten.terms() if coefficient)
        if not terms:
            raise ValueError(f"h_{level}: rewritten polynomial is zero")
        for _coefficient, exponents in terms:
            actual = sum(w * n for w, n in zip(weights, exponents))
            if actual != expected_weight:
                raise ValueError(f"h_{level}: rewritten weight {actual} != {expected_weight}")
        for window in WINDOWS:
            caps = (*window.var_caps, window.e_cap)
            maximum = max(sum(c * n for c, n in zip(caps, exponents))
                          for _coefficient, exponents in terms)
            if maximum != window.h_cap(level):
                raise ValueError(f"{window.name} h_{level}: cap {maximum} != "
                                 f"{window.h_cap(level)}")
        result[level] = terms
    return result


def g_caps(window: Window, branch: str, a_t: int) -> dict[int, int]:
    """Derive the stated min(forward, backward) cap at every level."""
    if not 0 <= a_t <= 10:
        raise ValueError(f"nonstandard a={a_t}")
    terminal = TERMINAL[branch]
    if window.name == "sub2":
        return {level: 10 + 3 * a_t for level in range(1, terminal + 1)}
    v, ehat_cap = 30 - 3 * a_t, window.e_cap - a_t
    forward = {1: window.h_cap(0) - v}
    for level in range(1, terminal):
        forward[level + 1] = (
            max(3 * ehat_cap + forward[level],
                4 * level + window.h_cap(level)) - v)
    cap = window.terminal_caps[0 if branch == "T1" else 1]
    backward = {terminal: cap}
    for level in range(terminal - 1, 0, -1):
        backward[level] = max(v + backward[level + 1],
                              4 * level + window.h_cap(level))
    return {level: min(forward[level], backward[level])
            for level in range(1, terminal + 1)}


def term_label(term) -> str:
    coefficient, (i, j, k, m) = term
    return f"{coefficient}*d2^{i}*d1^{j}*sigma^{k}*e^{m}"


def active_terms(terms, degrees):
    active = []
    for term in terms:
        total = 0
        for exponent, degree in zip(term[1], degrees):
            if exponent and degree is NEG_INF:
                break
            if exponent:
                total += exponent * int(degree)
        else:
            active.append((total, term))
    return active


def norm_obligation(item: dict) -> dict:
    return {key: item.get(key) for key in ("kind", "level", "depth", "tied")}


def exact_infinity(branch: str, a_t: int, degrees, deg_g, recorded,
                   monomials) -> list[str]:
    """Reconstruct the anchor, all max-plus equations, and obligations."""
    errors, expected = [], []
    terminal, v, e_degree = TERMINAL[branch], 30 - 3 * a_t, degrees[3]
    by_level = {}
    for raw in recorded:
        item = norm_obligation(raw)
        if not isinstance(item["level"], int):
            errors.append(f"malformed obligation {raw!r}")
            continue
        by_level.setdefault(item["level"], []).append(item)

    def choose_h(level: int):
        active = active_terms(monomials[level], degrees)
        if not active:
            return NEG_INF, [], False
        maximum = max(value for value, _term in active)
        tied = [term_label(term) for value, term in active if value == maximum]
        obligations = [item for item in by_level.get(level, [])
                       if item["kind"] in ("degree_tie_drop",
                                            "identical_vanishing")]
        if len(obligations) > 1:
            errors.append(f"level {level}: multiple h obligations")
            return maximum, tied, False
        if not obligations:
            return maximum, tied, False
        item = obligations[0]
        if len(tied) < 2:
            errors.append(f"level {level}: {item['kind']} without maximum tie")
        if item["tied"] != tied:
            errors.append(f"level {level}: tied monomials recorded={item['tied']!r}, "
                          f"expected={tied!r}")
        expected.append(item)
        if item["kind"] == "identical_vanishing":
            if item["depth"] != 0:
                errors.append(f"level {level}: identical_vanishing depth != 0")
            return NEG_INF, tied, False
        depth = item["depth"]
        if isinstance(depth, bool) or not isinstance(depth, int) or depth <= 0:
            errors.append(f"level {level}: drop depth is not positive")
            return maximum, tied, False
        if maximum - depth < 0:
            errors.append(f"level {level}: drop below degree zero")
        return maximum - depth, tied, True

    target = NEG_INF if deg_g[1] is NEG_INF else v + int(deg_g[1])
    h0, _top, _dropped = choose_h(0)
    if h0 != target:
        errors.append(f"anchor deg(h0)={h0}, expected {target}")

    for level in range(1, terminal):
        h_degree, h_top, h_dropped = choose_h(level)
        left = NEG_INF if deg_g[level + 1] is NEG_INF else v + int(deg_g[level + 1])
        first = (NEG_INF if deg_g[level] is NEG_INF
                 else 3 * (e_degree - a_t) + int(deg_g[level]))
        second = NEG_INF if h_degree is NEG_INF else 4 * level + h_degree
        sum_obs = [item for item in by_level.get(level, [])
                   if item["kind"] in ("leading_cancellation", "exact_identity")]
        wanted = None
        h_form = (f"lc(h{level}@deg={h_degree})" if h_dropped
                  else f"[{' + '.join(h_top)}]")
        sides = [f"lc(ehat)^3*lc(g{level})",
                 f"(-1024/3315)^{level}*{h_form}"]
        if left is NEG_INF:
            if first != second:
                errors.append(f"level {level}: zero left, side degrees {first},{second}")
            elif first is not NEG_INF:
                wanted = {"kind": "exact_identity", "level": level,
                          "depth": 0, "tied": sides}
        elif first is NEG_INF:
            if second != left:
                errors.append(f"level {level}: h-side {second} != left {left}")
        elif second is NEG_INF:
            if first != left:
                errors.append(f"level {level}: ehat-side {first} != left {left}")
        elif first != second:
            if left != max(first, second):
                errors.append(f"level {level}: unique max {max(first, second)} != {left}")
        elif left > first:
            errors.append(f"level {level}: left {left} exceeds tied sides {first}")
        elif left < first:
            wanted = {"kind": "leading_cancellation", "level": level,
                      "depth": first - left, "tied": sides}
        if wanted is not None:
            expected.append(wanted)
            if sum_obs != [wanted]:
                errors.append(f"level {level}: sum obligation recorded={sum_obs!r}, "
                              f"expected={[wanted]!r}")
        elif sum_obs:
            errors.append(f"level {level}: unexpected sum obligation {sum_obs!r}")

    normalized = [norm_obligation(item) for item in recorded]
    if sorted(normalized, key=repr) != sorted(expected, key=repr):
        errors.append(f"obligation set recorded={normalized!r}, expected={expected!r}")
    return errors


def squeeze_enabled(metadata: dict) -> bool:
    """The artifacts expose the squeeze through residue_kills=true."""
    return metadata.get("places") == "q+t+inf" and metadata.get("residue_kills") is True


def validate_survivor(window: Window, row: dict, case: dict, monomials,
                      squeeze: bool) -> list[str]:
    errors = []
    key, flags = branch_key(row), case_flags(case)
    a_t, bs, branch = key
    label = f"{branch_label(key)} {flags.label()}"
    terminal, caps = TERMINAL[branch], g_caps(window, branch, a_t)
    records = case.get("witness")
    if not isinstance(records, list) or len(records) != 6:
        return [f"{label}: witness does not have four q, t, infinity records"]
    q_records, t_record, inf_record = records[:4], records[4], records[5]
    if [record.get("place") for record in q_records] != ["q"] * 4:
        errors.append(f"{label}: first four records are not q places")
    if t_record.get("place") != "t" or inf_record.get("place") != "inf":
        errors.append(f"{label}: final records are not t, infinity")

    finite = []
    for index, record in enumerate(q_records + [t_record]):
        place = "q" if index < 4 else "t"
        parameter = bs[index] if index < 4 else a_t
        if record.get("b") != parameter:
            errors.append(f"{label}: {place} parameter {record.get('b')} != {parameter}")
        values = {}
        for position, field in enumerate(("v_d2", "v_d1", "v_sigma")):
            try:
                value = ext_int(record.get(field), f"{label} {place} {field}")
            except ValueError as exc:
                errors.append(str(exc)); value = NEG_INF
            values[field] = value
            if value is not NEG_INF and not 0 <= value <= window.var_caps[position]:
                errors.append(f"{label}: {place} {field}={value} outside cap")
        raw_g = record.get("v_g", {})
        if set(raw_g) != {str(level) for level in range(4, terminal + 1)}:
            errors.append(f"{label}: {place} has wrong g-level keys")
        for level in range(4, terminal + 1):
            try:
                value = ext_int(raw_g.get(str(level)),
                                f"{label} {place} v_g{level}")
            except ValueError as exc:
                errors.append(str(exc)); value = NEG_INF
            values[f"g{level}"] = value
            if value is not NEG_INF and not 0 <= value <= caps[level]:
                errors.append(f"{label}: {place} v_g{level}={value} outside cap")
        auxiliary = values["v_d1"] if branch == "T1" else values["v_sigma"]
        g_terminal = values[f"g{terminal}"]
        if auxiliary is NEG_INF or g_terminal is NEG_INF:
            errors.append(f"{label}: {place} terminal data are infinite")
        else:
            wanted = (terminal + 2 * auxiliary - 3 * parameter
                      if place == "q" else 2 * auxiliary)
            if g_terminal != wanted:
                errors.append(f"{label}: {place} terminal g={g_terminal}, expected {wanted}")
        finite.append(values)

    try:
        variable_degrees = tuple(ext_int(inf_record.get(field), f"{label} {field}")
                                 for field in ("deg_d2", "deg_d1", "deg_sigma"))
        e_degree = ext_int(inf_record.get("deg_e"), f"{label} deg_e")
    except ValueError as exc:
        return errors + [str(exc)]
    if e_degree is NEG_INF:
        errors.append(f"{label}: e is zero"); e_degree = 0
    expected_zero = (flags.d2_zero, branch == "T2", flags.sigma_zero)
    for name, degree, cap, is_zero in zip(("d2", "d1", "sigma"), variable_degrees,
                                          window.var_caps, expected_zero):
        if is_zero and degree is not NEG_INF:
            errors.append(f"{label}: zero {name} has degree {degree}")
        if not is_zero and (degree is NEG_INF or not 0 <= degree <= cap):
            errors.append(f"{label}: nonzero {name} degree {degree} outside cap")
    if not a_t + sum(bs) <= e_degree <= window.e_cap:
        errors.append(f"{label}: deg(e)={e_degree} outside sandwich")

    raw_g = inf_record.get("deg_g", {})
    if set(raw_g) != {str(level) for level in range(1, terminal + 1)}:
        errors.append(f"{label}: infinity has wrong g-level keys")
    deg_g = {}
    for level in range(1, terminal + 1):
        try:
            degree = ext_int(raw_g.get(str(level)), f"{label} deg_g{level}")
        except ValueError as exc:
            errors.append(str(exc)); degree = NEG_INF
        deg_g[level] = degree
        if level >= 4 and ((degree is NEG_INF) != (level in flags.g_zero)):
            errors.append(f"{label}: deg(g{level}) disagrees with zero flag")
        if degree is not NEG_INF and not 0 <= degree <= caps[level]:
            errors.append(f"{label}: deg(g{level})={degree} outside cap")
    if deg_g[terminal] is NEG_INF:
        errors.append(f"{label}: terminal g is zero")

    for name, degree in zip(("d2", "d1", "sigma"), variable_degrees):
        total = sum(value for values in finite
                    if (value := values[f"v_{name}"]) is not NEG_INF)
        if degree is not NEG_INF and total > degree:
            errors.append(f"{label}: finite {name} sum {total} > degree {degree}")
    for level in range(4, terminal + 1):
        total = sum(value for values in finite
                    if (value := values[f"g{level}"]) is not NEG_INF)
        if deg_g[level] is not NEG_INF and total > deg_g[level]:
            errors.append(f"{label}: finite g{level} sum {total} > degree {deg_g[level]}")

    auxiliary_degree = variable_degrees[1 if branch == "T1" else 2]
    if auxiliary_degree is not NEG_INF and deg_g[terminal] is not NEG_INF:
        lhs = 3 * (e_degree - a_t) + deg_g[terminal]
        rhs = 4 * terminal + 2 * auxiliary_degree
        if lhs != rhs:
            errors.append(f"{label}: terminal degree identity {lhs} != {rhs}")
    if squeeze and branch == "T2" and 2 not in bs and deg_g[6] is not NEG_INF:
        finite_g6 = sum(int(values["g6"]) for values in finite)
        bound = deg_g[6] - 2 * (e_degree - a_t - sum(bs))
        if finite_g6 > bound:
            errors.append(f"{label}: squeeze {finite_g6} > {bound}")

    errors.extend(f"{label}: {error}" for error in exact_infinity(
        branch, a_t, (*variable_degrees, e_degree), deg_g,
        inf_record.get("obligations", []), monomials))
    actual_count = sum(len(record.get("obligations", [])) for record in records)
    if case.get("obligation_count") != actual_count:
        errors.append(f"{label}: obligation_count={case.get('obligation_count')} != "
                      f"witness count {actual_count}")
    return errors


def terminal_profiles(window: Window, branch: str, a_t: int, bs: tuple[int, ...]):
    """Own finite enumeration, conservatively projected to terminal totals."""
    terminal = TERMINAL[branch]
    aux_cap = window.var_caps[1 if branch == "T1" else 2]
    g_cap = g_caps(window, branch, a_t)[terminal]
    states = {(0, 0)}
    for place, parameter in [("q", b) for b in bs] + [("t", a_t)]:
        options = []
        for auxiliary in range(aux_cap + 1):
            g_value = (terminal + 2 * auxiliary - 3 * parameter
                       if place == "q" else 2 * auxiliary)
            if 0 <= g_value <= g_cap:
                options.append((auxiliary, g_value))
        states = {(old_a + new_a, old_g + new_g)
                  for old_a, old_g in states for new_a, new_g in options
                  if old_a + new_a <= aux_cap and old_g + new_g <= g_cap}
        if not states:
            return ()
    frontier = []
    for candidate in sorted(states, key=lambda item: (sum(item), item)):
        if any(a <= candidate[0] and g <= candidate[1] for a, g in frontier):
            continue
        frontier.append(candidate)
    return tuple(frontier)


def local_h_domain(terms, values, cap: int):
    active = active_terms(terms, values)
    if not active:
        return None, None, True
    minimum = min(value for value, _term in active)
    if sum(value == minimum for value, _term in active) == 1:
        return (minimum, minimum, False) if minimum <= cap else (None, None, False)
    return ((minimum, cap, True) if minimum <= cap else (None, None, True))


def min_plus_domain_holds(left, first, shift: int, domain) -> bool:
    low, high, allow_infinity = domain
    if left is NEG_INF:
        if first is NEG_INF:
            return allow_infinity
        return low is not None and shift + low <= first <= shift + high
    if first is NEG_INF:
        return low is not None and shift + low <= left <= shift + high
    if allow_infinity and first == left:
        return True
    if low is None:
        return False
    low, high = shift + low, shift + high
    if low <= first <= high and left >= first:
        return True
    if first == left and high > first:
        return True
    return left < first and low <= left <= high


def pareto(vectors):
    frontier = []
    value = lambda item: tuple(0 if entry is NEG_INF else entry for entry in item)
    for vector in sorted(set(vectors), key=lambda item: (sum(value(item)), value(item))):
        if any(all(old <= new for old, new in zip(value(saved), value(vector)))
               for saved in frontier):
            continue
        frontier.append(vector)
    return tuple(frontier)


def compact_partial_states(states):
    """Pareto-reduce prefixes only when their current g-value is identical."""
    groups = {}
    for state in states:
        groups.setdefault(state[-1], []).append(state)
    return {state for group in groups.values() for state in pareto(group)}


_LOCAL_CACHE = {}
_JOIN_CACHE = {}


def local_profiles(window: Window, branch: str, a_t: int, parameter: int,
                   flags: Flags, monomials, place: str):
    cache_key = (window.name, branch, a_t, parameter, flags, place)
    if cache_key in _LOCAL_CACHE:
        return _LOCAL_CACHE[cache_key]
    terminal, caps = TERMINAL[branch], g_caps(window, branch, a_t)
    d2s = (NEG_INF,) if flags.d2_zero else range(window.var_caps[0] + 1)
    d1s = (NEG_INF,) if branch == "T2" else range(window.var_caps[1] + 1)
    sigmas = ((NEG_INF,) if flags.sigma_zero
              else range(window.var_caps[2] + 1))
    profiles = []
    for vd2, vd1, vsigma in itertools.product(d2s, d1s, sigmas):
        auxiliary = vd1 if branch == "T1" else vsigma
        if auxiliary is NEG_INF:
            continue
        g_terminal = (terminal + 2 * auxiliary - 3 * parameter
                      if place == "q" else 2 * auxiliary)
        if not 0 <= g_terminal <= caps[terminal]:
            continue
        states = {(g_terminal,)}
        h_values = (vd2, vd1, vsigma, parameter if place == "q" else a_t)
        for level in range(terminal - 1, 3, -1):
            h_domain = local_h_domain(monomials[level], h_values,
                                      window.h_cap(level))
            if h_domain == (None, None, False):
                states = set(); break
            candidates = ((NEG_INF,) if level in flags.g_zero
                          else range(caps[level] + 1))
            valid_by_next = {}
            for next_value in {state[-1] for state in states}:
                left = (next_value if place == "q"
                        else (NEG_INF if next_value is NEG_INF
                              else 30 - 3 * a_t + next_value))
                valid = []
                for here in candidates:
                    first = (NEG_INF if here is NEG_INF else
                             ((3 * parameter + here) if place == "q" else here))
                    if min_plus_domain_holds(
                            left, first, level if place == "q" else 0, h_domain):
                        valid.append(here)
                valid_by_next[next_value] = valid
            expanded = {state + (here,) for state in states
                        for here in valid_by_next[state[-1]]}
            states = compact_partial_states(expanded)
            if not states:
                break
        for state in states:
            g_values = dict(zip(range(terminal, 3, -1), state))
            profiles.append((0 if vd2 is NEG_INF else vd2,
                             0 if vd1 is NEG_INF else vd1,
                             0 if vsigma is NEG_INF else vsigma,
                             *(0 if g_values[level] is NEG_INF else g_values[level]
                               for level in range(4, terminal + 1))))
    result = pareto(profiles)
    _LOCAL_CACHE[cache_key] = result
    return result


def finite_frontier(window: Window, key, flags: Flags, monomials):
    cache_key = (window.name, key, flags)
    if cache_key in _JOIN_CACHE:
        return _JOIN_CACHE[cache_key]
    a_t, bs, branch = key
    terminal, caps = TERMINAL[branch], g_caps(window, branch, a_t)
    candidates = [local_profiles(window, branch, a_t, b, flags, monomials, "q")
                  for b in bs]
    candidates.append(local_profiles(window, branch, a_t, a_t, flags, monomials, "t"))
    if any(not profiles for profiles in candidates):
        _JOIN_CACHE[cache_key] = ()
        return ()
    budget = (*window.var_caps, *(caps[level] for level in range(4, terminal + 1)))
    totals = {(0,) * len(budget)}
    for profiles in sorted(candidates, key=len):
        totals = pareto(tuple(tuple(a + b for a, b in zip(old, profile))
                              for old, profile in itertools.product(totals, profiles)
                              if all(a + b <= cap
                                     for a, b, cap in zip(old, profile, budget))))
        if not totals:
            break
    _JOIN_CACHE[cache_key] = totals
    return totals


def relaxed_h_values(terms, degrees):
    active = active_terms(terms, degrees)
    if not active:
        return (NEG_INF,)
    maximum = max(value for value, _term in active)
    if sum(value == maximum for value, _term in active) == 1:
        return (maximum,)
    # Sound relaxation: every tied maximum may drop arbitrarily or vanish.
    return (*range(maximum + 1), NEG_INF)


def max_plus_holds(left, first, second) -> bool:
    if left is NEG_INF:
        return first == second
    if first is NEG_INF:
        return second == left
    if second is NEG_INF:
        return first == left
    return left <= first if first == second else left == max(first, second)


def relaxed_survives(window: Window, key, flags: Flags, monomials,
                     squeeze: bool, resources=None) -> bool:
    """Search an over-approximation of all finite+infinity witnesses."""
    a_t, bs, branch = key
    terminal, caps = TERMINAL[branch], g_caps(window, branch, a_t)
    profiles = (terminal_profiles(window, branch, a_t, bs)
                if resources is None else resources)
    if not profiles:
        return False
    v, e_min = 30 - 3 * a_t, a_t + sum(bs)
    d2s = (NEG_INF,) if flags.d2_zero else range(window.var_caps[0] + 1)
    d1s = (NEG_INF,) if branch == "T2" else range(window.var_caps[1] + 1)
    sigmas = ((NEG_INF,) if flags.sigma_zero
              else range(window.var_caps[2] + 1))
    for d2, d1, sigma, e_degree in itertools.product(
            d2s, d1s, sigmas, range(e_min, window.e_cap + 1)):
        auxiliary = d1 if branch == "T1" else sigma
        if auxiliary is NEG_INF:
            continue
        terminal_degree = (4 * terminal + 2 * auxiliary
                           - 3 * (e_degree - a_t))
        if not 0 <= terminal_degree <= caps[terminal]:
            continue
        fitting = []
        for profile in profiles:
            if resources is None:
                finite_aux, finite_terminal = profile
                if finite_aux > auxiliary or finite_terminal > terminal_degree:
                    continue
            else:
                chosen = (d2, d1, sigma,
                          *(0 if level in flags.g_zero else
                            (terminal_degree if level == terminal else caps[level])
                            for level in range(4, terminal + 1)))
                if any(degree is NEG_INF and used or
                       degree is not NEG_INF and used > degree
                       for used, degree in zip(profile, chosen)):
                    continue
                finite_terminal = profile[-1]
            if squeeze and branch == "T2" and 2 not in bs:
                margin = e_degree - a_t - sum(bs)
                if finite_terminal > terminal_degree - 2 * margin:
                    continue
            fitting.append(profile)
        if not fitting:
            continue
        degrees = (d2, d1, sigma, e_degree)
        for profile in fitting:
            lower_g = ({level: 0 for level in range(4, terminal + 1)}
                       if resources is None else
                       dict(zip(range(4, terminal + 1), profile[3:])))
            next_options = {terminal_degree}
            for level in range(terminal - 1, 0, -1):
                if level >= 4 and level in flags.g_zero:
                    here_domain = (NEG_INF,)
                else:
                    here_domain = range(lower_g.get(level, 0), caps[level] + 1)
                h_values = relaxed_h_values(monomials[level], degrees)
                here_options = set()
                for next_degree in next_options:
                    left = NEG_INF if next_degree is NEG_INF else v + next_degree
                    for here in here_domain:
                        first = (NEG_INF if here is NEG_INF
                                 else 3 * (e_degree - a_t) + here)
                        if any(max_plus_holds(
                                left, first,
                                NEG_INF if h is NEG_INF else 4 * level + h)
                               for h in h_values):
                            here_options.add(here)
                if not here_options:
                    break
                next_options = here_options
            else:
                h0_values = relaxed_h_values(monomials[0], degrees)
                if any((NEG_INF if g1 is NEG_INF else v + g1) in h0_values
                       for g1 in next_options):
                    return True
    return False


def validate_metadata(window: Window, artifact: dict, finite: dict) -> list[str]:
    errors = []
    if artifact.get("window") != window.name:
        errors.append(f"{window.name}: wrong window metadata")
    if artifact.get("places") != "q+t+inf" or artifact.get("depth") != 4:
        errors.append(f"{window.name}: wrong places/depth metadata")
    rows = artifact.get("branches")
    if not isinstance(rows, list):
        return errors + [f"{window.name}: branches is not a list"]
    keys = [branch_key(row) for row in rows]
    if len(keys) != len(set(keys)):
        errors.append(f"{window.name}: duplicate branches")
    finite_keys = {branch_key(row) for row in finite.get("branches", [])}
    for key in set(keys) - finite_keys:
        errors.append(f"{branch_label(key)}: absent from finite artifact")
    summary = artifact.get("summary", {})
    survives = sum(row.get("status") == "survives" for row in rows)
    expected = {"open_branches_processed": len(rows),
                "surviving_branches": survives,
                "engine_killed_pending_audit": len(rows) - survives}
    for field, value in expected.items():
        if summary.get(field) != value:
            errors.append(f"{window.name}: summary {field}={summary.get(field)!r}, "
                          f"expected {value}")
    for row in rows:
        cases = row.get("survivor_cases")
        if not isinstance(cases, list):
            errors.append(f"{branch_label(branch_key(row))}: cases is not a list")
            continue
        if row.get("survivor_case_count") != len(cases):
            errors.append(f"{branch_label(branch_key(row))}: case count mismatch")
        flags = [case_flags(case) for case in cases]
        if len(flags) != len(set(flags)):
            errors.append(f"{branch_label(branch_key(row))}: duplicate flag case")
        wanted = "survives" if cases else "engine_killed_pending_audit"
        if row.get("status") != wanted:
            errors.append(f"{branch_label(branch_key(row))}: status/cases mismatch")
    return errors


def branch_record(window: Window, key, row: dict, n_inf: int, n_qt: int,
                  removed_here: int, branch_clean: bool) -> dict:
    """One joinable per-branch verdict for the infinity layer.

    The ONLY verdict this auditor is entitled to assert about a branch is what
    happens at the infinity layer, on top of the q+t_rl survivor set it reads as
    given.  So:

      n_inf > 0                     -> 'survives'      (branch is still open at inf)
      n_inf == 0 and n_qt > 0       -> 'killed'        (kill_layer='inf': every one
                                       of the branch's q+t survivor cases was
                                       removed at infinity, and this auditor
                                       re-derived every one of those removals)
      n_inf == 0 and n_qt == 0      -> 'not_covered'   (kill_layer='pre_inf': the
                                       branch was already empty BEFORE the
                                       infinity layer.  This auditor re-derives
                                       NOTHING about such a branch and must not
                                       claim to -- it is the q- and t-layer
                                       auditors' business.)

    A branch that produced any disagreement is 'disagreement', never 'killed'.
    `agreement` is None for 'not_covered' so a consumer cannot read the absence
    of a disagreement as a confirmation.
    """
    claim = row.get("status")
    if n_inf:
        audit, layer = "survives", None
        agreement = (claim == "survives")
    elif not branch_clean:
        audit, layer, agreement = "disagreement", "inf", False
    elif n_qt:
        audit, layer = "killed", "inf"
        agreement = (claim != "survives")
    else:
        audit, layer, agreement = "not_covered", "pre_inf", None
    return {
        "window": window.name,
        "a_t": key[0], "b": list(key[1]), "branch": key[2],
        "claim": claim, "audit": audit, "agreement": agreement,
        "kill_layer": layer,
        "inf_survivor_cases": n_inf,
        "qt_survivor_cases": n_qt,
        "removed_cases_confirmed": removed_here if branch_clean else 0,
    }


def emit_artifact(path: str, records: list, stats: list) -> None:
    """Write the per-branch infinity-layer audit verdicts so the coverage
    proof-DAG can machine-join this independent audit, exactly as
    audit_cascade_kills{,_sub1}.json are joined.  Deterministic (sorted keys, no
    timestamps).  Both windows go in ONE file; each record names its window.

    Only records with audit=='killed' and kill_layer=='inf' support anything on
    the consumer side: those are the branches the q-cascade auditor verdicts
    'survives' (it sees no q-level kill) and which this auditor re-derives empty
    at the infinity layer.
    """
    recs = sorted(records, key=lambda r: (r["window"], r["a_t"], r["b"], r["branch"]))
    summary, shortfalls = {}, []
    for item in stats:
        if item.get("skipped"):
            continue
        name = item["window"]
        rows = [r for r in recs if r["window"] == name]
        killed = sum(1 for r in rows if r["audit"] == "killed")
        summary[name] = {
            "total": len(rows),
            "audit_survives": sum(1 for r in rows if r["audit"] == "survives"),
            "audit_killed_at_inf": killed,
            "not_covered_killed_before_inf":
                sum(1 for r in rows if r["audit"] == "not_covered"),
            "disagreements": sum(1 for r in rows if r["audit"] == "disagreement"),
            "removed_cases_confirmed":
                sum(r["removed_cases_confirmed"] for r in rows),
            "inf_survivor_cases": item["survivors"],
            "removed_cases_seen": item["removed"],
        }
        want = EXPECTED_INF_ONLY_KILLS.get(name)
        if want is not None and killed != want:
            shortfalls.append(f"{name}: {killed} inf-only branch kills, "
                              f"expected {want}")
        if killed and summary[name]["removed_cases_confirmed"] <= 0:
            shortfalls.append(f"{name}: {killed} branches marked killed at "
                              "infinity but ZERO removed cases were confirmed")
    if shortfalls:
        raise RuntimeError(
            "refusing to emit a join artifact that does not match the pinned "
            "expectation (EXPECTED_INF_ONLY_KILLS): " + "; ".join(shortfalls))
    out = {
        "schema": 1,
        "generator": Path(__file__).name,
        "generator_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "audits": "cascade_cones_qt_inf_rl.json / "
                  "cascade_cones_sub1_qt_inf_rl.json (the infinity layer, on top "
                  "of the q+t_rl survivor sets it reads as given)",
        "windows": sorted(summary),
        "summary": summary,
        "branches": recs,
    }
    target = path if Path(path).is_absolute() else str(HERE / path)
    with open(target, "w") as handle:
        json.dump(out, handle, indent=1, sort_keys=True)
    print(f"emitted inf-audit artifact: {target} ({len(recs)} branches; "
          + ", ".join(f"{w}={summary[w]['audit_killed_at_inf']} inf-only kills"
                      for w in sorted(summary)) + ")")


def audit_window(window: Window, monomials, quiet: bool):
    started = time.perf_counter()
    path = HERE / window.artifact
    if not path.exists():
        return [], {"window": window.name, "skipped": True, "seconds": 0.0}, []
    artifact = json.loads(path.read_text(encoding="utf-8"))
    finite = json.loads((HERE / window.finite_artifact).read_text(encoding="utf-8"))
    errors = validate_metadata(window, artifact, finite)
    finite_rows = {branch_key(row): row for row in finite["branches"]}
    squeeze = squeeze_enabled(artifact)
    survivors = removed = cited = 0
    records = []
    for number, row in enumerate(artifact["branches"], 1):
        key = branch_key(row)
        if key not in finite_rows:
            continue
        before = len(errors)
        target_cases = {case_flags(case): case
                        for case in row.get("survivor_cases", [])}
        finite_cases = {case_flags(case): case
                        for case in finite_rows[key].get("survivor_cases", [])}
        for flags in set(target_cases) - set(finite_cases):
            errors.append(f"{branch_label(key)} {flags.label()}: infinity survivor "
                          "was not a finite q+t survivor")
        for case in target_cases.values():
            survivors += 1
            errors.extend(validate_survivor(window, row, case, monomials, squeeze))
        removed_here = 0
        for flags in set(finite_cases) - set(target_cases):
            removed += 1
            removed_here += 1
            if relaxed_survives(window, key, flags, monomials, squeeze):
                resources = finite_frontier(window, key, flags, monomials)
                if resources and relaxed_survives(
                        window, key, flags, monomials, squeeze, resources):
                    errors.append(f"{branch_label(key)} {flags.label()}: conservative "
                                  "relaxed semantics admits a witness")
        terminal = TERMINAL[key[2]]
        legal_count = 2 * (2 if key[2] == "T1" else 1) * (1 << (terminal - 4))
        cited += legal_count - len(finite_cases)
        records.append(branch_record(window, key, row, len(target_cases),
                                     len(finite_cases), removed_here,
                                     len(errors) == before))
        if not quiet:
            print(f"{window.name} {number:04d} {branch_label(key):<38} "
                  f"survivor_cases={len(target_cases)}")
    elapsed = time.perf_counter() - started
    return errors, {
        "window": window.name, "skipped": False,
        "branches": len(artifact["branches"]), "survivors": survivors,
        "removed": removed, "cited": cited,
        "partial": bool(artifact.get("partial_checkpoint")), "seconds": elapsed,
    }, records


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quiet", action="store_true",
                        help="omit the per-branch progress table")
    parser.add_argument("--emit-artifact", nargs="?", const="audit_inf_kills.json",
                        default=None, metavar="PATH",
                        help="write per-branch infinity-layer audit verdicts as "
                             "JSON (for the proof-DAG join); not written if the "
                             "audit itself disagrees")
    args = parser.parse_args()
    started = time.perf_counter()
    try:
        monomials = load_h_monomials()
        errors, stats, records = [], [], []
        for window in WINDOWS:
            window_errors, window_stats, window_records = audit_window(
                window, monomials, args.quiet)
            errors.extend(window_errors)
            stats.append(window_stats)
            records.extend(window_records)
            if window_stats["skipped"] and not args.quiet:
                print(f"{window.name}: SKIP ({window.artifact} absent)")
        elapsed = time.perf_counter() - started
    except Exception as exc:
        print(f"AUDIT_ERROR {type(exc).__name__}: {exc}")
        return 2
    if errors:
        for error in errors:
            print("DISAGREEMENT", error)
        print(f"audit_inf_cases: FAIL; disagreements={len(errors)}; "
              f"runtime_seconds={elapsed:.3f}")
        if args.emit_artifact:
            # A failing audit emits NOTHING.  Anything else would leave a
            # joinable artifact on disk that the DAG would read as support.
            print("NOT emitting a join artifact: the audit disagreed")
        return 1
    details = []
    for item in stats:
        if item["skipped"]:
            details.append(f"{item['window']}=SKIP(absent)")
        else:
            partial = ",partial" if item["partial"] else ""
            details.append(f"{item['window']}=PASS({item['branches']} branches,"
                           f"{item['survivors']} survivor cases,"
                           f"{item['removed']} removed cases{partial},"
                           f"{item['seconds']:.3f}s)")
    print("audit_inf_cases: PASS; weighted_homogeneity=PASS; " +
          "; ".join(details) + f"; runtime_seconds={elapsed:.3f}")
    if args.emit_artifact:
        emit_artifact(args.emit_artifact, records, stats)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
