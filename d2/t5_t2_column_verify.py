#!/usr/bin/env python3
"""Source-linked exact checks for T5_T2_COLUMN.md."""

from __future__ import annotations

import json
import re
from pathlib import Path

import sympy as sp

def _require(_cond, _msg):
    """Proof-critical check: fails loudly and exits nonzero, unaffected by python -O."""
    if not _cond:
        import sys as _sys
        print("FAIL: " + str(_msg))
        _sys.exit(1)


ROOT = Path(__file__).resolve().parent
y = sp.Symbol("y")
d0, d1, d2, e = sp.symbols("d0 d1 d2 dm1")
sigma = 4 * d0 - d2**2
t = y + 1
q = 2048 * y**4 - 512 * y**3 + 320 * y**2 - 240 * y + 195
c = sp.Rational(-1, 6630)
Phi = c * t**30 * q


def load_h() -> dict[int, sp.Expr]:
    """Parse h_l; no h_l coefficient is copied into this verifier."""
    pattern = re.compile(r"h_(\d+)\s*\([^)]*\)\s*=\s*(.+)$")
    result = {}
    for line in (ROOT / "f31_graded.txt").read_text(encoding="utf-8").splitlines():
        match = pattern.match(line.strip())
        if match:
            result[int(match.group(1))] = sp.sympify(
                match.group(2),
                locals={"d0": d0, "d1": d1, "d2": d2, "dm1": e},
            )
    _require(sorted(result) == list(range(8)), "sorted(result) == list(range(8))")
    return result


TARGETS = (
    (5, (1, 0, 0, 0)),
    (6, (1, 0, 0, 0)), (6, (1, 1, 0, 0)), (6, (1, 1, 1, 0)),
    (7, (1, 0, 0, 0)), (7, (1, 1, 0, 0)),
    (7, (1, 1, 1, 0)), (7, (3, 0, 0, 0)),
    (8, (0, 0, 0, 0)), (8, (1, 0, 0, 0)), (8, (1, 1, 0, 0)),
    (9, (1, 0, 0, 0)),
)

FLAGS_2 = ((False, ()), (False, (4,)))
FLAGS_3 = FLAGS_2 + ((True, ()),)
EXPECTED_FLAGS = {
    TARGETS[0]: FLAGS_2,
    TARGETS[1]: FLAGS_2, TARGETS[2]: FLAGS_2, TARGETS[3]: FLAGS_2,
    TARGETS[4]: FLAGS_3, TARGETS[5]: FLAGS_3, TARGETS[6]: FLAGS_3,
    TARGETS[7]: FLAGS_2,
    TARGETS[8]: FLAGS_3, TARGETS[9]: FLAGS_3, TARGETS[10]: FLAGS_3,
    TARGETS[11]: ((False, ()), (False, (5,)), (False, (4,)), (True, ())),
}


def finite(value: int | str) -> int | None:
    return None if value == "inf" else int(value)


def case_flag(case: dict) -> tuple[bool, tuple[int, ...]]:
    return bool(case["d2_zero"]), tuple(case["g_zero_levels"])


def load_cells() -> dict[tuple[int, tuple[int, ...]], dict]:
    data = json.loads((ROOT / "cascade_cones_qt.json").read_text(encoding="utf-8"))
    cells = {
        (row["a_t"], tuple(row["b"])): row
        for row in data["branches"]
        if row["branch"] == "T2" and (row["a_t"], tuple(row["b"])) in TARGETS
    }
    _require(len(cells) == 12 and all(key in cells for key in TARGETS), "len(cells) == 12 and all(key in cells for key in TARGETS)")
    _require(sum(row["survivor_case_count"] for row in cells.values()) == 32, "sum(row[\"survivor_case_count\"] for row in cells.values()) == 32")
    for key, row in cells.items():
        _require(row["status"] == "survives", "row[\"status\"] == \"survives\"")
        _require(row["survivor_case_count"] == len(row["survivor_cases"]), "row[\"survivor_case_count\"] == len(row[\"survivor_cases\"])")
        _require(tuple(case_flag(case) for case in row["survivor_cases"]) == EXPECTED_FLAGS[key], "tuple(case_flag(case) for case in row[\"survivor_cases\"]) == EXPECTED_FLAGS[key]")
        _require(all(not case["sigma_zero"] for case in row["survivor_cases"]), "all(not case[\"sigma_zero\"] for case in row[\"survivor_cases\"])")
        _require(all(len(case["witness"]) == 5 for case in row["survivor_cases"]), "all(len(case[\"witness\"]) == 5 for case in row[\"survivor_cases\"])")
    return cells


def q_profile(case: dict) -> tuple[tuple[int, ...], tuple[int, ...], tuple[int, ...]]:
    places = [place for place in case["witness"] if place["place"] == "q"]
    _require(len(places) == 4, "len(places) == 4")
    return (
        tuple(int(place["b"]) for place in places),
        tuple(int(place["v_sigma"]) for place in places),
        tuple(int(place["v_g"]["6"]) for place in places),
    )


def h_raw_orders(level: int, vd2: int | None, vs: int, ve: int) -> list[int]:
    if level == 5:
        return [2 * ve] + ([] if vd2 is None else [vd2 + 2 * vs])
    _require(level == 4, "level == 4")
    return [3 * vs] + ([] if vd2 is None else [2 * vd2 + 2 * vs, vd2 + 2 * ve])


def check_tropical_line(left: int | None, outer: int | None, horders: list[int]) -> None:
    """Necessary ultrametric check, allowing a tied h_l minimum to rise."""
    hmin = min(horders)
    htied = horders.count(hmin) >= 2
    if outer is None:
        _require((left == hmin) if not htied else (left is None or left >= hmin), "(left == hmin) if not htied else (left is None or left >= hmin)")
    elif outer < hmin:
        _require(left == outer, "left == outer")
    elif outer > hmin and not htied:
        _require(left == hmin, "left == hmin")
    else:
        _require(left is None or left >= min(outer, hmin), "left is None or left >= min(outer, hmin)")


def verify_witnesses(cells: dict) -> None:
    for (a, bexpected), row in cells.items():
        v = 30 - 3 * a
        for case in row["survivor_cases"]:
            b, s, m = q_profile(case)
            _require(b == bexpected, "b == bexpected")
            _require(all(3 * bi + mi == 6 + 2 * si for bi, si, mi in zip(b, s, m)), "all(3 * bi + mi == 6 + 2 * si for bi, si, mi in zip(b, s, m))")
            tp = next(place for place in case["witness"] if place["place"] == "t")
            _require(int(tp["b"]) == a and int(tp["v_sigma"]) == 0, "int(tp[\"b\"]) == a and int(tp[\"v_sigma\"]) == 0")
            _require(int(tp["v_g"]["6"]) == 0, "int(tp[\"v_g\"][\"6\"]) == 0")
            vd2 = finite(tp["v_d2"])
            s4, s5, s6 = (finite(tp["v_g"][str(level)]) for level in (4, 5, 6))
            check_tropical_line(None if s6 is None else v + s6, s5,
                                h_raw_orders(5, vd2, 0, a))
            check_tropical_line(None if s5 is None else v + s5, s4,
                                h_raw_orders(4, vd2, 0, a))


def invariants(row: dict) -> tuple[int, int, int, int, int, int]:
    a = int(row["a_t"])
    b, s, m = q_profile(row["survivor_cases"][0])
    _require(all(q_profile(case) == (b, s, m) for case in row["survivor_cases"]), "all(q_profile(case) == (b, s, m) for case in row[\"survivor_cases\"])")
    B, S, M = sum(b), sum(s), sum(m)
    _require(all(mi >= 1 for mi in m), "all(mi >= 1 for mi in m)")
    return B, S, M, 10 - a - B, 10 + 3 * a - M, 30 - 3 * a


def terminal_states(a: int, B: int, S: int, fcap: int, gcap: int) -> tuple:
    """(deg F,deg Z,deg G,deg e,deg sigma), using F^2|G and F^3G=const*Z^2."""
    states = []
    for fdeg in range(fcap + 1):
        for zdeg in range(9 - S):
            gdeg = 2 * zdeg - 3 * fdeg
            if 2 * fdeg <= gdeg <= gcap:
                states.append((fdeg, zdeg, gdeg, a + B + fdeg, S + zdeg))
    return tuple(states)


def hcap(poly: sp.Poly, D: int, S: int, delta: int | None) -> int:
    values = []
    for kd, ks, ke in poly.monoms():
        if delta is None and kd:
            continue
        values.append((0 if delta is None else delta * kd) + S * ks + D * ke)
    return max(values)


def killed_at_infinity(H: dict[int, sp.Poly], D: int, S: int, d2zero: bool) -> bool:
    modes = (None,) if d2zero else (0, 1, 2, 3, 4)
    for delta in modes:
        degrees = [34 * level + (21 - 3 * level) * D + hcap(H[level], D, S, delta)
                   for level in range(6)]
        degrees.append(204 + 3 * D + 2 * S)
        tops = [level for level, degree in enumerate(degrees) if degree == max(degrees)]
        exact_unique = tops == [6]
        if tops == [5]:
            # H_5=-9216*d2*sigma^2+2048*e^2.
            exact_unique = delta is None or delta + 2 * S != 2 * D
        if not exact_unique:
            return False
    return True


KILLED_STATES = {
    TARGETS[0]: ((0, 0),),
    TARGETS[1]: ((0, 0), (0, 1), (1, 3)),
    TARGETS[2]: ((0, 0), (0, 1)),
    TARGETS[3]: ((0, 0),),
}
RESIDUAL = {
    TARGETS[4]: ((0, 1), (2, 5), (2, 6)),
    TARGETS[5]: ((1, 3), (1, 4)),
    TARGETS[6]: ((0, 0), (0, 1), (0, 2)),
    TARGETS[7]: ((0, 0), (0, 1)),
    TARGETS[8]: ((0, 3), (2, 5), (2, 6), (2, 7), (2, 8)),
    TARGETS[9]: ((1, 3), (1, 4), (1, 5), (1, 6)),
    TARGETS[10]: ((0, 0), (0, 1), (0, 2), (0, 3), (0, 4)),
    TARGETS[11]: tuple((0, z) for z in range(7)),
}


def main() -> None:
    h = load_h()
    cells = load_cells()

    # C1: every displayed collapse is checked against f31_graded.txt.
    _require(sp.expand(h[7].subs(d1, 0)) == 0, "sp.expand(h[7].subs(d1, 0)) == 0")
    _require(sp.expand(h[6].subs(d1, 0) + 3072 * sigma**2) == 0, "sp.expand(h[6].subs(d1, 0) + 3072 * sigma**2) == 0")
    _require(sp.expand(h[5].subs(d1, 0) - (-9216 * d2 * sigma**2 + 2048 * e**2)) == 0, "sp.expand(h[5].subs(d1, 0) - (-9216 * d2 * sigma**2 + 2048 * e**2)) == 0")
    h4 = -5184 * d2**2 * sigma**2 + 5632 * d2 * e**2 - 12096 * sigma**3
    _require(sp.expand(h[4].subs(d1, 0) - h4) == 0, "sp.expand(h[4].subs(d1, 0) - h4) == 0")
    _require(sp.expand(h[4].subs({d1: 0, d2: 0}) + 12096 * sigma.subs(d2, 0)**3) == 0, "sp.expand(h[4].subs({d1: 0, d2: 0}) + 12096 * sigma.subs(d2, 0)**3) == 0")
    _require(sp.degree(Phi, y) == 34, "sp.degree(Phi, y) == 34")
    print("C1. parsed h4--h7 collapses and deg(Phi)=34                         OK")

    # C2: generalized split-support rearrangement. Q=q*Qbar.
    qs, R, F, Qbar, G, g5, dd2, tv, ta2 = sp.symbols(
        "qs R F Qbar G g5 dd2 tv ta2", nonzero=True)
    sig2 = R**3 * F**3 * Qbar * G / (3072 * c**6 * qs**5)
    lhs = tv * qs * Qbar * G
    rhs = R**3 * F**3 * g5 + c**5 * qs**5 * (
        -9216 * dd2 * sig2 + 2048 * ta2 * R**2 * F**2)
    wanted = (tv * qs * Qbar * G - 2048 * c**5 * qs**5 * ta2 * R**2 * F**2
              - R**3 * F**3 * (g5 + 19890 * dd2 * Qbar * G))
    _require(sp.cancel(lhs - rhs - wanted) == 0, "sp.cancel(lhs - rhs - wanted) == 0")
    _require(sp.Rational(-9216, 3072) / c == 19890, "sp.Rational(-9216, 3072) / c == 19890")
    print("C2. terminal absorption and split-support level-5 rearrangement    OK")

    verify_witnesses(cells)
    print("C3. 12 cells / 32 cases: flags, q terminals, t coupling            OK")

    # C4: source-derived infinity cap engine and every verdict.
    s = sp.Symbol("s")
    Hexpr = {level: sp.expand(h[level].subs({d1: 0, d0: (d2**2 + s) / 4}))
             for level in range(8)}
    _require(Hexpr[7] == 0 and sp.expand(Hexpr[6] + 3072 * s**2) == 0, "Hexpr[7] == 0 and sp.expand(Hexpr[6] + 3072 * s**2) == 0")
    _require(sp.expand(Hexpr[5] - (-9216 * d2 * s**2 + 2048 * e**2)) == 0, "sp.expand(Hexpr[5] - (-9216 * d2 * s**2 + 2048 * e**2)) == 0")
    H = {level: sp.Poly(Hexpr[level], d2, s, e) for level in range(7)}
    killed_cases = open_cases = 0
    for key in TARGETS:
        row = cells[key]
        a = key[0]
        B, S, M, fcap, gcap, v = invariants(row)
        _require(M + gcap == 10 + 3 * a and v == 30 - 3 * a, "M + gcap == 10 + 3 * a and v == 30 - 3 * a")
        states = terminal_states(a, B, S, fcap, gcap)
        pairs = tuple((state[0], state[1]) for state in states)
        if key in KILLED_STATES:
            _require(pairs == KILLED_STATES[key], "pairs == KILLED_STATES[key]")
            for case in row["survivor_cases"]:
                _require(all(killed_at_infinity(H, D, sd, bool(case["d2_zero"]))
                           for _, _, _, D, sd in states), "all(killed_at_infinity(H, D, sd, bool(case[\"d2_zero\"])) for _, _, _, D, sd in states)")
                killed_cases += 1
        else:
            for case in row["survivor_cases"]:
                residual = tuple((fdeg, zdeg) for fdeg, zdeg, _, D, sd in states
                                 if not killed_at_infinity(H, D, sd, bool(case["d2_zero"])))
                _require(residual == RESIDUAL[key], "residual == RESIDUAL[key]")
                if tuple(case["g_zero_levels"]) == (5,):
                    # Exact g5=0 line at a=9: t^3*g6=c^5*q^5*h5.
                    tp = next(p for p in case["witness"] if p["place"] == "t")
                    _require(key == TARGETS[11] and tp["v_d2"] == 3, "key == TARGETS[11] and tp[\"v_d2\"] == 3")
                    narrowed = []
                    for fdeg, zdeg in residual:
                        D, sd = a + B + fdeg, S + zdeg
                        lhsdeg = v + M + (2 * zdeg - 3 * fdeg)
                        if sd < 8:
                            _require(max(4 + 2 * sd, 2 * D) == 20, "max(4 + 2 * sd, 2 * D) == 20")
                            _require(lhsdeg != 5 * 4 + 20, "lhsdeg != 5 * 4 + 20")
                        else:
                            _require(lhsdeg == 5 * 4 + 20, "lhsdeg == 5 * 4 + 20")
                            narrowed.append((fdeg, zdeg))
                    _require(tuple(narrowed) == ((0, 6),), "tuple(narrowed) == ((0, 6),)")
                open_cases += 1
    _require((killed_cases, open_cases) == (8, 24), "(killed_cases, open_cases) == (8, 24)")
    print("C4. squeeze states, source-derived infinity caps, all verdicts     OK")

    # C5: strict full-chain margins for all seven killed degree states.
    margins = []
    for key in KILLED_STATES:
        row = cells[key]
        a = key[0]
        B, S, _, fcap, gcap, _ = invariants(row)
        for fdeg, zdeg, _, D, sd in terminal_states(a, B, S, fcap, gcap):
            other = max(34 * level + (21 - 3 * level) * D
                        + hcap(H[level], D, sd, mode)
                        for mode in (None, 0, 1, 2, 3, 4) for level in range(6))
            exact6 = 204 + 3 * D + 2 * sd
            _require(exact6 > other, "exact6 > other")
            margins.append((key, fdeg, zdeg, other, exact6))
    _require(tuple(top - other for _, _, _, other, top in margins) == (8, 3, 5, 4, 2, 4, 1), "tuple(top - other for _, _, _, other, top in margins) == (8, 3, 5, 4, 2, 4, 1)")
    print("C5. killed-state margins:")
    for key, fdeg, zdeg, other, exact6 in margins:
        print(f"    a={key[0]} b={''.join(map(str, key[1]))} (f,z)=({fdeg},{zdeg}): "
              f"max(T0..T5)={other} < deg(T6)={exact6}")
    print("\nALL T2-COLUMN CHECKS PASS: 4 cells / 8 cases killed; 8 cells / 24 cases open")


if __name__ == "__main__":
    main()
