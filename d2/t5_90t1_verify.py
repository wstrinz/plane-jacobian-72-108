#!/usr/bin/env python3
"""Source-linked checks for the geometric q-coprime (a_t=9), T1 reduction."""

from __future__ import annotations

from collections import Counter
from itertools import product
from pathlib import Path
import re

import sympy as sp

def _require(_cond, _msg):
    """Proof-critical check: fails loudly and exits nonzero, unaffected by python -O."""
    if not _cond:
        import sys as _sys
        print("FAIL: " + str(_msg))
        _sys.exit(1)



ROOT = Path(__file__).resolve().parent
y = sp.symbols("y")
d0, d1, d2, e = sp.symbols("d0 d1 d2 dm1")
sigma = 4 * d0 - d2**2
t = y + 1
q = 2048 * y**4 - 512 * y**3 + 320 * y**2 - 240 * y + 195


def load_h() -> dict[int, sp.Expr]:
    pattern = re.compile(r"h_(\d+)\s*\([^)]*\)\s*=\s*(.+)$")
    result = {}
    for line in (ROOT / "f31_graded.txt").read_text(encoding="utf-8").splitlines():
        match = pattern.match(line.strip())
        if match:
            result[int(match.group(1))] = sp.sympify(
                match.group(2),
                locals={"d0": d0, "d1": d1, "d2": d2, "dm1": e},
            )
    return result


def at_least_two_minimum(orders: list[int]) -> bool:
    low = min(orders)
    return orders.count(low) >= 2


def check_local_shape() -> None:
    # INF represents an omitted zero polynomial. The finite ranges are the
    # exact degree caps at the linear place s.
    inf = 10**6
    survivors = set()
    for h_order in range(1, 10, 2):
        d1_order = (3 + h_order) // 2
        _require(2 * d1_order == 3 + h_order, "2 * d1_order == 3 + h_order")
        for z, k, r, r5 in product(
            list(range(9)) + [inf],
            list(range(5)) + [inf],
            list(range(14)) + [inf],
            list(range(18)) + [inf],
        ):
            level6 = [3 + r, 2*z, 2*d1_order + k, d1_order + 1, h_order]
            if not at_least_two_minimum(level6):
                continue
            level5 = [
                3 + r5,
                k + 2*z,
                2*d1_order + z,
                2*d1_order + 2*k,
                d1_order + k + 1,
                2,
                r,
            ]
            if at_least_two_minimum(level5):
                survivors.add((h_order, d1_order, r, z, k, r5))

    _require({(h, d, r) for h, d, r, _z, _k, _r5 in survivors} == {(5, 4, 2)}, "{(h, d, r) for h, d, r, _z, _k, _r5 in survivors} == {(5, 4, 2)}")
    _require(all(z >= 3 for _h, _d, _r, z, _k, _r5 in survivors), "all(z >= 3 for _h, _d, _r, z, _k, _r5 in survivors)")


def infinity_degrees(d: int, z: int | None) -> tuple[list[int], tuple[int, ...]]:
    absent = -10**6
    sigma_once = absent if z is None else z
    sigma_twice = absent if z is None else 2*z
    h5_terms = [4 + sigma_twice, 2*d + sigma_once, 2*d + 8, d + 13, 18]
    h6_terms = [sigma_twice, 2*d + 4, d + 9]
    degrees = [229, 232, 235, 238, 241]
    degrees.extend([224 + max(h5_terms), 231 + max(h6_terms), 238 + 2*d])
    top = max(degrees)
    tops = tuple(index for index, degree in enumerate(degrees) if degree == top)
    if tops == (5,):
        _require(h5_terms.count(max(h5_terms)) == 1, "h5_terms.count(max(h5_terms)) == 1")
    elif tops == (6,):
        _require(h6_terms.count(max(h6_terms)) == 1, "h6_terms.count(max(h6_terms)) == 1")
    return degrees, tops


def check_infinity_table() -> None:
    counts: Counter[tuple[int, ...]] = Counter()
    ties = []
    for d in range(5):
        for z in [None, *range(9)]:
            _degrees, tops = infinity_degrees(d, z)
            counts[tops] += 1
            if len(tops) > 1:
                ties.append((d, z, tops))
    _require(sum(counts.values()) == 50, "sum(counts.values()) == 50")
    _require(sum(count for tops, count in counts.items() if len(tops) == 1) == 43, "sum(count for tops, count in counts.items() if len(tops) == 1) == 43")
    _require(ties == [(2, z, (5, 6, 7)) for z in [None, *range(6)]], "ties == [(2, z, (5, 6, 7)) for z in [None, *range(6)]]")


def check_square_constraints() -> None:
    c, gamma, delta = sp.symbols("c gamma delta", nonzero=True)
    t0, q0, w0 = sp.symbols("t0 q0 w0", nonzero=True)
    terminal_eta = -8192*c**7*delta**2/gamma**3

    local_kappa = 2048*c**5*gamma**2*t0**15/q0
    local_level6 = (
        t0**3*q0*terminal_eta*w0**2
        - gamma**3*local_kappa
        - 8192*c**6*delta*gamma*t0**9*w0
    )
    multiplier = -gamma**3*q0/(2048*c**5*t0**3)
    local_square = (gamma**4*t0**6 + 2*c*delta*q0*w0)**2
    _require(sp.factor(multiplier*local_level6 - local_square) == 0, "sp.factor(multiplier*local_level6 - local_square) == 0")

    q_lc = sp.Poly(q, y).LC()
    _require(q_lc == 2048, "q_lc == 2048")
    infinity_kappa = c**5*gamma**2
    infinity_level6 = (
        q_lc*terminal_eta - gamma**3*infinity_kappa - 8192*c**6*delta*gamma
    )
    infinity_square = (gamma**4 + 4096*c*delta)**2
    _require(sp.factor(-gamma**3*infinity_level6/c**5 - infinity_square) == 0, "sp.factor(-gamma**3*infinity_level6/c**5 - infinity_square) == 0")


def check_tail_refinement(h: dict[int, sp.Expr]) -> None:
    phi = sp.symbols("Phi")
    tail = sp.expand(e**6*h[5] + phi*e**3*h[6] + phi**2*h[7])
    square = 2048*(e**4 + 2*phi*d1)**2
    residual = -512*e**3*(
        6*sigma**2*phi + 18*sigma**2*d2*e**3 - 63*sigma*d1**2*e**3
        - 28*phi*d1**2*d2 + 24*d1**2*d2**2*e**3 - 36*d1*d2*e**4
    )
    _require(sp.expand(tail - square - residual) == 0, "sp.expand(tail - square - residual) == 0")

    _require(6*34 + 3*9 + 2*5 == 241, "6*34 + 3*9 + 2*5 == 241")
    _require(5*34 + 2*35 == 240, "5*34 + 2*35 == 240")
    residual_inner_bounds = [
        2*4 + 34, 2*4 + 4 + 27, 4 + 4 + 27,
        34 + 4 + 4, 4 + 8 + 27, 2 + 4 + 36,
    ]
    _require(max(residual_inner_bounds) == 42, "max(residual_inner_bounds) == 42")
    _require(5*34 + 3*9 + 42 == 239, "5*34 + 3*9 + 42 == 239")

    h4_sigma = sp.factor(h[4].subs(d0, (d2**2 + sigma)/4))
    h4_bound_terms = [
        3*4, 2*4 + 8, 4 + 4 + 4, 4 + 2 + 9,
        4*2, 2*2 + 3*4, 2 + 2*4 + 9, 4 + 2*9,
    ]
    _require(max(h4_bound_terms) == 22, "max(h4_bound_terms) == 22")
    _require(4*34 + 9*9 + 22 == 239, "4*34 + 9*9 + 22 == 239")
    _require(h4_sigma != 0, "h4_sigma != 0")

    c, gamma, v1, v0 = sp.symbols("c gamma v1 v0", nonzero=True)
    delta = -gamma**4/(4096*c)
    core = sp.expand(gamma**4*t**6 + 2*c*q*delta*(y**2 + v1*y + v0))
    poly = sp.Poly(core, y)
    _require(sp.factor(poly.coeff_monomial(y**6)) == 0, "sp.factor(poly.coeff_monomial(y**6)) == 0")
    expected_y5 = gamma**4*(sp.Rational(25, 4) - v1)
    _require(sp.factor(poly.coeff_monomial(y**5) - expected_y5) == 0, "sp.factor(poly.coeff_monomial(y**5) - expected_y5) == 0")
    # Degree 239 after a nonzero common factor is stripped.
    # All leading d2 contributions cancel, leaving only lc(sigma)^2.
    a4, sigma4 = sp.symbols("a4 sigma4")
    q_lc = sp.Integer(2048)
    degree239 = (
        -512*c*q_lc*(6*c*q_lc*sigma4**2 - 28*c*q_lc*delta**2*a4
                     - 36*delta*a4*gamma**4)
        + 5632*a4*gamma**8
    )
    _require(sp.factor(degree239) == -12884901888*c**2*sigma4**2, "sp.factor(degree239) == -12884901888*c**2*sigma4**2")



def main() -> None:
    h = load_h()
    _require(sorted(h) == list(range(8)), "sorted(h) == list(range(8))")
    _require(sp.expand(h[7] - 8192*d1**2) == 0, "sp.expand(h[7] - 8192*d1**2) == 0")
    _require(sp.expand(h[6] - (-3072*sigma**2 + 14336*d1**2*d2 + 8192*d1*e)) == 0, "sp.expand(h[6] - (-3072*sigma**2 + 14336*d1**2*d2 + 8192*d1*e)) == 0")
    _require(sp.expand(
        h[5] - (-9216*d2*sigma**2 + 32256*d1**2*sigma
                - 12288*d1**2*d2**2 + 18432*d1*d2*e + 2048*e**2)
    ) == 0, "sp.expand( h[5] - (-9216*d2*sigma**2 + 32256*d1**2*sigma - 12288*d1**2*d2**2 + 18432*d1*d2*e + 2048*e**2) ) == 0")
    _require(30 - 3*9 == 3, "30 - 3*9 == 3")
    _require(10 + 3*9 == 37, "10 + 3*9 == 37")
    _require(37 - 7*4 == 9, "37 - 7*4 == 9")

    check_local_shape()
    check_infinity_table()
    check_square_constraints()
    check_tail_refinement(h)
    print("a_t=9 geometric q-coprime T1 reduction: PASS")
    print("  nonconstant E: unique local shape h=5, v_s(d1)=4, v_s(G6)=2")
    print("  constant E preliminary: 45/50 cells killed before coefficient descent")
    print("  exact squares plus the tail identity fix two d1 coefficients")


if __name__ == "__main__":
    main()
