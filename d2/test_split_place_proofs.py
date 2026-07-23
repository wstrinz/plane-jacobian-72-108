#!/usr/bin/env python3
"""Source-linked checks for FIELD_SPLIT_AUDIT.md.

Unlike the compact audit-bundle checker, this file parses f31_graded.txt and
checks that the formulas and degree bounds used by the a_t=7 proof are the
actual repository formulas.
"""

from __future__ import annotations

from itertools import product
from pathlib import Path
import re

import sympy as sp


ROOT = Path(__file__).resolve().parent
d0, d1, d2, e = sp.symbols("d0 d1 d2 dm1")
sigma = 4 * d0 - d2**2


def load_h() -> dict[int, sp.Expr]:
    out: dict[int, sp.Expr] = {}
    pattern = re.compile(r"h_(\d+)\s*\([^)]*\)\s*=\s*(.+)$")
    for line in (ROOT / "f31_graded.txt").read_text(encoding="utf-8").splitlines():
        match = pattern.match(line.strip())
        if match:
            out[int(match.group(1))] = sp.sympify(
                match.group(2),
                locals={"d0": d0, "d1": d1, "d2": d2, "dm1": e},
            )
    assert sorted(out) == list(range(8))
    return out


H = load_h()


def weighted_degree_bound(expr: sp.Expr, degrees: dict[sp.Symbol, int]) -> int:
    poly = sp.Poly(expr, d2, d1, d0, e)
    return max(
        sum(power * degrees[var] for power, var in zip(monom, (d2, d1, d0, e)))
        for monom, _coeff in poly.terms()
    )


def check_source_formulas() -> None:
    h6_expected = (
        -3072 * sigma**2 + 14336 * d1**2 * d2 + 8192 * d1 * e
    )
    h5_expected = (
        -9216 * d2 * sigma**2
        + 32256 * d1**2 * sigma
        - 12288 * d1**2 * d2**2
        + 18432 * d1 * d2 * e
        + 2048 * e**2
    )
    assert sp.expand(H[7] - 8192 * d1**2) == 0
    assert sp.expand(H[6] - h6_expected) == 0
    assert sp.expand(H[5] - h5_expected) == 0
    assert sp.expand(H[6].subs(d1, 0) + 3072 * sigma**2) == 0
    assert sp.expand(
        H[5].subs(d1, 0) - (-9216 * d2 * sigma**2 + 2048 * e**2)
    ) == 0

    generic_caps = {d2: 4, d1: 6, d0: 8, e: 10}
    for f, expr in H.items():
        assert weighted_degree_bound(expr, generic_caps) <= 40 - 4 * f


def local_parity_escape(include_sigma: bool) -> list[tuple[int, int, int]]:
    escapes = []
    for e_order in (1, 3):
        for h_order in (1, 3):
            d1_order = (3 * e_order + h_order) // 2
            assert 2 * d1_order == 3 * e_order + h_order
            sigma_orders = range(9) if include_sigma else (None,)
            possible = False
            for sigma_order in sigma_orders:
                for g6_order in range(8):
                    orders = [
                        3 * e_order + g6_order,
                        2 * d1_order,
                        d1_order + e_order,
                        h_order,
                    ]
                    if sigma_order is not None:
                        orders.append(2 * sigma_order)
                    low = min(orders)
                    if sum(order == low for order in orders) >= 2:
                        possible = True
                        break
                if possible:
                    break
            if possible:
                escapes.append((e_order, h_order, d1_order))
    return escapes


def term_bounds(
    deg_u: int, deg_v: int, deg_sigma: int
) -> tuple[list[int], int]:
    deg_e = 7 + 2 * deg_u
    deg_d1 = 3 * deg_u + deg_v
    bounds = [34 * f + (21 - 3 * f) * deg_e + (40 - 4 * f) for f in range(6)]

    if deg_u == 1 and deg_v == 0 and deg_sigma <= 6:
        h5_terms = (
            4 + 2 * deg_sigma,
            2 * deg_d1 + deg_sigma,
            2 * deg_d1 + 8,
            deg_d1 + 4 + deg_e,
            2 * deg_e,
        )
        assert max(h5_terms) == 18
        bounds[5] = 5 * 34 + 6 * deg_e + 18

    h6_terms = (
        2 * deg_sigma,
        2 * deg_d1 + 4,
        deg_d1 + deg_e,
    )
    h6_bound = max(h6_terms)
    bounds.append(6 * 34 + 3 * deg_e + h6_bound)
    bounds.append(7 * 34 + 2 * deg_d1)
    top = max(bounds)
    top_index = bounds.index(top)

    if top_index == 6:
        # The sigma^2 term must strictly lead h6, so the asserted top degree
        # cannot disappear through cancellation.
        assert h6_terms[0] > max(h6_terms[1:])
    else:
        assert top_index == 7
        assert H[7] == 8192 * d1**2
    return bounds, top_index


def check_a7_degree_table() -> None:
    assert local_parity_escape(include_sigma=True) == [(1, 3, 3)]
    # sigma=0 is an omitted term, not a finite valuation; it has the same sole
    # escape and is checked separately to close that logical edge case.
    assert local_parity_escape(include_sigma=False) == [(1, 3, 3)]

    for deg_u, deg_v, deg_sigma in product((0, 1), (0, 1), range(9)):
        bounds, top_index = term_bounds(deg_u, deg_v, deg_sigma)
        assert bounds.count(max(bounds)) == 1
        assert top_index in (6, 7)


def main() -> None:
    check_source_formulas()
    check_a7_degree_table()
    print("source-linked split-place proof checks: PASS")
    print("  parsed h_0,...,h_7 from f31_graded.txt")
    print("  h7/h6/h5 formulas and all generic degree caps verified")
    print("  a=7 parity escape checked with sigma nonzero and sigma=0")
    print("  all 36 infinity cases have a noncancelling unique top term")


if __name__ == "__main__":
    main()
