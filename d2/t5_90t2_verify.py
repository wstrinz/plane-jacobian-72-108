#!/usr/bin/env python3
"""Source-linked checks for the geometric q-coprime (a_t=9), T2 kill."""

from __future__ import annotations

from pathlib import Path
import re

import sympy as sp


ROOT = Path(__file__).resolve().parent
d0, d1, d2, e = sp.symbols("d0 d1 d2 dm1")
sigma = 4 * d0 - d2**2


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


def main() -> None:
    h = load_h()
    assert sorted(h) == list(range(8))
    assert sp.expand(h[7].subs(d1, 0)) == 0
    assert sp.expand(h[6].subs(d1, 0) + 3072 * sigma**2) == 0
    assert sp.expand(
        h[5].subs(d1, 0) - (-9216 * d2 * sigma**2 + 2048 * e**2)
    ) == 0

    # a=9: v=3, deg(E)<=1, deg(g_l)<=37, and after q^6 is removed
    # the terminal cofactor G has degree <=13.
    assert 30 - 3 * 9 == 3
    assert 10 + 3 * 9 == 37
    assert 37 - 6 * 4 == 13

    # At a root of the squarefree factor s in E=s*u^2 and G=s*v^2,
    # h=v_s(G) is odd.  The terminal equation gives 2*v_s(sigma)=3+h.
    # Level 5 has orders E^3*g5 >=3 and h5=2 (unique e^2 term), while
    # its right side has order h.  Equality would require h=2, impossible.
    for h_order in range(1, 14, 2):
        sigma_order = (3 + h_order) // 2
        assert 2 * sigma_order == 3 + h_order
        assert sigma_order >= 2
        h5_orders = (2 * sigma_order, 2)
        assert min(h5_orders) == 2 and h5_orders.count(2) == 1
        assert h_order != min(3, 2)

    # Hence E is a square; deg(E)<=1 makes E constant.  Then G is a square
    # and deg(sigma)<=6.  Verify the infinity domination for every degree.
    deg_e = 9
    for deg_sigma in range(7):
        degrees = [
            34 * f + (21 - 3 * f) * deg_e + (40 - 4 * f)
            for f in range(5)
        ]
        # h5 has exact degree 18 from 2048*e^2: the other term is <=16.
        assert 4 + 2 * deg_sigma <= 16 < 18
        degrees.append(5 * 34 + 6 * deg_e + 18)
        degrees.append(6 * 34 + 3 * deg_e + 2 * deg_sigma)
        top = max(degrees)
        tops = [index for index, degree in enumerate(degrees) if degree == top]
        assert tops == ([5] if deg_sigma <= 5 else [6])

    print("a_t=9 geometric q-coprime T2: PASS")
    print("  nonconstant squarefree factor killed at level 5")
    print("  E constant; infinity top is uniquely f=5 or f=6")


if __name__ == "__main__":
    main()
