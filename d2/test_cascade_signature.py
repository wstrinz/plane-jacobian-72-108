#!/usr/bin/env python3
"""Independent exact regression checks for cascade_signature.py."""

from __future__ import annotations

import sympy as sp

import cascade_signature as signature


def main() -> None:
    levels = signature.load_levels()

    assert sorted(levels) == list(range(8))
    assert [levels[index].degree_cap for index in range(8)] == [
        40,
        36,
        32,
        28,
        24,
        20,
        16,
        12,
    ]

    d2, d1, sigma, e = signature.SIGNATURE_VARIABLES
    expected = {
        7: 8192 * d1**2,
        6: -3072 * sigma**2 + 14336 * d1**2 * d2 + 8192 * d1 * e,
        5: (
            -9216 * d2 * sigma**2
            + 32256 * d1**2 * sigma
            - 12288 * d1**2 * d2**2
            + 18432 * d1 * d2 * e
            + 2048 * e**2
        ),
        4: -16
        * (
            756 * sigma**3
            + 324 * sigma**2 * d2**2
            + 1476 * sigma * d1**2 * d2
            - 2160 * sigma * d1 * e
            + 13797 * d1**4
            + 1952 * d1**2 * d2**3
            + 192 * d1 * d2**2 * e
            - 352 * d2 * e**2
        ),
    }
    for index, expression in expected.items():
        assert sp.expand(levels[index].sigma_expression - expression) == 0, index

    encoded = signature.build_signature()
    assert encoded["source"] == "f31_graded.txt"
    assert [level["term_count"] for level in encoded["levels"][-4:]] == [8, 5, 3, 1]
    assert all(level["monomials"] for level in encoded["levels"])

    print("cascade signature: PASS")
    print("  source-linked h_0,...,h_7 parse and sigma round trip")
    print("  exact h_7/h_6/h_5/h_4 terminal ladder")
    print("  machine-readable monomial tables and subcase-(2) degree caps")


if __name__ == "__main__":
    main()
