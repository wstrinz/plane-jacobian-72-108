#!/usr/bin/env python3
"""Local SymPy self-test for the independent f37 replay construction.

The left-hand side is rebuilt here from the mathematical formal-series
definition. Only after that construction is complete is generators.json
loaded, as a reference for exact comparison. The Macaulay2 and Sage replay
scripts do not read generators.json.
"""
from __future__ import annotations

import sympy as sp

from system_generators import load_generators


Phi = sp.Symbol("Phi")
d2, d1, d0 = sp.symbols("d2 d1 d0")
dm = {k: sp.Symbol(f"dm{k}") for k in range(1, 14)}


def s_coeff(n: int) -> sp.Expr:
    """Coefficient of u^n in S for the t=4 construction."""
    if n == 0:
        return sp.Integer(1)
    if n == 2:
        return d2
    if n == 3:
        return d1
    if n == 4:
        return d0
    if 5 <= n <= 17:
        return dm[n - 4]
    return sp.Integer(0)


def square_coeff(n: int) -> sp.Expr:
    return sp.expand(sum(s_coeff(i) * s_coeff(n - i) for i in range(n + 1)))


def cube_coeff(n: int) -> sp.Expr:
    return sp.expand(
        sum(
            s_coeff(i) * s_coeff(j) * s_coeff(n - i - j)
            for i in range(n + 1)
            for j in range(n - i + 1)
        )
    )


def build_from_definition() -> dict[str, sp.Expr]:
    """Apply the eight linear D^2 substitutions, then take four D^3 slices."""
    substitutions: dict[sp.Symbol, sp.Expr] = {}
    linear_steps = (
        (1, dm[5]),
        (2, dm[6]),
        (3, dm[7]),
        (4, dm[8]),
        (5, dm[9]),
        (6, dm[10]),
        (7, dm[11]),
        (9, dm[13]),
    )
    for k, target in linear_steps:
        equation = sp.expand(square_coeff(8 + k).subs(substitutions))
        coefficient = sp.diff(equation, target)
        assert coefficient != 0 and not coefficient.has(target)
        assert sp.diff(equation, target, 2) == 0
        substitutions[target] = sp.cancel(
            -equation.subs(target, 0) / coefficient
        )

    G1, G2, G3, G5body = (
        sp.expand(cube_coeff(12 + j).subs(substitutions))
        for j in (1, 2, 3, 5)
    )
    sol4 = sp.cancel(-G1.subs(dm[4], 0) / sp.diff(G1, dm[4]))

    # These are the cleared numerators produced after substituting sol4.
    H2 = sp.expand(2 * (dm[1] * G2 - dm[2] * G1))
    H3 = sp.expand(2 * (dm[1] * G3 - dm[3] * G1))
    return {
        "G1": G1,
        "G2": G2,
        "G3": G3,
        "G5body": G5body,
        "H2": H2,
        "H3": H3,
        "sol4": sol4,
    }


def check_t3_template() -> bool:
    """The published t=3 checkpoint used to validate the indexing convention."""

    def c3(n: int) -> sp.Expr:
        if n == 0:
            return sp.Integer(1)
        if n == 2:
            return d1
        if n == 3:
            return d0
        if 4 <= n <= 13:
            return dm[n - 3]
        return sp.Integer(0)

    coeff_u7 = sp.expand(sum(c3(i) * c3(7 - i) for i in range(8)))
    expected = 2 * d0 * dm[1] + 2 * d1 * dm[2] + 2 * dm[4]
    return sp.expand(coeff_u7 - expected) == 0


def main() -> None:
    ours = build_from_definition()
    reference = load_generators()

    checks: list[tuple[str, bool]] = [
        ("published t=3 indexing checkpoint", check_t3_template())
    ]
    for name in ("G1", "G2", "G3", "G5body", "H2", "H3"):
        checks.append(
            (
                f"independent construction equals generators.json: {name}",
                sp.expand(ours[name] - reference[name]) == 0,
            )
        )
    checks.append(
        (
            "independent construction equals generators.json: sol4",
            sp.cancel(ours["sol4"] - reference["sol4"]) == 0,
        )
    )

    for label, ok in checks:
        print(("PASS" if ok else "FAIL") + f": {label}")
    if not all(ok for _, ok in checks):
        raise SystemExit(1)
    print("PASS: all 8 replay-construction self-tests")


if __name__ == "__main__":
    main()
