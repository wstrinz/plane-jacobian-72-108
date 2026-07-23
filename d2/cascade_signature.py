#!/usr/bin/env python3
"""Extract the exact f31 lower-cascade signature from ``f31_graded.txt``.

The module deliberately keeps source parsing separate from valuation logic.
It supplies the monomial data that a tropical transition engine consumes,
and proves that rewriting with sigma = 4*d0 - d2**2 is reversible.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import re

import sympy as sp


ROOT = Path(__file__).resolve().parent
DEFAULT_SOURCE = ROOT / "f31_graded.txt"

d2, d1, d0, e, sigma = sp.symbols("d2 d1 d0 e sigma")
SOURCE_VARIABLES = (d2, d1, d0, e)
SIGNATURE_VARIABLES = (d2, d1, sigma, e)
SOURCE_WEIGHTS = {d2: 2, d1: 3, d0: 4, e: 5}
SIGNATURE_WEIGHTS = {d2: 2, d1: 3, sigma: 4, e: 5}

LINE_PATTERN = re.compile(
    r"h_(\d+)\s*\(weight\s+(\d+),\s*dm1-power\s+(\d+)\)\s*=\s*(.+)$"
)


@dataclass(frozen=True)
class CascadeLevel:
    index: int
    weight: int
    dm1_power: int
    source_expression: sp.Expr
    sigma_expression: sp.Expr

    @property
    def degree_cap(self) -> int:
        """Subcase-(2) stripped-window cap: twice the weighted degree."""

        return 2 * self.weight

    def monomial_records(self) -> list[dict[str, object]]:
        records = []
        for exponents, coefficient in sp.Poly(
            self.sigma_expression, *SIGNATURE_VARIABLES
        ).terms():
            records.append(
                {
                    "exponents": dict(
                        zip(("d2", "d1", "sigma", "e"), exponents)
                    ),
                    "coefficient": str(coefficient),
                }
            )
        return records


def rewrite_in_sigma(expression: sp.Expr) -> sp.Expr:
    """Eliminate d0 exactly using d0=(sigma+d2**2)/4."""

    return sp.expand(expression.subs(d0, (sigma + d2**2) / 4))


def load_levels(path: Path = DEFAULT_SOURCE) -> dict[int, CascadeLevel]:
    levels: dict[int, CascadeLevel] = {}
    local_symbols = {"d2": d2, "d1": d1, "d0": d0, "dm1": e}
    for line_number, line in enumerate(
        path.read_text(encoding="utf-8").splitlines(), start=1
    ):
        match = LINE_PATTERN.match(line.strip())
        if not match:
            continue
        index, weight, dm1_power = map(int, match.group(1, 2, 3))
        if index in levels:
            raise ValueError(f"duplicate h_{index} at {path}:{line_number}")
        source_expression = sp.sympify(match.group(4), locals=local_symbols)
        levels[index] = CascadeLevel(
            index=index,
            weight=weight,
            dm1_power=dm1_power,
            source_expression=source_expression,
            sigma_expression=rewrite_in_sigma(source_expression),
        )

    if sorted(levels) != list(range(8)):
        raise ValueError(f"expected h_0,...,h_7 in {path}; found {sorted(levels)}")
    validate_levels(levels)
    return levels


def _check_weight(expression: sp.Expr, variables, weights, expected: int) -> None:
    for monomial in sp.Poly(expression, *variables).monoms():
        actual = sum(weights[var] * exponent for var, exponent in zip(variables, monomial))
        if actual != expected:
            raise ValueError(
                f"nonhomogeneous monomial {monomial}: weight {actual}, expected {expected}"
            )


def validate_levels(levels: dict[int, CascadeLevel]) -> None:
    """Check source metadata, homogeneity, and the sigma round trip."""

    for index, level in levels.items():
        if level.weight != 20 - 2 * index:
            raise ValueError(f"h_{index}: wrong weight metadata")
        if level.dm1_power != 21 - 3 * index:
            raise ValueError(f"h_{index}: wrong dm1-power metadata")
        _check_weight(
            level.source_expression,
            SOURCE_VARIABLES,
            SOURCE_WEIGHTS,
            level.weight,
        )
        _check_weight(
            level.sigma_expression,
            SIGNATURE_VARIABLES,
            SIGNATURE_WEIGHTS,
            level.weight,
        )
        round_trip = level.sigma_expression.subs(sigma, 4 * d0 - d2**2)
        if sp.expand(round_trip - level.source_expression) != 0:
            raise ValueError(f"h_{index}: sigma rewrite is not reversible")


def build_signature(path: Path = DEFAULT_SOURCE) -> dict[str, object]:
    levels = load_levels(path)
    return {
        "component": "f31",
        "source": path.name,
        "variables": ["d2", "d1", "sigma", "e"],
        "weights": {str(var): weight for var, weight in SIGNATURE_WEIGHTS.items()},
        "sigma_definition": "sigma = 4*d0 - d2**2",
        "levels": [
            {
                "index": level.index,
                "weight": level.weight,
                "dm1_power": level.dm1_power,
                "degree_cap": level.degree_cap,
                "term_count": len(level.monomial_records()),
                "expression": str(level.sigma_expression),
                "monomials": level.monomial_records(),
            }
            for level in (levels[index] for index in range(8))
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", action="store_true", help="emit machine-readable JSON")
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    args = parser.parse_args()
    signature = build_signature(args.source)
    if args.json:
        print(json.dumps(signature, indent=2, sort_keys=True))
        return

    print(f"{signature['component']} cascade signature from {signature['source']}")
    print("sigma = 4*d0 - d2**2; variables=(d2,d1,sigma,e)")
    for level in signature["levels"]:
        print(
            f"h_{level['index']}: weight={level['weight']}, "
            f"degree_cap={level['degree_cap']}, terms={level['term_count']}"
        )
        print(f"  {level['expression']}")


if __name__ == "__main__":
    main()
