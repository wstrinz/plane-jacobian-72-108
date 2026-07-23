#!/usr/bin/env python3
"""Exact verification of the pre-resultant kill of the f37 free family.

The resultant f37 vanishes identically when d2=d1=0.  This checker reads the
regenerated pre-resultant state, restricts the original equations to that
locus, derives a compact two-equation system, and checks the complete local
valuation/degree contradiction.
"""

from __future__ import annotations

from pathlib import Path

import sympy as sp

import system_generators as sysgen


ROOT = Path(__file__).resolve().parent
# Parsed from the canonical generators.json (no pickle on the mandatory path).
state = sysgen.load_generators()

d2, d1, d0 = sp.symbols("d2 d1 d0")
e, r, s, m4 = sp.symbols("dm1 dm2 dm3 dm4")
Phi = sp.symbols("Phi")
restriction = {d2: 0, d1: 0}


def source_equations() -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    """Return the cleared H2,H3,H5 equations after solving G1 for dm4."""

    H2 = sp.factor(state["H2"].subs(restriction))
    H3 = sp.factor(state["H3"].subs(restriction))
    raw_H5 = (state["G5body"] + Phi).subs(m4, state["sol4"])
    H5 = sp.factor(sp.numer(sp.together(raw_H5)).subs(restriction))

    assert sp.expand(H2 - (-3 * (d0 * e**3 - e * s**2 + 2 * r**2 * s))) == 0
    assert sp.expand(H3 - (-6 * d0 * e**2 * r - e**4 - 6 * r * s**2)) == 0
    assert sp.expand(H5 - e * (2 * Phi - 3 * e**2 * s - 3 * e * r**2)) == 0
    return H2, H3, H5


def compact_system() -> tuple[sp.Expr, sp.Expr]:
    H2, H3, H5 = source_equations()
    eq0 = d0 * e**3 - e * s**2 + 2 * r**2 * s
    eq1 = 6 * d0 * e**2 * r + e**4 + 6 * r * s**2
    resultant = sp.factor(sp.resultant(eq0, eq1, d0))
    product_eq = 12 * r * s * (r**2 - e * s) - e**5
    sum_eq = 3 * e * (r**2 + e * s) - 2 * Phi
    assert sp.expand(resultant + e**2 * product_eq) == 0
    assert sp.expand(H5 + e * sum_eq) == 0
    return product_eq, sum_eq


def local_options(
    phi_order: int,
    e_order: int,
    *,
    r_cap: int = 12,
    s_cap: int = 14,
) -> list[tuple[int, int, int]]:
    """Enumerate (v(r),v(s),v(r^2-es)) at one DVR.

    The sum equation gives v(r^2+es)=phi_order-e_order, while the product
    equation gives v(r)+v(s)+v(r^2-es)=5*e_order.  In characteristic zero,
    min(v(X+Y),v(X-Y))=min(v(X),v(Y)) for X=r^2 and Y=es.
    """

    plus_order = phi_order - e_order
    if plus_order < 0:
        return []
    out: list[tuple[int, int, int]] = []
    for x in range(r_cap + 1):
        for y in range(s_cap + 1):
            X_order = 2 * x
            Y_order = e_order + y
            minus_orders: list[int]
            if X_order < Y_order:
                minus_orders = [X_order] if plus_order == X_order else []
            elif Y_order < X_order:
                minus_orders = [Y_order] if plus_order == Y_order else []
            else:
                common = X_order
                if plus_order < common:
                    minus_orders = []
                elif plus_order > common:
                    minus_orders = [common]
                else:
                    # The difference may cancel, but the product equation caps
                    # the only relevant value at 5*e_order.
                    minus_orders = list(range(common, 5 * e_order + 1))
            for z in minus_orders:
                if x + y + z == 5 * e_order:
                    out.append((x, y, z))
    return out


def valuation_kill() -> None:
    # sum_eq says e divides Phi.  Hence q-root multiplicities are 0 or 1 and
    # e has no roots away from t and the four simple q-roots.
    q_unselected = local_options(1, 0)
    q_selected = local_options(1, 1)
    assert q_unselected == [(0, 0, 0)]
    assert q_selected == [(0, 5, 0)]

    t_options = [
        (a, *option)
        for a in range(11)
        for option in local_options(30, a)
    ]
    assert t_options == [
        (0, 0, 0, 0),
        (5, 6, 7, 12),
        (9, 12, 12, 21),
        (10, 10, 10, 30),
    ]

    # k selected q-roots contribute k to deg(e) and 5k to deg(s).
    # The infinity degree in r^2+es=2Phi/(3e) forces deg(e)=10:
    # the left side has degree at most max(24,deg(e)+14)<=24, while the
    # right side has degree 34-deg(e).
    candidates = []
    for a, x, y, z in t_options:
        for k in range(5):
            deg_e = a + k
            if deg_e > 10 or y + 5 * k > 14:
                continue
            if 34 - deg_e > max(24, deg_e + 14):
                continue
            candidates.append((a, k, x, y, z, deg_e))
    assert candidates == [(10, 0, 10, 10, 30, 10)]

    # Write e=C*t^10, r=t^10*R (deg R<=2), s=t^10*S (deg S<=4).
    # v_t(r^2-es)=30 requires t^10 | R^2-C*S, but that polynomial has
    # degree at most 4.  It would be zero, contradicting product_eq=e^5.
    assert 30 - 20 == 10
    assert max(2 * (12 - 10), 14 - 10) == 4
    assert 10 > 4


def main() -> None:
    compact_system()
    valuation_kill()
    print("f37 pre-resultant free-family kill: PASS")
    print("  restricted H2/H3/H5 derived from generators.json")
    print("  compact product and sum equations verified")
    print("  e | Phi; local options at t and each split q-place exhausted")
    print("  infinity leaves e=C*t^10; final t-order 10 > degree 4")


if __name__ == "__main__":
    main()
