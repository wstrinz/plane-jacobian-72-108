#!/usr/bin/env python3
"""Exact coefficient certificate killing the constant-E part of a_t=9 T1.

Only requested y-coefficients are convolved.  This avoids expanding the full
degree-242 master polynomial while still deriving every value from the actual
repository h_0,...,h_7 formulas.
"""

from __future__ import annotations

import sympy as sp

import t5_90t1_verify as base

def _require(_cond, _msg):
    """Proof-critical check: fails loudly and exits nonzero, unaffected by python -O."""
    if not _cond:
        import sys as _sys
        print("FAIL: " + str(_msg))
        _sys.exit(1)



y = base.y
H = base.load_h()
c, gamma, v0 = sp.symbols("c gamma v0", nonzero=True)
a0, a1, a2, a3, a4 = sp.symbols("a0:5")
s0, s1, s2, s3 = sp.symbols("s0:4")


def from_expr(expr: sp.Expr) -> dict[int, sp.Expr]:
    poly = sp.Poly(sp.expand(expr), y)
    return {monomial[0]: coefficient for monomial, coefficient in poly.terms()}


def add(left: dict[int, sp.Expr], right: dict[int, sp.Expr]) -> dict[int, sp.Expr]:
    out = dict(left)
    for degree, coefficient in right.items():
        out[degree] = out.get(degree, 0) + coefficient
    return {degree: coefficient for degree, coefficient in out.items() if coefficient != 0}


def multiply(left: dict[int, sp.Expr], right: dict[int, sp.Expr]) -> dict[int, sp.Expr]:
    out: dict[int, sp.Expr] = {}
    for i, left_coefficient in left.items():
        for j, right_coefficient in right.items():
            out[i+j] = out.get(i+j, 0) + left_coefficient*right_coefficient
    return out


def power(poly: dict[int, sp.Expr], exponent: int) -> dict[int, sp.Expr]:
    out = {0: sp.Integer(1)}
    for _ in range(exponent):
        out = multiply(out, poly)
    return out


d2_expr = a4*y**4 + a3*y**3 + a2*y**2 + a1*y + a0
sigma_expr = s3*y**3 + s2*y**2 + s1*y + s0
d0_expr = (d2_expr**2 + sigma_expr)/4
delta = -gamma**4/(4096*c)
d1_expr = delta*(y**2 + sp.Rational(25, 4)*y + v0)
e_expr = gamma*(y+1)**9
phi_expr = c*(y+1)**30*base.q

substitution_polynomials = {
    base.d0: from_expr(d0_expr),
    base.d1: from_expr(d1_expr),
    base.d2: from_expr(d2_expr),
    base.e: from_expr(e_expr),
}


def evaluate_h(expr: sp.Expr) -> dict[int, sp.Expr]:
    source = sp.Poly(expr, base.d0, base.d1, base.d2, base.e)
    out: dict[int, sp.Expr] = {}
    for monomial, coefficient in source.terms():
        term = {0: coefficient}
        for symbol, exponent in zip((base.d0, base.d1, base.d2, base.e), monomial):
            term = multiply(term, power(substitution_polynomials[symbol], exponent))
        out = add(out, term)
    return out


phi = from_expr(phi_expr)
e_poly = from_expr(e_expr)
term_cache: dict[int, tuple[dict[int, sp.Expr], ...]] = {}


def term_coefficient(f: int, target: int) -> sp.Expr:
    if f not in term_cache:
        term_cache[f] = (
            evaluate_h(H[f]),
            power(phi, f),
            power(e_poly, 21-3*f),
        )
    h_poly, phi_power, e_power = term_cache[f]
    total = 0
    for i, phi_coefficient in phi_power.items():
        for j, e_coefficient in e_power.items():
            h_degree = target-i-j
            if h_degree in h_poly:
                total += phi_coefficient*e_coefficient*h_poly[h_degree]
    return total


def master_coefficient(target: int) -> sp.Expr:
    return sum(term_coefficient(f, target) for f in range(8))


def check(target: int, forced: dict[sp.Symbol, sp.Expr], expected: sp.Expr) -> None:
    actual = master_coefficient(target).subs(forced)
    _require(sp.factor(actual-expected) == 0, target)
    print(f"  degree {target}: checked")


def main() -> None:
    forced: dict[sp.Symbol, sp.Expr] = {}
    check(238, forced, 72057594037927936*c**5*gamma**8*(32*v0-525)**2)
    forced[v0] = sp.Rational(525, 32)

    check(237, forced, -1125899906842624*201326592*c**6*gamma**3*s3**2)
    forced[s3] = 0
    check(236, forced, 4398046511104*c**3*gamma**8*(a4*gamma**3+95200*c)**2)
    forced[a4] = -sp.Integer(95200)*c/gamma**3

    check(235, forced, -1099511627776*206158430208*c**6*gamma**3*s2**2)
    forced[s2] = 0
    check(234, forced, 137438953472*32*c**3*gamma**8*(a3*gamma**3-255850*c)**2)
    forced[a3] = sp.Integer(255850)*c/gamma**3

    check(233, forced, -226673591177742970257408*c**6*gamma**3*s1**2)
    forced[s1] = 0
    check(232, forced, 4398046511104*c**3*gamma**8*(a2*gamma**3+513451*c)**2)
    forced[a2] = -sp.Integer(513451)*c/gamma**3

    check(231, forced, -1099511627776*206158430208*c**6*gamma**3*s0**2)
    forced[s0] = 0
    check(230, forced, 68719476736*c**3*gamma**8*(8*a1*gamma**3+10656467*c)**2)
    forced[a1] = -sp.Rational(10656467, 8)*c/gamma**3

    check(229, forced, sp.Integer(0))
    check(228, forced, 68719476736*c**3*gamma**8*(8*a0*gamma**3-132899897*c)**2)
    forced[a0] = sp.Rational(132899897, 8)*c/gamma**3
    check(227, forced, sp.Integer(0))
    final = 29570349989420274657771126784*c**5*gamma**8
    check(226, forced, final)
    _require(final != 0, "final != 0")

    print("a_t=9 geometric q-coprime T1, constant E: INFEASIBLE")
    print("  exact coefficients 238 down to 226 checked from f31_graded.txt")
    print("  final degree-226 coefficient is nonzero")


if __name__ == "__main__":
    main()
