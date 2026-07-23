"""
D3 — an explicit counterexample to the Jacobian Conjecture in dimension 3.

The Jacobian Conjecture JC(n) asserts: a polynomial map F: C^n -> C^n whose
Jacobian determinant is a nonzero constant is invertible (hence bijective, hence
injective). This map has constant Jacobian determinant -2 yet is NOT injective:
three distinct points share the image (-1/4, 0, 0). That refutes JC in dim 3.

Source: L. Alpoge, announced 2026-07-19. This file is our independent
reconstruction and exact re-verification of the stated map and its properties.

Run `verify.py` for the full exact check (used as a test; exits nonzero on
failure). Import `F`, `variables`, `jacobian_det`, and `collision` from here.
"""
import sympy as sp

x, y, z = sp.symbols('x y z')
variables = (x, y, z)

# F: C^3 -> C^3
F1 = (1 + x*y)**3 * z + y**2 * (1 + x*y) * (4 + 3*x*y)
F2 = y + 3*x*(1 + x*y)**2 * z + 3*x*y**2 * (4 + 3*x*y)
F3 = 2*x - 3*x**2*y - x**3*z
F = (F1, F2, F3)

# The three distinct source points that collide, and their common image.
collision = {
    'sources': [
        (sp.Integer(0),  sp.Integer(0),  sp.Rational(-1, 4)),
        (sp.Integer(1),  sp.Rational(-3, 2), sp.Rational(13, 2)),
        (sp.Integer(-1), sp.Rational(3, 2),  sp.Rational(13, 2)),
    ],
    'image': (sp.Rational(-1, 4), sp.Integer(0), sp.Integer(0)),
}


def jacobian_det():
    """Return the (expanded) Jacobian determinant of F. Should be the constant -2."""
    M = sp.Matrix(F).jacobian(variables)
    return sp.expand(M.det())


def evaluate(point):
    """Evaluate F at a point (3-tuple), returning a tuple of simplified values."""
    sub = dict(zip(variables, point))
    return tuple(sp.simplify(c.subs(sub)) for c in F)
