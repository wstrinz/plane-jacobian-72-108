"""
Exact verification of the D3 dimension-3 Jacobian Conjecture counterexample.

Checks, all in exact arithmetic:
  1. The Jacobian determinant of F is the nonzero constant -2 (F is etale).
  2. The three listed source points are pairwise distinct.
  3. All three map to the common image (-1/4, 0, 0)  =>  F is not injective.
  4. The full fiber F^{-1}(-1/4, 0, 0) consists of exactly those three points.

A non-injective polynomial map with nonzero-constant Jacobian refutes JC(3).
Exits 0 on success, 1 on any failure.
"""
import sys
import sympy as sp
from counterexample import F, variables, jacobian_det, collision, evaluate

x, y, z = variables
ok = True


def check(label, condition):
    global ok
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    ok = ok and bool(condition)


print("D3 counterexample — exact verification")

# 1. constant Jacobian -2
det = jacobian_det()
check(f"Jacobian determinant == -2   (got {det})", det == -2)

# 2. distinct sources
srcs = collision['sources']
check("three source points are pairwise distinct", len(set(srcs)) == 3)

# 3. common image
img = collision['image']
for p in srcs:
    check(f"F{tuple(p)} == {tuple(img)}", evaluate(p) == tuple(img))

# 4. exact fiber — via the lex Groebner triangular basis (primary argument).
#    The fiber ideal I = <F1-img0, F2-img1, F3-img2> in Q[x,y,z].  Its reduced
#    lex Groebner basis (order z > y > x) is triangular:
#         z + 1/4 - 27/4 x^2,   y + 3/2 x,   x^3 - x.
#    We DERIVE the basis and assert it has exactly this shape.  A triangular
#    basis whose x-polynomial x^3-x = x(x-1)(x+1) has 3 simple roots, each
#    determining y and z uniquely, proves the fiber is EXACTLY three points
#    (a complete argument: no reliance on solve() enumerating every branch).
G = sp.groebner([F[0] - img[0], F[1] - img[1], F[2] - img[2]],
                z, y, x, order='lex')
gb = [sp.expand(p.as_expr()) for p in G.polys]
expected_gb = [
    sp.expand(z + sp.Rational(1, 4) - sp.Rational(27, 4) * x**2),
    sp.expand(y + sp.Rational(3, 2) * x),
    sp.expand(x**3 - x),
]
check("fiber ideal lex Groebner basis is the expected triangular form "
      "[z+1/4-27/4 x^2, y+3/2 x, x^3-x]",
      gb == expected_gb)
# x^3 - x factors into three distinct linear factors => three x-values, each
# back-substituted to a unique (y, z); enumerate them from the triangular basis.
x_roots = sp.solve(sp.Eq(x**3 - x, 0), x)
check("elimination polynomial x^3-x has 3 distinct roots", len(set(x_roots)) == 3)
gb_fiber_pts = set()
for xv in x_roots:
    yv = -sp.Rational(3, 2) * xv
    zv = -sp.Rational(1, 4) + sp.Rational(27, 4) * xv**2
    gb_fiber_pts.add((sp.nsimplify(xv), sp.nsimplify(yv), sp.nsimplify(zv)))
check(f"triangular basis yields exactly the 3 collision points "
      f"(found {len(gb_fiber_pts)})",
      gb_fiber_pts == set(srcs))

# Secondary assertion: the direct sympy.solve enumeration agrees.
eqs = [sp.Eq(F[0], img[0]), sp.Eq(F[1], img[1]), sp.Eq(F[2], img[2])]
fiber = sp.solve(eqs, [x, y, z], dict=True)
fiber_pts = {(s[x], s[y], s[z]) for s in fiber}
check(f"(secondary) sympy.solve fiber over {tuple(img)} agrees, 3 points "
      f"(found {len(fiber_pts)})",
      fiber_pts == set(srcs))

print("\nRESULT:", "all checks passed — JC(3) is refuted by this map."
      if ok else "FAILURES above.")
sys.exit(0 if ok else 1)
