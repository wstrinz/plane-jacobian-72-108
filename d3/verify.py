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

# 4. exact fiber
eqs = [sp.Eq(F[0], img[0]), sp.Eq(F[1], img[1]), sp.Eq(F[2], img[2])]
fiber = sp.solve(eqs, [x, y, z], dict=True)
fiber_pts = {(s[x], s[y], s[z]) for s in fiber}
check(f"fiber over {tuple(img)} has exactly the 3 collision points "
      f"(found {len(fiber_pts)})",
      fiber_pts == set(srcs))

print("\nRESULT:", "all checks passed — JC(3) is refuted by this map."
      if ok else "FAILURES above.")
sys.exit(0 if ok else 1)
