#!/usr/bin/env python3
"""marked_polynomial.py  --  GGHV's marked-polynomial generator + the corner law

GGHV22 prints exactly ONE explicit marked polynomial (the (66,99) / A_0=(9,24)
case, tex line 1594).  Our own home case (72,108) / A_0=(8,28) has another,
recorded in AUDIT.md item 4.  Both come from the SAME forcing ODE, which this
module reconstructs and solves for arbitrary parameters -- turning "read the
paper for a polynomial" into "generate the polynomial for any case".

THE FORCING ODE.  With C(y) = y^(a-1)(y+1) the shape polynomial and
f1 := C^3 * F_{-v}:

    2c * C * f1'  -  2s * C' * f1  =  C^2 ,        s = c + e + 1

Instantiations (both reproduced byte-exact by this module):
    (a,c,e) = (9,3,1) ->  6*C*f' - 10*C'*f = C^2   (GGHV tex line 1588)
    (a,c,e) = (8,4,2) ->  8*C*f' - 14*C'*f = C^2   (AUDIT.md item 4)

Structural facts recovered from it:
  * a polynomial solution exists iff c | a, and it is then unique;
  * the solution always has the form f1 = y^a (y+1)^2 g, so f = f1/C = y(y+1) g;
  * hence GGHV's "deg f = 6, separable, y(y+1) | f" forces deg g = 4;
  * with n := deg g = (c-1)(a/c - 1) + (c-3), the equation n = 4 has EXACTLY
    four integer solutions (a,c) = (12,2), (9,3), (8,4), (7,7).
    Two are GGHV's published cases; two are new.

THE CORNER LAW (see CORNER_RESOLVENT.md).  With C_int := (a-1)(c-1) - c and
n = deg g EVEN:

    disc(g)  =  (-1)^(n/2) * sqfree(C_int) * (perfect square)

i.e. the quadratic resolvent of the marked polynomial is
Q(sqrt((-1)^(n/2) * C_int)).  FALSE for odd n (9/9 counterexamples).

Read-only.  Usage:
    python marked_polynomial.py             # anchors + the four deg-4 cases
    python marked_polynomial.py --sweep     # law check over a range of (a,c)
    python marked_polynomial.py --quiet     # self-check, exit 0 iff all pass
"""
from __future__ import annotations

import argparse

import sympy as sp

y = sp.Symbol("y")

# The two published anchors, verbatim (primitive, positive leading coefficient).
ANCHORS = {
    (9, 3): "243*y**4 - 81*y**3 + 54*y**2 - 42*y + 35",       # GGHV tex 1594, (66,99)
    (8, 4): "2048*y**4 - 512*y**3 + 320*y**2 - 240*y + 195",  # AUDIT.md item 4, (72,108)
}


def corner_integer(a, c):
    """C = q(kappa+1) - t under the dictionary (t,kappa,q) = (c, c-2, a-1).

    NOTE: this dictionary is FITTED from the two anchors, not read off GGHV's
    Newton-polygon data.  The arithmetic law below is solid; the identification
    of this quantity with the census corner integer rests on the dictionary.
    See CORNER_RESOLVENT.md sec.5 -- this is the weakest link.
    """
    return (a - 1) * (c - 1) - c


def marked(a, c, e=None):
    """Primitive marked polynomial g for parameters (a, c); None if none exists."""
    if e is None:
        e = c - 2                      # kappa = t - 2, the standard class
    s = c + e + 1
    C = y**(a - 1) * (y + 1)
    ansatz_deg = a + 8
    co = sp.symbols("u0:%d" % (ansatz_deg + 1))
    f1 = sum(co[i] * y**i for i in range(ansatz_deg + 1))
    residual = sp.expand(2 * c * C * sp.diff(f1, y) - 2 * s * sp.diff(C, y) * f1 - C**2)
    sol = sp.solve([sp.Eq(k, 0) for k in sp.Poly(residual, y).all_coeffs()], co, dict=True)
    if not sol:
        return None
    f1v = sp.expand(f1.subs(sol[0]))
    if f1v == 0 or f1v.free_symbols - {y}:
        return None
    quotient = sp.cancel(sp.together(f1v / (y**a * (y + 1)**2)))
    if quotient.has(y**-1) or sp.denom(quotient).has(y):
        return None
    g = sp.Poly(sp.expand(quotient), y)
    g = sp.Poly(g.as_expr() * sp.lcm([t.q for t in g.all_coeffs()]), y)
    if g.LC() < 0:
        g = sp.Poly(-g.as_expr(), y)
    return sp.Poly(g.as_expr() / sp.gcd([abs(t) for t in g.all_coeffs()]), y)


def squarefree_part(value):
    value = sp.Integer(value)
    if value == 0:
        return sp.Integer(0)
    out = sp.Integer(-1) if value < 0 else sp.Integer(1)
    for p, m in sp.factorint(abs(value)).items():
        if m % 2:
            out *= p
    return out


def law_check(a, c):
    """Return a row testing the corner law at (a, c)."""
    g = marked(a, c)
    if g is None:
        return None
    n = g.degree()
    disc = sp.discriminant(g.as_expr(), y)
    lhs = squarefree_part(disc)
    Ci = corner_integer(a, c)
    rhs = squarefree_part((-1)**(n // 2) * Ci) if n % 2 == 0 else None
    return {"a": a, "c": c, "n": n, "C": Ci, "disc_sqfree": lhs,
            "predicted": rhs, "even": n % 2 == 0,
            "holds": (lhs == rhs) if rhs is not None else None, "g": g.as_expr()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sweep", action="store_true")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()
    failures = []
    checks = 0

    # --- 1. anchors reproduce byte-exact -------------------------------------
    for (a, c), text in ANCHORS.items():
        g = marked(a, c)
        checks += 1
        ok = g is not None and sp.expand(g.as_expr() - sp.sympify(text)) == 0
        if not ok:
            failures.append("anchor (a=%d,c=%d) NOT reproduced: got %s" % (a, c, g))
        if not args.quiet:
            print("anchor (a=%2d,c=%d): %s  %s"
                  % (a, c, g.as_expr() if g else None, "OK" if ok else "FAIL"))

    # --- 2. the four deg-4 (GGHV-shaped) cases -------------------------------
    if not args.quiet:
        print("\nthe four deg-g=4 solutions (GGHV shape: deg f = 6):")
    for (a, c) in ((12, 2), (9, 3), (8, 4), (7, 7)):
        row = law_check(a, c)
        checks += 1
        if row is None or row["n"] != 4:
            failures.append("(a=%d,c=%d) expected deg g = 4" % (a, c))
            continue
        if not row["holds"]:
            failures.append("(a=%d,c=%d) corner law FAILS: %s vs %s"
                            % (a, c, row["disc_sqfree"], row["predicted"]))
        if not args.quiet:
            print("  (a=%2d,c=%d)  C=%-3d  disc_sqfree=%-4s  predicted=%-4s  %s"
                  % (a, c, row["C"], row["disc_sqfree"], row["predicted"],
                     "OK" if row["holds"] else "FAIL"))

    # --- 3. optional sweep ---------------------------------------------------
    if args.sweep and not args.quiet:
        print("\nsweep (even n = law applies; odd n = law is FALSE by design):")
        for c in range(2, 8):
            for a in range(c, 25):
                if a % c:
                    continue
                row = law_check(a, c)
                if row is None:
                    continue
                tag = ("holds" if row["holds"] else "FAIL") if row["even"] else "odd-n (n/a)"
                print("   (a=%2d,c=%d) n=%-3d C=%-4d disc_sqfree=%-10s %s"
                      % (a, c, row["n"], row["C"], row["disc_sqfree"], tag))

    if failures:
        print("\nFAILURES (%d of %d):" % (len(failures), checks))
        for f in failures:
            print("  -", f)
        return 1
    print("\nALL %d MARKED-POLYNOMIAL CHECKS PASSED" % checks)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
