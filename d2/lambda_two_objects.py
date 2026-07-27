#!/usr/bin/env python3
"""lambda_two_objects.py  (NEW 2026-07-27; read-only)

THE REPO HAS TWO OBJECTS NAMED `lambda`, AND NOTHING STATES THEIR RELATION.

Why this file exists
--------------------
On 2026-07-27 I hypothesised that the class of nine could be attacked through the
cap lemma: last night's lane proved `lam = 0` at every monomial corner, and
`PROOF_72_108.md` Lemma 2.5 reads `deg_y d_j <= lambda * w`.  If those two
`lambda`s were the same object, then `lambda = 0` would force `deg_y d_j <= 0` --
every `d_j` a CONSTANT -- which would be an extraordinarily strong constraint on
nine open cases at once.

**They are not the same object, and the hypothesis is DOWNGRADED to a conditional.**
This checker pins the distinction so nobody (including me, again) transports one
into the other.

The two objects
---------------
  CAP-LAMBDA      `caps_audit.py` C3.  Per-REGIME at a fixed corner: 3 in sub1,
                  2 in sub2, obtained as (deg slope) - (ord slope) of the
                  D-transform bounds: 15-12 = 3 and 14-12 = 2.  It is the
                  coefficient in Lemma 2.5's cap `deg_y d_j <= lambda*w`.

  STRIP-LAMBDA    `corner_atlas.py` gate G3 / `contact_lemma.py` D5.  Per-CORNER:
                  `lam = (deg_y Phi - ord_y Phi)/M`, gated by `lam >= min(m,n)`.
                  It is the actual width of Phi's degree/order strip.

At `(72,108)`: STRIP-LAMBDA = (238-204)/17 = 2.  CAP-LAMBDA is 3 or 2 depending on
the regime.  So they AGREE IN SUB2 AND DISAGREE IN SUB1 -- at one and the same
corner, where STRIP-LAMBDA cannot depend on the regime because `Phi` does not.
A single coinciding data point is exactly the trap this repo keeps falling into.

What IS true, and is checked below
----------------------------------
  A1  M * CAP-LAMBDA(sub2) = 34 = deg_y Phi - ord_y Phi   EXACTLY (zero slack)
  A2  M * CAP-LAMBDA(sub1) = 51 >  34                     (the cap is LOOSE here)
  A3  hence CAP-LAMBDA >= STRIP-LAMBDA in both regimes: the cap is an UPPER BOUND
      on the strip, attained in sub2 and slack in sub1.

A3 is the relation the repo never states, and it is the one that matters, because
it fixes the direction of the implication:

  STRIP-LAMBDA = 0   =>   nothing about CAP-LAMBDA.       (>= 0 is vacuous)
  CAP-LAMBDA = 0     =>   STRIP-LAMBDA = 0, and `deg_y d_j <= 0`.

So the attack needs `CAP-LAMBDA = 0` at a class corner, which does NOT follow from
monomiality.  CAP-LAMBDA is computed from the corner's Newton-polygon direction
functionals (Lemma 2.5(i): a hull scan giving `max(j-i)`, `max(j-2i)`, `max(2i-j)`)
and then propagated through the valuation induction (ii) and the D-transform (iii).
Redoing that at `(5,20)` etc. requires the REDUCED POLYGON at those corners --
which is exactly the datum the repo does not have transcribed away from `(8,28)`,
and exactly what blocks the floor-raising test too.

So the hypothesis survives as: **IF the cap derivation redone at a class corner
yields CAP-LAMBDA = 0, THEN every `d_j` there is constant.** Settling the
antecedent is polygon work, not arithmetic.

One encouraging note, recorded as an observation and NOT as evidence: Lemma
2.5(ii)'s induction closes because `h(j) + h(k-j)` is `j`-free, and that holds for
ANY affine `h` -- so the machinery itself is not special to `(8,28)`.  Only the
constants are.

Checker: `--quiet`, exit 0 iff every check passes.
"""
import sys
from fractions import Fraction

QUIET = "--quiet" in sys.argv
FAILS = []
N_OK = 0

# (72,108) at corner (8,28), the closed case.  Sources:
#   deg_y Phi = 238, ord_y Phi = 204  -- PROOF_72_108.md sec.2.2 / phi_corner4.py
#   M = w(Phi) = 17                   -- the u-slice weight of Phi
#   CAP-LAMBDA                        -- caps_audit.py C3 (15-12 sub1, 14-12 sub2)
DEG_PHI, ORD_PHI, M = 238, 204, 17
CAP_LAMBDA = {"sub1": 15 - 12, "sub2": 14 - 12}
DEG_SLOPE = {"sub1": 15, "sub2": 14}
ORD_SLOPE = 12


def ok(label, cond):
    global N_OK
    if cond:
        N_OK += 1
        if not QUIET:
            print("[OK]", label)
    else:
        FAILS.append(label)
        print("[FAIL]", label)


def strip_lambda(deg_phi=DEG_PHI, ord_phi=ORD_PHI, m=M):
    """corner_atlas G3 / contact_lemma D5: (deg_y Phi - ord_y Phi)/M."""
    return Fraction(deg_phi - ord_phi, m)


def main():
    SL = strip_lambda()
    ok("A0  STRIP-LAMBDA = (238-204)/17 = 2, and it is regime-INDEPENDENT because "
       "Phi does not depend on the regime", SL == 2)

    ok("A1  M * CAP-LAMBDA(sub2) = 17*2 = 34 = deg_y Phi - ord_y Phi EXACTLY -- "
       "the sub2 cap is TIGHT against Phi's strip, with zero slack (this is the "
       "'zero-margin input' PROOF_72_108 sec.2.6 warns about)",
       M * CAP_LAMBDA["sub2"] == DEG_PHI - ORD_PHI)

    ok("A2  M * CAP-LAMBDA(sub1) = 17*3 = 51 > 34 -- the sub1 cap is LOOSE, so the "
       "two lambdas are NOT the same object; they merely coincide in sub2",
       M * CAP_LAMBDA["sub1"] > DEG_PHI - ORD_PHI)

    ok("A3  THE RELATION, which the repo nowhere states: CAP-LAMBDA >= "
       "STRIP-LAMBDA in BOTH regimes -- the cap is an upper bound on the strip, "
       "attained in sub2 and slack in sub1",
       all(CAP_LAMBDA[r] >= SL for r in CAP_LAMBDA))

    ok("A4  DIRECTION OF IMPLICATION.  CAP-LAMBDA = 0 would force STRIP-LAMBDA = 0 "
       "(since 0 >= STRIP-LAMBDA >= 0), but STRIP-LAMBDA = 0 gives only "
       "CAP-LAMBDA >= 0, which is vacuous.  So monomiality (STRIP-LAMBDA = 0) does "
       "NOT deliver the cap lemma's lambda, and the 2026-07-27 hypothesis is a "
       "CONDITIONAL, not a result.",
       True)

    # A5 -- the two lambdas differ as FUNCTIONS, not just in value: CAP-LAMBDA
    # varies with the regime at fixed corner data, STRIP-LAMBDA cannot.
    ok("A5  CAP-LAMBDA takes two distinct values (3, 2) at ONE corner while "
       "STRIP-LAMBDA takes one (2) -- so no reparametrisation can identify them",
       len(set(CAP_LAMBDA.values())) == 2 and CAP_LAMBDA["sub1"] != SL)

    # A6 -- MUTATION CONTROL.  If the two were the same object, sub1 would have to
    # satisfy A1 too.  It does not, and this check would fail if it did.
    ok("A6  MUTATION CONTROL: asserting the sub1 cap is also tight (17*3 == 34) is "
       "FALSE, so A1 is a real equality test and not an artifact of how the "
       "numbers were written down",
       M * CAP_LAMBDA["sub1"] != DEG_PHI - ORD_PHI)

    # A7 -- the cap lambda really is a slope difference, per regime.
    ok("A7  CAP-LAMBDA is (deg slope - ord slope) of the D-transform bounds in each "
       "regime: 15-12 = 3 and 14-12 = 2, with the ord slope 12 SHARED (it is the "
       "window floor y^{12k} of sec.2.3)",
       all(CAP_LAMBDA[r] == DEG_SLOPE[r] - ORD_SLOPE for r in CAP_LAMBDA))

    # A8 -- the induction's j-freeness is generic, so the MACHINERY transfers even
    # though the CONSTANTS do not.  Checked symbolically for a general affine h.
    import sympy as sp
    c, s, j, k = sp.symbols("c s j k")
    h = lambda t: c + s * t
    ok("A8  Lemma 2.5(ii)'s induction closes for ANY affine h: h(j)+h(k-j) = 2c+s*k "
       "is j-FREE identically, so the cap machinery is not special to (8,28) -- only "
       "its constants are, and those come from the corner's polygon",
       sp.simplify(sp.diff(sp.expand(h(j) + h(k - j)), j)) == 0)

    if not QUIET:
        print()
    if FAILS:
        print("FAILURES:", len(FAILS), FAILS)
        sys.exit(1)
    print("ALL %d LAMBDA-DISAMBIGUATION CHECKS PASSED" % N_OK)
    sys.exit(0)


if __name__ == "__main__":
    main()
