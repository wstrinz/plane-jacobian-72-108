#!/usr/bin/env python3
"""window_functions_75_125_verify.py  (NEW; read-only over all existing artifacts)

Independent EXACT checker for window_functions_75_125.py / WINDOW_FUNCTIONS_75_125.md
(the period-12 window functions for (75,125)).  Pure sympy; --quiet; exit 0 iff
every check passes.  Structured after the three deliverables:

  CHECK 1  the floor/ceiling window-cap functions (alpha, q, beta_m, deg_slope)
  CHECK 2  the class-interaction table (mod-q composition + beta 1-cocycle)
  CHECK 3  consistency:  (a) Phi at its caps exactly;
                         (b) a=2 period-7 control against f2_tower.py's window
                             table + the (72,108) affine (integral) limit;
                         (c) the q_window = 5a-3 law.

All landed ground-truth constants are cross-checked against the published
signatures (C_SERIES_75_125.md, f2_family_verify.py, WINDOW_CAPS.md) so the
derivation is anchored, not self-referential.
"""
import sys
from sympy import Rational, gcd, ceiling, floor, denom

sys.path.insert(0, __file__.rsplit("window_functions_75_125_verify.py", 1)[0] or ".")
from window_functions_75_125 import (          # noqa: E402
    family, window_law, L, L_ceil, U, q_window, target_75_125,
    generator_uweights, state_uweights, spare_uweights, cls, compose, carry,
)

QUIET = "--quiet" in sys.argv
FAILS = []
N_OK = 0


def ok(label, cond):
    global N_OK
    if cond:
        N_OK += 1
        if not QUIET:
            print("[OK]", label)
    else:
        FAILS.append(label)
        print("[FAIL]", label)


# ===========================================================================
# CHECK 1 -- the window-cap functions
# ===========================================================================
def check1():
    r = target_75_125()
    alpha, q, beta, dsl = r["alpha"], r["q"], r["beta"], r["deg_slope"]

    # 1.0 landed anchors: the built (75,125) Phi slice (C_SERIES_75_125.md)
    ok("1.0 landed M=36, ord_y(Phi)=201, deg_y(Phi)=504 reproduced by family(3)",
       (r["M"], r["ordPhi"], r["degPhi"]) == (36, 201, 504))
    ok("1.0 Phi signature (deg,ord,mult,cof)=(504,201,101,202): deg=ord+mult+cof",
       r["degPhi"] == 201 + 101 + 202 and r["ordPhi"] == 201)

    # 1.1 the derived constants
    ok("1.1 alpha = 10a^2-8a+1 = 67", alpha == 10 * 9 - 8 * 3 + 1 == 67)
    ok("1.1 q = 5a-3 = 12 (quasi-period)", q == 5 * 3 - 3 == 12)
    ok("1.1 deg_slope = 5a-1 = 14", dsl == 5 * 3 - 1 == 14)
    ok("1.1 W_step = 67/12 non-integral (=> ord cap quasipolynomial)",
       r["W_step"] == Rational(67, 12) and r["W_step"].q != 1)

    # 1.2 the twelve beta_m, derived as (-alpha m) mod q
    ok("1.2 beta_m = (-alpha m) mod q = [0,5,10,3,8,1,6,11,4,9,2,7]",
       beta == [0, 5, 10, 3, 8, 1, 6, 11, 4, 9, 2, 7]
       and beta == [(-alpha * m) % q for m in range(q)])

    # 1.3 floor form == ceil form == tight integer lower bound, for all w
    ok("1.3 L(w)=floor((alpha w+beta_m)/q)=ceil(alpha w/q) for w=0..200",
       all(L(w, alpha, q, beta) == L_ceil(w, alpha, q)
           == int(ceiling(Rational(alpha * w, q))) for w in range(0, 201)))
    # 1.3 the tight-lower-bound property: L(w) >= alpha*w/q > L(w)-1
    ok("1.3 L(w) is the least integer >= (alpha/q) w (tight ceiling)",
       all(Rational(alpha * w, q) <= L(w, alpha, q, beta) < Rational(alpha * w, q) + 1
           for w in range(0, 201)))

    # 1.4 quasi-periodicity: L(w+q) = L(w) + alpha ; U affine
    ok("1.4 period law L(w+12) = L(w) + 67",
       all(L(w + q, alpha, q, beta) == L(w, alpha, q, beta) + alpha for w in range(0, 60)))
    ok("1.4 U(w) = 14 w affine, integral slope (deg cap NOT quasipolynomial)",
       all(U(w, dsl) == 14 * w for w in range(0, 60)))

    # 1.5 caps are ordered 0 <= L(w) <= U(w) (an object's ord <= its deg)
    ok("1.5 0 <= L(w) <= U(w) for all w>=0",
       all(0 <= L(w, alpha, q, beta) <= U(w, dsl) for w in range(0, 200)))


# ===========================================================================
# CHECK 2 -- the class-interaction (composition) table
# ===========================================================================
def check2():
    r = target_75_125()
    alpha, q, beta = r["alpha"], r["q"], r["beta"]

    # 2.1 generator / Phi / state / spare residue-class occupancy
    gcls = sorted(set(cls(w, q) for w in generator_uweights(3)))
    ok("2.1 G-generators occupy classes {0,2,3,4,5,6,7,8,9,10} (all but {1,11})",
       gcls == [0, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    ok("2.1 Phi (w=M=36) sits in class 0; skipped G10 (w=35) in class 11",
       cls(36, q) == 0 and cls(35, q) == 11)
    ok("2.1 state d3..d0,e occupy classes {2,3,4,5,6}",
       sorted(set(cls(w, q) for w in state_uweights(3))) == [2, 3, 4, 5, 6])
    ok("2.1 nine spares dm2..dm10 occupy classes [7,8,9,10,11,0,1,2,3]",
       [cls(w, q) for w in spare_uweights(3)] == [7, 8, 9, 10, 11, 0, 1, 2, 3])

    # 2.2 the 12x12 composition table: class(x*y) = (class x + class y) mod q,
    #     because multiplication ADDS u-weights.
    ok("2.2 composition is additive mod 12 for all class pairs",
       all(compose(m1, m2, q) == (m1 + m2) % q for m1 in range(q) for m2 in range(q)))
    ok("2.2 composition closed & consistent with u-weight addition",
       all(cls(w1 + w2, q) == compose(cls(w1, q), cls(w2, q), q)
           for w1 in range(0, 40) for w2 in range(0, 40)))

    # 2.3 the beta offsets form a 1-cocycle: beta_{m1}+beta_{m2}
    #     = beta_{(m1+m2)%q} + q*carry, carry in {0,1}; the carry is exactly
    #     the ceil-superadditivity defect  L(w1)+L(w2) - L(w1+w2) in {0,1}.
    ok("2.3 beta is a 1-cocycle: carry(m1,m2) in {0,1} for all pairs",
       all(carry(m1, m2, beta, q) in (0, 1) for m1 in range(q) for m2 in range(q)))
    ok("2.3 carry(m1,m2) = ceil-defect L(w1)+L(w2)-L(w1+w2) (window superadditivity)",
       all(carry(cls(w1, q), cls(w2, q), beta, q)
           == L(w1, alpha, q, beta) + L(w2, alpha, q, beta) - L(w1 + w2, alpha, q, beta)
           for w1 in range(0, 30) for w2 in range(0, 30)))


# ===========================================================================
# CHECK 3 -- consistency (Phi at caps; a=2 & (72,108) controls; q_window law)
# ===========================================================================
def check3():
    r = target_75_125()
    alpha, q, beta, dsl = r["alpha"], r["q"], r["beta"], r["deg_slope"]

    # 3(a) the known Phi point sits AT both caps, at EQUALITY (the way (72,108)'s
    #      tight rows do).  M = 36 = 3q is a multiple of q, so the floor is exact.
    ok("3a M = 36 is a multiple of q=12 (=> quasipolynomial floor is exact at Phi)",
       r["M"] % q == 0)
    ok("3a L(M) = ord_y(Phi) = 201 exactly (lower cap attained at equality)",
       L(r["M"], alpha, q, beta) == r["ordPhi"] == 201)
    ok("3a U(M) = deg_y(Phi) = 504 exactly (upper cap attained at equality)",
       U(r["M"], dsl) == r["degPhi"] == 504)
    ok("3a stripped Phi degree U(M)-L(M) = 504-201 = 303 = deg of (y^3+1)^101 cofactor part",
       U(r["M"], dsl) - L(r["M"], alpha, q, beta) == 303)

    # 3(b) CONTROL 1: reduce to the a=2 period-7 analogue (the (50,75) rung).
    #      Landed Phi_2 signature (deg,ord)=(189,75) at M=21 (f2_family_verify.py).
    f2 = family(2)
    ok("3b(a=2) landed M=21, ord_y(Phi_2)=75, deg_y(Phi_2)=189 reproduced",
       (f2["M"], f2["ordPhi"], f2["degPhi"]) == (21, 75, 189))
    law2 = window_law(f2["ordPhi"], f2["M"], f2["degPhi"])
    ok("3b(a=2) alpha=25, q=7, deg_slope=9 (period-7 analogue)",
       (law2["alpha"], law2["q"], law2["deg_slope"]) == (25, 7, 9))
    ok("3b(a=2) beta_m = (-25 m) mod 7 = [0,3,6,2,5,1,4]",
       law2["beta"] == [0, 3, 6, 2, 5, 1, 4] == [(-25 * m) % 7 for m in range(7)])
    ok("3b(a=2) Phi_2 at caps: L(21)=75, U(21)=189 (M=21=3q exact)",
       L(f2["M"], law2["alpha"], law2["q"], law2["beta"]) == 75
       and U(f2["M"], law2["deg_slope"]) == 189 and f2["M"] % law2["q"] == 0)
    # tie to f2_tower.py's window table: naive physical-order fractional
    # denominators of the forcing slices (F2_TOWER.md sec.2b / f2_tower_verify §G).
    def fracdenoms(a):
        f = family(a)
        b, t, jphi, M, ordphi = f["b"], f["t"], f["jphi"], f["M"], f["ordPhi"]
        js = [j for j in range(1, jphi + 1) if j != jphi - 1]
        return sorted(set(denom(Rational(ordphi * (b * t + j), M)) for j in js))
    ok("3b control vs f2_tower window table: frac-denoms {1,7} (a=2) -> "
       "{1,2,3,4,6,12} (a=3), the divisor lattice of the period",
       fracdenoms(2) == [1, 7] and fracdenoms(3) == [1, 2, 3, 4, 6, 12])
    ok("3b these frac-denoms are exactly the divisors of q that occur "
       "(period 7: {1,7}; period 12: divisors of 12 except {8,12->..})",
       all(d == 1 or (q % d == 0) for d in fracdenoms(3)))

    # 3(b) CONTROL 2: the (72,108) INTEGRAL limit (q_window=1).  Landed Phi
    #      signature (238,204) at M=17 (WINDOW_CAPS.md / STATE.md); the caps
    #      degenerate to the affine ord>=12w, deg<=14w recited there.
    law_72 = window_law(204, 17, 238)
    ok("3b(72,108) integral limit: W_step=12, alpha=12, q_window=1, deg_slope=14",
       law_72["W_step"] == 12 and law_72["alpha"] == 12
       and law_72["q"] == 1 and law_72["deg_slope"] == 14)
    ok("3b(72,108) caps affine: L(w)=12w, U(w)=14w (matches WINDOW_CAPS.md ord>=12k, "
       "deg<=14k sub2)",
       all(L(w, 12, 1, [0]) == 12 * w and U(w, 14) == 14 * w for w in range(0, 20)))
    ok("3b(72,108) Phi at caps: L(17)=204, U(17)=238 (deg 238=14*17, ord 204=12*17)",
       L(17, 12, 1, [0]) == 204 and U(17, 14) == 238)

    # 3(c) the q_window = 5a-3 law with gcd(alpha,q)=1 (F2 window-denominator law)
    ok("3c q_window(a) = denom(W_step) = 5a-3 for a=2..6",
       all(q_window(a) == 5 * a - 3 for a in range(2, 7)))
    ok("3c gcd(alpha, q_window) = 1 for a=2..6 (reduced; incommensurate lattices)",
       all(gcd(10 * a ** 2 - 8 * a + 1, 5 * a - 3) == 1 for a in range(2, 7)))
    ok("3c gcd(q_window(2), q_window(3)) = gcd(7,12) = 1 (a=2 -> a=3 incommensurate)",
       gcd(q_window(2), q_window(3)) == 1)

    # 3 general: family formulas reproduce Phi at caps at EVERY rung a=2..6
    ok("3 general: L(M)=ord_y(Phi) and U(M)=deg_y(Phi) exactly for a=2..6",
       all((lambda f, lw: L(f["M"], lw["alpha"], lw["q"], lw["beta"]) == f["ordPhi"]
            and U(f["M"], lw["deg_slope"]) == f["degPhi"])(
               family(a), window_law(family(a)["ordPhi"], family(a)["M"], family(a)["degPhi"]))
           for a in range(2, 7)))


def main():
    check1()
    check2()
    check3()
    if not QUIET:
        print()
    if FAILS:
        print("FAILURES:", len(FAILS), FAILS)
        sys.exit(1)
    print("ALL %d WINDOW-FUNCTION CHECKS PASSED" % N_OK)
    sys.exit(0)


if __name__ == "__main__":
    main()
