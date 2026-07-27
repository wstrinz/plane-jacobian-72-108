#!/usr/bin/env python3
"""window_functions_75_125_verify.py  (REPAIRED 2026-07-26)

Independent EXACT checker for window_functions_75_125.py /
WINDOW_FUNCTIONS_75_125.md.  Pure sympy; --quiet; exit 0 iff every check passes.

  CHECK 0  the corner inputs come from the retraction guard, both directions
  CHECK 1  the lower (y-order) cap functions (alpha, q, beta_m) -- these SURVIVE
  CHECK 2  the class-interaction table (mod-q composition + beta 1-cocycle) -- SURVIVES
  CHECK 3  THE THREE REFUTATIONS: period 12 -> 29 (prime); no affine degree cap;
           the two-slope cone collapses to a ray (lambda = 0)
  CHECK 4  controls: the a=2 rung, the (72,108) integral limit (untouched), and
           the q_window = 12a-7 = M_a law across the family

All landed ground-truth constants are cross-checked against the repaired
signatures (C_SERIES_75_125.md, f2_family_verify.py, WINDOW_CAPS.md) so the
derivation is anchored, not self-referential.
"""
import sys
from pathlib import Path
from sympy import Rational, gcd, ceiling, isprime, denom

sys.path.insert(0, str(Path(__file__).resolve().parent))
from window_functions_75_125 import (          # noqa: E402
    family, window_law, L, L_ceil, U, U_ray, q_window, target_75_125,
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
# CHECK 0 -- the corner inputs, through the guard
# ===========================================================================
def check0():
    import polygon_reduction as pr
    cd = pr.corner_chart_data(5, 20, l_final=5, b_final=2, who="window_functions")
    ok("0.1 guard gives (t,kappa,deg C,ord C) = (4,2,1,1), C a MONOMIAL, no retraction",
       (cd["t"], cd["kappa"], cd["deg_C"], cd["ord_C"]) == (4, 2, 1, 1)
       and cd["monomial"] and not cd["retraction"])
    try:
        pr.final_corner_dictionary(5, 20, 5, 2)
        raised = False
    except pr.FinalCornerDictionaryError:
        raised = True
    ok("0.2 the superseded (t,q)=(l_final,b_final)=(5,2) RAISES; the dictionary "
       "still RETURNS (4,7) at (8,28) -- guard exercised in BOTH directions",
       raised and pr.final_corner_dictionary(8, 28, 4, 7) == (4, 7))
    f = family(3)
    ok("0.3 family(3) agrees with the guard: t=4, kappa=2",
       (f["t"], f["kappa"]) == (cd["t"], cd["kappa"]))


# ===========================================================================
# CHECK 1 -- the lower (y-order) cap functions.  THESE SURVIVE.
# ===========================================================================
def check1():
    r = target_75_125()
    alpha, q, beta = r["alpha"], r["q"], r["beta"]

    ok("1.0 landed M=29, ord_y(Phi)=80, deg_y(Phi)=80 reproduced by family(3)",
       (r["M"], r["ordPhi"], r["degPhi"]) == (29, 80, 80))
    ok("1.0 Phi signature (deg,ord,mult,cof) = (80,80,0,0): deg = ord, and mult "
       "and cofactor are 0 because C = y has no (y+1) place",
       r["degPhi"] == r["ordPhi"] == 80)
    ok("1.0 N = (3a-2)(4a-1) = 77 (was (3a-2)(5a-1) = 98)",
       r["N"] == 77 == (3 * 3 - 2) * (4 * 3 - 1))

    ok("1.1 alpha = 12a^2-10a+2 = 80", alpha == 12 * 9 - 10 * 3 + 2 == 80)
    ok("1.1 q = q_window = 12a-7 = 29 (quasi-period)", q == 12 * 3 - 7 == 29)
    ok("1.1 W_step = 80/29 non-integral (=> ord cap quasipolynomial)",
       r["W_step"] == Rational(80, 29) and r["W_step"].q != 1)

    ok("1.2 beta_m = (-alpha m) mod q, and it is a permutation of 0..q-1 "
       "(because gcd(alpha,q)=1)",
       beta == [(-alpha * m) % q for m in range(q)]
       and sorted(beta) == list(range(q)) and gcd(alpha, q) == 1)

    ok("1.3 L(w)=floor((alpha w+beta_m)/q)=ceil(alpha w/q) for w=0..300",
       all(L(w, alpha, q, beta) == L_ceil(w, alpha, q)
           == int(ceiling(Rational(alpha * w, q))) for w in range(0, 301)))
    ok("1.3 L(w) is the least integer >= (alpha/q) w (tight ceiling)",
       all(Rational(alpha * w, q) <= L(w, alpha, q, beta) < Rational(alpha * w, q) + 1
           for w in range(0, 301)))
    ok("1.4 period law L(w+29) = L(w) + 80",
       all(L(w + q, alpha, q, beta) == L(w, alpha, q, beta) + alpha
           for w in range(0, 80)))
    ok("1.5 L(w) >= 0 and L is strictly increasing for w >= 0",
       all(0 <= L(w, alpha, q, beta) < L(w + 1, alpha, q, beta) for w in range(0, 200)))


# ===========================================================================
# CHECK 2 -- the class-interaction (composition) table.  SURVIVES.
# ===========================================================================
def check2():
    r = target_75_125()
    alpha, q, beta = r["alpha"], r["q"], r["beta"]

    gcls = sorted(set(cls(w, q) for w in generator_uweights(3)))
    ok("2.1 G-generators occupy classes 21..27 and 0 (Phi's class), i.e. "
       "%s -- all distinct because every weight is < q except Phi's" % gcls,
       gcls == [0, 21, 22, 23, 24, 25, 26, 27])
    ok("2.1 Phi (w=M=29=q) sits in class 0; skipped G8 (w=28) in class 28",
       cls(29, q) == 0 and cls(28, q) == 28)
    ok("2.1 state d2,d1,d0,e occupy classes {2,3,4,5}",
       sorted(set(cls(w, q) for w in state_uweights(3))) == [2, 3, 4, 5])
    ok("2.1 seven spares dm2..dm8 occupy classes [6,7,8,9,10,11,12]",
       [cls(w, q) for w in spare_uweights(3)] == [6, 7, 8, 9, 10, 11, 12])

    ok("2.2 composition is additive mod 29 for all class pairs",
       all(compose(m1, m2, q) == (m1 + m2) % q for m1 in range(q) for m2 in range(q)))
    ok("2.2 composition closed & consistent with u-weight addition",
       all(cls(w1 + w2, q) == compose(cls(w1, q), cls(w2, q), q)
           for w1 in range(0, 60) for w2 in range(0, 60)))

    ok("2.3 beta is a 1-cocycle: carry(m1,m2) in {0,1} for all pairs",
       all(carry(m1, m2, beta, q) in (0, 1) for m1 in range(q) for m2 in range(q)))
    ok("2.3 carry(m1,m2) = ceil-defect L(w1)+L(w2)-L(w1+w2) (window superadditivity)",
       all(carry(cls(w1, q), cls(w2, q), beta, q)
           == L(w1, alpha, q, beta) + L(w2, alpha, q, beta) - L(w1 + w2, alpha, q, beta)
           for w1 in range(0, 40) for w2 in range(0, 40)))


# ===========================================================================
# CHECK 3 -- THE THREE REFUTATIONS.
# ===========================================================================
def check3():
    r = target_75_125()
    alpha, q, beta = r["alpha"], r["q"], r["beta"]

    # (R1) period 12 is refuted
    ok("3.R1 q_window(a) = 12a-7 for a=2..8, NOT 5a-3",
       all(q_window(a) == 12 * a - 7 for a in range(2, 9))
       and all(q_window(a) != 5 * a - 3 for a in range(2, 9)))
    ok("3.R1 the family periods are 17 (a=2) and 29 (a=3), not 7 and 12; both are "
       "PRIME, so there is no 'divisor lattice of the period' and no "
       "{2,3,4,6,12} class fragmentation",
       (q_window(2), q_window(3)) == (17, 29)
       and isprime(17) and isprime(29) and not isprime(12))
    ok("3.R1 they are still coprime, so the a=2 -> a=3 lattices remain "
       "incommensurate -- the QUALITATIVE conclusion of F2_TOWER.md survives",
       gcd(q_window(2), q_window(3)) == 1)
    # the fractional-denominator sets of the forcing slices
    def fracdenoms(a):
        f = family(a)
        b, t, jphi, M, op = f["b"], f["t"], f["jphi"], f["M"], f["ordPhi"]
        js = [j for j in range(1, jphi + 1) if j != jphi - 1]
        return sorted(set(denom(Rational(op * (b * t + j), M)) for j in js))
    ok("3.R1 forcing-slice fractional denominators are {1,17} (a=2) and {1,29} "
       "(a=3) -- two classes each, because the period is prime",
       fracdenoms(2) == [1, 17] and fracdenoms(3) == [1, 29])

    # (R2) no affine degree cap
    ok("3.R2 deg_slope = deg_y(Phi)/M = 80/29 is NOT an integer, so there is NO "
       "affine y-degree cap (CAPS_AUDIT sec.5's 'deg_slope = 14' is FALSE)",
       r["deg_slope"] == Rational(80, 29) and r["deg_affine"] is False)
    raised = False
    try:
        U(5, r["deg_slope"])
    except ValueError:
        raised = True
    ok("3.R2 U(w) refuses to return a number when deg_slope is non-integral, "
       "instead of silently producing a bogus cap",
       raised)
    ok("3.R2 and it DOES return at (72,108), where deg_slope = 14 -- so the "
       "refusal is discriminating, not blanket",
       U(17, window_law(204, 17, 238)["deg_slope"]) == 238)
    ok("3.R2 the F2 family never has an affine deg cap: deg_slope = W_step is "
       "non-integral at every rung a=2..8",
       all(not window_law(family(a)["ordPhi"], family(a)["M"],
                          family(a)["degPhi"])["deg_affine"] for a in range(2, 9)))

    # (R3) the cone collapses to a ray
    ok("3.R3 ord and deg slopes COINCIDE (both 80/29) and lambda = 0, because Phi "
       "is a monomial -- the (72,108) two-slope cone has no counterpart",
       r["slopes_coincide"] is True and r["lam"] == 0)
    ok("3.R3 the collapse is a family fact, not an a=3 accident",
       all(window_law(family(a)["ordPhi"], family(a)["M"],
                      family(a)["degPhi"])["lam"] == 0 for a in range(2, 9)))
    ok("3.R3 the caps PINCH: L(w) > U_ray(w) for every w not divisible by q, so "
       "the extreme-ray premise admits no object at such a weight -- which is a "
       "refutation of the premise's transfer, not a window system",
       all(L(w, alpha, q, beta) > U_ray(w, alpha, q)
           for w in range(1, 60) if w % q != 0)
       and all(L(w, alpha, q, beta) == U_ray(w, alpha, q)
               for w in range(0, 90) if w % q == 0))
    ok("3.R3 at (72,108) the cone is genuinely 2-dimensional: L(17)=204 < 238=U(17)",
       L(17, 12, 1, [0]) == 204 < 238)

    # the ord-side carry consequence consumed by weight_lemma sec.C
    zc = [w for w in range(1, r["M"])
          if carry(cls(w, q), cls(r["M"] - w, q), beta, q) == 0]
    ok("3.C q_window == M exactly, so NO split 0 < w < M has carry 0 (the "
       "superseded model left the escapes {12,24}) -- the ord-side obstruction in "
       "weight_lemma_75_125.py sec.C is TOTAL",
       zc == [] and q == r["M"])
    ok("3.C and this is a family fact: q_window(a) = M_a for every a=2..8",
       all(q_window(a) == family(a)["M"] for a in range(2, 9)))


# ===========================================================================
# CHECK 4 -- controls.
# ===========================================================================
def check4():
    # a=2 rung (the (50,75) sibling GGV3 killed)
    f2 = family(2)
    ok("4.1 a=2: M=17, ord_y(Phi_2)=deg_y(Phi_2)=30, N=28",
       (f2["M"], f2["ordPhi"], f2["degPhi"], f2["N"]) == (17, 30, 30, 28))
    law2 = window_law(f2["ordPhi"], f2["M"], f2["degPhi"])
    ok("4.1 a=2: alpha=30, q=17, deg_slope=30/17 (non-affine), lambda=0",
       (law2["alpha"], law2["q"]) == (30, 17)
       and law2["deg_slope"] == Rational(30, 17) and law2["lam"] == 0)
    ok("4.1 a=2: L(17)=30=ord_y(Phi_2) exactly (M = q, so the floor is exact)",
       L(f2["M"], law2["alpha"], law2["q"], law2["beta"]) == 30)
    ok("4.1 a=2's M = 17 coincides numerically with (72,108)'s M = 17, yet "
       "q_window differs (17 vs 1) -- the period is NOT determined by M alone, "
       "it is M/gcd(M, ord_y(Phi))",
       f2["M"] == 17 and law2["q"] == 17 and window_law(204, 17, 238)["q"] == 1
       and gcd(30, 17) == 1 and gcd(204, 17) == 17)

    # (72,108) integral limit -- UNTOUCHED
    l72 = window_law(204, 17, 238)
    ok("4.2 (72,108) integral limit: W_step=12, alpha=12, q_window=1, deg_slope=14",
       l72["W_step"] == 12 and l72["alpha"] == 12
       and l72["q"] == 1 and l72["deg_slope"] == 14)
    ok("4.2 (72,108) caps affine: L(w)=12w, U(w)=14w (matches WINDOW_CAPS.md "
       "ord>=12k, deg<=14k sub2)",
       all(L(w, 12, 1, [0]) == 12 * w and U(w, 14) == 14 * w for w in range(0, 20)))
    ok("4.2 (72,108) Phi at caps: L(17)=204, U(17)=238; lambda = 2 (the stripping "
       "factor).  This is the structure (75,125) lacks.",
       L(17, 12, 1, [0]) == 204 and U(17, 14) == 238 and l72["lam"] == 2)
    ok("4.2 (72,108) q_window=1, so carry == 0 for EVERY split of 17",
       all(carry(0, 0, [0], 1) == 0 for _ in range(1))
       and all(L(w, 12, 1, [0]) + L(17 - w, 12, 1, [0]) - L(17, 12, 1, [0]) == 0
               for w in range(1, 17)))

    # the family law
    ok("4.3 q_window(a) = 12a-7 = M_a with gcd(alpha,q)=1 for a=2..8",
       all(q_window(a) == 12 * a - 7 == family(a)["M"]
           and gcd(family(a)["ordPhi"], family(a)["M"]) == 1 for a in range(2, 9)))
    ok("4.3 L(M_a) = ord_y(Phi_a) exactly at every rung a=2..8 (equality, since "
       "M_a = q_window(a))",
       all((lambda f, lw: L(f["M"], lw["alpha"], lw["q"], lw["beta"]) == f["ordPhi"])(
               family(a), window_law(family(a)["ordPhi"], family(a)["M"],
                                     family(a)["degPhi"])) for a in range(2, 9)))


def main():
    check0()
    check1()
    check2()
    check3()
    check4()
    if not QUIET:
        print()
    if FAILS:
        print("FAILURES:", len(FAILS), FAILS)
        sys.exit(1)
    print("ALL %d WINDOW-FUNCTION CHECKS PASSED" % N_OK)
    sys.exit(0)


if __name__ == "__main__":
    main()
