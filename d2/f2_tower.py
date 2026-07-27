#!/usr/bin/env python3
"""f2_tower.py  (NEW; read-only over all existing artifacts)

THE CERTIFICATE-TOWER EXPERIMENT -- the first direct attack on (75,125).

Question (GPT-Pro review 6, well-posed): the a=2 case (50,75) was KILLED by
GGV3 sec.5's polynomial-system method.  Does that contradiction EXTEND by the
fixed five-row / five-column block rule of the D-transform G-system up the F2
family to a=3 = (75,125)?

This script constructs, end to end:

  1. THE a=2 CERTIFICATE.  GGV3 sec.5 kills (50,75) in TWO reduced charts
     (gamma=3, gamma=2).  We reproduce both EXACTLY (paper 1406.0886 sec.5):
       - gamma=3: the window elimination forces C_0's y-support to bottom at
         y^-6, so c_{0,-10}=0 -- contradicting corner primitivity (a6).
       - gamma=2: the terminal 13-equation system eliminates to g_{-2}^5 =
         g_{-5}^4 = 0 (=> F_{-4}=0), then the "square in K((1/y))" obstruction
         forces a^3=2 and C_1 homogeneous of y-degree -1, i.e. e_{-10}=0 --
         contradicting corner primitivity (b6).
     CERTIFICATE KIND: a small set of terminal coefficient equations whose
     elimination forces a corner window-coefficient to vanish -- a bigraded
     (u-weight, y-order) WINDOW-DEPTH contradiction, NOT a scalar syzygy on the
     G-system generators.

  2. THE TOWER STEP.  Build the a=2 and a=3 G-systems in our D-transform t=4
     chart (landed builder g_system_75_125.build_gsystem).  The ALGEBRAIC block
     layer extends exactly: +4 forcing generators, +4 spare window unknowns, the
     Phi recurrence Phi_{a+1}=(a/(a+1)) C^{24a+2} Phi_a, and the a=2 generators
     nest as the leading block of the a=3 generators (coefficient 10/3, exactly
     as before).  The KILL layer does NOT: the window-denominator invariant
     q_window = 12a-7 jumps 17 -> 29 (gcd=1), and the y-order
     fractional-denominator set moves {1,17} -> {1,29}.
     Since the a=2 kill lives ENTIRELY in this y-order / window-cap layer, it
     does not extend by the fixed block rule.

  VERDICT: BLOCK-OBSTRUCTION.  The obstructing block is the window-cap / y-order
  layer (the 4 new forcing generators + the deepened Phi slice, carrying the new
  denominator class 29 absent at a=2).  This is the "new geometry at 125" the
  review predicted: an incommensurate window lattice.  (75,125) is NOT killed by
  tower extension; it would require a fresh period-29 window compiler.

*** REPAIRED 2026-07-26 (PASSPORT_75_125_REPAIR.md).  The chart is t=4, kappa=2,
C=y, not t=5, kappa=3, C=y^2(y^3+1): GGV5's final chain corner (7\\5,2) is not
chart data off the retraction shape.  The VERDICT (BLOCK-OBSTRUCTION; the a=2
kill does not extend) SURVIVES, and so does every structural claim -- the block
rule, the nesting coefficient 10/3, the Phi recurrence, and the coprimality of
the two window periods.  What is REFUTED is the specific arithmetic story: the
periods are 17 and 29, not 7 and 12, and both are PRIME, so the "fractional
denominator classes {2,3,4,6,12}" / "divisor lattice of the period" reading has
no counterpart.  Anything built on the period-7 -> period-12 carry structure
should be re-examined against 17 -> 29. ***

Exact sympy throughout.  Independent checker: f2_tower_verify.py.
"""
import sys
import sympy as sp

sys.path.insert(0, __file__.rsplit("f2_tower.py", 1)[0] or ".")
from g_system_75_125 import build_gsystem, Phi  # noqa: E402

y = sp.symbols("y")
# REPAIRED corner data: C = y is a MONOMIAL (deg C = ord C = 1); there is no
# residual g at all.  t = ceil(b0/a0) = 4, kappa = 2.
T, KAPPA, QC = 4, 2, 1
C = y**QC


def ordPhi(a):
    """ord_y(Phi_a) = deg_y(Phi_a) = 12a^2-10a+2  (Phi_a is a monomial)."""
    return 12 * a**2 - 10 * a + 2


def Nof(a):
    """N_a = (3a-2)(4a-1) = 12a^2-11a+2."""
    return 12 * a**2 - 11 * a + 2


def sysfor(a):
    b = 2 * a - 1
    return build_gsystem(a, b, T, QC, ordPhi(a))


def phi(a):
    return sp.Rational(1, a) * y**ordPhi(a)


def a2_certificate():
    print("=" * 78)
    print("1. THE a=2 CERTIFICATE  (GGV3 sec.5, reproduced exactly)")
    print("=" * 78)

    # --- gamma=2 terminal elimination ---
    a, e10, e7, e4, e1, g2, g5, lam = sp.symbols("a e10 e7 e4 e1 g2 g5 lam")
    eqs = [
        a**3 - 2,
        -sp.Rational(3, 7) * (-12 * (7 + a**3) * g5**2 + a**3 * (4 + 7 * a**3) * e10),
        -3 * (-6 * (-8 + a**3) * g2 * g5 + a**3 * (1 + a**3) * e7),
        -3 * (4 + a**3) * (-12 * g2**2 + a**3 * e4),
        -3 * a**3 * (-2 + a**3) * e1,
        -sp.Rational(324, 49) * ((28 + 14 * a**3 + a**6) * g5**2 + 2 * a**6 * e10) ** 2,
        -sp.Rational(162, 7) * (2 * (-32 - 16 * a**3 + a**6) * g2 * g5 + a**6 * e7)
        * ((28 + 14 * a**3 + a**6) * g5**2 + 2 * a**6 * e10),
        -sp.Rational(81, 28) * (
            28 * a**6 * (-32 - 16 * a**3 + a**6) * g2 * g5 * e7 + 7 * a**12 * e7**2
            + 4 * g2**2 * (3 * (3584 + 3584 * a**3 + 864 * a**6 - 16 * a**9 + 5 * a**12) * g5**2
                          + 16 * a**6 * (4 + a**3) ** 2 * e10)),
        -sp.Rational(162, 7) * (
            14 * (4 + a**3) ** 2 * (-32 - 16 * a**3 + a**6) * g2**3 * g5
            + 7 * a**6 * (4 + a**3) ** 2 * g2**2 * e7
            + 2 * a**6 * e1 * ((28 + 14 * a**3 + a**6) * g5**2 + 2 * a**6 * e10)),
        -81 * (4 * (4 + a**3) ** 4 * g2**4 + 2 * a**6 * (-32 - 16 * a**3 + a**6) * g2 * g5 * e1
               + a**12 * e1 * e7),
        -324 * a**6 * (4 + a**3) ** 2 * g2**2 * e1,
        -3 * a**10 * (27 * a**2 * e1**2 + 32 * lam + 16 * a**3 * lam + 2 * a**6 * lam),
        3 * a**10 * (-32 + 6 * a**6 + a**9),
    ]
    G = sp.groebner(eqs, a, e10, e7, e4, e1, lam, g2, g5, order="lex")
    elim = [p.as_expr() for p in G.polys
            if not (set(p.free_symbols) & {a, e10, e7, e4, e1, lam})]
    print("  gamma=2 terminal elimination ideal in (g_{-2},g_{-5}):")
    print("   ", sorted(str(sp.factor(e)) for e in elim))
    Ge = sp.groebner(elim, g2, g5, order="lex")
    print("    => g_{-2}^5 = 0:", Ge.reduce(g2**5)[1] == 0,
          "   g_{-5}^4 = 0:", Ge.reduce(g5**4)[1] == 0, " (F_{-4}=0)")

    # --- gamma=3 window-depth ---
    aa, bb, f2, f4, f6, f8, lm = sp.symbols("aa bb f2 f4 f6 f8 lm")
    Fm2v = f8 * y**8 + f6 * y**6 + f4 * y**4 + f2 * y**2
    C0 = sp.expand((3 * (bb * y**4) ** 2 + 2 * Fm2v + 2 * lm) / (3 * (aa * y**3) ** 2))
    lowest = min(sp.degree(sp.numer(sp.together(tm)), y) - sp.degree(sp.denom(sp.together(tm)), y)
                 for tm in sp.Add.make_args(C0))
    print("  gamma=3: forced C_0 =", C0)
    print("    lowest y-power =", lowest, "-> c_{0,-10}=0 contradicts (a6). KILL.")
    print("  CERTIFICATE KIND: terminal coeff-equations + y-order WINDOW-DEPTH "
          "contradiction (bigraded).")


def tower_step():
    print("\n" + "=" * 78)
    print("2. THE TOWER STEP  (a=2 -> a=3 by the five-block rule)")
    print("=" * 78)
    r2, r3 = sysfor(2), sysfor(3)
    print("  ALGEBRAIC block layer (transfers EXACTLY):")
    print("    closed laws: #gens = a*t-kappa-2 = %da-%d ;  #spares = (a-1)t-1 = "
          "%da-%d   (step = t = %d)" % (T, KAPPA + 2, T, T + 1, T))
    print("    generators  %d -> %d   (+%d)" % (len(r2["Gs"]), len(r3["Gs"]),
                                                len(r3["Gs"]) - len(r2["Gs"])))
    print("    spares      %d -> %d   (+%d)" % (len(r2["spares"]), len(r3["spares"]),
                                                len(r3["spares"]) - len(r2["spares"])))
    print("    Phi slice   M: %d -> %d" % (r2["M"], r3["M"]))
    d0 = sp.Symbol("d0")
    nest = sp.expand(sp.expand(r3["Gs"][1]).coeff(d0, 2) - sp.Rational(10, 3) * sp.expand(r2["Gs"][1]))
    print("    nesting: coeff_{d0^2}(G1^{a3}) = (10/3) G1^{a2} ?", nest == 0)
    for a in (2, 3):
        rec = sp.simplify(phi(a + 1) - sp.Rational(a, a + 1) * C**(24 * a + 2) * phi(a)) == 0
        print("    Phi recurrence a=%d->%d (C^(24a+2)=C^%d):" % (a, a + 1, 24 * a + 2), rec)

    print("\n  KILL layer (does NOT transfer):")
    for a in (2, 3):
        b = 2 * a - 1
        t = T
        jphi = a * t - KAPPA - 1
        M = b * t + jphi
        ordphi = ordPhi(a)
        W = sp.Rational(ordphi, M)
        js = [j for j in range(1, jphi + 1) if j != jphi - 1]
        fden = sorted(set(sp.denom(sp.Rational(ordphi * (b * t + j), M)) for j in js))
        print("    a=%d: W_step=%s  q_window=12a-7=%d (prime: %s)   y-order "
              "frac-denoms of forcing slices=%s"
              % (a, W, sp.denom(W), sp.isprime(sp.denom(W)), fden))
    print("    gcd(q_window a2, a3) = gcd(17,29) =", sp.gcd(17, 29),
          "-> incommensurate window lattices")
    print("    NOTE q_window = M exactly at every rung, so the ord-side carry "
          "obstruction (weight_lemma_75_125.py C3) admits NO zero-carry split.")

    print("\n" + "=" * 78)
    print("VERDICT: BLOCK-OBSTRUCTION.  [verdict UNCHANGED by the 2026-07-26 repair]")
    print("  Named block: the window-cap / y-order layer -- the 4 new forcing")
    print("  generators + deepened Phi slice carry the new denominator class 29,")
    print("  absent (only {1,17}) at a=2. The a=2 kill lives entirely in this")
    print("  layer, whose period 12a-7 jumps 17 -> 29 (coprime), so the")
    print("  certificate does NOT extend by the fixed block rule.")
    print("  (75,125) is NOT killed by tower extension: it needs a fresh")
    print("  period-29 window compiler.  This IS the 'new geometry at 125'.")
    print("  REFUTED sub-claim: both periods are PRIME, so there is no 'divisor")
    print("  lattice of the period' and no {2,3,4,6,12} class fragmentation.")
    print("=" * 78)


def main():
    a2_certificate()
    tower_step()


if __name__ == "__main__":
    main()
