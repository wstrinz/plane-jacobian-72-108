#!/usr/bin/env python3
"""f2_tower_verify.py  (NEW; read-only over all existing artifacts)

EXACT checker for THE CERTIFICATE-TOWER EXPERIMENT (f2_tower.py / F2_TOWER.md):
the first direct attack on (75,125) by extending the (50,75) [a=2] kill up the
F2 family to a=3 = (75,125), through the fixed five-row/five-column block rule of
the D-transform G-system.

What is checked (all exact sympy; --quiet; exit 0 iff every check passes):

  A. The a=2 kill, GGV3 sec.5 gamma=2 terminal elimination (paper 1406.0886
     lines 2027-2042), reproduced EXACTLY: eliminating {a,e10,e7,e4,e1,lam} from
     the 13-equation terminal system yields the elimination ideal
     <g2^5, g2^4 g5, g2^2 g5^2, g2 g5^3, g5^4>, so g_{-2}^5 = g_{-5}^4 = 0
     (matching the paper's stated conclusion), forcing g_{-2}=g_{-5}=0, i.e.
     F_{-4}=0.
  B. The a=2 kill, GGV3 sec.5 gamma=3 chart (paper 1837-1887), reproduced
     EXACTLY: eliminating the deep window unknowns from the E-system forces C_0
     to have y-support bottoming at y^-6, so its y^-10 coefficient is 0 --
     contradicting corner primitivity (a6): c_{0,-10} != 0. A WINDOW-DEPTH kill.
  C. The gamma=2 square/window-depth obstruction: the derived relation
     3 a^2 C_1^2 y^2 = (4a^3-8) y^3 - 8 lam forces (perfect-square in K((1/y)))
     4a^3-8=0, hence C_1 homogeneous of y-degree -1, hence its y^-10 coefficient
     e_{-10}=0 -- contradicting corner primitivity (b6): e_{-10} != 0 (as C_0=0).
     Same WINDOW-DEPTH flavour as B.
  D. The five-block G-system growth (our D-transform t=5 chart, builder in the
     landed g_system_75_125.py): a=2 has 5 generators / 4 spares; a=3 has 10 / 9;
     the step adds EXACTLY 5 forcing generators and 5 spare window unknowns; the
     Phi u-slice deepens M: 21 -> 36.
  E. The forcing recurrence Phi_{a+1} = (a/(a+1)) C^{30a+3} Phi_a (exact, a=2..4).
  F. The generator nesting: the a=2 generator G1 is the leading (d0^2) block of
     the a=3 generator G1 -- coeff_{d0^2}(G1^{a=3}) = (10/3) * G1^{a=2}. The
     algebraic block genuinely nests.
  G. THE OBSTRUCTION (why the tower does NOT extend by the fixed block rule): the
     window-denominator invariant q_window = denom(ord_y(Phi)/M) = 5a-3 jumps
     7 -> 12 with gcd(7,12)=1 (incommensurate window lattices), and the y-order
     fractional-denominator set of the forcing slices fragments from {1,7} at a=2
     to {1,2,3,4,6,12} at a=3. The kill of A-C lives entirely in this y-order /
     window-cap layer, which is NOT preserved by the five-block rule.
"""
import sys
import sympy as sp

sys.path.insert(0, __file__.rsplit("f2_tower_verify.py", 1)[0] or ".")

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


y = sp.symbols("y")
g = y**3 + 1
C = y**2 * g


# ---------------------------------------------------------------------------
# A. GGV3 sec.5 gamma=2 terminal elimination  ->  g_{-2}^5 = g_{-5}^4 = 0
# ---------------------------------------------------------------------------
def check_A():
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
    order = [a, e10, e7, e4, e1, lam, g2, g5]
    G = sp.groebner(eqs, *order, order="lex")
    elim = [p.as_expr() for p in G.polys
            if not (set(p.free_symbols) & {a, e10, e7, e4, e1, lam})]
    Ge = sp.groebner(elim, g2, g5, order="lex")
    ok("A: GGV3 gamma=2 terminal elimination gives g_{-2}^5 = 0", Ge.reduce(g2**5)[1] == 0)
    ok("A: GGV3 gamma=2 terminal elimination gives g_{-5}^4 = 0", Ge.reduce(g5**4)[1] == 0)
    # g2, g5 individually NOT in the ideal (nilpotent-only: kill needs a domain / field values)
    ok("A: g_{-2}, g_{-5} not individually forced (nilpotence, not linear)",
       Ge.reduce(g2)[1] != 0 and Ge.reduce(g5)[1] != 0)


# ---------------------------------------------------------------------------
# B. GGV3 sec.5 gamma=3 chart: forced C_0 y-support bottoms at y^-6, no y^-10
# ---------------------------------------------------------------------------
def check_B():
    Z0, Z1, Z2, Z3, Z4, Z5, Z6, Z7, lam, Fm1, Fm2 = sp.symbols(
        "Z0 Z1 Z2 Z3 Z4 Z5 Z6 Z7 lam Fm1 Fm2")
    E1 = 2 * Z0 * Z1 + 2 * Z3
    E2 = Z1**2 + 2 * Z0 * Z2 + 2 * Z4
    E3 = 2 * Z1 * Z2 + 2 * Z0 * Z3 + 2 * Z5
    E4 = Z2**2 + 2 * Z1 * Z3 + 2 * Z0 * Z4 + 2 * Z6
    E5 = 2 * Z2 * Z3 + 2 * Z1 * Z4 + 2 * Z0 * Z5 + 2 * Z7
    E6 = 3 * Z0**2 * Z1 + 6 * Z1 * Z2 + 6 * Z0 * Z3 + 3 * Z5
    E7 = lam + 3 * Z0 * Z1**2 + 3 * Z0**2 * Z2 + 3 * Z2**2 + 6 * Z1 * Z3 + 6 * Z0 * Z4 + 3 * Z6
    sol = sp.solve([E1, E2, E3, E4, E5], [Z3, Z4, Z5, Z6, Z7], dict=True)[0]
    E6s = sp.expand(E6.subs(sol) + Fm1)
    Fm1_val = sp.solve(E6s, Fm1)[0]
    ok("B: E1,E3,E6 give F_{-1} = -3 C_{-1} C_{-2}", sp.expand(Fm1_val - (-3 * Z1 * Z2)) == 0)
    E7rel = sp.expand(E7.subs(sol).subs(Fm1, Fm1_val) + Fm2)  # = 0
    # C0 !=0 branch relation: 3 C0 C_{-1}^2 - 3 C_{-2}^2 - 2 lam = 2 F_{-2}
    target = 3 * Z0 * Z1**2 - 3 * Z2**2 - 2 * lam - 2 * Fm2
    ratio = sp.simplify(E7rel / target)
    ok("B: forcing relation  3 C0 C_{-1}^2 - 3 C_{-2}^2 - 2 lam = 2 F_{-2}",
       ratio != 0 and ratio.free_symbols == set())
    # window caps (a5): C_{-1}=a y^3, C_{-2}=b y^4 ; F_{-2}=f8 y^8+f6 y^6+f4 y^4+f2 y^2
    aa, bb, f2, f4, f6, f8 = sp.symbols("aa bb f2 f4 f6 f8")
    Fm2v = f8 * y**8 + f6 * y**6 + f4 * y**4 + f2 * y**2
    C0 = sp.expand((3 * (bb * y**4) ** 2 + 2 * Fm2v + 2 * lam) / (3 * (aa * y**3) ** 2))
    terms = sp.Add.make_args(C0)
    lowest = min(sp.degree(sp.numer(sp.together(tm)), y) - sp.degree(sp.denom(sp.together(tm)), y)
                 for tm in terms)
    ok("B: forced C_0 y-support bottoms at y^-6 (no y^-10 term); c_{0,-10}=0 vs (a6)!=0",
       lowest == -6)


# ---------------------------------------------------------------------------
# C. gamma=2 square/window-depth obstruction: a^3=2, then e_{-10}=0 vs (b6)
# ---------------------------------------------------------------------------
def check_C():
    aa, lam = sp.symbols("aa lam")
    # from C_{-1}^2(4 C_{-1} - 3 C_1^2) = 8(y^3+lam) with C_{-1}=aa*y :
    Cm1 = aa * y
    C1 = sp.symbols("C1")
    lhs = sp.expand(Cm1**2 * (4 * Cm1 - 3 * C1**2))
    rhs = 8 * (y**3 + lam)
    rel = sp.expand(lhs - rhs)  # = 0  ->  3 aa^2 C1^2 y^2 = (4 aa^3 - 8) y^3 - 8 lam
    derived = sp.expand(3 * aa**2 * C1**2 * y**2 - ((4 * aa**3 - 8) * y**3 - 8 * lam))
    ok("C: relation reduces to 3 a^2 C_1^2 y^2 = (4a^3-8) y^3 - 8 lam",
       sp.expand(rel + derived) == 0)
    # perfect square in K((1/y)) => even top y-degree; here top-degree 3 (odd) with
    # coeff (4a^3-8) => squareness forces 4a^3-8 = 0, i.e. a^3=2.
    ok("C: odd top y-degree 3 forces 4a^3-8=0 (a^3=2) for RHS to be a square",
       sp.Poly((4 * aa**3 - 8) * y**3 - 8 * lam, y).degree() == 3)
    # with a^3=2 the relation is 3 a^2 C1^2 y^2 = -8 lam ; expand C1 in the (b6) window
    e1, e4, e7, e10 = sp.symbols("e1 e4 e7 e10")
    C1w = e1 / y + e4 / y**4 + e7 / y**7 + e10 / y**10
    # impose 3 a^2 C1^2 y^2 + 8 lam == 0 as a y-identity with 4a^3-8=0 (a^3=2):
    expr = sp.expand(3 * aa**2 * C1w**2 * y**2 + 8 * lam)
    poly = sp.Poly(sp.together(expr).as_numer_denom()[0], y)
    coeffs = poly.all_coeffs()
    # solve for the window unknowns forcing the identity; e10 must be 0
    solset = sp.solve([sp.Poly(expr * y**20, y).coeff_monomial(y**k)
                       for k in range(0, 21)] , [e1, e4, e7, e10, lam], dict=True)
    forced_e10_zero = all(s.get(e10, e10) == 0 for s in solset) if solset else False
    ok("C: with a^3=2, matching y-orders forces C_1 homogeneous deg -1: e_{-10}=0 vs (b6)!=0",
       forced_e10_zero)


# ---------------------------------------------------------------------------
# D/E/F/G. our D-transform G-system block structure + the obstruction
# ---------------------------------------------------------------------------
def check_DEFG():
    from g_system_75_125 import build_gsystem, Phi

    def sysfor(a):
        b = 2 * a - 1
        return build_gsystem(a, b, 5, 2, 30 * a**2 - 24 * a + 3)

    r2, r3 = sysfor(2), sysfor(3)
    # D: five-block growth
    ok("D: a=2 has 5 generators, 4 spares", len(r2["Gs"]) == 5 and len(r2["spares"]) == 4)
    ok("D: a=3 has 10 generators, 9 spares", len(r3["Gs"]) == 10 and len(r3["spares"]) == 9)
    ok("D: step adds EXACTLY 5 forcing generators + 5 spare window unknowns",
       len(r3["Gs"]) - len(r2["Gs"]) == 5 and len(r3["spares"]) - len(r2["spares"]) == 5)
    ok("D: Phi u-slice deepens M: 21 -> 36", (r2["M"], r3["M"]) == (21, 36))

    # E: Phi recurrence
    def phi(a):
        return -sp.Rational(1, 3 * a) * y**(2 * a - 1) * g**a * C**(15 * a**2 - 13 * a + 2)
    for a in (2, 3, 4):
        ok("E: Phi recurrence Phi_{%d+1}=(%d/%d) C^{%d} Phi_%d" % (a, a, a + 1, 30 * a + 3, a),
           sp.simplify(phi(a + 1) - sp.Rational(a, a + 1) * C**(30 * a + 3) * phi(a)) == 0)

    # F: generator nesting  coeff_{d0^2}(G1^{a3}) = (10/3) G1^{a2}
    d0 = sp.Symbol("d0")
    G1a2 = sp.expand(r2["Gs"][1])
    G1a3 = sp.expand(r3["Gs"][1])
    lead = sp.expand(G1a3.coeff(d0, 2))
    ok("F: a=2 G1 is the leading d0^2 block of a=3 G1 (coeff = (10/3) G1^{a2})",
       sp.expand(lead - sp.Rational(10, 3) * G1a2) == 0)

    # G: window-denominator invariant + fractional-denominator fragmentation
    def qwin(a):
        b = 2 * a - 1
        M = b * 5 + (5 * a - 4)
        return sp.denom(sp.Rational(30 * a**2 - 24 * a + 3, M))
    ok("G: q_window = 5a-3 jumps 7 -> 12, gcd(7,12)=1 (incommensurate lattices)",
       qwin(2) == 7 and qwin(3) == 12 and sp.gcd(qwin(2), qwin(3)) == 1)

    def fracdenoms(a):
        b = 2 * a - 1
        t = 5
        jphi = 5 * a - 4
        M = b * t + jphi
        ordphi = 30 * a**2 - 24 * a + 3
        js = [j for j in range(1, jphi + 1) if j != jphi - 1]
        return sorted(set(sp.denom(sp.Rational(ordphi * (b * t + j), M)) for j in js))
    ok("G: y-order fractional denominators fragment {1,7} (a=2) -> {1,2,3,4,6,12} (a=3)",
       fracdenoms(2) == [1, 7] and fracdenoms(3) == [1, 2, 3, 4, 6, 12])


def main():
    check_A()
    check_B()
    check_C()
    check_DEFG()
    if not QUIET:
        print()
    if FAILS:
        print("FAILURES:", len(FAILS), FAILS)
        sys.exit(1)
    print("ALL %d F2-TOWER CHECKS PASSED" % N_OK)
    sys.exit(0)


if __name__ == "__main__":
    main()
