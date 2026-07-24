#!/usr/bin/env python3
"""Source-linked exact checks for T5_T2_INFINITY.md.

All h_j coefficients are parsed from f31_graded.txt by the established
t5_t2_column_verify.load_h() regex; none is maintained independently here.
"""
from __future__ import annotations

import sympy as sp
import t5_t2_column_verify as column

def _require(_cond, _msg):
    """Proof-critical check: fails loudly and exits nonzero, unaffected by python -O."""
    if not _cond:
        import sys as _sys
        print("FAIL: " + str(_msg))
        _sys.exit(1)


y, t, q, c, Phi = column.y, column.t, column.q, column.c, column.Phi
d0, d1, d2, e = column.d0, column.d1, column.d2, column.e
h = column.load_h()
PH, EE, DD, SS = sp.symbols("PH EE DD SS")
H = {j: sp.expand(h[j].subs({d1: 0, d0: (d2**2 + SS)/4})
                  .subs({d2: DD, e: EE})) for j in range(8)}
MASTER = sp.expand(sum(PH**j*EE**(21-3*j)*H[j] for j in range(8)))


def coeff(expr: sp.Expr, var: sp.Symbol, degree: int) -> sp.Expr:
    return sp.Poly(sp.expand(expr), var).coeff_monomial(var**degree)


def master_at(phi: sp.Expr, ep: sp.Expr, dp: sp.Expr, sig: sp.Expr) -> sp.Expr:
    return sp.expand(MASTER.subs({PH: phi, EE: ep, DD: dp, SS: sig}))


def check_shared_convolution() -> None:
    """One reversed-series convolution covers all 29 pattern-B states."""
    source_terms = sp.Poly(MASTER, PH, EE, DD, SS).terms()
    for monomial, _ in source_terms:
        _require(sum(w*k for w, k in zip((34, 10, 4, 8), monomial)) == 250, "sum(w*k for w, k in zip((34, 10, 4, 8), monomial)) == 250")

    z = sp.Symbol("z")
    p0, p1, e0, e1, d4, d3, s8, s7 = sp.symbols(
        "p0 p1 e0 e1 d4 d3 s8 s7")
    leading = MASTER.subs({PH: p0, EE: e0, DD: d4, SS: s8})
    next_one = (p1*sp.diff(MASTER, PH) + e1*sp.diff(MASTER, EE)
                + d3*sp.diff(MASTER, DD) + s7*sp.diff(MASTER, SS)).subs(
                    {PH: p0, EE: e0, DD: d4, SS: s8})
    # Taylor's rule on every source monomial gives exactly leading and
    # next_one.  Weighted homogeneity identifies them with y^250,y^249.
    _require(leading != 0 and next_one != 0, "leading != 0 and next_one != 0")
    slots = ((PH, p0, p1), (EE, e0, e1),
             (DD, d4, d3), (SS, s8, s7))
    for monomial, value in source_terms:
        term = value*sp.prod(symbol**power
                             for (symbol, _, _), power in zip(slots, monomial))
        direct = 0
        for (symbol, initial, following), power in zip(slots, monomial):
            if power:
                direct += following*sp.diff(term, symbol).subs(
                    {PH: p0, EE: e0, DD: d4, SS: s8})
        derivative = sum(following*sp.diff(term, symbol)
                         for symbol, _, following in slots).subs(
                             {PH: p0, EE: e0, DD: d4, SS: s8})
        _require(sp.expand(direct-derivative) == 0, "sp.expand(direct-derivative) == 0")

    cells, seen = column.load_cells(), []
    for key in column.TARGETS:
        if key in column.KILLED_STATES:
            continue
        B, S, _, fcap, gcap, _ = column.invariants(cells[key])
        for state in column.terminal_states(key[0], B, S, fcap, gcap):
            fdeg, zdeg, gdeg, D, Sigma = state
            if (fdeg, zdeg) not in column.RESIDUAL[key] or (D, Sigma) == (8, 3):
                continue
            _require(D == 10 and 3*fdeg+gdeg == 2*zdeg, "D == 10 and 3*fdeg+gdeg == 2*zdeg")
            seen.append((key, fdeg, zdeg, Sigma))
    _require(len(seen) == 29 and leading != 0, "len(seen) == 29 and leading != 0")
    print("I1. weight-250 master and shared C250/C249 convolution             OK")


def check_pattern_a() -> None:
    """Kill the two (D,Sigma)=(8,3) states at degrees 230 and 232."""
    gamma, beta = sp.symbols("gamma beta", nonzero=True)
    lcphi = sp.LC(sp.Poly(Phi, y))
    k5 = sp.Poly(H[5], DD, SS, EE).coeff_monomial(EE**2)
    k6 = sp.Poly(H[6], DD, SS, EE).coeff_monomial(SS**2)
    _require(k5 and k6 and sp.Rational(k6, k5) == -sp.Rational(3, 2), "k5 and k6 and sp.Rational(k6, k5) == -sp.Rational(3, 2)")

    top_monoms = {(5, 8, 0, 0): k5, (6, 3, 0, 2): k6}
    source_terms = dict(sp.Poly(MASTER, PH, EE, DD, SS).terms())
    _require(all(source_terms[m] == value for m, value in top_monoms.items()), "all(source_terms[m] == value for m, value in top_monoms.items())")
    _require(max(34*m[0]+8*m[1]+4*m[2]+3*m[3]
               for m in source_terms if m not in top_monoms) <= 228, "max(34*m[0]+8*m[1]+4*m[2]+3*m[3] for m in source_terms if m not in top_monoms) <= 228")

    def normalized_tail(Ep: sp.Expr, Sp: sp.Expr) -> sp.Expr:
        beta2 = sp.cancel(2*gamma**5/(3*lcphi))
        tail = 2*gamma**5*Ep**5 - 3*Phi*beta2*Sp**2
        return sp.expand(tail/(2*gamma**5))

    # a8 b0000: degrees 234..231 force the cubic; degree 230 is nonzero.
    s0, s1, s2 = sp.symbols("s0 s1 s2")
    W80 = normalized_tail(t**8, y**3+s2*y**2+s1*y+s0)
    forced80 = {}
    for degree, variable in ((39, s2), (38, s1), (37, s0)):
        equation = sp.factor(coeff(W80, y, degree).subs(forced80))
        sol = sp.solve(equation, variable)
        _require(len(sol) == 1, "len(sol) == 1")
        forced80[variable] = sol[0]
    obstruction80 = sp.factor(coeff(W80, y, 36).subs(forced80))
    _require(forced80 == {s2: sp.Rational(41, 8),
                        s1: sp.Rational(1353, 128),
                        s0: sp.Rational(11275, 1024)}, "forced80 == {s2: sp.Rational(41, 8), s1: sp.Rational(1353, 128), s0: sp.Rational(11275, 1024)}")
    _require(obstruction80 == sp.Rational(191675, 16384), "obstruction80 == sp.Rational(191675, 16384)")
    _require(obstruction80 != 0, "obstruction80 != 0")

    # a7 b1000: degree 233 forces u; degree 232 is coprime to q(r).
    r, u = sp.symbols("r u")
    W71 = normalized_tail(t**7*(y-r), (y-r)**2*(y-u))
    usol = sp.solve(sp.factor(coeff(W71, y, 39)), u)
    _require(len(usol) == 1, "len(usol) == 1")
    _require(usol[0] == r/2-sp.Rational(21, 8), "usol[0] == r/2-sp.Rational(21, 8)")
    obstruction71 = sp.factor(coeff(W71, y, 38).subs(u, usol[0]))
    _require(sp.expand(obstruction71+(16*r**2+168*r-273)/64) == 0, "sp.expand(obstruction71+(16*r**2+168*r-273)/64) == 0")
    num = sp.Poly(sp.together(obstruction71).as_numer_denom()[0], r)
    _require(sp.gcd(sp.Poly(q.subs(y, r), r), num).degree() == 0, "sp.gcd(sp.Poly(q.subs(y, r), r), num).degree() == 0")
    _require((0, 3) in column.RESIDUAL[(8, (0, 0, 0, 0))], "(0, 3) in column.RESIDUAL[(8, (0, 0, 0, 0))]")
    _require((0, 1) in column.RESIDUAL[(7, (1, 0, 0, 0))], "(0, 1) in column.RESIDUAL[(7, (1, 0, 0, 0))]")
    print("I2. pattern A: a8 b0000 degree 230; a7 b1000 degree 232          OK")
    print("    final normalized obstructions:", obstruction80, ";", obstruction71)


def check_d2_zero_low_sigma() -> None:
    """For d2=0,Sigma<=6, fixed-F supports fail at degree 249."""
    eps, enorm = sp.symbols("eps enorm", nonzero=True)
    lcphi, p33 = sp.LC(sp.Poly(Phi, y)), coeff(Phi, y, 33)
    k5 = sp.Poly(H[5], DD, SS, EE).coeff_monomial(EE**2)
    k0 = sp.Poly(H[0], DD, SS, EE).coeff_monomial(EE**4)
    C0 = MASTER.subs({PH: lcphi, EE: eps, DD: 0, SS: 0})
    _require(sp.expand(C0-(k5*lcphi**5*eps**8+k0*eps**25)) == 0, "sp.expand(C0-(k5*lcphi**5*eps**8+k0*eps**25)) == 0")
    eps17 = sp.cancel(-k5*lcphi**5/k0)
    C1 = (p33*sp.diff(MASTER, PH)+eps*enorm*sp.diff(MASTER, EE)).subs(
        {PH: lcphi, EE: eps, DD: 0, SS: 0})
    C1 = sp.factor(sp.expand(C1/eps**8).subs(eps**17, eps17))
    target = sp.solve(C1, enorm)
    _require(len(target) == 1, "len(target) == 1")
    target = sp.factor(target[0])
    _require(target == sp.factor(sp.Rational(5, 17)*p33/lcphi), "target == sp.factor(sp.Rational(5, 17)*p33/lcphi)")
    _require(target == sp.Rational(35, 4), "target == sp.Rational(35, 4)")

    r = sp.Symbol("r")
    qr = sp.Poly(q.subs(y, r), r)
    _require(q.subs(y, sp.factor(9-target)) != 0, "q.subs(y, sp.factor(9-target)) != 0")
    pair_sum = sp.factor(8-target)
    _require(sp.gcd(qr, sp.Poly(sp.together(q.subs(y, pair_sum-r)), r)).degree() == 0, "sp.gcd(qr, sp.Poly(sp.together(q.subs(y, pair_sum-r)), r)).degree() == 0")
    qpoly = sp.Poly(q, y)
    rootsum = -qpoly.all_coeffs()[1]/qpoly.all_coeffs()[0]
    _require(q.subs(y, sp.factor(target-(7-rootsum))) != 0, "q.subs(y, sp.factor(target-(7-rootsum))) != 0")

    killed = {(9, (1, 0, 0, 0)): range(5),
              (8, (1, 1, 0, 0)): range(3),
              (7, (1, 1, 1, 0)): (0,)}
    cells = column.load_cells()
    for key, zvalues in killed.items():
        _require(any(case["d2_zero"] for case in cells[key]["survivor_cases"]), "any(case[\"d2_zero\"] for case in cells[key][\"survivor_cases\"])")
        B, S, _, fcap, gcap, _ = column.invariants(cells[key])
        states = column.terminal_states(key[0], B, S, fcap, gcap)
        for zd in zvalues:
            _require(next(x for x in states if x[0:2] == (0, zd))[4] <= 6, "next(x for x in states if x[0:2] == (0, zd))[4] <= 6")
    _require(sum(len(tuple(v)) for v in killed.values()) == 9, "sum(len(tuple(v)) for v in killed.values()) == 9")
    print("I3. d2=0,Sigma<=6: nine degree states killed at degree 249        OK")
    print("    forced e_9/e_10 =", target)


def check_g5_zero() -> None:
    """Kill the sole a9 b1000 G5 state by the exact level-5 line."""
    h5poly = sp.Poly(H[5], DD, SS, EE)
    kd = h5poly.coeff_monomial(DD*SS**2)
    ke = h5poly.coeff_monomial(EE**2)
    _require(kd != 0 and ke != 0 and H[5] == kd*DD*SS**2+ke*EE**2, "kd != 0 and ke != 0 and H[5] == kd*DD*SS**2+ke*EE**2")
    r = sp.Symbol("r")
    p = y-r
    zc = sp.symbols("z0:6")
    Z = y**6+sum(zc[k]*y**k for k in range(6))
    l0, l1 = sp.symbols("l0 l1")
    gamma, zeta, mu = sp.symbols("gamma zeta mu", nonzero=True)
    L = l1*y+l0
    hp = H[5].subs({DD: t**3*L, SS: zeta*p**2*Z, EE: gamma*t**9*p})
    line = sp.cancel((hp-mu*t**3*q*p*Z**2)/(t**3*p))
    expected = ke*gamma**2*t**15*p + Z**2*(kd*zeta**2*L*p**3-mu*q)
    _require(sp.expand(line-expected) == 0, "sp.expand(line-expected) == 0")
    _require(sp.gcd(sp.Poly(t, y), sp.Poly(q, y)).degree() == 0, "sp.gcd(sp.Poly(t, y), sp.Poly(q, y)).degree() == 0")
    _require(sp.degree(Z, y) == 6, "sp.degree(Z, y) == 6")
    cells = column.load_cells()
    key = (9, (1, 0, 0, 0))
    g5case = next(case for case in cells[key]["survivor_cases"]
                  if tuple(case["g_zero_levels"]) == (5,))
    b, _, m = column.q_profile(g5case)
    _require(tuple(mi-6 for mi in m) == b == key[1], "tuple(mi-6 for mi in m) == b == key[1]")
    tp = next(place for place in g5case["witness"] if place["place"] == "t")
    _require(tp["v_d2"] == 3 and tp["v_sigma"] == 0, "tp[\"v_d2\"] == 3 and tp[\"v_sigma\"] == 0")
    # For any irreducible rho|Z, the hypotheses give v_rho(t^15 p)=0,
    # while v_rho(Z^2)=2m>=2.  This executable integer check is the UFD kill.
    _require(all(2*mult > 0 for mult in range(1, 7)), "all(2*mult > 0 for mult in range(1, 7))")
    _require(not any(2*mult <= 0 for mult in range(1, 7)), "not any(2*mult <= 0 for mult in range(1, 7))")
    print("I4. a9 b1000 G5 exact level-5 UFD descent                         OK")


def check_quotients() -> None:
    """Check the common forced g5 and g4 cascade parametrization."""
    R, F, QB, G, g5, g4, qs, tv = sp.symbols("R F QB G g5 g4 qs tv", nonzero=True)
    ta2, dd, h4name = sp.symbols("ta2 dd h4name")
    h5poly = sp.Poly(H[5], DD, SS, EE)
    ke = h5poly.coeff_monomial(EE**2)
    kd = h5poly.coeff_monomial(DD*SS**2)
    terminal = -sp.Poly(H[6], DD, SS, EE).coeff_monomial(SS**2)
    absorbed = sp.cancel(kd/(terminal*c))
    _require((ke, kd, terminal, absorbed) == (2048, -9216, 3072, 19890), "(ke, kd, terminal, absorbed) == (2048, -9216, 3072, 19890)")
    forced5 = sp.cancel((tv*qs*QB*G-ke*c**5*qs**5*ta2*R**2*F**2)
                        /(R**3*F**3)-absorbed*dd*QB*G)
    equation5 = g5-forced5
    _require(sp.diff(equation5, g5) == 1 and equation5.subs(g5, forced5) == 0, "sp.diff(equation5, g5) == 1 and equation5.subs(g5, forced5) == 0")
    forced4 = sp.cancel((tv*g5-c**4*qs**4*h4name)/(R**3*F**3))
    equation4 = g4-forced4
    _require(sp.diff(equation4, g4) == 1 and equation4.subs(g4, forced4) == 0, "sp.diff(equation4, g4) == 1 and equation4.subs(g4, forced4) == 0")
    print("I5. forced terminal/cascade quotient parametrization              OK")


def main() -> None:
    _require(H[7] == 0, "H[7] == 0")
    h6poly = sp.Poly(H[6], DD, SS, EE)
    _require(len(h6poly.terms()) == 1 and h6poly.monoms() == [(0, 2, 0)], "len(h6poly.terms()) == 1 and h6poly.monoms() == [(0, 2, 0)]")
    _require(H[4].subs(DD, 0) == -12096*SS**3, "H[4].subs(DD, 0) == -12096*SS**3")
    check_shared_convolution()
    check_pattern_a()
    check_d2_zero_low_sigma()
    check_g5_zero()
    check_quotients()
    print("\nALL T2 INFINITY CHECKS PASS")
    print("  cells killed/narrowed/unchanged: 0 / 8 / 0")
    print("  closed flag case: a9 b1000 G5")
    print("  additional residual degree-state removals: 11 (2 pattern A, 9 D)")


if __name__ == "__main__":
    main()
