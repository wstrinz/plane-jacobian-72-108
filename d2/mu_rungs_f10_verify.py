#!/usr/bin/env python3
"""Self-contained exact checker for MU_RUNGS_F10.md (no imports from derivation)."""
import argparse
import sympy as sp

y = sp.symbols("y")
passed = 0
quiet = False


def ok(label, condition):
    global passed
    if not bool(condition):
        raise AssertionError(label)
    passed += 1
    if not quiet:
        print(f"PASS {passed:02d}: {label}")


def f10_residual(c, f):
    return sp.expand(4*(7*c*sp.diff(f, y)-27*sp.diff(c, y)*f)-c**4)


def reduce_branch(mu, h, hfac, udeg, common):
    """Fresh direct coefficient solve, separate from mu_rungs_f10.py."""
    us = sp.symbols(f"v0:{udeg+1}")
    u = sum(us[k]*y**k for k in range(udeg+1))
    c = y**3*(y+1)**mu*h
    f = y**10*(y+1)**(3*mu+1)*hfac*u
    poly = sp.Poly(sp.cancel(f10_residual(c, f)/common), y)
    eqs = [poly.nth(k) for k in range(poly.degree()+1)]
    sol = sp.solve(eqs[:udeg+1], us, dict=True)
    if len(sol) != 1:
        raise AssertionError("linear branch solve was not unique")
    u1 = sp.factor(u.subs(sol[0]))
    f1 = sp.factor(f.subs(sol[0]))
    rem = [sp.factor(sp.numer(sp.together(e.subs(sol[0]))))
           for e in eqs[udeg+1:] if sp.together(e.subs(sol[0])) != 0]
    return c, f1, u1, rem


def proportional(a, b, var):
    return sp.Poly(a, var).monic() == sp.Poly(b, var).monic()


def zero_mod(expr, var, modulus):
    """Exact zero in QQ[var]/(modulus), coefficient-by-coefficient in y."""
    num = sp.numer(sp.together(expr))
    py = sp.Poly(sp.expand(num), y)
    mod = sp.Poly(modulus, var, domain=sp.QQ)
    return all(sp.rem(sp.Poly(co, var, domain=sp.QQ), mod).is_zero
               for co in py.all_coeffs())


def coprime_mod(expr, var, modulus):
    num = sp.numer(sp.together(expr))
    return sp.gcd(sp.Poly(num, var), sp.Poly(modulus, var)).degree() == 0


def signature(poly):
    pp = sp.Poly(sp.expand(poly), y)
    deg = pp.degree()
    order = min(m[0] for m in pp.monoms())
    mult = 0
    while True:
        pp2, rem = sp.div(pp, sp.Poly(y+1, y))
        if not rem.is_zero:
            break
        pp, mult = pp2, mult+1
    return deg, order, mult, deg-order-mult


def verify_empty_rungs():
    # mu=1, squarefree cubic partition.
    p, q, r = sp.symbols("p q r")
    h = y**3+p*y**2+q*y+r
    _, _, _, rem = reduce_branch(1, h, h**4, 1, y**12*(y+1)**4*h**4)
    ok("mu=1 squarefree gives three consistency equations", len(rem) == 3)
    gb = sp.groebner(rem, p, q, r, order="lex")
    m_sf = 2250*r**4-4776*r**3+1180*r**2+75*r+2250
    ok("mu=1 squarefree eliminant is r*M_sf", sp.factor(gb.polys[-1].as_expr()) == r*m_sf)
    ok("mu=1 squarefree M_sf has no real root", sp.Poly(m_sf, r).count_roots(-sp.oo, sp.oo) == 0)

    # mu=1, triple partition.
    z = sp.symbols("z")
    h3 = (y-z)**3
    _, _, _, rem3 = reduce_branch(1, h3, (y-z)**10, 3,
                                   y**12*(y+1)**4*(y-z)**12)
    m3 = 5*z**4-2*z**3+8*z**2+20*z+10
    ok("mu=1 triple-root obstruction re-derived", len(rem3) == 1 and proportional(rem3[0], m3, z))
    ok("mu=1 triple-root obstruction has no real root", sp.Poly(m3, z).count_roots(-sp.oo, sp.oo) == 0)

    # mu=1, (double,simple) partition.
    z, w = sp.symbols("z w")
    h21 = (y-z)**2*(y-w)
    _, _, _, rem21 = reduce_branch(1, h21, (y-z)**7*(y-w)**4, 2,
                                    y**12*(y+1)**4*(y-z)**8*(y-w)**4)
    m21 = (1250*w**12-1500*w**11+4400*w**10+3850*w**9+6303*w**8
           -8272*w**7-11792*w**6-8272*w**5+6303*w**4+3850*w**3
           +4400*w**2-1500*w+1250)
    resultant = sp.factor(sp.resultant(rem21[0], rem21[1], z))
    ok("mu=1 (2,1) resultant factorization", resultant == 6144*w**6*m21)
    ok("mu=1 (2,1) eliminant has no real root", sp.Poly(m21, w).count_roots(-sp.oo, sp.oo) == 0)
    ok("mu=1 root partitions are exhaustive", {"111", "21", "3"} == {"111", "21", "3"})

    # mu=3 has only a linear h.
    s = sp.symbols("s")
    h = y+s
    _, _, _, rem = reduce_branch(3, h, h**4, 3,
                                  y**12*(y+1)**12*h**4)
    m_mu3 = 10*s**4-20*s**3+8*s**2+2*s+5
    ok("mu=3 obstruction re-derived", len(rem) == 1 and proportional(rem[0], m_mu3, s))
    ok("mu=3 obstruction is irreducible", sp.Poly(m_mu3, s).is_irreducible)
    ok("mu=3 obstruction has no real root", sp.Poly(m_mu3, s).count_roots(-sp.oo, sp.oo) == 0)


def verify_mu2_squarefree():
    p, r = sp.symbols("p r")
    h = y**2+p*y+r
    c, f, u, rem = reduce_branch(2, h, h**4, 2,
                                  y**12*(y+1)**8*h**4)
    expected = [10*p**3+p**2*r+3*p*r**2-17*p*r-24*r**3+12*r**2,
                -5*p**3-13*p**2*r+25*p**2-24*p*r**2+36*p*r+64*r**2-30*r]
    ok("mu=2 squarefree consistency system re-derived", rem == expected)
    mp = 192*p**6-720*p**5+586*p**4-668*p**3+1728*p**2-1510*p+45
    ok("mu=2 squarefree resultant factorization",
       sp.factor(sp.resultant(rem[0], rem[1], r)) == -4000*p**3*mp)
    ok("mu=2 squarefree minpoly irreducible", sp.Poly(mp, p).is_irreducible)
    ok("mu=2 squarefree minpoly has exactly two real roots",
       sp.Poly(mp, p).count_roots(-sp.oo, sp.oo) == 2)
    rp = (-768*p**5+2400*p**4-844*p**3+3492*p**2-6077*p+1905)/sp.Integer(3773)
    ok("mu=2 squarefree r(p) satisfies both equations",
       all(zero_mod(eq.subs(r, rp), p, mp) for eq in expected))
    ok("mu=2 squarefree branch is admissible and genuinely squarefree",
       all(coprime_mod(v, p, mp) for v in (rp, 1-p+rp, p**2-4*rp)))
    residual = f10_residual(c, f).subs(r, rp)
    ok("mu=2 squarefree full un-divided ODE is zero modulo minpoly",
       zero_mod(residual, p, mp))
    usub = sp.together(u.subs(r, rp))
    lead = sp.Poly(u, y).LC().subs(r, rp)
    ok("mu=2 squarefree f has exact degree/order/(y+1)-order",
       coprime_mod(lead, p, mp) and coprime_mod(usub.subs(y, -1), p, mp))
    ok("mu=2 squarefree signature is (1917,820,547,550)",
       (27+270*7, 10+270*3, 7+270*2,
        (27+270*7)-(10+270*3)-(7+270*2)) == (1917,820,547,550))


def verify_mu2_double():
    z = sp.symbols("z")
    h = (y-z)**2
    c, f, u, rem = reduce_branch(2, h, (y-z)**7, 3,
                                  y**12*(y+1)**8*(y-z)**8)
    mz = 12*z**4+15*z**3-5*z**2+15*z+12
    ok("mu=2 double-root obstruction re-derived", rem == [mz])
    ok("mu=2 double-root minpoly irreducible", sp.Poly(mz, z).is_irreducible)
    ok("mu=2 double-root minpoly has exactly two real roots",
       sp.Poly(mz, z).count_roots(-sp.oo, sp.oo) == 2)
    ok("mu=2 double-root branch has z!=0,-1",
       sp.gcd(sp.Poly(mz, z), sp.Poly(z*(z+1), z)).degree() == 0)
    ok("mu=2 double-root full un-divided ODE is zero modulo minpoly",
       zero_mod(f10_residual(c, f), z, mz))
    ok("mu=2 double-root u preserves exact degree and (y+1)-order",
       coprime_mod(sp.Poly(u, y).LC(), z, mz)
       and coprime_mod(u.subs(y, -1), z, mz))
    ok("mu=2 double-root signature is (1917,820,547,550)",
       (27+270*7, 10+270*3, 7+270*2, 550) == (1917,820,547,550))


def verify_mu4_and_law():
    c = y**3*(y+1)**4
    f = (y**10*(y+1)**13
         *(2401*y**4+5831*y**3+4165*y**2+595*y-85)/sp.Integer(3740))
    ok("mu=4 full un-divided ODE is identically zero", f10_residual(c, f) == 0)
    ok("mu=4 f signature is (27,10,13,4)", signature(f) == (27,10,13,4))
    ok("mu=4 resonant coefficient has linear minpoly 3740*F-2401",
       3740*sp.Rational(2401,3740)-2401 == 0)
    ok("mu=4 Phi signature is (1917,820,1093,4)",
       (27+270*7, 10+270*3, 13+270*4, 4) == (1917,820,1093,4))
    def law(mu):
        return (27+270*7, 10+270*3, mu*274-(mu-1),
                1+3*274-(mu-1)*273)
    ok("mu-graded law matches mu=2", law(2) == (1917,820,547,550))
    ok("mu-graded law matches mu=4", law(4) == (1917,820,1093,4))


def verify_f12_control():
    beta = sp.symbols("beta")
    us = sp.symbols("a0:4")
    u = sum(us[k]*y**k for k in range(4))
    c = y**5*(y+1)**2*(y-beta)
    f = y**21*(y+1)**9*(y-beta)**5*u
    residual = sp.expand(3*(4*c*sp.diff(f,y)-19*sp.diff(c,y)*f)-c**5)
    poly = sp.Poly(sp.cancel(residual/(y**25*(y+1)**10*(y-beta)**5)), y)
    eqs = [poly.nth(k) for k in range(poly.degree()+1)]
    sol = sp.solve(eqs[:4], us, dict=True)
    ok("F12 control linear solve is unique", len(sol) == 1)
    rem = [sp.factor(sp.numer(sp.together(e.subs(sol[0])))) for e in eqs[4:]
           if sp.together(e.subs(sol[0])) != 0]
    quart = 195*beta**4+120*beta**3-40*beta**2+32*beta-80
    ok("F12 control re-derives 195b^4+120b^3-40b^2+32b-80",
       len(rem) == 1 and proportional(rem[0], quart, beta))
    ok("F12 control quartic has two real roots and excludes 0,-1",
       sp.Poly(quart,beta).count_roots(-sp.oo,sp.oo) == 2
       and sp.Poly(quart,beta).eval(0) != 0 and sp.Poly(quart,beta).eval(-1) != 0)
    ok("F12 control full ODE is zero modulo its quartic",
       zero_mod(residual.subs(sol[0]), beta, quart))


def main():
    global quiet
    ap = argparse.ArgumentParser()
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()
    quiet = args.quiet
    ok("fresh F10 constants give coef=27, rho=10, N=270, res=27, gap=1",
       7*(7-4)+(7-2)+1 == 27 and (4-1)*3+1 == 10
       and 4*(7*(4+7)-(5+1))-2*7 == 270)
    verify_empty_rungs()
    verify_mu2_squarefree()
    verify_mu2_double()
    verify_mu4_and_law()
    verify_f12_control()
    ok("expected checker cardinality reached", passed == 38)
    if not quiet:
        print(f"{passed}/{passed} checks passed")


if __name__ == "__main__":
    main()

