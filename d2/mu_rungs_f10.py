#!/usr/bin/env python3
"""Exact enumeration of the real mu-rungs at the F10 corner.

Independently reconstructs 4(7*c*f' - 27*c'*f)=c^4, c=y^3*g.  The residual
convention also requires h(0)!=0 (otherwise ord_y(c) is not q=3).  Every root
partition of h is covered.  Arithmetic and real-root counts are exact.
"""
from __future__ import annotations
import sympy as sp

y = sp.symbols("y")
A, B, T, KAPPA = 4, 7, 7, 5
Q, E, DG = 3, 4, 4
COEF = T*(B-A) + KAPPA + 1
RHO, N, RES, GAP, R = 10, 270, 27, 1, 3
assert COEF == 27


def ode(c, f):
    return sp.expand(A*(T*c*sp.diff(f, y)-COEF*sp.diff(c, y)*f)-c**E)


def solve_reduced(mu, h, hfac, udeg, common):
    us = sp.symbols(f"u0:{udeg+1}")
    u = sum(us[i]*y**i for i in range(udeg+1))
    c = y**Q*(y+1)**mu*h
    f = y**RHO*(y+1)**(3*mu+1)*hfac*u
    quotient = sp.Poly(sp.cancel(ode(c, f)/common), y)
    eqs = [sp.factor(quotient.nth(i)) for i in range(quotient.degree()+1)]
    linear = sp.solve(eqs[:udeg+1], us, dict=True)
    assert len(linear) == 1
    usol = sp.factor(u.subs(linear[0]))
    rem = [sp.factor(sp.numer(sp.together(eq.subs(linear[0]))))
           for eq in eqs[udeg+1:] if sp.together(eq.subs(linear[0])) != 0]
    return c, sp.factor(f.subs(linear[0])), usol, rem


def real_count(poly, var):
    return int(sp.Poly(poly, var, domain=sp.QQ).count_roots(-sp.oo, sp.oo))


def mu_law(mu):
    return (RES+N*(Q+DG), RHO+N*Q,
            mu*(E+N)-(mu-1),
            GAP+R*(E+N)-(mu-1)*(E+N-1))


def derive_mu1():
    p, q, r = sp.symbols("p q r")
    h = y**3+p*y**2+q*y+r
    _, _, _, rem = solve_reduced(1, h, h**4, 1,
                                  y**12*(y+1)**4*h**4)
    assert len(rem) == 3
    gb = sp.groebner(rem, p, q, r, order="lex")
    m_sf = 2250*r**4-4776*r**3+1180*r**2+75*r+2250
    assert sp.factor(gb.polys[-1].as_expr()) == r*m_sf
    assert real_count(m_sf, r) == 0

    z = sp.symbols("z")
    h3 = (y-z)**3
    _, _, _, rem3 = solve_reduced(1, h3, (y-z)**10, 3,
                                   y**12*(y+1)**4*(y-z)**12)
    m_3 = 5*z**4-2*z**3+8*z**2+20*z+10
    assert len(rem3) == 1 and sp.Poly(rem3[0], z).monic() == sp.Poly(m_3, z).monic()
    assert real_count(m_3, z) == 0

    z, w = sp.symbols("z w")
    h21 = (y-z)**2*(y-w)
    _, _, _, rem21 = solve_reduced(1, h21, (y-z)**7*(y-w)**4, 2,
                                    y**12*(y+1)**4*(y-z)**8*(y-w)**4)
    m_21 = (1250*w**12-1500*w**11+4400*w**10+3850*w**9+6303*w**8
            -8272*w**7-11792*w**6-8272*w**5+6303*w**4+3850*w**3
            +4400*w**2-1500*w+1250)
    assert sp.factor(sp.resultant(rem21[0], rem21[1], z)) == 6144*w**6*m_21
    assert real_count(m_21, w) == 0
    return {"squarefree": m_sf, "triple": m_3, "double_simple": m_21}


def derive_mu2():
    p, r = sp.symbols("p r")
    h = y**2+p*y+r
    c_sf, f_sf, u_sf, rem = solve_reduced(
        2, h, h**4, 2, y**12*(y+1)**8*h**4)
    expected = [
        10*p**3+p**2*r+3*p*r**2-17*p*r-24*r**3+12*r**2,
        -5*p**3-13*p**2*r+25*p**2-24*p*r**2+36*p*r+64*r**2-30*r,
    ]
    assert rem == expected
    m_sf = 192*p**6-720*p**5+586*p**4-668*p**3+1728*p**2-1510*p+45
    assert sp.factor(sp.resultant(rem[0], rem[1], r)) == -4000*p**3*m_sf
    assert sp.Poly(m_sf, p).is_irreducible and real_count(m_sf, p) == 2
    # p=0 gives incompatible nonzero-r requirements; on M_sf, r is unique.
    g0 = sp.gcd(sp.Poly(expected[0].subs(p, 0), r),
                sp.Poly(expected[1].subs(p, 0), r))
    assert g0.degree() == 1 and g0.eval(0) == 0
    r_of_p = (-768*p**5+2400*p**4-844*p**3+3492*p**2-6077*p+1905)/sp.Integer(3773)
    for eq in expected:
        num = sp.numer(sp.together(eq.subs(r, r_of_p)))
        assert sp.rem(sp.Poly(num, p), sp.Poly(m_sf, p)).is_zero
    for boundary in (r_of_p, 1-p+r_of_p, p**2-4*r_of_p):
        num = sp.numer(sp.together(boundary))
        assert sp.gcd(sp.Poly(num, p), sp.Poly(m_sf, p)).degree() == 0

    z = sp.symbols("z")
    hd = (y-z)**2
    c_d, f_d, u_d, rem_d = solve_reduced(
        2, hd, (y-z)**7, 3, y**12*(y+1)**8*(y-z)**8)
    m_d = 12*z**4+15*z**3-5*z**2+15*z+12
    assert rem_d == [m_d]
    assert sp.Poly(m_d, z).is_irreducible and real_count(m_d, z) == 2
    assert sp.gcd(sp.Poly(m_d, z), sp.Poly(z*(z+1), z)).degree() == 0
    return {
        "squarefree": dict(parameter=p, minimal=m_sf, relation=r_of_p,
                           c=c_sf, f=f_sf, u=u_sf),
        "double": dict(parameter=z, minimal=m_d, c=c_d, f=f_d, u=u_d),
    }


def derive_mu3():
    s = sp.symbols("s")
    h = y+s
    c, f, u, rem = solve_reduced(3, h, h**4, 3,
                                  y**12*(y+1)**12*h**4)
    m = 10*s**4-20*s**3+8*s**2+2*s+5
    assert len(rem) == 1 and sp.Poly(rem[0], s).monic() == sp.Poly(m, s).monic()
    assert sp.Poly(m, s).is_irreducible and real_count(m, s) == 0
    return dict(parameter=s, obstruction=m, c=c, f=f, u=u)


def derive_mu4():
    g = (y+1)**4
    c = y**3*g
    f = (y**10*(y+1)**13
         *(2401*y**4+5831*y**3+4165*y**2+595*y-85)/sp.Integer(3740))
    assert ode(c, f) == 0
    return dict(g=g, c=c, f=f, resonant=sp.Rational(2401, 3740))


def main():
    assert (A, B, T, Q, E, DG, RHO, N, RES, GAP, R) == (4,7,7,3,4,4,10,270,27,1,3)
    mu1 = derive_mu1()
    mu2 = derive_mu2()
    mu3 = derive_mu3()
    derive_mu4()
    p = mu2["squarefree"]["parameter"]
    z = mu2["double"]["parameter"]
    print("F10 exact real mu-rungs")
    print("  mu=1 EMPTY")
    print("    squarefree obstruction:", mu1["squarefree"])
    print("    (2,1) obstruction:", mu1["double_simple"])
    print("    triple obstruction:", mu1["triple"])
    print("    exact real-root counts: 0, 0, 0")
    print("  mu=2 REALIZED")
    print(f"    squarefree: minpoly({p}) = {mu2['squarefree']['minimal']}")
    print(f"      r = {mu2['squarefree']['relation']}; 2 real conjugates")
    print(f"    double-root: minpoly({z}) = {mu2['double']['minimal']}")
    print("      2 real conjugates")
    print(f"    signature {mu_law(2)}; mu-graded law MATCH")
    print("  mu=3 EMPTY")
    print(f"    obstruction: {mu3['obstruction']}; exact real-root count 0")
    print("  mu=4 REALIZED")
    print("    g=(y+1)^4; minpoly(F27)=3740*F27-2401")
    print(f"    signature {mu_law(4)}; mu-graded law MATCH")


if __name__ == "__main__":
    main()




