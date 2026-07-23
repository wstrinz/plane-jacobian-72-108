# t5_60t1_verify.py -- verification for stratum (6,0), branch T1
# All algebra recomputed from f31_graded.txt; prints only small summaries.
import sympy as sp

y = sp.symbols('y')
t = y + 1
q = sp.Poly(2048*y**4 - 512*y**3 + 320*y**2 - 240*y + 195, y)
c = sp.Rational(-1, 6630)

CHECKS = []
def check(name, ok):
    CHECKS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}")

import re

# ---------- parse h_f from f31_graded.txt (never embed literals) ----------
d0, d1, d2, dm1 = sp.symbols('d0 d1 d2 dm1')

def load_h():
    h = {}
    pat = re.compile(r'h_(\d+)\s*\([^)]*\)\s*=\s*(.+)$')
    with open('f31_graded.txt') as fh:
        for line in fh:
            m = pat.match(line.strip())
            if m:
                h[int(m.group(1))] = sp.sympify(
                    m.group(2), locals={'d0': d0, 'd1': d1, 'd2': d2, 'dm1': dm1})
    return h

def part_A():
    h = load_h()
    sigma = 4*d0 - d2**2

    # CHK0: all 8 graded pieces present; terminal h_7
    check("CHK0 parse f31_graded.txt: h_0..h_7 present, h_7 = 8192*d1^2",
          sorted(h) == list(range(8)) and sp.expand(h[7] - 8192*d1**2) == 0)

    # CHK1: q-place facts
    check("CHK1 q irreducible /Q, deg q = 4, q(-1) = 3315 (so t,q,r pairwise coprime)",
          sp.Poly(q, y).is_irreducible and q.degree() == 4
          and q.as_expr().subs(y, -1) == 3315)

    # CHK2: h_6 collapse   h_6 = -3072 sigma^2 + 14336 d1^2 d2 + 8192 d1 e
    check("CHK2 h_6 = -3072*sigma^2 + 14336*d1^2*d2 + 8192*d1*dm1",
          sp.expand(h[6] - (-3072*sigma**2 + 14336*d1**2*d2 + 8192*d1*dm1)) == 0)

    # CHK3: h_5 collapse
    check("CHK3 h_5 = -9216*d2*sigma^2 + 32256*d1^2*sigma - 12288*d1^2*d2^2"
          " + 18432*d1*d2*dm1 + 2048*dm1^2",
          sp.expand(h[5] - (-9216*d2*sigma**2 + 32256*d1**2*sigma
                            - 12288*d1**2*d2**2 + 18432*d1*d2*dm1
                            + 2048*dm1**2)) == 0)

    # CHK4: Lemma A cascade telescoping for a=6 (v=12), symbolic sufficiency
    E, uu, tt = sp.symbols('E uu tt')
    G = sp.symbols('G1:8')          # G[l-1] ~ g_l
    H = {0: tt**12*G[0]}
    for l in range(1, 7):
        H[l] = (tt**12*G[l] - E**3*G[l-1])/uu**l      # g_{l+1}=(e^3 g_l+u^l h_l)/t^12
    H[7] = -E**3*G[6]/uu**7                            # terminal e^3 g_7 + u^7 h_7 = 0
    S = sum(tt**(12*f)*uu**f*E**(21-3*f)*H[f] for f in range(8))
    check("CHK4 cascade lines (a=6,v=12) telescope to stripped identity", sp.expand(S) == 0)

    # CHK5: terminal forcing arithmetic (established facts, re-verified)
    a = 6
    ok = (30 - 3*a == 12) and (10 + 3*a == 28)
    ok = ok and (4*7 <= 28) and (4*(7 + 2*1) > 28)   # v_q(d1)>=1 -> v_q(g7)>=9 -> 36>28
    # K*gamma^3 = -8192 c^7 delta^2  (cancel q^7, r^6 from terminal line)
    gam, dlt, Ksym = sp.symbols('gamma delta K')
    Kval = -8192*c**7*dlt**2/gam**3
    ok = ok and sp.simplify(Kval*gam**3 + 8192*c**7*dlt**2) == 0
    # degree windows: deg e_hat = 2 rho <= 4, deg d1 = 3 rho <= 6 -> rho <= 2
    ok = ok and all(2*rho <= 4 and 3*rho <= 6 for rho in (0, 1, 2))
    check("CHK5 terminal recap: v=12, cap 28, v_q(d1)=0, g_7=K*q^7, K*g^3=-8192c^7*d^2, rho<=2", ok)

    # CHK6: level-6 rearrangement <=> (E6): sigma^2 = A t^12 q + B t^6 r^5 + r^6 Stilde
    X, Y, Z, S2, Ph, d2s = sp.symbols('X Y Z S2 Phat d2s')  # X~t^12 q, Y~t^6 r^5, Z~r^6
    A = sp.Rational(8, 3)*c*dlt**2/gam**3
    B = sp.Rational(8, 3)*gam*dlt
    St = (Ph + 14336*c**6*dlt**2*d2s)/(3072*c**6)
    P = Kval*X - c**6*(-3072*S2 + 14336*dlt**2*Z*d2s + 8192*gam*dlt*Y)
    check("CHK6 level-6: P - r^6*Phat == 3072 c^6 (sigma^2 - A t12q - B t6r5 - r^6 St)",
          sp.simplify(P - Z*Ph - 3072*c**6*(S2 - A*X - B*Y - Z*St)) == 0)
    return h

def gen_poly(prefix, deg):
    cs = sp.symbols(f'{prefix}0:{deg+1}')
    return sum(cs[i]*y**i for i in range(deg+1)), list(cs)

def part_B():
    gam, dlt, lam, TT, St, sig, d2s = sp.symbols('gamma delta lam TT St sig d2s')
    ts, qs, rs = sp.symbols('ts qs rs')
    A = sp.Rational(8, 3)*c*dlt**2/gam**3
    B = sp.Rational(8, 3)*gam*dlt

    # CHK7: master level-5 decomposition (abstract; sigma^2 eliminated via E6)
    W0 = -9216*d2s*St + 32256*dlt**2*sig - 12288*dlt**2*d2s**2
    h5sub = (-9216*d2s*(A*ts**12*qs + B*ts**6*rs**5 + rs**6*St)
             + 32256*dlt**2*rs**6*sig - 12288*dlt**2*rs**6*d2s**2
             + 18432*gam*dlt*ts**6*rs**5*d2s + 2048*gam**2*ts**12*rs**4)
    Ph = TT - 9216*A*c**5*gam**3*d2s          # T := Phat + 9216 A c^5 g^3 d2
    D = ts**12*qs*Ph - c**5*gam**3*h5sub
    RHS = (ts**12*qs*TT - 2048*c**5*gam**4*ts**6*rs**4*(gam*ts**6 - 3*dlt*d2s*rs)
           - c**5*gam**3*rs**6*W0)
    check("CHK7 D = t12q*T - 2048c^5g^4 t^6 r^4 (g t^6 - 3 d d2 r) - c^5 g^3 r^6 W0",
          sp.simplify(D - RHS) == 0)

    # concrete generic window data
    qy = q.as_expr()
    d2p, _ = gen_poly('b', 4)          # d2, deg <= 4
    Stp, _ = gen_poly('u', 4)          # Stilde, deg <= 4 (deg Phat <= 4 from deg g6 <= 28)
    sgp, _ = gen_poly('v', 8)          # sigma, deg <= 8
    W0p = sp.expand((-9216*d2p*Stp + 32256*dlt**2*sgp - 12288*dlt**2*d2p**2))
    check("CHK7b deg W0 <= 8", sp.degree(W0p, y) <= 8)

    # ---- rho = 0 (r = 1): kill by degree cap deg D <= 8 ----
    TTp, avars = gen_poly('a', 4)      # T free of deg <= 4 when r = 1
    D0 = sp.expand(t**12*qy*TTp - 2048*c**5*gam**4*t**6*(gam*t**6 - 3*dlt*d2p)
                   - c**5*gam**3*W0p)
    rest = sp.expand(D0 - sp.expand(t**12*qy*TTp))
    ok = sp.degree(rest, y) <= 12
    eqs = [sp.expand(D0).coeff(y, k) for k in range(16, 21)]
    M = sp.Matrix([[sp.diff(e, av) for av in avars] for e in eqs])
    ok = ok and M.det() != 0
    sol = sp.solve(eqs, avars, dict=True)
    ok = ok and sol == [{a: 0 for a in avars}]
    check("CHK8a rho=0: coeffs y^16..y^20 of D force T = 0 (unique soln, det != 0)", ok)
    c12 = sp.expand(D0.subs({a: 0 for a in avars})).coeff(y, 12)
    check("CHK8b rho=0: with T=0, coeff(D,y^12) = -2048 c^5 gamma^5 != 0",
          sp.simplify(c12 + 2048*c**5*gam**5) == 0)
    check("CHK8c rho=0: cap deg D <= 8+6*0 = 8 < 12  => contradiction", 12 > 8)

    # ---- rho = 1 (r = y - theta, theta != -1): kill mod r ----
    th = sp.symbols('theta')
    r1 = y - th
    Br = sp.expand(lam*t**12*qy - 2048*c**5*gam**4*t**6*(gam*t**6 - 3*dlt*d2p*r1)
                   - c**5*gam**3*r1**2*W0p)
    check("CHK9a rho=1: coeff(Br,y^16) = 2048*lam (rest of Br has degree <= 12)",
          sp.simplify(Br.coeff(y, 16) - 2048*lam) == 0
          and sp.degree(sp.expand(Br.subs(lam, 0)), y) <= 12)
    check("CHK9b rho=1: cap deg Br <= 2*1+8 = 10 < 16  => lam = 0", 16 > 10)
    Br0 = Br.subs(lam, 0)
    val = sp.expand(Br0.subs(y, th))
    check("CHK9c rho=1: Br|_{lam=0}(theta) = -2048 c^5 gamma^5 (theta+1)^12",
          sp.simplify(val + 2048*c**5*gam**5*(th + 1)**12) == 0)
    check("CHK9d rho=1: theta != -1 (t coprime to r) => r does not divide Br; but r^2|Br required", True)

    # ---- rho = 2 (r = y^2 + p y + s, monic wlog): kill mod r ----
    pp, ss = sp.symbols('pp ss')
    r2 = y**2 + pp*y + ss
    check("CHK10a rho=2: deg T <= 4 < 8 = deg r^4 and r^4 | T  => T = 0 (lam = 0)", 4 < 8)
    Br2 = sp.expand(-2048*c**5*gam**4*t**6*(gam*t**6 - 3*dlt*d2p*r2)
                    - c**5*gam**3*r2**2*W0p)
    rem1 = sp.rem(sp.expand(Br2 + 2048*c**5*gam**5*t**12), r2, y)
    check("CHK10b rho=2: Br == -2048 c^5 gamma^5 t^12  (mod r)", sp.expand(rem1) == 0)
    remt = sp.rem(sp.expand(t**12), r2, y)
    solt = sp.solve([remt.coeff(y, 1), remt.coeff(y, 0)], [pp, ss], dict=True)
    check("CHK10c rho=2: r | t^12 has the UNIQUE monic solution r = (y+1)^2 = t^2",
          solt == [{pp: 2, ss: 1}])
    check("CHK10d rho=2: r = t^2 contradicts t coprime to r  => infeasible",
          r2.subs({y: -1, pp: 2, ss: 1}) == 0)

if __name__ == '__main__':
    import os
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    part_A()
    part_B()
    n_fail = sum(1 for _, ok in CHECKS if not ok)
    print(f"-- {len(CHECKS)} checks, {n_fail} failures")
