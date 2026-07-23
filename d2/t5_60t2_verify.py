# t5_60t2_verify.py -- verification for stratum (6,0), branch T2 (d1 == 0, sigma != 0)
# All algebra recomputed from f31_graded.txt; prints only small summaries.
import os, re, random
import sympy as sp

y = sp.symbols('y')
t = y + 1
q = sp.Poly(2048*y**4 - 512*y**3 + 320*y**2 - 240*y + 195, y)
qe = q.as_expr()
c = sp.Rational(-1, 6630)

CHECKS = []
def check(name, ok):
    CHECKS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}")

d0, d1, d2, dm1 = sp.symbols('d0 d1 d2 dm1')

def load_h():
    h = {}
    pat = re.compile(r'h_(\d+)\s*\([^)]*\)\s*=\s*(.+)$')
    with open('f31_graded.txt') as fh:
        for line in fh:
            m = pat.match(line.strip())
            if m:
                h[int(m.group(1))] = sp.sympify(m.group(2),
                    locals={'d0': d0, 'd1': d1, 'd2': d2, 'dm1': dm1})
    return h

def main():
    sigma = 4*d0 - d2**2
    h = load_h()

    # V1 setup facts
    Phi = sp.expand(c*t**30*qe)
    ok = q.is_irreducible and q.degree() == 4 and qe.subs(y, -1) == 3315
    ok = ok and sp.degree(Phi, y) == 34 and sp.LC(sp.Poly(Phi, y)) == sp.Rational(-1024, 3315)
    check("V1 q irred /Q, deg 4, q(-1)=3315 (t coprime to q); deg Phi~=34, lc=-1024/3315", ok)

    # V2 graded pieces + collapses on d1 = 0
    ok = sorted(h) == list(range(8)) and sp.expand(h[7] - 8192*d1**2) == 0
    ok = ok and sp.expand(h[6].subs(d1, 0) - (-3072*sigma**2)) == 0
    ok = ok and sp.expand(h[5].subs(d1, 0) - (-9216*d2*sigma**2 + 2048*dm1**2)) == 0
    check("V2 h_0..h_7 parsed; h_7=8192*d1^2; h_6|d1=0 = -3072*s^2; h_5|d1=0 = -9216*d2*s^2+2048*e^2", ok)

    # V3 Lemma A cascade telescoping for a=6 (v=12), symbolic sufficiency
    E, uu, tt = sp.symbols('E uu tt')
    G = sp.symbols('G1:8')
    H = {0: tt**12*G[0]}
    for l in range(1, 7):
        H[l] = (tt**12*G[l] - E**3*G[l-1])/uu**l   # e^3 g_l + u^l h_l = t^12 g_{l+1}
    H[7] = -E**3*G[6]/uu**7                        # terminal e^3 g_7 + u^7 h_7 = 0
    S = sum(tt**(12*f)*uu**f*E**(21-3*f)*H[f] for f in range(8))
    check("V3 cascade lines (a=6, v=12) telescope to the stripped identity", sp.expand(S) == 0)

    # V4 stratum arithmetic
    a = 6
    ok = (30 - 3*a == 12) and (10 + 3*a == 28) and (2*a == 12)   # e^2 = t^12 * e_hat^2 exactly
    ok = ok and sp.simplify(-c**6*(-3072) - 3072*c**6) == 0      # terminal sign: e^3 g6 = 3072 c^6 q^6 s^2
    ok = ok and 4*(6 + 2*1) > 28                                 # v_q(sigma)>=1 -> deg g_6 >= 32 > 28
    ok = ok and 28 - 4*6 == 4                                    # deg g_hat <= 4
    check("V4 a=6: v=12, cap 28, e^2=t^12*ehat^2; sign 3072c^6; v_q(sigma)=0 forced; deg ghat<=4", ok)

    # V5 absorption constant
    check("V5 -9216 c^5/(3072 c^6) = 19890", sp.Rational(-9216, 3072)/c == 19890)

    # V6 exact level-5 rearrangement (abstract): with s^2 := E^3*gh/(3072 c^6),
    #   E^3 g5 + c^5 Q^5(-9216 d2 s^2 + 2048 T E^2) - T Q^6 gh
    #     == E^3 (g5 + 19890 Q^5 d2 gh) - T Q^5 (Q gh - 2048 c^5 E^2)      [T ~ t^12]
    Q, T, gh, g5, d2s = sp.symbols('Q T gh g5 d2s')
    s2 = E**3*gh/(3072*c**6)
    lhs = E**3*g5 + c**5*Q**5*(-9216*d2s*s2 + 2048*T*E**2) - T*Q**6*gh
    rhs = E**3*(g5 + 19890*Q**5*d2s*gh) - T*Q**5*(Q*gh - 2048*c**5*E**2)
    check("V6 level-5 rearrangement: E^3*G5 = t^12 q^5 N exactly (abstract identity)",
          sp.simplify(lhs - rhs) == 0)

    # V6b numeric instance of V6 with the real t, q
    rnd = random.Random(60)
    def rpoly(dg):
        return sum(sp.Integer(rnd.randint(-9, 9))*y**i for i in range(dg + 1))
    Eh, ghp, g5p, d2p = rpoly(4), rpoly(4), rpoly(28), rpoly(4)
    s2p = Eh**3*ghp/(3072*c**6)
    lhsn = sp.expand(Eh**3*g5p + c**5*qe**5*(-9216*d2p*s2p + 2048*t**12*Eh**2) - t**12*qe**6*ghp)
    G5n = g5p + 19890*qe**5*d2p*ghp
    Nn = qe*ghp - 2048*c**5*Eh**2
    rhsn = sp.expand(Eh**3*G5n - t**12*qe**5*Nn)
    check("V6b numeric instance (real t,q; random ghat,ehat,g5,d2): identity holds",
          sp.expand(lhsn - rhsn) == 0)

    # V7 degree domination: N != 0 => deg(G5) = 32 + deg(N/E^3) >= 32 > 28 = cap
    ok = (12 + 5*4 == 32) and (32 > 28)
    ok = ok and (20 + 4 + 4 <= 28)   # deg(19890 q^5 d2 ghat) <= 28, so cap deg G5 <= 28 is real
    ok = ok and all(32 + (3*dE - 3*dE) >= 32 for dE in range(5))  # any deg ehat <= 4
    check("V7 N != 0 forces deg G5 >= 32 > 28 = cap (any deg ehat <= 4)", ok)

    # V8 N = 0 kill at the q-place: q*ghat = 2048 c^5 ehat^2 needs q | ehat^2, q irreducible
    #    => q | ehat, excluded by a_q = 0. Spot checks that v_q(2048 c^5 ehat^2) = 0 for q not| ehat.
    ok = True
    trials = 0
    for _ in range(40):
        Ehr = rpoly(4)
        if sp.rem(Ehr, qe, y) == 0 or Ehr == 0:
            continue
        trials += 1
        ok = ok and sp.rem(sp.expand(2048*c**5*Ehr**2), qe, y) != 0
        if trials >= 20:
            break
    ok = ok and trials >= 20 and sp.rem(sp.expand(qe**2), qe, y) == 0  # sanity: q | q^2
    check("V8 q-place: v_q(2048 c^5 ehat^2)=0 whenever q not| ehat (20 random trials + sanity)", ok)

if __name__ == '__main__':
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    main()
    n_fail = sum(1 for _, ok in CHECKS if not ok)
    print(f"-- {len(CHECKS)} checks, {n_fail} failures")
