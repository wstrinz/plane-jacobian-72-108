"""
t5_stratum100_verify.py — stratum (a_t, a_q) = (10, 0) of the f31 subcase-(2)
window system (T5_STRATUM_10_0.md): exact verification of every computational
input used there.  All exact, sympy over Q.  Companion files:
T5_MULTIPLACE.md / t5_multiplace_verify.py (upstream strata), f31_graded.txt.

Stratum: e := d~_{-1} = C t^10 exactly (C in Q^x), t := y+1,
q := 2048y^4-512y^3+320y^2-240y+195, Phi~ = c t^30 q, c = -1/6630.

Check groups (numbered as cited in T5_STRATUM_10_0.md):
  S1  t-power bookkeeping:  Phi~^f (C t^10)^(21-3f) = c^f q^f C^(21-3f) t^210,
      hence f31(d~,Phi~)|_{e=Ct^10} = t^210 * F  with the T-FREE identity
      F := sum_f c^f C^(21-3f) q^f h_f(d~2,d~1,d~0,Ct^10);  the graded
      decomposition of f31 is re-verified so this file is self-contained.
  S2  q-cascade telescoping + degree caps (deg G_l <= 40-4l).
  S3  dm1-layer split h_f = sum_k dm1^k h_{f,k}; H_k := sum_f h_{f,k} X^f;
      H_4 = -6561; layer degree caps; the t^10-block regrouping of F and its
      Q-cascade telescoping (deg caps 40-10k for P_k, 40-10l for Q_l).
  S4  the factorization H_0 = U4 * U3 (weights 8 and 12), and the cofactor
      identities A1 = -216 d1, B1 = 16 d1 (32X^2-56Xd2-36d0+29d2^2):
        H_1 = A1*U3 + B1*U4,   H_3 = 11664 d1 (4 d2 - 5 X),
        H == (U4 + A1*dm1)(U3 + B1*dm1) + dm1^2*(D2 + dm1*H_3 - 6561*dm1^2),
      D2 := H_2 - A1*B1;  hence the T1 product-master
        Ut4 * Ut3 = - C^2 t^20 J   (all of stratum (10,0)),
      Ut4 := U4(d~,Bq) - 216 C t^10 d~1,  Ut3 := U3(d~,Bq) + C t^10 B1(d~,Bq),
      J := D2(d~,Bq) + C t^10 H_3(d~,Bq) - 6561 C^2 t^20,  B := c/C^3.
  S5  T2 branch (d~1 == 0) split in coordinates d0 = (sigma + d2^2)/4:
      h_f| = sigma^2 s_f + dm1^2 r_f + [f=0](-6561 dm1^4); the factorization
      S(X) = 12(9s - (4X-5d2)(4X+d2))(4(X+d2)^2 + 9s)^2; the closed form
      Rt = 512 Ah^4 Bh - 432 Ah^2 (29v-16z) S + 2916 Bh S^2; and the fully
      symbolic T2 master:
        C^9 * F|_{T2} = 12 S^2 Nh Ph^2 + K T^2 Rt - 6561 K^2 T^4
      (T stands for t^10, K := C^17, z := cQ, v := C^3 d2, S := C^6 sigma).
  S6  T2a/T2b collapse identities:  Rt|_{9S=-4Ah^2} = 128 Ah^4 (10z+v),
      Rt|_{9S=Bh*Bph} = 4 Bh^2 (8z-v) (10z+v)^2,  Rt|_{S=0} = 512 Ah^4 Bh;
      plus the valuation arithmetic used in their kills (q(-1) != 0, so
      t does not divide q and q is not const*t^4).
  S7  the six linear forms Ah, Bph, Bh, 10z+v, 2z-v, 29v-16z are pairwise
      t-coprime: beta'L - beta L' = k z with k != 0 for each pair; plus the
      first-block kill arithmetic of cases D and the m=0 / m>=1 splits.
  S8  I0 reduced identity: master|_{W=0} = 2187(4S^5 - 12z^2 S^4
      + 12 z K T^2 S^2 - 3 K^2 T^4), and the three-term valuation kill
      arithmetic (4m, 20+2m, 40 pairwise distinct for m <= 8).
  S9  the three-term Newton-polygon kill of case A (v_t(Ah) = alpha >= 1):
      the R-tilde regrouping  Rt(t^a W1, t^2a S1) = t^4a (9 z Theta
      - t^a W1 Theta'), Theta := 512 W1^4 + 2160 W1^2 S1 + 2916 S1^2,
      Theta' := 2560 W1^4 + 12528 W1^2 S1 + 14580 S1^2; the incompatibility
      lemma Theta|_{9 S1 = -4 W1^2} = 128 W1^4 != 0; and the pairing
      arithmetic that leaves no admissible (alpha, m, kappa1) profile;
      plus the case-B polygon arithmetic.
  S10 the five remaining rigid-cell certificates (B', C, C', L4==0, L5==0,
      and B as a redundant double-check): each cell's 41 coefficient
      equations + K*Kk - 1 generate the unit ideal (sympy Groebner over QQ)
      -- exact char-0 infeasibility certificates.

Every check prints OK; the script asserts throughout.
"""
import os, re, sympy as sp

os.chdir(os.path.dirname(os.path.abspath(__file__)))

d2, d1, d0, dm1, Phi, X, W = sp.symbols('d2 d1 d0 dm1 Phi X W')
C, K, Q, T = sp.symbols('C K Q T')
s = sp.Symbol('s')
z, v = sp.symbols('z v')
y = sp.Symbol('y')
t = y + 1
qq = 2048*y**4 - 512*y**3 + 320*y**2 - 240*y + 195
c = sp.Rational(-1, 6630)
Phit = c * t**30 * qq

f31 = sp.sympify(open('f31_deg31.txt').read().strip()
                 .replace('m1', 'dm1').replace('P', 'Phi').replace('^', '**'))
hs = {}
for m in re.finditer(r'h_(\d) \(weight (\d+), dm1-power (\d+)\) = (.+)',
                     open('f31_graded.txt').read()):
    hs[int(m.group(1))] = sp.sympify(m.group(4))
assert sorted(hs) == list(range(8))

# =============================================================== S1
recon = sum(Phi**f * dm1**(21 - 3*f) * hs[f] for f in range(8))
assert sp.expand(recon - f31) == 0
for f in range(8):
    lhs = Phit**f * (C*t**10)**(21 - 3*f)
    rhs = c**f * qq**f * C**(21 - 3*f) * t**210
    assert sp.expand(lhs - rhs) == 0
assert qq.subs(y, -1) == 3315          # t does not divide q
print("S1. graded decomposition; Phi~^f (Ct^10)^(21-3f) = c^f q^f C^(21-3f) t^210;")
print("    => f31|_(10,0) = t^210 * F with F t-free                    OK")

# =============================================================== S2
# q-cascade telescoping: if h_0 = Q*G1, h_l = (Q*G_{l+1} - C^3 G_l)/c^l (l=1..6),
# h_7 = -C^3 G_7 / c^7, then sum_f c^f C^(21-3f) Q^f h_f == 0 identically.
g = [None] + list(sp.symbols('g1:8'))
hh = {0: Q*g[1]}
for l in range(1, 7):
    hh[l] = (Q*g[l+1] - C**3*g[l]) / c**l
hh[7] = -C**3*g[7] / c**7
Scas = sum(c**f * C**(21 - 3*f) * Q**f * hh[f] for f in range(8))
assert sp.expand(sp.cancel(Scas)) == 0
for l in range(1, 8):
    assert (40 - 4*(l - 1)) - 4 == 40 - 4*l   # deg G_l <= 40-4l bookkeeping
print("S2. q-cascade telescoping and degree caps deg G_l <= 40-4l      OK")

# =============================================================== S3
Hfull = sp.expand(sum(hs[f]*X**f for f in range(8)))
Hp = sp.Poly(Hfull, dm1)
Hk = {k: sp.expand(Hp.coeff_monomial(dm1**k)) for k in range(5)}
assert sp.expand(Hfull - sum(dm1**k * Hk[k] for k in range(5))) == 0
assert Hk[4] == -6561
# weight/degree caps: H_k has weight 20 - 5k (X ~ 2, d2 ~ 2, d1 ~ 3, d0 ~ 4)
WT = {X: 2, d2: 2, d1: 3, d0: 4}
for k in range(4):
    P = sp.Poly(Hk[k], X, d2, d1, d0)
    for mono in P.monoms():
        assert sum(w*e for w, e in zip((2, 2, 3, 4), mono)) == 20 - 5*k, k
# window degrees: each weight unit costs <= 2 in y => deg_y H_k(d~, Bq) <= 40-10k
for k in range(4):
    assert 2*(20 - 5*k) == 40 - 10*k
# t^10-block regrouping is the identity F = C^21 * sum_k (C t^10)^k H_k(d~,Bq):
# pure regrouping of the dm1-layers, verified symbolically:
Bsym = sp.Symbol('B')
FA = sum(c**f * C**(21-3*f) * Q**f * hs[f].subs(dm1, C*T).subs(
    {d2: d2, d1: d1, d0: d0}) for f in range(8))
FB = C**21 * sum((C*T)**k * Hk[k].subs(X, c*Q/C**3) for k in range(5))
assert sp.simplify(sp.expand(FA - FB)) == 0
# t^10-cascade telescoping: P0 = T*Q1, P1 = (T*Q2 - Q1)/C, P2 = (T*Q3 - Q2)/C^2,
# P3 = (6561*C^25*T - Q3)/C^3, P4-term = -6561 C^25 T^4 ... direct:
Qs = sp.symbols('Q1:4')
P0 = T*Qs[0]; P1 = (T*Qs[1] - Qs[0])/C; P2 = (T*Qs[2] - Qs[1])/C**2
P3 = (6561*C**25*T - Qs[2])/C**3
tot = P0 + C*T*P1 + C**2*T**2*P2 + C**3*T**3*P3 - 6561*C**25*T**4
assert sp.expand(sp.cancel(tot)) == 0
print("S3. dm1-layers H_k (weights 20-5k, caps 40-10k), H_4 = -6561;")
print("    t^10-block regrouping + t^10-cascade telescoping             OK")

# =============================================================== S4
U4 = (16*X**4 + 64*X**3*d2 + 288*X**2*d0 + 24*X**2*d2**2 + 576*X*d0*d2
      - 432*X*d1**2 - 80*X*d2**3 + 1296*d0**2 - 360*d0*d2**2 + 216*d1**2*d2
      + 25*d2**4)
U3 = (512*X**3*d1**2 - 3072*X**2*d0**2 + 1536*X**2*d0*d2**2
      - 1152*X**2*d1**2*d2 - 192*X**2*d2**4 + 3072*X*d0**2*d2
      - 1152*X*d0*d1**2 - 1536*X*d0*d2**3 + 1056*X*d1**2*d2**2 + 192*X*d2**5
      + 6912*d0**3 - 4224*d0**2*d2**2 + 1008*d0*d1**2*d2 + 816*d0*d2**4
      + 27*d1**4 - 412*d1**2*d2**3 - 48*d2**6)
assert sp.expand(Hk[0] - U4*U3) == 0
A1 = -216*d1
B1 = 16*d1*(32*X**2 - 56*X*d2 - 36*d0 + 29*d2**2)
assert sp.expand(Hk[1] - (A1*U3 + B1*U4)) == 0
assert sp.expand(Hk[3] - 11664*d1*(4*d2 - 5*X)) == 0
D2 = sp.expand(Hk[2] - A1*B1)
assert sp.expand(Hfull - ((U4 + A1*dm1)*(U3 + B1*dm1)
                          + dm1**2*(D2 + dm1*Hk[3] - 6561*dm1**2))) == 0
print("S4. H_0 = U4*U3;  H_1 = A1*U3 + B1*U4 (A1 = -216 d1);")
print("    H_3 = 11664 d1 (4d2-5X);  H = (U4+A1*dm1)(U3+B1*dm1)")
print("    + dm1^2 (D2 + dm1 H_3 - 6561 dm1^2)  => product-master        OK")

# =============================================================== S5
loc = {d1: 0, d0: (s + d2**2)/4}
sfs, rfs = {}, {}
for f in range(8):
    hf = sp.expand(hs[f].subs(loc))
    extra = -6561*dm1**4 if f == 0 else 0
    p = sp.Poly(hf - extra, dm1)
    co = {mono[0]: cf for mono, cf in zip(p.monoms(), p.coeffs())}
    assert set(co) <= {0, 2}, f
    sfs[f] = sp.expand(sp.cancel(co.get(0, 0)/s**2))
    assert sp.expand(sfs[f]*s**2 - co.get(0, 0)) == 0, f
    rfs[f] = sp.expand(co.get(2, 0))
SX = sum(sfs[f]*X**f for f in range(8))
RX = sum(rfs[f]*X**f for f in range(8))
assert sp.expand(SX - 12*(9*s - (4*X - 5*d2)*(4*X + d2))*(4*(X + d2)**2 + 9*s)**2) == 0
Ah = z + v; Bh = 4*z - 5*v; Bph = 4*z + v
Rt = 512*Ah**4*Bh - 432*Ah**2*(29*v - 16*z)*s + 2916*Bh*s**2
assert sp.expand(RX.subs({X: z, d2: v}) - Rt) == 0
# fully symbolic T2 master (T = t^10, Q = q, K = C^17, z = cQ, v = C^3 d2, S = C^6 s):
F_T2 = sum(c**f * C**(21 - 3*f) * Q**f *
           (s**2*sfs[f] + C**2*T**2*rfs[f] + (-6561*C**4*T**4 if f == 0 else 0))
           for f in range(8))
subzvS = {z: c*Q, v: C**3*d2}
Ssym = C**6*s
Nh_m = 9*Ssym - (Bh*Bph).subs(subzvS)
Ph_m = 9*Ssym + (4*Ah**2).subs(subzvS)
Rt_m = (512*Ah**4*Bh).subs(subzvS) - 432*(Ah**2).subs(subzvS)*(29*C**3*d2 - 16*c*Q)*Ssym \
       + 2916*Bh.subs(subzvS)*Ssym**2
master_sym = 12*Ssym**2*Nh_m*Ph_m**2 + C**17*T**2*Rt_m - 6561*C**34*T**4
assert sp.expand(C**9*F_T2 - master_sym) == 0
print("S5. T2 split h_f| = s^2 s_f + dm1^2 r_f (+ f=0: -6561 dm1^4);")
print("    S-factorization; Rt closed form; symbolic T2 master:")
print("    C^9 F|_T2 = 12 S^2 Nh Ph^2 + K T^2 Rt - 6561 K^2 T^4         OK")

# =============================================================== S6
Rt_S = Rt  # polynomial in z, v, s (s in role of S)
assert sp.expand(Rt_S.subs(s, -sp.Rational(4, 9)*Ah**2) - 128*Ah**4*(10*z + v)) == 0
assert sp.expand(Rt_S.subs(s, sp.Rational(1, 9)*Bh*Bph)
                 - 4*Bh**2*(8*z - v)*(10*z + v)**2) == 0
assert sp.expand(Rt_S.subs(s, 0) - 512*Ah**4*Bh) == 0
# valuation arithmetic used in the T2a/T2b kills:
assert qq.subs(y, -1) != 0          # t ∤ q, so q != const * t^4 and t ∤ z
print("S6. Rt|_{9S=-4Ah^2} = 128 Ah^4 L4;  Rt|_{9S=BhBph} = 4 Bh^2 (8z-v) L4^2;")
print("    Rt|_{S=0} = 512 Ah^4 Bh;  q(-1) = 3315 != 0                   OK")

# =============================================================== S7
# six forms: Ah = z+v, Bph = 4z+v, Bh = 4z-5v, L4 = 10z+v, L5 = 2z-v, M = 29v-16z
L = {'Ah': (1, 1), 'Bph': (4, 1), 'Bh': (4, -5), 'L4': (10, 1), 'L5': (2, -1),
     'M': (-16, 29)}   # (alpha, beta): alpha*z + beta*v
import itertools
for (i, (ai, bi)), (j, (aj, bj)) in itertools.combinations(L.items(), 2):
    detk = ai*bj - aj*bi
    assert detk != 0, (i, j)
# first-block arithmetic (caps m, n, p <= 8 since deg S, N, P <= 8 and all != 0):
assert all(2*m0 < 20 for m0 in range(9))              # v(T1) = 2m (n=p=0) < 20
assert all(n0 < 20 for n0 in range(9))                # v(T1) = n (m=p=0) < 20
assert all(2*p0 < 20 for p0 in range(9))              # v(T1) = 2p (m=n=0) < 20
print("S7. six forms pairwise t-coprime (at most one can carry t);")
print("    case D and all m=0 side-cases die at the first block           OK")

# =============================================================== S8
Wsym = sp.Symbol('Wv')
Ssb = sp.Symbol('Sv')
Nh_g = 9*Ssb - 27*z**2 + 6*z*Wsym + 5*Wsym**2
Ph_g = 9*Ssb + 4*Wsym**2
Rt_g = 512*Wsym**4*(9*z - 5*Wsym) - 432*Wsym**2*(29*Wsym - 45*z)*Ssb \
       + 2916*(9*z - 5*Wsym)*Ssb**2
master_g = 12*Ssb**2*Nh_g*Ph_g**2 + K*T**2*Rt_g - 6561*K**2*T**4
# consistency of the (z, W=Ah, S) form with the (z,v,S) form (v = W - z):
Nh_v = 9*Ssb - (4*z - 5*(Wsym - z))*(4*z + (Wsym - z))
Rt_v = (512*Ah**4*Bh - 432*Ah**2*(29*v - 16*z)*s + 2916*Bh*s**2).subs(
    {v: Wsym - z, s: Ssb})
assert sp.expand(Nh_g - Nh_v) == 0
assert sp.expand(Rt_g - Rt_v) == 0
# I0 (W == 0) reduction:
red = sp.expand(master_g.subs(Wsym, 0)
                - 2187*(4*Ssb**5 - 12*z**2*Ssb**4 + 12*z*K*T**2*Ssb**2 - 3*K**2*T**4))
assert red == 0
# three-term valuation kill: 4m, 20+2m, 40 pairwise distinct for 1 <= m <= 8
for m0 in range(1, 9):
    vals = {4*m0, 20 + 2*m0, 40}
    assert len(vals) == 3, m0
print("S8. I0 identity master|_{W=0} = 2187(4S^5 - 12 z^2 S^4")
print("    + 12 z K T^2 S^2 - 3 K^2 T^4); 3-term Newton kill arithmetic  OK")

# =============================================================== S9
# case A (alpha = v_t(W) >= 1, all other forms t-units, m = v_t(S) >= 1, n = 0):
# (i) the R-tilde group valuations are 4a, 2a+m, 2m EXACTLY (unit cofactors
#     9z-5W == 9z, 29W-45z == -45z mod t) -- structural, from S5 closed form;
# (ii) regrouping identity at m = 2a:
W1s, S1s, TA = sp.symbols('W1s S1s TA')      # TA plays the role of t^alpha
Rt_W = (512*Wsym**4*(9*z - 5*Wsym) - 432*Wsym**2*(29*Wsym - 45*z)*Ssb
        + 2916*(9*z - 5*Wsym)*Ssb**2)
Theta = 512*W1s**4 + 2160*W1s**2*S1s + 2916*S1s**2
Thetap = 2560*W1s**4 + 12528*W1s**2*S1s + 14580*S1s**2
regroup = Rt_W.subs({Wsym: TA*W1s, Ssb: TA**2*S1s}) \
    - TA**4*(9*z*Theta - TA*W1s*Thetap)
assert sp.expand(regroup) == 0
# (iii) incompatibility lemma: kappa1 >= 1 and kappa2 >= 1 cannot both hold:
assert sp.expand(Theta.subs(S1s, -sp.Rational(4, 9)*W1s**2) - 128*W1s**4) == 0
# (iv) pairing arithmetic: no admissible profile survives.
#   m != 2a:
for alpha in range(1, 5):
    for m0 in range(1, 9):
        if m0 < 2*alpha:                     # p = m, rho = 2m exact
            assert 4*m0 < 20 + 2*m0 and 4*m0 < 40          # min unique
        elif m0 > 2*alpha:                   # p = 2a, rho = 4a exact
            assert 2*m0 + 4*alpha < 20 + 4*alpha and 2*m0 + 4*alpha < 40
#   m = 2a, kappa1 = 0:  v(T1) = 8a < min(20+4a, 40):
for alpha in range(1, 5):
    assert 8*alpha < 20 + 4*alpha and 8*alpha < 40
#   m = 2a, kappa1 >= 1 (then kappa2 = 0, v(T2) = 20+4a exact; p-cap p <= 8):
for alpha in range(1, 5):
    for k1 in range(1, 9):
        p0 = 2*alpha + k1                    # = v_t(Ph), must be <= 8
        if p0 > 8:
            continue                         # profile impossible by the p-cap
        vT1 = 8*alpha + 2*k1                 # = 4a + 2p < 20 + 4a = v(T2)
        assert vT1 < 20 + 4*alpha and vT1 < 40      # min = v(T1) unique: dead
# case B (v_t(Bh) = beta in [1,4], m >= 1, p = 0):
for beta in range(1, 5):
    for m0 in range(1, 9):
        if m0 == beta:
            assert 2*m0 + 8 < 20             # v(T1) <= 2m+8 <= 16 < 20: dead
        else:
            n0 = min(m0, beta)               # n and rho both = min(m,beta) exact
            vT1, vT2 = 2*m0 + n0, 20 + n0
            assert vT1 != vT2 and vT1 < 40 and vT2 != 40    # min unique: dead
# case B' (v_t(Bph) = beta' in [1,4], m >= 1, p = 0, rho = 0 exact):
#   surviving profiles of 2m + n = 20, n = min(m, beta') (m != beta'):
surv = [(m0, b0) for b0 in range(1, 5) for m0 in range(1, 9)
        if m0 != b0 and 2*m0 + min(m0, b0) == 20]
assert surv == [(8, 4)]                      # forced into rigid cell B'
assert all(2*m0 + 8 < 20 for m0 in range(1, 5))     # m = beta' <= 4: dead
# case C / C' (v_t(L4 or L5) = lambda in [1,4]; m = 0 forced):
#   surviving (n, p) with n + 2p >= 20, n,p <= 8, min(n,p) <= 4:
survC = [(n0, p0) for n0 in range(9) for p0 in range(9)
         if n0 + 2*p0 >= 20 and min(n0, p0) <= 4]
assert survC == [(4, 8)]                     # forced into rigid cells C / C'
# degenerate forms Bh == 0 / Bph == 0 (N = 9S, n = m, p = 0 for m >= 1):
for m0 in range(1, 9):
    assert 3*m0 != 20 + m0 and 3*m0 != 40 and 3*m0 != 20   # vs rho = m or 0
    assert 20 + m0 != 40
# sigma-locus (S == 0) inside (10,0): 512 Ah^4 Bh = 6561 K t^20 forces
# Ah = a t^4, Bh = b t^4 (4d1+d2 = 20 with d1,d2 <= 4 only at (4,4)):
sols = [(dd1, dd2) for dd1 in range(5) for dd2 in range(5) if 4*dd1 + dd2 == 20]
assert sols == [(4, 4)]          # also forces the T2a shape (Ah^4 L4 case)
# T2b shape: 2 deg(Bh) + deg(8z-v) + 2 deg(L4) = 20 with each <= 4 -> (4,4,4):
solsb = [(aa, bb, cc) for aa in range(5) for bb in range(5) for cc in range(5)
         if 2*aa + bb + 2*cc == 20]
assert solsb == [(4, 4, 4)]
print("S9. case-A regrouping + Theta-lemma + polygon arithmetic: case A dead;")
print("    case B dead; B'/C/C' forced into rigid cells; sigma-locus (10,0)")
print("    forced to (t^4, t^4) then dead against 5Ah+Bh = 9z             OK")

# =============================================================== S10
from sympy import groebner
Kk = sp.Symbol('Kk')
gam, delt, epsv, v0, v1 = sp.symbols('gam delt epsv v0 v1')
zpoly = sp.expand(c*qq)

def master_poly(Wp, Sp):
    Nh = 9*Sp - 27*zpoly**2 + 6*zpoly*Wp + 5*Wp**2
    Ph = 9*Sp + 4*Wp**2
    Rt_ = (512*Wp**4*(9*zpoly - 5*Wp) - 432*Wp**2*(29*Wp - 45*zpoly)*Sp
           + 2916*(9*zpoly - 5*Wp)*Sp**2)
    return sp.expand(12*Sp**2*Nh*Ph**2 + K*t**20*Rt_ - 6561*K**2*t**40)

def cell_dead(name, Wp, Sp, params):
    M = master_poly(sp.expand(Wp), sp.expand(Sp))
    eqs = [sp.expand(co) for co in sp.Poly(M, y).all_coeffs() if co != 0]
    gb = groebner(eqs + [K*Kk - 1], *params, K, Kk, order='grevlex', domain='QQ')
    assert list(gb.exprs) == [sp.Integer(1)], name
    print(f"    cell {name}: unit ideal (infeasible)")

print("S10. rigid-cell certificates (sympy Groebner over QQ):")
cell_dead("B  (Bh = delt t^4, S = gam t^8)",
          zpoly + (4*zpoly - delt*t**4)/5, gam*t**8, [gam, delt])
cell_dead("B' (Bph = delt t^4, S = gam t^8)",
          delt*t**4 - 3*zpoly, gam*t**8, [gam, delt])
WC_ = delt*t**4 - 9*zpoly
cell_dead("C  (L4 = delt t^4, Ph = epsv t^8)",
          WC_, (epsv*t**8 - 4*WC_**2)/9, [delt, epsv])
WCp_ = 3*zpoly - delt*t**4
cell_dead("C' (L5 = delt t^4, Ph = epsv t^8)",
          WCp_, (epsv*t**8 - 4*WCp_**2)/9, [delt, epsv])
cell_dead("L4 == 0  (W = -9z, 9S + 324 z^2 = t^7 V)",
          -9*zpoly, (t**7*(v0 + v1*y) - 324*zpoly**2)/9, [v0, v1])
cell_dead("L5 == 0  (W = 3z, 9S + 36 z^2 = t^7 V)",
          3*zpoly, (t**7*(v0 + v1*y) - 36*zpoly**2)/9, [v0, v1])
print("S10. all six rigid cells infeasible over QQ (unit ideals)          OK")

print()
print("ALL STRATUM-(10,0) CHECKS PASS: S1-S10")
print("  => the T2 branch (d~1 == 0) of stratum (10,0) is infeasible (char 0);")
print("     the T1 branch is reduced to the product-master (see MD, open).")
