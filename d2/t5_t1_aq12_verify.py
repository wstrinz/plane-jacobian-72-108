"""
t5_t1_aq12_verify.py -- exact sympy verification of every computational input
used in T5_T1_AQ12.md: the a_q >= 1 surviving strata of the f31 branch,
subcase (2): (a, a_q) in {(2,1),(3,1),(4,1),(5,1),(6,1),(0,2),(1,2),(2,2)},
both terminal branches T1 (d~1 != 0) and T2 (d~1 = 0, sigma != 0).

Notation (T5_NP.md / T5_MULTIPLACE.md / T5_STRATA_50_11.md): t := y+1,
q := 2048y^4 - 512y^3 + 320y^2 - 240y + 195 (irreducible/Q, Galois group S4),
c := -1/6630, Phi~ = c t^30 q, u := c q, e := d~_{-1} (deg <= 10, e != 0),
windows deg d~2 <= 4, d~1 <= 6, d~0 <= 8; sigma := 4 d~0 - d~2^2;
stratum (a, a_q) := (v_t(e), v_q(e)); e = t^a q^(a_q) ebar, gcd(ebar, tq) = 1.
Identity under attack:  0 = sum_{f=0}^{7} Phi~^f e^(21-3f) h_f(d~).

The q-adic Newton polygon lemma (T5_MULTIPLACE.md section 6, proven there):
w := Phi~/e^3 is a root of H(d~,W) = sum_f h_f(d~) W^f (nonzero by T5_NP
Lemma 1), so  min_f over the support of (delta_f + (1-3a_q) f)  must be
attained at least twice, where delta_f := v_q(h_f(d~)).

Checks (all exact, sympy over Q):
  X1. setup: q irreducible/Q, Galois group S4 (transitive), q(-1)=3315,
      graded decomposition + weight caps; the collapse identities WITH the
      d1-parts exact:
        h_7 = 8192 d1^2
        h_6 = -3072 sigma^2 + 14336 d2 d1^2 + 8192 d1 dm1
        h_5 = -9216 d2 sigma^2 + 2048 dm1^2
              + 129024 d0 d1^2 - 44544 d1^2 d2^2 + 18432 d1 d2 dm1
      and: every h_f monomial containing d1 contains d1^i dm1^j with
      i + j >= 2 and i >= 1 (so d1 = q lam, q | dm1  ==>  q^2 divides every
      d1-part); every d1-free monomial has dm1-power even, and if the
      dm1-power is 0 the monomial is divisible by sigma after
      d0 = (d2^2+s)/4 (no pure-d2 terms) -- the W3 structure, re-verified,
      plus H_4 := h_4|_{d1=0, d0=(d2^2+s)/4} = -5184 d2^2 s^2 - 12096 s^3
      + 5632 d2 dm1^2  (NO d2^4 s term -- used in X7).
  X2. a_q = 2 strata (0,2), (1,2), (2,2) are DEAD, both branches, by the
      q-NP alone: v_q(w) = -5; window caps give v_q(d1) <= 1, v_q(sigma) <= 2;
      T1: mu_7 = 2 v_q(d1) - 35 <= -33 < -30 <= mu_f (f <= 6): unique min.
      T2: h_7 = 0, h_6 = -3072 sigma^2 != 0,
          mu_6 = 2 v_q(sigma) - 30 <= -26 < -25 <= mu_f (f <= 5): unique min.
      (sigma = 0 with d1 = 0 is the sigma-locus, dead by T5_MULTIPLACE Thm 2.)
  X3. a_q = 1, branch T1: v_q(w) = -2; mu_7 = 2 v_q(d1) - 14; v_q(d1) = 0
      gives unique min -14 < -12 <= mu_f: so v_q(d1) = 1 (cap floor(6/4)=1),
      d1 = q lam, 0 != lam, deg lam <= 2, and the only possible partner for
      mu_7 = -12 is f = 6 with delta_6 = 0; since the d1-part of h_6 is
      divisible by q^2, delta_6 = 0 forces v_q(sigma) = 0 (sigma != 0).
      Terminal (Lemma A): ebar^3 g_7 = -8192 c^7 q^6 lam^2, so
      v_q(g_7) = 6 and 24 = 4*6 <= deg g_7 <= 10+3a  ==>  a >= 5:
      T1 of (2,1), (3,1), (4,1) DEAD.
  X4. a_q = 1, branch T2: v_q(w) = -2, h_7 = 0, h_6 = -3072 sigma^2,
      w_s := v_q(sigma) <= 2:
      w_s = 0: mu_6 = -12 unique min: DEAD.
      w_s = 1: mu_6 = -10; only possible partner f = 5 with delta_5 = 0;
        but h_5|_{d1=0} = -9216 d2 sigma^2 + 2048 e^2 with q^2 | sigma^2 and
        q^2 | e^2 (a_q = 1): delta_5 >= 2: no partner: DEAD.
      w_s = 2: cascade level 6: ebar^3 g_6 = 3072 c^6 q^(3+2*2) sigmahat^2,
        g_6 != 0 (sigma != 0), so 28 = 4*7 <= deg g_6 <= 10+3a ==> a = 6.
      T2 of (2,1), (3,1), (4,1), (5,1) DEAD; (6,1) T2 reduced to
      v_q(sigma) = 2, i.e. sigma = s q^2 with s in Q^x (deg sigma <= 8).
  X5. (5,1) T1, deg ebar = 1: terminal ebar^3 ghat7 = -8192 c^7 lam^2 with
      deg ghat7 <= 1 forces (v_ebar parity) ghat7 = beta ebar, lam = alpha
      ebar^2 (alpha, beta != 0); the d1-part of h_6 is then divisible by
      ebar^3 (symbolic identity), and the level-6 line divided by q^6 reads
        3072 c^6 sigma^2 = ebar^3 (m + c^6 q^2 Jt) - beta t^15 ebar ,
      whose right side has v_ebar = 1 (odd) while the left side has even
      v_ebar (or is 0, still != 1-valuation): DEAD.
  X6. (5,1) T1, ebar = E constant: then lam = alpha constant (deg ghat7 <= 1),
      and reducing the identity divided by q^9 modulo q gives EXACTLY
        sigma^2 = (8 c alpha^2 / (3 E^3)) t^15   in F := Q[y]/(q) ,
      i.e. t * (rational) is a square in F; taking N_{F/Q}:
      N(t) = q(-1)/2048 = 3315/2048 must be a square in Q -- but
      3315 * 2048 = 6789120 is not a perfect square: DEAD.
      (Symbolic layer identity verified with generic windows; norm exact.)
  X7. (6,1) T2 (the only surviving T2 cell): e = C t^6 q (deg ebar <= 0),
      sigma = s q^2, d0 = (d2^2 + s q^2)/4, d2 generic quartic.
      Layer q^13 of the identity: 2048 C^5 = 3072 c s^2, i.e.
      C^5 = -s^2/4420;  layer q^15: q | d2 * (5632 C^5 - 9216 c s^2) and
      5632(3c/2) - 9216 c = -768 c != 0, so q | d2: d2 = kappa q.
      Then d0 = mu q^2 (mu = (kappa^2+s)/4) and the identity is equivalent
      to W_n(kappa, mu, C) = 0, n = 0..6 (t-order separation): W_6 gives
      4420 C^5 = -(kappa^2-4mu)^2, W_5 gives kappa (kappa^2-4mu)^2/4 = 0 so
      kappa = 0, and then W_4 forces C^5 = -1008 mu^2/238680 while W_6
      forces C^5 = -16 mu^2/4420; 1008/238680 != 16/4420: DEAD.
  X8. (6,1) T1 REDUCTION (open cell): e = C t^6 q, d1 = q lam (X3),
      v_q(sigma) = 0; layers q^9..q^10 of the identity give exactly
        3072 C^3 sigma^2 = 8192 c t^12 lam^2  mod q^2 ,
      so (sigma/(t^6 lam))^2 = 8c/(3C^3) in F; since Gal(q) = S4, F has no
      quadratic subfield, so 8c/(3C^3) = r^2 with r rational, and then
      sigma = r t^6 lam mod q^2 (sign absorbed into r):  3 r^2 C^3 = 8 c.
      (Layer identity verified symbolically; the norm test is void here:
      N(8c/3C^3) is automatically a 4th-power-square.)  Further layer
      constraints are recorded in T5_T1_AQ12.md; the cell is NOT closed in
      this file unless X8b below is present and passing.
"""
import re
import sympy as sp

d2, d1, d0, dm1, y = sp.symbols('d2 d1 d0 dm1 y')
V4 = (d2, d1, d0, dm1)
WTS = {d2: 2, d1: 3, d0: 4, dm1: 5}
t = y + 1
qq = 2048*y**4 - 512*y**3 + 320*y**2 - 240*y + 195
c = sp.Rational(-1, 6630)
Phit = c * t**30 * qq
sig = 4*d0 - d2**2
s = sp.Symbol('s')

f31 = sp.sympify(open('f31_deg31.txt').read().strip()
                 .replace('m1', 'dm1').replace('P', 'Phi').replace('^', '**'))
Phi = sp.Symbol('Phi')
hs = {}
for m in re.finditer(r'h_(\d) \(weight (\d+), dm1-power (\d+)\) = (.+)',
                     open('f31_graded.txt').read()):
    hs[int(m.group(1))] = sp.sympify(m.group(4))
assert sorted(hs) == list(range(8))

# ---------------------------------------------------------------- check X1
fl = sp.factor_list(qq, y)
assert fl[1] == [(qq, 1)], "q NOT irreducible"
assert qq.subs(y, -1) == 3315
from sympy.polys.numberfields.galoisgroups import galois_group
qmonic = sp.Poly(sp.expand(qq/2048), y)
Gname, _alt = galois_group(qmonic, by_name=True)
assert 'S4' in str(Gname), Gname
recon = sum(Phi**f * dm1**(21 - 3*f) * hs[f] for f in range(8))
assert sp.expand(recon - f31) == 0
for f in range(8):
    for mono in sp.Poly(hs[f], *V4).monoms():
        assert sum(WTS[v]*k for v, k in zip(V4, mono)) == 20 - 2*f
# exact collapse identities WITH d1-parts:
assert sp.expand(hs[7] - 8192*d1**2) == 0
assert sp.expand(hs[6] - (-3072*sig**2 + 14336*d2*d1**2 + 8192*d1*dm1)) == 0
assert sp.expand(hs[5] - (-9216*d2*sig**2 + 2048*dm1**2 + 129024*d0*d1**2
                          - 44544*d1**2*d2**2 + 18432*d1*d2*dm1)) == 0
# structure of d1-parts and d1-free parts, all f:
for f in range(8):
    P = sp.Poly(hs[f], *V4)
    for (k, i, mm, j), co in zip(P.monoms(), P.coeffs()):
        if i >= 1:
            assert i + j >= 2          # d1-monomials get q^2 when d1=q lam, q|dm1
        else:
            assert j % 2 == 0          # d1-free: even dm1-powers
# no pure-d2 monomials after sigma-substitution (W3 fact, re-verified):
H = {f: sp.expand(hs[f].subs({d1: 0, d0: (d2**2 + s)/4})) for f in range(8)}
for f in range(8):
    if H[f] == 0:
        continue
    for (k, mm, j) in sp.Poly(H[f], d2, s, dm1).monoms():
        assert 2*k + 4*mm + 5*j == 20 - 2*f
        assert j % 2 == 0
        if j == 0:
            assert mm >= 1
assert sp.expand(H[4] - (-5184*d2**2*s**2 - 12096*s**3 + 5632*d2*dm1**2)) == 0
assert sp.expand(H[6] + 3072*s**2) == 0 and H[7] == 0
assert sp.expand(H[5] - (-9216*d2*s**2 + 2048*dm1**2)) == 0
print("X1. setup: q irreducible/Q, Gal = S4, q(-1)=3315; graded decomposition")
print("    + weights; exact h_7/h_6/h_5 with d1-parts; d1-monomials have")
print("    i+j>=2; H_4 = -5184 d2^2 s^2 - 12096 s^3 + 5632 d2 e^2          OK")

# ---------------------------------------------------------------- check X2
# caps: v_q(d1) <= 1, v_q(sigma) <= 2  (deg d1 <= 6, deg sigma <= 8, deg q = 4)
assert 6 // 4 == 1 and 8 // 4 == 2
# T1 unique-min: mu_7 = 2 v_q(d1) - 35 vs mu_f >= -5f >= -30 (f <= 6)
for vd1 in (0, 1):
    assert 2*vd1 - 35 < -30
# T2 unique-min: mu_6 = 2 v_q(sigma) - 30 vs mu_f >= -5f >= -25 (f <= 5)
for vs in (0, 1, 2):
    assert 2*vs - 30 < -25
print("X2. a_q = 2: T1 min 2v-35 <= -33 < -30, T2 min 2v-30 <= -26 < -25,")
print("    both unique => (0,2), (1,2), (2,2) dead in BOTH branches        OK")

# ---------------------------------------------------------------- check X3
# a_q = 1 T1: v_q(d1)=0 gives mu_7 = -14 < -12 <= mu_f: unique min
assert -14 < -12
for f in range(7):
    assert -2*f >= -12
# partner for mu_7 = -12: delta_f = 2f - 12 >= 0 only f = 6 (delta_6 = 0)
assert [f for f in range(7) if 2*f - 12 >= 0] == [6]
# terminal: ebar^3 g_7 = -8192 c^7 q^(7-3+2) lam^2, v_q(g_7) = 6:
for a in (2, 3, 4):
    assert 4*6 > 10 + 3*a          # g_7 != 0 impossible: T1 dead
for a in (5, 6):
    assert 4*6 <= 10 + 3*a         # (5,1), (6,1) T1 stay open here
print("X3. a_q = 1 T1: v_q(d1) = 1 and delta_6 = 0 (so v_q(sigma) = 0)")
print("    forced; v_q(g_7) = 6, 24 > 10+3a for a = 2,3,4:")
print("    T1 of (2,1), (3,1), (4,1) DEAD; (5,1), (6,1) survive to X5/X6/X8 OK")

# ---------------------------------------------------------------- check X4
# a_q = 1 T2: mu_6 = 2 w_s - 12 vs mu_f >= -2f >= -10 (f <= 5)
assert 2*0 - 12 < -10                        # w_s = 0 unique min: dead
assert [f for f in range(6) if 2*f - 10 >= 0] == [5]   # w_s=1 partner: f=5 only
# delta_5 >= 2 always (q^2 | sigma^2, q^2 | e^2): no partner: w_s = 1 dead
# w_s = 2: v_q(g_6) = 3 + 2*2 = 7 (from ebar^3 g_6 = 3072 c^6 q^7 sighat^2)
for a in (2, 3, 4, 5):
    assert 4*7 > 10 + 3*a          # T2 dead
assert 4*7 <= 10 + 3*6             # (6,1) T2 survives with w_s = 2
assert 8 - 4*2 < 4                 # deg sigma <= 8: sigma = s q^2, s constant
print("X4. a_q = 1 T2: w_s = 0 unique min; w_s = 1 needs delta_5 = 0 but")
print("    delta_5 >= 2; w_s = 2 needs 28 <= 10+3a: T2 of (2,1)-(5,1) DEAD;")
print("    (6,1) T2 reduced to sigma = s q^2, s in Q^x                     OK")

# ---------------------------------------------------------------- check X5
# (5,1) T1, deg ebar = 1.  Terminal: ebar^3 ghat7 = -8192 c^7 lam^2 with
# deg ghat7 <= 1 (cap 25 - 24), q !| lam, lam != 0, deg lam <= 2.
assert (10 + 3*5) - 4*6 == 1
# v_ebar parity: 3 + v_ebar(ghat7) in {3,4} must equal 2 v_ebar(lam):
# forced ghat7 = beta ebar, lam = alpha ebar^2 (deg lam = 2 <= 2 OK).
# d1-part of h_6 under d1 = alpha q ebar^2, dm1 = e = t^5 q ebar is q^2 ebar^3 Jt:
al, be, eb0, eb1 = sp.symbols('al be eb0 eb1')
EB = eb1*y + eb0                       # generic linear ebar
D1v = al*qq*EB**2
DM1v = t**5*qq*EB
Jt = 14336*al**2*d2*EB + 8192*al*t**5
assert sp.expand(14336*d2*D1v**2 + 8192*D1v*DM1v - qq**2*EB**3*Jt) == 0
# level-6 line divided by q^6 rearranges exactly (symbolic m, sigma):
mM, sg = sp.symbols('mM sg')
lhs = t**15*qq**6*(be*EB)                                  # t^15 g_7, g_7=q^6 be ebar
rhs = qq**3*EB**3*(qq**3*mM) + c**6*qq**6*(-3072*sg**2 + qq**2*EB**3*Jt)
assert sp.expand((lhs - rhs) - qq**6*(be*t**15*EB - EB**3*(mM + c**6*qq**2*Jt)
                                      + 3072*c**6*sg**2)) == 0
# so 3072 c^6 sigma^2 = ebar^3(m + c^6 q^2 Jt) - beta t^15 ebar: v_ebar(RHS) = 1
# (min(>=3, exactly 1) -- t, beta coprime to ebar), v_ebar(LHS) even: DEAD.
print("X5. (5,1) T1 deg ebar = 1: ghat7 = beta ebar, lam = alpha ebar^2;")
print("    d1-part of h_6 = q^2 ebar^3 Jt exactly; level-6 rearrangement")
print("    3072 c^6 sig^2 = ebar^3(...) - beta t^15 ebar: v_ebar 1 vs even:")
print("    DEAD                                                            OK")

# ---------------------------------------------------------------- check X6
# (5,1) T1, ebar = E const: lam = alpha const (deg lam^2 <= deg ghat7 <= 1).
# With generic windows d2 (deg 4), d0 (deg 8) and d1 = alpha q, e = E t^5 q,
# the term of the identity at f has the EXACT shape
#   term_f = c^f E^(21-3f) t^(105+15f) q^(21-2f) h_f(d~),
# and (exact closed forms, verified symbolically below)
#   h_7(d~) = 8192 alpha^2 q^2 ,
#   h_6(d~) = -3072 sigma^2 + q^2 J6,  J6 := 14336 alpha^2 d2 + 8192 alpha E t^5,
# while every h_f monomial has v_q >= 0 and terms f <= 5 carry q^(21-2f),
# 21 - 2f >= 11.  Hence
#   P == q^9 [ 8192 c^7 alpha^2 t^210 - 3072 c^6 E^3 t^195 sigma^2 ]  mod q^10
# and P = 0 forces  sigma^2 == (8 c alpha^2 / 3E^3) t^15  mod q.
aa = sp.symbols('aa0:5'); bb = sp.symbols('bb0:9')
E_, alc = sp.symbols('E_ alc')
D2g = sum(aa[i]*y**i for i in range(5))
D0g = sum(bb[i]*y**i for i in range(9))
SIGg = sp.expand(4*D0g - D2g**2)
subs61 = {d2: D2g, d1: sp.expand(alc*qq), d0: D0g, dm1: sp.expand(E_*t**5*qq)}
assert sp.expand(hs[7].subs(subs61) - 8192*alc**2*qq**2) == 0
J6 = 14336*alc**2*D2g + 8192*alc*E_*t**5
assert sp.expand(hs[6].subs(subs61) - (-3072*SIGg**2 + qq**2*J6)) == 0
for f in range(6):
    assert 21 - 2*f >= 11          # terms f <= 5 sit at q^11 or deeper
# 105+15*7 = 210, 105+15*6 = 195, and 21-2*7 = 7 (+2 from h_7) = 9 = 21-2*6:
assert 105 + 15*7 == 210 and 105 + 15*6 == 195 and 21 - 14 + 2 == 9 == 21 - 12
# norm obstruction: sigma^2 = (8 c alpha^2/3E^3) t^15 in F, t^15 = (t^7)^2 t,
# so t * (8c alpha^2/3E^3) in (F^x)^2; N(t rho) = rho^4 N(t) needs N(t) square:
Nt = sp.Rational(qq.subs(y, -1), 2048)
assert Nt == sp.Rational(3315, 2048)
assert not sp.sqrt(3315 * 2048).is_rational
print("X6. (5,1) T1 ebar = E: exact q-adic shapes of h_7, h_6; layer q^9")
print("    forces sigma^2 = (8c alpha^2/3E^3) t^15 in F = Q[y]/(q), i.e.")
print("    t*(rational) a square in F; N(t) = 3315/2048, 3315*2048 not a")
print("    square in Q: DEAD                                               OK")

# ---------------------------------------------------------------- check X7
# (6,1) T2: e = C t^6 q, sigma = s q^2, d0 = (d2^2 + s q^2)/4, d2 generic.
# EXACT q-adic shapes of the substituted h_f (verified symbolically):
#   h_6 = -3072 s^2 q^4                     (exactly),
#   h_5 = q^2 (2048 C^2 t^12 - 9216 s^2 q^2 d2),
#   h_4 = q^2 (5632 C^2 t^12 d2 - 5184 s^2 q^2 d2^2 - 12096 s^3 q^4),
# and every h_f, f <= 5, has v_q >= 2 (each monomial carries sigma or e^2),
# so term_f = c^f C^(21-3f) t^(126+12f) q^(21-2f) h_f sits at q^(23-2f) for
# f <= 5 and the whole identity reads EXACTLY
#   P = q^13 t^198 K2 + q^15 t^186 L d2 + q^17 R ,
#   K2 := 2048 c^5 C^8 - 3072 c^6 C^3 s^2,
#   L  := 5632 c^4 C^11 - 9216 c^5 C^6 s^2 ,
# with R a polynomial (f = 4 tail and the f <= 3 terms, floors >= 17).
C_, kap, mu = sp.symbols('C_ kap mu')
subs62 = {d2: D2g, d1: 0, d0: sp.expand((D2g**2 + s*qq**2)/4),
          dm1: sp.expand(C_*t**6*qq)}
assert hs[7].subs(subs62) == 0
assert sp.expand(hs[6].subs(subs62) + 3072*s**2*qq**4) == 0
assert sp.expand(hs[5].subs(subs62)
                 - qq**2*(2048*C_**2*t**12 - 9216*s**2*qq**2*D2g)) == 0
assert sp.expand(hs[4].subs(subs62)
                 - qq**2*(5632*C_**2*t**12*D2g - 5184*s**2*qq**2*D2g**2
                          - 12096*s**3*qq**4)) == 0
# q-exponent bookkeeping: f=6: (21-12)+4 = 13; f=5: (21-10)+2 = 13 and the
# second piece 11+4 = 15; f=4: 13+2 = 15 (t^12-piece), 13+4 = 17, 13+6 = 19;
# floors f <= 3: 21-2f+2 >= 17:
assert 9 + 4 == 13 and 11 + 2 == 13 and 11 + 4 == 15 and 13 + 2 == 15
for f in range(4):
    assert 21 - 2*f + 2 >= 17
# t-exponents: f=6: 126+72 = 198; f=5: 126+60+12 = 198 (t^12-piece) and
# 126+60 = 186 (d2-piece); f=4: 126+48+12 = 186 (d2-piece):
assert 126+72 == 198 and 126+60+12 == 198 and 126+60 == 186 and 126+48+12 == 186
K2 = 2048*c**5*C_**8 - 3072*c**6*C_**3*s**2
# layer q^13 => q | t^198 K2 => K2 = 0 <=> C^5 = (3c/2) s^2 = -s^2/4420:
assert sp.expand(K2 - 2048*c**5*C_**3*(C_**5 - sp.Rational(3,2)*c*s**2)) == 0
assert sp.Rational(3, 2)*c == sp.Rational(-1, 4420)
# layer q^15 => q | t^186 L d2; with C^5 = 3c s^2/2 (substituting C^11 =
# C (C^5)^2, C^6 = C C^5):
L = sp.expand(5632*c**4*C_*(sp.Rational(3,2)*c*s**2)**2
              - 9216*c**5*C_*(sp.Rational(3,2)*c*s**2)*s**2)
assert sp.expand(L - (-768)*sp.Rational(3,2)*c**6*C_*s**4) == 0
assert sp.expand(5632*sp.Rational(3,2)*c - 9216*c) == -768*c and -768*c != 0
# => q | d2:  d2 = kappa q.  Then d0 = mu q^2, mu = (kap^2+s)/4.
# (iii) Case A: T, Q formal; the identity collapses to W_n(kap, mu, C):
T_, Q_ = sp.symbols('T_ Q_')
subsA = {d2: kap*Q_, d1: 0, d0: mu*Q_**2, dm1: C_*T_**6*Q_}
EA = sp.expand(sum((c*T_**30*Q_)**f * (C_*T_**6*Q_)**(21-3*f)
                   * hs[f].subs(subsA) for f in range(8)))
PA = sp.Poly(EA, T_, Q_)
Ws = {}
for (et, eq), coeff in zip(PA.monoms(), PA.coeffs()):
    n = (et - 126) // 12
    assert et == 126 + 12*n and eq == 31 - 3*n   # t-orders 126+12n distinct
    Ws[n] = sp.expand(Ws.get(n, 0) + coeff)
assert sorted(Ws) == list(range(7))
sA = 4*mu - kap**2          # sigma = sA q^2
# W_6 <=> 4420 C^5 + sA^2 = 0:
assert sp.expand(Ws[6] + sp.Rational(16, 442364168356546921875)*C_**3
                 * (4420*C_**5 + sA**2)) == 0
# W_5 <=> kap (12155 C^5 + 3 sA^2) = 0; with C^5 = -sA^2/4420 this is
# kap * sA^2 / 4 * (nonzero rational) = 0  => kap = 0 (sA != 0):
w5core = 12155*(-sA**2/4420) + 3*sA**2
assert sp.simplify(w5core - sA**2*sp.Rational(1105, 4420)) == 0
assert sp.expand(Ws[5] - sp.Rational(32, 133443188041190625)*C_**6*kap
                 * (12155*C_**5 + 3*sA**2)) == 0
# W_4 at kap = 0 <=> mu (238680 C^5 + 1008 mu^2) = 0 => C^5 = -1008 mu^2/238680,
# but W_6 at kap = 0 gives C^5 = -16 mu^2/4420; the two differ:
w4k0 = sp.expand(Ws[4].subs(kap, 0))
assert sp.expand(w4k0 - sp.Rational(16, 40254355366875)*C_**9
                 * (-238680*C_**5*mu - 1008*mu**3)) == 0
assert sp.Rational(1008, 238680) != sp.Rational(16, 4420)
print("X7. (6,1) T2: layer q^13 gives C^5 = -s^2/4420; layer q^15 gives")
print("    q | d2 (coefficient -768 c s^2 != 0); then W_6, W_5 => kappa = 0,")
print("    W_4 vs W_6 => 1008/238680 != 16/4420: DEAD                      OK")

# ---------------------------------------------------------------- check X8
# (6,1) T1 reduction: e = C t^6 q, d1 = q lam (lam generic quadratic),
# d2, d0 generic.  Exact shapes:
lam0, lam1, lam2 = sp.symbols('lam0:3')
LAMg = lam2*y**2 + lam1*y + lam0
subs63 = {d2: D2g, d1: sp.expand(qq*LAMg), d0: D0g,
          dm1: sp.expand(C_*t**6*qq)}
assert sp.expand(hs[7].subs(subs63) - 8192*qq**2*LAMg**2) == 0
J6g = 14336*D2g*LAMg**2 + 8192*C_*LAMg*t**6
assert sp.expand(hs[6].subs(subs63) - (-3072*SIGg**2 + qq**2*J6g)) == 0
for f in range(6):
    assert 21 - 2*f >= 11             # f <= 5 terms sit at q^11 or deeper
# => P == q^9 c^6 t^198 N1 mod q^11,  N1 := 8192 c t^12 lam^2 - 3072 C^3 sig^2
# (t^210 = t^198 t^12), so q^2 | N1:  sigma^2 == (8c/3C^3) t^12 lam^2 mod q^2.
assert 126 + 12*7 == 210 and 126 + 12*6 == 198
# Gal(q) = S4 (X1) double-checked: resolvent cubic irreducible + disc nonsquare
b3, b2, b1, b0 = sp.Rational(-1,4), sp.Rational(5,32), sp.Rational(-15,128), \
    sp.Rational(195, 2048)
x = sp.Symbol('x')
resolvent = (x**3 - b2*x**2 + (b3*b1 - 4*b0)*x
             - (b3**2*b0 - 4*b2*b0 + b1**2))
assert sp.factor_list(resolvent, x)[1][0][1] == 1 and \
    len(sp.factor_list(resolvent, x)[1]) == 1        # irreducible cubic
disc = sp.discriminant(qq, y)
assert not sp.sqrt(disc).is_rational
# S4 root field has no proper subfield (S3 maximal in S4), so a rational
# that is a square in F is a square in Q:  8c/(3 C^3) = r^2, r in Q^x, i.e.
#     3 r^2 C^3 = 8 c ,    and    sigma = r t^6 lam + nu q^2,  nu in Q
# (sign of r absorbed; deg sigma <= 8 = deg q^2 pins the q^2-cofactor to Q).
# Sanity of the parametrization: with these shapes N1 is exactly divisible
# by q^2 when 3 r^2 C^3 = 8c:
rr, nn = sp.symbols('rr nn')
SIGp = sp.expand(rr*t**6*LAMg + nn*qq**2)
N1p = sp.expand(8192*c*t**12*LAMg**2 - 3072*C_**3*SIGp**2)
N1p = N1p.subs(rr**2, 8*c/(3*C_**3))     # impose the relation
N1p = sp.expand(N1p)
quo2, rem2 = sp.div(N1p, sp.expand(qq**2), y)
assert sp.expand(rem2) == 0
assert sp.expand(quo2 - (-3072*C_**3)*(2*rr*nn*t**6*LAMg + nn**2*qq**2)) == 0
print("X8. (6,1) T1 REDUCTION: layers q^9-q^10 give sigma^2 = (8c/3C^3)")
print("    t^12 lam^2 mod q^2; Gal(q) = S4 (resolvent cubic irreducible,")
print("    disc nonsquare) => 8c/(3C^3) = r^2 with r rational, and")
print("    sigma = r t^6 lam + nu q^2, 3 r^2 C^3 = 8c, nu in Q            OK")

# --------------------------------------------------------------- check X8b
# (6,1) T1, layered verification with the sigma-parametrization: compute
# P mod q^16 exactly (Poly ring over QQ[a*,l*,C,r,nu], y the generator),
# check q^9 | P, layers q^9/q^10 vanish identically modulo 3 r^2 C^3 = 8c,
# and the layer-q^11 class equals
#     c^6 C^3 t^198 lam (-10240 d2 lam + (8192 C - 6144 r nu) t^6)  mod q,
# i.e. (q !| lam):  10240 d2 lam == (8192 C - 6144 r nu) t^6  mod q.
av = sp.symbols('xa0:5'); lv = sp.symbols('xl0:3')
Cv, rv, nv = sp.symbols('xC xr xnu')
DOMX = sp.QQ[list(av)+list(lv)+[Cv, rv, nv]]
def PolyY(e):
    return sp.Poly(e, y, domain=DOMX)
D2x = sum(av[i]*y**i for i in range(5))
LAMx = sum(lv[i]*y**i for i in range(3))
SIGx = sp.expand(rv*t**6*LAMx + nv*qq**2)
D0x = sp.expand((SIGx + D2x**2)/4)
VALSx = {d2: PolyY(D2x), d1: PolyY(sp.expand(qq*LAMx)), d0: PolyY(D0x),
         dm1: PolyY(sp.expand(Cv*t**6*qq))}
KX = 12                       # enough for layers 9, 10, 11
qpx = {k: PolyY(sp.expand(qq**k)) for k in range(1, KX+1)}
Ptot = PolyY(0)
for f in range(8):
    base = 21 - 2*f
    need = KX - base
    if need <= 0:
        continue
    Qn = qpx[need]
    red = {v: VALSx[v].rem(Qn) for v in VALSx}
    Ph = sp.Poly(hs[f], d2, d1, d0, dm1)
    pw = {v: {0: PolyY(1)} for v in VALSx}
    def getpow(v, n):
        dd = pw[v]
        if n not in dd:
            dd[n] = (getpow(v, n-1)*red[v]).rem(Qn)
        return dd[n]
    tot = PolyY(0)
    for mono, co in zip(Ph.monoms(), Ph.coeffs()):
        term = PolyY(co)
        for v, n in zip((d2, d1, d0, dm1), mono):
            if n:
                term = (term*getpow(v, n)).rem(Qn)
        tot = tot + term
    pref = PolyY(sp.expand(c**f * t**(126+12*f))).rem(Qn) * PolyY(Cv**(21-3*f))
    Ptot = Ptot + qpx[base]*((pref*tot).rem(Qn))
Ptot = Ptot.rem(qpx[KX])
W9, rem9 = Ptot.div(qpx[9])
assert rem9.is_zero                                     # q^9 | P
def redr(e):
    """reduce even r-powers via r^2 = 8c/(3C^3), clear C-denominators."""
    p = sp.Poly(sp.expand(e), rv)
    out = 0
    for (n,), co in zip(p.monoms(), p.coeffs()):
        out += co * (rv if n % 2 else 1) * (8*c/(3*Cv**3))**(n//2)
    return sp.expand(sp.cancel(sp.together(out)))
QP1 = PolyY(qq)
Wc = W9
for layer in (9, 10):
    Wc, remq = Wc.div(QP1)
    assert redr(remq.as_expr()) == 0                    # layers 9, 10 vanish
Wc, rem11 = Wc.div(QP1)
lhs11 = redr(rem11.as_expr())
target11 = sp.rem(sp.expand(c**6*Cv**3*t**198*LAMx
                            * (-10240*D2x*LAMx + (8192*Cv - 6144*rv*nv)*t**6)),
                  qq, y)
assert sp.expand(lhs11 - sp.expand(target11)) == 0
print("X8b. layered check: q^9 | P; layers q^9, q^10 vanish iff 3r^2C^3 = 8c;")
print("     layer q^11 class = c^6 C^3 t^198 lam(-10240 d2 lam +")
print("     (8192C - 6144 r nu) t^6) mod q, forcing")
print("     10240 d2 lam == (8192C - 6144 r nu) t^6 mod q                  OK")
