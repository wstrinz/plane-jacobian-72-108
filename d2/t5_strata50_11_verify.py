"""
t5_strata50_11_verify.py -- exact sympy verification of every computational
input used in T5_STRATA_50_11.md: the kill of the joint (t,q)-adic strata
(a_t, a_q) = (5,0) and (1,1) for the f31 branch, subcase (2).

Notation (T5_NP.md / T5_MULTIPLACE.md): t := y+1,
q := 2048y^4 - 512y^3 + 320y^2 - 240y + 195 (irreducible/Q), c := -1/6630,
Phi~ = c t^30 q, u := Phi~/t^30 = c q, e := d~_{-1} (deg <= 10, e != 0),
windows deg d~2 <= 4, d~1 <= 6, d~0 <= 8; sigma := 4 d~0 - d~2^2;
stratum (a, a_q) := (v_t(e), v_q(e)).  Identity under attack:
    0 = sum_{f=0}^{7} Phi~^f e^(21-3f) h_f(d~).

Checks (all exact, sympy over Q):
  W1. setup re-verified (self-containedness): q irreducible/Q, q(-1) = 3315,
      Phi~ = c t^30 q with v_t = 30 and v_q = 1 exact, deg Phi~ = 34,
      lc(Phi~) = -1024/3315; graded decomposition
      f31 = sum_f Phi^f dm1^(21-3f) h_f with the weight/degree caps.
  W2. collapse identities used here:
        h_7 = 8192 d1^2,
        h_6|_{d1=0} = -3072 sigma^2,
        h_5|_{d1=0} = -9216 d2 sigma^2 + 2048 dm1^2.
  W3. structure of H_f(d2, s, e) := h_f(d2, 0, (d2^2+s)/4, e) (i.e. h_f with
      d1 = 0 and sigma = s substituted): every monomial d2^k s^m e^j has
      weight 2k + 4m + 5j = 20 - 2f, j even, j = 4 only at f = 0, and
      j = 0 ==> m >= 1 (no pure-d2 monomials); consequently, with
      deg d2 <= 4, deg e = 5, s constant:  deg_y H_f <= 32 - 4f  (f <= 5);
      H_6 = -3072 s^2 exactly, H_7 = 0; s^0-parts match the sigma-locus
      collapse constants (c_0..c_5) = (-2560,-8192,-7168,2048,5632,2048).
  W4. stratum arithmetic for (a,aq) in {(5,0),(1,1)}: v = 30-3a, cascade cap
      10+3a; C1: 4(7-3aq) > 10+3a (forces d~1 = 0 and g_7 = 0, Prop 1 of
      T5_MULTIPLACE.md); level-6 valuation squeeze: 4(6-3aq) <= 10+3a
      (T2 open) but 4((6-3aq)+2) > 10+3a (forces v_q(sigma) = 0), and
      (10+3a) - 4(6-3aq) = 1 (so deg ghat <= 1).
  W5. cascade telescoping re-verified (Lemma A sufficiency, symbolic), plus
      the constant  -9216 c^5 / (3072 c^6) = -3/c = 19890  and the exact
      level-5 rearrangement identities:
      (5,0):  t^15 q^6 gh - [eh^3 g5 + c^5 q^5 (-9216 d2 s2 + 2048 t^10 eh^2)]
              = q^5 t^10 (t^5 q gh - 2048 c^5 eh^2) - eh^3 (g5 + 19890 q^5 d2 gh)
              where s2 := eh^3 gh / (3072 c^6)   [terminal relation],
      (1,1):  t^27 q^3 gh - [(q eb)^3 g5 + c^5 q^5 (-9216 d2 s2' + 2048 t^2 q^2 eb^2)]
              = q^3 [ t^2 (t^25 gh - 2048 c^5 q^4 eb^2) - eb^3 (g5 + 19890 q^2 d2 gh) ]
              where s2' := eb^3 gh / (3072 c^6).
  W6. endgame (5,0), symbolic coefficients: d2 = a4 y^4 + ... + a0, e = C t^5,
      H_f := h_f(d2, 0, (d2^2+s)/4, e):  H_7 = 0;
      T_6 := Phi~^6 e^3 H_6 = -3072 s^2 c^6 C^3 t^195 q^6, deg 219,
      lc = -3072 (1024/3315)^6 s^2 C^3;  and for f = 0..5 the formal degree
      34 f + 5(21-3f) + deg_y H_f <= 212 < 219  (valid for every
      specialization: formal generic degree is an upper bound).
  W7. endgame (1,1): same with e = C t q:  T_6 = -3072 s^2 c^6 C^3 t^183 q^9,
      deg 219, lc = -3072 (2048^9/6630^6) s^2 C^3; per-term caps <= 212.
  W8. numeric sanity for both strata: a concrete instance (random d2, C, s)
      gives P := sum_f Phi~^f e^(21-3f) H_f  of degree exactly 219 with the
      predicted leading coefficient (hence P != 0), and P agrees with direct
      evaluation of the original graded sum at the substituted window tuple.
"""
import re
import sympy as sp

d2, d1, d0, dm1, Phi, y, X = sp.symbols('d2 d1 d0 dm1 Phi y X')
V4 = (d2, d1, d0, dm1)
WTS = {d2: 2, d1: 3, d0: 4, dm1: 5}
t = y + 1
qq = 2048*y**4 - 512*y**3 + 320*y**2 - 240*y + 195
c = sp.Rational(-1, 6630)
Phit = c * t**30 * qq
s = sp.Symbol('s')          # the constant value of sigma

f31 = sp.sympify(open('f31_deg31.txt').read().strip()
                 .replace('m1', 'dm1').replace('P', 'Phi').replace('^', '**'))
hs = {}
for m in re.finditer(r'h_(\d) \(weight (\d+), dm1-power (\d+)\) = (.+)',
                     open('f31_graded.txt').read()):
    hs[int(m.group(1))] = sp.sympify(m.group(4))
assert sorted(hs) == list(range(8))

# ---------------------------------------------------------------- check W1
fl = sp.factor_list(qq, y)
assert fl[1] == [(qq, 1)], "q NOT irreducible"
assert qq.subs(y, -1) == 3315
assert sp.degree(Phit, y) == 34
assert sp.LC(sp.Poly(sp.expand(Phit), y)) == sp.Rational(-1024, 3315)
assert sp.rem(sp.expand(t**30), qq, y) != 0            # q !| t^30
Pq, Rq = sp.div(sp.expand(Phit), qq, y)
assert sp.expand(Rq) == 0 and sp.rem(Pq, qq, y) != 0   # q || Phi~
recon = sum(Phi**f * dm1**(21 - 3*f) * hs[f] for f in range(8))
assert sp.expand(recon - f31) == 0
for f in range(8):
    for mono in sp.Poly(hs[f], *V4).monoms():
        assert sum(WTS[v]*k for v, k in zip(V4, mono)) == 20 - 2*f
print("W1. setup: q irreducible/Q, q(-1)=3315; Phi~ = c t^30 q, v_t=30,")
print("    v_q=1, deg=34, lc=-1024/3315; graded decomposition + weights  OK")

# ---------------------------------------------------------------- check W2
sig = 4*d0 - d2**2
assert sp.expand(hs[7] - 8192*d1**2) == 0
assert sp.expand(hs[6].subs(d1, 0) + 3072*sig**2) == 0
assert sp.expand(hs[5].subs(d1, 0) - (-9216*d2*sig**2 + 2048*dm1**2)) == 0
print("W2. h_7 = 8192 d1^2;  h_6|_{d1=0} = -3072 sigma^2;")
print("    h_5|_{d1=0} = -9216 d2 sigma^2 + 2048 dm1^2                   OK")

# ---------------------------------------------------------------- check W3
CF = {0: -2560, 1: -8192, 2: -7168, 3: 2048, 4: 5632, 5: 2048}
H = {f: sp.expand(hs[f].subs({d1: 0, d0: (d2**2 + s)/4})) for f in range(8)}
assert H[7] == 0
assert sp.expand(H[6] + 3072*s**2) == 0
assert sp.expand(H[5] - (-9216*d2*s**2 + 2048*dm1**2)) == 0
for f in range(7):
    if H[f] == 0:
        continue
    P = sp.Poly(H[f], d2, s, dm1)
    degs = []
    for (k, mm, j) in P.monoms():
        assert 2*k + 4*mm + 5*j == 20 - 2*f          # weight
        assert j % 2 == 0 and j <= 4                 # e-powers even, <= 4
        if j == 4:
            assert f == 0                            # e^4 only in H_0
        if j == 0:
            assert mm >= 1                           # no pure-d2 monomials
        degs.append(4*k + 5*j)      # deg_y bound with deg d2 <= 4, deg e = 5
    if f <= 5:
        assert max(degs) <= 32 - 4*f, (f, max(degs))
# s^0-parts = sigma-locus collapse:
for f in range(1, 6):
    assert sp.expand(H[f].subs(s, 0) - CF[f]*d2**(5-f)*dm1**2) == 0
assert sp.expand(H[0].subs(s, 0) - (-2560*d2**5*dm1**2 - 6561*dm1**4)) == 0
print("W3. H_f(d2,s,e) monomial structure: weights OK, e-powers even,")
print("    e^4 only f=0, no pure-d2 terms; deg_y H_f <= 32-4f (f<=5,")
print("    deg d2<=4, deg e=5); H_6 = -3072 s^2, H_7 = 0; s^0-parts = c_f  OK")

# ---------------------------------------------------------------- check W4
for (a, aq) in [(5, 0), (1, 1)]:
    v = 30 - 3*a
    cap = 10 + 3*a
    assert a + 4*aq <= 10
    assert 4*(7 - 3*aq) > cap          # C1: T1 closed => d~1 = 0, g_7 = 0
    base = 6 - 3*aq                    # v_q(g_6) = base + 2 v_q(sigma)
    assert 4*base <= cap               # T2 not closed by degrees alone
    assert 4*(base + 2) > cap          # v_q(sigma) >= 1 impossible
    assert cap - 4*base == 1           # g_6 = q^base ghat with deg ghat <= 1
print("W4. (5,0): v=15, cap 25, C1 28>25, 24<=25<4*8, deg ghat<=1;")
print("    (1,1): v=27, cap 13, C1 16>13, 12<=13<4*5, deg ghat<=1        OK")

# ---------------------------------------------------------------- check W5
# (i) Lemma A telescoping (sufficiency), symbolic:
T_, U_, E_ = sp.symbols('T_ U_ E_')
g = [None] + list(sp.symbols('g1:8'))
hh = {0: T_*g[1]}
for l in range(1, 7):
    hh[l] = (T_*g[l+1] - E_**3*g[l]) / U_**l
hh[7] = -E_**3*g[7] / U_**7
S = sum(T_**f * U_**f * E_**(21 - 3*f) * hh[f] for f in range(8))
assert sp.expand(sp.cancel(S)) == 0
# (ii) the level-5 constant:
assert sp.Rational(-9216, 3072)/c == 19890
# (iii) exact rearrangement identities (symbolic eh, eb, gh, g5, dd2):
eh, eb, gh, g5s, dd2 = sp.symbols('eh eb gh g5s dd2')
# stratum (5,0): terminal relation sigma^2 = eh^3 gh/(3072 c^6)
s2 = eh**3*gh/(3072*c**6)
lhs = t**15 * qq**6 * gh
rhs = eh**3*g5s + c**5*qq**5*(-9216*dd2*s2 + 2048*t**10*eh**2)
X50 = t**5*qq*gh - 2048*c**5*eh**2
G50 = g5s + 19890*qq**5*dd2*gh
assert sp.simplify(lhs - rhs - (qq**5*t**10*X50 - eh**3*G50)) == 0
# stratum (1,1): ehat = q*eb, e^2 = t^2 q^2 eb^2, sigma^2 = eb^3 gh/(3072 c^6)
s2b = eb**3*gh/(3072*c**6)
lhs2 = t**27 * qq**3 * gh
rhs2 = (qq*eb)**3*g5s + c**5*qq**5*(-9216*dd2*s2b + 2048*t**2*qq**2*eb**2)
X11 = t**25*gh - 2048*c**5*qq**4*eb**2
G11 = g5s + 19890*qq**2*dd2*gh
assert sp.simplify(lhs2 - rhs2 - qq**3*(t**2*X11 - eb**3*G11)) == 0
print("W5. Lemma A telescoping OK; -9216 c^5/(3072 c^6) = 19890; level-5")
print("    rearrangements exact for (5,0) and (1,1)                      OK")

# ---------------------------------------------------------------- check W6
a0, a1, a2, a3, a4, C = sp.symbols('a0 a1 a2 a3 a4 C')
D2gen = a4*y**4 + a3*y**3 + a2*y**2 + a1*y + a0

def endgame(E, name, T6_target, lc_target):
    """E = the forced e; verify T_6 exactly and per-term degree caps."""
    Hval = {f: sp.expand(hs[f].subs({d1: 0, d2: D2gen,
                                     d0: (D2gen**2 + s)/4, dm1: E}))
            for f in range(8)}
    assert Hval[7] == 0
    assert sp.expand(Hval[6] + 3072*s**2) == 0     # = -3072 s^2, deg 0
    T6 = sp.expand(Phit**6 * E**3 * (-3072*s**2))
    assert sp.expand(T6 - T6_target) == 0
    P6 = sp.Poly(T6, y)
    assert P6.degree() == 219
    assert sp.expand(P6.LC() - lc_target) == 0
    degs = {}
    for f in range(6):
        dH = sp.degree(Hval[f], y)     # formal degree: upper bound for
        degs[f] = 34*f + 5*(21 - 3*f) + dH  # every specialization
        assert degs[f] <= 212, (name, f, degs[f])
    return Hval, degs

E50 = sp.expand(C * t**5)
T6_50 = sp.expand(-3072*s**2*c**6*C**3 * t**195 * qq**6)
lc50 = sp.Rational(-3072) * sp.Rational(1024, 3315)**6 * s**2 * C**3
Hval50, degs50 = endgame(E50, "(5,0)", T6_50, lc50)
print("W6. (5,0) endgame: T_6 = -3072 s^2 c^6 C^3 t^195 q^6, deg 219,")
print("    lc = -3072 (1024/3315)^6 s^2 C^3;  per-term degrees f=0..5:")
print("    ", degs50, " all <= 212 < 219                                 OK")

# ---------------------------------------------------------------- check W7
E11 = sp.expand(C * t * qq)
T6_11 = sp.expand(-3072*s**2*c**6*C**3 * t**183 * qq**9)
lc11 = sp.Rational(-3072) * sp.Rational(2048)**9 / sp.Rational(6630)**6 \
    * s**2 * C**3
Hval11, degs11 = endgame(E11, "(1,1)", T6_11, lc11)
print("W7. (1,1) endgame: T_6 = -3072 s^2 c^6 C^3 t^183 q^9, deg 219,")
print("    lc = -3072 (2048^9/6630^6) s^2 C^3;  per-term degrees f=0..5:")
print("    ", degs11, " all <= 212 < 219                                 OK")

# ---------------------------------------------------------------- check W8
import random
random.seed(11)
sub = {a4: 2, a3: -1, a2: 3, a1: 1, a0: -2, C: 5, s: 7}
for (name, E, lct) in [("(5,0)", E50, lc50), ("(1,1)", E11, lc11)]:
    D2n = D2gen.subs(sub)
    En = E.subs(sub)
    Hn = {f: sp.expand(hs[f].subs({d1: 0, d2: D2n,
                                   d0: (D2n**2 + sub[s])/4, dm1: En}))
          for f in range(8)}
    P = sp.expand(sum(Phit**f * En**(21 - 3*f) * Hn[f] for f in range(8)))
    PP = sp.Poly(P, y)
    assert PP.degree() == 219
    assert PP.LC() == lct.subs(sub)
    assert P != 0
    # cross-check against direct evaluation of the graded sum at a point:
    y0 = sp.Rational(3, 2)
    direct = sum((Phit**f * En**(21 - 3*f)
                  * hs[f].subs({d1: 0, d2: D2n, d0: (D2n**2 + sub[s])/4,
                                dm1: En})).subs(y, y0) for f in range(8))
    assert sp.nsimplify(P.subs(y, y0) - direct) == 0
print("W8. numeric instances: P has degree exactly 219 with the predicted")
print("    nonzero lc in both strata; matches direct graded evaluation    OK")

print("\nALL STRATA-(5,0)/(1,1) CHECKS PASS")
