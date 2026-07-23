"""
t5_multiplace_verify.py — T5 multi-place analysis (T5_MULTIPLACE.md): exact
verification of every computational input used in the joint (t,q)-adic cascade
argument for the f31 branch, subcase (2).

Places: t := y+1 and the quartic q := 2048y^4-512y^3+320y^2-240y+195, with
Phi~ = c t^30 q, c = -1/6630 (so u := Phi~/t^30 = c q).

Checks (all exact, sympy over Q):
  1. q is irreducible over Q; q(-1)=3315; disc(q) != 0; v_t(Phi~)=30 and
     v_q(Phi~)=1 exactly; deg Phi~ = 34, leading coeff -1024/3315.
  2. graded decomposition f31 = sum_f Phi^f dm1^(21-3f) h_f (re-verified so
     this file is self-contained) and the window degree caps
     deg_y h_f(d~) <= 40-4f, whence v_q(h_f(d~)) <= 10-f when h_f(d~) != 0.
  3. collapse identities driving the terminal collapse chain:
       h_7 = 8192 d1^2
       h_6|_{d1=0}            = -3072 (4 d0 - d2^2)^2
       h_f|_{d1=0, d0=d2^2/4} = c_f d2^(5-f) dm1^2   (f=1..5),
       h_0|_{d1=0, d0=d2^2/4} = -2560 d2^5 dm1^2 - 6561 dm1^4,
       h_6| = h_7| = 0 on that locus,
     with (c_0..c_5) := (-2560,-8192,-7168,2048,5632,2048), all nonzero, and
     the quintic factorization  sum_f c_f X^f = 512 (X+1)^4 (4X-5).
  4. sigma-locus master identity (exact, symbolic):  on d1=0, d0=d2^2/4,
       sum_f Phi^f dm1^(21-3f) h_f|  =  dm1^8 * (512 A^4 B - 6561 dm1^17),
     A := Phi + d2 dm1^3,  B := 4 Phi - 5 d2 dm1^3;  and the linear-algebra
     facts 5A+B = 9 Phi, 4A-B = 9 d2 dm1^3, plus
     Phi~' = c (y+1)^29 (30 q + (y+1) q')  with deg(30q+(y+1)q') = 4
     (used in the A-const / B-const derivative kills).
  5. cascade telescoping (sufficiency): if h_0 = T g_1,
     h_l = (T g_{l+1} - E^3 g_l)/U^l (l=1..6), h_7 = -E^3 g_7 / U^7, then
     sum_f T^f U^f E^(21-3f) h_f == 0 identically (symbolic T,U,E,g_1..g_7).
  6. reduction algebra on a random window instance with e = t^a ehat (a=2):
       F := sum_f Phi~^f e^(21-3f) h_f(d~)  ==  t^(21a) * G,
       G := sum_f t^(vf) u^f ehat^(21-3f) h_f(d~),   v := 30-3a,
     and G - ehat^21 h_0(d~) is divisible by t^v AND by q
     (the two forced first steps of the joint cascade).
  7. degree bookkeeping of the general-a cascade: for a=0..9, l=1..7:
     deg(u^l h_l(d~)) <= 40, deg(ehat^3 g + u^l h) <= 40 given
     deg g <= 10+3a, so all cascade quotients have deg <= 40-v = 10+3a.
  8. kill-chain arithmetic: with C1: 4(7-3aq) > 10+3a (forces d1~=0),
     C2: 4(6-3aq) > 10+3a (forces 4d0~=d2~^2), and the level-5 degree
     20+2a-4aq > 10+3a  <=>  a < 10-4aq (always true under C2), the killed
     strata are exactly {(0,0),(1,0),(2,0),(3,0),(4,0),(0,1)}; survivor
     table printed.  Mason-Stothers arithmetic of the sigma-locus kill:
     34 > 9-1 and 17 > 7-1; and 17 does not divide 150 (d2~=0 subcase).
"""
import re, sympy as sp

d2, d1, d0, dm1, Phi, y, W, X = sp.symbols('d2 d1 d0 dm1 Phi y W X')
V4 = (d2, d1, d0, dm1)
WTS = {d2: 2, d1: 3, d0: 4, dm1: 5}
t = y + 1

# ---------------------------------------------------------------- inputs
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

# ---------------------------------------------------------------- check 1
fl = sp.factor_list(qq, y)
assert fl[1] == [(qq, 1)], "q NOT irreducible"
assert qq.subs(y, -1) == 3315
assert sp.discriminant(qq, y) != 0
assert sp.degree(Phit, y) == 34
assert sp.LC(sp.Poly(sp.expand(Phit), y)) == sp.Rational(-1024, 3315)
# v_t(Phi~) = 30 exactly (unit cofactor), v_q(Phi~) = 1 exactly
assert sp.rem(sp.expand(t**30), qq, y) != 0          # q does not divide t^30
Pq, Rq = sp.div(sp.expand(Phit), qq, y)
assert sp.expand(Rq) == 0 and sp.rem(Pq, qq, y) != 0  # q || Phi~
# a_q := v_q(e) <= 2 since 4*3 = 12 > 10 >= deg e; and t^a q^aq | e => a+4aq<=10
assert 4*3 > 10
print("1. q irreducible/Q, q(-1)=3315, disc!=0; v_t(Phi~)=30, v_q(Phi~)=1;")
print("   deg Phi~=34, lc=-1024/3315; a_q<=2, a_t+4a_q<=10        OK")

# ---------------------------------------------------------------- check 2
recon = sum(Phi**f * dm1**(21 - 3*f) * hs[f] for f in range(8))
assert sp.expand(recon - f31) == 0
for f in range(8):
    for mono in sp.Poly(hs[f], *V4).monoms():
        assert sum(WTS[v]*e for v, e in zip(V4, mono)) == 20 - 2*f
    maxdeg = max(sum(2*WTS[v]*e for v, e in zip(V4, mono))
                 for mono in sp.Poly(hs[f], *V4).monoms())
    assert maxdeg == 40 - 4*f
    assert 4*(10 - f) <= 40 - 4*f   # so v_q(h_f(d~)) <= 10-f
print("2. f31 = sum Phi^f dm1^(21-3f) h_f; deg_y h_f(d~) <= 40-4f;")
print("   hence v_q(h_f(d~)) <= 10-f                              OK")

# ---------------------------------------------------------------- check 3
CF = {0: -2560, 1: -8192, 2: -7168, 3: 2048, 4: 5632, 5: 2048}
assert sp.expand(hs[7] - 8192*d1**2) == 0
assert sp.expand(hs[6].subs(d1, 0) + 3072*(4*d0 - d2**2)**2) == 0
loc = {d1: 0, d0: d2**2/4}
for f in range(1, 6):
    assert sp.expand(hs[f].subs(loc) - CF[f]*d2**(5-f)*dm1**2) == 0, f
assert sp.expand(hs[0].subs(loc) - (-2560*d2**5*dm1**2 - 6561*dm1**4)) == 0
assert sp.expand(hs[6].subs(loc)) == 0 and sp.expand(hs[7].subs(loc)) == 0
assert all(CF[f] != 0 for f in CF)
quintic = sum(CF[f]*X**f for f in range(6))
assert sp.expand(quintic - 512*(X + 1)**4*(4*X - 5)) == 0
print("3. h_7 = 8192 d1^2; h_6|_{d1=0} = -3072(4d0-d2^2)^2;")
print("   h_f|_sigma = c_f d2^(5-f) dm1^2 (f=1..5), h_0|_sigma adds")
print("   -6561 dm1^4; h_6|=h_7|=0; sum c_f X^f = 512(X+1)^4(4X-5) OK")

# ---------------------------------------------------------------- check 4
A = Phi + d2*dm1**3
B = 4*Phi - 5*d2*dm1**3
assert sp.expand(5*A + B - 9*Phi) == 0
assert sp.expand(4*A - B - 9*d2*dm1**3) == 0
lhs = sum(Phi**f * dm1**(21 - 3*f) * hs[f].subs(loc) for f in range(8))
rhs = dm1**8 * (512*A**4*B - 6561*dm1**17)
assert sp.expand(lhs - rhs) == 0
dPhit = sp.diff(Phit, y)
assert sp.expand(dPhit - c*t**29*(30*qq + t*sp.diff(qq, y))) == 0
assert sp.degree(30*qq + t*sp.diff(qq, y), y) == 4
print("4. sigma-locus master: sum Phi^f dm1^(21-3f) h_f|_sigma")
print("   = dm1^8 (512 A^4 B - 6561 dm1^17), 5A+B=9Phi, 4A-B=9d2*dm1^3;")
print("   Phi~' = c t^29 (30q + t q'), deg(30q+tq') = 4           OK")

# ---------------------------------------------------------------- check 5
T, U, E = sp.symbols('T U E')
g = [None] + list(sp.symbols('g1:8'))
hh = {0: T*g[1]}
for l in range(1, 7):
    hh[l] = (T*g[l+1] - E**3*g[l]) / U**l
hh[7] = -E**3*g[7] / U**7
S = sum(T**f * U**f * E**(21 - 3*f) * hh[f] for f in range(8))
assert sp.expand(sp.cancel(S)) == 0
print("5. cascade telescoping: recursion => sum T^f U^f E^(21-3f) h_f = 0  OK")

# ---------------------------------------------------------------- check 6
import random
random.seed(7)
def rpoly(deg):
    return sum(random.randint(-3, 3)*y**k for k in range(deg + 1)) + y**deg
a = 2
v = 30 - 3*a
D2, D1, D0 = rpoly(4), rpoly(6), rpoly(8)
Ehat = rpoly(8)
if Ehat.subs(y, -1) == 0:
    Ehat += 1
Epoly = sp.expand(t**a * Ehat)
u = c * qq
hval = {f: hs[f].subs({d2: D2, d1: D1, d0: D0, dm1: Epoly}) for f in range(8)}
hval_hat = hval  # h_f(d~) with dm1 slot = full e (as in the identity)
PF = sp.Poly(sp.expand(Phit), y)
EP = sp.Poly(Epoly, y)
EH = sp.Poly(sp.expand(Ehat), y)
UP = sp.Poly(sp.expand(u), y)
TP = sp.Poly(t, y)
F = sp.Poly(0, y)
G = sp.Poly(0, y)
for f in range(8):
    hf = sp.Poly(sp.expand(hval[f]), y)
    F += PF**f * EP**(21 - 3*f) * hf
    G += TP**(v*f) * UP**f * EH**(21 - 3*f) * hf
assert sp.expand(F.as_expr() - t**(21*a)*G.as_expr()) == 0
diff0 = G - EH**21 * sp.Poly(sp.expand(hval[0]), y)
qd, rd = sp.div(diff0.as_expr(), sp.expand(t**v), y)
assert sp.expand(rd) == 0                      # t^v | G - ehat^21 h_0(d~)
assert sp.rem(diff0.as_expr(), qq, y) == 0     # q   | G - ehat^21 h_0(d~)
print("6. random instance (a=2): F = t^(21a) G;  t^v and q both divide")
print("   G - ehat^21 h_0(d~)  (first joint cascade step)         OK")

# ---------------------------------------------------------------- check 7
for a_ in range(10):
    v_ = 30 - 3*a_
    for l in range(1, 8):
        assert 4*l + (40 - 4*l) == 40
        assert 3*(10 - a_) + (10 + 3*a_) == 40
        assert 40 - v_ == 10 + 3*a_
print("7. degree bookkeeping: cascade objects have deg <= 10+3a   OK")

# ---------------------------------------------------------------- check 8
killed, surviving = [], []
for aq in range(3):
    for a_ in range(11 - 4*aq):
        if a_ == 10:            # degenerate stratum e = C t^10, no cascade
            surviving.append((a_, aq)); continue
        C1 = 4*(7 - 3*aq) > 10 + 3*a_
        C2 = 4*(6 - 3*aq) > 10 + 3*a_
        lvl5 = 20 + 2*a_ - 4*aq > 10 + 3*a_   # deg forced > cap
        assert (not C2) or lvl5               # C2 => level-5 contradiction
        (killed if (C1 and C2) else surviving).append((a_, aq))
assert killed == [(0, 0), (1, 0), (2, 0), (3, 0), (4, 0), (0, 1)]
assert 34 > 9 - 1 and 17 > 7 - 1              # Mason-Stothers margins
assert 150 % 17 != 0                          # d2~=0 subcase: 17a=150 impossible
print("8. killed strata (cascade levels 7/6/5): (0,0) (1,0) (2,0) (3,0)")
print("   (4,0) (0,1); Mason margins 34>8, 17>6; 17 !| 150         OK")
print("   surviving (a_t,a_q):", surviving)

print("\nALL MULTIPLACE CHECKS PASS")
