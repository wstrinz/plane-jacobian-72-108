"""
t5_f37_verify.py — T5: exact verification of the bigraded (two-edge) decomposition
of f37 and every computational fact used in T5_F37_GRADED.md.

Check groups (all exact, sympy over Q; run from this directory):
  V1  Phi~ facts reused from the f31 campaign: q irreducible, q(-1)=3315,
      v_t(Phi~)=30 exactly, v_q(Phi~)=1 exactly, deg Phi~ = 34.
  V2  decomposition f37 = sum_{f=0}^7 Phi^f dm1^{p_f} h_f with
      p_f = (18,15,12,9,6,4,2,0), dm1 does NOT divide h_f (p_f maximal),
      h_f weighted-homogeneous of weight 134-17f-5p_f = (44,42,40,38,36,29,22,15),
      term counts (145,124,106,88,78,51,25,1), h_0 = f37|_{Phi=0}/dm1^18 (= h37).
  V3  Newton polygon of {(f,p_f)}: lower hull has vertices (0,18),(4,6),(7,0),
      edge slopes -3 and -2, and ALL eight points lie on the hull boundary.
  V4  the exact two-substitution form
        f37 = dm1^18 * A(d,w) + Phi^5 dm1^4 * B(d,z),
        A = sum_{f=0}^4 h_f W^f,  B = h_5 + h_6 Z + h_7 Z^2,
        w = Phi/dm1^3,  z = Phi/dm1^2,
      and the degree caps deg_y h_f(d~) <= 2*wt_f (sub2) / 3*wt_f (sub1).
  V5  cascade telescoping (Lemma A'): for every a = 0..9, with
      Delta = (3,3,3,3,2,2,2), delta_l = 30 - a*Delta_l, the recursion
        h_0 = t^{delta_1} g_1,
        u^l h_l = t^{delta_{l+1}} g_{l+1} - e^{Delta_l} g_l   (l = 1..6),
        u^7 h_7 = -e^{Delta_7} g_7
      telescopes to sum_f t^{eps_f} u^f e^{p_f} h_f = 0, eps_f = sum delta.
  V6  terminal collapse: h_7 = 221184 d1^5; h_6|_{d1=0} = -82944 d2 sigma^5
      (sigma = 4d0-d2^2); h_5|_{d1=0,sigma=0} = 131072 d2^2 dm1^5;
      d2 | h_f|_{d1=0} for EVERY f;  and the FREE FAMILY:
      f37|_{d1=0,d2=0} == 0 identically (d0, dm1, Phi all free).
  V7  sigma-locus master identity:
      f37|_{d1=0,d0=d2^2/4} = 64 d2^2 dm1^9 (32 A'^4 B' - 27 dm1^17),
      A' = 2Phi + 3 d2 dm1^3,  B' = 4Phi + 3 d2 dm1^3,
      B'-A' = 2Phi,  2A'-B' = 3 d2 dm1^3,  and the edge-quintic factorization
      2048X^5+13824X^4+36864X^3+48384X^2+31104X+7776 = 32(2X+3)^4(4X+3).
  V8  degree bookkeeping and the NO-KILL audit: cascade quotient caps
      deg g_l <= (58+3a, 58+3a, 58+3a, 58+3a, 58+2a, 48+2a, 38+2a),
      terminal degree balance 2(10-a) + (38+2a) = 58 = 28 + 30, and the f31
      starvation conditions C1'/C2' fail on ALL 21 joint strata (a + 4 a_q <= 10):
      4(7-2 a_q) <= 38+2a and 4(6-2 a_q) <= 48+2a always.  q-side caps
      delta_f <= floor(2 wt_f/4) = (22,21,20,19,18,14,11,7).
  V9  mod-q reduction on a random exact instance:
      f37(d~,Phi~) == e^18 h_0(d~)  (mod q)   [Lemma B' seed], and
      the t-strip exponents: with e = t^a e^, term f of f37(d~,Phi~) has
      t-valuation 18a + eps_f + v_t(u^f e^^{p_f} h_f)  (checked via the
      algebraic identity t^{30f+a p_f} = t^{18a} t^{eps_f} for all f, a).
"""
import re
import sympy as sp

d2, d1, d0, dm1, Phi, y, W, Z, X, t = sp.symbols('d2 d1 d0 dm1 Phi y W Z X t')
V4 = (d2, d1, d0, dm1)
WTS = {d2: 2, d1: 3, d0: 4, dm1: 5}

P_F = [18, 15, 12, 9, 6, 4, 2, 0]
WT_F = [44, 42, 40, 38, 36, 29, 22, 15]
NTERMS = [145, 124, 106, 88, 78, 51, 25, 1]
DELTA = [3, 3, 3, 3, 2, 2, 2]            # Delta_l = p_{l-1} - p_l, l = 1..7

s = open('f37_deg37.txt').read().strip()
f37 = sp.sympify(s.replace('m1', 'dm1').replace('P', 'Phi').replace('^', '**'))

hs = {}
for m in re.finditer(r'h_(\d) \(weight (\d+), dm1-power (\d+)\) = (.+)',
                     open('f37_graded.txt').read()):
    f, wt, pw, expr = int(m.group(1)), int(m.group(2)), int(m.group(3)), m.group(4)
    hs[f] = (sp.sympify(expr), wt, pw)
assert sorted(hs) == list(range(8)), "expected h_0..h_7"

# ---------------------------------------------------------------- V1
quartic = 2048*y**4 - 512*y**3 + 320*y**2 - 240*y + 195
Phit = -(y + 1)**30 * quartic / 6630
assert sp.Poly(quartic, y).is_irreducible
assert quartic.subs(y, -1) == 3315
Pt = sp.Poly(sp.expand(Phit.subs(y, t - 1)), t)
assert min(m[0] for m in Pt.monoms()) == 30           # v_t = 30 exactly
assert sp.degree(Phit, y) == 34
assert sp.rem(sp.expand(sp.numer(sp.together(Phit / (y+1)**30))), quartic, y) == 0
print("V1. Phi~ = c t^30 q, q irreducible, q(-1)=3315, deg 34      OK")

# ---------------------------------------------------------------- V2
recon = sum(Phi**f * dm1**hs[f][2] * hs[f][0] for f in range(8))
assert sp.expand(recon - f37) == 0, "DECOMPOSITION FAILS"
for f in range(8):
    h, wt, pw = hs[f]
    assert pw == P_F[f] and wt == WT_F[f] == 134 - 17*f - 5*pw
    hp = sp.Poly(h, *V4)
    assert len(hp.terms()) == NTERMS[f]
    assert min(mo[3] for mo in hp.monoms()) == 0      # dm1 does not divide h_f
    for mono in hp.monoms():
        assert sum(WTS[v]*e for v, e in zip(V4, mono)) == wt, (f, mono)
q0 = sp.cancel(f37.subs(Phi, 0) / dm1**18)
assert sp.expand(q0 - hs[0][0]) == 0
print("V2. f37 == sum Phi^f dm1^{p_f} h_f, p_f=(18,15,12,9,6,4,2,0),")
print("    h_f weighted-homog wt 134-17f-5p_f, h_0 = h37           OK")

# ---------------------------------------------------------------- V3
pts = [(f, P_F[f]) for f in range(8)]
hull = []
for pt in pts:                                        # lower hull, monotone chain
    while len(hull) >= 2:
        (x1, y1), (x2, y2) = hull[-2], hull[-1]
        if (x2 - x1)*(pt[1] - y1) - (y2 - y1)*(pt[0] - x1) <= 0:
            hull.pop()
        else:
            break
    hull.append(pt)
assert hull == [(0, 18), (4, 6), (7, 0)], hull
for f in range(5):                                    # all points ON the hull
    assert P_F[f] == 18 - 3*f
for f in range(4, 8):
    assert P_F[f] == 14 - 2*f
print("V3. Newton hull: vertices (0,18),(4,6),(7,0); slopes -3,-2;")
print("    all eight points on the boundary                        OK")

# ---------------------------------------------------------------- V4
A_head = sum(hs[f][0]*W**f for f in range(5))
B_tail = hs[5][0] + hs[6][0]*Z + hs[7][0]*Z**2
lhs = dm1**18*A_head.subs(W, Phi/dm1**3) + Phi**5*dm1**4*B_tail.subs(Z, Phi/dm1**2)
assert sp.expand(sp.together(lhs - f37)) == 0
for f in range(8):
    hp = sp.Poly(hs[f][0], *V4)
    cap2 = max(sum(2*WTS[v]*e for v, e in zip(V4, mono)) for mono in hp.monoms())
    cap3 = max(sum(3*WTS[v]*e for v, e in zip(V4, mono)) for mono in hp.monoms())
    assert cap2 == 2*WT_F[f] and cap3 == 3*WT_F[f], (f, cap2, cap3)
print("V4. f37 == dm1^18 A(w) + Phi^5 dm1^4 B(z)  (w=Phi/dm1^3, z=Phi/dm1^2);")
print("    caps deg h_f(d~) = 2wt_f (sub2) / 3wt_f (sub1)          OK")

# ---------------------------------------------------------------- V5
us, es = sp.symbols('us es')                          # u, e-hat as free symbols
gsym = sp.symbols('g1:8')
hsym = list(sp.symbols('H0:8'))
for a in range(10):
    delta = [30 - a*D for D in DELTA]                 # delta_1..delta_7
    eps = [0]
    for dl in delta:
        eps.append(eps[-1] + dl)                      # eps_0..eps_7
    # u^f h_f =: N_f (polynomial in the g's); then sum_f t^eps_f e^{p_f} N_f = 0
    N = {0: t**delta[0]*gsym[0]}
    for l in range(1, 7):
        N[l] = t**delta[l]*gsym[l] - es**DELTA[l-1]*gsym[l-1]
    N[7] = -es**DELTA[6]*gsym[6]
    total = sum(t**eps[f]*es**P_F[f]*N[f] for f in range(8))
    assert sp.expand(total) == 0, ("telescoping fails", a)
print("V5. cascade telescoping holds for every a = 0..9            OK")

# ---------------------------------------------------------------- V6
sigma = 4*d0 - d2**2
assert sp.expand(hs[7][0] - 221184*d1**5) == 0
assert sp.expand(hs[6][0].subs(d1, 0) + 82944*d2*sigma**5) == 0
assert sp.expand(hs[5][0].subs(d1, 0).subs(d0, d2**2/4) - 131072*d2**2*dm1**5) == 0
for f in range(8):
    r = sp.expand(hs[f][0].subs(d1, 0))
    assert sp.expand(r.subs(d2, 0)) == 0              # d2 | h_f|_{d1=0}
free = sp.expand(f37.subs({d1: 0, d2: 0}))
assert free == 0, "FREE FAMILY FAILS"
# numeric spot-check of the free family
import random
random.seed(37)
for _ in range(3):
    vals = {d1: 0, d2: 0, d0: sp.Rational(random.randint(-99, 99), 7),
            dm1: sp.Rational(random.randint(-99, 99), 5),
            Phi: sp.Rational(random.randint(-99, 99), 3)}
    assert f37.subs(vals) == 0
print("V6. h_7 = 221184 d1^5;  h_6|_{d1=0} = -82944 d2 sigma^5;")
print("    h_5|_{d1=0,sigma} = 131072 d2^2 dm1^5;  d2 | h_f|_{d1=0} all f;")
print("    FREE FAMILY: f37|_{d1=d2=0} == 0 identically            OK")

# ---------------------------------------------------------------- V7
Dq = d2*dm1**3
Ap, Bp = 2*Phi + 3*Dq, 4*Phi + 3*Dq
Ssig = sp.expand(f37.subs({d1: 0, d0: d2**2/4}))
assert sp.expand(Ssig - 64*d2**2*dm1**9*(32*Ap**4*Bp - 27*dm1**17)) == 0
assert sp.expand(Bp - Ap - 2*Phi) == 0
assert sp.expand(2*Ap - Bp - 3*Dq) == 0
quint = 2048*X**5 + 13824*X**4 + 36864*X**3 + 48384*X**2 + 31104*X + 7776
assert sp.expand(quint - 32*(2*X + 3)**4*(4*X + 3)) == 0
print("V7. sigma-locus: f37| = 64 d2^2 dm1^9 (32 A'^4 B' - 27 dm1^17),")
print("    A' = 2Phi~+3d2e^3, B' = 4Phi~+3d2e^3; B'-A' = 2Phi~;")
print("    edge quintic = 32(2X+3)^4(4X+3)                         OK")

# ---------------------------------------------------------------- V8
for a in range(10):
    ecap = 10 - a
    gcap = []
    prev = None
    delta = [30 - a*D for D in DELTA]
    caps_uh = [4*l + 2*WT_F[l] for l in range(8)]     # deg u^l h_l
    g = 2*WT_F[0] - delta[0]
    gcap.append(g)
    for l in range(1, 7):
        num = max(DELTA[l-1]*ecap + gcap[-1], caps_uh[l])
        gcap.append(num - delta[l])
    expect = [58+3*a, 58+3*a, 58+3*a, 58+3*a, 58+2*a, 48+2*a, 38+2*a]
    assert gcap == expect, (a, gcap)
    # terminal degree balance: deg e^^2 g_7 <= 2(10-a) + 38+2a = 58 = 28+30
    assert DELTA[6]*ecap + gcap[6] == 58 == caps_uh[7]
strata = [(a, aq) for a in range(11) for aq in range(3) if a + 4*aq <= 10]
assert len(strata) == 21
killed = [(a, aq) for (a, aq) in strata if a <= 9 and
          (4*(7 - 2*aq) > 38 + 2*a or 4*(6 - 2*aq) > 48 + 2*a)]
assert killed == [], "unexpected starvation kill"
assert [2*w//4 for w in WT_F] == [22, 21, 20, 19, 18, 14, 11, 7]
print("V8. cascade caps deg g_l = (58+3a x4, 58+2a, 48+2a, 38+2a);")
print("    terminal balance 2(10-a)+38+2a = 58 = 28+30;")
print("    NO stratum dies by C1'/C2' starvation (kill set empty)  OK")

# ---------------------------------------------------------------- V9
# strip-exponent identity: 30f + a p_f = 18a + eps_f for all f, a
asym = sp.symbols('a_')
for f in range(8):
    eps_f = sum(30 - asym*DELTA[l] for l in range(f))
    assert sp.expand(30*f + asym*P_F[f] - (18*asym + eps_f)) == 0, f
# mod-q reduction on a random exact instance (a_q = 0 generic windows):
#   6630^7 * f37(d~,Phi~) = sum_f (6630*Phi~)^f * 6630^{7-f} * (dm1^{p_f} h_f)(d~)
# and every f >= 1 summand is divisible by q (since q | Phi~), so
#   6630^7 * f37(d~,Phi~) == 6630^7 * e^18 h_0(d~)  (mod q),
# while q does NOT divide e^18 h_0(d~) generically  ->  the identity
# f37(d~,Phi~) == 0 genuinely FORCES q | h_0(d~) when q does not divide e.
random.seed(203)
def rnd(deg):
    return sum(random.randint(-9, 9)*y**i for i in range(deg + 1))
inst = {d2: rnd(4), d1: rnd(6), d0: rnd(8), dm1: rnd(10)}
P6 = sp.Poly(sp.expand(6630*Phit), y)                 # = -(y+1)^30 q, integer
qpoly = sp.Poly(quartic, y)
Tf = [sp.Poly(sp.expand((dm1**hs[f][2]*hs[f][0]).subs(inst)), y) for f in range(8)]
S = sp.Poly(0, y)
for f in range(8):
    S = S + P6**f * Tf[f] * sp.Integer(6630)**(7 - f)
assert sp.rem(S - sp.Integer(6630)**7*Tf[0], qpoly).is_zero      # == mod q
assert not sp.rem(Tf[0], qpoly).is_zero               # q does not divide e^18 h_0
assert not sp.rem(S, qpoly).is_zero                   # so f37(d~,Phi~) != 0 mod q
print("V9. 30f + a p_f == 18a + eps_f (all f, symbolic a);")
print("    f37(d~,Phi~) == e^18 h_0(d~) (mod q) exactly, and a random")
print("    instance has q coprime to h_0(d~)  ->  mod-q forcing is real  OK")

print("\nALL f37 GRADED-STRUCTURE CHECKS PASS")
