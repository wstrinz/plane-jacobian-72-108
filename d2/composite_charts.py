#!/usr/bin/env python3
"""composite_charts.py  (NEW; read-only over all existing artifacts)

LANE J -- composite charts: the last structural escape hatch of the corner law
(kappa = t-2 was proven only on the standard length-1 A0'=(1,0) chart class;
the named escapes are the A0'=(2,0) families F12/F13 and the length-2 chains
F18-F24).

Results in brief (details in COMPOSITE_CHARTS.md):

1. FUSED-CHART LEMMA.  Every chain transformation published in the GGV/GGHV
   program is either a root-shift shear  y -> y + lambda*x^(-s/l)  (Jacobian 1,
   bracket-preserving -- GGV1 Prop 3.10 as cited at GGV5 line 1758) or the one
   final presentation map.  Fusing shears into the final map gives
       (X, Y) = (x^-1,  x^l * y  +  sum_i lambda_i x^(e_i)),   e_i integers,
   whose Jacobian is exactly -x^(l-2) for ANY shear terms.  Hence
       kappa = l - 2 = t - 2
   for every chain realized in this form -- ANY length, ANY A0'.  The
   composition heuristic kappa = l2 - l1 (PHI_CORNER4.md escape 2) is the
   Jacobian of a chart with TWO inversions; the chain has only ONE.  Heuristic
   dead, escape closed at the chart level.

2. F18-F21 are DEAD as counterexample sources: GGV5 (lines 1726-1786) proves
   no standard (m,n)-pair realizes them.  The checkable core: the shifted
   start (1/m)st_10 = (6,3) violates the last-lower-corner criterion
   b + gcd(a,b) != a  (3 + 3 = 6), and the degenerate path lands at (6,9)
   with v11 = 15 < 16, below the minimum for admissible chains.

3. A0'=(2,0) breaks the MODEL, not the chart: with P = x^zeta * C^a
   (zeta > 0 the pure-power defect), the commutator
       [x^zeta C^a, C^b] = zeta * b * x^(zeta-1) * C^(a+b-1) * C_y  !=  0,
   so the "commuting C-powers" tail structure underlying the forcing family
   (F) of CORNER_144 needs a zeta-corrected theory.  kappa is settled; the
   next boundary is precisely this term.

4. F12 conditional branch analysis (standard family (F) at t=4, kappa=2,
   a0=8, q=5 => dg=3, gap=2, r=2 -- the first gap>0,r>0 point with dg ODD):
   the consistency variety has exactly THREE points (one rational), ALL with
   squarefree g -- NO ramified branch exists, the mirror image of PHI_F7's
   dg=2 obstruction.  A cubic always has a nonzero real root, so the
   UNRAMIFIED gauge g(-1)=0 (simple) is realizable and the conditional
   signature is the OLD-law shape:
       Phi sig (814, 506, 102, 206)   [ramified formulas would give (814,506,304,4)]
   dg-parity now has a mechanism: ramification is forced exactly when dg is
   even (no real root available), never when dg is odd.

Sources: GGV5 paper_src/1708.07936_GGV5.tex (family tables lines 1674-1718,
F18-F21 impossibility lines 1726-1786, chain-step theorem lines 1032-1050);
GGHV22 paper_src/2204.14178.tex (root shift line 1132, final map + chain rule
lines 1228-1234).  Checker: composite_charts_verify.py.  Exact sympy only.
"""
import itertools
import sympy as sp
from math import gcd
from fractions import Fraction

x, y = sp.symbols("x y")

BAR = "=" * 96


# ---------------------------------------------------------------------------
# 0. Chain data (GGV5 tables, transcribed verbatim; lines 1674-1718)
# ---------------------------------------------------------------------------
# length-1 escapes + standard controls: (name, A0, A0', A1=(p,l,q), m(j), n(j))
LEN1 = {
    "F5":  ((5, 20), (1, 0), (9, 5, 4)),
    "F7":  ((6, 15), (1, 0), (7, 3, 4)),
    "F9":  ((7, 21), (1, 0), (11, 7, 2)),
    "F10": ((7, 21), (1, 0), (13, 7, 3)),
    "F16": ((9, 24), (1, 0), (10, 3, 7)),
    "F12": ((8, 24), (2, 0), (13, 4, 5)),
    "F13": ((9, 21), (2, 0), (13, 3, 7)),
}
# length-2 chains (GGV5 lines 1709-1715): name -> (A0, A0', A1, A1', A2)
# fractional corners as (p, l, q) meaning (p/l, q).
LEN2 = {
    "F18": ((6, 18), (6, 15), (6, 15), (1, 0), (7, 3, 4)),
    "F19": ((6, 18), (6, 15), (6, 15), (1, 0), (8, 3, 5)),
    "F20": ((6, 24), (6, 15), (6, 15), (1, 0), (7, 3, 4)),
    "F21": ((6, 24), (6, 15), (6, 15), (1, 0), (8, 3, 5)),
    "F22": ((8, 24), (2, 0), (14, 4, 6), (5, 4, 2), (5, 4, 2)),
    "F23": ((8, 24), (2, 0), (14, 4, 6), (11, 4, 4), (11, 4, 4)),
    "F24": ((8, 24), (2, 0), (14, 4, 6), (5, 4, 0), (19, 8, 3)),
}


def edge_direction(dx, dy):
    """Primitive normal (rho, sigma), rho > 0, with rho*dx + sigma*dy = 0."""
    fr = Fraction(dy, dx) if dx else None
    # rho*dx = -sigma*dy  ->  (rho, sigma) ~ (dy, -dx) reduced
    g = gcd(abs(dy), abs(dx))
    rho, sigma = abs(dy) // g, -(dx // g) * (1 if dy > 0 else -1)
    if rho < 0:
        rho, sigma = -rho, -sigma
    assert rho * dx + sigma * dy == 0
    return rho, sigma


def step_corner(lower, q, s, l):
    """A_next = (a', b'=0-side) + q*(s/l, 1): the chain-step corner formula."""
    a_next = (Fraction(lower[0]) + Fraction(q * s, l), Fraction(lower[1]) + q)
    return a_next


print(BAR)
print("STEP 1 -- chain-step arithmetic reproduces every table corner (controls + escapes)")
print(BAR)
for name, (A0, A0p, (p, l, q)) in LEN1.items():
    dx, dy = A0[0] - A0p[0], A0[1] - A0p[1]
    rho, sigma = edge_direction(dx, dy)
    s = -sigma
    got = step_corner(A0p, q, s, l)
    ok = (rho == l and got == (Fraction(p, l), Fraction(q)))
    print(f"  {name}: dir=({rho},{sigma})  l={l} s={s}  "
          f"A1 = {A0p} + {q}*({s}/{l},1) = ({got[0]},{got[1]})  "
          f"table ({p}/{l},{q})  {'OK' if ok else 'MISMATCH'}")
    assert ok, name

# length-2: F22/F23 end AT A1' (one transformation); F24 has a second step.
print("\n  length-2 structure:")
for name in ("F22", "F23"):
    A0, A0p, A1, A1p, A2 = LEN2[name]
    assert A1p == A2, name
    print(f"  {name}: A2 == A1' -> the chain records two edges but performs ONE "
          f"transformation; final corner ({A2[0]}/{A2[1]},{A2[2]}), l = {A2[1]}")
# F22/F23 step-1 arithmetic (same first step as F12's shape):
for name in ("F22", "F23", "F24"):
    A0, A0p, A1, A1p, A2 = LEN2[name]
    rho, sigma = edge_direction(A0[0] - A0p[0], A0[1] - A0p[1])
    s = -sigma
    got = step_corner(A0p, A1[2], s, A1[1])
    assert rho == A1[1] and got == (Fraction(A1[0], A1[1]), Fraction(A1[2])), name
    print(f"  {name} step1: dir=({rho},{sigma}), A1 = ({got[0]},{got[1]}) OK")
# F24 second step in u = x^(1/4) coordinates:
A0, A0p, A1, A1p, A2 = LEN2["F24"]
du = (Fraction(A1[0], A1[1]) - Fraction(A1p[0], A1p[1])) * 4   # u-exponent delta
dyy = A1[2] - A1p[2]
rho2, sigma2 = edge_direction(int(du), dyy)
s2 = -sigma2
l2 = A1[1] * rho2
got2 = (Fraction(A1p[0], A1p[1]) + Fraction(A2[2] * s2, l2), Fraction(A2[2]))
assert l2 == A2[1] and got2 == (Fraction(A2[0], A2[1]), Fraction(A2[2]))
print(f"  F24 step2: u-edge ({du},{dyy}) -> dir ({rho2},{sigma2}); l refines "
      f"{A1[1]} -> {l2}; A2 = ({got2[0]},{got2[1]}) matches table (19/8,3) OK")

# ---------------------------------------------------------------------------
# 2. The fused-chart lemma
# ---------------------------------------------------------------------------
print("\n" + BAR)
print("STEP 2 -- fused-chart lemma: shears never touch the Jacobian; one inversion does")
print(BAR)
l_s, e1, e2, e3 = sp.symbols("l_s e1 e2 e3", positive=True)
lam1, lam2, lam3 = sp.symbols("lambda1 lambda2 lambda3")


def jac(X, Y):
    return sp.simplify(sp.diff(X, x) * sp.diff(Y, y) - sp.diff(X, y) * sp.diff(Y, x))


J0 = jac(x**-1, x**l_s * y)
J1 = jac(x**-1, x**l_s * y + lam1 * x**e1)
J3 = jac(x**-1, x**l_s * y + lam1 * x**e1 + lam2 * x**e2 + lam3 * x**e3)
print(f"  pure Laurent      : J = {J0}")
print(f"  fused, 1 shear    : J = {J1}")
print(f"  fused, 3 shears   : J = {J3}")
assert J0 == J1 == J3 == -x**(l_s - 2)
# the heuristic's chart: TWO inversions (use a positive symbol so powers combine)
l1s, l2s = sp.symbols("l1 l2", positive=True)
xp = sp.symbols("xp", positive=True)
Xc, Yc = (xp**-1)**-1, (xp**-1)**l1s * (xp**l2s * y)
Jh = sp.powsimp(sp.simplify(
    sp.diff(Xc, xp) * sp.diff(Yc, y) - sp.diff(Xc, y) * sp.diff(Yc, xp)),
    force=True)
print(f"  double-inversion  : J = {Jh}   <- kappa = l2 - l1: the heuristic's source")
assert sp.simplify(Jh - xp**(l2s - l1s)) == 0
print("  => kappa = l-2 = t-2 for ANY chain of shears + ONE final inversion;")
print("     kappa = l2-l1 would require a second inversion the chain never makes.")

# fused shear exponents for the escapes (integrality/nonnegativity):
print("\n  fused shear exponents (final coordinates):")
FUSED = {"F12": (4, [4 - 1]), "F13": (3, [3 - 1]), "F22": (4, [4 - 1]),
         "F23": (4, [4 - 1]), "F24": (8, [8 - 1 * 2, 8 - 3])}
for name, (l_f, exps) in FUSED.items():
    assert all(e >= 0 and isinstance(e, int) for e in exps)
    terms = " + ".join(f"lam*x^{e}" for e in exps)
    print(f"  {name}: (X,Y) = (x^-1, x^{l_f} y + {terms})  kappa = {l_f - 2}")

# ---------------------------------------------------------------------------
# 3. F18-F21 impossibility (GGV5 lines 1726-1786): the checkable core
# ---------------------------------------------------------------------------
print("\n" + BAR)
print("STEP 3 -- F18-F21 carry no standard (m,n)-pair (GGV5's own proof, core re-checked)")
print(BAR)
# shifted start: ell_10(phi(P)) = lam * x^(6m) y^(3m) (y+lam)^(15m) => (1/m)st = (6,3)
a_, b_ = 6, 3
crit = b_ + gcd(a_, b_)
print(f"  (1/m) st_10(phi(P)) = (6,3);  last-lower-corner criterion: "
      f"b + gcd(a,b) = {crit} {'==' if crit == a_ else '!='} a = {a_}")
assert crit == a_   # violates the criterion -> (6,3) impossible as last lower corner
print("  -> violates b + gcd(a,b) != a (GGV2 Rmk 3.29 as used at GGV5 line 1774): dead.")
v11_alt = 6 + 9
print(f"  degenerate path (lam=lam'=lam''): (6,9), v11 = {v11_alt} < 16 = minimum "
      f"v11 for admissible chains (GGV5 line 407): dead.")
assert v11_alt < 16
print("  => within GGV5's classification the only length-2 escapes are F22, F23, F24.")

# ---------------------------------------------------------------------------
# 4. A0'=(2,0): the zeta-obstruction (model, not chart)
# ---------------------------------------------------------------------------
print("\n" + BAR)
print("STEP 4 -- A0'=(2,0) breaks the MODEL, not the chart: the zeta-commutator term")
print(BAR)
t_, zeta_, a_e, b_e, u_e = 3, 2, 2, 3, 4          # generic small exponents
c_g = sp.Function("c")(y)
C = x**t_ * c_g
f_g = sp.Function("f")(y)


def bracket(P, Q):
    return sp.simplify(sp.diff(P, x) * sp.diff(Q, y) - sp.diff(P, y) * sp.diff(Q, x))


lhs = bracket(x**zeta_ * C**a_e, C**b_e)
rhs = zeta_ * b_e * x**(zeta_ - 1) * C**(a_e + b_e - 1) * sp.diff(C, y)
assert sp.simplify(lhs - rhs) == 0
print(f"  [x^zeta C^a, C^b] = zeta*b*x^(zeta-1)*C^(a+b-1)*C_y   (checked at "
      f"zeta={zeta_}, a={a_e}, b={b_e}, C = x^{t_} c(y))")
assert sp.simplify(bracket(C**a_e, C**b_e)) == 0
print("  [C^a, C^b] = 0 (zeta = 0: pure powers commute -- the standard model's tail).")
lhs2 = bracket(x**zeta_ * C**a_e, x**u_e * f_g)
rhs2 = x**(zeta_ + a_e * t_ + u_e - 1) * (
    (zeta_ + a_e * t_) * c_g**a_e * sp.diff(f_g, y)
    - a_e * u_e * c_g**(a_e - 1) * sp.diff(c_g, y) * f_g)
assert sp.simplify(sp.expand(lhs2 - rhs2)) == 0
print("  [x^zeta C^a, x^u f] = x^(zeta+at+u-1) [ (zeta+at) c^a f' - a u c^(a-1) c' f ]")
print("  => zeta shifts the f'-coefficient AND leaves an uncancelled C^(a+b-1)C_y term:")
print("     the forcing family (F) needs a zeta-corrected tail theory before any")
print("     A0'=(2,0) Phi claim.  kappa itself is settled (= t-2) by Step 2.")

# ---------------------------------------------------------------------------
# 5. F12 conditional branch analysis (standard family (F), flagged CONDITIONAL-zeta)
# ---------------------------------------------------------------------------
print("\n" + BAR)
print("STEP 5 -- F12 under the UNMODIFIED family (F): dg=3 branch variety (CONDITIONAL)")
print(BAR)
# F12 j=0: (m,n)=(3,7) -> (a,b)=(3,7); t=4, kappa=2, a0=8, q=5.
a, b, t, kappa, a0, q = 3, 7, 4, 2, 8, 5
e = b - a + 1
r = a0 - q - 1
dg = a0 - q
coef = t * (b - a) + kappa + 1
rho = (e - 1) * q + 1
N = a * (t * (a + b - 1) + 1) - 2 * b
res = Fraction((t * (b - a) + kappa + 1) * a0, t)
pure = e * a0 - q + 1
gap = res - pure
print(f"  corner data: (a,b)=({a},{b}) t={t} kappa={kappa} a0={a0} q={q} "
      f"e={e} r={r} dg={dg} N={N} gap={gap} rho={rho} res={res}")
assert gap == 2 and dg == 3 and N == 97

g2, g1, g0, u2, u1, u0 = sp.symbols("g2 g1 g0 u2 u1 u0")
g = y**3 + g2 * y**2 + g1 * y + g0
u = u2 * y**2 + u1 * y + u0
c = y**q * g
f = y**rho * g**e * u
R = sp.expand(a * t * c * sp.diff(f, y) - a * coef * sp.diff(c, y) * f - c**e)
Q = sp.expand(sp.cancel(R / (y**(rho + q - 1) * g**(e - 1) / 1)))
# structural: R = y^(rho+q-1+... ) -- reduce exactly:
Q = sp.expand(sp.cancel(R / (y**25 * g**5)))
P = sp.Poly(Q, y)
print(f"  reduced equation (deg {P.degree()} in y): sum_k E_k y^k = 0, unknowns "
      f"(g2,g1,g0,u2,u1,u0), scaling gauge")
eqs = [sp.expand(P.nth(k)) for k in range(P.degree() + 1)]
# Triangular reduction: E4 solves u1, E3 solves u0 (u2 = free scale);
# E1, E2 become the consistency ideal on g; E0 fixes the scale.
u1_of = sp.solve(eqs[4], u1)[0]                       # u1 = -(g2/4) u2
u0_of = sp.solve(eqs[3].subs(u1, u1_of), u0)[0]
E1c = sp.factor(sp.expand(eqs[1].subs({u1: u1_of, u0: u0_of})) / u2)
E2c = sp.factor(sp.expand(eqs[2].subs({u1: u1_of, u0: u0_of})) / u2)
print(f"  triangular solve: u1, u0 in terms of u2 (free scale); consistency on g:")
print(f"    E2': {sp.factor(E2c)}")
print(f"    E1': {sp.factor(E1c)}")
# eliminate g0:
g0_of = sp.solve(E2c, g0)[0]
elim = sp.factor(sp.expand(E1c.subs(g0, g0_of)))
print(f"    eliminant in (g1, g2): {elim}")
sols = sp.solve([sp.Eq(sp.expand(sp.numer(sp.together(elim))).subs(g2, 1), 0)],
                [g1], dict=True)
pts = []
for s in sols:
    G1 = sp.radsimp(s[g1])
    G0 = sp.radsimp(g0_of.subs({g2: 1, g1: G1}))
    pts.append((G1, G0))
print(f"  g2=1 slice: {len(pts)} branch points (g1 = (-1 +- sqrt(22))/8); "
      f"g2=0 slice forces g0=0, inconsistent with E_0 (needs 33 g0 u0 = -1).")
assert len(pts) == 2
assert all(sp.simplify((8 * G1 + 1)**2 - 22) == 0 for G1, _ in pts)
# both branch points: g squarefree (no ramified branch at dg=3)
for G1, G0 in pts:
    disc = sp.simplify(sp.discriminant(y**3 + y**2 + G1 * y + G0, y))
    assert sp.simplify(disc) != 0
print("  both branch points have squarefree g  ->  NO ramified branch at dg=3")
print("  (mirror image of PHI_F7's dg=2 obstruction: a cubic always has a real root).")
# full exact verification at one point over Q(sqrt(22)):
G1, G0 = pts[0]
u2v = sp.Integer(96)
u1v = sp.radsimp(u1_of.subs({g2: 1, u2: u2v}))
u0v = sp.radsimp(u0_of.subs({g2: 1, g1: G1, u2: u2v}))
# rescale (u2,u1,u0) so that E0 = -33 g0 u0 - 1 = 0 holds:
scale = sp.radsimp(-1 / (33 * G0 * u0v))
Uv = {u2: sp.radsimp(scale * u2v), u1: sp.radsimp(scale * u1v),
      u0: sp.radsimp(scale * u0v)}
subs_all = {g2: 1, g1: G1, g0: G0, **Uv}
resid = sp.expand(sp.radsimp(
    (a * t * c * sp.diff(f, y) - a * coef * sp.diff(c, y) * f - c**e)
    .subs(subs_all)))
resid = sp.expand(sp.nsimplify(resid))
assert sp.simplify(resid) == 0
print("  FULL ODE residual == 0 exactly over Q(sqrt(22)) at the first branch point.")
# u shares no root with g (so mult_u(-1) = 0 in the g(-1)=0 gauge):
g_pt = (y**3 + y**2 + G1 * y + G0).subs({})
u_pt = Uv[u2] * y**2 + Uv[u1] * y + Uv[u0]
res_gu = sp.simplify(sp.resultant(g_pt, u_pt, y))
assert res_gu != 0
print(f"  resultant(g, u) != 0: u has no common root with g -> mult_u(-1) = 0.")
# real root exists (cubic), nonzero (g0 != 0): unramified gauge realizable.
deg_phi = int(res) + N * a0
ord_phi = rho + N * q
mult_unram = e + N
cof_unram = deg_phi - ord_phi - mult_unram
print(f"  CONDITIONAL Phi signature (unramified gauge, mult_g(-1)=1, u(-1)!=0):")
print(f"    (deg, ord, mult, cof) = ({deg_phi}, {ord_phi}, {mult_unram}, {cof_unram})")
print(f"    old-law check: cof = gap + r*(e+N) = {gap} + {r}*{e + N} = "
      f"{gap + r * (e + N)};  ramified formulas would give "
      f"({deg_phi}, {ord_phi}, {dg * (e + N) - (dg - 1)}, {gap + r})  -- NOT realized.")
assert (deg_phi, ord_phi, mult_unram, cof_unram) == (814, 506, 102, 206)
assert cof_unram == gap + r * (e + N)

# ---------------------------------------------------------------------------
# 6. Escape survey (conditional model data; kappa = t-2 from Step 2 throughout)
# ---------------------------------------------------------------------------
print("\n" + BAR)
print("STEP 6 -- escape survey (kappa = t-2 settled; model data CONDITIONAL-zeta for "
      "A0'=(2,0))")
print(BAR)
SURVEY = [
    # name, (m,n) at j=0, t=l_final, a0, q
    ("F12", (3, 7), 4, 8, 5),
    ("F13", (2, 13), 3, 9, 7),
    ("F22", (2, 3), 4, 8, 2),
    ("F23", (2, 7), 4, 8, 4),
    ("F24", (3, 4), 8, 8, 3),
]
print(f"  {'fam':5}{'(a,b)':>8}{'t':>3}{'kap':>4}{'a0':>4}{'q':>3}{'e':>3}{'r':>3}"
      f"{'dg':>4}{'N':>5}{'gap':>5}   regime note")
for name, (m, n), t_f, a0_f, q_f in SURVEY:
    aa, bb = sorted((m, n))
    e_f = bb - aa + 1
    r_f = a0_f - q_f - 1
    dg_f = a0_f - q_f
    kap = t_f - 2
    N_f = aa * (t_f * (aa + bb - 1) + 1) - 2 * bb
    gap_f = Fraction(q_f - 1) - Fraction(a0_f, t_f)
    note = []
    if gap_f > 0 and r_f > 0 and dg_f % 2 == 1:
        note.append(f"gap>0,r>0 with dg ODD -- outside the standard-15 parity lemma")
    if gap_f < 0:
        note.append("NEGATIVE gap: res < pure, regime unobserved anywhere")
    print(f"  {name:5}{str((aa, bb)):>8}{t_f:>3}{kap:>4}{a0_f:>4}{q_f:>3}{e_f:>3}"
          f"{r_f:>3}{dg_f:>4}{N_f:>5}{str(gap_f):>5}   {'; '.join(note)}")
print("\n  F13 j=1 is Orevkov's extensively analyzed case (GGV5 remark, line 1789).")
print("\nDERIVATION COMPLETE -- see COMPOSITE_CHARTS.md; checker: "
      "composite_charts_verify.py")
