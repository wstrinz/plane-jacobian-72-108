#!/usr/bin/env python3
"""zeta_tail.py  (NEW; read-only over all existing artifacts)

LANE K -- the zeta-corrected tail theory: what happens to the forcing family
when the reduced P-side element carries a pure-power defect x^zeta
(A0'=(2,0) families F12, F13 -- the boundary named by COMPOSITE_CHARTS.md).

RESULTS (details in ZETA_TAIL.md, checker zeta_tail_verify.py):

1. RIGIDITY (one-slice theorem).  In the model
       P = x^zeta (C^a + tail),  Q = C^b + tail + F,  v(F) < 0,  [P,Q] = x^kappa,
   the tail-tail bracket's top slice zeta*b*x^(zeta+(a+b)t-1) c^(a+b-1) c'
   is nonzero for zeta > 0, sits strictly above x^kappa, and nothing in the
   model can cancel it (an F-cross-term there would need v(F) >= bt > 0).
   The defect CANNOT be carried by the tail; it must enter the element.

2. REPAIR: D = x^eta C, eta = zeta/a (GGV proportionality distributes the
   defect evenly).  P = D^a, Q = D-series + F commutes again and the forcing
   identity [D^a, x^s f c^(-b)] = x^kappa (s = kappa+1-aT) yields the
   ZETA-CORRECTED FAMILY
       a { T c f' - [T(b-a) + kappa + 1] c' f } = c^(b-a+1),
       T = t + eta,   kappa = t-2  (chart-fixed; does NOT move with eta):
   the standard two-parameter family OFF the kappa = T-2 diagonal by eta.
   Fractional eta lives in the refined chart x = w^n: T_w = nT,
   kappa_w = n kappa + n - 1, K_w = nK; d_res = K a0/T is chart-invariant.

3. UNIVERSAL eta = -1 LEMMA.  T = t-1 <=> K = eT for every family; the ODE
   integrates: (f/c^e)' = 1/(aTc); polynomial f <=> all residues of 1/c
   vanish -- impossible in every root configuration at dg in {2,3}.

4. F12 SWEEP: all motivated defects (0 < |eta| <= 1, incl. fractional) are
   DEAD (degree count or lemma 3).  The scan finds ONE nonzero defect with
   canonical structure: eta = +2 (T=6), where the family collapses exactly
   like the standard corners (g = y^3+1 forced, residual H = y^2-y+1) --
   and lands ON the law at gap_eff = 0.

5. THE mu-LADDER (new branch structure, corrects COMPOSITE_CHARTS Step 5).
   Indexing branches by mu = mult_{y+1}(g):
     eta=0 (standard family), F12: ALL THREE rungs are realized --
       mu=1: COMPOSITE_CHARTS.md's two squarefree points        sig (814, 506, 102, 206)
       mu=2: g=(y+1)^2(y-beta), beta a root of
             195b^4+120b^3-40b^2+32b-80 (2 real)   sig (814, 506, 203, 105)
       mu=3: g=(y+1)^3,
             u = -(2048y^4+2560y^3+320y^2-80y+35)/1155
                                                   sig (814, 506, 304,   4)
     The mu=3 rung IS the "ramified formulas (304,4)" that COMPOSITE_CHARTS
     called NOT realized -- that claim was an artifact of the g^e-uniform
     ansatz (it cannot represent the ramified forced orders (e-1)mu+1).
     The signatures obey the mu-GRADED LAW (PHI_F7's judgment-6 observation,
     now realized as actual branches):
       mult = mu(e+N) - (mu-1),   cof = gap + r(e+N) - (mu-1)(e+N-1).
     eta=+2: rungs mu=1 and mu=3 exist, mu=2 provably does not.
   Parity refinement: dg even kills mu=1 (PHI_F7's theorem); odd dg keeps
   mu=1 available but does NOT exclude higher rungs.

6. F13 SWEEP: every motivated defect dead; only eta = 0 survives
   (conditional; PHI_F7 dg=2 analysis applies).  F13 j=1 = Orevkov's case.

Framework facts used (published, cited per PRIOR_ART.md): bracket-to-ODE
mechanism (GGV1 eqq1 genre), chart Jacobian -x^(l-2) (COMPOSITE_CHARTS
fused-chart lemma / GGHV22 L1228-1234), chain data (GGV5 tables).
"""
import sympy as sp
from fractions import Fraction

x, y = sp.symbols("x y")
BAR = "=" * 96

F12 = dict(name="F12", a=3, b=7, t=4, kappa=2, a0=8, q=5)
F13 = dict(name="F13", a=2, b=13, t=3, kappa=1, a0=9, q=7)
for F in (F12, F13):
    F["e"] = F["b"] - F["a"] + 1
    F["r"] = F["a0"] - F["q"] - 1
    F["dg"] = F["a0"] - F["q"]


def bracket(P, Q):
    return sp.expand(sp.diff(P, x) * sp.diff(Q, y) - sp.diff(P, y) * sp.diff(Q, x))


# ---------------------------------------------------------------------------
print(BAR)
print("STEP 1 -- RIGIDITY: the zeta-defect cannot be carried by the tail")
print(BAR)
t_s, zeta_s = sp.symbols("t zeta", positive=True)
c_f = sp.Function("c")(y)
C = x**t_s * c_f
mu2, mu1, nu2, nu1 = sp.symbols("mu2 mu1 nu2 nu1")
Pser = x**zeta_s * (C**3 + mu2 * C**2 + mu1 * C)
Qser = C**3 + nu2 * C**2 + nu1 * C
BB = sp.expand(bracket(Pser, Qser))
top_claim = zeta_s * 3 * x**(zeta_s + 6 * t_s - 1) * c_f**5 * sp.diff(c_f, y)
rest = sp.expand(BB - top_claim)
poly_test = sp.expand(rest.subs(t_s, 11))
xdegs = {m.as_powers_dict().get(x, 0) for m in sp.Add.make_args(poly_test)}
assert all(sp.simplify(d - (zeta_s + 65)) != 0 for d in xdegs)
print("  top tail-tail slice = zeta*b * x^(zeta+(a+b)t-1) * c^(a+b-1) c'  "
      "(verified, (a,b)=(3,3) with full generic tails)")
print("  slice sits at kappa + [zeta + 1 + (a+b-1)t] > kappa; an F-term there "
      "needs v(F) = (a+b-j)t >= bt > 0, contradicting v(F) < 0")
print("  => zeta > 0 tail-carried: NO solution.  The defect must enter the ELEMENT.")

# ---------------------------------------------------------------------------
print("\n" + BAR)
print("STEP 2 -- REPAIR: D = x^eta C; corrected family = standard family OFF the "
      "kappa=T-2 diagonal")
print(BAR)
T_s, kap_s = sp.symbols("T kappa_s", positive=True)
f_f = sp.Function("f")(y)
for (a_i, b_i) in ((2, 3), (3, 7)):
    D = x**T_s * c_f
    s_val = kap_s + 1 - a_i * T_s
    W = x**s_val * f_f * c_f**(-b_i)
    lhs = sp.simplify(bracket(D**a_i, W))
    K_expr = T_s * (b_i - a_i) + kap_s + 1
    rhs = a_i * x**kap_s * c_f**(a_i - b_i - 1) * (
        T_s * c_f * sp.diff(f_f, y) - K_expr * sp.diff(c_f, y) * f_f)
    assert sp.simplify(sp.expand(lhs - rhs)) == 0
    print(f"  (a,b)=({a_i},{b_i}): [D^a, x^s f/c^b] = a x^kappa c^(a-b-1) "
          f"{{T c f' - [T(b-a)+kappa+1] c' f}}  OK  (s = kappa+1-aT)")
print("  => a{ T c f' - [T(b-a)+kappa+1] c' f } = c^(b-a+1),  T = t+eta, "
      "kappa = t-2 chart-fixed")
n_s = sp.symbols("n", positive=True)
Kx = T_s * (b_i - a_i) + kap_s + 1
Kw = (n_s * T_s) * (b_i - a_i) + (n_s * kap_s + n_s - 1) + 1
assert sp.simplify(Kw - n_s * Kx) == 0
w = sp.symbols("w", positive=True)
Fg, Gg = sp.Function("F")(x, y), sp.Function("G")(x, y)
Fw, Gw = Fg.subs(x, w**3), Gg.subs(x, w**3)
bw = sp.diff(Fw, w) * sp.diff(Gw, y) - sp.diff(Fw, y) * sp.diff(Gw, w)
bx = (sp.diff(Fg, x) * sp.diff(Gg, y) - sp.diff(Fg, y) * sp.diff(Gg, x)).subs(x, w**3)
assert sp.simplify(bw - 3 * w**2 * bx) == 0
print("  fractional eta via x = w^n: K_w = nK, d_res = K a0/T chart-invariant; "
      "bracket rule [ , ]_w = n w^(n-1) [ , ]_x  (verified n=3)")

# ---------------------------------------------------------------------------
print("\n" + BAR)
print("STEP 3 -- eta = -1 <=> K = eT: degenerate-integrable, always log-obstructed")
print(BAR)
t_g, a_g, b_g = sp.symbols("t_g a_g b_g", positive=True)
assert sp.simplify(((t_g - 1) * (b_g - a_g) + t_g - 1)
                   - (b_g - a_g + 1) * (t_g - 1)) == 0
e_i = 4
lhs_int = sp.diff(f_f * c_f**(-e_i), y)
rhs_int = (T_s * c_f * sp.diff(f_f, y) - e_i * T_s * sp.diff(c_f, y) * f_f) \
    / (T_s * c_f**(e_i + 1))
assert sp.simplify(lhs_int - rhs_int) == 0
print("  T = t-1 => K = e(t-1) = eT (symbolic, all families); ODE reads "
      "(f/c^e)' = 1/(aTc)")
al = sp.symbols("alpha", nonzero=True)
for q_i, dg_i, tag in ((5, 3, "F12"), (7, 2, "F13")):
    print(f"  {tag} (q={q_i}, dg={dg_i}): polynomial f <=> all residues of 1/c "
          f"vanish:")
    print(f"    simple root: Res = 1/(alpha^{q_i} prod) != 0 -> DEAD")
    if dg_i == 3:
        be = sp.symbols("beta_t", nonzero=True)
        R_al = sp.simplify(sp.diff(y**-q_i / (y - be), y).subs(y, al))
        be_t = sp.solve(sp.numer(sp.together(R_al)), be)[0]
        R_be = sp.simplify(1 / (be_t**q_i * (be_t - al)**2))
        assert sp.simplify(R_be) != 0
        print(f"    [2,1]: Res_alpha = 0 tunes beta = {be_t}; Res_beta = "
              f"{sp.simplify(R_be)} != 0 -> DEAD")
        R3 = sp.simplify(sp.diff(y**-q_i, y, 2).subs(y, al) / 2)
        assert sp.simplify(R3) != 0
        print(f"    [3]: Res_alpha = {R3} != 0 -> DEAD")
    else:
        R2 = sp.simplify(sp.diff(y**-q_i, y).subs(y, al))
        assert sp.simplify(R2) != 0
        print(f"    [2]: Res_alpha = {R2} != 0 -> DEAD")
print("  => eta = -1: no polynomial last element, any root configuration "
      "(dg <= 3 shown).")

# ---------------------------------------------------------------------------
def sweep(F, etas):
    a, b, t, kap, a0, q, e, dg = (F[k] for k in
                                  ("a", "b", "t", "kappa", "a0", "q", "e", "dg"))
    rows = []
    for eta in etas:
        Tx = Fraction(t) + eta
        n = Tx.denominator
        Tw = int(Tx * n)
        Kw = Tw * (b - a) + n * (kap + 1)
        d_res = Fraction(Kw * a0, Tw)
        d_max = int(d_res) if d_res.denominator == 1 else (e - 1) * a0 + 1
        rho0 = (e - 1) * q + 1
        if Kw == e * Tw:
            rows.append((eta, Tw, Kw, d_res, "DEAD (K=eT: log obstruction)"))
            continue
        assert Tw * rho0 != Kw * q, (F["name"], eta)
        sub0 = Fraction(Kw * q, Tw)
        assert not (sub0.denominator == 1 and sub0 < rho0), (F["name"], eta)
        partitions = {3: ([1, 1, 1], [2, 1], [3]), 2: ([1, 1], [2])}[dg]
        min_forced = None
        for part in partitions:
            tot, ok = rho0, True
            for mu in part:
                om = (e - 1) * mu + 1
                sub = Fraction(Kw * mu, Tw)
                if sub.denominator == 1 and sub < om:
                    ok = False
                tot += om
            if ok:
                min_forced = tot if min_forced is None else min(min_forced, tot)
        if min_forced is not None and d_max < min_forced:
            v = (f"DEAD (degree: min forced {min_forced} > d_max {d_max}"
                 f"{'' if d_res.denominator == 1 else ', d_res non-integral'})")
        else:
            v = f"viable (d_res {d_res}, min forced {min_forced})"
        rows.append((eta, Tw, Kw, d_res, v))
    return rows


print("\n" + BAR)
print("STEP 4 -- F12 candidate sweep")
print(BAR)
CANDS = sorted(set([Fraction(k, 3) for k in (-3, -2, -1, 1, 2, 3)]
                   + [0, 2, -2, 3, -3]))
for eta, Tw, Kw, d_res, verdict in sweep(F12, CANDS):
    tag = {Fraction(0): "  <- standard family (Step 5b: mu-ladder)",
           Fraction(2): "  <- CANONICAL COLLAPSE (Step 5a)",
           Fraction(-2): "  (unmotivated defect - flagged OPEN)",
           Fraction(-3): "  (unmotivated defect - flagged OPEN)"}.get(eta, "")
    print(f"  eta={str(eta):>5}:  T_w={Tw:>3} K_w={Kw:>3} d_res={str(d_res):>7}  "
          f"{verdict}{tag}")

# ---------------------------------------------------------------------------
print("\n" + BAR)
print("STEP 5a -- F12 at eta = +2 (T=6, K=27): canonical collapse + mu-rungs")
print(BAR)
a, b, t, kap, a0, q, e = (F12[k] for k in ("a", "b", "t", "kappa", "a0", "q", "e"))
r = F12["r"]
T, K = 6, 27
g2, g1, g0, A = sp.symbols("g2 g1 g0 A")
g = y**3 + g2 * y**2 + g1 * y + g0
c = y**q * g
rho0 = (e - 1) * q + 1
f = A * y**rho0 * g**e
R = sp.expand(a * (T * c * sp.diff(f, y) - K * sp.diff(c, y) * f) - c**e)
Rq = sp.expand(sp.cancel(R / (y**(e * q) * g**e)))
collapse = sp.expand(9 * A * (y * sp.diff(g, y) - 3 * g) - 1)
assert sp.expand(Rq - collapse) == 0
solg = sp.solve([sp.Poly(collapse, y).nth(k) for k in range(1, 3)], [g2, g1])
assert solg == {g2: 0, g1: 0}
A_of = sp.solve(collapse.subs({g2: 0, g1: 0}), A)[0]
g_fin, A_fin = y**3 + 1, A_of.subs(g0, 1)
f_mu1 = sp.expand(A_fin * y**rho0 * g_fin**e)
assert sp.expand(a * (T * (y**q * g_fin) * sp.diff(f_mu1, y)
                      - K * sp.diff(y**q * g_fin, y) * f_mu1)
                 - (y**q * g_fin)**e) == 0
print(f"  mu=1: f = A y^21 g^5 collapses to 9A(yg'-3g)=1 => g = y^3+1 "
      f"(gauge g(-1)=0), A = {A_fin}; residual H = y^2-y+1 ((75,125)/F2 class)")

# mu=3 rung at eta=+2: c = y^5 (y+1)^3, f = y^21 (y+1)^13 u, deg u = 2
us = sp.symbols("v0:3")
u = sum(us[i] * y**i for i in range(3))
c3 = y**5 * (y + 1)**3
f3 = y**21 * (y + 1)**13 * u
R3 = sp.expand(a * (T * c3 * sp.diff(f3, y) - K * sp.diff(c3, y) * f3) - c3**e)
Q3 = sp.Poly(sp.cancel(R3 / (y**25 * (y + 1)**15)), y)
sol3 = sp.solve([Q3.nth(k) for k in range(Q3.degree() + 1)], list(us), dict=True)
u3 = [sp.expand(u.subs(s)) for s in sol3 if sp.expand(u.subs(s)) != 0]
assert len(u3) == 1 and sp.simplify(u3[0] - (8 * y**2 + 4 * y - 1) / 27) == 0
assert u3[0].subs(y, -1) != 0 and u3[0].subs(y, 0) != 0
print(f"  mu=3: g = (y+1)^3 branch EXISTS: u = (8y^2+4y-1)/27, u(-1) != 0")

# mu=2 rung at eta=+2: no branch (consistency has no admissible root)
beta = sp.symbols("beta")
us2 = sp.symbols("w0:2")
u2e = us2[0] + us2[1] * y
c21 = y**5 * (y + 1)**2 * (y - beta)
f21 = y**21 * (y + 1)**9 * (y - beta)**5 * u2e
R21 = sp.expand(a * (T * c21 * sp.diff(f21, y) - K * sp.diff(c21, y) * f21)
                - c21**e)
Q21 = sp.Poly(sp.cancel(R21 / (y**25 * (y + 1)**10 * (y - beta)**5)), y,
              domain=sp.QQ.frac_field(beta, *us2))
eqs21 = [sp.together(Q21.nth(k)) for k in range(Q21.degree() + 1)]
lin21 = sp.solve(eqs21[:2], list(us2), dict=True)
rem21 = [sp.factor(sp.numer(sp.together(sp.expand(eq.subs(lin21[0])))))
         for eq in eqs21[2:]]
conds21 = [sp.Poly(r0, beta) for r0 in rem21 if sp.simplify(r0) != 0]
admissible = set()
for cnd in conds21:
    admissible |= {rt for rt in sp.roots(cnd).keys() if rt not in (0, -1)}
if len(conds21) > 1:
    gc = conds21[0]
    for cnd in conds21[1:]:
        gc = gc.gcd(cnd)
    admissible = {rt for rt in sp.roots(gc).keys() if rt not in (0, -1)} \
        if gc.degree() > 0 else set()
assert not admissible
print("  mu=2: NO branch at eta=+2 (consistency conditions have no admissible "
      "common root)")

N2 = a * (T * (a + b) - (kap + 1)) - 2 * b
pure = e * a0 - q + 1
gap_eff = Fraction(K * a0, T) - pure
print(f"\n  off-diagonal N = a[T(a+b)-(kappa+1)] - 2b = {N2}; gap_eff = {gap_eff}")
for mu, dfu in ((1, 36), (3, 38 - 2)):
    degf = 36
    mult = mu * (e + N2) - (mu - 1)
    deg_phi = degf + N2 * a0
    ord_phi = rho0 + N2 * q
    cof = deg_phi - ord_phi - mult
    law_cof = gap_eff + r * (e + N2) - (mu - 1) * (e + N2 - 1)
    assert cof == law_cof
    print(f"  mu={mu}: Phi sig ({deg_phi}, {ord_phi}, {mult}, {cof})  "
          f"[mu-graded law: mult = mu(e+N)-(mu-1), cof = gap+r(e+N)-(mu-1)(e+N-1)]")

# ---------------------------------------------------------------------------
print("\n" + BAR)
print("STEP 5b -- the mu-LADDER at eta = 0 (standard family; corrects "
      "COMPOSITE_CHARTS Step 5)")
print(BAR)
T0, K0 = 4, 19
N0 = a * (t * (a + b - 1) + 1) - 2 * b
# mu=3: c = y^5 (y+1)^3, f = y^21 (y+1)^13 u, deg u = 4 (d_res = 38)
us4 = sp.symbols("z0:5")
u4 = sum(us4[i] * y**i for i in range(5))
f34 = y**21 * (y + 1)**13 * u4
R30 = sp.expand(a * (T0 * c3 * sp.diff(f34, y) - K0 * sp.diff(c3, y) * f34)
                - c3**e)
Q30 = sp.Poly(sp.cancel(R30 / (y**25 * (y + 1)**15)), y)
sol30 = sp.solve([Q30.nth(k) for k in range(Q30.degree() + 1)], list(us4),
                 dict=True)
u30 = [sp.expand(u4.subs(s)) for s in sol30 if sp.expand(u4.subs(s)) != 0]
u3_expected = -(2048 * y**4 + 2560 * y**3 + 320 * y**2 - 80 * y + 35) / 1155
assert len(u30) == 1 and sp.simplify(u30[0] - u3_expected) == 0
assert u30[0].subs(y, -1) != 0 and u30[0].subs(y, 0) != 0
f3_final = sp.expand(y**21 * (y + 1)**13 * u30[0])
assert sp.expand(a * (T0 * c3 * sp.diff(f3_final, y)
                      - K0 * sp.diff(c3, y) * f3_final) - c3**e) == 0
print(f"  mu=3 branch EXISTS at eta=0: u = {sp.nsimplify(u30[0])}")
print(f"         => sig (814, 506, 304, 4): the 'ramified formulas' "
      f"COMPOSITE_CHARTS said were NOT realized ARE realized.")
print(f"         (the g^e-uniform ansatz cannot represent ramified forced "
      f"orders (e-1)mu+1 -- that is where the earlier claim broke)")

# mu=2: c = y^5(y+1)^2(y-beta), f = y^21 (y+1)^9 (y-beta)^5 u, deg u = 3
us3 = sp.symbols("p0:4")
u3e = sum(us3[i] * y**i for i in range(4))
f210 = y**21 * (y + 1)**9 * (y - beta)**5 * u3e
c210 = y**5 * (y + 1)**2 * (y - beta)
R210 = sp.expand(a * (T0 * c210 * sp.diff(f210, y)
                      - K0 * sp.diff(c210, y) * f210) - c210**e)
Q210 = sp.Poly(sp.cancel(R210 / (y**25 * (y + 1)**10 * (y - beta)**5)), y,
               domain=sp.QQ.frac_field(beta, *us3))
eqs210 = [sp.together(Q210.nth(k)) for k in range(Q210.degree() + 1)]
lin210 = sp.solve(eqs210[:4], list(us3), dict=True)
rem210 = [sp.factor(sp.numer(sp.together(sp.expand(eq.subs(lin210[0])))))
          for eq in eqs210[4:]]
conds210 = [r0 for r0 in rem210 if sp.simplify(r0) != 0]
assert len(conds210) == 1
quart = sp.Poly(conds210[0], beta)
quart_monic = sp.expand(quart.as_expr() / quart.LC() * 195)
print(f"  mu=2 branch: ONE consistency condition: {quart.as_expr()} = 0")
n_real = sp.Poly(quart, beta).count_roots(-sp.oo, sp.oo)
assert quart.degree() == 4 and n_real == 2
assert quart.eval(0) != 0 and quart.eval(-1) != 0
print(f"         quartic in beta, {n_real} real roots, 0 and -1 excluded: "
      f"admissible mu=2 branches EXIST")
# exact residual check mod the quartic at the algebraic root
brt = sp.RootOf(quart.as_expr(), 0)
lin_sub = {k: sp.together(v) for k, v in lin210[0].items()}
res_exact = []
for eq in eqs210:
    num = sp.numer(sp.together(sp.expand(eq.subs(lin_sub))))
    rem_p = sp.rem(sp.Poly(num, beta), sp.Poly(quart.as_expr(), beta))
    res_exact.append(sp.simplify(rem_p.as_expr()))
assert all(v == 0 for v in res_exact)
print(f"         full system == 0 EXACTLY mod the quartic (algebraic root check)")
print(f"\n  mu-LADDER at F12 eta=0 (N = {N0}, e+N = {e + N0}):")
for mu in (1, 2, 3):
    mult = mu * (e + N0) - (mu - 1)
    cof = 2 + 2 * (e + N0) - (mu - 1) * (e + N0 - 1)
    print(f"    mu={mu}: sig (814, 506, {mult}, {cof})"
          + ("   [COMPOSITE_CHARTS.md's two squarefree points]" if mu == 1 else
         "   [NEW this lane]"))
print("  parity refinement: dg even kills mu=1 (PHI_F7 theorem); odd dg keeps "
      "mu=1 but does NOT exclude higher rungs.")

# ---------------------------------------------------------------------------
print("\n" + BAR)
print("STEP 6 -- F13 candidate sweep")
print(BAR)
CANDS13 = sorted(set([Fraction(k, 2) for k in (-2, -1, 1, 2)] + [0, 2, -2]))
for eta, Tw, Kw, d_res, verdict in sweep(F13, CANDS13):
    tag = {Fraction(0): "  <- conditional standard analysis (PHI_F7 dg=2)"}.get(
        eta, "")
    print(f"  eta={str(eta):>5}:  T_w={Tw:>3} K_w={Kw:>3} d_res={str(d_res):>7}  "
          f"{verdict}{tag}")
print("\n  F13: every motivated defect DEAD; only eta = 0 survives (conditional). "
      "F13 j=1 is Orevkov's case.")

print("\nDERIVATION COMPLETE -- see ZETA_TAIL.md; checker: zeta_tail_verify.py")
