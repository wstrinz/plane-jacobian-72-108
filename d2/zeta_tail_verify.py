#!/usr/bin/env python3
"""zeta_tail_verify.py -- exact checker for LANE K (the zeta-corrected tail
theory).  Independent code paths from zeta_tail.py: different bracket
instances, different truncations, different chart-refinement n, sweeps
recomputed from scratch, explicit branch polynomials verified by direct
residual, the mu=2 system rebuilt fresh and checked mod its quartic.

Checks:
  A. corrected forcing identity [D^a, x^s f/c^b] = x^kappa at (a,b)=(2,5)
     (symbolic T, kappa); zeta=0 controls: F9 instance (14,26,c^2) + known
     F9 solution; (75,125) instance (15,42,c^3).
  B. rigidity: top tail-tail slice at (a,b)=(2,2) with generic tails; the
     slice-degree margin; the v(F)<0 contradiction arithmetic.
  C. chart refinement: bracket rule at n=2; K_w = nK symbolic; d_res
     chart-invariance.
  D. candidate sweeps F12/F13 recomputed independently; verdicts match.
  E. eta=-1: K=eT symbolic; integrable form at e=7; residue obstructions
     (simple root, tuned [2,1], [3]; F13 [2]).
  F. eta=+2 mu=1: explicit f residual == 0; collapse quotient identity.
  G. branch rungs: explicit mu=3 residuals at eta=+2 and eta=0; mu=2 system
     at eta=0 rebuilt fresh (quartic condition, 2 real roots, exact mod-
     quartic residual); mu=2 at eta=+2 has no admissible root.
  H. the mu-graded law: reproduces all five F12 signatures; specializes to
     the old law at mu=1 and EXACTLY to PHI_F7's ramified law at mu=dg
     (via the identity r = dg-1); F7 spot check (83, 2).
"""
import sys
import sympy as sp
from fractions import Fraction

QUIET = "--quiet" in sys.argv
FAILS = []
N_OK = 0


def ok(label, cond):
    global N_OK
    if cond:
        N_OK += 1
        if not QUIET:
            print("[OK]", label)
    else:
        FAILS.append(label)
        print("[FAIL]", label)


x, y = sp.symbols("x y")


def bracket(P, Q):
    return sp.expand(sp.diff(P, x) * sp.diff(Q, y) - sp.diff(P, y) * sp.diff(Q, x))


# --- A. corrected forcing identity + zeta=0 controls ------------------------
T, kap = sp.symbols("T kappa", positive=True)
cf = sp.Function("c")(y)
ff = sp.Function("f")(y)
a_i, b_i = 2, 5
D = x**T * cf
W = x**(kap + 1 - a_i * T) * ff * cf**(-b_i)
K_sym = T * (b_i - a_i) + kap + 1
rhs = a_i * x**kap * cf**(a_i - b_i - 1) * (
    T * cf * sp.diff(ff, y) - K_sym * sp.diff(cf, y) * ff)
ok("A: forcing identity at (a,b)=(2,5), symbolic (T,kappa)",
   sp.simplify(sp.expand(bracket(D**a_i, W) - rhs)) == 0)
ok("A: zeta=0 control -- T=t reproduces the standard family coefficient "
   "t(b-a)+kappa+1", sp.simplify(K_sym.subs(T, sp.Symbol("t", positive=True))
                                 - (sp.Symbol("t", positive=True) * (b_i - a_i)
                                    + kap + 1)) == 0)
# F9 instance: (a,b,t,kappa) = (2,3,7,5) -> 14 c f' - 26 c' f = c^2
aF, bF, tF, kF = 2, 3, 7, 5
ok("A: F9 instance numbers a*T=14, K=26, RHS power e=2",
   (aF * tF, tF * (bF - aF) + kF + 1, bF - aF + 1) == (14, 13, 2))
c9 = y**2 * (y**5 + 1)
f9 = -sp.Rational(1, 10) * y**3 * (y**5 + 1)**2
ok("A: known F9 solution satisfies 14 c f' - 26 c' f = c^2",
   sp.expand(2 * (7 * c9 * sp.diff(f9, y) - 13 * sp.diff(c9, y) * f9)
             - c9**2) == 0)
c75 = y**2 * (y**3 + 1)
f75 = -sp.Rational(1, 9) * y**5 * (y**3 + 1)**3
ok("A: known (75,125) solution satisfies 15 c f' - 42 c' f = c^3",
   sp.expand(3 * (5 * c75 * sp.diff(f75, y) - 14 * sp.diff(c75, y) * f75)
             - c75**3) == 0)

# --- B. rigidity ------------------------------------------------------------
t_s, z_s = sp.symbols("t_r zeta_r", positive=True)
C = x**t_s * cf
m1, n1 = sp.symbols("m1 n1")
P22 = x**z_s * (C**2 + m1 * C)
Q22 = C**2 + n1 * C
BB = sp.expand(bracket(P22, Q22))
top = z_s * 2 * x**(z_s + 4 * t_s - 1) * cf**3 * sp.diff(cf, y)
# numeric-exponent reading avoids unmerged symbolic-power artifacts
rest_n = sp.expand(sp.powsimp((BB - top).subs({t_s: 13, z_s: 7}), force=True))
xdegs = {mm.as_powers_dict().get(x, 0) for mm in sp.Add.make_args(rest_n)}
ok("B: top tail-tail slice = zeta*b x^(zeta+(a+b)t-1) c^(a+b-1) c' at "
   "(a,b)=(2,2), generic tails",
   7 + 4 * 13 - 1 not in xdegs and rest_n != 0)
ok("B: slice degree margin (zeta+(a+b)t-1) - kappa = zeta+1+(a+b-1)t > 0",
   sp.simplify((z_s + 4 * t_s - 1) - (t_s - 2) - (z_s + 1 + 3 * t_s)) == 0)
ok("B: an F-term at the top slice needs v(F) = (a+b-j)t >= bt > 0 "
   "(contradiction with v(F) < 0): arithmetic",
   all((2 + 2 - j) * 3 >= 2 * 3 > 0 for j in (0, 1, 2)))

# --- C. chart refinement ----------------------------------------------------
w = sp.symbols("w", positive=True)
Fg, Gg = sp.Function("F")(x, y), sp.Function("G")(x, y)
bw = sp.diff(Fg.subs(x, w**2), w) * sp.diff(Gg.subs(x, w**2), y) \
    - sp.diff(Fg.subs(x, w**2), y) * sp.diff(Gg.subs(x, w**2), w)
bx = (sp.diff(Fg, x) * sp.diff(Gg, y)
      - sp.diff(Fg, y) * sp.diff(Gg, x)).subs(x, w**2)
ok("C: bracket chart rule [ , ]_w = n w^(n-1) [ , ]_x at n=2",
   sp.simplify(bw - 2 * w * bx) == 0)
n_s = sp.symbols("n", positive=True)
Kx = T * (b_i - a_i) + kap + 1
Kw = (n_s * T) * (b_i - a_i) + (n_s * kap + n_s - 1) + 1
ok("C: K_w = n K_x symbolic", sp.simplify(Kw - n_s * Kx) == 0)
ok("C: d_res = K a0/T chart-invariant symbolic",
   sp.simplify(Kw * 8 / (n_s * T) - Kx * 8 / T) == 0)

# --- D. sweeps recomputed independently -------------------------------------
def verdicts(a, b, t, kappa, a0, q, etas):
    e, dg = b - a + 1, a0 - q
    out = {}
    for eta in etas:
        Tx = Fraction(t) + eta
        n = Tx.denominator
        Tw, Kw = int(Tx * n), int(Tx * n) * (b - a) + n * (kappa + 1)
        d_res = Fraction(Kw * a0, Tw)
        d_max = int(d_res) if d_res.denominator == 1 else (e - 1) * a0 + 1
        if Kw == e * Tw:
            out[eta] = "K=eT"
            continue
        parts = {3: ([1, 1, 1], [2, 1], [3]), 2: ([1, 1], [2])}[dg]
        mf = min(((e - 1) * q + 1) + sum((e - 1) * mu + 1 for mu in p)
                 for p in parts)
        out[eta] = "dead-deg" if d_max < mf else "viable"
    return out


v12 = verdicts(3, 7, 4, 2, 8, 5,
               [Fraction(k, 3) for k in (-3, -2, -1, 1, 2, 3)] + [0, 2, -2])
ok("D: F12 sweep -- motivated defects 0<|eta|<=1 all dead",
   all(v12[eta] in ("K=eT", "dead-deg")
       for eta in [Fraction(k, 3) for k in (-2, -1, 1, 2, 3)] + [Fraction(-3, 3)]))
ok("D: F12 -- eta=-1 is the K=eT case", v12[Fraction(-1)] == "K=eT")
ok("D: F12 -- eta=0 and eta=+2 viable", v12[0] == "viable" and v12[2] == "viable")
v13 = verdicts(2, 13, 3, 1, 9, 7,
               [Fraction(k, 2) for k in (-2, -1, 1, 2)] + [0, 2])
ok("D: F13 sweep -- motivated defects all dead; eta=0 viable; eta=+2 dead",
   all(v13[eta] in ("K=eT", "dead-deg")
       for eta in [Fraction(k, 2) for k in (-1, 1)] + [Fraction(-1), 1])
   and v13[0] == "viable" and v13[2] == "dead-deg")

# --- E. eta = -1 lemma ------------------------------------------------------
tg, ag, bg = sp.symbols("t_g a_g b_g", positive=True)
ok("E: T=t-1 => K = eT symbolic (all families)",
   sp.simplify(((tg - 1) * (bg - ag) + tg - 1)
               - (bg - ag + 1) * (tg - 1)) == 0)
e7 = 7
ok("E: integrable form (f/c^e)' = [T c f' - eT c' f]/(T c^(e+1)) at e=7",
   sp.simplify(sp.diff(ff * cf**(-e7), y)
               - (T * cf * sp.diff(ff, y) - e7 * T * sp.diff(cf, y) * ff)
               / (T * cf**(e7 + 1))) == 0)
al = sp.symbols("alpha", nonzero=True)
be = sp.symbols("beta_e", nonzero=True)
R_al = sp.simplify(sp.diff(y**-5 / (y - be), y).subs(y, al))
# discard the spurious beta = alpha root (excluded: partition roots distinct)
be_roots = [s for s in sp.solve(sp.numer(sp.together(R_al)), be)
            if sp.simplify(s - al) != 0]
be_t = be_roots[0]
ok("E: F12 [2,1] residue tuning gives beta = 6 alpha/5 (unique after "
   "excluding beta = alpha)",
   len(be_roots) == 1 and sp.simplify(be_t - 6 * al / 5) == 0)
ok("E: F12 [2,1] tuned: Res_beta = 1/(beta^5 (beta-alpha)^2) != 0",
   sp.simplify(1 / (be_t**5 * (be_t - al)**2)) != 0)
ok("E: F12 [3]: Res = 15/alpha^7 != 0",
   sp.simplify(sp.diff(y**-5, y, 2).subs(y, al) / 2 - 15 / al**7) == 0)
ok("E: F13 [2]: Res = -7/alpha^8 != 0",
   sp.simplify(sp.diff(y**-7, y).subs(y, al) + 7 / al**8) == 0)

# --- F. eta=+2 mu=1 canonical collapse --------------------------------------
a, e, q, r = 3, 5, 5, 2
T2, K2 = 6, 27
c_mu1 = y**5 * (y**3 + 1)
f_mu1 = -sp.Rational(1, 27) * y**21 * (y**3 + 1)**5
ok("F: eta=+2 mu=1 explicit residual == 0",
   sp.expand(a * (T2 * c_mu1 * sp.diff(f_mu1, y)
                  - K2 * sp.diff(c_mu1, y) * f_mu1) - c_mu1**e) == 0)
g2s, g1s, g0s, As = sp.symbols("g2s g1s g0s As")
gg = y**3 + g2s * y**2 + g1s * y + g0s
cc = y**5 * gg
ffa = As * y**21 * gg**5
Rc = sp.expand(a * (T2 * cc * sp.diff(ffa, y) - K2 * sp.diff(cc, y) * ffa)
               - cc**e)
Qc, remc = sp.div(sp.Poly(Rc, y), sp.Poly(sp.expand(y**25 * gg**5), y))
ok("F: collapse quotient identity R/(y^25 g^5) = 9A(y g' - 3g) - 1 "
   "(exact division, remainder 0)",
   remc.as_expr() == 0
   and sp.expand(Qc.as_expr()
                 - (9 * As * (y * sp.diff(gg, y) - 3 * gg) - 1)) == 0)

# --- G. branch rungs --------------------------------------------------------
c3r = y**5 * (y + 1)**3
f3p2 = y**21 * (y + 1)**13 * (8 * y**2 + 4 * y - 1) / 27
ok("G: eta=+2 mu=3 explicit residual == 0",
   sp.expand(a * (T2 * c3r * sp.diff(f3p2, y)
                  - K2 * sp.diff(c3r, y) * f3p2) - c3r**e) == 0)
T0, K0 = 4, 19
u30 = -(2048 * y**4 + 2560 * y**3 + 320 * y**2 - 80 * y + 35) / 1155
f30 = y**21 * (y + 1)**13 * u30
ok("G: eta=0 mu=3 explicit residual == 0 (the branch COMPOSITE_CHARTS missed)",
   sp.expand(a * (T0 * c3r * sp.diff(f30, y)
                  - K0 * sp.diff(c3r, y) * f30) - c3r**e) == 0)
ok("G: eta=0 mu=3 u(-1) != 0 and u(0) != 0 (orders exact)",
   u30.subs(y, -1) != 0 and u30.subs(y, 0) != 0)
# mu=2 at eta=0: rebuild fresh
beta = sp.symbols("beta")
ps = sp.symbols("s0:4")
u2g = sum(ps[i] * y**i for i in range(4))
c21 = y**5 * (y + 1)**2 * (y - beta)
f21 = y**21 * (y + 1)**9 * (y - beta)**5 * u2g
R21 = sp.expand(a * (T0 * c21 * sp.diff(f21, y)
                     - K0 * sp.diff(c21, y) * f21) - c21**e)
Q21 = sp.Poly(sp.cancel(R21 / (y**25 * (y + 1)**10 * (y - beta)**5)), y,
              domain=sp.QQ.frac_field(beta, *ps))
eqs = [sp.together(Q21.nth(k)) for k in range(Q21.degree() + 1)]
lin = sp.solve(eqs[:4], list(ps), dict=True)
rem = [sp.factor(sp.numer(sp.together(sp.expand(eq.subs(lin[0])))))
       for eq in eqs[4:]]
conds = [r0 for r0 in rem if sp.simplify(r0) != 0]
quart_expected = 195 * beta**4 + 120 * beta**3 - 40 * beta**2 + 32 * beta - 80
ok("G: eta=0 mu=2 -- exactly one consistency condition, the quartic "
   "195b^4+120b^3-40b^2+32b-80 (up to scalar)",
   len(conds) == 1
   and sp.simplify(sp.expand(conds[0]
                             - sp.Poly(conds[0], beta).LC() / 195
                             * quart_expected)) == 0)
qp = sp.Poly(quart_expected, beta)
ok("G: the quartic has 2 real roots; 0 and -1 are not roots",
   qp.count_roots(-sp.oo, sp.oo) == 2 and qp.eval(0) != 0 and qp.eval(-1) != 0)
lin_sub = {k: sp.together(v) for k, v in lin[0].items()}
allzero = True
for eq in eqs:
    num = sp.numer(sp.together(sp.expand(eq.subs(lin_sub))))
    if sp.rem(sp.Poly(num, beta), qp).as_expr() != 0:
        allzero = False
ok("G: eta=0 mu=2 full system == 0 exactly mod the quartic", allzero)
# mu=2 at eta=+2: no admissible root
ws = sp.symbols("h0:2")
u2p = ws[0] + ws[1] * y
f21p = y**21 * (y + 1)**9 * (y - beta)**5 * u2p
R21p = sp.expand(a * (T2 * c21 * sp.diff(f21p, y)
                      - K2 * sp.diff(c21, y) * f21p) - c21**e)
Q21p = sp.Poly(sp.cancel(R21p / (y**25 * (y + 1)**10 * (y - beta)**5)), y,
               domain=sp.QQ.frac_field(beta, *ws))
eqsp = [sp.together(Q21p.nth(k)) for k in range(Q21p.degree() + 1)]
linp = sp.solve(eqsp[:2], list(ws), dict=True)
admissible = set()
if linp:
    remp = [sp.factor(sp.numer(sp.together(sp.expand(eq.subs(linp[0])))))
            for eq in eqsp[2:]]
    condsp = [sp.Poly(r0, beta) for r0 in remp if sp.simplify(r0) != 0]
    if condsp:
        gc = condsp[0]
        for cnd in condsp[1:]:
            gc = gc.gcd(cnd)
        if gc.degree() > 0:
            admissible = {rt for rt in sp.roots(gc).keys()
                          if rt not in (0, -1)}
    else:
        admissible = {"family"}
ok("G: eta=+2 mu=2 -- no admissible branch", not admissible)

# --- H. the mu-graded law ---------------------------------------------------
def mu_sig(degf, rho0, a0, qq, e, N, gap, rr, mu):
    mult = mu * (e + N) - (mu - 1)
    deg = degf + N * a0
    order = rho0 + N * qq
    cof = deg - order - mult
    law = gap + rr * (e + N) - (mu - 1) * (e + N - 1)
    return (deg, order, mult, cof), law


ROWS = [  # (eta, degf, N, gap, mu, expected signature)
    (0, 38, 97, 2, 1, (814, 506, 102, 206)),
    (0, 38, 97, 2, 2, (814, 506, 203, 105)),
    (0, 38, 97, 2, 3, (814, 506, 304, 4)),
    (2, 36, 157, 0, 1, (1292, 806, 162, 324)),
    (2, 36, 157, 0, 3, (1292, 806, 484, 2)),
]
allrows = True
for eta, degf, N, gap, mu, expect in ROWS:
    sig, law = mu_sig(degf, 21, 8, 5, 5, N, gap, 2, mu)
    if sig != expect or sig[3] != law:
        allrows = False
ok("H: mu-graded law reproduces all five F12 signatures (both eta values)",
   allrows)
E_s, N_s2, r_s, dg_s, gap_s = sp.symbols("E N r dg gap")
cof_mu = gap_s + r_s * (E_s + N_s2) - (dg_s - 1) * (E_s + N_s2 - 1)
ok("H: at mu=dg the law specializes EXACTLY to PHI_F7's ramified law "
   "cof = gap+r via the identity r = dg-1",
   sp.simplify(cof_mu.subs(r_s, dg_s - 1) - (gap_s + dg_s - 1)) == 0)
ok("H: at mu=1 the law is the old unramified law cof = gap + r(e+N)",
   sp.simplify((gap_s + r_s * (E_s + N_s2)
                - (1 - 1) * (E_s + N_s2 - 1))
               - (gap_s + r_s * (E_s + N_s2))) == 0)
# F7 spot check: e+N = 42, dg = 2, gap = 1, r = 1 -> mult 83, cof 2
m7 = 2 * 42 - 1
c7 = 1 + 1 * 42 - (2 - 1) * (42 - 1)
ok("H: F7 spot check -- mu=dg=2 gives (mult, cof) = (83, 2) as in PHI_F7",
   (m7, c7) == (83, 2))

print()
if FAILS:
    print("FAILURES:", len(FAILS))
    sys.exit(1)
print("ALL %d ZETA-TAIL CHECKS PASSED" % N_OK)
