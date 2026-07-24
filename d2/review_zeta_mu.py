#!/usr/bin/env python3
"""review_zeta_mu.py -- LANE R (skeptic): independent verification battery for
the adversarial review of ZETA_TAIL.md / the mu-ladder (REVIEW_ZETA_MU.md).

Independence notes: item D re-derives the F12 mu=2 quartic by a DIFFERENT
elimination (augmented-matrix maximal-minor rank conditions, then gcd), not
the author's solve-first-k-coefficients route.  Exit 0 iff every
confirmation below holds.  --quiet supported.
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

# ---------------------------------------------------------------------------
# A. Rigidity top slice: [x^zeta C^a, C^b] head-head term, three (a,b), two t
# ---------------------------------------------------------------------------
zeta = sp.Symbol("zeta", positive=True)
cf = sp.Function("c")(y)

def bracket(P, Q):
    return sp.expand(sp.diff(P, x) * sp.diff(Q, y) - sp.diff(P, y) * sp.diff(Q, x))

for (a_i, b_i) in ((2, 3), (3, 7), (2, 5)):
    for t_i in (101, 103):
        C = x**t_i * cf
        # full generic scalar tails on both sides
        mus = sp.symbols(f"m1:{a_i}")
        nus = sp.symbols(f"n1:{b_i}")
        P = x**zeta * (C**a_i + sum(mus[j - 1] * C**j for j in range(1, a_i)))
        Q = C**b_i + sum(nus[j - 1] * C**j for j in range(1, b_i))
        BB = bracket(P, Q)
        top = zeta * b_i * x**(zeta + (a_i + b_i) * t_i - 1) * cf**(a_i + b_i - 1) * sp.diff(cf, y)
        rest = sp.expand(BB - top)
        deg_target = (a_i + b_i) * t_i - 1
        bad = [m for m in sp.Add.make_args(rest)
               if (m.as_powers_dict().get(x, 0) - zeta - deg_target).is_zero]
        ok(f"A: (a,b)=({a_i},{b_i}), t={t_i}: top slice is exactly "
           f"zeta*b*x^(zeta+(a+b)t-1)*c^(a+b-1)*c' and nothing else lands there",
           len(bad) == 0)
# slice height above kappa: zeta+(a+b)t-1 - (t-2) = zeta+1+(a+b-1)t > 0
a_s, b_s, t_s = sp.symbols("a_s b_s t_s", positive=True)
ok("A: slice sits kappa + [zeta+1+(a+b-1)t] (symbolic bookkeeping)",
   sp.simplify((zeta + (a_s + b_s) * t_s - 1) - (t_s - 2)
               - (zeta + 1 + (a_s + b_s - 1) * t_s)) == 0)
# max cross x-degree from head_P x tail_Q is zeta+(a+b-1)t-1 < slice
ok("A: tails are strictly below the slice (head x tail max index a+b-1)",
   sp.simplify((zeta + (a_s + b_s) * t_s - 1)
               - (zeta + (a_s + b_s - 1) * t_s - 1)) == t_s)

# ---------------------------------------------------------------------------
# B. eta = -1 lemma: K = eT identity + integration + residues
# ---------------------------------------------------------------------------
T, tg, ag, bg = sp.symbols("T t_g a_g b_g", positive=True)
K_expr = T * (bg - ag) + (tg - 2) + 1
e_expr = bg - ag + 1
ok("B: K = eT  <=>  T = t-1  (both directions, symbolic)",
   sp.simplify(K_expr.subs(T, tg - 1) - e_expr * (tg - 1)) == 0
   and sp.solve(sp.Eq(K_expr, e_expr * T), T) == [tg - 1])
ff = sp.Function("f")(y)
for e_i in (5, 12):  # F12, F13 values of e
    lhs = sp.diff(ff * cf**(-e_i), y)
    rhs = (T * cf * sp.diff(ff, y) - e_i * T * sp.diff(cf, y) * ff) / (T * cf**(e_i + 1))
    ok(f"B: ODE integrates to (f/c^e)' = 1/(aTc) at e={e_i}",
       sp.simplify(lhs - rhs) == 0)
al = sp.Symbol("alpha", nonzero=True)
be = sp.Symbol("beta")
# simple root residue never vanishes: c = y^q (y-al) * h(y), h(al) != 0
# Res_{al} 1/c = 1/c'(al) with c'(al) = al^q h(al) != 0  -- structural.
for q_i, dg_i in ((5, 3), (7, 2)):
    if dg_i == 3:
        c21 = y**q_i * (y - al)**2 * (y - be)
        R_al = sp.residue(1 / c21, y, al)
        be_tuned = sp.solve(sp.numer(sp.together(R_al)), be)
        ok(f"B: q={q_i} [2,1]: Res_alpha = 0 has a unique tuning beta(alpha)",
           len(be_tuned) == 1)
        R_be = sp.simplify(sp.residue(1 / c21, y, be).subs(be, be_tuned[0]))
        ok(f"B: q={q_i} [2,1]: tuned Res_beta = {R_be} is a nonzero multiple "
           f"of a power of alpha", sp.simplify(R_be) != 0
           and sp.simplify(R_be * al**7).is_constant())
        c3 = y**q_i * (y - al)**3
        R3 = sp.residue(1 / c3, y, al)
        ok(f"B: q={q_i} [3]: Res_alpha = {sp.simplify(R3)} != 0",
           sp.simplify(R3) != 0)
    else:
        c2 = y**q_i * (y - al)**2
        R2 = sp.residue(1 / c2, y, al)
        ok(f"B: q={q_i} [2]: Res_alpha = {sp.simplify(R2)} != 0",
           sp.simplify(R2) != 0)

# ---------------------------------------------------------------------------
# C. Specialization algebra: r = dg-1 identity; cof-law redundancy
# ---------------------------------------------------------------------------
a0_s, q_s = sp.symbols("a0 q", positive=True)
ok("C: r = dg - 1 is definitional (r = a0-q-1, dg = a0-q)",
   sp.simplify((a0_s - q_s - 1) - ((a0_s - q_s) - 1)) == 0)
mu, e_s, N_s, gap_s, dg_s = sp.symbols("mu e N gap dg", positive=True)
r_s = dg_s - 1
mult_law = mu * (e_s + N_s) - (mu - 1)
cof_law = gap_s + r_s * (e_s + N_s) - (mu - 1) * (e_s + N_s - 1)
ok("C: mu=1 specializes to the unramified law (mult = e+N, cof = gap+r(e+N))",
   sp.simplify(mult_law.subs(mu, 1) - (e_s + N_s)) == 0
   and sp.simplify(cof_law.subs(mu, 1) - (gap_s + r_s * (e_s + N_s))) == 0)
ok("C: mu=dg specializes to PHI_F7's ramified law (mult = dg(e+N)-(dg-1), "
   "cof = gap+r)",
   sp.simplify(mult_law.subs(mu, dg_s) - (dg_s * (e_s + N_s) - (dg_s - 1))) == 0
   and sp.simplify(cof_law.subs(mu, dg_s) - (gap_s + r_s)) == 0)
# cof-law is REDUNDANT given deg/ord: deg-ord = (e+N)dg + gap identically.
# Work entirely in (a0, q, e, N, gap, mu): dg := a0-q, r := a0-q-1.
rho0 = (e_s - 1) * q_s + 1
res_deg = (e_s * a0_s - q_s + 1) + gap_s
deg_m_ord = sp.expand((res_deg + N_s * a0_s) - (rho0 + N_s * q_s))
mult_l = mu * (e_s + N_s) - (mu - 1)
cof_l = gap_s + (a0_s - q_s - 1) * (e_s + N_s) - (mu - 1) * (e_s + N_s - 1)
ok("C: deg - ord = (e+N)dg + gap identically, so cof-law == (deg-ord) - mult-law "
   "(the cof formula carries no independent content)",
   sp.simplify(deg_m_ord - ((e_s + N_s) * (a0_s - q_s) + gap_s)) == 0
   and sp.simplify((deg_m_ord - mult_l) - cof_l) == 0)

# ---------------------------------------------------------------------------
# D. F12 mu=2 at eta=0: INDEPENDENT elimination (augmented-minor route)
# ---------------------------------------------------------------------------
a_f, b_f, t_f, kap_f, a0_f, q_f = 3, 7, 4, 2, 8, 5
e_f = b_f - a_f + 1
T0, K0 = t_f, t_f * (b_f - a_f) + kap_f + 1          # 4, 19
N0 = a_f * (t_f * (a_f + b_f - 1) + 1) - 2 * b_f     # 97
ok("D: F12 eta=0 bookkeeping T=4, K=19, N=97, e+N=102",
   (T0, K0, N0, e_f + N0) == (4, 19, 97, 102))
beta = sp.Symbol("beta")
p0, p1, p2, p3 = sp.symbols("p0:4")
u = p0 + p1 * y + p2 * y**2 + p3 * y**3
c = y**q_f * (y + 1)**2 * (y - beta)
f = y**21 * (y + 1)**9 * (y - beta)**5 * u
R = sp.expand(a_f * (T0 * c * sp.diff(f, y) - K0 * sp.diff(c, y) * f) - c**e_f)
Qp = sp.Poly(sp.cancel(R / (y**25 * (y + 1)**10 * (y - beta)**5)), y,
             domain=sp.QQ.frac_field(beta, p0, p1, p2, p3))
eqs = [sp.expand(sp.numer(sp.together(Qp.nth(k)))) for k in range(Qp.degree() + 1)]
# Naive top degree of the quotient is 5, but deg f = 38 = d_res makes the
# top coefficient cancel identically (T*38 = K*8 = 152) -- the resonance that
# defines d_res.  So the quotient has degree 4: FIVE coefficient equations.
ok("D: residual quotient degree is 4 (resonant top-coefficient cancellation "
   "T*38 = K*8), giving 5 coefficient equations",
   len(eqs) == 5 and T0 * 38 == K0 * 8)
# linear system M(beta) * (p0..p3) = v(beta): rank condition via 5x5 minors of [M|v]
M = sp.zeros(len(eqs), 4)
v = sp.zeros(len(eqs), 1)
for i, eq in enumerate(eqs):
    peq = sp.Poly(eq, p0, p1, p2, p3)
    for j, pj in enumerate((p0, p1, p2, p3)):
        M[i, j] = peq.diff(pj).as_expr()
    v[i] = -peq.subs({p0: 0, p1: 0, p2: 0, p3: 0})
ok("D: system is affine-linear in the u-coefficients",
   all(sp.degree(sp.Poly(eq, p0, p1, p2, p3), gen=g) <= 1
       for eq in eqs for g in (p0, p1, p2, p3)))
Aug = M.row_join(v)
import itertools
minors = []
for rows in itertools.combinations(range(len(eqs)), 5):
    d = Aug[list(rows), :].det()
    if sp.simplify(d) != 0:
        minors.append(sp.Poly(sp.factor(d), beta))
g = minors[0]
for m in minors[1:]:
    g = g.gcd(m)
# strip beta / (beta+1) unit factors (excluded root positions)
g_core = g.as_expr()
for excl in (beta, beta + 1):
    while sp.simplify(sp.rem(sp.Poly(g_core, beta), sp.Poly(excl, beta)).as_expr()) == 0:
        g_core = sp.cancel(g_core / excl)
quart_claimed = 195 * beta**4 + 120 * beta**3 - 40 * beta**2 + 32 * beta - 80
ratio = sp.simplify(sp.cancel(g_core / quart_claimed))
ok("D: augmented-minor gcd reproduces the claimed quartic "
   "195b^4+120b^3-40b^2+32b-80 up to a nonzero rational scalar (and unit "
   "factors at the excluded positions 0, -1)",
   ratio.is_constant() and ratio != 0)
qP = sp.Poly(quart_claimed, beta)
ok("D: quartic has exactly 2 real roots; 0 and -1 are not roots",
   qP.count_roots(-sp.oo, sp.oo) == 2 and qP.eval(0) != 0 and qP.eval(-1) != 0)
# exact branch check at the algebraic root
brt = sp.RootOf(quart_claimed, 0)
minp = sp.Poly(quart_claimed, beta)
sol_u = sp.solve(eqs[:4], [p0, p1, p2, p3], dict=True)
ok("D: first four coefficients determine u uniquely (rational in beta)",
   len(sol_u) == 1)
subs_u = sol_u[0]
resid_mod = []
for eq in eqs:
    num = sp.numer(sp.together(sp.expand(eq.subs(subs_u))))
    resid_mod.append(sp.rem(sp.Poly(num, beta), minp).as_expr())
ok("D: FULL system vanishes exactly mod the quartic (all 6 equations)",
   all(sp.simplify(r0) == 0 for r0 in resid_mod))
u_at = sp.expand(u.subs(subs_u))
vals = {}
for tag, point in (("u(-1)", -1), ("u(0)", 0)):
    val_num = sp.numer(sp.together(u_at.subs(y, point)))
    vals[tag] = sp.rem(sp.Poly(sp.expand(val_num), beta), minp).as_expr()
lc_u = sp.numer(sp.together(u_at.coeff(y, 3)))
vals["lc(u)"] = sp.rem(sp.Poly(sp.expand(lc_u), beta), minp).as_expr()
ok("D: u(-1), u(0), lc(u) all nonzero mod the quartic (orders are exact, "
   "deg f = 38)", all(sp.simplify(v0) != 0 for v0 in vals.values()))
# signature bookkeeping: mult = 9 + 2N = 203; deg = 38+8N = 814; ord = 21+5N = 506
ok("D: signature (814, 506, 203, 105) is exact bookkeeping for the mu=2 branch",
   (38 + 8 * N0, 21 + 5 * N0, 9 + 2 * N0,
    (38 + 8 * N0) - (21 + 5 * N0) - (9 + 2 * N0)) == (814, 506, 203, 105))

# ---------------------------------------------------------------------------
# E. F12 mu=3 at eta=0: independent rebuild
# ---------------------------------------------------------------------------
z = sp.symbols("z0:5")
u4 = sum(z[i] * y**i for i in range(5))
c3 = y**5 * (y + 1)**3
f3 = y**21 * (y + 1)**13 * u4
R3 = sp.expand(a_f * (T0 * c3 * sp.diff(f3, y) - K0 * sp.diff(c3, y) * f3) - c3**e_f)
Q3 = sp.Poly(sp.cancel(R3 / (y**25 * (y + 1)**15)), y)
sol3 = sp.solve([Q3.nth(k) for k in range(Q3.degree() + 1)], list(z), dict=True)
u3_sols = [sp.expand(u4.subs(s)) for s in sol3 if sp.expand(u4.subs(s)) != 0]
u3_claim = -(2048 * y**4 + 2560 * y**3 + 320 * y**2 - 80 * y + 35) / 1155
ok("E: mu=3 branch at eta=0 exists, is unique, and matches the claimed u",
   len(u3_sols) == 1 and sp.simplify(u3_sols[0] - u3_claim) == 0)
ok("E: mu=3 u is a unit at 0 and -1 (orders exact -> sig (814,506,304,4))",
   u3_claim.subs(y, -1) != 0 and u3_claim.subs(y, 0) != 0
   and (38 + 8 * N0, 21 + 5 * N0, 13 + 3 * N0,
        (38 + 8 * N0) - (21 + 5 * N0) - (13 + 3 * N0)) == (814, 506, 304, 4))

# ---------------------------------------------------------------------------
# F. Latent sweep-logic probe: the sub-resonant "ok=False" evasion branch
#    must never fire on the actually-swept rows (else DEAD verdicts unsound)
# ---------------------------------------------------------------------------
def sweep_evasions(fam, etas):
    a, b, t, kap, a0, q = (fam[k] for k in ("a", "b", "t", "kappa", "a0", "q"))
    e, dg = b - a + 1, a0 - q
    hits = []
    for eta in etas:
        Tx = Fraction(t) + eta
        n = Tx.denominator
        Tw = int(Tx * n)
        Kw = Tw * (b - a) + n * (kap + 1)
        if Kw == e * Tw:
            continue
        partitions = {3: ([1, 1, 1], [2, 1], [3]), 2: ([1, 1], [2])}[dg]
        for part in partitions:
            for mu_i in part:
                om = (e - 1) * mu_i + 1
                sub = Fraction(Kw * mu_i, Tw)
                if sub.denominator == 1 and sub < om:
                    hits.append((fam, eta, part, mu_i))
    return hits

F12d = dict(a=3, b=7, t=4, kappa=2, a0=8, q=5)
F13d = dict(a=2, b=13, t=3, kappa=1, a0=9, q=7)
ev = (sweep_evasions(F12d, sorted(set([Fraction(k, 3) for k in (-3, -2, -1, 1, 2, 3)]
                                      + [0, 2, -2, 3, -3])))
      + sweep_evasions(F13d, sorted(set([Fraction(k, 2) for k in (-2, -1, 1, 2)]
                                        + [0, 2, -2]))))
ok("F: the latent 'sub-resonant evasion' branch (ok=False) never fires on any "
   "swept (family, eta, partition) row -- current DEAD verdicts unaffected",
   len(ev) == 0)

print()
if FAILS:
    print("REVIEW CONFIRMATION FAILURES:", len(FAILS))
    for f0 in FAILS:
        print("  -", f0)
    sys.exit(1)
print(f"ALL {N_OK} REVIEW CONFIRMATIONS HOLD (review_zeta_mu.py)")
