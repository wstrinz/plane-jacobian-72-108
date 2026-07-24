#!/usr/bin/env python3
"""phi_f7.py  (NEW; read-only over all existing artifacts)

The last unprobed regime of the corner law: gap>0 WITH r>0 (F3, F7, F10, F15,
F16) -- where PHI_F14.md flagged the unified cofactor formula gap + r(e+N) as
a CONJECTURE.  Verdict: **DIFFERS -- the conjecture is REFUTED, and a sharper
regime law replaces it.**

Result in brief (all exact, ODE residual == 0 checked for every solution):

  * The deg and ord_y components of the unified law hold on EVERY branch:
        deg Phi = res + N*a0   (res = pure + gap),   ord_y Phi = rho + N*q.
  * The mult and cofactor components FAIL.  Structural reason, proven at
    dg = 2 by exact elimination: the ODE's polynomial-solvability obstruction
    factors as
        F7 :  E ~ g0^27 * (3 g0 - 2 g1^2)   * (4 g0 - g1^2)^6
        F3 :  E ~ g0^6  * (5 g0 - 3 g1^2)   * (4 g0 - g1^2)^2
        F16:  E ~ g0^18 * (54 g0^2 - 126 g0 g1^2 + 35 g1^4) * (4 g0 - g1^2)^3
    so the residual g (monic, deg 2) is either RAMIFIED (double root; the
    discriminant factor) or has a COMPLEX-CONJUGATE pair (the ratio factor(s)
    force disc < 0 over R).  A simple real root at -1 -- what the old
    mult = e+N formula needs -- is impossible.  (dg odd never faces this:
    g = y^dg + 1 has a simple root at -1; dg is even on every gap>0,r>0 row
    of the survey and odd elsewhere -- checked.)
  * On the ramified branch g = (y+1)^dg (the branch continuous with the
    audited pattern (y+1) | C), four fresh exact points obey the AMENDED law
        mult  = dg*(e+N) - (dg-1),      cofactor = gap + r,
    which also retro-explains the audited (72,108) quartic (gap+r = 4+0 = 4):

      F7  (42,147) : f = (1/10)  y^21 (y+1)^11 (9y^2+3y-1)
                     Phi sig (250, 165,  83, 2)   [old law said (250,165,42,43)]
      F3  (75,50)  : f = (1/42)  y^4  (y+1)^3  (25y^2+15y-3)
                     Phi sig (189, 112,  75, 2)   [old law said (189,112,38,39)]
      F10 (196,112): f = (1/3740) y^10 (y+1)^13 (2401y^4+5831y^3+4165y^2+595y-85)
                     Phi sig (1917, 820, 1093, 4)
      F16 (99,165) : f = (1/330) y^15 (y+1)^5  (243y^4+81y^3-27y^2+15y-10)
                     Phi sig (528, 407, 117, 4)

    F16 (gap=3) separates deg u = gap+r from deg u = dg (4 vs 2): gap+r wins.
  * On the complex-pair branch the signature is (deg, ord, 0, deg-ord) --
    no (y+1) place at all; recorded for F7 and F3 with rational scale reps.
  * Control: the same machinery re-derives F14 (gap=0) exactly, and shows the
    F1 (dg=1, r=0) obstruction s^6(15 F s^2 - 4) is solvable for EVERY root
    position s != 0 -- confirming PHI_F14.md's judgment that dg=1 root
    placement is gauge, while dg=2 root RATIOS are forced.

Branch selection is a judgment item (see PHI_F7.md): the ramified branch is
what continuity with all seven previous points selects, but the actual tower
C-series for these families is built in no paper.

Sources: GGV5 family tables as transcribed in phi_corner4.py (Diophantine
re-checked here); method template PHI_F14.md / PHI_CORNER4.md.  The
independent PASS/FAIL checker is phi_f7_verify.py.  Everything is exact sympy.
"""
import sympy as sp
from fractions import Fraction
from math import gcd

y = sp.symbols("y")

# ---------------------------------------------------------------------------
# 0. Corner rows (GGV5 v11<=35 tables; same transcription as phi_corner4.py)
#    (name, A0, p, l, q, k, (m0,dm), (n0,dn))
# ---------------------------------------------------------------------------
ROWS = {
    "F7":  ((6, 15), 7, 3, 4, 1, (2, 1), (7, 4)),
    "F3":  ((5, 20), 8, 5, 3, 1, (3, 4), (2, 3)),
    "F10": ((7, 21), 13, 7, 3, 1, (7, 5), (4, 3)),
    "F16": ((9, 24), 10, 3, 7, 1, (3, 4), (5, 7)),
    "F14": ((9, 24), 7, 3, 4, 1, (2, 1), (7, 4)),
    "F1":  ((4, 12), 7, 4, 3, 1, (3, 2), (4, 3)),
}

def corner_params(name):
    A0, p, l, q, k, (m0, dm), (n0, dn) = ROWS[name]
    j = 0
    while gcd(m0 + dm * j, n0 + dn * j) != 1:
        j += 1
    m, n = m0 + dm * j, n0 + dn * j
    assert (m + n) * q * k - n * (q * l - p) == k, f"{name}: Diophantine failed"
    v11 = A0[0] + A0[1]
    a, b = sorted((m, n))
    t, kappa = l, l - 2               # standard chart (A0'=(1,0), PHI_CORNER4 sec.2)
    a0 = A0[0]
    e, r, dg = b - a + 1, a0 - q - 1, a0 - q
    coef = t * (b - a) + kappa + 1
    rho = (e - 1) * q + 1
    N = a * (t * (a + b) - (kappa + 1)) - 2 * b
    res = Fraction(coef * a0, t)
    assert res.denominator == 1
    res = int(res)
    gap = res - (e * a0 - q + 1)
    return dict(name=name, j=j, m=m, n=n, degs=(v11 * m, v11 * n), a=a, b=b,
                t=t, kappa=kappa, a0=a0, q=q, e=e, r=r, dg=dg, coef=coef,
                rho=rho, N=N, res=res, gap=gap)

def old_law_sig(P):
    """PHI_F14's unified law -- the conjecture under test in this regime."""
    deg = (P["e"] * P["a0"] - P["q"] + 1) + P["gap"] + P["N"] * P["a0"]
    return (deg, P["rho"] + P["N"] * P["q"], P["e"] + P["N"],
            P["gap"] + P["r"] * (P["e"] + P["N"]))

def ram_law_sig(P):
    """The amended ramified-branch law established here."""
    deg = P["res"] + P["N"] * P["a0"]
    mult = P["dg"] * (P["e"] + P["N"]) - (P["dg"] - 1)
    return (deg, P["rho"] + P["N"] * P["q"], mult, P["gap"] + P["r"])

# ---------------------------------------------------------------------------
# 1. The ODE in coefficient form + exact triangular solver
#    L(f) = a*t*c*f' - a*coef*c'*f,  c = y^q g,  g = sum g_i y^i monic deg dg.
#    L(f) = a * sum_{d,i} [t*d - coef*(q+i)] g_i f_d y^{d+i+q-1}, so with
#    shifted index k (y-power K = k+q-1) the equation reads
#        a * sum_i g_i [t(k-i) - coef(q+i)] f_{k-i} = [g^e]_{k+q-1-q*e}.
#    Top-down pivot f_{k-dg} (g_dg = 1) has integer coefficient
#    a*t*(k-dg-res); it vanishes exactly at the resonant degree res.
# ---------------------------------------------------------------------------
def solve_corner(P, g_coeffs, slack=2):
    a, t, coef, q, e, dg, res = (P[k] for k in
                                 ("a", "t", "coef", "q", "e", "dg", "res"))
    g_list = list(g_coeffs) + [sp.Integer(1)]
    ge = sp.Poly(sp.expand(sum(g_list[i] * y**i for i in range(dg + 1))**e), y)
    def rhs(kk):
        idx = kk + q - 1 - q * e
        if idx < 0 or idx > e * dg:
            return sp.Integer(0)
        return ge.coeff_monomial(y**idx) or sp.Integer(0)
    Dmax = res + slack
    f, frees, conditions = {}, [], []
    for k in range(Dmax + dg, -1, -1):
        expr = -rhs(k)
        piv_idx, piv_coeff = k - dg, a * (t * (k - dg) - coef * (q + dg))
        assert piv_coeff == a * t * (k - dg - res)      # resonance bookkeeping
        for i in range(dg + 1):
            d = k - i
            if d < 0 or d > Dmax or (d == piv_idx and piv_coeff != 0 and d not in f):
                continue
            if d not in f:                              # resonant free coefficient
                fs = sp.symbols(f"FREE{d}")
                frees.append((d, fs))
                f[d] = fs
            expr += a * (t * d - coef * (q + i)) * g_list[i] * f[d]
        if 0 <= piv_idx <= Dmax and piv_idx not in f:
            if piv_coeff != 0:
                f[piv_idx] = sp.expand(sp.cancel(-expr / piv_coeff))
            else:
                fs = sp.symbols(f"FREE{piv_idx}")
                frees.append((piv_idx, fs))
                f[piv_idx] = fs
                if sp.expand(expr) != 0:
                    conditions.append((k, sp.expand(expr)))
        elif sp.expand(expr) != 0:
            conditions.append((k, sp.expand(expr)))
    return f, dict(frees), conditions, Dmax

def ode_residual(P, g_coeffs, fdict):
    g = sum((list(g_coeffs) + [sp.Integer(1)])[i] * y**i
            for i in range(P["dg"] + 1))
    c = y**P["q"] * g
    fpoly = sum(v * y**d for d, v in fdict.items())
    return sp.expand(P["a"] * P["t"] * c * sp.diff(fpoly, y)
                     - P["a"] * P["coef"] * sp.diff(c, y) * fpoly - c**P["e"])

def solve_branch(P, g_coeffs):
    """Fix g, solve the two consistency conditions for the resonant free
    coefficient, return (f_poly, signature_of_f)."""
    fdict, frees, conds, _ = solve_corner(P, g_coeffs)
    Fs = frees[P["res"]]
    sol = sp.solve([c for _, c in conds], Fs, dict=True)
    assert len(sol) == 1, f"{P['name']}: branch not uniquely solvable"
    fnum = {d: sp.expand(v.subs(sol[0])) for d, v in fdict.items()}
    assert ode_residual(P, g_coeffs, fnum) == 0
    return sum(v * y**d for d, v in fnum.items())

def signature_poly(p):
    P2 = sp.Poly(sp.expand(p), y)
    deg = P2.degree()
    ordy = min(mm[0] for mm in P2.monoms())
    mult, Q = 0, P2
    while True:
        q2, r2 = sp.div(Q, sp.Poly(y + 1, y))
        if not r2.is_zero:
            break
        Q, mult = q2, mult + 1
    return (deg, ordy, mult, deg - ordy - mult)

def signature_Phi(P, g_coeffs, fpoly):
    """Signature of Phi = f * c^N without expanding c^N: exponent arithmetic
    on top of f's trial-division signature (validated against a full
    expansion at F3 in phi_f7_verify.py)."""
    df, of, mf, cf = signature_poly(fpoly)
    g = sum((list(g_coeffs) + [sp.Integer(1)])[i] * y**i
            for i in range(P["dg"] + 1))
    dgs, ogs, mgs, cgs = signature_poly(y**P["q"] * g)
    N = P["N"]
    return (df + N * dgs, of + N * ogs, mf + N * mgs, cf + N * cgs)

# ---------------------------------------------------------------------------
# 2. F7: pure-ansatz failure, obstruction, both branches
# ---------------------------------------------------------------------------
def run_dg2(name, ratio_reps):
    P = corner_params(name)
    print("=" * 96)
    print(f"{name} j={P['j']}, case {P['degs']}: (a,b)=({P['a']},{P['b']}) "
          f"t={P['t']} kappa={P['kappa']} a0={P['a0']} q={P['q']} e={P['e']} "
          f"r={P['r']} dg={P['dg']} coef={P['coef']} rho={P['rho']} N={P['N']} "
          f"res={P['res']} gap={P['gap']}")
    print("=" * 96)
    tr, te = P["t"] * P["rho"] - P["coef"] * P["q"], P["t"] * P["e"] - P["coef"]
    print(f"  pure ansatz f = A y^{P['rho']} g^{P['e']}: bracket multipliers "
          f"{tr}+i*{te} for i=1..{P['dg']} are "
          f"{[tr + i * te for i in range(1, P['dg'] + 1)]} -- all nonzero, so "
          f"g_1..g_{P['dg']} = 0 is forced: CONTRADICTION with monic deg "
          f"{P['dg']}.  Resonance is broken (gap={P['gap']}); generic f needed.")
    g1, g0 = sp.symbols("g1 g0")
    fdict, frees, conds, _ = solve_corner(P, [g0, g1])
    Fs = frees[P["res"]]
    assert len(conds) == 2
    pa, pb = (sp.Poly(c, Fs) for _, c in conds)
    E = sp.factor(sp.expand(pa.all_coeffs()[0] * pb.all_coeffs()[1]
                            - pb.all_coeffs()[0] * pa.all_coeffs()[1]))
    print(f"  generic g = y^2+g1*y+g0: 2 consistency conditions, linear in the")
    print(f"  resonant free f_{P['res']}; eliminant E(g1,g0) factors as")
    print(f"    E = {E}")
    L = sp.symbols("L", positive=True)
    ratio = sp.simplify(sp.expand(E.subs({g1: L * g1, g0: L**2 * g0},
                                         simultaneous=True)) / sp.expand(E))
    assert ratio.is_Pow or ratio.is_Symbol or ratio.is_number
    print(f"  (weighted-homogeneous: scale ratio {ratio} -- root RATIOS are "
          f"forced, unlike dg=1)")
    print(f"  Real branches: RAMIFIED g=(y+1)^2 (disc factor) or COMPLEX PAIR")
    print(f"  (ratio factor(s) force disc<0).  Simple real root at -1: IMPOSSIBLE.")
    out = {}
    # ramified branch
    f_ram = solve_branch(P, [sp.Integer(1), sp.Integer(2)])
    sig_ram = signature_Phi(P, [sp.Integer(1), sp.Integer(2)], f_ram)
    print(f"  RAMIFIED branch: f = {sp.factor(f_ram)}")
    print(f"    Phi sig = {sig_ram};  old law {old_law_sig(P)}  "
          f"ram law {ram_law_sig(P)}")
    assert sig_ram == ram_law_sig(P)
    out["ram"] = (f_ram, sig_ram)
    # complex-pair branch, rational scale representative(s)
    for g1v, g0v in ratio_reps:
        assert sp.expand(E.subs({g1: g1v, g0: g0v})) == 0
        f_cx = solve_branch(P, [g0v, g1v])
        sig_cx = signature_Phi(P, [g0v, g1v], f_cx)
        print(f"  COMPLEX-PAIR branch (rep g1={g1v}, g0={g0v}): "
              f"f = {sp.factor(f_cx)}")
        print(f"    Phi sig = {sig_cx}  (no (y+1) place)")
        assert sig_cx[2] == 0
        assert sig_cx[0] == ram_law_sig(P)[0] and sig_cx[1] == ram_law_sig(P)[1]
        out["cx"] = (f_cx, sig_cx)
    return P, out

P7, out7 = run_dg2("F7", [(sp.Integer(3), sp.Integer(6))])     # 3g0=2g1^2
print()
P3, out3 = run_dg2("F3", [(sp.Integer(3), sp.Rational(27, 5))])  # 5g0=3g1^2
print()
P16, out16 = run_dg2("F16", [])   # ratio factor quadratic in g0/g1^2: skip reps

# ---------------------------------------------------------------------------
# 3. dg=4: F10 ramified-branch point (branch completeness NOT claimed at dg=4)
# ---------------------------------------------------------------------------
print()
print("=" * 96)
P10 = corner_params("F10")
print(f"F10 j={P10['j']}, case {P10['degs']}: (a,b)=({P10['a']},{P10['b']}) "
      f"t={P10['t']} a0={P10['a0']} q={P10['q']} e={P10['e']} r={P10['r']} "
      f"dg={P10['dg']} rho={P10['rho']} N={P10['N']} res={P10['res']} "
      f"gap={P10['gap']}")
print("=" * 96)
gc10 = [sp.Integer(1), sp.Integer(4), sp.Integer(6), sp.Integer(4)]  # (y+1)^4
f10 = solve_branch(P10, gc10)
sig10 = signature_Phi(P10, gc10, f10)
print(f"  RAMIFIED branch g=(y+1)^4: f = {sp.factor(f10)}")
print(f"  Phi sig = {sig10};  old law {old_law_sig(P10)}  ram law {ram_law_sig(P10)}")
assert sig10 == ram_law_sig(P10)

# ---------------------------------------------------------------------------
# 4. Controls: F14 (gap=0) re-derivation; F1 (dg=1) root-gauge freedom
# ---------------------------------------------------------------------------
print()
print("=" * 96)
print("CONTROLS")
print("=" * 96)
P14 = corner_params("F14")
f14 = solve_branch(P14, [sp.Integer(1), 0, 0, 0, 0])           # g = y^5+1
assert sp.expand(f14 - (sp.Rational(-1, 10) * y**21 * (y**5 + 1)**6)) == 0
print("  F14 (gap=0): machinery re-derives f = -(1/10) y^21 (y^5+1)^6 exactly")
P1 = corner_params("F1")
s = sp.symbols("s")
_, frees1, conds1, _ = solve_corner(P1, [s])
assert len(conds1) == 1
c1 = sp.factor(conds1[0][1])
print(f"  F1 (dg=1): single condition {c1} -- solvable for EVERY s != 0:")
print("  root position is GAUGE at dg=1 (PHI_F14 judgment confirmed);")
print("  at dg=2 the ratio is forced -- the new phenomenon of this regime.")

# dg parity observation over the whole survey
SURVEY = ["F1", "F3", "F7", "F10", "F14", "F16"]
for nm in SURVEY:
    Pn = corner_params(nm)
    even = Pn["dg"] % 2 == 0
    in_regime = Pn["gap"] > 0 and Pn["r"] > 0
    assert even == in_regime, nm
print("  dg-parity: dg even <=> (gap>0 and r>0) on every family used here")

# ---------------------------------------------------------------------------
# 5. Verdict
# ---------------------------------------------------------------------------
print()
print("=" * 96)
print("VERDICT")
print("=" * 96)
for Pn, on in ((P7, out7), (P3, out3), (P16, out16)):
    print(f"  {Pn['name']}: ramified sig {on['ram'][1]}  vs old-law "
          f"{old_law_sig(Pn)}  -> DIFFERS (deg, ord match; mult, cof do not)")
print(f"  F10: ramified sig {sig10}  vs old-law {old_law_sig(P10)}  -> DIFFERS")
print("""  The gap>0,r>0 regime REFUTES the old unified mult/cofactor conjecture and
  obeys the amended ramified law
      deg  = res + N*a0          ord = rho + N*q        (unchanged)
      mult = dg*(e+N) - (dg-1)   cofactor = gap + r     (new; 4 exact points)
  with cofactor = gap + r also retro-explaining the audited (72,108) quartic
  (gap+r = 4).  Structural driver, proven at dg=2: the ODE obstruction only
  admits ramified or complex-pair residuals -- a simple real root at -1 is
  impossible, so mult = e+N cannot occur.""")
