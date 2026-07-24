#!/usr/bin/env python3
"""c_series_75_125.py  (NEW; read-only over all existing artifacts)

THE TRANSFER TEST, phase 1 -- BUILD the actual C-series for the (75,125) case
from its polygon/chain data, then DERIVE the tower length N from the built
D-transform tower (not from the N-formula), and read off Phi = f * C^N.

This is the pending item flagged in PHI_75_125.md judgment item 3: N = 98 was a
FORMULA-based extrapolation of the corner-144 clearing exponent to (t,kappa) =
(5,3); the concern was that the corner-144 derivation used a per-term "forcing
slice index" k = t - b*t + s + (b-1)/a that is INTEGER only when a | (b-1), and
for (75,125) (b-1)/a = 4/3 is not integral.  This script builds the tower and
shows N is a slice-SUM invariant that does not see that per-term index at all.

Case (family F_2, j=1): corner A_0 = (5,20), final corner (7/5,2), (m,n)=(3,5),
reduced C-power pair (a,b)=(3,5), t=5, kappa=3, a0=5, q=2.

Independent PASS/FAIL checker: c_series_75_125_verify.py (--quiet, exit 0).
Everything below is exact sympy.  Judgment items are called out inline.
"""
import sympy as sp

y, u, c = sp.symbols("y u c")


def order(e):
    return min(m[0] for m in sp.Poly(sp.expand(e), y).monoms())


def multiplicity(e, fac):
    cnt, qq, d = 0, sp.Poly(sp.expand(e), y), sp.Poly(fac, y)
    while True:
        qq, rem = sp.div(qq, d)
        if not rem.is_zero:
            return cnt
        cnt += 1


def c_exponent(term):
    _, cpart = term.as_independent(c, as_Add=False)
    if cpart == 1:
        return 0
    base, expo = cpart.as_base_exp()
    assert base == c
    return int(expo)


# ===========================================================================
# 0. Corner geometry -> reduction parameters.
#    [judgment 1: chain data] read from GGV5 line 1679 (F_2 j=1 table row).
#    [judgment 2: unreduced polygon] the (5,20) reduction is in no paper; we
#    assume the standard type-II.b root shift + final Laurent chart, exactly as
#    CORNER_144_COMPARISON.md does for its own conditional boundary.
# ===========================================================================
a, b = 3, 5              # P = C^a, Q = C^b + (commuting C-powers) + F
a0, l, q = 5, 5, 2       # deg C = a0 ; Laurent denominator l ; selected mult q
t = l                    # ell(C) = x^t c
kappa = l - 2            # [P,Q] = x^kappa ; Jacobian of (x^-1, x^l y) is -x^(l-2)
e = b - a + 1            # forcing RHS exponent c^e  (= 3; non-adjacent pair)
r = a0 - q - 1           # residual degree in the leading form
s = kappa + 1 - a * t    # v(F), the obstruction valuation
j = -s                   # forcing-slice u-power offset

print("=" * 78)
print("(75,125)  F_2 j=1   corner (5,20) -> (7/5,2)   (m,n)=(3,5)  (a,b)=(3,5)")
print("=" * 78)
print(f"  t=l={t}   kappa=l-2={kappa}   e=b-a+1={e}   residual deg r={r}   deg C=a0={a0}")
print(f"  s=v(F)=kappa+1-a*t={s}   j=-s={j}")
print("  [judgment 1] chain data: GGV5 line 1679 (F_2 j=1), Diophantine-checked")
print("  [judgment 2] unreduced polygon: standard type-II.b shift + Laurent chart assumed")

# ===========================================================================
# 1. Build the C-series leading polynomial by solving the forcing ODE.
#    The operator identity [P, x^s f/c^b] = x^kappa gives
#        a { t c f' - [t(b-a)+kappa+1] c' f } = c^(b-a+1),   f = c^b F_s.
#    Here:  15 c f' - 42 c' f = c^3,   c = y^q g,  deg g = a0-q = 3.
# ===========================================================================
coef = t * (b - a) + kappa + 1              # = 42/a ... actually a*coef below
rho = (e - 1) * q + 1                        # local order of f at y=0  (= 5)
dg = a0 - q                                  # deg g  (= 3)
print("\n--- 1. C-series ODE ---------------------------------------------------")
print(f"  forcing ODE:  {a*t} c f' - {a*coef} c' f = c^{e},   c = y^{q} g,  deg g = {dg}")

gc = sp.symbols(f"g0:{dg+1}")
A = sp.symbols("A")
g = sum(gc[i] * y**i for i in range(dg + 1))
cc = y**q * g
f = A * y**rho * g**e
resid = sp.expand(a * t * cc * sp.diff(f, y) - a * coef * sp.diff(cc, y) * f - cc**e)
quo = sp.expand(sp.factor(resid) / (y**(e * q) * g**(e - 1)))
print(f"  ansatz f = A y^{rho} g^{e} collapses the ODE to:  3A(y g' - {dg} g) = 1")
print(f"    sympy residual / (y^{e*q} g^{e-1}) = {sp.expand(quo)}")

# forced solution: g1 = g2 = 0, g3 resonant (free), g(-1)=0 forces g0=g3, monic.
g_sol = y**dg + 1
subs = {gc[i]: sp.Poly(g_sol, y).coeff_monomial(y**i) for i in range(dg + 1)}
A_sol = sp.solve(sp.expand(quo.subs(subs)).coeff(y, 0), A)[0]
C = y**q * g_sol
f = sp.expand(A_sol * y**rho * g_sol**e)
H2 = sp.factor(g_sol / (y + 1))
assert sp.expand(a * t * C * sp.diff(f, y) - a * coef * sp.diff(C, y) * f - C**e) == 0
print(f"  g = y^3 + 1 = (y+1)(y^2-y+1),   H2 = {H2}   (separable, avoids 0,-1)")
print(f"  C  = y^2 (y^3+1)            deg C = {sp.degree(C, y)}  (= a0)")
print(f"  f  = -(1/9) y^5 (y^3+1)^3   deg f = {sp.degree(f, y)}   A = {A_sol}")

# ===========================================================================
# 2. THE TOWER -- build the D-transform tower and DERIVE N (the decisive step).
#
#    Normalised C-series:  S = sum_k d_k u^(t-k),   d_t = 1, d_{t-1} = 0 (shift),
#    with the D-transform d_k = c_k * c^(a(t-k)-1)  (c := leading poly C).
#    P = C^a lives in S^a (linear window); the forcing/Phi lives in the C^b
#    tower S^b, at the (D~^b)_{-j} slice = coeff of u^M, M = b*t + j.
#
#    KEY LEMMA (slice-sum invariant).  For any monomial d_{k_1}..d_{k_b} of the
#    u^M coefficient of S^b we have sum_i (t - k_i) = M, hence the total
#    c-exponent is
#        sum_i (a(t-k_i) - 1) = a * (sum_i (t-k_i)) - b = a*M - b  =: clear,
#    depending ONLY on the (integer) slice u-power M -- never on the individual
#    k_i.  So Phi = c^clear * F_s = f * c^(clear - b), and N := clear - b.
# ===========================================================================
print("\n--- 2. D-transform tower and the DERIVED tower length N ---------------")
M = b * t + j
Kwin = 8
dvars = {}
for k in range(t - Kwin, t + 1):
    if k == t:
        dvars[k] = sp.Integer(1)
    elif k == t - 1:
        dvars[k] = sp.Integer(0)
    else:
        dvars[k] = sp.symbols(f"c_{t-k}")     # generic C-series coefficient c_k
S = sum(dvars[k] * c**(a * (t - k) - 1) * u**(t - k) for k in dvars)
Sb = sp.Poly(sp.expand(S**b), u)

# verify the slice-sum invariant on every reachable slice, incl. the forcing one
homog = {}
for MM in range(0, b * Kwin + 1):
    expr = sp.expand(Sb.coeff_monomial(u**MM))
    if expr == 0:
        continue
    exps = {c_exponent(tm) for tm in sp.Add.make_args(expr)}
    assert len(exps) == 1 and next(iter(exps)) == a * MM - b, (MM, exps)
    homog[MM] = next(iter(exps))

clear = homog[M]
N = clear - b
assert clear == a * M - b
print(f"  S^{b} u-slices are c-homogeneous with c-exponent = a*M - b  (verified {len(homog)} slices)")
print(f"  forcing slice:  M = b*t + j = {b}*{t} + {j} = {M}   (reached directly: {M in homog})")
print(f"  clearing exponent read off the tower:  clear = a*M - b = {clear}")
print(f"  DERIVED tower length:  N = clear - b = {N}")

N_formula = a * (t * (a + b) - (kappa + 1)) - 2 * b
print(f"  corner-144 N-formula a[t(a+b)-(kappa+1)]-2b = {N_formula}   (agrees: {N == N_formula})")

# the judgment-item-3 concern, made precise and dissolved
slice_index = sp.Rational(t) - b * t + s + sp.Rational(b - 1, a)
print(f"\n  [judgment 3 -> DERIVED]  the per-term slice index k = t-bt+s+(b-1)/a = {slice_index}")
print(f"     is NON-integral because (b-1)/a = {sp.Rational(b-1,a)} (unlike (72,108)/(108,144) where it is 1).")
print(f"     But Phi sits at the u-power SLICE M = {M} (integral), and clear = a*M-b is a")
print("     slice-SUM invariant, independent of that per-term index.  The non-integral")
print("     slice index therefore cannot move Phi or change N.  N = 98 is DERIVED, not")
print("     merely extrapolated -- PHI_75_125 judgment item 3 is upgraded.")

# ===========================================================================
# 3. Phi = f * C^N and its divisor signature.
# ===========================================================================
print("\n--- 3. Phi = f * C^98 and its signature -------------------------------")
Phi = sp.expand(f * C**N)
deg = sp.degree(Phi, y)
ordy = order(Phi)
m1 = multiplicity(Phi, y + 1)
cof = sp.cancel(Phi / (y**ordy * (y + 1)**m1))
cofdeg = sp.degree(cof, y)
assert sp.expand(Phi + sp.Rational(1, 9) * y**201 * (y**3 + 1)**101) == 0
print(f"  Phi = f * C^{N} = -(1/9) y^{ordy} (y^3+1)^{m1}")
print(f"      = -(1/9) y^{ordy} (y+1)^{m1} (y^2-y+1)^{m1}")
print(f"  SIGNATURE (deg, ord_y, mult_(y+1), cofactor) = ({deg}, {ordy}, {m1}, {cofdeg})")
print(f"  cofactor = -(1/9)(y^2-y+1)^{m1}  (residual H2 rides INSIDE C, not a new unit place)")
print(f"  target signature (504,201,101,202): {(deg,ordy,m1,cofdeg) == (504,201,101,202)}")

# ===========================================================================
# 4. Controls -- the same machinery on two landed checkpoints.
# ===========================================================================
print("\n--- 4. Controls (same tower+ODE machinery) ----------------------------")
def tower_N(a, b, t, kappa, Kwin):
    s = kappa + 1 - a * t
    M = b * t - s
    dvars = {k: (sp.Integer(1) if k == t else sp.Integer(0) if k == t - 1
                 else sp.symbols(f"z_{t-k}")) for k in range(t - Kwin, t + 1)}
    S = sum(dvars[k] * c**(a * (t - k) - 1) * u**(t - k) for k in dvars)
    Sb = sp.Poly(sp.expand(S**b), u)
    expr = sp.expand(Sb.coeff_monomial(u**M))
    exps = {c_exponent(tm) for tm in sp.Add.make_args(expr)}
    assert exps == {a * M - b}
    return M, a * M - b, (a * M - b) - b

for tag, (aa, bb, tt, kk, Cpoly, fpoly, want) in {
    "(108,144)": (3, 4, 4, 2, y**3 * (y**5 + 1),
                  -y**4 * (y**5 + 1)**2 / sp.Integer(15), (550, 205, 69, 276)),
    "(72,108)":  (2, 3, 4, 2, y**7 * (y + 1),
                  -y**8 * (y + 1)**2 * (2048 * y**4 - 512 * y**3 + 320 * y**2
                                        - 240 * y + 195) / sp.Integer(6630),
                  (238, 204, 30, 4)),
}.items():
    MM, cl, NN = tower_N(aa, bb, tt, kk, 7)
    Ph = sp.expand(fpoly * Cpoly**NN)
    o = order(Ph); mm = multiplicity(Ph, y + 1)
    sg = (sp.degree(Ph, y), o, mm, sp.degree(sp.cancel(Ph / (y**o * (y + 1)**mm)), y))
    print(f"  {tag}: tower M={MM} clear={cl} DERIVED N={NN}; Phi signature {sg}  "
          f"(want {want}: {sg == want})")

print("\n" + "=" * 78)
print("VERDICT: N = 98 CONFIRMED -- derived from the built D-transform tower.")
print("Phi = -(1/9) y^201 (y^3+1)^101 emerges with signature (504,201,101,202).")
print("=" * 78)

if __name__ == "__main__":
    pass
