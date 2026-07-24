#!/usr/bin/env python3
"""phi_f14.py  (NEW; read-only over all existing artifacts)

SIXTH corner point: F14 at t=3  --  and the first fresh derivation in the
gap>0 / r=0 resonance regime (F1), testing whether the (72,108) extra-unit
story generalizes.

Result in brief:
  * F14 j=0 (case (66,231), corner (9,24) -> (7\\3,4), (a,b)=(2,7), t=3):
        Phi = -(1/10) y^165 (y^5+1)^42,  signature (375, 165, 42, 168)
    -- the five-parameter corner law's parameter-free prediction, exactly.
    Coverage after this point: t in {3,4,5,7}, e in {2,3,6}, a0 in {5,7,8,9},
    q in {2,3,4}.  (F14 is the first e=6 and first q=4 point.)
  * F1 j=0 (case (48,64), corner (4,12) -> (7\\4,3), (a,b)=(3,4), t=4,
    r=0, gap=1): the forced pure ansatz FAILS at the top (resonance broken,
    exactly as in (72,108)) and the unique polynomial solution of the ODE is
        f = (1/15) y^4 (y+1)^2 (4y-1)
    -- the pure-ansatz shape times a UNIT cofactor of degree exactly gap=1
    (u(0) != 0, u(-1) != 0), the degree-1 analogue of the (72,108) quartic.
        Phi = (1/15) y^205 (y+1)^69 (4y-1),  signature (275, 205, 69, 1).
    The r=0-amended law   deg = (e*a0-q+1) + gap + N*a0,  cofactor = gap
    now has TWO exact points ((72,108) audited + F1 fresh); with gap=0 it is
    the ordinary law, so ONE unified statement covers all seven known points:
        signature = ( (e*a0-q+1) + gap + N*a0,  (e-1)q+1 + N*q,
                      e + N,                    gap + r*(e+N) ).
  * Two structural mini-lemmas (checked over the whole GGV5 survey):
        gap = (q-1) - a0/t          (so gap = 0  <=>  a0 = t(q-1)),
        residual index dg = a0 - q  (g = y^dg + 1; the 10th-cyclotomic
        residual of (108,144)/F9 recurs at F14 because all three have dg=5).
    UNTESTED regime, named honestly: gap>0 with r>0 (F3, F7, F10, F15, F16)
    -- the unified cofactor formula gap + r(e+N) is a CONJECTURE there.

Sources: GGV5 family tables as transcribed in phi_corner4.py (Diophantine
re-checked here); method template PHI_75_125.md / PHI_CORNER4.md.  The
independent PASS/FAIL checker is phi_f14_verify.py.  Everything is exact sympy.
"""
import sympy as sp
from fractions import Fraction
from math import gcd

y = sp.symbols("y")

# ---------------------------------------------------------------------------
# 0. Corner rows (GGV5 v11<=35 tables; same transcription as phi_corner4.py)
#    (name, A0, p, l, q, k, (m0,dm), (n0,dn))
# ---------------------------------------------------------------------------
F14 = ("F14", (9, 24), 7, 3, 4, 1, (2, 1), (7, 4))
F1  = ("F1",  (4, 12), 7, 4, 3, 1, (3, 2), (4, 3))

def corner_params(row):
    name, A0, p, l, q, k, (m0, dm), (n0, dn) = row
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
    gap = res - (e * a0 - q + 1)
    return dict(name=name, j=j, m=m, n=n, degs=(v11 * m, v11 * n), a=a, b=b,
                t=t, kappa=kappa, a0=a0, q=q, e=e, r=r, dg=dg, coef=coef,
                rho=rho, N=N, gap=gap)

def law_sig(P):
    """Unified law (gap term included; reduces to the old law when gap=0)."""
    deg = (P["e"] * P["a0"] - P["q"] + 1) + P["gap"] + P["N"] * P["a0"]
    ordy = P["rho"] + P["N"] * P["q"]
    mult = P["e"] + P["N"]
    cof = P["gap"] + P["r"] * (P["e"] + P["N"])
    return (int(deg), int(ordy), int(mult), int(cof))

def signature(Phi):
    P = sp.Poly(sp.expand(Phi), y)
    deg = P.degree()
    ordy = min(mm[0] for mm in P.monoms())
    mult, Q = 0, P
    while True:
        q2, r2 = sp.div(Q, sp.Poly(y + 1, y))
        if not r2.is_zero:
            break
        Q, mult = q2, mult + 1
    return (deg, ordy, mult, deg - ordy - mult)

# ---------------------------------------------------------------------------
# 1. F14: generic-g collapse (gap = 0, so the pure ansatz is exact)
# ---------------------------------------------------------------------------
P = corner_params(F14)
print("=" * 96)
print(f"STEP 1 -- F14 j={P['j']}, case {P['degs']}: (a,b)=({P['a']},{P['b']}) "
      f"t={P['t']} kappa={P['kappa']} a0={P['a0']} q={P['q']} e={P['e']} "
      f"r={P['r']} dg={P['dg']} rho={P['rho']} N={P['N']} gap={P['gap']}")
print("=" * 96)
a, t, coef, q, e, rho, dg, N = (P[k] for k in
                                ("a", "t", "coef", "q", "e", "rho", "dg", "N"))
print(f"  forcing ODE:  {a*t} c f' - {a*coef} c' f = c^{e},   c = y^{q} g,  deg g = {dg}")

gc = sp.symbols(f"g0:{dg+1}")
A = sp.symbols("A")
g = sum(gc[i] * y**i for i in range(dg + 1))
c = y**q * g
f = A * y**rho * g**e
resid = sp.expand(a * t * c * sp.diff(f, y) - a * coef * sp.diff(c, y) * f - c**e)
bracket = sp.expand(A * a * ((t * rho - coef * q) * g + (t * e - coef) * y * sp.diff(g, y)) - 1)
assert sp.expand(resid - y**(e * q) * g**e * bracket) == 0
print(f"  resid = y^{e*q} g^{e} * [A*a*(({t*rho - coef*q}) g + ({t*e - coef}) y g') - 1]")
for i in range(1, dg):
    mult_i = (t * rho - coef * q) + i * (t * e - coef)
    assert mult_i != 0
    print(f"    y^{i}: multiplier {mult_i} != 0  =>  g_{i} = 0")
assert (t * rho - coef * q) + dg * (t * e - coef) == 0
print(f"    y^{dg}: multiplier 0  =>  top coefficient RESONANT (free; monic normalization)")
g_sol = y**dg + 1                       # monic + g(-1)=0 (mult_(y+1) Phi > 0)
A_sol = sp.Rational(1, a * (t * rho - coef * q))
c_sol = y**q * g_sol
f_sol = sp.expand(A_sol * y**rho * g_sol**e)
assert sp.expand(a * t * c_sol * sp.diff(f_sol, y)
                 - a * coef * sp.diff(c_sol, y) * f_sol - c_sol**e) == 0
print(f"  g = y^{dg}+1,  H = {sp.factor(g_sol / (y+1))},  A = {A_sol}")
print(f"  f = {A_sol} y^{rho} (y^{dg}+1)^{e}   (deg {sp.degree(f_sol, y)})")

Phi14 = f_sol * c_sol**N
sig14 = signature(Phi14)
fit14 = law_sig(P)
v14 = "MATCHES" if sig14 == fit14 else "DIFFERS"
print(f"  Phi = f C^{N} = {A_sol} y^{rho + N*q} (y^{dg}+1)^{e+N}")
print(f"  SIGNATURE (deg, ord_y, mult_(y+1), cofactor) = {sig14}")
print(f"  law prediction (parameter-free)              = {fit14}   ==> {v14}")

# ---------------------------------------------------------------------------
# 2. F1: the gap>0 / r=0 probe -- solve the ODE with FULLY GENERIC f
# ---------------------------------------------------------------------------
P1 = corner_params(F1)
print("\n" + "=" * 96)
print(f"STEP 2 -- F1 j={P1['j']}, case {P1['degs']}: (a,b)=({P1['a']},{P1['b']}) "
      f"t={P1['t']} kappa={P1['kappa']} a0={P1['a0']} q={P1['q']} e={P1['e']} "
      f"r={P1['r']} rho={P1['rho']} N={P1['N']} gap={P1['gap']}  (r=0 regime)")
print("=" * 96)
a, t, coef, q, e, rho, N = (P1[k] for k in ("a", "t", "coef", "q", "e", "rho", "N"))
res_deg = P1["rho"] + P1["e"] * P1["dg"] + int(P1["gap"])   # = e*a0-q+1 + gap
print(f"  ODE:  {a*t} c f' - {a*coef} c' f = c^{e},   c = y^{q} (y+1)")
print(f"  pure-ansatz degree {P1['rho'] + P1['e']*P1['dg']} < resonant degree {res_deg}"
      f"  (gap {P1['gap']}) -- resonance broken, exactly as (72,108)")
c1 = y**q * (y + 1)
fc = sp.symbols(f"f0:{res_deg + 3}")            # allow 2 beyond resonant degree
f_gen = sum(fc[i] * y**i for i in range(res_deg + 3))
resid1 = sp.expand(a * t * c1 * sp.diff(f_gen, y)
                   - a * coef * sp.diff(c1, y) * f_gen - c1**e)
sols = sp.solve(sp.Poly(resid1, y).all_coeffs(), fc, dict=True)
assert len(sols) == 1, "polynomial solution not unique"
f1_sol = sp.expand(f_gen.subs(sols[0]))
assert sp.degree(f1_sol, y) == res_deg
u = sp.cancel(f1_sol / (y**rho * (y + 1)**e))
assert sp.denom(sp.together(u)) == 1 or sp.Poly(sp.together(u), y)  # polynomial
u = sp.expand(u)
print(f"  UNIQUE polynomial solution (generic solve, deg allowed up to {res_deg+2}):")
print(f"    f = {sp.factor(f1_sol)}")
print(f"    = (pure-ansatz shape) * u,   u = {u}   (deg u = {sp.degree(u, y)} = gap;"
      f" u(0) = {u.subs(y,0)}, u(-1) = {u.subs(y,-1)} -- both nonzero: UNIT)")

Phi1 = f1_sol * c1**N
sig1 = signature(Phi1)
fit1 = law_sig(P1)
v1 = "MATCHES" if sig1 == fit1 else "DIFFERS"
print(f"  Phi = f C^{N},  SIGNATURE = {sig1}")
print(f"  r=0-amended law prediction  = {fit1}   ==> {v1}")
P_GGHV = dict(e=2, a0=8, q=7, gap=Fraction(4), N=28, rho=8, r=0)
print(f"  (72,108) control: stored (238,204,30,4) vs amended law {law_sig(P_GGHV)} "
      f"-- the audited case obeys the same amended law")

# ---------------------------------------------------------------------------
# 3. Structural mini-lemmas over the whole length-1 survey
# ---------------------------------------------------------------------------
print("\n" + "=" * 96)
print("STEP 3 -- mini-lemmas:  gap = (q-1) - a0/t   and   residual index dg = a0-q")
print("=" * 96)
ROWS = [
    ("F1", (4,12), 7,4,3,1,(3,2),(4,3)),  ("F2", (5,20), 7,5,2,1,(2,1),(3,2)),
    ("F3", (5,20), 8,5,3,1,(3,4),(2,3)),  ("F4", (5,20), 8,5,3,2,(3,2),(16,12)),
    ("F5", (5,20), 9,5,4,1,(9,7),(5,4)),  ("F6", (5,20), 9,5,4,2,(4,3),(10,8)),
    ("F7", (6,15), 7,3,4,1,(2,1),(7,4)),  ("F8", (6,15), 8,3,5,1,(3,2),(7,5)),
    ("F9", (7,21), 11,7,2,1,(2,1),(3,2)), ("F10",(7,21), 13,7,3,1,(7,5),(4,3)),
    ("F11",(7,21), 13,7,3,2,(2,1),(5,3)), ("F14",(9,24), 7,3,4,1,(2,1),(7,4)),
    ("F15",(9,24), 8,3,5,1,(3,2),(7,5)),  ("F16",(9,24), 10,3,7,1,(3,4),(5,7)),
    ("F17",(9,24), 11,3,8,1,(2,5),(3,8)),
]
for row in ROWS:
    Pr = corner_params((row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7]))
    lem = Fraction(Pr["q"] - 1) - Fraction(Pr["a0"], Pr["t"])
    assert Pr["gap"] == lem, row[0]
print("  gap = (q-1) - a0/t verified on all 15 standard-chart families;")
print("  gap = 0  <=>  a0 = t(q-1).  dg = a0-q throughout (g = y^(a0-q)+1), so the")
print("  10th-cyclotomic residual recurs at every dg=5 corner: (108,144), F9, F14.")
print("  UNTESTED (named): gap>0 with r>0 (F3,F7,F10,F15,F16) -- cofactor formula")
print("  gap + r(e+N) is a conjecture in that regime until one is derived.")

# ---------------------------------------------------------------------------
# 4. Verdict
# ---------------------------------------------------------------------------
print("\n" + "=" * 96)
print("VERDICT")
print("=" * 96)
print(f"  sixth point  F14 (66,231), t=3 : derived {sig14}  law {fit14}  -> {v14}")
print(f"  seventh      F1  (48,64),  r=0 : derived {sig1}  law {fit1}  -> {v1}")
print("""  The unified law   ( pure+gap+N*a0, rho+N*q, e+N, gap+r(e+N) )   now
  reproduces SEVEN exact points at t in {3,4,5,7}: (72,108) r=0*, (108,144),
  (75,125), (56,84), (50,75), (66,231), (48,64) -- the last a fresh derivation
  in the resonance-gap regime, confirming the (72,108) extra-unit story
  generalizes (unit cofactor of degree exactly gap).""")
