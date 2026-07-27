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
  * F17 j=0 (case (66,99), corner (9,24) -> (11\\3,8), (a,b)=(2,3), t=3,
    r=0, gap=4): the forced pure ansatz FAILS at the top (resonance broken,
    exactly as in (72,108)) and the unique polynomial solution of the ODE is
        f = -(1/910) y^9 (y+1)^2 (243y^4 - 81y^3 + 54y^2 - 42y + 35)
    -- the pure-ansatz shape times a UNIT cofactor of degree exactly gap=4
    (u(0) != 0, u(-1) != 0), the exact analogue of the (72,108) quartic AND at
    the same gap value 4.
        Phi = f * C^20,  signature (195, 169, 22, 4).
    The r=0-amended law   deg = (e*degC - ordC + 1) + gap + N*degC,
    cofactor = gap   now has THREE exact points ((72,108) audited, F17 gap=4,
    F8 gap=2); with gap=0 it is the ordinary law, so ONE unified statement
    covers every known point:
        signature = ( (e*degC-ordC+1) + gap + N*degC,  (e-1)ordC+1 + N*ordC,
                      e + N,                           gap + r*(e+N) ).
  * Two structural mini-lemmas (checked over the whole GGV5 survey):
        gap = (ordC-1) - degC/t     (so gap = 0  <=>  degC = t(ordC-1)),
        residual index dg = degC - ordC  (g = y^dg + 1; the 10th-cyclotomic
        residual of (108,144) recurs at F14 because both have dg=5).
    UNTESTED regime, named honestly: gap>0 with r>0 (F7, F15, F16) -- the
    unified cofactor formula gap + r(e+N) is a CONJECTURE there, and PHI_F7.md
    REFUTES it on the ramified branch.

*** 2026-07-27 CHART REPAIR -- WHY THIS FILE'S SECOND POINT MOVED. ***

The gap>0 / r=0 probe used to be F1 j=0 (48,64) at the corner (4,12).  GGV5's
final chain corner (7\\4,3) was read as chart data (t, deg C, ord C) = (4,4,3),
giving dg=1, r=0, gap=1.  That dictionary is valid only on the retraction shape
b0 == ceil(b0/a0)*(a0-1), and (4,12) FAILS it (ceil(12/4) = 3, 3*3 = 9 != 12), so
polygon_reduction.final_corner_dictionary() now RAISES there.  The repaired
chart is (t, deg C, ord C) = (3,1,1): C is the MONOMIAL y, dg = 0, gap = -1/3.
So (4,12) is NOT in the gap>0 / r=0 regime at all -- the probe had no subject.

REPLACEMENTS, on corners that DO retract ((6,15) and (9,24)):
  * F17 (66,99), gap = 4, r = 0  -- the same gap as the audited (72,108);
  * F8  (63,147), gap = 2, r = 0 -- a second gap value, so the regime is not
    tested at one gap only.
Both are derived in STEP 2 below.  The regime therefore ends up with THREE exact
points rather than the two it had, and none of them rests on a refused corner.

F14 itself is UNAFFECTED: its corner (9,24) retracts (24 = 3*(9-1)), so
(t, deg C, ord C) = (3,9,4) is exactly what the guard returns.

INDEPENDENT TARGET: every ord_y(Phi) here is checked against
ord_y(Phi) = a*q*M - H, PROVED in BRIDGE_GENERALITY.md and computed by neither
this file nor its checker.

Sources: GGV5 family tables as transcribed in phi_corner4.py (Diophantine
re-checked here); chart data through polygon_reduction.corner_chart_data; method
template PHI_75_125.md / PHI_CORNER4.md; bridge BRIDGE_GENERALITY.md.  The
independent PASS/FAIL checker is phi_f14_verify.py.  Everything is exact sympy.
"""
import sympy as sp
from fractions import Fraction
from math import gcd

import polygon_reduction as pr

y = sp.symbols("y")

# ---------------------------------------------------------------------------
# 0. Corner rows (GGV5 v11<=35 tables; same transcription as phi_corner4.py)
#    (name, A0, p, l_final, b_final, k, (m0,dm), (n0,dn))
#
# l_final / b_final are CHAIN data (they enter the Diophantine identity).  The
# CHART data (t, kappa, deg C, ord C) is DERIVED from A0 through the guard.
# ---------------------------------------------------------------------------
F14 = ("F14", (9, 24), 7, 3, 4, 1, (2, 1), (7, 4))
F17 = ("F17", (9, 24), 11, 3, 8, 1, (2, 5), (3, 8))   # gap=4 r=0, RETRACTS
F8  = ("F8",  (6, 15), 8, 3, 5, 1, (3, 2), (7, 5))    # gap=2 r=0, RETRACTS
F1  = ("F1",  (4, 12), 7, 4, 3, 1, (3, 2), (4, 3))    # REFUSED corner (retired)

# SUPERSEDED, kept LABELLED so the retirement stays falsifiable (see the docstring
# and phi_f14_verify.py sec. H).  (t_stale, degC_stale, ordC_stale, N, sig, f)
SUPERSEDED = {
    "F1": (4, 4, 3, 67, (275, 205, 69, 1),
           sp.Rational(1, 15) * y**4 * (y + 1)**2 * (4 * y - 1)),
}


def bridge_ord(a, b, t, kappa, ordC):
    """ord_y(Phi) = a*q*M - H.  PROVED in BRIDGE_GENERALITY.md, not fitted here."""
    s = a + b
    return a * ordC * (t * s - (kappa + 1)) - (ordC * s - 1)


def corner_params(row):
    name, A0, p, l, q, k, (m0, dm), (n0, dn) = row
    j = 0
    while gcd(m0 + dm * j, n0 + dn * j) != 1:
        j += 1
    m, n = m0 + dm * j, n0 + dn * j
    assert (m + n) * q * k - n * (q * l - p) == k, f"{name}: Diophantine failed"
    v11 = A0[0] + A0[1]
    a, b = sorted((m, n))
    # ---- CHART data through the guard; RAISES off the retraction shape --------
    cd = pr.corner_chart_data(A0[0], A0[1], l_final=l, b_final=q,
                              who="phi_f14 %s" % name)
    t, kappa, degC, ordC = cd["t"], cd["kappa"], cd["deg_C"], cd["ord_C"]
    e, dg = b - a + 1, degC - ordC
    r = dg - 1
    coef = t * (b - a) + kappa + 1
    rho = (e - 1) * ordC + 1
    N = a * (t * (a + b) - (kappa + 1)) - 2 * b
    res = Fraction(coef * degC, t)
    gap = res - (e * degC - ordC + 1)
    return dict(name=name, j=j, m=m, n=n, degs=(v11 * m, v11 * n), a=a, b=b,
                t=t, kappa=kappa, a0=degC, q=ordC, degC=degC, ordC=ordC,
                e=e, r=r, dg=dg, coef=coef, rho=rho, N=N, gap=gap,
                retracts=cd["retraction"], A0=A0, l_final=l, b_final=q,
                bridge=bridge_ord(a, b, t, kappa, ordC))

# ---------------------------------------------------------------------------
# 2026-07-26 corner-law generalizations (PASSPORT_75_125_REPAIR.md).  See the
# patch note in the repair doc; both are forced by the (5,20) corner, which is
# the first case with C a MONOMIAL (deg g = 0) and with a non-integral gap.
# ---------------------------------------------------------------------------
def gap_effective(gap):
    """The resonance gap AS AN EXTRA UNIT-FACTOR DEGREE.

    gap = (q-1) - a0/t is the resonant degree minus the pure-ansatz degree of f.
    It contributes an extra factor only when it is a POSITIVE INTEGER ((72,108):
    gap = 4).  Negative or non-integral gap means the resonance does not sit at
    or above the ansatz degree, no extra factor appears, and the pure ansatz is
    exact -- effective gap 0.  At the repaired (5,20) corner gap = -1/4 and the
    independent ODE solve confirms deg f = rho exactly.
    """
    from fractions import Fraction as _F
    g = _F(gap)
    return int(g) if (g.denominator == 1 and g > 0) else 0


def mult_and_cofactor(e, N, degC, ordC, gap):
    """(mult_(y+1), cofactor_deg) with the residual-free branch dg = degC-ordC = 0.

    2026-07-27: the third and fourth arguments are deg_y(C) and ord_y(C), NOT the
    corner coordinate a0 and the chain datum b_final -- those coincide only under
    GGV5's final-corner dictionary, i.e. only on the retraction shape.

    dg > 0 : g = y^dg + 1 contributes (y+1)^(e+N) and a residual H2 = g/(y+1) of
             degree dg-1, so mult = e+N and cof = gap + (dg-1)*(e+N).
    dg == 0: C is a MONOMIAL.  There is NO g, hence no (y+1) place and no
             residual: mult = 0 and cof = gap.
    """
    ge = gap_effective(gap)
    dg = degC - ordC
    if dg == 0:
        return 0, ge
    return e + N, ge + (dg - 1) * (e + N)


def law_sig(P):
    """Unified law (gap term included; reduces to the old law when gap=0).

    2026-07-26: gap enters through gap_effective (an extra unit factor only when
    it is a POSITIVE INTEGER) and mult/cofactor through the residual-free branch.
    """
    ge = gap_effective(P["gap"])
    deg = (P["e"] * P["a0"] - P["q"] + 1) + ge + P["N"] * P["a0"]
    ordy = P["rho"] + P["N"] * P["q"]
    mult, cof = mult_and_cofactor(P["e"], P["N"], P["a0"], P["q"], P["gap"])
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
# 2. The gap>0 / r=0 probe -- solve the ODE with FULLY GENERIC f.
#    2026-07-27: run at F17 (gap=4) and F8 (gap=2), both on RETRACTING corners.
#    F1 (48,64) used to be the single point here; its corner (4,12) is refused.
# ---------------------------------------------------------------------------
def probe_r0(row, label):
    P1 = corner_params(row)
    assert P1["retracts"], (P1["name"], "gap>0/r=0 probe needs a real residual")
    assert P1["r"] == 0 and P1["gap"] > 0 and P1["dg"] == 1, P1
    print("\n" + "=" * 96)
    print(f"STEP 2{label} -- {P1['name']} j={P1['j']}, case {P1['degs']}: "
          f"(a,b)=({P1['a']},{P1['b']}) t={P1['t']} kappa={P1['kappa']} "
          f"deg C={P1['degC']} ord C={P1['ordC']} e={P1['e']} r={P1['r']} "
          f"rho={P1['rho']} N={P1['N']} gap={P1['gap']}  (r=0 regime; corner "
          f"{P1['A0']} RETRACTS)")
    print("=" * 96)
    a, t, coef, q, e, rho, N = (P1[k] for k in
                                ("a", "t", "coef", "q", "e", "rho", "N"))
    res_deg = P1["rho"] + P1["e"] * P1["dg"] + int(P1["gap"])  # = e*degC-ordC+1+gap
    print(f"  ODE:  {a*t} c f' - {a*coef} c' f = c^{e},   c = y^{q} (y+1)")
    print(f"  pure-ansatz degree {P1['rho'] + P1['e']*P1['dg']} < resonant degree "
          f"{res_deg}  (gap {P1['gap']}) -- resonance broken, exactly as (72,108)")
    c1 = y**q * (y + 1)
    fc = sp.symbols(f"f{label}0:{res_deg + 3}")   # allow 2 beyond resonant degree
    f_gen = sum(fc[i] * y**i for i in range(res_deg + 3))
    resid1 = sp.expand(a * t * c1 * sp.diff(f_gen, y)
                       - a * coef * sp.diff(c1, y) * f_gen - c1**e)
    sols = sp.solve(sp.Poly(resid1, y).all_coeffs(), fc, dict=True)
    assert len(sols) == 1, "polynomial solution not unique"
    f1_sol = sp.expand(f_gen.subs(sols[0]))
    assert sp.degree(f1_sol, y) == res_deg
    u = sp.expand(sp.cancel(f1_sol / (y**rho * (y + 1)**e)))
    assert sp.Poly(u, y).total_degree() == int(P1["gap"])
    print(f"  UNIQUE polynomial solution (generic solve, deg allowed up to "
          f"{res_deg+2}):")
    print(f"    f = {sp.factor(f1_sol)}")
    print(f"    = (pure-ansatz shape) * u,   deg u = {sp.degree(u, y)} = gap;"
          f" u(0) = {u.subs(y,0)}, u(-1) = {u.subs(y,-1)} -- both nonzero: UNIT")
    Phi1 = f1_sol * c1**N
    sig1 = signature(Phi1)
    fit1 = law_sig(P1)
    v1 = "MATCHES" if sig1 == fit1 else "DIFFERS"
    print(f"  Phi = f C^{N},  SIGNATURE = {sig1}")
    print(f"  r=0-amended law prediction  = {fit1}   ==> {v1}")
    assert sig1[1] == P1["bridge"], (P1["name"], sig1, P1["bridge"])
    print(f"  bridge identity a*q*M - H = {P1['bridge']}  ==> ord_y AGREES "
          f"[INDEPENDENT, PROVED: BRIDGE_GENERALITY.md]")
    return P1, sig1, fit1, v1

P17, sig17, fit17, v17 = probe_r0(F17, "a")
P8, sig8, fit8, v8 = probe_r0(F8, "b")
P_GGHV = dict(e=2, a0=8, q=7, degC=8, ordC=7, gap=Fraction(4), N=28, rho=8, r=0)
print(f"\n  (72,108) control: stored (238,204,30,4) vs amended law "
      f"{law_sig(P_GGHV)} -- the audited case obeys the same amended law")
print(f"  gap coverage in this regime is now {{2, 4}} on two DIFFERENT corners "
      f"((6,15) and (9,24)) plus the audited gap=4 at (8,28): three points, none")
print(f"  of them on a corner the retraction guard refuses.")

# ---- RETIRED: F1 (48,64).  Right ODE, wrong corner. -------------------------
print("\n" + "=" * 96)
print("RETIRED 2026-07-27 -- F1 (48,64): corner (4,12) is guard-REFUSED")
print("=" * 96)
P1 = corner_params(F1)                       # the REPAIRED chart, via the guard
st, sdC, soC, sN, ssig, sf = SUPERSEDED["F1"]
print(f"  GGV5's chain row gives (l_final,b_final) = ({F1[3]},{F1[4]}); the "
      f"pre-repair file used that AS chart data,")
print(f"  i.e. (t,deg C,ord C) = ({st},{sdC},{soC}) -> dg=1, r=0, gap=1, N={sN}, "
      f"sig={ssig}, f = (1/15) y^4 (y+1)^2 (4y-1).")
print(f"  Corner (4,12) is REFUSED: ceil(12/4) = {pr.chart_exponent(4,12)} and "
      f"{pr.chart_exponent(4,12)}*(4-1) = {pr.chart_exponent(4,12)*3} != 12.")
print(f"  Repaired chart: (t,deg C,ord C) = ({P1['t']},{P1['degC']},{P1['ordC']}), "
      f"dg = {P1['dg']}, gap = {P1['gap']} < 0.")
print(f"  So (4,12) is NOT in the gap>0/r=0 regime: C = y is a MONOMIAL, there is "
      f"no residual (y+1)")
print(f"  and no unit cofactor.  Repaired N = {P1['N']}, ord_y(Phi) = "
      f"{P1['bridge']} (bridge identity) vs {ssig[1]} before.")
# the retired f still solves the ODE it was computed for: retired as a claim
# about the corner, not withdrawn as arithmetic.
_a, _b = P1["a"], P1["b"]
_e = _b - _a + 1
_coef = st * (_b - _a) + (st - 2) + 1
_c = y**soC * (y + 1)
assert sp.expand(_a * st * _c * sp.diff(sf, y)
                 - _a * _coef * sp.diff(_c, y) * sf - _c**_e) == 0
_crep = y**P1["ordC"]
assert sp.expand(_a * P1["t"] * _crep * sp.diff(sf, y)
                 - _a * P1["coef"] * sp.diff(_crep, y) * sf - _crep**_e) != 0
print(f"  The retired f DOES still solve the ODE at the stale parameters "
      f"({st},{st-2},{soC},1) and does NOT")
print(f"  solve the repaired corner's ODE -- retired as a claim about (4,12), not "
      f"as arithmetic.")
_frep = sp.Rational(1, _a) * y**_e
assert sp.expand(_a * P1["t"] * _crep * sp.diff(_frep, y)
                 - _a * P1["coef"] * sp.diff(_crep, y) * _frep - _crep**_e) == 0
print(f"  The repaired corner's unique solution is f = (1/{_a}) y^{_e}, so "
      f"Phi = (1/{_a}) y^{_e + P1['N']}: signature "
      f"({_e + P1['N']}, {_e + P1['N']}, 0, 0).")
sig1, fit1, v1 = signature(_frep * _crep**P1["N"]), law_sig(P1), None
v1 = "MATCHES" if sig1 == fit1 else "DIFFERS"
print(f"  law prediction {fit1}  ==> {v1}")

# ---------------------------------------------------------------------------
# 3. Structural mini-lemmas over the whole length-1 survey
# ---------------------------------------------------------------------------
print("\n" + "=" * 96)
print("STEP 3 -- mini-lemmas:  gap = (ordC-1) - degC/t   and   dg = degC - ordC")
print("         (2026-07-27: in terms of deg_y(C) and ord_y(C) from the guard, NOT")
print("          the corner coordinate a0 and the chain datum b_final)")
print("=" * 96)
ROWS = [
    ("F1", (4,12), 7,4,3,1,(3,2),(4,3)),  ("F2", (5,20), 7,5,2,1,(2,1),(3,2)),
    ("F3", (5,20), 8,5,3,1,(3,4),(2,3)),  ("F4", (5,20), 8,5,3,2,(3,2),(16,12)),
    ("F5", (5,20), 9,5,4,1,(9,7),(5,4)),  ("F6", (5,20), 9,5,4,2,(7,6),(18,16)),
    ("F7", (6,15), 7,3,4,1,(2,1),(7,4)),  ("F8", (6,15), 8,3,5,1,(3,2),(7,5)),
    ("F9", (7,21), 11,7,2,1,(2,1),(3,2)), ("F10",(7,21), 13,7,3,1,(7,5),(4,3)),
    ("F11",(7,21), 13,7,3,2,(2,1),(5,3)), ("F14",(9,24), 7,3,4,1,(2,1),(7,4)),
    ("F15",(9,24), 8,3,5,1,(3,2),(7,5)),  ("F16",(9,24), 10,3,7,1,(3,4),(5,7)),
    ("F17",(9,24), 11,3,8,1,(2,5),(3,8)),
]
_refused, _retr, _bridge_ok = [], [], True
for row in ROWS:
    Pr = corner_params(row)
    lem = Fraction(Pr["ordC"] - 1) - Fraction(Pr["degC"], Pr["t"])
    assert Pr["gap"] == lem, row[0]
    assert Pr["dg"] == Pr["degC"] - Pr["ordC"], row[0]
    (_retr if Pr["retracts"] else _refused).append(row[0])
    # the law's ord component must equal the PROVED bridge value on EVERY row
    if law_sig(Pr)[1] != Pr["bridge"]:
        _bridge_ok = False
assert _bridge_ok, "the law's ord_y disagrees with the bridge identity somewhere"
print("  gap = (ordC-1) - degC/t verified on all 15 standard-chart families;")
print("  gap = 0  <=>  degC = t(ordC-1).  dg = degC-ordC throughout, so the")
print("  10th-cyclotomic residual recurs at every dg=5 corner: (108,144) and F14.")
print(f"  ord_y(Phi) == a*q*M - H on all 15 rows (bridge identity, PROVED).")
print(f"  RETRACTING rows (dictionary VALID, dg >= 1): {', '.join(_retr)}")
print(f"  GUARD-REFUSED rows (C = y, dg = 0, gap < 0): {', '.join(_refused)}")
print("  NOTE: pre-repair this lemma was stated with a0 and b_final in place of")
print("  deg C and ord C, and it 'held' on all 15 rows because BOTH sides used the")
print("  same substitution.  It is a real lemma only in the guarded variables.")
print("  UNTESTED (named): gap>0 with r>0 (F7, F15, F16 -- the retracting members")
print("  of that regime) -- cofactor formula gap + r(e+N) is a conjecture there,")
print("  and PHI_F7.md REFUTES it on the ramified branch.  F3 and F10 are NOT in")
print("  the regime after the repair: their corners are refused and dg = 0.")

# ---------------------------------------------------------------------------
# 4. Verdict
# ---------------------------------------------------------------------------
print("\n" + "=" * 96)
print("VERDICT")
print("=" * 96)
print(f"  F14 (66,231), t=3, gap=0     : derived {sig14}  law {fit14}  -> {v14}")
print(f"  F17 (66,99),  t=3, gap=4 r=0 : derived {sig17}  law {fit17}  -> {v17}")
print(f"  F8  (63,147), t=3, gap=2 r=0 : derived {sig8}   law {fit8}   -> {v8}")
print(f"  F1  (48,64)  RETIRED (corner (4,12) refused): repaired {sig1}  law "
      f"{fit1}  -> {v1}")
print("""  The unified law   ( pure+gap+N*degC, rho+N*ordC, e+N, gap+r(e+N) )
  reproduces every exact point on a corner the retraction guard ACCEPTS:
    (72,108) r=0 gap=4* audited, (108,144), F14 (66,231), F17 (66,99) r=0 gap=4,
    F8 (63,147) r=0 gap=2,
  and, in the monomial regime the guard forces at refused corners, the repaired
  (75,125), (50,75), (56,84) and (48,64) -- all of shape (D, D, 0, 0).

  2026-07-27 CHART REPAIR, evidence ledger for a reader of the public tree:
    * t in {3,4,5,7} was an artefact of reading t = l_final.  The DERIVED chart
      exponents of GGV5's v11<=35 length-1 tables are t in {3,4} only
      (phi_corner4.py STEP 1b).  This file's census sentence is withdrawn.
    * the gap>0 / r=0 regime is UNCHANGED as a law but its evidence moved: it had
      the audited (72,108) plus F1 (48,64); F1's corner (4,12) is refused, so F1
      is RETIRED and REPLACED by TWO fresh points on retracting corners,
      F17 (gap=4) and F8 (gap=2).  Net: two evidence points became three, and
      the regime is now tested at two distinct gap values on two distinct
      corners rather than at one corner.
    * (72,108), (108,144) and F14 are at (8,28) and (9,24), both of which
      retract, so the audited home case, the corner-144 template and this file's
      original sixth point are untouched.
    * every ord_y above is checked against the PROVED bridge identity
      ord_y(Phi) = a*q*M - H (BRIDGE_GENERALITY.md), which neither this file nor
      its checker computes.  Before the repair the law was validated only
      against targets produced by the same superseded chart dictionary.""")
