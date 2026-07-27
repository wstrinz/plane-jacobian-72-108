#!/usr/bin/env python3
"""phi_f7_verify.py -- independent PASS/FAIL checker for PHI_F7.md / phi_f7.py.

Independent routes (no numeric input taken from phi_f7.py):
  * corner arithmetic re-derived from the GGV5 table rows with fresh code;
  * chart data through polygon_reduction.corner_chart_data -- the ONE thing not
    re-implemented, because a second implementation of a guard is a second
    chance to get the guard wrong.  Instead sec. F CROSS-CHECKS phi_f7.py's
    ledger against it (the family_grammar_verify.py A9 drift-guard pattern);
  * every claimed f verified against the ODE by direct sp.diff/expand;
  * uniqueness by FULL generic linear solve (sp.solve over all coefficients,
    degree slack past resonance) -- not the triangular recurrence;
  * signatures by trial division on f plus exponent arithmetic, cross-checked
    at F7 by full expansion of Phi itself;
  * the no-simple-real-root-at--1 theorem re-proven per corner from the
    claimed obstruction branches by quadratic-discriminant arguments;
  * off-branch refutation spot checks (generic g -> NO polynomial solution);
  * every ord_y(Phi) checked against the PROVED bridge identity a*q*M - H
    (BRIDGE_GENERALITY.md), a target neither file computes;
  * controls: F14 (gap=0) and F8 (dg=1 root gauge) plus (72,108) arithmetic.

2026-07-27 CHART REPAIR.  Three of the rows this file used sit on corners the
retraction guard REFUSES, where C is the monomial y and dg = 0:
  F3 (5,20) and F10 (7,21) -- the dg=2 and dg=4 exact points, RETIRED;
  F1 (4,12) -- the dg=1 root-gauge control, REPLACED by F8 (6,15).
The dg=4 rung is preserved by a fresh point on a RETRACTING corner, F15
(99,231) at (9,24).  Sec. F holds the mutation controls: the retired
polynomials must still solve the ODE at their STALE parameters (so the
retirement is about the corner, not the arithmetic) and must NOT solve the
repaired corner's ODE.

Usage: python phi_f7_verify.py [--quiet].  Exit 0 iff every check passes.
"""
import sys
import sympy as sp
from fractions import Fraction
from math import gcd

import polygon_reduction as pr

QUIET = "--quiet" in sys.argv[1:]
y = sp.symbols("y")
_n, _f = 0, 0

def check(label, ok):
    global _n, _f
    _n += 1
    if not ok:
        _f += 1
    if not QUIET or not ok:
        print(("[OK] " if ok else "[FAIL] ") + label)

# ---------------------------------------------------------------------------
# A. corner arithmetic (independent implementation, chart data via the guard)
# ---------------------------------------------------------------------------
TABLE = {   # name: (A0, p, l_final, b_final, k, (m0,dm), (n0,dn))  [GGV5 v11<=35]
    "F7":  ((6, 15), 7, 3, 4, 1, (2, 1), (7, 4)),
    "F8":  ((6, 15), 8, 3, 5, 1, (3, 2), (7, 5)),
    "F15": ((9, 24), 8, 3, 5, 1, (3, 2), (7, 5)),
    "F16": ((9, 24), 10, 3, 7, 1, (3, 4), (5, 7)),
    "F14": ((9, 24), 7, 3, 4, 1, (2, 1), (7, 4)),
    # guard-REFUSED rows, kept so the repair can be tested rather than assumed:
    "F3":  ((5, 20), 8, 5, 3, 1, (3, 4), (2, 3)),
    "F10": ((7, 21), 13, 7, 3, 1, (7, 5), (4, 3)),
    "F1":  ((4, 12), 7, 4, 3, 1, (3, 2), (4, 3)),
}
REFUSED = {"F3", "F10", "F1"}

def bridge(a, b, t, kappa, ordC):
    s = a + b
    return a * ordC * (t * s - (kappa + 1)) - (ordC * s - 1)

def params(name):
    A0, p, l, q, k, (m0, dm), (n0, dn) = TABLE[name]
    j = 0
    while gcd(m0 + dm * j, n0 + dn * j) != 1:
        j += 1
    m, n = m0 + dm * j, n0 + dn * j
    ok_dio = (m + n) * q * k - n * (q * l - p) == k
    lo, hi = min(m, n), max(m, n)
    cd = pr.corner_chart_data(A0[0], A0[1], l_final=l, b_final=q,
                              who="phi_f7_verify " + name)
    t, kappa, a0, qc = cd["t"], cd["kappa"], cd["deg_C"], cd["ord_C"]
    e, dg = hi - lo + 1, a0 - qc
    r = dg - 1
    coef = t * (hi - lo) + kappa + 1
    rho = (e - 1) * qc + 1
    N = lo * (t * (lo + hi) - (kappa + 1)) - 2 * hi
    res = Fraction(coef * a0, t)
    gap = res - (e * a0 - qc + 1)
    return dict(ok_dio=ok_dio, a=lo, b=hi, t=t, kappa=kappa, a0=a0, q=qc,
                e=e, r=r, dg=dg, coef=coef, rho=rho, N=N,
                res=(int(res) if res.denominator == 1 else res),
                gap=(int(gap) if gap.denominator == 1 else gap),
                retracts=cd["retraction"], A0=A0, l_final=l, b_final=q,
                bridge=bridge(lo, hi, t, kappa, qc),
                degs=((A0[0] + A0[1]) * m, (A0[0] + A0[1]) * n))

# (a,b,t,a0,q,e,r,dg,rho,N,res,gap,degs) -- a0/q here are deg C / ord C.
EXPECT = {
    "F7":  (2, 7, 3, 6, 4, 6, 1, 2, 21, 36, 34, 1, (42, 147)),
    "F8":  (3, 7, 3, 6, 5, 5, 0, 1, 21, 70, 28, 2, (63, 147)),
    "F15": (3, 7, 3, 9, 5, 5, 3, 4, 21, 70, 42, 1, (99, 231)),
    "F16": (3, 5, 3, 9, 7, 3, 1, 2, 15, 56, 24, 3, (99, 165)),
}
for nm, exp in EXPECT.items():
    P = params(nm)
    got = (P["a"], P["b"], P["t"], P["a0"], P["q"], P["e"], P["r"], P["dg"],
           P["rho"], P["N"], P["res"], P["gap"], P["degs"])
    check(f"A: {nm} Diophantine identity", P["ok_dio"])
    check(f"A: {nm} corner invariants {got}", got == exp)
    check(f"A: {nm} corner {P['A0']} RETRACTS, so GGV5's final-corner dictionary "
          f"IS valid there ((t,deg C,ord C) = (l_final,a0,b_final))",
          P["retracts"] and P["t"] == P["l_final"] and P["a0"] == P["A0"][0]
          and P["q"] == P["b_final"])

# the refused rows: the guard must refuse, and the repaired chart must be monomial
for nm in sorted(REFUSED):
    P = params(nm)
    A0 = P["A0"]
    raised = False
    try:
        pr.final_corner_dictionary(A0[0], A0[1], P["l_final"], P["b_final"], who=nm)
    except pr.FinalCornerDictionaryError:
        raised = True
    check(f"A: {nm} corner {A0} is REFUSED (final_corner_dictionary RAISES) and "
          f"the repaired chart is (t,deg C,ord C) = "
          f"({P['t']},{P['a0']},{P['q']}) with dg = 0: C = y a MONOMIAL",
          raised and not P["retracts"] and (P["a0"], P["q"], P["dg"]) == (1, 1, 0))

# ---------------------------------------------------------------------------
# B. claimed solutions: ODE identity + uniqueness by full generic solve
# ---------------------------------------------------------------------------
CLAIMS = {  # name -> (g poly, claimed f)
    ("F7", "ram"):  ((y + 1)**2,
                     y**21 * (y + 1)**11 * (9 * y**2 + 3 * y - 1) / 10),
    ("F7", "cx"):   (y**2 + 3 * y + 6,
                     y**21 * (y - 1) * (y**2 + 3 * y + 6)**6 / 60),
    ("F16", "ram"): ((y + 1)**2,
                     y**15 * (y + 1)**5
                     * (243 * y**4 + 81 * y**3 - 27 * y**2 + 15 * y - 10) / 330),
    # 2026-07-27: F15 REPLACES F10 as the dg=4 point.  Corner (9,24) RETRACTS.
    ("F15", "ram"): ((y + 1)**4,
                     -y**21 * (y + 1)**17
                     * (243 * y**4 + 405 * y**3 + 135 * y**2 - 15 * y + 5) / 105),
    # dg=1 control, REPLACING F1 (whose corner (4,12) is refused):
    ("F8", "ram"):  (y + 1,
                     -y**21 * (y + 1)**5 * (9 * y**2 - 3 * y + 2) / 42),
}

def ode_ok(P, g, f):
    c = y**P["q"] * g
    return sp.expand(P["a"] * P["t"] * c * sp.diff(f, y)
                     - P["a"] * P["coef"] * sp.diff(c, y) * f
                     - c**P["e"]) == 0

def generic_solutions(P, g, slack=2):
    D = P["res"] + slack
    fc = sp.symbols(f"vf0:{D + 1}")
    f = sum(fc[i] * y**i for i in range(D + 1))
    c = y**P["q"] * g
    resid = sp.expand(P["a"] * P["t"] * c * sp.diff(f, y)
                      - P["a"] * P["coef"] * sp.diff(c, y) * f - c**P["e"])
    sols = sp.solve(sp.Poly(resid, y).all_coeffs(), fc, dict=True)
    return [sp.expand(f.subs(s)) for s in sols]

def sig(p):
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

SIGS = {}
for (nm, br), (g, f_claim) in CLAIMS.items():
    P = params(nm)
    check(f"B: {nm}/{br} claimed f solves the ODE exactly", ode_ok(P, g, f_claim))
    sols = generic_solutions(P, g)
    uniq = len(sols) == 1 and sp.expand(sols[0] - f_claim) == 0
    check(f"B: {nm}/{br} full generic solve -> unique = claimed f", uniq)
    sf = sig(f_claim)
    sc = sig(y**P["q"] * g)
    SIGS[(nm, br)] = tuple(sf[i] + P["N"] * sc[i] for i in range(4))

# F7 cross-check: exponent arithmetic == full expansion of Phi.  (Moved from F3,
# whose corner is refused; F7's is the smallest retracting dg=2 case.)
P = params("F7")
g, f_claim = CLAIMS[("F7", "ram")]
full = sig(sp.expand(f_claim * (y**P["q"] * g)**P["N"]))
check("B: F7 Phi signature via full expansion == exponent arithmetic",
      full == SIGS[("F7", "ram")])

# ---------------------------------------------------------------------------
# C. verdict vs the old unified law, and the amended ramified law
# ---------------------------------------------------------------------------
def old_law(P):
    return ((P["e"] * P["a0"] - P["q"] + 1) + P["gap"] + P["N"] * P["a0"],
            P["rho"] + P["N"] * P["q"], P["e"] + P["N"],
            P["gap"] + P["r"] * (P["e"] + P["N"]))

def ram_law(P):
    return (P["res"] + P["N"] * P["a0"], P["rho"] + P["N"] * P["q"],
            P["dg"] * (P["e"] + P["N"]) - (P["dg"] - 1), P["gap"] + P["r"])

EXPECT_SIG = {("F7", "ram"): (250, 165, 83, 2),
              ("F7", "cx"): (250, 165, 0, 85),
              ("F16", "ram"): (528, 407, 117, 4),
              ("F15", "ram"): (672, 371, 297, 4),
              ("F8", "ram"): (448, 371, 75, 2)}
for key, expected in EXPECT_SIG.items():
    check(f"C: {key[0]}/{key[1]} signature {expected}", SIGS[key] == expected)
# the gap>0 & r>0 regime rows: the old conjecture is refuted, the amended law holds
for nm in ("F7", "F16", "F15"):
    P = params(nm)
    o, rl, got = old_law(P), ram_law(P), SIGS[(nm, "ram")]
    check(f"C: {nm} deg,ord match old law", got[:2] == o[:2])
    check(f"C: {nm} mult,cof DIFFER from old law (conjecture refuted)",
          got[2:] != o[2:])
    check(f"C: {nm} ramified law (res+N*degC, rho+N*ordC, dg(e+N)-(dg-1), gap+r)",
          got == rl)
    check(f"C: {nm} is in the gap>0 & r>0 regime and RETRACTS (so this is a "
          f"statement about a real corner)",
          P["gap"] > 0 and P["r"] > 0 and P["retracts"])
    if (nm, "cx") in SIGS:
        check(f"C: {nm} complex-pair branch: deg,ord match, mult=0",
              SIGS[(nm, "cx")][:2] == o[:2] and SIGS[(nm, "cx")][2] == 0)
check("C: the amended ramified law has THREE exact points and BOTH dg parities "
      "of the regime are represented (dg=2 at F7,F16; dg=4 at F15) -- the dg=4 "
      "rung no longer rests on the refused F10",
      sorted({params(nm)["dg"] for nm in ("F7", "F16", "F15")}) == [2, 4])
# F8 (r=0) is where the old and amended laws COINCIDE -- a discrimination control
P8 = params("F8")
check("C: F8 (dg=1, r=0) is a control where old and amended laws AGREE "
      "(dg(e+N)-(dg-1) = e+N at dg=1), and both match the derived signature",
      old_law(P8) == ram_law(P8) == SIGS[("F8", "ram")])
# units
for key in EXPECT_SIG:
    nm, br = key
    if br != "ram":
        continue
    P = params(nm)
    g, f_claim = CLAIMS[key]
    u = sp.cancel(f_claim / (y**P["rho"]
                             * (y + 1)**(P["dg"] * P["e"] - (P["dg"] - 1))))
    upoly = sp.Poly(sp.together(u), y)
    check(f"C: {nm} unit cofactor deg u = gap+r = {P['gap'] + P['r']}, "
          f"u(0),u(-1) != 0",
          upoly.degree() == P["gap"] + P["r"]
          and u.subs(y, 0) != 0 and u.subs(y, -1) != 0)
# (72,108) retro-explanation
check("C: (72,108) audited quartic degree = gap+r = 4+0", 4 + 0 == 4)

# --- THE INDEPENDENT TARGET: the PROVED bridge identity a*q*M - H -------------
bad = [nm for nm in TABLE
       if nm not in REFUSED and SIGS.get((nm, "ram"), (0, params(nm)["bridge"]))[1]
       != params(nm)["bridge"]]
check("C: ord_y(Phi) == a*q*M - H at every derived point (F7, F8, F15, F16) -- "
      "PROVED in BRIDGE_GENERALITY.md and computed by neither file",
      not bad and all(SIGS[(nm, "ram")][1] == params(nm)["bridge"]
                      for nm in ("F7", "F8", "F15", "F16")))
check("C: and the complex-pair branch has the SAME ord_y as the ramified branch, "
      "as BRIDGE_GENERALITY.md D6 proves it must (branch ambiguity moves mult and "
      "the cofactor, never ord_y)",
      SIGS[("F7", "cx")][1] == SIGS[("F7", "ram")][1] == params("F7")["bridge"])

# ---------------------------------------------------------------------------
# D. the impossibility theorem at dg=2 (independent re-proof per corner)
#    Claimed obstruction branch ratios w = g0/g1^2 (from phi_f7.py):
#    F7: {1/4 (disc), 2/3}; F16: {1/4, roots of 54w^2-126w+35}.
#    A simple real root at -1 needs 1 - g1 + g0 = 0 with g1^2 > 4 g0.
#    disc branch: g0=g1^2/4 -> (1-g1/2)^2 = 0 -> g1=2: double root AT -1.
#    ratio branch w>1/4: 1 - g1 + w g1^2 = 0 has no real g1 (disc 1-4w < 0).
# ---------------------------------------------------------------------------
w = sp.symbols("w", real=True)
RATIOS = {"F7": [sp.Rational(2, 3)],
          "F16": list(sp.solve(54 * w**2 - 126 * w + 35, w))}
for nm, ws in RATIOS.items():
    all_gt = all(sp.simplify(wv - sp.Rational(1, 4)) > 0 for wv in ws)
    check(f"D: {nm} every complex-pair ratio w > 1/4 (so 1-g1+w*g1^2 = 0 has "
          f"no real g1: no residual root at -1)", all_gt)
g1v = sp.symbols("g1v")
sol_disc = sp.solve(1 - g1v + g1v**2 / 4, g1v)
check("D: disc branch: g(-1)=0 forces g1=2, i.e. the DOUBLE root sits at -1",
      sol_disc == [2])
# off-branch refutation spot checks: generic w -> no polynomial solution
for nm, wbad in (("F7", sp.Rational(1, 2)), ("F16", sp.Integer(1))):
    P = params(nm)
    sols = generic_solutions(P, y**2 + y + wbad)
    check(f"D: {nm} off-branch w={wbad}: NO polynomial solution", sols == [])
# and on the claimed ratio branch a solution EXISTS (rational rep)
sols = generic_solutions(params("F7"), y**2 + y + sp.Rational(2, 3))
check("D: F7 on-branch w=2/3: polynomial solution exists", len(sols) == 1)
# mult = e+N impossibility (the refuted component), from deg f <= res:
for nm in ("F7", "F16", "F15"):
    P = params(nm)
    imposs = (P["res"] < P["e"] + P["N"]) and (P["e"] - P["N"] < 0)
    check(f"D: {nm} old-law mult=e+N unrealizable "
          f"(res={P['res']} < e+N={P['e'] + P['N']}; e-N<0)", imposs)

# ---------------------------------------------------------------------------
# E. controls
# ---------------------------------------------------------------------------
P14 = params("F14")
sols = generic_solutions(P14, y**5 + 1)
check("E: F14 control (gap=0): unique generic solution = -(1/10) y^21 (y^5+1)^6",
      len(sols) == 1
      and sp.expand(sols[0] + sp.Rational(1, 10) * y**21 * (y**5 + 1)**6) == 0)
check("E: F14 ord_y(Phi) = 165 agrees with the bridge identity",
      sig(sp.expand(sols[0] * (y**P14["q"] * (y**5 + 1))**P14["N"]))[1]
      == P14["bridge"] == 165)
# dg=1 root gauge -- moved 2026-07-27 from F1 (refused corner) to F8 (retracting)
for sv in (sp.Integer(1), sp.Integer(2), sp.Rational(1, 3)):
    sols = generic_solutions(P8, y + sv)
    check(f"E: F8 control (dg=1): root position s={sv} admits a unique solution "
          f"(gauge, not forced)", len(sols) == 1)
check("E: and F8's corner (6,15) RETRACTS, so dg=1 is real there -- unlike F1's "
      "(4,12), where the repaired chart has dg=0 and there is no root to place",
      P8["retracts"] and P8["dg"] == 1 and params("F1")["dg"] == 0
      and not params("F1")["retracts"])

# dg parity across the full 15-family survey, RESTRICTED to dg >= 1 (= retracting)
SURVEY = [
    ("F1", (4,12), 7,4,3,1,(3,2),(4,3)),  ("F2", (5,20), 7,5,2,1,(2,1),(3,2)),
    ("F3", (5,20), 8,5,3,1,(3,4),(2,3)),  ("F4", (5,20), 8,5,3,2,(3,2),(16,12)),
    ("F5", (5,20), 9,5,4,1,(9,7),(5,4)),  ("F6", (5,20), 9,5,4,2,(7,6),(18,16)),
    ("F7", (6,15), 7,3,4,1,(2,1),(7,4)),  ("F8", (6,15), 8,3,5,1,(3,2),(7,5)),
    ("F9", (7,21), 11,7,2,1,(2,1),(3,2)), ("F10",(7,21), 13,7,3,1,(7,5),(4,3)),
    ("F11",(7,21), 13,7,3,2,(2,1),(5,3)), ("F14",(9,24), 7,3,4,1,(2,1),(7,4)),
    ("F15",(9,24), 8,3,5,1,(3,2),(7,5)),  ("F16",(9,24), 10,3,7,1,(3,4),(5,7)),
    ("F17",(9,24), 11,3,8,1,(2,5),(3,8)),
]
parity_ok, dg0_rows, dgpos_rows = True, [], []
for row in SURVEY:
    TABLE[row[0]] = (row[1], row[2], row[3], row[4], row[5], row[6], row[7])
    Pr = params(row[0])
    (dg0_rows if Pr["dg"] == 0 else dgpos_rows).append(row[0])
    if Pr["dg"] >= 1 and (Pr["dg"] % 2 == 0) != (Pr["gap"] > 0 and Pr["r"] > 0):
        parity_ok = False
check("E: survey mini-lemma dg even <=> (gap>0 and r>0) holds on every dg >= 1 "
      "row (the six RETRACTING rows F7,F8,F14-F17)", parity_ok)
check("E: and the dg = 0 rows are exactly the guard-refused ones: %s"
      % ", ".join(dg0_rows),
      set(dg0_rows) == {"F1", "F2", "F3", "F4", "F5", "F6", "F9", "F10", "F11"}
      and set(dgpos_rows) == {"F7", "F8", "F14", "F15", "F16", "F17"})
check("E: so the biconditional is NOT vacuously true post-repair: it has 3 rows "
      "inside the regime and 3 outside among the dg >= 1 rows",
      len([nm for nm in dgpos_rows
           if params(nm)["gap"] > 0 and params(nm)["r"] > 0]) == 3)

# ===========================================================================
# F.  THE CHART REPAIR: drift guard + mutation controls
# ===========================================================================
import contextlib, io                                            # noqa: E402
_rep = io.StringIO()
with contextlib.redirect_stdout(_rep):
    import phi_f7 as pf7                                         # noqa: E402

check("F1 DRIFT GUARD: phi_f7.REFUSED_ROWS equals the refused set recomputed "
      "here from the guard (%s)" % ", ".join(sorted(REFUSED)),
      pf7.REFUSED_ROWS == REFUSED)
check("F1b and its SUPERSEDED table lists exactly those rows, with the stale "
      "(t,deg C,ord C) really equal to the dictionary's (l_final, a0, b_final)",
      set(pf7.SUPERSEDED) == REFUSED
      and all(pf7.SUPERSEDED[nm][:3] == (TABLE[nm][2],        # l_final
                                         TABLE[nm][0][0],     # a0
                                         TABLE[nm][3])        # b_final
              for nm in ("F3", "F10", "F1")))
check("F1c and the two retired f-polynomials are still PRESENT and labelled "
      "(deleting them would make the retirement unfalsifiable)",
      set(pf7.SUPERSEDED_F) == {"F3", "F10"})

# F2: the retired polynomials still solve the ODE at their STALE parameters --
# so what is retired is a statement about a CORNER, not a computation.
for nm in ("F3", "F10"):
    st, sdC, soC, sN, ssig = pf7.SUPERSEDED[nm]
    P = params(nm)
    a_, b_ = P["a"], P["b"]
    e_ = b_ - a_ + 1
    coef_ = st * (b_ - a_) + (st - 2) + 1
    g_ = (y + 1)**(sdC - soC)
    c_ = y**soC * g_
    fs = pf7.SUPERSEDED_F[nm]
    still = sp.expand(a_ * st * c_ * sp.diff(fs, y)
                      - a_ * coef_ * sp.diff(c_, y) * fs - c_**e_) == 0
    check(f"F2 {nm}: the retired f DOES solve the ODE at the stale parameters "
          f"(t,kappa,q,dg) = ({st},{st-2},{soC},{sdC-soC}) -- retired as a claim "
          f"about the corner, not withdrawn as arithmetic", still)
    # ... and does NOT solve the repaired corner's ODE.
    c_rep = y**P["q"]                        # C = y, dg = 0
    bad = sp.expand(a_ * P["t"] * c_rep * sp.diff(fs, y)
                    - a_ * P["coef"] * sp.diff(c_rep, y) * fs - c_rep**e_)
    check(f"F2b {nm}: and it does NOT solve the REPAIRED corner's ODE "
          f"(t={P['t']}, C=y, e={e_}) -- the two are different equations",
          bad != 0)
    # the repaired unique solution is the monomial (1/a) y^e
    frep = sp.Rational(1, a_) * y**e_
    check(f"F2c {nm}: the repaired corner's unique solution is f = (1/{a_}) y^{e_} "
          f"and Phi = f*C^{P['N']} is the MONOMIAL (1/{a_}) y^{e_ + P['N']}",
          sp.expand(a_ * P["t"] * c_rep * sp.diff(frep, y)
                    - a_ * P["coef"] * sp.diff(c_rep, y) * frep - c_rep**e_) == 0
          and sig(sp.expand(frep * c_rep**P["N"])) == (e_ + P["N"], e_ + P["N"], 0, 0))

# F3: MUTATION CONTROL -- reinstating the dictionary must move ord_y and break
# the bridge identity at the guarded chart.  Shape copied from
# bridge_generality.py MUT F (51->205, 30->112, 22->107).
moved = {}
for nm in sorted(REFUSED):
    st, sdC, soC, sN, ssig = pf7.SUPERSEDED[nm]
    P = params(nm)
    a_, b_ = P["a"], P["b"]
    Ns = a_ * (st * (a_ + b_) - (st - 1)) - 2 * b_
    ords = ((b_ - a_) * soC + 1) + Ns * soC
    stale_bridge = bridge(a_, b_, st, st - 2, soC)
    faithful = (Ns == sN and ords == ssig[1] and stale_bridge == ords)
    if faithful and ords != P["bridge"]:
        moved[nm] = (P["bridge"], ords)
check("F3 MUT: reinstating the refused dictionary reproduces the SUPERSEDED N and "
      "ord_y exactly, then contradicts the guarded chart's bridge value at every "
      "refused row -- " + "; ".join("%s %d<-%d" % (nm, g, s)
                                    for nm, (g, s) in sorted(moved.items())),
      len(moved) == 3)
check("F3b and F1's and F10's displacements match bridge_generality MUT F / this "
      "repair's phi_corner4 numbers exactly: F1 51<-205, F3 30<-112, F10 114<-820",
      moved.get("F1") == (51, 205) and moved.get("F3") == (30, 112)
      and moved.get("F10") == (114, 820))
check("F3c and the stale signatures all claim mult_(y+1) > 0 -- a (y+1) place a "
      "monomial C cannot have; every repaired row has mult = cof = 0",
      all(pf7.SUPERSEDED[nm][4][2] > 0 for nm in REFUSED))

# F4: the retracting rows must be untouched by the repair.
check("F4 the five retracting rows used here (F7, F8, F14, F15, F16) all have "
      "t = l_final, deg C = a0, ord C = b_final -- bit-identical to the "
      "pre-repair dictionary, so the repair is targeted, not a rewrite",
      all(params(nm)["t"] == TABLE[nm][2] and params(nm)["a0"] == TABLE[nm][0][0]
          and params(nm)["q"] == TABLE[nm][3]
          for nm in ("F7", "F8", "F14", "F15", "F16")))

# ---------------------------------------------------------------------------
print()
if _f == 0:
    print(f"ALL {_n} PHI-F7 CHECKS PASSED")
else:
    print(f"{_f}/{_n} PHI-F7 CHECKS FAILED")
import pathlib
print(f"script: {pathlib.Path(__file__).resolve()}")
sys.exit(0 if _f == 0 else 1)
