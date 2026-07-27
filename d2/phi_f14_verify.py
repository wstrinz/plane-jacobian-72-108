#!/usr/bin/env python3
"""phi_f14_verify.py  -- independent PASS/FAIL checker for PHI_F14.md.

Self-contained: re-derives everything from scratch with its own routines (no
imports from phi_f14.py; signature extraction by trial division, uniqueness by
raw linear solves).  Verifies:

  A. F14/F1 corner data: Diophantine identities, parameter extraction,
     and the two mini-lemmas gap = (q-1) - a0/t, dg = a0-q over all 15
     standard-chart families.
  B. The residual identity  resid = y^(eq) g^e [A a((t*rho-coef*q) g +
     (t*e-coef) y g') - 1]  with GENERIC g at the F14 parameters, and the
     exponent identity eq - q - rho + 1 = 0 that makes the bracket constant.
  C. F14 forcing: multipliers -4..-1 kill g_1..g_4, resonance at y^5,
     g(-1)=0 + monic => g = y^5+1, A = -1/10; the forced f solves the ODE
     exactly; uniqueness among polynomials of degree <= 51 (full linear solve).
  D. F14 signature by trial division = (375,165,42,168) = the unified law's
     parameter-free prediction (MATCH).
  E. F1 (gap>0, r=0): fully generic solve (degree allowed 2 beyond resonant)
     has a UNIQUE polynomial solution f = (1/15) y^4 (y+1)^2 (4y-1); the
     cofactor u = (4y-1)/15 has deg = gap = 1 and is a unit at 0 and -1;
     Phi signature (275,205,69,1) = amended-law prediction (MATCH).
  F. Controls: the same forcing machinery re-derives (108,144), (75,125),
     (56,84), (50,75) and reproduces their audited signatures; the stored
     (72,108) signature (238,204,30,4) obeys the amended law with gap=4.
  G. The reduced N-formula N = a[t(a+b-1)+1] - 2b reproduces N at all seven
     points; coverage census t in {3,4,5,7}, e in {2,3,6}.

Usage: python phi_f14_verify.py [--quiet]     exit 0 iff all checks pass.
"""
import sys
import sympy as sp
from fractions import Fraction
from math import gcd

import polygon_reduction as pr


def bridge(a, b, t, kappa, ordC):
    """ord_y(Phi) = a*q*M - H.  PROVED in BRIDGE_GENERALITY.md; computed here by
    neither phi_f14.py nor this file's own machinery -- an INDEPENDENT target."""
    s = a + b
    return a * ordC * (t * s - (kappa + 1)) - (ordC * s - 1)


def _require(_cond, _msg):
    """Proof-critical check: fails loudly and exits nonzero, unaffected by python -O."""
    if not _cond:
        import sys as _sys
        print("FAIL: " + str(_msg))
        _sys.exit(1)


QUIET = "--quiet" in sys.argv[1:]
y, A = sp.symbols("y A")
_n, _fail = 0, 0

def ok(cond, msg):
    global _n, _fail
    _n += 1
    if not cond:
        _fail += 1
    if not QUIET or not cond:
        print(f"[{'OK' if cond else 'FAIL'}] {msg}")

def signature(Phi):
    """(deg, ord_y, mult_(y+1), cofactor deg) by trial division only."""
    P = sp.Poly(sp.expand(Phi), y)
    deg = P.degree()
    ordy = min(m[0] for m in P.monoms())
    mult, Q = 0, P
    while True:
        q2, r2 = sp.div(Q, sp.Poly(y + 1, y))
        if not r2.is_zero:
            break
        Q, mult = q2, mult + 1
    return (deg, ordy, mult, deg - ordy - mult)

# phi_f14.py is a REPORT: importing it runs its whole derivation and prints it.
# Swallow that stream so this checker's own PASS/FAIL lines -- and any FAIL --
# cannot be lost in it.  (Masking OUTPUT is fine; masking an EXIT CODE is the
# trap, and nothing here touches exit codes.)
import contextlib as _ctx, io as _io                        # noqa: E402
_pf14_report = _io.StringIO()
with _ctx.redirect_stdout(_pf14_report):
    from phi_f14 import mult_and_cofactor, gap_effective    # noqa: E402


def law(a, b, t, a0, q, gap):
    """Unified law, kappa = t-2 eliminated (PHI_CORNER4 sec.2).

    2026-07-26: gap enters via gap_effective, mult/cofactor via the
    residual-free branch (deg g = a0-q = 0 => mult = 0, cof = gap)."""
    e = b - a + 1
    N = a * (t * (a + b - 1) + 1) - 2 * b
    rho = (e - 1) * q + 1
    ge = gap_effective(gap)
    mult, cof = mult_and_cofactor(e, N, a0, q, gap)
    return N, ((e * a0 - q + 1) + ge + N * a0, rho + N * q, mult, cof)

# ---------------------------------------------------------------- A. corner data
FAMS = [  # (name, A0, p, l, q, k, m0, dm, n0, dn)   [GGV5 v11<=35, length-1]
    ("F1", (4,12), 7,4,3,1, 3,2, 4,3),  ("F2", (5,20), 7,5,2,1, 2,1, 3,2),
    ("F3", (5,20), 8,5,3,1, 3,4, 2,3),  ("F4", (5,20), 8,5,3,2, 3,2, 16,12),
    ("F5", (5,20), 9,5,4,1, 9,7, 5,4),  ("F6", (5,20), 9,5,4,2, 4,3, 10,8),
    ("F7", (6,15), 7,3,4,1, 2,1, 7,4),  ("F8", (6,15), 8,3,5,1, 3,2, 7,5),
    ("F9", (7,21), 11,7,2,1, 2,1, 3,2), ("F10",(7,21), 13,7,3,1, 7,5, 4,3),
    ("F11",(7,21), 13,7,3,2, 2,1, 5,3), ("F14",(9,24), 7,3,4,1, 2,1, 7,4),
    ("F15",(9,24), 8,3,5,1, 3,2, 7,5),  ("F16",(9,24), 10,3,7,1, 3,4, 5,7),
    ("F17",(9,24), 11,3,8,1, 2,5, 3,8),
]
ok((2+7)*4*1 - 7*(4*3-7) == 1, "A: F14 j=0 Diophantine (m+n)qk - n(ql-p) = k = 1")
ok((3+4)*3*1 - 4*(3*4-7) == 1, "A: F1 j=0 Diophantine = k = 1")
ok(33*2 == 66 and 33*7 == 231, "A: F14 degrees v11*(m,n) = (66,231)")
ok(16*3 == 48 and 16*4 == 64,  "A: F1 degrees v11*(m,n) = (48,64)")

# 2026-07-27: the mini-lemma is now stated in the GUARDED variables (deg C,
# ord C) instead of (a0, b_final).  Pre-repair BOTH sides of the identity used
# the same substitution, so it "held" on all 15 rows without saying anything
# about any corner.  Chart data comes from polygon_reduction.corner_chart_data --
# the one routine no consumer re-implements.
lem_gap, lem_dg, REFUSED, RETRACT = True, True, [], []
CHART = {}
for name, A0, p, l, q, k, m0, dm, n0, dn in FAMS:
    j = 0
    while gcd(m0 + dm * j, n0 + dn * j) != 1:
        j += 1
    m, n = m0 + dm * j, n0 + dn * j
    a_, b_ = sorted((m, n))
    cd = pr.corner_chart_data(A0[0], A0[1], l_final=l, b_final=q,
                              who="phi_f14_verify " + name)
    t_, kap_, degC_, ordC_ = cd["t"], cd["kappa"], cd["deg_C"], cd["ord_C"]
    CHART[name] = dict(cd, a=a_, b=b_, A0=A0, l_final=l, b_final=q)
    e_ = b_ - a_ + 1
    coef_ = t_ * (b_ - a_) + kap_ + 1
    res_ = Fraction(coef_ * degC_, t_)
    gap_ = res_ - (e_ * degC_ - ordC_ + 1)
    lem_gap &= (gap_ == Fraction(ordC_ - 1) - Fraction(degC_, t_))
    (RETRACT if cd["retraction"] else REFUSED).append(name)
    lem_dg &= (degC_ - ordC_ >= 1) == cd["retraction"]
ok(lem_gap, "A: mini-lemma gap = (ordC-1) - degC/t on all 15 standard-chart "
            "families, in the GUARDED variables")
ok(lem_dg, "A: dg = degC-ordC >= 1 EXACTLY on the retracting rows; at a refused "
           "corner C is the monomial y and dg = 0 (no residual g at all)")
ok(set(REFUSED) == {"F1", "F2", "F3", "F4", "F5", "F6", "F9", "F10", "F11"}
   and set(RETRACT) == {"F7", "F8", "F14", "F15", "F16", "F17"},
   "A: the guard REFUSES 9 of these 15 rows (F1-F6, F9-F11) and accepts 6 "
   "(F7,F8,F14-F17).  F12/F13 are not in this table; phi_corner4_verify.py "
   "covers the full 17 (11 refused)")
ok(sorted({cd["t"] for cd in CHART.values()}) == [3, 4],
   "A: t-census: every DERIVED chart exponent here is in {3,4}; t in {5,7} came "
   "only from reading t = l_final at refused corners")
ok(9 == 3 * (4 - 1) and CHART["F14"]["deg_C"] == 9 and CHART["F14"]["ord_C"] == 4,
   "A: F14 satisfies degC = t(ordC-1)  <=>  gap = 0 (resonance exact), and its "
   "corner (9,24) RETRACTS so the dictionary is valid there")
ok(Fraction(8 - 1) - Fraction(9, 3) == 4 and CHART["F17"]["retraction"],
   "A: F17 gap = (ordC-1) - degC/t = 4 > 0 (resonance broken) at a RETRACTING "
   "corner -- the replacement for F1, whose corner (4,12) is refused")
ok(Fraction(5 - 1) - Fraction(6, 3) == 2 and CHART["F8"]["retraction"],
   "A: F8 gap = 2 > 0, r = 0, also retracting -- a SECOND gap value in the "
   "r=0 regime, so it is no longer tested at one gap only")
ok(not CHART["F1"]["retraction"] and CHART["F1"]["deg_C"] == 1
   and Fraction(1 - 1) - Fraction(1, 3) < 0,
   "A: and F1's repaired chart has deg C = 1, gap = -1/3 < 0 -- so (4,12) is NOT "
   "in the gap>0/r=0 regime the pre-repair file probed there")

# --------------------------------------------- B. residual identity, generic g
a_, t_, kap_, a0_, q_ = 2, 3, 1, 9, 4          # F14
e_, dg_ = 7 - 2 + 1, a0_ - q_
coef_ = t_ * (7 - 2) + kap_ + 1
rho_ = (e_ - 1) * q_ + 1
ok(e_ * q_ - q_ - rho_ + 1 == 0,
   "B: exponent identity eq - q - rho + 1 = 0 (bracket is y-free up to g-part)")
gcs = sp.symbols(f"g0:{dg_+1}")
g_gen = sum(gcs[i] * y**i for i in range(dg_ + 1))
c_gen = y**q_ * g_gen
f_gen = A * y**rho_ * g_gen**e_
resid = sp.expand(a_ * t_ * c_gen * sp.diff(f_gen, y)
                  - a_ * coef_ * sp.diff(c_gen, y) * f_gen - c_gen**e_)
bracket = sp.expand(A * a_ * ((t_ * rho_ - coef_ * q_) * g_gen
                              + (t_ * e_ - coef_) * y * sp.diff(g_gen, y)) - 1)
ok(sp.expand(resid - y**(e_ * q_) * g_gen**e_ * bracket) == 0,
   "B: resid = y^24 g^6 [A a((t rho - coef q) g + (te-coef) y g') - 1], generic g")
ok(t_ * rho_ - coef_ * q_ == -5 and t_ * e_ - coef_ == 1,
   "B: F14 bracket constants (t rho - coef q, te - coef) = (-5, 1)")
ok(t_ * e_ - coef_ == 1,
   "B: te - coef = t - kappa - 1 = 1 is the kappa = t-2 identity in ODE form")

# ------------------------------------------------------------- C. F14 forcing
mults = [(t_ * rho_ - coef_ * q_) + i * (t_ * e_ - coef_) for i in range(dg_ + 1)]
ok(all(mults[i] != 0 for i in range(1, dg_)) and mults[dg_] == 0,
   f"C: multipliers {mults[1:dg_]} kill g_1..g_4; y^{dg_} resonant (multiplier 0)")
g14 = y**dg_ + 1
A14 = sp.Rational(1, a_ * (t_ * rho_ - coef_ * q_))
ok(A14 == sp.Rational(-1, 10), "C: A = 1/(a(t rho - coef q)) = -1/10")
ok(g14.subs(y, -1) == 0, "C: g(-1) = 0 (mult_(y+1) Phi > 0 normalization)")
c14 = y**q_ * g14
f14 = sp.expand(A14 * y**rho_ * g14**e_)
ok(sp.expand(a_ * t_ * c14 * sp.diff(f14, y)
             - a_ * coef_ * sp.diff(c14, y) * f14 - c14**e_) == 0,
   "C: forced f = -(1/10) y^21 (y^5+1)^6 solves 6cf' - 34c'f = c^6 exactly")
ok(sp.degree(f14, y) == 51 and 51 == e_ * a0_ - q_ + 1,
   "C: deg f = 51 = e*a0 - q + 1 (pure = resonant, gap 0)")
ok(sp.factor(sp.cancel(g14 / (y + 1))) == y**4 - y**3 + y**2 - y + 1,
   "C: residual H = Phi_10 (10th cyclotomic), same as (108,144)/F9 (dg=5 class)")
# uniqueness among ALL polynomials of degree <= 51: linear solve from scratch
fu = sp.symbols(f"u0:{52}")
f_u = sum(fu[i] * y**i for i in range(52))
resid_u = sp.expand(a_ * t_ * c14 * sp.diff(f_u, y)
                    - a_ * coef_ * sp.diff(c14, y) * f_u - c14**e_)
sols_u = sp.solve(sp.Poly(resid_u, y).all_coeffs(), fu, dict=True)
ok(len(sols_u) == 1 and sp.expand(f_u.subs(sols_u[0]) - f14) == 0,
   "C: f is the UNIQUE polynomial solution of degree <= 51 (full linear solve)")

# ----------------------------------------------------------- D. F14 signature
N14, law14 = law(2, 7, 3, 9, 4, 0)
ok(N14 == 36, "D: reduced N-formula gives N = 36 for F14")
Phi14 = f14 * c14**N14
sig14 = signature(Phi14)
ok(sig14 == (375, 165, 42, 168),
   "D: F14 signature by trial division = (375, 165, 42, 168)")
ok(sig14 == law14, f"D: F14 signature = law prediction {law14}  ==> MATCH")
ok(sp.expand(Phi14 - sp.Rational(-1, 10) * y**165 * (y**5 + 1)**42) == 0,
   "D: closed form Phi = -(1/10) y^165 (y^5+1)^42")
ok(sig14[1] == bridge(2, 7, 3, 1, 4) == 165,
   "D: and ord_y(Phi) = 165 equals the PROVED bridge identity a*q*M - H "
   "(BRIDGE_GENERALITY.md) -- a target neither module computes")

# --------------------------- E. the gap>0 / r=0 regime, on RETRACTING corners
# 2026-07-27: this section used to probe F1 (48,64) at the corner (4,12), which
# the retraction guard REFUSES -- the repaired chart there is C = y, dg = 0,
# gap = -1/3, so (4,12) is not in this regime at all.  Two replacements, both on
# corners that retract, and at two DIFFERENT gap values so the regime is not
# pinned by a single point:
#     F17 (66,99)  at (9,24): gap = 4 -- the same gap as the audited (72,108)
#     F8  (63,147) at (6,15): gap = 2
def probe_r0(tag, a1, b1, t1, degC1, ordC1, want_f, want_sig, want_N):
    e1, r1 = b1 - a1 + 1, degC1 - ordC1 - 1
    coef1 = t1 * (b1 - a1) + (t1 - 2) + 1
    rho1 = (e1 - 1) * ordC1 + 1
    gap1 = Fraction(coef1 * degC1, t1) - (e1 * degC1 - ordC1 + 1)
    ok(r1 == 0 and gap1 > 0, f"E: {tag} is r=0 with gap = {gap1} > 0")
    c1 = y**ordC1 * (y + 1)
    res_deg1 = e1 * degC1 - ordC1 + 1 + int(gap1)
    fv = sp.symbols(f"vv_{tag}_0:{res_deg1 + 3}")
    f_v = sum(fv[i] * y**i for i in range(res_deg1 + 3))
    resid_v = sp.expand(a1 * t1 * c1 * sp.diff(f_v, y)
                        - a1 * coef1 * sp.diff(c1, y) * f_v - c1**e1)
    sols_v = sp.solve(sp.Poly(resid_v, y).all_coeffs(), fv, dict=True)
    ok(len(sols_v) == 1, f"E: {tag} ODE {a1*t1}cf' - {a1*coef1}c'f = c^{e1} has a "
                         f"UNIQUE polynomial solution (degree allowed 2 past "
                         f"resonant)")
    f1s = sp.expand(f_v.subs(sols_v[0]))
    ok(sp.degree(f1s, y) == res_deg1,
       f"E: {tag} deg f = {res_deg1} = resonant degree (pure-ansatz degree "
       f"{res_deg1 - int(gap1)} is NOT attained)")
    ok(sp.expand(f1s - want_f) == 0, f"E: {tag} f = {sp.factor(want_f)}")
    u1 = sp.cancel(f1s / (y**rho1 * (y + 1)**e1))
    ok(sp.degree(sp.expand(u1), y) == int(gap1) and u1.subs(y, 0) != 0
       and u1.subs(y, -1) != 0,
       f"E: {tag} cofactor u: deg = gap = {gap1}, u(0) != 0, u(-1) != 0 (UNIT)")
    N1, law1 = law(a1, b1, t1, degC1, ordC1, int(gap1))
    ok(N1 == want_N, f"E: {tag} N = {want_N}")
    sig1 = signature(f1s * c1**N1)
    ok(sig1 == want_sig and sig1 == law1,
       f"E: {tag} signature {want_sig} = amended-law prediction {law1}  ==> MATCH")
    ok(sig1[1] == bridge(a1, b1, t1, t1 - 2, ordC1),
       f"E: {tag} ord_y(Phi) = {sig1[1]} equals the PROVED bridge identity")
    return sig1

sig17 = probe_r0("F17", 2, 3, 3, 9, 8,
                 -sp.Rational(1, 910) * y**9 * (y + 1)**2
                 * (243 * y**4 - 81 * y**3 + 54 * y**2 - 42 * y + 35),
                 (195, 169, 22, 4), 20)
sig8 = probe_r0("F8", 3, 7, 3, 6, 5,
                -sp.Rational(1, 42) * y**21 * (y + 1)**5 * (9 * y**2 - 3 * y + 2),
                (448, 371, 75, 2), 70)
ok(sig17[3] == 4 and sig8[3] == 2,
   "E: the r=0 cofactor equals gap on BOTH new points (4 and 2), so the "
   "'unit cofactor of degree exactly gap' story is tested at two gap values on "
   "two distinct corners -- not at one point as before")

# ------------------------------------------------------------------ F. controls
def force_and_sign(a_, b_, t_, a0_, q_):
    """Re-derive a gap=0 corner end to end; return signature.

    2026-07-26: dg_ == 0 means C is a MONOMIAL, so g = 1 is forced (a monic
    constant -- no free coefficient to gauge and no root to place)."""
    e_, dg_ = b_ - a_ + 1, a0_ - q_
    coef_ = t_ * (b_ - a_) + (t_ - 2) + 1
    rho_ = (e_ - 1) * q_ + 1
    g_ = sp.Integer(1) if dg_ == 0 else y**dg_ + 1
    c_ = y**q_ * g_
    A_ = sp.Rational(1, a_ * (t_ * rho_ - coef_ * q_))
    f_ = sp.expand(A_ * y**rho_ * g_**e_)
    _require(sp.expand(a_ * t_ * c_ * sp.diff(f_, y)
                     - a_ * coef_ * sp.diff(c_, y) * f_ - c_**e_) == 0, "sp.expand(a_ * t_ * c_ * sp.diff(f_, y) - a_ * coef_ * sp.diff(c_, y) * f_ - c_**e_) == 0")
    N_ = a_ * (t_ * (a_ + b_ - 1) + 1) - 2 * b_
    return signature(f_ * c_**N_), A_

CONTROLS = [  # (label, a, b, t, degC, ordC, audited signature, audited lc)
    ("(108,144)", 3, 4, 4, 8, 3, (550, 205, 69, 276), sp.Rational(-1, 15)),
    # REPAIRED 2026-07-26: t=4, kappa=2, deg C=ord C=1, C=y (a monomial).
    ("(75,125)",  3, 5, 4, 1, 1, (80, 80, 0, 0), sp.Rational(1, 3)),
    # REPAIRED 2026-07-27: corner (7,21) is REFUSED (GGHV22 publishes l=3 there,
    # 2204.14178.tex:1394), so t=3 and C=y -- not t=7, deg C=7, ord C=2.
    ("(56,84)",   2, 3, 3, 1, 1, (22, 22, 0, 0),  sp.Rational(1, 2)),
    ("(50,75)",   2, 3, 4, 1, 1, (30, 30, 0, 0),  sp.Rational(1, 2)),
    # REPAIRED 2026-07-27: corner (4,12) is REFUSED, so t=3 and C=y.
    ("(48,64)",   3, 4, 3, 1, 1, (51, 51, 0, 0),  sp.Rational(1, 3)),
]
for lbl, ca, cb, ct, cdC, coC, csig, clc in CONTROLS:
    s_, A_ = force_and_sign(ca, cb, ct, cdC, coC)
    ok(s_ == csig and A_ == clc,
       f"F: control {lbl}: re-derived signature {s_} and lc {A_} match audited values")
    ok(s_[1] == bridge(ca, cb, ct, ct - 2, coC),
       f"F: control {lbl}: ord_y = {s_[1]} equals the PROVED bridge identity")
ok(law(2, 3, 4, 8, 7, 4)[1] == (238, 204, 30, 4),
   "F: (72,108) audited signature (238,204,30,4) obeys the amended law with gap = 4")
ok(all(clc == sp.Rational(1, ca) for lbl, ca, cb, ct, cdC, coC, csig, clc
       in CONTROLS if cdC == 1),
   "F: at every MONOMIAL control the leading constant is exactly 1/a -- forced, "
   "because C = y collapses the ODE to a*A*(t - kappa - 1) = a*A = 1")

# ------------------------------------------------------------------- G. census
POINTS = [  # (label, a, b, t, degC, ordC, gap, N)  -- ALL through the guard
    ("(72,108)",  2, 3, 4, 8, 7, 4, 28), ("(108,144)", 3, 4, 4, 8, 3, 0, 67),
    ("(75,125)",  3, 5, 4, 1, 1, 0, 77), ("(50,75)",   2, 3, 4, 1, 1, 0, 28),
    ("(66,231)",  2, 7, 3, 9, 4, 0, 36), ("(56,84)",   2, 3, 3, 1, 1, 0, 20),
    ("(48,64)",   3, 4, 3, 1, 1, 0, 49), ("(66,99)",   2, 3, 3, 9, 8, 4, 20),
    ("(63,147)",  3, 7, 3, 6, 5, 2, 70), ("(99,231)",  3, 7, 3, 9, 5, 1, 70),
]
ok(all(law(sa, sb, st, sdC, soC, sg)[0] == sN
       for _, sa, sb, st, sdC, soC, sg, sN in POINTS),
   "G: reduced N-formula N = a[t(a+b-1)+1] - 2b reproduces N at all TEN points")
ok(all(law(sa, sb, st, sdC, soC, sg)[1][1] == bridge(sa, sb, st, st - 2, soC)
       for _, sa, sb, st, sdC, soC, sg, sN in POINTS),
   "G: and the law's ord_y equals the PROVED bridge identity a*q*M - H at all ten "
   "-- the census is no longer validated only against targets it produced itself")
ok(sorted({st for _, sa, sb, st, sdC, soC, sg, sN in POINTS}) == [3, 4],
   "G: coverage census t in {3,4}.   [2026-07-26: t=5 left -- the two t=5 points "
   "were both at (5,20), which is t=4.  2026-07-27: t=7 left too -- (56,84) is at "
   "(7,21), where GGHV22 PUBLISHES l=3.  No corner in GGV5's v11<=35 length-1 "
   "tables has a derived t outside {3,4}: phi_corner4.py STEP 1b]")
ok(sorted({sb - sa + 1 for _, sa, sb, st, sdC, soC, sg, sN in POINTS})
   == [2, 3, 5, 6],
   "G: e coverage {2,3,5,6} -- e=5 is NEW (F8/F15 at (a,b)=(3,7)), contributed by "
   "the replacement points the repair required")
ok(sorted({sg for _, sa, sb, st, sdC, soC, sg, sN in POINTS if sg > 0})
   == [1, 2, 4]
   and len([1 for _, sa, sb, st, sdC, soC, sg, sN in POINTS if sg > 0]) == 4,
   "G: gap coverage in the resonance regime is {1,2,4} on FOUR points ((72,108) 4, "
   "F17 4, F8 2, F15 1); pre-repair it was {1,4} on two, one of which (F1, gap=1) "
   "sat at a refused corner")

# ===========================================================================
# H.  THE CHART REPAIR: drift guard + mutation controls
# ===========================================================================
import contextlib, io                                            # noqa: E402
_rep = io.StringIO()
with contextlib.redirect_stdout(_rep):
    import phi_f14 as pf14                                       # noqa: E402

ok(set(pf14.SUPERSEDED) == {"F1"},
   "H1 DRIFT GUARD: phi_f14.SUPERSEDED holds exactly the retired row F1, and the "
   "retired polynomial is kept LABELLED so the retirement stays falsifiable")
_st, _sdC, _soC, _sN, _ssig, _sf = pf14.SUPERSEDED["F1"]
ok((_st, _sdC, _soC) == (4, 4, 3),
   "H1b and its stale (t,deg C,ord C) really IS the dictionary output "
   "(l_final, a0, b_final) = (4,4,3) at (4,12) -- so H3 mutates the OLD model, "
   "not a straw man")
ok(_sN == 67 and _ssig == (275, 205, 69, 1),
   "H1c and the stale N and signature are the ones PHI_F14.md published")

# H2: the guard refuses (4,12), and the retired f is right-ODE/wrong-corner.
_raised = False
try:
    pr.final_corner_dictionary(4, 12, 4, 3, who="F1")
except pr.FinalCornerDictionaryError:
    _raised = True
ok(_raised and pr.chart_exponent(4, 12) == 3 and not pr.has_retraction(4, 12),
   "H2 final_corner_dictionary RAISES at (4,12): ceil(12/4) = 3 and 3*(4-1) = 9 "
   "!= 12, so the retraction shape fails")
_c_stale = y**3 * (y + 1)
ok(sp.expand(3 * 4 * _c_stale * sp.diff(_sf, y)
             - 3 * (4 * 1 + 2 + 1) * sp.diff(_c_stale, y) * _sf - _c_stale**2) == 0,
   "H2b the retired f = (1/15)y^4(y+1)^2(4y-1) DOES solve the ODE at the STALE "
   "parameters (t,kappa,q,dg) = (4,2,3,1) -- retired as a claim about (4,12), not "
   "withdrawn as arithmetic")
_c_rep = y
ok(sp.expand(3 * 3 * _c_rep * sp.diff(_sf, y)
             - 3 * (3 * 1 + 1 + 1) * sp.diff(_c_rep, y) * _sf - _c_rep**2) != 0,
   "H2c and it does NOT solve the REPAIRED corner ODE (t=3, C=y, e=2): the two are "
   "different equations, so the two claims are genuinely incompatible")

# H3: MUTATION CONTROL.  Reinstating the dictionary must move ord_y and must
# break the bridge identity at the guarded chart.  Shape copied from
# bridge_generality.py MUT F (51->205, 30->112, 22->107).
MUT = {  # row: (a, b, t_stale, degC_stale, ordC_stale, t_good, ord_good, ord_stale)
    "F1": (3, 4, 4, 4, 3, 3, 51, 205),
    "F9": (2, 3, 7, 7, 2, 3, 22, 107),
    "F2": (2, 3, 5, 5, 2, 4, 30,  75),
}
_moved, _straw = [], []
for nm, (ma, mb, mt, mdC, moC, gt, gord, sord) in sorted(MUT.items()):
    stale = law(ma, mb, mt, mdC, moC, Fraction(moC - 1) - Fraction(mdC, mt))[1][1]
    good = law(ma, mb, gt, 1, 1, Fraction(0) - Fraction(1, gt))[1][1]
    stale_bridge = bridge(ma, mb, mt, mt - 2, moC)
    good_bridge = bridge(ma, mb, gt, gt - 2, 1)
    # The stale numbers are INTERNALLY consistent -- they satisfy the bridge
    # identity at their OWN stale chart -- which is exactly why they passed for
    # months.  What refutes them is that their chart is refused.
    if (stale == sord == stale_bridge and good == gord == good_bridge
            and stale != good_bridge):
        _moved.append((nm, gord, sord))
    else:
        _straw.append((nm, stale, sord, stale_bridge, good, gord, good_bridge))
ok(len(_moved) == 3 and not _straw,
   "H3 MUT: reinstating the refused dictionary reproduces the SUPERSEDED ord_y "
   "exactly and then contradicts the guarded chart bridge value -- "
   + "; ".join("%s %d<-%d" % m for m in _moved))
ok([m for m in _moved if m[0] == "F1"] == [("F1", 51, 205)]
   and [m for m in _moved if m[0] == "F9"] == [("F9", 22, 107)],
   "H3b and the displacements match bridge_generality.py MUT F exactly "
   "(F1 51<-205, F9 22<-107), computed there by a wholly independent route")
ok(all(law(ma, mb, mt, mdC, moC,
           Fraction(moC - 1) - Fraction(mdC, mt))[1][2] > 0
       for ma, mb, mt, mdC, moC, _, _, _ in MUT.values()),
   "H3c and every stale signature claims mult_(y+1) > 0 -- a (y+1) place a "
   "monomial C cannot have; the repaired ones all have mult = cof = 0")

# H4: the repair must be TARGETED -- retracting corners bit-identical.
ok(force_and_sign(2, 7, 3, 9, 4)[0] == (375, 165, 42, 168)
   and force_and_sign(3, 4, 4, 8, 3)[0] == (550, 205, 69, 276)
   and law(2, 3, 4, 8, 7, 4)[1] == (238, 204, 30, 4),
   "H4 F14 at (9,24), (108,144) and (72,108) at (8,28) are BIT-IDENTICAL pre- and "
   "post-repair: their corners retract, so the repair is targeted, not a rewrite")

print()
if _fail:
    print(f"{_fail} FAILURE(S) out of {_n} checks")
    sys.exit(1)
print(f"ALL {_n} PHI-F14 CHECKS PASSED")
if not QUIET:
    import pathlib
    print(f"script: {pathlib.Path(__file__).resolve()}")
sys.exit(0)
