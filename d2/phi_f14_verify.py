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

def law(a, b, t, a0, q, gap):
    """Unified law, kappa = t-2 eliminated (PHI_CORNER4 sec.2)."""
    e, r = b - a + 1, a0 - q - 1
    N = a * (t * (a + b - 1) + 1) - 2 * b
    rho = (e - 1) * q + 1
    return N, ((e * a0 - q + 1) + gap + N * a0, rho + N * q,
               e + N, gap + r * (e + N))

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

lem_gap, lem_dg = True, True
for name, A0, p, l, q, k, m0, dm, n0, dn in FAMS:
    j = 0
    while gcd(m0 + dm * j, n0 + dn * j) != 1:
        j += 1
    m, n = m0 + dm * j, n0 + dn * j
    a_, b_ = sorted((m, n))
    t_, kap_, a0_ = l, l - 2, A0[0]
    e_ = b_ - a_ + 1
    coef_ = t_ * (b_ - a_) + kap_ + 1
    res_ = Fraction(coef_ * a0_, t_)
    gap_ = res_ - (e_ * a0_ - q + 1)
    lem_gap &= (gap_ == Fraction(q - 1) - Fraction(a0_, t_))
    lem_dg &= (a0_ - q >= 1)
ok(lem_gap, "A: mini-lemma gap = (q-1) - a0/t on all 15 standard-chart families")
ok(lem_dg, "A: dg = a0-q >= 1 on all 15 (residual g = y^(a0-q)+1 well-defined)")
ok(9 == 3 * (4 - 1), "A: F14 satisfies a0 = t(q-1)  <=>  gap = 0 (resonance exact)")
ok(Fraction(3-1) - Fraction(4, 4) == 1, "A: F1 gap = (q-1) - a0/t = 1 > 0 (resonance broken)")

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

# --------------------------------------------- E. F1: gap>0 / r=0 fresh solve
a1, b1, t1, a01, q1 = 3, 4, 4, 4, 3
e1, r1 = b1 - a1 + 1, a01 - q1 - 1
coef1 = t1 * (b1 - a1) + (t1 - 2) + 1
rho1 = (e1 - 1) * q1 + 1
gap1 = Fraction(coef1 * a01, t1) - (e1 * a01 - q1 + 1)
ok(r1 == 0 and gap1 == 1, "E: F1 is r=0 with gap = 1")
c1 = y**q1 * (y + 1)
res_deg1 = e1 * a01 - q1 + 1 + int(gap1)
fv = sp.symbols(f"v0:{res_deg1 + 3}")
f_v = sum(fv[i] * y**i for i in range(res_deg1 + 3))
resid_v = sp.expand(a1 * t1 * c1 * sp.diff(f_v, y)
                    - a1 * coef1 * sp.diff(c1, y) * f_v - c1**e1)
sols_v = sp.solve(sp.Poly(resid_v, y).all_coeffs(), fv, dict=True)
ok(len(sols_v) == 1, "E: F1 ODE 12cf' - 21c'f = c^2 has a UNIQUE polynomial solution"
                     " (degree allowed 2 past resonant)")
f1s = sp.expand(f_v.subs(sols_v[0]))
ok(sp.degree(f1s, y) == res_deg1 == 7,
   "E: deg f = 7 = resonant degree (pure ansatz degree 6 is NOT attained)")
ok(sp.expand(f1s - sp.Rational(1, 15) * y**4 * (y + 1)**2 * (4 * y - 1)) == 0,
   "E: f = (1/15) y^4 (y+1)^2 (4y-1)")
u1 = sp.cancel(f1s / (y**rho1 * (y + 1)**e1))
ok(sp.degree(u1, y) == 1 and u1.subs(y, 0) != 0 and u1.subs(y, -1) != 0,
   "E: cofactor u = (4y-1)/15: deg = gap = 1, u(0) != 0, u(-1) != 0 (UNIT)")
N1, law1 = law(a1, b1, t1, a01, q1, int(gap1))
ok(N1 == 67, "E: N = 67 for F1 (same N as (108,144): same (a,b,t,q), different a0)")
sig1 = signature(f1s * c1**N1)
ok(sig1 == (275, 205, 69, 1) and sig1 == law1,
   f"E: F1 signature (275,205,69,1) = amended-law prediction {law1}  ==> MATCH")

# ------------------------------------------------------------------ F. controls
def force_and_sign(a_, b_, t_, a0_, q_):
    """Re-derive a gap=0 corner end to end; return signature."""
    e_, dg_ = b_ - a_ + 1, a0_ - q_
    coef_ = t_ * (b_ - a_) + (t_ - 2) + 1
    rho_ = (e_ - 1) * q_ + 1
    g_ = y**dg_ + 1
    c_ = y**q_ * g_
    A_ = sp.Rational(1, a_ * (t_ * rho_ - coef_ * q_))
    f_ = sp.expand(A_ * y**rho_ * g_**e_)
    _require(sp.expand(a_ * t_ * c_ * sp.diff(f_, y)
                     - a_ * coef_ * sp.diff(c_, y) * f_ - c_**e_) == 0, "sp.expand(a_ * t_ * c_ * sp.diff(f_, y) - a_ * coef_ * sp.diff(c_, y) * f_ - c_**e_) == 0")
    N_ = a_ * (t_ * (a_ + b_ - 1) + 1) - 2 * b_
    return signature(f_ * c_**N_), A_

CONTROLS = [  # (label, a, b, t, a0, q, audited signature, audited lc)
    ("(108,144)", 3, 4, 4, 8, 3, (550, 205, 69, 276), sp.Rational(-1, 15)),
    ("(75,125)",  3, 5, 5, 5, 2, (504, 201, 101, 202), sp.Rational(-1, 9)),
    ("(56,84)",   2, 3, 7, 7, 2, (377, 107, 54, 216), sp.Rational(-1, 10)),
    ("(50,75)",   2, 3, 5, 5, 2, (189, 75, 38, 76),  sp.Rational(-1, 6)),
]
for lbl, ca, cb, ct, ca0, cq, csig, clc in CONTROLS:
    s_, A_ = force_and_sign(ca, cb, ct, ca0, cq)
    ok(s_ == csig and A_ == clc,
       f"F: control {lbl}: re-derived signature {s_} and lc {A_} match audited values")
ok(law(2, 3, 4, 8, 7, 4)[1] == (238, 204, 30, 4),
   "F: (72,108) audited signature (238,204,30,4) obeys the amended law with gap = 4")

# ------------------------------------------------------------------- G. census
SEVEN = [("(72,108)", 2, 3, 4, 8, 7, 4, 28), ("(108,144)", 3, 4, 4, 8, 3, 0, 67),
         ("(75,125)", 3, 5, 5, 5, 2, 0, 98), ("(56,84)", 2, 3, 7, 7, 2, 0, 52),
         ("(50,75)", 2, 3, 5, 5, 2, 0, 36), ("(66,231)", 2, 7, 3, 9, 4, 0, 36),
         ("(48,64)", 3, 4, 4, 4, 3, 1, 67)]
ok(all(law(sa, sb, st, sa0, sq, sg)[0] == sN
       for _, sa, sb, st, sa0, sq, sg, sN in SEVEN),
   "G: reduced N-formula N = a[t(a+b-1)+1] - 2b reproduces N at all SEVEN points")
ok(sorted({st for _, sa, sb, st, sa0, sq, sg, sN in SEVEN}) == [3, 4, 5, 7]
   and sorted({sb - sa + 1 for _, sa, sb, st, sa0, sq, sg, sN in SEVEN}) == [2, 3, 6],
   "G: coverage census t in {3,4,5,7}, e in {2,3,6}")

print()
if _fail:
    print(f"{_fail} FAILURE(S) out of {_n} checks")
    sys.exit(1)
print(f"ALL {_n} PHI-F14 CHECKS PASSED")
if not QUIET:
    import pathlib
    print(f"script: {pathlib.Path(__file__).resolve()}")
sys.exit(0)
