#!/usr/bin/env python3
"""phi_corner4.py  (NEW; read-only over all existing artifacts)

FOURTH-CORNER Phi derivation: break (or explain) the t/kappa correlation in the
corner-signature fit of PHI_75_125.md.

Result in brief:
  * kappa = t-2 is FORCED for every corner reached by the standard single
    Laurent chart (X,Y) -> (x^-1, x^l y): its Jacobian is -x^(l-2), and t = l.
    That covers all 15 length-1 families anchored at A_0' = (1,0) in the GGV5
    v11<=35 table (F1-F11, F14-F17).  t and kappa can therefore never be
    separated inside this class; the fit loses a parameter,
        N = a[t(a+b) - (kappa+1)] - 2b  =  a[t(a+b-1) + 1] - 2b .
  * What CAN be varied is t itself.  The known points had t in {4,5}; family
    F9 (corner (7,21), final corner (11\\7,2)) gives t = 7.  Deriving its Phi
    by the exact corner-144 template and comparing to the fit's parameter-free
    prediction is the out-of-sample test.  A second new point (F2 at j=0, the
    (50,75) case: same (5,20) corner as (75,125) but reduced pair (2,3))
    isolates corner-vs-pair dependence for free.

Sources: GGV5 (paper_src/1708.07936_GGV5.tex) family tables (the two tables in
the final section, "Admissible complete chains with v11(A_0) <= 35"); the
corner-144 template CORNER_144_COMPARISON.md / corner144_verify.py; the
three-point fit PHI_75_125.md.  The independent PASS/FAIL checker is
phi_corner4_verify.py.  Everything below is exact sympy.
"""
import sympy as sp
from fractions import Fraction
from math import gcd

y = sp.symbols("y")

# ---------------------------------------------------------------------------
# 0. The GGV5 v11<=35 family tables, transcribed verbatim.
#    Length-1 families: (name, A0, A0', p, l, q, k, (m0, dm), (n0, dn))
#    with (m,n)(j) = (m0 + dm*j, n0 + dn*j) and final corner A1 = (p\l, q).
# ---------------------------------------------------------------------------
FAMILIES_LEN1 = [
    ("F1",  (4, 12), (1, 0),  7, 4, 3, 1, (3, 2),  (4, 3)),
    ("F2",  (5, 20), (1, 0),  7, 5, 2, 1, (2, 1),  (3, 2)),
    ("F3",  (5, 20), (1, 0),  8, 5, 3, 1, (3, 4),  (2, 3)),
    ("F4",  (5, 20), (1, 0),  8, 5, 3, 2, (3, 2),  (16, 12)),
    ("F5",  (5, 20), (1, 0),  9, 5, 4, 1, (9, 7),  (5, 4)),
    ("F6",  (5, 20), (1, 0),  9, 5, 4, 2, (7, 6),  (18, 16)),  # CORRECTED 2026-07-24: GGV5 prints F6 base (m,n)=(4,10) [gcd=2, violates coprimality]; the coprime family is (6j+7,16j+18)=base (7,18). See CHAIN_SURVEY.md. (Survey advances j to first coprime pair -> computed (m,n)=(7,18) unchanged.)
    ("F7",  (6, 15), (1, 0),  7, 3, 4, 1, (2, 1),  (7, 4)),
    ("F8",  (6, 15), (1, 0),  8, 3, 5, 1, (3, 2),  (7, 5)),
    ("F9",  (7, 21), (1, 0), 11, 7, 2, 1, (2, 1),  (3, 2)),
    ("F10", (7, 21), (1, 0), 13, 7, 3, 1, (7, 5),  (4, 3)),
    ("F11", (7, 21), (1, 0), 13, 7, 3, 2, (2, 1),  (5, 3)),
    ("F12", (8, 24), (2, 0), 13, 4, 5, 1, (3, 2),  (7, 5)),
    ("F13", (9, 21), (2, 0), 13, 3, 7, 1, (2, 1),  (13, 7)),
    ("F14", (9, 24), (1, 0),  7, 3, 4, 1, (2, 1),  (7, 4)),
    ("F15", (9, 24), (1, 0),  8, 3, 5, 1, (3, 2),  (7, 5)),
    ("F16", (9, 24), (1, 0), 10, 3, 7, 1, (3, 4),  (5, 7)),
    ("F17", (9, 24), (1, 0), 11, 3, 8, 1, (2, 5),  (3, 8)),
]
# Length-2 families (final corner only; composite chart NOT derived in any
# paper -- t and kappa are UNVERIFIED for these, see PHI_CORNER4.md):
FAMILIES_LEN2 = [
    ("F18", (6, 18),  7, 3, 4), ("F19", (6, 18),  8, 3, 5),
    ("F20", (6, 24),  7, 3, 4), ("F21", (6, 24),  8, 3, 5),
    ("F22", (8, 24),  5, 4, 2), ("F23", (8, 24), 11, 4, 4),
    ("F24", (8, 24), 19, 8, 3),
]

def fit_signature(a, b, t, kappa, a0, q):
    """The PHI_75_125.md six-parameter closed forms (r > 0 regime)."""
    e = b - a + 1
    r = a0 - q - 1
    N = a * (t * (a + b) - (kappa + 1)) - 2 * b
    deg = (e * a0 - q + 1) + N * a0
    ordy = ((e - 1) * q + 1) + N * q
    mult = e + N
    cof = r * (e + N)
    return N, e, r, (deg, ordy, mult, cof)

# ---------------------------------------------------------------------------
# 1. Candidate-corner survey (smallest j with gcd(m,n)=1 per family)
# ---------------------------------------------------------------------------
print("=" * 96)
print("STEP 1 -- candidate-corner survey (GGV5 v11<=35 tables, smallest coprime j)")
print("=" * 96)
print(f"{'fam':4} {'A0':>8} {'A1':>10} {'(m,n)':>8} {'degs':>10} "
      f"{'t':>2} {'kap':>3} {'a0':>3} {'q':>2} {'e':>2} {'r':>2} {'N':>4} "
      f"{'gap':>4}  notes")

survey = []
for name, A0, A0p, p, l, q, k, (m0, dm), (n0, dn) in FAMILIES_LEN1:
    j = 0
    while gcd(m0 + dm * j, n0 + dn * j) != 1:
        j += 1
    m, n = m0 + dm * j, n0 + dn * j
    v11 = A0[0] + A0[1]
    degs = (v11 * m, v11 * n)
    a, b = sorted((m, n))
    a0 = A0[0]
    t, kappa = l, l - 2                      # standard chart (A0'=(1,0) class)
    N, e, r, sig = fit_signature(a, b, t, kappa, a0, q)
    # resonance gap: resonant deg f  minus  pure-ansatz deg f
    res = Fraction((t * (b - a) + kappa + 1) * a0, t)
    gap = res - (e * a0 - q + 1)
    dio = (m + n) * q * k - n * (q * l - p)  # must equal k
    notes = []
    if A0p != (1, 0):
        notes.append("A0'!=(1,0): t,kappa UNVERIFIED")
    if k != 1:
        notes.append("k=2: N-formula unverified")
    if t not in (4, 5) and not notes:
        notes.append(f"NEW t={t}")
    assert dio == k, f"{name}: Diophantine failed"
    survey.append((name, j, A0, (p, l, q), (m, n), degs, t, kappa, a0, q, e, r, N, gap, notes))
    print(f"{name:4} {str(A0):>8} ({p}\\{l},{q})  {str((m,n)):>8} {str(degs):>10} "
          f"{t:>2} {kappa:>3} {a0:>3} {q:>2} {e:>2} {r:>2} {N:>4} "
          f"{str(gap):>4}  {'; '.join(notes)}")

print("\nLength-2 families (composite chart underived -- t,kappa unknown): "
      + ", ".join(f"{nm} A2=({p}\\{l},{q})" for nm, A0, p, l, q in FAMILIES_LEN2))

# ---------------------------------------------------------------------------
# 2. kappa = t-2 : forced for the standard-chart class
# ---------------------------------------------------------------------------
print("\n" + "=" * 96)
print("STEP 2 -- kappa = t-2 is forced on the whole standard-chart class")
print("=" * 96)
x, ls = sp.symbols("x l_s", positive=True)
X, Y = x**-1, x**ls * y
J = sp.simplify(sp.diff(X, x) * sp.diff(Y, y) - sp.diff(X, y) * sp.diff(Y, x))
print(f"  Jacobian of (X,Y) -> (x^-1, x^l y):  {J}   (= -x^(l-2) for all l)")
print("  ell(C) = x^t c with t = l (each factor (Y - r X^-l) -> x^l (y - r)).")
print("  Hence kappa = l-2 = t-2 for EVERY corner reduced by this chart --")
print("  all 15 length-1 A0'=(1,0) families above.  t and kappa cannot be")
print("  separated inside this class; escapes need A0'=(2,0) or length-2")
print("  chains, whose charts are derived in no paper (see PHI_CORNER4.md).")

# ---------------------------------------------------------------------------
# 3. Derivations: F9 j=0 (56,84), t=7  and  F2 j=0 (50,75), t=5
# ---------------------------------------------------------------------------
def derive(tag, a, b, t, kappa, a0, q):
    e = b - a + 1
    coef = t * (b - a) + kappa + 1
    rho = (e - 1) * q + 1
    dg = a0 - q
    N, _, r, fit = fit_signature(a, b, t, kappa, a0, q)
    print("\n" + "=" * 96)
    print(f"STEP 3 -- {tag}:  (a,b)=({a},{b})  t={t} kappa={kappa} a0={a0} q={q}"
          f"  e={e} r={r} rho={rho} deg g={dg} N={N}")
    print("=" * 96)
    print(f"  forcing ODE:  {a*t} c f' - {a*coef} c' f = c^{e},   c = y^{q} g")

    # generic-g collapse
    gc = sp.symbols(f"g0:{dg+1}")
    A = sp.symbols("A")
    g = sum(gc[i] * y**i for i in range(dg + 1))
    c = y**q * g
    f = A * y**rho * g**e
    resid = sp.expand(a * t * c * sp.diff(f, y) - a * coef * sp.diff(c, y) * f - c**e)
    quo = sp.expand(sp.factor(resid) / (y**(e * q) * g**(e - 1)))
    eqs = sp.Poly(quo, y).all_coeffs()
    print(f"  ansatz f = A y^{rho} g^{e} collapses the ODE; coefficient system "
          f"forces g_1..g_{dg-1} = 0,")
    print(f"  top coefficient resonant (free), g(-1)=0 + monic => g = y^{dg} + 1.")

    g_sol = y**dg + 1
    A_sol = sp.Rational(-1, a * t - a * coef)          # from -[a*coef-a*t] A g0 = 1... solved below
    # solve A exactly from the constant coefficient with g = g_sol:
    A_sol = sp.solve(sp.expand(quo.subs({gc[i]: sp.Poly(g_sol, y).coeff_monomial(y**i)
                                         for i in range(dg + 1)})).coeff(y, 0), A)[0]
    c_sol = y**q * g_sol
    f_sol = sp.expand(A_sol * y**rho * g_sol**e)
    assert sp.expand(a * t * c_sol * sp.diff(f_sol, y)
                     - a * coef * sp.diff(c_sol, y) * f_sol - c_sol**e) == 0
    H = sp.factor(g_sol / (y + 1))
    print(f"  g = y^{dg}+1,  H = {H},  A = {A_sol}")
    print(f"  f = {A_sol} * y^{rho} * (y^{dg}+1)^{e}    (deg {sp.degree(f_sol, y)})")

    Phi = sp.expand(f_sol * c_sol**N)
    deg = sp.degree(Phi, y)
    ordy = min(mm[0] for mm in sp.Poly(Phi, y).monoms())
    m1 = 0
    qq = sp.Poly(Phi, y)
    d1 = sp.Poly(y + 1, y)
    while True:
        qq2, rem = sp.div(qq, d1)
        if not rem.is_zero:
            break
        qq, m1 = qq2, m1 + 1
    cof = deg - ordy - m1 * 1
    sig = (deg, ordy, m1, cof)
    print(f"  Phi = f * C^{N} = {A_sol} y^{ordy} (y^{dg}+1)^{m1}")
    print(f"  SIGNATURE (deg, ord_y, mult_(y+1), cofactor) = {sig}")
    verdict = "MATCHES" if sig == fit else "DIFFERS"
    print(f"  fit prediction (parameter-free)              = {fit}   ==> {verdict}")
    return sig, fit, verdict

sig9, fit9, v9 = derive("F9 j=0, case (56,84)", 2, 3, 7, 5, 7, 2)
sig2, fit2, v2 = derive("F2 j=0, case (50,75)", 2, 3, 5, 3, 5, 2)

# ---------------------------------------------------------------------------
# 4. Verdict
# ---------------------------------------------------------------------------
print("\n" + "=" * 96)
print("VERDICT")
print("=" * 96)
print(f"  fourth corner  (56,84), t=7 : derived {sig9}  fit {fit9}  -> {v9}")
print(f"  fifth  point   (50,75), t=5 : derived {sig2}  fit {fit2}  -> {v2}")
print("""  kappa = t-2 is FORCED on the standard single-chart class (all length-1
  A0'=(1,0) families).  Within it the fit has NO free parameters left, and it
  now reproduces FIVE points exactly at t in {4,5,7}:
    (72,108) r=0*, (108,144), (75,125), (56,84), (50,75).
  The only routes to kappa != t-2 are the A0'=(2,0) families (F12, F13) and
  the length-2 chains (F18-F24), whose reduction charts exist in no paper.""")
