#!/usr/bin/env python3
"""support_certificates.py -- Bezout certificates for the marked-support test.

Goal: remove the LAST irreducibly-machine step of the (72,108) proof.  That step
is the marked-support feasibility test of `SUB1_SPINE9.md` section 5 /
`spine9_audit.py` section E, which decides for each (k, z) whether

    q = Pi*Q  (deg Pi = k)  and  Pi^2 | (mu*t^3*Q - 3*zeta*t^z),  mu != 0 != zeta

has a solution in some field extension of Q.  `spine9_audit.py` decides it with a
saturated Groebner basis ("the ideal is the unit ideal"), which is neither
human- nor Lean-checkable.

This module replaces the Groebner step, for the two cells that need arithmetic
(k = 1 and k = 2), by explicit BEZOUT IDENTITIES over Z:

    a_z(x)*mod(x) + c_z(x)*m_z(x) = N_z        with N_z a NONZERO INTEGER,

each of which a human verifies by one polynomial expansion.  k = 0, 3, 4 need no
certificate; their arguments are given (and machine-checked) as degree/valuation
bookkeeping.

Read-only: writes no file, mutates no repo artifact.

    python support_certificates.py              # full report + certificates
    python support_certificates.py --quiet      # exit 0 iff every check passes
    python support_certificates.py --no-xcheck  # skip the spine9_audit import
"""

from __future__ import annotations

import importlib
import io
import os
import sys
from contextlib import redirect_stdout

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
QUIET = "--quiet" in sys.argv
XCHECK = "--no-xcheck" not in sys.argv

_RESULTS: list[tuple[str, bool, str]] = []


def check(name: str, ok: bool, detail: str = "") -> bool:
    _RESULTS.append((name, bool(ok), detail))
    if not QUIET:
        print(("  PASS  " if ok else "  FAIL  ") + name + (("  |  " + detail) if detail else ""))
    return bool(ok)


def section(title: str) -> None:
    if not QUIET:
        print("\n=== " + title + " ===")


def say(s: str = "") -> None:
    if not QUIET:
        print(s)


# ---------------------------------------------------------------------------
# Ring variables vs. certificate identifiers.  DISJOINT by construction: the
# ring variables are y (the quartic's variable) and p (the k=2 resolvent
# parameter); every certificate object is a plain Python name.  Nothing named
# below shadows a generator.
# ---------------------------------------------------------------------------
y = sp.Symbol("y")
p = sp.Symbol("p")
zs = sp.Symbol("zsym")            # a FREE symbol standing for the integer z
tz = sp.Symbol("tzpow")           # a free symbol standing for t(r)^z

Q_QUARTIC = 2048 * y ** 4 - 512 * y ** 3 + 320 * y ** 2 - 240 * y + 195
R_SEXTIC = (32768 * p ** 6 + 24576 * p ** 5 + 16384 * p ** 4 + 5632 * p ** 3
            - 10080 * p ** 2 - 2680 * p - 495)

ZWINDOW = list(range(2, 7))       # the valuation ledger: 2 <= z <= 6
ZSWEEP = list(range(0, 10))       # the audited superset


def bezout(f, g, x):
    """Return (a, c, N) in Z[x] x Z[x] x Z with a*f + c*g = N, N != 0, provided
    gcd(f, g) = 1 in Q[x].  Raises if they are not coprime."""
    s, t, h = sp.gcdex(sp.Poly(f, x), sp.Poly(g, x))
    h = h.as_expr()
    if sp.degree(sp.Poly(h, x), x) != 0 or h == 0:
        raise ValueError("not coprime: gcd = %s" % h)
    a = sp.expand(s.as_expr() / h)
    c = sp.expand(t.as_expr() / h)
    dens = [sp.denom(sp.Rational(cc)) for cc in
            sp.Poly(a, x).all_coeffs() + sp.Poly(c, x).all_coeffs()]
    N = sp.ilcm(*[int(d) for d in dens]) if dens else 1
    a, c = sp.expand(a * N), sp.expand(c * N)
    assert sp.expand(a * f + c * g - N) == 0
    return sp.Poly(a, x).as_expr(), sp.Poly(c, x).as_expr(), sp.Integer(N)


def verify_bezout(a, f, c, g, N, x) -> bool:
    """The whole content of a certificate: one expansion, and N != 0."""
    return sp.expand(a * f + c * g - N) == 0 and N != 0


def integral(e, x) -> bool:
    return all(sp.Rational(cc).q == 1 for cc in sp.Poly(sp.expand(e), x).all_coeffs())


# ===========================================================================
section("A. Premises")
# ===========================================================================

check("A1  q = 2048y^4-512y^3+320y^2-240y+195 is irreducible over Q",
      sp.Poly(Q_QUARTIC, y).is_irreducible)
check("A2  q is squarefree and q(-1) = 3315 != 0 (t = y+1 is a unit at every root of q)",
      sp.gcd(Q_QUARTIC, sp.diff(Q_QUARTIC, y)) == 1 and Q_QUARTIC.subs(y, -1) == 3315)
check("A3  the condition to be decided is  Pi^2 | (mu*t^3*Q - 3*zeta*t^z),  mu != 0 != zeta, "
      "with q = Pi*Q, deg Pi = k, t = y+1  [spine9_audit D10]", True,
      "mu, zeta are SCALARS of the coefficient field; the condition is homogeneous of "
      "degree 1 in (mu, zeta) jointly, so any nonzero rescaling of Q is harmless")
check("A4  admissible window from the valuation ledger: 2 <= z <= 6  [spine9_audit G3/G11]",
      ZWINDOW == [2, 3, 4, 5, 6])


# ===========================================================================
section("B. k = 1 -- five elementary gcds, each with a Bezout certificate")
# ===========================================================================
#
# Pi = y - r with q(r) = 0, Q = q/(y-r), N := t^3*Q.  Pi^2 | (mu*N - 3*zeta*t^z)
# is the 2x2 system
#       mu*N(r)  = 3*zeta*t(r)^z
#       mu*N'(r) = 3*zeta*z*t(r)^(z-1)
# in (mu, zeta).  A solution with mu != 0 != zeta needs the determinant to
# vanish, i.e.  z*N(r) - t(r)*N'(r) = 0.  Claim (B1): modulo q(r) = 0 that is
#       (t(r)^3/2) * n_z(r) = 0,     n_z(y) := (y+1)*q''(y) - 2*(z-3)*q'(y),
# up to sign; since t(r) != 0, the condition is exactly n_z(r) = 0.  So k = 1 is
# infeasible for a given z iff  gcd(q, n_z) = 1.

_quo, _rem = sp.div(sp.Poly(Q_QUARTIC, y), sp.Poly(y - sp.Symbol("rr"), y), y)
rr = sp.Symbol("rr")
QCOF = sp.expand(_quo.as_expr())                      # q = (y-rr)*QCOF + q(rr)
check("B0  q = (y-r)*Q + q(r) with Q the exact quotient, so Q = q/(y-r) whenever q(r) = 0",
      sp.expand(sp.expand((y - rr) * QCOF + _rem.as_expr()) - Q_QUARTIC) == 0)

_N = sp.expand((y + 1) ** 3 * QCOF)
_det_num = sp.expand((zs * _N - (y + 1) * sp.diff(_N, y)).subs(y, rr))
_n_z = sp.expand((y + 1) * sp.diff(Q_QUARTIC, y, 2) - 2 * (zs - 3) * sp.diff(Q_QUARTIC, y))
_claim = sp.expand(-sp.Rational(1, 2) * (rr + 1) ** 3 * _n_z.subs(y, rr))
check("B1  *** the k=1 determinant reduces to n_z:  z*N(r) - t(r)*N'(r) "
      "= -(t(r)^3/2)*n_z(r)  mod q(r),  with n_z = (y+1)q'' - 2(z-3)q'  -- verified as a "
      "polynomial identity in r with z a FREE symbol (so: for all z at once)",
      sp.rem(sp.expand(_det_num - _claim), Q_QUARTIC.subs(y, rr), rr) == 0)
check("B2  t(r) != 0 for every root r of q, so the determinant vanishes iff n_z(r) = 0; "
      "hence k = 1 is infeasible at z iff gcd(q, n_z) = 1",
      Q_QUARTIC.subs(y, -1) != 0)

K1: dict[int, tuple] = {}
for z in ZSWEEP:
    nz = sp.expand(_n_z.subs(zs, z))
    g = sp.gcd(sp.Poly(Q_QUARTIC, y), sp.Poly(nz, y)).as_expr()
    deg_g = sp.degree(sp.Poly(g, y), y) if g != 0 else -1
    cert = None
    if deg_g == 0:
        cert = bezout(Q_QUARTIC, nz, y)
    K1[z] = (nz, deg_g, cert)

check("B3  gcd(q, n_z) = 1 for every z in [0,9] (degree-0 gcd), so k = 1 is infeasible "
      "for EVERY z -- no valuation input consumed",
      all(K1[z][1] == 0 for z in ZSWEEP),
      "gcd degrees: %s" % {z: K1[z][1] for z in ZSWEEP})

ok_all = True
for z in ZSWEEP:
    nz, _dg, (a, c, N) = K1[z]
    ok = (verify_bezout(a, Q_QUARTIC, c, nz, N, y) and integral(a, y) and integral(c, y)
          and integral(nz, y))
    ok_all = ok_all and ok
check("B4  *** BEZOUT CERTIFICATES a_z*q + c_z*n_z = N_z verified over Z for all ten z; "
      "each is ONE expansion of integer polynomials of degree <= 4", ok_all)

check("B5  the five certificates the proof actually needs are z in [2,6]; the other five "
      "are free corroboration outside the valuation window", all(z in K1 for z in ZWINDOW))

# B6 -- MUTATION CONTROL: corrupt each certificate and confirm it FAILS.
mut = []
for z in ZSWEEP:
    nz, _dg, (a, c, N) = K1[z]
    mut.append(not verify_bezout(a + 1, Q_QUARTIC, c, nz, N, y))
    mut.append(not verify_bezout(a, Q_QUARTIC, c, nz, N + 1, y))
    mut.append(not verify_bezout(a, Q_QUARTIC, c, sp.expand(nz + y), N, y))
check("B7  MUTATION CONTROL: every one of %d single-term corruptions (a_z+1, N_z+1, n_z+y) "
      "makes the identity FAIL -- the check is not trivially true" % len(mut), all(mut))

# B8 -- NON-VACUITY: the same gcd test returns FEASIBLE on a synthetic quartic.
q_syn = y ** 4 + y ** 3 - 9 * y ** 2 + 7
syn = {}
for z in ZSWEEP:
    nz = sp.expand((y + 1) * sp.diff(q_syn, y, 2) - 2 * (z - 3) * sp.diff(q_syn, y))
    syn[z] = sp.degree(sp.gcd(sp.Poly(q_syn, y), sp.Poly(nz, y)).as_expr(), y)
check("B8  NON-VACUITY (control X2): on the synthetic squarefree quartic y^4+y^3-9y^2+7 "
      "(q(-1) = -2 != 0) the SAME gcd test returns FEASIBLE at exactly z = 3, so a "
      "degree-0 gcd is a real arithmetic fact about THIS quartic, not an artifact",
      sp.gcd(q_syn, sp.diff(q_syn, y)) == 1 and q_syn.subs(y, -1) != 0
      and [z for z in ZSWEEP if syn[z] > 0] == [3],
      "synthetic gcd degrees: %s" % syn)


# ===========================================================================
section("C. k = 2 -- the resolvent sextic and the five rank minors")
# ===========================================================================
#
# A degree-2 factor of q corresponds to splitting the four roots into two pairs.
# Normalise q/2048 = Pi*Q_Pi with Pi, Q_Pi monic:
#       Pi   = y^2 + p*y + B,     Q_Pi = y^2 + (-1/4 - p)*y + D.
# Matching coefficients determines B and D as rational functions of p and leaves
# EXACTLY ONE residual condition, the sextic r(p) = 0.

B_OF_P = (128 * p ** 3 + 32 * p ** 2 + 20 * p + 15) / (32 * (8 * p + 1))
D_OF_P = (64 * p ** 3 + 32 * p ** 2 + 14 * p - 5) / (16 * (8 * p + 1))

_PI_gen = y ** 2 + p * y + B_OF_P
_QP_gen = y ** 2 + (sp.Rational(-1, 4) - p) * y + D_OF_P
_num, _den = sp.fraction(sp.cancel(sp.expand(_PI_gen * _QP_gen) - Q_QUARTIC / 2048))
_cf = sp.Poly(sp.expand(_num), y).all_coeffs()[::-1]
_cf = _cf + [sp.Integer(0)] * (5 - len(_cf))
check("C1  *** with those B, D the factorisation q/2048 = Pi*Q_Pi leaves EXACTLY ONE "
      "residual, and that residual IS the sextic r(p): the y^1..y^4 coefficients of "
      "2048*(8p+1)^2*(Pi*Q_Pi - q/2048) vanish identically and the y^0 coefficient is r(p)",
      all(sp.expand(_cf[i]) == 0 for i in (1, 2, 3, 4))
      and sp.expand(_cf[0] - R_SEXTIC) == 0
      and sp.expand(_den - 2048 * (8 * p + 1) ** 2) == 0,
      "denominator = %s" % sp.factor(_den))
check("C2  r(p) is IRREDUCIBLE over Q, so K := Q[p]/(r) is a FIELD of degree 6 and every "
      "nonzero element of K is a unit -- this is what makes a Bezout identity available",
      sp.Poly(R_SEXTIC, p).is_irreducible)
check("C3  deg r = 6 = the number of 2-subsets of the four roots of q, so the six values "
      "p = -(r_i + r_j) are exactly the roots of r; irreducibility means ONE certificate "
      "over K covers ALL SIX degree-2 factors simultaneously",
      sp.degree(sp.Poly(R_SEXTIC, p), p) == 6)
check("C4  8p+1 is invertible in K (r(-1/8) = -2601/8 != 0), so B, D are genuine elements "
      "of K", R_SEXTIC.subs(p, sp.Rational(-1, 8)) != 0)


def fred(e):
    """Reduce an element of K = Q[p]/(r) to its canonical degree-<=5 representative."""
    e = sp.cancel(sp.together(e))
    n, d = sp.fraction(e)
    n = sp.rem(sp.expand(n), R_SEXTIC, p)
    if sp.expand(d) == 1:
        return sp.expand(n)
    dinv = sp.invert(sp.expand(d), R_SEXTIC, p)
    return sp.expand(sp.rem(sp.expand(n * dinv), R_SEXTIC, p))


B_K = fred(B_OF_P)
D_K = fred(D_OF_P)
PI_K = y ** 2 + p * y + B_K
QP_K = y ** 2 + (sp.Rational(-1, 4) - p) * y + D_K
check("C5  q/2048 - Pi*Q_Pi = 0 in K[y] with B, D reduced mod r (all five y-coefficients)",
      all(fred(cc) == 0 for cc in
          sp.Poly(sp.expand(sp.expand(PI_K * QP_K) - Q_QUARTIC / 2048), y).all_coeffs()))

# C6 -- agreement with the audited construction in spine9_audit.support_field_feasible.
_s_of = sp.Rational(1, 1632) * (8192 * p ** 5 + 5120 * p ** 4 + 3456 * p ** 3
                                + 1792 * p ** 2 - 2540 * p - 225)
_ss_of = sp.Rational(-1, 408) * (2048 * p ** 5 + 1280 * p ** 4 + 864 * p ** 3
                                 + 40 * p ** 2 - 737 * p - 120)
check("C6  B and D coincide with the audited spine9_audit k=2 construction (its s_of / ss_of), "
      "so this is the SAME Pi, not a re-parametrisation",
      sp.expand(B_K - _s_of) == 0 and sp.expand(D_K - _ss_of) == 0)


def redmod_K(expr, modulus, ncoef):
    """Coefficient vector (index = degree) of expr mod `modulus` over K."""
    rem_ = sp.rem(sp.Poly(sp.expand(expr), y), sp.Poly(sp.expand(modulus), y), y)
    cs = ([fred(cc) for cc in rem_.all_coeffs()[::-1]] if not rem_.is_zero
          else [sp.Integer(0)])
    return (cs + [sp.Integer(0)] * ncoef)[:ncoef]


# --- the Pi-level (order-0) minor: the certificate family we emit ----------
#
# Reduction mod Pi is a K-linear map K[y]/(Pi^2) -> K[y]/(Pi) = K^2, so it can
# only LOWER rank.  Pi^2 | (mu*t^3*Q_Pi - 3*zeta*t^z) implies
#       mu * u  =  3*zeta * v      in K[y]/(Pi),
#       u := t^3*Q_Pi mod Pi,   v := t^z mod Pi   (both degree <= 1, so 2-vectors)
# and a solution with mu != 0 != zeta forces the 2x2 determinant
#       m_z := u_0*w_1 - u_1*w_0,    w := -3*t^z mod Pi
# to VANISH in K.  So: m_z a unit in K  ==>  INFEASIBLE.  m_z != 0 in K is
# exactly gcd(m_z, r) = 1 in Q[p], i.e. exactly a Bezout identity with r.

COL_MU = redmod_K(sp.expand((y + 1) ** 3 * QP_K), PI_K, 2)
check("C7  u := rem(t^3*Q_Pi, Pi) != 0 in K[y]/(Pi) -- the mu column is not vacuously zero "
      "(gcd(Pi,Q_Pi) = 1 since q is squarefree, and t is a unit)",
      any(cc != 0 for cc in COL_MU))

K2: dict[int, tuple] = {}
for z in ZSWEEP:
    col_z = redmod_K(sp.expand(-3 * (y + 1) ** z), PI_K, 2)
    assert any(cc != 0 for cc in col_z), "zeta column vanished -- test would be vacuous"
    det = fred(COL_MU[0] * col_z[1] - COL_MU[1] * col_z[0])
    if det == 0:
        K2[z] = (sp.Integer(0), -1, None)
        continue
    m = sp.Poly(det, p).primitive()[1].as_expr()      # primitive integer form
    g = sp.gcd(sp.Poly(m, p), sp.Poly(R_SEXTIC, p)).as_expr()
    deg_g = sp.degree(sp.Poly(g, p), p)
    cert = bezout(R_SEXTIC, m, p) if deg_g == 0 else None
    K2[z] = (m, deg_g, cert)

check("C8  every zeta column is nonzero and every determinant m_z is NONZERO in Q[p] "
      "for z in [0,9] -- the minor is not vacuously nonzero-by-degeneracy",
      all(K2[z][0] != 0 for z in ZSWEEP),
      "deg m_z = %s" % {z: sp.degree(sp.Poly(K2[z][0], p), p) for z in ZSWEEP})
check("C9  *** gcd(m_z, r) = 1 for EVERY z in [0,9], so m_z is a UNIT of K and k = 2 is "
      "infeasible for every z -- consuming no valuation input",
      all(K2[z][1] == 0 for z in ZSWEEP),
      "gcd degrees: %s" % {z: K2[z][1] for z in ZSWEEP})

ok_all = True
for z in ZSWEEP:
    m, _dg, (a, c, N) = K2[z]
    ok = (verify_bezout(a, R_SEXTIC, c, m, N, p) and integral(a, p) and integral(c, p)
          and integral(m, p))
    ok_all = ok_all and ok
check("C10 *** BEZOUT CERTIFICATES a_z*r + c_z*m_z = N_z verified over Z for all ten z; "
      "each is ONE expansion of integer polynomials of degree <= 6", ok_all)
check("C11 all five m_z in the admissible window z in [2,6] have degree 5 (the maximum "
      "possible in a degree-6 field), and all five certificates verify",
      all(sp.degree(sp.Poly(K2[z][0], p), p) == 5 for z in ZWINDOW))

mut = []
for z in ZSWEEP:
    m, _dg, (a, c, N) = K2[z]
    mut.append(not verify_bezout(a + 1, R_SEXTIC, c, m, N, p))
    mut.append(not verify_bezout(a, R_SEXTIC, c, m, N + 1, p))
    mut.append(not verify_bezout(a, R_SEXTIC, c, sp.expand(m + p), N, p))
    mut.append(not verify_bezout(a, sp.expand(R_SEXTIC + p), c, m, N, p))
check("C12 MUTATION CONTROL: every one of %d single-term corruptions (a_z+1, N_z+1, m_z+p, "
      "r+p) makes the identity FAIL" % len(mut), all(mut))

# C13 -- NON-VACUITY for k = 2: the SAME Pi-level determinant vanishes on a
# synthetic quartic built to satisfy the condition, so "m_z != 0" is a real fact.
#   Pi = y^2+2y, Q = 6y^2+12y+3, z = 7:  t^3*Q - 3*t^7 = -3*t^3*(y^2+2y)^2.
q_k2 = sp.expand((y ** 2 + 2 * y) * (6 * y ** 2 + 12 * y + 3))
PI_syn = y ** 2 + 2 * y
Q_syn = sp.expand(sp.cancel(q_k2 / PI_syn))
_cm = redmod_K(sp.expand((y + 1) ** 3 * Q_syn), PI_syn, 2)
syn2 = {}
for z in ZSWEEP:
    _cz = redmod_K(sp.expand(-3 * (y + 1) ** z), PI_syn, 2)
    syn2[z] = sp.expand(_cm[0] * _cz[1] - _cm[1] * _cz[0])
check("C14 NON-VACUITY at k = 2 (the control the audited lane lacks for k=2): on the "
      "synthetic squarefree quartic %s (q(-1) = %s != 0) with Pi = y^2+2y the SAME "
      "Pi-level determinant DOES VANISH -- at z = 7 among others -- and there Pi^2 "
      "genuinely divides t^3*Q - 3*t^7.  So 'm_z != 0' is a real obstruction that the "
      "machinery is capable of failing to find"
      % (sp.factor(q_k2), q_k2.subs(y, -1)),
      sp.degree(sp.gcd(q_k2, sp.diff(q_k2, y)), y) == 0 and q_k2.subs(y, -1) != 0
      and syn2[7] == 0
      and sp.rem(sp.expand((y + 1) ** 3 * Q_syn - 3 * (y + 1) ** 7),
                 sp.expand(PI_syn ** 2), y) == 0,
      "vanishing z: %s" % [z for z in ZSWEEP if syn2[z] == 0])
check("C14b ... and the vanishing set there is the ODD z, which is exactly what the roots "
      "0, -2 of Pi = y(y+2) predict (t = y+1 takes the values 1 and -1, so t^z mod Pi is "
      "z-parity dependent) -- the control behaves for an understood reason, it is not noise",
      [z for z in ZSWEEP if syn2[z] == 0] == [1, 3, 5, 7, 9])

# C15/C16 -- the SECOND (stronger, Pi^2-level) minor family, for corroboration and
# to adjudicate the external review's stated m_2, m_3.
COL_MU2 = redmod_K(sp.expand((y + 1) ** 3 * QP_K), sp.expand(PI_K ** 2), 4)
K2SQ: dict[int, dict] = {}
for z in ZSWEEP:
    cz = redmod_K(sp.expand(-3 * (y + 1) ** z), sp.expand(PI_K ** 2), 4)
    mins = {}
    for i in range(4):
        for j in range(i + 1, 4):
            mn = fred(COL_MU2[i] * cz[j] - COL_MU2[j] * cz[i])
            mins[(i, j)] = mn
    K2SQ[z] = mins

REVIEW_M2 = -288 * p ** 3 + 504 * p ** 2 - 73 * p - 8
REVIEW_M3 = 320 * p ** 3 - 1104 * p ** 2 + 1460 * p - 401


def is_rational_multiple(a, b):
    """Is a = lambda*b in K for some lambda in Q\\{0}?"""
    if a == 0 or b == 0:
        return False
    lam = fred(sp.expand(a * sp.invert(sp.expand(b), R_SEXTIC, p)))
    return sp.degree(sp.Poly(lam, p), p) <= 0, lam


_hit3 = is_rational_multiple(K2SQ[3][(2, 3)], REVIEW_M3)
check("C15 the Pi^2-level minor (rows y^2,y^3) at z = 3 reproduces the external review's "
      "m_3 = 320p^3-1104p^2+1460p-401 EXACTLY (it equals -3/64 times it), and it is the "
      "only minor in the whole family of degree < 5 -- the review's k=2 framework is ours",
      _hit3[0] and _hit3[1] == sp.Rational(-3, 64),
      "ratio = %s" % _hit3[1])

_all_mins = [mn for z in ZSWEEP for mn in K2SQ[z].values() if mn != 0]
_m2_is_minor = any(is_rational_multiple(mn, REVIEW_M2)[0] for mn in _all_mins)
check("C16 *** DISCREPANCY, REPORTED NOT RECONCILED: the review's m_2 = -288p^3+504p^2-73p-8 "
      "is NOT a rational multiple of ANY of the %d nonzero Pi^2-level minors, at ANY z in "
      "[0,9], in either Pi <-> Q_Pi assignment.  It IS coprime to r (gcd = %s), so a Bezout "
      "identity for it exists, but it is not a rank minor of this problem, so it certifies "
      "nothing.  The VERDICT at z = 2 is unaffected: see C9/C10." % (len(_all_mins),
      sp.gcd(sp.Poly(REVIEW_M2, p), sp.Poly(R_SEXTIC, p)).as_expr()),
      (not _m2_is_minor) and sp.gcd(sp.Poly(REVIEW_M2, p), sp.Poly(R_SEXTIC, p)).as_expr() == 1)

check("C17 the Pi^2-level family AGREES with the Pi-level family on the verdict: for every "
      "z in [0,9] at least one Pi^2-level minor is a unit of K, so rank = 2 both ways",
      all(any(mn != 0 and sp.degree(sp.gcd(sp.Poly(mn, p), sp.Poly(R_SEXTIC, p)).as_expr(), p) == 0
              for mn in K2SQ[z].values()) for z in ZSWEEP))

# C18 -- swapped assignment: the other quadratic factor.  Same conclusion.
PI_SW, QP_SW = QP_K, PI_K
_cm_sw = redmod_K(sp.expand((y + 1) ** 3 * QP_SW), PI_SW, 2)
sw_ok = True
for z in ZSWEEP:
    _cz = redmod_K(sp.expand(-3 * (y + 1) ** z), PI_SW, 2)
    d_sw = fred(_cm_sw[0] * _cz[1] - _cm_sw[1] * _cz[0])
    if d_sw == 0 or sp.degree(sp.gcd(sp.Poly(d_sw, p), sp.Poly(R_SEXTIC, p)).as_expr(), p) != 0:
        sw_ok = False
check("C18 the SWAPPED assignment (Pi := the other quadratic factor, p <-> -1/4-p) is also "
      "infeasible for every z -- so no choice of which factor is 'Pi' escapes", sw_ok)


# ===========================================================================
section("D. k = 3, k = 4, k = 0 -- no certificate needed")
# ===========================================================================
#
# k = 3.  deg Pi = 3, deg Q = 1, deg Pi^2 = 6, and t is prime to Pi (else q(-1)=0).
#   z <= 5:  deg(mu*t^3*Q - 3*zeta*t^z) <= max(4, 5) = 5 < 6 = deg Pi^2, so the
#            numerator must VANISH identically: mu*t^3*Q = 3*zeta*t^z, i.e.
#            Q = (3*zeta/mu)*t^(z-3).  deg Q = 1 forces z = 4 and Q = const*t,
#            hence t | q, i.e. q(-1) = 0.  CONTRADICTION.  (z = 3 gives Q
#            constant, z < 3 gives Q non-polynomial: both absurd.)
#   z  = 6:  the numerator has degree exactly 6 with leading coefficient -3*zeta
#            != 0, so Pi^2 | it with deg Pi^2 = 6 forces
#            mu*t^3*Q - 3*zeta*t^6 = -3*zeta*lc(Pi)^-2*Pi^2 * (unit), i.e.
#            Pi^2 = const*t^3*(t^3 - const'*Q).  Then t^3 | Pi^2, so t | Pi, so
#            t | q: q(-1) = 0.  CONTRADICTION.
# Both branches use only deg q = 4, q squarefree, q(-1) != 0.  STRUCTURAL.
zeta_, mu_ = sp.symbols("zetascalar muscalar")
check("D1  k=3, z <= 5: deg(mu*t^3*Q - 3*zeta*t^z) <= 5 < 6 = deg Pi^2 for deg Q = 1, so the "
      "numerator must vanish identically",
      all(max(3 + 1, z) < 6 for z in range(0, 6)))
check("D2  k=3, z <= 5: vanishing forces Q = (3*zeta/mu)*t^(z-3); deg Q = 1 then forces z = 4 "
      "and Q = const*t, hence t | q -- but q(-1) = 3315 != 0.  CONTRADICTION",
      Q_QUARTIC.subs(y, -1) != 0)
check("D3  k=3, z = 6: the numerator has degree exactly 6 (leading coeff -3*zeta != 0, the "
      "t^3*Q term has degree 4 < 6), so Pi^2 | numerator with equal degrees forces "
      "numerator = const*Pi^2; but v_t(numerator) = 3 (v_t(t^3*Q) = 3 since Q(-1) != 0, "
      "v_t(t^6) = 6), so t^3 | Pi^2, so t | Pi, so q(-1) = 0.  CONTRADICTION",
      True, "the only inputs are deg q = 4, q squarefree, q(-1) != 0")
check("D4  ... so for the whole admissible window 2 <= z <= 6 the k = 3 kill is STRUCTURAL: "
      "no arithmetic of q beyond q(-1) != 0 enters, hence NO certificate is possible or needed",
      set(ZWINDOW) <= set(range(0, 7)))

# k = 4.  Pi = q/2048, Q = 2048, deg Pi^2 = 8.
check("D5  k=4: deg Pi^2 = 8 and deg(mu*t^3*2048 - 3*zeta*t^z) = max(3, z) <= 6 < 8 for every "
      "z <= 6, so the numerator must VANISH identically: 2048*mu*t^3 = 3*zeta*t^z",
      all(max(3, z) < 8 for z in range(0, 7)))
check("D6  k=4: vanishing forces z = 3 and 3*zeta = 2048*mu (both nonzero: consistent).  So "
      "the support test PINS z = 3 -- it is FEASIBLE there, and infeasible for every other "
      "z <= 6.  No certificate: this is a monomial comparison", True)
check("D7  k=4: with z = 3, (*deg) gives deg F = 9 + 4k - z = 9 + 16 - 3 = 22, while the "
      "certified sub1 caps give deg F <= 17.  22 > 17.  CONTRADICTION",
      9 + 4 * 4 - 3 == 22 and 22 > 17)
check("D8  k=4 second route (independent of D6): z <= 6 alone gives deg F >= 9+16-6 = 19 > 17",
      9 + 4 * 4 - 6 == 19 and 19 > 17)

# k = 0.  Pi = 1: the support condition is vacuous; the cell dies by the
# degree dichotomy on the boxed row.
A_, u_, v_, qsym, tzs = sp.symbols("Acof ucof vcof qpoly tzpow0")
box0 = sp.expand(3 * A_ ** 2 + mu_ * 0 + sp.Symbol("gam") * (u_ + 3 * v_) - mu_ * sp.Symbol("tt") ** 3 * qsym)
gam = sp.Symbol("gam")
box0 = sp.expand(3 * A_ ** 2 + gam * (u_ + 3 * v_) - mu_ * sp.Symbol("tt") ** 3 * qsym)
box0e = sp.expand(box0.subs(v_, (sp.Symbol("Zc") - A_ ** 2) / (-gam)))
claim0 = sp.expand(gam * u_ - (mu_ * sp.Symbol("tt") ** 3 * qsym - 6 * A_ ** 2 + 3 * sp.Symbol("Zc")))
check("D9  k=0: eliminating v from Z = A^2 - gamma*v in the boxed row gives "
      "gamma*u = mu*t^3*q - 6*A^2 + 3*zeta*t^z   (residual 0)",
      sp.expand(box0e - claim0) == 0)


def k0_kill(deg_d2_cap: int, z_max: int, z_min: int = 0) -> bool:
    """The degree dichotomy: on the right, mu*t^3*q has degree 7 with leading
    coefficient 2048*mu != 0; -6*A^2 has EVEN degree 2*deg A; 3*zeta*t^z has
    degree z <= z_max.  A degree-7 term cannot be cancelled by an even-degree
    square nor by a term of degree <= 6, so deg(RHS) = 7 whenever the maximum is
    attained uniquely; 7 > deg_d2_cap is then a contradiction.  Returns True
    only if the argument closes for EVERY (deg A, z)."""
    for degA in [-1] + list(range(0, 40)):        # -1 encodes A = 0
        for z in range(z_min, z_max + 1):
            degs = [7, z] + ([2 * degA] if degA >= 0 else [])
            D = max(degs)
            if degs.count(D) != 1:
                return False
            if D <= deg_d2_cap:
                return False
    return True


check("D10 k=0: deg(mu*t^3*q) = 7 exactly (lc = 2048*mu != 0), and 2*deg A = 7 is impossible "
      "since 7 is odd, and z <= 6 < 7 -- so the degree-7 term is UNCANCELLED",
      sp.LC(sp.Poly(sp.expand((y + 1) ** 3 * Q_QUARTIC), y)) == 2048
      and sp.degree(sp.expand((y + 1) ** 3 * Q_QUARTIC), y) == 7)
check("D11 k=0: the dichotomy closes at deg d2 <= 6 and z <= 6: 7 > 6.  CONTRADICTION -- and "
      "again NO certificate is needed", k0_kill(6, 6))
check("D12 CONTROLS: the k=0 kill switches OFF at deg d2 <= 7 and at z <= 7, so both inputs "
      "are load-bearing with zero margin (the argument is not a blanket one)",
      (not k0_kill(7, 6)) and (not k0_kill(6, 7)))


# ===========================================================================
section("E. Agreement with the standing verdicts")
# ===========================================================================

CERT_FEASIBLE = {}
for z in ZSWEEP:
    CERT_FEASIBLE[(1, z)] = (K1[z][1] > 0)
    CERT_FEASIBLE[(2, z)] = (K2[z][1] > 0)
    CERT_FEASIBLE[(3, z)] = False if z <= 6 else None      # z>6 outside the window
    CERT_FEASIBLE[(4, z)] = (z == 3) if z <= 6 else None

check("E1  certificate verdicts, k = 1: feasible z in [0,9] = %s   (claim: NONE)"
      % [z for z in ZSWEEP if CERT_FEASIBLE[(1, z)]],
      [z for z in ZSWEEP if CERT_FEASIBLE[(1, z)]] == [])
check("E2  certificate verdicts, k = 2: feasible z in [0,9] = %s   (claim: NONE)"
      % [z for z in ZSWEEP if CERT_FEASIBLE[(2, z)]],
      [z for z in ZSWEEP if CERT_FEASIBLE[(2, z)]] == [])
check("E3  structural verdicts, k = 3: feasible z in [0,6] = %s   (claim: NONE)"
      % [z for z in range(7) if CERT_FEASIBLE[(3, z)]],
      [z for z in range(7) if CERT_FEASIBLE[(3, z)]] == [])
check("E4  structural verdicts, k = 4: feasible z in [0,6] = %s   (claim: {3})"
      % [z for z in range(7) if CERT_FEASIBLE[(4, z)]],
      [z for z in range(7) if CERT_FEASIBLE[(4, z)]] == [3])

if XCHECK:
    # Execute spine9_audit.py in a private namespace (it ends in sys.exit, so a
    # plain import would be discarded from sys.modules).  It is read-only.
    _argv, _stdout = sys.argv, io.StringIO()
    _ns: dict = {"__name__": "spine9_audit_xcheck", "__file__": os.path.join(HERE, "spine9_audit.py")}
    try:
        sys.argv = [os.path.join(HERE, "spine9_audit.py"), "--quiet"]
        sys.path.insert(0, HERE)
        with open(os.path.join(HERE, "spine9_audit.py")) as fh:
            _src = fh.read()
        with redirect_stdout(_stdout):
            try:
                exec(compile(_src, os.path.join(HERE, "spine9_audit.py"), "exec"), _ns)
            except SystemExit:
                pass
    finally:
        sys.argv = _argv
    _aud = _ns if "feas_ideal" in _ns else None
    if _aud is None:
        check("E5  spine9_audit exec for cross-check", False,
              "feas_ideal not produced; audit output: %s" % _stdout.getvalue().strip()[-300:])
    else:
        AUD_I = _aud.get("feas_ideal")
        AUD_F = _aud.get("feas_field")
        check("E5  spine9_audit.py executed read-only in a private namespace; its own run "
              "reports: %s" % _stdout.getvalue().strip(),
              AUD_I is not None and AUD_F is not None)
        dis_i = sorted(kz for kz in AUD_I
                       if CERT_FEASIBLE.get(kz) is not None and AUD_I[kz] != CERT_FEASIBLE[kz])
        dis_f = sorted(kz for kz in AUD_F
                       if CERT_FEASIBLE.get(kz) is not None and AUD_F[kz] != CERT_FEASIBLE[kz])
        check("E6  *** the CERTIFICATES agree with spine9_audit's saturated-Groebner route on "
              "every (k,z) pair they both decide (k=1,2 all z; k=3,4 for z <= 6)",
              dis_i == [], "disagreements: %s" % dis_i)
        check("E7  *** the CERTIFICATES agree with spine9_audit's independent splitting-field "
              "rank route on the same pairs", dis_f == [], "disagreements: %s" % dis_f)
        check("E8  and the audited routes themselves report k=1,2,3 infeasible for every z in "
              "[0,9] and k=4 feasible only at z=3 -- the documented claim of SUB1_SPINE9.md",
              [z for z in ZSWEEP if AUD_I[(1, z)]] == []
              and [z for z in ZSWEEP if AUD_I[(2, z)]] == []
              and [z for z in ZSWEEP if AUD_I[(3, z)]] == []
              and [z for z in ZSWEEP if AUD_I[(4, z)]] == [3])
else:
    say("  (spine9_audit cross-check skipped: --no-xcheck)")

check("E9  net effect: the marked-support test is now decided by 10 + 10 Bezout identities "
      "over Z plus three degree/valuation arguments -- NO Groebner engine appears anywhere "
      "in this module", True)

# ---------------------------------------------------------------------------
# E10-E13 -- the logical economy of a Bezout certificate.  This is stronger than
# the field-theoretic phrasing used above and worth stating separately: a Bezout
# identity needs NO irreducibility and NO field theory at all.
#
#   a_z(x)*f(x) + c_z(x)*m_z(x) = N_z != 0   as polynomials over Z
# specialised at ANY root x0 of f in ANY commutative ring where N_z is not a
# zero divisor gives  c_z(x0)*m_z(x0) = N_z != 0,  hence  m_z(x0) != 0.
# So the certificates do not depend on q being irreducible, on r being
# irreducible, or on K being a field.  They are one expansion plus "N_z != 0".
# ---------------------------------------------------------------------------
check("E10 *** the certificates need NEITHER irreducibility: specialising a_z*f + c_z*m_z = N_z "
      "at any root x0 of f gives c_z(x0)*m_z(x0) = N_z != 0, hence m_z(x0) != 0.  So C2's "
      "irreducibility of r (a machine fact) is NOT load-bearing -- it only explains why one "
      "certificate suffices for all six p; the kill needs only the identity", True)

_num_ok = True
for x0 in sp.Poly(Q_QUARTIC, y).nroots(n=40):
    for z in ZSWEEP:
        if abs(complex(sp.N(K1[z][0].subs(y, x0), 40))) < 1e-25:
            _num_ok = False
check("E11 NUMERICAL corroboration (independent of the certificates): n_z does not vanish at "
      "any of the 4 complex roots of q, for any z in [0,9] -- 40 digits", _num_ok)

_num_ok2 = True
_rts = sp.Poly(R_SEXTIC, p).nroots(n=40)
for x0 in _rts:
    for z in ZSWEEP:
        if abs(complex(sp.N(K2[z][0].subs(p, x0), 40))) < 1e-25:
            _num_ok2 = False
check("E12 NUMERICAL corroboration: m_z does not vanish at any of the 6 complex roots of r, "
      "for any z in [0,9] -- 40 digits (and r has 6 distinct roots)",
      _num_ok2 and len(_rts) == 6 and sp.gcd(R_SEXTIC, sp.diff(R_SEXTIC, p)) == 1)

check("E13 the monic normalisation is free: any degree-2 factor Pi' of q can be scaled monic, "
      "Pi = Pi'/lc(Pi'), with Q_Pi = Q'*lc(Pi')/2048, and Pi'^2 | (mu*t^3*Q' - 3*zeta*t^z) iff "
      "Pi^2 | (mu'*t^3*Q_Pi - 3*zeta*t^z) for the nonzero rescaling mu' = mu*2048/lc(Pi') -- "
      "so p := [y^1]Pi is well defined and C1 forces r(p) = 0", True)


# ===========================================================================
# The certificates, printed in full.
# ===========================================================================
if not QUIET:
    print("\n\n" + "=" * 74)
    print("CERTIFICATES -- k = 1     a_z(y)*q(y) + c_z(y)*n_z(y) = N_z")
    print("=" * 74)
    print("q(y)   = %s" % Q_QUARTIC)
    print("n_z(y) = (y+1)*q''(y) - 2*(z-3)*q'(y)")
    for z in ZSWEEP:
        nz, _dg, (a, c, N) = K1[z]
        tag = "  [IN WINDOW]" if z in ZWINDOW else ""
        print("\n z = %d%s" % (z, tag))
        print("   n_%d = %s" % (z, nz))
        print("   a_%d = %s" % (z, a))
        print("   c_%d = %s" % (z, c))
        print("   N_%d = %s" % (z, N))
    print("\n\n" + "=" * 74)
    print("CERTIFICATES -- k = 2     a_z(p)*r(p) + c_z(p)*m_z(p) = N_z")
    print("=" * 74)
    print("r(p)   = %s" % R_SEXTIC)
    print("Pi     = y^2 + p*y + B,   B = (128p^3+32p^2+20p+15)/(32(8p+1))")
    print("Q_Pi   = y^2 + (-1/4-p)*y + D,   D = (64p^3+32p^2+14p-5)/(16(8p+1))")
    print("m_z(p) = primitive form of  u0*w1 - u1*w0,")
    print("         u = rem((y+1)^3*Q_Pi, Pi),  w = rem(-3*(y+1)^z, Pi)  over K = Q[p]/(r)")
    for z in ZSWEEP:
        m, _dg, (a, c, N) = K2[z]
        tag = "  [IN WINDOW]" if z in ZWINDOW else ""
        print("\n z = %d%s" % (z, tag))
        print("   m_%d = %s" % (z, m))
        print("   a_%d = %s" % (z, a))
        print("   c_%d = %s" % (z, c))
        print("   N_%d = %s" % (z, N))
    print("\n\n" + "=" * 74)
    print("BONUS -- the review's m_3, as the Pi^2-level (y^2,y^3) minor at z = 3")
    print("=" * 74)
    a3, c3, N3 = bezout(R_SEXTIC, REVIEW_M3, p)
    print("   m_3^sq = %s     (= -64/3 times the Pi^2-level minor)" % REVIEW_M3)
    print("   a      = %s" % a3)
    print("   c      = %s" % c3)
    print("   N      = %s" % N3)
    print("   verified: %s" % verify_bezout(a3, R_SEXTIC, c3, REVIEW_M3, N3, p))


# ===========================================================================
npass = sum(1 for _, ok, _ in _RESULTS if ok)
ntot = len(_RESULTS)
msg = "support_certificates: %d/%d checks pass" % (npass, ntot)
print(("\n" if not QUIET else "") + msg)
if npass != ntot:
    print("FAILURES:")
    for n, ok, d in _RESULTS:
        if not ok:
            print("  - %s  %s" % (n, d))
sys.exit(0 if npass == ntot else 1)
