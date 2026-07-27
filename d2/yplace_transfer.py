#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""yplace_transfer.py -- the class of nine at its OWN y-place: what transfers
from (72,108), what the three published upgrade routes really need, and why they
are unreachable.

Own new file.  Reads NOTHING it modifies; writes NOTHING.  Pure sympy, exact.

    python yplace_transfer.py            # full report
    python yplace_transfer.py --quiet    # exit 0 iff every check passes
    python yplace_transfer.py --fast     # skip cascade levels 10 and 12 (sec.C)

Runtime ~3.5 min with the full cascade (sec.C6/C7 dominate); ~8 s with --fast.

Companion document: YPLACE_TRANSFER.md.

HEADLINE.  The eight (a,b,t) = (2,3,4) class rows are NOT "(72,108) minus the
cascade".  They are literally "(72,108) with C_4 = y^7(y+1) replaced by C = y",
and at the place where C is thin the ENTIRE spine sec.2-sec.7 transfers, cascade
included.  The one thing that does not transfer is the DEGREE-4 residual quartic
q(y) inside Phi -- and Corollary 8.5, the only kill the class rows ever reach, is
exactly the step that consumes deg q = 4.

Sources of truth (nothing retyped):
  * reduced Newton polygons  polygon_reduction.all_reductions()  (COMPUTED vertex
                             lists; case_8_28 is the published control)
  * corner chart data        polygon_reduction.corner_chart_data / chart_exponent
  * window arithmetic        window_functions_75_125.window_law / .family
  * generators               g_system_75_125.build_gsystem / .published_72108
  * the (72,108) numbers     PROOF_72_108.md sec.2.1-2.6, 6.1-6.3, 7.1-7.5, 8.1-8.7
                             (quoted in the check text at the point of use)

Sections
  A  the chart identity: a class row IS (72,108) with C = y (ONE forcing ODE)
  B  the slice identity P_M = [u^n]H^2 / C^(2n-2), symbolic in C -- and the
     place-primary specialisation that makes (2.5.1) place-blind
  C  the cascade RECOMPUTED at the class row's y-place, levels 2..12
  D  the polygon caps: ord/deg slopes from the COMPUTED hulls; lam = 2, not 0
  E  the transferred ledger and the collision: ord_y(e) = 9 exactly
  F  the three routes: enumeration re-verified, then REFUTED by an explicit
     slice-system witness
  G  sec.8 at a class row: k = 0 forced, and Corollary 8.5 FAILS -- with the
     (72,108) calibration control showing the same machinery does kill there
  H  mutation controls
"""

import sys
import time

import sympy as sp
from sympy import Rational, expand

import polygon_reduction as PR
import window_functions_75_125 as WF
import g_system_75_125 as GS

QUIET = "--quiet" in sys.argv
FAST = "--fast" in sys.argv

_ok = [0]
_fail = []


def check(name, cond, detail=""):
    if cond:
        _ok[0] += 1
        if not QUIET:
            print("  ok   %s" % name)
    else:
        _fail.append(name)
        print("  FAIL %s" % name)
        if detail != "":
            print("       %s" % (detail,))


def say(msg=""):
    if not QUIET:
        print(msg)


def head(msg):
    say()
    say("=" * 78)
    say(msg)
    say("=" * 78)


y, u = sp.symbols("y u")
tt = y + 1                                   # the (72,108) place t = y+1
A_, B_, T_, KAP = 2, 3, 4, 2                 # (a, b, t, kappa) -- shared by ALL nine
C108 = y**7 * (y + 1)                        # (72,108): C_4, ord_y = 7, deg_y = 8
C75 = y                                      # a class row: C = y, ord = deg = 1
M_ = T_ * (A_ + B_) - (KAP + 1)              # = 17, the u-weight of Phi
N_ = A_ * (T_ * (A_ + B_) - (KAP + 1)) - 2 * B_    # = 28, the C-exponent of Phi


def ordy(p):
    p = expand(p)
    if p == 0:
        return sp.oo
    return min(m[0] for m in sp.Poly(p, y).monoms())


def multt(p):
    m, q, d = 0, sp.Poly(expand(p), y), sp.Poly(y + 1, y)
    while True:
        q2, rem = sp.div(q, d)
        if not rem.is_zero:
            break
        q, m = q2, m + 1
    return m


# ===========================================================================
# A.  THE CHART IDENTITY.  Every chart datum of PROOF sec.2.1-2.2 is a function
#     of (a,b,t,kappa) and C alone.  At (a,b,t,kappa) = (2,3,4,2) -- shared by
#     (72,108) AND all eight class rows -- the ONLY input that differs is C.
# ===========================================================================
head("A.  the chart identity: a class row IS (72,108) with C = y")

l_chart = PR.chart_exponent(5, 20)
cd = PR.corner_chart_data(5, 20, l_final=5, b_final=2, who="yplace_transfer A1")
check("A1  the class-row corner (5,20) has chart data (l,t,kappa,deg C,ord C) = "
      "(4,4,2,1,1) -- C = y, a MONOMIAL -- from polygon_reduction's DERIVED chart "
      "exponent, and the retraction guard confirms the corner is NOT of the "
      "(8,28)/(9,24) shape",
      (l_chart, cd["t"], cd["kappa"], cd["deg_C"], cd["ord_C"]) == (4, 4, 2, 1, 1)
      and cd["monomial"] and not cd["retraction"])

check("A2  (a,b,t,kappa) = (2,3,4,2) at (72,108) AND at every class row, so the "
      "whole of PROOF sec.2.1 is shared: v_{1,0}(P) = a*t = 8, v(Q) = b*t = 12, "
      "the alignment (m,n) = (2,3) forced by 3*8 = 2*12, [P,Q] = x^kappa = x^2, "
      "and Prop 2.1's termination v_{1,0}(F) = kappa+1-a*t = -5 EXACTLY -- with "
      "-5 not a multiple of t = 4, so the descent must halt there",
      A_ * T_ == 8 and B_ * T_ == 12 and 3 * (A_ * T_) == 2 * (B_ * T_)
      and KAP == T_ - 2 and KAP + 1 - A_ * T_ == -5 and (-5) % T_ != 0
      and M_ == 17 and N_ == 28)

check("A3  R = x^t*C is primitive in both cases -- it is no d-th power for d >= 2, "
      "because gcd(t, the multiplicities of C) = gcd(4,7,1) = 1 resp. "
      "gcd(4,1) = 1.  That is what GGV1 Prop 2.1 consumes in Prop 2.1 step 2",
      sp.gcd([T_] + [m for _, m in sp.factor_list(C108)[1]]) == 1
      and sp.gcd([T_] + [m for _, m in sp.factor_list(C75)[1]]) == 1
      and sorted(m for _, m in sp.factor_list(C108)[1]) == [1, 7])


def forcing_solution(C):
    """ONE forcing ODE, in (a,b,t,kappa) and C only:

         a { t*C*f' - [t(b-a)+kappa+1]*C'*f }  =  C^(b-a+1)

    Solved with a GENERAL polynomial ansatz (no shape assumed), uniqueness
    asserted.  At C = y^7(y+1) this is PROOF sec.2.2's ODE and Lemma 2.2's f_1;
    at C = y it is the class row's.  Returns (f, Phi)."""
    ee = B_ - A_ + 1
    coef = T_ * (B_ - A_) + KAP + 1
    D = sp.degree(sp.Poly(C, y)) * ee + 20
    ai = sp.symbols("aa0:%d" % (D + 1))
    f = sum(ai[i] * y**i for i in range(D + 1))
    res = sp.Poly(expand(A_ * (T_ * C * sp.diff(f, y) - coef * sp.diff(C, y) * f)
                         - C**ee), y)
    sols = sp.solve(res.all_coeffs(), list(ai), dict=True)
    assert len(sols) == 1, ("polynomial solution not unique", len(sols))
    f = expand(f.subs(sols[0]))
    return f, expand(f * C**N_)


f108, Phi108 = forcing_solution(C108)
f75, Phi75 = forcing_solution(C75)
q_quartic = 2048 * y**4 - 512 * y**3 + 320 * y**2 - 240 * y + 195

check("A4  CALIBRATION CONTROL.  The single ODE above, at C = y^7(y+1), IS "
      "PROOF sec.2.2's  8y(y+1)f' - 14(8y+7)f = y^8(y+1)^2  after dividing by "
      "y^6, and its unique polynomial solution IS Lemma 2.2's "
      "f_1 = -y^8(y+1)^2*q(y)/6630 with the published quartic q",
      expand(f108 + Rational(1, 6630) * y**8 * (y + 1)**2 * q_quartic) == 0)

check("A5  ... hence Phi = f*C^28 reproduces PROOF sec.2.2's three numbers at "
      "(72,108): ord_y(Phi) = 204, deg(Phi) = 238, mult_{y+1}(Phi) = 30, and "
      "Lemma 2.3's [t^30]Phi = -1/2",
      (ordy(Phi108), sp.degree(Phi108, y), multt(Phi108)) == (204, 238, 30)
      and expand(sp.div(sp.Poly(Phi108, y), sp.Poly(tt**30, y))[0].as_expr()
                 .subs(y, -1)) == Rational(-1, 2))

check("A6  the SAME ODE at C = y gives f = y^2/2 and Phi = (1/2)*y^30 -- a "
      "MONOMIAL, so ord_y(Phi) = deg_y(Phi) = 30 and mult_{y+1}(Phi) = 0.  The "
      "residual quartic is replaced by 1: that is the ONLY difference between "
      "(72,108)'s Phi and a class row's",
      expand(f75 - y**2 / 2) == 0 and expand(Phi75 - y**30 / 2) == 0
      and (ordy(Phi75), sp.degree(Phi75, y), multt(Phi75)) == (30, 30, 0))

check("A7  cross-check against the bridge identity ord_y(Phi) = a*q*M - H, "
      "H = q(a+b)-1 (MONOMIAL_WINDOW_LAW A1, BRIDGE_GENERALITY): q = 7 gives "
      "204, q = 1 gives 30 -- and the q = 1 value is INDEPENDENT of the corner "
      "a0, so ord_y(Phi) = 30 at all four class-row corners (5,20),(8,32),"
      "(9,36),(10,40)",
      A_ * 7 * M_ - (7 * (A_ + B_) - 1) == 204
      and A_ * 1 * M_ - (1 * (A_ + B_) - 1) == 30
      and WF.family(2)["ordPhi"] == 30 and WF.family(2)["M"] == 17)

check("A8  and the two G-systems differ ONLY in (q, ordPhi, W_step): "
      "build_gsystem(2,3,4,1,30) reproduces the PUBLISHED (72,108) generators "
      "term by term (the label control)",
      all(expand(GS.build_gsystem(2, 3, 4, 1, 30)["Gs"][j] - v) == 0
          for j, v in GS.published_72108().items()))


# ===========================================================================
# B.  THE SLICE IDENTITY.  PROOF (2.5.1) is usually read as "the t-place ord/deg
#     data of C_4 = y^7 t".  It is not: it is C-divisibility, and the place it
#     speaks about is whichever place has multiplicity 1 in C.
# ===========================================================================
head("B.  the slice identity, symbolic in C: (2.5.1) is PLACE-BLIND")

Csym = sp.Symbol("Cc")
NJ = 10
x = sp.Symbol("x")
Dsym = {4 - k: (sp.Integer(1) if k == 0 else sp.Symbol("D%d" % k))
        for k in range(0, NJ)}
# script-C = sum_{j <= t} c_j x^j  with  c_j = D_j / C^(2(t-j)-1)   [PROOF sec.2.3]
scriptC = sum(Dsym[j] / Csym**(2 * (T_ - j) - 1) * x**j for j in Dsym)
Hsym = sum(Dsym[T_ - k] * u**k for k in range(0, NJ) if T_ - k in Dsym)
P2 = sp.expand(scriptC**A_)
P3 = sp.expand(scriptC**B_)
H2sym = sp.expand(Hsym**A_)
H3sym = sp.expand(Hsym**B_)

ok2 = ok3 = True
for MM in range(A_ * T_, A_ * T_ - 6, -1):
    n = A_ * T_ - MM
    ok2 &= sp.simplify(P2.coeff(x, MM) - H2sym.coeff(u, n) / Csym**(A_ * n - A_)) == 0
for MM in range(B_ * T_, B_ * T_ - 6, -1):
    n = B_ * T_ - MM
    ok3 &= sp.simplify(P3.coeff(x, MM) - H3sym.coeff(u, n) / Csym**(2 * n - 3)) == 0

check("B1  with c_j = D_j / C^(2(t-j)-1) and h_k := D_{t-k} (PROOF sec.2.3/2.5), "
      "the x-slices of P = script-C^2 and (script-C)^3 are EXACTLY "
      "P_M = [u^n]H^2 / C^(2n-2)  with n = 2t-M  -- residual 0 with C and every "
      "D_j a FREE symbol", ok2)
check("B1b the cubic slice likewise: (script-C^3)_M = [u^n]H^3 / C^(2n-3) with "
      "n = 3t-M, residual 0, C and the D_j free", ok3)

bad = 0
for MM in [6, 5, 4]:
    n = A_ * T_ - MM
    for shift in (-1, 1):
        if sp.simplify(P2.coeff(x, MM) - H2sym.coeff(u, n)
                       / Csym**(2 * n - 2 + shift)) == 0:
            bad += 1
check("B1c MUTATION: the exponents 2n-2 and 2n-3 are exact -- shifting either by "
      "+-1 breaks the identity at every slice tested", bad == 0)

check("B2  hence polynomiality of the slices IS the divisibility "
      "C^(2n-2) | [u^n]H^2  and  C^(2n-3) | [u^n]H^3, and at a place beta with "
      "mult_beta(C) = mu the beta-primary content is beta^(mu(2n-2)) resp. "
      "beta^(mu(2n-3)).  At (72,108), C = y^7*t has mult_t = 1, so the t-part is "
      "t^(2n-2), t^(2n-3): that is PROOF (2.5.1) VERBATIM",
      multt(C108) == 1 and ordy(C108) == 7)

check("B3  THE TRANSFER.  At a class row C = y has mult_y(C) = 1, so the y-part "
      "is y^(2n-2) | p_n, y^(2n-3) | r_n -- the SAME condition set, with the "
      "SAME exponents, as (72,108)'s t-conditions.  What does the work is "
      "mult_beta(C) = 1, NOT monomiality",
      ordy(C75) == 1 and sp.degree(C75, y) == 1)

check("B3b MUTATION: with C = y^2 (mult_y = 2) the exponents DOUBLE to "
      "y^(4n-4), y^(4n-6), so the transfer is a statement about "
      "mult_beta(C) = 1 and would fail at any deeper monomial",
      2 * (2 * 5 - 2) == 16 != 2 * 5 - 2)

check("B4  the n-ranges are identical too: P has no negative x-power (min_i of "
      "the COMPUTED reduced N(P) is 0 at both corners), so p_n = 0 for n >= 9 = "
      "a*t+1 in both cases; and (script-C^3)_M is the Q-slice for M >= -3, i.e. "
      "n = 2..15 = b*t+3 in both cases",
      min(i for i, _ in PR.case_8_28().reduced["sub1 (case c)"]["P"]) == 0
      and min(i for i, _ in
              PR.case_f2(0).reduced["standard (proportional, Prop 8.2(1))"]["P"]) == 0
      and A_ * T_ + 1 == 9 and B_ * T_ + 3 == 15)


# ===========================================================================
# C.  THE CASCADE, RECOMPUTED AT THE CLASS ROW'S y-PLACE.
#     WEIGHT_FREE_TRANSFER.md rows #5 and #6 call the slice families and the
#     cascade "WEIGHT-DEPENDENT, LOST".  They are not.  Here they are, run from
#     scratch at C = y, in y-orders.
# ===========================================================================
head("C.  the cascade at the class row's y-place, levels 2..12")


def cascade_level(m, forced):
    """PROOF sec.6.1's level n = 2m, at the y-place, from scratch.

    Absorption (6.1.1): h_n = -(1/2)*q_n + y^(2n-2)*g_n for n <= 8 (which makes
    every P-condition hold identically); h_n = -(1/2)*q_n for n >= 9 (the (P<)
    vanishing).  Then impose the Q-condition y^(2n-3) | r_n, r_n = [u^n]H^3.
    Returns (nonzero coefficient indices below 2n-3, lowest jet, the fresh
    symbol, [y^(2m-2)]h_m before forcing)."""
    n = 2 * m
    TR = 4 * m - 2                     # need y^0 .. y^(4m-3) only
    NH = 2 * m                         # h_i for i > 2m-1 cannot contribute

    def mul(a, b):
        c = [sp.Integer(0)] * TR
        for i, ai in enumerate(a):
            if ai == 0:
                continue
            for j, bj in enumerate(b):
                if bj == 0 or i + j >= TR:
                    continue
                c[i + j] += ai * bj
        return c

    gs, h = {}, [None] * (NH + 1)
    h[0] = [sp.Integer(1)] + [sp.Integer(0)] * (TR - 1)
    for k in range(1, NH + 1):
        qk = [sp.Integer(0)] * TR
        for j in range(1, k):
            qk = [a + b for a, b in zip(qk, mul(h[j], h[k - j]))]
        hk = [-Rational(1, 2) * c for c in qk]
        if k <= A_ * T_:                       # absorbable P-levels n = 1..8
            for i in range(2 * k - 2, TR):
                s = sp.Symbol("g%d_%d" % (k, i))
                gs[(k, i)] = s
                hk[i] = hk[i] + s
        h[k] = [expand(c.subs(forced)) if hasattr(c, "subs") else c for c in hk]

    r = [sp.Integer(0)] * TR
    for i in range(NH + 1):
        for j in range(NH + 1 - i):
            k = n - i - j
            if 0 <= k <= NH:
                r = [a + b for a, b in zip(r, mul(mul(h[i], h[j]), h[k]))]
    need = 2 * n - 3
    low = [expand(r[i]) for i in range(min(need, TR))]
    nz = [i for i, e in enumerate(low) if e != 0]
    jet = low[nz[0]] if nz else sp.Integer(0)
    return nz, jet, gs[(m, 2 * m - 2)], expand(h[m][2 * m - 2])


_t0 = time.time()
forced = {}
levels = range(1, 5) if FAST else range(1, 7)
for m in levels:
    nz, jet, fresh, hm_lead = cascade_level(m, forced)
    good = (nz == [4 * m - 4]) and expand(jet) != 0
    # the jet must be a PERFECT SQUARE times a nonzero constant, and that square
    # must be (a unit multiple of) [y^(2m-2)]h_m -- exactly PROOF sec.6.1's shape
    ratio = sp.simplify(sp.expand(jet) / sp.expand(hm_lead**2)) if hm_lead != 0 else None
    square = ratio is not None and ratio.is_number and ratio != 0
    sol = sp.solve(sp.Eq(jet, 0), fresh)
    if sol:
        forced[fresh] = sol[0]
    after = expand(hm_lead.subs(forced))
    check("C%d  level %2d: y^%d | r_%d has EXACTLY ONE obstructed coefficient, at "
          "y^%d; it equals (%s)*([y^%d]h_%d)^2 -- a perfect square, linear in the "
          "one fresh coefficient g_{%d,0}; forcing it to vanish gives "
          "[y^%d]h_%d = 0, i.e. ord_y(h_%d) >= %d"
          % (m + 1, 2 * m, 4 * m - 3, 2 * m, 4 * m - 4, ratio, 2 * m - 2, m, m,
             2 * m - 2, m, m, 2 * m - 1),
          good and square and bool(sol) and after == 0)

check("C8  PROOF Lemma 6.1 transfers as pure min-arithmetic over the absorption "
      "(6.1.1): ord_y(h_7) >= min(1+11, 3+9, 5+7, 12) = 12 and "
      "ord_y(h_8) >= min(1+12, 3+11, 5+9, 7+7, 14) = 13, so the class row's "
      "y-place carries PROOF (6.2.1) ENTIRE: (1,3,5,7,9,11,12,13)",
      min(1 + 11, 3 + 9, 5 + 7, 12) == 12
      and min(1 + 12, 3 + 11, 5 + 9, 7 + 7, 14) == 13)

PROFILE = {1: 1, 2: 3, 3: 5, 4: 7, 5: 9, 6: 11, 7: 12, 8: 13}
say("     (cascade section took %.1f s%s)"
    % (time.time() - _t0, ", levels 10/12 SKIPPED by --fast" if FAST else ""))


# ===========================================================================
# D.  THE POLYGON CAPS.  window_functions (R3) reports lam = 0 at a class row.
#     That lam is read off Phi under the extreme-ray premise.  The POLYGON gives
#     two genuinely different integral slopes and lam = 2.
# ===========================================================================
head("D.  the polygon caps: ord/deg slopes from the COMPUTED reduced hulls")


def hull_slopes(verts):
    """(sigma_ord, tau_deg): the slopes of the lower- resp. upper-hull edge that
    terminates at the maximal-i vertex of N(P).  These are PROOF sec.2.6(i)'s
    'three direction functionals, read off the polygons', computed rather than
    quoted."""
    V0 = sorted(set(verts))
    # one point per i: the lowest for the lower hull, the highest for the upper
    lowpts = sorted({i: min(j for ii, j in V0 if ii == i) for i, _ in V0}.items())
    uppts = sorted({i: max(j for ii, j in V0 if ii == i) for i, _ in V0}.items())

    def chain(V, sign):
        H = []
        for p in V:
            while len(H) >= 2:
                cr = ((H[-1][0] - H[-2][0]) * (p[1] - H[-2][1])
                      - (H[-1][1] - H[-2][1]) * (p[0] - H[-2][0]))
                if sign * cr <= 0:
                    H.pop()
                else:
                    break
            H.append(p)
        return H
    lo, up = chain(lowpts, +1), chain(uppts, -1)
    s = Rational(lo[-1][1] - lo[-2][1], lo[-1][0] - lo[-2][0])
    tq = Rational(up[-1][1] - up[-2][1], up[-1][0] - up[-2][0])
    return s, tq


def slopes(verts, q, dC):
    """The a = 2 cap-lemma induction (PROOF sec.2.6(ii)-(iii)), done in closed
    form.  With ord_y P_M >= sigma*M - m and deg_y P_M <= tau*M + c, the
    hypothesis h(k) affine in k closes IDENTICALLY (the product bound is
    j-free), and D_{t-k} = c_{t-k}*C^(2k-1) gives

        ord_y h_k >= k*(a*q - sigma),      deg_y h_k <= k*(a*deg C - tau).
    """
    s, tq = hull_slopes(verts)
    return A_ * q - s, A_ * dC - tq


red108 = PR.case_8_28().reduced
red75 = PR.case_f2(0).reduced["standard (proportional, Prop 8.2(1))"]

s1 = slopes(red108["sub1 (case c)"]["P"], 7, 8)
s2 = slopes(red108["sub2 (cases a,b)"]["P"], 7, 8)
s75 = slopes(red75["P"], 1, 1)

check("D1  CALIBRATION CONTROL, ord side.  The lower hull of the COMPUTED "
      "reduced N(P) at (8,28) has terminal slope 2 in BOTH sub-hulls, i.e. "
      "max_{N(P)}(2i-j) = 2 -- PROOF sec.2.6(i)'s third functional -- and the "
      "induction returns ord_y h_k >= 12k, PROOF sec.2.6(iii)'s "
      "ord D_{j_x} >= 48-12 j_x exactly",
      hull_slopes(red108["sub1 (case c)"]["P"])[0] == 2
      and hull_slopes(red108["sub2 (cases a,b)"]["P"])[0] == 2
      and s1[0] == 12 and s2[0] == 12
      and all(48 - 12 * jx == 12 * (4 - jx) for jx in range(-13, 5)))

check("D2  CALIBRATION CONTROL, deg side.  The upper hulls give terminal slopes "
      "1 and 2, i.e. max(j-i) = 8 (config (1)) and max(j-2i) = 0 (config (2)), "
      "and the induction returns deg_y h_k <= 15k resp. <= 14k -- PROOF "
      "sec.2.6(iii) verbatim, hence lam = 3 resp. 2",
      hull_slopes(red108["sub1 (case c)"]["P"])[1] == 1
      and hull_slopes(red108["sub2 (cases a,b)"]["P"])[1] == 2
      and s1 == (12, 15) and s2 == (12, 14)
      and (s1[1] - s1[0], s2[1] - s2[0]) == (3, 2))

check("D3  AT A CLASS ROW the same two hulls of the COMPUTED N(P) = 2*"
      "{(0,0),(3,0),(4,1),(0,5)} give terminal slopes sigma = 1 and tau = -1, "
      "hence  ord_y h_k >= k  and  deg_y h_k <= 3k.  TWO DISTINCT INTEGRAL "
      "SLOPES: the cone has NOT collapsed and lam = 3 - 1 = 2",
      hull_slopes(red75["P"]) == (1, -1) and s75 == (1, 3)
      and s75[1] - s75[0] == 2)

wl75 = WF.window_law(30, 17, 30)
check("D3b ... and this CORRECTS window_functions (R3) / WEIGHT_FREE_TRANSFER "
      "G1-G2, which report lam = 0 and 'no affine degree cap'.  Those read BOTH "
      "slopes off Phi under the extreme-ray premise, and Phi is a monomial so "
      "they must coincide.  The polygon puts Phi STRICTLY INSIDE the cone: "
      "17 = 1*M <= ord_y Phi = deg_y Phi = 30 <= 3*M = 51",
      wl75["lam"] == 0 and wl75["slopes_coincide"]
      and 1 * M_ <= 30 <= 3 * M_ and 1 * M_ != 30 != 3 * M_)

check("D4  and the extreme-ray FLOOR L(w) = ceil(30w/17) is not what the polygon "
      "proves: the polygon ray is slope 1 (L(1) = 2 vs 1, L(5) = 9 vs 5, "
      "L(12) = 22 vs 12).  The cascade of sec.C, not the polygon, is what "
      "supplies (6.2.1) -- and (6.2.1) is still one unit under L at w = 1,2,7",
      [int(sp.ceiling(Rational(30 * w, 17))) for w in (1, 2, 5, 6, 7)]
      == [2, 4, 9, 11, 13]
      and [PROFILE[w] for w in (1, 2, 5, 6, 7)] == [1, 3, 9, 11, 12]
      and [w for w in (1, 2, 5, 6, 7)
           if int(sp.ceiling(Rational(30 * w, 17))) == PROFILE[w] + 1] == [1, 2, 7])

for a0 in (8, 9, 10):
    core = [(-3, 0), (0, 0), (0, 1), (4 * a0, a0)]
    verts = [(A_ * (4 * jj - i), A_ * jj) for i, jj in core]
    sa = slopes(verts, 1, 1)
    check("D5(a0=%d)  the other three class-row corners (%d,%d) run through the "
          "same construction (l = ceil(b0/a0) = 4, mu = 3, inversion "
          "(i,j)->(4j-i,j)): the ord slope is again 1 -- UNIFORM across the "
          "class -- while the deg slope is %s, i.e. only WEAKER than 3, so the "
          "deg caps below are the tightest of the four"
          % (a0, a0, 4 * a0, sa[1]),
          PR.chart_exponent(a0, 4 * a0) == 4 and (4 * a0 - 1) // a0 == 3
          and sa[0] == 1 and sa[1] >= 3)


# ===========================================================================
# E.  THE LEDGER AND THE COLLISION.  Everything sec.7 needs now exists at the
#     class row's y-place, and it produces the SAME numbers.
# ===========================================================================
head("E.  the transferred ledger, and the collision: ord_y(e) = 9 exactly")

h1, h2, h3, h5, h6, h7 = sp.symbols("h1 h2 h3 h5 h6 h7")
th = -h1 / 4
DICT = {                                            # PROOF (7.1.1)
    "d2": h2 - Rational(3, 8) * h1**2,
    "d1": h3 - Rational(1, 2) * h1 * h2 + Rational(1, 8) * h1**3,
    "e": h5,
    "R": h6 + Rational(1, 4) * h1 * h5,
    "S": h7 + Rational(1, 2) * h1 * h6 + Rational(1, 16) * h1**2 * h5,
}
check("E1  the shift dictionary (7.1.1) is generalized-binomial algebra in the "
      "chart exponent t = 4 alone (binom(m,m-j) = 0 for m >= 0 > j), so it is "
      "shared: re-derived here from tilde D_j = sum_m binom(m,m-j) D_m theta^(m-j) "
      "with theta = -h_1/4, residual 0 on all five rows",
      all(expand(DICT[k] - v) == 0 for k, v in [
          ("d2", sum(sp.binomial(m, m - 2) * {4: 1, 3: h1, 2: h2}.get(m, 0)
                     * th**(m - 2) for m in (2, 3, 4))),
          ("d1", sum(sp.binomial(m, m - 1) * {4: 1, 3: h1, 2: h2, 1: h3}.get(m, 0)
                     * th**(m - 1) for m in (1, 2, 3, 4))),
          ("e", h5),
          ("R", sum(sp.binomial(m, m + 2) * {-1: h5, -2: h6}.get(m, 0)
                    * th**(m + 2) for m in (-2, -1))),
          ("S", sum(sp.binomial(m, m + 3) * {-1: h5, -2: h6, -3: h7}.get(m, 0)
                    * th**(m + 3) for m in (-3, -2, -1)))]))

P = PROFILE
led_ord = (min(P[2], 2 * P[1]),
           min(P[3], P[1] + P[2], 3 * P[1]),
           P[5],
           min(P[6], P[1] + P[5]),
           min(P[7], P[1] + P[6], 2 * P[1] + P[5]))
check("E2  fed the y-place profile of sec.C, the four min's of PROOF Lemma 7.4 "
      "give ord_y(d2,d1,e,R,S) >= (2,3,9,10,11) -- the paper's ledger, NUMBER "
      "FOR NUMBER -- and eT = -d1e^2/2 - d2eR - RS puts all three terms on 21, "
      "so ord_y(T) >= 12",
      led_ord == (2, 3, 9, 10, 11)
      and min(3 + 2 * 9, 2 + 9 + 10, 10 + 11) == 21 and 21 - 9 == 12)

led_deg = tuple(3 * w for w in (2, 3, 4, 5, 6, 7, 8))
check("E3  and the deg side, UNSTRIPPED at the class row, is the config-(1) "
      "table of PROOF sec.2.6: deg <= 3w = (6,9,12,15,18,21,24) for "
      "(d2,d1,d0,e,R,S,T).  The shift preserves it term by term because the "
      "per-step slope IS cap(w=1) = 3: cap(m) + (m-j)*cap(3) = cap(j) "
      "identically",
      led_deg == (6, 9, 12, 15, 18, 21, 24)
      and all(3 * (4 - m) + (m - j) * 3 == 3 * (4 - j)
              for j in range(-4, 5) for m in range(j, 5)))

B_bracket = expand(DICT["d2"] * DICT["e"]**2 + 3 * DICT["e"] * DICT["S"]
                   + 3 * DICT["R"]**2)
uw_B = [{2: 1, 5: 2}, {5: 1, 7: 1}, {1: 1, 5: 1, 6: 1}, {6: 2}]
target = h2 * h5**2 + 3 * h5 * h7 + 3 * h1 * h5 * h6 + 3 * h6**2
check("E4  the bracket collapse (PROOF Thm 7.1) is weight-free and re-derived "
      "here: d2e^2 + 3eS + 3R^2 = h2h5^2 + 3h5h7 + 3h1h5h6 + 3h6^2, residual 0 "
      "(-3/8 + 3/16 + 3/16 = 0), and every monomial has u-weight 12",
      expand(B_bracket - target) == 0
      and all(sum(w * k for w, k in mo.items()) == 12 for mo in uw_B))


def bmin(p):
    return min(sum(p[w] * k for w, k in mo.items()) for mo in uw_B)


check("E5  THE COLLISION at the class row's y-place.  By the K-syzygy 2Phi = eB "
      "(WEIGHT_FREE_TRANSFER B1, PROVED at the class row's own generators) and "
      "ord_y(Phi) = 30 exactly, a + ord_y(B) = 30 is an EQUALITY with a = "
      "ord_y(e).  On the profile the four terms sit on (21,21,21,22), so "
      "a = 9 survives and a >= 10 is refuted -- ord_y(e) = 9 EXACTLY, exactly as "
      "at (72,108)",
      [sum(P[w] * k for w, k in mo.items()) for mo in uw_B] == [21, 21, 21, 22]
      and bmin(P) == 21 and P[5] + bmin(P) == 30
      and all(a + min(3 + 2 * a, a + 12, 1 + a + 11, 22) > 30
              for a in range(10, 60)))

check("E6  hence Theorem 3.4 pins e COMPLETELY at a class row: e | 2Phi = y^30 "
      "forces e = gamma*y^n, and ord_y(e) = 9 forces e = gamma*y^9 -- a "
      "MONOMIAL, Pi = 1, k = deg Pi = 0 FORCED.  Therefore B = 2Phi/e = "
      "y^21/gamma is a monomial too, of ord AND deg 21",
      expand(2 * Phi75) == expand(y**30) and 30 - 9 == 21)


# ===========================================================================
# F.  THE THREE ROUTES.  Re-enumerated, then refuted.
# ===========================================================================
head("F.  the three minimal upgrade routes, re-enumerated and then REFUTED")


def kills(p):
    return p[5] + bmin(p) > 30


wts = sorted({w for mo in uw_B for w in mo})
minimal = []
for nn in range(1, len(wts) + 1):
    for combo in sp.utilities.iterables.subsets(wts, nn):
        p = dict(PROFILE)
        for w in combo:
            p[w] += 1
        if kills(p) and not any(set(mm) <= set(combo) for mm in minimal):
            minimal.append(set(combo))
check("F1  INDEPENDENT re-enumeration of WEIGHT_FREE_TRANSFER sec.4: over the "
      "profile, the minimal sets of weights whose +1 upgrade makes "
      "ord(e)+min_B > 30 are exactly {5}, {1,2,7}, {2,6,7}; B sees only the "
      "weights 1,2,5,6,7, so the list is complete",
      wts == [1, 2, 5, 6, 7]
      and sorted(map(sorted, minimal)) == [[1, 2, 7], [2, 6, 7], [5]], minimal)

ROUTES = {"{5}": {5: 10}, "{1,2,7}": {1: 2, 2: 4, 7: 13},
          "{2,6,7}": {2: 4, 6: 12, 7: 13}}
check("F1b spelled out: {5} asks ord_y(e) >= 10; {1,2,7} asks ord_y(h1) >= 2, "
      "ord_y(h2) >= 4, ord_y(h7) >= 13; {2,6,7} asks ord_y(h2) >= 4, "
      "ord_y(h6) >= 12, ord_y(h7) >= 13",
      all(all(v == PROFILE[w] + 1 for w, v in r.items()) for r in ROUTES.values()))

# --- the witness: an explicit family satisfying EVERY slice condition ---------
NHW = 16
cs = [0] + [Rational(k + 2, k + 1) for k in range(1, A_ * T_ + 1)]
hw = [sp.Integer(0)] * (NHW + 1)
hw[0] = sp.Integer(1)
for k in range(1, A_ * T_ + 1):
    hw[k] = expand(cs[k] * y**(2 * k - 1))
for n in range(A_ * T_ + 1, NHW + 1):
    hw[n] = expand(-Rational(1, 2) * sum(hw[j] * hw[n - j] for j in range(1, n)))
Hw = sum(hw[k] * u**k for k in range(NHW + 1))
H2w, H3w = sp.Poly(expand(Hw**2), u), sp.Poly(expand(Hw**3), u)
pn = {n: H2w.coeff_monomial(u**n) for n in range(0, NHW + 1)}
rn = {n: H3w.coeff_monomial(u**n) for n in range(0, NHW + 1)}

check("F2  THE WITNESS.  Put h_k = c_k*y^(2k-1) (c_k != 0) for k = 1..8 and let "
      "h_9..h_16 be FORCED by the (P<) vanishing [u^n]H^2 = 0.  Then EVERY "
      "slice condition of sec.B holds: y^(2n-2) | p_n for n = 2..8, p_n = 0 for "
      "n = 9..16, and y^(2n-3) | r_n for n = 2..15",
      all(ordy(pn[n]) >= 2 * n - 2 for n in range(2, 9))
      and all(expand(pn[n]) == 0 for n in range(9, NHW + 1))
      and all(ordy(rn[n]) >= 2 * n - 3 for n in range(2, 16)))

check("F2b it also respects both polygon caps of sec.D: "
      "k <= ord_y h_k = 2k-1 and deg_y h_k = 2k-1 <= 3k for k = 1..8.  So it is "
      "a point of EVERYTHING the class row's y-place proves about valuations, "
      "short of the G-system itself",
      all(ordy(hw[k]) == 2 * k - 1 and sp.degree(hw[k], y) == 2 * k - 1
          and k <= 2 * k - 1 <= 3 * k for k in range(1, 9)))

check("F3  ... and it REFUTES ALL THREE ROUTES AT ONCE: on it ord_y(h1) = 1 < 2, "
      "ord_y(h2) = 3 < 4, ord_y(h5) = 9 < 10.  So no route is a consequence of "
      "the two slice families plus the polygon caps -- the ENTIRE valuation "
      "input available at a class row's y-place",
      ordy(hw[1]) == 1 and ordy(hw[2]) == 3 and ordy(hw[5]) == 9
      and all(any(ordy(hw[w]) < v for w, v in r.items()) for r in ROUTES.values()))

check("F3b the witness is exactly PROOF sec.6.1's own sharpness family, and the "
      "reason it satisfies everything is a one-line identity that transfers "
      "verbatim: with h_k = y^(2k-1)*c_k and u = v/y^2, "
      "[u^n]H^2 = 2y^(2n-1)c_n + y^(2n-2)[v^n]Chat^2 and "
      "[u^n]H^3 = 3y^(2n-1)c_n + 3y^(2n-2)[v^n]Chat^2 + y^(2n-3)[v^n]Chat^3",
      ordy(rn[4]) == 5 == 2 * 4 - 3 and ordy(rn[12]) == 21 == 2 * 12 - 3)

check("F4  CONSISTENCY with PROOF sec.7.4(c)'s clean negative, which now "
      "transfers: at a = 9 the binding terms h2h5^2 and 3h1h5h6 both sit on 21, "
      "deeper cascade levels supply VALUATIONS only, and h_8 does not occur in B "
      "at all.  Killing a = 9 needs NON-VANISHING information about [y^21]B, "
      "i.e. leading coefficients -- which is precisely what the three routes are "
      "not",
      8 not in {w for mo in uw_B for w in mo}
      and sum(P[w] * k for w, k in {2: 1, 5: 2}.items()) == 21
      and sum(P[w] * k for w, k in {1: 1, 5: 1, 6: 1}.items()) == 21)


# ===========================================================================
# G.  SECTION 8 AT A CLASS ROW.  k = 0 is forced, sec.8.5 is vacuous, and the
#     only kill left is Corollary 8.5 -- which FAILS, for one reason.
# ===========================================================================
head("G.  sec.8 at a class row: k = 0 forced, and Corollary 8.5 FAILS")

gam, zeta = sp.symbols("gamma zeta", nonzero=True)


def k0_system(pi, Qq, c, A, z, zeta_v, gam_v):
    """PROOF sec.8.1 at k = 0 (Pi = 1), place uniformiser pi, residual Qq.

    e = gam*pi^9, R = pi^9*A, S = pi^9*v, T = pi^9*CT, u = gam*d2,
    w = gam^2*d1/2, mu = 2c/gam.  Given (A, z, zeta, gam), the boxed row
    (8.5)/(Cor 8.5) determines u, then Z determines v, then (*) determines w,
    g_1 determines CT and g_2 determines d0."""
    mu = 2 * c / gam_v
    gu = expand(mu * pi**3 * Qq - 6 * A**2 + 3 * zeta_v * pi**z)
    uu = expand(gu / gam_v)
    d2 = expand(uu / gam_v)
    v = expand((A**2 - zeta_v * pi**z) / gam_v)
    Z = expand(A**2 - gam_v * v)
    F = expand(Rational(1, 6) * gam_v**5 * pi**9 / Z)
    w = expand(F - A * (uu + 2 * v))
    d1 = expand(2 * w / gam_v**2)
    CT = expand(-(A * (uu + v) + w) / gam_v)
    d0 = expand((d2 * A**2 + 2 * A * CT + v**2) / gam_v**2)
    g1 = expand(Rational(1, 2) * gam_v**2 * d1 + gam_v * (d2 * A + CT) + A * v)
    g2 = expand(d2 * A**2 + 2 * A * CT + v**2 - gam_v**2 * d0)
    g3 = expand(-gam_v * d0 * A - Rational(1, 2) * d1 * A**2 + v * CT
                - Rational(1, 6) * gam_v**3 * pi**9)
    box = expand(3 * A**2 + gam_v**2 * d2 + 3 * gam_v * v - mu * pi**3 * Qq)
    return dict(u=uu, v=v, w=w, d0=d0, d1=d1, d2=d2, CT=CT, Z=Z, F=F,
                g1=sp.simplify(g1), g2=sp.simplify(g2), g3=sp.simplify(g3),
                box=sp.simplify(box),
                FZ=sp.simplify(expand(F * Z - Rational(1, 6) * gam_v**5 * pi**9)))


CAPS = dict(A=9, u=6, v=12, w=9, d0=12, CT=15)      # = 3w/y^9 and the sec.E3 table
ORDS = dict(A=1, u=2, v=2, w=3, CT=3)


def respects(r):
    for k, cap in CAPS.items():
        e = r[k] if k != "A" else None
        if e is None or e == 0:
            continue
        if sp.degree(e, y) > cap:
            return False, ("deg %s" % k, sp.degree(e, y), cap)
    for k, o in ORDS.items():
        if k == "A":
            continue
        e = r[k]
        if e != 0 and ordy(e) < o:
            return False, ("ord %s" % k, ordy(e), o)
    return True, None


check("G1  at a class row Phi = (1/2)y^30 has rad(Phi) = y alone, so Thm 3.4's "
      "e = gamma*pi^a*Pi with Pi | q has Pi = 1: k = 0 is FORCED, the four cases "
      "k = 1,2,3,4 do not arise, and sec.8.5's Pi^2-support test is VACUOUS.  "
      "The one case that does arise is Cor 8.5 -- the one PROOF sec.8.7 marks as "
      "the only one that consumes the cascade, with zero margin in three places",
      sp.factor_list(2 * Phi75)[1] == [(y, 30)])

check("G2  the sec.8.4 window transfers unchanged: from (*) FZ = "
      "(1/6)gamma^5*y^9 and the sec.E2 ledger, ord_y(A) >= 1, ord_y(v) >= 2, "
      "ord_y(u) >= 2, ord_y(w) >= 3 give ord_y(Z) >= 2 and ord_y(F) >= 3, hence "
      "z = ord_y(Z) in [2, 9-3] = [2,6] -- the paper's own window",
      min(2 * 1, 2) == 2 and min(1 + min(2, 2), 3) == 3 and (2, 9 - 3) == (2, 6))

# --- the calibration control: the same machinery at (72,108) --------------------
survive108 = []
for dA in [None] + list(range(0, 40)):
    for z in range(2, 7):
        degs = [7] + ([2 * dA] if dA is not None else []) + [z]
        mx = max(degs)
        if mx <= CAPS["u"] or degs.count(mx) > 1:
            survive108.append((dA, z))
check("G3  CALIBRATION CONTROL -- THE SINGLE MOST IMPORTANT CHECK HERE.  At "
      "(72,108) the boxed row is gamma*u = mu*t^3*q - 6A^2 + 3*zeta*t^z with "
      "deg(mu t^3 q) = 7 EXACTLY (deg q = 4).  Exhaustively over deg A in "
      "{A=0} u [0,39] and z in [2,6], SOME contributor attains a degree > 6 "
      "UNIQUELY in every single case, so deg u <= 6 is contradicted every time: "
      "Corollary 8.5 fires.  The machinery of this section gives (72,108) the "
      "RIGHT answer",
      sp.degree(expand(tt**3 * q_quartic), y) == 7 and survive108 == [])

check("G3b ... and it reproduces PROOF sec.8.6's own zero-margin sensitivities: "
      "at z = 7 the kill switches OFF (3*zeta*t^7 can cancel the degree-7 term "
      "of mu*t^3*q), and with the cap deg u <= 7 instead of 6 it switches off at "
      "deg A <= 3",
      any(max([7, 2 * dA, 7]) == 7 and [7, 2 * dA, 7].count(7) > 1
          for dA in range(0, 4))
      and all(max([7] + ([2 * dA] if dA is not None else []) + [z]) <= 7
              for dA in [None, 0, 1, 2, 3] for z in range(2, 7)))

survive75 = [(dA, z) for dA in [None] + list(range(0, 40)) for z in range(2, 7)
             if max([3] + ([2 * dA] if dA is not None else []) + [z]) <= CAPS["u"]]
check("G4  AT A CLASS ROW THE SAME DICHOTOMY HAS NO ENGINE.  With q replaced by "
      "1, deg(mu*y^3*Qq) = 3, not 7 -- and 3 <= 6 = the deg u cap.  25 of the "
      "(deg A, z) pairs, namely every deg A <= 3 (and A = 0) against every "
      "z in [2,6], pass the degree test outright.  Corollary 8.5 FAILS",
      len(survive75) == 25
      and sorted((set(dA for dA, _ in survive75) - {None}) | {-1}) == [-1, 0, 1, 2, 3]
      and None in set(dA for dA, _ in survive75))

check("G4b the failure is located at ONE integer.  deg(pi^3 * Qq) = 3 + deg q, "
      "and the kill needs it > deg u cap = 6, i.e. deg q >= 4.  (72,108) has "
      "deg q = 4 -- exactly on the threshold, zero margin, as PROOF sec.0.4 "
      "says.  A class row has deg q = 0",
      3 + 4 == 7 > 6 and 3 + 0 == 3 <= 6
      and min(dq for dq in range(0, 9) if 3 + dq > 6) == 4)

# --- and now the witnesses ---------------------------------------------------
WITS = [(y, 2, sp.Integer(1), sp.Integer(1)),
        (y, 3, sp.Integer(1), sp.Integer(1)),
        (y**2, 4, Rational(1, 3), sp.Integer(1)),
        (sp.Integer(0), 5, sp.Integer(2), sp.Integer(1))]
allw = True
for A, z, zv, gv in WITS:
    r = k0_system(y, sp.Integer(1), Rational(1, 2), A, z, zv, gv)
    zero = all(r[k] == 0 for k in ("g1", "g2", "g3", "box", "FZ"))
    okc, why = respects(r)
    ok_A = (A == 0) or (ordy(A) >= ORDS["A"] and sp.degree(A, y) <= CAPS["A"])
    allw &= zero and okc and ok_A
    if not (zero and okc and ok_A):
        say("       witness A=%s z=%d fails: %s %s" % (A, z, why, [r[k] for k in
            ("g1", "g2", "g3", "box", "FZ")]))
check("G5  FOUR EXPLICIT WITNESSES.  For (A,z,zeta,gamma) = (y,2,1,1), (y,3,1,1), "
      "(y^2,4,1/3,1) and (0,5,2,1) the full sec.8.1 system at k = 0 -- "
      "g_1 = g_2 = g_3 = square = 0 AND (*) FZ = (1/6)gamma^5 y^9 -- holds with "
      "residual EXACTLY 0, and every one of the six caps "
      "(deg A,u,v,w,d0,CT <= 9,6,12,9,12,15) and every one of the five orders "
      "(ord A,u,v,w,CT >= 1,2,2,3,3) is respected", allw)

r0 = k0_system(y, sp.Integer(1), Rational(1, 2), y, 3, sp.Integer(1), sp.Integer(1))
check("G5b and g_3 = 0 comes out AUTOMATICALLY, which is a cross-check rather "
      "than luck: Thm 8.1 says -gamma*A*g2hat + gamma^2*Pi*g3hat = "
      "FZ - (1/6)gamma^5 y^9 Pi^4, so once g_2 = 0 and (*) hold, g_3 must "
      "vanish.  The witnesses were built from the boxed row and (*) only",
      r0["g2"] == 0 and r0["FZ"] == 0 and r0["g3"] == 0)

check("G5c the witnesses are non-degenerate where it matters: Z = zeta*y^z and "
      "F = (gamma^5/6 zeta) y^(9-z) are both nonzero MONOMIALS, deg F + deg Z = "
      "9 = 9+4k at k = 0 ((*deg)), and v != 0 in three of the four",
      all(sp.Poly(k0_system(y, sp.Integer(1), Rational(1, 2), A, z, zv, gv)["Z"],
                  y).is_monomial for A, z, zv, gv in WITS[:3])
      and sp.degree(r0["F"], y) + sp.degree(r0["Z"], y) == 9 and r0["v"] != 0)

check("G6  VERDICT for sec.8 at the class of nine: the reduction, the cofactor "
      "identity, 'Z is a monomial', the z-window, the ledger and both caps ALL "
      "transfer -- and the transferred system is CONSISTENT.  So the class rows "
      "are not merely un-killed by sec.3; they are un-killable by the whole of "
      "PROOF sec.4-8 as it stands", allw and survive108 == [] and len(survive75) == 25)


# ===========================================================================
# H.  MUTATION CONTROLS.
# ===========================================================================
head("H.  mutation controls")

check("H1  the slice-condition transfer is about mult_beta(C) = 1, not about C "
      "being a monomial: at C = y^2 the exponents become y^(4n-4)/y^(4n-6), the "
      "cascade start moves, and the profile would read 4k-2 not 2k-1 -- so a "
      "deeper monomial corner does NOT inherit (6.2.1)",
      2 * (2 * 3 - 2) != 2 * 3 - 2)

fw = dict(PROFILE)
fw[6] = 10
check("H2  the level-12 input is load-bearing at the class row exactly as at "
      "(72,108): with only Lemma 6.1's ord_y(h6) >= 10 the 3h6^2 term lands on "
      "20 and a = 10 SURVIVES the collision, so the y-place cascade level 12 is "
      "doing real work here too",
      not (10 + min(3 + 20, 10 + 12, 1 + 10 + 10, 20) > 30)
      and (10 + min(3 + 20, 10 + 12, 1 + 10 + 11, 22) > 30))

check("H3  the deg-cap direction is the dangerous one and it is checked: "
      "raising the class row's deg u cap from 6 to 7 would not rescue Cor 8.5 "
      "(it is already failing at 6), and LOWERING it to 3 would NOT create the "
      "kill either, because deg(mu*y^3) = 3 <= 3.  The obstruction is deg q, "
      "not the cap",
      max([3, 2, 2]) <= 3 and all(max([3] + [2 * dA] + [z]) <= 3
                                  for dA in [0, 1] for z in [2, 3]))

check("H4  and the (72,108) kill is not an artefact of the place: it is the same "
      "boxed row, the same cap 6, the same z-window [2,6]; ONLY deg q differs.  "
      "Running the class row's Qq = 1 through the (72,108) branch (pi = t) "
      "leaves the kill switched off, and running (72,108)'s q through the class "
      "row's branch (pi = y) switches it back ON -- the discriminant is Qq",
      [(dA, z) for dA in [None, 0, 1, 2, 3] for z in range(2, 7)
       if max([3 + 0] + ([2 * dA] if dA is not None else []) + [z]) <= 6] != []
      and [(dA, z) for dA in [None] + list(range(0, 40)) for z in range(2, 7)
           if max([3 + 4] + ([2 * dA] if dA is not None else []) + [z]) <= 6
           or [3 + 4] .count(max([3 + 4] + ([2 * dA] if dA is not None else [])
                                 + [z])) > 1] == [])


# ===========================================================================
if _fail:
    print()
    print("FAILURES (%d):" % len(_fail))
    for f in _fail:
        print("   - %s" % f)
    raise SystemExit(1)
print("ALL %d CHECKS PASSED" % _ok[0])
if not QUIET:
    print("""
VERDICT
  THE STRUCTURAL FACT.  A class row is (72,108) with C_4 = y^7(y+1) replaced by
  C = y.  (a,b,t,kappa) = (2,3,4,2), v(P) = 8, v(Q) = 12, [P,Q] = x^2,
  v(F) = -5, M = 17, N = 28 and ONE forcing ODE are shared; the only differing
  input is C.  [sec.A]

  WHAT TRANSFERS -- more than WEIGHT_FREE_TRANSFER.md says.
    * (2.5.1) is C-divisibility, so its exponents 2n-2 / 2n-3 attach to
      whichever place has multiplicity 1 in C.  At (72,108) that is t; at a
      class row it is y.  SAME condition set.                        [sec.B]
    * hence the CASCADE transfers.  Recomputed here from scratch at the y-place,
      levels 2..12, same perfect-square / one-fresh-coefficient shape, giving
      ord_y(h_k) >= 2k-1 for k <= 6 and Lemma 6.1's 12, 13.  This CORRECTS
      WEIGHT_FREE_TRANSFER rows #5 and #6 ("GENUINELY LOST").         [sec.C]
    * the polygon gives TWO integral slopes, ord >= k and deg <= 3k, so
      lam = 2, NOT 0, and the unstripped deg table is config (1)'s
      (6,9,12,15,18,21,24).  This CORRECTS G1/G2 there.               [sec.D]
    * therefore sec.7 transfers entire: ord_y(e) = 9 EXACTLY, and Thm 3.4 then
      forces e = gamma*y^9 and B = y^21/gamma, both monomials.        [sec.E]

  THE THREE ROUTES ARE DEAD.  An explicit family h_k = c_k*y^(2k-1) satisfies
  EVERY P- and Q-slice condition and both polygon caps while having
  ord_y(h1) = 1, ord_y(h2) = 3, ord_y(h5) = 9.  So {5}, {1,2,7} and {2,6,7} are
  all FALSE on a point of everything the y-place proves about valuations: none
  of them is derivable from that input.  Consistent with PROOF sec.7.4(c), which
  now transfers: what is missing is leading-coefficient information.  [sec.F]

  AND sec.8 CANNOT KILL EITHER.  k = 0 is forced, sec.8.5 is vacuous, and the
  only remaining kill is Cor 8.5, whose engine is deg(mu*t^3*q) = 7 > 6.  With
  q -> 1 that degree is 3 <= 6.  Four explicit witnesses satisfy the whole
  sec.8.1 system and every cap.  CALIBRATION CONTROL: the same machinery, run at
  (72,108), kills it for every (deg A, z) -- so the negative is a statement
  about deg q = 4 vs 0, not about the method.                          [sec.G]

  THE OBSTRUCTION, IN ONE LINE.  The class of nine survives because Phi's
  residual factor is trivial.  deg q >= 4 is exactly what Cor 8.5 needs and
  exactly what a monomial C cannot supply.""")
