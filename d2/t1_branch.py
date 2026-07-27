#!/usr/bin/env python3
"""t1_branch.py -- settle `t^a | R` on the T1 branch (d1 != 0), redo the pole /
valuation analysis on the double cover, and test which T2-only results transfer.

READ-ONLY over every pre-existing artifact.  Writes nothing.

=====================================================================
HEADLINE
=====================================================================

  `t^a | R` HOLDS ON T1.  More: the proof is BRANCH-INDEPENDENT and needs
  only two facts that are themselves branch-independent --

     (K)  the K-syzygy      2*Phi = d2*e^3 + 3*e^2*S + 3*e*R^2
     (D)  e | S             (SYZYGY_SWEEP.md sec.4; re-derived here from
                             scratch by Sylvester resultant, check C3)

  Given (D), write S = e*s with s polynomial.  Then at ANY place beta,
  with m = v_beta(e), rho = v_beta(R), P_b = v_beta(Phi):

     v(d2*e^3) >= 3m ,   v(3*e^2*S) = 2m + v(S) >= 3m ,   v(3*e*R^2) = m + 2*rho

  so if rho < m the R-term is the STRICT unique minimum and

     P_b  =  m + 2*rho   <   3m .                                   (*)

  THEOREM (place trichotomy).  On every genuine lift, at every place,
     EITHER  v_beta(R) >= v_beta(e)   OR   v_beta(Phi) = v_beta(e) + 2*v_beta(R).

  At beta = -1: P_b = 30 EXACTLY (q(-1) = 3315 != 0), m = a.  So either
  v_t(R) >= a, or a + 2*rho = 30 with rho <= a-1, which forces
  30 <= 3a - 2, i.e. a >= 11, AND a even (rho = (30-a)/2 must be an integer).

     ==>  a <= 10  or  a odd     ==>   t^a | R .

  Every sub2 state has deg e = 10, hence a <= 10.  So `t^a | R` holds on the
  WHOLE sub2 universe, T1 included.  The T2-only status is removed.

  Why the old argument was branch-bound and this one is not: the T2 proof
  (`tpower_divisibility.verify_tpower`) allowed v_t(S) to range freely over
  0..cap_S and therefore NEEDED the second relation `H3 = 0`, whose shape
  changes when d1 != 0.  `e | S` removes that need entirely: the K-syzygy
  alone does the work, and the K-syzygy has no branch hypothesis.

  SECOND, INDEPENDENT ROUTE (check C8): the spine lane's dm4-KEPT chain
  K -> G1 -> G3 reaches the same conclusion WITHOUT using e|S at all.  It is
  adjudicated here from scratch and CONFIRMED (with one recorded ablation
  disagreement -- G2 is not redundant, it is interchangeable with G1).  Two
  routes, disjoint hypotheses, same theorem.

=====================================================================
THE DOUBLE COVER
=====================================================================

  `POLE_THEOREM.md` sec.7.1 named the T1 obstruction: with d1 != 0,
  2*(H5 + d2*H3) = dm1*K5gen + 6*d1*dm2^2*dm3 with K5gen QUADRATIC in dm3,
  so dm3 is two-valued.  Computed here (checks C8/C9):

  * The discriminant of that quadratic, REDUCED MODULO THE G-IDEAL, is
        9*(2*d1*(e*S - R^2) + e^3)^2 ,
    a PERFECT SQUARE.  The K5gen cover is REDUCIBLE -- the two-valuedness is
    an artifact of not having used the K-syzygy.  On the lift the quadratic
    is exactly  (3*d1*e) * P2, where P2 is the G2 relation.

  * The genuine double cover is P2 itself, in s = S/e:
        e^2*s^2 - 2*R^2*s - (d0*e^2 + d1*e*R + d2*R^2) = 0
    with discriminant (over 4)
        DELTA  =  R^4 + d2*e^2*R^2 + d1*e^3*R + d0*e^4  =  e^4 * P(R/e),
        P(X)   =  X^4 + d2*X^2 + d1*X + d0 .
    Setting u = R/e and v = s - u^2 = (e*S - R^2)/e^2 this is EXACTLY
        v^2  =  P(u)  =  u^4 + d2*u^2 + d1*u + d0 ,
    a genus-<=1 quartic model.  Every G-system lift is a Q(y)-point of it.

  * FIELD-SCOPE-SAFE COROLLARY (no square classes over Q are used):
    W := e*S - R^2 is a polynomial with W^2 = DELTA, so
        every root of  R^4 + d2*e^2*R^2 + d1*e^3*R + d0*e^4
        has EVEN multiplicity,
    and e^2 | (R^2 + W).  This is new, branch-independent, and cap-free.

  * At a marked root beta (root of q dividing e): m = 1 and rho_b = 0 are
    forced by (*), so DELTA(beta) = R(beta)^4 != 0 -- the cover is UNRAMIFIED
    at every marked root, and the branch is pinned by e^2 | R^2 + W.

=====================================================================
WHAT TRANSFERS / WHAT DOES NOT  (details in T1_BRANCH.md)
=====================================================================
  TRANSFERS (branch-independent):  Theorem 2A, Theorem 2B (incl. the
    square relation 2*Phi'(beta) = 3*e'(beta)*R(beta)^2), Theorem 2C
    (v_t(dm2) >= a AND v_t(dm3) >= a), the sec.4 spare collapse 28 -> 13-a.
  DOES NOT TRANSFER:  R | e^2  and  e*R | Phi.  And the failure is PERMANENT
    for the divisor method: the explicit point (check C12)
        e = 1, R = 0, S = 1, T = 1/6, d0 = 1, d1 = -1/3, d2 = 0, Phi = 3/2
    lies on V(G1,G2,G3) with R = 0, e != 0, Phi != 0 and d1 != 0, so no power
    of e and no power of Phi lies in I + (R) or I + (e*R).  On T2 no such
    point exists (R = 0 & d1 = 0 forces e = 0), which is exactly why those
    two results are T2-only.

COLON-IDEAL DISCIPLINE: this file computes NO colon ideal of any kind.  Every
certificate is an explicit polynomial identity verified by sympy expansion.

Usage:  python -u t1_branch.py            # full report
        python -u t1_branch.py --quiet    # self-check, exit 0/1
"""
from __future__ import annotations

import argparse
import itertools
import json
import os
import sys

import sympy as sp

import face_kill_sweep as fks
import full_system_bridge as fsb

HERE = os.path.dirname(os.path.abspath(__file__)) or "."
INF = 10 ** 9                       # valuation of the zero polynomial

Y = sp.Symbol("y")
Q_POLY = 2048 * Y ** 4 - 512 * Y ** 3 + 320 * Y ** 2 - 240 * Y + 195
V_T_PHI = 30                        # v_{y+1}(Phi), exact: q(-1) = 3315 != 0

d0, d1, d2, e, R, S, T, Phi = fks._gsystem_symbols()
s_ = sp.Symbol("s_")                # s = S/e
W_ = sp.Symbol("W_")                # W = e*S - R^2

RESULTS: list[tuple[str, bool, str]] = []


def record(tag, ok, msg=""):
    RESULTS.append((tag, bool(ok), msg))
    return bool(ok)


def G():
    return {k: sp.expand(v[0]) for k, v in fks.canonical_G_generators().items()}


# =====================================================================
#  C1 / C2 -- the standing guard and the K-syzygy
# =====================================================================
def c1_guard(v=True):
    g = G()
    ok = sp.expand(g["G5"]).coeff(Phi) == 1
    if v:
        print("[C1] standing guard  coeff(G5, Phi) == 1 :", ok)
    return record("C1_G5_normalisation", ok)


def c2_ksyzygy(v=True):
    g = G()
    lhs = sp.expand(2 * (g["G5"] + d2 * g["G3"] + d1 * g["G2"] + d0 * g["G1"]))
    rhs = sp.expand(2 * Phi - e * (d2 * e ** 2 + 3 * e * S + 3 * R ** 2))
    resid = sp.expand(lhs - rhs)
    ok = resid == 0
    if v:
        print("[C2] K-syzygy 2*(G5+d2*G3+d1*G2+d0*G1) == 2*Phi - e*(d2*e^2+3*e*S+3*R^2)")
        print("     residual =", resid, "-> EXACT" if ok else "*** NONZERO ***")
        print("     => on any lift  2*Phi = d2*e^3 + 3*e^2*S + 3*e*R^2   (branch-free)")
    return record("C2_K_syzygy", ok)


def c2b_phi(v=True):
    from bigrade_annotator import _phi_stripped
    ph = sp.expand(_phi_stripped())
    ratio = sp.cancel(ph / sp.expand((Y + 1) ** 30 * Q_POLY))
    sqfree = sp.gcd(sp.Poly(Q_POLY, Y), sp.Poly(sp.diff(Q_POLY, Y), Y)).degree() == 0
    qm1 = Q_POLY.subs(Y, -1)
    ok = (ratio == sp.Rational(-1, 6630)) and sqfree and qm1 != 0
    if v:
        print("[C2b] Phi = %s * (y+1)^30 * q ; q squarefree: %s ; q(-1) = %s"
              % (ratio, sqfree, qm1))
        print("      => v_t(Phi) = 30 EXACTLY, v_beta(Phi) = 1 at each root of q")
    return record("C2b_phi_shape", ok)


# =====================================================================
#  C3 -- `e | S` re-derived from scratch (Sylvester resultant + adjugate)
# =====================================================================
def c3_e_divides_S(v=True):
    """Independent re-derivation of e | S.  NO colon ideal is used."""
    g = G()
    A = sp.expand(R * g["G1"] - e * g["G2"])
    B = sp.expand(S * g["G1"] - e * g["G3"])
    ok_T = (not A.has(T)) and (not B.has(T))

    pa, pb = sp.Poly(A, R), sp.Poly(B, R)
    a2, a1, a0 = pa.all_coeffs()
    b2, b1, b0 = pb.all_coeffs()
    M = sp.Matrix([[a2, a1, a0, 0], [0, a2, a1, a0],
                   [b2, b1, b0, 0], [0, b2, b1, b0]])
    Res = sp.expand(M.det())
    # adjugate: last row of adj(M) gives Res = c1*(R*A) + c2*A + c3*(R*B) + c4*B
    adj = M.adjugate()
    c1, c2, c3, c4 = [sp.expand(adj[3, j]) for j in range(4)]
    u = sp.expand(c1 * R + c2)
    w = sp.expand(c3 * R + c4)
    cof_ok = sp.expand(u * A + w * B - Res) == 0

    # ... hence Res = (u*R + w*S)*G1 + (-u*e)*G2 + (-w*e)*G3  in I
    ideal_ok = sp.expand((u * R + w * S) * g["G1"] - u * e * g["G2"]
                         - w * e * g["G3"] - Res) == 0

    Qq = sp.expand(sp.cancel(Res / (sp.Rational(729, 16) * e)))
    div_ok = sp.expand(Res - sp.Rational(729, 16) * e * Qq) == 0
    pq = sp.Poly(Qq, S)
    deg_ok = pq.degree() == 7 and pq.all_coeffs()[0] == sp.Rational(-8, 9)

    # normalise to monic in S and read off the alpha_i; each must be e^i * poly
    Qm = sp.expand(sp.Rational(-9, 8) * Qq)
    cs = sp.Poly(Qm, S).all_coeffs()
    ok_monic = cs[0] == 1
    alphas, poly_ok = [], True
    for i in range(1, 8):
        ci = sp.expand(cs[i])
        # exact polynomial division by e^i in Q[d0,d1,d2,e]
        quo, rem = sp.div(sp.Poly(ci, e), sp.Poly(e ** i, e))
        if sp.expand(rem.as_expr()) != 0:
            poly_ok = False
        ai = sp.expand(quo.as_expr())
        if sp.together(ai).has(sp.Pow(e, -1)) or sp.expand(ai * e ** i - ci) != 0:
            poly_ok = False
        alphas.append(ai)
    poly_ok = poly_ok and ok_monic
    ok = ok_T and cof_ok and ideal_ok and div_ok and deg_ok and poly_ok
    if v:
        print("[C3] e | S, re-derived here (no colon ideal, no trust in resultant theory)")
        print("     A := R*G1 - e*G2 , B := S*G1 - e*G3 : dm4-free =", ok_T,
              ", deg_R = 2 each")
        print("     Sylvester Res_R(A,B) = u*A + w*B with u,w from the ADJUGATE:", cof_ok)
        print("     Res = (u*R+w*S)*G1 - u*e*G2 - w*e*G3   (so Res in I):", ideal_ok)
        print("     Res = (729/16)*e*Q ,  deg_S Q = 7, lc = -8/9 :", div_ok and deg_ok)
        print("     monic form  S^7 + sum_i e^i*alpha_i*S^(7-i) = 0, alpha_i polynomial:",
              poly_ok)
        for i, a in enumerate(alphas, 1):
            print("        alpha_%d = %s" % (i, sp.factor(a)))
        print("     => at any root p of e: 7*sigma >= min_i(i*mu + (7-i)*sigma)")
        print("        forces sigma >= mu, i.e. e | S.  BOTH BRANCHES, cap-free.")
    return record("C3_e_divides_S", ok)


# =====================================================================
#  C4 -- the reduced system after S = e*s and the T-formula
# =====================================================================
def reduced_system():
    g = G()
    sub = {S: e * s_, T: -R * (s_ + d2) - d1 * e / 2}
    P1 = sp.expand(g["G1"].subs(sub))
    P2 = sp.expand(sp.Rational(-2, 3) * g["G2"].subs(sub))
    P3 = sp.expand(-2 * g["G3"].subs(sub))
    P5 = sp.expand(-2 * g["G5"].subs(sub))
    return P1, P2, P3, P5


def c4_reduction(v=True):
    P1, P2, P3, P5 = reduced_system()
    ok1 = P1 == 0
    P2_want = sp.expand(d0 * e ** 2 + d1 * e * R + d2 * R ** 2 - e ** 2 * s_ ** 2
                        + 2 * R ** 2 * s_)
    P3_want = sp.expand(6 * d0 * e * R + 3 * d1 * e ** 2 * s_ + 3 * d1 * R ** 2
                        + 6 * d2 * e * R * s_ + e ** 3 + 6 * e * R * s_ ** 2)
    ok2 = sp.expand(P2 - P2_want) == 0
    ok3 = sp.expand(P3 - P3_want) == 0
    # the K-syzygy in reduced form: 2*Phi = e^3*(d2 + 3*s) + 3*e*R^2
    Kred = sp.expand(2 * Phi - e ** 3 * (d2 + 3 * s_) - 3 * e * R ** 2)
    okK = sp.expand(Kred - (2 * Phi - e * (d2 * e ** 2 + 3 * e * (e * s_) + 3 * R ** 2))) == 0
    ok = ok1 and ok2 and ok3 and okK
    if v:
        print("[C4] reduced system (S = e*s, T = -R*(s+d2) - d1*e/2) -- branch-free")
        print("     G1 vanishes IDENTICALLY:", ok1, "  (dm4 is determined, not a spare)")
        print("     P2 := -(2/3)*G2| =", P2)
        print("     P3 := -2*G3|     =", P3)
        print("     K   : 2*Phi = e^3*(d2+3*s) + 3*e*R^2 :", okK)
    return record("C4_reduced_system", ok)


# =====================================================================
#  C5 -- THE THEOREM: the place trichotomy and t^a | R
# =====================================================================
def place_trichotomy(m, rho, v_s, delta2, P_b):
    """Is the K-syzygy 2*Phi = d2*e^3 + 3*e^2*S + 3*e*R^2 satisfiable at these
    orders?  v(S) = m + v_s with v_s >= 0 because e | S."""
    orders = [delta2 + 3 * m, 3 * m + v_s, m + 2 * rho, P_b]
    orders = [o for o in orders if o < INF]
    if not orders:
        return True
    mn = min(orders)
    return orders.count(mn) >= 2


def c5_theorem(v=True, amax=10):
    """Machine-check the closed-form argument at the t-place for a = 0..amax
    (and report the exceptional even a >= 11), by exhaustive enumeration."""
    bad = {}
    for a in range(0, amax + 1):
        surviving = []
        for rho in range(0, a):                      # rho < a is what we refute
            for v_s in list(range(0, 3 * a + 40)) + [INF]:
                for delta2 in list(range(0, 3 * a + 40)) + [INF]:
                    if place_trichotomy(a, rho, v_s, delta2, V_T_PHI):
                        surviving.append((rho, v_s, delta2))
                        break
                if surviving:
                    break
            if surviving:
                break
        if surviving:
            bad[a] = surviving[0]
    ok = not bad

    # the exceptional regime, closed form: a + 2*rho = 30, rho <= a-1
    exc = [a for a in range(11, 40)
           if (30 - a) % 2 == 0 and 0 <= (30 - a) // 2 <= a - 1]
    exc_check = sorted(a for a in range(11, 40)
                       if any(place_trichotomy(a, rho, v_s, delta2, V_T_PHI)
                              for rho in range(0, a)
                              for v_s in [0, 1, 5, 100, INF]
                              for delta2 in [0, 1, 5, 100, INF]))
    ok_exc = exc == exc_check
    # ADMISSIBILITY CONTROL: rho >= a MUST remain feasible, or the test is vacuous
    ctrl = [(a, rho, v_s, dl2)
            for a in (7, 9, 10) for rho in (a, a + 1)
            for v_s in (0, 1, 2, 30 - 3 * a) for dl2 in (0, 30 - 3 * a)
            if v_s >= 0 and dl2 >= 0 and place_trichotomy(a, rho, v_s, dl2, V_T_PHI)]
    ok_ctrl = len(ctrl) > 0
    if v:
        print("[C5] THEOREM  t^a | R  --  exhaustive refutation of v_t(R) < a")
        print("     inputs: K-syzygy + (e | S, so v_t(S) >= a) + d2 polynomial")
        print("             + v_t(Phi) = 30 exactly.   NO branch hypothesis.")
        print("     a = 0..%d : v_t(R) >= a  forced :" % amax, ok,
              "" if ok else "  *** survivor %s ***" % bad)
        print("     exceptional regime (closed form a+2rho=30, rho<a): a in", exc)
        print("     machine agrees:", ok_exc)
        print("     -> a even and 12 <= a <= 30 are the ONLY escapes; every odd a")
        print("        and every a <= 10 is forced.  deg e <= 10 on sub2 => a <= 10.")
        print("     admissibility control (rho >= a must stay feasible):",
              "%d feasible tuples, e.g. (a,rho,v_s,delta2) = %s"
              % (len(ctrl), ctrl[0]) if ok_ctrl else "*** VACUOUS TEST ***")
    return record("C5_tpower_R_on_T1", ok and ok_exc and ok_ctrl)


def _p2_ok(a, rho, v_s, e0, e1, e2):
    """P2 = d0*e^2 + d1*e*R + d2*R^2 - e^2*s^2 + 2*R^2*s = 0 ; min attained >= 2x."""
    o = [e0 + 2 * a, e1 + a + rho, e2 + 2 * rho, 2 * a + 2 * v_s, 2 * rho + v_s]
    o = [x for x in o if x < INF]
    if not o:
        return True
    mn = min(o)
    return o.count(mn) >= 2


def _p3_ok(a, rho, v_s, e0, e1, e2):
    """P3 = 6d0eR + 3d1e^2 s + 3d1R^2 + 6d2eRs + e^3 + 6eRs^2 = 0 ; min >= 2x.
    The e^3 term is always present (e != 0), so the all-infinite case cannot occur."""
    o = [e0 + a + rho, e1 + 2 * a + v_s, e1 + 2 * rho,
         e2 + a + rho + v_s, 3 * a, a + rho + 2 * v_s]
    o = [x for x in o if x < INF]
    mn = min(o)
    return o.count(mn) >= 2


def c6_t2_consistency(v=True):
    """Cross-check against the RECORDED T2 route.  `tpower_divisibility` proves
    t^a|R on T2 WITHOUT e|S, using H3; this file proves it on BOTH branches
    WITH e|S, using the K-syzygy alone.  Both must agree on T2 at a = 7, 8, 9."""
    import tpower_divisibility as tpd
    old = {a: bool(tpd.verify_tpower(a, verbose=False)[0]) for a in (7, 8, 9)}
    new = {a: not any(place_trichotomy(a, rho, v_s, dl2, V_T_PHI)
                      for rho in range(0, a)
                      for v_s in list(range(0, 3 * a + 40)) + [INF]
                      for dl2 in list(range(0, 3 * a + 40)) + [INF])
           for a in (7, 8, 9)}
    ok = all(old[a] and new[a] for a in (7, 8, 9))
    if v:
        print("[C6] consistency with the recorded T2 route (tpower_divisibility)")
        print("     old (T2 only, H3-based, v_t(S) free)  v_t(R) >= a :", old)
        print("     new (both branches, K-only, uses e|S)             :", new)
        print("     the new route uses a STRICT SUBSET of the hypotheses and")
        print("     drops the d1 = 0 assumption entirely.")
    return record("C6_T2_consistency", ok)


def c7_escapes(v=True):
    """The FRONTIER_REBUILD.md:273 lead: does a_t <= 10 now follow?

    K alone: if rho < a then a+2rho <= 3a-2 < 3a <= min(delta2+3a, 3a+v_s), so
    the R-term is the strict unique minimum and a + 2*rho = v_t(2*Phi) = 30.
    The ONLY candidate rho is (30-a)/2, and only when that is an integer < a."""
    rng = list(range(0, 61)) + [INF]
    rows = []
    for a in range(11, 16):
        cand = [rho for rho in range(0, a)
                if any(place_trichotomy(a, rho, v_s, dl2, V_T_PHI)
                       for v_s in rng for dl2 in rng)]
        alive = []
        for rho in cand:
            hit = None
            for v_s in rng:
                for e0 in rng:
                    for e1 in rng:            # e1 = INF is the T2 sub-case
                        for e2 in rng:
                            if (place_trichotomy(a, rho, v_s, e2, V_T_PHI)
                                    and _p2_ok(a, rho, v_s, e0, e1, e2)
                                    and _p3_ok(a, rho, v_s, e0, e1, e2)):
                                hit = (v_s, e0, e1, e2)
                                break
                        if hit:
                            break
                    if hit:
                        break
                if hit:
                    break
            if hit:
                alive.append((rho, hit))
        rows.append((a, cand, alive))
    closed = {11: [], 12: [9], 13: [], 14: [8], 15: []}
    ok = all([r for r, _ in al] == closed[a] for a, cand, al in rows) \
        and all(cand == closed[a] for a, cand, al in rows)
    if v:
        print("[C7] the a_t <= 10 lead (FRONTIER_REBUILD.md:273)")
        print("     if v_t(R) >= a then 30 = a + v_t(bracket) >= 3a, so a <= 10.")
        print("     the only alternative is a + 2*v_t(R) = 30 with v_t(R) < a:")
        for a, cand, al in rows:
            print("        a = %2d : K-candidate rho %-6s ; survives K+P2+P3 -> %s"
                  % (a, cand or "NONE", [r for r, _ in al] or "NONE"))
        print("     => a_t in {11, 13, 15} DELETED (parity: (30-a)/2 not an integer)")
        print("        a_t in {12, 14} SURVIVE, with v_t(R) = 9, 8 respectively.")
        print("     Lead PARTIALLY unblocked: 3 of the 5 alternate-regime values die.")
    return record("C7_at_le_10_lead", ok)


# =====================================================================
#  C8 -- ADJUDICATION of the spine lane's dm4-KEPT argument
# =====================================================================
#  The spine lane claims `t^a | dm2, dm3, dm4` on BOTH branches, a = 6..10,
#  from G1 and G3 with dm4 KEPT (never from the dm4-eliminated H3, which is
#  where d1 = 0 entered `tpower_divisibility`).  Its chain:
#      rho < a  =>  K forces  s = 2*rho - a
#               =>  G1 forces tau = 3*rho - 2*a
#               =>  G3's 3*S*T term is a STRICT UNIQUE minimum, so G3 != 0.
#  Adjudicated here from scratch, two ways.
# =====================================================================
def _orders_K(a, rho, s, tau, e0, e1, e2, P_b=V_T_PHI):
    """2*Phi = d2*e^3 + 3*e^2*S + 3*e*R^2  as a 4-term vanishing sum."""
    return [e2 + 3 * a, 2 * a + s, a + 2 * rho, P_b]


def _orders_G1(a, rho, s, tau, e0, e1, e2):
    """G1 = (3/2)*d1*e^2 + 3*d2*e*R + 3*e*T + 3*R*S   (dm4 KEPT)."""
    return [e1 + 2 * a, e2 + a + rho, a + tau, rho + s]


def _orders_G2(a, rho, s, tau, e0, e1, e2):
    """G2 = -(3/2)*d0*e^2 + (3/2)*d2*R^2 + 3*R*T + (3/2)*S^2."""
    return [e0 + 2 * a, e2 + 2 * rho, rho + tau, 2 * s]


def _orders_G3(a, rho, s, tau, e0, e1, e2):
    """G3 = -3*d0*e*R - (3/2)*d1*R^2 - (1/2)*e^3 + 3*S*T   (dm4 KEPT).
    The e^3 term is always present (e != 0)."""
    return [e0 + a + rho, e1 + 2 * rho, 3 * a, s + tau]


def _sum_ok(orders):
    o = [x for x in orders if x < INF]
    if not o:
        return True
    mn = min(o)
    return o.count(mn) >= 2


ROWS = {"K": _orders_K, "G1": _orders_G1, "G2": _orders_G2, "G3": _orders_G3}


def _spine_hand_proof(amax=10):
    """The spine chain as pure arithmetic assertions, for every (a, rho) with
    rho < a <= amax and every delta_i >= 0.  Returns (ok, failures)."""
    fails = []
    for a in range(1, amax + 1):
        for rho in range(0, a):
            # (i) K forces s = 2*rho - a
            if not (a + 2 * rho < 3 * a):                    # < delta2 + 3a
                fails.append(("i-A", a, rho))
            if not (a + 2 * rho < V_T_PHI):                  # < v_t(2*Phi)
                fails.append(("i-B", a, rho))
            s = 2 * rho - a
            # v(S) >= 0 is required; s < 0 is itself a refutation, so skip
            if s < 0:
                continue
            # (ii) G1 forces tau = 3*rho - 2*a; rho+s must beat the d1 and d2 terms
            if not (rho + s < 2 * a):                        # < delta1 + 2a  (d1!)
                fails.append(("ii-d1", a, rho))
            if not (rho + s < a + rho):                      # < delta2 + a + rho
                fails.append(("ii-d2", a, rho))
            tau = 3 * rho - 2 * a
            if tau < 0:
                continue
            # (iii) G3: s + tau = 5*rho - 3*a is the STRICT unique minimum
            if not (s + tau < a + rho):                      # < delta0 + a + rho
                fails.append(("iii-d0", a, rho))
            if not (s + tau < 2 * rho):                      # < delta1 + 2*rho  (d1!)
                fails.append(("iii-d1", a, rho))
            if not (s + tau < 3 * a):                        # < e^3
                fails.append(("iii-e3", a, rho))
    return (not fails), fails


def _spine_enumerate(amax=10, LIM=36, rows=("K", "G1", "G2", "G3"), t2=False):
    """Exhaustive dm4-KEPT valuation enumeration.  Returns the first surviving
    configuration with v_t(R) < a, or None."""
    grid = list(range(0, LIM + 1)) + [INF]
    d1rng = [INF] if t2 else grid
    fs = [ROWS[r] for r in rows]
    for a in range(1, amax + 1):
        for rho in range(0, a):
            for s in grid:
                for e2 in grid:
                    if "K" in rows and not _sum_ok(_orders_K(a, rho, s, 0, 0, 0, e2)):
                        continue
                    for e1 in d1rng:
                        for tau in grid:
                            if "G1" in rows and not _sum_ok(
                                    _orders_G1(a, rho, s, tau, 0, e1, e2)):
                                continue
                            for e0 in grid:
                                if all(_sum_ok(f(a, rho, s, tau, e0, e1, e2))
                                       for f in fs):
                                    return dict(a=a, rho=rho, s=s, tau=tau,
                                                d0=e0, d1=e1, d2=e2)
    return None


def c8_spine_adjudication(v=True):
    hand_ok, fails = _spine_hand_proof(10)
    enum_t1 = _spine_enumerate(10, 36, t2=False)
    enum_t2 = _spine_enumerate(10, 36, t2=True)
    enum_cap = _spine_enumerate(10, 44, t2=False)          # cap independence
    # ablation: drop one row at a time (and the {G1,G2} pair) and see whether a
    # survivor appears
    abl = {}
    for drop in ("K", "G1", "G2", "G3"):
        rows = tuple(r for r in ("K", "G1", "G2", "G3") if r != drop)
        abl[drop] = _spine_enumerate(10, 24, rows=rows, t2=False)
    abl["G1+G2"] = _spine_enumerate(10, 24, rows=("K", "G3"), t2=False)
    # admissibility control: rho >= a MUST be feasible (else the test is vacuous)
    ctrl = None
    grid = list(range(0, 37)) + [INF]
    for a in (7, 9, 10):
        for rho in range(a, a + 3):
            for s in grid:
                for tau in grid:
                    if all(_sum_ok(f(a, rho, s, tau, 0, 0, 0))
                           for f in ROWS.values()):
                        ctrl = dict(a=a, rho=rho, s=s, tau=tau)
                        break
                if ctrl:
                    break
            if ctrl:
                break
        if ctrl:
            break
    ok = (hand_ok and enum_t1 is None and enum_t2 is None and enum_cap is None
          and ctrl is not None
          and abl["K"] is not None and abl["G3"] is not None
          and abl["G1"] is None and abl["G2"] is None
          and abl["G1+G2"] is not None)
    if v:
        print("[C8] ADJUDICATION -- the spine lane's dm4-KEPT chain, re-derived here")
        print("     hand chain as arithmetic assertions, all (a,rho) with rho<a<=10:",
              hand_ok, "" if hand_ok else fails[:4])
        print("       (i)   K  : a+2rho < 3a and a+2rho < 30  =>  v(S) = 2rho - a")
        print("       (ii)  G1 : rho+s < delta1+2a and < delta2+a+rho => v(T) = 3rho-2a")
        print("             ^^ the d1 term enters ONLY as delta1 >= 0, and the")
        print("                inequality 3rho-a < 2a is strict already at delta1 = 0")
        print("       (iii) G3 : s+tau = 5rho-3a < min(delta0+a+rho, delta1+2rho, 3a),")
        print("                all three strict for rho < a  =>  G3 != 0.  CONFIRMED.")
        print("     exhaustive dm4-KEPT enumeration, T1 (d1 free), a<=10, LIM=36:",
              "no survivor" if enum_t1 is None else "*** %s ***" % enum_t1)
        print("     same with d1 = 0 (T2):",
              "no survivor" if enum_t2 is None else "*** %s ***" % enum_t2)
        print("     cap independence (LIM=44):",
              "no survivor" if enum_cap is None else "*** %s ***" % enum_cap)
        print("     admissibility control (rho >= a MUST be feasible):", ctrl)
        print("     row ablation (a survivor appearing = that row was load-bearing):")
        for k, sv in abl.items():
            print("        drop %-5s -> %s"
                  % (k, "survivor %s" % sv if sv else "still empty"))
        print("     DISAGREEMENT WITH THE SPINE LANE, recorded loudly: it reports")
        print("       `G1, G3, K each load-bearing; G2 redundant`.  Here G1 and G2")
        print("       are EACH individually redundant and only their PAIR is")
        print("       load-bearing -- G2 pins v(T) exactly as G1 does, because")
        print("       2*s = 4rho-2a is strictly below both delta0+2a and delta2+2rho")
        print("       when rho < a, so G2 too forces rho + tau = 2*s, tau = 3rho-2a.")
        print("       This STRENGTHENS the lane's conclusion; it does not weaken it.")
        print("     VERDICT: the spine lane's argument is CORRECT.  d1 never enters")
        print("       a step; both d1 occurrences are lower bounds that the strict")
        print("       inequalities clear at delta1 = 0.  Independent of this file's")
        print("       own route (K + e|S), which reaches the same conclusion.")
    return record("C8_spine_adjudication", ok)


# =====================================================================
#  C8 / C9 -- the double cover
# =====================================================================
def c8_double_cover(v=True):
    K5gen = (2 * Phi + 3 * d0 * d1 * e ** 2 + 3 * d1 ** 2 * e * R + 3 * d1 * d2 * R ** 2
             - 3 * d1 * S ** 2 - d2 * e ** 3 - 3 * e ** 2 * S - 3 * e * R ** 2)
    g = G()
    H3 = sp.expand(e * g["G3"] - S * g["G1"])
    H5 = sp.expand(e * g["G5"] + (d0 * e + d1 * R + d2 * S) * g["G1"])
    ok_thm1p = sp.expand(2 * (H5 + d2 * H3) - (e * K5gen + 6 * d1 * R ** 2 * S)) == 0

    quad = sp.expand(e * K5gen + 6 * d1 * R ** 2 * S)
    pq = sp.Poly(quad, S)
    c2_, c1_, c0_ = pq.all_coeffs()
    disc = sp.expand(c1_ ** 2 - 4 * c2_ * c0_)
    # reduce modulo the K-syzygy (2*Phi -> d2 e^3 + 3 e^2 S + 3 e R^2) and S -> e*s
    lift = {Phi: (d2 * e ** 3 + 3 * e ** 2 * (e * s_) + 3 * e * R ** 2) / 2, S: e * s_}
    disc_l = sp.expand(disc.subs(lift))
    Wl = sp.expand(e ** 2 * s_ - R ** 2)
    P2 = reduced_system()[1]
    ok_sq = sp.expand(disc_l - 9 * (2 * d1 * Wl + e ** 3) ** 2 - 36 * d1 ** 2 * e ** 2 * P2) == 0
    # and the quadratic itself is -3*d1*e*P2 on the lift
    quad_l = sp.expand(quad.subs(lift))
    ok_quad = sp.expand(quad_l - 3 * d1 * e * P2) == 0
    ok = ok_thm1p and ok_sq and ok_quad
    if v:
        print("[C9] the K5gen double cover SPLITS")
        print("     Thm 1' 2*(H5+d2*H3) = e*K5gen + 6*d1*R^2*S :", ok_thm1p)
        print("     quadratic in dm3: %s*S^2 + %s*S + (...)"
              % (sp.factor(c2_), sp.factor(c1_)))
        print("     on the lift (K-syzygy used) it equals (3*d1*e)*P2 :", ok_quad)
        print("     its discriminant = 9*(2*d1*(e*S-R^2) + e^3)^2 + 36*d1^2*e^2*P2")
        print("       i.e. a PERFECT SQUARE modulo the ideal :", ok_sq)
        print("     => dm3 is NOT genuinely two-valued; the cover is reducible.")
    return record("C9_K5gen_cover_splits", ok)


def c10_genus1(v=True):
    P2 = reduced_system()[1]
    disc = sp.expand(sp.discriminant(P2, s_))
    Delta = sp.expand(R ** 4 + d2 * e ** 2 * R ** 2 + d1 * e ** 3 * R + d0 * e ** 4)
    ok_disc = sp.expand(disc - 4 * Delta) == 0
    Wl = sp.expand(e ** 2 * s_ - R ** 2)
    ok_W = sp.expand(Wl ** 2 - Delta + e ** 2 * P2) == 0
    u, vv = sp.symbols("u v")
    chk = sp.expand((P2 / e ** 2).subs(s_, vv + (R / e) ** 2).subs(R, u * e))
    ok_curve = sp.simplify(chk - (d0 + d1 * u + d2 * u ** 2 + u ** 4 - vv ** 2)) == 0
    ok = ok_disc and ok_W and ok_curve
    if v:
        print("[C10] the GENUINE double cover")
        print("     P2 is quadratic in s;  disc_s(P2)/4 = DELTA =")
        print("        R^4 + d2*e^2*R^2 + d1*e^3*R + d0*e^4  =  e^4 * P(R/e) :", ok_disc)
        print("        P(X) = X^4 + d2*X^2 + d1*X + d0")
        print("     W := e*S - R^2 satisfies W^2 = DELTA modulo P2 :", ok_W)
        print("     u = R/e, v = s - u^2  =>  v^2 = P(u) exactly :", ok_curve)
        print("     COROLLARY (field-scope safe): every root of DELTA has EVEN")
        print("        multiplicity, and e^2 | (R^2 + W).")
    return record("C10_genus1_cover", ok)


def c11_places(v=True):
    """Theorems 2A / 2B, re-derived branch-independently from (*)"""
    # 2A: P_b = 0, m >= 1  -> impossible
    a2 = not any(place_trichotomy(m, rho, v_s, dl2, 0)
                 for m in range(1, 8) for rho in range(0, 12)
                 for v_s in list(range(0, 12)) + [INF]
                 for dl2 in list(range(0, 12)) + [INF])
    # 2B: P_b = 1 -> m = 1 and rho = 0 are the only survivors
    surv = set()
    for m in range(1, 8):
        for rho in list(range(0, 12)) + [INF]:
            for v_s in list(range(0, 20)) + [INF]:
                for dl2 in list(range(0, 20)) + [INF]:
                    if place_trichotomy(m, rho, v_s, dl2, 1):
                        surv.add((m, rho))
    b2 = surv == {(1, 0)}
    # the square relation, symbolically: differentiate
    #   2*Phi = d2*e^3 + 3*e^3*s + 3*e*R^2   and set e(beta) = 0.
    yv = sp.Symbol("y")
    ee, RR, ss, dd2, PP = (sp.Function(f)(yv) for f in ("e", "R", "s", "d2", "Phi"))
    expr = 2 * PP - dd2 * ee ** 3 - 3 * ee ** 3 * ss - 3 * ee * RR ** 2
    dexpr = sp.expand(sp.diff(expr, yv))
    ep, Php = sp.symbols("ep Php")
    at_beta = {sp.Derivative(ee, yv): ep, sp.Derivative(PP, yv): Php}
    dexpr0 = sp.expand(dexpr.doit().subs(at_beta).subs(ee, 0))
    ok_rel = sp.expand(dexpr0 - (2 * Php - 3 * ep * RR ** 2)) == 0
    ok = a2 and b2 and ok_rel
    if v:
        print("[C11] Theorems 2A / 2B transfer to T1 (branch-free, from (*) alone)")
        print("     2A  root of e with v_beta(Phi) = 0 is impossible :", a2)
        print("         => rad(e) | (y+1)*q  on BOTH branches")
        print("     2B  at a root of q: only (m, rho) = (1, 0) survives :", b2)
        print("         => v_beta(e) <= 1, and R(beta) != 0")
        print("     2B' differentiating 2*Phi = d2*e^3+3*e^2*S+3*e*R^2 at e(beta)=0:")
        print("         2*Phi'(beta) = 3*e'(beta)*R(beta)^2 :", ok_rel)
    return record("C11_place_theorems", ok)


# =====================================================================
#  C12 -- the clean NEGATIVES: R | e^2 and e*R | Phi do NOT transfer
# =====================================================================
def c12_negatives(v=True):
    g = G()
    pt = {e: sp.Integer(1), R: sp.Integer(0), S: sp.Integer(1),
          T: sp.Rational(1, 6), d0: sp.Integer(1), d1: sp.Rational(-1, 3),
          d2: sp.Integer(0)}
    vals = {k: sp.expand(g[k].subs(pt)) for k in ("G1", "G2", "G3")}
    phi_val = sp.solve(sp.Eq(sp.expand(g["G5"].subs(pt)), 0), Phi)[0]
    on_V = all(vv == 0 for vv in vals.values())
    ok_T1 = pt[d1] != 0 and pt[e] != 0 and pt[R] == 0 and phi_val != 0

    # the T2 analogue does NOT exist: R = 0 and d1 = 0 force e = 0
    sols = sp.solve([sp.expand(g[k].subs({R: 0, d1: 0})) for k in ("G1", "G3")],
                    [T, e], dict=True)
    t2_forces = all(sol.get(e, None) == 0 for sol in sols) if sols else None
    # direct: G1|_{R=0,d1=0} = 3*e*T ; G3|_{R=0,d1=0} = -e^3/2 + 3*S*T
    # e != 0 => T = 0 => e^3 = 0.  Assert symbolically.
    g1r = sp.expand(g["G1"].subs({R: 0, d1: 0}))
    g3r = sp.expand(g["G3"].subs({R: 0, d1: 0}))
    ok_t2 = sp.expand(g1r - 3 * e * T) == 0 and sp.expand(g3r - (-e ** 3 / 2 + 3 * S * T)) == 0

    # the T2-only identities, and their exact T1 replacements
    P2, P3 = reduced_system()[1], reduced_system()[2]
    t2_Re2 = sp.expand(P3.subs(d1, 0) - (6 * R * (d0 * e + e * s_ * (s_ + d2)) + e ** 3))
    ok_Re2 = t2_Re2 == 0
    # T1 replacement: P3 = R*(6*d0*e + 3*d1*R + 6*e*s*(s+d2)) + e^2*(e + 3*d1*s)
    t1_repl = sp.expand(P3 - (R * (6 * d0 * e + 3 * d1 * R + 6 * e * s_ * (s_ + d2))
                              + e ** 2 * (e + 3 * d1 * s_)))
    ok_repl = t1_repl == 0
    # and the P2 replacement: R | e^2*(s^2 - d0)
    t1_repl2 = sp.expand(P2 - (e ** 2 * (d0 - s_ ** 2) + R * (d1 * e + d2 * R + 2 * R * s_)))
    ok_repl2 = t1_repl2 == 0

    ok = on_V and ok_T1 and ok_t2 and ok_Re2 and ok_repl and ok_repl2
    if v:
        print("[C12] NEGATIVE: R | e^2 and e*R | Phi are T2-ONLY, permanently")
        print("     point  e=1, R=0, S=1, T=1/6, d0=1, d1=-1/3, d2=0, Phi=%s" % phi_val)
        print("     lies on V(G1,G2,G3,G5):", on_V, " with R=0, e!=0, Phi!=0, d1!=0:", ok_T1)
        print("     => no power of e lies in I+(R) and no power of Phi in I+(e*R):")
        print("        R | e^2 and e*R | Phi have NO divisor certificate on T1.")
        print("     T2 has no such point: G1|_{R=0,d1=0} = 3*e*T and")
        print("        G3|_{R=0,d1=0} = -e^3/2 + 3*S*T, so e != 0 => T = 0 => e = 0 :", ok_t2)
        print("     exact T2 identity   6*R*(d0 + s^2 + d2*s) = -e^2 :", ok_Re2)
        print("     exact T1 replacement  R*(6*d0*e+3*d1*R+6*e*s*(s+d2))")
        print("                            = -e^2*(e + 3*d1*s)  => R | e^2*(e+3*d1*s) :", ok_repl)
        print("     exact T1 replacement  R*(d1*e+d2*R+2*R*s) = e^2*(s^2-d0)")
        print("                            => R | e^2*(s^2-d0) :", ok_repl2)
    return record("C12_negatives", ok)


# =====================================================================
#  C13 -- the T1 spare collapse
# =====================================================================
def spare_collapse(window="sub2", deg_e=10, a=10):
    caps = fsb.STRIP_DEGCAP[window]
    full = sum(caps[n] + 1 for n in ("dm2", "dm3", "dm4"))
    deg_s = caps["dm3"] - deg_e
    after_eS = (caps["dm2"] + 1) + max(deg_s + 1, 0) + 0          # e|S + T-formula
    after_tR = max(caps["dm2"] - a + 1, 0) + max(deg_s + 1, 0) + 0
    after_div = max(caps["dm2"] - a + 1, 0)                        # K-syzygy division
    return full, after_eS, after_tR, after_div


def c13_collapse(v=True):
    rows = []
    for a in (6, 7, 8, 9, 10):
        rows.append((a,) + spare_collapse("sub2", 10, a))
    ok = all(r[4] == 13 - r[0] for r in rows)
    if v:
        print("[C13] the T1 spare collapse on sub2 (deg e = 10, caps 12/14/16)")
        print("      a  |  n | full | after e|S,T | + t^a|R | + K-division")
        for a, f, b1, b2, b3 in rows:
            print("      %-2d | %2d | %3d  |   %3d       |  %3d    | %3d"
                  % (a, 10 - a, f, b1, b2, b3))
        print("      col `+ t^a|R` = 18 - a = n + 8  -- EXACTLY the spine lane's")
        print("        honest `45 -> n+8` collapse (n = 10 - a).  AGREES.")
        print("      col `+ K-division` = 13 - a  -- POLE_THEOREM sec.4's `28 -> 13-a`,")
        print("        which additionally eliminates s by the exact division")
        print("        3*gamma^3*Tm^3*s = W, at the cost of 2k remainder rows.")
    return record("C13_spare_collapse", ok)


# =====================================================================
#  C14 -- the census (READ-ONLY)
# =====================================================================
def census(v=True):
    out = {}
    for win, fn in (("sub2", "phase_d_states_sub2.json"),
                    ("sub1", "phase_d_states_sub1.json")):
        data = json.load(open(os.path.join(HERE, fn), encoding="utf-8"))
        cases = data["cases"]
        ns = lambda L: sum(c["state_count"] for c in L)
        t1 = [c for c in cases if c["branch"] == "T1"]
        t2 = [c for c in cases if c["branch"] == "T2"]
        bad_b = [c for c in cases if any(x > 1 for x in c["b"])]
        bad_b_t1 = [c for c in bad_b if c["branch"] == "T1"]
        # the authoritative post-filter universe, from the frontier lane's own
        # filter (imported READ-ONLY; nothing is written)
        import divisor_filter as dfilt
        filtered = dfilt.filter_universe(data, win)
        alive_t1 = [c for c in filtered["cases"] if c["branch"] == "T1"]
        for c in alive_t1:
            c = c
        alive_t1 = [dict(c, state_count=len(c["states"])) for c in alive_t1]
        amax = max(c["a_t"] for c in cases)
        dege = max(st["deg_e"] for c in cases for st in c["states"])
        out[win] = dict(
            flagcases=len(cases), states=ns(cases),
            t1_flagcases=len(t1), t1_states=ns(t1),
            t2_flagcases=len(t2), t2_states=ns(t2),
            bge2_t1_flagcases=len(bad_b_t1), bge2_t1_states=ns(bad_b_t1),
            alive_t1_flagcases=len(alive_t1), alive_t1_states=ns(alive_t1),
            a_t_max=amax, deg_e_max=dege,
        )
        if v:
            print("[C14] %s  (READ-ONLY from %s)" % (win, fn))
            print("      total            : %4d flag cases  %6d states"
                  % (len(cases), ns(cases)))
            print("      T1               : %4d flag cases  %6d states"
                  % (len(t1), ns(t1)))
            print("      T1 with b_j >= 2 : %4d flag cases  %6d states"
                  "   <- POLE_THEOREM sec.7.1's named prize"
                  % (len(bad_b_t1), ns(bad_b_t1)))
            print("      T1 after divisor_filter.filter_universe (READ-ONLY):"
                  " %4d flag cases  %6d states" % (len(alive_t1), ns(alive_t1)))
            print("        <- t^a|R now applies to EVERY one of these")
            print("      max a_t = %d, max deg_e = %d  (a = a_t <= 10 throughout)"
                  % (amax, dege))
            cols = {}
            for c in alive_t1:
                k = c["a_t"]
                cols.setdefault(k, [0, 0])
                cols[k][0] += 1
                cols[k][1] += c["state_count"]
            print("      surviving T1 by a_t%s:"
                  % (" (and its sub2 spare collapse)" if win == "sub2" else ""))
            for a, (nc, nst) in sorted(cols.items()):
                extra = ("   45 -> %d -> %d" % (18 - a, 13 - a)) if win == "sub2" else ""
                print("         a_t=%-2d  %4d flag cases %6d states%s"
                      % (a, nc, nst, extra))
    ok = out["sub2"]["a_t_max"] <= 10 and out["sub2"]["deg_e_max"] <= 10 \
        and out["sub1"]["a_t_max"] <= 10
    return record("C14_census", ok), out


# =====================================================================
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()
    v = not args.quiet
    if v:
        print("=" * 78)
        print("T1 BRANCH -- t^a | R, the double cover, and what transfers")
        print("=" * 78)
    c1_guard(v)
    c2_ksyzygy(v)
    c2b_phi(v)
    c3_e_divides_S(v)
    c4_reduction(v)
    c5_theorem(v)
    c6_t2_consistency(v)
    c7_escapes(v)
    c8_spine_adjudication(v)
    c8_double_cover(v)
    c10_genus1(v)
    c11_places(v)
    c12_negatives(v)
    c13_collapse(v)
    census(v)

    npass = sum(1 for _, ok, _ in RESULTS if ok)
    if v:
        print("\n" + "=" * 78)
        for tag, ok, msg in RESULTS:
            print("  %-28s %s %s" % (tag, "PASS" if ok else "**FAIL**", msg))
    print("t1_branch: %d/%d checks pass" % (npass, len(RESULTS)))
    return 0 if npass == len(RESULTS) else 1


if __name__ == "__main__":
    raise SystemExit(main())
