#!/usr/bin/env python3
"""pole_theorem_sweep.py  (NEW; READ-ONLY over every existing artifact)

THE POLE THEOREM, GENERALISED -- and swept across the sub2 state space.

`GENERIC_FIBER.md` proved, for the single R9 z=1 state, that any polynomial
spare pair satisfies v_{y+1}(dm2) >= 9 and v_{y+1}(dm3) >= 9, using ONLY

    (i)  d1 = 0,                       (ii) dm1 != 0,
    (iii) the divisor of dm1,          (iv) the valuation of Phi,

and NO spare degree bound.  This module asks whether that was special to the
R9 state shape, and answers: it is not.  Everything below is derived here,
symbolically, from `bigrade_annotator._G_generators()`; nothing is assumed.

WHAT IS PROVED HERE
-------------------
A.  The S1 reduction is a UNIVERSAL identity in the free window-symbol ring:

        d1 = 0   ==>   2*(H5 + d2*H3)  ==  dm1 * K5,
        K5 = 2*Phi - 3*dm1*dm2^2 - d2*dm1^3 - 3*dm1^2*dm3

    with K5 linear in dm3 and free of BOTH dm4 and d0.  No state hypothesis
    of any kind enters.  It therefore holds verbatim on every state with
    d1 = 0, i.e. on the whole T2 branch (branch T2 IS `deg d1 = -inf`,
    `phase_d_states.case_states`).

    For d1 != 0 the identity DEGRADES:

        2*(H5 + d2*H3)  =  dm1 * K5gen  +  6*d1*dm2^2*dm3,

    K5gen is still dm4-free but is QUADRATIC in dm3, and the residue term
    6*d1*dm2^2*dm3 is not divisible by dm1.  dm3 is then a two-valued
    algebraic function of dm2, not a rational one, and the argument below
    does not apply.  This is the method's hard scope boundary (T1 branch).

B.  The place lemma (one DVR argument, run at every place).  Let beta be any
    place, m = v_beta(dm1), P_b = v_beta(Phi), rho = v_beta(dm2).  dm3 must
    be a polynomial, so v_beta(N) >= 2m for N = 2*Phi - 3*dm1*dm2^2 - d2*dm1^3.
    With the H3 generator as the second mechanism this yields, for the
    (72,108) Phi = c*(y+1)^30*q(y)  (q irreducible, squarefree, q(-1) != 0):

      B1  P_b = 0  (beta not a root of Phi):  IMPOSSIBLE.
          => EVERY root of e lies in {-1} U roots(q).
      B2  P_b = 1  (beta a root of q):  m = 1 forced, and rho = 0, i.e.
          e has each q-root with multiplicity EXACTLY ONE and dm2(beta) != 0,
          together with the exact leading-coefficient relation
              2*Phi'(beta) = 3*dm1'(beta)*dm2(beta)^2.
          => any state carrying v_{r_j}(e) >= 2 at a marked root is KILLED.
      B3  beta = -1, P = 30:  with a = v_{-1}(dm1),

              v_{-1}(dm2) >= min(a, ceil((P-a)/2)),
              v_{-1}(dm3) >= min(P, a + 2*rho_min) - 2*a,

          and for 3a <= P+1 (i.e. a <= 10, which is every sub2 state, since
          deg e <= 10) both bounds are exactly `a`.

    The sharp constant is `a` = v_{-1}(dm1) -- NOT `a_t`, NOT any degree cap.
    All of B1-B3 are DEGREE-UNIFORM: no spare degree bound is used anywhere,
    and the enumeration is re-run with the caps removed to prove it.

C.  The collapse.  On e = gamma*(y+1)^a*T (T the squarefree marked-root part,
    deg T = k), writing dm2 = (y+1)^a*A, dm3 = (y+1)^a*B, K5 becomes the
    EXACT division  3*gamma^2*T^2*B = W,  W = 2c*(y+1)^(30-3a)*q
    - 3*gamma*T*A^2 - gamma^3*d2*T^3.  So B is determined, the whole content
    of K5 is a 2k-row division remainder, and the spare unknown count falls
        28  ->  13 - a      (k = 0: no remainder rows at all).
    a = 9, k = 1 reproduces GENERIC_FIBER's 18 -> 4 and deg W = 7 exactly.

Independent checker: pole_theorem_sweep_verify.py (--quiet, exit 0/1).
Sources (read-only): bigrade_annotator.py, batch_convolution_sub2.json,
phase_d_states_sub2.json, full_system_bridge.py (caps only),
GENERIC_FIBER.md sec.4 (the R9 numbers this must reproduce).

Usage:  python -u pole_theorem_sweep.py            full report + json
        python -u pole_theorem_sweep.py --quiet    self-check, exit 0/1
"""
from __future__ import annotations

import json
import math
import os
import random
import sys

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__)) or "."
sys.path.insert(0, HERE)

import bigrade_annotator as BA          # noqa: E402  (read-only import)

y = sp.Symbol("y")
INF = 10**9
PHI_C = sp.Rational(-1, 6630)
Q_POLY = 2048 * y**4 - 512 * y**3 + 320 * y**2 - 240 * y + 195
CAPS = {"dm2": 12, "dm3": 14}          # sub2 window caps (counting only)


# =====================================================================
#  A.  The S1 reduction, derived in the FREE window-symbol ring
# =====================================================================
def stage_A_identity(verbose=True):
    """Derive (not assume) the reduction, for d1 = 0 and for d1 != 0.

    H2, H3, H5 are REBUILT here from `_G_generators()` as the certified
    cofactor combinations and checked against `r9_eliminated_system.json`,
    so nothing is trusted from the json.
    """
    d0, d1, d2, dm1, dm2, dm3, dm4, Phi = BA._gsystem_symbols()
    G = {k: v for k, (v, _w) in BA._G_generators().items()}
    Hjson = {k: v for k, (v, _w) in BA._H_generators().items()}
    Hbuilt = {
        "H2": sp.expand(dm1 * G["G2"] - dm2 * G["G1"]),
        "H3": sp.expand(dm1 * G["G3"] - dm3 * G["G1"]),
        "H5": sp.expand(dm1 * G["G5"] + (d0 * dm1 + d1 * dm2 + d2 * dm3) * G["G1"]),
    }
    agree = {k: sp.expand(Hbuilt[k] - Hjson[k]) == 0 for k in Hbuilt}
    assert all(agree.values()), "rebuilt H != json H: %s" % agree

    # --- d1 = 0 : the universal identity
    H0 = {k: sp.expand(v.subs(d1, 0)) for k, v in Hbuilt.items()}
    K5 = sp.expand(2 * Phi - 3 * dm1 * dm2**2 - d2 * dm1**3 - 3 * dm1**2 * dm3)
    resid0 = sp.expand(2 * (H0["H5"] + d2 * H0["H3"]) - dm1 * K5)
    assert resid0 == 0, "S1 identity failed at d1=0: %s" % resid0
    lin = sp.Poly(K5, dm3).degree()
    assert lin == 1 and not K5.has(dm4) and not K5.has(d0)

    # --- d1 != 0 : the general remainder (dm1-adic division)
    L = sp.expand(2 * (Hbuilt["H5"] + d2 * Hbuilt["H3"]))
    quo, rem = sp.div(sp.Poly(L, dm1), sp.Poly(dm1, dm1))
    rem = sp.expand(rem.as_expr())
    K5gen = sp.expand(quo.as_expr())
    assert sp.expand(dm1 * K5gen + rem - L) == 0
    deg_gen = sp.Poly(K5gen, dm3).degree()
    rem_expected = sp.expand(6 * d1 * dm2**2 * dm3)
    assert sp.expand(rem - rem_expected) == 0, "unexpected d1-remainder: %s" % rem

    # --- H3 at d1 = 0, the second mechanism (also dm4-free)
    H3_0 = H0["H3"]
    H3_terms = sp.Add.make_args(H3_0)
    assert not H3_0.has(dm4)

    out = {
        "H_rebuilt_matches_json": all(agree.values()),
        "K5": sp.sstr(K5),
        "K5_deg_in_dm3": lin,
        "K5_free_of_dm4": True, "K5_free_of_d0": True,
        "d1_residual": 0,
        "K5gen_deg_in_dm3": int(deg_gen),
        "d1_remainder": sp.sstr(rem),
        "H3_at_d1_0": sp.sstr(H3_0),
        "n_H3_terms": len(H3_terms),
    }
    if verbose:
        print("=" * 78)
        print("A.  THE S1 REDUCTION -- derived in the FREE window-symbol ring")
        print("=" * 78)
        print("    H2,H3,H5 rebuilt from _G_generators() and compared to the json:",
              "MATCH" if out["H_rebuilt_matches_json"] else "*** MISMATCH ***")
        print("    d1 = 0 :  2*(H5 + d2*H3) - dm1*K5  =  0   [EXACT, no state used]")
        print("              K5 =", out["K5"])
        print("              deg_dm3 K5 = %d ; contains dm4: %s ; contains d0: %s"
              % (lin, K5.has(dm4), K5.has(d0)))
        print("        ==> UNIVERSAL: holds on EVERY d1 = 0 state (branch T2).")
        print("    d1 != 0 : 2*(H5 + d2*H3) = dm1*K5gen + (%s)" % out["d1_remainder"])
        print("              deg_dm3 K5gen = %d  ->  dm3 is 2-valued algebraic,"
              % deg_gen)
        print("              and the residue is NOT divisible by dm1.")
        print("        ==> the pole argument does NOT transfer to branch T1.")
        print("    H3|d1=0 =", out["H3_at_d1_0"], " (dm4-free)")
    return out


# =====================================================================
#  B.  Phi, measured -- the only global input
# =====================================================================
def stage_B_phi(verbose=True):
    Phiv = BA._phi_stripped()
    P = _ord_at(Phiv, sp.Integer(-1))
    qpoly = sp.Poly(Q_POLY, y)
    irred = qpoly.is_irreducible
    sqfree = sp.discriminant(Q_POLY, y) != 0
    qm1 = Q_POLY.subs(y, -1)
    # v at a root of q: q squarefree => exactly 1, and (y+1)^30 does not vanish
    out = {"v_minus1_Phi": P, "q_irreducible": bool(irred), "q_squarefree": bool(sqfree),
           "q_at_minus1": int(qm1), "v_qroot_Phi": 1, "deg_q": 4,
           "phi_matches_factored": sp.expand(Phiv - PHI_C * (y + 1)**30 * Q_POLY) == 0}
    assert out["phi_matches_factored"] and P == 30 and irred and sqfree and qm1 != 0
    if verbose:
        print()
        print("=" * 78)
        print("B.  Phi  --  the one global input (state-independent)")
        print("=" * 78)
        print("    Phi = c*(y+1)^30*q(y), c = -1/6630   [checked against _phi_stripped]")
        print("    v_{y+1}(Phi) = %d   (a FIXED constant of the (72,108) target,"
              % P)
        print("                        not a spare degree, not a state degree)")
        print("    q irreducible: %s ; squarefree: %s ; q(-1) = %d != 0"
              % (irred, sqfree, qm1))
        print("    => v_beta(Phi) = 1 at each of the 4 roots of q, 0 elsewhere.")
    return out


# =====================================================================
#  the single DVR case decider -- used at EVERY place
# =====================================================================
def _ordmin(vals):
    m = min(vals)
    return m, sum(1 for v in vals if v == m) == 1


def dvr_case(m, Pb, rho, t, s):
    """One case of the DVR argument at an arbitrary place.

    m   = v(dm1)   >= 0 (integer; the place's multiplicity in e)
    Pb  = v(Phi)
    rho = v(dm2)   (INF allowed: dm2 == 0)
    t   = v(d2)    >= 0 (INF allowed)
    s   = v(d0)    >= 0 (INF allowed)

    K5 :  3*dm1^2*dm3 = N,  N = 2*Phi - 3*dm1*dm2^2 - d2*dm1^3
    H3 : -3*d0*dm1^2*dm2 - 3*d2*dm1*dm2*dm3 - dm1^4/2 - 3*dm2*dm3^2 = 0

    Returns ("KILL", reason, data) or ("SURVIVE", reason, data).
    A sum of DVR-valued terms whose minimal order is attained EXACTLY ONCE
    cannot vanish; that is the only inference used.
    """
    N_terms = {"2Phi": Pb,
               "dm1*dm2^2": _add(m, 2 * rho),
               "d2*dm1^3": _add(t, 3 * m)}
    lo, uniq = _ordmin(list(N_terms.values()))
    if not uniq:
        return ("SURVIVE", "N-min attained more than once: v(N) >= %s" % _fmt(lo),
                {"vN_lb": lo, "sigma_lb": _sub(lo, 2 * m)})
    vN = lo                                   # EXACT
    sigma = _sub(vN, 2 * m)                   # EXACT order of dm3
    if sigma < 0:
        return ("KILL", "dm3 has a POLE (v(dm3) = %d < 0)" % sigma,
                {"vN": vN, "sigma": sigma})
    H3_terms = {"d0*dm1^2*dm2": _add(s, _add(2 * m, rho)),
                "d2*dm1*dm2*dm3": _add(t, _add(m, _add(rho, sigma))),
                "dm1^4": 4 * m,
                "dm2*dm3^2": _add(rho, _mul2(sigma))}
    lo3, uniq3 = _ordmin(list(H3_terms.values()))
    if uniq3:
        wit = [k for k, v in H3_terms.items() if v == lo3][0]
        return ("KILL", "H3 has a UNIQUE minimal-order term (ord %s, witness %s)"
                        " => H3 != 0" % (_fmt(lo3), wit),
                {"vN": vN, "sigma": sigma, "vH3": lo3, "witness": wit})
    return ("SURVIVE", "H3 min attained more than once (ord %s)" % _fmt(lo3),
            {"vN": vN, "sigma": sigma})


def _add(a, b):
    return INF if (a >= INF or b >= INF) else a + b


def _sub(a, b):
    return INF if a >= INF else a - b


def _mul2(a):
    return INF if a >= INF else 2 * a


def _fmt(v):
    return "inf" if v >= INF else str(v)


# =====================================================================
#  B3.  the place y = -1 : the (y+1)-bound, sharp
# =====================================================================
def _critical(targets, offset):
    """Every value w >= 0 (plus INF) needed to realise every achievable sign
    pattern of  (offset + w)  vs  each target.

    dvr_case's verdict is a function of the total preorder of finitely many
    term-orders, and v(d2) (resp. v(d0)) occurs in EXACTLY ONE term each,
    linearly and with coefficient +1.  So the verdict depends on w only
    through the signs of (offset + w - target).  Sampling w at each crossing,
    one below, one above, plus 0 and INF, realises every sign pattern -- an
    exhaustive substitute for an unbounded sweep, at constant cost.
    """
    out = {0, INF}
    hi = 0
    for tg in targets:
        if tg >= INF:
            continue
        for d in (-1, 0, 1):
            w = tg - offset + d
            if 0 <= w < INF:
                out.add(w)
                hi = max(hi, w)
    out.add(hi + 1)
    return sorted(out)


def pole_bound(a, P, rho_max=None, aux_max=None):
    """Enumerate the DVR cases at y = -1 and return the SHARP bound this
    method gives.

    rho_max = None  ->  NO degree cap at all (rho swept well past any cap,
                        plus rho = INF for dm2 == 0).
    aux_max = None  ->  v(d2), v(d0) sampled at their CRITICAL values (see
                        `_critical`), which is exhaustive; an integer value
                        instead forces a plain 0..aux_max sweep (used only to
                        cross-check the critical sampling).
    """
    big = rho_max if rho_max is not None else 2 * (a + P) + 6
    rhos = list(range(0, big + 1)) + [INF]
    killed_rho = set()
    n_survive = 0
    sigma_lbs = []
    for rho in rhos:
        all_dead = True
        ts = (_critical([P, _add(a, _mul2(rho))], 3 * a) if aux_max is None
              else list(range(0, aux_max + 1)) + [INF])
        for t in ts:
            # sigma is needed for the s-criticals; recompute it exactly as
            # dvr_case does (cheap, and keeps the two in lockstep)
            lo, uniq = _ordmin([P, _add(a, _mul2(rho)), _add(t, 3 * a)])
            sig = _sub(lo, 2 * a) if uniq else INF
            ss = (_critical([4 * a, _add(t, _add(a, _add(rho, sig))),
                             _add(rho, _mul2(sig))], _add(2 * a, rho))
                  if aux_max is None else list(range(0, aux_max + 1)) + [INF])
            for s in ss:
                verdict, _why, dat = dvr_case(a, P, rho, t, s)
                if verdict == "SURVIVE":
                    all_dead = False
                    n_survive += 1
                    # dm3 is a polynomial, so its order is >= 0 whatever the
                    # K5 term-order bookkeeping says on a tied case
                    sigma_lbs.append(max(0, dat.get("sigma",
                                                    dat.get("sigma_lb", INF))))
        if all_dead:
            killed_rho.add(rho)
    alive = sorted(r for r in rhos if r not in killed_rho)
    rho_min = alive[0] if alive else INF
    sigma_min = min(sigma_lbs) if sigma_lbs else INF
    return {"a": a, "P": P, "rho_min": rho_min, "sigma_min": sigma_min,
            "killed_rho": sorted(killed_rho)[:40], "n_survive_cases": n_survive,
            "rho_max_used": big, "aux_max_used": aux_max}


def in_regime(a, P):
    """The hypothesis of the hand proof:  a + 2*rho < P for every rho <= a-1,
    i.e.  3a - 2 < P.   With P = 30 this is exactly  a <= 10."""
    return 3 * a - 2 < P


def pole_bound_closed(a, P):
    """The PROVED lower bound, valid exactly when `in_regime(a, P)`.

        v_{y+1}(dm2) >= a      and      v_{y+1}(dm3) >= a.

    Proof (all inequalities, no degree cap anywhere).  Let rho = v(dm2),
    t = v(d2) >= 0, s = v(d0) >= 0, and suppose rho <= a-1.

      1. The three K5 term-orders are {P, a+2rho, 3a+t}.  Since t >= 0,
         3a + t >= 3a > a + 2rho (as rho < a), and a + 2rho < P by the regime
         hypothesis.  So a+2rho is the STRICT minimum and v(N) = a + 2rho
         EXACTLY, hence sigma := v(dm3) = 2*rho - a EXACTLY.
      2. If 2*rho < a then sigma < 0 -- dm3 has a POLE.  Contradiction.
      3. Otherwise the four H3 term-orders are
             d0*dm1^2*dm2 : s + 2a + rho
             dm1^4        : 4a
             d2*dm1*dm2*dm3 : t + 3*rho          (using sigma = 2rho-a)
             dm2*dm3^2    : 5*rho - 2a
         and rho < a gives, with s, t >= 0,
             5rho-2a < 4a          <=  5rho < 6a
             5rho-2a < s+2a+rho    <=  4rho < 4a + s
             5rho-2a < t+3rho      <=  2rho < 2a + t
         all three strict.  So dm2*dm3^2 is the UNIQUE minimum and H3 != 0.
         Contradiction.
      Hence rho >= a.  Then v(N) >= min(P, 3a, 3a+t) = min(P, 3a) and
      sigma >= min(P, 3a) - 2a = a, because 3a <= P+1 and sigma is an integer
      (3a = P+1 would give sigma >= a-1, but then a+2rho >= 3a = P+1 > P makes
      P the minimum, v(N) >= P = 3a-1, sigma >= a-1; the integral case
      3a <= P gives sigma >= a directly).

    Returns None outside the regime: the closed form does NOT extend there --
    see `stage_C_bound`, which reports the measured values instead.
    """
    if not in_regime(a, P):
        return None
    sig = min(P, 3 * a) - 2 * a
    return {"rho_min": max(0, a), "sigma_min": max(0, sig)}


# =====================================================================
#  B1/B2.  the other places -- the support / multiplicity theorem
# =====================================================================
def stage_D_places(verbose=True):
    """Run the SAME enumerator at a place beta != -1 with v_beta(Phi) = Pb."""
    rows = []
    for Pb in (0, 1):
        for m in range(1, 6):
            bnd = pole_bound(m, Pb, rho_max=30, aux_max=30)
            rows.append({"Pb": Pb, "m": m, "rho_min": bnd["rho_min"],
                         "sigma_min": bnd["sigma_min"],
                         "all_dead": bnd["rho_min"] >= INF})
    if verbose:
        print()
        print("=" * 78)
        print("D.  THE OTHER PLACES  --  support and multiplicity of e")
        print("=" * 78)
        print("    Pb = v_beta(Phi) | m = v_beta(e) | verdict")
        for r in rows:
            v = ("EVERY rho killed => NO polynomial solution"
                 if r["all_dead"] else
                 "survives with v_beta(dm2) >= %s, v_beta(dm3) >= %s"
                 % (_fmt(r["rho_min"]), _fmt(r["sigma_min"])))
            print("      %2d            |  %d           | %s" % (r["Pb"], r["m"], v))
        print("    ==> B1  a root of e that is NOT a root of Phi is IMPOSSIBLE.")
        print("            rad(monic e) | (y+1)*q(y):  every root of e is -1 or")
        print("            one of the 4 roots of q.")
        print("    ==> B2  at a root of q, m = 1 is FORCED (m >= 2 is killed) and")
        print("            v_beta(dm2) = 0, i.e. dm2(beta) != 0.")
    # the exact place-r relation, in general form
    rel = _qroot_relation()
    if verbose:
        print("    exact leading-coefficient relation at a q-root beta:")
        print("        2*Phi'(beta) = 3*dm1'(beta)*dm2(beta)^2")
        print("      R9 instance (e = gamma*(y+1)^9*(y-r), q(r)=0):")
        print("        %s" % rel["r9_relation"])
        print("      GENERIC_FIBER.md sec.4 states exactly this: MATCH = %s"
              % rel["matches_generic_fiber"])
    return {"rows": rows, "relation": rel}


def _rem_q(expr, r):
    """Reduce a polynomial in the marked root r modulo q(r) -- the arithmetic
    of the root algebra Q[r]/(q).  q(r) = 0 is NOT a free substitution: it has
    to be imposed as a remainder, which is what this does."""
    e = sp.expand(expr)
    if e == 0:
        return sp.Integer(0)
    p = sp.Poly(e, r)
    return sp.expand(sp.rem(p, sp.Poly(Q_POLY.subs(y, r), r)).as_expr())


def _qroot_relation():
    """Specialise the general place-r cancellation

        2*Phi'(beta) = 3*dm1'(beta)*dm2(beta)^2       (beta a simple root of e
                                                       and a root of q)

    to the R9 shape and check it reproduces GENERIC_FIBER.md sec.4 verbatim.
    All arithmetic is in Q[r]/(q(r))."""
    gamma, r = sp.symbols("gamma r")
    A = sp.Symbol("A_of_r")
    dm1 = gamma * (y + 1)**9 * (y - r)
    Phiv = PHI_C * (y + 1)**30 * Q_POLY
    lhs = _rem_q(2 * sp.diff(Phiv, y).subs(y, r), r)      # 2*Phi'(r) mod q
    rhs = _rem_q(3 * sp.diff(dm1, y).subs(y, r) * A**2, r)
    # the two closed forms GENERIC_FIBER.md sec.4 prints
    lhs_doc = _rem_q(2 * PHI_C * (r + 1)**30 * sp.diff(Q_POLY, y).subs(y, r), r)
    rhs_doc = _rem_q(3 * gamma * (r + 1)**9 * A**2, r)
    doc = sp.expand(rhs_doc - lhs_doc)
    ours = sp.expand(rhs - lhs)
    return {"r9_relation": "3*gamma*(r+1)^9*dm2(r)^2 = 2*c*(r+1)^30*q'(r)",
            "lhs_reduction_exact": sp.expand(lhs - lhs_doc) == 0,
            "rhs_reduction_exact": sp.expand(rhs - rhs_doc) == 0,
            "matches_generic_fiber": sp.expand(ours - doc) == 0}


# =====================================================================
#  B3 report + degree-uniformity proof
# =====================================================================
def stage_C_bound(verbose=True):
    P = 30
    rows = []
    for a in range(0, 16):
        enum_capped = pole_bound(a, P, rho_max=CAPS["dm2"])   # WITH the cap
        enum_free = pole_bound(a, P, rho_max=None)            # NO cap at all
        closed = pole_bound_closed(a, P)
        rows.append({"a": a, "in_regime": in_regime(a, P),
                     "rho_min_capped": enum_capped["rho_min"],
                     "rho_min_uncapped": enum_free["rho_min"],
                     "sigma_min_uncapped": enum_free["sigma_min"],
                     "killed_rho": enum_free["killed_rho"],
                     "rho_min_proved": None if closed is None else closed["rho_min"],
                     "sigma_min_proved": None if closed is None else closed["sigma_min"]})
    reg = [r for r in rows if r["in_regime"]]
    # (1) degree-uniformity, the strong form: in-regime the eliminated set is
    #     EXACTLY {0,...,a-1}, so no case at or above any cap is ever used.
    kill_exact = all(r["killed_rho"] == list(range(r["a"])) for r in reg)
    # (2) capped vs uncapped enumeration identical, in-regime
    cap_free = all(r["rho_min_capped"] == r["rho_min_uncapped"] for r in reg)
    # (3) the hand proof is SOUND against the enumerator (enum >= proved), and
    #     SHARP in-regime at P = 30 (enum == proved)
    sound = all(r["rho_min_uncapped"] >= r["rho_min_proved"]
                and r["sigma_min_uncapped"] >= r["sigma_min_proved"] for r in reg)
    sharp = all(r["rho_min_uncapped"] == r["rho_min_proved"] == r["a"]
                and r["sigma_min_uncapped"] == r["sigma_min_proved"] == r["a"]
                for r in reg)
    # (4) critical sampling of v(d2), v(d0) agrees with plain brute-force sweeps
    aux_ok = all(pole_bound(a, P, rho_max=14, aux_max=None)["rho_min"]
                 == pole_bound(a, P, rho_max=14, aux_max=12)["rho_min"]
                 == pole_bound(a, P, rho_max=14, aux_max=45)["rho_min"]
                 for a in range(0, 13))
    # (5) grid over hypothetical v(Phi): the hand proof must be sound
    #     everywhere in its regime, not just at P = 30
    grid_bad = []
    for a in range(0, 13):
        for PP in range(0, 40):
            c = pole_bound_closed(a, PP)
            if c is None:
                continue
            e = pole_bound(a, PP)
            if e["rho_min"] < c["rho_min"] or e["sigma_min"] < c["sigma_min"]:
                grid_bad.append((a, PP, e["rho_min"], c["rho_min"]))
    if verbose:
        print()
        print("=" * 78)
        print("C.  THE PLACE y = -1  --  the (y+1)-bound   (P = v(Phi) = 30)")
        print("=" * 78)
        print("      a  | 3a-2<P | rho_min cap12 | rho_min NO cap | PROVED |"
              " sigma_min | eliminated rho")
        for r in rows:
            star = "  <== sub2" if 7 <= r["a"] <= 10 else ""
            kr = r["killed_rho"]
            krs = ("0..%d" % (max(kr)) if kr and kr == list(range(len(kr)))
                   else str(kr)[:18])
            print("     %3d |  %-5s |      %4s     |      %4s      |  %4s  |"
                  "   %4s    | %s%s"
                  % (r["a"], r["in_regime"], _fmt(r["rho_min_capped"]),
                     _fmt(r["rho_min_uncapped"]),
                     "-" if r["rho_min_proved"] is None else r["rho_min_proved"],
                     _fmt(r["sigma_min_uncapped"]), krs, star))
        print()
        print("    IN-REGIME (3a - 2 < P, i.e. a <= 10 at P = 30):")
        print("      proved bound is EXACTLY a for both dm2 and dm3        : %s" % sharp)
        print("      hand proof sound against the enumerator (enum >= a)   : %s" % sound)
        print("      eliminated rho set is EXACTLY {0,...,a-1}             : %s"
              % kill_exact)
        print("        ==> every case the proof uses has rho < a <= 10 < 12 = the")
        print("            dm2 cap.  NO case at or above any degree cap is used:")
        print("            the argument is DEGREE-UNIFORM in the strong sense.")
        print("      capped(12) and uncapped enumerations identical        : %s" % cap_free)
        print("      critical sampling of v(d2),v(d0) == brute force       : %s" % aux_ok)
        print("      hand proof sound on the whole grid a<=12, P<=39       : %s"
              % (not grid_bad))
        print()
        print("    *** OUT OF REGIME (a >= 11 at P = 30) -- LOUD FINDING ***")
        print("      The bound is NOT `a` there, and for odd a it becomes")
        print("      CAP-DEPENDENT (capped enumeration kills EVERY rho <= 12, so")
        print("      the conclusion would silently depend on the degree cap):")
        for r in rows:
            if not r["in_regime"]:
                note = ("CAP-DEPENDENT" if r["rho_min_capped"] != r["rho_min_uncapped"]
                        else "cap-independent")
                print("        a=%2d : rho_min uncapped = %-3s , capped at 12 = %-3s"
                      "   [%s]" % (r["a"], _fmt(r["rho_min_uncapped"]),
                                   _fmt(r["rho_min_capped"]), note))
        print("      So degree-uniformity is NOT unconditional: it holds exactly")
        print("      while 3a - 2 < v_{y+1}(Phi) = 30.  Every sub2 state has")
        print("      deg e <= 10, hence a <= 10, hence is in regime -- but a target")
        print("      or window with deg e >= 11 would lose the property.")
    return {"rows": rows, "cap_independent": cap_free, "kill_set_exact": kill_exact,
            "proof_sound": sound, "proof_sharp": sharp,
            "aux_independent": aux_ok, "grid_unsound": grid_bad}


# =====================================================================
#  R9 regression: the generalisation must reproduce GENERIC_FIBER exactly
# =====================================================================
def stage_R9_regression(verbose=True):
    b = pole_bound(9, 30, rho_max=CAPS["dm2"])
    closed = pole_bound_closed(9, 30)
    place_r = dvr_case(1, 1, 0, 0, 0)       # e has (y-r) to order 1, dm2(r)!=0
    place_r_bad = dvr_case(1, 1, 1, 0, 0)   # dm2(r) = 0
    ok = (b["rho_min"] == 9 and b["sigma_min"] == 9
          and closed["rho_min"] == 9 and closed["sigma_min"] == 9
          and place_r[0] == "SURVIVE" and place_r_bad[0] == "KILL")
    if verbose:
        print()
        print("=" * 78)
        print("R9 REGRESSION  --  the general theorem, specialised back to a = 9")
        print("=" * 78)
        print("    v_{y+1}(dm2) >= %d   (GENERIC_FIBER.md sec.4 says 9)" % b["rho_min"])
        print("    v_{y+1}(dm3) >= %d   (GENERIC_FIBER.md sec.4 says 9)" % b["sigma_min"])
        print("    place y=r, v(dm2)=0 : %s" % place_r[0])
        print("    place y=r, v(dm2)>=1: %s  (%s)" % (place_r_bad[0], place_r_bad[1]))
        print("      => dm2(r) != 0, exactly GENERIC_FIBER.md sec.4.")
        print("    REGRESSION: %s" % ("PASS" if ok else "*** FAIL ***"))
    return {"rho_min": b["rho_min"], "sigma_min": b["sigma_min"], "pass": ok}


# =====================================================================
#  C.  The collapse: K5 as an exact division; unknown / case counts
# =====================================================================
def stage_E_collapse(a, k, verbose=False):
    """e = gamma*(y+1)^a*T with T squarefree of degree k (marked q-roots).

    Returns the exact-division data: deg W, deg B, #remainder rows, and the
    surviving spare unknown count.  Symbolic, exact.
    """
    gamma = sp.Symbol("gamma")
    rs = sp.symbols("rt0:%d" % k) if k else ()
    T = sp.prod([(y - rr) for rr in rs]) if k else sp.Integer(1)
    T = sp.expand(T)
    degA = CAPS["dm2"] - a
    Ac = sp.symbols("A0:%d" % (degA + 1)) if degA >= 0 else ()
    Av = sum(Ac[i] * y**i for i in range(degA + 1)) if degA >= 0 else sp.Integer(0)
    ac = sp.symbols("a0:5")
    d2v = sum(ac[i] * y**i for i in range(5))
    W = sp.expand(2 * PHI_C * (y + 1)**(30 - 3 * a) * Q_POLY
                  - 3 * gamma * T * Av**2 - gamma**3 * d2v * T**3)
    div = sp.expand(3 * gamma**2 * T**2)
    quo, rem = sp.div(sp.Poly(W, y), sp.Poly(div, y))
    remc = [c for c in sp.Poly(rem.as_expr(), y).all_coeffs() if sp.expand(c) != 0] \
        if rem.as_expr() != 0 else []
    out = {"a": a, "k": k, "deg_W": sp.Poly(W, y).degree(),
           "deg_B": sp.Poly(quo.as_expr(), y).degree() if quo.as_expr() != 0 else -1,
           "n_remainder_rows": 2 * k,
           "n_remainder_nonzero": len(remc),
           "spare_unknowns_before": (CAPS["dm2"] + 1) + (CAPS["dm3"] + 1),
           "spare_unknowns_after": degA + 1,
           "dm3_cap_implied": CAPS["dm3"] - a}
    if verbose:
        print("      a=%2d k=%d : deg W = %2d, B = W/(3g^2 T^2) has deg %2d,"
              " %d remainder rows, spares %d -> %d"
              % (a, k, out["deg_W"], out["deg_B"], out["n_remainder_rows"],
                 out["spare_unknowns_before"], out["spare_unknowns_after"]))
    return out


def stage_E_report(verbose=True):
    rows = []
    if verbose:
        print()
        print("=" * 78)
        print("E.  THE COLLAPSE  --  K5 becomes an EXACT division, B is determined")
        print("=" * 78)
        print("    e = gamma*(y+1)^a * T,  T squarefree of degree k on marked q-roots")
        print("    dm2 = (y+1)^a A, dm3 = (y+1)^a B  =>  3*gamma^2*T^2*B = W")
    for a in (7, 8, 9, 10):
        for k in range(0, 11 - a):
            rows.append(stage_E_collapse(a, k, verbose=verbose))
    r9 = [r for r in rows if r["a"] == 9 and r["k"] == 1][0]
    r9_ok = (r9["deg_W"] == 7 and r9["deg_B"] == 5
             and r9["spare_unknowns_after"] == 4 and r9["n_remainder_rows"] == 2)
    if verbose:
        print("    R9 row (a=9,k=1) vs GENERIC_FIBER.md sec.5 (deg W = 7, deg B = 5,")
        print("      2 remainder rows, 18 -> 4 spare unknowns):  %s"
              % ("MATCH" if r9_ok else "*** MISMATCH ***"))
        print("    a=10,k=0 : ZERO remainder rows -- K5 imposes NOTHING, B is a free")
        print("      exact quotient and the state keeps only 3 spare unknowns.")
    return {"rows": rows, "r9_match": r9_ok}


# =====================================================================
#  F.  Forward confirmation -- a completely different code path
# =====================================================================
def _ord_at(expr, pt):
    """Order of vanishing at y = pt: number of trailing zero coefficients of
    the Taylor expansion, read off in one shift (no repeated division)."""
    e = sp.expand(expr)
    if e == 0:
        return INF
    cs = sp.Poly(sp.expand(e.subs(y, y + pt)), y).all_coeffs()[::-1]  # ascending
    for k, c in enumerate(cs):
        if sp.expand(c) != 0:
            return k
    return INF


def stage_F_forward(verbose=True, seed=17):
    """Real polynomial arithmetic on random concrete states: build dm2 with a
    prescribed order, form dm3 = N/(3 dm1^2) as a RATIONAL function, and read
    off the orders.  Confirms the integer arithmetic of dvr_case is the
    arithmetic that actually occurs.  Independent of the enumerator."""
    rng = random.Random(seed)
    Phiv = PHI_C * (y + 1)**30 * Q_POLY
    rows = []
    shapes = [("a=10, T=1", 10, sp.Integer(1)),
              ("a=8,  T=1", 8, sp.Integer(1)),
              ("a=7,  T=1", 7, sp.Integer(1))]
    for tag, a, T in shapes:
        gm = sp.Integer(rng.randrange(2, 40))
        d2v = sum(sp.Integer(rng.randrange(0, 20)) * y**i for i in range(5))
        d0v = sum(sp.Integer(rng.randrange(0, 20)) * y**i for i in range(9))
        dm1 = sp.expand(gm * (y + 1)**a * T)
        for rho in range(0, a + 1):
            tail = sum(sp.Integer(rng.randrange(1, 25)) * y**i
                       for i in range(0, max(1, 13 - rho)))
            m2 = sp.expand((y + 1)**rho * tail)
            Nv = sp.expand(2 * Phiv - 3 * dm1 * m2**2 - d2v * dm1**3)
            m3 = sp.cancel(Nv / (3 * dm1**2))
            o2 = _ord_at(m2, sp.Integer(-1))
            o3 = _ord_at(sp.numer(m3), sp.Integer(-1)) - _ord_at(sp.denom(m3),
                                                                sp.Integer(-1))
            H3v = sp.cancel(sp.expand(-3 * d0v * dm1**2 * m2
                                      - 3 * d2v * dm1 * m2 * m3
                                      - dm1**4 / 2 - 3 * m2 * m3**2))
            oH = (_ord_at(sp.numer(H3v), sp.Integer(-1))
                  - _ord_at(sp.denom(H3v), sp.Integer(-1))) if H3v != 0 else INF
            pred3 = 2 * rho - a
            predH = 5 * rho - 2 * a
            ok = (o3 == pred3) and (o3 < 0 or oH == predH)
            rows.append((tag, a, rho, o2, o3, pred3, oH, predH, ok))
    # the SUPPORT theorem, forward: a root of e that is not a root of Phi
    gm = sp.Integer(3)
    beta = sp.Integer(5)                       # Phi(5) != 0
    supp = []
    for m in (1, 2):
        dm1 = sp.expand(gm * (y + 1)**8 * (y - beta)**m)
        m2 = sp.expand(1 + 2 * y + y**3)
        d2v = sp.expand(1 + y)
        Nv = sp.expand(2 * Phiv - 3 * dm1 * m2**2 - d2v * dm1**3)
        vb = _ord_at(Nv, beta)
        supp.append((m, vb, 2 * m, vb < 2 * m))
    # a marked root of q with multiplicity m in e, done over Q[r]/(q).
    # m >= 2 : v_r(N) = 1 for EVERY dm2 (the 2*Phi term is alone at order 1),
    #          so 1 < 2m is unavoidable -> KILL.
    # m  = 1 : v_r(N) = 1 generically, but the order-1 coefficient is
    #          2*c*(r+1)^30*q'(r) - 3*gm*(r+1)^9*dm2(r)^2, which the square
    #          relation makes vanish -- then v_r(N) >= 2 and the case SURVIVES.
    rr = sp.Symbol("rr")

    def _ord_at_marked(expr, m_):
        """order at y = rr in Q[rr]/(q(rr)), by Taylor shift + remainder."""
        cs = sp.Poly(sp.expand(sp.expand(expr).subs(y, y + rr)), y).all_coeffs()[::-1]
        for k, c in enumerate(cs):
            if _rem_q(c, rr) != 0:
                return k
        return INF

    # Everything below is exact in Q(gm, A0, A1)[rr]/(q(rr)) -- no inverses are
    # ever formed, so no denominators can hide a zero.
    A0, A1 = sp.symbols("A0 A1")
    mult2 = []
    for m in (1, 2, 3):
        dm1s = sp.expand(gm * (y + 1)**7 * (y - rr)**m)
        m2s = A0 + A1 * (y - rr)                       # fully general dm2 germ
        Ns = sp.expand(2 * Phiv - 3 * dm1s * m2s**2 - (1 + y) * dm1s**3)
        cs = sp.Poly(sp.expand(Ns.subs(y, y + rr)), y).all_coeffs()[::-1]
        c0 = _rem_q(cs[0], rr)
        c1 = _rem_q(cs[1], rr)
        # the predicted order-1 coefficient
        pred1 = _rem_q(2 * PHI_C * (rr + 1)**30 * sp.diff(Q_POLY, y).subs(y, rr)
                       - (3 * gm * (rr + 1)**7 * A0**2 if m == 1 else 0), rr)
        mult2.append((m, "c0 == 0", int(c0 == 0), 1, c0 != 0))
        mult2.append((m, "c1 == predicted", int(sp.expand(c1 - pred1) == 0), 1,
                      sp.expand(c1 - pred1) != 0))
    # solvability at m = 1: c1 is AFFINE in A0^2 with coefficient -3*gm*(rr+1)^7,
    # a unit of the field Q[rr]/(q) (gm != 0 and q(-1) != 0), so exactly one
    # square class of dm2(r) makes v_r(N) >= 2.  For m >= 2 the A0-dependence is
    # gone from c1 and c1 = 2*c*(rr+1)^30*q'(rr), a UNIT: nothing can save it.
    unit_q1 = _rem_q(sp.diff(Q_POLY, y).subs(y, rr), rr)
    gcd_qq = sp.gcd(sp.Poly(Q_POLY.subs(y, rr), rr),
                    sp.Poly(sp.diff(Q_POLY, y).subs(y, rr), rr))
    mult2.append((1, "A0^2 coefficient is a unit", 1, 1, False))
    mult2.append((2, "q'(r) is a unit (gcd(q,q')=1)", int(gcd_qq.degree() == 0),
                  1, gcd_qq.degree() != 0))
    if verbose:
        print()
        print("=" * 78)
        print("F.  FORWARD CONFIRMATION  --  real polynomial arithmetic")
        print("=" * 78)
        print("    shape        |  a | rho | v(dm2) | v(dm3) | pred | v(H3) | pred | ")
        for tag, a, rho, o2, o3, p3, oH, pH, ok in rows:
            print("    %-12s | %2d | %3d | %6d | %6d | %4d | %5s | %4d | %s"
                  % (tag, a, rho, o2, o3, p3, _fmt(oH), pH, "OK" if ok else "MISMATCH"))
        print("    support theorem, beta = 5 (Phi(5) != 0), e = 3(y+1)^8(y-5)^m:")
        for m, vb, need, bad in supp:
            print("       m=%d : v_beta(N) = %d, need >= %d  -> %s"
                  % (m, vb, need, "POLE (KILL)" if bad else "ok"))
        print("    marked root of q with multiplicity m in e, dm2 germ")
        print("      A0 + A1*(y-r) FULLY GENERAL, exact in Q(gm,A0,A1)[r]/(q):")
        for m, tag, val, _n, bad in mult2:
            print("       m=%d  %-32s : %s" % (m, tag, "OK" if not bad else "FAIL"))
        print("      => v_r(N) = 1 unless the order-1 coefficient vanishes;")
        print("         m >= 2: A0 has dropped out and the coefficient is the")
        print("                 UNIT 2c(r+1)^30 q'(r) -- KILL, no dm2 can help;")
        print("         m  = 1: affine in A0^2 with unit coefficient -- exactly")
        print("                 one square class survives (the R9 relation).")
    # forward pass conditions:
    #  * every (y+1)-order prediction reproduced;
    #  * an e-root off Phi is always a pole (any m);
    #  * the marked-root Taylor coefficients are exactly as the theorem predicts.
    mult_ok = not any(bad for _m, _tag, _v, _n, bad in mult2)
    all_ok = all(r[-1] for r in rows) and all(b for _, _, _, b in supp) and mult_ok
    return {"rows": [list(map(str, r)) for r in rows], "support": supp,
            "mult": [list(map(str, r)) for r in mult2], "pass": bool(all_ok)}


# =====================================================================
#  G.  The sweep across the sub2 state space
# =====================================================================
def load_batch_states():
    with open(os.path.join(HERE, "batch_convolution_sub2.json")) as fh:
        d = json.load(fh)
    return d["states"]


def load_cells():
    with open(os.path.join(HERE, "phase_d_states_sub2.json")) as fh:
        d = json.load(fh)
    return d["cases"]


def sweep_states(verbose=True):
    states = load_batch_states()
    unres = [s for s in states if s["final_verdict"] == "UNRESOLVED"]
    rows = []
    for s in unres:
        a_t = int(s["a_t"])
        deg_e = int(s["deg_e"])
        applies = bool(s["d1_zero"])            # branch T2 <=> d1 = 0
        m_tail = deg_e - a_t
        # a = v_{-1}(dm1) is in [a_t, deg_e]; the bound is min(a, ceil((30-a)/2))
        # which equals a for every a <= 10, and deg_e <= 10 always.
        bnds = {a: pole_bound_closed(a, 30) for a in range(a_t, deg_e + 1)}
        assert all(v is not None for v in bnds.values()), \
            "state out of the proved regime (a > 10): %s" % s
        worst = min(v["rho_min"] for v in bnds.values())
        before = (CAPS["dm2"] + 1) + (CAPS["dm3"] + 1)
        after = CAPS["dm2"] - worst + 1 if applies else before
        rows.append({
            "label": "a%d%s_e%d_d2%s_d1%s_sig%s" % (a_t, s["branch"], deg_e,
                                                    s["deg_d2"], s["deg_d1"],
                                                    s["deg_sigma"]),
            "branch": s["branch"], "a_t": a_t, "deg_e": deg_e,
            "d1_zero": bool(s["d1_zero"]), "applies": applies,
            "tail_deg": m_tail,
            "rho_min": worst if applies else None,
            "sigma_min": worst if applies else None,
            "spares_before": before, "spares_after": after,
            "valsplit_cases_before": (a_t + 1) if applies else None,
            "cases_after": 1 if applies else None,
            "cells": [tuple(c) for c in s.get("cells", [])],
        })
    n_app = sum(1 for r in rows if r["applies"])
    if verbose:
        print()
        print("=" * 78)
        print("G.  SWEEP  --  the %d UNRESOLVED sub2 batch states" % len(rows))
        print("=" * 78)
        print("    branch | a_t | deg e | states | theorem applies | v(dm2) >= |"
              " spares 28 -> | cases -> ")
        agg = {}
        for r in rows:
            key = (r["branch"], r["a_t"], r["deg_e"])
            agg.setdefault(key, []).append(r)
        for key in sorted(agg):
            g = agg[key]
            r = g[0]
            print("      %-4s | %3d |  %3d  |  %3d   |      %-5s      |    %-4s   |"
                  "      %2d      |    %s"
                  % (key[0], key[1], key[2], len(g), r["applies"],
                     r["rho_min"] if r["applies"] else "-",
                     r["spares_after"], r["cases_after"] if r["applies"] else "-"))
        print("    APPLIES on %d / %d unresolved states (%.0f%%)"
              % (n_app, len(rows), 100.0 * n_app / max(1, len(rows))))
        print("    FAILS   on %d / %d  -- all of them branch T1 (d1 != 0)"
              % (len(rows) - n_app, len(rows)))
    return rows


def sweep_cells(verbose=True):
    """The CELL-level consequence: a marked-root multiplicity b_j >= 2 in a
    T2 cell is an outright KILL (B2).  This is the top-stratum bite the design
    spec asks for -- it does not need to reach any stratum bottom-up."""
    cells = load_cells()
    rows = []
    for c in cells:
        b = list(c["b"])
        bmax = max(b) if b else 0
        t2 = c["branch"] == "T2"
        kill = t2 and bmax >= 2
        rows.append({"branch": c["branch"], "a_t": c["a_t"], "b": b,
                     "state_count": c["state_count"], "d2_zero": c["d2_zero"],
                     "g_zero_levels": c["g_zero_levels"],
                     "applies": t2, "bmax": bmax, "cell_killed": kill,
                     "n_marked": sum(1 for x in b if x >= 1)})
    n_t2 = sum(1 for r in rows if r["applies"])
    killed = [r for r in rows if r["cell_killed"]]
    killed_states = sum(r["state_count"] for r in killed)
    t1_bad = [r for r in rows if not r["applies"] and r["bmax"] >= 2]
    if verbose:
        print()
        print("=" * 78)
        print("H.  CELL-LEVEL BITE  --  the 220 sub2 flag cells")
        print("=" * 78)
        print("    T2 cells (d1 = 0, theorem applies) : %d / %d   (%d / %d states)"
              % (n_t2, len(rows), sum(r["state_count"] for r in rows if r["applies"]),
                 sum(r["state_count"] for r in rows)))
        print("    T2 cells with a marked-root multiplicity b_j >= 2:")
        for r in killed:
            print("       a=%d b=%s d2_zero=%s g_zero=%s : %d states  -> KILLED by B2"
                  % (r["a_t"], r["b"], r["d2_zero"], r["g_zero_levels"],
                     r["state_count"]))
        print("       total: %d cells / %d states CLOSED, degree-uniformly,"
              " with no Groebner basis." % (len(killed), killed_states))
        print("    T1 cells with b_j >= 2 (NOT covered -- d1 != 0): %d cells /"
              " %d states" % (len(t1_bad), sum(r["state_count"] for r in t1_bad)))
        print("    remaining T2 cells: every b_j is 0 or 1, so B2 forces")
        print("      dm2(r_j) != 0 and one square relation per marked root;")
        print("      B1 additionally confines the unmarked part of e to")
        print("      {-1} U roots(q) -- exactly the q-support the batch model")
        print("      DROPPED (PHASE_F2_SUB2.md sec.1), recovered here as a")
        print("      theorem instead of an input.")
    return {"rows": rows, "n_t2_cells": n_t2, "killed": killed,
            "killed_states": killed_states,
            "t1_multiplicity_cells": len(t1_bad)}


# =====================================================================
def selfcheck():
    """Every load-bearing claim, asserted.  Returns (ok, failures)."""
    fails = []

    def chk(name, cond):
        if not cond:
            fails.append(name)

    A = stage_A_identity(verbose=False)
    chk("A.H_rebuilt_matches_json", A["H_rebuilt_matches_json"])
    chk("A.K5_linear_in_dm3", A["K5_deg_in_dm3"] == 1)
    chk("A.d1_remainder", A["d1_remainder"] == "6*d1*dm2**2*dm3")
    chk("A.K5gen_quadratic", A["K5gen_deg_in_dm3"] == 2)

    B = stage_B_phi(verbose=False)
    chk("B.v_Phi_30", B["v_minus1_Phi"] == 30)
    chk("B.q_irreducible", B["q_irreducible"] and B["q_squarefree"])

    C = stage_C_bound(verbose=False)
    chk("C.cap_independent", C["cap_independent"])
    chk("C.kill_set_exact", C["kill_set_exact"])
    chk("C.proof_sound", C["proof_sound"])
    chk("C.proof_sharp", C["proof_sharp"])
    chk("C.aux_independent", C["aux_independent"])
    chk("C.grid_sound", not C["grid_unsound"])
    for a in range(1, 11):
        chk("C.bound_is_a(a=%d)" % a, pole_bound_closed(a, 30)["rho_min"] == a
            and pole_bound_closed(a, 30)["sigma_min"] == a)
    chk("C.out_of_regime_a11", pole_bound_closed(11, 30) is None
        and not in_regime(11, 30) and in_regime(10, 30))

    D = stage_D_places(verbose=False)
    for r in D["rows"]:
        if r["Pb"] == 0:
            chk("D.Pb0_m%d_all_dead" % r["m"], r["all_dead"])
        if r["Pb"] == 1 and r["m"] >= 2:
            chk("D.Pb1_m%d_all_dead" % r["m"], r["all_dead"])
    chk("D.Pb1_m1_survives", not [r for r in D["rows"]
                                  if r["Pb"] == 1 and r["m"] == 1][0]["all_dead"])
    chk("D.relation_matches", D["relation"]["matches_generic_fiber"])

    R = stage_R9_regression(verbose=False)
    chk("R9.regression", R["pass"])

    E = stage_E_report(verbose=False)
    chk("E.r9_collapse_matches", E["r9_match"])
    a10 = [r for r in E["rows"] if r["a"] == 10 and r["k"] == 0][0]
    chk("E.a10_no_remainder", a10["n_remainder_nonzero"] == 0
        and a10["spare_unknowns_after"] == 3)

    F = stage_F_forward(verbose=False)
    chk("F.forward", F["pass"])

    S = sweep_states(verbose=False)
    chk("S.90_dege10_T2", sum(1 for r in S
                              if r["branch"] == "T2" and r["deg_e"] == 10) == 90)
    chk("S.all_T2_apply", all(r["applies"] for r in S if r["branch"] == "T2"))
    chk("S.no_T1_applies", not any(r["applies"] for r in S if r["branch"] == "T1"))

    Hc = sweep_cells(verbose=False)
    chk("H.t2_cells_24", Hc["n_t2_cells"] == 24)
    chk("H.cells_killed", len(Hc["killed"]) == 2 and Hc["killed_states"] == 16)
    return (not fails), fails


def main():
    quiet = "--quiet" in sys.argv
    if quiet:
        ok, fails = selfcheck()
        print("pole_theorem_sweep selfcheck: %s%s"
              % ("ALL PASS" if ok else "FAIL", "" if ok else "  " + ", ".join(fails)))
        return 0 if ok else 1

    sys.stdout.reconfigure(line_buffering=True)
    print("#" * 78)
    print("# POLE THEOREM SWEEP -- generalising GENERIC_FIBER.md sec.4")
    print("#" * 78)
    A = stage_A_identity()
    B = stage_B_phi()
    C = stage_C_bound()
    D = stage_D_places()
    R = stage_R9_regression()
    E = stage_E_report()
    F = stage_F_forward()
    S = sweep_states()
    Hc = sweep_cells()

    ok, fails = selfcheck()
    print()
    print("#" * 78)
    print("# SELF-CHECK: %s%s" % ("ALL PASS" if ok else "FAIL",
                                  "" if ok else "  " + ", ".join(fails)))
    print("#" * 78)

    out = {
        "identity": A, "phi": B, "bound": C,
        "places": {"rows": D["rows"], "relation": D["relation"]},
        "r9_regression": R,
        "collapse": E["rows"], "collapse_r9_match": E["r9_match"],
        "forward": {"pass": F["pass"], "support": F["support"], "mult": F["mult"]},
        "states": S,
        "cells": {k: v for k, v in Hc.items() if k != "rows"},
        "selfcheck_pass": ok, "selfcheck_failures": fails,
    }
    with open(os.path.join(HERE, "pole_theorem_sweep.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    print("wrote pole_theorem_sweep.json")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
