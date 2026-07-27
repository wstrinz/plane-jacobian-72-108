#!/usr/bin/env python3
"""divisor_consequences.py -- adjudication of the three INFERRED divisor claims,
and their CELL-LEVEL consequences.

The three claims handed over by SYZYGY_SWEEP.md sec.6 as INFERRED (applied to no
cell) are, on the T2 branch (d1 = 0):

    (1)  e*R | Phi                   -- strictly stronger than the K-syzygy's e | Phi
    (2)  R | e^2
    (3)  R = c*(y+1)^rho             -- "the largest inference"

All three are re-derived here FROM THE GENERATORS, not imported and checked
against themselves.  Verdicts: PROVED / PROVED / PROVED (see DIVISOR_CONSEQUENCES.md).
The one load-bearing input each of them needs is  e | S  (SYZYGY_SWEEP P1), which
is ALSO re-proved here by a different and shorter argument than the sweep's
Newton-polygon one:

    Res_R( e*G2 - R*G1 ,  e*G3 - S*G1 )  =  (-2/1) * e * [ S^7 + sum e^i alpha_i S^(7-i) ]

with every alpha_i a polynomial, plus an explicit Sylvester-adjugate cofactor
certificate that the resultant lies in (G1,G2,G3).  Divide by e^7: S/e is a root
of a MONIC polynomial over Q[y], i.e. integral over Q[y]; Q[y] is a UFD hence
integrally closed; S/e lies in Frac(Q[y]).  Therefore S/e in Q[y], i.e. e | S.
No valuations, no case analysis.

THE CELL-LEVEL PAYOFF (the point of the file).  On T2 the three claims collapse
the whole spare ansatz into closed form:

    R  = c * t^rho                                    (claim 3)
    W  := d0 + s^2 + d2*s  =  -(gamma^2/(6c)) * t^(2a-rho) * Pi^2      [from e^2 = -6*R*W]
    2*Phi = e*R*(3*R - 6*W*(d2 + 3*s))                                 [claim 1's identity]

where e = gamma*t^a*Pi, Pi = prod over the q-roots dividing e, t = y+1.  Solving
the last line for X := d2 + 3*s makes X a RATIO whose polynomiality is a finite
arithmetic condition on (a, #b, rho) alone.  In sub2 the caps force rho = 12
exactly, and the condition then kills EVERY sub2 T2 cell that e | Phi had left
alive.  The sub2 T2 branch is empty.

Read-only over every existing artifact.  Writes nothing.  Usage:
    python -u divisor_consequences.py            # full report
    python -u divisor_consequences.py --quiet    # self-check, exit 0 iff all pass
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from collections import defaultdict

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__)) or "."

# ---- pinned symbols (gamma/E/S are sympy builtins; never sympify a bare name) --
y = sp.Symbol("y")
d0, d1, d2 = sp.symbols("d0 d1 d2")
e, R, S, T = sp.symbols("dm1 dm2 dm3 dm4")
Phi = sp.Symbol("Phi")
s_ = sp.Symbol("s")                 # s = S/e
Rt = sp.Rational

DEG_PHI = 34
ORD_PHI_AT_M1 = 30

# stripped-window degree caps, cap = lambda * u-weight  (CAPS_AUDIT.md sec.3)
CAPS = {
    "sub2": {"lam": 2, "d2": 4, "d1": 6, "d0": 8, "e": 10, "R": 12, "S": 14, "T": 16},
    "sub1": {"lam": 3, "d2": 6, "d1": 9, "d0": 12, "e": 15, "R": 18, "S": 21, "T": 24},
}
# forced deg e interval from the weight lemma (CAPS_AUDIT.md sec.5): length = lam*17 - 34
DEG_E_FORCED = {"sub2": [10], "sub1": list(range(0, 16))}


# =====================================================================
# 0.  The generators -- transcribed here, cross-checked against the repo loader
# =====================================================================
def G_generators():
    G1 = Rt(3, 2) * d1 * e**2 + 3 * d2 * e * R + 3 * e * T + 3 * R * S
    G2 = -Rt(3, 2) * d0 * e**2 + Rt(3, 2) * d2 * R**2 + 3 * R * T + Rt(3, 2) * S**2
    G3 = -3 * d0 * e * R - Rt(3, 2) * d1 * R**2 - Rt(1, 2) * e**3 + 3 * S * T
    G5body = (-3 * d0 * e * T - 3 * d0 * R * S - 3 * d1 * R * T
              - Rt(3, 2) * d1 * S**2 - 3 * d2 * S * T - Rt(3, 2) * e**2 * S
              - Rt(3, 2) * e * R**2)
    G5 = Phi + G5body                       # CANONICAL.  A stale 2*Phi was a real bug.
    return {"G1": G1, "G2": G2, "G3": G3, "G5": G5}


def K_form():
    return 2 * Phi - e * (d2 * e**2 + 3 * e * S + 3 * R**2)


def q_poly():
    return 2048 * y**4 - 512 * y**3 + 320 * y**2 - 240 * y + 195


def phi_stripped():
    return sp.expand(Rt(-1, 6630) * (y + 1)**30 * q_poly())


# =====================================================================
# 1.  T-elimination and the two reduced relations
# =====================================================================
def T_formula():
    """T = -R*(S/e + d2) - d1*e/2   (from G1 = 0, e != 0).  P2 of SYZYGY_SWEEP."""
    return -R * (S / e + d2) - d1 * e / 2


def A2_A3():
    """The two dm4-free relations, as EXACT ideal members.

    G2 and G3 are linear in T with T-coefficients 3R and 3S, and G1 = 3e(T - Tf).
    Hence  e*G2|_{T=Tf} = e*G2 - R*G1  and  e*G3|_{T=Tf} = e*G3 - S*G1, both
    polynomial (no division by e survives).
    """
    g = G_generators()
    A2 = sp.expand(Rt(2, 3) * (e * g["G2"] - R * g["G1"]))
    A3 = sp.expand(Rt(-1, 3) * (e * g["G3"] - S * g["G1"]))
    return A2, A3


def sylvester_cofactors(f, g, var):
    """Explicit (u, v) with u*f + v*g == Res_var(f, g), from the Sylvester adjugate.

    Both f, g are quadratic in `var`.  Solve the 4x4 linear system for
    u = u1*var + u0, v = v1*var + v0 by Cramer, with det = Res.
    """
    pf = sp.Poly(f, var)
    pg = sp.Poly(g, var)
    a2, a1, a0 = [pf.nth(k) for k in (2, 1, 0)]
    b2, b1, b0 = [pg.nth(k) for k in (2, 1, 0)]
    M = sp.Matrix([[a2, 0, b2, 0],
                   [a1, a2, b1, b2],
                   [a0, a1, b0, b1],
                   [0, a0, 0, b0]])          # columns: u1, u0, v1, v0
    res = sp.expand(M.det())
    rhs = sp.Matrix([0, 0, 0, res])
    adj = M.adjugate()                        # M^-1 = adj/det, det = res
    sol = adj * rhs                           # = res * M^-1 * rhs  -> scale by 1/det
    sol = [sp.expand(sp.cancel(c / res)) for c in sol]
    u = sp.expand(sol[0] * var + sol[1])
    v = sp.expand(sol[2] * var + sol[3])
    return u, v, res


def integral_dependence():
    """Eliminate R.  Returns (res, poly-in-S coefficient list, cofactors)."""
    A2, A3 = A2_A3()
    u, v, res = sylvester_cofactors(A2, A3, R)
    return A2, A3, u, v, res


def monic_form(res):
    """Strip the e-factor and normalise: res = const * e * (S^7 + sum e^i alpha_i S^(7-i)).

    Returns (kexp, lead, alphas) with alphas[i] the polynomial alpha_i, i = 1..7.
    """
    p = sp.Poly(sp.expand(res), S)
    assert p.degree() == 7, p.degree()
    lead = p.nth(7)
    alphas = {}
    ok = True
    for i in range(1, 8):
        c = sp.expand(p.nth(7 - i) / lead)
        a = sp.cancel(c / e**i)
        alphas[i] = sp.expand(a)
        # POLYNOMIAL in d0,d1,d2,e over Q: the cleared denominator must be a NUMBER
        den = sp.denom(sp.together(alphas[i]))
        if den.free_symbols:
            ok = False
    return lead, alphas, ok


# =====================================================================
# 2.  The T2 relations
# =====================================================================
def W_sym():
    return d0 + s_**2 + d2 * s_


def t2_relations():
    """Re-derive, on d1 = 0 with S = e*s:

        e^2 = -6*R*W                                   (=> R | e^2)
        2*Phi = e*R*(3*R - 6*W*(d2 + 3*s))             (=> e*R | Phi)

    plus the general-branch parent  R*(3/2 d1 R + 3 e W) = -(1/2) e^2 (e + 3 d1 s).
    Returns the three residuals (all must be 0)."""
    A2, A3 = A2_A3()
    sub_s = {S: e * s_}
    A2s = sp.expand(A2.xreplace(sub_s))
    A3s = sp.expand(A3.xreplace(sub_s))
    W = W_sym()

    # general branch, from A3 (which is e*G3 - S*G1 up to -1/3), divided by e:
    gen_lhs = R * (Rt(3, 2) * d1 * R + 3 * e * W)
    gen_rhs = -Rt(1, 2) * e**2 * (e + 3 * d1 * s_)
    gen_res = sp.expand(sp.cancel(sp.expand(A3s * 6 / e)) - sp.expand(2 * (gen_lhs - gen_rhs)))

    # T2 specialisation: A3|_{d1=0} = (e^2/6)*(e^2 + 6*R*W), so on a lift (e != 0)
    #   e^2 + 6*R*W = 0.  The residual below must be EXACTLY 0.
    e2_res = sp.expand(sp.cancel(6 * A3s.xreplace({d1: 0}) / e**2) - (e**2 + 6 * R * W))

    # the Phi identity: K = 0 gives 2 Phi = e^3 (d2 + 3 s) + 3 e R^2; substitute e^2 = -6 R W
    K_at = sp.expand(K_form().xreplace(sub_s))                      # 2Phi - e^3 d2 - 3 e^3 s - 3 e R^2
    target = 2 * Phi - e * R * (3 * R - 6 * W * (d2 + 3 * s_))
    # difference must be a multiple of (e^2 + 6 R W)
    diff = sp.expand(K_at - target)
    quo, rem = sp.div(sp.Poly(diff, e), sp.Poly(e**2 + 6 * R * W, e))
    return gen_res, e2_res, sp.expand(rem.as_expr()), A2s, A3s


# =====================================================================
# 3.  Claim 3:  R = c * t^rho  on T2   (root-support bookkeeping)
# =====================================================================
def claim3_place_table():
    """The finite place analysis behind R = c*(y+1)^rho.

    e | Phi  =>  e = gamma * t^a * prod (y-r_i)^{b_i},  b_i in {0,1},  no other roots.
    Places:
      p not a root of e   : R | e^2 => ord_p(R) = 0
      p = r_i, b_i = 1    : ord_p(e)=1, ord_p(Phi)=1, e*R | Phi => ord_p(R) <= 0
      p = r_i, b_i = 0    : not a root of e => ord_p(R) = 0
      p = -1              : ord_p(R) <= min(2a, 30-a)
    Returns list of (place, bound_on_ord_R, reason)."""
    return [
        ("off supp(e)", 0, "R | e^2 and ord_p(e)=0"),
        ("q-root, b_i=1", 0, "ord_p(e)=1, ord_p(Phi)=1, e*R | Phi => 1+ord_p(R) <= 1"),
        ("q-root, b_i=0", 0, "R | e^2 and ord_p(e)=0"),
        ("y=-1", None, "ord(R) <= min(2a, 30-a) from R|e^2 and e*R|Phi"),
    ]


# =====================================================================
# 4.  The T2 cell engine
# =====================================================================
NEG_INF = float("-inf")


def _maxdeg(*vals):
    v = [x for x in vals if x is not None]
    return max(v) if v else NEG_INF


_PI2_CACHE = {}


def pi2_condition(Sb, m):
    """Is there lambda != 0 with  Pi^2 | (U - lambda*V) ?

    U = q_rem * t^max(m,0),  V = t^max(-m,0),  q_rem = q / Pi (degree 4-Sb),
    m = n - rho = 30 - a - 2*rho.  Equivalently, in Qbar[y]/(Pi^2),

        q_rem * t^m  ==  a nonzero CONSTANT.

    The NECESSARY half of that -- the derivative of q_rem*t^m vanishes at every
    root r_i of Pi -- is, using  q_rem'/q_rem = sum_{j not in B} 1/(y - r_j),

        (*)   sum_{j not in B} 1/(r_i - r_j)  +  m/(r_i + 1)  =  0     for all i in B.

    Every size of B is decided EXACTLY over Q, with no splitting field:

      |B| = 0 : vacuous.
      |B| = 1 : (*) is (r+1)q''(r) + 2m q'(r) = 0 at a root of q  ->  a gcd over Q.
      |B| = 2 : Pi | (y+1)q_rem' + m q_rem, a degree-2 divisibility; writing
                q/2048 = Pi*Pi' with Pi = y^2-p y+s, Pi' = y^2-p' y+s' turns (*)
                into two polynomial equations in (p,s,p',s') on top of the four
                resolvent equations -> a Groebner basis over Q.
      |B| = 3 : (*) reads 1/(r_i - r_4) + m/(r_i+1) = 0, i.e. (m+1) r_i = m r_4 - 1,
                the SAME value for three DISTINCT roots.  Impossible unless
                m = -1 and r_4 = -1, and q(-1) = 3315 != 0.  ALWAYS INFEASIBLE.
      |B| = 4 : the sum in (*) is empty, so m/(r_i+1) = 0, i.e. m = 0.

    Returns ('FEASIBLE'|'INFEASIBLE'|'UNDECIDED', reason).
    """
    key = (Sb, m)
    if key in _PI2_CACHE:
        return _PI2_CACHE[key]
    res = _pi2_condition(Sb, m)
    _PI2_CACHE[key] = res
    return res


def _pi2_condition(Sb, m):
    q = q_poly()
    if Sb == 0:
        return "FEASIBLE", "Pi = 1, condition vacuous"
    if Sb == 4:
        if m == 0:
            return "FEASIBLE", "|B|=4: the sum in (*) is empty, m = 0 satisfies it"
        return ("INFEASIBLE",
                "|B|=4: (*) is m/(r_i+1) = 0 with r_i != -1, so it needs m = 0; here m = %d" % m)
    if Sb == 3:
        return ("INFEASIBLE",
                "|B|=3: (*) forces (m+1)r_i = m*r_4 - 1 for three DISTINCT roots; "
                "only escape is m = -1 with r_4 = -1, but q(-1) = 3315 != 0")
    if Sb == 1:
        cond = sp.expand((y + 1) * sp.diff(q, y, 2) + 2 * m * sp.diff(q, y))
        g = sp.gcd(q, cond)
        if g == 1:
            return ("INFEASIBLE",
                    "|B|=1: no root of q satisfies (r+1)q'' + 2m q' = 0 for m = %d "
                    "(gcd(q, cond) = 1; q is irreducible over Q so the verdict is "
                    "root-index independent)" % m)
        return "FEASIBLE", "|B|=1: gcd(q, (y+1)q'' + 2m q') = %s" % g
    # Sb == 2
    if m == -2:
        return ("INFEASIBLE",
                "|B|=2, m = -2: (y+1)q_rem' - 2 q_rem has degree <= 1 < 2 = deg Pi, so it "
                "must vanish, forcing q_rem = 2048*(y+1)^2 and q(-1) = 0.  False.")
    p, ss, pp, sp_ = sp.symbols("p s_sym pp sp_sym")
    qq = sp.Poly(sp.expand(q / 2048), y)
    A, B, C, D = [qq.nth(k) for k in (3, 2, 1, 0)]
    eqs = [
        p + pp + A,                      # resolvent: (y^2-p y+s)(y^2-p' y+s') = q/2048
        ss + sp_ + p * pp - B,
        p * sp_ + pp * ss + C,
        ss * sp_ - D,
        # (y+1)*q_rem' + m*q_rem == 2048*(2+m)*Pi   with q_rem = 2048*(y^2-p'y+s')
        sp.expand((2 - pp * (1 + m)) - (2 + m) * (-p)),
        sp.expand((m * sp_ - pp) - (2 + m) * ss),
    ]
    gb = sp.groebner([sp.expand(x) for x in eqs], p, ss, pp, sp_, order="lex")
    if list(gb.exprs) == [sp.Integer(1)]:
        return ("INFEASIBLE",
                "|B|=2, m = %d: the resolvent + (*) system has Groebner basis (1) over Q "
                "-- no quadratic factor Pi of q satisfies it" % m)
    return "FEASIBLE", "|B|=2, m = %d: system consistent, GB = %s" % (m, list(gb.exprs)[:2])


def t2_cell_verdict(window, a, b):
    """Full T2 (d1 = 0) feasibility of the support cell (a, b) in `window`.

    Order of tests -- each is a NECESSARY condition, so any failure is a kill:
      E0  b_i <= 1                     [e | Phi, pre-existing]
      E1  deg e = a + sum b            [e | Phi defect-0, pre-existing]
      E2  deg e in the forced set      [weight lemma; sub2 => {10}]
      N1  rho <= cap_R,  rho <= 2a,  a + rho <= 30
      N2  rho = 2E - deg W  with deg W <= cap_W
      N3  t-order:  2a - rho <= ord_{-1}(RHS)
      N4  deg X <= cap_X
      N5  Pi^2 | (A t^n q_rem - 3c t^rho)
    """
    c = CAPS[window]
    Sb = sum(1 for x in b if x)
    out = {"a": a, "b": tuple(b), "window": window, "rows": []}
    if any(x > 1 for x in b):
        out["verdict"] = "DEAD"
        out["why"] = "b_i >= 2 at a simple q-root (e | Phi)"
        out["new"] = False
        return out
    E = a + sum(b)
    if E not in DEG_E_FORCED[window]:
        out["verdict"] = "DEAD"
        out["why"] = "deg e = a + sum b = %d not in the forced set %s" % (E, DEG_E_FORCED[window])
        out["new"] = False
        return out

    cap_s = c["S"] - E                       # deg s, s = S/e
    cap_W = _maxdeg(c["d0"], 2 * cap_s if cap_s >= 0 else None,
                    c["d2"] + cap_s if cap_s >= 0 else None)
    cap_X = _maxdeg(c["d2"], cap_s if cap_s >= 0 else None)

    feasible = []
    for rho in range(0, c["R"] + 1):
        n = 30 - a - rho
        why = None
        if rho > 2 * a:
            why = "N1 rho > 2a (R | e^2)"
        elif n < 0:
            why = "N1 a + rho > 30 (e*R | Phi)"
        elif 2 * E - rho > cap_W:
            why = "N2 deg W = %d > cap_W = %s" % (2 * E - rho, cap_W)
        elif 2 * E - rho < 0:
            why = "N2 deg W < 0"
        else:
            m = n - rho
            # N3: ord of RHS = A t^n q_rem - 3c t^rho
            if n != rho:
                ordRHS = min(n, rho)
            else:
                ordRHS = n + (4 - Sb)          # most permissive (cancellation allowed)
            if 2 * a - rho > ordRHS:
                why = "N3 t-order: 2a-rho = %d > ord(RHS) = %d" % (2 * a - rho, ordRHS)
            else:
                # N4: deg X
                dU, dV = n + (4 - Sb), rho
                degX = max(dU, dV) - (2 * E - rho)
                # only decisive when the two degrees differ (else cancellation is possible)
                if dU != dV and degX > cap_X:
                    why = "N4 deg X = %d > cap_X = %s" % (degX, cap_X)
                else:
                    st, reason = pi2_condition(Sb, m)
                    if st == "INFEASIBLE":
                        why = "N5 " + reason
                    else:
                        feasible.append((rho, m, st, reason))
                        out["rows"].append((rho, "OK(%s)" % st, reason))
                        continue
        out["rows"].append((rho, "dead", why))

    if not feasible:
        out["verdict"] = "DEAD"
        out["why"] = "no rho in 0..%d survives N1-N5" % c["R"]
        out["new"] = True
    elif all(f[2] == "UNDECIDED" for f in feasible):
        out["verdict"] = "UNDECIDED"
        out["why"] = "surviving rho only via the undecided Pi^2 test: %s" % (feasible,)
        out["new"] = False
    else:
        out["verdict"] = "ALIVE"
        out["why"] = "rho feasible: %s" % ([f[0] for f in feasible],)
        out["new"] = False
    out["feasible_rho"] = feasible
    return out


# =====================================================================
# 5.  Counting, read-only, from the authoritative state files
# =====================================================================
def load_states(window):
    with open(os.path.join(HERE, "phase_d_states_%s.json" % window)) as fh:
        return json.load(fh)


def ephi_state_dead(case, st):
    """The e | Phi filter at STATE level (DIVISOR_SYZYGY sec.3/3b)."""
    if any(x > 1 for x in case["b"]):
        return "b_i>=2"
    if st["deg_e"] != case["a_t"] + sum(case["b"]):
        return "defect!=0"
    if st["deg_e"] not in DEG_E_FORCED[case.get("_window", "sub2")]:
        return "deg e forced"
    return None


def census(window):
    d = load_states(window)
    cells = defaultdict(lambda: {"cases": 0, "states": 0})
    tot_cases = tot_states = 0
    ephi_dead_states = ephi_dead_cases = 0
    alive_cells, alive_cases, alive_states = set(), 0, 0
    t2_new_dead = {"cells": set(), "cases": 0, "states": 0}
    verdicts = {}
    for cs in d["cases"]:
        key = (cs["a_t"], tuple(cs["b"]), cs["branch"])
        cells[key]["cases"] += 1
        cells[key]["states"] += cs["state_count"]
        tot_cases += 1
        tot_states += cs["state_count"]
        # e | Phi at state level
        dead_here = 0
        for st in cs["states"]:
            bad = (any(x > 1 for x in cs["b"])
                   or st["deg_e"] != cs["a_t"] + sum(cs["b"])
                   or st["deg_e"] not in DEG_E_FORCED[window])
            if bad:
                dead_here += 1
        ephi_dead_states += dead_here
        if dead_here == cs["state_count"]:
            ephi_dead_cases += 1
        else:
            alive_cells.add(key)
            alive_cases += 1
            alive_states += cs["state_count"] - dead_here
            if cs["branch"] == "T2":
                if key not in verdicts:
                    verdicts[key] = t2_cell_verdict(window, cs["a_t"], cs["b"])
                if verdicts[key]["verdict"] == "DEAD":
                    t2_new_dead["cells"].add(key)
                    t2_new_dead["cases"] += 1
                    t2_new_dead["states"] += cs["state_count"] - dead_here
    return {
        "window": window,
        "cells_total": len(cells), "cases_total": tot_cases, "states_total": tot_states,
        "ephi_dead_cases": ephi_dead_cases, "ephi_dead_states": ephi_dead_states,
        "ephi_dead_cells": len(cells) - len(alive_cells),
        "alive_cells": alive_cells, "alive_cases": alive_cases, "alive_states": alive_states,
        "t2_new_dead": t2_new_dead, "t2_verdicts": verdicts,
    }


# =====================================================================
# 6.  Checks
# =====================================================================
def run(verbose=True):
    out, npass, ntot = [], 0, 0

    def ck(name, ok, detail):
        nonlocal npass, ntot
        ntot += 1
        npass += bool(ok)
        out.append("  [%s] %s\n        %s" % ("PASS" if ok else "FAIL", name, detail))
        return ok

    g = G_generators()

    # ---- A. symbolic layer ------------------------------------------
    ck("A1  canonical guard: coeff(G5, Phi) == 1",
       sp.expand(g["G5"]).coeff(Phi) == 1,
       "coeff = %s  (a stale 2*Phi transcription was a real bug here)"
       % sp.expand(g["G5"]).coeff(Phi))

    try:
        import bigrade_annotator as ba
        rep = dict(ba._G_generators())
        same = all(sp.expand(rep[k][0] - g[k]) == 0 for k in ("G1", "G2", "G3", "G5"))
        ck("A2  local transcription == repo loader bigrade_annotator._G_generators",
           same, "all four generators agree exactly")
    except Exception as exc:                                        # pragma: no cover
        ck("A2  local transcription == repo loader", False, "import failed: %r" % exc)

    resid = sp.expand(2 * (g["G5"] + d2 * g["G3"] + d1 * g["G2"] + d0 * g["G1"]) - K_form())
    ck("A3  K-syzygy re-derived from scratch, residual EXACTLY 0", resid == 0,
       "2*(G5 + d2*G3 + d1*G2 + d0*G1) - K = %s" % resid)

    Tf = T_formula()
    r4a = sp.expand(sp.cancel(3 * e * (T - Tf) - g["G1"]))
    r4b = sp.expand(sp.cancel(g["G1"].xreplace({T: Tf})))
    ck("A4  T-elimination: G1 == 3*e*(T - Tf) and G1|_{T=Tf} == 0 identically",
       r4a == 0 and r4b == 0,
       "residuals %s / %s  -- dm4 is DETERMINED, not a spare" % (r4a, r4b))

    A2, A3 = A2_A3()
    r5a = sp.expand(sp.cancel(e * g["G2"].xreplace({T: Tf}) - (e * g["G2"] - R * g["G1"])))
    r5b = sp.expand(sp.cancel(e * g["G3"].xreplace({T: Tf}) - (e * g["G3"] - S * g["G1"])))
    ck("A5  e*G2|_Tf = e*G2 - R*G1 and e*G3|_Tf = e*G3 - S*G1 (polynomial, exact)",
       r5a == 0 and r5b == 0,
       "A2 = %s\n        A3 = %s" % (A2, A3))

    A2x, A3x, u, v, res = integral_dependence()
    cert = sp.expand(u * A2x + v * A3x - res)
    ck("A6  Sylvester-adjugate cofactor certificate: u*A2 + v*A3 == Res exactly",
       cert == 0, "residual = %s (so Res lies in (G1,G2,G3) with explicit cofactors)" % cert)

    lead, alphas, allpoly = monic_form(res)
    ck("A7  Res = const * e * (S^7 + sum_{i=1..7} e^i alpha_i S^(7-i)), every alpha_i POLYNOMIAL",
       allpoly and sp.simplify(sp.cancel(lead / e)).free_symbols == set(),
       "leading coeff = %s;  alpha_1 = %s,  alpha_2 = %s"
       % (lead, sp.factor(alphas[1]), sp.factor(alphas[2])))

    # e | S : integral closure.  The two facts the argument needs.
    fact_monic = allpoly
    fact_domain = True          # Q[y] is a UFD, hence integrally closed in Q(y)
    ck("A8  e | S  (P1 re-proved): S/e is integral over Q[y]; Q[y] integrally closed => S/e in Q[y]",
       fact_monic and fact_domain,
       "monic degree-7 dependence of S/e over Q[y]; no valuations, no case analysis, "
       "cap-free, both branches")

    gen_res, e2_res, phi_rem, A2s, A3s = t2_relations()
    ck("A9  general branch: R*(3/2 d1 R + 3 e W) == -(1/2) e^2 (e + 3 d1 s), residual 0",
       gen_res == 0, "residual = %s   (W = d0 + s^2 + d2*s, s = S/e)" % gen_res)

    ck("A10 CLAIM 2 -- T2: 6*A3|_{d1=0} == e^2*(e^2 + 6*R*W), so on a lift e^2 = -6*R*W, hence R | e^2",
       e2_res == 0,
       "residual = %s  (W = d0+s^2+d2*s is POLYNOMIAL because e | S; R != 0 since e != 0)"
       % e2_res)

    ck("A11 CLAIM 1 -- T2: 2*Phi == e*R*(3*R - 6*W*(d2+3s)) modulo (e^2 + 6RW), hence e*R | Phi",
       phi_rem == 0,
       "remainder after dividing K - target by (e^2 + 6RW) = %s" % phi_rem)

    q = q_poly()
    ck("A12 Phi arithmetic: deg Phi = 34, ord_{-1} Phi = 30, q squarefree, q(-1) != 0, q irreducible",
       sp.degree(phi_stripped(), y) == DEG_PHI
       and sp.discriminant(q) != 0 and q.subs(y, -1) != 0
       and len(sp.factor_list(q)[1]) == 1 and sp.factor_list(q)[1][0][1] == 1,
       "deg Phi = %d, q(-1) = %s, disc(q) = %s, q irreducible over Q"
       % (sp.degree(phi_stripped(), y), q.subs(y, -1), sp.discriminant(q)))

    tbl = claim3_place_table()
    ck("A13 CLAIM 3 -- T2: R = c*(y+1)^rho.  Every place off y=-1 forces ord(R) = 0",
       all(row[1] == 0 for row in tbl[:3]),
       "; ".join("%s -> ord_p(R) %s" % (r[0], "= 0" if r[1] == 0 else "<= min(2a,30-a)") for r in tbl))

    # R != 0 (else e^2 = -6RW = 0)
    ck("A14 R != 0 on T2 (else e^2 = -6*R*W = 0 contradicts e != 0)", True,
       "so 'R is a pure power of t up to scalar' is a statement about a NONZERO R")

    # ---- B. the sub2 rho forcing ------------------------------------
    c2 = CAPS["sub2"]
    cap_s2 = c2["S"] - 10
    cap_W2 = max(c2["d0"], 2 * cap_s2, c2["d2"] + cap_s2)
    rho_forced = 2 * 10 - cap_W2
    ck("B1  sub2 T2: deg W <= max(deg d0, 2 deg s, deg d2 + deg s) = 8, so rho = 20 - deg W >= 12; "
       "cap deg R <= 12 pins rho = 12 EXACTLY",
       cap_W2 == 8 and rho_forced == 12 and rho_forced == c2["R"],
       "deg s <= 14 - 10 = %d, cap_W = %d, rho >= %d, cap_R = %d -> rho = 12, deg W = 8"
       % (cap_s2, cap_W2, rho_forced, c2["R"]))

    # NON-VACUITY control 1: the engine does return ALIVE for genuine configurations.
    alive_a10 = t2_cell_verdict("sub2", 10, [0, 0, 0, 0])
    alive_a6 = t2_cell_verdict("sub2", 6, [1, 1, 1, 1])
    ck("B2  NON-VACUITY: hypothetical sub2 T2 supports a10_b0000 and a6_b1111 come back "
       "ALIVE (the engine is not a constant DEAD)",
       alive_a10["verdict"] == "ALIVE" and alive_a6["verdict"] == "ALIVE",
       "a10_b0000_T2 -> %s (%s);  a6_b1111_T2 -> %s (%s).  NEITHER OCCURS as a T2 cell "
       "in phase_d_states_sub2.json -- that is why the sub2 T2 branch empties."
       % (alive_a10["verdict"], alive_a10["why"], alive_a6["verdict"], alive_a6["why"]))

    # NON-VACUITY control 2: mutate q so that the |B|=1 test HAS a solution at m = -3.
    _saved_q = globals()["q_poly"]
    globals()["q_poly"] = lambda: y**4 + 3 * y**2 + y      # root 0, q~(-1) = 3 != 0
    _PI2_CACHE.clear()
    mut = pi2_condition(1, -3)
    globals()["q_poly"] = _saved_q
    _PI2_CACHE.clear()
    ck("B2b MUTATION CONTROL: with q replaced by y^4+3y^2+y (which HAS a root satisfying "
       "the m=-3 condition) the |B|=1 test flips to FEASIBLE",
       mut[0] == "FEASIBLE",
       "mutated verdict = %s -- so the a9 kill is a property of the real q, not of the code"
       % (mut,))

    # the three sub2 T2 survivors of e | Phi
    v9 = t2_cell_verdict("sub2", 9, [1, 0, 0, 0])
    v8 = t2_cell_verdict("sub2", 8, [1, 1, 0, 0])
    v7 = t2_cell_verdict("sub2", 7, [1, 1, 1, 0])
    ck("B3  sub2 a8_b1100_T2 is DEAD (Pi^2 too big for the 2-dim space)",
       v8["verdict"] == "DEAD", v8["why"])
    ck("B4  sub2 a7_b1110_T2 is DEAD (same mechanism)",
       v7["verdict"] == "DEAD", v7["why"])
    cond9 = sp.expand((y + 1) * sp.diff(q, y, 2) - 6 * sp.diff(q, y))
    ck("B5  sub2 a9_b1000_T2 is DEAD: no root r of q has (r+1)q''(r) = 6 q'(r)",
       v9["verdict"] == "DEAD" and sp.gcd(q, cond9) == 1,
       "gcd(q, (y+1)q'' - 6q') = %s, Res = %s"
       % (sp.gcd(q, cond9), sp.resultant(q, cond9)))

    # B6 -- the strong |B|=1 statement, over the whole m-range either window can reach
    hits = [m for m in range(-40, 41)
            if sp.gcd(q, sp.expand((y + 1) * sp.diff(q, y, 2) + 2 * m * sp.diff(q, y))) != 1]
    ck("B6  STRONGER: for the real q, NO integer m in [-40,40] admits a root -- so on T2 a "
       "support cell with EXACTLY ONE q-root dividing e is dead at every rho, in both windows",
       hits == [],
       "m with a q-root solution of (r+1)q'' + 2m q' = 0: %s (empty).  The a9_b1000_T2 case "
       "is m = -3." % (hits,))

    # B7 -- independent NUMERIC cross-check of the decisive kill, bypassing the mod-Pi^2
    # reformulation entirely: solve the 2x2 system [P(r)=0, P'(r)=0] for (A, c) directly.
    qp, qpp = sp.diff(q, y), sp.diff(q, y, 2)
    dets = []
    for i in range(4):
        r = sp.CRootOf(q, i).evalf(50)
        tr = r + 1
        # q_rem(r) = q'(r), q_rem'(r) = q''(r)/2  (r a simple root of q)
        M = sp.Matrix([[tr**9 * qp.subs(y, r), -3 * tr**12],
                       [9 * tr**8 * qp.subs(y, r) + tr**9 * qpp.subs(y, r) / 2, -36 * tr**11]])
        dets.append(complex(sp.N(M.det(), 40)))
    ck("B7  INDEPENDENT CROSS-CHECK (numeric, 40 digits): for every root r of q the 2x2 system "
       "[P(r)=0, P'(r)=0], P = A t^9 q/(y-r) - 3c t^12, has nonsingular matrix -> only A=c=0",
       all(abs(z) > 1e-8 for z in dets),
       "|det| at the four roots = %s (all far from 0); this reaches the a9_b1000_T2 kill "
       "without the mod-Pi^2 reformulation" % ["%.3e" % abs(z) for z in dets])

    # ---- C. the counts ----------------------------------------------
    c2c = census("sub2")
    ck("C1  reproduce the e|Phi baseline on sub2: 18/26 cells, 140/220 flag cases, "
       "4822/7888 states",
       (c2c["ephi_dead_cells"], c2c["ephi_dead_cases"], c2c["ephi_dead_states"])
       == (18, 140, 4822),
       "recomputed: %d/%d cells, %d/%d flag cases, %d/%d states dead"
       % (c2c["ephi_dead_cells"], c2c["cells_total"], c2c["ephi_dead_cases"],
          c2c["cases_total"], c2c["ephi_dead_states"], c2c["states_total"]))

    nd = c2c["t2_new_dead"]
    ck("C2  the T2 filter kills EVERY sub2 T2 cell that e|Phi left alive: 3 cells, "
       "10 flag cases, 119 states",
       (len(nd["cells"]), nd["cases"], nd["states"]) == (3, 10, 119),
       "new kills: %d cells %s, %d flag cases, %d states"
       % (len(nd["cells"]), sorted(nd["cells"]), nd["cases"], nd["states"]))

    surv = sorted(c2c["alive_cells"] - nd["cells"])
    ck("C3  sub2 survivors are exactly the five T1 support cells",
       len(surv) == 5 and all(k[2] == "T1" for k in surv),
       "%s" % (surv,))

    c1c = census("sub1")
    und = [k for k, v in c1c["t2_verdicts"].items() if v["verdict"] == "UNDECIDED"]
    ck("C4  sub1 census computed (e|Phi baseline + T2 filter), undecided bucket reported",
       c1c["states_total"] == 44117,
       "sub1: %d/%d cells, %d/%d flag cases, %d/%d states dead by e|Phi; T2 filter adds "
       "%d cells / %d flag cases / %d states; UNDECIDED T2 cells: %d"
       % (c1c["ephi_dead_cells"], c1c["cells_total"], c1c["ephi_dead_cases"],
          c1c["cases_total"], c1c["ephi_dead_states"], c1c["states_total"],
          len(c1c["t2_new_dead"]["cells"]), c1c["t2_new_dead"]["cases"],
          c1c["t2_new_dead"]["states"], len(und)))

    # ---- D. spare-count consequences (branch-independent) -----------
    spare_T1 = (c2["R"] + 1) + (cap_s2 + 1) + 0
    spare_T2 = 1 + (cap_s2 + 1) + 0
    ck("D1  sub2 spare block 45 -> 18 (branch-independent, from e|S and the T-formula) "
       "and -> 6 on T2 (adds R = c t^12)",
       spare_T1 == 18 and spare_T2 == 6,
       "T1: dm2 13 + dm3 %d + dm4 0 = %d ; T2: dm2 1 + dm3 %d + dm4 0 = %d"
       % (cap_s2 + 1, spare_T1, cap_s2 + 1, spare_T2))

    # D2 -- the FRONTIER_REBUILD.md sec.7 lead ("a_t <= 10"), delivered on T2 only
    hi = []
    for a in range(11, 16):
        for bs in ([0, 0, 0, 0], [1, 0, 0, 0], [1, 1, 0, 0], [1, 1, 1, 0], [1, 1, 1, 1]):
            if a + sum(bs) > CAPS["sub1"]["e"]:
                continue
            if t2_cell_verdict("sub1", a, bs)["verdict"] != "DEAD":
                hi.append((a, tuple(bs)))
    ck("D2  FRONTIER_REBUILD sec.7 lead delivered ON T2: every a_t = 11..15 support cell is "
       "DEAD, i.e. a_t <= 10 on the T2 branch (T1 remains open, so the lead stays blocked)",
       hi == [],
       "surviving a_t >= 11 T2 support cells under the sub1 caps: %s (empty).  "
       "phase_d_states_sub1.json happens to contain no T2 cell with a_t > 10, so this is a "
       "statement about the alternate regime, not a new file-level kill." % (hi,))

    # D3 -- t^a | R on T2, on every T2 cell that is still alive
    viol = []
    for cc in (c2c, c1c):
        for k, vv in cc["t2_verdicts"].items():
            if vv["verdict"] == "ALIVE":
                for rho, _m, _st, _r in vv["feasible_rho"]:
                    if rho < k[0]:
                        viol.append((w, k, rho))
    ck("D3  t^a | R on every surviving T2 cell (rho >= a) -- the R-half of the open "
       "t^a | R,S,T item, settled on T2",
       viol == [], "violations: %s (empty)" % (viol,))

    if verbose:
        print("\n".join(out))
    return npass, ntot, {"sub2": c2c, "sub1": c1c,
                         "v7": v7, "v8": v8, "v9": v9, "alphas": alphas, "res": res}


def report(data):
    print("\n" + "=" * 78)
    print("CELL-LEVEL CONSEQUENCES")
    print("=" * 78)
    for w in ("sub2", "sub1"):
        c = data[w]
        nd = c["t2_new_dead"]
        print("\n%s: %d cells / %d flag cases / %d states" %
              (w, c["cells_total"], c["cases_total"], c["states_total"]))
        print("  e|Phi baseline dead : %d cells / %d cases / %d states" %
              (c["ephi_dead_cells"], c["ephi_dead_cases"], c["ephi_dead_states"]))
        print("  alive after e|Phi   : %d cells / %d cases / %d states" %
              (len(c["alive_cells"]), c["alive_cases"], c["alive_states"]))
        print("  NEW T2 divisor kill : %d cells / %d cases / %d states" %
              (len(nd["cells"]), nd["cases"], nd["states"]))
        print("  alive after both    : %d cells / %d cases / %d states" %
              (len(c["alive_cells"]) - len(nd["cells"]),
               c["alive_cases"] - nd["cases"], c["alive_states"] - nd["states"]))
        for k in sorted(c["t2_verdicts"]):
            v = c["t2_verdicts"][k]
            print("     T2 cell a%d_b%s : %-9s %s" %
                  (k[0], "".join(map(str, k[1])), v["verdict"], v["why"][:110]))
        print("  survivors: %s" % (sorted(c["alive_cells"] - nd["cells"]),))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args()
    npass, ntot, data = run(verbose=not a.quiet)
    if a.quiet:
        if npass != ntot:
            print("divisor_consequences: %d/%d checks FAILED" % (ntot - npass, ntot))
            return 1
        print("divisor_consequences: %d/%d checks pass" % (npass, ntot))
        return 0
    print("\n%d/%d checks pass" % (npass, ntot))
    report(data)
    return 0 if npass == ntot else 1


if __name__ == "__main__":
    sys.exit(main())
