#!/usr/bin/env python3
"""tpower_divisibility.py -- verify the `t^a | dm2,dm3,dm4` reduction, and the
cell-level DIVISIBILITY KILL it comes from, from the exact G-system syzygy.

THE SYZYGY (verified here, exactly, against `face_kill_sweep.canonical_G_generators`,
which itself asserts the canonical normalisation `G5 = G5body + Phi`):

    2*(G5 + d2*G3 + d1*G2 + d0*G1)  ==  2*Phi - dm1*(d2*dm1^2 + 3*dm1*dm3 + 3*dm2^2)

so on ANY common zero of G1,G2,G3,G5 -- i.e. on any genuine G-system lift --

    2*Phi = e*(d2*e^2 + 3*e*S + 3*R^2),      e = dm1, R = dm2, S = dm3, T = dm4.

Since `Phi = -(1/6630)*(y+1)^30*q(y)` with `q` the fixed quartic, this is a
DIVISIBILITY constraint `e | 2*Phi` in `Q[y]`, and it is checkable by hand.

CONSEQUENCES, in increasing depth (each verified by a function below):

  [A]  q is SQUAREFREE and coprime to (y+1), so every divisor of `(y+1)^30*q`
       has multiplicity <= 1 at each root of q.  With `v_t(e) = a_t` and
       `v_{r_j}(e) = b_j` the flag coordinates, this forces  b_j in {0,1}.

  [B]  deg e = 10 EXACTLY on sub2, hence  a_t + sum(b_j) = 10.
       (>= 10 from the identity's degree count; <= 10 from the sub2 window cap
       for `dm1`, k = 5, cap 2k = 10 -- `full_system_bridge.WEIGHT`/`STRIP_DEGCAP`.)

  [C]  On the T2 branch (d1 = 0), with t = y+1, the identity TOGETHER WITH the
       dm4-free relation `H3 := dm1*G3|_{dm4 from G1}` forces
            t^a | dm2,   t^a | dm3,   t^a | dm4,
       so substituting dm_i = t^a * dm_i_bar and dividing exactly collapses the
       stripped spare ansatz from 45 coefficients to 45 - 3a.

[A] and [B] are cell-level and need NO Groebner basis: they kill flag cases
outright.  [C] is a reduction, not a kill.

The [C] proof is a valuation case analysis, and it is done here by EXHAUSTIVE
ENUMERATION rather than by prose, so that no case is waved through.  Orders are
taken at t = y+1:  a = v(e), rho = v(R), s = v(S), and d2/d0 contribute unknown
orders delta2, delta0 >= 0 which are enumerated too.  A configuration survives
only if BOTH relations can hold:

  * the identity `d2*e^3 + 3*e^2*S + 3*e*R^2 = 2*Phi` with v(2*Phi) = 30:
    if the minimum order is attained ONCE it must equal 30; otherwise it must
    be <= 30 (cancellation);
  * `H3 = -6*d0*e^2*R - 6*d2*e*R*S - e^4 - 6*R*S^2 = 0`:
    the minimum order must be attained at least TWICE (a unique minimum is a
    nonzero lowest term, so H3 != 0).

Usage:  python tpower_divisibility.py            # everything
        python tpower_divisibility.py --cell a9_b1000_T2
"""
from __future__ import annotations

import argparse
import itertools
import json
import os

import sympy as sp

import face_kill_sweep as fks
import full_system_bridge as fsb

HERE = os.path.dirname(os.path.abspath(__file__)) or "."
INF = 10 ** 9                      # order of the zero polynomial

Y = sp.Symbol("y")
Q_POLY = 2048 * Y ** 4 - 512 * Y ** 3 + 320 * Y ** 2 - 240 * Y + 195


# =====================================================================
#  0.  The syzygy itself
# =====================================================================
def check_syzygy(verbose=True):
    d0, d1, d2, dm1, dm2, dm3, dm4, Phi = fks._gsystem_symbols()
    G = {k: sp.expand(v[0]) for k, v in fks.canonical_G_generators().items()}
    lhs = sp.expand(2 * (G["G5"] + d2 * G["G3"] + d1 * G["G2"] + d0 * G["G1"]))
    rhs = sp.expand(2 * Phi - dm1 * (d2 * dm1 ** 2 + 3 * dm1 * dm3 + 3 * dm2 ** 2))
    resid = sp.expand(lhs - rhs)
    if verbose:
        print("[0] SYZYGY  2*(G5 + d2*G3 + d1*G2 + d0*G1)")
        print("            == 2*Phi - dm1*(d2*dm1^2 + 3*dm1*dm3 + 3*dm2^2)")
        print("    residual =", resid, "  -> EXACT" if resid == 0 else "  *** NONZERO ***")
        print("    G5 Phi-coefficient =", sp.expand(G["G5"]).coeff(Phi), "(canonical form requires 1)")
    return resid == 0


def check_phi(verbose=True):
    """Phi = -(1/6630)*(y+1)^30*q, q squarefree, q(-1) != 0."""
    from bigrade_annotator import _phi_stripped
    ph = sp.expand(_phi_stripped())
    ratio = sp.cancel(ph / sp.expand((Y + 1) ** 30 * Q_POLY))
    sqfree = sp.gcd(sp.Poly(Q_POLY, Y), sp.Poly(sp.diff(Q_POLY, Y), Y)).degree() == 0
    coprime = Q_POLY.subs(Y, -1) != 0
    irred = len(sp.factor_list(Q_POLY)[1]) == 1 and sp.factor_list(Q_POLY)[1][0][1] == 1
    if verbose:
        print("\n[A] Phi / ((y+1)^30 * q) =", ratio, " (expected -1/6630)")
        print("    deg Phi =", sp.Poly(ph, Y).degree(), " (sub2 cap for Phi: k=17 -> 2k = 34)")
        print("    q squarefree:", sqfree, "  q irreducible/Q:", irred,
              "  q(-1) =", Q_POLY.subs(Y, -1), "!= 0 ->", coprime)
        print("    => every divisor of 2*Phi has multiplicity <= 1 at each root of q")
        print("    => b_j = v_{r_j}(e) in {0,1}   [CELL-LEVEL KILL]")
    return ratio == sp.Rational(-1, 6630) and sqfree and coprime


def check_deg_e(verbose=True):
    """deg e = 10 exactly on sub2."""
    caps = fsb.STRIP_DEGCAP["sub2"]
    # window cap for dm1: WEIGHT[dm1] = 12k -> k, stripped cap 2k on sub2
    k_dm1 = fsb.WEIGHT["dm1"] // 12
    cap_e = 2 * k_dm1
    k_phi = fsb.WEIGHT["Phi"] // 12
    deg_phi = 2 * k_phi
    lo = None
    for E in range(0, 40):
        # deg of e*(d2*e^2 + 3*e*S + 3*R^2)  <=  E + max(4+2E, E+cap_S, 2*cap_R)
        ub = E + max(4 + 2 * E, E + caps["dm3"], 2 * caps["dm2"])
        if ub >= deg_phi and lo is None:
            lo = E
    if verbose:
        print("\n[B] sub2 caps: deg d2<=4, deg R<=%d, deg S<=%d;  deg Phi = %d"
              % (caps["dm2"], caps["dm3"], deg_phi))
        print("    smallest E with  E + max(4+2E, E+%d, %d) >= %d  is E = %s"
              % (caps["dm3"], 2 * caps["dm2"], deg_phi, lo))
        print("    window cap for dm1 (k=%d, sub2 -> 2k): deg e <= %d" % (k_dm1, cap_e))
        print("    => deg e = %d EXACTLY  =>  a_t + sum(b_j) = %d   [CELL-LEVEL KILL]"
              % (cap_e, cap_e))
    return lo == cap_e


# =====================================================================
#  1.  The cell-level census
# =====================================================================
def census(universe="phase_d_states_sub2.json", verbose=True):
    data = json.load(open(os.path.join(HERE, universe), encoding="utf-8"))
    cases = data["cases"]
    bad_b, bad_deg, alive = [], [], []
    for c in cases:
        b = tuple(c["b"])
        if any(x > 1 for x in b):
            bad_b.append(c)
        elif c["a_t"] + sum(b) != 10:
            bad_deg.append(c)
        else:
            alive.append(c)
    ns = lambda L: sum(c["state_count"] for c in L)
    if verbose:
        print("\n[CENSUS] %s: %d cells, %d states" % (universe, len(cases), ns(cases)))
        print("    killed by  b_j <= 1        : %3d cells  %5d states"
              % (len(bad_b), ns(bad_b)))
        print("    killed by  a + sum(b) = 10 : %3d cells  %5d states"
              % (len(bad_deg), ns(bad_deg)))
        print("    SURVIVING                  : %3d cells  %5d states"
              % (len(alive), ns(alive)))
        cols = sorted({(c["a_t"], tuple(c["b"]), c["branch"]) for c in alive})
        print("    surviving (a_t, b, branch) columns:")
        for a, b, br in cols:
            n = sum(c["state_count"] for c in alive
                    if c["a_t"] == a and tuple(c["b"]) == b and c["branch"] == br)
            print("        a%-2d b%s %-3s  %d cells  %d states"
                  % (a, "".join(str(x) for x in b), br,
                     len([c for c in alive if c["a_t"] == a
                          and tuple(c["b"]) == b and c["branch"] == br]), n))
    return bad_b, bad_deg, alive


# =====================================================================
#  2.  [C] the t^a divisibility, by exhaustive valuation enumeration
# =====================================================================
def _identity_ok(a, rho, s, delta2, v_phi=30):
    """`d2*e^3 + 3*e^2*S + 3*e*R^2 = 2*Phi` possible at these orders?"""
    orders = [delta2 + 3 * a, 2 * a + s, a + 2 * rho]
    orders = [o for o in orders if o < INF]
    if not orders:
        return v_phi >= INF
    m = min(orders)
    if orders.count(m) == 1:
        return m == v_phi
    return m <= v_phi


def _h3_ok(a, rho, s, delta0, delta2):
    """`-6*d0*e^2*R - 6*d2*e*R*S - e^4 - 6*R*S^2 = 0` possible at these orders?

    A unique minimum is a surviving lowest-order term, so H3 != 0.  (e != 0, so
    the e^4 term is always present: the all-zero degenerate case cannot occur.)
    """
    orders = [delta0 + 2 * a + rho, delta2 + a + rho + s, 4 * a, rho + 2 * s]
    orders = [o for o in orders if o < INF]
    m = min(orders)
    return orders.count(m) >= 2


def verify_tpower(a, cap_R=12, cap_S=14, dmax=60, verbose=True):
    """Exhaustively refute `v_t(R) < a`, then `v_t(S) < a`, then conclude for T."""
    # -- R = 0 is impossible: H3 degenerates to -e^4 = 0, and e != 0.
    # -- so 0 <= rho <= cap_R.
    bad_rho = []
    for rho in range(0, a):
        for s in list(range(0, cap_S + 1)) + [INF]:
            for delta2 in range(0, dmax):
                if not _identity_ok(a, rho, s, delta2):
                    continue
                for delta0 in range(0, dmax):
                    if _h3_ok(a, rho, s, delta0, delta2):
                        bad_rho.append((rho, s, delta2, delta0))
                        break
                if bad_rho:
                    break
            if bad_rho:
                break
        if bad_rho:
            break
    rho_ok = not bad_rho

    # -- given rho >= a, refute v_t(S) < a from the identity alone
    bad_s = []
    for s in range(0, a):
        for rho in range(a, cap_R + 1):
            for delta2 in range(0, dmax):
                if _identity_ok(a, rho, s, delta2):
                    bad_s.append((s, rho, delta2))
                    break
            if bad_s:
                break
        if bad_s:
            break
    s_ok = not bad_s

    # -- G1 (d1=0):  3*d2*e*R + 3*e*T + 3*R*S = 0  =>  T = -(d2*R + R*S/e)
    #    v(T) >= min(delta2 + rho, rho + s - a) >= min(a, a + a - a) = a
    t_ok = rho_ok and s_ok

    if verbose:
        print("\n[C] t^a divisibility for a = %d  (T2 branch, d1 = 0; caps R<=%d, S<=%d)"
              % (a, cap_R, cap_S))
        print("    v_t(dm2) >= a : %s%s"
              % (rho_ok, "" if rho_ok else "   *** surviving config %s ***" % bad_rho[:1]))
        print("    v_t(dm3) >= a : %s%s"
              % (s_ok, "" if s_ok else "   *** surviving config %s ***" % bad_s[:1]))
        print("    v_t(dm4) >= a : %s   (from G1: T = -(d2*R + R*S/e), "
              "v >= min(v(d2)+rho, rho+s-a) >= a)" % t_ok)
    return rho_ok, s_ok, t_ok


def reduced_count(a, window="sub2", verbose=True):
    caps = fsb.STRIP_DEGCAP[window]
    full = sum(caps[n] + 1 for n in ("dm2", "dm3", "dm4"))
    per = {n: caps[n] - a + 1 for n in ("dm2", "dm3", "dm4")}
    red = sum(max(v, 0) for v in per.values())
    if verbose:
        print("\n[REDUCTION] window %s, a = %d" % (window, a))
        print("    full spare ansatz : %s  = %d coefficients"
              % ({n: caps[n] + 1 for n in ("dm2", "dm3", "dm4")}, full))
        print("    after dm_i = t^%d * bar : %s  = %d coefficients" % (a, per, red))
        print("    %d - 3a = %d - %d = %d   -> matches: %s"
              % (full, full, 3 * a, full - 3 * a, red == full - 3 * a))
    return full, red


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", type=int, default=9)
    args = ap.parse_args()
    print("=" * 78)
    print("G-SYSTEM DIVISIBILITY CONSTRAINT  e | 2*Phi   and the t^a reduction")
    print("=" * 78)
    ok = check_syzygy()
    ok &= check_phi()
    ok &= check_deg_e()
    census()
    for a in (7, 8, 9):
        verify_tpower(a)
    reduced_count(args.a)
    print("\nsyzygy + Phi factorisation + deg-e forcing all verified:", ok)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
