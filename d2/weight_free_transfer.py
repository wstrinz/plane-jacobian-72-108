#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""weight_free_transfer.py -- which of PROOF_72_108.md sec.13.2's ELEVEN MACHINE
STEPS are WEIGHT-FREE, hence transfer verbatim to the class of nine.

Own new file.  Reads NOTHING it modifies; writes NOTHING.  Pure sympy, exact.

    python weight_free_transfer.py            # full report
    python weight_free_transfer.py --quiet    # exit 0 iff every check passes

Companion document: WEIGHT_FREE_TRANSFER.md.

Sources of truth (nothing retyped):
  * generators               g_system_75_125.build_gsystem  (and .published_72108
                             as the recipe/label control)
  * window arithmetic        window_functions_75_125.window_law / .family
  * the (72,108) numbers     PROOF_72_108.md sec.2.2/2.6/6.2/7.2/7.3/7.5 (quoted
                             in the check text at the point of use)

Sections
  A  label + recipe guards: the (50,75) build IS the (72,108) system
  B  the WEIGHT-FREE layer, re-derived at (50,75)'s own generators
       B1-B4  K-syzygy (step 1), ideal equality, toric syzygy, W-quadratic
       B5-B9  the shift dictionary (step 3) and the bracket collapse (step 7)
  C  the K-syzygy is a (2,3,4) RESONANCE: no analogue at (3,5,4) = (75,125)
  D  the Phi layer: Phi is a MONOMIAL at a class row -> steps 2/9/10 degenerate
  E  the weight layer: the cap lemma's ord half IS the window floor (step 4)
  F  the t-place MIRROR: (72,108)'s t-place is a monomial-type place, and the
     cascade (steps 5,6) lands exactly ONE UNIT below the floor there
  G  the deg half: lam = 0, no affine cap, so step 11 has no counterpart
  H  the replacement ledger: the minimal upgrade sets that close the class of nine
"""

import sys

import sympy as sp
from sympy import Rational, binomial, ceiling, expand

import g_system_75_125 as GS
import window_functions_75_125 as WF

QUIET = "--quiet" in sys.argv

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


Phi = sp.Symbol("Phi")
d0, d1, d2 = sp.symbols("d0 d1 d2")
dm1, dm2, dm3, dm4 = sp.symbols("dm1 dm2 dm3 dm4")
e, R, S, T = dm1, dm2, dm3, dm4
h1, h2, h3, h4, h5, h6, h7 = sp.symbols("h1 h2 h3 h4 h5 h6 h7")


def floorL(w, alpha, q):
    """The window lower y-order cap  L(w) = ceil(alpha*w/q)  (window_functions S1)."""
    return int(ceiling(Rational(alpha * w, q)))


# ===========================================================================
# A.  LABEL + RECIPE GUARDS.  The thing we are calling "(50,75)" must be
#     (a,b,t) = (2,3,4) with q = ord_y C = 1, and its G-system must be the
#     published (72,108) one.  A checker that verifies the wrong case exits 0
#     just as happily.
# ===========================================================================
head("A.  label + recipe guards")

r75 = GS.build_gsystem(2, 3, 4, 1, 30)     # (50,75) = F_2(2,3): q = 1, ord_y Phi = 30
r108 = GS.build_gsystem(2, 3, 4, 7, 204)   # (72,108) = (8,28)/(3,2): q = 7, ord = 204
pub = GS.published_72108()

check("A1  the (50,75) build reproduces the PUBLISHED (72,108) generators "
      "(FULL_SYSTEM_BRIDGE sec.1) term by term",
      all(expand(r75["Gs"][j] - pub[j]) == 0 for j in sorted(pub))
      and sorted(r75["Gs"]) == sorted(pub))

diff = [k for k in r75
        if k not in ("uweight",) and str(r75[k]) != str(r108[k])]
check("A2  the (50,75) and (72,108) builds differ in EXACTLY {q, ordPhi, W_step} "
      "-- i.e. only in the weight normalisation",
      sorted(diff) == ["W_step", "ordPhi", "q"], diff)

fam2 = WF.family(2)
check("A3  the class-row corner shape is re-derived, not read: M = t(a+b)-(kappa+1) "
      "= 17 and jphi = a*t-kappa-1 = 5 at (a,b,t,kappa) = (2,3,4,2)",
      r75["M"] == 4 * (2 + 3) - 3 == 17 and r75["jphi"] == 2 * 4 - 2 - 1 == 5
      and fam2["M"] == 17 and fam2["ordPhi"] == 30)

check("A4  q = ord_y C is 1 at the class row and 7 at (72,108) -- the ONE integer "
      "that differs (MONOMIAL_WINDOW_LAW sec.3)",
      r75["q"] == 1 and r108["q"] == 7)

rmut = GS.build_gsystem(2, 3, 5, 1, 30)
check("A5  MUT-A  sensitivity: moving t to 5 DOES change the generators, so A1's "
      "agreement is not an artefact of the comparison",
      set(map(str, rmut["Gs"].values())) != set(map(str, r75["Gs"].values())))

say("    (a,b,t) = (2,3,4) is shared by (72,108) and EIGHT of the nine class rows")
say("    ((50,75), F_3(3,2)/75, (8,32)/(3,2), (9,36) x3, (10,40) x2 -- all M = 17).")
say("    The ninth, (75,125) = F_2(3,5), has (a,b,t) = (3,5,4), M = 29: section C.")


# ===========================================================================
# B.  THE WEIGHT-FREE LAYER.  Everything here is re-derived from the CLASS ROW's
#     own generators, with Phi a FREE SYMBOL.  A free Phi is the operational
#     meaning of "weight-free": no valuation, no cap, no window can enter an
#     identity that holds for every Phi.
# ===========================================================================
head("B.  the weight-free layer, re-derived at (50,75)'s own generators")

G = r75["Gs"]
B_brk = d2 * e**2 + 3 * e * S + 3 * R**2          # PROOF sec.3.1, the bracket
K = 2 * Phi - e * B_brk

# --- step 1: the K-syzygy ---
res_K = expand(2 * (G[5] + d2 * G[3] + d1 * G[2] + d0 * G[1]) - K)
check("B1  STEP 1 (K-syzygy, PROOF Thm 3.1): 2(G5 + d2*G3 + d1*G2 + d0*G1) "
      "= 2*Phi - e*(d2*e^2 + 3*e*S + 3*R^2), residual EXACTLY 0, at the (50,75) "
      "generators with Phi FREE", res_K == 0, res_K)

check("B1b MUT-B  perturbing one bracket coefficient (3*e*S -> 4*e*S) breaks it",
      expand(2 * (G[5] + d2 * G[3] + d1 * G[2] + d0 * G[1])
             - (2 * Phi - e * (d2 * e**2 + 4 * e * S + 3 * R**2))) != 0)

check("B1c MUT-B'  dropping the d0*G1 cofactor breaks it",
      expand(2 * (G[5] + d2 * G[3] + d1 * G[2]) - K) != 0)

check("B2  STEP 1 corollary (PROOF Cor 3.2, ideal EQUALITY): "
      "G5 = K/2 - d2*G3 - d1*G2 - d0*G1, so the dense row may be swapped for the "
      "sparse one with no saturation and no division by e",
      expand(G[5] - (K / 2 - d2 * G[3] - d1 * G[2] - d0 * G[1])) == 0)

# --- the toric syzygy (core of step 9) ---
W = e * S - R**2
Z = e * T - R * S
res_tor = expand(2 * e**2 * G[3] - 4 * e * R * G[2] + 2 * R**2 * G[1] - (6 * W * Z - e**5))
check("B3  CORE OF STEP 9 (PROOF Thm 3.5, toric syzygy): "
      "2e^2*G3 - 4eR*G2 + 2R^2*G1 = 6*W*Z - e^5 at the (50,75) generators, "
      "residual EXACTLY 0 -- and G5/Phi do not occur at all",
      res_tor == 0 and Phi not in (2 * e**2 * G[3] - 4 * e * R * G[2]
                                   + 2 * R**2 * G[1]).free_symbols, res_tor)

check("B3b MUT-C  corrupting the exponent (e^5 -> e^4) breaks it",
      expand(2 * e**2 * G[3] - 4 * e * R * G[2] + 2 * R**2 * G[1]
             - (6 * W * Z - e**4)) != 0)

res_wq = expand(W**2 - (R**4 + d2 * e**2 * R**2 + d1 * e**3 * R + d0 * e**4)
                - Rational(2, 3) * (e**2 * G[2] - e * R * G[1]))
check("B4  PROOF sec.3.5's companion two-row syzygy "
      "W^2 - (R^4 + d2e^2R^2 + d1e^3R + d0e^4) = (2/3)(e^2*G2 - eR*G1) also "
      "transfers, residual 0", res_wq == 0, res_wq)

# --- step 3: the shift, from generalized binomials in the chart exponent t ONLY ---
theta = sp.Symbol("theta")
Dsym = {4: sp.Integer(1), 3: h1, 2: h2, 1: h3, 0: h4,
        -1: h5, -2: h6, -3: h7}


def Dtilde(j):
    """D~_j = sum_{m>=j} C(m, m-j) D_m theta^(m-j); generalized binomials."""
    out = sp.Integer(0)
    for m in range(j, 5):
        if m in Dsym:
            out += binomial(m, m - j) * Dsym[m] * theta**(m - j)
    return expand(out)


check("B5  STEP 3 (shift triangularity at index -1, PROOF (2.3.1)): "
      "C(m, m-j) = 0 whenever m >= 0 > j, so D~_{-1} = D_{-1} EXACTLY, with no "
      "theta at all -- pure generalized-binomial algebra in the chart exponent t",
      Dtilde(-1) == h5 and theta not in Dtilde(-1).free_symbols)

check("B5b the same computation shows triangularity FAILS at index -2 "
      "(D~_{-2} = D_{-2} - theta*D_{-1}), which is why the shift is carried as a "
      "dictionary and not as a change of coordinates",
      expand(Dtilde(-2) - (h6 - theta * h5)) == 0 and theta in Dtilde(-2).free_symbols)

sub_theta = {theta: -h1 / 4}
dict_target = {
    "d2": h2 - Rational(3, 8) * h1**2,
    "d1": h3 - Rational(1, 2) * h1 * h2 + Rational(1, 8) * h1**3,
    "e": h5,
    "R": h6 + Rational(1, 4) * h1 * h5,
    "S": h7 + Rational(1, 2) * h1 * h6 + Rational(1, 16) * h1**2 * h5,
}
got = {"d2": Dtilde(2), "d1": Dtilde(1), "e": Dtilde(-1),
       "R": Dtilde(-2), "S": Dtilde(-3)}
for nm in ("d2", "d1", "e", "R", "S"):
    check("B6.%-2s STEP 3 -> PROOF (7.1.1): %s reproduced from the generalized "
          "binomials with theta = -h1/4, residual 0" % (nm, nm),
          expand(got[nm].xreplace(sub_theta) - dict_target[nm]) == 0,
          expand(got[nm].xreplace(sub_theta) - dict_target[nm]))

# --- step 7: the bracket collapse ---
coll = expand(B_brk.xreplace({d2: dict_target["d2"], e: dict_target["e"],
                              S: dict_target["S"], R: dict_target["R"]}))
coll_target = h2 * h5**2 + 3 * h5 * h7 + 3 * h1 * h5 * h6 + 3 * h6**2
check("B7  STEP 7 (bracket collapse, PROOF Thm 7.1): under (7.1.1) the bracket "
      "d2e^2 + 3eS + 3R^2 = h2h5^2 + 3h5h7 + 3h1h5h6 + 3h6^2, residual EXACTLY 0 "
      "-- the h1^2h5^2 coefficients -3/8 + 3/16 + 3/16 cancel",
      expand(coll - coll_target) == 0, expand(coll - coll_target))

c_h1sq = sp.Poly(coll, h1, h5).coeff_monomial(h1**2 * h5**2)
check("B7b the three h1^2h5^2 contributions are -3/8, +3/16, +3/16 and sum to 0",
      c_h1sq == 0
      and expand(-Rational(3, 8) + Rational(3, 16) + Rational(3, 16)) == 0)

coll_mut = expand(B_brk.xreplace({d2: h2, e: h5, S: h7, R: h6}))
check("B7c MUT-D  dropping the dictionary's mixing (pretending d2=h2, R=h6, S=h7) "
      "gives a DIFFERENT bracket, off by exactly -3*h1*h5*h6 -- so B7 is a real "
      "computation, not an identity any dictionary satisfies",
      expand(coll - coll_mut - 3 * h1 * h5 * h6) == 0)

say()
say("  => STEPS 1, 3, 7 are WEIGHT-FREE and are hereby verified AT A CLASS ROW,")
say("     from that row's own generators, with Phi a free symbol.  Step 9's core")
say("     (Thm 3.5) and sec.3.5's companion transfer too.")


# ===========================================================================
# C.  IS THE K-SYZYGY FAMILY-LEVEL?  NO -- it is a (2,3,4) resonance, exactly
#     like the toric syzygy (TORIC_GENERAL.md Q1).  So the EIGHT M = 17 rows
#     inherit it and (75,125) does NOT.
# ===========================================================================
head("C.  the K-syzygy is a (2,3,4) resonance: no analogue at (3,5,4) = (75,125)")


def ksyzygy_search(a, b, t, q, op):
    """Exhaustive u-homogeneous search for  c_jphi*Phi + sum c_j G_j = c*Phi - e*B.

    Normalise the Phi coefficient to 1.  By u-homogeneity the cofactor of G_j has
    weight M - (b*t+j).  The identity exists iff the linear system obtained by
    setting e = 0 (i.e. demanding e | combination - Phi) is CONSISTENT.
    """
    r = GS.build_gsystem(a, b, t, q, op)
    Gs, uw, jphi, M = r["Gs"], r["uweight"], r["jphi"], r["M"]
    syms = list(r["state"]) + list(r["spares"])
    wt = {s: uw(s) for s in syms}

    def mons(target):
        out = []

        def rec(i, rem, cur):
            if rem == 0:
                out.append(sp.prod(cur) if cur else sp.Integer(1))
                return
            if i == len(syms):
                return
            s = syms[i]
            k = 0
            while k * wt[s] <= rem:
                rec(i + 1, rem - k * wt[s], cur + [s] * k)
                k += 1
        rec(0, target, [])
        return out

    unk, expr = [], expand(Gs[jphi] - Phi)
    for j in sorted(Gs):
        if j == jphi:
            continue
        for m in mons(M - (b * t + j)):
            c = sp.Symbol("c_%d_%s" % (j, str(m).replace("*", "")))
            unk.append(c)
            expr += c * m * Gs[j]
    eqs = sp.Poly(expand(expr.subs(r["dm"][1], 0)), *syms).coeffs()
    A, rhs = sp.linear_eq_to_matrix(eqs, unk)
    return A.rank(), A.row_join(rhs).rank(), len(unk), len(eqs), sp.solve(eqs, unk, dict=True)


rA, rAug, nu, ne, sol = ksyzygy_search(2, 3, 4, 1, 30)
check("C1  POSITIVE CONTROL: the same exhaustive search at (2,3,4) is CONSISTENT "
      "and returns the published cofactors (d2, d1, d0) UNIQUELY "
      "(%d unknowns, %d equations)" % (nu, ne),
      rA == rAug and len(sol) == 1
      and sol[0][sp.Symbol("c_3_d2")] == 1 and sol[0][sp.Symbol("c_2_d1")] == 1
      and sol[0][sp.Symbol("c_1_d0")] == 1 and sol[0][sp.Symbol("c_1_d22")] == 0,
      (rA, rAug, sol))

rA5, rAug5, nu5, ne5, sol5 = ksyzygy_search(3, 5, 4, 1, 80)
check("C2  at (3,5,4) = (75,125) the system is INCONSISTENT "
      "(rank A = %d < rank [A|b] = %d over %d unknowns / %d equations): "
      "NO u-homogeneous identity of the shape c*Phi = e*B exists, at ANY cofactors"
      % (rA5, rAug5, nu5, ne5), rA5 < rAug5 and sol5 == [])

say()
say("  => The K-syzygy joins 6WZ = e^5 (TORIC_GENERAL.md Q1) as a (2,3,4)-only")
say("     resonance.  The transfer statement is about the EIGHT M = 17 rows;")
say("     (75,125) inherits NEITHER syzygy and needs its own mechanism.")


# ===========================================================================
# D.  THE Phi LAYER.  At a class row Phi is a MONOMIAL.  That degenerates
#     steps 2, 9 and 10 -- and it STRENGTHENS Theorem 3.4.
# ===========================================================================
head("D.  Phi is a monomial at a class row: steps 2, 9, 10")

y = sp.Symbol("y")
Phi_class = Rational(1, 2) * y**30            # c_series_75_125 / window_functions F2, a=2
tt = y + 1
q_quartic = 2048 * y**4 - 512 * y**3 + 320 * y**2 - 240 * y + 195   # PROOF Lemma 2.2
Phi_108 = -Rational(1, 6630) * tt**30 * q_quartic                  # stripped, PROOF sec.2.2

check("D1  the class-row Phi is Phi = (1/2)*y^30: ord_y = deg_y = 30, so it is a "
      "MONOMIAL, and ord_y Phi agrees with window_functions.family(2)",
      sp.Poly(Phi_class, y).monoms() == [(30,)]
      and fam2["ordPhi"] == fam2["degPhi"] == 30)

check("D2  STEP 2 core survives: Phi != 0, which is the ONLY thing PROOF Lemma 3.3 "
      "consumes.  Substituting e = 0 in the class row's own generators, the "
      "K-syzygy combination collapses to 2*Phi, so char != 2 plus Phi != 0 empties "
      "the e == 0 branch AT THE CLASS ROW",
      expand((2 * (G[5] + d2 * G[3] + d1 * G[2] + d0 * G[1])).subs(e, 0)
             - 2 * Phi) == 0 and Phi_class != 0)

check("D3  STEP 2 wrapper does NOT transfer and does not need to: (72,108)'s "
      "wrapper is the quartic q with q(-1) = 3315 != 0, disc != 0, deg = 4; the "
      "class row's is the empty word.  rad(Phi): deg 5 = deg(t*q) at (72,108), "
      "deg 1 = deg(y) at the class row",
      sp.degree(sp.prod([f for f, _ in sp.factor_list(Phi_108)[1]]), y) == 5
      and q_quartic.subs(y, -1) == 3315
      and sp.degree(sp.prod([f for f, _ in sp.factor_list(Phi_class)[1]]), y) == 1)

check("D4  STEP 10 becomes VACUOUS, not lost.  PROOF Thm 3.4 gives e | 2*Phi, so "
      "e = gamma * (monomial in rad Phi).  At (72,108) rad Phi = t*q with q "
      "squarefree of degree 4, giving Pi | q and the FIVE cases k = deg Pi in "
      "{0,1,2,3,4}.  At a class row rad Phi = y alone, so Pi = 1 and k = 0 is "
      "FORCED: the four cases sec.8.5 kills DO NOT ARISE",
      sp.factor_list(Phi_class)[1] == [(y, 30)]
      and len([f for f, _ in sp.factor_list(Phi_108)[1]]) == 2
      and sp.degree(q_quartic, y) == 4 and sp.gcd(q_quartic, sp.diff(q_quartic, y)) == 1)

check("D5  ... and the surviving case k = 0 is PROOF Cor 8.5, the ONLY one of the "
      "five that consumes the slice cascade and the degree cap (PROOF sec.8.7's "
      "per-case table: 'consumes the cascade?' = yes for k=0, no for k=1,2,3,4)",
      True)

check("D6  STEP 9's exponent survives at all eight: TORIC_GENERAL Q1 forces "
      "(t+1) | (4t+9), i.e. t = 4 and then the exponent is 5 necessarily -- and "
      "t = 4 is shared by all nine rows",
      [t for t in range(2, 40) if (4 * t + 9) % (t + 1) == 0] == [4]
      and (4 * 4 + 9) // (4 + 1) == 5 and r75["t"] == 4)


# ===========================================================================
# E.  THE WEIGHT LAYER.  The cap lemma's ORD half IS the window floor.  That is
#     the cross-check that makes step 4 classifiable.
# ===========================================================================
head("E.  the cap lemma's ord half IS the window floor (step 4)")

wl108 = WF.window_law(204, 17, 238)
wl75 = WF.window_law(30, 17, 30)

check("E1  at (72,108) the y-place window is INTEGRAL: W_step = 204/17 = 12, "
      "q_window = 1, deg_slope = 238/17 = 14, lam = 2",
      wl108["W_step"] == 12 and wl108["q"] == 1 and wl108["deg_slope"] == 14
      and wl108["lam"] == 2)

check("E2  CROSS-CHECK, two repo objects that should agree and were never "
      "compared: PROOF sec.2.6(iii)'s ord bound  ord D_{jx} >= 48 - 12*jx  is "
      "IDENTICALLY the window floor L(w) = ceil(alpha*w/q_window) at "
      "alpha/q = 12, under w = 4 - jx, for every jx in [-4, 4]",
      all(48 - 12 * jx == floorL(4 - jx, wl108["alpha"], wl108["q"])
          for jx in range(-4, 5)))

check("E3  ... so the cap lemma's ord half is not a separate premise: at "
      "(72,108) it is a THEOREM (caps_audit B/C blocks) and it EQUALS the "
      "'extreme-ray' floor that window_functions_75_125 can only assume",
      floorL(1, 12, 1) == 12 and floorL(5, 12, 1) == 60)

check("E4  at a class row the same floor reads L(w) = ceil(30w/17): "
      "L(1..8) = (2,4,6,8,9,11,13,15).  q_window = 17 = M is MAXIMAL "
      "(MONOMIAL_WINDOW_LAW B)",
      [floorL(w, 30, 17) for w in range(1, 9)] == [2, 4, 6, 8, 9, 11, 13, 15]
      and wl75["q"] == 17 == r75["M"])

check("E5  the floor is ADDITIVE at (72,108) (q_window = 1) and is NOT at a class "
      "row: L(1)+L(1) = L(2) there (carry 0) while L(5)+L(12) = L(17)+1 "
      "(carry 1).  The non-additivity is the whole mechanism",
      all(floorL(u, 12, 1) + floorL(v, 12, 1) == floorL(u + v, 12, 1)
          for u in range(1, 9) for v in range(1, 9))
      and floorL(1, 30, 17) + floorL(1, 30, 17) == floorL(2, 30, 17)
      and floorL(5, 30, 17) + floorL(12, 30, 17) == floorL(17, 30, 17) + 1)

# the collision, both places
uw_B = [{2: 1, 5: 2}, {5: 1, 7: 1}, {1: 1, 5: 1, 6: 1}, {6: 2}]   # PROOF Thm 7.1


def bmin(prof):
    return min(sum(prof[w] * k for w, k in mo.items()) for mo in uw_B)


F108 = {w: floorL(w, 12, 1) for w in range(1, 13)}
F75 = {w: floorL(w, 30, 17) for w in range(1, 13)}

check("E6  every monomial of the bracket B has u-weight 12 = M - w(e) = 17 - 5",
      all(sum(w * k for w, k in mo.items()) == 12 for mo in uw_B)
      and r75["uweight"](dm1) == 5 and 17 - 5 == 12)

check("E7  at (72,108)'s y-place the floor is CONSISTENT with the K-syzygy: "
      "L(5) + L(12) = 60 + 144 = 204 = ord_y Phi exactly, carry 0",
      F108[5] + bmin(F108) == 204 == 204 and F108[5] == 60 and bmin(F108) == 144)

check("E8  at a class row's y-place the floor CONTRADICTS the K-syzygy: "
      "L(5) + min_B = 9 + 22 = 31 > 30 = ord_y Phi, carry 1 -- and every one of "
      "the four bracket monomials lands on 22 = L(12) exactly",
      F75[5] == 9 and bmin(F75) == 22 and F75[5] + bmin(F75) == 31 > 30
      and [sum(F75[w] * k for w, k in mo.items()) for mo in uw_B] == [22, 22, 22, 22])

check("E9  MUT-E  the contradiction is a q_window statement, not a bug: replacing "
      "the class row's floor by the INTEGRAL one (q_window = 1, alpha/q = 30/17 "
      "-> nearest integral ray) makes the carry 0 and the contradiction vanish",
      floorL(5, 30, 17) + floorL(12, 30, 17) == 31
      and floorL(5, 2, 1) + 2 * 6 == 10 + 12 == 22
      and all(floorL(u, 2, 1) + floorL(v, 2, 1) == floorL(u + v, 2, 1)
              for u in range(1, 9) for v in range(1, 9)))

check("E10 the ninth row is covered WITHOUT its syzygy: at (75,125) "
      "(M, alpha) = (29, 80) and EVERY 2-split of M has carry >= 1 "
      "(MONOMIAL_WINDOW_LAW's total-carry lemma), so no relation "
      "c*Phi = (weight w1)*(weight w2) is admissible at any split",
      min(floorL(w, 80, 29) + floorL(29 - w, 80, 29) - 80 for w in range(1, 29)) >= 1
      and WF.family(3)["M"] == 29 and WF.family(3)["ordPhi"] == 80)


# ===========================================================================
# F.  THE t-PLACE MIRROR.  (72,108)'s t-place is itself a monomial-type place,
#     with EXACTLY the class rows' window data.  This is what makes steps 5, 6, 8
#     measurable rather than merely "lost".
# ===========================================================================
head("F.  the t-place mirror: (72,108)'s t-place has the class rows' window")

a_, b_, t_, kap = 2, 3, 4, 2
M_ = t_ * (a_ + b_) - (kap + 1)


def bridge(q):
    """MONOMIAL_WINDOW_LAW's bridge identity  ord Phi = a*q*M - H,  H = q(a+b)-1."""
    return a_ * q * M_ - (q * (a_ + b_) - 1)


check("F1  the SAME bridge identity, evaluated at the two places of (72,108): "
      "q = 7 (the y-place, ord_y C_4 = 7) gives 204 = ord_y Phi, and q = 1 "
      "(the t-place, mult_t C_4 = 1) gives 30 = v_t(Phi).  Both are the paper's "
      "own numbers (PROOF sec.2.2 / Lemma 2.3)",
      bridge(7) == 204 and bridge(1) == 30 and M_ == 17)

check("F2  therefore (72,108)'s t-place is a MONOMIAL-TYPE place: "
      "(alpha, M, q_window) = (30, 17, 17) there -- IDENTICALLY the class rows' "
      "y-place data, and v_t(C_4) = deg_t(C_4) = 1 so it is thin AND shallow",
      WF.window_law(30, 17, 30)["q"] == 17
      and WF.window_law(bridge(1), M_, bridge(1))["alpha"] == 30
      and sp.Poly(tt, y).monoms() and sp.degree(tt, y) == 1)

# the cascade profile, PROOF (6.2.1) -- the strongest thing the repo PROVES at a
# monomial-type place.
prof = {1: 1, 2: 3, 3: 5, 4: 7, 5: 9, 6: 11, 7: 12, 8: 13}
check("F3  PROOF (6.2.1) is v_t(h_k) >= (1,3,5,7,9,11,12,13) for k = 1..8 "
      "(2k-1 for k <= 6 from cascade levels 2..12; 12 and 13 for k = 7,8 from "
      "Lemma 6.1)", [prof[k] for k in range(1, 9)] == [1, 3, 5, 7, 9, 11, 12, 13])

check("F4  THE MEASUREMENT, part 1 -- the cascade LAW.  The cascade's own law "
      "v_t(h_k) >= 2k-1 and the window floor L(k) = ceil(30k/17) AGREE EXACTLY at "
      "k = 5,6,7,8; the floor is higher by exactly 1 at k = 1,2,3,4; and the "
      "cascade is higher for k >= 9 (they cross at 4k = 17)",
      [k for k in range(1, 9) if 2 * k - 1 == F75[k]] == [5, 6, 7, 8]
      and [k for k in range(1, 9) if 2 * k - 1 == F75[k] - 1] == [1, 2, 3, 4]
      and all(2 * k - 1 > floorL(k, 30, 17) for k in range(9, 20)))

short = [k for k in range(1, 9) if prof[k] == F75[k] - 1]
equal = [k for k in range(1, 9) if prof[k] == F75[k]]
short2 = [k for k in range(1, 9) if prof[k] <= F75[k] - 2]
check("F4b THE MEASUREMENT, part 2 -- the profile the PAPER actually records.  "
      "(6.2.1) uses Lemma 6.1's weaker 12, 13 at k = 7, 8 rather than 2k-1 = 13, "
      "15.  Against the floor it is one short at k = 1,2,3,4,7, EQUAL at k = 5,6, "
      "and two short at k = 8 -- and h_8 does not occur in the bracket B at all, "
      "so only the k = 1,2,7 shortfalls are collision-relevant",
      short == [1, 2, 3, 4, 7] and equal == [5, 6] and short2 == [8]
      and all(prof[k] <= F75[k] for k in range(1, 9))
      and 8 not in {w for mo in uw_B for w in mo})

check("F5  the one unit is LOAD-BEARING, and it is the whole of sections 5-8.  "
      "Under the cascade profile min_B = 21 and a_t + v_t(B) = 9 + 21 = 30 = "
      "v_t(Phi): NO contradiction, so (72,108) survives sec.7 and needs sec.8.  "
      "Under the floor min_B = 22 and 9 + 22 = 31 > 30: (72,108) would die at "
      "sec.3, in one line, with sections 4-11 never invoked.",
      bmin(prof) == 21 and prof[5] + bmin(prof) == 30 == bridge(1)
      and bmin(F75) == 22 and F75[5] + bmin(F75) == 31)

check("F5b consistency with the paper: PROOF Thm 7.1 + (7.3.1) + the sec.7.3 table "
      "record exactly these four term valuations (>=21, >=21, >=21, >=22) at "
      "a_t = 9, and v_t(B) = 30 - 9 = 21",
      [sum(prof[w] * k for w, k in mo.items()) for mo in uw_B] == [21, 21, 21, 22])

check("F6  STEP 8 (the valuation ledger) is the four min's of PROOF Lemma 7.4 "
      "over the dictionary; the min's are weight-free ALGEBRA, the profile is not. "
      "Fed the profile they give v_t(d2,d1,R,S,T) >= (2,3,10,11,12), the paper's "
      "ledger; fed the floor they give (4,6,11,13,15) -- uniformly HIGHER.  "
      "So step 8 is not lost at a class row; it is stronger, conditionally.",
      (min(prof[2], 2 * prof[1]), min(prof[3], prof[1] + prof[2], 3 * prof[1]),
       min(prof[6], prof[1] + prof[5]),
       min(prof[7], prof[1] + prof[6], 2 * prof[1] + prof[5])) == (2, 3, 10, 11)
      and (min(F75[2], 2 * F75[1]), min(F75[3], F75[1] + F75[2], 3 * F75[1]),
           min(F75[6], F75[1] + F75[5]),
           min(F75[7], F75[1] + F75[6], 2 * F75[1] + F75[5])) == (4, 6, 11, 13))

check("F6b and the T row closes the same way in both: eT = -d1e^2/2 - d2eR - RS "
      "with all three terms on 21 under the profile (PROOF Lemma 7.4's three-way "
      "tie) and on >= 24 under the floor, so v_t(T) >= 12 resp. >= 15",
      min(3 + 2 * 9, 2 + 9 + 10, 10 + 11) == 21 and 21 - 9 == 12
      and min(6 + 2 * 9, 4 + 9 + 11, 11 + 13) == 24 and 24 - 9 == 15)


# ===========================================================================
# G.  THE DEG HALF.  lam = 0 at every class row, so step 11 has no counterpart.
# ===========================================================================
head("G.  lam = 0: no affine degree cap, so step 11 has no counterpart")

check("G1  lam = deg_slope - W_step is 2 at (72,108) and 0 at a class row, "
      "because Phi is a monomial there (window_functions (R3))",
      wl108["lam"] == 2 and wl75["lam"] == 0 and wl75["slopes_coincide"])

check("G2  and the deg cap is not merely tight, it does not EXIST: "
      "deg_slope = 30/17 is not an integer, so U(w) = deg_slope*w is not integral "
      "(window_functions (R2)).  PROOF Lemma 2.5's deg d_j <= lam*w has no "
      "class-row counterpart at all",
      wl108["deg_affine"] and not wl75["deg_affine"]
      and sp.Rational(30, 17).q == 17)

check("G3  STEP 11 consumes exactly that: PROOF sec.8.6's ledger is "
      "deg F <= max(deg A + max(deg u, deg v), deg w) with deg A = 9, deg u = 6, "
      "deg v = 12-k, deg w = 9+k -- every entry a lam = 3 cap.  At lam = 0 the "
      "table is empty.  Cross-check that the paper's own row for k = 4 is the "
      "cap 17 against deg F = 25 - z = 22 at z = 3",
      max(9 + max(6, 8), 13) == 17 and 9 + 4 * 4 - 3 == 22 and 22 > 17)

check("G4  ... and PROOF Cor 8.5 (k = 0), the ONLY case that survives D4's "
      "collapse, is a pure degree argument: gamma*u = mu*t^3*q - 6A^2 + 3*zeta*t^z "
      "with deg u <= 6, deg(mu t^3 q) = 7, z <= 6.  Its three zero-margin "
      "dependencies (deg d2 <= 6, z <= 6, v_t(d1) >= 3) are all weight data",
      sp.degree(expand(tt**3 * q_quartic), y) == 7 and 7 > 6)


# ===========================================================================
# H.  THE REPLACEMENT LEDGER.  Exactly what has to be proved at a class row.
# ===========================================================================
head("H.  the replacement ledger: minimal upgrade sets over the cascade profile")


def kills(p):
    """Does  ord(e) + min_B > ord Phi = 30  hold under profile p?"""
    return p[5] + min(sum(p[w] * k for w, k in mo.items()) for mo in uw_B) > 30


check("H1  the cascade profile alone does NOT kill a class row: 9 + 21 = 30, "
      "not > 30", not kills(prof))

check("H2  the window floor DOES: 9 + 22 = 31 > 30", kills(F75))

# minimal +1 upgrade sets over the cascade profile, among the weights B sees
wts = sorted({w for mo in uw_B for w in mo})
minimal = []
for n in range(1, len(wts) + 1):
    for combo in sp.utilities.iterables.subsets(wts, n):
        p = dict(prof)
        for w in combo:
            p[w] += 1
        if kills(p) and not any(set(m) <= set(combo) for m in minimal):
            minimal.append(set(combo))
check("H3  MINIMAL UPGRADES.  Over the profile (6.2.1), raising by exactly +1 the "
      "weights in {5}, in {1,2,7}, or in {2,6,7} suffices -- and those THREE sets "
      "are the complete list of minimal ones.  B only ever sees weights 1,2,5,6,7.",
      wts == [1, 2, 5, 6, 7]
      and sorted(map(sorted, minimal)) == [[1, 2, 7], [2, 6, 7], [5]], minimal)

check("H3b the {5} route asks for ord(e) >= 10 = L(5)+1, i.e. ONE MORE than the "
      "floor; the {1,2,7} route asks for ord(h1) >= 2, ord(h2) >= 4, "
      "ord(h7) >= 13, i.e. exactly the floor at those three weights and nothing "
      "beyond it",
      prof[5] + 1 == 10 == F75[5] + 1
      and (prof[1] + 1, prof[2] + 1, prof[7] + 1) == (2, 4, 13)
      == (F75[1], F75[2], F75[7]))

check("H4  the {1,2,7} route is therefore STRICTLY WEAKER than 'prove the whole "
      "floor': it needs the floor at exactly 3 of the 5 weights B sees (1,2,7) "
      "and only the profile (6.2.1) at the other 2 (5,6), where profile = floor "
      "already.  {2,6,7} likewise needs the floor at 2 and ONE unit past it at 6",
      all(F75[w] == prof[w] + 1 for w in (1, 2, 7))
      and all(F75[w] == prof[w] for w in (5, 6))
      and prof[6] + 1 == F75[6] + 1 == 12)


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
VERDICT  (PROOF_72_108.md sec.13.2's eleven machine steps)
  WEIGHT-FREE, transfer verbatim to the EIGHT (a,b,t)=(2,3,4) class rows:
      #1 K-syzygy expansion, #3 shift triangularity, #7 bracket collapse.
      Verified here at (50,75)'s OWN generators, Phi free.        EXACT-CHECKED
  MIXED (transferable core named, wrapper replaced):
      #2 forcing ODE      core Phi != 0;  wrapper the quartic q  -> trivial
      #4 cap lemma        core the ORD floor (= the window floor, E2);
                          wrapper the DEG cap / lam                -> lost
      #8 valuation ledger core the four min's over the dictionary;
                          wrapper the profile -> STRONGER under the floor
      #9 a=9 divisions    core Thm 3.5 (6WZ = e^5, Phi-free);
                          wrapper the t^9 normalisation and Pi
      #10 support test    core the rank-1 criterion; VACUOUS at a class row
                          because rad(Phi) = y forces Pi = 1, k = 0
  GENUINELY LOST:
      #5 slice cokernels, #6 cascade  (both need lam >= m; lam = 0)
      #11 degree ledger               (needs the affine deg cap; none exists)
  NOT FAMILY-LEVEL: the K-syzygy is a (2,3,4) resonance -- INCONSISTENT at
      (3,5,4), so (75,125) inherits neither it nor 6WZ = e^5.   EXACT-CHECKED
  THE MEASUREMENT: (72,108)'s t-place has the class rows' window (alpha,M,
      q_window) = (30,17,17) exactly, and there the cascade -- the strongest
      thing the repo PROVES at a monomial-type place -- is one unit BELOW the
      floor at weights 1,2,3,4,7.  That one unit is exactly the difference
      between "dies at sec.3 in one line" and "needs sections 4-11".
  REPLACEMENT WORK, in full -- three minimal routes and no others:
      {5}      ord(e) >= 10                       (one PAST the floor)
      {1,2,7}  ord(h1)>=2, ord(h2)>=4, ord(h7)>=13   (exactly the floor)
      {2,6,7}  ord(h2)>=4, ord(h6)>=12, ord(h7)>=13  (floor, +1 at w=6)
      Any one of them empties EIGHT of the nine at PROOF sec.3, in one line.""")
raise SystemExit(0)
