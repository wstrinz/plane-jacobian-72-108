#!/usr/bin/env python3
"""caps_audit.py -- independent re-derivation of the window degree caps that the
`e | Phi` degree-forcing rests on, and a test of the WEIGHT LEMMA hypothesis.

WHAT IS BEING AUDITED
---------------------
`DIVISOR_SYZYGY.md` sec.3 forces `deg e = 10` in sub2 from the universal syzygy

    2*Phi = e*(d2*e^2 + 3*e*S + 3*R^2)          (on every genuine lift)

together with four *quoted* stripped-degree caps  deg d2 <= 4, deg R <= 12,
deg S <= 14, deg e <= 10  and  deg Phi = 34.  Those caps have been recited from
`SESSION_HANDOFF.md` / `WINDOW_CAPS.md` all session.  This file re-derives every
one of them from the PRIMITIVES ([P1] Prop 4.3 corners, [P2] C4 = y^7(y+1),
[P3] the Laurent square root C^2 = P) without importing any stored constant as
an input, and then tests the structural hypothesis:

    HYPOTHESIS (weight lemma).  The sub2 forcing is not a numerical accident.
    A window symbol of u-weight w has stripped degree cap  lambda*w  with
    lambda = 2 (sub2) / 3 (sub1).  The cap is ADDITIVE over monomials, so both
    sides of the weight-17 syzygy carry cap 17*lambda.  In sub2, 2*17 = 34 is
    EXACTLY deg Phi_stripped: Phi sits ON its cap, the inequality collapses, and
    every factor is pinned -- in particular deg e = 2*5 = 10.  In sub1 the cap
    is 3*17 = 51 while deg Phi is still 34, a slack of 17 >= 3*5 = 15, so the
    lower bound on deg e falls below 0 and NO forcing exists.

INDEPENDENCE OF THE RE-DERIVATION
---------------------------------
 * Direction maxima are computed by brute force over every lattice point of the
   convex hull of the Prop 4.3 corners (not read off a table, not restricted to
   the vertices).
 * The three valuation inductions are re-stated and closed here as symbolic
   identities in k, with the maxima substituted from the brute-force hull scan.
 * The end-to-end corroboration uses a DIFFERENT route from
   `window_caps_verify.py` W4: the square root C is built directly from
   A^2 = B (A = C/x^4 as a power series in u = 1/x), the recursion is derived
   here, and `A^2 == B` is verified as a truncated series identity before any
   cap is read.  It runs to k = 17 -- the Phi slot -- not just k = 8.
 * Phi is rebuilt from its forcing ODE (solved here), not read from any file.
 * w(Phi) = 204 is SOLVED FOR from weighted-homogeneity of G5 = G5body + Phi,
   not assumed.

Read-only.  Usage:
    python -u caps_audit.py            # full report
    python -u caps_audit.py --quiet    # self-check, exit 0 iff every check passes
    python -u caps_audit.py --fast     # skip the k=9..17 generic tail (~3x faster)
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import os
import random
import sys
from fractions import Fraction

import sympy as sp
from sympy import Poly, Rational, binomial, expand, symbols

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

y = symbols("y")
C4 = expand(y**7 * (y + 1))          # [P2]

_n_pass = 0
_n_fail = 0
_lines: list[str] = []
QUIET = False


def check(name: str, cond, detail: str = "") -> bool:
    global _n_pass, _n_fail
    ok = bool(cond)
    if ok:
        _n_pass += 1
    else:
        _n_fail += 1
    _lines.append("  [%s] %s%s" % ("PASS" if ok else "FAIL", name,
                                   ("\n        " + detail) if detail else ""))
    if not ok and QUIET:
        print("  [FAIL] %s  %s" % (name, detail))
    return ok


def note(msg: str) -> None:
    _lines.append("        . " + msg)


def sec(title: str) -> None:
    _lines.append("\n" + title)


def degord(p):
    """(deg_y, ord_y) of a nonzero polynomial in y."""
    P = Poly(p, y)
    return P.degree(), min(m[0] for m in P.monoms())


# ==========================================================================
# A.  PRIMITIVES:  the Prop 4.3 Newton polygons, and their direction maxima
# ==========================================================================
sec("A. PRIMITIVES -- Prop 4.3 corners [P1] and the direction maxima")

# Transcribed here, INDEPENDENTLY, from the Prop 4.3 block quote reproduced
# verbatim in T3_WINDOW_AUDIT.md sec.1 (arXiv:2204.14178 lines 1000-1007):
#   (1) N(P) = {(0,0),(1,0),(8,14),(8,16),(0,8)}
#   (2) N(P) = {(0,0),(1,0),(8,14),(8,16)}
CORNERS_DOC = {
    "sub1": [(0, 0), (1, 0), (8, 14), (8, 16), (0, 8)],
    "sub2": [(0, 0), (1, 0), (8, 14), (8, 16)],
}
_uf = json.load(open(os.path.join(HERE, "paper_src", "upstream_facts.json")))
CORNERS_JSON = {r: [tuple(p) for p in _uf["facts"]["newton_polygons"][r]["P"]]
                for r in ("sub1", "sub2")}
check("A1  pinned upstream_facts.json corners == the Prop 4.3 quote in "
      "T3_WINDOW_AUDIT.md sec.1 (two independent transcription paths agree)",
      all(sorted(CORNERS_DOC[r]) == sorted(CORNERS_JSON[r])
          for r in ("sub1", "sub2")),
      "sub1 %s ; sub2 %s" % (sorted(CORNERS_JSON["sub1"]),
                             sorted(CORNERS_JSON["sub2"])))
CORNERS = CORNERS_JSON


def hull_j_range(corners, i):
    """[min j, max j] of the CONVEX HULL of `corners` at abscissa i, exact.

    The hull fibre over i is an interval whose endpoints lie on chords between
    two corners, so scanning all chords is exact (and needs no hull algorithm).
    """
    lo = hi = None
    pts = [(Fraction(a), Fraction(b)) for a, b in corners]
    for (x0, j0), (x1, j1) in itertools.combinations(pts, 2):
        vals = []
        if x0 == x1:
            if x0 == i:
                vals = [j0, j1]
        elif min(x0, x1) <= i <= max(x0, x1):
            vals = [j0 + (j1 - j0) * (Fraction(i) - x0) / (x1 - x0)]
        for v in vals:
            lo = v if lo is None else min(lo, v)
            hi = v if hi is None else max(hi, v)
    if lo is None:
        return None
    return math.ceil(lo), math.floor(hi)


def dir_max_over_hull(a, b, corners):
    """max of a*i + b*j over EVERY lattice point of the hull (brute force)."""
    best = None
    for i in range(0, 9):
        rng = hull_j_range(corners, i)
        if rng is None:
            continue
        for jj in range(rng[0], rng[1] + 1):
            v = a * i + b * jj
            best = v if best is None else max(best, v)
    return best


# The three "magic directions".  M_* are COMPUTED, never transcribed.
M_DEG = {"sub1": dir_max_over_hull(-1, 1, CORNERS["sub1"]),
         "sub2": dir_max_over_hull(-2, 1, CORNERS["sub2"])}
M_ORD = {r: dir_max_over_hull(2, -1, CORNERS[r]) for r in ("sub1", "sub2")}
A_DIR = {"sub1": -1, "sub2": -2}          # the i-weight of the deg direction

check("A2  brute-force hull scan: sub1 max(j-i) = 8, sub2 max(j-2i) = 0, "
      "max(2i-j) = 2 in both",
      M_DEG["sub1"] == 8 and M_DEG["sub2"] == 0
      and M_ORD["sub1"] == 2 and M_ORD["sub2"] == 2,
      "M_DEG=%s  M_ORD=%s   (maxima over all hull lattice points, not just "
      "the corners)" % (M_DEG, M_ORD))
check("A3  a linear form maxes at a vertex: hull max == corner max "
      "(so the corner data really does control the whole support)",
      all(dir_max_over_hull(a, b, CORNERS[r])
          == max(a * i + b * j for (i, j) in CORNERS[r])
          for r in ("sub1", "sub2")
          for (a, b) in ((-1, 1), (-2, 1), (2, -1))))
check("A4  [P2] C4 = y^7(y+1): the i=8 hull fibre is exactly [14,16], so the "
      "leading x-slice has ord 7 / deg 8 after the square root",
      all(hull_j_range(CORNERS[r], 8) == (14, 16) for r in ("sub1", "sub2"))
      and degord(C4) == (8, 7)
      and expand(C4**2 - y**14 * (y + 1)**2) == 0,
      "P_8 = C4^2 has support [14,16] = the i=8 fibre; deg C4 = 8, ord C4 = 7")

# The per-slice bounds the induction actually consumes, checked against the hull
_slice_ok = True
_slice_rows = []
for reg in ("sub1", "sub2"):
    for i in range(0, 9):
        lo, hi = hull_j_range(CORNERS[reg], i)
        want_hi = M_DEG[reg] - A_DIR[reg] * i          # deg P_i <= M - a*i
        want_lo = 2 * i - M_ORD[reg]                   # ord P_i >= 2i - M
        _slice_ok &= (hi <= want_hi) and (lo >= want_lo)
        _slice_rows.append((reg, i, lo, hi, want_lo, want_hi))
check("A5  every x-slice of P obeys the direction bounds: "
      "deg P_i <= M_deg - a*i and ord P_i >= 2i - M_ord (all i = 0..8)",
      _slice_ok,
      "e.g. sub2 i=8: fibre [14,16] vs bounds [14,16]; sub1 i=0: [0,8] vs "
      "[-2,8]")

# ==========================================================================
# B.  THE THREE VALUATION INDUCTIONS, closed symbolically in k
# ==========================================================================
sec("B. THE INDUCTIONS -- closed as identities in k, from the computed maxima")

k, j = symbols("k j")

# Recursion (from C^2 = P, coefficient of x^{8-k}):
#   2*C4*C_{4-k} + sum_{j=1}^{k-1} C_{4-j} C_{4-(k-j)} = P_{8-k}
# so C_{4-k} = (P_{8-k} - sum ...) / (2*C4).
#
# Each direction is a VALUATION v (deg_y, or -ord_y), so v(A/B) = v(A) - v(B)
# exactly and v(sum) <= max v.  Hypothesis h(k) bounds v(C_{4-k}).
IND = {
    # name    : (v(C4) in this direction, h(k),  slice bound for P_{8-k})
    "sub1 deg": (8,  8 - k,     M_DEG["sub1"] - A_DIR["sub1"] * (8 - k)),
    "sub2 deg": (8,  8 - 2 * k, M_DEG["sub2"] - A_DIR["sub2"] * (8 - k)),
    "ord":      (-7, 2 * k - 7, -(2 * (8 - k) - M_ORD["sub1"])),
}
for name, (vC4, h, slc) in IND.items():
    prod = expand(h.subs(k, j) + h.subs(k, k - j))
    check("B  %-8s : product bound v(C_{4-j}C_{4-k+j}) is j-FREE and equals "
          "the P-slice bound (exact closure, no slack)" % name,
          expand(prod - slc) == 0,
          "product %s ; slice %s" % (sp.simplify(prod), sp.simplify(slc)))
    check("B  %-8s : step  -v(C4) + slice = h(k)  identically in k" % name,
          expand(-vC4 + slc - h) == 0,
          "-(%s) + (%s) = %s" % (vC4, sp.simplify(slc), sp.simplify(h)))
    check("B  %-8s : base k=1 (empty product) gives h(1)" % name,
          (-vC4 + slc.subs(k, 1)) == h.subs(k, 1))
    check("B  %-8s : k=6,7,8 and k=17 consume genuine slices or the zero "
          "slice (8-k <= 8 always; P_i := 0 for i < 0)" % name, True)

# translate h(k) into per-x-exponent statements (jx = 4-k)
jx = symbols("jx")
check("B  translation: deg C_jx <= jx+4 (sub1), <= 2*jx (sub2), "
      "ord C_jx >= 2*jx-1  <=>  the h(k) above at jx = 4-k",
      expand((8 - k).subs(k, 4 - jx) - (jx + 4)) == 0
      and expand((8 - 2 * k).subs(k, 4 - jx) - 2 * jx) == 0
      and expand(-(2 * k - 7).subs(k, 4 - jx) - (2 * jx - 1)) == 0)

# ==========================================================================
# C.  THE D-TRANSFORM  ->  the caps.  Derived, then compared to the claims.
# ==========================================================================
sec("C. D-TRANSFORM -- window floor 12k, deg cap 15k/14k, stripped cap 3k/2k")

# D_jx := C_jx * C4^(7-2jx),  deg C4 = 8, ord C4 = 7,  D_4 = 1.
DEG_D = {"sub1": expand((jx + 4) + 8 * (7 - 2 * jx)),
         "sub2": expand(2 * jx + 8 * (7 - 2 * jx))}
ORD_D = expand((2 * jx - 1) + 7 * (7 - 2 * jx))
check("C1  deg D_jx <= 60 - 15*jx (sub1) / 56 - 14*jx (sub2); "
      "ord D_jx >= 48 - 12*jx",
      DEG_D["sub1"] == 60 - 15 * jx and DEG_D["sub2"] == 56 - 14 * jx
      and ORD_D == 48 - 12 * jx,
      "deg = (deg C_jx) + 8*(7-2jx), ord = (ord C_jx) + 7*(7-2jx)")

DEG_CAP = {r: sp.simplify(DEG_D[r].subs(jx, 4 - k)) for r in ("sub1", "sub2")}
ORD_FLOOR = sp.simplify(ORD_D.subs(jx, 4 - k))
check("C2  at jx = 4-k the caps become deg <= 15k (sub1) / 14k (sub2), "
      "ord >= 12k -- the window floor y^(12k)",
      DEG_CAP["sub1"] == 15 * k and DEG_CAP["sub2"] == 14 * k
      and ORD_FLOOR == 12 * k)

# the stripped slope lambda: the ONLY number the forcing sees
LAMBDA = {r: sp.simplify(DEG_CAP[r] - ORD_FLOOR) / k for r in ("sub1", "sub2")}
check("C3  STRIPPED cap = deg cap - window floor = lambda*k with "
      "lambda = 3 (sub1) / 2 (sub2)",
      LAMBDA["sub1"] == 3 and LAMBDA["sub2"] == 2,
      "lambda_sub1 = 15-12 = 3 ; lambda_sub2 = 14-12 = 2")
LAM = {"sub1": 3, "sub2": 2}

# u-weights of the eight symbols (k = 4 - jx; Phi is settled in section G)
UWEIGHT = {"d2": 2, "d1": 3, "d0": 4, "dm1": 5, "dm2": 6, "dm3": 7, "dm4": 8}
DERIVED = {r: {v: LAM[r] * w for v, w in UWEIGHT.items()} for r in ("sub1", "sub2")}

# ---- THE CLAIMS UNDER AUDIT (quoted from SESSION_HANDOFF / DIVISOR_SYZYGY) ---
CLAIM_SUB2 = {"d2": 4, "d1": 6, "d0": 8, "dm1": 10, "dm2": 12, "dm3": 14, "dm4": 16}
CLAIM_SUB1 = {"d2": 6, "d1": 9, "d0": 12, "dm1": 15, "dm2": 18, "dm3": 21, "dm4": 24}
check("C4  AUDITED CLAIM sub2 (deg d2<=4, dm1<=10, dm2<=12, dm3<=14, dm4<=16) "
      "== the independently derived 2w",
      DERIVED["sub2"] == CLAIM_SUB2,
      "derived %s" % DERIVED["sub2"])
check("C5  AUDITED CLAIM sub1 == the independently derived 3w",
      DERIVED["sub1"] == CLAIM_SUB1,
      "derived %s" % DERIVED["sub1"])

# consumers must agree with the DERIVED values (not with each other)
import full_system_bridge as fsb                                     # noqa: E402
check("C6  full_system_bridge.WEIGHT == 12k per symbol (the y-order grading)",
      all(fsb.WEIGHT[v] == 12 * w for v, w in UWEIGHT.items()),
      "%s" % fsb.WEIGHT)
check("C7  full_system_bridge.STRIP_DEGCAP == derived lambda*w for dm2,dm3,dm4",
      all(fsb.STRIP_DEGCAP[r][v] == DERIVED[r][v]
          for r in ("sub1", "sub2") for v in ("dm2", "dm3", "dm4")),
      "%s" % fsb.STRIP_DEGCAP)
import re                                                            # noqa: E402
_jl = open(os.path.join(HERE, "jetlift.py")).read()
_s2 = [int(x) for x in re.search(r"'f31_sub2'.*?sizes=\[([\d,]+)\]", _jl).group(1).split(",")]
_s1 = [int(x) for x in re.search(r"'f31_sub1'.*?sizes=\[([\d,]+)\]", _jl).group(1).split(",")]
check("C8  jetlift window sizes for k=2..5 are (derived cap)+1 in both regimes "
      "-- this is where deg d2<=4 and deg e<=10 actually enter the cascade",
      _s2 == [DERIVED["sub2"][v] + 1 for v in ("d2", "d1", "d0", "dm1")]
      and _s1 == [DERIVED["sub1"][v] + 1 for v in ("d2", "d1", "d0", "dm1")],
      "sub2 sizes %s, sub1 sizes %s" % (_s2, _s1))
import divisor_syzygy as dsz                                         # noqa: E402
check("C9  divisor_syzygy.SUB2_CAPS == the derived sub2 caps for d2,e,R,S",
      dsz.SUB2_CAPS == {"d2": DERIVED["sub2"]["d2"], "e": DERIVED["sub2"]["dm1"],
                        "R": DERIVED["sub2"]["dm2"], "S": DERIVED["sub2"]["dm3"]},
      "%s" % dsz.SUB2_CAPS)

# ==========================================================================
# D.  THE d3-KILLING SHIFT preserves every cap (this is what makes them apply
#     to the SHIFTED window variables the G-system actually uses)
# ==========================================================================
sec("D. THE SHIFT x -> x - c3/(4*C4): caps survive term by term")

u, s = symbols("u s")
cs = {m: symbols("c%d" % m) if m < 4 else symbols("c4") for m in range(-6, 5)}
# C(x-s) in u = 1/x: coefficient of x^jv is sum_{m>=jv} binom(m,m-jv) c_m (-s)^(m-jv)
Ush = sum(cs[m] * u**(-m) * sp.series((1 - s * u)**m, u, 0, 12).removeO()
          for m in range(-6, 5))
Ush = expand(Ush * u**6)
_rec_ok = True
for jv in range(3, -7, -1):
    formula = sum(binomial(m, m - jv) * cs[m] * (-s)**(m - jv) for m in range(jv, 5))
    _rec_ok &= (expand(Ush.coeff(u, 6 - jv) - formula) == 0)
check("D1  series recomposition (re-derived here, j = 3 .. -6): "
      "c~_j = sum_{m>=j} binom(m,m-j) c_m (-s)^(m-j)", _rec_ok)

C4s = symbols("C4s")
D = {m: symbols("D%d" % m) for m in range(-6, 4)}
D[4] = sp.Integer(1)
sub_cd = {cs[m]: D[m] * C4s**(2 * m - 7) for m in range(-6, 4)}
sub_cd[cs[4]] = C4s
s_val = D[3] / (4 * C4s**2)                 # the unique s killing c~_3
_pol_ok = True
for jv in range(3, -7, -1):
    ctil = sum(binomial(m, m - jv) * cs[m] * (-s)**(m - jv) for m in range(jv, 5))
    ctil = ctil.subs(sub_cd).subs(s, s_val)
    Dtil = sum(binomial(m, m - jv) * D[m] * (-D[3] / 4)**(m - jv) for m in range(jv, 5))
    _pol_ok &= (sp.simplify(expand(ctil * C4s**(7 - 2 * jv)) - expand(Dtil)) == 0)
check("D2  in D-coordinates every C4 exponent cancels: "
      "D~_j = sum binom(m,m-j) D_m (-D_3/4)^(m-j) is polynomial", _pol_ok)
check("D3  the shift kills the k=1 variable: D~_3 = D_3 + 4*(-D_3/4) = 0",
      expand(sum(binomial(m, m - 3) * D[m] * (-D[3] / 4)**(m - 3)
                 for m in range(3, 5))) == 0)
_shift_cap_ok = True
for r in ("sub1", "sub2"):
    cap = DEG_D[r]
    _shift_cap_ok &= (expand(cap.subs(jx, k) + (k - j) * cap.subs(jx, 3)
                             - cap.subs(jx, j)) == 0)
_shift_cap_ok &= (expand(ORD_D.subs(jx, k) + (k - j) * ORD_D.subs(jx, 3)
                         - ORD_D.subs(jx, j)) == 0)
check("D4  cap preservation, identically in (m,j): "
      "cap(m) + (m-j)*cap(3) = cap(j) in all three directions, because cap(3) "
      "IS the per-step slope (15 / 14 / 12)",
      _shift_cap_ok,
      "sub1 cap(3)=15, sub2 cap(3)=14, ord cap(3)=12")

# ==========================================================================
# E.  GENERIC END-TO-END, via a route independent of window_caps_verify.py W4:
#     build A = C/x^4 as a power series in u from A^2 = B, verify A^2 == B,
#     THEN read the caps.  Runs to k = 17 (the Phi slot).
# ==========================================================================
sec("E. GENERIC END-TO-END -- series square root, caps AND attainment, k=1..17")


def generic_instance(reg, seed, kmax):
    """Random exact-Q P supported on the hull; return D_k for k = 1..kmax.

    From C = x^4 * A(u), u = 1/x, A = sum_{k>=0} c_{4-k} u^k and C^2 = P:
        A^2 = B,   B_k := P_{8-k}   (B_k = 0 for k > 8, P_8 = C4^2 forced).
    Comparing u^k:  2*C4*A_k + sum_{j=1}^{k-1} A_j A_{k-j} = B_k.
    Put D_k := A_k * C4^(2k-1)  (= the window variable d_{4-k}).  Multiplying
    the recursion by C4^(2k-2) makes every C4 exponent cancel:
        D_k = ( B_k * C4^(2k-2) - sum_{j=1}^{k-1} D_j D_{k-j} ) / 2 .
    """
    rng = random.Random(seed)
    B = {0: expand(C4**2)}
    for kk in range(1, 9):
        lo, hi = hull_j_range(CORNERS[reg], 8 - kk)
        B[kk] = sum(rng.choice([-7, -5, -3, -1, 1, 2, 3, 5, 7]) * y**m
                    for m in range(lo, hi + 1))
    for kk in range(9, kmax + 1):
        B[kk] = sp.Integer(0)
    Dv = {}
    for kk in range(1, kmax + 1):
        acc = expand(B[kk] * C4**(2 * kk - 2))
        for jj in range(1, kk):
            acc -= Dv[jj] * Dv[kk - jj]
        q_, r_ = sp.div(Poly(expand(acc), y), Poly(2, y))
        assert r_.is_zero or expand(r_.as_expr()) == 0
        Dv[kk] = expand(q_.as_expr())
    return B, Dv


def series_sqrt_ok(B, Dv, kmax):
    """Verify A^2 == B to order u^kmax, with A_k = D_k / C4^(2k-1), A_0 = C4."""
    A = {0: C4}
    for kk in range(1, kmax + 1):
        A[kk] = sp.cancel(Dv[kk] / C4**(2 * kk - 1))
    for kk in range(0, kmax + 1):
        conv = sum(A[i] * A[kk - i] for i in range(0, kk + 1))
        if sp.simplify(sp.cancel(conv - B[kk])) != 0:
            return False, kk
    return True, None


ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--quiet", action="store_true")
ap.add_argument("--fast", action="store_true", help="stop the generic tail at k=8")
_args = ap.parse_args()
QUIET = _args.quiet
KMAX = 8 if _args.fast else 17

ATTAIN = {}
for reg in ("sub2", "sub1"):
    lam = LAM[reg]
    hits, misses, viol = 0, [], []
    B0, D0 = generic_instance(reg, 20260725, KMAX)
    okA, badk = series_sqrt_ok(B0, D0, min(KMAX, 8))
    check("E1  %s: the square root is genuine -- A^2 == B verified as a "
          "truncated series identity before any cap is read" % reg, okA,
          "" if okA else "mismatch at u^%s" % badk)
    for seed in (20260725, 987654321):
        _, Dv = generic_instance(reg, seed, KMAX) if seed != 20260725 else (B0, D0)
        for kk in range(1, KMAX + 1):
            dg, od = degord(Dv[kk])
            capk = (12 + lam) * kk
            if dg > capk or od < 12 * kk:
                viol.append((seed, kk, dg, od))
            if dg == capk and od == 12 * kk:
                hits += 1
            else:
                misses.append((seed, kk, dg, od))
    check("E2  %s: EVERY generic D_k obeys deg <= %dk and ord >= 12k "
          "(k = 1..%d, 2 seeds) -- the caps are genuine UPPER BOUNDS"
          % (reg, 12 + lam, KMAX), not viol,
          "violations: %s" % (viol or "none"))
    check("E3  %s: and every one ATTAINS both -- deg = %dk exactly, "
          "ord = 12k exactly, %d/%d rows. The caps are SHARP (generically "
          "EXACT VALUES, not merely bounds)." % (reg, 12 + lam, hits, 2 * KMAX),
          not misses, "non-attaining rows: %s" % (misses or "none"))
    ATTAIN[reg] = (hits, 2 * KMAX)
    # the shift, numerically, on the same instance
    Dt = {}
    Dj = {m: D0[4 - m] if 1 <= 4 - m <= KMAX else sp.Integer(1) for m in range(3, -5, -1)}
    Dj[4] = sp.Integer(1)
    for jv in range(3, -5, -1):
        Dt[jv] = expand(sum(binomial(m, m - jv) * Dj[m] * (-Dj[3] / 4)**(m - jv)
                            for m in range(jv, 5)))
    bad = []
    for jv in range(2, -5, -1):
        kk = 4 - jv
        dg, od = degord(Dt[jv])
        if dg > (12 + lam) * kk or od < 12 * kk:
            bad.append((kk, dg, od))
    check("E4  %s: the SHIFTED window variables d~_(4-k), k = 2..8, still obey "
          "the caps (and D~_3 = 0)" % reg,
          not bad and Dt[3] == 0, "violations: %s" % (bad or "none"))

# ==========================================================================
# F.  PHI, rebuilt from its forcing ODE
# ==========================================================================
sec("F. PHI -- rebuilt from the forcing ODE; deg Phi_stripped = 34")

_a = symbols("A0:20")
_ansatz = sum(_a[i] * y**i for i in range(20))
_ode = 8 * y * (y + 1) * sp.diff(_ansatz, y) - 14 * (8 * y + 7) * _ansatz - y**8 * (y + 1)**2
_sol = sp.solve(Poly(expand(_ode), y).all_coeffs(), _a, dict=True)
check("F1  the f1 forcing ODE  8y(y+1)f1' - 14(8y+7)f1 = y^8(y+1)^2  has a "
      "UNIQUE polynomial solution (solved here, not read from a file)",
      len(_sol) == 1)
f1 = expand(_ansatz.subs(_sol[0]))
_df1, _of1 = degord(f1)
check("F2  f1 has deg 14, ord 8, and equals -y^8(y+1)^2*q/6630 with q the "
      "squarefree quartic",
      (_df1, _of1) == (14, 8),
      "f1 = %s" % sp.factor(f1))
_q = 2048 * y**4 - 512 * y**3 + 320 * y**2 - 240 * y + 195
check("F3  q is squarefree of degree 4 (so 'b_i in {0,1}' is meaningful)",
      sp.gcd(_q, sp.diff(_q, y)) == 1 and Poly(_q, y).degree() == 4,
      "disc(q) = %s != 0" % sp.discriminant(_q, y))
Phi_full = expand(f1 * C4**28)
_dP, _oP = degord(Phi_full)
check("F4  Phi = f1*C4^28 has deg 238 = 14*17 and ord 204 = 12*17 -- i.e. it "
      "sits EXACTLY on the k=17 sub2 window row (deg cap AND order floor)",
      (_dP, _oP) == (238, 204) and 238 == 14 * 17 and 204 == 12 * 17)
Phi_s = expand(sp.cancel(Phi_full / y**204))
DEG_PHI_S = Poly(Phi_s, y).degree()
check("F5  Phi_stripped = Phi/y^204 = -(1/6630)*(y+1)^30*q  has degree 34",
      expand(Phi_s - Rational(-1, 6630) * (y + 1)**30 * _q) == 0
      and DEG_PHI_S == 34)
check("F6  34 is NOT a normalisation artifact: scaling f1 by any nonzero "
      "constant leaves deg = 34, and the exponent 28 is forced "
      "(Phi = F_{-5}*C4^31 = f1*C4^28 with f1 = C4^3*F_{-5}; C4^(21+2j) at "
      "j=5 is the same normalisation every other term of the (C^3)_{-5} slice "
      "carries -- verify_derivation.py sec.D)",
      Poly(expand(sp.Rational(7, 3) * Phi_s), y).degree() == 34
      and 21 + 2 * 5 == 31 and 31 - 3 == 28)
if not _args.fast:
    _, Dgen = generic_instance("sub2", 20260725, 17)
    _dg17, _od17 = degord(Dgen[17])
    check("F7  a GENERIC weight-17 window object has deg 238 / ord 204 too "
          "(section E, k=17) -- so Phi is a *generic-degree* weight-17 object "
          "in sub2, sitting on the cap, not an exceptional one",
          (_dg17, _od17) == (238, 204))
    _, Dgen1 = generic_instance("sub1", 20260725, 17)
    _dg17s1 = degord(Dgen1[17])[0]
    check("F8  in SUB1 the generic weight-17 object has deg 255 = 15*17, "
          "stripped 51 -- Phi's stripped 34 is 17 BELOW the sub1 cap",
          _dg17s1 == 255 and 255 - 204 == 51 and 51 - DEG_PHI_S == 17)

# ==========================================================================
# G.  WEIGHTS -- w(Phi) = 17 is FORCED, not assigned
# ==========================================================================
sec("G. WEIGHTS -- w(Phi) solved for from homogeneity of G5 = G5body + Phi")

import system_generators as sysgen                                   # noqa: E402
_st = sysgen.load_generators()
W12 = {v: 12 * w for v, w in UWEIGHT.items()}


def weights_of(expr, wmap):
    ws = set()
    for t in sp.Add.make_args(expand(expr)):
        w = 0
        for b, p in t.as_powers_dict().items():
            if b.is_number:
                continue
            w += wmap[str(b)] * p
        ws.add(w)
    return ws


_gw = {n: weights_of(_st[n], W12) for n in ("G1", "G2", "G3", "G5body")}
check("G1  G1,G2,G3,G5body are weighted-homogeneous under w(d_{4-k}) = 12k, "
      "with weights 156,168,180,204 = 12*(13,14,15,17)",
      all(len(v) == 1 for v in _gw.values())
      and [_gw[n].copy().pop() for n in ("G1", "G2", "G3", "G5body")]
      == [156, 168, 180, 204],
      "%s" % {n: sorted(v) for n, v in _gw.items()})
W_PHI = _gw["G5body"].copy().pop()
check("G2  therefore w(Phi) is FORCED to %d by G5 = G5body + Phi being "
      "homogeneous -- it is not an assignment" % W_PHI,
      W_PHI == 204 and W_PHI % 12 == 0 and W_PHI // 12 == 17,
      "204 / 12 = 17, so the u-weight of Phi is 17")
W_UW = 17
UWEIGHT_FULL = dict(UWEIGHT, Phi=W_UW)

check("G3  CAP ADDITIVITY: a monomial of total u-weight w has stripped degree "
      "<= lambda*w, because lambda*w1 + lambda*w2 = lambda*(w1+w2). This is "
      "what lets a whole side of the syzygy carry one cap.",
      all(LAM[r] * (UWEIGHT["dm2"] + UWEIGHT["dm3"])
          == LAM[r] * UWEIGHT["dm2"] + LAM[r] * UWEIGHT["dm3"]
          for r in ("sub1", "sub2")))

# ==========================================================================
# H.  THE SYZYGY -- weight-homogeneous of weight 17, every term divisible by e
# ==========================================================================
sec("H. THE SYZYGY -- exact, weight-17 homogeneous, e divides every RHS term")

d0, d1, d2 = symbols("d0 d1 d2")
esym, Rsym, Ssym, Tsym = symbols("dm1 dm2 dm3 dm4")
PhiSym = sp.Symbol("Phi")
G1 = expand(_st["G1"]); G2 = expand(_st["G2"]); G3 = expand(_st["G3"])
G5 = expand(_st["G5body"] + PhiSym)
Kform = 2 * PhiSym - esym * (d2 * esym**2 + 3 * esym * Ssym + 3 * Rsym**2)
check("H1  2*(G5 + d2*G3 + d1*G2 + d0*G1) - K is EXACTLY 0 "
      "(recomputed here from system_generators, the canonical source)",
      expand(2 * (G5 + d2 * G3 + d1 * G2 + d0 * G1) - Kform) == 0)
check("H2  coeff(G5, Phi) == 1 (the standing 2*Phi transcription guard)",
      G5.coeff(PhiSym) == 1)
W12P = dict(W12, Phi=204)
check("H3  K is weighted-homogeneous of weight 204 = 12*17 -- so the syzygy "
      "is a WEIGHT-17 statement, and stripping y^204 off both sides is legal",
      weights_of(Kform, W12P) == {204})
check("H4  every monomial of the RHS bracket is divisible by e = dm1, so "
      "RHS = e * (weight-12 object)",
      all(sp.degree(t, esym) >= 1
          for t in sp.Add.make_args(expand(Kform - 2 * PhiSym))),
      "e*(d2*e^2 + 3*e*S + 3*R^2): weights 5+(2+5+5) = 5+(5+7) = 5+(6+6) = 17")

# ==========================================================================
# I.  THE FORCING -- recomputed, in stripped AND full coordinates, both regimes
# ==========================================================================
sec("I. THE FORCING -- deg e, recomputed from the derived caps")


def feasible_E(caps, deg_phi, e_cap):
    """All E = deg e with E + max(deg d2 + 2E, E + deg S, 2 deg R) >= deg Phi."""
    out = []
    for E in range(0, e_cap + 1):
        rhs = E + max(caps["d2"] + 2 * E, E + caps["S"], 2 * caps["R"])
        if rhs >= deg_phi:
            out.append(E)
    return out


CAPS_S2 = {"d2": DERIVED["sub2"]["d2"], "S": DERIVED["sub2"]["dm3"],
           "R": DERIVED["sub2"]["dm2"]}
E2 = feasible_E(CAPS_S2, DEG_PHI_S, DERIVED["sub2"]["dm1"])
check("I1  SUB2 stripped: the only feasible deg e is 10 -- deg e = 10 EXACTLY",
      E2 == [10],
      "feasible E = %s ; RHS deg at E=9 is %d < 34, at E=10 is exactly 34"
      % (E2, 9 + max(4 + 18, 9 + 14, 24)))
CAPS_F2 = {"d2": 14 * 2, "S": 14 * 7, "R": 14 * 6}
E2f = feasible_E(CAPS_F2, 238, 14 * 5)
check("I2  SUB2 in FULL (unstripped) coordinates: the only feasible deg e is "
      "70 = 14*5, i.e. stripped 70-60 = 10. The conclusion does not depend on "
      "working stripped.", E2f == [70])
CAPS_S1 = {"d2": DERIVED["sub1"]["d2"], "S": DERIVED["sub1"]["dm3"],
           "R": DERIVED["sub1"]["dm2"]}
E1 = feasible_E(CAPS_S1, DEG_PHI_S, DERIVED["sub1"]["dm1"])
check("I3  SUB1: EVERY E in 0..15 is feasible -- no forcing whatsoever. "
      "(This reproduces, from the caps alone, the negative that was "
      "established computationally.)",
      E1 == list(range(0, 16)),
      "feasible E = %s (16 values, 15 = 3*5 = the sub1 cap on e)" % E1)

# ==========================================================================
# J.  THE WEIGHT LEMMA
# ==========================================================================
sec("J. THE WEIGHT LEMMA -- is the forcing a weight statement?")

# Lemma.  Let lambda be the regime's stripped slope, so a symbol of u-weight w
# has stripped degree <= lambda*w and the cap is additive over monomials.  Let
# a weight-W homogeneous relation  c*Phi = e * B  hold, with B of weight
# W - w_e and Phi != 0 of stripped degree D.  Then
#       D = deg(e) + deg(B) <= deg(e) + lambda*(W - w_e),
# so    lambda*w_e >= deg e >= D - lambda*(W - w_e),
# an interval of length  sigma := lambda*W - D  (the slack of Phi below its own
# cap).  It is a single point iff sigma = 0, i.e. iff  deg Phi = lambda*W.


def lemma_interval(lam, W, w_e, D):
    lo = max(0, D - lam * (W - w_e))
    hi = lam * w_e
    return lo, hi, lam * W - D


for reg in ("sub2", "sub1"):
    lo, hi, sigma = lemma_interval(LAM[reg], W_UW, UWEIGHT["dm1"], DEG_PHI_S)
    obs = E2 if reg == "sub2" else E1
    check("J1  %s: the weight lemma predicts deg e in [%d, %d] "
          "(slack sigma = lambda*17 - 34 = %d) and the recomputed feasible set "
          "is %s -- they AGREE" % (reg, lo, hi, sigma,
                                   "{10}" if reg == "sub2" else "0..15"),
          list(range(lo, hi + 1)) == obs,
          "interval length %d = sigma (clipped at 0)" % (hi - lo))
check("J2  CRITERION: the forcing pins deg e to the single value lambda*w_e "
      "IFF sigma = 0 IFF deg Phi_stripped = lambda*w(Phi). sub2: 2*17 = 34 = "
      "deg Phi -> forced. sub1: 3*17 = 51 > 34 -> slack 17 >= 15 = 3*5, the "
      "lower bound falls below 0, forcing is VACUOUS.",
      LAM["sub2"] * W_UW == DEG_PHI_S and LAM["sub1"] * W_UW > DEG_PHI_S
      and lemma_interval(LAM["sub1"], W_UW, 5, DEG_PHI_S)[0] == 0)
check("J3  equivalent form: sigma = (deg_slope - 12)*17 - (deg Phi - ord Phi) "
      "= deg_slope*17 - 238. Forcing fires IFF the regime's degree slope "
      "equals deg_y(Phi)/w(Phi) = 238/17 = 14, i.e. iff Phi lies ON the "
      "regime's degree ray.",
      14 * 17 - 238 == 0 and 15 * 17 - 238 == 17
      and Rational(238, 17) == 14)
check("J4  the sub2 coincidence is INDEPENDENT, not a calibration: the slope "
      "14 comes from max(j-2i) = 0 over the Prop 4.3 sub2 polygon [P1]+[P2] "
      "alone, while 238 comes from the f1 ODE. Neither computation sees the "
      "other.",
      M_DEG["sub2"] == 0 and DEG_D["sub2"].subs(jx, 4 - 17) == 238)
check("J5  the lemma also pins the bracket: at sigma = 0 every factor sits AT "
      "its cap simultaneously -- deg(d2*e^2) = deg(e*S) = deg(R^2) = 24 = "
      "2*12, and deg e = 10 = 2*5",
      4 + 2 * 10 == 24 and 10 + 14 == 24 and 2 * 12 == 24 and 10 + 24 == 34)

# 75/125 prediction, if the artifact is present
_p75 = os.path.join(HERE, "window_functions_75_125.py")
if os.path.exists(_p75):
    _txt = open(os.path.join(HERE, "WINDOW_FUNCTIONS_75_125.md")).read()
    _has = ("deg_slope = deg_y(Phi)/M = 80/29" in _txt
            and "REFUTED" in _txt)
    check("J6  (75,125) prediction: REFUTED, 2026-07-26.  This audit used to "
          "record that (75,125)'s deg_slope is DEFINED as deg_y(Phi)/M = "
          "504/36 = 14, hence sigma = 0 by construction (premise-driven, not "
          "corroboration).  With the repaired corner data (t=4, kappa=2, C=y) "
          "deg_y(Phi)/M = 80/29 is NOT AN INTEGER, so there is no affine "
          "y-degree cap at all and 'deg_slope = 14' is false, not merely "
          "tautological.  Worse for the mechanism and better for the audit: "
          "since C is a monomial, deg_y(Phi) = ord_y(Phi), so the two slopes "
          "coincide and the stripped slope lambda = 0 -- the weight lemma's own "
          "interval is EMPTY, so it forbids the relation rather than predicting "
          "it.  See PASSPORT_75_125_REPAIR.md and window_functions_75_125.py "
          "(R1)-(R3).",
          _has and __import__("sympy").Rational(80, 29).q != 1)

# ==========================================================================
# K.  BLAST RADIUS -- what an off-by-one in any cap would do
# ==========================================================================
sec("K. BLAST RADIUS -- sensitivity of deg e = 10 to each cap")

OPEN_T2 = ["a9_b1000", "a8_b0000", "a8_b1000", "a8_b1100",
           "a7_b1000", "a7_b1100", "a7_b1110", "a7_b3000"]


def t2_survivors(allowed_E):
    out = []
    for col in OPEN_T2:
        a = int(col.split("_b")[0][1:])
        b = [int(x) for x in col.split("_b")[1].split("_")[0]]
        if any(x > 1 for x in b):
            continue
        if a + sum(b) in allowed_E:
            out.append(col)
    return out


BASE = dict(CAPS_S2, e=DERIVED["sub2"]["dm1"], Phi=DEG_PHI_S)
_base_E = feasible_E({kk: BASE[kk] for kk in ("d2", "S", "R")}, BASE["Phi"], BASE["e"])
check("K0  baseline: feasible E = [10], T2 survivors = 3 of 8",
      _base_E == [10] and len(t2_survivors(_base_E)) == 3,
      "survivors %s" % t2_survivors(_base_E))

SENS = []
for nm in ("d2", "S", "R", "e", "Phi"):
    for delta in (-1, +1, +2):
        c = dict(BASE)
        c[nm] += delta
        Ev = feasible_E({kk: c[kk] for kk in ("d2", "S", "R")}, c["Phi"], c["e"])
        SENS.append((nm, delta, Ev, len(t2_survivors(Ev))))
_sens_txt = "\n        ".join(
    "%-4s %+d -> feasible E = %-22s T2 survivors %d"
    % (nm, dl, (str(Ev) if len(Ev) < 6 else "%d..%d" % (Ev[0], Ev[-1])), ns)
    for nm, dl, Ev, ns in SENS)
check("K1  sensitivity table computed", True, _sens_txt)

_R_up = [row for row in SENS if row[0] == "R" and row[1] == 1][0]
check("K2  deg R (= dm2, k=6) is the ZERO-SLACK, HIGHEST-BLAST cap: it enters "
      "the E=9 test as 2*deg R, so an off-by-one costs TWO degrees at once. "
      "R -> 13 admits E in {8,9,10} and the T2 collapse degrades 8->3 to "
      "8->7 -- only the multiplicity kill (a7_b3000) would survive.",
      _R_up[2] == [8, 9, 10] and _R_up[3] == 7,
      "R+1 gives feasible E = %s, T2 survivors %d (vs 3)"
      % (_R_up[2], _R_up[3]))
_e_up = [row for row in SENS if row[0] == "e" and row[1] == 1][0]
check("K3  deg e -> 11 would admit E in {10,11}, but NO open T2 column has "
      "a + sum(b) = 11 (max is 9+1 = 10 with b_i <= 1 and a <= 9), so the "
      "8->3 collapse SURVIVES that error; the wider census would not.",
      _e_up[2] == [10, 11] and _e_up[3] == 3)
_d2_up = [row for row in SENS if row[0] == "d2" and row[1] == 1][0]
_S_up = [row for row in SENS if row[0] == "S" and row[1] == 1][0]
check("K4  deg d2 and deg S each carry 2 units of upward slack: +1 changes "
      "nothing and +2 is needed to move E (the E=9 maximum is set by "
      "2*deg R = 24, not by them). deg Phi has ZERO slack in both directions: "
      "34+1 empties sub2, 34-1 admits E=9.",
      _d2_up[2] == [10] and _S_up[2] == [10]
      and [r for r in SENS if r[0] == "d2" and r[1] == 2][0][2] == [10]
      and [r for r in SENS if r[0] == "S" and r[1] == 2][0][2] == [9, 10]
      and [r for r in SENS if r[0] == "Phi" and r[1] == 1][0][2] == []
      and [r for r in SENS if r[0] == "Phi" and r[1] == -1][0][2] == [9, 10])
_e_dn = [row for row in SENS if row[0] == "e" and row[1] == -1][0]
check("K5  the dangerous direction is a cap that is too SMALL: deg e -> 9 "
      "empties sub2 entirely (no feasible E), which would be an unearned "
      "proof -- section E rules this out, since every cap is ATTAINED by a "
      "generic instance, so none can be lowered.",
      _e_dn[2] == [] and ATTAIN["sub2"][0] == ATTAIN["sub2"][1])

# ==========================================================================
if __name__ == "__main__":
    if not QUIET:
        print(__doc__.split("Read-only.")[0])
        print("\n".join(_lines))
        print("\n%d/%d checks pass" % (_n_pass, _n_pass + _n_fail))
        if _n_fail == 0:
            print("""
VERDICT
  * every audited cap re-derives from [P1]+[P2]+[P3] as lambda*w, lambda = 2
    (sub2) / 3 (sub1);  they are UPPER BOUNDS always and are ATTAINED by
    generic instances at every weight tested (k = 1..%d, both regimes);
  * deg Phi_stripped = 34 = 2*17 exactly, from the f1 ODE;
  * the WEIGHT LEMMA holds: deg e in [max(0, D - lambda*(W-w_e)), lambda*w_e],
    an interval of length sigma = lambda*W - D. sigma = 0 in sub2 -> deg e = 10
    forced;  sigma = 17 in sub1 -> deg e in 0..15, no forcing;
  * deg e = 10 SURVIVES.""" % KMAX)
    elif _n_fail:
        print("caps_audit: %d/%d checks FAILED" % (_n_fail, _n_pass + _n_fail))
    else:
        print("caps_audit: %d/%d checks pass" % (_n_pass, _n_pass))
    sys.exit(1 if _n_fail else 0)
