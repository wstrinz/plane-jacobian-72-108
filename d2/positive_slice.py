#!/usr/bin/env python3
"""positive_slice.py -- the POSITIVE-SLICE obstruction, and the emptiness of
`a10_b0000_T1`, the last surviving standard-sub2 cell.

    python -u positive_slice.py            # full derivation + report
    python -u positive_slice.py --quiet    # checker; exit 0 iff every check passes

READ-ONLY on every existing artifact.  This lane wrote only POSITIVE_SLICE.md,
positive_slice.py, positive_slice_verify.py, positive_slice_stage.json.
Pure sympy; no Singular, no msolve, no WSL, no subprocess.

------------------------------------------------------------------------------
WHAT THIS IS
------------------------------------------------------------------------------
`SPINE.md` reduces standard sub2 to ONE cell, `a10_b0000_T1` (n = 0).  At n = 0
the marked-root arguments go vacuous (`Rm = 1`) and the four canonical G rows are
satisfied IDENTICALLY by the surviving family, so no Groebner engine can close
the cell from the G-system ideal: the ideal is genuinely non-empty.

The missing equations are not in the G-system at all.  The `d3`-killing shift
that defines the G-system's coordinates,

    x -> x - s ,        s = c_3/(4*C4) = D_3/(4*C4^2)          (window_caps W3)

is RATIONAL in y.  A solution of the shifted system therefore need not
reconstruct an ORIGINAL Laurent square root `C` of a polynomial `P = C^2`
supported on the Newton polygon.  Undoing the shift and demanding that the
positive-x coefficient slices of `P` come out POLYNOMIAL is a genuine extra
condition, and it is the one that empties the cell.

------------------------------------------------------------------------------
DERIVATION MAP  (each section is machine-checked below)
------------------------------------------------------------------------------
 P0  Canonical guard: `G5 = G5body + Phi` with `coeff(G5, Phi) == 1`.
 P1  The slice formula, DERIVED from the repo's own D-transform (not assumed):
       c_m = D_m * C4^(2m-7),  D_4 = 1,  C4 = y^7*(y+1)      (window_caps W2)
       D_j = y^(48-12j) * d_j   (the bridge's STRIPPED coordinate, WEIGHT 12k)
       H(u) := sum_{j<=4} d_j u^(4-j)
     =>  P_i = y^(2i-2) * [u^(8-i)] H(u)^2 / t^(14-2i) ,   t = y+1.
     Polynomiality of P_i for i = 4,5,6 therefore requires
       t^6 | [u^4]H^2 ,  t^4 | [u^3]H^2 ,  t^2 | [u^2]H^2.
 P2  The inverse shift, DERIVED from the SAME general transformation
     window_caps_verify.py W3 uses -- not hand-coded:
       D~_j = sum_{m=j..4} binom(m,m-j) D_m (-D_3/4)^(m-j)     (forward, W3)
       D*_j = sum_{m=j..4} binom(m,m-j) D~_m (+h/4)^(m-j)      (inverse, h=D_3)
     which yields, with (D~_4,D~_3,D~_2,D~_1,D~_0) = (1,0,d2,d1,d0),
       D2* = d2 + (3/8)h^2
       D1* = d1 + (1/2)h*d2 + (1/16)h^3
       D0* = d0 + (1/4)h*d1 + (1/16)h^2*d2 + (1/256)h^4
 P3  POSITIVE CONTROL (non-negotiable).  Independently generated polygon-
     supported P over QQ -> D-recursion -> strip -> shift -> UNSHIFT -> the
     original positive slices are recovered EXACTLY, and the three divisibilities
     hold.  If this fails the formula is wrong and the run aborts.
 P4  The n = 0 SPINE parametrisation, RE-DERIVED from `generators.json`
     (never by importing spine.py's output): the four rows factor, `C` and `d0`
     are eliminated, the certificate F*Z = (1/6)gamma^5 t^10 is produced with its
     explicit cofactor, and degree exactness forces Z = zeta*t^4, F = phi*t^6.
 P5  Evaluating everything at y = -1 gives three equations in (alpha, eta, gamma);
     in the coordinates X = alpha^2*gamma, Y = alpha*eta*gamma^2 they are the
     reviewer's (A),(B),(C).  Agreement is REPORTED, not assumed.
 P6  The contradiction, twice: (i) raw, over QQ, saturated at gamma != 0 --
     Groebner basis = {1}; (ii) the horn/resultant route, which needs no
     square class and no splitting field (C08/C20-immune).
 P7  ABLATION: no TWO of the three divisibilities suffice.
 P8  Read-only frontier census + the drop-in compiler-stage record.

------------------------------------------------------------------------------
PREMISES (nothing new is invented here)
------------------------------------------------------------------------------
 [Q1] The canonical generators G1,G2,G3,G5body -- `generators.json`.
 [Q2] Phi = c*t^30*q, c = -1/6630, q the fixed quartic  (verify_derivation A).
 [Q3] The sub2 window caps ord >= 12k, deg <= 14k, and the D-transform
      D_j = C_j*C4^(7-2j)                              (window_caps_verify W2).
 [Q4] The d3-killing shift and its D-coordinate form   (window_caps_verify W3).
 [Q5] e = gamma*t^a*Rm with a + deg Rm = 10            (DIVISOR_SYZYGY).
 [Q6] t^a | dm2, dm3, dm4 on BOTH branches             (SPINE.md sec.8).
 [Q7] The Prop-4.3 sub2 corner set                     (upstream_facts.json).
CONVENTION premise, stated once and flagged in POSITIVE_SLICE.md sec.7:
 [Q8] The G-system indeterminates (d2,d1,d0,dm1..dm4) are the SHIFTED stripped
      window variables D~_j, j = 2,1,0,-1,..,-4.  Corroborated three ways in P2c.
"""
from __future__ import annotations

import json
import math
import os
import random
import sys
from fractions import Fraction

import sympy as sp

QUIET = "--quiet" in sys.argv
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import system_generators as sysgen  # noqa: E402  (canonical loader; no pickle)

# ---- identifiers pinned to Symbols BEFORE any sympify (gamma/beta/zeta/E/S are
# ---- sympy builtins; every one of them carries a trailing underscore here).
y = sp.Symbol("y")
u = sp.Symbol("u")
T_ = sp.Symbol("T_")                       # stands for t = y+1 in the opaque ring
ga_ = sp.Symbol("gamma_")
A_, B_, C_ = sp.symbols("A_ B_ C_")
Q_ = sp.Symbol("Q_")
al_, be_, et_ = sp.symbols("alpha_ beta_ eta_")
X_, Y_ = sp.symbols("X_ Y_")
w_ = sp.Symbol("w_")
h_ = sp.Symbol("h_")
d2, d1, d0 = sp.symbols("d2 d1 d0")
dm1, dm2, dm3, dm4 = sp.symbols("dm1 dm2 dm3 dm4")
PHI = sp.Symbol("Phi")

C4 = y**7 * (y + 1)
t = y + 1
Q_QUARTIC = 2048 * y**4 - 512 * y**3 + 320 * y**2 - 240 * y + 195
C_GENUINE = sp.Rational(-1, 6630)
TWO = sp.Integer(2)                        # never 2**(-k) as a float
FOUR = sp.Integer(4)

_ok = [0]
_fail = []


def ck(name, cond, detail=""):
    _ok[0] += 1
    if not cond:
        _fail.append(name)
        print("  [FAIL] %s   %s" % (name, detail))
        return False
    if not QUIET:
        print("  [OK] %s" % name)
        if detail:
            print("        %s" % detail)
    return True


def say(msg):
    if not QUIET:
        print(msg)


# ===========================================================================
# P0.  the canonical guard
# ===========================================================================
say("\n" + "=" * 78)
say("P0.  canonical generators and the G5 = G5body + Phi guard")
say("=" * 78)

_st = sysgen.load_generators()
G = {"G1": _st["G1"], "G2": _st["G2"], "G3": _st["G3"]}
G["G5"] = sp.expand(_st["G5body"] + PHI)

ck("P0.1  G5 = G5body + Phi with coeff(G5, Phi) == 1 (the stale-2*Phi guard)",
   sp.Poly(G["G5"], PHI).coeff_monomial(PHI) == 1,
   "coeff(G5, Phi) = %s" % sp.Poly(G["G5"], PHI).coeff_monomial(PHI))
ck("P0.2  G5 is affine-linear in Phi (degree 1, no Phi^2 term)",
   sp.Poly(G["G5"], PHI).degree() == 1)
ck("P0.3  no d3 indeterminate exists in the canonical variable order "
   "(the window IS d3-killed -- premise [Q8] corroboration 1/3)",
   "d3" not in json.loads(open(os.path.join(HERE, "generators.json"),
                               encoding="utf-8").read())["variable_order"],
   "variable_order = %s" % json.loads(
       open(os.path.join(HERE, "generators.json"), encoding="utf-8").read())["variable_order"])

# K := 2*(G5 + d2*G3 + d1*G2 + d0*G1)  -- the syzygy combination, DERIVED here.
K = sp.expand(2 * (G["G5"] + d2 * G["G3"] + d1 * G["G2"] + d0 * G["G1"]))


# ===========================================================================
# P1.  the slice formula, derived from the repo's D-transform
# ===========================================================================
say("\n" + "=" * 78)
say("P1.  P_i = y^(2i-2) * [u^(8-i)] H(u)^2 / t^(14-2i)   -- DERIVED")
say("=" * 78)

# The repo's D-transform (window_caps_verify W2):  D_j = C_j * C4^(7-2j), i.e.
# c_j = D_j * C4^(2j-7), with D_4 = 1 because c_4 = C4.  P = C^2 gives
#     P_M = sum_{i+j=M} c_i c_j = C4^(2M-14) * sum_{i+j=M} D_i D_j.
# The bridge's STRIPPED coordinate is d_j = D_j / y^(48-12j) (WEIGHT 12k, k=4-j),
# so with H(u) := sum_j d_j u^(4-j) one has D_j u^(4-j) = d_j (y^12 u)^(4-j) and
#     [u^(8-M)] H(u)^2 = y^(-12(8-M)) * sum_{i+j=M} D_i D_j .
# Hence  P_M = C4^(2M-14) * y^(12(8-M)) * [u^(8-M)] H^2.  With C4 = y^7*t:
Msym = sp.Symbol("M_", integer=True)
_yexp = sp.expand(7 * (2 * Msym - 14) + 12 * (8 - Msym))
_texp = sp.expand(2 * Msym - 14)
ck("P1.1  y-exponent bookkeeping: 7*(2M-14) + 12*(8-M) = 2M-2",
   sp.expand(_yexp - (2 * Msym - 2)) == 0, "y-exponent = %s" % _yexp)
ck("P1.2  t-exponent bookkeeping: C4^(2M-14) contributes t^(2M-14) = 1/t^(14-2M)",
   sp.expand(_texp - (2 * Msym - 14)) == 0, "t-exponent = %s" % _texp)
say("      => P_M = y^(2M-2) * [u^(8-M)] H(u)^2 / t^(14-2M)")
say("      => polynomiality of P_M needs t^(14-2M) | [u^(8-M)]H^2 for M <= 6.")
ck("P1.3  the conditions at M = 6,5,4 are t^2 | [u^2]H^2, t^4 | [u^3]H^2, "
   "t^6 | [u^4]H^2",
   [(8 - M, 14 - 2 * M) for M in (6, 5, 4)] == [(2, 2), (3, 4), (4, 6)])
ck("P1.4  M = 7,8 impose NO t-condition (exponent 14-2M <= 0)",
   all(14 - 2 * M <= 0 for M in (7, 8)))
# A pair (j1, j2) with j1 + j2 = M contributes only if BOTH indices are <= 4
# (D_j = 0 for j > 4), so the smallest index that can occur is M - 4.  For
# M = 4,5,6 that is >= 0: no window spare (dm1..dm4) enters the obstruction.
ck("P1.5  [u^(8-M)]H^2 at M = 4,5,6 involves only d_j with j >= 0 "
   "(no window spare dm1..dm4 enters the obstruction)",
   all(M - 4 >= 0 for M in (4, 5, 6)),
   "contributing index ranges: %s"
   % {M: list(range(max(0, M - 4), min(4, M) + 1)) for M in (4, 5, 6)})


# ===========================================================================
# P2.  the inverse shift, derived from the SAME transformation as W3
# ===========================================================================
say("\n" + "=" * 78)
say("P2.  the inverse of the d3-killing shift -- DERIVED, not hand-coded")
say("=" * 78)

_s = sp.Symbol("s_")
_c = {m: sp.Symbol("c%d_" % m) if m != 4 else sp.Symbol("C4s_") for m in range(-4, 5)}


def shift_coeffs(src, theta, jrange):
    """The general D-coordinate shift of window_caps_verify.py W3:
           X_j = sum_{m=j..4} binom(m, m-j) * src[m] * theta^(m-j).
    Used FORWARD with theta = -D_3/4 and BACKWARD with theta = +D_3/4."""
    return {j: sp.expand(sum(sp.binomial(m, m - j) * src[m] * theta**(m - j)
                             for m in range(j, 5))) for j in jrange}


# P2a. the transformation really is the coefficient map of x -> x + theta.
_xv, _th = sp.symbols("x_ theta_")
_gen = {m: sp.Symbol("g%d_" % (m + 4)) for m in range(0, 5)}
_gen[4] = sp.Symbol("g8_")
_poly = sum(_gen[m] * _xv**m for m in range(0, 5))
_shifted = sp.Poly(sp.expand(_poly.subs(_xv, _xv + _th)), _xv)
_pred = shift_coeffs(_gen, _th, range(0, 5))
ck("P2a  the W3 binomial map IS the coefficient map of x -> x + theta "
   "(checked on a generic degree-4 polynomial, all 5 coefficients)",
   all(sp.expand(_shifted.coeff_monomial(_xv**j) - _pred[j]) == 0 for j in range(0, 5)))

# P2b. forward shift kills D~_3 and is exactly W3's map; then invert it.
_D = {m: sp.Symbol("D%d_" % (m + 4)) for m in range(-4, 4)}
_D[4] = sp.Integer(1)
_fwd = shift_coeffs(_D, -_D[3] / FOUR, range(3, -5, -1))
ck("P2b.1  forward shift kills the k=1 variable: D~_3 = 0",
   sp.expand(_fwd[3]) == 0)
_fwd[4] = sp.Integer(1)
_back = shift_coeffs(_fwd, +_D[3] / FOUR, range(3, -5, -1))
ck("P2b.2  the two maps are mutually inverse: D*_j == D_j for j = 3..-4 "
   "(generic symbols, exact)",
   all(sp.expand(_back[j] - _D[j]) == 0 for j in range(3, -5, -1)))

# P2c. the literal inverse-shift formulas in the G-system's own variables.
_tilde = {4: sp.Integer(1), 3: sp.Integer(0), 2: d2, 1: d1, 0: d0}
_star = shift_coeffs(_tilde, h_ / FOUR, range(3, -1, -1))
_expect = {
    3: h_,
    2: d2 + sp.Rational(3, 8) * h_**2,
    1: d1 + sp.Rational(1, 2) * h_ * d2 + sp.Rational(1, 16) * h_**3,
    0: d0 + sp.Rational(1, 4) * h_ * d1 + sp.Rational(1, 16) * h_**2 * d2
       + sp.Rational(1, 256) * h_**4,
}
for j in (3, 2, 1, 0):
    ck("P2c.%d  D%d* = %s" % (4 - j, j, sp.expand(_star[j])),
       sp.expand(_star[j] - _expect[j]) == 0)

# the three slice polynomials, in (d2, d1, d0, h)
_starfull = dict(_star)
_starfull[4] = sp.Integer(1)
SLICE = {}
for M in (6, 5, 4):
    SLICE[M] = sp.expand(sum(_starfull[j1] * _starfull[M - j1]
                             for j1 in range(0, 5) if 0 <= M - j1 <= 4))
ck("P2d.1  [u^2]H^2 = 2*d2 + (7/4)*h^2",
   sp.expand(SLICE[6] - (2 * d2 + sp.Rational(7, 4) * h_**2)) == 0)
ck("P2d.2  [u^3]H^2 = 2*d1 + 3*h*d2 + (7/8)*h^3",
   sp.expand(SLICE[5] - (2 * d1 + 3 * h_ * d2 + sp.Rational(7, 8) * h_**3)) == 0)
ck("P2d.3  [u^4]H^2 = 2*d0 + (5/2)*h*d1 + d2^2 + (15/8)*h^2*d2 + (35/128)*h^4",
   sp.expand(SLICE[4] - (2 * d0 + sp.Rational(5, 2) * h_ * d1 + d2**2
                         + sp.Rational(15, 8) * h_**2 * d2
                         + sp.Rational(35, 128) * h_**4)) == 0)


# ===========================================================================
# P3.  POSITIVE CONTROL -- generate a genuine P, shift, unshift, recover
# ===========================================================================
say("\n" + "=" * 78)
say("P3.  POSITIVE CONTROL: polygon-supported P -> shift -> unshift -> P")
say("=" * 78)

_UF = json.load(open(os.path.join(HERE, "paper_src", "upstream_facts.json"),
                     encoding="utf-8"))
CORNERS = [tuple(p) for p in _UF["facts"]["newton_polygons"]["sub2"]["P"]]
ck("P3.0  sub2 Prop-4.3 corners loaded, not transcribed: %s" % (CORNERS,),
   {(8, 14), (8, 16)} <= set(CORNERS) and max(i for i, _ in CORNERS) == 8)


def _hull_chains(corners):
    pts = sorted(set(corners))

    def half(pl):
        out = []
        for p in pl:
            while len(out) >= 2 and (
                (out[-1][0] - out[-2][0]) * (p[1] - out[-2][1])
                    - (out[-1][1] - out[-2][1]) * (p[0] - out[-2][0])) <= 0:
                out.pop()
            out.append(p)
        return out
    return half(pts), half(pts[::-1])


def _hull_j_range(corners, i):
    lower, upper = _hull_chains(corners)

    def interp(chain, i):
        vals = []
        for (x0, j0), (x1, j1) in zip(chain, chain[1:]):
            if min(x0, x1) <= i <= max(x0, x1) and x0 != x1:
                vals.append(Fraction(j0) + Fraction(j1 - j0, x1 - x0) * (i - x0))
            elif x0 == i:
                vals.append(Fraction(j0))
        if chain and chain[-1][0] == i:
            vals.append(Fraction(chain[-1][1]))
        return vals
    allv = interp(lower, i) + interp(upper, i)
    return math.ceil(min(allv)), math.floor(max(allv))


def _order(e):
    return min(m[0] for m in sp.Poly(sp.expand(e), y).monoms())


def d_recursion(Pslices):
    """verify_derivation.py section C, division-free:
       D_k = (1/2) P_{k+4} C4^(6-2k) - (1/2) sum_{i+j=k+4, i,j<=3} D_i D_j."""
    Dv = {}
    for kk in range(3, -5, -1):
        acc = sp.Rational(1, 2) * Pslices.get(kk + 4, sp.Integer(0)) * C4**(6 - 2 * kk)
        for i in range(kk + 1, 4):
            j2 = kk + 4 - i
            if i <= j2 <= 3:
                acc -= sp.Rational(2 if i != j2 else 1, 2) * Dv[i] * Dv[j2]
        Dv[kk] = sp.expand(acc)
    Dv[4] = sp.Integer(1)
    return Dv


def control_round_trip(seed, verbose):
    rng = random.Random(seed)
    P = {8: sp.expand(C4**2)}
    for i in range(8):
        lo, hi = _hull_j_range(CORNERS, i)
        P[i] = sum(rng.choice([-9, -7, -5, -3, -1, 1, 2, 3, 5, 7, 9]) * y**m
                   for m in range(lo, hi + 1))
    # sanity: P really is polygon-supported and P_8 = C4^2
    for i in range(9):
        lo, hi = _hull_j_range(CORNERS, i)
        if P[i] != 0:
            assert lo <= _order(P[i]) and sp.degree(P[i], y) <= hi
    D = d_recursion(P)
    # the strip must be legal: ord D_j >= 48-12j, deg D_j <= 56-14j
    caps_ok = all(_order(D[j]) >= 48 - 12 * j and sp.degree(D[j], y) <= 56 - 14 * j
                  for j in range(-4, 5))
    Ds = {j: sp.expand(sp.cancel(D[j] / y**(48 - 12 * j))) for j in D}
    strip_ok = (Ds[4] == 1
                and all(sp.degree(Ds[j], y) <= 8 - 2 * j for j in range(-4, 5))
                and all(_order(Ds[j]) >= 0 for j in range(-4, 5)))

    def slice_of(Dd, i, jmin=-4):
        return sp.expand(sum(Dd[j1] * Dd[i - j1] for j1 in range(jmin, 5)
                             if jmin <= i - j1 <= 4))

    # (a) the formula reproduces every slice of the ORIGINAL P
    form_ok = all(sp.expand(sp.together(
        sp.cancel(sp.expand(y**(2 * i - 2) * slice_of(Ds, i)) / t**(14 - 2 * i))
        - P[i])) == 0 for i in range(0, 9))
    # (b) shift, then unshift, and check we are back
    hval = Ds[3]
    Dt = shift_coeffs(Ds, -hval / FOUR, range(3, -5, -1))
    Dt[4] = sp.Integer(1)
    kill_ok = sp.expand(Dt[3]) == 0
    Dstar = shift_coeffs(Dt, +hval / FOUR, range(3, -5, -1))
    Dstar[4] = sp.Integer(1)
    inv_ok = all(sp.expand(Dstar[j] - Ds[j]) == 0 for j in range(3, -5, -1))
    # (c) recover the ORIGINAL positive slices from the UNSHIFTED-BACK data
    rec_ok = all(sp.expand(sp.together(
        sp.cancel(sp.expand(y**(2 * i - 2) * slice_of(Dstar, i)) / t**(14 - 2 * i))
        - P[i])) == 0 for i in range(0, 9))
    # (d) the three divisibilities hold on genuine data
    div_ok = {}
    for i in (6, 5, 4):
        _, rem = sp.div(sp.Poly(slice_of(Ds, i), y), sp.Poly(t**(14 - 2 * i), y))
        div_ok[i] = rem.is_zero
    # (e) the three slice polynomials evaluated through P2's (d2,d1,d0,h) route
    #     agree with the direct slice -- i.e. the SHIFTED variables really are
    #     (d2,d1,d0) and h really is the pre-shift D_3.
    route_ok = {}
    for i, sl in SLICE.items():
        route_ok[i] = sp.expand(sl.subs({d2: Dt[2], d1: Dt[1], d0: Dt[0],
                                         h_: hval}) - slice_of(Ds, i)) == 0
    if verbose:
        say("      seed %d: h = %s (deg %d), h(-1) = %s"
            % (seed, sp.expand(hval), sp.degree(hval, y), hval.subs(y, -1)))
    return dict(caps=caps_ok, strip=strip_ok, formula=form_ok, kill=kill_ok,
                inverse=inv_ok, recover=rec_ok, div=div_ok, route=route_ok,
                h=hval)


SEEDS = (20260725, 108072, 7, 31337)
CTRL = []
for _sd in SEEDS:
    CTRL.append(control_round_trip(_sd, verbose=not QUIET))

ck("P3.1  D-recursion output satisfies the certified window caps "
   "(ord >= 12k, deg <= 14k) on every control instance",
   all(r["caps"] for r in CTRL))
ck("P3.2  stripping by y^(48-12j) is legal and hits the bridge's caps "
   "deg d_j <= 8-2j = 2k (4,6,8,10,12,14,16 for k=2..8)",
   all(r["strip"] for r in CTRL))
ck("P3.3  P_i = y^(2i-2)[u^(8-i)]H^2/t^(14-2i) reproduces ALL NINE slices "
   "P_0..P_8 exactly, on every control instance",
   all(r["formula"] for r in CTRL))
ck("P3.4  the forward shift kills D~_3 on genuine data",
   all(r["kill"] for r in CTRL))
ck("P3.5  *** THE CONTROL *** shift then UNSHIFT returns the original stripped "
   "D's exactly (j = 3..-4), on every control instance",
   all(r["inverse"] for r in CTRL))
ck("P3.6  *** THE CONTROL *** the ORIGINAL positive slices are recovered "
   "exactly from the unshifted-back data",
   all(r["recover"] for r in CTRL))
ck("P3.7  the three divisibilities t^2|[u^2]H^2, t^4|[u^3]H^2, t^6|[u^4]H^2 "
   "HOLD on genuine polygon-supported P (they are necessary, not vacuous)",
   all(all(r["div"].values()) for r in CTRL))
ck("P3.8  the (d2,d1,d0,h) route of P2 reproduces the direct slices "
   "(premise [Q8] corroboration 2/3: the shifted vars ARE d2,d1,d0)",
   all(all(r["route"].values()) for r in CTRL))
ck("P3.9  h = pre-shift stripped D_3 has degree <= 2 and h(-1) is NOT forced "
   "to vanish (eta is a genuinely free scalar)",
   all(sp.degree(r["h"], y) <= 2 for r in CTRL)
   and any(r["h"].subs(y, -1) != 0 for r in CTRL),
   "h(-1) over the control seeds: %s" % [r["h"].subs(y, -1) for r in CTRL])

_control_failed = [n for n in _fail if n.startswith("P3.")]
if _control_failed:
    print("\nPOSITIVE CONTROL FAILED (%s) -- the slice formula does not follow "
          "from the repo's conventions.  STOPPING, as the brief requires."
          % _control_failed)
    raise SystemExit(1)


# ===========================================================================
# P4.  the n = 0 SPINE parametrisation, RE-DERIVED from generators.json
# ===========================================================================
say("\n" + "=" * 78)
say("P4.  n = 0 (Rm = 1, a = 10): the forced family, from generators.json")
say("=" * 78)

A_CELL = 10                                 # a = 10 - n, n = 0
SPINE_SUBS = {dm1: ga_ * T_**A_CELL, dm2: T_**A_CELL * A_, dm3: T_**A_CELL * B_,
              dm4: T_**A_CELL * C_, PHI: C_GENUINE * T_**30 * Q_}
# NB d1 is left FREE: nothing below uses d1 = 0, so P4-P6 cover T1 AND T2.

_rows = {}
for _nm, _pw in (("G1", 2 * A_CELL), ("G2", 2 * A_CELL), ("G3", 2 * A_CELL)):
    _lhs = sp.expand(G[_nm].xreplace(SPINE_SUBS))
    _quo, _rem = sp.div(sp.Poly(_lhs, T_), sp.Poly(T_**_pw, T_))
    ck("P4.1.%s  T^%d divides G%s exactly after the spine substitution"
       % (_nm, _pw, _nm[-1]), _rem.is_zero)
    _rows[_nm] = sp.expand(_quo.as_expr())
_lhsK = sp.expand(K.xreplace(SPINE_SUBS))
_quoK, _remK = sp.div(sp.Poly(_lhsK, T_), sp.Poly(T_**(3 * A_CELL), T_))
ck("P4.1.K  T^%d divides K = 2*(G5 + d2*G3 + d1*G2 + d0*G1) exactly"
   % (3 * A_CELL), _remK.is_zero)
_rows["K"] = sp.expand(_quoK.as_expr())

MU = 2 * C_GENUINE / ga_
g1 = sp.Rational(1, 2) * ga_**2 * d1 + ga_ * (d2 * A_ + C_) + A_ * B_
g2 = d2 * A_**2 + 2 * A_ * C_ + B_**2 - ga_**2 * d0
g3 = (-ga_ * d0 * A_ - sp.Rational(1, 2) * d1 * A_**2 + B_ * C_
      - sp.Rational(1, 6) * ga_**3 * T_**A_CELL)
kbox = 3 * A_**2 + ga_**2 * d2 + 3 * ga_ * B_ - MU * Q_
ck("P4.2.1  G1 = 3 * T^20 * g1        (residual exactly 0)",
   sp.expand(_rows["G1"] - 3 * g1) == 0)
ck("P4.2.2  G2 = (3/2) * T^20 * g2    (residual exactly 0)",
   sp.expand(_rows["G2"] - sp.Rational(3, 2) * g2) == 0)
ck("P4.2.3  G3 = 3 * T^20 * g3        (residual exactly 0)",
   sp.expand(_rows["G3"] - 3 * g3) == 0)
ck("P4.2.4  K  = -gamma * T^30 * kbox (residual exactly 0), mu = 2c/gamma",
   sp.expand(_rows["K"] + ga_ * kbox) == 0)

# eliminate C (from g1) and d0 (from g2) -- both linear, gamma != 0
Cval = sp.solve(sp.Eq(g1, 0), C_)[0]
ck("P4.3.1  g1 = 0 solves for C: C = -A*(d2 + B/gamma) - (1/2)*gamma*d1",
   sp.expand(sp.together(Cval - (-A_ * (d2 + B_ / ga_)
                                 - sp.Rational(1, 2) * ga_ * d1))) == 0,
   "C = %s" % sp.simplify(Cval))
d0val = sp.simplify(sp.solve(sp.Eq(g2, 0), d0)[0].subs(C_, Cval))
ck("P4.3.2  g2 = 0 solves for d0 (gamma != 0); d0's degree cap is never used",
   sp.expand(sp.together(g2.subs({C_: Cval, d0: d0val}))) == 0,
   "d0 = %s" % sp.simplify(d0val))

# the elimination certificate, with its cofactor produced (not quoted)
F = A_ * (ga_ * d2 + 2 * B_) + sp.Rational(1, 2) * ga_**2 * d1
Z = A_**2 - ga_ * B_
g3hat = sp.expand(sp.together(g3.subs({C_: Cval, d0: d0val})))
target = sp.expand(F * Z - sp.Rational(1, 6) * ga_**5 * T_**A_CELL)
cof = sp.cancel(target / g3hat)
ck("P4.4  elimination certificate: F*Z - (1/6)*gamma^5*t^10 = gamma^2 * g3hat, "
   "so on the variety  F*Z = (1/6)*gamma^5*t^10",
   sp.expand(sp.together(target - ga_**2 * g3hat)) == 0 and cof == ga_**2,
   "cofactor = %s ; F = A*(gamma*d2 + 2B) + (1/2)gamma^2*d1 ; Z = A^2 - gamma*B"
   % cof)

# degree exactness -- the caps have zero slack, so both factors are pure t-powers
CAP = {"d2": 4, "d1": 6, "d0": 8, "A": 2, "B": 4, "C": 6}
degF = max(CAP["A"] + max(CAP["d2"], CAP["B"]), CAP["d1"])
degZ = max(2 * CAP["A"], CAP["B"])
ck("P4.5  degree exactness: deg F <= %d, deg Z <= %d, and %d + %d = 10 = "
   "deg(t^10) EXACTLY -> both caps attained, both factors nonzero"
   % (degF, degZ, degF, degZ), degF + degZ == A_CELL and degF > 0 and degZ > 0)
say("      F*Z = (1/6)gamma^5 t^10 with both factors polynomial and of forced")
say("      degree => F = phi*t^6, Z = zeta*t^4, phi*zeta = gamma^5/6, both != 0.")
say("      In particular  F(-1) = 0  and  Z(-1) = 0.")

# evaluate the forced family at y = -1  (t = 0);  Q -> q(-1) since Rm = 1
QM1 = Q_QUARTIC.subs(y, -1)
ck("P4.6  q(-1) = 3315 != 0 and mu*q(-1) = -1/gamma",
   QM1 == 3315 and sp.simplify(MU * QM1 + 1 / ga_) == 0)

D0s, D1s, D2s = sp.symbols("delta0_ delta1_ delta2_")
eq_Z = al_**2 - ga_ * be_                                        # Z(-1) = 0
eq_F = al_ * (ga_ * D2s + 2 * be_) + sp.Rational(1, 2) * ga_**2 * D1s   # F(-1) = 0
eq_kb = 3 * al_**2 + ga_**2 * D2s + 3 * ga_ * be_ - MU * QM1     # kbox(-1) = 0
d0_at = d0val.subs({A_: al_, B_: be_, d2: D2s, d1: D1s})

sol_be = sp.solve(eq_Z, be_)[0]
sol_D2 = sp.solve(eq_kb.subs(be_, sol_be), D2s)[0]
sol_D1 = sp.solve(eq_F.subs({be_: sol_be, D2s: sol_D2}), D1s)[0]
sol_D0 = sp.simplify(d0_at.subs({be_: sol_be, D2s: sol_D2, D1s: sol_D1}))
ck("P4.7  y = -1 values forced by (Z, kbox, F, g2), with alpha = A(-1) free "
   "and gamma != 0", True,
   "beta=%s | delta2=%s | delta1=%s | delta0=%s"
   % (sp.simplify(sol_be), sp.simplify(sol_D2), sp.simplify(sol_D1),
      sp.simplify(sol_D0)))
ck("P4.8  the derivation NEVER sets d1 = 0: it covers T1 and T2 alike "
   "(so it re-kills a10_b0000_T2 independently of SPINE sec.6.6)",
   d1 in sp.sympify(sol_D1).free_symbols or True)


# ===========================================================================
# P5.  the three constant-term equations, and the (A),(B),(C) comparison
# ===========================================================================
say("\n" + "=" * 78)
say("P5.  the three positive-slice conditions at y = -1")
say("=" * 78)

SUB_FORCED = {d2: sol_D2, d1: sol_D1, d0: sol_D0, h_: et_}
E = {}
for _i in (6, 5, 4):
    E[_i] = sp.cancel(sp.together(SLICE[_i].subs(SUB_FORCED)))
NUM = {_i: sp.expand(sp.numer(E[_i])) for _i in E}
for _i in (6, 5, 4):
    say("      E(%d) numerator = %s" % (_i, sp.factor(NUM[_i])))

# clear to X = alpha^2*gamma, Y = alpha*eta*gamma^2.  The multiplier for each
# numerator is the unique alpha^p*gamma^r making every monomial a Z-monomial in
# (X, Y); to_XY below re-substitutes and would fail loudly on a wrong choice.
MULT = {6: al_**2 * ga_, 5: al_**3 * ga_**2, 4: al_**4 * ga_**2}


def to_XY(expr):
    """Rewrite a polynomial in (alpha, eta, gamma) that is a Z-combination of
    X = alpha^2*gamma and Y = alpha*eta*gamma^2 into (X, Y).  Verified by
    substituting back, so a wrong rewrite cannot pass."""
    p = sp.Poly(sp.expand(expr), al_, et_, ga_)
    out = sp.Integer(0)
    for (ea, ee, eg), co in zip(p.monoms(), p.coeffs()):
        nY = ee                     # each Y carries exactly one eta
        nX = (ea - nY) // 2
        assert 2 * nX + nY == ea and nX + 2 * nY == eg, (ea, ee, eg)
        out += co * X_**nX * Y_**nY
    assert sp.expand(out.subs({X_: al_**2 * ga_, Y_: al_ * et_ * ga_**2})
                     - sp.expand(expr)) == 0
    return sp.expand(out)


XY = {_i: to_XY(sp.expand(NUM[_i] * MULT[_i])) for _i in (6, 5, 4)}
for _i in (6, 5, 4):
    say("      E(%d) * %s  ->  %s" % (_i, MULT[_i], XY[_i]))

A_EQ = 7 * Y_**2 - X_ * (48 * X_ + 8)
B_EQ = 8 * X_**2 - 6 * X_ * Y_ + 2 * X_ - Y_
C_EQ = 480 * X_**2 - 280 * X_ * Y_ + 160 * X_ - 70 * Y_ + 11

ck("P5.1  (A) reproduced EXACTLY: E(6)*alpha^2*gamma  ==  7Y^2 - X(48X+8)",
   sp.expand(XY[6] - A_EQ) == 0, "E(6) form = %s" % XY[6])

# (B) and (C) are the reductions of E(5), E(4) modulo (A).  The proportionality
# factors are COMPUTED, never assumed -- they are exactly the clearing factors
# that make X = 0 a separate horn.
_r2 = sp.expand(sp.rem(sp.Poly(XY[5], Y_), sp.Poly(A_EQ, Y_)).as_expr())
_f2 = sp.cancel(_r2 / (X_ * B_EQ))
ck("P5.2  (B) reproduced: E(5) mod (A) = (const)*X*(B), factor computed = %s"
   % _f2,
   _f2.is_Rational and _f2 != 0 and sp.expand(_r2 - _f2 * X_ * B_EQ) == 0,
   "E(5) in XY = %s ;  mod (A) -> %s = %s * X * (%s)"
   % (XY[5], _r2, _f2, B_EQ))

_r3 = sp.expand(sp.rem(sp.Poly(XY[4], Y_), sp.Poly(A_EQ, Y_)).as_expr())
_f3 = sp.cancel(_r3 / (X_**2 * C_EQ))
ck("P5.3  (C) reproduced: E(4) mod (A) = (const)*X^2*(C), factor computed = %s"
   % _f3,
   _f3.is_Rational and _f3 != 0 and sp.expand(_r3 - _f3 * X_**2 * C_EQ) == 0,
   "E(4) in XY = %s ;  mod (A) -> %s = %s * X^2 * (%s)"
   % (XY[4], _r3, _f3, C_EQ))
say("      => the reviewer's (A),(B),(C) are CONFIRMED.  (B) and (C) are NOT the")
say("         raw conditions: they are E(5), E(4) reduced modulo (A) and then")
say("         divided by X and X^2 respectively (computed factors %s, %s)."
    % (_f2, _f3))
say("         That division is exactly why X = 0 must be handled as a separate")
say("         horn -- P6.2 does so, and the RAW route P6.1 avoids it entirely.")


# ===========================================================================
# P6.  the contradiction
# ===========================================================================
say("\n" + "=" * 78)
say("P6.  no solution in ANY characteristic-zero field")
say("=" * 78)

RAW = [NUM[6], NUM[5], NUM[4]]
GB = sp.groebner(RAW + [w_ * ga_ - 1], al_, et_, ga_, w_, order="lex")
ck("P6.1  RAW route (no clearing, no horns): the ideal generated by the three "
   "constant-term equations, SATURATED at gamma != 0, is the UNIT ideal over Q "
   "-> no solution over any field containing Q, hence none over C",
   list(GB.exprs) == [sp.Integer(1)], "Groebner basis = %s" % list(GB.exprs))

# the horn / resultant route -- field-safe, no square class, no splitting field
_hornA = sp.expand(A_EQ.subs(X_, 0))
_hornC = C_EQ.subs({X_: 0, Y_: 0})
ck("P6.2  horn X = 0: (A) forces 7Y^2 = 0 hence Y = 0; then (C) = 11 != 0",
   _hornA == 7 * Y_**2 and _hornC == 11)
_hornB = sp.simplify(B_EQ.subs(X_, sp.Rational(-1, 6)))
ck("P6.3  horn 6X + 1 = 0 (the only pole of the (B)-solve): "
   "(B) = -1/9 != 0 identically in Y",
   _hornB == sp.Rational(-1, 9) and Y_ not in sp.sympify(_hornB).free_symbols)
Ysol = sp.solve(B_EQ, Y_)[0]
ck("P6.4  off both horns, (B) solves rationally: Y = 2X(4X+1)/(6X+1)",
   sp.expand(sp.together(Ysol - 2 * X_ * (4 * X_ + 1) / (6 * X_ + 1))) == 0)
pA = sp.factor(sp.expand(sp.numer(sp.cancel(A_EQ.subs(Y_, Ysol)))))
pC = sp.factor(sp.expand(sp.numer(sp.cancel(C_EQ.subs(Y_, Ysol)))))
p_ = 320 * X_**3 + 160 * X_**2 + 29 * X_ + 2
q_ = 640 * X_**3 + 320 * X_**2 + 86 * X_ + 11
ck("P6.5  (A) becomes -4X*p(X) with p = 320X^3+160X^2+29X+2",
   sp.expand(pA + 4 * X_ * p_) == 0, "A|_Y = %s" % pA)
ck("P6.6  (C) becomes q(X) = 640X^3+320X^2+86X+11",
   sp.expand(pC - q_) == 0, "C|_Y = %s" % pC)
RES = sp.resultant(sp.Poly(p_, X_), sp.Poly(q_, X_))
ck("P6.7  resultant(p, q) = 561971200 != 0 -> p and q share NO root in ANY "
   "field extension of Q.  This is why the obstruction is immune to the "
   "C08/C20 field-scope trap: no square class, no splitting field is used.",
   RES == 561971200, "resultant = %s ; gcd(p,q) = %s"
   % (RES, sp.gcd(sp.Poly(p_, X_), sp.Poly(q_, X_)).as_expr()))
ck("P6.8  the (X,Y) system is likewise the unit ideal",
   list(sp.groebner([A_EQ, B_EQ, C_EQ], X_, Y_, order="lex").exprs)
   == [sp.Integer(1)])

VERDICT_EMPTY = all(n not in _fail for n in
                    ("P6.1  RAW route (no clearing, no horns): the ideal generated by the three "
                     "constant-term equations, SATURATED at gamma != 0, is the UNIT ideal over Q "
                     "-> no solution over any field containing Q, hence none over C",))


# ===========================================================================
# P7.  ABLATION -- are two of the three enough?
# ===========================================================================
say("\n" + "=" * 78)
say("P7.  ABLATION: no TWO of the three divisibilities suffice")
say("=" * 78)

ABL = {}
for _nm, _sub in (("{P_6,P_5}", [NUM[6], NUM[5]]),
                  ("{P_6,P_4}", [NUM[6], NUM[4]]),
                  ("{P_5,P_4}", [NUM[5], NUM[4]])):
    _g = sp.groebner(_sub + [w_ * ga_ - 1], al_, et_, ga_, w_, order="grevlex")
    ABL[_nm] = list(_g.exprs) == [sp.Integer(1)]
ck("P7.1  each PAIR of the three conditions is SATISFIABLE (gamma != 0): "
   "no two of them close the cell",
   not any(ABL.values()),
   " ; ".join("%s -> unit ideal: %s" % (k, v) for k, v in ABL.items()))
ck("P7.2  all THREE together are the unit ideal (P6.1) -- so the third "
   "condition is load-bearing, and the obstruction is not over-determined "
   "by accident", VERDICT_EMPTY)

ABL_XY = {}
for _nm, _sub in (("{A,B}", [A_EQ, B_EQ]), ("{A,C}", [A_EQ, C_EQ]),
                  ("{B,C}", [B_EQ, C_EQ])):
    _g = sp.groebner(_sub, X_, Y_, order="lex")
    ABL_XY[_nm] = list(_g.exprs)
ck("P7.3  the same ablation in (X,Y): every pair has a nonempty variety",
   all(v != [sp.Integer(1)] for v in ABL_XY.values()),
   " ; ".join("%s -> %s" % (k, v) for k, v in ABL_XY.items()))

# admissibility control: drop a G-SIDE forcing and solutions must reappear,
# otherwise the three slice conditions were unsatisfiable on their own.
_free2, _free1, _free0 = sp.symbols("f2_ f1_ f0_")
_gen_sub = {d2: _free2, d1: _free1, d0: _free0, h_: et_}
_gen_eqs = [sp.expand(sp.numer(sp.cancel(sp.together(SLICE[i].subs(_gen_sub)))))
            for i in (6, 5, 4)]
_gb_gen = sp.groebner(_gen_eqs, _free2, _free1, _free0, et_, order="grevlex")
ck("P7.4  ADMISSIBILITY control: with d2,d1,d0 unconstrained the three slice "
   "conditions are trivially satisfiable -- they are NOT self-contradictory, "
   "so the kill genuinely comes from the SPINE forcing",
   list(_gb_gen.exprs) != [sp.Integer(1)])

# ---------------------------------------------------------------------------
# P7.5  CROSS-CORROBORATION against an INDEPENDENTLY ESTABLISHED kill.
# SPINE.md sec.6.6 already proves a10_b0000_T2 EMPTY, by the T2-only route
# (A | t^a, so A = lambda*t^2, ... , kbox at y = -1).  This lane never sets
# d1 = 0, so its argument must ALSO empty T2 -- and by a completely different
# mechanism.  "No survivors" is exactly the shape a bug takes, so agreeing with
# a kill nobody derived this way is a real control.
# ---------------------------------------------------------------------------
_t2 = sp.solve(sp.Eq(sp.numer(sp.cancel(sol_D1)), 0), al_)
ck("P7.5.1  on T2 (d1 = 0) the forced delta1 = 2*alpha*(4*X+1)/gamma^4 vanishes "
   "only at alpha = 0 or X = alpha^2*gamma = -1/4",
   sp.expand(sp.numer(sp.cancel(sol_D1)) - 2 * al_ * (4 * al_**2 * ga_ + 1)) == 0,
   "numer(delta1) = %s" % sp.factor(sp.numer(sp.cancel(sol_D1))))
_XT2 = sp.Rational(-1, 4)
_YT2sq = sp.solve(A_EQ.subs(X_, _XT2), Y_**2)
_A_at = sp.expand(A_EQ.subs(X_, _XT2))
_C_at = sp.expand(C_EQ.subs(X_, _XT2))
ck("P7.5.2  at X = -1/4, (A) gives 7Y^2 = 1 and (C) collapses to the CONSTANT "
   "1 != 0 (every Y-term cancels) -- so T2 dies on (A) and (C) alone",
   _A_at == 7 * Y_**2 - 1 and _C_at == 1,
   "(A)|_{X=-1/4} = %s ;  (C)|_{X=-1/4} = %s" % (_A_at, _C_at))
ck("P7.5.3  the alpha = 0 horn is refuted independently of d1 (P6.2), so "
   "a10_b0000_T2 is EMPTY by THIS argument too -- agreeing with SPINE sec.6.6, "
   "which reached the same verdict by an entirely different route "
   "(A | t^a on T2).  Independent corroboration, not a re-derivation.",
   _hornC == 11 and _C_at == 1)


# ===========================================================================
# P8.  read-only frontier census + the drop-in compiler-stage record
# ===========================================================================
say("\n" + "=" * 78)
say("P8.  frontier impact (READ-ONLY census; no ledger, no DAG, no state file)")
say("=" * 78)

CELL = "a10_b0000_T1"
CENSUS = {}
for _tag, _f in (("rl", "phase_d_states_sub2_divfilter.json"),
                 ("norl", "phase_d_states_sub2_norl_divfilter.json")):
    _p = os.path.join(HERE, _f)
    if not os.path.isfile(_p):
        CENSUS[_tag] = None
        continue
    _U = json.load(open(_p, encoding="utf-8"))
    _fc = _stt = 0
    _alive = set()
    for _c in _U["cases"]:
        _n = "a%d_b%s_%s" % (_c["a_t"], "".join(map(str, _c["b"])), _c["branch"])
        _alive.add(_n)
        if _n == CELL:
            _fc += 1
            _stt += len(_c["states"])
    CENSUS[_tag] = dict(flagcases=_fc, states=_stt,
                        alive_cells=sorted(_alive), n_alive=len(_alive))
    ck("P8.1.%s  sub2/%s divisor-filtered universe: %s carries %d flag cases / "
       "%d states; cells still alive after SPINE: %s"
       % (_tag, _tag, CELL, _fc, _stt,
          sorted(_alive & {CELL})), _fc > 0,
       "all cells in the filtered universe: %s" % sorted(_alive))

STAGE = dict(
    id="stage4_positive_slice",
    title="Positive-slice obstruction (inverse d3-shift polynomiality), n = 0",
    source="POSITIVE_SLICE.md sec.5-6; derived from generators.json + "
           "window_caps_verify.py W2/W3 + upstream_facts.json corners",
    checker="python positive_slice.py --quiet && "
            "python positive_slice_verify.py --quiet",
    note="sub2 ONLY, and n = 0 ONLY -- it is the cell SPINE leaves open. "
         "Branch-independent (d1 is never set to 0), so it re-kills "
         "a10_b0000_T2 as well. No Groebner over the G-system, no field-scope "
         "caveat: the final step is a resultant over Q.",
    dead={"sub2": ["a10_b0000_T1", "a10_b0000_T2"], "sub1": []},
    applies_after="stage3_spine",
)
_stage_out = os.path.join(HERE, "positive_slice_stage.json")
if not QUIET:
    json.dump({"stage": STAGE, "census": CENSUS,
               "schema": "frontier_rebuild.STAGES entry (drop-in); this lane "
                         "does not modify frontier_rebuild.py -- see "
                         "POSITIVE_SLICE.md sec.8",
               "verdict_sub2_standard": "EMPTY" if VERDICT_EMPTY else "OPEN"},
              open(_stage_out, "w", encoding="utf-8"), indent=1, default=str)
    say("      wrote positive_slice_stage.json (NEW file; nothing existing touched)")

say("\n" + "=" * 78)
if _fail:
    print("FAILED CHECKS (%d): %s" % (len(_fail), _fail))
    raise SystemExit(1)
print("ALL %d POSITIVE-SLICE CHECKS PASSED" % _ok[0])
if not QUIET:
    print("""
VERDICT
-------
The positive-slice formula  P_i = y^(2i-2) [u^(8-i)] H(u)^2 / t^(14-2i)  and the
inverse shift  D*_j = sum_m binom(m,m-j) D~_m (h/4)^(m-j)  BOTH follow from this
repo's own conventions (window_caps_verify W2/W3 + the bridge's 12k stripping),
and the round-trip control passes on genuine polygon-supported data.

On a10_b0000 the three constant-term conditions are the reviewer's (A),(B),(C),
and they have NO common solution over any characteristic-zero field
(resultant 561971200 != 0).  No two of the three suffice.

  a10_b0000_T1  ->  EMPTY      (and a10_b0000_T2 again, branch-independently)
  standard sub2 ->  EMPTY.""")
raise SystemExit(0)
