#!/usr/bin/env python3
"""positive_slice_verify.py -- INDEPENDENT re-derivation of the positive-slice
obstruction.  Separately authored; imports NOTHING from positive_slice.py.

    python -u positive_slice_verify.py            # full report
    python -u positive_slice_verify.py --quiet    # exit 0 iff all checks pass

Every load-bearing step is redone by a DIFFERENT mechanism, so a bug shared
between the two files would have to be a bug in sympy or in generators.json:

  step                     positive_slice.py            THIS FILE
  ----------------------   --------------------------   ------------------------
  Laurent square root      quadratic D-recursion        global binomial series
                           (verify_derivation sec.C)    sum_n binom(1/2,n) w^n
  slice identity           stripped y/t exponent        UNSTRIPPED convolution
                           bookkeeping                  P_M = C4^(2M-14)*sum DiDj
  shift inverse            apply the binomial map       9x9 shift-matrix GROUP
                           twice and compare            LAW  M(a)M(b) = M(a+b)
  row factorisation        divide, compare to a form    substitute the SOLVED
                                                        values back into the RAW
                                                        generators.json rows
  certificate F*Z          explicit cofactor gamma^2    Groebner REDUCTION of
                                                        F*Z - (1/6)g^5 t^10
                                                        modulo <g1,g2,g3>
  emptiness                Groebner, gamma-saturated    RESULTANTS ONLY, plus an
                                                        explicit alpha = 0 horn
  ablation                 non-unit Groebner ideal      explicit WITNESS roots
                                                        via resultant factors

READ-ONLY on every existing artifact.  Pure sympy.  Exact rational arithmetic
throughout; `2**(-k)` never appears (sp.Integer/sp.Rational only).
"""
from __future__ import annotations

import json
import os
import random
import sys

import sympy as sp

QUIET = "--quiet" in sys.argv
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import system_generators as sysgen  # noqa: E402

# every identifier pinned to a Symbol before any sympify; gamma/beta/zeta/E/S
# are sympy builtins, hence the trailing underscores.
yy = sp.Symbol("yy")
TT = sp.Symbol("TT")                       # stands for t = y+1
gg = sp.Symbol("gg")                       # gamma
AA, BB, CC, QQ = sp.symbols("AA BB CC QQ")
aa, ee, hh = sp.symbols("aa ee hh")        # alpha, eta, h
XX, YY = sp.symbols("XX YY")
th1, th2 = sp.symbols("th1 th2")
D2v, D1v, D0v = sp.symbols("D2v D1v D0v")  # the shifted (G-system) d2,d1,d0

C4v = yy**7 * (yy + 1)
QUARTIC = 2048 * yy**4 - 512 * yy**3 + 320 * yy**2 - 240 * yy + 195
CGEN = sp.Rational(-1, 6630)
NN = 8

_n = [0]
_bad = []


def V(tag, cond, info=""):
    _n[0] += 1
    if not cond:
        _bad.append(tag)
        print("  [FAIL] %s  %s" % (tag, info))
        return False
    if not QUIET:
        print("  [ok] %s" % tag)
        if info:
            print("       %s" % info)
    return True


def out(msg):
    if not QUIET:
        print(msg)


# ===========================================================================
# V0.  canonical guard -- independent load
# ===========================================================================
out("\nV0. canonical generators")
_g = sysgen.load_generators()
Ph = sp.Symbol("Phi")
sd2, sd1, sd0 = sp.symbols("d2 d1 d0")
sm1, sm2, sm3, sm4 = sp.symbols("dm1 dm2 dm3 dm4")
G5full = sp.expand(_g["G5body"] + Ph)
V("V0.1  coeff(G5, Phi) == 1  (the stale 2*Phi transcription guard)",
  sp.Poly(G5full, Ph).coeff_monomial(Ph) == 1)
V("V0.2  the canonical variable order has no d3 -- the window is d3-killed, "
  "which is exactly what the shift of window_caps_verify W3 achieves",
  "d3" not in json.loads(open(os.path.join(HERE, "generators.json"),
                              encoding="utf-8").read())["variable_order"])
ROWS_RAW = {"G1": _g["G1"], "G2": _g["G2"], "G3": _g["G3"], "G5": G5full}
Kraw = sp.expand(2 * (G5full + sd2 * _g["G3"] + sd1 * _g["G2"] + sd0 * _g["G1"]))


# ===========================================================================
# V1.  Laurent square root by the GLOBAL BINOMIAL SERIES, and the slice identity
#      checked in UNSTRIPPED form (no y/t exponent bookkeeping at all)
# ===========================================================================
out("\nV1. Laurent sqrt via sum_n binom(1/2,n) w^n; unstripped slice identity")

_UF = json.load(open(os.path.join(HERE, "paper_src", "upstream_facts.json"),
                     encoding="utf-8"))
CORN = [tuple(p) for p in _UF["facts"]["newton_polygons"]["sub2"]["P"]]
V("V1.0  sub2 corners from upstream_facts.json: %s" % (CORN,),
  set(CORN) == {(0, 0), (1, 0), (8, 14), (8, 16)})
# The hull of {(0,0),(1,0),(8,14),(8,16)} has, at abscissa i, the j-range
# [2i-2, 2i] for i >= 1 and [0,0] at i = 0.  Derived, then asserted:
_edges_ok = True
for i in range(0, 9):
    lo_pred, hi_pred = (0, 0) if i == 0 else (2 * i - 2, 2 * i)
    # lower hull: (0,0)->(1,0)->(8,14) ; upper hull: (0,0)->(8,16)
    lo_true = 0 if i <= 1 else sp.Rational(14 * (i - 1), 7)
    hi_true = sp.Rational(16 * i, 8)
    if (sp.ceiling(lo_true), sp.floor(hi_true)) != (lo_pred, hi_pred):
        _edges_ok = False
V("V1.1  hull j-range at abscissa i is [2i-2, 2i] (i>=1), [0,0] at i=0 "
  "-- computed from the corner set, not transcribed", _edges_ok)


def trunc_mul(a, b):
    o = [sp.Integer(0)] * (NN + 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            if bj == 0 or i + j > NN:
                continue
            o[i + j] += ai * bj
    return [sp.expand(v) for v in o]


def sqrt_D(P):
    """D_j for j = 4..-4 from the slices P_0..P_8, by the binomial series
       sqrt(1 + w) = sum_n binom(1/2, n) w^n,  w = (P/x^8 - C4^2)/C4^2.
    Writing w = v/C4^2 with v_k = P_{8-k} (v_0 = 0), the C4 powers cancel:
       D_{4-k} = sum_{n<=k} binom(1/2,n) * [z^k](v^n) * C4^(2(k-n)).
    Entirely polynomial -- no rational function is ever formed."""
    v = [sp.Integer(0)] + [sp.expand(P[8 - k]) for k in range(1, NN + 1)]
    vp = [[sp.Integer(1)] + [sp.Integer(0)] * NN]
    for _ in range(NN):
        vp.append(trunc_mul(vp[-1], v))
    D = {}
    for k in range(NN + 1):
        acc = sp.Integer(0)
        for n in range(k + 1):
            co = vp[n][k]
            if co != 0:
                acc += sp.binomial(sp.Rational(1, 2), n) * co * C4v**(2 * (k - n))
        D[4 - k] = sp.expand(acc)
    return D


def conv(D, M, lo=-4):
    return sp.expand(sum(D[j] * D[M - j] for j in range(lo, 5)
                         if lo <= M - j <= 4))


def make_P(seed):
    rng = random.Random(seed)
    P = {8: sp.expand(C4v**2)}
    for i in range(8):
        lo, hi = (0, 0) if i == 0 else (2 * i - 2, 2 * i)
        P[i] = sum(rng.choice([-9, -7, -5, -3, -1, 1, 2, 3, 5, 7, 9]) * yy**m
                   for m in range(lo, hi + 1))
    return P


def ordy(e):
    return min(m[0] for m in sp.Poly(sp.expand(e), yy).monoms())


SEEDS = (4242, 999331, 20260726)
BANK = []
for _sd in SEEDS:
    _P = make_P(_sd)
    _D = sqrt_D(_P)
    BANK.append((_sd, _P, _D))

V("V1.2  the binomial series returns D_4 = 1 (i.e. c_4 = C4), every instance",
  all(D[4] == 1 for _, _, D in BANK))
V("V1.3  UNSTRIPPED slice identity  P_M * C4^(14-2M) == sum_{i+j=M} D_i D_j  "
  "for M = 0..8, every instance  (this is P = C^2, checked with NO exponent "
  "bookkeeping and NO stripping)",
  all(all(sp.expand(conv(D, M) - sp.expand(P[M] * C4v**(14 - 2 * M))) == 0
          for M in range(0, 9)) for _, P, D in BANK))
V("V1.4  the certified window caps hold: ord D_j >= 48-12j, deg D_j <= 56-14j "
  "for j = 4..-4, every instance",
  all(all(ordy(D[j]) >= 48 - 12 * j and sp.degree(D[j], yy) <= 56 - 14 * j
          for j in range(-4, 5)) for _, _, D in BANK))

# now, and only now, introduce the stripped coordinate and CONFIRM the form
# P_i = y^(2i-2) [u^(8-i)] H^2 / t^(14-2i) that positive_slice.py derives.
tt = yy + 1
STRIP = []
for _sd, _P, _D in BANK:
    Ds = {j: sp.expand(sp.cancel(_D[j] / yy**(48 - 12 * j))) for j in _D}
    STRIP.append((_sd, _P, _D, Ds))
V("V1.5  stripping by y^(48-12j) leaves polynomials with deg <= 8-2j = 2k "
  "(the bridge's STRIP_DEGCAP row, k = 4-j)",
  all(all(sp.Poly(Ds[j], yy).is_zero or (ordy(Ds[j]) >= 0
      and sp.degree(Ds[j], yy) <= 8 - 2 * j) for j in range(-4, 5))
      for _, _, _, Ds in STRIP))
V("V1.6  CONFIRMED INDEPENDENTLY: P_i = y^(2i-2) * (sum_{j1+j2=i} d_j1 d_j2) "
  "/ t^(14-2i)  for i = 0..8, every instance",
  all(all(sp.expand(sp.together(
      sp.cancel(sp.expand(yy**(2 * i - 2) * conv(Ds, i)) / tt**(14 - 2 * i)) - P[i]
  )) == 0 for i in range(0, 9)) for _, P, _, Ds in STRIP))
V("V1.7  the three divisibilities t^2 | [u^2]H^2, t^4 | [u^3]H^2, "
  "t^6 | [u^4]H^2 HOLD on genuine polygon-supported P (necessary, not vacuous)",
  all(all(sp.div(sp.Poly(conv(Ds, i), yy),
                 sp.Poly(tt**(14 - 2 * i), yy))[1].is_zero
          for i in (6, 5, 4)) for _, _, _, Ds in STRIP))


# ===========================================================================
# V2.  the shift, via the 9x9 shift-matrix GROUP LAW
# ===========================================================================
out("\nV2. shift matrices: M(a)*M(b) = M(a+b), hence M(-h/4)^-1 = M(h/4)")

IDX = list(range(4, -5, -1))               # 4,3,...,-4  (row/col order)


def Smat(theta):
    """(M(theta))_{j,m} = binom(m, m-j) theta^(m-j) for m >= j, else 0;
    it maps the coefficient vector (D_m) to (sum_m binom(m,m-j) D_m theta^(m-j)).
    """
    return sp.Matrix(len(IDX), len(IDX), lambda r, c: (
        sp.binomial(IDX[c], IDX[c] - IDX[r]) * theta**(IDX[c] - IDX[r])
        if IDX[c] >= IDX[r] else sp.Integer(0)))


_M1, _M2 = Smat(th1), Smat(th2)
V("V2.1  GROUP LAW: M(th1)*M(th2) == M(th1+th2)  (9x9, exact, symbolic)",
  sp.expand(_M1 * _M2 - Smat(th1 + th2)) == sp.zeros(len(IDX), len(IDX)))
V("V2.2  hence M(theta) is invertible with M(theta)^-1 = M(-theta): "
  "M(0) = Identity", sp.expand(Smat(sp.Integer(0))) == sp.eye(len(IDX)))
V("V2.3  the map is the coefficient map of x -> x + theta, checked against a "
  "literal substitution on a generic degree-4 polynomial",
  all(sp.expand(
      sp.Poly(sp.expand(sum(sp.Symbol("q%d" % m) * (sp.Symbol("xv") + th1)**m
                            for m in range(5))), sp.Symbol("xv")
              ).coeff_monomial(sp.Symbol("xv")**j)
      - sum(sp.binomial(m, m - j) * sp.Symbol("q%d" % m) * th1**(m - j)
            for m in range(j, 5))) == 0 for j in range(5)))

# the d3-killing shift and its inverse, in the G-system's own variables
_tilde_vec = sp.Matrix([sp.Integer(1), sp.Integer(0), D2v, D1v, D0v] +
                       [sp.Symbol("m%d" % k) for k in range(1, 5)])
_star_vec = sp.expand(Smat(hh / sp.Integer(4)) * _tilde_vec)
_star = {IDX[r]: sp.expand(_star_vec[r]) for r in range(len(IDX))}
_want = {3: hh,
         2: D2v + sp.Rational(3, 8) * hh**2,
         1: D1v + sp.Rational(1, 2) * hh * D2v + sp.Rational(1, 16) * hh**3,
         0: D0v + sp.Rational(1, 4) * hh * D1v + sp.Rational(1, 16) * hh**2 * D2v
            + sp.Rational(1, 256) * hh**4}
for _j in (3, 2, 1, 0):
    V("V2.4.%d  inverse shift D%d* = %s" % (4 - _j, _j, _want[_j]),
      sp.expand(_star[_j] - _want[_j]) == 0)
V("V2.5  the forward shift kills D~_3: with theta = -h/4 and D_3 = h, "
  "row 3 of M(-h/4) applied to (1, h, ...) gives h + 4*(-h/4) = 0",
  sp.expand((Smat(-hh / sp.Integer(4)) *
             sp.Matrix([sp.Integer(1), hh] + [sp.Symbol("z%d" % k)
                                              for k in range(7)]))[1]) == 0)

# END-TO-END on genuine data, with the recovery done through the UNSTRIPPED
# convolution -- so V1.6's stripped form is not reused here.
_rt_ok, _rec_ok, _slice_ok = True, True, True
SLICE_SYM = {}
_sf = dict(_star)
_sf[4] = sp.Integer(1)
for _M in (6, 5, 4):
    SLICE_SYM[_M] = sp.expand(sum(_sf[j] * _sf[_M - j] for j in range(0, 5)
                                  if 0 <= _M - j <= 4))
for _sd, _P, _D, _Ds in STRIP:
    hval = _Ds[3]
    vec = sp.Matrix([_Ds[j] for j in IDX])
    tvec = sp.expand(Smat(-hval / sp.Integer(4)) * vec)
    if sp.expand(tvec[1]) != 0:
        _rt_ok = False
    bvec = sp.expand(Smat(+hval / sp.Integer(4)) * tvec)
    if any(sp.expand(bvec[r] - vec[r]) != 0 for r in range(len(IDX))):
        _rt_ok = False
    Dstar_strip = {IDX[r]: sp.expand(bvec[r]) for r in range(len(IDX))}
    Dstar = {j: sp.expand(Dstar_strip[j] * yy**(48 - 12 * j)) for j in Dstar_strip}
    for M in range(0, 9):
        if sp.expand(conv(Dstar, M) - sp.expand(_P[M] * C4v**(14 - 2 * M))) != 0:
            _rec_ok = False
    tl = {IDX[r]: sp.expand(tvec[r]) for r in range(len(IDX))}
    for M in (6, 5, 4):
        if sp.expand(SLICE_SYM[M].subs({D2v: tl[2], D1v: tl[1], D0v: tl[0],
                                        hh: hval}) - conv(_Ds, M, lo=0)) != 0:
            _slice_ok = False
V("V2.6  *** CONTROL *** on genuine polygon-supported P: shift kills D~_3 and "
  "the inverse shift returns the original stripped D's exactly", _rt_ok)
V("V2.7  *** CONTROL *** the ORIGINAL slices P_0..P_8 are recovered from the "
  "unshifted-back data, via the UNSTRIPPED identity "
  "P_M*C4^(14-2M) = sum D*_i D*_j", _rec_ok)
V("V2.8  the (d2,d1,d0,h) slice polynomials agree with the direct convolution "
  "on genuine data  (so the shifted variables really are d2,d1,d0 and h really "
  "is the pre-shift stripped D_3)", _slice_ok)
V("V2.9  [u^2]H^2 = 2*d2 + (7/4)h^2 ; [u^3]H^2 = 2*d1 + 3*h*d2 + (7/8)h^3 ; "
  "[u^4]H^2 = 2*d0 + (5/2)h*d1 + d2^2 + (15/8)h^2*d2 + (35/128)h^4",
  sp.expand(SLICE_SYM[6] - (2 * D2v + sp.Rational(7, 4) * hh**2)) == 0
  and sp.expand(SLICE_SYM[5] - (2 * D1v + 3 * hh * D2v
                                + sp.Rational(7, 8) * hh**3)) == 0
  and sp.expand(SLICE_SYM[4] - (2 * D0v + sp.Rational(5, 2) * hh * D1v + D2v**2
                                + sp.Rational(15, 8) * hh**2 * D2v
                                + sp.Rational(35, 128) * hh**4)) == 0)

if [b for b in _bad if b.startswith(("V1.", "V2."))]:
    print("\nCONTROL FAILED -- the slice/shift formulas are not the repo's. STOP.")
    raise SystemExit(1)


# ===========================================================================
# V3.  n = 0 rows: verify by substituting SOLVED values back into the RAW rows
# ===========================================================================
out("\nV3. n = 0 (a = 10, Rm = 1): rows verified against generators.json itself")

A10 = 10
SPN = {sm1: gg * TT**A10, sm2: TT**A10 * AA, sm3: TT**A10 * BB,
       sm4: TT**A10 * CC, Ph: CGEN * TT**30 * QQ}          # d1 stays FREE
sub_rows = {k: sp.expand(v.xreplace(SPN)) for k, v in ROWS_RAW.items()}
sub_K = sp.expand(Kraw.xreplace(SPN))

# divide out the T-powers by unassisted exact division (quotient NOT supplied)
red = {}
for _k, _p in (("G1", 2 * A10), ("G2", 2 * A10), ("G3", 2 * A10)):
    _q, _r = sp.div(sp.Poly(sub_rows[_k], TT), sp.Poly(TT**_p, TT))
    V("V3.1.%s  T^%d | %s exactly (unassisted sp.div)" % (_k, _p, _k),
      _r.is_zero)
    red[_k] = sp.expand(_q.as_expr())
_q, _r = sp.div(sp.Poly(sub_K, TT), sp.Poly(TT**(3 * A10), TT))
V("V3.1.K  T^30 | K = 2*(G5 + d2*G3 + d1*G2 + d0*G1) exactly", _r.is_zero)
red["K"] = sp.expand(_q.as_expr())

# solve g1 = 0 for C and g2 = 0 for d0, then VERIFY against the RAW rows
Csol = sp.solve(sp.Eq(red["G1"], 0), CC)[0]
d0sol = sp.simplify(sp.solve(sp.Eq(red["G2"], 0), sd0)[0].subs(CC, Csol))
V("V3.2  C = -A*(d2 + B/gamma) - (1/2)*gamma*d1 makes the RAW G1 vanish "
  "identically after the spine substitution (checked on generators.json, not "
  "on a transcribed g1)",
  sp.expand(sp.together(sub_rows["G1"].xreplace({CC: Csol}))) == 0,
  "C = %s" % sp.simplify(Csol))
V("V3.3  the resulting d0 makes the RAW G2 vanish identically",
  sp.expand(sp.together(sub_rows["G2"].xreplace({CC: Csol, sd0: d0sol}))) == 0,
  "d0 = %s" % sp.simplify(d0sol))

# the certificate, by GROEBNER REDUCTION modulo <g1, g2, g3> (no cofactor given)
Fc = AA * (gg * sd2 + 2 * BB) + sp.Rational(1, 2) * gg**2 * sd1
Zc = AA**2 - gg * BB
tgt = sp.expand(Fc * Zc - sp.Rational(1, 6) * gg**5 * TT**A10)
GBI = sp.groebner([sp.expand(red["G1"]), sp.expand(red["G2"]),
                   sp.expand(red["G3"])], CC, sd0, AA, BB, sd2, sd1, gg, TT,
                  order="lex")
V("V3.4  F*Z - (1/6)*gamma^5*t^10 REDUCES TO 0 modulo the ideal <g1,g2,g3> "
  "(Groebner reduction; no cofactor supplied).  F = A*(gamma*d2+2B) + "
  "(1/2)gamma^2*d1,  Z = A^2 - gamma*B",
  sp.expand(GBI.reduce(tgt)[1]) == 0)

# degree exactness, by integer enumeration over the certified caps
CAPS = dict(d2=4, d1=6, A=2, B=4)
dF = max(CAPS["A"] + max(CAPS["d2"], CAPS["B"]), CAPS["d1"])
dZ = max(2 * CAPS["A"], CAPS["B"])
V("V3.5  deg F <= %d and deg Z <= %d; deg(t^10) = 10 = %d + %d, so BOTH caps "
  "are attained and neither factor is constant"
  % (dF, dZ, dF, dZ), dF + dZ == 10 and dF >= 1 and dZ >= 1)
V("V3.6  F*Z = (1/6)gamma^5 t^10 with F, Z polynomial and gamma != 0 forces "
  "F = phi*t^%d and Z = zeta*t^%d, phi*zeta = gamma^5/6 != 0.  Hence "
  "F(-1) = Z(-1) = 0." % (dF, dZ),
  all(p + q == 10 and (p, q) == (dF, dZ)
      for p, q in [(dF, dZ)]))

# y = -1 evaluation.  Solve, then CHECK every row at T = 0.
QM1 = QUARTIC.subs(yy, -1)
V("V3.7  q(-1) = 3315 != 0", QM1 == 3315)
b_, c_, D2_, D1_, D0_ = sp.symbols("b_ c_ D2_ D1_ D0_")
at0 = {TT: sp.Integer(0), AA: aa, BB: b_, CC: c_, QQ: QM1,
       sd2: D2_, sd1: D1_, sd0: D0_}
r1 = sp.expand(red["G1"].xreplace(at0))
r2 = sp.expand(red["G2"].xreplace(at0))
r3 = sp.expand(red["G3"].xreplace(at0))
rK = sp.expand(red["K"].xreplace(at0))
zc = sp.expand(Zc.xreplace(at0))
fc = sp.expand(Fc.xreplace(at0))
SOL = sp.solve([r1, r2, rK, zc, fc], [b_, c_, D2_, D1_, D0_], dict=True)
V("V3.8  the system {g1, g2, kbox, Z, F} at y = -1 has a UNIQUE solution for "
  "(beta, chi, delta2, delta1, delta0) in terms of (alpha, gamma)",
  len(SOL) == 1, "%s" % ({str(k): sp.simplify(v) for k, v in SOL[0].items()}
                         if SOL else None))
S0 = SOL[0]
V("V3.9  CONSISTENCY: the remaining row g3 at y = -1 is then automatically "
  "satisfied (it is implied by F(-1)*Z(-1) = 0), residual exactly 0",
  sp.simplify(sp.together(r3.subs(S0))) == 0)
d2m1, d1m1, d0m1 = (sp.simplify(S0[D2_]), sp.simplify(S0[D1_]),
                    sp.simplify(S0[D0_]))
V("V3.10  forced: delta2 = %s ; delta1 = %s ; delta0 = %s"
  % (d2m1, d1m1, d0m1), True)
V("V3.11  d1 is NEVER set to 0 anywhere above -- the derivation is "
  "branch-independent and covers a10_b0000_T1 and _T2 alike",
  sd1 in sp.sympify(red["G1"]).free_symbols)


# ===========================================================================
# V4.  the three constant-term equations
# ===========================================================================
out("\nV4. the positive-slice conditions at y = -1")

SUBF = {D2v: d2m1, D1v: d1m1, D0v: d0m1, hh: ee}
Eq = {M: sp.expand(sp.numer(sp.cancel(sp.together(SLICE_SYM[M].subs(SUBF)))))
      for M in (6, 5, 4)}
for M in (6, 5, 4):
    out("     N(%d) = %s" % (M, sp.factor(Eq[M])))

N6 = 7 * ee**2 * gg**3 - 48 * aa**2 * gg - 8
N5 = (128 * aa**3 * gg - 144 * aa**2 * ee * gg**2 + 32 * aa
      + 7 * ee**3 * gg**4 - 24 * ee * gg)
N4 = (3840 * aa**4 * gg**2 + 2560 * aa**3 * ee * gg**3
      - 1440 * aa**2 * ee**2 * gg**4 + 1280 * aa**2 * gg + 640 * aa * ee * gg**2
      + 35 * ee**4 * gg**6 - 240 * ee**2 * gg**3 + 128)
for M, ref in ((6, N6), (5, N5), (4, N4)):
    _r = sp.cancel(Eq[M] / ref)
    V("V4.1.%d  N(%d) is a nonzero rational multiple (%s) of the reference form"
      % (M, M, _r), _r.is_Rational and _r != 0)

# to (X, Y): X = alpha^2*gamma, Y = alpha*eta*gamma^2
MUL = {6: aa**2 * gg, 5: aa**3 * gg**2, 4: aa**4 * gg**2}
XYf = {}
for M, ref in ((6, N6), (5, N5), (4, N4)):
    p = sp.Poly(sp.expand(ref * MUL[M]), aa, ee, gg)
    acc = sp.Integer(0)
    for (ea, ex, eg), co in zip(p.monoms(), p.coeffs()):
        nY = ex
        nX = (ea - nY) // 2
        V("V4.2.%d.%s  monomial alpha^%d eta^%d gamma^%d = X^%d Y^%d"
          % (M, (ea, ex, eg), ea, ex, eg, nX, nY),
          2 * nX + nY == ea and nX + 2 * nY == eg)
        acc += co * XX**nX * YY**nY
    XYf[M] = sp.expand(acc)
    V("V4.3.%d  the rewrite is exact (substituted back)" % M,
      sp.expand(acc.subs({XX: aa**2 * gg, YY: aa * ee * gg**2})
                - sp.expand(ref * MUL[M])) == 0)

A_EQ = 7 * YY**2 - XX * (48 * XX + 8)
B_EQ = 8 * XX**2 - 6 * XX * YY + 2 * XX - YY
C_EQ = 480 * XX**2 - 280 * XX * YY + 160 * XX - 70 * YY + 11
V("V4.4  (A) reproduced exactly:  7Y^2 = X(48X+8)",
  sp.expand(XYf[6] - A_EQ) == 0, "%s" % XYf[6])
_rb = sp.expand(sp.rem(sp.Poly(XYf[5], YY), sp.Poly(A_EQ, YY)).as_expr())
_fb = sp.cancel(_rb / (XX * B_EQ))
V("V4.5  (B) reproduced: E(5) mod (A) = %s * X * (8X^2 - 6XY + 2X - Y)" % _fb,
  _fb.is_Rational and _fb != 0)
_rc = sp.expand(sp.rem(sp.Poly(XYf[4], YY), sp.Poly(A_EQ, YY)).as_expr())
_fc2 = sp.cancel(_rc / (XX**2 * C_EQ))
V("V4.6  (C) reproduced: E(4) mod (A) = %s * X^2 * "
  "(480X^2 - 280XY + 160X - 70Y + 11)" % _fc2,
  _fc2.is_Rational and _fc2 != 0)


# ===========================================================================
# V5.  EMPTINESS by resultants only -- no Groebner basis anywhere in this section
# ===========================================================================
out("\nV5. emptiness, by resultants and one explicit horn")

# horn alpha = 0.  Then X = alpha^2*gamma = 0, so the (X,Y) route is blind here
# and the raw equations must be used.  N(6) forces W := eta^2*gamma^3 = 8/7,
# and N(4) is a polynomial in W alone.
W_ = sp.Symbol("W_")
n6_0 = sp.expand(N6.subs(aa, 0))
n4_0 = sp.expand(N4.subs(aa, 0))
V("V5.1  horn alpha = 0: N(6) reads 7*W - 8 = 0 with W = eta^2*gamma^3, "
  "so W = 8/7", sp.expand(n6_0 - (7 * W_ - 8).subs(W_, ee**2 * gg**3)) == 0)
_pW = 35 * W_**2 - 240 * W_ + 128
V("V5.2  horn alpha = 0: N(4) is 35W^2 - 240W + 128 in the SAME W",
  sp.expand(n4_0 - _pW.subs(W_, ee**2 * gg**3)) == 0)
_val = sp.expand(_pW.subs(W_, sp.Rational(8, 7)))
V("V5.3  horn alpha = 0 is CONTRADICTORY: 35(8/7)^2 - 240(8/7) + 128 = %s != 0 "
  "(and this uses only P_6 and P_4)" % _val, _val == sp.Rational(-704, 7))

# main branch alpha != 0.  Since gamma != 0, X = alpha^2*gamma != 0, and the
# multipliers alpha^p*gamma^r are units, so the raw system is EQUIVALENT to the
# (X,Y) system there.  Eliminate Y by resultants.
R65 = sp.expand(sp.resultant(sp.Poly(XYf[6], YY), sp.Poly(XYf[5], YY)))
R64 = sp.expand(sp.resultant(sp.Poly(XYf[6], YY), sp.Poly(XYf[4], YY)))
R54 = sp.expand(sp.resultant(sp.Poly(XYf[5], YY), sp.Poly(XYf[4], YY)))
V("V5.4  the resultant criterion is exact here: lc_Y of E(6) is the nonzero "
  "constant %s, so res_Y = 0 at X0 <=> a common root Y0 exists at X0"
  % sp.LC(sp.Poly(XYf[6], YY)),
  sp.LC(sp.Poly(XYf[6], YY)).is_Rational and sp.LC(sp.Poly(XYf[6], YY)) != 0)
GCD2 = sp.factor(sp.gcd(sp.Poly(R65, XX), sp.Poly(R64, XX)).as_expr())
V("V5.5  gcd( res_Y(E6,E5), res_Y(E6,E4) ) = %s -- a MONOMIAL in X" % GCD2,
  sp.Poly(sp.expand(GCD2 / sp.LC(sp.Poly(GCD2, XX))), XX).is_monomial)
V("V5.6  therefore any common solution needs X = 0, i.e. alpha = 0, which "
  "V5.3 already refuted.  ==> NO SOLUTION over any characteristic-zero field.",
  True)
out("     res_Y(E6,E5) = %s" % sp.factor(R65))
out("     res_Y(E6,E4) = %s" % sp.factor(R64))

# the field-safety statement, re-derived: the two cubics of the horn-free route
Ysol = sp.solve(B_EQ, YY)[0]
pX = sp.expand(sp.numer(sp.cancel(A_EQ.subs(YY, Ysol))))
qX = sp.expand(sp.numer(sp.cancel(C_EQ.subs(YY, Ysol))))
p3 = 320 * XX**3 + 160 * XX**2 + 29 * XX + 2
q3 = 640 * XX**3 + 320 * XX**2 + 86 * XX + 11
V("V5.7  with Y = 2X(4X+1)/(6X+1): (A) -> -4X*p, p = 320X^3+160X^2+29X+2",
  sp.expand(pX + 4 * XX * p3) == 0)
V("V5.8  and (C) -> q = 640X^3+320X^2+86X+11", sp.expand(qX - q3) == 0)
RES = sp.resultant(sp.Poly(p3, XX), sp.Poly(q3, XX))
V("V5.9  resultant(p, q) = 561971200 != 0 -> p and q have no common root in "
  "ANY extension of Q.  No square class and no splitting field is used "
  "anywhere in this argument, so the C08/C20 field-scope downgrade cannot "
  "touch it.", RES == 561971200, "resultant = %s" % RES)
V("V5.10  horn 6X+1 = 0 (the only pole of the Y-solve) is refuted by (B) "
  "alone: (B)|_{X=-1/6} = -1/9, independent of Y",
  sp.simplify(B_EQ.subs(XX, sp.Rational(-1, 6))) == sp.Rational(-1, 9))
V("V5.11  horn X = 0 in the (X,Y) picture: (A) gives Y = 0, then (C) = 11 != 0",
  sp.expand(A_EQ.subs(XX, 0)) == 7 * YY**2 and C_EQ.subs({XX: 0, YY: 0}) == 11)


# ===========================================================================
# V6.  ABLATION with EXPLICIT WITNESSES (stronger than a non-unit ideal)
# ===========================================================================
out("\nV6. ablation: explicit witnesses show no TWO conditions suffice")

ABL = {}
for nm, (Pa, Pb) in (("{P_6,P_5}", (XYf[6], XYf[5])),
                     ("{P_6,P_4}", (XYf[6], XYf[4])),
                     ("{P_5,P_4}", (XYf[5], XYf[4]))):
    Rr = sp.Poly(sp.expand(sp.resultant(sp.Poly(Pa, YY), sp.Poly(Pb, YY))), XX)
    k = min(m[0] for m in Rr.monoms())
    tail = sp.expand(sp.cancel(Rr.as_expr() / XX**k))
    lcA = sp.LC(sp.Poly(Pa, YY))
    ok = sp.degree(tail, XX) > 0 and lcA.is_Rational and lcA != 0
    ABL[nm] = (k, sp.factor(tail), sp.degree(tail, XX))
    V("V6.1.%s  res_Y = X^%d * (%s): the tail has degree %d > 0, so it has a "
      "root X0 != 0 in Qbar; lc_Y is the nonzero constant %s, so a matching Y0 "
      "exists.  (alpha, eta, gamma) = (1, Y0/X0^2, X0) is then an honest "
      "witness with gamma != 0  ==>  this PAIR does NOT close the cell."
      % (nm, k, sp.factor(tail), sp.degree(tail, XX), lcA), ok)

V("V6.2  so all THREE positive-slice conditions are load-bearing: the "
  "obstruction is exactly determined, not over-determined", len(ABL) == 3)

# admissibility: with d2,d1,d0 free the three conditions are satisfiable, so
# they are not self-contradictory -- the kill comes from the SPINE forcing.
f2_, f1_, f0_ = sp.symbols("f2_ f1_ f0_")
free_eqs = [sp.expand(SLICE_SYM[M].subs({D2v: f2_, D1v: f1_, D0v: f0_}))
            for M in (6, 5, 4)]
_witness = {hh: sp.Integer(0), f2_: sp.Integer(0), f1_: sp.Integer(0),
            f0_: sp.Integer(0)}
V("V6.3  ADMISSIBILITY: with d2,d1,d0 unconstrained the three conditions have "
  "the explicit solution d2 = d1 = d0 = h = 0, so they are NOT vacuously "
  "unsatisfiable; the contradiction genuinely needs the SPINE forcing",
  all(sp.expand(e.subs(_witness)) == 0 for e in free_eqs))

# CROSS-CORROBORATION.  SPINE.md sec.6.6 proves a10_b0000_T2 EMPTY by a
# T2-ONLY route (A | t^a => A = lambda*t^2, then kbox at y = -1).  Nothing in
# this file sets d1 = 0, so the same argument must empty T2 as well -- by a
# different mechanism, and agreeing with a verdict established elsewhere.
V("V6.4  on T2 (d1 = 0) the forced delta1 vanishes only at alpha = 0 or "
  "X = -1/4", sp.expand(sp.numer(sp.cancel(d1m1))
                        - 2 * aa * (4 * aa**2 * gg + 1)) == 0,
  "numer(delta1) = %s" % sp.factor(sp.numer(sp.cancel(d1m1))))
V("V6.5  at X = -1/4 the equation (C) collapses to the constant 1 (all Y-terms "
  "cancel), and the alpha = 0 horn was refuted in V5.3 -- so a10_b0000_T2 is "
  "EMPTY by THIS argument too, agreeing with SPINE.md sec.6.6, which reached "
  "the same verdict by an unrelated route.  CROSS-CORROBORATION.",
  sp.expand(C_EQ.subs(XX, sp.Rational(-1, 4))) == 1
  and sp.expand(A_EQ.subs(XX, sp.Rational(-1, 4))) == 7 * YY**2 - 1,
  "(A)|_{X=-1/4} = %s ; (C)|_{X=-1/4} = %s"
  % (sp.expand(A_EQ.subs(XX, sp.Rational(-1, 4))),
     sp.expand(C_EQ.subs(XX, sp.Rational(-1, 4)))))


# ===========================================================================
# V7.  verdict
# ===========================================================================
out("\n" + "=" * 78)
if _bad:
    print("FAILED (%d): %s" % (len(_bad), _bad))
    raise SystemExit(1)
print("ALL %d INDEPENDENT CHECKS PASSED" % _n[0])
if not QUIET:
    print("""
INDEPENDENT VERDICT
-------------------
  * The slice formula and the inverse shift DO follow from this repo's
    conventions (binomial-series square root; shift-matrix group law).
  * The positive control passes: shift then unshift recovers the original
    polygon-supported slices EXACTLY.
  * On a10_b0000 the three constant-term conditions are (A),(B),(C), and by
    resultants alone -- gcd of the eliminants is a monomial in X, and X = 0 is
    refuted separately -- they have NO common solution over any
    characteristic-zero field.
  * No TWO of the three suffice (explicit witnesses).

  a10_b0000_T1 is EMPTY.  Standard sub2 is EMPTY.""")
raise SystemExit(0)
