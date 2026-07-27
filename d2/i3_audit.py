#!/usr/bin/env python3
"""i3_audit.py -- FOUNDATIONAL audit of premise [I3] / [Q8].

    python -u i3_audit.py            # full report
    python -u i3_audit.py --quiet    # exit 0 iff every check passes

TARGET (verbatim, as the two dependent lanes state it):

    the G-system indeterminates ARE the SHIFTED stripped coefficients
        dm1 = D~_{-1} = e ,  dm2 = D~_{-2} = R ,  dm3 = D~_{-3} = S ,
        dm4 = D~_{-4} = T .

    `POSITIVE_SLICE.md` sec.3.3 / sec.7 [Q8]  -- "a convention, not a theorem"
    `SLICE_OBSTRUCTION.md` sec.7 [Q8]          -- "convention"
    `ALT_LEVEL12.md` sec.6 [I3]                -- "the weakest link"

Everything below is re-derived from primitives.  NOTHING is imported from
`window_caps_verify.py`, `slice_obstruction_basis.py`, `positive_slice.py`,
`alt_level12.py`, `spine.py`, `divisor_syzygy.py` or `full_system_bridge.py`.
The only artifacts READ are `generators.json` (parsed by hand, not via
`system_generators`), and three source files read as TEXT by regex for their
declared constants.  Nothing is written.  Pure sympy; no Singular, no msolve,
no subprocess, no solver -- there are no exit codes to misread.

Read-only.  Creates no files.  Pinned against commit 3739c77.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import sympy as sp

QUIET = "--quiet" in sys.argv
ROOT = Path(__file__).resolve().parent

_n_ok = [0]
_fails: list[str] = []


def say(*a):
    if not QUIET:
        print(*a)


def head(t):
    say("\n" + "=" * 78 + "\n" + t + "\n" + "=" * 78)


def ck(name, cond):
    """A check.  `cond` MUST be a bool computed from an exact symbolic fact."""
    if cond is True:
        _n_ok[0] += 1
        say(f"  [OK]   {name}")
    else:
        _fails.append(name)
        print(f"  [FAIL] {name}")


y, t, u, x, th = sp.symbols("y t u x theta")
C4s = sp.Symbol("C4s")           # formal C4, for exponent bookkeeping
C4y = y**7 * (y + 1)             # the literal C4 = y^7*(y+1)

# =============================================================================
head("A.  The shift coefficient map, re-derived from scratch")
# =============================================================================
# A1.  generalized binomial from the FALLING FACTORIAL -- not sympy.binomial,
#      so no library convention is inherited.


def gbinom(m, r):
    """binom(m, r) for integer m (any sign) and integer r >= 0."""
    if r < 0:
        return sp.Integer(0)
    num = sp.Integer(1)
    for i in range(r):
        num *= sp.Integer(m - i)
    return sp.Rational(num, sp.factorial(r))


ck("A1  gbinom sanity: C(4,1)=4, C(4,4)=1, C(4,5)=0, C(-1,1)=-1, C(-1,2)=1, "
   "C(-2,1)=-2",
   gbinom(4, 1) == 4 and gbinom(4, 4) == 1 and gbinom(4, 5) == 0
   and gbinom(-1, 1) == -1 and gbinom(-1, 2) == 1 and gbinom(-2, 1) == -2)

# A2.  TRIANGULARITY ACROSS ZERO, and its exact scope.
_tri = all(gbinom(m, m - j) == 0 for m in range(0, 9) for j in range(-1, -9, -1))
ck("A2a TRIANGULARITY: gbinom(m, m-j) = 0 for every m in 0..8 and every "
   "j in -1..-8  (no NON-NEGATIVE source coefficient feeds any spare)", _tri)
# NON-VACUITY: the negative rows must genuinely survive, or A2a would be
# saying "everything vanishes" and the whole map would be trivial.
_surv = [(m, j) for m in range(-1, -6, -1) for j in range(m - 1, -9, -1)
         if gbinom(m, m - j) != 0]
ck("A2b NON-VACUITY of A2a: negative-index rows do NOT vanish -- "
   f"{len(_surv)} nonzero gbinom(m,m-j) with m<0 (e.g. gbinom(-1,1) = -1)",
   len(_surv) >= 10 and gbinom(-1, 1) != 0)
# NON-VACUITY 2: the map is not the identity -- non-negative rows DO mix.
ck("A2c NON-VACUITY: the map is NOT the identity -- gbinom(4,4-2) = 6 != 0, "
   "so D_4 feeds D~_2", gbinom(4, 4 - 2) == 6)

# A3.  Literal substitution x -> x + theta in a Laurent series, then compare.
NJ = 8                                    # keep c_4 .. c_{-8}
cc = {j: sp.Symbol(f"c{j}" if j >= 0 else f"cm{-j}") for j in range(4, -NJ - 1, -1)}
# (x+theta)^j = u^{-j} * (1+theta*u)^j        with u = 1/x
ORD = NJ + 6
_shifted = sp.Integer(0)
for j in range(4, -NJ - 1, -1):
    ser = sp.series((1 + th * u) ** j, u, 0, ORD).removeO()
    _shifted += cc[j] * u ** (-j) * sp.expand(ser)
_shifted = sp.expand(_shifted * u ** 4)   # coeff of u^(4-j) is now c~_j
_ok_A3 = True
for jv in range(4, -5, -1):
    lhs = sp.expand(_shifted.coeff(u, 4 - jv))
    rhs = sp.expand(sum(gbinom(m, m - jv) * cc[m] * th ** (m - jv)
                        for m in range(jv, 5)))
    if sp.expand(lhs - rhs) != 0:
        _ok_A3 = False
ck("A3  literal substitution x -> x+theta reproduces "
   "c~_j = sum_{m>=j} gbinom(m,m-j) c_m theta^(m-j)  for j = 4..-4", _ok_A3)


def shift_c(src, theta, jrange):
    """c~_j = sum_{m=j..4} gbinom(m,m-j) src[m] theta^(m-j)."""
    return {j: sp.expand(sum(gbinom(m, m - j) * src.get(m, 0) * theta ** (m - j)
                             for m in range(j, 5)))
            for j in jrange}


# A4.  D-coordinates: c_m = D_m * C4^(2m-7), and the x-shift that kills c~_3.
D = {m: sp.Symbol(f"D{m}" if m >= 0 else f"Dm{-m}") for m in range(3, -5, -1)}
D[4] = sp.Integer(1)
theta_x = -D[3] / (4 * C4s ** 2)          # = -c_3/(4*c_4);  see A4b
_c_of_D = {m: D[m] * C4s ** (2 * m - 7) for m in range(-4, 5)}
ck("A4a the x-shift that kills c~_3 is theta_x = -c_3/(4*c_4) = -D_3/(4*C4^2)",
   sp.simplify(-_c_of_D[3] / (4 * _c_of_D[4]) - theta_x) == 0)
# NON-VACUITY / notational finding: theta_x is NOT -D_3/4.
ck("A4b NOTE (documentation finding): theta_x != -D_3/4 -- the frequently "
   "written 'x -> x - D_3/4' is the D-COORDINATE parameter, not the x-shift",
   sp.simplify(theta_x - (-D[3] / 4)) != 0)

_ctil = shift_c(_c_of_D, theta_x, range(3, -5, -1))
_ok_A4 = True
_ok_A4poly = True
for jv in range(3, -5, -1):
    lhs = sp.expand(sp.cancel(sp.together(_ctil[jv] * C4s ** (7 - 2 * jv))))
    rhs = sp.expand(sum(gbinom(m, m - jv) * D[m] * (-D[3] / 4) ** (m - jv)
                        for m in range(jv, 5)))
    if sp.simplify(lhs - rhs) != 0:
        _ok_A4 = False
    if lhs.has(C4s):
        _ok_A4poly = False
ck("A4c D-COORDINATE FORM: c~_j*C4^(7-2j) = sum gbinom(m,m-j)*D_m*(-D_3/4)^(m-j) "
   "for j = 3..-4  -- the x-shift is RATIONAL, its D-coordinate form is not",
   _ok_A4)
ck("A4d every C4 exponent cancels: D~_j is a POLYNOMIAL in the D's "
   "(this is what makes the shifted window variables polynomial at all)",
   _ok_A4poly)

theta_D = -D[3] / 4
Dt = shift_c(D, theta_D, range(3, -5, -1))
ck("A5  the shift kills the k=1 variable: D~_3 = 0 identically", Dt[3] == 0)

# A6.  group law (the inverse shift is the same map at -theta).
_gen = {m: sp.Symbol(f"g{m}" if m >= 0 else f"gm{-m}") for m in range(4, -5, -1)}
_gen[4] = sp.Symbol("g4")
_fwd = shift_c(_gen, th, range(4, -5, -1))
_back = shift_c(_fwd, -th, range(4, -5, -1))
ck("A6  group law M(theta)*M(-theta) = I on rows j = 4..-4 (so the inverse "
   "shift needs no separate derivation)",
   all(sp.expand(_back[j] - _gen[j]) == 0 for j in range(4, -5, -1)))

# =============================================================================
head("B.  WHERE the identification is exact and where it MIXES\n"
     "    (with theta SPECIALIZED to -D_3/4 -- the audit's precision point)")
# =============================================================================
ck("B1a D~_{-1} = D_{-1} EXACTLY, with theta held INDEPENDENT",
   sp.expand(shift_c(D, th, [-1])[-1] - D[-1]) == 0)
ck("B1b D~_{-1} = D_{-1} EXACTLY, AFTER substituting theta = -D_3/4 "
   "(the specialization does NOT reintroduce D_3 at index -1)",
   sp.expand(Dt[-1] - D[-1]) == 0)
ck("B1c D~_{-1} is FREE of D_3 after specialization: d/dD_3 == 0",
   sp.expand(sp.diff(Dt[-1], D[3])) == 0)
# NON-VACUITY of B1c: index -2 is NOT free of D_3.
ck("B1d NON-VACUITY of B1c: d(D~_{-2})/dD_3 != 0 -- D_3 DOES reappear at "
   "index -2 after the specialization",
   sp.expand(sp.diff(Dt[-2], D[3])) != 0)

_mix = {}
for k in (2, 3, 4):
    _mix[k] = sp.expand(Dt[-k] - D[-k])
    say(f"       exact mixing term at index -{k}:  D~_(-{k}) - D_(-{k}) = "
        f"{sp.factor(_mix[k])}")
ck("B2a exact mixing at -2:  D~_{-2} = D_{-2} + (D_3/4)*D_{-1}",
   sp.expand(Dt[-2] - (D[-2] + D[3] * D[-1] / 4)) == 0)
ck("B2b exact mixing at -3:  D~_{-3} = D_{-3} + (1/2)*D_3*D_{-2} "
   "+ (1/16)*D_3^2*D_{-1}",
   sp.expand(Dt[-3] - (D[-3] + D[3] * D[-2] / 2
                       + sp.Rational(1, 16) * D[3] ** 2 * D[-1])) == 0)
ck("B2c exact mixing at -4:  D~_{-4} = D_{-4} + (3/4)*D_3*D_{-3} "
   "+ (3/16)*D_3^2*D_{-2} + (1/64)*D_3^3*D_{-1}",
   sp.expand(Dt[-4] - (D[-4] + sp.Rational(3, 4) * D[3] * D[-3]
                       + sp.Rational(3, 16) * D[3] ** 2 * D[-2]
                       + sp.Rational(1, 64) * D[3] ** 3 * D[-1])) == 0)
ck("B2d all three mixing terms are NONZERO (so 'exact at -1, mixes below' is "
   "a real dichotomy, not a vacuous one)",
   all(_mix[k] != 0 for k in (2, 3, 4)))

# B3.  the spares are fed ONLY by spares and by D_3; never by D_0,D_1,D_2.
_clean = all(sp.expand(sp.diff(Dt[-k], D[m])) == 0
             for k in (1, 2, 3, 4) for m in (0, 1, 2))
ck("B3a NO non-negative source coefficient (D_0,D_1,D_2) feeds ANY spare "
   "D~_{-1..-4}  -- triangularity, in the specialized map", _clean)
ck("B3b NON-VACUITY of B3a: the NON-NEGATIVE rows DO mix -- "
   "d(D~_2)/dD_3 != 0", sp.expand(sp.diff(Dt[2], D[3])) != 0)

# B4.  stripping commutes with the shift (ALT_LEVEL12 L4.6).
_d = {m: sp.Symbol(f"dd{m}" if m >= 0 else f"ddm{-m}") for m in range(3, -5, -1)}
_d[4] = sp.Integer(1)
_D_from_d = {m: (_d[m] * y ** (12 * (4 - m)) if m != 4 else sp.Integer(1))
             for m in range(4, -5, -1)}
_Dt_from_d = shift_c(_D_from_d, -_D_from_d[3] / 4, range(3, -5, -1))
_dt = shift_c(_d, -_d[3] / 4, range(3, -5, -1))
ck("B4a STRIPPING COMMUTES WITH THE SHIFT: (D~_j)/y^(12(4-j)) computed from "
   "unstripped D equals the same map applied to the stripped d, j = 3..-4",
   all(sp.expand(sp.cancel(_Dt_from_d[j] / y ** (12 * (4 - j))) - _dt[j]) == 0
       for j in range(3, -5, -1)))
ck("B4b stripping cannot move v_t: y is a UNIT at the place t = y+1 "
   "(y(-1) = -1 != 0), so v_t(y^N * f) = v_t(f)",
   sp.Integer(-1) != 0 and sp.simplify((y).subs(y, -1)) == -1)

# B5.  the INVERSE shift -- POSITIVE_SLICE.md sec.3.2's D2*, D1*, D0*.
_d2s, _d1s, _d0s, _h = sp.symbols("d2 d1 d0 h")
_tilde_src = {4: sp.Integer(1), 3: sp.Integer(0), 2: _d2s, 1: _d1s, 0: _d0s}
_star = shift_c(_tilde_src, _h / 4, range(3, -1, -1))   # inverse: theta = +h/4
ck("B5a inverse shift at theta = +D_3/4 reproduces POSITIVE_SLICE sec.3.2: "
   "D2* = d2 + (3/8)h^2, D1* = d1 + (1/2)h*d2 + (1/16)h^3, "
   "D0* = d0 + (1/4)h*d1 + (1/16)h^2*d2 + (1/256)h^4",
   sp.expand(_star[3] - _h) == 0
   and sp.expand(_star[2] - (_d2s + sp.Rational(3, 8) * _h ** 2)) == 0
   and sp.expand(_star[1] - (_d1s + _h * _d2s / 2
                             + sp.Rational(1, 16) * _h ** 3)) == 0
   and sp.expand(_star[0] - (_d0s + _h * _d1s / 4
                             + sp.Rational(1, 16) * _h ** 2 * _d2s
                             + sp.Rational(1, 256) * _h ** 4)) == 0)
ck("B5b LOAD-BEARING at indices 2,1,0: the h-corrections are NONZERO, so "
   "reading the G-system's d2,d1,d0 as UNSHIFTED would change the "
   "POSITIVE_SLICE equations",
   sp.expand(_star[2] - _d2s) != 0 and sp.expand(_star[1] - _d1s) != 0
   and sp.expand(_star[0] - _d0s) != 0)

# =============================================================================
head("C.  Is the shift ADMISSIBLE?  Cap preservation, from the caps alone")
# =============================================================================
jx, mm, jj = sp.symbols("jx mm jj")
for nm, cap in (("sub1 deg (60-15j)", 60 - 15 * jx),
                ("sub2 deg (56-14j)", 56 - 14 * jx),
                ("ord      (48-12j)", 48 - 12 * jx)):
    ck(f"C1  {nm}: cap(m) + (m-j)*cap(3) - cap(j) == 0 identically in (m,j) "
       "-- so deg[D_m*(-D_3/4)^(m-j)] <= cap(j) term by term",
       sp.expand(cap.subs(jx, mm) + (mm - jj) * cap.subs(jx, 3)
                 - cap.subs(jx, jj)) == 0)
_bad = 60 - 16 * jx
ck("C1z NON-VACUITY of C1: the identity FAILS for a wrong slope (60-16j), so "
   "C1 is a real property of the certified caps, not an algebraic triviality",
   sp.expand(_bad.subs(jx, mm) + (mm - jj) * _bad.subs(jx, 3)
             - _bad.subs(jx, jj)) != 0)

# C2.  numeric: worst-case D's sitting exactly ON the caps.
import random as _random


def _rand_poly(rng, ordf, degc):
    """random poly with ord exactly ordf and deg exactly degc."""
    if degc < ordf:
        return sp.Integer(0)
    co = [rng.choice([-3, -1, 1, 2, 3]) for _ in range(degc - ordf + 1)]
    co[0] = rng.choice([-3, -1, 1, 3])
    co[-1] = rng.choice([-3, -1, 1, 3])
    return sp.expand(sum(co[i] * y ** (ordf + i) for i in range(len(co))))


def _degord(p):
    if p == 0:
        return (-10 ** 6, 10 ** 6)
    P = sp.Poly(p, y)
    return P.degree(), min(m[0] for m in P.monoms())


for reg, slope in (("sub1", 15), ("sub2", 14)):
    rng = _random.Random(3739 if reg == "sub1" else 77)
    Dv = {4: sp.Integer(1)}
    for m in range(3, -5, -1):
        k = 4 - m
        Dv[m] = _rand_poly(rng, 12 * k, slope * k)      # ord = 12k, deg = 15k/14k
    Dtv = shift_c(Dv, -Dv[3] / 4, range(3, -5, -1))
    _capok = all(_degord(Dtv[m])[0] <= slope * (4 - m)
                 and _degord(Dtv[m])[1] >= 12 * (4 - m)
                 for m in range(3, -5, -1))
    ck(f"C2  {reg}: with EVERY unshifted D_m sitting exactly ON its certified "
       f"caps (ord = 12k, deg = {slope}k), the SHIFTED D~_m still obey "
       f"ord >= 12k and deg <= {slope}k  for k = 1..8 (worst case, zero slack)",
       _capok)
    # NON-VACUITY: inflate D_3 by one degree and the shifted caps must BREAK.
    Dbad = dict(Dv)
    Dbad[3] = sp.expand(Dv[3] + y ** (slope + 1))
    Dtb = shift_c(Dbad, -Dbad[3] / 4, range(3, -5, -1))
    _broke = any(_degord(Dtb[m])[0] > slope * (4 - m) for m in range(2, -5, -1))
    ck(f"C2z {reg} NON-VACUITY: giving D_3 ONE extra degree beyond its k=1 cap "
       "BREAKS the shifted degree caps -- C2 is a real consequence of the "
       "k=1 cap, not automatic", _broke)

# =============================================================================
head("D.  Does the PIPELINE actually use this convention?")
# =============================================================================
_gj = json.loads((ROOT / "generators.json").read_text(encoding="utf-8"))
_order = _gj["variable_order"]
ck("D1a generators.json variable_order = "
   "[d2,d1,d0,dm1,dm2,dm3,dm4,Phi] -- there is NO d3 row",
   _order == ["d2", "d1", "d0", "dm1", "dm2", "dm3", "dm4", "Phi"]
   and "d3" not in _order)

_vs = [sp.Symbol(n) for n in _order]


def _poly_of(terms):
    acc = sp.Integer(0)
    for cs, exps in terms:
        mono = sp.Rational(cs)
        for v, ex in zip(_vs, exps):
            mono *= v ** ex
        acc += mono
    return sp.expand(acc)


GEN = {k: _poly_of(v) for k, v in _gj["polynomials"].items()}
d2, d1, d0, dm1, dm2, dm3, dm4, Phi = _vs

# D2.  rebuild the G-system MYSELF from a series with NO u^1 term.
NSP = 13
dmv = {k: sp.Symbol(f"dm{k}") for k in range(1, NSP + 1)}


def _build(with_d3):
    d3sym = sp.Symbol("d3")
    S = (1 + (d3sym * u if with_d3 else 0) + d2 * u ** 2 + d1 * u ** 3
         + d0 * u ** 4 + sum(dmv[k] * u ** (4 + k) for k in range(1, NSP + 1)))
    S2 = sp.Poly(sp.expand(S * S), u)
    S3 = sp.Poly(sp.expand(S2.as_expr() * S), u)
    E2 = lambda k: S2.coeff_monomial(u ** (8 + k))
    E3 = lambda j: S3.coeff_monomial(u ** (12 + j))
    sub = {}
    for k, fresh in [(1, dmv[5]), (2, dmv[6]), (3, dmv[7]), (4, dmv[8]),
                     (5, dmv[9]), (6, dmv[10]), (7, dmv[11]), (9, dmv[13])]:
        sub[fresh] = sp.expand(sp.solve(E2(k).subs(sub), fresh)[0])
    return {"G1": sp.expand(E3(1).subs(sub)), "G2": sp.expand(E3(2).subs(sub)),
            "G3": sp.expand(E3(3).subs(sub)),
            "G5body": sp.expand(E3(5).subs(sub)), "E2": E2, "E3": E3}


_shifted_sys = _build(with_d3=False)
ck("D2a INDEPENDENT REBUILD: the four canonical rows built here from "
   "S = 1 + d2*u^2 + d1*u^3 + d0*u^4 + sum dm_k*u^(4+k)  (NO u^1 term, "
   "i.e. d3 = 0) reproduce generators.json's G1,G2,G3,G5body EXACTLY",
   all(sp.expand(_shifted_sys[g] - GEN[g]) == 0
       for g in ("G1", "G2", "G3", "G5body")))

_unshifted_sys = _build(with_d3=True)
_d3 = sp.Symbol("d3")
ck("D2b DISCRIMINATOR: rebuilding the SAME rows from a series that KEEPS the "
   "u^1 term gives polynomials that genuinely CONTAIN d3 -- so 'no d3 in "
   "generators.json' is informative, not an accident of notation",
   all(_unshifted_sys[g].has(_d3) for g in ("G1", "G2", "G3", "G5body")))
ck("D2c ... and setting d3 = 0 in those recovers the committed rows exactly. "
   "The G-system IS the d3 = 0 (i.e. SHIFTED) system, on the nose",
   all(sp.expand(_unshifted_sys[g].subs(_d3, 0) - GEN[g]) == 0
       for g in ("G1", "G2", "G3", "G5body")))

# D3.  index alignment: dm_k sits at u^(4+k), i.e. dm_k = d_{-k} = h_{4+k}.
_S = (1 + d2 * u ** 2 + d1 * u ** 3 + d0 * u ** 4
      + sum(dmv[k] * u ** (4 + k) for k in range(1, NSP + 1)))
ck("D3a INDEX ALIGNMENT: in the series the pipeline is built from, dm_k is "
   "the coefficient of u^(4+k) = u^(4-j) at j = -k.  So dm1 = d_{-1}, "
   "dm2 = d_{-2}, dm3 = d_{-3}, dm4 = d_{-4}",
   all(sp.expand(sp.Poly(_S, u).coeff_monomial(u ** (4 + k)) - dmv[k]) == 0
       for k in range(1, 5)))
_dgen = {j: sp.Symbol(f"dg{j}" if j >= 0 else f"dgm{-j}") for j in range(4, -5, -1)}
_Hgen = sum(_dgen[j] * u ** (4 - j) for j in range(4, -5, -1))
ck("D3b and in the cascade's level coordinates h_n := [u^n]H with "
   "H(u) = sum_j d_j u^(4-j): h_n = d_(4-n) for n = 0..8, so h_5 = d_{-1}, "
   "h_6 = d_{-2}, h_7 = d_{-3}, h_8 = d_{-4}.  With D3a: h_5 = dm1, "
   "h_6 = dm2, h_7 = dm3, h_8 = dm4",
   all(sp.expand(sp.Poly(_Hgen, u).coeff_monomial(u ** n) - _dgen[4 - n]) == 0
       for n in range(0, 9)))

# D4.  the consumers' constants, read as TEXT (never imported).
_fsb = (ROOT / "full_system_bridge.py").read_text(encoding="utf-8")
_wm = re.search(r"WEIGHT\s*=\s*\{(.*?)\}", _fsb, re.S)
_W = dict(re.findall(r'"(\w+)"\s*:\s*(\d+)', _wm.group(1))) if _wm else {}
ck("D4a full_system_bridge.WEIGHT (read by regex) gives dm_k the stripping "
   "weight 12*(4+k): dm1=60, dm2=72, dm3=84, dm4=96 -- consistent with "
   "dm_k = d_{-k} at window index k' = 4+k",
   all(int(_W.get(f"dm{k}", -1)) == 12 * (4 + k) for k in range(1, 5)))
_sd = re.search(r"STRIP_DEGCAP\s*=\s*\{(.*?)\n\}", _fsb, re.S)
_txt = _sd.group(1) if _sd else ""
_s1 = dict(re.findall(r'"(dm\d)"\s*:\s*(\d+)',
                      _txt.split('"sub1"')[1].split("}")[0])) if _sd else {}
_s2 = dict(re.findall(r'"(dm\d)"\s*:\s*(\d+)',
                      _txt.split('"sub2"')[1].split("}")[0])) if _sd else {}
ck("D4b full_system_bridge.STRIP_DEGCAP = lam*k with lam = 3 (sub1) / 2 "
   "(sub2) and k = 4+|j| -- exactly cap - 12k for the caps section C uses",
   all(int(_s1[f"dm{k}"]) == 3 * (4 + k) for k in range(2, 5))
   and all(int(_s2[f"dm{k}"]) == 2 * (4 + k) for k in range(2, 5)))

_ce = (ROOT / "cascade_engine.py").read_text(encoding="utf-8")
_c2 = re.search(r"SUB2\s*=\s*WindowConfig\(.*?aux_caps=\((\d+),\s*(\d+),\s*(\d+)\).*?e_cap=(\d+)", _ce)
_c1 = re.search(r"SUB1\s*=\s*WindowConfig\(.*?aux_caps=\((\d+),\s*(\d+),\s*(\d+)\).*?e_cap=(\d+)", _ce)
ck("D4c cascade_engine's OWN caps (read by regex) are lam*k at the right k: "
   "e_cap = 10 = 2*5 (sub2) / 15 = 3*5 (sub1); d2 cap = 2*2 / 3*2; "
   "d1 cap = 2*3 / 3*3.  So the CASCADE's (d2, d1, e) sit at window indices "
   "k = 2, 3, 5 -- i.e. e is d_{-1}, the SAME slot as dm1",
   _c2 is not None and _c1 is not None
   and (int(_c2.group(4)), int(_c2.group(3)), int(_c2.group(1))) == (2 * 5, 2 * 2, 2 * 3)
   and (int(_c1.group(4)), int(_c1.group(3)), int(_c1.group(1))) == (3 * 5, 3 * 2, 3 * 3))

# D5.  e, R, S, T are ALIASES, not a second identification.
_dsy = (ROOT / "divisor_syzygy.py").read_text(encoding="utf-8")
_alias = re.search(r"e,\s*R,\s*S,\s*T\s*=\s*sp\.symbols\(\s*[\"']dm1 dm2 dm3 dm4[\"']\s*\)",
                   _dsy)
ck("D5a divisor_syzygy.py DEFINES  e, R, S, T = symbols('dm1 dm2 dm3 dm4') -- "
   "so 'dm1 = e, dm2 = R, dm3 = S, dm4 = T' is a NAMING convention with NO "
   "mathematical content; the only substantive part of [I3] is shifted-vs-"
   "unshifted", _alias is not None)
_dfi = (ROOT / "divisor_filter.py").read_text(encoding="utf-8")
ck("D5b divisor_filter.py defines the cell coordinate by factoring THAT e: "
   "'e = gamma * t^a * prod_i (y - r_i)^{b_i} * (off-support factor)', so "
   "a_t = v_t(dm1) BY DEFINITION",
   "e = gamma * t^a * prod_i (y - r_i)^{b_i}" in _dfi)

# D6.  the K-syzygy, recomputed from generators.json alone.
_K = sp.expand(2 * (GEN["G5body"] + Phi + d2 * GEN["G3"] + d1 * GEN["G2"]
                    + d0 * GEN["G1"]))
_Krhs = sp.expand(2 * Phi - dm1 * (d2 * dm1 ** 2 + 3 * dm1 * dm3 + 3 * dm2 ** 2))
ck("D6a THE (K) SYZYGY, recomputed from generators.json: "
   "2*(G5 + d2*G3 + d1*G2 + d0*G1) == 2*Phi - e*(d2*e^2 + 3*e*S + 3*R^2) "
   "EXACTLY -- a pure identity in the G-system indeterminates",
   sp.expand(_K - _Krhs) == 0)
ck("D6b NON-VACUITY of D6a: the two sides are not both zero", _K != 0)
ck("D6c (K) mentions NO shift, NO theta, NO d3, and NO window coordinate: it "
   "is INTERNAL to the G-system, hence INDEPENDENT of [I3]'s "
   "shifted-vs-unshifted content",
   not _K.has(sp.Symbol("d3")) and not _K.has(th)
   and set(_K.free_symbols) <= set(_vs))

# D6d.  A FIFTH corroboration of [I3], not among POSITIVE_SLICE 3.3's four.
_Kbad = sp.expand(2 * (_unshifted_sys["G5body"] + Phi
                       + d2 * _unshifted_sys["G3"] + d1 * _unshifted_sys["G2"]
                       + d0 * _unshifted_sys["G1"]))
_res = sp.expand(_Kbad - _Krhs)
ck("D6d NEW CORROBORATION (not among POSITIVE_SLICE sec.3.3's four).  The (K) "
   "syzygy FAILS for the d3-present rows -- residual != 0 -- and the residual "
   "is EXACTLY divisible by d3, so it vanishes precisely on d3 = 0.  The "
   "whole divisor layer built on (K) (e | Phi, b_i in {0,1}, deg e = 10 in "
   "sub2, the place trichotomy) is therefore an identity of the SHIFTED "
   "system specifically.  Independent evidence that the pipeline is "
   "coherently d3-killed",
   _res != 0 and sp.expand(_res.subs(sp.Symbol("d3"), 0)) == 0
   and sp.rem(_res, sp.Symbol("d3"), sp.Symbol("d3")) == 0)

# =============================================================================
head("E.  Admissibility of the shift for the G-system DERIVATION\n"
     "    (the steps that are prose-only elsewhere in the repo)")
# =============================================================================
# E1.  Phi is SHIFT-INVARIANT:  F~_{-5} = F_{-5}.
Fm = {m: sp.Symbol(f"F{-m}") for m in range(-5, -14, -1)}
_Fsh = sp.Integer(0)
for m in range(-5, -14, -1):
    _Fsh += Fm[m] * u ** (-m) * sp.expand(sp.series((1 + th * u) ** m, u, 0, 12).removeO())
_Fsh = sp.expand(_Fsh)
ck("E1a Phi IS SHIFT-INVARIANT: for F with top x-degree -5 (premise [QQ1]), "
   "the shifted F~ has the same top coefficient F~_{-5} = F_{-5}.  Hence "
   "f1 = C4^3*F_{-5} and Phi = f1*C4^28 are unchanged by the shift",
   sp.expand(_Fsh.coeff(u, 5) - Fm[-5]) == 0)
ck("E1b NON-VACUITY of E1a: the NEXT coefficient does move -- "
   "F~_{-6} = F_{-6} - 5*theta*F_{-5} != F_{-6}",
   sp.expand(_Fsh.coeff(u, 6) - Fm[-6]) != 0
   and sp.expand(_Fsh.coeff(u, 6) - (Fm[-6] - 5 * th * Fm[-5])) == 0)
ck("E1c the shift cannot RAISE the top degree: F~ has no term of x-degree "
   "above -5", all(_Fsh.coeff(u, r) == 0 for r in range(0, 5)))

# E2.  no negative x-powers is preserved (the soundness of every G row).
_pc = sp.symbols("p0:9")
# expand P(x+theta) as a Laurent series in u = 1/x and demand every u^k, k >= 1,
# coefficient vanishes.  (x+theta)^i = u^(-i)*(1+theta*u)^i.
_Psh = sp.expand(sum(_pc[i] * u ** (-i)
                     * sp.series((1 + th * u) ** i, u, 0, 10).removeO()
                     for i in range(9)))
ck("E2a SOUNDNESS OF THE ROWS: P is a polynomial in x, so P(x+theta) has NO "
   "negative x-power for ANY theta -- including a theta RATIONAL in y.  Same "
   "for Q.  Every (D~^2)_{-k} = 0 and (D~^3)_{-j} = 0 row is therefore a "
   "genuine necessary condition on the SHIFTED coefficients",
   all(sp.expand(_Psh.coeff(u, k)) == 0 for k in range(1, 9)))
_Plau = sp.expand(_Psh + sp.Symbol("pm1") * u
                  * sp.series((1 + th * u) ** (-1), u, 0, 10).removeO())
ck("E2b NON-VACUITY of E2a: a genuinely LAURENT P (one x^-1 term) does NOT "
   "have that property after the shift, so E2a is a real consequence of "
   "P being polynomial in x",
   any(sp.expand(_Plau.coeff(u, k)) != 0 for k in range(1, 9)))

# E3.  the shift is USED, not cosmetic: it is what removes lambda from j = 5.
_cser = {4: C4s, 3: sp.Symbol("c3"), 2: sp.Symbol("c2"), 1: sp.Symbol("c1"),
         0: sp.Symbol("c0")}
_U = C4s + sum(_cser[3 - i] * u ** (i + 1) for i in range(4))
_inv = sp.series(1 / _U, u, 0, 4).removeO()
ck("E3a lambda-isolation: (C^-1)_{-4} = 1/C4 (a unit) and "
   "(C^-1)_{-5} = -c_3/C4^2",
   sp.simplify(_inv.coeff(u, 0) - 1 / C4s) == 0
   and sp.simplify(_inv.coeff(u, 1) + _cser[3] / C4s ** 2) == 0)
ck("E3b THE SHIFT IS LOAD-BEARING, NOT COSMETIC: (C^-1)_{-5} vanishes iff "
   "c_3 = 0.  Killing c_3 is exactly what makes the j=5 slice "
   "(C^3)_{-5} + F_{-5} = 0 free of lambda -- i.e. what makes the G5 row "
   "exist in its committed form",
   sp.simplify(_inv.coeff(u, 1).subs(_cser[3], 0)) == 0
   and sp.simplify(_inv.coeff(u, 1)) != 0)

# E4.  the shift is NOT a plane automorphism -- the divergence from upstream.
#      s = D_3/(4*C4^2) is polynomial iff C4^2 | D_3; under the k=1 caps
#      (ord >= 12, deg <= 15 sub1 / 14 sub2) that forces D_3 = 0.
_C4sq = sp.expand(C4y ** 2)
ck("E4a deg(C4^2) = 16 and ord(C4^2) = 14", _degord(_C4sq) == (16, 14))
_viol = []
for reg, cap in (("sub1", 15), ("sub2", 14)):
    # any nonzero multiple of C4^2 has degree >= 16 > cap
    _viol.append(16 > cap)
ck("E4b THE STATED HYPOTHESIS, MADE PRECISE.  s = D_3/(4*C4^2) is a "
   "POLYNOMIAL only if C4^2 | D_3; every nonzero multiple of C4^2 has "
   "degree >= 16, which exceeds the k=1 degree cap (15 sub1 / 14 sub2).  So "
   "D_3 = 0 is the only case: the d3-killing shift is NOT an automorphism of "
   "K[x,y].  It is an automorphism of K(y)[x] / K[y]((x^-1)) only",
   all(_viol))
ck("E4c NON-VACUITY of E4b: relax the k=1 cap to 16 and a nonzero admissible "
   "multiple of C4^2 EXISTS (namely C4^2 itself, ord 14 >= 12, deg 16)",
   _degord(_C4sq)[0] <= 16 and _degord(_C4sq)[1] >= 12)

# E4d.  The one-sided-ness, DEMONSTRATED rather than asserted.  Un-shifting a
#       point of the shifted system must reproduce POLYNOMIAL P-slices; that is
#       an extra demand the G-system does not make.  Exhibit a violation.
def _unshift_slice2(d2v, d1v, d0v, hv):
    """[u^2] of (the un-shifted H)^2 = 2*d2 + (7/4)*h^2  (POSITIVE_SLICE 3.2)."""
    st = shift_c({4: sp.Integer(1), 3: sp.Integer(0), 2: d2v, 1: d1v, 0: d0v},
                 hv / 4, range(3, -1, -1))
    Hs = 1 * u ** 0 + st[3] * u + st[2] * u ** 2 + st[1] * u ** 3 + st[0] * u ** 4
    return sp.expand(sp.expand(Hs * Hs).coeff(u, 2))


_viol_slice = _unshift_slice2(0, 0, 0, sp.Integer(1))
ck("E4d ONE-SIDED, DEMONSTRATED.  Take the shifted point d2 = d1 = d0 = 0 with "
   "h = 1 (h is invisible to the G-system: there is no d3 row).  The "
   "un-shifted slice [u^2]H*^2 = 2*d2 + (7/4)*h^2 = 7/4, so t does NOT divide "
   "it and P_6 = y^10*[u^2]H*^2/t^2 is NOT a polynomial.  The G-system objects "
   "to none of this => V(G-system) is STRICTLY larger than the germ image; "
   "emptiness still implies no germ, the converse fails",
   _viol_slice == sp.Rational(7, 4)
   and sp.simplify(_viol_slice.subs(y, -1)) != 0)
ck("E4e NON-VACUITY of E4d: with h = t instead, [u^2]H*^2 = (7/4)*t^2 IS "
   "divisible by t^2 -- so the condition E4d violates is a real, "
   "satisfiable one, not a contradiction manufactured for the demo",
   sp.expand(_unshift_slice2(0, 0, 0, t) - sp.Rational(7, 4) * t ** 2) == 0)

# E5.  The CASCADE's series is the UNSHIFTED one -- the pivot of [I3].
_dst = {j: sp.Symbol(f"e{j}" if j >= 0 else f"em{-j}") for j in range(3, -5, -1)}
_dst[4] = sp.Integer(1)
_C4 = y ** 7 * t
_cj = {j: sp.expand(_dst[j] * y ** (12 * (4 - j)) * _C4 ** (2 * j - 7))
       for j in range(4, -5, -1)}
ck("E5a the stripped substitution: c_j = d_j*y^(2j-1)*t^(2j-7), from "
   "D_j = c_j*C4^(7-2j), C4 = y^7*t and d_j = D_j/y^(12(4-j))",
   all(sp.simplify(_cj[j] - _dst[j] * y ** (2 * j - 1) * t ** (2 * j - 7)) == 0
       for j in range(4, -5, -1)))
_H = sum(_dst[j] * u ** (4 - j) for j in range(4, -5, -1))
_H2 = sp.expand(_H * _H)
_slice_ok = True
for M in range(8, -1, -1):
    PM = sp.expand(sum(_cj[i] * _cj[M - i] for i in range(-4, 5)
                       if -4 <= M - i <= 4))
    rhs = sp.expand(y ** (2 * M - 2) * sp.Poly(_H2, u).coeff_monomial(u ** (8 - M))
                    / t ** (14 - 2 * M))
    if sp.simplify(PM - rhs) != 0:
        _slice_ok = False
ck("E5b THE CASCADE'S H IS THE UNSHIFTED SERIES.  With H(u) = sum_j d_j u^(4-j) "
   "over j = 4..-4 INCLUDING j = 3, the slice identity "
   "P_M = y^(2M-2)*[u^(8-M)]H^2/t^(14-2M) holds exactly for M = 8..0.  So "
   "h_n := [u^n]H has h_1 = d_3 -- a coordinate that does NOT exist in the "
   "G-system (D1a)", _slice_ok)
ck("E5c THE PIVOT, STATED SHARPLY.  h_1 = d_3 is a live cascade coordinate; "
   "there is no d3 indeterminate in generators.json.  The two series are "
   "therefore DIFFERENT objects and [I3] is a substantive claim, not a "
   "relabelling", "d3" not in _order and sp.expand(_H.coeff(u, 1) - _dst[3]) == 0)
ck("E5d ... and the difference is visible in the data: the UNSHIFTED slice "
   "P_7 = 2*C4*c_3 is nonzero for c_3 != 0, while the SHIFTED P~_7 = "
   "2*C4*c~_3 = 0 identically",
   sp.simplify(sp.expand(sum(_cj[i] * _cj[7 - i] for i in range(3, 5)))
               - 2 * _C4 * _cj[3]) == 0
   and sp.expand(_cj[3]) != 0)

# =============================================================================
head("G.  Documentation hazards found while auditing (pinned to the source)")
# =============================================================================
_dict_block = _fsb[_fsb.index("VARIABLE DICTIONARY"):_fsb.index("STRIPPED COORDINATES")]
ck("G1  HAZARD: full_system_bridge.py's VARIABLE DICTIONARY writes the tilde "
   "on d~2, d~1, d~0 but NOT on d_-1..d_-4 -- i.e. the tilde is dropped at "
   "exactly the four indices where shifted and unshifted DIFFER.  The only "
   "unqualified tilde in the file is about the GENERATORS, (D~^3)_{-1,..}, "
   "not about the indeterminates",
   "d~2" in _dict_block and "d~1" in _dict_block and "d~0" in _dict_block
   and "d~-1" not in _dict_block and "d_-1" in _dict_block
   and "(D~^3)_{-1,-2,-3,-5}" in _fsb)
ck("G2  HAZARD: cascade_engine.py uses h_6, h_7 for its TERMINAL polynomials "
   "(h_7 = 8192*d1^2), a namespace collision with the slice lane's window "
   "coefficients h_6 = d_{-2}, h_7 = d_{-3}.  Same symbol, unrelated objects, "
   "in files that are read together",
   "h_7 = 8192*d1^2" in _ce and "h_slope" in _ce)
ck("G3  HAZARD: dm4 is named T in divisor_syzygy.py / divisor_consequences.py "
   "but M in full_system_bridge.SPARE_PREFIX and g4_row.py",
   'SPARE_PREFIX = {"dm2": "R", "dm3": "S", "dm4": "M"}' in _fsb
   and _alias is not None)

# =============================================================================
head("F.  BLAST RADIUS: what actually depends on [I3]")
# =============================================================================
# F1.  a_t >= 9.
ck("F1a a_t >= 9 IS [I3]-INVARIANT.  The bound is v_t(h_5) >= 9 with "
   "h_5 = d_{-1}; the census reads a_t = v_t(dm1).  Since D~_{-1} = D_{-1} "
   "IDENTICALLY (B1a/B1b), the shifted and unshifted readings of dm1 are the "
   "SAME polynomial, so the bound is unchanged under either reading",
   sp.expand(Dt[-1] - D[-1]) == 0 and sp.expand(shift_c(D, th, [-1])[-1] - D[-1]) == 0)

# F2.  the alternate regime, both readings.
A_T = {12: 9, 14: 8}       # horn-2 pinned v_t(R) per ALT_LEVEL12 sec.3
VH6 = 10                   # v_t(h_6) >= 10 from the LEVEL-10 cascade
VH1 = 1                    # v_t(h_1) >= 1 from level 2
_readA = {a: min(VH6, VH1 + a) for a in A_T}          # R = h6 + (h1/4)*h5
_readB = {a: VH6 for a in A_T}                        # R = h6
ck("F2a ALT six-branch kill under reading A (dm2 = D~_{-2}, the convention): "
   "R = h_6 + (h_1/4)*h_5 so v_t(R) >= min(10, 1+a_t) = 10 for a_t in "
   "{12,14}; horn 2 pins v_t(R) = 9 / 8.  10 > 9 and 10 > 8: CONTRADICTION",
   all(_readA[a] == 10 and _readA[a] > A_T[a] for a in A_T))
ck("F2b ALT six-branch kill under reading B (dm2 = D_{-2}, unshifted): "
   "v_t(R) = v_t(h_6) >= 10 directly; same contradiction",
   all(_readB[a] == 10 and _readB[a] > A_T[a] for a in A_T))
ck("F2c => the ALT kill is [I3]-INVARIANT: BOTH readings give v_t(R) >= 10 "
   "against a pinned 9 / 8", all(_readA[a] == _readB[a] for a in A_T))
# NON-VACUITY: the collision must FAIL where the frontier says it should.
_std = {9: min(VH6, VH1 + 9), 10: min(VH6, VH1 + 10)}
ck("F2d NON-VACUITY of F2a: the same arithmetic does NOT fire in the standard "
   "regime -- at a_t = 9 the bound degrades to min(10,10) = 10 while horn 1 "
   "is FEASIBLE there (3a <= 30), so no contradiction is manufactured",
   _std[9] == 10 and 3 * 9 <= 30 and 3 * 10 <= 30)
# and the symbolic version of F2a, not just the integer arithmetic:
_U1, _U5, _U6 = sp.symbols("U1 U5 U6")
for a in (12, 14):
    _h1 = t ** VH1 * _U1
    _h5 = t ** a * _U5
    _h6 = t ** VH6 * _U6
    _R = sp.expand(_h6 + _h1 * _h5 / 4)
    _ordR = min(sp.Poly(_R, t).monoms())[0] if _R != 0 else 10 ** 6
    ck(f"F2e symbolic, a_t = {a}: R = h_6 + (h_1/4)*h_5 has t-order exactly "
       f"{VH6} on generic units (the h_1*h_5 term sits at t^{1+a}, far above), "
       f"so v_t(R) = 10 != {A_T[a]}", _ordR == VH6)

# F3.  horn 1 exclusion -- from (K) ALONE, no slice calculus, no shift.
INF = 10 ** 9


def _min_twice(vals):
    """the minimum of the multiset is attained at least twice"""
    m = min(vals)
    return sum(1 for v in vals if v == m) >= 2


def _feasible(a, rho, PB=30, grid=range(0, 41)):
    """is there (v_s >= 0, delta2 >= 0) making (K)'s four orders attain their
    minimum at least twice?   orders: delta2+3m, 3m+v_s, m+2rho, P_b"""
    for vs in list(grid) + [INF]:
        for de in list(grid) + [INF]:
            o = [min(de + 3 * a, INF), min(3 * a + vs, INF),
                 (INF if rho >= INF else a + 2 * rho), PB]
            if _min_twice(o):
                return True
    return False


_h1_feas = {a: any(_feasible(a, rho) for rho in list(range(a, 3 * a + 1)) + [INF])
            for a in range(0, 21)}
ck("F3a HORN 1 (v_t(R) >= a_t) is INFEASIBLE for every a_t >= 11 -- "
   "3*a > 30 = v_t(Phi).  Exhaustive over rho in [a,3a] and rho = INF, "
   "v_s and delta2 over 0..40 and INF",
   all(not _h1_feas[a] for a in range(11, 21)))
ck("F3b NON-VACUITY of F3a: horn 1 is FEASIBLE for every a_t = 0..10, so the "
   "test is a threshold and not a constant refutation (a constant refutation "
   "would wrongly empty the standard regime too)",
   all(_h1_feas[a] for a in range(0, 11)))
_h2 = {a: sorted({rho for rho in range(0, a)
                  if _feasible(a, rho)}) for a in range(0, 21)}
_closed = {a: ([(30 - a) // 2] if (30 - a) % 2 == 0 and 0 <= (30 - a) // 2 < a
               else []) for a in range(0, 21)}
ck("F3c HORN 2 (v_t(R) = rho < a_t) forces 30 = a_t + 2*rho: the machine's "
   "feasible set equals the closed form {(30-a)/2} exactly, for a = 0..20",
   all(_h2[a] == _closed[a] for a in range(0, 21)))
ck("F3d => v_t(R) is PINNED: exactly 9 at a_t = 12, exactly 8 at a_t = 14",
   _h2[12] == [9] and _h2[14] == [8])
_Kmine = sp.expand(2 * (_shifted_sys["G5body"] + Phi + d2 * _shifted_sys["G3"]
                        + d1 * _shifted_sys["G2"] + d0 * _shifted_sys["G1"]))
ck("F3e THE SALVAGE CLAIM (ALT_LEVEL12.md sec.3) IS CORRECT.  Re-run: the (K) "
   "syzygy also holds for the rows I rebuilt MYSELF at D2a, and its free "
   "symbols are exactly the G-system indeterminates.  F3a-F3d consume only "
   "(K), v_t(Phi) = 30, e | S and d2 polynomial -- no shift, no theta, no "
   "slice calculus, no cascade.  Horn 1's exclusion SURVIVES a total failure "
   "of [I3]",
   sp.expand(_Kmine - _Krhs) == 0
   and set(_Kmine.free_symbols) <= set(_vs)
   and not _Kmine.has(sp.Symbol("d3")))

# F4.  what does NOT survive.
ck("F4  POSITIVE_SLICE.md (standard sub2 EMPTY) IS [I3]-DEPENDENT, by B5b: "
   "its three equations are built from D2* = d2 + (3/8)h^2 etc., which are "
   "the INVERSE shift.  Read the G-system's d2,d1,d0 as unshifted and those "
   "corrections vanish, changing (A),(B),(C).  Unlike a_t >= 9 and the ALT "
   "kill, this one has no [I3]-free reading",
   sp.expand(_star[2] - _d2s) != 0)

# F5.  v_t(Phi) = 30 exactly -- the other input to F3.
_q = 2048 * y ** 4 - 512 * y ** 3 + 320 * y ** 2 - 240 * y + 195
_Phi_str = sp.expand(-(y + 1) ** 30 * _q / 6630)
ck("F5  v_t(Phi) = 30 EXACTLY: Phi_stripped = -(1/6630)*t^30*q with "
   "q(-1) = 3315 != 0 and q separable",
   sp.simplify(_q.subs(y, -1)) == 3315 and sp.discriminant(_q, y) != 0
   and sp.simplify(sp.expand(_Phi_str / (y + 1) ** 30).subs(y, -1)) != 0)

# =============================================================================
head("SUMMARY")
# =============================================================================
say(f"checks passed : {_n_ok[0]}")
say(f"checks failed : {len(_fails)}")
for f in _fails:
    say("   FAILED:", f)
if _fails:
    print(f"i3_audit: {len(_fails)} FAILED of {_n_ok[0] + len(_fails)}")
    sys.exit(1)
print(f"i3_audit: ALL {_n_ok[0]} CHECKS PASSED")
sys.exit(0)
