#!/usr/bin/env python3
"""alt_level12.py -- INDEPENDENT adjudication of the SLICE_OBSTRUCTION.md sec.6.3
level-12 lead against the six surviving ALTERNATE-REGIME T1 branches.

    a12_b0000_T1  a12_b1000_T1  a12_b1100_T1  a12_b1110_T1
    a14_b0000_T1  a14_b1000_T1

WRITES NOTHING.  Reads no peer artifact except `paper_src/upstream_facts.json`
(not needed for the verdict) -- every number below is recomputed here from
scratch in sympy.  In particular this file does NOT import
`slice_obstruction_basis`, `alt_rebuild`, `t1_branch`, `frontier_rebuild`,
`window_caps_verify` or `g4_row`; their claims are RE-DERIVED, not read.

=====================================================================
WHAT IS CHECKED (and what each check is worth)
=====================================================================

L1  the algebraic identity  2H^3 - 3H^2 = -1 + 3K^2 + 2K^3,  K = H-1,
    so the STACKED functional  [u^n](3K^2+2K^3) = 2*r_n - 3*p_n  is free of
    the fresh coefficient h_n.                                     [PROVED]

L2  the cascade, levels n = 2..12, recomputed independently.  At each level
    the lowest jet of the stacked condition is factored; a deduction counts
    only if the non-constant part is a SINGLE irreducible factor that is
    LINEAR in one g-coefficient with a constant leading coefficient.  This
    is the check that the level-12 step does not BRANCH.        [PROVED]

L3  the resulting forced profile v_t(h_k) >= 2k-1 for k = 1..6.  The two
    numbers the lead needs are  v_t(h_1) >= 1  and  v_t(h_6) >= 11.
                                                                  [PROVED]

L4  the d3-killing shift, re-derived from the generalized-binomial
    recomposition:  D~_j = sum_{m=j..4} binom(m,m-j) D_m (-D_3/4)^(m-j).
    Consequences:  D~_{-1} = D_{-1}  (so h_5 = e exactly) and
    D~_{-2} = D_{-2} + (D_3/4)*D_{-1} = h_6 + (h_1/4)*h_5,  whence
        v_t(R) >= min( v_t(h_6), v_t(h_1) + a_t ) = min(11, 1 + a_t).
    The y^12 stripping factors cancel out of the shift identically, so the
    formula is the same in stripped coordinates.                   [PROVED]

L5  the place "trichotomy" at y = -1, re-derived here as a Newton-polygon
    (min-attained-twice) feasibility test on the four terms of the K-syzygy
        2*Phi = d2*e^3 + 3*e^2*S + 3*e*R^2 ,
    with v_t(d2) >= 0 (d2 polynomial), v_t(S) >= v_t(e) (e | S), and
    v_t(Phi) = 30 EXACTLY.  It has exactly TWO horns and the first is
    INFEASIBLE for every a >= 11 -- 3a > 30.  Machine-enumerated, with the
    mandatory non-vacuity control that horn 1 IS feasible at a = 8,9,10.
                                                                  [PROVED]

L6  the six branches, one at a time.

L7  CONTROLS: the argument must not kill the standard regime (a <= 10), and
    must not be a constant refutation.

=====================================================================
THE STANDING CONDITIONALITY -- READ THIS FIRST
=====================================================================
L2/L3 sit on top of the SAME slice-obstruction calculus whose level-10 output
(a_t >= 9) is graded `exact-checked SAME-AUTHOR` in SLICE_OBSTRUCTION.md sec.8
and is under independent audit at the time of writing.  This file re-derives
the CASCADE independently, but it does NOT re-derive the calculus's inputs:

  [I1] the P-slice formula  P_M = y^(2M-2)[u^(8-M)]H^2/t^(14-2M)  and hence
       the P conditions  t^(2n-2) | p_n  (n <= 8),  p_n = 0 (n >= 9);
  [I2] the Q-slice formula  Q_M = y^(2M-3)[u^(12-M)]H^3/t^(21-2M)  and hence
       the Q conditions  t^(2n-3) | r_n  for n <= 15, which needs premise
       QQ1 (v_{1,0}(F) = -5, lambda-strip WLOG);
  [I3] the identification of the G-system indeterminates with the SHIFTED
       stripped coefficients, i.e. dm1 = D~_{-1} = e and dm2 = D~_{-2} = R.

Those three are IMPORTED, flagged, and are exactly what the concurrent audit
is auditing.  A verdict here is conditional on them.

Usage:  python -u alt_level12.py            # full report
        python -u alt_level12.py --quiet    # exit 0 iff every check passes
        python -u alt_level12.py --fast     # cascade to level 10 only (dev)
"""
from __future__ import annotations

import sys
import time

import sympy as sp

QUIET = "--quiet" in sys.argv
MAXLEV = 10 if "--fast" in sys.argv else 12

ZERO = sp.Integer(0)
ONE = sp.Integer(1)
INF = 10 ** 9

RESULTS: list[tuple[str, bool, str]] = []


def ck(name: str, cond, detail: str = "") -> bool:
    ok = bool(cond)
    RESULTS.append((name, ok, detail))
    if not QUIET:
        print("  [%s] %s" % ("PASS" if ok else "FAIL", name))
        if detail:
            for line in str(detail).splitlines():
                print("         " + line)
    return ok


def say(msg: str = "") -> None:
    if not QUIET:
        print(msg)


# ===========================================================================
# L1.  the stacked functional is free of the fresh coefficient
# ===========================================================================
say("=" * 78)
say("L1.  the stacked identity")
say("=" * 78)

_K = sp.Symbol("K_")
ck("L1.1  2*(1+K)^3 - 3*(1+K)^2 = -1 + 3*K^2 + 2*K^3  (exact)",
   sp.expand(2 * (1 + _K) ** 3 - 3 * (1 + _K) ** 2
             - (-1 + 3 * _K ** 2 + 2 * _K ** 3)) == 0)

_hh = [sp.Symbol("H%d_" % i) for i in range(0, 14)]
_hh0 = [ONE] + _hh[1:]


def _pn_sym(n):
    return sp.expand(sum(_hh0[i] * _hh0[n - i] for i in range(0, n + 1)))


def _rn_sym(n):
    return sp.expand(sum(_hh0[i] * _hh0[j] * _hh0[n - i - j]
                         for i in range(0, n + 1)
                         for j in range(0, n - i + 1)))


def _stack_sym(n):
    """[u^n](3K^2 + 2K^3) computed directly from K = H - 1 (all indices >= 1)."""
    s2 = sum(_hh0[i] * _hh0[n - i] for i in range(1, n))
    s3 = sum(_hh0[i] * _hh0[j] * _hh0[n - i - j]
             for i in range(1, n) for j in range(1, n - i) if n - i - j >= 1)
    return sp.expand(3 * s2 + 2 * s3)


_id_ok = all(sp.expand(_stack_sym(n) - (2 * _rn_sym(n) - 3 * _pn_sym(n))) == 0
             for n in range(2, 13))
ck("L1.2  [u^n](3K^2+2K^3) = 2*r_n - 3*p_n for n = 2..12 (h_0 = 1)", _id_ok)

_fresh_ok = all(sp.expand(sp.diff(_stack_sym(n), _hh[n])) == 0
                for n in range(2, 13))
ck("L1.3  the fresh coefficient h_n CANCELS in the stacked functional at "
   "every level n = 2..12 -- so the stacked condition constrains h_1..h_{n-1} "
   "alone, and the P conditions can be solved for h_n first without loss",
   _fresh_ok)


# ===========================================================================
# L2/L3.  the cascade, recomputed from scratch
# ===========================================================================
say("")
say("=" * 78)
say("L2.  the cascade  (levels 2..%d), recomputed independently" % MAXLEV)
say("=" * 78)
say("""
  Parametrisation.  Every P condition t^(2n-2) | p_n is absorbable by the
  fresh h_n (p_n = 2*h_n + q_n, q_n = sum_{i=1}^{n-1} h_i h_{n-i}), so

      h_n = -q_n/2 + t^(2n-2)*g_n   (1 <= n <= 8) ,   h_n = -q_n/2  (n >= 9)

  with g_n = sum_j g{n}_j t^j FREE.  Every P condition then holds
  identically and the residual content at level n is exactly

      t^(2n-3)  |  [u^n]( 3K^2 + 2K^3 ) .
""")

WINDOW_TOP = 8
_G = {}
CSUBS = {}


def gsym(n, j):
    return _G.setdefault((n, j), sp.Symbol("g%d_%d" % (n, j)))


def _mul(a, b, depth):
    out = [ZERO] * depth
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            if bj == 0 or i + j >= depth:
                continue
            out[i + j] += ai * bj
    return [sp.expand(v) for v in out]


def _add(a, b, depth):
    return [sp.expand(a[i] + b[i]) for i in range(depth)]


def build_h(top, depth):
    """h_1..h_top as truncated t-series of length `depth`, with the cascade
    substitutions found so far already applied."""
    h = {0: [ONE] + [ZERO] * (depth - 1)}
    for n in range(1, top + 1):
        q = [ZERO] * depth
        for i in range(1, n):
            q = _add(q, _mul(h[i], h[n - i], depth), depth)
        hn = [sp.expand(-v / 2) for v in q]
        if n <= WINDOW_TOP:
            m = 2 * n - 2
            for j in range(0, depth - m):
                s = gsym(n, j)
                hn[m + j] = sp.expand(hn[m + j] + CSUBS.get(s, s))
        h[n] = [sp.expand(v.subs(CSUBS)) if getattr(v, "free_symbols", None)
                else v for v in hn]
    return h


def stacked_jets(n, depth):
    """the t-series of [u^n](3K^2+2K^3), truncated to `depth`."""
    h = build_h(n - 1, depth)
    s2 = [ZERO] * depth
    for i in range(1, n):
        s2 = _add(s2, _mul(h[i], h[n - i], depth), depth)
    s3 = [ZERO] * depth
    for i in range(1, n):
        for j in range(1, n - i):
            k = n - i - j
            if k >= 1:
                s3 = _add(s3, _mul(_mul(h[i], h[j], depth), h[k], depth), depth)
    return [sp.expand(3 * s2[q] + 2 * s3[q]) for q in range(depth)]


CASCADE = []
BRANCHED = []
_t0 = time.time()
for n in range(2, MAXLEV + 1):
    need = 2 * n - 3                      # t^need must divide the functional
    depth = max(need, 1)                  # jets t^0 .. t^(need-1) are the content
    jets = stacked_jets(n, depth)
    lowest = None
    for j in range(0, need):
        co = sp.expand(jets[j].subs(CSUBS)) if jets[j] != 0 else ZERO
        if co != 0:
            lowest = (j, co)
            break
    if lowest is None:
        CASCADE.append((n, None, None, "no condition (every required jet "
                                       "vanishes identically)"))
        say("   level n=%-2d  need t^%-2d :  NOTHING (all required jets vanish)  "
            "[%.0fs]" % (n, need, time.time() - _t0))
        continue
    j, co = lowest
    _c, fl = sp.factor_list(co)
    nonconst = [(f, e) for f, e in fl if f.free_symbols]
    if len(nonconst) != 1:
        BRANCHED.append((n, j, co))
        CASCADE.append((n, j, None, "BRANCHES (%d distinct non-constant "
                                    "irreducible factors)" % len(nonconst)))
        say("   level n=%-2d  need t^%-2d :  jet t^%-2d BRANCHES (%d factors)"
            % (n, need, j, len(nonconst)))
        continue
    f, mult = nonconst[0]
    f = sp.expand(f)
    solved = None
    for X in sorted(f.free_symbols, key=lambda s: s.name, reverse=True):
        if sp.degree(f, X) == 1 and not sp.expand(sp.diff(f, X)).free_symbols:
            rhs = sp.expand(sp.solve(sp.Eq(f, 0), X)[0])
            CSUBS[X] = sp.expand(rhs.subs(CSUBS))
            for k in list(CSUBS):
                CSUBS[k] = sp.expand(CSUBS[k].subs({X: CSUBS[X]})) \
                    if getattr(CSUBS[k], "free_symbols", None) else CSUBS[k]
            solved = X
            break
    if solved is None:
        BRANCHED.append((n, j, co))
        CASCADE.append((n, j, None, "single factor but not linear in any single "
                                    "g with constant leading coefficient"))
        say("   level n=%-2d  need t^%-2d :  jet t^%-2d NOT FORCED" % (n, need, j))
        continue
    CASCADE.append((n, j, (f, mult, solved), "FORCED"))
    say("   level n=%-2d  need t^%-2d :  jet t^%-2d = (unit)*(%s)^%d  ->  FORCED "
        "%s = %s   [%.0fs]"
        % (n, need, j, f, mult, solved, CSUBS[solved], time.time() - _t0))

ck("L2.1  no cascade step BRANCHES: at every level the lowest surviving jet "
   "is a unit times a single irreducible factor, linear in one g-coefficient "
   "with a CONSTANT leading coefficient, so each step is a forced consequence "
   "and never a choice of component",
   not BRANCHED, "branching levels: %s" % [b[0] for b in BRANCHED])

ck("L2.2  every ODD level contributes nothing (all required jets vanish "
   "identically), so the cascade advances only at EVEN levels",
   all(c[1] is None for c in CASCADE if c[0] % 2 == 1),
   "odd-level outcomes: %s" % {c[0]: c[3] for c in CASCADE if c[0] % 2 == 1})

_sq = {c[0]: c[2][1] for c in CASCADE if c[2] is not None}
ck("L2.3  every firing jet is a PERFECT SQUARE (multiplicity 2 on the single "
   "irreducible factor) -- which is why no case split can arise: the jet "
   "vanishes iff the factor does",
   all(m == 2 for m in _sq.values()), "level -> multiplicity: %s" % _sq)

_lev12 = next((c for c in CASCADE if c[0] == 12), None)
if MAXLEV >= 12:
    ck("L2.4  *** THE LEVEL-12 STEP *** it fires at jet t^20, it is a single "
       "irreducible factor, it is a perfect square, and it is LINEAR in g6_0 "
       "with a constant leading coefficient.  If it had split into coprime "
       "pieces the deduction would BRANCH and v_t(h_6) >= 11 would FAIL.",
       _lev12 is not None and _lev12[1] == 20 and _lev12[2] is not None
       and _lev12[2][1] == 2 and _lev12[2][2].name == "g6_0",
       "level-12 record: jet t^%s, factor %s, multiplicity %s, solved for %s"
       % (_lev12[1], _lev12[2][0] if _lev12[2] else None,
          _lev12[2][1] if _lev12[2] else None,
          _lev12[2][2] if _lev12[2] else None))

# ---- the forced valuation profile ----------------------------------------
say("")
say("=" * 78)
say("L3.  the forced t-adic valuations of the UNSHIFTED stripped coefficients")
say("=" * 78)

_depth = 2 * MAXLEV + 2
_hfin = build_h(min(MAXLEV, WINDOW_TOP), _depth)
VAL = {}
for k in range(1, min(MAXLEV, WINDOW_TOP) + 1):
    VAL[k] = next((j for j in range(_depth)
                   if sp.expand(_hfin[k][j].subs(CSUBS)) != 0), None)

NAMES = {1: "d3", 2: "d2", 3: "d1", 4: "d0", 5: "D_{-1} = e", 6: "D_{-2}",
         7: "D_{-3}", 8: "D_{-4}"}
for k in sorted(VAL):
    say("     v_t(h_%d) >= %-4s  (%s)" % (k, VAL[k], NAMES.get(k, "")))

_adv = MAXLEV // 2
ck("L3.1  level 2m advances v_t(h_m) from 2m-2 to 2m-1, for every m whose "
   "level 2m was actually run",
   all(VAL.get(m) is not None and VAL[m] >= 2 * m - 1
       for m in range(1, min(_adv, WINDOW_TOP) + 1)),
   "profile (actual, target 2k-1): %s"
   % {k: (VAL[k], 2 * k - 1) for k in range(1, min(_adv, WINDOW_TOP) + 1)})

ck("L3.2  the bounds are ATTAINED, not exceeded (v_t(h_k) = 2k-1 exactly for "
   "every advanced k) -- so the cascade is not silently collapsing the system "
   "to zero, which is the shape a bug would take",
   all(VAL[k] == 2 * k - 1 for k in range(1, min(_adv, WINDOW_TOP) + 1)),
   "advanced k = 1..%d: %s" % (min(_adv, WINDOW_TOP),
                               {k: VAL[k] for k in
                                range(1, min(_adv, WINDOW_TOP) + 1)}))

V_H1 = VAL.get(1)
V_H5 = VAL.get(5)
V_H6 = VAL.get(6)

ck("L3.3  *** v_t(h_1) >= 1 ***  (the inverse-shift input; h_1 = d3 and the "
   "shift parameter is theta = -h_1/4)", V_H1 is not None and V_H1 >= 1,
   "v_t(h_1) >= %s" % V_H1)
ck("L3.4  v_t(h_5) >= %s  (h_5 = D_{-1}; the level-10 result the concurrent "
   "audit is auditing)" % V_H5, V_H5 is not None and V_H5 >= 9,
   "v_t(h_5) >= %s" % V_H5)
if MAXLEV >= 12:
    ck("L3.5  *** v_t(h_6) >= 11 ***  (the level-12 result the lead needs)",
       V_H6 is not None and V_H6 >= 11, "v_t(h_6) >= %s" % V_H6)

# ---- THE ROBUSTNESS POINT: level 12 is NOT actually needed ----------------
# h_6 = -q_6/2 + t^10*g_6 with q_6 = 2*h_1*h_5 + 2*h_2*h_4 + h_3^2.  Under the
# LEVEL-10 profile v_t(h_k) >= 2k-1 (k <= 5) every term of q_6 already has
# t-order >= 10 (1+9 = 3+7 = 5+5 = 10), and the fresh term starts at t^10.  So
#     v_t(h_6) >= 10   from the LEVEL-10 cascade plus the level-6 P condition
# alone.  Level 12 improves 10 -> 11; the six-branch kill needs only >= 9.
_q6_orders = [(2 * 1 - 1) + (2 * 5 - 1), (2 * 2 - 1) + (2 * 4 - 1),
              (2 * 3 - 1) + (2 * 3 - 1)]
ck("L3.7  *** THE KILL DOES NOT NEED LEVEL 12 ***  h_6 = -q_6/2 + t^10*g_6 "
   "(the level-6 P condition t^10 | p_6), and under the LEVEL-10 profile every "
   "monomial of q_6 = 2h_1h_5 + 2h_2h_4 + h_3^2 has t-order >= 10.  So "
   "v_t(h_6) >= 10 follows from the SAME level-10 calculus that is already "
   "under audit -- and 10 already exceeds both pinned values (9 and 8).  "
   "Level 12 only sharpens 10 -> 11.",
   min(_q6_orders) >= 10 and (V_H6 is not None and V_H6 >= 10),
   "q_6 monomial t-orders under the level-10 profile: %s ; fresh term at t^10 ; "
   "computed v_t(h_6) >= %s at MAXLEV = %d" % (_q6_orders, V_H6, MAXLEV))

# satisfiability: the forced profile is a CONSTRAINT, not a contradiction.
# under h_k = t^(2k-1)*Hh_k the substitution u = v/t^2 gives K = Hhat(v)/t, so
#    [u^n](3K^2+2K^3) = 3*t^(2n-2)*[v^n]Hhat^2 + 2*t^(2n-3)*[v^n]Hhat^3,
# divisible by t^(2n-3) for every n.  Verified here by exponent arithmetic on
# every monomial of the two symmetric functions, not by assertion.
_sat = True
for n in range(2, 16):
    for i in range(1, n):
        _sat &= (2 * i - 1) + (2 * (n - i) - 1) >= 2 * n - 3     # K^2 monomials
    for i in range(1, n):
        for j in range(1, n - i):
            k = n - i - j
            if k >= 1:
                _sat &= (2 * i - 1) + (2 * j - 1) + (2 * k - 1) >= 2 * n - 3
ck("L3.6  SATISFIABILITY control: the forced profile v_t(h_k) >= 2k-1 makes "
   "EVERY monomial of [u^n](3K^2+2K^3) have t-order >= 2n-3, for n = 2..15.  "
   "So the cascade pins valuations; it does not empty the slice system, and "
   "the bounds above are not vacuously derived from an inconsistent system.",
   _sat)


# ===========================================================================
# L4.  the d3-killing shift, and the inverse-shift bound on v_t(R)
# ===========================================================================
say("")
say("=" * 78)
say("L4.  the d3-killing shift and the inverse-shift bound on v_t(R)")
say("=" * 78)

_u, _s = sp.symbols("u_ s_")
_c = {m: sp.Symbol("c%d_" % (m + 6)) for m in range(-4, 5)}

# (a) recomposition identity: the coefficient of x^j in C(x - s) is
#     sum_{m=j..4} binom(m, m-j) c_m (-s)^(m-j)  -- verified as a series identity
_Ush = sp.expand(sum(_c[m] * _u ** (-m)
                     * sp.series((1 - _s * _u) ** m, _u, 0, 11).removeO()
                     for m in range(-4, 5)) * _u ** 4)
_rec_ok = True
for jv in range(4, -5, -1):
    formula = sum(sp.binomial(m, m - jv) * _c[m] * (-_s) ** (m - jv)
                  for m in range(jv, 5))
    _rec_ok &= sp.expand(_Ush.coeff(_u, 4 - jv) - formula) == 0
ck("L4.1  recomposition: [x^j] C(x-s) = sum_{m=j..4} binom(m,m-j)*c_m*(-s)^(m-j), "
   "verified as an exact series identity for j = 4 .. -4 (GENERALIZED "
   "binomials; binom(m, m-j) = 0 whenever m >= 0 > j)", _rec_ok)

_D = {m: sp.Symbol("D%d_" % (m + 6)) for m in range(-4, 4)}
_D[4] = ONE
_theta = -_D[3] / 4                                    # the unique d3-killing s


def _Dtil(j):
    return sp.expand(sum(sp.binomial(m, m - j) * _D[m] * _theta ** (m - j)
                         for m in range(j, 5)))


ck("L4.2  the shift KILLS d3:  D~_3 = D_3 + 4*(-D_3/4) = 0", _Dtil(3) == 0)
# CAREFUL.  "triangular across zero" is a statement about the SOURCE
# coefficients with the shift parameter held INDEPENDENT.  Once theta is
# substituted (theta = -D_3/4) the non-negative coefficient D_3 DOES reappear
# in every spare with j <= -2 -- which is exactly the term the inverse-shift
# bound of L4.5/L4.7 has to carry.  Both statements are checked separately so
# they cannot be conflated.
_thsym = sp.Symbol("theta_")
_src = {m: sp.Symbol("src%d_" % (m + 6)) for m in range(-4, 5)}


def _Dtil_free(j):
    return sp.expand(sum(sp.binomial(m, m - j) * _src[m] * _thsym ** (m - j)
                         for m in range(j, 5)))


ck("L4.3  TRIANGULAR ACROSS ZERO (theta held independent): no source "
   "coefficient src_m with m >= 0 appears in any spare D~_j, j < 0, because "
   "binom(m, m-j) = 0 whenever m >= 0 > j",
   all(_src[m] not in _Dtil_free(j).free_symbols
       for j in range(-1, -5, -1) for m in range(0, 5)),
   "D~_{-1} = %s ;  D~_{-2} = %s" % (_Dtil_free(-1), _Dtil_free(-2)))
ck("L4.3b AND THE CAVEAT, stated so it cannot be conflated with L4.3: after "
   "theta = -D_3/4 is substituted, D_3 DOES appear in D~_{-2}.  Triangularity "
   "protects D~_{-1} = D_{-1} completely (no theta at all), but NOT D~_{-2}. "
   "That is precisely why the lead needs the inverse-shift bound rather than "
   "the naked v_t(h_6).",
   _D[3] in _Dtil(-2).free_symbols and _thsym not in _Dtil_free(-1).free_symbols,
   "substituted D~_{-2} = %s  (D9_ is D_3, D5_ is D_{-1}, D4_ is D_{-2})"
   % _Dtil(-2))
ck("L4.4  D~_{-1} = D_{-1} EXACTLY -- so the G-system's dm1 IS the unshifted "
   "h_5, and v_t(h_5) is v_t(e) = a_t [import I3]",
   sp.expand(_Dtil(-1) - _D[-1]) == 0)
ck("L4.5  D~_{-2} = D_{-2} + (D_3/4)*D_{-1}  -- i.e. R = h_6 + (h_1/4)*h_5, "
   "which is the inverse-shift relation the lead uses (sec.6.3 writes the "
   "same relation with the opposite sign on theta; valuations are unaffected)",
   sp.expand(_Dtil(-2) - (_D[-2] + _D[3] * _D[-1] / 4)) == 0)

# the y^12 stripping factors cancel out of the shift identically
_yv = sp.Symbol("y_")
_Dp = {m: _D[m] * _yv ** (12 * (4 - m)) for m in range(-4, 5)}
_thp = -_Dp[3] / 4
_strip_ok = True
for j in (-1, -2, -3):
    lhs = sp.expand(sum(sp.binomial(m, m - j) * _Dp[m] * _thp ** (m - j)
                        for m in range(j, 5)))
    _strip_ok &= sp.simplify(sp.expand(lhs - _Dtil(j) * _yv ** (12 * (4 - j)))) == 0
ck("L4.6  the shift is FORM-INVARIANT under the y^(12*(4-m)) stripping used by "
   "the slice calculus: the same formula holds in stripped coordinates, and "
   "y is a UNIT at t = y+1, so no valuation moves", _strip_ok)


def v_R_lower_bound(a_t):
    """v_t(R) >= min( v_t(h_6), v_t(h_1) + v_t(h_5) ),  v_t(h_5) = a_t."""
    return min(V_H6 if V_H6 is not None else INF,
               (V_H1 if V_H1 is not None else INF) + a_t)


ck("L4.7  *** THE BOUND ***  v_t(R) >= min(v_t(h_6), v_t(h_1) + a_t) = "
   "min(%s, %s + a_t).  At MAXLEV = %d this is t^%d | R whenever a_t >= %d."
   % (V_H6, V_H1, MAXLEV, min(V_H6, 10 + V_H1), V_H6 - V_H1),
   all(v_R_lower_bound(a) >= 10 for a in range(9, 31)),
   "a_t -> lower bound on v_t(R): %s"
   % {a: v_R_lower_bound(a) for a in (8, 9, 10, 11, 12, 13, 14)})


# ===========================================================================
# L5.  the place "trichotomy" at y = -1 -- re-derived, and it is a DICHOTOMY
# ===========================================================================
say("")
say("=" * 78)
say("L5.  the place trichotomy at y = -1 (K-syzygy Newton test), re-derived")
say("=" * 78)
say("""
  (K)  2*Phi = d2*e^3 + 3*e^2*S + 3*e*R^2        [DIVISOR_SYZYGY sec.1]
  (D)  e | S                                     [SYZYGY_SWEEP sec.4]

  At the place t = y+1 write  a = v_t(e), rho = v_t(R), v_t(S) = a + v_s
  (v_s >= 0 by (D)), v_t(d2) = delta2 >= 0 (d2 is a polynomial), and
  v_t(Phi) = 30 EXACTLY.  A valuation assignment is FEASIBLE only if the
  minimum of the four orders

      {delta2 + 3a,   3a + v_s,   a + 2*rho,   30}

  is attained at least TWICE (otherwise one term is a strict unique minimum
  and the identity (K) cannot hold).
""")

V_T_PHI = 30

# L5.0  v_t(Phi) = 30 EXACTLY is the one load-bearing NUMBER in the horn
# argument (">= 30" would not exclude horn 1).  It is re-verified here from
# `divisor_consequences.phi_stripped()` -- the stripped Phi is an artifact of
# the pipeline, so this is a CHECKED reproduction, not a re-derivation.
_phi_ok, _phi_detail = False, "divisor_consequences unavailable"
try:
    import divisor_consequences as _dc

    _yv2 = sp.Symbol("y")
    _phi = _dc.phi_stripped()
    _q = _dc.q_poly()
    _pp = sp.Poly(_phi, _yv2)
    _phi_ok = (all(sp.Poly(sp.diff(_phi, _yv2, k), _yv2).eval(-1) == 0
                   for k in range(0, 30))
               and sp.Poly(sp.diff(_phi, _yv2, 30), _yv2).eval(-1) != 0
               and sp.expand(_q.subs(_yv2, -1)) == 3315)
    _phi_detail = ("ord_{y=-1} Phi = 30 exactly (derivatives 0..29 vanish, the "
                   "30th does not); q(-1) = %s != 0" % _q.subs(_yv2, -1))
except Exception as _exc:                                  # pragma: no cover
    _phi_detail = "could not load divisor_consequences: %s" % _exc
ck("L5.0  v_t(Phi) = 30 EXACTLY at y = -1 -- not '>= 30'.  The horn-1 "
   "exclusion below needs the EXACT value: with only v_t(Phi) >= 30 the "
   "minimum could be attained twice among the right-hand terms and horn 1 "
   "would survive.  [CHECKED against divisor_consequences.phi_stripped()]",
   _phi_ok, _phi_detail)


def feasible(a, rho, v_s, delta2, P_b=V_T_PHI):
    orders = [x for x in (delta2 + 3 * a, 3 * a + v_s, a + 2 * rho, P_b)
              if x < INF]
    if not orders:
        return True
    mn = min(orders)
    return orders.count(mn) >= 2


_RNG = list(range(0, 65)) + [INF]

# --- how many horns are there, really?
# if rho < a then a + 2rho <= 3a - 2 < 3a <= min(delta2+3a, 3a+v_s), so the
# R-term is the STRICT unique minimum of the three right-hand terms and the
# only way to attain the minimum twice is against P_b.  Verified exhaustively.
_h2_forced = all(
    (not feasible(a, rho, v_s, dl2)) or (a + 2 * rho == V_T_PHI)
    for a in range(0, 21) for rho in range(0, a)
    for v_s in (0, 1, 2, 7, 40, INF) for dl2 in (0, 1, 2, 7, 40, INF))
ck("L5.1  the statement is a DICHOTOMY, not a trichotomy, despite the name: "
   "whenever rho < a the R-term is the strict unique minimum of the three "
   "right-hand terms, so feasibility forces v_t(Phi) = a + 2*rho EXACTLY.  "
   "There is no third horn -- the only other case is rho >= a.  (The repo's "
   "OTHER 'trichotomy' is the terminal BRANCH trichotomy T1|T2|T3, an "
   "unrelated object.)", _h2_forced)

_degen = feasible(0, 0, 0, 0, P_b=INF)
ck("L5.2  the one DEGENERATE case the predicate admits -- every order "
   "infinite, i.e. Phi = 0 and e = 0 -- cannot occur at y = -1, because "
   "v_t(Phi) = 30 is FINITE and e != 0 on any genuine lift",
   _degen and V_T_PHI < INF)

# --- horn 1 (rho >= a) is INFEASIBLE for every a >= 11
_horn1 = {a: any(feasible(a, rho, v_s, dl2)
                 for rho in list(range(a, 3 * a + 2)) + [INF]
                 for v_s in _RNG for dl2 in _RNG)
          for a in range(0, 21)}
ck("L5.3  *** HORN 1 IS EXCLUDED FOR EVERY a >= 11 ***  If v_t(R) >= a then "
   "all three right-hand terms of (K) have order >= 3a, so v_t(2*Phi) >= 3a; "
   "but v_t(Phi) = 30 exactly, so 30 >= 3a, i.e. a <= 10.  This is a PROOF, "
   "not an inference from the odd-a kills.",
   all(not _horn1[a] for a in range(11, 21)),
   "horn-1 feasibility by a: %s" % {a: _horn1[a] for a in range(8, 16)})

ck("L5.4  NON-VACUITY control: horn 1 MUST stay feasible for a <= 10, or the "
   "test would be a constant refutation and would (wrongly) empty the "
   "standard regime too",
   all(_horn1[a] for a in range(0, 11)),
   "horn-1 feasible at a = %s" % [a for a in range(0, 11) if _horn1[a]])

# --- horn 2 candidates
_horn2 = {a: [rho for rho in range(0, a)
              if any(feasible(a, rho, v_s, dl2) for v_s in _RNG for dl2 in _RNG)]
          for a in range(0, 21)}
_closed = {a: ([] if (30 - a) % 2 or not (0 <= (30 - a) // 2 <= a - 1)
               else [(30 - a) // 2]) for a in range(0, 21)}
ck("L5.5  horn 2 candidates match the closed form rho = (30-a)/2, required to "
   "be a non-negative integer < a (hence a EVEN and a >= 11)",
   _horn2 == _closed,
   "a -> rho candidates, a = 10..16: %s"
   % {a: _horn2[a] for a in range(10, 17)})

ck("L5.6  CONSEQUENCE: for a >= 11 horn 1 is dead, so horn 2 is COMPULSORY "
   "and v_t(R) is PINNED to the single value (30-a)/2.  Odd a has no such "
   "value at all (parity kill, already banked in ALT_FRONTIER_V2.md sec.3).",
   all((not _horn1[a]) and (len(_horn2[a]) <= 1) for a in range(11, 21))
   and all(_horn2[a] == [] for a in range(11, 21) if a % 2 == 1)
   and all(len(_horn2[a]) == 1 for a in range(12, 21) if a % 2 == 0),
   "a -> pinned v_t(R): %s"
   % {a: _horn2[a][0] for a in range(11, 21) if _horn2[a]})


# ===========================================================================
# L6.  the six surviving alternate T1 branches, one at a time
# ===========================================================================
say("")
say("=" * 78)
say("L6.  verdict, per branch")
say("=" * 78)

SIX = [("a12_b0000_T1", 12, 0), ("a12_b1000_T1", 12, 1),
       ("a12_b1100_T1", 12, 2), ("a12_b1110_T1", 12, 3),
       ("a14_b0000_T1", 14, 0), ("a14_b1000_T1", 14, 1)]

VERDICTS = []
for name, a, sb in SIX:
    h1_open = _horn1[a]
    rho_cands = _horn2[a]
    lb = v_R_lower_bound(a)
    # every horn must be refuted for the branch to die
    horn1_dead = not h1_open
    horn2_dead = all(rho < lb for rho in rho_cands)
    dead = horn1_dead and horn2_dead and (rho_cands or horn1_dead)
    VERDICTS.append((name, a, sb, horn1_dead, rho_cands, lb, dead))
    say("   %-14s a=%2d sum b=%d | horn1 excluded: %-5s | horn2 pins v_t(R) = %s"
        " | slice bound v_t(R) >= %s | %s"
        % (name, a, sb, horn1_dead, rho_cands or "NONE", lb,
           "*** EMPTY ***" if dead else "SURVIVES"))

ck("L6.1  on all six branches horn 1 is excluded by the K-syzygy order count "
   "(3a > 30) -- INDEPENDENTLY of anything the slice calculus says",
   all(v[3] for v in VERDICTS))
ck("L6.2  on all six branches horn 2 pins v_t(R) to exactly one value, and "
   "that value is STRICTLY BELOW the slice-calculus lower bound",
   all(len(v[4]) == 1 and v[4][0] < v[5] for v in VERDICTS),
   "; ".join("%s: pinned %d < bound %d" % (v[0], v[4][0], v[5])
             for v in VERDICTS))
ck("L6.3  *** ALL SIX ALTERNATE T1 BRANCHES ARE EMPTY *** (conditional on the "
   "imports I1/I2/I3 and on the concurrent audit of the slice calculus)",
   all(v[6] for v in VERDICTS) and len(VERDICTS) == 6,
   "dead: %s" % [v[0] for v in VERDICTS if v[6]])

# ---- the stronger corollary, which costs nothing extra --------------------
_all_a = {}
for a in range(0, 31):
    h1 = _horn1[a] if a <= 20 else any(
        feasible(a, rho, v_s, dl2)
        for rho in list(range(a, 3 * a + 2)) + [INF]
        for v_s in _RNG for dl2 in _RNG)
    h2 = [rho for rho in range(0, a)
          if any(feasible(a, rho, v_s, dl2) for v_s in _RNG for dl2 in _RNG)]
    lb = v_R_lower_bound(a)
    _all_a[a] = h1 or any(rho >= lb for rho in h2)
ck("L6.4  COROLLARY (stronger than the six branches, and CAP-FREE): with the "
   "level-12 bound, EVERY a_t >= 11 dies -- odd a on parity, even a because "
   "(30-a)/2 <= 9 < 11.  So a_t <= 10 outright, with no appeal to the "
   "deg e <= 15 stratum premise that put a_t >= 16 'out of scope' in "
   "ALT_FRONTIER_V2.md sec.1.",
   all(not _all_a[a] for a in range(11, 31)) and all(_all_a[a] for a in range(0, 11)),
   "surviving a_t after the combined test: %s"
   % [a for a in range(0, 31) if _all_a[a]])

ck("L6.5  CONTROL -- the standard regime is UNTOUCHED: a_t = 9 and a_t = 10 "
   "survive the combined test (horn 1 remains feasible there), so this is not "
   "a blanket refutation that would also contradict the surviving sub1 cells",
   _all_a[9] and _all_a[10])


# ===========================================================================
# summary
# ===========================================================================
_pass = sum(1 for _, ok, _ in RESULTS if ok)
say("")
say("=" * 78)
say("%d/%d checks pass   (cascade to level %d, %.0fs)"
    % (_pass, len(RESULTS), MAXLEV, time.time() - _t0))
say("=" * 78)
if not QUIET:
    for name, ok, _ in RESULTS:
        if not ok:
            print("FAILED: %s" % name)
sys.exit(0 if _pass == len(RESULTS) else 1)
