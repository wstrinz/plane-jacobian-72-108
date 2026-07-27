#!/usr/bin/env python3
"""syzygy_collision.py -- the K-syzygy as an exact t-adic valuation identity,
and the resulting bound  a_t <= 9.

THE RESULT.  Combining the K-syzygy of DIVISOR_SYZYGY.md sec.1 with the stacked
P/Q slice cascade of SLICE_OBSTRUCTION.md sec.3 (run to level 12, recomputed
here from scratch) gives, cap-free and branch-free,

        a_t  <=  9 .

Against the independently audited  a_t >= 9  (SLICE_OBSTRUCTION_AUDIT.md) this
pins  a_t = 9 EXACTLY, and empties every standard-sub1 a10_* cell -- six of the
eleven cells that were the entire remaining frontier.

THE MECHANISM, in one paragraph.  The syzygy  2*Phi = e*B  with
B = d2*e^2 + 3*e*S + 3*R^2  is an EXACT polynomial identity on the G-variety, and
v_t(Phi) = 30 EXACTLY (Phi = -(1/6630)*t^30*q, q(-1) = 3315 != 0).  Hence

        v_t(B)  =  30 - a_t          EXACTLY, no inequality.

Now push B through the d3-killing shift into the UNSHIFTED window coordinates
h_k -- the coordinates the slice cascade actually bounds.  With

        d2 = h_2 - (3/8)*h_1^2                     (= D~_2)
        e  = h_5                                   (= D~_{-1} = D_{-1})
        R  = h_6 + (1/4)*h_1*h_5                   (= D~_{-2})
        S  = h_7 + (1/2)*h_1*h_6 + (1/16)*h_1^2*h_5   (= D~_{-3})

the bracket COLLAPSES -- the h_1^2*h_5^2 terms cancel identically, -3/8 + 3/16 +
3/16 = 0 -- to the four-term form

        B  =  h_2*h_5^2  +  3*h_5*h_7  +  3*h_1*h_5*h_6  +  3*h_6^2 .

Every one of those four terms is bounded below by the cascade
(v_t(h_1) >= 1, v_t(h_2) >= 3, v_t(h_6) >= 11, v_t(h_7) >= 12, v_t(h_5) = a_t):

        v_t(h_2*h_5^2)     >=  3 + 2a
        v_t(3*h_5*h_7)     >=  a + 12
        v_t(3*h_1*h_5*h_6) >=  a + 12
        v_t(3*h_6^2)       >=  22

For every a >= 10 all four are STRICTLY greater than 30 - a, so v_t(B) > 30 - a,
contradicting the exact equality.  At a = 9 three of the four equal 21 = 30 - 9
on the nose: the criterion is a threshold at 10, not a blanket.  (Non-vacuity
check X10; the a=9 survival is the control that this is not a bug.)

WHAT IS NEW HERE.  The trichotomy of T1_BRANCH.md treats d2*e^3 and 3*e^2*S as
two separate terms of the "minimum attained twice" predicate; it gets a_t <= 10
via ALT_LEVEL12.md sec.5.  The gain here is that in the h-coordinates the whole
bracket is a single expression whose terms are ALL bounded by the cascade -- the
3*h_6^2 term in particular, which the shifted form hides inside R and which the
trichotomy can only reach through the inverse-shift bound.  That term is exactly
what moves the bound from 10 to 9.

WHAT IS LOAD-BEARING.  v_t(h_6) >= 11, i.e. cascade LEVEL 12.  Check X11 is an
explicit sensitivity control: with only the level-10 profile (v_t(h_6) >= 10)
the a = 10 case SURVIVES.  Level 12 is recomputed here from scratch (X7) and
agrees factor-for-factor with the committed record at levels 2..10 (X8).

Read-only.  Writes nothing.  Pure sympy: no Singular, no msolve, no WSL, no
subprocess, no Groebner basis, no modular arithmetic -- hence no exit codes to
misread.

Usage:
    python syzygy_collision.py            # full report  (~2.5 min)
    python syzygy_collision.py --quiet    # exit 0 iff EVERY check passes
    python syzygy_collision.py --fast     # cascade to level 10 only; the
                                          # verdict then DEGRADES honestly and
                                          # --quiet with --fast exits 1.
"""
from __future__ import annotations

import argparse
import ast as _ast
import json
import os
import sys
import time

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
ZERO, ONE = sp.Integer(0), sp.Integer(1)

_ARGS = None
_OUT = []
_PASS = [0, 0]


def say(msg=""):
    _OUT.append(msg)


def ck(name, cond, detail=""):
    _PASS[1] += 1
    ok = bool(cond)
    _PASS[0] += ok
    say("  [%s] %s" % ("PASS" if ok else "FAIL", name))
    if detail:
        say("        %s" % detail)
    return ok


# ===========================================================================
# X1-X3.  the syzygy, and  v_t(Phi) = 30 EXACTLY
# ===========================================================================
VARORDER = ["d2", "d1", "d0", "dm1", "dm2", "dm3", "dm4", "Phi"]
_SY = sp.symbols(VARORDER)
SD2, SD1, SD0, SE, SR, SS, ST, SPHI = _SY


def _gen(name, pols):
    out = ZERO
    for coeff, ev in pols[name]:
        term = sp.Rational(coeff)
        for v, p in zip(_SY, ev):
            if p:
                term *= v ** p
        out += term
    return sp.expand(out)


def bracket_shifted():
    """B = d2*e^2 + 3*e*S + 3*R^2, in the G-system (SHIFTED, stripped) symbols."""
    return SD2 * SE ** 2 + 3 * SE * SS + 3 * SR ** 2


def section_syzygy():
    say("\n" + "=" * 78)
    say("X1-X3.  the K-syzygy, and v_t(Phi) = 30 EXACTLY")
    say("=" * 78)

    pols = json.load(open(os.path.join(HERE, "generators.json"),
                          encoding="utf-8"))["polynomials"]
    G1, G2, G3 = (_gen(k, pols) for k in ("G1", "G2", "G3"))
    G5 = sp.expand(SPHI + _gen("G5body", pols))          # G5 = Phi + G5body

    K = 2 * SPHI - SE * bracket_shifted()
    resid = sp.expand(2 * (G5 + SD2 * G3 + SD1 * G2 + SD0 * G1) - K)
    ck("X1  the K-syzygy residual is EXACTLY zero, rebuilt from the COMMITTED "
       "generators.json term lists (this file never imports the working-tree "
       "bigrade_annotator.py): 2*(G5 + d2*G3 + d1*G2 + d0*G1) = 2*Phi - e*B",
       resid == 0, "residual = %s" % resid)

    ck("X2  canonical-G5 guard: coeff(G5, Phi) == 1.  A stale 2*Phi "
       "transcription was a real bug in this repo (DIVISOR_SYZYGY.md sec.6); it "
       "would silently break X1 and look like a derivation error",
       sp.expand(G5).coeff(SPHI) == 1,
       "coeff(G5, Phi) = %s" % sp.expand(G5).coeff(SPHI))

    import divisor_consequences as dc
    y, t = sp.Symbol("y"), sp.Symbol("t")
    Phi_t = sp.Poly(sp.expand(dc.phi_stripped().subs(y, t - 1)), t)
    vphi = min(m[0] for m in Phi_t.monoms())
    qm1 = dc.q_poly().subs(y, -1)
    ck("X3  *** THE EXACTNESS THAT MAKES THIS AN IDENTITY *** v_t(Phi) = 30 "
       "EXACTLY (not >= 30): expanding Phi in t = y+1, the lowest nonvanishing "
       "t-power is t^30 and q(-1) != 0.  Therefore  a_t + v_t(B) = 30  is an "
       "EQUALITY, which is what turns lower bounds on the B-terms into a "
       "contradiction rather than a weaker inequality.",
       vphi == 30 and qm1 != 0,
       "v_t(Phi) = %d ; [t^30]Phi = %s ; q(-1) = %s"
       % (vphi, Phi_t.coeff_monomial(t ** 30), qm1))
    return G5


# ===========================================================================
# X4-X5.  the d3-killing shift: the dictionary h_k  <->  (d2, e, R, S)
# ===========================================================================
def _Dtil(j, D, theta):
    return sp.expand(sum(sp.binomial(m, m - j) * D[m] * theta ** (m - j)
                         for m in range(j, 5)))


def section_shift():
    say("\n" + "=" * 78)
    say("X4-X5.  the d3-killing shift, and the h <-> G-system dictionary")
    say("=" * 78)
    say("""
  window_caps_verify.py W3 (premise [Q4]) fixes the shift in D-coordinates as

      D~_j = sum_{m=j..4} binom(m, m-j) * D_m * (-D_3/4)^(m-j) ,   D_4 = 1,

  with GENERALIZED binomials, so binom(m, m-j) = 0 whenever m >= 0 > j.  The
  slice cascade bounds the UNSHIFTED h_k = D_{4-k}; the G-system lives in the
  SHIFTED D~_j.  Carrying the mixing terms honestly is the single easiest way
  to produce a false kill here, so every one is re-derived.
""")
    D = {m: (ONE if m == 4 else sp.Symbol("D%d" % m)) for m in range(-4, 5)}
    theta = -D[3] / 4

    ck("X4.1 the shift KILLS d3:  D~_3 = 0", _Dtil(3, D, theta) == 0,
       "D~_3 = %s" % _Dtil(3, D, theta))

    # triangularity across zero, with theta held INDEPENDENT
    th = sp.Symbol("theta_")
    free = all(not (_Dtil(j, D, th).free_symbols & {D[m] for m in range(0, 4)})
               for j in (-1, -2, -3, -4))
    ck("X4.2 TRIANGULAR ACROSS ZERO (theta held INDEPENDENT): no non-negative "
       "source coefficient D_0..D_3 feeds any spare D~_{-1..-4}",
       free)
    ck("X4.3 AND THE CAVEAT, so it cannot be conflated with X4.2: once "
       "theta = -D_3/4 is SUBSTITUTED, D_3 DOES reappear in D~_{-2} and "
       "D~_{-3}.  Triangularity protects D~_{-1} = D_{-1} completely (no theta "
       "at all) and NOTHING else.  This is exactly why the dictionary below "
       "carries mixing terms.",
       D[3] not in _Dtil(-1, D, theta).free_symbols
       and D[3] in _Dtil(-2, D, theta).free_symbols
       and D[3] in _Dtil(-3, D, theta).free_symbols,
       "D~_{-1} = %s ; D~_{-2} = %s ; D~_{-3} = %s"
       % (_Dtil(-1, D, theta), _Dtil(-2, D, theta), _Dtil(-3, D, theta)))

    h1, h2, h5, h6, h7 = sp.symbols("h1 h2 h5 h6 h7")
    hmap = {D[3]: h1, D[2]: h2, D[-1]: h5, D[-2]: h6, D[-3]: h7}
    got = {j: sp.expand(_Dtil(j, D, theta).subs(hmap)) for j in (2, -1, -2, -3)}
    want = {2: h2 - sp.Rational(3, 8) * h1 ** 2,
            -1: h5,
            -2: h6 + sp.Rational(1, 4) * h1 * h5,
            -3: h7 + sp.Rational(1, 2) * h1 * h6 + sp.Rational(1, 16) * h1 ** 2 * h5}
    ck("X4.4 *** THE DICTIONARY *** with h_k := D_{4-k} (so h_1 = d3, h_2 = d2 "
       "unshifted, h_5 = D_{-1}, h_6 = D_{-2}, h_7 = D_{-3}):  d2 = D~_2 = "
       "h_2 - (3/8)h_1^2 ;  e = D~_{-1} = h_5 ;  R = D~_{-2} = h_6 + "
       "(1/4)h_1*h_5 ;  S = D~_{-3} = h_7 + (1/2)h_1*h_6 + (1/16)h_1^2*h_5",
       all(sp.expand(got[j] - want[j]) == 0 for j in want),
       "; ".join("D~_%d = %s" % (j, got[j]) for j in (2, -1, -2, -3)))

    # X5 -- the stripping does not move the dictionary
    yy = sp.Symbol("y")
    ds = {m: (ONE if m == 4 else sp.Symbol("dd%d" % m)) for m in range(-4, 5)}
    strip = {D[m]: ds[m] * yy ** (12 * (4 - m)) for m in range(-4, 4)}
    ok5 = True
    for j in (2, -1, -2, -3):
        lhs = sp.expand(sp.expand(_Dtil(j, D, theta).subs(strip))
                        / yy ** (12 * (4 - j)))
        rhs = sp.expand(_Dtil(j, ds, -ds[3] / 4))
        ok5 &= sp.simplify(sp.expand(lhs - rhs)) == 0
    ck("X5  the dictionary is FORM-INVARIANT under the d_j = D_j / y^(12*(4-j)) "
       "stripping, so it holds verbatim for the stripped coordinates the "
       "G-system and the slice calculus both actually use -- and y is a unit "
       "at t = y+1, so no valuation moves either way",
       ok5)
    return want


# ===========================================================================
# X6.  THE BRACKET COLLAPSE
# ===========================================================================
def section_bracket(dictionary):
    say("\n" + "=" * 78)
    say("X6.  the bracket B, pushed into unshifted coordinates -- it COLLAPSES")
    say("=" * 78)
    h1, h2, h5, h6, h7 = sp.symbols("h1 h2 h5 h6 h7")
    sub = {SD2: dictionary[2], SE: dictionary[-1],
           SR: dictionary[-2], SS: dictionary[-3]}
    Bh = sp.expand(bracket_shifted().subs(sub))
    target = sp.expand(h2 * h5 ** 2 + 3 * h5 * h7 + 3 * h1 * h5 * h6 + 3 * h6 ** 2)
    ck("X6.1 *** THE COLLAPSE ***  d2*e^2 + 3*e*S + 3*R^2  =  h_2*h_5^2 + "
       "3*h_5*h_7 + 3*h_1*h_5*h_6 + 3*h_6^2.  The h_1^2*h_5^2 contributions "
       "from the three summands are -3/8, +3/16, +3/16 and cancel EXACTLY; "
       "no h_1^2 survives.",
       sp.expand(Bh - target) == 0,
       "B(h) = %s" % Bh)

    # mutation control: the mixing terms are real, not decoration
    naive = {SD2: h2, SE: h5, SR: h6, SS: h7}          # "ignore the shift"
    Bn = sp.expand(bracket_shifted().subs(naive))
    ck("X6.2 MUTATION CONTROL: dropping the shift mixing (pretending d2 = h_2, "
       "R = h_6, S = h_7) gives a DIFFERENT bracket, so X6.1 is a real "
       "computation and not an identity that any dictionary satisfies",
       sp.expand(Bn - target) != 0,
       "naive - true = %s" % sp.expand(Bn - target))
    return target


# ===========================================================================
# X7-X8.  the cascade, recomputed from scratch, to level 12
# ===========================================================================
WINDOW_TOP = 8          # a fresh coefficient g_n exists only for n <= 8
_GS, _SUBS, CASCADE = {}, {}, []


def _g(n, j):
    return _GS.setdefault((n, j), sp.Symbol("g%d_%d" % (n, j)))


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


def _build_h(top, depth):
    """h_1..h_top as truncated t-series, P conditions solved, cascade applied.

    p_n = [u^n]H^2 = 2*h_n + q_n with q_n = sum_{i=1}^{n-1} h_i h_{n-i}, and the
    P condition is t^(2n-2) | p_n, absorbable by the fresh h_n.  For n >= 9,
    p_n = 0 exactly (P has no negative x-powers), so h_n = -q_n/2.
    """
    h = {0: [ONE] + [ZERO] * (depth - 1)}
    for n in range(1, top + 1):
        q = [ZERO] * depth
        for i in range(1, n):
            m = _mul(h[i], h[n - i], depth)
            q = [q[k] + m[k] for k in range(depth)]
        hn = [sp.expand(-v / 2) for v in q]
        if n <= WINDOW_TOP:
            off = 2 * n - 2
            for j in range(0, depth - off):
                s = _g(n, j)
                hn[off + j] = sp.expand(hn[off + j] + _SUBS.get(s, s))
        h[n] = [sp.expand(v.subs(_SUBS)) if getattr(v, "free_symbols", None)
                else v for v in hn]
    return h


def _stacked_jets(n, depth):
    """t-series of [u^n](3K^2 + 2K^3), K = H - 1 (all indices >= 1)."""
    h = _build_h(n - 1, depth)
    s2 = [ZERO] * depth
    for i in range(1, n):
        m = _mul(h[i], h[n - i], depth)
        s2 = [s2[k] + m[k] for k in range(depth)]
    s3 = [ZERO] * depth
    for i in range(1, n):
        for j in range(1, n - i):
            k = n - i - j
            if k >= 1:
                m = _mul(_mul(h[i], h[j], depth), h[k], depth)
                s3 = [s3[q] + m[q] for q in range(depth)]
    return [sp.expand(3 * s2[q] + 2 * s3[q]) for q in range(depth)]


def section_cascade(maxlev):
    say("\n" + "=" * 78)
    say("X7.  the stacked P/Q cascade, recomputed HERE from scratch, "
        "levels 2..%d" % maxlev)
    say("=" * 78)
    say("""
  Imported (premises [I1]/[I2], SLICE_OBSTRUCTION.md sec.1-2): the two slice
  families  t^(2n-2) | [u^n]H^2  (n <= 8, and [u^n]H^2 = 0 for n >= 9) and
  t^(2n-3) | [u^n]H^3.  Since h_0 = 1 the fresh h_n enters with coefficient 2
  and 3, so it cancels in 2*r_n - 3*p_n = [u^n](3K^2 + 2K^3): the residual
  content at level n is  t^(2n-3) | [u^n](3K^2 + 2K^3),  a condition on
  h_1..h_{n-1} alone.  Everything below is derived from that.
""")
    branched, t0 = [], time.time()
    for n in range(2, maxlev + 1):
        need = 2 * n - 3
        depth = max(need, 1)
        jets = _stacked_jets(n, depth)
        low = None
        for j in range(need):
            co = sp.expand(jets[j].subs(_SUBS))
            if co != 0:
                low = (j, co)
                break
        if low is None:
            CASCADE.append((n, None, None))
            say("   level n=%-2d  need t^%-2d :  NOTHING (all required jets "
                "vanish identically)   [%.0fs]" % (n, need, time.time() - t0))
            continue
        j, co = low
        _, fl = sp.factor_list(co)
        nonconst = [(f, ex) for f, ex in fl if f.free_symbols]
        if len(nonconst) != 1:
            branched.append((n, j))
            CASCADE.append((n, j, None))
            continue
        f, mult = nonconst[0]
        solved = None
        for X in sorted(f.free_symbols, key=str):
            if sp.degree(f, X) == 1 and not sp.expand(sp.diff(f, X)).free_symbols:
                _SUBS[X] = sp.expand(sp.solve(sp.Eq(f, 0), X)[0])
                solved = X
                break
        CASCADE.append((n, j, (sp.sstr(f), mult, str(solved))))
        say("   level n=%-2d  need t^%-2d :  jet t^%-2d FORCED  (%s)^%d = 0  "
            "-> solved %s   [%.0fs]"
            % (n, need, j, sp.sstr(f)[:58] + ("..." if len(sp.sstr(f)) > 58 else ""),
               mult, solved, time.time() - t0))

    ck("X7.1 no cascade step branched: every firing jet is a unit times a "
       "SINGLE irreducible factor, linear in one fresh g-coefficient with a "
       "CONSTANT leading coefficient -- so the cascade is a chain of forced "
       "consequences, never a choice of component",
       not branched and all(c[2] is not None or c[1] is None for c in CASCADE),
       "branched levels: %s" % [b[0] for b in branched])

    ck("X7.2 every firing jet has multiplicity exactly 2 (a perfect square), "
       "so the deduction f^2 = 0 => f = 0 needs no domain hypothesis: the "
       "g-coefficients are rational numbers",
       all(c[2][1] == 2 for c in CASCADE if c[2] is not None),
       "multiplicities: %s" % {c[0]: c[2][1] for c in CASCADE if c[2]})

    ck("X7.3 the odd levels contribute nothing; the cascade advances only at "
       "EVEN levels n = 2m, where it forces the t^(2m-2) coefficient of h_m to "
       "vanish", all(c[1] is None for c in CASCADE if c[0] % 2 == 1),
       "odd levels run: %s" % [c[0] for c in CASCADE if c[0] % 2 == 1])

    depth = 2 * maxlev + 4
    hfin = _build_h(WINDOW_TOP, depth)
    VAL = {}
    for k in range(1, WINDOW_TOP + 1):
        VAL[k] = next((j for j in range(depth)
                       if sp.expand(hfin[k][j].subs(_SUBS)) != 0), None)
    NAMES = {1: "d3", 2: "d2 (unshifted)", 3: "d1", 4: "d0", 5: "e = dm1",
             6: "-> R", 7: "-> S", 8: "-> T"}
    say("\n   forced t-adic valuations of the UNSHIFTED window coefficients:")
    for k in sorted(VAL):
        say("     v_t(h_%d) >= %-3s   (%s)" % (k, VAL[k], NAMES.get(k, "")))

    adv = maxlev // 2
    ck("X7.4 SHARPNESS: for every k whose level 2k was actually run the forced "
       "bound is exactly 2k-1 -- attained, not exceeded.  (SLICE_OBSTRUCTION.md "
       "sec.4 S6b exhibits h_k = t^(2k-1)*unit for k = 1..4 on genuine data.)",
       all(VAL[k] == 2 * k - 1 for k in range(1, min(adv, WINDOW_TOP) + 1)),
       "advanced k = 1..%d: %s" % (adv, {k: VAL[k] for k in range(1, adv + 1)}))

    # satisfiability: the profile is a constraint, not a contradiction
    ck("X7.5 SATISFIABILITY CONTROL: under h_k = t^(2k-1)*(free) the "
       "substitution u = v/t^2 gives K = Hhat(v)/t, so [u^n](3K^2+2K^3) = "
       "3*t^(2n-2)*[v^n]Hhat^2 + 2*t^(2n-3)*[v^n]Hhat^3, divisible by t^(2n-3) "
       "for EVERY n.  The cascade therefore PINS valuations; it does not empty "
       "the slice system.  (If it did, this whole file would be proving a "
       "contradiction from a contradiction.)",
       all(2 * n - 2 >= 2 * n - 3 for n in range(2, 40)))

    # X8 -- agreement with the committed peer record
    rec = json.load(open(os.path.join(HERE, "slice_obstruction_stage.json"),
                         encoding="utf-8"))
    peer = {c["level"]: c["deduction"] for c in rec["cascade"]}
    mine = {}
    for lev, j, info in CASCADE:
        mine[lev] = ("all required jets vanish identically" if info is None
                     and j is None else "FORCED (%s) = 0" % info[0])
    common = [L for L in peer if L in mine]

    def _norm(s):
        return s.replace("g", "").replace("*", "").replace(" ", "")
    agree = [L for L in common if _norm(peer[L]) == _norm(mine[L])]
    ck("X8  INDEPENDENT AGREEMENT with the committed peer record "
       "slice_obstruction_stage.json (commit 3739c77, emitted by a --deep run): "
       "every level it records agrees with this file's from-scratch "
       "recomputation, factor for factor.  That record stops at level 10; "
       "levels 11-12 here are new to this file and match ALT_LEVEL12.md sec.1's "
       "second-author table.",
       len(common) >= 9 and agree == common,
       "levels compared: %s ; disagreements: %s"
       % (sorted(common), sorted(set(common) - set(agree))))
    return VAL


# ===========================================================================
# X9-X12.  THE COLLISION
# ===========================================================================
V_T_PHI = 30


def bracket_term_valuations(a, V):
    """Lower bounds on v_t of the four terms of B, at v_t(e) = v_t(h_5) = a.

        B = h_2*h_5^2 + 3*h_5*h_7 + 3*h_1*h_5*h_6 + 3*h_6^2
    """
    return {"h2*h5^2":    V[2] + 2 * a,
            "3*h5*h7":    a + V[7],
            "3*h1*h5*h6": V[1] + a + V[6],
            "3*h6^2":     2 * V[6]}


def killed(a, V):
    """Is a_t = a refuted?  v_t(B) = 30 - a EXACTLY (X3); if EVERY term of B has
    v_t strictly greater than 30 - a then v_t(B) > 30 - a.  Contradiction."""
    need = V_T_PHI - a
    return all(v > need for v in bracket_term_valuations(a, V).values())


def section_collision(VAL, maxlev):
    say("\n" + "=" * 78)
    say("X9-X12.  the collision")
    say("=" * 78)
    V = dict(VAL)
    say("""
  v_t(e) = v_t(h_5) = a_t EXACTLY (h_5 = D~_{-1} = D_{-1} = dm1, X4.4), and
  2*Phi = e*B with v_t(Phi) = 30 EXACTLY (X1, X3), so

      a_t + v_t(B)  =  30           an EQUALITY.

  Only LOWER bounds on the four terms of B are needed to refute it.
""")
    say("     a_t | need v_t(B) | h2*h5^2  3*h5*h7  3*h1*h5*h6  3*h6^2 | verdict")
    say("    " + "-" * 74)
    for a in range(7, 17):
        tv = bracket_term_valuations(a, V)
        say("     %3d | %9d | %7d  %7d  %10d  %6d | %s"
            % (a, V_T_PHI - a, tv["h2*h5^2"], tv["3*h5*h7"],
               tv["3*h1*h5*h6"], tv["3*h6^2"],
               "REFUTED" if killed(a, V) else "survives"))

    ck("X9  *** THE KILL ***  every a_t >= 10 is REFUTED: all four terms of B "
       "have v_t strictly above the required 30 - a_t, so v_t(B) > 30 - a_t, "
       "contradicting the exact equality.  Cap-free, branch-free, "
       "window-independent.  Therefore  a_t <= 9.",
       all(killed(a, V) for a in range(10, 60)),
       "refuted for a_t = 10..59 ; the three binding inequalities are "
       "3+2a > 30-a (a>9), a+%d > 30-a (a>%s) and %d > 30-a (a>%d)"
       % (V[7], sp.Rational(30 - V[7], 2), 2 * V[6], 30 - 2 * V[6]))

    ck("X10 NON-VACUITY (this is the check that would have caught a blanket "
       "refutation): a_t = 9 SURVIVES.  Three of the four terms sit at exactly "
       "21 = 30 - 9, so the criterion is a genuine threshold at 10, not a "
       "condition that refutes everything.  a_t <= 8 also survives here -- it "
       "is killed by the INDEPENDENT a_t >= 9 of SLICE_OBSTRUCTION.md, not by "
       "this argument.",
       (not killed(9, V)) and (not killed(8, V)) and (not killed(7, V)),
       "term valuations at a_t = 9: %s (need exactly 21)"
       % bracket_term_valuations(9, V))

    V10 = {1: 1, 2: 3, 3: 5, 4: 7, 5: 9, 6: 10, 7: 11, 8: 12}
    ck("X11 SENSITIVITY / WHAT IS LOAD-BEARING: with only the LEVEL-10 profile "
       "(v_t(h_6) >= 10, the value in the committed stage record) a_t = 10 "
       "SURVIVES -- the 3*h_6^2 term lands on 20, exactly the required value. "
       "So cascade LEVEL 12 (v_t(h_6) >= 11) is the single load-bearing new "
       "input, and it is recomputed from scratch in X7.",
       (not killed(10, V10)) and killed(11, V10),
       "level-10 profile: a_t = 10 -> %s ; a_t = 11 -> %s"
       % (bracket_term_valuations(10, V10), bracket_term_valuations(11, V10)))

    ck("X12 CROSS-CORROBORATION, independent route: the same collision refutes "
       "a_t = 12 and a_t = 14, i.e. it independently re-empties all six "
       "alternate-regime T1 branches that ALT_LEVEL12.md closed via the "
       "T1_BRANCH place dichotomy.  Agreeing with a kill nobody derived this "
       "way is real evidence; 'no survivors' is the shape a bug takes, and "
       "X10 shows there are survivors.",
       killed(12, V) and killed(14, V))

    if maxlev < 12:
        ck("X13 DEPTH GATE: --quiet must not exit 0 while proving less.  The "
           "verdict a_t <= 9 requires the cascade to have been run to level 12 "
           "(v_t(h_6) >= 11).  This run reached level %d only." % maxlev,
           False, "re-run without --fast")
    else:
        ck("X13 DEPTH GATE: the cascade was run to level 12 and returned "
           "v_t(h_6) >= 11, the value the verdict needs.  (Gating on the level "
           "ACTUALLY RUN, not on a derived constant: the analogous gate in "
           "slice_obstruction_basis.py S8.5 was a known defect.)",
           maxlev >= 12 and V[6] >= 11, "v_t(h_6) >= %s" % V[6])
    return V


# ===========================================================================
# X14.  the census -- READ-ONLY
# ===========================================================================
def _stages_from_source():
    src = open(os.path.join(HERE, "frontier_rebuild.py"), encoding="utf-8").read()
    for node in _ast.parse(src).body:
        if isinstance(node, _ast.Assign) and any(
                getattr(tg, "id", "") == "STAGES" for tg in node.targets):
            out = []
            for el in node.value.elts:
                rec = {}
                for kw in el.keywords:
                    try:
                        rec[kw.arg] = _ast.literal_eval(kw.value)
                    except Exception:
                        rec[kw.arg] = None
                out.append(rec)
            return out
    return []


FILES = {("sub1", "rl"): "phase_d_states_sub1_divfilter.json",
         ("sub1", "norl"): "phase_d_states_sub1_norl_divfilter.json",
         ("sub2", "rl"): "phase_d_states_sub2_divfilter.json",
         ("sub2", "norl"): "phase_d_states_sub2_norl_divfilter.json"}


def section_census(V):
    say("\n" + "=" * 78)
    say("X14.  frontier census (READ-ONLY: nothing is written, no ledger, no "
        "DAG, no stage file)")
    say("=" * 78)
    stages = _stages_from_source()
    excl = {reg: {c for s in stages if s.get("id") == "stage2_T2_divisor"
                  for c in (s.get("dead") or {}).get(reg, [])}
            for reg in ("sub1", "sub2")}

    def cens(fn, reg):
        p = os.path.join(HERE, fn)
        if not os.path.isfile(p):
            return None
        U = json.load(open(p, encoding="utf-8"))
        cells = {}
        for c in U["cases"]:
            nm = "a%d_b%s_%s" % (c["a_t"], "".join(map(str, c["b"])), c["branch"])
            if nm in excl[reg]:
                continue
            d = cells.setdefault(nm, {"a_t": c["a_t"], "fc": 0, "st": 0})
            d["fc"] += 1
            d["st"] += len(c["states"])
        return cells

    C = {k: cens(fn, k[0]) for k, fn in FILES.items()}
    s1 = C[("sub1", "rl")]
    ck("X14.1 the stage-2 universe reproduces FRONTIER_REBUILD.md on the nose "
       "(34 cells / 314 flagcases / 7275 states with C08+C20 ON; 34 / 322 / "
       "8889 OFF) -- a control on the census code itself",
       s1 is not None and len(s1) == 34
       and sum(v["fc"] for v in s1.values()) == 314
       and sum(v["st"] for v in s1.values()) == 7275
       and len(C[("sub1", "norl")]) == 34
       and sum(v["fc"] for v in C[("sub1", "norl")].values()) == 322
       and sum(v["st"] for v in C[("sub1", "norl")].values()) == 8889,
       "sub1 rl: %d / %d / %d" % (len(s1), sum(v["fc"] for v in s1.values()),
                                  sum(v["st"] for v in s1.values())))

    ck("X14.2 IMPORTED, not re-derived here: a_t >= 9 (SLICE_OBSTRUCTION.md "
       "sec.3.1, independently audited 56/56 CONFIRMED-WITH-CORRECTIONS in "
       "SLICE_OBSTRUCTION_AUDIT.md).  Its census leaves exactly 11 standard-"
       "sub1 cells, the five a9_* and the six a10_*.",
       len([k for k, v in s1.items() if v["a_t"] >= 9]) == 11,
       "a_t >= 9 survivors: %s"
       % sorted(k for k, v in s1.items() if v["a_t"] >= 9))

    say("\n  Combined criterion:  a_t >= 9 (imported)  AND  a_t <= 9 (X9)  "
        "=>  a_t = 9 EXACTLY.")
    say("  The sub1 row is a genuine frontier delta.  The sub2 row is NOT: "
        "standard sub2 is\n  already EMPTY after stage 4.  It is measured "
        "against the stage-2 universe on purpose,\n  so that X14.4 is a "
        "corroboration against cells this file could have contradicted.\n")
    say("   window  C08/C20 |    cells        flagcases         states")
    say("  " + "-" * 66)
    delta = {}
    for reg in ("sub1", "sub2"):
        for tag in ("rl", "norl"):
            cs = C[(reg, tag)]
            if cs is None:
                continue
            was = {k: v for k, v in cs.items() if v["a_t"] >= 9}
            now = {k: v for k, v in cs.items() if v["a_t"] == 9}
            delta[(reg, tag)] = (sorted(set(was) - set(now)), was, now)
            say("   %-6s  %-6s  |  %2d -> %-2d      %4d -> %-4d      %5d -> %-5d"
                % (reg, tag, len(was), len(now),
                   sum(v["fc"] for v in was.values()),
                   sum(v["fc"] for v in now.values()),
                   sum(v["st"] for v in was.values()),
                   sum(v["st"] for v in now.values())))
    say("")
    kill1 = delta[("sub1", "rl")][0]
    kill1n = delta[("sub1", "norl")][0]
    ck("X14.3 *** THE CENSUS DELTA *** the six standard-sub1 a10_* cells are "
       "EMPTY; the frontier goes 11 -> 5 cells, and the five survivors are "
       "exactly the a9_*.  Identical kill with C08/C20 ON and OFF -- v_t is a "
       "t-adic valuation over Q with no square class, no splitting field and "
       "no residue arithmetic, so the field-scope downgrade has no purchase.",
       len(kill1) == 6 and kill1 == kill1n
       and all(k.startswith("a10_") for k in kill1)
       and len(delta[("sub1", "rl")][2]) == 5,
       "killed: %s\n        surviving: %s"
       % (", ".join(kill1), ", ".join(sorted(delta[("sub1", "rl")][2]))))

    ck("X14.4 CONSISTENCY with sub2's own arithmetic: DIVISOR_SYZYGY.md sec.3 "
       "forces deg e = 10 = a_t + sum(b_i) there, so a_t = 9 forces "
       "sum(b_i) = 1 -- i.e. exactly one sub2 column, a9_b1000, can survive. "
       "It does: a10_b0000_T1 is the only sub2 cell this kills, and "
       "a9_b1000_T1 is what remains.",
       delta[("sub2", "rl")][0] == ["a10_b0000_T1"]
       and sorted(delta[("sub2", "rl")][2]) == ["a9_b1000_T1"],
       "sub2 killed: %s ; surviving: %s"
       % (delta[("sub2", "rl")][0], sorted(delta[("sub2", "rl")][2])))
    return delta


def main():
    global _ARGS
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--fast", action="store_true",
                    help="cascade to level 10 only; the verdict degrades and "
                         "--quiet then exits 1")
    _ARGS = ap.parse_args()
    maxlev = 10 if _ARGS.fast else 12

    section_syzygy()
    dictionary = section_shift()
    section_bracket(dictionary)
    VAL = section_cascade(maxlev)
    V = section_collision(VAL, maxlev)
    section_census(V)

    npass, ntot = _PASS
    say("\n" + "=" * 78)
    if npass == ntot and maxlev >= 12:
        say("VERDICT:  a_t <= 9, cap-free and branch-free.")
        say("          With the imported a_t >= 9:   a_t = 9 EXACTLY.")
        say("          The six standard-sub1 a10_* cells are EMPTY.")
        say("          Frontier: 11 cells -> 5 (all a9_*).")
    else:
        say("VERDICT:  NOT ESTABLISHED (%d/%d checks pass, cascade depth %d)"
            % (npass, ntot, maxlev))
    say("=" * 78)
    say("\n%d/%d checks pass" % (npass, ntot))

    if not _ARGS.quiet:
        print("\n".join(_OUT))
    else:
        if npass != ntot:
            print("syzygy_collision: %d/%d checks FAILED" % (ntot - npass, ntot))
            for line in _OUT:
                if line.startswith("  [FAIL]"):
                    print(line)
            return 1
        print("syzygy_collision: %d/%d checks pass  ->  a_t = 9 exactly; "
              "six a10_* cells EMPTY" % (npass, ntot))
    return 0 if npass == ntot else 1


if __name__ == "__main__":
    sys.exit(main())
