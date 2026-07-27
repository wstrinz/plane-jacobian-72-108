#!/usr/bin/env python3
"""
prop43_audit.py -- mechanized audit of the *exhaustiveness* of GGHV22 Proposition 4.3.

Target
------
Guccione-Guccione-Horruitiner-Valqui, "Increasing the degree of a possible
counterexample to the Jacobian Conjecture from 100 to 108" (arXiv:2204.14178),
Proposition 4.3 [Case (8,28)].  Local source copy: paper_src/2204.14178.tex,
Proposition statement at L1000-1008, proof at L1009-1290.

Both of the (72,108) exclusion programs -- ours (repo HEAD 1e2d99b) and Helali's
(doi:10.5281/zenodo.21479814, adjudicated in HELALI_ADJUDICATION.md) -- are
conditional on Prop 4.3 being an EXHAUSTIVE reduction: on a (72,108)
counterexample necessarily landing in one of its two Newton configurations.
This script mechanizes the finitely-checkable part of that claim.

What is mechanized
------------------
The proof of Prop 4.3 is Newton-polygon combinatorics over a lattice, driven by
four upstream results.  Their statements are transcribed here from the local
arXiv sources (see CITATIONS below), and every step of the proof that is a
finite lattice computation is re-run from those statements.

Design rule: no check is allowed to be trivially true.  Two checks are
"calibration" checks -- they re-derive tables that the paper itself PRINTS, so a
transcription error in the upstream statement shows up as a mismatch:

  * check_calibration_ggv2_algorithm  reproduces the two sets that GGHV22 states
    for the *neighbouring* Proposition 4.2 (case (9,24)) from the algorithm the
    paper prints there.
  * check_calibration_prop41_table    reproduces, row for row and column for
    column, the divisibility table that GGHV22 PRINTS inside Proposition 4.1
    (case (9,27)) -- 9 candidate points, both numeric columns, both survivors.

Everything downstream then reuses exactly those calibrated rules on the (8,28)
data, where the paper prints no table at all.

CITATIONS (local arXiv sources under paper_src/)
------------------------------------------------
  GGV1 = arXiv:1401.1784   Cor 7.4  (paper_src/1401.1784_GGV1.tex L4238-4270)
                           Prop 8.2 (paper_src/1401.1784_GGV1.tex L5320-5358)
                           Prop 7.3 (L4192), v/st/en definitions (L322-380)
  GGV2 = arXiv:1605.09430  Prop 3.12  (used via GGV6's restatement)
  GGV6 = arXiv:1708.09367  "b = 2" proposition (L316-405), cited by GGHV22 as
                           "GGV6 Proposition 2.5"; the local arXiv version has
                           coarser section numbering than the published
                           Pro Mathematica version, so it is matched BY CONTENT
                           (it is the unique result in GGV6 whose hypothesis is
                           "set b := 2", which is exactly en(R) = (7,2)).
  GGV5 = arXiv:1708.07936  case table L1831-1832: A_0 = (8,28), A_1 = (11/4,7),
                           (m,n) = (3,2), max deg = 108.

Usage:  python prop43_audit.py [--quiet]
Exit 0 iff every check passes.
"""

from __future__ import annotations

import json
import os
import sys
from fractions import Fraction
from math import gcd

HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------------------
# lattice / direction helpers (GGV1 L322-380: v_{rho,sigma}(i,j) = rho*i+sigma*j,
# v_{rho,sigma}(P) is the MAX over Supp(P), ell is the leading form, st/en are the
# endpoints of the face, ordered counterclockwise)
# ---------------------------------------------------------------------------


def v(d, p):
    return d[0] * p[0] + d[1] * p[1]


def prim(w):
    g = gcd(abs(w[0]), abs(w[1]))
    return (w[0] // g, w[1] // g)


def normal_pos_rho(w):
    """Primitive (rho,sigma) with rho>0 and v_{rho,sigma}(w)=0 (lower-side form)."""
    n = prim((w[1], -w[0]))
    if n[0] < 0 or (n[0] == 0 and n[1] < 0):
        n = (-n[0], -n[1])
    return n


def normal_pos_sum(w):
    """Primitive (rho,sigma) with v(w)=0 and rho+sigma>0 (GGV1's V_{>0} form)."""
    n = prim((w[1], -w[0]))
    if n[0] + n[1] < 0:
        n = (-n[0], -n[1])
    return n


def cross(a, b):
    return a[0] * b[1] - a[1] * b[0]


def hull(points):
    """Convex hull vertices (counterclockwise), exact integer arithmetic."""
    pts = sorted(set(points))
    if len(pts) <= 2:
        return pts

    def build(seq):
        out = []
        for p in seq:
            while len(out) >= 2 and cross(
                (out[-1][0] - out[-2][0], out[-1][1] - out[-2][1]),
                (p[0] - out[-1][0], p[1] - out[-1][1]),
            ) <= 0:
                out.pop()
            out.append(p)
        return out

    lower = build(pts)
    upper = build(reversed(pts))
    return lower[:-1] + upper[:-1]


# ---------------------------------------------------------------------------
# Laurent polynomial arithmetic in K[x, x^-1, y] over Q (for the phi_3 step)
# ---------------------------------------------------------------------------


def pmul(a, b):
    out = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            k = (ka[0] + kb[0], ka[1] + kb[1])
            out[k] = out.get(k, Fraction(0)) + va * vb
    return {k: c for k, c in out.items() if c != 0}


def padd(a, b):
    out = dict(a)
    for k, c in b.items():
        out[k] = out.get(k, Fraction(0)) + c
    return {k: c for k, c in out.items() if c != 0}


def ppow(a, n):
    out = {(0, 0): Fraction(1)}
    for _ in range(n):
        out = pmul(out, a)
    return out


def psubst_y(p, shift):
    """y -> y + shift, where shift is a Laurent polynomial (here alpha*x^-4)."""
    ybase = padd({(0, 1): Fraction(1)}, shift)
    out = {}
    for (i, j), c in p.items():
        term = pmul({(i, 0): c}, ppow(ybase, j))
        out = padd(out, term)
    return out


# ---------------------------------------------------------------------------
# UPSTREAM RULE 1 -- GGV2 Prop 3.12, in the algorithmic form GGHV22 PRINTS
# (paper_src/2204.14178.tex L698-719, "Algorithm: Possible starting points").
# Necessary conditions on st_{rho,sigma}(R) given en_{rho,sigma}(R) = (a/l, b),
# under the "R is not a monomial times a power of a linear form" hypothesis
# (GGV6 eq. "no potencia de un lineal", paper_src/1708.09367.tex L286-291).
# ---------------------------------------------------------------------------


def possible_starting_points(a, b, l=1):
    out = []
    for d in range(0, b):
        c_lo = (d * a) // b + 1
        c_hi = l * d + a - b * l - 1
        for c in range(c_lo, c_hi + 1):
            N1 = gcd(abs(a - c), abs(b - d))
            N2 = gcd(abs(c), abs(d))
            rho, sigma = normal_pos_rho((a - c, l * (b - d)))
            num = rho * a + sigma * b * l
            if num <= 0:
                # v_{rho,sigma}(R) > 0 is a standing hypothesis of GGV2 Prop 3.12
                # (num = 0 is the terminal face, handled separately).
                continue
            den = l * (rho + sigma)
            g = gcd(abs(num), abs(den))
            s = abs(num) // g if g else 0
            ok = (d > 0 and N2 % s == 0 if s else False) or (s <= N1)
            if ok:
                out.append(((c, d), (rho, sigma), s, N1, N2))
    return out


# ---------------------------------------------------------------------------
# UPSTREAM RULE 2 -- GGV6 "b = 2" proposition, criterion (4)
# (paper_src/1708.09367.tex L316-345):
#   equivalent to the existence of R,G with en(R) = (a/l, 2) and the split
#   ("no potencia de un lineal") hypothesis, is:
#       exists Delta in N with l < Delta < a/2 and (a - 2*Delta) | (Delta - l);
#   and then (rho,sigma) ~ (l, -Delta).
# ---------------------------------------------------------------------------


def ggv6_b2_deltas(a, l=1):
    out = []
    for Delta in range(l + 1, (a + 1) // 2 + 1):
        if 2 * Delta >= a:
            continue
        if (Delta - l) % (a - 2 * Delta) == 0:
            out.append(Delta)
    return out


# ---------------------------------------------------------------------------
# UPSTREAM RULE 3 -- the excluded ("non-split") shape.
# If R = lambda * x^u * h(z)^j with h linear in z = x^Delta*y, then with
# en(R) = (A, B) we need j = B (the y-degree of R), st(R) = (A - Delta*B, 0),
# and L^{(1)}-integrality forces rho | l = 1, i.e. (rho,sigma) = (1,-Delta).
# v_{rho,sigma}(R) > 0 forces A - Delta*B > 0.
# ---------------------------------------------------------------------------


def nonsplit_deltas(en, l=1):
    A, B = en
    return [D for D in range(l + 1, A + 1) if A - D * B > 0]


# ---------------------------------------------------------------------------
# UPSTREAM RULE 4 -- GGV1 Prop 8.2 (paper_src/1401.1784_GGV1.tex L5320-5358).
# (1) aligned case: (a',b') in Z x N_0 with
#        v_{rho1,sigma1}(a',b') < v_{rho1,sigma1}(a,b)   and   a*b' - b*a' > 0.
# (2) non-aligned case: exists k in N with (k+1)*b < a and
#        {en(P), en(Q)} = {(-k,0), (k+1,1)}.
# The extra filter GGHV22 applies inside Prop 4.1 (L540-560): the corner of the
# GGV1-Thm-2.6 element F other than (1,1) must satisfy
#        (1,1) + c*(a-a', b-b')/g = (p/q)*(a,b),  g = gcd(a-a', b-b'),
# and evaluating v_{-b,a} on both sides gives
#        (a - b) * g + c * (b*a' - a*b') = 0,   hence   (a*b' - b*a') | (a-b)*g.
# Points on the diagonal Z(1,1) are excluded (dir is undefined there).
# ---------------------------------------------------------------------------


def prop82_candidates(corner, d1):
    """Aligned-case candidate opposite vertices (GGV1 Prop 8.2(1) inequalities).

    b' >= 0 and a' in Z with v_{d1}(a',b') < v_{d1}(a,b) and a*b' - b*a' > 0.
    Both constraints are linear in (a',b') and the feasible region is a bounded
    wedge, so a generous finite box is exhaustive."""
    a, b = corner
    lim = v(d1, corner)
    R = 8 * (abs(a) + abs(b)) + 16
    out = []
    for bp in range(0, R):
        for ap in range(-R, R):
            if v(d1, (ap, bp)) < lim and a * bp - b * ap > 0:
                out.append((ap, bp))
    return sorted(set(out), key=lambda p: (p[1], p[0]))


def prop82_table(corner, d1):
    """The full filter: returns rows (point, a*b'-b*a', (a-b)*gcd, survives)."""
    a, b = corner
    rows = []
    for (ap, bp) in prop82_candidates(corner, d1):
        if ap == bp:  # diagonal Z(1,1): dir undefined, discarded in GGHV22 L556
            rows.append(((ap, bp), None, None, False, "diagonal"))
            continue
        lhs = a * bp - b * ap
        g = gcd(abs(a - ap), abs(b - bp))
        rhs = (a - b) * g
        rows.append(((ap, bp), lhs, rhs, rhs % lhs == 0 if lhs else False, ""))
    return rows


# ===========================================================================
# CHECKS
# ===========================================================================

CHECKS = []


def check(name):
    def deco(fn):
        CHECKS.append((name, fn))
        return fn

    return deco


# --- calibration -----------------------------------------------------------


@check("CAL-1  GGV2/Prop-3.12 algorithm reproduces GGHV22's own Prop 4.2 output")
def _cal1():
    """GGHV22 L695: 'The following algorithm with l=1, a=8 and b=3 shows that
    Pred_P(1,0) in {(1,-2),(2,-5)}'.  GGHV22 L727: 'st_{1,-2}(R) in {(2,0),(4,1)}'
    (after separately excluding (6,2) at L723-726)."""
    res = possible_starting_points(8, 3, 1)
    dirs = sorted({r[1] for r in res})
    pts12 = sorted({r[0] for r in res if r[1] == (1, -2)})
    assert dirs == [(1, -3), (1, -2), (2, -5)] or dirs == [(1, -2), (2, -5)], dirs
    assert set(dirs) == {(1, -2), (2, -5)}, f"directions {dirs} != paper's {{(1,-2),(2,-5)}}"
    assert pts12 == [(2, 0), (4, 1), (6, 2)], pts12
    return f"directions={sorted(set(dirs))} (paper: {{(1,-2),(2,-5)}}); (1,-2)-points={pts12} (paper keeps {{(2,0),(4,1)}} after excluding (6,2))"


@check("CAL-2  GGV1/Prop-8.2 filter reproduces GGHV22's PRINTED table in Prop 4.1")
def _cal2():
    """GGHV22 L544-562 prints, for corner (21,8) and direction (-1,3), a 9-row
    table with columns 21b'-8a' and 13*gcd(21-a',8-b'), survivors (5,2),(13,5)."""
    printed = {
        (-2, 0): (16, 13),
        (-1, 0): (8, 26),
        (2, 1): (5, 13),
        (4, 2): (10, 13),
        (5, 2): (2, 26),
        (7, 3): (7, 13),
        (10, 4): (4, 13),
        (13, 5): (1, 13),
    }
    rows = prop82_table((21, 8), (-1, 3))
    got_pts = [r[0] for r in rows]
    expect_pts = sorted(list(printed) + [(1, 1)], key=lambda p: (p[1], p[0]))
    assert got_pts == expect_pts, f"candidate set {got_pts} != paper's {expect_pts}"
    for (pt, lhs, rhs, ok, note) in rows:
        if pt == (1, 1):
            assert note == "diagonal"
            continue
        assert (lhs, rhs) == printed[pt], f"row {pt}: got {(lhs,rhs)}, paper prints {printed[pt]}"
    surv = sorted(r[0] for r in rows if r[3])
    assert surv == [(5, 2), (13, 5)], surv
    return f"9/9 candidate points, 16/16 printed table entries, survivors {surv} == paper's"


@check("CAL-3  GGV6 'b=2' criterion is consistent with the GGV2 algorithm for all a<=80")
def _cal3():
    """Two independently-transcribed upstream statements: GGV6's b=2 proposition
    (an EQUIVALENCE, paper_src/1708.09367.tex L316-345 criterion (4)) and the
    GGV2-Prop-3.12 algorithm GGHV22 prints at L698-719 (NECESSARY conditions only
    -- GGV6 L296-299 states explicitly that sufficiency is unproven).  So the
    criterion's Delta-set must be CONTAINED in the algorithm's for every a.  A
    transcription slip in either direction breaks containment."""
    bad, tight = [], []
    l = 1
    for a in range(5, 81):
        crit = set(ggv6_b2_deltas(a, l))
        alg = {(-d[1]) for (_, d, _, _, _) in possible_starting_points(a, 2, l)}
        if not crit <= alg:
            bad.append((a, sorted(crit), sorted(alg)))
        if crit == alg and crit:
            tight.append(a)
    assert not bad, f"containment fails at: {bad[:6]}"
    a7 = (set(ggv6_b2_deltas(7, 1)),
          {(-d[1]) for (_, d, _, _, _) in possible_starting_points(7, 2, 1)})
    assert a7[0] == a7[1] == {3}, a7
    return (f"76/76 values of a: criterion(4) subset of algorithm; the two coincide for "
            f"a in {tight[:12]}{'...' if len(tight)>12 else ''} -- in particular for the case at hand, "
            f"a=7: both give Delta={{3}}")


# --- Prop 4.3, first half: the a)/b)/c) trichotomy -------------------------

EN_R = (7, 2)  # GGHV22 L1017: 'Thus en_{rho,sigma}(R) = (7,2)'


@check("P43-1  the split branch at en(R)=(7,2) forces Delta=3, i.e. (1,-3)")
def _p1():
    d = ggv6_b2_deltas(7, 1)
    assert d == [3], d
    alg = sorted({dd for (_, dd, _, _, _) in possible_starting_points(7, 2, 1)})
    assert alg == [(1, -3)], alg
    st = sorted({p for (p, dd, _, _, _) in possible_starting_points(7, 2, 1)})
    assert st == [(1, 0), (4, 1)], st
    return f"Delta in {d} -> direction (1,-3); admissible st(R) in {st} (d=0 form and its normalised d=1 form)"


@check("P43-2  the non-split branch at en(R)=(7,2) allows exactly Delta in {2,3}")
def _p2():
    d = nonsplit_deltas(EN_R)
    assert d == [2, 3], d
    return "R = x^(7-2D)*(x^D y - lam)^2 with 7-2D>0 and D>1  =>  D in {2,3}"


@check("P43-3  union == GGHV22's claim 'Pred_P(1,0) in {(1,-2),(1,-3)}' (L1019)")
def _p3():
    union = sorted({(1, -d) for d in set(ggv6_b2_deltas(7, 1)) | set(nonsplit_deltas(EN_R))},
                   key=lambda t: -t[1])
    assert union == [(1, -2), (1, -3)], union
    return f"{union} == paper's {{(1,-2),(1,-3)}}  [split contributes only (1,-3); (1,-2) is non-split only]"


@check("P43-4  after the (1,-2) shift the next direction is in {(1,-3),(2,-7)} (L1073)")
def _p4():
    """A direction (1,-D) strictly below (1,-2) needs D>2; v_{1,-D}(7,2)=7-2D>0
    needs D<3.5.  So D=3 is the only live edge; otherwise v=0, which (since
    (0,0) in N(P)) forces the terminal face dir((7,2)) = (2,-7) with st=(0,0)."""
    live = [d for d in nonsplit_deltas(EN_R) + ggv6_b2_deltas(7, 1) if d > 2]
    assert sorted(set(live)) == [3], live
    term = normal_pos_rho(EN_R)
    assert term == (2, -7), term
    return f"live edge directions below (1,-2): {{(1,-3)}}; terminal face direction dir((7,2)) = {term} -> paper's {{(1,-3),(2,-7)}}"


@check("P43-5  the (1,-3) edge carries at most 2 distinct linear factors (L1077)")
def _p5():
    """en(R)-st(R) has y-difference 2 in every admissible case, so the z-degree of
    R is 2: at most two distinct roots.  st=(1,0) gives R=x(z-a1)(z-a2); st=(4,1)
    gives R=x^4 y(z-a) = x*z*(z-a), i.e. the same shape with a1=0."""
    shapes = []
    for (st, d, _, _, _) in possible_starting_points(7, 2, 1):
        shapes.append((st, EN_R[1] - st[1]))
    assert all(deg <= 2 for _, deg in shapes), shapes
    assert sorted(shapes) == [((1, 0), 2), ((4, 1), 1)], shapes
    return "z-degree of R is 2 in the d=0 form and 1+monomial in the d=1 form: 'one or two different linear factors', never three"


@check("P43-6  case c): the corner below (16,4) is forced, giving (0,0) via (1,-4)")
def _p6():
    """After the split shift, st_{1,-3}(P) = 4m(4,1), so the next R' has
    en(R')=(4,1).  Every direction below (1,-3) with v>0 is impossible, hence the
    terminal face dir((4,1)) = (1,-4) with st = (0,0)."""
    en2 = (4, 1)
    cand = set(nonsplit_deltas(en2)) | set(ggv6_b2_deltas(4, 1) if en2[1] == 2 else [])
    alg = {(-d[1]) for (_, d, _, _, _) in possible_starting_points(4, 1, 1)}
    live = sorted({d for d in (cand | alg) if d > 3})
    assert live == [], live
    term = normal_pos_rho(en2)
    assert term == (1, -4), term
    return f"no live direction below (1,-3) at en(R')=(4,1); terminal dir((4,1)) = {term}, st=(0,0)  ->  lower side {{(-3,0),(0,0),(16,4),(28,8)}}"


@check("P43-7  the trichotomy a)/b)/c) is COMPLETE (full leaf enumeration)")
def _p7():
    """Enumerate every root-to-terminal path of the direction recursion from
    en(R)=(7,2), using RULE 2 (split) + RULE 3 (non-split) + the terminal rule.
    Each leaf yields a lower-side corner set at the 1/(m,n) scale."""
    leaves = set()

    def walk(en, prev_delta, neg_corner, extra):
        # terminal option: v=0 face, st=(0,0)
        base = [(0, 0), (28, 8)] + extra
        if neg_corner is not None:
            base = [(-neg_corner, 0)] + base
        leaves.add(tuple(sorted(base)))
        # live options
        for D in sorted(set(nonsplit_deltas(en)) | (set(ggv6_b2_deltas(en[0], 1)) if en[1] == 2 else set())):
            if prev_delta is not None and D <= prev_delta:
                continue
            split_ok = (en[1] == 2 and D in ggv6_b2_deltas(en[0], 1))
            nonsplit_ok = D in nonsplit_deltas(en)
            if nonsplit_ok:
                # single linear factor: the shift collapses the whole edge
                walk(en, D, D, extra)
            if split_ok:
                # two distinct factors: shift normalises st(R) to (en_x - D, 1)
                st = (en[0] - D, 1)
                walk(st, D, D, extra + [(4 * st[0], 4 * st[1])])

    walk(EN_R, None, None, [])
    paper = {
        tuple(sorted([(-2, 0), (0, 0), (28, 8)])): "a)",
        tuple(sorted([(-3, 0), (0, 0), (28, 8)])): "b)",
        tuple(sorted([(-3, 0), (0, 0), (16, 4), (28, 8)])): "c)",
    }
    degenerate = tuple(sorted([(0, 0), (28, 8)]))
    unexplained = {L for L in leaves if L not in paper and L != degenerate}
    assert not unexplained, f"leaf not covered by a)/b)/c): {unexplained}"
    assert set(paper) <= leaves, f"paper case not reached: {set(paper)-leaves}"
    return (f"{len(leaves)} leaves: a),b),c) all reached, no uncovered leaf. "
            f"One extra degenerate leaf {degenerate} (no shift applied, hence no negative corner) "
            f"is subsumed: the negative corner is washed out by phi_3 before it is ever used.")


@check("P43-8  the shifts y->y+lam*x^-D (D=2,3) leave the (-1,4) edge intact; D=4 does not")
def _p8():
    """v_{-1,4}(y) = 4.  A shift by lam*x^-D perturbs ell_{-1,4} iff
    v_{-1,4}(x^-D) = D >= 4.  So phi (D=2,3) preserves the corner m(0,1) and the
    edge {(0,1),(28,8)}, while phi_3 (D=4) is exactly the one that moves it."""
    got = {D: v((-1, 4), (-D, 0)) for D in (2, 3, 4)}
    assert got[2] < 4 and got[3] < 4 and got[4] == 4, got
    return f"v_(-1,4)(x^-D) = {got}, v_(-1,4)(y) = 4"


# --- Prop 4.3, second half --------------------------------------------------


@check("P43-9  phi_3 on y*(x^4 y - alpha)^7 reduces {(0,1),(28,8)} to {(24,7),(28,8)}")
def _p9():
    """GGHV22 L1132.  Exact Laurent arithmetic, alpha = 1, m = 1,2,3."""
    res = {}
    for m in (1, 2, 3):
        R = pmul({(0, 1): Fraction(1)}, ppow(padd({(4, 1): Fraction(1)}, {(0, 0): Fraction(-1)}), 7))
        Rm = ppow(R, m)
        out = psubst_y(Rm, {(-4, 0): Fraction(1)})
        vs = hull(list(out.keys()))
        res[m] = sorted(vs)
        assert sorted(vs) == sorted([(24 * m, 7 * m), (28 * m, 8 * m)]), (m, sorted(vs))
    return f"m=1,2,3 all reduce to m*{{(24,7),(28,8)}}: {res[1]}, {res[2]}, {res[3]}"


@check("P43-10 the Prop-8.2 filter at (24,7) yields exactly GGHV22's list (L1136)")
def _p10():
    """GGHV22 L1132 says only 'as in Proposition 4.1, one can analyze the possibilities
    ... and obtain Succ_P(-1,4)=Succ_Q(-1,4)=(-2,7)', and later lists
    (a,b) in {(24,7),(17,5),(10,3),(3,1)}.  We re-run the calibrated filter."""
    rows = prop82_table((24, 7), (-1, 4))
    cands = [r[0] for r in rows]
    surv = sorted((r[0] for r in rows if r[3]), key=lambda p: -p[1])
    assert len(cands) == 12, cands
    assert surv == [(17, 5), (10, 3), (3, 1)], surv
    full = [(24, 7)] + surv
    assert full == [(24, 7), (17, 5), (10, 3), (3, 1)], full
    # all four collinear along (7,2) -> a single edge, direction (-2,7)
    diffs = {prim((full[i][0] - full[i + 1][0], full[i][1] - full[i + 1][1])) for i in range(3)}
    assert diffs == {(7, 2)}, diffs
    d = normal_pos_sum((7, 2))
    assert d == (-2, 7), d
    return (f"{len(cands)} candidates -> survivors {surv}; with (24,7) gives {full} == paper's list; "
            f"all collinear along (7,2) so Succ_P(-1,4)=Succ_Q(-1,4)={d} == paper's (-2,7)")


@check("P43-11 GGV1 Prop 8.2(2) bound (k+1)b<a over the four corners gives k in {1,2}")
def _p11():
    corners = [(24, 7), (17, 5), (10, 3), (3, 1)]
    ks = set()
    per = {}
    for (a, b) in corners:
        kk = [k for k in range(1, 20) if (k + 1) * b < a]
        per[(a, b)] = kk
        ks |= set(kk)
    assert sorted(ks) == [1, 2], sorted(ks)
    return f"{per} -> k in {sorted(ks)} == paper's 'obtaining k in {{1,2}}'"


@check("P43-12 k=2 is impossible by parallelism; k=1 forces (en Q, en P)=((2,1),(-1,0))")
def _p12():
    """GGHV22 L1136: 'The case k=2 is impossible, as the edges of P and Q would
    have no way of being parallel.'  st_{rho,sigma}(P)=m(a,b), st(Q)=n(a,b), and
    the two edges must share the direction (rho,sigma)."""
    corners = [(24, 7), (17, 5), (10, 3), (3, 1)]
    m, n = 2, 3  # NB: the figure captions read (1/2)N(P) = (1/3)N(Q); see MN-1
    ok2, ok1 = [], []
    for (a, b) in corners:
        stP, stQ = (m * a, m * b), (n * a, n * b)
        for k in (1, 2):
            if (k + 1) * b >= a:
                continue
            for (eP, eQ) in (((-k, 0), (k + 1, 1)), ((k + 1, 1), (-k, 0))):
                wP = (eP[0] - stP[0], eP[1] - stP[1])
                wQ = (eQ[0] - stQ[0], eQ[1] - stQ[1])
                if cross(wP, wQ) == 0:
                    (ok1 if k == 1 else ok2).append(((a, b), eP, eQ, normal_pos_sum(wP)))
    assert ok2 == [], f"k=2 survived parallelism at {ok2}"
    assert ok1, "k=1 must survive"
    assigns = {(r[1], r[2]) for r in ok1}
    dirs = {r[3] for r in ok1}
    assert assigns == {((-1, 0), (2, 1))}, assigns
    assert dirs == {(-2, 7)}, dirs
    return (f"k=2: 0/{2*len(corners)*2 - 2} parallel-compatible assignments (impossible, as the paper says); "
            f"k=1: unique assignment en(P)=(-1,0), en(Q)=(2,1), common direction {dirs.pop()}")


@check("MN-1  the literal '(m,n)=(3,2)' at L1011 is inconsistent; (2,3) is forced")
def _mn():
    """GGHV22 L1009-1010 say the corners are '(0,0),(1,0),(8,28),(0,4) multiplied by
    (m,n)=(3,2)' (matching GGV5's table, arXiv:1708.07936 L1832), but the figure
    captions inside the same proof read (1/2)N(P) = (1/3)N(Q) and the assembled
    polygons at L1139-1140 read N(P)=2(28,8),..., N(Q)=3(28,8),....  Only the
    latter is compatible with Prop 8.2's output."""
    printed_NP = sorted([(-1, 0), (0, 0), (2 * 28, 2 * 8), (2 * 24, 2 * 7)])  # L1139
    out = {}
    for (m, n) in ((2, 3), (3, 2)):
        forced = set()
        for (a, b) in [(24, 7), (17, 5), (10, 3), (3, 1)]:
            stP, stQ = (m * a, m * b), (n * a, n * b)
            for (eP, eQ) in (((-1, 0), (2, 1)), ((2, 1), (-1, 0))):
                wP = (eP[0] - stP[0], eP[1] - stP[1])
                wQ = (eQ[0] - stQ[0], eQ[1] - stQ[1])
                if cross(wP, wQ) == 0:
                    forced.add((eP, eQ))
        assert len(forced) == 1, (m, n, forced)
        eP = next(iter(forced))[0]
        NP = sorted([eP, (0, 0), (m * 28, m * 8), (m * 24, m * 7)])
        out[(m, n)] = (eP, NP == printed_NP)
    assert out[(2, 3)][1] and not out[(3, 2)][1], out
    return (f"under (m_P,n_Q)=(2,3): en(P)={out[(2,3)][0]} and N(P) reproduces the paper's printed "
            f"L1139 list; under the literal (3,2): en(P)={out[(3,2)][0]} and N(P) would be the "
            f"3-multiple, i.e. P and Q swap roles relative to the Proposition statement. "
            f"L1009-1010's '(3,2)' (which matches GGV5 L1832) is therefore a P/Q labelling slip -- "
            f"mathematically harmless, the pair being symmetric, but it contradicts this proof's "
            f"own figure captions '(1/2)N(P)=(1/3)N(Q)' and its assembly at L1139-1140.")


@check("P43-13 assembled polygons in cases a)/b) and c) match GGHV22 L1139-1140 / L1184-1185")
def _p13():
    m, n = 2, 3
    NPab = hull([(-1, 0), (0, 0), (m * 28, m * 8), (m * 24, m * 7)])
    NQab = hull([(2, 1), (0, 0), (n * 28, n * 8), (n * 24, n * 7)])
    NPc = hull([(-1, 0), (0, 0), (m * 16, m * 4), (m * 28, m * 8), (m * 24, m * 7)])
    NQc = hull([(2, 1), (0, 0), (n * 16, n * 4), (n * 28, n * 8), (n * 24, n * 7)])
    assert len(NPab) == 4 and len(NQab) == 4, (NPab, NQab)
    assert len(NPc) == 5 and len(NQc) == 5, (NPc, NQc)
    assert set(NPc) == {(-1, 0), (0, 0), (32, 8), (56, 16), (48, 14)}, NPc
    return (f"a)/b): N(P) hull = {sorted(NPab)}, N(Q) hull = {sorted(NQab)} (4 vertices each); "
            f"c): 5 vertices each, every listed point is a genuine hull vertex")


@check("P43-14 final morphism (i,j)->(4j-i,j) lands on exactly the two stated configurations")
def _p14():
    with open(os.path.join(HERE, "paper_src", "upstream_facts.json"), encoding="utf-8") as fh:
        facts = json.load(fh)["facts"]["newton_polygons"]
    m, n = 2, 3
    f = lambda p: (4 * p[1] - p[0], p[1])
    NPab = [f(p) for p in [(-1, 0), (0, 0), (m * 28, m * 8), (m * 24, m * 7)]]
    NQab = [f(p) for p in [(2, 1), (0, 0), (n * 28, n * 8), (n * 24, n * 7)]]
    NPc = [f(p) for p in [(-1, 0), (0, 0), (m * 16, m * 4), (m * 28, m * 8), (m * 24, m * 7)]]
    NQc = [f(p) for p in [(2, 1), (0, 0), (n * 16, n * 4), (n * 28, n * 8), (n * 24, n * 7)]]
    sub2 = (sorted(map(tuple, facts["sub2"]["P"])), sorted(map(tuple, facts["sub2"]["Q"])))
    sub1 = (sorted(map(tuple, facts["sub1"]["P"])), sorted(map(tuple, facts["sub1"]["Q"])))
    assert (sorted(NPab), sorted(NQab)) == sub2, (sorted(NPab), sorted(NQab), sub2)
    assert (sorted(NPc), sorted(NQc)) == sub1, (sorted(NPc), sorted(NQc), sub1)
    return ("cases a),b) -> configuration (2) = sub2 = alok CASE_2; "
            "case c) -> configuration (1) = sub1 = alok CASE_1 "
            "(vertex-for-vertex against paper_src/upstream_facts.json)")


@check("P43-15 the chain-rule factor of x->x^-1, y->x^4 y is exactly -x^2")
def _p15():
    """d(x^-1)/dx * d(x^4 y)/dy - d(x^-1)/dy * d(x^4 y)/dx
       = (-x^-2)(x^4) - 0 = -x^2, so [phi(P),phi(Q)] = -[P,Q] x^2 (GGHV22 L1229)."""
    dXdx = {(-2, 0): Fraction(-1)}
    dXdy = {}
    dYdx = {(3, 1): Fraction(4)}
    dYdy = {(4, 0): Fraction(1)}
    det = padd(pmul(dXdx, dYdy), {k: -c for k, c in pmul(dXdy, dYdx).items()})
    assert det == {(2, 0): Fraction(-1)}, det
    return "det = -x^2  =>  [phi(P),phi(Q)] = -[P,Q]*x^2, and [P,Q] in K^x gives [P,Q]=x^2 after scaling"


# --- sensitivity of the ONE load-bearing unverified step -------------------


@check("SENS-1 sensitivity: what the second half would give if the (-1,4) edge were NOT a 7th power")
def _sens():
    """GGHV22 L1132 asserts the edge {(0,1),(28,8)} 'must be of the form
    y(x^4y-alpha)^7, corresponding to its form before the transformations' -- the
    ONE step of the proof for which no result is cited.  If instead the maximal
    multiplicity of a root were k < 7, phi_3 would reduce the edge only to
    {(4k-4,k),(28,8)}.  We re-run the calibrated Prop-8.2 machinery for each k."""
    m, n = 2, 3
    f = lambda p: (4 * p[1] - p[0], p[1])
    lines, top_viable = [], {}
    for k in range(2, 8):
        corner = (4 * k - 4, k)
        rows = prop82_table(corner, (-1, 4))
        surv = sorted((r[0] for r in rows if r[3]), key=lambda p: -p[1])
        full = [corner] + surv
        # (i) does the paper's own endpoint pattern close at the TOP corner?
        stP, stQ = (m * corner[0], m * corner[1]), (n * corner[0], n * corner[1])
        hit = []
        for kk in range(1, 20):
            if (kk + 1) * corner[1] >= corner[0]:
                continue
            for (eP, eQ) in (((-kk, 0), (kk + 1, 1)), ((kk + 1, 1), (-kk, 0))):
                if cross((eP[0] - stP[0], eP[1] - stP[1]),
                         (eQ[0] - stQ[0], eQ[1] - stQ[1])) == 0:
                    hit.append((kk, eP, eQ))
        # (ii) does anything at all survive, at any corner of the chain?
        anywhere = False
        for (a, b) in full:
            for kk in range(1, 20):
                if (kk + 1) * b >= a:
                    continue
                sP, sQ = (m * a, m * b), (n * a, n * b)
                for (eP, eQ) in (((-kk, 0), (kk + 1, 1)), ((kk + 1, 1), (-kk, 0))):
                    if cross((eP[0] - sP[0], eP[1] - sP[1]),
                             (eQ[0] - sQ[0], eQ[1] - sQ[1])) == 0:
                        anywhere = True
        if hit:
            kk, eP, eQ = hit[0]
            NP = sorted(f(p) for p in [eP, (0, 0), (m * 28, m * 8), stP])
            NQ = sorted(f(p) for p in [eQ, (0, 0), (n * 28, n * 8), stQ])
            top_viable[k] = (NP, NQ)
            lines.append(f"k={k}: closes at the top corner -> N(P)={NP}, N(Q)={NQ}")
        else:
            lines.append(f"k={k}: does NOT close at the top corner "
                         f"(survivors {surv}; survives further down the chain: {anywhere})")
    assert 7 in top_viable, "the paper's k=7 must close"
    assert top_viable[7][0] == [(0, 0), (1, 0), (8, 14), (8, 16)], top_viable[7]
    extra = sorted(set(top_viable) - {7})
    assert extra == [], extra
    return (" | ".join(lines) +
            "  ==>  k=7 is the ONLY multiplicity for which the GGV1-Prop-8.2 endpoint pattern "
            "closes directly at the reduced corner; for k<7 the argument would have to be "
            "continued down the chain of collinear corners, which the paper does not do. "
            "So the 7th-power claim is load-bearing but the failure mode is 'more work needed', "
            "not 'a configuration silently escapes at this step'.")


# ===========================================================================


def main(argv):
    quiet = "--quiet" in argv
    failures = []
    for name, fn in CHECKS:
        try:
            detail = fn()
            if not quiet:
                print(f"[PASS] {name}")
                if detail:
                    print(f"       {detail}")
        except AssertionError as exc:
            failures.append((name, str(exc)))
            print(f"[FAIL] {name}\n       {exc}", file=sys.stderr)
        except Exception as exc:  # noqa: BLE001
            failures.append((name, repr(exc)))
            print(f"[ERROR] {name}\n       {exc!r}", file=sys.stderr)
    if not quiet:
        print(f"\n{len(CHECKS)-len(failures)}/{len(CHECKS)} checks passed.")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
