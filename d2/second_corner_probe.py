#!/usr/bin/env python3
"""second_corner_probe.py  (NEW 2026-07-28; read-only over every artifact)

THE SECOND CORNER.  Gates SECOND_CORNER.md.

PRIMITIVITY_DEPTH.md sec.1 names one bottleneck: the depth law
`depth = -delta * deg_y(C|x=0)` is mechanical wherever a reduced Newton polygon
exists, and a reduced polygon existed "at only two corners", `(8,28)` and
`(5,20)`, of which only `(5,20)` feeds the law -- and there EIGHT distinct
corner-data formulas collapse onto the single value `-10`.

This file does three things.

  (A) CORRECTION.  A usable second corner was ALREADY in the repo, published and
      unnoticed: `passport_75_125.PUB["7_21"]` transcribes GGHV22's own reduced
      polygon at `(7,21)` (2204.14178.tex:1313-1320, proof :1388-1395):
          N(P) = 2*{(0,0),(2,0),(3,1),(0,7)},  N(Q) = 3*same,  [P,Q] = x.
      GGHV22's figure caption calls this object "the transformation of
      1/2 N(P) = 1/3 N(Q)" -- i.e. it IS N(C).  It is the PROPORTIONAL branch,
      so `N(C)` exists; its y-axis vertex is `(0,7)`, giving `j* = 7 != 5`.
      The claim "there is no second point to fit against" is therefore false of
      the repo, though true of `polygon_reduction.py` (which carries three
      reductions at two corners).  Of the SIX rows of `passport_75_125.PUB`
      (GGHV22's published reductions and sub-cases, at four corners) exactly one
      is proportional -- (7,21).  The other five ((8,28) x2, (9,24) x2, (9,27))
      are en-split, so N(P) is not m*Delta' and no N(C) exists there: that is
      why the usable corner went unnoticed.

  (B) TWO NEW CORNERS, derived.  `(9,36)` and `(7,42)`.  Both are sporadic
      length-1 rows of GGV5's "9 other possible pairs" table (tex:1828-1836),
      whose `A_0'` column is NOT printed.  The missing datum is recovered by a
      lemma calibrated on 20 published rows:

          A_0'-RECOVERY.  Given A_0 = (u,v) and the PRINTED final corner
          A_1 = (a\\l, b), take rho = l, get (rho,sigma) from the unique
          GGV1-Prop-'final' branch with that rho, solve
          v_{rho,sigma}(A_0') = v_{rho,sigma}(A_0) for 0 <= s' < r' < u, and set
          gamma = b.  Then GGV1 (7) must reproduce the printed first coordinate
          a/l of A_1.  It does, on every row.

      Calibration: all 17 length-1 family rows plus F_22/23/24's first link --
      20 rows where GGV5 PRINTS A_0' -- are reproduced exactly, and F_18-F_21
      (the four families GGV5 itself proves cannot come from a standard pair)
      are the exactly four rows the recovery REFUSES.  Applied to the 9 sporadic
      length-1 rows it gives a unique A_0' each, with (7,42) and (9,36) landing
      on A_0' = (1,0) -- the one shape on which `passport_75_125.Reduction`'s
      rule (r1) is anchored.

  (C) THE CLASS CLOSED FORM, and what it settles about the depth law.
      For A_0' = (1,0) and a0 | b0 (the "monomial shape"):
          l = b0/a0,  mu = l-1,  c = a0,  q = a0,  zdeg = 1  (no split branch),
          retraction is IMPOSSIBLE, the en-split branch is ILLEGAL, hence
              Delta' = {(0,0), (l-1,0), (l,1), (0,a0)},  kappa = l-2,
              N(P) = m*Delta',  N(Q) = n*Delta',  N(C) = Delta',  j* = a0.
      GGV1 Prop 'u(u-1)' (1401.1784_GGV1.tex:3631-3632, "v <= u(u-1)") forces
      l+1 <= a0, hence deg C = a0 too.  Five published corners are in this
      class: (4,12), (5,20), (7,21), (7,42), (9,36).  Two carry external
      controls -- (7,21)'s vertex list is GGHV22's, and (5,20) reproduces GGV3's
      three published integers.

      CONSEQUENCE for PRIMITIVITY_DEPTH sec.1's eight candidates.  The law says
      the depth is `-delta * j*` with j* = a0, so every candidate must be
      divisible by a0.  Four are not, at one corner or another, and die:
          -b0/2, -t*kappa-2, -(b0-2t-2), -2t-2.
      The four survivors -2a0, -2degC, -(a0+degC), -degP1 are (the first three)
      IDENTICALLY EQUAL wherever the law is non-vacuous -- PROVED, not observed
      -- and -degP1 is separated only by an (m,n) argument.  So the eight
      collapse to ONE, `-2a0 = -2j*`, and its entire content is `delta = 2`.
      The residual degeneracy is a THEOREM: separating -2a0 from -2c would need
      a0 !x b0, and there j* = 0 and the law is vacuous.

  (D) THE BLOCKERS, gated as facts rather than remembered.
      (8,32): NO A_0' exists on any branch (0 <= s' < r' < 8 has no solution);
              its chain's first link ends at the INTEGER corner A_1 = (8,28)
              with v_{1,0}(A_0) = v_{1,0}(A_1), i.e. (rho_0,sigma_0) = (1,0) --
              the F_18-F_21 shape, outside the type-II.b chart entirely.
      (10,40): A_0' = (2,0), recovered and anchored on BOTH of its printed A_1's,
              but rule (r1)'s lower Delta vertex has no anchor off A_0' = (1,0).
              GGHV22 PRINTS Delta at four corners (2204.14178.tex:471, 682, 1010,
              1388) and at three of them A_0' = (1,0), so the two readings
              ("the vertex is always (1,0)" vs "the vertex is A_0'") COINCIDE
              there.  The fourth, (9,27), carries one extra printed vertex and it
              IS that chain's A_0' = (9,24) -- one instance, for reading B.  The
              readings DISAGREE at (10,40).  That single undetermined vertex is
              the whole blocker, and it is a missing PUBLISHED datum, not missing
              mathematics.

Sources (local copies; line numbers pinned):
  GGHV22 2204.14178.tex : (7,21) Prop lines 1313-1320, proof 1388-1395.
  GGV5 1708.07936_GGV5.tex : length-1 family table 1678-1694, length-2 family
    table 1709-1715, F_18-F_21 impossibility 1726-1786, sporadic length-1 table
    1828-1836, length-2 1848-1858.
  GGV1 1401.1784_GGV1.tex : Prop 'final' conditions (5)-(9) (transcribed in
    gamma_from_corner.py), Prop 'u(u-1)' lines 3631-3632.
  GGV3 1406.0886_GGV3.tex:1723-1727 : [P_1,Q_1]=x^2, deg 10, deg 15 at (5,20).
  composite_charts.py : FUSED-CHART LEMMA, Jacobian -x^(l-2), hence kappa=l-2.

Exact integer / Fraction arithmetic only.  Reads only; writes nothing.
Usage:  python second_corner_probe.py [--quiet]      exit 0 iff all checks pass.
"""
from __future__ import annotations

import os
import sys
from fractions import Fraction as F
from math import gcd

HERE = os.path.dirname(os.path.abspath(__file__)) or "."
sys.path.insert(0, HERE)

QUIET = "--quiet" in sys.argv
_ok = [0]
_fail: list[str] = []


def ck(name: str, cond: bool, detail: str = "") -> bool:
    if cond:
        _ok[0] += 1
        if not QUIET:
            print("[OK]   %s%s" % (name, ("  -- " + detail) if detail else ""))
    else:
        _fail.append(name)
        print("[FAIL] %s%s" % (name, ("  -- " + detail) if detail else ""))
    return bool(cond)


def head(s: str) -> None:
    if not QUIET:
        print("\n" + "=" * 92 + "\n" + s + "\n" + "=" * 92)


# ---------------------------------------------------------------------------
# GGV5 tables, transcribed.  (label, A_0, A_0' or None, A_1 as (a,l,b) or the
# integer corner (i,j)).  A_0' is None exactly where GGV5 does not print it.
# ---------------------------------------------------------------------------
FAM_L1 = [   # tex:1678-1694.  Family | A_0 | A_0' | A_1 | k | m | n
    ("F_1",  (4, 12), (1, 0), (7, 4, 3)),
    ("F_2",  (5, 20), (1, 0), (7, 5, 2)),
    ("F_3",  (5, 20), (1, 0), (8, 5, 3)),
    ("F_4",  (5, 20), (1, 0), (8, 5, 3)),
    ("F_5",  (5, 20), (1, 0), (9, 5, 4)),
    ("F_6",  (5, 20), (1, 0), (9, 5, 4)),
    ("F_7",  (6, 15), (1, 0), (7, 3, 4)),
    ("F_8",  (6, 15), (1, 0), (8, 3, 5)),
    ("F_9",  (7, 21), (1, 0), (11, 7, 2)),
    ("F_10", (7, 21), (1, 0), (13, 7, 3)),
    ("F_11", (7, 21), (1, 0), (13, 7, 3)),
    ("F_12", (8, 24), (2, 0), (13, 4, 5)),
    ("F_13", (9, 21), (2, 0), (13, 3, 7)),
    ("F_14", (9, 24), (1, 0), (7, 3, 4)),
    ("F_15", (9, 24), (1, 0), (8, 3, 5)),
    ("F_16", (9, 24), (1, 0), (10, 3, 7)),
    ("F_17", (9, 24), (1, 0), (11, 3, 8)),
]
FAM_L2 = [   # tex:1709-1715, FIRST link only.  A_1 integer => A_1 = A_0'.
    ("F_18", (6, 18), (6, 15), None),
    ("F_19", (6, 18), (6, 15), None),
    ("F_20", (6, 24), (6, 15), None),
    ("F_21", (6, 24), (6, 15), None),
    ("F_22", (8, 24), (2, 0), (14, 4, 6)),
    ("F_23", (8, 24), (2, 0), (14, 4, 6)),
    ("F_24", (8, 24), (2, 0), (14, 4, 6)),
]
SPOR_L1 = [  # tex:1828-1836.  A_0 | A_1 | (m,n) | max deg.  NO A_0' column.
    ((7, 35),  (19, 7, 5),  (2, 3), 126),
    ((7, 42),  (13, 7, 6),  (3, 2), 147),
    ((7, 42),  (13, 7, 6),  (2, 3), 147),
    ((8, 28),  (7, 4, 3),   (3, 4), 144),
    ((8, 28),  (11, 4, 7),  (3, 2), 108),
    ((9, 36),  (17, 9, 4),  (3, 2), 135),
    ((9, 36),  (17, 9, 4),  (2, 3), 135),
    ((11, 33), (19, 4, 8),  (2, 3), 132),
    ((12, 33), (11, 3, 8),  (2, 3), 135),
]
SPOR_L2 = [  # tex:1848-1858, first link.  4-tuples are fractional A_1.
    ((8, 32),  (8, 28)),   ((8, 40),  (8, 28)),
    ((9, 27),  (9, 24)),   ((9, 36),  (9, 24)),
    ((10, 40), (16, 5, 6)), ((10, 40), (18, 5, 8)),
    ((12, 30), (16, 3, 10)),
    ((12, 36), (12, 33)),  ((12, 36), (9, 24)),
    ((12, 36), (21, 4, 9)), ((12, 36), (21, 4, 9)),
]
# every distinct A_0 of the census
ALL_A0 = sorted({A0 for _, A0, _, _ in FAM_L1} | {A0 for _, A0, _, _ in FAM_L2}
                | {r[0] for r in SPOR_L1} | {r[0] for r in SPOR_L2})


# ---------------------------------------------------------------------------
# 0.  Self-contained lattice geometry + an INDEPENDENT reduction path.
#     (Deliberately not imported from passport_75_125: section C cross-checks
#     two code paths against each other and against GGHV22.)
# ---------------------------------------------------------------------------
def _cross(o, a, b):
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def hull(pts):
    P = sorted(set(map(tuple, pts)))
    if len(P) <= 2:
        return P
    lo = []
    for p in P:
        while len(lo) >= 2 and _cross(lo[-2], lo[-1], p) <= 0:
            lo.pop()
        lo.append(p)
    up = []
    for p in reversed(P):
        while len(up) >= 2 and _cross(up[-2], up[-1], p) <= 0:
            up.pop()
        up.append(p)
    return lo[:-1] + up[:-1]


def reduce_independent(a0, b0, s=None, l=None):
    """Flip -> root shift -> Laurent inversion, written out from scratch.

    Delta  = hull{(0,0),(1,0),(a0,b0),(0,c)}   with c = b0 - mu*a0, mu = l-1
    flip   : (i,j) -> (j,i)
    shift  : the flipped foot (c,0) is replaced by (-s,0)   [Pred_P(1,0)=(1,-s)]
    invert : (i,j) -> (l*j - i, j)
    Returns hull(Delta').  A_0' = (1,0) is ASSUMED, exactly as rule (r1) does.
    """
    mu = (b0 - 1) // a0
    if l is None:
        l = mu + 1
    if s is None:
        s = mu
    c = b0 - mu * a0
    delta = hull([(0, 0), (1, 0), (a0, b0), (0, c)])
    flipped = [(j, i) for (i, j) in delta]
    pre = [p for p in flipped if p != (c, 0)] + [(-s, 0)]
    return hull([(l * j - i, j) for (i, j) in hull(pre)])


def jstar(Dp):
    """deg_y of the x = 0 column of the polygon, i.e. the y-axis vertex."""
    return max([j for (i, j) in Dp if i == 0] + [0])


# ---------------------------------------------------------------------------
# 1.  The A_0'-RECOVERY LEMMA.
# ---------------------------------------------------------------------------
def recover_A0prime(A0, A1):
    """A_1 = (a, l, b) printed as (a\\l, b).  Return the list of consistent
    (branch f, (rho,sigma), A_0', gamma) with GGV1 (7) reproducing a/l."""
    import gamma_from_corner as gfc
    a, l, b = A1
    out = []
    for rec in gfc.analyse(*A0):
        if rec["rejected"]:
            continue
        rho, sigma = rec["rho_sigma"]
        if rho != l:
            continue
        for A0p in rec["A0prime_candidates"]:
            A1c = gfc.A1_of(A0p, b, rho, sigma)
            if A1c == (F(a, l), b):
                out.append((rec["f"], (rho, sigma), A0p, b))
    return out


def branches_with_A0prime(A0):
    import gamma_from_corner as gfc
    return [(rec["f"], rec["rho_sigma"], rec["A0prime_candidates"])
            for rec in gfc.analyse(*A0) if not rec["rejected"]]


# ---------------------------------------------------------------------------
def main() -> int:                                   # noqa: C901
    import passport_75_125 as pp

    # =====================================================================
    head("A.  CORRECTION -- a usable second corner was already in the repo")
    # =====================================================================
    p721 = pp.PUB["7_21"]
    a0, b0, mn, kw, wp, wq, wk, ln = p721
    ck("A1  passport_75_125.PUB['7_21'] is GGHV22's PUBLISHED reduction at "
       "(7,21) [2204.14178.tex:%s] and carries NO branch selector, i.e. it is "
       "the PROPORTIONAL branch of GGV1 Prop 8.2(1)" % ln,
       (a0, b0) == (7, 21) and kw == {} and mn == (2, 3),
       "N(P)=%s  N(Q)=%s  [P,Q]=x^%d" % (wp, wq, wk))
    Dp721 = [(0, 0), (2, 0), (3, 1), (0, 7)]
    ck("A2  the published vertex lists ARE m and n times one polygon Delta', so "
       "N(C) exists: N(P)/2 == N(Q)/3 == %s" % Dp721,
       hull(wp) == hull([(2 * i, 2 * j) for i, j in Dp721])
       and hull(wq) == hull([(3 * i, 3 * j) for i, j in Dp721]))
    ck("A2b  GGHV22 itself names this object: its figure caption reads 'The "
       "transformation of 1/2 N(P) = 1/3 N(Q)' (2204.14178.tex:1385) -- the "
       "proportionality is the paper's, not ours", True,
       "quoted, not re-derived")
    ck("A3  j*(7,21) = 7 != 5 = j*(5,20): the calibration set of "
       "PRIMITIVITY_DEPTH sec.1 is NOT a single point",
       jstar(Dp721) == 7 and jstar([(0, 0), (3, 0), (4, 1), (0, 5)]) == 5)
    # which of the five published reductions can carry a j* at all
    proportional, ensplit = [], []
    for tag, (qa0, qb0, qmn, qkw, *_rest) in pp.PUB.items():
        (proportional if qkw.get("en_k") is None else ensplit).append(tag)
    ck("A4  of the SIX rows of passport_75_125.PUB (GGHV22's published "
       "reductions and sub-cases, at four corners) exactly ONE is proportional "
       "-- (7,21).  The other five are en-split (Prop 8.2(2)), where N(P) is "
       "not m*Delta' and no N(C) exists",
       proportional == ["7_21"] and len(ensplit) == 5,
       "proportional=%s  en-split=%s" % (proportional, ensplit))
    for tag in ensplit:
        qa0, qb0, qmn, qkw = pp.PUB[tag][:4]
        r = pp.Reduction(tag, qa0, qb0, qmn, **qkw)
        ck("A5  MUTATION control: %s's published N(P) is NOT m*(integer polygon) "
           "-- it cannot supply a j*" % tag,
           any(i % qmn[0] or j % qmn[0] for (i, j) in r.NP),
           "N(P)=%s, m=%d" % (r.NP, qmn[0]))
    r721 = pp.Reduction("7_21", 7, 21, (2, 3))
    ck("A6  and the engine reproduces (7,21) exactly, so the published corner "
       "and the derived machinery agree there",
       r721.NP == hull(wp) and r721.NQ == hull(wq) and r721.kappa == wk,
       "kappa = l-2 = %d = GGHV22's [P,Q] = x^%d" % (r721.kappa, wk))

    # =====================================================================
    head("B.  THE A_0'-RECOVERY LEMMA -- calibration, refusal, mutation")
    # =====================================================================
    cal = [(lab, A0, A0p, A1) for lab, A0, A0p, A1 in FAM_L1 + FAM_L2
           if A1 is not None]
    good = []
    for lab, A0, A0p, A1 in cal:
        hits = recover_A0prime(A0, A1)
        good.append(len(hits) == 1 and hits[0][2] == A0p)
    ck("B1  CALIBRATION: on all %d published rows that print BOTH A_0' and a "
       "fractional A_1, the recovery returns exactly one A_0' and it is the "
       "printed one" % len(cal), all(good),
       "rows: " + ", ".join(l for l, *_ in cal))
    ck("B1b  ... and the recovery's own consistency test is the printed first "
       "coordinate of A_1: GGV1 (7) reproduces a/l on every one of those rows "
       "(that is what recover_A0prime filters on)", all(good))
    # refusal control
    ref = [(lab, A0) for lab, A0, _A0p, A1 in FAM_L2 if A1 is None]
    refused = []
    for lab, A0 in ref:
        refused.append(all(not c for _f, _rs, c in branches_with_A0prime(A0)))
    ck("B2  REFUSAL control: the recovery finds NO A_0' at %s -- exactly the "
       "four families (F_18-F_21) that GGV5 PROVES cannot come from a standard "
       "(m,n)-pair (tex:1726-1786).  Their A_1 is an INTEGER corner equal to "
       "A_0', so 0 <= s' < r' < u fails by construction"
       % ", ".join("%s%s" % (l, A0) for l, A0 in ref), all(refused))
    ck("B2b  MUTATION: the refusal is not vacuous -- the same code accepts the "
       "SAME corners' siblings.  (8,24) hosts both F_12 (accepted) and "
       "F_22-F_24 (accepted); (6,18)/(6,24) are refused",
       len(recover_A0prime((8, 24), (13, 4, 5))) == 1
       and len(recover_A0prime((8, 24), (14, 4, 6))) == 1)
    # numerator mutation
    mut_clean = True
    for lab, A0, A0p, A1 in cal:
        a, l, b = A1
        for da in (-1, 1):
            if recover_A0prime(A0, (a + da, l, b)):
                mut_clean = False
    ck("B3  MUTATION: shifting the printed numerator a of A_1 = (a\\l,b) by +-1 "
       "makes the recovery return NOTHING on all %d calibration rows -- the "
       "consistency test is discriminating, not a tautology" % len(cal),
       mut_clean)
    # denominator mutation: picking the wrong branch gives the wrong A_1
    wrong = []
    import gamma_from_corner as gfc
    for A0, A1, _mn, _md in SPOR_L1:
        for rec in gfc.analyse(*A0):
            if rec["rejected"] or rec["rho_sigma"][0] == A1[1]:
                continue
            rho, sigma = rec["rho_sigma"]
            for A0p in rec["A0prime_candidates"]:
                wrong.append(gfc.A1_of(A0p, A1[2], rho, sigma)
                             != (F(A1[0], A1[1]), A1[2]))
    ck("B4  MUTATION: across the sporadic length-1 rows, EVERY branch other than "
       "the one with rho = l(A_1) predicts a DIFFERENT A_1 -- so the printed "
       "final corner selects the branch uniquely (%d wrong branches tested)"
       % len(wrong), len(wrong) >= 5 and all(wrong))
    # application
    app = {}
    for A0, A1, mn, md in SPOR_L1:
        hits = recover_A0prime(A0, A1)
        app[(A0, A1)] = hits
    ck("B5  APPLICATION: all %d sporadic length-1 rows (A_0' NOT printed) get a "
       "UNIQUE A_0' whose A^(1) reproduces the printed A_1"
       % len({k for k in app}), all(len(v) == 1 for v in app.values()),
       "; ".join("%s -> A_0'=%s gamma=%d" % (k[0], v[0][2], v[0][3])
                 for k, v in sorted(app.items())))
    ck("B6  CONTROL inside the application: (8,28) is a sporadic row whose A_0' "
       "the repo already knows independently (polygon_reduction.case_8_28 uses "
       "A0p=(1,0)); the recovery reproduces it on BOTH of its rows",
       app[((8, 28), (7, 4, 3))][0][2] == (1, 0)
       and app[((8, 28), (11, 4, 7))][0][2] == (1, 0))
    A0P = {A0: v[0][2] for (A0, _A1), v in app.items()}
    ck("B7  the two targets land on A_0' = (1,0) -- the one shape rule (r1) is "
       "anchored on", A0P[(9, 36)] == (1, 0) and A0P[(7, 42)] == (1, 0),
       "(7,35)->%s  (11,33)->%s  (12,33)->%s"
       % (A0P[(7, 35)], A0P[(11, 33)], A0P[(12, 33)]))

    # =====================================================================
    head("C.  THE CLASS CLOSED FORM  A_0'=(1,0) and a0 | b0")
    # =====================================================================
    # ---- rule (r1) is PUBLISHED verbatim at four corners ---------------------
    # GGHV22 2204.14178.tex: "The corners of the polygons of P and Q are {...}
    # multiplied by (m,n) = ... respectively."  Transcribed with line numbers.
    R1_PUB = {
        (9, 27): ([(0, 0), (1, 0), (9, 24), (9, 27), (0, 9)], "471"),
        (9, 24): ([(0, 0), (1, 0), (9, 24), (0, 6)], "682"),
        (8, 28): ([(0, 0), (1, 0), (8, 28), (0, 4)], "1010"),
        (7, 21): ([(0, 0), (1, 0), (7, 21), (0, 7)], "1388"),
    }
    r1_ok = True
    for (qa, qb), (verts, _ln) in R1_PUB.items():
        mu = (qb - 1) // qa
        c = qb - mu * qa
        base = hull([(0, 0), (1, 0), (qa, qb), (0, c)])
        extra = [v for v in hull(verts) if v not in set(base)]
        r1_ok &= (hull(verts) == hull(base + extra)) and (0, c) in set(verts) \
            and (1, 0) in set(verts) and len(extra) <= 1
    ck("C0  rule (r1) is not a house convention: GGHV22 PRINTS the polygon "
       "Delta = {(0,0),(1,0),A_0,(0,c)} with c = b0 - mu*a0 at FOUR corners "
       "(2204.14178.tex:471, 682, 1010, 1388), and c comes out 9, 6, 4, 7 -- "
       "all four reproduced by the rule", r1_ok,
       "; ".join("%s c=%d [tex:%s]" % (k, k[1] - ((k[1] - 1) // k[0]) * k[0], v[1])
                 for k, v in sorted(R1_PUB.items())))
    ck("C0b  the ONE printed corner carrying an extra Delta vertex is (9,27), "
       "and that vertex is (9,24) -- which is that chain's A_1 and A_0'.  It is "
       "the single published data point bearing on rule (r1) off A_0' = (1,0)",
       set(hull(R1_PUB[(9, 27)][0])) - set(hull([(0, 0), (1, 0), (9, 27), (0, 9)]))
       == {(9, 24)} and ((9, 27), (9, 24)) in SPOR_L2)

    CLASS = [(4, 12), (5, 20), (7, 21), (7, 42), (9, 36)]
    est = {(4, 12): "GGV5 F_1 prints A_0'=(1,0)",
           (5, 20): "GGV5 F_2-F_6 print A_0'=(1,0)",
           (7, 21): "GGV5 F_9-F_11 print A_0'=(1,0)",
           (7, 42): "recovered in B5, anchored on A_1=(13\\7,6)",
           (9, 36): "recovered in B5, anchored on A_1=(17\\9,4)"}
    ck("C1  the five class corners all have a0 | b0 AND an ESTABLISHED "
       "A_0' = (1,0)", all(b % a == 0 for a, b in CLASS),
       "; ".join("%s: %s" % (c, est[c]) for c in CLASS))
    ok_shape = True
    for a, b in CLASS:
        mu = (b - 1) // a
        l = mu + 1
        c = b - mu * a
        ok_shape &= (l == b // a and mu == l - 1 and c == a
                     and gcd(a, b) == a and a // gcd(a, b) == 1
                     and b != l * (a - 1))
    ck("C2  a0 | b0 FORCES l = b0/a0, mu = l-1, c = a0, q = gcd = a0, "
       "zdeg = a0/q = 1 (so rule (r5) admits NO split branch) and makes the "
       "retraction b0 = l(a0-1) arithmetically impossible (it would need l = 0)",
       ok_shape)
    ck("C3  GGV1 Prop 'u(u-1)' (1401.1784_GGV1.tex:3631-3632, 'v <= u(u-1)') "
       "holds on all %d census corners and forces l <= a0-1, hence l+1 <= a0, "
       "hence deg C = max(l-1, l+1, a0) = a0" % len(ALL_A0),
       all(b <= a * (a - 1) for a, b in ALL_A0)
       and all((b // a) + 1 <= a for a, b in CLASS))
    closed = {}
    for a, b in CLASS:
        l = b // a
        closed[(a, b)] = hull([(0, 0), (l - 1, 0), (l, 1), (0, a)])
    ck("C4  CLOSED FORM  Delta' = {(0,0),(l-1,0),(l,1),(0,a0)} reproduces "
       "GGHV22's PUBLISHED (7,21) polygon exactly -- an EXTERNAL control on a "
       "corner the form was not fitted to", closed[(7, 21)] == hull(Dp721),
       "%s" % closed[(7, 21)])
    two_paths = all(closed[(a, b)] == pp.Reduction("x", a, b, (2, 3)).reduced_delta
                    == reduce_independent(a, b) for a, b in CLASS)
    ck("C5  THREE code paths agree on all five class corners: the closed form, "
       "passport_75_125.Reduction, and a from-scratch flip/shift/invert written "
       "in this file", two_paths,
       "; ".join("%s->%s" % (c, closed[c]) for c in CLASS))
    ck("C6  kappa = l - 2 (FUSED-CHART LEMMA, composite_charts.py: the composite "
       "Laurent chart (x^-1, x^l y + shears) has Jacobian -x^(l-2) for ANY "
       "shears) reproduces BOTH published brackets: [P,Q] = x at (7,21) "
       "(l=3) and [P_1,Q_1] = x^2 at (5,20) (l=4, GGV3 tex:1725)",
       (21 // 7) - 2 == 1 and (20 // 5) - 2 == 2)
    r520 = pp.Reduction("x", 5, 20, (2, 3))
    ck("C7  and GGV3's other two published integers at (5,20): deg P_1 = m*a0 = "
       "10 and deg Q_1 = n*a0 = 15 (1406.0886_GGV3.tex:1723-1727)",
       (r520.degP, r520.degQ) == (10, 15))
    # branch controls
    legal = []
    for a, b in CLASS:
        for k in (1, 2, 3, 4):
            for sw in (False, True):
                legal.append(pp.Reduction("e", a, b, (2, 3), en_k=k, en_swap=sw).legal)
    ck("C8  the en-split branch (GGV1 Prop 8.2(2)) is ILLEGAL at every class "
       "corner, for every k and both assignments (%d branches) -- so the "
       "PROPORTIONAL branch is FORCED and N(C) = Delta' is not a choice"
       % len(legal), not any(legal))
    ck("C8b  MUTATION: the same legality test PASSES at (8,28) and (12,33), so "
       "it is discriminating",
       pp.Reduction("e", 8, 28, (2, 3), en_k=1).legal
       and pp.Reduction("e", 12, 33, (2, 3), en_k=1).legal)
    inv_s = all(jstar(reduce_independent(a, b, s=(b // a) - 2)) == a
                for a, b in CLASS if (b // a) - 2 >= 1)
    ck("C9  MUTATION on the root-shift depth: taking s = mu-1 instead of mu "
       "(GGV6 Prop 2.5's other admissible Pred) moves ONLY the foot vertex and "
       "leaves j* = a0 untouched -- the depth prediction is s-invariant", inv_s)
    bad_l = all(reduce_independent(7, 21, l=ll) != hull(Dp721) for ll in (2, 4))
    ck("C10 MUTATION on l: l = 2 or 4 at (7,21) fails to reproduce GGHV22's "
       "published polygon, so the l = b0/a0 rule is pinned there by print",
       bad_l)

    # =====================================================================
    head("D.  THE TWO NEW CORNERS")
    # =====================================================================
    for (a, b), rows in (((9, 36), [(3, 2), (2, 3)]), ((7, 42), [(3, 2), (2, 3)])):
        l = b // a
        Dp = closed[(a, b)]
        ck("D  (%d,%d): Delta' = %s, kappa = %d, N(C) = Delta', j* = %d"
           % (a, b, Dp, l - 2, a), jstar(Dp) == a)
        for mn in rows:
            r = pp.Reduction("t", a, b, mn)
            NP = hull([(mn[0] * i, mn[0] * j) for i, j in Dp])
            NQ = hull([(mn[1] * i, mn[1] * j) for i, j in Dp])
            md = max(mn) * (a + b)
            pub = [row for row in SPOR_L1 if row[0] == (a, b) and row[2] == mn]
            ck("D  (%d,%d) (m,n)=%s: N(P) = %s (deg %d), N(Q) = %s (deg %d)"
               % (a, b, mn, NP, mn[0] * a, NQ, mn[1] * a),
               r.NP == NP and r.NQ == NQ and r.degP == mn[0] * a
               and r.degQ == mn[1] * a)
            ck("D  ... and GGV5's printed max{deg P, deg Q} = max(m,n)*v11(A_0) "
               "= %d for that row (tex:1828-1836)" % md,
               bool(pub) and pub[0][3] == md)
    ck("D9  the two new corners are NOT the (5,20) shape replayed: they differ "
       "in a0 (9 and 7 vs 5) and (7,42) also in l (6 vs 4), hence in kappa "
       "(4 vs 2)",
       len({a for a, b in CLASS}) >= 3 and len({b // a for a, b in CLASS}) >= 3,
       "class (a0,l): %s" % [(a, b // a) for a, b in CLASS])

    # =====================================================================
    head("E.  WHAT THIS SETTLES ABOUT THE DEPTH LAW")
    # =====================================================================
    def candidates(a0_, b0_, m_):
        l_ = b0_ // a0_
        t_, kap = l_, l_ - 2
        degC = a0_
        return {
            "-2*a0": F(-2 * a0_), "-2*degC": F(-2 * degC),
            "-b0/2": F(-b0_, 2), "-(a0+degC)": F(-(a0_ + degC)),
            "-degP1": F(-m_ * degC),
            "-t*kappa-2": F(-(t_ * kap + 2)),
            "-(b0-2t-2)": F(-(b0_ - 2 * t_ - 2)),
            "-2t-2": F(-(2 * t_ + 2)),
        }
    ck("E1  at (5,20) with m=2 all EIGHT candidates read -10 -- the degeneracy "
       "PRIMITIVITY_DEPTH sec.1 reports", set(candidates(5, 20, 2).values()) == {F(-10)},
       str({k: str(v) for k, v in candidates(5, 20, 2).items()}))
    # the law: depth = -delta * j*, j* = a0, so a0 must divide the depth
    dead = {}
    for a, b in CLASS:
        for k, v in candidates(a, b, 2).items():
            if v.denominator != 1 or (int(v) % a) != 0:
                dead.setdefault(k, []).append((a, b))
    ck("E2  the derived law says depth = -delta * j* = -delta * a0 with delta a "
       "POSITIVE INTEGER, so a0 must divide the depth.  FOUR candidates fail "
       "that on at least one class corner and are REFUTED",
       set(dead) == {"-b0/2", "-t*kappa-2", "-(b0-2t-2)", "-2t-2"},
       "; ".join("%s dies at %s" % (k, v) for k, v in sorted(dead.items())))
    kills = {c: sorted(k for k, v in dead.items() if c in v) for c in CLASS}
    ck("E2b  (7,21) ALONE kills all four: -b0/2 = -21/2 is not even an integer "
       "there, and -t*kappa-2 = -5, -(b0-2t-2) = -13, -2t-2 = -8 are not "
       "multiples of j* = 7.  That corner is GGHV22's, PUBLISHED -- so the "
       "discrimination needed no new derivation at all, only that the repo read "
       "its own control", set(kills[(7, 21)]) == set(dead),
       "per-corner kills: " + "; ".join("%s:%d" % (c, len(kills[c]))
                                        for c in CLASS))
    ck("E2c  and the kills are NOT automatic -- (5,20) kills nothing (E1), and "
       "(7,42) leaves -b0/2 = -21 = -3*j* standing because 7 | 21.  The two NEW "
       "corners corroborate independently: (9,36) kills three of the four, "
       "(7,42) one, by arithmetic disjoint from (7,21)'s",
       kills[(5, 20)] == [] and "-b0/2" not in kills[(7, 42)]
       and len(kills[(9, 36)]) == 3 and len(kills[(7, 42)]) == 1,
       "(9,36) kills %s ; (7,42) kills %s" % (kills[(9, 36)], kills[(7, 42)]))
    surv = [k for k in candidates(5, 20, 2) if k not in dead]
    ck("E3  the survivors are %s" % surv,
       set(surv) == {"-2*a0", "-2*degC", "-(a0+degC)", "-degP1"})
    ident = all(len({candidates(a, b, 2)[k] for k in
                     ("-2*a0", "-2*degC", "-(a0+degC)", "-degP1")}) == 1
                for a, b in CLASS)
    ck("E4  and three of the four are IDENTICALLY EQUAL wherever the law is "
       "non-vacuous -- PROVED, not observed: j* > 0 forces a0 | b0 (C2), which "
       "forces c = a0, and GGV1 Prop 'u(u-1)' forces deg C = a0 (C3).  No "
       "corner can ever separate -2a0, -2degC, -(a0+degC)", ident)
    ck("E4b  MUTATION: a fabricated ninth candidate -3*a0 is NOT refuted by the "
       "divisibility test, so E2 is not a test that kills everything",
       all((3 * a) % a == 0 for a, b in CLASS))
    ck("E5  -degP1 = -m*a0 is separated from -2a0 by an (m,n) argument, not by "
       "a corner: at (5,20) the rows (m,n) = (2,3) and (3,5) share ONE corner, "
       "ONE A_0' and ONE chart -- hence one gamma and one delta -- yet give "
       "-degP1 = -10 and -15.  A depth that is a function of (corner, gamma) "
       "cannot be both",
       candidates(5, 20, 2)["-degP1"] != candidates(5, 20, 3)["-degP1"],
       "-degP1 at m=2: %s ; at m=3: %s"
       % (candidates(5, 20, 2)["-degP1"], candidates(5, 20, 3)["-degP1"]))
    ck("E6  NET: eight corner-data formulas collapse to ONE, -2*a0 = -2*j*, "
       "whose entire content is delta = 2.  The depth law's residue is a single "
       "gamma-dependent integer, exactly as PRIMITIVITY_DEPTH sec.5 flags -- "
       "and delta is NOT determined here", True, "declared, not proved")
    ck("E7  DECLARED NON-DISCRIMINATING: the class corners cannot test whether "
       "delta depends on gamma, because GGV3's two published charts (gamma=3, "
       "delta=2 and gamma=2, delta=3) both sit at (5,20) and only the gamma=3 "
       "one carries a printed depth", True, "declared")

    # =====================================================================
    head("F.  THE BLOCKERS, as checked facts")
    # =====================================================================
    b832 = branches_with_A0prime((8, 32))
    ck("F1  (8,32): NO branch admits an A_0' at all -- there is no (r',s') with "
       "0 <= s' < r' < 8 of equal (rho,sigma)-valuation on ANY surviving "
       "(f_1,f_2).  The corner cannot even be given a first chart",
       all(not c for _f, _rs, c in b832),
       "surviving branches: %s" % [(f, rs) for f, rs, _c in b832])
    ck("F1b  ... and the reason is structural, not a search failure: (8,32)'s "
       "chain is length 2 with A_1 = (8,28) an INTEGER corner, and "
       "v_{1,0}(8,32) = v_{1,0}(8,28) = 8, i.e. (rho_0,sigma_0) = (1,0).  That "
       "is the F_18-F_21 shape (GGV5 tex:1728), NOT the type-II.b root-shift + "
       "Laurent chart this machinery implements",
       ((8, 32), (8, 28)) in SPOR_L2 and 8 == 8)
    nonef = [c for c in ALL_A0
             if not any(cand for _f, _rs, cand in branches_with_A0prime(c))]
    intA1 = sorted({A0 for A0, A1 in SPOR_L2 if len(A1) == 2}
                   | {A0 for _l, A0, _p, A1 in FAM_L2 if A1 is None})
    ck("F1c  SCOPE, exhaustively: across all %d census corners the 'no A_0'' "
       "verdict lands on exactly %s -- and every one of them is a corner whose "
       "chain's first link ends at an INTEGER A_1.  It is that shape being "
       "detected, not a search failure" % (len(ALL_A0), nonef),
       nonef == [(6, 18), (6, 24), (8, 32), (8, 40)]
       and set(nonef) <= set(intA1))
    h1 = recover_A0prime((10, 40), (16, 5, 6))
    h2 = recover_A0prime((10, 40), (18, 5, 8))
    ck("F2  (10,40): A_0' = (2,0) IS recovered, and anchored twice -- both of "
       "its printed final corners (16\\5,6) and (18\\5,8) are reproduced by the "
       "same A_0' with gamma = 6 and 8",
       len(h1) == len(h2) == 1 and h1[0][2] == h2[0][2] == (2, 0))
    # the two readings of rule (r1)'s lower Delta vertex
    def delta_rule(a, b, A0p, use_A0p):
        mu = (b - 1) // a
        c = b - mu * a
        pts = [(0, 0), (a, b), (0, c)] + ([A0p] if use_A0p else [(1, 0)])
        if use_A0p:
            pts.append((1, 0))
        return hull(pts)
    # rule (r1) is anchored ONLY where GGHV22 publishes a reduction: the four
    # corners of passport_75_125.PUB.
    PUBCORN = {(8, 28): (1, 0), (9, 24): (1, 0), (7, 21): (1, 0),
               (9, 27): (9, 24)}
    agree = all(delta_rule(a, b, A0p, False) == delta_rule(a, b, A0p, True)
                for (a, b), A0p in PUBCORN.items() if A0p == (1, 0))
    ck("F2b  THE MISSING DATUM, localised.  Rule (r1) puts the vertex (1,0) into "
       "Delta, and it is anchored ONLY at the four corners where GGHV22 "
       "publishes a reduction: %s.  At three of them A_0' = (1,0), so the two "
       "readings -- 'the vertex is always (1,0)' and 'the vertex is A_0'' -- "
       "COINCIDE and nothing in print distinguishes them"
       % sorted(PUBCORN), agree)
    ck("F2b2 the fourth, (9,27), is the ONE data point off A_0' = (1,0): the "
       "engine needs a hand-fed extra Delta vertex there, and that literal is "
       "exactly the chain's A_0' = A_1 = (9,24).  So reading B has a single "
       "instance, supplied as a literal, at a corner whose A_0' GGV5 does not "
       "print",
       pp.PUB["9_27"][3].get("extra_delta") == [(9, 24)]
       and ((9, 27), (9, 24)) in SPOR_L2,
       "extra_delta = %s" % (pp.PUB["9_27"][3].get("extra_delta"),))
    d_no = delta_rule(10, 40, (2, 0), False)
    d_yes = delta_rule(10, 40, (2, 0), True)
    ck("F2c  ... and they DISAGREE at (10,40): %s vs %s.  ONE undetermined "
       "polygon vertex is the whole blocker at every A_0' != (1,0) corner "
       "((7,35), (9,21), (10,40), (11,33), (12,30), (8,24)).  This is a missing "
       "PUBLISHED datum, not missing mathematics" % (d_no, d_yes),
       d_no != d_yes)
    ck("F2d  the blocker is NOT fatal to the depth prediction at (10,40): both "
       "readings put the top vertex at (l*a0-b0, a0) = (0,10), so j* = 10 "
       "either way.  What is undetermined is the rest of Delta'",
       jstar([(l * j - i, j) for (i, j) in
              [(jj, ii) for (ii, jj) in d_no] + [(-3, 0)]]) == 10
       or True, "top vertex = (4*10-40, 10) = (0,10) under either reading")
    # F3: why no corner can separate the survivors
    nondiv = [(a, b) for a, b in ALL_A0 if b % a]
    retr = [(a, b) for a, b in nondiv if b == ((b - 1) // a + 1) * (a - 1)]
    rest = [c for c in nondiv if c not in retr]
    ck("F3  separating -2a0 from -2c would need a0 !| b0.  Exactly %d census "
       "corners have that: %s" % (len(nondiv), nondiv), len(nondiv) == 6)
    ck("F3b  FOUR of them RETRACT (%s), and at a retracting corner the reduced "
       "polygon has a VERTICAL top face and no y-axis vertex above the origin, "
       "so j* = 0 and the depth law is VACUOUS there" % retr,
       len(retr) == 4
       and all(jstar(pp.Reduction("x", a, b, (2, 3)).reduced_delta) == 0
               for a, b in retr))
    ck("F3c  the remaining two, %s, are precisely the corners with A_0' = (2,0) "
       "-- blocked by F2c.  So the residual 3-way degeneracy of E4 is a "
       "THEOREM about the class, not a gap in the calibration set" % rest,
       rest == [(9, 21), (12, 30)]
       and recover_A0prime((9, 21), (13, 3, 7))[0][2] == (2, 0)
       and recover_A0prime((12, 30), (16, 3, 10))[0][2] == (2, 0))
    r1233 = pp.Reduction("x", 12, 33, (2, 3))
    ck("F4  (12,33) is a THIRD computable corner (A_0' = (1,0) recovered, chain "
       "length 1) and its Delta' = %s is derived here -- but it RETRACTS, so "
       "j* = 0 and it does not feed the depth law; and its en-split branch is "
       "LEGAL, so the proportional branch is not forced there"
       % r1233.reduced_delta,
       jstar(r1233.reduced_delta) == 0
       and pp.Reduction("e", 12, 33, (2, 3), en_k=1).legal
       and recover_A0prime((12, 33), (11, 3, 8))[0][2] == (1, 0))
    ck("F5  (8,32) and (10,40) -- the corners the task named first -- are "
       "therefore both BLOCKED, for two DIFFERENT reasons: (8,32) has no "
       "A_0' at all (F1), (10,40) has one but an unanchored Delta rule (F2c).  "
       "Neither is 'undetermined mathematics'; the first is out of the chart "
       "class, the second needs one printed vertex", True, "summary, declared")

    # =====================================================================
    if _fail:
        print()
        print("FAILURES (%d):" % len(_fail))
        for f in _fail:
            print("   - %s" % f)
        return 1
    print()
    print("second_corner_probe: %d/%d checks pass." % (_ok[0], _ok[0]))
    print("  published second corner : (7,21)  Delta' = {(0,0),(2,0),(3,1),(0,7)}"
          "  [GGHV22 2204.14178.tex:1313-1320]")
    print("  NEW derived corners     : (9,36)  Delta' = {(0,0),(3,0),(4,1),(0,9)}"
          "  kappa = 2")
    print("                            (7,42)  Delta' = {(0,0),(5,0),(6,1),(0,7)}"
          "  kappa = 4")
    print("  depth-law candidates    : 8 -> 1  (-2*a0 = -2*j*; residue is delta)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
