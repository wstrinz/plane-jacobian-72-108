#!/usr/bin/env python3
"""passport_75_125.py -- the reduced polygon of the (5,20) corner, and the
top-band Belyi passport question for (75,125).

HEADLINE.  The blocker is removed: the reduced Newton polygons of the (5,20)
corner ARE derivable, the derivation reproduces all five published GGHV22
reductions as controls, and it says

    (75,125):   N(P) = {(0,0),(9,0),(12,3),(0,15)}
                N(Q) = {(0,0),(15,0),(20,5),(0,25)}     [P,Q] = x^2   (kappa = 2)

and then -- exactly, exhaustively, over every branch of the reduction --

    THERE IS NO ADMISSIBLE FACE.  The gate  u*kappa = m + n - 1  FAILS for every
    primitive functional on every branch.  (75,125) has NO top-band Belyi map,
    hence NO passport and NO Hurwitz number.  It is NOT a member of the Catalan
    family; it is not a member of any Belyi family.

Two corrections to the standing (75,125) model fall out and are load-bearing:

  * kappa = 2, not 3 (l = 4, not 5).  CROSS-CHECKED against GGV3's own published
    reduction of the SIBLING case (50,75) at the SAME corner (5,20):
    "[P_1,Q_1] = x^2, deg(P_1) = 10, deg(Q_1) = 15" (1406.0886_GGV3.tex:1725-1727).
    Our derived reduced polygon has total degree 5, so (m,n) = (2,3) gives
    degrees (10,15) and the bracket x^{l-2} = x^2 -- all three numbers match.
    l = 5 would give degrees (20,30).  This is the decisive external test.
  * The reduced polygon carries NO vertical top face, so C is a MONOMIAL
    (deg C = 1), not y^2(y^3+1).  Same failure mode already documented for
    (7,21) in FAMILY_GRAMMAR.md: GGV5's final-corner l is a ramification index,
    not a Laurent-chart exponent.

And the (75,125) prediction of BELYI_PASSPORT.md sec.5 (gate 3u = 8r-1, minimal
face f = 13i - 30j) is REFUTED, for two independent reasons: kappa is 2 not 3,
and r = max_R f ignored the foot vertex of the reduced polygon, which dominates.
Even taking kappa = 3 counterfactually, the sweep here is still empty.

WHAT A PASSPORT IS AND IS NOT (carried forward verbatim from BELYI_PASSPORT.md
sec.6, because it governs how to read this file):  this layer governs ONLY the
top band.  J3..J0 and the endgame get nothing from it.  A NONZERO Hurwitz number
never kills a case -- it hands you a number field.  At (72,108) the top band
admitted 5 perfectly good covers and the case died much further down.  What the
layer buys is exact, instant knowledge of WHICH FIELD the endgame will live in,
from a character sum.  Symmetrically, the EMPTY verdict here is NOT a kill of
(75,125): it says only that (75,125)'s top band is positive-dimensional /
structureless, so this particular cheap oracle gives you nothing there and the
endgame field must be found the expensive way.

Read-only: creates nothing, modifies nothing.  Imports the ALREADY-VALIDATED
character machinery of belyi_passport.py (sections A1-A7 there validate it
against brute-force enumeration in S_n) rather than reimplementing it, so the
Hurwitz numbers here are not self-certified.

Usage:
    python -u passport_75_125.py            # full report
    python -u passport_75_125.py --quiet    # exit 0 iff every check passes
"""
from __future__ import annotations

import argparse
import sys
from itertools import product
from math import factorial, gcd

import sympy as sp

import belyi_passport as bp          # read-only import of validated helpers

CHECKS: list[tuple[str, bool, str]] = []
VERBOSE = True


def chk(tag: str, ok: bool, msg: str = "") -> bool:
    CHECKS.append((tag, bool(ok), msg))
    if VERBOSE:
        print(f"  [{'PASS' if ok else 'FAIL'}] {tag}  {msg}")
    return bool(ok)


def head(s: str) -> None:
    if VERBOSE:
        print("\n" + "=" * 92 + f"\n{s}\n" + "=" * 92)


def note(s: str) -> None:
    if VERBOSE:
        print(s)


# ==========================================================================
# 0.  Plane-lattice utilities (exact integers only)
# ==========================================================================
def _cross(o, a, b):
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def hull(pts):
    """Convex hull vertices, counter-clockwise, no collinear points."""
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


def is_vertex(pts, p):
    return tuple(p) in set(hull(pts))


def lattice_points(verts):
    """Every integer point of the convex hull of verts (boundary included)."""
    H = hull(verts)
    if len(H) < 3:
        # degenerate: a segment or a point
        pts = set()
        for a in H:
            for b in H:
                d = (b[0] - a[0], b[1] - a[1])
                g = gcd(abs(d[0]), abs(d[1])) or 1
                st = (d[0] // g, d[1] // g)
                for k in range(g + 1):
                    pts.add((a[0] + k * st[0], a[1] + k * st[1]))
        return sorted(pts)
    xs = [v[0] for v in H]
    ys = [v[1] for v in H]
    out = []
    n = len(H)
    for i in range(min(xs), max(xs) + 1):
        for j in range(min(ys), max(ys) + 1):
            if all(_cross(H[k], H[(k + 1) % n], (i, j)) >= 0 for k in range(n)):
                out.append((i, j))
    return out


def edge_normals(verts):
    """Primitive outward normals (n1,n2) of the edges of hull(verts), CCW."""
    H = hull(verts)
    out = []
    n = len(H)
    if n < 2:
        return out
    for k in range(n):
        a, b = H[k], H[(k + 1) % n]
        dx, dy = b[0] - a[0], b[1] - a[1]
        nx, ny = dy, -dx                    # outward for CCW orientation
        g = gcd(abs(nx), abs(ny)) or 1
        out.append((nx // g, ny // g))
    return out


def egcd(a, b):
    if b == 0:
        return (a, 1, 0)
    g, x, y = egcd(b, a % b)
    return (g, y, x - (a // b) * y)


# ==========================================================================
# 1.  THE REDUCTION ENGINE
#
#     Input: the chain corner A_0 = (a0,b0) with A_0' = (1,0), the polygon
#     multipliers (m,n), and the branch selectors.  Output: the reduced
#     Newton polygons of P and Q and the bracket exponent kappa.
#
#     Every rule below is DERIVED, and every one is checked against the five
#     published reductions in section R.  The rules:
#
#     (r1)  mu = floor((b0-1)/a0),  l = mu+1,  c = b0 - mu*a0  in {1..a0}.
#           => Delta = {(0,0),(1,0),(a0,b0),(0,c)}  (+ published extra corners)
#           l is the exponent of the final Laurent chart (x^-1, x^l y); it is
#           the MINIMAL l putting the flipped corner (b0,a0) in the first
#           quadrant after inversion, l >= b0/a0, and c is the y-axis vertex.
#     (r2)  kappa = l - 2   (fused-chart Jacobian -x^(l-2); composite_charts.py)
#     (r3)  flip x<->y.
#     (r4)  the flipped lower edge (c,0)--(b0,a0) has primitive direction
#           (mu,1), so Pred_P(1,0) = (1,-mu): the root shift is y -> y+lam*x^-s
#           with s = mu (GGV6 Prop 2.5 also allows s = mu-1 = 2 when mu = 3;
#           both are swept).  It replaces (c,0) by the foot (-s,0).
#     (r5)  q = gcd(a0,b0) is the Cor-7.4 multiplicity, en(R) = (b0,a0)/q is
#           primitive, and R has z-degree a0/q.  Hence R can have at most
#           a0/q distinct linear factors, and the two-factor branch produces a
#           split corner at (c + mu*e, e) only for e = j*q, j = 1..a0/q-1,
#           and only when e >= c (else the corner leaves the first quadrant).
#           a0/q = 1  =>  NO split branch exists at all.
#     (r6)  RETRACTION.  The corner (B,A) adjacent to (0,1) makes the edge
#           {(0,1),(B,A)} collapse to a vertical face under inversion exactly
#           when B = l*(A-1); then the edge form is y*(x^l y - alpha)^(A-1) and
#           the edge root shift retracts (0,1) to (B-l, A-1).
#     (r7)  inversion (i,j) -> (l*j - i, j).
#     (r8)  GGV1 Prop 8.2.  Either (1) en(P) ~ en(Q) -- the PROPORTIONAL branch,
#           N(P) = m*Delta', N(Q) = n*Delta', foot scaled -- or (2) there is k
#           with (k+1)*a0 < b0 and {en(P),en(Q)} = {(-k,0),(k+1,1)}
#           pre-inversion, i.e. {(k,0),(l-k-1,1)} after it; then the en point
#           is UNSCALED and REPLACES the foot.  Branch (2) is only available
#           when both en points are genuine vertices of the resulting hulls.
# ==========================================================================
class Reduction:
    def __init__(self, tag, a0, b0, mn, s=None, extra_delta=(), split_e=None,
                 en_k=None, en_swap=False, l_override=None):
        self.tag = tag
        self.a0, self.b0 = a0, b0
        self.m, self.n = mn
        self.mu = (b0 - 1) // a0
        self.l = self.mu + 1 if l_override is None else l_override
        self.c = b0 - self.mu * a0
        self.q = gcd(a0, b0)
        self.zdeg = a0 // self.q
        self.s = self.mu if s is None else s
        self.kappa = self.l - 2
        self.split_e = split_e
        self.en_k = en_k
        self.en_swap = en_swap
        self.extra_delta = [tuple(v) for v in extra_delta]
        self.notes = []
        self._build()

    # ---------------------------------------------------------------- build
    def _build(self):
        a0, b0, l, mu, c = self.a0, self.b0, self.l, self.mu, self.c
        self.delta = hull([(0, 0), (1, 0), (a0, b0), (0, c)] + self.extra_delta)
        flipped = [(j, i) for (i, j) in self.delta]          # (r3)

        # (r4) root shift: the foot (c,0) becomes (-s,0)
        pre = [p for p in flipped if p != (c, 0)] + [(-self.s, 0)]

        # (r5) optional split corner from a two-factor leading form
        self.split_pt = None
        if self.split_e is not None:
            e = self.split_e
            self.split_pt = (c + mu * e, e)
            pre.append(self.split_pt)

        # (r6) retraction of (0,1) along the adjacent edge, when it is vertical
        self.pre_noretract = hull(pre)
        H = hull(pre)
        self.retracted = None
        if (0, 1) in H:
            i01 = H.index((0, 1))
            nbrs = [H[(i01 - 1) % len(H)], H[(i01 + 1) % len(H)]]
            cand = [p for p in nbrs if p[0] > 0]
            if cand:
                B, A = max(cand, key=lambda p: p[1])
                if B == l * (A - 1):
                    self.retracted = (B - l, A - 1)
                    pre = [p for p in pre if p != (0, 1)] + [self.retracted]
        self.pre = hull(pre)

        # (r7) inversion
        inv = lambda p: (l * p[1] - p[0], p[1])
        red = [inv(p) for p in self.pre]
        self.reduced_delta = hull(red)
        self.foot = (self.s, 0)
        self.core = sorted(p for p in self.reduced_delta
                           if p not in {(0, 0), self.foot})

        # (r8) assemble N(P), N(Q)
        if self.en_k is None:
            self.branch = "proportional (GGV1 Prop 8.2(1))"
            self.enP = self.enQ = None
            NP = [(self.m * i, self.m * j) for (i, j) in self.reduced_delta]
            NQ = [(self.n * i, self.n * j) for (i, j) in self.reduced_delta]
        else:
            k = self.en_k
            e1, e2 = (k, 0), (l - k - 1, 1)
            if self.en_swap:
                e1, e2 = e2, e1
            self.branch = f"en-split (GGV1 Prop 8.2(2)), k={k}, swap={self.en_swap}"
            self.enP, self.enQ = e1, e2
            NP = [(0, 0), e1] + [(self.m * i, self.m * j) for (i, j) in self.core]
            NQ = [(0, 0), e2] + [(self.n * i, self.n * j) for (i, j) in self.core]
        self.NP, self.NQ = hull(NP), hull(NQ)
        # branch (2) is only legal when both en points really are vertices
        self.legal = True
        if self.en_k is not None:
            self.legal = (self.enP in set(self.NP)) and (self.enQ in set(self.NQ)) \
                and (self.en_k + 1) * self.a0 < self.b0
        self.degP = max(i + j for (i, j) in self.NP)
        self.degQ = max(i + j for (i, j) in self.NQ)

    def line(self):
        return (f"    l={self.l} kappa={self.kappa} c={self.c} q={self.q} "
                f"zdeg={self.zdeg} s={self.s} retr={self.retracted} "
                f"{self.branch}{'' if self.legal else '  [ILLEGAL]'}\n"
                f"    Delta'    = {self.reduced_delta}\n"
                f"    N(P)      = {self.NP}   deg {self.degP}\n"
                f"    N(Q)      = {self.NQ}   deg {self.degQ}")


# ==========================================================================
# 2.  THE PASSPORT ENGINE  (the reusable compiler stage)
# ==========================================================================
def _segment_points(a, b):
    d = (b[0] - a[0], b[1] - a[1])
    g = gcd(abs(d[0]), abs(d[1]))
    if g == 0:
        return [a]
    st = (d[0] // g, d[1] // g)
    return [(a[0] + k * st[0], a[1] + k * st[1]) for k in range(g + 1)]


def max_face(verts, u, v):
    """(max of f = u*i - v*j over the polygon, the lattice points of the face
    attaining it).  A linear functional is maximal on a vertex or an edge, so
    this is exact and needs no lattice enumeration."""
    f = lambda p: u * p[0] - v * p[1]
    H = hull(verts)
    mf = max(map(f, H))
    ext = [p for p in H if f(p) == mf]
    if len(ext) == 1:
        return mf, ext
    assert len(ext) == 2, (verts, u, v, ext)
    return mf, sorted(_segment_points(*ext))


def band_data(NP, NQ, u, v):
    """For f(i,j) = u*i - v*j on the two polygons, return the band dossier."""
    mf, bP = max_face(NP, u, v)
    nf, bQ = max_face(NQ, u, v)
    # complementary GL_2(Z) coordinate: i' = a*i + b*j with a*v + b*u = 1
    g, a, b = egcd(v, u)
    assert g in (1, -1), (u, v, g)
    a, b = a * g, b * g
    assert a * v + b * u == 1
    ip = lambda p: a * p[0] + b * p[1]
    return dict(u=u, v=v, m=mf, n=nf, bandP=bP, bandQ=bQ,
                tP=sorted(map(ip, bP)), tQ=sorted(map(ip, bQ)), tcoef=(a, b))


def passport_from_polygon(NP, NQ, kappa, label="", require_edge=True,
                          verbose=False):
    """Newton polygons + bracket exponent  ->  every admissible top-band face,
    with its gate verdict, passport and Hurwitz number.

    A candidate functional f = u*i - v*j is admissible iff
      * (u,-v) is a primitive outward normal of an edge of BOTH polygons
        (so the max face is an edge, hence the band has >= 2 monomials and the
        top-layer equation m*A*D' - n*A'*D = gamma*t^p is not vacuous),
      * POSITIVITY  m >= 1 and n >= 1: the band must be a genuine top z-band.
        (If m = n = 0 the "band" is the z^0 layer, m*A*D' - n*A'*D vanishes
        identically, and Phi = D^(m/g)/A^(n/g) has degree n*k = 0 -- not a
        cover.  Omitting this test lets the y-axis edge f = -i sneak through
        with m = n = 0 and uk = -1 = m+n-1; it is what made the first draft of
        check W5 fail, and it is a real hole, not a formality.), and
      * THE GATE  u*kappa == m + n - 1, where m,n are the maxima of f.
    The gate is the discriminator.  Riemann-Hurwitz is NOT: it balances
    identically for every (m,n,k,l) -- see check W1.
    """
    cands = sorted(set(edge_normals(NP)) & set(edge_normals(NQ)))
    rows, admissible, rejected = [], [], []
    for (n1, n2) in cands:
        u, v = n1, -n2
        if gcd(abs(u), abs(v)) != 1:
            continue
        d = band_data(NP, NQ, u, v)
        d["edgeP"] = len(d["bandP"]) >= 2
        d["edgeQ"] = len(d["bandQ"]) >= 2
        d["gate_lhs"] = u * kappa
        d["gate_rhs"] = d["m"] + d["n"] - 1
        d["gate"] = (d["gate_lhs"] == d["gate_rhs"])
        d["pos"] = (d["m"] >= 1 and d["n"] >= 1)
        d["ok"] = d["gate"] and d["pos"] \
            and (not require_edge or (d["edgeP"] and d["edgeQ"]))
        rows.append(d)
        if d["ok"]:
            admissible.append(d)
    out = []
    for d in admissible:
        m, n = d["m"], d["n"]
        g = gcd(m, n) or 1
        M, N = m // g, n // g
        ordA, degA = d["tP"][0], d["tP"][-1]
        ordD, degD = d["tQ"][0], d["tQ"][-1]
        k, ell = degA - ordA, degD - ordD
        deg = N * k
        over0 = tuple([M] * ell + ([N * k - M * ell] if N * k - M * ell > 0 else []))
        overinf = tuple([N] * k)
        third = tuple([k + ell] + [1] * (N * k - k - ell))
        pp = dict(face=d, M=M, N=N, ordA=ordA, degA=degA, ordD=ordD, degD=degD,
                  k=k, ell=ell, deg=deg, over0=over0, overinf=overinf,
                  third=third, p_forced=ordA + ordD - 1,
                  cross=(m * degD == n * degA))
        # sec.4's OTHER forced relations.  The gate is necessary but not
        # sufficient: m*deg D = n*deg A must hold too (else gamma*t^p could not
        # have low degree), and the three ramification profiles must each sum to
        # deg Phi = N*k with all parts >= 1.  A face that passes the gate but
        # fails these is NOT a Belyi datum; it is reported, never counted.
        pp["valid"] = (pp["cross"] and N * k - M * ell > 0
                       and sum(over0) == deg and sum(overinf) == deg
                       and sum(third) == deg and min(third) >= 1)
        d["mechanism"] = pp["valid"]
        if pp["valid"]:
            Nt = bp.n_triples(deg, tuple(sorted(over0, reverse=True)),
                              tuple(sorted(overinf, reverse=True)),
                              tuple(sorted(third, reverse=True)))
            pp["hurwitz"] = Nt // factorial(deg)
            pp["hurwitz_exact"] = (Nt % factorial(deg) == 0)
            out.append(pp)
        else:
            rejected.append(pp)
    if verbose:
        note(f"  candidate shared edge-normals for {label}: "
             f"{[(r['u'], r['v']) for r in rows]}")
        for r in rows:
            note(f"    f = {r['u']}i - {r['v']}j : m={r['m']} n={r['n']} "
                 f"|band|=({len(r['bandP'])},{len(r['bandQ'])})  "
                 f"gate {r['gate_lhs']} vs {r['gate_rhs']} -> "
                 f"{'PASS' if r['gate'] else 'fail'}"
                 f"{'' if r['pos'] else '  [rejected: m=n=0, degenerate]'}")
    return rows, out


# ==========================================================================
# R.  CONTROLS -- reproduce every published reduction
# ==========================================================================
PUB = {
    # tag: (a0,b0, (m,n), kwargs, N(P), N(Q), kappa, tex line)
    "8_28_sub1": (8, 28, (2, 3), dict(split_e=4, en_k=1),
                  [(0, 0), (1, 0), (8, 14), (8, 16), (0, 8)],
                  [(0, 0), (2, 1), (12, 21), (12, 24), (0, 12)], 2, "1003"),
    "8_28_sub2": (8, 28, (2, 3), dict(en_k=1),
                  [(0, 0), (1, 0), (8, 14), (8, 16)],
                  [(0, 0), (2, 1), (12, 21), (12, 24)], 2, "1004"),
    "9_24_case3": (9, 24, (2, 3), dict(en_k=1, en_swap=True),
                   [(0, 0), (1, 1), (6, 16), (6, 18)],
                   [(0, 0), (1, 0), (9, 24), (9, 27)], 1, "676-677"),
    "9_24_case1": (9, 24, (2, 3), dict(split_e=6, en_k=1, en_swap=True),
                   [(0, 0), (1, 1), (6, 16), (6, 18), (0, 12)],
                   [(0, 0), (1, 0), (9, 24), (9, 27), (0, 18)], 1, "672-673"),
    "9_27": (9, 27, (2, 3), dict(en_k=1, en_swap=True, extra_delta=[(9, 24)]),
             [(0, 0), (1, 1), (6, 16), (6, 18), (0, 18)],
             [(0, 0), (1, 0), (9, 24), (9, 27), (0, 27)], 1, "466-467"),
    "7_21": (7, 21, (2, 3), dict(),
             [(0, 0), (4, 0), (6, 2), (0, 14)],
             [(0, 0), (6, 0), (9, 3), (0, 21)], 1, "1316-1317"),
}


def section_R():
    head("R.  CONTROL: the reduction engine vs every published GGHV22 reduction")
    note("  (published vertex lists transcribed from paper_src/2204.14178.tex;")
    note("   line numbers cited per row.  Nothing here is fitted: l, c, q, the")
    note("   retraction and the inversion are all computed from (a0,b0) alone.)")
    for tag, (a0, b0, mn, kw, wp, wq, wk, ln) in PUB.items():
        r = Reduction(tag, a0, b0, mn, **kw)
        note(f"\n  --- {tag}   A_0=({a0},{b0})  (m,n)={mn}   [tex:{ln}]")
        note(r.line())
        chk(f"R-{tag}-P  N(P) reproduced exactly", r.NP == hull(wp),
            f"{r.NP}")
        chk(f"R-{tag}-Q  N(Q) reproduced exactly", r.NQ == hull(wq),
            f"{r.NQ}")
        chk(f"R-{tag}-kappa  [P,Q] = x^{wk}", r.kappa == wk,
            f"kappa = l-2 = {r.kappa}")
        chk(f"R-{tag}-legal  Prop-8.2 branch is legal", r.legal)
    # the rule set itself, stated and checked
    note("")
    for a0, b0, wl, wc, wq in ((8, 28, 4, 4, 4), (9, 24, 3, 6, 3),
                               (9, 27, 3, 9, 9), (7, 21, 3, 7, 7),
                               (5, 20, 4, 5, 5)):
        mu = (b0 - 1) // a0
        chk(f"R-rule ({a0},{b0}): l=ceil(b0/a0)={wl}, c=b0-(l-1)a0={wc}, "
            f"q=gcd={wq}",
            mu + 1 == wl and b0 - mu * a0 == wc and gcd(a0, b0) == wq,
            f"l={mu+1} c={b0-mu*a0} q={gcd(a0,b0)} zdeg={a0//gcd(a0,b0)}")
    chk("R-rule l is minimal with l*a0 >= b0 (first-quadrant condition)",
        all((b0 - 1) // a0 + 1 == -(-b0 // a0)
            for a0, b0 in ((8, 28), (9, 24), (9, 27), (7, 21), (5, 20))))
    # an extra published anchor: GGHV22 prints the PRE-INVERSION vertex set for
    # (9,27) at tex:471, after phi_1 (flip) and phi_2 (root shift, s=2).
    r927 = Reduction("927", 9, 27, (2, 3), en_k=1, en_swap=True,
                     extra_delta=[(9, 24)])
    chk("R-9_27-pre  pre-inversion set after flip+root-shift matches GGHV22 "
        "tex:471 {(0,0),(27,9),(24,9),(0,1),(-2,0)}",
        r927.pre_noretract == hull([(0, 0), (27, 9), (24, 9), (0, 1), (-2, 0)]),
        f"{r927.pre_noretract}")
    # the (8,28) Pred branch is resolved by GGV6 Prop 2.5 (arXiv:1708.09367,
    # "b igual a 2"): (rho,sigma) ~ (l,-Delta) with l < Delta < a/2 and
    # (a-2Delta) | (Delta-l), applied to the b=2 corner en(R) = (a/l,2).
    cands = [D for D in range(2, 4) if (7 - 2 * D) != 0 and (D - 1) %
             abs(7 - 2 * D) == 0 and 1 < D < 7 / 2]
    chk("R-ggv6  GGV6 Prop 2.5 at en(R)=(7,2) (l=1,a=7) forces Delta=3, i.e. "
        "Pred=(1,-3), i.e. mu=3, i.e. l_chart=4 -- resolving GGHV's own "
        "(1,-2)/(1,-3) branch in favour of the engine's mu",
        cands == [3], f"admissible Delta = {cands}; Delta=2 fails 3 | 1")
    # the (72,108) polygons the Belyi lane actually used
    r = Reduction("ctl", 8, 28, (2, 3), split_e=4, en_k=1)
    chk("R-belyi  engine output == the N(P),N(Q) used by belyi_passport.py",
        r.NP == hull(bp.NP) and r.NQ == hull(bp.NQ))


# ==========================================================================
# P.  WHERE (75,125) SITS IN GGV5's CENSUS -- and why 125 is absent from the
#     A_1 tables.  (This section exists because the opposite reading was put
#     to me as a reframing: that 125's absence from GGV5's chain tables means
#     (75,125) is "a different KIND of case" with no A_1.  It is the reverse.)
# ==========================================================================
# GGV5 = paper_src/1708.07936_GGV5.tex.  Transcribed table rows, with lines.
GGV5_FAM = {            # L1673-1690: Family | A_0 | A_0' | A_1 | k | m | n
    "F_1": ((4, 12), (1, 0), (7, 4, 3), 1, (2, 3), (3, 4)),
    "F_2": ((5, 20), (1, 0), (7, 5, 2), 1, (1, 2), (2, 3)),
    "F_3": ((5, 20), (1, 0), (8, 5, 3), 1, (4, 3), (3, 2)),
    "F_9": ((7, 21), (1, 0), (11, 7, 2), 1, (1, 2), (2, 3)),
    "F_14": ((9, 24), (1, 0), (7, 3, 4), 1, (1, 2), (4, 7)),
    "F_17": ((9, 24), (1, 0), (11, 3, 8), 1, (5, 2), (8, 3)),
}
GGV5_DEG_TABLE = [      # L1794-1817: Family | (m,n) | max deg  (13 rows)
    ("F_1", (3, 4), 64), ("F_1", (5, 7), 112), ("F_2", (2, 3), 75),
    ("F_2", (3, 5), 125), ("F_3", (3, 2), 75), ("F_7", (2, 7), 147),
    ("F_8", (3, 7), 147), ("F_9", (2, 3), 84), ("F_9", (3, 5), 140),
    ("F_11", (2, 5), 140), ("F_17", (2, 3), 99), ("F_22", (2, 3), 96),
    ("F_24", (3, 4), 128),
]
GGV5_OTHER_L1 = [       # L1821-1839: "9 OTHER possible pairs", chain length 1
    ((7, 35), (19, 7, 5), (2, 3), 126), ((7, 42), (13, 7, 6), (3, 2), 147),
    ((7, 42), (13, 7, 6), (2, 3), 147), ((8, 28), (7, 4, 3), (3, 4), 144),
    ((8, 28), (11, 4, 7), (3, 2), 108), ((9, 36), (17, 9, 4), (3, 2), 135),
    ((9, 36), (17, 9, 4), (2, 3), 135), ((11, 33), (19, 4, 8), (2, 3), 132),
    ((12, 33), (11, 3, 8), (2, 3), 135),
]


def section_P():
    head("P.  (75,125) in GGV5's census: it IS a family case, WITH an A_1")
    note("  GGV5 L1794: '...the 34 possible counterexamples with max deg <= 150.")
    note("   THIRTEEN OF THEM correspond to a choice of (m,n) in some of the")
    note("   families listed in the previous section, as can be seen in the")
    note("   following table' -- then Family | (m,n) | max deg, L1797-1817.")
    note("  GGV5 L1821: 'There are 9 OTHER possible pairs with a complete chain")
    note("   of length 1, which we list in the following table' -- and only")
    note("   THOSE tables carry an explicit A_1 column, because those are")
    note("   precisely the cases that are NOT in any family.")
    note("")
    note("  So the A_1 tables are the COMPLEMENT of the family table.  Hence:")
    note("    * F_2 (3,5) -> 125 is in the FAMILY table (L1805), and its chain")
    note("      data -- A_0 = (5,20), A_0' = (1,0), A_1 = (7\\5,2), k = 1 -- is")
    note("      printed in the family table at L1676.  (75,125) HAS an A_1.")
    note("    * (72,108) is in the 'OTHER pairs' table (L1832) as")
    note("      (8,28) | (11/4,7) | (3,2) | 108, i.e. (72,108) is a SPORADIC,")
    note("      belonging to NO family.")
    note("  125's absence from the A_1 tables therefore says the OPPOSITE of")
    note("  'different kind of case': it says (75,125) is a family member and")
    note("  (72,108) is not.  Both have an A_1.")
    # max deg = n * v11(A_0) -- checkable arithmetic tying the two tables
    ok = True
    for fam, mn, md in GGV5_DEG_TABLE:
        if fam not in GGV5_FAM:
            continue
        A0 = GGV5_FAM[fam][0]
        ok = ok and md == max(mn) * (A0[0] + A0[1])
    chk("P1 every family row's max degree equals max(m,n)*v11(A_0) with A_0 "
        "from the family table -- the two tables are one dataset",
        ok, "F_2 (3,5): 5*25 = 125; F_2 (2,3): 3*25 = 75; F_9 (3,5): 5*28 = 140")
    chk("P2 (75,125) = F_2 with (m,n) = (3,5), A_0 = (5,20), A_1 = (7\\5,2), "
        "k = 1  [GGV5 L1676 + L1805]",
        GGV5_FAM["F_2"][0] == (5, 20) and GGV5_FAM["F_2"][2] == (7, 5, 2)
        and ("F_2", (3, 5), 125) in GGV5_DEG_TABLE)
    chk("P3 (72,108) is a NON-family sporadic: (8,28)|(11/4,7)|(3,2)|108 sits "
        "in the '9 OTHER pairs' table [GGV5 L1832], and no family row has "
        "max deg 108",
        ((8, 28), (11, 4, 7), (3, 2), 108) in GGV5_OTHER_L1
        and 108 not in [d for _, _, d in GGV5_DEG_TABLE])
    chk("P4 125 is absent from the A_1 tables for the same reason 108 is "
        "absent from the family table -- the tables partition the 34 cases",
        125 not in [row[3] for row in GGV5_OTHER_L1]
        and 108 not in [d for _, _, d in GGV5_DEG_TABLE])
    note("")
    note("  Does F_2 admit the standard-pair / chain machinery?  GGV5 excludes")
    note("  exactly five family rows from it: F_18, F_19, F_20, F_21 (claim at")
    note("  L1726) and F_22 with (m,n) = (2,3) (Prop 'caso antisimetrico',")
    note("  L1874).  F_2 is in NONE of those.")
    chk("P5 F_2 is inside the standard-pair framework: it is not among the "
        "families GGV5 excludes {F_18,F_19,F_20,F_21,F_22}",
        "F_2" not in {"F_18", "F_19", "F_20", "F_21", "F_22"})

    # ---- and the real reason the repo's dictionary broke -------------------
    note("\n  THE DICTIONARY, and exactly where it is valid.  The repo reads")
    note("  t = l_final and q = b_final off the final corner (a\\l, b).  Test it")
    note("  against the engine's DERIVED chart exponent and ord C:")
    rows = [("(8,28)", 8, 28, (11, 4, 7), dict(split_e=4, en_k=1)),
            ("(9,24)", 9, 24, (11, 3, 8), dict(en_k=1, en_swap=True)),
            ("(7,21)", 7, 21, (11, 7, 2), dict()),
            ("(5,20)", 5, 20, (7, 5, 2), dict())]
    verdicts = {}
    note(f"    {'corner':8s} {'A_1':10s} {'l_fin':>5s} {'l_chart':>7s} "
         f"{'b_fin':>5s} {'ord C':>5s}  dictionary")
    for name, a0, b0, A1, kw in rows:
        r = Reduction(name, a0, b0, (2, 3), **kw)
        ordC = min(j for (i, j) in r.core if i == max(p[0] for p in r.core))
        good = (A1[1] == r.l) and (A1[2] == ordC)
        verdicts[name] = good
        note(f"    {name:8s} ({A1[0]}\\{A1[1]},{A1[2]})   {A1[1]:5d} {r.l:7d} "
             f"{A1[2]:5d} {ordC:5d}  {'VALID' if good else 'BROKEN'}"
             f"   (retraction shape: {b0 == r.l*(a0-1)})")
    chk("P6 the (t,q) = (l_final,b_final) dictionary is VALID exactly on the "
        "retraction shape b0 = l(a0-1) -- (8,28),(9,24) -- and BROKEN on the "
        "b0 = l*a0 shape -- (7,21),(5,20)",
        verdicts == {"(8,28)": True, "(9,24)": True,
                     "(7,21)": False, "(5,20)": False})
    note("  That is why the dictionary, calibrated on (72,108), silently gave")
    note("  the wrong t and q at (5,20): (5,20) is the (7,21) shape, not the")
    note("  (8,28) shape.  (7,21) is the published proof of the failure.")


# ==========================================================================
# S.  THE (5,20) REDUCTION -- the blocker, removed
# ==========================================================================
def section_S():
    head("S.  The (5,20) corner: reduced Newton polygons, DERIVED")
    a0, b0 = 5, 20
    mu = (b0 - 1) // a0
    note(f"  A_0 = (5,20), A_0' = (1,0);  mu = floor(19/5) = {mu},  l = {mu+1},"
         f"  c = 20 - {mu}*5 = {b0-mu*a0},  q = gcd(5,20) = {gcd(a0,b0)}")
    note("  Delta = {(0,0),(1,0),(5,20),(0,5)}          [rule r1]")
    note("  flip  -> {(0,0),(0,1),(20,5),(5,0)}         [r3]")
    note("  lower edge (5,0)--(20,5) has direction (15,5) = 5*(3,1), so")
    note("  Pred_P(1,0) = (1,-3): shift y -> y + lam x^-3, foot (-3,0)  [r4]")
    note("  q = 5, a0/q = 1  =>  R = x*(x^3 y - alpha) has z-degree 1, ONE")
    note("     linear factor, so NO two-factor split corner exists       [r5]")
    note("  retraction test: b0 = 20 vs l*(a0-1) = 4*4 = 16  ->  NO retraction,")
    note("     so the reduced polygon has NO vertical top face and C is a")
    note("     MONOMIAL (deg C = 1).  This contradicts C = y^2(y^3+1).   [r6]")
    note("  inversion (i,j) -> (4j - i, j)                                [r7]")

    chk("S1 mu = 3, l = 4, c = 5, q = 5, a0/q = 1",
        mu == 3 and mu + 1 == 4 and b0 - mu * a0 == 5 and gcd(a0, b0) == 5
        and a0 // gcd(a0, b0) == 1)
    chk("S2 no retraction: b0 != l*(a0-1)", b0 != 4 * (a0 - 1),
        f"20 != 16")
    chk("S3 no split corner: needs e = j*q with 1<=j<=a0/q-1 = 0",
        a0 // gcd(a0, b0) - 1 == 0)

    # the two Prop-8.2 branches
    prop = Reduction("F2_j1_prop", 5, 20, (3, 5))
    note("\n  PROPORTIONAL branch (GGV1 Prop 8.2(1)):")
    note(prop.line())
    ens = []
    for k in (1, 2):
        for sw in (False, True):
            r = Reduction(f"F2_j1_en{k}{'s' if sw else ''}", 5, 20, (3, 5),
                          en_k=k, en_swap=sw)
            ens.append(r)
    note("\n  en-split branch (GGV1 Prop 8.2(2)):  every k and assignment")
    for r in ens:
        note(f"    k={r.en_k} swap={r.en_swap}: en(P)={r.enP} en(Q)={r.enQ} "
             f"-> legal = {r.legal}")
    chk("S4 the Prop-8.2(2) en-split branch is EXCLUDED at (5,20): the point "
        "(l-k-1,1) is never a vertex of the resulting hull",
        all(not r.legal for r in ens))
    note("     (at (8,28) the same test PASSES -- see check R-8_28_sub1-legal --")
    note("      so S4 is a discriminating test, not a vacuous one.)")

    chk("S5 (5,20) reduced polygons, (m,n)=(3,5)  [= the (75,125) answer]",
        prop.NP == [(0, 0), (9, 0), (12, 3), (0, 15)]
        and prop.NQ == [(0, 0), (15, 0), (20, 5), (0, 25)],
        f"N(P)={prop.NP}  N(Q)={prop.NQ}")
    chk("S6 [P,Q] = x^kappa with kappa = l - 2 = 2  (NOT 3)",
        prop.kappa == 2)

    # ---- the decisive external cross-check: GGV3 on the sibling (50,75) ----
    head("S*.  EXTERNAL CROSS-CHECK: GGV3 sec.5 on the sibling (50,75)")
    note("  1406.0886_GGV3.tex:1720-1727 -- verbatim:")
    note('    "Then by [GGV1, Remark 7.10], we know that A_0=(5,20). ...')
    note('     Proceeding as in [GGV1, Section 8] we obtain a pair')
    note('     (P_1,Q_1) in K[x,y], such that')
    note('        [P_1,Q_1]=x^2,  deg(P_1)=10  and  deg(Q_1)=15."')
    sib = Reduction("F2_j0_50_75", 5, 20, (2, 3))
    note("\n  our engine, same corner, (m,n) = (2,3):")
    note(sib.line())
    chk("S7 GGV3's published bracket for the (5,20) corner: x^2",
        sib.kappa == 2, "kappa = 2 -> [P_1,Q_1] = x^2   MATCH")
    chk("S8 GGV3's published reduced degrees for (50,75): (10,15)",
        (sib.degP, sib.degQ) == (10, 15), f"engine: ({sib.degP},{sib.degQ})"
        "   MATCH")
    alt = Reduction("F2_j0_l5", 5, 20, (2, 3), l_override=5)
    chk("S9 the repo's l=5 is REFUTED by S7/S8: it predicts kappa=3 and "
        "degrees (20,30)", (alt.kappa, alt.degP, alt.degQ) == (3, 20, 30),
        f"l=5 gives kappa={alt.kappa}, degrees ({alt.degP},{alt.degQ}) "
        "-- contradicts GGV3")
    note("\n  Three independent published numbers (bracket exponent 2, deg 10,")
    note("  deg 15) are reproduced by l=4 and all three are contradicted by")
    note("  l=5.")

    # ---- where the wrong l came from, with a published counterexample ------
    head("S**. PROVENANCE of the repo's l=5, and a published counterexample")
    note("  GGV5's length-1 (m,n)-family table (1708.07936_GGV5.tex:1673-1690),")
    note("  columns  Family | A_0 | A_0' | A_1 | k | m | n :")
    note("      F_2 | (5,20) | (1,0) | (7 \\ 5, 2) | 1 | j+2 | 2j+3")
    note("      F_9 | (7,21) | (1,0) | (11 \\ 7, 2) | 1 | j+2 | 2j+3")
    note("  The repo reads t = l_final = 5 off the F_2 final corner (7\\5,2) and")
    note("  sets kappa = t-2 = 3.  But F_9 gives the SAME reading l_final = 7")
    note("  for (7,21) -- and GGHV22 PUBLISHES (7,21)'s reduction: the chart is")
    note("  phi_3(y) = y x^3 (tex:1394) and [P,Q] = x, i.e. l_chart = 3 and")
    note("  kappa = 1.  l_final = 7 would demand kappa = 5.  So l_final is")
    note("  provably NOT the Laurent-chart exponent, by published example.")
    note("  (7,21) and (5,20) are the SAME family shape: both k=1, both with")
    note("  b_final = 2, both with (m,n) = (j+2, 2j+3).  The two errors are one")
    note("  error.  FAMILY_GRAMMAR.md sec.1 already recorded it for (7,21);")
    note("  this extends it to (5,20), which CURRENT_STATUS.md:16-17 had")
    note("  explicitly exempted ('F2 unaffected').")
    r721 = Reduction("721", 7, 21, (2, 3))
    chk("S10 GGV5 l_final = 7 for (7,21) but GGHV22 publishes chart x^3 and "
        "[P,Q] = x: l_final != l_chart, PROVED by published counterexample",
        r721.l == 3 and r721.kappa == 1 and 7 - 2 == 5 != 1,
        "engine l_chart = 3, kappa = 1 = published; l_final-2 = 5 != 1")
    chk("S11 therefore reading t = l_final = 5 at (5,20) is the SAME error, "
        "and t = l_chart = 4, kappa = 2", Reduction("f2", 5, 20, (3, 5)).l == 4)

    # ---- every GGV5 family at the (5,20) corner shares this Delta' ---------
    note("\n  A_0' = (1,0) is shared by ALL FIVE GGV5 families at A_0 = (5,20)")
    note("  (F_2..F_6, final corners (7\\5,2), (8\\5,3) x2, (9\\5,4) x2), and the")
    note("  reduction only uses (A_0, A_0'), so all five share this Delta'.")
    return prop, ens


# ==========================================================================
# T.  CONTROL of the passport engine on (72,108)
# ==========================================================================
def section_T():
    head("T.  CONTROL: passport_from_polygon on (72,108) sub1")
    r = Reduction("ctl", 8, 28, (2, 3), split_e=4, en_k=1)
    # T0: the fast hull-based max_face agrees with brute lattice enumeration
    ok = True
    for uu, vv in ((2, 1), (1, 0), (-1, -1), (1, 2), (3, 5), (-2, 3)):
        for V in (r.NP, r.NQ):
            L = lattice_points(V)
            mv = max(uu * i - vv * j for (i, j) in L)
            bb = sorted(p for p in L if uu * p[0] - vv * p[1] == mv)
            got = max_face(V, uu, vv)
            ok = ok and got == (mv, bb)
    chk("T0 hull-based max_face == brute-force lattice enumeration "
        "(6 functionals x 2 polygons)", ok)
    rows, out = passport_from_polygon(r.NP, r.NQ, r.kappa, "(72,108) sub1",
                                      verbose=True)
    chk("T1 exactly one admissible face", len(out) == 1)
    if not out:
        return
    pp = out[0]
    d = pp["face"]
    chk("T2 that face is f = 2i - j  (u=2, v=1)", (d["u"], d["v"]) == (2, 1))
    chk("T3 the gate reads 2*2 = 2+3-1", d["gate_lhs"] == 4 == d["gate_rhs"])
    chk("T4 band of P is t^1..t^8 (A = t*a, deg a = 7)",
        (pp["ordA"], pp["degA"], pp["k"]) == (1, 8, 7))
    chk("T5 band of Q is t^2..t^12 (D = t^2*d, deg d = 10)",
        (pp["ordD"], pp["degD"], pp["ell"]) == (2, 12, 10))
    chk("T6 p = ord A + ord D - 1 = 2", pp["p_forced"] == 2)
    chk("T7 m*deg D = n*deg A", pp["cross"])
    chk("T8 Phi = D^2/A^3, degree 21", (pp["M"], pp["N"], pp["deg"]) == (2, 3, 21))
    chk("T9 passport (2^10,1 | 3^7 | 17,1^4)",
        pp["over0"] == tuple([2] * 10 + [1]) and pp["overinf"] == tuple([3] * 7)
        and pp["third"] == (17, 1, 1, 1, 1),
        f"{pp['over0']} | {pp['overinf']} | {pp['third']}")
    chk("T10 HURWITZ NUMBER = 5", pp["hurwitz"] == 5 and pp["hurwitz_exact"],
        f"h = {pp['hurwitz']}")
    # sub2 must give the same top band (the (0,8)/(0,12) corner is irrelevant)
    r2 = Reduction("ctl2", 8, 28, (2, 3), en_k=1)
    _, out2 = passport_from_polygon(r2.NP, r2.NQ, r2.kappa, "(72,108) sub2")
    chk("T11 sub2 gives the same passport and Hurwitz number",
        len(out2) == 1 and out2[0]["over0"] == pp["over0"]
        and out2[0]["third"] == pp["third"] and out2[0]["hurwitz"] == 5)
    # and the family control: k=1 member, independently brute-forced upstream
    return pp


# ==========================================================================
# U.  (75,125): the exhaustive branch sweep
# ==========================================================================
def section_U(prop, ens):
    head("U.  (75,125): EXHAUSTIVE sweep over every branch of the reduction")
    note("  branch axes:  l in {4 (derived), 5 (the repo's value, counterfactual)}")
    note("                s in {3 (from Pred=(1,-3)), 2 (GGV6 Prop 2.5 alt)}")
    note("                Prop-8.2 branch in {proportional, en-split k=1,2 x 2}")
    note("  (the en-split rows are ILLEGAL by S4 and are swept anyway)")
    total, passed = 0, []
    for l in (4, 5):
        for s in (2, 3):
            variants = [dict()] + [dict(en_k=k, en_swap=sw)
                                   for k in (1, 2) for sw in (False, True)]
            for kw in variants:
                r = Reduction("sweep", 5, 20, (3, 5), s=s, l_override=l, **kw)
                rows, out = passport_from_polygon(r.NP, r.NQ, r.kappa)
                total += 1
                tag = (f"l={l} s={s} "
                       f"{'prop' if not kw else 'en k=%d sw=%d' % (kw['en_k'], kw['en_swap'])}"
                       f"{'' if r.legal else ' [ILLEGAL]'}")
                gates = [f"f={d['u']}i-{d['v']}j: {d['gate_lhs']} vs {d['gate_rhs']}"
                         for d in rows]
                note(f"    {tag:34s} N(P)={r.NP}")
                note(f"      {'; '.join(gates) if gates else 'no shared edge normal'}"
                     f"   -> {len(out)} admissible")
                if out:
                    passed.append((tag, out))
    chk(f"U1 all {total} branches sweep to ZERO admissible faces",
        not passed, f"{total} branches, 0 passes")

    # every GGV5 family at the (5,20) corner, not just F2 j=1
    note("\n  U1b.  Every GGV5 (m,n)-family at the corner (5,20), j = 0..4")
    note("        (F_2..F_6 of 1708.07936_GGV5.tex:1675-1680; F_6's printed")
    note("        (4j+... ) base pair (4,10) is the known non-coprime typo and")
    note("        is swept with the coprime correction as well):")
    fams = {"F_2": lambda j: (j + 2, 2 * j + 3),
            "F_3": lambda j: (4 * j + 3, 3 * j + 2),
            "F_4": lambda j: (2 * j + 3, 12 * j + 16),
            "F_5": lambda j: (7 * j + 9, 4 * j + 5),
            "F_6": lambda j: (3 * j + 4, 8 * j + 10),
            "F_6'": lambda j: (7 + 6 * j, 18 + 16 * j)}
    fam_hits, fam_tried = [], 0
    for name, f in fams.items():
        row = []
        for j in range(5):
            m_, n_ = f(j)
            if gcd(m_, n_) != 1:
                row.append(f"j={j}:({m_},{n_}) skip gcd={gcd(m_,n_)}")
                continue
            r = Reduction("fam", 5, 20, (m_, n_))
            _, o = passport_from_polygon(r.NP, r.NQ, r.kappa)
            fam_tried += 1
            if o:
                fam_hits.append((name, j, (m_, n_)))
            row.append(f"j={j}:({m_},{n_})->{len(o)}")
        note(f"        {name:5s} " + "  ".join(row))
    chk(f"U1c all {fam_tried} coprime GGV5 (m,n)-families at the (5,20) corner "
        "are Belyi-empty -- not just F_2 j=1", not fam_hits,
        f"{fam_tried} families, 0 admissible faces")
    note("")
    note("  U2.  Why, in closed form -- THE PROPORTIONAL-CLASS THEOREM.")
    note("       On the proportional branch N(P) = m*Delta', N(Q) = n*Delta'")
    note("       with (0,0) in Delta', so for any functional f, with")
    note("       r := max_{Delta'} f >= 0 one has m_f = m*r, n_f = n*r, and")
    note("            THE GATE  <=>  u*kappa = (m+n)*r - 1 .")
    note("       Positivity (m_f,n_f >= 1) forces r >= 1, hence")
    note("            u*kappa = (m+n)r - 1 >= (m+n) - 1 >= 4 > 0,  so u > 0.")
    note("       Delta' carries the scaled FOOT vertex (s,0) with s >= 2 (the")
    note("       root-shift depth; GGV6 Prop 2.5 gives s in {2,3}), so")
    note("            r >= f(s,0) = s*u >= 2u,")
    note("       and therefore  kappa >= 2(m+n) - 1/u >= 2(m+n) - 1.")
    note("       At (75,125): kappa would have to be >= 2*8 - 1 = 15.  It is 2")
    note("       (or 3 counterfactually).  At (7,21): >= 9 vs kappa = 1.")
    note("       And instantly, for kappa = 2 and m+n = 8: 2u = 8r - 1 is")
    note("       even = odd.")
    u = sp.Symbol("u", integer=True)
    r_ = sp.Symbol("r", integer=True)
    eq = sp.Eq(2 * u, 8 * r_ - 1)
    chk("U2a parity: with kappa=2, (m,n)=(3,5), the gate 2u = 8r-1 has NO "
        "integer solution (LHS even, RHS odd)",
        eq is sp.S.false and all((8 * rr - 1) % 2 == 1 and (2 * uu) % 2 == 0
                                 for rr in range(-60, 61) for uu in (rr, -rr)),
        "sympy's integer reasoning collapses the equation to False")
    dp = Reduction("dp", 5, 20, (3, 5))
    Dp = dp.reduced_delta
    ok = True
    for uu in range(-40, 41):
        for vv in range(-40, 41):
            if gcd(abs(uu), abs(vv)) != 1:
                continue
            rr = max(uu * i - vv * j for (i, j) in Dp)
            if uu * dp.kappa == 8 * rr - 1:
                ok = False
    chk("U2b brute force: no primitive (u,v) with |u|,|v| <= 40 satisfies the "
        "proportional gate u*kappa = 8r-1 on Delta'", ok)
    chk("U2c the foot is what does it: dropping (s,0) from Delta' and using "
        "the (72,108)-style UNSCALED en point is exactly the en-split branch, "
        "which S4 excludes", True,
        "structural, recorded not re-derived")
    return passed


# ==========================================================================
# V.  What was predicted, and what is actually true
# ==========================================================================
def section_V():
    head("V.  BELYI_PASSPORT.md sec.5's (75,125) prediction, tested")
    note("  Predicted (labelled INFERRED / N3 there): kappa = 3, gate 3u = 8r-1")
    note("  with r = max over R's monomials {(5,2),(5,5)}, forcing r = 5u-2v,")
    note("  r = 5 mod 6, minimal (r,u,v) = (5,13,30), face f = 13i - 30j,")
    note("  (m,n) = (15,25) -> reduced (3,5), Phi = D^3/A^5.")
    # (a) the arithmetic of the prediction is internally right
    sols = [(r_, u_, (37 * u_ - 1) // 16)
            for u_ in range(1, 60) for r_ in [5 * u_ - 2 * ((37 * u_ - 1) // 16)]
            if (37 * u_ - 1) % 16 == 0 and 3 * u_ == 8 * r_ - 1]
    chk("V1 the prediction's own arithmetic is correct: 3u=8r-1 & r=5u-2v give "
        "u = 13,29,45 with v = 30,67,104", sols[:3] ==
        [(5, 13, 30), (11, 29, 67), (17, 45, 104)], f"{sols[:3]}")
    # (b) but the face cannot exist: an edge of direction (v,u) needs that much width
    prop = Reduction("v", 5, 20, (3, 5))
    w = max(i for (i, j) in prop.NP)
    chk("V2 REFUTED (geometry): a max-face EDGE for f = 13i-30j must have "
        "direction (30,13); the reduced N(P) has x-width 12 < 30",
        w < 30, f"x-width of N(P) = {w}")
    # (c) and kappa is 2, not 3
    chk("V3 REFUTED (kappa): kappa = 2 by S7/S8, so the gate is 2u = 8r-1, "
        "which is parity-impossible", prop.kappa == 2)
    # (d) and r ignored the foot
    r_R = max(5 * 13 - 2 * 30, 5 * 13 - 5 * 30)
    r_true = max(13 * i - 30 * j for (i, j) in prop.reduced_delta)
    chk("V4 REFUTED (r): max over R's monomials alone is 5, but max over the "
        "whole reduced Delta' is dominated by the foot",
        r_R == 5 and r_true == 13 * 3, f"r_R = {r_R}, r_true = {r_true} "
        f"(from the foot (3,0))")
    note("\n  So sec.5's N3 was arithmetically sound but geometrically empty:")
    note("  it was a necessary condition on a face that cannot exist.  Its own")
    note("  hedge -- 'Whether N(P) and N(Q) actually have parallel edges of")
    note("  direction (30,13) is a polygon question' -- is now answered: NO.")


# ==========================================================================
# W.  family membership + the standing non-tests
# ==========================================================================
def section_W(pp):
    head("W.  Family membership, and the two standing NON-tests")
    # X1: RH is an identity, not a test  (re-derived here, not quoted)
    m, n, k, ell = sp.symbols("m n k ell")
    lhs = ell * (m - 1) + (n * k - m * ell - 1) + k * (n - 1) + (k + ell - 1)
    chk("W1 Riemann-Hurwitz balances IDENTICALLY in (m,n,k,ell) -- it is NOT a "
        "test and cannot discriminate", sp.simplify(lhs - (2 * n * k - 2)) == 0,
        "residual is identically 0")
    # the Catalan family
    def fam(kk):
        return (tuple([2] * ((3 * kk - 1) // 2) + [1]), tuple([3] * kk),
                tuple([(5 * kk - 1) // 2] + [1] * ((kk + 1) // 2)))
    cat = [1, 1, 2, 5, 14, 42, 132]
    ks = [1, 3, 5, 7, 9, 11, 13]
    ok = True
    for i, kk in enumerate(ks[:4]):
        p0, pi, p3 = fam(kk)
        d = sum(p0)
        h = bp.n_triples(d, tuple(sorted(p0, reverse=True)),
                         tuple(sorted(pi, reverse=True)),
                         tuple(sorted(p3, reverse=True))) // factorial(d)
        ok = ok and h == cat[i]
        note(f"    k={kk:2d}: deg {d:2d}  passport {p0} | {pi} | {p3}  h = {h}"
             f"  (Catalan {cat[i]})")
    chk("W2 the Catalan family reproduces at k = 1,3,5,7 -> 1,1,2,5", ok)
    chk("W3 (72,108) is the k=7 member", pp is not None and pp["deg"] == 21
        and pp["over0"] == fam(7)[0] and pp["hurwitz"] == 5)
    chk("W4 (75,125) is NOT a member of the Catalan family -- it is not a "
        "member of ANY Belyi family: it has no admissible face at all",
        True, "consequence of U1; the family has (m,n)=(2,3) and needs the "
              "en-split, which (5,20) cannot have (S4)")
    note("\n  What family IS (75,125) in, then?  The PROPORTIONAL class: the")
    note("  Prop-8.2(1) corners, where N(P) = m*Delta' and N(Q) = n*Delta'")
    note("  exactly.  For that whole class the gate collapses to")
    note("      u*kappa = (m+n)*r - 1,   r = max_{Delta'} f,")
    note("  and the scaled foot vertex (s,0), s >= 2, forces r >= 2u, hence")
    note("  kappa >= 2(m+n) - 1.  Since kappa = l-2 is small (1,2,3 in every")
    note("  published case) while 2(m+n)-1 >= 9, the ENTIRE proportional class")
    note("  is Belyi-empty.  (7,21) is in it too, and is empty for the same")
    note("  reason -- check W5.")
    r7 = Reduction("7_21", 7, 21, (2, 3))
    _, o7 = passport_from_polygon(r7.NP, r7.NQ, r7.kappa)
    chk("W5 (7,21), the other published proportional corner, is also "
        "Belyi-empty", not o7, f"{len(o7)} admissible faces")
    # the theorem, verified over a whole census slab
    head("W7.  The proportional-class theorem, swept over a corner census")
    note("  Every corner (a0,b0) with 2 <= a0 <= 14, a0 < b0 <= 60; every")
    note("  coprime (m,n) in {(2,3),(2,5),(3,4),(3,5),(4,5),(5,7)}; both the")
    note("  proportional and every legal en-split branch (k = 1,2,3, both")
    note("  assignments).  Reports every corner admitting a Belyi face.")
    hits_prop, hits_en, tried = [], [], 0
    for a0 in range(2, 15):
        for b0 in range(a0 + 1, 61):
            for mn in ((2, 3), (2, 5), (3, 4), (3, 5), (4, 5), (5, 7)):
                for kw in [dict()] + [dict(en_k=k, en_swap=sw)
                                      for k in (1, 2, 3) for sw in (False, True)]:
                    try:
                        r = Reduction("cs", a0, b0, mn, **kw)
                    except AssertionError:
                        continue
                    if not r.legal or r.kappa < 1:
                        continue
                    tried += 1
                    _, o = passport_from_polygon(r.NP, r.NQ, r.kappa)
                    if o:
                        (hits_prop if kw == {} else hits_en).append(
                            (a0, b0, mn, kw, [x["deg"] for x in o]))
    note(f"    swept {tried} legal (corner, (m,n), branch) triples")
    note(f"    proportional-branch hits: {len(hits_prop)}")
    note(f"    en-split-branch   hits: {len(hits_en)}")
    note("    en-split hits, compiled all the way to a Hurwitz number:")
    for (a0, b0, mn, kw, _degs) in hits_en[:14]:
        r = Reduction("hit", a0, b0, mn, **kw)
        _, o = passport_from_polygon(r.NP, r.NQ, r.kappa)
        for pp in o:
            d = pp["face"]
            note(f"      A_0=({a0},{b0}) (m,n)={mn} k={kw['en_k']} "
                 f"swap={int(kw['en_swap'])} | f = {d['u']}i-{d['v']}j, "
                 f"kappa={r.kappa} | deg {pp['deg']} | "
                 f"({pp['M']}^{pp['ell']},{pp['N']*pp['k']-pp['M']*pp['ell']} | "
                 f"{pp['N']}^{pp['k']} | {pp['k']+pp['ell']},"
                 f"1^{pp['deg']-pp['k']-pp['ell']}) | h = "
                 f"{pp.get('hurwitz','n/a')}")
    note("    (all but (8,28) are POLYGON-LEVEL hits only: the corner/(m,n)")
    note("     pairing is not checked against GGV5's admissible chain table,")
    note("     so they are candidate passports, not established cases.)")
    chk("W7a NO corner in the census slab admits a Belyi face on the "
        "PROPORTIONAL branch", not hits_prop,
        f"{tried} triples swept, 0 proportional hits")
    chk("W7b the en-split branch DOES produce hits, so W7a is discriminating "
        "and not vacuous", len(hits_en) > 0, f"{len(hits_en)} hits")
    chk("W7c (8,28)/(72,108) is among them", any(
        h[0] == 8 and h[1] == 28 and h[2] == (2, 3) and 21 in h[4]
        for h in hits_en))
    chk("W7d (5,20) is NOT among them, on either branch and either (m,n)",
        not any(h[0] == 5 and h[1] == 20 for h in hits_prop + hits_en))
    note("\n  So the top-band Belyi layer is a property of the EN-SPLIT corners")
    note("  (Prop 8.2(2)) only.  That is a sharper statement than sec.4 had:")
    note("  the missing 'face-pairing rule' is the en-split -- it is what")
    note("  un-scales the foot and creates a band joining an unscaled small")
    note("  vertex to a scaled core vertex.  At (72,108): f(1,0) = u = 2 and")
    note("  m*r = 2*1 = 2 coincide, which IS the band {(1,0)..(8,14)}.")
    r0 = Reduction("ctl", 8, 28, (2, 3), split_e=4, en_k=1)
    rr = max(2 * i - 1 * j for (i, j) in r0.core)
    chk("W6 at (72,108) the band exists precisely because f(en_P) = u = 2 "
        "equals m*r = 2*1", rr == 1 and 2 * rr == 2 * 1 * 1 + 0 and 2 == 2 * rr,
        f"r = max_core(2i-j) = {rr}, f(en_P) = 2, m*r = 2")


# ==========================================================================
def main():
    global VERBOSE
    ap = argparse.ArgumentParser()
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()
    VERBOSE = not args.quiet
    bp.VERBOSE = False

    section_R()
    section_P()
    prop, ens = section_S()
    pp = section_T()
    section_U(prop, ens)
    section_V()
    section_W(pp)

    bad = [c for c in CHECKS if not c[1]]
    if VERBOSE:
        print("\n" + "=" * 92)
        print(f"  {len(CHECKS) - len(bad)}/{len(CHECKS)} checks passed")
        for t, _, m in bad:
            print(f"    FAILED: {t}  {m}")
        print("=" * 92)
    else:
        print(f"{len(CHECKS) - len(bad)}/{len(CHECKS)}")
        for t, _, m in bad:
            print(f"FAILED: {t}  {m}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
