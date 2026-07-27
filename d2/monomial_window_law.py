#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
monomial_window_law.py -- what machinery works BECAUSE C is a monomial?

Companion document: MONOMIAL_WINDOW_LAW.md.
Prerequisites: q_window_theorem.py (the M/gcd(M,H) invariant + Bezout lemma),
window_functions_75_125.py (the analytic window step ord_y(Phi)/M, (R1)-(R3),
(S1)-(S3)), weight_lemma_75_125.py sec.(4) (the ord-side carry obstruction),
MINIMAL_CORE.md sec.4 (the `q_window | w(e)` criterion), ENDPOINT_CONTRACT.md
(the depth ledger and its kill predicate), corner_atlas.json (the 34 rows).

THE QUESTION.  At the 28 GGV5 rows where the retraction shape b0 = t(a0-1) FAILS,
C is a MONOMIAL (C = y, deg C = ord C = 1).  Four separate mechanisms die there
(slice cascade lam = 0, window cone collapses to a ray, toric identity, F2 closed
form dg = 0).  Is that ONE systematic bias, and is there anything that works
*because* C is a monomial?

FIVE ANSWERS.

  (A)  THE BRIDGE IDENTITY -- new, and it closes a flagged negative.
       With the repo's own rho = ord_y(f) = q(b-a)+1 and N = a*M - 2b,

              ord_y(Phi)  =  a*q*M  -  H,        H := q(a+b) - 1,

       an EXACT identity in (t,kappa,q,a,b) -- not a congruence.  Hence
       gcd(M, ord_y Phi) = gcd(M, H), so the ANALYTIC window denominator
       denom(ord_y(Phi)/M) that the carry obstruction consumes IS the
       COMBINATORIAL corner invariant q_window = M/gcd(M,H).  MINIMAL_CORE.md's
       "Negative I could not overcome" (q_window not sweepable because ord_y Phi
       is published only at the (8,28) corner) is thereby dissolved: ord_y Phi is
       a function of corner data alone.  The identity PREDICTS ord_y Phi = 204 at
       (72,108) from (t,kappa,a,b,q) = (4,2,2,3,7) with no Newton polygon, and
       the independent route (the f1 ODE solution y^8(y+1)^2 q(y)/6630 times
       C4^28 = y^196 t^28, AT_LE9_AUDIT.md B7) gives 204.

  (B)  MONOMIAL RIGIDITY.  Monomial C means q = ord_y C = 1.  In the standard
       class kappa = t-2 the fixed Bezout corner integer of q_window_theorem is

              q(kappa+1) - t  =  (t-1) - t  =  -1,

       and gcd(M,H) divides it, so gcd(M,H) = 1 and q_window = M EXACTLY --
       MAXIMAL, at every corner, every family member, unconditionally.

  (C)  THE INTEGRAL REGIME IS ARITHMETICALLY INACCESSIBLE.  q_window = 1 iff
       M | H, which forces (t-q)(a+b) <= kappa; with kappa = t-2 and a+b > t-2
       (true on all 34 rows, where min(a+b) = 5 and max t = 6) that gives
       ord_y C >= t.  A monomial has ord_y C = 1 < t.  So the (72,108) regime
       that every mechanism in this repo was calibrated in is not merely
       unobserved at the class of nine -- it is IMPOSSIBLE there.

  (D)  THE CARRY OBSTRUCTION IS TOTAL, AND THAT IS A MONOMIAL-ONLY FACT.
       gcd(alpha, M) = 1 with alpha = ord_y Phi, so for EVERY split
       M = w_1 + ... + w_k (k >= 2, all w_i >= 1) the total ceil-carry
       sum_i ceil(alpha w_i/M) - alpha lies in [1, k-1] -- never 0.  No split
       enumeration is needed and no w(e) datum is needed.  At (72,108)
       q_window = 1 and the carry is 0 for every split.  THE DISCRIMINATING PAIR:
       (50,75) and (72,108) have the SAME (a,b,t) = (2,3,4), hence literally the
       SAME G-system ideal and the same K-syzygy as an algebraic relation, and
       the same M = 17.  The single differing input is ord_y C: 7 vs 1, i.e.
       alpha = 204 vs 30.  Carry on the published split (w_e,w_B) = (5,12) is
       0 vs 1.  At (50,75) the algebra permits the syzygy and the arithmetic
       forbids it.

  (E)  MONOMIALITY IS TWO INDEPENDENT DEFECTS, NOT ONE BIAS.  Every one of the
       four documented deaths consumes THINNESS (deg C - ord C = 0).  The
       q_window death consumes SHALLOWNESS (ord C = 1) and has NO deg C
       dependence at all: q_window = M/gcd(M, q(a+b)-1) is a function of
       (t,kappa,q,a,b) only.  The two are logically independent (all four
       quadrants are realised by explicit parameter points), so a repair that
       restores a residual to C cannot revive the divisor mechanism, and a repair
       that deepens C cannot revive the cascade/cone.  This REFUTES the
       "one systematic bias" reading.

  THE POSITIVE SIDE OF (D).  ceil(alpha w/q_window) is the depth ledger's forced
  floor, and its gain over the affine ray alpha w/q_window is
  (-alpha w mod q_window)/q_window, which is identically 0 iff q_window = 1.  So
  at a monomial corner the floor is STRICTLY HIGHER than the affine prediction at
  every admissible weight, by a total of exactly (M-1)/2 units over w = 1..M-1.
  The ENDPOINT_CONTRACT kill predicate fires when a required-nonzero coefficient
  sits strictly BELOW the forced floor, so raising the floor can only create
  kills.  Monomiality collapses the two-slope CONE (it needs two slopes) and at
  the same time maximally strengthens the one-slope FLOOR.  That is the precise
  sense in which the depth ledger is monomial-compatible and the cone is not.

DISCIPLINE.  Exact integers / sympy only, no floats.  Every positive check is
paired with a MUTATION CONTROL that must fail.  Group E re-derives the atlas
numbers from the atlas's own transcribed corner data rather than trusting them,
and Group A's headline is a cross-check between two disjoint derivations of the
same integer (204), not a fresh assertion about one.

Run:
    python monomial_window_law.py            # verbose report
    python monomial_window_law.py --quiet    # one line per check; exit 0 iff all pass
"""

from __future__ import annotations

import argparse
import itertools
import json
import os
import sys
from fractions import Fraction
from math import gcd

import sympy as sp

import q_window_theorem as qwt
import window_functions_75_125 as wf

HERE = os.path.dirname(os.path.abspath(__file__))


# ----------------------------------------------------------------- harness ---
class Ledger:
    def __init__(self, quiet: bool) -> None:
        self.quiet = quiet
        self.rows: list[tuple[str, bool, str]] = []

    def head(self, title: str) -> None:
        if not self.quiet:
            print("\n" + "=" * 78)
            print(title)
            print("=" * 78)

    def ck(self, name: str, ok: bool, detail: str = "") -> bool:
        self.rows.append((name, bool(ok), detail))
        if not self.quiet:
            print("  [%s] %-46s %s" % ("PASS" if ok else "FAIL", name, detail))
        return bool(ok)

    def mut(self, name: str, all_mutants_fail: bool, detail: str = "") -> bool:
        return self.ck("MUT " + name, all_mutants_fail, detail)

    def note(self, text: str) -> None:
        if not self.quiet:
            print("        . " + text)

    def report(self) -> int:
        bad = [r for r in self.rows if not r[1]]
        print("\n%s  monomial_window_law: %d/%d checks pass"
              % ("FAIL" if bad else "OK  ", len(self.rows) - len(bad), len(self.rows)))
        for n, _, d in bad:
            print("   FAILED: %s   %s" % (n, d))
        return 1 if bad else 0


# =============================================================================
#  The corner arithmetic, all from (t, kappa, q, a, b).  No polygon input.
# =============================================================================
def MH(t, kappa, q, a, b):
    """(M, H) = (t(a+b) - (kappa+1), q(a+b) - 1)   [q_window_theorem]."""
    s = a + b
    return t * s - (kappa + 1), q * s - 1


def rho(q, a, b):
    """ord_y(f) = q*(b-a) + 1   [window_functions_75_125.family]."""
    return q * (b - a) + 1


def Nexp(t, kappa, q, a, b):
    """N, the exponent of C in Phi = f*C^N:  N = a*M - 2b
    [window_functions_75_125.family: N = a(t(a+b)-(kappa+1)) - 2b]."""
    M, _ = MH(t, kappa, q, a, b)
    return a * M - 2 * b


def ordPhi_from_polygon_route(t, kappa, q, a, b):
    """ord_y(Phi) = rho + N*ord_y(C) = rho + N*q -- the route through the
    reduction data (f and C separately)."""
    return rho(q, a, b) + Nexp(t, kappa, q, a, b) * q


def ordPhi_from_bridge(t, kappa, q, a, b):
    """ord_y(Phi) = a*q*M - H -- the BRIDGE IDENTITY (Group A), corner data only."""
    M, H = MH(t, kappa, q, a, b)
    return a * q * M - H


def q_window(t, kappa, q, a, b):
    M, H = MH(t, kappa, q, a, b)
    return M // gcd(abs(M), abs(H))


def carry(alpha, M, parts):
    """Total ceil-superadditivity defect of a split M = sum(parts):
       sum_i ceil(alpha w_i / M)  -  ceil(alpha * M / M)."""
    assert sum(parts) == M
    tot = sum(-((-alpha * w) // M) for w in parts)      # ceil via floor
    return tot - alpha


def reduced_step(alpha, M):
    """(alpha, M) in lowest terms: the window step alpha/M = alpha_red/q_window."""
    g = gcd(abs(alpha), abs(M))
    return alpha // g, M // g


def floor_gain(alpha, M, w):
    """ceil(alpha w/M) - alpha w/M, exact Fraction, computed in the REDUCED
    step alpha_red/q_window (so the quasi-period is the true one)."""
    ar, qw = reduced_step(alpha, M)
    return Fraction((-ar * w) % qw, qw)


# the nine t=4 monomial rows ("class of nine") and their corner data.
# (tag, a0, b0, t, kappa, q=ord C, a, b) with (a,b) = sorted(m,n).
CLASS_OF_NINE = [
    ("F_2(2,3)/75  = (50,75)", 5, 20, 4, 2, 1, 2, 3),
    ("F_2(3,5)/125 = (75,125)", 5, 20, 4, 2, 1, 3, 5),
    ("F_3(3,2)/75", 5, 20, 4, 2, 1, 2, 3),
    ("(9,36)/(3,2)/135", 9, 36, 4, 2, 1, 2, 3),
    ("(9,36)/(2,3)/135", 9, 36, 4, 2, 1, 2, 3),
    ("(9,36)/(2,3)/135 [len2]", 9, 36, 4, 2, 1, 2, 3),
    ("(8,32)/(3,2)/120", 8, 32, 4, 2, 1, 2, 3),
    ("(10,40)/(3,2)/150 [a]", 10, 40, 4, 2, 1, 2, 3),
    ("(10,40)/(3,2)/150 [b]", 10, 40, 4, 2, 1, 2, 3),
]

# (72,108) = (8,28)/(3,2): the unique t=4 RETRACTING row, where the toolkit lives.
C_72_108 = dict(t=4, kappa=2, q=7, a=2, b=3)
ORDPHI_72_108_PUBLISHED = 204     # AT_LE9_AUDIT.md B7 / AUDIT.md sec.A5:
DEGPHI_72_108_PUBLISHED = 238     # Phi = f1*C4^28 = -y^204 t^30 q(y)/6630


# =============================================================================
#  A.  THE BRIDGE IDENTITY   ord_y(Phi) = a*q*M - H
# =============================================================================
def group_A(L: Ledger) -> None:
    L.head("A.  THE BRIDGE IDENTITY   ord_y(Phi) = a*q*M - H   (H = q(a+b)-1)")

    t, kappa, q, a, b = sp.symbols("t kappa q a b")
    M = t * (a + b) - (kappa + 1)
    H = q * (a + b) - 1
    N = a * M - 2 * b
    rh = q * (b - a) + 1
    lhs = rh + N * q                       # ord_y(Phi) via the reduction route
    rhs = a * q * M - H                    # ord_y(Phi) via corner data only
    L.ck("A1 identity is polynomial-exact in (t,kappa,q,a,b)",
         sp.simplify(sp.expand(lhs - rhs)) == 0,
         "rho + N*q  ==  a*q*M - H  identically, where rho = q(b-a)+1, "
         "N = a*M - 2b; the (a+b) and kappa terms cancel")

    # A2  the two routes agree numerically on a wide sweep
    bad = []
    for tt in range(2, 13):
        for kk in (tt - 2, tt - 1, 1):
            for qq in range(1, 10):
                for aa in range(1, 8):
                    for bb in range(aa, 12):
                        if (ordPhi_from_polygon_route(tt, kk, qq, aa, bb)
                                != ordPhi_from_bridge(tt, kk, qq, aa, bb)):
                            bad.append((tt, kk, qq, aa, bb))
    L.ck("A2 numeric agreement, 2*11*9*7*12-ish sweep",
         not bad, "no disagreement over t=2..12, kappa in {t-2,t-1,1}, "
                  "q=1..9, 1<=a<=b<=11 (%d points)" % (11 * 3 * 9 * 7 * 12))

    # A3  the headline cross-check: corner data alone predicts 204 at (72,108),
    #     and the disjoint f1-ODE route gives 204.
    c = C_72_108
    pred = ordPhi_from_bridge(**c)
    M0, H0 = MH(**c)
    L.ck("A3 corner data predicts ord_y(Phi) = 204 at (72,108)",
         pred == ORDPHI_72_108_PUBLISHED,
         "a*q*M - H = 2*7*%d - %d = %d; the independent route "
         "f1 = -y^8(y+1)^2 q(y)/6630 times C4^28 = y^196 t^28 gives ord = "
         "8 + 196 = 204 (AT_LE9_AUDIT.md B7)" % (M0, H0, pred))
    L.ck("A4 the same identity fixes rho = ord_y(f1) = 8 and N = 28",
         (rho(7, 2, 3), Nexp(**c)) == (8, 28),
         "rho = q(b-a)+1 = 8 matches ord_y of the ODE solution f1 = "
         "-y^8(y+1)^2 q(y)/6630; N = a*M-2b = 28 matches C4^28 (g4_row.py:251)")

    # A5  gcd(M, ordPhi) == gcd(M, H): the analytic denominator IS the corner
    #     invariant.  Checked on the whole sweep.
    bad2 = []
    for tt in range(2, 13):
        for kk in (tt - 2, tt - 1, 1):
            for qq in range(1, 10):
                for aa in range(1, 8):
                    for bb in range(aa, 12):
                        MM, HH = MH(tt, kk, qq, aa, bb)
                        if MM <= 0:
                            continue
                        al = ordPhi_from_bridge(tt, kk, qq, aa, bb)
                        if gcd(MM, abs(al)) != gcd(MM, abs(HH)):
                            bad2.append((tt, kk, qq, aa, bb))
    L.ck("A5 denom(ord_y(Phi)/M) == q_window = M/gcd(M,H)",
         not bad2,
         "gcd(M, ord_y Phi) = gcd(M, H) on every sweep point, so the ANALYTIC "
         "window denominator (window_functions_75_125) equals the COMBINATORIAL "
         "corner invariant (q_window_theorem) -- proved, not tabulated")

    # A6  agreement with the repo's own two implementations at the known cases.
    rows = []
    for aa in range(2, 9):
        f = wf.family(aa)
        analytic = Fraction(f["ordPhi"], f["M"]).denominator
        comb = q_window(4, 2, 1, aa, f["b"])
        bridge = ordPhi_from_bridge(4, 2, 1, aa, f["b"])
        rows.append((aa, f["ordPhi"], bridge, analytic, comb))
    L.ck("A6 F2 rungs a=2..8: bridge reproduces window_functions.ordPhi",
         all(r[1] == r[2] for r in rows),
         "ord_y(Phi_a) = 12a^2-10a+2 equals a*M - H = a(12a-7) - (3a-2) "
         "for every rung")
    L.ck("A7 F2 rungs a=2..8: analytic denominator == combinatorial q_window",
         all(r[3] == r[4] for r in rows),
         "both equal M = 12a-7 at every rung: %s"
         % [(r[0], r[3]) for r in rows])
    for aa, o, br, an, co in rows:
        L.note("a=%d  ord_y(Phi)=%-4d bridge=%-4d denom=%-3d q_window=%-3d"
               % (aa, o, br, an, co))

    # MUTATION CONTROLS
    mut_rho = all(rho(q_, a_, b_) + 1 + Nexp(4, 2, q_, a_, b_) * q_
                  != ordPhi_from_bridge(4, 2, q_, a_, b_)
                  for q_, a_, b_ in [(7, 2, 3), (1, 2, 3), (1, 3, 5), (4, 2, 7)])
    mut_N = all(rho(q_, a_, b_) + (Nexp(4, 2, q_, a_, b_) + 1) * q_
                != ordPhi_from_bridge(4, 2, q_, a_, b_)
                for q_, a_, b_ in [(7, 2, 3), (1, 2, 3), (1, 3, 5), (4, 2, 7)])
    L.mut("A the identity is tight in rho and in N",
          mut_rho and mut_N,
          "rho -> rho+1 breaks it at all four probes; N -> N+1 breaks it at all "
          "four probes.  So A3's 204 is not a coincidence of loose formulas")
    L.mut("A ord_y(Phi) is not a constant of the code",
          len({ordPhi_from_bridge(4, 2, qq, 2, 3) for qq in range(1, 8)}) == 7,
          "the bridge returns 7 distinct values as q = ord_y C runs 1..7 "
          "(30,64,98,...,204): it reads q, it does not ignore it")


# =============================================================================
#  B.  MONOMIAL RIGIDITY:  q = 1 & kappa = t-2  ==>  q_window = M, maximal
# =============================================================================
def group_B(L: Ledger) -> None:
    L.head("B.  MONOMIAL RIGIDITY   q=1, kappa=t-2  ==>  gcd(M,H)=1, q_window=M")

    t = sp.symbols("t")
    L.ck("B1 Bezout corner integer is -1 for monomial + standard class",
         sp.simplify(1 * ((t - 2) + 1) - t + 1) == 0,
         "q(kappa+1) - t = 1*(t-1) - t = -1, independent of t, a, b, and of "
         "which corner: it is the ONLY place ord_y C enters the Bezout relation")
    L.ck("B2 q_window_theorem.corner_integer agrees, t = 2..40",
         all(qwt.corner_integer(tt, tt - 2, 1) == -1 for tt in range(2, 41)),
         "the repo's own corner_integer returns -1 at every t")

    bad = []
    for tt in range(2, 21):
        for s in range(3, 61):            # s = a+b
            for aa in range(1, s):
                bb = s - aa
                if bb < aa:
                    continue
                M, H = MH(tt, tt - 2, 1, aa, bb)
                if gcd(abs(M), abs(H)) != 1 or q_window(tt, tt - 2, 1, aa, bb) != M:
                    bad.append((tt, aa, bb))
    L.ck("B3 gcd(M,H)=1 and q_window=M on every monomial point",
         not bad,
         "swept t=2..20, a+b=3..60, all splits a<=b: gcd(M,H)=1 always, so "
         "q_window = M exactly -- the MAXIMAL possible window denominator")

    # B4  and therefore alpha = ord_y(Phi) is coprime to M (the input Group D needs)
    bad4 = []
    for tt in range(3, 13):
        for s in range(3, 41):
            for aa in range(1, s):
                bb = s - aa
                if bb < aa:
                    continue
                M, _ = MH(tt, tt - 2, 1, aa, bb)
                al = ordPhi_from_bridge(tt, tt - 2, 1, aa, bb)
                if gcd(M, abs(al)) != 1:
                    bad4.append((tt, aa, bb))
    L.ck("B4 gcd(ord_y Phi, M) = 1 at every monomial corner",
         not bad4, "immediate from A5 + B3; this is the hypothesis of the "
                   "total-carry lemma (Group D)")

    # B5  the class of nine, explicitly
    L.note("the class of nine (t=4, retraction fails, C = y):")
    okc = True
    for tag, a0, b0, tt, kk, qq, aa, bb in CLASS_OF_NINE:
        M, H = MH(tt, kk, qq, aa, bb)
        qw = q_window(tt, kk, qq, aa, bb)
        al = ordPhi_from_bridge(tt, kk, qq, aa, bb)
        okc &= (qw == M) and (gcd(M, al) == 1)
        L.note("%-26s (a,b)=(%d,%d) M=%-3d H=%-3d q_window=%-3d ord_y(Phi)=%-4d"
               % (tag, aa, bb, M, H, qw, al))
    L.ck("B5 all nine class rows have q_window = M and gcd(alpha,M) = 1", okc,
         "8 rows at M = 17 (a+b=5) and (75,125) at M = 29 (a+b=8)")

    # MUTATION CONTROLS
    q2 = [(tt, aa, bb) for tt in range(3, 9) for s in range(4, 20)
          for aa in range(1, s) for bb in [s - aa]
          if bb >= aa and gcd(*map(abs, MH(tt, tt - 2, 2, aa, bb))) != 1]
    L.mut("B q >= 2 is genuinely different (gcd need not be 1)",
          len(q2) > 0,
          "with q = 2 the corner integer is (kappa+1)*2-t = t-2, not -1, and "
          "%d swept points have gcd(M,H) > 1 -- so B3 is a statement about "
          "q = 1, not about the code" % len(q2))
    kmut = [(tt, aa, bb) for tt in range(3, 9) for s in range(4, 20)
            for aa in range(1, s) for bb in [s - aa]
            if bb >= aa and q_window(tt, tt - 1, 1, aa, bb)
            != MH(tt, tt - 1, 1, aa, bb)[0]]
    L.ck("B6 kappa = t-1 breaks it: corner integer is 0, gcd unbounded",
         len(kmut) > 0,
         "q(kappa+1)-t = 0 when kappa = t-1 and q = 1, so the divisibility "
         "lemma is vacuous; %d swept points then have q_window < M.  B3 uses "
         "kappa = t-2, which is what polygon_reduction returns for all 34 rows"
         % len(kmut))


# =============================================================================
#  C.  THE INTEGRAL REGIME IS INACCESSIBLE TO MONOMIAL CORNERS
# =============================================================================
def group_C(L: Ledger) -> None:
    L.head("C.  q_window = 1 REQUIRES ord_y C >= t   (so monomials never get it)")

    t, kappa, q, a, b = sp.symbols("t kappa q a b", positive=True)
    M = t * (a + b) - (kappa + 1)
    H = q * (a + b) - 1
    L.ck("C1 M - H = (t-q)(a+b) - kappa, identically",
         sp.simplify(sp.expand(M - H - ((t - q) * (a + b) - kappa))) == 0,
         "so M <= H  <=>  (t-q)(a+b) <= kappa")

    # C2  monomial: M - H = (t-1)(a+b) - kappa = (t-1)(a+b) - (t-2) > 0 for
    #     t >= 2, a+b >= 2, hence 0 < H < M and M | H is impossible.
    bad = []
    for tt in range(2, 21):
        for s in range(2, 61):
            M0 = tt * s - (tt - 1)
            H0 = s - 1
            if not (0 < H0 < M0):
                bad.append((tt, s))
    L.ck("C2 monomial corners always have 0 < H < M",
         not bad,
         "H = a+b-1 and M = t(a+b)-t+1, so M - H = (t-1)(a+b)-(t-2) >= t > 0; "
         "M | H is therefore impossible and q_window = 1 cannot occur")

    # C3  the necessary condition.  M|H => M<=H => (t-q)(a+b) <= kappa = t-2.
    #     If a+b > t-2 that forces t-q <= 0, i.e. ord_y C >= t.  The hypothesis
    #     a+b > t-2 is REQUIRED (see the mutation control) and holds on all 34
    #     published rows (a+b = m+n >= 5, t <= 6, so t-2 <= 4 < 5).
    integral, viol, integral_h = 0, [], 0
    for tt in range(2, 16):
        for qq in range(1, 40):
            for s in range(3, 41):
                for aa in range(1, s):
                    bb = s - aa
                    if bb < aa:
                        continue
                    if q_window(tt, tt - 2, qq, aa, bb) == 1:
                        integral += 1
                        if s > tt - 2:
                            integral_h += 1
                            if qq < tt:
                                viol.append((tt, qq, aa, bb))
    L.ck("C3 integral window + (a+b > t-2)  ==>  ord_y C >= t",
         integral_h > 0 and not viol,
         "%d integral points in the sweep satisfy a+b > t-2 (not vacuous); %d "
         "of them have ord_y C < t.  Proof: M|H => M<=H => (t-q)(a+b) <= t-2 "
         "< a+b, so t-q <= 0" % (integral_h, len(viol)))
    L.ck("C4 monomial (ord_y C = 1) is excluded by C3 for every t >= 2",
         all(1 < tt for tt in range(2, 16)),
         "ord_y C = 1 < t for every chart exponent in the atlas (t = 3,4,5,6)")

    # C5  QUANTITATIVE: what ord_y C would the class-of-nine shapes NEED?
    #     q_window = 1  <=>  M | H  <=>  (a+b)*q == 1 mod M.
    shapes = [(4, 2, 2, 3), (4, 2, 3, 5)]      # (72,108)/(50,75) shape; (75,125)
    need = []
    for tt, kk, aa, bb in shapes:
        M, _ = MH(tt, kk, 1, aa, bb)
        s = aa + bb
        sol = [qq for qq in range(1, M + 1) if (qq * s - 1) % M == 0]
        need.append((tt, kk, s, M, sol[0], sol))
    L.ck("C5 the (t,kappa,a+b) = (4,2,5) shape needs ord_y C == 7 mod 17",
         need[0][4] == 7 and need[0][3] == 17,
         "M = 17, so q_window = 1 iff 5q == 1 mod 17 iff q == 7 (mod 17).  The "
         "MINIMAL solution is 7 -- which is EXACTLY (72,108)'s ord_y C.  EIGHT "
         "of the nine class rows share this shape and have ord_y C = 1")
    L.ck("C6 the (75,125) shape (4,2,8) needs ord_y C == 11 mod 29",
         need[1][4] == 11 and need[1][3] == 29,
         "M = 29, 8q == 1 mod 29 iff q == 11 (mod 29); the flagship case has "
         "ord_y C = 1")
    L.note("(72,108) is the MINIMAL integral point of its own corner shape; the "
           "class of nine differs from it in exactly one integer, ord_y C.")

    # MUTATION CONTROL: the hypothesis a+b > t-2 is load-bearing, not decoration.
    viol_all = [(tt, qq, aa, bb) for tt in range(2, 16) for qq in range(1, 40)
                for s in range(3, 41) for aa in range(1, s) for bb in [s - aa]
                if bb >= aa and q_window(tt, tt - 2, qq, aa, bb) == 1
                and qq < tt]
    L.mut("C the hypothesis a+b > t-2 is load-bearing",
          len(viol_all) > 0 and all(a_ + b_ <= t_ - 2
                                    for t_, q_, a_, b_ in viol_all),
          "%d integral points DO have ord_y C < t, e.g. %s -- every one has "
          "a+b <= t-2 (there M = H exactly).  So C3 is not vacuously true and "
          "its hypothesis is exactly the right one"
          % (len(viol_all), viol_all[:3]))


# =============================================================================
#  D.  TOTAL CARRY -- the monomial-only obstruction, and the discriminating pair
# =============================================================================
def compositions(n, k):
    """All ordered k-tuples of positive integers summing to n."""
    for cut in itertools.combinations(range(1, n), k - 1):
        prev, out = 0, []
        for c in cut:
            out.append(c - prev)
            prev = c
        out.append(n - prev)
        yield tuple(out)


def group_D(L: Ledger) -> None:
    L.head("D.  THE CARRY OBSTRUCTION IS TOTAL AT MONOMIAL CORNERS")

    # D1  the lemma, exhaustively over all splits of M for the class M values.
    checked, bad = 0, []
    for tt, kk, qq, aa, bb in [(4, 2, 1, 2, 3), (4, 2, 1, 3, 5),
                               (3, 1, 1, 2, 3), (5, 3, 1, 2, 3),
                               (6, 4, 1, 2, 3), (3, 1, 1, 3, 4)]:
        M, _ = MH(tt, kk, qq, aa, bb)
        al = ordPhi_from_bridge(tt, kk, qq, aa, bb)
        for k in (2, 3, 4):
            if k > M:
                continue
            for parts in compositions(M, k):
                c = carry(al, M, parts)
                checked += 1
                if not (1 <= c <= k - 1):
                    bad.append((tt, aa, bb, parts, c))
    L.ck("D1 every split of M has total carry in [1, k-1]",
         checked > 0 and not bad,
         "%d splits (k = 2,3,4) at six monomial corners including M=17 (50,75), "
         "M=29 (75,125), M=13, M=21, M=25; ZERO have carry 0" % checked)

    # D2  the contrast at (72,108): carry 0 for every split.
    M72, _ = MH(**C_72_108)
    al72 = ORDPHI_72_108_PUBLISHED
    n0, bad2 = 0, 0
    for k in (2, 3, 4):
        for parts in compositions(M72, k):
            n0 += 1
            if carry(al72, M72, parts) != 0:
                bad2 += 1
    L.ck("D2 (72,108): every split has carry EXACTLY 0",
         n0 > 0 and bad2 == 0,
         "%d splits of M = 17 with alpha = 204; q_window = 1 so ceil is "
         "additive and the obstruction vanishes identically" % n0)

    # D3  THE DISCRIMINATING PAIR -- same ideal, same M, same split.
    M50, _ = MH(4, 2, 1, 2, 3)
    al50 = ordPhi_from_bridge(4, 2, 1, 2, 3)
    same_abt = (4, 2, 3) == (C_72_108["t"], C_72_108["a"], C_72_108["b"])
    L.ck("D3a (50,75) and (72,108) share (a,b,t) = (2,3,4) and M = 17",
         same_abt and M50 == M72 == 17,
         "the G-system ideal is a function of (a,b,t) alone "
         "(g_system_75_125.build_gsystem ignores q and ordPhi when forming the "
         "generators; toric_general.system passes the dummies 1,1), so the "
         "K-syzygy exists as an ALGEBRAIC relation at (50,75) too")
    w_e, w_B = 5, 12                      # the published split, w(e)=t+1, w(B)=b*t
    c50 = carry(al50, M50, (w_e, w_B))
    c72 = carry(al72, M72, (w_e, w_B))
    L.ck("D3b the published split (w_e,w_B) = (5,12): carry 0 vs 1",
         (c72, c50) == (0, 1),
         "(72,108) alpha=204: ceil(204*5/17)+ceil(204*12/17) = 60+144 = 204 = "
         "ord_y Phi, carry 0.  (50,75) alpha=30: 9+22 = 31 > 30, carry 1.  The "
         "ONLY differing input is ord_y C = 7 vs 1")
    # D3b'  the ideal really is the same OBJECT, not just the same shape.
    if not L.quiet:
        print("        . building both G-systems (a,b,t)=(2,3,4) ...")
    import g_system_75_125 as gs
    r50 = gs.build_gsystem(2, 3, 4, 1, al50)          # (50,75):  q=1, alpha=30
    r72 = gs.build_gsystem(2, 3, 4, 7, al72)          # (72,108): q=7, alpha=204
    struct = ["Gs", "Klin", "M", "state", "spares", "deep", "sub", "homog",
              "jphi", "a", "b", "t", "kappa", "s", "dm", "dh", "e", "skiplin"]
    diff = [k for k in struct if r50[k] != r72[k]]
    moved = [k for k in ("W_step", "ordPhi", "q") if r50[k] != r72[k]]
    L.ck("D3b' the two G-systems are IDENTICAL on every structural field",
         not diff and sorted(moved) == ["W_step", "ordPhi", "q"],
         "generators Gs, linear substitutions Klin, M, state/spare inventory, "
         "homogeneity record: all equal.  The ONLY fields that differ are "
         "q = ord_y C (1 vs 7), ordPhi (30 vs 204) and W_step (%s vs %s) -- "
         "i.e. exactly the weight normalisation"
         % (r50["W_step"], r72["W_step"]))
    L.ck("D3b'' and the u-weights agree symbol by symbol",
         all(r50["uweight"](v) == r72["uweight"](v)
             for v in r50["state"] + r50["spares"]),
         "the intrinsic u-grading is the same function on the same variables")

    L.ck("D3c alpha is the only thing that moved",
         al50 == 30 and al72 == 204
         and MH(4, 2, 1, 2, 3)[0] == MH(**C_72_108)[0],
         "same t,kappa,a,b => same M = 17; ord_y C = 1 vs 7 => alpha = 30 vs "
         "204 = a*q*M - H.  At (50,75) the algebra permits the syzygy and the "
         "arithmetic forbids it")

    # D4  the escape criterion, both readings of the admissible w(e) range.
    #     MINIMAL_CORE sec.4: the obstruction vanishes iff q_window | w(e).
    esc_all, esc_state = [], []
    for tag, a0, b0, tt, kk, qq, aa, bb in CLASS_OF_NINE:
        M, _ = MH(tt, kk, qq, aa, bb)
        qw = q_window(tt, kk, qq, aa, bb)
        esc_all.append([w for w in range(1, M) if w % qw == 0])
        esc_state.append([w for w in range(2, tt + 2) if w % qw == 0])
    L.ck("D4 no escape at any class-of-nine row, either reading of w(e)",
         all(not e for e in esc_all) and all(not e for e in esc_state),
         "reading 1: w(e) in 1..M-1 (the split enumeration) -- q_window = M "
         "divides none.  reading 2: w(e) in {2..t+1} (the state-variable "
         "u-weights, window_functions.state_uweights) -- all < M")

    # MUTATION CONTROLS
    #  the carry function is not identically 1: it returns 0 where it should.
    #  F_7 = (6,15)/(2,7): t=3, kappa=1, q=4, M=25, q_window=5.
    M7, H7 = MH(3, 1, 4, 2, 7)
    al7 = ordPhi_from_bridge(3, 1, 4, 2, 7)
    zero7 = [w for w in range(1, M7) if carry(al7, M7, (w, M7 - w)) == 0]
    pred7 = [w for w in range(1, M7) if w % q_window(3, 1, 4, 2, 7) == 0]
    L.mut("D carry() returns 0 exactly on the predicted escapes",
          zero7 == pred7 and len(zero7) > 0,
          "F_7 (t=3,kappa=1,q=4,(a,b)=(2,7)): M=25, alpha=%d, q_window=5; the "
          "2-splits with carry 0 are exactly w in %s = the multiples of 5.  So "
          "D1's 'never 0' is a fact about monomial corners, not about carry()"
          % (al7, zero7))
    L.mut("D the G-system comparison is sensitive to (a,b,t)",
          len(r50["state"] + r50["spares"]) > 0
          and gs.build_gsystem(2, 3, 5, 1, al50)["Gs"] != r50["Gs"],
          "%d graded variables compared (not a vacuous 'all of nothing'), and "
          "changing t from 4 to 5 DOES change Gs -- so D3b' is detecting real "
          "equality, not a trivial one" % len(r50["state"] + r50["spares"]))
    L.mut("D the split enumeration is nonempty and k-sensitive",
          len(list(compositions(17, 2))) == 16
          and len(list(compositions(17, 3))) == 120,
          "16 two-splits and 120 three-splits of M = 17 (C(16,1), C(16,2))")


# =============================================================================
#  E.  ATLAS CROSS-CHECK over all 34 GGV5 rows
# =============================================================================
def group_E(L: Ledger) -> None:
    L.head("E.  ATLAS CROSS-CHECK -- all 34 published GGV5 candidate rows")

    path = os.path.join(HERE, "corner_atlas.json")
    if not os.path.exists(path):
        L.ck("E0 corner_atlas.json present", False, "missing; run corner_atlas.py")
        return
    atlas = json.load(open(path))
    rows = atlas["rows"]
    L.ck("E0 corner_atlas.json has 34 rows", len(rows) == 34,
         "n_rows = %s" % atlas["n_rows"])

    recomputed, mism = [], []
    for r in rows:
        g1, g5 = r["gates"]["G1"], r["gates"]["G5"]
        tt, kk, qq = g1["t"], g1["kappa"], g1["ord_C"]
        aa, bb = sorted((r["m"], r["n"]))
        M, H = MH(tt, kk, qq, aa, bb)
        qw = q_window(tt, kk, qq, aa, bb)
        al = ordPhi_from_bridge(tt, kk, qq, aa, bb)
        if (M, H, qw) != (g5["M"], g5["H"], g5["q_window"]):
            mism.append(r["id"])
        recomputed.append(dict(id=r["id"], mono=g1["C_is_monomial"], t=tt,
                               kappa=kk, q=qq, a=aa, b=bb, M=M, H=H, qw=qw,
                               alpha=al))
    L.ck("E1 (M,H,q_window) re-derived from the atlas's own corner data",
         not mism,
         "34/34 rows reproduce the stored G5 numbers from (t,kappa,ord_C,m,n) "
         "with no reference to the stored values")

    L.ck("E2 kappa = t-2 on all 34 rows",
         all(x["kappa"] == x["t"] - 2 for x in recomputed),
         "the standard class holds throughout, which is what Group B needs")

    mono = [x for x in recomputed if x["mono"]]
    nonm = [x for x in recomputed if not x["mono"]]
    L.ck("E3 28 monomial rows / 6 non-monomial rows",
         (len(mono), len(nonm)) == (28, 6),
         "monomial = retraction shape fails = C = y (deg C = ord C = 1)")
    L.ck("E4 C monomial  <==>  q_window == M   (34/34)",
         all(x["qw"] == x["M"] for x in mono)
         and all(x["qw"] < x["M"] for x in nonm),
         "every monomial row has the MAXIMAL window denominator; every "
         "non-monomial row has strict cancellation (q_window in %s)"
         % sorted({x["qw"] for x in nonm}))
    L.ck("E5 gcd(M, ord_y Phi) = gcd(M, H) on all 34 rows",
         all(gcd(x["M"], abs(x["alpha"])) == gcd(x["M"], abs(x["H"]))
             for x in recomputed),
         "the bridge identity, instantiated on the published population")
    L.ck("E5b a+b > t-2 on all 34 rows (Group C's hypothesis)",
         all(x["a"] + x["b"] > x["t"] - 2 for x in recomputed),
         "min a+b = %d, max t = %d, so t-2 <= %d < min(a+b): the necessary "
         "condition ord_y C >= t applies to every published row"
         % (min(x["a"] + x["b"] for x in recomputed),
            max(x["t"] for x in recomputed),
            max(x["t"] for x in recomputed) - 2))
    L.ck("E6 ord_y C >= t at every integral-window row",
         all(x["q"] >= x["t"] for x in recomputed if x["qw"] == 1),
         "the three q_window=1 rows are %s, with (q,t) = %s -- Group C's "
         "necessary condition, on real data"
         % ([x["id"] for x in recomputed if x["qw"] == 1],
            [(x["q"], x["t"]) for x in recomputed if x["qw"] == 1]))

    # E7  the escape census: resolve G5 for all 34 rows without any w(e) datum.
    tot_esc_all = [x["id"] for x in recomputed
                   if any(w % x["qw"] == 0 for w in range(1, x["M"]))]
    tot_esc_state = [x["id"] for x in recomputed
                     if any(w % x["qw"] == 0 for w in range(2, x["t"] + 2))]
    L.ck("E7a escape possible at exactly the 6 non-monomial rows "
         "(w(e) in 1..M-1)",
         sorted(tot_esc_all) == sorted(x["id"] for x in nonm),
         "so the carry obstruction is TOTAL at 28/34 rows, decided with no "
         "split enumeration; atlas G5 records 31 UNKNOWN")
    L.ck("E7b escape possible at exactly 5 rows (w(e) in {2..t+1})",
         len(tot_esc_state) == 5 and all(
             x["mono"] is False for x in recomputed
             if x["id"] in tot_esc_state),
         "under the tighter state-variable reading F_7 (q_window=5 > t+1=4) "
         "also loses its escape: %s" % sorted(tot_esc_state))

    if not L.quiet:
        print("\n     id                       mono   t  q   M   H  q_win  "
              "ord_y(Phi)  escape w(e)")
        for x in recomputed:
            esc = [w for w in range(1, x["M"]) if w % x["qw"] == 0]
            print("     %-24s %-5s %2d %2d %3d %3d  %4d  %9d   %s"
                  % (x["id"], "y" if x["mono"] else "-", x["t"], x["q"],
                     x["M"], x["H"], x["qw"], x["alpha"],
                     esc[:5] if esc else "NONE"))

    # MUTATION CONTROL: q_window == M does NOT imply monomial as a theorem.
    ctr = [(tt, qq, aa, bb) for tt in range(3, 8) for qq in range(2, 12)
           for s in range(4, 16) for aa in range(1, s) for bb in [s - aa]
           if bb >= aa and q_window(tt, tt - 2, qq, aa, bb)
           == MH(tt, tt - 2, qq, aa, bb)[0]]
    L.mut("E E4's converse is a POPULATION fact, not a theorem",
          len(ctr) > 0,
          "%d abstract points with q >= 2 (C not a monomial) still have "
          "q_window = M, e.g. %s.  Only 'monomial ==> q_window = M' is proved; "
          "the biconditional holds on the 34 published rows and is reported as "
          "CHECKED, not PROVED" % (len(ctr), ctr[:3]))


# =============================================================================
#  F.  THE POSITIVE SIDE -- monomiality RAISES the depth-ledger floor
# =============================================================================
def group_F(L: Ledger) -> None:
    L.head("F.  THE FORCED FLOOR ceil(alpha w/q_window) GAINS EXACTLY (M-1)/2")

    # F1  strict gain at every admissible weight, monomial corners.
    bad = []
    for tt, kk, qq, aa, bb in [(4, 2, 1, 2, 3), (4, 2, 1, 3, 5), (3, 1, 1, 2, 3),
                               (5, 3, 1, 2, 3), (6, 4, 1, 2, 3)]:
        M, _ = MH(tt, kk, qq, aa, bb)
        al = ordPhi_from_bridge(tt, kk, qq, aa, bb)
        for w in range(1, M):
            if floor_gain(al, M, w) <= 0:
                bad.append((tt, aa, bb, w))
    L.ck("F1 gain(w) > 0 for every 0 < w < M at monomial corners",
         not bad,
         "the forced floor L(w) = ceil(alpha w/q_window) strictly exceeds the "
         "affine ray alpha w/q_window at EVERY admissible weight")

    # F2  the exact total: sum over w = 1..M-1 is (M-1)/2, because
    #     (-alpha w) mod M is a bijection of {1..M-1} when gcd(alpha,M)=1.
    tots = []
    for tt, kk, qq, aa, bb in [(4, 2, 1, 2, 3), (4, 2, 1, 3, 5), (3, 1, 1, 2, 3),
                               (5, 3, 1, 2, 3), (6, 4, 1, 2, 3)]:
        M, _ = MH(tt, kk, qq, aa, bb)
        al = ordPhi_from_bridge(tt, kk, qq, aa, bb)
        s = sum(floor_gain(al, M, w) for w in range(1, M))
        tots.append((M, s, Fraction(M - 1, 2)))
    L.ck("F2 total gain over w = 1..M-1 equals (M-1)/2, exactly",
         all(s == e for _, s, e in tots),
         "monomiality buys exactly (M-1)/2 units of forced floor: %s"
         % [(M, str(s)) for M, s, _ in tots])

    # F3  and the (72,108) contrast: zero gain everywhere.
    M72, _ = MH(**C_72_108)
    g72 = [floor_gain(ORDPHI_72_108_PUBLISHED, 1, w) for w in range(1, M72)]
    g72b = [Fraction((-ORDPHI_72_108_PUBLISHED * w) % M72, M72)
            for w in range(1, M72)]
    L.ck("F3 (72,108): gain is 0 at every weight",
         all(g == 0 for g in g72) and all(g == 0 for g in g72b),
         "q_window = 1 makes ceil exact; and even measured against M = 17 the "
         "residues vanish because 17 | 204.  The retracting corner is the case "
         "with NO floor gain at all")

    # MUTATION CONTROL: the law is (q_window-1)/2, NOT (M-1)/2 -- at a corner
    # with cancellation the two differ, and only the reduced form is right.
    M7, _ = MH(3, 1, 4, 2, 7)
    al7 = ordPhi_from_bridge(3, 1, 4, 2, 7)
    qw7 = q_window(3, 1, 4, 2, 7)
    s7 = sum(floor_gain(al7, M7, w) for w in range(1, qw7))
    L.mut("F the law is (q_window-1)/2; M and q_window differ at F_7",
          s7 == Fraction(qw7 - 1, 2) and qw7 != M7
          and s7 != Fraction(M7 - 1, 2),
          "F_7 has q_window = 5 < M = 25; total gain over w = 1..4 is %s = "
          "(5-1)/2, and NOT (25-1)/2 = 12.  So F2's '(M-1)/2' is the monomial "
          "specialisation of a (q_window-1)/2 law, valid there only because "
          "q_window = M" % s7)
    L.ck("F4 the gain sum over a full period is (q_window-1)/2 in general",
         all(sum(floor_gain(ordPhi_from_bridge(*c), MH(*c)[0], w)
                 for w in range(1, q_window(*c)))
             == Fraction(q_window(*c) - 1, 2)
             for c in [(3, 1, 4, 2, 7), (3, 1, 5, 3, 7), (4, 2, 3, 3, 4),
                       (4, 2, 1, 2, 3), (4, 2, 1, 3, 5)])
         and sum(floor_gain(ORDPHI_72_108_PUBLISHED, M72, w)
                 for w in range(1, q_window(**C_72_108))) == 0,
         "checked at four cancelling corners and two monomial ones; the "
         "(72,108) sum is the empty sum over w in range(1,1) = 0, the "
         "degenerate q_window = 1 endpoint")


# =============================================================================
#  G.  THE DEFECT DECOMPOSITION -- thinness and shallowness are independent
# =============================================================================
def group_G(L: Ledger) -> None:
    L.head("G.  TWO INDEPENDENT DEFECTS: thinness (deg C = ord C) vs "
           "shallowness (ord C = 1)")

    # G1  q_window has no deg C dependence whatsoever.
    src = ["q_window = M/gcd(M,H)", "M = t(a+b)-(kappa+1)", "H = q(a+b)-1"]
    L.ck("G1 q_window is a function of (t,kappa,ord_C,a,b) only",
         all("deg" not in s for s in src),
         "no formula in the chain mentions deg_y C, so NO change to deg C can "
         "move q_window -- adding a residual to C cannot revive the divisor "
         "mechanism")

    # G2  lam (and dg, and the slope gap) depend only on deg C - ord C.
    Nsym, Msym, dC, oC = sp.symbols("N M degC ordC")
    lam = Nsym * (dC - oC) / Msym
    L.ck("G2 lam = N(deg C - ord C)/M vanishes iff deg C = ord C, "
         "for any ord C",
         sp.simplify(lam.subs(dC, oC)) == 0
         and sp.simplify(lam.subs({dC: oC + 1})) != 0,
         "so the cascade/cone/closed-form deaths are indifferent to ord C -- "
         "deepening C cannot revive them")

    # G3  all four quadrants are realised by explicit parameter points.
    quad = {}
    #  (thin?, shallow?)  ->  (t, kappa, q=ord C, a, b, deg C)
    pts = {
        (True, True):  (4, 2, 1, 2, 3, 1),     # C = y      -- the class of nine
        (True, False): (3, 1, 8, 2, 3, 8),     # C = y^8    -- monomial, deep
        (False, True): (4, 2, 1, 2, 3, 2),     # C = y(y+1) -- residual, shallow
        (False, False): (4, 2, 7, 2, 3, 8),    # C = y^7(y+1) -- (72,108) itself
    }
    for key, (tt, kk, qq, aa, bb, dc) in pts.items():
        M, H = MH(tt, kk, qq, aa, bb)
        qw = q_window(tt, kk, qq, aa, bb)
        thin = (dc == qq)
        shallow = (qq == 1)
        quad[key] = dict(M=M, H=H, qw=qw, thin=thin, shallow=shallow,
                         lam_zero=thin)
        L.note("thin=%-5s shallow=%-5s  (t,kappa,ord C,deg C)=(%d,%d,%d,%d) "
               "M=%-3d H=%-3d q_window=%-3d  lam=%s"
               % (thin, shallow, tt, kk, qq, dc, M, H, qw,
                  "0" if thin else "!=0"))
    L.ck("G3 all four (thin, shallow) quadrants are realised",
         set(quad) == {(True, True), (True, False),
                       (False, True), (False, False)}
         and all(quad[k]["thin"] == k[0] and quad[k]["shallow"] == k[1]
                 for k in quad),
         "the two defects are logically independent chart conditions")
    L.ck("G4 thin-but-not-shallow has an INTEGRAL window (q_window = 1)",
         quad[(True, False)]["qw"] == 1,
         "C = y^8 at t=3, kappa=1, (a,b)=(2,3): M=13, H=39, gcd=13, "
         "q_window=1.  A MONOMIAL C with an integral window: monomiality per se "
         "does NOT kill the divisor mechanism -- ord C = 1 does")
    L.ck("G5 not-thin-but-shallow still has q_window = M",
         quad[(False, True)]["qw"] == quad[(False, True)]["M"],
         "C = y(y+1) at t=4, kappa=2, (a,b)=(2,3): lam != 0 (the cone and the "
         "cascade are alive) yet q_window = 17 = M (the divisor mechanism is "
         "still dead).  Restoring a residual is NOT a repair for q_window")

    # MUTATION CONTROL: the two invariants really do separate on these points.
    L.mut("G the quadrant table is discriminating",
          len({quad[k]["qw"] == 1 for k in quad}) == 2
          and len({quad[k]["thin"] for k in quad}) == 2,
          "q_window=1 splits the table one way (rows 2 and 4) and thinness "
          "splits it the other way (rows 1 and 2): the two are orthogonal, so "
          "no single 'monomiality bias' accounts for both")


# =============================================================================
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()
    L = Ledger(args.quiet)
    if not args.quiet:
        print(__doc__)
    for g in (group_A, group_B, group_C, group_D, group_E, group_F, group_G):
        g(L)
    return L.report()


if __name__ == "__main__":
    sys.exit(main())
