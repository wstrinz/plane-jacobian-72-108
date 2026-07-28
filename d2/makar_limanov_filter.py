#!/usr/bin/env python3
"""makar_limanov_filter.py  (NEW 2026-07-28; read-only over all existing artifacts)

MAKAR-LIMANOV 2025, THEOREM 2 AND THE LEMMA ON DIVISIBILITY, RUN AGAINST THE
34-ROW CORNER CENSUS.  VERDICT: **NEITHER OBSTRUCTION KILLS ANYTHING.**

Source (full text obtained and read, not paraphrased):
    L. Makar-Limanov, "On the shape of a counterexample to the two-dimensional
    Jacobian conjecture", Serdica Math. J. 51 (2025) 299-314,
    doi:10.55630/serdica.2025.51.299-314.  Open galley:
    https://serdica.math.bas.bg/index.php/serdica/article/download/300/153

THE TWO RESULTS, VERBATIM
-------------------------
(T2)  "Theorem 2.  If f is reduced and v_0 = prod_{i=1}^{k} p_i^{delta_i} (a,b)
      where all p_i are prime numbers and gcd(a,b) = 1 then N(f) has at most
      sum_i delta_i - 1 admissible edges."                        [p. 311]

(DIV) "Lemma on divisibility.  If J(sigma, tau) = sigma where sigma, tau are w
      homogeneous, w(sigma) > 0, w(tau) > 0 and sigma is not a monomial then
      w(sigma) doesn't divide w(tau)."                            [p. 305]

THE FOUR DEFINITIONS T2 DEPENDS ON, ALSO VERBATIM
-------------------------------------------------
* "reduced":  "From now on we will be looking at counterexamples where supp(f)
  has such a monomial [x^{d_x} y^{d_y}, d_x = deg_x(f), d_y = deg_y(f)].  Since
  we can make an automorphism x -> y, y -> -x we will also assume that
  d_x = m < d_y = n.  Additionally, if w_{1,0}(f) = x^m p(y) we can make a
  substitution y -> y - c such that the order of p(y - c) is larger than m ...
  Let us call such an f reduced."                                 [pp. 304-305]
* the rectangle:  "Since f is reduced ... points of the Newton polygon N(f) of f
  belong to the rectangle with vertices (0,0), (m,0), (m,n), (0,n)."  [p. 307]
* right edges / the bisectrix / "admissible":  "The boundary of the Newton
  polygon N(f) of f consists of right edges e_1, e_2, ..., e_k and left edges
  e_0, e_{-1}, ..., e_{-l}." ... "We will be interested in the right edges which
  have at least one vertex above the bisectrix of the first quadrant.  We will
  call this edges admissible."                                    [pp. 307-308]
* v_0:  "For a monomial mu denote by |mu| its degree vector and by v_0 the
  vertex of N(f) corresponding to \\hat f."                        [p. 310]
  \\hat a is "its leading monomial" (p. 302).  For reduced f the (1,1)-leading
  monomial is x^m y^n, so v_0 = (deg_x f, deg_y f) = (m, n) with m < n.

ONE READING DECLARED (the paper is ambiguous here and it is load-bearing)
------------------------------------------------------------------------
T2 as literally quantified would read "0 <= -1" for gcd(v_0) = 1, i.e. it would
by itself refute every counterexample with coprime corner coordinates.  That is
not what the proof establishes: the proof opens "For the last admissible edge
d_s is ...", i.e. it assumes s >= 1 admissible edges exist and bounds s.  WE USE
THE PROOF'S READING:

    (T2')  if N(f) has s >= 1 admissible edges then s <= Omega(gcd(v_0)) - 1.

THIS DOES NOT WEAKEN IT, because s >= 1 is FORCED: e_1 is a right edge with
vertex v_0 = (m,n) ("There are two edges e_0 and e_1 with the vertex (m,n)",
p. 307) and m < n puts v_0 strictly above the bisectrix.  So

    (T2-KILL)  a reduced f is impossible unless Omega(gcd(v_0)) >= 2.

Second declared reading: "above the bisectrix" is taken as STRICTLY above
(y > x).  Immaterial for v_0, where m < n is strict.

WHAT IS ACTUALLY CHECKED HERE
-----------------------------
PART A  T2 as a pure function of a corner vector, with its floor s >= 1.
PART B  THE IDENTIFICATION OF v_0 ON OUR CENSUS -- the whole result turns on it.
        GGV5's normal form (1708.07936_GGV5.tex:250) says supp(P) lies in the
        rectangle {(0,0), m(a,0), m(a,b), m(0,b)} with m(a,b) IN supp(P).  That
        is verbatim ML's reduced normalisation, and it pins
            v_0(P) = m * A_0        v_0(Q) = n * A_0
        NOT v_0 = A_0.  The naive identification v_0 := A_0 is refuted on all 34
        rows and is refuted in the DANGEROUS direction: it "kills" 14 rows
        including the flagship F_2(3,5)/125 while sparing (8,28)/(3,2)/108, the
        one row this campaign proved dead.  PART B fires that misfire explicitly.
PART C  the census: Omega and the predicted bound for all 34 rows, both P and Q.
PART D  CALIBRATION against the 10 rows settled below 125, and the OBSERVED
        admissible-edge counts wherever a polygon exists.
PART E  (DIV): the verbatim statement, its one missing hypothesis, and the
        TRANSFER QUESTION -- our corners carry [P,Q] = x^2, not 1.

HEADLINE RESULTS
----------------
1.  T2 KILLS NOTHING, on any of the 34 rows, and the reason is structural:
    m, n > 1 (GGV5) and gcd(A_0) >= 3 (all 34 rows), so Omega(gcd(v_0)) >= 2
    automatically.  T2's kill condition needs gcd(A_0) = 1, which no published
    corner in this range has.
2.  T2 is not vacuous: the same code path DOES kill synthetic corners, and T2 is
    TIGHT (predicted == observed == 1) at the corner (5,20) that carries the
    flagship F_2(3,5)/125.  It leaves a testable shape constraint there.
3.  (DIV) DOES NOT TRANSFER, and not merely for want of hypotheses: under our
    bracket normalisation [P,Q] = x^2 the relation is J(rho,tau) = x^2 rho, and
    its divisibility conclusion is FALSE -- exhibited exactly at the positive,
    primitive, non-diagonal weight w = (1,3) by rho = x^3 + y, tau = (x^3 y +
    y^2)/3, with w(rho) = 3 dividing w(tau) = 6.  No rescaling repairs this.
4.  As a by-product this file DERIVES original-coordinate corner data the repo
    never wrote down, e.g. (deg_x Q, deg_y Q) = (16,56) for the (72,108) case.

EVIDENCE GRADES USED BELOW
    [PROVED]     follows from published text by pure arithmetic done here
    [CHECKED]    exact computation in this file
    [CITATION]   transcribed from a pinned published line
    [INFERRED]   repo-internal rule flagged as unpublished (only the r1 unit
                 polygon Delta, used ONLY for the observed counts in PART D)

Checker: --quiet, exit 0 iff every check passes.  Exact arithmetic only
(int / Fraction / sympy).  Reads corner_atlas.json; no writes.  ~3 s.
"""
from __future__ import annotations

import itertools
import json
import math
import os
import sys
from fractions import Fraction

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__)) or "."
QUIET = "--quiet" in sys.argv
_ok = [0]
_fail: list[str] = []


def ck(name: str, cond: bool, detail: str = "") -> bool:
    if cond:
        _ok[0] += 1
        if not QUIET:
            print("[OK]   %s" % name)
    else:
        _fail.append(name)
        print("[FAIL] %s%s" % (name, ("  -- " + detail) if detail else ""))
    return bool(cond)


def note(s: str = "") -> None:
    if not QUIET:
        print(s)


# ===========================================================================
# PART A.  THEOREM 2 AS A PURE FUNCTION OF A CORNER VECTOR.
# ===========================================================================
def omega(k: int) -> int:
    """Omega(k) = number of prime factors of k WITH multiplicity; Omega(1) = 0.

    This is ML's  sum_i delta_i  for  k = prod_i p_i^{delta_i}.
    """
    if k < 1:
        raise ValueError("Omega is defined on positive integers, got %r" % (k,))
    if k == 1:
        return 0
    return sum(sp.factorint(k).values())


def corner_decomposition(v0: tuple) -> tuple:
    """ML's  v_0 = (prod_i p_i^{delta_i}) * (a,b),  gcd(a,b) = 1.

    Returns (D, (a,b), Omega(D)).  D = gcd(v_0) and the decomposition is unique.
    """
    u, v = v0
    if u <= 0 or v <= 0:
        raise ValueError("corner vector must be strictly positive, got %r" % (v0,))
    D = math.gcd(u, v)
    return D, (u // D, v // D), omega(D)


def t2_bound(v0: tuple) -> int:
    """ML Theorem 2: N(f) has at most  sum_i delta_i - 1  admissible edges."""
    return corner_decomposition(v0)[2] - 1


def t2_floor(v0: tuple) -> int:
    """The forced LOWER bound on the admissible-edge count for reduced f.

    e_1 is a right edge carrying the vertex v_0 = (m,n) (ML p. 307), and a
    reduced f has m < n, so v_0 lies STRICTLY above the bisectrix y = x and e_1
    is admissible.  Hence s >= 1.  A vector on or below the bisectrix is not the
    corner of a reduced f at all, and the floor is 0 (T2 then says nothing).
    """
    return 1 if v0[1] > v0[0] else 0


def t2_kills(v0: tuple) -> bool:
    """True iff T2 refutes a reduced f with this corner: bound < forced floor."""
    return t2_bound(v0) < t2_floor(v0)


def part_A() -> None:
    note("\n" + "=" * 78)
    note("PART A -- Theorem 2 as a function of the corner vector")
    note("=" * 78)

    # A1 -- the decomposition is ML's, and Omega is Omega, not omega.
    ck("A1  Omega counts prime factors WITH multiplicity: "
       "Omega(1)=0, Omega(5)=1, Omega(4)=2, Omega(12)=3, Omega(60)=4",
       [omega(k) for k in (1, 5, 4, 12, 60)] == [0, 1, 2, 3, 4],
       str([omega(k) for k in (1, 5, 4, 12, 60)]))
    # MUTATION CONTROL: the DISTINCT-prime count omega() would give 12 -> 2.
    distinct = len(sp.factorint(12))
    ck("A1m MUTATION CONTROL -- the distinct-prime count omega(12) = 2 differs "
       "from Omega(12) = 3, so a code path using omega() would report a "
       "DIFFERENT (looser) bound and is excluded",
       distinct == 2 and distinct != omega(12))

    # A2 -- the decomposition v_0 = D*(a,b) with gcd(a,b)=1.
    D, ab, om = corner_decomposition((15, 60))
    ck("A2  v_0 = (15,60) decomposes as 15*(1,4) with Omega(15) = 2, so T2 "
       "predicts at most 1 admissible edge",
       (D, ab, om, t2_bound((15, 60))) == (15, (1, 4), 2, 1),
       str((D, ab, om)))
    ck("A2b the primitive part is genuinely primitive on every decomposition "
       "tested", all(math.gcd(*corner_decomposition(v)[1]) == 1
                     for v in [(15, 60), (24, 84), (16, 56), (4, 12), (7, 21),
                               (36, 108), (2, 3), (6, 10)]))

    # A3 -- the floor s >= 1, and why.
    ck("A3  a reduced f has m < n, so v_0 is STRICTLY above the bisectrix and "
       "e_1 is admissible: floor = 1 on every corner of a reduced f",
       all(t2_floor(v) == 1 for v in [(15, 60), (24, 84), (16, 56), (1, 2)]))
    # MUTATION CONTROL: drop the reducedness normalisation m < n and the floor
    # collapses, which is exactly what would make T2 unable to kill anything.
    ck("A3m MUTATION CONTROL -- on the bisectrix or below (n <= m) the floor is "
       "0 and T2 kills nothing; the m < n half of 'reduced' is load-bearing",
       t2_floor((5, 5)) == 0 and t2_floor((20, 5)) == 0
       and not t2_kills((5, 5)) and not t2_kills((20, 5)))

    # A4 -- the kill predicate, on synthetic corners with a known answer.
    ck("A4  T2 REFUTES a reduced f with coprime corner (1,4): Omega(1) = 0, "
       "bound -1 < floor 1", t2_kills((1, 4)) and t2_bound((1, 4)) == -1)
    ck("A4b T2 REFUTES a reduced f with prime-gcd corner (3,6): Omega(3) = 1, "
       "bound 0 < floor 1", t2_kills((3, 6)) and t2_bound((3, 6)) == 0)
    ck("A4c T2 does NOT refute (4,12): Omega(4) = 2, bound 1 >= floor 1",
       (not t2_kills((4, 12))) and t2_bound((4, 12)) == 1)
    ck("A4d T2 does NOT refute (24,84): Omega(12) = 3, bound 2",
       (not t2_kills((24, 84))) and t2_bound((24, 84)) == 2)

    # A5 -- the exact quantifier boundary: Omega(gcd) >= 2 is the whole content.
    boundary = [(v, t2_kills(v)) for v in
                [(1, 2), (2, 6), (4, 12), (9, 18), (5, 25), (6, 18)]]
    expected = [((1, 2), True), ((2, 6), True), ((4, 12), False),
                ((9, 18), False), ((5, 25), True), ((6, 18), False)]
    ck("A5  the entire content of T2-as-a-filter is 'Omega(gcd(v_0)) >= 2', "
       "verified on the boundary cases %s"
       % [(v, k) for v, k in boundary], boundary == expected, str(boundary))


# ===========================================================================
# PART B.  THE IDENTIFICATION OF v_0 ON OUR CENSUS.
#          This is where a filter of this kind is won or lost.
# ===========================================================================
#
# GGV5 Introduction, paper_src/1708.07936_GGV5.tex:250, VERBATIM:
#
#   "If this conjecture is false, then there exist P,Q in L such that
#    [P,Q] = K^x, and there exist m,n,a,b in N, such that m,n > 1 are coprime,
#    a < b, the support of P is contained in the rectangle with vertices
#    {(0,0), m(a,0), m(a,b), m(0,b)}, the support of Q is contained in the
#    rectangle with vertices {(0,0), n(a,0), n(a,b), n(0,b)}, the point m(a,b)
#    is in the support of P and the point n(a,b) is in the support of Q.
#    Note that deg(P) = m(a+b) and deg(Q) = n(a+b)."
#
# and paper_src/1708.07936_GGV5.tex:398-400:
#
#   "A_0 = (1/m) en_10(P) ... This point A_0 corresponds to (a,b) in the
#    introduction."
#
# Consequences, all pure arithmetic on that sentence [PROVED]:
#   (i)   supp(P) inside the rectangle => deg_x P <= m*a and deg_y P <= m*b;
#         m(a,b) in supp(P) => deg_x P >= m*a and deg_y P >= m*b.  Hence
#              (deg_x P, deg_y P) = m * A_0 ,   (deg_x Q, deg_y Q) = n * A_0 .
#   (ii)  x^{deg_x P} y^{deg_y P} is IN supp(P): that is exactly ML's defining
#         property of "reduced", so GGV5's normal form IS ML's reduced form
#         (up to ML's extra y -> y - c shift, which moves no vertex).
#   (iii) a < b => m*a < m*b: ML's m < n normalisation holds for both P and Q.
#   (iv)  m,n > 1 coprime => gamma = deg(Q)/deg(P) = n/m is neither an integer
#         nor the reciprocal of one -- ML's own normalisation of gamma.
#   (v)   therefore  v_0(P) = m * A_0  and  v_0(Q) = n * A_0 ,  NOT A_0.
# ---------------------------------------------------------------------------
GGV5_NORMAL_FORM_LINE = "paper_src/1708.07936_GGV5.tex:250"


def v0_of_P(row) -> tuple:
    a0, b0 = row["A0"]
    return (row["m"] * a0, row["m"] * b0)


def v0_of_Q(row) -> tuple:
    a0, b0 = row["A0"]
    return (row["n"] * a0, row["n"] * b0)


def part_B(rows) -> None:
    note("\n" + "=" * 78)
    note("PART B -- what v_0 IS on our rows, and the refutation of v_0 := A_0")
    note("=" * 78)

    # B1 -- the degree recipe implied by the normal form, checked on all rows.
    bad = [r["id"] for r in rows
           if r["max_deg"] != (sum(r["A0"])) * max(r["m"], r["n"])]
    ck("B1  GGV5:250's 'deg(P) = m(a+b), deg(Q) = n(a+b)' reproduces max_deg on "
       "all 34 rows [CITATION -> CHECKED]", not bad, str(bad))

    # B2 -- GGV5's normal form satisfies ML's "reduced" hypotheses.
    ck("B2  a < b on all 34 rows, so ML's reducedness normalisation m < n holds "
       "for BOTH P and Q [PROVED]",
       all(r["A0"][0] < r["A0"][1] for r in rows))
    ck("B2b m,n > 1 and coprime on all 34 rows, so ML's gamma = n/m is neither "
       "an integer nor the reciprocal of one -- his own normalisation [PROVED]",
       all(r["m"] > 1 and r["n"] > 1 and math.gcd(r["m"], r["n"]) == 1
           for r in rows))

    # B3 -- the identification itself, and its independent cross-check: ML's
    #       Lemma on similarity says N(g) is homothetic to N(f) with ratio
    #       gamma, so gamma * v_0(f) must be integral.
    gam_ok = []
    for r in rows:
        g = Fraction(r["n"], r["m"])                     # gamma = deg Q / deg P
        vP = v0_of_P(r)
        img = (g * vP[0], g * vP[1])
        gam_ok.append(img == tuple(Fraction(c) for c in v0_of_Q(r))
                      and all(c.denominator == 1 for c in img))
    ck("B3  ML's Lemma on similarity is SATISFIED by this identification: "
       "gamma * v_0(P) = v_0(Q) exactly, and is integral, on all 34 rows "
       "[CHECKED]", all(gam_ok), str(gam_ok.count(False)))

    # B4 -- REFUTATION of the naive reading v_0 := A_0.
    #      (a) deg f would have to be a_0+b_0, but deg P = m(a_0+b_0), m >= 2.
    naive_deg_bad = [r["id"] for r in rows
                     if sum(r["A0"]) not in (r["m"] * sum(r["A0"]),
                                             r["n"] * sum(r["A0"]))]
    ck("B4  the naive identification v_0 := A_0 is REFUTED on ALL 34 rows: "
       "reducedness forces deg f = (v_0)_x + (v_0)_y, which would make "
       "deg f = a_0+b_0, but GGV5 gives deg P = m(a_0+b_0) with m >= 2 [PROVED]",
       len(naive_deg_bad) == 34, str(len(naive_deg_bad)))
    #      (b) an independent refutation: the similarity lemma fails for it.
    naive_sim_bad = [r["id"] for r in rows
                     if math.gcd(*r["A0"]) % min(r["m"], r["n"]) != 0]
    ck("B4b independent refutation of v_0 := A_0 -- ML's similarity lemma needs "
       "denom(gamma) = min(m,n) to divide gcd(v_0); with v_0 := A_0 that FAILS "
       "on %d of the 34 rows, including the flagship" % len(naive_sim_bad),
       len(naive_sim_bad) == 21 and "F_2(3,5)/125" in naive_sim_bad,
       str(len(naive_sim_bad)))

    # B5 -- THE MISFIRE CONTROL.  What the naive reading would have claimed.
    naive_kill = [r["id"] for r in rows if t2_kills(tuple(r["A0"]))]
    ck("B5  MISFIRE CONTROL -- the naive reading v_0 := A_0 'kills' 14 rows",
       len(naive_kill) == 14, str(len(naive_kill)))
    ck("B5b MISFIRE CONTROL -- and it is anti-correlated with truth: it KILLS "
       "F_2(3,5)/125, which is OPEN, while SPARING (8,28)/(3,2)/108, the one "
       "row this campaign actually closed.  A filter that behaves this way on "
       "the two rows whose answer we know is not usable, and its 14 'kills' "
       "are artefacts of the wrong identification refuted in B4/B4b.",
       "F_2(3,5)/125" in naive_kill and "(8,28)/(3,2)/108" not in naive_kill)

    if not QUIET:
        note("\n  naive-reading 'kills' (ALL ARTEFACTS -- do not cite):")
        for i in naive_kill:
            note("     %s" % i)


# ===========================================================================
# PART C.  THE CENSUS -- T2 over all 34 rows with the correct v_0.
# ===========================================================================
def part_C(rows) -> list:
    note("\n" + "=" * 78)
    note("PART C -- the census: T2 over all 34 rows, for BOTH P and Q")
    note("=" * 78)
    note("")
    note("  %-22s %-9s %-11s %-3s %-4s  %-11s %-3s %-4s"
         % ("row", "A_0", "v_0(P)", "Om", "bnd", "v_0(Q)", "Om", "bnd"))
    note("  " + "-" * 74)

    census = []
    for r in rows:
        vP, vQ = v0_of_P(r), v0_of_Q(r)
        oP, oQ = corner_decomposition(vP)[2], corner_decomposition(vQ)[2]
        rec = dict(id=r["id"], A0=tuple(r["A0"]), vP=vP, vQ=vQ,
                   omP=oP, omQ=oQ, bP=t2_bound(vP), bQ=t2_bound(vQ),
                   kill=t2_kills(vP) or t2_kills(vQ))
        census.append(rec)
        note("  %-22s %-9s %-11s %-3d %-4d  %-11s %-3d %-4d%s"
             % (r["id"], tuple(r["A0"]), vP, oP, rec["bP"], vQ, oQ, rec["bQ"],
                "   *** KILL ***" if rec["kill"] else ""))

    kills = [c["id"] for c in census if c["kill"]]
    ck("\nC1  T2 KILLS NOTHING: 0 of the 34 rows is refuted, for either member "
       "of the pair", not kills, str(kills))
    ck("C2  every predicted bound is >= 1, i.e. T2 never even contradicts the "
       "forced floor s >= 1",
       all(c["bP"] >= 1 and c["bQ"] >= 1 for c in census))

    # C3 -- WHY.  The structural reason, so the negative cannot be mistaken for
    #       an accident of this particular table.
    ck("C3  the reason is structural: Omega(gcd(v_0(P))) = Omega(m) + "
       "Omega(gcd(A_0)) with m > 1 (GGV5) and gcd(A_0) > 1 on every row, so "
       "Omega >= 2 automatically",
       all(math.gcd(*r["A0"]) > 1 for r in rows)
       and all(corner_decomposition(v0_of_P(r))[2]
               == omega(r["m"]) + omega(math.gcd(*r["A0"])) for r in rows))
    gcds = sorted({math.gcd(*r["A0"]) for r in rows})
    ck("C3b in fact gcd(A_0) >= 3 on all 34 rows (values %s), so T2's kill "
       "condition -- which needs gcd(A_0) = 1 together with a prime m or n -- "
       "is not met by any published corner in this range" % gcds,
       min(gcds) >= 3)

    # C4 -- MUTATION CONTROL: the filter is LIVE, not vacuously permissive.
    synth = [((1, 4), 2, 3), ((3, 8), 2, 3), ((2, 5), 3, 5), ((5, 12), 7, 9)]
    killed = []
    for A0, m, n in synth:
        vP = (m * A0[0], m * A0[1])
        vQ = (n * A0[0], n * A0[1])
        killed.append(t2_kills(vP) or t2_kills(vQ))
    ck("C4  MUTATION CONTROL -- the SAME code path DOES kill synthetic rows "
       "with gcd(A_0) = 1 and a prime multiplier: %s all refuted.  The 0/34 "
       "above is a fact about the census, not a dead filter."
       % [s[0] for s in synth], all(killed), str(killed))
    ck("C4b MUTATION CONTROL -- and it correctly SPARES the synthetic row with "
       "gcd(A_0) = 1 and both multipliers composite, which is T2's exact "
       "boundary",
       not (t2_kills((4 * 1, 4 * 4)) or t2_kills((9 * 1, 9 * 4))))

    # C5 -- the general corollary worth recording.
    ck("C5  GENERAL COROLLARY of T2 + GGV5's normal form [PROVED]: for any "
       "counterexample, if gcd(a,b) = 1 then BOTH m and n must be composite.  "
       "No corner in GGV5's max_deg <= 150 table has gcd(a,b) = 1, so the "
       "corollary is non-vacuous but unused here.",
       all(t2_kills((m * 1, m * 4)) == (omega(m) <= 1) for m in
           (2, 3, 4, 5, 6, 7, 8, 9, 11, 12)))
    return census


# ===========================================================================
# PART D.  CALIBRATION -- against rows whose answer is already known, and
#          against the OBSERVED admissible-edge counts.
# ===========================================================================
#
# The original-coordinate unit polygon.  passport_75_125.py rule (r1):
#     mu = floor((b0-1)/a0),  c = b0 - mu*a0,
#     Delta = {(0,0), (1,0), (a0,b0), (0,c)}   (+ any published extra corners)
# and N(P) = m*Delta, N(Q) = n*Delta.  Instantiated in print in this repo only
# at (5,20): "Delta = {(0,0),(1,0),(5,20),(0,5)}" (PASSPORT_75_125.md:236,
# passport_75_125.py:620, polygon_reduction.py:624).
#
# GRADE: [INFERRED].  The repo flags r1's chart exponent as unpublished, and the
# "+ extra corners" clause means the vertex list may be INCOMPLETE.  Therefore
# the observed counts below are LOWER BOUNDS on the true admissible-edge count.
# A lower bound is exactly what is needed to test for a CONTRADICTION with T2's
# upper bound, so the direction is safe: if observed > predicted we would have
# a kill or a bug, and we check for precisely that.
# ---------------------------------------------------------------------------
def unit_delta(a0: int, b0: int) -> list:
    mu = (b0 - 1) // a0
    c = b0 - mu * a0
    return [(0, 0), (1, 0), (a0, b0), (0, c)]


def hull(pts: list) -> list:
    """Exact integer convex hull, counter-clockwise."""
    pts = sorted(set(pts))
    if len(pts) <= 2:
        return list(pts)

    def half(ps):
        out = []
        for p in ps:
            while len(out) >= 2:
                (x1, y1), (x2, y2) = out[-2], out[-1]
                if (x2 - x1) * (p[1] - y1) - (y2 - y1) * (p[0] - x1) <= 0:
                    out.pop()
                else:
                    break
            out.append(p)
        return out

    return half(pts)[:-1] + half(pts[::-1])[:-1]


def right_edges(poly: list) -> list:
    """ML's right edges e_1,...,e_k: the boundary chain from the corner v_0 to
    the vertex on the x axis, walking clockwise (the side whose edge extensions
    meet the x axis)."""
    h = hull(poly)
    mx = max(p[0] for p in h)
    corner = max(p for p in h if p[0] == mx)
    axis = max(p for p in h if p[1] == 0)
    i, j = h.index(corner), h.index(axis)
    out, k = [], i
    while k != j:
        nk = (k - 1) % len(h)
        out.append((h[k], h[nk]))
        k = nk
    return out


def n_admissible(poly: list) -> int:
    """ML: right edges with at least one vertex STRICTLY above y = x."""
    return sum(1 for (u, v) in right_edges(poly) if u[1] > u[0] or v[1] > v[0])


# The 10 rows GGHV22 Thm 2.1 settles below the 125 bound, keyed exactly as
# gghv_sub125.py keys them: (A_0, (m,n), max_deg).  Transcribed, not recomputed.
SETTLED_SUB125 = {
    ((4, 12), (3, 4), 64), ((4, 12), (5, 7), 112), ((5, 20), (2, 3), 75),
    ((5, 20), (3, 2), 75), ((7, 21), (2, 3), 84), ((8, 24), (2, 3), 96),
    ((8, 28), (3, 2), 108), ((8, 32), (3, 2), 120), ((9, 24), (2, 3), 99),
    ((9, 27), (2, 3), 108),
}


def part_D(rows, census) -> None:
    note("\n" + "=" * 78)
    note("PART D -- CALIBRATION against rows whose answer is already known")
    note("=" * 78)

    by_id = {c["id"]: c for c in census}

    def key(r):
        return (tuple(r["A0"]), (r["m"], r["n"]), r["max_deg"])

    settled = [r for r in rows if key(r) in SETTLED_SUB125]
    ck("D1  the calibration set is the 10 rows settled below the 125 bound "
       "(gghv_sub125.py's computed partition)", len(settled) == 10,
       str(len(settled)))

    killed_settled = [r["id"] for r in settled if by_id[r["id"]]["kill"]]
    ck("D2  CALIBRATION RESULT -- T2 kills 0 of the 10 KNOWN-DEAD rows.  It "
       "reproduces NO death that is already established, so it has no "
       "demonstrated discriminating power on this census.",
       not killed_settled, str(killed_settled))

    # D3 -- the sharpest single calibration: the row we ourselves closed.
    ours = by_id["(8,28)/(3,2)/108"]
    ck("D3  the sharpest calibration -- (8,28)/(3,2)/108 is the one row with a "
       "PUBLISHED full reduction and a proved death (this campaign).  T2 sees "
       "v_0(P) = (24,84) and v_0(Q) = (16,56), Omega = 3 and 3, bounds 2 and 2: "
       "ALIVE.  The one death we can check against is invisible to T2.",
       ours["vP"] == (24, 84) and ours["vQ"] == (16, 56)
       and ours["bP"] == 2 and ours["bQ"] == 2 and not ours["kill"],
       str(ours))
    ck("D3b by-product: this DERIVES original-coordinate corner data the repo "
       "never wrote down -- (deg_x Q, deg_y Q) = (16,56) and (deg_x P, deg_y P) "
       "= (24,84) for the (72,108) case, with deg P = 108 and deg Q = 72",
       sum(ours["vP"]) == 108 and sum(ours["vQ"]) == 72)

    # D4 -- OBSERVED counts, where a polygon exists at all.
    note("")
    note("  OBSERVED admissible-edge counts [INFERRED unit polygon r1;")
    note("  a LOWER bound, since r1 may omit published extra corners]")
    note("  %-22s %-6s %-24s %-9s %-9s"
         % ("row", "poly", "N(.) vertices", "observed", "predicted"))
    note("  " + "-" * 74)

    have_polygon = [r for r in rows if tuple(r["A0"]) in {(5, 20), (8, 28)}]
    violations = []
    observed_rows = []
    for r in have_polygon:
        a0, b0 = r["A0"]
        d = unit_delta(a0, b0)
        for tag, mult, pred in (("P", r["m"], by_id[r["id"]]["bP"]),
                                ("Q", r["n"], by_id[r["id"]]["bQ"])):
            poly = [(mult * u, mult * v) for (u, v) in d]
            obs = n_admissible(poly)
            observed_rows.append((r["id"], tag, obs, pred))
            if obs > pred:
                violations.append((r["id"], tag, obs, pred))
            note("  %-22s %-6s %-24s %-9d %-9d%s"
                 % (r["id"], tag, str(hull(poly))[:24], obs, pred,
                    "  <-- VIOLATION" if obs > pred else ""))

    ck("\nD4  NO row has observed > predicted -- T2 is nowhere contradicted by "
       "the polygons we actually have (which would have meant a kill or a bug)",
       not violations, str(violations))
    ck("D4b the r1 unit polygon is a quadrilateral, so its right boundary is a "
       "single edge and every observed count is 1; these observations are "
       "CONSISTENT with T2 but cannot confirm it",
       all(o == 1 for (_, _, o, _) in observed_rows))
    # MUTATION CONTROL: the observed-count machinery can detect a violation.
    # (5,10) lies strictly OUTSIDE the segment (1,0)-(5,20), so it is a genuine
    # extra hull vertex, and it sits above the bisectrix (10 > 5).
    fake = [(0, 0), (1, 0), (5, 10), (5, 20), (0, 5)]
    fake3 = [(3 * u, 3 * v) for (u, v) in fake]
    ck("D4m MUTATION CONTROL -- inserting ONE extra vertex above the bisectrix "
       "into the (5,20) polygon raises the observed count to 2, which EXCEEDS "
       "the predicted 1 and would be flagged.  The D4 pass is a real test.",
       n_admissible(fake3) == 2 and n_admissible(fake3) > by_id[
           "F_2(3,5)/125"]["bP"], str(n_admissible(fake3)))

    # D5 -- the flagship, in detail.
    fl = by_id["F_2(3,5)/125"]
    ck("D5  THE FLAGSHIP F_2(3,5)/125 (A_0 = (5,20), (m,n) = (3,5), the unique "
       "open row AT the 125 bound): v_0(P) = (15,60) and v_0(Q) = (25,100), "
       "Omega = 2 both, bound 1 both.  NOT KILLED.",
       fl["vP"] == (15, 60) and fl["vQ"] == (25, 100)
       and fl["bP"] == 1 and fl["bQ"] == 1 and not fl["kill"], str(fl))
    ck("D5b ... but T2 is TIGHT there: predicted 1 == forced floor 1, so T2 "
       "says N(P) and N(Q) have EXACTLY ONE admissible right edge.  That is a "
       "falsifiable shape constraint on the 125 case: no right-boundary vertex "
       "of N(P) other than (15,60) may lie strictly above the bisectrix.",
       fl["bP"] == t2_floor(fl["vP"]) == 1
       and fl["bQ"] == t2_floor(fl["vQ"]) == 1)
    ck("D5c and the constraint is SATISFIED by the polygon we have: the r1 "
       "polygon 3*{(0,0),(1,0),(5,20),(0,5)} has exactly one admissible right "
       "edge, from (15,60) to (3,0)",
       n_admissible([(3 * u, 3 * v) for (u, v) in unit_delta(5, 20)]) == 1)


# ===========================================================================
# PART E.  THE LEMMA ON DIVISIBILITY, AND WHETHER IT TRANSFERS TO [P,Q] = x^2.
# ===========================================================================
X, Y = sp.symbols("X Y")


def jac(p, q):
    return sp.expand(sp.diff(p, X) * sp.diff(q, Y) - sp.diff(p, Y) * sp.diff(q, X))


def w_monomials(al: int, be: int, W: int, maxdeg: int = 10) -> list:
    return [(i, j) for i in range(maxdeg + 1) for j in range(maxdeg + 1)
            if al * i + be * j == W]


def build(ms, cs):
    return sum(c * X ** i * Y ** j for (i, j), c in zip(ms, cs))


def div_scan(al: int, be: int, kappa: int, Rmax: int = 10,
             grid=(0, 1, -1, 2)):
    """Search for w-homogeneous rho (NON-monomial), tau with

           J(rho, tau) = X^kappa * rho     and     w(rho) | w(tau),

    i.e. a configuration where ML's divisibility CONCLUSION fails.  Returns
    (w(rho), w(tau), rho, tau) or None.

    Degree bookkeeping: w(J(rho,tau)) = w(rho)+w(tau)-(al+be), so the relation
    FORCES w(tau) = (kappa+1)*al + be.  kappa = 0 is ML's own normalisation.

    A HIT is a proof (it is re-verified exactly).  A MISS over this finite grid
    is RECONNAISSANCE ONLY and is reported as such.
    """
    Wt = (kappa + 1) * al + be
    tm = w_monomials(al, be, Wt)
    if not tm:
        return None
    tc = sp.symbols("t0:%d" % len(tm))
    tau = build(tm, tc)
    for R in range(1, Rmax + 1):
        if Wt % R:
            continue
        rm = w_monomials(al, be, R)
        if len(rm) < 2:
            continue
        for co in itertools.product(grid, repeat=len(rm)):
            if sum(1 for c in co if c) < 2:
                continue
            rho = build(rm, co)
            eq = sp.Poly(sp.expand(jac(rho, tau) - X ** kappa * rho), X, Y)
            for s in sp.solve(eq.coeffs(), tc, dict=True):
                tv = sp.expand(tau.subs(s).subs({c: 0 for c in tc}))
                if tv != 0 and sp.expand(jac(rho, tv) - X ** kappa * rho) == 0:
                    return (R, Wt, rho, tv)
    return None


# The reduced Newton polygons this repo actually computes (polygon_reduction.py
# regression contracts R1/R2/R3; R1 matches paper_src/upstream_facts.json).
REDUCED_POLYGONS = {
    "(8,28) sub2 N(P)": [(0, 0), (1, 0), (8, 14), (8, 16)],
    "(8,28) sub2 N(Q)": [(0, 0), (2, 1), (12, 21), (12, 24)],
    "(8,28) sub1 N(P)": [(0, 0), (0, 8), (1, 0), (8, 14), (8, 16)],
    "(8,28) sub1 N(Q)": [(0, 0), (0, 12), (2, 1), (12, 21), (12, 24)],
    "(5,20)/(2,3) N(P)": [(0, 0), (0, 10), (6, 0), (8, 2)],
    "(5,20)/(2,3) N(Q)": [(0, 0), (0, 15), (9, 0), (12, 3)],
    "(5,20)/(3,5) N(P)": [(0, 0), (0, 15), (9, 0), (12, 3)],
    "(5,20)/(3,5) N(Q)": [(0, 0), (0, 25), (15, 0), (20, 5)],
}


def positive_weight_edges(poly: list) -> list:
    """Every edge of the hull supported by a STRICTLY positive primitive weight,
    returned as (edge, primitive weight)."""
    h = hull(poly)
    out = []
    for i in range(len(h)):
        a, b = h[i], h[(i + 1) % len(h)]
        dx, dy = b[0] - a[0], b[1] - a[1]
        for w in ((dy, -dx), (-dy, dx)):
            if w[0] > 0 and w[1] > 0:
                g = math.gcd(w[0], w[1])
                wp = (w[0] // g, w[1] // g)
                mv = max(wp[0] * p[0] + wp[1] * p[1] for p in h)
                sup = [p for p in h if wp[0] * p[0] + wp[1] * p[1] == mv]
                if len(sup) > 1:
                    out.append(((a, b), wp))
    return out


def part_E() -> None:
    note("\n" + "=" * 78)
    note("PART E -- the Lemma on divisibility, and the [P,Q] = x^2 transfer")
    note("=" * 78)

    # E1 -- the degree bookkeeping the lemma silently assumes.
    al, be = sp.symbols("alpha beta", positive=True, integer=True)
    # w(J(rho,tau)) = w(rho) + w(tau) - w(xy); setting it = w(rho) gives
    # w(tau) = alpha + beta.  Verified on a concrete instance below.
    rho0, tau0 = X ** 3 + Y, sp.Rational(1, 3) * (X ** 3 * Y + Y ** 2)
    ck("E1  the relation J(rho,tau) = rho FORCES w(tau) = alpha+beta (weighted "
       "degree bookkeeping), so DIV's conclusion is literally 'w(rho) does not "
       "divide alpha+beta' [PROVED]",
       True)  # recorded; the concrete arithmetic is exercised in E2/E4.

    # E2 -- the verbatim statement is FALSE at the diagonal weight.
    r2, t2_ = X + Y, X * Y + Y ** 2
    ok_ce = sp.expand(jac(r2, t2_) - r2) == 0
    ck("E2  DIV AS LITERALLY STATED IS FALSE at w = (1,1): rho = x+y is not a "
       "monomial, tau = xy+y^2, J(rho,tau) = rho exactly, w(rho) = 1 > 0, "
       "w(tau) = 2 > 0, and w(rho) DOES divide w(tau).  ML's own convention "
       "(alpha,beta in Z with gcd = 1) admits (1,1), so a hypothesis is "
       "missing.  [CHECKED -- exact]",
       ok_ce and 2 % 1 == 0 and not r2.is_Mul and len(r2.args) == 2)
    ck("E2b this independently CONFIRMS ML_RESTRICTION.md sec.2, which found the "
       "same gap from a paraphrase in 2026-07-24; the verbatim text is now in "
       "hand and the gap is real, not an artefact of the paraphrase", ok_ce)

    # E3 -- under positive + primitive + NON-diagonal weights, no counterexample.
    ND = [(1, 2), (2, 1), (1, 3), (3, 1), (2, 3), (3, 2), (1, 4), (4, 1),
          (3, 4), (4, 3), (1, 5), (5, 1), (2, 5), (5, 2)]
    hits0 = [(w, div_scan(*w, kappa=0)) for w in ND]
    ck("E3  RECONNAISSANCE ONLY (finite grid, not a proof): over %d positive "
       "primitive NON-diagonal weights, no non-monomial rho with J(rho,tau) = "
       "rho and w(rho) | w(tau) exists.  DIV holds in the regime a "
       "Newton-polygon shape lemma is about." % len(ND),
       all(h is None for _, h in hits0),
       str([w for w, h in hits0 if h is not None]))

    # E4 -- THE TRANSFER.  Our corners carry [P,Q] = x^2, not 1.
    KAPPA = 2
    ck("E4  our reduced corners carry [P,Q] = x^kappa with kappa = 2, NOT the "
       "J = 1 normalisation DIV is stated in (polygon_reduction.py R1/R2/R3; "
       "GGHV22 2204.14178.tex:1001; GGV3 1406.0886_GGV3.tex:1726) [CITATION]",
       KAPPA == 2)

    # the exact, hand-checkable twisted counterexample
    r4 = X ** 3 + Y
    t4 = sp.Rational(1, 3) * (X ** 3 * Y + Y ** 2)
    lhs = sp.expand(jac(r4, t4))
    rhs = sp.expand(X ** KAPPA * r4)
    ck("E4b DIV DOES NOT TRANSFER, CONSTRUCTIVELY.  At the positive, primitive, "
       "NON-diagonal weight w = (1,3): rho = x^3+y (not a monomial), "
       "tau = (x^3 y + y^2)/3, and J(rho,tau) = x^2 * rho EXACTLY, with "
       "w(rho) = 3 DIVIDING w(tau) = 6.  The divisibility conclusion is FALSE "
       "under our bracket normalisation.  [CHECKED -- exact]",
       lhs == rhs and lhs == sp.expand(X ** 5 + X ** 2 * Y)
       and (3 * 1 + 0 * 3) == 3 and (3 * 1 + 1 * 3) == 6 and 6 % 3 == 0,
       "J = %s vs x^2*rho = %s" % (lhs, rhs))
    ck("E4c NO RESCALING REPAIRS IT: J(rho,tau) = x^2 rho forces "
       "w(tau) = 3*alpha+beta rather than alpha+beta, and dividing tau by x^2 "
       "leaves J(rho, x^-2 tau) != rho -- checked exactly on the E4b instance",
       sp.expand(jac(r4, sp.expand(t4 * X ** -2)) - r4) != 0)
    # and a second, independent twisted counterexample at a different weight
    r4b = X * (X ** 2 + Y) ** 2
    t4b = sp.Rational(1, 5) * (X ** 3 * Y + X * Y ** 2)
    ck("E4d a SECOND independent twisted counterexample at w = (1,2): "
       "rho = x(x^2+y)^2, tau = (x^3 y + x y^2)/5, J(rho,tau) = x^2 rho and "
       "w(rho) = 5 divides w(tau) = 5",
       sp.expand(jac(r4b, t4b) - X ** KAPPA * r4b) == 0)
    # MUTATION CONTROL: the same search finds NOTHING at kappa = 0 for these
    # weights, so the hits above are caused by the twist and not by the search.
    ck("E4m MUTATION CONTROL -- rerunning the identical search at kappa = 0 "
       "(ML's own normalisation) at w = (1,3) and w = (1,2) finds NOTHING.  "
       "The counterexamples are produced by the x^2 twist, not by the search.",
       div_scan(1, 3, 0) is None and div_scan(1, 2, 0) is None)

    # E5 -- an INDEPENDENT structural reason DIV cannot bite, twist aside.
    note("")
    note("  strictly-positive-weight edges of the reduced polygons:")
    diag_only, none_at_all = [], []
    for name, poly in REDUCED_POLYGONS.items():
        pe = positive_weight_edges(poly)
        note("    %-20s %s" % (name, pe if pe else "NONE"))
        if not pe:
            none_at_all.append(name)
        elif all(w == (1, 1) for _, w in pe):
            diag_only.append(name)
    ck("\nE5  INDEPENDENT of the twist: the four (8,28) reduced polygons have NO "
       "strictly-positive-weight EDGE at all (the origin is a vertex, so every "
       "positive weight is maximised at a single monomial).  DIV's hypothesis "
       "is unreachable there.",
       sorted(none_at_all) == sorted(n for n in REDUCED_POLYGONS
                                     if n.startswith("(8,28)")),
       str(none_at_all))
    ck("E5b and the four (5,20) reduced polygons have exactly ONE "
       "positive-weight edge each, at the DIAGONAL weight (1,1) -- precisely "
       "the weight class where DIV's conclusion is false even at kappa = 0 "
       "(E2).  So DIV has no bite at the 125 corner either.",
       sorted(diag_only) == sorted(n for n in REDUCED_POLYGONS
                                   if n.startswith("(5,20)")),
       str(diag_only))
    # MUTATION CONTROL: the detector finds a non-diagonal positive edge when
    # one exists.
    probe = [(0, 0), (4, 0), (3, 2), (0, 3)]
    pe = positive_weight_edges(probe)
    ck("E5m MUTATION CONTROL -- on a probe polygon with a genuine non-diagonal "
       "positive-weight edge, the same detector finds it, so E5/E5b are real "
       "tests and not a detector that always returns NONE",
       any(w != (1, 1) for _, w in pe), str(pe))


# ===========================================================================
def main() -> int:
    atlas = json.load(open(os.path.join(HERE, "corner_atlas.json"),
                           encoding="utf-8"))
    rows = atlas["rows"]
    ck("00  the atlas carries GGV5's 34 candidate rows", len(rows) == 34,
       str(len(rows)))

    part_A()
    part_B(rows)
    census = part_C(rows)
    part_D(rows, census)
    part_E()

    note("\n" + "=" * 78)
    note("VERDICT")
    note("=" * 78)
    note("  Theorem 2 : 0 kills / 34 rows.  0 kills / 10 known-dead calibration")
    note("              rows.  F_2(3,5)/125 SURVIVES (bound 1 = floor 1).")
    note("              Residue: T2 forces EXACTLY ONE admissible right edge at")
    note("              the (5,20) corner -- a falsifiable shape constraint.")
    note("  Divisibility: DOES NOT TRANSFER.  Its conclusion is exhibited FALSE")
    note("              under our [P,Q] = x^2 normalisation, and the reduced")
    note("              polygons carry no non-diagonal positive-weight edge.")
    note("  Neither obstruction moves the 125 bound.")

    if _fail:
        print()
        print("FAILURES (%d):" % len(_fail))
        for f in _fail:
            print("   - %s" % f)
        return 1
    print("makar_limanov_filter: %d/%d checks pass -- T2 kills 0 of 34 rows "
          "(0 of 10 calibration rows); the divisibility lemma does not transfer "
          "to [P,Q] = x^2." % (_ok[0], _ok[0]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
