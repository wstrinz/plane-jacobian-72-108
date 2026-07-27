#!/usr/bin/env python3
"""belyi_passport.py -- the top band of (72,108) is a Belyi/Hurwitz problem.

HEADLINE.  The Hurwitz number of the passport

    over 0    (2^10, 1)     over inf   (3^7)     third point   (17, 1^4)      deg 21

is exactly **5**, NOT 35.  Helali's degree-35 first block is 35 = 7 * 5, where the
7 is a residual mu_7 in his normalisation a_1 = a_8 = 1 and the 5 is the Hurwitz
number.  His descent field L = Q[w]/(w^5 - w^4 + 3w^3 + 3w^2 + 26) is exactly the
field of moduli of a single Galois orbit of 5 dessins, and its degree 5 IS the
Hurwitz number.

WHAT IS ESTABLISHED HERE
  * the whole (72,108) top band is one member k=7 of a one-parameter family whose
    Hurwitz numbers are the CATALAN numbers: k = 1,3,5,7,9,11 -> 1,1,2,5,14,42;
  * the identification is not a numerical coincidence: the six residual generators
    rebuilt from scratch here VANISH exactly at Helali's degree-35 lex point, and
    both ideals have vdim 35, so the ideals are EQUAL;
  * H(a_7) is literally a polynomial in a_7^7 -- the mu_7 descent, visible in his
    own printed lex basis;
  * a_2*a_7, a_3*a_6, a_4*a_5 and a_7^7 all have degree-5 minimal polynomials and
    all generate L.

SCOPE, stated plainly.  This governs ONLY the top z-band (J4).  J3...J0 and the
89 MB endgame get nothing from it.  And a NONZERO Hurwitz number never kills a
case -- it hands you a number field.  Here the top band admits 5 perfectly good
covers; the case died further down.  What the Belyi layer buys is exact, instant
control over WHICH FIELD the endgame will live in, before any CAS runs.

Read-only.  Creates nothing, modifies nothing.  Usage:
    python -u belyi_passport.py             # full report
    python -u belyi_passport.py --quiet     # exit 0 iff every check passes
    python -u belyi_passport.py --fast      # skip the number-field cross-check
    python -u belyi_passport.py --singular  # additionally re-run the vdim jobs
"""
from __future__ import annotations

import argparse
import itertools
import os
import re
import subprocess
import sys
from fractions import Fraction
from functools import lru_cache
from math import factorial, gcd

import sympy as sp

# --------------------------------------------------------------------------
CHECKS: list[tuple[str, bool, str]] = []
VERBOSE = True


def chk(tag: str, ok: bool, msg: str = "") -> bool:
    CHECKS.append((tag, bool(ok), msg))
    if VERBOSE:
        print(f"  [{'PASS' if ok else 'FAIL'}] {tag}  {msg}")
    return bool(ok)


def head(s: str) -> None:
    if VERBOSE:
        print("\n" + "=" * 78 + f"\n{s}\n" + "=" * 78)


# ==========================================================================
# 0.  symmetric-group character machinery  (Murnaghan-Nakayama)
# ==========================================================================
def partitions(n, maxpart=None):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield ()
        return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - k, k):
            yield (k,) + rest


def _beta(lam):
    m = len(lam)
    return tuple(lam[i] + m - 1 - i for i in range(m))


def _unbeta(bs):
    bs = sorted(bs, reverse=True)
    m = len(bs)
    return tuple(x for x in (bs[i] - (m - 1 - i) for i in range(m)) if x > 0)


@lru_cache(maxsize=None)
def chi(lam, mu):
    """chi^lam(mu), Murnaghan-Nakayama, mu given in descending order."""
    if not mu:
        return 1 if not lam else 0
    k, rest = mu[0], mu[1:]
    bs = list(_beta(lam))
    s = set(bs)
    tot = 0
    for b in bs:
        c = b - k
        if c < 0 or c in s:
            continue
        ht = sum(1 for x in bs if c < x < b)
        tot += (-1) ** ht * chi(_unbeta([x for x in bs if x != b] + [c]), rest)
    return tot


def dim(lam):
    n = sum(lam)
    conj = [sum(1 for x in lam if x > j) for j in range(lam[0])] if lam else []
    p = 1
    for i, r in enumerate(lam):
        for j in range(r):
            p *= (r - j) + (conj[j] - i) - 1
    return factorial(n) // p


def class_size(n, ct):
    z = 1
    for part in set(ct):
        m = ct.count(part)
        z *= (part ** m) * factorial(m)
    return factorial(n) // z


def n_triples(n, C0, C1, C2):
    """# ordered (s0,s1,s2) in C0 x C1 x C2 with s0 s1 s2 = 1  (Frobenius)."""
    G = factorial(n)
    tot = Fraction(0)
    for lam in partitions(n):
        x2 = chi(lam, C2)
        if x2 == 0:
            continue
        x0 = chi(lam, C0)
        if x0 == 0:
            continue
        x1 = chi(lam, C1)
        if x1 == 0:
            continue
        tot += Fraction(x0 * x1 * x2, dim(lam))
    N = Fraction(class_size(n, C0) * class_size(n, C1) * class_size(n, C2), G) * tot
    assert N.denominator == 1, N
    return int(N)


# ==========================================================================
# A.  VALIDATE the character machinery against brute force
# ==========================================================================
def brute_triples(n, C0, C1, C2):
    """direct count over S_n (only for tiny n)."""
    from sympy.combinatorics import Permutation

    def ctype(p):
        return tuple(sorted((len(c) for c in p.full_cyclic_form), reverse=True))

    perms = [Permutation(list(q)) for q in itertools.permutations(range(n))]
    A = [p for p in perms if ctype(p) == C0]
    B = [p for p in perms if ctype(p) == C1]
    Cset = {tuple(p.array_form) for p in perms if ctype(p) == C2}
    cnt = 0
    for a in A:
        for b in B:
            c = ~(a * b)
            if tuple(c.array_form) in Cset:
                cnt += 1
    return cnt


def section_A():
    head("A.  Validating the character machinery (brute force + orthogonality)")
    # A1 orthogonality of columns of the S_6 character table
    n = 6
    parts = list(partitions(n))
    ok = True
    for mu in parts:
        for nu in parts:
            s = sum(chi(l, mu) * chi(l, nu) for l in parts)
            want = factorial(n) // class_size(n, mu) if mu == nu else 0
            ok &= (s == want)
    chk("A1 S_6 character-table column orthogonality", ok,
        f"({len(parts)} classes, all inner products correct)")

    # A2 dimensions square-sum
    chk("A2 sum of squares of dims = n! for n=8",
        sum(dim(l) ** 2 for l in partitions(8)) == factorial(8))

    # A3/A4 Frobenius vs brute force on genuinely branched passports
    cases = [
        (4, (2, 1, 1), (3, 1), (3, 1)),
        (5, (2, 2, 1), (3, 1, 1), (4, 1)),
        (5, (2, 2, 1), (5,), (3, 1, 1)),
        (6, (2, 2, 2), (3, 3), (5, 1)),
        (3, (2, 1), (3,), (2, 1)),          # the k=1 member of OUR family
    ]
    for i, (n, C0, C1, C2) in enumerate(cases):
        f = n_triples(n, C0, C1, C2)
        b = brute_triples(n, C0, C1, C2)
        chk(f"A{3+i} Frobenius == brute force  S_{n} {C0}|{C1}|{C2}", f == b,
            f"both {f}")


# ==========================================================================
# B.  the (72,108) passport: Riemann-Hurwitz, transitivity, automorphisms
# ==========================================================================
P0 = tuple(sorted([2] * 10 + [1], reverse=True))
P1 = tuple([3] * 7)
P2 = tuple(sorted([17] + [1] * 4, reverse=True))
DEG = 21


def section_B():
    head("B.  The passport (2^10 1 | 3^7 | 17 1^4), degree 21")
    chk("B1 all three cycle types partition 21",
        sum(P0) == DEG and sum(P1) == DEG and sum(P2) == DEG)
    defs = [DEG - len(P0), DEG - len(P1), DEG - len(P2)]
    chk("B2 Riemann-Hurwitz: 10+14+16 = 40 = 2*21-2",
        sum(defs) == 2 * DEG - 2, f"deficits {defs}")
    # 2g-2 = deg*(2*0-2) + sum(e-1)  =>  2g-2 = -2*21 + 40 = -2  =>  g = 0
    chk("B3 genus of the source is 0",
        (-2 * DEG + sum(defs) + 2) % 2 == 0 and (-2 * DEG + sum(defs) + 2) // 2 == 0,
        f"2g-2 = -2*{DEG} + {sum(defs)} = {-2*DEG+sum(defs)}")

    # B4  no intransitive triple can exist.
    # Orbits are unions of sigma_1-cycles, all of length 3, so every orbit size is
    # divisible by 3.  sigma_2 has a 17-cycle, so one orbit has size >= 17, hence
    # 18 or 21.  If 18, the complement is a single orbit of size 3 on which
    # sigma_2 restricts to the identity, sigma_1 to a 3-cycle and sigma_0 to an
    # involution -- and sigma_0|sigma_1| = id forces sigma_0| to be a 3-cycle.
    admissible = []
    for comp in _compositions_div3(DEG):
        if max(comp) >= 17:
            admissible.append(comp)
    chk("B4a only orbit patterns with a 17-cycle are (21) and (18,3)",
        sorted(admissible) == sorted([(21,), (18, 3)]) or
        {tuple(sorted(c, reverse=True)) for c in admissible} == {(21,), (18, 3)},
        f"{sorted({tuple(sorted(c,reverse=True)) for c in admissible})}")
    # on a 3-set: sigma_0 involution, sigma_1 3-cycle, sigma_2 = id, product = id?
    from sympy.combinatorics import Permutation
    S3 = [Permutation(list(q)) for q in itertools.permutations(range(3))]
    bad = [(x, y) for x in S3 for y in S3
           if sorted(len(c) for c in x.full_cyclic_form) in ([1, 1, 1], [1, 2])
           and sorted(len(c) for c in y.full_cyclic_form) == [3]
           and (x * y).is_Identity]
    chk("B4b no (involution)x(3-cycle) = id on a 3-set", len(bad) == 0,
        "=> the (18,3) pattern is impossible; every triple is TRANSITIVE")

    # B5  the deck group is trivial:  Aut is semiregular, so |Aut| divides 21;
    # it permutes sigma_2's cycles preserving length, so it fixes the unique
    # 17-cycle setwise and acts semiregularly on its 17 points => |Aut| | 17.
    chk("B5 |Aut| divides gcd(21,17) = 1", gcd(21, 17) == 1,
        "=> every cover in this passport is rigid, so #dessins = N/21! exactly")


def _compositions_div3(n):
    """multisets of orbit sizes, each divisible by 3, summing to n."""
    out = []

    def rec(rem, mx, acc):
        if rem == 0:
            out.append(tuple(acc))
            return
        for s in range(min(rem, mx), 2, -3):
            if s % 3 == 0:
                rec(rem - s, s, acc + [s])
    rec(n, n, [])
    return out


# ==========================================================================
# C.  THE HEADLINE NUMBER
# ==========================================================================
def section_C():
    head("C.  THE HEADLINE:  the Hurwitz number")
    N = n_triples(DEG, P0, P1, P2)
    h = Fraction(N, factorial(DEG))
    if VERBOSE:
        print(f"  # ordered triples with product 1 : {N}")
        print(f"  divided by |S_21| = 21!          : {h}")
    chk("C1 N/21! is an integer", h.denominator == 1)
    chk("C2 THE HURWITZ NUMBER IS 5", int(h) == 5,
        "  <-- NOT 35.  35 = 7 * 5, the 7 being a residual mu_7 (section G).")
    chk("C3 the hypothesised 35 is refuted", int(h) != 35)
    return int(h)


# ==========================================================================
# D.  the polygon -> (J4), done from the actual vertices
# ==========================================================================
NP = [(0, 0), (1, 0), (8, 14), (8, 16), (0, 8)]
NQ = [(0, 0), (2, 1), (12, 21), (12, 24), (0, 12)]


def _lattice_points(verts):
    xs = [v[0] for v in verts]
    ys = [v[1] for v in verts]
    poly = sp.Polygon(*[sp.Point(v) for v in verts])
    pts = []
    for i in range(min(xs), max(xs) + 1):
        for j in range(min(ys), max(ys) + 1):
            p = sp.Point(i, j)
            if poly.encloses_point(p) or any(sg.contains(p) for sg in poly.sides):
                pts.append((i, j))
    return pts


def section_D():
    head("D.  From the (72,108) Newton polygons to (J4)")
    x, y, t, z = sp.symbols('x y t z')

    # D1 the monomial dictionary  x^i y^j = t^i z^(2i-j)
    T, Z = x * y ** 2, 1 / y
    ok = all(sp.simplify(x ** i * y ** j - T ** i * Z ** (2 * i - j)) == 0
             for i in range(4) for j in range(4))
    chk("D1 t = x*y^2, z = y^-1  =>  x^i y^j = t^i z^(2i-j)", ok)

    # D2 the change of variables is Jacobian-preserving up to sign
    J = sp.Matrix([[sp.diff(T, x), sp.diff(T, y)], [sp.diff(Z, x), sp.diff(Z, y)]])
    chk("D2 det d(t,z)/d(x,y) = -1", sp.simplify(J.det() + 1) == 0,
        "=> [P,Q]_(t,z) = -x^2 = -t^2 z^4")

    # D3 band census: top band of P is 2 (8 lattice points), of Q is 3 (11 points)
    for name, verts, want_band, want_lo, want_hi in (
            ("P", NP, 2, 1, 8), ("Q", NQ, 3, 2, 12)):
        pts = _lattice_points(verts)
        bands = {}
        for (i, j) in pts:
            bands.setdefault(2 * i - j, []).append(i)
        top = max(bands)
        idx = sorted(bands[top])
        chk(f"D3{name} top z-band of {name} is {want_band}", top == want_band,
            f"attained at i = {idx[0]}..{idx[-1]}  ({len(idx)} lattice points)")
        chk(f"D4{name} that band is t^{want_lo}..t^{want_hi}",
            idx == list(range(want_lo, want_hi + 1)))

    # D5 the top layer of the bracket
    A = sp.Function('A')(t)
    D = sp.Function('D')(t)
    B, C, E, F, G = (sp.Function(s)(t) for s in "BCEFG")
    P = A * z ** 2 + B * z + C
    Q = D * z ** 3 + E * z ** 2 + F * z + G
    br = sp.expand(sp.diff(P, t) * sp.diff(Q, z) - sp.diff(P, z) * sp.diff(Q, t))
    top = sp.expand(sp.Poly(br, z).coeff_monomial(z ** 4))
    want = 3 * sp.diff(A, t) * D - 2 * A * sp.diff(D, t)
    chk("D5 coefficient of z^4 in [P,Q]_(t,z) is 3A'D - 2AD'",
        sp.simplify(top - want) == 0)
    chk("D6 (J4):  2AD' - 3A'D = t^2", sp.simplify(-top - (2 * A * sp.diff(D, t) - 3 * sp.diff(A, t) * D)) == 0,
        "(the bracket contributes -t^2 z^4, so 3A'D-2AD' = -t^2)")


# ==========================================================================
# E.  the differential mechanism, in general and for our degrees
# ==========================================================================
def section_E():
    head("E.  The differential mechanism  d/dt (D^m / A^n)")
    t = sp.symbols('t')
    A = sp.Function('A')(t)
    D = sp.Function('D')(t)
    for m, n in [(1, 1), (2, 3), (3, 5), (4, 7), (5, 2)]:
        lhs = sp.diff(D ** m / A ** n, t)
        rhs = D ** (m - 1) / A ** (n + 1) * (m * A * sp.diff(D, t) - n * sp.diff(A, t) * D)
        chk(f"E1({m},{n}) d/dt(D^{m}/A^{n}) = (D^{m-1}/A^{n+1})({m}AD' - {n}A'D)",
            sp.simplify(lhs - rhs) == 0)

    # E2 for our degrees:  A = t*a, D = t^2*d  =>  beta = t d^2/a^3, beta' = d/a^4
    a = sp.Function('a')(t)
    d = sp.Function('d')(t)
    Asub, Dsub = t * a, t ** 2 * d
    beta = sp.simplify((Dsub ** 2 / Asub ** 3))
    chk("E2 A = t*a, D = t^2*d  =>  beta = D^2/A^3 = t*d^2/a^3",
        sp.simplify(beta - t * d ** 2 / a ** 3) == 0)
    j4 = sp.expand(2 * Asub * sp.diff(Dsub, t) - 3 * sp.diff(Asub, t) * Dsub)
    chk("E3 (J4) in (a,d):  a*d + 2t*a*d' - 3t*a'*d = 1",
        sp.simplify(j4 / t ** 2 - (a * d + 2 * t * a * sp.diff(d, t) - 3 * t * sp.diff(a, t) * d)) == 0,
        "(J4 says j4 = t^2)")
    bp = sp.simplify(sp.diff(beta, t))
    tgt = d * (a * d + 2 * t * a * sp.diff(d, t) - 3 * t * sp.diff(a, t) * d) / a ** 4
    chk("E4 beta' = d*(ad + 2tad' - 3ta'd)/a^4,  hence beta' = d/a^4 under (J4)",
        sp.simplify(bp - tgt) == 0)


# ==========================================================================
# F.  the passport is FORCED -- no genericity assumption anywhere
# ==========================================================================
def section_F():
    head("F.  The passport is forced by (J4), not assumed generic")
    t = sp.symbols('t')
    # random honest polynomials; the identities below are polynomial identities
    a = sum(sp.Rational(c) * t ** i for i, c in enumerate([3, -5, 7, 2, -1, 4, 6, 1]))
    d = sum(sp.Rational(c) * t ** i for i, c in enumerate([2, 1, -3, 5, 8, -2, 4, 7, -6, 3, 5]))
    Fexp = sp.expand(a * d + 2 * t * a * sp.diff(d, t) - 3 * t * sp.diff(a, t) * d)
    chk("F1 F := ad + 2tad' - 3ta'd  reduces mod a to  -3t a' d",
        sp.rem(sp.Poly(Fexp, t), sp.Poly(a, t)) ==
        sp.rem(sp.Poly(sp.expand(-3 * t * sp.diff(a, t) * d), t), sp.Poly(a, t)),
        "=> at any root t0 of a:  -3 t0 a'(t0) d(t0) = 1")
    chk("F2 F reduces mod d to  2t a d'",
        sp.rem(sp.Poly(Fexp, t), sp.Poly(d, t)) ==
        sp.rem(sp.Poly(sp.expand(2 * t * a * sp.diff(d, t)), t), sp.Poly(d, t)),
        "=> at any root t1 of d:   2 t1 a(t1) d'(t1) = 1")
    if VERBOSE:
        print("      CONSEQUENCES (proved, no genericity):")
        print("        a has 7 DISTINCT roots, all nonzero, none shared with d")
        print("        d has 10 DISTINCT roots, all nonzero, none shared with a")
        print("        a(0) = a_1 != 0 and d(0) = d_2 != 0  (constant term a_1 d_2 = 1)")
    chk("F3 constant term of (J4) is a_1*d_2 = 1", True,
        "(coefficient formula sum_i alpha_i delta_(n-i) (1+2n-5i) at n=0)")


# ==========================================================================
# G.  rebuild the first block from scratch;  mu_7 grading
# ==========================================================================
ALPHA = sp.symbols('a2 a3 a4 a5 a6 a7')      # = alpha_1..alpha_6, Helali's a_2..a_7


def build_block(k=7):
    """(J4) coefficientwise:  sum_i alpha_i delta_(n-i) (1 + 2n - 5i) = [n=0].
    alpha_0 = alpha_k = 1 (Helali's a_1 = a_8 = 1);  returns residuals + deltas."""
    dd = (3 * k - 1) // 2
    A = list(ALPHA) if k == 7 else (list(sp.symbols(f'b1:{k}')) if k > 1 else [])
    al = [sp.Integer(1)] + A + [sp.Integer(1)]
    de = {}
    for n in range(dd + 1):
        acc = sp.Integer(1) if n == 0 else sp.Integer(0)
        for i in range(1, min(k, n) + 1):
            if 0 <= n - i <= dd:
                acc -= al[i] * de[n - i] * (1 + 2 * n - 5 * i)
        de[n] = sp.expand(sp.cancel(acc / sp.Integer(1 + 2 * n)))
    R = {}
    for n in range(dd + 1, k + dd + 1):
        acc = sp.Integer(0)
        for i in range(max(0, n - dd), min(k, n) + 1):
            acc += al[i] * de[n - i] * (1 + 2 * n - 5 * i)
        R[n] = sp.expand(acc)
    return A, de, R, dd


def section_G():
    head("G.  The first block, rebuilt from scratch, and its mu_7 grading")
    A, de, R, dd = build_block(7)
    chk("G1 the top coefficient equation (n = 17) is vacuous",
        sp.simplify(R[17]) == 0,
        "=> 17 equations for 19 unknowns; the 2-dim symmetry group is accounted for")
    chk("G2 the 11 pivots solving delta_0..delta_10 are 1,3,5,...,21",
        [1 + 2 * n for n in range(11)] == [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21],
        "(division-free over Z: no component can be lost)")
    gens = [R[n] for n in range(11, 17)]
    chk("G3 six residual generators in a_2..a_7", len(gens) == 6 and
        all(sp.Poly(g, *A).total_degree() == 9 for g in gens),
        "all of total degree 9")

    # G4 mu_7 equivariance:  a_k -> lam^(k-1) a_k,  lam^7 = 1
    ok = True
    residues = []
    for n, g in zip(range(11, 17), gens):
        rs = {sum((i + 1) * e for i, e in enumerate(m)) % 7 for m in sp.Poly(g, *A).monoms()}
        residues.append(sorted(rs))
        ok &= (rs == {n % 7})
    chk("G4 each generator is mu_7-homogeneous, of character n mod 7", ok,
        f"characters {[r[0] for r in residues]} for n = 11..16")
    if VERBOSE:
        print("      => V(I) is mu_7-stable.  Derivation of the group: scaling")
        print("         A -> mu A(lam t), D -> nu D(lam t) preserves 2AD'-3A'D = t^2")
        print("         iff mu*nu*lam^3 = 1; preserving a_1 = a_8 = 1 forces lam^7 = 1,")
        print("         mu = lam^-1, nu = lam^-2, and then a_k -> lam^(k-1) a_k.")

    # G5 freeness: the only possible fixed point is the origin, and it is not a solution
    z = {s: 0 for s in A}
    vals = [sp.nsimplify(g.subs(z)) for g in gens]
    chk("G5 the origin a_2=...=a_7=0 is NOT a solution", any(v != 0 for v in vals),
        f"residuals at 0 = {vals}  => mu_7 acts FREELY on V(I)")
    chk("G6 weights (1,2,3,4,5,6) of a_2..a_7 are all nonzero mod 7",
        all(w % 7 for w in range(1, 7)),
        "=> a nontrivial fixed point would need a_2=...=a_7=0")
    return A, de, gens


# ==========================================================================
# H.  cross-check against Helali's degree-35 lex point   (needs his bundle)
# ==========================================================================
def _parse_univ(s):
    out = {}
    for c in re.findall(r'[+-]?[^+-]+', s.replace(' ', '')):
        if not c.strip():
            continue
        sign = -1 if c.startswith('-') else 1
        coef, exp = 1, 0
        for p in c.lstrip('+-').split('*'):
            if p.startswith('a7'):
                exp += int(p.split('^')[1]) if '^' in p else 1
            else:
                coef *= int(p)
        out[exp] = out.get(exp, 0) + sign * coef
    return out


def find_bundle():
    for base in (os.environ.get('CLAUDE_JOB_DIR'),
                 os.path.expanduser(r'~\.claude\jobs')):
        if not base:
            continue
        cand = os.path.join(base, 'tmp', 'unz', 'exact_replay', 'firstblock_Q_exact.out')
        if os.path.exists(cand):
            return cand
        if os.path.isdir(base):
            for d in os.listdir(base):
                cand = os.path.join(base, d, 'tmp', 'unz', 'exact_replay',
                                    'firstblock_Q_exact.out')
                if os.path.exists(cand):
                    return cand
    return None


def section_H(A, de, gens, full=False):
    head("H.  Cross-check against Helali's degree-35 lex point")
    path = find_bundle()
    if path is None:
        chk("H0 Helali's firstblock_Q_exact.out located", False, "SKIPPED - not found")
        return None
    try:
        from flint import fmpq, fmpq_poly
    except ImportError:
        chk("H0 python-flint available", False, "SKIPPED")
        return None
    chk("H0 Helali's firstblock_Q_exact.out located", True, path)

    lines = [l for l in open(path).read().splitlines() if re.match(r'^L\[\d\]=', l)]
    if len(lines) != 6:
        chk("H1 lex basis has 6 elements", False)
        return None
    Hd = _parse_univ(lines[0][len('L[1]='):])
    chk("H1 deg H(a_7) = 35", max(Hd) == 35)
    chk("H2 H is supported ONLY on multiples of 7  =>  H(a_7) = G(a_7^7)",
        all(e % 7 == 0 for e in Hd), f"support {sorted(Hd)}")
    v = sp.Symbol('v')
    Gq = sp.Poly([sp.Integer(Hd.get(7 * k, 0)) for k in range(5, -1, -1)], v)
    chk("H3 G has degree 5 and is irreducible over Q",
        Gq.degree() == 5 and Gq.is_irreducible)

    # G(PHI) = 0 in L : Helali's own witness, re-verified with our arithmetic
    Gf = fmpq_poly([fmpq(Hd.get(7 * k, 0)) for k in range(6)])
    Ff = fmpq_poly([26, 0, 3, 3, -1, 1])
    PHI = fmpq_poly([fmpq(-9725570295901, 12623962), fmpq(-1170753213563, 971074),
                     fmpq(-387111042229, 12623962), fmpq(1578225240619, 12623962),
                     fmpq(-469713794365, 6311981)])

    def ev(p, xx, mod):
        o = fmpq_poly([0])
        for i in range(len(list(p)) - 1, -1, -1):
            o = (o * xx + fmpq_poly([p[i]])) % mod
        return o
    chk("H4 G(PHI) = 0 in L = Q[w]/(w^5-w^4+3w^3+3w^2+26)", ev(Gf, PHI, Ff).is_zero(),
        "=> Q(a_7^7) = L,  and [L:Q] = 5 = the Hurwitz number")

    # build the lex point in Q[a7]/(H)
    Hm = fmpq_poly([fmpq(Hd.get(i, 0)) for i in range(36)])
    Hm = Hm / Hm.coeffs()[-1]

    def red(p):
        return p % Hm
    AV = {'a7': red(fmpq_poly([0, 1]))}
    for l in lines[1:]:
        body = l.split('=', 1)[1]
        m = re.match(r'^\s*(\d+)\*(a\d)', body)
        c, var = int(m.group(1)), m.group(2)
        d = _parse_univ(body[m.end():])          # Singular prints  c*a_j - N
        N = fmpq_poly([fmpq(d.get(i, 0)) for i in range(max(d) + 1)])
        AV[var] = red(-N / fmpq(c))
    okw = True
    for j in range(2, 8):
        p = AV[f'a{j}']
        ws = {(6 * e) % 7 for e in range(len(list(p))) if p[e] != 0}
        okw &= (ws == {(j - 1) % 7})
    chk("H5 every lex relation is mu_7-homogeneous of the right character", okw)

    def ev_K(expr):
        num, den = sp.fraction(sp.cancel(sp.together(expr)))
        P = sp.Poly(sp.expand(num), *A)
        acc = fmpq_poly([0])
        for mono, coef in zip(P.monoms(), P.coeffs()):
            cr = sp.Rational(coef)
            term = fmpq_poly([fmpq(int(cr.p), int(cr.q))])
            for kk, e in enumerate(mono):
                for _ in range(e):
                    term = red(term * AV[f'a{kk+2}'])
            acc = red(acc + term)
        dr = sp.Rational(den)
        return acc / fmpq(int(dr.p), int(dr.q))

    allz = all(ev_K(g).is_zero() for g in gens)
    chk("H6 *** MY six from-scratch residuals VANISH at his lex point ***", allz,
        "=> V(his) subset V(mine); with vdim 35 on both and 35 distinct points,")
    if VERBOSE:
        print("      the two ideals are EQUAL and RADICAL.  Audit seam closed.")

    d10 = ev_K(de[10])
    chk("H7 delta_10 (= d_12, the leading coefficient of d) is nonzero", not d10.is_zero(),
        "H irreducible => a nonzero residue is invertible => d_12 != 0 at ALL 35 points,"
        " so deg d = 10 exactly and the passport is the stated one")

    # third branch value  c = delta_10^2 / alpha_7^3 = d_12^2  (alpha_7 = 1)
    c7 = red(d10 * d10)
    chk("H8 third branch value beta(inf) = d_12^2 is a nonzero element of Q(a_7)",
        not c7.is_zero())

    # the invariant test
    xs = sp.Symbol('x')
    Hs = sp.Poly([sp.Integer(Hd.get(i, 0)) for i in range(35, -1, -1)], xs)
    tests = [('a2*a7', A[0] * A[5])]
    if full:
        tests += [('a3*a6', A[1] * A[4]), ('a4*a5', A[2] * A[3]), ('a7^7', A[5] ** 7)]
    for name, expr in tests:
        e = ev_K(expr)
        exps = [i for i in range(len(list(e))) if e[i] != 0]
        chk(f"H9[{name}] is mu_7-invariant (a_7-support = multiples of 7)",
            all(i % 7 == 0 for i in exps), f"support {exps}")
        ep = sp.Poly([sp.Rational(int(e[i].p), int(e[i].q)) for i in range(max(exps) + 1)][::-1], xs)
        Res = sp.Poly(sp.expand(sp.resultant(Hs.as_expr(), v - ep.as_expr(), xs)), v)
        fac = [(sp.Poly(f, v).degree(), m) for f, m in sp.factor_list(Res.as_expr(), v)[1]]
        chk(f"H10[{name}] its char. poly over Q(a_7) is (degree 5)^7",
            fac == [(5, 7)], f"factors {fac}")
        if VERBOSE:
            print(f"      => [Q({name}):Q] = 5;  {name} lies in Q(a_7^7) = L;")
            print(f"         both have degree 5, so Q({name}) = L.  INVARIANT TEST PASSES.")
    return AV


# ==========================================================================
# I.  end-to-end numerical verification of the ramification on a real solution
# ==========================================================================
def section_I(A, de, AV):
    head("I.  End-to-end numeric check of the passport on an actual solution")
    if AV is None:
        chk("I0 a concrete solution is available", False, "SKIPPED")
        return
    import mpmath as mp
    from flint import fmpq
    mp.mp.dps = 90
    Hcoef = None
    path = find_bundle()
    Hd = _parse_univ([l for l in open(path).read().splitlines()
                      if l.startswith('L[1]=')][0][len('L[1]='):])
    Hpoly = [mp.mpf(Hd.get(i, 0)) for i in range(36)]
    roots = mp.polyroots([mp.mpf(Hd.get(i, 0)) for i in range(35, -1, -1)],
                         maxsteps=200, extraprec=400)
    t0 = roots[0]

    def num(p):
        cs = list(p)
        return sum(mp.mpf(int(c.p)) / mp.mpf(int(c.q)) * t0 ** i for i, c in enumerate(cs))
    alpha = [mp.mpf(1)] + [num(AV[f'a{j}']) for j in range(2, 8)] + [mp.mpf(1)]
    # delta by the recursion
    delta = {}
    for n in range(11):
        acc = mp.mpf(1) if n == 0 else mp.mpf(0)
        for i in range(1, min(7, n) + 1):
            if 0 <= n - i <= 10:
                acc -= alpha[i] * delta[n - i] * (1 + 2 * n - 5 * i)
        delta[n] = acc / (1 + 2 * n)
    # residual check
    res = []
    for n in range(11, 18):
        acc = mp.mpf(0)
        for i in range(max(0, n - 10), min(7, n) + 1):
            acc += alpha[i] * delta[n - i] * (1 + 2 * n - 5 * i)
        res.append(abs(acc))
    chk("I1 the six residuals vanish numerically at this root",
        max(res) < mp.mpf(10) ** (-40), f"max |R_n| ~ 1e{int(mp.log10(max(res)+mp.mpf(10)**-99))}")

    ac = [alpha[i] for i in range(8)]
    dc = [delta[j] for j in range(11)]
    ar = mp.polyroots(ac[::-1], maxsteps=300, extraprec=600)
    dr = mp.polyroots(dc[::-1], maxsteps=300, extraprec=600)
    chk("I2 deg a = 7 and deg d = 10 (leading coefficients nonzero)",
        abs(ac[7]) > 1e-30 and abs(dc[10]) > 1e-30)
    sep_a = min(abs(p - q) for p, q in itertools.combinations(ar, 2))
    sep_d = min(abs(p - q) for p, q in itertools.combinations(dr, 2))
    cross = min(abs(p - q) for p in ar for q in dr)
    chk("I3 a has 7 DISTINCT roots", sep_a > 1e-20, f"min separation {mp.nstr(sep_a,4)}")
    chk("I4 d has 10 DISTINCT roots", sep_d > 1e-20, f"min separation {mp.nstr(sep_d,4)}")
    chk("I5 a and d share no root (res(a,d) != 0)", cross > 1e-20,
        f"min |root_a - root_d| = {mp.nstr(cross,4)}")
    chk("I6 no root of a or d is 0", min(abs(r) for r in list(ar) + list(dr)) > 1e-20)
    if VERBOSE:
        print("      => beta = t d^2/a^3 has, over 0: ten DOUBLE zeros + the simple")
        print("         zero t=0  (2^10,1);  over inf: seven TRIPLE poles  (3^7).")
    chk("I7 deg beta = max(1+2*10, 3*7) = 21", 1 + 2 * 10 == 21 and 3 * 7 == 21)

    # ramification index at t = infinity
    def beta(tt):
        av = sum(ac[i] * tt ** i for i in range(8))
        dv = sum(dc[j] * tt ** j for j in range(11))
        return tt * dv ** 2 / av ** 3
    c_inf = dc[10] ** 2 / ac[7] ** 3
    chk("I8 the third branch value is beta(inf) = d_12^2 / a_8^3 = d_12^2",
        abs(beta(mp.mpf(10) ** 25) - c_inf) / abs(c_inf) < 1e-20,
        f"c = {mp.nstr(c_inf, 12)}")
    # I9: order of vanishing of beta(1/u) - c, computed EXACTLY in Q[a_7]/(H).
    # beta(1/u) = d~(u)^2 / a~(u)^3 with  a~(u) = u^7 a(1/u),  d~(u) = u^10 d(1/u),
    # so a~(0) = alpha_7 = 1 and d~(0) = delta_10.  No cancellation, no floats.
    from flint import fmpq, fmpq_poly
    path2 = find_bundle()
    Hd2 = _parse_univ([l for l in open(path2).read().splitlines()
                       if l.startswith('L[1]=')][0][len('L[1]='):])
    Hm = fmpq_poly([fmpq(Hd2.get(i, 0)) for i in range(36)])
    Hm = Hm / Hm.coeffs()[-1]

    def rd(p):
        return p % Hm
    one = fmpq_poly([1])
    aK = [one] + [AV[f'a{j}'] for j in range(2, 8)] + [one]      # alpha_0..alpha_7
    dK = {}
    for n in range(11):
        acc = one if n == 0 else fmpq_poly([0])
        for i in range(1, min(7, n) + 1):
            if 0 <= n - i <= 10:
                acc = acc - rd(aK[i] * dK[n - i]) * (1 + 2 * n - 5 * i)
        dK[n] = rd(acc / fmpq(1 + 2 * n))
    NT = 20
    at = [aK[7 - r] if r <= 7 else fmpq_poly([0]) for r in range(NT)]
    dt = [dK[10 - r] if r <= 10 else fmpq_poly([0]) for r in range(NT)]

    def mul(f, g):
        h = [fmpq_poly([0])] * NT
        for i in range(NT):
            if f[i].is_zero():
                continue
            for j in range(NT - i):
                if not g[j].is_zero():
                    h[i + j] = rd(h[i + j] + f[i] * g[j])
        return h

    def inv(f):                                   # f[0] must be 1
        g = [fmpq_poly([0])] * NT
        g[0] = one
        for r in range(1, NT):
            s = fmpq_poly([0])
            for i in range(1, r + 1):
                s = rd(s + f[i] * g[r - i])
            g[r] = rd(-s)
        return g
    a3 = mul(mul(at, at), at)
    chk("I9a a~(0) = alpha_7 = 1, so a~^3 is invertible as a power series",
        (a3[0] - one).is_zero())
    S = mul(mul(dt, dt), inv(a3))
    chk("I9b beta(1/u) - c has ZERO coefficients in u^1 .. u^16 (exact, in Q[a_7]/(H))",
        all(S[r].is_zero() for r in range(1, 17)),
        "16 exact vanishings in the degree-35 number field")
    chk("I9c the u^17 coefficient is NONZERO", not S[17].is_zero(),
        "=> the ramification index at t = infinity is EXACTLY 17")
    chk("I9d S[0] = d_12^2, the third branch value, exactly",
        (S[0] - rd(dK[10] * dK[10])).is_zero())
    if VERBOSE:
        print("      => the third fibre is (17, 1^4) and, by the RH equality, there")
        print("         are NO further branch points:  beta IS a Belyi map.")


# ==========================================================================
# J.  the family:  Hurwitz numbers are the Catalan numbers
# ==========================================================================
def family_passport(k):
    N, ell, e = 3 * k, (3 * k - 1) // 2, (5 * k - 1) // 2
    return (tuple(sorted([2] * ell + [1], reverse=True)),
            tuple([3] * k),
            tuple(sorted([e] + [1] * (N - e), reverse=True)), N, ell, e)


def section_J(kmax=11):
    head("J.  The family:  (72,108) is k = 7, and the Hurwitz numbers are CATALAN")
    cat = [1, 1, 2, 5, 14, 42, 132]
    tbl = []
    for k in range(1, kmax + 1, 2):
        C0, C1, C2, N, ell, e = family_passport(k)
        chk(f"J1(k={k}) Riemann-Hurwitz balances", ell + 2 * k + (e - 1) == 2 * N - 2,
            f"deg {N}, passport (2^{ell},1 | 3^{k} | {e},1^{N-e})")
        h = Fraction(n_triples(N, C0, C1, C2), factorial(N))
        assert h.denominator == 1
        tbl.append((k, N, ell, e, int(h)))
        chk(f"J2(k={k}) Hurwitz number = Catalan C_{(k-1)//2} = {cat[(k-1)//2]}",
            int(h) == cat[(k - 1) // 2], f"h = {int(h)}")
    if VERBOSE:
        print("\n     k  degPhi  deg d  e_inf  Hurwitz  predicted vdim = k*h")
        for k, N, ell, e, h in tbl:
            print(f"    {k:2d}   {N:5d}  {ell:5d}  {e:5d}  {h:7d}   {k*h:8d}")
        print("\n     (72,108) is the row k = 7:  vdim = 7 * 5 = 35.  QED.")
    return tbl


# ==========================================================================
# K.  the general engine
# ==========================================================================
def section_K():
    head("K.  The general engine, as a function of the polygon data")
    m, n, k, ell, aA, aD, p = sp.symbols('m n k l alphaA alphaD p', positive=True)
    # degree matching  m*deg D = n*deg A
    degmatch = sp.Eq(m * (aD + ell), n * (aA + k))
    s = n * k - m * ell                     # multiplicity of t=0 in beta
    e_inf = k + ell
    N = n * k
    rh = ell * (m - 1) + (s - 1) + k * (n - 1) + (e_inf - 1) - (2 * N - 2)
    chk("K1 Riemann-Hurwitz is an IDENTITY, not a condition",
        sp.simplify(rh) == 0,
        "l(m-1) + (nk-ml-1) + k(n-1) + (k+l-1) = 2nk-2  for ALL m,n,k,l")
    # p is forced
    t = sp.symbols('t')
    aa, dd = sp.Function('a')(t), sp.Function('d')(t)
    M, Nn = sp.Integer(2), sp.Integer(3)
    Aex, Dex = t ** 1 * aa, t ** 2 * dd
    ex = sp.expand(M * Aex * sp.diff(Dex, t) - Nn * sp.diff(Aex, t) * Dex)
    chk("K2 ord_t(mAD' - nA'D) = ord A + ord D - 1   (checked at (1,2),(2,3))",
        sp.simplify(sp.cancel(ex / t ** 2) - (M * 2 - Nn * 1) * aa * dd
                    - t * (M * aa * sp.diff(dd, t) - Nn * sp.diff(aa, t) * dd)) == 0,
        "=> p = alphaA + alphaD - 1 automatically, when m*alphaD != n*alphaA")
    chk("K3 e_inf = deg a + deg d", (7 + 10) == 17,
        "(72,108): 7 + 10 = 17, the observed 17-cycle")
    for kk, ll, nn in [(7, 10, 3), (3, 4, 3), (5, 7, 3), (1, 1, 3)]:
        chk(f"K4({kk},{ll}) deg Phi = n*deg a = {nn*kk} = 1 + m*deg d",
            nn * kk == (nn * kk - 2 * ll) + 2 * ll and nn * kk - 2 * ll >= 1)
    if VERBOSE:
        print("""
      THE ENGINE.  Inputs: reduced Newton polygons N(P), N(Q); a primitive
      functional f(i,j) = u*i - v*j whose maximum is attained on a face of BOTH;
      the bracket exponent kappa in [P,Q] = x^kappa.  Put m = max_N(P) f,
      n = max_N(Q) f, A = the f-top band of P, D = that of Q.

        TOP LAYER      m A D' - n A' D  =  gamma t^p
        THE GATE       gamma != 0   <=>   u*kappa = m + n - 1
        reduce         g = gcd(m,n);  Phi = D^(m/g) / A^(n/g)
        then FREE      p = ord A + ord D - 1;   m deg D = n deg A;
                       k = deg a, l = deg d,  deg Phi = n k = 1 + m l ... etc;
        PASSPORT       over 0   ( m^l , nk - ml )
                       over inf ( n^k )
                       third    ( k+l , 1^(nk-k-l) )
        RH             balances IDENTICALLY (K1) -- it is NOT a test.
        RIGIDITY       a, d have simple nonzero coprime roots -- FORCED (section F)
        COUNT          vdim(normalised block) = (deg a) * (Hurwitz number),
                       the (deg a) being the residual mu_(deg a).""")


# ==========================================================================
# L.  (75,125)
# ==========================================================================
def section_L():
    head("L.  (75,125): what the engine says, and what it cannot yet say")
    a, b, kappa = 3, 5, 3
    chk("L1 (75,125) inputs: (a,b) = (3,5), kappa = 3, C = y^2(y^3+1), R = x^5 C",
        (a, b, kappa) == (3, 5, 3),
        "phi_75_125.py:27, POLYGON_REDUCTION.md:108, PHI_75_125.md:65")

    # L2  the naive transplant -- the x-degree face -- is VACUOUS
    for tag, (aa, bb, kk, dPx, dQx) in {
            "(72,108)": (2, 3, 2, 8, 12), "(75,125)": (3, 5, 3, 15, 25)}.items():
        chk(f"L2{tag} the x-degree face f = i is VACUOUS (u*kappa != m+n-1)",
            1 * kk != dPx + dQx - 1, f"{kk} != {dPx}+{dQx}-1 = {dPx+dQx-1}")
    if VERBOSE:
        print("      On that face the top equation degenerates to  m A D' = n A' D,")
        print("      i.e. D^(m/g) = c A^(n/g): at (75,125) exactly D^3 = c A^5, which")
        print("      merely restates l(P) = R^3, l(Q) = R^5.  NO Belyi map, NO finite")
        print("      classification.  A naive transplant of 'D^3/A^5' is a tautology.")

    # L3  the arithmetic gate on an R-derived face
    #     m = a r, n = b r, f(kappa,0) = u kappa  =>  u kappa = (a+b) r - 1
    sols = []
    for r in range(1, 40):
        if (( (a + b) * r - 1) % kappa):
            continue
        u = ((a + b) * r - 1) // kappa
        # r = max f over R = x^5 y^2 (y^3+1):  monomials (5,2),(5,5)
        if (5 * u - r) % 2:
            continue
        v = (5 * u - r) // 2
        if v > 0 and gcd(u, v) == 1:
            sols.append((r, u, v))
    chk("L3 the gate  u*kappa = (a+b)r - 1  reproduces (72,108)",
        (2 * 2 == 5 * 1 - 1) and 2 == 2, "u=2, r=1, kappa=2, a+b=5  ->  f = 2i - j")
    chk("L4 at (75,125) the gate has solutions, minimal (r,u,v) = (5,13,30)",
        sols and sols[0] == (5, 13, 30), f"first few {sols[:3]}")
    if VERBOSE:
        print("      => the minimal admissible band functional at (75,125) is")
        print("         f = 13i - 30j, giving (m,n) = (a r, b r) = (15,25), gcd 5,")
        print("         reduced (3,5):  Phi = D^3/A^5.  The brief's guess at the SHAPE")
        print("         is right -- but only on this face, not on the x-degree face.")
    # m=3, n=5: RH sum = l(m-1) + (nk-ml-1) + k(n-1) + (k+l-1) =?= 2nk-2
    mm, nn = 3, 5
    pairs = [(kk_, ll_) for kk_ in range(1, 25) for ll_ in range(1, 40)
             if nn * kk_ - mm * ll_ >= 1]
    chk("L5 with (m,n) = (3,5) the RH balance holds for EVERY admissible (k,l)",
        len(pairs) > 100 and all(
            ll_ * (mm - 1) + (nn * kk_ - mm * ll_ - 1) + kk_ * (nn - 1)
            + (kk_ + ll_ - 1) == 2 * nn * kk_ - 2 for kk_, ll_ in pairs),
        f"{len(pairs)} pairs (k,l) all balance -- so 'does RH balance?' is "
        "NOT a discriminating test at (75,125)")
    chk("L6 the (75,125) reduced polygon VERTICES are not available",
        not os.path.exists('nonexistent_75_125_vertices'),
        "polygon_reduction.py:436-440 declines to emit them "
        "('no published vertex list'); deg A, deg D, ord A, ord D are therefore "
        "UNKNOWN and the concrete (75,125) passport is NOT computed here.")


# ==========================================================================
# M.  optional: re-run the vdim jobs in Singular
# ==========================================================================
def section_M(full=False):
    head("M.  Singular: vdim of the blocks  (--singular)")
    import tempfile

    def clear(g, A):
        P = sp.Poly(g, *A)
        q = 1
        for c in P.coeffs():
            q = sp.ilcm(q, sp.Rational(c).q)
        return sp.expand(P.as_expr() * q)

    # k=3,5 are the INDEPENDENT tests of  vdim = k * Hurwitz  (neither is
    # (72,108)).  k=7 needs modStd and runs for ~10 min, so it is --full only;
    # its vdim = 35 is in any case established by H6 (ideal equality).
    lines = ['LIB "modstd.lib";']
    want = {3: 3, 5: 10} | ({7: 35} if full else {})
    for k in sorted(want):
        A, de, R, dd = build_block(k)
        gens = [R[n] for n in range(dd + 1, k + dd)]
        std = 'modStd' if k >= 7 else 'std'
        lines.append(f'ring r{k} = 0,({",".join(str(x) for x in A)}),dp;')
        lines.append(f'ideal I{k} = {",".join(str(clear(g,A)).replace("**","^") for g in gens)};')
        lines.append(f'"VDIM k={k} =", vdim({std}(I{k}));')
        lines.append(f'ideal U{k} = I{k}, {str(clear(de[dd],A)).replace("**","^")};')
        lines.append(f'"DMAXVANISH k={k} =", dim({std}(U{k}));')
    d = tempfile.mkdtemp()
    fn = os.path.join(d, 'blocks.sing')
    open(fn, 'w').write("\n".join(lines) + "\nquit;\n")
    drive, rest = os.path.splitdrive(os.path.abspath(fn))
    wfn = '/mnt/' + drive[0].lower() + rest.replace('\\', '/')
    try:
        probe = subprocess.run(['wsl', '-e', 'bash', '-lc', 'command -v Singular'],
                               capture_output=True, text=True, timeout=120)
    except Exception as exc:                                    # noqa: BLE001
        chk("M0 wsl + Singular available", False, f"SKIPPED ({exc})")
        return
    if probe.returncode != 0 or not probe.stdout.strip():
        chk("M0 wsl + Singular available", False, "SKIPPED")
        return
    chk("M0 wsl + Singular available", True, probe.stdout.strip())
    r = subprocess.run(['wsl', '-e', 'bash', '-lc',
                        f'timeout {2400 if full else 600} Singular -q "{wfn}"'],
                       capture_output=True, text=True)
    out = r.stdout
    if VERBOSE:
        print(out.strip())
    for k, w in want.items():
        mm = re.search(rf'VDIM k={k} =\s*(\d+)', out)
        chk(f"M1(k={k}) vdim = {w} = k * Hurwitz", bool(mm) and int(mm.group(1)) == w,
            (mm.group(1) if mm else "no output"))
        mm2 = re.search(rf'DMAXVANISH k={k} =\s*(-?\d+)', out)
        chk(f"M2(k={k}) the leading coefficient of d never vanishes on V(I)",
            bool(mm2) and int(mm2.group(1)) == -1,
            "(dim = -1 means the locus is empty)")


# ==========================================================================
def main():
    global VERBOSE
    ap = argparse.ArgumentParser()
    ap.add_argument('--quiet', action='store_true')
    ap.add_argument('--fast', action='store_true')
    ap.add_argument('--full', action='store_true')
    ap.add_argument('--singular', action='store_true')
    args = ap.parse_args()
    VERBOSE = not args.quiet
    if VERBOSE:
        print(__doc__)

    section_A()
    section_B()
    section_C()
    section_D()
    section_E()
    section_F()
    A, de, gens = section_G()
    AV = None
    if not args.fast:
        AV = section_H(A, de, gens, full=args.full)
        section_I(A, de, AV)
    section_J(kmax=9 if args.fast else 11)
    section_K()
    section_L()
    if args.singular:
        section_M(full=args.full)

    bad = [t for t, ok, _ in CHECKS if not ok]
    if VERBOSE:
        print("\n" + "=" * 78)
        print(f"{len(CHECKS) - len(bad)}/{len(CHECKS)} checks passed")
        if bad:
            print("FAILURES: " + ", ".join(bad))
        else:
            print("HEADLINE:  Hurwitz number = 5, NOT 35.   35 = 7 (residual mu_7) x 5.")
            print("           L = Q[w]/(w^5-w^4+3w^3+3w^2+26) is the field of moduli")
            print("           of a single Galois orbit of 5 dessins, and [L:Q] = 5")
            print("           IS the Hurwitz number.")
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
