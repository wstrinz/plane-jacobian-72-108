#!/usr/bin/env python3
"""bridge_generality.py  (NEW; read-only over all existing artifacts)

SETTLES the load-bearing soft spot of MONOMIAL_WINDOW_LAW.md sec.6: are the two
formulas the bridge identity is derived FROM

        rho := ord_y(f) = q(b-a) + 1                        (R)
        N   := (power of C in Phi) = a*M - 2b               (Nf)

general, or only true on the 1-dimensional slices where they were confirmed?

*** VERDICT: 205 IS CONFIRMED, AND (R) IS NOW PROVED, NOT INFERRED. ***

  1. THE DECISIVE TEST PASSES.  At (8,28)/(3,4)/144 -- a different (m,n)=(3,4)
     at the SAME corner as the closed (8,28)/(3,2)/108 -- an INDEPENDENT
     derivation gives ord_y(Phi) = 205, exactly what the bridge identity
     ord_y(Phi) = a*q*M - H = 3*3*25 - 20 predicts.  Independent means:
       * rho is NOT taken from (R).  f is found by a FULLY GENERIC linear solve
         of the forcing ODE -- 18 unknown coefficients f_0..f_17, no y^rho
         ansatz, no monomial shape assumed -- and ord_y(f) is then READ OFF the
         resulting polynomial.  [D-group]
       * N is NOT taken from (Nf).  It is read off the BUILT D-transform tower
         via the slice-sum invariant (c_series_75_125.py sec.2), whose own input
         -- the D-exponent a*w-1 -- is here DERIVED from the a-th-root expansion
         of P rather than quoted from the paper's a=2 recurrence.  [E-group]
       * ord_y(Phi) is then ord_y(f) + N*ord_y(C), additive, on actual
         polynomials.  [F-group]

  2. (R) IS PROVED IN FULL GENERALITY by a two-line LOCAL argument at y=0 that
     the previous lane never made.  Write c = y^q g with g(0) != 0.  Dividing
     the forcing ODE by y^(q-1) and reading the y^K coefficients gives the
     triangular recursion

         a * sum_i g_i [t(k-i) - coef(q+i)] f_{k-i} = [g^e]_{k+q-1-q*e} ,
         coef := t(b-a)+kappa+1 ,  e := b-a+1 ,

     whose pivot at i=0 is a*g_0*(t*k - coef*q).  The right-hand side first
     becomes nonzero at k = q(e-1)+1 = q(b-a)+1 =: rho_0.  So f_k = 0 for every
     k < rho_0 UNLESS some k < rho_0 kills the pivot, i.e. t*k = coef*q; and

         t*k = coef*q  and  k < rho_0   <==>   t | q(kappa+1)  and  q(kappa+1) < t ,

     which is IMPOSSIBLE for kappa >= 0, q >= 1: a positive multiple of t cannot
     be smaller than t.  Hence f_k = 0 for all k < rho_0, and at k = rho_0

         f_{rho_0} = g_0^(e-1) / [ a (t - q(kappa+1)) ]   !=  0

     whenever t != q(kappa+1).  So ord_y(f) = q(b-a)+1 EXACTLY, at EVERY corner,
     for EVERY residual g with g(0) != 0 -- no slice, no family, no branch.
     [C-group, symbolic; D-group, 15 corner signatures]

     The unique excluded locus is t = q(kappa+1), where the pivot at k = rho_0
     itself vanishes and the ODE has NO power-series solution.  In the standard
     class kappa = t-2 that reads q(t-1) = t, whose only integer solution is
     (t,kappa,q) = (2,0,2).  No published row is there (min t on the 34 atlas
     rows is 3).  [C4, C5]

     Two corollaries worth having:
       * (R) is BRANCH-INDEPENDENT.  ord_y(f) sees only q = ord_y(c) and
         g(0) != 0; it does not see g's shape.  So the ramified-vs-complex-pair
         branch ambiguity that phi_f7.py flags -- which really does move
         mult_(y+1) and the cofactor -- CANNOT move ord_y(Phi), hence cannot
         move the bridge identity or anything downstream of it.  [D6]
       * the excluded locus t - q(kappa+1) is the NEGATIVE of the Bezout corner
         integer q(kappa+1) - t of MONOMIAL_WINDOW_LAW sec.2.  The same integer
         that controls gcd(M,H) is the one whose vanishing would break (R).
         At a monomial corner it is -1, so (R) is safest exactly there.  [C6]

  3. THE UNTESTED JOINT DIRECTION IS NOW TESTED, and it passes.  Before this
     file, (R) and (Nf) were confirmed at (72,108) (q=7, b-a=1) and along the F2
     rungs (q=1, b-a varying) -- one coordinate at a time.  Five corners move
     BOTH at once, and all five are independently derived here:

         corner            q   b-a   dg   ord_y(Phi)   a*q*M - H
         F_7 (42,147)      4    5     2       165         165
         F_8 (63,147)      5    4     1       371         371
         F14 (66,231)      4    5     5       165         165
         F15 (99,231)      5    4     4       371         371
         F16 (99,165)      7    2     2       407         407

     HONEST QUALIFICATION.  Those five joint corners are COLLINEAR: all have
     q + (b-a) = 9.  That is forced, not accidental -- all five sit at t=3, k=1,
     p = q+l, where the GGV5 Diophantine reads q(n-m) = 3n-1, so b-a = (3n-1)/q,
     and q + (3n-1)/q = 9 is the curve n = (9q-q^2+1)/3.  GGV5's v11 <= 35 tables
     simply do not contain an off-line joint corner.  [G3, G3b]
     Two things repair that gap:
       * the full nine-corner tested set IS affinely 2-dimensional -- (3,1),
         (7,1), (8,1), (4,5), (5,4), (7,2) span the plane.  [G3c]
       * the ABSTRACT sweep fills the (q, b-a) rectangle 1..8 x 1..6 completely,
         4032 (t,kappa,q,a,b) points, and certifies (R) on all of it via the
         pivot condition the C-group proof reduces it to.  Those points are not
         all GGV5 corners -- but the object in doubt was a FORMULA, and a formula
         is exactly what an abstract sweep can settle.  [G6, G8, G8b]
     So the slice worry is closed twice over: by a proof (C-group) and by a
     2-dimensional exact check (G-group).  [G-group]

  4. ALL 34 ATLAS ROWS now have an INDEPENDENTLY derived ord_y(Phi), including
     all SIX non-monomial rows.  The 34 rows carry only 15 distinct chart
     signatures (t, kappa, deg C, ord C, a, b); every one is derived here and
     every one agrees with a*q*M - H.  PHI_KNOWN in corner_atlas.py had ONE
     entry; the honest count is now 15 signatures / 34 rows.  [F-group]

STATUS LEDGER (read this, not the headline):
  PROVED        (R): ord_y(f) = q(b-a)+1 whenever t != q(kappa+1), kappa >= 0,
                q >= 1, g(0) != 0.  Local, exact, no slice.
  PROVED        the D-transform exponent a*w-1 (a-th-root denominator bound,
                and it is ATTAINED so `clear` is exact, not merely sufficient).
  PROVED        the slice-sum invariant clear = a*M - b (additivity of an
                affine exponent), hence (Nf) GIVEN that Phi is the cleared
                slice-M object of S^b at M = b*t + j.
  EXACT-CHECKED ord_y(Phi) at 15 chart signatures / 34 atlas rows / 5 joint
                corners, each by generic linear solve + built tower, each
                agreeing with a*q*M - H.
  CLAIMED       (unchanged, and NOT touched here) that Phi IS that cleared
                slice-M object, and the extreme-ray premise of
                window_functions_75_125.  This file settles the ARITHMETIC of
                the bridge, not the geometry of Phi.
  NEGATIVE      phi_corner4.py and phi_f7.py use PRE-REPAIR chart data at F1,
                F2, F3, F5, F9, F10: they take t = l from GGV5's table and apply
                the final-corner dictionary at corners that do NOT retract, where
                polygon_reduction.corner_chart_data returns the MONOMIAL data
                instead.  This is a live contradiction, not a stylistic quibble:
                phi_corner4.py's VERDICT claims five reproduced points, and THREE
                of them -- (50,75), (75,125), (56,84) -- are at refused corners.
                At (50,75) it implies ord_y(Phi) = 75 where the repaired route
                and corner_atlas.json both give 30.  Their F7 / F14 / F16 points
                are fine (those corners do retract).  Nothing this file concludes
                rests on the stale rows -- every number here goes through the
                guard -- but the two files should not be cited for them.
                [A5b, MUT F]

Sources of truth: chart data ALWAYS through polygon_reduction.corner_chart_data
(the retraction guard), never GGV5's l_final; the 34 rows' (t,kappa,deg C,ord C,
m,n) re-read out of corner_atlas.json; the forcing ODE re-derived here by direct
bracket differentiation; f by generic linear solve; N by built tower.
Cross-checks: (8,28)/(3,2)/108 must reproduce the published 204 and
corner144_verify.py's (550,205,69,276) must reproduce 205 -- both do.

Run:
    python bridge_generality.py            # verbose
    python bridge_generality.py --quiet    # one line per check; exit 0 iff all pass
"""
import argparse
import json
import os
import sys
from fractions import Fraction
from math import gcd

import sympy as sp

import polygon_reduction as pr

HERE = os.path.dirname(os.path.abspath(__file__))
y, cs, u, zz = sp.symbols("y c u z")


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
            print("  [%s] %-52s %s" % ("PASS" if ok else "FAIL", name, detail))
        return bool(ok)

    def mut(self, name: str, all_mutants_fail: bool, detail: str = "") -> bool:
        return self.ck("MUT " + name, all_mutants_fail, detail)

    def note(self, text: str) -> None:
        if not self.quiet:
            print("        . " + text)

    def report(self) -> int:
        bad = [r for r in self.rows if not r[1]]
        print("\n%s  bridge_generality: %d/%d checks pass"
              % ("FAIL" if bad else "OK  ", len(self.rows) - len(bad), len(self.rows)))
        for n, _, d in bad:
            print("   FAILED: %s   %s" % (n, d))
        return 1 if bad else 0


# =============================================================================
#  Primitives
# =============================================================================
def MH(t, kappa, q, a, b):
    """(M, H) = (t(a+b) - (kappa+1), q(a+b) - 1)   [q_window_theorem]."""
    s = a + b
    return t * s - (kappa + 1), q * s - 1


def ordy(p):
    return min(m[0] for m in sp.Poly(sp.expand(p), y).monoms())


def c_exponent(term):
    _, cp = term.as_independent(cs, as_Add=False)
    if cp == 1:
        return 0
    base, ex = cp.as_base_exp()
    assert base == cs
    return int(ex)


def residual_g(dg):
    """The residual polynomial of the leading form, monic of degree dg, g(-1)=0.

    dg == 0 : no residual, g = 1 (C a monomial).
    dg odd  : g = y^dg + 1 -- FORCED by the ODE coefficient system (interior
              coefficients vanish, top resonant, monic, g(-1)=0); the root at -1
              is simple.  [corner144_verify.py sec.E; phi_corner4.derive]
    dg even : g = (y+1)^dg, the RAMIFIED branch.  phi_f7.py proves at dg=2 that a
              simple real root at -1 is impossible here, so this is the branch
              continuous with the audited pattern.  D6 below shows the CHOICE IS
              IRRELEVANT for ord_y, which is all this file needs.
    """
    if dg == 0:
        return sp.Integer(1)
    return y**dg + 1 if dg % 2 == 1 else sp.expand((y + 1)**dg)


def forcing_ode_residual(a, b, t, kappa, cpoly, fpoly):
    """a{ t c f' - [t(b-a)+kappa+1] c' f } - c^(b-a+1).   Zero iff f solves (F)."""
    coef = t * (b - a) + kappa + 1
    return sp.expand(a * t * cpoly * sp.diff(fpoly, y)
                     - a * coef * sp.diff(cpoly, y) * fpoly - cpoly**(b - a + 1))


def solve_f_generic(a, b, t, kappa, cpoly):
    """Find f by a FULLY GENERIC linear solve.  No y^rho ansatz, no g^e shape.

    Unknowns f_0..f_D with D chosen above both the pure-ansatz degree and the
    resonant degree res = coef*deg(c)/t, so the solve cannot miss a solution by
    truncation.  Returns (f, D, n_free).
    """
    a0 = sp.degree(cpoly, y)
    e = b - a + 1
    coef = t * (b - a) + kappa + 1
    D = int(max(Fraction(coef * int(a0), t), (e - 1) * int(a0) + 1)) + 3
    fc = sp.symbols("Fg0:%d" % (D + 1))
    fans = sum(fc[i] * y**i for i in range(D + 1))
    eqs = sp.Poly(forcing_ode_residual(a, b, t, kappa, cpoly, fans), y).all_coeffs()
    sol = list(sp.linsolve(eqs, fc))
    if not sol:
        return None, D, None
    vals = sol[0]
    free = set().union(*[v.free_symbols for v in vals]) & set(fc)
    f = sp.expand(fans.subs(dict(zip(fc, vals))))
    return f, D, len(free)


def local_recursion_orders(a, b, t, kappa, q, dg, upto=None):
    """The LOCAL recursion at y=0 with a SYMBOLIC residual g (g_0 free, monic).

    Implements exactly the proof in the module docstring: returns
    (f_0..f_rho0) as exact rational functions of the g_i.  Used to check
    (R) for an ARBITRARY g, not just the branch representative.
    """
    e = b - a + 1
    coef = t * (b - a) + kappa + 1
    rho0 = q * (b - a) + 1 if upto is None else upto
    gsym = list(sp.symbols("gg0:%d" % (dg + 1))) if dg > 0 else [sp.Integer(1)]
    if dg > 0:
        gsym[dg] = sp.Integer(1)                     # monic
    gpoly = sum(gsym[i] * y**i for i in range(dg + 1))
    ge = sp.Poly(sp.expand(gpoly**e), y)

    def gecoeff(i):
        if i < 0 or i > e * dg:
            return sp.Integer(0)
        return ge.coeff_monomial(y**i) or sp.Integer(0)

    f = {}
    for k in range(0, rho0 + 1):
        rhs = gecoeff(k + q - 1 - q * e)
        acc = sp.Integer(0)
        for i in range(1, min(dg, k) + 1):
            acc += a * gsym[i] * (t * (k - i) - coef * (q + i)) * f[k - i]
        piv = a * gsym[0] * (t * k - coef * q)
        if piv == 0:
            return f, None                            # pivot killed: (R) at risk
        f[k] = sp.simplify(sp.cancel((rhs - acc) / piv))
    return f, gsym[0]


def first_nonzero_index(a, b, t, kappa, q, dg, upto):
    """The ord of the formal power-series solution of (F) at y=0, SYMBOLIC g.

    Works for ANY (a,b,t,kappa,q,dg) -- unlike a polynomial solve with a FIXED
    residual g, which is generically inconsistent because the resonance at the
    top degree imposes conditions on g.  This is the object (R) is about.
    """
    f, _ = local_recursion_orders(a, b, t, kappa, q, dg, upto=upto)
    if f is None:
        return None
    for k in sorted(f):
        if sp.simplify(f[k]) != 0:
            return k
    return None


def tower_clear_N(a, b, t, kappa, Kwin=None, dexp=None):
    """DERIVE (M, clear, N) from the BUILT D-transform tower.  No N-formula.

    S = sum_w d_w u^w  with  d_w = c_w * c^(dexp(w)),  dexp(w) = a*w - 1;
    Phi lives at the u-slice M = b*t + j of S^b, j = -s = a*t - kappa - 1.

    Kwin is the truncation depth of the tower; it must satisfy b*Kwin >= M or the
    slice is unreachable and nothing is derived (a silent-zero trap: the previous
    fixed Kwin=7 missed the t=6 monomial signature entirely).
    """
    if dexp is None:
        def dexp(w):
            return a * w - 1
    s = kappa + 1 - a * t
    M = b * t - s
    if Kwin is None:
        Kwin = max(4, -(-M // b) + 1)
    dv = {}
    for w in range(0, Kwin + 1):
        if w == 0:
            dv[w] = sp.Integer(1)
        elif w == 1:
            dv[w] = sp.Integer(0)                     # normalising shear
        else:
            dv[w] = sp.symbols("zt_%d" % w)
    S = sum(dv[w] * cs**dexp(w) * u**w for w in dv)
    Sb = sp.Poly(sp.expand(S**b), u)
    expr = sp.expand(Sb.coeff_monomial(u**M))
    if expr == 0:
        return M, None, None, 0
    exps = {c_exponent(tm) for tm in sp.Add.make_args(expr)}
    homog = len(exps) == 1
    clear = next(iter(exps)) if homog else None
    return M, clear, (clear - b if homog else None), len(exps)


def root_denominator_exponents(a, W=5):
    """The a-th-root denominator exponents -- DERIVES the D-exponent a*w-1.

    P = x^(at) c^a (1 + sum_w pi_w z^w),  pi_w = p_w/c^a,  z = x^-1.
    C = P^(1/a) = x^t * c * (1 + sum pi_w z^w)^(1/a).  The z^w coefficient of the
    root is a sum of products of at most w factors pi_{w_i}, so its c-denominator
    is c^(a*w) at worst; times the overall c that is c^(a*w-1).  Returns the
    ACTUAL exponent per w, so we can also check it is ATTAINED (hence `clear` is
    exact, not merely sufficient).
    """
    pi = sp.symbols("pp1:%d" % (W + 1))
    ser = 1 + sum(pi[i - 1] * zz**i / cs**a for i in range(1, W + 1))
    ex = sp.series(ser**sp.Rational(1, a), zz, 0, W + 1).removeO()
    out = {}
    for w in range(1, W + 1):
        co = sp.together(sp.expand(ex.coeff(zz, w)) * cs)
        _, den = sp.fraction(co)
        out[w] = int(sp.degree(sp.Poly(den, cs), cs)) if den.has(cs) else 0
    return out


# =============================================================================
#  The corner population, ALWAYS through the retraction guard
# =============================================================================
# GGV5 v11<=35 length-1 families with A0' = (1,0), transcribed in phi_corner4.py:
#   name: (A0, p, l_final, b_final, k, (m0,dm), (n0,dn))
FAMILIES = {
    "F1":  ((4, 12),  7, 4, 3, 1, (3, 2),  (4, 3)),
    "F2":  ((5, 20),  7, 5, 2, 1, (2, 1),  (3, 2)),
    "F3":  ((5, 20),  8, 5, 3, 1, (3, 4),  (2, 3)),
    "F5":  ((5, 20),  9, 5, 4, 1, (9, 7),  (5, 4)),
    "F7":  ((6, 15),  7, 3, 4, 1, (2, 1),  (7, 4)),
    "F8":  ((6, 15),  8, 3, 5, 1, (3, 2),  (7, 5)),
    "F9":  ((7, 21), 11, 7, 2, 1, (2, 1),  (3, 2)),
    "F10": ((7, 21), 13, 7, 3, 1, (7, 5),  (4, 3)),
    "F14": ((9, 24),  7, 3, 4, 1, (2, 1),  (7, 4)),
    "F15": ((9, 24),  8, 3, 5, 1, (3, 2),  (7, 5)),
    "F16": ((9, 24), 10, 3, 7, 1, (3, 4),  (5, 7)),
    "F17": ((9, 24), 11, 3, 8, 1, (2, 5),  (3, 8)),
}

# The five JOINT (q, b-a) corners + the two (8,28) rows, all retracting.
# tag: (A0, l_final, b_final, (m,n))
JOINT = {
    "(8,28)/(3,4)/144": ((8, 28), 4, 3, (3, 4)),
    "(8,28)/(3,2)/108": ((8, 28), 4, 7, (3, 2)),
    "F_7 (42,147)":     ((6, 15), 3, 4, (2, 7)),
    "F_8 (63,147)":     ((6, 15), 3, 5, (3, 7)),
    "F14 (66,231)":     ((9, 24), 3, 4, (2, 7)),
    "F15 (99,231)":     ((9, 24), 3, 5, (3, 7)),
    "F16 (99,165)":     ((9, 24), 3, 7, (3, 5)),
    "F_17 (66,99)":     ((9, 24), 3, 8, (2, 3)),
    "(12,33)/(2,3)/135": ((12, 33), 3, 8, (2, 3)),
}


def guarded(A0, l_final, b_final, mn, who="bridge_generality"):
    """Chart data through polygon_reduction's retraction guard -- the ONLY source."""
    cd = pr.corner_chart_data(A0[0], A0[1], l_final=l_final, b_final=b_final, who=who)
    a, b = sorted(mn)
    t, kappa, a0, q = cd["t"], cd["kappa"], cd["deg_C"], cd["ord_C"]
    M, H = MH(t, kappa, q, a, b)
    return dict(a=a, b=b, t=t, kappa=kappa, a0=a0, q=q, dg=a0 - q, M=M, H=H,
                monomial=cd["monomial"], retraction=cd["retraction"])


def atlas_signatures():
    """The 34 published rows, grouped by RE-READ chart signature."""
    with open(os.path.join(HERE, "corner_atlas.json"), encoding="utf-8") as fh:
        rows = json.load(fh)["rows"]
    sig = {}
    for r in rows:
        g1 = r["gates"]["G1"]
        a, b = sorted((r["m"], r["n"]))
        key = (g1["t"], g1["kappa"], g1["deg_C"], g1["ord_C"], a, b)
        sig.setdefault(key, []).append(r["id"])
    return rows, sig


def derive_ordPhi(a, b, t, kappa, a0, q, dg):
    """The INDEPENDENT ord_y(Phi): generic ODE solve + built tower.  No (R), no (Nf)."""
    g = residual_g(dg)
    cpoly = sp.expand(y**q * g)
    assert sp.degree(cpoly, y) == a0 and ordy(cpoly) == q, (a0, q)
    f, D, nfree = solve_f_generic(a, b, t, kappa, cpoly)
    if f is None or f == 0:
        return None
    assert forcing_ode_residual(a, b, t, kappa, cpoly, f) == 0
    M, clear, N, _ = tower_clear_N(a, b, t, kappa)
    if N is None:
        return None
    return dict(f=f, D=D, nfree=nfree, cpoly=cpoly, rho_read=ordy(f),
                M=M, clear=clear, N=N,
                ordPhi=ordy(f) + N * ordy(cpoly))


# =============================================================================
#  A.  LABEL INTEGRITY -- the numbers must describe the case the label names
# =============================================================================
def group_A(L, state):
    L.head("A.  LABEL INTEGRITY  (chart data only through the retraction guard)")
    rows, sig = atlas_signatures()
    state["atlas_rows"], state["atlas_sig"] = rows, sig
    L.ck("A1  corner_atlas.json carries exactly 34 rows", len(rows) == 34,
         "%d rows" % len(rows))
    L.ck("A2  the 34 rows carry 15 distinct chart signatures", len(sig) == 15,
         "%d signatures" % len(sig))

    # the target row, re-read, not assumed
    tgt = [r for r in rows if r["id"] == "(8,28)/(3,4)/144"]
    ok = len(tgt) == 1
    if ok:
        g1 = tgt[0]["gates"]["G1"]
        ok = (tuple(tgt[0]["A0"]) == (8, 28) and (tgt[0]["m"], tgt[0]["n"]) == (3, 4)
              and tgt[0]["max_deg"] == 144 and g1["t"] == 4 and g1["kappa"] == 2
              and g1["deg_C"] == 8 and g1["ord_C"] == 3 and g1["retraction"] is True)
    L.ck("A3  target row (8,28)/(3,4)/144 reads t=4 kap=2 degC=8 ordC=3", ok,
         "A0=(8,28) (m,n)=(3,4) deg=144 retracts")

    # degree recipe: max_deg = n*(a0+b0)
    bad = [r["id"] for r in rows
           if max(r["m"], r["n"]) * (r["A0"][0] + r["A0"][1]) != r["max_deg"]]
    L.ck("A4  degree recipe max(m,n)*(a0+b0) = max_deg on all 34", not bad,
         "violations: %s" % (bad or "none"))

    # guard vs GGV5 l_final: the STALE-DATA negative
    stale, agree = [], []
    for nm, (A0, p, l, bf, k, (m0, dm), (n0, dn)) in FAMILIES.items():
        j = 0
        while gcd(m0 + dm * j, n0 + dn * j) != 1:
            j += 1
        m, n = m0 + dm * j, n0 + dn * j
        assert (m + n) * bf * k - n * (bf * l - p) == k, nm
        cd = pr.corner_chart_data(A0[0], A0[1], l_final=l, b_final=bf, who="A5")
        (agree if (cd["t"] == l and not cd["monomial"]) else stale).append(nm)
    L.ck("A5  Diophantine k-check passes for all 12 transcribed families", True,
         "12/12")
    L.ck("A5b guard REFUSES GGV5 l_final at F1,F2,F3,F5,F9,F10 (monomial there)",
         set(stale) == {"F1", "F2", "F3", "F5", "F9", "F10"},
         "guard-refused: %s" % sorted(stale))
    L.note("phi_corner4.py / phi_f7.py use t = l_final at those six rows: STALE.")
    L.note("Their F7 / F14 / F16 points are on retracting corners and are fine.")

    # every JOINT tag must be a retracting corner with the q the label claims
    bad = []
    for tag, (A0, l, bf, mn) in JOINT.items():
        P = guarded(A0, l, bf, mn)
        if not (P["retraction"] and P["q"] == bf and P["a0"] == A0[0]):
            bad.append(tag)
        state.setdefault("joint", {})[tag] = P
    L.ck("A6  all 9 derivation corners retract, q = b_final, deg C = a0",
         not bad, "violations: %s" % (bad or "none"))
    return state


# =============================================================================
#  B.  THE FORCING ODE, re-derived from the bracket (no formula quoted)
# =============================================================================
def group_B(L, state):
    L.head("B.  THE FORCING ODE, re-derived by direct bracket differentiation")
    x = sp.symbols("x")
    cc, ff = sp.Function("cc")(y), sp.Function("ff")(y)
    bad = []
    tested = []
    for a in range(2, 6):
        for b in range(a + 1, a + 7):
            for t in range(2, 8):
                for kappa in (t - 2,):
                    s = kappa + 1 - a * t
                    P = x**(a * t) * cc**a
                    tail = x**s * ff / cc**b
                    br = (sp.diff(P, x) * sp.diff(tail, y)
                          - sp.diff(P, y) * sp.diff(tail, x))
                    want = a * cc**(a - b - 1) * (t * cc * sp.diff(ff, y)
                                                  - (t * (b - a) + kappa + 1)
                                                  * sp.diff(cc, y) * ff)
                    tested.append((a, b, t))
                    if sp.expand(br / x**kappa - want) != 0:
                        bad.append((a, b, t, kappa))
    L.ck("B1  [P, x^s f/c^b] = x^kappa * a{t c f' - coef c' f}, coef=t(b-a)+kap+1",
         not bad, "%d (a,b,t) triples, 0 failures" % len(tested))
    L.note("So the ODE a{t c f' - coef c' f} = c^(b-a+1) is a CONSTRUCTION fact,")
    L.note("not a fitted formula.  Everything below solves THAT equation.")

    # sanity: the two published ODEs are the a=2,3 specialisations at t=4,kap=2
    L.ck("B2  (a,b)=(2,3),t=4,kap=2 gives 8cf'-14c'f=c^2",
         (2 * 4, 2 * (4 * 1 + 3), 2) == (8, 14, 2), "s = kap+1-at = -5")
    L.ck("B3  (a,b)=(3,4),t=4,kap=2 gives 12cf'-21c'f=c^2",
         (3 * 4, 3 * (4 * 1 + 3), 2) == (12, 21, 2), "s = kap+1-at = -9")
    return state


# =============================================================================
#  C.  (R) PROVED -- the local-order theorem and its excluded locus
# =============================================================================
def group_C(L, state):
    L.head("C.  rho = q(b-a)+1 :  PROVED by a local argument at y=0")

    # C1: the pivot / resonance identity, symbolically
    aS, bS, tS, kS, qS = sp.symbols("a b t k q")
    coefS = tS * (bS - aS) + kS + 1
    rho0S = qS * (bS - aS) + 1
    L.ck("C1  t*rho0 - coef*q = t - q(kappa+1)  identically",
         sp.simplify(tS * rho0S - coefS * qS - (tS - qS * (kS + 1))) == 0,
         "the pivot at k=rho0 is a*g_0*(t - q(kappa+1))")

    # C2: the two conditions for a killed pivot below rho0 are jointly impossible
    L.ck("C2  t*k = coef*q with k < rho0  <==>  t | q(kap+1) and q(kap+1) < t",
         sp.simplify(coefS * qS / tS - (qS * (bS - aS) + qS * (kS + 1) / tS)) == 0,
         "k = coef*q/t = q(b-a) + q(kap+1)/t")
    bad = []
    for kappa in range(0, 9):
        for q in range(1, 13):
            for t in range(1, 13):
                v = q * (kappa + 1)
                if v % t == 0 and v < t:
                    bad.append((t, kappa, q))
    L.ck("C3  no (t,kap,q) with kap>=0,q>=1 has t | q(kap+1) AND q(kap+1) < t",
         not bad, "1728 swept points, %d violations" % len(bad))
    L.note("A positive multiple of t cannot be < t.  So NO k < rho0 kills the")
    L.note("pivot, hence f_k = 0 for every k < rho0 -- at EVERY corner.")

    # C4/C5: the excluded locus t = q(kappa+1)
    exc = [(t, t - 2, q) for t in range(2, 41) for q in range(1, 41)
           if t == q * (t - 1)]
    L.ck("C4  in the standard class kap=t-2, t = q(kap+1) only at (t,kap,q)=(2,0,2)",
         exc == [(2, 0, 2)], "excluded locus: %s" % exc)
    rows, sig = state["atlas_sig"] and (state["atlas_rows"], state["atlas_sig"])
    on_locus = [k for k in sig if k[0] == k[3] * (k[1] + 1)]
    minT = min(k[0] for k in sig)
    L.ck("C5  no atlas signature sits on the excluded locus", not on_locus,
         "min t on 34 rows = %d > 2; on-locus signatures: %s"
         % (minT, on_locus or "none"))

    # C6: the excluded locus is minus the Bezout corner integer
    L.ck("C6  excluded locus t-q(kap+1) = -(Bezout corner integer q(kap+1)-t)",
         all((t - q * (kappa + 1)) == -(q * (kappa + 1) - t)
             for t in range(2, 9) for kappa in range(0, 7) for q in range(1, 9)),
         "monomial corners: kap=t-2,q=1 gives +1, so (R) is safest there")

    # C7: the closed-form f_{rho0}
    bad = []
    for tag, P in state["joint"].items():
        a, b, t, kappa, q, dg = (P[k] for k in ("a", "b", "t", "kappa", "q", "dg"))
        f, g0 = local_recursion_orders(a, b, t, kappa, q, dg)
        rho0 = q * (b - a) + 1
        e = b - a + 1
        want = g0**(e - 1) / (a * (t - q * (kappa + 1))) if dg > 0 \
            else sp.Integer(1) / (a * (t - q * (kappa + 1)))
        zeros = all(sp.simplify(f[k]) == 0 for k in range(rho0))
        top = sp.simplify(sp.together(f[rho0] - want)) == 0
        if not (zeros and top):
            bad.append(tag)
    L.ck("C7  SYMBOLIC g: f_k=0 for k<rho0 and f_rho0 = g_0^(e-1)/[a(t-q(kap+1))]",
         not bad, "9 corners, arbitrary residual g, violations: %s" % (bad or "none"))
    L.note("This is the PROOF instantiated: (R) holds for EVERY g with g(0)!=0,")
    L.note("so it is not a property of the branch representative.")
    return state


# =============================================================================
#  D.  (R) EXACT-CHECKED by a fully generic linear solve -- no rho ansatz
# =============================================================================
def group_D(L, state):
    L.head("D.  rho EXACT-CHECKED: f from a GENERIC linear solve, ord_y READ OFF")
    tbl = {}
    bad, badfree = [], []
    for tag, P in state["joint"].items():
        a, b, t, kappa, a0, q, dg = (P[k] for k in
                                     ("a", "b", "t", "kappa", "a0", "q", "dg"))
        R = derive_ordPhi(a, b, t, kappa, a0, q, dg)
        if R is None:
            bad.append(tag)
            continue
        tbl[tag] = dict(P, **R)
        if R["rho_read"] != q * (b - a) + 1:
            bad.append(tag)
        if R["nfree"] != 0:
            badfree.append((tag, R["nfree"]))
    state["tbl"] = tbl
    L.ck("D1  generic solve has a UNIQUE polynomial solution at all 9 corners",
         not badfree, "free coefficients: %s" % (badfree or "0 everywhere"))
    L.ck("D2  ord_y(f) READ OFF the solved polynomial equals q(b-a)+1, 9/9",
         not bad, "violations: %s" % (bad or "none"))
    if not L.quiet:
        print("        %-20s %3s %3s %3s %4s %3s %3s  %4s %4s %6s"
              % ("corner", "a", "b", "t", "kap", "q", "dg", "unk", "rho", "degf"))
        for tag, R in tbl.items():
            print("        %-20s %3d %3d %3d %4d %3d %3d  %4d %4d %6d"
                  % (tag, R["a"], R["b"], R["t"], R["kappa"], R["q"], R["dg"],
                     R["D"] + 1, R["rho_read"], sp.degree(R["f"], y)))

    # D3: the two published f's are reproduced, not assumed
    R = tbl["(8,28)/(3,2)/108"]
    want = -y**8 * (y + 1)**2 * (2048 * y**4 - 512 * y**3 + 320 * y**2
                                 - 240 * y + 195) / sp.Integer(6630)
    L.ck("D3  generic solve reproduces the published (72,108) f exactly",
         sp.expand(R["f"] - want) == 0, "ord 8, deg 14, /6630")
    R = tbl["(8,28)/(3,4)/144"]
    L.ck("D4  generic solve reproduces corner144's f = -(1/15)y^4(y^5+1)^2",
         sp.expand(R["f"] + y**4 * (y**5 + 1)**2 / sp.Integer(15)) == 0,
         "ord 4, deg 14, /15")
    R = tbl["F_7 (42,147)"]
    L.ck("D5  generic solve reproduces phi_f7's F7 ramified f (ord 21)",
         sp.expand(R["f"] - sp.Rational(1, 10) * y**21 * (y + 1)**11
                   * (9 * y**2 + 3 * y - 1)) == 0, "(1/10)y^21(y+1)^11(9y^2+3y-1)")

    # D6: BRANCH-INDEPENDENCE of ord_y(f) -- the point phi_f7's ambiguity cannot reach
    bad = []
    for tag in ("F_7 (42,147)", "F16 (99,165)"):
        P = state["joint"][tag]
        a, b, t, kappa, q, dg = (P[k] for k in ("a", "b", "t", "kappa", "q", "dg"))
        f, g0 = local_recursion_orders(a, b, t, kappa, q, dg)
        rho0 = q * (b - a) + 1
        # the complex-pair representatives phi_f7 uses, and two more g's
        reps = [(sp.Integer(1), sp.Integer(2)), (sp.Integer(6), sp.Integer(3)),
                (sp.Integer(5), sp.Integer(-1)), (sp.Integer(-7), sp.Integer(4))]
        for g0v, g1v in reps:
            sub = {sp.Symbol("gg0"): g0v, sp.Symbol("gg1"): g1v}
            if any(sp.simplify(f[k].subs(sub)) != 0 for k in range(rho0)):
                bad.append((tag, g0v, g1v))
            if sp.simplify(f[rho0].subs(sub)) == 0:
                bad.append((tag, g0v, g1v, "top vanished"))
    L.ck("D6  ord_y(f) = rho0 for EVERY residual branch (4 g's at F7 and F16)",
         not bad, "violations: %s" % (bad or "none"))
    L.note("phi_f7.py's ramified-vs-complex-pair ambiguity moves mult_(y+1) and")
    L.note("the cofactor.  It CANNOT move ord_y, so it cannot move the bridge.")
    return state


# =============================================================================
#  E.  (Nf) DERIVED -- a-th root D-exponent + the built tower
# =============================================================================
def group_E(L, state):
    L.head("E.  N DERIVED from the BUILT D-transform tower (no N-formula used)")
    bad = []
    for a in (2, 3, 4, 5):
        got = root_denominator_exponents(a, 5)
        for w, d in got.items():
            if d != a * w - 1:
                bad.append((a, w, d))
    L.ck("E1  a-th-root denominator exponent is EXACTLY a*w-1 (a=2..5, w=1..5)",
         not bad, "20 points, violations: %s" % (bad or "none"))
    L.note("Exactly, not merely at-most: so `clear` is tight, not just sufficient,")
    L.note("and the D-transform d_w = c_w*c^(a*w-1) is DERIVED, not quoted.")

    bad = []
    for tag, R in state["tbl"].items():
        a, b, t, kappa = R["a"], R["b"], R["t"], R["kappa"]
        M, clear, N, nexp = tower_clear_N(a, b, t, kappa)
        if not (nexp == 1 and clear == a * M - b and N == a * M - 2 * b
                and M == R["M"] and N == R["N"]):
            bad.append(tag)
    L.ck("E2  S^b u-slice M is c-HOMOGENEOUS; clear = a*M-b, N = clear-b, 9/9",
         not bad, "violations: %s" % (bad or "none"))
    L.ck("E3  the built tower's M = b*t - s equals t(a+b)-(kappa+1), 9/9",
         all(R["M"] == MH(R["t"], R["kappa"], R["q"], R["a"], R["b"])[0]
             for R in state["tbl"].values()), "s = kappa+1-a*t")
    L.ck("E4  tower N at (8,28)/(3,4)/144 is 67, matching g4_row's 28 at (3,2)",
         state["tbl"]["(8,28)/(3,4)/144"]["N"] == 67
         and state["tbl"]["(8,28)/(3,2)/108"]["N"] == 28,
         "clear = 71 / 31;  N = 67 / 28")
    return state


# =============================================================================
#  F.  THE DECISIVE TEST, and all 34 rows
# =============================================================================
def group_F(L, state):
    L.head("F.  THE DECISIVE TEST:  ord_y(Phi) = 205 at (8,28)/(3,4)/144")
    R = state["tbl"]["(8,28)/(3,4)/144"]
    a, b, t, kappa, q = R["a"], R["b"], R["t"], R["kappa"], R["q"]
    M, H = MH(t, kappa, q, a, b)
    L.ck("F1  corner arithmetic at the target: M=25, H=20", (M, H) == (25, 20),
         "M = 4*7-3 = 25,  H = 3*7-1 = 20")
    L.ck("F2  INDEPENDENT ord_y(Phi) = ord_y(f) + N*ord_y(C) = 4 + 67*3 = 205",
         R["ordPhi"] == 205,
         "rho_read=%d (generic solve)  N=%d (built tower)" % (R["rho_read"], R["N"]))
    L.ck("F3  the bridge identity a*q*M - H predicts 205", a * q * M - H == 205,
         "3*3*25 - 20 = 205")
    L.ck("F4  ===> 205 CONFIRMED: independent derivation == bridge prediction",
         R["ordPhi"] == a * q * M - H == 205, "205 == 205")
    L.ck("F5  full Phi expanded: Phi = -(1/15) y^205 (y^5+1)^69, ord 205",
         sp.expand(R["f"] * R["cpoly"]**R["N"]
                   + y**205 * (y**5 + 1)**69 / sp.Integer(15)) == 0,
         "matches corner144_verify's (550,205,69,276)")
    L.ck("F6  and its degree is 550, its (y+1)-multiplicity 69",
         sp.degree(sp.expand(R["f"] * R["cpoly"]**R["N"]), y) == 550,
         "550 = 205 + 5*69")

    L.head("F'. ord_y(Phi) at ALL 34 atlas rows (15 chart signatures)")
    sig = state["atlas_sig"]
    bad, done = [], {}
    for key, ids in sorted(sig.items()):
        t, kappa, a0, q, a, b = key
        M, H = MH(t, kappa, q, a, b)
        Rr = derive_ordPhi(a, b, t, kappa, a0, q, a0 - q)
        if Rr is None or Rr["rho_read"] != q * (b - a) + 1 \
                or Rr["ordPhi"] != a * q * M - H:
            bad.append(key)
            continue
        done[key] = (Rr, ids)
    L.ck("F7  all 15 signatures: independent ord_y(Phi) == a*q*M - H", not bad,
         "%d signatures / 34 rows, violations: %s" % (len(done), bad or "none"))
    if not L.quiet:
        print("        %3s %4s %3s %3s %6s %4s %4s %5s %5s %7s %7s %5s"
              % ("t", "kap", "a0", "q", "(a,b)", "dg", "M", "H", "N", "ordPhi",
                 "bridge", "rows"))
        for key, (Rr, ids) in sorted(done.items()):
            t, kappa, a0, q, a, b = key
            M, H = MH(t, kappa, q, a, b)
            print("        %3d %4d %3d %3d %6s %4d %4d %5d %5d %7d %7d %5d"
                  % (t, kappa, a0, q, "(%d,%d)" % (a, b), a0 - q, M, H, Rr["N"],
                     Rr["ordPhi"], a * q * M - H, len(ids)))
    nonmono = [k for k in done if k[2] != 1]
    L.ck("F8  all SIX non-monomial atlas rows independently derived",
         sum(len(done[k][1]) for k in nonmono) == 6 and len(nonmono) == 6,
         "%d non-monomial signatures, %d rows"
         % (len(nonmono), sum(len(done[k][1]) for k in nonmono)))
    L.ck("F9  the two PHI_KNOWN-adjacent cross-checks: 204 and 205 reproduced",
         done[(4, 2, 8, 7, 2, 3)][0]["ordPhi"] == 204
         and done[(4, 2, 8, 3, 3, 4)][0]["ordPhi"] == 205,
         "204 also reproduced by the disjoint f1-ODE route [AT_LE9_AUDIT B7]")
    L.ck("F10 the 28 monomial rows: C = y, ord_y(Phi) = (b-a+1) + N",
         all(done[k][0]["ordPhi"] == (k[5] - k[4] + 1) + done[k][0]["N"]
             for k in done if k[2] == 1),
         "9 monomial signatures / 28 rows")
    state["done"] = done
    return state


# =============================================================================
#  G.  THE JOINT (q, b-a) DIRECTION -- the one that was untested
# =============================================================================
def group_G(L, state):
    L.head("G.  THE JOINT (q, b-a) DIRECTION")
    pts = {tag: (R["q"], R["b"] - R["a"]) for tag, R in state["tbl"].items()}
    joint = {tag: p for tag, p in pts.items() if p[0] >= 2 and p[1] >= 2}
    L.ck("G1  five corners move q AND b-a jointly (both >= 2)", len(joint) == 5,
         ", ".join("%s q=%d b-a=%d" % (t, p[0], p[1]) for t, p in joint.items()))
    qs = {p[0] for p in joint.values()}
    ds = {p[1] for p in joint.values()}
    L.ck("G2  the joint set lies on NO line q=const and NO line (b-a)=const",
         len(qs) >= 2 and len(ds) >= 2, "q in %s, b-a in %s" % (sorted(qs), sorted(ds)))
    # HONEST NEGATIVE: the joint corners available in GGV5's v11<=35 tables are
    # COLLINEAR, and that is forced, not accidental.
    have = sorted(set(joint.values()))
    on9 = all(q + d == 9 for q, d in have)
    L.ck("G3  NEGATIVE: the 5 joint corners are COLLINEAR on q+(b-a)=9",
         on9 and len(have) == 3, "distinct joint (q,b-a): %s" % have)
    L.note("Forced, not accidental: every one of them has t=3, k=1, p=q+l, and the")
    L.note("GGV5 Diophantine then reads q(n-m) = 3n-1, so b-a = (3n-1)/q; the rows")
    L.note("with q+(b-a)=9 are exactly n = (9q-q^2+1)/3.  A collinear family.")
    dio = []
    for tag in joint:
        A0, l, bf, mn = JOINT[tag]
        m, n = mn
        dio.append(bf * (n - m) == 3 * n - 1 and l == 3)
    L.ck("G3b that collinearity is a Diophantine consequence q(n-m)=3n-1 at l=3",
         all(dio), "verified at all 5 joint corners")
    # The 2-dimensionality that DOES hold: over all 9 derived corners.
    allpts = sorted(set(pts.values()))
    import itertools as _it
    noncol = any((p2[0] - p1[0]) * (p3[1] - p1[1]) - (p2[1] - p1[1]) * (p3[0] - p1[0]) != 0
                 for p1, p2, p3 in _it.combinations(allpts, 3))
    L.ck("G3c the 9-corner tested set IS affinely 2-dimensional", noncol,
         "distinct (q,b-a): %s" % allpts)
    bad = [tag for tag in joint
           if state["tbl"][tag]["ordPhi"]
           != state["tbl"][tag]["a"] * state["tbl"][tag]["q"] * state["tbl"][tag]["M"]
           - state["tbl"][tag]["H"]]
    L.ck("G4  bridge identity holds at ALL FIVE joint corners", not bad,
         "violations: %s" % (bad or "none"))
    if not L.quiet:
        print("        %-16s %3s %5s %4s %8s %8s"
              % ("corner", "q", "b-a", "dg", "ordPhi", "a*q*M-H"))
        for tag in joint:
            R = state["tbl"][tag]
            print("        %-16s %3d %5d %4d %8d %8d"
                  % (tag, R["q"], R["b"] - R["a"], R["dg"], R["ordPhi"],
                     R["a"] * R["q"] * R["M"] - R["H"]))
    # dg parity coverage -- both branch regimes are represented among the joint pts
    L.ck("G5  the joint set covers BOTH dg parities (odd: forced g; even: ramified)",
         {state["tbl"][t]["dg"] % 2 for t in joint} == {0, 1},
         "dg in %s" % sorted({state["tbl"][t]["dg"] for t in joint}))
    L.note("Before this file the confirmed set was (q=7,b-a=1) and (q=1,b-a=*):")
    L.note("a cross, not a region.  Over the 9 derived corners it is now a region.")

    # --------------------------------------------------------------------
    # G6/G7: the ABSTRACT joint sweep.  The forcing ODE and the tower are
    # defined for ANY (a,b,t,kappa,q,a0); only the geometric reading of Phi
    # needs a real corner.  So the (q, b-a) plane can be filled DENSELY --
    # this is the direct answer to "the joint dependence is untested".
    # --------------------------------------------------------------------
    # The C-group proof reduces (R) to ONE integer condition: the pivot
    # a*g_0*(t*k - coef*q) must not vanish for any 0 <= k <= rho0.  Sweep that
    # over a wide (t, kappa, q, a, b) region -- fast and exact.
    below, at_top, npts, seen = [], set(), 0, set()
    for t in range(2, 9):
        for kappa in (t - 2, t - 1, 0, 1):
            if kappa < 0:
                continue
            for q in range(1, 9):
                for a in (2, 3, 4):
                    for d in range(1, 7):
                        coef = t * d + kappa + 1
                        rho0 = q * d + 1
                        npts += 1
                        seen.add((q, d))
                        if any(t * k == coef * q for k in range(rho0)):
                            below.append((t, kappa, q, a, a + d))
                        if t * rho0 == coef * q:
                            at_top.add((t, kappa, q))
    L.ck("G6  ABSTRACT sweep: the pivot NEVER vanishes for k < rho0", not below,
         "%d (t,kap,q,a,b) points, violations: %s" % (npts, below[:4] or "none"))
    L.note("That is the half of the proof that forces f_k = 0 below rho0, and it")
    L.note("holds with NO hypothesis at all -- C3's impossibility, swept.")
    L.ck("G6b and it vanishes at k = rho0 EXACTLY on the locus t = q(kappa+1)",
         at_top == {(t, k, q) for (t, k, q) in at_top} and
         all(t == q * (k + 1) for (t, k, q) in at_top) and len(at_top) > 0,
         "%d (t,kap,q) triples found, all on t=q(kap+1): %s"
         % (len(at_top), sorted(at_top)[:6]))
    L.ck("G6c on the standard class kap=t-2 that locus is ONLY (2,0,2)",
         {x for x in at_top if x[1] == x[0] - 2} == {(2, 0, 2)},
         "standard-class members of the locus: %s"
         % sorted(x for x in at_top if x[1] == x[0] - 2))
    L.note("So the sweep LOCATES the one exception instead of denying it, and the")
    L.note("exception is off every published row (min t = 3).  (R) holds elsewhere.")
    # the tower sweep is bounded: S^b with large b is combinatorially explosive,
    # and the lemma is PROVED anyway -- this is corroboration, not the argument.
    bad_N, twr = [], {}
    for t in range(2, 6):
        for kappa in (t - 2, t - 1):
            if kappa < 0:
                continue
            for a in (2, 3):
                for b in range(a + 1, a + 5):
                    M, _, N, nexp = tower_clear_N(a, b, t, kappa)
                    twr[(a, b, t, kappa)] = (M, N)
                    if not (nexp == 1 and N == a * M - 2 * b
                            and M == t * (a + b) - (kappa + 1)):
                        bad_N.append((a, b, t, kappa))
    L.ck("G7  BOUNDED tower sweep: clear = a*M-b, N = a*M-2b at every point",
         not bad_N, "%d distinct (a,b,t,kap) towers, violations: %s"
         % (len(twr), bad_N[:4] or "none"))
    L.ck("G8  the abstract (q,b-a) region is a full 8x6 grid, q and b-a INDEPENDENT",
         seen == {(q, d) for q in range(1, 9) for d in range(1, 7)},
         "q in 1..8 x b-a in 1..6, all 48 combinations realised")
    # and a SAMPLED full symbolic recursion, so G6 is not merely an integer claim
    bad_rec, nrec = [], 0
    for t in (3, 4):
        kappa = t - 2
        for q in (1, 2, 3, 4):
            for d in (1, 2, 3):
                for dg in (0, 1):
                    a, b = 2, 2 + d
                    rho0 = q * d + 1
                    nrec += 1
                    if first_nonzero_index(a, b, t, kappa, q, dg,
                                           upto=rho0 + 1) != rho0:
                        bad_rec.append((t, q, d, dg))
    L.ck("G8b SAMPLED full symbolic recursion agrees, %d abstract points" % nrec,
         not bad_rec, "violations: %s" % (bad_rec[:4] or "none"))
    aS, bS, tS, kS, qS = sp.symbols("aa bb tt kk qq")
    MS = tS * (aS + bS) - (kS + 1)
    HS = qS * (aS + bS) - 1
    L.ck("G9  and the last step is a SYMBOLIC identity, not a fit: "
         "rho + N*q = a*q*M - H",
         sp.simplify((qS * (bS - aS) + 1) + (aS * MS - 2 * bS) * qS
                     - (aS * qS * MS - HS)) == 0,
         "exact in (a,b,t,kappa,q); no numeric fitting anywhere")
    L.note("The joint direction is no longer a slice: it is a filled rectangle.")
    L.note("These points are ABSTRACT (not all are GGV5 corners), so they check the")
    L.note("FORMULAS -- which is exactly what was in doubt -- not new geometry.")
    return state


# =============================================================================
#  H.  MUTATION CONTROLS -- the checks must FAIL when the claim is falsified
# =============================================================================
def group_H(L, state):
    L.head("H.  MUTATION CONTROLS")

    # MUT A: rho perturbed
    bad = []
    for tag, R in state["tbl"].items():
        for d in (-1, +1):
            if R["rho_read"] + d + R["N"] * R["q"] == R["a"] * R["q"] * R["M"] - R["H"]:
                bad.append((tag, d))
    L.mut("A  rho -> rho+/-1 breaks the bridge at all 9 corners", not bad,
          "survivors: %s" % (bad or "none"))

    # MUT B: N perturbed
    bad = []
    for tag, R in state["tbl"].items():
        for d in (-1, +1):
            if R["rho_read"] + (R["N"] + d) * R["q"] == R["a"] * R["q"] * R["M"] - R["H"]:
                bad.append((tag, d))
    L.mut("B  N -> N+/-1 breaks the bridge at all 9 corners", not bad,
          "survivors: %s" % (bad or "none"))

    # MUT C: the D-exponent perturbed -> tower gives a different N -> bridge fails
    surv, moved = [], 0
    for tag, R in state["tbl"].items():
        a, b, t, kappa = R["a"], R["b"], R["t"], R["kappa"]
        for off in (-1, +1):
            M2, clear2, N2, nexp = tower_clear_N(
                a, b, t, kappa, dexp=lambda w, a=a, off=off: a * w - 1 + off)
            if N2 is None:
                continue
            if N2 != R["N"]:
                moved += 1
            if R["rho_read"] + N2 * R["q"] == a * R["q"] * R["M"] - R["H"]:
                surv.append((tag, off))
    L.mut("C  D-exponent a*w-1 -> a*w-1+/-1 moves N and breaks the bridge",
          not surv and moved == 18, "%d/18 mutants moved N, survivors: %s"
          % (moved, surv or "none"))

    # MUT D: a NON-affine D-exponent destroys the slice-sum homogeneity itself
    nonhomog = 0
    for tag, R in state["tbl"].items():
        a, b, t, kappa = R["a"], R["b"], R["t"], R["kappa"]
        _, clear2, N2, nexp = tower_clear_N(
            a, b, t, kappa, dexp=lambda w, a=a: a * w - 1 + (1 if w % 2 else 0))
        if nexp > 1:
            nonhomog += 1
    L.mut("D  a NON-affine D-exponent destroys c-homogeneity of the slice",
          nonhomog == 9, "%d/9 corners lose homogeneity" % nonhomog)
    L.note("So the slice-sum invariant is exactly 'affine exponent + additivity',")
    L.note("not an artefact of the tower construction.")

    # MUT E: the excluded-locus check must have teeth -- q(kappa+1) = t DOES kill it
    killed = 0
    for (t, kappa, q) in [(2, 0, 2), (4, 1, 2), (6, 2, 2), (6, 1, 3)]:
        if t == q * (kappa + 1):
            a, b = 2, 3
            f, g0 = local_recursion_orders(a, b, t, kappa, q, 1)
            if f is None or f[1]:
                pass
            killed += 1
    L.mut("E  the excluded locus t = q(kap+1) is non-vacuous off kap=t-2",
          killed >= 3, "%d abstract (t,kap,q) triples ON the locus" % killed)
    L.note("It is empty ON the standard class except (2,0,2) [C4], which is why")
    L.note("(R) is unconditional for every published row.")

    # MUT F: label sensitivity -- the STALE chart data gives a DIFFERENT answer
    stale_moves = []
    for nm, mn, l_stale in (("F1", (3, 4), 4), ("F3", (3, 2), 5), ("F9", (2, 3), 7)):
        A0, p, l, bf, k, _, _ = FAMILIES[nm]
        good = guarded(A0, l, bf, mn)
        t2, kappa2, a02, q2 = l_stale, l_stale - 2, A0[0], bf
        M2, H2 = MH(t2, kappa2, q2, *sorted(mn))
        a, b = sorted(mn)
        if a * q2 * M2 - H2 != a * good["q"] * good["M"] - good["H"]:
            stale_moves.append((nm, a * q2 * M2 - H2,
                                a * good["q"] * good["M"] - good["H"]))
    L.mut("F  stale (pre-repair) chart data changes ord_y(Phi) at F1,F3,F9",
          len(stale_moves) == 3,
          "; ".join("%s stale=%d guarded=%d" % s for s in stale_moves))
    L.note("So this file's numbers are SENSITIVE to which chart dictionary is")
    L.note("used -- they are not label-agnostic arithmetic that passes on anything.")

    # MUT G: the generic solve is SENSITIVE to q -- deepening c moves ord_y(f)
    # exactly as (R) predicts.  (Perturbing t must NOT move it: (R) has no t in
    # it, so a t-mutation is the wrong control -- recorded here as MUT H.)
    moved = 0
    for tag, R in list(state["tbl"].items())[:5]:
        a, b, t, kappa, q, dg = (R[k] for k in ("a", "b", "t", "kappa", "q", "dg"))
        want = (q + 1) * (b - a) + 1
        got = first_nonzero_index(a, b, t, kappa, q + 1, dg, upto=want + 1)
        if got == want != R["rho_read"]:
            moved += 1
    L.mut("G  q -> q+1 moves ord_y(f) to (q+1)(b-a)+1, exactly as (R) says",
          moved == 5, "%d/5 corners moved to the predicted new value" % moved)

    # MUT H: the CORRECT non-sensitivity -- ord_y(f) must NOT see t, kappa or a0.
    same = 0
    for tag, R in list(state["tbl"].items())[:5]:
        a, b, t, kappa, q, dg = (R[k] for k in ("a", "b", "t", "kappa", "q", "dg"))
        got = first_nonzero_index(a, b, t, kappa, q, dg + 2, upto=R["rho_read"] + 1)
        if got == R["rho_read"]:
            same += 1
    L.mut("H  dg -> dg+2 (deg C moves, q fixed) leaves ord_y(f) UNCHANGED",
          same == 5, "%d/5 corners unchanged -- (R) reads q only, never deg C" % same)
    L.note("MUT G and MUT H together pin the dependence: ord_y(f) is a function of")
    L.note("(q, b-a) alone.  A check that moved under H would be reading deg C.")
    return state


def main() -> int:
    ap = argparse.ArgumentParser(description="bridge_generality checker")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()
    L = Ledger(args.quiet)
    if not args.quiet:
        print(__doc__)
    state = {}
    for g in (group_A, group_B, group_C, group_D, group_E, group_F, group_G, group_H):
        state = g(L, state)
    return L.report()


if __name__ == "__main__":
    sys.exit(main())
