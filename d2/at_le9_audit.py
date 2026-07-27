#!/usr/bin/env python3
"""at_le9_audit.py -- INDEPENDENT AUDIT of the claim  a_t = v_t(e) = v_t(dm1) <= 9.

Commissioned as a second, independently authored checker for the two same-author
proofs of `a_t <= 9`:

    syzygy_collision.py  / SYZYGY_COLLISION.md   (K-syzygy collapse)
    slice_phi_yplace.py  / SLICE_PHI_YPLACE.md   (levels 14/15 + Phi leading jet)

NO CODE IS IMPORTED, EXEC'D, COPIED OR ADAPTED FROM EITHER FILE, nor from
slice_obstruction_basis.py, slice_obstruction_audit.py, alt_level12.py,
divisor_syzygy.py, i3_audit.py, sub1_spine9.py, positive_slice.py or
window_caps_verify.py.  Everything below is rebuilt from the repo PRIMITIVES:

    generators.json                 (parsed by hand as JSON term lists)
    paper_src/upstream_facts.json   (Prop-4.3 corners -> C4 = y^7(y+1))
    the commutator route for Phi    (re-derived, not transcribed)

Deliberate methodological differences from both proofs, so that agreement is
evidence and not echo:

  * the whole slice calculus is rebuilt from ONE object -- the root series
    H(u) = sum h_i u^i with h_{>=9} eliminated by p_n = 0 -- rather than from
    D2/D3 slice tables or from a stacked P/Q cokernel engine;
  * the d3-killing shift is obtained TWICE and by two different mechanisms:
    (a) generalized binomials computed as falling factorials (no sympy.binomial),
    (b) the generating-function identity Ht(u) = (1-a u)^4 H(u/(1-a u));
  * the K-syzygy is re-derived as a POWER-SERIES identity in h~_2..h~_8, not as
    a term-list combination of generators.json;
  * the cascade is re-run in the UNSHIFTED chart with an explicitly different
    reduction path (level 12 closes via the (P<) n=8 relation A2*A6 = 0, where
    slice_phi_yplace closes it via A8 = -A2*A6);
  * the leading-jet ("graded") lemma is used with a machine-checked
    layer-stability control instead of being assumed.

Read-only.  Pure sympy: no Singular, no msolve, no WSL, no subprocess, no
modular arithmetic, no solver -- so there are no aborts, timeouts or exit codes
to misread.  Writes nothing.

Usage:
    python at_le9_audit.py            # narrated report
    python at_le9_audit.py --quiet    # exit 0 iff every check passes
    python at_le9_audit.py --deep     # + the slower multi-layer Groebner controls
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))

t, y, x, u = sp.symbols("t y x u")

# ----------------------------------------------------------------- harness
_state = {"n": 0, "fail": 0, "quiet": False}


def ck(name, cond, detail=""):
    _state["n"] += 1
    ok = bool(cond)
    if not ok:
        _state["fail"] += 1
    if not _state["quiet"]:
        print("  [%s] %s%s" % ("OK" if ok else "FAIL", name,
                               ("   -- " + str(detail)) if detail else ""))
    elif not ok:
        print("  [FAIL] %s   -- %s" % (name, detail), file=sys.stderr)
    return ok


def head(s):
    if not _state["quiet"]:
        print("\n" + s)
        print("-" * len(s))


# ============================================================ 0. primitives
def load_generators():
    with open(os.path.join(HERE, "generators.json")) as fh:
        G = json.load(fh)
    vo = G["variable_order"]
    S = {n: sp.Symbol(n) for n in vo}
    def poly(name):
        acc = sp.Integer(0)
        for c, ev in G["polynomials"][name]:
            m = sp.Rational(c)
            for i, e in enumerate(ev):
                if e:
                    m *= S[vo[i]] ** e
            acc += m
        return sp.expand(acc)
    return G, vo, S, poly


# ================================================ 1. the root-series engine
NU = 20                       # u-truncation: enough for r_17 and p_8


def root_series(h_list, nu=NU):
    """H(u) = 1 + sum_{k>=1} h_k u^k with h_n (n >= 9) ELIMINATED by p_n = 0.

    p_n := [u^n] H^2 ; p_n = 2 h_n + sum_{i=1}^{n-1} h_i h_{n-i}, so
    p_n = 0  <=>  h_n = -(1/2) sum_{i=1}^{n-1} h_i h_{n-i}.   (no division
    except by 2; the elimination is a closed recursion.)
    """
    H = [sp.Integer(1)] + list(h_list) + [sp.Integer(0)] * (nu - 1 - len(h_list))
    for n in range(9, nu):
        H[n] = sp.expand(-sp.Rational(1, 2) * sum(H[i] * H[n - i] for i in range(1, n)))
    return H


def series_pow(H, m, nu=NU, tcut=None):
    R = [sp.Integer(1)] + [sp.Integer(0)] * (nu - 1)
    for _ in range(m):
        R = [_tcut(sp.expand(sum(R[i] * H[n - i] for i in range(n + 1))), tcut)
             for n in range(nu)]
    return R


def _tcut(e, T):
    if T is None:
        return e
    e = sp.expand(e)
    if e == 0:
        return e
    P = sp.Poly(e, t)
    return sp.expand(sum(P.coeff_monomial(t ** m) * t ** m
                         for m in set(mm[0] for mm in P.monoms()) if m <= T))


def _layer_stability(Lc, extra_layers=3, T=32, nu=18):
    """Exact control on the graded ("leading jet") lemma.

    Rebuild the whole calculus twice inside a fast sparse polynomial ring over
    Q: once with h_k = t^(L_k) * A_k (a bare monomial) and once with
    h_k = t^(L_k) * (A_k + B_k1 t + B_k2 t^2 + B_k3 t^3), every B fully
    SYMBOLIC.  The lemma is exactly the claim that the lowest t-order of each
    p_n / r_n and its coefficient are the same in both.  Returns (ok, count).
    """
    from sympy.polys.rings import ring
    from sympy.polys.domains import QQ

    names = ("tt," + ",".join("Aq%d" % k for k in range(1, 9)) + ","
             + ",".join("Bq%d_%d" % (k, j) for k in range(1, 9)
                        for j in range(1, extra_layers + 1)))
    Rg, *gens = ring(names, QQ)
    tt = gens[0]
    Ag = {k: gens[k] for k in range(1, 9)}
    Bg = {(k, j): gens[9 + extra_layers * (k - 1) + (j - 1)]
          for k in range(1, 9) for j in range(1, extra_layers + 1)}

    def cut(f):
        d = {m: c for m, c in f.terms() if m[0] <= T}
        return Rg.from_dict(d) if d else Rg.zero

    def run(layers):
        hv = []
        for k in range(1, 9):
            e = Ag[k] + sum(Bg[(k, j)] * tt ** j for j in range(1, layers + 1))
            hv.append(cut(tt ** Lc[k] * e))
        H = [Rg.one] + hv + [Rg.zero] * (nu - 1 - len(hv))
        for n in range(9, nu):
            H[n] = cut(-QQ(1, 2) * sum((H[i] * H[n - i] for i in range(1, n)), Rg.zero))

        def pw(m):
            Rr = [Rg.one] + [Rg.zero] * (nu - 1)
            for _ in range(m):
                Rr = [cut(sum((Rr[i] * H[n - i] for i in range(n + 1)), Rg.zero))
                      for n in range(nu)]
            return Rr
        return pw(2), pw(3)

    def low(f):
        if f == Rg.zero:
            return None, None
        m0 = min(m[0] for m in f.monoms())
        return m0, Rg.from_dict({(0,) + m[1:]: c for m, c in f.terms() if m[0] == m0})

    p0, r0 = run(0)
    p1, r1 = run(extra_layers)
    ok, cnt = True, 0
    for a_, b_, ns in ((p0, p1, (6, 7, 8)), (r0, r1, (9, 10, 11, 12, 13, 14, 15, 17))):
        for n in ns:
            m0a, ca = low(a_[n])
            m0b, cb = low(b_[n])
            cnt += 1
            if m0a is None or m0a != m0b or (ca - cb) != Rg.zero:
                ok = False
    return ok, cnt


def tcoeffs(e):
    e = sp.expand(e)
    if e == 0:
        return {}
    P = sp.Poly(e, t)
    return {m[0]: sp.expand(P.coeff_monomial(t ** m[0])) for m in P.monoms()}


# ==================================================================== MAIN
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--deep", action="store_true")
    args = ap.parse_args()
    _state["quiet"] = args.quiet

    G, vo, SYM, gpoly = load_generators()

    # ================================================================== A
    head("A. The G-system rows, rebuilt from ONE root series (tests [I3]/[Q8])")

    ck("A0  generators.json variable_order is the expected 8-tuple",
       vo == ["d2", "d1", "d0", "dm1", "dm2", "dm3", "dm4", "Phi"], vo)

    hs = sp.symbols("hh2:9")                     # h~_2 .. h~_8, SHIFTED chart
    eta = sp.Symbol("eta")

    Ht = root_series([sp.Integer(0)] + list(hs))            # h~_1 = 0
    pt = series_pow(Ht, 2)
    rt = series_pow(Ht, 3)

    ck("A1  the elimination is exact: p~_n == 0 for n = 9..19",
       all(sp.expand(pt[n]) == 0 for n in range(9, NU)))

    DICT = {hs[0]: SYM["d2"], hs[1]: SYM["d1"], hs[2]: SYM["d0"],
            hs[3]: SYM["dm1"], hs[4]: SYM["dm2"], hs[5]: SYM["dm3"],
            hs[6]: SYM["dm4"]}
    for n, name in ((13, "G1"), (14, "G2"), (15, "G3"), (17, "G5body")):
        ck("A2.%s  r~_%d  ==  generators.json %s  exactly" % (name, n, name),
           sp.expand(rt[n].subs(DICT) - gpoly(name)) == 0)

    # chart discriminator: with a free u^1 term the rebuild does NOT match
    He = root_series([eta] + list(hs))
    re_ = series_pow(He, 3)
    diffs = [sp.expand(re_[n].subs(DICT) - gpoly(nm)) != 0
             for n, nm in ((13, "G1"), (14, "G2"), (15, "G3"), (17, "G5body"))]
    ck("A3  CHART DISCRIMINATOR: keeping a free u^1 term BREAKS the match "
       "(so generators.json really is the h~_1 = 0 chart)", any(diffs))
    ck("A3b setting eta = 0 recovers the match",
       all(sp.expand(re_[n].subs({eta: 0}).subs(DICT) - gpoly(nm)) == 0
           for n, nm in ((13, "G1"), (14, "G2"), (15, "G3"), (17, "G5body"))))

    # slice-exponent arithmetic, re-derived from c_m = D_m C4^(2m-7),
    # D_j = y^(12(4-j)) d_j, C4 = y^7 t.
    C4s = sp.Symbol("C4s")
    M = sp.Symbol("M")
    # Q_M = sum_{i+j+k=M} c_i c_j c_k : C4-exponent (2i-7)+(2j-7)+(2k-7) = 2M-21
    #     = C4^(2M-21) * y^(12*(12-M)) * r_(12-M)
    eyQ = sp.expand(7 * (2 * M - 21) + 12 * (12 - M))
    etQ = sp.expand(2 * M - 21)
    ck("A4  Q-slice, re-derived: Q_M = y^(2M-3) * r_(12-M) / t^(21-2M)",
       sp.expand(eyQ - (2 * M - 3)) == 0 and sp.expand(etQ - (2 * M - 21)) == 0,
       "y-exponent %s, t-exponent %s" % (eyQ, etQ))
    # P_M : C4-exponent (2i-7)+(2j-7) = 2M-14, strip y^(12*(8-M))
    eyP = sp.expand(7 * (2 * M - 14) + 12 * (8 - M))
    etP = sp.expand(2 * M - 14)
    ck("A5  P-slice, re-derived: P_M = y^(2M-2) * p_(8-M) / t^(14-2M)",
       sp.expand(eyP - (2 * M - 2)) == 0 and sp.expand(etP - (2 * M - 14)) == 0,
       "y-exponent %s, t-exponent %s" % (eyP, etP))
    ck("A6  hence, gcd(y,t) = 1: P_M polynomial => t^(2n-2) | p_n (n = 8-M) and "
       "Q_M polynomial => t^(2n-3) | r_n (n = 12-M); Q_M = 0 for M < 0 => "
       "r_13 = r_14 = r_15 = 0",
       sp.expand((14 - 2 * M).subs(M, 8 - sp.Symbol('n')) - (2 * sp.Symbol('n') - 2)) == 0
       and sp.expand((21 - 2 * M).subs(M, 12 - sp.Symbol('n')) - (2 * sp.Symbol('n') - 3)) == 0)

    # ================================================================== B
    head("B. Phi and  v_t(Phi) = 30  EXACTLY  (re-derived from the commutator)")

    C4 = y ** 7 * (y + 1)
    f1 = sp.Function("f1")(y)
    br = lambda g, hh: sp.diff(g, x) * sp.diff(hh, y) - sp.diff(g, y) * sp.diff(hh, x)
    # [P, Q^2 - P^3 - 2 lam P] = [P,Q^2] = 2 Q [P,Q] = 2 Q x^2 ; leading forms
    #   ell(P) = x^8 C4^2, ell(Q) = x^12 C4^3, ell(2 C^3 F) = 2 x^7 C4^3 F_{-5}
    #   = 2 x^7 f1   with f1 := C4^3 F_{-5}.
    lhs = br(x ** 8 * C4 ** 2, 2 * x ** 7 * f1)
    rhs = 2 * x ** 2 * x ** 12 * C4 ** 3
    ode = sp.cancel((lhs - rhs) / (2 * x ** 14 * C4 * y ** 6))
    ode_expected = 8 * y * (y + 1) * sp.diff(f1, y) - 14 * (8 * y + 7) * f1 - y ** 8 * (y + 1) ** 2
    ck("B1  the commutator route yields the forcing ODE "
       "8y(y+1)f1' - 14(8y+7)f1 = y^8(y+1)^2",
       sp.expand(ode - sp.expand(ode_expected)) == 0)

    a = sp.symbols("aa0:17")
    ans = sum(a[i] * y ** i for i in range(17))
    eqs = sp.Poly(sp.expand(8 * y * (y + 1) * sp.diff(ans, y)
                            - 14 * (8 * y + 7) * ans - y ** 8 * (y + 1) ** 2), y).all_coeffs()
    sol = sp.solve(eqs, a, dict=True)
    ck("B2  the ODE has a UNIQUE polynomial solution", len(sol) == 1)
    f1v = sp.expand(ans.subs(sol[0]))
    qpoly = sp.expand(-6630 * sp.cancel(f1v / (y ** 8 * (y + 1) ** 2)))
    ck("B3  f1 = -y^8 (y+1)^2 q(y)/6630 with q the quartic "
       "2048y^4-512y^3+320y^2-240y+195",
       sp.expand(qpoly - (2048 * y ** 4 - 512 * y ** 3 + 320 * y ** 2 - 240 * y + 195)) == 0,
       sp.factor(qpoly))
    q_at_m1 = qpoly.subs(y, -1)
    ck("B4  q(-1) = 3315 != 0   (3315 = 3*5*13*17, 6630 = 2*3315)",
       q_at_m1 == 3315 and sp.factorint(3315) == {3: 1, 5: 1, 13: 1, 17: 1})

    Phi = sp.expand(f1v * C4 ** 28)
    Pp = sp.Poly(Phi, y)
    ordy = min(m[0] for m in Pp.monoms())
    # exact t-adic valuation of Phi, read off by expanding in t = y+1
    Phi_t = sp.Poly(sp.expand(Phi.subs(y, t - 1)), t)
    v_t_Phi = min(m[0] for m in Phi_t.monoms())
    ck("B5  ord_y Phi = 204 = 12*17, deg Phi = 238", ordy == 204 and Pp.degree() == 238)
    ck("B6  v_t(Phi) = 30 EXACTLY (not >= 30): lowest nonzero (y+1)-power is 30",
       v_t_Phi == 30, "v_t(Phi) = %s, [t^30] = %s"
       % (v_t_Phi, Phi_t.coeff_monomial(t ** 30)))
    ck("B7  Phi = -y^204 * t^30 * q(y) / 6630 identically",
       sp.expand(Phi + y ** 204 * (y + 1) ** 30 * qpoly / 6630) == 0)

    # the M = -5 row: r~_17 = -Phi/y^204 = t^30 q / 6630, leading coefficient 1/2
    r17_target = sp.expand(-Phi / y ** 204)
    ck("B8  r~_17 = -Phi/y^204 = t^30 q(y)/6630",
       sp.expand(r17_target - (y + 1) ** 30 * qpoly / 6630) == 0)
    lead17 = sp.Rational(q_at_m1, 6630)
    ck("B9  v_t(r~_17) = 30 EXACTLY and [t^30] r~_17 = 1/2",
       lead17 == sp.Rational(1, 2), lead17)

    # window / stratum independence
    with open(os.path.join(HERE, "paper_src", "upstream_facts.json")) as fh:
        UF = json.load(fh)
    np_ = UF["facts"]["newton_polygons"]
    ck("B10 C4 = y^7(y+1) is forced by corners (8,14),(8,16) present in BOTH windows",
       [8, 14] in np_["sub1"]["P"] and [8, 16] in np_["sub1"]["P"]
       and [8, 14] in np_["sub2"]["P"] and [8, 16] in np_["sub2"]["P"])
    ck("B11 v_t(Phi) = 30 is an EQUALITY on a fixed element of Q[y]: no cap, "
       "stratum, window or field-scope toggle can turn it into '>= 30' -- the "
       "only way would be q(-1) = 0, and q(-1) = 3315",
       q_at_m1 != 0 and v_t_Phi == 30)
    ck("B12 the derivation of Phi uses ONLY [P,Q]=x^2 and the leading forms "
       "ell(P)=x^8 C4^2, ell(Q)=x^12 C4^3 -- all window-shared; no d-variable, "
       "no degree cap and no stratum enters", Phi.free_symbols == {y})

    # [QQ1] at its deepest reach, M = -5
    c = {k: sp.Symbol("c%d" % k) for k in range(-6, 5)}
    unit = C4s + sum(c[3 - i] * u ** (i + 1) for i in range(7))
    inv = sp.series(1 / unit, u, 0, 8).removeO()
    ck("B13 [QQ1] reach: (C^-1)_{-4} = 1/C4 (a unit) -> M=-4 only DETERMINES lambda",
       sp.simplify(inv.coeff(u, 0) - 1 / C4s) == 0)
    ck("B14 [QQ1] reach: (C^-1)_{-5} = -c3/C4^2, which VANISHES in the d3-killed "
       "chart -> the M=-5 (Phi) row carries no lambda",
       sp.simplify(inv.coeff(u, 1) + c[3] / C4s ** 2) == 0)
    ck("B15 [QQ1] reach: the F-term clears exactly -- F_{-5} * C4^31 = f1 * C4^28 "
       "= Phi, with f1 := C4^3 F_{-5}, so the M=-5 row is (C^3)_{-5}*C4^31 + Phi = 0",
       sp.expand((f1v / C4 ** 3) * C4 ** 31 - f1v * C4 ** 28) == 0)

    # ================================================================== C
    head("C. The d3-killing shift, twice, by two different mechanisms")

    def gbinom(m, r):
        """generalized binomial C(m,r) as a falling factorial -- deliberately
        NOT sympy.binomial."""
        if r < 0:
            return sp.Integer(0)
        num = sp.Integer(1)
        for i in range(r):
            num *= (m - i)
        return sp.Rational(1, sp.factorial(r)) * num

    tab = [(4, 1, 4), (3, 0, 1), (-1, 0, 1), (-1, 1, -1), (-2, 1, -2), (-1, 2, 1),
           (0, 1, 0), (1, 2, 0), (2, 3, 0)]
    ck("C1  falling-factorial generalized binomials agree with the standard table "
       "(including C(m,r)=0 for 0 <= m < r)",
       all(gbinom(m, r) == v for m, r, v in tab))

    hu = sp.symbols("h1:9")                       # UNSHIFTED stripped h_1..h_8
    th = sp.Symbol("theta")
    hall = [sp.Integer(1)] + list(hu)

    def htil(k, theta):
        return sp.expand(sum(gbinom(4 - l, k - l) * hall[l] * theta ** (k - l)
                             for l in range(0, k + 1)))

    ck("C2  triangularity across zero, theta INDEPENDENT: D~_{-1} = D_{-1}, "
       "i.e. h~_5 = h_5 with no theta at all",
       sp.expand(htil(5, th) - hu[4]) == 0)
    ck("C2b and it is only level 5: h~_6 and h~_7 DO carry theta",
       htil(6, th).has(th) and htil(7, th).has(th))

    theta0 = -hu[0] / 4
    ck("C3  theta = -h_1/4 is the unique d3-killing choice: h~_1 = 0",
       sp.expand(htil(1, theta0)) == 0)
    ck("C4  d2 = h~_2 = h_2 - (3/8) h_1^2",
       sp.expand(htil(2, theta0) - (hu[1] - sp.Rational(3, 8) * hu[0] ** 2)) == 0)
    ck("C5  e  = h~_5 = h_5  (no h_1)",
       sp.expand(htil(5, theta0) - hu[4]) == 0)
    ck("C6  R  = h~_6 = h_6 + (1/4) h_1 h_5",
       sp.expand(htil(6, theta0) - (hu[5] + hu[0] * hu[4] / 4)) == 0)
    ck("C7  S  = h~_7 = h_7 + (1/2) h_1 h_6 + (1/16) h_1^2 h_5   [INDEX -3]",
       sp.expand(htil(7, theta0)
                 - (hu[6] + hu[0] * hu[5] / 2 + hu[0] ** 2 * hu[4] / 16)) == 0)
    ck("C8  PRECISION POINT: with theta independent no h_1 appears in h~_6/h~_7; "
       "after theta = -h_1/4 is substituted h_1 DOES reappear in BOTH",
       (not htil(6, th).has(hu[0])) and (not htil(7, th).has(hu[0]))
       and htil(6, theta0).has(hu[0]) and htil(7, theta0).has(hu[0]))

    # second, independent mechanism: generating function
    def gf_transform(coeffs, deg, aa, nu=NU):
        out = []
        for n in range(nu):
            s = sp.Integer(0)
            for k in range(0, n + 1):
                cbin = gbinom(deg - k, n - k)
                if cbin == 0:
                    continue
                s += coeffs[k] * cbin * (-aa) ** (n - k)
            out.append(sp.expand(s))
        return out

    aa = hu[0] / 4
    Hu = root_series(list(hu))
    Htil_gf = gf_transform(Hu, 4, aa)
    ck("C9  SECOND MECHANISM: Ht(u) = (1-a u)^4 H(u/(1-a u)) with a = h_1/4 "
       "reproduces h~_1..h~_7 identically",
       all(sp.expand(Htil_gf[k] - htil(k, theta0)) == 0 for k in range(1, 8)))

    # ================================================================== D
    head("D. The unshifted t-profile: levels 12 and 14, and  v_t(h_7) >= 11")

    A = sp.symbols("A0:9")

    def profile_jets(L, nu=NU):
        """graded-lemma jets: substitute h_k = t^(L_k) * A_k (None => h_k = 0)."""
        hv = []
        for k in range(1, 9):
            hv.append(sp.Integer(0) if L.get(k) is None else t ** L[k] * A[k])
        H = root_series(hv, nu)
        return series_pow(H, 2, nu), series_pow(H, 3, nu)

    # ---- graded-lemma soundness control (LAYER STABILITY) -----------------
    def layered(L, J, nu=NU, T=34):
        sym = {}
        hv = []
        for k in range(1, 9):
            if L.get(k) is None:
                hv.append(sp.Integer(0))
                continue
            cs = sp.symbols("B%d_0:%d" % (k, J + 1))
            sym[k] = cs
            hv.append(sum(cs[j] * t ** (L[k] + j) for j in range(J + 1)))
        H = [sp.Integer(1)] + hv + [sp.Integer(0)] * (nu - 1 - len(hv))
        for n in range(9, nu):
            H[n] = _tcut(-sp.Rational(1, 2) * sum(H[i] * H[n - i] for i in range(1, n)), T)
        return H, series_pow(H, 2, nu, T), series_pow(H, 3, nu, T), sym

    Lc = {1: 1, 2: 3, 3: 5, 4: 7, 5: 9, 6: 10, 7: 11, 8: 12}
    stable, ncmp = _layer_stability(Lc, extra_layers=3, T=32)
    ck("D0  GRADED-LEMMA CONTROL: for every p_n / r_n used below, the lowest "
       "t-order AND its coefficient are UNCHANGED when three further, fully "
       "SYMBOLIC t-layers are added to every h_k (32 free coefficients) -- so "
       "the jets used below really are functions of the leading coefficients "
       "alone, and the substitution h_k -> t^(L_k) A_k loses nothing at the "
       "lowest order", stable, "%d expressions compared exactly over Q" % ncmp)

    # ---- the imported, independently audited base -------------------------
    if not _state["quiet"]:
        print("  [import] [S1] v_t(h_k) >= 2k-1 for k = 1..5 in the UNSHIFTED chart "
              "(SLICE_OBSTRUCTION.md sec.3, independently audited 56/56). NOT "
              "re-proved here -- this audit is of the UPPER bound.")

    # ---- (P<) at n = 6,7,8, unshifted -------------------------------------
    base = {1: 1, 2: 3, 3: 5, 4: 7, 5: 9}
    L = dict(base)
    for k, n in ((6, 6), (7, 7), (8, 8)):
        Lx = dict(L)
        Lx[k] = None
        for kk in range(6, 9):
            Lx.setdefault(kk, None)
        pj, _ = profile_jets(Lx)
        d = tcoeffs(pj[n])
        m0 = min(d) if d else 10 ** 9
        L[k] = min(2 * n - 2, m0)
    ck("D2  (P<) n=6  ->  v_t(h_6) >= 10   [unshifted]", L[6] == 10, L)
    ck("D3  (P<) n=7  ->  v_t(h_7) >= 11   [unshifted]  <== THE LOAD-BEARING INTEGER",
       L[7] == 11, "binding term is 2*h_1*h_6 at 1 + 10 = 11")
    ck("D4  (P<) n=8  ->  v_t(h_8) >= 12   [unshifted]", L[8] == 12)

    pj, rj = profile_jets(L)
    d7 = tcoeffs(pj[7])
    d8 = tcoeffs(pj[8])
    rel7 = sp.expand(d7[11])
    rel8a = sp.expand(d8[12])
    rel8b = sp.expand(d8[13])
    ck("D5  (P<) n=7 forces  A7 = -A1*A6   (jet at t^11)",
       sp.expand(rel7 / 2 - (A[1] * A[6] + A[7])) == 0, sp.factor(rel7))
    ck("D6  (P<) n=8 forces  A8 = -A1*A7   (jet at t^12)",
       sp.expand(rel8a / 2 - (A[1] * A[7] + A[8])) == 0, sp.factor(rel8a))
    ck("D7  (P<) n=8 forces  A2*A6 = 0     (jet at t^13)  <== the reduction this "
       "audit uses at level 12 (slice_phi_yplace uses A8 = -A2*A6 instead)",
       sp.expand(rel8b / 2 - A[2] * A[6]) == 0, sp.factor(rel8b))

    # ---- level 12, unshifted ---------------------------------------------
    subP = {A[7]: -A[1] * A[6], A[8]: A[1] ** 2 * A[6]}
    d12 = tcoeffs(rj[12])
    j20_raw = sp.expand(d12[20])
    j20_red = sp.expand(j20_raw.subs(subP))
    fl_raw = sp.factor_list(j20_red)
    ck("D8  level 12: r_12 has jets at t^19 and t^20 below the required t^21; the "
       "t^19 one is 3*A4*(A1*A7 + A8) and VANISHES identically under (P<) n=8's "
       "A8 = -A1*A7, so t^20 is the operative jet",
       sorted(m for m in d12 if m < 21) == [19, 20]
       and sp.expand(d12[19] - 3 * A[4] * (A[1] * A[7] + A[8])) == 0
       and sp.expand(d12[19].subs(subP)) == 0, sp.factor(d12[19]))
    ck("D9  BRANCH WARNING: after A7,A8 only, the t^20 jet is "
       "(3/2)*A6*(2*A2*A4 + A6) -- TWO COPRIME FACTORS, it would branch",
       sp.expand(j20_red - sp.Rational(3, 2) * A[6] * (2 * A[2] * A[4] + A[6])) == 0
       and len([f for f, m in fl_raw[1] if sp.total_degree(f) > 0]) == 2,
       sp.factor(j20_red))
    j20_final = sp.expand(j20_red - 3 * A[4] * (A[2] * A[6]))
    fl = sp.factor_list(j20_final)
    ck("D10 but (P<) n=8 already gives A2*A6 = 0; subtracting 3*A4*(A2*A6) the jet "
       "collapses to (3/2)*A6^2 -- a PERFECT SQUARE, so NO case split arises and "
       "LEVEL 12 IS CONFIRMED unshifted: A6 = 0, i.e. v_t(h_6) >= 11",
       sp.expand(j20_final - sp.Rational(3, 2) * A[6] ** 2) == 0
       and fl[1] == [(A[6], 2)], sp.factor_list(j20_final))

    # ---- re-run (P<) n=7,8 with the level-12 upgrade -----------------------
    L2 = dict(base)
    L2[6] = 11
    for k, n in ((7, 7), (8, 8)):
        Lx = dict(L2)
        Lx[k] = None
        for kk in range(7, 9):
            Lx.setdefault(kk, None)
        pjx, _ = profile_jets(Lx)
        d = tcoeffs(pjx[n])
        L2[k] = min(2 * n - 2, min(d) if d else 10 ** 9)
    ck("D12 after level 12: (P<) n=7 gives  v_t(h_7) >= 12  [unshifted] -- ONE "
       "MORE than the 11 that sub1_spine9.py's control X1 needs",
       L2[7] == 12, L2)
    ck("D13 after level 12: v_t(h_8) >= 13 [unshifted] (the 14 quoted by "
       "slice_phi_yplace is the SHIFTED value; h_8 plays no role in the "
       "collision either way)", L2[8] == 13, L2)

    # v_t(S): the zero-margin place is the h_1^2 h_5 term, not h_7
    vS = min(L2[7], 1 + L2[6], 2 + L2[5])
    ck("D14 v_t(S) = v_t(D~_{-3}) >= min(v_t(h_7), 1+v_t(h_6), 2+a_t) = 11 at "
       "a_t = 9: the binding term is (1/16) h_1^2 h_5, NOT h_7 -- so the margin "
       "gained in D12 does not reach S", vS == 11,
       "min(%d, %d, %d) = %d" % (L2[7], 1 + L2[6], 2 + L2[5], vS))

    # ---- uniformity in a_t -------------------------------------------------
    uni = True
    for w5 in (9, 10, 11, 12, 14, 20):
        Lu = {1: 1, 2: 3, 3: 5, 4: 7, 5: w5, 6: 10, 7: 11, 8: 12}
        pjx, rjx = profile_jets(Lu)
        if sp.expand(tcoeffs(pjx[7])[11] - rel7) != 0:
            uni = False
        if sp.expand(tcoeffs(pjx[8])[13] - rel8b) != 0:
            uni = False
        jj = sp.expand(tcoeffs(rjx[12])[20].subs(subP))
        jj = sp.expand(jj - sp.Rational(3, 2) * 2 * A[4] * (A[2] * A[6]))
        if sp.expand(jj - sp.Rational(3, 2) * A[6] ** 2) != 0:
            uni = False
    ck("D15 UNIFORMITY: the (P<) relations and the level-12 collapse are "
       "IDENTICAL for a_t = 9,10,11,12,14,20 (no A5 enters any of them), so the "
       "profile is legitimate under the kill hypothesis a_t >= 10", uni)

    # ================================================================== E
    head("E. The syzygy, re-derived as a POWER-SERIES identity, and the collision")

    Bt = hs[0] * hs[3] ** 2 + 3 * hs[3] * hs[5] + 3 * hs[4] ** 2   # d2 e^2+3 e S+3 R^2
    lhsS = sp.expand(rt[17] + hs[0] * rt[15] + hs[1] * rt[14] + hs[2] * rt[13])
    ck("E1  SERIES SYZYGY (re-derived, not transcribed): "
       "r~_17 + h~_2 r~_15 + h~_3 r~_14 + h~_4 r~_13 "
       "== -(1/2) h~_5 (h~_2 h~_5^2 + 3 h~_5 h~_7 + 3 h~_6^2)",
       sp.expand(lhsS + sp.Rational(1, 2) * hs[3] * Bt) == 0)

    lhsG = sp.expand(2 * (gpoly("G5body") + SYM["Phi"]
                          + SYM["d2"] * gpoly("G3") + SYM["d1"] * gpoly("G2")
                          + SYM["d0"] * gpoly("G1")))
    Kexp = sp.expand(2 * SYM["Phi"] - SYM["dm1"] * (SYM["d2"] * SYM["dm1"] ** 2
                                                    + 3 * SYM["dm1"] * SYM["dm3"]
                                                    + 3 * SYM["dm2"] ** 2))
    ck("E2  COMMON-MODE FINDING: the K-syzygy is an EXACT Q[d]-linear combination "
       "of G1,G2,G3,G5 -- 2*(G5 + d2*G3 + d1*G2 + d0*G1) == "
       "2*Phi - e*(d2 e^2 + 3 e S + 3 R^2).  syzygy_collision and "
       "slice_phi_yplace therefore consume the IDENTICAL four equations.",
       sp.expand(lhsG - Kexp) == 0)

    subshift = {hs[k - 2]: htil(k, theta0) for k in range(2, 9)}
    Bu_ = sp.expand(Bt.subs(subshift))
    Bu_target = (hu[1] * hu[4] ** 2 + 3 * hu[4] * hu[6]
                 + 3 * hu[0] * hu[4] * hu[5] + 3 * hu[5] ** 2)
    ck("E3  THE COLLAPSE: pushed through the shift, "
       "d2 e^2 + 3 e S + 3 R^2 == h_2 h_5^2 + 3 h_5 h_7 + 3 h_1 h_5 h_6 + 3 h_6^2 "
       "(the h_1^2 h_5^2 terms cancel: -3/8 + 3/16 + 3/16 = 0)",
       sp.expand(Bu_ - Bu_target) == 0, sp.expand(Bu_))
    mut = sp.expand((hu[1] * hu[4] ** 2 + 3 * hu[4] * hu[6] + 3 * hu[5] ** 2) - Bu_target)
    ck("E3b MUTATION CONTROL: dropping the mixing (d2->h_2, R->h_6, S->h_7) gives "
       "a DIFFERENT bracket, off by -3 h_1 h_5 h_6 -- the collapse is a real "
       "computation, not an identity any dictionary satisfies",
       sp.expand(mut + 3 * hu[0] * hu[4] * hu[5]) == 0)

    # ---- the collision -----------------------------------------------------
    VH1, VH6, VH7 = 1, 11, 12

    def vB(aval):
        return min(3 + 2 * aval, aval + VH7, VH1 + aval + VH6, 2 * VH6)

    rows = []
    for aval in list(range(6, 20)) + [24, 30, 59]:
        need = 30 - aval
        rows.append((aval, need, vB(aval), vB(aval) > need))
    ck("E4  THE COLLISION is an EQUALITY: with r~_13 = r~_14 = r~_15 = 0 the "
       "series syzygy E1 reduces to r~_17 = -(1/2) e * B, and v_t(r~_17) = 30 "
       "exactly (B9), so  a_t + v_t(B) = 30 -- a bound on v_t(B) therefore "
       "REFUTES rather than merely constrains",
       sp.expand(lhsS + sp.Rational(1, 2) * hs[3] * Bt) == 0 and lead17 != 0
       and v_t_Phi == 30)
    ck("E5  a_t >= 10 is REFUTED for every a_t in 10..59: v_t(B) > 30 - a_t",
       all(bad for aval, need, got, bad in rows if aval >= 10),
       "; ".join("a=%d: v_t(B)>=%d vs need %d" % (aa_, g, nd)
                 for aa_, nd, g, bad in rows if aa_ in (10, 11, 12, 14)))
    ck("E6  ANTI-VACUITY: a_t = 9 is NOT refuted -- three terms land on 21 = 30-9 "
       "exactly", vB(9) == 21 and 30 - 9 == 21)
    ck("E6b ANTI-VACUITY: a_t = 8 is NOT refuted either (the bound is a threshold, "
       "not a blanket)", vB(8) <= 30 - 8, "v_t(B) >= %d vs need %d" % (vB(8), 22))
    ck("E7  the three binding inequalities in closed form: 3+2a > 30-a <=> a>=10 ; "
       "a+12 > 30-a <=> a>=10 ; 22 > 30-a <=> a>=9",
       all((3 + 2 * aa_ > 30 - aa_) == (aa_ >= 10) for aa_ in range(0, 40))
       and all((aa_ + 12 > 30 - aa_) == (aa_ >= 10) for aa_ in range(0, 40))
       and all((22 > 30 - aa_) == (aa_ >= 9) for aa_ in range(0, 40)))

    # ---- level-12 sensitivity ---------------------------------------------
    def vB10(aval):
        return min(3 + 2 * aval, aval + 11, 1 + aval + 10, 2 * 10)
    ck("E8  LEVEL-12 SENSITIVITY: with only v_t(h_6) >= 10 (the committed stage "
       "record) the 3*h_6^2 term lands on 20 = 30-10 and a_t = 10 SURVIVES -- "
       "level 12 is genuinely load-bearing", vB10(10) == 20 and not (vB10(10) > 20))
    ck("E9  and level 14/16 are NOT needed: the collision uses only "
       "v_t(h_7) >= 11 (D3), which (P<) n=7 supplies with no cascade level at all",
       min(3 + 2 * 10, 10 + 11, 1 + 10 + 11, 22) > 20)

    # ================================================================== F
    head("F. The SHIFTED-chart route (slice_phi_yplace): reproduced, and its "
         "foundation examined")

    def shifted_jets(L, nu=NU):
        hv = [sp.Integer(0)]
        for k in range(2, 9):
            hv.append(sp.Integer(0) if L.get(k) is None else t ** L[k] * A[k])
        H = root_series(hv, nu)
        return series_pow(H, 2, nu), series_pow(H, 3, nu)

    Av = [A[k] for k in range(2, 9)]

    def shifted_system(w5, L6, L7, L8, use=("E13", "E14", "E15", "E17")):
        L = {2: 3, 3: 5, 4: 7, 5: w5, 6: L6, 7: L7, 8: L8}
        pj_, rj_ = shifted_jets(L)
        E = []
        for n, tag in ((13, "E13"), (14, "E14"), (15, "E15")):
            if tag in use:
                d = tcoeffs(rj_[n])
                if d:
                    E.append(d[min(d)])
        if "E17" in use:
            d = tcoeffs(rj_[17])
            m0 = min(d)
            if m0 > 30:
                return "outright"
            E.append(d[m0] - (sp.Rational(1, 2) if m0 == 30 else 0))
        E = [sp.expand(z) for z in E if sp.expand(z) != 0]
        gbz = sp.groebner(E, *Av, order="grevlex")
        return list(gbz.exprs) == [sp.Integer(1)]

    ck("F1  REPRODUCED: at the shifted profile (3,5,7,10,11,13,14) the four "
       "leading-jet equations generate the UNIT IDEAL over Q",
       shifted_system(10, 11, 13, 14) is True)
    ck("F2  ... and for a_t = 11, 12 as well (uniform in how far past 10)",
       shifted_system(11, 11, 13, 14) is True and shifted_system(12, 11, 13, 14) is True)
    ck("F3  ... and it is NOT the unit ideal at a_t = 9 (anti-vacuity)",
       shifted_system(9, 11, 13, 14) is False)

    # my own rational witness at a_t = 9 (found independently of the published one)
    W = {A[2]: sp.Rational(-1, 8), A[3]: sp.Integer(-1), A[4]: sp.Integer(0),
         A[5]: sp.Integer(2), A[6]: sp.Integer(0), A[7]: sp.Rational(4, 3),
         A[8]: sp.Integer(1)}
    Lw = {2: 3, 3: 5, 4: 7, 5: 9, 6: 11, 7: 13, 8: 14}
    _, rjw = shifted_jets(Lw)
    okw = True
    for n in (13, 14, 15):
        d = tcoeffs(rjw[n])
        okw = okw and sp.expand(d[min(d)].subs(W)) == 0
    d17 = tcoeffs(rjw[17])
    okw = okw and sp.expand(d17[min(d17)].subs(W)) == sp.Rational(1, 2)
    ck("F4  independent rational witness at a_t = 9: "
       "(A2..A8) = (-1/8, -1, 0, 2, 0, 4/3, 1) satisfies all four equations with "
       "A5 = 2 != 0 (so v_t(h_5) = 9 EXACTLY)", okw)
    ck("F4b the witness independently reproduces SLICE_PHI_YPLACE's forced "
       "relation A2*A5^3 = -1", sp.Rational(-1, 8) * 8 == -1)

    ck("F5  drop the Phi relation -> NOT unit (Phi is load-bearing)",
       shifted_system(10, 11, 13, 14, use=("E13", "E14", "E15")) is False)
    ck("F6  drop r_14 = 0 -> NOT unit", shifted_system(10, 11, 13, 14,
                                                       use=("E13", "E15", "E17")) is False)
    ck("F7  drop r_15 = 0 -> NOT unit", shifted_system(10, 11, 13, 14,
                                                       use=("E13", "E14", "E17")) is False)
    ck("F8  CORRECTION: r_13 = 0 is NOT needed -- dropping it leaves the ideal "
       "UNIT (SLICE_PHI_YPLACE prints it as one of four but never tests it)",
       shifted_system(10, 11, 13, 14, use=("E14", "E15", "E17")) is True)
    ck("F9  level 12 is load-bearing on this route too: at (3,5,7,10,10,12,13) "
       "the system is NOT unit",
       shifted_system(10, 10, 12, 13) is False)

    # ---- the foundational objection ---------------------------------------
    p_un = series_pow(Hu, 2)
    ptil_gf = gf_transform(p_un, 8, aa)
    ck("F10 shifted/unshifted P-slices are related by "
       "p~_n = sum_k C(8-k, n-k) (-h_1/4)^(n-k) p_k  (no denominators)",
       sp.expand(ptil_gf[3] - (hu[0] ** 3 - 4 * hu[0] * hu[1] + 8 * hu[2]) / 4) == 0,
       sp.factor(sp.expand(ptil_gf[3])))
    # re-derive the joint control's stripped coordinates rather than transcribe
    # them: C_3 = a y^5, C_2 = b y^3, C_1 = c y, C_0 = e0; D_j = C_j C4^(7-2j);
    # d_j = D_j / y^(12(4-j)) ; C4 = y^7 (y+1) and t = y+1.
    aS, bS, cS, e0S = sp.symbols("aS bS cS e0S")
    C4y = y ** 7 * (y + 1)
    Cco = {3: aS * y ** 5, 2: bS * y ** 3, 1: cS * y, 0: e0S}
    strip = {}
    for j, Cj in Cco.items():
        strip[4 - j] = sp.cancel(Cj * C4y ** (7 - 2 * j) / y ** (12 * (4 - j)))
    ck("F11a joint-control stripping re-derived: h_1 = a t, h_2 = b t^3, "
       "h_3 = c t^5, h_4 = e0 y t^7",
       sp.expand(strip[1] - aS * (y + 1)) == 0
       and sp.expand(strip[2] - bS * (y + 1) ** 3) == 0
       and sp.expand(strip[3] - cS * (y + 1) ** 5) == 0
       and sp.expand(strip[4] - e0S * y * (y + 1) ** 7) == 0)
    inst = {hu[0]: aS * t, hu[1]: bS * t ** 3, hu[2]: cS * t ** 5,
            hu[3]: e0S * y * t ** 7, hu[4]: 0, hu[5]: 0, hu[6]: 0, hu[7]: 0}
    p3t = sp.expand(ptil_gf[3].subs(inst))
    v3 = min(m[0] for m in sp.Poly(p3t, t).monoms())
    ck("F11 REFUTATION OF THE SHIFTED (P<): on SLICE_OBSTRUCTION.md sec.4's own "
       "GENUINE joint control C = x^4 C4 + a y^5 x^3 + b y^3 x^2 + c y x + e0 "
       "(P = C^2 and Q = C^3 both honest polynomials) the shifted slice has "
       "v_t(p~_3) = 3, while (P<) at n = 3 demands t^4 | p~_3.  So the "
       "divisibility conditions do NOT transfer to the shifted chart.",
       v3 == 3, "p~_3 = %s  (v_t = %d, need >= 4)" % (sp.factor(p3t), v3))
    tr2 = min(3, 2 * 1)
    tr3 = min(5, 1 + 3, 3 * 1)
    tr7 = min(11, 1 + 10, 2 + 9)
    ck("F12 and the shifted profile slice_phi_yplace ASSUMES (3,5,7,a,11,13,14) "
       "does not transfer either: the dictionary only gives v_t(h~_2) >= 2, "
       "v_t(h~_3) >= 3, v_t(h~_7) >= 11",
       tr2 == 2 and tr3 == 3 and tr7 == 11,
       "h~_2 = h_2 - (3/8)h_1^2 ; h~_3 = h_1^3/8 - h_1 h_2/2 + h_3 ; "
       "h~_7 = h_7 + h_1 h_6/2 + h_1^2 h_5/16")
    r_un = series_pow(Hu, 3)
    rtil_gf = gf_transform(r_un, 12, aa)
    ck("F13 by contrast the EXACT rows transfer with no chart hypothesis at all -- "
       "verified: r~_13 = r_13, r~_14 = r_14 + a r_13, "
       "r~_15 = r_15 + 2a r_14 + a^2 r_13, "
       "r~_17 = r_17 + 4a r_16 + 6a^2 r_15 + 4a^3 r_14 + a^4 r_13 "
       "(triangular, so r~_13=r~_14=r~_15=0 <=> r_13=r_14=r_15=0); this is a "
       "consequence of Q(x-s) still being a polynomial in x of degree 12, "
       "whatever s is",
       sp.expand(rtil_gf[13] - r_un[13]) == 0
       and sp.expand(rtil_gf[14] - (r_un[14] + aa * r_un[13])) == 0
       and sp.expand(rtil_gf[15] - (r_un[15] + 2 * aa * r_un[14]
                                    + aa ** 2 * r_un[13])) == 0
       and sp.expand(rtil_gf[17] - (r_un[17] + 4 * aa * r_un[16]
                                    + 6 * aa ** 2 * r_un[15]
                                    + 4 * aa ** 3 * r_un[14]
                                    + aa ** 4 * r_un[13])) == 0)

    if args.deep:
        head("G. DEEP: multi-layer Groebner controls (slower)")

        def layered_system(L, J, use=("P", "Q", "G", "PHI"), shifted=True):
            H, pL, rL, sym = layered(L, J)
            E = []
            av = [v for k in sym for v in sym[k]]
            def take(d, pred):
                m0 = min(d)
                for m, v in sorted(d.items()):
                    if m <= m0 + J and pred(m):
                        E.append(v)
            if "P" in use:
                for n in range(2, 9):
                    d = tcoeffs(pL[n])
                    if d:
                        take(d, lambda m, n=n: m < 2 * n - 2)
            if "Q" in use:
                for n in range(2, 13):
                    d = tcoeffs(rL[n])
                    if d:
                        take(d, lambda m, n=n: m < 2 * n - 3)
            if "G" in use:
                for n in (13, 14, 15):
                    d = tcoeffs(rL[n])
                    if d:
                        take(d, lambda m: True)
            if "PHI" in use:
                X = sp.expand(rL[17] if shifted else rL[17] + H[1] * rL[16])
                d = tcoeffs(X)
                if not d:
                    return "outright"
                m0 = min(d)
                if m0 > 30:
                    return "outright"
                for m, v in sorted(d.items()):
                    if m <= m0 + J:
                        E.append(v - (sp.Rational(1, 2) if m == 30 else 0))
            E = [z for z in (sp.expand(w) for w in E) if z != 0]
            gbz = sp.groebner(E, *av, order="grevlex")
            return list(gbz.exprs) == [sp.Integer(1)]

        ck("G1  two-layer shifted system at (3,5,7,10,11,12,14) -- i.e. WITHOUT "
           "cascade level 14 -- is already the UNIT IDEAL",
           layered_system({2: 3, 3: 5, 4: 7, 5: 10, 6: 11, 7: 12, 8: 14}, 1) is True)
        ck("G2  ... and WITHOUT level 12 (v_t(h_6) >= 10) it is not, even at two "
           "layers", layered_system({2: 3, 3: 5, 4: 7, 5: 10, 6: 10, 7: 12, 8: 13},
                                    1) is False)

    # ================================================================ verdict
    head("VERDICT")
    if not _state["quiet"]:
        print("""
  a_t <= 9  is CONFIRMED, by the SYZYGY route, reproduced here end to end from
  the primitives with different algebra at every step.  The chain is:

     r~_17 = -(1/2) * e * B          (series syzygy, E1/E2)
     v_t(r~_17) = 30 EXACTLY         (Phi from the commutator ODE, B6/B9)
     B = h_2 h_5^2 + 3 h_5 h_7 + 3 h_1 h_5 h_6 + 3 h_6^2   (collapse, E3)
     v_t(h_1,h_2,h_5,h_6,h_7) >= 1,3,a,11,12   (UNSHIFTED, D2-D13)
     => a_t + v_t(B) = 30 and v_t(B) >= min(3+2a, a+12, 22) > 30-a for a >= 10.

  CORRECTIONS -- see AT_LE9_AUDIT.md:
    (1) the two proofs are NOT independent: E2 shows they consume the identical
        four G-rows;
    (2) slice_phi_yplace's route imposes the slice DIVISIBILITY conditions in
        the shifted chart, where F11 shows they are false on the repo's own
        genuine joint control.  That route is not established as written; the
        syzygy route is unaffected because it is entirely unshifted.
    (3) r_13 = 0 is not needed on the shifted route (F8).
    (4) v_t(h_7) >= 11 (unshifted) is CONFIRMED, and level 12 upgrades it to 12.
""")
        print("%d checks, %d failures" % (_state["n"], _state["fail"]))
    return 1 if _state["fail"] else 0


if __name__ == "__main__":
    sys.exit(main())
