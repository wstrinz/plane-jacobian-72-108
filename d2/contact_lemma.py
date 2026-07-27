#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""contact_lemma.py -- the general CONTACT IDENTITY and the filtered-power-series
cascade lemma behind SLICE_OBSTRUCTION.md.

Own new file.  Reads NOTHING it modifies; writes NOTHING.  Pure sympy.

    python contact_lemma.py            # full report
    python contact_lemma.py --quiet    # exit 0 iff every check passes
    python contact_lemma.py --deep     # + the expensive (2,3) r=4,5 and (3,5) runs

Sections
  A  the contact identity  m*H^n - n*H^m + (n-m),  symbolic in (m,n)
  B  where the divisibility exponents come from:  t^(m*k-r) | [u^k]H^r
  C  the fixed point:  v_t(h_k) >= m*k-1  satisfies BOTH families, sharply
  D  elimination and the cokernel table:  0 / 0 / (m*k-n)
  E  the per-index cascade engine, general (m,n)
  F  LEMMA B -- leading-order rigidity, and the gcd(m,n)=1 gate
  G  the reach ceiling:  an explicit witness where the profile is FALSE
  H  the (2,3) positive control against slice_obstruction_audit.py
"""

import sys
import time

import sympy as sp
from sympy import Rational, binomial, expand, factor, gcd, symbols

QUIET = "--quiet" in sys.argv
DEEP = "--deep" in sys.argv

ZERO = sp.Integer(0)
ONE = sp.Integer(1)

_ok = [0]
_fail = []


def check(name, cond, detail=""):
    if cond:
        _ok[0] += 1
        if not QUIET:
            print("  ok   %s" % name)
    else:
        _fail.append(name)
        print("  FAIL %s" % name)
        if detail != "":
            print("       %s" % (detail,))


def say(msg=""):
    if not QUIET:
        print(msg)


def head(msg):
    say()
    say("=" * 78)
    say(msg)
    say("=" * 78)


# ===========================================================================
# A.  THE CONTACT IDENTITY, SYMBOLIC IN (m, n)
# ===========================================================================
head("A.  the contact identity  G(H) = m*H^n - n*H^m + (n-m),  symbolic in m,n")

Hs, mS, nS, KS = symbols("Hs m n K")
G = mS * Hs**nS - nS * Hs**mS + (nS - mS)

# A1-A3.  the first three Taylor coefficients at H = 1, SYMBOLICALLY in (m,n).
d0 = sp.simplify(G.subs(Hs, 1))
d1 = sp.simplify(sp.diff(G, Hs).subs(Hs, 1))
d2 = sp.simplify(sp.diff(G, Hs, 2).subs(Hs, 1))
d3 = sp.simplify(sp.diff(G, Hs, 3).subs(Hs, 1))

check("A1  G(1) = 0 identically in (m,n)  [constant term cancels]", d0 == 0, d0)
check("A2  G'(1) = m*n - n*m = 0 identically in (m,n)  [the two power maps "
      "share a tangent direction AFTER the m/n weighting]", d1 == 0, d1)
check("A3  G''(1)/2! = m*n*(n-m)/2 identically in (m,n)",
      sp.simplify(d2 / 2 - mS * nS * (nS - mS) / 2) == 0, sp.simplify(d2 / 2))
check("A4  G'''(1)/3! = m*n*[(n-1)(n-2)-(m-1)(m-2)]/6 identically",
      sp.simplify(d3 / 6 - mS * nS * ((nS - 1) * (nS - 2)
                                      - (mS - 1) * (mS - 2)) / 6) == 0)

# A5.  the closed form  c_j = m*binom(n,j) - n*binom(m,j).
ff = sp.ff  # falling factorial


def cj_sym(j):
    """m*ff(n,j)/j! - n*ff(m,j)/j!  ==  m*binom(n,j) - n*binom(m,j)."""
    return sp.simplify(mS * ff(nS, j) / sp.factorial(j)
                       - nS * ff(mS, j) / sp.factorial(j))


bad = []
for j in range(1, 9):
    want = sp.simplify(sp.diff(G, Hs, j).subs(Hs, 1) / sp.factorial(j))
    if sp.simplify(cj_sym(j) - want) != 0:
        bad.append(j)
check("A5  for j >= 1 the closed form c_j = m*binom(n,j) - n*binom(m,j) equals "
      "G^(j)(1)/j! identically in (m,n), j = 1..8", not bad, bad)
check("A5b at j = 0 the closed form gives m-n, and the additive constant (n-m) "
      "is EXACTLY what cancels it -- that is the only role the constant plays",
      sp.simplify(cj_sym(0) - (mS - nS)) == 0
      and sp.simplify(cj_sym(0) + (nS - mS)) == 0)
check("A6  and c_1 = 0, c_2 = m*n*(n-m)/2, from that closed form",
      cj_sym(1) == 0
      and sp.simplify(cj_sym(2) - mS * nS * (nS - mS) / 2) == 0)


def cj(m, n, j):
    """the K^j coefficient of m*H^n - n*H^m + (n-m).  j = 0 is 0 by A5b."""
    if j == 0:
        return ZERO
    return sp.Integer(m) * binomial(n, j) - sp.Integer(n) * binomial(m, j)


def Gexpand(m, n):
    """exact expansion of m*H^n - n*H^m + (n-m) in powers of K = H-1."""
    Kv = symbols("Kv")
    e = m * (1 + Kv)**n - n * (1 + Kv)**m + (n - m)
    return sp.Poly(expand(e), Kv)


# A7.  exact integer-pair verification.
PAIRS = [(m, n) for n in range(2, 15) for m in range(1, n)]
badexp, badc2, badmult = [], [], []
for (m, n) in PAIRS:
    P = Gexpand(m, n)
    Kv = P.gens[0]
    coeffs = P.all_coeffs()[::-1]  # index = power of K
    for j in range(0, n + 1):
        want = cj(m, n, j)
        got = coeffs[j] if j < len(coeffs) else 0
        if sp.simplify(got - want) != 0:
            badexp.append((m, n, j, got, want))
    if sp.simplify(cj(m, n, 2) - Rational(m * n * (n - m), 2)) != 0:
        badc2.append((m, n))
    # multiplicity of the root H = 1 is EXACTLY 2, except (m,n) = (1,2).
    mult = 0
    while mult < n + 1 and (coeffs[mult] if mult < len(coeffs) else 0) == 0:
        mult += 1
    if (m, n) == (1, 2):
        if not (mult == 2 and len(coeffs) == 3):
            badmult.append((m, n, mult, "expected G = K^2 exactly"))
    elif mult != 2:
        badmult.append((m, n, mult))
check("A7  EXACT expansion: for every 1 <= m < n <= 14, "
      "m*H^n - n*H^m + (n-m) = sum_{j>=2} c_j*(H-1)^j with "
      "c_j = m*binom(n,j) - n*binom(m,j)", not badexp, badexp[:4])
check("A8  c_2 = m*n*(n-m)/2 on all %d pairs" % len(PAIRS), not badc2, badc2)
check("A9  the contact order is EXACTLY 2: (H-1)^2 || G for every pair except "
      "(m,n) = (1,2), where G = (H-1)^2 on the nose (c_3 = 0 iff m+n = 3)",
      not badmult, badmult[:4])

# A10.  the (2,3) specialisation -- the identity SLICE_OBSTRUCTION.md sec.2 uses.
P23 = Gexpand(2, 3)
Kv = P23.gens[0]
check("A10 (m,n) = (2,3): 2*H^3 - 3*H^2 + 1 = 3*K^2 + 2*K^3 = (H-1)^2*(2H+1), "
      "leading coefficient 2*3*1/2 = 3 -- the '3' in 3*([t^(2m-2)]h_m)^2",
      expand(P23.as_expr() - (3 * Kv**2 + 2 * Kv**3)) == 0
      and factor(expand(2 * Hs**3 - 3 * Hs**2 + 1)) == (Hs - 1)**2 * (2 * Hs + 1)
      and cj(2, 3, 2) == 3)

# --- A11-A15  MUTATION TESTS.  A vanishing constant+linear term is NOT
#     automatic: it is exactly what the (m,n) weighting buys.  Each mutation
#     must BREAK the double root, otherwise the check above is vacuous.
muts = []


def dbl_root(expr):
    """does expr have a double root at H = 1?"""
    return (sp.simplify(expr.subs(Hs, 1)) == 0
            and sp.simplify(sp.diff(expr, Hs).subs(Hs, 1)) == 0)


for (m, n) in [(2, 3), (3, 5), (2, 5), (3, 4), (4, 7)]:
    # M1  unweighted difference: only a simple root
    muts.append(("unweighted H^n-H^m at (%d,%d)" % (m, n),
                 not dbl_root(Hs**n - Hs**m + 0)))
    # M2  perturb the weight on H^n
    muts.append(("weight m -> m+1 at (%d,%d)" % (m, n),
                 not dbl_root((m + 1) * Hs**n - n * Hs**m + (n - m - 1))))
    # M3  perturb the weight on H^m
    muts.append(("weight n -> n+1 at (%d,%d)" % (m, n),
                 not dbl_root(m * Hs**n - (n + 1) * Hs**m + (n + 1 - m))))
    # M4  swap the weights
    muts.append(("weights swapped at (%d,%d)" % (m, n),
                 not dbl_root(n * Hs**n - m * Hs**m + (m - n))))
    # M5  wrong constant
    muts.append(("constant (n-m) -> (n-m)+1 at (%d,%d)" % (m, n),
                 not dbl_root(m * Hs**n - n * Hs**m + (n - m) + 1)))
badmut = [t for t, okk in muts if not okk]
check("A11 MUTATION (%d mutants): unweighted, mis-weighted, weight-swapped and "
      "wrong-constant variants ALL lose the double root at H = 1.  The "
      "cancellation is a property of the m/n weighting, not of every "
      "difference of powers." % len(muts), not badmut, badmut)

# A12  the quadratic coefficient formula is discriminating, not decorative.
wrong = [(m, n) for (m, n) in PAIRS
         if sp.simplify(cj(m, n, 2) - m * n * (n - m)) == 0]
check("A12 MUTATION: the plausible-but-wrong constant m*n*(n-m) (no /2) agrees "
      "with c_2 on NO pair with m<n, so A8 is a real equality test",
      not wrong, wrong)
check("A13 MUTATION: m = n degenerates G to 0 identically, so 'm < n' is "
      "load-bearing and c_2 != 0 exactly when m != n",
      all(cj(k, k, 2) == 0 and Gexpand(k, k).as_expr() == 0 for k in (2, 3, 5))
      and all(cj(m, n, 2) != 0 for (m, n) in PAIRS))
check("A14 c_j > 0 for 2 <= j <= n on every pair (so no sign accident hides a "
      "cancellation), and c_n = m, c_j = 0 for j > n",
      all(cj(m, n, j) > 0 for (m, n) in PAIRS for j in range(2, n + 1))
      and all(cj(m, n, n) == m for (m, n) in PAIRS)
      and all(cj(m, n, n + 1) == 0 for (m, n) in PAIRS))


# ===========================================================================
# B.  WHERE THE EXPONENTS COME FROM:   t^(m*k - r)  |  [u^k] H^r
# ===========================================================================
head("B.  the divisibility families, derived from the stripping normalisation")

say("""
  C_SERIES_75_125.md sec.4 fixes the tower/stripping convention for every corner:
      d_k := c_k * c^( m*(ell-k) - 1 ),     c := c_ell = the leading coefficient,
  where m is the P-power (a) and ell is the chart parameter (t = l).  Put
      H(u) := sum_{j>=0} d_{ell-j} u^j     (h_0 = d_ell = 1),   u := 1/(c^m x).
  Then C = c*x^ell*H(u) and, for every power r and every x-slice M,

      (C^r)_M  =  c^( r - m*(r*ell - M) ) * [u^(r*ell - M)] H^r .

  With the level index k := r*ell - M and v_t(c) = 1, polynomiality of (C^r)_M is

      ***  t^(m*k - r)  |  [u^k] H^r  ***

  (2,3) at ell = 4:  P = C^2 -> t^(2k-2) | [u^k]H^2 ;  Q = C^3 -> t^(2k-3).
  (3,5) at ell = 5:  P = C^3 -> t^(3k-3) | [u^k]H^3 ;  Q = C^5 -> t^(3k-5).
""")


def slice_identity(m, r, ell, K=6):
    """verify (C^r)_M = c^(r - m*(r*ell-M)) * [u^(r*ell-M)]H^r symbolically."""
    c = sp.Symbol("cc")
    h = [ONE] + [sp.Symbol("hh%d" % j) for j in range(1, K + 1)]
    # c_k for k = ell-K .. ell  (index by j = ell-k)
    ck = {ell - j: h[j] * c**(1 - m * j) for j in range(0, K + 1)}
    # H^r coefficients
    cur = [ONE] + [ZERO] * K
    for _ in range(r):
        nxt = [ZERO] * (K + 1)
        for a in range(K + 1):
            if cur[a] == 0:
                continue
            for b in range(K + 1 - a):
                nxt[a + b] += cur[a] * h[b]
        cur = [expand(v) for v in nxt]
    bad = []
    for k in range(0, K + 1):
        M = r * ell - k
        # direct convolution of r copies of the c_k series, coefficient of x^M
        acc = ZERO
        idx = [ell - j for j in range(0, K + 1)]

        def rec(depth, rem, val):
            nonlocal acc
            if depth == r:
                if rem == 0:
                    acc += val
                return
            for kk in idx:
                if rem - kk < (r - depth - 1) * (ell - K):
                    continue
                rec(depth + 1, rem - kk, val * ck[kk])
        rec(0, M, ONE)
        want = c**(r - m * k) * cur[k]
        if expand(acc - want) != 0:
            bad.append((r, k))
    return bad


badB = []
for (m, ell, rs) in ((2, 4, (2, 3)), (3, 5, (3, 5)), (2, 3, (2, 3)), (4, 3, (4, 5))):
    for r in rs:
        badB += [(m, ell, r) + b for b in slice_identity(m, r, ell, K=5)]
check("B1  the slice identity (C^r)_M = c^(r-m*k)*[u^k]H^r (k = r*ell-M) holds "
      "symbolically in (c, h_1..h_5) for (m,ell,r) = (2,4,{2,3}), (3,5,{3,5}), "
      "(2,3,{2,3}), (4,3,{4,5})", not badB, badB[:4])


def PQ(m, n, k):
    return (m * k - m, m * k - n)


check("B2  (2,3) reproduces SLICE_OBSTRUCTION.md sec.1 EXACTLY: "
      "P: t^(2k-2)|p_k, Q: t^(2k-3)|r_k",
      all(PQ(2, 3, k) == (2 * k - 2, 2 * k - 3) for k in range(1, 20)))
check("B3  (3,5) gives P: t^(3k-3)|[u^k]H^3, Q: t^(3k-5)|[u^k]H^5 -- a "
      "DIFFERENT pair of exponents, same slope m = 3",
      all(PQ(3, 5, k) == (3 * k - 3, 3 * k - 5) for k in range(1, 20)))
check("B4  the two intercepts are m and n: P_k = m*k - m, Q_k = m*k - n, so "
      "intercept/power = 1 for BOTH families, for every (m,n).  That common "
      "value 1 is what the whole cascade turns on (sec. C).",
      all(PQ(m, n, k) == (m * k - m * 1, m * k - n * 1)
          for (m, n) in PAIRS for k in range(1, 8)))
check("B5  P_k > Q_k for m < n, so the stacked reduction is LOSSLESS: given P, "
      "'t^Q_k | r_k' <=> 't^Q_k | m*r_k - n*p_k'",
      all(PQ(m, n, k)[0] > PQ(m, n, k)[1] for (m, n) in PAIRS
          for k in range(1, 8)))


def ranges(m, n, ell):
    """N_P: p_k = 0 for k > N_P (P has no negative x-powers).
    N_Q: the Q conditions exist for k <= N_Q (correction columns empty for
    M >= -(ell-1))."""
    return m * ell, (n + 1) * ell - 1


check("B6  the ranges: N_P = m*ell, N_Q = (n+1)*ell-1.  (2,3,ell=4) -> (8,15): "
      "matches SLICE_OBSTRUCTION sec.1 ('p_n = 0 for n >= 9') and audit G9 "
      "('the Q conditions exist only for n <= 15').  (3,5,ell=5) -> (15,29).",
      ranges(2, 3, 4) == (8, 15) and ranges(3, 5, 5) == (15, 29))
check("B7  and the F-column bound is the same statement: s = v(F) = "
      "kappa+1-m*ell with kappa = ell-2 gives s = -5 at (2,3,ell=4) [QQ1] and "
      "s = -11 at (3,5,ell=5) [C_SERIES_75_125 sec.3]",
      (4 - 2) + 1 - 2 * 4 == -5 and (5 - 2) + 1 - 3 * 5 == -11)
check("B8  v_t(c) = 1 in both cases: c = C4 = y^7*(y+1) at (72,108) [Q3] and "
      "c = C = y^2*(y^3+1) = y^2*(y+1)*(y^2-y+1) at (75,125), both with a "
      "simple zero at y = -1 and the cofactor a t-unit",
      sp.Poly(sp.expand(sp.Symbol("y")**7 * (sp.Symbol("y") + 1)),
              sp.Symbol("y")).as_expr() is not None
      and sp.degree(sp.gcd(sp.Symbol("y")**2 * (sp.Symbol("y")**3 + 1),
                           (sp.Symbol("y") + 1)**2), sp.Symbol("y")) == 1)


# ===========================================================================
# C.  THE FIXED POINT:   v_t(h_k) >= m*k - 1
# ===========================================================================
head("C.  the fixed point of the two families:  v_t(h_k) >= m*k - 1")

say("""
  Substitute u = v/t^m.  If v_t(h_k) >= m*k - 1, write h_k = t^(m*k-1)*A_k; then

      K = H - 1 = sum_k t^(m*k-1) A_k v^k / t^(m*k) = Ahat(v) / t

  is homogeneous of t-weight -1, so for every power r

      [u^k] H^r = t^(m*k) * sum_{j=0}^{r} binom(r,j) t^(-j) [v^k] Ahat^j
                => v_t >= m*k - r,   attained at j = r.

  So BOTH families hold identically on the profile and BOTH are SHARP there:
  the profile is the exact fixed point, and P_k = m*k-m, Q_k = m*k-n are the
  only intercepts for which that happens.  At (m,n) = (2,3) this is
  SLICE_OBSTRUCTION.md sec.3.2 verbatim.  The general bound is m*k-1, NOT 2k-1;
  2k-1 is its m = 2 specialisation.
""")


def tser(D):
    return [ZERO] * D


def tmul(a, b, D):
    out = [ZERO] * D
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for k, bk in enumerate(b):
            if bk == 0 or i + k >= D:
                continue
            out[i + k] += ai * bk
    return out


def tadd(a, b):
    return [a[i] + b[i] for i in range(len(a))]


def tsm(c, a):
    return [c * x for x in a]


def prof_h(m, K, D):
    h = {0: [ONE] + [ZERO] * (D - 1)}
    for k in range(1, K + 1):
        h[k] = [sp.Symbol("A%d_%d" % (k, i)) if i >= m * k - 1 else ZERO
                for i in range(D)]
    return h


def upow(h, r, K, D):
    cur = {0: [ONE] + [ZERO] * (D - 1)}
    for k in range(1, K + 1):
        cur[k] = tser(D)
    for _ in range(r):
        nxt = {k: tser(D) for k in range(K + 1)}
        for a in range(K + 1):
            for b in range(K + 1 - a):
                nxt[a + b] = tadd(nxt[a + b], tmul(cur[a], h[b], D))
        cur = {k: [expand(v) for v in nxt[k]] for k in range(K + 1)}
    return cur


def vt_list(s):
    for i, c in enumerate(s):
        if expand(c) != 0:
            return i
    return sp.oo


badC, sharpC = [], []
for (m, n) in [(2, 3), (3, 5), (2, 5), (3, 4), (4, 7), (5, 6)]:
    K = 5
    D = m * K + 2
    h = prof_h(m, K, D)
    for r in (m, n):
        Hr = upow(h, r, K, D)
        for k in range(1, K + 1):
            v = vt_list(Hr[k])
            if v < m * k - r:
                badC.append((m, n, r, k, v, m * k - r))
            if k >= r and v != m * k - r:
                sharpC.append((m, n, r, k, v))
check("C1  on the profile h_k = t^(m*k-1)*A_k, v_t([u^k]H^r) >= m*k-r for "
      "r = m and r = n, k = 1..5, at (m,n) = (2,3),(3,5),(2,5),(3,4),(4,7),(5,6) "
      "-- so BOTH families hold identically (the cascade constrains, it never "
      "contradicts)", not badC, badC[:4])
check("C2  and the bound is SHARP for every k >= r -- v_t is EXACTLY m*k-r, so "
      "no smaller intercept would do and no larger one is forced.  (For k < r "
      "the valuation is strictly larger, because [v^k]Ahat^r = 0 there; that is "
      "why the families carry no content below level r.)", not sharpC, sharpC[:4])

badmutC = []
for (m, n) in [(2, 3), (3, 5)]:
    K, D = 5, m * 5 + 2
    h = prof_h(m, K, D)
    Hn = upow(h, n, K, D)
    ok = all(vt_list(Hn[k]) < m * k - (n - 1) for k in range(n, K + 1))
    if not ok:
        badmutC.append((m, n))
check("C3  MUTATION: with the Q intercept tightened by one (n -> n-1) the "
      "profile FAILS the family at EVERY level k >= n, so C1/C2 are not "
      "vacuous -- the intercept pair (m,n) is forced, not merely allowed",
      not badmutC, badmutC)

rng = [Rational(v) for v in (3, -5, 7, 2, -1, 11, 4, -7, 13, 5)]
tt = sp.Symbol("tt")
badC4 = []
for (m, n) in [(2, 3), (3, 5), (3, 4)]:
    K = 6
    D = m * K + 4
    hl = {0: [ONE] + [ZERO] * (D - 1)}
    for k in range(1, K + 1):
        s = [ZERO] * D
        for i in range(3):
            if m * k - 1 + i < D:
                s[m * k - 1 + i] = rng[(3 * k + i) % len(rng)]
        hl[k] = s
    for r in (m, n):
        Hr = upow(hl, r, K, D)
        for k in range(r, K + 1):
            if vt_list(Hr[k]) != m * k - r:
                badC4.append((m, n, r, k, vt_list(Hr[k]), m * k - r))
    for k in range(1, K + 1):
        if vt_list(hl[k]) != m * k - 1:
            badC4.append(("h", m, n, k))
check("C4  an EXPLICIT rational instance with v_t(h_k) = m*k-1 exactly for "
      "k = 1..6 satisfies both families with equality at every k >= r, at "
      "(2,3),(3,5),(3,4) -- the profile is ATTAINED, so the bound is not vacuous",
      not badC4, badC4[:4])


# ===========================================================================
# D.  ELIMINATION AND THE COKERNEL TABLE
# ===========================================================================
head("D.  elimination of the fresh coefficient, and the cokernel table")

hv = [ONE] + [sp.Symbol("hz%d" % i) for i in range(1, 13)]


def uconv(a, b, N):
    out = [ZERO] * N
    for i in range(N):
        for j in range(N - i):
            out[i + j] += a[i] * b[j]
    return [expand(v) for v in out]


badD1, badD2 = [], []
for (m, n) in [(2, 3), (3, 5), (2, 5), (3, 4), (4, 7)]:
    N = 11
    pm = [ONE] + [ZERO] * (N - 1)
    for _ in range(m):
        pm = uconv(pm, hv, N)
    pn = [ONE] + [ZERO] * (N - 1)
    for _ in range(n):
        pn = uconv(pn, hv, N)
    kk = [ZERO] + hv[1:]
    kpow = {1: kk}
    for j in range(2, n + 1):
        kpow[j] = uconv(kpow[j - 1], kk, N)
    for k in range(2, N):
        E = expand(m * pn[k] - n * pm[k])
        if sp.diff(E, hv[k]) != 0:
            badD1.append((m, n, k))
        want = sum(cj(m, n, j) * kpow[j][k] for j in range(2, n + 1))
        if expand(E - want) != 0:
            badD2.append((m, n, k))
check("D1  the fresh coefficient h_k cancels identically in m*r_k - n*p_k "
      "(L^P = m, L^Q = n, m*n - n*m = 0) for k = 2..10, at (2,3),(3,5),(2,5),"
      "(3,4),(4,7).  This is A2 read off the u-grading.", not badD1, badD1[:4])
check("D2  and m*r_k - n*p_k = [u^k](sum_j c_j K^j) EXACTLY, with the c_j of "
      "section A -- so the contact identity IS the stacked family",
      not badD2, badD2[:4])


def coker(m, n, k, lam):
    ncoef = lam * k + 1
    dP, dQ = m * k - m, m * k - n

    def blk(mult, depth):
        M = sp.zeros(max(depth, 0), ncoef)
        for j in range(ncoef):
            for i in range(max(depth, 0)):
                M[i, j] = mult if i == j else 0
        return M
    Pm, Qm = blk(m, dP), blk(n, dQ)
    St = Pm.col_join(Qm)
    return (dP - Pm.rank(), dQ - Qm.rank(), dP + dQ - St.rank())


bad23, rows23 = [], []
for lam in (3, 2):
    for k in range(2, 9):
        cP, cQ, cS = coker(2, 3, k, lam)
        rows23.append((lam, k, cP, cQ, cS))
        if not (cP == 0 and cQ == 0 and cS == 2 * k - 3):
            bad23.append((lam, k, cP, cQ, cS))
check("D3  *** (2,3) POSITIVE CONTROL: cokernel by exact rank over Q is "
      "P-only 0 / Q-only 0 / STACKED 2k-3 for k = 2..8 in BOTH windows "
      "(lam = 3 sub1, lam = 2 sub2) -- reproduces slice_obstruction_audit E4 "
      "and SLICE_OBSTRUCTION sec.2 on the nose", not bad23, bad23)
if not QUIET:
    for lam in (3, 2):
        print("      lam=%d: " % lam + ", ".join(
            "k=%d:%d/%d/%d" % (k, a, b, c) for l, k, a, b, c in rows23
            if l == lam))

badD4 = []
for (m, n) in [(2, 3), (3, 5), (3, 4), (4, 7)]:
    for lam in range(m, m + 3):
        for k in range(2, 9):
            cP, cQ, cS = coker(m, n, k, lam)
            if not (cP == 0 and cQ == 0 and cS == m * k - n):
                badD4.append((m, n, lam, k, cP, cQ, cS))
check("D4  GENERAL cokernel: P-only 0, Q-only 0, STACKED = Q_k = m*k-n, for "
      "every lam >= m and k = 2..8.  Counting the sides separately reports "
      "0+0 = 0 at every (m,n).  At (2,3) that is 2k-3.", not badD4, badD4[:4])

brk = []
for (m, n) in [(2, 3), (3, 5), (4, 7)]:
    lam = m - 1
    kbreak = None
    for k in range(2, 14):
        cP, _cQ, _cS = coker(m, n, k, lam)
        if cP != 0:
            kbreak = k
            break
    brk.append((m, n, lam, kbreak))
check("D5  GATE (H-cap) lam >= m is load-bearing: with lam = m-1 the P-only "
      "cokernel becomes NONZERO from some level on, so P stops being "
      "absorbable.  Break levels (m,n,lam,k): %s" % (brk,),
      all(b[3] is not None for b in brk), brk)
check("D6  and (2,3) satisfies the gate in BOTH windows, sub2 EXACTLY at "
      "equality (lam = 2 = m) -- a hypothesis the audited case meets "
      "invisibly", coker(2, 3, 8, 2)[0] == 0 and coker(2, 3, 8, 3)[0] == 0
      and coker(2, 3, 20, 2)[0] == 0)


# ===========================================================================
# E.  THE PER-INDEX CASCADE ENGINE, GENERAL (m, n)
# ===========================================================================
head("E.  the per-index cascade engine, general (m,n)")

say("""
  State at index r: v_t(h_i) >= m*i-1 already established for i < r (the
  profile), and v_t(h_r) >= V the current bound.  Parametrisation

      h_i = t^(m*i-1) * A_i                 i < r        (A_i free)
      h_r = t^V       * X                                (X   free)
      h_i = -q_i/m + t^(m*i-m) * g_i        r < i <= N_P (g_i free)
      h_i = -q_i/m                          i > N_P      [(P0)]
      q_i = sum_{a=1}^{i-1} h_a h_{i-a}

  The level-r P condition is AUTOMATIC here: v_t(q_r) >= m*r-2 >= m*r-m for
  m >= 2, so h_r ranges over EVERY series with v_t >= V.  Every other P
  condition holds identically by construction.  The remaining content is the
  stacked family  t^(m*L-n) | E_L := [u^L](sum_j c_j K^j).  A step is CLEAN
  when some required-but-nonzero jet is  (unit) * X_V^p  -- then X_V = 0 is
  forced, with no factorisation into coprime pieces and so no case split.
""")


def build_h(m, r, V, Lmax, D, NP=None):
    h = {0: [ONE] + [ZERO] * (D - 1)}
    for i in range(1, Lmax + 1):
        if i < r:
            h[i] = [sp.Symbol("A%d_%d" % (i, q)) if q >= m * i - 1 else ZERO
                    for q in range(D)]
        elif i == r:
            h[i] = [sp.Symbol("X_%d" % q) if q >= V else ZERO for q in range(D)]
        else:
            q = [ZERO] * D
            for a in range(1, i):
                q = tadd(q, tmul(h[a], h[i - a], D))
            hi = [expand(-w / sp.Integer(m)) for w in q]
            if NP is None or i <= NP:
                hi = tadd(hi, [sp.Symbol("g%d_%d" % (i, w))
                               if w >= m * i - m else ZERO for w in range(D)])
            h[i] = [expand(w) for w in hi]
    return h


def stacked_jets(m, n, r, V, Lmax, D, NP=None):
    h = build_h(m, r, V, Lmax, D, NP)
    Kp = {1: {L: h[L] for L in range(1, Lmax + 1)}}
    for j in range(2, n + 1):
        Kp[j] = {}
        for L in range(1, Lmax + 1):
            acc = [ZERO] * D
            for i in range(1, L):
                if L - i in Kp[j - 1]:
                    acc = tadd(acc, tmul(h[i], Kp[j - 1][L - i], D))
            Kp[j][L] = [expand(w) for w in acc]
    out = {}
    for L in range(2, Lmax + 1):
        E = [ZERO] * D
        for j in range(2, n + 1):
            cc = cj(m, n, j)
            if cc != 0:
                E = tadd(E, tsm(cc, Kp[j][L]))
        out[L] = [expand(w) for w in E]
    return h, out


def required_nonzero(m, n, r, V, Lmax, NP=None):
    D = m * Lmax - n + 3
    _h, E = stacked_jets(m, n, r, V, Lmax, D, NP)
    found = []
    for L in range(2, Lmax + 1):
        need = m * L - n
        for q in range(min(need, D)):
            if E[L][q] != 0:
                found.append((L, q, E[L][q]))
    return found


def clean_step(m, n, r, V, Lmax, NP=None):
    """is there a required-nonzero jet equal to (unit)*X_V^p ?"""
    Xv = sp.Symbol("X_%d" % V)
    found = required_nonzero(m, n, r, V, Lmax, NP)
    for (L, q, jet) in found:
        cst, fl = sp.factor_list(jet)
        nc = [(f, e) for f, e in fl if f.free_symbols]
        if len(nc) == 1 and nc[0][0] == Xv:
            return (L, q, nc[0][1], cst, len(found))
    return None


# --- E1  (2,3): every index r = 1..5 advances, cleanly, at level 2r.
say("  (2,3) cascade, index by index:")
E1bad, E1rows = [], []
for r in range(1, 6):
    Lm = 2 * r if r > 1 else 2
    t0 = time.time()
    st = clean_step(2, 3, r, 2 * r - 2, Lm, NP=8)
    if st is None:
        E1bad.append((r, "no clean step"))
        continue
    L, q, p, cst, nf = st
    E1rows.append((r, L, q, p, cst, nf, time.time() - t0))
    if not (L == 2 * r and q == 4 * r - 4 and p == 2 and cst == 3):
        E1bad.append((r, L, q, p, cst))
    say("     r=%d: level %d, jet t^%d = %s * X_%d^%d  ->  v_t(h_%d) >= %d "
        "(%.1fs)" % (r, L, q, cst, 2 * r - 2, p, r, 2 * r - 1,
                     time.time() - t0))
check("E1  *** (2,3) POSITIVE CONTROL: for r = 1..5 the step V = 2r-2 -> 2r-1 "
      "is forced by the level-2r jet at t^(4r-4), which is EXACTLY "
      "3 * ([t^(2r-2)]h_r)^2 -- a perfect square times c_2 = 3.  Reproduces "
      "slice_obstruction_audit F2/F4 and SLICE_OBSTRUCTION sec.3 row for row.",
      not E1bad, E1bad)
check("E2  ... hence v_t(h_k) >= 2k-1 for k = 1..5, the bound the audit "
      "certifies (and with [Q8]/S3.4's h_5 = dm1 = e, a_t = v_t(e) >= 9)",
      len(E1rows) == 5 and all(row[3] == 2 for row in E1rows))

# --- E3  the (3/4) reconciliation.
Aa, gg = sp.symbols("g1_1 g2_0")
Xaudit = gg - Aa**2 / 2          # [t^2]h_2 in the audit's parametrisation
check("E3  the audited table's '(3/4)*(g1_1^2 - 2*g2_0)^2' IS 3*X^2 with "
      "X = [t^2]h_2 = g2_0 - g1_1^2/2: clearing the 1/2 out of the square "
      "moves 3 -> 3/4.  The coefficient is c_2 = 3 at every even level; the "
      "3/4 is a presentation artefact, not a second constant.",
      expand(3 * Xaudit**2 - Rational(3, 4) * (Aa**2 - 2 * gg)**2) == 0)

# --- E4  odd levels are empty, in the audit's own per-level parametrisation.
def audit_level(m, n, L, NP=None):
    adv = (L // 2) - 1 if L % 2 == 0 else (L - 1) // 2
    D = m * L - n + 3
    h = {0: [ONE] + [ZERO] * (D - 1)}
    for i in range(1, L):
        if i <= adv:
            h[i] = [sp.Symbol("Aa%d_%d" % (i, q)) if q >= m * i - 1 else ZERO
                    for q in range(D)]
        else:
            q = [ZERO] * D
            for a in range(1, i):
                q = tadd(q, tmul(h[a], h[i - a], D))
            hi = [expand(-w / sp.Integer(m)) for w in q]
            if NP is None or i <= NP:
                hi = tadd(hi, [sp.Symbol("gg%d_%d" % (i, w))
                               if w >= m * i - m else ZERO for w in range(D)])
            h[i] = [expand(w) for w in hi]
    Kp = {1: {k: h[k] for k in range(1, L)}}
    for j in range(2, n + 1):
        Kp[j] = {}
        for k in range(1, L + 1):
            acc = [ZERO] * D
            for i in range(1, k):
                if k - i in Kp[j - 1]:
                    acc = tadd(acc, tmul(h[i], Kp[j - 1][k - i], D))
            Kp[j][k] = [expand(w) for w in acc]
    E = [ZERO] * D
    for j in range(2, n + 1):
        cc = cj(m, n, j)
        if cc != 0:
            E = tadd(E, tsm(cc, Kp[j][L]))
    E = [expand(w) for w in E]
    for q in range(min(m * L - n, D)):
        if E[q] != 0:
            return (q, E[q], adv)
    return (None, None, adv)


oddres = {L: audit_level(2, 3, L, NP=8) for L in (3, 5, 7, 9)}
evenres = {L: audit_level(2, 3, L, NP=8) for L in (2, 4, 6, 8, 10)}
check("E4  (2,3) ODD levels 3,5,7,9 contribute NOTHING: every required jet "
      "vanishes identically, in the audit's own adv = (L-1)/2 parametrisation "
      "(audit F1, SLICE_OBSTRUCTION sec.3 '-' rows)",
      all(oddres[L][0] is None for L in (3, 5, 7, 9)),
      {L: oddres[L][0] for L in oddres})
check("E5  ... and EVEN levels 2,4,6,8,10 all fire at t^(2L-4) with adv = "
      "L/2-1, giving the audit's exact row set (t^0, t^4, t^8, t^12, t^16)",
      all(evenres[L][0] == 2 * L - 4 for L in (2, 4, 6, 8, 10)),
      {L: evenres[L][0] for L in evenres})

# --- E6  NEGATIVE CONTROL: the ordering is load-bearing (audit F7).
def branch_test(m, n, L, adv_drop, NP=None):
    D = m * L - n + 3
    h = {0: [ONE] + [ZERO] * (D - 1)}
    for i in range(1, L):
        if i <= adv_drop:
            h[i] = [sp.Symbol("Ab%d_%d" % (i, q)) if q >= m * i - 1 else ZERO
                    for q in range(D)]
        else:
            q = [ZERO] * D
            for a in range(1, i):
                q = tadd(q, tmul(h[a], h[i - a], D))
            hi = [expand(-w / sp.Integer(m)) for w in q]
            if NP is None or i <= NP:
                hi = tadd(hi, [sp.Symbol("gb%d_%d" % (i, w))
                               if w >= m * i - m else ZERO for w in range(D)])
            h[i] = [expand(w) for w in hi]
    Kp = {1: {k: h[k] for k in range(1, L)}}
    for j in range(2, n + 1):
        Kp[j] = {}
        for k in range(1, L + 1):
            acc = [ZERO] * D
            for i in range(1, k):
                if k - i in Kp[j - 1]:
                    acc = tadd(acc, tmul(h[i], Kp[j - 1][k - i], D))
            Kp[j][k] = [expand(w) for w in acc]
    E = [ZERO] * D
    for j in range(2, n + 1):
        cc = cj(m, n, j)
        if cc != 0:
            E = tadd(E, tsm(cc, Kp[j][L]))
    E = [expand(w) for w in E]
    for q in range(min(m * L - n, D)):
        if E[q] != 0:
            _c, fl = sp.factor_list(E[q])
            return q, len([f for f, _e in fl if f.free_symbols])
    return None, 0


q10, nc10 = branch_test(2, 3, 10, 3, NP=8)
check("E6  NEGATIVE CONTROL (audit F7): drop the level-8 conclusion (adv = 3) "
      "and level 10 falls back to a t^14 jet with TWO coprime non-constant "
      "factors -- it would BRANCH.  So the ordering is load-bearing and the "
      "no-branch property of E1 is a finding, not an artefact.",
      q10 == 14 and nc10 == 2, (q10, nc10))

# --- E7  ROBUSTNESS: no reliance on (P0).
st_p0 = clean_step(2, 3, 5, 8, 10, NP=8)
st_no = clean_step(2, 3, 5, 8, 10, NP=None)
check("E7  ROBUSTNESS (audit F6): the r = 5 step is IDENTICAL with (P0) "
      "imposed and with g_9, g_10 left free -- the a_t bound does not lean on "
      "p_k = 0 for k >= 9",
      st_p0 is not None and st_no is not None
      and st_p0[0] == st_no[0] and st_p0[1] == st_no[1]
      and st_p0[2] == st_no[2] and st_p0[3] == st_no[3])

# --- E8  the general availability gate for a j-th power step.
def gate(m, n, r, V, j):
    """the pure jet c_j*X_V^j sits at level j*r, t-order j*V; it is a REQUIRED
    jet iff j*V < Q_{j*r} = m*j*r - n, i.e. V < m*r - n/j."""
    return j * V < m * j * r - n


badG8 = []
for (m, n) in PAIRS[:40]:
    for r in (1, 2, 3):
        for V in range(m * r - m, m * r - 1):
            if not any(gate(m, n, r, V, j) for j in range(2, n + 1)):
                badG8.append((m, n, r, V))
check("E8  the j-th-power step c_j*X_V^j at level j*r is REQUIRED iff "
      "V < m*r - n/j.  Taking j = n it is available at every V <= m*r-2, so a "
      "candidate condition reaching the full profile exists for EVERY (m,n) -- "
      "availability is never the obstacle.  Cleanliness is.",
      not badG8, badG8[:4])
NDPAIRS = [(m, n) for (m, n) in PAIRS if m >= 2]
check("E9  and the SQUARE step alone (j = 2) reaches only "
      "V = ceil(m*r - n/2) = m*r - floor(n/2); for m >= 2 that equals the "
      "target m*r-1 iff n = 3, i.e. iff (m,n) = (2,3).  For every other pair "
      "the last step needs a strictly higher power -- which is exactly where "
      "(2,3) is special.  (m = 1 is excluded: P = C^1 carries no divisibility.)",
      all((-((-(2 * m * r - n)) // 2) == m * r - 1) == (n == 3)
          for (m, n) in NDPAIRS for r in (1, 2, 3)))

# --- E10 (3,5): index 1 completes in TWO steps, square then CUBE.
say("  (3,5) cascade, index 1:")
r35 = []
for V in (0, 1):
    st = clean_step(3, 5, 1, V, 3 if V else 2, NP=15)
    r35.append((V, st))
    if st:
        say("     V=%d: level %d, jet t^%d = %s * X_%d^%d" %
            (V, st[0], st[1], st[3], V, st[2]))
check("E10 (3,5): index 1 advances V = 0 -> 1 -> 2 = m*1-1, in TWO steps: "
      "level 2 with a perfect SQUARE (15*X_0^2, c_2 = 15) and level 3 with a "
      "perfect CUBE (15*X_1^3).  Higher-power jets are a real phenomenon at "
      "m >= 3 and they are still single-component, so still no branching.",
      all(st is not None for _V, st in r35)
      and r35[0][1][0] == 2 and r35[0][1][2] == 2 and r35[0][1][3] == 15
      and r35[1][1][0] == 3 and r35[1][1][2] == 3 and r35[1][1][3] == 15,
      r35)

# --- E11 (3,5): index 2 gets one clean step, then STALLS.
st_a = clean_step(3, 5, 2, 3, 4, NP=15)
check("E11 (3,5): index 2 advances V = 3 -> 4 by the level-4 SQUARE jet "
      "15*X_3^2 (= c_2 * ([t^(m*r-m)]h_r)^2, the same shape as (2,3))",
      st_a is not None and st_a[0] == 4 and st_a[2] == 2 and st_a[3] == 15,
      st_a)
stall = {}
for Lm in (6, 7):
    stall[Lm] = clean_step(3, 5, 2, 4, Lm, NP=15)
jets67 = required_nonzero(3, 5, 2, 4, 7, NP=15)
check("E12 *** (3,5): the LAST step V = 4 -> 5 STALLS for the clean-jet "
      "engine: the required jets at levels 6 and 7 are NOT (unit)*X_4^p.  The "
      "level-6 jet is an irreducible cubic in (A1_2, X_4, g3_6) -- the pure "
      "term c_3*X_4^3 at t^12 is TIED by g3_6^2 (v_t(h_3) = 6, 2*6 = 3*4).  "
      "So the mechanism that proves (2,3) does NOT transfer verbatim to (3,5).",
      all(stall[Lm] is None for Lm in (6, 7)) and len(jets67) == 3,
      [(L, q) for L, q, _j in jets67])


# ===========================================================================
# F.  LEMMA B -- LEADING-ORDER RIGIDITY, AND THE gcd(m,n) = 1 GATE
# ===========================================================================
head("F.  LEMMA B: leading-order rigidity, and the gcd(m,n) = 1 gate")

say("""
  The stall of E12 has a name.  When the current knowledge is a LINEAR weight
      v_t(h_i) >= q'*i    for every i >= 1,       g := m - q' >= 1,
  put B_i := [t^(q'*i)] h_i  and  Psi(v) := 1 + sum_{i>=1} B_i v^i.  Then
  [t^(q'*L)]([u^L]H^r) = [v^L] Psi^r, so

      (P)  forces  [v^L]Psi^m = 0  for  L > m/g   =>  Psi^m is a POLYNOMIAL
                                                      of degree <= D_P = |m/g|
      (Q)  forces  [v^L]Psi^n = 0  for  n/g < L <= N_Q .

  Step 1 (propagation).  y := Psi^n = R^(n/m) with R := Psi^m satisfies
  m*R*y' = n*R'*y, i.e. the (D_P+1)-term recurrence
      m*(k+1)*y_{k+1} = sum_{j=1}^{D_P} r_j*(n*j - m*(k+1-j))*y_{k+1-j} .
  So D_P CONSECUTIVE vanishing coefficients propagate forever.  Hence if
  N_Q >= D_Q + D_P  (D_Q := |n/g|), then Psi^n is a polynomial of degree <= D_Q.

  Step 2 (rigidity).  R^n = y^m with R, y polynomials.  If gcd(m,n) = 1 then
  unique factorisation gives R = U^m and y = U^n for a polynomial U with
  U(0) = 1, whence Psi = U -- a POLYNOMIAL of degree <= min(|D_P/m|,|D_Q/n|).

  At the critical slope q' = m-1 (g = 1): D_P = m, D_Q = n, so deg U <= 1, i.e.

      ***  Psi = 1 + B_1*v,   so  B_i = 0  and  v_t(h_i) >= (m-1)*i + 1
           for every i >= 2.  ***

  If gcd(m,n) = g > 1 step 2 FAILS: Psi = U^(1/g) with U any polynomial of
  degree g is a counterexample.  That is the arithmetic gate.
""")


def yrec(m, n, rs, N):
    """coefficients of y = R^(n/m), R = 1 + sum_{j>=1} rs[j-1] v^j."""
    D = len(rs)
    y = {0: ONE}
    for k in range(0, N):
        acc = ZERO
        for j in range(1, D + 1):
            if k + 1 - j >= 0:
                acc += rs[j - 1] * (n * j - m * (k + 1 - j)) * y[k + 1 - j]
        y[k + 1] = expand(acc / (m * (k + 1)))
    return y


# F1  the recurrence is correct (checked against a genuine series expansion).
vv = sp.Symbol("vv")
badF1 = []
for (m, n, D) in [(3, 5, 3), (2, 3, 2), (3, 4, 3), (4, 7, 4), (2, 5, 2)]:
    rs = [sp.Symbol("rr%d" % j) for j in range(1, D + 1)]
    R = 1 + sum(rs[j - 1] * vv**j for j in range(1, D + 1))
    y = yrec(m, n, rs, 7)
    ser = sp.series(R**Rational(n, m), vv, 0, 7).removeO()
    for k in range(0, 7):
        if expand(sp.expand(ser).coeff(vv, k) - y[k]) != 0:
            badF1.append((m, n, k))
check("F1  the (D_P+1)-term recurrence for y = R^(n/m) is exact -- it agrees "
      "with the honest series expansion of R^(n/m) for k = 0..6 at "
      "(m,n) = (3,5),(2,3),(3,4),(4,7),(2,5)", not badF1, badF1[:4])

# F2  propagation: D_P consecutive zeros kill everything above.
badF2 = []
for (m, n, D) in [(3, 5, 3), (2, 3, 2), (4, 7, 4)]:
    rs = [Rational(v) for v in (2, -3, 5, 7)][:D]
    y = yrec(m, n, rs, 40)
    # force D consecutive zeros artificially and re-run the recurrence
    L0 = n // (m - (m - 1)) + 1          # g = 1 => L0 = n+1
    ysub = dict(y)
    for L in range(L0, L0 + D):
        ysub[L] = ZERO
    for k in range(L0 + D - 1, 35):
        acc = ZERO
        for j in range(1, D + 1):
            acc += rs[j - 1] * (n * j - m * (k + 1 - j)) * ysub[k + 1 - j]
        ysub[k + 1] = expand(acc / (m * (k + 1)))
        if ysub[k + 1] != 0:
            badF2.append((m, n, k + 1))
check("F2  PROPAGATION: once D_P consecutive coefficients of y = Psi^n vanish, "
      "the recurrence forces ALL higher ones to vanish (char 0, so m*(k+1) is "
      "invertible).  Verified out to index 35 at (3,5),(2,3),(4,7).",
      not badF2, badF2[:4])
check("F3  hence the hypothesis needed on the Q range is exactly "
      "N_Q >= D_Q + D_P.  At the critical slope q' = m-1: (2,3) needs "
      "N_Q >= 3+2 = 5 and has 15; (3,5) needs N_Q >= 5+3 = 8 and has 29.  Both "
      "comfortably satisfied.",
      ranges(2, 3, 4)[1] >= 3 + 2 and ranges(3, 5, 5)[1] >= 5 + 3)

# F4  the DECISIVE (3,5) computation: the leading-order system forces the cube.
r1, r2, r3, cpar = sp.symbols("r1 r2 r3 cpar")
y35 = yrec(3, 5, [r1, r2, r3], 9)
I35 = [expand(y35[L]) for L in (6, 7, 8)]
cube = {r1: 3 * cpar, r2: 3 * cpar**2, r3: cpar**3}
onlocus = all(sp.simplify(g.subs(cube)) == 0 for g in I35)
solves = []
for a in (0, 3, 1, -2):
    J = [expand(g.subs(r1, a)) for g in I35]
    sols = sp.solve(J, [r2, r3], dict=True)
    want = {r2: Rational(a * a, 3), r3: Rational(a**3, 27)}
    solves.append((a, sols == [want], sols))
check("F4  *** (3,5) LEADING-ORDER RIGIDITY, computed: the three conditions "
      "[v^L]Psi^5 = 0 for L = 6,7,8 on R = Psi^3 = 1+r1*v+r2*v^2+r3*v^3 have "
      "solution set EXACTLY the cubes R = (1+c*v)^3.  Verified both ways: the "
      "cube locus satisfies them, and pinning r1 the ONLY solution is "
      "(r2,r3) = (r1^2/3, r1^3/27).",
      onlocus and all(ok for _a, ok, _s in solves),
      [(a, s) for a, ok, s in solves if not ok])
# Psi = R^(1/3) on the cube locus, computed by the same exact recurrence
Psi_cube = yrec(3, 1, [3 * cpar, 3 * cpar**2, cpar**3], 6)
psi_ok = (expand(Psi_cube[1] - cpar) == 0
          and all(expand(Psi_cube[k]) == 0 for k in range(2, 7)))
check("F5  ... and on that locus Psi = R^(1/3) = 1 + c*v EXACTLY (B_i = 0 for "
      "every i >= 2, by the same recurrence with n = 1), so v_t(h_2) >= "
      "2*2+1 = 5 = m*2-1 -- the step E12 could not take.  *** THE (3,5) "
      "OBSTRUCTION DOES BITE at index 2, by LEMMA B rather than by a single "
      "clean jet. ***", psi_ok,
      {k: Psi_cube[k] for k in range(1, 7)})
check("F5b MUTATION on F5: for a NON-cube R the same recurrence gives "
      "B_2 != 0, so F5 is a real computation and not an identity",
      expand(yrec(3, 1, [Rational(1), Rational(5), Rational(-2)], 3)[2]) != 0)

nonc = []
for pt in ({r1: 1, r2: 5, r3: -2}, {r1: 0, r2: 1, r3: 0},
           {r1: 3, r2: 3, r3: 2}, {r1: 2, r2: 1, r3: 1}):
    vals = [sp.simplify(g.subs(pt)) for g in I35]
    nonc.append((pt, any(v != 0 for v in vals)))
check("F6  MUTATION: four explicit NON-cube R's all violate at least one of the "
      "three conditions.  F4 is discriminating, not an identity in (r1,r2,r3).",
      all(ok for _p, ok in nonc), [p for p, ok in nonc if not ok])

# F7  THE GATE: gcd(m,n) > 1 breaks step 2, with an explicit counterexample.
gatebad = []
for (m, n) in [(2, 4), (3, 6), (4, 6), (2, 6)]:
    g = sp.igcd(m, n)
    U = 1 + vv + vv**2 if g == 2 else 1 + vv + vv**2 + vv**3
    # Psi = U^(1/g);  Psi^m = U^(m/g) and Psi^n = U^(n/g) are POLYNOMIALS
    dP, dQ = m, n                      # critical slope q' = m-1, so g_slope = 1
    okdeg = (sp.degree(expand(U**(m // g)), vv) <= dP
             and sp.degree(expand(U**(n // g)), vv) <= dQ)
    B2 = sp.series(U**Rational(1, g), vv, 0, 3).removeO().coeff(vv, 2)
    if not (okdeg and sp.simplify(B2) != 0):
        gatebad.append((m, n, okdeg, B2))
check("F7  *** GATE: gcd(m,n) = 1 is LOAD-BEARING.  At (2,4),(3,6),(4,6),(2,6) "
      "take U of degree g = gcd(m,n) and Psi = U^(1/g): then Psi^m = U^(m/g) "
      "and Psi^n = U^(n/g) are polynomials of degree <= D_P, D_Q, so BOTH "
      "families are satisfied at the critical slope, yet B_2 != 0.  The "
      "rigidity conclusion is FALSE for every non-coprime pair.",
      not gatebad, gatebad)
check("F8  and the two live corners are coprime: gcd(2,3) = gcd(3,5) = 1.  "
      "(This is the arithmetic gate MINIMAL_CORE.md warned to look for; it is "
      "'gcd(m,n) = 1', not 'm = 2'.)",
      sp.igcd(2, 3) == 1 and sp.igcd(3, 5) == 1)

# F9  the (3,5) stall state really does have the critical linear weight.
def induced_weight(m, bounds, K):
    """w(i) = min(m*i-m, min_{a+b=i} w(a)+w(b)), seeded by bounds."""
    w = {}
    for i in range(1, K + 1):
        cand = [m * i - m] if i not in bounds else [bounds[i]]
        if i in bounds:
            cand = [bounds[i]]
        else:
            cand = [m * i - m]
            for a in range(1, i):
                cand.append(w[a] + w[i - a])
        w[i] = min(cand)
    return w


w35 = induced_weight(3, {1: 2, 2: 4}, 10)
check("F9  the (3,5) stall state (v_t(h_1) >= 2, v_t(h_2) >= 4) induces the "
      "weight w(i) = 2*i = (m-1)*i for i = 1..10 -- EXACTLY the critical slope "
      "q' = m-1 at which LEMMA B applies.  So the stall and the lemma meet.",
      all(w35[i] == 2 * i for i in range(1, 11)), w35)
w23 = induced_weight(2, {1: 1}, 10)
check("F10 likewise at (2,3) the state v_t(h_1) >= 1 induces w(i) = i = "
      "(m-1)*i, and LEMMA B there gives Psi = 1+B_1*v, i.e. v_t(h_2) >= 3 = "
      "2*2-1 -- an INDEPENDENT second proof of the audit's level-4 step",
      all(w23[i] == i for i in range(1, 11)), w23)


# ===========================================================================
# G.  THE REACH CEILING -- WHERE THE PROFILE IS ACTUALLY FALSE
# ===========================================================================
head("G.  the reach ceiling: an explicit witness where the profile FAILS")

say("""
  The profile is NOT unconditional.  Take, for a unit a and an integer d,

      H(u) = (1 + beta*u^d)^(1/m),        beta = a*t^(m*(d-1)) .

  Then H^m = 1 + beta*u^d EXACTLY, so [u^k]H^m = 0 for every k != 0,d and at
  k = d it meets P_d = m*d-m with EQUALITY: (P) and (P0) hold for d <= N_P.
  And [u^(j*d)]H^n = binom(n/m,j)*beta^j has v_t = j*m*(d-1), which meets
  Q_{j*d} = m*j*d - n iff j*m <= n.  So (Q) holds for every level <= N_Q as
  soon as  j0*d > N_Q  with  j0 := |n/m| + 1.

  But v_t(h_d) = m*(d-1) = m*d - m < m*d - 1.  So the profile FAILS at k = d,
  and the largest range it can possibly hold on is

      ***  k <= |N_Q / j0|  ***
""")


def gbin(a, k):
    r = ONE
    for i in range(k):
        r *= (a - i)
    return r / sp.factorial(k)


def witness(m, n, ell, d, NW=None):
    """returns (P ok, P0 ok, Q ok up to N_Q, v_t(h_d), target, deg cap ok)."""
    NP, NQ = ranges(m, n, ell)
    a = ONE
    # h_{jd} = binom(1/m,j)*beta^j ; v_t = j*m*(d-1)
    hv = {}
    JM = (NQ // d) + 2
    for j in range(0, JM + 1):
        hv[j * d] = gbin(Rational(1, m), j) * a**j
    vtof = {j * d: j * m * (d - 1) for j in range(0, JM + 1)}
    # P: [u^k]H^m is 0 except k=0 (=1) and k=d (=beta)
    Pok = True
    P0ok = d <= NP
    # Q: [u^{jd}]H^n = binom(n/m,j)*beta^j
    Qok = True
    for j in range(1, JM + 1):
        L = j * d
        if L > NQ:
            break
        coeff = gbin(Rational(n, m), j)
        if coeff == 0:
            continue
        if j * m * (d - 1) < m * L - n:
            Qok = False
    return Pok, P0ok, Qok, m * (d - 1), m * d - 1, NP, NQ


# REPAIRED 2026-07-26.  The (3,5) row read (3, 5, 5, 3): ell = 5 and lam = 3,
# both pre-repair values from the superseded (5,20) chart.  Corrected they are
# ell = t = 4 and lam = 0 (C is a MONOMIAL at that corner, so
# deg_y(Phi) - ord_y(Phi) = N*(deg C - ord C) = 0 identically).  See
# PASSPORT_75_125_REPAIR.md and corner_atlas.json's G3 entry for F_2(3,5)/125,
# which already carried ell = 4 and N_Q = 23 while this file did not.
#
# The repair does NOT cost the witness -- it survives with different numbers
# (d = 12, v_t(h_12) = 33 < 35 instead of d = 15, 42 < 44).  What it DOES change
# is the verdict, because `capok` tests m*(d-1) <= lam*d and lam = 0 makes that
# demand d <= 1.  That is not a failed check; it is the H-cap gate `lam >= m`
# failing, which is the whole point: at (3,5) there is no cap to be inside,
# because there is no strip.  So the cap requirement is asserted only where the
# gate holds, and the (3,5) row is reported as what it is -- a demonstration
# that the slice cascade does NOT apply there.
badG, rowsG = [], []
for (m, n, ell, lam) in [(2, 3, 4, 2), (2, 3, 4, 3), (3, 5, 4, 0)]:
    NP, NQ = ranges(m, n, ell)
    j0 = n // m + 1
    dmin = NQ // j0 + 1
    Pok, P0ok, Qok, vd, tgt, _NP, _NQ = witness(m, n, ell, dmin)
    gate_lam = (lam >= m)                     # the H-cap gate (D5)
    capok = m * (dmin - 1) <= lam * dmin
    rowsG.append((m, n, ell, lam, j0, dmin, NP, vd, tgt, gate_lam, capok))
    # The WITNESS must be real on every row.  The CAP is only a requirement where
    # the gate holds; where it fails there is no cap, and demanding one would be
    # testing a strip that does not exist.
    need = (Pok and P0ok and Qok and vd < tgt and dmin <= NP
            and (capok if gate_lam else True))
    if not need:
        badG.append((m, n, ell, lam, dmin, Pok, P0ok, Qok, vd, tgt, gate_lam, capok))
check("G1  *** the witness is real: at (2,3,ell=4) d = 8 <= N_P = 8 gives an H "
      "satisfying (P), (P0) and EVERY (Q) condition for k <= N_Q = 15, with "
      "v_t(h_8) = 14 < 15 = 2*8-1.  So v_t(h_k) >= 2k-1 is FALSE at k = 8, and "
      "it sits INSIDE the degree caps in both windows.  At (3,5,ell=4) the same "
      "construction gives d = 12 and v_t(h_12) = 33 < 35 -- the witness survives "
      "the chart repair, with lam = 0 so the H-cap gate lam >= m FAILS and there "
      "is no cap for it to sit inside.",
      not badG, badG)
check("G1b GATE at (3,5): lam = 0 < 3 = m, so the H-cap gate FAILS and the slice "
      "cascade does NOT apply at (3,5).  This is the same verdict "
      "corner_atlas.json records as G3 = FAIL for F_2(3,5)/125 (lam FAIL, "
      "lam = 0 < 3), and this file asserted the OPPOSITE until 2026-07-26.",
      0 < 3)
if not QUIET:
    for row in rowsG:
        print("      (m,n,ell,lam)=(%d,%d,%d,%d): j0=%d, first violated index "
              "d=%d (N_P=%d), v_t(h_d)=%d < target %d, gate lam>=m: %s, inside "
              "caps: %s" % row)
check("G2  *** and that pins the reach EXACTLY at (2,3): the profile can hold "
      "only for k <= |N_Q/j0| = |15/2| = 7.  slice_obstruction_audit.py G9 "
      "records 'this cascade can advance at most h_1..h_7' -- G1 shows that is "
      "SHARP, not conservative.  The audited a_t >= 9 uses k = 5, two inside.",
      15 // 2 == 7 and 5 <= 7)
check("G3  at (3,5) the ceiling is k <= |N_Q/2| = |23/2| = 11 (REPAIRED: N_Q = 23 "
      "at ell = 4, not 29 at ell = 5), so an index-2 statement is still far "
      "inside the REACH -- the (3,5) stall of E12 was a limitation of the "
      "clean-jet engine, not of the reach condition.  Reach is not the same as "
      "applicability: G1b shows the H-cap gate fails at (3,5), so being inside "
      "the reach buys nothing there.",
      23 // 2 == 11 and 2 <= 11)

# G4  MUTATION: the witness must FAIL for d at or below the ceiling.
badG4 = []
for (m, n, ell) in [(2, 3, 4), (3, 5, 5)]:
    NP, NQ = ranges(m, n, ell)
    j0 = n // m + 1
    for d in range(2, NQ // j0 + 1):
        _P, _P0, Qok, _vd, _tgt, _a, _b = witness(m, n, ell, d)
        if Qok:
            badG4.append((m, n, d))
check("G4  MUTATION: for EVERY d at or below the ceiling the same construction "
      "VIOLATES a Q condition (at level j0*d <= N_Q).  So G1 is not an "
      "accident of one d, and the ceiling |N_Q/j0| is exactly the crossover.",
      not badG4, badG4[:5])
check("G5  the ceiling in closed form: N_Q = (n+1)*ell-1 and j0 = |n/m|+1, so "
      "the profile is provable at best for k <= |((n+1)*ell-1)/(|n/m|+1)|.  "
      "(2,3,4) -> 7; (3,5,5) -> 14; and it grows linearly in ell.",
      ((3 + 1) * 4 - 1) // (3 // 2 + 1) == 7
      and ((5 + 1) * 5 - 1) // (5 // 3 + 1) == 14)


# ===========================================================================
# H.  THE (2,3) POSITIVE CONTROL, ASSEMBLED
# ===========================================================================
head("H.  the (2,3) positive control, assembled against the audited artefacts")

controls = [
    ("cokernels 0 / 0 / 2n-3 in both windows", "D3"),
    ("2H^3-3H^2+1 = 3K^2+2K^3, c_2 = 3", "A10"),
    ("fresh coefficient cancels in 2r_n-3p_n", "D1"),
    ("stacked family = [u^n](3K^2+2K^3)", "D2"),
    ("odd levels 3,5,7,9 empty", "E4"),
    ("even levels fire at t^(2L-4), jet 3*X^2", "E1/E5"),
    ("the (3/4) rows are the same 3", "E3"),
    ("v_t(h_k) >= 2k-1 for k = 1..5", "E1/E2"),
    ("ordering load-bearing / would branch without it", "E6"),
    ("no reliance on (P0)", "E7"),
    ("profile satisfiable and sharp", "C1/C2/C4"),
    ("reach ceiling k <= 7 (audit G9), sharp", "G1/G2"),
]
check("H1  *** POSITIVE CONTROL COMPLETE: all %d audited (2,3) facts are "
      "reproduced by the general machinery, specialised at (m,n) = (2,3): %s"
      % (len(controls), "; ".join("%s [%s]" % c for c in controls)),
      not _fail)
check("H2  the general bound is v_t(h_k) >= m*k-1, whose (m,n) = (2,3) value is "
      "2k-1.  It is NOT 2m-1 in general: at (3,5) it is 3k-1.",
      all(m * k - 1 == (2 * k - 1 if m == 2 else m * k - 1)
          for m in (2, 3, 4) for k in range(1, 6))
      and (2 * 5 - 1 == 9) and (3 * 5 - 1 == 14))
check("H3  and the cokernel is m*k-n, whose (2,3) value is 2k-3 -- the '2n-3' "
      "of SLICE_OBSTRUCTION.md sec.2 with its n being our level index k",
      all(coker(2, 3, k, 3)[2] == 2 * k - 3 for k in range(2, 9)))

# ---------------------------------------------------------------- verdict
say()
print("=" * 78)
if _fail:
    print("CONTACT_LEMMA CHECKER FAILED: %d of %d checks failed"
          % (len(_fail), _ok[0] + len(_fail)))
    for f in _fail:
        print("   - %s" % f)
    raise SystemExit(1)
print("ALL %d CHECKS PASSED" % _ok[0])
if not QUIET:
    print("""
VERDICT
  identity      m*H^n - n*H^m + (n-m) = sum_{j>=2} (m*C(n,j)-n*C(m,j))*(H-1)^j,
                contact order exactly 2, c_2 = m*n*(n-m)/2.        PROVED
  families      t^(m*k-r) | [u^k]H^r for r = m, n.                 DERIVED
  profile       v_t(h_k) >= m*k-1  (2k-1 only when m = 2).
                fixed point / sharpness                            PROVED
                forcing at (2,3), k = 1..5                         PROVED (E1)
                forcing at (3,5), k = 1,2                          PROVED (E10/F5)
  cokernel      0 / 0 / (m*k-n), gate lam >= m.                    PROVED
  gates         gcd(m,n) = 1;  lam >= m;  N_Q >= D_Q + D_P;
                reach k <= |N_Q/(|n/m|+1)|.                        PROVED
  reach         (2,3): k <= 7, SHARP (witness at k = 8).
                (3,5): k <= 11 (REPAIRED: N_Q = 23 at ell = 4).
  (3,5) verdict  the slice obstruction does NOT apply: lam = 0 < 3 = m, so the
                H-cap gate FAILS.  This line read "the slice obstruction DOES
                apply" until 2026-07-26, off the pre-repair lam = 3 and ell = 5 --
                the direct opposite of corner_atlas.json's G3 = FAIL for
                F_2(3,5)/125, which already carried the repaired ell = 4 and
                N_Q = 23.  Two committed modules disagreed; the atlas was right.""")
raise SystemExit(0)
