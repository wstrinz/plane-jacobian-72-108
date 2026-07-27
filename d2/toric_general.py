#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
toric_general.py -- is the toric syzygy `6*W*Z = e^5` a FAMILY-LEVEL fact?

Companion document: TORIC_GENERAL.md.
Prerequisite: TORIC_SYZYGY.md / toric_syzygy.py (the (72,108) identity itself),
CONTACT_LEMMA.md (the `lam >= m` gate), PASSPORT_75_125_REPAIR.md (the repaired
(75,125) inputs).

THREE QUESTIONS, THREE ANSWERS (all reproduced below):

  Q1  is the exponent `m + n`?                       ->  NO.  It is forced by the
      chart exponent `t` alone, through a weight-divisibility condition that has
      a solution only at `t = 4`.  At (72,108) FIVE different formulas evaluate
      to 5 (`m+n`, `t+1`, `(4t+9)/(t+1)`, `2n-1`, `mn-1`); they separate at
      (3,5) and at `t = 5,6`, and the survivor is the `t` one.
  Q2  does (3,5) admit the identity?                 ->  NO, and COMPLETELY, not
      merely below a search bound: at (75,125) the weight-admissible exponents
      for a product of two window-Hankel minors are exactly k in {5,6,7,8}
      (m+n = 8 IS among them, so the prediction got a fair test), and exact
      rational linear algebra kills all four.
  Q3  what is `lam` at (75,125)?                     ->  `lam = 0`.  The (H-cap)
      gate `lam >= m = 3` FAILS.  (72,108) has lam = 2 = m, exactly at equality.

Everything is EXACT (sympy over Q, or Fraction linear algebra over Q).  Nothing
is retyped: the generators come from `g_system_75_125.build_gsystem` (the
committed parametric builder, whose (72,108) output is checked against
`g_system_75_125.published_72108`), and the `lam` arithmetic comes from
`window_functions_75_125.window_law` / `.family`.

Run:
    python toric_general.py            # verbose
    python toric_general.py --quiet    # one line per group; exit 0 iff all pass
    python toric_general.py --fast     # skip the two slowest groups (D, and the
                                       # (2,3,6) rung of C)

DISCIPLINE.  Every positive check is paired with a MUTATION CONTROL that must
fail, and every negative check is paired with a POSITIVE CONTROL run through the
same code path (so a "no identity found" cannot be a broken search).
"""

from __future__ import annotations

import argparse
import math
import sys
import time
from fractions import Fraction

import sympy as sp

import g_system_75_125 as gsys
import window_functions_75_125 as wf

PHI = sp.Symbol("Phi")


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
            print("  [%s] %-40s %s" % ("PASS" if ok else "FAIL", name, detail))
        return bool(ok)

    def mut(self, name: str, all_mutants_fail: bool, detail: str = "") -> bool:
        return self.ck("MUT " + name, all_mutants_fail, detail)

    def note(self, text: str) -> None:
        if not self.quiet:
            print("       . " + text)

    def report(self) -> int:
        bad = [r for r in self.rows if not r[1]]
        print("\n%s  toric_general: %d/%d checks pass"
              % ("FAIL" if bad else "OK  ", len(self.rows) - len(bad), len(self.rows)))
        for n, _, d in bad:
            print("   FAILED: %s   %s" % (n, d))
        return 1 if bad else 0


# =============================================================================
#  The parametric G-system, its weights, and exact graded linear algebra
# =============================================================================
_SYS: dict = {}


def system(a: int, b: int, t: int) -> dict:
    """Build (or reuse) the (a,b,t) G-system.  a = m (linear window S^a),
    b = n (forcing window S^b), t = chart exponent."""
    key = (a, b, t)
    if key in _SYS:
        return _SYS[key]
    r = gsys.build_gsystem(a, b, t, 1, 1)
    D = (a - 1) * t                                   # spare inventory dm1..dmD
    dvars = [sp.Symbol("d%d" % i) for i in range(t - 2, -1, -1)]
    mvars = [sp.Symbol("dm%d" % k) for k in range(1, D + 1)]
    vs = dvars + mvars
    # the u-grading: w(s_i) = t - i for EVERY coefficient, positive or negative
    w = [t - int(str(v)[1:]) for v in dvars] + [t + k for k in range(1, D + 1)]
    gens = {j: g for j, g in r["Gs"].items() if j != r["jphi"]}   # Phi-free rows
    out = dict(a=a, b=b, t=t, D=D, vars=vs, w=w, gens=gens, dvars=dvars,
               mvars=mvars, gw={j: b * t + j for j in gens},
               jphi=r["jphi"], M=r["M"], raw=r)
    _SYS[key] = out
    return out


def monos(weights: list[int], W: int) -> list[tuple]:
    """Every exponent tuple over `weights` of exact total weight W."""
    if W < 0:
        return []
    n = len(weights)
    res, cur = [], [0] * n

    def rec(i, rem):
        if rem == 0:
            res.append(tuple(cur))
            return
        if i == n:
            return
        for k in range(rem // weights[i] + 1):
            cur[i] = k
            rec(i + 1, rem - k * weights[i])
        cur[i] = 0
    rec(0, W)
    return res


def to_dict(expr, vs) -> dict:
    P = sp.Poly(sp.expand(expr), *vs)
    out = {}
    for mon, c in P.terms():
        c = sp.Rational(c)
        out[tuple(mon)] = Fraction(int(c.p), int(c.q))
    return out


def shift(d: dict, expo: tuple) -> dict:
    return {tuple(x + y for x, y in zip(m, expo)): c for m, c in d.items()}


def mul(d1: dict, d2: dict) -> dict:
    out = {}
    for m1, c1 in d1.items():
        for m2, c2 in d2.items():
            k = tuple(x + y for x, y in zip(m1, m2))
            out[k] = out.get(k, Fraction(0)) + c1 * c2
    return {k: v for k, v in out.items() if v}


class Reducer:
    """Incremental exact Gaussian elimination over Q on sparse dict-vectors."""

    def __init__(self) -> None:
        self.piv: dict = {}

    def reduce(self, vec: dict):
        vec = {k: v for k, v in vec.items() if v}
        while vec:
            k = max(vec)
            if k not in self.piv:
                return vec, k
            row = self.piv[k]
            f = vec[k] / row[k]
            for kk, vv in row.items():
                nv = vec.get(kk, Fraction(0)) - f * vv
                if nv:
                    vec[kk] = nv
                elif kk in vec:
                    del vec[kk]
        return {}, None

    def add(self, vec: dict) -> bool:
        red, k = self.reduce(vec)
        if k is None:
            return False
        self.piv[k] = red
        return True


def ideal_piece(S: dict, W: int) -> list:
    """A spanning set for I_W, I = (Phi-free generators).  Exact, complete:
    the ideal is u-weight homogeneous, so its weight-W piece is spanned by
    {monomial * G_j} with w(monomial) = W - w(G_j)."""
    vs, w = S["vars"], S["w"]
    if "gdict" not in S:
        S["gdict"] = {j: to_dict(g, vs) for j, g in S["gens"].items()}
    cols = []
    for j, gdj in sorted(S["gdict"].items()):
        for expo in monos(w, W - S["gw"][j]):
            cols.append(shift(gdj, expo))
    return cols


def hankel_minors(S: dict, extra_state: int = 0) -> tuple:
    """2x2 minors s_i*s_j - s_k*s_l (i+j=k+l) of the Hankel/catalecticant matrix
    on the coefficient sequence of s.  `extra_state = 0` is the window-only
    Hankel [[e,R,S,..],[R,S,..]] used at (72,108); extra_state = r also admits
    the r deepest state coefficients d_0, d_1, ... as entries."""
    t, D = S["t"], S["D"]
    # nu = -i indexes the coefficient s_{-nu}: nu >= 1 is the window (dm_nu),
    # nu = 0 is d_0, nu = -1 is d_1, ...  extra_state = 0 is window-ONLY.
    nus = [nu for nu in range(1 - extra_state, D + 1) if nu > 0 or -nu <= t - 2]
    var = {nu: (sp.Symbol("d%d" % (-nu)) if nu <= 0 else sp.Symbol("dm%d" % nu))
           for nu in nus}
    mn, mw = {}, {}
    for i in nus:
        for j in nus:
            if j < i + 2:
                continue
            for k in nus:
                if k <= i:
                    continue
                l = i + j - k
                if l not in var or l < k or l >= j:
                    continue
                mn[(i, j, k, l)] = var[i] * var[j] - var[k] * var[l]
                mw[(i, j, k, l)] = 2 * t + i + j
    return mn, mw


def e_power(S: dict, k: int) -> dict:
    nd = len(S["dvars"])
    return {tuple(k if i == nd else 0 for i in range(len(S["vars"]))): Fraction(1)}


def product_search(S: dict, extra_state: int = 0, kmax: int = 14) -> dict:
    """For every u-weight-admissible k, decide EXACTLY whether some Q-linear
    combination of {minor_1 * minor_2} is congruent to e^k modulo the Phi-free
    generators.  Returns {k: (n_pairs, n_Icols, hit)} over admissible k only."""
    vs, w, t = S["vars"], S["w"], S["t"]
    mn, mw = hankel_minors(S, extra_state)
    mind = {k: to_dict(v, vs) for k, v in mn.items()}
    ks = sorted(mn)
    res = {}
    for k in range(1, kmax + 1):
        W = k * (t + 1)
        pairs = [(ks[i], ks[j]) for i in range(len(ks)) for j in range(i, len(ks))
                 if mw[ks[i]] + mw[ks[j]] == W]
        if not pairs:
            continue
        cols = ideal_piece(S, W)
        R = Reducer()
        for c in cols:
            R.add(c)
        for p1, p2 in pairs:
            R.add(mul(mind[p1], mind[p2]))
        _, kk = R.reduce(e_power(S, k))
        res[k] = (len(pairs), len(cols), kk is None)
    return res


def toric_ideal_search(S: dict, klist) -> dict:
    """Is e^k in I + J, J = the ideal of the 2x2 window-Hankel minors?
    Two coefficient rings: `win` (coefficients in the window variables only --
    this is the shape of `6*W*Z = e^5`, whose cofactor W is itself a minor) and
    `all` (coefficients allowed to involve the state variables d_i too)."""
    vs, w, t = S["vars"], S["w"], S["t"]
    nd = len(S["dvars"])
    mw_all = w[nd:]
    mn, mw = hankel_minors(S, 0)
    mind = {k: to_dict(v, vs) for k, v in mn.items()}
    out = {}
    for k in klist:
        W = k * (t + 1)
        cols = ideal_piece(S, W)
        row = {}
        for tag, weights, pad in (("win", mw_all, lambda ex: (0,) * nd + ex),
                                  ("all", w, lambda ex: ex)):
            R = Reducer()
            for c in cols:
                R.add(c)
            n = 0
            for key, md in mind.items():
                for expo in monos(weights, W - mw[key]):
                    R.add(shift(md, pad(expo)))
                    n += 1
            _, kk = R.reduce(e_power(S, k))
            row[tag] = (n, kk is None)
        out[k] = row
    return out


# =============================================================================
#  A.  What the G-system IS, parametrically:  G_j = [x^-j] p^(b/a)
# =============================================================================
def laurent_pow(coeffs: dict, top: int, alpha, depth: int) -> dict:
    """(sum coeffs[i] x^i)^alpha as a Laurent series at x = infinity, for a
    MONIC top term x^top, down to x^(top*alpha - depth)."""
    f = {top - pw: c for pw, c in coeffs.items() if pw != top}
    assert coeffs[top] == 1
    cur, out = {0: sp.Integer(1)}, {0: sp.Integer(1)}
    for k in range(1, depth + 1):
        nxt = {}
        for i, ci in cur.items():
            for j, cj in f.items():
                if i + j <= depth:
                    nxt[i + j] = nxt.get(i + j, 0) + ci * cj
        cur = nxt
        if not cur:
            break
        bk = sp.binomial(alpha, k)
        for i, ci in cur.items():
            out[i] = out.get(i, 0) + bk * ci
    return {sp.Rational(top) * alpha - i: sp.expand(c) for i, c in out.items()}


def root_presentation(a, b, t, pvals=None, alpha_num=None):
    """Test  G_j == [x^-j] p^(b/a)  with the window variables read off p^(1/a).

    p is monic of degree a*t with [x^(a t - 1)]p = 0.  `pvals` (a list of
    rationals) specialises p's free coefficients; None keeps them symbolic.
    `alpha_num` overrides the numerator b of the forcing exponent (mutation).
    Returns (all_match, per-generator residual-is-zero list)."""
    S = system(a, b, t)
    at, D, jphi = a * t, S["D"], S["jphi"]
    ps = sp.symbols("pp0:%d" % at)
    coeffs = {at: sp.Integer(1)}
    for i in range(at - 1):
        coeffs[i] = ps[i] if pvals is None else sp.Rational(pvals[i])
    depth = max(S["M"], at) + 2
    srt = laurent_pow(coeffs, at, sp.Rational(1, a), depth)
    num = b if alpha_num is None else alpha_num
    sb = laurent_pow(coeffs, at, sp.Rational(num, a), depth)
    sub = {sp.Symbol("d%d" % (t - i)): sp.expand(srt.get(t - i, sp.Integer(0)))
           for i in range(2, t + 1)}
    for k in range(1, D + 1):
        sub[sp.Symbol("dm%d" % k)] = sp.expand(srt.get(-k, sp.Integer(0)))
    shift_zero = sp.expand(srt.get(t - 1, sp.Integer(0)))
    oks = []
    for j in sorted(S["raw"]["Gs"]):
        body = S["raw"]["Gs"][j] - (PHI if j == jphi else 0)
        oks.append(sp.expand(sp.expand(body.xreplace(sub))
                             - sp.expand(sb.get(-j, sp.Integer(0)))) == 0)
    return all(oks), oks, shift_zero


def group_A(L: Ledger) -> None:
    L.head("A.  The general G-system:  G_j = [x^-j] p^(b/a),  window = p^(1/a)")

    # ---- A1: the builder's (72,108) control -------------------------------
    r = gsys.build_gsystem(2, 3, 4, 1, 1)
    pub = gsys.published_72108()
    L.ck("A1 builder reproduces (72,108)",
         all(sp.expand(r["Gs"][j] - pub[j]) == 0 for j in (1, 2, 3, 5)),
         "g_system_75_125.build_gsystem(2,3,4) == published_72108(); "
         "(m,n) = (a,b) = (2,3), chart exponent t = 4")

    # ---- A2: the closed form, symbolically --------------------------------
    ok234, _, sh234 = root_presentation(2, 3, 4)
    ok235, _, sh235 = root_presentation(2, 3, 5)
    L.ck("A2 closed form exact at (2,3,4),(2,3,5)", ok234 and ok235,
         "p monic of degree a*t, [x^(at-1)]p = 0; s := p^(1/a); "
         "d_i = [x^i]s, dm_k = [x^-k]s; G_j = [x^-j] p^(b/a) -- residual 0")
    L.ck("A2 the x-shift d_{t-1} = 0 is automatic",
         sh234 == 0 and sh235 == 0,
         "[x^(t-1)] p^(1/a) = p_{at-1}/a = 0")
    bad = []
    for (a, b, t) in ((2, 3, 4), (2, 3, 5)):
        for num in (b + 1, b - 1, b + a):
            okm, _, _ = root_presentation(a, b, t, alpha_num=num)
            bad.append(not okm)
    L.mut("A2 wrong forcing exponent dies", all(bad),
          "p^((b+1)/a), p^((b-1)/a), p^((b+a)/a) all leave a residual "
          "(%d/%d mutants die)" % (sum(bad), len(bad)))

    # ---- A3: the same closed form at (75,125), at exact rational points ----
    pts = [[sp.Rational(n * n % 17 - 8, (n % 5) + 1) for n in range(i, i + 11)]
           for i in (1, 7, 23, 41)]
    ok354 = [root_presentation(3, 5, 4, pvals=p)[0] for p in pts]
    L.ck("A3 closed form at (3,5,4), 4 exact points", all(ok354),
         "all eight generators G1..G7,G9 match [x^-j] p^(5/3) exactly over Q "
         "at 4 independent rational p-vectors  [CHECKED-at-points, not symbolic:"
         " the 12-symbol expansion is minutes-scale]")
    badpt = [not root_presentation(3, 5, 4, pvals=pts[0], alpha_num=n)[0]
             for n in (4, 6, 8)]
    L.mut("A3 wrong exponent dies at the points too", all(badpt),
          "p^(4/3), p^(6/3), p^(8/3) all fail at the same point")

    # ---- A4: the structural counts ----------------------------------------
    rows, homog, hetero = [], [], []
    for (a, b, t) in ((2, 3, 4), (2, 3, 5), (3, 4, 4), (3, 5, 4)):
        S = system(a, b, t)
        D = (a - 1) * t
        rows.append(S["D"] == D and S["jphi"] == D + 1
                    and sorted(S["gens"]) == list(range(1, D)))
        # VERIFY homogeneity: every monomial of G_j has u-weight exactly b*t+j,
        # under w(s_i) = t - i.  (This is what P1's weight law rests on.)
        for j, g in S["gens"].items():
            ws = {sum(ww * ex for ww, ex in zip(S["w"], mon))
                  for mon, _ in sp.Poly(g, *S["vars"]).terms()}
            homog.append(ws == {b * t + j})
            # control: a DIFFERENT weight assignment must NOT homogenise them
            w2 = [x + 1 for x in S["w"]]
            ws2 = {sum(ww * ex for ww, ex in zip(w2, mon))
                   for mon, _ in sp.Poly(g, *S["vars"]).terms()}
            hetero.append(len(ws2) > 1)
    L.ck("A4 counts + VERIFIED u-homogeneity",
         all(rows) and all(homog) and len(homog) == 3 + 4 + 7 + 7,
         "D=(a-1)t spares, j_Phi = D+1, Phi-free rows are G_1..G_{D-1}; and "
         "every monomial of every one of the %d generators checked has u-weight "
         "exactly b*t+j under w(s_i) = t-i" % len(homog))
    L.mut("A4 the grading is the RIGHT one", all(hetero),
          "shifting every weight by +1 destroys homogeneity in all %d "
          "generators -- so w(s_i) = t-i is forced, not one of many gradings"
          % len(hetero))
    L.note("so `b` = n enters ONLY through the exponent b/a and the weight "
           "shift b*t; the NUMBER of equations depends on (a,t) alone")


# =============================================================================
#  B.  The (2,3) positive control, and the weight law that fixes the exponent
# =============================================================================
def group_B(L: Ledger) -> None:
    L.head("B.  (2,3,4) control: 6*W*Z = e^5 re-derived, constant and all")

    S = system(2, 3, 4)
    vs = S["vars"]
    e, R, Sv, T = S["mvars"]
    W, Z = e * Sv - R ** 2, e * T - R * Sv

    # ---- B1: the exponent 5 is the ONLY weight-admissible one --------------
    mn, mw = hankel_minors(S, 0)
    prods = sorted({mw[i] + mw[j] for i in mn for j in mn})
    adm = sorted({p for p in prods if p % (S["t"] + 1) == 0})
    L.ck("B1 unique admissible exponent at (2,3,4)",
         adm == [25] and [p // 5 for p in adm] == [5],
         "minor weights %s; products span %d..%d; only 25 = 5*w(e) is a "
         "multiple of w(e) = 5, so k = 5 is forced BEFORE any algebra"
         % (sorted(set(mw.values())), min(prods), max(prods)))

    # ---- B2: the constant is DETERMINED, and comes out 6 -------------------
    cols = ideal_piece(S, 25)
    Rd = Reducer()
    for c in cols:
        Rd.add(c)
    nf_e5, k1 = Rd.reduce(e_power(S, 5))
    nf_wz, k2 = Rd.reduce(to_dict(sp.expand(W * Z), vs))
    L.ck("B2 e^5 is NOT in I alone", k1 is not None,
         "so the relation is not vacuous: e^5 needs the minor product")
    const = None
    if k2 is not None and nf_wz:
        ratios = {nf_e5.get(kk, Fraction(0)) / vv for kk, vv in nf_wz.items()}
        if len(ratios) == 1 and all(nf_e5.get(kk, Fraction(0)) == r * vv
                                    for kk, vv in nf_wz.items() for r in ratios) \
                and set(nf_e5) <= set(nf_wz):
            const = ratios.pop()
    L.ck("B2 the constant is FORCED and equals 6", const == Fraction(6),
         "e^5 == c * W*Z modulo I has a unique solution c, computed by exact "
         "normal-form comparison: c = %s  (NOT typed in)" % const)

    # ---- B3: recover the cofactors, also without typing them ---------------
    wts = {v: ww for v, ww in zip(vs, S["w"])}
    cof = {}
    unk = []
    for j, g in sorted(S["gens"].items()):
        need = 25 - S["gw"][j]
        terms = [e ** i * R ** ((need - 5 * i) // 6)
                 for i in range(need // 5 + 1) if (need - 5 * i) % 6 == 0]
        cs = sp.symbols("q%d_0:%d" % (j, max(1, len(terms))))
        cof[j] = sum(c * tm for c, tm in zip(cs, terms))
        unk += list(cs[:len(terms)])
    expr = sp.expand(sum(cof[j] * g for j, g in S["gens"].items())
                     - (6 * W * Z - e ** 5))
    sol = sp.solve([c for _, c in sp.Poly(expr, *vs).terms()], unk, dict=True)
    got = ({j: sp.expand(cof[j].xreplace(sol[0]).xreplace({u: 0 for u in unk}))
            for j in cof} if sol else None)
    ok = got is not None and sp.expand(
        sum(got[j] * g for j, g in S["gens"].items()) - (6 * W * Z - e ** 5)) == 0
    L.ck("B3 cofactors recovered by solving, not typing", ok,
         "restricting cofactors to monomials in (e,R) gives "
         "%s" % ({("G%d" % j): got[j] for j in sorted(got)} if got else None))
    L.note("that is exactly TORIC_SYZYGY.md R2:  2e^2 G3 - 4eR G2 + 2R^2 G1 "
           "= 6WZ - e^5")

    # ---- B3 mutation: drop a row, and perturb the target -------------------
    S3 = dict(S)
    S3["gens"] = {j: g for j, g in S["gens"].items() if j != 3}
    S3.pop("gdict", None)
    R3 = Reducer()
    for c in ideal_piece(S3, 25):
        R3.add(c)
    R3.add(to_dict(sp.expand(W * Z), vs))
    _, kk3 = R3.reduce(e_power(S, 5))
    Rb = Reducer()
    for c in cols:
        Rb.add(c)
    Rb.add(to_dict(sp.expand(W * W), vs))
    Rb.add(to_dict(sp.expand(Z * Z), vs))
    _, kk4 = Rb.reduce(e_power(S, 5))
    L.mut("B3 dropping G3 / swapping the minors kills it",
          kk3 is not None and kk4 is not None,
          "without G3 there is no certificate; and 6W^2/6Z^2 in place of 6WZ "
          "fails -- the search is discriminating, not a rank accident")

    # ---- B4: THE WEIGHT LAW -- the exponent is a function of t -------------
    # w(W) = 2t+4, w(Z) = 2t+5, w(e) = t+1.  c*W*Z = e^k needs (t+1) | (4t+9).
    sols = [tt for tt in range(2, 200) if (4 * tt + 9) % (tt + 1) == 0]
    L.ck("B4 (t+1) | (4t+9) has the unique solution t = 4",
         sols == [4] and (4 * 4 + 9) // 5 == 5,
         "4t+9 = 4(t+1)+5, so (t+1)|5, so t+1 = 5.  The exponent of the "
         "(72,108) identity is (4t+9)/(t+1) = 5 -- a function of t ALONE. "
         "m and n do not appear.")
    many = [tt for tt in range(2, 200) if (4 * tt + 16) % (tt + 1) == 0]
    none_ = [tt for tt in range(2, 200) if (4 * tt + 4) % (tt + 1) == 0]
    L.mut("B4 the divisibility test is not vacuous",
          len(many) >= 3 and len(none_) == 198,
          "the same test on numerator 4t+16 has %d solutions %s, and on 4t+4 "
          "it is satisfied by EVERY t -- so `unique solution at t=4` is a "
          "property of the numerator 4t+9, not of the test"
          % (len(many), many))

    # ---- B5: the coincidence ledger ---------------------------------------
    m, n, t = 2, 3, 4
    coincide = {"m+n": m + n, "t+1": t + 1, "(4t+9)/(t+1)": (4 * t + 9) // (t + 1),
                "2n-1": 2 * n - 1, "m*n-1": m * n - 1}
    L.ck("B5 five formulas all give 5 at (72,108)",
         set(coincide.values()) == {5},
         "%s -- (72,108) cannot distinguish them" % coincide)
    m2, n2, t2 = 3, 5, 4
    sep = {"m+n": m2 + n2, "t+1": t2 + 1,
           "(4t+9)/(t+1)": sp.Rational(4 * t2 + 9, t2 + 1),
           "2n-1": 2 * n2 - 1, "m*n-1": m2 * n2 - 1}
    L.ck("B5 they separate at (3,5)", len(set(sep.values())) >= 4,
         "%s -- so (75,125) IS a discriminating test" % sep)


# =============================================================================
#  C.  The (3,5) test -- and the (2,3,t) t-sweep that isolates the cause
# =============================================================================
def group_C(L: Ledger, fast: bool) -> None:
    L.head("C.  (3,5) at t=4: NO product-of-minors identity, at ANY exponent")

    S = system(3, 5, 4)
    mn, mw = hankel_minors(S, 0)
    prods = sorted({mw[i] + mw[j] for i in mn for j in mn})
    adm = sorted({p // 5 for p in prods if p % 5 == 0})
    L.ck("C1 admissible exponents at (3,5,4) are exactly {5,6,7,8}",
         adm == [5, 6, 7, 8],
         "minor weights %d..%d, products %d..%d; k*w(e) = 5k lands in range "
         "only for k = 5..8.  m+n = 8 IS admissible -- the prediction gets a "
         "fair test" % (min(mw.values()), max(mw.values()),
                        min(prods), max(prods)))

    res = product_search(S, extra_state=0)
    L.ck("C2 every admissible k FAILS at (3,5,4)",
         set(res) == {5, 6, 7, 8} and not any(v[2] for v in res.values()),
         "exact rank over Q: %s"
         % {k: "%d pairs, |I_W|=%d -> %s" % (v[0], v[1], "HIT" if v[2] else "no")
            for k, v in sorted(res.items())})
    L.ck("C2a in particular  X*Y = e^(m+n) = e^8 is FALSE",
         res.get(8, (0, 0, True))[2] is False,
         "all %d weight-admissible minor pairs at u-weight 40 tested "
         "simultaneously (a single linear-span test, so no pair is missed)"
         % res.get(8, (0,))[0])
    L.ck("C2b and the direct analogue  c*W*Z = e^5 is FALSE",
         res.get(5, (0, 0, True))[2] is False,
         "W*Z is the UNIQUE weight-25 pair, so this is the one candidate that "
         "the (72,108) shape maps onto at t = 4")

    # POSITIVE CONTROL through the identical code path
    res234 = product_search(system(2, 3, 4), extra_state=0)
    L.ck("C3 POSITIVE CONTROL: same routine finds (2,3,4)",
         [k for k, v in res234.items() if v[2]] == [5],
         "product_search(2,3,4) -> admissible %s, HIT at k=5 (%d pair, "
         "|I_W|=%d).  So C2's negative is not a broken search"
         % (sorted(res234), res234[5][0], res234[5][1]))

    # ---- C4: robustness -- widen the minor family -------------------------
    wide = {}
    for extra in (1, 2):
        wide[extra] = product_search(S, extra_state=extra)
    ctrl = {extra: product_search(system(2, 3, 4), extra_state=extra)
            for extra in (1, 2)}
    L.ck("C4 negative survives widening the Hankel",
         all(not v[2] for r in wide.values() for v in r.values()),
         "admitting the state coefficient d_0 (and d_0,d_1) as Hankel entries "
         "raises the candidate count to %d/%d pairs and still finds nothing"
         % (sum(v[0] for v in wide[1].values()),
            sum(v[0] for v in wide[2].values())))
    L.ck("C4 widened control still finds (2,3,4)",
         all(any(v[2] for v in ctrl[x].values()) for x in (1, 2)),
         "the widened routine still HITs at (2,3,4) k=5")

    # ---- C5: the sharper and the looser toric shapes ----------------------
    tor = toric_ideal_search(S, [3, 4, 5, 6, 7])
    tor234 = toric_ideal_search(system(2, 3, 4), [3, 4, 5])
    L.ck("C5 no window-coefficient toric relation at (3,5,4), k<=7",
         all(not v["win"][1] for v in tor.values()),
         "e^k = (toric-ideal element with coefficients in dm1..dm8) fails for "
         "k = 5,6,7  [%s]"
         % {k: v["win"][0] for k, v in sorted(tor.items())})
    L.ck("C5 control: it DOES hold at (2,3,4), k=5",
         tor234[5]["win"][1] and not tor234[4]["win"][1],
         "the (72,108) identity is exactly this shape (its cofactor W is "
         "itself a minor); minimal k = 5")
    L.ck("C5 what (3,5,4) DOES satisfy: e^6, with state coefficients",
         tor[6]["all"][1] and not any(tor[k]["all"][1] for k in (3, 4, 5)),
         "e^6 in I + J once the toric cofactors may involve d0,d1,d2 -- and "
         "k = 3,4,5 all fail, so 6 is MINIMAL, and 6 != 8 = m+n")
    L.ck("C5 (2,3,4) counterpart is e^4, not e^5",
         tor234[4]["all"][1] and not tor234[3]["all"][1],
         "so even this looser invariant is 4 at (2,3) and 6 at (3,5): "
         "neither m+n (5,8) nor any shared formula is reproduced")

    # ---- C6: the t-sweep at fixed (m,n) = (2,3) ---------------------------
    tlist = (3, 4, 5) if fast else (3, 4, 5, 6)
    sweep = {}
    for tt in tlist:
        sweep[tt] = product_search(system(2, 3, tt), extra_state=0)
    hits = {tt: [k for k, v in r.items() if v[2]] for tt, r in sweep.items()}
    L.ck("C6 at fixed (m,n)=(2,3) the identity exists ONLY at t=4",
         hits.get(4) == [5] and all(not hits[tt] for tt in tlist if tt != 4),
         "t -> hits: %s.  Same m, same n, identity gone.  The exponent -- and "
         "the identity itself -- is a fact about t, not about (m,n)" % hits)
    L.mut("C6 the sweep is not vacuously empty",
          all(sum(v[0] for v in sweep[tt].values()) > 0 for tt in tlist),
          "every rung had weight-admissible candidates to test: %s"
          % {tt: sum(v[0] for v in sweep[tt].values()) for tt in tlist})


# =============================================================================
#  D.  The (m,n) sweep at t = 4:  (2,3) is alone
# =============================================================================
def group_D(L: Ledger) -> None:
    L.head("D.  Sweep over (m,n) at t = 4 -- (2,3) is the only member")

    cases = [(2, 3), (2, 5), (3, 4), (3, 5)]
    table = {}
    for (a, b) in cases:
        if math.gcd(a, b) != 1 or b <= a:
            continue
        r = product_search(system(a, b, 4), extra_state=0)
        table[(a, b)] = sorted(k for k, v in r.items() if v[2])
    L.ck("D1 only (m,n) = (2,3) admits the identity",
         table.get((2, 3)) == [5]
         and all(not v for kk, v in table.items() if kk != (2, 3)),
         "t = 4, gcd(m,n)=1: %s" % {("m=%d,n=%d" % k): (v or "none")
                                    for k, v in sorted(table.items())})
    L.ck("D2 the constant `6` has no second data point",
         sum(1 for v in table.values() if v) == 1,
         "exactly one case in the whole (m,n,t) sweep carries an identity, so "
         "any formula for the constant -- m*n*(n-m) = 6, or 2*c_2 = 6 -- fits "
         "one point and is UNFALSIFIABLE.  Reported as such, not as a law.")


# =============================================================================
#  E.  lam at (75,125):  the (H-cap) gate
# =============================================================================
def group_E(L: Ledger) -> None:
    L.head("E.  lam at (75,125) -- the CONTACT_LEMMA (H-cap) gate `lam >= m`")

    l72 = wf.window_law(204, 17, 238)
    L.ck("E1 (72,108): lam = 2 = m, exactly at equality",
         l72["lam"] == 2,
         "lam := deg_slope - W_step = (deg_y(Phi) - ord_y(Phi))/M "
         "= (238-204)/17 = %s; m = 2, so (H-cap) holds with ZERO margin"
         % l72["lam"])

    fam = wf.family(3)                       # the repaired F2 rung a = 3
    l75 = wf.window_law(fam["ordPhi"], fam["M"], fam["degPhi"])
    L.ck("E2 repaired (75,125) inputs come through the guard",
         (fam["t"], fam["kappa"], fam["b"], fam["M"], fam["ordPhi"],
          fam["degPhi"]) == (4, 2, 5, 29, 80, 80),
         "t=4, kappa=2, n=b=5, M=29, ord_y(Phi)=deg_y(Phi)=80 "
         "(Phi = (1/3) y^80 is a MONOMIAL because C = y is)")
    L.ck("E3 (75,125): lam = 0", l75["lam"] == 0,
         "(80 - 80)/29 = 0.  deg_C = ord_C = 1, so deg_y(Phi) - ord_y(Phi) "
         "= N*(deg_C - ord_C) = 0 identically")
    L.ck("E4 the (H-cap) gate `lam >= m` FAILS at (75,125)",
         not (l75["lam"] >= 3) and (l72["lam"] >= 2),
         "lam = 0 < 3 = m.  CONTACT_LEMMA.md sec.4.3's open flag closes "
         "NEGATIVE: the slice-cascade section is void at (75,125)")
    L.mut("E1-E3 the lam formula is not identically zero",
          l72["lam"] != 0 and wf.window_law(204, 17, 255)["lam"] != 0,
          "it returns 2 on the (72,108) signature and 3 on a perturbed one; "
          "the 0 at (75,125) is a property of Phi, not of the code")
    lams = [wf.window_law(wf.family(aa)["ordPhi"], wf.family(aa)["M"],
                          wf.family(aa)["degPhi"])["lam"] for aa in range(2, 9)]
    L.ck("E5 lam = 0 for the whole F2 family, not just a = 3",
         all(x == 0 for x in lams),
         "rungs a = 2..8 all give lam = 0: a family fact (C monomial), not an "
         "a = 3 accident")


# ------------------------------------------------------------------- main ----
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--fast", action="store_true",
                    help="skip group D and the t=6 rung of C6")
    args = ap.parse_args()
    L = Ledger(args.quiet)
    t0 = time.time()
    group_A(L)
    group_B(L)
    group_C(L, args.fast)
    if not args.fast:
        group_D(L)
    group_E(L)
    if not args.quiet:
        print("\n(%.1f s)" % (time.time() - t0))
    return L.report()


if __name__ == "__main__":
    sys.exit(main())
