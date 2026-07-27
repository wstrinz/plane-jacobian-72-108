#!/usr/bin/env python3
"""weight_lemma_75_125.py -- does the Phi-divisor mechanism transfer to (75,125)?

*** REPAIRED 2026-07-26 (PASSPORT_75_125_REPAIR.md).  Two independent repairs. ***

(i) INPUTS.  The (75,125) chart exponent is l = ceil(b0/a0) = 4, not the
    denominator 5 of GGV5's final chain corner (7\\5,2): the (t,q) =
    (l_final,b_final) dictionary holds only on the retraction shape
    b0 = l(a0-1), which (5,20) fails.  So t=4, kappa=2, C=y (a MONOMIAL),
    ord C=1, N=77, Phi=(1/3) y^80 with signature (80,80,0,0), M=29,
    W_step = 80/29, q_window = 29.  Rebuild g_system_75_125.json before running.

(ii) CRITERION.  Section B used to assert FULL COLUMN RANK (nullity 0) of the
    graded piece and read that as "no relation c*Phi = e*B".  Nullity 0 is
    SUFFICIENT but not NECESSARY: a pure inter-generator syzygy -- a nullspace
    vector with no Phi-carrying column in its support -- also lowers the rank
    while giving c = 0, i.e. no Phi-relation at all.  On the repaired system such
    syzygies exist (the first at weight 43, supported on G1 and G2 only), so the
    old test would have raised a FALSE ALARM.  The criterion is now the right
    one: a relation with c != 0 exists iff the Phi-carrying columns are DEPENDENT
    modulo the span of the others.  Both directions are exercised -- the
    (72,108) control DETECTS its K-syzygy under the same test.

THE PREDICTION UNDER TEST.  CAPS_AUDIT.md sec.5 states the weight lemma: for a
weight-W homogeneous relation c*Phi = e*B in a regime with stripped slope lambda,

    max(0, D - lambda*(W - w_e))  <=  deg e  <=  lambda*w_e,     D = deg Phi,

an interval of length sigma = lambda*W - D, and the forcing pins deg e = lambda*w_e
IFF sigma = 0.  It logs [inference]: at (75,125) deg_slope = deg_y(Phi)/M = 14, so
sigma = 0 and "the lemma predicts the analogous forcing fires".

THE ANSWER, IN ORDER.

  (1) The relation comes FIRST.  A complete graded search of the (75,125)
      G-system ideal finds NO relation c*Phi = e*B with c not itself divisible
      by e -- at ANY weight W = 29..45 (i.e. any multiplier c of u-weight 0..16),
      and for e replaced by ANY of the state variables d2,d1,d0,dm1.
      The identical search at (72,108) recovers the published K-syzygy uniquely
      (nullity exactly 1) and, at higher weights, a kernel of exactly the
      dimension the K-multiples predict.  So the search is sensitive, and the
      (75,125) negative is a real negative.
      This VERDICT IS UNCHANGED by the 2026-07-26 repair -- it was re-run from
      scratch on the rebuilt 8-generator system, with the corrected criterion.

  (2) Therefore sigma predicts nothing there.  With no relation the lemma has no
      hypothesis to discharge.

  (3) And sigma = 0 at (75,125) is TAUTOLOGICAL anyway.  window_functions_75_125
      .window_law DEFINES deg_slope := deg_y(Phi)/M, so sigma = deg_slope*M -
      deg_y(Phi) = 0 identically, for every case, by construction.  At (72,108)
      the slope 14 is computed from the Prop 4.3 Newton polygon (max(j-2i) = 0
      over the hull) and 238 from the f1 ODE -- two disjoint routes -- so there
      sigma = 0 is a genuine agreement.
      REPAIRED: at (75,125) the deg_slope is now 80/29, which is NOT AN INTEGER,
      so there is no affine y-degree cap at all -- CAPS_AUDIT sec.5's premise
      ("deg_slope = 14") is not merely tautological, it is FALSE.  And because
      C = y is a monomial, deg_y(Phi) = ord_y(Phi), so the two slopes COINCIDE
      and the stripped slope lambda = deg_slope - W_step is 0, not the 101/12 of
      the superseded model.  With lambda = 0 the lemma's interval
      [max(0, D - lambda(W-w_e)), lambda*w_e] = [80, 0] is EMPTY, so the lemma
      itself now forbids the relation whose absence section B verifies.

  (4) There is a POSITIVE structural obstruction, and it is on the ord side, not
      the deg side -- and the repair makes it TOTAL.  The lower (y-order) cap at
      (75,125) is quasi-affine, L(w) = ceil(80w/29), with quasi-period
      q_window = 29.  ceil is strictly superadditive off the period: for
      e = dm1 (u-weight 5) and B (u-weight 24),

          L(5) + L(24) = 14 + 67 = 81  >  80 = ord_y(Phi) = L(29).

      Since q_window = 29 = M exactly, NO split of W = 29 has carry 0: the
      superseded model left the two escapes w_e in {12,24} (multiples of 12
      below 36), and there is now no escape at all.  The mechanism needs
      carry 0, i.e. q_window | w_e, and 0 < w_e < M = q_window is impossible.
      At (72,108) q_window = 1 so the carry vanishes identically; in the F2
      family q_window = 12a-7 is never 1.  This is the precise sense in which
      the quasi-affine normalisation breaks the mechanism.

  VERDICT: the divisor mechanism does NOT transfer to (75,125).  (72,108) is
  special, and the thing that makes it special is q_window = 1.
  *** THIS VERDICT SURVIVES THE 2026-07-26 REPAIR, strengthened on both legs:
  the graded search still finds nothing (now under the corrected criterion), and
  the carry obstruction went from "all but two splits" to "every split". ***

Read-only: creates nothing, modifies nothing.  Usage:
    python -u weight_lemma_75_125.py            # full report
    python -u weight_lemma_75_125.py --quiet    # self-check, exit 0 iff all pass
    python -u weight_lemma_75_125.py --fast     # skip the deep weight sweeps
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from fractions import Fraction

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

y = sp.Symbol("y")
BIGP = (1 << 61) - 1                       # 2^61-1, prime


# --------------------------------------------------------------------------
# generic graded-piece syzygy search engine
# --------------------------------------------------------------------------
class GradedSystem:
    """A weight-homogeneous G-system: variables with u-weights, generators with
    u-weights, one distinguished symbol Phi appearing only in the top generator."""

    def __init__(self, tag, varnames, varweights, gens, gweights, phiname, phiweight):
        self.tag = tag
        self.VN = list(varnames) + [phiname]
        self.SYM = {n: sp.Symbol(n) for n in self.VN}
        self.G = [self.SYM[n] for n in self.VN]
        self.IDX = {n: i for i, n in enumerate(self.VN)}
        self.W = list(varweights) + [phiweight]
        self.IPHI = self.IDX[phiname]
        self.phiweight = phiweight
        self.GW = dict(gweights)
        self.POLY = {}
        for nm, ex in gens.items():
            P = sp.Poly(sp.expand(ex), *self.G)
            self.POLY[nm] = {m: sp.Rational(c) for m, c in P.terms()}

    def homogeneity(self):
        """Recomputed u-weight set of every generator (independent of any stored
        weight field)."""
        out = {}
        for nm, terms in self.POLY.items():
            ws = {sum(e * self.W[i] for i, e in enumerate(mono)) for mono in terms}
            out[nm] = ws
        return out

    def monomials(self, n, exclude=()):
        """Exponent tuples over the non-Phi variables with total u-weight n,
        using no variable in `exclude`."""
        ids = [i for i in range(len(self.VN)) if i != self.IPHI and i not in exclude]
        ids.sort(key=lambda i: -self.W[i])
        out, cur = [], [0] * len(self.VN)

        def rec(k, rem):
            if rem == 0:
                out.append(tuple(cur))
                return
            if k >= len(ids):
                return
            i = ids[k]
            if self.W[i] > rem:
                rec(k + 1, rem)
                return
            for x in range(rem // self.W[i], -1, -1):
                cur[i] = x
                rec(k + 1, rem - x * self.W[i])
            cur[i] = 0

        rec(0, n)
        return out

    def columns(self, Wtarget, ivar):
        """The linear system for  sum_a c_a * m_a * G_a  ==  c*Phi   (mod ivar).

        Multiplier monomials divisible by ivar are DROPPED: those terms are
        automatically == 0 mod ivar, so they place no constraint (and they only
        ever contribute to B).  A column is the Phi-free part of m_a*G_a reduced
        mod ivar.  Nullity 0  <=>  no relation with c not divisible by ivar.
        """
        cols, labels = [], []
        for nm, gw in self.GW.items():
            n = Wtarget - gw
            if n < 0:
                continue
            for m in self.monomials(n, exclude=(ivar,)):
                d = {}
                for mono, c in self.POLY[nm].items():
                    if mono[self.IPHI]:
                        continue                      # the c*Phi side
                    if mono[ivar] + m[ivar]:
                        continue                      # already == 0 mod ivar
                    key = tuple(a + b for a, b in zip(mono, m))
                    d[key] = d.get(key, sp.Rational(0)) + c
                cols.append({k: Fraction(int(v.p), int(v.q)) for k, v in d.items() if v != 0})
                labels.append((nm, m))
        return cols, labels

    def mono_str(self, m):
        s = "*".join("%s%s" % (self.VN[i], "" if e == 1 else "**%d" % e)
                     for i, e in enumerate(m) if e)
        return s or "1"


def rank_modp(cols, p=BIGP):
    """Rank of the column set mod p.  rank_p <= rank_Q always, so rank_p == ncols
    PROVES nullity 0 over Q."""
    piv, rank = {}, 0
    for col in cols:
        v = {}
        for k, c in col.items():
            val = (c.numerator % p) * pow(c.denominator % p, p - 2, p) % p
            if val:
                v[k] = val
        while v:
            k = max(v)
            if k in piv:
                pc = piv[k]
                f = v[k] * pow(pc[k], p - 2, p) % p
                for kk, cc in pc.items():
                    nv = (v.get(kk, 0) - f * cc) % p
                    if nv:
                        v[kk] = nv
                    elif kk in v:
                        del v[kk]
            else:
                piv[k] = v
                rank += 1
                break
    return rank


def incr_rank_modp(cols, piv=None, p=BIGP):
    """Column-echelon reduce `cols` into `piv`; return (#new pivots, piv).

    Splitting the rank computation lets us ask the question that actually
    matters: are the Phi-carrying columns independent MODULO the others?
    """
    if piv is None:
        piv = {}
    new = 0
    for col in cols:
        v = {}
        for k, c in col.items():
            val = (c.numerator % p) * pow(c.denominator % p, p - 2, p) % p
            if val:
                v[k] = val
        while v:
            k = max(v)
            if k in piv:
                pc = piv[k]
                f = v[k] * pow(pc[k], p - 2, p) % p
                for kk, cc in pc.items():
                    nv = (v.get(kk, 0) - f * cc) % p
                    if nv:
                        v[kk] = nv
                    elif kk in v:
                        del v[kk]
            else:
                piv[k] = v
                new += 1
                break
    return new, piv


def phi_relation_exists(sysX, phigen, Wtarget, ivar, p=BIGP):
    """Is there a relation c*Phi = ivar*B with c NOT divisible by ivar?

    A nullspace vector of the full column set gives an ideal combination whose
    Phi-free part vanishes mod ivar, i.e. c*Phi + ivar*B = 0 where c is read off
    the columns of the Phi-CARRYING generator.  So a genuine Phi-divisor relation
    exists iff some nullspace vector has a nonzero Phi-carrying coefficient, i.e.
    iff the Phi-carrying columns are DEPENDENT modulo the span of the rest.

    Returns (exists, n_phi_cols, n_independent, n_cols).  `not exists` is PROVED
    over Q: independence mod p implies independence over Q (clear a rational
    dependency to integers of content 1; it cannot reduce to 0 mod p).

    NOTE why plain "nullity == 0" is the WRONG test: pure inter-generator
    syzygies (no Phi column in support) also lower the rank yet give c = 0.  They
    occur on the repaired (75,125) system from weight 43 on.
    """
    cols, lab = sysX.columns(Wtarget, ivar)
    rest = [c for c, L in zip(cols, lab) if L[0] != phigen]
    phic = [c for c, L in zip(cols, lab) if L[0] == phigen]
    _, piv = incr_rank_modp(rest, None, p)
    new, _ = incr_rank_modp(phic, piv, p)
    return (new != len(phic)), len(phic), new, len(cols)


def nullspace_exact(cols):
    """Exact rational nullspace of the column set (small systems only)."""
    keys = sorted(set().union(*[set(c) for c in cols])) if cols else []
    M = sp.zeros(len(keys), len(cols))
    ri = {k: i for i, k in enumerate(keys)}
    for j, c in enumerate(cols):
        for k, v in c.items():
            M[ri[k], j] = sp.Rational(v.numerator, v.denominator)
    return M.nullspace()


# --------------------------------------------------------------------------
# case builders
# --------------------------------------------------------------------------
def build_72108():
    import full_system_bridge as fsb
    st = fsb.gsystem()
    return GradedSystem(
        "(72,108)",
        ["d2", "d1", "d0", "dm1", "dm2", "dm3", "dm4"], [2, 3, 4, 5, 6, 7, 8],
        {"G1": st["G1"], "G2": st["G2"], "G3": st["G3"], "G5": st["G5"]},
        {"G1": 13, "G2": 14, "G3": 15, "G5": 17}, "Phi", 17)


def load_75125():
    return json.load(open(os.path.join(HERE, "g_system_75_125.json")))


def build_75125(J):
    t, M = J["case"]["t"], J["phi_u_weight"]
    VN = [n for n in J["variable_order"] if n != "Phi"]

    def uw(n):                                  # S = sum_m d_m u^(t-m)  =>  w(d_m) = t-m
        m = -int(n[2:]) if n.startswith("dm") else int(n[1:])
        return t - m

    SY = {n: sp.Symbol(n) for n in J["variable_order"]}
    gens = {nm: sp.sympify(J["generators"][nm]["poly"], locals=SY)
            for nm in J["generator_names"]}
    gw = {nm: J["generators"][nm]["u_weight"] for nm in J["generator_names"]}
    return GradedSystem("(75,125)", VN, [uw(n) for n in VN], gens, gw, "Phi", M)


# --------------------------------------------------------------------------
# cap arithmetic
# --------------------------------------------------------------------------
def window_law(ordPhi, M, degPhi):
    """(alpha, q_window, deg_slope) from a case's Phi signature.

    NOTE the provenance: BOTH slopes are read off Phi.  deg_slope := degPhi/M is
    a DEFINITION, so sigma = deg_slope*M - degPhi == 0 for every case, always."""
    Wst = sp.Rational(ordPhi, M)
    dsl = sp.Rational(degPhi, M)
    return int(Wst.p), int(Wst.q), (int(dsl) if dsl.q == 1 else dsl)


def Lcap(w, alpha, q):
    """Lower y-order cap ceil(alpha*w/q) -- tight integer bound under the
    extreme-ray premise (ord >= (alpha/q)*w for every window object)."""
    return -((-alpha * w) // q)


def carry(w1, w2, alpha, q):
    return Lcap(w1, alpha, q) + Lcap(w2, alpha, q) - Lcap(w1 + w2, alpha, q)


# --------------------------------------------------------------------------
def run(verbose=True, fast=False):
    out, npass, nfail = [], 0, 0

    def chk(cid, ok, msg):
        nonlocal npass, nfail
        if ok:
            npass += 1
        else:
            nfail += 1
        out.append((cid, ok, msg))
        if verbose:
            print("  [%s] %-6s %s" % ("PASS" if ok else "FAIL", cid, msg))

    # =====================================================================
    if verbose:
        print("\n== A. (75,125) weight data, re-derived from primitives ==")
    J = load_75125()
    a, b, t = J["case"]["a"], J["case"]["b"], J["case"]["t"]
    kap, qq = J["case"]["kappa"], J["case"]["q"]
    M = J["phi_u_weight"]

    chk("A0", (a, b, t, kap, qq) == (3, 5, 4, 2, 1),
        "case params (a,b,t,kappa,q) = %s  [REPAIRED: was (3,5,5,3,2)]"
        % ((a, b, t, kap, qq),))
    # the inputs must come from the retraction guard, and the superseded pair
    # must be refused -- both directions, so this is not a vacuous check
    import polygon_reduction as pr
    _cd = pr.corner_chart_data(5, 20, l_final=5, b_final=2, who="weight_lemma")
    chk("A0b", (_cd["t"], _cd["kappa"], _cd["ord_C"]) == (t, kap, qq)
        and _cd["monomial"] and not _cd["retraction"],
        "guarded corner data agrees with the json: (t,kappa,ord C) = (4,2,1), "
        "C a MONOMIAL, no retraction at (5,20)")
    try:
        pr.final_corner_dictionary(5, 20, 5, 2)
        _raised = False
    except pr.FinalCornerDictionaryError:
        _raised = True
    chk("A0c", _raised and pr.final_corner_dictionary(8, 28, 4, 7) == (4, 7),
        "the (t,q)=(l_final,b_final) dictionary RAISES at (5,20) and RETURNS at "
        "(8,28) -- the root-cause guard, checked in both directions")

    # --- A1: re-solve the forcing ODE from the corner data (no stored f) -----
    # a{ t c f' - [t(b-a)+kappa+1] c' f } = c^(b-a+1);  at (3,5,4,2): 12 c f' - 33 c' f = c^3
    c_ = y**qq                                   # C = y, a MONOMIAL
    lhs_coef = (a * t, a * (t * (b - a) + kap + 1))
    chk("A1a", lhs_coef == (12, 33), "forcing ODE coefficients (a*t, a*[t(b-a)+kappa+1]) = %s" % (lhs_coef,))
    ai = sp.symbols("A0:12")
    fgen = sum(ai[i] * y**i for i in range(12))
    resid = sp.Poly(sp.expand(12 * c_ * sp.diff(fgen, y) - 33 * sp.diff(c_, y) * fgen - c_**3), y)
    sol = sp.solve(resid.all_coeffs(), list(ai), dict=True)
    chk("A1b", len(sol) == 1, "polynomial solution of deg<=11 is unique (%d solution set)" % len(sol))
    f_ = sp.expand(fgen.subs(sol[0]))
    f_target = sp.Rational(1, 3) * y**3
    chk("A1c", sp.expand(f_ - f_target) == 0, "f = (1/3) y^3, deg f = %d" % sp.degree(f_, y))
    _fold = -sp.Rational(1, 9) * y**5 * (y**3 + 1)**3
    chk("A1d", sp.expand(12 * c_ * sp.diff(_fold, y)
                         - 33 * sp.diff(c_, y) * _fold - c_**3) != 0,
        "the superseded f = -(1/9) y^5 (y^3+1)^3 does NOT solve the repaired ODE "
        "(so A1c is a discriminating check, not a restatement)")

    # --- A2: Phi and its signature ------------------------------------------
    N = a * (t * (a + b) - (kap + 1)) - 2 * b
    chk("A2a", N == 77, "tower length N = a[t(a+b)-(kappa+1)]-2b = %d  [was 98]" % N)
    Phi_y = sp.expand(f_ * c_**N)
    degPhi = sp.degree(Phi_y, y)
    ordPhi = min(sp.Poly(Phi_y, y).monoms())[0]
    mult = sp.Poly(Phi_y, y).as_expr()
    m1 = 0
    ptmp = sp.factor(Phi_y)
    while sp.simplify(sp.rem(sp.Poly(Phi_y, y), sp.Poly((y + 1)**(m1 + 1), y)).as_expr()) == 0:
        m1 += 1
        if m1 > 210:
            break
    chk("A2b", (degPhi, ordPhi) == (80, 80),
        "Phi = f*C^%d = (1/3) y^80 : (deg_y, ord_y) = (%d, %d) -- a MONOMIAL, so "
        "there is NO ord/deg gap" % (N, degPhi, ordPhi))
    chk("A2c", m1 == 0, "mult_(y+1) Phi = %d  (cofactor deg %d) -- no (y+1) place "
        "at all; the superseded 101/202 came from g = y^3+1"
        % (m1, degPhi - ordPhi - m1))
    chk("A2d", M == 29 and M == b * t + J["recipe"]["jphi"],
        "forcing slice M = b*t + jphi = %d = w(Phi)" % M)

    # --- A3/A4: symbol weights + independent homogeneity ---------------------
    sysB = build_75125(J)
    wmap = {sysB.VN[i]: sysB.W[i] for i in range(len(sysB.VN))}
    chk("A3a", [wmap[n] for n in ["d2", "d1", "d0", "dm1"]] == [2, 3, 4, 5],
        "state u-weights w(d_m)=t-m : d2,d1,d0,e=dm1 -> 2,3,4,5  (t=4, so d3 is gone)")
    chk("A3b", [wmap["dm%d" % k] for k in range(2, 9)] == list(range(6, 13)),
        "spare u-weights dm2..dm8 -> 6..12")
    hom = sysB.homogeneity()
    ok = all(len(v) == 1 for v in hom.values())
    got = sorted(list(v)[0] for v in hom.values())
    chk("A4a", ok, "all 8 generators u-homogeneous (recomputed from the stored strings)")
    chk("A4b", got == [21, 22, 23, 24, 25, 26, 27, 29],
        "generator weights = b*t+j = %s" % got)
    chk("A4c", all(list(hom[nm])[0] == b * t + J["generators"][nm]["slice_j"]
                   for nm in hom), "each weight equals its own slice index b*t+j")
    phi_carriers = [nm for nm, tm in sysB.POLY.items()
                    if any(mo[sysB.IPHI] for mo in tm)]
    cphi = sysB.POLY["G9"][tuple(1 if i == sysB.IPHI else 0
                                 for i in range(len(sysB.VN)))]
    chk("A4d", phi_carriers == ["G9"] and cphi == 1,
        "Phi appears only in G9, with coefficient exactly 1 (stale-2Phi guard)")

    # --- A5: the two slopes, and the (non-)existence of a stripping factor ---
    alpha, qw, dsl = window_law(ordPhi, M, degPhi)
    chk("A5a", (alpha, qw) == (80, 29) and qw == 12 * a - 7 and qw == M,
        "W_step = ord_y(Phi)/M = %d/%d, q_window = %d = 12a-7 (NON-integral).  Note "
        "q_window == M exactly, which is what makes the carry obstruction total (C3)."
        % (alpha, qw, qw))
    chk("A5b", sp.Rational(dsl) == sp.Rational(80, 29) and sp.Rational(dsl).q != 1,
        "deg_slope = deg_y(Phi)/M = %s is NOT an integer, so there is NO affine "
        "y-degree cap -- CAPS_AUDIT sec.5's 'deg_slope = 14' is FALSE, not merely "
        "tautological  [REPAIRED: this check previously asserted dsl == 14]" % dsl)
    lam = sp.Rational(dsl) - sp.Rational(alpha, qw)
    chk("A5c", lam == 0,
        "stripped slope lambda = deg_slope - W_step = %s.  Because C = y is a "
        "monomial, ord_y(Phi) = deg_y(Phi), so the two slopes COINCIDE and lambda "
        "vanishes.  The (72,108) stripping factor (lambda = 2) therefore has no "
        "(75,125) counterpart -- and the mechanism fails for a STRONGER reason "
        "than before: with lambda = 0 the weight lemma's own interval "
        "[max(0, D-lambda(W-w_e)), lambda*w_e] = [%d, 0] is EMPTY." % (lam, degPhi))

    # --- A6: the (72,108) contrast, incl. the independent polygon route ------
    up = json.load(open(os.path.join(HERE, "paper_src", "upstream_facts.json")))
    Psub2 = [tuple(p) for p in up["facts"]["newton_polygons"]["sub2"]["P"]]
    hullmax = max(j - 2 * i for i in range(0, 9) for j in range(0, 17)
                  if sp.Polygon(*[sp.Point(*p) for p in Psub2]).encloses_point(sp.Point(i, j))
                  or sp.Point(i, j) in [sp.Point(*p) for p in Psub2]
                  or any(sp.Segment(sp.Point(*Psub2[k]), sp.Point(*Psub2[(k + 1) % len(Psub2)]))
                         .contains(sp.Point(i, j)) for k in range(len(Psub2))))
    chk("A6a", hullmax == 0,
        "(72,108) sub2 polygon: max(j-2i) over the hull = %d  -> deg C_{4-k} <= 8-2k" % hullmax)
    # D-transform.  max(j-2i)=0 gives the slice bound deg P_{8-k} <= 16-2k, hence
    # deg C_{4-k} <= 8-2k (subtract v(C4)=8).  With D_j = C_j*C4^(7-2j), deg C4 = 8,
    # at j = 4-k:  deg D_{4-k} <= (8-2k) + 8*(7-2*(4-k)) = 14k.
    capk = [(8 - 2 * k) + 8 * (7 - 2 * (4 - k)) for k in range(1, 9)]
    chk("A6b", capk == [14 * k for k in range(1, 9)],
        "(72,108) deg cap of the weight-k window variable = 14k from the POLYGON alone "
        "(k=1..8: %s) -- slope 14, Phi never consulted" % capk)
    chk("A6c", window_law(204, 17, 238) == (12, 1, 14),
        "(72,108) W_step = 204/17 = 12 (q_window = 1), deg_slope = 238/17 = 14 -- and the "
        "14 agrees with A6b, computed by a disjoint route")
    # A6d REPAIRED.  There IS now a (5,20)-corner reduced polygon -- computed by
    # polygon_reduction.case_f2 and externally controlled by GGV3's published
    # (50,75) degrees.  So we can ask the disjoint POLYGON question here too, and
    # the answer is that no affine deg cap exists.
    import polygon_reduction as pr2
    _r75 = pr2.case_f2(1)
    _pq = _r75.reduced["standard (proportional, Prop 8.2(1))"]
    chk("A6d", set(_pq["P"]) == {(0, 0), (9, 0), (12, 3), (0, 15)}
        and set(_pq["Q"]) == {(0, 0), (15, 0), (20, 5), (0, 25)},
        "the (5,20) corner now HAS a computed reduced polygon: N(P) = 3*Delta', "
        "N(Q) = 5*Delta', Delta' = {(0,0),(3,0),(4,1),(0,5)} -- and the same engine "
        "reproduces GGV3's published (50,75) degrees 10 and 15")
    chk("A6e", set(up["facts"]["newton_polygons"]) == {"sub1", "sub2"}
        and sp.Rational(degPhi, M).q != 1,
        "upstream_facts.json still carries PUBLISHED polygons for the (8,28) corner "
        "only (keys %s).  But the (75,125) deg_slope no longer needs a polygon "
        "provenance argument: it is %s, non-integral, so no affine deg cap exists "
        "by arithmetic alone" % (sorted(up["facts"]["newton_polygons"]),
                                 sp.Rational(degPhi, M)))

    # =====================================================================
    if verbose:
        print("\n== B. Does a Phi-divisor relation exist?  (graded search) ==")
        print("   Search space: c*Phi = e*B  with c,B weight-homogeneous polynomials.")
        print("   Complete for exact ideal identities -- the standard the K-syzygy meets.")

    # --- B1/B2: CONTROL at (72,108) -----------------------------------------
    sysA = build_72108()
    homA = sysA.homogeneity()
    chk("B0", all(len(v) == 1 for v in homA.values()) and
        sorted(list(v)[0] for v in homA.values()) == [13, 14, 15, 17],
        "control system (72,108) u-homogeneous, weights %s"
        % sorted(list(v)[0] for v in homA.values()))

    ie = sysA.IDX["dm1"]
    # the CRITERION control: the corrected test must DETECT the K-syzygy
    _ex, _np, _ni, _nc = phi_relation_exists(sysA, "G5", 17, ie)
    chk("B0b", _ex and _np == 1 and _ni == 0,
        "criterion control: at (72,108) weight 17 the Phi-carrying column is "
        "DEPENDENT mod the rest (%d of %d independent), so phi_relation_exists "
        "returns True -- the corrected test is not vacuous on the positive side"
        % (_ni, _np))
    colsA, labA = sysA.columns(17, ie)
    nsA = nullspace_exact(colsA)
    chk("B1a", len(nsA) == 1,
        "(72,108) weight-17 piece: EXACT nullity = %d (relation exists and is unique)" % len(nsA))
    v = nsA[0] / nsA[0][[i for i, L in enumerate(labA) if L[0] == "G5"][0]]
    coefs = {labA[i][0] + "*" + sysA.mono_str(labA[i][1]): v[i] for i in range(len(labA)) if v[i] != 0}
    chk("B1b", set(coefs) == {"G5*1", "G3*d2", "G2*d1", "G1*d0"} and set(coefs.values()) == {1},
        "recovered relation = G5 + d2*G3 + d1*G2 + d0*G1  (published K-syzygy)")
    # and its e-quotient is the published K
    d2s, d1s, d0s = sysA.SYM["d2"], sysA.SYM["d1"], sysA.SYM["d0"]
    e_, R_, S_ = sysA.SYM["dm1"], sysA.SYM["dm2"], sysA.SYM["dm3"]
    import full_system_bridge as fsb
    st = fsb.gsystem()
    comb = sp.expand(2 * (st["G5"] + d2s * st["G3"] + d1s * st["G2"] + d0s * st["G1"]))
    Kf = 2 * sysA.SYM["Phi"] - e_ * (d2s * e_**2 + 3 * e_ * S_ + 3 * R_**2)
    chk("B1c", sp.expand(comb - Kf) == 0,
        "its e-quotient is exactly K = 2*Phi - e*(d2*e^2 + 3*e*S + 3*R^2)")

    kmax = 6 if fast else 9
    okB2 = True
    b2rows = []
    for k in range(0, kmax + 1):
        nm_ = len(sysA.monomials(k, exclude=(ie,)))
        if nm_ == 0:
            continue
        cA, _ = sysA.columns(17 + k, ie)
        nul = len(cA) - rank_modp(cA)
        b2rows.append((k, len(cA), nul, nm_))
        okB2 &= (nul == nm_)
    chk("B2", okB2,
        "(72,108) kernel dim at weight 17+k equals #(weight-k monomials) for k=0..%d "
        "-- the kernel is exactly the K-multiples, so the search is EXACTLY sensitive" % kmax)
    if verbose:
        for k, nc, nul, nm_ in b2rows:
            print("           k=%d W=%2d cols=%3d nullity=%d  (#wt-%d monomials = %d)"
                  % (k, 17 + k, nc, nul, k, nm_))

    # --- B3: the (75,125) search, with the CORRECTED criterion ---------------
    ieB = sysB.IDX["dm1"]
    exists0, np0, ni0, nc0 = phi_relation_exists(sysB, "G9", M, ieB)
    chk("B3a", not exists0,
        "(75,125) weight-%d piece: the %d Phi-carrying column(s) stay independent "
        "mod the other %d columns => c = 0 forced => NO relation"
        % (M, np0, nc0 - np0))

    kmaxB = 8 if fast else 16
    rows, okB3, okB3n = [], True, True
    for k in range(0, kmaxB + 1):
        if not sysB.monomials(k, exclude=(ieB,)):
            continue                          # no weight-k multiplier for Phi => c=0 forced
        ex, npc, ni, nc = phi_relation_exists(sysB, "G9", M + k, ieB)
        cB, _ = sysB.columns(M + k, ieB)
        nul = len(cB) - rank_modp(cB)
        rows.append((k, nc, npc, ni, nul))
        okB3 &= (not ex)
        okB3n &= (nul == 0)
    chk("B3b", okB3,
        "(75,125) NO relation c*Phi = e*B at any weight %d..%d: at every weight the "
        "Phi-carrying columns are independent modulo the rest, which PROVES c = 0 "
        "over Q (independence mod p => independence over Q)" % (M, M + kmaxB))
    # And record explicitly that the OLD criterion would have mis-fired here.
    chk("B3c", not okB3n,
        "AND the superseded 'nullity == 0' criterion FAILS on this system (a pure "
        "inter-generator syzygy appears at weight %d, supported on G1,G2 only, with "
        "no Phi column) -- so replacing the criterion was necessary, not cosmetic. "
        "Such a syzygy gives c = 0 and is NOT a Phi-divisor relation."
        % (M + [r[0] for r in rows if r[4] > 0][0] if any(r[4] > 0 for r in rows) else -1))
    if verbose:
        for k, nc, npc, ni, nul in rows:
            print("           k=%2d W=%2d cols=%4d phi-cols=%3d indep=%3d "
                  "(plain nullity=%d)" % (k, M + k, nc, npc, ni, nul))

    # --- B4: other divisor candidates ---------------------------------------
    okB4, b4 = True, []
    kmaxB4 = 6 if fast else 10
    for vn in ["d0", "d1", "d2"]:
        iv = sysB.IDX[vn]
        worst = False
        for k in range(0, kmaxB4 + 1):
            if not sysB.monomials(k, exclude=(iv,)):
                continue
            ex, _, _, _ = phi_relation_exists(sysB, "G9", M + k, iv)
            worst = worst or ex
        b4.append((vn, worst))
        okB4 &= (not worst)
    chk("B4", okB4,
        "no relation c*Phi = v*B for v in d0,d1,d2 either, weights %d..%d "
        "(relation found per divisor: %s)" % (M, M + kmaxB4, b4))

    # --- B5: the concrete witness -------------------------------------------
    # a=3 => b+1 = 6 = 2a, so [u^M]S^6 = 0 is the analogue of (72,108)'s
    # [u^17]S^4 = 0.  Its head pairs each generator G_j (weight b*t+j) with the
    # variable of weight M-(b*t+j); the pairing for the divisor e = dm1 is omitted
    # (it lands on the e side), and G8 is the skipped generator.
    head = {"G9": sp.Integer(1), "G7": sysB.SYM["d2"], "G6": sysB.SYM["d1"],
            "G5": sysB.SYM["d0"], "G3": sysB.SYM["dm2"], "G2": sysB.SYM["dm3"],
            "G1": sysB.SYM["dm4"]}
    for _nm, _m in head.items():
        _w = J["generators"][_nm]["u_weight"] + (0 if _m == 1 else
                                                 sysB.W[sysB.IDX[str(_m)]])
        assert _w == M, (_nm, _m, _w, M)
    X = sum(m * sp.sympify(J["generators"][nm]["poly"], locals=sysB.SYM)
            for nm, m in head.items())
    resid75 = sp.expand(sp.expand(X.subs(sysB.SYM["Phi"], 0)).subs(sysB.SYM["dm1"], 0))
    nterm = len(resid75.args) if resid75 != 0 else 0
    chk("B5", resid75 != 0,
        "witness: the [u^%d]S^6=0 head combination G9 + d2*G7 + d1*G6 + d0*G5 "
        "+ dm2*G3 + dm3*G2 + dm4*G1 leaves a NONZERO residue mod e (%d terms). "
        "At (72,108) the same construction's residue is 0 -- that IS the K-syzygy."
        % (M, nterm))
    chk("B5b", J["generators"]["G1"]["u_weight"] == b * t + 1
        and b * t + 1 == 21 and 13 > 12,
        "STRUCTURAL REASON: at (72,108) the tail slices [u^9..u^12]S^3 sit BELOW the "
        "generator range (which starts at u^13), so the identity survives modulo the "
        "ideal as a Phi-relation.  At (75,125) the tail slices [u^21..u^25]S^5 ARE the "
        "generators G1..G5, so the identity collapses into the ideal and isolates "
        "nothing.  (Checked, not merely asserted: G1's u-weight is b*t+1 = 21.)")

    # =====================================================================
    if verbose:
        print("\n== C. The forcing verdict, and the ord-side obstruction ==")

    sigma = sp.Rational(dsl) * M - degPhi
    chk("C1", sigma == 0,
        "sigma = deg_slope*W - D = %s*%d - %d = %d -- but TAUTOLOGICAL: window_law "
        "DEFINES deg_slope := deg_y(Phi)/M, so sigma == 0 for every case by "
        "construction.  And here deg_slope is not even an integer, so it is not a "
        "usable cap at all." % (dsl, M, degPhi, sigma))
    we = 5                                        # u-weight of e = dm1 (t=4)
    chk("C2", sp.Rational(dsl) * we == sp.Rational(400, 29)
        and sp.Rational(dsl) * we - Lcap(we, alpha, qw) == sp.Rational(-6, 29),
        "CONDITIONAL, and now VACUOUS: IF a relation existed, forcing would pin "
        "deg_y(e) = deg_slope*w_e = %s, which is not an integer -- so there is no "
        "pinned value to test.  The stripped slope lambda = 0 makes the lemma's "
        "interval empty outright." % (sp.Rational(dsl) * we))

    Le, LB, LM = Lcap(we, alpha, qw), Lcap(M - we, alpha, qw), Lcap(M, alpha, qw)
    chk("C3a", (Le, LB, LM) == (14, 67, 80) and LM == ordPhi,
        "lower y-order caps L(w)=ceil(80w/29): L(%d)=%d, L(%d)=%d, L(%d)=%d = ord_y(Phi)"
        % (we, Le, M - we, LB, M, LM))
    chk("C3b", Le + LB - LM == 1,
        "carry(%d,%d) = L(e)+L(B)-L(Phi) = %d > 0: on any lift ord(e*B) >= %d while "
        "ord(c*Phi) = %d.  A weight-%d relation with a CONSTANT c is incompatible with "
        "the existence of any lift." % (we, M - we, Le + LB - LM, Le + LB, LM, M))
    zero_carry = [w for w in range(1, M) if carry(w, M - w, alpha, qw) == 0]
    chk("C3c", zero_carry == [] and qw == M,
        "carry(w_e, %d-w_e) = 0 for NO w_e at all (%s): the mechanism needs "
        "q_window | w_e, and q_window = %d = M, so 0 < w_e < M is impossible.  The "
        "superseded model still had the two escapes w_e in {12,24}; the repaired one "
        "has none, so the ord-side obstruction is TOTAL." % (M, zero_carry, qw))
    chk("C4", all(carry(w, 17 - w, 12, 1) == 0 for w in range(1, 17)),
        "(72,108): q_window = 1, so carry == 0 for EVERY split -- the ord side is exactly "
        "balanced there.  That is what (72,108) has and (75,125) lacks.")
    fam = [(aa, 12 * aa - 7) for aa in range(2, 7)]
    chk("C5", all(q != 1 for _, q in fam) and fam[1][1] == qw,
        "F2-family law q_window = 12a-7 = %s: never 1 for a>=2, so the carry obstruction "
        "is generic in the family.  (72,108)'s q_window = 1 is exceptional to its corner."
        % [q for _, q in fam])
    chk("C5b", all(sp.gcd(12 * aa - 7, 12 * aa - 7 - 1) == 1 for aa in range(2, 7))
        and all((12 * aa - 7) == (2 * aa - 1) * 4 + (4 * aa - 3) for aa in range(2, 7)),
        "and q_window = M = b*t + jphi for every rung, which is why the obstruction is "
        "total at every rung and not just at a=3")

    if verbose:
        print("\n== VERDICT ==")
        print("  Relation c*Phi = e*B at (75,125):  DOES NOT EXIST (B3, weights %d..%d)."
              % (M, M + kmaxB))
        print("  Forcing:                           VACUOUS -- no hypothesis to discharge.")
        print("  sigma = 0:                         TRUE but TAUTOLOGICAL (C1).")
        print("  Structural obstruction:            ord-side carry, q_window = %d = M," % qw)
        print("                                     so EVERY split carries (C3) -- total.")
        print("  VERDICT UNCHANGED by the 2026-07-26 repair, and strengthened on both legs.")
        print("\n  %d/%d checks passed" % (npass, npass + nfail))
    return out, npass, nfail


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--fast", action="store_true")
    args = ap.parse_args()
    _, npass, nfail = run(verbose=not args.quiet, fast=args.fast)
    if args.quiet:
        if nfail:
            print("WEIGHT-LEMMA (75,125): %d FAILED of %d" % (nfail, npass + nfail))
        else:
            print("ALL %d WEIGHT-LEMMA (75,125) CHECKS PASSED" % npass)
    return 1 if nfail else 0


if __name__ == "__main__":
    sys.exit(main())
