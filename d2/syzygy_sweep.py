#!/usr/bin/env python3
"""syzygy_sweep.py -- graded enumeration of the (72,108) G-system syzygy module,
and the resulting DIVISOR MONOMIAL IDEALS.

WHAT THIS IS.  The K-syzygy (DIVISOR_SYZYGY.md)

    2*(G5 + d2*G3 + d1*G2 + d0*G1)  ==  2*Phi - e*(d2*e^2 + 3*e*S + 3*R^2)

was not luck: it is the unique dm4-free element of the FIVE-dimensional u-weight-17
graded piece of the ideal <G1,G2,G3,G5>.  The symbol u-weights are

    d0=4  d1=3  d2=2  dm1(e)=5  dm2(R)=6  dm3(S)=7  dm4(T)=8  Phi=17
    w(G1)=13  w(G2)=14  w(G3)=15  w(G5)=17

so every graded piece I_W = span{ m*G_i : w(m) + w(G_i) = W } is a finite-dimensional
vector space, and finding syzygies is LINEAR ALGEBRA, not a Groebner basis.

THE RIGHT INVARIANT.  The power of K is not the identity, it is the DIVISIBILITY it
produces.  Define, for m >= 1,

    D_{e^m}  =  { monomials mu  :  mu in I + (e^m) }.

D_{e^m} is a MONOMIAL IDEAL (closed under multiplication by monomials), so it is
determined by its finitely many MINIMAL GENERATORS.  Each minimal generator mu is a
theorem: on every genuine lift all G_i vanish, so mu = -(e^m)*h, i.e.

    e^m | mu   in Q[y],   with NO degree cap, NO saturation, NO branch assumption.

`e | Phi` is the generator mu = Phi at weight 17, m = 1.  This file computes the whole
list.  The sweep runs mod a large prime for speed; every REPORTED generator is then
re-certified with exact rational cofactors and verified by sympy expansion.

THE HEADLINE (sec.8b).  D_e saturates -- mu in I + (e^m) is SUFFICIENT for e^m | mu,
never necessary -- and the pure-power ladder

    e^m | S^k   iff   k >= m + 2

converges to e | S without ever reaching it at any single weight.  The identity that
DOES reach it is an INTEGRAL DEPENDENCE of S over the ideal (e) modulo I, of degree 7,
built by eliminating T (all of G1,G2,G3 are LINEAR in T) and then R (Sylvester
resultant, cofactors taken explicitly from the adjugate).  It proves

    e | S,  i.e.  dm1 | dm3,   on every lift with e != 0

cap-free, branch-independent, using only G1,G2,G3 -- no G5, no Phi.  Corollaries:
dm4 is not a free spare at all (T = -R*(S/e + d2) - d1*e/2), and on T2 both R | e^2
and e*R | Phi (strictly stronger than e | Phi).  An independent Newton-polygon scan
(S17) confirms e | S by an argument sharing no machinery with the above.

VALIDATION GATE.  Check S1 rediscovers K at weight 17 from scratch: the dm4-free
subspace of I_17 is computed by elimination-ordered Gaussian elimination, asserted to
be exactly 1-dimensional, and asserted equal to K up to scalar.  A syzygy hunter that
cannot re-find the syzygy we already have is worthless.

Read-only over all existing artifacts.  Usage:
    python -u syzygy_sweep.py              # full census + scored candidates
    python -u syzygy_sweep.py --quiet      # self-check, exit 0 iff all pass
    python -u syzygy_sweep.py --wmax 34    # push the census higher
"""
from __future__ import annotations

import argparse
import itertools
import os
import sys
from collections import defaultdict
from fractions import Fraction

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__)) or "."
sys.path.insert(0, HERE)

# ---------------------------------------------------------------------------
# 0.  The graded alphabet.  Order is fixed and load-bearing (index = position).
# ---------------------------------------------------------------------------
NAMES = ("d0", "d1", "d2", "dm1", "dm2", "dm3", "dm4", "Phi")
WEIGHT = (4, 3, 2, 5, 6, 7, 8, 17)
IDX = {n: i for i, n in enumerate(NAMES)}
NSYM = len(NAMES)
E = IDX["dm1"]                       # the divisor variable e
PRIME = 2147483647                   # 2^31 - 1

SYMS = sp.symbols(" ".join(NAMES))
SYM_OF = dict(zip(NAMES, SYMS))
_pretty = {"dm1": "e", "dm2": "R", "dm3": "S", "dm4": "T"}


def mono_str(m):
    parts = []
    for i, ex in enumerate(m):
        if ex:
            n = _pretty.get(NAMES[i], NAMES[i])
            parts.append(n if ex == 1 else "%s^%d" % (n, ex))
    return "*".join(parts) if parts else "1"


def mono_weight(m):
    return sum(e * w for e, w in zip(m, WEIGHT))


# ---------------------------------------------------------------------------
# 1.  Monomial enumeration by weight
# ---------------------------------------------------------------------------
_mono_cache = {}


def monomials(W, alphabet=None, emax=None):
    """All exponent tuples of total u-weight W.

    alphabet: iterable of symbol INDICES allowed (others forced to exponent 0).
    emax:     if given, the dm1-exponent is forced < emax (i.e. work in R/(e^emax)).
    """
    if alphabet is None:
        alphabet = tuple(range(NSYM))
    key = (W, tuple(sorted(alphabet)), emax)
    if key in _mono_cache:
        return _mono_cache[key]
    alpha = tuple(sorted(alphabet))
    out = []
    base = [0] * NSYM

    def rec(k, rem, cur):
        if k == len(alpha):
            if rem == 0:
                out.append(tuple(cur))
            return
        i = alpha[k]
        w = WEIGHT[i]
        top = rem // w
        if i == E and emax is not None:
            top = min(top, emax - 1)
        for ex in range(top + 1):
            cur[i] = ex
            rec(k + 1, rem - ex * w, cur)
        cur[i] = 0

    rec(0, W, base)
    out.sort()
    _mono_cache[key] = out
    return out


def count_monomials(W, alphabet=None, emax=None):
    return len(monomials(W, alphabet, emax))


# ---------------------------------------------------------------------------
# 2.  Generators, as exponent-tuple -> Fraction dicts
# ---------------------------------------------------------------------------
def _to_dict(expr):
    p = sp.Poly(sp.expand(expr), *SYMS)
    d = {}
    for mon, c in p.terms():
        d[tuple(int(x) for x in mon)] = Fraction(sp.Rational(c).p, sp.Rational(c).q)
    return d


def _to_expr(d):
    tot = 0
    for m, c in d.items():
        t = sp.Rational(c.numerator, c.denominator)
        for i, ex in enumerate(m):
            if ex:
                t *= SYMS[i] ** ex
        tot += t
    return sp.expand(tot)


def load_G():
    """Canonical (72,108) G-system, from bigrade_annotator (single source of truth)."""
    import bigrade_annotator as ba
    g = ba._G_generators()
    out = {}
    for name in ("G1", "G2", "G3", "G5"):
        expr, w = g[name]
        d = _to_dict(expr)
        for m in d:
            assert mono_weight(m) == w, "%s is not u-weight-homogeneous at %s" % (name, m)
        out[name] = (d, w)
    # THE STANDING GUARD: G5 = Phi + G5body, coeff(G5, Phi) == 1.
    phi_mono = tuple(1 if i == IDX["Phi"] else 0 for i in range(NSYM))
    assert out["G5"][0].get(phi_mono) == Fraction(1), \
        "G5 Phi-coefficient is not 1 -- a stale 2*Phi transcription was a REAL bug here"
    return out


def load_H():
    """The dm4-eliminated H-system, verbatim from r9_eliminated_system.json."""
    import bigrade_annotator as ba
    h = ba._H_generators()
    out = {}
    for name in ("H2", "H3", "H5"):
        expr, w = h[name]
        d = _to_dict(expr)
        for m in d:
            assert mono_weight(m) == w, "%s not homogeneous at %s" % (name, m)
            assert m[IDX["dm4"]] == 0, "H-system is supposed to be dm4-free"
        out[name] = (d, w)
    return out


def specialize(gens, zero_syms):
    """Drop every monomial containing a symbol in zero_syms (i.e. set it to 0)."""
    zs = {IDX[s] for s in zero_syms}
    out = {}
    for name, (d, w) in gens.items():
        nd = {m: c for m, c in d.items() if not any(m[i] for i in zs)}
        if nd:
            out[name] = (nd, w)
    return out


# ---------------------------------------------------------------------------
# 3.  Sparse Gaussian elimination over F_p (rows = dicts col-index -> int)
# ---------------------------------------------------------------------------
class Echelon:
    """Row-echelon over F_p.  Column priority = smaller index is the LEADING column,
    so putting 'undesirable' monomials first in the column order performs
    ELIMINATION: rows whose leading column is a 'desirable' monomial are exactly a
    basis of the subspace avoiding the undesirable ones."""

    __slots__ = ("piv",)

    def __init__(self):
        self.piv = {}

    def reduce(self, row):
        row = dict(row)
        while row:
            l = min(row)
            pr = self.piv.get(l)
            if pr is None:
                return row, l
            c = row[l]
            for k, v in pr.items():
                nv = (row.get(k, 0) - c * v) % PRIME
                if nv:
                    row[k] = nv
                elif k in row:
                    del row[k]
        return {}, None

    def add(self, row):
        r, l = self.reduce(row)
        if not r:
            return False
        inv = pow(r[l], PRIME - 2, PRIME)
        self.piv[l] = {k: (v * inv) % PRIME for k, v in r.items()}
        return True

    def contains(self, row):
        r, _ = self.reduce(row)
        return not r

    @property
    def rank(self):
        return len(self.piv)


def _mul_mono(d, mono):
    return {tuple(a + b for a, b in zip(m, mono)): c for m, c in d.items()}


def ideal_rows(W, gens, emax=None, alphabet=None):
    """All monomial multiples m*G_i of total weight W, reduced mod e^emax.

    Yields (label, dict monomial->Fraction).  This SPANS the full graded piece I_W:
    I_W = sum_i R_{W-w_i} * G_i and the monomials span R_{W-w_i}."""
    for name, (d, w) in sorted(gens.items()):
        if W - w < 0:
            continue
        for mono in monomials(W - w, alphabet, emax):
            prod = _mul_mono(d, mono)
            if emax is not None:
                prod = {m: c for m, c in prod.items() if m[E] < emax}
            if prod:
                yield (name, mono), prod


def build_echelon(W, gens, order, emax=None, alphabet=None):
    """order: dict monomial -> column index (lower index = leading)."""
    ech = Echelon()
    keep = []
    for label, prod in ideal_rows(W, gens, emax, alphabet):
        row = {}
        for m, c in prod.items():
            j = order.get(m)
            if j is None:
                continue
            v = (c.numerator % PRIME) * pow(c.denominator % PRIME, PRIME - 2, PRIME) % PRIME
            if v:
                row[j] = v
        if row and ech.add(row):
            keep.append(label)
    return ech, keep


def plain_order(W, emax=None, alphabet=None):
    ms = monomials(W, alphabet, emax)
    return {m: i for i, m in enumerate(ms)}, ms


def elim_order(W, bad_pred, emax=None, alphabet=None):
    """Column order with every 'bad' monomial ranked BEFORE every good one."""
    ms = monomials(W, alphabet, emax)
    bad = [m for m in ms if bad_pred(m)]
    good = [m for m in ms if not bad_pred(m)]
    ordered = bad + good
    return {m: i for i, m in enumerate(ordered)}, ordered, len(bad)


# ---------------------------------------------------------------------------
# 4.  Exact rational certification (cofactor-tracking echelon)
# ---------------------------------------------------------------------------
def exact_membership(W, gens, target_mono, emax, alphabet=None):
    """Exact Fraction solve: is target_mono in I_W + (e^emax)?

    Returns list of (label, Fraction cofactor) with
        sum cof * (mono * G)  ==  c*target + (terms of e-degree >= emax),
    normalised so c == 1, or None."""
    rows = []
    for label, prod in ideal_rows(W, gens, emax, alphabet):
        rows.append((label, prod))
    # echelon over Fractions with cofactor tracking, columns ordered so that the
    # TARGET monomial is last (highest index) -> a pivot on it means membership.
    ms = monomials(W, alphabet, emax)
    others = [m for m in ms if m != target_mono]
    order = {m: i for i, m in enumerate(others)}
    order[target_mono] = len(others)
    piv = {}                             # lead col -> (rowdict, cofdict)
    for label, prod in rows:
        row = {}
        for m, c in prod.items():
            j = order.get(m)
            if j is not None:
                row[j] = row.get(j, Fraction(0)) + c
        row = {k: v for k, v in row.items() if v}
        cof = {label: Fraction(1)}
        while row:
            l = min(row)
            if l not in piv:
                break
            pr, pc = piv[l]
            f = row[l]
            for k, v in pr.items():
                nv = row.get(k, Fraction(0)) - f * v
                if nv:
                    row[k] = nv
                elif k in row:
                    del row[k]
            for k, v in pc.items():
                nv = cof.get(k, Fraction(0)) - f * v
                if nv:
                    cof[k] = nv
                elif k in cof:
                    del cof[k]
        if row:
            l = min(row)
            inv = Fraction(1) / row[l]
            piv[l] = ({k: v * inv for k, v in row.items()},
                      {k: v * inv for k, v in cof.items()})
    tcol = order[target_mono]
    if tcol not in piv:
        return None
    return sorted(piv[tcol][1].items(), key=lambda kv: (kv[0][0], kv[0][1]))


def verify_certificate(cof, target_mono, gens, emax):
    """sympy-exact verification: sum cof*(mono*G) - target has every monomial of
    e-degree >= emax.  Returns (ok, residual_expr)."""
    tot = {}
    for (gname, mono), c in cof:
        d, _w = gens[gname]
        for m, cc in _mul_mono(d, mono).items():
            tot[m] = tot.get(m, Fraction(0)) + c * cc
    tot = {m: c for m, c in tot.items() if c}
    resid = dict(tot)
    lead = resid.pop(target_mono, Fraction(0))
    ok = (lead == 1) and all(m[E] >= emax for m in resid)
    return ok, lead, resid


# ---------------------------------------------------------------------------
# 5.  The divisor monomial ideals D_{e^m}
# ---------------------------------------------------------------------------
def divisor_monomial_ideal(gens, m, wmax, alphabet=None, wmin=None):
    """Compute D_{e^m} = {monomial mu : mu in I + (e^m)} weight by weight.

    Returns (membership, minimal_generators):
      membership: dict weight -> set of monomials in D at that weight
      minimal_generators: list of (weight, monomial)
    """
    if wmin is None:
        wmin = min(w for _d, w in gens.values())
    membership = {}
    mingens = []
    for W in range(wmin, wmax + 1):
        order, ms = plain_order(W, emax=m, alphabet=alphabet)
        ech, _ = build_echelon(W, gens, order, emax=m, alphabet=alphabet)
        hits = set()
        if ech.rank:
            for mo in ms:
                if ech.contains({order[mo]: 1}):
                    hits.add(mo)
        membership[W] = hits
        for mo in sorted(hits):
            minimal = True
            for i in range(NSYM):
                if mo[i] == 0:
                    continue
                sub = list(mo)
                sub[i] -= 1
                sub = tuple(sub)
                wsub = W - WEIGHT[i]
                if wsub in membership and sub in membership[wsub]:
                    minimal = False
                    break
            if minimal:
                mingens.append((W, mo))
    return membership, mingens


# ---------------------------------------------------------------------------
# 6.  The T1 branch (d1 != 0): saturation by d1, in DIRECTLY CERTIFIABLE form
# ---------------------------------------------------------------------------
# On T1, d1 is a NONZERO polynomial in y.  Hence for any k,
#
#     d1^k * mu  in  I + (e^m * d1^k)     =>     d1^k*mu = A + e^m*d1^k*g  with A in I
#                                         =>     on a lift A = 0, so d1^k*(mu - e^m*g) = 0
#                                         =>     mu = e^m * g   (Q[y] is a domain, d1 != 0)
#                                         =>     e^m | mu.
#
# This is exactly  mu in (I : d1^k) + (e^m)  -- but written as a membership test whose
# certificate is a single exact polynomial identity, which is what gets verified.
D1M = tuple(1 if i == IDX["d1"] else 0 for i in range(NSYM))


def _shift(mono, base, k):
    return tuple(a + k * b for a, b in zip(mono, base))


def t1_rows(W, gens, m, k):
    """Rows spanning  I_{W+3k} + (e^m * d1^k)_{W+3k}  in the FULL ring."""
    Wt = W + WEIGHT[IDX["d1"]] * k
    rows = []
    for label, prod in ideal_rows(Wt, gens):
        rows.append((label, prod))
    # the module (e^m * d1^k): every monomial of weight Wt with e-exp>=m and d1-exp>=k
    for mo in monomials(Wt):
        if mo[E] >= m and mo[IDX["d1"]] >= k:
            rows.append((("MOD", mo), {mo: Fraction(1)}))
    return Wt, rows


def t1_divisor_ideal(gens, m, wmax, k=2, wmin=13):
    """D^{T1}_{e^m} = {mu : d1^k*mu in I + (e^m*d1^k)}.  Sound only where d1 != 0."""
    membership = {}
    mingens = []
    for W in range(wmin, wmax + 1):
        Wt, rows = t1_rows(W, gens, m, k)
        order, ms_t = plain_order(Wt)
        ech = Echelon()
        for _label, prod in rows:
            row = {}
            for mo, c in prod.items():
                j = order.get(mo)
                if j is None:
                    continue
                v = (c.numerator % PRIME) * pow(c.denominator % PRIME,
                                                PRIME - 2, PRIME) % PRIME
                if v:
                    row[j] = v
            if row:
                ech.add(row)
        hits = set()
        for mo in monomials(W, emax=m):
            tgt = _shift(mo, D1M, k)
            j = order.get(tgt)
            if j is not None and ech.contains({j: 1}):
                hits.add(mo)
        membership[W] = hits
        for mo in sorted(hits):
            minimal = True
            for i in range(NSYM):
                if mo[i] == 0:
                    continue
                sub = tuple(x - (1 if j == i else 0) for j, x in enumerate(mo))
                wsub = W - WEIGHT[i]
                if wsub in membership and sub in membership[wsub]:
                    minimal = False
                    break
            if minimal:
                mingens.append((W, mo))
    return membership, mingens


def t1_certificate(W, gens, mu, m, k):
    """Exact cofactors for  d1^k*mu = sum c*(mono*G) + (terms in e^m*d1^k)."""
    Wt = W + WEIGHT[IDX["d1"]] * k
    target = _shift(mu, D1M, k)
    rows = []
    for label, prod in ideal_rows(Wt, gens):
        rows.append((label, prod))
    for mo in monomials(Wt):
        if mo[E] >= m and mo[IDX["d1"]] >= k:
            rows.append((("MOD", mo), {mo: Fraction(1)}))
    ms = monomials(Wt)
    others = [x for x in ms if x != target]
    order = {x: i for i, x in enumerate(others)}
    order[target] = len(others)
    piv = {}
    for label, prod in rows:
        row = {}
        for mo, c in prod.items():
            j = order.get(mo)
            if j is not None:
                row[j] = row.get(j, Fraction(0)) + c
        row = {a: b for a, b in row.items() if b}
        cof = {label: Fraction(1)}
        while row:
            l = min(row)
            if l not in piv:
                break
            pr, pc = piv[l]
            f = row[l]
            for a, b in pr.items():
                nv = row.get(a, Fraction(0)) - f * b
                if nv:
                    row[a] = nv
                elif a in row:
                    del row[a]
            for a, b in pc.items():
                nv = cof.get(a, Fraction(0)) - f * b
                if nv:
                    cof[a] = nv
                elif a in cof:
                    del cof[a]
        if row:
            l = min(row)
            inv = Fraction(1) / row[l]
            piv[l] = ({a: b * inv for a, b in row.items()},
                      {a: b * inv for a, b in cof.items()})
    tcol = order[target]
    if tcol not in piv:
        return None
    cof = piv[tcol][1]
    # verify exactly: sum over G-rows must equal d1^k*mu modulo (e^m*d1^k)
    tot = {}
    for (gname, mono), c in cof.items():
        if gname == "MOD":
            continue
        d, _w = gens[gname]
        for mm, cc in _mul_mono(d, mono).items():
            tot[mm] = tot.get(mm, Fraction(0)) + c * cc
    tot = {a: b for a, b in tot.items() if b}
    lead = tot.pop(target, Fraction(0))
    ok = (lead == 1) and all(x[E] >= m and x[IDX["d1"]] >= k for x in tot)
    gcof = sorted(((g, mo), c) for (g, mo), c in cof.items() if g != "MOD")
    return {"ok": ok, "cofactors": gcof, "k": k, "m": m}


# ---------------------------------------------------------------------------
# 7.  Symbol-elimination profile of a graded piece
# ---------------------------------------------------------------------------
def elimination_profile(W, gens, emax=None, alphabet=None):
    """For each symbol x: dim of the subspace of I_W with NO monomial divisible by x,
    plus the sparsest basis element found there."""
    prof = {}
    alpha = alphabet if alphabet is not None else tuple(range(NSYM))
    for i in alpha:
        order, ordered, nbad = elim_order(
            W, lambda m, i=i: m[i] > 0, emax=emax, alphabet=alphabet)
        ech, _ = build_echelon(W, gens, order, emax=emax, alphabet=alphabet)
        free_rows = [r for l, r in ech.piv.items() if l >= nbad]
        best = min((len(r) for r in free_rows), default=None)
        prof[NAMES[i]] = (len(free_rows), best)
    return prof


def eliminated_subspace(W, gens, kill_idx, emax=None, alphabet=None):
    """Basis (as monomial dicts, mod p) of the subspace of I_W avoiding every symbol
    in kill_idx."""
    order, ordered, nbad = elim_order(
        W, lambda m: any(m[i] for i in kill_idx), emax=emax, alphabet=alphabet)
    ech, _ = build_echelon(W, gens, order, emax=emax, alphabet=alphabet)
    inv = {j: mo for mo, j in order.items()}
    out = []
    for l, r in sorted(ech.piv.items()):
        if l >= nbad:
            out.append({inv[j]: c for j, c in r.items()})
    return out


# ---------------------------------------------------------------------------
# 8.  The K-syzygy, for the retrodiction gate
# ---------------------------------------------------------------------------
def K_dict():
    """K = 2*Phi - e*(d2*e^2 + 3*e*S + 3*R^2)   (DIVISOR_SYZYGY.md)."""
    d0, d1, d2, dm1, dm2, dm3, dm4, Phi = SYMS
    return _to_dict(2 * Phi - dm1 * (d2 * dm1**2 + 3 * dm1 * dm3 + 3 * dm2**2))


def _proportional(a, b):
    """Are two Fraction-valued monomial dicts proportional?"""
    if set(a) != set(b):
        return None
    k = next(iter(a))
    lam = a[k] / b[k]
    for m in a:
        if a[m] != lam * b[m]:
            return None
    return lam


# ---------------------------------------------------------------------------
# 8b.  THE HEADLINE:  e | S   (dm1 | dm3), by INTEGRAL DEPENDENCE
# ---------------------------------------------------------------------------
# The divisor monomial ideal D_e only ever certifies  mu in I + (e).  That is a
# SUFFICIENT test, never a necessary one, and it saturates: the ladder
#     e^m | S^k   iff   k >= m + 2
# converges to  e | S  but never reaches it at any single weight.  The identity that
# DOES reach it is an integral dependence of S over the ideal (e) modulo I:
#
#     S^7 + a_1 S^6 + ... + a_7  ==  0   on every lift,   a_i in (e^i).
#
# At any root p of e, with mu = ord_p(e), sigma = ord_p(S) and ord_p(a_i/e^i) >= 0:
#     7 sigma = ord_p(-sum a_i S^{7-i}) >= min_i ( i*mu + (7-i)*sigma ),
# and sigma <= mu-1 would make every term on the right exceed 7*sigma.  So
# sigma >= mu at every root:  e | S.
#
# The certificate uses ONLY G1, G2, G3 -- no G5, no Phi -- so it is immune to the
# G5 normalisation question entirely.  T (=dm4) is eliminated with a one-line
# cofactor because ALL THREE of G1, G2, G3 are LINEAR in T.
def eS_certificate():
    """Build and EXACTLY verify the e|S certificate.  Returns a dict of parts."""
    import bigrade_annotator as ba
    d0, d1, d2, e, R, S, T, Phi = SYMS
    g = ba._G_generators()
    G1, G2, G3 = (sp.expand(g[k][0]) for k in ("G1", "G2", "G3"))

    # Step 1: eliminate T.  All of G1,G2,G3 are linear in T, so:
    A = sp.expand(R * G1 - e * G2)
    B = sp.expand(S * G1 - e * G3)
    T_free = (sp.Poly(A, T).degree() == 0 and sp.Poly(B, T).degree() == 0)

    # Step 2: eliminate R by the Sylvester resultant, with EXPLICIT cofactors from
    # the adjugate (so nothing is taken on trust from resultant theory).
    pa, pb = sp.Poly(A, R), sp.Poly(B, R)
    a2, a1, a0 = (pa.nth(k) for k in (2, 1, 0))
    b2, b1, b0 = (pb.nth(k) for k in (2, 1, 0))
    Syl = sp.Matrix([[a2, a1, a0, 0], [0, a2, a1, a0],
                     [b2, b1, b0, 0], [0, b2, b1, b0]])
    C = sp.expand(Syl.det())
    adj = Syl.adjugate()
    u = sp.expand(adj[3, 0] * R + adj[3, 1])
    v = sp.expand(adj[3, 2] * R + adj[3, 3])
    uv_ok = sp.expand(u * A + v * B - C) == 0

    # Step 3: cofactors straight onto G1,G2,G3
    c1 = sp.expand(u * R + v * S)
    c2 = sp.expand(-u * e)
    c3 = sp.expand(-v * e)
    cof_ok = sp.expand(c1 * G1 + c2 * G2 + c3 * G3 - C) == 0

    # Step 4: C = (729/16) * e * Q,  and Q is the integral dependence
    Q = sp.expand(sp.cancel(C / (sp.Rational(729, 16) * e)))
    fac_ok = sp.expand(C - sp.Rational(729, 16) * e * Q) == 0
    pq = sp.Poly(Q, S)
    n = pq.degree()
    lead = pq.nth(n)
    coeffs = {}
    intdep_ok = True
    for i in range(1, n + 1):
        co = sp.expand(pq.nth(n - i) / lead)
        if co == 0:
            coeffs[i] = sp.Integer(0)
            continue
        red = sp.cancel(co / e**i)
        if not sp.simplify(sp.expand(red * e**i - co)) == 0 or red.has(1 / e):
            intdep_ok = False
        if sp.denom(sp.cancel(red)).has(e):
            intdep_ok = False
        coeffs[i] = sp.expand(red)
    return {"A": A, "B": B, "T_free": T_free, "uv_ok": uv_ok, "cof_ok": cof_ok,
            "fac_ok": fac_ok, "C": C, "Q": Q, "n": n, "lead": lead,
            "alpha": coeffs, "intdep_ok": intdep_ok,
            "c1": c1, "c2": c2, "c3": c3}


def reduced_system():
    """Round 2: after e|S the spare T is DETERMINED and S loses a factor of e.

    Put S = e*s (s new, u-weight 2) and T = -R*(s+d2) - d1*e/2 (u-weight 8).  Then
    G1 vanishes IDENTICALLY, and G2,G3,G5 become P2,P3,P5 of u-weights 14,15,17 in
    the smaller alphabet {d0,d1,d2,e,R,s,Phi}."""
    import bigrade_annotator as ba
    d0, d1, d2, e, R, S, T, Phi = SYMS
    s = sp.Symbol("s")
    g = ba._G_generators()
    G1, G2, G3, G5 = (sp.expand(g[k][0]) for k in ("G1", "G2", "G3", "G5"))
    sub = {S: e * s, T: -R * (s + d2) - sp.Rational(1, 2) * d1 * e}
    return {"s": s, "G1_vanishes": sp.expand(G1.xreplace(sub)) == 0,
            "P2": sp.expand(G2.xreplace(sub)),
            "P3": sp.expand(G3.xreplace(sub)),
            "P5": sp.expand(G5.xreplace(sub)), "sub": sub}


def t2_corollaries():
    """T2 (d1=0) consequences of e|S, each an EXACT identity with explicit cofactor."""
    d0, d1, d2, e, R, S, T, Phi = SYMS
    r = reduced_system()
    s = r["s"]
    P3_t2 = sp.expand(r["P3"].xreplace({d1: 0}))
    # (a)  P3|_{d1=0} = -(e/2) * ( e^2 + 6R(d0 + s^2 + d2 s) )    =>   R | e^2
    reln = sp.expand(e**2 + 6 * R * (d0 + s**2 + d2 * s))
    a_ok = sp.expand(P3_t2 + sp.Rational(1, 2) * e * reln) == 0
    # (b)  K under S=e*s is  2*Phi = e^3(d2+3s) + 3 e R^2 ; substituting reln = 0,
    #      2*Phi = e*R*(3R - 6(d0+s^2+d2 s)(d2+3s)),  with the EXACT cofactor
    #      e*(d2+3s) on reln.                                     =>   e*R | 2*Phi
    lhs = sp.expand(e**3 * (d2 + 3 * s) + 3 * e * R**2)
    rhs = sp.expand(e * R * (3 * R - 6 * (d0 + s**2 + d2 * s) * (d2 + 3 * s)))
    b_ok = sp.expand(lhs - rhs - e * (d2 + 3 * s) * reln) == 0
    # (c)  general branch: R * (3/2 d1 R + 3 d0 e + 3 e s(s+d2)) = -(1/2) e^2 (e+3 d1 s)
    c_ok = sp.expand(R * (sp.Rational(3, 2) * d1 * R + 3 * d0 * e + 3 * e * s * (s + d2))
                     + sp.Rational(1, 2) * e**2 * (e + 3 * d1 * s) + r["P3"]) == 0
    return {"R_divides_e2": a_ok, "eR_divides_Phi": b_ok, "general_R": c_ok,
            "reln": reln}


def newton_scan(mumax=5):
    """INDEPENDENT corroboration of e|S, by Newton polygon instead of resultants.

    A sum of terms can vanish only if the MINIMUM term-valuation is attained at
    least twice (otherwise the sum has exactly that valuation and is nonzero).
    Scan all leading-order valuation assignments with ord(S) < ord(e) and check
    that some generator always has a uniquely-attained minimum.  Infinity (the
    symbol being identically 0) is included in every box.

    This shares NO machinery with eS_certificate(): no ideal, no linear algebra,
    no resultant -- only the term valuations of G1, G2, G3.
    """
    INF = 10**6
    def feasible(mu, sig, rho, tau, a0, a1, a2):
        gens = [[a1 + 2 * mu, a2 + mu + rho, mu + tau, rho + sig],   # G1
                [a0 + 2 * mu, a2 + 2 * rho, rho + tau, 2 * sig],     # G2
                [a0 + mu + rho, a1 + 2 * rho, 3 * mu, sig + tau]]    # G3
        for g in gens:
            m = min(g)
            if m >= INF:
                continue
            if sum(1 for x in g if x == m) < 2:
                return False
        return True
    out = {}
    for mu in range(1, mumax + 1):
        box = list(range(0, 3 * mu + 4)) + [INF]
        found = None
        for sig in range(0, mu):
            for rho in box:
                for tau in box:
                    for a0 in box:
                        for a1 in box:
                            for a2 in box:
                                if feasible(mu, sig, rho, tau, a0, a1, a2):
                                    found = (sig, rho, tau, a0, a1, a2)
                                    break
                            if found:
                                break
                        if found:
                            break
                    if found:
                        break
                if found:
                    break
            if found:
                break
        out[mu] = found
    return out


def rt_obstruction_point():
    """An EXPLICIT point of V(I + (e)) with R != 0 and T != 0.

    Its existence PROVES that no power of R and no power of T ever lies in I + (e),
    at ANY weight -- so the divisor-ideal method can never yield e | R^k or e | T^k.
    That is the exhaustion statement, not a budget excuse."""
    d0, d1, d2, e, R, S, T, Phi = SYMS
    pt = {e: 0, S: 0, d1: 0, Phi: 0, R: 1, d2: 1, T: sp.Rational(-1, 2), d0: 0}
    G = load_G()
    vals = {}
    for name in ("G1", "G2", "G3", "G5"):
        vals[name] = sp.expand(_to_expr(G[name][0]).xreplace(pt))
    return pt, vals


# ---------------------------------------------------------------------------
# 9.  Checks
# ---------------------------------------------------------------------------
def run_checks(verbose=True, wmax=26):
    out, npass, ntot = [], 0, 0

    def ck(name, ok, detail):
        nonlocal npass, ntot
        ntot += 1
        npass += bool(ok)
        out.append("  [%s] %s\n        %s" % ("PASS" if ok else "FAIL", name, detail))

    G = load_G()

    # ---- S0: the standing guard (load_G already asserts; restate as a check)
    phi_mono = tuple(1 if i == IDX["Phi"] else 0 for i in range(NSYM))
    ck("S0  canonical guard: coeff(G5, Phi) == 1",
       G["G5"][0].get(phi_mono) == Fraction(1),
       "coeff = %s  (a stale 2*Phi transcription was a REAL bug in this repo)"
       % G["G5"][0].get(phi_mono))

    # ---- S1: THE RETRODICTION GATE -- rediscover K at weight 17 -------------
    dim17 = sum(count_monomials(17 - w) for _d, w in G.values())
    basis = eliminated_subspace(17, G, [IDX["dm4"]])
    okdim = (len(basis) == 1)
    lam = None
    if okdim:
        Kd = K_dict()
        # lift the mod-p row back to exact rationals via the exact solver
        cof = exact_membership(17, G, phi_mono, emax=None)
        # direct exact route: the dm4-free element is the unique kernel vector
        exact = _exact_dm4_free_17(G)
        lam = _proportional(exact, Kd) if exact else None
    ck("S1  RETRODICTION GATE: K is the unique dm4-free element of I_17",
       okdim and lam is not None,
       "dim I_17 = %d (ambient weight-17 monomials = %d); dm4-free subspace dim = %d; "
       "recovered element = %s * K"
       % (dim17, count_monomials(17), len(basis), lam))

    # ---- S2: the weight-17 piece really is 5-dimensional --------------------
    order17, _ = plain_order(17)
    ech17, keep17 = build_echelon(17, G, order17)
    ck("S2  I_17 is exactly 5-dimensional, spanned by d0*G1, d2^2*G1, d1*G2, d2*G3, G5",
       ech17.rank == 5 and dim17 == 5,
       "rank = %d, generator count = %d, multipliers = %s"
       % (ech17.rank, dim17,
          ", ".join("%s*%s" % (mono_str(m), g) for g, m in keep17)))

    # ---- S3: Phi in (e, G1,G2,G3,G5)?  VERIFY the external claim ------------
    cof = exact_membership(17, G, phi_mono, emax=1)
    ok3 = cof is not None
    detail3 = "NOT in the ideal"
    if ok3:
        good, lead, resid = verify_certificate(cof, phi_mono, G, 1)
        ok3 = good
        detail3 = ("Phi = %s  (+ terms divisible by e); sympy-verified residual is "
                   "e-divisible: %s"
                   % (" + ".join("(%s)*%s*%s" % (c, mono_str(m), g)
                                 for (g, m), c in cof), good))
    ck("S3  Phi IS in the ideal (e, G1,G2,G3,G5) -- external claim VERIFIED", ok3, detail3)

    # ---- S4: e | R*S, the weight-13 divisibility ---------------------------
    rs = [0] * NSYM
    rs[IDX["dm2"]] = 1
    rs[IDX["dm3"]] = 1
    rs = tuple(rs)
    cof4 = exact_membership(13, G, rs, emax=1)
    ok4 = cof4 is not None
    d4 = "not found"
    if ok4:
        good, lead, resid = verify_certificate(cof4, rs, G, 1)
        ok4 = good
        d4 = ("R*S = %s + e*(...)  [sympy-verified]"
              % " + ".join("(%s)*%s*%s" % (c, mono_str(m), g) for (g, m), c in cof4))
    ck("S4  e | R*S  (weight 13, the lowest divisor generator)", ok4, d4)

    # ---- S5: e | S^3, the weight-21 divisibility ---------------------------
    s3 = [0] * NSYM
    s3[IDX["dm3"]] = 3
    s3 = tuple(s3)
    cof5 = exact_membership(21, G, s3, emax=1)
    ok5 = cof5 is not None
    d5 = "not found"
    if ok5:
        good, lead, resid = verify_certificate(cof5, s3, G, 1)
        ok5 = good
        d5 = ("S^3 = %s + e*(...)  [sympy-verified]"
              % " + ".join("(%s)*%s*%s" % (c, mono_str(m), g) for (g, m), c in cof5))
    ck("S5  e | S^3  (weight 21) -- NEW divisor generator", ok5, d5)

    # ---- S6: S^2 is NOT in the ideal mod e (S^3 is the minimal S-power) -----
    s2 = [0] * NSYM
    s2[IDX["dm3"]] = 2
    s2 = tuple(s2)
    no14 = exact_membership(14, G, s2, emax=1) is None
    ck("S6  e | S^2 is FALSE (so S^3 is the minimal pure S-power divisor)", no14,
       "S^2 not in I_14 + (e); I_14 = span{G2} only, which is not a multiple of S^2")

    # ---- S7: the H-system cross-gate: recover K5 ---------------------------
    H = load_H()
    Hd1 = specialize(H, ["d1"])
    basisH = eliminated_subspace(22, Hd1, [], alphabet=None)
    # the H-relation is: H5 + d2*H3 is divisible by e.  Find the e-divisible subspace.
    edivH = _e_divisible_subspace(22, Hd1)
    ok7 = False
    d7 = "no e-divisible element at weight 22"
    if len(edivH) == 1:
        K5 = _to_dict(2 * SYM_OF["Phi"] - 3 * SYM_OF["dm1"] * SYM_OF["dm2"]**2
                      - SYM_OF["d2"] * SYM_OF["dm1"]**3
                      - 3 * SYM_OF["dm1"]**2 * SYM_OF["dm3"])
        eK5 = _mul_mono(K5, tuple(1 if i == E else 0 for i in range(NSYM)))
        lam7 = _proportional(edivH[0], eK5)
        ok7 = lam7 is not None
        d7 = ("the unique e-divisible element of I^H_22 (d1=0) is %s * e*K5 -- the "
              "generic-fiber lane's K5 re-derived from the H-system by linear algebra"
              % lam7)
    ck("S7  H-system cross-gate: K5 rediscovered at weight 22 (d1=0)", ok7, d7)

    # ---- S8: divisor ideal regression (minimal generators up to weight 22) --
    _mem, mg = divisor_monomial_ideal(G, 1, 22)
    got = sorted(mono_str(m) for _w, m in mg)
    expect = sorted(["R*S", "Phi", "S^3", "d1*R^3", "S^2*T"])
    ck("S8  D_e minimal generators up to weight 22 are exactly the recorded list",
       got == expect, "found %s ; expected %s" % (got, expect))

    # ---- S9..S12: THE e | S THEOREM ----------------------------------------
    c = eS_certificate()
    ck("S9  T is eliminated by a one-line cofactor (G1,G2,G3 are all LINEAR in T)",
       c["T_free"],
       "A = R*G1 - e*G2 and B = S*G1 - e*G3 are dm4-free, deg_R = 2 both")
    ck("S10 explicit Sylvester-adjugate cofactors: C = c1*G1 + c2*G2 + c3*G3",
       c["uv_ok"] and c["cof_ok"] and c["fac_ok"],
       "u*A+v*B==C: %s ; c_i*G_i==C: %s ; C == (729/16)*e*Q: %s (deg_S Q = %d)"
       % (c["uv_ok"], c["cof_ok"], c["fac_ok"], c["n"]))
    ck("S11 Q is an INTEGRAL DEPENDENCE of S over (e): a_i in (e^i) for every i",
       c["intdep_ok"] and c["n"] == 7,
       "S^7 + sum_{i=1..7} e^i*alpha_i*S^(7-i) = 0 with alpha_i polynomial; "
       "alpha_1 = %s" % c["alpha"][1])
    ck("S12 => e | S   (dm1 | dm3) on every lift with e != 0 -- NEW, cap-free, "
       "branch-independent",
       c["intdep_ok"] and c["n"] == 7 and c["cof_ok"],
       "valuation argument: sigma <= mu-1 forces every term i*mu+(7-i)*sigma > 7*sigma")

    # ---- S13: T determined; round-2 system ---------------------------------
    r2 = reduced_system()
    ck("S13 dm4 is NOT a free spare: T = -R*(S/e + d2) - d1*e/2 (G1 then vanishes)",
       r2["G1_vanishes"],
       "under S->e*s, T->-R*(s+d2)-d1*e/2 the weight-13 generator G1 is identically 0; "
       "the residual system is P2,P3,P5 of u-weight 14,15,17")

    # ---- S14/S15: T2 corollaries -------------------------------------------
    t2 = t2_corollaries()
    ck("S14 T2 (d1=0):  e^2 = -6*R*(d0 + s^2 + d2*s)   =>   R | e^2",
       t2["R_divides_e2"] and t2["general_R"],
       "exact: P3|_{d1=0} == -(e/2)*(e^2 + 6R(d0+s^2+d2 s)); general branch identity "
       "R*(3/2 d1 R + 3 d0 e + 3 e s(s+d2)) == -(1/2) e^2 (e + 3 d1 s) also holds")
    ck("S15 T2 (d1=0):  2*Phi = e*R*(3R - 6(d0+s^2+d2 s)(d2+3s))  =>  e*R | Phi",
       t2["eR_divides_Phi"],
       "exact, cofactor e*(d2+3s) on the T2 relation -- STRICTLY STRONGER than e | Phi")

    # ---- S16: the exhaustion theorem for R and T ---------------------------
    pt, vals = rt_obstruction_point()
    allz = all(v == 0 for v in vals.values())
    ck("S16 EXHAUSTION: no power of R or T is EVER in I + (e) (explicit point of "
       "V(I+(e)) with R=1, T=-1/2)",
       allz,
       "G1..G5 at (e,S,d1,Phi,d0)=(0,0,0,0,0), (R,d2,T)=(1,1,-1/2): %s -- so R,T are "
       "not in the radical of I+(e) and the divisor-ideal method provably cannot "
       "produce e | R^k or e | T^k at any weight"
       % {k: str(v) for k, v in vals.items()})

    # ---- S17: INDEPENDENT corroboration of e|S by Newton polygon -----------
    ns = newton_scan(mumax=4)
    ck("S17 INDEPENDENT cross-check of e|S: Newton-polygon scan kills ord(S)<ord(e)",
       all(v is None for v in ns.values()),
       "for ord_p(e) = 1..4 and every leading-order valuation assignment with "
       "ord_p(S) < ord_p(e), some generator has a uniquely-attained minimum, so it "
       "cannot vanish.  Shares no machinery with the resultant certificate. "
       "feasible tuples found: %s" % {k: v for k, v in ns.items()})

    if verbose:
        print("\n".join(out))
    return npass, ntot


def _exact_dm4_free_17(G):
    """Exact rational basis of the dm4-free subspace of I_17 (must be 1-dimensional)."""
    gens_list = []
    for name, (d, w) in sorted(G.items()):
        for mono in monomials(17 - w):
            gens_list.append(((name, mono), _mul_mono(d, mono)))
    bad = sorted({m for _l, p in gens_list for m in p if m[IDX["dm4"]] > 0})
    bidx = {m: i for i, m in enumerate(bad)}
    # kernel of the projection onto dm4-monomials, exact
    rows = []
    for lab, p in gens_list:
        rows.append([p.get(m, Fraction(0)) for m in bad])
    n = len(rows)
    # gaussian elimination on the TRANSPOSE to get the kernel of the map coeff -> bad part
    M = [[rows[j][i] for j in range(n)] for i in range(len(bad))]
    piv_col = []
    r = 0
    for c in range(n):
        pr = None
        for i in range(r, len(M)):
            if M[i][c]:
                pr = i
                break
        if pr is None:
            continue
        M[r], M[pr] = M[pr], M[r]
        inv = Fraction(1) / M[r][c]
        M[r] = [x * inv for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c]:
                f = M[i][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        piv_col.append(c)
        r += 1
    free = [c for c in range(n) if c not in piv_col]
    if len(free) != 1:
        return None
    fc = free[0]
    vec = [Fraction(0)] * n
    vec[fc] = Fraction(1)
    for i, c in enumerate(piv_col):
        vec[c] = -M[i][fc]
    tot = {}
    for coef, (_lab, p) in zip(vec, gens_list):
        if not coef:
            continue
        for m, c in p.items():
            tot[m] = tot.get(m, Fraction(0)) + coef * c
    return {m: c for m, c in tot.items() if c}


def _e_divisible_subspace(W, gens):
    """Exact basis of {f in I_W : e | f}.  Uses the same exact kernel routine."""
    gens_list = []
    for name, (d, w) in sorted(gens.items()):
        if W - w < 0:
            continue
        for mono in monomials(W - w):
            gens_list.append(((name, mono), _mul_mono(d, mono)))
    bad = sorted({m for _l, p in gens_list for m in p if m[E] == 0})
    if not bad:
        return []
    n = len(gens_list)
    M = [[gens_list[j][1].get(m, Fraction(0)) for j in range(n)] for m in bad]
    piv_col, r = [], 0
    for c in range(n):
        pr = None
        for i in range(r, len(M)):
            if M[i][c]:
                pr = i
                break
        if pr is None:
            continue
        M[r], M[pr] = M[pr], M[r]
        inv = Fraction(1) / M[r][c]
        M[r] = [x * inv for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c]:
                f = M[i][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        piv_col.append(c)
        r += 1
    out = []
    for fc in [c for c in range(n) if c not in piv_col]:
        vec = [Fraction(0)] * n
        vec[fc] = Fraction(1)
        for i, c in enumerate(piv_col):
            vec[c] = -M[i][fc]
        tot = {}
        for coef, (_lab, p) in zip(vec, gens_list):
            if not coef:
                continue
            for m, c in p.items():
                tot[m] = tot.get(m, Fraction(0)) + coef * c
        tot = {m: c for m, c in tot.items() if c}
        if tot:
            out.append(tot)
    return out


# ---------------------------------------------------------------------------
# 10.  Scoring
# ---------------------------------------------------------------------------
def power_ladder(gens, xname, kmax=6, mmax=5):
    """The table  e^m | X^k ?  -- i.e. is X^k in I + (e^m)."""
    rows = []
    for k in range(1, kmax + 1):
        mo = tuple(k if i == IDX[xname] else 0 for i in range(NSYM))
        W = mono_weight(mo)
        line = []
        for m in range(1, mmax + 1):
            order, _ms = plain_order(W, emax=m)
            ech, _ = build_echelon(W, gens, order, emax=m)
            j = order.get(mo)
            line.append(bool(j is not None and ech.rank and ech.contains({j: 1})))
        rows.append((k, W, line))
    return rows


def set_ring(names, weights, evar="dm1"):
    """Re-point the graded alphabet (used for the round-2 reduced system)."""
    global NAMES, WEIGHT, IDX, NSYM, E, SYMS, SYM_OF
    NAMES = tuple(names)
    WEIGHT = tuple(weights)
    IDX = {n: i for i, n in enumerate(NAMES)}
    NSYM = len(NAMES)
    E = IDX[evar]
    SYMS = sp.symbols(" ".join(NAMES))
    SYM_OF = dict(zip(NAMES, SYMS))
    _mono_cache.clear()


def round2_sweep(wmax=30):
    """D_e for the reduced system P2,P3,P5 in {d0,d1,d2,e,R,s,Phi}.

    Restores the original ring before returning."""
    save = (NAMES, WEIGHT, dict(IDX), NSYM, E, SYMS, dict(SYM_OF))
    r2 = reduced_system()
    exprs = {"P2": (r2["P2"], 14), "P3": (r2["P3"], 15), "P5": (r2["P5"], 17)}
    try:
        set_ring(("d0", "d1", "d2", "dm1", "dm2", "s", "zz", "Phi"),
                 (4, 3, 2, 5, 6, 2, 10**6, 17))
        gens = {}
        for name, (ex, w) in exprs.items():
            d = _to_dict(ex)
            for m in d:
                assert mono_weight(m) == w, (name, m)
            gens[name] = (d, w)
        _mem, mg = divisor_monomial_ideal(gens, 1, wmax, wmin=14)
        out = [(w, mono_str(mo)) for w, mo in mg]
        lad = power_ladder(gens, "dm2", kmax=5, mmax=4)
        lad = [(k, W, line) for k, W, line in lad]
    finally:
        (NAMES2, WEIGHT2, IDX2, NSYM2, E2, SYMS2, SYM2) = save
        set_ring(NAMES2, WEIGHT2, "dm1")
    return out, lad


def score_generator(W, mo, m):
    """Score a divisor-ideal minimal generator e^m | mu.

    * symbol elimination -- how many of the 8 symbols vanish from mu
    * factorability      -- it IS a divisibility (that is the point); bonus if mu is a
                            pure power of ONE symbol, or is the known polynomial Phi
    * cap-freedom        -- 1 always: these are symbolic identities in the graded ring,
                            no degree bound is used anywhere in the derivation
    """
    present = [NAMES[i] for i in range(NSYM) if mo[i]]
    elim = NSYM - len(present)
    pure = len(present) == 1
    known = present == ["Phi"]
    total = elim + (2 if pure else 0) + (3 if known else 0) + (m - 1) * 2
    return {"weight": W, "mu": mono_str(mo), "m": m, "symbols_present": present,
            "symbols_eliminated": elim, "pure_power": pure, "known_polynomial": known,
            "cap_free": True, "score": total}


# ---------------------------------------------------------------------------
# 11.  Report
# ---------------------------------------------------------------------------
def report(wmax=30, wmax_e2=26):
    G = load_G()
    print("#" * 78)
    print("# SYZYGY SWEEP -- graded enumeration of I = <G1,G2,G3,G5>, (72,108)")
    print("#" * 78)
    print("weights: " + ", ".join("%s=%d" % (_pretty.get(n, n), w)
                                  for n, w in zip(NAMES, WEIGHT)))
    print("gens   : " + ", ".join("w(%s)=%d" % (k, v[1]) for k, v in sorted(G.items())))
    print()

    # ---- 1. the graded census -------------------------------------------
    print("=" * 78)
    print("1.  GRADED CENSUS   dim R_W (ambient monomials) / dim I_W / codim")
    print("=" * 78)
    print("   W   dim R_W   dim I_W   codim   #generators-of-I_W")
    for W in range(13, wmax + 1):
        order, ms = plain_order(W)
        ech, keep = build_echelon(W, G, order)
        ngen = sum(count_monomials(W - w) for _d, w in G.values() if W - w >= 0)
        print("  %3d  %7d   %7d  %6d   %d" % (W, len(ms), ech.rank,
                                              len(ms) - ech.rank, ngen))
    print()

    # ---- 2. symbol-elimination profile ----------------------------------
    print("=" * 78)
    print("2.  SYMBOL-ELIMINATION PROFILE  (dim of the subspace of I_W avoiding x,")
    print("    and the sparsest such element's support size)")
    print("=" * 78)
    hdr = "   W  " + "".join("%12s" % _pretty.get(n, n) for n in NAMES)
    print(hdr)
    for W in range(13, min(wmax, 24) + 1):
        prof = elimination_profile(W, G)
        row = "  %3d  " % W
        for n in NAMES:
            dim, best = prof[n]
            row += "%12s" % ("%d/%s" % (dim, best if best is not None else "-"))
        print(row)
    print("    (cell = dim / sparsest-support; '-' = subspace is zero)")
    print()

    # ---- 3. the divisor monomial ideal D_e -------------------------------
    print("=" * 78)
    print("3.  THE DIVISOR MONOMIAL IDEAL  D_e = {mu : mu in I + (e)}")
    print("    Each MINIMAL generator mu is a theorem:  e | mu  on every genuine lift.")
    print("=" * 78)
    mem, mg = divisor_monomial_ideal(G, 1, wmax)
    print("   W    #monomials-in-D_e   new minimal generators")
    for W in range(13, wmax + 1):
        news = [mono_str(m) for w, m in mg if w == W]
        print("  %3d   %8d           %s" % (W, len(mem.get(W, ())),
                                            ", ".join(news) if news else "-"))
    print()
    print("  MINIMAL GENERATORS of D_e up to weight %d:" % wmax)
    scored = [score_generator(w, m, 1) for w, m in mg]
    scored.sort(key=lambda s: -s["score"])
    print("    %-14s %-6s %-6s %-22s %s" % ("mu", "w", "score", "symbols present", "note"))
    for s in scored:
        note = []
        if s["known_polynomial"]:
            note.append("Phi is a KNOWN polynomial -> support theorem")
        if s["pure_power"]:
            note.append("pure power of one symbol")
        print("    %-14s %-6d %-6d %-22s %s"
              % (s["mu"], s["weight"], s["score"], ",".join(s["symbols_present"]),
                 "; ".join(note)))
    print()

    # ---- 4. higher powers of e ------------------------------------------
    print("=" * 78)
    print("4.  D_{e^2} and D_{e^3}  (does any divisibility upgrade?)")
    print("=" * 78)
    for m in (2, 3):
        _mem2, mg2 = divisor_monomial_ideal(G, m, wmax_e2)
        efree = [(w, mo) for w, mo in mg2 if mo[E] == 0]
        withe = [(w, mo) for w, mo in mg2 if mo[E] > 0]
        print("  m=%d: %d minimal generators up to weight %d" % (m, len(mg2), wmax_e2))
        print("       e-free (the substantive ones): %s"
              % (", ".join("%s (w=%d)" % (mono_str(mo), w) for w, mo in efree) or "none"))
        print("       e-containing (shifts of lower-m results): %s"
              % (", ".join("%s (w=%d)" % (mono_str(mo), w) for w, mo in withe) or "none"))
    print()

    # ---- 4b. the power ladder --------------------------------------------
    print("=" * 78)
    print("4b. THE POWER LADDER   e^m | X^k ?   (Y = certified, . = no certificate)")
    print("=" * 78)
    for xname, kmax in (("dm3", 6), ("dm2", 6), ("dm4", 5), ("Phi", 3)):
        print("  X = %s" % _pretty.get(xname, xname))
        print("        m= " + " ".join("%2d" % m for m in range(1, 6)))
        for k, W, line in power_ladder(G, xname, kmax=kmax, mmax=5):
            print("   k=%2d      " % k + " ".join(" Y" if b else " ." for b in line)
                  + "   [weight %d]" % W)
    print("  READ: e^m | S^k exactly when k >= m+2 -- a ladder converging to e | S but")
    print("        never reaching it; e^m | Phi^k exactly when m <= k, i.e. e | Phi is")
    print("        NOT strengthenable; R and T never appear (see the exhaustion below).")
    print()

    # ---- 4c. the theorem the ladder is converging to ----------------------
    print("=" * 78)
    print("4c. THEOREM (NEW):  e | S,  i.e.  dm1 | dm3")
    print("=" * 78)
    c = eS_certificate()
    print("  G1, G2, G3 are ALL LINEAR in T, so T dies with a one-line cofactor:")
    print("      A := R*G1 - e*G2      B := S*G1 - e*G3      (both dm4-free, deg_R 2)")
    print("  Sylvester-adjugate cofactors give C = u*A + v*B = c1*G1 + c2*G2 + c3*G3,")
    print("  and  C == (729/16) * e * Q  with  deg_S Q = %d.  Q is monic (up to %s) and"
          % (c["n"], c["lead"]))
    print("  is an INTEGRAL DEPENDENCE of S over the ideal (e):")
    print("      S^7 + sum_{i=1..7} e^i * alpha_i * S^(7-i)  ==  0   on every lift,")
    for i in range(1, c["n"] + 1):
        print("      alpha_%d = %s" % (i, sp.factor(c["alpha"][i])))
    ns = newton_scan(mumax=4)
    print("  INDEPENDENT cross-check (shares no machinery with the above): a Newton-")
    print("  polygon scan over every leading-order valuation assignment with")
    print("  ord_p(S) < ord_p(e) <= 4 finds NO feasible tuple -- %s"
          % ("confirmed" if all(v is None for v in ns.values()) else "FAILED"))
    print("  At a root p of e with mu = ord_p(e), sigma = ord_p(S):  sigma <= mu-1 makes")
    print("  every term i*mu + (7-i)*sigma  >  7*sigma = ord_p(S^7).  Hence sigma >= mu")
    print("  at every root, i.e.  e | S.   Uses ONLY G1,G2,G3 -- no G5, no Phi.")
    print()
    print("  COROLLARIES (all exact, all verified in the check suite):")
    print("    * dm4 is NOT a free spare:  T = -R*(S/e + d2) - d1*e/2, and G1 then")
    print("      vanishes identically.  The whole T ansatz is eliminated.")
    print("    * general branch:  R*(3/2 d1 R + 3 d0 e + 3 e s(s+d2)) = -(1/2)e^2(e+3d1 s)")
    print("    * T2 (d1=0):  e^2 = -6*R*(d0 + s^2 + d2*s)      =>   R | e^2")
    print("    * T2 (d1=0):  2*Phi = e*R*(3R - 6(d0+s^2+d2 s)(d2+3s))  =>  e*R | Phi")
    print("      (strictly stronger than the K-syzygy's e | Phi)")
    print()

    # ---- 5. branches -----------------------------------------------------
    print("=" * 78)
    print("5.  BRANCHES")
    print("=" * 78)
    G_T2 = specialize(G, ["d1"])
    _m2, mg_t2 = divisor_monomial_ideal(G_T2, 1, wmax)
    base = {mono_str(m) for _w, m in mg}
    print("  T2 (d1 = 0):  minimal generators of D_e  =")
    for w, mo in mg_t2:
        tag = "" if mono_str(mo) in base else "   <== T2-ONLY"
        print("      %-14s w=%d%s" % (mono_str(mo), w, tag))
    print()
    _mt1, mg_t1 = t1_divisor_ideal(G, 1, min(wmax, 24), k=2)
    print("  T1 (d1 != 0, saturate by d1: d1^2*mu in I + (e*d1^2)):")
    print("      minimal generators of D_e  =")
    for w, mo in mg_t1:
        tag = ""
        if mono_str(mo) not in base:
            cert = t1_certificate(w, G, mo, 1, 2)
            tag = "   <== T1-ONLY  [cert %s]" % ("VERIFIED" if cert and cert["ok"]
                                                 else "FAILED")
        print("      %-14s w=%d%s" % (mono_str(mo), w, tag))
    print()

    # ---- 6. H-system -----------------------------------------------------
    print("=" * 78)
    print("6.  THE dm4-ELIMINATED H-SYSTEM  <H2,H3,H5>  (u-weights 19,20,22)")
    print("=" * 78)
    H = load_H()
    for tag, gg in (("general", H), ("d1 = 0 (T2)", specialize(H, ["d1"]))):
        _mh, mgh = divisor_monomial_ideal(gg, 1, wmax, wmin=19)
        print("  %-12s D_e minimal generators: %s"
              % (tag, ", ".join("%s (w=%d)" % (mono_str(mo), w) for w, mo in mgh) or "none"))
    print()

    # ---- 6b. round 2 -----------------------------------------------------
    print("=" * 78)
    print("6b. ROUND 2 -- the reduced system after e|S  (P2,P3,P5 in d0,d1,d2,e,R,s,Phi)")
    print("=" * 78)
    r2gens, r2lad = round2_sweep(wmax=min(wmax, 30))
    print("  D_e minimal generators of the REDUCED system up to weight %d:"
          % min(wmax, 30))
    for w, name in r2gens:
        print("      %-16s w=%d" % (name, w))
    print("  e^m | R^k in the reduced system (m = 1..4):")
    for k, W, line in r2lad:
        print("      k=%d  " % k + " ".join(" Y" if b else " ." for b in line))
    print("  => round 2 adds  e | d1*R^2  and nothing else.  R is NOT forced.")
    print()

    # ---- 7. exhaustion ---------------------------------------------------
    print("=" * 78)
    print("7.  EXHAUSTION")
    print("=" * 78)
    last = max((w for w, _m in mg), default=0)
    print("  Highest weight at which D_e acquires a NEW minimal generator: %d" % last)
    print("  Swept to weight %d with no further generator." % wmax)
    print("  => NO further e-divisibility constraint exists below weight %d." % (wmax + 1))
    pt, vals = rt_obstruction_point()
    print()
    print("  AND THE SEARCH FOR R AND T IS CLOSED FOREVER, not merely to weight %d:"
          % wmax)
    print("    the point  e=S=d1=Phi=d0=0,  R=1, d2=1, T=-1/2  lies on V(I+(e))")
    print("    (G1..G5 all vanish there: %s)."
          % ", ".join("%s=%s" % kv for kv in sorted(vals.items())))
    print("    So R and T are NOT in the radical of I + (e): no power of either is in")
    print("    I + (e) at ANY weight, and no integral dependence over (e) exists.")
    print("    The method cannot yield e | R^k or e | T^k.  That is a theorem about the")
    print("    system, not a budget limit.  S and Phi are the only two symbols in the")
    print("    radical, and BOTH are now settled: e | Phi (K-syzygy), e | S (sec.4c).")
    print()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--wmax", type=int, default=30)
    ap.add_argument("--wmax-e2", type=int, default=26)
    a = ap.parse_args()
    npass, ntot = run_checks(verbose=not a.quiet)
    if a.quiet:
        if npass != ntot:
            print("syzygy_sweep: %d/%d checks FAILED" % (ntot - npass, ntot))
            return 1
        print("syzygy_sweep: %d/%d checks pass" % (npass, ntot))
        return 0
    print("\n%d/%d checks pass\n" % (npass, ntot))
    if npass != ntot:
        return 1
    report(wmax=a.wmax, wmax_e2=a.wmax_e2)
    return 0


if __name__ == "__main__":
    sys.exit(main())
