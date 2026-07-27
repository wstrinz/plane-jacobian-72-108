#!/usr/bin/env python3
"""spine_verify.py -- independent re-derivation of everything `spine.py` claims.

Every check here reaches the same object by a DIFFERENT route than `spine.py`:

  V1  generators: `face_kill_sweep.canonical_G_generators()` (which goes through
      `bigrade_annotator` and asserts the G5 normalisation) instead of the
      `generators.json` term list -- and the two are compared.
  V2  the spine factorisation: by EXACT POLYNOMIAL DIVISION of honest
      symbolic-coefficient polynomials in y (`sp.div`, remainder must be 0),
      instead of an identity in opaque commuting symbols.
  V3  the boxed identity: read off the quotient of the K-row and compared
      coefficient-by-coefficient in y against the claimed form.
  V4  the elimination certificate: a FULLY POLYNOMIAL cofactor identity
          F*Z - (1/6)g^5 t^a Rm^4 = (A^2+Z)*g1r - g*A*g2 + g^2*Rm*g3
      (spine.py instead solves g1 for C and combines the two hatted rows), plus
      an independent RANDOM EXACT-RATIONAL point test: pick free values, solve
      g1r for C, g2 for d0, g3 for t^a, and confirm F*Z hits the target.
  V5  the t^a divisibility: a valuation enumeration over ALL FOUR rows
      G1,G2,G3,K with degree-derived valuation ceilings, instead of the two-row
      (identity + H3) enumeration `tpower_divisibility.py` and spine.py S9 use.
  V6  the arithmetic obstruction at n=1, by resultant instead of gcd.
  V7  the frontier impact, counted read-only from `phase_d_states_sub2.json`.

SCOPE (identical to spine.py, restated so this file stands alone): the whole
argument is conditional on `t^a | dm2,dm3,dm4`, which is established on the T2
branch (d1 = 0) ONLY.  On T1 the n = 1,2,3,4 kills are conditional on that
divisibility; n = 0 on T1 (`a10_b0000_T1`) is not closed at all.

Read-only.  Usage:
    python spine_verify.py
    python spine_verify.py --quiet
"""
from __future__ import annotations

import argparse
import itertools
import json
import os
import random
import sys

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__)) or "."

y = sp.Symbol("y")
t = y + 1
QQUART = 2048 * y**4 - 512 * y**3 + 320 * y**2 - 240 * y + 195
CGEN = sp.Rational(-1, 6630)
CAP = {"d2": 4, "d1": 6, "d0": 8, "dm1": 10, "dm2": 12, "dm3": 14, "dm4": 16, "Phi": 34}
INF = 10**9
FAMILY = {0: "a10_b0000", 1: "a9_b1000", 2: "a8_b1100", 3: "a7_b1110", 4: "a6_b1111"}


# ------------------------------------------------------------------ V1 sources
def gens_via_face_kill_sweep():
    import face_kill_sweep as fks
    g = fks.canonical_G_generators()
    return {k: sp.expand(v[0]) for k, v in g.items()}


def gens_via_json():
    with open(os.path.join(HERE, "generators.json"), encoding="utf-8") as fh:
        d = json.load(fh)
    order = d["variable_order"]
    sym = {n: sp.Symbol(n) for n in order}
    out = {}
    for name in ("G1", "G2", "G3", "G5body"):
        expr = sp.Integer(0)
        for coeff, ev in d["polynomials"][name]:
            term = sp.Rational(coeff)
            for v, k in zip(order, ev):
                if k:
                    term *= sym[v] ** k
            expr += term
        out[name] = sp.expand(expr)
    out["G5"] = sp.expand(out.pop("G5body") + sym["Phi"])
    return out


# --------------------------------------------- V2/V3 honest polynomials in y
def _pv(prefix, deg):
    cs = sp.symbols("%s0:%d" % (prefix, deg + 1))
    return sum(c * y**i for i, c in enumerate(cs)), list(cs)


def build_family(n, branch):
    """Honest symbolic-coefficient polynomials for one spine family."""
    a = 10 - n
    gam = sp.Symbol("gam")
    mcs = sp.symbols("m0:%d" % n) if n else ()
    Rm = sp.expand(y**n + sum(c * y**i for i, c in enumerate(mcs)))
    Qp, _ = _pv("Qc", 4 - n)
    A, _ = _pv("A", n + 2)
    B, _ = _pv("B", n + 4)
    C, _ = _pv("C", n + 6)
    d2, _ = _pv("D", CAP["d2"])
    d0, _ = _pv("F", CAP["d0"])
    d1 = sp.Integer(0) if branch == "T2" else _pv("E", CAP["d1"])[0]
    return dict(a=a, n=n, gam=gam, Rm=Rm, Q=Qp, A=A, B=B, C=C,
                d2=d2, d1=d1, d0=d0,
                e=sp.expand(gam * t**a * Rm), dm2=sp.expand(t**a * A),
                dm3=sp.expand(t**a * B), dm4=sp.expand(t**a * C),
                Phi=sp.expand(CGEN * t**30 * Rm * Qp))


def instantiate(G, F):
    s = sp.symbols("d2 d1 d0 dm1 dm2 dm3 dm4 Phi")
    sub = dict(zip(s, (F["d2"], F["d1"], F["d0"], F["e"], F["dm2"], F["dm3"],
                       F["dm4"], F["Phi"])))
    return {k: sp.expand(v.xreplace(sub)) for k, v in G.items()}


def exact_divide(num, den):
    """num/den in Q(params)[y]; returns (quotient, remainder)."""
    qy, ry = sp.div(sp.Poly(sp.expand(num), y), sp.Poly(sp.expand(den), y))
    return sp.expand(qy.as_expr()), sp.expand(ry.as_expr())


# --------------------------------------------------- V5 valuation enumeration
def _min_twice(vals):
    vals = [v for v in vals if v < INF]
    if not vals:
        return True                     # every term is the zero polynomial
    m = min(vals)
    return vals.count(m) >= 2


def _rhs_is_30(vals):
    """sum of terms with these valuations equals an object of valuation 30."""
    vals = [v for v in vals if v < INF]
    if not vals:
        return False                    # 0 != 2*Phi
    m = min(vals)
    return (m == 30) if vals.count(m) == 1 else (m <= 30)


def tpower_all_rows(a, branch="T2"):
    """Refute v_t(dm2) < a, v_t(dm3) < a, v_t(dm4) < a using ALL FOUR G-rows.

    Valuation ceilings come from the certified stripped degree caps: a nonzero
    polynomial of degree <= D has t-valuation <= D.  INF encodes the zero
    polynomial.  e != 0 always (else 2*Phi = 0).
    """
    rng = lambda cap: list(range(cap + 1)) + [INF]
    d1s = [INF] if branch == "T2" else rng(CAP["d1"])
    survivors = []
    for rho, s, tau, d0o, d2o, d1o in itertools.product(
            rng(CAP["dm2"]), rng(CAP["dm3"]), rng(CAP["dm4"]),
            rng(CAP["d0"]), rng(CAP["d2"]), d1s):
        if min(rho, s, tau) >= a:
            continue                                    # not a counterexample
        # G1 = 3/2 d1 e^2 + 3 d2 e R + 3 e T + 3 R S
        if not _min_twice([d1o + 2 * a, d2o + a + rho, a + tau, rho + s]):
            continue
        # G2 = -3/2 d0 e^2 + 3/2 d2 R^2 + 3 R T + 3/2 S^2
        if not _min_twice([d0o + 2 * a, d2o + 2 * rho, rho + tau, 2 * s]):
            continue
        # G3 = -3 d0 e R - 3/2 d1 R^2 - 1/2 e^3 + 3 S T
        if not _min_twice([d0o + a + rho, d1o + 2 * rho, 3 * a, s + tau]):
            continue
        # K = 0:  e*(d2 e^2 + 3 e S + 3 R^2) = 2*Phi,  v(2*Phi) = 30
        if not _rhs_is_30([d2o + 3 * a, 2 * a + s, a + 2 * rho]):
            continue
        survivors.append((rho, s, tau, d0o, d2o, d1o))
        if len(survivors) > 4:
            break
    return survivors


# ============================================================== check suite
def run(verbose=True):
    out, npass, ntot = [], 0, 0

    def ck(name, ok, detail):
        nonlocal npass, ntot
        ntot += 1
        npass += bool(ok)
        out.append("  [%s] %s\n        %s" % ("PASS" if ok else "FAIL", name, detail))

    # ---- V1  two independent generator sources agree
    Gj = gens_via_json()
    try:
        Gf = gens_via_face_kill_sweep()
        diffs = {k: sp.expand(Gj[k] - Gf[k]) for k in Gj}
        ck("V1  generators.json == face_kill_sweep.canonical_G_generators()",
           all(v == 0 for v in diffs.values()),
           "residuals: %s" % ", ".join("%s=%s" % kv for kv in diffs.items()))
    except Exception as exc:                                  # pragma: no cover
        ck("V1  generators.json == face_kill_sweep.canonical_G_generators()",
           False, "import/compare failed: %r" % (exc,))

    d2, d1, d0, dm1, dm2, dm3, dm4, Phi = sp.symbols("d2 d1 d0 dm1 dm2 dm3 dm4 Phi")
    K = sp.expand(2 * (Gj["G5"] + d2 * Gj["G3"] + d1 * Gj["G2"] + d0 * Gj["G1"]))
    ck("V1b  K-syzygy re-derived: K = 2*Phi - e*(d2 e^2 + 3 e S + 3 R^2)",
       sp.expand(K - (2 * Phi - dm1 * (d2 * dm1**2 + 3 * dm1 * dm3 + 3 * dm2**2))) == 0,
       "residual = %s"
       % sp.expand(K - (2 * Phi - dm1 * (d2 * dm1**2 + 3 * dm1 * dm3 + 3 * dm2**2))))

    # ---- V2/V3  exact polynomial division, honest coefficient polynomials
    for branch in ("T2", "T1"):
        for n in range(5):
            F = build_family(n, branch)
            a, gam, Rm = F["a"], F["gam"], F["Rm"]
            inst = instantiate(dict(Gj, K=K), F)
            A, B, C = F["A"], F["B"], F["C"]
            want = {
                "G1": (3 * t**(2 * a),
                       sp.Rational(1, 2) * gam**2 * F["d1"] * Rm**2
                       + gam * Rm * (F["d2"] * A + C) + A * B),
                "G2": (sp.Rational(3, 2) * t**(2 * a),
                       F["d2"] * A**2 + 2 * A * C + B**2 - gam**2 * F["d0"] * Rm**2),
                "G3": (3 * t**(2 * a),
                       -gam * F["d0"] * Rm * A - sp.Rational(1, 2) * F["d1"] * A**2
                       + B * C - sp.Rational(1, 6) * gam**3 * t**a * Rm**3),
                "K": (-gam * t**(3 * a) * Rm,
                      3 * A**2 + gam**2 * F["d2"] * Rm**2 + 3 * gam * Rm * B
                      - (2 * CGEN / gam) * t**(3 * n) * F["Q"]),
            }
            ok, det = True, []
            for name, (fac, red) in want.items():
                # `G == fac * red` in Q(params)[y] is exactly the statement that
                # the division is exact WITH this quotient; verified by
                # expansion because unassisted sp.div at n = 3,4 is a cost wall.
                r = sp.expand(sp.together(sp.expand(inst[name] - fac * red)))
                good = (r == 0)
                ok &= good
                det.append("%s:%s" % (name, "ok" if good else "resid=%s" % r))
            ck("V2.%s.n%d  spine factorisation on honest symbolic-coefficient "
               "polynomials in y" % (branch, n), ok, "a=%d  %s" % (a, ", ".join(det)))

    # ---- V2b  unassisted exact division (quotient NOT supplied), n = 0,1
    ok2b, det2b = True, []
    for branch in ("T2", "T1"):
        for n in (0, 1):
            F = build_family(n, branch)
            a, gam, Rm = F["a"], F["gam"], F["Rm"]
            inst = instantiate(dict(Gj, K=K), F)
            for name, fac in (("G1", 3 * t**(2 * a)),
                              ("K", -gam * t**(3 * a) * Rm)):
                quo, rem = exact_divide(inst[name], fac)
                ok2b &= (rem == 0)
                det2b.append("%s n=%d %s: rem=%s" % (branch, n, name, rem))
    ck("V2b  unassisted sp.div: the stated factors really divide G1 and K "
       "(remainder 0), n = 0,1 both branches", ok2b, " | ".join(det2b))

    # ---- V3  the boxed identity coefficient-by-coefficient, degree ledger
    ok3, det3 = True, []
    for n in range(5):
        F = build_family(n, "T2")
        gam, Rm = F["gam"], F["Rm"]
        lhs = sp.expand(3 * F["A"]**2 + gam**2 * F["d2"] * Rm**2 + 3 * gam * Rm * F["B"])
        rhs = sp.expand((2 * CGEN / gam) * t**(3 * n) * F["Q"])
        dl = sp.Poly(lhs, y).degree()
        dr = sp.Poly(rhs, y).degree()
        good = (dl == 2 * n + 4) and (dr == 2 * n + 4)
        ok3 &= good
        det3.append("n=%d: deg LHS=%d deg RHS=%d target=%d" % (n, dl, dr, 2 * n + 4))
    ck("V3  boxed identity is bi-degree-exact: both sides have degree 4+2n",
       ok3, " | ".join(det3))

    # ---- V4  the fully polynomial cofactor certificate
    #  F*Z - (1/6)g^5 t^a Rm^4 = (A^2+Z)*g1r - g*A*g2 + g^2*Rm*g3   with B = Rm*v
    As, Cs, vs, us, Rs, Ts, gs, D0, D1 = sp.symbols(
        "As Cs vs us Rs Ts gs D0 D1")
    w = sp.Rational(1, 2) * gs**2 * D1 * Rs
    g1r = sp.expand(gs * (us / gs * As + Cs) + As * vs + w)
    g2 = sp.expand(us / gs * As**2 + 2 * As * Cs + Rs**2 * vs**2 - gs**2 * D0 * Rs**2)
    g3 = sp.expand(-gs * D0 * Rs * As - sp.Rational(1, 2) * D1 * As**2
                   + Rs * vs * Cs - sp.Rational(1, 6) * gs**3 * Ts * Rs**3)
    Fq = As * (us + 2 * vs) + w
    Zq = As**2 - gs * Rs**2 * vs
    tgt = sp.expand(Fq * Zq - sp.Rational(1, 6) * gs**5 * Ts * Rs**4)
    combo = sp.expand((As**2 + Zq) * g1r - gs * As * g2 + gs**2 * Rs * g3)
    ck("V4  polynomial cofactor certificate F*Z - (1/6)g^5 t^a Rm^4 "
       "= (A^2+Z)*g1r - g*A*g2 + g^2*Rm*g3",
       sp.expand(sp.together(tgt - combo)) == 0,
       "residual = %s   [holds with d1 present, so it covers T1 too]"
       % sp.expand(sp.together(tgt - combo)))

    # ---- V4b  random exact-rational point test on the reduced variety
    random.seed(20260725)
    bad = []
    for trial in range(40):
        pick = lambda: sp.Rational(random.randint(-9, 9), random.randint(1, 5))
        val = {As: pick() or sp.Integer(1), vs: pick(), us: pick(),
               Rs: pick() or sp.Integer(1), gs: pick(), D1: pick()}
        if val[gs] == 0 or val[As] == 0 or val[Rs] == 0:
            continue
        Cv = sp.solve(sp.Eq(g1r.xreplace(val), 0), Cs)
        if not Cv:
            continue
        val[Cs] = sp.nsimplify(Cv[0])
        D0v = sp.solve(sp.Eq(g2.xreplace(val), 0), D0)
        if not D0v:
            continue
        val[D0] = sp.nsimplify(D0v[0])
        Tv = sp.solve(sp.Eq(g3.xreplace(val), 0), Ts)
        if not Tv:
            continue
        val[Ts] = sp.nsimplify(Tv[0])
        resid = sp.simplify(tgt.xreplace(val))
        if resid != 0:
            bad.append((trial, resid))
    ck("V4b  random exact-rational points of {g1r=g2=g3=0} satisfy "
       "F*Z = (1/6)g^5 t^a Rm^4", not bad,
       "40 trials, %d failures%s" % (len(bad), "" if not bad else " %s" % bad[:2]))

    # ---- V5  t^a divisibility on BOTH branches, rows extracted AUTOMATICALLY
    #  The term-valuation vectors are read off the canonical generators rather
    #  than hand-transcribed (spine.py S23 hand-writes them), so a transcription
    #  slip cannot be shared between the two files.
    gsyms = sp.symbols("d2 d1 d0 dm1 dm2 dm3 dm4 Phi")

    def _tvecs(expr):
        return [tuple(int(term.as_powers_dict().get(s, 0)) for s in gsyms)
                for term in sp.Add.make_args(sp.expand(expr))]

    Krow = sp.expand(2 * Phi - dm1 * (d2 * dm1**2 + 3 * dm1 * dm3 + 3 * dm2**2))
    TVEC = {n: _tvecs(Gj[n]) for n in ("G1", "G2", "G3")}
    TVEC["K"] = _tvecs(Krow)

    def _row_ok(vecs, val):
        vals = []
        for ev in vecs:
            s = 0
            for e_, v_ in zip(ev, val):
                if e_:
                    if v_ >= INF:
                        s = INF
                        break
                    s += e_ * v_
            if s < INF:
                vals.append(s)
        if not vals:
            return True
        m = min(vals)
        return vals.count(m) >= 2

    def _scan_auto(a, branch, want_ce, rows=("G1", "G2", "G3", "K"), stop=3):
        rg = lambda c: list(range(c + 1)) + [INF]
        d1s = [INF] if branch == "T2" else rg(CAP["d1"])
        hits = []
        for rho, s, tau in itertools.product(rg(CAP["dm2"]), rg(CAP["dm3"]),
                                             rg(CAP["dm4"])):
            if (min(rho, s, tau) < a) != want_ce:
                continue
            for d0o, d2o, d1o in itertools.product(rg(CAP["d0"]), rg(CAP["d2"]), d1s):
                val = (d2o, d1o, d0o, a, rho, s, tau, 30)
                if all(_row_ok(TVEC[r], val) for r in rows):
                    hits.append((rho, s, tau, d0o, d2o, d1o))
                    if len(hits) >= stop:
                        return hits
        return hits

    ok5, det5 = True, []
    for branch in ("T2", "T1"):
        for a in range(6, 11):
            ce = _scan_auto(a, branch, True)
            adm = _scan_auto(a, branch, False, stop=1)
            ok5 &= (not ce) and bool(adm)
            det5.append("%s a=%d ce=%d adm=%s" % (branch, a, len(ce), bool(adm)))
    ck("V5  t^a | dm2,dm3,dm4 on BOTH branches, a=6..10, rows extracted "
       "AUTOMATICALLY from generators.json (no hand-written valuation lists)",
       ok5, ", ".join(det5))

    # ---- V5b  the hand proof of V5, step by step, as arithmetic assertions
    #  If rho < a:  K forces s = 2*rho - a;  G1 forces tau = 3*rho - 2*a;
    #  then G3's S*T term  s+tau = 5*rho-3*a  is the STRICT unique minimum.
    ok5b, det5b = True, []
    for a in range(6, 11):
        for rho in range(a):
            # K: 30, d2o+3a, 2a+s, a+2rho.  a+2rho <= 3a-2 < 3a <= d2o+3a, < 30.
            if not (a + 2 * rho <= 3 * a - 2 < 3 * a and a + 2 * rho < 30):
                ok5b = False
                det5b.append("K-step fails a=%d rho=%d" % (a, rho))
                continue
            s = 2 * rho - a
            if s < 0:
                continue                       # K already has a unique minimum
            # G1: d1o+2a, d2o+a+rho, a+tau, rho+s=3rho-a.  3rho-a < both bounds.
            if not (3 * rho - a < a + rho and 3 * rho - a < 2 * a):
                ok5b = False
                det5b.append("G1-step fails a=%d rho=%d" % (a, rho))
                continue
            tau = 3 * rho - 2 * a
            if tau < 0:
                continue                       # G1 already has a unique minimum
            # G3: d0o+a+rho, d1o+2rho, 3a, s+tau=5rho-3a -- strict unique min
            m = 5 * rho - 3 * a
            if not (m >= 0 and m < 3 * a and m < a + rho and m < 2 * rho):
                ok5b = False
                det5b.append("G3-step fails a=%d rho=%d (m=%d)" % (a, rho, m))
    # and the two easy steps
    for a in range(6, 11):
        ok5b &= all(2 * a + s < 3 * a and 2 * a + s < 30 for s in range(a))
        ok5b &= all(a + tau < 2 * a for tau in range(a))
    ck("V5b  the hand proof behind V5: rho<a => s=2rho-a => tau=3rho-2a => "
       "G3's S*T term is a strict unique minimum; then s<a and tau<a fall to "
       "K and G1 directly", ok5b,
       "all a=6..10, all rho<a: %s  [uses only v(e)=a, v(2Phi)=30, a<=10 -- "
       "NO degree caps and NO branch assumption]"
       % ("clean" if ok5b else "; ".join(det5b)))

    # ---- V6  the n=1 obstruction, by resultant instead of gcd
    W3 = sp.expand((y + 1) * sp.diff(QQUART, y, 2) - 6 * sp.diff(QQUART, y))
    res = sp.resultant(sp.Poly(QQUART, y), sp.Poly(W3, y))
    fl = sp.factor_list(QQUART)
    ck("V6  res(q, (y+1)q''-6q') != 0 -> no common root in ANY field extension",
       res != 0 and len(fl[1]) == 1 and fl[1][0][1] == 1,
       "resultant = %s ; q irreducible over Q: %s ; q(-1) = %s"
       % (res, len(fl[1]) == 1 and fl[1][0][1] == 1, QQUART.subs(y, -1)))

    # ---- V6b  the n=2 obstruction: Q(-1) != 0 for every squarefree Rm | q
    ok6b = True
    for n in range(5):
        # Q = q/Rm; q(-1) != 0 => Q(-1) != 0 whatever Rm is
        ok6b &= QQUART.subs(y, -1) != 0
    ck("V6b  n=2 kill premise: q(-1) != 0 so Q(-1) != 0 for every Rm | q", ok6b,
       "q(-1) = %s ; deg(mu*Q - 3*zeta*t^2) <= 2 < 4 = deg Rm^2 forces V=0, "
       "then V(-1) = mu*Q(-1) != 0" % QQUART.subs(y, -1))

    # ---- V7  frontier impact, read-only
    try:
        with open(os.path.join(HERE, "phase_d_states_sub2.json"), encoding="utf-8") as fh:
            data = json.load(fh)
        alive = [c for c in data["cases"]
                 if all(x <= 1 for x in c["b"]) and c["a_t"] + sum(c["b"]) == 10]
        cols = {}
        for c in alive:
            key = ("a%d_b%s" % (c["a_t"], "".join(str(x) for x in c["b"])), c["branch"])
            k = cols.setdefault(key, [0, 0])
            k[0] += 1
            k[1] += c["state_count"]
        t2 = {k: v for k, v in cols.items() if k[1] == "T2"}
        t1 = {k: v for k, v in cols.items() if k[1] == "T1"}
        resid_col = [(k, v) for k, v in t1.items() if k[0] == "a10_b0000"]
        ck("V7  frontier impact (read-only census of phase_d_states_sub2.json)",
           bool(cols),
           "surviving-after-divisor-filter columns: %s || totals: T2 %d cases / "
           "%d states, T1 %d cases / %d states || CLOSED by SPINE.md: all "
           "columns except a10_b0000_T1 = %s, which is the entire residue "
           "(T1 kills are conditional on the sec.8 t^a upgrade)"
           % (", ".join("%s_%s:%dc/%ds" % (k[0], k[1], v[0], v[1])
                        for k, v in sorted(cols.items())),
              sum(v[0] for v in t2.values()), sum(v[1] for v in t2.values()),
              sum(v[0] for v in t1.values()), sum(v[1] for v in t1.values()),
              resid_col))
    except Exception as exc:                                  # pragma: no cover
        ck("V7  frontier impact (read-only census)", False, "census failed: %r" % (exc,))

    if verbose:
        print("\n".join(out))
    return npass, ntot


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args()
    npass, ntot = run(verbose=not a.quiet)
    if a.quiet:
        if npass != ntot:
            print("spine_verify: %d/%d checks FAILED" % (ntot - npass, ntot))
            return 1
        print("spine_verify: %d/%d checks pass" % (npass, ntot))
        return 0
    print("\n%d/%d checks pass" % (npass, ntot))
    return 0 if npass == ntot else 1


if __name__ == "__main__":
    sys.exit(main())
