#!/usr/bin/env python3
"""spine.py -- the forcing-divisor spine: the five-family reduction of the sub2
frontier, and the emptiness argument it carries.

CONTEXT.  `DIVISOR_SYZYGY.md` established the universal K-syzygy

    2*(G5 + d2*G3 + d1*G2 + d0*G1)  ==  2*Phi - e*(d2*e^2 + 3*e*S + 3*R^2)

hence `e | 2*Phi`, hence (sub2) `deg e = 10`, `e = gamma*t^a*Rm` with
`t = y+1`, `Rm` a squarefree divisor of the quartic `q`, and `a + deg Rm = 10`.
The survivors of that filter are precisely the cells the lemma cannot see:
`a + sum(b_i) = 10`, `b_i in {0,1}`, i.e. FIVE support families indexed by

    n := deg Rm = 10 - a   in   {0,1,2,3,4}
      a10_b0000 (n=0)  a9_b1000 (n=1)  a8_b1100 (n=2)  a7_b1110 (n=3)  a6_b1111 (n=4)

WHAT THIS FILE DOES.

(1) DERIVES, from the canonical generators in `generators.json` (never from a
    quoted reduced form), the spine parametrisation

        e = gamma*t^a*Rm,  dm2 = t^a*A,  dm3 = t^a*B,  dm4 = t^a*C,   a = 10-n
        deg A <= n+2,  deg B <= n+4,  deg C <= n+6,   q = Rm*Q

    and verifies the EXACT factorisations

        G1 =  3   * t^(2a) * g1,          K  = -gamma * t^(3a) * Rm * kbox
        G2 = (3/2)* t^(2a) * g2
        G3 =  3   * t^(2a) * g3

    with, in particular, the BOXED IDENTITY  kbox = 0, i.e.

        3*A^2 + gamma^2*d2*Rm^2 + 3*gamma*Rm*B  ==  (2c/gamma) * t^(3n) * Q

    every term of degree exactly <= 4 + 2n.  Residual 0 for every n = 0..4.

(2) Confirms the spare collapse the brief asked about, 45 -> 45 - 3a =
    15/18/21/24/27 (S6), and then SUPERSEDES it: once dm4 is eliminated (it is
    determined, S19) and dm3 = e*Sbar (S20), the honest free-spare count is
    45 -> n+8 = 8/9/10/11/12 (S22).  Both results are reported, not conflated.

(3) Runs the STRUCTURAL argument on the reduced system.  Eliminating C (from
    g1) and d0 (from g2, g3) gives the exact certificate

        F * Z  -  (1/6)*gamma^5*t^a*Rm^4  =  -gamma*A*g2hat + gamma^2*Rm*g3hat
        F := A*(u + 2*v) + (1/2)*gamma^2*d1*Rm,      Z := A^2 - gamma*Rm^2*v
        u := gamma*d2,   B = Rm*v   (v exists: see below),   deg u,v <= 4

    so on the variety  F*Z = (1/6)*gamma^5*t^a*Rm^4.  The three degree caps
    deg A <= n+2, deg F <= n+6, deg Z <= 2n+4 sum to exactly deg(t^a*Rm^4) =
    10+3n, so every one of them is ATTAINED.  Since kbox forces A(r_i) != 0 at
    each marked root, Z(r_i) = A(r_i)^2 != 0, so gcd(Z, Rm) = 1 and Z | t^a:

        2n+4 = deg Z <= a = 10-n   =>   n <= 2      [kills n = 3, 4]

    and Z = zeta*t^(2n+4).  Feeding that back into kbox gives

        gamma*Rm^2*(u+6v) = t^(3n) * ( mu*Q - 3*zeta*t^(4-n) ),   mu := 2c/gamma

    so Rm^2 divides  V_n := mu*Q - 3*zeta*t^(4-n)  (deg V_n <= 4-n):

        n=2: deg V < deg Rm^2 => V=0 => mu*Q(-1) = 0    CONTRADICTION
        n=1: V = Rm^2*(linear); the 2x2 system in (mu,zeta) has a nonzero
             solution only if  (r+1)*q''(r) = 6*q'(r);  but
             gcd(q, (y+1)*q'' - 6*q') = 1 and q is irreducible   CONTRADICTION
        n=0: Rm = 1 -- no condition.  n=0 needs the T2-only step below.

    On T2 (d1 = 0) F = A*(u+2v), so A | t^a*Rm^4 too, hence A = lambda*t^(n+2),
    hence Z | t^(8-2n) and 2n+4 <= 8-2n forces n <= 1; and at n = 0 the
    resulting A, B, d2 are all multiples of t^4, so kbox at y = -1 reads
    mu*q(-1) = 0.   CONTRADICTION.

SCOPE -- READ THIS BEFORE QUOTING ANY VERDICT.

  * The whole argument rests on  t^a | dm2, dm3, dm4.  The RECORDED status of
    that lemma (GSYSTEM_CELL.md sec.7.9, DIVISOR_SYZYGY.md sec.4) is: proved on
    the T2 branch (d1 = 0) only, open on T1.  Check S23 UPGRADES it to BOTH
    branches, a = 6..10, by a four-row valuation enumeration (G1,G3,K
    load-bearing; G2 redundant) with an admissibility control and a row
    ablation; SPINE.md sec.8.1 gives the hand proof, and spine_verify.py V5/V5b
    re-derive it with the rows extracted automatically from generators.json.
    THIS CONTRADICTS A RECORDED STATUS and is flagged for adjudication.  The
    whole T1 column of the verdict table depends on it.
  * The n = 0 kill is T2-only by construction (it uses d1 = 0), so
    `a10_b0000_T1` is NOT closed by anything here -- it is the entire residue.
  * Everything is exact rational arithmetic in char 0.  No Groebner basis, no
    modular reconnaissance, no ansatz for sigma -- only the certified sub2
    window caps deg d2 <= 4, deg d1 <= 6, deg dm2 <= 12, deg dm3 <= 14,
    deg e <= 10 (`full_system_bridge.WEIGHT`/`STRIP_DEGCAP`).  deg d0 is not
    used at all: d0 is eliminated.  The t^a lemma itself uses NO caps.

Read-only.  Usage:
    python spine.py            # full report
    python spine.py --quiet    # self-check, exit 0 iff all checks pass
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__)) or "."

# ---------------------------------------------------------------- primitives
y = sp.Symbol("y")
t = y + 1
Q_QUARTIC = 2048 * y**4 - 512 * y**3 + 320 * y**2 - 240 * y + 195
C_GENUINE = sp.Rational(-1, 6630)          # Phi = c * t^30 * q
DEG_E_SUB2 = 10
# certified sub2 stripped window caps (weights 12k, stripped cap 2k)
WEIGHT = {"d2": 24, "d1": 36, "d0": 48, "dm1": 60,
          "dm2": 72, "dm3": 84, "dm4": 96, "Phi": 204}
CAP = {k: 2 * (w // 12) for k, w in WEIGHT.items()}      # d2:4 d1:6 d0:8 dm1:10 ...

FAMILY = {0: "a10_b0000", 1: "a9_b1000", 2: "a8_b1100",
          3: "a7_b1110", 4: "a6_b1111"}


def load_canonical_G():
    """G1,G2,G3,G5 = G5body + Phi, parsed from the canonical generators.json.

    This is the same artifact `full_system_bridge.gsystem()` reads; it is the
    de-pickled term list, not a hand transcription.  The `G5body + Phi`
    normalisation is asserted (a stale `2*Phi` form was a real bug here).
    """
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
    return out, sym


def K_from_syzygy(G, sym):
    """K := 2*(G5 + d2*G3 + d1*G2 + d0*G1), expanded -- DERIVED, not quoted."""
    return sp.expand(2 * (G["G5"] + sym["d2"] * G["G3"]
                          + sym["d1"] * G["G2"] + sym["d0"] * G["G1"]))


# ------------------------------------------------------- the spine substitution
#
# A, B, C, Rm, Q, d0, d1, d2 are treated as OPAQUE commuting symbols.  Every
# factorisation below is then a polynomial identity in the free commutative ring
# Q[gamma, T, A, B, C, Rm, Q, d0, d1, d2] with T standing for t = y+1, and holds
# for any polynomial values whatsoever.  Degree bookkeeping is done separately
# from the certified caps (check S4) and re-instantiated concretely in S5.
Tsym = sp.Symbol("T_")            # stands for t = y+1
gam = sp.Symbol("gamma_")         # sympy's `gamma` is the Gamma function
Asym, Bsym, Csym = sp.symbols("A_ B_ C_")
Rmsym, Qsym = sp.symbols("Rm_ Q_")
vsym = sp.Symbol("v_")            # B = Rm * v
usym = sp.Symbol("u_")            # u = gamma*d2


def spine_subs(n, sym, branch="T2"):
    """dict substituting the spine parametrisation into the canonical symbols."""
    a = DEG_E_SUB2 - n
    return {
        sym["dm1"]: gam * Tsym**a * Rmsym,
        sym["dm2"]: Tsym**a * Asym,
        sym["dm3"]: Tsym**a * Bsym,
        sym["dm4"]: Tsym**a * Csym,
        sym["Phi"]: C_GENUINE * Tsym**30 * Rmsym * Qsym,
        sym["d1"]: sp.Integer(0) if branch == "T2" else sym["d1"],
    }


def reduced_generators(n, sym, branch="T2"):
    """The REDUCED rows g1,g2,g3,kbox for family n, written out for comparison.

    These are the forms the factorisation checks target; the checks derive the
    left-hand sides from `generators.json` and confirm the quotient is exactly
    this, with zero residual.  Nothing here is assumed.
    """
    a = DEG_E_SUB2 - n
    d0, d1, d2 = sym["d0"], (sp.Integer(0) if branch == "T2" else sym["d1"]), sym["d2"]
    mu = 2 * C_GENUINE / gam
    g1 = (sp.Rational(1, 2) * gam**2 * d1 * Rmsym**2
          + gam * Rmsym * (d2 * Asym + Csym) + Asym * Bsym)
    g2 = d2 * Asym**2 + 2 * Asym * Csym + Bsym**2 - gam**2 * d0 * Rmsym**2
    g3 = (-gam * d0 * Rmsym * Asym - sp.Rational(1, 2) * d1 * Asym**2
          + Bsym * Csym - sp.Rational(1, 6) * gam**3 * Tsym**a * Rmsym**3)
    kbox = (3 * Asym**2 + gam**2 * d2 * Rmsym**2 + 3 * gam * Rmsym * Bsym
            - mu * Tsym**(3 * n) * Qsym)
    return {"g1": sp.expand(g1), "g2": sp.expand(g2),
            "g3": sp.expand(g3), "kbox": sp.expand(kbox)}


FACTOR = {"G1": (3, 2), "G2": (sp.Rational(3, 2), 2), "G3": (3, 2)}   # (coeff, a-power)


def factorisation_residuals(n, branch="T2"):
    """G_i - coeff*T^(2a)*g_i  and  K + gamma*T^(3a)*Rm*kbox.  All must be 0."""
    G, sym = load_canonical_G()
    K = K_from_syzygy(G, sym)
    a = DEG_E_SUB2 - n
    subs = spine_subs(n, sym, branch)
    red = reduced_generators(n, sym, branch)
    res = {}
    for name in ("G1", "G2", "G3"):
        coeff, _ = FACTOR[name]
        lhs = sp.expand(G[name].xreplace(subs))
        res[name] = sp.expand(lhs - coeff * Tsym**(2 * a) * red[name.replace("G", "g")])
    lhsK = sp.expand(K.xreplace(subs))
    res["K"] = sp.expand(lhsK + gam * Tsym**(3 * a) * Rmsym * red["kbox"])
    return res, red


# ------------------------------------------------------- the elimination certificate
def elimination_certificate(n, branch="T2"):
    """Verify  F*Z - (1/6)*gamma^5*T^a*Rm^4  ==  -gamma*A*g2hat + gamma^2*Rm*g3hat.

    g1 is linear in C with leading coefficient gamma*Rm, so on the variety
    (gamma != 0, Rm != 0)  C = -(A*(u+v) + w)/gamma  with w = (1/2)*gamma^2*d1*Rm.
    g2hat, g3hat are g2, g3 after B -> Rm*v and that substitution for C.
    """
    a = DEG_E_SUB2 - n
    _, sym = load_canonical_G()
    d0 = sym["d0"]
    d1 = sp.Integer(0) if branch == "T2" else sym["d1"]
    red = reduced_generators(n, sym, branch)
    w = sp.Rational(1, 2) * gam**2 * d1 * Rmsym
    Cval = -(Asym * (usym + vsym) + w) / gam
    sub = {Bsym: Rmsym * vsym, Csym: Cval, sym["d2"]: usym / gam}
    g1hat = sp.expand(sp.together(red["g1"].xreplace(sub)))
    g2hat = sp.expand(sp.together(red["g2"].xreplace(sub)))
    g3hat = sp.expand(sp.together(red["g3"].xreplace(sub)))
    F = Asym * (usym + 2 * vsym) + w
    Z = Asym**2 - gam * Rmsym**2 * vsym
    target = sp.expand(F * Z - sp.Rational(1, 6) * gam**5 * Tsym**a * Rmsym**4)
    combo = sp.expand(-gam * Asym * g2hat + gam**2 * Rmsym * g3hat)
    return {"g1hat": sp.simplify(g1hat), "resid": sp.expand(sp.together(target - combo)),
            "F": F, "Z": Z, "d0_present": d0 in (target - combo).free_symbols}


# ------------------------------------------------------------------ degree ledger
def degree_ledger(n, branch="T2"):
    """Every degree used by the structural argument, from the certified caps."""
    a = DEG_E_SUB2 - n
    degA = CAP["dm2"] - a                 # dm2 = t^a * A
    degB = CAP["dm3"] - a
    degC = CAP["dm4"] - a
    degv = degB - n                       # B = Rm * v
    degu = CAP["d2"]                      # u = gamma*d2
    degF = max(degA + max(degu, degv),
               (0 if branch == "T2" else CAP["d1"]) + n)
    degZ = max(2 * degA, n * 2 + degv)
    rhs = a + 4 * n                       # deg( t^a * Rm^4 )
    return {"a": a, "n": n, "degA": degA, "degB": degB, "degC": degC,
            "degv": degv, "degu": degu, "degF": degF, "degZ": degZ,
            "deg_rhs": rhs, "cap_sum": degF + degZ, "kbox_deg": 2 * n + 4}


# ------------------------------------------------ the n=1 arithmetic obstruction
def n1_obstruction():
    """(r+1)*q''(r) = 6*q'(r) with q(r) = 0 is impossible: gcd(q, W3) = 1."""
    qq = sp.Poly(Q_QUARTIC, y)
    W3 = sp.Poly(sp.expand((y + 1) * sp.diff(Q_QUARTIC, y, 2)
                           - 6 * sp.diff(Q_QUARTIC, y)), y)
    g = sp.gcd(qq, W3)
    fl = sp.factor_list(Q_QUARTIC)
    irred = len(fl[1]) == 1 and fl[1][0][1] == 1
    return {"W3": W3.as_expr(), "gcd": g.as_expr(), "trivial": g.degree() == 0,
            "q_irreducible": irred, "q_at_-1": Q_QUARTIC.subs(y, -1),
            "q_squarefree": sp.gcd(qq, sp.Poly(sp.diff(Q_QUARTIC, y), y)).degree() == 0}


# ------------------------------------------------------------ the spare collapse
def spare_collapse():
    rows = []
    full = sum(CAP[k] + 1 for k in ("dm2", "dm3", "dm4"))
    for n in range(5):
        a = DEG_E_SUB2 - n
        per = {k: CAP[k] - a + 1 for k in ("dm2", "dm3", "dm4")}
        red = sum(per.values())
        rows.append((n, a, FAMILY[n], full, per, red, full - 3 * a, red == full - 3 * a))
    return rows


# ------------------------------ independent re-implementation of the t^a lemma
INF = 10**9


def _tpow_ok(a, cap_R=12, cap_S=14, dmax=64):
    """Exhaustive t-valuation refutation of v(dm2) < a and v(dm3) < a on d1 = 0.

    Written from the two relations directly, not imported:
      identity   d2*e^3 + 3*e^2*S + 3*e*R^2 = 2*Phi,  v_t(2*Phi) = 30
      H3 = 0     -6*d0*e^2*R - 6*d2*e*R*S - e^4 - 6*R*S^2 = 0   (d1 = 0)
    A sum of terms with a UNIQUE minimal valuation cannot vanish, and if it
    equals a nonzero object of valuation 30 the unique minimum must BE 30.
    """
    def ident(a_, rho, s, d2o):
        o = [x for x in (d2o + 3 * a_, 2 * a_ + s, a_ + 2 * rho) if x < INF]
        if not o:
            return False
        m = min(o)
        return (m == 30) if o.count(m) == 1 else (m <= 30)

    def h3(a_, rho, s, d0o, d2o):
        o = [x for x in (d0o + 2 * a_ + rho, d2o + a_ + rho + s, 4 * a_, rho + 2 * s)
             if x < INF]
        m = min(o)
        return o.count(m) >= 2

    for rho in range(a):
        for s in list(range(cap_S + 1)) + [INF]:
            for d2o in range(dmax):
                if not ident(a, rho, s, d2o):
                    continue
                for d0o in range(dmax):
                    if h3(a, rho, s, d0o, d2o):
                        return False, ("rho", rho, s, d2o, d0o)
    for s in range(a):
        for rho in range(a, cap_R + 1):
            for d2o in range(dmax):
                if ident(a, rho, s, d2o):
                    return False, ("s", s, rho, d2o)
    return True, None


# ============================================================== the check suite
def run(verbose=True):
    out, npass, ntot = [], 0, 0

    def ck(name, ok, detail):
        nonlocal npass, ntot
        ntot += 1
        npass += bool(ok)
        out.append("  [%s] %s\n        %s" % ("PASS" if ok else "FAIL", name, detail))

    G, sym = load_canonical_G()

    # ---- S0  the primitives themselves
    ck("S0  canonical G5 normalisation is G5body + Phi",
       sp.expand(G["G5"]).coeff(sym["Phi"]) == 1,
       "coeff(G5, Phi) = %s  (a stale 2*Phi form breaks every check below)"
       % sp.expand(G["G5"]).coeff(sym["Phi"]))

    # ---- S1  the syzygy, DERIVED from generators.json
    K = K_from_syzygy(G, sym)
    e, R_, S_ = sym["dm1"], sym["dm2"], sym["dm3"]
    Kshape = sp.expand(2 * sym["Phi"] - e * (sym["d2"] * e**2 + 3 * e * S_ + 3 * R_**2))
    ck("S1  K := 2*(G5+d2*G3+d1*G2+d0*G1) = 2*Phi - e*(d2*e^2+3*e*S+3*R^2)",
       sp.expand(K - Kshape) == 0,
       "residual = %s   [derived from generators.json, not quoted]"
       % sp.expand(K - Kshape))

    # ---- S2  the spine factorisation, per family, per branch
    for branch in ("T2", "T1"):
        for n in range(5):
            res, _ = factorisation_residuals(n, branch)
            allz = all(v == 0 for v in res.values())
            ck("S2.%s.n%d  spine factorisation exact (G1,G2,G3,K)" % (branch, n), allz,
               "n=%d a=%d branch=%s  residuals: %s"
               % (n, DEG_E_SUB2 - n, branch,
                  ", ".join("%s=%s" % (k, v) for k, v in res.items())))

    # ---- S3  the boxed identity, spelled out
    for n in range(5):
        _, red = factorisation_residuals(n, "T2")
        want = sp.expand(3 * Asym**2 + gam**2 * sym["d2"] * Rmsym**2
                         + 3 * gam * Rmsym * Bsym
                         - (2 * C_GENUINE / gam) * Tsym**(3 * n) * Qsym)
        ck("S3.n%d  boxed identity 3A^2 + g^2*d2*Rm^2 + 3g*Rm*B = (2c/g)*t^(3n)*Q" % n,
           sp.expand(red["kbox"] - want) == 0,
           "residual = %s ; 2c/gamma = %s ; t-power = t^%d"
           % (sp.expand(red["kbox"] - want), 2 * C_GENUINE / gam, 3 * n))

    # ---- S4  the degree ledger: every boxed term has degree exactly <= 2n+4
    ok4, rows4 = True, []
    for n in range(5):
        L = degree_ledger(n, "T2")
        terms = [2 * L["degA"], CAP["d2"] + 2 * n, n + L["degB"]]     # 3A^2, d2*Rm^2, Rm*B
        rhs = 3 * n + (4 - n)
        good = all(x == L["kbox_deg"] for x in terms) and rhs == L["kbox_deg"]
        ok4 &= good
        rows4.append("n=%d: terms %s  rhs %d  target %d %s"
                     % (n, terms, rhs, L["kbox_deg"], "ok" if good else "MISMATCH"))
    ck("S4  every boxed term has degree exactly 4+2n (deg t^(3n)Q = 3n+(4-n))",
       ok4, " | ".join(rows4))

    # ---- S5  concrete polynomial instantiation (belt and braces)
    ok5, det5 = True, []
    for n in range(5):
        a = DEG_E_SUB2 - n
        rng = [sp.Rational(k) for k in (3, -5, 7, -2, 11, 13, -17, 19, 23, -29,
                                        31, 37, -41, 43, 47, -53, 59, 61, -67, 71)]
        def poly(deg, off):
            return sum(rng[(off + i) % len(rng)] * y**i for i in range(deg + 1))
        Rmp = sp.expand(y**n + sum(rng[i] * y**i for i in range(n)))
        Qp = poly(4 - n, 5)
        Ap, Bp, Cp = poly(n + 2, 1), poly(n + 4, 3), poly(n + 6, 7)
        d2p, d1p, d0p = poly(4, 11), poly(6, 2), poly(8, 13)
        gv = sp.Rational(5, 3)
        real = {sym["dm1"]: gv * t**a * Rmp, sym["dm2"]: t**a * Ap,
                sym["dm3"]: t**a * Bp, sym["dm4"]: t**a * Cp,
                sym["Phi"]: C_GENUINE * t**30 * Rmp * Qp,
                sym["d0"]: d0p, sym["d1"]: d1p, sym["d2"]: d2p}
        redp = reduced_generators(n, sym, "T1")
        conc = {gam: gv, Tsym: t, Asym: Ap, Bsym: Bp, Csym: Cp,
                Rmsym: Rmp, Qsym: Qp, sym["d0"]: d0p, sym["d1"]: d1p,
                sym["d2"]: d2p}
        for name in ("G1", "G2", "G3"):
            coeff, _ = FACTOR[name]
            lhs = sp.expand(G[name].xreplace(real))
            rhs = sp.expand(coeff * t**(2 * a) * redp[name.replace("G", "g")].xreplace(conc))
            r = sp.expand(lhs - rhs)
            ok5 &= (r == 0)
        Kl = sp.expand(K.xreplace(real))
        Kr = sp.expand(-gv * t**(3 * a) * Rmp * redp["kbox"].xreplace(conc))
        r = sp.expand(Kl - Kr)
        ok5 &= (r == 0)
        det5.append("n=%d ok" % n if r == 0 else "n=%d K-resid %s" % (n, r))
    ck("S5  concrete exact-rational instantiation of all 5 families reproduces "
       "the same factorisation", ok5, " | ".join(det5))

    # ---- S6  spare collapse 45 -> 45-3a
    rows = spare_collapse()
    ok6 = all(r[-1] for r in rows) and [r[5] for r in rows] == [15, 18, 21, 24, 27]
    ck("S6  spare collapse 45 -> 45-3a = 15/18/21/24/27 for a=10/9/8/7/6", ok6,
       " | ".join("%s n=%d a=%d: %d->%d (45-3a=%d)" % (r[2], r[0], r[1], r[3], r[5], r[6])
                  for r in rows))

    # ---- S7  the elimination certificate F*Z = (1/6) g^5 t^a Rm^4
    for branch in ("T2", "T1"):
        for n in range(5):
            cert = elimination_certificate(n, branch)
            ck("S7.%s.n%d  F*Z - (1/6)g^5 t^a Rm^4 = -g*A*g2hat + g^2*Rm*g3hat"
               % (branch, n), cert["resid"] == 0,
               "residual = %s ; d0 eliminated: %s" % (cert["resid"], not cert["d0_present"]))

    # ---- S8  degree exactness: deg F + deg Z == deg(t^a Rm^4)
    ok8, rows8 = True, []
    for branch in ("T2", "T1"):
        for n in range(5):
            L = degree_ledger(n, branch)
            good = (L["cap_sum"] == L["deg_rhs"]) and L["degZ"] == 2 * n + 4
            ok8 &= good
            rows8.append("%s n=%d: degF=%d degZ=%d sum=%d rhs=%d %s"
                         % (branch, n, L["degF"], L["degZ"], L["cap_sum"],
                            L["deg_rhs"], "ok" if good else "MISMATCH"))
    ck("S8  cap sum deg F + deg Z equals deg(t^a Rm^4) exactly -> both attained",
       ok8, " | ".join(rows8))

    # ---- S9  the t^a divisibility, re-implemented independently (T2 only)
    ok9, det9 = True, []
    for a in (6, 7, 8, 9, 10):
        good, why = _tpow_ok(a)
        ok9 &= good
        det9.append("a=%d:%s" % (a, "OK" if good else "SURVIVOR %s" % (why,)))
    ck("S9  t^a | dm2,dm3,dm4 by independent exhaustive valuation enumeration "
       "(T2, d1=0), a = 6..10", ok9,
       "%s  [T1 is NOT covered -- the d1 terms change the case analysis]"
       % ", ".join(det9))

    # ---- S10  Z | t^a  =>  n <= 2   (both branches)
    ok10 = all((2 * n + 4 <= DEG_E_SUB2 - n) == (n <= 2) for n in range(5))
    ck("S10  deg Z = 2n+4 and Z | t^a force 2n+4 <= 10-n, i.e. n <= 2", ok10,
       " | ".join("n=%d: 2n+4=%d vs a=%d -> %s"
                  % (n, 2 * n + 4, DEG_E_SUB2 - n,
                     "OK" if 2 * n + 4 <= DEG_E_SUB2 - n else "IMPOSSIBLE (kill)")
                  for n in range(5)))

    # ---- S11  the kbox feedback  gamma*Rm^2*(u+6v) = t^(3n)*(mu*Q - 3*zeta*t^(4-n))
    zeta = sp.Symbol("zeta_")
    ok11, rows11 = True, []
    for n in range(5):
        kb = reduced_generators(n, sym, "T2")["kbox"]
        # substitute A^2 = zeta*t^(2n+4) + gamma*Rm^2*v  and  B = Rm*v, d2 = u/gamma
        lhs = sp.expand(kb.xreplace({Asym**2: zeta * Tsym**(2 * n + 4)
                                     + gam * Rmsym**2 * vsym,
                                     Bsym: Rmsym * vsym,
                                     sym["d2"]: usym / gam}))
        want = sp.expand(3 * zeta * Tsym**(2 * n + 4) + gam * Rmsym**2 * (usym + 6 * vsym)
                         - (2 * C_GENUINE / gam) * Tsym**(3 * n) * Qsym)
        r = sp.expand(sp.together(lhs - want))
        ok11 &= (r == 0)
        rows11.append("n=%d resid=%s" % (n, r))
    ck("S11  kbox with A^2 = zeta*t^(2n+4)+g*Rm^2*v becomes "
       "g*Rm^2*(u+6v) = t^(3n)*(mu*Q - 3*zeta*t^(4-n))", ok11, " | ".join(rows11))

    # ---- S12  n=2: deg V_n < deg Rm^2 forces V=0, contradicted at y=-1
    ok12 = (4 - 2) < 2 * 2 and Q_QUARTIC.subs(y, -1) != 0
    ck("S12  n=2 kill: deg(mu*Q - 3*zeta*t^2) <= 2 < 4 = deg Rm^2 so V=0, but "
       "V(-1) = mu*Q(-1) != 0", ok12,
       "deg V_2 <= %d, deg Rm^2 = %d, q(-1) = %s (so Q(-1) != 0 for any Rm | q)"
       % (4 - 2, 4, Q_QUARTIC.subs(y, -1)))

    # ---- S13  n=1: the arithmetic obstruction
    ob = n1_obstruction()
    ck("S13  n=1 kill: (r+1)q''(r) = 6q'(r) has no root in common with q",
       ob["trivial"] and ob["q_irreducible"] and ob["q_squarefree"]
       and ob["q_at_-1"] != 0,
       "W3 = (y+1)q'' - 6q' = %s ; gcd(q, W3) = %s ; q irreducible = %s ; "
       "q squarefree = %s ; q(-1) = %s"
       % (ob["W3"], ob["gcd"], ob["q_irreducible"], ob["q_squarefree"], ob["q_at_-1"]))

    # ---- S14  n=1 determinant: the 2x2 system in (mu, zeta) is what S13 answers
    r = sp.Symbol("r_")
    Qr, Qpr = sp.symbols("Qr_ Qpr_")           # Q(r), Q'(r)
    M = sp.Matrix([[Qr, -3 * (r + 1)**3], [Qpr, -9 * (r + 1)**2]])
    det = sp.factor(sp.expand(M.det()))
    ck("S14  n=1: V = Rm^2*(linear) forces det[[Q(r),-3(r+1)^3],[Q'(r),-9(r+1)^2]] = 0",
       sp.expand(det - 3 * (r + 1)**2 * ((r + 1) * Qpr - 3 * Qr)) == 0,
       "det = %s ; with Q(r)=q'(r), Q'(r)=q''(r)/2 this is (r+1)q''(r) = 6q'(r)" % det)

    # ---- S15  T2-only: A | t^a, so A = lambda*t^(n+2), and Z | t^(8-2n) => n <= 1
    ok15 = all((2 * n + 4 <= 8 - 2 * n) == (n <= 1) for n in range(5))
    ck("S15  T2 only (d1=0, F = A*(u+2v)): A | t^a gives A = lam*t^(n+2), then "
       "Z | t^(8-2n) forces n <= 1", ok15,
       " | ".join("n=%d: degZ=%d vs 8-2n=%d -> %s"
                  % (n, 2 * n + 4, 8 - 2 * n,
                     "OK" if 2 * n + 4 <= 8 - 2 * n else "IMPOSSIBLE (kill)")
                  for n in range(5)))

    # ---- S16  n=0 on T2: everything is a multiple of t^4 and kbox(-1) = mu*q(-1)
    lam, ze = sp.symbols("lam_ zeta_")
    Ap = lam * Tsym**2
    Zp = ze * Tsym**4
    Bp = (lam**2 - ze) / gam * Tsym**4                 # from Z = A^2 - gamma*B
    ok16a = sp.expand(Ap**2 - gam * Bp - Zp) == 0
    kb0 = reduced_generators(0, sym, "T2")["kbox"]
    tau = sp.Symbol("tau_")
    kb0s = kb0.xreplace({Asym: Ap, Bsym: Bp, Rmsym: sp.Integer(1),
                         sym["d2"]: tau * Tsym**4 / gam})
    # y = -1 is t = 0, and there Q = q (Rm = 1 at n = 0), so Q -> q(-1)
    at_m1 = sp.expand(sp.together(kb0s.xreplace({Tsym: sp.Integer(0),
                                                 Qsym: Q_QUARTIC.subs(y, -1)})))
    ok16 = ok16a and sp.expand(at_m1 + (2 * C_GENUINE / gam) * Q_QUARTIC.subs(y, -1)) == 0
    ck("S16  n=0 on T2: A=lam*t^2, B=((lam^2-z)/g)*t^4, u=tau*t^4 -> kbox at y=-1 "
       "is -mu*q(-1) != 0", ok16,
       "Z = A^2 - g*B check: %s ; kbox|_{t=0} = %s ; mu*q(-1) = %s"
       % (ok16a, at_m1, sp.expand((2 * C_GENUINE / gam) * Q_QUARTIC.subs(y, -1))))

    # ---- S17  n=1 on T2: v=0, u = kappa*Rm^4, kbox at y=-1 is g*kappa*Rm(-1)^6
    kap = sp.Symbol("kappa_")
    kb1 = reduced_generators(1, sym, "T2")["kbox"]
    kb1s = kb1.xreplace({Asym: lam * Tsym**3, Bsym: sp.Integer(0),
                         sym["d2"]: kap * Rmsym**4 / gam})
    at_m1_1 = sp.expand(kb1s.xreplace({Tsym: sp.Integer(0)}))
    ok17 = sp.expand(at_m1_1 - gam * kap * Rmsym**6) == 0
    ck("S17  n=1 on T2: v=0 and u=kappa*Rm^4 (kappa = g^5/(6 lam^3) != 0) -> "
       "kbox at y=-1 is g*kappa*Rm(-1)^6 != 0", ok17,
       "kbox|_{t=0} = %s   [Rm(-1) != 0 because q(-1) = %s != 0]"
       % (at_m1_1, Q_QUARTIC.subs(y, -1)))

    # ---- S18  a9_b1000_T2 TOP STRATUM: the complete forced family, and its death
    #  A = lam*t^3, v = 0 (so dm3 = 0), u = kappa*Rm^4, and then C, d0 follow.
    #  g3 = 0 pins kappa = gamma^5/(6*lam^3).  g1 = g2 = g3 = 0 EXACTLY; only
    #  kbox fails, and it fails at y = -1.
    lamS, kapS = sp.symbols("lam_ kappa_")
    RmS = Rmsym
    Aq = lamS * Tsym**3
    d2q = kapS * RmS**4 / gam
    Cq = -Aq * (kapS * RmS**4) / gam
    d0q = -lamS**2 * kapS * Tsym**6 * RmS**2 / gam**3
    redn1 = reduced_generators(1, sym, "T2")
    forced = {Asym: Aq, Bsym: sp.Integer(0), Csym: Cq,
              sym["d2"]: d2q, sym["d0"]: d0q}
    g1v = sp.expand(sp.together(redn1["g1"].xreplace(forced)))
    g2v = sp.expand(sp.together(redn1["g2"].xreplace(forced)))
    g3v = sp.factor(sp.expand(sp.together(redn1["g3"].xreplace(forced))))
    kbv = sp.expand(sp.together(redn1["kbox"].xreplace(forced)))
    kap_star = gam**5 / (6 * lamS**3)
    g3_at_kap = sp.expand(sp.together(g3v.xreplace({kapS: kap_star})))
    kb_at_m1 = sp.expand(sp.together(
        kbv.xreplace({kapS: kap_star, Tsym: sp.Integer(0)})))
    ok18 = (g1v == 0 and g2v == 0 and g3_at_kap == 0
            and sp.expand(kb_at_m1 - gam**6 * RmS**6 / (6 * lamS**3)) == 0)
    ck("S18  a9_b1000_T2 top stratum: the forced family satisfies g1=g2=g3=0 "
       "exactly and dies ONLY on kbox, at y=-1", ok18,
       "g1=%s g2=%s ; g3 = %s -> 0 iff kappa = g^5/(6 lam^3) ; then "
       "kbox|_{y=-1} = %s != 0 (Rm(-1) != 0, lam != 0, gamma != 0). "
       "deg d2 = 4 and deg d0 = 8 are FORCED: the family lives only at the TOP "
       "stratum, and it is empty there."
       % (g1v, g2v, g3v, kb_at_m1))

    # =================================================================
    #  S19-S24: the coordinator's mid-lane results, adjudicated, and the
    #  reduction re-expressed in the (A, Sbar) coordinates they imply.
    # =================================================================
    Sbar = sp.Symbol("Sbar_")

    # ---- S19  dm4 is NOT a free spare: T = -R*(S/e + d2) - d1*e/2 kills G1
    Tval = -sym["dm2"] * (sym["dm3"] / sym["dm1"] + sym["d2"]) \
        - sym["d1"] * sym["dm1"] / 2
    g1_after_T = sp.simplify(sp.together(G["G1"].xreplace({sym["dm4"]: Tval})))
    # and it is the SAME elimination this file already does through g1
    Cmine = -(Asym * (gam * sym["d2"] + vsym)
              + sp.Rational(1, 2) * gam**2 * sym["d1"] * Rmsym) / gam
    a9 = 9
    theirs = sp.expand(sp.together(Tval.xreplace(
        {sym["dm1"]: gam * Tsym**a9 * Rmsym, sym["dm2"]: Tsym**a9 * Asym,
         sym["dm3"]: Tsym**a9 * Rmsym * vsym})))
    ck("S19  dm4 is determined: T = -R*(S/e + d2) - d1*e/2 annihilates G1, and "
       "equals this file's own g1-elimination of C",
       g1_after_T == 0 and sp.expand(Tsym**a9 * Cmine - theirs) == 0,
       "G1|_{dm4=T} = %s ; (this file's C) - (coordinator's T)/t^a = %s "
       "-> the two eliminations are the SAME map, reached independently"
       % (g1_after_T, sp.expand(Tsym**a9 * Cmine - theirs)))

    # ---- S20  e | dm3 falls out of this lane's own marked-root step
    #  B = Rm*v (S6.1) means dm3 = t^a*Rm*v = e*(v/gamma), i.e. Sbar = v/gamma.
    ck("S20  e | dm3 (the coordinator's `e | S`) is already implied here: "
       "B = Rm*v gives dm3 = e*(v/gamma)",
       sp.simplify((Tsym**a9 * Rmsym * vsym) / (gam * Tsym**a9 * Rmsym) - vsym / gam) == 0,
       "Sbar = v/gamma.  Derived here from kbox + g1 at the marked roots "
       "(S6.1), i.e. by a DIFFERENT mechanism than the Sylvester-resultant "
       "integral-dependence route -- two lanes, one object.")

    # ---- S21  the boxed identity in the (A, Sbar) coordinates
    ok21, det21 = True, []
    for n in range(5):
        kb = reduced_generators(n, sym, "T2")["kbox"]
        lhs = sp.expand(kb.xreplace({Bsym: gam * Rmsym * Sbar}))
        want = sp.expand(3 * Asym**2 + gam**2 * Rmsym**2 * (sym["d2"] + 3 * Sbar)
                         - (2 * C_GENUINE / gam) * Tsym**(3 * n) * Qsym)
        r = sp.expand(lhs - want)
        ok21 &= (r == 0)
        det21.append("n=%d resid=%s" % (n, r))
    ck("S21  boxed identity in (A,Sbar): 3A^2 + g^2*Rm^2*(d2 + 3*Sbar) "
       "= (2c/g)*t^(3n)*Q", ok21, " | ".join(det21))

    # ---- S22  the CORRECTED spare count: 45 -> n+8, not 45-3a
    ok22, rows22 = True, []
    for n in range(5):
        a = DEG_E_SUB2 - n
        nA = (CAP["dm2"] - a) + 1                    # deg A = n+2
        nS = (CAP["dm3"] - CAP["dm1"]) + 1           # deg Sbar = 14-10 = 4
        nC = 0                                       # dm4 determined
        tot = nA + nS + nC
        good = (tot == n + 8)
        ok22 &= good
        rows22.append("%s n=%d a=%d: A:%d + Sbar:%d + dm4:0 = %d (n+8=%d) "
                      "[t^a-only count was %d]"
                      % (FAMILY[n], n, a, nA, nS, tot, n + 8, 45 - 3 * a))
    ck("S22  true free-spare count after e|S and the dm4 elimination: "
       "45 -> n+8 = 8/9/10/11/12, superseding 45-3a = 15/18/21/24/27",
       ok22, " | ".join(rows22))

    # ---- S23  t^a | dm2,dm3,dm4 on BOTH branches -- stronger than the record
    def _rows_ok(a, rho, s, tau, d0o, d2o, d1o, rows=("G1", "G2", "G3", "K")):
        def mt(vals):
            vals = [v for v in vals if v < INF]
            if not vals:
                return True
            m = min(vals)
            return vals.count(m) >= 2
        R = {"G1": [d1o + 2 * a, d2o + a + rho, a + tau, rho + s],
             "G2": [d0o + 2 * a, d2o + 2 * rho, rho + tau, 2 * s],
             "G3": [d0o + a + rho, d1o + 2 * rho, 3 * a, s + tau],
             "K":  [30, d2o + 3 * a, 2 * a + s, a + 2 * rho]}
        return all(mt(R[r]) for r in rows)

    def _scan(a, branch, want_ce, rows=("G1", "G2", "G3", "K")):
        rng = lambda c: list(range(c + 1)) + [INF]
        d1s = [INF] if branch == "T2" else rng(CAP["d1"])
        hits = []
        for rho in rng(CAP["dm2"]):
            for s in rng(CAP["dm3"]):
                for tau in rng(CAP["dm4"]):
                    if (min(rho, s, tau) < a) != want_ce:
                        continue
                    for d0o in rng(CAP["d0"]):
                        for d2o in rng(CAP["d2"]):
                            for d1o in d1s:
                                if _rows_ok(a, rho, s, tau, d0o, d2o, d1o, rows):
                                    hits.append((rho, s, tau, d0o, d2o, d1o))
                                    if len(hits) > 2:
                                        return hits
        return hits

    ok23, det23 = True, []
    for branch in ("T2", "T1"):
        for a in (6, 7, 8, 9, 10):
            ce = _scan(a, branch, True)
            adm = _scan(a, branch, False)
            good = (not ce) and bool(adm)          # refuted AND not vacuous
            ok23 &= good
            det23.append("%s a=%d: ce=%d adm>0=%s" % (branch, a, len(ce), bool(adm)))
    abl = {d: len(_scan(9, "T1", True, tuple(r for r in ("G1", "G2", "G3", "K")
                                             if r != d)))
           for d in ("G1", "G2", "G3", "K")}
    # G1, G3 and K must each be load-bearing; G2 turns out to be REDUNDANT,
    # which is exactly what the hand proof in SPINE.md sec.8.1 shows (it uses
    # only K, G1, G3).  Assert the three, and record the fourth as redundant.
    ok23 &= all(abl[d] > 0 for d in ("G1", "G3", "K")) and abl["G2"] == 0
    ck("S23  t^a | dm2,dm3,dm4 on BOTH branches, a=6..10, from the rows "
       "G1,G3,K (+G2, redundant) -- admissibility control and row ablation",
       ok23,
       "%s || ablation at a=9 T1 (counterexamples per dropped row): %s -> "
       "G1,G3,K are each load-bearing; G2 is REDUNDANT, matching the hand "
       "proof which uses only those three. Admissible non-counterexample "
       "configurations exist throughout, so the enumeration is not vacuously "
       "restrictive." % (", ".join(det23), abl))

    # ---- S24  on T2, R = dm2 is forced to lambda*(y+1)^12 for every n
    ok24 = all((DEG_E_SUB2 - n) + (n + 2) == 12 for n in range(5))
    ck("S24  on T2, A = lam*t^(n+2) gives dm2 = t^a*A = lam*(y+1)^12 for EVERY "
       "n -- the flagged inference `R = c*(y+1)^rho` is a theorem here, rho = 12",
       ok24, " | ".join("n=%d: a+(n+2) = %d + %d = 12"
                        % (n, DEG_E_SUB2 - n, n + 2) for n in range(5)))

    if verbose:
        print("\n".join(out))
    return npass, ntot


VERDICT = """
VERDICT TABLE  (n = deg Rm = 10 - a)

  n  family      branch T2 (d1 = 0)                    branch T1 (d1 != 0)
  -  ----------  ------------------------------------  --------------------------
  0  a10_b0000   EMPTY  (S15 -> A=lam*t^2, S16)         OPEN   <-- the residue
  1  a9_b1000    EMPTY  (S10/S11/S13/S14; also S15/S17) EMPTY  (S11+S13+S14)
  2  a8_b1100    EMPTY  (S10/S11/S12)                   EMPTY  (S11+S12)
  3  a7_b1110    EMPTY  (S10)                           EMPTY  (S10)
  4  a6_b1111    EMPTY  (S10)                           EMPTY  (S10)

  The T1 column is unconditional ONLY because S23 upgrades  t^a | dm2,dm3,dm4
  from T2-only (the recorded status, GSYSTEM_CELL.md sec.7.9) to BOTH branches.
  That upgrade is the single contested premise here and is flagged for
  adjudication in SPINE.md sec.8 -- it has a hand proof, two independent
  machine checks, an admissibility control and a row ablation, but it does
  contradict a recorded status and should be second-partied before the T1
  verdicts are entered anywhere.

  The argument never touches the stratum coordinates (deg d2, deg sigma) except
  through the certified CAP deg d2 <= 4 -- i.e. it is evaluated AT the top
  stratum and holds a fortiori below it.  d0 is eliminated; no sigma ansatz is
  used; no Groebner basis is computed; nothing is mod p.

  Free-spare count after the dm4 elimination and e | dm3:  45 -> n+8
  (8/9/10/11/12), superseding the brief's 45-3a = 15/18/21/24/27 (S22).
  On T2, dm2 = lambda*(y+1)^12 is FORCED for every n (S24).
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args()
    npass, ntot = run(verbose=not a.quiet)
    if a.quiet:
        if npass != ntot:
            print("spine: %d/%d checks FAILED" % (ntot - npass, ntot))
            return 1
        print("spine: %d/%d checks pass" % (npass, ntot))
        return 0
    print("\n%d/%d checks pass" % (npass, ntot))
    if npass == ntot:
        print(VERDICT)
    return 0 if npass == ntot else 1


if __name__ == "__main__":
    sys.exit(main())
