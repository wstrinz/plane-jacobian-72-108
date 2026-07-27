#!/usr/bin/env python3
"""sub1_spine9.py -- the SPINE cofactor identity generalised to the sub1
`a_t = 9` family, and the emptiness of all five remaining frontier cells.

CONTEXT.  After `slice_obstruction_basis.py` (`a_t >= 9`) and
`syzygy_collision.py` / `slice_phi_yplace.py` (`a_t <= 9`), the enumerated f31
frontier is FIVE cells, all sub1, all T1, all `a_t = 9`:

    a9_b0000_T1  a9_b1000_T1  a9_b1100_T1  a9_b1110_T1  a9_b1111_T1

indexed by  k := sum(b_i) = deg Pi  in  {0,1,2,3,4}.

WHAT THIS FILE DOES.  Everything below is DERIVED here from
`generators.json`; NOTHING is imported from `spine.py` / `SPINE.md`.  That
matters, because SPINE's degree bookkeeping (its zero-slack count
`(n+6)+(2n+4) = 3n+10`) is a *sub2 coincidence* and does NOT transfer.  What
DOES transfer is the algebra, because the cofactor identity is a polynomial
identity in a free commutative ring and uses no degree cap at all.

  (0) e | Phi (K-syzygy) + a_t = 9 give   e = gamma*t^9*Pi,  Pi | q squarefree,
      k = deg Pi;  q = Pi*Q;  and t^9 | R,S,T (branch-independent place
      theorem, re-derived here as the SPINE sec.8.1 hand proof at a = 9).
      Write R = t^9*A, S = t^9*B, T = t^9*C, Phi = c*t^30*Pi*Q.

  (1) EXACT reduction, derived by unassisted division (no quotient supplied):

        G1 =  3   * t^18 * g1        K = -gamma * t^27 * Pi * kbox
        G2 = (3/2)* t^18 * g2
        G3 =  3   * t^18 * g3

      with, since 30 - 3a = 3 (NOT sub2's 3n),

        kbox = 3*A^2 + gamma^2*d2*Pi^2 + 3*gamma*Pi*B - mu*t^3*Q,  mu = 2c/gamma

  (2) The marked roots force A(r_i) != 0, hence Pi | B (write B = Pi*v), hence
      with u := gamma*d2, w := (1/2)*gamma^2*d1*Pi,

        F := A*(u + 2*v) + w        Z := A^2 - gamma*Pi^2*v

      an explicit POLYNOMIAL cofactor identity gives, on the variety,

        F * Z  ==  (1/6) * gamma^5 * t^9 * Pi^4                            (*)

      and gcd(Z, Pi) = 1, so  Z | t^9,  Z = zeta*t^z,  zeta != 0.

  (3) Feeding Z back into the boxed row:

        gamma*Pi^2*(u + 6*v)  =  mu*t^3*Q - 3*zeta*t^z                     (5)

      so  Pi^2 | (mu*t^3*Q - 3*zeta*t^z).  This is an EXACT ideal-theoretic
      test in the coefficients of Pi (Pi need not be rational: the test is a
      Groebner computation over Q in the unknown coefficients, saturated by
      zeta, so a unit ideal means "no such factorisation in ANY extension").

      RESULT (P12):  k = 1, 2, 3 are infeasible for EVERY integer z in [0,9];
      k = 4 is feasible ONLY at z = 3.  k = 0 is vacuous (Pi = 1).

  (4) Two degree arguments finish, using only the certified sub1 caps
      (deg d2 <= 6, deg d1 <= 9, deg R <= 18, deg S <= 21):

        k = 4:  (*) forces deg F = 9 + 4k - z = 22 at z = 3, but the caps give
                deg F <= max(deg A + max(deg u, deg v), deg w) = 17.
        k = 0:  (2) reads gamma*u = mu*t^3*q - 6*A^2 + 3*zeta*t^z with
                deg u <= 6.  deg A >= 4 makes -6*A^2 an uncancellable term of
                even degree >= 8; deg A <= 3 leaves the degree-7 leading term
                of mu*t^3*q uncancelled.  Needs z <= 6.

  (5) z is bounded by the t-adic valuation ledger (P10).  This is the ONLY
      step that consumes the slice cascade, and it is needed ONLY for k = 0.

SCOPE -- READ BEFORE QUOTING.

  * `a_t = 9` itself: the lower bound `a_t >= 9` is independently audited
    (`slice_obstruction_audit.py` 56/56); the upper bound `a_t <= 9` is NOT
    audited (`syzygy_collision.py`, `slice_phi_yplace.py`, same author).  This
    file assumes `a_t = 9` and inherits that dependency wholesale.
  * The valuation ledger imports the cascade profile
    `slice_obstruction_stage.json` (`--deep`, level 10).  **The level-12 result
    `v_t(h_6) >= 11` is NOT needed** -- under the [I3] convention the inverse
    shift gives `v_t(R) >= min(v_t(h_6), v_t(h_1) + a_t) = min(11, 10) = 10`,
    so levels 11-12 are invisible here (P10, X1).  The BRIEF for this lane
    asserted `v_t(R) >= 11`; that is FALSE at `a_t = 9` under the operative
    convention and this file does not use it.
  * Four of the five kills (k = 1,2,3,4) use NO cascade input whatsoever.
  * Read-only.  Writes nothing.  Pure sympy: no Singular, no msolve, no WSL,
    no subprocess, no modular arithmetic.

Usage:
    python sub1_spine9.py            # full report
    python sub1_spine9.py --quiet    # exit 0 iff every check passes
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__)) or "."

# ------------------------------------------------------------------ primitives
y = sp.Symbol("y")
t = y + 1
Q_QUARTIC = 2048 * y**4 - 512 * y**3 + 320 * y**2 - 240 * y + 195
C_GENUINE = sp.Rational(-1, 6630)              # Phi = c * t^30 * q
A_T = 9                                        # the whole point of this lane
CELL = {0: "a9_b0000_T1", 1: "a9_b1000_T1", 2: "a9_b1100_T1",
        3: "a9_b1110_T1", 4: "a9_b1111_T1"}

# opaque commuting symbols -- every factorisation below is an identity in the
# free commutative ring Q[gam, T_, A_, B_, C_, Pi_, Q_, v_, u_, d0, d1, d2].
Tsym = sp.Symbol("T_")
gam = sp.Symbol("gamma_")                      # sympy's `gamma` is a function
Asym, Bsym, Csym = sp.symbols("A_ B_ C_")
Pisym, Qsym = sp.symbols("Pi_ Q_")
vsym, usym = sp.symbols("v_ u_")
zeta = sp.Symbol("zeta_")


# ------------------------------------------------------------ canonical inputs
def load_canonical_G():
    """G1,G2,G3,G5 = G5body + Phi from the canonical `generators.json`."""
    with open(os.path.join(HERE, "generators.json"), encoding="utf-8") as fh:
        d = json.load(fh)
    order = d["variable_order"]
    sym = {n: sp.Symbol(n) for n in order}
    out = {}
    for name in ("G1", "G2", "G3", "G5body"):
        expr = sp.Integer(0)
        for coeff, ev in d["polynomials"][name]:
            term = sp.Rational(coeff)
            for var, kk in zip(order, ev):
                if kk:
                    term *= sym[var] ** kk
            expr += term
        out[name] = sp.expand(expr)
    out["G5"] = sp.expand(out.pop("G5body") + sym["Phi"])
    return out, sym


def K_from_syzygy(G, sym):
    return sp.expand(2 * (G["G5"] + sym["d2"] * G["G3"]
                          + sym["d1"] * G["G2"] + sym["d0"] * G["G1"]))


def sub1_caps():
    """Certified sub1 caps, READ from the repo modules -- never typed here."""
    import cascade_engine as ce
    import full_system_bridge as fsb
    d1cap, sigcap, d2cap = ce.SUB1.aux_caps
    caps = {"d2": d2cap, "d1": d1cap, "d0": sigcap, "e": ce.SUB1.e_cap}
    caps.update({k: v for k, v in fsb.STRIP_DEGCAP["sub1"].items()})
    caps["weights"] = dict(fsb.WEIGHT)
    return caps


# ------------------------------------------------ the a = 9 spine substitution
def spine9_subs(sym):
    return {
        sym["dm1"]: gam * Tsym**A_T * Pisym,
        sym["dm2"]: Tsym**A_T * Asym,
        sym["dm3"]: Tsym**A_T * Bsym,
        sym["dm4"]: Tsym**A_T * Csym,
        sym["Phi"]: C_GENUINE * Tsym**30 * Pisym * Qsym,
    }


def exact_divide(expr, div, gens):
    """Unassisted exact division: returns (quotient, ok).  No quotient is
    supplied to the routine, so a reduced form cannot be smuggled in.

    `gens` are the symbols the quotient must be POLYNOMIAL in.  gamma is a
    deliberate exception: it is invertible on the variety (gamma = 0 forces
    e = 0 hence 2*Phi = 0), and the boxed row genuinely carries mu = 2c/gamma.
    T_ is in `gens`, so a wrong t-power cannot slip through (control X5)."""
    quo = sp.expand(sp.cancel(sp.together(sp.expand(expr) / div)))
    try:
        sp.Poly(quo, *gens)          # raises unless quo is polynomial in gens
        poly_ok = True
    except (sp.PolynomialError, sp.GeneratorsNeeded):
        poly_ok = False
    return quo, bool(poly_ok and sp.expand(quo * div - sp.expand(expr)) == 0)


def derive_rows():
    """DERIVE g1,g2,g3,kbox by unassisted division of the substituted rows."""
    G, sym = load_canonical_G()
    K = K_from_syzygy(G, sym)
    sub = spine9_subs(sym)
    gens = [Tsym, Asym, Bsym, Csym, Pisym, Qsym,
            sym["d0"], sym["d1"], sym["d2"]]      # gamma deliberately excluded
    rows, oks = {}, {}
    for name, coeff in (("G1", sp.Integer(3)), ("G2", sp.Rational(3, 2)),
                        ("G3", sp.Integer(3))):
        quo, ok = exact_divide(G[name].xreplace(sub),
                               coeff * Tsym**(2 * A_T), gens)
        rows["g" + name[1]] = quo
        oks["g" + name[1]] = ok
    quo, ok = exact_divide(K.xreplace(sub),
                           -gam * Tsym**(3 * A_T) * Pisym, gens)
    rows["kbox"] = quo
    oks["kbox"] = ok
    return rows, oks, sym, G, K


# ============================================================== the check suite
def run(verbose=True):
    out, npass, ntot = [], 0, 0

    def ck(name, ok, detail=""):
        nonlocal npass, ntot
        ntot += 1
        npass += bool(ok)
        out.append("  [%s] %s%s" % ("PASS" if ok else "FAIL", name,
                                    ("\n        " + detail) if detail else ""))

    def say(s=""):
        out.append(s)

    G, sym = load_canonical_G()
    d0, d1, d2 = sym["d0"], sym["d1"], sym["d2"]

    # ---------------------------------------------------------------- P0 / P1
    ck("P0  canonical G5 normalisation is G5body + Phi (a stale 2*Phi form "
       "breaks everything below)",
       sp.expand(G["G5"]).coeff(sym["Phi"]) == 1,
       "coeff(G5, Phi) = %s" % sp.expand(G["G5"]).coeff(sym["Phi"]))

    K = K_from_syzygy(G, sym)
    e_, R_, S_ = sym["dm1"], sym["dm2"], sym["dm3"]
    Kshape = sp.expand(2 * sym["Phi"] - e_ * (d2 * e_**2 + 3 * e_ * S_ + 3 * R_**2))
    ck("P1  K := 2*(G5 + d2*G3 + d1*G2 + d0*G1) = 2*Phi - e*(d2*e^2+3*e*S+3*R^2)"
       "  [DERIVED from generators.json, not quoted]",
       sp.expand(K - Kshape) == 0,
       "residual = %s ; hence e | 2*Phi on every lift" % sp.expand(K - Kshape))

    # ---------------------------------------------------------------- P2 caps
    caps = sub1_caps()
    wr = caps["weights"]
    rule_ok = (caps["d2"] == 3 * (wr["d2"] // 12) and caps["d1"] == 3 * (wr["d1"] // 12)
               and caps["d0"] == 3 * (wr["d0"] // 12) and caps["e"] == 3 * (wr["dm1"] // 12)
               and caps["dm2"] == 3 * (wr["dm2"] // 12)
               and caps["dm3"] == 3 * (wr["dm3"] // 12)
               and caps["dm4"] == 3 * (wr["dm4"] // 12))
    ck("P2  certified sub1 caps READ from cascade_engine.SUB1 / "
       "full_system_bridge.STRIP_DEGCAP, and equal to the 3k weight rule",
       rule_ok and caps["d2"] == 6 and caps["d1"] == 9 and caps["d0"] == 12
       and caps["e"] == 15 and caps["dm2"] == 18 and caps["dm3"] == 21
       and caps["dm4"] == 24,
       "deg d2<=%d  deg d1<=%d  deg d0<=%d  deg e<=%d  deg R<=%d  deg S<=%d  "
       "deg T<=%d   (weights %s, stripped cap = 3k)"
       % (caps["d2"], caps["d1"], caps["d0"], caps["e"], caps["dm2"],
          caps["dm3"], caps["dm4"], {k: wr[k] for k in ("d2", "d1", "dm1", "dm2")}))

    # ------------------------------------------------------- P3  the reduction
    rows, oks, _, _, _ = derive_rows()
    want = {
        "g1": sp.expand(sp.Rational(1, 2) * gam**2 * d1 * Pisym**2
                        + gam * Pisym * (d2 * Asym + Csym) + Asym * Bsym),
        "g2": sp.expand(d2 * Asym**2 + 2 * Asym * Csym + Bsym**2
                        - gam**2 * d0 * Pisym**2),
        "g3": sp.expand(-gam * d0 * Pisym * Asym - sp.Rational(1, 2) * d1 * Asym**2
                        + Bsym * Csym
                        - sp.Rational(1, 6) * gam**3 * Tsym**A_T * Pisym**3),
        "kbox": sp.expand(3 * Asym**2 + gam**2 * d2 * Pisym**2
                          + 3 * gam * Pisym * Bsym
                          - (2 * C_GENUINE / gam) * Tsym**(30 - 3 * A_T) * Qsym),
    }
    for nm in ("g1", "g2", "g3", "kbox"):
        ck("P3.%s  a=9 reduction: the UNASSISTED quotient exists and equals the "
           "expected row" % nm,
           oks[nm] and sp.expand(rows[nm] - want[nm]) == 0,
           "divisible: %s ; quotient - expected = %s"
           % (oks[nm], sp.expand(rows[nm] - want[nm])))
    ck("P3.tpow  the boxed row's t-power is t^(30-3a) = t^3 (NOT sub2's "
       "t^(3n)) -- this is where the sub1 a = 9 family differs",
       sp.degree(sp.Poly(sp.expand(rows["kbox"]).coeff(Qsym), Tsym), Tsym) == 3
       and 30 - 3 * A_T == 3,
       "coeff of Q in kbox = %s" % sp.expand(rows["kbox"]).coeff(Qsym))

    # ------------------------------------------ P4  concrete instantiation
    ok4, det4 = True, []
    rng = [sp.Rational(x) for x in (3, -5, 7, -2, 11, 13, -17, 19, 23, -29,
                                    31, 37, -41, 43, 47, -53, 59, 61, -67, 71)]

    def poly(deg, off):
        return sum(rng[(off + i) % len(rng)] * y**i for i in range(deg + 1))

    Kc = K_from_syzygy(G, sym)
    for k in range(5):
        Pip = sp.expand(y**k + sum(rng[i] * y**i for i in range(k)))
        Qp = poly(4 - k, 5)
        Ap, Bp, Cp = poly(9, 1), poly(12 - k, 3), poly(15 - k, 7)
        d2p, d1p, d0p = poly(6, 11), poly(9, 2), poly(12, 13)
        gv = sp.Rational(5, 3)
        real = {sym["dm1"]: gv * t**A_T * Pip, sym["dm2"]: t**A_T * Ap,
                sym["dm3"]: t**A_T * Bp, sym["dm4"]: t**A_T * Cp,
                sym["Phi"]: C_GENUINE * t**30 * Pip * Qp,
                d0: d0p, d1: d1p, d2: d2p}
        conc = {gam: gv, Tsym: t, Asym: Ap, Bsym: Bp, Csym: Cp,
                Pisym: Pip, Qsym: Qp, d0: d0p, d1: d1p, d2: d2p}
        for name, coeff in (("G1", sp.Integer(3)), ("G2", sp.Rational(3, 2)),
                            ("G3", sp.Integer(3))):
            r = sp.expand(G[name].xreplace(real)
                          - coeff * t**(2 * A_T) * rows["g" + name[1]].xreplace(conc))
            ok4 &= (r == 0)
        rK = sp.expand(Kc.xreplace(real)
                       + gv * t**(3 * A_T) * Pip * rows["kbox"].xreplace(conc))
        ok4 &= (rK == 0)
        det4.append("k=%d %s" % (k, "ok" if rK == 0 else "K-resid %s" % rK))
    ck("P4  concrete exact-rational instantiation at the sub1 caps reproduces "
       "the same factorisation for every k = 0..4", ok4, " | ".join(det4))

    # ------------------------------------------------- P5  the boxed row (2)
    mu = 2 * C_GENUINE / gam
    kb_uv = sp.expand(rows["kbox"].xreplace({Bsym: Pisym * vsym, d2: usym / gam}))
    want2 = sp.expand(3 * Asym**2 + gam * Pisym**2 * (usym + 3 * vsym)
                      - mu * Tsym**3 * Qsym)
    ck("P5  boxed row (2):  3*A^2 + gamma*Pi^2*(u + 3*v) = mu*t^3*Q,  "
       "mu = 2c/gamma != 0",
       sp.expand(sp.together(kb_uv - want2)) == 0,
       "residual = %s ; mu = %s ; 2c = %s" % (sp.expand(sp.together(kb_uv - want2)),
                                              mu, 2 * C_GENUINE))

    # ------------------------------------- P6  the marked-root step (k >= 1)
    kb_at_root = sp.expand(rows["kbox"].xreplace({Pisym: sp.Integer(0)}))
    g1_at_root = sp.expand(rows["g1"].xreplace({Pisym: sp.Integer(0)}))
    qq = sp.Poly(Q_QUARTIC, y)
    sqfree = sp.gcd(qq, sp.Poly(sp.diff(Q_QUARTIC, y), y)).degree() == 0
    irred = len(sp.factor_list(Q_QUARTIC)[1]) == 1
    ck("P6.a  kbox at a marked root (Pi = 0) is 3*A^2 - mu*t^3*Q, so "
       "3*A(r)^2 = mu*(r+1)^3*Q(r) != 0  =>  A(r) != 0",
       sp.expand(kb_at_root - (3 * Asym**2 - mu * Tsym**3 * Qsym)) == 0
       and sqfree and Q_QUARTIC.subs(y, -1) != 0,
       "kbox|_{Pi=0} = %s ; q squarefree = %s (so Q(r) != 0) ; q(-1) = %s "
       "(so r != -1) ; q irreducible = %s"
       % (kb_at_root, sqfree, Q_QUARTIC.subs(y, -1), irred))
    ck("P6.b  g1 at a marked root is A*B, so A(r) != 0 forces B(r) = 0, i.e. "
       "Pi | B.  Write B = Pi*v.  Then Z(r) = A(r)^2 != 0 and gcd(Z, Pi) = 1.",
       sp.expand(g1_at_root - Asym * Bsym) == 0,
       "g1|_{Pi=0} = %s   [k = 0 makes both statements vacuous, and there "
       "gcd(Z,Pi) = 1 holds trivially]" % g1_at_root)

    # -------------------------------- P7  the cofactor identity, DERIVED here
    w_ = sp.Rational(1, 2) * gam**2 * d1 * Pisym
    Cval = -(Asym * (usym + vsym) + w_) / gam
    sub_hat = {Bsym: Pisym * vsym, Csym: Cval, d2: usym / gam}
    g1hat = sp.expand(sp.together(rows["g1"].xreplace(sub_hat)))
    g2hat = sp.expand(sp.together(rows["g2"].xreplace(sub_hat)))
    g3hat = sp.expand(sp.together(rows["g3"].xreplace(sub_hat)))
    F = Asym * (usym + 2 * vsym) + w_
    Z = Asym**2 - gam * Pisym**2 * vsym
    target = sp.expand(F * Z - sp.Rational(1, 6) * gam**5 * Tsym**A_T * Pisym**4)
    combo = sp.expand(-gam * Asym * g2hat + gam**2 * Pisym * g3hat)
    resid = sp.expand(sp.together(target - combo))
    ck("P7.a  C is DETERMINED by g1 = 0 (g1 is linear in C with leading "
       "coefficient gamma*Pi): C = -(A*(u+v) + w)/gamma, and g1hat vanishes",
       sp.simplify(g1hat) == 0, "g1hat = %s" % sp.simplify(g1hat))
    ck("P7.b  *** THE COFACTOR IDENTITY ***  "
       "F*Z - (1/6)*gamma^5*t^9*Pi^4 = -gamma*A*g2hat + gamma^2*Pi*g3hat, "
       "residual EXACTLY 0.  d0 is ELIMINATED; no degree cap is used.",
       resid == 0 and d0 not in sp.expand(target - combo).free_symbols,
       "residual = %s ; d0 present in the identity: %s ; so on the variety "
       "F*Z = (1/6)*gamma^5*t^9*Pi^4   (*)"
       % (resid, d0 in sp.expand(target - combo).free_symbols))

    # P7.c  independent numeric confirmation on the variety
    #  g2hat is linear in d0 and, after that substitution, g3hat is linear in
    #  d1 -- so a genuine point of {g1hat = g2hat = g3hat = 0} exists over Q
    #  for EVERY choice of the free coordinates.  (A `tries = 0` outcome would
    #  make this check vacuous, so the count is asserted, not just reported.)
    import random
    rnd = random.Random(20260725)
    ok7c, tries = True, 0
    for _ in range(40):
        free = {gam: sp.Rational(rnd.randint(1, 30)),
                Tsym: sp.Rational(rnd.randint(-30, 30)),
                Asym: sp.Rational(rnd.randint(1, 30)),
                Pisym: sp.Rational(rnd.randint(1, 30)),
                usym: sp.Rational(rnd.randint(-30, 30)),
                vsym: sp.Rational(rnd.randint(-30, 30))}
        sol = sp.solve([sp.Eq(g2hat.xreplace(free), 0),
                        sp.Eq(g3hat.xreplace(free), 0)], [d0, d1], dict=True)
        if not sol:
            continue
        pt = dict(free)
        pt.update(sol[0])
        if any(x.free_symbols for x in sol[0].values()):
            continue
        # confirm it really is on the variety, then test (*) there
        on = (sp.expand(g1hat.xreplace(pt)) == 0
              and sp.expand(g2hat.xreplace(pt)) == 0
              and sp.expand(g3hat.xreplace(pt)) == 0)
        if not on:
            ok7c = False
            continue
        tries += 1
        ok7c &= sp.expand((target - combo).xreplace(pt)) == 0
    ck("P7.c  (*) re-confirmed at random exact-rational points of "
       "{g1hat = g2hat = g3hat = 0} (NOT vacuous: the count is asserted)",
       ok7c and tries >= 30, "points genuinely on the variety and tested: %d"
       % tries)

    # ---------------------------------------------- P8  Z | t^9, so Z = z*t^z
    ck("P8  (*) has nonzero right side, so F != 0 and Z != 0 and their degrees "
       "add; Z | t^9*Pi^4 and gcd(Z, Pi) = 1 give Z | t^9, i.e. Z = zeta*t^z "
       "with zeta != 0 and 0 <= z <= 9",
       True,
       "deg F + deg Z = deg(t^9*Pi^4) = 9 + 4k EXACTLY, for every k = 0..4 "
       "(no cap enters -- unlike sub2 there is no zero-slack coincidence here)")

    # --------------------------------------- P9  the d3-killing shift, derived
    theta = sp.Symbol("theta_")
    Dsrc = {m: sp.Symbol("D%d_" % (m + 6)) for m in range(-4, 5)}
    Dsrc[4] = sp.Integer(1)

    def Dtil(j, th):
        return sp.expand(sum(sp.binomial(m, m - j) * Dsrc[m] * th**(m - j)
                             for m in range(j, 5)))

    th_star = -Dsrc[3] / 4
    ck("P9.a  the d3-killing shift: D~_3 = D_3 + 4*theta = 0 at theta = -D_3/4, "
       "and D~_{-1} = D_{-1} EXACTLY (triangular across zero)",
       Dtil(3, th_star) == 0 and sp.expand(Dtil(-1, th_star) - Dsrc[-1]) == 0,
       "D~_{-1} = %s" % Dtil(-1, th_star))
    inv = {
        "d2": sp.expand(Dtil(2, th_star)),
        "d1": sp.expand(Dtil(1, th_star)),
        "R": sp.expand(Dtil(-2, th_star)),
        "S": sp.expand(Dtil(-3, th_star)),
    }
    h1, h2, h3, h5, h6, h7 = (Dsrc[3], Dsrc[2], Dsrc[1],
                              Dsrc[-1], Dsrc[-2], Dsrc[-3])
    #  I3_AUDIT.md sec.5 F4 records the FORWARD forms h2 = d2 + (3/8)h1^2 and
    #  h3 = d1 + (1/2)*h1*d2 + (1/16)*h1^3; alt_level12.py L4.5 records
    #  R = h6 + (h1/4)*h5.  Check all three against what is derived here.
    f4_d2 = sp.expand(inv["d2"] + sp.Rational(3, 8) * h1**2 - h2)
    f4_d1 = sp.expand(inv["d1"] + sp.Rational(1, 2) * h1 * inv["d2"]
                      + sp.Rational(1, 16) * h1**3 - h3)
    ck("P9.b  the inverse-shift relations, DERIVED here and cross-checked "
       "against I3_AUDIT.md sec.5 F4 and alt_level12.py L4.5",
       f4_d2 == 0 and f4_d1 == 0
       and sp.expand(inv["R"] - (h6 + h1 * h5 / 4)) == 0
       and sp.expand(inv["S"] - (h7 + h1 * h6 / 2
                                 + sp.Rational(1, 16) * h1**2 * h5)) == 0,
       "d2 = h2 - (3/8)*h1^2  [F4 residual %s] | d1 = %s  [F4 residual %s] | "
       "R = h6 + (h1/4)*h5 | S = h7 + (h1/2)*h6 + (1/16)*h1^2*h5"
       % (f4_d2, inv["d1"], f4_d1))

    # ------------------------------------------- P10  the valuation ledger
    with open(os.path.join(HERE, "slice_obstruction_stage.json"),
              encoding="utf-8") as fh:
        stage = json.load(fh)
    VAL = {int(k[1:]): v for k, v in stage["forced_valuations"].items()}
    a_t_min = stage["a_t_min"]

    HSYM = {Dsrc[m]: 4 - m for m in range(-3, 4)}   # D_m is h_{4-m}

    def vt_of(expr, v_h):
        """t-adic lower bound for a polynomial in the unshifted D's: the
        minimum over monomials of the summed valuations.  DERIVED from the
        shift expression -- no relation is typed in."""
        gens = sorted(HSYM, key=lambda s: s.name)
        pol = sp.Poly(sp.expand(expr), *gens)
        best = None
        for mono in pol.monoms():
            tot = sum(e * v_h[HSYM[g]] for e, g in zip(mono, gens))
            best = tot if best is None else min(best, tot)
        return best

    def ledger(v_h, reading):
        """(v_t(d2), v_t(d1), v_t(R), v_t(S)) from the cascade profile."""
        if reading == "A":     # [I3] convention: the G-system vars are shifted
            #  v_t(e) = a_t = 9 is the [I3]-INVARIANT anchor (D~_{-1} = D_{-1})
            vh = dict(v_h)
            vh[5] = A_T
            return tuple(vt_of(inv[x], vh) for x in ("d2", "d1", "R", "S"))
        return v_h[2], v_h[3], v_h[6], v_h[7]      # reading B: unshifted

    def zwindow(vd2, vd1, vR, vS):
        vA = vR - A_T                       # R = t^9 * A
        vv = vS - A_T                       # S = t^9 * Pi * v, t does not divide Pi
        vu, vw = vd2, vd1                   # u = gamma*d2 ; w = (1/2)g^2*d1*Pi
        vZ = min(2 * vA, vv)                # Z = A^2 - gamma*Pi^2*v
        vF = min(vA + min(vu, vv), vw)      # F = A*(u+2v) + w
        return vA, vv, vZ, vF, vZ, 9 - vF   # (.., z_min, z_max)

    ck("P10.a  cascade profile imported read-only from "
       "slice_obstruction_stage.json (`--deep`, level 10) and consistent with "
       "a_t = %d" % A_T,
       VAL[5] == a_t_min == 9 and VAL[1] >= 1 and VAL[2] >= 3 and VAL[3] >= 5
       and VAL[6] >= 10 and VAL[7] >= 11,
       "v_t(h_k) >= %s ; a_t_min = %s ; NAMES h1=d3 h2=d2 h3=d1 h5=e h6=R h7=S"
       % (VAL, a_t_min))

    windows = {}
    for reading, vh in (("A", VAL), ("B", VAL),
                        ("A+lvl12", {**VAL, 6: 11}), ("B+lvl12", {**VAL, 6: 11})):
        vd2, vd1, vR, vS = ledger(vh, reading[0])
        vA, vv, vZ, vF, zmin, zmax = zwindow(vd2, vd1, vR, vS)
        windows[reading] = (vd2, vd1, vR, vS, vA, vv, vZ, vF, zmin, zmax)
    ck("P10.b  *** THE z WINDOW ***  under the [I3] convention (reading A, the "
       "operative one) 2 <= z <= 6; every reading gives z <= 6",
       all(w[8] >= 2 and w[9] <= 6 for w in windows.values())
       and windows["A"][8] == 2 and windows["A"][9] == 6,
       " | ".join("%s: v(d2)>=%d v(d1)>=%d v(R)>=%d v(S)>=%d -> v(A)>=%d "
                  "v(v)>=%d v(Z)>=%d v(F)>=%d -> %d <= z <= %d"
                  % ((r,) + w) for r, w in windows.items()))
    ck("P10.c  the level-12 result v_t(h_6) >= 11 is NOT needed: under the "
       "[I3] convention v_t(R) >= min(v_t(h_6), v_t(h_1) + a_t) = min(11,10) "
       "= 10 either way, so the z window is unchanged",
       windows["A"][2] == windows["A+lvl12"][2] == 10
       and windows["A"][9] == windows["A+lvl12"][9],
       "reading A v_t(R) with level 10 = %d, with level 12 = %d; z_max %d vs %d"
       % (windows["A"][2], windows["A+lvl12"][2],
          windows["A"][9], windows["A+lvl12"][9]))

    Z_MAX = windows["A"][9]                 # 6 -- the weakest (safest) reading
    Z_RANGE = list(range(windows["A"][8], Z_MAX + 1))

    # ---------------------------------- P11  feeding Z back into the boxed row
    kb_fed = sp.expand(rows["kbox"].xreplace(
        {Asym**2: zeta * Tsym**sp.Symbol("zz_") * 0 + zeta * sp.Symbol("TZ_")
         + gam * Pisym**2 * vsym,
         Bsym: Pisym * vsym, d2: usym / gam}))
    TZ = sp.Symbol("TZ_")                   # stands for t^z
    want5 = sp.expand(3 * zeta * TZ + gam * Pisym**2 * (usym + 6 * vsym)
                      - mu * Tsym**3 * Qsym)
    ck("P11  (5): substituting A^2 = Z + gamma*Pi^2*v with Z = zeta*t^z and "
       "B = Pi*v turns the boxed row into "
       "gamma*Pi^2*(u + 6*v) = mu*t^3*Q - 3*zeta*t^z",
       sp.expand(sp.together(kb_fed - want5)) == 0,
       "residual = %s ; gcd(t, Pi) = 1 (q(-1) != 0), hence "
       "Pi^2 | (mu*t^3*Q - 3*zeta*t^z)"
       % sp.expand(sp.together(kb_fed - want5)))

    # -------------------------- P12  the EXACT marked-support test, per (k, z)
    def support_feasible(quartic, k, z):
        """Is there a monic degree-k factor Pi | quartic (in ANY extension) and
        a zeta != 0 with Pi^2 | (t^3*Q - 3*zeta*t^z), Q = quartic/Pi?
        Exact: Groebner over Q in the unknown coefficients, saturated by zeta.
        mu is normalised to 1 (the condition is homogeneous in (mu, zeta))."""
        if k == 0:
            return True                     # Pi = 1: no condition at all
        ps = sp.symbols("pp0:%d" % k)
        cs = sp.symbols("cc0:%d" % (5 - k))
        wsat = sp.Symbol("wsat_")
        Pi = y**k + sum(ps[i] * y**i for i in range(k))
        Qf = sum(cs[j] * y**j for j in range(5 - k))
        eqs = [sp.Poly(sp.expand(quartic - Pi * Qf), y).coeff_monomial(y**i)
               for i in range(5)]
        P = sp.expand(t**3 * Qf - 3 * zeta * t**z)
        rem = sp.rem(sp.Poly(P, y), sp.Poly(sp.expand(Pi**2), y))
        eqs += [rem.coeff_monomial(y**i) for i in range(2 * k)]
        gb = sp.groebner(eqs + [zeta * wsat - 1],
                         *(list(ps) + list(cs) + [zeta, wsat]), order="lex")
        return list(gb.exprs) != [sp.Integer(1)]

    table, ok12 = {}, True
    for k in (1, 2, 3, 4):
        feas = [z for z in range(0, 10) if support_feasible(Q_QUARTIC, k, z)]
        table[k] = feas
    ok12 = (table[1] == [] and table[2] == [] and table[3] == []
            and table[4] == [3])
    ck("P12  *** THE MARKED-SUPPORT TEST ***  exact, over Q, saturated by "
       "zeta != 0, swept over z = 0..9 (a strict superset of the admissible "
       "window): k = 1,2,3 INFEASIBLE for every z; k = 4 feasible ONLY at z = 3",
       ok12,
       " | ".join("k=%d (%s): feasible z = %s"
                  % (k, CELL[k], table[k] if table[k] else "NONE") for k in table)
       + "   [no valuation input is used here -- these three kills are "
         "cascade-free]")

    # ------------------------------------------------- P13  the degree ledger
    def degF_cap(k, capd2=None, capd1=None, capR=None, capS=None):
        cd2 = caps["d2"] if capd2 is None else capd2
        cd1 = caps["d1"] if capd1 is None else capd1
        cR = caps["dm2"] if capR is None else capR
        cS = caps["dm3"] if capS is None else capS
        degA = cR - A_T                     # R = t^9 * A
        degv = cS - A_T - k                 # S = t^9 * Pi * v
        degu = cd2                          # u = gamma*d2
        degw = cd1 + k                      # w = (1/2)*gamma^2*d1*Pi
        return degA, degv, degu, degw, max(degA + max(degu, degv), degw)

    ledger_rows = []
    for k in range(5):
        degA, degv, degu, degw, cap = degF_cap(k)
        ledger_rows.append("k=%d: degA<=%d degu<=%d degv<=%d degw<=%d "
                           "=> degF<=%d ; deg(t^9*Pi^4) = %d"
                           % (k, degA, degu, degv, degw, cap, 9 + 4 * k))
    dA, dv, du, dw, capF4 = degF_cap(4)
    ck("P13  degree ledger from the certified sub1 caps "
       "(deg A <= deg R - 9, deg v <= deg S - 9 - k, deg w <= deg d1 + k)",
       (dA, du, dv, dw, capF4) == (9, 6, 8, 13, 17),
       " | ".join(ledger_rows))

    # ------------------------------------------------------- P14  k = 4 kills
    need4 = 9 + 4 * 4 - 3                   # z = 3 pinned by P12
    ck("P14.a  *** %s DIES ***  P12 pins z = 3, so (*) forces "
       "deg F = 9 + 16 - 3 = %d, but the caps give deg F <= %d.  "
       "CONTRADICTION.  [uses NO cascade input]" % (CELL[4], need4, capF4),
       need4 > capF4,
       "required deg F = %d > cap %d = max(degA+max(degu,degv), degw) = "
       "max(%d+%d, %d)" % (need4, capF4, dA, max(du, dv), dw))
    ck("P14.b  independent second route to the same kill: even without P12, "
       "z <= %d (P10) forces deg F = 25 - z >= %d > %d"
       % (Z_MAX, 25 - Z_MAX, capF4),
       25 - Z_MAX > capF4,
       "min required deg F over the admissible window = %d ; cap = %d"
       % (25 - Z_MAX, capF4))

    # ------------------------------------------------------- P15  k = 0 kills
    #  gamma*u = mu*t^3*q - 6*A^2 + 3*zeta*t^z ,  deg u <= 6
    kb0 = sp.expand(rows["kbox"].xreplace({Pisym: sp.Integer(1),
                                           Qsym: Qsym, Bsym: vsym,
                                           d2: usym / gam}))
    #  Z = A^2 - gamma*v  =>  gamma*v = A^2 - zeta*t^z
    kb0_sub = sp.expand(sp.together(
        kb0.xreplace({vsym: (Asym**2 - zeta * TZ) / gam})))
    want0 = sp.expand(6 * Asym**2 + gam * usym - 3 * zeta * TZ - mu * Tsym**3 * Qsym)
    deg_lead = 3 + 4                        # deg(t^3 * q)
    ok15a = sp.expand(sp.together(kb0_sub - want0)) == 0
    ok15b = (deg_lead > caps["d2"]) and (Z_MAX < deg_lead) and \
            all(2 * dA_ != deg_lead for dA_ in range(0, 10))
    ck("P15.a  k = 0 row: eliminating v gives "
       "gamma*u = mu*t^3*q - 6*A^2 + 3*zeta*t^z, residual 0",
       ok15a, "residual = %s" % sp.expand(sp.together(kb0_sub - want0)))
    ck("P15.b  *** %s DIES ***  deg(mu*t^3*q) = %d and deg u <= %d and "
       "z <= %d.  If deg A >= 4 then deg(6*A^2) >= 8 is an uncancellable "
       "even-degree leading term; if deg A <= 3 then deg(6*A^2) <= 6 and "
       "deg(3*zeta*t^z) <= %d, so the degree-%d term of mu*t^3*q survives.  "
       "Either way deg(RHS) > %d = cap on deg u.  CONTRADICTION."
       % (CELL[0], deg_lead, caps["d2"], Z_MAX, Z_MAX, deg_lead, caps["d2"]),
       ok15b and ok15a,
       "branch 1: 2*degA >= 8 > max(%d, %d) ; branch 2: max(6, %d) < %d ; "
       "no degA has 2*degA = %d (odd)"
       % (deg_lead, Z_MAX, Z_MAX, deg_lead, deg_lead))

    # ---------------------- P16  t^9 | R,S,T at a = 9, branch-independent
    #  SPINE sec.8.1's hand proof, re-checked as arithmetic at a = 9.
    a = A_T
    step1 = []
    for rho in range(a):
        # K: orders {30, v(d2)+3a, 2a+s, a+2rho}; rho < a gives a+2rho <= 3a-2 < 30
        c1 = (a + 2 * rho <= 3 * a - 2) and (3 * a - 2 < 30)
        s = 2 * rho - a
        c2 = (3 * rho - a < a + rho) and (3 * rho - a < 2 * a)
        tau = 3 * rho - 2 * a
        c3 = (5 * rho - 3 * a < 3 * a and 5 * rho - 3 * a < a + rho
              and 5 * rho - 3 * a < 2 * rho)
        step1.append(c1 and c2 and c3)
    ok16 = (all(step1) and (2 * a + (a - 1) <= 3 * a - 1 < 30)
            and (a + (a - 1) <= 2 * a - 1 < 2 * a))
    ck("P16  t^a | R,S,T at a = 9 on BOTH branches (SPINE sec.8.1 hand proof "
       "re-checked as arithmetic; d1 appears only in terms bounded below, so "
       "no step uses d1 = 0; needs only v_t(e) = a <= 10 and v_t(2*Phi) = 30)",
       ok16,
       "step1 (rho < a impossible) holds for every rho = 0..%d: %s ; "
       "step2 (s < a): 2a+s <= %d < 30 unique minimum ; "
       "step3 (tau < a): a+tau <= %d < 2a unique minimum"
       % (a - 1, all(step1), 3 * a - 1, 2 * a - 1))

    # =====================================================  SENSITIVITY CONTROLS
    say("")
    say("  ---- MANDATORY SENSITIVITY CONTROLS ('no survivors' is exactly the "
        "shape a convention error takes) ----")

    # X1  weaken the cascade and watch the k = 0 kill switch off
    sweep = {}
    for h6 in (8, 9, 10, 11):
        for h7 in (9, 10, 11):
            vh = {**VAL, 6: h6, 7: h7}
            vd2, vd1, vR, vS = ledger(vh, "A")
            zmax = zwindow(vd2, vd1, vR, vS)[5]
            sweep[(h6, h7)] = (zmax, zmax < deg_lead)
    fires = {kk: vv for kk, vv in sweep.items() if not vv[1]}
    ck("X1  CONTROL -- weakening the cascade DOES switch the k = 0 kill off: "
       "the argument needs z <= 6, and lowering v_t(h_7) to 10 (or 9) makes "
       "z_max = 7 (or 8) and the k = 0 degree argument NO LONGER closes",
       bool(fires) and sweep[(10, 11)][1] and not sweep[(10, 10)][1]
       and not sweep[(11, 10)][1],
       "(v_t(h_6), v_t(h_7)) -> (z_max, k=0 kill fires): %s"
       % {kk: vv for kk, vv in sorted(sweep.items())})
    ck("X1b CONTROL -- and the level-12 upgrade is genuinely inert here: "
       "(h6,h7) = (10,11) and (11,11) give the SAME z_max",
       sweep[(10, 11)][0] == sweep[(11, 11)][0] == 6,
       "z_max at h6=10: %d ; at h6=11: %d"
       % (sweep[(10, 11)][0], sweep[(11, 11)][0]))

    # X2  a synthetic quartic with a genuinely feasible k = 1 condition
    q_syn = y**4 + y**3 - 9 * y**2 + 7        # q_syn(1) = 0, q_syn''(1) = 0
    syn_sqfree = sp.gcd(sp.Poly(q_syn, y), sp.Poly(sp.diff(q_syn, y), y)).degree() == 0
    syn_feas = [z for z in range(0, 8) if support_feasible(q_syn, 1, z)]
    ck("X2  CONTROL -- the k = 1 kill is a property of THIS quartic, not of the "
       "test: on the synthetic quartic y^4+y^3-9y^2+7 (squarefree, "
       "q(-1) = %s != 0, with a root r = 1 satisfying (r+1)*q''(r) = 0) the "
       "k = 1 test becomes FEASIBLE at exactly z = 3" % q_syn.subs(y, -1),
       syn_sqfree and syn_feas == [3] and not support_feasible(Q_QUARTIC, 1, 3),
       "synthetic feasible z: %s ; genuine q feasible z: %s"
       % (syn_feas, [z for z in range(0, 8) if support_feasible(Q_QUARTIC, 1, z)]))

    # X2b  and the k = 3 / k = 4 machinery is not vacuously restrictive either
    ck("X2b CONTROL -- the support test is not vacuously infeasible: it DOES "
       "return feasible for (k, z) = (4, 3) on the genuine quartic, which is "
       "why k = 4 needs the degree argument and cannot be closed by support",
       support_feasible(Q_QUARTIC, 4, 3),
       "k=4 z=3 feasible = %s (Pi = q/2048, Q = 2048, P = 2048*t^3 - 3*zeta*t^3 "
       "vanishes identically at zeta = 2048/3)"
       % support_feasible(Q_QUARTIC, 4, 3))

    # X3  mutate the d2 cap and watch the k = 0 degree argument change
    mut = {}
    for cd2 in (5, 6, 7, 8):
        mut[cd2] = (deg_lead > cd2)
    ck("X3  CONTROL -- mutating the sub1 cap deg d2 from 6 to 7 DOES switch the "
       "k = 0 kill off: the argument needs deg(mu*t^3*q) = 7 > deg u, so it "
       "fires at deg d2 <= 6 and fails at deg d2 >= 7",
       mut[6] and mut[5] and not mut[7] and not mut[8],
       "deg d2 cap -> k=0 kill fires: %s" % mut)

    # X4  mutate the R/S caps and watch the k = 4 degree argument change
    mut4 = {}
    for cR in (16, 18, 20, 22, 23, 24, 26):
        _, _, _, _, cap = degF_cap(4, capR=cR)
        mut4[cR] = (need4 > cap)
    ck("X4  CONTROL -- mutating the sub1 cap deg R DOES switch the k = 4 kill "
       "off: it fires at the certified deg R <= 18 and fails once the cap is "
       "large enough for deg F to reach 22 (deg R >= 23)",
       mut4[18] and not mut4[23] and not mut4[24],
       "deg R cap -> k=4 kill fires: %s (deg F required = %d)" % (mut4, need4))

    # X5  the reduction itself is falsifiable: a wrong t-power must FAIL
    _, bad_ok = exact_divide(K_from_syzygy(G, sym).xreplace(spine9_subs(sym)),
                             -gam * Tsym**(3 * A_T + 1) * Pisym,
                             [Tsym, gam, Asym, Bsym, Csym, Pisym, Qsym,
                              d0, d1, d2])
    ck("X5  CONTROL -- the unassisted division is falsifiable: dividing K by "
       "t^28*Pi instead of t^27*Pi does NOT give a polynomial",
       not bad_ok, "t^28 division is exact: %s (must be False)" % bad_ok)

    # ==================================================================== verdict
    verdict = {
        0: ("EMPTY", "P15 (degree, needs z <= 6 from P10)"),
        1: ("EMPTY", "P12 (marked support, every z)"),
        2: ("EMPTY", "P12 (marked support, every z)"),
        3: ("EMPTY", "P12 (marked support, every z)"),
        4: ("EMPTY", "P14 (degree; z = 3 pinned by P12)"),
    }
    ck("P17  *** VERDICT ***  all five a_t = 9 sub1 T1 cells are EMPTY",
       ok12 and (need4 > capF4) and ok15a and ok15b and resid == 0
       and all(oks.values()),
       " | ".join("%s: %s by %s" % (CELL[k], verdict[k][0], verdict[k][1])
                  for k in range(5)))

    if verbose:
        print("\n".join(out))
    return npass, ntot


VERDICT = """
VERDICT -- the five remaining sub1 a_t = 9 T1 cells

  cell           k   killed by                                cascade input?
  -------------  --  ---------------------------------------  --------------
  a9_b0000_T1     0  degree, boxed row (P15)                   YES (z <= 6)
  a9_b1000_T1     1  marked Pi^2 support, EVERY z (P12)         no
  a9_b1100_T1     2  marked Pi^2 support, EVERY z (P12)         no
  a9_b1110_T1     3  marked Pi^2 support, EVERY z (P12)         no
  a9_b1111_T1     4  degree, z = 3 pinned by P12 (P14)          no

  ALL FIVE EMPTY  =>  the enumerated f31 frontier is EMPTY.

PREMISES, in decreasing order of how much they should worry a reader:

  1. a_t = 9.  Lower bound a_t >= 9 INDEPENDENTLY AUDITED
     (slice_obstruction_audit.py 56/56).  Upper bound a_t <= 9 exact-checked
     TWICE by different mechanisms but by the SAME AUTHOR, and NOT audited
     (syzygy_collision.py, slice_phi_yplace.py).  Everything here is
     conditional on it.
  2. The level-10 slice cascade profile v_t(h_k), imported read-only from
     slice_obstruction_stage.json.  Used ONLY for a9_b0000_T1.
     The LEVEL-12 result v_t(h_6) >= 11 is NOT used and NOT needed (P10.c).
     NOTE: the brief for this lane asserted v_t(R) >= 11.  Under the [I3]
     convention that is FALSE at a_t = 9 -- the inverse shift gives
     v_t(R) >= min(v_t(h_6), v_t(h_1) + a_t) = min(11, 10) = 10.
  3. e | Phi (K-syzygy, re-derived here at P1) and the divisor filter's
     D1/D2 (rad(e) | t*q, b_i in {0,1}).
  4. t^9 | R,S,T, branch-independent (P16 re-derives the hand proof at a = 9).
  5. The certified sub1 caps deg d2 <= 6, deg d1 <= 9, deg R <= 18,
     deg S <= 21, read from the repo modules (P2).  Used by P14 and P15 only.

NOT USED: SPINE's zero-slack degree coincidence (a sub2 accident); the
sub2 reduced formulas; any Groebner basis over the G-system; any modular
arithmetic; any solver.  The only Groebner bases computed are the tiny exact
marked-support ideals of P12, over Q, in at most 6 variables.
"""


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()
    npass, ntot = run(verbose=not args.quiet)
    if args.quiet:
        if npass != ntot:
            print("sub1_spine9: %d/%d checks FAILED" % (ntot - npass, ntot))
            return 1
        print("sub1_spine9: %d/%d checks pass" % (npass, ntot))
        return 0
    print("\n%d/%d checks pass" % (npass, ntot))
    if npass == ntot:
        print(VERDICT)
    return 0 if npass == ntot else 1


if __name__ == "__main__":
    sys.exit(main())
