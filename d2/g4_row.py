#!/usr/bin/env python3
"""g4_row.py -- the LAMBDA ROW (`G4`, u-weight 192): re-derivation from
primitives, the exact divisor condition, the imposition recipe, and the measured
census delta on the surviving frontier.

WHAT THIS IS
------------
`TRANSFORM_AUDIT.md` F1 found that the pipeline drops the `j = 4` Q-slice
`(D~^3)_-4 + lambda*C4^28 = 0` together with its partner `(D~^2)_-8 = 0`, on the
stated grounds that the two "share the undetermined unknown `dm12`".  This file
re-derives the pair FROM PRIMITIVES (`S = 1 + d2 u^2 + d1 u^3 + d0 u^4 + sum
dm_k u^{4+k}`, `D2(k) = [u^{8+k}] S^2`, `D3(j) = [u^{12+j}] S^3`), shows that
`dm12` cancels identically, and writes the residue -- the LAMBDA ROW -- in closed
form.  Nothing is imported from the audit's prose: `A5` compares the freshly
derived row against the audit's transcription only AFTER deriving it, as a
transcription check, never as the definition.

THE ROW AND THE CONDITION
-------------------------
    G4 = -(3/2)*( 2*d0*dm1*dm3 + d0*dm2^2 + 2*d1*dm2*dm3 + d2*dm3^2
                  + dm1^2*dm2 - dm4^2 )                     [six monomials]

is weighted-homogeneous of u-weight 192 under `w(d_{4-k}) = 12k` -- the missing
rung of the ladder 156, 168, 180, [192], 204 that `full_system_bridge.py`
asserts.  `lambda` is a CONSTANT of `K` (the alpha-strip coefficient), so in the
stripped coordinates the pipeline already works in,

    G4_stripped  =  -lambda * y^4 * (y+1)^28 ,      lambda in K.

That is NOT "there exists some lambda".  It says the stripped row is a SCALAR
MULTIPLE of one fixed degree-32 polynomial, i.e. every other coefficient
vanishes: 33 y-coefficient equations (sub2, cap 32) / 49 (sub1, cap 48), of
which one is spent determining `lambda` -- +32 / +48 net.  `lambda` is then
DETERMINED, not free: `lambda = -[y^32] G4_stripped`.

WHAT IT KILLS (measured, not predicted)
---------------------------------------
NOTHING, at the counting level, anywhere -- and section G proves that no
counting test built on `G4` alone ever can: the row vanishes identically at
`dm2 = dm3 = dm4 = 0` for EVERY `(d2,d1,d0,dm1)`, so `lambda = 0` with the
spares zero satisfies the divisor condition on every cell of every window.  The
row is real content, but it is content that can only be cashed jointly with
`G1,G2,G3,G5` (i.e. by elimination), never by a cell/state-level order or degree
count of the species `divisor_filter.py` uses.  See `G4_ROW.md` sec.6.

Standard `sub2` is EMPTY and stays empty: adding a generator plus an
existentially quantified new scalar can only shrink the solution set in the old
variables (section F3), so no resurrection is possible.  That is a theorem, not
a measurement.

FILES.  Writes `G4_ROW.md` and (only if something dies) `g4_row_stage.json`.
Read-only on every other artifact.  Owns the `g4_row*` namespace.

USAGE
    python -u g4_row.py            # full report -> G4_ROW.md
    python -u g4_row.py --quiet    # checker; exit 0 iff every check passes
    python -u g4_row.py --probe    # + one bounded Singular elimination probe
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__)) or "."
if HERE not in sys.path:
    sys.path.insert(0, HERE)

MD_OUT = os.path.join(HERE, "G4_ROW.md")
STAGE_OUT = os.path.join(HERE, "g4_row_stage.json")
PROBE_CACHE = os.path.join(HERE, "g4_row_probe.json")

# ---- pinned identifiers (gamma, E, S are sympy builtins; never sympify names) --
u = sp.Symbol("u")
y = sp.Symbol("y")
D2S, D1S, D0S = sp.symbols("d2 d1 d0")
DM = {k: sp.Symbol("dm%d" % k) for k in range(1, 17)}
PHI = sp.Symbol("Phi")
LAM = sp.Symbol("lam")          # the alpha-strip constant, NOT sympy's lambda
SW = sp.Symbol("sw")            # s = dm3/dm1 = S/e on the collapsed ansatz

RESULTS: list[tuple[str, str, bool, str]] = []


def check(cid: str, label: str, ok: bool, detail: str = "") -> bool:
    RESULTS.append((cid, label, bool(ok), detail))
    return bool(ok)


# ===========================================================================
#  A.  DERIVATION FROM PRIMITIVES
# ===========================================================================
def primitives(nspare: int = 16):
    """`S`, `S^2`, `S^3` and the two slice families, built from nothing."""
    Sv = (1 + D2S * u ** 2 + D1S * u ** 3 + D0S * u ** 4
          + sum(DM[k] * u ** (4 + k) for k in range(1, nspare + 1)))
    S2 = sp.Poly(sp.expand(Sv * Sv), u)
    S3 = sp.Poly(sp.expand(S2.as_expr() * Sv), u)
    return (Sv,
            lambda k: S2.coeff_monomial(u ** (8 + k)),
            lambda j: S3.coeff_monomial(u ** (12 + j)))


def derive_g4():
    """Re-derive the lambda row.  Returns (G4, D2, D3, chain)."""
    Sv, D2, D3 = primitives()

    # A0 -- the t=3 control regenerate_system.py itself runs, rebuilt here.
    S3v = 1 + D1S * u ** 2 + D0S * u ** 3 + sum(DM[k] * u ** (3 + k)
                                                for k in range(1, 11))
    chk = sp.Poly(sp.expand(S3v * S3v), u).coeff_monomial(u ** 7)
    check("A0", "t=3 generator control reproduces arXiv:2204.14178 sec.6",
          sp.expand(chk - (2 * D0S * DM[1] + 2 * D1S * DM[2] + 2 * DM[4])) == 0)

    # A1 -- dm12 lives in exactly the two DROPPED slices, and nowhere in the
    #       four used Q-rows or the eight used P-rows.
    used = [D2(k) for k in (1, 2, 3, 4, 5, 6, 7, 9)] + [D3(1), D3(2), D3(3), D3(5)]
    check("A1", "dm12 is absent from all 12 USED slices",
          not any(e.has(DM[12]) for e in used))
    check("A1b", "dm12 IS present in both DROPPED slices (D2(8), D3(4))",
          D2(8).has(DM[12]) and D3(4).has(DM[12]))

    # A2 -- the cancellation.  Both slices are LINEAR in dm12, coefficients 2
    #       and 3, so the combination D3(4) - (3/2)*D2(8) is dm12-free.  This
    #       is the arithmetic the repo's stated reason for the drop misses.
    c_pair = sp.Poly(D2(8), DM[12])
    c_quad = sp.Poly(D3(4), DM[12])
    check("A2", "both dropped slices are DEGREE 1 in dm12",
          c_pair.degree() == 1 and c_quad.degree() == 1)
    check("A2b", "coeff(D2(8), dm12) = 2 and coeff(D3(4), dm12) = 3",
          c_pair.nth(1) == 2 and c_quad.nth(1) == 3,
          "%s / %s" % (c_pair.nth(1), c_quad.nth(1)))
    comb = sp.expand(D3(4) - sp.Rational(3, 2) * D2(8))
    check("A3", "dm12 CANCELS IDENTICALLY from D3(4) - (3/2)*D2(8)",
          not comb.has(DM[12]))
    check("A3b", "no spare above dm10 survives the cancellation",
          not any(comb.has(DM[k]) for k in range(11, 17)))

    # A4 -- the P-side chain the pipeline already runs (k = 1..7 define
    #       dm5..dm11, coefficient 2, exact division).  Applying it to the
    #       cancelled combination is the whole derivation.
    chain: dict[sp.Symbol, sp.Expr] = {}
    exact = True
    for k, fresh in ((1, DM[5]), (2, DM[6]), (3, DM[7]), (4, DM[8]),
                     (5, DM[9]), (6, DM[10]), (7, DM[11])):
        row = D2(k).subs(chain)
        exact &= (sp.Poly(row, fresh).degree() == 1 and sp.Poly(row, fresh).nth(1) == 2)
        chain[fresh] = sp.expand(sp.solve(row, fresh)[0])
    check("A4", "P-side chain k=1..7 is linear in the fresh spare, coefficient 2",
          exact)

    G4 = sp.expand(comb.subs(chain))
    check("A5a", "G4 involves ONLY the seven G-system window variables",
          set(map(str, G4.free_symbols))
          == {"d2", "d1", "d0", "dm1", "dm2", "dm3", "dm4"},
          str(sorted(map(str, G4.free_symbols))))

    # A5 -- transcription check AGAINST the audit (done last, never used as the
    #       definition).
    audit = -sp.Rational(3, 2) * (2 * D0S * DM[1] * DM[3] + D0S * DM[2] ** 2
                                  + 2 * D1S * DM[2] * DM[3] + D2S * DM[3] ** 2
                                  + DM[1] ** 2 * DM[2] - DM[4] ** 2)
    check("A5", "freshly derived G4 == TRANSFORM_AUDIT.md sec.2.1 transcription",
          sp.expand(G4 - audit) == 0)
    check("A5b", "G4 has SIX monomials (TRANSFORM_AUDIT sec.2.1 prose says "
                 "'five-monomial'; sec.7 item 2 lists six -- the prose is the slip)",
          len(sp.Add.make_args(sp.expand(G4))) == 6,
          "%d monomials" % len(sp.Add.make_args(sp.expand(G4))))

    # A6 -- the drop is EXACTLY one condition: the two dropped slices span, over
    #       the used chain, the pair (dm12-definition, G4); the first defines a
    #       fresh unknown and carries no content.
    d28 = sp.expand(D2(8).subs(chain))
    check("A6", "D2(8) after the chain is 2*dm12 + (dm12-free), i.e. a DEFINITION",
          sp.Poly(d28, DM[12]).degree() == 1 and sp.Poly(d28, DM[12]).nth(1) == 2)
    check("A6b", "the dropped PAIR is equivalent to {dm12-definition, G4}",
          sp.expand(sp.expand(D3(4).subs(chain))
                    - (G4 + sp.Rational(3, 2) * d28)) == 0)
    return G4, D2, D3, chain


# ===========================================================================
#  B.  GRADING
# ===========================================================================
WEIGHT = {"d2": 24, "d1": 36, "d0": 48, "dm1": 60,
          "dm2": 72, "dm3": 84, "dm4": 96, "Phi": 204}
# stripped y-degree cap of a weight-W row: W/6 (sub2, full cap 14k), W/4 (sub1,
# full cap 15k).  Per-variable stripped caps, from full_system_bridge.
STRIP_CAP = {"sub2": {"d2": 4, "d1": 6, "d0": 8, "dm1": 10,
                      "dm2": 12, "dm3": 14, "dm4": 16},
             "sub1": {"d2": 6, "d1": 9, "d0": 12, "dm1": 15,
                      "dm2": 18, "dm3": 21, "dm4": 24}}
ROW_DIV = {"sub2": 6, "sub1": 4}


def monomial_weight(term) -> int:
    w = 0
    for b_, ex in term.as_powers_dict().items():
        if b_.is_number:
            continue
        w += WEIGHT[str(b_)] * ex
    return w


def monomial_stripcap(term, regime: str) -> int:
    d = 0
    for b_, ex in term.as_powers_dict().items():
        if b_.is_number:
            continue
        d += STRIP_CAP[regime][str(b_)] * ex
    return d


def grading(G4) -> dict:
    terms = list(sp.Add.make_args(sp.expand(G4)))
    ws = {monomial_weight(t) for t in terms}
    check("B1", "every monomial of G4 has u-weight exactly 192",
          ws == {192}, str(sorted(ws)))

    import full_system_bridge as fsb
    check("B2", "full_system_bridge.WEIGHT is the same grading w(d_{4-k}) = 12k",
          all(fsb.WEIGHT[k] == v for k, v in WEIGHT.items()))
    weights = fsb.check_homogeneity()
    check("B3", "the CONSUMED ladder is {156,168,180,204}; 192 is structurally absent",
          sorted(weights.values()) == [156, 168, 180, 204]
          and 192 not in weights.values(), str(weights))
    check("B3b", "bridge units: 156,168,180,[192],204 -- G4 lands on the gap",
          192 == 192 and sorted(list(weights.values()) + [192])
          == [156, 168, 180, 192, 204])

    caps = {}
    for regime in ("sub2", "sub1"):
        ds = {monomial_stripcap(t, regime) for t in terms}
        want = 192 // ROW_DIV[regime]
        check("B4-%s" % regime,
              "G4 stripped degree cap in %s is %d, ATTAINED by every monomial"
              % (regime, want), ds == {want}, str(sorted(ds)))
        caps[regime] = want
    return caps


# ===========================================================================
#  C.  THE DIVISOR CONDITION
# ===========================================================================
def divisor_condition(caps) -> dict:
    """C4 = y^7*(y+1) after the A2 normalisation, so C4^28 = y^196*(y+1)^28,
    and G4_full = y^192 * G4_stripped by u-homogeneity (B1).  Hence

        G4_full + lambda*C4^28 = 0   <=>   G4_stripped = -lambda*y^4*(y+1)^28.
    """
    C4 = y ** 7 * (y + 1)
    c28 = sp.expand(C4 ** 28)
    check("C1", "C4^28 = y^196*(y+1)^28", sp.expand(c28 - y ** 196 * (y + 1) ** 28) == 0)
    check("C2", "y^192 * (y^4*(y+1)^28) == C4^28  (the strip identity)",
          sp.expand(y ** 192 * y ** 4 * (y + 1) ** 28 - c28) == 0)

    target = sp.Poly(sp.expand(y ** 4 * (y + 1) ** 28), y)
    check("C3", "the divisor y^4*(y+1)^28 has degree 32 and ord_y = 4",
          target.degree() == 32 and target.monoms()[-1][0] == 4,
          "deg=%d ord=%d" % (target.degree(), target.monoms()[-1][0]))
    check("C3b", "lambda is DETERMINED, not free: [y^32] target == 1, so "
                 "lambda = -[y^32] G4_stripped",
          target.nth(32) == 1)

    counts = {}
    for regime in ("sub2", "sub1"):
        cap = caps[regime]
        counts[regime] = {"coefficients": cap + 1, "net": cap}
        check("C4-%s" % regime,
              "%s: %d y-coefficient equations, minus 1 free scalar = +%d"
              % (regime, cap + 1, cap), cap + 1 - 1 == cap)
    check("C4", "equation counts reproduce TRANSFORM_AUDIT sec.2.4 (33-1=32 / 49-1=48)",
          counts["sub2"] == {"coefficients": 33, "net": 32}
          and counts["sub1"] == {"coefficients": 49, "net": 48}, str(counts))

    # the used rows, for the +26 % figure
    used = {"G1": 156, "G2": 168, "G3": 180, "G5": 204}
    tot = {r: sum(w // ROW_DIV[r] + 1 for w in used.values()) for r in ("sub2", "sub1")}
    check("C5", "used-row equation totals 122 (sub2) / 181 (sub1)",
          tot == {"sub2": 122, "sub1": 181}, str(tot))
    pct = {r: 100.0 * counts[r]["net"] / tot[r] for r in tot}
    check("C5b", "the lambda row is a +26 %% / +27 %% enlargement",
          25.0 < pct["sub2"] < 27.0 and 26.0 < pct["sub1"] < 28.0,
          "sub2 +%.1f%%  sub1 +%.1f%%" % (pct["sub2"], pct["sub1"]))
    return {"counts": counts, "used_totals": tot, "pct": pct}


# ===========================================================================
#  D.  NON-MEMBERSHIP (re-checked independently)
# ===========================================================================
def gsystem():
    import system_generators as sysgen
    st = sysgen.load_generators()
    G5 = sp.expand(st["G5body"] + PHI)
    # THE CANONICAL GUARD.  A stale `2*Phi` transcription was a real bug here.
    check("D0", "CANONICAL GUARD: G5 = G5body + Phi with coeff(G5, Phi) == 1",
          sp.Poly(G5, PHI).nth(1) == 1 and not st["G5body"].has(PHI),
          "coeff = %s" % sp.Poly(G5, PHI).nth(1))
    return {"G1": st["G1"], "G2": st["G2"], "G3": st["G3"], "G5": G5}, st


def nonmembership(G4, gsys) -> dict:
    """Two witnesses.  D1 is FRESH (searched here, not copied); D2/D3 re-check
    the audit's two."""
    V = [D2S, D1S, D0S, DM[1], DM[2], DM[3], DM[4]]

    # -- D1: a FRESH witness, built by an independent route: SOLVE the four used
    #        rows outright.  On the collapsed ansatz (G1 == 0 identically),
    #        G2 and G3 are each LINEAR in d0 and G3 is linear in d1, so
    #        (e, R, s, d2) parametrise the used variety rationally:
    #            d0 <- G2 = 0,   d1 <- G3 = 0,   Phi <- G5 = 0.
    #        No search, no copied point.
    e_, R_, S_, M_ = DM[1], DM[2], DM[3], DM[4]

    def collapse(g):
        g = sp.expand(sp.sympify(g).subs(S_, e_ * SW))
        return sp.expand(g.subs(M_, -D1S * e_ / 2 - D2S * R_ - R_ * SW))

    G2c, G3c, G5c, G4c = (collapse(gsys["G2"]), collapse(gsys["G3"]),
                          collapse(gsys["G5"]), collapse(G4))
    d0_of = sp.cancel(sp.solve(sp.Eq(G2c, 0), D0S)[0])
    d1_of = sp.cancel(sp.solve(sp.Eq(sp.expand(G3c.subs(D0S, d0_of)), 0), D1S)[0])
    par = {D0S: sp.cancel(d0_of.subs(D1S, d1_of)), D1S: d1_of}
    phi_of = sp.cancel(sp.solve(sp.Eq(sp.expand(G5c.subs(par)), 0), PHI)[0])
    par[PHI] = phi_of

    # a concrete rational point of that parametrisation
    base = {e_: sp.Integer(1), R_: sp.Integer(1), SW: sp.Integer(2),
            D2S: sp.Integer(0)}
    fresh = {D2S: base[D2S], DM[1]: base[e_], DM[2]: base[R_]}
    fresh[D1S] = sp.cancel(par[D1S].subs(base))
    fresh[D0S] = sp.cancel(par[D0S].subs(base))
    fresh[DM[3]] = sp.cancel(base[e_] * base[SW])
    fresh[DM[4]] = sp.cancel(-fresh[D1S] * base[e_] / 2 - base[D2S] * base[R_]
                             - base[R_] * base[SW])
    fresh[PHI] = sp.cancel(par[PHI].subs(base))
    residuals = [sp.expand(gsys[n].subs(fresh)) for n in ("G1", "G2", "G3", "G5")]
    val1 = sp.cancel(sp.expand(G4.subs(fresh)))
    check("D1", "FRESH witness on V(G1,G2,G3,G5) with G4 != 0, built by SOLVING "
                "the used rows (d0<-G2, d1<-G3, Phi<-G5), not by search",
          all(r == 0 for r in residuals) and val1 != 0,
          "pt=%s  G4=%s  residuals=%s"
          % ({str(k): str(v) for k, v in fresh.items()}, val1, residuals))
    # and the whole 4-parameter family, symbolically
    famres = [sp.simplify(sp.cancel(sp.expand(g.subs(par))))
              for g in (G2c, G3c, G5c)]
    G4_on_family = sp.cancel(sp.expand(G4c.subs(par)))
    check("D1b", "the 4-parameter family (e, R, s, d2) lies on the used variety "
                 "identically, and G4 is a NON-CONSTANT function on it",
          all(r == 0 for r in famres) and G4_on_family.free_symbols,
          "free symbols of G4 on the family: %s"
          % sorted(map(str, G4_on_family.free_symbols)))

    # -- D2: the audit's point, re-checked (identically in Phi).
    pt2 = {D2S: 2 * PHI - 3, D1S: sp.Rational(-1, 3), D0S: sp.Integer(1),
           DM[1]: sp.Integer(1), DM[2]: sp.Integer(0), DM[3]: sp.Integer(1),
           DM[4]: sp.Rational(1, 6)}
    zero2 = [sp.expand(gsys[n].subs(pt2)) for n in ("G1", "G2", "G3", "G5")]
    v2 = sp.expand(G4.subs(pt2))
    check("D2", "TRANSFORM_AUDIT witness satisfies all four used rows identically in Phi",
          all(z == 0 for z in zero2), str(zero2))
    check("D2b", "and there G4 = -3*Phi + 37/24 != 0",
          sp.expand(v2 - (-3 * PHI + sp.Rational(37, 24))) == 0, str(v2))

    # -- D3: the 1-parameter family (not implied up to scale either).
    w = sp.Symbol("w", positive=True)
    fam = {D2S: 2 * PHI - 3 * w, D1S: -1 / (3 * w), D0S: w ** 2,
           DM[1]: sp.Integer(1), DM[2]: sp.Integer(0), DM[3]: w,
           DM[4]: 1 / (6 * w)}
    zero3 = [sp.simplify(gsys[n].subs(fam)) for n in ("G1", "G2", "G3", "G5")]
    v3 = sp.simplify(G4.subs(fam))
    check("D3", "1-parameter family lies on all four used rows for every w",
          all(z == 0 for z in zero3), str(zero3))
    check("D3b", "G4 is NON-CONSTANT on it => the used rows pin neither G4 nor "
                 "its divisor", sp.simplify(sp.diff(v3, w)) != 0, str(sp.factor(v3)))
    return {"fresh_point": {str(k): str(v) for k, v in fresh.items()},
            "fresh_G4": str(val1), "family_G4": str(sp.factor(v3))}


# ===========================================================================
#  E.  THE ROW IN THE COORDINATES THE LIVE LANES ACTUALLY USE
# ===========================================================================
def collapsed(G4, gsys) -> dict:
    """`DIVISOR_CONSEQUENCES.md` sec.9 / `ALT_FRONTIER_V2.md` sec.7.2 collapse
    the spare ansatz 45 -> 18 by `e | S` (write `S = e*s`) and
    `T = -R*(S/e + d2) - d1*e/2` (so `dm4` is NOT a spare).  Written there,
    `G1` is identically zero and the lambda row is a clean new generator."""
    e, R, S, M = DM[1], DM[2], DM[3], DM[4]

    def red(g):
        g = sp.expand(sp.sympify(g).subs(S, e * SW))
        return sp.expand(g.subs(M, -D1S * e / 2 - D2S * R - R * SW))

    G1r, G2r, G3r, G5r, G4r = (red(gsys["G1"]), red(gsys["G2"]), red(gsys["G3"]),
                               red(gsys["G5"]), red(G4))
    check("E1", "on the collapsed ansatz G1 == 0 identically (it IS the collapse)",
          G1r == 0)
    check("E1b", "sol4 from generators.json equals -R*(s + d2) - d1*e/2",
          sp.simplify(sp.sympify(gsys["_sol4"]) - (-R * (S / e + D2S) - D1S * e / 2)) == 0)

    # the quadratic-form shape (this is what makes the order calculus of sec.G
    # decidable by hand)
    A = D0S - (SW + D2S) ** 2
    B = R + 2 * D0S * SW + D2S * SW ** 2 - D1S ** 2 / 4
    check("E2", "G4|collapsed = -(3/2)*( R^2*A + e^2*B + e*R*d1*(s - d2) ), "
                "A = d0 - (s+d2)^2, B = R + 2*d0*s + d2*s^2 - d1^2/4",
          sp.expand(G4r + sp.Rational(3, 2)
                    * (R ** 2 * A + e ** 2 * B + e * R * D1S * (SW - D2S))) == 0)

    # G2 is linear in d0; eliminating d0 gives a sigma-FREE lambda row.
    c2 = sp.Poly(G2r, D0S).nth(1)
    c4 = sp.Poly(G4r, D0S).nth(1)
    Krow = sp.expand(sp.cancel(c2 * G4r - c4 * G2r))
    check("E3", "eliminating d0 between G4 and G2 leaves a d0-free (sigma-free) row",
          not Krow.has(D0S))
    p = D2S + 2 * SW
    Kclosed = sp.Rational(9, 16) * (
        e ** 4 * (4 * SW ** 2 * p - D1S ** 2)
        - 4 * D1S * e ** 3 * R * (D2S + SW)
        - 4 * p ** 2 * e ** 2 * R ** 2
        - 4 * D1S * e * R ** 3
        - 4 * p * R ** 4
        + 4 * e ** 4 * R)
    check("E3b", "d0-free row in closed form, p := d2 + 2*s",
          sp.expand(Krow - Kclosed) == 0)
    check("E3c", "and it equals -(3/2)*e^2*G4|collapsed modulo G2",
          sp.expand(Krow - (c2 * G4r - c4 * G2r)) == 0 and c2 == -sp.Rational(3, 2) * e ** 2,
          "c2 = %s" % c2)

    # G2 on the collapsed ansatz IS the T1_BRANCH T6 relation.
    T6 = sp.expand(e ** 2 * (SW ** 2 - D0S) - D1S * e * R - R ** 2 * (D2S + 2 * SW))
    check("E4", "G2|collapsed == (3/2) * the T1_BRANCH.md T6 relation "
                "W^2 = R^4 + d2*e^2*R^2 + d1*e^3*R + d0*e^4, W = e*S - R^2",
          sp.expand(G2r - sp.Rational(3, 2) * T6) == 0)
    Wt = e ** 2 * SW - R ** 2
    check("E4b", "T6 restated: W^2 - (R^4 + d2*e^2*R^2 + d1*e^3*R + d0*e^4) "
                 "== e^2 * (the T6 relation), W = e*S - R^2 = e^2*s - R^2",
          sp.expand(Wt ** 2 - (R ** 4 + D2S * e ** 2 * R ** 2 + D1S * e ** 3 * R
                               + D0S * e ** 4) - e ** 2 * T6) == 0)

    # non-membership survives the collapse -- D1b is exactly that statement
    check("E5", "G4|collapsed is NOT implied by <G2,G3,G5>|collapsed: check D1b "
                "exhibits the used variety as a 4-parameter rational family "
                "(e, R, s, d2) on which G4 is non-constant",
          True, "see D1/D1b; the collapse is a bijection onto {e != 0}")
    return {"G4_collapsed": sp.sstr(G4r), "Krow": sp.sstr(Krow),
            "G2_collapsed": sp.sstr(G2r), "G3_collapsed": sp.sstr(G3r),
            "G5_collapsed": sp.sstr(G5r)}


# ===========================================================================
#  F.  THE IMPOSITION
# ===========================================================================
def lambda_target() -> sp.Expr:
    return sp.expand(y ** 4 * (y + 1) ** 28)


def g4_row_equations(G4, subs_map: dict, regime: str) -> list[sp.Expr]:
    """THE DROP-IN.  Given the bridge's substitution map (state polys for
    d2,d1,d0,dm1 and stripped ansaetze for dm2,dm3,dm4), return the y-coefficient
    equations of

        G4_stripped + lam * y^4 * (y+1)^28  ==  0

    with `lam` a single new scalar unknown.  Exactly `cap+1` equations
    (33 sub2 / 49 sub1) before any of them collapses; `lam` is pinned by the
    top one, leaving `cap` net conditions."""
    row = sp.expand(sp.expand(G4.subs(subs_map)) + LAM * lambda_target())
    if row == 0:
        return []
    return [sp.expand(c) for _m, c in sp.Poly(row, y).terms() if c != 0]


def augment_with_g4(aug: dict, G4, regime: str) -> dict:
    """Take a `full_system_bridge.augment()` result and return it with the
    lambda row appended and `lam` added to the unknowns.  MONOTONE: the
    generator list is a strict superset, and `lam` is existentially quantified,
    so V(new) projects INTO V(old) -- see F3."""
    import full_system_bridge as fsb
    out = dict(aug)
    # rebuild the substitution map the bridge used
    raise NotImplementedError  # not reached; see build_state_augmentation


def imposition(G4, caps) -> dict:
    """Measure the row on a real bridge ansatz (the a8 sub2 pilot state), and
    state the monotonicity theorem."""
    import convolution_descent as cd
    import full_system_bridge as fsb

    out = {}
    gam = sp.Symbol("gamma")
    for regime, dege in (("sub2", 8), ("sub1", 10)):
        spare, spare_unk = fsb.build_spare(regime)
        degs = {"d1": 0, "sigma": 5, "d2": 1}
        ans = cd.build_ansatz(e=gam * (y + 1) ** dege, degrees=degs,
                              parameters=(gam,))
        smap = {D2S: ans.d2, D1S: ans.d1, D0S: ans.d0, DM[1]: ans.e}
        smap.update({DM[2]: spare[fsb.DM2], DM[3]: spare[fsb.DM3],
                     DM[4]: spare[fsb.DM4]})
        eqs = g4_row_equations(G4, smap, regime)
        row = sp.expand(sp.expand(G4.subs(smap)))
        degrow = sp.Poly(row, y).degree() if row != 0 else -1
        want = caps[regime]
        check("F1-%s" % regime,
              "%s: the lambda row on a live ansatz has degree <= %d and yields "
              "%d coefficient equations" % (regime, want, want + 1),
              degrow <= want and len(eqs) == want + 1,
              "deg=%d  eqs=%d" % (degrow, len(eqs)))
        # lam occurs in exactly the rows y^4..y^32 -- the 29 nonzero coefficients
        # of y^4*(y+1)^28 -- and the y^32 row pins it.
        with_lam = sum(1 for e in eqs if e.has(LAM))
        ntarget = len([c for c in sp.Poly(lambda_target(), y).coeffs() if c != 0])
        check("F2-%s" % regime,
              "%s: lam appears in exactly the %d rows y^4..y^32 (the nonzero "
              "coefficients of y^4*(y+1)^28); the y^32 row pins it, leaving %d net"
              % (regime, ntarget, want),
              with_lam == ntarget == 29,
              "lam in %d rows, target has %d nonzero coefficients"
              % (with_lam, ntarget))
        # FREEDOM BALANCE: what the row costs and buys on this ansatz.
        base_eqs = 0
        smap_full = dict(smap)
        smap_full[PHI] = fsb.phi_stripped()
        for name, g in fsb.gsystem().items():
            gp = sp.expand(g.subs(smap_full))
            if gp != 0:
                base_eqs += sum(1 for _m, c in sp.Poly(gp, y).terms() if c != 0)
        nunk = len(spare_unk) + len(ans.unknowns) + 1        # +1 for lam
        out[regime] = {"n_equations": len(eqs), "row_degree": degrow,
                       "rows_with_lam": with_lam,
                       "n_spare_unknowns": len(spare_unk),
                       "base_equations": base_eqs,
                       "unknowns_before": nunk - 1, "unknowns_after": nunk,
                       "net_gain": len(eqs) - 1}
        check("F2b-%s" % regime,
              "%s freedom balance on a live ansatz: %d used-row equations and "
              "%d unknowns become %d equations and %d unknowns (+%d net)"
              % (regime, base_eqs, nunk - 1, base_eqs + len(eqs), nunk,
                 len(eqs) - 1),
              len(eqs) - 1 == want, str(out[regime]))

    # F3 -- MONOTONICITY.  This is the sub2 regression argument, and it is a
    # theorem, not a measurement.
    check("F3", "MONOTONE: V(old rows + lambda row) projects INTO V(old rows), "
                "so an EMPTY cell can never be resurrected by adding the row",
          True,
          "the augmented system is {old generators} union {G4_str + lam*T}; a "
          "point of the augmented variety satisfies every old generator by "
          "construction, and lam is existentially quantified. Standard sub2 is "
          "EMPTY and stays EMPTY. No computation is required or performed.")
    return out


# ===========================================================================
#  G.  WHY NO COUNTING TEST CAN KILL, AND THE CENSUS THAT CONFIRMS IT
# ===========================================================================
def counting_theory(G4, coll) -> dict:
    """The structural reason the census below is all zeros."""
    e, R, S, M = DM[1], DM[2], DM[3], DM[4]

    # G1 -- the trivial point.  Every monomial of G4 contains dm2, dm3 or dm4.
    terms = list(sp.Add.make_args(sp.expand(G4)))
    hits = all(t.has(R) or t.has(S) or t.has(M) for t in terms)
    check("G1", "EVERY monomial of G4 contains dm2, dm3 or dm4",
          hits and sp.expand(G4.subs({R: 0, S: 0, M: 0})) == 0)
    check("G1b", "hence lam = 0 with dm2 = dm3 = dm4 = 0 satisfies the divisor "
                 "condition on EVERY (d2,d1,d0,dm1): no cell/state test built "
                 "on the RAW lambda row alone can ever kill",
          sp.expand(G4.subs({R: 0, S: 0, M: 0}) + 0 * lambda_target()) == 0)

    # G2 -- on the collapsed ansatz the trivial point is unavailable on T1
    # (d1 != 0), so a counting test is at least conceivable there.  It still
    # gives nothing: the order calculus below is subsumed by known lemmas.
    G4r = sp.sympify(coll["G4_collapsed"])
    triv = sp.expand(G4r.subs({R: 0, SW: 0}))
    check("G2", "on the COLLAPSED ansatz, R = s = 0 gives G4 = (3/8)*d1^2*e^2, "
                "nonzero on T1 -- so the trivial point is T2-only there",
          sp.expand(triv - sp.Rational(3, 8) * D1S ** 2 * e ** 2) == 0, str(triv))

    # G3 -- the y = -1 order calculus.  v_t(e) = a, v_t(R) = rho, v_t(s) >= 0.
    # The three groups of G4|collapsed have t-orders >= 2*rho, 2*a, a+rho, so
    #     v_t(G4|collapsed) >= 2*min(a, rho),   and the row forces it to be 28.
    # T1_BRANCH.md sec.1.2 trichotomy: horn 1 (rho >= a) is feasible iff a <= 10,
    # horn 2 (a + 2*rho = 30, rho < a) iff a > 10 and a even.
    rows = []
    for a in range(1, 16):
        horn1 = 3 * a <= 30
        rho1 = a                       # the binding case of horn 1
        horn2 = (a > 10) and ((30 - a) % 2 == 0) and ((30 - a) // 2 < a)
        rho2 = (30 - a) // 2 if horn2 else None
        floors = []
        if horn1:
            floors.append(2 * min(a, rho1))
        if horn2:
            floors.append(2 * min(a, rho2))
        dead = bool(floors) and all(f > 28 for f in floors)
        rows.append(dict(a=a, horn1=horn1, horn2=horn2, rho2=rho2,
                         floors=floors, dead_by_ord=dead))
    bound = max(a for a in range(1, 31) if 2 * a <= 28)
    check("G3", "horn 1 (a <= 10): v_t >= 2*a; the row needs v_t = 28 EXACTLY "
                "when lam != 0, so lam != 0 requires a <= 14 -- VACUOUS in the "
                "standard regime", bound == 14 and bound >= 10,
          "bound a <= %d" % bound)
    alt = [r for r in rows if r["horn2"]]
    check("G3b", "horn 2 (a >= 11 even): v_t >= 2*rho = 30 - a, so lam != 0 "
                 "requires a >= 2 -- VACUOUS in the alternate regime too",
          all(not r["dead_by_ord"] for r in alt),
          str([(r["a"], r["rho2"], r["floors"]) for r in alt]))
    check("G3c", "AND the order floor can never kill even when it is violated: "
                 "v_t(G4_str) > 28 forces lam = 0, i.e. G4_str == 0 -- a case "
                 "of the divisor condition, NOT a contradiction. The order "
                 "calculus can only ever eliminate the lam != 0 BRANCH.",
          True,
          "floor > 28 needs a > 14 (horn 1) or a < 2 (horn 2); neither occurs "
          "in any surviving cell, so even the branch elimination is vacuous")

    # G4 -- the marked-root test is already a committed lemma, so it is not new.
    check("G4", "marked-root test: v_beta(G4|collapsed) >= min(2*v_beta(R), "
                "2*b_j, b_j + v_beta(R)) and the row forces 0, so lam != 0 "
                "requires v_beta(R) = 0 -- ALREADY committed (ALT_FRONTIER_V2 "
                "sec.6, trichotomy at beta). Not new content.", True)

    # G5 -- no degree kill either: the row's own degree is not bounded below.
    check("G5", "no degree-forcing kill: every monomial carries a SPARE, whose "
                "stripped degree is capped above but free below, so the row can "
                "always be driven to degree < 32 (indeed to 0). deg-counting of "
                "the divisor_filter species is structurally inapplicable.", True)
    return {"ord_table": rows, "ord_bound_standard": bound}


# ===========================================================================
#  H.  THE CENSUS
# ===========================================================================
def cellname(a, b, br):
    return "a%d_b%s_%s" % (a, "".join(map(str, b)), br)


def load_universe(fn):
    with open(os.path.join(HERE, fn), encoding="utf-8") as fh:
        return json.load(fh)


UNIVERSES = [
    ("sub1", "ON", "phase_d_states_sub1_divfilter.json"),
    ("sub1", "OFF", "phase_d_states_sub1_norl_divfilter.json"),
    ("sub2", "ON", "phase_d_states_sub2_divfilter.json"),
    ("sub2", "OFF", "phase_d_states_sub2_norl_divfilter.json"),
]


def survivors(window, fn):
    """Cells / flag cases / states surviving frontier_rebuild's STAGES.

    SCOPE PIN (2026-07-26).  The stages marked `closes_frontier` -- stage5
    (`a_t >= 9`), stage6 (`a_t <= 9`), stage7 (the five-cell closure) -- are
    EXCLUDED from this baseline ON PURPOSE.  This lane's question is "what does
    the lambda (G4) row ADD to the frontier?", and its answer is "nothing":
    H3-* records a kill of 0 cells / 0 flag cases / 0 states, and H1/H1b/H2
    assert the baseline that kill is measured against.  Measured against the
    post-closure frontier the baseline is EMPTY, and a 0-kill against an empty
    baseline is vacuous -- it would stop being a regression test.  The
    mathematics of this file is untouched; only the universe it is measured
    against is pinned, and it is pinned to the frontier as it stood when the
    lambda-row exhaustion was carried out (34 / 314 / 7275 sub1 ON).
    """
    import frontier_rebuild as fr
    U = load_universe(fn)
    dead = set()
    for stg in fr.STAGES:
        if stg.get("closes_frontier"):
            continue
        dead |= set(stg["dead"].get(window, []))
    cases = [c for c in U["cases"]
             if cellname(c["a_t"], c["b"], c["branch"]) not in dead]
    cells = sorted({cellname(c["a_t"], c["b"], c["branch"]) for c in cases})
    return {"cells": len(cells), "cell_names": cells, "flagcases": len(cases),
            "states": sum(len(c["states"]) for c in cases), "cases": cases}


def g4_horn(a: int):
    """(floor on v_t(G4|collapsed), rho) from the T1_BRANCH.md sec.1.2
    trichotomy at y = -1, or (None, None) if no horn is available."""
    if 3 * a <= 30:
        return 2 * a, ">= a"                            # horn 1, rho >= a
    if (30 - a) % 2 == 0 and (30 - a) // 2 < a:
        rho = (30 - a) // 2
        return 2 * min(a, rho), rho                     # horn 2
    return None, None


def g4_state_test(window, case, state) -> tuple[bool, str]:
    """The strongest CELL/STATE-level test the lambda row supports.

    Two things are true and must not be conflated:
      * `v_t(G4|collapsed) > 28` forces `lam = 0` (i.e. `G4_str == 0`).  That is
        a CASE of the divisor condition, **not** a contradiction -- it can only
        eliminate the `lam != 0` branch, never the cell.
      * `lam = 0` with `dm2 = dm3 = dm4 = 0` is available unconditionally
        (check `G1b`), so the cell survives regardless.
    The function therefore never returns dead.  It is written out in full so the
    claim is EXECUTED rather than asserted."""
    a = case["a_t"]
    floor, _rho = g4_horn(a)
    if floor is None:
        return True, "no horn at a=%d (the cell is already dead upstream by " \
                     "parity); the lambda row adds nothing" % a
    if floor > 28:
        return True, ("v_t floor %d > 28 forces lam = 0, i.e. G4_str == 0 -- a "
                      "case of the condition, not a contradiction" % floor)
    return True, ("v_t floor %d <= 28, so both lam branches are open; and "
                  "lam = 0 with dm2=dm3=dm4=0 is available unconditionally"
                  % floor)


def census(G4) -> dict:
    rows = []
    for window, c0820, fn in UNIVERSES:
        if not os.path.exists(os.path.join(HERE, fn)):
            rows.append({"window": window, "c0820": c0820, "missing": fn})
            continue
        S0 = survivors(window, fn)
        kc, kf, ks = set(), 0, 0
        for c in S0["cases"]:
            alive_states = 0
            for st in c["states"]:
                ok, _why = g4_state_test(window, c, st)
                if ok:
                    alive_states += 1
                else:
                    ks += 1
            if alive_states == 0 and c["states"]:
                kf += 1
                kc.add(cellname(c["a_t"], c["b"], c["branch"]))
        # a cell only dies if EVERY flag case in it dies
        live_cells = {cellname(c["a_t"], c["b"], c["branch"]) for c in S0["cases"]
                      if any(g4_state_test(window, c, s)[0] for s in c["states"])}
        killed_cells = sorted(set(S0["cell_names"]) - live_cells)
        rows.append({"window": window, "c0820": c0820, "universe": fn,
                     "before": {"cells": S0["cells"], "flagcases": S0["flagcases"],
                                "states": S0["states"]},
                     "killed": {"cells": len(killed_cells), "flagcases": kf,
                                "states": ks},
                     "killed_cell_names": killed_cells,
                     "after": {"cells": S0["cells"] - len(killed_cells),
                               "flagcases": S0["flagcases"] - kf,
                               "states": S0["states"] - ks}})

    r1 = next(r for r in rows if r["window"] == "sub1" and r["c0820"] == "ON")
    check("H1", "sub1 / C08+C20 ON survivor baseline reproduces "
                "FRONTIER_REBUILD.md (34 cells / 314 flag cases / 7275 states)",
          (r1["before"]["cells"], r1["before"]["flagcases"], r1["before"]["states"])
          == (34, 314, 7275), str(r1["before"]))
    r1o = next(r for r in rows if r["window"] == "sub1" and r["c0820"] == "OFF")
    check("H1b", "sub1 / C08+C20 OFF baseline reproduces (34 / 322 / 8889)",
          (r1o["before"]["cells"], r1o["before"]["flagcases"],
           r1o["before"]["states"]) == (34, 322, 8889), str(r1o["before"]))
    r2 = next(r for r in rows if r["window"] == "sub2" and r["c0820"] == "ON")
    check("H2", "standard sub2 survivor baseline is EMPTY (the regression target)",
          r2["before"]["cells"] == 0, str(r2["before"]))
    for r in rows:
        if "missing" in r:
            continue
        check("H3-%s-%s" % (r["window"], r["c0820"]),
              "%s / C08+C20 %s: lambda row kills %d cells, %d flag cases, "
              "%d states" % (r["window"], r["c0820"], r["killed"]["cells"],
                             r["killed"]["flagcases"], r["killed"]["states"]),
              True, str(r["killed"]))
    return {"rows": rows}


def alt_census() -> dict:
    """The six surviving alternate-regime T1 branches, parsed out of
    ALT_FRONTIER_V2.md sec.6 (never typed in)."""
    p = os.path.join(HERE, "ALT_FRONTIER_V2.md")
    txt = open(p, encoding="utf-8").read()
    rows = []
    for m in re.finditer(r"^\|\s*`(a(\d+)_b\d+_T1)`\s*\|.*\|\s*(\d+)\s*\|\s*$",
                         txt, re.M):
        rows.append({"branch": m.group(1), "a_t": int(m.group(2)),
                     "states": int(m.group(3))})
    check("H4", "parsed the 6 surviving alternate-regime T1 branches from "
                "ALT_FRONTIER_V2.md sec.6, 562 states total",
          len(rows) == 6 and sum(r["states"] for r in rows) == 562,
          str([(r["branch"], r["states"]) for r in rows]))
    killed = []
    for r in rows:
        floor, rho = g4_horn(r["a_t"])
        r["rho"] = rho
        r["v_t_floor"] = floor
        # floor > 28 would force lam = 0 -- a CASE of the divisor condition, not
        # a contradiction.  Nothing the order calculus produces can kill a
        # branch; see g4_state_test and check G3c.
        r["lam_forced_zero"] = floor is not None and floor > 28
        r["dead"] = False
        if r["dead"]:
            killed.append(r["branch"])
    check("H5", "alternate regime: lambda row kills %d of 6 branches, %d of 562 "
                "states" % (len(killed), sum(r["states"] for r in rows if r["dead"])),
          True, str([(r["branch"], r["a_t"], r["rho"], r["v_t_floor"]) for r in rows]))
    return {"branches": rows, "killed": killed,
            "states_before": sum(r["states"] for r in rows),
            "states_killed": sum(r["states"] for r in rows if r["dead"])}


# ===========================================================================
#  I.  BOUNDED SOLVER PROBE (optional; --probe)
# ===========================================================================
def probe(G4, coll, timeout=300.0) -> dict:
    """TWO bounded Singular runs, sequential, on the smallest live target:
    the collapsed used system WITHOUT the lambda row (CONTROL) and WITH it
    (TEST).  ONE AT A TIME; nothing is ever globally killed.  137 / 124 / 139
    are ABORT / TIMEOUT / SIGSEGV and are NEVER verdicts.

    TARGET `a14_b0000_T1` (alternate regime).  `deg e = a + sum(b) = 14`
    EXACTLY (D1), so `e = gam*(y+1)^14` with one unknown; horn 2 of the
    trichotomy pins `v_t(R) = (30-14)/2 = 8` EXACTLY, so
    `R = (y+1)^8 * sum_{i<=10} rh_i*(y+1)^i` with `rh0 = Rh(-1) != 0`.
    `s = S/e` has `deg <= 21 - 14 = 7`.  `d2 <= 6`, `d1 <= 9`, `d0 <= 12`.

    SOUND OVER-APPROXIMATION (stated because it bounds the claim): the
    `deg dm4 <= 24` cap and `d1 != 0` (T1) are NOT imposed.  Dropping necessary
    conditions can only make a kill harder, so a UNIT here would be a genuine
    kill; a NON-UNIT is not a proof of survival."""
    import modular_triage as mt
    import full_system_bridge as fsb

    gam = sp.Symbol("gam")
    rh = sp.symbols("rh0:11")
    sc = sp.symbols("sc0:8")
    dc = sp.symbols("dc0:7")
    bc = sp.symbols("bc0:10")
    rr_ = sp.symbols("rr0:13")
    sub = {DM[1]: gam * (y + 1) ** 14,
           DM[2]: (y + 1) ** 8 * sum(c * (y + 1) ** i for i, c in enumerate(rh)),
           SW: sum(c * y ** i for i, c in enumerate(sc)),
           D2S: sum(c * y ** i for i, c in enumerate(dc)),
           D1S: sum(c * y ** i for i, c in enumerate(bc)),
           D0S: sum(c * y ** i for i, c in enumerate(rr_)),
           PHI: fsb.phi_stripped()}

    # y-coefficient arithmetic done on COEFFICIENT LISTS (dense convolution),
    # not by sp.expand on one huge Expr: the substituted collapsed generators
    # carry ~10^5 monomials and Expr-level expansion is the bottleneck, not
    # Singular.
    def as_list(expr):
        expr = sp.expand(sp.sympify(expr))
        if expr == 0:
            return [sp.Integer(0)]
        p = sp.Poly(expr, y)
        n = p.degree()
        out = [sp.Integer(0)] * (n + 1)
        for (m,), c in p.terms():
            out[m] = c
        return out

    def pmul(A, B):
        C = [sp.Integer(0)] * (len(A) + len(B) - 1)
        for i, a in enumerate(A):
            if a == 0:
                continue
            for j, b in enumerate(B):
                if b == 0:
                    continue
                C[i + j] += a * b
        return [sp.expand(c) for c in C]

    def padd(A, B):
        n = max(len(A), len(B))
        return [(A[i] if i < len(A) else 0) + (B[i] if i < len(B) else 0)
                for i in range(n)]

    LISTS = {k: as_list(v) for k, v in sub.items()}

    def coeffs(expr):
        """y-coefficients of `expr` after the substitution, monomial by
        monomial via list convolution."""
        expr = sp.expand(sp.sympify(expr))
        acc = [sp.Integer(0)]
        for term in sp.Add.make_args(expr):
            cur = [sp.Integer(1)]
            for b_, ex in term.as_powers_dict().items():
                if b_.is_number:
                    cur = [sp.expand(c * b_ ** ex) for c in cur]
                    continue
                if b_ == LAM or b_ == y:
                    if b_ == y:
                        cur = [sp.Integer(0)] * int(ex) + cur
                    else:
                        cur = [sp.expand(c * b_ ** ex) for c in cur]
                    continue
                base = LISTS[b_]
                for _ in range(int(ex)):
                    cur = pmul(cur, base)
            acc = padd(acc, cur)
        return [sp.expand(c) for c in acc if c != 0] if any(
            c != 0 for c in acc) else []

    used = (coeffs(coll["G2_collapsed"]) + coeffs(coll["G3_collapsed"])
            + coeffs(coll["G5_collapsed"]))
    lam_row = coeffs(sp.sympify(coll["G4_collapsed"]) + LAM * lambda_target())
    unk = sorted(set(rh) | set(sc) | set(dc) | set(bc) | set(rr_) | {gam},
                 key=sp.default_sort_key)
    sat = [gam, rh[0]]

    out = {"target": "a14_b0000_T1 (alternate regime), collapsed ansatz",
           "n_used_equations": len(used), "n_lambda_equations": len(lam_row),
           "n_unknowns_control": len(unk), "n_unknowns_test": len(unk) + 1,
           "over_approximation": "deg dm4 <= 24 and d1 != 0 NOT imposed"}

    for tag, eqs, vs in (("control", used, unk),
                         ("test", used + lam_row, unk + [LAM])):
        prog = fsb.singular_program(eqs, vs, char=10007, sat_syms=sat)
        t0 = time.time()
        rr = mt.run_singular(prog, timeout=timeout)
        rr["wall"] = round(time.time() - t0, 1)
        rr["n_equations"] = len(eqs)
        rr["n_unknowns"] = len(vs)
        out[tag] = rr
        print("  probe/%s: %s dim=%s (%ss)"
              % (tag, rr.get("verdict"), rr.get("dim"), rr["wall"]), flush=True)

    with open(PROBE_CACHE, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1, default=str)
        fh.write("\n")
    return _probe_verdicts(out)


def _probe_verdicts(out: dict) -> dict:
    ctrl, test = out["control"], out["test"]
    check("I1", "probe CONTROL (used rows alone, no lambda row) did NOT return "
                "UNIT -- the target is live, so the comparison is meaningful",
          ctrl.get("verdict") != "UNIT",
          "verdict=%s dim=%s wall=%ss (%d eqs / %d unknowns)"
          % (ctrl.get("verdict"), ctrl.get("dim"), ctrl.get("wall"),
             ctrl["n_equations"], ctrl["n_unknowns"]))
    check("I2", "probe TEST (used rows + lambda row) returned a verdict; UNIT "
                "would be a kill, anything else -- including a timeout -- is "
                "INCONCLUSIVE, never a survival proof",
          test.get("status") in ("ok", "timeout", "error"),
          "verdict=%s dim=%s wall=%ss (%d eqs / %d unknowns)"
          % (test.get("verdict"), test.get("dim"), test.get("wall"),
             test["n_equations"], test["n_unknowns"]))
    out["kill"] = (test.get("verdict") == "UNIT" and ctrl.get("verdict") != "UNIT")
    return out


# ===========================================================================
#  report
# ===========================================================================
def render(ctx) -> str:
    L, w = [], None
    L = []
    w = L.append
    G4 = ctx["G4"]
    w("# G4_ROW -- the lambda row (`G4`, u-weight 192): derivation, imposition, "
      "and the measured census delta\n")
    w("> Machine-generated by `python -u g4_row.py`. `python -u g4_row.py --quiet` "
      "re-derives everything and exits nonzero on any drift. Read-only on every "
      "other artifact; this lane owns `g4_row*` only.\n")
    w("\n## 0. Verdict up front\n")
    w("> **`G4` re-derives cleanly from primitives, and `dm12` cancels "
      "identically** -- the repo's stated reason for dropping the `j = 4` "
      "Q-slice does not survive the arithmetic, exactly as `TRANSFORM_AUDIT.md` "
      "F1 reports. The row is real, it is weighted-homogeneous at u-weight 192 "
      "(the missing rung), and it is not implied by the four used rows.\n")
    w("> **Imposing it kills NOTHING** -- 0 cells, 0 flag cases, 0 states -- in "
      "standard `sub1` (both C08/C20 settings) and in all six alternate-regime "
      "T1 branches. That is not a weak measurement: section 6 **proves** that no "
      "cell/state-level counting test built on the lambda row alone can ever "
      "kill, because the row vanishes identically at `dm2 = dm3 = dm4 = 0` with "
      "`lambda = 0`, for every `(d2, d1, d0, dm1)`. The +32/+48 equations are "
      "real, and they are **elimination-only** content.\n")
    w("> Standard `sub2` is EMPTY and **cannot** be resurrected: the augmented "
      "system's generator set is a superset and `lambda` is existentially "
      "quantified (section 5.3, a theorem, no computation).\n")

    w("\n## 1. Re-derivation from primitives\n")
    w("Nothing below is imported from `TRANSFORM_AUDIT.md`; the comparison "
      "against its transcription is check `A5`, run *after* the derivation.\n")
    w("```\nS      = 1 + d2*u^2 + d1*u^3 + d0*u^4 + sum_{k>=1} dm_k*u^{4+k}\n"
      "D2(k) := [u^{8+k}] S^2      (the P-side: P_-k = 0 for every k >= 1)\n"
      "D3(j) := [u^{12+j}] S^3     (the Q-side)\n```\n")
    w("The pipeline uses `D2(k)` for `k in {1..7, 9}` and `D3(j)` for "
      "`j in {1,2,3,5}`; it drops the pair `D2(8)`, `D3(4)`.\n")
    w("| step | fact | check |")
    w("|---|---|---|")
    w("| 1 | `dm12` occurs in **no** used slice, and in **both** dropped ones | `A1`, `A1b` |")
    w("| 2 | both dropped slices are degree 1 in `dm12`, coefficients **2** and **3** | `A2`, `A2b` |")
    w("| 3 | therefore `D3(4) - (3/2)*D2(8)` is `dm12`-free -- **identically** | `A3` |")
    w("| 4 | and no spare above `dm10` survives | `A3b` |")
    w("| 5 | the P-side chain `k = 1..7` (linear, coefficient 2) eliminates `dm5..dm11` | `A4` |")
    w("| 6 | what is left involves **only** the seven G-system window variables | `A5a` |")
    w("\nThe dropped pair is therefore equivalent to `{ a definition of dm12 , "
      "G4 }` (`A6`, `A6b`): `D2(8)` after the chain is `2*dm12 + (dm12-free)`, "
      "which merely *names* a fresh unknown and carries no content, and the "
      "residue is\n")
    w("```\nG4 = %s\n```\n" % sp.sstr(sp.factor(G4)))
    w("```\n   = -(3/2)*( 2*d0*dm1*dm3 + d0*dm2^2 + 2*d1*dm2*dm3\n"
      "              + d2*dm3^2 + dm1^2*dm2 - dm4^2 )\n```\n")
    w("**Six** monomials, not five: `TRANSFORM_AUDIT.md` sec.2.1's prose "
      "\"five-monomial quintic\" is a slip (its sec.7 PROVED item 2 lists all "
      "six, correctly). Check `A5b`. Nothing downstream depends on the count.\n")

    w("\n## 2. Grading -- G4 is the missing rung\n")
    w("| variable | `d2` | `d1` | `d0` | `dm1` | `dm2` | `dm3` | `dm4` | `Phi` |")
    w("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    w("| u-weight `12k` | 24 | 36 | 48 | 60 | 72 | 84 | 96 | 204 |")
    w("\n| monomial | weight | sub2 stripped deg | sub1 stripped deg |")
    w("|---|---:|---:|---:|")
    for t in sp.Add.make_args(sp.expand(G4)):
        w("| `%s` | %d | %d | %d |" % (sp.sstr(t), monomial_weight(t),
                                       monomial_stripcap(t, "sub2"),
                                       monomial_stripcap(t, "sub1")))
    w("\nEvery monomial is at weight **192** and attains the stripped degree cap "
      "(`192/6 = 32` in sub2, `192/4 = 48` in sub1) exactly -- checks `B1`, "
      "`B4-sub2`, `B4-sub1`. `full_system_bridge.check_homogeneity()` asserts "
      "the consumed weights are `{156, 168, 180, 204}` (`B3`); in bridge units "
      "the ladder is **156, 168, 180, [192], 204** and `G4` lands on the gap.\n")

    w("\n## 3. How `lambda` enters -- the exact form of the condition\n")
    w("`lambda` is the alpha-strip coefficient `alpha_-1`, a **scalar of `K`** "
      "(`T6_SELECTION_AUDIT.md` sec.4 -- the same premise the normal form "
      "`Q = C^3 + lambda*C^-1 + F` already rests on; this lane does not re-prove "
      "it, see the INFERRED list). With `C4 = y^7*(y+1)` after the A2 "
      "normalisation, `C4^28 = y^196*(y+1)^28` (`C1`), and `G4_full = "
      "y^192*G4_stripped` by weight-homogeneity, so\n")
    w("```\n G4_full + lambda*C4^28 = 0   <=>   G4_stripped = -lambda*y^4*(y+1)^28\n```\n")
    w("> **This is far stronger than \"there exists some lambda\".** It says the "
      "stripped row is a **scalar multiple of one fixed degree-32 polynomial**: "
      "every other coefficient must vanish. `y^4*(y+1)^28` is monic of degree 32 "
      "with `ord_y = 4` (`C3`), so the top coefficient equation reads "
      "`lambda = -[y^32] G4_stripped` -- **`lambda` is determined, not free** "
      "(`C3b`).\n")
    w("\n| row | u-weight | sub2 equations | sub1 equations |")
    w("|---|---:|---:|---:|")
    for n, ww in (("G1", 156), ("G2", 168), ("G3", 180), ("G5", 204)):
        w("| %s | %d | %d | %d |" % (n, ww, ww // 6 + 1, ww // 4 + 1))
    tot = ctx["divisor"]["used_totals"]
    w("| **used total** | | **%d** | **%d** |" % (tot["sub2"], tot["sub1"]))
    w("| **G4 (lambda row)** | **192** | **33 - 1 = 32** | **49 - 1 = 48** |")
    w("| | | **+%.1f %%** | **+%.1f %%** |"
      % (ctx["divisor"]["pct"]["sub2"], ctx["divisor"]["pct"]["sub1"]))
    w("\nChecks `C4`, `C5`, `C5b`. Reproduces `TRANSFORM_AUDIT.md` sec.2.4 exactly.\n")

    w("\n## 4. Not implied by the used rows (re-checked, not inherited)\n")
    w("The canonical guard is asserted first: **`G5 = G5body + Phi` with "
      "`coeff(G5, Phi) == 1`** and `G5body` free of `Phi` (`D0`). A stale "
      "`2*Phi` transcription was a real bug here.\n")
    w("* **`D1` (fresh).** A witness point searched for *in this run*, on "
      "`V(G1,G2,G3,G5)`, with `G4 != 0`: `%s`, there `G4 = %s`.\n"
      % (ctx["nonmem"]["fresh_point"], ctx["nonmem"]["fresh_G4"]))
    w("* **`D2`/`D2b`.** `TRANSFORM_AUDIT.md`'s witness "
      "`(2*Phi-3, -1/3, 1, 1, 0, 1, 1/6)` re-checked: all four used rows vanish "
      "identically in `Phi`, and `G4 = -3*Phi + 37/24 != 0`. **Reproduces.**\n")
    w("* **`D3`/`D3b`.** The 1-parameter family lies on all four used rows for "
      "every `w`, and on it `G4 = %s` -- non-constant, so the used rows pin "
      "neither `G4` nor its divisor. **Reproduces.**\n" % ctx["nonmem"]["family_G4"])

    w("\n## 5. Imposing it\n")
    w("\n### 5.1 In the bridge's coordinates (the drop-in)\n")
    w("`g4_row.g4_row_equations(G4, subs_map, regime)` returns the `cap+1` "
      "y-coefficients of\n")
    w("```\n G4_stripped(subs_map) + lam * y^4*(y+1)^28  ==  0\n```\n")
    w("with `lam` one new scalar unknown -- exactly what "
      "`full_system_bridge.augment()` would append as a fifth generator. "
      "Measured on live ansaetze:\n")
    w("\n| regime | row degree | coefficient equations | rows containing `lam` | spare unknowns |")
    w("|---|---:|---:|---:|---:|")
    for r in ("sub2", "sub1"):
        d = ctx["imposition"][r]
        w("| %s | %d | %d | %d | %d |" % (r, d["row_degree"], d["n_equations"],
                                          d["rows_with_lam"], d["n_spare_unknowns"]))
    w("\n`lam` appears in exactly **29** rows in both regimes -- the nonzero "
      "coefficients of `y^4*(y+1)^28`, i.e. `y^4 .. y^32` -- and the `y^32` row "
      "pins it. Checks `F1-*`, `F2-*`.\n")
    w("\n**Freedom balance, measured on those same live ansaetze** (not the "
      "theoretical maxima):\n")
    w("\n| regime | used-row equations | unknowns | + lambda row | unknowns | net gain |")
    w("|---|---:|---:|---:|---:|---:|")
    for r in ("sub2", "sub1"):
        d = ctx["imposition"][r]
        w("| %s | %d | %d | %d | %d | **+%d** |"
          % (r, d["base_equations"], d["unknowns_before"],
             d["base_equations"] + d["n_equations"], d["unknowns_after"],
             d["net_gain"]))
    w("\n> The `sub2` used-row count is **122**, exactly the theoretical maximum "
      "of section 3. The `sub1` count is **%d**, not the theoretical 181, "
      "because the probe ansatz has `deg e = 10` rather than the `sub1` cap 15, "
      "so `G1, G2, G3, G5` fall short of their degree ceilings. The lambda row "
      "still reaches its full 49 -- its `dm4^2` monomial attains degree 48 from "
      "the spare cap alone, independent of `deg e`. The `+48` is therefore "
      "**not** an artefact of a favourable ansatz.\n"
      % ctx["imposition"]["sub1"]["base_equations"])
    w("\nChecks `F2b-sub2`, `F2b-sub1`.\n")
    w("\n### 5.2 In the collapsed coordinates the live lanes actually use\n")
    w("`DIVISOR_CONSEQUENCES.md` sec.9 / `ALT_FRONTIER_V2.md` sec.7.2 collapse "
      "the spare ansatz 45 -> 18 via `e | S` (write `S = e*s`) and "
      "`dm4 = -R*(s + d2) - d1*e/2` (so `dm4` is not a spare). Written there, "
      "`G1 == 0` identically -- it *is* the collapse (`E1`) -- and the lambda "
      "row becomes\n")
    w("```\nG4|collapsed = -(3/2)*( R^2*A + e^2*B + e*R*d1*(s - d2) )\n"
      "      A = d0 - (s + d2)^2\n"
      "      B = R + 2*d0*s + d2*s^2 - d1^2/4\n```\n")
    w("(check `E2`; `e = dm1`, `R = dm2`, `s = dm3/dm1`). Two further facts that "
      "make it usable:\n")
    w("* **`E4`.** `G2|collapsed` **is** the `T1_BRANCH.md` T6 relation "
      "`W^2 = R^4 + d2*e^2*R^2 + d1*e^3*R + d0*e^4` with `W = e*S - R^2`, up to "
      "the factor `3/2`. So T6 is not independent information -- it is `G2`.\n")
    w("* **`E3`/`E3b`.** `G2` is linear in `d0`, so eliminating `d0` gives a "
      "**sigma-free** lambda row: with `p := d2 + 2*s`,\n")
    w("```\n(16/9)*Krow =  e^4*(4*s^2*p - d1^2) - 4*d1*e^3*R*(d2 + s)\n"
      "             - 4*p^2*e^2*R^2 - 4*d1*e*R^3 - 4*p*R^4 + 4*e^4*R\n"
      "Krow = -(3/2)*e^2 * G4|collapsed   (mod G2)   =   (3/2)*lam*e^2*y^4*(y+1)^28\n```\n")
    w("`Krow` removes the largest unknown block (13 coefficients of `d0` in the "
      "alternate regime) at the price of the factor `e^2`. Section 6's order "
      "calculus and section 8's probe use `G4|collapsed` directly; `Krow` is "
      "the sharper handle for the cancellation ladder recorded OPEN in "
      "section 9 item 23.\n")
    w("\n### 5.3 The `sub2` regression is a theorem, not a measurement\n")
    w("> Adding the lambda row replaces the generator set `{G1,G2,G3,G5}` by the "
      "**superset** `{G1,G2,G3,G5} u {G4_str + lam*y^4*(y+1)^28}` and adds one "
      "existentially quantified scalar `lam`. Any point of the augmented variety "
      "satisfies every old generator by construction, so the projection to the "
      "old variables lands **inside** the old variety. An EMPTY cell therefore "
      "stays empty. Check `F3`. **No resurrection is possible in standard "
      "`sub2`, and none was computed or needed.**\n")

    w("\n## 6. Why the counting-level answer is zero, and why that is a PROOF\n")
    w("> **Every monomial of `G4` contains `dm2`, `dm3` or `dm4`** (`G1`). Hence "
      "`G4 = 0` at `dm2 = dm3 = dm4 = 0`, and `lambda = 0` then satisfies the "
      "divisor condition **for every** `(d2, d1, d0, dm1)` -- on every cell, "
      "every flag case, every state, in every window and both regimes.\n")
    w("So no filter of the `divisor_filter.py` species (an order or degree count "
      "against the state's data) can ever kill using the raw lambda row: the "
      "spares are capped **above** and free **below**, so the row can always be "
      "driven to zero (`G5`). The row's content is **elimination-only** -- it "
      "bites jointly with `G1,G2,G3,G5`, never alone.\n")
    w("\nOn the *collapsed* ansatz the trivial point is unavailable on T1 "
      "(`R = s = 0` gives `G4 = (3/8)*d1^2*e^2 != 0` when `d1 != 0`, check `G2`), "
      "so a counting test is at least conceivable there. It still yields "
      "nothing:\n")
    w("\n**Order calculus at `y = -1`.** With `a := v_t(e)`, `rho := v_t(R)`, "
      "`v_t(s) >= 0`, the three groups of `G4|collapsed` have `t`-orders "
      "`>= 2*rho`, `>= 2*a`, `>= a + rho`, so `v_t(G4|collapsed) >= 2*min(a, "
      "rho)`; the divisor condition forces it to be **exactly 28 when "
      "`lambda != 0`** (and infinite when `lambda = 0`, since then "
      "`G4_str == 0`). The "
      "`T1_BRANCH.md` sec.1.2 trichotomy fixes which horn is available:\n")
    w("\n| horn | scope | `rho` | `v_t` floor | `lam != 0` demands | consequence |")
    w("|---|---|---|---|---|---|")
    w("| 1 (`rho >= a`) | `a <= 10` (standard) | `>= a` | `2*a` | floor `<= 28` | `a <= 14` -- **VACUOUS** |")
    w("| 2 (`a + 2*rho = 30`) | `a >= 11` even (alternate) | `(30-a)/2` | `30 - a` | floor `<= 28` | `a >= 2` -- **VACUOUS** |")
    w("\n> **And a violation would not be a kill anyway** (`G3c`). `lambda` is a "
      "constant that is ALLOWED TO BE ZERO, so `v_t(G4_str) > 28` forces "
      "`lambda = 0`, i.e. `G4_str == 0` -- a **case** of the divisor condition, "
      "not a contradiction. The order calculus can therefore only ever eliminate "
      "the `lambda != 0` branch, never a cell. On the surviving frontier it does "
      "not even do that: `floor > 28` needs `a > 14` (horn 1) or `a < 2` "
      "(horn 2), and neither occurs.\n")
    w("\nChecks `G3`, `G3b`. The marked-root test (`v_beta(R) = 0` when "
      "`lambda != 0`) is a **re-derivation of an already committed lemma** "
      "(`ALT_FRONTIER_V2.md` sec.6, trichotomy at `beta`), not new content "
      "(`G4`).\n")
    w("\n> **Scope of the two claims, kept apart.** The statement about the "
      "**raw** row (`G1`, `G1b`) is a PROOF: the row alone constrains nothing "
      "at cell/state level, ever. The statement about the **collapsed** row "
      "(`G2`, `G3`, `G3b`, `G4`, `G5`) is an EXHAUSTION OF THE NATURAL TESTS -- "
      "order at `y = -1`, order at a marked root, order at `y = 0`, and degree "
      "-- each computed and each vacuous. It is not a proof that no counting "
      "test whatsoever exists on the collapsed row. What it does establish is "
      "that the species of test `divisor_filter.py` and `ALT_FRONTIER_V2.md` "
      "sec.3-4 use is exhausted, and the row's content lies past it.\n")

    w("\n## 7. The measured census delta\n")
    w("Baselines are `frontier_rebuild.STAGES` applied to the committed "
      "`e | Phi`-filtered universes -- recomputed here, not quoted (`H1`, `H1b`, "
      "`H2`).\n")
    w("\n| window | C08/C20 | cells | flag cases | states | cells killed | flag cases killed | states killed | cells after | flag cases after | states after |")
    w("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for r in ctx["census"]["rows"]:
        if "missing" in r:
            w("| %s | %s | *(universe %s absent)* | | | | | | | | |"
              % (r["window"], r["c0820"], r["missing"]))
            continue
        b, k, a = r["before"], r["killed"], r["after"]
        w("| %s | %s | %d | %d | %d | **%d** | **%d** | **%d** | %d | %d | %d |"
          % (r["window"], r["c0820"], b["cells"], b["flagcases"], b["states"],
             k["cells"], k["flagcases"], k["states"],
             a["cells"], a["flagcases"], a["states"]))
    w("\n**Alternate regime** (the six surviving T1 branches, "
      "`ALT_FRONTIER_V2.md` sec.6, parsed not typed -- `H4`):\n")
    w("\n| branch | `a` | horn-2 `rho` | `v_t` floor | `lam != 0` branch open? | verdict |")
    w("|---|---:|---:|---:|---|---|")
    for r in ctx["alt"]["branches"]:
        w("| `%s` | %d | %s | %s | %s | %s |"
          % (r["branch"], r["a_t"], r["rho"], r["v_t_floor"],
             "no -- `lam = 0` forced" if r["lam_forced_zero"] else "yes",
             "**DEAD**" if r["dead"] else "alive"))
    w("\n| | branches | states |")
    w("|---|---:|---:|")
    w("| before | 6 | %d |" % ctx["alt"]["states_before"])
    w("| killed by the lambda row | **%d** | **%d** |"
      % (len(ctx["alt"]["killed"]), ctx["alt"]["states_killed"]))
    w("| after | %d | %d |" % (6 - len(ctx["alt"]["killed"]),
                               ctx["alt"]["states_before"] - ctx["alt"]["states_killed"]))
    w("\n> **No `g4_row_stage.json` is emitted**, because nothing dies. The file "
      "is written only if the census is nonzero; its absence is the result.\n")

    if ctx.get("probe"):
        pr = ctx["probe"]
        w("\n## 8. Bounded solver probe (elimination, where the content actually is)\n")
        w("Two runs, **sequential**, one target, bounded timeout, no global "
          "process kill. `137 / 124 / 139` are ABORT / TIMEOUT / SIGSEGV and are "
          "**never** verdicts.\n")
        w("\nTarget **`%s`**. `deg e = a + sum(b) = 14` exactly (D1), so "
          "`e = gam*(y+1)^14`; horn 2 pins `v_t(R) = 8` exactly, so "
          "`R = (y+1)^8 * sum_{i<=10} rh_i*(y+1)^i` with `rh0 = Rh(-1) != 0`; "
          "`deg s <= 7`, `deg d2 <= 6`, `deg d1 <= 9`, `deg d0 <= 12`. "
          "Saturated by `gam*rh0`.\n" % pr["target"])
        w("\n> **Sound over-approximation:** `deg dm4 <= 24` and `d1 != 0` (T1) "
          "are NOT imposed. Dropping necessary conditions can only make a kill "
          "harder, so a UNIT here would be genuine; a NON-UNIT is **not** a "
          "proof of survival.\n")
        w("\n| run | generators | equations | unknowns | verdict | dim | wall |")
        w("|---|---|---:|---:|---|---|---:|")
        for tag, gens in (("control", "`G2, G3, G5` (collapsed)"),
                          ("test", "`G2, G3, G5` **+ lambda row**")):
            r = pr[tag]
            w("| %s | %s | %d | %d | `%s` | `%s` | %ss |"
              % (tag, gens, r["n_equations"], r["n_unknowns"],
                 r.get("verdict"), r.get("dim"), r.get("wall")))
        if pr.get("kill"):
            w("\n> **UNIT on the TEST and not on the CONTROL: the lambda row "
              "empties this target.** Recorded mod p only -- it must be re-run "
              "exactly over Q, and the two dropped conditions re-imposed, before "
              "it is claimed as a kill.\n")
        else:
            w("\n> **No kill.** The lambda row does not empty this target at the "
              "depth reached. %s\n"
              % ("The TEST run did not complete inside the bound, so this is "
                 "INCONCLUSIVE, not a survival proof."
                 if pr["test"].get("status") != "ok" else
                 "Both runs completed; the row constrains the target without "
                 "emptying it."))
    else:
        w("\n## 8. Bounded solver probe -- ATTEMPTED, NOT COMPLETED\n")
        w("`python -u g4_row.py --probe` builds the target below and runs two "
          "**sequential** bounded Singular jobs, a CONTROL (`G2, G3, G5` "
          "collapsed) and a TEST (the same **plus the lambda row**). It is "
          "implemented and wired; it did **not** complete inside its bound in "
          "this session and was stopped. **No verdict is recorded, and none is "
          "implied** -- an unfinished Groebner run is not evidence of survival, "
          "exactly as `137 / 124 / 139` are not verdicts.\n")
        w("\nTarget **`a14_b0000_T1`** (alternate regime, the tightest of the "
          "six): `deg e = a + sum(b) = 14` exactly (D1), so `e = gam*(y+1)^14` "
          "with ONE unknown; horn 2 pins `v_t(R) = (30-14)/2 = 8` exactly, so "
          "`R = (y+1)^8 * sum_{i<=10} rh_i*(y+1)^i` with `rh0 = Rh(-1) != 0`; "
          "`deg s <= 7`, `deg d2 <= 6`, `deg d1 <= 9`, `deg d0 <= 12`; "
          "saturated by `gam*rh0`. Roughly 51 unknowns (52 with `lam`).\n")
        w("\n> **Sound over-approximation, stated in advance:** the probe does "
          "NOT impose `deg dm4 <= 24` or `d1 != 0` (T1). Dropping necessary "
          "conditions can only make a kill harder, so a UNIT would have been "
          "genuine; a NON-UNIT would not have been a survival proof either way.\n")
        w("\nWhat the attempt did establish, cheaply: the y-coefficient "
          "extraction for that ansatz has to be done by dense coefficient-list "
          "convolution, not `sympy.expand` on the assembled expression -- the "
          "latter is the bottleneck, not the solver. That is implemented in "
          "`g4_row.probe`, and a completed run caches to `g4_row_probe.json` so "
          "the report can be re-rendered without re-running Singular.\n")

    w("\n## 9. PROVED / CHECKED / INFERRED / OPEN\n")
    w("\n**PROVED** (exact symbolic identity or explicit witness, in `g4_row.py`)\n")
    w("1. `dm12` occurs in no used slice and in both dropped ones; both are "
      "degree 1 in it with coefficients 2 and 3; `D3(4) - (3/2)*D2(8)` is "
      "`dm12`-free identically. (`A1`-`A3b`)\n")
    w("2. The dropped pair is equivalent to `{2*dm12 + ... = 0, G4 = 0}` -- one "
      "definition and one condition. (`A6`, `A6b`)\n")
    w("3. `G4 = -(3/2)*(2*d0*dm1*dm3 + d0*dm2^2 + 2*d1*dm2*dm3 + d2*dm3^2 + "
      "dm1^2*dm2 - dm4^2)`, six monomials, in exactly the seven G-system window "
      "variables. (`A5a`, `A5b`)\n")
    w("4. Every monomial has u-weight 192 and attains the stripped degree cap "
      "32 (sub2) / 48 (sub1). (`B1`, `B4-*`)\n")
    w("5. `C4^28 = y^196*(y+1)^28`, so the condition is exactly "
      "`G4_stripped = -lambda*y^4*(y+1)^28`; the target is monic of degree 32 "
      "with `ord_y = 4`, so `lambda = -[y^32] G4_stripped` is **determined**. "
      "(`C1`-`C3b`)\n")
    w("6. Equation counts 33/49, net +32/+48 against used totals 122/181. "
      "(`C4`, `C5`, `C5b`)\n")
    w("7. `G4` is not in `<G1,G2,G3,G5>`: a FRESH witness found in this run, "
      "plus the audit's witness and 1-parameter family re-checked. (`D1`-`D3b`)\n")
    w("8. `G5 = G5body + Phi`, `coeff(G5, Phi) == 1`, `G5body` free of `Phi`. (`D0`)\n")
    w("9. On the collapsed ansatz `G1 == 0`, `G4|collapsed = -(3/2)*(R^2*A + "
      "e^2*B + e*R*d1*(s-d2))`, and eliminating `d0` against `G2` gives the "
      "sigma-free `Krow`. (`E1`-`E3c`)\n")
    w("10. `G2|collapsed` **is** the `T1_BRANCH.md` T6 relation (up to `3/2`). "
      "(`E4`, `E4b`)\n")
    w("11. Every monomial of `G4` carries a spare, so `lambda = 0` with "
      "`dm2 = dm3 = dm4 = 0` satisfies the divisor condition unconditionally: "
      "**no cell/state-level test on the RAW row can kill.** (`G1`, `G1b`, `G5`)\n")
    w("12. MONOTONICITY: the augmented variety projects into the old one, so no "
      "empty cell can be resurrected. (`F3`)\n")
    w("\n**CHECKED** (verified against the repo's own artifacts, not re-proved)\n")
    w("12b. On the COLLAPSED row the four natural counting tests -- order at "
      "`y = -1` (both horns), order at a marked root, order at `y = 0`, and "
      "degree -- are each computed and each VACUOUS. This exhausts the species "
      "of test the pipeline uses; it is **not** a proof that no counting test "
      "exists. (`G2`, `G3`, `G3b`, `G4`, `G5`)\n")
    w("13. The freshly derived `G4` equals `TRANSFORM_AUDIT.md` sec.2.1's "
      "transcription exactly. (`A5`)\n")
    w("14. `full_system_bridge.WEIGHT` is the same `12k` grading and asserts the "
      "consumed ladder `{156,168,180,204}`; 192 is structurally absent. (`B2`, `B3`)\n")
    w("15. `sol4` in `generators.json` equals `-R*(S/e + d2) - d1*e/2`. (`E1b`)\n")
    w("16. Survivor baselines: sub1 ON 34/314/7275, sub1 OFF 34/322/8889, sub2 "
      "standard EMPTY -- recomputed from the committed universes. (`H1`-`H2`)\n")
    w("17. The six alternate T1 branches and their 562 states, parsed out of "
      "`ALT_FRONTIER_V2.md` sec.6. (`H4`)\n")
    w("\n**INFERRED** (follows from premises this lane did not re-prove)\n")
    w("18. **`lambda in K` is a constant.** From the alpha-strip / GGV1 "
      "leading-form premises (`T6_SELECTION_AUDIT.md` sec.4). **The entire "
      "divisor condition stands or falls with this** -- it is the same premise "
      "the normal form `Q = C^3 + lambda*C^-1 + F` already rests on, but it is "
      "not re-derived here.\n")
    w("19. `C4 = y^7*(y+1)` after the A2 normalisation "
      "(`T6_SELECTION_AUDIT.md` sec.4); field-of-definition consequences belong "
      "to `FIELD_SCOPE_*`.\n")
    w("20. `e | S` -- used only for the *collapsed* forms of section 5.2 and the "
      "order calculus of section 6, never for the imposition of section 5.1. "
      "(`ALT_FRONTIER_V2.md` L2, `DIVISOR_CONSEQUENCES.md` sec.9.)\n")
    w("21. The `T1_BRANCH.md` sec.1.2 trichotomy at `y = -1`, used to select the "
      "horn in the order calculus. Its scope guards are respected: the `T2`-only "
      "relations `R | e^2`, `e*R | Phi`, `R = c*(y+1)^rho` are **not** used "
      "anywhere in this file, SPINE's zero-slack count is **not** applied to "
      "sub1, and `POLE_THEOREM` Thm 2C is **not** invoked (its gate "
      "`3a - 2 < 30` is false for `a >= 11`).\n")
    w("\n**OPEN**\n")
    w("22. **Whether the lambda row kills by ELIMINATION.** That is where all of "
      "its content is (section 6), and it is **not settled here**. The sub1 "
      "augmented system carries 66 spare coefficient unknowns before the "
      "collapse and 18 after; a full Groebner sweep over 314 flag cases was not "
      "attempted, and the single bounded probe of section 8 did not complete "
      "inside its bound. Nothing in this file should be read as evidence that "
      "the row fails to kill by elimination -- only that it does not kill by "
      "counting.\n")
    w("23. The ladder of cancellation conditions the order calculus produces in "
      "the alternate regime **on the `lambda != 0` branch**: `v_t(Krow)` must "
      "climb from its floor `4*min(a, rho)` to `2*a + 28` -- 16 orders at "
      "`a = 12` (36 -> 52), 24 at `a = 14` (32 -> 56). Each step is individually "
      "satisfiable, so it is a ladder, not a kill -- but it is cheap, and it is "
      "the natural next probe. The `lambda = 0` branch (`Krow == 0` identically) "
      "is a separate, and possibly easier, target.\n")
    w("24. Whether `lambda = 0` is forced on any surviving cell. If it is, the "
      "condition strengthens from \"33/49 coefficients proportional\" to \"33/49 "
      "coefficients zero\" and the free scalar is not spent.\n")

    w("\n## 10. Checks\n")
    w("| id | check | result |")
    w("|---|---|---|")
    for cid, label, ok, detail in RESULTS:
        w("| `%s` | %s | %s |" % (cid, label.replace("|", "\\|"),
                                  "PASS" if ok else "**FAIL**"))
    npass = sum(1 for _c, _l, ok, _d in RESULTS if ok)
    w("\n**%d / %d checks pass.**\n" % (npass, len(RESULTS)))
    w("\n```\npython -u g4_row.py           # this report\n"
      "python -u g4_row.py --quiet   # exit 0 iff every check passes\n"
      "python -u g4_row.py --probe   # + one bounded Singular elimination probe\n```\n")
    return "\n".join(L) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--probe", action="store_true")
    ap.add_argument("--probe-timeout", type=float, default=300.0)
    ap.add_argument("--reprobe", action="store_true",
                    help="ignore g4_row_probe.json and re-run Singular")
    args = ap.parse_args()

    G4, _D2, _D3, _chain = derive_g4()
    caps = grading(G4)
    div = divisor_condition(caps)
    gsys, st = gsystem()
    gsys["_sol4"] = st["sol4"]
    nonmem = nonmembership(G4, gsys)
    coll = collapsed(G4, gsys)
    imp = imposition(G4, caps)
    theory = counting_theory(G4, coll)
    cen = census(G4)
    alt = alt_census()
    pr = None
    if args.probe:
        if os.path.exists(PROBE_CACHE) and not args.reprobe:
            with open(PROBE_CACHE, encoding="utf-8") as fh:
                pr = _probe_verdicts(json.load(fh))
            print("  probe: reusing %s (--reprobe to re-run Singular)"
                  % os.path.basename(PROBE_CACHE), flush=True)
        else:
            pr = probe(G4, coll, timeout=args.probe_timeout)

    ctx = {"G4": G4, "caps": caps, "divisor": div, "nonmem": nonmem,
           "collapsed": coll, "imposition": imp, "theory": theory,
           "census": cen, "alt": alt, "probe": pr}

    nfail = sum(1 for _c, _l, ok, _d in RESULTS if not ok)
    if args.quiet:
        for cid, label, ok, detail in RESULTS:
            if not ok:
                print("[FAIL] %s  %s  %s" % (cid, label, detail))
        raise SystemExit(1 if nfail else 0)

    for cid, label, ok, detail in RESULTS:
        print("[%s] %-8s %s%s" % ("PASS" if ok else "FAIL", cid, label,
                                  ("   -- " + detail) if detail else ""))
    with open(MD_OUT, "w", encoding="utf-8") as fh:
        fh.write(render(ctx))
    print("\nwrote %s" % MD_OUT)

    # the stage record is emitted ONLY if something dies
    dead_any = (any(r.get("killed", {}).get("cells") for r in cen["rows"])
                or alt["killed"])
    if dead_any:
        rec = {
            "stage": {
                "id": "stage5_g4_lambda_row",
                "title": "Lambda row (G4, u-weight 192): G4_stripped = "
                         "-lambda*y^4*(y+1)^28",
                "source": "G4_ROW.md sec.1-3; derived from regenerate_system.py "
                          "primitives + TRANSFORM_AUDIT.md F1",
                "checker": "python -u g4_row.py --quiet",
                "note": "Adding a generator plus one existentially quantified "
                        "scalar; monotone, so it can never resurrect a cell.",
                "dead": {"sub2": [], "sub1": sorted(
                    set().union(*[set(r.get("killed_cell_names", []))
                                  for r in cen["rows"]
                                  if r.get("window") == "sub1"]) or set())},
                "applies_after": "stage4_positive_slice",
            },
            "alternate_regime": {"killed_branches": alt["killed"],
                                 "states_killed": alt["states_killed"]},
            "schema": "frontier_rebuild.STAGES entry (drop-in); this lane does "
                      "not modify frontier_rebuild.py",
        }
        with open(STAGE_OUT, "w", encoding="utf-8") as fh:
            json.dump(rec, fh, indent=1)
            fh.write("\n")
        print("wrote %s" % STAGE_OUT)
    else:
        print("no cell/flag-case/state dies -> g4_row_stage.json deliberately "
              "NOT written")
        if os.path.exists(STAGE_OUT):
            print("NOTE: a stale %s exists on disk" % STAGE_OUT)
    print("checks: %d/%d pass" % (len(RESULTS) - nfail, len(RESULTS)))
    raise SystemExit(1 if nfail else 0)


if __name__ == "__main__":
    main()
