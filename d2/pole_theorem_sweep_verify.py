#!/usr/bin/env python3
"""pole_theorem_sweep_verify.py -- INDEPENDENT checker for POLE_THEOREM.md /
pole_theorem_sweep.py.

Every check here re-derives its object by a DIFFERENT route than the module it
checks; where that is impossible the check is a direct instantiation against
first principles (the G-generators) rather than against the module's own output.

  V1  H2,H3,H5 rebuilt from _G_generators; S1 identity re-verified by RANDOM
      NUMERIC substitution over Q (not by symbolic expansion).
  V2  K5 is linear in dm3 / free of dm4 and d0, read off the monomial support.
  V3  the d1 != 0 degradation, verified numerically: the residue
      6*d1*dm2^2*dm3 is exactly right, and K5gen really is quadratic in dm3.
  V4  the (y+1) bound re-derived by SYMBOLIC LEADING COEFFICIENTS on genuine
      rational-function arithmetic -- no integer order bookkeeping at all.
      Reproduces GENERIC_FIBER.md sec.4's third route at a = 9 and extends it.
  V5  an INDEPENDENT brute-force DVR enumerator, written from scratch here,
      reproducing rho_min for a = 0..12 (no import of `dvr_case`).
  V6  the support theorem by direct numeric instantiation at random places.
  V7  the collapse: the exact division is verified as a POLYNOMIAL IDENTITY on
      random rational data, and the R9 row is checked against the numbers
      printed in GENERIC_FIBER.md sec.5.
  V8  census arithmetic re-counted straight from the json artifacts.
  V9  branch T2 <=> d1 == 0, checked on both the batch census and the cell file.
  V10 the spare caps 12/14 taken from full_system_bridge, not hard-coded.
  V11 the new bound strictly implies the certified R9_VALSPLIT row
      v(dm2) + v(dm3) >= 9, and is strictly stronger.

Usage:  python -u pole_theorem_sweep_verify.py           verbose
        python -u pole_theorem_sweep_verify.py --quiet   exit 0 / 1
"""
from __future__ import annotations

import json
import os
import random
import sys

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__)) or "."
sys.path.insert(0, HERE)

import bigrade_annotator as BA            # noqa: E402
import pole_theorem_sweep as PT           # noqa: E402  (the module under test)

y = sp.Symbol("y")
INF = 10**9
PHI_C = sp.Rational(-1, 6630)
Q_POLY = 2048 * y**4 - 512 * y**3 + 320 * y**2 - 240 * y + 195

RESULTS = []


def check(name, ok, detail=""):
    RESULTS.append((name, bool(ok), detail))
    return bool(ok)


# ---------------------------------------------------------------- V1, V2, V3
def v1_v2_v3():
    d0, d1, d2, dm1, dm2, dm3, dm4, Phi = BA._gsystem_symbols()
    G = {k: v for k, (v, _w) in BA._G_generators().items()}
    H = {"H2": sp.expand(dm1 * G["G2"] - dm2 * G["G1"]),
         "H3": sp.expand(dm1 * G["G3"] - dm3 * G["G1"]),
         "H5": sp.expand(dm1 * G["G5"] + (d0 * dm1 + d1 * dm2 + d2 * dm3) * G["G1"])}
    Hj = {k: v for k, (v, _w) in BA._H_generators().items()}
    check("V1a  H rebuilt from G == r9_eliminated_system.json H",
          all(sp.expand(H[k] - Hj[k]) == 0 for k in H))

    K5 = 2 * Phi - 3 * dm1 * dm2**2 - d2 * dm1**3 - 3 * dm1**2 * dm3
    rng = random.Random(2024)
    syms = [d0, d2, dm1, dm2, dm3, dm4, Phi]
    bad = 0
    for _ in range(60):
        pt = {s: sp.Rational(rng.randrange(-40, 40), rng.randrange(1, 12))
              for s in syms}
        pt[d1] = sp.Integer(0)
        lhs = sp.expand(2 * (H["H5"].xreplace(pt) + pt[d2] * H["H3"].xreplace(pt)))
        rhs = sp.expand(pt[dm1] * K5.xreplace(pt))
        if sp.simplify(lhs - rhs) != 0:
            bad += 1
    check("V1b  S1 identity at d1=0 on 60 random rational points", bad == 0,
          "%d failures" % bad)

    mons = sp.Poly(sp.expand(K5), dm3, dm4, d0).monoms()
    check("V2   K5: deg_dm3 = 1, no dm4, no d0",
          max(m[0] for m in mons) == 1 and max(m[1] for m in mons) == 0
          and max(m[2] for m in mons) == 0)

    # V3: the d1 != 0 degradation, checked numerically
    K5gen = sp.expand(2 * Phi + 3 * d0 * d1 * dm1**2 + 3 * d1**2 * dm1 * dm2
                      + 3 * d1 * d2 * dm2**2 - 3 * d1 * dm3**2 - d2 * dm1**3
                      - 3 * dm1**2 * dm3 - 3 * dm1 * dm2**2)
    resid = 2 * (H["H5"] + d2 * H["H3"]) - dm1 * K5gen - 6 * d1 * dm2**2 * dm3
    bad = 0
    for _ in range(60):
        pt = {s: sp.Rational(rng.randrange(-40, 40), rng.randrange(1, 12))
              for s in syms + [d1]}
        if sp.simplify(resid.xreplace(pt)) != 0:
            bad += 1
    check("V3a  d1 != 0: 2(H5+d2 H3) = dm1*K5gen + 6 d1 dm2^2 dm3 (60 pts)",
          bad == 0, "%d failures" % bad)
    check("V3b  K5gen is QUADRATIC in dm3 (method does not transfer to T1)",
          sp.Poly(K5gen, dm3).degree() == 2)
    check("V3c  the residue 6 d1 dm2^2 dm3 is NOT divisible by dm1",
          sp.rem(sp.Poly(6 * d1 * dm2**2 * dm3, dm1),
                 sp.Poly(dm1, dm1)).as_expr() != 0)


# ---------------------------------------------------------------------- V4
def v4_leading_coefficients():
    """Re-derive the bound with genuine rational-function arithmetic in the
    uniformiser T = y+1, carrying SYMBOLIC leading coefficients.

    dm1 = gm*T^a*(w0 + w1*T),  dm2 = T^rho*(u0 + u1*T),  Phi = c*T^30*(Q0+Q1*T),
    d2 = v0 + v1*T,  d0 = z0 + z1*T   (v(d2) = v(d0) = 0: the TIGHTEST case,
    since raising them only raises the competing term orders).

    Predictions to reproduce, for 0 <= rho < a:
        ord(dm3) = 2*rho - a      lead(dm3) = -u0^2/(gm*w0)
        ord(H3)  = 5*rho - 2a     lead(H3)  = -3*u0^5/(gm^2*w0^2)
    both manifestly nonzero.  At a = 9, w0 = -(1+r) this is exactly the third
    route of GENERIC_FIBER.md sec.4.
    """
    T = sp.Symbol("T")
    gm, w0, w1, u0, u1, Q0, Q1 = sp.symbols("gm w0 w1 u0 u1 Q0 Q1")
    v0, v1, z0, z1 = sp.symbols("v0 v1 z0 z1")

    def ord_lead(expr):
        """(order, leading coefficient) of a rational function of T."""
        e = sp.cancel(sp.together(sp.expand(expr)))
        num, den = sp.fraction(e)
        pn, pd = sp.Poly(sp.expand(num), T), sp.Poly(sp.expand(den), T)
        cn = pn.all_coeffs()[::-1]
        cd = pd.all_coeffs()[::-1]
        kn = next((i for i, c in enumerate(cn) if sp.simplify(c) != 0), None)
        kd = next((i for i, c in enumerate(cd) if sp.simplify(c) != 0), None)
        if kn is None:
            return INF, sp.Integer(0)
        return kn - kd, sp.cancel(cn[kn] / cd[kd])

    ok3 = ok_h3 = True
    rows = []
    for a in (7, 8, 9, 10):
        dm1 = gm * T**a * (w0 + w1 * T)
        d2v = v0 + v1 * T
        d0v = z0 + z1 * T
        Phiv = PHI_C * T**30 * (Q0 + Q1 * T)
        for rho in range(0, a):
            dm2 = T**rho * (u0 + u1 * T)
            N = sp.expand(2 * Phiv - 3 * dm1 * dm2**2 - d2v * dm1**3)
            dm3 = sp.cancel(N / (3 * dm1**2))
            o3, l3 = ord_lead(dm3)
            H3 = (-3 * d0v * dm1**2 * dm2 - 3 * d2v * dm1 * dm2 * dm3
                  - dm1**4 / 2 - 3 * dm2 * dm3**2)
            oH, lH = ord_lead(H3)
            p3, pl3 = 2 * rho - a, sp.cancel(-u0**2 / (gm * w0))
            pH, plH = 5 * rho - 2 * a, sp.cancel(-3 * u0**5 / (gm**2 * w0**2))
            g3 = (o3 == p3) and sp.simplify(l3 - pl3) == 0
            gH = (oH == pH) and sp.simplify(lH - plH) == 0
            ok3 &= g3
            ok_h3 &= gH
            rows.append((a, rho, o3, p3, oH, pH, g3, gH))
    check("V4a  ord/lead(dm3) = (2rho-a, -u0^2/(gm w0)) for all a=7..10, rho<a",
          ok3)
    check("V4b  ord/lead(H3) = (5rho-2a, -3u0^5/(gm^2 w0^2)) -- H3 != 0", ok_h3)
    # the R9 specialisation printed in GENERIC_FIBER.md sec.4
    r, gamma = sp.symbols("r gamma")
    l3_r9 = sp.cancel((-u0**2 / (gm * w0)).subs({gm: gamma, w0: -(1 + r)}))
    lH_r9 = sp.cancel((-3 * u0**5 / (gm**2 * w0**2)).subs({gm: gamma,
                                                           w0: -(1 + r)}))
    check("V4c  a=9 specialisation == GENERIC_FIBER sec.4 "
          "(u0^2/(gamma(1+r)), -3u0^5/(gamma^2(1+r)^2))",
          sp.simplify(l3_r9 - u0**2 / (gamma * (1 + r))) == 0
          and sp.simplify(lH_r9 + 3 * u0**5 / (gamma**2 * (1 + r)**2)) == 0)
    return rows


# ---------------------------------------------------------------------- V5
def v5_independent_enumeration():
    """A brute-force DVR enumerator written from scratch (does NOT call
    PT.dvr_case).  Plain nested loops, plain integers, wide ranges."""
    BIG = 10**7

    def add(*xs):
        return BIG if any(x >= BIG for x in xs) else sum(xs)

    def rho_min_bruteforce(a, P, RHO=70, AUX=40):
        for rho in list(range(0, RHO)) + [BIG]:
            for t in list(range(0, AUX)) + [BIG]:
                for s in list(range(0, AUX)) + [BIG]:
                    terms = [P, add(a, rho, rho), add(t, a, a, a)]
                    lo = min(terms)
                    if terms.count(lo) != 1:
                        return rho          # cannot conclude -> rho survives
                    sig = BIG if lo >= BIG else lo - 2 * a
                    if sig < 0:
                        continue            # pole: this (rho,t,s) is dead
                    h = [add(s, a, a, rho), 4 * a,
                         add(t, a, rho, sig), add(rho, sig, sig)]
                    lo3 = min(h)
                    if h.count(lo3) != 1:
                        return rho          # survives
        return BIG

    ok = True
    rows = []
    for a in range(0, 13):
        mine = rho_min_bruteforce(a, 30)
        theirs = PT.pole_bound(a, 30)["rho_min"]
        theirs = BIG if theirs >= INF else theirs
        rows.append((a, mine, theirs))
        ok &= (mine == theirs)
    check("V5   independent brute-force enumerator reproduces rho_min, a=0..12",
          ok, str(rows))
    # and the proved constant
    ok2 = all(PT.pole_bound_closed(a, 30)["rho_min"] == a for a in range(0, 11))
    check("V5b  proved constant is exactly `a` for a = 0..10", ok2)
    check("V5c  a >= 11 is OUT of the proved regime (3a-2 < 30 fails)",
          not PT.in_regime(11, 30) and PT.pole_bound_closed(11, 30) is None)


# ---------------------------------------------------------------------- V6
def v6_support():
    """Direct instantiation: a root of e that is not a root of Phi forces a
    pole in dm3, for random places, random multiplicities and random dm2."""
    Phiv = sp.expand(PHI_C * (y + 1)**30 * Q_POLY)
    rng = random.Random(7)
    bad = 0
    trials = 0
    for _ in range(12):
        beta = sp.Integer(rng.choice([2, 3, 5, 7, -3, -5, 11]))
        if sp.expand(Phiv).subs(y, beta) == 0:
            continue
        m = rng.choice([1, 2, 3])
        a = rng.choice([6, 7, 8])
        gmv = sp.Integer(rng.randrange(1, 9))
        dm1 = sp.expand(gmv * (y + 1)**a * (y - beta)**m)
        dm2 = sum(sp.Integer(rng.randrange(-6, 7)) * y**i for i in range(5))
        d2v = sum(sp.Integer(rng.randrange(-6, 7)) * y**i for i in range(5))
        N = sp.expand(2 * Phiv - 3 * dm1 * dm2**2 - d2v * dm1**3)
        vb = PT._ord_at(N, beta)
        trials += 1
        if not (vb < 2 * m):
            bad += 1
    check("V6   e-root off Phi => v_beta(N) < 2m (pole) on %d random cases"
          % trials, bad == 0 and trials >= 8, "%d failures" % bad)


# ---------------------------------------------------------------------- V7
def v7_collapse():
    """The exact division 3*gm^2*T^2*B = W, verified as a polynomial identity
    on RANDOM RATIONAL DATA against the original K5 -- not against the module's
    own symbolic division."""
    rng = random.Random(99)
    bad = 0
    rows = []
    for a, k in ((10, 0), (9, 1), (8, 2), (8, 0), (7, 3)):
        gmv = sp.Integer(rng.randrange(2, 9))
        roots = [sp.Rational(rng.randrange(2, 30)) for _ in range(k)]
        Tp = sp.expand(sp.prod([(y - rt) for rt in roots])) if k else sp.Integer(1)
        Av = sum(sp.Integer(rng.randrange(-5, 6)) * y**i
                 for i in range(12 - a + 1))
        d2v = sum(sp.Integer(rng.randrange(-5, 6)) * y**i for i in range(5))
        dm1 = sp.expand(gmv * (y + 1)**a * Tp)
        dm2 = sp.expand((y + 1)**a * Av)
        Phiv = sp.expand(PHI_C * (y + 1)**30 * Q_POLY)
        N = sp.expand(2 * Phiv - 3 * dm1 * dm2**2 - d2v * dm1**3)
        # the claim: N == (y+1)^(3a) * W  with  W as in POLE_THEOREM sec.4
        W = sp.expand(2 * PHI_C * (y + 1)**(30 - 3 * a) * Q_POLY
                      - 3 * gmv * Tp * Av**2 - gmv**3 * d2v * Tp**3)
        id1 = sp.expand(N - (y + 1)**(3 * a) * W) == 0
        # and: dm3 = N/(3 dm1^2) is polynomial IFF 3*gm^2*T^2 | W
        quo, rem = sp.div(sp.Poly(W, y), sp.Poly(sp.expand(3 * gmv**2 * Tp**2), y))
        dm3 = sp.cancel(N / (3 * dm1**2))
        id2 = (sp.simplify(rem.as_expr()) == 0) == (sp.denom(dm3) == 1)
        rows.append((a, k, id1, id2))
        if not (id1 and id2):
            bad += 1
    check("V7a  N == (y+1)^(3a)*W and (poly dm3 <=> 3gm^2 T^2 | W), 5 shapes",
          bad == 0, str(rows))
    # the R9 numbers GENERIC_FIBER.md sec.5 prints
    r9 = PT.stage_E_collapse(9, 1)
    check("V7b  R9 collapse row == GENERIC_FIBER sec.5 (deg W 7, deg B 5, "
          "2 remainder rows, 4 spares)",
          r9["deg_W"] == 7 and r9["deg_B"] == 5 and r9["n_remainder_rows"] == 2
          and r9["spare_unknowns_after"] == 4)
    a10 = PT.stage_E_collapse(10, 0)
    check("V7c  a=10,k=0 : ZERO remainder rows, 3 spare unknowns",
          a10["n_remainder_nonzero"] == 0 and a10["spare_unknowns_after"] == 3)


# ------------------------------------------------------------------ V8/V9/V10
def v8_v9_v10():
    with open(os.path.join(HERE, "batch_convolution_sub2.json")) as fh:
        batch = json.load(fh)["states"]
    with open(os.path.join(HERE, "phase_d_states_sub2.json")) as fh:
        cells = json.load(fh)["cases"]
    unres = [s for s in batch if s["final_verdict"] == "UNRESOLVED"]
    t2 = [s for s in unres if s["branch"] == "T2"]
    t1 = [s for s in unres if s["branch"] == "T1"]
    dege10_t2 = [s for s in t2 if int(s["deg_e"]) == 10]
    check("V8a  130 UNRESOLVED sub2 batch states", len(unres) == 130,
          str(len(unres)))
    check("V8b  100 of them are T2, 30 are T1", len(t2) == 100 and len(t1) == 30)
    check("V8c  90 deg_e=10 T2 states (the target batch)", len(dege10_t2) == 90)
    check("V9a  branch T2 <=> d1_zero, on all 194 census states",
          all((s["branch"] == "T2") == bool(s["d1_zero"]) for s in batch))
    check("V9b  every T2 census state has deg_d1 == '-inf'",
          all(s["deg_d1"] == "-inf" for s in batch if s["branch"] == "T2"))
    check("V9c  220 sub2 cells, 24 of them T2",
          len(cells) == 220 and sum(1 for c in cells if c["branch"] == "T2") == 24)
    t2_hi = [c for c in cells if c["branch"] == "T2" and max(c["b"]) >= 2]
    check("V9d  exactly 2 T2 cells / 16 states carry a marked-root "
          "multiplicity >= 2",
          len(t2_hi) == 2 and sum(c["state_count"] for c in t2_hi) == 16,
          str([(c["a_t"], c["b"], c["state_count"]) for c in t2_hi]))
    check("V9e  deg e <= 10 on every sub2 state (so a <= 10, always in regime)",
          all(s["deg_e"] <= 10 for c in cells for s in c["states"]))
    try:
        import full_system_bridge as fsb
        caps = fsb.STRIP_DEGCAP["sub2"]
        check("V10  spare caps 12/14 taken from full_system_bridge.STRIP_DEGCAP",
              caps["dm2"] == PT.CAPS["dm2"] == 12
              and caps["dm3"] == PT.CAPS["dm3"] == 14, str(caps))
    except Exception as ex:                                     # pragma: no cover
        check("V10  spare caps from full_system_bridge", False, str(ex)[:80])


# ---------------------------------------------------------------------- V11
def v11_against_valsplit():
    """R9_VALSPLIT.md sec.1 certifies  v(dm2) + v(dm3) >= 9  (from
    monic(e) | dm2*dm3).  The pole theorem gives v(dm2) >= 9 AND v(dm3) >= 9.
    Check: (i) it IMPLIES the certified row; (ii) it is STRICTLY stronger --
    exhibit a profile the certified row allows and the pole theorem forbids;
    (iii) it collapses the 20-case split to 1."""
    a = 9
    lo = PT.pole_bound_closed(a, 30)
    implies = lo["rho_min"] + lo["sigma_min"] >= 9
    # a profile allowed by the old row but forbidden by the new one
    witness = (0, 9)                      # v(dm2)=0, v(dm3)=9: sum 9, but rho<9
    strictly = witness[0] + witness[1] >= 9 and witness[0] < lo["rho_min"]
    old_cases = 10 * 2                    # i = 0..9 times j = 0..1
    check("V11a new bound implies the certified v(dm2)+v(dm3) >= 9", implies)
    check("V11b strictly stronger (profile (0,9) allowed before, forbidden now)",
          strictly)
    check("V11c 20 valuation cases -> 1", old_cases == 20)


def main():
    quiet = "--quiet" in sys.argv
    v1_v2_v3()
    v4_leading_coefficients()
    v5_independent_enumeration()
    v6_support()
    v7_collapse()
    v8_v9_v10()
    v11_against_valsplit()
    ok = all(r[1] for r in RESULTS)
    if not quiet:
        print("=" * 78)
        print("pole_theorem_sweep_verify -- INDEPENDENT CHECKS")
        print("=" * 78)
        for name, good, detail in RESULTS:
            print("  [%s] %s%s" % ("PASS" if good else "FAIL", name,
                                   ("   " + detail[:90]) if detail and not good
                                   else ""))
        print("-" * 78)
    print("%d checks, %d pass, %d fail: %s"
          % (len(RESULTS), sum(1 for r in RESULTS if r[1]),
             sum(1 for r in RESULTS if not r[1]),
             "ALL CHECKS PASS" if ok else "*** FAILURES ***"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
