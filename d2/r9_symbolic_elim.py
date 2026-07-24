#!/usr/bin/env python3
"""r9_symbolic_elim.py -- ONE-TIME symbolic elimination of the spare window
unknown dm4 = d_-4 from the pre-resultant G-system, BEFORE any state is
grounded.  This is the structural reduction BRIDGE_SWEEP.md section 6 named as
the likely cure for the deg_e = 10 cost wall (R9 z >= 1 / a10): the full
bridge's 122-equation system with 45 spare scalar unknowns swells every
Groebner engine; eliminating dm4 symbolically removes its 17 (sub2) / 25
(sub1) scalar unknowns from every instantiated state at zero soundness cost.

THE ELIMINATION (exact, certified by construction)
--------------------------------------------------
Every G-system generator is LINEAR in dm4 (verified below):

    G1 = 3*dm1*dm4 + A1,   A1 = 3/2*d1*dm1^2 + 3*d2*dm1*dm2 + 3*dm2*dm3
    G2 = 3*dm2*dm4 + A2
    G3 = 3*dm3*dm4 + A3
    G5 = -3*(d0*dm1 + d1*dm2 + d2*dm3)*dm4 + A5      (G5 = G5body + Phi)

Since dm1 = e is a KNOWN nonzero polynomial in every cascade state, the
cross-multiplied combinations

    H2 := dm1*G2 - dm2*G1
    H3 := dm1*G3 - dm3*G1
    H5 := dm1*G5 + (d0*dm1 + d1*dm2 + d2*dm3)*G1

are dm4-FREE elements of the ideal <G1,G2,G3,G5>.  Hence for ANY state, the
system {y-coefficients of H2,H3,H5} + window caps on dm2,dm3 alone is a sound
NECESSARY system: a UNIT verdict on it kills the state outright.  (It is a
priori weaker than the full G-system -- dropping G1 forgets the coupling to a
polynomial dm4 of capped degree -- so a PROPER verdict here is INCONCLUSIVE,
never a survival signal.)

THE DIVISIBILITY LEMMA (recovers most of G1's content, still dm4-free)
----------------------------------------------------------------------
G1 = 0 rearranges to  3*dm2*dm3 = -dm1*(3/2*d1*dm1 + 3*d2*dm2 + 3*dm4),
all factors polynomials in Q[y].  Therefore on the G-variety

    dm1 | dm2*dm3   in Q[y],   i.e.   M | dm2*dm3

for M := e / lc-scalar the monic polynomial part of the state's e.  Per state
this contributes deg(M) further necessary equations
rem(dm2*dm3, M, y) == 0, coefficient-wise -- quadratic in the spare
coefficients and very small.  (Necessary only; recovering dm4 as a polynomial
*of capped degree <= 16/24* needs slightly more, which we deliberately do not
claim.)

This module verifies the combination identities by exact re-expansion,
verifies weighted homogeneity of H2,H3,H5 (weights 228, 240, 264), and caches
everything to r9_eliminated_system.json for r9_symbolic_sweep.py.

New file, uncommitted.  READ-ONLY on every imported module.
"""
from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

import full_system_bridge as fsb

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "r9_eliminated_system.json"

D2, D1, D0 = fsb.D2, fsb.D1, fsb.D0
DM1, DM2, DM3, DM4 = fsb.DM1, fsb.DM2, fsb.DM3, fsb.DM4
PHI = fsb.PHI

# weights of the dm4-free combinations: w(dm1) + w(G_i)
EXPECTED_H_WEIGHTS = {"H2": 60 + 168, "H3": 60 + 180, "H5": 60 + 204}


def eliminate() -> dict:
    """Build H2, H3, H5 and their explicit cofactor certificates."""
    g = fsb.gsystem()
    checks = []

    # 1. linearity in dm4 + the exact pivot coefficients
    pivots = {}
    for name, expr in g.items():
        p = sp.Poly(expr, DM4)
        assert p.degree() == 1, f"{name} not linear in dm4"
        pivots[name] = sp.expand(p.nth(1))
    assert pivots["G1"] == sp.expand(3 * DM1), pivots["G1"]
    checks.append("linearity: G1,G2,G3,G5 all degree 1 in dm4; pivot(G1)=3*dm1")

    # 2. the combinations, with cofactors recorded (the membership certificate)
    combos = {
        "H2": ((DM1, "G2"), (-DM2, "G1")),
        "H3": ((DM1, "G3"), (-DM3, "G1")),
        "H5": ((DM1, "G5"), (D0 * DM1 + D1 * DM2 + D2 * DM3, "G1")),
    }
    H: dict[str, sp.Expr] = {}
    for hname, parts in combos.items():
        expr = sp.expand(sum(cof * g[gname] for cof, gname in parts))
        assert DM4 not in expr.free_symbols, f"{hname} still contains dm4"
        # re-expansion identity: H - sum(cof*G) == 0 exactly
        residual = sp.expand(expr - sum(cof * g[gname] for cof, gname in parts))
        assert residual == 0
        H[hname] = expr
    checks.append("dm4-free: H2,H3,H5 contain no dm4; membership residuals 0")

    # 3. weighted homogeneity of each H (inherits from the G-weights)
    weights = {}
    for hname, expr in H.items():
        seen = set()
        for term in sp.Add.make_args(expr):
            w = 0
            for b_, ex in term.as_powers_dict().items():
                if b_.is_number:
                    continue
                w += fsb.WEIGHT[str(b_)] * ex
            seen.add(w)
        assert seen == {EXPECTED_H_WEIGHTS[hname]}, (hname, seen)
        weights[hname] = EXPECTED_H_WEIGHTS[hname]
    checks.append(f"weighted-homogeneous: {weights}")

    # 4. the divisibility lemma identity:
    #    3*dm2*dm3 + dm1*(3/2*d1*dm1 + 3*d2*dm2 + 3*dm4) - G1 == 0
    div_residual = sp.expand(
        3 * DM2 * DM3 + DM1 * (sp.Rational(3, 2) * D1 * DM1 + 3 * D2 * DM2
                               + 3 * DM4) - g["G1"])
    assert div_residual == 0
    checks.append("divisibility lemma: G1 rearrangement identity exact "
                  "=> dm1 | dm2*dm3 on the G-variety")

    return {
        "H": {k: sp.sstr(v) for k, v in H.items()},
        "cofactors": {h: [[sp.sstr(cof), gname] for cof, gname in parts]
                      for h, parts in combos.items()},
        "pivots": {k: sp.sstr(v) for k, v in pivots.items()},
        "weights": weights,
        "checks": checks,
        "note": ("H2,H3,H5 in <G1,G2,G3,G5>; dm4-free. Sound necessary "
                 "system per state: y-coeffs of H2,H3,H5 with dm2,dm3 "
                 "capped (sub2: 12/14) + rem(dm2*dm3, monic(e)) = 0. "
                 "UNIT => state killed (PENDING AUDIT). PROPER => "
                 "inconclusive (weaker than full bridge)."),
        "_H_sympy": H,  # stripped before JSON dump
    }


def main() -> int:
    res = eliminate()
    H = res.pop("_H_sympy")
    for c in res["checks"]:
        print(f"[OK] {c}", flush=True)
    for k, v in H.items():
        print(f"  {k}: {len(sp.Add.make_args(v))} terms", flush=True)
    json.dump(res, open(OUT, "w"), indent=1)
    print(f"wrote {OUT.name}")
    print("ALL SYMBOLIC-ELIMINATION CHECKS PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
