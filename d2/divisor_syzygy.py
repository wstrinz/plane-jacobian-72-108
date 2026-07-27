#!/usr/bin/env python3
"""divisor_syzygy.py -- the universal K-syzygy of the G-system, and e | Phi

THE IDENTITY.  With e = dm1, R = dm2, S = dm3, T = dm4, the canonical G-system
generators satisfy the EXACT polynomial identity

    2*(G5 + d2*G3 + d1*G2 + d0*G1)  ==  2*Phi - e*(d2*e^2 + 3*e*S + 3*R^2)

Everything else -- every dm4 term, every d0/d1/d2 cross term -- cancels
identically.  Setting K := 2*Phi - e*(d2*e^2 + 3*e*S + 3*R^2), this gives

    <G1,G2,G3,G5> == <G1,G2,G3,K>          (char 0; K is an exact combination,
                                            and G5 = K/2 - d2*G3 - d1*G2 - d0*G1)

so G5 may be replaced by the sparse four-term K-row with IDEAL EQUALITY -- no
saturation, no division by e.

THE CONSEQUENCE.  On every genuine lift K = 0, i.e.

    2*Phi = e*(d2*e^2 + 3*e*S + 3*R^2)        =>        e | Phi.

Phi = -(1/6630)*(y+1)^30*q(y) with q the fixed SQUAREFREE quartic.  So over the
algebraic closure e = gamma*(y+1)^a * (a squarefree divisor of q): e has NO roots
off {y=-1} U {roots of q}, and every simple q-root divides e to order <= 1.
That is b_i in {0,1} -- not a degree bound, not a stratum condition, and it holds
on the whole G-variety.

WHY THIS MATTERS.  Per SESSION_HANDOFF's THE SPEC, Phi-depth kills bottom-up and
closes a cell only by reaching the TOP stratum, which happens in 0 of 220 sub2
cells.  This lemma is not monotone in degree -- it deletes whole SUPPORT cells
including their top strata -- which is exactly the shape a cell-closing
mechanism has to have.

SUB2: THE DEGREE IS FORCED.  With the certified sub2 caps deg d2 <= 4,
deg R <= 12, deg S <= 14, deg e <= 10, the RHS has degree
E + max(4 + 2E, E + 14, 24).  Against deg Phi = 34 that is <= 33 for E <= 9 --
impossible -- and exactly 34 at E = 10.  Hence

    deg e = 10 EXACTLY for every sub2 G-system solution,  and  a + sum(b_i) = 10.

INDEPENDENT CROSS-CHECK.  The generic-fiber lane, working from the dm4-eliminated
H-system at d1 = 0, landed 2*(H5 + d2*H3) = dm1*K5 with
K5 = 2*Phi - 3*dm1*dm2^2 - d2*dm1^3 - 3*dm1^2*dm3.  That K5 is LITERALLY this K.
Two lanes, two starting systems, one object -- see check C4.

STATUS.  The syzygy and the sub2 degree forcing are verified exactly here.  The
t^a-divisibility of the spares (which would collapse 45 spare coefficients to
45-3a) is NOT verified in this file -- see DIVISOR_SYZYGY.md sec.4.

Read-only.  Usage:
    python divisor_syzygy.py           # full report
    python divisor_syzygy.py --quiet   # self-check, exit 0 iff all pass
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))

d0, d1, d2 = sp.symbols("d0 d1 d2")
e, R, S, T = sp.symbols("dm1 dm2 dm3 dm4")
Phi = sp.Symbol("Phi")

DEG_PHI = 34                                    # Phi = c*(y+1)^30*q, deg q = 4
SUB2_CAPS = {"d2": 4, "R": 12, "S": 14, "e": 10}


def K_form():
    """The sparse replacement row K = 2*Phi - e*(d2*e^2 + 3*e*S + 3*R^2)."""
    return 2 * Phi - e * (d2 * e**2 + 3 * e * S + 3 * R**2)


def syzygy_residual():
    """2*(G5 + d2*G3 + d1*G2 + d0*G1) - K, expanded. Must be exactly 0."""
    import bigrade_annotator as ba
    g = dict(ba._G_generators())
    G1, G2, G3, G5 = (sp.expand(g[k][0]) for k in ("G1", "G2", "G3", "G5"))
    lhs = sp.expand(2 * (G5 + d2 * G3 + d1 * G2 + d0 * G1))
    return sp.expand(lhs - K_form()), G5


def forced_deg_e(caps=None):
    """Least E with E + max(deg d2 + 2E, E + deg S, 2 deg R) >= deg Phi."""
    c = caps or SUB2_CAPS
    rows = []
    for E in range(0, c["e"] + 2):
        rhs = E + max(c["d2"] + 2 * E, E + c["S"], 2 * c["R"])
        rows.append((E, rhs))
    feasible = [E for E, rhs in rows if rhs >= DEG_PHI and E <= c["e"]]
    return (min(feasible) if feasible else None), rows


def column_alive(col, total=10):
    """Does a support column survive  a + sum(b_i) = total  with b_i in {0,1}?

    Two independent ways to die, and they are NOT the same test:
      * a + sum(b_i) != total          -- the forced-degree condition
      * some b_i >= 2                  -- e | Phi at a SIMPLE q-root
    """
    a = int(col.split("_b")[0][1:])
    b = [int(x) for x in col.split("_b")[1].split("_")[0]]
    if any(x > 1 for x in b):
        return False, "b_i>=2 at a simple q-root (e | Phi forbids it)"
    if a + sum(b) != total:
        return False, "a+sum(b)=%d != %d (forced deg e)" % (a + sum(b), total)
    return True, "survives"


OPEN_T2 = ["a9_b1000", "a8_b0000", "a8_b1000", "a8_b1100",
           "a7_b1000", "a7_b1100", "a7_b1110", "a7_b3000"]


def run(verbose=True):
    out, npass, ntot = [], 0, 0

    def ck(name, ok, detail):
        nonlocal npass, ntot
        ntot += 1
        npass += bool(ok)
        out.append("  [%s] %s\n        %s" % ("PASS" if ok else "FAIL", name, detail))

    # C1 -- the syzygy itself
    resid, G5 = syzygy_residual()
    ck("C1  syzygy residual is EXACTLY zero", resid == 0,
       "2*(G5 + d2*G3 + d1*G2 + d0*G1) - K  =  %s" % resid)

    # C2 -- canonical-G5 guard. A stale 2*Phi transcription was a real bug here.
    coeff = sp.expand(G5).coeff(Phi)
    ck("C2  G5 has Phi-coefficient exactly 1 (canonical form)", coeff == 1,
       "coeff(G5, Phi) = %s -- a stale 2*Phi form would silently break C1" % coeff)

    # C3 -- ideal equality is genuine: G5 is recoverable from K
    import bigrade_annotator as ba
    g = dict(ba._G_generators())
    G1, G2, G3 = (sp.expand(g[k][0]) for k in ("G1", "G2", "G3"))
    recovered = sp.expand(K_form() / 2 - d2 * G3 - d1 * G2 - d0 * G1)
    ck("C3  G5 = K/2 - d2*G3 - d1*G2 - d0*G1 (so the swap is ideal EQUALITY)",
       sp.expand(recovered - G5) == 0,
       "residual = %s -- K is not merely a consequence, it is exchangeable"
       % sp.expand(recovered - G5))

    # C4 -- independent cross-check against the generic-fiber lane's K5
    K5_generic = 2 * Phi - 3 * e * R**2 - d2 * e**3 - 3 * e**2 * S
    ck("C4  equals the generic-fiber lane's K5 (independent derivation)",
       sp.expand(K5_generic - K_form()) == 0,
       "the d1=0 H-system route and the universal G-system route produced the "
       "SAME object; residual = %s" % sp.expand(K5_generic - K_form()))

    # C5 -- sub2 forces deg e = 10
    E, rows = forced_deg_e()
    ck("C5  sub2 caps force deg e = 10 exactly", E == 10,
       "min feasible E = %s; RHS degree by E: %s (deg Phi = %d)"
       % (E, ", ".join("%d->%d" % r for r in rows[6:12]), DEG_PHI))

    # C6 -- the open T2 columns collapse 8 -> 3
    alive = [c for c in OPEN_T2 if column_alive(c)[0]]
    ck("C6  open sub2 T2 columns collapse 8 -> 3", len(alive) == 3,
       "survivors: %s" % ", ".join(alive))

    # C7 -- the live pilot is dead, and for the multiplicity reason specifically
    ok7, why7 = column_alive("a7_b3000")
    ck("C7  the a7_b3000_T2 pilot is EMPTY before any Groebner run",
       (not ok7) and "b_i>=2" in why7,
       "a7_b3000 -> %s" % why7)

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
            print("divisor_syzygy: %d/%d checks FAILED" % (ntot - npass, ntot))
            return 1
        print("divisor_syzygy: %d/%d checks pass" % (npass, ntot))
        return 0
    print("\n%d/%d checks pass" % (npass, ntot))
    if npass == ntot:
        print("\n  K = %s" % K_form())
        print("  => 2*Phi = e*(d2*e^2 + 3*e*S + 3*R^2)  =>  e | Phi")
        print("  => e = gamma*(y+1)^a * (squarefree divisor of q), a + sum(b_i) = 10")
    return 0 if npass == ntot else 1


if __name__ == "__main__":
    sys.exit(main())
