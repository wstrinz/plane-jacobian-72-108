#!/usr/bin/env python3
"""Galois-descent sweep of the full residue library C01..C23 -- exact checker.

Companion: GALOIS_LIBRARY.md.  Extends GALOIS_DESCENT_PILOT.md (C08) to the
whole 23-shape library and proves the C20 kill through the same
S4/C4/A4/D4-V4 branching.

Source linking: every residue equation is re-derived from f31_graded.txt via
the audited machinery of residue_lemmas_verify.py (load_h / source_terms /
inventory / residue_shapes).  No coefficient in this file is hand-copied from
RESIDUE_LEMMAS.md; the three witness quartics below are synthetic test data
whose every claimed property (irreducibility, Galois type, discriminant
class, membership identities) is verified here from scratch.

WHAT IS PROVED (see GALOIS_LIBRARY.md for the prose):

A. Support-geometry classification.  For each shape, the exponent support is
   either collinear (rank of the difference lattice = 1) or not.  Collinear
   shapes reduce ON THE COEFFICIENT TORUS (all variables nonzero) to a
   one-variable polynomial P(rho) in a primitive monomial character rho;
   primitivity gives an integer vector u with <g,u> = 1, so the character is
   SURJECTIVE on F*-points over every field F and torus solvability is
   exactly "P has a nonzero root in F".  The verifier checks the reduction
   identity and the surjectivity certificate exactly, then:
     - deg P = 1: root -c0/c1 in Q*, explicit torus point built from u and
       substituted back into the source equation => solvable over EVERY
       field => CONSTRAINT for every forcing quartic (12 shapes).
     - deg P = 2: constant discriminant square class Delta.  Delta != 1 =>
       QUADRATIC-OBSTRUCTION: solvable over F iff sqrt(Delta) in F.
       Exactly C08 (Delta=105) and C20 (Delta=170).
     - deg P >= 3: would need a higher resolvent; ZERO such shapes occur.
   Non-collinear (rank >= 2) shapes: rational witness substituted exactly =>
   solvable over every field (9 shapes).
   So the library's kill/constraint split is PREDICTED by support geometry:
   a shape can carry an arithmetic obstruction iff its support is collinear
   with ratio-degree 2, and both such shapes have non-square classes.

B. The C20 kill through the Galois branching (mirroring the C08 pilot):
   - fixed q (S4, disc class 17): Obs = squarefree(170*17) = 10 != 1 => KILL;
   - sharpness: qs = y^4+7y^2-8y+10 is S4 with disc = 170*52^2, so
     Obs(qs) = 170*170 = square => obstruction VANISHES and C20 is solvable
     over split(qs) via r = (3+2*sqrt(170))/11 (checked exactly);
   - branching necessity: qv = y^4-344y^2+28224 (minimal polynomial of
     sqrt(170)+sqrt(2)) is V4 with square discriminant; the naive cyclic
     obstruction 170*1 = 170 would predict a kill, but
     (512x-x^3)^2 = 170*336^2 mod qv exactly, so sqrt(170) in L and C20 is
     solvable there;
   - A4 branch live test: qa = y^4+8y+12 is A4 (order 12, square disc);
     no quadratic subfield => C20 kills over split(qa).
   - joint obstruction: the subgroup <105,170> = {1,105,170,714} of
     Q*/(Q*)^2 meets the fixed q's split-field square-class subgroup {1,17}
     trivially -- C08 and C20 kill simultaneously in the same field.

D4 is decided by explicit membership witness only (same honest boundary as
pilot judgment J4); no D4 case is claimed here.

Run: python galois_library_verify.py [--quiet].  Exit 0 iff all checks pass.
"""

import sys
from math import gcd

import sympy as sp
from sympy.polys.numberfields import galois_group

import residue_lemmas_verify as rl

def _require(_cond, _msg):
    """Proof-critical check: fails loudly and exits nonzero, unaffected by python -O."""
    if not _cond:
        import sys as _sys
        print("FAIL: " + str(_msg))
        _sys.exit(1)


y = sp.Symbol("y")
QUIET = "--quiet" in sys.argv
CHECKS = []


def check(label, condition):
    CHECKS.append((label, bool(condition)))
    if not QUIET:
        print(f"[{'OK' if condition else 'FAIL'}] {label}")
    if not condition:
        raise AssertionError(label)


def squarefree_part(value):
    """Signed squarefree part = the square class of a nonzero rational."""
    value = sp.nsimplify(value)
    num, den = sp.fraction(sp.together(value))
    value = int(num) * int(den)
    _require(value != 0, "value != 0")
    sign = -1 if value < 0 else 1
    return sign * int(sp.prod(p for p, e in sp.factorint(abs(value)).items()
                              if e % 2))


# ---------------------------------------------------------------- source load
def load_library():
    h = rl.load_h()
    terms = rl.source_terms(h)
    rows = {label: rl.survivors(path) for label, path in rl.WINDOWS.items()}
    inventories = {label: rl.inventory(rows[label], terms) for label in rl.WINDOWS}
    ordered, equations = rl.residue_shapes(inventories, terms)
    return equations


# --------------------------------------------- support-geometry classification
def primitive_direction(diffs):
    """The primitive integer direction g with every diff an integer multiple,
    plus the multiples t_i.  Returns None if the diffs are not collinear."""
    nonzero = [d for d in diffs if any(x != 0 for x in d)]
    v = nonzero[0]
    g0 = 0
    for x in v:
        g0 = gcd(g0, abs(int(x)))
    g = tuple(int(x) // g0 for x in v)
    ts = []
    for d in diffs:
        pivot = next((i for i, x in enumerate(g) if x != 0))
        if int(d[pivot]) % g[pivot] != 0:
            return None
        t = int(d[pivot]) // g[pivot]
        if tuple(t * x for x in g) != tuple(int(z) for z in d):
            return None
        ts.append(t)
    return g, ts


def bezout_vector(g):
    """Integer u with <g,u> = 1 (exists iff gcd(g) = 1)."""
    # iterative extended gcd across the 4 coordinates
    u = [0, 0, 0, 0]
    cur_g, cur_u = 0, [0, 0, 0, 0]
    for i, gi in enumerate(g):
        old = cur_g
        cur_g_new = gcd(cur_g, abs(gi))
        if cur_g_new == 0:
            continue
        # find a,b with a*old + b*gi = cur_g_new
        a, b = sp.gcdex(old, gi)[:2] if old != 0 else (0, sp.sign(gi))
        cur_u = [int(a) * x for x in cur_u]
        cur_u[i] += int(b)
        cur_g = cur_g_new
    _require(sum(ui * gi for ui, gi in zip(cur_u, g)) == cur_g, "sum(ui * gi for ui, gi in zip(cur_u, g)) == cur_g")
    _require(cur_g == 1, "direction must be primitive")
    return cur_u


def classify(equations):
    """Classify all shapes; verify every reduction/solvability claim exactly."""
    D, X, S, E = rl.LEADS
    gens = (D, X, S, E)
    out = {}
    witnesses = rl.rational_witnesses()
    for name in sorted(equations):
        eq = equations[name]
        P = sp.Poly(eq, *gens, domain=sp.QQ)
        monoms = [m for m, _ in P.terms()]
        coeffs = [c for _, c in P.terms()]
        m0 = monoms[0]
        diffs = [tuple(a - b for a, b in zip(m, m0)) for m in monoms[1:]]
        line = primitive_direction(diffs) if diffs else None
        if line is not None:
            g, ts = line
            ts = [0] + ts
            shift = min(ts)
            ts = [t - shift for t in ts]
            deg = max(ts)
            rho = sp.Symbol("rho")
            Prho = sum(c * rho**t for c, t in zip(coeffs, ts))
            # exact reduction identity on the torus:
            #   eq / monom(m0) == rho^shift * P(rho) with rho = chi_g
            chi = sp.prod(v**e for v, e in zip(gens, g))
            lead = sp.prod(v**e for v, e in zip(gens, m0))
            reduction = sp.cancel(eq / lead - chi**shift * Prho.subs(rho, chi))
            check(f"{name}: collinear support reduces exactly to P(chi_g), "
                  f"deg {deg}", reduction == 0)
            # surjectivity certificate: <g,u> = 1
            u = bezout_vector(g)
            tau = sp.Symbol("tau", nonzero=True)
            point = {v: tau**ui for v, ui in zip(gens, u)}
            chi_val = sp.cancel(chi.subs(point))
            check(f"{name}: character surjectivity certificate <g,u>=1",
                  chi_val == tau)
            if deg == 1:
                lam = sp.Rational(-coeffs[ts.index(0)], coeffs[ts.index(1)])
                pt = {v: lam**ui for v, ui in zip(gens, u)}
                check(f"{name}: LINEAR -- explicit rational torus point "
                      f"(root {lam})",
                      lam != 0 and sp.cancel(eq.subs(pt)) == 0)
                out[name] = ("LINEAR", lam)
            elif deg == 2:
                disc = sp.Poly(Prho, rho).discriminant()
                delta = squarefree_part(disc)
                const_nonzero = Prho.subs(rho, 0) != 0
                check(f"{name}: QUADRATIC in chi_g -- constant residue "
                      f"discriminant class {delta}, roots nonzero",
                      const_nonzero)
                out[name] = ("QUADRATIC", delta)
            else:
                out[name] = ("HIGHER", deg)
        else:
            point = witnesses[name]
            used = eq.free_symbols
            check(f"{name}: MULTIDIM -- audited rational witness re-verified",
                  used <= set(point) and all(point[v] != 0 for v in used)
                  and sp.expand(eq.subs(point)) == 0)
            out[name] = ("MULTIDIM", None)
    return out


# --------------------------------------------------- Galois branching machine
def galois_label(q_poly):
    P = sp.Poly(q_poly, y, domain=sp.QQ)
    _require(P.is_irreducible, "P.is_irreducible")
    group, alt = galois_group(q_poly, y)
    order = group.order()
    disc_square = squarefree_part(sp.discriminant(q_poly, y)) == 1
    if order == 24:
        return "S4"
    if order == 12:
        return "A4"
    if order == 8:
        return "D4"
    if order == 4:
        return "V4" if disc_square else "C4"
    raise ValueError(order)


def decide_kill(q_poly, delta, membership_witness=None):
    """Master criterion: KILL <=> sqrt(delta) not in split(q).
    S4/C4: Obs = delta*disc class test.  A4: always kills (delta != 1).
    D4/V4: decided only by an exact membership witness w with
    w^2 == delta mod q (witness present => solvable => no kill); V4 with no
    witness supplied returns None (undecided here), same for D4."""
    _require(delta != 1, "delta != 1")
    label = galois_label(q_poly)
    disc_class = squarefree_part(sp.discriminant(q_poly, y))
    if label in ("S4", "C4"):
        obs = squarefree_part(delta * disc_class)
        return {"label": label, "kills": obs != 1, "obs": obs}
    if label == "A4":
        return {"label": label, "kills": True, "obs": None}
    if membership_witness is not None:
        rem = sp.rem(sp.expand(membership_witness**2 - delta), q_poly, y)
        _require(sp.expand(rem) == 0, "membership witness failed")
        return {"label": label, "kills": False, "obs": None}
    return {"label": label, "kills": None, "obs": None}


# -------------------------------------------------------------------- checks
def main():
    equations = load_library()
    check("23 source-derived residue shapes rebuilt from f31_graded.txt",
          len(equations) == 23)

    out = classify(equations)
    quad = {k: v[1] for k, v in out.items() if v[0] == "QUADRATIC"}
    linear = [k for k, v in out.items() if v[0] == "LINEAR"]
    multi = [k for k, v in out.items() if v[0] == "MULTIDIM"]
    higher = [k for k, v in out.items() if v[0] == "HIGHER"]

    check("census: exactly C08 and C20 are quadratic-obstruction shapes",
          quad == {"C08": 105, "C20": 170})
    check("census: 12 linear-in-character shapes (solvable over every field)",
          len(linear) == 12)
    check("census: 9 multidimensional shapes with rational witnesses",
          len(multi) == 9)
    check("census: ZERO shapes need a higher resolvent (no deg>=3 collinear)",
          higher == [])

    # ---- C20 (and C08 regression) through the Galois branching
    q_fixed = 2048*y**4 - 512*y**3 + 320*y**2 - 240*y + 195
    check("fixed q is S4 with discriminant class 17",
          galois_label(q_fixed) == "S4"
          and squarefree_part(sp.discriminant(q_fixed, y)) == 17)

    d08 = decide_kill(q_fixed, 105)
    check("C08 fixed-q verdict: Obs class 1785, KILL (pilot regression)",
          d08["kills"] and d08["obs"] == 1785)
    d20 = decide_kill(q_fixed, 170)
    check("C20 fixed-q verdict: Obs class 10, KILL",
          d20["kills"] and d20["obs"] == 10)

    # joint obstruction subgroup <105,170> vs split-field classes {1,17}
    joint = {1, 105, 170, squarefree_part(105 * 170)}
    check("joint obstruction subgroup {1,105,170,714} misses {1,17} "
          "away from 1", squarefree_part(105 * 170) == 714
          and joint & {17} == set())

    # ---- sharpness: S4 quartic with disc class 170 -> obstruction vanishes
    qs = y**4 + 7*y**2 - 8*y + 10
    disc_qs = sp.discriminant(qs, y)
    check("qs = y^4+7y^2-8y+10 is S4 with disc = 170*52^2 exactly",
          galois_label(qs) == "S4" and disc_qs == 170 * 52**2)
    ds = decide_kill(qs, 170)
    check("C20 obstruction VANISHES for qs (Obs = 170*170 = square)",
          ds["kills"] is False and ds["obs"] == 1)
    # constructive: sqrt(170) = Vandermonde/52 in split(qs); numeric confirm
    roots = sp.Poly(qs, y).all_roots()
    vand = sp.prod(roots[j] - roots[i] for i in range(4) for j in range(i+1, 4))
    check("Vandermonde(qs)^2 = disc(qs) numerically (so sqrt170 in L)",
          abs(complex(sp.N(vand**2, 30)) - float(disc_qs)) < 1e-6)
    r20 = (3 + 2*sp.sqrt(170)) / 11
    c20_eq = equations["C20"]
    D, X, S, E = rl.LEADS
    check("r = (3+2*sqrt170)/11 solves source C20 with X=D=1 exactly",
          sp.simplify(c20_eq.subs({D: 1, X: 1, E: r20})) == 0 and r20 != 0)

    # ---- branching necessity: V4 witness containing sqrt(170)
    qv = y**4 - 344*y**2 + 28224
    check("qv = minpoly(sqrt170+sqrt2) is irreducible V4 with square disc",
          galois_label(qv) == "V4")
    naive = squarefree_part(170 * sp.discriminant(qv, y))
    check("naive cyclic-regime obstruction for qv is 170 (would claim kill)",
          naive == 170)
    dv = decide_kill(qv, 170, membership_witness=(512*y - y**3) / 336)
    check("exact witness ((512x-x^3)/336)^2 = 170 mod qv => sqrt170 in L, "
          "NO kill: branching is necessary", dv["kills"] is False)

    # ---- A4 branch live test
    qa = y**4 + 8*y + 12
    check("qa = y^4+8y+12 is A4 (square disc, order 12)",
          galois_label(qa) == "A4")
    da = decide_kill(qa, 170)
    check("C20 kills over split(qa): A4 has no quadratic subfield",
          da["kills"] is True)

    n = len(CHECKS)
    print(f"\nALL {n} GALOIS-LIBRARY CHECKS PASSED")


if __name__ == "__main__":
    try:
        main()
    except AssertionError as exc:
        print(f"\nFAILED: {exc}")
        sys.exit(1)
    sys.exit(0)
