#!/usr/bin/env python3
"""Exact source-linked checks for RESIDUE_LEMMAS_DEPTH.md.

Extends ``residue_lemmas_verify.py`` from the depth-one initial forms to the
a-quantified *jet* structure on the surviving f31 subcase-(1) frontier
``cascade_cones_sub1_qt_inf_rl.json`` (171 branches / 1145 flag cases,
a in [2,10]).

The mathematical content proved here (all re-derived from ``f31_graded.txt``;
no h_l coefficient is copied into this file):

  1. Census.  The 1145 surviving flag cases are inventoried by
     (place, level, kind, tied-support).  Every t-place tied monomial string
     is checked to be an exact term of the source h_l (0 mismatches).  The
     affine t-depth law ``depth = 30 - 3a`` is verified against the survivor
     term_cancellation depths.

  2. Jet decomposition.  For a tied set T with common valuation-weight m,
     substituting local jet series d2=D+D1*pi+D2*pi^2, ... into the tie
     polynomial F=sum(T) and factoring pi^m gives F = C0 + C1*pi + C2*pi^2+...
     A depth-delta obligation is exactly C0=...=C_{delta-1}=0.  We prove the
     exact identities
         C0 = IF                       (the depth-1 initial form),
         C1 = grad(IF) . jet1,
         C2 = grad(IF) . jet2 + (1/2) jet1^T Hess(IF) jet1,
     for the TOP THREE highest-incidence a-growing t-place window patterns
     (C09 @ L5, C02 @ L6, C22 @ L4).

  3. Verdicts.  Because the t-place y=-1 is a *finite* place, the residues
     (D,X,S,E) and every jet are free rational Taylor coefficients.  At a
     rational nonzero-leading point where grad(IF) != 0 (a smooth point) the
     jet tower C1=0, C2=0, ... is solvable to every order: each C_delta is
     linear in the order-delta jet with the fixed nonzero coefficient
     grad(IF).  Hence the extended (depth-1 hypersurface + all jets) system is
     a CONSTRAINT for EVERY depth, i.e. for every a in [2,10] at once.  We
     exhibit an explicit exact depth-3 solution for each of the three, then
     prove the general no-kill theorem: every t-place monomial_tie_rise
     support in the frontier is smooth at a rational nonzero-leading point.

  4. Kill accounting.  The two program kills C08/C20 are q-place, depth-1 and
     a-independent (RESIDUE_LEMMAS.md).  We confirm their supports do not
     occur anywhere on the 1145 survivors (their carrying branches were
     already removed upstream by the residue-kill feedback), so the depth
     extension introduces no new kill and removes no additional survivor.

Run:  python residue_lemmas_depth_verify.py
"""

from __future__ import annotations

import json
from collections import Counter, defaultdict
from pathlib import Path

import sympy as sp

from cascade_signature import load_levels

def _require(_cond, _msg):
    """Proof-critical check: fails loudly and exits nonzero, unaffected by python -O."""
    if not _cond:
        import sys as _sys
        print("FAIL: " + str(_msg))
        _sys.exit(1)


ROOT = Path(__file__).resolve().parent
FRONTIER = ROOT / "cascade_cones_sub1_qt_inf_rl.json"

# Signature variables and their leading-coefficient / jet symbols.
d2, d1, sigma, e = sp.symbols("d2 d1 sigma e")
SIGVARS = (d2, d1, sigma, e)
D, X, S, E = sp.symbols("D X S E")
LEAD = (D, X, S, E)
J1 = sp.symbols("D1 X1 S1 E1")          # order-1 jets
J2 = sp.symbols("D2 X2 S2 E2")          # order-2 jets
pi = sp.Symbol("pi")
LEAD_OF = dict(zip(SIGVARS, LEAD))
VALW = {d2: 2, d1: 3, sigma: 4, e: 5}   # signature weights (info only)


# --------------------------------------------------------------------------- #
# Source parsing (no coefficient is copied; all pulled from f31_graded.txt).
# --------------------------------------------------------------------------- #
def source_h():
    lv = load_levels()
    return {l: sp.expand(lv[l].sigma_expression) for l in (4, 5, 6)}


def source_terms(h):
    return {l: dict(sp.Poly(h[l], *SIGVARS, domain=sp.QQ).terms())
            for l in (4, 5, 6)}


def tie_polynomial(level, support, terms):
    """Rebuild sum_{M in support} coeff_source(M) * M from the parsed source."""
    F = sp.Integer(0)
    for mono in support:
        _require(mono in terms[level], "support monomial absent from source "
                                      f"h_{level}: {mono}")
        F += terms[level][mono] * sp.prod(v**n for v, n in zip(SIGVARS, mono))
    return sp.expand(F)


def initial_form(F):
    """IF(D,X,S,E): replace each source variable by its leading coefficient."""
    return sp.expand(F.subs(LEAD_OF))


# --------------------------------------------------------------------------- #
# Jet expansion.
# --------------------------------------------------------------------------- #
def jet_coefficients(F):
    """Return [C0, C1, C2]: coefficients of pi^0,pi^1,pi^2 of F(jet series).

    Valid because every monomial of a *tied* support shares one weight m, so
    the pi^m prefactor divides out uniformly and the bare series substitution
    reproduces the convolution coefficients exactly.
    """
    series = {v: LEAD_OF[v] + J1[i] * pi + J2[i] * pi**2
              for i, v in enumerate(SIGVARS)}
    P = sp.Poly(sp.expand(F.subs(series)), pi)
    return [sp.expand(P.coeff_monomial(pi**j)) for j in range(3)]


def grad_dot(IF, jets):
    return sp.expand(sum(sp.diff(IF, v) * j for v, j in zip(LEAD, jets)))


def hess_quadratic(IF, jets):
    return sp.expand(sp.Rational(1, 2) * sum(
        sp.diff(IF, a, b) * ja * jb
        for a, ja in zip(LEAD, jets) for b, jb in zip(LEAD, jets)))


def tie_weight_exists(support):
    """A nonnegative integer weight (k,x,z,b) equating all support weights?"""
    k, x, z, b = sp.symbols("k x z b", integer=True, nonnegative=True)
    wv = {d2: k, d1: x, sigma: z, e: b}
    weights = [sum(wv[v] * n for v, n in zip(SIGVARS, mono)) for mono in support]
    eqs = [sp.Eq(weights[0], w) for w in weights[1:]]
    return sp.solve(eqs, [k, x, z, b], dict=True) != [] or len(support) == 1


# --------------------------------------------------------------------------- #
# Frontier census.
# --------------------------------------------------------------------------- #
def survivors():
    data = json.loads(FRONTIER.read_text(encoding="utf-8"))
    surv = [b for b in data["branches"] if b.get("survivor_cases")]
    return data, surv


def cell_id(b):
    return "a%d %s %s" % (b["a_t"], "".join(map(str, b["b"])), b["branch"])


def parse_support(tied):
    support = []
    for t in tied:
        ex = sp.sympify(t.replace("^", "**"),
                        locals={"d2": d2, "d1": d1, "sigma": sigma, "e": e})
        (mono,) = sp.Poly(ex, *SIGVARS).monoms()
        support.append((mono, sp.Poly(ex, *SIGVARS).coeffs()[0]))
    return support


def census(surv, terms):
    """(place,level,kind,support) -> freq/cells/avals/depths, with source check."""
    freq = Counter()
    cells = defaultdict(set)
    avals = defaultdict(Counter)
    depths = defaultdict(Counter)
    mism = 0
    for b in surv:
        cid = cell_id(b)
        for case in b["survivor_cases"]:
            for w in case["witness"]:
                # Only the q/t places carry monomial residue ties; the inf
                # place records bracketed leading-cancellation polynomials.
                if w["place"] not in ("q", "t"):
                    continue
                for ob in w["obligations"]:
                    support = []
                    tied = (ob.get("tied", [])
                            if ob["kind"] in ("monomial_tie_rise",
                                              "identical_vanishing") else [])
                    for mono, coeff in parse_support(tied):
                        if (mono not in terms.get(ob["level"], {})
                                or terms[ob["level"]][mono] != coeff):
                            mism += 1
                        support.append(mono)
                    key = (w["place"], ob["level"], ob["kind"],
                           frozenset(support))
                    freq[key] += 1
                    cells[key].add(cid)
                    avals[key][b["a_t"]] += 1
                    depths[key][ob["depth"]] += 1
    return freq, cells, avals, depths, mism


# --------------------------------------------------------------------------- #
# Smooth rational point finder (for the general no-kill theorem).
# --------------------------------------------------------------------------- #
# Rational nonzero-leading certificates (RESIDUE_LEMMAS.md section 4), keyed by
# the tied support.  Each point is VERIFIED below against the source-derived
# initial form -- the table supplies only candidate certificates.
r73 = sp.Rational(73, 4)
u22 = sp.Rational(152, 511)
WIT_BY_SUPPORT = {
    frozenset({(0, 0, 2, 0), (0, 1, 0, 1)}): {X: -6, S: -4, E: -1},          # C01
    frozenset({(0, 0, 2, 0), (1, 2, 0, 0)}): {D: sp.Rational(3, 14), X: 1, S: 1},  # C02
    frozenset({(1, 2, 0, 0), (0, 1, 0, 1)}): {D: 1, X: 1, E: sp.Rational(-7, 4)},  # C03
    frozenset({(2, 2, 0, 0), (1, 0, 2, 0), (0, 2, 1, 0)}):
        {D: sp.Rational(13, 6), X: sp.Rational(13, 6), S: sp.Rational(169, 36)},   # C09
    frozenset({(0, 0, 3, 0), (0, 1, 1, 1)}): {X: 1, S: 10, E: 35},           # C11
    frozenset({(0, 0, 3, 0), (0, 4, 0, 0)}): {X: r73**2, S: -r73**3},        # C12
    frozenset({(0, 0, 3, 0), (2, 0, 2, 0)}): {D: 7, S: -21},                 # C14
    frozenset({(0, 4, 0, 0), (3, 2, 0, 0)}):
        {D: -sp.Rational(13797, 1952), X: sp.Rational(13797, 1952)},         # C16
    frozenset({(0, 0, 3, 0), (0, 4, 0, 0), (0, 1, 1, 1)}):
        {X: 1, S: 1, E: sp.Rational(539, 80)},                               # C18
    frozenset({(0, 0, 3, 0), (0, 4, 0, 0), (1, 2, 1, 0), (3, 2, 0, 0),
               (2, 0, 2, 0)}):
        {D: 1 / u22, X: 1 / u22, S: sp.Rational(-4, 3) / u22**2},            # C22
}


def smooth_certificate(support, IF):
    """Verify the tabled point kills the source IF and IF is smooth there."""
    pt = WIT_BY_SUPPORT.get(support)
    if pt is None:
        return None
    present = sorted(IF.free_symbols, key=sp.default_sort_key)
    if any(pt.get(v, 0) == 0 for v in present):
        return None
    if sp.expand(IF.subs(pt)) != 0:
        return None
    g = {v: sp.expand(sp.diff(IF, v).subs(pt)) for v in present}
    return (pt, g) if any(g[v] != 0 for v in present) else None


# --------------------------------------------------------------------------- #
# The three top a-growing t-place window patterns.
# --------------------------------------------------------------------------- #
# supports as exponent tuples (d2,d1,sigma,e); identified with RESIDUE_LEMMAS ids
SUPP_A = frozenset({(2, 2, 0, 0), (1, 0, 2, 0), (0, 2, 1, 0)})   # L5  = C09
SUPP_B = frozenset({(0, 0, 2, 0), (1, 2, 0, 0)})                 # L6  = C02
SUPP_C = frozenset({(0, 0, 3, 0), (0, 4, 0, 0), (1, 2, 1, 0),
                    (3, 2, 0, 0), (2, 0, 2, 0)})                 # L4  = C22
TOP3 = [("C09", 5, SUPP_A, 201), ("C02", 6, SUPP_B, 196),
        ("C22", 4, SUPP_C, 160)]

# RESIDUE_LEMMAS.md primitive relations (recomputed, only used as a cross-check
# target up to a nonzero rational multiple).
PRIMITIVE = {
    "C09": 8 * X**2 * D**2 - 21 * X**2 * S + 6 * D * S**2,
    "C02": 14 * X**2 * D - 3 * S**2,
    "C22": (13797 * X**4 + 1952 * X**2 * D**3 + 1476 * X**2 * D * S
            + 324 * D**2 * S**2 + 756 * S**3),
}
def verify_top3(terms):
    results = {}
    for name, level, support, expect_freq in TOP3:
        F = tie_polynomial(level, support, terms)
        IF = initial_form(F)

        # (i) IF matches the RESIDUE_LEMMAS primitive relation up to nonzero const
        ratio = sp.simplify(IF / PRIMITIVE[name])
        _require(ratio.free_symbols == set() and ratio != 0, name)

        # (ii) jet decomposition identities
        C0, C1, C2 = jet_coefficients(F)
        _require(sp.expand(C0 - IF) == 0, f"{name}: C0 != IF")
        _require(sp.expand(C1 - grad_dot(IF, J1)) == 0, f"{name}: C1")
        _require(sp.expand(C2 - grad_dot(IF, J2) - hess_quadratic(IF, J1)) == 0, f"{name}: C2")

        # (iii) smooth rational nonzero-leading point (verified against source)
        cert = smooth_certificate(support, IF)
        _require(cert is not None, f"{name}: no smooth certificate")
        pt, grad = cert
        present = sorted(IF.free_symbols, key=sp.default_sort_key)

        # (iv) solve the depth-2 and depth-3 jet equations at that point with a
        #      NONZERO forced lower jet, exhibiting a nontrivial depth-3 point.
        piv = next(v for v in present if grad[v] != 0)
        piv1 = J1[LEAD.index(piv)]
        piv2 = J2[LEAD.index(piv)]
        base = {v: pt.get(v, 0) for v in LEAD}
        # force present non-pivot lower jets to distinct nonzero values (a
        # nonzero inhomogeneity that avoids accidental cancellation).
        vals = [2, 3, 5, 7]
        nonpiv = [v for v in present if v is not piv]
        forced1 = {J1[LEAD.index(v)]: vals[i] for i, v in enumerate(nonpiv)}
        zero1 = {j: 0 for j in J1 if j not in forced1 and j is not piv1}
        c1 = C1.subs(base).subs(forced1).subs(zero1)
        j1val = sp.solve(sp.Eq(c1, 0), piv1)[0]
        j1sub = {**forced1, **zero1, piv1: j1val}
        _require(sp.expand(C1.subs(base).subs(j1sub)) == 0, "sp.expand(C1.subs(base).subs(j1sub)) == 0")
        # order 2: force present non-pivot order-2 jets, solve C2=0
        forced2 = {J2[LEAD.index(v)]: vals[i] for i, v in enumerate(nonpiv)}
        zero2 = {j: 0 for j in J2 if j not in forced2 and j is not piv2}
        c2 = C2.subs(base).subs(j1sub).subs(forced2).subs(zero2)
        j2val = sp.solve(sp.Eq(c2, 0), piv2)[0]
        j2sub = {**forced2, **zero2, piv2: j2val}
        full = {**base, **j1sub, **j2sub}
        _require(sp.expand(C0.subs(full)) == 0, "sp.expand(C0.subs(full)) == 0")
        _require(sp.expand(C1.subs(full)) == 0, "sp.expand(C1.subs(full)) == 0")
        _require(sp.expand(C2.subs(full)) == 0, "sp.expand(C2.subs(full)) == 0")

        results[name] = {
            "IF": IF, "pt": pt, "pivot": piv,
            "j1": (piv1, j1val), "j2": (piv2, j2val),
        }
    return results


def main():
    print("Residue-lemma DEPTH verifier (a-quantified jet structure)\n")
    h = source_h()
    terms = source_terms(h)
    _require(tuple(len(terms[l]) for l in (6, 5, 4)) == (3, 5, 8), "tuple(len(terms[l]) for l in (6, 5, 4)) == (3, 5, 8)")
    print("V1. parsed f31_graded.txt; h6/h5/h4 term counts 3/5/8            OK")

    data, surv = survivors()
    ncases = sum(len(b["survivor_cases"]) for b in surv)
    avals = sorted({b["a_t"] for b in surv})
    _require(len(surv) == 171 and ncases == 1145, "len(surv) == 171 and ncases == 1145")
    _require(avals == [2, 3, 4, 5, 6, 7, 8, 9, 10], "avals == [2, 3, 4, 5, 6, 7, 8, 9, 10]")
    freq, cells, av, dep, mism = census(surv, terms)
    _require(mism == 0, "mism == 0")
    print(f"V2. frontier: {len(surv)} branches / {ncases} cases, a in "
          f"{avals[0]}..{avals[-1]}; 0 tied-vs-source mismatches         OK")

    # affine law: dominant t-place term_cancellation depth == 30-3a
    tdep = defaultdict(Counter)
    for b in surv:
        for case in b["survivor_cases"]:
            for w in case["witness"]:
                if w["place"] == "t":
                    for ob in w["obligations"]:
                        if ob["kind"] == "term_cancellation":
                            tdep[b["a_t"]][ob["depth"]] += 1
    for a in range(2, 10):
        _require(tdep[a].most_common(1)[0][0] == 30 - 3 * a, a)
    print("V3. affine t-depth law depth = 30-3a holds (a=2..9)             OK")

    # top-3 identification and frequencies
    tkeys = {k: freq[k] for k in freq
             if k[0] == "t" and k[2] == "monomial_tie_rise" and k[3]}
    top = sorted(tkeys, key=lambda k: -tkeys[k])[:3]
    got = {(k[1], k[3]): tkeys[k] for k in top}
    for name, level, support, ef in TOP3:
        _require(got.get((level, support)) == ef, (name, got.get((level, support))))
    print("V4. top-3 a-growing t window patterns = C09/C02/C22 "
          "(201/196/160)   OK")

    res = verify_top3(terms)
    for name in ("C09", "C02", "C22"):
        r = res[name]
        print(f"    {name}: IF={sp.srepr(r['IF'])[:0] or r['IF']}")
        print(f"        depth-2 eqn  C1 = grad(IF).jet1 = 0")
        print(f"        depth-3 eqn  C2 = grad(IF).jet2 + (1/2)jet1^T H jet1 = 0")
        print(f"        smooth point {r['pt']}; pivot {r['pivot']}")
        print(f"        exact jets   {r['j1'][0]}={r['j1'][1]}, "
              f"{r['j2'][0]}={r['j2'][1]}")
    print("V5. jet identities C0=IF, C1, C2 exact; depth-3 solution exhibited OK")

    # general no-kill theorem: every t-place monomial_tie_rise support smooth
    tsupports = {}
    for k in freq:
        if k[0] == "t" and k[2] == "monomial_tie_rise" and k[3]:
            tsupports.setdefault((k[1], k[3]), 0)
            tsupports[(k[1], k[3])] += freq[k]
    for (level, support), _ in sorted(tsupports.items()):
        F = tie_polynomial(level, support, terms)
        _require(tie_weight_exists(support), (level, support))
        IF = initial_form(F)
        cert = smooth_certificate(support, IF)
        _require(cert is not None, ("no smooth certificate", level, support))
    print(f"V6. all {len(tsupports)} t-place tie supports smooth at a rational "
          f"point   OK")
    print("    => every t-place monomial_tie_rise is a CONSTRAINT at all "
          "depths;")
    print("       no jet obstruction (KILL) exists at the t-place, all a.")

    # kill accounting: C08/C20 supports absent from the survivors
    C08 = frozenset({(2, 2, 0, 0), (1, 0, 0, 1), (0, 0, 0, 2)})
    C20 = frozenset({(3, 2, 0, 0), (2, 1, 0, 1), (1, 0, 0, 2)})
    kill_hits = [k for k in freq if k[3] in (C08, C20)]
    _require(kill_hits == [], "kill_hits == []")
    print("V7. C08/C20 kill supports occur 0x on the 1145 survivors "
          "(a-indep, q, d1) OK")

    print("\nALL DEPTH RESIDUE-LEMMA CHECKS PASS")
    print("Verdict: top-3 (C09/C02/C22) and all t-place ties are CONSTRAINT "
          "for every a in [2,10]; depth extension yields no new KILL.")


if __name__ == "__main__":
    main()
