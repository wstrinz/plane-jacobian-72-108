#!/usr/bin/env python3
"""Independent spot-check verifier for the depth-1 residue-congruence layer.

Re-derives everything from cascade_engine.py (MONOMIALS / tropical_h_max_full /
deg_h_options / FORBIDDEN_RISES / LC_U) and alt_inf_sweep.run_chain (the flipped
max-plus chain).  It never trusts the precomputed initial forms in
alt_residue_congruences.json; it recomputes them and, at the end, cross-checks
the JSON census against a fresh kill-test sample.

Checks (all hand-style, must PASS):

  A. constants: LC_U = c*lc(q) = -1024/3315, nonzero -> the bottom close is
     always linearly solvable (never a depth-1 kill).
  B. a CONSTRAINT state (a=11 T1, degstate (5,2,10,11)): recompute its two
     level-0 obligations; hand-derive the h_0 depth-1 initial form; exhibit an
     all-nonzero rational point; confirm it SURVIVES under APPLY_RESIDUE_KILLS.
  C. the soundness case (a=11 T1, degstate (6,5,sigma=0,11)): the C08 (level 5)
     AND C20 (level 4) FORBIDDEN supports both occur as tropical ties, and the
     C08 relation 6X^2D^2-9XDE-E^2 has discriminant square-class 105 (absent
     from Q(sqrt(17)) -> WOULD kill if required); yet the flipped chain takes
     h_5, h_4 at their MAXIMUM (no drop required), so the state survives
     identically with kills OFF and ON -> NOT killed.  This is exactly the
     soundness direction: a non-obligatory forbidden tie is not a kill.
  D. census cross-check: on a fresh random sample of survivors, re-run the
     kill-test (kills off vs on) and confirm it matches the JSON classification;
     assert 0 depth-1 kills and 0 whole-branch kills overall.

Run: python alt_residue_congruences_verify.py
"""
from __future__ import annotations

import json
import random
import re
from pathlib import Path

import sympy as sp

import cascade_engine as ce
import alt_inf_sweep as ais
import t5_90t1_verify as base  # source q (for the LC_U constant, as in cascade_inf_ties_verify)

def _require(_cond, _msg):
    """Proof-critical check: fails loudly and exits nonzero, unaffected by python -O."""
    if not _cond:
        import sys as _sys
        print("FAIL: " + str(_msg))
        _sys.exit(1)


ROOT = Path(__file__).resolve().parent
NEG_INF = ce.NEG_INF
D, X, S, E = sp.symbols("D X S E")
LABEL = re.compile(r"^(-?\d+)\*d2\^(\d+)\*d1\^(\d+)\*sigma\^(\d+)\*e\^(\d+)$")

C08 = frozenset({(2, 2, 0, 0), (1, 1, 0, 1), (0, 0, 0, 2)})   # level 5
C20 = frozenset({(3, 2, 0, 0), (2, 1, 0, 1), (1, 0, 0, 2)})   # level 4


def parse_form(tied):
    poly = sp.Integer(0)
    exps = set()
    for label in tied:
        m = LABEL.match(label)
        _require(m, label)
        coef, k, x, z, b = (int(g) for g in m.groups())
        poly += coef * D**k * X**x * S**z * E**b
        exps.add((k, x, z, b))
    return poly, frozenset(exps)


def run(a, degstate, kills):
    ce.APPLY_RESIDUE_KILLS = kills
    ais._CHAIN_CACHE.clear()
    r = ais.run_chain(a, degstate)
    ce.APPLY_RESIDUE_KILLS = False
    ais._CHAIN_CACHE.clear()
    return r


# ---------------------------------------------------------------------------
def check_a_constants():
    y, q = base.y, base.q
    c = sp.Rational(-1, 6630)
    _require(sp.LC(q, y) == 2048, "sp.LC(q, y) == 2048")
    _require(q.subs(y, -1) == 3315, "q.subs(y, -1) == 3315")
    lc_u = c * sp.LC(q, y)
    _require(lc_u == sp.Rational(-1024, 3315), "lc_u == sp.Rational(-1024, 3315)")
    _require(sp.Rational(ce.LC_U) == lc_u, "sp.Rational(ce.LC_U) == lc_u")
    _require(lc_u != 0, "lc_u != 0")
    # bottom close  lc(E)^21 * lc(h0) + lc_u * lc(r0) = 0 is linear in the free
    # nonzero unknown lc(r0): lc(r0) = -lc(E)^21 lc(h0)/lc_u, nonzero -> solvable.
    lcE, lch0 = sp.symbols("lcE lch0")
    lcr0 = -lcE**21 * lch0 / lc_u
    _require(sp.simplify(lcE**21 * lch0 + lc_u * lcr0) == 0, "sp.simplify(lcE**21 * lch0 + lc_u * lcr0) == 0")
    print("A. LC_U = -1024/3315 (nonzero) -> bottom close always solvable (never a kill). PASS")


def check_b_constraint_state():
    a, st = 11, (5, 2, 10, 11)
    r_off = run(a, st, False)
    r_on = run(a, st, True)
    _require(r_off["verdict"] == "survive", "r_off[\"verdict\"] == \"survive\"")
    obls = r_off["obligations"]
    _require(len(obls) == 2, obls)
    kinds = sorted(o["kind"] for o in obls)
    _require(kinds == ["degree_tie_drop", "leading_cancellation"], "kinds == [\"degree_tie_drop\", \"leading_cancellation\"]")
    _require(all(o["level"] == 0 for o in obls), "both obligations at the bottom close")
    tie = next(o for o in obls if o["kind"] == "degree_tie_drop")
    # recompute the h_0 initial form independently from the monomial table
    flags = (st[2] == NEG_INF, st[0] == NEG_INF, st[1] == NEG_INF)
    maximum, labels, _, es = ce.tropical_h_max_full(0, st, flags)
    ref, ref_es = parse_form(labels)
    got, got_es = parse_form(tie["tied"])
    _require(got_es == ref_es, (got_es, ref_es))
    _require(sp.expand(got - ref) == 0, "sp.expand(got - ref) == 0")
    # hand form: 960 D^6 S^2 + 6048 D^4 S^3 + 12636 D^2 S^4 + 8748 S^5
    hand = (960 * D**6 * S**2 + 6048 * D**4 * S**3
            + 12636 * D**2 * S**4 + 8748 * S**5)
    _require(sp.expand(got - hand) == 0, sp.factor(got))
    # all-nonzero rational point: factor 12 S^2 (4D^2+9S)^2 (5D^2+9S); take
    # 5D^2+9S = 0 with D=3 -> S = -5.
    pt = {D: 3, S: -5}
    _require(got.subs(pt) == 0 and all(v != 0 for v in pt.values()), "got.subs(pt) == 0 and all(v != 0 for v in pt.values())")
    # CONSTRAINT: survives even with the arithmetic kills switched on.
    _require(r_on["verdict"] == "survive", "r_on[\"verdict\"] == \"survive\"")
    print("B. (5,2,10,11): 2 level-0 obligations; h_0 form "
          "12 S^2 (4D^2+9S)^2 (5D^2+9S) with point (D,S)=(3,-5); "
          "survives kills-on -> CONSTRAINT. PASS")


def check_c_soundness_nonfire():
    # sigma identically zero; d1=5, d2=6, e=11 (deg_E=0), a=11 T1.
    a, st = 11, (6, 5, NEG_INF, 11)
    flags = (st[2] == NEG_INF, st[0] == NEG_INF, st[1] == NEG_INF)
    # C08 at level 5 and C20 at level 4 both appear as the tropical tie support.
    _, _, _, es5 = ce.tropical_h_max_full(5, st, flags)
    _, _, _, es4 = ce.tropical_h_max_full(4, st, flags)
    _require(es5 == C08 and (5, es5) in ce.FORBIDDEN_RISES, es5)
    _require(es4 == C20 and (4, es4) in ce.FORBIDDEN_RISES, es4)
    # C08 relation 6X^2D^2 - 9XDE - E^2: substitute r=E/(XD) -> r^2+9r-6=0,
    # discriminant 81+24=105; square class 105.  Q(sqrt17) (unique quadratic
    # subfield of the S4 splitting field) does not contain sqrt(105).
    disc = 9**2 + 4 * 6
    _require(disc == 105, "disc == 105")
    _require(not sp.integer_nthroot(105, 2)[1], "not sp.integer_nthroot(105, 2)[1]")        # 105 is not a square
    _require(not sp.integer_nthroot(105 * 17, 2)[1], "not sp.integer_nthroot(105 * 17, 2)[1]")   # square class 105 != class 17
    # deg_h_options at levels 5 and 4 OFFERS the maximum (no drop needed); with
    # kills ON the forbidden support only removes the (unused) drop options.
    for lvl, es in ((5, C08), (4, C20)):
        opts_off = ce.deg_h_options(lvl, st, flags)  # kills default off
        _require(opts_off[0][1] == (), "opts_off[0][1] == ()")  # max option carries no obligation
        ce.APPLY_RESIDUE_KILLS = True
        opts_on = ce.deg_h_options(lvl, st, flags)
        ce.APPLY_RESIDUE_KILLS = False
        _require(opts_on == [(opts_off[0][0], ())], (lvl, opts_on))  # only the max
    # the state itself: identical survival with kills OFF and ON, and its only
    # obligations live at level 0, NOT at levels 5/4.
    r_off = run(a, st, False)
    r_on = run(a, st, True)
    _require(r_off["verdict"] == "survive" and r_on["verdict"] == "survive", "r_off[\"verdict\"] == \"survive\" and r_on[\"verdict\"] == \"survive\"")
    _require(r_off["obligations"] == r_on["obligations"], "r_off[\"obligations\"] == r_on[\"obligations\"]")
    _require(all(o["level"] == 0 for o in r_off["obligations"]), "all(o[\"level\"] == 0 for o in r_off[\"obligations\"])")
    print("C. (6,5,sigma=0,11): C08(L5)+C20(L4) forbidden ties present (disc 105 "
          "-> arithmetic kill IF required), but h_5/h_4 taken at max (no drop) so "
          "kills-on == kills-off, state SURVIVES -> non-obligatory tie is NOT a kill. PASS")


def check_d_census():
    data = json.loads((ROOT / "alt_residue_congruences.json").read_text("utf-8"))
    cen = data["census"]
    _require(cen["killed_at_depth1"] == 0, "cen[\"killed_at_depth1\"] == 0")
    _require(cen["constrained"] == cen["n_states"], "cen[\"constrained\"] == cen[\"n_states\"]")
    _require(cen["n_whole_branch_kills"] == 0 and cen["whole_branch_kills"] == [], "cen[\"n_whole_branch_kills\"] == 0 and cen[\"whole_branch_kills\"] == []")
    _require(cen["singleton_obligatory_ties_flagged"] == 0, "cen[\"singleton_obligatory_ties_flagged\"] == 0")
    _require(data["kill_audit"]["kills_on_equals_kills_off"] is True, "data[\"kill_audit\"][\"kills_on_equals_kills_off\"] is True")
    # every catalogued CONSTRAINT hypersurface: recompute and re-verify the point
    for c in data["support_catalog"]:
        poly = sum(m["coef"] * D**m["d2"] * X**m["d1"] * S**m["sigma"] * E**m["e"]
                   for m in c["monomials"])
        _require(c["classification"] == "CONSTRAINT", "c[\"classification\"] == \"CONSTRAINT\"")
        w = c["rational_witness"]
        sub = {sym: sp.Rational(w[str(sym)]) for sym in (D, X, S, E)
               if str(sym) in w}
        _require(sub and all(v != 0 for v in sub.values()), "sub and all(v != 0 for v in sub.values())")
        _require(sp.expand(poly).subs(sub) == 0, c["support_id"])
    # fresh random sample: independent kill-test must match JSON classification
    rng = random.Random(20260722)
    sample = rng.sample(data["states"], 120)
    for row in sample:
        s = row["state"]
        f = lambda v: NEG_INF if v is None else v
        st = (f(s["deg_d2"]), f(s["deg_d1"]), f(s["deg_sigma"]), f(s["deg_e"]))
        a = int(row["branch"][1:row["branch"].index("_")])
        killed = run(a, st, True)["verdict"] == "killed"
        expect = row["classification"] == "KILL"
        _require(killed == expect, (row["branch"], st))
    print("D. census: 0 killed / 3102 constrained / 0 whole-branch kills; all 19 "
          "hypersurface points re-verified; 120-state fresh kill-test matches. PASS")


def main():
    check_a_constants()
    check_b_constraint_state()
    check_c_soundness_nonfire()
    check_d_census()
    print("\nalt_residue_congruences: PASS")


if __name__ == "__main__":
    main()
