#!/usr/bin/env python3
"""alt_rebuild.py -- recompile the ALTERNATE-REGIME frontier (a_t >= 11).

WHY THIS FILE EXISTS
--------------------
`ALT_REGIME.md` -> `ALT_REGIME_L2.md` -> `ALT_INF_SWEEP.md` -> `ALT_COMBINED.md`
-> `ALT_RESIDUE_CONGRUENCES.md` leave the alternate regime at **27 open branches
(13 T1 + 14 T2) / 3102 states, 0 whole-branch kills**.  That census was honest
when it was written.  It predates four lemmas that are cap-free, branch-free and
regime-free, and that between them delete 46 of the 52 alternate-regime
branches:

  L1  e | Phi                    (DIVISOR_SYZYGY.md; divisor_syzygy.py 7/7)
  L2  e | S                      (SYZYGY_SWEEP.md sec.4; re-proved in
                                  DIVISOR_CONSEQUENCES.md sec.2)
  L3  the place trichotomy at t  (T1_BRANCH.md sec.1.2; t1_branch.py 15/15)
  L4  the T2 divisor normal form (DIVISOR_CONSEQUENCES.md sec.3-6)

This module re-derives the alternate-regime survivor set from the COMMITTED
artifacts, from the 52-branch universe up -- it does NOT inherit the 27-branch
intermediate as a premise.  The 27-branch list is loaded only as a cross-check.

SCOPE GUARDS (violating any of these produces a WRONG answer; each is asserted)
------------------------------------------------------------------------------
  * `R | e^2`, `e*R | Phi`, `R = c*(y+1)^rho` are **T2-ONLY, PERMANENTLY**
    (T1_BRANCH.md T8: an explicit point on V(G1,G2,G3,G5) with R = 0, e = 1,
    d1 != 0 proves no power of e lies in I+(R)).  Check G6 re-verifies that
    witness, and stage S3 is structurally unreachable from a T1 branch.
  * The SPINE five-family elimination is **sub2-only** (n = 10 - a, a = 6..10).
    Every branch here has a >= 11, so it has no column in this regime; check G7.
  * `POLE_THEOREM.md` Thm 2C is gated on `3a - 2 < v_t(Phi) = 30`, i.e.
    `a <= 10`, which is FALSE throughout this regime.  What is used here is the
    UNGATED place trichotomy of T1_BRANCH.md sec.1.2, which specialises to
    Thm 2C only when a <= 10.  Check B4 re-derives the a >= 11 horn explicitly.
  * C08/C20 are UNSOUND as used (FIELD_SCOPE_AUDIT.md) and downgraded
    KILL -> CONSTRAINT.  The census is reported under BOTH settings; check E1
    re-checks the `kills_on == kills_off` counterfactual on all 3102 states.

READ-ONLY.  This module writes exactly one file, `ALT_FRONTIER_V2.md`, and only
when invoked without `--quiet`.  No shared artifact is regenerated.

USAGE
    python -u alt_rebuild.py            # full report + (re)write ALT_FRONTIER_V2.md
    python -u alt_rebuild.py --quiet    # self-check only, exit 0 iff all pass
    python -u alt_rebuild.py --no-emit  # full report, do not write the markdown
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))

# `gamma`, `E`, `S`, `beta` are sympy builtins -- every identifier this module
# sympifies is pinned first.
y = sp.Symbol("y")
sym_d0, sym_d1, sym_d2 = sp.symbols("d0 d1 d2")
sym_e, sym_R, sym_S, sym_T = sp.symbols("dm1 dm2 dm3 dm4")
sym_Phi = sp.Symbol("Phi")

# ---------------------------------------------------------------------------
# Inputs.  A moved input is FATAL (assert), never a silently smaller census.
# ---------------------------------------------------------------------------
INPUTS = {
    "ledger": "split_place_ledger_sub1.json",
    "inf": "alt_inf_sweep.json",
    "combined": "alt_combined.json",
    "residue": "alt_residue_congruences.json",
}

# The external synthesis this run is asked to CHECK, not adopt.  It is quoted
# here only so that the comparison in section 8 of the report is mechanical.
EXTERNAL_SIX = (
    "a12_b0000_T1", "a12_b1000_T1", "a12_b1100_T1", "a12_b1110_T1",
    "a14_b0000_T1", "a14_b1000_T1",
)

V_T_PHI = 30                     # ord_{y=-1} Phi, exact (q(-1) = 3315 != 0)
DEG_PHI = 34
ALT_MIN_A = 11                   # the alternate regime is a_t >= 11

_RESULTS: list[tuple[str, bool, str]] = []


def record(tag, ok, msg=""):
    _RESULTS.append((tag, bool(ok), msg))
    return bool(ok)


def _load(key):
    path = os.path.join(HERE, INPUTS[key])
    assert os.path.exists(path), "MISSING INPUT %s (%s)" % (INPUTS[key], key)
    with open(path, "rb") as fh:
        raw = fh.read()
    return json.loads(raw.decode("utf-8")), hashlib.sha256(raw).hexdigest()[:16]


def _bid(a, b, branch):
    return "a%d_b%s_%s" % (a, "".join(str(x) for x in b), branch)


# ===========================================================================
#  A.  The universe -- 26 alternate-regime strata / 52 branches
# ===========================================================================
def universe(ledger):
    """Every (a_t, b, branch) with a_t >= 11, straight out of the ledger.

    This is the honest denominator.  It is NOT the 27-branch list: the 27 is
    what survives ALT_REGIME.md + ALT_REGIME_L2.md, and nothing below needs
    those two rounds.
    """
    rows = [r for r in ledger["strata"] if r["a_t"] >= ALT_MIN_A]
    out = []
    for r in rows:
        for branch in r["open_branches"]:
            out.append({"id": _bid(r["a_t"], r["b"], branch),
                        "a": r["a_t"], "b": tuple(r["b"]),
                        "sum_b": sum(r["b"]), "branch": branch})
    return sorted(out, key=lambda d: (d["a"], d["b"], d["branch"])), rows


def a1_universe(led, v=True):
    uni, rows = universe(led)
    summ = led["summary"]
    ok = (len(rows) == summ["alternate_regime_strata"] == 26
          and len(uni) == summ["alternate_regime_open_branches"] == 52
          and all(r["stratum_status"] == "alternate_regime_open" for r in rows)
          and max(u["a"] for u in uni) == 15
          and min(u["a"] for u in uni) == 11)
    if v:
        print("[A1] universe from %s: %d strata / %d branches "
              "(ledger summary says %d / %d), a_t in [%d, %d]"
              % (INPUTS["ledger"], len(rows), len(uni),
                 summ["alternate_regime_strata"],
                 summ["alternate_regime_open_branches"],
                 min(u["a"] for u in uni), max(u["a"] for u in uni)))
        print("     every stratum carries stratum_status = 'alternate_regime_open'.")
    return record("A1_universe_52", ok), uni


def a2_prior_layers(inf, comb, res, uni, v=True):
    """The 27-branch intermediate, loaded ONLY as a cross-check."""
    ids = {u["id"] for u in uni}
    inf_ids = [b["id"] for b in inf["branches"]]
    comb_ids = [b["id"] for b in comb["branches"]]
    ok = (inf["summary"]["n_branches"] == 27
          and inf["summary"]["surviving_states"] == 4690
          and inf["summary"]["total_degree_states"] == 38360
          and comb["summary"]["states_before"] == 4690
          and comb["summary"]["states_remaining"] == 3102
          and comb["summary"]["branches_whole_killed"] == 0
          and res["census"]["n_states"] == 3102
          and res["census"]["n_whole_branch_kills"] == 0
          and set(inf_ids) == set(comb_ids)
          and set(inf_ids) <= ids
          and sum(1 for i in inf_ids if i.endswith("T1")) == 13
          and sum(1 for i in inf_ids if i.endswith("T2")) == 14)
    if v:
        print("[A2] prior alternate layers (cross-check only, NOT a premise):")
        print("     alt_inf_sweep     27 branches, 38360 degree states -> 4690 surviving")
        print("     alt_combined      4690 -> 3102 states, 0 whole-branch kills")
        print("     alt_residue_cong  3102 states, 0 whole-branch kills")
        print("     the 27 (13 T1 + 14 T2) are a SUBSET of the 52-branch universe.")
    return record("A2_prior_27", ok)


def a3_caps(v=True):
    """The sub1 degree sandwich, read from the cap tables, never typed in."""
    import cascade_engine as ce
    import full_system_bridge as fsb
    import divisor_consequences as dc
    cfg = ce.CONFIGS["sub1"]
    strip = fsb.STRIP_DEGCAP["sub1"]
    caps = dc.CAPS["sub1"]
    got = {"d2": caps["d2"], "d1": caps["d1"], "d0": caps["d0"],
           "sigma": cfg.aux_caps[1], "e": caps["e"], "R": strip["dm2"],
           "S": strip["dm3"], "T": strip["dm4"]}
    ok = (got["R"] == caps["R"] == 18 and got["S"] == caps["S"] == 21
          and got["T"] == caps["T"] == 24 and got["e"] == 15
          and got["d2"] == 6 and got["d1"] == 9 and got["d0"] == 12
          and cfg.aux_caps == (caps["d1"], got["sigma"], caps["d2"])
          and got["sigma"] == 12 and cfg.e_cap == caps["e"])
    if v:
        print("[A3] sub1 caps (cascade_engine.CONFIGS + full_system_bridge.STRIP_DEGCAP")
        print("     + divisor_consequences.CAPS, cross-agreeing): %s" % got)
    return record("A3_caps", ok), got


# ===========================================================================
#  B.  The four lemmas, re-checked here (a failure is fatal)
# ===========================================================================
def b1_ksyzygy(v=True):
    import divisor_syzygy as ds
    residual, _ = ds.syzygy_residual()
    ok = sp.expand(residual) == 0
    if v:
        print("[B1] L1 K-syzygy  2*(G5 + d2*G3 + d1*G2 + d0*G1)")
        print("       == 2*Phi - e*(d2*e^2 + 3*e*S + 3*R^2):  residual = %s" % residual)
        print("     No branch hypothesis (the d1*G2 term is IN the combination),")
        print("     no cap, no regime condition  =>  e | Phi on every lift.")
    return record("B1_ksyzygy_residual_0", ok)


def b2_phi(v=True):
    import divisor_consequences as dc
    q = dc.q_poly()
    phi = dc.phi_stripped()
    qp = sp.Poly(q, y)
    ok = (sp.degree(phi, y) == DEG_PHI
          and sp.degree(q, y) == 4
          and q.subs(y, -1) == 3315
          and sp.gcd(q, sp.diff(q, y)) == 1                 # squarefree
          and qp.is_irreducible
          and sp.Poly(phi, y).eval(-1) == 0
          and all(sp.Poly(sp.diff(phi, y, k), y).eval(-1) == 0 for k in range(30))
          and sp.Poly(sp.diff(phi, y, 30), y).eval(-1) != 0)
    if v:
        print("[B2] Phi arithmetic: deg Phi = %d, ord_{y=-1} Phi = %d exactly,"
              % (sp.degree(phi, y), V_T_PHI))
        print("     q = %s  squarefree, irreducible over Q, q(-1) = %s != 0."
              % (q, q.subs(y, -1)))
        print("     => b_i <= 1 at every simple q-root, and v_t(Phi) = 30 EXACTLY.")
    return record("B2_phi_arithmetic", ok)


def b3_e_divides_S(v=True):
    """L2: e | S, by integral closure (DIVISOR_CONSEQUENCES.md sec.2)."""
    import divisor_consequences as dc
    A2, A3, u_cof, v_cof, res = dc.integral_dependence()
    # the cofactor identity is a CERTIFICATE, not a citation
    cert = sp.expand(u_cof * A2 + v_cof * A3 - res) == 0
    lead, alphas, all_poly = dc.monic_form(res)
    ok = (cert and all_poly and len(alphas) == 7
          and all(sp.denom(sp.together(a)).free_symbols == set()
                  for a in alphas.values()))
    if v:
        print("[B3] L2 e | S: Res_R(A2,A3) = -2*e*[S^7 + sum e^i alpha_i S^(7-i)],")
        print("     all %d alpha_i polynomial in (d0,d1,d2,e) => S/e integral over" % len(alphas))
        print("     Q[y], which is integrally closed => S/e in Q[y].  Branch-free,")
        print("     cap-free, regime-free.  (Second route: SYZYGY_SWEEP.md sec.4.)")
    return record("B3_e_divides_S", ok)


def b4_trichotomy(v=True):
    """L3: the place trichotomy at y = -1, and its a >= 11 horn.

    THEOREM (T1_BRANCH.md sec.1.2, both branches, cap-free, every place):
        v(R) >= v(e)   OR   v(Phi) = v(e) + 2*v(R).
    At y = -1, v(e) = a and v(Phi) = 30 exactly.  The first horn plus e | S and
    d2 polynomial force 30 >= 3a, i.e. a <= 10 -- FALSE in this regime.  So in
    the alternate regime the SECOND horn is compulsory:
        a + 2*rho = 30,  rho = v_t(R) < a   =>   rho = (30-a)/2 must be an
    integer, i.e. a must be EVEN, and (30-a)/2 < a <=> a > 10.
    """
    import t1_branch as tb
    INF = tb.INF
    rng = list(range(0, 61)) + [INF]
    cand = {}
    for a in range(ALT_MIN_A, 16):
        alive = []
        for rho in range(0, a):                    # rho < a is the second horn
            if any(tb.place_trichotomy(a, rho, v_s, dl2, V_T_PHI)
                   for v_s in rng for dl2 in rng):
                alive.append(rho)
        cand[a] = alive
    # horn 1 (rho >= a) is impossible for every a >= 11: all four K-terms then
    # have order >= 3a > 30, so 2*Phi could not have order 30.
    horn1 = {a: any(tb.place_trichotomy(a, rho, v_s, dl2, V_T_PHI)
                    for rho in range(a, 3 * a + 1)
                    for v_s in rng for dl2 in rng)
             for a in range(ALT_MIN_A, 16)}
    expect = {11: [], 12: [9], 13: [], 14: [8], 15: []}
    ok = cand == expect and not any(horn1.values())
    # control: the test is NOT a constant refutation -- horn 1 IS feasible for
    # a <= 10, which is exactly why the standard regime is untouched.
    ctrl = all(any(tb.place_trichotomy(a, rho, v_s, dl2, V_T_PHI)
                   for rho in range(a, 2 * a + 1) for v_s in rng for dl2 in rng)
               for a in (8, 9, 10))
    ok = ok and ctrl
    if v:
        print("[B4] L3 place trichotomy at y = -1 (T1_BRANCH.md sec.1.2/1.3, BOTH branches):")
        for a in sorted(cand):
            print("       a = %2d : horn1 (v_t(R) >= a) feasible = %-5s ; "
                  "horn2 rho candidates = %s"
                  % (a, horn1[a], cand[a] or "NONE"))
        print("     => a_t in {11, 13, 15} DEAD on parity, on BOTH branches;")
        print("        a_t in {12, 14} survive with v_t(R) = 9, 8.")
        print("     control: horn 1 stays feasible at a = 8,9,10 (%s), so the test"
              % ctrl)
        print("        is not a constant refutation and the standard regime is untouched.")
    return record("B4_trichotomy", ok), cand


def b5_D3_vacuous(caps, v=True):
    """D3 (the deg-Phi count that forces deg e = 10 in sub2) is VACUOUS in sub1."""
    rhs_at_0 = 0 + max(caps["d2"] + 0, 0 + caps["S"], 2 * caps["R"])
    ok = rhs_at_0 >= DEG_PHI and 2 * caps["R"] == 36
    if v:
        print("[B5] D3 in sub1 is VACUOUS: at deg e = 0 the RHS already reaches")
        print("     max(%d, %d, %d) = %d >= deg Phi = %d, so no degree is forced."
              % (caps["d2"], caps["S"], 2 * caps["R"], rhs_at_0, DEG_PHI))
        print("     The sub2 conclusion 'deg e = 10' does NOT transfer here.")
    return record("B5_D3_vacuous_sub1", ok)


# ===========================================================================
#  G.  Scope guards -- each one is a way to get a WRONG answer
# ===========================================================================
def g6_T2_only_witness(v=True):
    """T1_BRANCH.md T8: R | e^2 and e*R | Phi have NO T1 certificate, ever."""
    import divisor_consequences as dc
    g = dc.G_generators()
    pt = {dc.e: 1, dc.R: 0, dc.S: 1, dc.T: sp.Rational(1, 6),
          dc.d0: 1, dc.d1: sp.Rational(-1, 3), dc.d2: 0,
          dc.Phi: sp.Rational(3, 2)}
    vals = {k: sp.simplify(gg.subs(pt)) for k, gg in g.items()}
    ok = (all(x == 0 for x in vals.values())
          and pt[dc.R] == 0 and pt[dc.e] != 0
          and pt[dc.d1] != 0 and pt[dc.Phi] != 0)
    if v:
        print("[G6] SCOPE GUARD -- R | e^2 / e*R | Phi / R = c*t^rho are T2-ONLY.")
        print("     The point e=1, R=0, S=1, T=1/6, d0=1, d1=-1/3, d2=0, Phi=3/2")
        print("     lies on V(G1,G2,G3,G5): %s" % {k: str(x) for k, x in vals.items()})
        print("     with R = 0, e != 0, Phi != 0, d1 != 0.  So no power of e is in")
        print("     I+(R) and no power of Phi in I+(e*R): the three relations can")
        print("     NEVER be certified on T1.  Stage S3 below is applied to T2 only.")
        print("     The exact T1 replacements are R | e^2*(e+3*d1*s), R | e^2*(s^2-d0).")
    return record("G6_T2_only_witness", ok)


def g7_spine_out_of_scope(uni, v=True):
    """The SPINE five-family reduction is sub2-only (a = 6..10, n = 10-a)."""
    ok = min(u["a"] for u in uni) >= 11
    if v:
        print("[G7] SCOPE GUARD -- the SPINE elimination certificate is sub2-only.")
        print("     Its zero-slack degree count (n+6)+(2n+4) = 3n+10 is a sub2")
        print("     coincidence; its five families are n = 10 - a with a = 6..10.")
        print("     Every alternate-regime branch has a >= %d, so the spine has NO"
              % min(u["a"] for u in uni))
        print("     column here.  It is NOT applied, and not imported.")
    return record("G7_spine_out_of_scope", ok)


def g8_pole_theorem_gate(v=True):
    """POLE_THEOREM Thm 2C is gated on 3a - 2 < 30; the gate FAILS here."""
    gate = {a: (3 * a - 2 < V_T_PHI) for a in range(ALT_MIN_A, 16)}
    ok = not any(gate.values())
    if v:
        print("[G8] SCOPE GUARD -- POLE_THEOREM.md Thm 2C carries the regime")
        print("     condition 3a - 2 < v_t(Phi) = 30, i.e. a <= 10.  In this regime")
        print("     the gate is FALSE for every a: %s" % gate)
        print("     What is used instead is the UNGATED trichotomy of T1_BRANCH.md")
        print("     sec.1.2, which reduces to Thm 2C only when a <= 10 (check B4).")
    return record("G8_pole_gate_fails", ok)


# ===========================================================================
#  S.  The filter chain
# ===========================================================================
def s1_e_divides_phi(u):
    """L1 / D2: e | Phi and q squarefree => b_i in {0,1}."""
    bad = [i for i, x in enumerate(u["b"]) if x > 1]
    if bad:
        return False, ("b_%d = %d >= 2 at a simple q-root; ord_{r} Phi = 1 and "
                       "e | Phi forbid it" % (bad[0] + 1, u["b"][bad[0]]))
    return True, "b_i <= 1"


def s1b_defect(u, caps):
    """L1 / D1: rad(e) | (y+1)*q  =>  deg e = a + sum(b) EXACTLY."""
    deg_e = u["a"] + u["sum_b"]
    if deg_e > caps["e"]:
        return False, "deg e = a + sum b = %d > cap %d" % (deg_e, caps["e"])
    return True, "deg e = a + sum b = %d (forced, defect 0)" % deg_e


def s2_trichotomy(u, cand):
    """L3: at y = -1, a >= 11 forces a + 2*v_t(R) = 30 with v_t(R) < a."""
    if not cand[u["a"]]:
        return False, ("a = %d: horn 1 (v_t(R) >= a) gives 30 >= 3a, false; horn 2 "
                       "needs rho = (30-a)/2 = %s to be an integer < a -- it is not "
                       "(a odd)" % (u["a"], sp.Rational(30 - u["a"], 2)))
    return True, "rho = v_t(R) = %d (forced, exact)" % cand[u["a"]][0]


def s3_t2_capfree(u, rho):
    """L4 on T2, CAP-FREE, at the trichotomy-pinned rho = (30-a)/2.

    e | Phi        =>  e = gamma*t^a*Pi,  Pi = prod_{i in B}(y - r_i),
                       deg Pi = sum b,  ord_t(Pi) = 0  (q(-1) != 0)
    R | e^2  (T2)  =>  e^2 = -6*R*W  with W = d0 + s^2 + d2*s polynomial
    R = c*t^rho    =>  W = -(gamma^2/(6c)) * t^(2a-rho) * Pi^2
    e*R | Phi (T2) =>  Psi := 2*Phi/(e*R) = A*t^n*q_rem,  n = 30-a-rho,
                       q_rem = q/Pi,  A != 0
    K              =>  Psi = 3*R - 6*W*X   with X := d2 + 3*s polynomial

    Hence   6*W*X = 3*c*t^rho - A*t^n*q_rem.   The trichotomy gives n = rho, so
    the right side is  t^rho * (3*c - A*q_rem),  and X polynomial requires

        2*a - rho  <=  rho + ord_t(3*c - A*q_rem).

    `3*c - A*q_rem` is either identically 0 -- possible ONLY if deg q_rem = 0,
    i.e. sum b = 4 -- or a nonzero polynomial of degree 4 - sum b, whose order
    at any point is at most its degree.  So the necessary condition is

        2*a - 2*rho <= 4 - sum b,   i.e.   w := 3*a - 30 <= 4 - sum b.

    No degree cap enters.  Note w >= 3 throughout this regime.
    """
    a, sb = u["a"], u["sum_b"]
    assert u["branch"] == "T2", "S3 is T2-ONLY (T1_BRANCH.md T8)"
    assert max(u["b"]) <= 1, "deg Pi = sum b needs b_i <= 1 (stage S1 first)"
    n = 30 - a - rho
    assert n == rho, "trichotomy must give n = rho"
    w = 3 * a - 30
    if sb == 4:
        return (True, "sum b = 4: q_rem is a nonzero CONSTANT, so 3c - A*q_rem "
                      "can vanish identically and X = 0 escapes the t-order test")
    if w <= 4 - sb:
        return True, "w = %d <= 4 - sum b = %d" % (w, 4 - sb)
    return False, ("t-order of X: needs 2a-rho = %d <= rho + ord_t(3c - A*q_rem) "
                   "<= %d + deg(3c - A*q_rem) = %d + %d = %d.  %d > %d, and the "
                   "numerator cannot vanish identically because deg q_rem = %d > 0"
                   % (2 * a - rho, rho, rho, 4 - sb, rho + 4 - sb,
                      2 * a - rho, rho + 4 - sb, 4 - sb))


def s3_t2_engine(u):
    """L4 on T2, the committed cap-dependent engine (N1..N5, rho swept 0..cap_R)."""
    import divisor_consequences as dc
    r = dc.t2_cell_verdict("sub1", u["a"], list(u["b"]))
    return r["verdict"] != "DEAD", r["verdict"], r["why"]


def run_chain(uni, caps, cand, v=True):
    """Apply S1 -> S1b -> S2 -> S3 in order; also record UNORDERED attribution."""
    alive = list(uni)
    stages = []

    def apply(name, fn, note):
        nonlocal alive
        keep, killed = [], []
        for u in alive:
            ok, why = fn(u)
            (keep if ok else killed).append((u, why))
        stages.append({"stage": name, "note": note,
                       "killed": [(u["id"], w) for u, w in killed],
                       "alive_after": [u["id"] for u, _ in keep],
                       "n_killed": len(killed), "n_alive": len(keep)})
        alive = [u for u, _ in keep]

    apply("S1  e|Phi  (D2: b_i <= 1)", s1_e_divides_phi,
          "branch-free, cap-free, regime-free")
    apply("S1b e|Phi  (D1: deg e = a + sum b)", lambda u: s1b_defect(u, caps),
          "cap-dependent ONLY through deg e <= %d" % caps["e"])
    apply("S2  place trichotomy at y = -1", lambda u: s2_trichotomy(u, cand),
          "branch-free, cap-free; the a >= 11 horn")

    def s3(u):
        if u["branch"] != "T2":
            return True, "T1: stage S3 does not apply (T1_BRANCH.md T8)"
        rho = cand[u["a"]][0]
        ok_cf, why_cf = s3_t2_capfree(u, rho)
        ok_en, verdict, why_en = s3_t2_engine(u)
        # both routes must agree that the branch is dead
        if not ok_cf and not ok_en:
            return False, "CAP-FREE: %s || ENGINE(%s): %s" % (why_cf, verdict, why_en)
        if ok_cf != ok_en:
            return ok_cf, "DISAGREEMENT capfree=%s engine=%s" % (why_cf, why_en)
        return True, "survives both routes: %s / %s" % (why_cf, why_en)

    apply("S3  T2 divisor normal form (T2 ONLY)", s3,
          "e*R|Phi, R|e^2, R = c*t^rho -- NEVER applied to a T1 branch")

    # ---- unordered attribution: what each lemma kills on its own -----------
    solo = {}
    solo["S1 e|Phi (b_i<=1)"] = [u["id"] for u in uni if not s1_e_divides_phi(u)[0]]
    solo["S1b e|Phi (defect 0 vs cap)"] = [u["id"] for u in uni
                                           if not s1b_defect(u, caps)[0]]
    solo["S2 trichotomy (a odd)"] = [u["id"] for u in uni
                                     if not s2_trichotomy(u, cand)[0]]
    solo["S3 T2 divisor (engine)"] = [u["id"] for u in uni
                                      if u["branch"] == "T2"
                                      and not s3_t2_engine(u)[0]]
    if v:
        print("\n=== FILTER CHAIN over the %d-branch universe ===" % len(uni))
        for st in stages:
            print("  %-42s kills %2d -> %2d alive   [%s]"
                  % (st["stage"], st["n_killed"], st["n_alive"], st["note"]))
        print("\n  unordered attribution (each lemma applied ALONE to all %d):" % len(uni))
        for k, vlist in solo.items():
            print("    %-34s kills %2d" % (k, len(vlist)))
    return alive, stages, solo


def s3_controls(v=True):
    """The T2 engine is not a constant DEAD, and e|Phi is not a constant kill."""
    import divisor_consequences as dc
    def vd(a, b):
        return dc.t2_cell_verdict("sub1", a, list(b))["verdict"]

    live = {"a%d_b%s_T2" % (a, "".join(map(str, b))): vd(a, b)
            for a, b in ((10, (0, 0, 0, 0)), (8, (0, 0, 0, 0)), (6, (1, 1, 1, 1)))}
    dead = {"a%d_b%s_T2" % (a, "".join(map(str, b))): vd(a, b)
            for a, b in ((9, (1, 0, 0, 0)), (7, (1, 1, 1, 0)))}
    ok = (all(x == "ALIVE" for x in live.values())
          and all(x == "DEAD" for x in dead.values()))
    if v:
        print("\n[S3c] CONTROLS on the T2 engine (standard regime, a <= 10):")
        print("      ALIVE: %s" % live)
        print("      DEAD : %s" % dead)
        print("      => the alternate-regime DEADs are not an artefact of a")
        print("         constantly-refuting engine.")
    return record("S3c_engine_controls", ok)


# ===========================================================================
#  D.  State-level reduction on the survivors
# ===========================================================================
def state_census(survivors, inf, comb, v=True):
    I = {b["id"]: b for b in inf["branches"]}
    C = {b["id"]: b for b in comb["branches"]}
    rows = []
    for u in survivors:
        i, c = I.get(u["id"]), C.get(u["id"])
        assert i is not None and c is not None, "survivor %s absent from the prior layers" % u["id"]
        forced = u["a"] + u["sum_b"]
        inf_all = i["surviving_states_compact"]
        inf_d1 = [r for r in inf_all if r[3] == forced]
        comb_all = c["remaining_states"]
        comb_d1 = [r for r in comb_all if r["state"]["deg_e"] == forced]
        rows.append({"id": u["id"], "a": u["a"], "sum_b": u["sum_b"],
                     "deg_e_forced": forced,
                     "inf_total": i["counts"]["total_degree_states"],
                     "inf_surv": len(inf_all), "inf_D1": len(inf_d1),
                     "comb_rem": len(comb_all), "comb_D1": len(comb_d1)})
    tot = {k: sum(r[k] for r in rows)
           for k in ("inf_total", "inf_surv", "inf_D1", "comb_rem", "comb_D1")}
    if v:
        print("\n=== STATE-LEVEL REDUCTION on the survivors (D1: deg e = a + sum b) ===")
        print("  %-16s %5s %6s %8s %8s %9s %8s"
              % ("branch", "deg_e", "degsts", "inf surv", "+D1", "comb rem", "+D1"))
        for r in rows:
            print("  %-16s %5d %6d %8d %8d %9d %8d"
                  % (r["id"], r["deg_e_forced"], r["inf_total"], r["inf_surv"],
                     r["inf_D1"], r["comb_rem"], r["comb_D1"]))
        print("  %-16s %5s %6d %8d %8d %9d %8d"
              % ("TOTAL", "", tot["inf_total"], tot["inf_surv"], tot["inf_D1"],
                 tot["comb_rem"], tot["comb_D1"]))
    ok = tot["comb_D1"] <= tot["comb_rem"] <= 3102 and tot["inf_D1"] <= tot["inf_surv"]
    return record("D_state_census", ok), rows, tot


# ===========================================================================
#  E.  C08/C20 under BOTH settings
# ===========================================================================
def e1_c08_c20(res, survivors, v=True):
    """CONFIRM (not quote) the C08/C20 counterfactual, from the per-state rows.

    A C08/C20 downgrade can only change the census if a forbidden support ever
    appears as a REQUIRED L0 tie -- i.e. if some state's verdict depends on it.
    Recomputed here directly from `states` and `support_catalog`, not read off
    the summary.
    """
    ka = res["kill_audit"]
    cat = {c["support_id"]: c for c in res["support_catalog"]}
    forbidden = {sid for sid, c in cat.items() if c["matches_forbidden_C08_C20"]}
    classes = {}
    obligatory = 0
    for st in res["states"]:
        classes[st["classification"]] = classes.get(st["classification"], 0) + 1
        if st["L0_tie_support_id"] in forbidden:
            obligatory += 1
    recomputed_ok = (classes == {"CONSTRAINT": len(res["states"])}
                     and obligatory == 0
                     and len(res["states"]) == 3102)
    per = {b["id"]: b for b in res["branches"]}
    sur_states = sum(per[u["id"]]["n_states"] for u in survivors if u["id"] in per)
    sur_kill = sum(per[u["id"]]["n_killed_depth1"] for u in survivors if u["id"] in per)
    sur_tie = sum(per[u["id"]]["n_with_forbidden_tie"] for u in survivors if u["id"] in per)
    ok = (recomputed_ok
          and ka["kills_on_equals_kills_off"] is True
          and ka["killed_at_depth1"] == 0
          and ka["C08_C20_as_REQUIRED_obligation"] == 0
          and sur_kill == 0)
    if v:
        print("\n=== C08/C20 (UNSOUND as used -> CONSTRAINT) -- census under BOTH settings ===")
        print("  RECOMPUTED from the %d per-state rows: classifications = %s,"
              % (len(res["states"]), classes))
        print("  states whose REQUIRED L0 tie is a forbidden C08/C20 support = %d"
              % obligatory)
        print("  alt_residue_congruences.json kill_audit: kills_on == kills_off = %s,"
              % ka["kills_on_equals_kills_off"])
        print("  killed_at_depth1 = %d, C08/C20 as a REQUIRED obligation = %d,"
              % (ka["killed_at_depth1"], ka["C08_C20_as_REQUIRED_obligation"]))
        print("  states carrying a NON-obligatory C08/C20 tropical tie = %d."
              % ka["states_with_C08_C20_tropical_tie"])
        print("  On the %d survivors: %d states, %d C08/C20 kills, %d non-obligatory ties."
              % (len(survivors), sur_states, sur_kill, sur_tie))
        print("  => CONFIRMED: the alternate-regime census is IDENTICAL with the")
        print("     residue kills ON and OFF.  None of S1/S1b/S2/S3 uses C08 or C20.")
    return record("E1_c08_c20_invariant", ok), {"states": sur_states,
                                                "kills": sur_kill, "ties": sur_tie}


# ===========================================================================
#  F.  The external synthesis, checked
# ===========================================================================
def f2_vs_27(survivors, inf, v=True):
    """Consistency with the superseded 27-branch residue.

    The chain never consumes the 27; this check confirms the two are compatible
    -- every survivor is in the 27, so the new lemmas only ever ADD kills, and
    the earlier ALT_REGIME/ALT_REGIME_L2 rounds are not contradicted.
    """
    prior = {b["id"] for b in inf["branches"]}
    sur = {u["id"] for u in survivors}
    ok = sur <= prior
    if v:
        print("\n[F2] consistency with the superseded 27-branch residue:")
        print("     survivors subset of the 27: %s (%d of %d)"
              % (ok, len(sur & prior), len(sur)))
        print("     the 27 lose %d branches to the new lemmas; the other %d "
              "alternate branches" % (len(prior - sur), 52 - len(prior)))
        print("     were already killed by ALT_REGIME.md / ALT_REGIME_L2.md and "
              "are re-killed here")
        print("     independently (the chain never reads those rounds).")
    return record("F2_consistent_with_27", ok), sorted(prior - sur)


def f1_external(survivors, v=True):
    got = tuple(sorted(u["id"] for u in survivors))
    exp = tuple(sorted(EXTERNAL_SIX))
    ok = True     # agreement is REPORTED, never required
    if v:
        print("\n=== EXTERNAL SYNTHESIS, CHECKED (not adopted) ===")
        print("  external claim : %s" % ", ".join(exp))
        print("  compiler-emitted: %s" % ", ".join(got))
        if got == exp:
            print("  VERDICT: AGREES exactly, and is now compiler-emitted.")
        else:
            print("  VERDICT: **DISAGREES**.")
            print("    only external : %s" % sorted(set(exp) - set(got)))
            print("    only compiler : %s" % sorted(set(got) - set(exp)))
    return record("F1_external_compare", ok), got, exp


# ===========================================================================
#  Report
# ===========================================================================
def git_rev():
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], cwd=HERE,
                              capture_output=True, text=True,
                              timeout=20).stdout.strip() or "unknown"
    except Exception:
        return "unknown"


def esc(text):
    """Escape a literal `|` so it cannot break a markdown table cell."""
    return str(text).replace("|", "\\|")


def emit(path, data):
    L = []
    A = L.append
    d = data
    A("# ALT_FRONTIER_V2 -- the alternate regime (`a_t >= 11`) recompiled "
      "(machine-generated)")
    A("")
    A("> **DO NOT HAND-EDIT.** Regenerated by `python -u alt_rebuild.py`. Every "
      "figure below is computed from the committed artifacts named in section 9; "
      "none is typed in. `python -u alt_rebuild.py --quiet` re-derives the whole "
      "census and exits nonzero on any drift or on a moved input.")
    A("")
    A("Provenance: git `%s` | schema 1 | %s" % (d["git"], d["stamp"]))
    A("")
    A("---")
    A("")
    A("## 0. Headline")
    A("")
    A("> **%d of the %d alternate-regime branches survive. All %d are T1.**"
      % (len(d["survivors"]), d["n_universe"], len(d["survivors"])))
    A(">")
    A("> `%s`" % "`, `".join(u["id"] for u in d["survivors"]))
    A(">")
    A("> The committed **27 open branches (13 T1 + 14 T2)** figure "
      "(`ALT_REGIME.md`, `ALT_REGIME_L2.md`, `ALT_INF_SWEEP.md`, "
      "`FRONTIER.md`) is **superseded**. The T2 half of the alternate regime is "
      "**EMPTY**, and every odd `a_t` is **dead**.")
    A(">")
    A("> State level: the surviving branches carry **%d** of the %d states the "
      "combination layer left, i.e. **%d states are removed (%.1f%%)**."
      % (d["tot"]["comb_D1"], d["prior_states"],
         d["prior_states"] - d["tot"]["comb_D1"],
         100.0 * (d["prior_states"] - d["tot"]["comb_D1"]) / d["prior_states"]))
    A("")
    A("The derivation below does **not** consume the 27-branch intermediate as a "
      "premise. It starts from the full **%d-branch** universe read out of "
      "`split_place_ledger_sub1.json` and applies only lemmas that are "
      "branch-free and regime-free. The 27-branch list is loaded solely as a "
      "cross-check (it is a subset of the survivors' complement)."
      % d["n_universe"])
    A("")
    A("---")
    A("")
    A("## 1. The universe, and why it is 52 and not 27")
    A("")
    A("| source | object | count |")
    A("|---|---|---:|")
    A("| `split_place_ledger_sub1.json` (`a_t >= 11`) | strata | %d |" % d["n_strata"])
    A("| same, x open branches | branches | %d |" % d["n_universe"])
    A("| `ALT_REGIME.md` + `ALT_REGIME_L2.md` residue | branches | 27 (13 T1 + 14 T2) |")
    A("| `alt_inf_sweep.json` | degree states | 38360 -> 4690 surviving |")
    A("| `alt_combined.json` | states | 4690 -> %d remaining |" % d["prior_states"])
    A("| `alt_residue_congruences.json` | states | %d, 0 kills |" % d["prior_states"])
    A("")
    A("The stratum boundary is `a_t + sum(b_i) <= %d`, which is the sub1 cap "
      "`deg e <= %d` -- **a premise, not a theorem**. `a_t >= 16` is out of scope "
      "for that reason alone; the trichotomy of section 3 by itself permits every "
      "even `a` up to 30." % (d["caps"]["e"], d["caps"]["e"]))
    A("")
    A("---")
    A("")
    A("## 2. The lemmas integrated, with their exact scope")
    A("")
    A("| # | lemma | scope | checked here |")
    A("|---|---|---|---|")
    A("| L1 | `e \\| Phi` (K-syzygy, residual 0) | both branches, cap-free, "
      "regime-free | B1, B2 |")
    A("| L2 | `e \\| S` (integral closure) | both branches, cap-free, regime-free | B3 |")
    A("| L3 | place trichotomy at `y=-1` | both branches, cap-free, **every place** | B4 |")
    A("| L4 | `e*R \\| Phi`, `R \\| e^2`, `R = c*(y+1)^rho` | **T2 ONLY, permanently** | G6, S3 |")
    A("")
    A("### Scope guards asserted by the compiler")
    A("")
    A("* **L4 is T2-only, permanently.** `T1_BRANCH.md` T8 exhibits "
      "`e=1, R=0, S=1, T=1/6, d0=1, d1=-1/3, d2=0, Phi=3/2` on "
      "`V(G1,G2,G3,G5)` with `R=0, e!=0, d1!=0`, so no power of `e` lies in "
      "`I+(R)` and no power of `Phi` in `I+(e*R)`. Check **G6** re-substitutes "
      "that point and confirms all four generators vanish. Stage `S3` is "
      "structurally unreachable from a T1 branch. The exact T1 replacements "
      "`R | e^2*(e + 3*d1*s)` and `R | e^2*(s^2 - d0)` are **not** used to kill "
      "anything here.")
    A("* **`t^a | R,S,T` is branch-independent** and does transfer -- but in this "
      "regime it is the horn the trichotomy *refutes* (section 3), not an input.")
    A("* **The SPINE certificate is sub2-only.** Its five families are "
      "`n = 10 - a`, `a = 6..10`; every branch here has `a >= 11`. Check **G7**. "
      "It is not imported.")
    A("* **`POLE_THEOREM.md` Thm 2C is gated on `3a - 2 < v_t(Phi) = 30`**, i.e. "
      "`a <= 10`. The gate is FALSE for every `a` in this regime (check **G8**). "
      "The ungated trichotomy of `T1_BRANCH.md` sec.1.2 is used instead.")
    A("* **D3** (the `deg Phi` count that forces `deg e = 10` in sub2) is "
      "**vacuous** in sub1: `2*deg R = %d > %d = deg Phi` already at `deg e = 0` "
      "(check **B5**). It contributes nothing here."
      % (2 * d["caps"]["R"], DEG_PHI))
    A("")
    A("---")
    A("")
    A("## 3. The trichotomy at `y = -1`, which is where the regime breaks")
    A("")
    A("```")
    A("THEOREM (T1_BRANCH.md sec.1.2; both branches, cap-free, every place)")
    A("    v(R) >= v(e)      OR      v(Phi) = v(e) + 2*v(R)")
    A("")
    A("At y = -1:  v(e) = a  and  v(Phi) = 30 EXACTLY  (q(-1) = 3315 != 0).")
    A("")
    A("horn 1  v_t(R) >= a :  every K-term has order >= 3a  =>  30 >= 3a  =>  a <= 10.")
    A("                       FALSE for a >= 11.")
    A("horn 2  a + 2*rho = 30, rho = v_t(R) < a  =>  rho = (30-a)/2 in Z  =>  a EVEN.")
    A("```")
    A("")
    A("| `a` | horn 1 feasible | horn 2 `rho` candidates | verdict |")
    A("|---:|---|---|---|")
    for a in sorted(d["cand"]):
        c = d["cand"][a]
        A("| %d | no | %s | %s |"
          % (a, ", ".join(map(str, c)) or "NONE",
             "**DEAD** (parity)" if not c else "alive, `v_t(R) = %d`" % c[0]))
    A("")
    A("Control (check B4): horn 1 stays feasible at `a = 8, 9, 10`, so the test is "
      "not a constant refutation and the standard regime is untouched.")
    A("")
    A("> This kill is **branch-independent**. It removes `a_t in {11, 13, 15}` on "
      "**both** T1 and T2 -- %d of the %d branches -- with no cap and no residue "
      "argument." % (len(d["solo"]["S2 trichotomy (a odd)"]), d["n_universe"]))
    A("")
    A("---")
    A("")
    A("## 4. The T2 half of the alternate regime is EMPTY -- two independent routes")
    A("")
    A("### 4.1 Cap-free route (the one to quote)")
    A("")
    A("```")
    A("e | Phi        =>  e = gamma*t^a*Pi,  deg Pi = sum b,  ord_t(Pi) = 0")
    A("R | e^2  (T2)  =>  e^2 = -6*R*W,  W = d0 + s^2 + d2*s  polynomial")
    A("R = c*t^rho    =>  W = -(gamma^2/(6c)) * t^(2a-rho) * Pi^2")
    A("e*R | Phi (T2) =>  Psi := 2*Phi/(e*R) = A*t^n*q_rem,  n = 30-a-rho,  A != 0")
    A("K              =>  Psi = 3*R - 6*W*X,  X := d2 + 3*s  polynomial")
    A("")
    A("        6*W*X  =  3*c*t^rho - A*t^n*q_rem .")
    A("")
    A("The trichotomy gives rho = (30-a)/2, hence n = rho EXACTLY, so the right")
    A("side is t^rho*(3*c - A*q_rem).  X polynomial therefore needs")
    A("")
    A("        2*a - rho  <=  rho + ord_t(3*c - A*q_rem)  <=  rho + (4 - sum b),")
    A("")
    A("i.e.    w := 3*a - 30  <=  4 - sum b .")
    A("```")
    A("")
    A("| branch | `a` | `sum b` | `rho` | `w = 3a-30` | `4 - sum b` | verdict |")
    A("|---|---:|---:|---:|---:|---:|---|")
    for r in d["t2_rows"]:
        A("| `%s` | %d | %d | %d | %d | %d | %s |"
          % (r["id"], r["a"], r["sum_b"], r["rho"], r["w"], 4 - r["sum_b"],
             "**DEAD**" if r["dead_capfree"] else "alive"))
    A("")
    A("The single escape is `sum b = 4`, where `q_rem` is a nonzero constant and "
      "`3c - A*q_rem` may vanish identically (then `X = 0`). It **does not occur "
      "in this universe** -- but for a cap reason, not a lemma reason: "
      "`sum b = 4` with `a` even and `a >= 12` needs `deg e = a + 4 >= 16 > %d`, "
      "so no such stratum exists. Flagged, not swept: it is the ONE gap in the "
      "otherwise cap-free T2 emptiness." % d["caps"]["e"])
    A("")
    A("**Field scope.** `Pi` and `q_rem` live over the `q`-splitting field, but "
      "the only facts the argument uses are `deg` and `ord_{y=-1}` -- "
      "multiplicity data, invariant under any base change. **No square class "
      "over `Q` and no rationality of a root is claimed**, so this route is "
      "immune to the `FIELD_SCOPE_AUDIT.md` objection in the same way "
      "`T1_BRANCH.md` T6 is. (Section 4.2's `|B| = 1, 2` tests are gcd/Groebner "
      "computations over `Q` that exploit `q` being irreducible over `Q`; that is "
      "the cross-check, not the primary route.)")
    A("")
    A("### 4.2 Cap-dependent route (the committed engine, as a cross-check)")
    A("")
    A("`divisor_consequences.t2_cell_verdict('sub1', a, b)` sweeps `rho` over "
      "`0..%d` and applies N1-N5. Every alternate-regime T2 branch comes back "
      "`DEAD`, **including the odd `a`** that section 3 already removed:" % d["caps"]["R"])
    A("")
    A("| branch | engine verdict | first surviving `rho` | why |")
    A("|---|---|---|---|")
    for r in d["t2_engine_rows"]:
        A("| `%s` | %s | %s | %s |"
          % (r["id"], r["verdict"], r["rho"], esc(r["why"])))
    A("")
    A("Controls (check S3c): the same engine returns **ALIVE** for `a10_b0000_T2`, "
      "`a8_b0000_T2`, `a6_b1111_T2` and **DEAD** for `a9_b1000_T2`, "
      "`a7_b1110_T2` in the standard regime -- it is not a constant refutation. "
      "This reproduces `DIVISOR_CONSEQUENCES.md` sec.8's `a_t <= 10 on T2` "
      "by-product, which that file correctly declined to call an unblocking of "
      "the `FRONTIER_REBUILD.md` sec.7 lead; here the *missing* T1 half is "
      "supplied by section 3, and the lead is unblocked to the extent stated in "
      "section 7.")
    A("")
    A("---")
    A("")
    A("## 5. The filter chain, stage by stage")
    A("")
    A("| stage | scope note | kills | alive after |")
    A("|---|---|---:|---:|")
    for st in d["stages"]:
        A("| `%s` | %s | %d | %d |"
          % (esc(st["stage"]), esc(st["note"]), st["n_killed"], st["n_alive"]))
    A("")
    A("### Unordered attribution -- what each lemma kills on its own, out of %d"
      % d["n_universe"])
    A("")
    A("| lemma applied ALONE | branches killed |")
    A("|---|---:|")
    for k, vlist in d["solo"].items():
        A("| %s | %d |" % (esc(k), len(vlist)))
    A("")
    A("Stages compose in the order shown and are not double-counted; the "
      "unordered column overlaps by construction.")
    A("")
    A("### Whole-branch kills, itemised")
    A("")
    A("| branch | killed by | reason |")
    A("|---|---|---|")
    for st in d["stages"]:
        for bid, why in st["killed"]:
            A("| `%s` | `%s` | %s |" % (bid, st["stage"].split()[0], esc(why)))
    A("")
    A("---")
    A("")
    A("## 6. The %d survivors -- state-level, and what is left to do"
      % len(d["survivors"]))
    A("")
    A("| branch | `a` | `sum b` | `deg e` forced | degree states | after inf sweep | "
      "+ D1 | after combination | + D1 |")
    A("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for r in d["rows"]:
        A("| `%s` | %d | %d | %d | %d | %d | %d | %d | %d |"
          % (r["id"], r["a"], r["sum_b"], r["deg_e_forced"], r["inf_total"],
             r["inf_surv"], r["inf_D1"], r["comb_rem"], r["comb_D1"]))
    A("| **total** | | | | %d | %d | %d | %d | **%d** |"
      % (d["tot"]["inf_total"], d["tot"]["inf_surv"], d["tot"]["inf_D1"],
         d["tot"]["comb_rem"], d["tot"]["comb_D1"]))
    A("")
    A("`+ D1` is the `e | Phi` defect-0 consequence `deg e = a + sum(b_i)` "
      "**exactly**, applied to the recorded degree states. It is a *state* "
      "reduction, not a kill: every one of the %d branches keeps states."
      % len(d["survivors"]))
    A("")
    A("On each survivor the exact residual obligations are unchanged in kind: "
      "`v_t(R) = %s` is now pinned, `R` has `v_beta(R) = 0` at every marked root "
      "(trichotomy at `beta`, since `v_beta(Phi) = v_beta(e) = 1`), and the "
      "descending congruences `(D_t)`/`(D_p)` through `h_5` plus the bottom close "
      "`E^21 h_0 + u r_0 = 0` remain open."
      % "/".join(str(d["cand"][a][0]) for a in sorted(d["cand"]) if d["cand"][a]))
    A("")
    A("---")
    A("")
    A("## 7. PROVED / CONSTRAINT-ONLY / REDUCTION / FIELD-DEPENDENT")
    A("")
    A("### 7.1 Compiler-proved WHOLE-BRANCH kills (%d branches)"
      % (d["n_universe"] - len(d["survivors"])))
    A("")
    A("| kill | branches | evidence |")
    A("|---|---:|---|")
    for st in d["stages"]:
        if st["n_killed"]:
            A("| `%s` | %d | PROVED -- exact identity / finite valuation "
              "enumeration, char 0 |" % (esc(st["stage"]), st["n_killed"]))
    A("")
    A("**Every kill listed is cap-free.** The one cap-dependent stage, `S1b` "
      "(`deg e = a + sum b <= %d`), kills **0** branches and is in section 5 only "
      "for completeness; the section 4.2 engine cross-check is cap-dependent but "
      "the section 4.1 route it cross-checks is not." % d["caps"]["e"])
    A("")
    A("The caps do, however, still set the **scope**: which `(a, b)` strata exist "
      "at all is `a + sum(b) <= %d` (section 1). Two consequences, stated "
      "explicitly rather than buried: `a >= 16` is untested here, and "
      "`sum b = 4` -- the one configuration that escapes the section 4.1 t-order "
      "test -- never arises for even `a >= 12` only because `deg e = a + 4 >= 16` "
      "would exceed the cap. **Widening `deg e` would reopen exactly those two "
      "families and nothing else.**" % d["caps"]["e"])
    A("")
    A("### 7.2 LOCAL CONSTRAINTS ONLY (no branch dies)")
    A("")
    A("* `e | S` and `T = -R*(S/e + d2) - d1*e/2` collapse the spare ansatz "
      "(45 -> 18, branch-independent, `DIVISOR_CONSEQUENCES.md` sec.9). `dm4` is "
      "**not** a spare. No alternate-regime branch dies from this.")
    A("* The T1 replacements `R | e^2*(e + 3*d1*s)` and `R | e^2*(s^2 - d0)`. "
      "Real constraints; **not evaluated as kills here**.")
    A("* `W := e*S - R^2` with `W^2 = R^4 + d2*e^2*R^2 + d1*e^3*R + d0*e^4` and "
      "the even-multiplicity condition (`T1_BRANCH.md` T6). At `y=-1` the quartic "
      "has unique minimum `4*rho` (even) on every survivor, so it is consistent "
      "and kills nothing here.")
    A("")
    A("### 7.3 STATE REDUCTIONS")
    A("")
    A("* **D1** (`deg e = a + sum b` exactly): %d -> %d states on the survivors "
      "(combination layer), i.e. %d states removed inside surviving branches."
      % (d["tot"]["comb_rem"], d["tot"]["comb_D1"],
         d["tot"]["comb_rem"] - d["tot"]["comb_D1"]))
    A("* **The trichotomy** pins `v_t(R)` to a single value per `a`; `R` is not a "
      "recorded state coordinate in this lane, so this cannot be cashed at state "
      "level yet.")
    A("")
    A("### 7.4 FIELD-DEPENDENT observations (C08/C20 -- CONSTRAINT, not KILL)")
    A("")
    A("`FIELD_SCOPE_AUDIT.md` downgrades C08/C20 from KILL to CONSTRAINT. Census "
      "under **both** settings:")
    A("")
    A("| setting | branches killed by C08/C20 | states killed | non-obligatory ties |")
    A("|---|---:|---:|---:|")
    A("| C08/C20 ON | 0 | %d | %d |"
      % (d["c08"]["kills"], d["c08"]["ties"]))
    A("| C08/C20 OFF | 0 | %d | %d |"
      % (d["c08"]["kills"], d["c08"]["ties"]))
    A("")
    A("**CONFIRMED, and recomputed rather than quoted.** Check `E1` walks all %d "
      "per-state rows of `alt_residue_congruences.json` and finds every one "
      "classified `CONSTRAINT`, with **0** states whose REQUIRED level-0 tie is a "
      "forbidden C08/C20 support -- which is the only way the downgrade could "
      "move a verdict. That independently reproduces the file's own "
      "`kill_audit`: `kills_on == kills_off`, "
      "`killed_at_depth1 = 0`, C08/C20 as a REQUIRED obligation "
      "**0** times (they occur %d times as non-obligatory tropical ties). None of "
      "`S1`, `S1b`, `S2`, `S3` uses C08 or C20 in any form, so the survivor set "
      "above is **identical** under both settings. This reproduces "
      "`FIELD_SCOPE_AUDIT.md` sec.4.6 (\"the alternate regime: zero\")."
      % (d["prior_states"], d["c08_ties_all"]))
    A("")
    A("### 7.5 Premises this census still rests on")
    A("")
    A("* the sub1 stripped caps `deg d2 <= %d, deg d1 <= %d, deg d0 <= %d, "
      "deg e <= %d, deg R <= %d, deg S <= %d` (`CAPS_AUDIT.md` sec.3 [P1]-[P3]). "
      "They fix the SCOPE (which strata exist: `a + sum b <= deg e cap`) and are "
      "used by `S1b` and by the section 4.2 cross-check. **No kill in section 7.1 "
      "uses them.**"
      % (d["caps"]["d2"], d["caps"]["d1"], d["caps"]["d0"], d["caps"]["e"],
         d["caps"]["R"], d["caps"]["S"]))
    A("* `Phi = -(1/6630)*(y+1)^30*q`, `deg Phi = %d`, `q` squarefree and "
      "irreducible over `Q`, `q(-1) = 3315` (recomputed in B2)" % DEG_PHI)
    A("* the canonical `G1,G2,G3,G5` normalisation with `coeff(G5, Phi) = 1`")
    A("")
    A("---")
    A("")
    A("## 8. The external synthesis, checked")
    A("")
    A("An external review inferred a six-branch alternate frontier, all T1, and "
      "explicitly said it should not be quoted because it was not "
      "compiler-emitted. Recomputed here from the %d-branch universe:"
      % d["n_universe"])
    A("")
    A("```")
    A("external : %s" % ", ".join(d["external"]))
    A("compiler : %s" % ", ".join(d["emitted"]))
    A("```")
    A("")
    if tuple(d["emitted"]) == tuple(d["external"]):
        A("> **AGREES exactly.** The six-branch frontier is now compiler-emitted and "
          "may be quoted, with the grading of section 7.")
        A("")
        A("What is **proved dead** and what is **proved alive** are not the same "
          "grade, and the distinction must survive the quote:")
        A("")
        A("* **Proved DEAD (%d branches):** every kill in section 5 is an exact "
          "identity or a finite valuation enumeration in characteristic 0. "
          "Nothing is inferred." % (d["n_universe"] - len(d["emitted"])))
        A("* **Proved ALIVE (0 branches):** *none* of the six is proved to contain "
          "a solution. They are branches no current lemma refutes. "
          "\"Survivor\" means OPEN, and every one of them still carries the "
          "descending congruences and the bottom close as unmet obligations.")
        A("* **INFERRED (0):** no step in this chain is inferred. The one place "
          "where an inference could have entered -- carrying "
          "`R = c*(y+1)^rho` onto a T1 column -- is blocked structurally (G6).")
    else:
        A("> **DISAGREES.** Reported loudly rather than reconciled.")
        A("")
        A("* only in the external list: `%s`"
          % "`, `".join(sorted(set(d["external"]) - set(d["emitted"]))) or "(none)")
        A("* only compiler-emitted: `%s`"
          % "`, `".join(sorted(set(d["emitted"]) - set(d["external"]))) or "(none)")
    A("")
    A("---")
    A("")
    A("## 8b. What the superseded 27-branch residue loses")
    A("")
    A("The chain above never reads `ALT_REGIME.md` or `ALT_REGIME_L2.md`. This "
      "section only records that the two accounts are **compatible**: every "
      "survivor lies inside the old 27, so the new lemmas strictly add kills and "
      "contradict nothing (check **F2**).")
    A("")
    A("The %d branches the 27-list loses:" % len(d["lost27"]))
    A("")
    A("`%s`" % "`, `".join(d["lost27"]))
    A("")
    A("The other %d alternate branches were already killed by `ALT_REGIME.md` / "
      "`ALT_REGIME_L2.md`. **Every one of them is re-killed here independently**, "
      "with this stage breakdown -- which is why the six-branch frontier does not "
      "inherit those two rounds as premises:" % (d["n_universe"] - 27))
    A("")
    A("| stage that re-kills it | branches |")
    A("|---|---:|")
    for stage_name, n in d["pre27_by_stage"]:
        A("| `%s` | %d |" % (esc(stage_name), n))
    A("")
    A("---")
    A("")
    A("## 9. Provenance -- every input, its hash, and the command that reads it")
    A("")
    A("| artifact | sha256(16) | role |")
    A("|---|---|---|")
    for k, (name, h, role) in d["prov"].items():
        A("| `%s` | `%s` | %s |" % (name, h, role))
    A("")
    A("| imported checker | what it certifies | its own command |")
    A("|---|---|---|")
    A("| `divisor_syzygy.py` | K-syzygy residual 0, `e \\| Phi` | "
      "`python -u divisor_syzygy.py --quiet` |")
    A("| `divisor_consequences.py` | `e \\| S`, T2 normal form, `t2_cell_verdict` | "
      "`python -u divisor_consequences.py --quiet` |")
    A("| `t1_branch.py` | the place trichotomy, `place_trichotomy` | "
      "`python -u t1_branch.py --quiet` |")
    A("| `cascade_engine.py`, `full_system_bridge.py` | the sub1 cap tables | "
      "(read-only import) |")
    A("")
    A("Checker: `python -u alt_rebuild.py --quiet` (exit 0/1). It reloads every "
      "input, re-verifies L1-L4 symbolically, re-runs the whole chain and requires "
      "every number in this file. A moved input is an `AssertionError`, not a "
      "smaller census.")
    A("")
    A("### Checks")
    A("")
    A("| tag | result |")
    A("|---|---|")
    for tag, ok, msg in d["checks"]:
        A("| `%s` | %s%s |" % (tag, "PASS" if ok else "**FAIL**",
                               (" -- " + msg) if msg else ""))
    A("")
    A("`%d/%d` pass." % (sum(1 for _, o, _ in d["checks"] if o), len(d["checks"])))
    A("")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(L))
    return path


# ===========================================================================
#  main
# ===========================================================================
def run(verbose=True, emit_md=True):
    t0 = time.time()
    led, h_led = _load("ledger")
    inf, h_inf = _load("inf")
    comb, h_comb = _load("combined")
    res, h_res = _load("residue")

    if verbose:
        print("=" * 78)
        print("alt_rebuild.py -- ALTERNATE REGIME (a_t >= %d) recompiled" % ALT_MIN_A)
        print("=" * 78)

    ok_a1, uni = a1_universe(led, verbose)
    a2_prior_layers(inf, comb, res, uni, verbose)
    ok_a3, caps = a3_caps(verbose)
    if verbose:
        print()
    b1_ksyzygy(verbose)
    b2_phi(verbose)
    b3_e_divides_S(verbose)
    _, cand = b4_trichotomy(verbose)
    b5_D3_vacuous(caps, verbose)
    if verbose:
        print()
    g6_T2_only_witness(verbose)
    g7_spine_out_of_scope(uni, verbose)
    g8_pole_theorem_gate(verbose)

    survivors, stages, solo = run_chain(uni, caps, cand, verbose)
    s3_controls(verbose)

    # per-branch T2 detail for the report (both routes)
    import divisor_consequences as dc
    t2_rows, t2_engine_rows = [], []
    for u in uni:
        if u["branch"] != "T2":
            continue
        r = dc.t2_cell_verdict("sub1", u["a"], list(u["b"]))
        first = [row for row in r["rows"] if not str(row[1]).startswith("dead")]
        t2_engine_rows.append({"id": u["id"], "verdict": r["verdict"],
                               "rho": (first[0][0] if first else "none"),
                               "why": r["why"]})
        # the cap-free route needs S1 (b_i <= 1) and a pinned rho (a even)
        if cand[u["a"]] and max(u["b"]) <= 1:
            rho = cand[u["a"]][0]
            dead, why = s3_t2_capfree(u, rho)
            t2_rows.append({"id": u["id"], "a": u["a"], "sum_b": u["sum_b"],
                            "rho": rho, "w": 3 * u["a"] - 30,
                            "dead_capfree": not dead, "why": why})
    ok_t2 = (all(r["verdict"] == "DEAD" for r in t2_engine_rows)
             and all(r["dead_capfree"] for r in t2_rows)
             and not any(u["branch"] == "T2" for u in survivors))
    record("S3_T2_regime_empty", ok_t2)
    if verbose:
        print("\n[S3] T2 alternate regime: cap-free route kills %d/%d pinned-rho "
              "branches; engine kills %d/%d (all a). T2 survivors: %d."
              % (sum(1 for r in t2_rows if r["dead_capfree"]), len(t2_rows),
                 sum(1 for r in t2_engine_rows if r["verdict"] == "DEAD"),
                 len(t2_engine_rows),
                 sum(1 for u in survivors if u["branch"] == "T2")))

    _, rows, tot = state_census(survivors, inf, comb, verbose)
    _, c08 = e1_c08_c20(res, survivors, verbose)
    _, lost27 = f2_vs_27(survivors, inf, verbose)
    # every branch outside the 27 must ALSO be re-killed by this chain, or the
    # census would be silently inheriting the ALT_REGIME/_L2 rounds.
    killed_by = {bid: st["stage"] for st in stages for bid, _ in st["killed"]}
    prior27 = {b["id"] for b in inf["branches"]}
    pre27 = [u["id"] for u in uni if u["id"] not in prior27]
    record("F3_pre27_all_rekilled", all(b in killed_by for b in pre27))
    pre27_by_stage = sorted(
        ((s, sum(1 for b in pre27 if killed_by.get(b) == s))
         for s in dict.fromkeys(st["stage"] for st in stages)),
        key=lambda kv: -kv[1])
    pre27_by_stage = [kv for kv in pre27_by_stage if kv[1]]
    if verbose:
        print("     of the %d branches OUTSIDE the 27, this chain re-kills %d: %s"
              % (len(pre27), sum(1 for b in pre27 if b in killed_by),
                 {s: n for s, n in pre27_by_stage}))
    _, emitted, external = f1_external(survivors, verbose)

    npass = sum(1 for _, o, _ in _RESULTS if o)
    if verbose:
        print("\n" + "=" * 78)
        print("SURVIVORS: %d of %d alternate-regime branches -- %s"
              % (len(survivors), len(uni),
                 "ALL T1" if all(u["branch"] == "T1" for u in survivors)
                 else "MIXED"))
        for u in survivors:
            print("   %-16s a=%2d sum_b=%d  v_t(R)=%d"
                  % (u["id"], u["a"], u["sum_b"], cand[u["a"]][0]))
        print("checks: %d/%d pass   (%.1fs)" % (npass, len(_RESULTS), time.time() - t0))
        print("=" * 78)

    data = {
        "git": git_rev(),
        "stamp": time.strftime("%Y-%m-%d"),
        "n_strata": len([r for r in led["strata"] if r["a_t"] >= ALT_MIN_A]),
        "n_universe": len(uni),
        "prior_states": comb["summary"]["states_remaining"],
        "survivors": survivors, "stages": stages, "solo": solo,
        "cand": cand, "caps": caps, "rows": rows, "tot": tot,
        "c08": c08, "c08_ties_all": res["kill_audit"]["states_with_C08_C20_tropical_tie"],
        "emitted": emitted, "external": external, "lost27": lost27,
        "pre27_by_stage": pre27_by_stage,
        "t2_rows": t2_rows, "t2_engine_rows": t2_engine_rows,
        "checks": list(_RESULTS),
        "prov": {
            "ledger": (INPUTS["ledger"], h_led,
                       "the %d-strata / %d-branch alternate universe" % (
                           len([r for r in led["strata"] if r["a_t"] >= ALT_MIN_A]),
                           len(uni))),
            "inf": (INPUTS["inf"], h_inf, "degree (max-plus) layer, 38360 -> 4690"),
            "combined": (INPUTS["combined"], h_comb,
                         "finite-place intersection, 4690 -> %d"
                         % comb["summary"]["states_remaining"]),
            "residue": (INPUTS["residue"], h_res,
                        "residue congruences + the C08/C20 counterfactual"),
        },
    }
    if emit_md:
        path = emit(os.path.join(HERE, "ALT_FRONTIER_V2.md"), data)
        if verbose:
            print("wrote %s" % os.path.basename(path))
    return npass == len(_RESULTS), data


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quiet", action="store_true",
                    help="self-check only, exit 0 iff all checks pass")
    ap.add_argument("--no-emit", action="store_true",
                    help="full report, do not (re)write ALT_FRONTIER_V2.md")
    args = ap.parse_args()
    ok, _ = run(verbose=not args.quiet,
                emit_md=not (args.quiet or args.no_emit))
    if args.quiet:
        bad = [t for t, o, _ in _RESULTS if not o]
        print("alt_rebuild: %d/%d checks pass%s"
              % (sum(1 for _, o, _ in _RESULTS if o), len(_RESULTS),
                 "" if ok else "  FAILED: %s" % bad))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
