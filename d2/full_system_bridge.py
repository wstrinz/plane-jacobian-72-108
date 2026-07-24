#!/usr/bin/env python3
"""full_system_bridge.py -- the bridge from a cascade STATE to the FULL
necessary polynomial system (the "G-system"), not merely the f31-alone
condition.

WHY THIS EXISTS
---------------
The cascade / convolution engine tests one necessary condition, f31 = 0.  But
f31 is only the *elimination-ideal* generator of the pre-resultant window system
in the variables (d~2, d~1, d~0, e, Phi) (F37_SATURATION_REPORT.md fact [5]:
E := <G-system> cap Q[d~2,d~1,d~0,e,Phi] = <f31>).  The pre-resultant G-system
carries STRICTLY MORE information: it constrains three further window unknowns
(dm2 = d_-2, dm3 = d_-3, dm4 = d_-4) that f31 has eliminated.  A state that is
solvable under f31-alone is NOT a counterexample germ until it lifts through the
G-system with those spare unknowns realised as honest polynomials of the bounded
degree the window forces.  This module builds exactly that augmented system.

VARIABLE DICTIONARY (see FULL_SYSTEM_BRIDGE.md for the full table)
------------------------------------------------------------------
cascade / window variable  |  G-system indeterminate  |  window meaning  |  k
  d2      (cascade d2)      |  d2                      |  d~2  = d_{4-2}  |  2
  d1      (cascade d1)      |  d1                      |  d~1  = d_{4-3}  |  3
  d0=(d2^2+sigma)/4         |  d0                      |  d~0  = d_{4-4}  |  4
  e       (cascade e)       |  dm1                     |  d_-1 = d_{4-5}  |  5
  --- spare window unknowns the bridge introduces ---
  r                         |  dm2                     |  d_-2 = d_{4-6}  |  6
  s                         |  dm3                     |  d_-3 = d_{4-7}  |  7
  (dm4)                     |  dm4                     |  d_-4 = d_{4-8}  |  8
  Phi = c t^30 q            |  Phi                     |  Phi/y^204       |  -

STRIPPED COORDINATES
--------------------
Every window variable d_{4-k} has order floor 12k (T3_WINDOW_AUDIT.md).  The
cascade works in *stripped* coordinates V_stripped := V_full / y^{12k}
(Phi_stripped := Phi_full / y^204 = c t^30 q, c=-1/6630 -- verify_derivation.py A).
Each G-system generator is weighted-homogeneous under w(d_{4-k}) = 12k,
w(Phi) = 204 (verified: G1,G2,G3,G5 have weights 156,168,180,204), so
    G_i(V_full) = y^{W_i} * G_i(V_stripped),
hence  G_i(full window polys) = 0  <=>  G_i(stripped polys) = 0.  The bridge
therefore substitutes STRIPPED ansaetze everywhere.  A window var of stripped
order-floor 12k and degree cap 15k (sub1) / 14k (sub2) becomes a stripped poly
of order >= 0 and degree <= 3k (sub1) / 2k (sub2).

SOUNDNESS (each added equation is a proven necessary condition)
---------------------------------------------------------------
 * G1,G2,G3,G5body = (D~^3)_{-1,-2,-3,-5} after the (D~^2) linear substitutions,
   i.e. exactly regenerate_system.py's system, validated by T6_SELECTION_AUDIT.md.
 * Phi = c t^30 q is the genuine instance's stripped Phi (verify_derivation.py A).
 * The degree caps for dm2,dm3,dm4 are the k=6,7,8 instances of T3_WINDOW_AUDIT's
   k-INDEPENDENT window induction (deg <= 15k/14k, ord >= 12k).  [JUDGMENT: T3
   states the caps explicitly only for the window k=2..5; the extension to
   k=6,7,8 reuses the identical magic-direction valuation bound -- which is
   k-independent by construction -- together with the D_k in K[y] polynomiality
   that verify_derivation.py section C checks for k down to -13.  Flagged.]

New file, uncommitted.  READ-ONLY on every imported module/artifact.
"""
from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

import convolution_descent as cd
import convolution_elim as ce
import t5_90t1_verify as base
import modular_triage as mt
import system_generators as sysgen

ROOT = Path(__file__).resolve().parent
y = base.y

# G-system indeterminate names (as stored in t4_state.pkl)
D2, D1, D0 = sp.symbols("d2 d1 d0")
DM1, DM2, DM3, DM4 = sp.symbols("dm1 dm2 dm3 dm4")
PHI = sp.symbols("Phi")

C_GENUINE = sp.Rational(-1, 6630)  # Phi = c t^30 q at the genuine (72,108) value

# order floor 12k per window variable d_{4-k}; used both for the homogeneity
# weights and to convert full<->stripped degree caps.
WEIGHT = {"d2": 24, "d1": 36, "d0": 48, "dm1": 60,
          "dm2": 72, "dm3": 84, "dm4": 96, "Phi": 204}

# stripped degree cap = (full deg cap) - (12k):  sub1 -> 3k, sub2 -> 2k.
# k(dm2)=6, k(dm3)=7, k(dm4)=8.
STRIP_DEGCAP = {
    "sub1": {"dm2": 18, "dm3": 21, "dm4": 24},
    "sub2": {"dm2": 12, "dm3": 14, "dm4": 16},
}
SPARE_PREFIX = {"dm2": "R", "dm3": "S", "dm4": "M"}


# --------------------------------------------------------------------------
#  the G-system (loaded, never hand-copied)
# --------------------------------------------------------------------------
def gsystem() -> dict[str, sp.Expr]:
    """The four pre-resultant generators G1,G2,G3,G5body+Phi parsed from the
    canonical generators.json (the same source f37_sat_verify.py reads; no
    pickle is touched)."""
    st = sysgen.load_generators()
    return {
        "G1": st["G1"],
        "G2": st["G2"],
        "G3": st["G3"],
        "G5": st["G5body"] + PHI,
    }


def phi_stripped(c: sp.Expr = C_GENUINE) -> sp.Expr:
    """Stripped Phi = c t^30 q (= Phi_full / y^204)."""
    return sp.expand(c * (y + 1) ** 30 * base.q)


def check_homogeneity() -> dict[str, int]:
    """Assert each generator is weighted-homogeneous; return its weight."""
    weights = {}
    for name, g in gsystem().items():
        seen = set()
        for term in sp.Add.make_args(sp.expand(g)):
            w = 0
            for b_, ex in term.as_powers_dict().items():
                if b_.is_number:
                    continue
                w += WEIGHT[str(b_)] * ex
            seen.add(w)
        assert len(seen) == 1, f"{name} not weighted-homogeneous: {seen}"
        weights[name] = seen.pop()
    assert weights == {"G1": 156, "G2": 168, "G3": 180, "G5": 204}, weights
    return weights


# --------------------------------------------------------------------------
#  spare-unknown ansaetze
# --------------------------------------------------------------------------
def build_spare(regime: str) -> tuple[dict[sp.Symbol, sp.Expr], list[sp.Symbol]]:
    """Generic stripped polynomial ansaetze for dm2,dm3,dm4 under the regime caps."""
    caps = STRIP_DEGCAP[regime]
    polys: dict[sp.Symbol, sp.Expr] = {}
    unknowns: list[sp.Symbol] = []
    for var, name in [(DM2, "dm2"), (DM3, "dm3"), (DM4, "dm4")]:
        deg = caps[name]
        cs = sp.symbols(f"{SPARE_PREFIX[name]}0:{deg + 1}")
        unknowns.extend(cs)
        polys[var] = sum(c * y ** i for i, c in enumerate(cs))
    return polys, unknowns


def ansatz_cost(regime: str) -> dict[str, int]:
    caps = STRIP_DEGCAP[regime]
    per = {n: caps[n] + 1 for n in ("dm2", "dm3", "dm4")}
    per["total_spare"] = sum(per.values())
    return per


# --------------------------------------------------------------------------
#  the augmentation recipe
# --------------------------------------------------------------------------
def augment(ansatz: cd.Ansatz, *, regime: str = "sub2", c: sp.Expr = C_GENUINE,
            nf31: int = 0) -> dict:
    """Given a cascade STATE (a convolution_descent.Ansatz for d2,d1,sigma,e),
    return the FULL necessary polynomial system:

        [ up to nf31 f31 master-coefficient equations of the state ]
      + [ every y-coefficient of G1,G2,G3,(G5body+Phi) on the stripped ansatz ]
      + [ dm2,dm3,dm4 introduced as bounded stripped polynomials ]

    Returns dict with 'equations' (sympy polys, = 0), 'unknowns', 'parameters',
    'spare_unknowns', 'Gpolys', 'f31_coeffs', and diagnostic sizes.
    """
    spare, spare_unk = build_spare(regime)
    subs = {D2: ansatz.d2, D1: ansatz.d1, D0: ansatz.d0, DM1: ansatz.e,
            PHI: phi_stripped(c)}
    subs.update(spare)

    Gpolys: dict[str, sp.Expr] = {}
    eqs: list[sp.Expr] = []
    for name, g in gsystem().items():
        gp = sp.expand(g.subs(subs))
        Gpolys[name] = gp
        if gp != 0:
            for _mono, coeff in sp.Poly(gp, y).terms():
                if coeff != 0:
                    eqs.append(sp.expand(coeff))

    # the state's own f31 master-coefficient equations (redundant -- f31 lies in
    # the G-ideal -- but they are the literal "state's equations" and sharpen GB)
    f31_coeffs: list[sp.Expr] = []
    if nf31 > 0:
        eng = cd.ConvolutionDescent(ansatz, c=c)
        # walk down from the state's top window degree, keep first nf31 nonzero
        top = 245
        found = 0
        d = top
        while found < nf31 and d > top - 40:
            mc = sp.expand(eng.master_coefficient(d))
            if mc != 0:
                f31_coeffs.append(mc)
                found += 1
            d -= 1
    eqs = list(f31_coeffs) + eqs

    params = tuple(sorted(ansatz.parameters, key=sp.default_sort_key))
    state_unk = tuple(ansatz.unknowns)
    unknowns = tuple(sorted(set(state_unk) | set(spare_unk),
                            key=sp.default_sort_key))
    return {
        "equations": eqs,
        "unknowns": unknowns,
        "state_unknowns": state_unk,
        "spare_unknowns": tuple(spare_unk),
        "parameters": params,
        "Gpolys": Gpolys,
        "f31_coeffs": f31_coeffs,
        "regime": regime,
        "c": c,
        "n_equations": len(eqs),
        "n_unknowns": len(unknowns),
    }


# --------------------------------------------------------------------------
#  sound linear pre-elimination (shrinks the variable count before Singular)
# --------------------------------------------------------------------------
def linear_presolve(equations, unknowns):
    """Repeatedly solve any equation that is degree 1 in some unknown with a
    NONZERO RATIONAL-CONSTANT coefficient (no parameter/other-unknown in the
    pivot), and substitute everywhere.  Purely an equivalent rewriting of the
    ideal generators; keeps every coefficient rational.  Returns the reduced
    equation list, the surviving unknowns, and the substitution map."""
    eqs = [sp.expand(e) for e in equations if e != 0]
    remaining = list(unknowns)
    subs_out: dict[sp.Symbol, sp.Expr] = {}
    changed = True
    while changed:
        changed = False
        for idx, eq in enumerate(eqs):
            fs = eq.free_symbols & set(remaining)
            for u in fs:
                p = sp.Poly(eq, u)
                if p.degree() != 1:
                    continue
                lead = p.nth(1)
                if lead.free_symbols or lead == 0:
                    continue  # need a nonzero constant pivot
                val = sp.expand(-p.nth(0) / lead)
                eqs = [sp.expand(e.subs(u, val)) for e in eqs]
                eqs = [e for e in eqs if e != 0]
                for k in list(subs_out):
                    subs_out[k] = sp.expand(subs_out[k].subs(u, val))
                subs_out[u] = val
                remaining.remove(u)
                changed = True
                break
            if changed:
                break
    # dedup
    uniq = []
    seen = set()
    for e in eqs:
        k = sp.srepr(e)
        if k not in seen:
            seen.add(k)
            uniq.append(e)
    return uniq, remaining, subs_out


# --------------------------------------------------------------------------
#  Singular emission (integer-cleared) + runners
# --------------------------------------------------------------------------
def _clear_int(expr: sp.Expr, gens):
    expr = sp.expand(expr)
    if expr == 0:
        return sp.Integer(0)
    p = sp.Poly(expr, *gens)
    L = 1
    for co in p.coeffs():
        L = sp.ilcm(L, sp.Rational(co).q)
    return sp.expand(expr * L)


def _to_singular(expr, gens):
    expr = _clear_int(expr, gens)
    s = sp.sstr(sp.expand(expr)).replace("**", "^").replace(" ", "")
    return s


def singular_program(equations, ring_vars, *, char: int, sat_syms):
    """Build a Singular script over F_char (char>0) or Q (char==0).  Saturates
    by the product of sat_syms (the parameters that must stay nonzero, via a
    Rabinowitsch variable w)."""
    rv = list(ring_vars) + [sp.Symbol("w")]
    var_txt = ",".join(v.name for v in rv)
    lines = ['LIB "elim.lib";', f"ring R = {char},({var_txt}),dp;"]
    members = []
    for i, g in enumerate(equations):
        if char > 0:
            s = mt.poly_to_singular_modp(g, rv, char)
        else:
            s = _to_singular(g, rv)
        if s in ("0", ""):
            continue
        lines.append(f"poly g{i} = {s};")
        members.append(f"g{i}")
    if not members:
        members = ["0"]
    lines.append(f"ideal I = {','.join(members)};")
    if sat_syms:
        prod = "*".join(s.name for s in sat_syms)
        lines.append(f"poly nz = w*{prod}-1;")
        lines.append("ideal Isat = I + ideal(nz);")
        lines.append("ideal G = std(Isat);")
    else:
        lines.append("ideal G = std(I);")
    lines.append("int u = (reduce(1,G)==0);")
    lines.append('"@@UNIT";')
    lines.append("u;")
    lines.append('"@@DIM";')
    lines.append("dim(G);")
    lines.append("quit;")
    return "\n".join(lines) + "\n"


def triage(aug, *, primes=(10007, 10009, 100019), presolve=False, timeout=300.0):
    """mod-p verdict of the augmented system across primes."""
    eqs = list(aug["equations"])
    unk = list(aug["unknowns"])
    params = list(aug["parameters"])
    subs_map = {}
    if presolve:
        eqs, unk, subs_map = linear_presolve(eqs, aug["unknowns"])
    ring_vars = sorted(set(unk) | set(params), key=sp.default_sort_key)
    out = []
    for p in primes:
        prog = singular_program(eqs, ring_vars, char=p, sat_syms=params)
        rr = mt.run_singular(prog, timeout=timeout)
        rr["prime"] = p
        out.append(rr)
        print(f"    p={p}: {rr['verdict']} dim={rr.get('dim')} "
              f"({rr.get('wall')}s)", flush=True)
    pred = mt.classify(out)
    return {"primes": out, "prediction": pred,
            "n_eq_reduced": len(eqs), "n_var_reduced": len(ring_vars),
            "eliminated": len(subs_map)}


def f31_alone_system(ansatz: cd.Ansatz, *, c: sp.Expr = C_GENUINE,
                     ncoeff: int = 16, top: int = 246) -> dict:
    """CONTROL: the f31-ALONE system -- only the state's f31 master
    coefficients (ncoeff of them) plus gamma-saturation.  No G-system, no spare
    unknowns.  This is what the cascade tests; used to show the full system
    kills where f31-alone (at reachable depth) does not."""
    eng = cd.ConvolutionDescent(ansatz, c=c)
    coeffs, d = [], top
    while len(coeffs) < ncoeff and d > top - 60:
        mc = sp.expand(eng.master_coefficient(d))
        if mc != 0:
            coeffs.append(mc)
        d -= 1
    params = tuple(sorted(ansatz.parameters, key=sp.default_sort_key))
    return {"equations": coeffs, "unknowns": tuple(ansatz.unknowns),
            "parameters": params, "n_equations": len(coeffs),
            "n_unknowns": len(ansatz.unknowns)}


def exact_kill(aug, *, presolve=False, timeout=420.0):
    """Attempt an exact CONTRADICTION over Q (UNIT ideal => state killed)."""
    eqs = list(aug["equations"])
    unk = list(aug["unknowns"])
    params = list(aug["parameters"])
    if presolve:
        eqs, unk, _ = linear_presolve(eqs, aug["unknowns"])
    ring_vars = sorted(set(unk) | set(params), key=sp.default_sort_key)
    prog = singular_program(eqs, ring_vars, char=0, sat_syms=params)
    rr = mt.run_singular(prog, timeout=timeout)
    return rr


# --------------------------------------------------------------------------
#  pilot: one a8 constant-E stall state
# --------------------------------------------------------------------------
def pilot_state(index: int = 0):
    """The a8 constant-E UNRESOLVED states of batch_convolution_sub2.json
    (System 2 of MODULAR_TRIAGE).  index 0 = the simplest (d2=0, deg d1=0,
    deg sigma=5)."""
    d = json.load(open(ROOT / "batch_convolution_sub2.json"))
    states = [s for s in d["states"]
              if s["a_t"] == 8 and s["branch"] == "T1"
              and s["deg_e"] == 8 and s["final_verdict"] == "UNRESOLVED"]
    s = states[index]
    gamma = sp.Symbol("gamma")
    e = gamma * (y + 1) ** 8
    degrees = {"d1": int(s["deg_d1"]), "sigma": int(s["deg_sigma"])}
    if s["d2_zero"]:
        d2arg = {"d2": sp.Integer(0)}
    else:
        degrees["d2"] = int(s["deg_d2"])
        d2arg = {}
    ans = cd.build_ansatz(e=e, degrees=degrees, parameters=(gamma,), **d2arg)
    label = f"a8_dd2{s['deg_d2']}_dd1{s['deg_d1']}_dsig{s['deg_sigma']}"
    return ans, s, label


def run_pilot(index: int = 0, regime: str = "sub2"):
    weights = check_homogeneity()
    print(f"G-system weights (homogeneity check): {weights}")
    print(f"ansatz cost sub1: {ansatz_cost('sub1')}")
    print(f"ansatz cost sub2: {ansatz_cost('sub2')}")
    ans, meta, label = pilot_state(index)
    print(f"\nPILOT state: {label}  (regime {regime})")
    print(f"  e={ans.e}  d1={ans.d1}  sigma={ans.sigma}  d2={ans.d2}")
    print(f"  f31-alone residual: {meta['gauge_detail']['residual_factored']}")

    # CONTROL: f31-alone (16 master coefficients) -- what the cascade sees
    print("\n  CONTROL: f31-alone (16 master coeffs) mod-p triage:")
    f31 = f31_alone_system(ans)
    tri_f31 = triage(f31)
    print(f"  => f31-alone: {tri_f31['prediction']}")

    aug = augment(ans, regime=regime)
    print(f"\n  augmented: {aug['n_equations']} eqs, {aug['n_unknowns']} unknowns "
          f"+ params {aug['parameters']}; spare {len(aug['spare_unknowns'])}")

    print("  mod-p triage of the FULL system:")
    tri = triage(aug)
    print(f"  => full-system: {tri['prediction']}")

    exact = None
    if tri["prediction"] == "LIKELY-EMPTY":
        print("\n  UNIT predicted -> attempting exact kill over Q:")
        exact = exact_kill(aug)
        print(f"  exact-Q verdict: {exact['verdict']} ({exact.get('wall')}s)")

    result = {
        "label": label, "regime": regime,
        "meta": {k: meta[k] for k in ("deg_d1", "deg_sigma", "deg_d2",
                                      "d2_zero", "final_verdict")},
        "f31_residual": meta["gauge_detail"]["residual_factored"],
        "weights": weights,
        "ansatz_cost": {"sub1": ansatz_cost("sub1"), "sub2": ansatz_cost("sub2")},
        "n_equations": aug["n_equations"], "n_unknowns": aug["n_unknowns"],
        "f31_alone_modp": tri_f31, "modp": tri, "exact": exact,
    }
    json.dump(result, open(ROOT / "full_system_bridge_pilot.json", "w"),
              indent=1, default=str)
    print("\nwrote full_system_bridge_pilot.json")
    return result


if __name__ == "__main__":
    import sys
    idx = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    reg = sys.argv[2] if len(sys.argv) > 2 else "sub2"
    run_pilot(idx, reg)
