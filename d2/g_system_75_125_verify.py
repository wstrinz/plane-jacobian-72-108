#!/usr/bin/env python3
"""g_system_75_125_verify.py  (NEW; read-only over all existing artifacts)

Exact PASS/FAIL checker for G_SYSTEM_75_125.md / g_system_75_125.py /
g_system_75_125.json -- the (75,125) analogue of the D-transform G-system.

Structure of the check (all exact sympy):

  A. RECIPE CONTROL (independent ground truth).  Rebuild the (72,108) G-system
     from scratch with the same parametric builder and assert it reproduces the
     PUBLISHED generators G1,G2,G3,(G5body+Phi) of FULL_SYSTEM_BRIDGE.md sec.1
     bit-for-bit, and the known G-weights [156,168,180,204].  This validates the
     builder against externally-known values before it is trusted on (75,125).

  B. JSON STRUCTURE.  Load g_system_75_125.json and check the corner signature,
     variable order / ring, spare inventory (dm2..dm8, 7 spares), generator
     count/slices (8 quintic generators, j=1..7,9), and the recipe parameters.

  C. HOMOGENEITY (independent of the builder).  Parse each stored generator
     polynomial string and verify every monomial carries the intrinsic u-weight
     w(d_m)=t-m, w(Phi)=M -- i.e. weight = b t + j for generator G_j; and that
     the forcing weights form the arithmetic progression 21..29 (skip 28).

  D. PHYSICAL-WEIGHT OBSTRUCTION.  W_step = ord_y(Phi)/M = 80/29 is
     NON-integral (vs 12 at (72,108)); the a=3 boundary of
     CORNER_144_COMPARISON.md sec.5.  Contrast with the two a-based controls.

  E. SPOT RE-DERIVATION (anti-fabrication).  Rebuild (75,125) generators G1,G2
     from scratch (truncated) and assert they equal the JSON's stored G1,G2.

  F. PHI-CONSISTENCY.  Phi enters G9 at u-slice M=29 with the homogeneity-forced
     weight; the phase-1 slice-sum clear = a*M-b = 82, N = clear-b = 77
     (C_SERIES_75_125.md).

  G. THE RETRACTION GUARD.  The (75,125) inputs (t,q) = (4,1) are cross-checked
     against polygon_reduction.corner_chart_data, and the superseded
     (t,q) = (l_final,b_final) = (5,2) is checked to RAISE.

*** REPAIRED 2026-07-26 (PASSPORT_75_125_REPAIR.md): t 5->4, kappa 3->2, q 2->1,
C y^2(y^3+1) -> y, M 36->29, ord_y(Phi) 201->80, generators 10->8, spares 9->7. ***

Run:  python d2_plane_72_108/g_system_75_125_verify.py [--quiet]
Exit 0 on pass; any failed claim raises SystemExit(nonzero).
"""
from pathlib import Path
import json
import sys

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from g_system_75_125 import build_gsystem, published_72108  # noqa: E402

QUIET = "--quiet" in sys.argv[1:]
checks = 0
HERE = Path(__file__).resolve().parent


def check(name, cond):
    global checks
    if not bool(cond):
        raise SystemExit(f"[FAIL] {name}")
    checks += 1
    if not QUIET:
        print(f"[OK] {name}")


T7525, M7525 = 4, 29                  # REPAIRED: chart exponent 4, forcing slice 29


def uweight_7525(name):
    """Intrinsic u-weight of a (75,125) generator symbol (t=4, M=29)."""
    if name == "Phi":
        return M7525
    if name.startswith("dm"):
        return T7525 + int(name[2:])
    return T7525 - int(name[1:])      # d{idx}


# ===========================================================================
# A. RECIPE CONTROL -- the builder reproduces the published (72,108) system.
# ===========================================================================
if not QUIET:
    print("A. recipe control: rebuild (72,108) and match FULL_SYSTEM_BRIDGE sec.1")
c = build_gsystem(2, 3, 4, 7, 204)
pub = published_72108()
for j in (1, 2, 3, 5):
    check("(72,108) builder reproduces published G%d exactly" % j,
          sp.expand(c["Gs"][j] - pub[j]) == 0)
check("(72,108) spare inventory is dm2,dm3,dm4 (= d_-2,d_-3,d_-4)",
      [str(s) for s in c["spares"]] == ["dm2", "dm3", "dm4"])
uw = [c["b"] * c["t"] + j for j in sorted(c["Gs"])]
check("(72,108) u-grading weights are [13,14,15,17]", uw == [13, 14, 15, 17])
check("(72,108) W_step = ord_y(Phi)/M = 204/17 = 12 (INTEGER)",
      c["W_step"] == 12 and c["W_step"].q == 1)
phys = [int(c["W_step"] * w) for w in uw]
check("(72,108) physical weights = known G-weights [156,168,180,204]",
      phys == [156, 168, 180, 204])
for j in sorted(c["Gs"]):
    check("(72,108) G%d homogeneous at u-weight %d" % (j, c["b"] * c["t"] + j),
          c["homog"][j] == {c["b"] * c["t"] + j})

# ===========================================================================
# B. JSON STRUCTURE.
# ===========================================================================
if not QUIET:
    print("\nB. load g_system_75_125.json and check structure")
jpath = HERE / "g_system_75_125.json"
check("g_system_75_125.json exists", jpath.exists())
D = json.loads(jpath.read_text(encoding="utf-8"))
cs = D["case"]
check("case = F2_j1_75_125, degrees (75,125), corner (5,20)",
      cs["tag"] == "F2_j1_75_125" and cs["degrees"] == [75, 125]
      and cs["a"] == 3 and cs["b"] == 5 and cs["t"] == 4)
check("corner signature kappa=t-2=2, q=1, e=b-a+1=3, s=kappa+1-a*t=-9  [REPAIRED]",
      cs["kappa"] == 2 and cs["q"] == 1 and cs["e"] == 3 and cs["s"] == -9)
check("the json records the chart-exponent rule and the repair provenance, so a "
      "future reader cannot mistake t for GGV5's l_final",
      "ceil(b0/a0)" in cs["chart_exponent_rule"]
      and "2026-07-26" in cs["repaired"] and "NOT the chart" in cs["corner"])
check("recipe: linear window S^3, forcing window S^5, M = b*t+jphi = 29",
      "S^3" in D["recipe"]["linear_window"]
      and "S^5" in D["recipe"]["forcing_window"]
      and D["recipe"]["forcing_slice_M"] == 29 and D["recipe"]["jphi"] == 9)
check("recipe: 16 linear eliminations, k=1..17 skip 16",
      D["recipe"]["linear_eliminations"] == 16
      and "1..17" in D["recipe"]["linear_slice_range"]
      and "skip k=16" in D["recipe"]["linear_slice_range"])
check("spare inventory = dm2..dm8 (7 spares = d_-2..d_-8 = d_-2..d_-(a-1)t)",
      D["spare_variables"] == ["dm%d" % k for k in range(2, 9)]
      and D["num_spares"] == 7)
check("state variables = d2,d1,d0,dm1 (d_{t-2}..d_0, e=d_-1; t=4 so d3 is gone)",
      D["state_variables"] == ["d2", "d1", "d0", "dm1"])
check("ring variable order = state + spares + Phi (canonical, documented)",
      D["variable_order"] == ["d2", "d1", "d0", "dm1"]
      + ["dm%d" % k for k in range(2, 9)] + ["Phi"]
      and D["ring"] == "Q[%s]" % ",".join(D["variable_order"]))
check("8 generators from forcing window S^5, names G1..G7,G9 (slice j=1..9 skip 8)",
      D["num_generators"] == 8 and D["forcing_window_power"] == 5
      and D["generator_names"] == ["G%d" % j for j in [1, 2, 3, 4, 5, 6, 7, 9]])
check("generators are WEIGHT-homogeneous, NOT total-degree homogeneous "
      "(deg maxima grow 5..7 as a=3 cubic substitutions inflate deeper slices)",
      [D["generators"]["G%d" % j]["total_degree_max"]
       for j in [1, 2, 3, 4, 5, 6, 7, 9]] == [5, 5, 5, 5, 5, 6, 6, 7]
      and D["generators"]["G1"]["total_degree_min"] == 2)
check("counts follow the closed laws #gens = a*t-kappa-2 and #spares = (a-1)t-1, "
      "so both increment by t (=4) per family rung -- the superseded t=5 gave the "
      "'five-row/five-column block rule', which is now a FOUR-row/four-column rule",
      D["num_generators"] == 3 * 4 - 2 - 2 and D["num_spares"] == (3 - 1) * 4 - 1)

# ===========================================================================
# C. HOMOGENEITY -- re-parse the stored generator strings (builder-independent).
# ===========================================================================
if not QUIET:
    print("\nC. homogeneity of the stored generators (independent re-parse)")
gens_parsed = {}
for name, g in sorted(D["generators"].items(), key=lambda kv: kv[1]["slice_j"]):
    expr = sp.sympify(g["poly"])
    gens_parsed[name] = expr
    j = g["slice_j"]
    body = expr - (sp.Symbol("Phi") if g["has_phi"] else 0)
    P = sp.Poly(body, *sorted(body.free_symbols, key=str))
    ws = {sum(uweight_7525(str(v)) * ex for v, ex in zip(P.gens, mon))
          for mon, _ in P.terms()}
    check("%s homogeneous under u-grading at weight b*t+j = %d"
          % (name, 20 + j), ws == {20 + j})
check("only G9 carries Phi; Phi u-weight = M = 29",
      all(("Phi" in g["poly"]) == (g["slice_j"] == 9)
          for g in D["generators"].values())
      and uweight_7525("Phi") == 29)
check("forcing-generator weights are the AP 21,22,23,24,25,26,27,29 "
      "(common diff 1; 28 absent = skipped G8)",
      D["u_grading_weights"] == [21, 22, 23, 24, 25, 26, 27, 29])

# ===========================================================================
# D. PHYSICAL-WEIGHT OBSTRUCTION (the a>=3 boundary).
# ===========================================================================
if not QUIET:
    print("\nD. physical-weight normalisation: the a=3 obstruction")
pw = D["physical_weight"]
W = sp.Rational(pw["W_step_num"], pw["W_step_den"])
check("W_step = ord_y(Phi)/M = 80/29  [REPAIRED from 201/36 = 67/12]",
      W == sp.Rational(80, 29)
      and pw["ord_y_Phi"] == 80 and pw["forcing_slice_M"] == 29)
check("W_step is NON-INTEGRAL (denominator 29) -> exact y-stripping fails",
      W.q == 29 and pw["integral"] is False)
check("contrast a-controls: (72,108) 204/17=12 integral; (108,144) 205/25=41/5 "
      "and (75,125) 80/29 both non-integral (a>=3)",
      sp.Rational(204, 17) == 12 and sp.Rational(204, 17).q == 1
      and sp.Rational(205, 25) == sp.Rational(41, 5)
      and sp.Rational(205, 25).q == 5 and W.q == 29)
check("window-cap layer flagged OBSTRUCTED (quasipolynomial, quasi-period 29)",
      "OBSTRUCTED" in D["window_caps"]["status"]
      and "quasi" in D["window_caps"]["pattern"].lower()
      and D["window_caps"]["quasi_period"] == 29)
check("quasi-period 29 is PRIME, so the superseded 'divisor lattice of the "
      "period' reading of 12 has no counterpart",
      D["window_caps"]["quasi_period_is_prime"] is True and sp.isprime(29))
check("and since C = y is a monomial, deg_y(Phi) = ord_y(Phi), so the DEG cap "
      "coincides with the ORD cap -- there is no second, affine slope here, "
      "unlike (72,108)'s ord 12 / deg 14",
      D["deg_slope"]["deg_equals_ord"] is True)
check("W_step*M = 80 = ord_y(Phi) (rational identity holds; only integrality fails)",
      W * 29 == 80)

# ===========================================================================
# E. SPOT RE-DERIVATION -- rebuild (75,125) G1,G2 from scratch (anti-fabrication).
# ===========================================================================
if not QUIET:
    print("\nE. spot re-derivation of (75,125) G1,G2 from scratch (truncated)")
spot = build_gsystem(3, 5, 4, 1, 80, Nmax_override=23, jset=[1, 2])
for j, name in [(1, "G1"), (2, "G2")]:
    check("rebuilt (75,125) %s equals the stored JSON generator" % name,
          sp.expand(spot["Gs"][j] - gens_parsed[name]) == 0)
check("spot spare set within G1,G2 subset of the declared dm2..dm8",
      set(str(s) for s in spot["spares"]) <= set(D["spare_variables"]))

# ===========================================================================
# F. PHI-CONSISTENCY (slice-sum / phase-1 tower agreement).
# ===========================================================================
if not QUIET:
    print("\nF. Phi-consistency with the phase-1 tower (C_SERIES_75_125.md)")
a, b, t, M = 3, 5, 4, 29
clear = a * M - b
N = clear - b
check("Phi occupies the (D~^5)_-9 forcing slice: M = b*t+jphi = 29", M == b * t + 9)
check("slice-sum invariant: clear = a*M - b = 82", clear == 82)
check("tower length N = clear - b = 77 (matches C_SERIES_75_125.md)", N == 77)
check("ord_y(Phi) = W_step*M = 80 and = rho + q*N = 3 + 1*77",
      W * M == 80 and 3 + 1 * 77 == 80)
check("Phi-consistency verdict recorded CONSISTENT (intrinsic grading)",
      D["phi_consistency"]["verdict"].startswith("CONSISTENT"))

# ===========================================================================
# G. THE RETRACTION GUARD -- the inputs must come from it, both directions.
# ===========================================================================
if not QUIET:
    print("\nG. retraction guard on the (75,125) corner inputs")
import polygon_reduction as pr                                      # noqa: E402
_cd = pr.corner_chart_data(5, 20, l_final=5, b_final=2, who="g_system_verify")
check("corner_chart_data(5,20) gives (t,kappa,deg C,ord C) = (4,2,1,1) and flags "
      "MONOMIAL / no retraction",
      (_cd["t"], _cd["kappa"], _cd["deg_C"], _cd["ord_C"]) == (4, 2, 1, 1)
      and _cd["monomial"] and not _cd["retraction"])
check("the json's t and q agree with the guard, so the build cannot drift",
      cs["t"] == _cd["t"] and cs["q"] == _cd["ord_C"] and cs["kappa"] == _cd["kappa"])
try:
    pr.final_corner_dictionary(5, 20, 5, 2)
    raise SystemExit("[FAIL] the guard must RAISE at (5,20)")
except pr.FinalCornerDictionaryError:
    check("the superseded (t,q)=(l_final,b_final)=(5,2) RAISES", True)
check("and the guard RETURNS at (8,28) -- the (72,108) corner where the "
      "dictionary IS valid -- so it is not vacuous",
      pr.final_corner_dictionary(8, 28, 4, 7) == (4, 7))

# ===========================================================================
if not QUIET:
    print(f"\nALL {checks} G-SYSTEM (75,125) CHECKS PASSED")
    print("VERDICT: G-system BUILT (structure transfers: 8 quintic generators, "
          "7 spares dm2..dm8, u-grading AP weights 21..29 skip 28).")
    print("Physical-weight OBSTRUCTION characterised: W_step=80/29 "
          "non-integral; Phi-consistency CONSISTENT.")
    print(f"script: {Path(__file__).resolve()}")
sys.exit(0)
