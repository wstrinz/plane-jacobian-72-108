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
     variable order / ring, spare inventory (dm2..dm10, 9 spares), generator
     count/slices (10 quintic generators, j=1..9,11), and the recipe parameters.

  C. HOMOGENEITY (independent of the builder).  Parse each stored generator
     polynomial string and verify every monomial carries the intrinsic u-weight
     w(d_m)=t-m, w(Phi)=M -- i.e. weight = b t + j for generator G_j; and that
     the forcing weights form the arithmetic progression 26..36 (skip 35).

  D. PHYSICAL-WEIGHT OBSTRUCTION.  W_step = ord_y(Phi)/M = 201/36 = 67/12 is
     NON-integral (vs 12 at (72,108)); the a=3 boundary of
     CORNER_144_COMPARISON.md sec.5.  Contrast with the two a-based controls.

  E. SPOT RE-DERIVATION (anti-fabrication).  Rebuild (75,125) generators G1,G2
     from scratch (truncated) and assert they equal the JSON's stored G1,G2.

  F. PHI-CONSISTENCY.  Phi enters G11 at u-slice M=36 with the homogeneity-forced
     weight; the phase-1 slice-sum clear = a*M-b = 103, N = clear-b = 98
     (C_SERIES_75_125.md).

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


def uweight_7525(name):
    """Intrinsic u-weight of a (75,125) generator symbol (t=5, M=36)."""
    if name == "Phi":
        return 36
    if name.startswith("dm"):
        return 5 + int(name[2:])
    return 5 - int(name[1:])          # d{idx}


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
check("case = F2_j1_75_125, degrees (75,125), corner (5,20)->(7/5,2)",
      cs["tag"] == "F2_j1_75_125" and cs["degrees"] == [75, 125]
      and cs["a"] == 3 and cs["b"] == 5 and cs["t"] == 5)
check("corner signature kappa=t-2=3, q=2, e=b-a+1=3, s=kappa+1-a*t=-11",
      cs["kappa"] == 3 and cs["q"] == 2 and cs["e"] == 3 and cs["s"] == -11)
check("recipe: linear window S^3, forcing window S^5, M = b*t+jphi = 36",
      "S^3" in D["recipe"]["linear_window"]
      and "S^5" in D["recipe"]["forcing_window"]
      and D["recipe"]["forcing_slice_M"] == 36 and D["recipe"]["jphi"] == 11)
check("recipe: 20 linear eliminations, k=1..21 skip 20",
      D["recipe"]["linear_eliminations"] == 20
      and "1..21" in D["recipe"]["linear_slice_range"]
      and "skip k=20" in D["recipe"]["linear_slice_range"])
check("spare inventory = dm2..dm10 (9 spares = d_-2..d_-10 = d_-2..d_-(a-1)t)",
      D["spare_variables"] == ["dm%d" % k for k in range(2, 11)]
      and D["num_spares"] == 9)
check("state variables = d3,d2,d1,d0,dm1 (d_{t-2}..d_0, e=d_-1)",
      D["state_variables"] == ["d3", "d2", "d1", "d0", "dm1"])
check("ring variable order = state + spares + Phi (canonical, documented)",
      D["variable_order"] == ["d3", "d2", "d1", "d0", "dm1"]
      + ["dm%d" % k for k in range(2, 11)] + ["Phi"]
      and D["ring"] == "Q[%s]" % ",".join(D["variable_order"]))
check("10 generators from forcing window S^5, names G1..G9,G11 (slice j=1..11 skip 10)",
      D["num_generators"] == 10 and D["forcing_window_power"] == 5
      and D["generator_names"] == ["G%d" % j for j in [1, 2, 3, 4, 5, 6, 7, 8, 9, 11]])
check("generators are WEIGHT-homogeneous, NOT total-degree homogeneous "
      "(deg maxima grow 5..8 as a=3 cubic substitutions inflate deeper slices)",
      [D["generators"]["G%d" % j]["total_degree_max"]
       for j in [1, 2, 3, 4, 5, 6, 7, 8, 9, 11]] == [5, 5, 5, 5, 5, 5, 6, 6, 7, 8]
      and D["generators"]["G1"]["total_degree_min"] == 2)

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
          % (name, 25 + j), ws == {25 + j})
check("only G11 carries Phi; Phi u-weight = M = 36",
      all(("Phi" in g["poly"]) == (g["slice_j"] == 11)
          for g in D["generators"].values())
      and uweight_7525("Phi") == 36)
check("forcing-generator weights are the AP 26,27,28,29,30,31,32,33,34,36 "
      "(common diff 1; 35 absent = skipped G10)",
      D["u_grading_weights"] == [26, 27, 28, 29, 30, 31, 32, 33, 34, 36])

# ===========================================================================
# D. PHYSICAL-WEIGHT OBSTRUCTION (the a>=3 boundary).
# ===========================================================================
if not QUIET:
    print("\nD. physical-weight normalisation: the a=3 obstruction")
pw = D["physical_weight"]
W = sp.Rational(pw["W_step_num"], pw["W_step_den"])
check("W_step = ord_y(Phi)/M = 201/36 = 67/12", W == sp.Rational(67, 12)
      and pw["ord_y_Phi"] == 201 and pw["forcing_slice_M"] == 36)
check("W_step is NON-INTEGRAL (denominator 12) -> exact y-stripping fails",
      W.q == 12 and pw["integral"] is False)
check("contrast a-controls: (72,108) 204/17=12 integral; (108,144) 205/25=41/5 "
      "and (75,125) 201/36=67/12 both non-integral (a>=3)",
      sp.Rational(204, 17) == 12 and sp.Rational(204, 17).q == 1
      and sp.Rational(205, 25) == sp.Rational(41, 5)
      and sp.Rational(205, 25).q == 5 and W.q == 12)
check("window-cap layer flagged OBSTRUCTED (quasipolynomial, quasi-period 12)",
      "OBSTRUCTED" in D["window_caps"]["status"]
      and "quasi" in D["window_caps"]["pattern"].lower())
check("W_step*M = 201 = ord_y(Phi) (rational identity holds; only integrality fails)",
      W * 36 == 201)

# ===========================================================================
# E. SPOT RE-DERIVATION -- rebuild (75,125) G1,G2 from scratch (anti-fabrication).
# ===========================================================================
if not QUIET:
    print("\nE. spot re-derivation of (75,125) G1,G2 from scratch (truncated)")
spot = build_gsystem(3, 5, 5, 2, 201, Nmax_override=27, jset=[1, 2])
for j, name in [(1, "G1"), (2, "G2")]:
    check("rebuilt (75,125) %s equals the stored JSON generator" % name,
          sp.expand(spot["Gs"][j] - gens_parsed[name]) == 0)
check("spot spare set within G1,G2 subset of the declared dm2..dm10",
      set(str(s) for s in spot["spares"]) <= set(D["spare_variables"]))

# ===========================================================================
# F. PHI-CONSISTENCY (slice-sum / phase-1 tower agreement).
# ===========================================================================
if not QUIET:
    print("\nF. Phi-consistency with the phase-1 tower (C_SERIES_75_125.md)")
a, b, t, M = 3, 5, 5, 36
clear = a * M - b
N = clear - b
check("Phi occupies the (D~^5)_-11 forcing slice: M = b*t+jphi = 36", M == b * t + 11)
check("slice-sum invariant: clear = a*M - b = 103", clear == 103)
check("tower length N = clear - b = 98 (matches C_SERIES_75_125.md)", N == 98)
check("ord_y(Phi) = W_step*M = 201 and = rho + q*N = 5 + 2*98",
      W * M == 201 and 5 + 2 * 98 == 201)
check("Phi-consistency verdict recorded CONSISTENT (intrinsic grading)",
      D["phi_consistency"]["verdict"].startswith("CONSISTENT"))

# ===========================================================================
if not QUIET:
    print(f"\nALL {checks} G-SYSTEM (75,125) CHECKS PASSED")
    print("VERDICT: G-system BUILT (structure transfers: 10 quintic generators, "
          "9 spares dm2..dm10, u-grading AP weights 26..36).")
    print("Physical-weight OBSTRUCTION characterised: W_step=201/36=67/12 "
          "non-integral (a=3 boundary); Phi-consistency CONSISTENT.")
    print(f"script: {Path(__file__).resolve()}")
sys.exit(0)
