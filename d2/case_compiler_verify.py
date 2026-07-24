#!/usr/bin/env python3
"""case_compiler_verify.py — PASS/FAIL checker for the Lane I case-compiler pilot.

Validates compiled dossiers against KNOWN LANDED VALUES hard-coded here from the
source documents (STATE.md audited facts, PHI_75_125.md, PHI_CORNER4.md,
PHI_F14.md, GALOIS_LIBRARY.md, FULL_SYSTEM_BRIDGE.md) — the verifier does not
trust the compiler's own registries for the expected constants.

Checks:
  A. the unified corner law reproduces all SEVEN landed points exactly;
  B/C/D. the three pilot dossiers ((72,108), (75,125)=F2 j=1, F9 j=0) field by
     field: corner signature, law signature, regime + conjectural flags,
     Galois transfer (label / disc class / verdicts / conditionality),
     presentation instantiation markers, master-identity parameters at (72,108);
  E. the Galois labeling routines on independent witnesses (S4 quartic with
     disc class 17; 10th cyclotomic -> C4 with a constructive order-4
     automorphism cross-check; V4 witness y^4+1; square-class arithmetic
     including squarefree(170*17) = 10);
  F. survey-wide sanity: Diophantine identity for every length-1 family at its
     smallest coprime j; conjectural flags fire exactly where they must
     (F12 A0'!=(1,0); F4 k=2), and F7 is grounded by the PHI_F7.md ramified
     law reproducing its derived signature;
  G. canonical JSON determinism + the three on-disk pilot dossiers match
     recomputation.

Usage: python case_compiler_verify.py [--quiet]     exit 0 iff all pass.
"""
import json
import sys
from math import gcd

import sympy as sp

import case_compiler as cc

y = sp.symbols("y")
QUIET = "--quiet" in sys.argv
count = 0

def ok(label, condition):
    global count
    if not condition:
        print("[FAIL] %s" % label)
        sys.exit(1)
    count += 1
    if not QUIET:
        print("[OK] %s" % label)

# ---------------------------------------------------------------------------
# A. the law reproduces all seven landed points (constants from the PHI docs)
# ---------------------------------------------------------------------------
SEVEN = [
    # (label, a, b, t, a0, q, expected signature)  — PHI_F14.md table
    ("(72,108)",  2, 3, 4, 8, 7, (238, 204, 30, 4)),
    ("(108,144)", 3, 4, 4, 8, 3, (550, 205, 69, 276)),
    ("(75,125)",  3, 5, 5, 5, 2, (504, 201, 101, 202)),
    ("(56,84)",   2, 3, 7, 7, 2, (377, 107, 54, 216)),
    ("(50,75)",   2, 3, 5, 5, 2, (189, 75, 38, 76)),
    ("(66,231)",  2, 7, 3, 9, 4, (375, 165, 42, 168)),
    ("(48,64)",   3, 4, 4, 4, 3, (275, 205, 69, 1)),
]
for label, a, b, t, a0, q, want in SEVEN:
    s = cc.corner_signature(a, b, t, a0, q)
    L = cc.law_signature(s)
    got = (L["deg"], L["ord_y"], L["mult_y_plus_1"], L["cofactor_deg"])
    ok("A: law reproduces %s -> %s" % (label, want), got == want)

# ---------------------------------------------------------------------------
# B. (72,108) dossier — audited home case
# ---------------------------------------------------------------------------
d72 = cc.compile_case("GGHV_72_108")
s = d72["corner_signature"]
ok("B: (72,108) corner signature",
   (s["a"], s["b"], s["t"], s["kappa"], s["a0"], s["q"]) == (2, 3, 4, 2, 8, 7)
   and (s["e"], s["r"], s["gap"], s["dg"], s["N"]) == (2, 0, 4, 1, 28))
sig = d72["phi_prediction"]["signature"]
ok("B: (72,108) law signature audited (238,204,30,4)",
   (sig["deg"], sig["ord_y"], sig["mult_y_plus_1"], sig["cofactor_deg"])
   == (238, 204, 30, 4))
ok("B: (72,108) regime resonance_gap_r0, grounded, NOT conjectural",
   d72["phi_prediction"]["regime"] == "resonance_gap_r0"
   and d72["phi_prediction"]["regime_grounded"] is True
   and d72["phi_prediction"]["conjectural"] is False)
g = d72["galois_transfer"]
ok("B: (72,108) forcing quartic is the audited 2048y^4-512y^3+320y^2-240y+195",
   sp.expand(sp.sympify(g["forcing_candidate"]["poly"])
             - (2048*y**4 - 512*y**3 + 320*y**2 - 240*y + 195)) == 0)
ok("B: (72,108) Galois label S4, disc class 17",
   g["galois"]["label"] == "S4" and g["galois"]["disc_class"] == 17)
ok("B: (72,108) C08 and C20 both KILL, status AUDITED",
   g["verdicts"]["C08"] == "KILLS" and g["verdicts"]["C20"] == "KILLS"
   and g["status"].startswith("AUDITED"))
pres = d72["presentations"]
ok("B: (72,108) BOTH presentations INSTANTIATED",
   pres["eliminated_f31_style"]["instantiated"] is True
   and pres["pre_resultant_G_system"]["instantiated"] is True)
ok("B: (72,108) master identity (D,F,E,s) = (31,7,21,3)",
   pres["eliminated_f31_style"]["parameters"] == dict(D=31, F=7, E=21, s=3))
ok("B: (72,108) G-weights 156/168/180/204 and 4 generators",
   pres["pre_resultant_G_system"]["G_weights"] == [156, 168, 180, 204]
   and pres["pre_resultant_G_system"]["generators"]
   == ["G1", "G2", "G3", "G5body+Phi"])
pn = pres["eliminated_f31_style"]["phi_normalized"]
ok("B: (72,108) Phi_stripped = t^30 * (deg-4 unit), lc -1/6630, total deg 34",
   pn["mult"] == 30 and pn["deg_u_gap"] == 4 and pn["H"] == "1"
   and pn["lc"] == "-1/6630" and pn["total_stripped_degree"] == 34)

# ---------------------------------------------------------------------------
# C. (75,125) = F2 j=1
# ---------------------------------------------------------------------------
d75 = cc.compile_case("F2", 1)
ok("C: F2 j=1 is the (75,125) case", d75["case"]["degrees"] == [75, 125]
   and (d75["case"]["m"], d75["case"]["n"]) == (3, 5))
s = d75["corner_signature"]
ok("C: (75,125) corner signature",
   (s["a"], s["b"], s["t"], s["a0"], s["q"], s["e"], s["r"], s["gap"],
    s["dg"], s["N"]) == (3, 5, 5, 5, 2, 3, 2, 0, 3, 98))
sig = d75["phi_prediction"]["signature"]
ok("C: (75,125) law signature (504,201,101,202)",
   (sig["deg"], sig["ord_y"], sig["mult_y_plus_1"], sig["cofactor_deg"])
   == (504, 201, 101, 202))
ok("C: (75,125) grounded regime, not conjectural, Diophantine OK",
   d75["phi_prediction"]["conjectural"] is False
   and d75["case"]["diophantine"].endswith("OK"))
g = d75["galois_transfer"]
ok("C: (75,125) forcing residual H = y^2-y+1 (PHI_75_125 H_2), label C2, disc class -3",
   sp.expand(sp.sympify(g["forcing_candidate"]["poly"]) - (y**2 - y + 1)) == 0
   and g["galois"]["label"] == "C2" and g["galois"]["disc_class"] == -3)
ok("C: (75,125) kills transfer (C08, C20) but CONDITIONAL",
   g["verdicts"]["C08"] == "KILLS" and g["verdicts"]["C20"] == "KILLS"
   and g["status"].startswith("CONDITIONAL"))
ok("C: (75,125) presentations schematic (no f31-analogue exists)",
   d75["presentations"]["eliminated_f31_style"]["instantiated"] is False
   and d75["presentations"]["pre_resultant_G_system"]["instantiated"] is False
   and d75["presentations"]["eliminated_f31_style"]["parameters"] is None)
ok("C: (75,125) N-formula judgment records the broken slice picture ((b-1)/a=4/3)",
   any("is NOT integral" in j for j in d75["judgment"]))

# ---------------------------------------------------------------------------
# D. F9 j=0 (56,84)
# ---------------------------------------------------------------------------
d56 = cc.compile_case("F9", 0)
ok("D: F9 j=0 is the (56,84) case", d56["case"]["degrees"] == [56, 84])
s = d56["corner_signature"]
ok("D: (56,84) corner signature",
   (s["a"], s["b"], s["t"], s["a0"], s["q"], s["e"], s["r"], s["gap"],
    s["dg"], s["N"]) == (2, 3, 7, 7, 2, 2, 4, 0, 5, 52))
sig = d56["phi_prediction"]["signature"]
ok("D: (56,84) law signature (377,107,54,216)",
   (sig["deg"], sig["ord_y"], sig["mult_y_plus_1"], sig["cofactor_deg"])
   == (377, 107, 54, 216))
g = d56["galois_transfer"]
ok("D: (56,84) forcing residual is the 10th cyclotomic",
   sp.expand(sp.sympify(g["forcing_candidate"]["poly"])
             - (y**4 - y**3 + y**2 - y + 1)) == 0)
ok("D: (56,84) label C4, disc class 5, kills transfer, CONDITIONAL",
   g["galois"]["label"] == "C4" and g["galois"]["disc_class"] == 5
   and g["verdicts"]["C08"] == "KILLS" and g["verdicts"]["C20"] == "KILLS"
   and g["status"].startswith("CONDITIONAL"))
ok("D: (56,84) N-formula judgment: slice picture transfers ((b-1)/a=1)",
   any("is integral" in j for j in d56["judgment"]))

# ---------------------------------------------------------------------------
# E. Galois routine witnesses (independent of the dossiers)
# ---------------------------------------------------------------------------
ok("E: square_class(disc of audited quartic) = 17",
   cc.square_class(sp.discriminant(
       2048*y**4 - 512*y**3 + 320*y**2 - 240*y + 195, y)) == 17)
ok("E: squarefree(105*17) = 1785 != 1 (C08 kill arithmetic)",
   cc.square_class(105 * 17) == 1785)
ok("E: squarefree(170*17) = 10 != 1 (C20 kill arithmetic, GALOIS_LIBRARY)",
   cc.square_class(170 * 17) == 10)
ok("E: squarefree(459680) = 170 (the qs sharpness witness disc, 170*52^2)",
   cc.square_class(459680) == 170)
lab = cc.galois_label(y**4 - y**3 + y**2 - y + 1)
ok("E: Phi_10 labeled C4 with disc class 5 by the generic router",
   lab["label"] == "C4" and lab["disc_class"] == 5)
# constructive cross-check that Gal(Phi_10) is cyclic of order 4: the map
# sigma: x -> x^3 permutes the roots (Phi10(x^3) = 0 mod Phi10) and has order 4.
P10 = sp.Poly(y**4 - y**3 + y**2 - y + 1, y)
def power_mod(k):
    return sp.rem(sp.Poly(y**k, y), P10)
ok("E: sigma(x)=x^3 is a root map of Phi_10 (Phi10(x^3) = 0 mod Phi10)",
   sp.rem(sp.Poly((y**3)**4 - (y**3)**3 + (y**3)**2 - y**3 + 1, y), P10)
   == sp.Poly(0, y))
ok("E: sigma has order 4 (x^(3^2), x^(3^3) != x; x^(3^4) = x mod Phi_10)",
   power_mod(9) != sp.rem(sp.Poly(y, y), P10)
   and power_mod(27) != sp.rem(sp.Poly(y, y), P10)
   and power_mod(81) == sp.rem(sp.Poly(y, y), P10))
lab = cc.galois_label(y**4 + 1)
ok("E: y^4+1 labeled V4 -> verdict WITNESS-NEEDED",
   lab["label"] == "V4"
   and cc.transfer_verdicts("V4", lab["disc_class"])["C08"].startswith("WITNESS-NEEDED"))
ok("E: a disc-class-105 field vetoes C08 but not C20 (rule sharpness)",
   cc.transfer_verdicts("C2", 105)["C08"].startswith("OBSTRUCTION-VANISHES")
   and cc.transfer_verdicts("C2", 105)["C20"] == "KILLS")
ok("E: A4 always kills (no quadratic subfield)",
   cc.transfer_verdicts("A4", 1)["C08"] == "KILLS")

# ---------------------------------------------------------------------------
# F. survey-wide sanity
# ---------------------------------------------------------------------------
for row in cc.FAMILIES_LEN1:
    name, A0, A0p, p, l, q, k, (m0, dm), (n0, dn) = row
    j = 0
    while gcd(m0 + dm * j, n0 + dn * j) != 1:
        j += 1
    m, n = m0 + dm * j, n0 + dn * j
    ok("F: %s Diophantine identity at smallest coprime j=%d" % (name, j),
       (m + n) * q * k - n * (q * l - p) == k)
d7 = cc.compile_case("F7", 0)
s7 = d7["phi_prediction"]["signature"]
ok("F: F7 j=0 (42,147) grounded via the RAMIFIED law (PHI_F7.md) and "
   "reproduces the derived signature (250,165,83,2), branch-annotated",
   d7["phi_prediction"]["conjectural"] is False
   and (s7["deg"], s7["ord_y"], s7["mult_y_plus_1"], s7["cofactor_deg"])
       == (250, 165, 83, 2)
   and "ramified" in s7["branch"])
d12 = cc.compile_case("F12", 0)
ok("F: F12 flagged CONJECTURAL (A0' != (1,0): chart unverified)",
   d12["phi_prediction"]["conjectural"] is True
   and any("A0'" in r for r in d12["phi_prediction"]["conjectural_reasons"]))
d4 = cc.compile_case("F4", 0)
ok("F: F4 flagged CONJECTURAL (k=2: N-formula unverified)",
   d4["phi_prediction"]["conjectural"] is True
   and any("k = 2" in r for r in d4["phi_prediction"]["conjectural_reasons"]))
d108 = cc.compile_case("GGHV_108_144")
sig = d108["phi_prediction"]["signature"]
ok("F: (108,144) special entry gives (550,205,69,276) and Phi_10 forcing",
   (sig["deg"], sig["ord_y"], sig["mult_y_plus_1"], sig["cofactor_deg"])
   == (550, 205, 69, 276)
   and sp.expand(sp.sympify(d108["galois_transfer"]["forcing_candidate"]["poly"])
                 - (y**4 - y**3 + y**2 - y + 1)) == 0)

# ---------------------------------------------------------------------------
# G. canonical determinism + on-disk pilot dossiers
# ---------------------------------------------------------------------------
for name, j, fname in [("GGHV_72_108", None, "case_dossier_GGHV_72_108.json"),
                       ("F2", 1, "case_dossier_F2_j1_75_125.json"),
                       ("F9", 0, "case_dossier_F9_j0_56_84.json")]:
    one = cc.canonical_json(cc.compile_case(name, j))
    two = cc.canonical_json(cc.compile_case(name, j))
    ok("G: %s canonical JSON deterministic" % (name if j is None else "%s j=%d" % (name, j)),
       one == two)
    with open(fname, encoding="utf-8") as fh:
        ok("G: on-disk %s matches recomputation" % fname, fh.read() == one)
    ok("G: %s parses and carries schema case-dossier-v1" % fname,
       json.loads(one)["schema"] == "case-dossier-v1")

print("\nALL %d CASE-COMPILER CHECKS PASSED" % count)
