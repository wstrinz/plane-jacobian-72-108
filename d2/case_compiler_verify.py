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
# A. the law reproduces every landed point (constants from the PHI docs)
#
# *** 2026-07-27: THIS TABLE WAS THE DUPLICATED-GROUND-TRUTH BUG. ***
#
# Until today the rows below read
#     ("(75,125)", 3, 5, 5, 5, 2, (504, 201, 101, 202))
#     ("(56,84)",  2, 3, 7, 7, 2, (377, 107, 54, 216))
#     ("(50,75)",  2, 3, 5, 5, 2, (189, 75,  38,  76))
# while cc.KNOWN_POINTS already held the REPAIRED (80,80,0,0) and (30,30,0,0).
# Both halves were green: check A fed the STALE chart data into the law and got
# the STALE signature back (self-consistent), and check C compiled the dossier,
# which routes through the guard, and got the REPAIRED one (also self-consistent).
# Two independent copies of the same ground truth, one repaired and one not,
# neither aware of the other -- the exact pattern family_grammar_verify.py's
# check A9 was written for after the F3 incident.
#
# So this table is repaired AND cross-checked against cc.KNOWN_POINTS below
# (check A4).  Independence of the two copies is only worth something if they are
# compared; without A4, repairing one of them again next time will go unnoticed.
#
# The (a,b,t,degC,ordC) here are the GUARDED chart data.  A3 re-derives them from
# polygon_reduction so this file does not merely transcribe the compiler's.
# ---------------------------------------------------------------------------
LANDED = [
    # (label, a, b, t, degC, ordC, expected signature, dossier tag)
    ("(72,108)",  2, 3, 4, 8, 7, (238, 204, 30, 4),    "GGHV_72_108"),
    ("(108,144)", 3, 4, 4, 8, 3, (550, 205, 69, 276),  "GGHV_108_144"),
    ("(75,125)",  3, 5, 4, 1, 1, (80, 80, 0, 0),       "F2_j1_75_125"),
    ("(56,84)",   2, 3, 3, 1, 1, (22, 22, 0, 0),       "F9_j0_56_84"),
    ("(50,75)",   2, 3, 4, 1, 1, (30, 30, 0, 0),       "F2_j0_50_75"),
    ("(66,231)",  2, 7, 3, 9, 4, (375, 165, 42, 168),  "F14_j0_66_231"),
    ("(48,64)",   3, 4, 3, 1, 1, (51, 51, 0, 0),       "F1_j0_48_64"),
    ("(66,99)",   2, 3, 3, 9, 8, (195, 169, 22, 4),    "F17_j0_66_99"),
    ("(63,147)",  3, 7, 3, 6, 5, (448, 371, 75, 2),    "F8_j0_63_147"),
    ("(99,231)",  3, 7, 3, 9, 5, (672, 371, 297, 4),   "F15_j0_99_231"),
]
for label, a, b, t, degC, ordC, want, _tag in LANDED:
    s = cc.corner_signature(a, b, t, degC, ordC)
    L = cc.law_signature(s)
    got = (L["deg"], L["ord_y"], L["mult_y_plus_1"], L["cofactor_deg"])
    ok("A: law reproduces %s -> %s" % (label, want), got == want)

# A2. THE INDEPENDENT TARGET.  ord_y(Phi) = a*q*M - H, PROVED in
# BRIDGE_GENERALITY.md, computed here from (a,b,t,kappa,ord C) alone.  Before the
# repair the law was checked only against targets the same chart dictionary had
# produced, so agreement carried no information.
def bridge(a, b, t, kappa, ordC):
    s_ = a + b
    return a * ordC * (t * s_ - (kappa + 1)) - (ordC * s_ - 1)

ok("A2: ord_y == a*q*M - H (PROVED bridge identity) at every landed point",
   all(cc.law_signature(cc.corner_signature(a, b, t, dC, oC))["ord_y"]
       == bridge(a, b, t, t - 2, oC)
       for _l, a, b, t, dC, oC, _w, _tg in LANDED))

# A3. the chart data in this table must come from the GUARD, not from GGV5's
# final chain corner.  Re-derived here from the corner alone.
import polygon_reduction as _pr                                  # noqa: E402
_CORNER_OF = {"(72,108)": ((8, 28), 4, 7), "(108,144)": ((8, 28), 4, 3),
              "(75,125)": ((5, 20), 5, 2), "(56,84)": ((7, 21), 7, 2),
              "(50,75)": ((5, 20), 5, 2),  "(66,231)": ((9, 24), 3, 4),
              "(48,64)": ((4, 12), 4, 3),  "(66,99)": ((9, 24), 3, 8),
              "(63,147)": ((6, 15), 3, 5), "(99,231)": ((9, 24), 3, 5)}
_bad3 = []
for label, a, b, t, degC, ordC, want, _tag in LANDED:
    A0, lf, bf = _CORNER_OF[label]
    cd = _pr.corner_chart_data(A0[0], A0[1], l_final=lf, b_final=bf,
                               who="case_compiler_verify " + label)
    if (cd["t"], cd["deg_C"], cd["ord_C"]) != (t, degC, ordC):
        _bad3.append((label, (cd["t"], cd["deg_C"], cd["ord_C"]), (t, degC, ordC)))
ok("A3: every landed row's (t, deg C, ord C) is what polygon_reduction's guard "
   "returns for its corner -- not GGV5's (l_final, a0, b_final).  Violations: %s"
   % (_bad3 or "none"), not _bad3)
_refused = [lbl for lbl in _CORNER_OF
            if not _pr.has_retraction(*_CORNER_OF[lbl][0])]
ok("A3b: and the refused corners among them are exactly (75,125), (56,84), "
   "(50,75), (48,64) -- (8,28), (6,15) and (9,24) all RETRACT, so the audited "
   "home case, the corner-144 template and F14/F17/F8/F15 are untouched",
   set(_refused) == {"(75,125)", "(56,84)", "(50,75)", "(48,64)"})

# A4. THE CROSS-CHECK.  This file's table and cc.KNOWN_POINTS are independent
# copies; they must AGREE.  Its absence is what let the two halves drift.
_theirs = {tag: cc.KNOWN_POINTS[tag][0] for tag in cc.KNOWN_POINTS}
_mine = {tag: want for _l, _a, _b, _t, _dC, _oC, want, tag in LANDED}
ok("A4: DRIFT GUARD -- this file's landed table and case_compiler.KNOWN_POINTS "
   "agree on every key.  (Added 2026-07-27: they did NOT, at (75,125), (56,84) "
   "and (50,75); each copy was self-consistent and both checkers were green.)",
   _theirs == _mine)
ok("A4b: and every monomial-corner leading constant registered in "
   "case_compiler.KNOWN_LC is exactly 1/a, which C = y forces",
   all(cc.KNOWN_LC[tag] == "1/%d" % a
       for _l, a, _b, _t, dC, _oC, _w, tag in LANDED
       if dC == 1 and tag in cc.KNOWN_LC))

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
# REPAIRED 2026-07-26 (PASSPORT_75_125_REPAIR.md): t 5->4, a0 (=deg C) 5->1,
# q (=ord C) 2->1, dg 3->0 (C = y is a MONOMIAL), N 98->77.  gap = (q-1)-a0/t is
# now -1/4, non-integral, so gap_effective = 0 and no extra unit factor appears
# -- confirmed independently by the ODE solve in phi_75_125_verify.py sec.E.
ok("C: (75,125) corner signature  [REPAIRED]",
   (s["a"], s["b"], s["t"], s["a0"], s["q"], s["e"], s["r"],
    s["dg"], s["N"]) == (3, 5, 4, 1, 1, 3, -1, 0, 77))
ok("C: (75,125) gap = -1/4 is NON-integral, so gap_effective = 0 (no resonance "
   "at or above the ansatz degree)",
   (s["gap_num"], s["gap_den"]) == (-1, 4) and s["gap"] is None)
sig = d75["phi_prediction"]["signature"]
ok("C: (75,125) law signature (80,80,0,0)  [REPAIRED]",
   (sig["deg"], sig["ord_y"], sig["mult_y_plus_1"], sig["cofactor_deg"])
   == (80, 80, 0, 0))
ok("C: (75,125) dossier records the REPAIRED judgment and cites GGV3's three "
   "published integers",
   any("REPAIRED 2026-07-2" in j and "deg(Q_1) = 15" in j
       for j in d75["judgment"]))
ok("C: (75,125) regime is monomial_corner, NOT resonance_gap_r0  [LABEL-INTEGRITY "
   "FIX 2026-07-27: the old regime key sent every monomial corner to "
   "'resonance_gap_r0', because gap = -1/t is non-integral (so the key's first "
   "component read True) and r = -1 (so the second read False).  A case with no "
   "resonance and no residual was labelled a resonance case]",
   d75["phi_prediction"]["regime"] == "monomial_corner")
ok("C: (75,125) dossier carries the bridge check and it AGREES with the law",
   d75["phi_prediction"]["bridge_check"]["ord_y_bridge"] == 80
   and "PROVED" in d75["phi_prediction"]["bridge_check"]["status"])
ok("C: (75,125) grounded regime, not conjectural, Diophantine OK",
   d75["phi_prediction"]["conjectural"] is False
   and d75["case"]["diophantine"].endswith("OK"))
g = d75["galois_transfer"]
# REPAIRED 2026-07-26.  The residual H2 = y^2-y+1 (label C2, disc class -3) does
# NOT exist: it was g/(y+1) with g = y^3+1, and g came from deg C = a0 = 5, which
# came from the retracted shape (5,20) does not have.  With C = y a monomial there
# is no residual polynomial, so the whole forcing-polynomial / Galois-descent
# layer is VACUOUS at this corner -- an ABSENT object, not an unresolved one.
ok("C: (75,125) has NO forcing residual: the Galois-transfer layer is VACUOUS "
   "(was H2 = y^2-y+1, label C2, disc class -3, all of which came from deg C=5)",
   g["forcing_candidate"]["poly"] is None
   and g["forcing_candidate"]["rationale"].startswith("VACUOUS")
   and g["galois"] is None and g["status"] == "VACUOUS")
ok("C: (75,125) kill classes report VACUOUS, not KILLS and not UNKNOWN -- so the "
   "superseded 'kills transfer (C08,C20)' claim is WITHDRAWN, and it is withdrawn "
   "for a structural reason rather than left pending",
   g["verdicts"]["C08"].startswith("VACUOUS")
   and g["verdicts"]["C20"].startswith("VACUOUS"))
# and the contrast case still works, so VACUOUS is discriminating.
# 2026-07-27: the contrast MOVED from F9 (56,84) to F14 (66,231).  F9's corner
# (7,21) is guard-refused, so it too is monomial now and its residual is VACUOUS
# as well -- using it as the contrast would have made this check vacuously true
# while still passing.  F14's corner (9,24) RETRACTS, deg C = 9, ord C = 4, so it
# genuinely has the residual g = y^5+1 and a real Galois layer.
_g231 = cc.compile_case("F14", 0)["galois_transfer"]
ok("C: contrast -- F14 (66,231), at the RETRACTING corner (9,24), still HAS a "
   "residual and still reports KILLS, so the VACUOUS verdict is discriminating, "
   "not blanket",
   _g231["forcing_candidate"]["poly"] is not None
   and _g231["verdicts"]["C08"] == "KILLS")
_g56 = cc.compile_case("F9", 0)["galois_transfer"]
ok("C: and F9 (56,84) is now VACUOUS TOO -- corner (7,21) is refused, so its "
   "residual and its C08/C20 'kills transfer' claim are WITHDRAWN for the same "
   "structural reason as (75,125)'s.  (This is why the contrast case had to "
   "move: with F9 as contrast the check above would pass vacuously.)",
   _g56["forcing_candidate"]["poly"] is None
   and _g56["verdicts"]["C08"].startswith("VACUOUS"))
ok("C: (75,125) presentations schematic (no f31-analogue exists)",
   d75["presentations"]["eliminated_f31_style"]["instantiated"] is False
   and d75["presentations"]["pre_resultant_G_system"]["instantiated"] is False
   and d75["presentations"]["eliminated_f31_style"]["parameters"] is None)
ok("C: (75,125) N-formula judgment records the broken slice picture ((b-1)/a=4/3)",
   any("is NOT integral" in j for j in d75["judgment"]))

# ---------------------------------------------------------------------------
# D. F9 j=0 (56,84)  --  REPAIRED 2026-07-27.
#
# This dossier used to read (t, deg C, ord C) = (7,7,2) with dg=5, N=52, signature
# (377,107,54,216), a 10th-cyclotomic forcing residual, Galois label C4, disc
# class 5 and "kills transfer (C08, C20), CONDITIONAL".  ALL OF THAT IS WITHDRAWN.
# Corner (7,21) does not satisfy the retraction shape (ceil(21/7) = 3 and
# 3*(7-1) = 18 != 21), and it is refuted IN PRINT: GGHV22 2204.14178.tex:1394
# publishes the chart phi_3(y) = y x^3 with [P,Q] = x there -- l = 3, kappa = 1.
# So the chart is (3,1,1,1): C = y a monomial, no residual, no Galois layer.
# The repaired ord_y(Phi) = 22 agrees with the PROVED bridge identity, and
# bridge_generality.py MUT F independently displaces 22 <- 107.
# ---------------------------------------------------------------------------
d56 = cc.compile_case("F9", 0)
ok("D: F9 j=0 is the (56,84) case", d56["case"]["degrees"] == [56, 84])
s = d56["corner_signature"]
ok("D: (56,84) corner signature  [REPAIRED 2026-07-27: t 7->3, deg C 7->1, "
   "ord C 2->1, dg 5->0, N 52->20]",
   (s["a"], s["b"], s["t"], s["a0"], s["q"], s["e"], s["r"],
    s["dg"], s["N"]) == (2, 3, 3, 1, 1, 2, -1, 0, 20))
ok("D: (56,84) gap = -1/3 is NON-integral and NEGATIVE, so gap_effective = 0",
   (s["gap_num"], s["gap_den"]) == (-1, 3) and s["gap"] is None)
sig = d56["phi_prediction"]["signature"]
ok("D: (56,84) law signature (22,22,0,0)  [was (377,107,54,216)]",
   (sig["deg"], sig["ord_y"], sig["mult_y_plus_1"], sig["cofactor_deg"])
   == (22, 22, 0, 0))
ok("D: (56,84) bridge check present and agreeing: ord_y = 22 = a*q*M - H "
   "(bridge_generality.py MUT F displaces this by exactly 22 <- 107)",
   d56["phi_prediction"]["bridge_check"]["ord_y_bridge"] == 22)
ok("D: (56,84) regime is monomial_corner",
   d56["phi_prediction"]["regime"] == "monomial_corner")
ok("D: (56,84) dossier records the REPAIRED judgment and cites GGHV22's PUBLISHED "
   "chart at (7,21) -- this row is refuted in print, not merely unproved",
   any("REPAIRED 2026-07-27" in j and "2204.14178.tex:1394" in j
       for j in d56["judgment"]))
g = d56["galois_transfer"]
ok("D: (56,84) forcing residual WITHDRAWN: the 10th cyclotomic came from "
   "g = y^5+1, which came from deg C = a0 = 7, which came from the retracted "
   "shape (7,21) does not have.  With C = y there is no residual at all",
   g["forcing_candidate"]["poly"] is None and g["galois"] is None)
ok("D: (56,84) 'label C4, disc class 5, kills transfer, CONDITIONAL' is WITHDRAWN "
   "and replaced by VACUOUS -- an ABSENT object, not an unresolved one",
   g["status"] == "VACUOUS" and g["verdicts"]["C08"].startswith("VACUOUS")
   and g["verdicts"]["C20"].startswith("VACUOUS"))
ok("D: (56,84) N-formula judgment: slice picture transfers ((b-1)/a=1)",
   any("is integral" in j for j in d56["judgment"]))
# and the 10th cyclotomic is still a real object of the program -- at F14, whose
# corner retracts and whose dg is also 5.  So what moved is the ATTRIBUTION.
_g231 = cc.compile_case("F14", 0)["galois_transfer"]
ok("D: the 10th cyclotomic residual is NOT lost to the program: it is F14 "
   "(66,231)'s residual, at the RETRACTING corner (9,24) where dg = 5 as well.  "
   "What the repair moved is the attribution, not the object",
   sp.expand(sp.sympify(_g231["forcing_candidate"]["poly"])
             - (y**4 - y**3 + y**2 - y + 1)) == 0
   and _g231["galois"]["label"] == "C4" and _g231["galois"]["disc_class"] == 5)

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
