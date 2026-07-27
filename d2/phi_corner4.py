#!/usr/bin/env python3
"""phi_corner4.py  (NEW; read-only over all existing artifacts)

FOURTH-CORNER Phi derivation: break (or explain) the t/kappa correlation in the
corner-signature fit of PHI_75_125.md.

*** 2026-07-27 CHART REPAIR -- READ THIS BEFORE ANY NUMBER BELOW. ***

Until today this file read its chart data (t, deg C, ord C) off GGV5's FINAL
CHAIN CORNER A_1 = (p\\l_final, b_final) by the dictionary
(t, deg C, ord C) = (l_final, a0, b_final).  Commit 2adb92a established that the
dictionary is valid ONLY on the RETRACTION SHAPE  b0 == l_chart*(a0 - 1)  with
l_chart = ceil(b0/a0), and polygon_reduction.final_corner_dictionary() now
RAISES off it.  Of the seven distinct corners in GGV5's v11<=35 length-1 tables
only TWO retract -- (6,15) and (9,24).  The other five -- (4,12), (5,20),
(7,21), (8,24), (9,21) -- do NOT, and there the corner has NO vertical top face,
so C is a MONOMIAL: deg C = ord C = 1, C = y.  Eleven of the seventeen family
rows are affected (F1-F6, F9-F13).  All chart data below now comes from
polygon_reduction.corner_chart_data(), which raises rather than guessing.

WHAT THAT COSTS THIS FILE'S ORIGINAL HEADLINE.  The out-of-sample test was
"F9's corner (7,21) gives t = 7".  It does not: chart_exponent(7,21) = 3, and
GGHV22 PUBLISHES the chart phi_3(y) = y x^3 with [P,Q] = x at exactly this
corner (2204.14178.tex:1394), i.e. l = 3, kappa = 1.  l_final = 7 was never the
chart exponent.  The claim "t can be varied" is therefore much weaker than
stated: after the repair EVERY corner in GGV5's v11<=35 length-1 tables has
t in {3,4} (STEP 1b proves this), so the tables contain no t-variation beyond
the two values already known from (8,28).

WHAT SURVIVES, UNCHANGED.  kappa = t-2 is still FORCED for every corner reached
by the standard single Laurent chart (X,Y) -> (x^-1, x^l y): its Jacobian is
-x^(l-2) and t = l.  That is a statement about the CHART, not about which l the
chart uses, so the repair does not touch it (STEP 2).  The fit still loses a
parameter, N = a[t(a+b) - (kappa+1)] - 2b = a[t(a+b-1) + 1] - 2b.

INDEPENDENT TARGET.  Every repaired row is checked against the bridge identity
        ord_y(Phi) = a*q*M - H,    M = t(a+b) - (kappa+1),  H = q(a+b) - 1,
PROVED in bridge_generality.py (59/59; rho = q(b-a)+1 by a local triangular
recursion at y=0, N = a*M-2b by the built D-transform tower).  Nothing here is
fitted to itself any more: STEP 3 derives Phi from the ODE and STEP 3b compares
it to a number this file does not produce.

Sources: GGV5 (paper_src/1708.07936_GGV5.tex) family tables (the two tables in
the final section, "Admissible complete chains with v11(A_0) <= 35"); the
corner-144 template CORNER_144_COMPARISON.md / corner144_verify.py; the
three-point fit PHI_75_125.md; the retraction guard polygon_reduction.py sec.0b
and PASSPORT_75_125_REPAIR.md; the bridge identity BRIDGE_GENERALITY.md.  The
independent PASS/FAIL checker is phi_corner4_verify.py.  Exact sympy throughout.
"""
import sympy as sp
from fractions import Fraction
from math import gcd

import polygon_reduction as pr

y = sp.symbols("y")

# ---------------------------------------------------------------------------
# 0. The GGV5 v11<=35 family tables, transcribed verbatim.
#    Length-1 families: (name, A0, A0', p, l, q, k, (m0, dm), (n0, dn))
#    with (m,n)(j) = (m0 + dm*j, n0 + dn*j) and final corner A1 = (p\l, q).
# ---------------------------------------------------------------------------
FAMILIES_LEN1 = [
    ("F1",  (4, 12), (1, 0),  7, 4, 3, 1, (3, 2),  (4, 3)),
    ("F2",  (5, 20), (1, 0),  7, 5, 2, 1, (2, 1),  (3, 2)),
    ("F3",  (5, 20), (1, 0),  8, 5, 3, 1, (3, 4),  (2, 3)),
    ("F4",  (5, 20), (1, 0),  8, 5, 3, 2, (3, 2),  (16, 12)),
    ("F5",  (5, 20), (1, 0),  9, 5, 4, 1, (9, 7),  (5, 4)),
    ("F6",  (5, 20), (1, 0),  9, 5, 4, 2, (7, 6),  (18, 16)),  # CORRECTED 2026-07-24: GGV5 prints F6 base (m,n)=(4,10) [gcd=2, violates coprimality]; the coprime family is (6j+7,16j+18)=base (7,18). See CHAIN_SURVEY.md. (Survey advances j to first coprime pair -> computed (m,n)=(7,18) unchanged.)
    ("F7",  (6, 15), (1, 0),  7, 3, 4, 1, (2, 1),  (7, 4)),
    ("F8",  (6, 15), (1, 0),  8, 3, 5, 1, (3, 2),  (7, 5)),
    ("F9",  (7, 21), (1, 0), 11, 7, 2, 1, (2, 1),  (3, 2)),
    ("F10", (7, 21), (1, 0), 13, 7, 3, 1, (7, 5),  (4, 3)),
    ("F11", (7, 21), (1, 0), 13, 7, 3, 2, (2, 1),  (5, 3)),
    ("F12", (8, 24), (2, 0), 13, 4, 5, 1, (3, 2),  (7, 5)),
    ("F13", (9, 21), (2, 0), 13, 3, 7, 1, (2, 1),  (13, 7)),
    ("F14", (9, 24), (1, 0),  7, 3, 4, 1, (2, 1),  (7, 4)),
    ("F15", (9, 24), (1, 0),  8, 3, 5, 1, (3, 2),  (7, 5)),
    ("F16", (9, 24), (1, 0), 10, 3, 7, 1, (3, 4),  (5, 7)),
    ("F17", (9, 24), (1, 0), 11, 3, 8, 1, (2, 5),  (3, 8)),
]
# Length-2 families (final corner only; composite chart NOT derived in any
# paper -- t and kappa are UNVERIFIED for these, see PHI_CORNER4.md):
FAMILIES_LEN2 = [
    ("F18", (6, 18),  7, 3, 4), ("F19", (6, 18),  8, 3, 5),
    ("F20", (6, 24),  7, 3, 4), ("F21", (6, 24),  8, 3, 5),
    ("F22", (8, 24),  5, 4, 2), ("F23", (8, 24), 11, 4, 4),
    ("F24", (8, 24), 19, 8, 3),
]

# ---------------------------------------------------------------------------
# 0b.  THE RETRACTION LEDGER  (2026-07-27 chart repair).
#
# GGV5's per-row final chain corner A_1 = (p\l_final, b_final) carries this
# corner's CHART data only on the retraction shape
#
#       b0 == l_chart * (a0 - 1),        l_chart = chart_exponent(a0,b0)
#                                                = ceil(b0/a0),
#
# i.e. exactly when the flipped edge {(0,1),(b0,a0)} collapses to a VERTICAL top
# face under the inversion (a,b) -> (l*b - a, b).  Off that shape there is no
# vertical top face, hence deg C = 1 and C is the MONOMIAL y.
#
# THE LEDGER IS NOT THE SOURCE OF THE NUMBERS -- polygon_reduction's guard is.
# The ledger records, per corner, WHICH answer the guard must give and WHY, and
# is asserted against the guard in chart_of_corner() below.  So a corner whose
# classification drifts, or a new corner added to the tables, RAISES instead of
# silently producing a wrong signature.  This is the whole point of the repair:
# the failure mode being fixed is a module that fitted its own superseded model
# to its own superseded targets and reported MATCHES.
#
#   corner   l_final  l_chart   retraction test          verdict
#   (4,12)      4        3      12 != 3*(4-1) = 9        REFUSED -> C = y
#   (5,20)      5        4      20 != 4*(5-1) = 16       REFUSED -> C = y
#   (6,15)      3        3      15 == 3*(6-1) = 15       RETRACTS -> deg C = 6
#   (7,21)      7        3      21 != 3*(7-1) = 18       REFUSED -> C = y
#   (8,24)      4        3      24 != 3*(8-1) = 21       REFUSED -> C = y
#   (9,21)      3        3      21 != 3*(9-1) = 24       REFUSED -> C = y
#   (9,24)      3        3      24 == 3*(9-1) = 24       RETRACTS -> deg C = 9
#
# TWO OF THE FIVE REFUSALS ARE REFUTED IN PRINT, not merely unproved:
#   (7,21) -- GGHV22 2204.14178.tex:1394 publishes phi_3(y) = y x^3 and
#             [P,Q] = x for this corner: l = 3, kappa = 1, against l_final = 7.
#   (5,20) -- GGV3 1406.0886.tex:1723-1727 publishes [P_1,Q_1] = x^2,
#             deg P_1 = 10, deg Q_1 = 15 for the sibling (50,75) AT THIS CORNER:
#             three integers, all reproduced by l = 4 and all contradicted by 5.
# The chart_exponent rule itself is [INFERRED] (polygon_reduction sec.0b): it is
# validated on all five published GGHV22 reductions and pinned at (5,20), but it
# appears in no published proposition.  Do not cite it as published.
# ---------------------------------------------------------------------------
RETRACTION_LEDGER = {
    (4, 12): (False, "l_chart = ceil(12/4) = 3;  12 != 3*(4-1) = 9"),
    (5, 20): (False, "l_chart = ceil(20/5) = 4;  20 != 4*(5-1) = 16  "
                     "[GGV3 pins l=4 via (50,75)]"),
    (6, 15): (True,  "l_chart = ceil(15/6) = 3;  15 == 3*(6-1)  -> vertical face"),
    (7, 21): (False, "l_chart = ceil(21/7) = 3;  21 != 3*(7-1) = 18  "
                     "[GGHV22 PUBLISHES l=3, [P,Q]=x at this corner]"),
    (8, 24): (False, "l_chart = ceil(24/8) = 3;  24 != 3*(8-1) = 21"),
    (9, 21): (False, "l_chart = ceil(21/9) = 3;  21 != 3*(9-1) = 24"),
    (9, 24): (True,  "l_chart = ceil(24/9) = 3;  24 == 3*(9-1)  -> vertical face"),
}

# Rows whose corner is guard-REFUSED, hence whose pre-2026-07-27 (t, deg C,
# ord C) were the broken dictionary's output.  DERIVED below from the ledger and
# asserted, never trusted as a literal.
AFFECTED_EXPECTED = {"F1", "F2", "F3", "F4", "F5", "F6",
                     "F9", "F10", "F11", "F12", "F13"}

# The A0' = (2,0) rows.  Their corner is refused too, but chart_exponent is
# validated only against A0' = (1,0) published reductions, so their repaired
# chart data carries ONE MORE conditional than the other nine.  F13's corner
# (9,21) is additionally the only refused corner with NO corroborating row in
# corner_atlas.json.  See the TRIPWIRE at the end of STEP 1.
A0P_20_ROWS = {"F12", "F13"}


def chart_of_corner(A0, l_final, b_final, who=""):
    """(t, kappa, deg C, ord C, retracts) for one corner, THROUGH THE GUARD.

    polygon_reduction.corner_chart_data does the work: on the retraction shape it
    takes ord C = b_final through final_corner_dictionary (which raises off the
    shape), otherwise it returns the monomial data (deg C = ord C = 1).  The
    RETRACTION_LEDGER entry is asserted against the guard, so drift RAISES.
    """
    a0, b0 = int(A0[0]), int(A0[1])
    retr = pr.has_retraction(a0, b0)
    assert (a0, b0) in RETRACTION_LEDGER, \
        ("corner %s is not in RETRACTION_LEDGER: classify it (and say why) "
         "before deriving anything from it" % ((a0, b0),))
    want, why = RETRACTION_LEDGER[(a0, b0)]
    assert retr is want, ("RETRACTION_LEDGER disagrees with the guard at %s: "
                          "ledger says %s (%s), guard says %s"
                          % ((a0, b0), want, why, retr))
    cd = pr.corner_chart_data(a0, b0, l_final=l_final, b_final=b_final, who=who)
    assert cd["retraction"] is retr and cd["kappa"] == cd["t"] - 2, cd
    return cd["t"], cd["kappa"], cd["deg_C"], cd["ord_C"], retr


# ---------------------------------------------------------------------------
# 0c.  SUPERSEDED (pre-repair) chart data and signatures, kept LABELLED.
#
# Deleting these would make the repair unfalsifiable: phi_corner4_verify.py's
# mutation control (sec. H3) reinstates them and must then FAIL.  Same discipline
# as family_grammar.SUPERSEDED_F.
#
# The recorded quantities are the MODULE-INDEPENDENT ones -- the chart data the
# refused dictionary returned, N, and ord_y(Phi).  The FULL signature is NOT
# recorded here because the three PHI_* files disagreed about its mult/cofactor
# components (phi_corner4's fit, phi_f14's unified law and phi_f7's amended
# ramified law give different mult at the same row); ord_y is the one component
# every route agreed on, and it is exactly the component the bridge identity
# pins.  Per-row provenance of the published full signature:
#
#   F1  (48,64)   phi_f14.py  law_sig -> (275, 205,   69,    1)
#   F2  (50,75)   phi_corner4 fit     -> (189,  75,   38,   76)
#   F3  (75,50)   phi_f7.py   ram law -> (189, 112,   75,    2)
#   F9  (56,84)   phi_corner4 fit     -> (377, 107,   54,  216)
#   F10 (196,112) phi_f7.py   ram law -> (1917, 820, 1093,   4)
#
#         row:  (t_stale, degC_stale, ordC_stale, N_stale, ordy_stale)
# ---------------------------------------------------------------------------
SUPERSEDED = {
    "F1":  (4, 4, 3,  67, 205),        # (48,64)
    "F2":  (5, 5, 2,  36,  75),        # (50,75)   [also pre-2026-07-26]
    "F3":  (5, 5, 3,  36, 112),        # (75,50)
    "F9":  (7, 7, 2,  52, 107),        # (56,84)
    "F10": (7, 7, 3, 270, 820),        # (196,112)
}


# ---------------------------------------------------------------------------
# 2026-07-26 corner-law generalizations (PASSPORT_75_125_REPAIR.md).  See the
# patch note in the repair doc; both are forced by the (5,20) corner, which is
# the first case with C a MONOMIAL (deg g = 0) and with a non-integral gap.
# ---------------------------------------------------------------------------
def gap_effective(gap):
    """The resonance gap AS AN EXTRA UNIT-FACTOR DEGREE.

    gap = (q-1) - a0/t is the resonant degree minus the pure-ansatz degree of f.
    It contributes an extra factor only when it is a POSITIVE INTEGER ((72,108):
    gap = 4).  Negative or non-integral gap means the resonance does not sit at
    or above the ansatz degree, no extra factor appears, and the pure ansatz is
    exact -- effective gap 0.  At the repaired (5,20) corner gap = -1/4 and the
    independent ODE solve confirms deg f = rho exactly.
    """
    from fractions import Fraction as _F
    g = _F(gap)
    return int(g) if (g.denominator == 1 and g > 0) else 0


def mult_and_cofactor(e, N, degC, ordC, gap):
    """(mult_(y+1), cofactor_deg) with the residual-free branch dg = degC-ordC = 0.

    2026-07-27: the third and fourth arguments are deg_y(C) and ord_y(C), NOT the
    corner coordinate a0 and the chain datum b_final.  Those coincide only under
    GGV5's final-corner dictionary, i.e. only on the retraction shape.

    dg > 0 : g = y^dg + 1 contributes (y+1)^(e+N) and a residual H2 = g/(y+1) of
             degree dg-1, so mult = e+N and cof = gap + (dg-1)*(e+N).
    dg == 0: C is a MONOMIAL.  There is NO g, hence no (y+1) place and no
             residual: mult = 0 and cof = gap.
    """
    ge = gap_effective(gap)
    dg = degC - ordC
    if dg == 0:
        return 0, ge
    return e + N, ge + (dg - 1) * (e + N)


def bridge_ord(a, b, t, kappa, ordC):
    """ord_y(Phi) = a*q*M - H  --  the INDEPENDENT target (BRIDGE_GENERALITY.md).

    PROVED there, not fitted here: rho := ord_y(f) = q(b-a)+1 by the local
    triangular recursion at y=0 (unique excluded locus t = q(kappa+1), whose only
    standard-class point is (t,kappa,q) = (2,0,2), off every published row), and
    N = a*M - 2b from the built D-transform tower.  q is ord_y(C).
    """
    s = a + b
    M = t * s - (kappa + 1)
    H = ordC * s - 1
    return a * ordC * M - H


def fit_signature(a, b, t, kappa, degC, ordC):
    """The PHI_75_125.md six-parameter closed forms.

    REPAIRED 2026-07-26: the old body was the r > 0 regime only and silently
    returned nonsense when C is a MONOMIAL (deg g = degC - ordC = 0): then
    r = -1, and mult = e+N / cof = -(e+N) are both wrong -- there is no (y+1)
    place and no residual cofactor at all.  Since the (5,20) corner IS in that
    regime (deg C = ord C = 1, C = y), the dg == 0 branch is required, not
    cosmetic.  Verified against the directly-derived (80,80,0,0) at (75,125) and
    (30,30,0,0) at (50,75) in phi_corner4_verify.py sec.G.

    REPAIRED 2026-07-27: the last two arguments are now (deg C, ord C) from
    polygon_reduction.corner_chart_data, never (a0, b_final) from GGV5's chain
    row.  The ord_y component is cross-checked against the PROVED bridge
    identity, so this is no longer a formula validated only against itself.
    """
    e = b - a + 1
    dg = degC - ordC             # degree of the residual g
    r = max(dg - 1, 0)           # degree of the residual H2 = g/(y+1); 0 if no g
    N = a * (t * (a + b) - (kappa + 1)) - 2 * b
    deg = (e * degC - ordC + 1) + N * degC
    ordy = ((e - 1) * ordC + 1) + N * ordC
    mult, cof = mult_and_cofactor(e, N, degC, ordC, 0)
    assert deg - ordy - mult == cof, (a, b, t, kappa, degC, ordC, deg, ordy, mult, cof)
    assert ordy == bridge_ord(a, b, t, kappa, ordC), \
        ("ord_y(Phi) disagrees with the PROVED bridge identity a*q*M - H",
         a, b, t, kappa, degC, ordC, ordy, bridge_ord(a, b, t, kappa, ordC))
    return N, e, r, (deg, ordy, mult, cof)

# ---------------------------------------------------------------------------
# 1. Candidate-corner survey (smallest j with gcd(m,n)=1 per family)
# ---------------------------------------------------------------------------
print("=" * 96)
print("STEP 1 -- candidate-corner survey (GGV5 v11<=35 tables, smallest coprime j)")
print("         CHART DATA THROUGH polygon_reduction.corner_chart_data -- NOT off")
print("         GGV5's final chain corner.  'chain' = the row's (l_final,b_final);")
print("         'chart' = (t, kappa, deg C, ord C) derived from the corner itself.")
print("=" * 96)
print(f"{'fam':4} {'A0':>8} {'A1':>10} {'(m,n)':>8} {'degs':>10} "
      f"{'t':>2} {'kap':>3} {'dC':>3} {'oC':>2} {'e':>2} {'r':>2} {'N':>4} "
      f"{'gap':>5} {'ordPhi':>7} {'bridge':>7}  notes")

survey = []
affected, bridge_bad = [], []
for name, A0, A0p, p, l, q, k, (m0, dm), (n0, dn) in FAMILIES_LEN1:
    j = 0
    while gcd(m0 + dm * j, n0 + dn * j) != 1:
        j += 1
    m, n = m0 + dm * j, n0 + dn * j
    v11 = A0[0] + A0[1]
    degs = (v11 * m, v11 * n)
    a, b = sorted((m, n))
    # ---- CHART data: from the corner, through the guard.  (l, q) stay CHAIN data.
    t, kappa, degC, ordC, retracts = chart_of_corner(
        A0, l_final=l, b_final=q, who="phi_corner4 %s" % name)
    if not (retracts and t == l and degC == A0[0] and ordC == q):
        affected.append(name)
    N, e, r, sig = fit_signature(a, b, t, kappa, degC, ordC)
    # resonance gap: resonant deg f  minus  pure-ansatz deg f.  Both are functions
    # of (deg C, ord C), NOT of (a0, b_final) -- that identification is precisely
    # the refused dictionary.
    res = Fraction((t * (b - a) + kappa + 1) * degC, t)
    gap = res - (e * degC - ordC + 1)
    gap_eff = gap_effective(gap)          # extra unit factor only if a POSITIVE int
    dio = (m + n) * q * k - n * (q * l - p)  # must equal k  (CHAIN identity)
    bord = bridge_ord(a, b, t, kappa, ordC)
    if sig[1] != bord:
        bridge_bad.append(name)
    notes = []
    if not retracts:
        notes.append("corner REFUSED -> C = y monomial")
    if A0p != (1, 0):
        notes.append("A0'!=(1,0): chart_exponent rule unvalidated here")
    if k != 1:
        notes.append("k=2: N-formula unverified")
    assert dio == k, f"{name}: Diophantine failed"
    survey.append((name, j, A0, (p, l, q), (m, n), degs, t, kappa, degC, ordC,
                   e, r, N, gap, retracts, sig, bord, notes))
    print(f"{name:4} {str(A0):>8} ({p}\\{l},{q})  {str((m,n)):>8} {str(degs):>10} "
          f"{t:>2} {kappa:>3} {degC:>3} {ordC:>2} {e:>2} {r:>2} {N:>4} "
          f"{str(gap):>5} {sig[1]:>7} {bord:>7}  {'; '.join(notes)}")

assert not bridge_bad, ("ord_y(Phi) disagrees with the bridge identity at %s"
                        % bridge_bad)
assert set(affected) == AFFECTED_EXPECTED, \
    ("the guard-refused set moved: got %s, expected %s"
     % (sorted(affected), sorted(AFFECTED_EXPECTED)))
print(f"\n  ALL {len(survey)} rows: derived ord_y(Phi) == a*q*M - H (bridge identity, "
      f"PROVED in BRIDGE_GENERALITY.md).  No row is fitted to itself.")
print(f"  GUARD-REFUSED rows ({len(affected)} of 17): {', '.join(sorted(affected, key=lambda s: int(s[1:])))}")
print("  These are exactly the rows whose pre-2026-07-27 (t, deg C, ord C) came")
print("  through the refused final-corner dictionary.  NOTE the previously")
print("  circulated affected set {F1,F2,F3,F5,F9,F10} is INCOMPLETE: it is the")
print("  refused subset of the TWELVE rows bridge_generality.py transcribes, not")
print("  of the SEVENTEEN here.  F4, F6, F11 (k=2 rows at the same (5,20)/(7,21)")
print("  corners) and F12, F13 (A0'=(2,0)) are refused too.")

# ---- TRIPWIRES: the two rows whose repair carries an extra conditional --------
_atlas_corners = {(4, 12), (5, 20), (6, 15), (7, 21), (8, 24), (9, 24)}
print("\n  TRIPWIRE A0'=(2,0):  F12 (8,24) and F13 (9,21) are guard-refused, and")
print("    their repaired chart data is therefore MORE correct than the dictionary's")
print("    -- but chart_exponent = ceil(b0/a0) is [INFERRED] and is validated only")
print("    against A0'=(1,0) published reductions.  Their t is CLAIMED, not")
print("    exact-checked.  F12's corner (8,24) is corroborated by corner_atlas.json")
print("    (rows F_22(2,3)/96 and F_24(3,4)/128: t=3, deg C = ord C = 1);")
print("    F13's corner (9,21) appears in NO atlas row and has NO second source.")
assert (9, 21) not in _atlas_corners, "if (9,21) gained an atlas row, retire this tripwire"
assert A0P_20_ROWS == {"F12", "F13"} and A0P_20_ROWS <= AFFECTED_EXPECTED

print("\nLength-2 families (composite chart underived -- t,kappa unknown): "
      + ", ".join(f"{nm} A2=({p}\\{l},{q})" for nm, A0, p, l, q in FAMILIES_LEN2))

# ---------------------------------------------------------------------------
# 1b.  THE t-CENSUS  --  what the repair does to this file's original headline.
# ---------------------------------------------------------------------------
print("\n" + "=" * 96)
print("STEP 1b -- t-CENSUS: the 'new t' out-of-sample test does NOT survive")
print("=" * 96)
_tvals = sorted({row[6] for row in survey})
print(f"  t over all 17 length-1 rows, derived from the corner: {_tvals}")
assert _tvals == [3, 4], _tvals
print("  Pre-repair this file read t = l_final and reported t in {3,4,5,7,8}, with")
print("  F9's t=7 as the out-of-sample point.  t = 7 was never a chart exponent:")
print("  chart_exponent(7,21) = 3, and GGHV22 2204.14178.tex:1394 PUBLISHES the")
print("  chart phi_3(y) = y x^3 with [P,Q] = x at that very corner.  So:")
print("    * GGV5's v11<=35 length-1 tables contain NO t outside {3,4};")
print("    * (8,28), the audited corner, is t=4 -- already in the census;")
print("    * therefore these tables offer NO t-variation beyond what (8,28) gave,")
print("      and the 'break the t/kappa correlation by varying t' programme is")
print("      CLOSED NEGATIVE on this data set, not confirmed.")
print("  A genuinely new t needs a corner with ceil(b0/a0) not in {3,4}; the atlas")
print("  has such rows ((7,35) t=5, (7,42) t=6) but they are not GGV5 family rows.")

# ---------------------------------------------------------------------------
# 2. kappa = t-2 : forced for the standard-chart class
# ---------------------------------------------------------------------------
print("\n" + "=" * 96)
print("STEP 2 -- kappa = t-2 is forced on the whole standard-chart class")
print("=" * 96)
x, ls = sp.symbols("x l_s", positive=True)
X, Y = x**-1, x**ls * y
J = sp.simplify(sp.diff(X, x) * sp.diff(Y, y) - sp.diff(X, y) * sp.diff(Y, x))
print(f"  Jacobian of (X,Y) -> (x^-1, x^l y):  {J}   (= -x^(l-2) for all l)")
print("  ell(C) = x^t c with t = l (each factor (Y - r X^-l) -> x^l (y - r)).")
print("  Hence kappa = l-2 = t-2 for EVERY corner reduced by this chart --")
print("  all 15 length-1 A0'=(1,0) families above.  t and kappa cannot be")
print("  separated inside this class; escapes need A0'=(2,0) or length-2")
print("  chains, whose charts are derived in no paper (see PHI_CORNER4.md).")

# ---------------------------------------------------------------------------
# 3. Derivations.  BOTH corners are guard-REFUSED, so both are re-derived at the
#    REPAIRED chart data (C = y a monomial) and both are checked against the
#    bridge identity -- a target this file does not produce.
#
#    F9's corner (7,21): t = 3 (GGHV22 PUBLISHES this chart), not l_final = 7.
#    F2's corner (5,20): t = 4 (GGV3 pins it via (50,75)),   not l_final = 5.
# ---------------------------------------------------------------------------
def derive(tag, a, b, A0, l_final, b_final, row=None):
    t, kappa, degC, ordC, retracts = chart_of_corner(
        A0, l_final=l_final, b_final=b_final, who="phi_corner4 derive %s" % tag)
    e = b - a + 1
    coef = t * (b - a) + kappa + 1
    rho = (e - 1) * ordC + 1
    dg = degC - ordC
    N, _, r, fit = fit_signature(a, b, t, kappa, degC, ordC)
    print("\n" + "=" * 96)
    print(f"STEP 3 -- {tag}:  corner A0={tuple(A0)}  (a,b)=({a},{b})  t={t} "
          f"kappa={kappa} deg C={degC} ord C={ordC}  e={e} r={r} rho={rho} "
          f"deg g={dg} N={N}")
    print("=" * 96)
    if not retracts:
        print(f"  corner REFUSED by the retraction guard: "
              f"{RETRACTION_LEDGER[tuple(A0)][1]}")
        print(f"  GGV5's chain row gives (l_final,b_final) = ({l_final},{b_final}); "
              f"that is NOT this corner's chart data.")
    print(f"  forcing ODE:  {a*t} c f' - {a*coef} c' f = c^{e},   c = y^{ordC} g")

    A = sp.symbols("A")
    if dg == 0:
        # C is a MONOMIAL (the guard-refused shape).  g = 1 is FORCED -- a monic
        # constant, no free coefficient to gauge and no root to place -- so the
        # generic-g collapse below is VACUOUS, not merely simpler.  The ODE
        # reduces to a*A*(t*e - coef) = 1 with t*e - coef = t - kappa - 1 = 1,
        # hence A = 1/a exactly and f = (1/a) y^e.
        g_sol = sp.Integer(1)
        print(f"  deg g = deg C - ord C = 0: C = y is a MONOMIAL, g = 1 FORCED.")
        print(f"  The common-root gauge branch is VACUOUS (no residual to place a")
        print(f"  root on).  ODE collapses to a*A*(t*e - coef) = a*A*(t-kappa-1) = "
              f"a*A = 1.")
        c_sol = y**ordC * g_sol
        f_sol = sp.expand(sp.Rational(1, a) * y**rho)
        A_sol = sp.Rational(1, a)
    else:
        # generic-g collapse
        gc = sp.symbols(f"g0:{dg+1}")
        g = sum(gc[i] * y**i for i in range(dg + 1))
        c = y**ordC * g
        f = A * y**rho * g**e
        resid = sp.expand(a * t * c * sp.diff(f, y)
                          - a * coef * sp.diff(c, y) * f - c**e)
        quo = sp.expand(sp.factor(resid) / (y**(e * ordC) * g**(e - 1)))
        print(f"  ansatz f = A y^{rho} g^{e} collapses the ODE; coefficient system "
              f"forces g_1..g_{dg-1} = 0,")
        print(f"  top coefficient resonant (free), g(-1)=0 + monic => g = y^{dg} + 1.")
        g_sol = y**dg + 1
        # solve A exactly from the constant coefficient with g = g_sol:
        A_sol = sp.solve(sp.expand(quo.subs({gc[i]: sp.Poly(g_sol, y).coeff_monomial(y**i)
                                             for i in range(dg + 1)})).coeff(y, 0), A)[0]
        c_sol = y**ordC * g_sol
        f_sol = sp.expand(A_sol * y**rho * g_sol**e)
        print(f"  g = y^{dg}+1,  H = {sp.factor(g_sol / (y + 1))},  A = {A_sol}")
    assert sp.expand(a * t * c_sol * sp.diff(f_sol, y)
                     - a * coef * sp.diff(c_sol, y) * f_sol - c_sol**e) == 0
    # UNIQUENESS, independently of the ansatz: a fully generic linear solve.
    Dmax = int(max(Fraction(coef * degC, t), rho)) + 3
    fc = sp.symbols("uf0:%d" % (Dmax + 1))
    fgen = sum(fc[i] * y**i for i in range(Dmax + 1))
    rg = sp.expand(a * t * c_sol * sp.diff(fgen, y)
                   - a * coef * sp.diff(c_sol, y) * fgen - c_sol**e)
    gsols = sp.solve(sp.Poly(rg, y).all_coeffs(), fc, dict=True)
    assert len(gsols) == 1 and sp.expand(fgen.subs(gsols[0]) - f_sol) == 0, \
        ("generic solve does not reproduce the ansatz solution", tag, len(gsols))
    print(f"  f = {sp.factor(f_sol)}    (deg {sp.degree(f_sol, y)}; UNIQUE polynomial "
          f"solution by a generic solve to degree {Dmax})")

    Phi = sp.expand(f_sol * c_sol**N)
    deg = sp.degree(Phi, y)
    ordy = min(mm[0] for mm in sp.Poly(Phi, y).monoms())
    m1 = 0
    qq = sp.Poly(Phi, y)
    d1 = sp.Poly(y + 1, y)
    while True:
        qq2, rem = sp.div(qq, d1)
        if not rem.is_zero:
            break
        qq, m1 = qq2, m1 + 1
    cof = deg - ordy - m1 * 1
    sig = (deg, ordy, m1, cof)
    print(f"  Phi = f * C^{N} = {sp.factor(Phi)}")
    print(f"  SIGNATURE (deg, ord_y, mult_(y+1), cofactor) = {sig}")
    verdict = "MATCHES" if sig == fit else "DIFFERS"
    print(f"  law prediction (parameter-free)              = {fit}   ==> {verdict}")
    # THE INDEPENDENT CHECK: a target this file does not produce.
    bord = bridge_ord(a, b, t, kappa, ordC)
    assert ordy == bord, (tag, ordy, bord)
    print(f"  bridge identity a*q*M - H = {a}*{ordC}*{t*(a+b)-(kappa+1)} - "
          f"{ordC*(a+b)-1} = {bord}   ==> ord_y AGREES  [INDEPENDENT: "
          f"BRIDGE_GENERALITY.md, PROVED]")
    if row is not None and row in SUPERSEDED:
        st, sdC, soC, sN, sordy = SUPERSEDED[row]
        print(f"  SUPERSEDED (pre-repair, do NOT cite): (t,deg C,ord C) = "
              f"({st},{sdC},{soC}), N = {sN}, ord_y = {sordy}")
        print(f"    -> N moved {sN} -> {N}, ord_y moved {sordy} -> {ordy}; the old "
              f"values are what the refused dictionary produced.")
    return sig, fit, verdict

sig9, fit9, v9 = derive("F9 j=0, case (56,84)", 2, 3, (7, 21), 7, 2, row="F9")
sig2, fit2, v2 = derive("F2 j=0, case (50,75)", 2, 3, (5, 20), 5, 2, row="F2")

# ---------------------------------------------------------------------------
# 4. Verdict
# ---------------------------------------------------------------------------
print("\n" + "=" * 96)
print("VERDICT")
print("=" * 96)
print(f"  F9 j=0 (56,84), corner (7,21), t=3 : derived {sig9}  law {fit9}  -> {v9}")
print(f"  F2 j=0 (50,75), corner (5,20), t=4 : derived {sig2}  law {fit2}  -> {v2}")
print("""  kappa = t-2 is FORCED on the standard single-chart class (all length-1
  A0'=(1,0) families), and that is UNTOUCHED by the repair: it follows from the
  chart's Jacobian -x^(l-2) for whatever l the chart uses.  Within the class the
  law has NO free parameters left.

  WHAT THE 2026-07-27 REPAIR CHANGED, stated plainly for a reader of the public
  tree:

    * The chart exponent is DERIVED from the corner (ceil(b0/a0)), never read
      off GGV5's final chain corner.  Eleven of the seventeen length-1 rows --
      F1-F6, F9-F13 -- sit on corners that do NOT retract, and there C is the
      MONOMIAL y (deg C = ord C = 1), not y^b_final * (residual of degree
      a0 - b_final).  Their old (t, deg C, ord C, N, signature) are SUPERSEDED.
    * This file's original headline is WITHDRAWN, not weakened.  'F9 gives the
      new value t = 7' is false at its premise: GGHV22 publishes l = 3 and
      [P,Q] = x at the corner (7,21) (2204.14178.tex:1394).  After the repair
      every GGV5 v11<=35 length-1 corner has t in {3,4} (STEP 1b), so these
      tables offer no t-variation beyond the audited (8,28) t=4.
    * The two derivations above are now checked against ord_y(Phi) = a*q*M - H,
      PROVED independently in BRIDGE_GENERALITY.md.  Before the repair this file
      compared its own superseded model to its own superseded targets and
      reported MATCHES; that is the specific failure the bridge check removes.
    * (72,108) and (108,144) are at (8,28), which DOES retract (it is the unique
      t=4 corner in GGV5's list with b0 = 4a0-4), so the audited home case and
      the corner-144 template are UNAFFECTED.  So are F7, F8, F14, F15, F16, F17
      -- corners (6,15) and (9,24) both retract.
    * The only routes to kappa != t-2 remain the length-2 chains (F18-F24), whose
      reduction charts exist in no paper.  A0'=(2,0) (F12, F13) is NOT such a
      route -- one final inversion still gives kappa = l-2 (composite_charts.py)
      -- but its l is CLAIMED, see the STEP 1 tripwire.""")
