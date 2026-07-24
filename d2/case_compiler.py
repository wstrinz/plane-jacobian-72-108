#!/usr/bin/env python3
"""case_compiler.py — pilot CASE COMPILER (Lane I; INDUCTIVE_PROGRAM.md architecture item).

Input : a GGV5 family case (family id + j) or a special externally-sourced case
        (GGHV22 (72,108)/(108,144)).
Output: a compiled CASE DOSSIER — canonical JSON + human-readable rendering —
        containing:
  (1) the corner signature (a,b,t,kappa=t-2,a0,q,e,r,gap,dg,N), Diophantine-checked;
  (2) the predicted Phi signature from the unified seven-point corner law
      (PHI_F14.md), flagged CONJECTURAL when the case sits in a regime without a
      derived point (data-driven REGIME_STATUS table below — update it when a
      new regime point lands, e.g. Lane H's F7 test of gap>0 & r>0);
  (3) the master-identity layout in BOTH presentations (INDUCTIVE_PROGRAM.md
      amendment): eliminated / f31-style for the tropical & infinity layers, and
      the pre-resultant G-system for terminal decisions.  STRUCTURE is emitted;
      full polynomial data exists only for (72,108) (there is no f31-analogue
      elsewhere), so every field carries an explicit instantiated-vs-schematic
      marker;
  (4) the Galois transfer check (GALOIS_LIBRARY.md sec. 4 two-line rule) applied
      to the case's candidate forcing polynomial;
  (5) the per-mechanism transfer inventory of the (72,108) machinery.

Honesty contract: nothing here is a new mathematical claim.  Instantiated
fields restate landed, audited/verified facts with sources; schematic fields
describe structure and say exactly what computation would instantiate them;
conjectural flags are raised, never silently dropped.

Checker: case_compiler_verify.py (validates the (72,108), (75,125)=F2 j=1 and
F9 j=0 dossiers against known landed values).  Sources: GGV5 family tables as
transcribed and Diophantine-checked in phi_corner4.py; PHI_75_125.md,
PHI_CORNER4.md, PHI_F14.md (corner law, 7 exact points); GALOIS_LIBRARY.md
(transfer rule); FULL_SYSTEM_BRIDGE.md (G-system, variable dictionary);
WINDOW_CAPS.md (window caps); STATE.md (audited (72,108) facts).
"""
import json
import sys
from fractions import Fraction
from math import gcd

import sympy as sp

y = sp.symbols("y")

# ---------------------------------------------------------------------------
# Case data.
# GGV5 v11<=35 length-1 family table, transcribed verbatim from the GGV5 final
# section exactly as in phi_corner4.py (same tuple layout:
# (name, A0, A0', p, l, q, k, (m0,dm), (n0,dn)), (m,n)(j) = (m0+dm*j, n0+dn*j),
# final corner A1 = (p\l, q)).  Diophantine identity (m+n)*q*k - n*(q*l-p) = k
# is enforced for every instantiated j.
# ---------------------------------------------------------------------------
FAMILIES_LEN1 = [
    ("F1",  (4, 12), (1, 0),  7, 4, 3, 1, (3, 2),  (4, 3)),
    ("F2",  (5, 20), (1, 0),  7, 5, 2, 1, (2, 1),  (3, 2)),
    ("F3",  (5, 20), (1, 0),  8, 5, 3, 1, (3, 4),  (2, 3)),
    ("F4",  (5, 20), (1, 0),  8, 5, 3, 2, (3, 2),  (16, 12)),
    ("F5",  (5, 20), (1, 0),  9, 5, 4, 1, (9, 7),  (5, 4)),
    ("F6",  (5, 20), (1, 0),  9, 5, 4, 2, (4, 3),  (10, 8)),
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
FAMILY_INDEX = {row[0]: row for row in FAMILIES_LEN1}

# Externally-sourced cases (chain data NOT from the GGV5 v11<=35 table; the
# GGV5 Diophantine identity is not asserted for them).
SPECIAL_CASES = {
    "GGHV_72_108": dict(
        A0=(8, 28), pair=(2, 3), t=4, a0=8, q=7, degrees=(72, 108),
        source="GGHV22 arXiv:2204.14178 Prop 4.3 case (8,28); l=4 chart at "
               "2204 lines 1228-1230 (PHI_CORNER4.md); audited campaign case",
        audited=True),
    "GGHV_108_144": dict(
        A0=(8, 28), pair=(3, 4), t=4, a0=8, q=3, degrees=(108, 144),
        source="GGHV22 (8,28) corner, second pair; CORNER_144_COMPARISON "
               "template case (PHI_75_125.md / PHI_F14.md table)",
        audited=False),
}

# ---------------------------------------------------------------------------
# Regime status — DATA-DRIVEN; update when a new regime point lands.
# A regime key is (gap_positive, r_positive).  Lane H is currently testing F7
# (42,147), the smallest gap>0 & r>0 case; if it lands MATCHES, flip that row
# to grounded and cite PHI_F7.
# ---------------------------------------------------------------------------
REGIME_STATUS = {
    (False, False): dict(name="gap0_r0", grounded=True,
                         evidence="no such length-1 survey row uses it alone; "
                                  "covered by the unified law's gap=0 limit"),
    (False, True):  dict(name="gap0_r_pos", grounded=True,
                         evidence="derived: (108,144), (75,125), F9 (56,84), "
                                  "F2 j=0 (50,75), F14 (66,231) — PHI_75_125/"
                                  "PHI_CORNER4/PHI_F14, all MATCH"),
    (True,  False): dict(name="resonance_gap_r0", grounded=True,
                         evidence="audited (72,108) + derived F1 (48,64) "
                                  "(PHI_F14.md: unit cofactor of degree gap)"),
    (True,  True):  dict(name="gap_pos_r_pos", grounded=True,
                         evidence="derived RAMIFIED law at F7 (42,147), F3, "
                                  "F10, F16 (PHI_F7.md): the old conjecture "
                                  "mult=e+N / cof=gap+r(e+N) is REFUTED "
                                  "(simple root at -1 impossible, proven at "
                                  "dg=2); amended: mult = dg(e+N)-(dg-1), "
                                  "cofactor = gap+r (retro-explains the "
                                  "(72,108) quartic, 4 = gap+r). "
                                  "BRANCH-CONDITIONAL: ramified branch "
                                  "selected by continuity with (y+1)|C; the "
                                  "complex-pair branch (mult 0) is recorded, "
                                  "not excluded (PHI_F7.md judgment)"),
}

# Known derived/audited corner-law points: tag -> (signature, source).
KNOWN_POINTS = {
    "GGHV_72_108":   ((238, 204, 30, 4),   "STATE.md (audited) / PHI_F14.md"),
    "GGHV_108_144":  ((550, 205, 69, 276), "PHI_75_125.md / corner-144"),
    "F2_j1_75_125":  ((504, 201, 101, 202), "PHI_75_125.md"),
    "F9_j0_56_84":   ((377, 107, 54, 216), "PHI_CORNER4.md"),
    "F2_j0_50_75":   ((189, 75, 38, 76),   "PHI_CORNER4.md"),
    "F14_j0_66_231": ((375, 165, 42, 168), "PHI_F14.md"),
    "F1_j0_48_64":   ((275, 205, 69, 1),   "PHI_F14.md"),
}

# Known forcing polynomials (gap-regime unit cofactors) and leading constants,
# from the landed derivations.  Registry, not derivation.
KNOWN_FORCING = {
    "GGHV_72_108": "2048*y**4 - 512*y**3 + 320*y**2 - 240*y + 195",
    "F1_j0_48_64": "4*y - 1",
}
KNOWN_LC = {
    "GGHV_72_108": "-1/6630", "F2_j1_75_125": "-1/9", "F9_j0_56_84": "-1/10",
    "F2_j0_50_75": "-1/6", "F14_j0_66_231": "-1/10", "F1_j0_48_64": "1/15",
}

KILL_CLASSES = {"C08": 105, "C20": 170}   # GALOIS_LIBRARY.md census

# ---------------------------------------------------------------------------
# Transfer inventory — per-mechanism status of the (72,108) machinery.
#   AS-IS       : transfers with no per-case re-derivation (stated caveats aside)
#   PARAMETRIC  : structure transfers; numeric data must be recomputed from case data
#   METHOD-ONLY : the technique transfers; the killing arithmetic is instance-specific
# ---------------------------------------------------------------------------
TRANSFER_INVENTORY = [
    dict(mechanism="tropical/infinity cascade engine", status="PARAMETRIC",
         rationale="engines run window-generically; kill laws affine in the "
                   "branch parameter (t-depth 30-3a at (72,108)); sub1 26-family "
                   "was a-independent until low-a infinity effects (understood)",
         source="INDUCTIVE_PROGRAM.md layer 1; cascade_engine.py"),
    dict(mechanism="T2 squeeze (F^2|G degree squeeze)", status="PARAMETRIC",
         rationale="the squeeze inequality is presentation-independent; the "
                   "terminal law m_i = 6+2s_i-3b_i gating it is case arithmetic",
         source="cascade engine C24; test_cascade_inf.py R5"),
    dict(mechanism="residue library + Galois descent", status="AS-IS",
         rationale="support-geometry classification of all 23 shapes is "
                   "corner-independent; only case inputs are the forcing "
                   "polynomial's Galois label + disc class (two-line check); "
                   "caveat G1: free-torus scope",
         source="GALOIS_LIBRARY.md sec. 2, 4"),
    dict(mechanism="corner law (Phi signature)", status="AS-IS",
         rationale="seven exact points, t in {3,4,5,7}; AS-IS within the "
                   "standard-chart length-1 class; CONJECTURAL for gap>0 & r>0, "
                   "A0'=(2,0) (F12,F13) and length-2 chains (F18-F24)",
         source="PHI_F14.md unified law"),
    dict(mechanism="full-system bridge (pre-resultant G-system)", status="PARAMETRIC",
         rationale="corner-144 recurrence: skeleton yes, numerics no; weights "
                   "and the ~122-equation bridge need the case's D-transform",
         source="FULL_SYSTEM_BRIDGE.md; INDUCTIVE_PROGRAM.md test T2"),
    dict(mechanism="window caps (ord/deg caps on D_k)", status="PARAMETRIC",
         rationale="the three T3 valuation inductions close as symbolic "
                   "identities in k; the slopes (12k/15k/14k at (72,108)) come "
                   "from case degree arithmetic",
         source="WINDOW_CAPS.md"),
    dict(mechanism="divisor-lemma engine (confluent Vandermonde)", status="PARAMETRIC",
         rationale="rank/dimension statements with resultant-product "
                   "determinants are polynomial in the case data",
         source="INDUCTIVE_PROGRAM.md layer 3; DIVISOR_LEMMAS"),
    dict(mechanism="S-unit / Mason-Stothers corner kills", status="METHOD-ONLY",
         rationale="the height gate transfers as a method; the 17-does-not-"
                   "divide-5 coprimality corner is instance arithmetic",
         source="S_UNIT_LAYER"),
    dict(mechanism="modular triage (GB mod p)", status="AS-IS",
         rationale="pipeline is case-agnostic; only the bad-prime set is "
                   "case-dependent",
         source="MODULAR_TRIAGE"),
]

# ---------------------------------------------------------------------------
# Arithmetic helpers
# ---------------------------------------------------------------------------

def square_class(rat):
    """Squarefree part (with sign) of a nonzero rational — its class in Q*/(Q*)^2."""
    r = sp.Rational(rat)
    assert r != 0
    n = int(r.p) * int(r.q)
    sign = -1 if n < 0 else 1
    cls = 1
    for p, e in sp.factorint(abs(n)).items():
        if e % 2:
            cls *= p
    return sign * cls


def galois_label(poly_expr):
    """Galois data of a forcing-polynomial candidate over Q (degree <= 4).

    Returns dict with: degree, irreducible, label, disc, disc_class,
    quadratic_subfield_classes (None if witness-decided), method.
    Quartic labeling: resolvent cubic + discriminant + factorization over
    Q(sqrt(disc_class)) for the C4/D4 split (standard criteria, computed
    exactly — no table lookups)."""
    P = sp.Poly(sp.expand(poly_expr), y)
    d = P.degree()
    if d == 0:
        return dict(degree=0, irreducible=False, label="TRIVIAL", disc=None,
                    disc_class=None, quadratic_subfield_classes=[],
                    method="constant")
    if d == 1:
        return dict(degree=1, irreducible=True, label="C1", disc=None,
                    disc_class=1, quadratic_subfield_classes=[],
                    method="linear: splitting field Q")
    D = sp.discriminant(P.as_expr(), y)
    Dc = square_class(D)
    irr = P.is_irreducible
    if d == 2:
        if not irr:
            return dict(degree=2, irreducible=False, label="C1", disc=str(D),
                        disc_class=1, quadratic_subfield_classes=[],
                        method="reducible quadratic: splits over Q")
        return dict(degree=2, irreducible=True, label="C2", disc=str(D),
                    disc_class=Dc, quadratic_subfield_classes=[Dc],
                    method="irreducible quadratic: splitting field Q(sqrt(%d))" % Dc)
    if d == 3:
        label = "C3" if Dc == 1 else "S3"
        return dict(degree=3, irreducible=bool(irr), label=label if irr else "REDUCIBLE",
                    disc=str(D), disc_class=Dc,
                    quadratic_subfield_classes=None if not irr else ([] if Dc == 1 else [Dc]),
                    method="cubic: disc square test")
    assert d == 4, "forcing-candidate degrees above 4 not needed by the pilot"
    if not irr:
        return dict(degree=4, irreducible=False, label="REDUCIBLE", disc=str(D),
                    disc_class=Dc, quadratic_subfield_classes=None,
                    method="reducible quartic: analyze factors instead")
    Pm = P.monic()
    c4, p3, p2, p1, p0 = [sp.Rational(c) for c in Pm.all_coeffs()]
    assert c4 == 1
    z = sp.symbols("z")
    resolvent = z**3 - p2 * z**2 + (p3 * p1 - 4 * p0) * z \
        - (p3**2 * p0 - 4 * p2 * p0 + p1**2)
    lin = [f for f, m in sp.factor_list(resolvent, z)[1] for _ in range(m)
           if sp.Poly(f, z).degree() == 1]
    n_rat = len(lin)
    if n_rat == 0:
        label = "A4" if Dc == 1 else "S4"
        subf = [] if Dc == 1 else [Dc]
        method = "irreducible resolvent cubic; disc class %s" % Dc
    elif n_rat == 3:
        label, subf = "V4", None
        method = "resolvent cubic splits; V4 (witness-decided per GALOIS_LIBRARY G2)"
    else:
        # one rational resolvent root, disc nonsquare: C4 iff the quartic
        # factors over Q(sqrt(disc_class)), else D4.
        assert Dc != 1
        factors = sp.factor_list(Pm.as_expr(), extension=sp.sqrt(Dc))[1]
        if sum(m for _, m in factors) > 1:
            label, subf = "C4", [Dc]
            method = "one rational resolvent root; factors over Q(sqrt(%d)) -> C4" % Dc
        else:
            label, subf = "D4", None
            method = "one rational resolvent root; irreducible over Q(sqrt(%d)) -> D4" % Dc
    return dict(degree=4, irreducible=True, label=label, disc=str(D),
                disc_class=Dc, quadratic_subfield_classes=subf, method=method)


def transfer_verdicts(label, disc_class):
    """GALOIS_LIBRARY.md sec.4 two-line rule per kill shape.

    S4/C4 (unique quadratic subfield Q(sqrt(disc))), C2, C1: kill iff the
    obstruction class survives, i.e. disc class != Delta.  A4: no quadratic
    subfield -> always kills.  D4/V4: witness-decided (rule insufficient)."""
    out = {}
    for name, delta in KILL_CLASSES.items():
        if label in ("S4", "C4", "C2", "C1", "C3", "S3", "TRIVIAL"):
            adjoined = {1}
            if label in ("S4", "C4", "C2") and disc_class is not None:
                adjoined.add(disc_class)
            if label == "S3" and disc_class is not None:
                adjoined.add(disc_class)
            out[name] = "KILLS" if delta not in adjoined else \
                "OBSTRUCTION-VANISHES (sqrt(%d) in splitting field)" % delta
        elif label == "A4":
            out[name] = "KILLS"
        else:  # D4, V4, REDUCIBLE
            out[name] = "WITNESS-NEEDED (D4/V4 or reducible: rule insufficient, "\
                        "GALOIS_LIBRARY.md G2)"
    return out


# ---------------------------------------------------------------------------
# Corner signature + law
# ---------------------------------------------------------------------------

def corner_signature(a, b, t, a0, q):
    """(a,b,t,kappa,a0,q,e,r,gap,dg,N) with kappa = t-2 (structural,
    PHI_CORNER4.md) and the PHI_F14.md mini-lemmas."""
    kappa = t - 2
    e = b - a + 1
    r = a0 - q - 1
    N = a * (t * (a + b - 1) + 1) - 2 * b
    gap = Fraction(q - 1, 1) - Fraction(a0, t)
    dg = a0 - q
    return dict(a=a, b=b, t=t, kappa=kappa, a0=a0, q=q, e=e, r=r,
                gap_num=gap.numerator, gap_den=gap.denominator,
                gap=(int(gap) if gap.denominator == 1 else None),
                dg=dg, N=N)


def law_signature(sig):
    """Corner law.  deg/ord hold on every branch (PHI_F7.md); mult/cofactor
    split by regime: unramified (PHI_F14.md, seven points) vs the ramified
    gap>0,r>0 regime (PHI_F7.md, four points, ramified-branch selection).
    Requires integral gap."""
    assert sig["gap_den"] == 1, "non-integral gap: law prediction undefined"
    g, e, N, a0, q, r = sig["gap"], sig["e"], sig["N"], sig["a0"], sig["q"], sig["r"]
    dg = sig["dg"]
    deg = (e * a0 - q + 1) + g + N * a0
    ordy = (e - 1) * q + 1 + N * q
    if g > 0 and r > 0:
        mult = dg * (e + N) - (dg - 1)
        cof = g + r
        branch = ("ramified (g has a forced double root at -1; selected by "
                  "continuity with (y+1)|C — complex-pair branch recorded in "
                  "PHI_F7.md, not excluded)")
    else:
        mult = e + N
        cof = g + r * mult
        branch = "unramified (simple root at -1; all seven derived points)"
    return dict(deg=deg, ord_y=ordy, mult_y_plus_1=mult, cofactor_deg=cof,
                branch=branch)


# ---------------------------------------------------------------------------
# Presentations (item 3): BOTH, with instantiated-vs-schematic markers
# ---------------------------------------------------------------------------

def build_presentations(tag, sig, law):
    inst = (tag == "GGHV_72_108")
    mult, cof = law["mult_y_plus_1"], law["cofactor_deg"]
    dg = sig["dg"]
    H = "1" if dg <= 1 else str(sp.expand(sp.quo(y**dg + 1, y + 1))) \
        if dg % 2 == 1 else None
    phi_norm = dict(
        form="Phi_stripped = Phi_full / (lc * y^ord) = t^mult * u_gap * H^mult, "
             "t = y+1, deg u_gap = gap, H = (y^dg+1)/(y+1)",
        mult=mult, deg_u_gap=sig["gap"], H=H,
        H_note=None if H is not None else
            "dg even: (y+1) does not divide y^dg+1, and PHI_F7.md proves a "
            "simple root at -1 is impossible here — the residual is the "
            "RAMIFIED shape (double root at -1, e.g. g=(y+1)^2 at dg=2), not "
            "(y^dg+1)/(y+1); see the gap_pos_r_pos regime law",
        lc=KNOWN_LC.get(tag),
        total_stripped_degree=mult + cof,
        note=("equals the audited master-identity Phi = c*t^30*q4 exactly "
              "(FULL_SYSTEM_BRIDGE.md sec.1: Phi_full/y^204 = c*t^30*q, verified)"
              if inst else
              "instantiated from the corner law; the underlying Phi_full = f*C^N"))
    eliminated = dict(
        purpose="tropical / infinity / discovery layers (window-generic engines)",
        instantiated=inst,
        master_identity_shape="f_D = sum_{f=0}^{F} Phi^f * e^(E - s*f) * h_f == 0 in Q[y]",
        parameters=(dict(D=31, F=7, E=21, s=3) if inst else None),
        schematic_note=(None if inst else
            "the tower data (D,F,E,s) requires the case's cascade construction; "
            "no f31-analogue exists off (72,108) (corner-144 recurrence: "
            "skeleton yes, numerics no — INDUCTIVE_PROGRAM.md T2). No numeric "
            "guess is offered; building the case's C-series/D-transform is the "
            "instantiating computation."),
        phi_normalized=phi_norm,
        cascade_chain_shape=dict(
            shape="t^v * g_(l+1) = ehat^3 * g_l + u^l * h_l; anchor t^v g_1 = h_0; "
                  "v = 30 - 3*a_branch at (72,108)",
            instantiated=inst,
            note="a_branch is the cascade branch parameter (NOT the corner 'a'); "
                 "the depth constant (30 = mult at (72,108)) and step are case "
                 "data — schematic off (72,108)"),
        infinity_layer="v_inf = -deg; unique-max forcing with recorded tie "
                       "obligations — engine transfers per INDUCTIVE_PROGRAM "
                       "layer 1 (PARAMETRIC)")
    g_system = dict(
        purpose="terminal decisions (the presentation holding strictly more "
                "information than the eliminated form)",
        instantiated=inst,
        generators=(["G1", "G2", "G3", "G5body+Phi"] if inst else None),
        generator_skeleton="4-generator pre-resultant ideal; elimination ideal "
                           "in the deep window variables is principal "
                           "(= <f_D>); corner-144: skeleton recurs, numerics do not",
        ring=("Q[d~2,d~1,d~0,e,r,s,dm4,Phi]" if inst else
              "deep window coefficients + e + spare unknowns "
              "(negative-index d's) + Phi"),
        spare_dictionary="r, s, dm4 = d_(-2), d_(-3), d_(-4) "
                         "(FULL_SYSTEM_BRIDGE.md sec.1, verified V1)",
        G_weights=([156, 168, 180, 204] if inst else None),
        bridge=("~122 quadratic equations, 45 sub2 spare coefficients" if inst
                else "schematic: the state->G-system bridge equations need the "
                     "case's D-transform"),
        window_caps=dict(
            pattern="ord >= 12k; deg <= 15k (sub1) / 14k (sub2); k = 2..8",
            status=("PROVEN (WINDOW_CAPS.md, 81 checks; k=6,7,8 flags retired)"
                    if inst else
                    "slopes are case degree arithmetic (schematic); the T3 "
                    "inductions close symbolically in k, so the METHOD "
                    "transfers (WINDOW_CAPS.md)")))
    return dict(
        amendment="BOTH presentations emitted (INDUCTIVE_PROGRAM.md adopted "
                  "architecture): eliminated form for tropical/discovery, "
                  "pre-resultant for terminal decisions — each layer wants its "
                  "own presentation",
        eliminated_f31_style=eliminated,
        pre_resultant_G_system=g_system)


# ---------------------------------------------------------------------------
# Galois transfer section (item 4)
# ---------------------------------------------------------------------------

def build_galois_section(tag, sig):
    dg, gap_int = sig["dg"], sig["gap"]
    conditional = (tag != "GGHV_72_108")
    if tag in KNOWN_FORCING:
        poly = sp.sympify(KNOWN_FORCING[tag])
        rationale = ("registry: audited resonance-gap unit cofactor (deg = gap)"
                     if tag == "GGHV_72_108" else
                     "registry: derived resonance-gap unit cofactor (deg = gap)")
    elif gap_int == 0 and dg % 2 == 1 and dg >= 2:
        poly = sp.expand(sp.quo(y**dg + 1, y + 1))
        rationale = ("computed: residual H = (y^dg+1)/(y+1) — the forcing-"
                     "polynomial correspondence of INDUCTIVE_PROGRAM.md "
                     "(disc-17 quartic at (72,108); 10th cyclotomic at (108,144))")
    elif gap_int == 0 and dg == 1:
        poly = sp.Integer(1)
        rationale = "trivial residual (dg=1) and gap=0: no forcing candidate"
    else:
        poly = None
        rationale = ("UNKNOWN: gap>0 unit cofactor not in registry (requires "
                     "the case's ODE solve)" if (gap_int is None or gap_int > 0)
                     else "UNKNOWN: dg even — residual mechanism underived")
    if poly is None:
        return dict(forcing_candidate=dict(poly=None, rationale=rationale),
                    galois=None, verdicts={k: "UNKNOWN (no forcing candidate)"
                                           for k in KILL_CLASSES},
                    rule="GALOIS_LIBRARY.md sec.4", status="UNKNOWN")
    glab = galois_label(poly)
    verd = transfer_verdicts(glab["label"], glab["disc_class"])
    status = ("AUDITED (this is the home case: kills C08+C20 are proved here)"
              if not conditional else
              "CONDITIONAL: rule assumes the case's residue layer reproduces "
              "the (72,108) shape library on its free torus (GALOIS_LIBRARY "
              "sec.4 + judgment G1); the residue analogue is underived off "
              "(72,108)")
    return dict(
        forcing_candidate=dict(poly=str(sp.expand(poly)),
                               degree=int(sp.Poly(poly, y).degree()) if poly != 1 else 0,
                               rationale=rationale),
        galois=glab, verdicts=verd,
        rule="kill transfers iff Galois label in {S4,C4,A4} (or any label whose "
             "splitting field misses sqrt(Delta)) and disc class avoids "
             "{105, 170}; D4/V4 witness-decided",
        status=status)


# ---------------------------------------------------------------------------
# The compiler
# ---------------------------------------------------------------------------

def compile_case(name, j=None):
    """Compile a case dossier.  name = family id (with j) or a SPECIAL_CASES tag."""
    judgment = []
    if name in SPECIAL_CASES:
        scd = SPECIAL_CASES[name]
        a, b = scd["pair"]
        t, a0, q = scd["t"], scd["a0"], scd["q"]
        tag = name
        case_block = dict(
            tag=tag, family=None, j=None, m=a, n=b, a=a, b=b,
            degrees=list(scd["degrees"]), v11=sum(scd["A0"]),
            A0=list(scd["A0"]), A1=None, A0prime=None, k=None,
            chain_length=1, source=scd["source"],
            diophantine="N/A (chain data externally sourced, not a GGV5 "
                        "v11<=35 table row)")
        if scd.get("audited"):
            judgment.append("[audited] this case's corner facts and Phi are "
                            "the campaign's proven ground truth (STATE.md)")
    else:
        assert name in FAMILY_INDEX, "unknown family %r" % name
        assert j is not None, "family case needs j"
        fam, A0, A0p, p, l, q, k, (m0, dm), (n0, dn) = FAMILY_INDEX[name]
        m, n = m0 + dm * j, n0 + dn * j
        assert gcd(m, n) == 1, \
            "%s j=%d gives non-coprime (m,n)=(%d,%d): not a reduced case" % (fam, j, m, n)
        v11 = A0[0] + A0[1]
        degrees = (v11 * m, v11 * n)
        a, b = sorted((m, n))
        a0, t = A0[0], l
        dio = (m + n) * q * k - n * (q * l - p)
        assert dio == k, "%s: Diophantine identity failed" % fam
        tag = "%s_j%d_%d_%d" % (fam, j, degrees[0], degrees[1])
        case_block = dict(
            tag=tag, family=fam, j=j, m=m, n=n, a=a, b=b,
            degrees=list(degrees), v11=v11, A0=list(A0),
            A1="(%d\\%d,%d)" % (p, l, q), A0prime=list(A0p), k=k,
            chain_length=1,
            source="GGV5 v11<=35 length-1 table (transcription Diophantine-"
                   "checked; phi_corner4.py judgment 1)",
            diophantine="(m+n)*q*k - n*(q*l - p) = %d = k  OK" % dio)
        if A0p != (1, 0):
            judgment.append("[conjectural] A0' != (1,0): chart settled "
                            "(COMPOSITE_CHARTS.md, kappa = t-2 extends); "
                            "zeta-corrected tail theory (ZETA_TAIL.md) "
                            "enumerates the surviving models — every "
                            "motivated defect 0<|eta|<=1 is "
                            "NO-POLYNOMIAL-SOLUTION; survivors are eta=0 "
                            "(mu-ladder rungs mu=1,2,3, signatures "
                            "(814,506,102,206)/(814,506,203,105)/"
                            "(814,506,304,4) for F12), the eta=+2 canonical "
                            "collapse ((1292,806,162,324), support {1,3} so "
                            "mu=2 absent there), and the viable-unmotivated "
                            "eta in {-2,-3} rows left OPEN in ZETA_TAIL's "
                            "sweep; model selection needs the actual polygon "
                            "reduction")
        if k != 1:
            judgment.append("[conjectural] k = %d != 1: N-formula unverified "
                            "at k=2 (phi_corner4.py survey note)" % k)

    sig = corner_signature(a, b, t, a0, q)
    regime_key = ((sig["gap"] is None or sig["gap"] > 0), sig["r"] > 0)
    regime = REGIME_STATUS[regime_key]
    conjectural_reasons = [jt for jt in judgment if jt.startswith("[conjectural]")]
    if not regime["grounded"]:
        conjectural_reasons.append(
            "[conjectural] regime %s has no derived point: %s"
            % (regime["name"], regime["evidence"]))
    if sig["gap"] is None:
        conjectural_reasons.append(
            "[conjectural] gap = %d/%d is non-integral: the law's degree "
            "formula is undefined as stated" % (sig["gap_num"], sig["gap_den"]))
        law = None
    else:
        law = law_signature(sig)
    known = KNOWN_POINTS.get(tag)
    if known is not None and law is not None:
        lawtuple = (law["deg"], law["ord_y"], law["mult_y_plus_1"],
                    law["cofactor_deg"])
        assert lawtuple == known[0], \
            "law disagrees with landed point %s: %s vs %s" % (tag, lawtuple, known[0])

    # standing judgment items (mirroring the PHI_* docs)
    if tag != "GGHV_72_108":
        judgment.append("[judgment] unreduced polygon: the corner's explicit "
                        "reduction is performed in no paper; standard type-II.b "
                        "root shift + Laurent chart assumed (t=l, kappa=l-2, "
                        "deg C=a0) — PHI_75_125/PHI_CORNER4/PHI_F14 shared "
                        "conditional boundary")
        slice_ok = (b - 1) % a == 0
        judgment.append("[judgment] N-formula: (b-1)/a %s integral, so the "
                        "corner-144 forcing-slice picture %s"
                        % ("is" if slice_ok else "is NOT",
                           "transfers verbatim (less-conditional class)"
                           if slice_ok else
                           "does not transfer directly (more-conditional, as "
                           "for (75,125))"))
    judgment.append("[judgment] forcing-polynomial identification for "
                    "non-audited cases follows the INDUCTIVE_PROGRAM.md "
                    "correspondence; only (72,108) is audited")

    phi_block = dict(
        signature=law,
        regime=regime["name"], regime_grounded=regime["grounded"],
        regime_evidence=regime["evidence"],
        conjectural=bool(conjectural_reasons),
        conjectural_reasons=conjectural_reasons,
        derived_reference=(dict(signature=list(known[0]), source=known[1])
                           if known else None))

    dossier = dict(
        schema="case-dossier-v1",
        case=case_block,
        corner_signature=sig,
        phi_prediction=phi_block,
        presentations=build_presentations(tag, sig, law) if law else dict(
            note="law undefined (non-integral gap); presentations skipped"),
        galois_transfer=build_galois_section(tag, sig),
        transfer_inventory=TRANSFER_INVENTORY,
        judgment=judgment)
    return dossier


def canonical_json(dossier):
    return json.dumps(dossier, sort_keys=True, indent=1) + "\n"


def render(dossier):
    """Human-readable dossier summary."""
    c, s = dossier["case"], dossier["corner_signature"]
    p = dossier["phi_prediction"]
    g = dossier["galois_transfer"]
    L = []
    L.append("=" * 78)
    L.append("CASE DOSSIER  %s   degrees %s   [%s]"
             % (c["tag"], tuple(c["degrees"]),
                "CONJECTURAL" if p["conjectural"] else "grounded regime"))
    L.append("=" * 78)
    L.append("corner: A0=%s A1=%s (a,b)=(%d,%d) t=%d kappa=%d a0=%d q=%d | "
             "e=%d r=%d gap=%s dg=%d N=%d"
             % (tuple(c["A0"]), c["A1"], s["a"], s["b"], s["t"], s["kappa"],
                s["a0"], s["q"], s["e"], s["r"],
                s["gap"] if s["gap"] is not None
                else "%d/%d" % (s["gap_num"], s["gap_den"]), s["dg"], s["N"]))
    L.append("diophantine: %s" % c["diophantine"])
    if p["signature"]:
        sg = p["signature"]
        L.append("law Phi signature (deg, ord, mult, cof) = (%d, %d, %d, %d)   "
                 "regime %s (%s)"
                 % (sg["deg"], sg["ord_y"], sg["mult_y_plus_1"],
                    sg["cofactor_deg"], p["regime"],
                    "grounded" if p["regime_grounded"] else "NO DERIVED POINT"))
        if p["derived_reference"]:
            L.append("   derived/audited reference: %s  [%s]"
                     % (tuple(p["derived_reference"]["signature"]),
                        p["derived_reference"]["source"]))
    for reason in p["conjectural_reasons"]:
        L.append("   ! %s" % reason)
    fc = g["forcing_candidate"]
    L.append("galois transfer: forcing = %s  (%s)" % (fc["poly"], fc["rationale"]))
    if g["galois"]:
        L.append("   label %s, disc class %s;  C08: %s;  C20: %s"
                 % (g["galois"]["label"], g["galois"]["disc_class"],
                    g["verdicts"]["C08"], g["verdicts"]["C20"]))
        L.append("   status: %s" % g["status"])
    pres = dossier["presentations"]
    if "eliminated_f31_style" in pres:
        el, gs = pres["eliminated_f31_style"], pres["pre_resultant_G_system"]
        L.append("presentations: eliminated %s | G-system %s"
                 % ("INSTANTIATED" if el["instantiated"] else "schematic",
                    "INSTANTIATED" if gs["instantiated"] else "schematic"))
        pn = el["phi_normalized"]
        L.append("   Phi_stripped = t^%d * u_gap(deg %s) * H^%d,  H = %s"
                 % (pn["mult"], pn["deg_u_gap"], pn["mult"], pn["H"]))
    L.append("transfer inventory: " + "; ".join(
        "%s=%s" % (t["mechanism"].split(" (")[0], t["status"])
        for t in dossier["transfer_inventory"]))
    for jt in dossier["judgment"]:
        L.append("judgment: %s" % jt)
    return "\n".join(L)


# ---------------------------------------------------------------------------
# Pilot main: compile the three validation dossiers (+ CLI for any case)
# ---------------------------------------------------------------------------
PILOT = [("GGHV_72_108", None), ("F2", 1), ("F9", 0)]

def main(argv):
    if len(argv) >= 2 and argv[1] not in ("--pilot",):
        name = argv[1]
        j = int(argv[2]) if len(argv) >= 3 else None
        d = compile_case(name, j)
        print(render(d))
        out = "case_dossier_%s.json" % d["case"]["tag"]
        with open(out, "w", encoding="utf-8") as fh:
            fh.write(canonical_json(d))
        print("\nwrote %s" % out)
        return
    for name, j in PILOT:
        d = compile_case(name, j)
        print(render(d))
        print()
        out = "case_dossier_%s.json" % d["case"]["tag"]
        with open(out, "w", encoding="utf-8") as fh:
            fh.write(canonical_json(d))
        print("wrote %s\n" % out)

if __name__ == "__main__":
    main(sys.argv)
