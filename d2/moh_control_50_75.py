#!/usr/bin/env python3
"""moh_control_50_75.py  (NEW 2026-07-27; read-only over all existing artifacts)

RESOLVING THE (50,75) TENSION:  our sec.8 witnesses DO survive there, and what
kills (50,75) is a condition our transferred machinery does not contain AT ALL.

THE TENSION
-----------
  half one   yplace_transfer.py sec.G (57/57):  at a class row's y-place the
             whole of PROOF_72_108 secs.3-7 transfers, the row lands at k = 0 /
             Cor 8.5, and FOUR EXPLICIT WITNESSES satisfy the entire sec.8.1
             k = 0 system with residual 0 inside every transferred cap.
  half two   moh_discards.py (21/21):  F_2(2,3)/75 IS the case (50,75), which
             GGV3 sec.5 kills outright in two gamma charts -- and (5,20) is a
             class corner (b0 = 4 a0).

So the witnesses must fail to lift, and something must kill them.

THE ANSWER, IN ONE SENTENCE
---------------------------
Every condition our transfer imposes is CLOSED -- an equation = 0, an order
floor, a degree cap.  GGV3's kill is an OPEN condition: corner primitivity
(a6)/(b6) requires one named window coefficient to be NONZERO, at a depth the
window equations force to vanish.  Our transferred sec.8 states no nonvanishing
requirement about any coefficient at any depth, so its solution set is strictly
larger than the germ locus and a point of it need not be a germ.  There is no
contradiction; the transfer is SOUND and INCOMPLETE, in a nameable way.

That is not a re-statement of yplace_transfer's own "NOT CLAIMED: the witnesses
are not germs".  It is the first EXTERNAL confirmation of which missing
ingredient does the work -- and it is exactly the one YPLACE_TRANSFER.md sec.8
opened as lead #2 ("leading-coefficient information at [y^21]B") and left
untouched.

WHAT IS DERIVED HERE VS REPLAYED
--------------------------------
f2_tower.a2_certificate() writes GGV3's gamma=2 system down as 13 LITERALS and
supplies a^3 = 2 as a GIVEN; it consumes no corner data (SESSION_HANDOFF's
"REPLAY TRAP").  Section C below does not do that.  It builds the Laurent series
Z from (a4)/(b4), computes E_k = (Z^2)_{-k} and (Z^3 + lam Z^{-1})_{-k} from the
definition, and reproduces GGV3's published E_1..E_8 (gamma=3) and E_1..E_13
(gamma=2) TERM FOR TERM -- then eliminates and derives, rather than transcribes,
    F_{-1} = -3 C_{-1}C_{-2},   3 C_0 C_{-1}^2 = 3 C_{-2}^2 + 2F_{-2} + 2 lam,
which are the paper's own two displayed conclusions.  The premises consumed are
(a1)-(a6) only.  GGV3 does not prove (a1)-(a6) ("We do not provide proofs for
this first part", tex:1716), so this is a REPRODUCTION of a published kill from
its stated premises, NOT an independent proof that (50,75) is dead.

SCOPE / EVIDENCE BOUNDARY -- read before citing
-----------------------------------------------
  PROVED         the E_k identities; the elimination; "a product of Laurent
                 polynomials is a unit iff both factors are"; the E-depth law
                 K + c(m-1) = L + c(n-1).
  EXACT-CHECKED  the sec.8.1 k=0 witnesses at (50,75)'s own caps; the reduced
                 polygons and their degrees; the gamma-admissibility tables.
  CITATION-LEVEL (a1)-(a6) and (b1)-(b6) themselves, and "Moh ruled (50,75)
                 out".  Both are GGV3/GGV5's word, not re-derived here.
  INFERRED       the shape of (a1) at (m,n) = (3,5) (P = C^3, Q = C^5 + ...).
                 One anchor only, at (2,3).  Flagged at its point of use.

Checker: --quiet, exit 0 iff every check passes.  ~20 s.  Exact sympy only.
"""
import json
import os
import re
import sys

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__)) or "."
sys.path.insert(0, HERE)
QUIET = "--quiet" in sys.argv
_ok = [0]
_fail = []
y = sp.Symbol("y")
expand = sp.expand
Rational = sp.Rational


def check(name, cond, detail=""):
    if cond:
        _ok[0] += 1
        if not QUIET:
            print("[OK]  %s" % name)
    else:
        _fail.append(name)
        print("[FAIL] %s %s" % (name, detail))


def head(msg):
    if not QUIET:
        print("\n" + "=" * 78 + "\n" + msg + "\n" + "=" * 78)


# GGV3's .tex is copyrighted and is NOT redistributed in the public release
# (PUBLICATION_AUDIT.md item 2), so every "does GGV3 say X" probe below answers
# from paper_src/upstream_quotes.json -- one sha256-pinned transcription -- and
# is re-derived from the source whenever the .tex is on disk (check Z1).
# Reading the .tex directly here used to make this checker die with
# FileNotFoundError on a clean public clone.
import upstream_quotes as uq

ATLAS = json.load(open(os.path.join(HERE, "corner_atlas.json"), encoding="utf-8"))
ROWS = {}
for _r in ATLAS["rows"]:
    ROWS.setdefault(_r["id"], []).append(_r)


# ===========================================================================
head("A.  IDENTITY -- the numbers describe the case the LABEL names")
# ===========================================================================
check("A1  GGV3 sec.5 says it verifies deg(P_0) = 50, deg(Q_0) = 75 -- so the "
      "label '(50,75)' is the paper's own (%s)" % uq.cite("ctl.deg_P0_50"),
      uq.present("ctl.deg_P0_50") and uq.present("ctl.deg_Q0_75"))

check("A2  GGV3 sec.5 says A_0 = (5,20) for that case (its cite: GGV1 Rmk 7.10) "
      "(%s)" % uq.cite("ctl.A0_5_20"), uq.present("ctl.A0_5_20"))

check("A3  GGV3 sec.5 says the reduced pair has [P_1,Q_1] = x^2, deg(P_1) = 10, "
      "deg(Q_1) = 15 (%s)" % uq.cite("ctl.bracket_P1Q1"),
      uq.present("ctl.bracket_P1Q1") and uq.present("ctl.deg_P1_10")
      and uq.present("ctl.deg_Q1_15"))

r5075 = ROWS["F_2(2,3)/75"][0]
r7550 = ROWS["F_3(3,2)/75"][0]
r75125 = ROWS["F_2(3,5)/125"][0]

check("A4  the atlas row F_2(2,3)/75 carries the SAME A_0 = (5,20) and max_deg "
      "= 75, and (m,n) = (2,3), so 50:75 = 2:3 -- the atlas row and GGV3's case "
      "are the same object",
      tuple(r5075["A0"]) == (5, 20) and r5075["max_deg"] == 75
      and (r5075["m"], r5075["n"]) == (2, 3))

check("A5  and (5,20) is a CLASS corner: b0 = 4*a0 exactly (20 = 4*5), the "
      "regime of the class of nine -- not (8,28)'s retraction shape "
      "b0 = 4*a0 - 4",
      20 == 4 * 5 and 28 == 4 * 8 - 4 and 28 != 4 * 8)

cd = None
import polygon_reduction as pr                                      # noqa: E402
cd = pr.corner_chart_data(5, 20)
check("A6  polygon_reduction's GUARDED corner_chart_data(5,20) gives t = 4, "
      "kappa = 2, C = y a monomial (deg C = ord C = 1), retraction FALSE -- "
      "i.e. exactly yplace_transfer's class-row chart (a,b,t,kappa) = (2,3,4,2)",
      (cd["t"], cd["kappa"], cd["deg_C"], cd["ord_C"], cd["monomial"],
       cd["retraction"]) == (4, 2, 1, 1, True, False))

red0 = pr.case_f2(0)
check("A7  polygon_reduction.case_f2(0) is itself LABELLED (50,75) "
      "(tag=%s, degs=%s) and carries A_0'=(1,0), reduced pair (2,3), l = 4, "
      "bracket x^2 -- so the in-repo reduction and GGV3 sec.5 agree on the "
      "chart, not merely on the corner"
      % (red0.tag, red0.signature["degs"]),
      red0.tag == "F2_j0_50_75" and red0.signature["degs"] == (50, 75)
      and red0.A0 == (5, 20) and red0.A0p == (1, 0) and red0.mn == (2, 3)
      and red0.l == 4 and red0.bracket == "x^2")

NP0 = red0.reduced["standard (proportional, Prop 8.2(1))"]["P"]
NQ0 = red0.reduced["standard (proportional, Prop 8.2(1))"]["Q"]
check("A8  CROSS-CHECK (two independent objects that must agree).  The COMPUTED "
      "reduced polygons give max_y N(P_1) = %d and max_y N(Q_1) = %d -- which "
      "ARE GGV3's published deg(P_1) = 10, deg(Q_1) = 15.  The reduction was "
      "not tuned to those numbers" % (max(p[1] for p in NP0),
                                      max(p[1] for p in NQ0)),
      max(p[1] for p in NP0) == 10 and max(p[1] for p in NQ0) == 15)

check("A9  and both agree with two independent closed forms: deg P_1 = "
      "m*(t+1) = 2*5 = 10 = deg P / a0 = 50/5, deg Q_1 = n*(t+1) = 15 = 75/5",
      2 * (4 + 1) == 10 == 50 // 5 and 3 * (4 + 1) == 15 == 75 // 5)

check("A10 F_3(3,2)/75 is the P<->Q SWAP of it: same A_0 = (5,20), same t = 4, "
      "kappa = 2, same max_deg 75, (m,n) = (3,2), and the atlas's own sorted "
      "reduced pair (D_P,D_Q) = (2,3) is IDENTICAL to F_2(2,3)/75's.  So every "
      "object yplace_transfer builds is bit-identical at the two rows",
      tuple(r7550["A0"]) == (5, 20) and (r7550["m"], r7550["n"]) == (3, 2)
      and r7550["max_deg"] == 75
      and (r7550["gates"]["G3"]["sub"]["N_Q"]["D_P"],
           r7550["gates"]["G3"]["sub"]["N_Q"]["D_Q"])
      == (r5075["gates"]["G3"]["sub"]["N_Q"]["D_P"],
          r5075["gates"]["G3"]["sub"]["N_Q"]["D_Q"]) == (2, 3))

check("A11 both /75 rows are among the atlas's moh_discarded list AND both are "
      "red_in_paper -- they are settled in the literature (CITATION-LEVEL, per "
      "moh_discards.py G1), which is exactly what makes them a CONTROL",
      {"F_2(2,3)/75", "F_3(3,2)/75"} <= set(ATLAS["moh_discarded"])
      and all(r["provenance"].get("red_in_paper") for r in (r5075, r7550)))

CLASS8 = [r for r in ATLAS["rows"] if r["A0"][1] == 4 * r["A0"][0]
          and (r["gates"]["G3"]["sub"]["N_Q"].get("D_P"),
               r["gates"]["G3"]["sub"]["N_Q"].get("D_Q")) == (2, 3)]
check("A12 the class of nine splits 8 + 1 exactly as yplace_transfer scopes it: "
      "EIGHT rows with reduced pair (2,3) (of which the two /75 rows are two) "
      "and ONE with (3,5), namely F_2(3,5)/125 = (75,125).  So (50,75) is "
      "INSIDE yplace_transfer's stated scope and (75,125) is outside it",
      len(CLASS8) == 8
      and sum(1 for r in ATLAS["rows"] if r["A0"][1] == 4 * r["A0"][0]) == 9
      and (r75125["gates"]["G3"]["sub"]["N_Q"]["D_P"],
           r75125["gates"]["G3"]["sub"]["N_Q"]["D_Q"]) == (3, 5)
      and "none of this applies" in open(
          os.path.join(HERE, "YPLACE_TRANSFER.md"), encoding="utf-8").read())


# ===========================================================================
head("B.  THE sec.8 WITNESSES, INSTANTIATED AT (50,75) -- they SURVIVE")
# ===========================================================================
# The caps are not generic: they are corner (5,20)'s own, recomputed here from
# the COMPUTED reduced polygon rather than imported.


def hull_slopes(verts):
    """(sigma, tau) = the slopes of the lower- resp. upper-hull EDGE terminating
    at the maximal-i vertex of N(P) -- PROOF sec.2.6(i)'s direction functionals.
    Same construction as yplace_transfer.hull_slopes (this is a deliberate
    re-implementation, so agreement with sec.D there is a cross-check)."""
    V0 = sorted(set(map(tuple, verts)))
    lowpts = sorted({i: min(j for ii, j in V0 if ii == i) for i, _ in V0}.items())
    uppts = sorted({i: max(j for ii, j in V0 if ii == i) for i, _ in V0}.items())

    def chain(V, sign):
        H = []
        for p in V:
            while len(H) >= 2:
                cr = ((H[-1][0] - H[-2][0]) * (p[1] - H[-2][1])
                      - (H[-1][1] - H[-2][1]) * (p[0] - H[-2][0]))
                if sign * cr <= 0:
                    H.pop()
                else:
                    break
            H.append(p)
        return H

    lo, up = chain(lowpts, +1), chain(uppts, -1)
    return (Rational(lo[-1][1] - lo[-2][1], lo[-1][0] - lo[-2][0]),
            Rational(up[-1][1] - up[-2][1], up[-1][0] - up[-2][0]))


r828 = pr.case_8_28().reduced
s1 = hull_slopes(r828["sub1 (case c)"]["P"])
s2 = hull_slopes(r828["sub2 (cases a,b)"]["P"])
check("B0  CALIBRATION of the slope reader on the PUBLISHED control.  Fed "
      "(8,28)'s computed reduced hulls it returns (sigma,tau) = %s (sub1) and "
      "%s (sub2), so the a=2 induction ord_y h_k >= k(a*q - sigma) = 12k and "
      "deg_y h_k <= k(a*degC - tau) = 15k / 14k -- which ARE PROOF sec.2.6's "
      "published ord D_{j_x} >= 48 - 12 j_x and lambda = 3 / 2" % (s1, s2),
      s1 == (2, 1) and s2 == (2, 2)
      and (2 * 7 - s1[0], 2 * 8 - s1[1]) == (12, 15)
      and (2 * 7 - s2[0], 2 * 8 - s2[1]) == (12, 14))

sig, tau = hull_slopes(NP0)
check("B1  corner (5,20)'s COMPUTED reduced polygon N(P_1) = %s has terminal "
      "hull slopes (sigma, tau) = (%s, %s), giving PROOF sec.2.6's affine "
      "induction ord_y h_k >= k*(a*q - sigma) = k and deg_y h_k <= "
      "k*(a*deg C - tau) = 3k, so cap-lambda = 3 - 1 = 2" % (NP0, sig, tau),
      (sig, tau) == (1, -1) and 2 * 1 - sig == 1 and 2 * 1 - tau == 3)

check("B2  which is what the shipped atlas records for THIS row (G3.lam = 2, "
      "verdict PASS) -- two independent readings of one number",
      r5075["gates"]["G3"]["sub"]["lam"]["lam"] == "2"
      and r5075["gates"]["G3"]["sub"]["lam"]["verdict"] == "PASS")

# --- PROOF (8.1.1), retyped from the paper, at k = 0 (Pi = 1), place pi = y ---
gam, zeta = sp.symbols("gamma zeta", nonzero=True)


def k0_system(A, z, zeta_v, gam_v, Qq=sp.Integer(1), c=Rational(1, 2),
              pi=None):
    """PROOF_72_108 (8.1.1) with Pi = 1, uniformiser pi, residual Q_Pi = Qq.

    Written from the paper's own display:
        g1 = 1/2 g^2 d1 P^2 + g P (d2 A + C) + A B
        g2 = d2 A^2 + 2 A C + B^2 - g^2 d0 P^2
        g3 = -g d0 P A - 1/2 d1 A^2 + B C - 1/6 g^3 pi^9 P^3
        box = 3 A^2 + g^2 d2 P^2 + 3 g P B - mu pi^3 Q_Pi
        u := g d2,  w := 1/2 g^2 d1 P,  F := A(u+2v)+w,  Z := A^2 - g P^2 v
        (*)  F Z = 1/6 g^5 pi^9 P^4
    At k = 0: Pi = 1 and B = Pi v = v.  Given (A, z, zeta, gamma) the boxed row
    determines u, Z determines v, (*) determines w, g1 determines C, g2 gives d0.
    """
    if pi is None:
        pi = y
    mu = 2 * c / gam_v
    uu = expand((mu * pi**3 * Qq - 6 * A**2 + 3 * zeta_v * pi**z) / gam_v)
    d2 = expand(uu / gam_v)
    v = expand((A**2 - zeta_v * pi**z) / gam_v)
    Z = expand(A**2 - gam_v * v)
    F = expand(Rational(1, 6) * gam_v**5 * pi**9 / Z)
    w = expand(F - A * (uu + 2 * v))
    d1 = expand(2 * w / gam_v**2)
    CT = expand(-(A * (uu + v) + w) / gam_v)
    d0 = expand((d2 * A**2 + 2 * A * CT + v**2) / gam_v**2)
    g1 = expand(Rational(1, 2) * gam_v**2 * d1 + gam_v * (d2 * A + CT) + A * v)
    g2 = expand(d2 * A**2 + 2 * A * CT + v**2 - gam_v**2 * d0)
    g3 = expand(-gam_v * d0 * A - Rational(1, 2) * d1 * A**2 + v * CT
                - Rational(1, 6) * gam_v**3 * pi**9)
    box = expand(3 * A**2 + gam_v**2 * d2 + 3 * gam_v * v - mu * pi**3 * Qq)
    FZ = expand(F * Z - Rational(1, 6) * gam_v**5 * pi**9)
    return dict(u=uu, v=v, w=w, d0=d0, d1=d1, d2=d2, CT=CT, Z=Z, F=F,
                g1=sp.simplify(g1), g2=sp.simplify(g2), g3=sp.simplify(g3),
                box=sp.simplify(box), FZ=sp.simplify(FZ))


def ordy(p):
    p = expand(p)
    if p == 0:
        return sp.oo
    return min(m[0] for m in sp.Poly(p, y).monoms())


# The caps ARE the sec.B1 slopes: deg <= 3w and ord >= w on the ledger
# (d2,d1,d0,e,R,S,T) at weights w = (2,3,4,5,6,7,8) stripped by y^9.
CAPS = dict(A=9, u=6, v=12, w=9, d0=12, CT=15)
ORDS = dict(A=1, u=2, v=2, w=3, CT=3)
check("B3  the six degree caps (A,u,v,w,d0,C) <= (9,6,12,9,12,15) ARE corner "
      "(5,20)'s: the unstripped ledger deg(d2,d1,d0,e,R,S,T) <= 3w = "
      "(6,9,12,15,18,21,24) is the slope-3 cap of B1 applied weight by weight, "
      "and stripping y^9 off (R,S,T) gives (A,v,C) <= (9,12,15)",
      [3 * w for w in (2, 3, 4, 5, 6, 7, 8)] == [6, 9, 12, 15, 18, 21, 24]
      and (18 - 9, 21 - 9, 24 - 9) == (CAPS["A"], CAPS["v"], CAPS["CT"])
      and CAPS["u"] == 6 and CAPS["w"] == 9 and CAPS["d0"] == 12)

WITS = [(y, 2, sp.Integer(1), sp.Integer(1)),
        (y, 3, sp.Integer(1), sp.Integer(1)),
        (y**2, 4, Rational(1, 3), sp.Integer(1)),
        (sp.Integer(0), 5, sp.Integer(2), sp.Integer(1))]
allw, why = True, []
for A, z, zv, gv in WITS:
    r = k0_system(A, z, zv, gv)
    zero = all(r[k] == 0 for k in ("g1", "g2", "g3", "box", "FZ"))
    capok = all(r[k] == 0 or sp.degree(r[k], y) <= CAPS[k] for k in CAPS
                if k != "A")
    ordok = all(r[k] == 0 or ordy(r[k]) >= ORDS[k] for k in ORDS if k != "A")
    aok = (A == 0) or (ordy(A) >= 1 and sp.degree(A, y) <= CAPS["A"])
    if not (zero and capok and ordok and aok):
        why.append((A, z, zero, capok, ordok, aok))
    allw &= zero and capok and ordok and aok
check("B4  ALL FOUR WITNESSES SURVIVE AT (50,75).  Rebuilt here from "
      "PROOF (8.1.1) as printed, not imported from yplace_transfer: for "
      "(A,z,zeta,gamma) = (y,2,1,1), (y,3,1,1), (y^2,4,1/3,1), (0,5,2,1) the "
      "full k=0 system g1 = g2 = g3 = box = 0 AND (*) FZ = (1/6)gamma^5 y^9 "
      "holds with residual EXACTLY 0, inside all six caps and all five orders",
      allw, why)

check("B5  and the caps they are tested against are the TIGHTEST in the class: "
      "yplace_transfer sec.D5 records deg slopes 15/4, 4, 17/4 at the other "
      "three class corners (8,32), (9,36), (10,40), all > 3.  So (50,75) is "
      "where the transferred machinery is MOST constrained, and the witnesses "
      "survive there",
      all(s > 3 for s in (Rational(15, 4), Rational(4), Rational(17, 4))))

check("B6  nothing in the transferred chain reads (deg P, deg Q) at all -- the "
      "witnesses depend only on (a,b,t,kappa,C) = (2,3,4,2,y), ord_y Phi = 30 "
      "and the (5,20) caps.  So they are the SAME four points at F_2(2,3)/75 "
      "and at its P<->Q swap F_3(3,2)/75, and no asymmetry can appear between "
      "the two rows",
      k0_system(*WITS[1])["Z"] == k0_system(*WITS[1])["Z"]
      and (r5075["gates"]["G1"]["t"], r5075["gates"]["G1"]["kappa"])
      == (r7550["gates"]["G1"]["t"], r7550["gates"]["G1"]["kappa"]))

# --- the FULL nonvanishing inventory of PROOF sec.8 at k = 0 -------------------
NONVANISHING_SEC8 = ["gamma != 0 (leading coeff of e, a scalar)  [sec.8.4]",
                     "mu != 0                                    [sec.8.5]",
                     "zeta != 0                                  [sec.8.4.1]",
                     "F != 0 and Z != 0                          [sec.8.3 (*)]",
                     "A(r) != 0 at every root r of Pi            [sec.8.2] "
                     "-- VACUOUS at k = 0, Pi = 1"]
nv_ok = True
for A, z, zv, gv in WITS:
    r = k0_system(A, z, zv, gv)
    nv_ok &= (gv != 0) and (zv != 0) and (r["F"] != 0) and (r["Z"] != 0)
check("B7  THE INVENTORY.  PROOF sec.8 at k = 0 states exactly FIVE "
      "nonvanishing conditions -- gamma, mu, zeta, F, Z -- and the fifth "
      "(A(r) != 0 at roots of Pi) is VACUOUS because Pi = 1.  All four "
      "witnesses satisfy every one of them.  There is no sixth: sec.8 never "
      "requires a NAMED COEFFICIENT at a NAMED depth to be nonzero",
      nv_ok and len(NONVANISHING_SEC8) == 5)


# ===========================================================================
head("C.  THE KILLER, DERIVED FROM (a1)-(a6) -- not replayed from literals")
# ===========================================================================
u_ = sp.Symbol("u"); lam = sp.Symbol("lam")


def E_system(c_lead, extra_pos, m, n, K, L, trunc=None):
    """E_k = (Z^m)_{-k}, k=1..K ;  E_{K+k} = (Z^n + lam Z^{m-n})_{-k}, k=1..L,
    for Z = x^c + (sum over extra_pos positive slots) + Z_0 + Z_{-1}x^{-1}+...

    Built from the DEFINITION.  Returns (E, Zs)."""
    c = c_lead
    if trunc is None:
        trunc = max(K + c * m, L + c * n) + 1
    depth = K + c * (m - 1) + 2
    Zs = {}
    for j in extra_pos:
        Zs[j] = sp.Symbol("Z%d" % j)
    Zs[0] = sp.Symbol("Z0")
    for k in range(1, depth + 1):
        Zs[-k] = sp.Symbol("Zm%d" % k)

    def tr(e):
        e = expand(e)
        if e == 0:
            return sp.Integer(0)
        p = sp.Poly(e, u_)
        return expand(sum(p.nth(i) * u_**i
                          for i in range(0, min(p.degree(), trunc) + 1)))

    # Zu := u^c * Z = 1 + sum_j Zs[j] u^{c-j}
    Zu = tr(1 + sum(Zs[j] * u_**(c - j) for j in sorted(Zs)))
    powm = sp.Integer(1)
    for _ in range(m):
        powm = tr(powm * Zu)
    pown = sp.Integer(1)
    for _ in range(n):
        pown = tr(pown * Zu)
    h = tr(Zu - 1)
    inv = sp.Integer(0)
    t = sp.Integer(1)
    for i in range(trunc + 2):
        inv = tr(inv + (-1)**i * t)
        t = tr(t * h)
    assert tr(inv * Zu) == 1
    e_mn = m - n                                  # negative
    powmn = sp.Integer(1)
    for _ in range(-e_mn):
        powmn = tr(powmn * inv)

    def co(e, k):
        if k < 0:
            return sp.Integer(0)
        e = expand(e)
        return sp.Integer(0) if e == 0 else expand(sp.Poly(e, u_).nth(k))

    E = {}
    for k in range(1, K + 1):
        E[k] = co(powm, k + c * m)                     # Z^m = u^{-cm} Zu^m
    for k in range(1, L + 1):
        E[K + k] = expand(co(pown, k + c * n)
                          + lam * co(powmn, k + c * e_mn))
    return E, Zs


# --- gamma = 3 chart: (a4) C = x^2 + C_0 + C_{-1}x^{-1} + ... ; (a1) P=C^2 ----
E3, Z3s = E_system(c_lead=2, extra_pos=[], m=2, n=3, K=5, L=3)
Z0, Zm1, Zm2, Zm3, Zm4, Zm5, Zm6, Zm7 = (Z3s[0], Z3s[-1], Z3s[-2], Z3s[-3],
                                         Z3s[-4], Z3s[-5], Z3s[-6], Z3s[-7])
PUB3 = {  # GGV3 tex:1861-1873, transcribed
    1: 2 * Z0 * Zm1 + 2 * Zm3,
    2: Zm1**2 + 2 * Z0 * Zm2 + 2 * Zm4,
    3: 2 * Zm1 * Zm2 + 2 * Z0 * Zm3 + 2 * Zm5,
    4: Zm2**2 + 2 * Zm1 * Zm3 + 2 * Z0 * Zm4 + 2 * Zm6,
    5: 2 * Zm2 * Zm3 + 2 * Zm1 * Zm4 + 2 * Z0 * Zm5 + 2 * Zm7,
    6: 3 * Z0**2 * Zm1 + 6 * Zm1 * Zm2 + 6 * Z0 * Zm3 + 3 * Zm5,
    7: (lam + 3 * Z0 * Zm1**2 + 3 * Z0**2 * Zm2 + 3 * Zm2**2 + 6 * Zm1 * Zm3
        + 6 * Z0 * Zm4 + 3 * Zm6),
    8: (Zm1**3 + 6 * Z0 * Zm1 * Zm2 + 3 * Z0**2 * Zm3 + 6 * Zm2 * Zm3
        + 6 * Zm1 * Zm4 + 6 * Z0 * Zm5 + 3 * Zm7),
}
check("C1  DERIVED, NOT TRANSCRIBED.  Forming E_k = (Z^2)_{-k} (k=1..5) and "
      "E_{5+k} = (Z^3 + lam Z^{-1})_{-k} (k=1..3) from Z = x^2 + Z_0 + "
      "Z_{-1}x^{-1}+... reproduces GGV3's published E_1..E_8 TERM FOR TERM.  "
      "This validates the series conventions before anything is eliminated",
      all(expand(E3[k] - PUB3[k]) == 0 for k in range(1, 9)),
      [k for k in range(1, 9) if expand(E3[k] - PUB3[k]) != 0])

check("C2  and the paper's own remark is reproduced: Z_{-7} is the DEEPEST "
      "coefficient occurring, and it occurs in E_5 and E_8 only",
      all(E3[k].has(Zm7) == (k in (5, 8)) for k in range(1, 9))
      and not any(E3[k].has(Z3s[-8]) for k in range(1, 9)))

sol3 = sp.solve([E3[k] for k in range(1, 6)],
                [Zm3, Zm4, Zm5, Zm6, Zm7], dict=True)[0]
E6s = sp.simplify(E3[6].subs(sol3))
E7s = sp.simplify(expand(E3[7].subs(sol3)))
check("C3  DERIVED: imposing E_1 = ... = E_5 = 0 collapses E_6 to exactly "
      "3*C_{-1}*C_{-2}.  With (a3)'s E_6 = -F_{-1} this IS the paper's "
      "displayed F_{-1} = -3 C_{-1} C_{-2}",
      expand(E6s - 3 * Zm1 * Zm2) == 0)

Fm2 = sp.Symbol("Fm2")
rel = expand(-2 * (E7s + Fm2))
check("C4  DERIVED: E_7 = -F_{-2} then reads 3 C_0 C_{-1}^2 - 3 C_{-2}^2 "
      "- 2 lam - 2 F_{-2} = 0.  Note this is STRONGER than GGV3's printed "
      "'C_0(3C_0C_{-1}^2 - 3C_{-2}^2 - 2 lam) = 2 C_0 F_{-2}': the spurious "
      "C_0 factor (and hence the 'either C_0 = 0 or ...' branch) is an "
      "artefact of their elimination, not of the system.  Multiplying by C_0 "
      "recovers the printed form exactly",
      expand(rel - (3 * Z0 * Zm1**2 - 3 * Zm2**2 - 2 * lam - 2 * Fm2)) == 0
      and expand(Z0 * rel
                 - (Z0 * (3 * Z0 * Zm1**2 - 3 * Zm2**2 - 2 * lam)
                    - 2 * Z0 * Fm2)) == 0)

check("C5  (a3) supplies F_{-1} = y^7, a UNIT of K[y,y^-1].  Since "
      "C_{-1}C_{-2} = -y^7/3 is a unit and K[y,y^-1] has unit group {c y^n}, "
      "BOTH factors are monomials -- this is where the kill is generated, and "
      "it needs no cap yet",
      uq.present("ctl.F_minus1_y7"))

check("C6  ZERO MARGIN.  (a5) is deg_y(C_{-k}) <= k+2, so deg_y C_{-1} <= 3 "
      "and deg_y C_{-2} <= 4.  Two monomials whose exponents sum to 7 with "
      "e1 <= 3, e2 <= 4 force e1 = 3, e2 = 4 EXACTLY -- one integer of slack "
      "in neither place.  Hence C_{-1} = a y^3, C_{-2} = b y^4, a,b != 0",
      uq.present("ctl.deg_y_C_cap")
      and [(e1, 7 - e1) for e1 in range(0, 4) if 7 - e1 <= 4] == [(3, 4)])

aa, bb, f2, f4, f6, f8, lm = sp.symbols("aa bb f2 f4 f6 f8 lm")
C0 = expand((3 * (bb * y**4)**2 + 2 * (f8 * y**8 + f6 * y**6 + f4 * y**4
                                       + f2 * y**2) + 2 * lm)
            / (3 * (aa * y**3)**2))
orders = sorted(set(sp.degree(sp.numer(sp.together(tm)), y)
                    - sp.degree(sp.denom(sp.together(tm)), y)
                    for tm in sp.Add.make_args(C0)))
check("C7  the forced C_0 = (3C_{-2}^2 + 2F_{-2} + 2 lam)/(3C_{-1}^2) has "
      "y-support %s -- FORCED FLOOR j_min(0) = -6, reached only by the lam "
      "term, and this reproduces GGV3's displayed C_0 term for term"
      % orders,
      min(orders) == -6 and orders == [-6, -4, -2, 0, 2]
      and expand(C0 - (2 * lm / (3 * aa**2 * y**6) + 2 * f2 / (3 * aa**2 * y**4)
                       + 2 * f4 / (3 * aa**2 * y**2) + 2 * f6 / (3 * aa**2)
                       + bb**2 * y**2 / aa**2
                       + 2 * f8 * y**2 / (3 * aa**2))) == 0)

check("C8  THE KILL.  (a6) declares C_0 = c_{0,2}y^2 + ... + c_{0,-10}y^{-10} "
      "with c_{0,-10} != 0 -- a REQUIRED-NONZERO at depth -10.  -10 < -6 = the "
      "forced floor, so c_{0,-10} = 0.  Contradiction.  Margin = 4 = two "
      "chart-steps (the gamma=3 chart is y -> y^-2, so C_0 lives on even "
      "y-orders)",
      uq.present("ctl.c0_minus10_ne0") and -10 < -6 and (-6) - (-10) == 4
      and all(o % 2 == 0 for o in orders))

check("C9  MUTATION on C6: the zero margin is load-bearing.  Had (a5) read "
      "k+3, the caps would be 4 and 5, 4+5 = 9 > 7, and (3,4),(4,3) would both "
      "survive -- no forced C_{-1} = a y^3, no forced floor, no kill",
      len([(e1, 7 - e1) for e1 in range(0, 5) if 7 - e1 <= 5]) == 3)

check("C10 MUTATION on C8: the kill predicate is not a tautology.  Run with the "
      "SAME required-nonzero (0,-10) against a hypothetical floor of -12 it "
      "does NOT fire; run with required (0,-4) against the true floor -6 it "
      "does NOT fire.  It fires exactly when required-depth < forced floor",
      not (-10 < -12) and not (-4 < -6) and (-10 < -6))

# --- the gamma = 2 chart: same schema, different slot -------------------------
E2s_, Z2s = E_system(c_lead=3, extra_pos=[1], m=2, n=3, K=8, L=5)
W1, W0 = Z2s[1], Z2s[0]
Wm = {k: Z2s[-k] for k in range(1, 12)}
PUB2 = {  # GGV3 tex:1905-1932, transcribed
    1: 2 * W0 * Wm[1] + 2 * W1 * Wm[2] + 2 * Wm[4],
    2: Wm[1]**2 + 2 * W0 * Wm[2] + 2 * W1 * Wm[3] + 2 * Wm[5],
    3: 2 * Wm[1] * Wm[2] + 2 * W0 * Wm[3] + 2 * W1 * Wm[4] + 2 * Wm[6],
    4: (Wm[2]**2 + 2 * Wm[1] * Wm[3] + 2 * W0 * Wm[4] + 2 * W1 * Wm[5]
        + 2 * Wm[7]),
    5: (2 * Wm[2] * Wm[3] + 2 * Wm[1] * Wm[4] + 2 * W0 * Wm[5]
        + 2 * W1 * Wm[6] + 2 * Wm[8]),
    6: (Wm[3]**2 + 2 * Wm[2] * Wm[4] + 2 * Wm[1] * Wm[5] + 2 * W0 * Wm[6]
        + 2 * W1 * Wm[7] + 2 * Wm[9]),
    7: (2 * Wm[3] * Wm[4] + 2 * Wm[2] * Wm[5] + 2 * Wm[1] * Wm[6]
        + 2 * W0 * Wm[7] + 2 * W1 * Wm[8] + 2 * Wm[10]),
    8: (Wm[4]**2 + 2 * Wm[3] * Wm[5] + 2 * Wm[2] * Wm[6] + 2 * Wm[1] * Wm[7]
        + 2 * W0 * Wm[8] + 2 * W1 * Wm[9] + 2 * Wm[11]),
    9: (3 * W0**2 * Wm[1] + 3 * W1 * Wm[1]**2 + 6 * W0 * W1 * Wm[2]
        + 3 * Wm[2]**2 + 3 * W1**2 * Wm[3] + 6 * Wm[1] * Wm[3]
        + 6 * W0 * Wm[4] + 6 * W1 * Wm[5] + 3 * Wm[7]),
    10: (3 * W0 * Wm[1]**2 + 3 * W0**2 * Wm[2] + 6 * W1 * Wm[1] * Wm[2]
         + 6 * W0 * W1 * Wm[3] + 6 * Wm[2] * Wm[3] + 3 * W1**2 * Wm[4]
         + 6 * Wm[1] * Wm[4] + 6 * W0 * Wm[5] + 6 * W1 * Wm[6] + 3 * Wm[8]),
    11: (lam + Wm[1]**3 + 6 * W0 * Wm[1] * Wm[2] + 3 * W1 * Wm[2]**2
         + 3 * W0**2 * Wm[3] + 6 * W1 * Wm[1] * Wm[3] + 3 * Wm[3]**2
         + 6 * W0 * W1 * Wm[4] + 6 * Wm[2] * Wm[4] + 3 * W1**2 * Wm[5]
         + 6 * Wm[1] * Wm[5] + 6 * W0 * Wm[6] + 6 * W1 * Wm[7] + 3 * Wm[9]),
    12: (3 * Wm[1]**2 * Wm[2] + 3 * W0 * Wm[2]**2 + 6 * W0 * Wm[1] * Wm[3]
         + 6 * W1 * Wm[2] * Wm[3] + 3 * W0**2 * Wm[4] + 6 * W1 * Wm[1] * Wm[4]
         + 6 * Wm[3] * Wm[4] + 6 * W0 * W1 * Wm[5] + 6 * Wm[2] * Wm[5]
         + 3 * W1**2 * Wm[6] + 6 * Wm[1] * Wm[6] + 6 * W0 * Wm[7]
         + 3 * Wm[10] + 6 * W1 * Wm[8]),
    13: (-lam * W1 + 3 * Wm[1] * Wm[2]**2 + 3 * Wm[1]**2 * Wm[3]
         + 6 * W0 * Wm[2] * Wm[3] + 3 * W1 * Wm[3]**2 + 6 * W0 * Wm[1] * Wm[4]
         + 6 * W1 * Wm[2] * Wm[4] + 3 * Wm[4]**2 + 3 * W0**2 * Wm[5]
         + 6 * W1 * Wm[1] * Wm[5] + 6 * Wm[3] * Wm[5] + 6 * W0 * W1 * Wm[6]
         + 6 * Wm[2] * Wm[6] + 3 * W1**2 * Wm[7] + 6 * Wm[1] * Wm[7]
         + 6 * W0 * Wm[8] + 6 * W1 * Wm[9] + 3 * Wm[11]),
}
check("C11 THE SAME SCHEMA IN THE SIBLING CHART.  From (b4) Z = x^3 + Z_1 x + "
      "Z_0 + ... the same construction reproduces GGV3's published gamma=2 "
      "list E_1..E_13 TERM FOR TERM, including the -lam Z_1 of E_13 and the "
      "deepest coefficient Z_{-11} appearing in E_8 and E_13 only",
      all(expand(E2s_[k] - PUB2[k]) == 0 for k in range(1, 14))
      and all(E2s_[k].has(Wm[11]) == (k in (8, 13)) for k in range(1, 14)),
      [k for k in range(1, 14) if expand(E2s_[k] - PUB2[k]) != 0])

check("C12 and its required-nonzero sits at a DIFFERENT slot: (b6) e_{-10} != 0 "
      "if C_0 = 0, and (b5) c_{-1,1} != 0.  Same contract, three different "
      "(series, y_order) slots across two charts -- so 'required-nonzero vs "
      "forced floor' is the chart-general SCHEMA, and (0,-10) is the (50,75) "
      "gamma=3 INSTANCE of it",
      uq.present("ctl.e_minus10_ne0") and uq.present("ctl.c_minus1_1_ne0"))


# ===========================================================================
head("D.  WHY OUR TRANSFER MISSES IT -- closed conditions cannot state (a6)")
# ===========================================================================
check("D1  every condition the transferred sec.8 imposes is CLOSED: five "
      "polynomial equations (g1,g2,g3,box,(*)) and eleven inequalities "
      "(six deg caps, five ord floors).  A degree cap and an order floor are "
      "both satisfied by MORE vanishing, never less -- so no combination of "
      "them can require a coefficient to be nonzero",
      len(("g1", "g2", "g3", "box", "FZ")) == 5
      and len(CAPS) + len(ORDS) == 11)

# The contract kill predicate (ENDPOINT_CONTRACT.md sec.2), applied to both sides.
def kill_predicate(required_nonzero, forced_floor):
    return [(s, j) for (s, j) in required_nonzero
            if s in forced_floor and j < forced_floor[s]]


THEIRS = kill_predicate([(-1, 3), (-2, 4), (0, -10)], {0: -6, -1: 3, -2: 4})
OURS = kill_predicate([], {0: 9})          # our transferred class-row data
check("D2  THE DECISIVE COMPARISON, one predicate, two inputs.  With GGV3's "
      "contract the ENDPOINT_CONTRACT kill predicate returns %s.  With OUR "
      "transferred class-row data it returns %s -- because our required_nonzero "
      "list is EMPTY.  The predicate is identical; the missing datum is the "
      "primitivity list" % (THEIRS, OURS),
      THEIRS == [(0, -10)] and OURS == [])

# Our side's strongest coefficient statement, made into a contract and tested.
# yplace_transfer sec.E6: e = gamma*y^9 and B = 2*Phi/e = y^21/gamma, both
# FORCED MONOMIALS.  Grant our side the two nonzero facts it does have.
OUR_FLOOR = {"e": 9, "B": 21}
OUR_REQ = [("e", 9), ("B", 21)]                 # the ONLY nonzero coefficients
check("D4  our required_nonzero list cannot fire, and that is structural rather "
      "than an oversight.  Grant the transfer its two strongest coefficient "
      "facts (yplace_transfer sec.E6: e = gamma*y^9 and B = y^21/gamma are "
      "FORCED MONOMIALS).  A forced monomial marks every slot ABOVE the bottom "
      "as forbidden and the bottom one as nonzero -- so its nonzero slot sits "
      "EXACTLY AT the forced floor, never below it, and the predicate returns "
      "%s.  A required-nonzero can only fire when it is DEEPER than the floor, "
      "which is a fact no monomial statement can produce"
      % kill_predicate(OUR_REQ, OUR_FLOOR),
      kill_predicate(OUR_REQ, OUR_FLOOR) == []
      and all(j == OUR_FLOOR[s] for s, j in OUR_REQ))

YP = open(os.path.join(HERE, "yplace_transfer.py"), encoding="utf-8").read()
hits = [w for w in ("(a6)", "c_{0,-10}", "e_{-10}", "gamma=3 chart",
                    "required-nonzero", "forced floor") if w in YP]
check("D5  MECHANICAL CONFIRMATION (the repo's 'read the module, not its "
      "output' discipline).  yplace_transfer.py references (a6), c_{0,-10}, "
      "e_{-10}, 'required-nonzero' and 'forced floor' ZERO times each: %s.  "
      "Its ONE use of the word 'primitive' is check A3, 'R = x^t C is "
      "primitive ... it is no d-th power for d >= 2' -- a statement about R "
      "being non-a-power, NOT a statement that a named window coefficient at a "
      "named depth is nonzero.  Different notion, and the only one we carry"
      % (hits or "none"),
      hits == [] and YP.count("primitive") == 1
      and "no d-th power" in YP)

check("D6  hence NO CONTRADICTION.  The witnesses are points of a strictly "
      "weaker condition set than 'is a germ'.  yplace_transfer said exactly "
      "this ('the sec.8 witnesses are points of the sec.8.1 reduced system "
      "with the caps, not germs') -- what (50,75) adds is WHICH missing "
      "ingredient does the work, and it is the same one PROOF sec.7.4(c) names "
      "(leading-coefficient non-vanishing) and YPLACE_TRANSFER sec.8 lead #2 "
      "opens and leaves untouched",
      "not germs" in open(os.path.join(HERE, "YPLACE_TRANSFER.md"),
                          encoding="utf-8").read())


# ===========================================================================
head("E.  MUTATION CONTROLS -- the killer must NOT kill (72,108)")
# ===========================================================================
tt = y + 1
q_quartic = None
survive108 = []
for dA in [None] + list(range(0, 40)):
    for z in range(2, 7):
        degs = [7] + ([2 * dA] if dA is not None else []) + [z]
        mx = max(degs)
        if mx <= 6 or degs.count(mx) > 1:
            survive108.append((dA, z))
check("E1  MUTATION CONTROL #1 (the mandatory one).  (72,108) is killed at "
      "Cor 8.5 by a pure DEGREE count: deg(mu t^3 q) = 3 + deg q = 7 > 6 = the "
      "deg u cap, uniquely attained, for every (deg A, z) with z in [2,6].  "
      "That kill consumes NO primitivity input whatsoever, so removing (a6) -- "
      "which is what our transfer does -- leaves (72,108) dead.  The diagnosis "
      "does not damage the closed case",
      3 + 4 == 7 > 6 and survive108 == [])

check("E2  MUTATION CONTROL #2.  Conversely the diagnosis is not vacuous: it "
      "does NOT predict that every case survives.  Feeding a class row's "
      "Q_Pi = 1 (deg q = 0) into the same count gives 3 <= 6 and 25 surviving "
      "(deg A, z) pairs, while (72,108)'s q gives 0.  So the same closed "
      "machinery separates the two cases already",
      len([(dA, z) for dA in [None] + list(range(0, 40)) for z in range(2, 7)
           if max([3] + ([2 * dA] if dA is not None else []) + [z]) <= 6])
      == 25)

import gamma_from_corner as GFC                                     # noqa: E402


def gammas_at(a0, b0):
    out = {}
    for rec in GFC.analyse(a0, b0):
        if rec.get("rejected") is None and "gamma_admissible" in rec:
            out[tuple(rec["f"])] = rec["gamma_admissible"]
    return out


g520, g828 = gammas_at(5, 20), gammas_at(8, 28)
check("E3  MUTATION CONTROL #3 -- the gamma-chart ROUTE is not even runnable at "
      "(72,108).  gamma_from_corner (43 checks, calibrated on 28 published GGV1 "
      "data points) gives at (5,20) ONE surviving branch f=(4,16) with "
      "A_0'=(1,0) and gamma in %s, but at (8,28) ONE branch f=(6,21) with "
      "gamma in %s -- SIX charts, not three.  GGV3's two-chart argument has no "
      "(8,28) counterpart to be spuriously fired"
      % (list(g520.values())[0] if g520 else None,
         list(g828.values())[0] if g828 else None),
      list(g520.values()) == [[2, 3, 4]]
      and len(list(g828.values())[0]) == 6)

check("E4  MUTATION CONTROL #4 -- the class corners are NOT interchangeable, so "
      "the (50,75) kill must not be assumed to spread across the eight.  "
      "(8,32): NO branch survives the corner conditions at all.  (9,36): "
      "THREE surviving branches (A_0' = (4,1), (2,1), (1,0)).  (10,40): one "
      "branch with SIX admissible gammas.  Only (5,20) has the unique "
      "A_0'=(1,0) / gamma in {2,3,4} shape GGV3 sec.5 analyses",
      gammas_at(8, 32) == {}
      and len(gammas_at(9, 36)) == 3
      and sorted(gammas_at(10, 40).values())[0] == [2, 4, 5, 6, 7, 8])

check("E5  HONEST NEGATIVE, stated so it is not inferred from absence.  Nothing "
      "here shows the (50,75) witnesses fail a test we can RUN.  The witnesses "
      "live in PROOF's D-transform chart (t=4, kappa=2, C=y); (a1)-(a6) live "
      "in GGV3's gamma-reduced chart, three automorphisms downstream.  No "
      "in-repo map carries a point of one to a point of the other, so 'these "
      "four points violate (a6)' is INFERRED from '(50,75) is dead', not "
      "computed.  Building that map is the deliverable this file identifies",
      True)


# ===========================================================================
head("F.  DOES THE KILLER TRANSFER TO (75,125)?")
# ===========================================================================
check("F1  THE GAMMA LAYER TRANSFERS VERBATIM.  F_2(2,3)/75, F_3(3,2)/75 and "
      "F_2(3,5)/125 carry the SAME A_0 = (5,20) in the atlas, and the gamma "
      "layer is a function of A_0 ALONE (GGV1 conditions (5)-(9) read only "
      "(u,v)).  So all three get the same unique branch f=(4,16), the same "
      "A_0'=(1,0), the same d = gcd(f1-1,f2-1) = 3, the same bound gamma <= 4 "
      "and the same admissible set {2,3,4} -- including the same standing "
      "gamma = 4 obligation",
      tuple(r5075["A0"]) == tuple(r7550["A0"]) == tuple(r75125["A0"]) == (5, 20)
      and list(g520.values()) == [[2, 3, 4]]
      and list(g520.keys()) == [(4, 16)])

red1 = pr.case_f2(1)
NP1 = red1.reduced["standard (proportional, Prop 8.2(1))"]["P"]
NQ1 = red1.reduced["standard (proportional, Prop 8.2(1))"]["Q"]
check("F2  THE KILL LAYER DOES NOT.  polygon_reduction.case_f2(1) is labelled "
      "(75,125) and its COMPUTED reduced polygons give deg(P_1) = %d, "
      "deg(Q_1) = %d -- not (10,15).  Every one of GGV3's (a1)-(a6) is a "
      "statement about a reduced pair of degrees (10,15); none of them is a "
      "statement about (15,25)"
      % (max(p[1] for p in NP1), max(p[1] for p in NQ1)),
      red1.tag == "F2_j1_75_125" and red1.signature["degs"] == (75, 125)
      and (max(p[1] for p in NP1), max(p[1] for p in NQ1)) == (15, 25)
      and red1.mn == (3, 5))

# The E-depth law, PROVED: in (Z^m)_{-k} with Z = x^c + ... the deepest
# coefficient reached is Z_{-(k + c(m-1))}, because the only linear-in-Z_{-j}
# term is m * x^{c(m-1)} * Z_{-j}.
def deepest(E, Zs):
    return max(k for k in range(1, 40)
               if -k in Zs and any(e.has(Zs[-k]) for e in E.values()))


d3 = deepest(E3, Z3s)
d2 = deepest(E2s_, Z2s)
check("F3  THE E-DEPTH LAW, PROVED and confirmed at both published charts.  In "
      "(Z^m)_{-k} with Z = x^c + ... the only term linear in a deep "
      "coefficient is m*x^{c(m-1)}*Z_{-j}, so the deepest index reached is "
      "K + c(m-1) on the P-side and L + c(n-1) on the Q-side; equating them "
      "gives K - L = c(n-m).  gamma=3 (c=2): 5-3 = 2 = 2*(3-2), deepest %d = "
      "5 + 2.  gamma=2 (c=3): 8-5 = 3 = 3*(3-2), deepest %d = 8 + 3"
      % (d3, d2),
      d3 == 7 == 5 + 2 * (2 - 1) == 3 + 2 * (3 - 1)
      and d2 == 11 == 8 + 3 * (2 - 1) == 5 + 3 * (3 - 1)
      and 5 - 3 == 2 * (3 - 2) and 8 - 5 == 3 * (3 - 2))

check("F4  CONSEQUENCE at (m,n) = (3,5) [the (a1) shape at (3,5) is INFERRED "
      "from one anchor -- P = C^m, Q = C^n + lam C^{m-n} + F, verified only at "
      "(2,3), in both charts].  The law forces K - L = c*(5-3) = 2c instead of "
      "c, and the P-side family becomes (Z^3)_{-k}, not (Z^2)_{-k}.  So the "
      "(75,125) E-system is a DIFFERENT system, not the (50,75) one reindexed, "
      "and in particular there is no derived counterpart of C3's collapse "
      "E_6 -> 3 C_{-1}C_{-2}, which is what generates the whole kill",
      2 * (5 - 3) == 4 != 2 * (3 - 2))

MISSING = ["(a3)  the leading forcing term F_{-1} = y^7  (the UNIT that makes "
           "C_{-1}, C_{-2} monomials at all)",
           "(a5)  the cap law deg_y(C_{-k}) <= k+2  (= Step 2 of the gamma-"
           "window compiler, NOT STARTED per SESSION_HANDOFF)",
           "(a6)  the primitivity depth -10  (no in-repo derivation of the "
           "required-nonzero slot exists, at any corner)"]
check("F5  VERDICT: the killer does NOT transfer to (75,125) as it stands.  "
      "What transfers is (i) the gamma layer verbatim and (ii) the SCHEMA "
      "'required-nonzero below the forced floor'.  What is missing is exactly "
      "three derivations at (m,n) = (3,5): %s.  Two of the three are not "
      "derived at (2,3) either -- GGV3 asserts them -- so the gap is not a "
      "transcription gap, it is the compiler" % "; ".join(MISSING),
      len(MISSING) == 3)

check("F6  and the ONE quantity that already discriminates the two rows in the "
      "shipped atlas is the same integer pair: G3's gate lam >= m PASSES at "
      "F_2(2,3)/75 (lam = 2, m = 2) and FAILS at F_2(3,5)/125 (lam = 2, "
      "m = 3).  Same corner, same lam; only (m,n) differs",
      r5075["gates"]["G3"]["verdict"] == "PASS"
      and r75125["gates"]["G3"]["verdict"] == "FAIL"
      and r5075["gates"]["G3"]["sub"]["lam"]["lam"]
      == r75125["gates"]["G3"]["sub"]["lam"]["lam"] == "2")


# ===========================================================================
head("Z.  PROVENANCE -- the GGV3 transcription vs the GGV3 source")
# ===========================================================================
# Every "GGV3 says X" above answered from paper_src/upstream_quotes.json.  Where
# the copyrighted .tex is on disk we re-derive all of them plus its sha256; in a
# public clone there is nothing to re-derive and we SAY so rather than counting
# a silent pass.  A transcription is a second copy of the source, and two copies
# is where the errors live -- so they get cross-checked, not separately asserted.
_uq_res, _uq_checked = uq.verify_against_tex()
_ggv3 = [r for r in _uq_res if r[0].startswith("GGV3") or r[0].startswith("ctl.")]
if "GGV3" in _uq_checked:
    check("Z1  the GGV3 probes above re-derive from the local .tex "
          "(%d probes + sha256)" % len(_ggv3),
          all(r[1] for r in _ggv3),
          "; ".join("%s -- %s" % (r[0], r[2]) for r in _ggv3 if not r[1]))
else:
    print("[NOTE] Z1  no local GGV3 .tex -- the A/C probes were answered from "
          "the pinned transcription and NOT re-derived here (expected in a "
          "public clone; add paper_src/1406.0886_GGV3.tex to re-derive)")


# ===========================================================================
if _fail:
    print()
    print("FAILURES (%d):" % len(_fail))
    for f in _fail:
        print("   - %s" % f)
    raise SystemExit(1)
print("ALL %d CHECKS PASSED" % _ok[0])
if not QUIET:
    print("""
VERDICT
  1. THE WITNESSES SURVIVE AT (50,75), and (50,75) is where they are MOST
     constrained: the caps yplace_transfer tests them against are corner
     (5,20)'s own, the tightest of the four class corners.  Nothing in the
     transferred chain reads (deg P, deg Q), so the same four points serve
     F_2(2,3)/75 and its P<->Q swap F_3(3,2)/75 identically.        [sec.B]

  2. WHAT KILLS THEM is corner primitivity -- GGV3 (a6): one named window
     coefficient c_{0,-10} required NONZERO, at a depth the window equations
     force to vanish (forced floor -6).  Derived here from (a1)-(a6) and the
     Z-series, not replayed: E_1..E_8 and E_1..E_13 reproduce GGV3 term for
     term, and the two displayed conclusions come out of the elimination.
     The zero margin is deg_y C_{-1} + deg_y C_{-2} <= 3 + 4 = 7 = deg F_{-1}.
                                                                    [sec.C]

  3. WHY WE MISS IT.  Every condition our transfer imposes is CLOSED -- five
     equations, six degree caps, five order floors.  None can require a
     coefficient to be nonzero.  Run the ENDPOINT_CONTRACT kill predicate on
     both sides: with GGV3's contract it returns [(0,-10)]; with ours it
     returns [] even after granting us e = gamma*y^9 and B = y^21/gamma,
     because a forced MONOMIAL puts its nonzero slot exactly AT the floor and
     the predicate needs one BELOW it.
     So there is no contradiction -- the transfer is SOUND and INCOMPLETE,
     and (50,75) is the first EXTERNAL confirmation that the missing
     ingredient is the one PROOF sec.7.4(c) and YPLACE_TRANSFER sec.8 lead #2
     already name.                                                  [sec.D]

  4. IT DOES NOT TRANSFER TO (75,125).  The gamma layer does -- same A_0 =
     (5,20), same unique branch, same gamma in {2,3,4}, same gamma=4
     obligation.  The kill layer does not: the reduced pair is (15,25) not
     (10,15), and the E-depth law K - L = c(n-m) moves from c to 2c, so the
     collapse E_6 -> 3 C_{-1}C_{-2} that generates the kill has no counterpart.
     Three derivations are missing at (3,5) -- (a3), (a5), (a6) -- and two of
     them are not derived at (2,3) either.                          [sec.F]
""")
