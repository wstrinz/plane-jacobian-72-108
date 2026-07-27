#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""corner_atlas.py -- the CROSS-CORNER ATLAS of GGV5's 34 possible counterexamples.

Companion document: CORNER_ATLAS.md.  Machine-readable output: corner_atlas.json.

WHAT THIS IS.  One row per published possible counterexample with
max(deg P, deg Q) <= 150 (GGV5 = paper_src/1708.07936_GGV5.tex sec. "Possible
counterexamples ... <= 150", the 13 family rows at tex:1802-1814 and the 21
sporadic chain rows at tex:1828-1836 / 1848-1858 / 1869).  Every published datum
carries its exact tex line.  Every row is then run through the five gates that
this campaign proved, and the rows are clustered by gate signature.

    G1  CHART / DICTIONARY   polygon_reduction.chart_exponent + has_retraction
    G2  TORIC ADMISSIBILITY  (t+1) | (4t+9)  <=>  t = 4        [toric_general B4]
    G3  SLICE CASCADE        gcd(m,n)=1, lam >= m, N_Q >= D_P+D_Q  [contact_lemma]
    G4  BELYI PASSPORT       u*kappa = m_f+n_f-1 over the branch sweep
                                                             [passport_75_125]
    G5  DIVISOR SYZYGY       q_window | w(e)                  [q_window_theorem]

DISCIPLINE (this repo has been bitten by vacuous checks repeatedly):

  * TRANSCRIPTION IS CROSS-CHECKED, NOT TRUSTED.  Two independent published
    identities are re-derived from the transcribed data and must hold on all 34
    rows / all 13 family rows respectively:
      (X1)  max{deg P, deg Q} == v11(A_0) * max(m,n)      -- all 34 rows
      (X2)  (m+n)*q*k - n*(q*l - p) == k                  -- all 13 family rows
    (X1) in particular validates the family A_0 column, which is NOT printed in
    the counterexample table and had to be joined in from the chain-data table
    at tex:1678-1694 / 1709-1715.  A wrong A_0 would break (X1) immediately.

  * EVERY GATE IS MUTATION-TESTED.  A gate that returns one verdict for the
    whole population is reported as NON-DISCRIMINATING, in the report and in the
    JSON, rather than being passed off as a filter.  Two of the eight sub-gates
    ARE non-discriminating on this population and are labelled as such.

  * MISSING INPUT != FALSE.  Where a gate's input is not published and not
    derivable in-repo, the verdict is UNKNOWN and the JSON names the exact
    missing quantity.

  * THE RETRACTION GUARD IS USED, NEVER BYPASSED.  t comes from
    polygon_reduction.chart_exponent (INFERRED, not citable as published); the
    chart data (deg C, ord C) comes from polygon_reduction.corner_chart_data,
    which refuses GGV5's final-corner dictionary off the retraction shape.

Run:
    python corner_atlas.py             # full report + write corner_atlas.json
    python corner_atlas.py --quiet     # one line; exit 0 iff every check passes
    python corner_atlas.py --no-write  # do not touch corner_atlas.json
"""

from __future__ import annotations

import argparse
import json
import sys
from fractions import Fraction
from math import gcd

import polygon_reduction as pr
import passport_75_125 as pas
import q_window_theorem as qwt

TEX = "paper_src/1708.07936_GGV5.tex"


# ===========================================================================
# 1.  THE PUBLISHED DATA.  Transcribed verbatim; every field carries its line.
# ===========================================================================
# The chain-data tables (GGV5 tex:1678-1694 length 1, tex:1709-1715 length 2).
# fields: A_0, A_0', final corner (p, l_final, b_final), k, m(j), n(j), chain_len
FAMILY_CHAIN = {
    #        A_0       A_0'    final (p,l,q)  k   m=(m0,dm)  n=(n0,dn) len line
    "F_1":  ((4, 12), (1, 0), (7, 4, 3),  1, (3, 2), (4, 3),  1, 1678),
    "F_2":  ((5, 20), (1, 0), (7, 5, 2),  1, (2, 1), (3, 2),  1, 1679),
    "F_3":  ((5, 20), (1, 0), (8, 5, 3),  1, (3, 4), (2, 3),  1, 1680),
    "F_7":  ((6, 15), (1, 0), (7, 3, 4),  1, (2, 1), (7, 4),  1, 1684),
    "F_8":  ((6, 15), (1, 0), (8, 3, 5),  1, (3, 2), (7, 5),  1, 1685),
    "F_9":  ((7, 21), (1, 0), (11, 7, 2), 1, (2, 1), (3, 2),  1, 1686),
    "F_11": ((7, 21), (1, 0), (13, 7, 3), 2, (2, 1), (5, 3),  1, 1688),
    "F_17": ((9, 24), (1, 0), (11, 3, 8), 1, (2, 5), (3, 8),  1, 1694),
    # length-2 families: chain is A_0 -> A_1=(14\4,6) -> A_2 = the final corner
    "F_22": ((8, 24), (2, 0), (5, 4, 2),  1, (2, 1), (3, 2),  2, 1713),
    "F_24": ((8, 24), (2, 0), (19, 8, 3), 1, (3, 2), (4, 3),  2, 1715),
}
# the intermediate corner of the two length-2 families (same for both)
FAMILY_A1 = {"F_22": "(14\\4,6)", "F_24": "(14\\4,6)"}

# The 13 family rows of the counterexample table, GGV5 tex:1802-1814.
FAMILY_ROWS = [
    ("F_1",  (3, 4),  64, 1802, True),
    ("F_1",  (5, 7), 112, 1803, False),
    ("F_2",  (2, 3),  75, 1804, True),
    ("F_2",  (3, 5), 125, 1805, False),
    ("F_3",  (3, 2),  75, 1806, True),
    ("F_7",  (2, 7), 147, 1807, False),
    ("F_8",  (3, 7), 147, 1808, False),
    ("F_9",  (2, 3),  84, 1809, True),
    ("F_9",  (3, 5), 140, 1810, False),
    ("F_11", (2, 5), 140, 1811, False),
    ("F_17", (2, 3),  99, 1812, True),
    ("F_22", (2, 3),  96, 1813, True),    # the starred row, see GGV5 tex:1818
    ("F_24", (3, 4), 128, 1814, False),
]

# The 21 sporadic chain rows.  fields:
#   A_0, [intermediate corners...], final (p,l,q), (m,n), maxdeg, chain_len, line
SPORADIC_ROWS = [
    ((7, 35),  [], (19, 7, 5),  (2, 3), 126, 1, 1828),
    ((7, 42),  [], (13, 7, 6),  (3, 2), 147, 1, 1829),
    ((7, 42),  [], (13, 7, 6),  (2, 3), 147, 1, 1830),
    ((8, 28),  [], (7, 4, 3),   (3, 4), 144, 1, 1831),
    ((8, 28),  [], (11, 4, 7),  (3, 2), 108, 1, 1832),   # *** (72,108), CLOSED
    ((9, 36),  [], (17, 9, 4),  (3, 2), 135, 1, 1833),
    ((9, 36),  [], (17, 9, 4),  (2, 3), 135, 1, 1834),
    ((11, 33), [], (19, 4, 8),  (2, 3), 132, 1, 1835),
    ((12, 33), [], (11, 3, 8),  (2, 3), 135, 1, 1836),
    ((8, 32),  ["(8,28)"],   (11, 4, 7), (3, 2), 120, 2, 1848),
    ((8, 40),  ["(8,28)"],   (11, 4, 7), (3, 2), 144, 2, 1849),
    ((9, 27),  ["(9,24)"],   (11, 3, 8), (2, 3), 108, 2, 1850),
    ((9, 36),  ["(9,24)"],   (11, 3, 8), (2, 3), 135, 2, 1851),
    ((10, 40), ["(16\\5,6)"],  (23, 10, 3), (3, 2), 150, 2, 1852),
    ((10, 40), ["(18\\5,8)"],  (8, 5, 3),   (3, 2), 150, 2, 1853),
    ((12, 30), ["(16\\3,10)"], (11, 6, 3),  (3, 2), 126, 2, 1854),
    ((12, 36), ["(12,33)"],  (11, 3, 8),  (2, 3), 144, 2, 1855),
    ((12, 36), ["(9,24)"],   (11, 3, 8),  (2, 3), 144, 2, 1856),
    ((12, 36), ["(21\\4,9)"],  (19, 4, 8),  (2, 3), 144, 2, 1857),
    ((12, 36), ["(21\\4,9)"],  (12, 4, 5),  (2, 3), 144, 2, 1858),
    ((12, 36), ["(12,30)", "(16\\3,10)"], (11, 6, 3), (3, 2), 144, 3, 1869),
]

# The one row GGV5 itself discards, tex:1818 + Proposition "caso antisimetrico".
GGV5_SELF_DISCARDED = {("F_22", (2, 3))}

# The one row this campaign has CLOSED.
CLOSED_ROWS = {((8, 28), (3, 2), 108)}

# The single case whose Phi signature (ord_y Phi, M, deg_y Phi) is published /
# derived in this repo.  Source: window_functions_75_125 / toric_general E1.
PHI_KNOWN = {((8, 28), (3, 2)): dict(ordPhi=204, M=17, degPhi=238, lam=2)}


# ===========================================================================
# 2.  HARNESS
# ===========================================================================
class Ledger:
    def __init__(self, quiet):
        self.quiet, self.rows = quiet, []

    def head(self, s):
        if not self.quiet:
            print("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)

    def ck(self, name, ok, detail=""):
        self.rows.append((name, bool(ok), detail))
        if not self.quiet:
            print("  [%s] %-52s %s" % ("PASS" if ok else "FAIL", name, detail))
        return bool(ok)

    def note(self, s):
        if not self.quiet:
            print("       . " + s)

    def report(self):
        bad = [r for r in self.rows if not r[1]]
        print("\n%s  corner_atlas: %d/%d checks pass"
              % ("FAIL" if bad else "OK  ", len(self.rows) - len(bad), len(self.rows)))
        for n, _, d in bad:
            print("   FAILED: %s   %s" % (n, d))
        return 1 if bad else 0


# ===========================================================================
# 3.  ROW CONSTRUCTION
# ===========================================================================
def solve_j(m0dm, n0dn, mn):
    """The unique j >= 0 with (m(j), n(j)) == mn, or with the pair swapped.

    GGV5 tex:1794 states explicitly that only the cases satisfying the first
    Diophantine equality are listed and 'the other cases can be obtained by
    swapping m with n', so a swapped match is a legitimate identification.
    """
    (m0, dm), (n0, dn) = m0dm, n0dn
    for j in range(0, 64):
        p = (m0 + dm * j, n0 + dn * j)
        if p == tuple(mn):
            return j, False
        if p == (mn[1], mn[0]):
            return j, True
    return None, None


def build_rows():
    rows = []
    for fam, mn, deg, line, red in FAMILY_ROWS:
        A0, A0p, fin, k, mj, nj, clen, cline = FAMILY_CHAIN[fam]
        j, swapped = solve_j(mj, nj, mn)
        rows.append(dict(
            id="%s(%d,%d)/%d" % (fam, mn[0], mn[1], deg),
            kind="FAMILY", family=fam, j=j, mn_swapped_vs_family_law=swapped,
            A0=list(A0), A0prime=list(A0p),
            intermediate=[FAMILY_A1[fam]] if fam in FAMILY_A1 else [],
            final_corner=dict(p=fin[0], l_final=fin[1], b_final=fin[2]),
            k=k, m=mn[0], n=mn[1], max_deg=deg, chain_length=clen,
            provenance=dict(
                counterexample_row="%s:%d" % (TEX, line),
                chain_row="%s:%d" % (TEX, cline),
                red_in_paper=red),
        ))
    for A0, mid, fin, mn, deg, clen, line in SPORADIC_ROWS:
        rows.append(dict(
            id="(%d,%d)/(%d,%d)/%d" % (A0[0], A0[1], mn[0], mn[1], deg),
            kind="SPORADIC", family=None, j=None,
            mn_swapped_vs_family_law=None,
            A0=list(A0), A0prime=None, intermediate=list(mid),
            final_corner=dict(p=fin[0], l_final=fin[1], b_final=fin[2]),
            k=None, m=mn[0], n=mn[1], max_deg=deg, chain_length=clen,
            provenance=dict(counterexample_row="%s:%d" % (TEX, line),
                            chain_row=None, red_in_paper=False),
        ))
    return rows


# ===========================================================================
# 4.  THE GATES
# ===========================================================================
def gate1_chart(r):
    """G1  chart exponent + dictionary validity.  polygon_reduction sec.0b."""
    a0, b0 = r["A0"]
    t = pr.chart_exponent(a0, b0)
    retr = pr.has_retraction(a0, b0, t)
    lf = r["final_corner"]["l_final"]
    cd = pr.corner_chart_data(a0, b0, l_final=lf,
                              b_final=r["final_corner"]["b_final"],
                              who=r["id"])
    verdict = "PASS" if (retr and lf == t) else "FAIL"
    return dict(
        gate="G1_chart_dictionary", verdict=verdict,
        t=t, kappa=t - 2, retraction=retr, l_final=lf, l_final_eq_t=(lf == t),
        deg_C=cd["deg_C"], ord_C=cd["ord_C"], C_is_monomial=cd["monomial"],
        dictionary_trust=("TRUSTED" if retr else "SUSPECT"),
        detail=("retraction shape b0 == t*(a0-1) holds (%d == %d*%d) and GGV5's "
                "l_final = %d agrees with t = ceil(b0/a0) = %d"
                % (b0, t, a0 - 1, lf, t)) if retr else
               ("retraction shape FAILS: %d != %d*%d = %d, so GGV5's final-corner "
                "dictionary is REFUSED here; C is a monomial (deg C = ord C = 1) "
                "and l_final = %d is %sequal to t = %d"
                % (b0, t, a0 - 1, t * (a0 - 1), lf,
                   "" if lf == t else "NOT ", t)),
    )


def gate2_toric(r, g1):
    """G2  toric admissibility.  toric_general.py B4: (t+1) | (4t+9) <=> t = 4."""
    t = g1["t"]
    ok = (4 * t + 9) % (t + 1) == 0
    return dict(gate="G2_toric_admissible", verdict="PASS" if ok else "FAIL",
                t=t, exponent=((4 * t + 9) // (t + 1)) if ok else None,
                detail="4t+9 = 4(t+1)+5 so (t+1)|5; t = %d %s the unique solution "
                       "t = 4" % (t, "IS" if ok else "is NOT"))


def gate3_cascade(r, g1):
    """G3  slice-cascade gates.  CONTACT_LEMMA.md / contact_lemma.py.

    gcd(m,n) = 1                                             [F7,  PROVED]
    lam >= a,  lam = (deg_y Phi - ord_y Phi)/M               [D5,  PROVED]
    N_Q >= D_P + D_Q,  N_Q = (b+1)*ell - 1, D_P = a, D_Q = b [F3,  PROVED]

    ell is the chart parameter and equals t (MINIMAL_CORE.md sec.2.0 line 150,
    "`ell` the chart parameter (`t = l`)").  NOTE contact_lemma.py:1115 still
    hardcodes the PRE-REPAIR row (m,n,ell,lam) = (3,5,5,3) for (75,125); with
    the repaired t = 4 that row should be (3,5,4,0).  REPORTED, not fixed --
    it does not block this atlas, which derives ell from the guarded t.
    """
    a, b = sorted((r["m"], r["n"]))
    sub = {}
    sub["gcd"] = dict(verdict="PASS" if gcd(a, b) == 1 else "FAIL",
                      value=gcd(a, b))
    # lam
    key = (tuple(r["A0"]), (r["m"], r["n"]))
    if key in PHI_KNOWN:
        lam = Fraction(PHI_KNOWN[key]["degPhi"] - PHI_KNOWN[key]["ordPhi"],
                       PHI_KNOWN[key]["M"])
        sub["lam"] = dict(verdict="PASS" if lam >= a else "FAIL",
                          lam=str(lam), source="published/derived Phi signature "
                          "(ord=%d, M=%d, deg=%d)" % (PHI_KNOWN[key]["ordPhi"],
                                                      PHI_KNOWN[key]["M"],
                                                      PHI_KNOWN[key]["degPhi"]))
    elif g1["C_is_monomial"]:
        # deg_y(Phi) - ord_y(Phi) = N*(deg C - ord C) = 0 identically when C is a
        # monomial -- the (75,125) argument (toric_general E3), which depends on
        # nothing but deg C = ord C.
        sub["lam"] = dict(verdict="FAIL", lam="0",
                          source="C is a MONOMIAL at this corner (retraction "
                          "fails), so deg_y(Phi) - ord_y(Phi) = N*(deg C - ord C) "
                          "= 0 identically and lam = 0 < %d = min(m,n).  Same "
                          "argument as toric_general E3 at (75,125)." % a)
    else:
        sub["lam"] = dict(verdict="UNKNOWN", lam=None,
                          missing="ord_y(Phi), deg_y(Phi), M at this corner -- "
                          "the corner retracts so C is NOT a monomial and lam is "
                          "not forced to 0, but Phi has never been derived here")
    ell = g1["t"]
    NQ = (b + 1) * ell - 1
    sub["N_Q"] = dict(verdict="PASS" if NQ >= a + b else "FAIL",
                      N_Q=NQ, D_P=a, D_Q=b, ell=ell)
    vs = [v["verdict"] for v in sub.values()]
    verdict = "FAIL" if "FAIL" in vs else ("UNKNOWN" if "UNKNOWN" in vs else "PASS")
    return dict(gate="G3_slice_cascade", verdict=verdict, sub=sub)


def _branches(a0, b0):
    """Every reduction branch legal at this corner (passport_75_125 rules r1-r8)."""
    mu = (b0 - 1) // a0
    c = b0 - mu * a0
    q = gcd(a0, b0)
    zdeg = a0 // q
    ss = [mu] + ([mu - 1] if mu - 1 >= 1 else [])
    splits = [None] + [j * q for j in range(1, zdeg) if j * q >= c]
    ens = [dict()] + [dict(en_k=k, en_swap=sw)
                      for k in range(0, b0 // a0 + 2) if (k + 1) * a0 < b0
                      for sw in (False, True)]
    for s in ss:
        for sp in splits:
            for en in ens:
                yield s, sp, en


def gate4_belyi(r, g1):
    """G4  Belyi passport.  passport_75_125.py: build N(P), N(Q) from (a0,b0,m,n)
    through the published reduction rules and sweep every branch for a face
    passing the gate u*kappa = m_f + n_f - 1 with a valid ramification datum.

    kappa = t - 2 comes from the GUARD, not from l_final.
    """
    a0, b0 = r["A0"]
    kappa = g1["kappa"]
    tried = legal = 0
    faces = mech = 0
    hits = []
    for mm, nn in sorted({(r["m"], r["n"]), (r["n"], r["m"])}):
        for s, sp, en in _branches(a0, b0):
            try:
                red = pas.Reduction("atlas", a0, b0, (mm, nn), s=s,
                                    split_e=sp, **en)
            except Exception:
                continue
            tried += 1
            if not red.legal:
                continue
            legal += 1
            rowsf, out = pas.passport_from_polygon(red.NP, red.NQ, kappa)
            if out:
                faces += len(out)
                nv = sum(1 for d in rowsf if d.get("mechanism"))
                mech += nv
                hits.append(dict(mn=[mm, nn], s=s, split_e=sp, en=en,
                                 n_faces=len(out), n_valid=nv))
    # the proportional class is provably empty unless kappa >= 2(m+n)-1
    need = 2 * (r["m"] + r["n"]) - 1
    return dict(gate="G4_belyi_passport",
                verdict="PASS" if faces else "FAIL",
                kappa=kappa, branches_tried=tried, branches_legal=legal,
                n_admissible_faces=faces, n_valid_mechanisms=mech,
                hits=hits,
                proportional_class=dict(
                    verdict="FAIL", kappa=kappa, required_kappa_at_least=need,
                    detail="the proportional branch needs kappa >= 2(m+n)-1 = %d "
                           "(passport_75_125 U2); kappa = %d" % (need, kappa)))


def gate5_syzygy(r, g1):
    """G5  divisor-syzygy criterion q_window | w(e).  MINIMAL_CORE.md sec.4.

    q_window = M/gcd(M,H), M = t(a+b)-(kappa+1), H = q(a+b)-1  [q_window_theorem,
    proved symbolically].  We feed it the GUARDED (t, kappa, q = ord C).

    q_window = 1 divides EVERY w(e), so the criterion holds unconditionally there
    -- that is the only verdict obtainable without the corner's split
    enumeration.  Otherwise UNKNOWN, missing w(e).  (MINIMAL_CORE.md sec."Negative
    I could not overcome" already records that the split data exists only at the
    (8,28) corner.)
    """
    a, b = sorted((r["m"], r["n"]))
    t, kap, q = g1["t"], g1["kappa"], g1["ord_C"]
    M, H, g, qw = qwt.q_window(t, kap, q, a, b)
    # the as-published (unguarded) reading, for comparison only
    lf = g1["l_final"]
    Mu, Hu, gu, qwu = qwt.q_window(lf, lf - 2, r["final_corner"]["b_final"], a, b)
    return dict(gate="G5_divisor_syzygy",
                verdict="PASS" if qw == 1 else "UNKNOWN",
                q_window=qw, M=M, H=H, gcd=g,
                q_window_unguarded=qwu, M_unguarded=Mu,
                missing=None if qw == 1 else
                        "w(e) -- the split-weight enumeration at this corner; "
                        "q_window = %s != 1 so the criterion is w(e)-dependent "
                        "and only decidable per split" % qw,
                detail="q_window = 1 divides every w(e), so the syzygy's carry "
                       "obstruction vanishes on EVERY split at this corner"
                       if qw == 1 else
                       "q_window = %s; the criterion holds exactly on the splits "
                       "with q_window | w(e)" % qw)


# ===========================================================================
# 5.  MAIN
# ===========================================================================
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--no-write", action="store_true")
    args = ap.parse_args()
    L = Ledger(args.quiet)

    rows = build_rows()

    # ---------------------------------------------------------------- census
    L.head("0.  CENSUS -- the published population")
    L.ck("A0 the atlas carries exactly GGV5's 34 rows", len(rows) == 34,
         "%d FAMILY + %d SPORADIC"
         % (sum(1 for r in rows if r["kind"] == "FAMILY"),
            sum(1 for r in rows if r["kind"] == "SPORADIC")))
    L.ck("A1 sporadic chain-length split is 9 / 11 / 1 as published",
         [sum(1 for r in rows if r["kind"] == "SPORADIC" and r["chain_length"] == c)
          for c in (1, 2, 3)] == [9, 11, 1],
         "GGV5 tex:1825 ('9 other possible pairs'), tex:1839 ('11 other'), "
         "tex:1862 ('another possible pair with a complete chain of length 3')")
    L.ck("A2 every family row identified a unique j >= 0 in its family law",
         all(r["j"] is not None for r in rows if r["kind"] == "FAMILY"),
         "j = %s" % [r["j"] for r in rows if r["kind"] == "FAMILY"])

    # ------------------------------------------- transcription cross-checks
    L.head("X.  TRANSCRIPTION CROSS-CHECKS -- published identities, re-derived")
    badX1 = [r["id"] for r in rows
             if r["max_deg"] != sum(r["A0"]) * max(r["m"], r["n"])]
    L.ck("X1 max{deg P, deg Q} == v11(A_0) * max(m,n) on all 34 rows",
         not badX1,
         "this is the check that validates the FAMILY A_0 column, which the "
         "counterexample table does NOT print -- it had to be joined in from "
         "the chain-data table.  A mis-joined A_0 breaks it at once.")
    if badX1:
        L.note("offenders: %s" % badX1)
    badX2 = []
    for r in rows:
        if r["kind"] != "FAMILY":
            continue
        f, m, n, k = r["final_corner"], r["m"], r["n"], r["k"]
        if (m + n) * f["b_final"] * k - n * (f["b_final"] * f["l_final"]
                                             - f["p"]) != k:
            badX2.append(r["id"])
    L.ck("X2 Diophantine identity (m+n)qk - n(ql - p) = k on all 13 family rows",
         not badX2, "GGV5's own admissibility identity; it also confirms that "
         "the FINAL corner of the two length-2 families is A_2, not A_1")
    if badX2:
        L.note("offenders: %s" % badX2)
    # X1 must be a real test, not an identity: perturb one A_0 and watch it break
    mut = [r for r in rows
           if r["max_deg"] == (sum(r["A0"]) + 1) * max(r["m"], r["n"])]
    L.ck("X1-MUT the degree identity is not vacuous", not mut,
         "shifting v11(A_0) by 1 breaks it on every one of the 34 rows")

    # --------------------------------------------------------------- gates
    L.head("1.  THE FIVE GATES")
    for r in rows:
        g1 = gate1_chart(r)
        g2 = gate2_toric(r, g1)
        g3 = gate3_cascade(r, g1)
        g4 = gate4_belyi(r, g1)
        g5 = gate5_syzygy(r, g1)
        r["gates"] = dict(G1=g1, G2=g2, G3=g3, G4=g4, G5=g5)
        r["signature"] = "%s|%s|%s|%s|%s" % (g1["verdict"], g2["verdict"],
                                             g3["verdict"], g4["verdict"],
                                             g5["verdict"])

    def vs(key):
        return [r["gates"][key]["verdict"] for r in rows]

    for key, name in (("G1", "chart/dictionary"), ("G2", "toric admissibility"),
                      ("G3", "slice cascade"), ("G4", "Belyi passport"),
                      ("G5", "divisor syzygy")):
        v = vs(key)
        dist = {x: v.count(x) for x in sorted(set(v))}
        L.ck("%s %s DISCRIMINATES across the 34" % (key, name),
             len(set(v)) >= 2, "verdict distribution %s" % dist)

    # sub-gate honesty: which sub-gates are non-discriminating?
    L.head("1b.  SUB-GATE MUTATION AUDIT -- which sub-tests actually separate")
    sub_gcd = {r["gates"]["G3"]["sub"]["gcd"]["verdict"] for r in rows}
    L.ck("G3a gcd(m,n)=1 is NON-DISCRIMINATING here, and is reported as such",
         sub_gcd == {"PASS"},
         "all 34 published pairs are coprime by construction (GGV5 Definition "
         "'mn families' requires it).  The gate is PROVED and load-bearing "
         "in general (contact_lemma F7 exhibits (2,4),(3,6),(4,6),(2,6) "
         "counterexamples) but it filters NOTHING on this population.")
    sub_nq = {r["gates"]["G3"]["sub"]["N_Q"]["verdict"] for r in rows}
    L.ck("G3c N_Q >= D_P+D_Q is NON-DISCRIMINATING here, and is reported as such",
         sub_nq == {"PASS"},
         "min slack = %d over the 34 rows; contact_lemma F3 notes it 'only "
         "binds at small ell'"
         % min(r["gates"]["G3"]["sub"]["N_Q"]["N_Q"]
               - r["gates"]["G3"]["sub"]["N_Q"]["D_P"]
               - r["gates"]["G3"]["sub"]["N_Q"]["D_Q"] for r in rows))
    sub_lam = [r["gates"]["G3"]["sub"]["lam"]["verdict"] for r in rows]
    L.ck("G3b lam >= m DOES discriminate", len(set(sub_lam)) >= 2,
         "%s" % {x: sub_lam.count(x) for x in sorted(set(sub_lam))})
    prop = {r["gates"]["G4"]["proportional_class"]["verdict"] for r in rows}
    L.ck("G4-prop the proportional class is empty on ALL 34 (non-discriminating)",
         prop == {"FAIL"},
         "min margin 2(m+n)-1 - kappa = %d, i.e. never close"
         % min(r["gates"]["G4"]["proportional_class"]["required_kappa_at_least"]
               - r["gates"]["G4"]["kappa"] for r in rows))

    # G1 mutation: the historic (5,20) bug -- "some l satisfies b0 = l(a0-1)"
    exists_l = [r["id"] for r in rows
                if r["A0"][0] > 1 and r["A0"][1] % (r["A0"][0] - 1) == 0]
    guarded = [r["id"] for r in rows if r["gates"]["G1"]["retraction"]]
    L.ck("G1-MUT the quantifier matters: 'SOME l' != 'the l that is used'",
         set(exists_l) != set(guarded) and len(exists_l) > len(guarded),
         "'exists l with b0 = l(a0-1)' admits %d rows, the guarded test admits "
         "%d.  The %d extra rows are exactly the historic l=5 trap: %s"
         % (len(exists_l), len(guarded), len(exists_l) - len(guarded),
            sorted(set(exists_l) - set(guarded))))
    L.ck("G2-MUT the divisibility test is not vacuous",
         len([t for t in range(2, 200) if (4 * t + 4) % (t + 1) == 0]) == 198
         and len([t for t in range(2, 200) if (4 * t + 9) % (t + 1) == 0]) == 1,
         "numerator 4t+4 passes for EVERY t; 4t+9 passes for exactly one")
    # cross-implementation control on t: polygon_reduction.chart_exponent
    # (ceil(b0/a0)) vs passport_75_125's independent rule (r1) l = mu+1,
    # mu = floor((b0-1)/a0).  Two files, two derivations, same number.
    bad_t = [r["id"] for r in rows
             if pas.Reduction("ctrl", r["A0"][0], r["A0"][1],
                              (r["m"], r["n"])).l != r["gates"]["G1"]["t"]]
    L.ck("G1-CTRL two independent implementations of t agree on all 34",
         not bad_t,
         "polygon_reduction.chart_exponent = ceil(b0/a0) vs passport_75_125 "
         "rule (r1) l = floor((b0-1)/a0)+1")

    # G5 control: the guarded q_window must reproduce the three values this
    # repo already knows, from an entirely different derivation path.
    known_qw = {"F_2(2,3)/75": 17,        # (50,75), q_window_theorem KNOWN_CASES
                "F_2(3,5)/125": 29,       # (75,125), window_functions q_window(3)
                "(8,28)/(3,2)/108": 1}    # (72,108), the home case
    got_qw = {r["id"]: r["gates"]["G5"]["q_window"] for r in rows
              if r["id"] in known_qw}
    L.ck("G5-CTRL guarded q_window reproduces the three repo-known values",
         got_qw == known_qw,
         "expected %s, got %s -- the closed form M/gcd(M,H) fed with the "
         "GUARDED (t,kappa,ord C) lands on q_window_theorem's KNOWN_CASES and "
         "on window_functions.q_window(3) = 12*3-7 = 29" % (known_qw, got_qw))

    # G4 mutation: the sweep is not empty-by-construction
    L.ck("G4-MUT the passport sweep is not empty by construction",
         all(r["gates"]["G4"]["branches_legal"] > 0 for r in rows)
         and any(r["gates"]["G4"]["n_admissible_faces"] > 0 for r in rows),
         "every row had >= %d legal reduction branches to test, and at least "
         "one row DOES produce admissible faces"
         % min(r["gates"]["G4"]["branches_legal"] for r in rows))

    # ------------------------------------------------- Deliverable 3 claims
    L.head("2.  THE (8,28) UNIQUENESS CLAIM")
    both = sorted(r["id"] for r in rows
                  if r["gates"]["G1"]["retraction"] and r["gates"]["G2"]["verdict"] == "PASS")
    spor_both = sorted(r["id"] for r in rows if r["kind"] == "SPORADIC"
                       and r["gates"]["G1"]["retraction"]
                       and r["gates"]["G2"]["verdict"] == "PASS")
    retr_corners = sorted({tuple(r["A0"]) for r in rows
                           if r["gates"]["G1"]["retraction"]})
    L.ck("D3a among the SPORADIC rows, (8,28) is the ONLY corner passing both "
         "retraction and t = 4",
         {tuple(r["A0"]) for r in rows if r["kind"] == "SPORADIC"
          and r["gates"]["G1"]["retraction"]
          and r["gates"]["G2"]["verdict"] == "PASS"} == {(8, 28)},
         "sporadic rows passing both: %s" % spor_both)
    L.ck("D3b extended to the FAMILY rows: still only (8,28)",
         {tuple(r["A0"]) for r in rows
          if r["gates"]["G1"]["retraction"]
          and r["gates"]["G2"]["verdict"] == "PASS"} == {(8, 28)},
         "all 34 rows passing both: %s" % both)
    L.ck("D3c the conjunction is doing work -- retraction alone admits 4 corners",
         len(retr_corners) == 4 and retr_corners == [(6, 15), (8, 28), (9, 24),
                                                     (12, 33)],
         "corners with the retraction shape: %s; of these only (8,28) has "
         "t = 4 (the other three all have t = 3), so neither condition alone "
         "isolates it" % (retr_corners,))
    L.ck("D3d at every retracting corner GGV5's l_final agrees with "
         "ceil(b0/a0)", all(r["gates"]["G1"]["l_final_eq_t"] for r in rows
                            if r["gates"]["G1"]["retraction"]),
         "6 rows / 4 distinct corners -- INDEPENDENT corroboration of the "
         "INFERRED chart_exponent rule beyond the 5 published reductions it "
         "was validated on.  The converse is false: %d non-retracting rows "
         "also happen to have l_final = t."
         % sum(1 for r in rows if not r["gates"]["G1"]["retraction"]
               and r["gates"]["G1"]["l_final_eq_t"]))

    # ------------------------------------------------------ the clustering
    L.head("3.  CLUSTERING BY GATE SIGNATURE  (G1|G2|G3|G4|G5)")
    clusters = {}
    for r in rows:
        clusters.setdefault(r["signature"], []).append(r["id"])
    for sig in sorted(clusters, key=lambda s: (-len(clusters[s]), s)):
        L.note("%-34s  n=%2d   %s" % (sig, len(clusters[sig]),
                                      ", ".join(sorted(clusters[sig]))))
    L.ck("C1 the 34 rows collapse to a SMALL number of gate signatures",
         len(clusters) <= 8,
         "%d distinct signatures over 34 rows -- the frontier is mechanism-"
         "structured, not 34 singletons" % len(clusters))
    L.ck("C2 the clustering is not trivial (more than one class)",
         len(clusters) >= 2, "%d classes" % len(clusters))
    sizes = sorted((len(v) for v in clusters.values()), reverse=True)
    L.ck("C3 the population is concentrated: the two largest signatures cover "
         ">= 80% of it", sum(sizes[:2]) >= 0.8 * len(rows),
         "class sizes %s; top two = %d/%d = %.0f%%"
         % (sizes, sum(sizes[:2]), len(rows),
            100.0 * sum(sizes[:2]) / len(rows)))
    allpass = sorted(r["id"] for r in rows
                     if set(r["signature"].split("|")) == {"PASS"})
    L.ck("C4 exactly ONE row passes all five gates, and it is the case this "
         "campaign closed", allpass == ["(8,28)/(3,2)/108"],
         "all-PASS rows: %s.  No other row in the published population has the "
         "full (72,108) mechanism stack." % allpass)
    nofail = sorted(r["id"] for r in rows if "FAIL" not in r["signature"])
    L.ck("C5 no row is 'all-UNKNOWN' -- every row is decided by at least one "
         "gate", all("PASS" in r["signature"] or "FAIL" in r["signature"]
                     for r in rows),
         "rows with no FAIL at all: %s" % nofail)

    # ------------------------------------------------------ evaluability
    L.head("4.  EVALUABILITY -- how much of the atlas is actually decided")
    n_unknown_any = sum(1 for r in rows
                        if "UNKNOWN" in r["signature"])
    per_gate_unknown = {k: sum(1 for r in rows
                               if r["gates"][k]["verdict"] == "UNKNOWN")
                        for k in ("G1", "G2", "G3", "G4", "G5")}
    L.ck("E1 G1, G2, G4 are fully evaluated on all 34 rows",
         all(per_gate_unknown[k] == 0 for k in ("G1", "G2", "G4")),
         "UNKNOWN counts per gate: %s" % per_gate_unknown)
    L.ck("E2 the residual UNKNOWNs are named, not silent",
         all(r["gates"]["G3"]["sub"]["lam"].get("missing")
             for r in rows if r["gates"]["G3"]["sub"]["lam"]["verdict"] == "UNKNOWN")
         and all(r["gates"]["G5"]["missing"]
                 for r in rows if r["gates"]["G5"]["verdict"] == "UNKNOWN"),
         "%d rows carry at least one UNKNOWN; every one names its missing input"
         % n_unknown_any)

    # ------------------------------------------------------------- output
    out = dict(
        schema="corner_atlas/1",
        source=dict(tex=TEX,
                    section="Possible counterexamples with max(deg P,deg Q) <= 150",
                    section_line=1792,
                    family_table="%s:1802-1814" % TEX,
                    chain_table_len1="%s:1678-1694" % TEX,
                    chain_table_len2="%s:1709-1715" % TEX,
                    sporadic_len1="%s:1828-1836" % TEX,
                    sporadic_len2="%s:1848-1858" % TEX,
                    sporadic_len3="%s:1869" % TEX),
        gate_status=dict(
            G1="INFERRED (chart_exponent rule is not published; the retraction "
               "shape and its two published counterexamples are CHECKED)",
            G2="PROVED (toric_general.py B4, 38/38)",
            G3="PROVED (contact_lemma.py, 64/64) -- but lam requires the "
               "corner's Phi, which is published nowhere",
            G4="PROVED gate + CHECKED reduction engine (passport_75_125.py "
               "81/81, engine reproduces 5 published GGHV22 reductions)",
            G5="PROVED criterion (MINIMAL_CORE.md, exhaustive over 1.5M "
               "splits) -- but w(e) requires the corner's split enumeration"),
        n_rows=len(rows),
        clusters={s: sorted(v) for s, v in clusters.items()},
        unknown_counts=per_gate_unknown,
        ggv5_self_discarded=["%s(%d,%d)" % (f, m, n)
                             for (f, (m, n)) in sorted(GGV5_SELF_DISCARDED)],
        closed_by_this_campaign=["(8,28)/(3,2)/108"],
        rows=rows,
    )
    if not args.no_write:
        with open("corner_atlas.json", "w") as fh:
            json.dump(out, fh, indent=1, sort_keys=False)
        if not args.quiet:
            print("\n  wrote corner_atlas.json (%d rows)" % len(rows))

    return L.report()


if __name__ == "__main__":
    sys.exit(main())
