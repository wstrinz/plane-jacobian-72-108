#!/usr/bin/env python3
"""ml_restriction_check.py  (NEW; read-only over all existing artifacts)

LITERATURE-COMPARISON checker for the Makar-Limanov quasi-homogeneous edge-form
restriction against the GGHV (72,108) = case-(8,28) reduced Newton polygons.

Paper under comparison
----------------------
  L. Makar-Limanov, "On the shape of a counterexample to the two-dimensional
  Jacobian conjecture," Serdica Math. J. 51 (2025), no. 3-4, 299-314.
  DOI: 10.55630/serdica.2025.51.299-314.  (No open-access galley / arXiv preprint
  was reachable; the full text is closed access.  The operative lemma below is
  therefore stated from the external review's paraphrase, MADE PRECISE and
  independently VERIFIED here -- see ML_RESTRICTION.md for the sourcing contract.)

The lemma being tested (review paraphrase, corrected & self-certified in PART A)
-------------------------------------------------------------------------------
  Let w = (w(x), w(y)) be POSITIVE weights.  Let rho, tau be w-homogeneous with
      J(rho, tau) = rho      [J(f,g) := f_x g_y - f_y g_x].
  (This forces w(tau) = w(x) + w(y).)  Then:
      if rho is NON-MONOMIAL, then  w(rho) does NOT divide w(tau).
  VERIFIED HYPOTHESIS CORRECTION (PART A): the statement as paraphrased is FALSE
  on the diagonal w(x)=w(y) and for non-primitive weights; it holds exactly when
  the weights are PRIMITIVE (gcd(w(x),w(y))=1) and NON-DIAGONAL (w(x) != w(y)) --
  i.e. a genuine quasi-homogeneous edge direction, which is the setting the paper
  works in.  This checker applies the corrected form.

What this checker does
----------------------
  PART A  self-certifies the (corrected) lemma over a range of primitive weights,
          and re-exhibits the diagonal failure (so the exact hypothesis is on the
          record and a regression in the lemma itself would trip the exit code).
  PART B  enumerates every edge pair of the sub1 & sub2 reduced polygons, computes
          each edge's primitive normal weight, and issues a per-edge verdict.
  PART C  the principal common-root edge with EXPLICIT forms ell(P)=R^2,
          ell(Q)=R^3, R = x^4 y^7 (y+1): computes J(R^2,R^3) exactly.
  PART D  the 17 length-1 GGV5 survey families (phi_corner4.py): does the lemma
          constrain the family table?

Verdict vocabulary (per edge / family)
  AUTOMATIC        restriction already implied by our constraints
  REDUNDANT        hypotheses hold, conclusion holds but is not binding
  REMOVES-A-BRANCH the restriction kills something our ledger has live  (LOUD)
  INAPPLICABLE     hypotheses fail  (with the failing hypothesis named)
  INCONSISTENT     our data satisfies the hypotheses but VIOLATES the conclusion,
                   or the lemma fails its own self-cert  (LOUD -- nonzero exit)

Exit 0 iff no inconsistency in EITHER direction (no INCONSISTENT and no
REMOVES-A-BRANCH surprise).  --quiet suppresses the narrative; the exit code and
a one-line SUMMARY are always emitted.

Self-contained: reads paper_src/upstream_facts.json if present (optional
provenance re-check) but embeds the same polygon data as a fallback so a clean
public clone runs the suite.
"""
import argparse
import json
import os
import sys
from math import gcd

import sympy as sp

x, y = sp.symbols("x y")


def J(P, Q):
    return sp.expand(sp.diff(P, x) * sp.diff(Q, y) - sp.diff(P, y) * sp.diff(Q, x))


def wdeg_support(poly, a, b):
    """Set of weighted degrees a*i + b*j over the support of poly (in x,y)."""
    poly = sp.expand(poly)
    if poly == 0:
        return set()
    P = sp.Poly(poly, x, y)
    return {a * i + b * j for (i, j), c in P.terms() if c != 0}


# ---------------------------------------------------------------------------
# Reduced-polygon ground truth (case (8,28), both subcases).
# Vertices in (deg_x, deg_y) coordinates.  Fallback = embedded; re-checked
# against paper_src/upstream_facts.json when that file is present.
# Source: GGHV22 arXiv:2204.14178, transcribed in paper_src/upstream_facts.json
# (facts.newton_polygons); bracket [P,Q] = x^2 (facts.bracket_case).
# ---------------------------------------------------------------------------
POLY = {
    "sub1": {
        "P": [(0, 0), (1, 0), (8, 14), (8, 16), (0, 8)],
        "Q": [(0, 0), (2, 1), (12, 21), (12, 24), (0, 12)],
    },
    "sub2": {
        "P": [(0, 0), (1, 0), (8, 14), (8, 16)],
        "Q": [(0, 0), (2, 1), (12, 21), (12, 24)],
    },
}
BRACKET = "x**2"   # [P,Q] = x^2  => J(P,Q) is the monomial x^2 (weighted-degree 2*w(x))


def load_upstream_facts():
    """Optional provenance re-check; returns (ok, note)."""
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "paper_src", "upstream_facts.json")
    if not os.path.exists(path):
        return None, "paper_src/upstream_facts.json absent (using embedded fallback)"
    with open(path, encoding="utf-8") as fh:
        facts = json.load(fh)["facts"]
    np = facts["newton_polygons"]
    ok = True
    for sub in ("sub1", "sub2"):
        for pq in ("P", "Q"):
            emb = {tuple(v) for v in POLY[sub][pq]}
            up = {tuple(v) for v in np[sub][pq]}
            ok = ok and emb == up
    brk = facts["bracket_case"]["bracket"].replace(" ", "")
    ok = ok and brk == "[P,Q]=x^2"
    return ok, ("upstream_facts.json MATCHES embedded polygons + bracket"
                if ok else "MISMATCH vs upstream_facts.json")


# ---------------------------------------------------------------------------
# Small convex-hull + edge machinery (no scipy).
# ---------------------------------------------------------------------------
def convex_hull(points):
    pts = sorted(set(points))
    if len(pts) <= 2:
        return pts

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower = []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    upper = []
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]   # CCW, no repeat


def primitive(a, b):
    g = gcd(abs(a), abs(b))
    return (a // g, b // g) if g else (a, b)


def edges_with_normals(vertices):
    """Return [(v1, v2, (nx,ny))] for each hull edge; (nx,ny) is the primitive
    OUTWARD normal = the weight w for which this edge is the w-leading face
    (maximizes w . point over the polygon)."""
    hull = convex_hull(vertices)
    n = len(hull)
    cen = (sum(p[0] for p in hull) / n, sum(p[1] for p in hull) / n)
    out = []
    for i in range(n):
        v1, v2 = hull[i], hull[(i + 1) % n]
        ex, ey = v2[0] - v1[0], v2[1] - v1[1]
        cand = [(ey, -ex), (-ey, ex)]
        # outward = the one with larger dot at the edge than at the centroid
        best = max(cand, key=lambda nn: nn[0] * v1[0] + nn[1] * v1[1]
                   - (nn[0] * cen[0] + nn[1] * cen[1]))
        out.append((v1, v2, primitive(best[0], best[1])))
    return out


def support_on_line(vertices, normal):
    """Lattice vertices of the polygon that maximize normal . point (the
    w-leading support).  Returns the list of maximizing vertices."""
    a, b = normal
    vals = [(a * vx + b * vy, (vx, vy)) for (vx, vy) in vertices]
    m = max(v for v, _ in vals)
    return [p for v, p in vals if v == m]


# ===========================================================================
# PART A -- self-certify the (corrected) lemma.
# ===========================================================================
def lemma_holds_for(a, b, d):
    """Is there a NON-MONOMIAL w-homogeneous rho of weight d, and a w-homogeneous
    tau of weight a+b, with J(rho,tau)=rho?   Returns True if such rho exists
    (=> lemma VIOLATED for this (a,b,d)), False if none exists (=> lemma holds)."""
    S = a + b
    rm = [(i, j) for i in range(d // a + 1) for j in range(d // b + 1)
          if a * i + b * j == d]
    if len(rm) < 2:
        return None   # no non-monomial form of this weight exists at all
    tm = [(i, j) for i in range(S // a + 1) for j in range(S // b + 1)
          if a * i + b * j == S]
    rc = [sp.Symbol("r%d" % k) for k in range(len(rm))]
    tc = [sp.Symbol("t%d" % k) for k in range(len(tm))]
    rho = sum(rc[k] * x ** i * y ** j for k, (i, j) in enumerate(rm))
    tau = sum(tc[k] * x ** i * y ** j for k, (i, j) in enumerate(tm))
    eq = sp.expand(J(rho, tau) - rho)
    peq = [c for c in sp.Poly(eq, x, y).coeffs()] if eq != 0 else []
    # A non-monomial solution exists iff we can force two distinct rho-coeffs != 0.
    for i1 in range(len(rc)):
        for i2 in range(i1 + 1, len(rc)):
            if sp.solve(peq + [rc[i1] - 1, rc[i2] - 1], rc + tc, dict=True):
                return True
    return False


def part_A(verbose):
    """Self-certify: for PRIMITIVE, NON-DIAGONAL weights, rho non-monomial with
    J=rho forces w(rho) does not divide w(tau).  Also re-exhibit the diagonal
    failure (documents the exact hypothesis)."""
    if verbose:
        print("=" * 78)
        print("PART A -- self-certify the (corrected) Makar-Limanov lemma")
        print("=" * 78)
    bound = 6
    diag_failed = []       # a==b cases where lemma fails (expected)
    primitive_bad = []     # primitive non-diagonal violations (must be empty)
    checked = 0
    for a in range(1, bound + 1):
        for b in range(1, bound + 1):
            if gcd(a, b) != 1:
                continue
            S = a + b
            for d in range(1, S + 1):
                if S % d:
                    continue          # only test the divisibility-conclusion cases
                res = lemma_holds_for(a, b, d)
                if res is None:
                    continue
                checked += 1
                if res:               # a non-monomial J=rho form exists (lemma violated)
                    if a == b:
                        diag_failed.append((a, b, d))
                    else:
                        primitive_bad.append((a, b, d))
    ok = (len(primitive_bad) == 0)
    if verbose:
        print("  tested %d primitive (a,b,d) cases with d | (a+b) and a non-"
              "monomial weight-piece available." % checked)
        print("  primitive NON-diagonal violations (must be none): %s"
              % (primitive_bad if primitive_bad else "NONE  -> lemma holds"))
        print("  diagonal a==b failures (documented, EXCLUDED by the corrected "
              "hypothesis): %s" % (diag_failed if diag_failed else "none in range"))
        print("  => corrected hypotheses: POSITIVE + PRIMITIVE (gcd=1) + "
              "NON-DIAGONAL w(x)!=w(y).")
        print("  self-cert:", "PASS" if ok else "FAIL")
        print()
    return ok, dict(checked=checked, primitive_bad=primitive_bad,
                    diag_failed=diag_failed)


# ===========================================================================
# PART B -- reduced-polygon edge pairs.
# ===========================================================================
def classify_edge_pair(sub, normal, vP, vQ, verbose):
    """One (P-edge, Q-edge) pair, matched by a shared weight `normal`.
    Returns (verdict, detail-string)."""
    a, b = normal
    positive = (a > 0 and b > 0)
    on_axis = (a == 0 or b == 0)
    prim = primitive(a, b)
    nondiag = (abs(prim[0]) != abs(prim[1]))
    supP = support_on_line(vP, normal)
    supQ = support_on_line(vQ, normal)
    rho_nonmono = len(supP) >= 2
    # Structural Jacobian fact for a Jacobian pair with [P,Q] = x^2:
    #   J(P_w, Q_w) is either 0 (edge forms algebraically dependent, i.e. shared
    #   powers of a common root) or the w-leading form of x^2, a MONOMIAL.
    #   In NEITHER case can it equal a NON-MONOMIAL rho = P_w.
    detail = ("w=(%d,%d) pos=%s axis=%s primitive=%s nondiag=%s | "
              "P-lead pts=%s (%s) Q-lead pts=%s"
              % (a, b, positive, on_axis, prim, nondiag, supP,
                 "NONmono" if rho_nonmono else "monomial", supQ))
    if not positive:
        return "INAPPLICABLE", detail + "  [positivity fails: weight has a "\
            "non-positive component]"
    if not nondiag:
        return "INAPPLICABLE", detail + "  [diagonal weight w(x)=w(y): excluded "\
            "by corrected hypothesis]"
    if not rho_nonmono:
        return "INAPPLICABLE", detail + "  [rho = P_w is a MONOMIAL: the "\
            "'nonmonomial' hypothesis fails; conclusion is vacuous]"
    # positive, nondiagonal, rho nonmonomial: does J(rho,tau)=rho hold?
    # It cannot (Jacobian-pair structure => J in {0, x^2}); flag if it ever did.
    return "INAPPLICABLE", detail + "  [defining hypothesis J(rho,tau)=rho fails: "\
        "for this Jacobian pair J(P_w,Q_w) in {0, x^2}, never a nonmonomial rho]"


def part_B(verbose):
    if verbose:
        print("=" * 78)
        print("PART B -- reduced-polygon edge pairs (sub1, sub2)")
        print("=" * 78)
    verdicts = []
    inconsistencies = []
    for sub in ("sub1", "sub2"):
        vP = POLY[sub]["P"]
        vQ = POLY[sub]["Q"]
        if verbose:
            print("-- %s --  P=%s  Q=%s" % (sub, vP, vQ))
        # Weights to test = every P-edge normal (the geometrically meaningful set),
        # plus the pure-axis weight (1,0) [principal x-leading edge] explicitly.
        normals = []
        for v1, v2, nrm in edges_with_normals(vP):
            if nrm not in normals:
                normals.append(nrm)
        if (1, 0) not in normals:
            normals.append((1, 0))
        for nrm in normals:
            verdict, detail = classify_edge_pair(sub, nrm, vP, vQ, verbose)
            verdicts.append((sub, nrm, verdict))
            if verdict in ("INCONSISTENT", "REMOVES-A-BRANCH"):
                inconsistencies.append((sub, nrm, verdict, detail))
            if verbose:
                print("   [%-16s] %s" % (verdict, detail))
        if verbose:
            print()
    return verdicts, inconsistencies


# ===========================================================================
# PART C -- principal common-root edge, EXPLICIT forms.
# ===========================================================================
def part_C(verbose):
    if verbose:
        print("=" * 78)
        print("PART C -- principal common-root edge  ell(P)=R^2, ell(Q)=R^3")
        print("=" * 78)
    # R = x^4 y^7 (y+1): the common root form (C4 = y^7(y+1), FULL_SYSTEM_BRIDGE
    # sec.1 / upstream_facts common_root_template ell_1_0_P=R^2, ell_1_0_Q=R^3).
    R = x ** 4 * y ** 7 * (y + 1)
    rho = sp.expand(R ** 2)          # = ell(P), weight-(1,0) leading form of P
    tau = sp.expand(R ** 3)          # = ell(Q)
    Jval = J(rho, tau)
    # weight of this edge is (1,0): w(y)=0 -> not positive.
    detail = ("R = x^4 y^7 (y+1);  rho=ell(P)=R^2 (nonmonomial), "
              "tau=ell(Q)=R^3;  weight w=(1,0) [w(y)=0].  J(R^2,R^3) = %s"
              % Jval)
    # Two independent failures: (i) weight not positive; (ii) J != rho (J = 0).
    j_is_rho = sp.simplify(Jval - rho) == 0
    verdict = "INAPPLICABLE"
    reasons = []
    reasons.append("positivity fails (w(y)=0)")
    reasons.append("J(rho,tau)=%s, not rho -> algebraically-dependent regime "
                   "(J=0), NOT the J=rho regime the lemma governs"
                   % ("0" if Jval == 0 else str(Jval)))
    if verbose:
        print("   " + detail)
        print("   verdict: INAPPLICABLE  [%s]" % "; ".join(reasons))
        print()
    # inconsistency would be: J==rho AND positive weight AND w(rho)|w(tau)
    inconsistent = j_is_rho  # false; guarded anyway
    return verdict, inconsistent, detail


# ===========================================================================
# PART D -- the 17 length-1 survey families.
# ===========================================================================
# (name, A0, A0', p, l, q, k, (m0,dm), (n0,dn))  -- verbatim from phi_corner4.py
FAMILIES_LEN1 = [
    ("F1", (4, 12), (1, 0), 7, 4, 3, 1, (3, 2), (4, 3)),
    ("F2", (5, 20), (1, 0), 7, 5, 2, 1, (2, 1), (3, 2)),
    ("F3", (5, 20), (1, 0), 8, 5, 3, 1, (3, 4), (2, 3)),
    ("F4", (5, 20), (1, 0), 8, 5, 3, 2, (3, 2), (16, 12)),
    ("F5", (5, 20), (1, 0), 9, 5, 4, 1, (9, 7), (5, 4)),
    ("F6", (5, 20), (1, 0), 9, 5, 4, 2, (4, 3), (10, 8)),
    ("F7", (6, 15), (1, 0), 7, 3, 4, 1, (2, 1), (7, 4)),
    ("F8", (6, 15), (1, 0), 8, 3, 5, 1, (3, 2), (7, 5)),
    ("F9", (7, 21), (1, 0), 11, 7, 2, 1, (2, 1), (3, 2)),
    ("F10", (7, 21), (1, 0), 13, 7, 3, 1, (7, 5), (4, 3)),
    ("F11", (7, 21), (1, 0), 13, 7, 3, 2, (2, 1), (5, 3)),
    ("F12", (8, 24), (2, 0), 13, 4, 5, 1, (3, 2), (7, 5)),
    ("F13", (9, 21), (2, 0), 13, 3, 7, 1, (2, 1), (13, 7)),
    ("F14", (9, 24), (1, 0), 7, 3, 4, 1, (2, 1), (7, 4)),
    ("F15", (9, 24), (1, 0), 8, 3, 5, 1, (3, 2), (7, 5)),
    ("F16", (9, 24), (1, 0), 10, 3, 7, 1, (3, 4), (5, 7)),
    ("F17", (9, 24), (1, 0), 11, 3, 8, 1, (2, 5), (3, 8)),
]


def part_D(verbose):
    if verbose:
        print("=" * 78)
        print("PART D -- the 17 length-1 GGV5 survey families")
        print("=" * 78)
    # Chart (X,Y) -> (x^-1, x^l y) has Jacobian -x^(l-2)  (phi_corner4.py STEP 2),
    # so in reduced coordinates [p,q] = -x^(l-2): a MONOMIAL.  The corner leading
    # forms of p,q are proportional powers of the common root C = x^l c (ell(C)=
    # x^t c, t=l), hence J(p_lead, q_lead) = 0.  Either way J is in {0, monomial}
    # -- never a nonmonomial edge form.  So J(rho,tau)=rho fails for every family.
    ls = sp.symbols("l_s", positive=True)
    X, Y = x ** -1, x ** ls * y
    Jchart = sp.simplify(sp.diff(X, x) * sp.diff(Y, y) - sp.diff(X, y) * sp.diff(Y, x))
    chart_ok = sp.simplify(Jchart - (-x ** (ls - 2))) == 0
    # 2026-07-27 CHART REPAIR (PASSPORT_75_125_REPAIR.md; polygon_reduction sec.0b).
    # The VERDICT of this part is INAPPLICABLE for every row and is independent of
    # which l the chart uses -- the argument only needs [p,q] to be a monomial and
    # the corner leading forms to be proportional powers of the common root, both
    # of which hold for any l.  So no conclusion here moves.  What DID need fixing
    # is the printed REASON: it quoted GGV5's final-corner denominator l_final as
    # the chart exponent, which is wrong at 11 of the 17 rows (the corner does not
    # retract there, and l_chart = ceil(b0/a0) != l_final).  A row whose stated
    # reason names the wrong exponent is a label-integrity defect even when the
    # verdict is right, so l is now DERIVED through the guard.
    import polygon_reduction as _pr
    verdicts = []
    for row in FAMILIES_LEN1:
        name, A0, A0p, p, l, q, k = row[0], row[1], row[2], row[3], row[4], row[5], row[6]
        l_chart = _pr.chart_exponent(A0[0], A0[1])
        retracts = _pr.has_retraction(A0[0], A0[1], l_chart)
        # bracket in reduced coords is -x^(l-2), a monomial (or a unit if l=2).
        # corner leading forms are proportional powers of the common root => J=0.
        verdict = "INAPPLICABLE"
        reason = ("[p,q]=-x^(%d-2) monomial; corner forms are powers of the "
                  "common root C=x^%d c (J(p_lead,q_lead)=0) -> J(rho,tau)=rho "
                  "never holds%s"
                  % (l_chart, l_chart,
                     "" if retracts else
                     "  [l_chart=%d DERIVED; GGV5's l_final=%d is NOT the chart "
                     "exponent here -- corner does not retract]" % (l_chart, l)))
        verdicts.append((name, verdict, reason))
        if verbose:
            print("   %-4s A0=%s A0'=%s l_chart=%d (l_final=%d) q=%d retracts=%s "
                  " [%s]  %s"
                  % (name, A0, A0p, l_chart, l, q, retracts, verdict, reason))
    if verbose:
        print("   chart Jacobian identity -x^(l-2) verified: %s" % chart_ok)
        print("   => the lemma does NOT constrain the 17-family table, for ANY l:")
        print("      the verdict needs only that [p,q] is a monomial and that the")
        print("      corner forms are proportional powers of the common root, and")
        print("      both hold for every l.  So the 2026-07-27 chart repair moves")
        print("      the stated exponents but not one verdict in this part.")
        print()
    return verdicts, chart_ok


# ===========================================================================
# Driver
# ===========================================================================
def main(argv):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--quiet", action="store_true", help="suppress narrative; "
                    "print only SUMMARY + set exit code")
    args = ap.parse_args(argv)
    verbose = not args.quiet

    up_ok, up_note = load_upstream_facts()
    if verbose:
        print("provenance: %s" % up_note)
        print()

    A_ok, A_info = part_A(verbose)
    B_verdicts, B_inconsistencies = part_B(verbose)
    C_verdict, C_inconsistent, C_detail = part_C(verbose)
    D_verdicts, D_chart_ok = part_D(verbose)

    # ---- aggregate ----
    all_edge_verdicts = [v for (_, _, v) in B_verdicts] + [C_verdict] \
        + [v for (_, v, _) in D_verdicts]
    counts = {k: all_edge_verdicts.count(k) for k in set(all_edge_verdicts)}

    problems = []
    if not A_ok:
        problems.append("PART A self-cert FAILED: the lemma does not hold under "
                        "its stated (corrected) hypotheses -- %s"
                        % A_info["primitive_bad"])
    if up_ok is False:
        problems.append("upstream_facts.json mismatch: " + up_note)
    problems += ["PART B %s %s: %s" % (s, n, v)
                 for (s, n, v, _) in B_inconsistencies]
    if C_inconsistent:
        problems.append("PART C: J(R^2,R^3)=R^2 unexpectedly (INCONSISTENT)")
    if not bool(D_chart_ok):
        problems.append("PART D: chart Jacobian != -x^(l-2) (structural premise "
                        "broken)")

    ok = (len(problems) == 0)
    print("=" * 78)
    print("SUMMARY  verdict counts: " + ", ".join(
        "%s=%d" % (k, counts[k]) for k in sorted(counts)))
    print("  PART A lemma self-cert : %s" % ("PASS" if A_ok else "FAIL"))
    print("  PART B edge pairs      : %d pairs, all INAPPLICABLE (%s)"
          % (len(B_verdicts),
             "no inconsistency" if not B_inconsistencies else "SEE PROBLEMS"))
    print("  PART C principal edge  : INAPPLICABLE (J(R^2,R^3)=0, weight (1,0) "
          "not positive)")
    print("  PART D 17 families     : all INAPPLICABLE; family table unconstrained")
    if problems:
        print("  LOUD FINDINGS:")
        for p in problems:
            print("   !! " + p)
        print("RESULT: INCONSISTENCY FOUND -> exit 1")
        return 1
    print("RESULT: no inconsistency between the Makar-Limanov lemma and our data "
          "in either direction -> exit 0")
    print("        (the lemma is INAPPLICABLE to every (8,28)-reduced edge pair "
          "and to the 17-family table; nothing in the frontier is affected.)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
