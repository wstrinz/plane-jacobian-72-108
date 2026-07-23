#!/usr/bin/env python3
"""validate_candidate.py -- exact checker for a claimed *plane* Jacobian
counterexample (JC(2)), staged in order of cheapness.

Given two polynomials P, Q in x, y this reports, in order:

  1. BASIC   -- is [P,Q] := P_x Q_y - P_y Q_x a nonzero *constant* (exact sympy)?
               Are P, Q non-invertible-looking (total degree > 1)?  Reports degrees.
  2. BOUND   -- is max(deg P, deg Q) < 125?  If so, GGHV22 (arXiv:2204.14178,
               Prop 4.3) forces the pair to be (72,108)-shaped up to symmetry.
               Any other shape CONTRADICTS the published theorem.
  3. SHAPE   -- (relevant for the (72,108) target) compute the Newton-polygon
               corners of P, Q and compare to the Prop 4.3 subcase polygons as
               transcribed in T3_WINDOW_AUDIT.md section 1 (read at runtime -- the
               numbers are sourced from that doc, not hand-copied here).
  4. NEXT    -- if 1-3 leave the claim standing, name (do NOT run) the deeper
               framework checks: the f31 necessary condition and the frontier
               degree-state lists.

SCOPE / HONESTY.  This tool checks *necessary* conditions only.  It can REFUTE a
claim (non-constant Jacobian, wrong degree shape, wrong Newton polygon) but it
CANNOT certify that a surviving pair is a genuine counterexample -- that requires
the full framework (f31 window infeasibility, frontier states, the T6 reduction
debt), none of which is implemented here.

Usage:
    python validate_candidate.py CANDIDATE.py      # a .py module defining P, Q
    python validate_candidate.py CANDIDATE.json    # {"P": "...", "Q": "..."}
    python validate_candidate.py --self-test       # run the built-in self-tests

See VALIDATE_CANDIDATE.md for the input formats and a verdict glossary.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys

import sympy
from sympy import symbols, Poly, expand, diff, sympify

X, Y = symbols("x y")

# Repo docs used as *sources* (read at runtime; nothing hand-copied from them).
_HERE = os.path.dirname(os.path.abspath(__file__))
T3_DOC = os.path.join(_HERE, "T3_WINDOW_AUDIT.md")
STATE_DOC = os.path.join(_HERE, "STATE.md")
PHASE_D_SUB2 = os.path.join(_HERE, "phase_d_states_sub2.json")
PHASE_D_SUB1 = os.path.join(_HERE, "phase_d_states_sub1.json")

BOUND = 125          # GGHV22 max-degree threshold below which (72,108) is forced
TARGET = (72, 108)   # the (deg P, deg Q) shape forced below the bound


# --------------------------------------------------------------------------- #
# Input loading
# --------------------------------------------------------------------------- #
def load_candidate(path):
    """Load (P, Q) as sympy expressions in x, y from a .py or .json file.

    .py  : an ordinary Python module that assigns `P` and `Q` (sympy exprs).
           `x` and `y` are pre-injected so the module can use them directly.
    .json: an object with string fields "P" and "Q" that sympy.sympify parses
           over the symbols x, y (optional "vars" defaults to ["x","y"]).
    """
    ext = os.path.splitext(path)[1].lower()
    if ext == ".py":
        ns = {"x": X, "y": Y, "sympy": sympy, "symbols": symbols}
        with open(path, "r", encoding="utf-8") as fh:
            code = fh.read()
        exec(compile(code, path, "exec"), ns)          # noqa: S102 (local tool)
        if "P" not in ns or "Q" not in ns:
            raise ValueError(f"{path}: module must define both `P` and `Q`.")
        P, Q = ns["P"], ns["Q"]
    elif ext == ".json":
        with open(path, "r", encoding="utf-8") as fh:
            obj = json.load(fh)
        for key in ("P", "Q"):
            if key not in obj:
                raise ValueError(f"{path}: JSON must contain a '{key}' field.")
        local = {"x": X, "y": Y}
        P = sympify(str(obj["P"]), locals=local)
        Q = sympify(str(obj["Q"]), locals=local)
    else:
        raise ValueError(f"Unsupported input extension '{ext}' (use .py or .json).")

    P, Q = sympify(P), sympify(Q)
    stray = (P.free_symbols | Q.free_symbols) - {X, Y}
    if stray:
        raise ValueError(
            f"P, Q must be polynomials in x, y only; found extra symbols {stray}."
        )
    return P, Q


# --------------------------------------------------------------------------- #
# Geometry / algebra helpers
# --------------------------------------------------------------------------- #
def total_degree(expr):
    expr = expand(expr)
    if expr == 0:
        return -1
    return int(Poly(expr, X, Y).total_degree())


def bracket(P, Q):
    """The Jacobian bracket [P,Q] = P_x Q_y - P_y Q_x, fully expanded."""
    return expand(diff(P, X) * diff(Q, Y) - diff(P, Y) * diff(Q, X))


def exponent_points(expr):
    """Exponent tuples (i, j) of the monomials actually present in expr."""
    expr = expand(expr)
    if expr == 0:
        return []
    return [(int(m[0]), int(m[1])) for m in Poly(expr, X, Y).monoms()]


def convex_hull(points):
    """Vertices of the convex hull of integer points (Andrew's monotone chain).

    Returns the hull vertices as a set of (i, j) tuples; collinear points that
    are not extreme are dropped, so the result is exactly the polygon corners.
    """
    pts = sorted(set(points))
    if len(pts) <= 1:
        return set(pts)

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
    return set(lower[:-1] + upper[:-1]) or set(pts)


def newton_corners(expr):
    """Newton-polygon corner set (hull vertices) of a polynomial in x, y."""
    return frozenset(convex_hull(exponent_points(expr)))


# --------------------------------------------------------------------------- #
# Prop 4.3 reference polygons -- parsed from T3_WINDOW_AUDIT.md section 1
# --------------------------------------------------------------------------- #
def parse_prop43_reference(doc=T3_DOC):
    """Read the Prop 4.3 subcase corner sets from T3_WINDOW_AUDIT.md section 1.

    Numbers are sourced from the doc (which sources them verbatim from
    arXiv:2204.14178 lines 1000-1007), not hand-copied into this file.  Returns
    {subcase_index: {"P": frozenset(pts), "Q": frozenset(pts)}} and raises if the
    doc is missing or its structure is not what we expect.
    """
    if not os.path.exists(doc):
        raise FileNotFoundError(
            f"Prop 4.3 source doc not found: {doc}\n"
            "SHAPE CHECK needs T3_WINDOW_AUDIT.md; cannot proceed without it."
        )
    with open(doc, "r", encoding="utf-8") as fh:
        text = fh.read()

    # Isolate section 1 so we don't accidentally pick up later restatements.
    m = re.search(r"##\s*1\.\s*Prop 4\.3 statement.*?(?=\n##\s)", text, re.S)
    section = m.group(0) if m else text

    pt_re = re.compile(r"\((\d+),\s*(\d+)\)")

    def pts_after(label, chunk):
        m2 = re.search(re.escape(label) + r"\s*=\s*\{([^}]*)\}", chunk)
        if not m2:
            return None
        return frozenset(
            (int(a), int(b)) for a, b in pt_re.findall(m2.group(1))
        )

    # Each subcase is introduced by "(1)" / "(2)" and carries N(P) then N(Q).
    ref = {}
    for idx_match in re.finditer(r"\((\d)\)\s*N\(P\)", section):
        idx = int(idx_match.group(1))
        chunk = section[idx_match.start():idx_match.start() + 400]
        setP = pts_after("N(P)", chunk)
        setQ = pts_after("N(Q)", chunk)
        if setP and setQ:
            ref[idx] = {"P": setP, "Q": setQ}

    # Structural sanity (validates the parse, does NOT re-supply coordinates):
    if set(ref) != {1, 2}:
        raise ValueError(
            f"Expected exactly subcases 1 and 2 in {doc}; parsed {sorted(ref)}."
        )
    expect_sizes = {1: (5, 5), 2: (4, 4)}
    for idx, (nP, nQ) in expect_sizes.items():
        got = (len(ref[idx]["P"]), len(ref[idx]["Q"]))
        if got != (nP, nQ):
            raise ValueError(
                f"Subcase {idx}: expected {(nP, nQ)} corners, parsed {got} "
                f"from {doc} -- parser or doc changed; refusing to guess."
            )
    return ref


def ref_shape_degree(ref):
    """(deg P, deg Q) of the reduced-frame reference polygons (max i+j)."""
    dP = max(i + j for i, j in ref[2]["P"])
    dQ = max(i + j for i, j in ref[2]["Q"])
    return dP, dQ


# --------------------------------------------------------------------------- #
# Reporting helpers
# --------------------------------------------------------------------------- #
class Report:
    def __init__(self):
        self.refuted = False        # a published fact is contradicted / hard fail
        self.out_of_scope = False   # valid map but not a counterexample claim
        self.lines = []

    def say(self, s=""):
        self.lines.append(s)
        print(s)

    def header(self, s):
        self.say("")
        self.say(s)
        self.say("-" * len(s))


def fmt_pts(pts):
    return "{" + ", ".join(f"({i},{j})" for i, j in sorted(pts)) + "}"


# --------------------------------------------------------------------------- #
# The staged validation
# --------------------------------------------------------------------------- #
def validate(P, Q, label="candidate"):
    r = Report()
    r.say(f"=== validate_candidate: {label} ===")
    r.say(f"P = {P}")
    r.say(f"Q = {Q}")

    # ---- Step 1: BASIC ---------------------------------------------------- #
    r.header("1. BASIC  (Jacobian bracket + non-invertibility)")
    degP, degQ = total_degree(P), total_degree(Q)
    J = bracket(P, Q)
    is_const = (J.free_symbols == set())
    is_zero = (J == 0)
    jac_ok = is_const and not is_zero

    r.say(f"deg P = {degP},  deg Q = {degQ}")
    r.say(f"[P,Q] = P_x Q_y - P_y Q_x = {J}")
    if is_zero:
        r.say("  -> VERDICT: FAIL -- [P,Q] is identically 0 (P, Q algebraically "
              "dependent). Not a Jacobian pair.")
        r.refuted = True
    elif not is_const:
        r.say("  -> VERDICT: FAIL -- [P,Q] is NON-CONSTANT. The Jacobian "
              "condition is violated; this is not a counterexample.")
        r.refuted = True
    else:
        r.say(f"  -> Jacobian condition PASSES: [P,Q] = {J} is a nonzero constant.")

    noninvertible = (degP > 1 and degQ > 1)
    if not noninvertible:
        r.say(f"  -> Degrees (deg P={degP}, deg Q={degQ}): at least one is <= 1, "
              "so the map looks INVERTIBLE (e.g. affine/triangular).")
        if jac_ok:
            r.say("     A constant Jacobian with such degrees is an automorphism, "
                  "not a counterexample -- OUT OF SCOPE.")
            r.out_of_scope = True
    else:
        r.say(f"  -> Non-invertibility OK: both degrees exceed 1.")

    if r.refuted:
        _finish(r)
        return r

    # ---- Step 2: BOUND CHECK --------------------------------------------- #
    r.header("2. BOUND CHECK  (GGHV22 degree constraint)")
    md = max(degP, degQ)
    r.say(f"max(deg P, deg Q) = {md}   (threshold {BOUND})")
    shape = (degP, degQ)
    shape_ok = shape in (TARGET, (TARGET[1], TARGET[0]))
    if not noninvertible:
        r.say("  -> Degrees put this outside the counterexample regime; the "
              "GGHV22 (72,108) constraint does not apply.")
    elif md >= BOUND:
        r.say(f"  -> max degree >= {BOUND}: OUTSIDE the range GGHV22 pins down. "
              "No published-theorem contradiction can be asserted; the (72,108) "
              "shape check below is informational only.")
    elif shape_ok:
        r.say(f"  -> max degree < {BOUND} and (deg P, deg Q) = {shape} matches "
              f"the forced shape {TARGET} (up to swap). CONSISTENT with GGHV22.")
    else:
        r.say(f"  -> *** FLAG *** max degree < {BOUND} but (deg P, deg Q) = "
              f"{shape} is NOT {TARGET} (or its swap).")
        r.say("      GGHV22 (arXiv:2204.14178, Prop 4.3) proves every plane JC "
              "counterexample with max degree < 125 has this shape up to symmetry.")
        r.say("      CAVEAT: this tool tests non-invertibility only by the "
              "degree>1 heuristic; it does NOT decide whether the map is a "
              "genuine automorphism.")
        r.say("      * If the pair is genuinely non-invertible, this is a hard "
              "CONTRADICTION: either the claimed counterexample is wrong or the "
              "published theorem is. Investigate BOTH the pair and the citation.")
        r.say("      * If the pair is actually an automorphism (invertible), it "
              "is simply NOT a counterexample -- fully consistent with GGHV22.")
        r.say("      Either way the claim 'this is a plane JC counterexample' "
              "does not stand as given.")
        r.refuted = True

    # ---- Step 3: SHAPE CHECK  (Newton polygons vs Prop 4.3) -------------- #
    r.header("3. SHAPE CHECK  (Newton polygon vs Prop 4.3 subcases)")
    try:
        ref = parse_prop43_reference()
    except (FileNotFoundError, ValueError) as e:
        r.say(f"  -> SHAPE CHECK could not run: {e}")
        _finish(r)
        return r

    cP, cQ = newton_corners(P), newton_corners(Q)
    refdeg = ref_shape_degree(ref)
    r.say(f"Newton corners of P: {fmt_pts(cP)}")
    r.say(f"Newton corners of Q: {fmt_pts(cQ)}")
    r.say("Reference (Prop 4.3, sourced from T3_WINDOW_AUDIT.md section 1,")
    r.say(" itself verbatim from arXiv:2204.14178 lines 1000-1007):")
    for idx in (1, 2):
        r.say(f"  subcase ({idx}): N(P)={fmt_pts(ref[idx]['P'])}  "
              f"N(Q)={fmt_pts(ref[idx]['Q'])}")
    r.say(f"Reference frame total degree = {refdeg}  (the reduced [P,Q]=x^2 "
          "representative).")

    # Try to match, in either P<->Q assignment, against either subcase.
    matched = None
    for idx in (1, 2):
        rP, rQ = ref[idx]["P"], ref[idx]["Q"]
        if (cP == rP and cQ == rQ) or (cP == rQ and cQ == rP):
            matched = idx
            break

    if matched:
        r.say(f"  -> Newton polygons EXACTLY match Prop 4.3 subcase ({matched}). "
              "The candidate is supplied in (or already reduced to) the "
              "normalized frame and passes the shape check.")
    else:
        # Is the candidate at full (72,108) scale rather than the reduced frame?
        scale_num = None
        if md and refdeg[1]:
            if degP * refdeg[1] == degQ * refdeg[0] and degP % refdeg[0] == 0:
                scale_num = degP // refdeg[0]
        if shape_ok or scale_num:
            r.say("  -> No exact match. NOTE: Prop 4.3's corners describe the "
                  f"REDUCED [P,Q]=x^2 representative of total degree {refdeg}. A "
                  f"raw (72,108) candidate is "
                  + (f"{scale_num}x that reduced frame" if scale_num
                     else "a scaled/un-normalized representative")
                  + ".")
            r.say("     Matching requires first running the paper's reduction "
                  "(inversion morphism phi(x)=x^-1, phi(y)=x^4 y, plus the "
                  "normalization chain). That reduction is NOT implemented here "
                  "(it is the open T6 debt -- see T3_WINDOW_AUDIT.md section 5).")
            r.say("     Therefore a non-match here is NOT a refutation; it only "
                  "means the candidate is not in reduced normalized form.")
        else:
            r.say("  -> No match, and the candidate is neither in the reduced "
                  f"frame {refdeg} nor at the (72,108) scale. Newton polygon is "
                  "inconsistent with the Prop 4.3 subcases as given.")

    _finish(r)

    # ---- Step 4: NEXT (pointers only) ------------------------------------ #
    if not r.refuted and not r.out_of_scope:
        r.header("4. NEXT  (deeper framework checks -- pointers only, not run)")
        r.say("This candidate survives the cheap necessary conditions. The "
              "repo framework would next apply (none implemented here):")
        r.say("  * f31 necessary condition -- the master elimination identity")
        r.say("      f31 * f37 * d_{-1}^21 == 0 in K[y]; the generic branch is")
        r.say("      f31 == 0, a weighted-homogeneous 102-term factor. See "
              f"STATE.md items 5-6 ({_relpath(STATE_DOC)}).")
        r.say("  * Frontier degree-states -- the surviving flag cases and their")
        r.say("      residual degree windows, in "
              f"{_relpath(PHASE_D_SUB2)} (sub2) and "
              f"{_relpath(PHASE_D_SUB1)} (sub1); overview in FRONTIER.md.")
        r.say("  * The T6 reduction debt (equation selection / normalization "
              "chain) -- T3_WINDOW_AUDIT.md section 5.")

    return r


def _relpath(p):
    try:
        return os.path.relpath(p, _HERE)
    except ValueError:
        return p


def _finish(r):
    r.header("OVERALL VERDICT")
    if r.refuted:
        r.say("REFUTED -- a necessary condition (or a published theorem) is "
              "violated. This is NOT a valid plane JC counterexample as stated.")
    elif r.out_of_scope:
        r.say("OUT OF SCOPE -- the Jacobian condition may hold, but the degrees "
              "describe an invertible map, not a counterexample claim.")
    else:
        r.say("SURVIVES the implemented necessary conditions. This does NOT "
              "certify a counterexample: the tool checks necessary conditions "
              "only and cannot prove sufficiency. See step 4 for what remains.")
    r.say("")


# --------------------------------------------------------------------------- #
# Self-test
# --------------------------------------------------------------------------- #
def self_test():
    print("##################################################################")
    print("# SELF-TEST (a): P = x, Q = y")
    print("#   Jacobian [P,Q] = 1 (constant) but the map is the identity --")
    print("#   invertible, degrees <= 1, so OUT OF SCOPE (not a counterexample).")
    print("##################################################################")
    ra = validate(X, Y, label="self-test (a): identity map P=x, Q=y")

    print("\n\n##################################################################")
    print("# SELF-TEST (b): P = x**3 + y, Q = x + y**2")
    print("#   [P,Q] = 6*x**2*y - 1 is NON-CONSTANT -> rejected at step 1.")
    print("##################################################################")
    rb = validate(X**3 + Y, X + Y**2,
                  label="self-test (b): non-constant Jacobian")

    print("\n\n##################################################################")
    print("# SHAPE-MATCHER UNIT DEMO (not one of the mandated cases):")
    print("#   synthetic polynomials whose Newton polygons EQUAL Prop 4.3")
    print("#   subcase (2), to exercise step-3 matching directly. (These are")
    print("#   NOT a Jacobian pair, so the full pipeline would reject them at")
    print("#   step 1 before reaching step 3 -- hence the direct call.)")
    print("##################################################################")
    ref = parse_prop43_reference()
    demoP = X**8 * Y**16 + X**8 * Y**14 + X + 1          # corners (8,16),(8,14),(1,0),(0,0)
    demoQ = X**12 * Y**24 + X**12 * Y**21 + X**2 * Y + 1  # corners (12,24),(12,21),(2,1),(0,0)
    cP, cQ = newton_corners(demoP), newton_corners(demoQ)
    print(f"  demo P corners: {fmt_pts(cP)}")
    print(f"  demo Q corners: {fmt_pts(cQ)}")
    match = (cP == ref[2]["P"] and cQ == ref[2]["Q"])
    print(f"  parsed subcase (2) N(P): {fmt_pts(ref[2]['P'])}")
    print(f"  parsed subcase (2) N(Q): {fmt_pts(ref[2]['Q'])}")
    print(f"  -> exact match to subcase (2): {match}")
    assert match, "shape-matcher demo failed -- parser/hull regression"

    # Assertions so the self-test doubles as a regression check.
    assert ra.out_of_scope and not ra.refuted, "self-test (a) verdict changed"
    assert rb.refuted, "self-test (b) verdict changed"
    print("\nAll self-test assertions passed.")


# --------------------------------------------------------------------------- #
def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Exact necessary-condition checker for a claimed plane "
                    "Jacobian (JC(2)) counterexample.")
    ap.add_argument("candidate", nargs="?",
                    help="path to a .py (defines P, Q) or .json ({'P':..,'Q':..}) file")
    ap.add_argument("--self-test", action="store_true",
                    help="run built-in self-tests and exit")
    args = ap.parse_args(argv)

    if args.self_test:
        self_test()
        return 0
    if not args.candidate:
        ap.print_help()
        return 2

    P, Q = load_candidate(args.candidate)
    r = validate(P, Q, label=os.path.basename(args.candidate))
    # Exit non-zero when the claim is refuted, so scripts can gate on it.
    return 1 if r.refuted else 0


if __name__ == "__main__":
    sys.exit(main())
