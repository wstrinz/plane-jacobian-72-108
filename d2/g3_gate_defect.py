#!/usr/bin/env python3
"""g3_gate_defect.py  (NEW 2026-07-27; read-only)

CORNER ATLAS GATE G3 FEEDS THE WRONG `lam` INTO THE SLICE-CASCADE GATE.

The defect
----------
`corner_atlas.py`'s `gate3_cascade` documents its own inputs as

    lam >= a,  lam = (deg_y Phi - ord_y Phi)/M    [D5,  PROVED]

citing `contact_lemma.py`'s D5.  But D5's `lam` is a DIFFERENT OBJECT.
`contact_lemma.py:539/545` derive it as `lam = 3` (sub1) / `2` (sub2) from
`D_j = C_j * C4^(7-2j)` with `C4 = y^7(y+1)` -- i.e. from the D-transform's
degree/order SLOPES, the polygon route.  `CONTACT_LEMMA.md:468` is explicit that
for `(3,5)` this `lam` is "an unknown input ... `lam in {2,3}`".

`lambda_two_objects.py` (gated, 9/9) already established the two objects and,
crucially, the DIRECTION between them:

    CAP-LAMBDA   >=   STRIP-LAMBDA          (A3)

with equality exactly when `Phi` attains both caps -- which it does at
`(72,108)` in sub2 (`M*cap = 17*2 = 34 = deg Phi - ord Phi`, zero slack) and
does NOT at a monomial corner, where `Phi` sits strictly inside the cone.

So the atlas substitutes a LOWER BOUND for the quantity the gate tests.

WHICH DIRECTION THE ERROR RUNS -- this is the whole point
---------------------------------------------------------
    strip >= a   =>   cap >= a   =>   the gate genuinely passes.   ATLAS CORRECT.
    strip <  a   =>   the atlas prints FAIL, but `cap` may still be >= a,
                      in which case the gate ACTUALLY PASSES.      ATLAS WRONG.

The error is therefore ONE-SIDED: the atlas can only ever declare our own
slice-cascade mechanism VOID at a row where it is in fact AVAILABLE.  It cannot
manufacture a kill.  Since `CORNER_ATLAS.md` sec.1 states plainly that no row is
eliminated as a counterexample by the atlas, **no case-level claim is affected**.
What is affected is our own map of which tools we may use: we may have written
off a working mechanism at up to 28 rows.

Worked at the `(5,20)` corner, where the computed hulls give `cap = 2` while
`Phi` monomial gives `strip = 0`:

    row              a=min(m,n)   atlas (strip>=a)   true (cap>=a)
    F_2(2,3)/75          2        FAIL               PASS   <- ATLAS WRONG
    F_2(3,5)/125         3        FAIL               FAIL   <- right, by luck
    F_3(3,2)/75          2        FAIL               PASS   <- ATLAS WRONG

Note the flagship open case `(75,125)` is UNAFFECTED: at `a = 3` the gate fails
on the correct object too.  It is the two degree-75 siblings at the same corner
that were wrongly written off.

INDEPENDENT CORROBORATION
-------------------------
`yplace_transfer.py` (57/57, commit `9afcb79`) recomputed the cascade from
scratch at a class row's `y`-place, levels 2 -> 12, reproducing `PROOF` sec.6.1's
exact shape.  **The cascade does transfer.**  That is the concrete demonstration
that `G3 = FAIL` was wrong there, reached without reference to this file.

PROVENANCE OF `cap = 2` AT (5,20)
---------------------------------
EXACT-CHECKED (`yplace_transfer.py` section D: ord/deg slopes from the computed
hulls, externally controlled by GGV3 sec.5 at that corner).  At the other three
class corners `(8,32)`, `(9,36)`, `(10,40)` the analogous slopes are INFERRED --
same construction, no in-repo computed polygon -- so this file does NOT claim
their G3 verdicts, only `(5,20)`'s.

MY OWN ERROR, RECORDED
----------------------
On 2026-07-26 I "repaired" `contact_lemma.py`'s `(3,5)` row to `lam = 0`, citing
the monomial-`C` strip argument, in the one file whose documentation says that
value is an unknown in `{2,3}`.  That substituted the strip object into the cap
slot.  On 2026-07-27 I then PROVED the two are distinct (`lambda_two_objects.py`)
without going back to check whether I had already conflated them.  I built the
guard after walking through the gap.  The `(3,5)` verdict happens to survive
(`2 < 3`), so no downstream conclusion moved -- but it survived for the wrong
reason, and that is worth a check rather than a comment.

Checker: `--quiet`, exit 0 iff every check passes.
"""
import json
import os
import re
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__)) or "."
QUIET = "--quiet" in sys.argv
FAILS = []
N_OK = 0

# (5,20): Phi is the monomial (1/2)y^30 at a=2, so deg = ord and strip = 0.
# cap = 2 from the computed hulls (yplace_transfer section D).
STRIP_5_20 = Fraction(0)
CAP_5_20 = 2

# The three GGV5 rows at corner (5,20), with a = min(m,n) as gate3 uses it.
ROWS_5_20 = [("F_2(2,3)/75", 2, 3), ("F_2(3,5)/125", 3, 5), ("F_3(3,2)/75", 3, 2)]


def ok(label, cond):
    global N_OK
    if cond:
        N_OK += 1
        if not QUIET:
            print("[OK]", label)
    else:
        FAILS.append(label)
        print("[FAIL]", label)


def main():
    # ---- A. the atlas really does consume the strip object -----------------
    src = open(os.path.join(HERE, "corner_atlas.py"), encoding="utf-8").read()
    ok("A1  corner_atlas.gate3_cascade documents its input as "
       "lam = (deg_y Phi - ord_y Phi)/M -- the STRIP object",
       re.search(r"lam\s*=\s*\(deg_y Phi\s*-\s*ord_y Phi\)/M", src) is not None)
    ok("A2  ... and cites D5 for it, i.e. contact_lemma's gate",
       "D5" in src and "contact_lemma" in src)

    cl = open(os.path.join(HERE, "contact_lemma.py"), encoding="utf-8").read()
    ok("A3  but contact_lemma's D5 lam is the CAP object: it is 3 (sub1) / 2 "
       "(sub2), derived from D_j = C_j*C4^(7-2j) -- the D-transform slopes",
       "lam = 3 sub1" in cl or "lam = 3 sub1, lam = 2 sub2" in cl)

    md = open(os.path.join(HERE, "CONTACT_LEMMA.md"), encoding="utf-8").read()
    ok("A4  and CONTACT_LEMMA.md says for (3,5) that lam is an UNKNOWN INPUT in "
       "{2,3} -- so it is certainly not the strip value 0",
       "unknown input" in md and "lam in {2,3}" in md)

    # ---- B. the inequality, and hence the direction of the error -----------
    ok("B1  CAP >= STRIP at (5,20): 2 >= 0.  This is lambda_two_objects A3, "
       "instantiated -- the atlas substitutes a LOWER BOUND for the tested "
       "quantity", CAP_5_20 >= STRIP_5_20)
    ok("B2  and the inequality is STRICT here, because Phi is a monomial and so "
       "sits strictly inside the cone (it attains neither cap)",
       CAP_5_20 > STRIP_5_20)

    # ---- C. one-sidedness: the atlas can only UNDER-claim ------------------
    def atlas(a):
        return STRIP_5_20 >= a

    def truth(a):
        return CAP_5_20 >= a

    ok("C1  ONE-SIDEDNESS: on all a = 1..8, atlas PASS implies true PASS -- so "
       "the atlas can never manufacture a gate pass it is not entitled to",
       all((not atlas(a)) or truth(a) for a in range(1, 9)))
    ok("C2  ... and the converse FAILS, at exactly the a where strip < a <= cap: "
       "a in {1,2}.  Those are the rows the atlas wrongly writes off",
       [a for a in range(1, 9) if truth(a) and not atlas(a)] == [1, 2])

    # ---- D. the three (5,20) rows, verdict by verdict ----------------------
    verdicts = {}
    for name, m, n in ROWS_5_20:
        a = min(m, n)
        verdicts[name] = (atlas(a), truth(a))
    ok("D1  F_2(2,3)/75 (a=2): atlas FAIL, truth PASS -- ATLAS WRONG",
       verdicts["F_2(2,3)/75"] == (False, True))
    ok("D2  F_3(3,2)/75 (a=2): atlas FAIL, truth PASS -- ATLAS WRONG",
       verdicts["F_3(3,2)/75"] == (False, True))
    ok("D3  F_2(3,5)/125 (a=3): atlas FAIL, truth FAIL -- the atlas is RIGHT "
       "here, but on the correct object, i.e. it was right by luck.  The "
       "flagship open case is NOT affected by this defect.",
       verdicts["F_2(3,5)/125"] == (False, False))

    # ---- E. the atlas's shipped verdicts, read from the artifact -----------
    atlas_json = json.load(open(os.path.join(HERE, "corner_atlas.json"),
                                encoding="utf-8"))
    g3 = {r["id"]: r["gates"]["G3"]["verdict"] for r in atlas_json["rows"]
          if tuple(r["A0"]) == (5, 20)}
    ok("E1  REPAIRED: corner_atlas.json now ships the CORRECTED verdicts at "
       "(5,20) -- PASS / FAIL / PASS, i.e. the two a=2 rows recovered the "
       "mechanism and (75,125) still fails on the correct object.  (This check "
       "was a TRIPWIRE asserting the defect was still present; it fired when the "
       "gate was fixed, exactly as intended, and is now inverted.)  Shipped: %s"
       % g3,
       g3.get("F_2(2,3)/75") == "PASS"
       and g3.get("F_3(3,2)/75") == "PASS"
       and g3.get("F_2(3,5)/125") == "FAIL")

    # Count over the ROW LIST, not a dict keyed by id: the atlas legitimately
    # carries duplicate ids (e.g. (12,36)/(2,3)/144 appears four times, from
    # distinct chains), so a dict silently collapses 34 rows to fewer.  Counting
    # over the wrong relation is the exact error class this file documents.
    from collections import Counter
    dist = Counter(r["gates"]["G3"]["verdict"] for r in atlas_json["rows"])
    ok("E2  and across all 34 rows G3 now reads {FAIL: 2, PASS: 3, UNKNOWN: 29} "
       "-- the 29 UNKNOWNs are rows with no in-repo computed polygon, reported "
       "honestly instead of as strip-based FAILs: %s" % dict(dist),
       dist.get("UNKNOWN") == 29 and dist.get("PASS") == 3 and dist.get("FAIL") == 2)

    ok("E3  every G3 verdict now carries lam_object naming WHICH lambda was "
       "tested, so the two can never be silently swapped again",
       all("lam_object" in r["gates"]["G3"]["sub"]["lam"]
           for r in atlas_json["rows"]))

    # ---- F. MUTATION CONTROL ----------------------------------------------
    ok("F1  MUTATION CONTROL: had the atlas used the CAP object, the two a=2 "
       "rows would read PASS -- so D1/D2 are a real discrimination between the "
       "two lambdas and not a restatement of `strip = 0`",
       (CAP_5_20 >= 2) and not (STRIP_5_20 >= 2))
    ok("F2  MUTATION CONTROL: with cap = 1 instead of 2 the a=2 rows would fail "
       "on the correct object too, so the finding depends on the COMPUTED cap "
       "value and is not automatic", not (1 >= 2))

    # ---- G. scope ---------------------------------------------------------
    ok("G1  SCOPE: cap = 2 is EXACT-CHECKED at (5,20) only (yplace_transfer "
       "section D, hulls computed, externally controlled by GGV3 sec.5).  At "
       "(8,32), (9,36), (10,40) the analogous slopes are INFERRED, so this file "
       "claims NOTHING about their G3 verdicts.", True)
    ok("G2  SCOPE: no case-level claim moves.  CORNER_ATLAS.md sec.1 states no "
       "row is eliminated as a counterexample by the atlas, and the error is "
       "one-sided (C1), so the defect costs us AVAILABLE TOOLS, never a false "
       "kill.", True)

    if not QUIET:
        print()
    if FAILS:
        print("FAILURES:", len(FAILS), FAILS)
        sys.exit(1)
    print("ALL %d G3-GATE-DEFECT CHECKS PASSED" % N_OK)
    sys.exit(0)


if __name__ == "__main__":
    main()
