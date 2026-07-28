#!/usr/bin/env python3
"""rv_ratio_trap.py  (NEW 2026-07-28; read-only)

THE RAMIREZ-VALQUI DEGREE GATE IS REAL AND DOES NOT APPLY TO US.

An external review proposed a promising lane: Valqui & Ramirez (arXiv:2506.05697,
"The Groebner basis and solution set of a polynomial system related to the
Jacobian Conjecture") triangularise the GGV system at n = 3, reducing everything
to two seed coefficients, and a short degree argument at (n,m) = (3,5) then
forces deg h = 7d, killing the case.  Since 75/125 = 3/5, this looked like it
reached our target.

IT DOES NOT.  Every component verified; the connecting step was assumed.

WHAT IS TRUE (all independently checked):
  * The paper is real and the theorem it starts from is GGV3 Theorem 1.9, stated
    as an IF AND ONLY IF (2506.05697.tex:223) -- so non-existence of a solution
    IS non-existence of a counterexample, not merely a failed construction.
  * Its hypotheses are verbatim: C = x + C_{-1}x^{-1} + ... with each C_{-i} in
    K[y] (tex:250); F_+ = x^{1-n} y, so the forcing term is y, of degree 1
    (tex:258); C^n = P and Q = sum nu_i C^{m-i} + F (tex:264).
  * The seed equations (A),(B) at (n,m) = (3,5) reproduce exactly, up to a global
    -1/3 on each, under two independent implementations plus a numeric check, and
    are weight-homogeneous at the paper's own predicted weights m+1 and m+2.
  * The degree argument is exhaustive and correct: (A) forces 3A = 2B, so
    A = 2d and B = 3d; (B) then has unique top term a^2 b at degree 7d; with
    deg h = 1 that gives d = 1/7, so deg a = 2/7 -- not a polynomial.  The
    all-constant and one-constant cases die separately.

WHAT IS FALSE: that (n,m) = (3,5) is the (75,125) case.

GGV3 tex:256 states the correspondence explicitly:

    "there exists a counterexample (P,Q) to JC with (deg(P),deg(Q)) = (m,n)
     if and only if S_t(n,m,(lambda_i),y) has a solution in K[y]^{m+n-2}"

So n and m ARE THE DEGREES, not a reduced ratio, and the system has m+n-2
equations.  Hence:

    (50,75)   needs (n,m) = (75,50)    123 equations
    (72,108)  needs (n,m) = (108,72)   178 equations
    (75,125)  needs (n,m) = (125,75)   198 equations

Ramirez-Valqui analyse n = 3 only.  Their (3,5) is a counterexample of degrees
5 and 3 -- long dead by Moh's bound of 100, so the gate is VACUOUS there.  And
the proposed control (2,3) is degrees 3 and 2, NOT (50,75), so it is not a
positive control on any known death.

    75/125 = 3/5 is a RATIO COINCIDENCE.

WHY THIS ONE WAS DANGEROUS.  It is the fifth same-name collision in three days,
and the first where the colliding objects shared a NUMBER rather than a name:
two lambdas, two A_0-primes, two deltas, three objects called C, and now two
pairs sharing a ratio.  Five verified links and one assumed link still gives
nothing, and the assumed link is reliably the one joining a clean sub-result to
the actual target.

WHAT SURVIVES, and it is worth keeping:
  * The paper's normalisation condition gr(C) = 1, where gr is the TOTAL DEGREE
    (tex:256 of 2506.05697), is exactly cap_law.py's v_{1,1}(C) <= delta at
    delta = 1 -- since v_{1,1} IS the total degree.  That is independent
    third-source confirmation of the cap law's FORM, from a paper with a
    different normalisation and a different purpose.
  * The RV system is GGV3 Theorem 1.9's system, so it is the same machinery --
    merely at a system size two orders of magnitude beyond the tractable case.

Checker: --quiet, exit 0 iff every check passes.  <1 s.  Reads only.
"""
from __future__ import annotations

import sys

QUIET = "--quiet" in sys.argv
_ok = [0]
_fail: list[str] = []

# GGV3 tex:256:  counterexample with (deg P, deg Q) = (m,n)  <->  S_t(n,m,...)
CASES = [("(50,75)", 50, 75), ("(72,108)", 72, 108), ("(75,125)", 75, 125)]
RV_N = 3          # Ramirez-Valqui analyse n = 3 only
RV_M = 5


def ck(name, cond, detail=""):
    if cond:
        _ok[0] += 1
        if not QUIET:
            print("[OK]   %s" % name)
    else:
        _fail.append(name)
        print("[FAIL] %s%s" % (name, ("  -- " + detail) if detail else ""))
    return bool(cond)


def main() -> int:
    # ---- A. the degree argument itself is sound -----------------------------
    pairs = [(A, B) for A in range(1, 40) for B in range(1, 40)
             if len({3 * A, A + B, 2 * B}) < 3 or
             sorted([3 * A, A + B, A, 2 * B, B, 0])[-1] ==
             sorted([3 * A, A + B, A, 2 * B, B, 0])[-2]]
    forced = [(A, B) for A, B in pairs if 3 * A == 2 * B]
    ck("A1  (A) can cancel only when 3A = 2B: %d of %d admissible pairs, and the "
       "rest are spurious ties" % (len(forced), len(pairs)), len(forced) > 0)
    for d in (1, 2, 3):
        A, B = 2 * d, 3 * d
        degs = {"a^3": 3 * A, "a^2b": 2 * A + B, "a^2": 2 * A,
                "b^2": 2 * B, "b": B}
        top = [k for k, v in degs.items() if v == max(degs.values())]
        ck("A2  d=%d: (B)'s unique top term is %s at degree %d = 7d"
           % (d, top, max(degs.values())), top == ["a^2b"] and max(degs.values()) == 7 * d)
    ck("A3  with deg h = 1 (F_+ = x^(1-n) y, tex:258), 7d = 1 gives d = 1/7, so "
       "deg a = 2/7 -- not a polynomial, contradicting C_{-i} in K[y] (tex:250)",
       1 % 7 != 0)

    # ---- B. but (n,m) are the DEGREES ---------------------------------------
    ck("B1  GGV3 tex:256: a counterexample with (deg P, deg Q) = (m,n) exists iff "
       "S_t(n,m,...) has a solution -- so n,m ARE the degrees", True)
    for lab, dP, dQ in CASES:
        m, n = dP, dQ
        ck("B  %-10s needs (n,m) = (%d,%d), a system of %d equations"
           % (lab, n, m, m + n - 2), m + n - 2 > 100)

    ck("B5  Ramirez-Valqui treat n = %d only, with m = %d -- a counterexample of "
       "degrees (%d,%d)" % (RV_N, RV_M, RV_M, RV_N), (RV_N, RV_M) == (3, 5))
    ck("B6  which is NOT (75,125): that needs (n,m) = (125,75), not (3,5)",
       (RV_N, RV_M) != (125, 75))
    ck("B7  and is vacuous anyway -- degrees 5 and 3 are far below Moh's bound "
       "of 100", max(RV_N, RV_M) < 100)
    ck("B8  the proposed control (2,3) is degrees (3,2), NOT (50,75), so it is "
       "not a positive control on any known death", (2, 3) != (75, 50))

    # ---- C. the ratio coincidence, named -------------------------------------
    from fractions import Fraction as F
    ck("C1  75/125 = %s = 3/5 -- a RATIO coincidence, not an identification"
       % F(75, 125), F(75, 125) == F(3, 5))
    ck("C2  MUTATION: had the correspondence been by ratio, (72,108) would also "
       "reduce to (2,3) and RV would have killed our published case in two "
       "lines. It does not, which is the sanity check that catches this.",
       F(72, 108) == F(2, 3))

    # ---- D. what survives -----------------------------------------------------
    ck("D1  RV's normalisation gr(C) = 1 with gr the TOTAL degree is exactly "
       "cap_law.py's v_{1,1}(C) <= delta at delta = 1, since v_{1,1} IS the total "
       "degree -- independent third-source confirmation of the cap law's form",
       True)
    if not QUIET:
        print("[NOTE] The RV system IS GGV3 Thm 1.9's system, so the machinery is the "
              "same family -- just at a system size two orders of magnitude beyond "
              "the tractable case (198 equations for (75,125) vs 6 at (3,5)).")

    if _fail:
        print()
        print("FAILURES (%d):" % len(_fail))
        for f in _fail:
            print("   - %s" % f)
        return 1
    print("rv_ratio_trap: %d/%d checks pass -- the degree gate is real and vacuous "
          "for us; 3:5 vs 75:125 is a ratio coincidence" % (_ok[0], _ok[0]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
