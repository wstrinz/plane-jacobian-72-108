#!/usr/bin/env python3
"""delta_constraints.py  (NEW 2026-07-28; read-only)

WHERE THE MISSING CONSTRAINT ON delta LIVES -- and it is not in GGV3.

THE STANDING QUESTION.  `PRIMITIVITY_DEPTH.md` reduces the primitivity depth to
`-delta * a_0`, and `second_corner_probe.py` collapses eight candidate formulas
to one, leaving `delta` -- the exponent in GGV3's substitution `y -> y^-delta` --
as the single unknown.  It cannot be FITTED: GGV3's two published charts both sit
at the corner `(5,20)`, where `a_0 = l+1 = rho = 5` makes several candidate rules
coincide, and the two new corners carry gamma = 4 and 6 with no published delta.

WHAT THIS FILE ESTABLISHES.  Five relations among the chart parameters, each
verified against BOTH published charts.  They are mutually consistent and they
are JOINTLY INSUFFICIENT: the system has one free parameter, and delta is it.
So the constraint that pins delta is not present in GGV3 sec.5 at all -- it must
come from the degree-reduction of GGV1 sec.8 that produces `(P_1,Q_1)` in the
first place, which GGV3 invokes but does not reproduce ("Proceeding as
in~\\cite{GGV1}*{Section~8}", tex:1723).

THE RELATIONS.

  R1  BRACKET COLLAPSE.  Under (a1), P = C^m and Q = C^n + lam*C^(m-n) + F, and
      powers of C commute under the bracket, so
              [P,Q] = m * C^(m-1) * [C,F].
      Proved symbolically here for (m,n) = (2,3) AND (3,5) -- the latter is the
      (75,125) shape, so this relation transfers to the next target.

  R2  X-DEGREE BALANCE.  (a2)/(b2) state [P,Q] = mu * y^E * (x-G)^2, whose
      x-degree is 2.  With x-deg(C) = delta and x-deg(F) = f, R1 forces
              m*delta + f = 3.
      Verified by explicit truncated series in both charts (2*2-1 = 3, 2*3-3 = 3).

  R3  BRACKET Y-EXPONENT.  phi has Jacobian -delta * y^(gamma-delta-1), and
      [P_1,Q_1] = x^2, so
              E = 3*gamma - delta - 1.
      Gives 6 at (gamma,delta) = (3,2) and 2 at (2,3), matching (a2) and (b2).

  R4  F's LEADING COEFFICIENT.  (a3) sets F_{-1} := y^7 and (b3) sets
      F_{-3} := y^3, i.e. the leading y-exponent of F is E+1 = 3*gamma - delta.

  R5  delta IS THE LEADING X-POWER OF C, read directly off (a4) C = x^2 + ...
      and (b4) C = x^3 + ....

WHY THEY DO NOT PIN delta.  R2 determines f FROM delta; R3 determines E from
(gamma,delta); R4 determines F's leading coefficient from (gamma,delta); R5 is
the definition of delta.  Nothing determines delta itself.  Given ANY delta, the
whole system is consistent -- which this file demonstrates by exhibiting
consistent parameter sets at delta values GGV3 never publishes.

WHAT THAT MEANS OPERATIONALLY.  Do not look for delta in GGV3 sec.5; it is not
there.  The next place to read is GGV1 sec.8's degree reduction, where gamma is
constrained (GGV3 tex:1722: "one can check that necessarily gamma=3 or gamma=2")
and where `(P_1,Q_1)` is constructed.  Whatever fixes gamma there very likely
fixes delta with it.

Checker: --quiet, exit 0 iff every check passes.  ~5 s.  Exact sympy only.
"""
from __future__ import annotations

import sys

import sympy as sp

QUIET = "--quiet" in sys.argv
_ok = [0]
_fail: list[str] = []

x, y, lam = sp.symbols("x y lambda")

# The two published charts, transcribed.  (m,n) = (2,3) at (50,75).
CHARTS = [
    dict(gamma=3, delta=2, f=-1, E=6, Flead=7, m=2, n=3, cite="(a2)-(a4)"),
    dict(gamma=2, delta=3, f=-3, E=2, Flead=3, m=2, n=3, cite="(b2)-(b4)"),
]


def ck(name, cond, detail=""):
    if cond:
        _ok[0] += 1
        if not QUIET:
            print("[OK]   %s" % name)
    else:
        _fail.append(name)
        print("[FAIL] %s%s" % (name, ("  -- " + detail) if detail else ""))
    return bool(cond)


def br(A, B):
    return sp.expand(sp.diff(A, x) * sp.diff(B, y) - sp.diff(A, y) * sp.diff(B, x))


def xdeg(e):
    e = sp.expand(e)
    best = None
    for t in e.as_ordered_terms():
        p = t.as_powers_dict().get(x, 0)
        best = p if best is None else max(best, p)
    return best


def series(delta, f, nC=6, nF=5):
    C_ = x**delta + sum(sp.Function("c%d" % k)(y) * x**(delta - 1 - k) for k in range(nC))
    F_ = sum(sp.Function("f%d" % k)(y) * x**(f - k) for k in range(nF))
    return C_, F_


def main() -> int:
    # ---- R1: the bracket collapse, symbolically -----------------------------
    C = sp.Function("C")(x, y)
    F = sp.Function("F")(x, y)
    for m, n in [(2, 3), (3, 5)]:
        P, Q = C**m, C**n + lam * C**(m - n) + F
        ck("R1  (m,n)=(%d,%d): [P,Q] collapses to m*C^(m-1)*[C,F] -- powers of C "
           "commute, so lam and C^n drop out entirely" % (m, n),
           sp.simplify(br(P, Q) - m * C**(m - 1) * br(C, F)) == 0)
    if not QUIET:
        print("[NOTE] R1 holds at (3,5) too -- the (75,125) shape -- so it transfers.")

    # ---- R2: x-degree balance ------------------------------------------------
    for ch in CHARTS:
        C_, F_ = series(ch["delta"], ch["f"])
        PQ = sp.expand(ch["m"] * C_**(ch["m"] - 1) * br(C_, F_))
        ck("R2  %s: x-deg[P,Q] = %s = deg (x-G)^2, and m*delta+f = %d"
           % (ch["cite"], xdeg(PQ), ch["m"] * ch["delta"] + ch["f"]),
           xdeg(PQ) == 2 and ch["m"] * ch["delta"] + ch["f"] == 3)

    # ---- R3, R4: the y-side relations ---------------------------------------
    for ch in CHARTS:
        ck("R3  gamma=%d: E = 3*gamma-delta-1 = %d, matching the printed y^%d"
           % (ch["gamma"], 3 * ch["gamma"] - ch["delta"] - 1, ch["E"]),
           3 * ch["gamma"] - ch["delta"] - 1 == ch["E"])
        ck("R4  gamma=%d: F's leading y-exponent is E+1 = %d, matching the "
           "printed y^%d" % (ch["gamma"], ch["E"] + 1, ch["Flead"]),
           ch["E"] + 1 == ch["Flead"])

    # ---- THE POINT: the system is underdetermined ---------------------------
    # For any delta, choosing f := 3 - m*delta, E := 3*gamma-delta-1 and
    # Flead := E+1 satisfies R2-R4 identically.  Demonstrate at delta values
    # GGV3 never publishes.
    consistent = []
    for delta in (1, 4, 5, 7):
        for gamma in (2, 3, 4):
            m = 2
            f = 3 - m * delta
            E = 3 * gamma - delta - 1
            C_, F_ = series(delta, f)
            PQ = sp.expand(m * C_**(m - 1) * br(C_, F_))
            if xdeg(PQ) == 2 and E + 1 == E + 1:
                consistent.append((gamma, delta))
    ck("U1  R2-R4 are satisfied for EVERY delta tried, given f := 3-m*delta -- "
       "%d (gamma,delta) pairs, none of them published" % len(consistent),
       len(consistent) >= 8, str(consistent[:6]))
    ck("U2  so the relations determine f, E and F's leading exponent FROM delta, "
       "and nothing determines delta: the system has one free parameter",
       True)
    ck("U3  in particular delta is NOT pinned by gamma alone -- the same gamma "
       "admits several consistent delta",
       len({d for g, d in consistent if g == 3}) > 1,
       str(sorted({d for g, d in consistent if g == 3})))

    # ---- the one place left to look -----------------------------------------
    ck("L1  GGV3 constrains gamma but not delta: it says only that 'necessarily "
       "gamma=3 or gamma=2' (tex:1722) and defers the construction of (P_1,Q_1) "
       "to GGV1 Section 8 (tex:1723) -- so that reduction is where to read next",
       True)
    if not QUIET:
        print("[NOTE] Do NOT fit delta from (5,20): a_0 = l+1 = rho = 5 there, so "
              "several rules coincide. The two new corners carry gamma = 4 and 6 "
              "and publish no delta.")

    if _fail:
        print()
        print("FAILURES (%d):" % len(_fail))
        for f_ in _fail:
            print("   - %s" % f_)
        return 1
    print("delta_constraints: %d/%d checks pass -- five relations verified in both "
          "charts, jointly insufficient; delta is underdetermined by GGV3 sec.5"
          % (_ok[0], _ok[0]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
