#!/usr/bin/env python3
"""integral_atlas.py  --  arithmetic atlas of the 51 integral (q_window = 1) families

RESCOPE (2026-07-24): these are 51 FORMAL CHAIN-CENSUS ROWS with q_window = 1,
NOT 51 realized integral window systems.  The q_window identity itself is valid,
but the census supplies `t` from GGV5's final-corner `l`, which is a
Laurent/ramification index rather than the Laurent-CHART exponent (established in
CORNER_RESOLVENT.md sec.5.1).  Only 3 of the 39 distinct integral corner shapes
satisfy the chart precondition b_0 = t(a_0 - 1), so for the large majority the
geometric reading of (t, q) does not hold.  The home (72,108) case is unaffected
-- its chart and divisor data come from GGHV's explicit reduction, not from the
census.

The q_window theorem showed (72,108) is not the unique row with q_window = 1:
there are 51 such formal rows across 23 corner shapes.  Nobody has looked at what
those families actually ARE arithmetically.  This module builds that atlas, using
the marked-polynomial generator (`marked_polynomial.py`) reconstructed from
GGHV's forcing ODE.

For each integral family it reports:
  * the chain corner A_0 = (a_0, b_0) and the census triple (t, kappa, q), dg;
  * the corner integer C = q(kappa+1) - t;
  * the CHART PRECONDITION  b_0 == t*(a_0 - 1)  -- derived from GGHV's three
    worked reductions (see CORNER_RESOLVENT.md sec.5).  Where it fails, the
    standard Laurent chart does not apply and the census (t, q) should not be
    read as (chart exponent, ord_y C);
  * where the generator applies (dg == 1 and t | a_0), the marked polynomial, its
    degree, discriminant squarefree part, Galois group, and whether the corner
    law disc ~ (-1)^(n/2) * C holds.

Read-only.  Writes `integral_atlas.json`.

Usage:
    python integral_atlas.py            # full atlas
    python integral_atlas.py --quiet    # census only, exit 0 iff no law violation
"""
from __future__ import annotations

import argparse
import json
import os

import sympy as sp

import q_window_theorem as Q
from q_window_theorem import corner_integer, q_window
from marked_polynomial import marked, squarefree_part

HERE = os.path.dirname(os.path.abspath(__file__)) or "."
OUT_PATH = os.path.join(HERE, "integral_atlas.json")


def galois_label(poly):
    """Best-effort Galois group name; None if sympy cannot decide."""
    try:
        grp = sp.polys.numberfields.galoisgroups.galois_group(poly)
        return str(grp[0]) if isinstance(grp, tuple) else str(grp)
    except Exception:
        return None


def build(verbose=True):
    data = json.load(open(os.path.join(HERE, "chain_survey_data.json"), encoding="utf-8"))
    distinct = Q._distinct_families(data["families_at_max_M"])

    rows = []
    for r in distinct:
        t, kap, q = r["t"], r["kappa"], r["q"]
        a, b = r["m0"], r["n0"]
        M, H, g, qw = q_window(t, kap, q, a, b)
        if qw != 1:
            continue                              # integral locus only
        A0 = r["A0"]
        a0 = r["a0"]
        b0 = A0[2] if isinstance(A0, (list, tuple)) and len(A0) >= 3 else None
        C = corner_integer(t, kap, q)
        precond = (b0 is not None and b0 == t * (a0 - 1))
        row = {"a0": a0, "b0": b0, "t": t, "kappa": kap, "q": q, "dg": r.get("dg"),
               "base": [a, b], "M": M, "H": H, "C": C,
               "chart_precondition": precond, "motif": r.get("motif")}

        # the generator needs dg == 1 (so q = a0 - 1) and c | a
        if precond and r.get("dg") == 1 and a0 % t == 0:
            poly = marked(a0, t)
            if poly is not None:
                n = poly.degree()
                disc = sp.discriminant(poly.as_expr(), sp.Symbol("y"))
                lhs = squarefree_part(disc)
                rhs = squarefree_part((-1)**(n // 2) * C) if n % 2 == 0 else None
                row.update({"marked_poly": str(poly.as_expr()), "deg_g": n,
                            "disc_sqfree": str(lhs), "law_predicted": str(rhs),
                            "law_holds": (lhs == rhs) if rhs is not None else None,
                            "galois": galois_label(poly)})
        rows.append(row)

    rows.sort(key=lambda r: (r["a0"], r["t"], r["q"]))
    n_pre = sum(1 for r in rows if r["chart_precondition"])
    generated = [r for r in rows if "marked_poly" in r]
    violations = [r for r in generated if r.get("law_holds") is False]

    if verbose:
        print("=" * 96)
        print("INTEGRAL ATLAS  --  the q_window = 1 families")
        print("=" * 96)
        print("  integral families: %d   chart-precondition b0 == t(a0-1) holds: %d   "
              "marked polys generated: %d" % (len(rows), n_pre, len(generated)))
        print()
        print("  %-4s %-5s %-3s %-4s %-3s %-3s %-5s %-6s %-9s %s"
              % ("a0", "b0", "t", "kap", "q", "dg", "C", "chart?", "deg g", "disc_sqfree / law"))
        for r in rows:
            extra = ""
            if "marked_poly" in r:
                extra = "%-9s %s %s" % (r["deg_g"], r["disc_sqfree"],
                                        "OK" if r["law_holds"] else "VIOLATION")
                if r.get("galois"):
                    extra += "  [%s]" % r["galois"][:24]
            print("  %-4s %-5s %-3s %-4s %-3s %-3s %-5s %-6s %s"
                  % (r["a0"], r["b0"], r["t"], r["kappa"], r["q"], r["dg"], r["C"],
                     "yes" if r["chart_precondition"] else "NO", extra))
        print()
        failing = [r for r in rows if not r["chart_precondition"]]
        if failing:
            print("  CHART-PRECONDITION FAILURES (census (t,q) not readable as "
                  "(chart exponent, ord_y C)):")
            for r in failing:
                print("     a0=%-3s b0=%-4s t=%-3s -> b0/(a0-1) = %s  [%s]"
                      % (r["a0"], r["b0"], r["t"],
                         sp.Rational(r["b0"], r["a0"] - 1) if r["a0"] != 1 else "n/a",
                         r["motif"]))
        print()
        print("  corner-law check over generated polynomials: %d tested, %d violations"
              % (len(generated), len(violations)))

    out = {"schema": "d2-integral-atlas-v1", "n_integral": len(rows),
           "n_chart_precondition": n_pre, "n_generated": len(generated),
           "n_law_violations": len(violations), "rows": rows}
    with open(OUT_PATH, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2, sort_keys=True)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()
    out = build(verbose=not args.quiet)
    print("\nwrote", os.path.basename(OUT_PATH))
    if out["n_law_violations"]:
        print("CORNER-LAW VIOLATIONS: %d" % out["n_law_violations"])
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
