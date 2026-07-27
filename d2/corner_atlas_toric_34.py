#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""corner_atlas_toric_34.py -- does the toric identity transfer to (8,28)/(3,4)/144?

THE QUESTION.  `corner_atlas.py` establishes that among GGV5's 34 possible
counterexamples with max deg <= 150, the corner A_0 = (8,28) is the ONLY one that
passes BOTH preconditions of the (72,108) toric mechanism:

    (P1)  the retraction shape  b0 == t*(a0-1)   (so C is not a monomial), and
    (P2)  the chart exponent  t = 4              (the unique solution of the
          weight condition (t+1) | (4t+9), toric_general.py B4).

Two of the 34 rows sit on that corner: `(8,28)/(11\\4,7)/(3,2)/108` -- which is
(72,108), the case this campaign closed -- and `(8,28)/(7\\4,3)/(3,4)/144`.
Since `toric_general.py` proved the weight arithmetic depends on `t` ALONE, the
exponent k = (4t+9)/(t+1) = 5 is weight-admissible at (3,4) too.  So the second
row is the sharpest available test of whether the toric identity `6*W*Z = e^5` is
a fact about the CORNER (t = 4, retraction) or about `(m,n) = (2,3)` specifically.

WHAT IS TESTED.  Both halves of the mechanism, exactly as `toric_general.py`
defines them for (75,125):

  T1  PRODUCT-OF-MINORS.  Is some Q-linear combination of products of two 2x2
      window-Hankel minors congruent to e^k modulo the Phi-free G-generators,
      for ANY weight-admissible k?  (`product_search`, exact rank over Q.)
  T2  THE SHARPER SHAPE -- "the cofactor of Z is W, itself a minor".  Is e^k in
      I + J with the toric cofactors restricted to the WINDOW coefficients
      (`toric_ideal_search` tag "win")?  This is the shape that turns the
      relation into a PRODUCT and yields the divisor law and the Pi^4 contact
      order.  T2 is the condition the atlas brief singles out.
  T3  the looser invariant (tag "all", cofactors may involve the state
      variables d_i) -- reported for calibration only.

DISCIPLINE.  Every negative is paired with a POSITIVE CONTROL run through the
identical code path at (2,3,4), where the identity is known to hold; a negative
that is not accompanied by a passing control is a broken search, not a result.
The admissible-exponent set is printed so that "no hit" cannot mean "nothing was
tried".

Run:
    python corner_atlas_toric_34.py            # ~3 min
    python corner_atlas_toric_34.py --quiet    # exit 0 iff all checks pass
"""

from __future__ import annotations

import argparse
import json
import sys
import time

import toric_general as tg

T = 4                      # the unique toric-admissible chart exponent
TARGET = (3, 4)            # (8,28)/(7\4,3)/(3,4)/144   -- GGV5 tex:1831
CONTROL = (2, 3)           # (8,28)/(11\4,7)/(3,2)/108  -- GGV5 tex:1832, CLOSED
KS = [3, 4, 5, 6, 7]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()
    rows = []

    def ck(name, ok, detail=""):
        rows.append((name, bool(ok), detail))
        if not args.quiet:
            print("  [%s] %-56s %s" % ("PASS" if ok else "FAIL", name, detail))
        return bool(ok)

    t0 = time.time()
    out = {}
    for tag, (a, b) in (("target", TARGET), ("control", CONTROL)):
        S = tg.system(a, b, T)
        mn, mw = tg.hankel_minors(S, 0)
        prods = sorted({mw[i] + mw[j] for i in mn for j in mn})
        adm = sorted({p // (T + 1) for p in prods if p % (T + 1) == 0})
        prod = tg.product_search(S, extra_state=0)
        tor = tg.toric_ideal_search(S, KS)
        out[tag] = dict(
            mn=(a, b), n_minors=len(mn),
            minor_weight_range=[min(mw.values()), max(mw.values())],
            admissible_k=adm,
            product_search={k: dict(n_pairs=v[0], dim_I_W=v[1], hit=v[2])
                            for k, v in sorted(prod.items())},
            toric_win={k: dict(n_cols=tor[k]["win"][0], hit=tor[k]["win"][1])
                       for k in KS},
            toric_all={k: dict(n_cols=tor[k]["all"][0], hit=tor[k]["all"][1])
                       for k in KS},
        )
        if not args.quiet:
            print("\n  (m,n) = (%d,%d), t = %d:  %d window-Hankel minors, "
                  "u-weights %d..%d" % (a, b, T, len(mn), min(mw.values()),
                                        max(mw.values())))
            print("    weight-admissible exponents k: %s" % adm)
            print("    T1 product-of-minors : %s"
                  % {k: ("HIT" if v[2] else "no") + " (%d pairs, dim I_W=%d)"
                     % (v[0], v[1]) for k, v in sorted(prod.items())})
            print("    T2 toric 'win'       : %s"
                  % {k: ("IN" if tor[k]["win"][1] else "no") for k in KS})
            print("    T3 toric 'all'       : %s"
                  % {k: ("IN" if tor[k]["all"][1] else "no") for k in KS})
    if not args.quiet:
        print()

    tgt, ctl = out["target"], out["control"]

    # -- the searches were real ------------------------------------------------
    ck("V0 the exponent k = 5 IS weight-admissible at (3,4,4)",
       5 in tgt["admissible_k"],
       "admissible k = %s -- the (72,108) exponent (4t+9)/(t+1) = 5 is among "
       "them, so the prediction gets a fair test" % tgt["admissible_k"])
    ck("V0b the search had candidates at every admissible k",
       all(v["n_pairs"] > 0 for v in tgt["product_search"].values()),
       "minor pairs per k: %s"
       % {k: v["n_pairs"] for k, v in tgt["product_search"].items()})

    # -- T1 --------------------------------------------------------------------
    ck("V1 (3,4,4) admits NO product-of-minors identity, at ANY admissible k",
       not any(v["hit"] for v in tgt["product_search"].values()),
       "exact rank over Q at k = %s: all miss"
       % sorted(tgt["product_search"]))
    ck("V1-CTRL positive control: the same routine FINDS (2,3,4) at k = 5",
       [k for k, v in ctl["product_search"].items() if v["hit"]] == [5],
       "so V1's negative is not a broken search")

    # -- T2, the condition the brief singles out -------------------------------
    ck("V2 *** (3,4,4) has NO window-coefficient toric relation e^k, k <= 7 -- "
       "the cofactor of Z is never itself a minor",
       not any(tgt["toric_win"][k]["hit"] for k in KS),
       "cofactor columns tried per k: %s"
       % {k: tgt["toric_win"][k]["n_cols"] for k in KS})
    ck("V2-CTRL positive control: (2,3,4) DOES have it, minimally at k = 5",
       ctl["toric_win"][5]["hit"]
       and not any(ctl["toric_win"][k]["hit"] for k in (3, 4)),
       "this is exactly the (72,108) identity 6*W*Z = e^5, whose cofactor W is "
       "itself a minor")

    # -- T3, calibration -------------------------------------------------------
    tmin = min([k for k in KS if tgt["toric_all"][k]["hit"]] or [None])
    cmin = min([k for k in KS if ctl["toric_all"][k]["hit"]] or [None])
    ck("V3 the LOOSER invariant is present but at a different exponent",
       tmin == 5 and cmin == 4 and tmin != cmin,
       "minimal k with e^k in I+J allowing state cofactors: %s at (3,4) vs %s "
       "at (2,3).  toric_general C5 records 6 at (3,5).  So the loose "
       "invariant is 4 / 5 / 6 at (2,3) / (3,4) / (3,5) -- it MOVES, and it is "
       "not the product mechanism." % (tmin, cmin))
    ck("V3-OBS the 4/5/6 pattern is reported as an OBSERVATION, not a law",
       True,
       "three data points fit max(m,n)+1 (4,5,6) and also fit m+n-1 only at "
       "(2,3).  Per toric_general D2's warning about one-point fits, this is "
       "NOT asserted as a formula; it is logged for the next corner to test.")

    # -- the conclusion --------------------------------------------------------
    ck("V4 *** VERDICT: the toric identity does NOT transfer to (8,28)/(3,4)",
       (not any(v["hit"] for v in tgt["product_search"].values()))
       and (not any(tgt["toric_win"][k]["hit"] for k in KS)),
       "(8,28)/(3,4)/144 shares BOTH preconditions with (72,108) -- same "
       "corner, same retraction, same t = 4, same admissible exponent 5 -- and "
       "still carries neither half of the mechanism.  So t = 4 and the "
       "retraction shape are NECESSARY but NOT SUFFICIENT: (m,n) = (2,3) is "
       "doing the work.")

    if not args.quiet:
        with open("corner_atlas_toric_34.json", "w") as fh:
            json.dump(dict(schema="corner_atlas_toric/1", t=T, ks=KS,
                           seconds=round(time.time() - t0, 1), data=out),
                      fh, indent=1)
        print("\n  wrote corner_atlas_toric_34.json   [%.1fs]"
              % (time.time() - t0))

    bad = [r for r in rows if not r[1]]
    print("\n%s  corner_atlas_toric_34: %d/%d checks pass"
          % ("FAIL" if bad else "OK  ", len(rows) - len(bad), len(rows)))
    for n, _, d in bad:
        print("   FAILED: %s   %s" % (n, d))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
