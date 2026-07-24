#!/usr/bin/env python3
"""Regression ladder for the Phase D max-plus infinity layer.

Gates cascade_engine.py's infinity place (v_inf = -deg) against the three
independently documented by-hand degree arguments, in the order prescribed
by CASCADE_ENGINE_PLAN.md ("Phase D implementation design"):

  R0  max-plus semantics unit checks (unique max forces; ties drop only
      with recorded obligations; descend cases (a)/(b)/(c) mirror the
      audited min-plus descent on a hand-worked level-6 example);
  R1  the a_t=9 T2 uniform infinity kill (T5_90_T2.md section 2);
  R2  the 43/50 constant-cell degree kills of T5_90_T1.md section 3,
      with the surviving cells exactly d=2 x {sigma=0, deg sigma<=5};
  R3  the T2-column degree dominations of T5_T2_COLUMN.md section 2
      (all seven strict margins, e.g. max(T0..T5)=234 < deg(T6)=236 at
      a6 b1100) and the documented pattern-A/B open states;
  R4  joint q+t+inf smoke on the a9 b1000 T2 cell: the infinity place
      may only remove survivor cases relative to the audited
      cascade_cones_qt.json, never add them, and the g5=0 flag case is
      narrowed to the single degree state (deg e, deg sigma) = (10, 8)
      (T5_T2_COLUMN.md section 4, G5 row).

The checks drive the engine's own tables (source-parsed from
f31_graded.txt); no h_l coefficient is entered here by hand.
"""

from __future__ import annotations

import json
from pathlib import Path

import cascade_engine as ce

def _require(_cond, _msg):
    """Proof-critical check: fails loudly and exits nonzero, unaffected by python -O."""
    if not _cond:
        import sys as _sys
        print("FAIL: " + str(_msg))
        _sys.exit(1)


NEG_INF = ce.NEG_INF
ROOT = Path(__file__).resolve().parent


def r0_semantics() -> None:
    # Unique maximum: T1 level 6 at (deg d2, deg d1, deg sigma, deg e)
    # = (1,2,3,9): monomial degrees 2X+K=5, X+E=11, 2Z=6 -- unique max 11.
    degstate = (1.0, 2.0, 3.0, 9.0)
    flags = (False, False, False)
    options = ce.deg_h_options(6, degstate, flags)
    _require(options == [(11.0, ())], options)

    # Descend level 6 (a=9: v=3, deg ehat=0) below deg g7 = 32:
    # target 35; case (b) gives deg g6 in [0,34], case (c) the depth-0 tie
    # at 35; case (a) is empty because the h-term equals the target.
    results = ce.descend_options_inf(6, 32.0, False, degstate, flags, 37, 0, 3)
    values = sorted(value for value, _ in results)
    _require(values == [float(k) for k in range(36)], values)
    _require(all(not obligations for _, obligations in results), "all(not obligations for _, obligations in results)")

    # Tie: (1,2,6,10) makes X+E = 2Z = 12 with 2X+K = 5 below: the maximum
    # may drop only with a degree_tie_drop obligation naming both tied
    # monomials, or vanish identically.
    degstate = (1.0, 2.0, 6.0, 10.0)
    options = ce.deg_h_options(6, degstate, flags)
    _require(options[0] == (12.0, ()), "options[0] == (12.0, ())")
    drops = [entry for entry in options if entry[0] not in (12.0, NEG_INF)]
    _require(len(drops) == 12, len(drops))
    for value, obligations in drops:
        (obligation,) = obligations
        _require(obligation.kind == "degree_tie_drop", "obligation.kind == \"degree_tie_drop\"")
        _require(obligation.depth == int(12 - value), "obligation.depth == int(12 - value)")
        _require(len(obligation.tied) == 2, "len(obligation.tied) == 2")
    _require(options[-1][0] == NEG_INF, "options[-1][0] == NEG_INF")
    _require(options[-1][1][0].kind == "identical_vanishing", "options[-1][1][0].kind == \"identical_vanishing\"")

    # Unique maxima cannot vanish or drop even when requested.
    _require(ce.deg_h_options(7, (0.0, 3.0, 0.0, 9.0), flags, required=5.0) == [], "ce.deg_h_options(7, (0.0, 3.0, 0.0, 9.0), flags, required=5.0) == []")
    _require(ce.deg_h_options(7, (0.0, 3.0, 0.0, 9.0), flags, required=NEG_INF)
        == [], "ce.deg_h_options(7, (0.0, 3.0, 0.0, 9.0), flags, required=NEG_INF) == []")
    print("R0 max-plus semantics: unique max forces, ties drop with obligations")


def _t2_state_survives(a: int, deg_e: int, deg_sigma: int) -> bool:
    """Any (deg d2, zero-flag) infinity chain for a T2 degree state."""

    for d2_zero in (False, True):
        k_domain = (NEG_INF,) if d2_zero else tuple(float(k) for k in range(5))
        for k_deg in k_domain:
            for mask in range(4):
                g_zero = {6: False, 5: bool(mask & 1), 4: bool(mask & 2)}
                if ce.inf_place_profiles(
                    a,
                    "T2",
                    10 + 3 * a,
                    4,
                    False,
                    d2_zero,
                    g_zero,
                    (k_deg, NEG_INF, float(deg_sigma), float(deg_e)),
                ):
                    return True
    return False


def r1_a9_t2() -> None:
    # T5_90_T2.md section 2: deg e = 9 (constant E) has no consistent
    # degree chain for ANY (deg sigma, deg d2, flags).
    for deg_sigma in range(0, 9):
        _require(not _t2_state_survives(9, 9, deg_sigma), deg_sigma)
    # deg e = 10 (linear E) must survive tropically -- section 1's UFD
    # argument is a residue argument, not a degree argument -- and every
    # witness must carry a cancellation obligation.
    surviving = [z for z in range(0, 9) if _t2_state_survives(9, 10, z)]
    _require(surviving, "linear-E states must remain tropically open")
    profile = ce.inf_place_profiles(
        9, "T2", 37, 4, False, False, {6: False, 5: False, 4: False},
        (0.0, NEG_INF, 0.0, 10.0),
    )[0]
    kinds = {obligation.kind for obligation in profile.obligations}
    _require(kinds & {"leading_cancellation", "degree_tie_drop"}, kinds)
    print("R1 a_t=9 T2 kill (T5_90_T2.md): constant-E dead, linear-E open")


def r2_a9_t1_constant() -> None:
    # T5_90_T1.md section 3: 50 constant-E cells (deg d1 = d in 0..4;
    # sigma == 0 or deg sigma = z in 0..8); 43 die on degrees alone and
    # the seven survivors are exactly d=2 with sigma=0 or z<=5.
    def cell_survives(d: int, sigma_zero: bool, z: int) -> bool:
        for d2_zero in (False, True):
            k_domain = (
                (NEG_INF,) if d2_zero else tuple(float(k) for k in range(5))
            )
            for k_deg in k_domain:
                for mask in range(8):
                    g_zero = {
                        7: False,
                        6: bool(mask & 1),
                        5: bool(mask & 2),
                        4: bool(mask & 4),
                    }
                    z_deg = NEG_INF if sigma_zero else float(z)
                    if ce.inf_place_profiles(
                        9, "T1", 37, 4, sigma_zero, d2_zero, g_zero,
                        (k_deg, float(d), z_deg, 9.0),
                    ):
                        return True
        return False

    survivors = set()
    for d in range(0, 5):
        for case in ["zero"] + list(range(9)):
            sigma_zero = case == "zero"
            z = 0 if sigma_zero else int(case)
            if cell_survives(d, sigma_zero, z):
                survivors.add((d, case))
    expected = {(2, "zero")} | {(2, z) for z in range(6)}
    _require(survivors == expected, survivors)
    print("R2 a_t=9 T1 constant cells (T5_90_T1.md): 43/50 dead, 7 ties open")


def r3_t2_column() -> None:
    # T5_T2_COLUMN.md section 2: the seven strict margins, as
    # (a, deg e, deg sigma) via deg e = a+B+f, deg sigma = S+z.
    kill_rows = [
        (5, 6, 2),   # a5 b1000 (0,0):  218 < 226
        (6, 7, 2),   # a6 b1000 (0,0):  226 < 229
        (6, 7, 3),   # a6 b1000 (0,1):  226 < 231
        (6, 8, 5),   # a6 b1000 (1,3):  234 < 238
        (6, 8, 4),   # a6 b1100 (0,0):  234 < 236
        (6, 8, 5),   # a6 b1100 (0,1):  234 < 238
        (6, 9, 6),   # a6 b1110 (0,0):  242 < 243
    ]
    for a, deg_e, deg_sigma in kill_rows:
        _require(not _t2_state_survives(a, deg_e, deg_sigma), (a, deg_e, deg_sigma))
    # Documented open residual states (patterns A and B) stay open.
    open_rows = [
        (9, 10, 2), (9, 10, 8),          # R9 endpoints
        (8, 8, 3), (8, 10, 5), (8, 10, 8),  # R80 pattern A + B
        (7, 8, 3), (7, 10, 7),           # R71 pattern A + B
    ]
    for a, deg_e, deg_sigma in open_rows:
        _require(_t2_state_survives(a, deg_e, deg_sigma), (a, deg_e, deg_sigma))
    print("R3 T2-column dominations (T5_T2_COLUMN.md): 7 margins dead, A/B open")


def r4_join_smoke() -> None:
    # Joint q+t+inf on the a9 b1000 T2 cell against the audited artifact.
    reference = json.loads(
        (ROOT / "cascade_cones_qt.json").read_text(encoding="utf-8")
    )
    row = next(
        record
        for record in reference["branches"]
        if record["a_t"] == 9
        and record["b"] == [1, 0, 0, 0]
        and record["branch"] == "T2"
    )
    _require(row["status"] == "survives", "row[\"status\"] == \"survives\"")
    reference_cases = {
        (case["sigma_zero"], case["d2_zero"], tuple(case["g_zero_levels"]))
        for case in row["survivor_cases"]
    }

    outcome = ce.analyze_branch(
        9, (1, 0, 0, 0), "T2", 4, include_t=True, include_inf=True
    )
    inf_cases = {
        (case["sigma_zero"], case["d2_zero"], tuple(case["g_zero_levels"]))
        for case in outcome["survivor_cases"]
    }
    _require(inf_cases <= reference_cases, inf_cases - reference_cases)

    # Every survivor witness carries an infinity record with a full chain.
    for case in outcome["survivor_cases"]:
        inf_records = [
            record
            for record in case["witness"]
            if record["place"] == "inf"
        ]
        _require(len(inf_records) == 1, "len(inf_records) == 1")
        _require(set(inf_records[0]["deg_g"]) == {str(l) for l in range(1, 7)}, "set(inf_records[0][\"deg_g\"]) == {str(l) for l in range(1, 7)}")

    # The g5=0 flag case is narrowed to (deg e, deg sigma) = (10, 8)
    # (T5_T2_COLUMN.md section 4, G5 row: state (0,6,12;10,8)).
    g5_cases = [
        case
        for case in outcome["survivor_cases"]
        if tuple(case["g_zero_levels"]) == (5,)
    ]
    for case in g5_cases:
        (record,) = [r for r in case["witness"] if r["place"] == "inf"]
        _require(record["deg_e"] == 10, record)
        _require(record["deg_sigma"] == 8, record)
    print(
        "R4 joint q+t+inf smoke (a9 b1000 T2): no new survivors vs "
        "cascade_cones_qt.json; g5=0 narrowed to (deg e, deg sigma)=(10,8)"
    )


def r5_t2_squeeze() -> None:
    # T5_T2_COLUMN.md final ledger: with the proven level-5 squeeze F^2|G
    # (C24) joined to the infinity layer, exactly the four documented
    # cells die and the eight documented cells stay open.
    column = {
        (5, (1, 0, 0, 0)): "dead",
        (6, (1, 0, 0, 0)): "dead",
        (6, (1, 1, 0, 0)): "dead",
        (6, (1, 1, 1, 0)): "dead",
        (7, (1, 0, 0, 0)): "open",
        (7, (1, 1, 0, 0)): "open",
        (7, (1, 1, 1, 0)): "open",
        (7, (3, 0, 0, 0)): "open",
        (8, (0, 0, 0, 0)): "open",
        (8, (1, 0, 0, 0)): "open",
        (8, (1, 1, 0, 0)): "open",
        (9, (1, 0, 0, 0)): "open",
    }
    for (a, b_vector), expected in column.items():
        outcome = ce.analyze_branch(
            a, b_vector, "T2", 4,
            include_t=True, include_inf=True, t2_squeeze=True,
        )
        actual = "open" if outcome["status"] == "survives" else "dead"
        _require(actual == expected, (a, b_vector, actual))
    print(
        "R5 T2 squeeze joined to infinity: the four T5_T2_COLUMN cells "
        "die, the eight open cells survive"
    )


def main() -> None:
    r0_semantics()
    r1_a9_t2()
    r2_a9_t1_constant()
    r3_t2_column()
    r4_join_smoke()
    r5_t2_squeeze()
    print("cascade infinity layer: PASS")


if __name__ == "__main__":
    main()
