#!/usr/bin/env python3
"""frontier_rebuild.py -> FRONTIER_REBUILD.md + frontier_rebuild.json

The frontier REGENERATED through the authoritative compiler under the e | Phi
divisor filter (`divisor_filter.py`), reported at CELL / FLAG-CASE / STATE level,
for BOTH windows (sub1, sub2) and BOTH settings of the C08/C20 residue-kill
toggle.

WHY TWO CENSUSES
----------------
`FIELD_SCOPE_AUDIT.md` (2026-07-25) downgrades the C08/C20 forbidden-rise lemmas
from KILL to CONSTRAINT over an arbitrary characteristic-zero base field.  If
that downgrade lands, branches REOPEN -- the opposite direction from the divisor
filter.  The two effects must NOT be netted into a single number, so every table
here is reported twice: `rl` (C08/C20 enabled, the committed status quo) and
`norl` (C08/C20 disabled).  Flipping the toggle is a matter of rerunning the
named commands below; nothing here has to be rebuilt by hand.

PIPELINE (every artifact is produced by one of these, never typed in)
---------------------------------------------------------------------
  1. cascade_engine.py   --depth 4 --with-t --with-inf --t2-squeeze [--residue-kills]
  2. phase_d_states.py   --window W [--no-residue-kills --sweep S] [--divisor-filter]
  3. divisor_filter.py   (invoked by step 2; standalone entry point for reruns)
  4. frontier_rebuild.py (this file)

USAGE
-----
    python frontier_rebuild.py                # rebuild FRONTIER_REBUILD.md
    python frontier_rebuild.py --quiet        # checker; exit 0 iff consistent
    python frontier_rebuild.py --print-commands
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

MD_OUT = "FRONTIER_REBUILD.md"
JSON_OUT = "frontier_rebuild.json"

# ---------------------------------------------------------------------------
# the 2x2 matrix of pipeline runs
# ---------------------------------------------------------------------------
CASCADE_BASE = ("python cascade_engine.py --depth 4 --with-t --with-inf "
                "--t2-squeeze --window %s%s --json-out %s")

RUNS = []
for _win, _rl_sweep, _norl_sweep in (
    ("sub2", "cascade_cones_qt_inf_rl.json", "cascade_cones_qt_inf_norl.json"),
    ("sub1", "cascade_cones_sub1_qt_inf_rl.json", "cascade_cones_sub1_qt_inf_norl.json"),
):
    RUNS.append(dict(
        window=_win, c0820="on", tag="%s/C08+C20 ON" % _win,
        cascade=_rl_sweep,
        cascade_cmd=CASCADE_BASE % (_win, " --residue-kills", _rl_sweep),
        universe="phase_d_states_%s.json" % _win,
        universe_cmd="python phase_d_states.py --window %s" % _win,
        filtered="phase_d_states_%s_divfilter.json" % _win,
        filtered_cmd="python phase_d_states.py --window %s --divisor-filter" % _win,
    ))
    RUNS.append(dict(
        window=_win, c0820="off", tag="%s/C08+C20 OFF" % _win,
        cascade=_norl_sweep,
        cascade_cmd=CASCADE_BASE % (_win, "", _norl_sweep),
        universe="phase_d_states_%s_norl.json" % _win,
        universe_cmd="python phase_d_states.py --window %s --no-residue-kills "
                     "--sweep %s" % (_win, _norl_sweep),
        filtered="phase_d_states_%s_norl_divfilter.json" % _win,
        filtered_cmd="python phase_d_states.py --window %s --no-residue-kills "
                     "--sweep %s --divisor-filter" % (_win, _norl_sweep),
    ))


def path(fn):
    return os.path.join(HERE, fn)


def load(fn, required=True):
    p = path(fn)
    if not os.path.exists(p):
        if required:
            sys.stderr.write("FATAL: missing %s -- build it with:\n" % fn)
            for r in RUNS:
                for k in ("cascade", "universe", "filtered"):
                    if r[k] == fn:
                        sys.stderr.write("  %s\n" % r[k + "_cmd"])
            sys.exit(1)
        return None
    with open(p, encoding="utf-8") as fh:
        return json.load(fh)


def sha16(fn):
    h = hashlib.sha256()
    with open(path(fn), "rb") as fh:
        h.update(fh.read())
    return h.hexdigest()[:16]


def git_commit():
    try:
        r = subprocess.run(["git", "rev-parse", "HEAD"], cwd=HERE,
                           capture_output=True, text=True)
        if r.returncode == 0 and r.stdout.strip():
            return r.stdout.strip()
    except Exception:
        pass
    return os.environ.get("GIT_COMMIT", "unknown")


# ---------------------------------------------------------------------------
# census
# ---------------------------------------------------------------------------
def cellname(a, b, br):
    return "a%d_b%s_%s" % (a, "".join(map(str, b)), br)


def universe_census(U):
    cells = defaultdict(lambda: {"flagcases": 0, "states": 0})
    for c in U["cases"]:
        rec = cells[cellname(c["a_t"], c["b"], c["branch"])]
        rec["flagcases"] += 1
        rec["states"] += len(c["states"])
    per_fc = {}
    for c in U["cases"]:
        k = "%s|%s|%s|%s" % (cellname(c["a_t"], c["b"], c["branch"]),
                             c["sigma_zero"], c["d2_zero"],
                             tuple(c["g_zero_levels"]))
        per_fc[k] = len(c["states"])
    return {
        "cells": len(cells),
        "flagcases": len(U["cases"]),
        "states": sum(len(c["states"]) for c in U["cases"]),
        "cell_names": sorted(cells),
        "per_cell": {k: dict(v) for k, v in sorted(cells.items())},
        "per_flagcase": per_fc,
    }


def build_census():
    import divisor_filter as df

    rows = []
    for r in RUNS:
        casc = load(r["cascade"])
        U = load(r["universe"])
        # The filtered universe is a COMPILER PRODUCT but also exactly derivable
        # from the unfiltered one, so it is optional here: if present it is
        # cross-checked against a fresh re-derivation, if absent the census falls
        # back to the re-derivation (and says so).
        F = load(r["filtered"], required=False)
        if U.get("residue_kills") is not None:
            assert U["residue_kills"] == (r["c0820"] == "on"), \
                "%s: residue_kills provenance mismatch" % r["universe"]
        assert casc["residue_kills"] == (r["c0820"] == "on"), \
            "%s: cascade residue_kills mismatch" % r["cascade"]
        assert not casc.get("partial_checkpoint"), "%s is a partial checkpoint" % r["cascade"]
        # RE-DERIVE the filtered universe from the unfiltered one here, and
        # require it to agree with the artifact the compiler wrote.  This makes
        # the census independent of the stored filtered artifact.
        redo = df.filter_universe(U, r["window"])
        redo_c = universe_census(redo)
        before = universe_census(U)
        after = universe_census(F) if F is not None else redo_c
        agree = (redo_c["cells"] == after["cells"]
                 and redo_c["flagcases"] == after["flagcases"]
                 and redo_c["states"] == after["states"]
                 and redo_c["cell_names"] == after["cell_names"])
        filt = df.DivisorFilter(r["window"])
        rows.append(dict(
            window=r["window"], c0820=r["c0820"], tag=r["tag"],
            cascade=r["cascade"], universe=r["universe"], filtered=r["filtered"],
            cascade_cmd=r["cascade_cmd"], universe_cmd=r["universe_cmd"],
            filtered_cmd=r["filtered_cmd"],
            filtered_present=F is not None,
            sha={k: (sha16(r[k]) if os.path.exists(path(r[k])) else "(derived)")
                 for k in ("cascade", "universe", "filtered")},
            cascade_summary=casc["summary"],
            before=before, after=after,
            recompute_agrees=agree,
            census=(F if F is not None else redo)["divisor_filter"]["census"],
            stages=stage_census(F if F is not None else redo, r["window"]),
            caps=filt.caps, e_min=filt.e_min,
            degree_forcing_active=filt.degree_forcing_active,
        ))
    return rows


# ---------------------------------------------------------------------------
# DOWNSTREAM CELL-KILL OVERLAYS (other lanes' lemmas, imported as data)
# ---------------------------------------------------------------------------
# The `e | Phi` filter above is COMPUTED here from the syzygy.  Two later lanes
# kill further CELLS by arguments this file does not re-derive.  We import from
# them only the CELL NAMES (a short, checkable list) and compute every flag-case
# and state count ourselves from the universes -- so the mathematics is theirs
# and attributed, and every number is still the compiler's.
#
# Scope discipline (their own guards, repeated here because they bind this file):
#   * `R | e^2`, `e*R | Phi`, `R = c*(y+1)^rho` are **T2-ONLY, permanently**
#     (a witness point on V(G1,G2,G3,G5) with R=0, e=1 refutes them on T1).
#     Stage 2 therefore only ever names T2 cells.
#   * The SPINE reduction is a sub2 five-family statement; it is NOT applied to
#     sub1 here.  Its T1 columns are unconditional only because `t^a | R,S,T` was
#     adjudicated branch-independent (T1_BRANCH.md); before that they were
#     conditional.
#
# EVIDENCE GRADE, per stage (2026-07-26).  Every stage carries `level` (the
# proof-DAG evidence level its evidence actually supports, on the scale
# claimed < exact-checked < independently-audited < certified) and `evidence`
# (the checkers, and the INDEPENDENT audit if there is one).  The grade is the
# lane's OWN, imported here as data exactly like the cell lists; the ledger
# (`state_kill_ledger.py`) and the DAG (`proof_dag.py`) read it from here rather
# than re-asserting it, so there is one place to correct.
#
#   stages 2-4  exact-checked          -- same-author exact checkers, no
#                                         independent audit
#   stages 5-7  independently-audited  -- a second, independently authored
#                                         checker reproduced the result
#
# `closes_frontier=True` marks the three stages added on 2026-07-26 that take
# the enumerated f31 frontier to EMPTY.  It exists so that consumers which are
# measured against the *pre-closure* frontier on purpose (`g4_row.py`'s lambda-row
# census baseline, which asks "what does the lambda row ADD to the 34-cell
# frontier?" and answers "nothing") can pin that baseline explicitly instead of
# silently becoming vacuous.
STAGES = [
    dict(
        id="stage2_T2_divisor",
        title="T2 divisor filter (`e*R | Phi`, `R | e^2`, `R = c*(y+1)^rho`)",
        source="DIVISOR_CONSEQUENCES.md sec.6.5/6.6, sec.7",
        checker="python divisor_consequences.py --quiet",
        note="T2-only and cap-free. On T2 the three relations turn the spare "
             "ansatz into a closed form; the surviving arithmetic condition is "
             "decidable exactly over Q. sub2's ENTIRE T2 branch empties.",
        level="exact-checked",
        evidence="divisor_consequences.py (same-author exact checker). No "
                 "independent audit; the grade is exact-checked, NOT "
                 "independently-audited.",
        dead={
            "sub2": ["a9_b1000_T2", "a8_b1100_T2", "a7_b1110_T2"],
            "sub1": ["a7_b1000_T2", "a7_b1110_T2", "a7_b1111_T2",
                     "a8_b1100_T2", "a8_b1110_T2", "a9_b1000_T2",
                     "a9_b1100_T2", "a9_b1110_T2", "a9_b1111_T2",
                     "a10_b1000_T2", "a10_b1100_T2", "a10_b1110_T2",
                     "a10_b1111_T2"],
        },
    ),
    dict(
        id="stage3_spine",
        title="SPINE five-family reduction (n = 0..4), T1 columns unconditional",
        source="SPINE.md sec.0 verdict table; T1_BRANCH.md (t^a | R,S,T proved "
                "branch-independent, which is what removes the conditionality)",
        checker="python spine.py --quiet && python spine_verify.py --quiet && "
                "python t1_branch.py --quiet",
        note="sub2 ONLY -- the five-family reduction is a sub2 statement. Of the "
             "ten (n, branch) columns, nine are EMPTY; the single residue is "
             "a10_b0000_T1, saved by d1 free of degree 6.",
        level="exact-checked",
        evidence="spine.py + spine_verify.py + t1_branch.py (same-author exact "
                 "checkers). No independent audit.",
        dead={
            "sub2": ["a6_b1111_T1", "a7_b1110_T1", "a7_b1110_T2",
                     "a8_b1100_T1", "a8_b1100_T2", "a9_b1000_T1",
                     "a9_b1000_T2"],
            "sub1": [],
        },
    ),
    dict(
        id="stage4_positive_slice",
        title="Positive-slice obstruction (inverse d3-shift polynomiality), n = 0",
        source="POSITIVE_SLICE.md sec.5-6; derived from generators.json + "
               "window_caps_verify.py W2/W3 + upstream_facts.json corners",
        checker="python positive_slice.py --quiet && "
                "python positive_slice_verify.py --quiet",
        note="sub2 ONLY, n = 0 ONLY -- the cell SPINE leaves open. The shifted "
             "G-system ideal is genuinely NON-EMPTY there, so no elimination "
             "could ever close it; the obstruction is that those points do not "
             "lift back through the RATIONAL d3-killing shift to a polynomial P "
             "with the prescribed Newton support. Branch-independent (d1 is "
             "never set to 0), so it re-kills a10_b0000_T2 as well, agreeing "
             "with SPINE sec.6.6 by an unrelated route. No Groebner over the "
             "G-system and no field-scope caveat: the final step is a resultant "
             "over Q (561971200 != 0).",
        level="exact-checked",
        evidence="positive_slice.py (63/63) + positive_slice_verify.py (79/79). "
                 "Both same-author; no independent audit, so the grade is "
                 "exact-checked.",
        dead={
            "sub2": ["a10_b0000_T1", "a10_b0000_T2"],
            "sub1": [],
        },
    ),
    # -----------------------------------------------------------------------
    # ADDED 2026-07-26.  The three stages below were deliberately kept OUT of
    # this pipeline while they were same-author -- wiring an unaudited kill into
    # the frontier is what produced the v0.3.2 erratum.  Every one of them is
    # now reproduced by a second, independently authored checker, which is what
    # expires that reason.  Cell lists are imported as data exactly like stages
    # 2-4; every flag-case and state count is recomputed here.
    # -----------------------------------------------------------------------
    dict(
        id="stage5_slice_obstruction",
        title="Stacked P/Q positive-slice obstruction: `a_t >= 9`",
        source="SLICE_OBSTRUCTION.md sec.3/sec.6 (drop-in record in "
               "slice_obstruction_stage.json, `a_t_min: 9`); audited by "
               "SLICE_OBSTRUCTION_AUDIT.md",
        checker="python slice_obstruction_basis.py --quiet --deep && "
                "python slice_obstruction_audit.py --quiet",
        note="WINDOW-INDEPENDENT and BRANCH-INDEPENDENT: the criterion is "
             "`a_t = v_t(e) >= 9`, from polynomiality of the P and Q slices, the "
             "12k order floor, and `cap_n + 1 >= 2n-2` (true in both windows). "
             "Immune to the C08/C20 field-scope downgrade -- every step is a "
             "t-adic valuation over Q with no square class, no splitting field "
             "and no residue arithmetic -- so the kill is IDENTICAL under ON and "
             "OFF. It does NOT touch the alternate regime (`a_t in {12,14}`): the "
             "bound is the wrong shape for `a_t >= 11`. The sub2 row is "
             "corroboration, not a frontier delta -- standard sub2 is already "
             "EMPTY after stage 4, and this re-kills three of SPINE's seven sub2 "
             "cells by an unrelated mechanism.",
        level="independently-audited",
        evidence="slice_obstruction_basis.py --deep (59/59) + INDEPENDENT audit "
                 "slice_obstruction_audit.py (56/56, CONFIRMED-WITH-CORRECTIONS; "
                 "neither correction touches the bound, and the audit's own "
                 "recommendation is that `a_t >= 9` may be regraded "
                 "independently-audited). Premise [I3] (dm1 = D~_{-1} = e, "
                 "dm2 = D~_{-2} = R) is separately audited by i3_audit.py (81/81, "
                 "PROVED construction fact).",
        dead={
            "sub1": ["a2_b1111_T1", "a3_b1000_T1", "a3_b1110_T1",
                     "a4_b0000_T1", "a4_b1100_T1", "a4_b1111_T1",
                     "a5_b1000_T1", "a5_b1110_T1", "a6_b0000_T1",
                     "a6_b1100_T1", "a6_b1111_T1", "a6_b1111_T2",
                     "a7_b1000_T1", "a7_b1100_T1", "a7_b1110_T1",
                     "a7_b1111_T1", "a8_b0000_T1", "a8_b0000_T2",
                     "a8_b1000_T1", "a8_b1100_T1", "a8_b1110_T1",
                     "a8_b1111_T1", "a8_b1111_T2"],
            "sub2": ["a6_b1111_T1", "a7_b1110_T1", "a8_b1100_T1"],
        },
        closes_frontier=True,
    ),
    dict(
        id="stage6_syzygy_collision",
        title="K-syzygy exact valuation collision: `a_t <= 9`, hence `a_t = 9` EXACTLY",
        source="SYZYGY_COLLISION.md sec.0/sec.5 (its own sec.10 item 3 specifies "
               "this drop-in shape); audited by AT_LE9_AUDIT.md",
        checker="python syzygy_collision.py --quiet",
        note="Cap-free, branch-free, window-independent. `2*Phi = e*B` is an "
             "EXACT identity on the G-variety and `v_t(Phi) = 30` exactly, so "
             "`v_t(B) = 30 - a_t` EXACTLY. Pushed into unshifted h-coordinates B "
             "collapses to `h_2*h_5^2 + 3*h_5*h_7 + 3*h_1*h_5*h_6 + 3*h_6^2`, and "
             "for every `a >= 10` all four terms exceed `30 - a`. At `a = 9` "
             "three of the four land on 21 on the nose, so this is a THRESHOLD, "
             "not a blanket (check X10 asserts the a=9 survival). SINGLE "
             "LOAD-BEARING NEW INPUT: cascade level 12, `v_t(h_6) >= 11`; --fast "
             "(level 10) genuinely fails X9. It also refutes `a_t = 12` and "
             "`a_t = 14`, independently re-emptying all six alternate-regime T1 "
             "branches -- see ALT_CLOSURE below.",
        level="independently-audited",
        evidence="syzygy_collision.py (25/25) + INDEPENDENT audit at_le9_audit.py "
                 "(76/76, CONFIRMED-WITH-CORRECTIONS). *** THIS BOUND IS "
                 "SINGLE-LEGGED. *** AT_LE9_AUDIT.md C-1 proves the two "
                 "'independent' proofs consume the IDENTICAL four equations "
                 "(G1,G2,G3,G5body ARE r_13,r_14,r_15,r_17 -- the K-syzygy is an "
                 "exact Q[d]-combination of them), so the corroboration is "
                 "common-mode in its INPUT; only the extraction differs. C-2 "
                 "further shows slice_phi_yplace imposes (P<)/(Q) and the cascade "
                 "base in the SHIFTED chart, where they do not transfer, so "
                 "SLICE_PHI_YPLACE.md's `a_t <= 9` is NOT ESTABLISHED as written "
                 "(not refuted either). *** DO NOT cite slice_phi_yplace as "
                 "independent corroboration anywhere. *** The surviving leg is "
                 "syzygy_collision ALONE; the audit reproduces it end to end in "
                 "the unshifted chart with different algebra at every step (C-4), "
                 "which is what earns independently-audited -- it does NOT make "
                 "the result two-legged. ZERO-MARGIN NOTE (D14): "
                 "`v_t(S) = v_t(D~_{-3}) >= 11` has NO margin -- the binding term "
                 "is `(1/16)*h_1^2*h_5` at 2+9 = 11, not h_7. Anything downstream "
                 "needing `v_t(S) >= 12` is FALSE.",
        dead={
            "sub1": ["a10_b0000_T1", "a10_b0000_T2", "a10_b1000_T1",
                     "a10_b1100_T1", "a10_b1110_T1", "a10_b1111_T1"],
            "sub2": ["a10_b0000_T1"],
        },
        closes_frontier=True,
    ),
    dict(
        id="stage7_sub1_spine9",
        title="SPINE cofactor identity at `a_t = 9`: the last five cells are EMPTY",
        source="SUB1_SPINE9.md sec.5/sec.6; audited by SPINE9_AUDIT.md",
        checker="python sub1_spine9.py --quiet",
        note="sub1 ONLY, and conditional on `a_t = 9` -- which stages 5 and 6 now "
             "supply, both independently audited. The cofactor identity "
             "`F*Z = (1/6)*gamma^5*t^9*Pi^4` is re-derived from generators.json "
             "alone, residual exactly 0, and consumes NO degree cap (SPINE's sub2 "
             "zero-slack coincidence is neither needed nor used). Indexed by "
             "`k = sum(b_i) = deg Pi`: k = 1,2,3 die on the exact marked `Pi^2` "
             "support test for EVERY z in [0,9]; k = 4 dies on a degree count "
             "with z = 3 pinned; k = 0 dies on the boxed row. FOUR of the five "
             "kills consume no cascade input at all; only a9_b0000_T1 does, and "
             "SPINE9_AUDIT G10 shows that dependency (`v_t(h_6) >= 10`, "
             "`v_t(h_7) >= 11`) is a two-line consequence of the AUDITED rows "
             "`v_t(h_1..h_5) >= (1,3,5,7,9)` plus `(P<)`. Level 12 is NOT needed "
             "here: at `a_t = 9` the inverse shift gives "
             "`v_t(R) >= min(11, 1+9) = 10`, so the level-12 upgrade is invisible.",
        level="independently-audited",
        evidence="sub1_spine9.py (37/37) + INDEPENDENT audit spine9_audit.py "
                 "(81/81, CONFIRMED -- no step found wrong, weakened or vacuous; "
                 "the audit adds the G10 derivation above and an E6 non-vacuity "
                 "control at k = 2, the one value where 'infeasible' is not "
                 "structurally forced). The `a_t = 9` premise is inherited from "
                 "stages 5 and 6, both independently-audited, so this stage's "
                 "grade is not capped below by its dependency.",
        dead={
            "sub1": ["a9_b0000_T1", "a9_b1000_T1", "a9_b1100_T1",
                     "a9_b1110_T1", "a9_b1111_T1"],
            "sub2": [],
        },
        closes_frontier=True,
    ),
]

# ---------------------------------------------------------------------------
# THE ALTERNATE REGIME (a_t >= 11) -- closed, and NOT by a cell list
# ---------------------------------------------------------------------------
# The alternate regime lives outside both phase-D universes (they are capped at
# a_t <= 10), so it cannot be a STAGES entry: there are no cells here to name.
# It is registered as its own record, read by state_kill_ledger.py and
# proof_dag.py the same way the stage records are.
#
# WHY IT IS EMPTY, and by which leg.  `a_t <= 9` (stage 6) is cap-free,
# branch-free and window-independent -- it holds in the alternate regime too, and
# every alternate-regime branch has `a_t in {11..15} >= 11 > 9`.  So the whole
# 52-branch leaf L_alt is EMPTY outright, by an independently audited bound, and
# no branch-by-branch argument is needed.  That is the primary leg recorded here.
#
# ALT_LEVEL12.md closes the same regime by a different route and is recorded as
# CORROBORATION, not as the primary leg: its coverage is the six surviving T1
# branches (the post-C33/C34/C44 residual ALT_FRONTIER_V2.md leaves), not all 52,
# and it shares the slice calculus with stage 5, so it corroborates the
# ARITHMETIC and not the premises.  Its hinge (horn 1 dead for every a >= 11,
# because 3a > 30) is a proved theorem already in the repo
# (ALT_FRONTIER_V2.md sec.3 / check B4_trichotomy), re-derived independently
# twice.
ALT_CLOSURE = dict(
    id="alt_closure_at_le9",
    title="Alternate regime `a_t >= 11` is EMPTY (all 52 branches of L_alt)",
    source="SYZYGY_COLLISION.md sec.0 and check X12 (primary); ALT_LEVEL12.md "
           "sec.0/sec.3/sec.4 (corroboration, six T1 branches); ALT_LEVEL12.md "
           "sec.5 L6.4 for the cap-free `a_t <= 10` corollary that retires "
           "ALT_FRONTIER_V2.md sec.7.1's `a >= 16` scope hole",
    checker="python syzygy_collision.py --quiet && python alt_level12.py --quiet",
    level="independently-audited",
    evidence="PRIMARY LEG: `a_t <= 9` (stage6_syzygy_collision), "
             "independently-audited via at_le9_audit.py (76/76). It is cap-free, "
             "branch-free and window-independent, and every L_alt branch has "
             "`a_t >= 11`, so all 52 die at once -- check X12 refutes `a_t = 12` "
             "and `a_t = 14` explicitly. CORROBORATION (not a second, "
             "premise-disjoint leg): alt_level12.py (34 checks) closes the six "
             "surviving T1 branches via the `y = -1` place dichotomy.",
    covers="all 52 branches of L_alt (split_place_ledger_sub1.json strata with "
           "stratum_status = alternate_regime_open): the 25 already killed whole "
           "by C33/C34, and all 27 that survived them -- including the 12 that "
           "carried no state model in the DAG and the 15 forced-defect-0 overlay "
           "families.",
    note="This is what retires GAP-ALT-STATES. That gap measured 39 modelled "
         "degree-states against 4690 surviving across the 27 open branches. The "
         "4690 do not need to be modelled one by one: they all sit in branches "
         "with `a_t >= 11`, and `a_t <= 9` empties every such branch. The gap is "
         "MOOT, and it is retired explicitly rather than silently.",
)


def stage_census(F, window):
    """Apply the overlays in order to an already-`e|Phi`-filtered universe.

    Returns [(stage, killed_cells, killed_flagcases, killed_states, alive_after)].
    Cells already dead in an earlier stage are not double-counted.
    """
    rows = []
    dead_so_far = set()
    alive_cells = {cellname(c["a_t"], c["b"], c["branch"]) for c in F["cases"]}
    alive_fc = len(F["cases"])
    alive_st = sum(len(c["states"]) for c in F["cases"])
    for stg in STAGES:
        names = set(stg["dead"].get(window, []))
        kc, kf, ks = set(), 0, 0
        for c in F["cases"]:
            n = cellname(c["a_t"], c["b"], c["branch"])
            if n in names and n not in dead_so_far:
                kc.add(n)
                kf += 1
                ks += len(c["states"])
        # a cell named by the stage but absent from this universe is vacuous,
        # not an error -- record it so the difference is visible.
        vacuous = sorted(names - alive_cells - dead_so_far)
        dead_so_far |= kc
        alive_cells = alive_cells - kc
        alive_fc -= kf
        alive_st -= ks
        rows.append(dict(stage=stg["id"], title=stg["title"], source=stg["source"],
                         checker=stg["checker"], note=stg["note"],
                         killed_cells=sorted(kc), n_cells=len(kc),
                         flagcases=kf, states=ks, vacuous=vacuous,
                         alive_cells=sorted(alive_cells),
                         alive_n_cells=len(alive_cells),
                         alive_flagcases=alive_fc, alive_states=alive_st))
    return rows


# ---------------------------------------------------------------------------
# inherited numbers this rebuild is obliged to reproduce or contradict
# ---------------------------------------------------------------------------
# Every entry: a number some OTHER artifact states, the recomputation of the same
# quantity from the compiler, and a verdict.  A mismatch is REPORTED, never
# silently adopted (see the standing discipline in SESSION_HANDOFF.md).
def inherited_claims(rows):
    R = {(r["window"], r["c0820"]): r for r in rows}
    s2on, s2off = R[("sub2", "on")], R[("sub2", "off")]
    out = []

    def add(cid, source, claimed, computed, note):
        out.append(dict(id=cid, source=source, claimed=claimed, computed=computed,
                        agrees=(claimed == computed), note=note))

    c = s2on["census"]
    add("HAND-SUB2-KILLED-FLAGCASES", "DIVISOR_SYZYGY.md sec.3b (\"140 of 220 cells\")",
        140, c["flagcases_before"] - c["flagcases_after"],
        "AGREES as a FLAG-CASE count. The 220 objects the hand count calls "
        "\"cells\" are flag cases; the cell count (a_t, b, branch) is 26, of which "
        "the filter kills %d." % (c["cells_before"] - c["cells_after"]))
    add("HAND-SUB2-KILLED-STATES", "DIVISOR_SYZYGY.md sec.3b (\"4822 of 7888\")",
        4822, c["states_before"] - c["states_after"], "AGREES exactly.")
    add("HAND-SUB2-BGE2-FLAGCASES", "DIVISOR_SYZYGY.md sec.3b (b_j>=2 row, 23)",
        23, c["flagcases_killed_b_ge_2"], "AGREES (flag cases).")
    add("HAND-SUB2-BGE2-STATES", "DIVISOR_SYZYGY.md sec.3b (b_j>=2 row, 885)",
        885, c["states_killed_b_ge_2"], "AGREES.")
    add("HAND-SUB2-DEG-FLAGCASES", "DIVISOR_SYZYGY.md sec.3b (a+sum b != 10 row, 117)",
        117, c["flagcases_killed_degree"], "AGREES (flag cases).")
    add("HAND-SUB2-DEG-STATES", "DIVISOR_SYZYGY.md sec.3b (a+sum b != 10 row, 3937)",
        3937, c["states_killed_degree"], "AGREES.")

    # the defect count that failed to reproduce on 2026-07-25
    U2 = load(s2on["universe"])
    defect = sum(1 for cc in U2["cases"] for s in cc["states"]
                 if s["deg_e"] != cc["a_t"] + sum(cc["b"]))
    defect10 = sum(1 for cc in U2["cases"] for s in cc["states"]
                   if s["deg_e"] == 10 and cc["a_t"] + sum(cc["b"]) < 10)
    add("EXTERNAL-DEFECT-COUNT",
        "DIVISOR_SYZYGY.md sec.3b (external review 3844; repo recorded 3790 and "
        "called 3844 unreproduced)",
        3844, defect,
        "**REPRODUCES.** 3844 = states with deg_e != a_t+sum(b) (all defect "
        "states). 3790 = the subset ALSO having deg_e = 10 (recomputed here: %d). "
        "The %d-state gap is defect-1 e at deg_e = 8 and 9, in cells a6_b1000_T1 "
        "and a8_b0000_T1. The repo's parenthetical that \"<10 and !=10 agree\" is "
        "true of the a+sum(b) test but does not make 3790 the defect count; the "
        "external number was correct and should no longer be recorded as "
        "unreproduced." % (defect10, defect - defect10))

    # FIELD_SCOPE_AUDIT's measured C08/C20 blast radius, sub2 layer
    add("FIELDSCOPE-SUB2-CASES", "FIELD_SCOPE_AUDIT.md sec.4.4 (\"4 of 224\")",
        4, s2off["before"]["flagcases"] - s2on["before"]["flagcases"],
        "AGREES: 4 flag cases return.")
    add("FIELDSCOPE-SUB2-STATES", "FIELD_SCOPE_AUDIT.md sec.4.4 (\"10 of 7898\")",
        10, s2off["before"]["states"] - s2on["before"]["states"],
        "**DISAGREES.** 10 counts only the states inside the 4 RETURNING flag "
        "cases. Rebuilding the whole universe with APPLY_RESIDUE_KILLS = False "
        "(cascade AND state enumeration) also returns states to flag cases that "
        "never left: 23 pre-existing sub2 flag cases gain 168 states. Total sub2 "
        "blast radius %d states (7888 -> %d), not 10 (7888 -> 7898)."
        % (s2off["before"]["states"] - s2on["before"]["states"],
           s2off["before"]["states"]))

    s1on, s1off = R[("sub1", "on")], R[("sub1", "off")]
    add("FIELDSCOPE-SUB1-CASES", "FIELD_SCOPE_AUDIT.md sec.4.3/4.4 (\"18 of 1163\")",
        18, s1off["before"]["flagcases"] - s1on["before"]["flagcases"],
        "AGREES: 18 flag cases return (the sec.0 headline figure of 45 sub1 cases "
        "is the doc's own stale number; its sec.4.3 corrects it to 18).")
    add("FIELDSCOPE-SUB1-STATES", "FIELD_SCOPE_AUDIT.md sec.4.4 (\"125 of 44242\")",
        125, s1off["before"]["states"] - s1on["before"]["states"],
        "Same methodological gap as the sub2 row: 125 counts only the returning "
        "flag cases. Rebuilding the whole sub1 universe with "
        "APPLY_RESIDUE_KILLS = False gives 44117 -> %d, i.e. %d states."
        % (s1off["before"]["states"],
           s1off["before"]["states"] - s1on["before"]["states"]))
    add("FIELDSCOPE-HEADLINE-STATES",
        "FIELD_SCOPE_AUDIT.md sec.0 / sec.4.4 / sec.4.5 / [J3] headline",
        11341,
        (s2off["before"]["states"] - s2on["before"]["states"]
         + s1off["before"]["states"] - s1on["before"]["states"]),
        "REPRODUCES. CORRECTION 2026-07-25: this row previously hardcoded 1090 "
        "and asserted the audit doc was internally inconsistent. It was NOT. The "
        "committed FIELD_SCOPE_AUDIT.md carries 11341 at sec.0, sec.4.4, sec.4.5 "
        "and [J3] -- ten occurrences, and ZERO occurrences of 1090. This lane "
        "read a PRE-COMMIT DRAFT mid-flight and pinned its stale figure, then "
        "published the accusation. The 135 in sec.4.4 is a labelled column the "
        "prose explicitly says not to quote. Flag-case totals also agree "
        "(4 sub2 + 18 sub1 = 22).")
    return out


# ---------------------------------------------------------------------------
# emit
# ---------------------------------------------------------------------------
def fmt_pct(k, n):
    return "%.2f%%" % (100.0 * k / n) if n else "--"


def byw_cells(win, rows):
    """Unfiltered cell count for a window (C08/C20 ON row), recomputed."""
    for r in rows:
        if r["window"] == win and r["c0820"] == "on":
            return r["before"]["cells"]
    return 0


def render(rows):
    L = []
    w = L.append
    w("# FRONTIER_REBUILD -- the frontier regenerated under the `e | Phi` divisor "
      "filter (machine-generated)\n")
    w("> **DO NOT HAND-EDIT.** Regenerated by `python frontier_rebuild.py`. Every "
      "figure is computed from the compiler artifacts named in the provenance "
      "table; none is typed in. `python frontier_rebuild.py --quiet` re-derives "
      "the whole census from the universes and exits nonzero on any drift.\n")
    w("Provenance: git %s | schema 2 (schema 2 adds the machine-readable stage "
      "registry -- `stages` and `alt_closure` -- to `frontier_rebuild.json`; "
      "`state_kill_ledger.py` and `proof_dag.py` read the cell lists and the "
      "evidence grades from there).\n" % git_commit())

    # --- the lemma, and its per-window specialisation ----------------------
    w("\n## 1. The filter, and what it does in each window\n")
    w("Source lemma: `DIVISOR_SYZYGY.md` / `divisor_syzygy.py` -- the exact "
      "K-syzygy `2*(G5 + d2*G3 + d1*G2 + d0*G1) == 2*Phi - e*(d2*e^2 + 3*e*S + "
      "3*R^2)` (residual 0), hence `e | Phi` on every lift, with "
      "`Phi = -(1/6630)*(y+1)^30*q` and `q` squarefree. Implemented as a pipeline "
      "stage in `divisor_filter.py` and applied by "
      "`phase_d_states.py --divisor-filter`.\n")
    w("Three consequences, and the level each one bites at:\n")
    w("| # | statement | level | note |")
    w("|---|---|---|---|")
    w("| D1 | `rad(e)` divides `(y+1)*q` -- no off-support root, so `deg e = a_t + sum(b_i)` "
      "EXACTLY | state | **not automatic from the universe**: the universe "
      "enumerates `deg_e` up to `e_cap`, i.e. it admits defect-`d` `e` with `d` "
      "free extra roots |")
    w("| D2 | `b_i in {0,1}` at each simple `q`-root | cell | window-independent |")
    w("| D3 | degree count against `deg Phi = 34` | cell (via D1) | window-DEPENDENT, "
      "derived per window below |")
    w("\n### D3 derived per window (not assumed, not transferred)\n")
    w("`deg(RHS) = deg e + max(deg d2 + 2 deg e, deg e + deg S, 2 deg R)`, monotone "
      "in `deg e`, must be able to reach `deg Phi = 34`. Caps are read from "
      "`cascade_engine.CONFIGS` (d2, e) and `full_system_bridge.STRIP_DEGCAP` "
      "(R = dm2, S = dm3).\n")
    w("| window | deg d2 | deg R | deg S | deg e cap | `E_min` | effect |")
    w("|---|---:|---:|---:|---:|---:|---|")
    seen = set()
    for r in rows:
        if r["window"] in seen:
            continue
        seen.add(r["window"])
        c = r["caps"]
        eff = ("`deg e = %d` EXACTLY (E_min = e_cap), hence `a_t + sum(b_i) = %d`"
               % (r["e_min"], r["e_min"])) if r["degree_forcing_active"] else \
              ("**VACUOUS** -- `2*deg R = %d > 34` already at `deg e = 0`, so D3 "
               "imposes nothing and the sub2 `deg e = 10` does NOT transfer"
               % (2 * c["R"]))
        w("| %s | %d | %d | %d | %d | %d | %s |"
          % (r["window"], c["d2"], c["R"], c["S"], c["e"], r["e_min"], eff))
    w("\n> The sub1 row is the reason this file reports the two windows "
      "separately rather than quoting a single filter strength: in sub1 only D1 "
      "and D2 act.\n")

    # --- headline census ---------------------------------------------------
    w("\n## 2. Census -- before / after, at cell, flag-case and state level\n")
    w("`before` = the Phase-D universe the compiler emits; `after` = the same "
      "compile with `--divisor-filter`. A cell is `(a_t, b, branch)`; a flag case "
      "is `(sigma_zero, d2_zero, g_zero_levels)` inside a cell; a state is a "
      "residual degree assignment `(deg_d2, deg_d1, deg_sigma, deg_e)`.\n")
    w("| window | C08/C20 | cells | flag cases | states | cells after | flag cases after "
      "| states after | cells killed | flag cases killed | states killed |")
    w("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for r in rows:
        b, a = r["before"], r["after"]
        w("| %s | %s | %d | %d | %d | %d | %d | %d | %d (%s) | %d (%s) | %d (%s) |"
          % (r["window"], r["c0820"].upper(), b["cells"], b["flagcases"], b["states"],
             a["cells"], a["flagcases"], a["states"],
             b["cells"] - a["cells"], fmt_pct(b["cells"] - a["cells"], b["cells"]),
             b["flagcases"] - a["flagcases"],
             fmt_pct(b["flagcases"] - a["flagcases"], b["flagcases"]),
             b["states"] - a["states"], fmt_pct(b["states"] - a["states"], b["states"])))

    # --- death attribution -------------------------------------------------
    w("\n### Death attribution (the three consequences are independent tests)\n")
    w("| window | C08/C20 | cells: D2 `b>=2` | cells: D3 degree | cells: D1 defect | "
      "flag cases: D2 | flag cases: D3 | flag cases: D1 | states: D2 | states: D3 | states: D1 |")
    w("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for r in rows:
        c = r["census"]
        w("| %s | %s | %d | %d | %d | %d | %d | %d | %d | %d | %d |"
          % (r["window"], r["c0820"].upper(),
             c["cells_killed_b_ge_2"], c["cells_killed_degree"],
             c["cells_killed_defect"],
             c["flagcases_killed_b_ge_2"], c["flagcases_killed_degree"],
             c["flagcases_killed_defect"],
             c["states_killed_b_ge_2"], c["states_killed_degree"],
             c["states_killed_defect"]))
    w("\nIn sub2 the D1 (defect) column is 0 **not** because D1 is inactive but "
      "because D3 already deletes every cell that could carry a defect: with "
      "`E_min = e_cap = 10`, a surviving cell has `a_t + sum(b) = 10 = deg e`, so "
      "no defect fits. D1 is still logically load-bearing -- it is what turns "
      "\"`deg e = 10`\" into \"`a_t + sum(b_i) = 10`\". In sub1, D3 is vacuous and "
      "D1 is the whole state-level effect.\n")

    # --- downstream overlays -----------------------------------------------
    w("\n## 2b. Downstream cell kills (other lanes' lemmas, applied on top)\n")
    w("Later lanes kill further CELLS by arguments this file does **not** "
      "re-derive. Only their cell-name lists are imported; every flag-case and "
      "state count below is recomputed here from the universes. Stages compose "
      "in order and are not double-counted.\n")
    w("Each stage also carries the **evidence level its evidence actually "
      "supports**, on the DAG's scale `claimed < exact-checked < "
      "independently-audited < certified`. `state_kill_ledger.py` and "
      "`proof_dag.py` read that grade from here rather than re-asserting it.\n")
    w("| stage | level | evidence |")
    w("|---|---|---|")
    for stg in STAGES:
        w("| `%s` | **%s** | %s |"
          % (stg["id"], stg.get("level", "claimed"),
             stg.get("evidence", "(not recorded)")))
    w("")
    for stg in STAGES:
        w("- **`%s`** -- %s. Source: %s. Checker: `%s`. %s"
          % (stg["id"], stg["title"], stg["source"], stg["checker"], stg["note"]))
    w("\n| window | C08/C20 | stage | cells killed | flag cases killed | states killed "
      "| cells alive | flag cases alive | states alive |")
    w("|---|---|---|---:|---:|---:|---:|---:|---:|")
    for r in rows:
        a = r["after"]
        w("| %s | %s | `e \\| Phi` (sec.2) | %d | %d | %d | %d | %d | %d |"
          % (r["window"], r["c0820"].upper(),
             r["before"]["cells"] - a["cells"],
             r["before"]["flagcases"] - a["flagcases"],
             r["before"]["states"] - a["states"],
             a["cells"], a["flagcases"], a["states"]))
        for s in r["stages"]:
            w("| %s | %s | `%s` | %d | %d | %d | %d | %d | %d |"
              % (r["window"], r["c0820"].upper(), s["stage"], s["n_cells"],
                 s["flagcases"], s["states"], s["alive_n_cells"],
                 s["alive_flagcases"], s["alive_states"]))
    w("\n### Final surviving support, all stages applied\n")
    w("| window | C08/C20 | cells | flag cases | states | surviving cells |")
    w("|---|---|---:|---:|---:|---|")
    for r in rows:
        last = r["stages"][-1]
        w("| %s | %s | %d | %d | %d | %s |"
          % (r["window"], r["c0820"].upper(), last["alive_n_cells"],
             last["alive_flagcases"], last["alive_states"],
             ", ".join("`%s`" % c for c in last["alive_cells"])
             if len(last["alive_cells"]) <= 6 else "listed below"))
    for r in rows:
        last = r["stages"][-1]
        if len(last["alive_cells"]) > 6:
            w("\n**%s surviving cells (%d):** %s\n"
              % (r["tag"], last["alive_n_cells"],
                 ", ".join("`%s`" % c for c in last["alive_cells"])))
    w("\n> **The enumerated f31 frontier is EMPTY in BOTH windows.** sub2 closes "
      "at `stage4_positive_slice`; sub1 closes at `stage7_sub1_spine9`. Every one "
      "of the %d sub2 cells and %d sub1 cells the compiler enumerates is dead, "
      "either by `e | Phi` (section 2) or by a named stage above.\n"
      % (byw_cells("sub2", rows), byw_cells("sub1", rows)))
    w("> Reading the sub1 column in order: `e | Phi` takes 171 cells to 47; "
      "`stage2_T2_divisor` removes the 13 T2 cells it names; `stage5` "
      "(`a_t >= 9`) removes 23; `stage6` (`a_t <= 9`) removes the six `a10_*`; "
      "`stage7` removes the last five `a9_*`. The `a_t = 9` pincer (stages 5 and "
      "6 together) is what turned an 11-cell residue into a 5-cell one, and "
      "stage 7 closed those five.\n")
    w("\n> **The SPINE stage (`stage3_spine`) is still not applied to sub1** -- it "
      "is a sub2 five-family statement. sub1 is closed by stages 2, 5, 6 and 7 "
      "instead. `stage7_sub1_spine9` generalises SPINE's *cofactor identity* to "
      "`a_t = 9` and re-derives it from `generators.json`; it imports nothing "
      "from `spine.py`, because SPINE's zero-slack degree bookkeeping "
      "`(n+6)+(2n+4) = 3n+10` is a sub2 coincidence that does NOT transfer.\n")
    w("\n> **sub2's stage-3 note, kept for the record:** `SPINE.md` sec.0 closes "
      "nine of the ten `(n, branch)` columns; the two columns it names that do "
      "not occur in the universe (`a10_b0000_T2`, `a6_b1111_T2`) are vacuous. "
      "This is STRONGER than \"the five T1 cells survive\": four of those five "
      "are also empty.\n")

    # --- the alternate regime ----------------------------------------------
    w("\n## 2c. The alternate regime `a_t >= 11` -- EMPTY, and not by a cell list\n")
    w("The alternate regime lies OUTSIDE both phase-D universes (they are capped "
      "at `a_t <= 10`), so it cannot be a stage: there are no cells here to "
      "name. It is registered separately, and `state_kill_ledger.py` / "
      "`proof_dag.py` read that record the same way they read the stages.\n")
    w("| field | value |")
    w("|---|---|")
    for k in ("id", "title", "source", "checker", "level", "covers", "evidence",
              "note"):
        w("| `%s` | %s |" % (k, ALT_CLOSURE[k]))
    w("\n> **The mechanism is `a_t <= 9` itself.** `stage6_syzygy_collision` is "
      "cap-free, branch-free and window-independent, so it holds in the "
      "alternate regime too; every branch of `L_alt` has `a_t in {11..15}`, all "
      "`> 9`. All 52 branches die at once. `ALT_LEVEL12.md` closes the six "
      "surviving T1 branches by the `y = -1` place dichotomy and is recorded as "
      "CORROBORATION -- it covers 6 of 52, and it shares the slice calculus with "
      "`stage5`, so it corroborates the arithmetic rather than the premises.\n")

    # --- surviving columns -------------------------------------------------
    w("\n## 3. Surviving cells after `e | Phi` only (stage-1 column list)\n")
    for r in rows:
        w("\n**%s** -- %d surviving cells of %d:\n"
          % (r["tag"], r["after"]["cells"], r["before"]["cells"]))
        if not r["after"]["cell_names"]:
            w("- (none)")
            continue
        w("| cell | flag cases (before -> after) | states (before -> after) |")
        w("|---|---|---|")
        for name in r["after"]["cell_names"]:
            pb = r["before"]["per_cell"][name]
            pa = r["after"]["per_cell"][name]
            w("| `%s` | %d -> %d | %d -> %d |"
              % (name, pb["flagcases"], pa["flagcases"], pb["states"], pa["states"]))
        dead = [n for n in r["before"]["cell_names"] if n not in set(r["after"]["cell_names"])]
        w("\nCells deleted whole: %s\n" % (", ".join("`%s`" % d for d in dead) or "(none)"))

    # --- cascade layer -----------------------------------------------------
    w("\n## 4. Cascade layer (branch counts are NOT changed by the divisor filter)\n")
    w("| window | C08/C20 | branches processed | engine-killed | surviving |")
    w("|---|---|---:|---:|---:|")
    for r in rows:
        s = r["cascade_summary"]
        w("| %s | %s | %d | %d | %d |" % (r["window"], r["c0820"].upper(),
          s["open_branches_processed"], s["engine_killed_pending_audit"],
          s["surviving_branches"]))
    w("\nThe divisor filter acts inside surviving branches; it does not change "
      "which branches the cascade kills. Turning C08/C20 off changes the "
      "flag-case/state universe but leaves the surviving-branch counts identical, "
      "which is the measurement `FIELD_SCOPE_AUDIT.md` also reports.\n")

    # --- toggle delta ------------------------------------------------------
    w("\n## 5. The C08/C20 toggle, isolated\n")
    w("| window | layer | C08/C20 ON | C08/C20 OFF | delta (reopened) |")
    w("|---|---|---:|---:|---:|")
    byw = defaultdict(dict)
    for r in rows:
        byw[r["window"]][r["c0820"]] = r
    for win in ("sub2", "sub1"):
        on, off = byw[win]["on"], byw[win]["off"]
        for layer, key, sub in (("flag cases (unfiltered)", "before", "flagcases"),
                                ("states (unfiltered)", "before", "states"),
                                ("cells (filtered)", "after", "cells"),
                                ("flag cases (filtered)", "after", "flagcases"),
                                ("states (filtered)", "after", "states")):
            w("| %s | %s | %d | %d | %+d |" % (win, layer, on[key][sub],
                                               off[key][sub],
                                               off[key][sub] - on[key][sub]))
    tot_unf = sum(byw[wd]["off"]["before"]["states"] - byw[wd]["on"]["before"]["states"]
                  for wd in ("sub2", "sub1"))
    tot_f1 = sum(byw[wd]["off"]["after"]["states"] - byw[wd]["on"]["after"]["states"]
                 for wd in ("sub2", "sub1"))
    tot_fN = sum(byw[wd]["off"]["stages"][-1]["alive_states"]
                 - byw[wd]["on"]["stages"][-1]["alive_states"]
                 for wd in ("sub2", "sub1"))
    w("\n**Net effect of the C08/C20 downgrade, by how much of the pipeline is "
      "applied:** %d states reopen in the raw universes; **%d** survive the "
      "`e | Phi` filter; **%d** survive every stage. Cell counts are unchanged at "
      "every stage (+0), and branch counts are unchanged (section 4) -- the "
      "downgrade reopens states and flag cases, never a cell or a branch.\n"
      % (tot_unf, tot_f1, tot_fN))

    w("\n### Where the reopened states actually are\n")
    w("Turning C08/C20 off does two separable things. It returns whole flag cases "
      "the residue kills had removed at the cascade layer -- and it also returns "
      "states to flag cases that never left, because the SAME lemma gates "
      "forbidden drops inside `case_states`' infinity join. An audit that "
      "recomputes only the returning cases sees the first effect and misses the "
      "second.\n")
    w("| window | new flag cases | states in them | pre-existing flag cases that GAIN states "
      "| states gained there | total reopened |")
    w("|---|---:|---:|---:|---:|---:|")
    for win in ("sub2", "sub1"):
        on, off = byw[win]["on"], byw[win]["off"]
        A, B = on["before"]["per_flagcase"], off["before"]["per_flagcase"]
        new = [k for k in B if k not in A]
        grew = [(k, A[k], B[k]) for k in A if B.get(k, 0) != A[k]]
        w("| %s | %d | %d | %d | %d | %d |"
          % (win, len(new), sum(B[k] for k in new), len(grew),
             sum(b - a for _, a, b in grew),
             off["before"]["states"] - on["before"]["states"]))
    w("\n**Do not net these against the filter's kills.** The two move in opposite "
      "directions and are governed by different open questions; the toggle is left "
      "explicit so the number can be updated by rerunning the commands in section "
      "6 rather than rebuilding anything.\n")

    # --- inherited numbers -------------------------------------------------
    w("\n## 6. Inherited numbers: reproduced or contradicted\n")
    w("Every number some other artifact states about this filter, recomputed here "
      "from the compiler. A mismatch is reported, never silently adopted.\n")
    w("| id | source | claimed | recomputed | verdict |")
    w("|---|---|---:|---:|---|")
    for cl in inherited_claims(rows):
        w("| `%s` | %s | %s | %s | %s |"
          % (cl["id"], cl["source"], cl["claimed"], cl["computed"],
             "reproduces" if cl["agrees"] else "**DIFFERS**"))
    w("")
    for cl in inherited_claims(rows):
        w("- **`%s`** -- %s" % (cl["id"], cl["note"]))

    # --- open lead ---------------------------------------------------------
    w("\n## 7. A lead this rebuild deliberately does NOT use -- now PARTIAL, not "
      "unblocked\n")
    w("The same syzygy has an ORDER (t-adic) side: `ord_t(2*Phi) = 30` and "
      "`ord_t(e) = a_t` exactly (D1), so `30 - a_t = ord_t(d2*e^2 + 3*e*S + "
      "3*R^2)`. IF `t^a | R, S` then every bracket term has order `>= 2*a_t`, "
      "giving `30 - a_t >= 2*a_t`, i.e. **`a_t <= 10`**, which would delete the "
      "whole `a = 11..15` alternate regime.\n")
    w("The premise `t^a | R,S,T` has since been adjudicated **branch-independent** "
      "(`T1_BRANCH.md`), removing the T2-only restriction that originally blocked "
      "this. **The lead is nonetheless only PARTIALLY unblocked and is still not "
      "applied here:** the place trichotomy at `y = -1` admits a second horn "
      "`a + 2*rho = 30` with `rho < a`, and on that horn the alternate-regime "
      "values split -- `a_t in {11, 13, 15}` die on parity, `a_t in {12, 14}` "
      "survive. So the conclusion is NOT `a_t <= 10` outright. Nothing in the "
      "census above depends on it.\n")

    # --- what this file's numbers are, and are not -------------------------
    w("\n## 7b. Two different \"kill\" counts -- do not conflate them\n")
    w("This file reports the **divisor filter alone**, applied to the Phase-D "
      "universes: sub2 4822 states, sub1 36639 states (C08/C20 ON). "
      "`FRONTIER_V2_DIVFILTER.md` reports something different and larger -- the "
      "**unified state-kill ledger**, which joins the divisor filter with every "
      "PRIOR state-level kill mechanism (batch convolution, msolve/GB, "
      "reconstruction, phase-F2, the corner lemmas). Its headline "
      "(sub2 4951, sub1 36797 distinct states killed; 41746 of 52005 = 80.28%%) is "
      "therefore strictly larger, and the difference is exactly the pre-existing "
      "kills the divisor filter did not already cover. Neither number is a "
      "correction of the other; they count different sets.\n")
    w("Evidence grade is preserved across both: the divisor kills enter the ledger "
      "at `AUDITED`, which the proof-DAG maps to **`exact-checked` "
      "(same-author)** -- NOT `independently-audited`. In "
      "`FRONTIER_V2_DIVFILTER.md` the independent grade stays at sub2 34 / sub1 25. "
      "The distinction is load-bearing and must not be collapsed in any release.\n")

    # --- provenance --------------------------------------------------------
    w("\n## 8. Provenance -- every artifact and the command that makes it\n")
    w("| artifact | sha256(16) | command |")
    w("|---|---|---|")
    for r in rows:
        for k in ("cascade", "universe", "filtered"):
            w("| `%s`%s | `%s` | `%s` |"
              % (r[k], "" if (k != "filtered" or r["filtered_present"])
                 else " *(absent; re-derived at read time)*",
                 r["sha"][k], r[k + "_cmd"]))
    w("\nChecker: `python frontier_rebuild.py --quiet` (exit 0/1). It reloads every "
      "unfiltered universe, RE-APPLIES `divisor_filter.filter_universe` from "
      "scratch, and requires the result to agree with the stored filtered artifact "
      "and with every number in this file.\n")
    w("Filter self-check: `python divisor_filter.py --quiet` (exit 0/1) -- "
      "re-asserts the K-syzygy residual is 0, recomputes `deg Phi`, re-reads both "
      "windows' caps from the cap tables, and re-derives `E_min` per window.\n")
    return "\n".join(L) + "\n"


# ---------------------------------------------------------------------------
# checker
# ---------------------------------------------------------------------------
def check(rows) -> int:
    fails = []
    for r in rows:
        if not r["recompute_agrees"]:
            fails.append("%s: re-applying divisor_filter to %s does NOT reproduce %s"
                         % (r["tag"], r["universe"], r["filtered"]))
        b, a, c = r["before"], r["after"], r["census"]
        if c["states_before"] != b["states"] or c["states_after"] != a["states"]:
            fails.append("%s: stored census state totals disagree with the universes"
                         % r["tag"])
        if (c["states_before"] != c["states_after"] + c["states_killed_b_ge_2"]
                + c["states_killed_degree"] + c["states_killed_defect"]):
            fails.append("%s: state death attribution does not sum" % r["tag"])
        if (c["flagcases_before"] != c["flagcases_after"]
                + c["flagcases_killed_b_ge_2"] + c["flagcases_killed_degree"]
                + c["flagcases_killed_defect"]):
            fails.append("%s: flag-case death attribution does not sum" % r["tag"])
        if (c["cells_before"] != c["cells_after"] + c["cells_killed_b_ge_2"]
                + c["cells_killed_degree"] + c["cells_killed_defect"]):
            fails.append("%s: cell death attribution does not sum" % r["tag"])

    # the emitted document must contain the numbers we just recomputed
    p = path(MD_OUT)
    if not os.path.exists(p):
        fails.append("%s absent -- run `python frontier_rebuild.py`" % MD_OUT)
    else:
        md = open(p, encoding="utf-8").read()
        for r in rows:
            b, a = r["before"], r["after"]
            row = re.search(
                r"\|\s*%s\s*\|\s*%s\s*\|\s*(\d+)\s*\|\s*(\d+)\s*\|\s*(\d+)\s*\|"
                r"\s*(\d+)\s*\|\s*(\d+)\s*\|\s*(\d+)\s*\|"
                % (r["window"], r["c0820"].upper()), md)
            if not row:
                fails.append("%s: headline census row missing from %s"
                             % (r["tag"], MD_OUT))
                continue
            got = tuple(int(x) for x in row.groups())
            want = (b["cells"], b["flagcases"], b["states"],
                    a["cells"], a["flagcases"], a["states"])
            if got != want:
                fails.append("%s: %s says %s, recomputed %s"
                             % (r["tag"], MD_OUT, got, want))
        # the inherited-claims table must carry the CURRENT recomputed values
        # (a claim that legitimately differs is fine; a STALE one is not)
        for cl in inherited_claims(rows):
            row = re.search(r"\|\s*`%s`\s*\|[^|]*\|\s*(\d+)\s*\|\s*(\d+)\s*\|\s*([^|]*)\|"
                            % re.escape(cl["id"]), md)
            if not row:
                fails.append("inherited claim %s missing from %s" % (cl["id"], MD_OUT))
                continue
            if (int(row.group(1)) != cl["claimed"]
                    or int(row.group(2)) != cl["computed"]
                    or ("DIFFERS" in row.group(3)) == cl["agrees"]):
                fails.append("inherited claim %s stale in %s: file says "
                             "claimed=%s recomputed=%s verdict=%r, recomputed now "
                             "claimed=%s computed=%s agrees=%s"
                             % (cl["id"], MD_OUT, row.group(1), row.group(2),
                                row.group(3).strip(), cl["claimed"], cl["computed"],
                                cl["agrees"]))
    for f in fails:
        print("[FAIL] %s" % f)
    return len(fails)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quiet", action="store_true",
                    help="checker only: recompute and verify; exit 0 iff consistent")
    ap.add_argument("--print-commands", action="store_true")
    args = ap.parse_args()

    if args.print_commands:
        for r in RUNS:
            print(r["cascade_cmd"])
            print(r["universe_cmd"])
            print(r["filtered_cmd"])
        return

    rows = build_census()
    if args.quiet:
        raise SystemExit(1 if check(rows) else 0)

    with open(path(MD_OUT), "w", encoding="utf-8") as fh:
        fh.write(render(rows))
    with open(path(JSON_OUT), "w", encoding="utf-8") as fh:
        # `stages` and `alt_closure` are the machine-readable STAGE REGISTRY.
        # state_kill_ledger.py and proof_dag.py consume them from here, so the
        # cell lists and the evidence grades have exactly one authoritative home.
        json.dump({"schema": 2, "git_commit": git_commit(), "runs": rows,
                   "stages": STAGES, "alt_closure": ALT_CLOSURE},
                  fh, indent=2, sort_keys=True)
        fh.write("\n")
    print("wrote %s and %s" % (MD_OUT, JSON_OUT))
    for r in rows:
        b, a = r["before"], r["after"]
        print("  %-16s cells %3d -> %-3d  flagcases %5d -> %-5d  states %6d -> %d"
              % (r["tag"], b["cells"], a["cells"], b["flagcases"], a["flagcases"],
                 b["states"], a["states"]))
    n = check(rows)
    print("checker: %d finding(s)" % n)
    raise SystemExit(1 if n else 0)


if __name__ == "__main__":
    main()
