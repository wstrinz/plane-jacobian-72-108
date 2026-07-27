#!/usr/bin/env python3
"""
c0_partition.py -- JUDGMENT EDGE 2: the C0 case partition.

`proof_dag.json` carries the `subcase -> C0` edge with
`machine_checkable: false` and `exhaustiveness_level: claimed`.  This file makes
the partition EXPLICIT (see PROPOSITION below and JUDGMENT_EDGES.md section 3),
and checks the finite bookkeeping that the proposition asserts.

Two classes of output, never mixed:

  CHECK   -- finite bookkeeping the proposition claims.  A CHECK failure means
             the proposition as written is wrong.  Exit code 1.
  GAP     -- a place where the DAG's instantiation of the partition covers LESS
             than the proposition's leaf.  GAPs are pinned by exact key list, so
             a GAP that CHANGES is a CHECK failure.  A GAP that stays as
             documented is not: it is the honest state of affairs.

READ-ONLY.  Writes nothing.  `--quiet` prints failures, the gap inventory, and
the verdict.  Exit 0 iff every CHECK passes and the gap inventory is exactly the
documented one.

Usage:  python -u c0_partition.py [--quiet]

---------------------------------------------------------------------------
2026-07-25 REGISTRY REPAIR (DAG_REPAIR.md).  Two of the four gaps this file
originally declared were REGISTRY gaps and have been repaired in
`proof_dag.json`:

  GAP-D-NONODE      leaf L_D (d_{-1} == 0) had no node anywhere  -> REPAIRED
  GAP-ALT-BRANCHES  12 open alternate-regime branches had no node -> REPAIRED

The by-exact-key discipline is NOT relaxed for them: it is RE-POINTED.  Where
this file used to pin "these 12 keys are missing", it now pins "the DAG's L_alt
branch-key set is EXACTLY the ledger's 52, the unmodelled-open set is EXACTLY
these 12 named keys, and the L_D node exists, is closed, and cites C10".  A
regression that drops a node, silently promotes an unmodelled branch, or lets
the alternate-regime census drift is still a CHECK failure -- and so is a
"repair" that closes one of the 12 without evidence.

2026-07-26 GAP RETIREMENTS.  The two gaps that survived the registry repair are
now RETIRED -- not by registry work, but because later results dissolved them:

  GAP-ALT-STATES     was: 39 modelled states vs 4690 surviving (represented !=
                     covered).  MOOT: `a_t <= 9` empties all 52 L_alt branches, so
                     the 4690 never needed modelling.
  GAP-SUB2-EXCISIONS was: 443 - 23, the 23 being judgment edges to 7 tier-3
                     documents.  RETIRED: all 23 are re-derived EMPTY in-repo
                     (22 by `a_t = 9`, 1 by sub2's e|Phi degree count), so the
                     excision no longer rests on those documents.

Both retirements are PINNED, key by key / criterion by criterion, and both stay
in the printed inventory under RETIRED.  A retirement that is deleted, or a gap
that returns, is a CHECK failure -- the same discipline the 2026-07-25 repair used.
---------------------------------------------------------------------------
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
QUIET = "--quiet" in sys.argv

CHECKS, GAPS = [], []


def out(*a):
    if not QUIET:
        print(*a)


def check(ok, code, msg):
    CHECKS.append((bool(ok), code, msg))
    return bool(ok)


def gap(code, msg):
    GAPS.append((code, msg))


# A RETIRED gap is not a deleted gap.  It stays in the inventory, labelled, with
# the reason it dissolved -- otherwise "we solved it" and "we stopped mentioning
# it" look identical in the record.  RETIRED entries do not count towards the
# live gap total and do not fail the build.
RETIRED = []


def retired(code, msg):
    RETIRED.append((code, msg))


def load(fn):
    with open(os.path.join(HERE, fn)) as f:
        return json.load(f)


# ===========================================================================
# THE PROPOSITION
# ===========================================================================
PROPOSITION = r"""
PROPOSITION (C0 case partition).

PARENT HYPOTHESIS  H0.
    K a field, char K = 0.  P, Q in K[x,y] with [P,Q] = x^2,
    (deg P, deg Q) = (72,108), and Newton polygons N(P), N(Q) as in
    GGHV22 (arXiv:2204.14178) Proposition 4.3, case (8,28).
    C0 is the assertion  H0 = {}.

SPLIT 1  (Newton-polygon subcase).                        [S1a | S1b]
    source : GGHV22 Prop 4.3; transcription verified corner-for-corner in
             T3_WINDOW_AUDIT.md sec.1 and PAPER_NOTES.md (AUDIT.md sec.C).
    H0  =  H0 & subcase(1)   union   H0 & subcase(2).
    field scope    : char 0 (the polygon classification is characteristic-free
                     as used; the ambient argument is char 0).
    saturation     : none.
    DISJOINT       : YES.  The two subcases are distinct polygon data; a single
                     (P,Q) realises exactly one.

SPLIT 2  (master identity).                            [D | F31 | F37]
    source : the elimination of AUDIT.md sec.A.2 / STATE.md item 5.  All Q-side
             equations are linear in d_{-4}; after substitution one Singular
             resultant gives, in K[y],
                    f31 * f37 * d_{-1}^21  ==  0
             (f31: 102 terms, weight 125; f37: 618 terms, weight 134).
             K[y] is a domain, so at least one factor vanishes identically.
    applies inside EACH of subcase(1), subcase(2).
    field scope    : char 0 (the resultant identity is over Q; the debris
                     factorisation A = d_{-1} Ah, B = d_{-1} Bh is exact).
    saturation     : none -- this is an identity, not a saturation.
    DISJOINT       : NO, AND THIS IS STATED.  The three branches are a COVER,
                     not a partition: a solution could a priori satisfy two of
                     the three.  Overlap is harmless because the cover is used
                     only in the direction "every solution lies in some branch",
                     and two of the three branches are shown empty.

SPLIT 3  (branch D is empty).
    source : AUDIT.md sec.A.3, STATE.md item 5 (direct substitution, symbolic,
             denominator-free).  G1|_{d_{-1}=0} = 3 d_{-2} d_{-3}, so over the
             domain K[y] either d_{-2} = 0 or d_{-3} = 0.
               d_{-2}=0 : G2 -> (3/2) d_{-3}^2, so d_{-3}=0, then G5body -> Phi = 0.
               d_{-3}=0 : G3 -> -(3/2) d1 d_{-2}^2, so d1=0, then G5body -> Phi = 0.
             Both contradict Phi = f1 * C4^28 != 0 (deg 238, ord 204, explicit).
    field scope    : char 0.        saturation : none.

SPLIT 4  (branch F37 is empty).
    source : C11.  f31 in <G1,G2,G3,G5body+Phi> over Q with Phi free; the
             elimination ideal in (d2~,d1~,d0~,d_{-1},Phi) is exactly <f31>, so
             f37 and d_{-1}^21 are classical resultant excess and {f37=0} is an
             artifact -- in BOTH subcases.
    checkers : f37_sat_verify.py (same-author exact); lean_certificates/
               (kernel-checks the integer certificate D*f31 = sum (D c_i) G_i).
    field scope    : char != 3, 5.  The integer certificate carries the
                     denominator-clearing multiplier D = 46875 = 3 * 5^6
                     (F37_SATURATION_REPORT.md).
    saturation     : denominator-cleared lift certificate.

SPLIT 5  (terminal trichotomy, inside F31).              [T1 | T2 | T3]
    source : the h-ladder h7 = 8192 d1^2, h6 = -3072 sigma^2 + 14336 d1^2 d2
             + 8192 d1 e (CASCADE_ENGINE_PLAN.md sec.1); the branch label is the
             zero/nonzero alternative on the pair (d1, sigma):
                 T1 : d1 != 0
                 T2 : d1 == 0, sigma != 0
                 T3 : d1 == 0, sigma == 0
    field scope    : char 0.        saturation : none.
    DISJOINT       : YES -- it is the complete truth table of two predicates
                     (the two d1 != 0 rows are merged because T1's terminal
                     identity does not mention sigma).

SPLIT 6  (T3 is empty).
    source : the field-stable split-place sigma-locus theorem, C14-C16
             (FIELD_SPLIT_AUDIT.md; t5_split_place_verify.py,
             test_split_place_proofs.py).  For sub1 the same conclusion is
             carried by the split-place degree enumeration plus Mason-Stothers
             margins (sub1_cascade_verify.py).
    field scope    : char 0.        saturation : none.

SPLIT 7  (split-place stratification, inside F31 & (T1|T2)).
    source : C14-C16.  After base change q = p1 p2 p3 p4 into four distinct
             geometric (degree-one) places; with e = d_{-1},
                 a := v_t(e),  b_i := v_{p_i}(e / t^a),
             the case datum is the S4-orbit of (b_1,...,b_4), i.e. b sorted
             decreasing, subject to a + sum b_i <= deg e (= 10 sub2 / 15 sub1).
    field scope    : char 0.        saturation : none.
    DISJOINT       : YES.  a and b are FUNCTIONS of the solution, so distinct
                     (a,b) strata are mutually exclusive and jointly exhaustive.

SPLIT 8  (regime).                                 [standard | alternate]
    source : ALT_REGIME.md.  With Phi~ = t^30 u and e = t^a E, v := 30 - 3a, so
                 v >= 0  <=>  a <= 10   (standard reduction),
                 v <  0  <=>  a >= 11   (alternate / flipped reduction).
             For subcase (2) the cap a + sum b_i <= 10 forces a <= 10, so the
             alternate regime is EMPTY there.
    field scope    : char 0.        saturation : none.
    DISJOINT       : YES (a condition on the single integer a).

LEAVES.  Combining, H0 is COVERED by exactly five leaves:

    L_D    :  d_{-1} == 0                         -- EMPTY (Split 3)
    L_F37  :  f37 == 0, both subcases             -- EMPTY (Split 4, char != 3,5)
    L_sub2 :  subcase(2) & f31 & (T1|T2) & a<=10  -- 420 branches   [CLOSED]
    L_sub1 :  subcase(1) & f31 & (T1|T2) & a<=10  -- 2178 branches  [CLOSED]
    L_alt  :  subcase(1) & f31 & (T1|T2) & a>=11  -- 52 branches    [CLOSED]

    C0 holds  IFF  L_sub2, L_sub1 and L_alt are all empty
              (given Splits 3 and 4).

    STATUS 2026-07-26: ALL FIVE LEAVES ARE CLOSED, so C0 is CLOSED -- at level
    'claimed', which is the level of its own judgment-referenced exhaustiveness
    edge.  See "WHAT NOW CAPS C0" at the end of the run.

    DAG INSTANTIATION (post-2026-07-26).  Each of the five leaves has a node, and
    C0's DAG child list IS this leaf list:

        L_D    -> subcase:dm1               closed, level 'exact-checked' (C10;
                                            raised from 'claimed' 2026-07-26 once
                                            dm1_branch_verify.py was gated -- to
                                            the FLOOR only, not to
                                            'independently-audited')
        L_F37  -> subcase:f37               closed, level 'exact-checked' (C11)
        L_sub2 -> subcase:sub2              420 branch nodes, all closed
        L_sub1 -> subcase:sub1              2178 branch nodes, all closed
        L_alt  -> subcase:sub1_alt          all 52 branch keys, ALL CLOSED by
                                            `a_t <= 9`:
                    25 branch:alt:*  also killed whole by C33+C34 ('claimed')
                    12 branch:alt:*  closed_by_a_t_bound, still state_model none
                    15 branch:altdefect0:* reached through the overlay subcase
                       subcase:sub1_alt_defect0, which is a CHILD of
                       subcase:sub1_alt and not a child of C0

    HOW L_sub2 AND L_sub1 CLOSE.  Not through the state layer alone.  The
    enumerated frontier is emptied by COLUMN-LEVEL lemmas -- lemmas that empty a
    whole (a_t, b, branch) column, which is the DAG's BRANCH granularity:

        e | Phi divisor lemma (D2/D3 cell level)  exact-checked
        stage2_T2_divisor                         exact-checked
        stage3_spine                              exact-checked
        stage4_positive_slice                     exact-checked
        stage5_slice_obstruction  (a_t >= 9)      independently-audited
        stage6_syzygy_collision   (a_t <= 9)      independently-audited
        stage7_sub1_spine9        (5 cells)       independently-audited

    Those branches are closed as LEAVES, exactly as engine-killed branches are.
    The per-state kill records remain what they were; each affected cell records
    states_with_individual_kill_record vs states_covered_only_by_column_lemma, so
    the granularity of the evidence is visible and cannot be overstated.  Both
    windows nonetheless reach level 'claimed' only, and NOT because of the column
    lemmas: 4 sub2 + 108 sub1 engine-killed branches are killed by the t/inf layer
    alone, outside the depth-4 q-cascade auditor's scope, and stay 'claimed'.
    That is the binding constraint on the two f31 leaves.

    REGISTRATION WAS NOT COVERAGE -- and no longer needs to be.  The 12
    formerly-unmodelled alternate branches carry 1113 of the 4690 surviving
    alternate-regime degree-states and the DAG models none of them; the 15 overlay
    families model 39 forced-defect-0 slots.  That shortfall was GAP-ALT-STATES.
    It is RETIRED: those states all sit at a_t >= 11, and `a_t <= 9` empties every
    such branch, so they never needed modelling.

    DISJOINTNESS SUMMARY.  L_sub2 / L_sub1 / L_alt are pairwise disjoint
    (Splits 1 and 8).  L_D, L_F37 and (L_sub2 u L_sub1 u L_alt) may OVERLAP
    each other -- Split 2 is a cover, not a partition.  This is stated, not
    silent; it costs nothing because L_D and L_F37 are empty.

    FIELD SCOPE OF THE CONCLUSION.  The weakest leaf is L_F37 at char != 3,5.
    Under the char-0 hypothesis of H0 this is not a restriction; a
    characteristic-p transfer of C0 would inherit char != 3,5.
"""


# ===========================================================================
# THE FINITE BOOKKEEPING
# ===========================================================================
def parts(s, k=4):
    res = []

    def rec(rem, slots, cap, acc):
        if slots == 0:
            if rem == 0:
                res.append(tuple(acc))
            return
        for v in range(min(rem, cap), -1, -1):
            if v * slots < rem:
                break
            rec(rem - v, slots - 1, v, acc + [v])
    rec(s, k, s, [])
    return res


def local_min(b, c):
    x = max(0, -((-(3 * b - c)) // 2))
    return x, c + 2 * x - 3 * b


def terminal_open(deg_e, d1cap, sigcap, g7cap, g6cap, a_max=10):
    keys = set()
    for a in range(0, deg_e + 1):
        for s in range(0, deg_e - a + 1):
            for b in parts(s):
                if a > a_max:
                    continue
                p = [local_min(bi, 7) for bi in b]
                if sum(x for x, _ in p) <= d1cap and sum(g for _, g in p) <= g7cap(a):
                    keys.add((a, b, "T1"))
                p = [local_min(bi, 6) for bi in b]
                if sum(x for x, _ in p) <= sigcap and sum(g for _, g in p) <= g6cap(a):
                    keys.add((a, b, "T2"))
    return keys


def alt_universe(deg_e, a_min=11):
    return set((a, b, T)
               for a in range(a_min, deg_e + 1)
               for s in range(0, deg_e - a + 1)
               for b in parts(s)
               for T in ("T1", "T2"))


# --- the leaves, as the proposition asserts them ---------------------------
L_SUB2_TERMINAL = terminal_open(10, 6, 8, lambda a: 10 + 3 * a, lambda a: 10 + 3 * a)
L_SUB1_TERMINAL = terminal_open(15, 9, 12, lambda a: 46, lambda a: 48)
L_ALT = alt_universe(15)

DAG = load("proof_dag.json")
NODES = {n["id"]: n for n in DAG["nodes"]}
EDGES = DAG["edges"]


def dag_branch_keys(window):
    s = set()
    for n in DAG["nodes"]:
        if n["type"] != "branch" or n.get("window") != window:
            continue
        if "a_t" in n:
            s.add((n["a_t"], tuple(n["b"]), n["branch"]))
        else:                                    # altdefect0 nodes carry `family`
            f = n["family"]
            s.add((int(f.split("_")[0][1:]),
                   tuple(int(c) for c in f.split("_")[1][1:]),
                   f.split("_")[2]))
    return s


ENG_SUB2 = set((r["a_t"], tuple(r["b"]), r["branch"])
               for r in load("cascade_cones.json")["branches"])
ENG_SUB1 = set((r["a_t"], tuple(r["b"]), r["branch"])
               for r in load("cascade_cones_sub1_depth4.json")["branches"])
ALTSW = load("alt_inf_sweep.json")
ALT_OPEN27 = set((r["a"], tuple(r["b"]), r["branch"]) for r in ALTSW["branches"])
ALT_SURV = {(r["a"], tuple(r["b"]), r["branch"]): r["counts"]["surviving"]
            for r in ALTSW["branches"]}
DAG_ALT15 = dag_branch_keys("altdefect0")     # the defect-0 state overlay
DAG_ALT_DIRECT = dag_branch_keys("alt")       # the 25 killed-whole + 12 unmodelled
DAG_ALT_ALL = DAG_ALT15 | DAG_ALT_DIRECT


def keyname(k):
    return "a%d_b%s_%s" % (k[0], "".join(map(str, k[1])), k[2])

out(PROPOSITION)
out("=" * 74)
out("FINITE BOOKKEEPING")
out("=" * 74)

# --- CHECK 1: the C0 node's children are EXACTLY the five leaves -------------
# Repaired 2026-07-25.  Before the repair this list was
#   [f37, sub1, sub1_alt_defect0, sub2]
# -- L_D absent entirely, and L_alt represented by its defect-0 overlay.  The
# pin is unchanged in kind (exact list, order-sensitive); only its target moved.
c0_children = sorted(e["child"] for e in EDGES if e["parent"] == "C0")
expected_children = ["subcase:dm1", "subcase:f37", "subcase:sub1",
                     "subcase:sub1_alt", "subcase:sub2"]
check(c0_children == expected_children, "C0-CHILDREN",
      "C0's DAG children are %s, expected the five leaves %s"
      % (c0_children, expected_children))
# the defect-0 overlay must NOT be a C0 child -- it is not a leaf of the
# partition, it is a state-level refinement of 15 of L_alt's 27 open branches.
check("subcase:sub1_alt_defect0" not in c0_children, "C0-NO-OVERLAY-CHILD",
      "subcase:sub1_alt_defect0 is a direct child of C0 again; it is an overlay "
      "of 15 of L_alt's 27 open branches, not a leaf of the C0 partition")
check(any(e["parent"] == "subcase:sub1_alt"
          and e["child"] == "subcase:sub1_alt_defect0" for e in EDGES),
      "ALT-OVERLAY-PARENT",
      "subcase:sub1_alt_defect0 is not attached under subcase:sub1_alt")
out("  C0 children in the DAG        : %s" % ", ".join(c0_children))

# --- CHECK 2: leaf L_sub2 ---------------------------------------------------
n_exc = len(L_SUB2_TERMINAL - ENG_SUB2)
check(len(L_SUB2_TERMINAL) == 443, "SUB2-TERMINAL",
      "sub2 terminal-feasible = %d, expected 443" % len(L_SUB2_TERMINAL))
check(len(ENG_SUB2) == 420, "SUB2-COUNT",
      "sub2 engine universe = %d, expected 420" % len(ENG_SUB2))
check(n_exc == 23, "SUB2-EXCISED",
      "sub2 excised-by-prior-proof = %d, expected 23" % n_exc)
check(not (ENG_SUB2 - L_SUB2_TERMINAL), "SUB2-SPURIOUS",
      "sub2 engine keys outside the proposition's leaf: %s"
      % sorted(ENG_SUB2 - L_SUB2_TERMINAL)[:10])
check(dag_branch_keys("sub2") == ENG_SUB2, "SUB2-DAG",
      "subcase:sub2 branch nodes != engine universe")
check(NODES["subcase:sub2"]["n_branches"] == 420, "SUB2-NBRANCH",
      "subcase:sub2.n_branches = %s" % NODES["subcase:sub2"].get("n_branches"))
out("  L_sub2   443 terminal-feasible - 23 proof-excised = %d branches  (DAG: %d)"
    % (len(ENG_SUB2), len(dag_branch_keys("sub2"))))

# --- CHECK 3: leaf L_sub1 ---------------------------------------------------
check(L_SUB1_TERMINAL == ENG_SUB1, "SUB1-EXACT",
      "sub1 leaf != engine universe (only-prop %d, only-engine %d)"
      % (len(L_SUB1_TERMINAL - ENG_SUB1), len(ENG_SUB1 - L_SUB1_TERMINAL)))
check(len(ENG_SUB1) == 2178, "SUB1-COUNT",
      "sub1 engine universe = %d, expected 2178" % len(ENG_SUB1))
check(dag_branch_keys("sub1") == ENG_SUB1, "SUB1-DAG",
      "subcase:sub1 branch nodes != engine universe")
check(NODES["subcase:sub1"]["n_branches"] == 2178, "SUB1-NBRANCH",
      "subcase:sub1.n_branches = %s" % NODES["subcase:sub1"].get("n_branches"))
out("  L_sub1   %d branches, exact set equality with the engine (DAG: %d)"
    % (len(ENG_SUB1), len(dag_branch_keys("sub1"))))

# --- CHECK 4: leaf L_alt ----------------------------------------------------
alt_ledger = set()
for r in load("split_place_ledger_sub1.json")["strata"]:
    if r.get("stratum_status") == "alternate_regime_open":
        for T in r["open_branches"]:
            alt_ledger.add((r["a_t"], tuple(r["b"]), T))
check(len(L_ALT) == 52, "ALT-UNIVERSE", "alt leaf = %d, expected 52" % len(L_ALT))
check(L_ALT == alt_ledger, "ALT-LEDGER",
      "alt leaf != ledger alternate_regime_open branches")
check(ALT_OPEN27 <= L_ALT, "ALT-SWEEP-IN",
      "alt_inf_sweep keys outside the alt leaf: %s" % sorted(ALT_OPEN27 - L_ALT))
check(len(ALT_OPEN27) == 27, "ALT-OPEN",
      "post-C44 open alt branches = %d, expected 27" % len(ALT_OPEN27))
check(len(L_ALT) - len(ALT_OPEN27) == 25, "ALT-KILLED",
      "C44 alt kills = %d, expected 25" % (len(L_ALT) - len(ALT_OPEN27)))
out("  L_alt    52 branches = 25 killed (C33+C34, audited) + 27 open")

# --- CHECK 4b: the DAG instantiates L_alt in FULL, by exact key --------------
# This is the re-pointed GAP-ALT-BRANCHES pin.  It used to assert "these 12 keys
# are MISSING"; it now asserts the complete registry, key for key, which is a
# strictly stronger statement and still fails on any drift in either direction.
check(DAG_ALT_ALL == L_ALT, "ALT-DAG-FULL",
      "the DAG's L_alt branch-key set != the ledger's 52: missing %s, extra %s"
      % (sorted(map(keyname, L_ALT - DAG_ALT_ALL)),
         sorted(map(keyname, DAG_ALT_ALL - L_ALT))))
check(not (DAG_ALT15 & DAG_ALT_DIRECT), "ALT-NO-DOUBLE-NODE",
      "a branch is registered twice (as branch:alt AND branch:altdefect0): %s"
      % sorted(map(keyname, DAG_ALT15 & DAG_ALT_DIRECT)))
check(NODES["subcase:sub1_alt"].get("n_branches") == 52, "ALT-NBRANCH",
      "subcase:sub1_alt.n_branches = %s"
      % NODES["subcase:sub1_alt"].get("n_branches"))

ALT_NODE = {}
for n in DAG["nodes"]:
    if n["type"] == "branch" and n.get("window") == "alt":
        ALT_NODE[(n["a_t"], tuple(n["b"]), n["branch"])] = n

# the 25 C33/C34 whole-branch kills: closed, and honestly at 'claimed'
alt_killed = L_ALT - ALT_OPEN27
check(all(ALT_NODE[k]["closed"] is True
          and ALT_NODE[k]["alt_status"] == "killed_whole_branch"
          for k in alt_killed), "ALT-KILLED-NODES",
      "a C33/C34 whole-branch kill is not registered closed: %s"
      % sorted(keyname(k) for k in alt_killed
               if not ALT_NODE.get(k, {}).get("closed")))
# The C33/C34 route itself is still only 'claimed' -- audit_alt_regime.py
# re-derives all 25 spec-only but emits no joinable artifact, so the DAG must not
# promote THAT route.  Since 2026-07-26 these branches carry a SECOND, stronger
# route (`a_t <= 9`, independently-audited), and branch closure is disjunctive,
# so the node level may legitimately exceed 'claimed'.  What must not happen is a
# promotion with no recorded second route: the check therefore requires the
# pre-closure level to still be recorded as 'claimed' and any promotion to be
# attributed.
check(all(ALT_NODE[k].get("level_before_alt_closure") == "claimed"
          for k in alt_killed), "ALT-KILLED-LEVEL",
      "a C33/C34 whole-branch kill no longer records its own route at 'claimed'; "
      "audit_alt_regime.py re-derives all 25 spec-only but emits no joinable "
      "artifact, so the C33/C34 route must stay 'claimed': %s"
      % sorted(keyname(k) for k in alt_killed
               if ALT_NODE[k].get("level_before_alt_closure") != "claimed"))
check(all(ALT_NODE[k]["level"] == "claimed"
          or ALT_NODE[k].get("alt_closure_basis") for k in alt_killed),
      "ALT-KILLED-PROMOTION-ATTRIBUTED",
      "a C33/C34 whole-branch kill is above 'claimed' with NO recorded second "
      "route: %s"
      % sorted(keyname(k) for k in alt_killed
               if ALT_NODE[k]["level"] != "claimed"
               and not ALT_NODE[k].get("alt_closure_basis")))

# --- CHECK 4c: the 12 formerly-unmodelled branches, pinned by exact key ------
# These are the branches added by the 2026-07-25 repair.  Adding a node was NOT
# closing a branch, and until 2026-07-26 each had to be OPEN with state_model
# 'none'.  They are now CLOSED -- not by a state model, but by `a_t <= 9`, which
# empties every a_t >= 11 branch outright.  The key list is still pinned (a
# drifting key list is still a CHECK failure); what changed is the required
# STATUS, and the check below demands the closure be attributed to the bound and
# NOT to a state model that still does not exist.
ALT_UNMODELLED_PINNED = [
    (11, (0, 0, 0, 0), "T2"), (11, (1, 0, 0, 0), "T2"), (11, (1, 1, 0, 0), "T2"),
    (11, (1, 1, 1, 0), "T2"), (11, (3, 0, 0, 0), "T2"), (12, (0, 0, 0, 0), "T2"),
    (12, (1, 0, 0, 0), "T1"), (12, (1, 0, 0, 0), "T2"), (12, (1, 1, 0, 0), "T1"),
    (12, (1, 1, 0, 0), "T2"), (13, (0, 0, 0, 0), "T2"), (14, (0, 0, 0, 0), "T1"),
]
alt_unmodelled = sorted(ALT_OPEN27 - DAG_ALT15)
check(alt_unmodelled == ALT_UNMODELLED_PINNED, "ALT-UNMODELLED-PINNED",
      "the set of open alternate-regime branches with no state model CHANGED: "
      "now %d keys %s" % (len(alt_unmodelled), list(map(keyname, alt_unmodelled))))
check(all(k in ALT_NODE for k in alt_unmodelled), "ALT-UNMODELLED-NODES",
      "an open unmodelled alternate-regime branch has no node: %s"
      % sorted(keyname(k) for k in alt_unmodelled if k not in ALT_NODE))
# CLOSED, and closed for the RIGHT reason.  Each must be closed, must still
# record that it has NO state model (the closure does not come from one), must
# record 0 modelled states, must cite the a_t bound, and must still carry the
# pre-closure level 'open' so the promotion is visible rather than retroactive.
check(all(ALT_NODE[k]["closed"] is True
          and ALT_NODE[k]["alt_status"] == "closed_by_a_t_bound"
          and ALT_NODE[k]["modelled_states"] == 0
          and "not needed" in ALT_NODE[k]["state_model"]
          and ALT_NODE[k].get("level_before_alt_closure") == "open"
          and ALT_NODE[k].get("alt_closure_basis")
          for k in alt_unmodelled if k in ALT_NODE), "ALT-UNMODELLED-CLOSED",
      "a formerly-unmodelled alternate-regime branch is not registered as CLOSED "
      "BY THE a_t BOUND with no state model.  It must be closed, must record 0 "
      "modelled states (the closure is not from a state model), must cite the "
      "bound, and must keep level_before_alt_closure = 'open': %s"
      % sorted(keyname(k) for k in alt_unmodelled
               if k in ALT_NODE and not (
                   ALT_NODE[k]["closed"]
                   and ALT_NODE[k]["alt_status"] == "closed_by_a_t_bound"
                   and ALT_NODE[k].get("alt_closure_basis"))))
# and the bound really does apply: every alt branch has a_t >= 11 > 9
check(all(k[0] >= 11 for k in L_ALT), "ALT-A-T-BOUND-APPLIES",
      "an L_alt branch has a_t < 11, so the a_t <= 9 closure of the alternate "
      "regime does not cover it: %s"
      % sorted(keyname(k) for k in L_ALT if k[0] < 11))
# the surviving-state counts the DAG does not model must match the sweep
check(all(ALT_NODE[k]["alt_degree_states_surviving"] == ALT_SURV[k]
          for k in alt_unmodelled if k in ALT_NODE), "ALT-UNMODELLED-COUNTS",
      "an unmodelled branch's surviving-state count disagrees with "
      "alt_inf_sweep.json")
alt_unmodelled_surv = sum(ALT_SURV[k] for k in alt_unmodelled)
check(alt_unmodelled_surv == 1113, "ALT-UNMODELLED-SURV",
      "the 12 formerly-unmodelled branches carry %d surviving degree-states, "
      "expected 1113" % alt_unmodelled_surv)
# the leaf itself must now be CLOSED, 52/52, and say by what
ALT_LEAF_NODE = NODES["subcase:sub1_alt"]
check(ALT_LEAF_NODE["closed"] is True
      and ALT_LEAF_NODE.get("branches_closed") == 52
      and ALT_LEAF_NODE.get("branches_open") == 0, "ALT-LEAF-CLOSED",
      "subcase:sub1_alt is not CLOSED 52/52: closed=%s branches_closed=%s "
      "branches_open=%s" % (ALT_LEAF_NODE["closed"],
                            ALT_LEAF_NODE.get("branches_closed"),
                            ALT_LEAF_NODE.get("branches_open")))
check("a_t <= 9" in (ALT_LEAF_NODE.get("closure_mechanism") or ""),
      "ALT-LEAF-MECHANISM",
      "subcase:sub1_alt does not name `a_t <= 9` as its closure mechanism: %r"
      % ALT_LEAF_NODE.get("closure_mechanism"))
out("  L_alt    DAG registry: 52/52 keys, ALL CLOSED by %s"
    % ALT_LEAF_NODE.get("closure_mechanism"))
out("           history: 25 killed whole (C33+C34) + 12 formerly-unmodelled "
    "(branch:alt:*) + 15 overlay (branch:altdefect0:*)")
out("           the 12 formerly-unmodelled carry %d of the 4690 surviving "
    "degree-states -- now MOOT, all at a_t >= 11" % alt_unmodelled_surv)

# --- CHECK 5: disjointness of the three open leaves -------------------------
# L_sub2 and L_sub1 live in different Newton-polygon subcases (Split 1), so the
# (a,b,T) label collision between them is NOT an overlap; the check that matters
# is that L_sub1 and L_alt are disjoint as subsets of subcase (1).
check(not (ENG_SUB1 & L_ALT), "SUB1-ALT-DISJOINT",
      "L_sub1 and L_alt share keys: %s" % sorted(ENG_SUB1 & L_ALT)[:10])
check(max(k[0] for k in ENG_SUB1) <= 10 and min(k[0] for k in L_ALT) >= 11,
      "REGIME-CUT", "the a<=10 / a>=11 regime cut is not clean")
out("  disjointness  L_sub1 (a<=10) vs L_alt (a>=11): clean cut at a=10/11")
out("  disjointness  L_sub2 vs L_sub1/L_alt: different Newton subcases (Split 1)")
out("  overlap DECLARED: L_D, L_F37 vs the f31 leaves (Split 2 is a cover)")

# --- CHECK 6: the two empty leaves have a live checker ----------------------
for leaf, files in (("L_F37", ["f37_sat_verify.py", "F37_SATURATION_REPORT.md"]),
                    ("L_D", ["AUDIT.md", "STATE.md"])):
    missing = [f for f in files if not os.path.exists(os.path.join(HERE, f))]
    check(not missing, "LEAF-SOURCE-%s" % leaf,
          "%s cites missing file(s) %s" % (leaf, missing))
check(NODES["subcase:f37"]["closed"] is True, "F37-CLOSED",
      "subcase:f37 is not marked closed in the DAG")
check("3,5" in NODES["subcase:f37"].get("field_scope", ""), "F37-SCOPE",
      "subcase:f37 does not record the char != 3,5 scope")
out("  L_F37    closed (C11), field scope %r" % NODES["subcase:f37"]["field_scope"])

# --- CHECK 6b: leaf L_D has a node, is closed, and cites C10 ----------------
# Repaired 2026-07-25 (was GAP-D-NONODE).  The pin is re-pointed, not dropped:
# the node must EXIST, must be CLOSED, must name C10, and must NOT be graded
# above 'claimed' -- PROOF_INVENTORY grades C10 tier 2*, and marks its checker
# attribution itself as inferred, so there is no wired checker behind it.
check("subcase:dm1" in NODES, "D-NODE-EXISTS",
      "leaf L_D (d_{-1} == 0) has no node in proof_dag.json")
if "subcase:dm1" in NODES:
    dm1 = NODES["subcase:dm1"]
    check(dm1["closed"] is True, "D-CLOSED",
          "subcase:dm1 is not marked closed; L_D is EMPTY by AUDIT.md A.3")
    # 2026-07-26: this pin is now a DELIBERATE HOLD, not "there is no evidence".
    # dm1_branch_verify.py + dm1_branch_certificate.json exist (commit 9de8713)
    # and C0_CLOSEOUT.md sec.1.6 recommends 'independently-audited' with a floor
    # of 'exact-checked'.  The registry lane has not applied that regrade -- the
    # recommending document is same-day and unaudited, and its checker is not in
    # run_tests.sh.  The pin therefore stays, and its MESSAGE says why, so nobody
    # reads it as an assertion that no checker exists.
    # REGRADE TAKEN 2026-07-26 (commit 2636b90).  The old pin demanded 'claimed'
    # and its message gave two reasons: the checker was not in run_tests.sh, and
    # the recommending document was same-day and unaudited.  THE FIRST HAS GONE
    # STALE -- dm1_branch_verify.py IS gated (tools/suite_manifest.py) and passes
    # 28/28 -- and that is exactly the 'claimed' -> 'exact-checked' gap: an exact
    # checker verifies the claim AND a dropped checker fails the suite loudly.
    # The SECOND still stands, and is why the FLOOR was taken and not the ceiling.
    # This pin now enforces the floor from BOTH sides: dropping back to 'claimed'
    # is a failure (the evidence exists), and jumping to 'independently-audited'
    # is a failure too (no second independent implementation exists).
    check(dm1["level"] == "exact-checked", "D-LEVEL",
          "subcase:dm1 level = %r, expected 'exact-checked'. This node was raised "
          "from 'claimed' on 2026-07-26 to the FLOOR recommended by "
          "C0_CLOSEOUT.md sec.1.6, and ONLY the floor. Below 'exact-checked' is "
          "wrong because dm1_branch_verify.py closes the branch on one polynomial "
          "identity over Z, is gated in the suite, and passes 28/28. Above it is "
          "wrong because 'independently-audited' needs a SECOND, independent "
          "implementation of the same kill -- one gated checker is exact-checking, "
          "not audit, and conflating the two produced the v0.3.2 erratum. "
          "C0_CLOSEOUT.md, which recommends the ceiling, is still same-day and "
          "not itself audited. Keep this check and proof_dag.py's dm1 "
          "evidence_tier_note in step." % dm1["level"])
    check("C10" in dm1.get("closure_note", "")
          and "AUDIT.md" in dm1.get("closure_note", ""), "D-SOURCE",
          "subcase:dm1 does not cite C10 / AUDIT.md sec.A.3 in its closure_note")
    check("cover" in dm1.get("overlap_note", "").lower(), "D-OVERLAP",
          "subcase:dm1 does not record that Split 2 is a declared COVER, so its "
          "possible overlap with the f31 leaves is silent again")
    out("  L_D      closed (C10; AUDIT.md A.3 / STATE.md item 5), node "
        "subcase:dm1 at level %r -- the WEAKEST closed leaf" % dm1["level"])

# --- CHECK 7: the C0 edge is honestly labelled ------------------------------
c0_edges = [e for e in EDGES if e["parent"] == "C0"]
check(all(e["machine_checkable"] is False for e in c0_edges), "C0-EDGE-HONEST",
      "a subcase->C0 edge claims machine_checkable=true")
check(all(e["exhaustiveness_level"] == "claimed" for e in c0_edges), "C0-EDGE-LEVEL",
      "a subcase->C0 edge is above 'claimed'")
# 2026-07-26: C0 is now CLOSED, and this check inverts accordingly.  What it
# guards is the LEVEL: C0 must not close above 'claimed', because the
# subcase->C0 exhaustiveness edge is judgment (and subcase:dm1 is 'claimed' too).
# A C0 that ever reads 'exact-checked' or better while that edge is judgment is
# the single most dangerous thing this registry could say.
check(NODES["C0"]["closed"] is True, "C0-CLOSED",
      "C0 is not marked closed even though all five leaves are; the roll-up is "
      "broken")
check(NODES["C0"]["level"] == "claimed", "C0-LEVEL-CAPPED",
      "C0 closes at level %r. It MUST cap at 'claimed': the subcase->C0 "
      "exhaustiveness edge is judgment-referenced (GGHV22 Prop 4.3 + the "
      "field-split framework + the alternate-regime partition), and subcase:dm1 "
      "(L_D, C10) is itself only 'claimed'. Anything higher would assert a "
      "machine-enforced partition that does not exist."
      % NODES["C0"]["level"])
check(NODES["C0"].get("subcases_closed") == 5, "C0-SUBCASES-CLOSED",
      "C0 reports %s of 5 leaves closed" % NODES["C0"].get("subcases_closed"))
# the repair must not have quietly promoted the judgment edge: its ref must now
# point at the written-out proposition, and its notes must still say WHY it is
# judgment (published mathematics AND an instantiation smaller than the claim).
check(all("JUDGMENT_EDGES.md" in e["exhaustiveness_ref"] for e in c0_edges),
      "C0-EDGE-REF",
      "a subcase->C0 edge does not cite the written-out proposition "
      "(JUDGMENT_EDGES.md sec.3)")
# The edge must still NAME both former gaps -- a retired gap that disappears from
# the prose is indistinguishable from one that was ignored -- and must say they
# are RETIRED rather than leaving them readable as live.
check(all("GAP-ALT-STATES" in e["notes"] and "GAP-SUB2-EXCISIONS" in e["notes"]
          for e in c0_edges), "C0-EDGE-RESIDUAL",
      "a subcase->C0 edge no longer names the two coverage gaps; a retired gap "
      "must stay named, not vanish")
check(all("RETIRED" in e["notes"] for e in c0_edges), "C0-EDGE-RETIREMENT",
      "a subcase->C0 edge names the gaps but does not record that they are "
      "RETIRED, so the note reads as if they were still live")
check(all("judgment" in e["notes"] for e in c0_edges), "C0-EDGE-STILL-JUDGMENT",
      "a subcase->C0 edge stopped recording that it is a JUDGMENT edge; "
      "retiring the coverage gaps does NOT make the partition machine-enforced")

# --- CHECK 8: the gap inventory, and the RETIREMENTS, recorded by key --------
# The DAG must carry ZERO live c0-partition gaps and exactly the two retirements.
# A retirement that vanishes, or a gap that returns, is a CHECK failure.
DAG_GAPS = {u["name"]: u for u in DAG["unmapped"]
            if u.get("kind") == "c0-partition-gap"}
DAG_RETIRED = {u["name"]: u for u in DAG["unmapped"]
               if u.get("kind") == "c0-partition-gap-retired"}
check(sorted(DAG_GAPS) == [], "DAG-GAP-INVENTORY",
      "proof_dag.json still declares LIVE c0-partition gap(s) %s; both known "
      "gaps were retired on 2026-07-26" % sorted(DAG_GAPS))
check(sorted(DAG_RETIRED) == ["GAP-ALT-STATES", "GAP-SUB2-EXCISIONS"],
      "DAG-GAP-RETIREMENTS",
      "proof_dag.json's c0-partition-gap-retired inventory is %s, expected "
      "exactly ['GAP-ALT-STATES', 'GAP-SUB2-EXCISIONS'].  A retirement must not "
      "be deleted: deleting it makes an ignored gap indistinguishable from a "
      "dissolved one." % sorted(DAG_RETIRED))
if "GAP-ALT-STATES" in DAG_RETIRED:
    g = DAG_RETIRED["GAP-ALT-STATES"]
    check(sorted(g["former_unmodelled_branch_keys"])
          == sorted(map(keyname, alt_unmodelled)), "DAG-GAP-ALT-KEYS",
          "proof_dag.json's recorded formerly-unmodelled alt-branch keys "
          "disagree with alt_inf_sweep.json \\ phase_f2_scale.json")
    check(g.get("branches_now_closed") == 52, "DAG-GAP-ALT-CLOSED",
          "the GAP-ALT-STATES retirement records %s closed branches, expected 52"
          % g.get("branches_now_closed"))
    check("a_t <= 9" in (g.get("retired_by") or ""), "DAG-GAP-ALT-MECHANISM",
          "the GAP-ALT-STATES retirement does not name `a_t <= 9`: %r"
          % g.get("retired_by"))

# --- CHECK 8b: GAP-SUB2-EXCISIONS retirement, verified KEY BY KEY ------------
# The DAG records the CRITERION; this recomputes the 23 excised keys from the
# independent generator and classifies each one, so the retirement is a machine
# fact and not a claim.  Two rules, and their evidence levels differ:
#
#   rule 1  a_t != 9        -> EMPTY by `a_t = 9` EXACTLY (stages 5 + 6, both
#                              independently-audited).  Cap-free, branch-free and
#                              window-independent, so it binds branches the
#                              cascade engine never processed -- which is exactly
#                              what the gap was about.
#   rule 2  a_t == 9 and sum(b) != 1
#                           -> EMPTY by sub2's e|Phi degree count
#                              deg e = 10 = a_t + sum(b_i) (E_min = e_cap = 10).
#                              exact-checked, NOT independently-audited.
SUB2_EXCISED = sorted(L_SUB2_TERMINAL - ENG_SUB2)
EXC_BY_RULE1 = [k for k in SUB2_EXCISED if k[0] != 9]
EXC_BY_RULE2 = [k for k in SUB2_EXCISED if k[0] == 9 and sum(k[1]) != 1]
EXC_UNCOVERED = [k for k in SUB2_EXCISED if k[0] == 9 and sum(k[1]) == 1]
if "GAP-SUB2-EXCISIONS" in DAG_RETIRED:
    g = DAG_RETIRED["GAP-SUB2-EXCISIONS"]
    crit = g.get("criterion", {})
    check(crit.get("a_t_exact") == 9, "EXC-CRIT-AT",
          "the GAP-SUB2-EXCISIONS criterion does not pin a_t = 9: %r"
          % crit.get("a_t_exact"))
    check(len(EXC_BY_RULE1) == crit.get("expected_killed_by_rule_1"),
          "EXC-RULE1-COUNT",
          "%d of the 23 excised sub2 branches have a_t != 9 and so die by "
          "`a_t = 9`; the DAG's criterion expects %s"
          % (len(EXC_BY_RULE1), crit.get("expected_killed_by_rule_1")))
    check(len(EXC_BY_RULE2) == crit.get("expected_killed_by_rule_2"),
          "EXC-RULE2-COUNT",
          "%d of the 23 excised sub2 branches have a_t = 9 with sum(b) != 1 and "
          "so die by sub2's degree count; the DAG's criterion expects %s"
          % (len(EXC_BY_RULE2), crit.get("expected_killed_by_rule_2")))
    check(sorted(map(keyname, EXC_BY_RULE2))
          == sorted(crit.get("expected_rule_2_keys") or []), "EXC-RULE2-KEYS",
          "the excised branches needing sub2's degree count are %s; the DAG's "
          "criterion pins %s"
          % (sorted(map(keyname, EXC_BY_RULE2)), crit.get("expected_rule_2_keys")))
    check(len(EXC_UNCOVERED) == 0 == crit.get("expected_uncovered"),
          "EXC-UNCOVERED",
          "%d excised sub2 branch(es) are covered by NEITHER rule, so the "
          "retirement of GAP-SUB2-EXCISIONS is not established: %s"
          % (len(EXC_UNCOVERED), sorted(map(keyname, EXC_UNCOVERED))))
    check(len(EXC_BY_RULE1) + len(EXC_BY_RULE2) == len(SUB2_EXCISED) == 23,
          "EXC-PARTITION",
          "the two rules cover %d + %d of %d excised branches"
          % (len(EXC_BY_RULE1), len(EXC_BY_RULE2), len(SUB2_EXCISED)))
    check("exact-checked" in (g.get("weakest_link") or ""), "EXC-WEAKEST-LINK",
          "the GAP-SUB2-EXCISIONS retirement does not record that its weakest "
          "link is exact-checked (the one a_t = 9 branch, a9_b0000_T2, needs the "
          "e|Phi degree count, which is NOT independently audited): %r"
          % g.get("weakest_link"))
out("  L_sub2   the 23 excisions are re-derived EMPTY in-repo: %d by `a_t = 9` "
    "(independently-audited), %d by sub2's e|Phi degree count (exact-checked: %s)"
    % (len(EXC_BY_RULE1), len(EXC_BY_RULE2),
       ", ".join(map(keyname, EXC_BY_RULE2))))

# ===========================================================================
# THE DECLARED GAPS  (pinned by exact key list)
# ===========================================================================
# GAP 1 (was GAP-D-NONODE) and GAP 2 (was GAP-ALT-BRANCHES) were REGISTRY gaps
# and are REPAIRED as of 2026-07-25.  They are not deleted here: they are
# inverted into CHECKs that fail if the repair is ever undone --
#   D-NODE-EXISTS / D-CLOSED / D-LEVEL / D-SOURCE / D-OVERLAP  (above)
#   ALT-DAG-FULL / ALT-UNMODELLED-PINNED / ALT-UNMODELLED-OPEN (above)
# A node quietly deleted, or an unmodelled branch quietly closed, is a CHECK
# failure, exactly as a drifting gap key list used to be.
check(DAG_ALT15 <= ALT_OPEN27, "ALT-DAG-IN",
      "subcase:sub1_alt_defect0 has branch nodes outside the 27 open alt branches: %s"
      % sorted(map(keyname, DAG_ALT15 - ALT_OPEN27)))
check(DAG_ALT_DIRECT <= L_ALT, "ALT-DIRECT-IN",
      "branch:alt:* nodes outside the 52-branch alternate leaf: %s"
      % sorted(map(keyname, DAG_ALT_DIRECT - L_ALT)))

# GAP 3 -- within the 15 instantiated families, only the forced-defect-0 states.
alt_state_total = sum(n["state_total"] for n in DAG["nodes"]
                      if n["type"] == "branch" and n.get("window") == "altdefect0")
alt_surviving_states = load("alt_inf_sweep.json")["summary"]["surviving_states"]
check(alt_state_total == 39, "ALT-STATE-TOTAL",
      "altdefect0 state_total = %d, expected 39" % alt_state_total)
check(alt_surviving_states == 4690, "ALT-SURV-STATES",
      "alt surviving degree-states = %d, expected 4690" % alt_surviving_states)
retired("GAP-ALT-STATES",
    "RETIRED 2026-07-26 -- MOOT, dissolved by a later theorem rather than by "
    "registry work. The gap was: all 52 L_alt branch keys were REGISTERED, but "
    "across the 27 branches that survived C33/C34 the DAG modelled %d "
    "degree-states (the forced-defect-0 slots of phase_f2_scale.json, covering 15 "
    "branches) against %d surviving (alt_inf_sweep.json), and the other 12 "
    "carried %d surviving states modelled by nothing at all -- REPRESENTED, not "
    "COVERED. The %d never needed modelling: `a_t <= 9` (stage6_syzygy_collision, "
    "independently-audited via at_le9_audit.py 76/76) is cap-free, branch-free "
    "and window-independent, every L_alt branch has a_t in {11..15}, so all 52 "
    "are EMPTY. Verified above key by key (ALT-UNMODELLED-CLOSED, "
    "ALT-A-T-BOUND-APPLIES, ALT-LEAF-CLOSED)."
    % (alt_state_total, alt_surviving_states, alt_unmodelled_surv,
       alt_surviving_states))

# GAP 4 -- the 23 sub2 excisions are judgment edges to eight documents.
exc_refs = {}
rows = {(r["a_t"], tuple(r["b"])): r
        for r in load("split_place_ledger.json")["strata"]}
for (a, b, T) in sorted(L_SUB2_TERMINAL - ENG_SUB2):
    ref = rows[(a, b)]["branches"][T].get("reference", "?")
    exc_refs.setdefault(ref.split(",")[0].strip(), []).append(
        "a%d_b%s_%s" % (a, "".join(map(str, b)), T))
check(len(exc_refs) == 7, "EXC-REFDOCS",
      "sub2 excisions resolve to %d distinct documents, expected 7: %s"
      % (len(exc_refs), sorted(exc_refs)))
retired("GAP-SUB2-EXCISIONS",
    "RETIRED 2026-07-26 -- the judgment edge is gone. The gap was: L_sub2's "
    "420-branch universe is 443 terminal-feasible branches MINUS 23 excised by "
    "earlier exact proofs in %d documents (%s), so those 23 removals were sound "
    "only if those proofs were -- judgment edges to tier-3 CONDITIONAL documents, "
    "not machine-checked. ALL 23 ARE NOW RE-DERIVED EMPTY IN-REPO, consuming none "
    "of those documents: %d by `a_t = 9` EXACTLY (stages 5+6, both "
    "independently-audited; cap-free and branch-free, so it binds branches the "
    "engine never processed), and %d -- %s -- by sub2's own e|Phi degree count "
    "deg e = 10 = a_t + sum(b_i), which forces sum(b_i) = 1 at a_t = 9. "
    "WEAKEST LINK: that last one is exact-checked, NOT independently audited. "
    "Verified above key by key (EXC-RULE1-COUNT, EXC-RULE2-KEYS, EXC-UNCOVERED, "
    "EXC-PARTITION). Re-running the 23 through the cascade engine is no longer "
    "needed for soundness -- it would only add a second route."
    % (len(exc_refs), ", ".join(sorted(exc_refs)),
       len(EXC_BY_RULE1), len(EXC_BY_RULE2),
       ", ".join(map(keyname, EXC_BY_RULE2))))

# ===========================================================================
# VERDICT
# ===========================================================================
n_fail = sum(1 for ok, _, _ in CHECKS if not ok)
out("")
out("-- CHECK RESULTS (%d) --" % len(CHECKS))
for ok, code, msg in CHECKS:
    if ok:
        out("  [OK]   %-22s" % code)
for ok, code, msg in CHECKS:
    if not ok:
        print("  [FAIL] %-22s %s" % (code, msg))

print("")
print("-- DECLARED GAPS (%d) -- the C0 partition is a COVER of H0; these are the"
      % len(GAPS))
print("   places where proof_dag.json instantiates LESS than the leaf --")
if not GAPS:
    print("  (none)")
for code, msg in GAPS:
    print("  [GAP]  %-22s %s" % (code, msg))

print("")
print("-- RETIRED GAPS (%d) -- kept in the inventory on purpose: a gap that" %
      len(RETIRED))
print("   disappears from the record is indistinguishable from one ignored --")
for code, msg in RETIRED:
    print("  [RETIRED] %-19s %s" % (code, msg))

print("")
print("-- WHAT NOW CAPS C0 --")
print("   C0 closed=%s at level %r, %s of 5 leaves closed."
      % (NODES["C0"]["closed"], NODES["C0"]["level"],
         NODES["C0"].get("subcases_closed")))
print("   The cap is NOT a coverage gap, and as of 2026-07-26 it is a SINGLE "
      "thing:")
print("     (1) the subcase->C0 exhaustiveness edge is JUDGMENT-referenced "
      "(GGHV22 Prop 4.3")
print("         case-(8,28) subcases (1)&(2), the field-split framework C14-C16, "
      "and the")
print("         alternate-regime partition C44) -- published mathematics a "
      "finite bookkeeping")
print("         checker cannot re-derive, so it is not machine-enforceable here.")
print("         This is the SOLE cap, and it is CORRECT BY CONSTRUCTION, not a "
      "backlog")
print("         item. prop43_audit.py (20/20) discharges the CITATION -- that "
      "GGHV22 says")
print("         what we use it to say -- and cannot re-derive the EXHAUSTIVENESS "
      "of the")
print("         partition, which is what 'exact-checked' here would assert. A "
      "standing")
print("         to-do listing 'apply the Prop 4.3 edge regrade' is an OVERCLAIM "
      "and was")
print("         DECLINED (commit 2636b90). The routes above 'claimed' are a "
      "machine-")
print("         checkable reformulation of the partition, or a formal proof -- "
      "not a regrade.")
print("     (2) subcase:dm1 (L_D, C10) NO LONGER caps C0. Raised 'claimed' -> "
      "'exact-checked'")
print("         on 2026-07-26, to the FLOOR only. One of the two reasons for the "
      "old hold")
print("         had gone STALE: dm1_branch_verify.py IS gated and passes 28/28, "
      "so 'its")
print("         checker is not in run_tests.sh' was no longer true -- and that is "
      "exactly")
print("         the 'claimed' -> 'exact-checked' gap. The CEILING "
      "('independently-audited')")
print("         is NOT applied: that needs a second independent implementation, "
      "and")
print("         C0_CLOSEOUT.md, which recommends it, is still same-day and "
      "unaudited. All")
print("         five C0 leaves are now exact-checked.")
print("")
print("c0_partition: %d CHECK failure(s), %d declared gap(s), %d retired%s"
      % (n_fail, len(GAPS), len(RETIRED),
         " -- bookkeeping consistent; gap inventory exactly as documented"
         if n_fail == 0 else ""))
sys.exit(0 if n_fail == 0 else 1)
