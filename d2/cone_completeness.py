#!/usr/bin/env python3
"""
cone_completeness.py -- JUDGMENT EDGE 1: cone-lemma branch completeness.

WHAT THIS IS
------------
`proof_dag.json` records the `branch -> subcase` exhaustiveness edge as
machine-checkable *for the count only* (`survivor + killed ==
open_branches_processed`).  The claim that the enumerated branch universe is
COMPLETE -- that no (a, b, T) case escaped it -- rests on `CASCADE_CONE_LEMMAS*.md`
and is JUDGMENT-referenced (`PROOF_DAG.md`, "Machine-checked vs
judgment-referenced edges").

This file closes that edge by building the branch universe a SECOND time,
directly from the stated mathematics (premises P1-P8 below), and comparing
branch-key SETS exactly against the engine's universe.

INDEPENDENCE STATEMENT (important -- this is what makes the check worth
anything).  The author of this file did NOT read `cascade_engine.py`,
`split_place_ledger.py`, `split_place_ledger_sub1.py`, `cone_lemmas.py`, or any
of the `audit_*.py` auditors.  The enumeration below was derived from:

  * SPLIT_PLACE_LEDGER.md / SPLIT_PLACE_LEDGER_SUB1.md  (the stated identities
    and caps, prose only -- NOT the JSON row data),
  * CASCADE_ENGINE_PLAN.md sections 1-3 (the h-ladder and the global degree
    coupling),
  * CAPS_AUDIT.md section 3 (the envelope caps),
  * ALT_REGIME.md (the v = 30 - 3a regime split),
  * AUDIT.md / CURRENT_STATUS.md (the terminal trichotomy and the T3 theorem).

The generated-artifact JSON files are read ONLY as the comparison target, never
as an input to the generator.  The generator's own numbers (327, 197, 246, 81,
654; 1333, 1007, 1171, 136, 2614, 2178) are produced before any artifact is
opened.

READ-ONLY.  Writes nothing.  `--quiet` prints only failures and the verdict;
exit 0 iff every check passes.

Usage:  python -u cone_completeness.py [--quiet] [--certs]
"""
import json
import os
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
QUIET = "--quiet" in sys.argv
SHOW_CERTS = "--certs" in sys.argv


def out(*a):
    if not QUIET:
        print(*a)


FAILURES = []


def fail(code, msg):
    FAILURES.append((code, msg))


# ===========================================================================
# PART A -- THE PREMISES  (each carries its source; none is derived here)
# ===========================================================================
#
# P1  FIELD SPLIT.  After the base change of C14-C16 the quartic q factors into
#     four distinct GEOMETRIC (degree-one) places q = p1 p2 p3 p4, and t is a
#     further, separate place.  Because each p_i has residue degree 1, a
#     valuation of v units at p_i costs exactly v units of polynomial degree.
#         source: CURRENT_STATUS.md C14-C16; t5_split_place_verify.py
#
# P2  STRATUM COORDINATES.  For the state polynomial e = d_{-1} put
#         a    := v_t(e)                >= 0
#         b_i  := v_{p_i}(e / t^a)      >= 0
#     The four places are interchangeable (S4 acts by permuting the roots of q),
#     so the case datum is the S4-ORBIT of (b_1,..,b_4), canonically the sorted
#     decreasing tuple.  a and b are determined by the solution, so distinct
#     (a,b) are mutually exclusive: the stratification is a PARTITION.
#         source: SPLIT_PLACE_LEDGER.md header; CASCADE_ENGINE_PLAN.md sec. 3
#
# P3  STRATUM BUDGET.  t and the p_i are distinct places, all of degree one, so
#         a + sum_i b_i  <=  deg e.
#     Envelope cap: deg e = deg d_{-1} <= 2*5 = 10 (sub2), <= 3*5 = 15 (sub1).
#         source: CAPS_AUDIT.md sec. 3 table (`dm1 = e`, k = 5, lambda = 2/3);
#                 CURRENT_STATUS.md C6 (envelope_bounds_verify.py)
#
# P4  TERMINAL TRICHOTOMY.  The terminal branch label is the zero/nonzero
#     alternative on the pair (d1, sigma):
#         T1 :  d1 != 0
#         T2 :  d1 == 0,  sigma != 0
#         T3 :  d1 == 0,  sigma == 0
#     Exhaustive and pairwise disjoint BY CONSTRUCTION -- it is the complete
#     truth table of two predicates, with the (d1 != 0) rows merged because the
#     T1 terminal identity does not use sigma.
#         source: SPLIT_PLACE_LEDGER.md ("T3 (d1=0, sigma=0)");
#                 ALT_INF_SWEEP.md judgment J2 ("For T2 (d1==0), h7 = 8192 d1^2 == 0")
#
# P5  T3 IS EMPTY.  The field-stable split-place sigma-locus theorem kills T3
#     globally in characteristic 0.
#         source: CURRENT_STATUS.md C14-C16; FIELD_SPLIT_AUDIT.md
#     [PREMISE -- consumed, not re-derived here.]
#
# P6  TOP-OF-LADDER IDENTITIES.  The cascade ladder extracted from
#     `f31_graded.txt` has
#         h7 = 8192 d1^2,
#         h6 = -3072 sigma^2 + 14336 d1^2 d2 + 8192 d1 e,   so  h6|_{d1=0} = -3072 sigma^2.
#     At a finite place p_i the top relation E^3 g_l = (c q)^l h_l gives, with
#     v_i(E) = b_i and v_i(q) = 1,
#         T1 (level 7):  3 b_i + v_i(g7) = 7 + 2 v_i(d1)
#         T2 (level 6):  3 b_i + v_i(g6) = 6 + 2 v_i(sigma)
#         source: CASCADE_ENGINE_PLAN.md sec. 1 and 2 (ladder + transition);
#                 SPLIT_PLACE_LEDGER{,_SUB1}.md "Exact terminal pruning"
#
# P7  GLOBAL DEGREE CAPS (the coupling of the four places).
#         sub2:  deg d1 <= 6,  deg sigma <= 8,  deg g_l <= 10 + 3a
#         sub1:  deg d1 <= 9,  deg sigma <= 12, deg g7 <= 46, deg g6 <= 48
#     sigma = 4 d0 - d2^2, hence deg sigma <= max(deg d0, 2 deg d2) = 8 (sub2)
#     / 12 (sub1), consistent with the CAPS_AUDIT lambda*k rows.
#         source: SPLIT_PLACE_LEDGER{,_SUB1}.md; CURRENT_STATUS.md C27
#                 ("independent correction of the 15+3a cap")
#
# P8  REGIME SPLIT.  With Phi~ = t^30 u, e = t^a E one has v = 30 - 3a, so
#         v >= 0  <=>  a <= 10   (STANDARD reduction; the ladder above applies)
#         v <  0  <=>  a >= 11   (ALTERNATE reduction; flipped/descending
#                                 cascade, different terminal condition)
#     For sub2 P3 already forces a <= 10, so sub2 has no alternate regime.
#     For sub1 the alternate regime is a in [11,15].
#         source: ALT_REGIME.md sec. "Flipped reduction"
#
# ---------------------------------------------------------------------------

WINDOWS = {
    "sub2": dict(deg_e=10, deg_d1=6, deg_sigma=8,
                 g7_cap=lambda a: 10 + 3 * a, g6_cap=lambda a: 10 + 3 * a),
    "sub1": dict(deg_e=15, deg_d1=9, deg_sigma=12,
                 g7_cap=lambda a: 46, g6_cap=lambda a: 48),
}
STANDARD_A_MAX = 10          # P8


# ===========================================================================
# PART B -- THE INDEPENDENT GENERATOR
# ===========================================================================

def sorted_compositions(s, k=4):
    """All non-increasing k-tuples of non-negative integers summing to s.
    (= the S4-orbit representatives of P2, i.e. partitions of s into <= 4 parts)."""
    res = []

    def rec(rem, slots, cap, acc):
        if slots == 0:
            if rem == 0:
                res.append(tuple(acc))
            return
        for v in range(min(rem, cap), -1, -1):
            if v * slots < rem:      # cannot reach rem with non-increasing tail
                break
            rec(rem - v, slots - 1, v, acc + [v])

    rec(s, k, s, [])
    return res


def stratum_universe(deg_e_cap):
    """P2 + P3: every (a, sorted b) with a + sum b <= deg_e_cap."""
    U = []
    for a in range(0, deg_e_cap + 1):
        for s in range(0, deg_e_cap - a + 1):
            for b in sorted_compositions(s):
                U.append((a, b))
    return U


def local_minima(b, c):
    """P6 at ONE place.

    The identity is  3b + G = c + 2x  with G = v_i(g), x = v_i(aux) and
    G, x non-negative integers.  Solving, G = c + 2x - 3b, so

        x >= x_min(b) := max(0, ceil((3b - c) / 2))

    and G is strictly increasing in x, so BOTH v_i(aux) and v_i(g) attain their
    minima at the same x.  Hence there is a well-defined pair of minimal local
    orders and no trade-off between the two budgets.
    """
    x_min = max(0, -((-(3 * b - c)) // 2))    # ceil((3b-c)/2) for ints
    g_min = c + 2 * x_min - 3 * b
    assert g_min >= 0 and (3 * b + g_min) == (c + 2 * x_min)
    return x_min, g_min


def terminal_certificate(a, b, T, win):
    """P6 + P7: decide (a, b, T) and return a LINEAR-BUDGET CERTIFICATE.

    Returns (feasible: bool, cert: dict).  When infeasible the certificate is a
    single violated linear inequality with all its numbers, so the exclusion is
    independently checkable by hand.
    """
    W = WINDOWS[win]
    if T == "T1":
        c, aux, aux_cap, g_name, g_cap = 7, "d1", W["deg_d1"], "g7", W["g7_cap"](a)
    elif T == "T2":
        c, aux, aux_cap, g_name, g_cap = 6, "sigma", W["deg_sigma"], "g6", W["g6_cap"](a)
    else:
        raise ValueError(T)

    pairs = [local_minima(bi, c) for bi in b]
    aux_sum = sum(p[0] for p in pairs)
    g_sum = sum(p[1] for p in pairs)

    if aux_sum > aux_cap:
        return False, dict(kind="budget", dim=aux, forced=aux_sum, cap=aux_cap,
                           orders=[p[0] for p in pairs],
                           ineq="sum_i v_i(%s) >= %d > %d = deg %s cap"
                                % (aux, aux_sum, aux_cap, aux))
    if g_sum > g_cap:
        return False, dict(kind="budget", dim=g_name, forced=g_sum, cap=g_cap,
                           orders=[p[1] for p in pairs],
                           ineq="sum_i v_i(%s) >= %d > %d = deg %s cap"
                                % (g_name, g_sum, g_cap, g_name))
    return True, dict(kind="feasible", aux=aux, aux_forced=aux_sum, aux_cap=aux_cap,
                      g=g_name, g_forced=g_sum, g_cap=g_cap)


def generate(win):
    """The independent branch universe for one window.

    Returns dict with:
      strata          : all (a,b) of P2/P3
      standard        : strata with a <= 10   (P8)
      alternate       : strata with a >= 11   (P8; sub1 only)
      terminal_open   : set of (a,b,T) surviving the P6/P7 terminal test
      alt_open        : set of (a,b,T) in the alternate regime (terminal test
                        NOT applicable -- P8)
      certs           : (a,b,T) -> certificate dict
    """
    W = WINDOWS[win]
    strata = stratum_universe(W["deg_e"])
    standard = [s for s in strata if s[0] <= STANDARD_A_MAX]
    alternate = [s for s in strata if s[0] > STANDARD_A_MAX]
    open_keys, certs = set(), {}
    for (a, b) in standard:
        for T in ("T1", "T2"):
            ok, cert = terminal_certificate(a, b, T, win)
            certs[(a, b, T)] = cert
            if ok:
                open_keys.add((a, b, T))
    alt_keys = set((a, b, T) for (a, b) in alternate for T in ("T1", "T2"))
    return dict(strata=strata, standard=standard, alternate=alternate,
                terminal_open=open_keys, alt_open=alt_keys, certs=certs)


# ===========================================================================
# PART C -- PINNED SELF-CONSISTENCY OF THE GENERATOR
# (numbers produced before any artifact is read; they are also the numbers the
#  ledger DOCS publish, so a drift on either side is loud)
# ===========================================================================
PINNED = {
    "sub2": dict(strata=327, T1=197, T2=246, dead_both=81, terminal_decisions=654,
                 alt_strata=0),
    "sub1": dict(strata=1333, T1=1007, T2=1171, dead_both=136, terminal_decisions=2614,
                 alt_strata=26),
}

GEN = {w: generate(w) for w in ("sub2", "sub1")}

out("=" * 74)
out("JUDGMENT EDGE 1 -- CONE-LEMMA BRANCH COMPLETENESS")
out("independent generator (premises P1-P8) vs the engine's branch universe")
out("=" * 74)

for w in ("sub2", "sub1"):
    g = GEN[w]
    n1 = sum(1 for k in g["terminal_open"] if k[2] == "T1")
    n2 = sum(1 for k in g["terminal_open"] if k[2] == "T2")
    dead = sum(1 for (a, b) in g["standard"]
               if (a, b, "T1") not in g["terminal_open"]
               and (a, b, "T2") not in g["terminal_open"])
    p = PINNED[w]
    got = dict(strata=len(g["strata"]), T1=n1, T2=n2, dead_both=dead,
               terminal_decisions=2 * len(g["standard"]),
               alt_strata=len(g["alternate"]))
    out("")
    out("-- %s : independently generated census --" % w)
    for k in ("strata", "T1", "T2", "dead_both", "terminal_decisions", "alt_strata"):
        mark = "OK " if got[k] == p[k] else "!! "
        out("   %s%-20s %5d   (published %d)" % (mark, k, got[k], p[k]))
        if got[k] != p[k]:
            fail("GEN-CENSUS-%s" % w.upper(),
                 "%s/%s: generator says %d, ledger doc publishes %d"
                 % (w, k, got[k], p[k]))
    # T1 subset T2 is a structural consequence worth stating out loud
    t1 = set(k[:2] for k in g["terminal_open"] if k[2] == "T1")
    t2 = set(k[:2] for k in g["terminal_open"] if k[2] == "T2")
    out("   %sT1-feasible strata subset of T2-feasible: %s"
        % ("OK " if t1 <= t2 else "!! ", t1 <= t2))


# ===========================================================================
# PART D -- THE ENGINE'S UNIVERSE (comparison target; artifacts read here only)
# ===========================================================================
def load(fn):
    with open(os.path.join(HERE, fn)) as f:
        return json.load(f)


def cones_keys(fn):
    return set((r["a_t"], tuple(r["b"]), r["branch"]) for r in load(fn)["branches"])


def ledger_open(fn):
    s = set()
    for r in load(fn)["strata"]:
        for T in r["open_branches"]:
            s.add((r["a_t"], tuple(r["b"]), T))
    return s


def ledger_rows(fn):
    return {(r["a_t"], tuple(r["b"])): r for r in load(fn)["strata"]}


ENGINE = {
    "sub2": cones_keys("cascade_cones.json"),                 # 420
    "sub1": cones_keys("cascade_cones_sub1_depth4.json"),     # 2178
}
LEDGER = {"sub2": ledger_open("split_place_ledger.json"),
          "sub1": ledger_open("split_place_ledger_sub1.json")}
ROWS = {"sub2": ledger_rows("split_place_ledger.json"),
        "sub1": ledger_rows("split_place_ledger_sub1.json")}

# proof_dag.json branch nodes (the object the judgment edge actually gates)
DAG = load("proof_dag.json")
DAGK = {}
for n in DAG["nodes"]:
    if n["type"] == "branch" and n.get("window") in ("sub2", "sub1"):
        DAGK.setdefault(n["window"], set()).add(
            (n["a_t"], tuple(n["b"]), n["branch"]))


# ===========================================================================
# PART E -- EXACT SET COMPARISON, BOTH DIRECTIONS
# ===========================================================================
out("")
out("=" * 74)
out("EXACT BRANCH-KEY SET COMPARISON")
out("=" * 74)

# The 23 sub2 branches the ledger removes BEYOND the terminal test.  These are
# not terminal-infeasible: each is excised by a NAMED earlier exact proof.  The
# check below does not take the removal on trust -- it requires (i) the row to
# be marked proven_infeasible, (ii) the reference to be a document reference
# rather than a restatement of the terminal test, and (iii) the referenced file
# to exist in the repository.
EXPECTED_EXCISIONS = 23

for w in ("sub2", "sub1"):
    g = GEN[w]
    indep_std = g["terminal_open"]
    eng = ENGINE[w]

    extra_in_engine = sorted(eng - indep_std)
    missing_from_engine = sorted(indep_std - eng)

    out("")
    out("-- %s --" % w)
    out("   independent terminal-feasible branches : %d" % len(indep_std))
    out("   engine branch universe                 : %d" % len(eng))
    out("   engine \\ independent  (SPURIOUS -- would be a soundness break): %d"
        % len(extra_in_engine))
    out("   independent \\ engine  (EXCISED -- must each carry a proof)   : %d"
        % len(missing_from_engine))

    if extra_in_engine:
        fail("SPURIOUS-BRANCH-%s" % w.upper(),
             "%d branch key(s) processed by the engine lie OUTSIDE the "
             "independently derived universe: %s"
             % (len(extra_in_engine), extra_in_engine[:20]))

    # every excised branch must resolve to a named, existing proof document
    for (a, b, T) in missing_from_engine:
        row = ROWS[w].get((a, b))
        if row is None:
            fail("EXCISION-NOROW-%s" % w.upper(),
                 "a%d_b%s_%s is in the independent universe but has no ledger row"
                 % (a, "".join(map(str, b)), T))
            continue
        br = row["branches"].get(T, {})
        status, ref = br.get("status"), br.get("reference", "")
        docs = [tok.strip().rstrip(",") for tok in ref.replace(",", " ").split()
                if tok.strip().rstrip(",").endswith(".md")]
        resolved = [d for d in docs if os.path.exists(os.path.join(HERE, d))]
        bad = (status != "proven_infeasible") or (not docs) or (len(resolved) != len(docs))
        if bad:
            fail("EXCISION-UNJUSTIFIED-%s" % w.upper(),
                 "a%d_b%s_%s excised from the engine universe with status=%r "
                 "reference=%r (docs found: %s)"
                 % (a, "".join(map(str, b)), T, status, ref, resolved))

    if missing_from_engine:
        out("")
        out("   EXCISION TABLE -- branches the mathematics admits but the engine")
        out("   universe drops, each with the proof that authorises the drop:")
        by_ref = {}
        for (a, b, T) in missing_from_engine:
            row = ROWS[w][(a, b)]
            ref = row["branches"][T].get("reference", "?")
            by_ref.setdefault(ref, []).append(
                "a%d_b%s_%s" % (a, "".join(map(str, b)), T))
        for ref in sorted(by_ref):
            out("     %-58s %s" % (ref, " ".join(sorted(by_ref[ref]))))

    # SCOPE COMPLIANCE of the excision layer.  PROOF_INVENTORY.md C35 declares
    # the pre-repair strata kills (T5_STRATA_50_11, T5_60_T1/T2, T5_T1_AQ12,
    # T5_STRATUM_10_0) valid ONLY as "geometrically-q-coprime / uniform-q^r"
    # statements after the field-split repair, and FIELD_SPLIT_AUDIT's a=7
    # theorem is likewise a geometrically-q-coprime statement.  Both scopes mean
    # exactly: b_1 = b_2 = b_3 = b_4 (b = 0000 is q-coprime, b = rrrr is
    # uniform q^r).  So the excisions are in scope IFF every excised branch has
    # uniform b -- a machine-checkable side condition on a judgment edge.
    nonuniform = [k for k in missing_from_engine if len(set(k[1])) != 1]
    if missing_from_engine:
        out("   excision scope check (C35: geometrically-q-coprime / uniform-q^r")
        out("   => b_1 = b_2 = b_3 = b_4): %d/%d excised branches have uniform b"
            % (len(missing_from_engine) - len(nonuniform), len(missing_from_engine)))
    if nonuniform:
        fail("EXCISION-OUT-OF-SCOPE-%s" % w.upper(),
             "%d excised branch(es) have NON-uniform b, outside the scope "
             "PROOF_INVENTORY C35 declares valid for the pre-repair strata "
             "proofs: %s" % (len(nonuniform), nonuniform))

    if w == "sub2" and len(missing_from_engine) != EXPECTED_EXCISIONS:
        fail("EXCISION-COUNT-SUB2",
             "sub2 excision list is %d branches, expected %d -- the excision "
             "layer changed and must be re-reviewed"
             % (len(missing_from_engine), EXPECTED_EXCISIONS))
    if w == "sub1" and missing_from_engine:
        fail("EXCISION-COUNT-SUB1",
             "sub1 was expected to need NO excisions beyond the terminal test; "
             "found %d" % len(missing_from_engine))

    # ledger vs engine vs DAG must be the same set of keys
    if LEDGER[w] - set(k for k in LEDGER[w] if k[0] > STANDARD_A_MAX) != eng:
        fail("LEDGER-ENGINE-%s" % w.upper(),
             "the ledger's standard-regime open branches and the engine's "
             "processed branches are different sets")
    if DAGK.get(w) != eng:
        fail("DAG-ENGINE-%s" % w.upper(),
             "proof_dag.json branch nodes for %s do not equal the engine "
             "universe (dag %d, engine %d)" % (w, len(DAGK.get(w, ())), len(eng)))


# ---------------------------------------------------------------------------
# the sub1 ALTERNATE regime: universe check only (P8 -- the terminal test does
# not apply there; the 52 -> 27 reduction is C44's, consumed not re-derived)
# ---------------------------------------------------------------------------
alt_indep = GEN["sub1"]["alt_open"]
alt_ledger = set(k for k in LEDGER["sub1"] if k[0] > STANDARD_A_MAX)
out("")
out("-- sub1 ALTERNATE regime (a in [11,15], v = 30-3a < 0) --")
out("   independent universe : %d branches over %d strata"
    % (len(alt_indep), len(GEN["sub1"]["alternate"])))
out("   ledger universe      : %d branches" % len(alt_ledger))
out("   sets equal           : %s" % (alt_indep == alt_ledger))
if alt_indep != alt_ledger:
    fail("ALT-UNIVERSE",
         "alternate-regime branch universes differ: only-mine=%s only-ledger=%s"
         % (sorted(alt_indep - alt_ledger), sorted(alt_ledger - alt_indep)))

alt_sweep = load("alt_inf_sweep.json")
alt_open27 = set((r["a"], tuple(r["b"]), r["branch"]) for r in alt_sweep["branches"])
out("   post-C44 open (alt_inf_sweep.json)      : %d" % len(alt_open27))
out("   killed by C44 (consumed, not re-derived): %d" % (len(alt_indep) - len(alt_open27)))
if not alt_open27 <= alt_indep:
    fail("ALT-SWEEP-OUTSIDE",
         "alt_inf_sweep branches outside the independent alternate universe: %s"
         % sorted(alt_open27 - alt_indep))


# ---------------------------------------------------------------------------
# T3: present in the raw case tree of P4, excised globally by P5
# ---------------------------------------------------------------------------
out("")
out("-- terminal branch T3 (d1 == 0 and sigma == 0) --")
t3_rows = 0
t3_ok = 0
for w in ("sub2", "sub1"):
    for row in ROWS[w].values():
        t3 = row.get("T3")
        if t3 is not None:
            t3_rows += 1
            if t3.get("status") == "proven_infeasible":
                t3_ok += 1
out("   ledger rows carrying an explicit T3 verdict : %d" % t3_rows)
out("   of which proven_infeasible                  : %d" % t3_ok)
out("   (P5: the split-place sigma-locus theorem, C14-C16 -- CONSUMED premise,")
out("    not re-derived here.)")
if t3_rows != t3_ok:
    fail("T3-NOT-CLOSED", "%d ledger rows carry a T3 verdict that is not "
         "proven_infeasible" % (t3_rows - t3_ok))
if t3_rows != len(ROWS["sub2"]) + len(ROWS["sub1"]):
    fail("T3-MISSING-ROWS",
         "only %d of %d strata rows carry an explicit T3 verdict"
         % (t3_rows, len(ROWS["sub2"]) + len(ROWS["sub1"])))


# ---------------------------------------------------------------------------
# optional: dump the exclusion certificates
# ---------------------------------------------------------------------------
if SHOW_CERTS and not QUIET:
    print("")
    print("=" * 74)
    print("EXCLUSION CERTIFICATES (one violated linear budget per excluded branch)")
    print("=" * 74)
    for w in ("sub2", "sub1"):
        g = GEN[w]
        hist = Counter()
        for k, c in g["certs"].items():
            if c["kind"] == "budget":
                hist[(k[2], c["dim"])] += 1
        print("  %s: %d branch keys excluded by a linear budget certificate"
              % (w, sum(hist.values())))
        for (T, dim), n in sorted(hist.items()):
            print("     %-4s %-8s %5d" % (T, dim, n))
        ex = [(k, c) for k, c in sorted(g["certs"].items()) if c["kind"] == "budget"][:4]
        for k, c in ex:
            print("     e.g. a%d_b%s_%s : %s   (min local orders %s)"
                  % (k[0], "".join(map(str, k[1])), k[2], c["ineq"], c["orders"]))


# ===========================================================================
# VERDICT
# ===========================================================================
out("")
out("=" * 74)
for code, msg in FAILURES:
    print("  [FAIL] %-28s %s" % (code, msg))
print("")
print("cone_completeness: %d failure(s)%s"
      % (len(FAILURES),
         " -- engine branch universe EQUALS the independently generated "
         "universe minus an enumerated, individually proof-referenced "
         "excision list" if not FAILURES else ""))
sys.exit(0 if not FAILURES else 1)
