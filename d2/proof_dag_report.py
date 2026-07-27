#!/usr/bin/env python3
"""
proof_dag_report.py  --  reads proof_dag.json, emits:

  1. CLOSURE CENSUS by level (per node type).
  2. WEAKEST-EDGE list: the exhaustiveness edges whose upgrade would promote the
     most downstream *closed* obligation mass -- plus a grouped node-evidence
     AUDIT-PRIORITY QUEUE (the actionable version).
  3. INCONSISTENCY findings: any numeric claim in CURRENT_STATUS.md or
     FRONTIER_V2.md that is STRONGER than the DAG supports.  Loud.

Exit 0 iff zero inconsistencies.  --quiet suppresses the census/edge sections
and prints only inconsistencies (and the final verdict line).

READ-ONLY.  Writes nothing.
"""
import json, os, re, sys
from collections import defaultdict, Counter

HERE = os.path.dirname(os.path.abspath(__file__))
QUIET = "--quiet" in sys.argv

# --- optional inputs (defaults unchanged: the committed DAG + FRONTIER_V2) ---
def _argval(flag, default):
    if flag in sys.argv:
        return sys.argv[sys.argv.index(flag) + 1]
    return default

DAG_FILE = _argval("--dag", "proof_dag.json")
FRONTIER_FILE = _argval("--frontier", "FRONTIER_V2.md")

def load_json(fn):
    with open(os.path.join(HERE, fn)) as f:
        return json.load(f)

def load_text(fn):
    p = os.path.join(HERE, fn)
    if not os.path.exists(p):
        return ""
    with open(p, encoding="utf-8") as f:
        return f.read()

DAG = load_json(DAG_FILE)
NODES = {n["id"]: n for n in DAG["nodes"]}
LEVELS = DAG["levels"]
LRANK = {l: i for i, l in enumerate(LEVELS)}

def out(*a):
    if not QUIET:
        print(*a)

# ---------------------------------------------------------------------------
# child map (structural edges only; certificate edges point state->cert and are
# excluded from the closure tree)
# ---------------------------------------------------------------------------
CHILDREN = defaultdict(list)
for e in DAG["edges"]:
    if NODES.get(e["child"], {}).get("type") == "certificate":
        continue
    CHILDREN[e["parent"]].append(e["child"])

def is_closed_leaf(node):
    if node["type"] == "state":
        return True
    if node["type"] == "branch" and node.get("cascade_status") not in (None, "survives") \
            and node.get("window") in ("sub2", "sub1"):
        return True   # engine-killed branch: a closed leaf with no state children
    return False

_MASS = {}
def downstream_mass(nid):
    """# of closed-obligation leaves in the subtree rooted at nid."""
    if nid in _MASS:
        return _MASS[nid]
    node = NODES[nid]
    kids = CHILDREN.get(nid, [])
    if is_closed_leaf(node) and not kids:
        _MASS[nid] = 1; return 1
    tot = sum(downstream_mass(c) for c in kids)
    if is_closed_leaf(node):
        tot += 1
    _MASS[nid] = tot
    return tot

for nid in NODES:
    downstream_mass(nid)

# ===========================================================================
# 1. CLOSURE CENSUS
# ===========================================================================
out("=" * 72)
out("COVERAGE PROOF-DAG REPORT")
out("=" * 72)
out("nodes: %d | edges: %d" % (DAG["counts"]["nodes"], DAG["counts"]["edges"]))
out("source_sha256: %s" % DAG["provenance"]["source_sha256"][:16])
out("")
out("-- CLOSURE CENSUS (node level = evidence strength) --")
for t in sorted(DAG["closure_census"]):
    cen = DAG["closure_census"][t]
    parts = ", ".join("%s=%d" % (l, cen[l]) for l in LEVELS if l in cen)
    out("  %-12s %s" % (t, parts))

c0 = NODES["C0"]
# C0's children are the leaves of the case partition (JUDGMENT_EDGES.md sec.3);
# read them off the graph rather than hardcoding, so a repair to the child list
# cannot leave this report quietly stale.
C0_CHILDREN = sorted(set(e["child"] for e in DAG["edges"] if e["parent"] == "C0"))
out("")
out("  TARGET C0: closed=%s level=%s (subcases closed %d/%d)" % (
    c0["closed"], c0["level"], c0.get("subcases_closed", 0), len(C0_CHILDREN)))
for sid in C0_CHILDREN + ["subcase:sub1_alt_defect0"]:
    s = NODES.get(sid)
    if s is None:
        continue
    out("    %-26s closed=%s level=%-14s %s" % (
        sid.split(":")[1], s["closed"], s["level"],
        ("branches_closed %d/%d" % (s.get("branches_closed", 0), s.get("n_branches", 0))
         if "n_branches" in s else s.get("closure_note", ""))))

# state level by window
byw = Counter()
for n in DAG["nodes"]:
    if n["type"] == "state":
        byw[(n["window"], n["level"])] += 1
out("")
out("  killed states by window x level:")
for w in ("sub2", "sub1", "altdefect0", "corner"):
    row = {l: byw.get((w, l), 0) for l in LEVELS}
    out("    %-11s %s" % (w, ", ".join("%s=%d" % (l, row[l]) for l in LEVELS if row[l])))

# ===========================================================================
# 2. WEAKEST EDGES + AUDIT-PRIORITY QUEUE
# ===========================================================================
out("")
out("-- WEAKEST EXHAUSTIVENESS EDGES (cap <= claimed), by downstream closed mass --")
weak = []
for e in DAG["edges"]:
    if NODES.get(e["child"], {}).get("type") == "certificate":
        continue
    exl = e["exhaustiveness_level"]
    if LRANK[exl] <= LRANK["claimed"]:
        weak.append((downstream_mass(e["child"]), e))
weak.sort(key=lambda x: (-x[0], x[1]["parent"], x[1]["child"]))
for mass, e in weak[:12]:
    out("  mass=%-5d %-14s %s -> %s" % (mass, "[" + e["exhaustiveness_level"] + "]",
        e["parent"], e["child"]))
    out("            exhaustiveness: %s" % e["exhaustiveness_ref"][:96])

out("")
out("-- AUDIT-PRIORITY QUEUE (grouped node-evidence upgrades, by mass) --")
# engine-killed branches split by evidence level after the cascade-audit JOIN
ek = defaultdict(int)          # residual: engine-killed still at 'claimed'
ek_aud = defaultdict(int)      # engine-killed promoted to >= independently-audited
for n in DAG["nodes"]:
    if n["type"] == "branch" and n.get("window") in ("sub2", "sub1") \
            and n.get("cascade_status") not in (None, "survives"):
        if LRANK[n["level"]] >= LRANK["independently-audited"]:
            ek_aud[n["window"]] += 1
        else:
            ek[n["window"]] += 1
prio = []
if ek:
    prio.append((sum(ek.values()),
        "engine-killed branches BELOW independently-audited (sub2 %d, sub1 %d): "
        "killed only by the t/inf layer, outside the depth-4 q-cascade auditor's "
        "scope.  The audit_inf_kills.json join (C43) already carries them to "
        "exact-checked; independently-audited additionally needs the q+t_rl "
        "narrowing audited -- audit_tplace_cases.py reads the kills-OFF q+t "
        "artifact, not the _rl one" %
        (ek.get("sub2", 0), ek.get("sub1", 0))))
# group 2: exact-checked states (ledger 'AUDITED' = same-author) -> independent audit
ec = Counter()
for n in DAG["nodes"]:
    if n["type"] == "state" and n["level"] == "exact-checked":
        ec[n["window"]] += 1
if ec:
    prio.append((sum(ec.values()),
        "exact-checked killed states (sub2 %d, sub1 %d): ledger-'AUDITED' but "
        "same-author -> independent audit raises to independently-audited"
        % (ec.get("sub2", 0), ec.get("sub1", 0))))
# group 3: claimed states (pending msolve/GB) -> exact check / cert
cl = Counter()
for n in DAG["nodes"]:
    if n["type"] == "state" and n["level"] == "claimed":
        cl[n["window"]] += 1
if cl:
    prio.append((sum(cl.values()),
        "claimed killed states (%s): msolve/GB/ambiguous pending -> need an exact "
        "re-check or object certificate" % ", ".join("%s %d" % kv for kv in sorted(cl.items()))))
# group 4: NOT-YET-CERTIFICATED certificates
nyc = sum(1 for n in DAG["nodes"] if n["type"] == "certificate" and not n["found"])
if nyc:
    prio.append((nyc, "object certificates NOT-YET-CERTIFICATED (%d): lift/GB "
                 "timeouts -> re-run consuming auditor to reach 'certified'" % nyc))
prio.sort(key=lambda x: -x[0])
for mass, desc in prio:
    out("  [%5d] %s" % (mass, desc))

# ===========================================================================
# 3. INCONSISTENCY FINDINGS  (doc/corpus claim STRONGER than DAG supports)
# ===========================================================================
FINDINGS = []       # (severity, code, message)
def inconsistency(code, msg): FINDINGS.append(("INCONSISTENCY", code, msg))
def gap(code, msg):          FINDINGS.append(("COVERAGE-GAP", code, msg))
def ok(code, msg):           FINDINGS.append(("OK", code, msg))

FR = load_text(FRONTIER_FILE)
CS = load_text("CURRENT_STATUS.md")

# DAG-supported counts
def state_ge(window, level):
    return sum(1 for n in DAG["nodes"] if n["type"] == "state"
               and n["window"] == window and LRANK[n["level"]] >= LRANK[level])

dag_audit_sub2 = state_ge("sub2", "independently-audited")
dag_audit_sub1 = state_ge("sub1", "independently-audited")

# --- I1: FRONTIER must SEPARATE same-author exact-checked from independently- ---
#     audited (the DAG levels), not lump both under a single "audited" column.
# (a) the misleading bare "Killed (audited)" header must be gone, replaced by the
#     honest "exact-checked, same-author" label; (b) FRONTIER's evidence-grade
#     table's ">= independently-audited" per-window counts must MATCH the DAG
#     (a staleness / consistency guard).
honest_label = ("Killed (exact-checked, same-author)" in FR
                and "Killed (audited)" not in FR)
gtab = FR.split("Killed-state evidence grade", 1)[-1] if \
    "Killed-state evidence grade" in FR else ""
grow = re.compile(r"\|\s*%s\s*\|\s*\d+\s*\|\s*\d+\s*\|\s*\d+\s*\|\s*\d+\s*\|\s*(\d+)\s*\|")
gm2 = re.search(grow.pattern % "sub2", gtab)
gm1 = re.search(grow.pattern % "sub1", gtab)
if not honest_label:
    inconsistency("FRONTIER-AUDITED-LABEL",
        "FRONTIER still labels same-author exact checks as 'Killed (audited)', "
        "which reads as an independent audit. It must be split into "
        "'exact-checked (same-author)' vs 'independently-audited' (the DAG levels).")
elif not (gm2 and gm1):
    inconsistency("FRONTIER-AUDITED-LABEL",
        "FRONTIER lacks a parseable evidence-grade table with the DAG's "
        "'>= independently-audited' per-window counts.")
elif int(gm2.group(1)) != dag_audit_sub2 or int(gm1.group(1)) != dag_audit_sub1:
    inconsistency("FRONTIER-AUDITED-LABEL",
        "FRONTIER evidence-grade '>= independently-audited' counts (sub2=%s, "
        "sub1=%s) disagree with the DAG (sub2=%d, sub1=%d) -- FRONTIER_V2.md is "
        "stale; regenerate `python frontier_rollup.py` after `python proof_dag.py`."
        % (gm2.group(1), gm1.group(1), dag_audit_sub2, dag_audit_sub1))
else:
    ok("FRONTIER-AUDITED-LABEL",
       "FRONTIER separates same-author exact-checked from independently-audited; "
       "its '>= independently-audited' counts (sub2=%d, sub1=%d) match the DAG."
       % (dag_audit_sub2, dag_audit_sub1))

# --- I2: CURRENT_STATUS S1a cascade INDEPENDENT-AUDIT claim vs the DAG join ----
# After the machine-join of audit_cascade_kills{,_sub1}.py, engine-killed branches
# the auditor confirms are 'independently-audited'.  The doc (C18 '390 killed',
# C29 '1899 killed', spec-only audited) is supported iff the DAG independently-
# audits at least those counts.  (Branches killed only by the t/inf layer stay
# 'claimed'; the doc does not count them under C18/C29, so they are not a shortfall.)
c18 = re.search(r"390 killed", CS)
c29 = re.search(r"1899 killed", CS)
doc_sub2 = 390 if c18 else 0
doc_sub1 = 1899 if c29 else 0
dag_ek_ia_sub2 = ek_aud.get("sub2", 0)
dag_ek_ia_sub1 = ek_aud.get("sub1", 0)
shortfall = []
if doc_sub2 and dag_ek_ia_sub2 < doc_sub2:
    shortfall.append("sub2 (C18): doc %d, DAG independently-audited %d"
                     % (doc_sub2, dag_ek_ia_sub2))
if doc_sub1 and dag_ek_ia_sub1 < doc_sub1:
    shortfall.append("sub1 (C29): doc %d, DAG independently-audited %d"
                     % (doc_sub1, dag_ek_ia_sub1))
if shortfall:
    inconsistency("CASCADE-BRANCH-AUDIT",
        "CURRENT_STATUS S1a claims spec-only INDEPENDENT AUDIT of the cascade "
        "branch kills, but the DAG join supports fewer: %s. Re-run the auditors "
        "with --emit-artifact and rebuild the DAG." % "; ".join(shortfall))
else:
    ok("CASCADE-BRANCH-AUDIT",
       "cascade branch kills machine-joined from audit_cascade_kills{,_sub1}.py: "
       "DAG independently-audits sub2=%d (C18 claims %d), sub1=%d (C29 claims %d); "
       "%d engine-killed branches are t/inf-layer-only kills: exact-checked via "
       "the audit_inf_kills.json join (C43), not independently-audited."
       % (dag_ek_ia_sub2, doc_sub2, dag_ek_ia_sub1, doc_sub1,
          sum(ek.values())))

# --- I3: surviving-branch counts --------------------------------------------
# THIS FILE'S CONTRACT (see the module docstring) is one-directional: it fails on
# a doc claim that is STRONGER than the DAG supports.  A DAG that is stronger than
# the doc is a STALE DOC, not an inconsistency, and must not fail the build --
# otherwise closing a branch breaks the suite, which is exactly backwards.
# Before 2026-07-26 this check demanded equality with CURRENT_STATUS's 26/171 and
# would have hard-failed the moment the frontier closed.
open_sub2 = NODES["subcase:sub2"].get("branches_open")
open_sub1 = NODES["subcase:sub1"].get("branches_open")
DOC_OPEN_SUB2, DOC_OPEN_SUB1 = 26, 171
if open_sub2 == DOC_OPEN_SUB2 and open_sub1 == DOC_OPEN_SUB1:
    ok("SURVIVING-BRANCHES",
       "surviving branches match: sub2=%d, sub1=%d (DAG == CURRENT_STATUS)"
       % (DOC_OPEN_SUB2, DOC_OPEN_SUB1))
elif open_sub2 <= DOC_OPEN_SUB2 and open_sub1 <= DOC_OPEN_SUB1:
    # the DAG closed MORE than the docs record: report loudly, do not fail
    gap("SURVIVING-BRANCHES",
        "the DAG is STRONGER than the docs here, which is not an inconsistency "
        "but IS staleness: DAG open branches sub2=%s sub1=%s vs CURRENT_STATUS "
        "%d/%d.%s Update CURRENT_STATUS.md / FRONTIER_V2.md; nothing in the DAG "
        "needs changing."
        % (open_sub2, open_sub1, DOC_OPEN_SUB2, DOC_OPEN_SUB1,
           "  Both f31 windows are now CLOSED (0 open branches) -- see the "
           "column_lemmas block in proof_dag.json counts."
           if open_sub2 == 0 and open_sub1 == 0 else ""))
else:
    inconsistency("SURVIVING-BRANCHES",
        "surviving-branch mismatch in the DANGEROUS direction: the DAG has MORE "
        "open branches than CURRENT_STATUS records -- DAG sub2=%s sub1=%s vs "
        "CURRENT_STATUS %d/%d. The doc claims more closure than the DAG supports."
        % (open_sub2, open_sub1, DOC_OPEN_SUB2, DOC_OPEN_SUB1))
# note the alt-regime 27 that the DAG models as 15 defect-0 families
if re.search(r"\b27\b", FR) or re.search(r"alt.*?27", CS):
    _alt = NODES.get("subcase:sub1_alt", {})
    if _alt.get("branches_open") == 0 and _alt.get("closed"):
        gap("ALT-REGIME-27",
            "CURRENT_STATUS/FRONTIER track a 27-branch alternate-regime (a11-15, "
            "v<0) sweep as OPEN. It is CLOSED: the DAG registers all %s L_alt "
            "branch keys and closes %s of them by %s. The 27 is the post-C33/C34 "
            "residual, a historical figure -- the docs are STALE, not wrong, and "
            "the direction is safe (the DAG is stronger). GAP-ALT-STATES (39 "
            "modelled states vs %s surviving) is RETIRED: those states all sit at "
            "a_t >= 11 and the bound empties every such branch, so they never "
            "needed modelling. Pinned in c0_partition.py under RETIRED."
            % (_alt.get("n_branches"), _alt.get("branches_closed"),
               _alt.get("closure_mechanism"), _alt.get("states_surviving")))
    else:
        gap("ALT-REGIME-27",
            "CURRENT_STATUS/FRONTIER track a 27-branch alternate-regime (a11-15, "
            "v<0) sweep. The DAG REGISTERS all %s L_alt branch keys (%s closed "
            "whole by C33+C34, %s open), but it COVERS far less: %s of the open "
            "branches carry only a forced-defect-0 state overlay and %s carry no "
            "state model at all, %s modelled states against %s surviving. "
            "Registration is not coverage -- scope gap, not a contradiction. "
            "Pinned by exact key in c0_partition.py (GAP-ALT-STATES)."
            % (_alt.get("n_branches"), _alt.get("branches_killed_whole"),
               _alt.get("branches_open"),
               _alt.get("surviving_C33_C34_with_state_overlay"),
               _alt.get("surviving_C33_C34_unmodelled"),
               _alt.get("states_modelled"), _alt.get("states_surviving")))

# --- I4: FRONTIER CLOSED defect-0 families vs DAG ---------------------------
dag_closed_fam = sorted(n["family"] for n in DAG["nodes"]
                        if n["type"] == "branch" and n.get("window") == "altdefect0"
                        and n["closed"])
fr_closed = sorted(set(re.findall(r"`(a1[0-9]_b\d+_T\d)`\s*--\s*\d+/\d+ CLOSED", FR)))
if fr_closed and set(fr_closed) - set(dag_closed_fam):
    inconsistency("DEFECT0-CLOSED",
        "FRONTIER lists CLOSED defect-0 families %s not closed in the DAG (%s)"
        % (fr_closed, dag_closed_fam))
else:
    ok("DEFECT0-CLOSED",
       "CLOSED defect-0 families agree: DAG=%s" % dag_closed_fam)

# --- I5: CERTIFICATE-FOUND kill with no state in the coverage graph ---------
found_unresolved = [u for u in DAG["unmapped"]
                    if u.get("kind") == "certificate-unresolved"
                    and u.get("status") == "CERTIFICATE-FOUND"]
if found_unresolved:
    inconsistency("ORPHAN-CERTIFICATE",
        "%d object certificate(s) are CERTIFICATE-FOUND yet map to NO state in the "
        "ledger/universe join (%s). A certified kill with no obligation node is a "
        "coverage gap: the certificate corpus is stronger than the ledger records."
        % (len(found_unresolved), ", ".join(u["kill_id"] for u in found_unresolved)))
else:
    ok("ORPHAN-CERTIFICATE", "every CERTIFICATE-FOUND maps to a graph node")

# --- report ----------------------------------------------------------------
out("")
out("-- CONSISTENCY / INCONSISTENCY FINDINGS --")
n_incon = 0
for sev, code, msg in FINDINGS:
    if sev == "OK":
        out("  [OK]            %-22s %s" % (code, msg))
for sev, code, msg in FINDINGS:
    if sev == "COVERAGE-GAP":
        # gaps are surfaced loudly but do not fail the build
        print("  [COVERAGE-GAP]  %-22s %s" % (code, msg))
for sev, code, msg in FINDINGS:
    if sev == "INCONSISTENCY":
        n_incon += 1
        print("  [INCONSISTENCY] %-22s %s" % (code, msg))

out("")
out("UNMAPPED bucket: %d" % len(DAG["unmapped"]))
uk = Counter(u.get("kind") for u in DAG["unmapped"])
for k, v in sorted(uk.items()):
    out("   %-32s %d" % (k, v))

print("")
print("VERDICT: %d inconsistency finding(s)%s"
      % (n_incon, " -- DAG is consistent with the docs" if n_incon == 0 else ""))
sys.exit(0 if n_incon == 0 else 1)
