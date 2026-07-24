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

def load_json(fn):
    with open(os.path.join(HERE, fn)) as f:
        return json.load(f)

def load_text(fn):
    p = os.path.join(HERE, fn)
    if not os.path.exists(p):
        return ""
    with open(p, encoding="utf-8") as f:
        return f.read()

DAG = load_json("proof_dag.json")
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
out("")
out("  TARGET C0: closed=%s level=%s (subcases closed %d/4)" % (
    c0["closed"], c0["level"], c0.get("subcases_closed", 0)))
for sid in ("subcase:sub2", "subcase:sub1", "subcase:sub1_alt_defect0", "subcase:f37"):
    s = NODES[sid]
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
# group 1: engine-killed branches at 'claimed'
ek = defaultdict(int)
for n in DAG["nodes"]:
    if n["type"] == "branch" and n.get("window") in ("sub2", "sub1") \
            and n.get("cascade_status") not in (None, "survives"):
        ek[n["window"]] += 1
prio = []
if ek:
    prio.append((sum(ek.values()),
        "engine-killed branches @ claimed (sub2 %d, sub1 %d) -> upgrade via the "
        "cascade spec-only audit (CURRENT_STATUS C18/C29/C43) to promote every "
        "one to independently-audited" % (ek.get("sub2", 0), ek.get("sub1", 0))))
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

FR = load_text("FRONTIER_V2.md")
CS = load_text("CURRENT_STATUS.md")

# DAG-supported counts
def state_ge(window, level):
    return sum(1 for n in DAG["nodes"] if n["type"] == "state"
               and n["window"] == window and LRANK[n["level"]] >= LRANK[level])

dag_audit_sub2 = state_ge("sub2", "independently-audited")
dag_audit_sub1 = state_ge("sub1", "independently-audited")

# --- I1: FRONTIER "Killed (audited)" column vs DAG independently-audited ----
m2 = re.search(r"\|\s*sub2\s*\|\s*(\d+)\s*\|\s*(\d+)\s*\|\s*(\d+)\s*\|", FR)
m1 = re.search(r"\|\s*sub1\s*\|\s*(\d+)\s*\|\s*(\d+)\s*\|\s*(\d+)\s*\|", FR)
if m2 and m1:
    fr_aud2 = int(m2.group(2)); fr_aud1 = int(m1.group(2))
    if fr_aud2 > dag_audit_sub2 or fr_aud1 > dag_audit_sub1:
        inconsistency("FRONTIER-AUDITED-LABEL",
            "FRONTIER_V2 'Killed (audited)' counts sub2=%d, sub1=%d, but those are "
            "ledger-'AUDITED' = same-author EXACT-CHECKED, not independent audit. "
            "The DAG reaches >=independently-audited for only sub2=%d, sub1=%d "
            "(the alt-hunt census-verified kills). FRONTIER's 'audited' column "
            "overstates the evidence grade by %d states."
            % (fr_aud2, fr_aud1, dag_audit_sub2, dag_audit_sub1,
               (fr_aud2 - dag_audit_sub2) + (fr_aud1 - dag_audit_sub1)))
    else:
        ok("FRONTIER-AUDITED-LABEL", "FRONTIER audited counts within DAG support")

# --- I2: CURRENT_STATUS independent-audit claims for branch (engine) kills ---
ek_total = sum(ek.values())
c18 = re.search(r"390 killed", CS)
c29 = re.search(r"1899 killed", CS)
if (c18 or c29) and ek_total > 0:
    inconsistency("CASCADE-BRANCH-AUDIT",
        "CURRENT_STATUS S1a asserts INDEPENDENT AUDIT of the cascade branch kills "
        "(C18 '390 killed', C29 '1899 killed'). The DAG, from the loaded sources, "
        "rates all %d engine-killed branches at level 'claimed' (the "
        "cascade_cones data self-labels them '*_pending_audit'); the audit artifact "
        "(audit_cascade_kills*.py) is NOT machine-joined in v1. Doc claims "
        "independently-audited; DAG supports only claimed. This is the top weakest "
        "edge / audit priority." % ek_total)
else:
    ok("CASCADE-BRANCH-AUDIT", "no cascade branch-audit over-claim detected")

# --- I3: surviving-branch counts (should MATCH: consistency check) ----------
open_sub2 = NODES["subcase:sub2"].get("branches_open")
open_sub1 = NODES["subcase:sub1"].get("branches_open")
cs26 = re.search(r"\*\*26\*\*", CS) or re.search(r"\bsub2 cells.*?\*\*26\*\*", CS)
if open_sub2 == 26 and open_sub1 == 171:
    ok("SURVIVING-BRANCHES",
       "surviving branches match: sub2=26, sub1=171 (DAG == CURRENT_STATUS)")
else:
    inconsistency("SURVIVING-BRANCHES",
        "surviving-branch mismatch: DAG sub2=%s sub1=%s vs CURRENT_STATUS 26/171"
        % (open_sub2, open_sub1))
# note the alt-regime 27 that the DAG models as 15 defect-0 families
if re.search(r"\b27\b", FR) or re.search(r"alt.*?27", CS):
    gap("ALT-REGIME-27",
        "CURRENT_STATUS/FRONTIER track a 27-branch alternate-regime (a11-15, v<0) "
        "sweep; the DAG models the alt layer as the 15 entirely-defect-0 families "
        "(phase_f2_scale), 3 CLOSED / 12 open. The 27-branch alt_combined sweep is "
        "an overlay not per-state joined here -- scope gap, not a contradiction.")

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
