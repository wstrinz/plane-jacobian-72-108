#!/usr/bin/env python3
"""
proof_dag.py  ->  proof_dag.json

The COVERAGE PROOF-DAG for the (72,108) program.

Purpose
-------
External review named the project's central missing artifact: a *machine-enforced*
object

    certificate -> state -> cell -> branch -> subcase -> target(C0)

in which "this branch is closed" is a COMPUTED fact, never a hand-maintained
roll-up.  This tool builds that DAG from the read-only kill/enumeration artifacts,
computes a per-node evidence level as the MINIMUM along each node's *required*
support (with a disjunctive exception at the state layer -- see below), and emits
`proof_dag.json` deterministically (sorted keys, no timestamps).

Node types
----------
    target      C0 (single root)
    subcase     coarse regime/window partition of C0 (sub2, sub1, sub1 alt
                defect-0, f37).  A `corner` bucket holds auxiliary lemma kills.
    branch      a cascade-cone branch (win, a_t, b, branch_label) OR an alt
                defect-0 family (bid).  Enumerated by cascade_cones_*.json.
    cell        a phase-D flag-case (d2_zero, sigma_zero, g_zero_levels) inside a
                surviving branch.  Enumerated by phase_d_states_*.json.
    state       a residual degree-state (the canonical kill key).  Only KILLED
                states are materialized as nodes; surviving states are carried as
                an aggregate `surviving` count on their cell (closure requires
                surviving==0).  This keeps the DAG ~4.5k nodes instead of ~52k
                while leaving closure a computed fact.
    certificate an object-level kill certificate (kill_certificates/*.json); a
                leaf that, when CERTIFICATE-FOUND, promotes its state to certified.

Evidence levels (ascending)
----------------------------
    open  <  claimed  <  exact-checked  <  independently-audited  <  certified

Aggregation
    * state  : DISJUNCTIVE.  A state needs only one valid kill, so its level is
               the MAX over its independent kill mechanisms / audits / certs.
    * cell / branch / subcase / target : CONJUNCTIVE.  Level = MIN over all
               required children AND the exhaustiveness edge.  A structural node
               is "closed" only if every required child is closed; its closed
               level is gated by its exhaustiveness edge (a branch closes at level
               X only if all states dead at >=X AND exhaustiveness >= exact-checked
               -- implemented by folding the exhaustiveness edge level into the
               MIN, judgment-only exhaustiveness being level 'claimed').

READ-ONLY on every committed artifact.  Writes only proof_dag.json.
"""
import json, os, re, sys, hashlib
from collections import defaultdict, Counter

HERE = os.path.dirname(os.path.abspath(__file__))
LOADED = []

def load(fn, required=True):
    p = os.path.join(HERE, fn)
    if not os.path.exists(p):
        if required:
            sys.stderr.write("FATAL: missing %s\n" % fn); sys.exit(1)
        return None
    LOADED.append(fn)
    with open(p) as f:
        return json.load(f)

def load_text(fn):
    p = os.path.join(HERE, fn)
    if not os.path.exists(p):
        return None
    LOADED.append(fn)
    with open(p, encoding="utf-8") as f:
        return f.read()

# ---------------------------------------------------------------------------
# Levels
# ---------------------------------------------------------------------------
LEVELS = ["open", "claimed", "exact-checked", "independently-audited", "certified"]
LRANK = {l: i for i, l in enumerate(LEVELS)}
def lmin(*xs): return min(xs, key=lambda l: LRANK[l])
def lmax(*xs): return max(xs, key=lambda l: LRANK[l])

# ledger attribution audit-field -> base evidence level.  The ledger's own
# "AUDITED" is the producing engine's exact self-check (same author, NOT an
# independent audit); it is deliberately mapped to exact-checked, not
# independently-audited.  Independent-audit / certified upgrades come only from
# the alt-hunt census and the object certificates.
AUDIT_LEVEL = {
    "AUDITED": "exact-checked",
    "TRANSFERRED-AUDITED": "exact-checked",
    "PENDING": "claimed",
    "PENDING-AMBIGUOUS-MAP": "claimed",
}

# ---------------------------------------------------------------------------
# Load sources
# ---------------------------------------------------------------------------
LED   = load("state_kill_ledger.json")
CC2   = load("cascade_cones_qt_inf_rl.json")
CC1   = load("cascade_cones_sub1_qt_inf_rl.json")
D2    = load("phase_d_states_sub2.json")
D1    = load("phase_d_states_sub1.json")
SCALE = load("phase_f2_scale.json")
F2S2  = load("phase_f2_sub2.json")
AH    = load("alt_hunt_results.json")
AUC   = load("audit_alt_hunt_census.json")

CERT_DIR = os.path.join(HERE, "kill_certificates")
CERTS = []
for fn in sorted(os.listdir(CERT_DIR)):
    if fn == "status_log.json" or not fn.endswith(".json"):
        continue
    with open(os.path.join(CERT_DIR, fn)) as f:
        CERTS.append(json.load(f))
    LOADED.append("kill_certificates/" + fn)
load("kill_certificates/status_log.json", required=False)

# ---------------------------------------------------------------------------
# phase-D universe indices (canonical raw-state key -> exists) and case index
# ---------------------------------------------------------------------------
def bstr(b): return "".join(str(x) for x in b)

UNIV = {"sub2": {}, "sub1": {}}          # canonical str -> flag-case id
CASE_STATES = {"sub2": {}, "sub1": {}}   # flag-case id -> [canonical str,...]
CASE_META = {}                            # flag-case id -> case dict meta
CELL_OF_CASE = {}                         # flag-case id -> branch id
CASE_BY_CELLID = {"sub2": {}, "sub1": {}}  # phase_f2 cellid -> flag-case id

def canon_str(win, c, s):
    return "|".join(str(x) for x in
        (win, c["a_t"], tuple(c["b"]), c["branch"], c["d2_zero"], c["sigma_zero"],
         tuple(c["g_zero_levels"]), s["deg_d1"], s["deg_d2"], s["deg_e"], s["deg_sigma"]))

def flagcase_id(win, c):
    gz = "".join(str(x) for x in c["g_zero_levels"]) or "-"
    return "cell:%s:a%d_b%s_%s:dz%d_sz%d_gz%s" % (
        win, c["a_t"], bstr(c["b"]), c["branch"],
        int(c["d2_zero"]), int(c["sigma_zero"]), gz)

def branch_id(win, a, b, br):
    return "branch:%s:a%d_b%s_%s" % (win, a, bstr(b), br)

def legacy_cellid(win, c):
    sz = "sz1" if c["sigma_zero"] else "sz0"
    dz = "dz1" if c["d2_zero"] else "dz0"
    gz = "gz" + ("-" if not c["g_zero_levels"] else "".join(map(str, c["g_zero_levels"])))
    return "%s:a%d_b%s_%s_%s_%s_%s" % (win, c["a_t"], bstr(c["b"]),
                                       c["branch"], sz, dz, gz)

for win, D in (("sub2", D2), ("sub1", D1)):
    for c in D["cases"]:
        fid = flagcase_id(win, c)
        bid = branch_id(win, c["a_t"], c["b"], c["branch"])
        CASE_META[fid] = dict(win=win, a_t=c["a_t"], b=list(c["b"]),
                              branch=c["branch"], d2_zero=c["d2_zero"],
                              sigma_zero=c["sigma_zero"],
                              g_zero_levels=list(c["g_zero_levels"]),
                              state_count=c["state_count"])
        CELL_OF_CASE[fid] = bid
        CASE_BY_CELLID[win][legacy_cellid(win, c)] = fid
        lst = CASE_STATES[win].setdefault(fid, [])
        for s in c["states"]:
            ck = canon_str(win, c, s)
            UNIV[win][ck] = fid
            lst.append(ck)

# ---------------------------------------------------------------------------
# STATE evidence accumulator
#   canonical str -> dict(window, flagcase, evidence[list], mechanisms[set],
#                         field_scope, notes[set], level)
# ---------------------------------------------------------------------------
STATE = {}          # canonical str -> record
UNMAPPED = []       # loud bucket, never silently dropped
CORNER = {}         # corner label -> record
DEFECT0 = {}        # (bid, dd2) -> record

def state_rec(win, ck, flagcase):
    r = STATE.get(ck)
    if r is None:
        r = dict(window=win, canonical_key=ck, flagcase=flagcase,
                 evidence=[], mechanisms=set(), field_scopes=set(), notes=set())
        STATE[ck] = r
    return r

def add_state_evidence(win, ck, level, mechanism, source, field_scope, note=None):
    fid = UNIV[win].get(ck)
    if fid is None:
        UNMAPPED.append(dict(kind="state-kill-off-universe", window=win,
                             canonical_key=ck, source=source, mechanism=mechanism,
                             reason="canonical key not present in phase-D universe"))
        return False
    r = state_rec(win, ck, fid)
    r["evidence"].append(dict(level=level, mechanism=mechanism, source=source,
                              field_scope=field_scope))
    r["mechanisms"].add(mechanism)
    r["field_scopes"].add(field_scope)
    if note:
        r["notes"].add(note)
    return True

# --- field scope from mechanism string -------------------------------------
def field_scope_of(mech):
    m = mech.lower()
    if "mod p" in m:
        return "field-split: characteristic p (mod-p certificate)"
    if "saturat" in m or "depth-cap" in m or "depth-8" in m:
        return "char-0 / Q-algebra (saturation/denominator-cleared)"
    return "char-0 / Q-algebra (all places)"

# ===========================================================================
# 1. Ingest the unified ledger kills (already normalized to canonical keys)
# ===========================================================================
n_ledger_phase_d = 0
for r in LED["kills"]:
    w = r["window"]
    if w in ("sub2", "sub1"):
        ck = r["canonical_key"]
        for a in r["attributions"]:
            lvl = AUDIT_LEVEL.get(a["audit"], "claimed")
            fs = field_scope_of(a["mechanism"])
            note = None
            if a["audit"] == "PENDING-AMBIGUOUS-MAP":
                note = "ambiguous secondary-source map (recorded, not a clean kill)"
            if add_state_evidence(w, ck, lvl, a["mechanism"], a["source"], fs, note):
                n_ledger_phase_d += 1
    elif w == "altdefect0":
        key = (r["bid"], r["deg_d2"])
        rec = DEFECT0.setdefault(key, dict(bid=r["bid"], deg_d2=r["deg_d2"],
                                           evidence=[], mechanisms=set(),
                                           field_scopes=set(), notes=set()))
        for a in r["attributions"]:
            lvl = AUDIT_LEVEL.get(a["audit"], "claimed")
            fs = field_scope_of(a["mechanism"])
            rec["evidence"].append(dict(level=lvl, mechanism=a["mechanism"],
                                        source=a["source"], field_scope=fs))
            rec["mechanisms"].add(a["mechanism"]); rec["field_scopes"].add(fs)
    elif w == "corner":
        lab = r.get("label") or r["canonical_key"]
        rec = CORNER.setdefault(lab, dict(label=lab, evidence=[], mechanisms=set(),
                                          field_scopes=set(), notes=set()))
        for a in r["attributions"]:
            lvl = AUDIT_LEVEL.get(a["audit"], "claimed")
            fs = field_scope_of(a["mechanism"])
            rec["evidence"].append(dict(level=lvl, mechanism=a["mechanism"],
                                        source=a["source"], field_scope=fs))
            rec["mechanisms"].add(a["mechanism"]); rec["field_scopes"].add(fs)

# surface the ledger's own ambiguous mappings into the loud unmapped bucket
for amb in LED.get("ambiguous", []):
    UNMAPPED.append(dict(kind="ledger-ambiguous-map",
                         source=amb.get("source"), name=amb.get("name"),
                         ncand=amb.get("ncand"),
                         reason=amb.get("reason", "ambiguous secondary-source map")))

# ===========================================================================
# 2. Ingest ALT-HUNT kills (NOT in the ledger) + attach the independent audit
#    alt_hunt state key:  <win>:aA_bBBBB_TX_szS_dzD_gzG#stateN
#    degs order = [deg_d1, deg_d2, deg_sigma, deg_e]  ->  canonical (d1,d2,e,sig)
# ===========================================================================
CELLID_RE = re.compile(r"(sub1|sub2):a(\d+)_b(\d+)_(T\d)_sz(\d)_dz(\d)_gz(-|\d+)")

# audit category by kill_id (independent spec-only auditor)
AUDIT_CAT = {res["kill_id"]: res["category"] for res in AUC["results"]}
AUDIT_LEVEL_BY_CAT = {
    "FULLY-VERIFIED": "independently-audited",
    "VERIFIED-DATA-ONLY": "exact-checked",
    "DISAGREEMENT": "open",
    "UNPARSEABLE": "claimed",
}

def alt_hunt_canonical(key, degs):
    m = CELLID_RE.match(key.split("#")[0])
    if not m:
        return None
    win = m.group(1); a = int(m.group(2)); b = tuple(int(x) for x in m.group(3))
    br = m.group(4); sz = bool(int(m.group(5))); dz = bool(int(m.group(6)))
    gz = [] if m.group(7) == "-" else [int(x) for x in m.group(7)]
    d1, sig, d2, e = degs                         # alt-hunt order [d1,sigma,d2,e]
    def v(x): return "-inf" if x in ("-inf", None) else int(x)
    ck = "|".join(str(x) for x in
        (win, a, b, br, dz, sz, tuple(gz), v(d1), v(d2), v(e), v(sig)))
    return win, ck

n_alt_hunt = 0
for s in AH["states"]:
    if s["verdict"] != "KILLED":
        continue
    parsed = alt_hunt_canonical(s["key"], s["degs"])
    if parsed is None:
        UNMAPPED.append(dict(kind="alt-hunt-unparseable", key=s["key"],
                             reason="cellid did not parse"))
        continue
    win, ck = parsed
    cat = AUDIT_CAT.get(s["key"], None)
    lvl = AUDIT_LEVEL_BY_CAT.get(cat, "claimed")
    note = "alt-hunt depth-2 residue kill; independent audit=%s" % (cat or "NONE")
    mech = "alt-hunt residue/Groebner (depth-2)"
    if add_state_evidence(win, ck, lvl, mech, "alt_hunt_results.json",
                          "char-0 / Q-algebra (all places)", note):
        n_alt_hunt += 1
    else:
        # add_state_evidence already logged it as off-universe; add richer reason
        UNMAPPED[-1]["reason"] = ("alt-hunt KILLED state did not join a phase-D "
                                  "canonical key (flag-case/deg signature absent)")

# j6_msolve audited kills that are not carried in alt_hunt_results.json
j6_ids = [res["kill_id"] for res in AUC["results"] if res.get("source") == "j6_msolve"]
for kid in j6_ids:
    UNMAPPED.append(dict(kind="j6-msolve-audited-nostate", key=kid,
                         audit=AUDIT_CAT.get(kid),
                         reason="j6_msolve kill independently audited "
                                "(census) but no degs in loaded sources to join "
                                "to a phase-D canonical key"))

# ===========================================================================
# 3. Ingest CERTIFICATES and resolve each to its state / family / corner.
#    A CERTIFICATE-FOUND cert promotes its target to 'certified'.
# ===========================================================================
# name -> canonical str, harvested from ledger attribution detail.name
NAME_TO_CK = {}
for r in LED["kills"]:
    if r["window"] not in ("sub2", "sub1"):
        continue
    for a in r["attributions"]:
        d = a.get("detail") or {}
        nm = d.get("name")
        if nm:
            NAME_TO_CK.setdefault(nm, set()).add(r["canonical_key"])

# phase_f2_sub2 key -> canonical str (via cellid + degs [d1,sigma,d2,e])
F2S2_KEY_TO_CK = {}
for st in F2S2["states"]:
    fid = CASE_BY_CELLID["sub2"].get(st["cellid"])
    if fid is None:
        continue
    d1, sg, d2, e = st["degs"]
    # locate matching canonical in that flag-case
    for ck in CASE_STATES["sub2"].get(fid, []):
        parts = ck.split("|")
        # parts[-4:] = deg_d1, deg_d2, deg_e, deg_sigma
        if (parts[-4] == str(d1) and parts[-3] == str(d2)
                and parts[-2] == str(e) and parts[-1] == str(sg)):
            F2S2_KEY_TO_CK[st["key"]] = ck
            break

# case-label -> canonical str, harvested from ledger msolve/blowup attributions
CASE_TO_CK = {}
for r in LED["kills"]:
    if r["window"] not in ("sub2", "sub1"):
        continue
    for a in r["attributions"]:
        d = a.get("detail") or {}
        c = d.get("case")
        if c and a["source"] in ("msolve_bridge_results.json",
                                 "blowup_sweep_results.json"):
            CASE_TO_CK.setdefault(c, set()).add(r["canonical_key"])

NOTE_SUB2 = re.compile(r"(sub2:\S+?)#state\d+ degs=\[([^\]]+)\]")

def cert_targets(cert):
    """Return list of ('state', canonical_str) / ('defect0', (bid,dd2)) /
    ('corner', label) targets, plus a human reason if unresolved."""
    kid = cert["kill_id"]
    cat = cert.get("category")
    body = kid.split(":", 1)[1] if ":" in kid else kid
    # phase_f2_sub2:<cellid-suffix>sN
    if cat == "phase_f2_sub2":
        full = "sub2:" + body                # -> sub2:...gz-s0  == key ...#stateN
        m = re.match(r"(sub2:.+_gz(?:-|\d+))s(\d+)$", full)
        if m:
            statekey = "%s#state%s" % (m.group(1), m.group(2))
            ck = F2S2_KEY_TO_CK.get(statekey)
            if ck:
                return [("state", ck)], None
        return [], "phase_f2_sub2 cert cellid/state-index unresolved"
    # harvest / bridge labelled kills: name matches a ledger detail.name
    if cat == "harvest":
        cks = NAME_TO_CK.get(body)
        if cks:
            return [("state", ck) for ck in sorted(cks)], None
        # alt defect-0 sup-index harvest (a<K>_b<...>_T<n>_sup<idx>)
        m = re.match(r"(a\d+_b\d+_T\d)_sup(\d+)$", body)
        if m:
            SUP = {"12": 5, "14": 6}          # per D2_THRESHOLD ladder
            dd2 = SUP.get(m.group(2))
            if dd2 is not None and (m.group(1), dd2) in DEFECT0:
                return [("defect0", (m.group(1), dd2))], None
            return [], "harvest sup-index not on tracked defect-0 deg_d2 ladder"
        return [], "harvest cert name did not match a ledger kill"
    # d2_threshold: a<..>_T<n>_d<dd2>
    if cat == "d2_threshold":
        m = re.match(r"(a\d+_b\d+_T\d)_d(\d+)$", body)
        if m and (m.group(1), int(m.group(2))) in DEFECT0:
            return [("defect0", (m.group(1), int(m.group(2))))], None
        return [], "d2_threshold cert family/deg_d2 not in defect-0 census"
    # msolve_blowup: sub2_sNN (phase-D) or a<..>_dK / a<..>_NN (defect-0)
    if cat == "msolve_blowup":
        # first try the ledger's own case->canonical map (sub2_sNN etc.)
        cks = CASE_TO_CK.get(body)
        if cks:
            return [("state", ck) for ck in sorted(cks)], None
        note = cert.get("reason", "") + " " + json.dumps(cert.get("manifest_recipe", {}))
        m = NOTE_SUB2.search(note)
        if m:
            degs = [None if t.strip() in ("-inf", "None") else int(t)
                    for t in m.group(2).split(",")]
            # note degs order [d1,d2,sig,e]
            fid = CASE_BY_CELLID["sub2"].get(m.group(1))
            if fid:
                d1, d2, sg, e = degs
                for ck in CASE_STATES["sub2"].get(fid, []):
                    p = ck.split("|")
                    def z(x): return "-inf" if x is None else str(x)
                    if (p[-4] == z(d1) and p[-3] == z(d2) and p[-2] == z(e)
                            and p[-1] == z(sg)):
                        return [("state", ck)], None
        md = re.match(r"(a\d+_b\d+_T\d)_d(\d+)$", body)
        if md and (md.group(1), int(md.group(2))) in DEFECT0:
            return [("defect0", (md.group(1), int(md.group(2))))], None
        mi = re.match(r"(a\d+_b\d+_T\d)_(\d+)$", body)
        if mi and (mi.group(1), 6) in DEFECT0:      # sup17 = deg_d2=6 last survivor
            return [("defect0", (mi.group(1), 6))], None
        return [], "msolve_blowup cert target unresolved from loaded sources"
    return [], "unknown certificate category"

CERT_NODES = []
n_cert_found = n_cert_resolved = 0
for cert in CERTS:
    kid = cert["kill_id"]
    status = cert["status"]
    found = (status == "CERTIFICATE-FOUND")
    n_cert_found += int(found)
    targets, reason = cert_targets(cert)
    tgt_ids = []
    for kind, val in targets:
        if kind == "state":
            tgt_ids.append("state:" + val)
        elif kind == "defect0":
            tgt_ids.append("state:altdefect0|%s|%s" % (val[0], val[1]))
        elif kind == "corner":
            tgt_ids.append("state:corner|%s" % val)
    node = dict(id="cert:" + kid, type="certificate", kill_id=kid,
                category=cert.get("category"), status=status,
                found=found, targets=tgt_ids,
                lift_method=cert.get("lift_method"),
                generator_sha256=cert.get("generator_sha256"),
                level="certified" if found else "claimed")
    if reason and not targets:
        node["unresolved_reason"] = reason
        UNMAPPED.append(dict(kind="certificate-unresolved", kill_id=kid,
                             status=status, reason=reason))
    CERT_NODES.append(node)
    if targets:
        n_cert_resolved += 1
        # promote target state(s) to certified when a certificate was found
        for kind, val in targets:
            if kind == "state" and found:
                r = STATE.get(val)
                if r is not None:
                    r["evidence"].append(dict(level="certified",
                        mechanism="object certificate (lift, kernel-checkable)",
                        source="kill_certificates/" + kid, field_scope=field_scope_of("")))
                    r["notes"].add("object certificate CERTIFICATE-FOUND")
            elif kind == "defect0" and found and val in DEFECT0:
                DEFECT0[val]["evidence"].append(dict(level="certified",
                    mechanism="object certificate (lift, kernel-checkable)",
                    source="kill_certificates/" + kid, field_scope=field_scope_of("")))
                DEFECT0[val]["notes"].add("object certificate CERTIFICATE-FOUND")

# ---------------------------------------------------------------------------
# finalize state levels (DISJUNCTIVE: best evidence closes a state)
# ---------------------------------------------------------------------------
for ck, r in STATE.items():
    r["level"] = "open"
    for e in r["evidence"]:
        r["level"] = lmax(r["level"], e["level"])
for key, r in DEFECT0.items():
    r["level"] = "open"
    for e in r["evidence"]:
        r["level"] = lmax(r["level"], e["level"])
for lab, r in CORNER.items():
    r["level"] = "open"
    for e in r["evidence"]:
        r["level"] = lmax(r["level"], e["level"])

# ===========================================================================
# 4. Build the node/edge graph bottom-up
# ===========================================================================
NODES = {}
EDGES = []

def add_node(node): NODES[node["id"]] = node
def add_edge(parent, child, predicate, exhaustiveness_ref, exhaustiveness_level,
             machine_checkable, field_scope, notes=None):
    EDGES.append(dict(parent=parent, child=child, predicate=predicate,
                      exhaustiveness_ref=exhaustiveness_ref,
                      exhaustiveness_level=exhaustiveness_level,
                      machine_checkable=machine_checkable,
                      field_scope=field_scope, notes=notes or ""))

# ---- target ---------------------------------------------------------------
add_node(dict(id="C0", type="target", label="C0: no [P,Q]=x^2 in the "
              "Prop-4.3 case-(8,28) subcases (1)&(2)", level="open", closed=False))

# ---- subcases -------------------------------------------------------------
SUBCASES = {
    "subcase:sub2": dict(label="f31 window sub2 (q+t+inf, depth 4)", window="sub2"),
    "subcase:sub1": dict(label="f31 window sub1 standard (q+t+inf, depth 4)", window="sub1"),
    "subcase:sub1_alt_defect0": dict(label="sub1 alternate regime, entirely-"
                                     "defect-0 families (a>=11)", window="altdefect0"),
    "subcase:f37": dict(label="f37 branch (closed, C11)", window="f37"),
    "subcase:corner": dict(label="auxiliary corner/lemma kills (not a C0 partition)",
                           window="corner"),
}
for sid, meta in SUBCASES.items():
    add_node(dict(id=sid, type="subcase", level="open", closed=False, **meta))

# subcase -> C0 exhaustiveness is a JUDGMENT reference (GGHV22 Prop 4.3 + the
# field-split framework C14-16 + alt-regime C44), not machine-checkable in v1.
for sid in ("subcase:sub2", "subcase:sub1", "subcase:sub1_alt_defect0", "subcase:f37"):
    add_edge("C0", sid,
             predicate="C0 holds iff every regime/window subcase is closed",
             exhaustiveness_ref="JUDGMENT: GGHV22 Prop 4.3 case-(8,28) subcases "
                                "(1)&(2); field-split framework CURRENT_STATUS "
                                "C14-C16; alternate-regime partition C44",
             exhaustiveness_level="claimed", machine_checkable=False,
             field_scope="char-0 / Q-algebra; f37 restricted to char != 3,5",
             notes="v1: subcase->C0 completeness is judgment-referenced, not "
                   "machine-enforced; this gates C0 closure at <= claimed.")

# f37 subcase: closed by C11 (judgment/cert reference; no state data here)
f37 = NODES["subcase:f37"]
f37["closed"] = True
f37["level"] = "exact-checked"
f37["field_scope"] = "char != 3,5 (integer cert D=46875=3*5^6; F37_SATURATION_REPORT)"
f37["closure_note"] = ("closed by C11 (f37_sat_verify.py, same-author exact "
                       "checker + Lean-kernel integer certificate); represented "
                       "as a judgment-referenced leaf, not recomputed in this DAG.")

# ---- branches from the cascade enumeration --------------------------------
def cascade_summary_ok(cc):
    s = cc["summary"]
    got = s["surviving_branches"] + s["engine_killed_pending_audit"]
    return got == s["open_branches_processed"], got, s["open_branches_processed"]

BRANCH_STATUS = {}   # branch id -> cascade status
for win, cc, sub in (("sub2", CC2, "subcase:sub2"), ("sub1", CC1, "subcase:sub1")):
    ok, got, exp = cascade_summary_ok(cc)
    exh_ref = ("MACHINE: cascade_cones summary survivor+killed==open_branches_"
               "processed (%d==%d) [%s]" % (got, exp, "PASS" if ok else "FAIL"))
    for b in cc["branches"]:
        bid = branch_id(win, b["a_t"], b["b"], b["branch"])
        BRANCH_STATUS[bid] = b["status"]
        killed = (b["status"] != "survives")
        node = dict(id=bid, type="branch", window=win, a_t=b["a_t"],
                    b=list(b["b"]), branch=b["branch"], cascade_status=b["status"],
                    survivor_case_count=b.get("survivor_case_count", 0),
                    level="open", closed=False)
        add_node(node)
        add_edge(sub, bid,
                 predicate="branch (a_t=%d, b=%s, %s) in the depth-4 %s cascade cone"
                           % (b["a_t"], bstr(b["b"]), b["branch"], win),
                 exhaustiveness_ref=exh_ref,
                 exhaustiveness_level="exact-checked" if ok else "claimed",
                 machine_checkable=True,
                 field_scope="char-0 / Q-algebra (all places)",
                 notes="cascade branch-enumeration completeness rests on the "
                       "cone lemmas (CASCADE_CONE_LEMMAS*.md).")

# ---- cells (phase-D flag-cases) inside surviving branches -----------------
for fid, meta in CASE_META.items():
    win = meta["win"]
    bid = CELL_OF_CASE[fid]
    add_node(dict(id=fid, type="cell", window=win, branch=bid,
                  a_t=meta["a_t"], b=meta["b"], branch_label=meta["branch"],
                  d2_zero=meta["d2_zero"], sigma_zero=meta["sigma_zero"],
                  g_zero_levels=meta["g_zero_levels"],
                  state_total=meta["state_count"], level="open", closed=False))
    add_edge(bid, fid,
             predicate="flag-case d2_zero=%s sigma_zero=%s g_zero_levels=%s"
                       % (meta["d2_zero"], meta["sigma_zero"], meta["g_zero_levels"]),
             exhaustiveness_ref="MACHINE: phase_d_states_%s.json case enumeration; "
                                "frontier_rollup state_total cross-check" % win,
             exhaustiveness_level="exact-checked", machine_checkable=True,
             field_scope="char-0 / Q-algebra (all places)")

# ---- state nodes (killed states only) -------------------------------------
for ck, r in STATE.items():
    fid = r["flagcase"]
    fs = sorted(r["field_scopes"])
    add_node(dict(id="state:" + ck, type="state", window=r["window"],
                  cell=fid, canonical_key=ck, level=r["level"], closed=True,
                  n_mechanisms=len(r["mechanisms"]),
                  mechanisms=sorted(r["mechanisms"]),
                  field_scope="; ".join(fs),
                  notes="; ".join(sorted(r["notes"])),
                  evidence=sorted(r["evidence"],
                                  key=lambda e: (LRANK[e["level"]], e["source"]))))
    add_edge(fid, "state:" + ck,
             predicate="residual degree-state (deg_d1,deg_d2,deg_e,deg_sigma) "
                       "= (%s)" % ",".join(ck.split("|")[-4:]),
             exhaustiveness_ref="MACHINE: state present in phase-D universe %s"
                                % r["window"],
             exhaustiveness_level="exact-checked", machine_checkable=True,
             field_scope="; ".join(fs))
    # certificate -> state edges
    for cert in CERT_NODES:
        if ("state:" + ck) in cert["targets"]:
            add_edge("state:" + ck, cert["id"],
                     predicate="object certificate for this state's kill (%s)"
                               % cert["kill_id"],
                     exhaustiveness_ref="single certificate (lift certificate)",
                     exhaustiveness_level=cert["level"], machine_checkable=True,
                     field_scope="char-0 / Q-algebra")

# ---- alt defect-0 family branches + their state nodes ---------------------
# ground-truth slots from phase_f2_scale census
FAM_SLOTS = defaultdict(dict)      # bid -> {dd2: census_verdict}
for s in SCALE["alt_states"]:
    FAM_SLOTS[s["bid"]][s["degs"][0]] = s["verdict"]

for bid, slots in FAM_SLOTS.items():
    fam_id = "branch:altdefect0:" + bid
    add_node(dict(id=fam_id, type="branch", window="altdefect0", family=bid,
                  state_total=len(slots), level="open", closed=False))
    add_edge("subcase:sub1_alt_defect0", fam_id,
             predicate="entirely-defect-0 family %s (d1,sigma,e forced; deg_d2 free)"
                       % bid,
             exhaustiveness_ref="MACHINE: phase_f2_scale.json alt_states census "
                                "(deg_d2 in {none,0..6})",
             exhaustiveness_level="exact-checked", machine_checkable=True,
             field_scope="char-0 / Q-algebra (all places)")
    for dd2 in sorted(slots, key=lambda x: (x is not None, -1 if x is None else x)):
        sid = "state:altdefect0|%s|%s" % (bid, dd2)
        rec = DEFECT0.get((bid, dd2))
        if rec is not None:
            fs = sorted(rec["field_scopes"])
            add_node(dict(id=sid, type="state", window="altdefect0", cell=fam_id,
                          family=bid, deg_d2=dd2, level=rec["level"], closed=True,
                          n_mechanisms=len(rec["mechanisms"]),
                          mechanisms=sorted(rec["mechanisms"]),
                          field_scope="; ".join(fs), notes="; ".join(sorted(rec["notes"])),
                          evidence=sorted(rec["evidence"],
                                          key=lambda e: (LRANK[e["level"]], e["source"]))))
            add_edge(fam_id, sid,
                     predicate="deg_d2=%s slot (census=%s)" % (dd2, slots[dd2]),
                     exhaustiveness_ref="MACHINE: phase_f2_scale slot census",
                     exhaustiveness_level="exact-checked", machine_checkable=True,
                     field_scope="; ".join(fs))
            for cert in CERT_NODES:
                if sid in cert["targets"]:
                    add_edge(sid, cert["id"],
                             predicate="object certificate for defect-0 kill (%s)"
                                       % cert["kill_id"],
                             exhaustiveness_ref="single certificate",
                             exhaustiveness_level=cert["level"],
                             machine_checkable=True, field_scope="char-0 / Q-algebra")
        # surviving slots are NOT materialized as nodes; they are counted below

# ---- corner / lemma kills (auxiliary) -------------------------------------
for lab, rec in CORNER.items():
    cid = "state:corner|%s" % lab
    add_node(dict(id=cid, type="state", window="corner", cell="subcase:corner",
                  label=lab, level=rec["level"], closed=True,
                  mechanisms=sorted(rec["mechanisms"]),
                  field_scope="; ".join(sorted(rec["field_scopes"])),
                  evidence=sorted(rec["evidence"], key=lambda e: LRANK[e["level"]])))
    add_edge("subcase:corner", cid,
             predicate="structural corner/lemma kill %s" % lab,
             exhaustiveness_ref="auxiliary lemma (not a C0 state partition)",
             exhaustiveness_level="claimed", machine_checkable=False,
             field_scope="; ".join(sorted(rec["field_scopes"])))

for cert in CERT_NODES:
    add_node(cert)

# ===========================================================================
# 5. Compute closure bottom-up (CONJUNCTIVE with exhaustiveness gate)
# ===========================================================================
# cells: closed iff surviving==0; level = min(state levels) folded with the
# state->cell exhaustiveness (exact-checked).
CELL_KILLED = defaultdict(list)     # cell id -> [state level,...]
for ck, r in STATE.items():
    CELL_KILLED[r["flagcase"]].append(r["level"])

for fid, node in list(NODES.items()):
    if node.get("type") != "cell":
        continue
    total = node["state_total"]
    killed_levels = CELL_KILLED.get(fid, [])
    killed = len(killed_levels)
    node["killed"] = killed
    node["surviving"] = total - killed
    # level census within the cell
    node["level_census"] = dict(Counter(killed_levels))
    if node["surviving"] == 0 and total > 0:
        node["closed"] = True
        node["level"] = lmin("exact-checked", *killed_levels)  # exhaustiveness gate
    else:
        node["closed"] = False
        node["level"] = "open"

# branches
CELLS_OF_BRANCH = defaultdict(list)
for fid, node in NODES.items():
    if node.get("type") == "cell":
        CELLS_OF_BRANCH[node["branch"]].append(fid)

for bid, node in NODES.items():
    if node.get("type") != "branch":
        continue
    if node["window"] == "altdefect0":
        # defect-0 family: closed iff every slot killed
        bidname = node["family"]
        slots = FAM_SLOTS[bidname]
        levels = []
        killed = 0
        for dd2 in slots:
            rec = DEFECT0.get((bidname, dd2))
            if rec is not None:
                killed += 1; levels.append(rec["level"])
        node["killed"] = killed
        node["surviving"] = len(slots) - killed
        node["level_census"] = dict(Counter(levels))
        if node["surviving"] == 0 and slots:
            node["closed"] = True
            node["level"] = lmin("exact-checked", *levels)
        else:
            node["closed"] = False; node["level"] = "open"
        continue
    st = node["cascade_status"]
    if st != "survives":
        # engine-killed branch: closed by the cascade engine; the data itself
        # labels it *_pending_audit -> level 'claimed'.  Exhaustiveness edge
        # survivor_case_count==0 is machine-checked (exact-checked) but does not
        # raise the KILL evidence above the pending-audit 'claimed'.
        node["closed"] = True
        node["level"] = "claimed"
        node["closure_basis"] = "cascade engine kill (pending independent audit)"
        node["exhaustiveness_ok"] = (node.get("survivor_case_count", 0) == 0)
    else:
        cells = CELLS_OF_BRANCH.get(bid, [])
        if cells and all(NODES[c]["closed"] for c in cells):
            node["closed"] = True
            node["level"] = lmin("exact-checked", *[NODES[c]["level"] for c in cells])
        else:
            node["closed"] = False
            node["level"] = "open"
        node["n_cells"] = len(cells)
        node["cells_closed"] = sum(1 for c in cells if NODES[c]["closed"])

# subcases
BR_OF_SUB = defaultdict(list)
for bid, node in NODES.items():
    if node.get("type") != "branch":
        continue
    if node["window"] == "sub2": BR_OF_SUB["subcase:sub2"].append(bid)
    elif node["window"] == "sub1": BR_OF_SUB["subcase:sub1"].append(bid)
    elif node["window"] == "altdefect0": BR_OF_SUB["subcase:sub1_alt_defect0"].append(bid)

for sid in ("subcase:sub2", "subcase:sub1", "subcase:sub1_alt_defect0"):
    brs = BR_OF_SUB.get(sid, [])
    node = NODES[sid]
    node["n_branches"] = len(brs)
    node["branches_closed"] = sum(1 for b in brs if NODES[b]["closed"])
    node["branches_open"] = len(brs) - node["branches_closed"]
    if brs and all(NODES[b]["closed"] for b in brs):
        # exhaustiveness subcase->branch (cascade count) exact-checked, but the
        # branch->subcase completeness rests on cone lemmas; keep the fold.
        node["closed"] = True
        node["level"] = lmin("exact-checked", *[NODES[b]["level"] for b in brs])
    else:
        node["closed"] = False; node["level"] = "open"

# C0
c0_subs = ["subcase:sub2", "subcase:sub1", "subcase:sub1_alt_defect0", "subcase:f37"]
c0 = NODES["C0"]
c0["subcases_closed"] = sum(1 for s in c0_subs if NODES[s]["closed"])
if all(NODES[s]["closed"] for s in c0_subs):
    # subcase->C0 exhaustiveness is judgment (claimed) -> gate
    c0["closed"] = True
    c0["level"] = lmin("claimed", *[NODES[s]["level"] for s in c0_subs])
else:
    c0["closed"] = False; c0["level"] = "open"

# ===========================================================================
# 6. Serialize state evidence sets -> sorted lists; emit
# ===========================================================================
def _clean(node):
    out = {}
    for k, v in node.items():
        if isinstance(v, set):
            v = sorted(v)
        out[k] = v
    return out

nodes_out = [_clean(NODES[i]) for i in sorted(NODES)]
edges_out = sorted(EDGES, key=lambda e: (e["parent"], e["child"]))

# census by level and type
level_census = defaultdict(lambda: Counter())
for n_ in nodes_out:
    level_census[n_["type"]][n_.get("level", "open")] += 1

def _sources_sha256(names):
    h = hashlib.sha256()
    for name in sorted(set(names)):
        h.update(name.encode("utf-8")); h.update(b"\0")
        try:
            with open(os.path.join(HERE, name), "rb") as fh:
                h.update(fh.read())
        except OSError:
            h.update(b"<MISSING>")
        h.update(b"\0")
    return h.hexdigest()

out = {
    "schema": 1,
    "generator": "proof_dag.py",
    "description": "Coverage proof-DAG: certificate->state->cell->branch->"
                   "subcase->C0 with computed per-node evidence levels. "
                   "Closure is a computed fact; no hand roll-up asserts it.",
    "levels": LEVELS,
    "level_semantics": {
        "state": "disjunctive: MAX over independent kills/audits/certificates",
        "cell/branch/subcase/target": "conjunctive: MIN over required children, "
            "folded with the exhaustiveness-edge level (>= exact-checked required "
            "to close above claimed)",
        "engine_killed_branch": "level 'claimed' -- cascade data self-labels "
            "*_pending_audit; independent audit (C18/C29/C43) is judgment-"
            "referenced, not machine-joined in v1",
        "f37": "judgment-referenced leaf (C11); not recomputed here",
    },
    "provenance": {
        "git_commit_note": "run `git rev-parse HEAD` to bind; omitted for byte-"
                           "determinism",
        "source_sha256": _sources_sha256(LOADED),
        "source_files": sorted(set(LOADED)),
    },
    "counts": {
        "nodes": len(nodes_out),
        "edges": len(edges_out),
        "nodes_by_type": {t: len(nc) and sum(nc.values()) for t, nc in
                          sorted(level_census.items())},
        "ledger_phase_d_attributions_ingested": n_ledger_phase_d,
        "alt_hunt_kills_ingested": n_alt_hunt,
        "distinct_killed_states": sum(1 for n_ in nodes_out
                                      if n_["type"] == "state" and n_["window"] in ("sub2", "sub1")),
        "distinct_killed_defect0": len(DEFECT0),
        "corner_kills": len(CORNER),
        "certificates_total": len(CERT_NODES),
        "certificates_found": n_cert_found,
        "certificates_resolved_to_target": n_cert_resolved,
        "unmapped": len(UNMAPPED),
    },
    "closure_census": {t: dict(nc) for t, nc in sorted(level_census.items())},
    "unmapped": sorted(UNMAPPED, key=lambda u: (u.get("kind", ""),
                        str(u.get("name") or u.get("key") or u.get("canonical_key") or ""))),
    "nodes": nodes_out,
    "edges": edges_out,
}

with open(os.path.join(HERE, "proof_dag.json"), "w") as f:
    json.dump(out, f, indent=1, sort_keys=True)

print("=== proof_dag.json written ===")
print("nodes:", len(nodes_out), "edges:", len(edges_out))
print("nodes by type:", out["counts"]["nodes_by_type"])
print("distinct killed states (phase-D):", out["counts"]["distinct_killed_states"],
      "| defect-0:", out["counts"]["distinct_killed_defect0"],
      "| corner:", out["counts"]["corner_kills"])
print("certificates: total %d found %d resolved %d" % (
    len(CERT_NODES), n_cert_found, n_cert_resolved))
print("alt-hunt kills ingested:", n_alt_hunt)
print("UNMAPPED bucket:", len(UNMAPPED))
uk = Counter(u.get("kind") for u in UNMAPPED)
for k, v in sorted(uk.items()):
    print("   ", k, v)
print("C0:", c0["closed"], c0["level"], "| subcases closed",
      c0["subcases_closed"], "/", len(c0_subs))
for sid in ("subcase:sub2", "subcase:sub1", "subcase:sub1_alt_defect0"):
    s = NODES[sid]
    print("  %-28s closed=%s level=%s branches_open=%d/%d" % (
        sid, s["closed"], s["level"], s.get("branches_open"), s.get("n_branches")))
