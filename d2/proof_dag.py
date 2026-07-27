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
import json, os, re, sys, hashlib, subprocess
from collections import defaultdict, Counter

HERE = os.path.dirname(os.path.abspath(__file__))
LOADED = []

# --- optional inputs/outputs (defaults keep proof_dag.json byte-identical) ---
# `--ledger F` builds the DAG from an alternative state-kill ledger (e.g. the
# divisor-lemma-augmented state_kill_ledger_divfilter.json); `--out F` writes
# elsewhere so the committed DAG is never clobbered by an experiment.
def _argval(flag, default):
    if flag in sys.argv:
        return sys.argv[sys.argv.index(flag) + 1]
    return default

LEDGER_FILE = _argval("--ledger", "state_kill_ledger.json")
DAG_OUT = _argval("--out", "proof_dag.json")

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
    # 2026-07-26: INDEPENDENTLY-AUDITED is the ledger's grade for a kill a
    # SECOND, independently authored checker reproduced.  It is the only ledger
    # grade that maps above exact-checked, and it is set from the frontier stage
    # registry's own `level` field, never guessed here.
    "INDEPENDENTLY-AUDITED": "independently-audited",
    "AUDITED": "exact-checked",
    "TRANSFERRED-AUDITED": "exact-checked",
    "PENDING": "claimed",
    "PENDING-AMBIGUOUS-MAP": "claimed",
}

# ---------------------------------------------------------------------------
# Load sources
# ---------------------------------------------------------------------------
LED   = load(LEDGER_FILE)
CC2   = load("cascade_cones_qt_inf_rl.json")
CC1   = load("cascade_cones_sub1_qt_inf_rl.json")
D2    = load("phase_d_states_sub2.json")
D1    = load("phase_d_states_sub1.json")
SCALE = load("phase_f2_scale.json")
F2S2  = load("phase_f2_sub2.json")
AH    = load("alt_hunt_results.json")
AUC   = load("audit_alt_hunt_census.json")
# The alternate-regime (L_alt) branch universe, for the C0-partition repair of
# 2026-07-25 (DAG_REPAIR.md).  SPL1 gives the full 52-branch alternate leaf;
# ALTSW gives the 27 that survive the C33/C34 first-/second-level kills.
SPL1  = load("split_place_ledger_sub1.json")
ALTSW = load("alt_inf_sweep.json")
# The frontier stage registry (frontier_rebuild.json schema 2, 2026-07-26): the
# authoritative home for the column-level stage lists AND their evidence grades,
# plus ALT_CLOSURE -- the alternate-regime emptiness record, which has no cell
# list because the alternate regime lies outside both phase-D universes.
FR = load("frontier_rebuild.json")
ALT_CLOSURE = FR.get("alt_closure")

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
# Independent cascade-audit artifacts (the branch-level spec-only audit).
# audit_cascade_kills{,_sub1}.py exhaustively re-derive each depth-4 q-cascade
# branch's kill/survival from scratch (no code shared with the engine) and emit
# a per-branch verdict artifact.  We JOIN those verdicts here so every engine-
# killed branch the auditor CONFIRMS (audit=killed, agreement) is promoted from
# 'claimed' to 'independently-audited', recording WHICH artifact supports it.
# If the artifact is absent we invoke the auditor to produce it (slow but exact);
# normally it is present and just read (adding it to the provenance digest).
# ---------------------------------------------------------------------------
def load_audit_artifact(fn, script):
    p = os.path.join(HERE, fn)
    if not os.path.exists(p):
        sys.stderr.write("proof_dag: %s missing; running %s --quiet "
                         "--emit-artifact (exact, may take minutes)...\n" % (fn, script))
        subprocess.run([sys.executable, os.path.join(HERE, script),
                        "--quiet", "--emit-artifact", fn], cwd=HERE, check=True)
    LOADED.append(fn)
    with open(p) as f:
        return json.load(f)

CASCADE_AUDIT = {}       # (win, a_t, b_tuple, branch) -> verdict record
AUDIT_ARTIFACTS = {}     # win -> dict(file, generator_sha256, summary)
for _fn, _script in (("audit_cascade_kills.json", "audit_cascade_kills.py"),
                     ("audit_cascade_kills_sub1.json", "audit_cascade_kills_sub1.py")):
    _art = load_audit_artifact(_fn, _script)
    AUDIT_ARTIFACTS[_art["window"]] = dict(
        file=_fn, generator_sha256=_art.get("generator_sha256"),
        summary=_art.get("summary"))
    for _r in _art["branches"]:
        CASCADE_AUDIT[(_art["window"], _r["a_t"], tuple(_r["b"]), _r["branch"])] = dict(
            audit=_r["audit"], claim=_r["claim"], agreement=_r["agreement"], artifact=_fn)

# ---------------------------------------------------------------------------
# Infinity-layer audit artifact (audit_inf_cases.py, C43), joined the same way.
# The depth-4 q-cascade auditor above verdicts a branch killed only by the t/inf
# layer as 'survives' -- correctly, it sees no q-level kill -- so those branches
# get NO support from it.  audit_inf_cases.py re-derives the infinity layer from
# CASCADE_INF_REPORT.md alone, taking the q+t_rl survivor set as GIVEN, and its
# 'killed'/kill_layer='inf' records say: every one of this branch's q+t survivor
# cases was removed at infinity, and every one of those removals was re-derived.
# That is an exact re-check of the branch's emptiness RELATIVE to the q+t_rl
# baseline, not an end-to-end independent audit of it -- the q+t_rl narrowing
# itself is an engine artifact (audit_tplace_cases.py audits the kills-OFF
# cascade_cones{,_sub1}_qt.json, not the _rl one).  So this join supports
# 'exact-checked' and DELIBERATELY NOT 'independently-audited'.
# ---------------------------------------------------------------------------
INF_AUDIT = {}           # (win, a_t, b_tuple, branch) -> verdict record
INF_AUDIT_ARTIFACT = None
_ifn = "audit_inf_kills.json"
_iart = load_audit_artifact(_ifn, "audit_inf_cases.py")
INF_AUDIT_ARTIFACT = dict(file=_ifn, generator_sha256=_iart.get("generator_sha256"),
                          summary=_iart.get("summary"))
for _r in _iart["branches"]:
    INF_AUDIT[(_r["window"], _r["a_t"], tuple(_r["b"]), _r["branch"])] = dict(
        audit=_r["audit"], claim=_r["claim"], agreement=_r["agreement"],
        kill_layer=_r.get("kill_layer"),
        removed=_r.get("removed_cases_confirmed"), artifact=_ifn)

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

# the sweep each phase-D universe was cut from (recorded on the exhaustiveness
# edges so the residue-kill field scope is traceable to its artifact)
PHASE_D_SOURCE = {"sub2": D2["source_artifact"], "sub1": D1["source_artifact"]}

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

# ===========================================================================
# 1b. COLUMN-LEVEL EMPTINESS LEMMAS  (added 2026-07-26)
# ===========================================================================
# Some lemmas empty a WHOLE (a_t, b, branch) COLUMN at once -- that is the DAG's
# BRANCH granularity, not its state granularity.  They are therefore recorded
# against branches, exactly as an engine-killed branch is: a closed leaf with no
# state children required.  Writing them out as one "state kill" per enumerated
# degree-state would assert a per-state certificate that does not exist.
#
# TWO SOURCES, kept apart because their evidence grades differ:
#
#   (1) the FRONTIER STAGE REGISTRY -- frontier_rebuild.STAGES, ingested by
#       state_kill_ledger.py into `column_kills`, carrying each stage's own
#       `level`.  Stages 2-4 are exact-checked (same-author); stages 5-7
#       (a_t >= 9, a_t <= 9, the five-cell closure) are independently-audited.
#
#   (2) the `e | Phi` DIVISOR LEMMA's CELL-level consequences -- D2 (b_i >= 2 at a
#       simple q-root) and D3 (the window degree count), plus cells all of whose
#       enumerated states are D1 defect states.  Recomputed here from
#       divisor_filter.py rather than read, so it does not depend on which ledger
#       variant this DAG was built from.  GRADE: exact-checked, NOT
#       independently-audited -- FRONTIER_REBUILD.md sec.7b is explicit that the
#       divisor kills are the producing lane's own exact check, and that the
#       distinction "is load-bearing and must not be collapsed in any release".
#
# The two are DISJUNCTIVE per column (one valid lemma suffices), so a column's
# level is the MAX over its column evidence -- the same rule the state layer uses.
COLUMN = {}          # (win, cellname) -> dict(evidence=[...], level=...)


def cellname_of(a, b, br):
    return "a%d_b%s_%s" % (a, bstr(b), br)


def add_column_evidence(win, cell, level, mechanism, source, detail=None):
    r = COLUMN.setdefault((win, cell), dict(window=win, cell=cell, evidence=[],
                                            mechanisms=set(), sources=set()))
    r["evidence"].append(dict(level=level, mechanism=mechanism, source=source,
                              detail=detail or {}))
    r["mechanisms"].add(mechanism)
    r["sources"].add(source)


# (1) the stage registry, via the ledger
n_column_stage = 0
for r in LED.get("column_kills", []):
    win = r["target_window"]
    for a in r["attributions"]:
        lvl = AUDIT_LEVEL.get(a["audit"], "claimed")
        det = a.get("detail") or {}
        # the stage registry's own recorded level is authoritative; the ledger's
        # audit word is a coarser encoding of it.  Require agreement.
        if det.get("level") and det["level"] != lvl:
            UNMAPPED.append(dict(
                kind="column-level-mismatch", window=win, cell=r["cell"],
                stage=det.get("stage"), ledger_audit=a["audit"],
                mapped_level=lvl, registry_level=det["level"],
                reason="frontier stage registry level disagrees with the "
                       "ledger audit grade it was encoded as"))
        add_column_evidence(win, r["cell"], lvl, a["mechanism"], a["source"],
                            detail=det)
        n_column_stage += 1

# (2) the e | Phi divisor lemma's cell-level deaths, recomputed
n_column_divisor = 0
try:
    import divisor_filter as _dfmod
    LOADED.append("divisor_filter.py")
    for _win, _D in (("sub2", D2), ("sub1", D1)):
        _filt = _dfmod.DivisorFilter(_win)
        _cells = defaultdict(lambda: [0, 0, Counter()])   # cell -> [tot, dead, why]
        for _c in _D["cases"]:
            _nm = cellname_of(_c["a_t"], _c["b"], _c["branch"])
            for _s in _c["states"]:
                _alive, _why = _filt.state_verdict(_c["a_t"], _c["b"], _s["deg_e"])
                _cells[_nm][0] += 1
                if not _alive:
                    _cells[_nm][1] += 1
                    if _why.startswith(_dfmod.DEATH_B_GE_2):
                        _cells[_nm][2]["D2 b_i>=2 at a simple q-root"] += 1
                    elif _why.startswith(_dfmod.DEATH_DEG):
                        _cells[_nm][2]["D3 forced deg e (window count)"] += 1
                    else:
                        _cells[_nm][2]["D1 defect-0 (no off-support root)"] += 1
        for _nm, (_tot, _dead, _why) in _cells.items():
            if _tot and _dead == _tot:
                add_column_evidence(
                    _win, _nm, "exact-checked",
                    "divisor lemma e|Phi, CELL level (%s)"
                    % "; ".join("%s x%d" % kv for kv in sorted(_why.items())),
                    "divisor_filter.py / DIVISOR_SYZYGY.md sec.1",
                    detail={"states_in_cell": _tot, "states_dead": _dead,
                            "E_min": _filt.e_min, "caps": _filt.caps,
                            "grade_note": "exact-checked, NOT "
                                          "independently-audited "
                                          "(FRONTIER_REBUILD.md sec.7b)"})
                n_column_divisor += 1
except Exception as _e:                       # loud, never silent
    UNMAPPED.append(dict(kind="column-divisor-unavailable", reason=str(_e),
                         detail="divisor_filter.py could not be consulted; the "
                                "e|Phi cell-level column closures are ABSENT "
                                "from this DAG"))

for _k, _r in COLUMN.items():
    _r["level"] = "open"
    for _e in _r["evidence"]:
        _r["level"] = lmax(_r["level"], _e["level"])

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

# ---------------------------------------------------------------------------
# FIELD SCOPE of the phase-D universe  (FIELD_SCOPE_AUDIT.md, 2026-07-25)
# ---------------------------------------------------------------------------
# cascade_cones_*_qt_inf_rl.json carry `residue_kills: true`.  The two C08/C20
# forbidden rises they apply are empty over Q and over the q-splitting field
# L = Q(sqrt 17) ONLY; both supports carry real (hence complex) torus points, so
# over any K containing sqrt(105) / sqrt(170) -- in particular over R and C --
# they are NONEMPTY.  The repo's default label "char-0 / Q-algebra" does NOT
# cover this, because C is a Q-algebra: that label is exactly what breaks.
#
# What moves is the FIELD QUANTIFIER, not the evidence grade (FIELD_SCOPE_AUDIT
# J6).  The kill VERDICTS are field-stable -- 0 branches, 0 cells and 0
# independently-audited kills change status (J2) -- so branch nodes, state nodes
# and their levels keep both their grade and their char-0 scope.  What is
# field-scoped is the ENUMERATION: the flag-case list and the degree-state list
# are complete only over a K that omits sqrt(105), sqrt(170).  Those are exactly
# the branch->cell and cell->state exhaustiveness edges, and only those.
RL_SCOPE = ("FIELD-SCOPED: char-0 K with sqrt(105), sqrt(170) NOT in K "
            "(C08/C20 residue kills; cascade_cones_*_qt_inf_rl.json "
            "residue_kills=true).  NOT valid over R or C. "
            "-- FIELD_SCOPE_AUDIT.md sec.0, sec.6.1(3)")
RL_NOTE = ("FIELD-SCOPE DOWNGRADE (2026-07-25, DAG_REPAIR.md): this "
           "exhaustiveness claim is inherited from the residue-kills-ON sweep. "
           "Over an arbitrary char-0 K the enumeration is INCOMPLETE -- flag "
           "cases 220->224 (sub2) / 1145->1163 (sub1), degree-states "
           "7888->8066 / 44117->55280, i.e. a 21.8% larger phase-D universe "
           "(FIELD_SCOPE_AUDIT.md sec.4).  The returning states were never "
           "attempted, so no recorded kill is invalidated (sec.5/J2): this is a "
           "recount, not a refutation.  Evidence GRADE is unchanged (J6).")
BRANCH_RL_NOTE = ("cascade branch-enumeration completeness rests on the cone "
                  "lemmas (CASCADE_CONE_LEMMAS*.md).  FIELD SCOPE: computed "
                  "from a residue-kills-ON sweep, but the branch layer is "
                  "MEASURED field-stable -- turning the C08/C20 kills off moves "
                  "0 branches in either window (FIELD_SCOPE_AUDIT.md sec.4.2, "
                  "J2), so this edge keeps its char-0 scope.  The field scope "
                  "bites one level down, on branch->cell.")

# ---- target ---------------------------------------------------------------
add_node(dict(id="C0", type="target", label="C0: no [P,Q]=x^2 in the "
              "Prop-4.3 case-(8,28) subcases (1)&(2)", level="open", closed=False))

# ---- subcases -------------------------------------------------------------
# The five C0 leaves of the partition written out in JUDGMENT_EDGES.md sec.3 are
# L_D, L_F37, L_sub2, L_sub1, L_alt.  Before the 2026-07-25 repair the DAG's C0
# child list was {f37, sub1, sub1_alt_defect0, sub2}: L_D had no node at all, and
# L_alt was represented by `sub1_alt_defect0`, an overlay of 15 of its 27 open
# branches.  The child list is now the leaf list (DAG_REPAIR.md).
SUBCASES = {
    "subcase:sub2": dict(label="f31 window sub2 (q+t+inf, depth 4)", window="sub2"),
    "subcase:sub1": dict(label="f31 window sub1 standard (q+t+inf, depth 4)", window="sub1"),
    "subcase:sub1_alt": dict(label="L_alt: sub1 alternate regime a>=11 "
                             "(52 branches = 25 closed by C33/C34 + 27 OPEN)",
                             window="alt"),
    "subcase:sub1_alt_defect0": dict(label="sub1 alternate regime, entirely-"
                                     "defect-0 families (a>=11): a state-level "
                                     "OVERLAY of 15 of L_alt's 27 open branches, "
                                     "modelling 39 of the 4690 surviving "
                                     "alternate-regime degree-states -- NOT the "
                                     "whole alternate regime",
                                     window="altdefect0"),
    "subcase:f37": dict(label="f37 branch (closed, C11)", window="f37"),
    "subcase:dm1": dict(label="L_D: the d_{-1} == 0 branch (closed, C10)",
                        window="dm1zero"),
    "subcase:corner": dict(label="auxiliary corner/lemma kills (not a C0 partition)",
                           window="corner"),
}
for sid, meta in SUBCASES.items():
    add_node(dict(id=sid, type="subcase", level="open", closed=False, **meta))

# subcase -> C0 exhaustiveness is a JUDGMENT reference (GGHV22 Prop 4.3 + the
# field-split framework C14-16 + alt-regime C44), not machine-checkable in v1.
# It STAYS judgment after the 2026-07-25 repair, for two independent reasons --
# see JUDGMENT_EDGES.md sec.3.3 and the notes= field below.
C0_SUBS = ("subcase:sub2", "subcase:sub1", "subcase:sub1_alt", "subcase:f37",
           "subcase:dm1")
for sid in C0_SUBS:
    add_edge("C0", sid,
             predicate="C0 holds iff every regime/window subcase is closed",
             exhaustiveness_ref="JUDGMENT: GGHV22 Prop 4.3 case-(8,28) subcases "
                                "(1)&(2); field-split framework CURRENT_STATUS "
                                "C14-C16; alternate-regime partition C44.  The "
                                "partition is written out in full, split by "
                                "split with source / field scope / disjointness, "
                                "in JUDGMENT_EDGES.md sec.3 and is re-asserted "
                                "by c0_partition.py (proposition + CHECKs).  "
                                "Split 2 (d_{-1} | f31 | f37) is a declared "
                                "COVER, not a partition.",
             exhaustiveness_level="claimed", machine_checkable=False,
             field_scope="char-0 / Q-algebra; f37 restricted to char != 3,5; the "
                         "sub1/sub2 flag-case + degree-state ENUMERATION is "
                         "field-scoped (C08/C20 residue kills, see RL_SCOPE on "
                         "the branch->cell edges)",
             notes="subcase->C0 completeness is judgment-referenced, not "
                   "machine-enforced; this gates C0 closure at <= claimed, and "
                   "that gate is now the BINDING one.  UPDATED 2026-07-26: "
                   "reason (2) has expired.  It REMAINS judgment for reason (1) "
                   "ALONE -- Splits 1, 2, 5, 6, 7 rest on published/audited "
                   "mathematics that a finite bookkeeping checker cannot "
                   "re-derive.  Reason (2) used to be that the DAG's "
                   "instantiation was provably smaller than the claim: L_alt's 27 "
                   "surviving branches carried 39 modelled states against 4690 "
                   "surviving (GAP-ALT-STATES), and L_sub2's universe was 443 "
                   "minus 23 branches excised by 7 tier-3 documents "
                   "(GAP-SUB2-EXCISIONS).  BOTH ARE NOW RETIRED and recorded as "
                   "unmapped[kind=c0-partition-gap-retired]: GAP-ALT-STATES is "
                   "MOOT because `a_t <= 9` empties all 52 L_alt branches, and "
                   "GAP-SUB2-EXCISIONS is retired because all 23 excised branches "
                   "are re-derived EMPTY by in-repo lemmas that consume none of "
                   "the 7 documents.  Earlier, the 2026-07-25 repair closed "
                   "GAP-D-NONODE (L_D had no node) and GAP-ALT-BRANCHES (12 open "
                   "alt branches had no node).  All four are pinned by exact key "
                   "or exact criterion in c0_partition.py, so a retirement that "
                   "is ever reverted is a CHECK failure.  WHAT STILL CAPS C0: "
                   "this judgment edge (level 'claimed'), and subcase:dm1 (L_D, "
                   "C10) which is closed but graded 'claimed' because "
                   "PROOF_INVENTORY marks its checker attribution itself "
                   "INFERRED and no checker is wired.  "
                   "UPDATED AGAIN 2026-07-26 (second pass): the dm1 half of that "
                   "sentence has EXPIRED -- dm1_branch_verify.py is gated and "
                   "passes, so subcase:dm1 is now 'exact-checked' and all five "
                   "C0 subcases are at that level.  THIS EDGE IS NOW THE SOLE "
                   "CAP ON C0, and it is a CORRECT one, not a backlog item.  "
                   "Recorded because the standing to-do list said otherwise: it "
                   "listed 'apply the Prop 4.3 edge regrade, taking C0 from "
                   "claimed to exact-checked'.  That would be an OVERCLAIM and is "
                   "DECLINED.  prop43_audit.py (20/20) discharges the GGHV22 "
                   "Prop 4.3 CITATION -- it establishes that the cited result "
                   "says what we use it to say.  It does not, and cannot, "
                   "re-derive the exhaustiveness of the partition, which is what "
                   "'exact-checked' on this edge would assert.  Splits 1, 2, 5, 6 "
                   "and 7 rest on published mathematics; a finite bookkeeping "
                   "checker cannot reproduce it.  C0 is therefore closed at "
                   "'claimed' BY CONSTRUCTION, and the only routes above it are a "
                   "machine-checkable reformulation of the partition or a formal "
                   "proof (see LEAN_FEASIBILITY.md) -- not a regrade.")

# f37 subcase: closed by C11 (judgment/cert reference; no state data here)
f37 = NODES["subcase:f37"]
f37["closed"] = True
f37["level"] = "exact-checked"
f37["field_scope"] = "char != 3,5 (integer cert D=46875=3*5^6; F37_SATURATION_REPORT)"
f37["closure_note"] = ("closed by C11 (f37_sat_verify.py, same-author exact "
                       "checker + Lean-kernel integer certificate); represented "
                       "as a judgment-referenced leaf, not recomputed in this DAG.")

# dm1 subcase (leaf L_D of the C0 partition): closed by C10.  ADDED 2026-07-25 --
# it was a child of C0 in the proposition with NO node anywhere in the DAG
# (JUDGMENT_EDGES.md GAP-D-NONODE).  It is CLOSED, so this is a registry repair,
# not new progress; the level is deliberately 'claimed', not 'exact-checked',
# because C10 is a tier-2* row -- an upstream derivation audit whose checker
# coverage PROOF_INVENTORY.md records as *inferred* ("within
# verify_derivation.py"), not wired.  That makes L_D the weakest CLOSED leaf of
# the partition, and saying so is the point of adding the node.
dm1 = NODES["subcase:dm1"]
dm1["closed"] = True
dm1["level"] = "exact-checked"
dm1["field_scope"] = "char 0 (symbolic, denominator-free; no saturation)"
dm1["closure_note"] = (
    "closed by C10: on the reduced equations G1|_{d_{-1}=0} = 3 d_{-2} d_{-3}, "
    "so over the domain K[y] either d_{-2}=0 (then G2 -> (3/2) d_{-3}^2 forces "
    "d_{-3}=0) or d_{-3}=0 (then G3 -> -(3/2) d1 d_{-2}^2 forces d1=0); in both "
    "legs G5body+Phi collapses to Phi, contradicting Phi = f1*C4^28 != 0 "
    "(deg 238, ord 204, explicit).  Symbolic and denominator-free.  "
    "AUDIT.md sec.A.3 / STATE.md item 5 / PROOF_INVENTORY.md C10.")
dm1["evidence_tier_note"] = (
    "PROOF_INVENTORY.md grades C10 tier 2* (checker-only, and the checker "
    "attribution is itself marked [inferred: within verify_derivation.py]).  "
    "No checker for this leaf is wired into run_tests.sh, so the DAG records it "
    "at 'claimed'.  "
    "*** STATUS CHANGED 2026-07-26, AND THIS LEVEL IS NOW A DELIBERATE HOLD, NOT "
    "AN ABSENCE OF EVIDENCE. *** A checker now EXISTS: dm1_branch_verify.py plus "
    "dm1_branch_certificate.json (landed at commit 9de8713 by the C0_CLOSEOUT "
    "lane), which replaces AUDIT.md A.3's two-leg case split with a single "
    "polynomial identity over Z -- no case split, no domain argument, no side "
    "condition -- and C0_CLOSEOUT.md sec.1.6 RECOMMENDS raising this node to "
    "'independently-audited' with a floor of 'exact-checked'.  The registry lane "
    "has NOT applied that regrade: C0_CLOSEOUT.md is same-day and was not itself "
    "audited, the checker is not in run_tests.sh, and wiring an unaudited regrade "
    "into the registry is exactly the mistake that produced the v0.3.2 erratum.  "
    "So this is the ONE remaining discretionary cap on C0, and applying it is a "
    "reviewed decision, not a mechanical one.  See also the level_semantics.dm1 "
    "entry.  "
    "*** REGRADE APPLIED 2026-07-26 (second pass), TO THE FLOOR ONLY. *** Of the "
    "three blockers recorded above, one has since become STALE: dm1_branch_verify.py "
    "IS now gated in the suite (tools/suite_manifest.py; it passes 28/28), so 'no "
    "checker for this leaf is wired into run_tests.sh' is no longer true.  That "
    "clears exactly the gap between 'claimed' and 'exact-checked': an exact checker "
    "verifies the claim, and a dropped checker would now fail the suite loudly.  "
    "Level raised 'claimed' -> 'exact-checked'.  The CEILING is NOT applied and this "
    "node is NOT 'independently-audited': that needs a second, independent "
    "implementation of the same kill, and C0_CLOSEOUT.md -- which recommends it -- "
    "remains same-day and unaudited.  One gated checker is exact-checking, not "
    "independent audit; conflating the two is precisely the v0.3.2 erratum.  So C0's "
    "remaining discretionary headroom is ONE rung, not two.")
dm1["overlap_note"] = (
    "Split 2 of the C0 partition (d_{-1}=0 | f31=0 | f37=0) is a declared COVER, "
    "not a partition: L_D may overlap the f31 leaves.  Harmless -- the cover is "
    "used only in the direction 'every solution lies in some branch', and L_D "
    "and L_F37 are both empty.")

# ---- branches from the cascade enumeration --------------------------------
def cascade_summary_ok(cc):
    s = cc["summary"]
    got = s["surviving_branches"] + s["engine_killed_pending_audit"]
    return got == s["open_branches_processed"], got, s["open_branches_processed"]

# --- branch-universe completeness: the exhaustiveness edge, per window -------
# 2026-07-25 (JUDGMENT_EDGES.md sec.5 (A)/(B), cone_completeness.py).  The
# residual judgment on `branch -> subcase` was never the COUNT (already machine-
# checked) nor the kill VERDICTS (already independently audited, C18/C29) -- it
# was the UNIVERSE: that no (a,b,T) case escaped enumeration.  An independently
# authored generator (cone_completeness.py), written from the prose mathematics
# only, regenerates that universe from the terminal identities + envelope caps.
#
#   sub1: EXACT SET EQUALITY, 2178 == 2178, zero asymmetric difference in either
#         direction -> the edge is promoted exact-checked -> independently-audited.
#   sub2: the independent generator admits 443; the engine processes 420.  The
#         difference is ONE-DIRECTIONAL (zero spurious engine branches -- a
#         spurious branch would have been a soundness break) and every one of the
#         23 missing branches is excised by a NAMED earlier exact proof.  The
#         level is HELD at exact-checked, deliberately: those 7 documents are
#         tier-3 CONDITIONAL (scoped) in PROOF_INVENTORY C35, backed by
#         same-author checkers.  Promotion waits on re-running the 23.
BRANCH_EXH = {
    "sub1": dict(
        level="independently-audited",
        ref=("MACHINE: cascade_cones summary survivor+killed==open_branches_"
             "processed (%d==%d) [%s]; UNIVERSE COMPLETENESS: exact branch-key "
             "set equality with an independently authored generator derived from "
             "the terminal identities T1/T2 (3b + G = c + 2x) plus the envelope "
             "caps, written from the prose mathematics with no access to "
             "cascade_engine.py / split_place_ledger*.py / cone_lemmas.py / "
             "audit_*.py (cone_completeness.py, SUB1-EXACT, 2178==2178, zero "
             "asymmetric difference in either direction, across a 1333-stratum "
             "universe; 436 excluded branch keys all carry an explicit violated "
             "linear-budget inequality)"),
        notes=BRANCH_RL_NOTE + "  UNIVERSE: independently regenerated, exact set "
              "equality (cone_completeness.py SUB1-EXACT).  Consumed premises "
              "NOT re-derived here: P1 (q splits into four distinct degree-one "
              "places) and P5 (T3 empty), both C14-C16; and the degree caps P7 "
              "(C6/C17/C27, CAPS_AUDIT.md).",
    ),
    "sub2": dict(
        level="exact-checked",
        ref=("MACHINE: cascade_cones summary survivor+killed==open_branches_"
             "processed (%d==%d) [%s]; UNIVERSE COMPLETENESS: the independently "
             "authored generator (cone_completeness.py) admits 443 terminal-"
             "feasible branches and the engine processes 420 = 443 - 23.  The "
             "difference is one-directional: engine \\ independent = 0 (zero "
             "spurious branches), independent \\ engine = 23.  Each of the 23 is "
             "excised by a NAMED earlier exact proof across 7 documents "
             "(FIELD_SPLIT_AUDIT.md a=7 geometric q-coprime theorem; T5_60_T1; "
             "T5_60_T2; T5_90_T2; T5_STRATA_50_11 Thm 3 & Thm 4; "
             "T5_STRATUM_10_0; T5_T1_AQ12), each verified to carry ledger status "
             "proven_infeasible with a resolving .md reference.  SCOPE CHECK: "
             "PROOF_INVENTORY.md C35 retracts those proofs to geometrically-q-"
             "coprime / uniform-q^r scope, i.e. b1=b2=b3=b4, and 23/23 excised "
             "branches have uniform b [MACHINE-VERIFIED]"),
        notes=BRANCH_RL_NOTE + "  UNIVERSE: 443 - 23, HELD at exact-checked, NOT "
              "promoted.  The 23 excisions are judgment edges to 7 tier-3 "
              "CONDITIONAL (scoped) documents (PROOF_INVENTORY C35) backed by "
              "same-author checkers; C35 also records that those documents carry "
              "no conditional banner in their own text (inventory issue I5).  "
              "The excision is in-scope (uniform-b side condition machine-"
              "verified 23/23) but NOT independently audited.  Promotion to "
              "independently-audited requires re-running the 23 branches through "
              "the cascade engine -- 23 branches, uniform b, the cheapest shape "
              "in the whole universe.  See JUDGMENT_EDGES.md sec.2.7.",
    ),
}

BRANCH_STATUS = {}   # branch id -> cascade status
for win, cc, sub in (("sub2", CC2, "subcase:sub2"), ("sub1", CC1, "subcase:sub1")):
    ok, got, exp = cascade_summary_ok(cc)
    exh_ref = BRANCH_EXH[win]["ref"] % (got, exp, "PASS" if ok else "FAIL")
    exh_lvl = BRANCH_EXH[win]["level"] if ok else "claimed"
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
                 exhaustiveness_level=exh_lvl,
                 machine_checkable=True,
                 field_scope="char-0 / Q-algebra (all places)",
                 notes=BRANCH_EXH[win]["notes"])

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
                                "frontier_rollup state_total cross-check "
                                "(source_artifact %s, residue_kills=true)"
                                % (win, PHASE_D_SOURCE[win]),
             exhaustiveness_level="exact-checked", machine_checkable=True,
             field_scope=RL_SCOPE, notes=RL_NOTE)

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
             exhaustiveness_ref="MACHINE: state present in phase-D universe %s "
                                "(%s, residue_kills=true)"
                                % (r["window"], PHASE_D_SOURCE[r["window"]]),
             exhaustiveness_level="exact-checked", machine_checkable=True,
             field_scope="kill: %s | ENUMERATION: %s" % ("; ".join(fs), RL_SCOPE),
             notes=RL_NOTE)
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
             predicate="entirely-defect-0 family %s (d1,sigma,e forced; deg_d2 "
                       "free) -- models only this branch's forced-defect-0 "
                       "slots, not its full degree-state fibre" % bid,
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

# ---- L_alt: the alternate-regime leaf, all 52 branches --------------------
# ADDED 2026-07-25 (JUDGMENT_EDGES.md GAP-ALT-BRANCHES, DAG_REPAIR.md).
#
# Before the repair the DAG's only alternate-regime node was
# `subcase:sub1_alt_defect0`, whose label read as the whole alternate regime but
# which instantiates 15 of L_alt's 27 open branches.  TWELVE open alternate-
# regime branches had no node in the DAG at all.
#
# The leaf L_alt is now instantiated in full, and every branch key is derived by
# JOINING two committed artifacts, not asserted:
#     52  = split_place_ledger_sub1.json strata with
#           stratum_status == "alternate_regime_open", open_branches expanded
#     27  = alt_inf_sweep.json branches (the post-C33/C34 residual, C44's input)
#     25  = 52 \ 27, whole-branch kills from C33 (19: first-level parity/degree)
#           + C34 (6: h6/h5 levels 3/2)
#     15  = the defect-0 families above, a SUBSET of the 27 (asserted below)
#     12  = 27 \ 15, the branches with NO state model in this DAG
#
# NOTHING here closes a branch that was not already closed.  The 12 are added as
# OPEN nodes with `state_model: "none"` and an explicit count of how many
# surviving degree-states they carry that the DAG does not model.
ALT_LEAF = set()
for _r in SPL1["strata"]:
    if _r.get("stratum_status") == "alternate_regime_open":
        for _T in _r["open_branches"]:
            ALT_LEAF.add((_r["a_t"], tuple(_r["b"]), _T))
ALT_SWEEP = {(_r["a"], tuple(_r["b"]), _r["branch"]): _r for _r in ALTSW["branches"]}
ALT_OPEN = set(ALT_SWEEP)
ALT_C44_KILLED = ALT_LEAF - ALT_OPEN
ALT_DEFECT0_KEYS = set()
for _bid in FAM_SLOTS:
    _p = _bid.split("_")
    ALT_DEFECT0_KEYS.add((int(_p[0][1:]), tuple(int(_c) for _c in _p[1][1:]), _p[2]))
ALT_UNMODELLED = ALT_OPEN - ALT_DEFECT0_KEYS

# structural invariants of the join -- a FATAL here means an input artifact moved
assert ALT_OPEN <= ALT_LEAF, "alt_inf_sweep branches escape the ledger's alt leaf"
assert ALT_DEFECT0_KEYS <= ALT_OPEN, "a defect-0 family is not an open alt branch"
assert (len(ALT_LEAF), len(ALT_OPEN), len(ALT_C44_KILLED),
        len(ALT_DEFECT0_KEYS), len(ALT_UNMODELLED)) == (52, 27, 25, 15, 12), \
    "alternate-regime census moved: %d/%d/%d/%d/%d" % (
        len(ALT_LEAF), len(ALT_OPEN), len(ALT_C44_KILLED),
        len(ALT_DEFECT0_KEYS), len(ALT_UNMODELLED))

def alt_key_name(k):
    return "a%d_b%s_%s" % (k[0], bstr(k[1]), k[2])

# --- THE ALTERNATE REGIME IS EMPTY (2026-07-26) -----------------------------
# `a_t <= 9` (frontier stage6_syzygy_collision) is cap-free, branch-free and
# window-independent, so it holds in the alternate regime too.  Every branch of
# L_alt has a_t in {11..15}, all > 9, so ALL 52 are EMPTY at once -- no
# branch-by-branch argument, and no state layer, is needed.  This is a
# BRANCH-level fact and is recorded as such, exactly like the column lemmas in
# the standard windows.
#
# THIS IS WHAT RETIRES GAP-ALT-STATES.  That gap measured 39 modelled
# degree-states against 4690 surviving across the 27 open branches.  The 4690 do
# not need modelling one by one: they all sit in branches with a_t >= 11.  The gap
# is MOOT, and it is retired EXPLICITLY below rather than silently dropped.
ALT_CLOSURE_LEVEL = (ALT_CLOSURE or {}).get("level", "open")
ALT_CLOSURE_APPLIES = ALT_CLOSURE is not None
ALT_A_T_BOUND = 9        # a_t <= 9; every L_alt branch has a_t >= 11
if ALT_CLOSURE_APPLIES:
    _bad = sorted(k for k in ALT_LEAF if k[0] <= ALT_A_T_BOUND)
    assert not _bad, ("alt_closure applied to a branch with a_t <= %d: %s"
                      % (ALT_A_T_BOUND, _bad))
    ALT_CLOSURE_BASIS = (
        "EMPTY by `a_t <= %d` (frontier %s): the K-syzygy `2*Phi = e*B` is exact "
        "on the G-variety and `v_t(Phi) = 30` exactly, so `v_t(B) = 30 - a_t` "
        "EXACTLY; in unshifted h-coordinates B collapses to four terms all of "
        "whose valuations exceed `30 - a_t` for every `a_t >= 10`.  This branch "
        "has a_t = %%d >= 11.  Cap-free, branch-free, window-independent -- it "
        "holds in the alternate regime with no extra premise.  %s"
        % (ALT_A_T_BOUND, "stage6_syzygy_collision", ALT_CLOSURE["evidence"]))
else:
    ALT_CLOSURE_BASIS = None
    UNMAPPED.append(dict(
        kind="alt-closure-unavailable",
        reason="frontier_rebuild.json carries no alt_closure record; the "
               "alternate-regime branches keep their pre-2026-07-26 status"))

ALT_DIRECT = []          # the 37 branches instantiated directly under L_alt
for _k in sorted(ALT_LEAF - ALT_DEFECT0_KEYS):
    _name = alt_key_name(_k)
    _id = "branch:alt:" + _name
    ALT_DIRECT.append(_id)
    if _k in ALT_C44_KILLED:
        add_node(dict(
            id=_id, type="branch", window="alt", family=_name,
            a_t=_k[0], b=list(_k[1]), branch=_k[2],
            alt_status="killed_whole_branch", closed=True,
            level=ALT_CLOSURE_LEVEL if ALT_CLOSURE_APPLIES else "claimed",
            level_before_alt_closure="claimed",
            state_model="none (whole branch killed; no residual state fibre)",
            closure_basis="whole-branch kill in the flipped (alternate) cascade: "
                          "C33 first-level parity/degree kills (19 of 52) + C34 "
                          "levels 3/2 h6/h5 kills (6 more), leaving 27 residual "
                          "branches (13 T1 + 14 T2).  ALT_REGIME.md, "
                          "ALT_REGIME_L2.md.",
            alt_closure_level=ALT_CLOSURE_LEVEL if ALT_CLOSURE_APPLIES else None,
            alt_closure_basis=(ALT_CLOSURE_BASIS % _k[0]) if ALT_CLOSURE_APPLIES
                              else None,
            alt_closure_note=(
                "TWO INDEPENDENT ROUTES close this branch.  (1) C33/C34's "
                "whole-branch kill in the flipped cascade -- tier 2, "
                "'claimed', because audit_alt_regime.py re-derives all 25 "
                "spec-only but emits no joinable artifact.  (2) `a_t <= 9`, "
                "independently-audited (at_le9_audit.py 76/76), which empties "
                "every a_t >= 11 branch outright.  Branch closure is "
                "DISJUNCTIVE -- one valid kill suffices -- so the node's level "
                "is route (2)'s, and route (1) is retained as the historical "
                "basis.  Route (2) does NOT depend on the flipped cascade, on "
                "the place dichotomy, or on any degree cap."
                if ALT_CLOSURE_APPLIES else None),
            audit_note="PROOF_INVENTORY.md grades C33/C34 tier 2 (same-author "
                       "checkers alt_regime_verify.py / alt_regime_l2_verify.py) "
                       "and C44 tier 1 -- audit_alt_regime.py is a Codex-authored "
                       "spec-only auditor that re-derives all 25 of these kills "
                       "by name and asserts the residual 27, and it runs in "
                       "run_tests.sh.  It emits NO joinable artifact, so this DAG "
                       "cannot machine-join the promotion and records the kill at "
                       "'claimed'.  This is the same discipline the engine-killed "
                       "sub1/sub2 branches follow.  UPGRADE: give "
                       "audit_alt_regime.py an --emit-artifact flag and join it "
                       "here, exactly as audit_cascade_kills{,_sub1}.json are "
                       "joined -> 'independently-audited'.",
            derivation="split_place_ledger_sub1.json alternate_regime_open (52) "
                       "MINUS alt_inf_sweep.json branches (27)"))
        _pred = ("alternate-regime branch %s: killed whole by the flipped "
                 "cascade's first/second level (C33+C34)" % _name)
    else:
        _sw = ALT_SWEEP[_k]
        add_node(dict(
            id=_id, type="branch", window="alt", family=_name,
            a_t=_k[0], b=list(_k[1]), branch=_k[2],
            alt_status=("closed_by_a_t_bound" if ALT_CLOSURE_APPLIES
                        else "open_unmodelled"),
            closed=bool(ALT_CLOSURE_APPLIES),
            level=ALT_CLOSURE_LEVEL if ALT_CLOSURE_APPLIES else "open",
            level_before_alt_closure="open",
            state_model="none (not needed: the whole branch is empty)"
                        if ALT_CLOSURE_APPLIES else "none",
            modelled_states=0,
            alt_closure_level=ALT_CLOSURE_LEVEL if ALT_CLOSURE_APPLIES else None,
            alt_closure_basis=(ALT_CLOSURE_BASIS % _k[0]) if ALT_CLOSURE_APPLIES
                              else None,
            alt_degree_states_total=_sw["counts"]["total_degree_states"],
            alt_degree_states_killed=_sw["counts"]["killed"],
            alt_degree_states_surviving=_sw["counts"]["surviving"],
            coverage_note=(
                ("CLOSED by `a_t <= %d`, and the state layer is no longer "
                 "needed.  Before 2026-07-26 this branch was OPEN with "
                 "state_model 'none': the C44 max-plus degree sweep killed %d of "
                 "its %d degree-states and left %d SURVIVING, none of them "
                 "modelled here.  Those %d survivors do not need to be modelled "
                 "one by one -- they all live at a_t = %d >= 11, and the bound "
                 "empties the branch.  This is a genuine CLOSURE, not a registry "
                 "repair: the obligation is discharged, not merely registered."
                 % (ALT_A_T_BOUND, _sw["counts"]["killed"],
                    _sw["counts"]["total_degree_states"],
                    _sw["counts"]["surviving"], _sw["counts"]["surviving"], _k[0]))
                if ALT_CLOSURE_APPLIES else
                ("OPEN.  The C44 max-plus degree sweep kills %d of this branch's "
                 "%d degree-states and leaves %d SURVIVING, none of which is "
                 "modelled as a state node in this DAG: this branch has no "
                 "phase-D cell layer and no defect-0 family.  The node exists so "
                 "that the obligation is REGISTERED, not because anything about "
                 "it is closed."
                 % (_sw["counts"]["killed"],
                    _sw["counts"]["total_degree_states"],
                    _sw["counts"]["surviving"]))),
            derivation="alt_inf_sweep.json branches (27) MINUS the "
                       "phase_f2_scale.json defect-0 families (15)"))
        _pred = ("alternate-regime branch %s: %s"
                 % (_name,
                    "EMPTY by a_t <= %d; no state model needed" % ALT_A_T_BOUND
                    if ALT_CLOSURE_APPLIES
                    else "OPEN after C33/C34/C44, with no state model in this DAG"))
    add_edge("subcase:sub1_alt", _id, predicate=_pred,
             exhaustiveness_ref="MACHINE: L_alt branch universe = "
                                "split_place_ledger_sub1.json strata with "
                                "stratum_status=alternate_regime_open, "
                                "open_branches expanded (52 branches over 26 "
                                "strata); independently regenerated with EXACT "
                                "SET EQUALITY by cone_completeness.py (alt: "
                                "independent 52 | ledger 52 | sets EQUAL) and "
                                "re-asserted by c0_partition.py ALT-UNIVERSE / "
                                "ALT-LEDGER.  The 52 = 25 + 27 split is "
                                "alt_inf_sweep.json's branch list.",
             exhaustiveness_level="independently-audited", machine_checkable=True,
             field_scope="char-0 / Q-algebra (all places)",
             notes="the alternate-regime UNIVERSE (52 branches, 26 strata) is "
                   "independently regenerated and exactly equal to the ledger's; "
                   "the 52 -> 27 REDUCTION is consumed from C33/C34/C44, not "
                   "re-derived here.")

# the defect-0 overlay is a refinement of 15 of L_alt's 27 open branches
add_edge("subcase:sub1_alt", "subcase:sub1_alt_defect0",
         predicate="the forced-defect-0 state overlay of 15 of L_alt's 27 open "
                   "branches (%s)"
                   % ", ".join(sorted(alt_key_name(k) for k in ALT_DEFECT0_KEYS)),
         exhaustiveness_ref="MACHINE: phase_f2_scale.json alt_states census; the "
                            "15 family keys are verified to be a SUBSET of "
                            "alt_inf_sweep.json's 27 open branches "
                            "(c0_partition.py ALT-DAG-IN)",
         exhaustiveness_level="claimed", machine_checkable=False,
         field_scope="char-0 / Q-algebra (all places)",
         notes="NOT exhaustive, and it never was: the overlay models 39 "
               "forced-defect-0 slots against 4690 surviving alternate-regime "
               "degree-states across the 27 branches that survived C33/C34.  "
               "Closing all 15 families would NOT close those 15 branches -- only "
               "their defect-0 slots.  That shortfall was GAP-ALT-STATES.  "
               "RETIRED 2026-07-26: the shortfall is MOOT, because `a_t <= 9` "
               "empties all 52 L_alt branches outright, so the 4690 never needed "
               "modelling.  The overlay is kept as the historical state-level "
               "record; it is no longer load-bearing for closure.  See "
               "unmapped[kind=c0-partition-gap-retired].")

_alt = NODES["subcase:sub1_alt"]
_alt.update(
    n_branches=len(ALT_LEAF),
    # branches_closed / branches_open are the LIVE roll-up.  Before 2026-07-26
    # they read 25 / 27 -- the C33/C34 split -- which became a stale live-looking
    # obligation the moment a_t <= 9 closed the leaf.  The C33/C34 split is kept
    # under its own explicit names below.
    branches_closed=(len(ALT_LEAF) if ALT_CLOSURE_APPLIES
                     else len(ALT_C44_KILLED)),
    branches_open=(0 if ALT_CLOSURE_APPLIES
                   else len(ALT_OPEN)),
    branches_killed_whole=len(ALT_C44_KILLED),
    branches_surviving_C33_C34=len(ALT_OPEN),
    surviving_C33_C34_with_state_overlay=len(ALT_DEFECT0_KEYS),
    surviving_C33_C34_unmodelled=len(ALT_UNMODELLED),
    unmodelled_branch_keys=sorted(alt_key_name(k) for k in ALT_UNMODELLED),
    states_modelled=sum(len(s) for s in FAM_SLOTS.values()),
    states_surviving=ALTSW["summary"]["surviving_states"],
    branches_closed_by_a_t_bound=len(ALT_LEAF) if ALT_CLOSURE_APPLIES else 0,
    branches_survived_C33_C34=len(ALT_OPEN),
    closure_mechanism=("a_t <= %d (stage6_syzygy_collision, "
                       "independently-audited)" % ALT_A_T_BOUND)
                      if ALT_CLOSURE_APPLIES else None,
    field_scope="char-0 / Q-algebra (all places)",
    coverage_note=(
        ("L_alt is CLOSED, 52/52 branches, as of 2026-07-26.  The mechanism is "
         "`a_t <= %d` -- cap-free, branch-free, window-independent -- and every "
         "L_alt branch has a_t in {11..15}, so all 52 die at once.  "
         "HISTORICALLY: %d were already closed whole by C33+C34 at 'claimed', "
         "and %d survived them, of which %d carried a forced-defect-0 state "
         "overlay (%d slots) and %d carried no state model at all.  Across "
         "those %d the DAG modelled %d states against %d surviving degree-states "
         "-- the shortfall declared as GAP-ALT-STATES.  That gap is now MOOT, "
         "not repaired: the %d survivors never needed modelling, because they "
         "all live at a_t >= 11.  DO NOT read branches_closed (whole-branch "
         "C33/C34 kills, %d) as the closure count; the closure count is 52 and "
         "the field is branches_closed_by_a_t_bound."
         % (ALT_A_T_BOUND, len(ALT_C44_KILLED), len(ALT_OPEN),
            len(ALT_DEFECT0_KEYS), sum(len(s) for s in FAM_SLOTS.values()),
            len(ALT_UNMODELLED), len(ALT_OPEN),
            sum(len(s) for s in FAM_SLOTS.values()),
            ALTSW["summary"]["surviving_states"],
            ALTSW["summary"]["surviving_states"], len(ALT_C44_KILLED)))
        if ALT_CLOSURE_APPLIES else
        ("L_alt is REGISTERED in full (52/52 branch keys) but COVERED only in "
         "part.  %d branches are closed whole (C33+C34, recorded at 'claimed' "
         "pending an emitted audit artifact); %d are OPEN, of which %d carry a "
         "forced-defect-0 state overlay and %d carry no state model at all.  "
         "Across the 27 open branches the DAG models %d states against %d "
         "surviving degree-states (GAP-ALT-STATES)."
         % (len(ALT_C44_KILLED), len(ALT_OPEN), len(ALT_DEFECT0_KEYS),
            len(ALT_UNMODELLED), sum(len(s) for s in FAM_SLOTS.values()),
            ALTSW["summary"]["surviving_states"]))))

_d0 = NODES["subcase:sub1_alt_defect0"]
_d0.update(
    models_open_branches_of_L_alt=len(ALT_DEFECT0_KEYS),
    of_open_branches=len(ALT_OPEN),
    states_modelled=sum(len(s) for s in FAM_SLOTS.values()),
    alt_states_surviving=ALTSW["summary"]["surviving_states"],
    scope_note=(
        "This node is an OVERLAY, not a partition member: it is the "
        "forced-defect-0 slot model of %d of the %d L_alt branches that survived "
        "C33/C34, %d slots against %d surviving degree-states.  Its parent is "
        "subcase:sub1_alt (the real L_alt leaf), not C0.  Closing it does NOT by "
        "itself close those %d branches -- only their defect-0 slots.  What "
        "closes those branches is `a_t <= %d`, recorded on the branch nodes "
        "themselves.  Every slot here now carries that kill too, which is why "
        "this node closes; it is a CONSEQUENCE of the branch-level bound, not an "
        "independent achievement of the overlay."
        % (len(ALT_DEFECT0_KEYS), len(ALT_OPEN),
           sum(len(s) for s in FAM_SLOTS.values()),
           ALTSW["summary"]["surviving_states"], len(ALT_DEFECT0_KEYS),
           ALT_A_T_BOUND)))

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

# ---- the RESIDUAL C0-partition coverage gaps, by key ----------------------
# JUDGMENT_EDGES.md sec.5 (C)(3): record the gaps in `unmapped` so the report
# surfaces them by KEY, not by prose.  Two of the four gaps found by that lane
# are closed by this repair (GAP-D-NONODE, GAP-ALT-BRANCHES) and are therefore
# absent here.  These two are NOT closed by adding nodes, and are not claimed to be.
# ---- GAP RETIREMENTS (2026-07-26) ------------------------------------------
# Both declared C0-partition gaps are retired here.  They are recorded as
# `c0-partition-gap-retired`, NOT deleted: a gap that quietly disappears is
# indistinguishable from a gap that was quietly ignored.  c0_partition.py pins
# these records and fails if either one vanishes or reverts.
#
UNMAPPED.append(dict(
    kind="c0-partition-gap-retired", name="GAP-ALT-STATES", leaf="L_alt",
    status="RETIRED 2026-07-26 -- MOOT, superseded by a stronger result",
    retired_by="a_t <= %d (frontier stage6_syzygy_collision, "
               "independently-audited via at_le9_audit.py 76/76)" % ALT_A_T_BOUND,
    former_modelled_states=sum(len(s) for s in FAM_SLOTS.values()),
    former_surviving_states=ALTSW["summary"]["surviving_states"],
    former_open_branches=len(ALT_OPEN),
    former_unmodelled_branch_keys=sorted(alt_key_name(k) for k in ALT_UNMODELLED),
    branches_now_closed=len(ALT_LEAF) if ALT_CLOSURE_APPLIES else 0,
    reason="The gap was: across L_alt's %d branches surviving C33/C34 the DAG "
           "modelled %d degree-states (the forced-defect-0 slots of "
           "phase_f2_scale.json) against %d surviving (alt_inf_sweep.json), so "
           "the leaf was REPRESENTED rather than COVERED, and the 2026-07-25 "
           "repair could not shrink it.  IT IS NOW MOOT.  `a_t <= %d` is "
           "cap-free, branch-free and window-independent; every L_alt branch has "
           "a_t in {11..15}, all > %d; so all 52 branches are EMPTY and the %d "
           "surviving degree-states never needed to be modelled one by one.  The "
           "gap is not repaired by registry work -- it is dissolved by a later "
           "theorem, which is why it is retired rather than closed."
           % (len(ALT_OPEN), sum(len(s) for s in FAM_SLOTS.values()),
              ALTSW["summary"]["surviving_states"], ALT_A_T_BOUND, ALT_A_T_BOUND,
              ALTSW["summary"]["surviving_states"])))
# GAP-SUB2-EXCISIONS.  The gap was that L_sub2's 420-branch engine universe is
# the independent generator's 443 terminal-feasible branches MINUS 23, and those
# 23 removals were judgment edges to 7 tier-3 CONDITIONAL documents.  It is
# retired because the 23 are now re-derived INSIDE this repo by machine-checked
# lemmas that consume none of those documents.  The CRITERION is recorded here;
# c0_partition.py recomputes the 23 keys from the independent generator and
# verifies the classification KEY BY KEY, so this record is pinned, not prose.
#
#   `a_t = 9` EXACTLY (stages 5 + 6, both independently-audited) kills every
#   excised branch with a_t != 9.  Exactly ONE excised branch has a_t = 9,
#   namely a9_b0000_T2, and it dies on sub2's own degree identity: the e|Phi
#   divisor lemma forces deg e = 10 = a_t + sum(b_i) in sub2 (E_min = e_cap = 10),
#   so a_t = 9 forces sum(b_i) = 1, and that branch has sum(b_i) = 0.  That last
#   step is exact-checked (divisor_syzygy.py / divisor_filter.py), NOT
#   independently-audited, and the record says so.
UNMAPPED.append(dict(
    kind="c0-partition-gap-retired", name="GAP-SUB2-EXCISIONS", leaf="L_sub2",
    status="RETIRED 2026-07-26 -- the 23 excisions are no longer judgment edges",
    terminal_feasible=443, engine_universe=420, excised=23, documents=7,
    retired_by="a_t = 9 EXACTLY (stage5_slice_obstruction + "
               "stage6_syzygy_collision, both independently-audited) for 22 of "
               "the 23; sub2's degree identity deg e = 10 = a_t + sum(b_i) "
               "(e|Phi divisor lemma D1+D3, exact-checked) for the 23rd",
    criterion=dict(
        a_t_exact=9,
        rule_1="every excised branch with a_t != 9 is EMPTY by a_t = 9 -- "
               "cap-free, branch-free and window-independent, so it applies to "
               "branches the engine never processed",
        rule_2="an excised branch with a_t == 9 must additionally satisfy "
               "sum(b_i) == 1, because sub2's e|Phi degree count forces "
               "deg e = 10 = a_t + sum(b_i); any excised branch with a_t == 9 "
               "and sum(b_i) != 1 is EMPTY by that identity",
        rule_1_level="independently-audited",
        rule_2_level="exact-checked",
        expected_killed_by_rule_1=22,
        expected_killed_by_rule_2=1,
        expected_rule_2_keys=["a9_b0000_T2"],
        expected_uncovered=0),
    weakest_link="exact-checked -- the single excised branch with a_t = 9 "
                 "(a9_b0000_T2) needs sub2's degree identity from the e|Phi "
                 "divisor lemma, which is the producing lane's own exact check "
                 "(FRONTIER_REBUILD.md sec.7b), NOT an independent audit.  The "
                 "other 22 are independently-audited.",
    reason="The gap was: those 23 removals are sound only if 7 tier-3 "
           "CONDITIONAL documents are (FIELD_SPLIT_AUDIT.md, T5_60_T1.md, "
           "T5_60_T2.md, T5_90_T2.md, T5_STRATA_50_11.md, T5_STRATUM_10_0.md, "
           "T5_T1_AQ12.md), each backed by a same-author checker, and none "
           "carrying a conditional banner in its own text (PROOF_INVENTORY C35, "
           "inventory issue I5).  IT IS NOW RETIRED.  All 23 excised branches are "
           "re-derived EMPTY by in-repo machine-checked lemmas that consume none "
           "of those documents, so the excision no longer rests on them -- "
           "whether or not the 23 are ever re-run through the cascade engine.  "
           "Note the direction of the argument: not that the 23 were "
           "re-processed, but that they are empty for a reason that does not care "
           "whether they were.  See JUDGMENT_EDGES.md sec.2.7 for the gap as "
           "originally declared."))
UNMAPPED.append(dict(
    kind="field-scope", name="PHASE-D-UNIVERSE-Q-SCOPED",
    status="OPEN -- field quantifier, not evidence grade",
    affects="branch->cell and cell->state exhaustiveness edges, both windows",
    reason="The phase-D universe is cut from cascade_cones_*_qt_inf_rl.json "
           "(residue_kills=true).  The C08/C20 forbidden rises those sweeps "
           "apply are empty over Q and over L=Q(sqrt 17) only; both supports "
           "carry real torus points, so over any K containing sqrt(105) or "
           "sqrt(170) -- in particular R and C -- they are nonempty.  The "
           "ENUMERATION is therefore incomplete over an arbitrary char-0 K: "
           "220->224 sub2 / 1145->1163 sub1 flag cases, 7888->8066 / "
           "44117->55280 degree-states (+21.8% universe).  No recorded kill is "
           "invalidated (the returning states were never attempted) and 0 "
           "branches, 0 cells and 0 independently-audited kills change status: "
           "this is a RECOUNT, not a refutation, and a FIELD-SCOPE downgrade, "
           "not an evidence-grade downgrade.  FIELD_SCOPE_AUDIT.md sec.0, "
           "sec.4, sec.6.1(3), J2/J6."))

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
    # COLUMN LEMMA (2026-07-26).  A lemma that empties the whole (a_t, b, branch)
    # column empties this flag-case with it.  This is DISJUNCTIVE with the
    # per-state roll-up above -- one valid route suffices -- so take the MAX, and
    # record honestly how many of the cell's states have an individual kill record
    # and how many are covered only by the column lemma.
    _col = COLUMN.get((node["window"],
                       cellname_of(node["a_t"], node["b"], node["branch_label"])))
    if _col is not None:
        node["column_lemma_level"] = _col["level"]
        node["column_lemma_mechanisms"] = sorted(_col["mechanisms"])
        node["states_with_individual_kill_record"] = killed
        node["states_covered_only_by_column_lemma"] = node["surviving"]
        node["closure_basis"] = (
            "the whole (a_t, b, branch) COLUMN is empty by %s.  %d of this "
            "flag-case's %d degree-states also carry an individual kill record; "
            "the other %d are covered by the column lemma alone -- there is no "
            "per-state certificate for them, and this field says so."
            % (", ".join(sorted(_col["mechanisms"])), killed, total,
               node["surviving"]))
        if node["closed"]:
            node["level"] = lmax(node["level"], _col["level"])
        else:
            node["closed"] = True
            node["level"] = _col["level"]

# branches
CELLS_OF_BRANCH = defaultdict(list)
for fid, node in NODES.items():
    if node.get("type") == "cell":
        CELLS_OF_BRANCH[node["branch"]].append(fid)

for bid, node in NODES.items():
    if node.get("type") != "branch":
        continue
    if node["window"] == "alt":
        # L_alt branch: closure was set at construction from the artifact join
        # (killed-whole by C33/C34, or OPEN with no state model).  Nothing to
        # roll up -- these nodes have no children by design, and that absence is
        # exactly what they record.
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
        # engine-killed branch.  JOIN the independent spec-only cascade audit:
        # if the auditor artifact confirms this branch killed (audit=killed AND
        # agreement), the kill is independently re-derived -> 'independently-
        # audited'.  Otherwise it stays 'claimed' (cascade data self-labels it
        # *_pending_audit).  Branches killed only by the t/inf layer are NOT
        # covered by the q-cascade auditor (it verdicts them 'survives' at the
        # q level) and honestly remain 'claimed'.
        node["closed"] = True
        akey = (node["window"], node["a_t"], tuple(node["b"]), node["branch"])
        aud = CASCADE_AUDIT.get(akey)
        if aud is not None and aud["audit"] == "killed" and aud["agreement"]:
            node["level"] = "independently-audited"
            node["closure_basis"] = (
                "cascade engine kill, INDEPENDENTLY AUDITED (spec-only "
                "re-derivation agrees) by %s" % aud["artifact"])
            node["audited_by"] = aud["artifact"]
            node["auditor_sha256"] = AUDIT_ARTIFACTS[node["window"]]["generator_sha256"]
        else:
            # Not covered by the q-cascade auditor.  JOIN the infinity-layer
            # audit: if audit_inf_cases.py re-derived the removal of EVERY one of
            # this branch's q+t survivor cases at infinity, the branch's
            # emptiness is exactly re-checked relative to the q+t_rl baseline ->
            # 'exact-checked'.  Not 'independently-audited': the _rl narrowing of
            # the q+t frontier that this auditor takes as given is not itself
            # independently audited (see the INF_AUDIT loader comment).
            inf = INF_AUDIT.get(akey)
            if (inf is not None and inf["audit"] == "killed" and inf["agreement"]
                    and inf["kill_layer"] == "inf" and (inf["removed"] or 0) > 0):
                node["level"] = "exact-checked"
                node["closure_basis"] = (
                    "cascade engine kill at the INFINITY layer, EXACT-CHECKED: "
                    "%s re-derives the removal of all %d of this branch's q+t "
                    "survivor cases at infinity (spec-only, conservative relaxed "
                    "semantics), leaving 0 survivor cases"
                    % (inf["artifact"], inf["removed"]))
                node["inf_audited_by"] = inf["artifact"]
                node["inf_auditor_sha256"] = INF_AUDIT_ARTIFACT["generator_sha256"]
                node["inf_removed_cases_confirmed"] = inf["removed"]
                node["audit_note"] = (
                    "q-cascade auditor verdict=survives (it sees no q-level kill, "
                    "correctly); the kill is at the infinity layer and is joined "
                    "from %s.  Capped at 'exact-checked' rather than "
                    "'independently-audited' because the q+t_rl survivor set that "
                    "auditor starts from is taken as given (audit_tplace_cases.py "
                    "audits the kills-OFF q+t frontier, not the _rl one)."
                    % inf["artifact"])
            else:
                node["level"] = "claimed"
                node["closure_basis"] = "cascade engine kill (pending independent audit)"
                if aud is not None:
                    node["audit_note"] = (
                        "q-cascade auditor verdict=%s (agreement=%s): this branch's "
                        "kill is not covered by the depth-4 q-cascade auditor "
                        "(killed only by the t/inf layer)%s"
                        % (aud["audit"], aud["agreement"],
                           "" if inf is None else
                           "; infinity auditor verdict=%s (kill_layer=%s), which "
                           "does not support this branch either"
                           % (inf["audit"], inf["kill_layer"])))
                else:
                    node["audit_note"] = "no q-cascade auditor record for this branch key"
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
        # COLUMN LEMMA (2026-07-26).  This branch IS a column, so a column lemma
        # closes it as a LEAF -- exactly as an engine-killed branch is closed,
        # with no state-child requirement.  It is NOT capped by the
        # lmin("exact-checked", ...) fold above, because that fold is the
        # branch->cell->state exhaustiveness gate and a column lemma does not go
        # through the state enumeration at all: it says the column is empty over
        # any K, for every degree assignment, enumerated or not.
        _col = COLUMN.get((node["window"],
                           cellname_of(node["a_t"], node["b"], node["branch"])))
        if _col is not None:
            node["column_lemma_level"] = _col["level"]
            node["column_lemma_mechanisms"] = sorted(_col["mechanisms"])
            node["column_lemma_sources"] = sorted(_col["sources"])
            node["column_lemma_evidence"] = sorted(
                _col["evidence"], key=lambda e: (LRANK[e["level"]], e["source"]))
            node["closure_basis"] = (
                "COLUMN-LEVEL EMPTINESS LEMMA: the whole (a_t=%d, b=%s, %s) column "
                "is empty.  Closed as a leaf -- the lemma is not a statement about "
                "the enumerated degree-states, so it is not gated by the "
                "state-enumeration exhaustiveness edge.  Mechanisms: %s"
                % (node["a_t"], bstr(node["b"]), node["branch"],
                   "; ".join(sorted(_col["mechanisms"]))))
            node["closed"] = True
            node["level"] = lmax(node["level"], _col["level"])

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

# L_alt: 37 branch children (25 killed whole + 12 open-unmodelled) PLUS the
# defect-0 overlay subcase, which carries the remaining 15 open branches.
# Closure stays a COMPUTED fact -- the 12 unmodelled branches are closed=False,
# so L_alt is open no matter what happens to the overlay.
_alt_children = sorted(ALT_DIRECT) + ["subcase:sub1_alt_defect0"]
_alt = NODES["subcase:sub1_alt"]
_alt["n_child_obligations"] = len(_alt_children)
_alt["child_obligations_closed"] = sum(1 for c in _alt_children if NODES[c]["closed"])
if all(NODES[c]["closed"] for c in _alt_children):
    _alt["closed"] = True
    _alt["level"] = lmin("exact-checked", *[NODES[c]["level"] for c in _alt_children])
else:
    _alt["closed"] = False; _alt["level"] = "open"

# C0 -- the child list is now the leaf list of JUDGMENT_EDGES.md sec.3:
#   L_sub2, L_sub1, L_alt, L_F37, L_D.
c0_subs = list(C0_SUBS)
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
                   "Closure is a computed fact; no hand roll-up asserts it. "
                   "C0's children are the five leaves of the case partition "
                   "written out in JUDGMENT_EDGES.md sec.3 (L_sub2, L_sub1, "
                   "L_alt, L_F37, L_D); the subcase->C0 edge itself stays "
                   "JUDGMENT.  Repaired 2026-07-25, see DAG_REPAIR.md.",
    "levels": LEVELS,
    "level_semantics": {
        "state": "disjunctive: MAX over independent kills/audits/certificates",
        "cell/branch/subcase/target": "conjunctive: MIN over required children, "
            "folded with the exhaustiveness-edge level (>= exact-checked required "
            "to close above claimed)",
        "engine_killed_branch": "level 'independently-audited' when the spec-only "
            "cascade auditor (audit_cascade_kills{,_sub1}.py) artifact confirms the "
            "kill (audit=killed, agreement); else 'claimed'.  Branches killed only "
            "by the t/inf layer are outside the q-cascade auditor's scope and stay "
            "'claimed' (see node.audit_note / audited_by).",
        "f37": "judgment-referenced leaf (C11); not recomputed here",
        "field_scope_recount": "the enumerated universe is residue-kills-ON. The "
            "11341 states a kills-OFF rebuild would add carry NO kill record -- "
            "they were never attempted, because they were never in the universe "
            "-- so they would enter at level 'open', closed=false, with no "
            "evidence entries.  They are NOT reopened kills, and no existing "
            "node's level drops (state_kill_ledger.py; FIELD_SCOPE_AUDIT.md J2; "
            "FIELD_SCOPE_REPAIR.md sec.4.1).",
        "dm1":"judgment-referenced leaf (C10, tier 2* -- checker attribution "
               "itself INFERRED); closed at level 'claimed', the weakest closed "
               "leaf of the C0 partition",
        "alt": "L_alt branch nodes carry NO state layer.  A 'killed_whole_branch' "
               "node is closed at 'claimed' (C33/C34; audit_alt_regime.py "
               "re-derives all 25 spec-only but emits no joinable artifact).  An "
               "'open_unmodelled' node is OPEN with state_model='none': the node "
               "registers an obligation the DAG does not model, and its "
               "alt_degree_states_surviving field says how much.",
    },
    # the key FIELD_SCOPE_REPAIR.md sec.4.1 specifies, with its wording
    "field_scope_note":
        "universe is Q-SCOPED: built from cascade_cones_{qt,sub1_qt}_inf_rl.json "
        "(residue_kills=true). C08/C20 are CONSTRAINTS, not kills "
        "(FIELD_SCOPE_AUDIT.md); the char-0 universe is 1387 flag cases / 63346 "
        "states, +22 cases / +11341 states, all entering as OPEN.",
    "field_scope": {
        "default": "char-0 / Q-algebra (all places)",
        "WARNING": "the default label does NOT imply 'over every char-0 K'. "
                   "C is a Q-algebra, and two enumeration layers are strictly "
                   "narrower than that -- see the entries below and "
                   "unmapped[kind=field-scope].",
        "phase_d_enumeration": RL_SCOPE,
        "f37_leaf": "char != 3,5 (integer certificate D = 46875 = 3*5^6)",
        "node_levels": "UNCHANGED, deliberately.  No node in this DAG uses C08 "
                       "or C20 as evidence -- the residue kills appear in no "
                       "node's evidence[].mechanism (enumerated in "
                       "FIELD_SCOPE_REPAIR.md sec.4.0).  They enter in exactly "
                       "one way: as part of the definition of the UNIVERSE this "
                       "DAG enumerates.  Lowering an evidence tier would be the "
                       "wrong repair (FIELD_SCOPE_AUDIT J6).  What is relabelled "
                       "is the branch->cell and cell->state EXHAUSTIVENESS "
                       "edges, i.e. the universe itself.",
        "audit": "FIELD_SCOPE_AUDIT.md and FIELD_SCOPE_REPAIR.md sec.4.1 "
                 "(2026-07-25); applied to this registry in DAG_REPAIR.md.  The "
                 "affected nodes keep their evidence GRADE; only the field "
                 "quantifier moves (FIELD_SCOPE_AUDIT J6).",
    },
    "provenance": {
        "git_commit_note": "run `git rev-parse HEAD` to bind; omitted for byte-"
                           "determinism",
        "source_sha256": _sources_sha256(LOADED),
        "source_files": sorted(set(LOADED)),
        "cascade_audit_join": {w: AUDIT_ARTIFACTS[w] for w in sorted(AUDIT_ARTIFACTS)},
        "inf_audit_join": INF_AUDIT_ARTIFACT,
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
        # REFRESHED 2026-07-26.  The old shape of this block ("branches_open: 27,
        # states_surviving: 4690, states_modelled: 39") went STALE the moment
        # a_t <= 9 landed: it read as a live obligation when the leaf was empty.
        # The historical numbers are kept but renamed `..._C33_C34` /
        # `former_...` so they cannot be misread as current status, and
        # `branches_open` is now a COMPUTED roll-up over the branch nodes.
        "alt_leaf": {
            "branches_total": len(ALT_LEAF),
            "branches_closed": sum(1 for n_ in NODES.values()
                                   if n_.get("type") == "branch"
                                   and n_.get("window") == "alt"
                                   and n_.get("closed")) + sum(
                1 for n_ in NODES.values()
                if n_.get("type") == "branch"
                and n_.get("window") == "altdefect0" and n_.get("closed")),
            "branches_open": sum(1 for n_ in NODES.values()
                                 if n_.get("type") == "branch"
                                 and n_.get("window") in ("alt", "altdefect0")
                                 and not n_.get("closed")),
            "closure_mechanism": ("a_t <= %d (stage6_syzygy_collision, "
                                  "independently-audited)" % ALT_A_T_BOUND)
                                 if ALT_CLOSURE_APPLIES else None,
            "closure_level": ALT_CLOSURE_LEVEL if ALT_CLOSURE_APPLIES else None,
            "branches_killed_whole_C33_C34": len(ALT_C44_KILLED),
            "branches_surviving_C33_C34": len(ALT_OPEN),
            "surviving_C33_C34_with_defect0_state_overlay": len(ALT_DEFECT0_KEYS),
            "surviving_C33_C34_unmodelled": len(ALT_UNMODELLED),
            "former_states_modelled": sum(len(s) for s in FAM_SLOTS.values()),
            "former_states_surviving": ALTSW["summary"]["surviving_states"],
            "states_note":
                "former_states_surviving (%d) is the C44 max-plus sweep's "
                "surviving degree-state count across the %d branches that "
                "survived C33/C34, and former_states_modelled (%d) is how many of "
                "them this DAG ever modelled.  The shortfall between them was "
                "GAP-ALT-STATES.  It is RETIRED: those states all live at "
                "a_t >= 11 and a_t <= %d empties every such branch, so the "
                "shortfall is MOOT.  Neither number is a live obligation."
                % (ALTSW["summary"]["surviving_states"], len(ALT_OPEN),
                   sum(len(s) for s in FAM_SLOTS.values()), ALT_A_T_BOUND),
        },
        # column-level emptiness lemmas (frontier stages + e|Phi cell deaths)
        "column_lemmas": {
            "columns": len(COLUMN),
            "from_frontier_stages": n_column_stage,
            "from_divisor_lemma_cells": n_column_divisor,
            "by_level": dict(Counter(r["level"] for r in COLUMN.values())),
            "branches_closed_by_column_lemma": sum(
                1 for n_ in NODES.values() if n_.get("type") == "branch"
                and n_.get("column_lemma_level")),
            "note": "a column lemma empties a whole (a_t, b, branch) column, "
                    "which is this DAG's BRANCH granularity.  Such a branch is "
                    "closed as a LEAF -- exactly like an engine-killed branch -- "
                    "and is NOT gated by the state-enumeration exhaustiveness "
                    "edge, because the lemma is not a statement about the "
                    "enumerated degree-states.  Each affected cell records "
                    "states_with_individual_kill_record vs "
                    "states_covered_only_by_column_lemma so the granularity of "
                    "the evidence is visible and cannot be overstated.",
        },
        "certificates_total": len(CERT_NODES),
        "certificates_found": n_cert_found,
        "certificates_resolved_to_target": n_cert_resolved,
        "unmapped": len(UNMAPPED),
        "engine_killed_branches_independently_audited": {
            w: sum(1 for n_ in nodes_out if n_["type"] == "branch"
                   and n_.get("window") == w
                   and n_.get("cascade_status") not in (None, "survives")
                   and n_.get("level") == "independently-audited")
            for w in ("sub2", "sub1")},
        "engine_killed_branches_claimed": {
            w: sum(1 for n_ in nodes_out if n_["type"] == "branch"
                   and n_.get("window") == w
                   and n_.get("cascade_status") not in (None, "survives")
                   and n_.get("level") == "claimed")
            for w in ("sub2", "sub1")},
    },
    "closure_census": {t: dict(nc) for t, nc in sorted(level_census.items())},
    "unmapped": sorted(UNMAPPED, key=lambda u: (u.get("kind", ""),
                        str(u.get("name") or u.get("key") or u.get("canonical_key") or ""))),
    "nodes": nodes_out,
    "edges": edges_out,
}

with open(os.path.join(HERE, DAG_OUT), "w") as f:
    json.dump(out, f, indent=1, sort_keys=True)

print("=== %s written ===" % DAG_OUT)
print("nodes:", len(nodes_out), "edges:", len(edges_out))
print("nodes by type:", out["counts"]["nodes_by_type"])
print("distinct killed states (phase-D):", out["counts"]["distinct_killed_states"],
      "| defect-0:", out["counts"]["distinct_killed_defect0"],
      "| corner:", out["counts"]["corner_kills"])
print("certificates: total %d found %d resolved %d" % (
    len(CERT_NODES), n_cert_found, n_cert_resolved))
print("cascade-audit JOIN -> engine-killed branches independently-audited:",
      out["counts"]["engine_killed_branches_independently_audited"],
      "| still claimed:", out["counts"]["engine_killed_branches_claimed"])
print("alt-hunt kills ingested:", n_alt_hunt)
print("UNMAPPED bucket:", len(UNMAPPED))
uk = Counter(u.get("kind") for u in UNMAPPED)
for k, v in sorted(uk.items()):
    print("   ", k, v)
print("C0:", c0["closed"], c0["level"], "| subcases closed",
      c0["subcases_closed"], "/", len(c0_subs))
for sid in ("subcase:sub2", "subcase:sub1", "subcase:sub1_alt",
            "subcase:sub1_alt_defect0"):
    s = NODES[sid]
    print("  %-28s closed=%s level=%s branches_open=%d/%d" % (
        sid, s["closed"], s["level"], s.get("branches_open"), s.get("n_branches")))
print("  %-28s closed=%s level=%s (L_D, C10)"
      % ("subcase:dm1", NODES["subcase:dm1"]["closed"], NODES["subcase:dm1"]["level"]))
print("L_alt: %d/%d branches CLOSED%s"
      % (NODES["subcase:sub1_alt"]["branches_closed"], len(ALT_LEAF),
         (" by %s" % out["counts"]["alt_leaf"]["closure_mechanism"])
         if ALT_CLOSURE_APPLIES else ""))
print("       history: %d killed whole by C33+C34, %d survived them "
      "(%d with a defect-0 overlay of %d slots, %d unmodelled); the %d surviving "
      "degree-states behind GAP-ALT-STATES are MOOT"
      % (len(ALT_C44_KILLED), len(ALT_OPEN), len(ALT_DEFECT0_KEYS),
         sum(len(s) for s in FAM_SLOTS.values()), len(ALT_UNMODELLED),
         ALTSW["summary"]["surviving_states"]))
print("column lemmas: %d columns (%d stage attributions, %d e|Phi cell deaths) "
      "-> %d branches closed as leaves; by level %s"
      % (out["counts"]["column_lemmas"]["columns"], n_column_stage,
         n_column_divisor,
         out["counts"]["column_lemmas"]["branches_closed_by_column_lemma"],
         out["counts"]["column_lemmas"]["by_level"]))
