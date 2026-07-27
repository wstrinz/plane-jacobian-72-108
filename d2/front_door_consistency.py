#!/usr/bin/env python3
"""front_door_consistency.py  (NEW 2026-07-27)

THE DEFECT THIS EXISTS TO MAKE IMPOSSIBLE.  On 2026-07-27 the public tree's
README.md still said

    "The (72,108) case is **not closed** ... the target theorem `C0` remains
     **open**: the frontier above is nonempty"

while `proof_dag.json` in the same commit recorded `C0: closed=true,
subcases_closed=5` and `PROOF_72_108.md` sat unreferenced in `d2/`.  A
mathematician landing on the repository learned the OPPOSITE of the result.  The
private tree's root README said the same thing.  Neither was caught by anything,
because 114 checkers gated the mathematics and zero gated the front door.

WHY A CHECKER RATHER THAN A CAREFUL EDIT.  Editing the prose fixes today.  The
repository's own standing trap is that "prose and code drift, and a supersession
banner does NOT stop it -- only editing the body works", and the durable form of
that lesson is "where possible make drift a TEST FAILURE".  This is that test.

WHAT IT CHECKS, AND WHY BIDIRECTIONALLY.  Every fix in this repository that stuck
was a CROSS-CHECK between two halves, not another assertion about one.  So this
does not merely assert "the README says closed".  It reads the registry, reads
the prose, and requires them to AGREE -- in whichever direction the registry
happens to point.  If `C0` were ever reopened, a README still claiming closure
would fail here just as loudly as today's defect does.

Detection is by PROXIMITY, not by a blocklist of exact phrasings: a blocklist
only ever catches the sentence you already fixed.  We find sentences that name
the target case, then flag those carrying an openness marker.  That is why
"27 open" (the wider 34-row frontier) and "JC(2) is a separate open problem"
do not trip it -- neither sentence names the target case.

Checker: --quiet, exit 0 iff every check passes.  <1 s.  Reads only.
"""
from __future__ import annotations

import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__)) or "."
QUIET = "--quiet" in sys.argv
_ok = [0]
_fail: list[str] = []
_notes: list[str] = []


def ck(name: str, cond: bool, detail: str = "") -> bool:
    if cond:
        _ok[0] += 1
        if not QUIET:
            print("[OK]   %s" % name)
    else:
        _fail.append(name)
        print("[FAIL] %s%s" % (name, ("  -- " + detail) if detail else ""))
    return bool(cond)


def note(msg: str) -> None:
    """A stated non-check.  Never counted as a pass."""
    _notes.append(msg)
    print("[NOTE] %s" % msg)


# --------------------------------------------------------------------------
# 1.  The registry half.
# --------------------------------------------------------------------------
DAG = json.load(open(os.path.join(HERE, "proof_dag.json"), encoding="utf-8"))
C0 = next((n for n in DAG["nodes"] if n.get("id") == "C0"), None)

ck("R1  proof_dag.json contains the target node C0", C0 is not None)
if C0 is None:
    raise SystemExit(1)

C0_CLOSED = bool(C0.get("closed"))
C0_LEVEL = C0.get("level")
ck("R2  C0 records an explicit closed flag and an evidence level (closed=%s, "
   "level=%s)" % (C0_CLOSED, C0_LEVEL),
   "closed" in C0 and isinstance(C0_LEVEL, str))

# The level cap is deliberate and documented; assert it is STILL the documented
# one, so a silent regrade cannot slip past while the prose keeps saying
# "claimed".  See PROOF_72_108.md sec.0.5 -- the only routes above `claimed` are a
# machine-checkable reformulation of the partition, or a formal proof.
ck("R3  C0's recorded level is one the front-door prose is allowed to describe "
   "(claimed | exact-checked | proved)",
   C0_LEVEL in ("claimed", "exact-checked", "proved"), str(C0_LEVEL))


# --------------------------------------------------------------------------
# 2.  The prose half.  Locate the front-door documents from the repo root.
# --------------------------------------------------------------------------
def repo_root(start: str) -> str:
    d = start
    for _ in range(4):
        if os.path.isdir(os.path.join(d, ".git")):
            return d
        d = os.path.dirname(d)
    return os.path.dirname(start)


ROOT = repo_root(HERE)
FRONT_DOORS = [p for p in (os.path.join(ROOT, "README.md"),
                           os.path.join(ROOT, "docs", "README.md"))
               if os.path.exists(p)]

ck("P1  at least one front-door document was located under the repo root (%s)"
   % ROOT, bool(FRONT_DOORS), "none found")

# Sentences that NAME the target case.  Anything not naming it is out of scope --
# the repository legitimately calls JC(2) open and legitimately reports 27 open
# rows in the wider 34-row atlas.
NAMES_TARGET = re.compile(r"\(72,\s*108\)|\bC0\b|\(8,\s*28\)", re.I)

# Markers asserting the thing named is NOT settled.
OPEN_MARKER = re.compile(
    r"not\s+closed|remains?\s+\*{0,2}open|still\s+open|is\s+\*{0,2}open\*{0,2}\b"
    r"|left\s+open|in\s+progress|unresolved|not\s+yet\s+closed", re.I)

# Markers asserting the thing named IS settled.
CLOSED_MARKER = re.compile(
    r"\bexcluded\b|\bclosed\b|\bempty\b|no\s+admissible\s+germ"
    r"|there\s+is\s+no\s+pair", re.I)

# Sentences that are NOT a present claim about the target's status: either they
# quote the history, or they are a directive to the reader about how to read
# some other document.  The second class matters because the honest way to flag
# a stale artifact is a sentence like "do not read FRONTIER.md as what is left
# open for (72,108)" -- which names the target and contains an openness phrase
# while asserting the exact opposite of openness.
NOT_A_PRESENT_CLAIM = re.compile(
    r"left\s+open\s+in\s+2022|was\s+the\s+GGV|GGHV22\s+state|used\s+to|"
    r"before\s+the\s+proof|at\s+one\s+time|earlier\s+draft|formerly|"
    r"prior\s+to|it\s+is\s+left\s+open|"
    r"do\s+not\s+read|should\s+not\s+be\s+read|do\s+not\s+treat|"
    r"is\s+not\s+a\s+claim|rather\s+than\s+left\s+to\s+be\s+discovered", re.I)

SENT = re.compile(r"[^.!?]*[.!?]")


def sentences(text: str):
    """Yield sentences, joining Markdown's hard line wrapping first.

    The first version of this split on newlines as well as terminators.  These
    documents are hard-wrapped at ~79 columns, so nearly every real sentence
    spans a line break and was silently never examined -- P2 reported green over
    a README that said "The (72,108) case is **not closed**".  A detector that
    cannot see the defect it was written for is worse than no detector, because
    it reads as evidence of absence.
    """
    # Fenced code blocks are transcripts, not claims.
    text = re.sub(r"```.*?```", " ", text, flags=re.S)
    for para in re.split(r"\n\s*\n", text):
        para = " ".join(para.split())          # unwrap
        if not para:
            continue
        for m in SENT.finditer(para):
            s = m.group(0).strip()
            if s:
                yield s
        tail = SENT.sub("", para).strip()      # unterminated final clause
        if tail:
            yield tail


violations = []
supporting = []
for path in FRONT_DOORS:
    text = open(path, encoding="utf-8").read()
    rel = os.path.relpath(path, ROOT)
    for s in sentences(text):
        if not NAMES_TARGET.search(s):
            continue
        if NOT_A_PRESENT_CLAIM.search(s):
            continue
        says_open = bool(OPEN_MARKER.search(s))
        # A negated settledness claim ("is **not closed**") contains the word
        # `closed`, so CLOSED_MARKER must be evaluated on the sentence with
        # negations neutralised first -- otherwise the openness assertion
        # SUPPRESSES its own violation.  This checker's first run did exactly
        # that on the private README and reported P2 green.
        neutral = re.sub(r"\b(?:not|never|isn't|is\s+not|nor)\s+"
                         r"\**(closed|excluded|empty|settled)\**",
                         " __OPENCLAIM__ ", s, flags=re.I)
        says_closed = bool(CLOSED_MARKER.search(neutral))
        if C0_CLOSED and says_open:
            violations.append((rel, s))
        elif (not C0_CLOSED) and says_closed:
            violations.append((rel, s))
        elif C0_CLOSED and says_closed:
            supporting.append((rel, s))

ck("P2  no front-door sentence naming the target case contradicts the registry "
   "(registry says closed=%s; %d contradicting sentence(s))"
   % (C0_CLOSED, len(violations)),
   not violations,
   " || ".join("%s: %s" % (r, s[:110]) for r, s in violations[:4]))

ck("P3  and the front door makes the registry's verdict POSITIVELY, rather than "
   "merely not contradicting it (%d supporting sentence(s))" % len(supporting),
   bool(supporting) if C0_CLOSED else True,
   "no sentence states the result")


# --------------------------------------------------------------------------
# 3.  The front door must route a reader to the proof and to the priority claim.
# --------------------------------------------------------------------------
ALL_FRONT = "\n".join(open(p, encoding="utf-8").read() for p in FRONT_DOORS)

ck("P4  the front door links PROOF_72_108.md -- the proof must be reachable "
   "without knowing its filename", "PROOF_72_108.md" in ALL_FRONT)

ck("P5  the front door records Helali's priority BY DOI, not merely by name "
   "(10.5281/zenodo.21479814)",
   "10.5281/zenodo.21479814" in ALL_FRONT)

ck("P6  ... and says in so many words that he was first",
   re.search(r"Helali\s+was\s+first|got\s+there\s+first|first\s+published"
             r"|the\s+first\s+published\s+exclusion", ALL_FRONT, re.I) is not None)

# The two conditions are the whole reason Corollary B is not a bare theorem.
ck("P7  the front door names BOTH standing conditions (Prop 4.3 exhaustiveness, "
   "and the alpha-strip WLOG [QQ1])",
   re.search(r"Prop(osition)?\.?\s*4\.3", ALL_FRONT) is not None
   and "QQ1" in ALL_FRONT)

# Backticked, because `claimed` as a bare word occurs in ordinary prose
# ("What is NOT claimed") and matching that would be a free pass.
ck("P8  the front door states the registry's recorded level as a code-quoted "
   "term, so a reader is not left to infer it (`%s`)" % C0_LEVEL,
   ("`%s`" % C0_LEVEL) in ALL_FRONT,
   "`%s` appears nowhere in the front door" % C0_LEVEL)


# --------------------------------------------------------------------------
# 4.  Cross-tree honesty: say which front doors were actually inspected.
# --------------------------------------------------------------------------
if len(FRONT_DOORS) < 2:
    note("only %d front-door document(s) found under %s (%s). This checker gates "
         "whatever front doors the tree HAS; a tree with fewer has fewer gated."
         % (len(FRONT_DOORS), ROOT,
            ", ".join(os.path.relpath(p, ROOT) for p in FRONT_DOORS)))

if _fail:
    print()
    print("FAILURES (%d):" % len(_fail))
    for f in _fail:
        print("   - %s" % f)
    raise SystemExit(1)

print("front_door_consistency: %d/%d checks pass (registry C0 closed=%s level=%s; "
      "front doors: %s)"
      % (_ok[0], _ok[0], C0_CLOSED, C0_LEVEL,
         ", ".join(os.path.relpath(p, ROOT) for p in FRONT_DOORS)))
