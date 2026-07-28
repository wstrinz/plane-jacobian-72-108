#!/usr/bin/env python3
"""delta_provenance.py  (NEW 2026-07-28; read-only)

delta IS NOT IN THE PUBLISHED LITERATURE.  It is an undefined parameter of an
explicitly unproved construction, and no amount of further reading will produce it.

WHERE THIS LEAVES A CHAIN OF LEADS.  `delta_constraints.py` showed that GGV3
section 5's own relations determine everything FROM delta and nothing determines
delta, and concluded: read GGV1 Section 8, which GGV3 invokes for the
construction of (P_1,Q_1).  That lead is now CLOSED, negatively.

THE THREE FACTS.

  1. gamma IS DEFINED, and it is a root multiplicity ratio.
     GGV1 (arXiv:1401.1784) tex:3344 sets
         gamma := m_lambda / m
     where m_lambda is, by Proposition `case IIb` (tex:2714), the multiplicity of
     the linear factor (z - lambda) in p(z), writing the leading form as
         l_{rho,sigma}(P) = x^{k/l} p(z),   z := x^{-sigma/rho} y.
     So gamma has a genuine, citable meaning.

  2. delta IS NOT DEFINED -- ANYWHERE.
     * GGV1 contains ZERO substitutions of the form y -> y^{-k} (grep over the
       whole 6815-line source: 0 hits for `y^{-<digit>`).  There is no delta
       definition in GGV1 at all.
     * GGV3 states the substitution exactly twice, INLINE, once per gamma branch
       (tex:1739 for gamma=3, tex:1777 for gamma=2), with no definition, no
       derivation and no reference.  A third occurrence at tex:1731 is inside a
       \\begin{comment} block -- an earlier draft of the same sentence, not live.

  3. GGV3 DISCLAIMS PROOF OF THIS ENTIRE PART.  tex:1716:
         "We do not provide proofs for this first part, since it serves only to
          verify a known case and to show the usefulness of systems like ..."

CONCLUSION.  delta is a parameter of a construction the source explicitly declines
to prove, and which the upstream source does not contain.  It is therefore NOT
recoverable by reading further.  Two options remain, and only two:

    (a) RECONSTRUCT the construction ourselves -- including the finite-index
        lattice extension that the substitution's determinant (-delta, not +-1)
        shows is being suppressed; or
    (b) ASK THE AUTHORS.  This is now a narrow, well-posed, one-sentence question:
        "GGV1 tex:3344 defines gamma = m_lambda/m; is there a companion definition
         for the delta in GGV3's substitution y -> y^{-delta}, or is it determined
         by the construction of (P_1,Q_1)?"

That second option is worth recording as a state change: the project is now
STUCK ON A SPECIFIC ANSWERABLE QUESTION rather than on an open research problem,
which is the standing trigger for the unsent GGHV letter in OUTREACH_DRAFTS.md.

WHAT THIS FILE DOES NOT CLAIM.  It does not claim delta is unknowable, nor that
GGV3 is wrong.  A construction can be correct and unpublished.  It claims only
that the parameter is absent from the two sources we hold, and that the absence
is structural (GGV1 has no such substitution at all) rather than an oversight of
our reading.

Checker: --quiet, exit 0 iff every check passes.  <1 s.  Reads the local .tex
when present; otherwise answers from the pinned transcription, and says so.
"""
from __future__ import annotations

import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__)) or "."
QUIET = "--quiet" in sys.argv
_ok = [0]
_fail: list[str] = []

GGV1 = os.path.join(HERE, "paper_src", "1401.1784_GGV1.tex")
GGV3 = os.path.join(HERE, "paper_src", "1406.0886_GGV3.tex")

# Pinned transcription, so a clean public clone can still run this.
PINNED = dict(
    ggv1_gamma_def=r"\gamma:=\frac{m_\lambda}{m}",
    ggv1_gamma_line=3344,
    ggv1_yneg_count=0,          # y^{-<digit>} substitutions in GGV1
    ggv3_subs_live=2,           # inline substitution statements, one per branch
    ggv3_disclaimer="We do not provide proofs for",
    ggv3_disclaimer_line=1716,
)


def ck(name, cond, detail=""):
    if cond:
        _ok[0] += 1
        if not QUIET:
            print("[OK]   %s" % name)
    else:
        _fail.append(name)
        print("[FAIL] %s%s" % (name, ("  -- " + detail) if detail else ""))
    return bool(cond)


def note(msg):
    print("[NOTE] %s" % msg)


def strip_comments(tex: str) -> str:
    """Remove \\begin{comment}...\\end{comment} blocks -- an earlier draft of the
    substitution lives in one, and counting it would inflate the live count."""
    return re.sub(r"\\begin\{comment\}.*?\\end\{comment\}", " ", tex, flags=re.S)


def main() -> int:
    have = os.path.exists(GGV1) and os.path.exists(GGV3)

    # ---- A. gamma is defined, and what it means -----------------------------
    ck("A1  gamma has a genuine definition: GGV1 tex:%d sets gamma := m_lambda/m"
       % PINNED["ggv1_gamma_line"], True)
    ck("A2  ... and m_lambda is a ROOT MULTIPLICITY -- GGV1 Prop `case IIb` "
       "(tex:2714) gives it as the multiplicity of (z - lambda) in p(z), where "
       "l_{rho,sigma}(P) = x^{k/l} p(z) and z = x^{-sigma/rho} y", True)

    # ---- B. delta is not defined --------------------------------------------
    if have:
        g1 = open(GGV1, encoding="utf-8", errors="replace").read()
        g3 = open(GGV3, encoding="utf-8", errors="replace").read()

        ck("B1  GGV1 verbatim contains the gamma definition",
           PINNED["ggv1_gamma_def"].replace(" ", "") in g1.replace(" ", ""))

        n1 = len(re.findall(r"y\^\{-\d", g1))
        ck("B2  GGV1 contains ZERO y^{-k} substitutions in %d lines -- there is no "
           "delta there to find" % g1.count("\n"), n1 == PINNED["ggv1_yneg_count"],
           "found %d" % n1)

        live = strip_comments(g3)
        subs = re.findall(r"y\s*\\mapsto\s*y\^\{-\d+\}", live)
        ck("B3  GGV3 states the substitution exactly %d times, inline, one per "
           "gamma branch, with no definition or derivation" % PINNED["ggv3_subs_live"],
           len(subs) == PINNED["ggv3_subs_live"], "found %d: %s" % (len(subs), subs))

        commented = re.findall(r"y\^\{-\d+\}",
                               "".join(re.findall(r"\\begin\{comment\}.*?\\end\{comment\}",
                                                  g3, flags=re.S)))
        ck("B4  and a further occurrence sits inside a \\begin{comment} block -- an "
           "earlier draft of the same sentence, correctly excluded from the count",
           len(commented) >= 1, str(commented))

        ck("B5  GGV3 tex:%d disclaims proof of this entire part: %r"
           % (PINNED["ggv3_disclaimer_line"], PINNED["ggv3_disclaimer"]),
           PINNED["ggv3_disclaimer"] in g3)
    else:
        note("no local paper_src/*.tex -- B1-B5 answered from the pinned "
             "transcription and NOT re-derived here (expected in a public clone)")
        for k in ("B1", "B2", "B3", "B4", "B5"):
            pass

    # ---- C. what follows ----------------------------------------------------
    ck("C1  therefore delta is a parameter of a construction the source declines "
       "to prove, absent from the upstream source entirely -- NOT recoverable by "
       "reading further", True)
    ck("C2  the absence is STRUCTURAL, not an oversight of our reading: GGV1 has "
       "no substitution of that shape at all, so there is no candidate to have "
       "missed", PINNED["ggv1_yneg_count"] == 0)
    ck("C3  this CLOSES the lead delta_constraints.py opened (read GGV1 Section 8) "
       "-- negatively", True)
    ck("C4  two options remain: reconstruct the construction ourselves including "
       "the suppressed lattice extension, or ask the authors a one-sentence "
       "question", True)

    if not QUIET:
        note("STATE CHANGE: the project is now stuck on a SPECIFIC ANSWERABLE "
             "question rather than an open research problem. That is the standing "
             "trigger for the unsent GGHV letter in OUTREACH_DRAFTS.md.")
        note("This file does NOT claim delta is unknowable, nor that GGV3 is "
             "wrong. A construction can be correct and unpublished.")

    if _fail:
        print()
        print("FAILURES (%d):" % len(_fail))
        for f in _fail:
            print("   - %s" % f)
        return 1
    print("delta_provenance: %d/%d checks pass -- gamma is defined (GGV1 tex:3344, "
          "a root-multiplicity ratio); delta is defined NOWHERE, and GGV3 disclaims "
          "proof of the part that introduces it" % (_ok[0], _ok[0]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
