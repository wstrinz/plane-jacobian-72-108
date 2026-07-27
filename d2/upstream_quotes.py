#!/usr/bin/env python3
"""Self-contained access to the upstream arXiv source strings the checkers test.

WHY THIS EXISTS.  Two checkers -- ``moh_discards.py`` and ``moh_control_50_75.py``
-- assert that specific literal strings appear in GGV5 (arXiv:1708.07936) and
GGV3 (arXiv:1406.0886).  Those ``.tex`` files are other authors' copyrighted
sources and are deliberately NOT redistributed in the public release
(``PUBLICATION_AUDIT.md`` item 2).  Reading them live therefore made both
checkers unrunnable on a clean public clone: ``moh_discards`` failed its own
B2-B6 and ``moh_control_50_75`` died with ``FileNotFoundError``.

THE CONTRACT, which is the same one ``paper_src/upstream_facts.json`` already
uses for GGHV22:

  * ``paper_src/upstream_quotes.json`` records each probe ONCE -- the literal
    string, its arXiv id, the 1-based line number it was transcribed from, and
    the sha256 of the ``.tex`` it was transcribed from.
  * ``present(key)`` answers from that JSON, so the check runs everywhere.
  * When the ``.tex`` IS on disk, ``verify_against_tex()`` re-derives every probe
    from the source.  A probe recorded present that is NOT in the source is a
    **FAILURE**, not a skip -- the whole point is that the transcription cannot
    silently drift away from the paper.

THE TRAP THIS IS SHAPED AROUND.  A transcription is a second copy of the source,
and "two objects sharing a name is where the errors live".  The mitigation is
that the two halves are cross-checked whenever both are available, rather than
each being separately asserted -- and that the sha256 pins WHICH version of the
source the transcription is true of.
"""
from __future__ import annotations

import hashlib
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__)) or "."
_JSON = os.path.join(HERE, "paper_src", "upstream_quotes.json")

with open(_JSON, encoding="utf-8") as _fh:
    DATA = json.load(_fh)

PROBES = DATA["probes"]
SOURCES = DATA["sources"]
COUNTS = DATA["counts"]


def present(key: str) -> bool:
    """True iff the probe was recorded present in the upstream source."""
    if key not in PROBES:
        raise KeyError("no such upstream probe: %r" % key)
    return PROBES[key]["occurrences"] > 0


def count(key: str) -> int:
    """A recorded numeric count read off the upstream source."""
    return COUNTS[key]["value"]


def cite(key: str) -> str:
    """Human-readable provenance for a probe, e.g. 'GGV3 tex:1720'."""
    p = PROBES[key]
    return "%s tex:%d" % (p["source"], p["tex_line"])


def tex_path(source: str) -> str:
    # local_tex is repo-relative and already includes the paper_src/ segment.
    return os.path.join(HERE, *SOURCES[source]["local_tex"].split("/"))


def available(source: str) -> bool:
    return os.path.exists(tex_path(source))


def verify_against_tex():
    """Re-derive every probe from the local .tex files, where present.

    Returns ``(results, checked_sources)`` where ``results`` is a list of
    ``(key, ok, detail)``.  Sources absent from disk contribute nothing -- they
    are not failures, because the public release ships without them by design.
    """
    results = []
    checked = []
    for source in sorted(SOURCES):
        path = tex_path(source)
        if not os.path.exists(path):
            continue
        checked.append(source)
        raw = open(path, encoding="utf-8", errors="replace").read()
        sha = hashlib.sha256(open(path, "rb").read()).hexdigest()
        recorded = SOURCES[source]["local_tex_sha256"]
        results.append((
            "%s sha256 matches the transcription's source" % source,
            sha == recorded,
            "on disk %s..., recorded %s..." % (sha[:16], recorded[:16]),
        ))
        for key, p in sorted(PROBES.items()):
            if p["source"] != source:
                continue
            hay = raw.replace(" ", "") if p["normalisation"] == "nospace" else raw
            needle = (p["probe"].replace(" ", "")
                      if p["normalisation"] == "nospace" else p["probe"])
            n = hay.count(needle)
            results.append((
                "%s (%s) occurs in source" % (key, cite(key)),
                (n > 0) == (p["occurrences"] > 0),
                "source n=%d, recorded n=%d" % (n, p["occurrences"]),
            ))
        for key, c in sorted(COUNTS.items()):
            if c["source"] != source:
                continue
            n = raw.count(r"{\color{red}")
            results.append((
                "%s count matches source" % key,
                n == c["value"],
                "source n=%d, recorded n=%d" % (n, c["value"]),
            ))
    return results, checked


if __name__ == "__main__":
    import sys

    quiet = "--quiet" in sys.argv
    res, checked = verify_against_tex()
    if not checked:
        print("upstream_quotes: %d probes recorded; no local .tex present, "
              "nothing to re-verify (this is the public-clone case)"
              % len(PROBES))
        raise SystemExit(0)
    bad = [r for r in res if not r[1]]
    if not quiet:
        for name, ok_, detail in res:
            print("  [%s] %s -- %s" % ("OK" if ok_ else "FAIL", name, detail))
    print("upstream_quotes: %d/%d re-verified against %s"
          % (len(res) - len(bad), len(res), ", ".join(checked)))
    if bad:
        for name, _, detail in bad:
            print("FAIL: %s -- %s" % (name, detail))
        raise SystemExit(1)
