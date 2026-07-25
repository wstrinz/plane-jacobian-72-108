#!/usr/bin/env python3
"""case_envelope.py  --  a thin provenance envelope for algebraic-system experiments

NOT an execution engine, NOT a case IR, NOT an abstraction layer.  It is the
smallest shared vocabulary that lets independent experiment lanes talk about the
same object, plus enough provenance to tell whether a result depended on a
particular backend version.

WHY IT EXISTS.  Four experiment lanes (jet-obstruction, monomial-lens, chordal
profiling, proof-pattern mining) each emit their own JSON.  Without a shared
envelope they invent four incompatible notions of "system", and then:
  * Monomial Lens cannot hand a proposed transform to the chordal lane;
  * the chordal lane cannot profile before/after that transform;
  * the obstruction lane cannot consume a transformed system;
  * the proof-DAG cannot ingest any of it later.

DESIGN RULE -- DESCRIPTIVE, NOT PRESCRIPTIVE.  `wrap()` takes an artifact a lane
has ALREADY produced and attaches an envelope.  Nothing here can block a lane or
force it to restructure its own output; a lane's payload is carried verbatim
under "payload".  If this module ever becomes a thing an experiment must satisfy
before it can run, that is a bug in how it is being used.

WHY BACKEND VERSIONS ARE FIRST-CLASS.  Local Macaulay2 is 1.19.1 with
VersalDeformations 3.0, while 4.0 (Jul 2025) changed the lifting/nested-
deformation functionality; Singular is 4.2.1; sympy 1.14.  A result that depends
on a backend version is a different result, and today's session already produced
two silent-model-drift bugs (the G5 normalisation, and a chart exponent confused
with a ramification index).  Recording the version is the cheapest possible guard
against a third.

RESULT KINDS -- the epistemic status of an answer, deliberately NOT collapsed to
pass/fail:
    PROOF               exact, self-contained, checkable
    COUNTEREXAMPLE      an explicit witness
    EXACT_CONSEQUENCE   exact, but a necessary condition / contraction only
    NUMERICAL_EVIDENCE  homotopy / floating point
    MODULAR_EVIDENCE    holds mod p, not lifted to Q
    COST                the computation did not finish (NOT a mathematical verdict)
    MODEL_GAP           the model does not faithfully represent the intended object

`COST` and `MODEL_GAP` are the two that matter most: a timeout is not emptiness,
and a wrong model is not a wrong theorem.

Read-only over the repo.  Usage:
    from case_envelope import wrap, backend_versions
    env = wrap(payload, system_id="j6-depth-3", question="lift_obstruction",
               result_kind="EXACT_CONSEQUENCE", equations=gens)
    json.dump(env, fh, indent=2, sort_keys=True)

    python case_envelope.py --versions     # print captured backend versions
    python case_envelope.py --validate F   # validate an envelope file
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys

SCHEMA = "d2-case-envelope-v1"

RESULT_KINDS = ("PROOF", "COUNTEREXAMPLE", "EXACT_CONSEQUENCE",
                "NUMERICAL_EVIDENCE", "MODULAR_EVIDENCE", "COST", "MODEL_GAP")

REQUIRED = ("schema", "system_id", "question", "result_kind", "backend", "scope")


def _sha256(text):
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def equations_digest(equations):
    """Order-independent digest of a system's equations.

    Sorting the per-equation digests makes the value independent of generator
    ORDER, so the same system presented in a different order digests the same --
    which is what you want when two lanes build it independently.
    """
    if equations is None:
        return None
    parts = sorted(_sha256(str(e)) for e in equations)
    return _sha256("|".join(parts))


def _run(argv, timeout=25):
    """Run a fixed argv (no shell). Returns output lines, [] on any failure."""
    try:
        out = subprocess.run(argv, capture_output=True, text=True, timeout=timeout)
        return (out.stdout or out.stderr or "").strip().splitlines()
    except Exception:
        return []


def backend_versions(probe=True):
    """Capture versions of every backend that could change an answer.

    `probe=False` skips the external calls (they cost a few seconds) and reports
    only the in-process ones.
    """
    info = {"python": sys.version.split()[0]}
    try:
        import sympy
        info["sympy"] = sympy.__version__
    except Exception:
        pass
    if not probe:
        return info

    lines = _run('wsl.exe -e bash -lc "M2 --version" 2>/dev/null')
    if lines:
        info["macaulay2"] = lines[0].strip()
    lines = _run('wsl.exe -e bash -lc "Singular --version" 2>/dev/null')
    if lines:
        info["singular"] = lines[0].strip()[:80]
    return info


# Package versions are NOT auto-probed (each costs an M2 start-up).  A lane that
# uses a package should pass its version explicitly -- e.g.
#   wrap(..., backend_packages={"VersalDeformations": "3.0", "Chordal": "0.2"})
# The known-local values, recorded 2026-07-24, are:
KNOWN_M2_PACKAGES = {
    "Jets": "1.0",
    "Chordal": "0.2",
    "NoetherianOperators": "2.2.1",
    "VersalDeformations": "3.0",      # upstream 4.0 (Jul 2025) differs -- see docstring
    "Saturation": "0.2",
    "NumericalAlgebraicGeometry": "1.17",
}


def wrap(payload, system_id, question, result_kind, *,
         parent_system_id=None, equations=None, variables=None, field="QQ",
         open_conditions=None, gradings=None, transforms=None,
         backend_name="python/sympy", backend_packages=None,
         artifacts=None, scope=None, probe_versions=True):
    """Attach a provenance envelope to a payload a lane has already produced."""
    if result_kind not in RESULT_KINDS:
        raise ValueError("result_kind must be one of %s, got %r"
                         % (", ".join(RESULT_KINDS), result_kind))
    env = {
        "schema": SCHEMA,
        "system_id": system_id,
        "parent_system_id": parent_system_id,
        "ring": {"field": field, "variables": [str(v) for v in (variables or [])]},
        "equations_sha256": equations_digest(equations),
        "n_equations": (len(equations) if equations is not None else None),
        "open_conditions": [str(c) for c in (open_conditions or [])],
        "gradings": gradings or {},
        "transforms": transforms or [],
        "backend": {
            "name": backend_name,
            "versions": backend_versions(probe=probe_versions),
            "packages": backend_packages or {},
        },
        "question": question,
        "result_kind": result_kind,
        "artifacts": artifacts or [],
        "scope": scope or {},
        "payload": payload,
    }
    return env


def validate(env):
    """Return a list of problems; empty means well-formed."""
    problems = []
    for key in REQUIRED:
        if key not in env:
            problems.append("missing required field: %s" % key)
    if env.get("schema") != SCHEMA:
        problems.append("schema is %r, expected %r" % (env.get("schema"), SCHEMA))
    if env.get("result_kind") not in RESULT_KINDS:
        problems.append("result_kind %r not in %s"
                        % (env.get("result_kind"), ", ".join(RESULT_KINDS)))
    if not env.get("scope"):
        problems.append("scope is empty -- every envelope must state what it does "
                        "NOT establish")
    backend = env.get("backend") or {}
    if not backend.get("versions"):
        problems.append("backend.versions is empty -- a result that depends on a "
                        "backend version is a different result")
    return problems


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--versions", action="store_true")
    ap.add_argument("--validate", metavar="FILE")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()

    if args.versions:
        print(json.dumps({"backends": backend_versions(),
                          "macaulay2_packages": KNOWN_M2_PACKAGES},
                         indent=2, sort_keys=True))
        return 0

    if args.validate:
        env = json.load(open(args.validate, encoding="utf-8"))
        problems = validate(env)
        if problems:
            print("INVALID (%d problems):" % len(problems))
            for p in problems:
                print("  -", p)
            return 1
        print("VALID envelope: %s / %s / %s"
              % (env["system_id"], env["question"], env["result_kind"]))
        return 0

    # default: self-test
    env = wrap({"demo": True}, system_id="selftest", question="smoke",
               result_kind="MODULAR_EVIDENCE",
               equations=["x^2 - 1", "y - x"], variables=["x", "y"],
               scope={"establishes": "nothing; this is a smoke test"},
               probe_versions=False)
    problems = validate(env)
    # order-independence of the digest is the one behaviour worth asserting
    d1 = equations_digest(["a", "b", "c"])
    d2 = equations_digest(["c", "a", "b"])
    if d1 != d2:
        problems.append("equations digest is order-DEPENDENT; it must not be")
    if equations_digest(["a", "b"]) == d1:
        problems.append("equations digest collides across different systems")
    if problems:
        print("SELFTEST FAILED:")
        for p in problems:
            print("  -", p)
        return 1
    print("case_envelope selftest OK (schema %s, %d result kinds)"
          % (SCHEMA, len(RESULT_KINDS)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
