#!/usr/bin/env python3
"""Consistency checks for the committed cone-certificate ledger.

1. Structural coverage: every certified kill has an L or B certificate in
   EVERY zero-flag case; survivors have at least one open case.
2. Verdict agreement with cascade_cones.json on all 420 open branches.
3. Sample regeneration: re-derive certificates for a deterministic sample
   of branches directly from the engine tables and compare.
"""

from __future__ import annotations

import json
from pathlib import Path

import cone_lemmas as C

def _require(_cond, _msg):
    """Proof-critical check: fails loudly and exits nonzero, unaffected by python -O."""
    if not _cond:
        import sys as _sys
        print("FAIL: " + str(_msg))
        _sys.exit(1)


ROOT = Path(__file__).resolve().parent


def main() -> None:
    certs = json.loads(
        (ROOT / "cascade_cone_certificates.json").read_text(encoding="utf-8")
    )
    cones = json.loads(
        (ROOT / "cascade_cones.json").read_text(encoding="utf-8")
    )
    engine_status = {
        (r["a_t"], tuple(r["b"]), r["branch"]): r["status"]
        for r in cones["branches"]
    }

    _require(len(certs["branches"]) == 420, "len(certs[\"branches\"]) == 420")
    kills = 0
    for record in certs["branches"]:
        key = (record["a_t"], tuple(record["b"]), record["branch"])
        open_cases = [c for c in record["cases"] if c["kind"] == "open"]
        if record["certified_kill"]:
            kills += 1
            _require(not open_cases, key)
            _require(engine_status[key] == "engine_killed_pending_audit", key)
        else:
            _require(open_cases, key)
            _require(engine_status[key] == "survives", key)
    _require(kills == 390, kills)
    _require(certs["summary"]["kills_certified"] == 390, "certs[\"summary\"][\"kills_certified\"] == 390")
    print("  ledger structure and engine agreement: 390 kills / 30 survivors")

    tables = C.Tables()
    sample = certs["branches"][::40]  # deterministic spread, ~11 records
    for record in sample:
        regenerated = C.certify_branch(
            tables, record["a_t"], tuple(record["b"]), record["branch"]
        )
        stored = record["cases"]
        _require(len(regenerated) == len(stored), "len(regenerated) == len(stored)")
        for fresh, kept in zip(regenerated, stored):
            if fresh is None:
                _require(kept["kind"] == "open", "kept[\"kind\"] == \"open\"")
            else:
                _require(fresh["kind"] == kept["kind"], "fresh[\"kind\"] == kept[\"kind\"]")
                if fresh["kind"] == "B":
                    _require(fresh["dimension"] == kept["dimension"], "fresh[\"dimension\"] == kept[\"dimension\"]")
                    _require(fresh["sum_of_minima"] == kept["sum_of_minima"], "fresh[\"sum_of_minima\"] == kept[\"sum_of_minima\"]")
    print(f"  regenerated {len(sample)} branch certificates match the ledger")

    unconditional = certs["summary"]["unconditional_local_kills"]
    live_betas = {0, 1, 3}
    for entry in unconditional:
        _require(entry["beta"] not in live_betas or (
            entry["branch"] == "T2" and entry["a_t"] == 0 and entry["beta"] == 3
        ), entry)
    print("  unconditional local kills avoid the live values {0,1,3}"
          " (single exception T2 a=0 beta=3)")
    print("cone lemmas: PASS")


if __name__ == "__main__":
    main()
