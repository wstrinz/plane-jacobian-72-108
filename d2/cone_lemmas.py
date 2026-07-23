#!/usr/bin/env python3
"""Cone-certificate extraction for the cascade-engine Phase B kills.

Compresses the 390 engine kills (cascade_cones.json, depth 4) into two
reusable lemma families and emits a row-by-row certificate ledger:

  L (local place kill): at a single split place with v_p(e) = beta, the
    chain of level identities 7/6 -> 4 admits NO consistent local valuation
    chain under the given zero flags.  This is independent of the other
    three places, so it kills every stratum whose multiplicity vector
    contains beta (for that a, branch, flag case).

  B (budget kill): every place admits local chains, but one single budget
    dimension already overflows: the sum over the four places of the
    minimum possible consumption exceeds the global degree cap.

A branch record is certificate-killed iff EVERY zero-flag case receives an
L or B certificate.  The script verifies that the certificate verdict
matches the engine verdict on all 420 open branches (kills AND survivors)
and refuses to emit otherwise.

Scope note: the tables are generated from the engine's own per-place
enumeration (cascade_engine.place_profiles), so this artifact provides
compression and reuse, not independent verification — that is the role of
the separately written audit checker.
"""

from __future__ import annotations

import json
from pathlib import Path

import cascade_engine as E

ROOT = Path(__file__).resolve().parent
CONES_PATH = ROOT / "cascade_cones.json"
JSON_OUT = ROOT / "cascade_cone_certificates.json"
MD_OUT = ROOT / "CASCADE_CONE_LEMMAS.md"

DEPTH = 4
DIM_NAMES = {
    "T1": ("d1", "sigma", "d2", "g7", "g6", "g5", "g4"),
    "T2": ("d1", "sigma", "d2", "g6", "g5", "g4"),
}


def flag_cases(branch: str):
    """Global zero-flag cases, exactly as enumerated by the engine."""

    terminal = 7 if branch == "T1" else 6
    below = list(range(terminal - 1, DEPTH - 1, -1))
    for sigma_zero in (False, True) if branch == "T1" else (False,):
        for d2_zero in (False, True):
            for mask in range(2 ** len(below)):
                g_zero = {
                    level: bool(mask >> i & 1) for i, level in enumerate(below)
                }
                g_zero[terminal] = False
                yield sigma_zero, d2_zero, g_zero


def flag_key(sigma_zero: bool, d2_zero: bool, g_zero: dict[int, bool]) -> str:
    parts = []
    if sigma_zero:
        parts.append("sigma=0")
    if d2_zero:
        parts.append("d2=0")
    zeros = sorted(level for level, flag in g_zero.items() if flag)
    if zeros:
        parts.append("g0:" + ",".join(map(str, zeros)))
    return "|".join(parts) if parts else "generic"


class Tables:
    """Per-place local tables, cached over (branch, a, beta, flags)."""

    def __init__(self, config=None) -> None:
        self.profiles: dict = {}
        self.config = config  # None => audited sub2 behavior

    def get(self, branch: str, a: int, beta: int, sz: bool, dz: bool, gz):
        key = (branch, a, beta, sz, dz, tuple(sorted(gz.items())))
        if key not in self.profiles:
            self.profiles[key] = E.place_profiles(
                beta, branch, 10 + 3 * a, DEPTH, sz, dz, gz,
                self.config, a,
            )
        return self.profiles[key]

    def minima(self, branch: str, a: int, beta: int, sz: bool, dz: bool, gz):
        terminal = 7 if branch == "T1" else 6
        levels = tuple(range(terminal, DEPTH - 1, -1))
        options = self.get(branch, a, beta, sz, dz, gz)
        if not options:
            return None
        vectors = [profile.budget_vector(levels) for profile in options]
        return tuple(
            min(vector[i] for vector in vectors) for i in range(len(vectors[0]))
        )


def certify_branch(tables: Tables, a: int, b, branch: str):
    """Return per-flag-case certificates, or None where a case survives."""

    r_cap = 10 + 3 * a
    terminal = 7 if branch == "T1" else 6
    levels = tuple(range(terminal, DEPTH - 1, -1))
    aux, g_caps, _ = E.resolve_caps(branch, r_cap, tables.config, a)
    caps = (
        aux["d1"],
        aux["sigma"],
        aux["d2"],
        *(g_caps[level] for level in levels),
    )
    names = DIM_NAMES[branch]
    cases = []
    for sz, dz, gz in flag_cases(branch):
        key = flag_key(sz, dz, gz)
        dead = [
            beta
            for beta in b
            if tables.minima(branch, a, beta, sz, dz, gz) is None
        ]
        if dead:
            cases.append(
                {"flags": key, "kind": "L", "dead_place_beta": dead[0]}
            )
            continue
        minima = [tables.minima(branch, a, beta, sz, dz, gz) for beta in b]
        violated = [
            i
            for i in range(len(caps))
            if sum(m[i] for m in minima) > caps[i]
        ]
        if violated:
            i = violated[0]
            cases.append(
                {
                    "flags": key,
                    "kind": "B",
                    "dimension": names[i],
                    "sum_of_minima": int(sum(m[i] for m in minima)),
                    "cap": int(caps[i]),
                    "per_place_minima": [int(m[i]) for m in minima],
                }
            )
        else:
            cases.append(None)
    return cases


def local_kill_table(tables: Tables):
    """(branch, a, beta) -> 'all' | list of flag cases where the place dies."""

    e_cap = 10 if tables.config is None else tables.config.e_cap
    result = {}
    for branch in ("T1", "T2"):
        for a in range(11):
            for beta in range(0, e_cap - a + 1):
                dead_keys = []
                total = 0
                for sz, dz, gz in flag_cases(branch):
                    total += 1
                    if tables.minima(branch, a, beta, sz, dz, gz) is None:
                        dead_keys.append(flag_key(sz, dz, gz))
                if dead_keys:
                    result[(branch, a, beta)] = (
                        "all" if len(dead_keys) == total else dead_keys
                    )
    return result


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--window", choices=("sub2", "sub1"), default="sub2")
    args = parser.parse_args()

    global CONES_PATH, JSON_OUT, MD_OUT
    if args.window == "sub1":
        CONES_PATH = ROOT / "cascade_cones_sub1_depth4.json"
        JSON_OUT = ROOT / "cascade_cone_certificates_sub1.json"
        MD_OUT = ROOT / "CASCADE_CONE_LEMMAS_SUB1.md"

    cones = json.loads(CONES_PATH.read_text(encoding="utf-8"))
    assert cones["depth"] == DEPTH
    tables = Tables(config=None if args.window == "sub2" else E.CONFIGS["sub1"])

    records = []
    mismatches = []
    for row in cones["branches"]:
        a, b, branch = row["a_t"], tuple(row["b"]), row["branch"]
        cases = certify_branch(tables, a, b, branch)
        certified_kill = all(case is not None for case in cases)
        engine_kill = row["status"] == "engine_killed_pending_audit"
        if certified_kill != engine_kill:
            mismatches.append((a, b, branch, engine_kill, certified_kill))
        records.append(
            {
                "a_t": a,
                "b": list(b),
                "branch": branch,
                "certified_kill": certified_kill,
                "cases": [
                    case if case is not None else {"kind": "open"}
                    for case in cases
                ],
            }
        )
    assert not mismatches, f"certificate/engine mismatches: {mismatches[:5]}"

    local_table = local_kill_table(tables)
    unconditional = sorted(
        (branch, a, beta)
        for (branch, a, beta), value in local_table.items()
        if value == "all"
    )

    kills = [r for r in records if r["certified_kill"]]
    kind_counts = {"L": 0, "B": 0}
    dim_counts: dict[str, int] = {}
    for record in kills:
        for case in record["cases"]:
            kind_counts[case["kind"]] += 1
            if case["kind"] == "B":
                key = f"{record['branch']}:{case['dimension']}"
                dim_counts[key] = dim_counts.get(key, 0) + 1

    payload = {
        "schema": 1,
        "description": (
            "Cone certificates covering the cascade-engine depth-4 kills; "
            "generated from the engine's per-place tables (compression, "
            "not independent audit)"
        ),
        "depth": DEPTH,
        "summary": {
            "kills_certified": len(kills),
            "survivors": len(records) - len(kills),
            "case_certificates": kind_counts,
            "budget_dimension_counts": dim_counts,
            "unconditional_local_kills": [
                {"branch": branch, "a_t": a, "beta": beta}
                for branch, a, beta in unconditional
            ],
        },
        "local_kill_table": [
            {
                "branch": branch,
                "a_t": a,
                "beta": beta,
                "flags": value,
            }
            for (branch, a, beta), value in sorted(local_table.items())
        ],
        "branches": records,
    }
    JSON_OUT.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    engine_kills = sum(
        row["status"] == "engine_killed_pending_audit"
        for row in cones["branches"]
    )
    if tables.config is None:
        caps_text = "(`deg d1<=6`, `deg sigma<=8`, `deg d2<=4`, `deg g_l<=10+3a`)."
    else:
        caps_text = (
            "(`deg d1<=9`, `deg sigma<=12`, `deg d2<=6`, per-level g caps"
            " from `sub1_cascade_verify.py`)."
        )
    lines = [
        f"# Cone lemmas for the depth-4 cascade kills ({args.window})",
        "",
        "**Generated by `cone_lemmas.py`; do not hand-edit.**  Certificates",
        "are extracted from the engine's per-place enumeration; independent",
        f"verification is the audit checker's role.  Every one of the {engine_kills}",
        "killed branches carries, for every zero-flag case, either an L or a",
        "B certificate below, and the certificate verdict matches the engine",
        f"verdict on all {len(records)} open branches (asserted at generation time).",
        "",
        "## Lemma family L (single-place kills)",
        "",
        "At one split place with `v_p(e)=beta`, the level identities",
        f"7/6 -> {DEPTH} admit no consistent local valuation chain, killing",
        "every stratum containing `beta` regardless of the other places.",
        "Pairs `(a, beta)` where this holds for EVERY zero-flag case:",
        "",
    ]
    for branch in ("T1", "T2"):
        rows = [(a, beta) for br, a, beta in unconditional if br == branch]
        by_a: dict[int, list[int]] = {}
        for a, beta in rows:
            by_a.setdefault(a, []).append(beta)
        lines.append(f"- **{branch}**: " + "; ".join(
            f"a={a}: beta in {sorted(betas)}" for a, betas in sorted(by_a.items())
        ))
    lines += [
        "",
        "Flag-conditional local kills (dead only under some zero-flag cases)",
        f"are enumerated in `{JSON_OUT.name}`.",
        "",
        "## Lemma family B (single-dimension budget kills)",
        "",
        "All places admit chains, but summing each place's minimum possible",
        "consumption in ONE dimension exceeds its global cap",
        caps_text,
        "Distribution of the certificate dimension over all flag cases:",
        "",
    ]
    for key, count in sorted(dim_counts.items(), key=lambda kv: -kv[1]):
        lines.append(f"- `{key}`: {count} cases")
    lines += [
        "",
        "## Counts",
        "",
        f"- Killed branches certified: **{len(kills)}** / {engine_kills}.",
        f"- Certificate cases: **{kind_counts['L']}** local (L), "
        f"**{kind_counts['B']}** budget (B).",
        f"- Unconditional local-kill pairs: **{len(unconditional)}**.",
        "",
        f"Row-level certificates: `{JSON_OUT.name}`.",
        "",
    ]
    MD_OUT.write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps(payload["summary"], indent=2)[:1200])
    print(f"wrote {JSON_OUT.name} and {MD_OUT.name}")


if __name__ == "__main__":
    main()
