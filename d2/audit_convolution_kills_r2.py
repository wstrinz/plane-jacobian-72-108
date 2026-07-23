#!/usr/bin/env python3
"""Audit the new convolution-kill batches with the round-1 spec-only checker.

Fresh kills run in isolated child interpreters. The 60 sub1 transfers are
only matched, one-for-one, by their window-independent state tuples against
batch_convolution_sub2.json; they are deliberately not symbolically re-proved.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import json
import os
from pathlib import Path
import subprocess
import sys
import time
from typing import Any

import audit_convolution_kills as round1

HERE = Path(__file__).resolve().parent
DEFAULT_TIMEOUT = 120.0
SUB1_FILE = "batch_convolution_sub1.json"
SUB2_ROUND1_FILE = "batch_convolution_sub2.json"
SUB2_ROUND2_FILE = "batch_convolution_sub2_round2.json"
OVERNIGHT_FILE = "batch_convolution_overnight.json"

FRESH_SPECS = {
    "sub1-fresh": (SUB1_FILE, "fresh_kills_pending_audit", 68, 2),
    "sub2-round2": (SUB2_ROUND2_FILE, "kills_pending_audit", 82, 2),
    "overnight": (OVERNIGHT_FILE, "kills_pending_audit", 30, 1),
}
CLASSIFICATIONS = ("CONFIRMED", "UNDECIDED-BY-AUDIT", "DISAGREEMENT")
KILL_VERDICTS = ("CONTRADICTION", "STATE_KILLED_BY_DEGREE_DROP")


class R2AuditError(RuntimeError):
    pass


def need(condition: bool, message: str) -> None:
    if not condition:
        raise R2AuditError(message)


def read_json(filename: str) -> dict[str, Any]:
    data = json.loads((HERE / filename).read_text(encoding="utf-8-sig"))
    need(isinstance(data, dict), f"{filename}: JSON root is not an object")
    return data


def load_fresh_rows(dataset: str) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    need(dataset in FRESH_SPECS, f"unknown fresh dataset {dataset!r}")
    filename, key, expected, schema = FRESH_SPECS[dataset]
    data = read_json(filename)
    rows = data.get(key)
    need(isinstance(rows, list), f"{filename}: {key} is not a list")
    need(len(rows) == expected,
         f"{filename}: expected {expected} {key} rows, found {len(rows)}")
    need(data.get("schema") == schema,
         f"{filename}: expected schema {schema}, found {data.get('schema')!r}")
    need(data.get("c") == "-1/6630", f"{filename}: unexpected c={data.get('c')!r}")
    need(data.get("floor_budget") == 14, f"{filename}: floor_budget is not 14")
    if dataset != "sub1-fresh":
        need(data.get("kill_count") == expected,
             f"{filename}: kill_count is not {expected}")
    return data, rows


def claimed_verdict(row: dict[str, Any]) -> str:
    """Normalize the schema-1/schema-2 fresh claim field conservatively."""
    has_verdict = "verdict" in row
    has_final = "final_verdict" in row
    need(has_verdict != has_final, "record must have exactly one kill-claim field")
    verdict = row["verdict"] if has_verdict else row["final_verdict"]
    need(verdict in KILL_VERDICTS, f"unsupported kill verdict {verdict!r}")
    return verdict


def normalized_degree(row: dict[str, Any], degree_key: str,
                      zero_key: str) -> int | None:
    need(zero_key in row, f"missing tuple flag {zero_key}")
    zero = row[zero_key]
    need(type(zero) is bool, f"{zero_key} is not Boolean")
    need(degree_key in row, f"missing tuple degree {degree_key}")
    value = row[degree_key]
    if zero:
        need(value in ("-inf", None),
             f"{degree_key} is finite while {zero_key}=true")
        return None
    need(value not in ("-inf", None),
         f"{degree_key} is absent while {zero_key}=false")
    try:
        degree = int(value)
    except (TypeError, ValueError) as exc:
        raise R2AuditError(f"{degree_key} is not an integer: {value!r}") from exc
    need(degree >= 0, f"{degree_key} is negative")
    return degree


def state_tuple(row: dict[str, Any]) -> tuple[Any, ...]:
    """Return (a_t, flags, degrees), the window-independent match key."""
    need("a_t" in row and "deg_e" in row, "tuple lacks a_t or deg_e")
    a_t, deg_e = int(row["a_t"]), int(row["deg_e"])
    flags = (row.get("d1_zero"), row.get("sigma_zero"), row.get("d2_zero"))
    need(all(type(flag) is bool for flag in flags),
         "tuple zero flags must be Booleans")
    deg_d2 = normalized_degree(row, "deg_d2", "d2_zero")
    deg_d1 = normalized_degree(row, "deg_d1", "d1_zero")
    deg_sigma = normalized_degree(row, "deg_sigma", "sigma_zero")
    return (a_t, *flags, deg_d2, deg_d1, deg_sigma, deg_e)


def format_tuple(value: tuple[Any, ...]) -> str:
    names = ("a_t", "d1_zero", "sigma_zero", "d2_zero",
             "deg_d2", "deg_d1", "deg_sigma", "deg_e")
    return "(" + ", ".join(
        f"{key}={item!r}" for key, item in zip(names, value)) + ")"


def audit_transfers() -> list[dict[str, Any]]:
    """Match transfers one-for-one; never invoke the symbolic prover here."""
    sub1 = read_json(SUB1_FILE)
    source = read_json(SUB2_ROUND1_FILE)
    transfers = sub1.get("transferred_kills")
    source_rows = source.get("kills_pending_audit")
    need(isinstance(transfers, list),
         f"{SUB1_FILE}: transferred_kills is not a list")
    need(isinstance(source_rows, list),
         f"{SUB2_ROUND1_FILE}: kills_pending_audit is not a list")
    need(len(transfers) == 60,
         f"{SUB1_FILE}: expected 60 transfers, found {len(transfers)}")
    need(len(source_rows) == 60,
         f"{SUB2_ROUND1_FILE}: expected 60 source kills, found {len(source_rows)}")
    need(sub1.get("sub2_source") == SUB2_ROUND1_FILE,
         f"{SUB1_FILE}: bad sub2_source")

    by_tuple: dict[tuple[Any, ...], list[tuple[int, dict[str, Any]]]] = defaultdict(list)
    source_errors: list[str] = []
    for source_index, row in enumerate(source_rows):
        try:
            key = state_tuple(row)
            verdict = row.get("final_verdict")
            need(verdict in KILL_VERDICTS,
                 f"unsupported source verdict {verdict!r}")
            by_tuple[key].append((source_index, row))
        except Exception as exc:
            source_errors.append(
                f"source K{source_index + 1:03d}: {type(exc).__name__}: {exc}")

    results: list[dict[str, Any]] = []
    used_source_indices: set[int] = set()
    for index, row in enumerate(transfers):
        identifier = f"SUB1-TRANSFER-{index + 1:03d}"
        started = time.perf_counter()
        try:
            need(not source_errors, "; ".join(source_errors))
            key = state_tuple(row)
            candidates = by_tuple.get(key, [])
            need(candidates,
                 f"no round-1 source tuple equals {format_tuple(key)}")
            need(len(candidates) == 1,
                 f"ambiguous source tuple has {len(candidates)} matches: "
                 f"{format_tuple(key)}")
            source_index, source_row = candidates[0]
            need(source_index not in used_source_indices,
                 f"round-1 source K{source_index + 1:03d} matched twice")
            transfer_verdict = row.get("transferred_verdict")
            source_verdict = source_row.get("final_verdict")
            need(transfer_verdict == source_verdict,
                 f"tuple matches but verdict differs: transfer={transfer_verdict!r}, "
                 f"source={source_verdict!r}")
            used_source_indices.add(source_index)
            result = {
                "identifier": identifier,
                "classification": "CONFIRMED",
                "mechanism": "TUPLE_IDENTITY",
                "detail": f"one-to-one match to round-1 K{source_index + 1:03d}",
                "tuple": format_tuple(key),
            }
        except Exception as exc:
            result = {
                "identifier": identifier,
                "classification": "DISAGREEMENT",
                "mechanism": "TUPLE_IDENTITY_MISMATCH",
                "detail": f"{type(exc).__name__}: {exc}",
            }
        result["runtime_seconds"] = time.perf_counter() - started
        results.append(result)

    unmatched = sorted(set(range(len(source_rows))) - used_source_indices)
    if unmatched and all(r["classification"] == "CONFIRMED" for r in results):
        results[-1]["classification"] = "DISAGREEMENT"
        results[-1]["mechanism"] = "TUPLE_IDENTITY_MISMATCH"
        results[-1]["detail"] = "unmatched round-1 sources: " + ", ".join(
            f"K{index + 1:03d}" for index in unmatched)
    return results


def fresh_identifier(dataset: str, index: int, row: dict[str, Any]) -> str:
    base = round1.identifier(index, row)
    suffix = base.split(":", 1)[1] if ":" in base else base
    return f"{dataset.upper()}-{index + 1:03d}:{suffix}"


def undecided_result(dataset: str, index: int, detail: str) -> dict[str, Any]:
    return {
        "identifier": f"{dataset.upper()}-{index + 1:03d}",
        "classification": "UNDECIDED-BY-AUDIT",
        "mechanism": "WORKER_ERROR",
        "detail": detail,
        "trace": [],
        "runtime_seconds": 0.0,
    }


def audit_fresh_worker(dataset: str, index: int) -> int:
    try:
        _, rows = load_fresh_rows(dataset)
        need(0 <= index < len(rows), f"worker index {index} is out of range")
        row = rows[index]
        need(isinstance(row, dict), "fresh record is not an object")
        expected = claimed_verdict(row)
        normalized = dict(row)
        normalized["final_verdict"] = expected
        graded = round1.parse_graded()
        result = round1.audit_record(index, normalized, graded)
        result["identifier"] = fresh_identifier(dataset, index, row)
    except Exception as exc:
        result = undecided_result(dataset, index,
                                  f"{type(exc).__name__}: {exc}")
    print(json.dumps(result, sort_keys=True))
    return 0


def run_fresh_child(dataset: str, index: int, row: dict[str, Any],
                    timeout: float) -> dict[str, Any]:
    command = [sys.executable, str(Path(__file__).resolve()),
               "--worker-dataset", dataset, "--worker-index", str(index)]
    environment = dict(os.environ)
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    started = time.perf_counter()
    try:
        done = subprocess.run(command, cwd=HERE, env=environment,
                              capture_output=True, text=True,
                              timeout=timeout, check=False)
    except subprocess.TimeoutExpired:
        result = undecided_result(
            dataset, index, f"child exceeded {timeout:g}s and was terminated")
        result["mechanism"] = "TIMEOUT"
        result["runtime_seconds"] = time.perf_counter() - started
        try:
            result["identifier"] = fresh_identifier(dataset, index, row)
        except Exception:
            pass
        return result
    if done.returncode:
        result = undecided_result(
            dataset, index,
            f"child exit {done.returncode}; stderr={done.stderr.strip()}")
        result["runtime_seconds"] = time.perf_counter() - started
        return result
    try:
        return json.loads(done.stdout)
    except json.JSONDecodeError as exc:
        result = undecided_result(
            dataset, index,
            f"invalid child JSON: {exc}; stdout={done.stdout[-1000:]!r}")
        result["runtime_seconds"] = time.perf_counter() - started
        return result


def print_fresh_record(result: dict[str, Any]) -> None:
    print(f"{result['identifier']} {result['classification']} "
          f"claim={result.get('expected')} mechanism={result.get('mechanism')} "
          f"degree={result.get('degree')} "
          f"seconds={result.get('runtime_seconds', 0.0):.3f}")
    print(f"  ansatz: {result.get('ansatz', '(worker did not reconstruct ansatz)')}")
    for step in result.get("trace", []):
        move = step.get("move")
        if move == "IDENTITY":
            continue
        fields = ", ".join(
            f"{key}={value}" for key, value in step.items() if key != "move")
        print(f"  {move}: {fields}")
    if result.get("classification") != "CONFIRMED":
        print(f"  DETAIL: {result.get('detail')}")
        if result.get("residual"):
            print(f"  RESIDUAL: {result['residual']}")
        if result.get("substitutions"):
            print(f"  SUBSTITUTIONS: {result['substitutions']}")
        if result.get("parameter_constraints"):
            print(f"  PARAMETER_CONSTRAINTS: {result['parameter_constraints']}")


def print_transfer_record(result: dict[str, Any]) -> None:
    print(f"{result['identifier']} {result['classification']} "
          f"mechanism={result['mechanism']} "
          f"seconds={result['runtime_seconds']:.6f}")
    print(f"  {result['detail']}")
    if result.get("tuple"):
        print(f"  tuple: {result['tuple']}")


def census(results: list[dict[str, Any]]) -> dict[str, int]:
    return {label: sum(row.get("classification") == label for row in results)
            for label in CLASSIFICATIONS}


def print_census(label: str, results: list[dict[str, Any]]) -> None:
    counts = census(results)
    print(f"{label}: TOTAL={len(results)} CONFIRMED={counts['CONFIRMED']} "
          f"UNDECIDED-BY-AUDIT={counts['UNDECIDED-BY-AUDIT']} "
          f"DISAGREEMENT={counts['DISAGREEMENT']}")
    undecided = [row["identifier"] for row in results
                 if row.get("classification") == "UNDECIDED-BY-AUDIT"]
    disagreements = [row["identifier"] for row in results
                     if row.get("classification") == "DISAGREEMENT"]
    print(f"{label} UNDECIDED: {', '.join(undecided) if undecided else 'none'}")
    print(f"{label} DISAGREEMENTS: "
          f"{', '.join(disagreements) if disagreements else 'none'}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quiet", action="store_true",
                        help="suppress per-kill output; print setup and censuses")
    parser.add_argument("--timeout", type=float, default=DEFAULT_TIMEOUT,
                        help=f"hard timeout per fresh kill (default {DEFAULT_TIMEOUT:g}s)")
    parser.add_argument("--worker-dataset", choices=tuple(FRESH_SPECS),
                        help=argparse.SUPPRESS)
    parser.add_argument("--worker-index", type=int, help=argparse.SUPPRESS)
    args = parser.parse_args(argv)

    if args.worker_dataset is not None or args.worker_index is not None:
        if args.worker_dataset is None or args.worker_index is None:
            print("worker dataset and index must be supplied together",
                  file=sys.stderr)
            return 2
        return audit_fresh_worker(args.worker_dataset, args.worker_index)
    if args.timeout <= 0:
        print("AUDIT SETUP FAILED: --timeout must be positive", file=sys.stderr)
        return 2

    started = time.perf_counter()
    try:
        graded = round1.parse_graded()
        loaded = {dataset: load_fresh_rows(dataset) for dataset in FRESH_SPECS}
        transfers = audit_transfers()
    except Exception as exc:
        print(f"AUDIT SETUP FAILED: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 2

    print("parser_and_homogeneity: PASS "
          f"(h_0..h_7; sigma term counts={graded.counts}; "
          "h_5/h_6/h_7 spot checks)")
    print("ansatz_semantics: round-1 degree-exact leading coefficients are NONZERO")
    print(f"process_isolation: PASS (hard timeout {args.timeout:g}s per fresh kill)")
    print("schema_notes: sub1 schema=2 uses fresh_kills_pending_audit + verdict; "
          "sub2 round2 schema=2 uses kills_pending_audit + verdict; overnight "
          "schema=1 retains the common state tuple/ansatz fields + verdict")
    print("transfer_semantics: tuple identity only; no transferred kill was "
          "symbolically reconstructed or re-proved")

    if not args.quiet:
        for result in transfers:
            print_transfer_record(result)

    fresh_results: dict[str, list[dict[str, Any]]] = {}
    for dataset in FRESH_SPECS:
        _, rows = loaded[dataset]
        results: list[dict[str, Any]] = []
        for index, row in enumerate(rows):
            if not isinstance(row, dict):
                result = undecided_result(
                    dataset, index, "fresh record is not an object")
            else:
                result = run_fresh_child(dataset, index, row, args.timeout)
            results.append(result)
            if not args.quiet:
                print_fresh_record(result)
        fresh_results[dataset] = results

    print_census("SUB1 FRESH CENSUS", fresh_results["sub1-fresh"])
    transfer_counts = census(transfers)
    print("SUB1 TRANSFER IDENTITY: "
          f"TOTAL={len(transfers)} CONFIRMED={transfer_counts['CONFIRMED']} "
          f"NOT_CONFIRMED={len(transfers) - transfer_counts['CONFIRMED']}")
    print_census("SUB1 FILE CENSUS (fresh + transferred identity)",
                 fresh_results["sub1-fresh"] + transfers)
    print_census("SUB2 ROUND2 FILE CENSUS", fresh_results["sub2-round2"])
    print_census("OVERNIGHT FILE CENSUS", fresh_results["overnight"])

    all_results = (transfers + fresh_results["sub1-fresh"]
                   + fresh_results["sub2-round2"]
                   + fresh_results["overnight"])
    print_census("GRAND CENSUS", all_results)
    claim_counts = Counter(
        result.get("expected")
        for dataset_results in fresh_results.values()
        for result in dataset_results)
    print("FRESH CLAIM CENSUS: "
          f"CONTRADICTION={claim_counts['CONTRADICTION']} "
          "STATE_KILLED_BY_DEGREE_DROP="
          f"{claim_counts['STATE_KILLED_BY_DEGREE_DROP']} "
          f"TOTAL={sum(claim_counts.values())}")
    elapsed = time.perf_counter() - started
    print(f"TOTAL RUNTIME SECONDS: {elapsed:.3f}")

    fully_confirmed = (
        len(all_results) == 240
        and all(row.get("classification") == "CONFIRMED"
                for row in all_results)
    )
    print(f"FINAL EXIT CONDITION: {'PASS' if fully_confirmed else 'FAIL'}")
    return 0 if fully_confirmed else 1


if __name__ == "__main__":
    raise SystemExit(main())
