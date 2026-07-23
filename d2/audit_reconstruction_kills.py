#!/usr/bin/env python3
"""Independent spec-only audit."""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import time
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import sympy as sp

BASE = Path(__file__).resolve().parent
Q_FILE = BASE / "f31_graded.txt"
ALT_COMBINED = BASE / "alt_combined.json"
ALT_RESIDUE = BASE / "alt_residue_congruences.json"
SCALE = BASE / "phase_f2_scale.json"
SOURCES = ("PHASE_F_PLAN.md", "PHASE_F2_PILOT.md", "PHASE_F2_SCALE.md",
           "ALT_REGIME_INF.md", "ALT_RESIDUE_CONGRUENCES.md",
           "alt_combined.json", "alt_residue_congruences.json",
           "phase_f2_scale.json", "f31_graded.txt")
FORBIDDEN = ("phase_f2_pilot.py", "phase_f2_scale.py",
             "alt_residue_congruences.py", "alt_inf_sweep.py",
             "cascade_engine.py")
Q_COEFFS = (2048, -512, 320, -240, 195)
C_U = sp.Rational(-1, 6630)
EXPECTED_LC_U = sp.Rational(-1024, 3315)
TOTAL_CLAIMS = 22


class AuditDisagreement(RuntimeError):
    pass


@dataclass(frozen=True)
class Case:
    label: str
    family: str
    deg_d2: int | None
    deg_d1: int | None
    deg_sigma: int
    deg_e: int
    tie_depth: int
    check_depth: int
    expected_j0: str


@dataclass
class Reconstruction:
    d2: sp.Expr
    d1: sp.Expr
    sigma: sp.Expr
    e: sp.Expr
    variables: list[sp.Symbol]
    nonzero: list[sp.Expr]
    field_equations: list[sp.Expr]
    valuations: dict[str, dict[str, int]]
    field: str


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def same_state(a: dict[str, Any], b: dict[str, Any]) -> bool:
    keys = ("deg_d2", "deg_d1", "deg_sigma", "deg_e", "deg_E")
    return all(a.get(k) == b.get(k) for k in keys)


def parse_f31_graded() -> tuple[dict[int, sp.Expr], sp.Expr]:
    """Own parser plus original/re-written homogeneity self-checks."""
    d0, d2, d1, dm1 = sp.symbols("d0 d2 d1 dm1")
    local = {"d0": d0, "d2": d2, "d1": d1, "dm1": dm1}
    parsed: dict[int, sp.Expr] = {}
    declared: dict[int, tuple[int, int]] = {}
    pat = re.compile(
        r"^h_(\d+)\s+\(weight\s+(\d+),\s+dm1-power\s+(\d+)\)\s*=\s*(.+)$")
    for raw in Q_FILE.read_text(encoding="utf-8").splitlines():
        m = pat.match(raw.strip())
        if m:
            f, weight, power = map(int, m.group(1, 2, 3))
            if f in parsed:
                raise AuditDisagreement(f"duplicate h_{f}")
            parsed[f] = sp.sympify(m.group(4), locals=local)
            declared[f] = weight, power
    if set(parsed) != set(range(8)):
        raise AuditDisagreement(f"own parser found levels {sorted(parsed)}")
    weights = (4, 2, 3, 5)
    for f in range(8):
        target = 20 - 2 * f
        if declared[f] != (target, 21 - 3 * f):
            raise AuditDisagreement(f"h_{f} header mismatch")
        for monom, _ in sp.Poly(parsed[f], d0, d2, d1, dm1).terms():
            if sum(w * e for w, e in zip(weights, monom)) != target:
                raise AuditDisagreement(f"h_{f} is not homogeneous")
    sigma = sp.Symbol("sigma")
    h0 = sp.Poly(sp.cancel(parsed[0].subs(d0, (sigma + d2**2) / 4)),
                 d2, d1, sigma, dm1).as_expr()
    if sp.expand(h0.subs(sigma, 4 * d0 - d2**2) - parsed[0]) != 0:
        raise AuditDisagreement("sigma rewrite does not round-trip")
    terms = sp.Poly(h0, d2, d1, sigma, dm1).terms()
    if len(terms) != 26:
        raise AuditDisagreement("rewritten h_0 is not the 26-term table")
    for monom, _ in terms:
        if sum(w * e for w, e in zip((2, 3, 4, 5), monom)) != 20:
            raise AuditDisagreement("rewritten h_0 is not homogeneous")
    return parsed, h0


def source_cases() -> tuple[list[Case], dict[str, Any]]:
    """Cross-link the 20 scale keys to both source data sets."""
    scale, residue, combined = map(read_json, (SCALE, ALT_RESIDUE, ALT_COMBINED))
    killed = [row for row in scale["alt_states"] if row.get("verdict") == "KILLED"]
    if len(killed) != 20:
        raise AuditDisagreement(f"scale contains {len(killed)} kills, not 20")
    if scale["census_summary"]["verdict_census"].get("KILLED") != 20:
        raise AuditDisagreement("scale census does not say 20 kills")
    branches = {b["id"]: b for b in combined["branches"]}
    cases = [
        Case("PILOT-A:a11_b3000_T1:(0,9,12,14)", "pilot_a",
             0, 9, 12, 14, 14, 3, "2187*S**2*(4*S**3 + X**4)"),
        Case("PILOT-B:a14_b0000_T2:(6,zero,12,14)", "pilot_b",
             6, None, 12, 14, 14, 9,
             "12*S**2*(4*D6**2 + 9*S)**2*(5*D6**2 + 9*S)"),
    ]
    family_map = {"a11_b1111_T1": "flagship",
                  "a11_b3100_T2": "two_root",
                  "a12_b1110_T2": "complement",
                  "a14_b0000_T2": "pilot_b"}
    counts: Counter[str] = Counter()
    for row in killed:
        m = re.fullmatch(r"(.+)#sup(\d+)#idx(\d+)", row["key"])
        if not m:
            raise AuditDisagreement(f"malformed scale key {row['key']}")
        bid, support, ri = m.group(1), int(m.group(2)), int(m.group(3))
        if bid != row["bid"] or support != row["support"] or bid not in family_map:
            raise AuditDisagreement(f"scale metadata mismatch for {row['key']}")
        state = dict(zip(("deg_d2", "deg_d1", "deg_sigma", "deg_e"),
                         row["degs"]))
        state["deg_E"] = state["deg_e"] - branches[bid]["a"]
        rr = residue["states"][ri]
        if rr["branch"] != bid or not same_state(rr["state"], state):
            raise AuditDisagreement(f"residue state mismatch for {row['key']}")
        if (rr["L0_tie_support_id"], rr["L0_tie_depth"]) != (
                support, row["tie_depth"]):
            raise AuditDisagreement(f"tie metadata mismatch for {row['key']}")
        matches = [s for s in branches[bid]["remaining_states"]
                   if same_state(s["state"], state)]
        if len(matches) != 1:
            raise AuditDisagreement(f"source state count mismatch for {row['key']}")
        family = family_map[bid]
        counts[family] += 1
        cases.append(Case("SCALE:" + row["key"], family, *row["degs"],
                          row["tie_depth"], row["kill_depth"], row["j0"]))
    if counts != Counter(flagship=7, two_root=6, complement=6, pilot_b=1):
        raise AuditDisagreement(f"unexpected family census {counts}")
    if len(cases) != TOTAL_CLAIMS:
        raise AuditDisagreement(f"built {len(cases)} claims")
    return cases, {"scale": scale, "residue": residue, "combined": combined,
                   "branches": branches}


def q_polynomial(y: sp.Symbol) -> sp.Expr:
    return sum(c * y ** (4 - i) for i, c in enumerate(Q_COEFFS))


def monic_q(y: sp.Symbol) -> sp.Expr:
    return sp.expand(q_polynomial(y) / 2048)


def complement_factor(y: sp.Symbol, r: sp.Symbol) -> sp.Expr:
    quotient, remainder = sp.div(monic_q(y), y - r, y)
    if sp.expand(remainder - monic_q(y).subs(y, r)) != 0:
        raise AuditDisagreement("synthetic division of q failed")
    if sp.Poly(quotient, y).degree() != 3 or sp.Poly(quotient, y).LC() != 1:
        raise AuditDisagreement("complement is not monic cubic")
    return sp.expand(quotient)


def reconstruction(case: Case, meta: dict[str, Any]) -> Reconstruction:
    """Rebuild every polynomial from documented per-place exponents."""
    y = sp.Symbol("y")
    X, S, E = sp.symbols("X S E")
    variables: list[sp.Symbol] = []
    nonzero: list[sp.Expr] = []
    field_eqs: list[sp.Expr] = []
    field = "Q"
    if case.deg_d2 is None:
        d2p = sp.Integer(0)
    else:
        ds = list(sp.symbols(f"D0:{case.deg_d2 + 1}"))
        d2p = sum(ds[i] * y**i for i in range(len(ds)))
        variables += ds
        nonzero.append(ds[-1])

    if case.family == "pilot_a":
        r = sp.Symbol("r")
        d1p = X * (y + 1)**5 * (y - r)**4
        sigp = S * (y + 1)**12
        ep = E * (y + 1)**11 * (y - r)**3
        variables += [X, S, E, r]
        nonzero += [X, S, E]
        field_eqs = [q_polynomial(r)]
        vals = {"d1": {"t": 5, "r": 4}, "sigma": {"t": 12},
                "e": {"t": 11, "r": 3}}
        field = "Q[r]/(q)"
    elif case.family == "flagship":
        d1p = X * (y + 1)**5 * monic_q(y)
        sigp = S * (y + 1)**3
        ep = E * (y + 1)**11 * monic_q(y)
        variables += [X, S, E]
        nonzero += [X, S, E]
        vals = {"d1": {"t": 5, "q1": 1, "q2": 1, "q3": 1, "q4": 1},
                "sigma": {"t": 3},
                "e": {"t": 11, "q1": 1, "q2": 1, "q3": 1, "q4": 1}}
    elif case.family == "two_root":
        r1, r2 = sp.symbols("r1 r2")
        d1p = sp.Integer(0)
        sigp = S * (y + 1)**3 * (y - r1)**7 * (y - r2)**2
        ep = E * (y + 1)**11 * (y - r1)**3 * (y - r2)
        variables += [S, E, r1, r2]
        nonzero += [S, E, r1 - r2]
        field_eqs = [q_polynomial(r1), q_polynomial(r2)]
        vals = {"sigma": {"t": 3, "r1": 7, "r2": 2},
                "e": {"t": 11, "r1": 3, "r2": 1}}
        field = "Q[r1,r2]/(q,q), r1!=r2"
    elif case.family == "complement":
        r = sp.Symbol("r")
        comp = complement_factor(y, r)
        d1p = sp.Integer(0)
        sigp = S * (y + 1)**6 * comp**2
        ep = E * (y + 1)**12 * comp
        variables += [S, E, r]
        nonzero += [S, E]
        field_eqs = [q_polynomial(r)]
        vals = {"sigma": {"t": 6, "q1": 2, "q2": 2, "q3": 2},
                "e": {"t": 12, "q1": 1, "q2": 1, "q3": 1}}
        field = "Q[r]/(q), r marks inactive root"
    elif case.family == "pilot_b":
        d1p = sp.Integer(0)
        sigp = S * (y + 1)**12
        ep = E * (y + 1)**14
        variables += [S, E]
        nonzero += [S, E]
        vals = {"sigma": {"t": 12}, "e": {"t": 14}}
    else:
        raise AuditDisagreement(f"unknown family {case.family}")

    expected = {"d2": case.deg_d2, "d1": case.deg_d1,
                "sigma": case.deg_sigma, "e": case.deg_e}
    actual = {"d2": None if d2p == 0 else sp.Poly(d2p, y).degree(),
              "d1": None if d1p == 0 else sp.Poly(d1p, y).degree(),
              "sigma": sp.Poly(sigp, y).degree(), "e": sp.Poly(ep, y).degree()}
    if actual != expected:
        raise AuditDisagreement(f"{case.label}: degrees {actual} != {expected}")
    for name, place_vals in vals.items():
        if sum(place_vals.values()) != expected[name]:
            raise AuditDisagreement(f"{case.label}: {name} is not defect zero")

    # The compact alt_combined witness records the placewise tight lower
    # valuations. Since their sum equals the degree, each coordinate is forced.
    if case.family == "pilot_a":
        bid, deg_e = "a11_b3000_T1", 3
    elif case.family == "pilot_b":
        bid, deg_e = "a14_b0000_T2", 0
    else:
        bid = {"flagship": "a11_b1111_T1", "two_root": "a11_b3100_T2",
               "complement": "a12_b1110_T2"}[case.family]
        deg_e = case.deg_e - meta["branches"][bid]["a"]
    state = {"deg_d2": case.deg_d2, "deg_d1": case.deg_d1,
             "deg_sigma": case.deg_sigma, "deg_e": case.deg_e, "deg_E": deg_e}
    records = [s for s in meta["branches"][bid]["remaining_states"]
               if same_state(s["state"], state)]
    if len(records) != 1:
        raise AuditDisagreement(f"{case.label}: source-state match is not unique")
    witness = records[0]["finite_place_witness"]
    if case.family in ("pilot_a", "flagship"):
        if witness.get("X") != case.deg_d1:
            raise AuditDisagreement(f"{case.label}: X witness mismatch")
        if case.family == "flagship" and witness.get("Z") != case.deg_sigma:
            raise AuditDisagreement(f"{case.label}: Z witness mismatch")
    else:
        details = {
            "two_root": [["t", 3], ["q(b=3)", 7], ["q(b=1)", 2]],
            "complement": [["t", 6], ["q(b=1)", 2], ["q(b=1)", 2],
                           ["q(b=1)", 2]],
            "pilot_b": [["t", 12]],
        }
        if (witness.get("detail") != details[case.family]
                or witness.get("Zmin") != case.deg_sigma):
            raise AuditDisagreement(f"{case.label}: tight split mismatch")
    return Reconstruction(*(sp.expand(p) for p in (d2p, d1p, sigp, ep)),
                          variables, nonzero, field_eqs, vals, field)


def reversed_coefficients(expr: sp.Expr, y: sp.Symbol, depth: int) -> list[sp.Expr]:
    if expr == 0:
        return [sp.Integer(0)] * depth
    p = sp.Poly(expr, y)
    return [p.nth(p.degree() - k) if k <= p.degree() else sp.Integer(0)
            for k in range(depth)]


def trunc_mul(a: list[sp.Expr], b: list[sp.Expr], depth: int) -> list[sp.Expr]:
    out = [sp.Integer(0)] * depth
    for i, ai in enumerate(a[:depth]):
        if ai:
            for j, bj in enumerate(b[:depth - i]):
                if bj:
                    out[i + j] += ai * bj
    return out


def trunc_pow(a: list[sp.Expr], power: int, depth: int) -> list[sp.Expr]:
    out = [sp.Integer(1)] + [sp.Integer(0)] * (depth - 1)
    for _ in range(power):
        out = trunc_mul(out, a, depth)
    return out


def tie_equations(h0: sp.Expr, recon: Reconstruction, depth: int) -> list[sp.Expr]:
    """Top coefficients of h_0, rebuilt by independent truncated convolution."""
    y = sp.Symbol("y")
    source_vars = sp.symbols("d2 d1 sigma dm1")
    polys = (recon.d2, recon.d1, recon.sigma, recon.e)
    degrees = [None if p == 0 else sp.Poly(p, y).degree() for p in polys]
    rev = [reversed_coefficients(p, y, depth) for p in polys]
    result = [sp.Integer(0)] * depth
    for exponents, coefficient in sp.Poly(h0, *source_vars).terms():
        if any(exponents[i] and polys[i] == 0 for i in range(4)):
            continue
        term_degree = sum(exponents[i] * degrees[i] for i in range(4)
                          if exponents[i])
        offset = 60 - term_degree
        if offset >= depth:
            continue
        series = [sp.Integer(1)] + [sp.Integer(0)] * (depth - 1)
        for i, exponent in enumerate(exponents):
            if exponent:
                series = trunc_mul(series, trunc_pow(rev[i], exponent, depth), depth)
        for k in range(offset, depth):
            result[k] += coefficient * series[k - offset]
    return [sp.expand(eq) for eq in result]


def groebner_unit(case_index: int) -> dict[str, Any]:
    started = time.perf_counter()
    cases, meta = source_cases()
    case = cases[case_index]
    _, h0 = parse_f31_graded()
    y = sp.Symbol("y")
    u = sp.expand(C_U * q_polynomial(y))
    if sp.Poly(u, y).degree() != 4 or sp.Poly(u, y).LC() != EXPECTED_LC_U:
        raise AuditDisagreement("u=c*q check failed")
    constants = meta["residue"]["schema"]["constants"]
    if (sp.Rational(constants["lc_u"]) != EXPECTED_LC_U
            or constants["bottom_close"] != "E^21 h_0 + u r_0 = 0"):
        raise AuditDisagreement("bottom-close tie semantics mismatch")
    recon = reconstruction(case, meta)
    equations = tie_equations(h0, recon, case.check_depth)
    names = {str(s): s for s in set().union(*(e.free_symbols for e in equations))}
    for name in ("D6", "X", "S", "E"):
        names.setdefault(name, sp.Symbol(name))
    expected = sp.sympify(case.expected_j0, locals=names)
    if not equations or sp.expand(equations[0] - expected) != 0:
        raise AuditDisagreement(f"{case.label}: rebuilt j0 mismatch")

    z = sp.Symbol("rabinowitsch")
    generators = equations + recon.field_equations + [z * sp.prod(recon.nonzero) - 1]
    variables: list[sp.Symbol] = []
    for var in recon.variables + [z]:
        if var not in variables:
            variables.append(var)
    missing = set().union(*(g.free_symbols for g in generators)) - set(variables)
    if missing:
        raise AuditDisagreement(f"unlisted Groebner variables {missing}")
    # Buchberger is much faster than SymPy F5B on these sparse saturations.
    basis = sp.groebner(generators, *variables, order="grevlex",
                        method="buchberger")
    unit = len(basis.polys) == 1 and basis.polys[0].as_expr() == 1
    return {"label": case.label, "unit": unit, "field": recon.field,
            "check_depth": case.check_depth, "tie_depth": case.tie_depth,
            "basis_size": len(basis.polys),
            "elapsed": time.perf_counter() - started}


def worker_main(index: int) -> int:
    try:
        print(json.dumps({"status": "ok", "result": groebner_unit(index)},
                         sort_keys=True))
        return 0
    except MemoryError as exc:
        print(json.dumps({"status": "budget", "error": str(exc)}))
        return 2
    except Exception as exc:
        print(json.dumps({"status": "error", "error_type": type(exc).__name__,
                          "error": str(exc)}, sort_keys=True))
        return 3


def run_audit(timeout: float, quiet: bool) -> int:
    wall_started = time.perf_counter()
    cases, _ = source_cases()
    parse_f31_graded()
    counts: Counter[str] = Counter()
    for index, case in enumerate(cases):
        started = time.perf_counter()
        try:
            done = subprocess.run(
                [sys.executable, "-I", str(Path(__file__).resolve()),
                 "--worker", str(index)], cwd=str(BASE), text=True,
                capture_output=True, timeout=timeout, check=False)
            lines = [line for line in done.stdout.splitlines() if line.strip()]
            payload = json.loads(lines[-1]) if lines else {
                "status": "error", "error": "worker produced no output"}
            if payload["status"] == "budget":
                verdict, detail = "UNDECIDED-BY-AUDIT(budget)", "memory budget"
            elif payload["status"] != "ok":
                verdict = "DISAGREEMENT"
                detail = f"{payload.get('error_type', 'worker')}: {payload['error']}"
            elif payload["result"]["unit"]:
                verdict = "CONFIRMED"
                r = payload["result"]
                detail = (f"unit ideal; depth={r['check_depth']}; "
                          f"field={r['field']}; GB={r['basis_size']}")
            else:
                verdict = "DISAGREEMENT"
                detail = f"non-unit basis size {payload['result']['basis_size']}"
        except subprocess.TimeoutExpired:
            verdict, detail = ("UNDECIDED-BY-AUDIT(budget)",
                               f"timeout after {timeout:g}s")
        except Exception as exc:
            verdict, detail = "DISAGREEMENT", f"{type(exc).__name__}: {exc}"
        elapsed = time.perf_counter() - started
        counts[verdict] += 1
        if not quiet:
            print(f"[{index + 1:02d}/{len(cases)}] {case.label}")
            print(f"  {verdict} ({elapsed:.3f}s) -- {detail}", flush=True)
    wall = time.perf_counter() - wall_started
    c, u, d = (counts["CONFIRMED"], counts["UNDECIDED-BY-AUDIT(budget)"],
               counts["DISAGREEMENT"])
    status = 0 if c == TOTAL_CLAIMS else 1
    print(f"FINAL CENSUS: CONFIRMED={c} / UNDECIDED-BY-AUDIT={u} / "
          f"DISAGREEMENT={d} / TOTAL={len(cases)}")
    if not quiet:
        print(f"WALL-CLOCK: {wall:.3f}s")
    print(f"EXIT STATUS: {status}")
    return status


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Independent spec-only audit of 22 reconstruction kills")
    parser.add_argument("--quiet", action="store_true",
                        help="print only final census and exit status")
    parser.add_argument("--timeout", type=float, default=120.0,
                        help="per-kill subprocess timeout in seconds")
    parser.add_argument("--worker", type=int, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.timeout <= 0:
        parser.error("--timeout must be positive")
    if args.worker is not None:
        return worker_main(args.worker)
    return run_audit(args.timeout, args.quiet)


if __name__ == "__main__":
    raise SystemExit(main())
