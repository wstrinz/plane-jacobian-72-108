#!/usr/bin/env python3
"""Produce cofactor certificates for 15 timed-out lifts and four J6 kills.

msolve has already decided these saturated systems; this script turns those
decisions into portable ``d2-kill-certificate-v1`` identities.  A verdict is
never treated as a certificate: success is emitted only after exact SymPy
expansion verifies ``sum(c_i*f_i) == 1`` over QQ.

Route A runs Singular's standard-basis transformation lift over prime fields,
groups answers by complete monomial-support signature (isolating unlucky
primes), CRT-combines compatible coefficients, and rationally reconstructs
them.  Reconstruction is tried after every prime and accepted only at the exact
expansion gate.  Generators are primitive integer polynomials before printing;
the emitted Singular source therefore contains no divisions (avoiding the
``gm^8/N`` parser trap).  Ring variables are renamed V0,V1,... inside Singular
to avoid SHORT output and generator-name shadowing.  Every WSL call has GNU
``timeout --kill-after`` plus a Python timeout that explicitly kills its child.

Route B first solves the exact sparse QQ coefficient system on modular supports
seen in Route A, then tries common cofactor total-degree bounds 0 through 3
under a 25,000-column cap.  Its result passes the same exact expansion gate.

Manifest recipes are reimplemented here rather than imported from
``kill_certificate_tools.py``.  J6 systems are replayed from
``j6_msolve_results.json`` and cross-checked against ``alt_hunt_results.json``.
The status log is merged by kill_id: this invocation replaces only its own
targets and preserves all unrelated entries and top-level fields.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess
import time
from typing import Any, Iterable

import sympy as sp
from sympy.polys.matrices import DomainMatrix

ROOT = Path(__file__).resolve().parent
MANIFEST_PATH = ROOT / "kill_manifest.json"
TIMEOUT_IDS_PATH = ROOT / "timeout_ids.txt"
J6_RESULTS_PATH = ROOT / "j6_msolve_results.json"
ALT_HUNT_PATH = ROOT / "alt_hunt_results.json"
Q_COEFFS = (2048, -512, 320, -240, 195)
WSL_PREFIX = ("wsl.exe", "-d", "Ubuntu", "--", "bash", "-lc")
LARGE_PRIMES = (
    2147483647, 2147483629, 2147483587, 2147483579,
    2147483563, 2147483549, 2147483543, 2147483497,
    2147483489, 2147483477, 2147483423, 2147483399,
    2147483353, 2147483323, 2147483269, 2147483249,
)
PROJECT_PRIMES = (100019, 10009, 10007)
MAX_LINEAR_COLUMNS = 25_000
MAX_LINEAR_DEGREE = 3


def qpoly(v: sp.Symbol) -> sp.Expr:
    return sum(c * v ** (4 - i) for i, c in enumerate(Q_COEFFS))


def primitive(expr: sp.Expr, variables: list[sp.Symbol]) -> sp.Expr:
    """Canonical primitive integer associate, matching the existing producer."""
    poly = sp.Poly(sp.expand(sp.cancel(expr)), *variables, domain=sp.QQ)
    if poly.is_zero:
        return sp.Integer(0)
    denominator = 1
    for coefficient in poly.coeffs():
        denominator = sp.ilcm(denominator, int(coefficient.q))
    numerators = [abs(int(c * denominator)) for c in poly.coeffs() if c]
    content = math.gcd(*numerators) if numerators else 1
    result = sp.expand(poly.as_expr() * sp.Rational(denominator, content))
    if sp.Poly(result, *variables, domain=sp.ZZ).LC() < 0:
        result = -result
    return sp.expand(result)


def poly_json(expr: sp.Expr, variables: list[sp.Symbol]) -> dict[str, Any]:
    poly = sp.Poly(sp.expand(expr), *variables, domain=sp.QQ)
    return {"terms": [
        {"coefficient": {"numerator": str(coefficient.p),
                          "denominator": str(coefficient.q)},
         "powers": list(powers)}
        for powers, coefficient in poly.terms()
    ]}


def expr_string(expr: sp.Expr) -> str:
    return sp.sstr(sp.expand(expr))


def ring_vars(expressions: Iterable[sp.Expr]) -> list[sp.Symbol]:
    symbols: set[sp.Symbol] = set()
    for expression in expressions:
        symbols.update(sp.sympify(expression).free_symbols)
    return sorted(symbols, key=lambda symbol: symbol.name)


def reduce_roots(expr: sp.Expr, roots: list[sp.Symbol]) -> sp.Expr:
    result = sp.expand(expr)
    for root in roots:
        result = sp.rem(sp.Poly(result, root), sp.Poly(qpoly(root), root)).as_expr()
    return sp.expand(result)


def generic_poly(prefix: str, degree: int, y: sp.Symbol) -> sp.Expr:
    coefficients = sp.symbols(f"{prefix}0:{degree + 1}")
    return sum(coefficients[index] * y**index for index in range(degree + 1))


# Recipe reconstruction mirrors kill_certificate_tools.py, but remains local.
def material_sys4(name: str) -> tuple[list[sp.Expr], dict[str, Any]]:
    match = re.fullmatch(r"sub2T2_a(\d+)_b(\d{4})_dd2(-?inf|\d+)_dsig(\d+)", name)
    if not match:
        raise ValueError(f"bad System-4 name {name}")
    a, btxt, d2txt, dsig = int(match[1]), match[2], match[3], int(match[4])
    bsum = sum(map(int, btxt))
    y, root = sp.symbols("y r")
    gamma, w = sp.symbols("gamma w")
    roots = [root] if bsum else []
    generic_g = generic_poly("g", dsig - 2 * bsum, y)
    polys = {
        "d2": sp.Integer(0) if d2txt == "-inf" else generic_poly("a", int(d2txt), y),
        "d1": sp.Integer(0),
        "sigma": sp.expand((y - root) ** (2 * bsum) * generic_g) if roots else generic_g,
        "e": sp.expand(gamma * (y + 1)**a * ((y - root)**bsum if roots else 1)),
    }
    import convolution_descent as cd
    import convolution_elim as ce
    ansatz = cd.build_ansatz(d2=polys["d2"], d1=0, e=polys["e"],
                             sigma=polys["sigma"], parameters=tuple(roots))
    engine = cd.ConvolutionDescent(ansatz, c=ce.DEFAULT_C)
    coefficients, targets = [], []
    for degree in range(260, 195, -1):
        coefficient = reduce_roots(engine.master_coefficient(degree), roots)
        if coefficient != 0:
            coefficients.append(coefficient)
            targets.append(degree)
        if len(coefficients) == 8:
            break
    saturation = [gamma, sp.LC(sp.Poly(generic_g, y))]
    if roots:
        saturation.append(reduce_roots(generic_g.subs(y, root), roots))
    members = coefficients + [reduce_roots(w * sp.prod(saturation) - 1, roots)]
    members += [qpoly(item) for item in roots]
    material = {
        "identity": "f31_master",
        "polynomials": {key: expr_string(value) for key, value in polys.items()},
        "targets": targets,
        "root_variables": [str(item) for item in roots],
        "saturation_factors": [expr_string(item) for item in saturation],
        "source_name": name,
    }
    return members, material


def material_a8(name: str) -> tuple[list[sp.Expr], dict[str, Any]]:
    import convolution_descent as cd
    import convolution_elim as ce
    data = json.loads((ROOT / "batch_convolution_sub2.json").read_text())
    wanted = None
    for state in data["states"]:
        candidate = f"a8_dd2{state['deg_d2']}_dd1{state['deg_d1']}_dsig{state['deg_sigma']}"
        if candidate == name:
            wanted = state
            break
    if wanted is None:
        raise KeyError(name)
    y = cd.y
    gamma, w = sp.symbols("gamma w")
    d2 = sp.Integer(0) if wanted["d2_zero"] else generic_poly("a", int(wanted["deg_d2"]), y)
    d1 = generic_poly("b", int(wanted["deg_d1"]), y)
    sigma = generic_poly("s", int(wanted["deg_sigma"]), y)
    e = gamma * (y + 1)**8
    polys = {"d2": d2, "d1": d1, "sigma": sigma, "e": e}
    ansatz = cd.build_ansatz(d2=d2, d1=d1, e=e, sigma=sigma, parameters=(gamma,))
    start = int(wanted["gauge_detail"]["start"])
    engine = ce.HighCoefficientEngine(ansatz, start_degree=start,
                                      target_count=40, c=ce.DEFAULT_C)
    coefficients, targets = [], []
    for degree in range(start, start - 40, -1):
        coefficient = sp.expand(engine.master_coefficient(degree))
        if coefficient != 0:
            coefficients.append(coefficient)
            targets.append(degree)
        if len(coefficients) == 16:
            break
    return coefficients + [w * gamma - 1], {
        "identity": "f31_master",
        "polynomials": {key: expr_string(value) for key, value in polys.items()},
        "targets": targets, "root_variables": [],
        "saturation_factors": ["gamma"], "source_name": name,
    }


def material_sys3(name: str) -> tuple[list[sp.Expr], dict[str, Any]]:
    import phase_f2_scale as f2
    bid, support_text = name.rsplit("_sup", 1)
    support = int(support_text)
    narrowed = {row["key"] for row in json.loads(
        (ROOT / "phase_f2_scale.json").read_text())["alt_states"]
        if str(row.get("verdict", "")).startswith("NARROWED")}
    matches = [target for target in f2.load_targets()
               if target["bid"] == bid and int(target["support"]) == support
               and f"{target['bid']}#sup{target['support']}#idx{target['idx']}" in narrowed]
    if len(matches) != 1:
        raise RuntimeError(f"expected one target for {name}, found {len(matches)}")
    target = matches[0]
    degrees = target["degs"]
    drop_d1, drop_sigma = target["branch"] == "T2", target["sz"]
    total_degree = f2.total_deg(degrees, drop_d1, drop_sigma)
    factors, roots, scalars, d2_coefficients = f2.reconstruct(
        target["a"], target["b"], target["split"], target["branch"],
        degrees, drop_d1, drop_sigma)
    depth = min(int(target["depth"]), 12)
    tower = f2.h0_top(factors, tuple(item or 0 for item in degrees),
                      total_degree, depth, drop_d1=drop_d1, drop_sig=drop_sigma)
    reducer = f2.make_reducer(roots)
    coefficients = [reducer(item) for item in tower if item != 0]
    saturation = list(scalars)
    if d2_coefficients is not None and f2.d2_in_window(
            degrees, total_degree, depth, drop_d1, drop_sigma):
        saturation.append(d2_coefficients[-1])
    w = sp.Symbol("w")
    members = coefficients + [w * sp.prod(saturation) - 1]
    members += [f2.qpoly(root) for root in roots]
    if len(roots) == 2:
        members.append(sp.Symbol("wd") * (roots[0] - roots[1]) - 1)
    polynomials = {key: (value[0] if isinstance(value, tuple) else value)
                   for key, value in factors.items()}
    material = {
        "identity": "h0_tower",
        "polynomials": {key: expr_string(value) for key, value in polynomials.items()},
        "degrees": list(degrees), "TD": total_degree, "depth": depth,
        "drop_d1": drop_d1, "drop_sigma": drop_sigma,
        "root_variables": [str(root) for root in roots],
        "saturation_factors": [expr_string(item) for item in saturation],
        "source_name": name,
        "source_key": f"{target['bid']}#sup{support}#idx{target['idx']}",
    }
    return members, material


def material_d2(branch: str, degree: int, depth: int) -> tuple[list[sp.Expr], dict[str, Any]]:
    import d2_threshold as dt
    d2, sigma, e, roots, d2_coefficients = dt.build_state(branch, degree)
    reducer = dt.reducer(roots)
    coefficients = [reducer(item) for item in dt.h0_top(d2, degree, sigma, e, depth)
                    if item != 0]
    saturation = [dt.S, dt.E, d2_coefficients[-1]]
    members = coefficients + [dt.w * sp.prod(saturation) - 1]
    members += [dt.qpoly(root) for root in roots]
    if len(roots) == 2:
        members.append(dt.wd * (roots[0] - roots[1]) - 1)
    material = {
        "identity": "h0_tower",
        "polynomials": {"d2": expr_string(d2), "sigma": expr_string(sigma),
                        "e": expr_string(e)},
        "degrees": [degree, None, 12, 15], "TD": 60, "depth": depth,
        "drop_d1": True, "drop_sigma": False,
        "root_variables": [str(root) for root in roots],
        "saturation_factors": [expr_string(item) for item in saturation],
        "source_name": branch,
    }
    return members, material


def material_phase(key: str, depth: int) -> tuple[list[sp.Expr], dict[str, Any]]:
    import phase_f2_sub2 as f2
    import convolution_descent as cd
    matches = []
    for cell, case, state, index, _maximal, pdelta in f2.load_targets(
            f2.TARGET_CELLS, max_defect=1):
        if f"{cell}#state{index}" == key:
            matches.append((case, state, pdelta))
    if len(matches) != 1:
        raise RuntimeError(f"expected one phase target for {key}")
    case, state, pdelta = matches[0]
    combo, why, _ = f2.unique_split(case, state, pdelta)
    if combo is None:
        raise RuntimeError(why)
    polys, scalars, marked, mode, cofactors = f2.reconstruct(case, state, combo, pdelta)
    parameters = ((f2.r,) if marked is not None else ()) + tuple(cofactors)
    ansatz = cd.build_ansatz(d2=polys["d2"], d1=polys["d1"], e=polys["e"],
                             sigma=polys["sigma"], parameters=parameters)
    engine = cd.ConvolutionDescent(ansatz, c=f2.C_VAL)
    top = f2.engine_top(engine)
    coefficients, targets = [], []
    for offset in range(depth):
        coefficient = f2.redq(engine.master_coefficient(top - offset), marked)
        if coefficient != 0:
            coefficients.append(coefficient)
            targets.append(top - offset)
    members = coefficients + [f2.w * sp.prod(scalars) - 1]
    roots = [f2.r] if marked is not None else []
    if roots:
        members.append(f2.QR_EXPR)
    material = {
        "identity": "f31_master",
        "polynomials": {name: expr_string(value) for name, value in polys.items()},
        "targets": targets, "root_variables": [str(root) for root in roots],
        "saturation_factors": [expr_string(item) for item in scalars],
        "source_key": key,
        "state_degrees": [state["deg_d1"], state["deg_sigma"],
                          state["deg_d2"], state["deg_e"]],
        "d2_mode": mode, "cofactor_variables": [str(item) for item in cofactors],
    }
    return members, material


def resolve_manifest(entry: dict[str, Any]) -> tuple[list[sp.Expr], dict[str, Any]]:
    recipe, builder = entry["recipe"], entry["recipe"]["builder"]
    if builder == "harvest_sys4":
        return material_sys4(recipe["name"])
    if builder == "harvest_a8":
        return material_a8(recipe["name"])
    if builder == "harvest_sys3":
        return material_sys3(recipe["name"])
    if builder == "d2_threshold":
        return material_d2(recipe["branch"], int(recipe["degree_d2"]), int(recipe["depth"]))
    if builder == "phase_f2_sub2":
        return material_phase(recipe["state_key"], int(recipe["depth"]))
    if builder == "blowup_case":
        case = recipe["case"]
        if case == "a12_b1110_T2_d6":
            members, material = material_d2("a12_b1110_T2", 6, 12)
        elif case == "a11_b1111_T1_17":
            members, material = material_sys3("a11_b1111_T1_sup17")
        else:
            keys = {"sub2_s14": "sub2:a9_b1000_T1_sz0_dz0_gz-#state14",
                    "sub2_s38": "sub2:a9_b1000_T1_sz0_dz0_gz-#state38",
                    "sub2_s94": "sub2:a9_b1000_T1_sz0_dz0_gz-#state94"}
            members, material = material_phase(keys[case], 6)
        material["archived_ms"] = {
            "path": recipe.get("archived_ms"), "reused": False,
            "note": "recipe regenerated locally; archived WSL input not required"}
        return members, material
    raise ValueError(f"unknown builder {builder}")


def parse_expr(text: str, names: Iterable[str]) -> sp.Expr:
    local = {name: sp.Symbol(name) for name in names}
    return sp.sympify(text.replace("^", "**"), locals=local)


def resolve_j6(entry: dict[str, Any]) -> tuple[list[sp.Expr], dict[str, Any]]:
    key = entry["recipe"]["state_key"]
    rows = json.loads(J6_RESULTS_PATH.read_text(encoding="utf-8"))["results"]
    matches = [row for row in rows if row["key"] == key]
    if len(matches) != 1 or matches[0].get("verdict") != "KILLED":
        raise RuntimeError(f"expected one KILLED J6 result for {key}")
    row = matches[0]
    depth, names = int(row["kill_depth"]), row["ring_vars"]
    coefficient_rows = row["gens"][:depth]
    if len(coefficient_rows) != depth:
        raise RuntimeError(f"J6 result for {key} lacks depth-{depth} generators")
    degrees = [int(item["degree"]) for item in coefficient_rows]
    expected = list(range(int(row["top_degree"]), int(row["top_degree"]) - depth, -1))
    if degrees != expected:
        raise RuntimeError(f"J6 generator degrees {degrees} do not match {expected}")
    coefficients = [parse_expr(item["coefficient"], names) for item in coefficient_rows]
    saturation = parse_expr(row["saturation"], names)
    class_relations = [parse_expr(text, names) for text in row["class_relations"]]
    states = json.loads(ALT_HUNT_PATH.read_text(encoding="utf-8"))["states"]
    state_matches = [state for state in states if state.get("key") == key]
    if len(state_matches) != 1:
        raise RuntimeError(f"expected one ALT_HUNT state for {key}")
    split_matches = [split for split in state_matches[0].get("splits", [])
                     if split.get("combo") == row.get("combo")]
    if len(split_matches) != 1:
        raise RuntimeError(f"J6 combo for {key} does not match ALT_HUNT")
    members = coefficients + [saturation] + class_relations
    material = {
        "identity": "recorded_j6_system",
        "polynomials": split_matches[0]["polys"], "targets": degrees,
        "root_variables": [], "saturation_factors": ["E", "X"],
        "recorded_saturation": row["saturation"],
        "recorded_class_relations": row["class_relations"],
        "source_key": key,
        "source_files": ["alt_hunt_results.json", "j6_msolve_results.json"],
        "kill_depth": depth, "combo": row["combo"],
    }
    return members, material


def resolve(entry: dict[str, Any]) -> tuple[list[sp.Expr], dict[str, Any]]:
    return resolve_j6(entry) if entry["category"] == "j6_msolve" else resolve_manifest(entry)


# Singular emission uses only integer coefficients and synthetic multi-letter variables.
def integer_poly_text(expr: sp.Expr, variables: list[sp.Symbol], names: list[str]) -> str:
    poly = sp.Poly(sp.expand(expr), *variables, domain=sp.ZZ)
    chunks = []
    for powers, coefficient in poly.terms():
        monomial = "*".join(name + (f"^{power}" if power != 1 else "")
                            for name, power in zip(names, powers) if power)
        value = int(coefficient)
        if monomial:
            term = monomial if value == 1 else "-" + monomial if value == -1 \
                else f"{value}*{monomial}"
        else:
            term = str(value)
        chunks.append(term)
    return "+".join(chunks).replace("+-", "-") or "0"


def emit_modular_program(generators: list[sp.Expr], variables: list[sp.Symbol],
                         prime: int) -> str:
    names = [f"V{index}" for index in range(len(variables))]
    lines = [f"ring R = {prime},({','.join(names)}),dp;"]
    for index, generator in enumerate(generators):
        text = integer_poly_text(generator, variables, names)
        if "/" in text or "**" in text:
            raise AssertionError("non-integer or Python power syntax in Singular input")
        lines.append(f"poly GEN{index} = {text};")
    lines.append("ideal I = " + ",".join(f"GEN{i}" for i in range(len(generators))) + ";")
    lines += ["ideal U = 1;", "ideal G = std(I);", "matrix T = lift(I,G);",
              "matrix A = lift(G,U);", "matrix L = T*A;"]
    for index in range(len(generators)):
        lines += [f'"@@COF_BEGIN_{index}";', f"L[{index + 1},1];",
                  f'"@@COF_END_{index}";']
    lines += ['"@@DONE";', "quit;"]
    return "\n".join(lines) + "\n"


def kill_process(process: subprocess.Popen[str]) -> None:
    try:
        process.kill()
    except OSError:
        pass
    try:
        process.communicate(timeout=5)
    except Exception:
        pass


def run_singular_program(program: str, timeout: float) -> dict[str, Any]:
    command = WSL_PREFIX + (
        f'cd "$HOME" && timeout --kill-after=5s {timeout:g}s Singular -q',)
    started = time.monotonic()
    try:
        process = subprocess.Popen(command, stdin=subprocess.PIPE,
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                   text=True, encoding="utf-8")
    except Exception as error:
        return {"ok": False, "wall_seconds": round(time.monotonic() - started, 3),
                "reason": f"Singular launch failed: {type(error).__name__}: {error}"}
    try:
        stdout, stderr = process.communicate(program, timeout=timeout + 10)
    except subprocess.TimeoutExpired:
        kill_process(process)
        return {"ok": False, "wall_seconds": round(time.monotonic() - started, 3),
                "reason": f"Singular subprocess timed out and was killed after {timeout:g}s"}
    wall = round(time.monotonic() - started, 3)
    combined = (stdout or "") + "\n" + (stderr or "")
    if process.returncode != 0 or "@@DONE" not in combined:
        tail = combined.replace("\x00", "").strip()[-800:]
        return {"ok": False, "wall_seconds": wall,
                "reason": f"Singular failed (exit {process.returncode}): {tail}"}
    return {"ok": True, "wall_seconds": wall, "output": combined}


def parse_modular_cofactors(output: str, count: int, variable_count: int,
                            prime: int) -> list[dict[tuple[int, ...], int]]:
    symbols = [sp.Symbol(f"V{index}") for index in range(variable_count)]
    local = {str(symbol): symbol for symbol in symbols}
    result = []
    for index in range(count):
        match = re.search(rf"@@COF_BEGIN_{index}\s*(.*?)\s*@@COF_END_{index}",
                          output, re.S)
        if not match:
            raise ValueError(f"missing modular cofactor {index}")
        lines = [line.strip() for line in match.group(1).splitlines()
                 if line.strip() and not line.lstrip().startswith("//")]
        expression = sp.sympify(("".join(lines) or "0").replace("^", "**"), locals=local)
        poly = sp.Poly(sp.expand(expression), *symbols, domain=sp.ZZ)
        result.append({tuple(powers): int(coefficient) % prime
                       for powers, coefficient in poly.terms()
                       if int(coefficient) % prime})
    return result


def modular_lift(generators: list[sp.Expr], variables: list[sp.Symbol], prime: int,
                 timeout: float) -> dict[str, Any]:
    run = run_singular_program(emit_modular_program(generators, variables, prime), timeout)
    if not run["ok"]:
        return run
    try:
        cofactors = parse_modular_cofactors(run["output"], len(generators),
                                            len(variables), prime)
    except Exception as error:
        return {"ok": False, "wall_seconds": run["wall_seconds"],
                "reason": f"modular output parse failed: {type(error).__name__}: {error}"}
    return {"ok": True, "wall_seconds": run["wall_seconds"], "cofactors": cofactors}


Support = tuple[tuple[tuple[int, ...], ...], ...]


def support_signature(cofactors: list[dict[tuple[int, ...], int]]) -> Support:
    return tuple(tuple(sorted(poly)) for poly in cofactors)


def crt_pair(residue: int, modulus: int, new: int, prime: int) -> int:
    step = ((new - residue) % prime) * pow(modulus, -1, prime) % prime
    return residue + modulus * step


def rational_reconstruct(residue: int, modulus: int) -> sp.Rational | None:
    """Symmetric rational reconstruction with sqrt(modulus/2) bounds."""
    residue %= modulus
    if residue == 0:
        return sp.Rational(0)
    bound = math.isqrt(modulus // 2)
    old_r, new_r, old_t, new_t = modulus, residue, 0, 1
    while new_r and abs(new_r) > bound:
        quotient = old_r // new_r
        old_r, new_r = new_r, old_r - quotient * new_r
        old_t, new_t = new_t, old_t - quotient * new_t
    numerator, denominator = new_r, new_t
    if denominator < 0:
        numerator, denominator = -numerator, -denominator
    if (denominator == 0 or denominator > bound or abs(numerator) > bound
            or math.gcd(abs(numerator), denominator) != 1
            or (residue * denominator - numerator) % modulus):
        return None
    return sp.Rational(numerator, denominator)


def update_crt_group(group: dict[str, Any], modular: list[dict[tuple[int, ...], int]],
                     prime: int) -> None:
    old_modulus = int(group["modulus"])
    if old_modulus == 1:
        group["residues"] = [dict(poly) for poly in modular]
    else:
        for accumulated, new_poly in zip(group["residues"], modular):
            for monomial in accumulated:
                accumulated[monomial] = crt_pair(accumulated[monomial], old_modulus,
                                                  new_poly.get(monomial, 0), prime)
    group["modulus"] = old_modulus * prime
    group["primes"].append(prime)


def reconstruct_group(group: dict[str, Any], variables: list[sp.Symbol]) -> list[sp.Expr] | None:
    output = []
    for polynomial in group["residues"]:
        expression = sp.Integer(0)
        for powers, residue in polynomial.items():
            coefficient = rational_reconstruct(residue, int(group["modulus"]))
            if coefficient is None:
                return None
            monomial = coefficient
            for variable, power in zip(variables, powers):
                if power:
                    monomial *= variable**power
            expression += monomial
        output.append(sp.expand(expression))
    return output


def verify_identity(cofactors: list[sp.Expr], generators: list[sp.Expr],
                    variables: list[sp.Symbol]) -> tuple[bool, str]:
    if len(cofactors) != len(generators):
        return False, "generator/cofactor length mismatch"
    try:
        residual = sp.Poly(-1, *variables, domain=sp.QQ)
        for cofactor, generator in zip(cofactors, generators):
            residual += sp.Poly(sp.expand(cofactor * generator), *variables, domain=sp.QQ)
    except Exception as error:
        return False, f"exact expansion raised {type(error).__name__}: {error}"
    if residual.is_zero:
        return True, "exact SymPy expansion verified"
    return False, f"exact expansion residual has {len(residual.terms())} nonzero terms"


def prime_sequence(count: int) -> list[int]:
    if count <= 0:
        raise ValueError("--primes must be positive")
    primes = list(LARGE_PRIMES) + list(PROJECT_PRIMES)
    candidate = LARGE_PRIMES[-1]
    while len(primes) < count:
        candidate = int(sp.prevprime(candidate))
        if candidate not in primes:
            primes.insert(len(primes) - len(PROJECT_PRIMES), candidate)
    return primes[:count]


def route_a(generators: list[sp.Expr], variables: list[sp.Symbol], prime_count: int,
            timeout: float, quiet: bool
            ) -> tuple[list[sp.Expr] | None, dict[str, Any], list[dict[str, Any]]]:
    groups: dict[Support, dict[str, Any]] = {}
    attempts, total_wall = [], 0.0
    for prime in prime_sequence(prime_count):
        if not quiet:
            print(f"    Route A prime {prime}", flush=True)
        result = modular_lift(generators, variables, prime, timeout)
        total_wall += float(result.get("wall_seconds", 0.0))
        attempt = {"prime": prime, "wall_seconds": result.get("wall_seconds"),
                   "status": "MODULAR-LIFT-FOUND" if result["ok"] else "FAILED"}
        if not result["ok"]:
            attempt["reason"] = result.get("reason", "modular lift failed")
            attempts.append(attempt)
            continue
        modular = result["cofactors"]
        signature = support_signature(modular)
        group = groups.setdefault(signature, {"signature": signature, "modulus": 1,
                                              "residues": [], "primes": []})
        update_crt_group(group, modular, prime)
        attempt.update({"support_group_size": len(group["primes"]),
                        "support_terms": sum(len(poly) for poly in signature)})
        attempts.append(attempt)
        candidate = reconstruct_group(group, variables)
        if candidate is None:
            attempt["reconstruction"] = "not yet unique within rational bound"
            continue
        verified, detail = verify_identity(candidate, generators, variables)
        attempt["reconstruction"] = detail
        if verified:
            return candidate, {
                "status": "CERTIFICATE-FOUND",
                "lift_method": "modular-Singular-CRT-rational-reconstruction",
                "lift_wall_seconds": round(total_wall, 3),
                "primes_used": list(group["primes"]),
                "crt_modulus_bits": int(group["modulus"]).bit_length(),
                "attempts": attempts}, list(groups.values())
    return None, {
        "status": "NOT-RECONSTRUCTED",
        "lift_method": "modular-Singular-CRT-rational-reconstruction",
        "lift_wall_seconds": round(total_wall, 3),
        "reason": f"no exact characteristic-zero identity after {prime_count} prime calls",
        "attempts": attempts}, list(groups.values())


# Exact sparse coefficient linear algebra.
def compositions(total: int, length: int) -> Iterable[tuple[int, ...]]:
    if length == 0:
        if total == 0:
            yield ()
        return
    if length == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in compositions(total - first, length - 1):
            yield (first,) + tail


def monomials_through_degree(variable_count: int, degree: int) -> tuple[tuple[int, ...], ...]:
    return tuple(powers for total in range(degree + 1)
                 for powers in compositions(total, variable_count))


def solve_on_support(generators: list[sp.Expr], variables: list[sp.Symbol],
                     supports: Support) -> tuple[list[sp.Expr] | None, dict[str, Any]]:
    columns = [(index, powers) for index, support in enumerate(supports)
               for powers in support]
    if not columns:
        return None, {"reason": "empty cofactor support", "columns": 0, "rows": 0}
    if len(columns) > MAX_LINEAR_COLUMNS:
        return None, {"reason": f"support exceeds {MAX_LINEAR_COLUMNS}-column cap",
                      "columns": len(columns), "rows": 0}
    generator_terms = [sp.Poly(item, *variables, domain=sp.QQ).terms()
                       for item in generators]
    equations: dict[tuple[int, ...], dict[int, sp.Rational]] = {}
    for column, (generator_index, cofactor_powers) in enumerate(columns):
        for generator_powers, coefficient in generator_terms[generator_index]:
            product = tuple(a + b for a, b in zip(cofactor_powers, generator_powers))
            row = equations.setdefault(product, {})
            value = row.get(column, sp.Rational(0)) + coefficient
            if value:
                row[column] = value
            else:
                row.pop(column, None)
    constant = (0,) * len(variables)
    equations.setdefault(constant, {})
    monomials = sorted(equations, reverse=True)
    augmented = {}
    for row_index, monomial in enumerate(monomials):
        row = dict(equations[monomial])
        if monomial == constant:
            row[len(columns)] = sp.Integer(1)
        if row:
            augmented[row_index] = row
    try:
        matrix = DomainMatrix.from_dict_sympy(len(monomials), len(columns) + 1, augmented)
        reduced, pivots = matrix.rref()
    except Exception as error:
        return None, {"reason": f"exact sparse rref failed: {type(error).__name__}: {error}",
                      "columns": len(columns), "rows": len(monomials)}
    if len(columns) in pivots:
        return None, {"reason": "bounded coefficient system is inconsistent over QQ",
                      "columns": len(columns), "rows": len(monomials),
                      "rank": len(pivots) - 1}
    entries = reduced.to_dok()
    solution = [sp.Rational(0)] * len(columns)
    for row_index, pivot in enumerate(pivots):
        if pivot < len(columns):
            value = entries.get((row_index, len(columns)), reduced.domain.zero)
            solution[pivot] = reduced.domain.to_sympy(value)
    cofactors = [sp.Integer(0) for _ in generators]
    for value, (generator_index, powers) in zip(solution, columns):
        monomial = value
        for variable, power in zip(variables, powers):
            if power:
                monomial *= variable**power
        cofactors[generator_index] += monomial
    cofactors = [sp.expand(item) for item in cofactors]
    verified, detail = verify_identity(cofactors, generators, variables)
    metadata = {"columns": len(columns), "rows": len(monomials),
                "rank": len(pivots), "verification": detail}
    return (cofactors if verified else None), metadata


def route_b(generators: list[sp.Expr], variables: list[sp.Symbol],
            modular_groups: list[dict[str, Any]], quiet: bool
            ) -> tuple[list[sp.Expr] | None, dict[str, Any]]:
    started, attempts = time.monotonic(), []
    groups = sorted(modular_groups,
                    key=lambda group: (len(group["primes"]),
                                       int(group["modulus"]).bit_length()), reverse=True)
    seen: set[Support] = set()
    for group in groups:
        support = group["signature"]
        if support in seen:
            continue
        seen.add(support)
        if not quiet:
            print(f"    Route B modular support ({sum(len(x) for x in support)} columns)",
                  flush=True)
        cofactors, detail = solve_on_support(generators, variables, support)
        detail.update({"ansatz": "modular-support", "support_primes": group["primes"]})
        attempts.append(detail)
        if cofactors is not None:
            return cofactors, {"status": "CERTIFICATE-FOUND",
                "lift_method": "degree-bounded-exact-linear-algebra",
                "linear_wall_seconds": round(time.monotonic() - started, 3),
                "attempts": attempts}
    for degree in range(MAX_LINEAR_DEGREE + 1):
        common = monomials_through_degree(len(variables), degree)
        support = tuple(common for _ in generators)
        columns = len(common) * len(generators)
        if not quiet:
            print(f"    Route B total degree {degree} ({columns} columns)", flush=True)
        if columns > MAX_LINEAR_COLUMNS:
            attempts.append({"ansatz": "common-total-degree", "degree": degree,
                             "columns": columns, "reason": "exceeds safety cap"})
            break
        cofactors, detail = solve_on_support(generators, variables, support)
        detail.update({"ansatz": "common-total-degree", "degree": degree})
        attempts.append(detail)
        if cofactors is not None:
            return cofactors, {"status": "CERTIFICATE-FOUND",
                "lift_method": "degree-bounded-exact-linear-algebra",
                "linear_wall_seconds": round(time.monotonic() - started, 3),
                "attempts": attempts}
    return None, {"status": "NOT-YET-CERTIFICATED",
        "lift_method": "degree-bounded-exact-linear-algebra",
        "linear_wall_seconds": round(time.monotonic() - started, 3),
        "reason": ("exact sparse linear algebra found no identity on observed modular "
                   f"supports or through total degree {MAX_LINEAR_DEGREE}"),
        "attempts": attempts}


def safe_name(kill_id: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "__", kill_id) + ".json"


def j6_entries() -> list[dict[str, Any]]:
    rows = json.loads(J6_RESULTS_PATH.read_text(encoding="utf-8"))["results"]
    entries = [{"id": f"j6_msolve:{row['key']}", "category": "j6_msolve",
                "recipe": {"builder": "j6_recorded", "state_key": row["key"],
                           "kill_depth": int(row["kill_depth"])}} for row in rows]
    if len(entries) != 4:
        raise RuntimeError(f"expected four J6 entries, found {len(entries)}")
    return entries


def target_entries() -> list[dict[str, Any]]:
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8-sig"))
    wanted = {line.strip() for line in TIMEOUT_IDS_PATH.read_text().splitlines()
              if line.strip()}
    entries = [entry for entry in manifest["entries"] if entry["id"] in wanted]
    found = {entry["id"] for entry in entries}
    if found != wanted:
        raise RuntimeError("timeout IDs absent from manifest: " + ", ".join(sorted(wanted - found)))
    if len(entries) != 15:
        raise RuntimeError(f"expected 15 timeout IDs, found {len(entries)}")
    return entries + j6_entries()


def generator_digest(items: list[dict[str, Any]]) -> str:
    packed = json.dumps(items, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(packed).hexdigest()


def status_summary(result: dict[str, Any]) -> dict[str, Any]:
    keys = ("kill_id", "category", "status", "reason", "lift_method",
            "lift_wall_seconds", "linear_wall_seconds", "attempt_wall_seconds")
    return {key: result[key] for key in keys if result.get(key) is not None}


def merge_status_log(path: Path, updates: list[dict[str, Any]], timeout: float) -> None:
    if path.exists():
        try:
            log = json.loads(path.read_text(encoding="utf-8"))
        except Exception as error:
            raise RuntimeError(f"refusing to clobber unreadable status log: {error}") from error
        if not isinstance(log, dict) or not isinstance(log.get("entries", []), list):
            raise RuntimeError("refusing to clobber malformed status log")
    else:
        log = {"schema": "d2-kill-status-log-v1", "manifest": MANIFEST_PATH.name,
               "timeout_seconds": timeout, "entries": []}
    update_map = {item["kill_id"]: item for item in updates}
    merged, replaced = [], set()
    for old in log.get("entries", []):
        kill_id = old.get("kill_id") if isinstance(old, dict) else None
        if kill_id in update_map:
            merged.append(update_map[kill_id])
            replaced.add(kill_id)
        else:
            merged.append(old)
    merged += [item for item in updates if item["kill_id"] not in replaced]
    log["entries"] = merged
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(log, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def produce_one(entry: dict[str, Any], prime_count: int, timeout: float,
                quiet: bool) -> dict[str, Any]:
    started = time.monotonic()
    base = {"schema": "d2-kill-certificate-v1", "kill_id": entry["id"],
            "category": entry["category"], "manifest_recipe": entry["recipe"],
            "timeout_seconds": timeout}
    try:
        raw, material = resolve(entry)
        variables = ring_vars(raw)
        generators = [primitive(item, variables) for item in raw if item != 0]
        if not variables or not generators:
            raise RuntimeError("resolved system has no variables or generators")
        for generator in generators:
            sp.Poly(generator, *variables, domain=sp.ZZ)
        generator_json = [poly_json(item, variables) for item in generators]
        common = {**base, "variable_order": [str(item) for item in variables],
            "generator_normalization": "primitive integer associate with positive leading coefficient",
            "generating_recipe": material, "generators": generator_json,
            "generator_sha256": generator_digest(generator_json)}
        cofactors, route_a_meta, groups = route_a(generators, variables, prime_count,
                                                  timeout, quiet)
        route_b_meta = None
        if cofactors is None:
            cofactors, route_b_meta = route_b(generators, variables, groups, quiet)
        if cofactors is not None:
            verified, detail = verify_identity(cofactors, generators, variables)
            if not verified:
                raise AssertionError("selected certificate failed final gate: " + detail)
            winner = route_a_meta if route_a_meta["status"] == "CERTIFICATE-FOUND" else route_b_meta
            result = {**common, "status": "CERTIFICATE-FOUND",
                "lift_method": winner["lift_method"],
                "lift_wall_seconds": route_a_meta.get("lift_wall_seconds"),
                "cofactors": [poly_json(item, variables) for item in cofactors],
                "exact_expansion_check": detail, "route_a": route_a_meta,
                "attempt_wall_seconds": round(time.monotonic() - started, 3)}
            if route_b_meta is not None:
                result["route_b"] = route_b_meta
                result["linear_wall_seconds"] = route_b_meta.get("linear_wall_seconds")
            return result
        assert route_b_meta is not None
        return {**common, "status": "NOT-YET-CERTIFICATED",
            "reason": route_b_meta["reason"], "lift_method": route_b_meta["lift_method"],
            "lift_wall_seconds": route_a_meta.get("lift_wall_seconds"),
            "linear_wall_seconds": route_b_meta.get("linear_wall_seconds"),
            "route_a": route_a_meta, "route_b": route_b_meta,
            "attempt_wall_seconds": round(time.monotonic() - started, 3)}
    except Exception as error:
        return {**base, "status": "NOT-YET-CERTIFICATED",
                "reason": f"recipe/build error: {type(error).__name__}: {error}",
                "attempt_wall_seconds": round(time.monotonic() - started, 3)}


def main() -> int:
    parser = argparse.ArgumentParser(description="Certify the 19 msolve-assisted kills")
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument("--all", action="store_true",
                           help="the 15 timeout IDs and all four J6 states")
    selection.add_argument("--id", action="append", default=[], metavar="KILL-ID",
                           help="one target; repeat for multiple targets")
    parser.add_argument("--primes", type=int, default=8, metavar="N",
                        help="Route-A prime calls per target (default: 8)")
    parser.add_argument("--timeout", type=float, metavar="SECONDS",
                        help="timeout per Singular call (default: manifest value)")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()
    if args.primes <= 0:
        parser.error("--primes must be positive")
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8-sig"))
    timeout = float(args.timeout if args.timeout is not None
                    else manifest.get("timeout_seconds", 300))
    if timeout <= 0:
        parser.error("--timeout must be positive")
    available = target_entries()
    by_id = {entry["id"]: entry for entry in available}
    for entry in available:
        if entry["category"] == "j6_msolve":
            by_id.setdefault(entry["recipe"]["state_key"], entry)
    if args.all:
        entries = available
    else:
        unknown = [kill_id for kill_id in args.id if kill_id not in by_id]
        if unknown:
            parser.error("unknown or non-target ID(s): " + ", ".join(unknown))
        entries, seen = [], set()
        for kill_id in args.id:
            entry = by_id[kill_id]
            if entry["id"] not in seen:
                entries.append(entry)
                seen.add(entry["id"])
    output_dir = ROOT / manifest.get("output_dir", "kill_certificates")
    output_dir.mkdir(parents=True, exist_ok=True)
    statuses = []
    for index, entry in enumerate(entries, 1):
        if not args.quiet:
            print(f"[{index}/{len(entries)}] {entry['id']}", flush=True)
        result = produce_one(entry, args.primes, timeout, args.quiet)
        (output_dir / safe_name(entry["id"])).write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        statuses.append(status_summary(result))
        if not args.quiet:
            message = result.get("reason", f"{result.get('lift_method')} in {result.get('attempt_wall_seconds')}s")
            print(f"  {result['status']}: {message}", flush=True)
    status_path = ROOT / manifest.get("status_log", "kill_certificates/status_log.json")
    merge_status_log(status_path, statuses, timeout)
    print(f"MSOLVE-ASSISTED PRODUCTION CENSUS {dict(Counter(x['status'] for x in statuses))}",
          flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
