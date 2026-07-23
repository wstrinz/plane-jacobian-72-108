#!/usr/bin/env python3
"""Exact source-linked checks for RESIDUE_LEMMAS.md.

All residue coefficients come from f31_graded.txt. The cone JSON files
select exponent supports and supply incidence only; no h_l coefficient is
copied into this verifier.
"""

from collections import Counter, defaultdict
import json
import re
from pathlib import Path

import sympy as sp
from sympy.polys.numberfields import galois_group

ROOT = Path(__file__).resolve().parent
WINDOWS = {"sub2": ROOT / "cascade_cones_qt.json",
           "sub1": ROOT / "cascade_cones_sub1_qt.json"}
d0, d2, d1, dm1, sigma = sp.symbols("d0 d2 d1 dm1 sigma")
D, X, S, E = sp.symbols("D X S E", nonzero=True)
y = sp.Symbol("y")
q = 2048*y**4 - 512*y**3 + 320*y**2 - 240*y + 195
VARS, LEADS = (d2, d1, sigma, dm1), (D, X, S, E)


def load_h():
    pattern = re.compile(r"h_(\d+)\s*\([^)]*\)\s*=\s*(.+)$")
    raw = {}
    for line in (ROOT / "f31_graded.txt").read_text(encoding="utf-8").splitlines():
        match = pattern.match(line.strip())
        if match:
            raw[int(match.group(1))] = sp.sympify(
                match.group(2), locals={"d0": d0, "d1": d1,
                                        "d2": d2, "dm1": dm1})
    assert sorted(raw) == list(range(8))
    return {level: sp.expand(expr.subs(d0, (sigma+d2**2)/4))
            for level, expr in raw.items()}


def source_terms(h):
    return {level: dict(sp.Poly(h[level], *VARS, domain=sp.QQ).terms())
            for level in (4, 5, 6)}


def survivors(path):
    data = json.loads(path.read_text(encoding="utf-8"))
    return [row for row in data["branches"] if row["status"] == "survives"]


def cell_id(row):
    return f"a{row['a_t']} b={''.join(map(str, row['b']))} {row['branch']}"


def tied_support(text):
    expr = sp.sympify(text, locals={"d2": d2, "d1": d1,
                                   "sigma": sigma, "e": dm1})
    parsed = sp.Poly(expr, *VARS, domain=sp.QQ).terms()
    assert len(parsed) == 1
    return parsed[0]


def inventory(rows, terms):
    freq, cells, depths = Counter(), defaultdict(set), defaultdict(Counter)
    for row in rows:
        for case in row["survivor_cases"]:
            for place in case["witness"]:
                for ob in place["obligations"]:
                    checked = []
                    for text in ob.get("tied", []):
                        support, json_coeff = tied_support(text)
                        assert support in terms[ob["level"]]
                        assert json_coeff == terms[ob["level"]][support]
                        checked.append(support)
                    key = (place["place"], int(ob["level"]), ob["kind"],
                           tuple(sorted(checked)))
                    freq[key] += 1
                    cells[key].add(cell_id(row))
                    depths[key][int(ob["depth"])] += 1
    return freq, cells, depths


def residue_shapes(inventories, terms):
    shapes = {(key[1], key[3]) for freq, _, _ in inventories.values()
              for key in freq if key[3] and key[2] in
              {"monomial_tie_rise", "identical_vanishing"}}
    ordered = sorted(shapes, key=lambda item: (-item[0], len(item[1]), item[1]))
    assert len(ordered) == 23
    equations = {}
    for number, (level, support) in enumerate(ordered, 1):
        equations[f"C{number:02}"] = sp.expand(sum(
            terms[level][monom] * sp.prod(v**n for v, n in zip(LEADS, monom))
            for monom in support))
    return ordered, equations


def Q(n, d=1):
    return sp.Rational(n, d)


def squarefree_part(value):
    value = int(value)
    return sp.prod(p for p, exponent in sp.factorint(value).items()
                   if exponent % 2)


def rational_witnesses():
    r73, r15 = Q(73, 4), Q(13797, 1952)
    u21, u22 = Q(-13797, 1792), Q(152, 511)
    return {
        "C01": {X: -6, S: -4, E: -1},
        "C02": {D: Q(3,14), X: 1, S: 1},
        "C03": {D: 1, X: 1, E: Q(-7,4)},
        "C04": {D: -2, X: -1, S: -2, E: -5},
        "C05": {X: 1, S: Q(-4,63), E: 1},
        "C06": {D: 2, S: -2, E: -6},
        "C07": {D: -6, X: -2, S: -1, E: -6},
        "C09": {D: Q(13,6), X: Q(13,6), S: Q(169,36)},
        "C10": {D: -6, X: -2, S: -4, E: 6},
        "C11": {X: 1, S: 10, E: 35},
        "C12": {X: r73**2, S: -r73**3},
        "C13": {D: Q(189,88), S: 1, E: 1},
        "C14": {D: 7, S: -21},
        "C15": {X: 1, S: 1, E: Q(511,80)},
        "C16": {D: -r15, X: r15},
        "C17": {D: Q(88,81), S: 1, E: 1},
        "C18": {X: 1, S: 1, E: Q(539,80)},
        "C19": {D: 1, S: Q(1,27), E: Q(1,27)},
        "C21": {D: u21, X: u21, E: u21**2},
        "C22": {D: 1/u22, X: 1/u22, S: Q(-4,3)/u22**2},
        "C23": {D: -54, X: -28, S: 100, E: 3730},
    }


def verify_constraints(equations):
    kills = {"C08", "C20"}
    witnesses = rational_witnesses()
    assert set(witnesses) == set(equations) - kills
    for name, point in witnesses.items():
        used = equations[name].free_symbols
        assert used <= set(point) and all(point[v] != 0 for v in used)
        assert sp.expand(equations[name].subs(point)) == 0

    # Derive the two torus quadratics from the source equations.
    r, quadratics = sp.Symbol("r"), {}
    for name in kills:
        reduced = equations[name].subs({D: 1, X: 1, E: r})
        quadratics[name] = sp.Poly(reduced, r, domain=sp.QQ).primitive()[1]
        assert quadratics[name].degree() == 2
    discriminants = {name: poly.discriminant()
                     for name, poly in quadratics.items()}
    assert {squarefree_part(value) for value in discriminants.values()} == {105, 170}
    for name, poly in quadratics.items():
        roots = sp.solve(poly.as_expr(), r)
        assert len(roots) == 2 and all(root.is_real is True and root != 0
                                       for root in roots)
        assert all(sp.expand(equations[name].subs({D: 1, X: 1, E: root})) == 0
                   for root in roots)

    # q has S4 splitting field. Its derived subgroup has order 12, hence the
    # unique quadratic subfield is its discriminant field Q(sqrt(17)).
    assert sp.Poly(q, y, domain=sp.QQ).is_irreducible
    group, alternating = galois_group(q, y)
    assert group.order() == 24 and not alternating
    assert group.derived_subgroup().order() == 12
    factors = sp.factorint(sp.discriminant(q, y))
    squarefree = sp.prod(p for p, exponent in factors.items() if exponent % 2)
    assert squarefree == 17
    assert all(squarefree_part(value) != squarefree
               for value in discriminants.values())


def verify_global_flag_cuts(h, rows_by_window):
    expected = {(False,False,False), (False,True,False),
                (False,False,True), (False,True,True),
                (True,False,False), (True,True,False)}
    for rows in rows_by_window.values():
        combos = {(row["branch"] == "T2", bool(case["d2_zero"]),
                   bool(case["sigma_zero"]))
                  for row in rows for case in row["survivor_cases"]}
        assert combos == expected
        for t2, d2_zero, sigma_zero in combos:
            subs = {}
            if t2:
                subs[d1] = 0
            if d2_zero:
                subs[d2] = 0
            if sigma_zero:
                subs[sigma] = 0
            for level in (4, 5, 6):
                cut = sp.Poly(sp.expand(h[level].subs(subs)), *VARS)
                assert cut.terms()
                if len(cut.terms()) == 1:
                    (_, coefficient), = cut.terms()
                    assert coefficient != 0


def verify_incidence(inventories, ordered):
    assert len(inventories["sub2"][0]) == 41
    assert len(inventories["sub1"][0]) == 67
    shape_by_id = {f"C{i:02}": shape for i, shape in enumerate(ordered, 1)}

    def uses(label, name):
        level, support = shape_by_id[name]
        freq, cells, _ = inventories[label]
        keys = [key for key in freq if key[1] == level and key[3] == support
                and key[2] in {"monomial_tie_rise", "identical_vanishing"}]
        if not keys:
            return 0, set()
        return sum(freq[key] for key in keys), set().union(*(cells[key] for key in keys))

    assert uses("sub2", "C04")[0] == 132 and uses("sub1", "C04")[0] == 621
    assert uses("sub2", "C10")[0] == 145 and uses("sub1", "C10")[0] == 910
    assert uses("sub2", "C23")[0] == 145 and uses("sub1", "C23")[0] == 937
    assert (uses("sub2", "C08")[0], len(uses("sub2", "C08")[1])) == (15, 3)
    assert (uses("sub1", "C08")[0], len(uses("sub1", "C08")[1])) == (304, 54)
    assert uses("sub2", "C20") == (0, set())
    assert (uses("sub1", "C20")[0], len(uses("sub1", "C20")[1])) == (17, 8)


def main():
    h = load_h()
    terms = source_terms(h)
    assert tuple(len(terms[level]) for level in (6,5,4)) == (3,5,8)
    print("V1. parsed source and recovered h6/h5/h4 term counts             OK")
    rows = {label: survivors(path) for label, path in WINDOWS.items()}
    inventories = {label: inventory(rows[label], terms) for label in WINDOWS}
    ordered, equations = residue_shapes(inventories, terms)
    print("V2. 41/67 inventories yield 23 source-derived residue shapes     OK")
    verify_constraints(equations)
    print("V3. 21 rational CONSTRAINT shapes; C08/C20 splitting-field KILL OK")
    verify_global_flag_cuts(h, rows)
    print("V4. all occurring global zero-flag cuts and singleton kills      OK")
    verify_incidence(inventories, ordered)
    print("V5. P6/P10/P11 and kill-shape per-window incidence              OK")
    print("\nALL RESIDUE-LEMMA CHECKS PASS")


if __name__ == "__main__":
    main()
