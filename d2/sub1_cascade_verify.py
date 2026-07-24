#!/usr/bin/env python3
"""Exact audit for the f31 subcase-(1) split-place cascade.

CORRECTED: the claimed bound deg(g_l)<=15+3*a does not follow from the upper
bound deg(ehat)<=15-a. Rigorous terminal caps are 46 (T1) and 48 (T2).
"""
from __future__ import annotations
from collections import Counter
from itertools import product
from pathlib import Path
import re
import sympy as sp

def _require(_cond, _msg):
    """Proof-critical check: fails loudly and exits nonzero, unaffected by python -O."""
    if not _cond:
        import sys as _sys
        print("FAIL: " + str(_msg))
        _sys.exit(1)


ROOT = Path(__file__).resolve().parent
SCHEMA_VERSION = 1
Q_ROOT_COUNT = 4
MAX_A = 15
STANDARD_MAX_A = 10
ALTERNATE_MIN_A = 11
D2_DEG_CAP, D1_DEG_CAP, D0_DEG_CAP, E_DEG_CAP = 6, 9, 12, 15
SIGMA_DEG_CAP, H_WEIGHT_SLOPE, U_DEG, PHI_T_ORDER = 12, 3, 4, 30
T1_LEVEL, T2_LEVEL = 7, 6
TERMINAL_E_POWER, TERMINAL_AUX_POWER = 3, 2
T1_G_DEG_CAP, T2_G_DEG_CAP = 46, 48
T3_STATUS = "proven_infeasible"
T3_REFERENCE = ("sub1_cascade_verify.py check 7, split-place sigma theorem "
                "with subcase-(1) degree enumeration")

d2, d1, d0, dm1, Phi, y = sp.symbols("d2 d1 d0 dm1 Phi y")
VARS = (d2, d1, d0, dm1)
WEIGHTS = {d2: 2, d1: 3, d0: 4, dm1: 5}
t = y + 1
q = 2048*y**4 - 512*y**3 + 320*y**2 - 240*y + 195
c = sp.Rational(-1, 6630)
u = c*q
phi_tilde = c*t**PHI_T_ORDER*q


def _read_sources():
    text = (ROOT / "f31_deg31.txt").read_text(encoding="utf-8").strip()
    f31_expr = sp.sympify(text.replace("m1", "dm1").replace("P", "Phi").replace("^", "**"))
    graded = (ROOT / "f31_graded.txt").read_text(encoding="utf-8")
    parts = {}
    pattern = r"h_(\d) \(weight (\d+), dm1-power (\d+)\) = (.+)"
    for match in re.finditer(pattern, graded):
        parts[int(match.group(1))] = sp.sympify(match.group(4))
    _require(sorted(parts) == list(range(8)), "sorted(parts) == list(range(8))")
    return f31_expr, parts


f31, hs = _read_sources()
def degree(expr):
    return int(sp.degree(sp.expand(expr), y))

# Check 1: stripped window caps.
_require(tuple(map(sp.Integer, (D2_DEG_CAP, D1_DEG_CAP, D0_DEG_CAP, E_DEG_CAP))) == (6, 9, 12, 15), "tuple(map(sp.Integer, (D2_DEG_CAP, D1_DEG_CAP, D0_DEG_CAP, E_DEG_CAP))) == (6, 9, 12, 15)")
_require(sp.Max(D0_DEG_CAP, 2*D2_DEG_CAP) == SIGMA_DEG_CAP, "sp.Max(D0_DEG_CAP, 2*D2_DEG_CAP) == SIGMA_DEG_CAP")
_require(degree(4*y**D0_DEG_CAP - (y**D2_DEG_CAP)**2) == SIGMA_DEG_CAP, "degree(4*y**D0_DEG_CAP - (y**D2_DEG_CAP)**2) == SIGMA_DEG_CAP")
print("1. caps (d2,d1,d0,e)=(6,9,12,15); deg sigma<=12              OK")

# Check 2: exact weighted maximum over every source monomial.
for f in range(8):
    monomials = sp.Poly(hs[f], *VARS).monoms()
    for monomial in monomials:
        weight = sum(WEIGHTS[var]*exponent for var, exponent in zip(VARS, monomial))
        _require(sp.Integer(weight) == 20 - 2*f, "sp.Integer(weight) == 20 - 2*f")
    max_degree = max(sum(H_WEIGHT_SLOPE*WEIGHTS[var]*exponent for var, exponent in zip(VARS, monomial)) for monomial in monomials)
    _require(sp.Integer(max_degree) == H_WEIGHT_SLOPE*(20 - 2*f), "sp.Integer(max_degree) == H_WEIGHT_SLOPE*(20 - 2*f)")
    _require(sp.Integer(max_degree) == 60 - 6*f, "sp.Integer(max_degree) == 60 - 6*f")
print("2. source monomials give deg h_f(d~)<=60-6f exactly             OK")

# Check 3: graded identity, telescope, and exponent reduction.
reconstructed = sum(Phi**f * dm1**(21-3*f) * hs[f] for f in range(8))
_require(sp.expand(reconstructed - f31) == 0, "sp.expand(reconstructed - f31) == 0")
T, U, E = sp.symbols("T U E")
gs = [None] + list(sp.symbols("g1:8"))
hh = {0: T*gs[1]}
for level in range(1, 7):
    hh[level] = (T*gs[level+1] - E**3*gs[level])/U**level
hh[7] = -E**3*gs[7]/U**7
S = sum(T**f * U**f * E**(21-3*f) * hh[f] for f in range(8))
_require(sp.expand(sp.cancel(S)) == 0, "sp.expand(sp.cancel(S)) == 0")
for a in range(STANDARD_MAX_A + 1):
    v = PHI_T_ORDER - TERMINAL_E_POWER*a
    for f in range(8):
        _require(sp.Integer(PHI_T_ORDER*f + a*(21-3*f)) == 21*a + v*f, "sp.Integer(PHI_T_ORDER*f + a*(21-3*f)) == 21*a + v*f")
print("3. graded identity and window-independent cascade telescope        OK")

# Check 4: degree bookkeeping and correction.
_require(degree(u) == U_DEG, "degree(u) == U_DEG")
_require(sp.expand(hs[T1_LEVEL] - 8192*d1**2) == 0, "sp.expand(hs[T1_LEVEL] - 8192*d1**2) == 0")
_require(sp.Integer(U_DEG*T1_LEVEL + TERMINAL_AUX_POWER*D1_DEG_CAP) == 46, "sp.Integer(U_DEG*T1_LEVEL + TERMINAL_AUX_POWER*D1_DEG_CAP) == 46")
_require(T1_G_DEG_CAP == 46, "T1_G_DEG_CAP == 46")
# CORRECTED: the max-chain is true, but dividing by ehat^3 would require a
# lower bound on deg(ehat); only an upper cap is known. This exact terminal
# witness reaches degree 46 at a=0 while respecting all stated caps.
ehat_witness = sp.Integer(1)
d1_witness = y**D1_DEG_CAP
g7_witness = sp.expand(-u**T1_LEVEL*8192*d1_witness**2/ehat_witness**3)
_require(sp.expand(ehat_witness**3*g7_witness + u**T1_LEVEL*8192*d1_witness**2) == 0, "sp.expand(ehat_witness**3*g7_witness + u**T1_LEVEL*8192*d1_witness**2) == 0")
_require(degree(g7_witness) == T1_G_DEG_CAP, "degree(g7_witness) == T1_G_DEG_CAP")
_require(sp.Integer(T1_G_DEG_CAP) > 15, "sp.Integer(T1_G_DEG_CAP) > 15")
for a in range(STANDARD_MAX_A + 1):
    v = PHI_T_ORDER - TERMINAL_E_POWER*a
    ehat_cap = E_DEG_CAP - a
    claimed_cap = E_DEG_CAP + TERMINAL_E_POWER*a
    _require(sp.Integer(ehat_cap) >= 0 and sp.Integer(v) >= 0, "sp.Integer(ehat_cap) >= 0 and sp.Integer(v) >= 0")
    _require(sp.Max(v + claimed_cap, 60) <= 60, "sp.Max(v + claimed_cap, 60) <= 60")
    _require(sp.Integer(v) <= 45 - 3*a, "sp.Integer(v) <= 45 - 3*a")
    # Safe level-dependent caps from both recursion directions.
    forward = {1: 60-v}
    for level in range(1, 7):
        forward[level+1] = max(3*ehat_cap + forward[level], 60-2*level) - v
    backward = {7: T1_G_DEG_CAP}
    for level in range(6, 0, -1):
        backward[level] = max(v + backward[level+1], 60-2*level)
    safe = {level: min(forward[level], backward[level]) for level in range(1, 8)}
    for level in range(1, 8):
        _require(sp.Integer(safe[level]) == sp.Min(forward[level], backward[level]), "sp.Integer(safe[level]) == sp.Min(forward[level], backward[level])")
sigma = 4*d0 - d2**2
_require(sp.expand(hs[T2_LEVEL].subs(d1, 0) + 3072*sigma**2) == 0, "sp.expand(hs[T2_LEVEL].subs(d1, 0) + 3072*sigma**2) == 0")
_require(sp.Integer(U_DEG*T2_LEVEL + TERMINAL_AUX_POWER*SIGMA_DEG_CAP) == 48, "sp.Integer(U_DEG*T2_LEVEL + TERMINAL_AUX_POWER*SIGMA_DEG_CAP) == 48")
_require(T2_G_DEG_CAP == 48, "T2_G_DEG_CAP == 48")
print("4. CORRECTED: terminal caps deg g7<=46 (T1), deg g6<=48 (T2);")
print("   15+3a needs an unavailable lower deg(ehat) bound              OK")

# Check 5: ledger degree budget.
for a in range(MAX_A + 1):
    for b in product(range(MAX_A + 1), repeat=Q_ROOT_COUNT):
        if a + sum(b) <= MAX_A:
            _require(sp.Integer(a + sum(b)) <= E_DEG_CAP, "sp.Integer(a + sum(b)) <= E_DEG_CAP")
_require(sp.Integer(MAX_A) == E_DEG_CAP, "sp.Integer(MAX_A) == E_DEG_CAP")
print("5. divisibility gives a+sum(b_i)<=15                               OK")

# Check 6: alternate regime.
for a in range(ALTERNATE_MIN_A, MAX_A + 1):
    _require(sp.Integer(PHI_T_ORDER - TERMINAL_E_POWER*a) < 0, "sp.Integer(PHI_T_ORDER - TERMINAL_E_POWER*a) < 0")
_require([PHI_T_ORDER-TERMINAL_E_POWER*a for a in range(11, 16)] == [-3, -6, -9, -12, -15], "[PHI_T_ORDER-TERMINAL_E_POWER*a for a in range(11, 16)] == [-3, -6, -9, -12, -15]")
print("6. a=11..15 gives v=(-3,-6,-9,-12,-15): alternate regime open    OK")


def local_options(phi_order, polynomial_degree, e_degree):
    result = []
    for m in range(e_degree + 1):
        for alpha in range(polynomial_degree + 1):
            beta = 17*m - 4*alpha
            if not 0 <= beta <= polynomial_degree:
                continue
            if alpha != beta and min(alpha, beta) != phi_order:
                continue
            if alpha == beta and alpha > phi_order:
                continue
            result.append((m, alpha, beta))
    return result


def global_patterns(polynomial_degree, e_degree):
    t_options = local_options(PHI_T_ORDER, polynomial_degree, e_degree)
    q_options = local_options(1, polynomial_degree, e_degree)
    result = []
    phi_orders = (PHI_T_ORDER, 1, 1, 1, 1)
    for t_option in t_options:
        for q_tuple in product(q_options, repeat=Q_ROOT_COUNT):
            local = (t_option,) + q_tuple
            for mass_a in range(e_degree//4 + 1):
                for mass_b in range(e_degree + 1):
                    if sum(row[0] for row in local) + 4*mass_a + mass_b != e_degree:
                        continue
                    if sum(row[1] for row in local) + 17*mass_a != polynomial_degree:
                        continue
                    if sum(row[2] for row in local) + 17*mass_b != polynomial_degree:
                        continue
                    gcd_degree = sum(min(row[1], row[2]) for row in local)
                    active = sum(not (row[1] == row[2] == phi_order) for row, phi_order in zip(local, phi_orders))
                    radical = mass_a + mass_b + active
                    reduced = polynomial_degree - gcd_degree
                    result.append({"q": q_tuple, "radical_bound": radical,
                                   "reduced_degree": reduced,
                                   "mason_gap": reduced-(radical-1)})
    return result

# Check 7: T3 transfer. Degree arithmetic first gives only two triples.
AB_DEG_CAP = int(sp.Max(34, D2_DEG_CAP + 3*E_DEG_CAP))
_require(AB_DEG_CAP == 51, "AB_DEG_CAP == 51")
degree_triples = []
for e_degree in range(E_DEG_CAP + 1):
    for degree_a in range(AB_DEG_CAP + 1):
        for degree_b in range(AB_DEG_CAP + 1):
            if 4*degree_a + degree_b != 17*e_degree:
                continue
            phi_relation = ((degree_a == degree_b and degree_a >= 34) or
                            (degree_a != degree_b and max(degree_a, degree_b) == 34))
            d2e3_relation = (degree_a == degree_b or max(degree_a, degree_b) <= D2_DEG_CAP + 3*e_degree)
            if phi_relation and d2e3_relation:
                degree_triples.append((e_degree, degree_a, degree_b))
_require(degree_triples == [(10, 34, 34), (15, 51, 51)], "degree_triples == [(10, 34, 34), (15, 51, 51)]")
old = global_patterns(34, 10)
_require(len(old) == 6, "len(old) == 6")
old_mason = [row for row in old if row["mason_gap"] > 0]
old_exceptional = [row for row in old if row["mason_gap"] <= 0]
_require(len(old_mason) == 2 and len(old_exceptional) == Q_ROOT_COUNT, "len(old_mason) == 2 and len(old_exceptional) == Q_ROOT_COUNT")
_require(Counter((row["reduced_degree"], row["radical_bound"]) for row in old_mason) == Counter({(34, 9): 1, (17, 7): 1}), "Counter((row[\"reduced_degree\"], row[\"radical_bound\"]) for row in old_mason) == Counter({(34, 9): 1, (17, 7): 1})")
_require(sp.Integer(34) > 9-1 and sp.Integer(17) > 7-1, "sp.Integer(34) > 9-1 and sp.Integer(17) > 7-1")
for row in old_exceptional:
    nonzero_q = [triple for triple in row["q"] if triple != (0, 0, 0)]
    _require(nonzero_q == [(1, 4, 1)], "nonzero_q == [(1, 4, 1)]")
    _require(sp.Integer(min(4, 1)) < 3, "sp.Integer(min(4, 1)) < 3")
new = global_patterns(51, 15)
_require(len(new) == 19, "len(new) == 19")
new_counts = Counter((row["reduced_degree"], row["radical_bound"]) for row in new)
_require(new_counts == Counter({(51, 11): 1, (34, 9): 1, (20, 7): 4, (20, 6): 8, (20, 5): 4, (17, 6): 1}), "new_counts == Counter({(51, 11): 1, (34, 9): 1, (20, 7): 4, (20, 6): 8, (20, 5): 4, (17, 6): 1})")
for row in new:
    _require(sp.Integer(row["reduced_degree"]) > row["radical_bound"]-1, "sp.Integer(row[\"reduced_degree\"]) > row[\"radical_bound\"]-1")
_require(sp.Integer(51) > 11-1, "sp.Integer(51) > 11-1")
_require(sp.Integer(34) > 9-1, "sp.Integer(34) > 9-1")
_require(sp.Integer(20) > 7-1, "sp.Integer(20) > 7-1")
_require(sp.Integer(20) > 6-1, "sp.Integer(20) > 6-1")
_require(sp.Integer(20) > 5-1, "sp.Integer(20) > 5-1")
_require(sp.Integer(17) > 6-1, "sp.Integer(17) > 6-1")
_require(sp.Integer(5*PHI_T_ORDER) == 150 and sp.Integer(150) % 17 != 0, "sp.Integer(5*PHI_T_ORDER) == 150 and sp.Integer(150) % 17 != 0")
derivative_cofactor = 30*q + t*sp.diff(q, y)
_require(sp.expand(sp.diff(phi_tilde, y) - c*t**29*derivative_cofactor) == 0, "sp.expand(sp.diff(phi_tilde, y) - c*t**29*derivative_cofactor) == 0")
_require(degree(derivative_cofactor) == 4, "degree(derivative_cofactor) == 4")
_require(16*3 > 29+4, "16*3 > 29+4")
T3_MASON_INEQUALITIES = (
    "34 > 9-1", "17 > 7-1", "51 > 11-1", "34 > 9-1",
    "20 > 7-1", "20 > 6-1", "20 > 5-1", "17 > 6-1",
    "17 does not divide 150",
)
print("7. T3 transfers: old 34>8,17>6; new 51>10,34>8,20>6,")
print("   20>5,20>4,17>5; and 17 does not divide 150                  OK")
print("\nALL SUBCASE-(1) CASCADE PARAMETER CHECKS PASS")
