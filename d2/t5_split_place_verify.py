#!/usr/bin/env python3
"""Exact checks for the split-place repair of the f31 T5 campaign.

This script is deliberately small.  It does not try to prove UFD or
Mason--Stothers; it checks every finite piece of arithmetic used in:

  1. the field-stable split-place repair of the sigma-locus argument; and
  2. the new (a_t=7, geometrically q-coprime) f31/subcase-2 proof.

Run with Python 3 and SymPy:

    python3 t5_split_place_verify.py
"""

from __future__ import annotations

from itertools import product
from typing import Iterable

import sympy as sp

def _require(_cond, _msg):
    """Proof-critical check: fails loudly and exits nonzero, unaffected by python -O."""
    if not _cond:
        import sys as _sys
        print("FAIL: " + str(_msg))
        _sys.exit(1)



# ---------------------------------------------------------------------------
# Part A. Split-place valuation ledger for A^4 B = const * e^17.


def local_options(phi_order: int) -> list[tuple[int, int, int]]:
    """Return (m, alpha, beta) at a root of Phi of order ``phi_order``.

    Here m=v(e), alpha=v(A), beta=v(B), and
        4*alpha + beta = 17*m.
    Since a nonzero linear combination of A and B equals Phi:
      * alpha != beta  => min(alpha,beta)=phi_order;
      * alpha == beta  => alpha <= phi_order (cancellation may raise order).
    Global degree bounds give alpha,beta <= 34 and m <= 10.
    """

    out: list[tuple[int, int, int]] = []
    for m in range(11):
        for alpha in range(35):
            beta = 17 * m - 4 * alpha
            if not 0 <= beta <= 34:
                continue
            if alpha != beta:
                if min(alpha, beta) != phi_order:
                    continue
            elif alpha > phi_order:
                continue
            out.append((m, alpha, beta))
    return out


T_OPTIONS = local_options(30)
Q_OPTIONS = local_options(1)

_require(T_OPTIONS == [(0, 0, 0), (5, 17, 17), (9, 30, 33)], "T_OPTIONS == [(0, 0, 0), (5, 17, 17), (9, 30, 33)]")
_require(Q_OPTIONS == [
    (0, 0, 0),
    (1, 1, 13),
    (1, 4, 1),
    (2, 1, 30),
    (5, 21, 1),
], "Q_OPTIONS == [ (0, 0, 0), (1, 1, 13), (1, 4, 1), (2, 1, 30), (5, 21, 1), ]")


def enumerate_global_patterns() -> list[dict[str, object]]:
    """Enumerate all degree-compatible local patterns after base change.

    Outside the five roots of Phi (t plus four simple roots of q), a root of e
    is either:
      * type A: e-multiplicity 4k and A-multiplicity 17k; or
      * type B: e-multiplicity m and B-multiplicity 17m.

    KA is the sum of the type-A k's; MB is the total type-B e-mass.
    """

    patterns: list[dict[str, object]] = []
    phis = (30, 1, 1, 1, 1)
    for t_opt in T_OPTIONS:
        for q_opts in product(Q_OPTIONS, repeat=4):
            local = (t_opt,) + q_opts
            local_e = sum(row[0] for row in local)
            local_a = sum(row[1] for row in local)
            local_b = sum(row[2] for row in local)
            for ka in range(3):  # 4*ka <= deg e = 10
                for mb in range(11):
                    if local_e + 4 * ka + mb != 10:
                        continue
                    if local_a + 17 * ka != 34:
                        continue
                    if local_b + 17 * mb != 34:
                        continue

                    gcd_degree = sum(min(row[1], row[2]) for row in local)
                    active_special_roots = sum(
                        not (row[1] == row[2] == phi)
                        for row, phi in zip(local, phis)
                    )
                    radical_bound = ka + mb + active_special_roots
                    reduced_degree = 34 - gcd_degree
                    mason_gap = reduced_degree - (radical_bound - 1)
                    patterns.append(
                        {
                            "t": t_opt,
                            "q": q_opts,
                            "KA": ka,
                            "MB": mb,
                            "gcd_degree": gcd_degree,
                            "radical_bound": radical_bound,
                            "reduced_degree": reduced_degree,
                            "mason_gap": mason_gap,
                        }
                    )
    return patterns


PATTERNS = enumerate_global_patterns()
_require(len(PATTERNS) == 6, "len(PATTERNS) == 6")

mason_killed = [p for p in PATTERNS if int(p["mason_gap"]) > 0]
exceptional = [p for p in PATTERNS if int(p["mason_gap"]) <= 0]
_require(len(mason_killed) == 2, "len(mason_killed) == 2")
_require(len(exceptional) == 4, "len(exceptional) == 4")  # same pattern, one choice of the q-root

# The exceptional pattern is e=t^9*p, with (v_p(A),v_p(B))=(4,1).
# The second linear relation has the form const*A + const*B = const*d*e^3.
# Its left side has order 1 at p; the right side has order >= 3.  Contradiction.
for row in exceptional:
    q_rows = row["q"]
    _require(isinstance(q_rows, tuple), "isinstance(q_rows, tuple)")
    nonzero = [triple for triple in q_rows if triple != (0, 0, 0)]
    _require(nonzero == [(1, 4, 1)], "nonzero == [(1, 4, 1)]")
    left_order = min(4, 1)
    right_order_lower_bound = 3 * 1
    _require(left_order == 1 < right_order_lower_bound, "left_order == 1 < right_order_lower_bound")


# Correct geometric stratum count: sorted multiplicity vectors at the four
# simple roots of q, with a+sum(b_i)<=10.
def sorted_q_vectors(total_cap: int) -> list[tuple[int, int, int, int]]:
    return [
        b
        for b in product(range(total_cap + 1), repeat=4)
        if sum(b) <= total_cap and tuple(sorted(b, reverse=True)) == b
    ]


_require(sum(len(sorted_q_vectors(10 - a)) for a in range(11)) == 327, "sum(len(sorted_q_vectors(10 - a)) for a in range(11)) == 327")


# ---------------------------------------------------------------------------
# Part B. Algebra for the a_t=7, gcd(E,q)=1 branch of f31/subcase 2.

E, G, g5, d2, sigma, t, q, c = sp.symbols(
    "E G g5 d2 sigma t q c", nonzero=True
)

# On d1=0:
#   h6=-3072*sigma^2,
#   h5=-9216*d2*sigma^2+2048*e^2,
# with e=t^7 E.
h5_t2 = -9216 * d2 * sigma**2 + 2048 * t**14 * E**2

# Terminal level: E^3 G = 3072*c^6*sigma^2 after g6=q^6 G.
raw_level5 = E**3 * g5 + c**5 * q**5 * h5_t2 - t**9 * q**6 * G
raw_level5 = sp.expand(raw_level5.subs(sigma**2, E**3 * G / (3072 * c**6)))

# Exact rearrangement used in the divisibility argument.
rearranged = E**3 * (g5 - (3 / c) * q**5 * d2 * G) - t**9 * q**5 * (
    q * G - 2048 * c**5 * t**5 * E**2
)
_require(sp.cancel(raw_level5 - rearranged) == 0, "sp.cancel(raw_level5 - rearranged) == 0")
_require(sp.simplify(-3 / sp.Rational(-1, 6630)) == 19890, "sp.simplify(-3 / sp.Rational(-1, 6630)) == 19890")

# T2 degree facts.  If N=qG-2048*c^5*t^5 E^2=E^3 S, then
# deg(t^9 q^5 S)<=31 forces deg S<=2.  The subsequent quotient L has deg<=1.
_require(9 + 5 * 4 == 29, "9 + 5 * 4 == 29")
_require(31 - 29 == 2, "31 - 29 == 2")
_require(max(4 + 7, 5 + 2 * 3) == 11, "max(4 + 7, 5 + 2 * 3) == 11")
_require(11 - 3 * 3 == 2, "11 - 3 * 3 == 2")
_require(5 - 4 == 1, "5 - 4 == 1")

# T1 UFD parameterization:
# E = gamma*s*u^2, H=eta*s*v^2, d1=delta*s^2*u^3*v.
# Enumerate the only local parity escape when deg(E),deg(H)<=3.
local_escapes: list[tuple[int, int, int]] = []
for e_order in (1, 3):
    for h_order in (1, 3):
        d1_order = (3 * e_order + h_order) // 2
        _require(2 * d1_order == 3 * e_order + h_order, "2 * d1_order == 3 * e_order + h_order")

        # Level 6 orders away from t and q:
        # E^3*G6, sigma^2, d1^2*d2, d1*e, RHS H.
        # We ask whether some even sigma^2 order and some G6 order can make
        # the minimum occur at least twice.
        possible = False
        for sigma_order in range(9):
            for g6_order in range(8):
                orders = (
                    3 * e_order + g6_order,
                    2 * sigma_order,
                    2 * d1_order,
                    d1_order + e_order,
                    h_order,
                )
                m = min(orders)
                if sum(order == m for order in orders) >= 2:
                    possible = True
                    break
            if possible:
                break
        if possible:
            local_escapes.append((e_order, h_order, d1_order))

_require(local_escapes == [(1, 3, 3)], "local_escapes == [(1, 3, 3)]")

# In that sole escape, level 6 forces v_s(sigma)>=2 and v_s(G6)=0.
# Then level 5 has orders: E^3 g5 >=3; h5 has unique e^2 order 2;
# RHS g6 has order 0.  Impossible.
exceptional_level5_orders = {
    "E3g5_lower": 3,
    "d2_sigma2_lower": 4,
    "sigma_d1sq_lower": 8,
    "d2sq_d1sq_lower": 6,
    "d1_d2_e_lower": 4,
    "e2": 2,
    "rhs_g6": 0,
}
_require(min(
    exceptional_level5_orders["E3g5_lower"],
    exceptional_level5_orders["d2_sigma2_lower"],
    exceptional_level5_orders["sigma_d1sq_lower"],
    exceptional_level5_orders["d2sq_d1sq_lower"],
    exceptional_level5_orders["d1_d2_e_lower"],
    exceptional_level5_orders["e2"],
) == 2, "min( exceptional_level5_orders[\"E3g5_lower\"], exceptional_level5_orders[\"d2_sigma2_lower\"], exceptional_level5_orders[\"sigma_d1sq_lower\"], exceptional_level5_orders[\"d2sq_d1sq_lower\"], exceptional_level5_orders[\"d1_d2_e_lower\"], exceptional_level5_orders[\"e2\"], ) == 2")
_require(exceptional_level5_orders["rhs_g6"] == 0, "exceptional_level5_orders[\"rhs_g6\"] == 0")

# With the squarefree factor s gone, deg u,deg v are 0 or 1.  The following
# table checks the degree-domination argument for every case and every possible
# degree of sigma.  The top term must be unique among T_0,...,T_7.

def top_term_indices(deg_u: int, deg_v: int, deg_sigma: int) -> tuple[int, ...]:
    deg_e = 7 + 2 * deg_u
    deg_d1 = 3 * deg_u + deg_v

    # Coarse but sharp enough for f<=5, with a refinement for (u,v)=(1,0)
    # and deg_sigma<=6 where the generic h5 bound would tie artificially.
    degrees: list[int] = []
    for f in range(6):
        h_bound = 40 - 4 * f
        if deg_u == 1 and deg_v == 0 and deg_sigma <= 6 and f == 5:
            # h5=max(d2*sigma^2, sigma*d1^2, d2^2*d1^2,
            #        d1*d2*e, e^2) = 18 here.
            h_bound = 18
        degrees.append(34 * f + (21 - 3 * f) * deg_e + h_bound)

    h6_bound = max(
        2 * deg_sigma,
        2 * deg_d1 + 4,
        deg_d1 + deg_e,
    )
    degrees.append(6 * 34 + 3 * deg_e + h6_bound)
    degrees.append(7 * 34 + 2 * deg_d1)

    top = max(degrees)
    return tuple(i for i, degree in enumerate(degrees) if degree == top)


for deg_u, deg_v, deg_sigma in product((0, 1), (0, 1), range(9)):
    tops = top_term_indices(deg_u, deg_v, deg_sigma)
    _require(len(tops) == 1, (deg_u, deg_v, deg_sigma, tops))


print("split-place sigma ledger: PASS")
print("  local t options:", T_OPTIONS)
print("  local simple-q options:", Q_OPTIONS)
print("  global patterns: 6 = 2 Mason-killed + 4 symmetric exceptional")
print("  exceptional pattern killed by v_p(d*e^3)>=3 > 1")
print("a_t=7 geometrically q-coprime f31 branch: PASS")
print("  T2 level-5 rearrangement and degree squeeze checked")
print("  T1 UFD local parity escape and infinity table checked")
