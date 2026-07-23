#!/usr/bin/env python3
"""Exact DEGREE-LAYER (max-plus at infinity) checks for the alternate regime.

Companion to ALT_REGIME_INF.md.  This is the infinity-place analogue of the
standard-regime layer in cascade_engine.py (deg_h_options / descend_options_inf
/ inf_place_profiles), transported to the FLIPPED regime a in [11,15]
(v = 30-3a < 0, w = |v| = 3a-30 > 0) whose polynomial reduction and descending
cascade are derived in ALT_REGIME.md and ALT_REGIME_L2.md:

    F = t^210 G',   G' = sum_{f=0..7} t^((7-f)w) u^f E^(21-3f) h_f(d~),

    (D_t)  T r_6 = h_7                       (top anchor, u r_7 = 0)
           T r_{f-1} = E^(3(7-f)) h_f + u r_f      for f = 6,...,1
           E^21 h_0 + u r_0 = 0               (bottom closing anchor)

with T = t^w, deg u = 4, E = ehat (e = t^a E, deg E = deg e - a).  Taking
degrees (deg t = 1 at infinity) gives the max-plus identity, for every level f,

    w + deg r_{f-1} = max( 3(7-f) deg_E + deg h_f ,  4 + deg r_f ),

the max attained (no drop) unless the two right-hand terms TIE, in which case a
drop below the max is permitted only as a recorded leading-coefficient
cancellation.  The closing anchor must be a tie (both terms equal) for the sum
to vanish; a unique maximum there is a contradiction (nonzero leading term).

METHOD (mirrors alt_regime_verify.py's random-window check 2, lifted to deg):
build a random admissible subcase-(1) window at a=12, form r_6..r_0 as EXACT
rational functions from (D_t) top-down, and verify every derived degree
identity.  The recursion is an exact rational-function identity for ANY window,
so deg(T r_{f-1}) = max(term1, term2) holds with a strict drop only on a tie;
for a generic (non-solution) window no leading coefficients cancel, so the max
is attained at every level and the closing residual is nonzero -- exactly the
degree-layer contradiction that kills a would-be counterexample.
"""
from __future__ import annotations
import random
import re
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parent
NEG_INF = float("-inf")
DEG_U = 4  # deg u = deg(c*q) = deg q, mirrors cascade_engine.DEG_U

d2, d1, d0, dm1, y = sp.symbols("d2 d1 d0 dm1 y")
t = y + 1
qpoly = 2048*y**4 - 512*y**3 + 320*y**2 - 240*y + 195
cconst = sp.Rational(-1, 6630)
u = cconst*qpoly
assert int(sp.degree(sp.Poly(u, y))) == DEG_U

text = (ROOT / "f31_graded.txt").read_text(encoding="utf-8")
pattern = r"h_(\d) \(weight \d+, dm1-power \d+\) = (.+)"
hs = {int(m.group(1)): sp.sympify(m.group(2)) for m in re.finditer(pattern, text)}
assert sorted(hs) == list(range(8))


def deg_rat(expr) -> float:
    """Degree at infinity of a rational function in y (deg num - deg den)."""
    expr = sp.cancel(sp.together(expr))
    num, den = sp.fraction(expr)
    if sp.expand(num) == 0:
        return NEG_INF
    return int(sp.degree(sp.Poly(num, y))) - int(sp.degree(sp.Poly(den, y)))


# --- random admissible window at a=12 (sub1 caps d2<=6,d1<=9,d0<=12,e<=15) ---
random.seed(21012)


def rpoly(degree: int):
    """Monic degree-`degree` polynomial with small random lower coefficients."""
    return sp.expand(y**degree + sum(random.randint(-2, 2)*y**j
                                     for j in range(degree)))


a = 12
v = 30 - 3*a          # = -6  (flipped: v < 0)
w = 3*a - 30          # = 6  = |v| = deg T
D2, D1, D0 = rpoly(6), rpoly(9), rpoly(12)      # caps deg d2<=6, d1<=9, d0<=12
E = rpoly(15 - a)                                # ehat, deg E <= 15-a = 3
if E.subs(y, -1) == 0:
    E += 1
deg_E = int(sp.degree(sp.Poly(E, y)))
assert deg_E == 15 - a and E.subs(y, -1) != 0    # E a unit at t (t = y+1)
e_full = sp.expand(t**a * E)                     # dm1 = e = t^a E, deg e = 15
assert int(sp.degree(sp.Poly(e_full, y))) == 15
sigma = sp.expand(4*D0 - D2**2)
assert int(sp.degree(sp.Poly(sigma, y))) == 12   # generic: no leading cancel

hval = {f: sp.expand(hs[f].subs({d2: D2, d1: D1, d0: D0, dm1: e_full}))
        for f in range(8)}
H = {f: (NEG_INF if sp.expand(hval[f]) == 0
         else int(sp.degree(sp.Poly(hval[f], y)))) for f in range(8)}

# ---------------------------------------------------------------------------
# 1. Flipped reduction orders and the surviving sub1 degree cap deg h_f<=60-6f.
# ---------------------------------------------------------------------------
orders = [30*f + a*(21-3*f) for f in range(8)]
assert v < 0 and w == -v and 21*a + 7*v == 210
assert orders[7] == min(orders) == 210
assert all(orders[f] - 210 == (7-f)*w for f in range(8))
assert all(H[f] <= 60 - 6*f for f in range(8))   # ALT_REGIME.md caps survive
print("1. flipped orders (unique min f=7 ->210) and deg h_f<=60-6f       OK")

# ---------------------------------------------------------------------------
# 2. Reduction F = t^210 G' on the random window (mirrors alt_regime_verify 2).
# ---------------------------------------------------------------------------
phi_tilde = sp.expand(cconst*t**30*qpoly)
assert sp.cancel(phi_tilde/(t**30*u)) == 1
for f in range(8):
    assert 30*f + a*(21-3*f) == 210 + (7-f)*w
    left = phi_tilde**f * e_full**(21-3*f) * hval[f]
    right = t**(210+(7-f)*w) * u**f * E**(21-3*f) * hval[f]
    for sample in (-2, 0, 1, 3):
        assert sp.cancel((left-right).subs(y, sample)) == 0
print("2. random a=12 window: F = t^210 G' term-by-term (deg t=1)        OK")

# ---------------------------------------------------------------------------
# 3. Build descending cofactors r_6..r_0 as exact rational functions (D_t).
#    Verify the exact recursion identity used by the degree layer.
# ---------------------------------------------------------------------------
T = sp.expand(t**w)
r: dict[int, sp.Expr] = {}
r[6] = sp.cancel(hval[7] / T)                              # T r_6 = h_7
for f in range(6, 0, -1):
    rhs = sp.expand(E**(3*(7-f))*hval[f]) + u*r[f]         # E^(3(7-f)) h_f + u r_f
    r[f-1] = sp.cancel(rhs / T)                            # = T r_{f-1}
assert sp.cancel(T*r[6] - hval[7]) == 0
for f in range(6, 0, -1):
    assert sp.cancel(T*r[f-1] - (E**(3*(7-f))*hval[f] + u*r[f])) == 0
print("3. exact descending cascade (D_t) rebuilt as rational cofactors   OK")

# ---------------------------------------------------------------------------
# 4. Top anchor  T r_6 = h_7:  single term forces  w + deg r_6 = H_7.
# ---------------------------------------------------------------------------
R = {f: deg_rat(r[f]) for f in range(7)}
assert deg_rat(T*r[6]) == H[7]           # unique achiever, no drop
assert w + R[6] == H[7]                  # derived: deg r_6 = H_7 - w
print(f"4. top anchor: w+deg r_6 = {int(w+R[6])} = H_7 (forced)              "
      "     OK")

# ---------------------------------------------------------------------------
# 5. Intermediate max-plus levels f = 6..1.
#    w + deg r_{f-1} = max(3(7-f) deg_E + H_f, 4 + deg r_f); drop only on tie.
# ---------------------------------------------------------------------------
derivedR = {6: H[7] - w}
dropcount = 0
for f in range(6, 0, -1):
    term1 = 3*(7-f)*deg_E + H[f]                 # deg(E^(3(7-f)) h_f)
    term2 = DEG_U + R[f]                         # deg(u r_f)
    assert term1 == deg_rat(E**(3*(7-f))*hval[f])
    assert term2 == deg_rat(u*r[f])
    lhs = deg_rat(T*r[f-1])                       # deg(T r_{f-1}) = w + deg r_{f-1}
    assert lhs == w + R[f-1]
    mx = max(term1, term2)
    assert lhs <= mx                              # sum degree never exceeds max
    if lhs < mx:
        assert term1 == term2                     # a drop REQUIRES a tie
        dropcount += 1
    if term1 != term2:
        assert lhs == mx                          # unique max -> forced, no drop
    derivedR[f-1] = mx - w                         # max-plus recursion
    if lhs == mx:
        assert derivedR[f-1] == R[f-1]            # derived formula reproduces deg
print(f"5. intermediate levels f=6..1: max-plus identities ({dropcount} drops) "
      "  OK")

# ---------------------------------------------------------------------------
# 6. Bottom closing anchor  E^21 h_0 + u r_0 = 0.
#    Tie forced for a solution; a generic window has a unique max -> residual
#    nonzero -> the degree layer alone contradicts a would-be counterexample.
# ---------------------------------------------------------------------------
term1_0 = 21*deg_E + H[0]
term2_0 = DEG_U + R[0]
assert term1_0 == deg_rat(E**21*hval[0])
assert term2_0 == deg_rat(u*r[0])
residual = sp.expand(E**21*hval[0] + u*r[0])
dres = deg_rat(residual)
mx0 = max(term1_0, term2_0)
assert dres <= mx0
if dres < mx0:
    assert term1_0 == term2_0
# generic non-solution window: strictly dominant term -> nonzero residual.
assert term1_0 != term2_0 and dres == mx0 and sp.expand(residual) != 0
print(f"6. closing anchor: unique max deg={int(mx0)} -> residual != 0 "
      "(contradiction) OK")

print("\nALL ALTERNATE-REGIME INFINITY (DEGREE-LAYER) CHECKS PASS")
