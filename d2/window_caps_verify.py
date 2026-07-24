"""window_caps_verify.py — full recitation of the k=6,7,8 window caps.

CLAIMS (the caps, per window variable d_{4-k} = dm(k-4), k = 6,7,8):

    var  | k |  ord >=  | deg <= (sub1) | deg <= (sub2)
    dm2  | 6 |   72     |     90        |     84
    dm3  | 7 |   84     |    105        |     98
    dm4  | 8 |   96     |    120        |    112

i.e. ord >= 12k, deg <= 15k (sub1) / 14k (sub2) — the three rows that
FULL_SYSTEM_BRIDGE.md flags as "[judgment] extension of T3_WINDOW_AUDIT §4"
and that sit under EVERY bridge kill (they bound the dm2,dm3,dm4 ansaetze).

This file recites the whole derivation mechanically, from the audited
premises, with NO hand-copied polygon data:

  W0  Newton-polygon corners loaded from paper_src/upstream_facts.json
      (sha-pinned transcription of Prop 4.3, audited in T3_WINDOW_AUDIT §1);
      the three direction maxima the induction uses are COMPUTED from them.
  W1  The three valuation inductions of T3 §3 close for ALL k >= 1 — the
      base, product-bound, slice-bound and step are verified as SYMBOLIC
      identities in k (not just for k = 2..5).
  W2  The D-transform arithmetic D_k := C_k*C4^(7-2k) converts the induction
      bounds into deg(d_{4-k}) <= 15k/14k, ord >= 12k — symbolic in k, plus
      the literal k = 6,7,8 rows.
  W3  The d3-killing shift x -> x - s (s = c_3/(4*C4), the unique choice
      killing c~_3) in D-coordinates is the POLYNOMIAL identity
          D~_j = sum_{k=j..4} binom(k,k-j) * D_k * (-D_3/4)^(k-j),
      verified (a) as an exact series-recomposition identity, (b) to kill
      D~_3, (c) to cancel all C4 exponents (shifted window vars stay
      polynomial), (d) to preserve all three caps term-by-term.
  W4  End-to-end generic corroboration: random exact-QQ polynomials P_i
      supported on the corner hull (computed, not transcribed), the
      division-free D-recursion of verify_derivation.py §C run down to
      D_{-4}, and every cap checked on the unshifted AND shifted variables
      for k = 2..8, both regimes.
  W5  Consumer cross-checks: full_system_bridge.WEIGHT / STRIP_DEGCAP and
      jetlift.CONFIGS window sizes agree with the caps; spare-unknown
      totals 66 (sub1) / 45 (sub2); Phi (weight 17) attains the sub2 caps
      (tightness witness), with f1 re-derived from its forcing ODE here.

INHERITED premises (same tier as the k=2..5 window, nothing NEW):
  [P1] Prop 4.3 corner sets (transcription audited verbatim, T3 §1).
  [P2] C4 = y^7(y+1) normalisation (forced by corners (8,14),(8,16)).
  [P3] Existence of the Laurent square root C of P (paper template).

Run:  python3 window_caps_verify.py [--quiet]     exit 0 iff all checks pass.
"""
import json
import os
import random
import sys
from fractions import Fraction

import sympy as sp
from sympy import Poly, Rational, binomial, expand, symbols

QUIET = "--quiet" in sys.argv
HERE = os.path.dirname(os.path.abspath(__file__))

y = symbols("y")
C4 = y**7 * (y + 1)
ok = [0]


def check(name, cond):
    if not cond:
        print(f"  [FAIL] {name}")
        raise SystemExit(1)
    ok[0] += 1
    if not QUIET:
        print(f"  [OK] {name}")


def vmax(a, b, corners):
    """max of the linear form a*i+b*j over the corner set (= over the whole
    polygon, since a linear form maxes at a vertex)."""
    return max(a * i + b * j for (i, j) in corners)


# ---------------------------------------------------------------- W0. corners
if not QUIET:
    print("W0. polygon inputs computed from the pinned Prop-4.3 corner data")
UF = json.load(open(os.path.join(HERE, "paper_src", "upstream_facts.json")))
NP = UF["facts"]["newton_polygons"]
CORNERS = {reg: [tuple(p) for p in NP[reg]["P"]] for reg in ("sub1", "sub2")}

check("x-degree of P is 8 in both subcases",
      all(max(i for i, _ in CORNERS[r]) == 8 for r in ("sub1", "sub2")))
check("C4-forcing corners (8,14),(8,16) present in both subcases",
      all({(8, 14), (8, 16)} <= set(CORNERS[r]) for r in ("sub1", "sub2")))
M_SUB1_DEG = vmax(-1, 1, CORNERS["sub1"])   # direction (-1,1): j - i
M_SUB2_DEG = vmax(-2, 1, CORNERS["sub2"])   # direction (-2,1): j - 2i
M_ORD = {r: vmax(2, -1, CORNERS[r]) for r in ("sub1", "sub2")}  # 2i - j
check("sub1 deg direction: max (j-i) over corners = 8", M_SUB1_DEG == 8)
check("sub2 deg direction: max (j-2i) over corners = 0", M_SUB2_DEG == 0)
check("ord direction: max (2i-j) over corners = 2, both subcases",
      M_ORD["sub1"] == 2 and M_ORD["sub2"] == 2)
check("sub1 deg direction attained at BOTH (8,16) and (0,8) (T3 remark)",
      {(8, 16), (0, 8)} <= {c for c in CORNERS["sub1"]
                            if -c[0] + c[1] == M_SUB1_DEG})

# ------------------------------------------- W1. the inductions, symbolic in k
if not QUIET:
    print("W1. the three T3 §3 inductions close for ALL k (symbolic)")
k, j = symbols("k j")

# Each induction: hypothesis h(m) bounds the direction-valuation of C_{4-m};
# recursion C_{4-k} = -(1/2C4) * (P_{8-k} + sum_{j=1..k-1} C_{4-j}C_{4-k+j}).
# v(C4^-1) = -v(C4) in each direction; slice bound = M - a*(8-k) for the
# y-part of P_{8-k} under direction (a,b=1) (resp. -ord under (2,-1)).
#   name          a   b   v(C4)          h(k)
#   sub1 deg     -1   1    8             8 - k
#   sub2 deg     -2   1    8             8 - 2k
#   ord (both)    2  -1   -7 (= -ord)    2k - 7
V_C4 = {"sub1": 8, "sub2": 8, "ord": -7}
H = {"sub1": 8 - k, "sub2": 8 - 2 * k, "ord": 2 * k - 7}
A_DIR = {"sub1": -1, "sub2": -2, "ord": 2}
M_DIR = {"sub1": M_SUB1_DEG, "sub2": M_SUB2_DEG, "ord": M_ORD["sub1"]}

for name in ("sub1", "sub2", "ord"):
    h = H[name]
    prod = expand(h.subs(k, j) + h.subs(k, k - j))          # v(C_{4-j}C_{4-(k-j)})
    slc = M_DIR[name] - A_DIR[name] * (8 - k)               # v of P_{8-k} y-part
    check(f"{name}: product bound is j-free and equals the slice bound "
          f"(exact closure, T3)", expand(prod - slc) == 0)
    check(f"{name}: induction step  -v(C4) + slice = h(k)  identically in k",
          expand(-V_C4[name] + slc - h) == 0)
    check(f"{name}: base k=1 (empty product sum) gives h(1)",
          (-V_C4[name] + slc.subs(k, 1)) == h.subs(k, 1))
# The slices P_{8-k} for k=6,7,8 are the genuine slices P_2,P_1,P_0 (8-k>=0):
check("k=6,7,8 use genuine polygon slices (8-k >= 0)", all(8 - K >= 0 for K in (6, 7, 8)))

# --------------------------------------- W2. D-transform arithmetic, symbolic
if not QUIET:
    print("W2. D-transform: caps 15k / 14k / 12k, symbolic + literal rows")
jx = symbols("jx")  # x-exponent of C_j; window variable d_{4-k} has jx = 4-k
# deg C_j <= jx+4 (sub1), <= 2jx (sub2); ord C_j >= 2jx-1; D_j = C_j*C4^(7-2jx)
check("sub1: deg D_j <= (jx+4) + 8(7-2jx) = 60 - 15jx",
      expand((jx + 4) + 8 * (7 - 2 * jx) - (60 - 15 * jx)) == 0)
check("sub2: deg D_j <= 2jx + 8(7-2jx) = 56 - 14jx",
      expand(2 * jx + 8 * (7 - 2 * jx) - (56 - 14 * jx)) == 0)
check("ord:  ord D_j >= (2jx-1) + 7(7-2jx) = 48 - 12jx",
      expand((2 * jx - 1) + 7 * (7 - 2 * jx) - (48 - 12 * jx)) == 0)
check("substituting jx = 4-k: caps become deg <= 15k / 14k, ord >= 12k",
      expand((60 - 15 * jx).subs(jx, 4 - k) - 15 * k) == 0
      and expand((56 - 14 * jx).subs(jx, 4 - k) - 14 * k) == 0
      and expand((48 - 12 * jx).subs(jx, 4 - k) - 12 * k) == 0)
ROWS = {6: (72, 90, 84), 7: (84, 105, 98), 8: (96, 120, 112)}
for K, (o, d1_, d2_) in ROWS.items():
    check(f"k={K} literal row: ord>={o}, deg<={d1_} (sub1) / {d2_} (sub2)",
          12 * K == o and 15 * K == d1_ and 14 * K == d2_)

# --------------------------- W3. the d3-killing shift in D-coordinates
if not QUIET:
    print("W3. shift x -> x - c3/(4 C4): polynomial D-coordinate identity")
u, s = symbols("u s")
cs = {m: symbols(f"c{m}") for m in range(-4, 4)}
cs[4] = symbols("c4")
# (a) series recomposition: coefficient of x^j in C(x - s) equals
#     sum_{m=j..4} binom(m, m-j) c_m (-s)^(m-j).   (x = 1/u; exact per j.)
U_shift = sum(cs[m] * u**(-m) * sp.series((1 - s * u)**m, u, 0, 9).removeO()
              for m in range(-4, 5))
U_shift = expand(U_shift * u**4)  # now coefficient of u^(4-j) is c~_j
for jv in range(3, -5, -1):
    formula = sum(binomial(m, m - jv) * cs[m] * (-s)**(m - jv)
                  for m in range(jv, 5))
    check(f"recomposition: c~_{jv} = sum binom(m,m-{jv}) c_m (-s)^(m-{jv})",
          expand(U_shift.coeff(u, 4 - jv) - formula) == 0)
# (b,c) in D-coordinates all C4 powers cancel: with c_m = D_m C4^(2m-7),
#     s = c_3/(4C4) = D_3/(4 C4^2)  (the unique s with c~_3 = 0):
C4s = symbols("C4s")
D = {m: symbols(f"D{m}") for m in range(-4, 4)}
D[4] = 1
sub_cd = {cs[m]: D[m] * C4s**(2 * m - 7) for m in range(-4, 4)}
sub_cd[cs[4]] = C4s  # c_4 = C4  <->  D_4 = 1
s_val = D[3] / (4 * C4s**2)
for jv in range(3, -5, -1):
    ctil = sum(binomial(m, m - jv) * cs[m] * (-s)**(m - jv) for m in range(jv, 5))
    ctil = ctil.subs(sub_cd).subs(s, s_val)
    Dtil = sum(binomial(m, m - jv) * D[m] * (-D[3] / 4)**(m - jv)
               for m in range(jv, 5))
    check(f"C4-exponent cancellation: c~_{jv} * C4^(7-2*{jv}) = D~_{jv} "
          f"(polynomial in the D's)",
          sp.simplify(expand(ctil * C4s**(7 - 2 * jv)) - expand(Dtil)) == 0)
Dtil3 = sum(binomial(m, m - 3) * D[m] * (-D[3] / 4)**(m - 3) for m in range(3, 5))
check("the shift kills the k=1 variable: D~_3 = D_3 + 4*(-D_3/4) = 0",
      expand(Dtil3) == 0)
# (d) cap preservation, term by term: deg(D_m (-D_3/4)^(m-j))
#     <= (cap at m) + (m-j)*(cap at 3) = cap at j, since cap(3) = slope:
for capname, cap in (("sub1", 60 - 15 * jx), ("sub2", 56 - 14 * jx)):
    step3 = cap.subs(jx, 3)
    check(f"{capname}: deg[D_m*(-D_3/4)^(m-j)] <= cap(m)+(m-j)*cap(3) = cap(j)",
          expand(cap.subs(jx, k) + (k - j) * step3 - cap.subs(jx, j)) == 0)
ordcap = 48 - 12 * jx
check("ord:  ord[D_m*(-D_3/4)^(m-j)] >= ordcap(m)+(m-j)*ordcap(3) = ordcap(j)",
      expand(ordcap.subs(jx, k) + (k - j) * ordcap.subs(jx, 3)
             - ordcap.subs(jx, j)) == 0)

# ------------------------- W4. end-to-end generic random-QQ corroboration
if not QUIET:
    print("W4. generic end-to-end: hull-supported random P, D-recursion to D_-4")


def hull_chains(corners):
    """Andrew monotone chain; returns (lower, upper) vertex chains."""
    pts = sorted(set(corners))
    def half(pl):
        out = []
        for p in pl:
            while len(out) >= 2 and (
                (out[-1][0] - out[-2][0]) * (p[1] - out[-2][1])
                - (out[-1][1] - out[-2][1]) * (p[0] - out[-2][0])) <= 0:
                out.pop()
            out.append(p)
        return out
    return half(pts), half(pts[::-1])  # lower, upper


def hull_j_range(corners, i):
    """[min j, max j] of the hull at abscissa i (exact, Fractions)."""
    lower, upper = hull_chains(corners)
    def interp(chain, i):
        vals = []
        for (x0, y0), (x1, y1) in zip(chain, chain[1:]):
            if min(x0, x1) <= i <= max(x0, x1) and x0 != x1:
                vals.append(Fraction(y0) + Fraction(y1 - y0, x1 - x0) * (i - x0))
            elif x0 == i:
                vals.append(Fraction(y0))
        if chain and chain[-1][0] == i:
            vals.append(Fraction(chain[-1][1]))
        return vals
    allv = interp(lower, i) + interp(upper, i)
    lo, hi = min(allv), max(allv)
    import math
    return math.ceil(lo), math.floor(hi)


def d_recursion(Pslices):
    """verify_derivation.py §C, division-free: D_3 = P_7/2;
    D_k = 1/2 P_{k+4} C4^(6-2k) - 1/2 sum_{pairs} D_i D_j  (i,j <= 3)."""
    Dv = {}
    for kk in range(3, -5, -1):
        pairs = [(i, kk + 4 - i) for i in range(kk + 1, 4)
                 if kk + 4 - i <= 3 and kk + 4 - i >= i]
        acc = Rational(1, 2) * Pslices.get(kk + 4, sp.Integer(0)) * C4**(6 - 2 * kk)
        for (i, jj2) in pairs:
            mult = 2 if i != jj2 else 1
            acc -= Rational(mult, 2) * Dv[i] * Dv[jj2]
        Dv[kk] = expand(acc)
    Dv[4] = sp.Integer(1)
    return Dv


def degord(p):
    P_ = Poly(p, y)
    return P_.degree(), min(m[0] for m in P_.monoms())


for reg, degcap in (("sub1", lambda m: 60 - 15 * m), ("sub2", lambda m: 56 - 14 * m)):
    rng = random.Random(72108 if reg == "sub1" else 108072)
    corners = CORNERS[reg]
    Pslices = {8: expand(C4**2)}
    d8, o8 = degord(Pslices[8])
    lo8, hi8 = hull_j_range(corners, 8)
    check(f"{reg}: P_8 = C4^2 support [{o8},{d8}] matches hull at i=8 "
          f"[{lo8},{hi8}]", (o8, d8) == (lo8, hi8))
    for i in range(0, 8):
        lo, hi = hull_j_range(corners, i)
        Pslices[i] = sum(rng.choice([-9, -7, -5, -3, -1, 1, 2, 3, 5, 7, 9])
                         * y**m for m in range(lo, hi + 1))
    Dv = d_recursion(Pslices)
    for kk in range(2, -5, -1):        # window k = 4-kk runs 2..8
        K = 4 - kk
        dg, od = degord(Dv[kk])
        check(f"{reg}: unshifted d_(4-{K}) = D_{kk}: deg {dg} <= {degcap(kk)}, "
              f"ord {od} >= {48 - 12 * kk}",
              dg <= degcap(kk) and od >= 48 - 12 * kk)
    # the shift, in D-coordinates (W3 identity), numeric:
    Dt = {}
    for jv in range(3, -5, -1):
        Dt[jv] = expand(sum(binomial(m, m - jv) * Dv[m] * (-Dv[3] / 4)**(m - jv)
                            for m in range(jv, 5)))
    check(f"{reg}: shifted D~_3 = 0 exactly", Dt[3] == 0)
    for jv in range(2, -5, -1):
        K = 4 - jv
        dg, od = degord(Dt[jv])
        check(f"{reg}: SHIFTED d~_(4-{K}): deg {dg} <= {degcap(jv)}, "
              f"ord {od} >= {48 - 12 * jv}  (k={K} cap holds after shift)",
              dg <= degcap(jv) and od >= 48 - 12 * jv)
    kk_attained = sum(1 for kk in range(2, -5, -1)
                      if degord(Dt[kk])[0] == degcap(kk))
    if not QUIET:
        print(f"    ({reg}: shifted deg cap attained on {kk_attained}/7 rows "
              f"— tightness, informational)")

# --------------------------------------------- W5. consumer cross-checks
if not QUIET:
    print("W5. the caps the consumers actually use")
sys.path.insert(0, HERE)
import full_system_bridge as fsb  # noqa: E402

check("full_system_bridge.WEIGHT is 12k for k=2..8 and 204 for Phi",
      fsb.WEIGHT == {"d2": 24, "d1": 36, "d0": 48, "dm1": 60,
                     "dm2": 72, "dm3": 84, "dm4": 96, "Phi": 204})
check("full_system_bridge.STRIP_DEGCAP is (15k-12k)=3k / (14k-12k)=2k",
      fsb.STRIP_DEGCAP == {"sub1": {"dm2": 18, "dm3": 21, "dm4": 24},
                           "sub2": {"dm2": 12, "dm3": 14, "dm4": 16}})
check("spare-unknown totals: sub1 19+22+25=66, sub2 13+15+17=45",
      sum(v + 1 for v in fsb.STRIP_DEGCAP["sub1"].values()) == 66
      and sum(v + 1 for v in fsb.STRIP_DEGCAP["sub2"].values()) == 45)
import re  # noqa: E402
jl = open(os.path.join(HERE, "jetlift.py")).read()
sizes2 = re.search(r"'f31_sub2'.*?sizes=\[([\d,]+)\]", jl).group(1)
sizes1 = re.search(r"'f31_sub1'.*?sizes=\[([\d,]+)\]", jl).group(1)
check("jetlift CONFIGS window sizes are 2k+1 (sub2) / 3k+1 (sub1), k=2..5",
      [int(x) for x in sizes2.split(",")] == [2 * K + 1 for K in range(2, 6)]
      and [int(x) for x in sizes1.split(",")] == [3 * K + 1 for K in range(2, 6)])

# tightness witness: Phi has weight 204 = 12*17, i.e. k = 17; re-derive f1
# from its forcing ODE (verify_derivation.py §A) and confirm Phi = f1*C4^28
# ATTAINS deg = 14*17 = 238 and ord = 12*17 = 204 (sub2 caps are sharp).
a = symbols("a0:16")
ansatz = sum(a[i] * y**i for i in range(16))
eqs = Poly(expand(8 * y * (y + 1) * sp.diff(ansatz, y) - 14 * (8 * y + 7) * ansatz
                  - y**8 * (y + 1)**2), y).all_coeffs()
sol = sp.solve(eqs, a, dict=True)
check("f1 forcing ODE has a unique polynomial solution", len(sol) == 1)
Phi = expand(ansatz.subs(sol[0]) * C4**28)
dgP, odP = degord(Phi)
check("Phi = f1*C4^28 attains the k=17 sub2 caps: deg 238 = 14*17, "
      "ord 204 = 12*17 (tightness witness)",
      dgP == 238 == 14 * 17 and odP == 204 == 12 * 17)

print(f"\nALL {ok[0]} WINDOW-CAP CHECKS PASSED"
      + ("" if QUIET else "  (k=6,7,8 caps: PROVEN from the audited premises;"
         " no NEW conditionality vs the k=2..5 window)"))
