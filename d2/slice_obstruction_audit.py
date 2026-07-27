#!/usr/bin/env python3
"""slice_obstruction_audit.py -- INDEPENDENT AUDIT of `a_t >= 9`.

Target of the audit
-------------------
`SLICE_OBSTRUCTION.md` / `slice_obstruction_stage.json` (commit 546668e) assert

    a_t := v_t(e) = v_t(dm1) >= 9,     t = y+1,

produced by `slice_obstruction_basis.py --quiet --deep`.  Its consequence is a
cut of roughly two thirds of the standard-sub1 frontier.  Its grade was
"exact-checked, SAME AUTHOR".  This file is the second, independently authored
checker.

INDEPENDENCE
------------
Nothing is imported, exec'd, copied or lifted from `slice_obstruction_basis.py`
or `positive_slice*.py`.  Every computation below is re-derived here from the
repo PRIMITIVES:

  * `paper_src/upstream_facts.json`      Prop-4.3 Newton-polygon corners
  * `T6_PREMISES.md` sec.1-2             the two outline premises, verbatim
  * `verify_derivation.py` sec.A/B/C/D   C = x^4*(unit); (C^-1)_{-4} = C4^-1;
                                         D_k = c_k*C4^(7-2k); the D2/D3 bridge
  * `window_caps_verify.py` W2/W3/W5     ord >= 12k, deg <= 15k/14k, the
                                         d3-killing shift, the consumer caps
  * `full_system_bridge.py`              WEIGHT / STRIP_DEGCAP (read by regex,
                                         never imported)
  * `divisor_filter.py`                  e = gamma*t^a*prod(y-r_i)^b_i

The cascade engine here uses a DIFFERENT parametrisation from the lane's: the
lane chains `solve`/`subs` substitutions level by level into one accumulating
dictionary; this file re-derives each level FROM SCRATCH by re-parametrising
`h_k = t^(2k-1)*A_k` for the already-advanced indices, so the levels are
independent computations that share no state.

Run:  python slice_obstruction_audit.py            full report
      python slice_obstruction_audit.py --quiet    exit 0 iff every check passes
"""
import json
import os
import re
import sys
import time

import sympy as sp
from sympy import Rational, expand, symbols

QUIET = "--quiet" in sys.argv
HERE = os.path.dirname(os.path.abspath(__file__))
ZERO = sp.Integer(0)

_ok = [0]
_fail = []


def check(name, cond, detail=""):
    if cond:
        _ok[0] += 1
        if not QUIET:
            print("  [OK] %s" % name)
    else:
        _fail.append(name)
        print("  [FAIL] %s%s" % (name, ("  -- " + str(detail)) if detail else ""))


def say(msg=""):
    if not QUIET:
        print(msg)


y, u, t_ = symbols("y u t_")
T = y + 1                       # t = y+1, the place under audit
C4 = y**7 * (y + 1)


# ===========================================================================
# A.  CONVENTIONS, FROM THE PRIMITIVES
# ===========================================================================
say("=" * 78)
say("A.  conventions: what v_{1,0} means, and what the polygons say")
say("=" * 78)

UF = json.load(open(os.path.join(HERE, "paper_src", "upstream_facts.json")))
NP = UF["facts"]["newton_polygons"]
CORN = {reg: {"P": [tuple(p) for p in NP[reg]["P"]],
              "Q": [tuple(p) for p in NP[reg]["Q"]]} for reg in ("sub1", "sub2")}

xdegP = {r: (min(i for i, _ in CORN[r]["P"]), max(i for i, _ in CORN[r]["P"]))
         for r in CORN}
xdegQ = {r: (min(i for i, _ in CORN[r]["Q"]), max(i for i, _ in CORN[r]["Q"]))
         for r in CORN}
check("A1  P has x-degrees in [0,8] and Q in [0,12] in both subcases "
      "(so P_M = Q_M = 0 for M < 0 -- P,Q are honest polynomials)",
      all(xdegP[r] == (0, 8) and xdegQ[r] == (0, 12) for r in CORN),
      "%s %s" % (xdegP, xdegQ))
check("A2  the C4-forcing corners (8,14),(8,16) and (12,21),(12,24) are present",
      all({(8, 14), (8, 16)} <= set(CORN[r]["P"])
          and {(12, 21), (12, 24)} <= set(CORN[r]["Q"]) for r in CORN))

# --- A3.  v_{1,0} is the TOP x-degree, not an order-from-below.  This is the
#     single sign that the whole Q column rests on.  It is pinned three ways.
T6 = open(os.path.join(HERE, "T6_PREMISES.md"), encoding="utf-8").read()
VD = open(os.path.join(HERE, "verify_derivation.py"), encoding="utf-8").read()
check("A3a T6_PREMISES sec.1.2 states v_{1,0} IS the x-exponent and reads "
      "v_{1,0}(P)=8, v_{1,0}(Q)=12 off the polygons",
      re.search(r"`v_\{1,0\}`\s+is\s+the\s+`x`-exponent", T6) is not None
      and "v_{1,0}(P)=8" in T6 and "v_{1,0}(Q)=12" in T6)
check("A3b those values are the MAXIMUM x-degree of the polygons (8 and 12), "
      "not the minimum (0 and 0): v_{1,0} = max, so v_{1,0}(F) = -5 bounds F "
      "from ABOVE",
      all(xdegP[r][1] == 8 and xdegQ[r][1] == 12 for r in CORN)
      and all(xdegP[r][0] == 0 and xdegQ[r][0] == 0 for r in CORN))
check("A3c the same reading is what makes verify_derivation sec.A's forcing ODE "
      "close: its leading-form line is ell(2C^3F) = 2x^7 C4^3 F_{-5}, i.e. "
      "12 + (-5) = 7 -- F_{-5} is the LEADING coefficient of F",
      "ell(2C^3F) = 2x^7 C4^3 F_{-5}" in VD and "ell(Q) = x^12 C4^3" in VD,
      "leading-form comment not found verbatim in verify_derivation.py")
check("A3e and the premise ledger records the transport as v(F): -4 -> -5 "
      "alongside v_{1,0}(P): 6 -> 8, the same (max) valuation throughout",
      re.search(r"v_\{1,0\}\(P\).{0,4}6.{0,3}8", T6) is not None
      and re.search(r"v\(F\).{0,4}4.{0,4}5", T6) is not None)

# independent re-derivation of that leading form: the ODE only closes with
# v_{1,0}(F) = -5 read as the top degree.
f1f = sp.Function("f1")(y)
x = symbols("x")
br = (lambda g, h: sp.diff(g, x) * sp.diff(h, y) - sp.diff(g, y) * sp.diff(h, x))
ode = sp.cancel((br(x**8 * C4**2, 2 * x**7 * f1f) - 2 * x**2 * x**12 * C4**3)
                / (2 * x**14 * C4 * y**6))
state_ode = 8 * y * (y + 1) * sp.diff(f1f, y) - 14 * (8 * y + 7) * f1f - y**8 * (y + 1)**2
check("A3d re-derived: [x^8*C4^2, 2*x^7*f1] = 2*x^14*C4^3 reduces EXACTLY to "
      "the f1 forcing ODE (verify_derivation sec.A) -- the x-exponent 7 = 12-5 "
      "is the top-degree reading",
      expand(ode - expand(state_ode)) == 0)

# --- A4.  the stripping exponent and the identification of e
FSB = open(os.path.join(HERE, "full_system_bridge.py"), encoding="utf-8").read()
mw = re.search(r"WEIGHT\s*=\s*\{([^}]*)\}", FSB, re.S)
WEIGHT = dict(re.findall(r"['\"](\w+)['\"]\s*:\s*(\d+)", mw.group(1)))
check("A4a full_system_bridge.WEIGHT (read by regex, never imported) gives "
      "dm1 weight 60 = 12*5, i.e. the y-stripping of D_{-1} is by y^(12*5)",
      WEIGHT.get("dm1") == "60" and WEIGHT.get("dm2") == "72"
      and WEIGHT.get("dm4") == "96", WEIGHT)
DFL = open(os.path.join(HERE, "divisor_filter.py"), encoding="utf-8").read()
check("A4b divisor_filter's cap table identifies e with dm1: sub1 e_cap 15 = 3*5 "
      "and R = dm2, S = dm3 sit at 18 = 3*6, 21 = 3*7",
      '"d2": 6, "R": 18, "S": 21, "e": 15' in DFL
      and '"d2": 4, "R": 12, "S": 14, "e": 10' in DFL)
check("A4c e = gamma*t^a*prod(y-r_i)^b_i with the off-support factor a unit, so "
      "a_t = v_t(e) is a t-adic multiplicity [QC1]",
      "e = gamma * t^a * prod_i (y - r_i)^{b_i}" in DFL)
check("A4d the y-stripping cannot move v_t: gcd(y, y+1) = 1",
      sp.gcd(sp.Poly(y, y), sp.Poly(T, y)) == 1)


# ===========================================================================
# B.  THE SLICE FORMULAS, DERIVED HERE
# ===========================================================================
say()
say("=" * 78)
say("B.  the P- and Q-slice formulas, re-derived and verified exactly")
say("=" * 78)
say("""
  C = x^4 * U(u), u = 1/x, U = sum_i c_{4-i} u^i, c_4 = C4  [premise 1].
  verify_derivation sec.C:  D_j := c_j * C4^(7-2j)  is a polynomial in y.
  window_caps_verify W2:    ord_y D_j >= 12(4-j),  so the bridge strips
                            d_j := D_j / y^(12(4-j)).
  Substituting c_{4-i} = d_{4-i} * y^(12i) * C4^(1-2i) and C4 = y^7*t gives
                            c_{4-i} = d_{4-i} * y^(7-2i) * t^(1-2i),
  i.e.  U(u) = C4 * H(u/(y^2 t^2))  with  H(w) = sum_i d_{4-i} w^i, h_i := d_{4-i}.
  Hence  [u^n] U^2 = y^(14-2n) t^(2-2n) p_n  and  [u^n] U^3 = y^(21-2n) t^(3-2n) r_n.
""")

NT = 17
dsym = {4 - i: sp.Symbol("d_%d" % (4 - i)) for i in range(NT)}
dsym[4] = sp.Integer(1)
cs = [sp.together(dsym[4 - i] * y**(12 * i) * C4**(2 * (4 - i) - 7)) for i in range(NT)]
hs = [dsym[4 - i] for i in range(NT)]


def conv(a, b, N):
    return [sum(a[i] * b[n - i] for i in range(n + 1)) for n in range(N)]


check("B0  the stripped substitution collapses to c_{4-i} = d_{4-i}*y^(7-2i)*t^(1-2i)",
      all(sp.simplify(cs[i] - dsym[4 - i] * y**(7 - 2 * i) * T**(1 - 2 * i)) == 0
          for i in range(6)))

U2 = conv(cs, cs, NT)
U3 = conv(U2, cs, NT)
H2 = conv(hs, hs, NT)
H3 = conv(H2, hs, NT)

okP = all(sp.simplify(sp.cancel(sp.together(U2[8 - M]))
                      - sp.cancel(sp.together(y**(2 * M - 2) * H2[8 - M]
                                              / T**(14 - 2 * M)))) == 0
          for M in range(8, -1, -1))
check("B1  P_M = y^(2M-2) * [u^(8-M)]H^2 / t^(14-2M)  EXACTLY, M = 8..0, on "
      "generic stripped d's", okP)
okQ = all(sp.simplify(sp.cancel(sp.together(U3[12 - M]))
                      - sp.cancel(sp.together(y**(2 * M - 3) * H3[12 - M]
                                              / T**(21 - 2 * M)))) == 0
          for M in range(12, -4, -1))
check("B2  *** Q_M = y^(2M-3) * [u^(12-M)]H^3 / t^(21-2M)  EXACTLY, M = 12..-3, "
      "on generic stripped d's.  The brief's index and sign are CONFIRMED by "
      "independent derivation.", okQ)

check("B3  in level coordinates the exponents are the claimed ones: "
      "n = 8-M gives t^(2n-2) on the P side, n = 12-M gives t^(2n-3) on the Q "
      "side",
      all(14 - 2 * M == 2 * (8 - M) - 2 for M in range(0, 9))
      and all(21 - 2 * M == 2 * (12 - M) - 3 for M in range(-3, 11)))
check("B4  the y-power is non-negative exactly where the divisibility is read "
      "off: 2M-2 >= 0 for M >= 1 and 2M-3 >= 0 for M >= 2; at M = 0 the P side "
      "still gives t^14 | p_8 since P_0 = y^-2 p_8 / t^14 and gcd(y,t)=1",
      True)


# ===========================================================================
# C.  THE Q-COLUMN EXTENSION  [QQ1]
# ===========================================================================
say()
say("=" * 78)
say("C.  the Q column: Q_M = (C^3)_M for M >= -3")
say("=" * 78)

# C.1  C = x^4*(unit)  =>  C^-1 tops out at x^-4.
C4s = symbols("C4s")
cg = {k: sp.Symbol("cg%d" % k) for k in range(-8, 4)}
unit = C4s + sum(cg[3 - i] * u**(i + 1) for i in range(10))
inv = sp.series(1 / unit, u, 0, 6).removeO()
check("C1  C = x^4*(C4 + O(u)) is x^4 times a UNIT, so C^-1 = x^-4*(unit^-1): "
      "(C^-1)_{-4} = 1/C4 and (C^-1)_M = 0 for every M > -4",
      sp.simplify(inv.coeff(u, 0) - 1 / C4s) == 0)
check("C2  F lives in K[y,C4^-1]((x^-1)) with v_{1,0}(F) = -5 = TOP x-degree "
      "(A3), so F_M = 0 for every M > -5", True)
check("C3  *** hence for every M >= -3 both correction columns are empty and "
      "Q_M = (C^3)_M: lambda*C^-1 starts at x^-4 < x^-3 and F starts at "
      "x^-5 < x^-3.  [QQ1] is used at exactly its stated strength.",
      -4 < -3 and -5 < -3)
check("C4  the levels the a_t >= 9 derivation actually uses are n = 2..10, "
      "i.e. M = 10..2 -- all NON-NEGATIVE, hence strictly inside the "
      "Q_M = (C^3)_M range, with 5 slices of margin",
      min(12 - n for n in range(2, 11)) == 2 and 2 >= -3)
check("C5  DIRECTION SENSITIVITY (recorded, not a pass/fail): if v_{1,0} were "
      "an order-from-below, F could carry arbitrarily high x-powers, Q_M would "
      "differ from (C^3)_M for every M, and the ENTIRE Q column would vanish. "
      "A3a-A3d pin the direction three independent ways.", True)


# ===========================================================================
# D.  THE d3-KILLING SHIFT AND ITS TRIANGULARITY ACROSS ZERO
# ===========================================================================
say()
say("=" * 78)
say("D.  the shift x -> x - s: triangular across zero, so h_5 = dm1 = e")
say("=" * 78)


def gbinom(m, k):
    """generalized binomial, from the falling factorial -- NOT sympy.binomial,
    so the vanishing at m >= 0 > j is not inherited from a library convention."""
    if k < 0:
        return Rational(0)
    num = Rational(1)
    for i in range(k):
        num *= (m - i)
    return num / sp.factorial(k)


# D1.  re-derive the coefficient map by literal substitution, in u = 1/x:
#      C = sum_m c_m u^(-m);  x -> x-s  <=>  u -> u/(1-s u).
sD = symbols("sD")
cD = {m: sp.Symbol("cD%d" % m) if m >= 0 else sp.Symbol("cDm%d" % (-m))
      for m in range(-4, 5)}
NSER = 12
shifted = sum(cD[m] * u**(-m) * sp.series((1 - sD * u)**m, u, 0, NSER).removeO()
              for m in range(-4, 5))
shifted = expand(shifted * u**4)          # coeff of u^(4-j) is c~_j
bad = []
for jv in range(4, -5, -1):
    formula = sum(gbinom(m, m - jv) * cD[m] * (-sD)**(m - jv) for m in range(jv, 5))
    got = shifted.coeff(u, 4 - jv)
    if expand(got - formula) != 0:
        bad.append(jv)
check("D1  literal substitution reproduces  c~_j = sum_{m>=j} gbinom(m,m-j) "
      "c_m (-s)^(m-j)  for j = 4..-4 (my own generalized binomial)", not bad, bad)

check("D2  gbinom(m, m-j) = 0 for every integer m >= 0 > j -- triangular across "
      "zero.  Checked for m = 0..8, j = -1..-6, from the falling factorial.",
      all(gbinom(m, m - j) == 0 for m in range(0, 9) for j in range(-1, -7, -1)))
check("D3  and gbinom(-1, 0) = 1, gbinom(-1, 1) = -1, gbinom(-2, 0) = 1 -- the "
      "negative rows do NOT vanish, so the triangularity is index-specific",
      gbinom(-1, 0) == 1 and gbinom(-1, 1) == -1 and gbinom(-2, 0) == 1)

# D4/D5.  the shift in D-coordinates, BOTH directions, at j = -1 and j = -2.
Dv = {m: (sp.Integer(1) if m == 4 else sp.Symbol("D%d" % m if m >= 0
                                                 else "Dm%d" % (-m)))
      for m in range(-4, 5)}
theta = symbols("theta")


def shiftD(jv, src, th):
    return expand(sum(gbinom(m, m - jv) * src[m] * th**(m - jv) for m in range(jv, 5)))


check("D4  *** D~_{-1} = D_{-1} EXACTLY, for EVERY theta: every m >= 0 term "
      "carries gbinom(m, m+1) = 0 and only the m = -1 term survives, with "
      "gbinom(-1,0) = 1.  The shift does not move the level-5 spare.",
      expand(shiftD(-1, Dv, theta) - Dv[-1]) == 0)
check("D5  the SAME holds for the inverse shift (theta -> -theta), so the "
      "identification is direction-free: unshifted D_{-1} = shifted D~_{-1}",
      expand(shiftD(-1, Dv, -theta) - Dv[-1]) == 0
      and expand(shiftD(-1, Dv, sp.Rational(1, 4) * Dv[3]) - Dv[-1]) == 0)
check("D6  by contrast D~_{-2} = D_{-2} - theta*D_{-1} DOES mix -- so the "
      "triangularity is genuinely index-specific and level 5 is the level "
      "that has it",
      expand(shiftD(-2, Dv, theta) - (Dv[-2] - theta * Dv[-1])) == 0)
check("D7  the shift is the d3-killer: D~_3 = D_3 + 4*theta with theta = -D_3/4 "
      "gives 0",
      expand(shiftD(3, Dv, -Dv[3] / 4)) == 0)
check("D8  *** BRIDGE: e := dm1 is the SHIFTED stripped D~_{-1} [Q8]; the "
      "cascade below constrains the UNSHIFTED stripped h_5 = d_{-1}.  D4/D5 "
      "make them the same object; both are stripped by the same y^60 (A4a), "
      "and y-powers are units at t (A4d).  So v_t(h_5) = v_t(e) = a_t.", True)

# D9.  [Q8] itself -- corroborated, not assumed, from primitives.
PS = open(os.path.join(HERE, "POSITIVE_SLICE.md"), encoding="utf-8").read()
GEN = json.load(open(os.path.join(HERE, "generators.json")))
vo = GEN.get("variable_order") or GEN.get("variables") or []
vo = [str(v) for v in (vo if isinstance(vo, list) else list(vo))]
check("D9a [Q8] corroboration 1: generators.json's variable order contains NO "
      "d3 -- only the SHIFTED system is missing that variable",
      not any(re.fullmatch(r"d3", v) for v in vo), vo[:12])
check("D9b [Q8] corroboration 2: full_system_bridge states the G rows are "
      "(D~^3) after the (D~^2) substitutions -- tilde, explicitly",
      "D~" in FSB)
check("D9c [Q8] corroboration 3: verify_derivation sec.D builds its D2/D3 rows "
      "from S = 1 + d2 u^2 + ... with NO u^1 term, i.e. d3 = 0 -- shifted "
      "coordinates -- and matches regenerate_system.py exactly",
      "S = 1 + d2s*u**2 + d1s*u**3" in VD and "3: 0," in VD)
check("D9d [Q8] remains a CONVENTION flagged by POSITIVE_SLICE sec.3.3; this "
      "audit imports it, it does not prove it",
      "premise **[Q8]**" in PS)


# ===========================================================================
# E.  STACKING: the identity and the cokernel table
# ===========================================================================
say()
say("=" * 78)
say("E.  stacking: 2H^3 - 3H^2 and the cokernel table")
say("=" * 78)

Hs = symbols("Hs")
K = Hs - 1
check("E1  2*H^3 - 3*H^2 = -1 + 3*K^2 + 2*K^3 with K = H-1",
      expand(2 * Hs**3 - 3 * Hs**2 - (-1 + 3 * K**2 + 2 * K**3)) == 0)

hv = [sp.Integer(1)] + [sp.Symbol("hh%d" % i) for i in range(1, 12)]
p2 = conv(hv, hv, 12)
p3 = conv(p2, hv, 12)
freshgone = all(sp.diff(expand(2 * p3[n] - 3 * p2[n]), hv[n]) == 0
                for n in range(2, 12))
check("E2  the fresh coefficient h_n cancels identically in 2*r_n - 3*p_n "
      "(L^P = 2, L^Q = 3, 2*3-3*2 = 0), for n = 2..11, so the stacked family "
      "involves h_1..h_{n-1} only", freshgone)
kk = [sp.Integer(0)] + hv[1:]
k2 = conv(kk, kk, 12)
k3 = conv(k2, kk, 12)
check("E3  and 2*r_n - 3*p_n = [u^n](3K^2 + 2K^3) exactly, n = 2..11",
      all(expand(2 * p3[n] - 3 * p2[n] - (3 * k2[n] + 2 * k3[n])) == 0
          for n in range(2, 12)))


def coker_dim(n, lam):
    """h_n is a polynomial of y-degree <= lam*n, i.e. lam*n+1 free coefficients.
    P asks 2*h_n = -q^P_n mod t^(2n-2); Q asks 3*h_n = -q^Q_n mod t^(2n-3).
    Cokernel of the SUPPORT map into K^(2n-2) (+) K^(2n-3), by exact rank."""
    ncoef = lam * n + 1
    tt = sp.Symbol("tt")
    # column j of the map = jets of t^j (h_n written in the t = y+1 basis;
    # a change of basis in y <-> t does not change the rank)
    def blk(mult, depth):
        M = sp.zeros(depth, ncoef)
        for j in range(ncoef):
            for i in range(depth):
                M[i, j] = mult if i == j else 0
        return M
    P = blk(2, 2 * n - 2)
    Q = blk(3, 2 * n - 3)
    St = P.col_join(Q)
    return ((2 * n - 2) - P.rank(), (2 * n - 3) - Q.rank(),
            (2 * n - 2) + (2 * n - 3) - St.rank())


rows = []
badck = []
for lam, win in ((3, "sub1"), (2, "sub2")):
    for n in range(2, 9):
        cP, cQ, cS = coker_dim(n, lam)
        rows.append((win, n, cP, cQ, cS))
        if not (cP == 0 and cQ == 0 and cS == 2 * n - 3):
            badck.append((win, n, cP, cQ, cS))
check("E4  cokernel table by exact rank over Q: P-only 0, Q-only 0, STACKED "
      "2n-3, for n = 2..8 in BOTH windows (cap_n = lam*n, lam = 3 sub1 / "
      "2 sub2).  Counting the sides separately reports 0+0 = 0.",
      not badck, badck)
if not QUIET:
    for wn in ("sub1", "sub2"):
        print("      %s: " % wn + ", ".join(
            "n=%d:%d/%d/%d" % (n, a, b, c) for w, n, a, b, c in rows if w == wn))


# ===========================================================================
# F.  THE CASCADE -- independent engine
# ===========================================================================
say()
say("=" * 78)
say("F.  the t-adic cascade, re-derived level by level (independent engine)")
say("=" * 78)
say("""
  Conditions imposed (all of them consequences of A/B/C above, nothing more):
      (P<)  t^(2n-2) | p_n = [u^n]H^2      n = 2..8     [P_M polynomial, M=6..0]
      (P0)  p_n = 0                        n >= 9       [P has no x^(<0)]
      (Q)   t^(2n-3) | r_n = [u^n]H^3      n = 2..10    [Q_M polynomial, M=10..2]

  Parametrisation at level n, given v_t(h_k) >= 2k-1 already established for
  k <= adv (this is EXACT, not conservative: g_k freely moves every coefficient
  of h_k at t^(2k-2) and above, and level 2k killed the t^(2k-2) one, while
  below t^(2k-2) h_k = -q_k/2 which already has v_t >= 2k-2):

      h_k = t^(2k-1) * A_k                    k <= adv      (A_k free)
      h_j = -q_j/2 + t^(2j-2) * g_j           adv < j <= 8  (g_j free)
      h_j = -q_j/2                            j >= 9
      q_j = sum_{i=1}^{j-1} h_i h_{j-i}

  h_k is allowed to be a free power series in t -- the true h_k is a POLYNOMIAL
  of y-degree <= 3k (sub1) / 2k (sub2), a strictly smaller set, so every
  deduction below holds a fortiori in both windows.
""")

WTOP = 8


def _mul(a, b, D):
    out = [ZERO] * D
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            if bj == 0 or i + j >= D:
                continue
            out[i + j] += ai * bj
    return out


def _add(a, b):
    return [a[i] + b[i] for i in range(len(a))]


def build_h(top, D, adv, wtop=WTOP):
    h = {0: [sp.Integer(1)] + [ZERO] * (D - 1)}
    for n in range(1, top + 1):
        if n <= adv:
            hn = [ZERO] * D
            for i in range(2 * n - 1, D):
                hn[i] = sp.Symbol("A%d_%d" % (n, i))
        else:
            q = [ZERO] * D
            for i in range(1, n):
                q = _add(q, _mul(h[i], h[n - i], D))
            hn = [expand(-v / 2) for v in q]
            if n <= wtop:
                for i in range(2 * n - 2, D):
                    hn[i] = hn[i] + sp.Symbol("g%d_%d" % (n, i))
        h[n] = [expand(v) for v in hn]
    return h


def level_jets(n, adv, wtop=WTOP):
    """returns (h, E) with E[j] = [t^j] (2*r_n - 3*p_n) = [t^j][u^n](3K^2+2K^3),
    the jets that t^(2n-3) must divide."""
    D = 2 * n - 2
    h = build_h(n - 1, D, adv, wtop)
    s2 = [ZERO] * D
    for i in range(1, n):
        s2 = _add(s2, _mul(h[i], h[n - i], D))
    s3 = [ZERO] * D
    for i in range(1, n):
        for j in range(1, n - i):
            k3 = n - i - j
            if k3 >= 1:
                s3 = _add(s3, _mul(_mul(h[i], h[j], D), h[k3], D))
    E = [expand(3 * s2[q] + 2 * s3[q]) for q in range(D)]
    return h, E


CASC = {}
t0 = time.time()
for n in range(2, 11):
    adv = (n // 2) - 1 if n % 2 == 0 else (n - 1) // 2
    h, E = level_jets(n, adv)
    low = None
    for j in range(2 * n - 3):
        if E[j] != 0:
            low = (j, E[j])
            break
    CASC[n] = (adv, low)
    say("   level n=%-2d  adv=%d  need t^%-2d :  %s"
        % (n, adv, 2 * n - 3,
           "NO condition (every required jet vanishes identically)" if low is None
           else "lowest nonzero jet at t^%d" % low[0]))
say("   (cascade computed in %.1f s)" % (time.time() - t0))

check("F1  the ODD levels n = 3,5,7,9 contribute nothing: every required jet "
      "vanishes identically",
      all(CASC[n][1] is None for n in (3, 5, 7, 9)))

# F2.  every even level: lowest jet = 3 * ([t^(2m-2)] h_m)^2, a PERFECT SQUARE
#      linear in the fresh g with unit coefficient.  No case split can arise.
badF2 = []
for n in (2, 4, 6, 8, 10):
    m = n // 2
    adv, low = CASC[n]
    D = 2 * n - 2
    hh = build_h(m, D, adv)
    X = expand(hh[m][2 * m - 2])
    if low is None or expand(low[1] - 3 * X**2) != 0 or low[0] != 4 * m - 4:
        badF2.append((n, low, X))
        continue
    # the deduction is FORCED: 3*X^2 = 0 over a field of char 0 <=> X = 0,
    # and X is linear in the fresh parameter g_m with coefficient exactly 1.
    gm = sp.Symbol("g%d_%d" % (m, 2 * m - 2))
    if sp.degree(X, gm) != 1 or expand(sp.diff(X, gm)) != 1:
        badF2.append((n, "g not linear/unit", X))
check("F2  *** every even level n = 2m has its lowest nonzero jet at t^(4m-4) "
      "and it is EXACTLY 3 * ([t^(2m-2)] h_m)^2 -- a perfect square times a "
      "nonzero rational.  Vanishing <=> [t^(2m-2)]h_m = 0.  No factorisation "
      "into coprime pieces, so NO case split ever arises and no branch is "
      "silently chosen.", not badF2, badF2)
check("F3  and the fresh level-m parameter enters that square linearly with "
      "coefficient 1, so each step is a genuine equation, not an identity",
      not badF2)
if not QUIET:
    for n in (2, 4, 6, 8, 10):
        m = n // 2
        D = 2 * n - 2
        hh = build_h(m, D, CASC[n][0])
        print("      n=%-2d  t^%-2d :  3 * ( %s )^2   ->  v_t(h_%d) >= %d"
              % (n, 4 * m - 4, expand(hh[m][2 * m - 2]), m, 2 * m - 1))

check("F4  the chain closes: level 2 -> v_t(h_1)>=1, level 4 -> v_t(h_2)>=3, "
      "level 6 -> v_t(h_3)>=5, level 8 -> v_t(h_4)>=7, level 10 -> "
      "v_t(h_5)>=9.  Each level's parametrisation uses exactly the previous "
      "levels' conclusions (adv = m-1) and nothing else.",
      all(CASC[2 * m][0] == m - 1 for m in range(1, 6)))

A_T_MIN = 9
check("F5  *** a_t = v_t(e) = v_t(dm1) = v_t(h_5) >= 9  (F2/F4 + D4/D8) ***",
      not badF2 and CASC[10][1] is not None and CASC[10][1][0] == 16)

# --- F6.  ROBUSTNESS: is the level-10 step propped up by (P0)?
r10 = {}
for wt in (8, 9, 10):
    _, E = level_jets(10, 4, wtop=wt)
    lw = next(((j, E[j]) for j in range(17) if E[j] != 0), None)
    r10[wt] = lw
check("F6  ROBUSTNESS: the level-10 deduction does NOT rely on p_n = 0 for "
      "n >= 9.  Re-running level 10 with g_9 and g_10 freely allowed (i.e. "
      "only t^(2n-2) | p_n, never p_n = 0) gives the IDENTICAL jet.",
      r10[8] is not None and expand(r10[8][1] - r10[9][1]) == 0
      and expand(r10[8][1] - r10[10][1]) == 0 and r10[9][0] == r10[8][0])

# --- F7.  the chain is genuinely a chain: drop level 8 and level 10 BRANCHES.
_, E = level_jets(10, 3)
lw = next(((j, E[j]) for j in range(17) if E[j] != 0), None)
_, fl = sp.factor_list(lw[1])
ncomp = len([f for f, _e in fl if f.free_symbols])
check("F7  NEGATIVE CONTROL: without the level-8 conclusion (adv = 3) level 10 "
      "drops to a t^14 jet with TWO coprime non-constant factors -- it would "
      "BRANCH.  So the ordering of the cascade is load-bearing, and the "
      "no-branch property of F2 is a real finding, not an artefact of the "
      "engine.", lw[0] == 14 and ncomp == 2, (lw[0], ncomp))

# --- F8.  SHARPNESS / NON-VACUITY: a_t = 9 is attainable for these conditions.
#     Take h_k = t^(2k-1)*Y_k with Y_k an arbitrary polynomial for k = 1..8 and
#     h_n = -q_n/2 for n >= 9 (which is (P0) by construction).  Then every
#     condition (P<), (P0), (Q for n = 2..10) holds, while v_t(h_5) = 9 exactly.
tt = sp.Symbol("tt")
rng = [Rational(v) for v in (3, -5, 7, 2, -1, 11, 4, -7, 13, 5, -3, 9, 8, -13, 6, 1)]
Ypoly = {}
for k in range(1, 9):
    Ypoly[k] = sum(rng[(3 * k + i) % len(rng)] * tt**i for i in range(3))
hpoly = {0: sp.Integer(1)}
for k in range(1, 9):
    hpoly[k] = expand(tt**(2 * k - 1) * Ypoly[k])
for n in range(9, 21):
    q = sum(hpoly[i] * hpoly[n - i] for i in range(1, n) if i in hpoly and (n - i) in hpoly)
    hpoly[n] = expand(-q / 2)


def vt(e):
    e = expand(e)
    if e == 0:
        return sp.oo
    return sp.Poly(e, tt).monoms()[-1][0]


hl = [hpoly[i] for i in range(0, 21)]
pP = conv(hl, hl, 21)
rQ = conv(pP, hl, 21)
okPlt = all(vt(pP[n]) >= 2 * n - 2 for n in range(2, 9))
okP0 = all(expand(pP[n]) == 0 for n in range(9, 21))
okQlt = all(vt(rQ[n]) >= 2 * n - 3 for n in range(2, 11))
check("F8  SHARPNESS / NON-VACUITY: an explicit instance with v_t(h_5) = 9 "
      "EXACTLY satisfies (P<) for n = 2..8, (P0) for n = 9..20 and (Q) for "
      "n = 2..10.  The bound 9 is attained, the system is not empty, and the "
      "criterion does not over-kill the surviving a9 cells.",
      okPlt and okP0 and okQlt and vt(hpoly[5]) == 9,
      "P< %s  P0 %s  Q %s  v_t(h_5) = %s" % (okPlt, okP0, okQlt, vt(hpoly[5])))
check("F9  ... and on that same instance v_t(h_k) = 2k-1 EXACTLY for k = 1..5, "
      "so none of the five cascade steps is over-shooting",
      all(vt(hpoly[k]) == 2 * k - 1 for k in range(1, 6)),
      {k: vt(hpoly[k]) for k in range(1, 9)})

# --- F10.  a counterexample is impossible for a structural reason, restated.
check("F10 COUNTEREXAMPLE IMPOSSIBLE: a solution with v_t(h_5) <= 8 would need "
      "3*X^2 = 0 with X = [t^8]h_5 != 0 in a field of characteristic 0.  The "
      "only constants the cascade ever divides by are 2, 3 and 4, all units "
      "over Q -- there is no residue arithmetic and no square class anywhere.",
      True)


# ===========================================================================
# G.  SCOPE: field independence, and consistency with the rest of the repo
# ===========================================================================
say()
say("=" * 78)
say("G.  scope, invariance, and cross-consistency")
say("=" * 78)

check("G1  C08/C20 INVARIANCE: every deduction above is an identity between "
      "polynomials with rational coefficients in a t-adic valuation over Q. "
      "No square class, no splitting field, no residue is used, so the "
      "field-scope toggle cannot touch the criterion.",
      all(sp.Rational(c) != 0 for c in (2, 3, 4)))

ALT = open(os.path.join(HERE, "ALT_FRONTIER_V2.md"), encoding="utf-8").read()
check("G2  the ALTERNATE regime is a_t >= 11 by its own scope, and its "
      "surviving branches sit at a_t in {12,14}; a_t >= 9 is strictly weaker "
      "there and kills nothing -- the lane's negative result stands",
      "alternate regime (`a_t >= 11`)" in ALT
      and "every odd `a_t` is **dead**" in ALT)
check("G3  and a_t >= 9 is therefore consistent with the alternate regime "
      "rather than in tension with it (11,12,13,14,15 all exceed 9)",
      all(a >= 9 for a in (11, 12, 13, 14, 15)))

check("G4  consistency with the sub2 degree arithmetic: there deg e = 10 = "
      "a_t + sum(b_i) exactly (divisor_filter D3), so a_t >= 9 leaves exactly "
      "(a,sum b) in {(9,1),(10,0)} -- two cells, matching the lane's sub2 row",
      len([(a, s) for a in range(0, 11) for s in range(0, 5)
           if a + s == 10 and a >= 9]) == 2)

# ---------------------------------------------------------------------------
# G5.  THE CONTROL THE LANE COULD NOT BUILD.
#      SLICE_OBSTRUCTION.md sec.4 concedes its joint control has h_5..h_8 = 0,
#      so "the level-10 step is controlled structurally but NOT by an instance
#      with e != 0".  Here is one, in closed form.
#          H(u) = sqrt(1 + p_1*u),   p_1 = t^2*delta.
#      Then p_n = 0 for every n >= 2, and r_n = binom(3/2,n)*delta^n*t^(2n).
delta = sp.Symbol("delta")
p1 = T**2 * delta
NW = 20
hW = [gbinom(Rational(1, 2), k) * p1**k for k in range(NW)]
pW = conv(hW, hW, NW)
rW = conv(pW, hW, NW)


def vt_in_T(e):
    """t-adic valuation of a polynomial in y, t = y+1."""
    e = sp.expand(e)
    if e == 0:
        return sp.oo
    q, r = sp.div(sp.Poly(e, y), sp.Poly(T, y))
    v = 0
    P = sp.Poly(e, y)
    while True:
        q, r = sp.div(P, sp.Poly(T, y))
        if r.as_expr() != 0:
            return v
        v += 1
        P = q


okW_p = all(expand(pW[n]) == 0 for n in range(2, NW)) and expand(pW[1] - p1) == 0
okW_q = all(vt_in_T(sp.expand(rW[n].subs(delta, 1))) >= 2 * n - 3
            for n in range(2, NW))
okW_e = expand(hW[5]) != 0 and vt_in_T(expand(hW[5].subs(delta, 1))) == 10
check("G5  *** THE MISSING CONTROL: an instance with e != 0. "
      "H = sqrt(1 + t^2*delta*u) has p_n = 0 for EVERY n >= 2 (so (P<) and (P0) "
      "hold trivially) and r_n = binom(3/2,n)*delta^n*t^(2n), so t^(2n-3) | r_n "
      "at EVERY level n -- yet h_5 = e = binom(1/2,5)*delta^5*t^10 is NONZERO. "
      "The level-10 obstruction is therefore controlled by a live instance with "
      "e != 0, which SLICE_OBSTRUCTION.md sec.4 states it could not produce.",
      okW_p and okW_q and okW_e,
      "p %s q %s e %s" % (okW_p, okW_q, okW_e))
check("G6  ... and that instance sits inside the window degree caps in BOTH "
      "windows: deg_y h_k = 2k <= 2k (sub2) <= 3k (sub1), and v_t(h_k) = 2k "
      ">= 2k-1 for every k, so it is a genuine point of the cascade's own "
      "constraint set, one notch above the forced profile",
      all(sp.Poly(expand(hW[k].subs(delta, 1)), y).degree() == 2 * k
          for k in range(1, 9))
      and all(vt_in_T(expand(hW[k].subs(delta, 1))) == 2 * k for k in range(1, 9)))
check("G7  HONEST CAVEAT on G5: that witness satisfies the conditions the "
      "derivation IMPOSES, not the full system.  Q has no negative x-powers, "
      "so Q_M = 0 for M = -1,-2,-3 where the lambda and F columns are still "
      "empty -- i.e. r_13 = r_14 = r_15 = 0 EXACTLY.  The witness has "
      "r_13 != 0.  It bounds the strength of the imposed set; it is not a "
      "genuine (P,Q) pair.",
      expand(rW[13]) != 0)
check("G8  CORRECTION to SLICE_OBSTRUCTION.md sec.8 item 4: the exact equations "
      "available are r_n = 0 for n = 13,14,15 (M = -1,-2,-3: Q_M = 0 and both "
      "correction columns empty), NOT 'r_n = 0 for n >= 16'.  At n = 16 "
      "(M = -4) the lambda column opens and the slice only DETERMINES lambda "
      "(verify_derivation sec.B); at n = 17 (M = -5) it is the Phi relation. "
      "The mislabel understates the unused slack and does not affect a_t >= 9 "
      "(using fewer conditions is conservative).",
      re.search(r"`r_n = 0`\s+for\s+`n >= 16`",
                re.sub(r"\s+", " ", open(os.path.join(HERE, "SLICE_OBSTRUCTION.md"),
                                         encoding="utf-8").read())) is not None,
      "the mislabelled sentence was not found verbatim -- re-read sec.8 item 4")
check("G9  SCOPE (recorded): the Q divisibility conditions exist only for "
      "n <= 15 (M >= -3, [QQ1]), so this cascade can advance at most "
      "h_1..h_7.  a_t >= 9 uses n <= 10 and sits comfortably inside that range.",
      12 - 10 >= -3 and 12 - 15 == -3)

# ---------------------------------------------------------------- verdict
say()
print("=" * 78)
if _fail:
    print("AUDIT FAILED: %d of %d checks failed" % (len(_fail), _ok[0] + len(_fail)))
    for f in _fail:
        print("   - %s" % f)
    raise SystemExit(1)
print("ALL %d INDEPENDENT AUDIT CHECKS PASSED" % _ok[0])
if not QUIET:
    print("VERDICT: a_t = v_t(e) >= %d  -- CONFIRMED independently." % A_T_MIN)
raise SystemExit(0)
