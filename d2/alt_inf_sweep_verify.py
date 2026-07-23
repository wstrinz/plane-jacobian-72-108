#!/usr/bin/env python3
"""Independent spot-checker for alt_inf_sweep.py (companion: ALT_INF_SWEEP.md).

This file does NOT import alt_inf_sweep.  It re-derives the flipped max-plus
degree chain from first principles: h_f degrees are recomputed directly from
the monomial table (cascade_signature.load_levels), and the (D_t) cofactors
r_f are, for the concrete kill windows, built as EXACT rational functions in
sympy (mirroring alt_regime_inf_verify.py) so the degree bookkeeping is
grounded in real polynomial degrees, not just tropical arithmetic.

It checks:
  A. A T1 top-anchor kill: 2*deg d1 < w forces deg r_6 < 0 (r_6 = h_7/T is not
     a polynomial).  [a=12, deg d1 = 2, w = 6]
  B. Two bottom-close kills, re-derived by an explicit hand-style degree chain
     with unique maxima at every level (so NO leading cancellation is possible,
     the drop rule being gated on a tie), ending in 21 deg_E + H_0 != 4 + R_0.
     Each is confirmed against an exact sympy window whose closing residual
     E^21 h_0 + u r_0 is nonzero of the dominant degree.
  C. One OPEN-branch survivor: an explicit witness chain (with its recorded
     drops) is checked against every level identity (I7)/(If)/(I0), including
     that each realized drop sits on a genuine tie and the close is a tie.

Exact integer / rational arithmetic only.  PASSES when the sweep is sound.
"""
from __future__ import annotations

import re
from pathlib import Path

import sympy as sp

import cascade_engine as ce
import cascade_signature as cs

ROOT = Path(__file__).resolve().parent
NEG_INF = ce.NEG_INF
DEG_U = ce.DEG_U

# --- symbols / u exactly as in alt_regime_inf_verify.py -------------------
d2s, d1s, d0s, dm1s, y = sp.symbols("d2 d1 d0 dm1 y")
t = y + 1
qpoly = 2048*y**4 - 512*y**3 + 320*y**2 - 240*y + 195
u_expr = sp.Rational(-1, 6630) * qpoly
assert int(sp.degree(sp.Poly(u_expr, y))) == DEG_U

# h_f in source variables (d2,d1,d0,dm1) from f31_graded.txt.
_txt = (ROOT / "f31_graded.txt").read_text(encoding="utf-8")
_pat = r"h_(\d) \(weight \d+, dm1-power \d+\) = (.+)"
H_SRC = {int(m.group(1)): sp.sympify(m.group(2))
         for m in re.finditer(_pat, _txt)}
assert sorted(H_SRC) == list(range(8))

LEVELS = cs.load_levels()   # for the independent monomial degree recompute


# ---------------------------------------------------------------------------
# Independent tropical max degree of h_f from the monomial table.
# ---------------------------------------------------------------------------
def h_deg_indep(f: int, st: tuple) -> tuple[float, int]:
    """(max degree, #achievers) of h_f under degree state st=(dd2,dd1,dsig,de).

    Recomputed straight from MONOMIALS with the zero-variable drops -- no call
    into the sweep or into deg_h_options.
    """
    dd2, dd1, dsig, de = st
    sigma_zero, d2_zero, d1_zero = dsig == NEG_INF, dd2 == NEG_INF, dd1 == NEG_INF
    best, ach = NEG_INF, 0
    for (k, x, z, ee), _c in ce.MONOMIALS[f]:
        if d1_zero and x:
            continue
        if sigma_zero and z:
            continue
        if d2_zero and k:
            continue
        val = ee * de
        ok = True
        for exp, dv in ((k, dd2), (x, dd1), (z, dsig)):
            if exp:
                if dv == NEG_INF:
                    ok = False
                    break
                val += exp * dv
        if not ok:
            continue
        if val > best:
            best, ach = val, 1
        elif val == best:
            ach += 1
    return best, ach


def window(a: int, st: tuple, seed: int = 7):
    """Concrete sympy window realizing degree state st (exact rationals).

    d2 -> deg dd2 poly (0 if d2_zero), sigma -> deg dsig poly (0 if sigma_zero),
    d1 -> deg dd1 poly (0 if d1_zero / T2), E -> deg (de-a) unit poly at t=y+1,
    e = t^a E, d0 = (sigma + d2^2)/4.
    """
    import random
    random.seed(seed)

    def rp(dg):
        if dg == NEG_INF:
            return sp.Integer(0)
        return sp.expand(y**dg + sum(random.randint(-2, 2)*y**j
                                     for j in range(dg))) if dg > 0 else sp.Integer(3)

    dd2, dd1, dsig, de = st
    D2 = rp(dd2)
    D1 = rp(dd1)
    S = rp(dsig)
    dE = de - a
    E = rp(dE)
    if E.subs(y, -1) == 0:
        E = E + 1
    e_full = sp.expand(t**a * E)
    d0_expr = sp.expand((S + D2**2) / 4)
    subs = {d2s: D2, d1s: D1, d0s: d0_expr, dm1s: e_full}
    hval = {f: sp.expand(H_SRC[f].subs(subs)) for f in range(8)}
    return hval, E, e_full


def deg_rat(expr) -> float:
    expr = sp.cancel(sp.together(expr))
    num, den = sp.fraction(expr)
    if sp.expand(num) == 0:
        return NEG_INF
    return int(sp.degree(sp.Poly(num, y))) - int(sp.degree(sp.Poly(den, y)))


PASS = True


def check(label: str, cond: bool) -> None:
    global PASS
    PASS = PASS and cond
    print(f"    [{'PASS' if cond else 'FAIL'}] {label}")


# ===========================================================================
print("PART 1. Ground tropical h-degrees against exact sympy polynomials.")
# a generic full-support window at a=12 (matches ALT_REGIME_INF worked table).
gst = (6, 9, 12, 15)
hval, E, e_full = window(12, gst)
for f in range(8):
    Hi, _ = h_deg_indep(f, gst)
    Hs = deg_rat(hval[f])
    check(f"deg h_{f}: tropical {int(Hi)} == sympy {int(Hs)} (<= {60-6*f})",
          Hi == Hs and Hi <= 60 - 6*f)


# ===========================================================================
print("\nPART 2 (kill A). T1 top anchor 2*deg d1 < w  =>  deg r_6 < 0.")
aA, dd1A = 12, 2
wA = 3*aA - 30
stA = (NEG_INF, dd1A, NEG_INF, aA)          # deg_E = 0 window, T1
hvalA, _, _ = window(aA, stA)
H7 = deg_rat(hvalA[7])
TA = sp.expand(t**wA)
r6 = sp.cancel(hvalA[7] / TA)               # r_6 = h_7 / T
_, denom = sp.fraction(sp.together(r6))
print(f"    w = 3a-30 = {wA};  H_7 = deg h_7 = 2*deg d1 = {int(H7)}")
print(f"    deg r_6 = H_7 - w = {int(H7)} - {wA} = {int(H7)-wA}")
check("2*deg d1 = 4 < w = 6", 2*dd1A < wA)
check("H_7 - w < 0", H7 - wA < 0)
check("r_6 = h_7/T is NOT a polynomial (t divides denominator)",
      sp.rem(sp.Poly(denom, y), sp.Poly(t, y)) == 0 or denom != 1)
check("deg_rat(r_6) < 0 (no polynomial r_6 exists)", deg_rat(r6) < 0)


# ===========================================================================
def hand_close_kill(a: int, st: tuple, label: str):
    """Rigid hand chain: unique maxima at every level => forced R_f, then the
    bottom close has a unique maximum (term1 != term2) => kill.  Confirmed on
    an exact sympy window."""
    print(f"\n{label}")
    w = 3*a - 30
    dE = st[3] - a
    print(f"    a={a}, w={w}, deg_E={dE}, state (dd2,dd1,dsig,de)={st}")
    # anchor
    H7, ac7 = h_deg_indep(7, st)
    check(f"anchor h_7 unique achiever (#={ac7})", ac7 == 1 and H7 != NEG_INF)
    R = H7 - w
    print(f"    (I7)  w + R_6 = H_7 = {int(H7)}  ->  R_6 = {int(R)}")
    check("R_6 >= 0", R >= 0)
    for f in range(6, 0, -1):
        Hf, acf = h_deg_indep(f, st)
        term1 = NEG_INF if Hf == NEG_INF else 3*(7-f)*dE + Hf
        term2 = 4 + R
        mx = max(term1, term2)
        # rigid requirement: strict unique max, and if it is the h-term the
        # h-achiever must be unique (else a tie-drop would be available).
        rigid = (term1 != term2) and not (term1 > term2 and acf > 1)
        check(f"(I{f}) level {f}: unique max "
              f"max({term1},{term2})={mx}, forced (rigid={rigid})", rigid)
        R = mx - w
        print(f"          w + R_{f-1} = {int(mx)}  ->  R_{f-1} = {int(R)}")
        check(f"R_{f-1} >= 0", R >= 0)
    H0, ac0 = h_deg_indep(0, st)
    term1 = 21*dE + H0
    term2 = 4 + R
    print(f"    (I0)  close: term1 = 21*deg_E + H_0 = {21*dE}+{int(H0)} = "
          f"{int(term1)};  term2 = 4 + R_0 = {int(term2)}")
    check(f"h_0 unique achiever (#={ac0}) so NO drop can lower term1", ac0 == 1)
    check("closing terms differ (unique maximum) => sum has nonzero leading "
          "term => E^21 h_0 + u r_0 != 0 : KILLED",
          term1 != term2)
    # exact sympy confirmation on a concrete window.
    hval, Ewin, _ = window(a, st)
    T = sp.expand(t**w)
    r = {6: sp.cancel(hval[7] / T)}
    for f in range(6, 0, -1):
        r[f-1] = sp.cancel((Ewin**(3*(7-f))*hval[f] + u_expr*r[f]) / T)
    residual = sp.expand(Ewin**21 * hval[0] + u_expr*r[0])
    dres = deg_rat(residual)
    check(f"sympy window: deg r_0 = {int(deg_rat(r[0]))} matches hand R_0 = "
          f"{int(R)}", deg_rat(r[0]) == R)
    check(f"sympy window: closing residual != 0, deg = {int(dres)} = "
          f"max(term1,term2) = {int(max(term1,term2))} (nonzero leading term)",
          sp.expand(residual) != 0 and dres == max(term1, term2))


hand_close_kill(12, (NEG_INF, 3, NEG_INF, 12),
                "PART 3 (kill B). Bottom-close kill, a=12 T1, sigma=0 & d2=0.")
hand_close_kill(11, (NEG_INF, 2, NEG_INF, 11),
                "PART 4 (kill C). Bottom-close kill, a=11 T1, sigma=0 & d2=0.")


# ===========================================================================
print("\nPART 5 (open survivor). Explicit witness chain satisfies every "
      "identity and the close ties.")
# independent forward reachability with parent pointers (own implementation).
def witness(a: int, st: tuple):
    w = 3*a - 30
    dE = st[3] - a
    d1_zero = st[1] == NEG_INF
    # per-level H options with #achievers (drops allowed only when >1 achiever)
    Hopt = {}
    for f in range(8):
        Hm, ac = h_deg_indep(f, st)
        Hopt[f] = (Hm, ac)
    # reach[f] : R_f -> (prevR_{f+1}, H_used, kind, depth)
    H7, ac7 = Hopt[7]
    if H7 == NEG_INF:
        reach = {NEG_INF: None}
    else:
        R6 = H7 - w
        if R6 < 0:
            return None
        reach = {R6: (None, H7, "max", 0)}
    chains = {6: reach}
    for f in range(6, 0, -1):
        Hm, ac = Hopt[f]
        gshift = 3*(7-f)*dE
        # h option values: max, and (if tie) drops down to 0 and NEG_INF
        hopts = [(Hm, "max", 0)]
        if ac > 1 and Hm != NEG_INF:
            for hv in range(int(Hm)-1, -1, -1):
                hopts.append((hv, "h_drop", int(Hm)-hv))
            hopts.append((NEG_INF, "h_vanish", 0))
        elif Hm == NEG_INF:
            hopts = [(NEG_INF, "max", 0)]
        new = {}
        for Rf, _p in chains[f].items():
            term2 = NEG_INF if Rf == NEG_INF else 4 + Rf
            for hv, hk, hd in hopts:
                term1 = NEG_INF if hv == NEG_INF else gshift + hv
                tie = term1 == term2
                mx = max(term1, term2)
                if mx == NEG_INF:
                    new.setdefault(NEG_INF, (Rf, hv, hk, hd))
                    continue
                base = mx - w
                if not tie:
                    if base >= 0:
                        new.setdefault(int(base), (Rf, hv, hk, hd))
                    continue
                if base >= 0:
                    new.setdefault(int(base), (Rf, hv, hk, hd))
                    for kk in range(0, int(base)):
                        new.setdefault(kk, (Rf, hv, "sum_drop", int(base)-kk))
                new.setdefault(NEG_INF, (Rf, hv, "sum_vanish", 0))
        chains[f-1] = new
    # close
    H0, ac0 = Hopt[0]
    h0opts = [(H0, 0)]
    if ac0 > 1 and H0 != NEG_INF:
        h0opts += [(hv, int(H0)-hv) for hv in range(int(H0)-1, -1, -1)]
        h0opts += [(NEG_INF, 0)]
    for R0, _p in chains[0].items():
        term2 = NEG_INF if R0 == NEG_INF else 4 + R0
        for hv, hd in h0opts:
            term1 = NEG_INF if hv == NEG_INF else 21*dE + hv
            if term1 == term2:
                return (a, st, w, dE, chains, R0, hv, hd, term1)
    return None


sst = (5, 2, 10, 11)          # survivor in branch a11_b0000_T1
res = witness(11, sst)
check("witness chain exists for survivor state (5,2,10,11)", res is not None)
if res:
    a, st, w, dE, chains, R0, H0u, H0d, closeval = res
    # rebuild the explicit R-chain by backtracking parents.
    order = []
    R = R0
    for f in range(0, 7):
        order.append((f, R, chains[f][R]))
        R = chains[f][R][0]
    order.reverse()   # from f=6 up... rebuild identity checks per level
    # (I7)
    H7, ac7 = h_deg_indep(7, st)
    R6 = chains[6]
    R6val = next(iter(R6))
    print(f"    (I7) w+R_6 = {w}+{R6val} = {w+R6val} = H_7 = {int(H7)}")
    check("(I7) top anchor identity", w + R6val == H7)
    # walk levels f=6..1 verifying (If)
    for f in range(6, 0, -1):
        Rf = None
        # find R_f used on the path: it is chains[f-1][R_{f-1}][0]
        # reconstruct forward from R6
        pass
    # simpler: reconstruct path forward
    path = {6: R6val}
    # choose the R_{f-1} on the path leading to R0: backtrack from R0
    cur = R0
    path[0] = R0
    for f in range(0, 6):
        prev = chains[f][cur][0]
        path[f+1] = prev
        cur = prev
    Hopt = {f: h_deg_indep(f, st) for f in range(8)}
    okall = True
    for f in range(6, 0, -1):
        Rf = path[f]
        Rprev = path[f-1]
        prevR, hused, hkind, hdepth = chains[f-1][Rprev]
        Hm, ac = Hopt[f]
        gshift = 3*(7-f)*dE
        term1 = NEG_INF if hused == NEG_INF else gshift + hused
        term2 = NEG_INF if Rf == NEG_INF else 4 + Rf
        mx = max(term1, term2)
        realized = w + Rprev
        tie = term1 == term2
        ok = (realized <= mx) and (realized == mx or tie)
        if hkind in ("h_drop", "h_vanish"):
            ok = ok and ac > 1              # h-drop needs a real h-tie
        print(f"    (I{f}) w+R_{f-1} = {w}+{Rprev} = {realized}; "
              f"max(term1={term1},term2={term2})={mx}; H_{f} used={hused}"
              f"({hkind})")
        okall = okall and ok
    check("(If) every level identity holds with drops gated on ties", okall)
    print(f"    (I0) close: 21*deg_E + H_0(used {H0u}) = {21*dE+ (0 if H0u==NEG_INF else H0u)}"
          f" = 4 + R_0 = {4+R0} : TIE")
    check("(I0) closing tie 21 deg_E + H_0 == 4 + R_0", 21*dE + H0u == 4 + R0)


# ===========================================================================
print(f"\n{'ALL ALT_INF_SWEEP SPOT CHECKS PASS' if PASS else 'SOME CHECKS FAILED'}")
raise SystemExit(0 if PASS else 1)
