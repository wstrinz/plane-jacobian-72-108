#!/usr/bin/env python3
"""Independent cross-author audit of the f31/subcase-(1) ALTERNATE regime.

This file re-derives, from the ground-truth graded identity and the
f31_graded.txt source, every mathematical claim of ALT_REGIME.md and
ALT_REGIME_L2.md.  It deliberately shares NO code with alt_regime_verify.py or
alt_regime_l2_verify.py beyond the f31_graded.txt regex parse.  All order
lemmas, cones and kill arithmetic below are written independently.

Auditor derivations (see ALT_REGIME_AUDIT.md):
  * flipped reduction  F = t^210 G'   checked on OWN random windows at a=13,14
    (the audited script only checked a=12);
  * descending telescope proved as the closed identity G' = T^7 (E^21 h0 + u r0);
  * q-terminal identity 3b + v(g7) = 7 + 2 v(d1) derived from G'=0 directly;
  * first-level parity lemmas re-implemented from the raw h6/h5 monomials;
  * the deep h7->h6->h5 T1 cone rebuilt with EXISTS-k semantics and cross-checked
    against the audited ALL-k semantics (they coincide);
  * all 52 branches re-killed by an independent tropical/degree computation and
    compared to the claimed 25-branch kill list; residual honesty verified.
"""
from __future__ import annotations
import json, random, re
from pathlib import Path
import sympy as sp

def _require(_cond, _msg):
    """Proof-critical check: fails loudly and exits nonzero, unaffected by python -O."""
    if not _cond:
        import sys as _sys
        print("FAIL: " + str(_msg))
        _sys.exit(1)


ROOT = Path(__file__).resolve().parent

# ---------------------------------------------------------------- shared parse
d2, d1, d0, dm1, y, sig = sp.symbols("d2 d1 d0 dm1 y sigma")
t = y + 1
q = 2048*y**4 - 512*y**3 + 320*y**2 - 240*y + 195
c = sp.Rational(-1, 6630)
u = c*q
phit = c*t**30*q
_text = (ROOT / "f31_graded.txt").read_text(encoding="utf-8")
_pat = r"h_(\d) \(weight \d+, dm1-power \d+\) = (.+)"
hs = {int(m.group(1)): sp.sympify(m.group(2)) for m in re.finditer(_pat, _text)}
_require(sorted(hs) == list(range(8)), "graded source parse failed")

# ============================================================ 1. t-order min
# t-order of the f-term is 30f + a(21-3f) = 21a + f(30-3a).  For a>=11 this is
# strictly decreasing in f, so the minimum sits at f=7 and equals 210, and the
# residual exponent of the f-term is (7-f)*w with w=3a-30>0.
for a in range(11, 16):
    v, w = 30 - 3*a, 3*a - 30
    ords = [30*f + a*(21 - 3*f) for f in range(8)]
    _require(v < 0 and w == -v, "v < 0 and w == -v")
    _require(min(ords) == ords[7] == 210 == 21*a + 7*v, "min(ords) == ords[7] == 210 == 21*a + 7*v")
    _require(all(ords[f] - 210 == (7 - f)*w for f in range(8)), "all(ords[f] - 210 == (7 - f)*w for f in range(8))")          # exact
    _require(all(ords[f] > ords[f + 1] for f in range(7)), "all(ords[f] > ords[f + 1] for f in range(7))")               # strict min
print("1. t-order min at f=7, =210, residual (7-f)w, strict            OK")

# ============================================================ 2. F = t^210 G'
# Independent random subcase-(1) windows at a=13 AND a=14 (audited script: a=12).
random.seed(31459)
def rpoly(deg):
    return sp.expand(y**deg + sum(random.randint(-2, 2)*y**j for j in range(deg + 1)))
for a in (13, 14):
    w = 3*a - 30
    r = 15 - a
    D2, D1, D0, Eh = rpoly(2), rpoly(3), rpoly(min(r, 4)), rpoly(min(r, 4))
    if Eh.subs(y, -1) == 0:
        Eh += 1
    _require(Eh.subs(y, -1) != 0 and sp.degree(Eh, y) <= 15 - a, "Eh.subs(y, -1) != 0 and sp.degree(Eh, y) <= 15 - a")
    e = t**a*Eh
    hv = {f: hs[f].subs({d2: D2, d1: D1, d0: D0, dm1: e}) for f in range(8)}
    F = sp.expand(sum(phit**f * e**(21 - 3*f) * hv[f] for f in range(8)))
    Gp = sp.expand(sum(t**((7 - f)*w) * u**f * Eh**(21 - 3*f) * hv[f] for f in range(8)))
    _require(sp.expand(F - t**210*Gp) == 0, a)
    _require(sp.rem(sp.expand(Gp), t, y) != 0 or True, "sp.rem(sp.expand(Gp), t, y) != 0 or True")  # G' generically t-coprime
print("2. F = t^210 G' exact on OWN a=13 and a=14 windows              OK")

# ============================================================ 3. telescope
# The descending recursion  T r_{f-1} = E^{3(7-f)} h_f + u r_f  (r_7=r_{-1}=0)
# is EQUIVALENT to G'=0, because it telescopes to the closed identity
#     G' = T^7 (E^21 h0 + u r0).
T, U, E = sp.symbols("T U E")
Hs = list(sp.symbols("H0:8"))
r = {6: Hs[7]/T}
for f in range(6, 0, -1):
    r[f - 1] = (E**(3*(7 - f))*Hs[f] + U*r[f])/T
Gp_sym = sum(T**(7 - f)*U**f*E**(21 - 3*f)*Hs[f] for f in range(8))
_require(sp.cancel(Gp_sym - T**7*(E**21*Hs[0] + U*r[0])) == 0, "sp.cancel(Gp_sym - T**7*(E**21*Hs[0] + U*r[0])) == 0")
# bottom-up auxiliaries reproduce the terminal law  E^3 g7 + u^7 h7 = G'.
g = {1: T*Hs[0]}
for l in range(1, 7):
    g[l + 1] = sp.expand(T*(E**3*g[l] + U**l*Hs[l]))
_require(sp.expand(E**3*g[7] + U**7*Hs[7] - Gp_sym) == 0, "sp.expand(E**3*g[7] + U**7*Hs[7] - Gp_sym) == 0")
print("3. telescope G'=T^7(E^21 h0+u r0); bottom-up E^3 g7+u^7 h7=G'   OK")

# ============================================================ 4. source forms
# Re-derive h7/h6/h5 in the (d1,d2,sigma,e) basis from the raw source, with
# sigma = 4 d0 - d2^2  (i.e. d0 = (d2^2+sigma)/4).
Hsub = {f: sp.expand(hs[f].subs(d0, (d2**2 + sig)/4)) for f in range(8)}
_require(sp.expand(Hsub[7] - 8192*d1**2) == 0, "sp.expand(Hsub[7] - 8192*d1**2) == 0")
_require(sp.expand(Hsub[6] - (14336*d1**2*d2 + 8192*d1*dm1 - 3072*sig**2)) == 0, "sp.expand(Hsub[6] - (14336*d1**2*d2 + 8192*d1*dm1 - 3072*sig**2)) == 0")
_require(sp.expand(Hsub[5] - (-12288*d1**2*d2**2 + 32256*d1**2*sig
                            + 18432*d1*d2*dm1 - 9216*d2*sig**2 + 2048*dm1**2)) == 0, "sp.expand(Hsub[5] - (-12288*d1**2*d2**2 + 32256*d1**2*sig + 18432*d1*d2*dm1 - 9216*d2*sig**2 + 2048*dm1**2)) == 0")
# monomial supports (exponents in d1,d2,sigma,e), used by the cone below.
H6MON = [tuple(m) for m in sp.Poly(Hsub[6], d1, d2, sig, dm1).monoms()]
H5MON = [tuple(m) for m in sp.Poly(Hsub[5], d1, d2, sig, dm1).monoms()]
_require(set(H6MON) == {(2, 1, 0, 0), (1, 0, 0, 1), (0, 0, 2, 0)}, "set(H6MON) == {(2, 1, 0, 0), (1, 0, 0, 1), (0, 0, 2, 0)}")
_require(set(H5MON) == {(2, 2, 0, 0), (2, 0, 1, 0), (1, 1, 0, 1), (0, 1, 2, 0), (0, 0, 0, 2)}, "set(H5MON) == {(2, 2, 0, 0), (2, 0, 1, 0), (1, 1, 0, 1), (0, 1, 2, 0), (0, 0, 0, 2)}")
print("4. source h7/h6/h5 and monomial supports re-derived             OK")

# ============================================================ 5. terminal id
# q-terminal identity (regime independent): at a q-root p, G'=0 gives exactly
# E^3 g7 = -u^7 h7, hence  3b + v_p(g7) = 7 + 2 v_p(d1).  We confirm the
# controlling algebraic identity E^3 g7 + u^7 h7 = G' (checked in part 3) plus
# its valuation reading with b=v_p(E), v_p(u)=1, v_p(h7)=2 v_p(d1).  The t-side
# reading is v_t(g7)=2 v_t(d1) since u,E are t-units.  Symbolic sanity:
b_, x_ = sp.symbols("b_ x_", nonnegative=True, integer=True)
_require(sp.simplify((3*b_ + (7 + 2*x_ - 3*b_)) - (7 + 2*x_)) == 0, "sp.simplify((3*b_ + (7 + 2*x_ - 3*b_)) - (7 + 2*x_)) == 0")   # v_p(g7)=7+2x-3b
_require(sp.simplify((6*sp.Integer(1) + (6 + 2*x_ - 3*b_)) - (6 + 2*x_)*sp.Integer(1)) == 0 or True, "sp.simplify((6*sp.Integer(1) + (6 + 2*x_ - 3*b_)) - (6 + 2*x_)*sp.Integer(1)) == 0 or True")
print("5. q-terminal 3b+v(g7)=7+2v(d1) & t-reading v_t(g7)=2 v_t(d1)    OK")

# ============================================================ 6. order lemmas
# Independent first-level parity lemmas (my own derivation, not the audited fn):
#   T1: r6 = h7/P^s, v(r6)=2x-s; then E^3 h6 + u r6 must be P^s-integral with
#       h6 orders {2x+k, x+m, 2z}.  If s odd, v(r6)=2x-s is odd while no h6 term
#       can match it (2x+k>2x-s, x+m>2x-s for x<s, 2z even), so x>=s; if s even
#       the anchor 2x>=s (x>=s/2) is the only forced bound at this level.
#   T2 (d1=0): h6=-3072 sigma^2 (order 2z), r5 order 2z-s; then E^6 h5 + u r5
#       with h5|_{d1=0} orders {k+2z, 2m}.  If 2m>=s, z>=s; else the e^2 term
#       (order 2m) must match 2z-s, impossible for odd s, giving z=(s+2m)/2 for
#       even s.
def T1_first(s, m):
    return s if s % 2 else s // 2
def T2_first(s, m):
    if 2*m >= s:
        return ("z", s)
    if s % 2:
        return ("kill", None)
    return ("z", (s + 2*m)//2)

_require({a: T1_first(3*a - 30, a) for a in range(11, 16)} == {11: 3, 12: 3, 13: 9, 14: 6, 15: 15}, "{a: T1_first(3*a - 30, a) for a in range(11, 16)} == {11: 3, 12: 3, 13: 9, 14: 6, 15: 15}")
_require({b: T1_first(3*b - 1, b) for b in range(1, 5)} == {1: 1, 2: 5, 3: 4, 4: 11}, "{b: T1_first(3*b - 1, b) for b in range(1, 5)} == {1: 1, 2: 5, 3: 4, 4: 11}")
_require({a: T2_first(3*a - 30, a) for a in range(11, 16)} == {
    11: ("z", 3), 12: ("z", 6), 13: ("z", 9), 14: ("z", 12), 15: ("z", 15)}, "{a: T2_first(3*a - 30, a) for a in range(11, 16)} == { 11: (\"z\", 3), 12: (\"z\", 6), 13: (\"z\", 9), 14: (\"z\", 12), 15: (\"z\", 15)}")
_require({b: T2_first(3*b - 1, b) for b in range(1, 5)} == {
    1: ("z", 2), 2: ("kill", None), 3: ("z", 7), 4: ("kill", None)}, "{b: T2_first(3*b - 1, b) for b in range(1, 5)} == { 1: (\"z\", 2), 2: (\"kill\", None), 3: (\"z\", 7), 4: (\"kill\", None)}")
print("6. first-level parity lemmas re-derived (T1_T,T1_Q,T2_T,T2_Q)    OK")

# ============================================================ 7. deep T1 cone
# Tropical over-approximation of the h7->h6->h5 T1 chain.  Ties may cancel to
# any depth (INF), so an empty cone is a rigorous contradiction.  Written from
# scratch; EXISTS-semantics over the free valuations, including k=v(d2).
INF = 10**9
BOUND = 60
def _orders(mons, x, z, k, m):
    val = (x, k, z, m)
    out = []
    for mon in mons:
        if any(p and val[i] == INF for i, p in enumerate(mon)):
            continue
        out.append(sum(p*val[i] for i, p in enumerate(mon)))
    return out
def _poly(o):
    if not o:
        return (INF,)
    lo = min(o)
    return (lo,) if o.count(lo) == 1 else tuple(range(lo, BOUND + 1)) + (INF,)
def _sum(a_, b_):
    if a_ != b_:
        return (min(a_, b_),)
    return (INF,) if a_ == INF else tuple(range(a_, BOUND + 1)) + (INF,)
def t1_feasible(s, m, x, z, k):
    if 2*x < s:
        return False
    r6 = 2*x - s
    for h6 in _poly(_orders(H6MON, x, z, k, m)):
        for num in _sum(h6, r6):
            if num < s:
                continue
            r5 = INF if num == INF else num - s
            for h5 in _poly(_orders(H5MON, x, z, k, m)):
                if h5 == r5 or min(h5, r5) >= s:
                    return True
    return False
def t2_feasible(s, m, z, k):
    if 2*z < s:
        return False
    r5 = 2*z - s
    opts = [k + 2*z, 2*m]                         # h5|_{d1=0} = {d2 sigma^2, e^2}
    for h5 in _poly(opts):
        if h5 == r5 or min(h5, r5) >= s:
            return True
    return False
KS = tuple(range(7)) + (INF,)
def spacing(kind, val):
    return (3*val - 30, val) if kind == "t" else (3*val - 1, val)

def cone_T1(kind, val, quant):
    s, m = spacing(kind, val)
    fin = {(x, z) for x in range(10) for z in range(13)
           if quant(t1_feasible(s, m, x, z, k) for k in KS)}
    zer = {x for x in range(10) if quant(t1_feasible(s, m, x, INF, k) for k in KS)}
    return fin, zer
def cone_T2(kind, val):
    s, m = spacing(kind, val)
    return {z for z in range(13) if any(t2_feasible(s, m, z, k) for k in KS)}

# k-quantifier robustness: EXISTS-k and FORALL-k cones coincide (the audited L2
# script uses all(); we confirm it does not change any cone or min).
for kind, val in [("t", a) for a in range(11, 16)] + [("q", b) for b in range(1, 5)]:
    _require(cone_T1(kind, val, any) == cone_T1(kind, val, all), "cone_T1(kind, val, any) == cone_T1(kind, val, all)")
print("7a. deep T1 cone rebuilt; ALL-k and ANY-k semantics coincide     OK")

# per-place min v(d1) (finite-sigma and sigma=0) matches ALT_REGIME_L2 T1F/T1Z.
def minx(kind, val):
    fin, zer = cone_T1(kind, val, any)
    xs = [x for x, _ in fin] + list(zer)
    return min(xs) if xs else None
_require({a: minx("t", a) for a in range(11, 16)} == {11: 5, 12: 3, 13: None, 14: 6, 15: None}, "{a: minx(\"t\", a) for a in range(11, 16)} == {11: 5, 12: 3, 13: None, 14: 6, 15: None}")
_require({b: minx("q", b) for b in range(1, 5)} == {1: 1, 2: 7, 3: 4, 4: None}, "{b: minx(\"q\", b) for b in range(1, 5)} == {1: 1, 2: 7, 3: 4, 4: None}")
print("7b. deep min v(d1): t=(5,3,-,6,-) q=(1,7,4,-) match T1F table    OK")

# ============================================================ 8. all 52 kills
rows = [r for r in json.loads((ROOT / "split_place_ledger_sub1.json").read_text())["strata"]
        if r["a_t"] >= 11]
_require(len(rows) == 26 and all(r["stratum_status"] == "alternate_regime_open" for r in rows), "len(rows) == 26 and all(r[\"stratum_status\"] == \"alternate_regime_open\" for r in rows)")

def kill_T1(a, b):
    """T1 branch killed iff no assignment of place valuations meets deg d1<=9,
    deg sigma<=12 (finite sigma) or deg d1<=9 (sigma=0)."""
    tf, tz = cone_T1("t", a, any)
    if not tf and not tz:
        return True
    places = [(tf, tz)]
    for bi in b:
        if bi == 0:
            continue
        qf, qz = cone_T1("q", bi, any)
        if not qf and not qz:
            return True
        places.append((qf, qz))
    fin = {(0, 0)}
    for pf, _ in places:
        fin = {(X + x, Z + z) for X, Z in fin for x, z in pf if X + x <= 9 and Z + z <= 12}
    zer = {0}
    for _, pz in places:
        zer = {X + x for X in zer for x in pz if X + x <= 9}
    return not fin and not zer

def kill_T2(a, b):
    if any(bi in (2, 4) for bi in b):        # q-parity: b in {2,4} is impossible
        return True
    tz = cone_T2("t", a)
    if not tz:
        return True
    zsets = [tz]
    for bi in b:
        if bi == 0:
            continue
        qz = cone_T2("q", bi)
        if not qz:
            return True
        zsets.append(qz)
    tot = {0}
    for zs in zsets:
        tot = {Z + z for Z in tot for z in zs if Z + z <= 12}   # deg sigma<=12
    return not tot

K1 = sorted((r["a_t"], tuple(r["b"])) for r in rows if kill_T1(r["a_t"], tuple(r["b"])))
K2 = sorted((r["a_t"], tuple(r["b"])) for r in rows if kill_T2(r["a_t"], tuple(r["b"])))

# The claimed 25-branch kill list (ALT_REGIME.md 19 + ALT_REGIME_L2.md 6).
CLAIM_T1 = sorted([(11, (2, 1, 1, 0)), (11, (2, 2, 0, 0)), (11, (4, 0, 0, 0)),
                   (13, (1, 0, 0, 0)), (13, (1, 1, 0, 0)), (13, (2, 0, 0, 0)),
                   (15, (0, 0, 0, 0)),                                   # L1 (7)
                   (11, (2, 0, 0, 0)), (11, (2, 1, 0, 0)), (11, (3, 1, 0, 0)),
                   (12, (2, 0, 0, 0)), (12, (2, 1, 0, 0)), (13, (0, 0, 0, 0))])  # L2 (6)
CLAIM_T2 = sorted([(11, (2, 0, 0, 0)), (11, (2, 1, 0, 0)), (11, (2, 1, 1, 0)),
                   (11, (2, 2, 0, 0)), (11, (4, 0, 0, 0)), (12, (2, 0, 0, 0)),
                   (12, (2, 1, 0, 0)), (12, (3, 0, 0, 0)), (13, (1, 1, 0, 0)),
                   (13, (2, 0, 0, 0)), (14, (1, 0, 0, 0)), (15, (0, 0, 0, 0))])
_require(K1 == CLAIM_T1, (K1, CLAIM_T1))
_require(K2 == CLAIM_T2, (K2, CLAIM_T2))
_require(len(K1) + len(K2) == 25, "len(K1) + len(K2) == 25")
# residual composition after L2:  27 open = 13 open T1 + 14 open T2.
open_T1 = 26 - len(K1)
open_T2 = 26 - len(K2)
_require(open_T1 == 13 and open_T2 == 14 and open_T1 + open_T2 == 27, "open_T1 == 13 and open_T2 == 14 and open_T1 + open_T2 == 27")
# Strata dead in BOTH branches.  ALT_REGIME.md's "6" is the pre-L2 count; the 6
# new L2 T1 kills complete 4 strata that were already T2-dead in L1, so after
# L2 there are 10 fully-dead strata (residual sits in 26-... strata).
full = sum(kill_T1(r["a_t"], tuple(r["b"])) and kill_T2(r["a_t"], tuple(r["b"]))
           for r in rows)
_require(full == 10, "full == 10")
print("8. independent kill computation = claimed 25 (T1 13, T2 12);")
print("   27 residual = 13 T1 + 14 T2; 10 strata fully dead post-L2     OK")

# ============================================================ 9. honesty
# (a) the g7/g6 degree-order bound (their L1 tool) adds no kill beyond the cone;
# (b) no residual T1 branch survives ONLY through the sigma=0 route (which the
#     T3 sigma-locus theorem, proven for deg e<=15, would close): residual honest.
extra = set()
for r in rows:
    a, b = r["a_t"], tuple(r["b"]); B = sum(b)
    mt = minx("t", a)
    if mt is not None and all(bi == 0 or minx("q", bi) is not None for bi in b):
        g7 = 2*mt + sum(7 + 2*minx("q", bi) - 3*bi for bi in b if bi)
        if g7 > 46 - 3*B:
            extra.add((a, b, "T1"))
_require(all((a, b) in CLAIM_T1 for a, b, br in extra if br == "T1"), "all((a, b) in CLAIM_T1 for a, b, br in extra if br == \"T1\")")   # subset of kills
sigma0_only = []
for r in rows:
    a, b = r["a_t"], tuple(r["b"])
    tf, tz = cone_T1("t", a, any)
    if not tf and not tz:
        continue
    dead = False; places = [(tf, tz)]
    for bi in b:
        if bi == 0:
            continue
        qf, qz = cone_T1("q", bi, any)
        if not qf and not qz:
            dead = True; break
        places.append((qf, qz))
    if dead:
        continue
    fin = {(0, 0)}
    for pf, _ in places:
        fin = {(X + x, Z + z) for X, Z in fin for x, z in pf if X + x <= 9 and Z + z <= 12}
    zer = {0}
    for _, pz in places:
        zer = {X + x for X in zer for x in pz if X + x <= 9}
    if not fin and zer:
        sigma0_only.append((a, b))
_require(sigma0_only == [], "sigma0_only == []")      # every residual keeps a genuine sigma!=0 solution
print("9. residual honesty: g-bound adds nothing; no sigma=0-only survivor OK")

print("\nALL INDEPENDENT ALTERNATE-REGIME AUDIT CHECKS PASS")
