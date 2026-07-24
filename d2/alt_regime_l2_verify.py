#!/usr/bin/env python3
"""Exact source-linked certificate for ALT_REGIME_L2.md.

Tied local minima are allowed arbitrary rise, so the search is a deliberate
over-approximation. An empty cone is therefore a rigorous contradiction.
"""
from __future__ import annotations
import json, re
from pathlib import Path
import sympy as sp

def _require(_cond, _msg):
    """Proof-critical check: fails loudly and exits nonzero, unaffected by python -O."""
    if not _cond:
        import sys as _sys
        print("FAIL: " + str(_msg))
        _sys.exit(1)


ROOT = Path(__file__).resolve().parent
d2, d1, d0, e, sigma = sp.symbols("d2 d1 d0 dm1 sigma")
source = (ROOT / "f31_graded.txt").read_text(encoding="utf-8")
pat = r"h_(\d) \(weight \d+, dm1-power \d+\) = (.+)"
hs = {int(m.group(1)): sp.sympify(m.group(2)) for m in re.finditer(pat, source)}
_require(sorted(hs) == list(range(8)), "sorted(hs) == list(range(8))")
H = {f: sp.expand(hs[f].subs(d0, (d2**2 + sigma)/4)) for f in range(8)}

# A. Source formulas.
_require(sp.expand(H[7] - 8192*d1**2) == 0, "sp.expand(H[7] - 8192*d1**2) == 0")
_require(sp.expand(H[6] - (14336*d1**2*d2 + 8192*d1*e - 3072*sigma**2)) == 0, "sp.expand(H[6] - (14336*d1**2*d2 + 8192*d1*e - 3072*sigma**2)) == 0")
_require(sp.expand(H[5] - (-12288*d1**2*d2**2 + 32256*d1**2*sigma
                         + 18432*d1*d2*e - 9216*d2*sigma**2 + 2048*e**2)) == 0, "sp.expand(H[5] - (-12288*d1**2*d2**2 + 32256*d1**2*sigma + 18432*d1*d2*e - 9216*d2*sigma**2 + 2048*e**2)) == 0")
print("A. source-linked h7/h6/h5 formulas                                OK")

# B. Exact flipped and q-unit transitions.
A, U, E = sp.symbols("A U E", nonzero=True)
r = list(sp.symbols("r0:7")); hh = {7: A*r[6]}
for f in range(6, 0, -1):
    hh[f] = (A*r[f-1] - U*r[f])/E**(21-3*f)
hh[0] = -U*r[0]/E**21
_require(sp.cancel(sum(A**(7-f)*U**f*E**(21-3*f)*hh[f] for f in range(8))) == 0, "sp.cancel(sum(A**(7-f)*U**f*E**(21-3*f)*hh[f] for f in range(8))) == 0")
for b in range(1, 5):
    s = 3*b-1
    _require(all(f+b*(21-3*f) == 7+(7-f)*s for f in range(8)), "all(f+b*(21-3*f) == 7+(7-f)*s for f in range(8))")
T = sp.symbols("T", nonzero=True); X = list(sp.symbols("X0:8")); g = {1: T*X[0]}
for level in range(1, 7):
    g[level+1] = sp.expand(T*(E**3*g[level] + U**level*X[level]))
_require(sp.expand(E**3*g[7] + U**7*X[7]
                 - sum(T**(7-f)*U**f*E**(21-3*f)*X[f] for f in range(8))) == 0, "sp.expand(E**3*g[7] + U**7*X[7] - sum(T**(7-f)*U**f*E**(21-3*f)*X[f] for f in range(8))) == 0")
print("B. exact flipped telescope and t/q transition normalizations       OK")

# C. Parse and tropicalize h6/h5.
VARS = (d1, d2, sigma, e)
MONS = {f: list(sp.Poly(H[f], *VARS).terms()) for f in (6, 5)}
_require([m for m, _ in MONS[6]] == [(2,1,0,0), (1,0,0,1), (0,0,2,0)], "[m for m, _ in MONS[6]] == [(2,1,0,0), (1,0,0,1), (0,0,2,0)]")
_require({m for m, _ in MONS[5]} == {
    (2,2,0,0), (2,0,1,0), (1,1,0,1), (0,1,2,0), (0,0,0,2)}, "{m for m, _ in MONS[5]} == { (2,2,0,0), (2,0,1,0), (1,1,0,1), (0,1,2,0), (0,0,0,2)}")
INF, BOUND = 1000, 100

def monorders(level, x, z, k, m):
    vals = (x, k, z, m); out = []
    for mon, _ in MONS[level]:
        if any(p and v == INF for p, v in zip(mon, vals)): continue
        out.append(sum(p*v for p, v in zip(mon, vals)))
    return out

def polyopts(orders):
    if not orders: return (INF,)
    least = min(orders)
    return ((least,) if orders.count(least) == 1
            else tuple(range(least, BOUND+1)) + (INF,))

def sumopts(left, right):
    if left != right: return (min(left, right),)
    if left == INF: return (INF,)
    return tuple(range(left, BOUND+1)) + (INF,)

def t1_local(s, m, x, z, k):
    """Safe over-approximation to the h7 -> h6 -> h5 local cone."""
    if 2*x < s: return False
    r6 = 2*x-s
    for h6 in polyopts(monorders(6, x, z, k, m)):
        for num in sumopts(h6, r6):
            if num < s: continue
            r5 = INF if num == INF else num-s
            for h5 in polyopts(monorders(5, x, z, k, m)):
                if h5 == r5 or min(h5, r5) >= s: return True
    return False

def t2_local(s, m, z, k):
    if 2*z < s: return False
    r5 = 2*z-s
    return any(h5 == r5 or min(h5, r5) >= s
               for h5 in polyopts(monorders(5, INF, z, k, m)))

PLACES = {"t11":(3,11), "t12":(6,12), "t13":(9,13), "t14":(12,14),
          "q1":(2,1), "q2":(5,2), "q3":(8,3), "q4":(11,4)}
KS = tuple(range(7)) + (INF,)
def project(place, zero=False):
    s, m = PLACES[place]
    if zero:
        return {x for x in range(10) if any(t1_local(s,m,x,INF,k) for k in KS)}
    return {(x,z) for x in range(10) for z in range(13)
            if any(t1_local(s,m,x,z,k) for k in KS)}

T1F = {p: project(p) for p in PLACES}; T1Z = {p: project(p, True) for p in PLACES}
_require(T1F == {
 "t11":{(x,z) for x in range(5,10) for z in range(3,13)},
 "t12":{(x,x-3) for x in range(3,9)}|{(9,z) for z in range(6,13)},
 "t13":set(), "t14":{(x,x-6) for x in range(6,10)},
 "q1":{(1,0),(2,1)}|{(x,z) for x in range(3,10) for z in range(2,13)},
 "q2":{(7,z) for z in range(5,13)},
 "q3":{(x,x-4) for x in range(4,10)}, "q4":set()}, "T1F == { \"t11\":{(x,z) for x in range(5,10) for z in range(3,13)}, \"t12\":{(x,x-3) for x in range(3,9)}|{(9,z) for z in range(6,13)}, \"t13\":set(), \"t14\":{(x,x-6) for x in range(6,10)}, \"q1\":{(1,0),(2,1)}|{(x,z) for x in range(3,10) for z in range(2,13)}, \"q2\":{(7,z) for z in range(5,13)}, \"q3\":{(x,x-4) for x in range(4,10)}, \"q4\":set()}")
_require(T1Z == {"t11":set(range(5,10)), "t12":{9}, "t13":set(), "t14":set(),
                "q1":set(range(3,10)), "q2":{7}, "q3":set(), "q4":set()}, "T1Z == {\"t11\":set(range(5,10)), \"t12\":{9}, \"t13\":set(), \"t14\":set(), \"q1\":set(range(3,10)), \"q2\":{7}, \"q3\":set(), \"q4\":set()}")
T2 = {p:{z for z in range(13) if any(t2_local(s,m,z,k) for k in KS)}
      for p,(s,m) in PLACES.items()}
_require(T2 == {"t11":set(range(3,13)), "t12":set(range(6,13)),
              "t13":set(range(9,13)), "t14":{12}, "q1":set(range(2,13)),
              "q2":set(), "q3":{7}, "q4":set()}, "T2 == {\"t11\":set(range(3,13)), \"t12\":set(range(6,13)), \"t13\":set(range(9,13)), \"t14\":{12}, \"q1\":set(range(2,13)), \"q2\":set(), \"q3\":{7}, \"q4\":set()}")
print("C. source-derived h6/h5 local order cones (ties over-allowed)       OK")

# D. Reproduce ALT_REGIME.md's 33-branch frontier, then audit it.
rows = [r for r in json.loads((ROOT/"split_place_ledger_sub1.json").read_text())["strata"]
        if r["a_t"] >= 11]
_require(len(rows) == 26, "len(rows) == 26")
T1T={11:3,12:3,13:9,14:6,15:15}; T1Q={0:0,1:1,2:5,3:4,4:11}
T2T={11:3,12:6,13:9,14:12,15:15}; T2Q={0:0,1:2,2:3,3:7,4:6}
def oldkill(a,b,branch):
    B=sum(b)
    if branch=="T1":
        tord,qord,cap,lev,raw,local=T1T[a],tuple(T1Q[x] for x in b),9,7,46,False
    else:
        tord,qord,cap,lev,raw,local=T2T[a],tuple(T2Q[x] for x in b),12,6,48,any(x in(2,4) for x in b)
    gos=tuple(lev+2*x-3*bi for x,bi in zip(qord,b)); _require(all(x>=0 for x in gos), "all(x>=0 for x in gos)")
    return local or tord+sum(qord)>cap or 2*tord+sum(gos)>raw-3*B
frontier=[]
for row in rows:
    a,b=row["a_t"],tuple(row["b"])
    frontier += [(a,b,br) for br in ("T1","T2") if not oldkill(a,b,br)]
_require(len(frontier)==33, "len(frontier)==33")

def global_states(a,b):
    keys=[f"t{a}"]+[f"q{x}" for x in b if x]; finite={(0,0)}
    for key in keys:
        finite={(X+x,Z+z) for X,Z in finite for x,z in T1F[key]
                if X+x<=9 and Z+z<=12}
    zero={0}
    for key in keys: zero={X+x for X in zero for x in T1Z[key] if X+x<=9}
    return finite,zero
newkills=set(); survivors=set()
for cell in frontier:
    a,b,br=cell
    if br=="T1" and not any(global_states(a,b)): newkills.add(cell)
    else: survivors.add(cell)
expected={(11,(2,0,0,0),"T1"),(11,(2,1,0,0),"T1"),(11,(3,1,0,0),"T1"),
          (12,(2,0,0,0),"T1"),(12,(2,1,0,0),"T1"),(13,(0,0,0,0),"T1")}
_require(newkills==expected and len(survivors)==27, "newkills==expected and len(survivors)==27")
print("D. exhaustive 33-branch audit: 6 new kills, 27 residual branches   OK")

# E. UFD residual normal form F=s*u^2, D=s^2*u^3*v, G=s*v^2.
aa,uu,vv=sp.symbols("aa uu vv", integer=True, nonnegative=True)
_require(sp.expand(3*(aa+2*uu)+(aa+2*vv)-2*(2*aa+3*uu+vv))==0, "sp.expand(3*(aa+2*uu)+(aa+2*vv)-2*(2*aa+3*uu+vv))==0")
for a,b,br in survivors:
    _require(15-a-sum(b)>=0, "15-a-sum(b)>=0")
    if br=="T1": _require(any(global_states(a,b)), "any(global_states(a,b))")
    else: _require(min(T2[f"t{a}"])+sum(min(T2[f"q{x}"]) for x in b if x)<=12, "min(T2[f\"t{a}\"])+sum(min(T2[f\"q{x}\"]) for x in b if x)<=12")
print("E. terminal UFD residual normal form and survivor caps             OK")
print("\nALL ALTERNATE-REGIME LEVEL-2/3 CHECKS PASS")
