#!/usr/bin/env python3
"""Exact checks for the f31/subcase-(1) alternate regime a=11,...,15."""
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
d2, d1, d0, dm1, Phi, y = sp.symbols("d2 d1 d0 dm1 Phi y")
t = y + 1
q = 2048*y**4 - 512*y**3 + 320*y**2 - 240*y + 195
c = sp.Rational(-1, 6630)
u, phi_tilde = c*q, c*t**30*q
text = (ROOT / "f31_graded.txt").read_text(encoding="utf-8")
pattern = r"h_(\d) \(weight \d+, dm1-power \d+\) = (.+)"
hs = {int(m.group(1)): sp.sympify(m.group(2)) for m in re.finditer(pattern, text)}
_require(sorted(hs) == list(range(8)), "sorted(hs) == list(range(8))")

# 1. Orders and constant common factor.
for a in range(11, 16):
    v, w = 30-3*a, 3*a-30
    orders = [30*f+a*(21-3*f) for f in range(8)]
    _require(v < 0 and w == -v and 21*a+7*v == 210, "v < 0 and w == -v and 21*a+7*v == 210")
    _require(orders[7] == min(orders) == 210, "orders[7] == min(orders) == 210")
    _require(all(orders[f]-210 == (7-f)*w for f in range(8)), "all(orders[f]-210 == (7-f)*w for f in range(8))")
print("1. alternate orders: min at f=7 and 21a+7v=210              OK")

# 2. Random-window mirror of t5_multiplace_verify.py check 6 (a=12).
random.seed(1207)
def rpoly(degree):
    return sp.expand(y**degree + sum(random.randint(-2, 2)*y**j for j in range(degree+1)))
a, w = 12, 6
D2, D1, D0, Ehat = rpoly(2), rpoly(3), rpoly(4), rpoly(2)
if Ehat.subs(y, -1) == 0: Ehat += 1
_require(Ehat.subs(y, -1) != 0 and sp.degree(Ehat, y) <= 15-a, "Ehat.subs(y, -1) != 0 and sp.degree(Ehat, y) <= 15-a")
Efull = t**a*Ehat
hval = {f: hs[f].subs({d2:D2, d1:D1, d0:D0, dm1:Efull}) for f in range(8)}
_require(sp.cancel(phi_tilde/(t**30*u)) == 1, "sp.cancel(phi_tilde/(t**30*u)) == 1")
for f in range(8):
    _require(30*f+a*(21-3*f) == 210+(7-f)*w, "30*f+a*(21-3*f) == 210+(7-f)*w")
    ratio = (phi_tilde**f*Efull**(21-3*f)
             / (t**(210+(7-f)*w)*u**f*Ehat**(21-3*f)))
    _require(sp.cancel(ratio) == 1, "sp.cancel(ratio) == 1")
    left = phi_tilde**f*Efull**(21-3*f)*hval[f]
    right = t**(210+(7-f)*w)*u**f*Ehat**(21-3*f)*hval[f]
    for sample in (-2, 0, 1):
        _require(sp.cancel((left-right).subs(y, sample)) == 0, "sp.cancel((left-right).subs(y, sample)) == 0")
print("2. random a=12 window: F=t^210 G' exactly                       OK")

# 3. Exact descending telescope.
T, U, E = sp.symbols("T U E", nonzero=True)
r = list(sp.symbols("r0:7")); hh = {7:T*r[6]}
for f in range(6, 0, -1): hh[f] = (T*r[f-1]-U*r[f])/E**(21-3*f)
hh[0] = -U*r[0]/E**21
_require(sp.cancel(sum(T**(7-f)*U**f*E**(21-3*f)*hh[f] for f in range(8))) == 0, "sp.cancel(sum(T**(7-f)*U**f*E**(21-3*f)*hh[f] for f in range(8))) == 0")
print("3. descending flipped recursion telescopes to G'=0                 OK")

# 4. Bottom-up polynomial g_l recover the old terminal identity.
H = list(sp.symbols("H0:8")); g = {1:T*H[0]}
for level in range(1, 7): g[level+1] = sp.expand(T*(E**3*g[level]+U**level*H[level]))
terminal = sp.expand(E**3*g[7]+U**7*H[7])
generic = sum(T**(7-f)*U**f*E**(21-3*f)*H[f] for f in range(8))
_require(sp.expand(terminal-generic) == 0, "sp.expand(terminal-generic) == 0")
print("4. polynomial bottom-up g_l gives the old q-terminal identity       OK")

# 5. Source identities for local first-level lemmas.
sigma = 4*d0-d2**2
_require(sp.expand(hs[7]-8192*d1**2) == 0, "sp.expand(hs[7]-8192*d1**2) == 0")
_require(sp.expand(hs[6]-(-3072*sigma**2+14336*d1**2*d2+8192*d1*dm1)) == 0, "sp.expand(hs[6]-(-3072*sigma**2+14336*d1**2*d2+8192*d1*dm1)) == 0")
_require(sp.expand(hs[5].subs(d1,0)-(-9216*d2*sigma**2+2048*dm1**2)) == 0, "sp.expand(hs[5].subs(d1,0)-(-9216*d2*sigma**2+2048*dm1**2)) == 0")
print("5. h7, h6, and h5|d1=0 local source identities                    OK")

def t1_order(s, m):
    anchor = (s+1)//2
    if s % 2:
        for x in range(anchor, s):
            k = 2*x-s
            _require(0 <= k < s and k%2 and 2*x > k and x+m > k, "0 <= k < s and k%2 and 2*x > k and x+m > k")
            _require(all(2*z != k for z in range(30)), "all(2*z != k for z in range(30))")
        return s
    return anchor

def t2_order(s, m):
    if 2*m >= s: return ("order", s)
    if s % 2:
        _require((2*m-s)%2, "(2*m-s)%2")
        return ("killed", None)
    z = (s+2*m)//2
    _require(2*z-s == 2*m < s, "2*z-s == 2*m < s")
    return ("order", z)

_require({a:t1_order(3*a-30,a) for a in range(11,16)} == {11:3,12:3,13:9,14:6,15:15}, "{a:t1_order(3*a-30,a) for a in range(11,16)} == {11:3,12:3,13:9,14:6,15:15}")
_require({b:t1_order(3*b-1,b) for b in range(1,5)} == {1:1,2:5,3:4,4:11}, "{b:t1_order(3*b-1,b) for b in range(1,5)} == {1:1,2:5,3:4,4:11}")
_require({a:t2_order(3*a-30,a) for a in range(11,16)} == {11:("order",3),12:("order",6),13:("order",9),14:("order",12),15:("order",15)}, "{a:t2_order(3*a-30,a) for a in range(11,16)} == {11:(\"order\",3),12:(\"order\",6),13:(\"order\",9),14:(\"order\",12),15:(\"order\",15)}")
_require({b:t2_order(3*b-1,b) for b in range(1,5)} == {1:("order",2),2:("killed",None),3:("order",7),4:("killed",None)}, "{b:t2_order(3*b-1,b) for b in range(1,5)} == {1:(\"order\",2),2:(\"killed\",None),3:(\"order\",7),4:(\"killed\",None)}")
print("6. terminal + first-level local order/parity lemmas                  OK")

# 7. Exhaustive audit of all 26 alternate strata / 52 branches.
payload = json.loads((ROOT/"split_place_ledger_sub1.json").read_text(encoding="utf-8"))
rows = [row for row in payload["strata"] if row["a_t"] >= 11]
_require(len(rows) == 26 and all(row["stratum_status"] == "alternate_regime_open" for row in rows), "len(rows) == 26 and all(row[\"stratum_status\"] == \"alternate_regime_open\" for row in rows)")
T1_T={11:3,12:3,13:9,14:6,15:15}; T1_Q={0:0,1:1,2:5,3:4,4:11}
T2_T={11:3,12:6,13:9,14:12,15:15}; T2_Q={0:0,1:2,2:3,3:7,4:6}

def audit(a, b, branch):
    B=sum(b); _require(B <= 15-a <= 4, "B <= 15-a <= 4")
    if branch=="T1":
        level,ac,raw,tord,qord,local=7,9,46,T1_T[a],tuple(T1_Q[x] for x in b),False
    else:
        level,ac,raw,tord=6,12,48,T2_T[a]
        local=any(x in (2,4) for x in b); qord=tuple(T2_Q[x] for x in b)
    al=tord+sum(qord)
    qg=tuple(level+2*x-3*bi for x,bi in zip(qord,b)); _require(all(x>=0 for x in qg), "all(x>=0 for x in qg)")
    gl=2*tord+sum(qg); gc=raw-3*B
    return {"killed":local or al>ac or gl>gc,"aux":al,"g":gl,"gcap":gc}

V={}
for row in rows:
    a,b=row["a_t"],tuple(row["b"]); _require(a+sum(b)<=15, "a+sum(b)<=15")
    for branch in ("T1","T2"): V[(a,b,branch)]=audit(a,b,branch)
K1=sorted((a,b) for (a,b,br),z in V.items() if br=="T1" and z["killed"])
K2=sorted((a,b) for (a,b,br),z in V.items() if br=="T2" and z["killed"])
_require(K1==[(11,(2,1,1,0)),(11,(2,2,0,0)),(11,(4,0,0,0)),(13,(1,0,0,0)),(13,(1,1,0,0)),(13,(2,0,0,0)),(15,(0,0,0,0))], "K1==[(11,(2,1,1,0)),(11,(2,2,0,0)),(11,(4,0,0,0)),(13,(1,0,0,0)),(13,(1,1,0,0)),(13,(2,0,0,0)),(15,(0,0,0,0))]")
_require(K2==[(11,(2,0,0,0)),(11,(2,1,0,0)),(11,(2,1,1,0)),(11,(2,2,0,0)),(11,(4,0,0,0)),(12,(2,0,0,0)),(12,(2,1,0,0)),(12,(3,0,0,0)),(13,(1,1,0,0)),(13,(2,0,0,0)),(14,(1,0,0,0)),(15,(0,0,0,0))], "K2==[(11,(2,0,0,0)),(11,(2,1,0,0)),(11,(2,1,1,0)),(11,(2,2,0,0)),(11,(4,0,0,0)),(12,(2,0,0,0)),(12,(2,1,0,0)),(12,(3,0,0,0)),(13,(1,1,0,0)),(13,(2,0,0,0)),(14,(1,0,0,0)),(15,(0,0,0,0))]")
_require(sum(z["killed"] for z in V.values())==19, "sum(z[\"killed\"] for z in V.values())==19")
_require(sum(not z["killed"] for z in V.values())==33, "sum(not z[\"killed\"] for z in V.values())==33")
full=sum(all(V[(row["a_t"],tuple(row["b"]),br)]["killed"] for br in ("T1","T2")) for row in rows)
_require(full==6, "full==6")
print("7. ledger audit: 19/52 branches killed; 33 open in 20 strata       OK")
print("\nALL ALTERNATE-REGIME CHECKS PASS")
