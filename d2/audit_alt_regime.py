#!/usr/bin/env python3
"""Spec-only independent audit of the sub1 alternate regime.

Uses only f31_graded.txt and the two documented data artifacts. It does not
import or read existing alternate-regime verifier/sweep implementations or
cascade_engine.py.
"""
from __future__ import annotations
import argparse, json, random, re, sys, time
from functools import lru_cache
from itertools import product
from pathlib import Path
import sympy as sp

HERE=Path(__file__).resolve().parent
NEG,INF=-10**9,10**9
CAPS=(6,9,12,15)
class AuditError(RuntimeError): pass
def need(ok,msg):
    if not ok: raise AuditError(msg)
def bid(a,b,branch): return f"a{a}_b{''.join(map(str,b))}_{branch}"
def load(name):
    try: return json.loads((HERE/name).read_text(encoding="utf-8-sig"))
    except Exception as exc: raise AuditError(f"cannot read {name}: {exc}") from exc

def parse_graded():
    """Own parser, with original and sigma-basis homogeneity checks."""
    text=(HERE/"f31_graded.txt").read_text(encoding="utf-8-sig")
    need("f31 = sum_{f=0}^7 Phi^f * dm1^(21-3f) * h_f" in text,"graded header missing")
    d0,d1,d2,dm1,sigma,e=sp.symbols("d0 d1 d2 dm1 sigma e")
    rows=re.findall(r"^h_(\d+) \(weight (\d+), dm1-power (\d+)\) = (.+)$",text,re.M)
    need(len(rows)==8,f"expected 8 h-levels, got {len(rows)}")
    exprs,tabs={},{}
    for fs,ws,ps,source in rows:
        f,w,p=int(fs),int(ws),int(ps); need(w==20-2*f and p==21-3*f,f"h_{f} metadata")
        try:
            q=sp.Poly(sp.sympify(source,locals={str(x):x for x in (d0,d1,d2,dm1)}),d0,d1,d2,dm1,domain=sp.QQ)
        except Exception as exc: raise AuditError(f"cannot parse h_{f}: {exc}") from exc
        for mon,c in q.terms(): need(c and sum(x*y for x,y in zip(mon,(4,3,2,5)))==w,f"h_{f} nonhomogeneous {mon}")
        conv=sp.Poly(sp.expand(q.as_expr().subs({d0:(sigma+d2**2)/4,dm1:e})),d2,d1,sigma,e,domain=sp.QQ)
        exps=[]
        for mon,c in conv.terms():
            need(c and sum(x*y for x,y in zip(mon,(2,3,4,5)))==w,f"sigma h_{f} nonhomogeneous {mon}")
            need(sum(x*y for x,y in zip(mon,CAPS))<=60-6*f,f"h_{f} degree cap at {mon}")
            exps.append(tuple(map(int,mon)))
        exprs[f],tabs[f]=conv.as_expr(),tuple(exps)
    need(set(exprs)==set(range(8)),"h indices are not 0..7")
    need(sp.expand(exprs[7]-8192*d1**2)==0,"h7 mismatch")
    h6=14336*d1**2*d2+8192*d1*e-3072*sigma**2
    h5=-12288*d1**2*d2**2+32256*d1**2*sigma+18432*d1*d2*e-9216*d2*sigma**2+2048*e**2
    need(sp.expand(exprs[6]-h6)==0 and sp.expand(exprs[5]-h5)==0,"h6/h5 mismatch")
    return tabs

def audit_algebra():
    for a in range(11,16):
        v,w=30-3*a,3*a-30; orders=[30*f+a*(21-3*f) for f in range(8)]
        need(v<0 and w==-v and orders.count(210)==1 and orders[-1]==min(orders),f"a={a}: flipped minimum")
        need([x-210 for x in orders]==[(7-f)*w for f in range(8)],f"a={a}: residual exponents")
    T,u,E=sp.symbols("T u E",nonzero=True); h=sp.symbols("h0:8"); r=h[7]/T
    for f in range(6,0,-1): r=(E**(3*(7-f))*h[f]+u*r)/T
    G=sum(T**(7-f)*u**f*E**(21-3*f)*h[f] for f in range(8))
    need(sp.cancel(G-T**7*(E**21*h[0]+u*r))==0,"closed telescope failed")
    for b in range(1,5):
        for f in range(8): need(f+b*(21-3*f)==7+(7-f)*(3*b-1),f"q identity b={b},f={f}")

def wanted_strata():
    return {(a,b) for a in range(11,16) for b in product(range(16-a),repeat=4)
            if b==tuple(sorted(b,reverse=True)) and sum(b)<=15-a}

def audit_ledger(data):
    rows=[x for x in data.get("strata",[]) if 11<=x.get("a_t",-1)<=15]
    got={(x["a_t"],tuple(x["b"])) for x in rows}; want=wanted_strata()
    need(len(rows)==len(got) and got==want and len(got)==26,
         f"strata mismatch: missing={sorted(want-got)}, extra={sorted(got-want)}")
    for x in rows:
        a,b=x["a_t"],tuple(x["b"])
        need(x.get("stratum_status")=="alternate_regime_open" and x.get("open_branches")==["T1","T2"],f"ledger {(a,b)} status")
        need(x.get("q_multiplicity_sum")==sum(b) and x.get("residual_degree_budget")==15-a-sum(b),f"ledger {(a,b)} budget")
    s=data.get("summary",{}); need(s.get("alternate_regime_strata")==26 and s.get("alternate_regime_open_branches")==52,"ledger headline")
    return sorted(want)

L1T1={"a11_b2110_T1","a11_b2200_T1","a11_b4000_T1","a13_b1000_T1",
      "a13_b1100_T1","a13_b2000_T1","a15_b0000_T1"}
L1T2={"a11_b2000_T2","a11_b2100_T2","a11_b2110_T2","a11_b2200_T2",
      "a11_b4000_T2","a12_b2000_T2","a12_b2100_T2","a12_b3000_T2",
      "a13_b1100_T2","a13_b2000_T2","a14_b1000_T2","a15_b0000_T2"}
def t1min(s): return s if s%2 else s//2
def t2cond(s,m):
    if 2*m>=s: return True,s
    return (False,0) if s%2 else (True,(s+2*m)//2)
def derive_l1(strata):
    k1,k2=set(),set()
    for a,bs in strata:
        if t1min(3*a-30)+sum(t1min(3*b-1) for b in bs if b)>9: k1.add(bid(a,bs,"T1"))
        ok,z=t2cond(3*a-30,a)
        for b in bs:
            if b and ok: yes,zz=t2cond(3*b-1,b); ok=ok and yes; z+=zz
        if not ok or z>12: k2.add(bid(a,bs,"T2"))
    need(k1==L1T1,f"L1 T1 disagreement: {sorted(k1)}"); need(k2==L1T2,f"L1 T2 disagreement: {sorted(k2)}")
    need(len(k1|k2)==19,"L1 total is not 19"); return k1,k2

def vorders(tab,vals):
    return [INF if any(p and v==INF for p,v in zip(mon,vals)) else sum(p*v for p,v in zip(mon,vals)) for mon in tab]
def divisible(orders,s):
    q=[x for x in orders if x!=INF]; return not q or min(q)>=s or q.count(min(q))>=2
def local_ok(s,m,x,z,k,tabs):
    r6=2*x-s
    if r6<0: return False
    first=vorders(tabs[6],(k,x,z,m))+[r6]; q=[x for x in first if x!=INF]; low=min(q)
    if q.count(low)==1:
        if low<s: return False
        r5s=[low-s]
    else:
        h5=vorders(tabs[5],(k,x,z,m)); ceiling=max([s]+[v for v in h5 if v!=INF])+s+1
        r5s=list(range(max(0,low-s),ceiling+1))+[INF]
    h5=vorders(tabs[5],(k,x,z,m)); return any(divisible(h5+[r],s) for r in r5s)
@lru_cache(None)
def cone(s,m,h6,h5):
    tabs={6:h6,5:h5}; return frozenset((x,z) for x in range(10) for z in list(range(13))+[INF]
        if any(local_ok(s,m,x,z,k,tabs) for k in list(range(7))+[INF]))

def doc_cone(kind,n):
    if (kind,n)==("t",11): return {(x,z) for x in range(5,10) for z in range(3,13)}|{(x,INF) for x in range(5,10)}
    if (kind,n)==("t",12): return {(x,x-3) for x in range(3,9)}|{(9,z) for z in range(6,13)}|{(9,INF)}
    if (kind,n)==("t",13): return set()
    if (kind,n)==("t",14): return {(x,x-6) for x in range(6,10)}
    if (kind,n)==("q",1): return {(1,0),(2,1)}|{(x,z) for x in range(3,10) for z in range(2,13)}|{(x,INF) for x in range(3,10)}
    if (kind,n)==("q",2): return {(7,z) for z in range(5,13)}|{(7,INF)}
    if (kind,n)==("q",3): return {(x,x-4) for x in range(4,10)}
    if (kind,n)==("q",4): return set()
    raise AuditError(f"missing documented cone {kind}={n}")
def combine(sets):
    totals={(0,0)}
    for choices in sets:
        finite={(x,z) for x,z in choices if z!=INF}
        totals={(X+x,Z+z) for X,Z in totals for x,z in finite if X+x<=9 and Z+z<=12}
    if totals: return True
    xs={0}
    for choices in sets:
        zs={x for x,z in choices if z==INF}; xs={X+x for X in xs for x in zs if X+x<=9}
    return bool(xs)
L2={"a11_b2000_T1","a11_b2100_T1","a11_b3100_T1","a12_b2000_T1","a12_b2100_T1","a13_b0000_T1"}
def derive_l2(strata,l1,tabs):
    h6,h5=tabs[6],tabs[5]
    for a in range(11,15): need(set(cone(3*a-30,a,h6,h5))==doc_cone("t",a),f"t cone a={a}")
    for b in range(1,5): need(set(cone(3*b-1,b,h6,h5))==doc_cone("q",b),f"q cone b={b}")
    killed=set()
    for a,bs in strata:
        name=bid(a,bs,"T1")
        if name in l1: continue
        places=[set(cone(3*a-30,a,h6,h5))]+[set(cone(3*b-1,b,h6,h5)) for b in bs if b]
        if not combine(places): killed.add(name)
    need(killed==L2,f"L2 disagreement: {sorted(killed)}"); return killed

def plus(n,d): return NEG if d==NEG else n+d
def mdeg(mon,degs): return NEG if any(p and d==NEG for p,d in zip(mon,degs)) else sum(p*d for p,d in zip(mon,degs))
def option_maker(tabs):
    @lru_cache(None)
    def opts(f,degs):
        finite=[x for x in (mdeg(mon,degs) for mon in tabs[f]) if x!=NEG]
        if not finite: return (NEG,)
        top=max(finite); return (top,) if finite.count(top)==1 else (NEG,)+tuple(range(top+1))
    return opts

def witness(a,degs,opts):
    """One relaxed max-plus witness, or None; trace is ({H_f},{R_f})."""
    w,E=3*a-30,degs[3]-a; cur={}
    for h7 in opts(7,degs):
        if h7==NEG: cur[NEG]=({7:NEG},{6:NEG})
        elif h7>=w: cur[h7-w]=({7:h7},{6:h7-w})
    for f in range(6,0,-1):
        nxt={}
        for rf,(Hs,Rs) in cur.items():
            for hf in opts(f,degs):
                x,y=plus(3*(7-f)*E,hf),plus(4,rf)
                if x==y==NEG: outs=[NEG]
                elif x==y: outs=[NEG]+list(range(max(0,x-w+1)))
                else: outs=[max(x,y)-w] if max(x,y)-w>=0 else []
                for rp in outs:
                    if rp not in nxt:
                        hh,rr=dict(Hs),dict(Rs); hh[f]=hf; rr[f-1]=rp; nxt[rp]=(hh,rr)
        cur=nxt
        if not cur: return None
    for r0,(Hs,Rs) in cur.items():
        for h0 in opts(0,degs):
            if plus(21*E,h0)==plus(4,r0): hh=dict(Hs); hh[0]=h0; return hh,Rs
    return None

def check_witness(a,degs,W,opts):
    H,R=W; w,E=3*a-30,degs[3]-a
    need(all(H[f] in opts(f,degs) for f in range(8)),"witness forbidden H")
    need((H[7]==R[6]==NEG) or w+R[6]==H[7],"witness top anchor")
    for f in range(6,0,-1):
        x,y=plus(3*(7-f)*E,H[f]),plus(4,R[f]); rp=R[f-1]
        if rp==NEG: need(x==y,f"vanishing without tie f={f}")
        else: need(w+rp<=max(x,y) and (w+rp==max(x,y) or x==y),f"illegal drop f={f}")
    need(plus(21*E,H[0])==plus(4,R[0]),"bottom close not tie")

def jstate(row): return tuple(NEG if x is None else int(x) for x in row[:4])
def sample_state(x):
    d=x["degstate"]; return tuple(NEG if d[k] is None else int(d[k]) for k in ("deg_d2","deg_d1","deg_sigma","deg_e"))
def states(a,b,branch):
    return product([NEG]+list(range(7)),list(range(10)) if branch=="T1" else [NEG],
                   [NEG]+list(range(13)) if branch=="T1" else list(range(13)),range(a+sum(b),16))

def audit_sweep(data,open_ids,tabs):
    rows=data.get("branches",[]); ids={x.get("id") for x in rows}
    need(len(rows)==len(ids)==27 and ids==open_ids,f"branch list mismatch: missing={sorted(open_ids-ids)}, extra={sorted(ids-open_ids)}")
    opts=option_maker(tabs); total=surv=killed=spots=0; killed_pool=[]
    for row in rows:
        a=int(row["a"]); b=tuple(map(int,row["b"])); br=row["branch"]; name=bid(a,b,br)
        need(row["id"]==name and row.get("sum_b")==sum(b) and row.get("w")==3*a-30,f"{name}: metadata")
        need(row.get("deg_E_range")==[sum(b),15-a],f"{name}: E range")
        AS={jstate(x) for x in row["surviving_states_compact"]}; AK={jstate(x) for x in row["killed_states_compact"]}
        need(len(AS)==len(row["surviving_states_compact"]) and len(AK)==len(row["killed_states_compact"]) and not AS&AK,f"{name}: compact partition")
        IS,IK=set(),set()
        for raw in states(a,b,br):
            degs=tuple(map(int,raw)); (IS if witness(a,degs,opts) else IK).add(degs)
        need(IS==AS and IK==AK,f"{name}: partition mismatch survivor {len(AS-IS)}/{len(IS-AS)}, killed {len(AK-IK)}/{len(IK-AK)}")
        c=row["counts"]; need(c.get("total_degree_states")==len(AS)+len(AK) and c.get("surviving")==len(AS) and c.get("killed")==len(AK),f"{name}: counts")
        need(row.get("verdict")==("OPEN" if IS else "KILLED"),f"{name}: verdict")
        samples=row.get("survive_samples",[]); need(samples,f"{name}: no OPEN samples")
        for x in samples:
            degs=sample_state(x); need(degs in AS,f"{name}: sample absent")
            W=witness(a,degs,opts); need(W is not None,f"{name}: sample dies"); check_witness(a,degs,W,opts); spots+=1
        killed_pool += [(a,b,br,d) for d in AK]
        total+=len(AS)+len(AK); surv+=len(AS); killed+=len(AK)
    need((total,surv,killed)==(38360,4690,33670),f"sweep totals {(total,surv,killed)}")
    s=data.get("summary",{}); want={"n_branches":27,"branches_OPEN":27,"branches_KILLED":0,
        "total_degree_states":38360,"surviving_states":4690,"killed_states":33670}
    need(all(s.get(k)==v for k,v in want.items()),"sweep summary")
    rng=random.Random(0xA17E2026); n=min(1000,len(killed_pool)); need(n>=500,"fewer than 500 killed states")
    for a,b,br,degs in rng.sample(killed_pool,n): need(witness(a,degs,opts) is None,f"killed state survives: {bid(a,b,br)} {degs}")
    return total,surv,killed,spots,n

def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__); ap.add_argument("--quiet",action="store_true"); args=ap.parse_args(argv)
    start=time.perf_counter()
    try:
        tabs=parse_graded(); audit_algebra(); strata=audit_ledger(load("split_place_ledger_sub1.json"))
        k1,k2=derive_l1(strata); k3=derive_l2(strata,k1,tabs); kills=k1|k2|k3
        need(len(kills)==25,"combined kill total is not 25")
        opens={bid(a,b,br) for a,b in strata for br in ("T1","T2") if bid(a,b,br) not in kills}
        need(len(opens)==27,"residual total is not 27")
        total,surv,dead,spots,sample=audit_sweep(load("alt_inf_sweep.json"),opens,tabs)
        elapsed=time.perf_counter()-start
        if not args.quiet: print(f"PASS alternate-regime independent audit: 26 strata; kills 19+6=25; 27 OPEN branches; degree states {total}={surv} surviving+{dead} killed; {spots} OPEN samples and {sample} killed states checked; {elapsed:.2f}s")
        return 0
    except AuditError as exc: print(f"FAIL alternate-regime independent audit: {exc}",file=sys.stderr); return 1
    except Exception as exc: print(f"FAIL alternate-regime independent audit: unexpected {type(exc).__name__}: {exc}",file=sys.stderr); return 2

if __name__=="__main__": raise SystemExit(main())
