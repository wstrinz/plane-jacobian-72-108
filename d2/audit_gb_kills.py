#!/usr/bin/env python3
"""Independent, spec-only consumer for ``d2-kill-certificate-v1`` JSON.

This file deliberately does not import the producer or any state builder.  It
implements its own canonical term-list decoder, exact SymPy polynomial
arithmetic, root-minimal-polynomial reduction, and a from-scratch parser for all
``h_f`` lines in ``f31_graded.txt``.  For a certificate it performs two gates:

1. independently rebuild the documented recipe generators from the recorded
   ansatz plus the repository state census, and compare their canonical term
   lists with the certificate's generators;
2. expand ``sum(c_i*f_i)-1`` exactly (then reduce adjoined root variables modulo
   q as an additional number-field check).

Any recipe mismatch or nonzero residual is EXPANSION-FAILED and is printed
prominently.  NOT-YET-CERTIFICATED records are UNCERTIFICATED and do not make
the process fail.  ``--quiet`` suppresses per-kill success/pending lines, never
failures or the final census.  No ideal-membership or basis computation occurs.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
import math
from pathlib import Path
import re
import sys
from typing import Any

import sympy as sp

ROOT=Path(__file__).resolve().parent
Q_COEFFS=(2048,-512,320,-240,195)
y=sp.Symbol("y")


def qpoly(v:sp.Symbol)->sp.Expr:
    return sum(c*v**(4-i) for i,c in enumerate(Q_COEFFS))


def decode(poly:dict[str,Any],variables:list[sp.Symbol])->sp.Expr:
    out=sp.Integer(0)
    for term in poly.get("terms",[]):
        powers=term["powers"]
        if len(powers)!=len(variables): raise ValueError("power vector length differs from variable_order")
        c=term["coefficient"]; den=int(c["denominator"]); num=int(c["numerator"])
        if den<=0: raise ValueError("coefficient denominator must be positive")
        mon=sp.Rational(num,den)
        for v,e in zip(variables,powers):
            if not isinstance(e,int) or e<0: raise ValueError("powers must be nonnegative integers")
            mon*=v**e
        out+=mon
    return sp.expand(out)


def primitive(expr:sp.Expr,variables:list[sp.Symbol])->sp.Expr:
    p=sp.Poly(sp.expand(sp.cancel(expr)),*variables,domain=sp.QQ)
    if p.is_zero:return sp.Integer(0)
    den=1
    for c in p.coeffs(): den=sp.ilcm(den,int(c.q))
    nums=[abs(int(c*den)) for c in p.coeffs() if c]
    content=math.gcd(*nums) if nums else 1
    out=sp.expand(p.as_expr()*sp.Rational(den, content))
    if sp.Poly(out,*variables,domain=sp.ZZ).LC()<0:out=-out
    return sp.expand(out)


def parse_f31()->dict[int,sp.Expr]:
    """Parse only the eight authoritative h_f assignments, without eval/imports."""
    text=(ROOT/"f31_graded.txt").read_text(encoding="utf-8")
    d0,d1,d2,dm1=sp.symbols("d0 d1 d2 dm1")
    loc={"d0":d0,"d1":d1,"d2":d2,"dm1":dm1}
    out={}
    for line in text.splitlines():
        m=re.match(r"h_(\d+) \([^)]*\) = (.+)$",line.strip())
        if m: out[int(m[1])]=sp.sympify(m[2],locals=loc)
    if sorted(out)!=list(range(8)):raise ValueError(f"expected h_0..h_7; found {sorted(out)}")
    return out


Sparse=dict[int,sp.Expr]


def sparse(expr:sp.Expr)->Sparse:
    return {m[0]:c for m,c in sp.Poly(sp.expand(expr),y).terms()}


def sadd(a:Sparse,b:Sparse)->Sparse:
    z=dict(a)
    for k,v in b.items():z[k]=z.get(k,0)+v
    return {k:sp.expand(v) for k,v in z.items() if v!=0}


def smul(a:Sparse,b:Sparse)->Sparse:
    z={}
    for i,x in a.items():
        for j,v in b.items():z[i+j]=z.get(i+j,0)+x*v
    return {k:sp.expand(v) for k,v in z.items() if v!=0}


def spow(a:Sparse,n:int)->Sparse:
    z={0:sp.Integer(1)};b=dict(a)
    while n:
        if n&1:z=smul(z,b)
        n//=2
        if n:b=smul(b,b)
    return z


def eval_sparse(expr:sp.Expr,subs:dict[sp.Symbol,Sparse])->Sparse:
    vars_=tuple(subs); z={}
    for powers,c in sp.Poly(expr,*vars_).terms():
        t={0:c}
        for v,n in zip(vars_,powers):
            if n:t=smul(t,spow(subs[v],n))
        z=sadd(z,t)
    return z


def reduce_roots(expr:sp.Expr,roots:list[sp.Symbol])->sp.Expr:
    out=sp.expand(expr)
    for r in roots:out=sp.rem(sp.Poly(out,r),sp.Poly(qpoly(r),r)).as_expr()
    return sp.expand(out)


def parse_material_polys(material:dict[str,Any])->tuple[dict[str,sp.Expr],dict[str,sp.Symbol]]:
    names=set(["y","d0","d1","d2","dm1","sigma","sig","e","w","wd","r","r1","r2"])
    for txt in material.get("polynomials",{}).values():names.update(re.findall(r"\b[A-Za-z_]\w*\b",txt))
    for txt in material.get("saturation_factors",[]):names.update(re.findall(r"\b[A-Za-z_]\w*\b",txt))
    loc={n:(y if n=="y" else sp.Symbol(n)) for n in names}
    polys={k:sp.sympify(v,locals=loc) for k,v in material.get("polynomials",{}).items()}
    return polys,loc


def independent_recipe(material:dict[str,Any],variables:list[sp.Symbol],h:dict[int,sp.Expr])->list[sp.Expr]:
    polys,loc=parse_material_polys(material); roots=[loc[n] for n in material.get("root_variables",[])]
    d0s,d1s,d2s,dm1s=sp.symbols("d0 d1 d2 dm1")
    sigma=polys.get("sigma",polys.get("sig",sp.Integer(0)))
    d2=polys.get("d2",sp.Integer(0));d1=polys.get("d1",sp.Integer(0));e=polys["e"]
    d0=sp.expand((d2**2+sigma)/4)
    source_sub={d0s:sparse(d0),d1s:sparse(d1),d2s:sparse(d2),dm1s:sparse(e)}
    hsp={f:eval_sparse(h[f],source_sub) for f in range(8)}
    coeffs=[]
    if material["identity"]=="f31_master":
        phi=sp.Rational(-1,6630)*(y+1)**30*qpoly(y); ph=sparse(phi); ep=sparse(e)
        phpow={f:spow(ph,f) for f in range(8)}; epow={f:spow(ep,21-3*f) for f in range(8)}
        for target in material["targets"]:
            total=sp.Integer(0)
            for f in range(8):
                pe=smul(phpow[f],epow[f])
                for deg,c in pe.items(): total+=c*hsp[f].get(int(target)-deg,0)
            c=reduce_roots(total,roots)
            if c!=0:coeffs.append(c)
    elif material["identity"]=="h0_tower":
        full=hsp[0]; TD=int(material["TD"])
        for i in range(int(material["depth"])):
            c=reduce_roots(full.get(TD-i,0),roots)
            if c!=0:coeffs.append(c)
    else:raise ValueError(f"unsupported identity {material['identity']}")
    sat=[sp.sympify(x,locals=loc) for x in material.get("saturation_factors",[])]
    wv=loc.get("w",sp.Symbol("w")); coeffs.append(wv*sp.prod(sat)-1)
    coeffs += [qpoly(r) for r in roots]
    if len(roots)==2:
        coeffs.append(loc.get("wd",sp.Symbol("wd"))*(roots[0]-roots[1])-1)
    # Map same-name independently created symbols to certificate symbols.
    cmap={str(v):v for v in variables}; subst={s:cmap[str(s)] for e0 in coeffs for s in e0.free_symbols if str(s) in cmap}
    return [primitive(e0.xreplace(subst),variables) for e0 in coeffs if e0!=0]


def documented_recipe_check(cert:dict[str,Any])->list[str]:
    """Cross-check IDs/depth/degrees against immutable census artifacts."""
    errors=[]; cat=cert["category"]; recipe=cert["manifest_recipe"]; mat=cert.get("generating_recipe",{})
    if cat=="harvest":
        data=json.loads((ROOT/"triage_harvest.json").read_text())
        names={x["name"] for group in (data.get("system3",[]),data.get("system4",[])) for x in group if x.get("exact_kill")}
        names.update(x["name"] for x in data.get("a8",{}).get("states",[]) if x.get("exact_kill"))
        if recipe.get("name") not in names:errors.append("state is not an exact_kill in triage_harvest.json")
        if recipe.get("builder")=="harvest_sys3":
            rows=json.loads((ROOT/"phase_f2_scale.json").read_text())["alt_states"]
            if mat.get("source_key") not in {x.get("key") for x in rows}:errors.append("source_key absent from phase_f2_scale.json")
    elif cat=="msolve_blowup":
        rows=json.loads((ROOT/"msolve_bridge_results.json").read_text())
        cases={x["result"].get("case") for x in rows if x["result"].get("char")==0 and x["result"].get("verdict")=="EMPTY(KILL)"}
        if recipe.get("case") not in cases:errors.append("case is not a rational EMPTY(KILL) in msolve_bridge_results.json")
    elif cat=="phase_f2_sub2":
        rows=json.loads((ROOT/"phase_f2_sub2.json").read_text())["states"]
        hit=[x for x in rows if x.get("key")==recipe.get("state_key")]
        if len(hit)!=1 or hit[0].get("verdict")!="KILLED":errors.append("state is not uniquely KILLED in phase_f2_sub2.json")
        elif int(hit[0].get("kill_depth",-1))!=int(recipe.get("depth",-2)):errors.append("manifest depth differs from documented kill_depth")
        elif mat.get("state_degrees")!=hit[0].get("degs"):errors.append("recorded recipe degrees differ from state census")
    elif cat=="d2_threshold":
        if recipe.get("branch") not in {"a12_b1110_T2","a11_b3100_T2"} or recipe.get("degree_d2")!=5 or recipe.get("depth")!=8:
            errors.append("entry differs from the two documented degree-5/depth-8 threshold kills")
    else:errors.append(f"unknown category {cat}")
    return errors


def audit_one(path:Path,h:dict[int,sp.Expr])->tuple[str,list[str]]:
    cert=json.loads(path.read_text())
    if cert.get("schema")!="d2-kill-certificate-v1":return "EXPANSION-FAILED",["wrong certificate schema"]
    if cert.get("status")!="CERTIFICATE-FOUND":return "UNCERTIFICATED",[cert.get("reason","no certificate")]
    try:
        variables=[sp.Symbol(n) for n in cert["variable_order"]]
        generators=[decode(x,variables) for x in cert["generators"]]
        cofactors=[decode(x,variables) for x in cert["cofactors"]]
        if len(generators)!=len(cofactors):return "EXPANSION-FAILED",["generator/cofactor length mismatch"]
        issues=documented_recipe_check(cert)
        derived=independent_recipe(cert["generating_recipe"],variables,h)
        if len(derived)!=len(generators):issues.append(f"generator count mismatch: certificate {len(generators)}, recipe {len(derived)}")
        else:
            bad=[i for i,(a,b) in enumerate(zip(generators,derived)) if sp.expand(a-b)!=0]
            if bad:issues.append("generator-recipe mismatch at indices "+",".join(map(str,bad)))
        residual=sp.Poly(-1,*variables,domain=sp.QQ)
        for c,g in zip(cofactors,generators):residual += sp.Poly(sp.expand(c*g),*variables,domain=sp.QQ)
        exact=residual.as_expr()
        roots=[sp.Symbol(n) for n in cert.get("generating_recipe",{}).get("root_variables",[])]
        reduced=reduce_roots(exact,roots)
        if exact!=0:issues.append("cofactor identity residual is nonzero")
        if reduced!=0:issues.append("number-field-reduced residual is nonzero")
        return ("EXPANSION-FAILED",issues) if issues else ("CERTIFIED",[])
    except Exception as ex:return "EXPANSION-FAILED",[f"auditor exception: {type(ex).__name__}: {ex}"]


def main()->int:
    ap=argparse.ArgumentParser();ap.add_argument("--manifest",default="kill_manifest.json");ap.add_argument("--quiet",action="store_true")
    ns=ap.parse_args();manifest=json.loads((ROOT/ns.manifest).read_text(encoding="utf-8-sig"));outdir=ROOT/manifest.get("output_dir","kill_certificates")
    h=parse_f31(); census=Counter(); mismatches=[]
    for entry in manifest["entries"]:
        fn=re.sub(r"[^A-Za-z0-9_.-]+","__",entry["id"])+".json";path=outdir/fn
        if not path.exists():status,issues="UNCERTIFICATED",["certificate output missing"]
        else:status,issues=audit_one(path,h)
        census[status]+=1
        if status=="EXPANSION-FAILED":
            mismatches.append((entry["id"],issues));print(f"!!! EXPANSION-FAILED {entry['id']}: {'; '.join(issues)}",file=sys.stderr)
        elif not ns.quiet:print(f"{status:16s} {entry['id']}"+(f" -- {issues[0]}" if issues else ""))
    print(f"AUDIT CENSUS CERTIFIED={census['CERTIFIED']} EXPANSION-FAILED={census['EXPANSION-FAILED']} UNCERTIFICATED={census['UNCERTIFICATED']}")
    if mismatches:
        print("GENERATOR-RECIPE MISMATCH / IDENTITY FINDINGS:",file=sys.stderr)
        for kid,issues in mismatches:print(f"  {kid}: {'; '.join(issues)}",file=sys.stderr)
    return 0 if census["EXPANSION-FAILED"]==0 else 1


if __name__=="__main__":raise SystemExit(main())

