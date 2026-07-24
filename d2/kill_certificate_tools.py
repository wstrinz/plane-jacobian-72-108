#!/usr/bin/env python3
"""Produce portable cofactor certificates for saturated Groebner kills.

Manifest ``d2-kill-manifest-v1`` entries select documented state recipes.  Each
output is ``d2-kill-certificate-v1`` JSON.  Polynomials are canonical term lists:
``terms=[{coefficient:{numerator,denominator}, powers:[...]}, ...]`` in the exact
``variable_order`` stored beside them.  The producer primitive-normalizes each
input generator (positive leading coefficient), adds number-field and
Rabinowitsch constraints as ordinary generators, and asks Singular for a lift of
1.  Thus a consumer only has to check ``sum(cofactor[i]*generator[i]) == 1``.

The primary Singular route is ``lift(I,ideal(1))``.  If that fails, a standard-
basis transformation ``lift(I,G)*lift(G,ideal(1))`` is tried.  Every subprocess
has the manifest timeout (default 300 seconds).  Missing recipes, archived msolve
files, WSL, or Singular become honest NOT-YET-CERTIFICATED records; the batch
continues.  Existing repository files are read only.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import time
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parent
Q_COEFFS = (2048, -512, 320, -240, 195)
WSL_SINGULAR = ("wsl.exe", "-d", "Ubuntu", "--", "bash", "-lc",
                "cd $HOME && Singular -q")


def qpoly(v: sp.Symbol) -> sp.Expr:
    return sum(c*v**(4-i) for i, c in enumerate(Q_COEFFS))


def primitive(expr: sp.Expr, variables: list[sp.Symbol]) -> sp.Expr:
    """Canonical primitive integer associate; ideal-preserving over Q."""
    p = sp.Poly(sp.expand(sp.cancel(expr)), *variables, domain=sp.QQ)
    if p.is_zero:
        return sp.Integer(0)
    den = 1
    for c in p.coeffs():
        den = sp.ilcm(den, int(c.q))
    nums = [abs(int(c*den)) for c in p.coeffs() if c]
    content = math.gcd(*nums) if nums else 1
    out = sp.expand(p.as_expr()*sp.Rational(den, content))
    pp = sp.Poly(out, *variables, domain=sp.ZZ)
    if pp.LC() < 0:
        out = -out
    return sp.expand(out)


def poly_json(expr: sp.Expr, variables: list[sp.Symbol]) -> dict[str, Any]:
    p = sp.Poly(sp.expand(expr), *variables, domain=sp.QQ)
    terms = []
    for powers, coeff in p.terms():
        terms.append({"coefficient": {"numerator": str(coeff.p),
                                       "denominator": str(coeff.q)},
                      "powers": list(powers)})
    return {"terms": terms}


def poly_text(expr: sp.Expr, variables: list[sp.Symbol]) -> str:
    p = sp.Poly(sp.expand(expr), *variables, domain=sp.QQ)
    chunks = []
    for powers, coeff in p.terms():
        mon = "*".join(v.name + (f"^{e}" if e != 1 else "")
                       for v, e in zip(variables, powers) if e)
        if mon:
            if coeff == 1:
                term = mon
            elif coeff == -1:
                term = "-" + mon
            else:
                term = f"({coeff})*{mon}"
        else:
            term = str(coeff)
        chunks.append(term)
    return "+".join(chunks).replace("+-", "-") or "0"


def expr_string(expr: sp.Expr) -> str:
    return sp.sstr(sp.expand(expr))


def reduce_roots(expr: sp.Expr, roots: list[sp.Symbol]) -> sp.Expr:
    out = sp.expand(expr)
    for r in roots:
        out = sp.rem(sp.Poly(out, r), sp.Poly(qpoly(r), r)).as_expr()
    return sp.expand(out)


def generic_poly(prefix: str, degree: int, y: sp.Symbol) -> sp.Expr:
    cs = sp.symbols(f"{prefix}0:{degree+1}")
    return sum(cs[i]*y**i for i in range(degree+1))


def ring_vars(expressions: list[sp.Expr], y: sp.Symbol | None = None) -> list[sp.Symbol]:
    symbols: set[sp.Symbol] = set()
    for e in expressions:
        symbols.update(sp.sympify(e).free_symbols)
    if y is not None:
        symbols.discard(y)
    return sorted(symbols, key=lambda s: s.name)


def full_master(polys: dict[str, sp.Expr], targets: list[int], roots: list[sp.Symbol]) -> list[sp.Expr]:
    import convolution_descent as cd
    import convolution_elim as ce
    params = tuple(roots) + tuple(sorted(
        (set().union(*(p.free_symbols for p in polys.values())) - {cd.y}) -
        set().union(*(p.free_symbols for p in polys.values() if False)), key=str))
    # Explicit expressions determine unknowns; parameters only stop the source
    # engine's forcing logic, which coefficient extraction does not use.
    ans = cd.build_ansatz(d2=polys["d2"], d1=polys["d1"], e=polys["e"],
                          sigma=polys["sigma"], parameters=tuple(roots))
    eng = cd.ConvolutionDescent(ans, c=ce.DEFAULT_C)
    return [reduce_roots(eng.master_coefficient(d), roots) for d in targets]


def material_sys4(name: str) -> tuple[list[sp.Expr], dict[str, Any]]:
    m = re.fullmatch(r"sub2T2_a(\d+)_b(\d{4})_dd2(-?inf|\d+)_dsig(\d+)", name)
    if not m:
        raise ValueError(f"bad System-4 name {name}")
    a, btxt, d2txt, dsig = int(m[1]), m[2], m[3], int(m[4])
    bsum = sum(map(int, btxt)); y, r = sp.symbols("y r"); gamma, w = sp.symbols("gamma w")
    roots = [r] if bsum else []
    deg_g = dsig - 2*bsum
    G = generic_poly("g", deg_g, y)
    polys = {"d2": sp.Integer(0) if d2txt == "-inf" else generic_poly("a", int(d2txt), y),
             "d1": sp.Integer(0),
             "sigma": sp.expand((y-r)**(2*bsum)*G) if roots else G,
             "e": sp.expand(gamma*(y+1)**a*((y-r)**bsum if roots else 1))}
    coeffs = []
    import convolution_descent as cd
    import convolution_elim as ce
    ans = cd.build_ansatz(d2=polys["d2"], d1=0, e=polys["e"], sigma=polys["sigma"],
                          parameters=tuple(roots))
    eng = cd.ConvolutionDescent(ans, c=ce.DEFAULT_C)
    targets = []
    for degree in range(260, 195, -1):
        c = reduce_roots(eng.master_coefficient(degree), roots)
        if c != 0:
            coeffs.append(c); targets.append(degree)
        if len(coeffs) == 8:
            break
    satf = [gamma, sp.LC(sp.Poly(G, y))]
    if roots:
        satf.append(reduce_roots(G.subs(y, r), roots))
    members = coeffs + [reduce_roots(w*sp.prod(satf)-1, roots)] + [qpoly(r) for r in roots]
    mat = {"identity":"f31_master", "polynomials":{k:expr_string(v) for k,v in polys.items()},
           "targets":targets, "root_variables":[str(x) for x in roots],
           "saturation_factors":[expr_string(x) for x in satf], "source_name":name}
    return members, mat


def material_a8(name: str) -> tuple[list[sp.Expr], dict[str, Any]]:
    import convolution_descent as cd
    import convolution_elim as ce
    data = json.loads((ROOT/"batch_convolution_sub2.json").read_text())
    wanted = None
    for st in data["states"]:
        candidate = f"a8_dd2{st['deg_d2']}_dd1{st['deg_d1']}_dsig{st['deg_sigma']}"
        if candidate == name:
            wanted = st; break
    if wanted is None:
        raise KeyError(name)
    y = cd.y; gamma, w = sp.symbols("gamma w")
    d2 = sp.Integer(0) if wanted["d2_zero"] else generic_poly("a", int(wanted["deg_d2"]), y)
    d1 = generic_poly("b", int(wanted["deg_d1"]), y)
    sigma = generic_poly("s", int(wanted["deg_sigma"]), y)
    e = gamma*(y+1)**8
    polys = {"d2":d2,"d1":d1,"sigma":sigma,"e":e}
    ans = cd.build_ansatz(d2=d2,d1=d1,e=e,sigma=sigma,parameters=(gamma,))
    eng = ce.HighCoefficientEngine(ans, start_degree=int(wanted["gauge_detail"]["start"]),
                                    target_count=40, c=ce.DEFAULT_C)
    coeffs=[]; targets=[]; start=int(wanted["gauge_detail"]["start"])
    for degree in range(start,start-40,-1):
        c=sp.expand(eng.master_coefficient(degree))
        if c != 0:
            coeffs.append(c); targets.append(degree)
        if len(coeffs)==16: break
    members=coeffs+[w*gamma-1]
    mat={"identity":"f31_master","polynomials":{k:expr_string(v) for k,v in polys.items()},
         "targets":targets,"root_variables":[],"saturation_factors":["gamma"],"source_name":name}
    return members,mat


def material_sys3(name: str) -> tuple[list[sp.Expr], dict[str, Any]]:
    import phase_f2_scale as f2
    bid, suptxt = name.rsplit("_sup",1); support=int(suptxt)
    narrowed={x["key"] for x in json.loads((ROOT/"phase_f2_scale.json").read_text())["alt_states"] if str(x.get("verdict","")).startswith("NARROWED")}
    matches=[t for t in f2.load_targets() if t["bid"]==bid and int(t["support"])==support and f"{t['bid']}#sup{t['support']}#idx{t['idx']}" in narrowed]
    if len(matches)!=1:
        raise RuntimeError(f"expected one target for {name}, found {len(matches)}")
    t=matches[0]; degs=t["degs"]; drop_d1=t["branch"]=="T2"; drop_sig=t["sz"]
    TD=f2.total_deg(degs,drop_d1,drop_sig)
    factors,roots,scalars,Dc=f2.reconstruct(t["a"],t["b"],t["split"],t["branch"],degs,drop_d1,drop_sig)
    depth=min(int(t["depth"]),12)
    C=f2.h0_top(factors,tuple(x or 0 for x in degs),TD,depth,drop_d1=drop_d1,drop_sig=drop_sig)
    red=f2.make_reducer(roots); coeffs=[red(c) for c in C if c!=0]
    w=sp.Symbol("w"); sat_scalars=list(scalars)
    if Dc is not None and f2.d2_in_window(degs,TD,depth,drop_d1,drop_sig): sat_scalars.append(Dc[-1])
    members=coeffs+[w*sp.prod(sat_scalars)-1]+[f2.qpoly(r) for r in roots]
    if len(roots)==2:
        wd=sp.Symbol("wd"); members.append(wd*(roots[0]-roots[1])-1)
    polyd={k:(v[0] if isinstance(v,tuple) else v) for k,v in factors.items()}
    mat={"identity":"h0_tower","polynomials":{k:expr_string(v) for k,v in polyd.items()},
         "degrees":[x for x in degs],"TD":TD,"depth":depth,"drop_d1":drop_d1,
         "drop_sigma":drop_sig,"root_variables":[str(x) for x in roots],
         "saturation_factors":[expr_string(x) for x in sat_scalars],"source_name":name,
         "source_key":f"{bid}#sup{support}#idx{t['idx']}"}
    return members,mat


def material_d2(branch: str, degree: int, depth: int) -> tuple[list[sp.Expr], dict[str, Any]]:
    import d2_threshold as dt
    d2,sigma,e,roots,Dc=dt.build_state(branch,degree); red=dt.reducer(roots)
    coeffs=[red(c) for c in dt.h0_top(d2,degree,sigma,e,depth) if c!=0]
    satf=[dt.S,dt.E,Dc[-1]]; members=coeffs+[dt.w*sp.prod(satf)-1]+[dt.qpoly(r) for r in roots]
    if len(roots)==2: members.append(dt.wd*(roots[0]-roots[1])-1)
    mat={"identity":"h0_tower","polynomials":{"d2":expr_string(d2),"sigma":expr_string(sigma),"e":expr_string(e)},
         "degrees":[degree,None,12,15],"TD":60,"depth":depth,"drop_d1":True,"drop_sigma":False,
         "root_variables":[str(x) for x in roots],"saturation_factors":[expr_string(x) for x in satf],
         "source_name":branch}
    return members,mat


def material_phase(key: str, depth: int) -> tuple[list[sp.Expr], dict[str, Any]]:
    import phase_f2_sub2 as f2
    matches=[]
    for cell,case,st,idx,mx,pdelta in f2.load_targets(f2.TARGET_CELLS,max_defect=1):
        if f"{cell}#state{idx}"==key: matches.append((case,st,pdelta))
    if len(matches)!=1: raise RuntimeError(f"expected one phase target for {key}")
    case,st,pdelta=matches[0]; combo,why,_=f2.unique_split(case,st,pdelta)
    if combo is None: raise RuntimeError(why)
    polys,scalars,marked,mode,cofactors=f2.reconstruct(case,st,combo,pdelta)
    params=((f2.r,) if marked is not None else ())+tuple(cofactors)
    import convolution_descent as cd
    ans=cd.build_ansatz(d2=polys["d2"],d1=polys["d1"],e=polys["e"],sigma=polys["sigma"],parameters=params)
    eng=cd.ConvolutionDescent(ans,c=f2.C_VAL); top=f2.engine_top(eng)
    coeffs=[]; targets=[]
    for n in range(depth):
        c=f2.redq(eng.master_coefficient(top-n),marked)
        if c!=0: coeffs.append(c); targets.append(top-n)
    members=coeffs+[f2.w*sp.prod(scalars)-1]
    roots=[f2.r] if marked is not None else []
    members += [f2.QR_EXPR] if roots else []
    mat={"identity":"f31_master","polynomials":{k:expr_string(v) for k,v in polys.items()},
         "targets":targets,"root_variables":[str(x) for x in roots],
         "saturation_factors":[expr_string(x) for x in scalars],"source_key":key,
         "state_degrees":[st["deg_d1"],st["deg_sigma"],st["deg_d2"],st["deg_e"]],
         "d2_mode":mode,"cofactor_variables":[str(x) for x in cofactors]}
    return members,mat


def archived_ms(name: str) -> tuple[list[str],list[str]]|None:
    cmd=("wsl.exe","-d","Ubuntu","--","bash","-lc",f"test -f \"$HOME/{name}\" && cat \"$HOME/{name}\"")
    try:
        cp=subprocess.run(cmd,text=True,encoding="utf-8",stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=10)
    except Exception:
        return None
    if cp.returncode or not cp.stdout.strip(): return None
    lines=[x.strip() for x in cp.stdout.splitlines() if x.strip()]
    if len(lines)<3 or lines[1] != "0": return None
    body="".join(lines[2:]).rstrip(",")
    return lines[0].split(","),[x.strip() for x in re.split(r",\s*",body) if x.strip()]


def resolve(entry: dict[str,Any]) -> tuple[list[sp.Expr],dict[str,Any]]:
    r=entry["recipe"]; b=r["builder"]
    if b=="harvest_sys4": return material_sys4(r["name"])
    if b=="harvest_a8": return material_a8(r["name"])
    if b=="harvest_sys3": return material_sys3(r["name"])
    if b=="d2_threshold": return material_d2(r["branch"],int(r["degree_d2"]),int(r["depth"]))
    if b=="phase_f2_sub2": return material_phase(r["state_key"],int(r["depth"]))
    if b=="blowup_case":
        case=r["case"]
        if case=="a12_b1110_T2_d6": members,mat=material_d2("a12_b1110_T2",6,12)
        elif case=="a11_b1111_T1_17": members,mat=material_sys3("a11_b1111_T1_sup17")
        else:
            keys={"sub2_s14":"sub2:a9_b1000_T1_sz0_dz0_gz-#state14",
                  "sub2_s38":"sub2:a9_b1000_T1_sz0_dz0_gz-#state38",
                  "sub2_s94":"sub2:a9_b1000_T1_sz0_dz0_gz-#state94"}
            members,mat=material_phase(keys[case],6)
        arc=archived_ms(r.get("archived_ms",""))
        mat["archived_ms"]={"path":r.get("archived_ms"),"reused":bool(arc)}
        if arc:
            names,texts=arc; loc={n:sp.Symbol(n) for n in names}; parsed=[sp.sympify(x.replace("^","**"),locals=loc) for x in texts]
            # Reuse only after ideal generators match the regenerated primitive list.
            av=[loc[n] for n in names]; pv=[primitive(x,av) for x in parsed if x!=0]
            rv=ring_vars(members); mv=[primitive(x,rv) for x in members if x!=0]
            if [expr_string(x) for x in pv]==[expr_string(x) for x in mv]: members=parsed
            else: mat["archived_ms"]["mismatch"]="archived program differed; regenerated recipe used"
        return members,mat
    raise ValueError(f"unknown builder {b}")


def emit_program(gens:list[sp.Expr],variables:list[sp.Symbol],fallback:bool=False)->str:
    # Generator identifiers must not shadow ring variables: Singular silently
    # (no warning) resolves `poly g1 = ...` ahead of ring variable g1, so later
    # generator bodies referencing g1 substitute the earlier polynomial and the
    # ideal being lifted is corrupted.  Root cause writeup: CERT_LIFT_DEBUG.md.
    pfx="GEN"
    while any(v.name.startswith(pfx) for v in variables): pfx+="Z"
    lines=[f"ring R = 0,({','.join(v.name for v in variables)}),dp;"]
    for i,g in enumerate(gens): lines.append(f"poly {pfx}{i} = {poly_text(g,variables)};")
    lines.append("ideal I = "+",".join(f"{pfx}{i}" for i in range(len(gens)))+";")
    lines.append("ideal U = 1;")
    if fallback:
        lines += ["ideal G=std(I);","matrix T=lift(I,G);","matrix A=lift(G,U);","matrix L=T*A;"]
    else: lines.append("matrix L=lift(I,U);")
    for i in range(len(gens)):
        lines += [f'"@@COF_BEGIN_{i}";',f"L[{i+1},1];",f'"@@COF_END_{i}";']
    lines += ['"@@DONE";','quit;']
    return "\n".join(lines)+"\n"


def run_lift(gens:list[sp.Expr],variables:list[sp.Symbol],timeout:float)->dict[str,Any]:
    last=""
    for method,fallback in (("lift(I,ideal(1))",False),("std-transform-lift",True)):
        started=time.monotonic()
        try:
            cp=subprocess.run(WSL_SINGULAR,input=emit_program(gens,variables,fallback),text=True,
                              encoding="utf-8",stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                              timeout=timeout,check=False)
        except subprocess.TimeoutExpired:
            return {"status":"NOT-YET-CERTIFICATED","reason":f"Singular {method} timed out after {timeout:g}s",
                    "lift_wall_seconds":round(time.monotonic()-started,3),"lift_method":method}
        except Exception as ex:
            last=f"Singular launch failed: {type(ex).__name__}: {ex}"; continue
        wall=round(time.monotonic()-started,3); combined=(cp.stdout or "")+"\n"+(cp.stderr or "")
        if cp.returncode or "@@DONE" not in combined:
            last=f"{method} failed (exit {cp.returncode}): "+combined.replace("\x00","").strip()[-800:]; continue
        cof=[]; loc={v.name:v for v in variables}
        # Singular emits SHORT monomial format (no '*' or '^': X8r4w8 means
        # X^8*r^4*w^8) exactly when every ring variable name is one character.
        # Expand it before sympify; the exact-expansion gate below still
        # guards against any mis-parse surviving.
        short_ok=all(len(v.name)==1 for v in variables)
        def _expand_short(s):
            out=[]; i=0
            while i<len(s):
                ch=s[i]
                if ch.isalpha():
                    j=i+1
                    while j<len(s) and s[j].isdigit(): j+=1
                    out.append(ch+("**"+s[i+1:j] if j>i+1 else ""))
                    if j<len(s) and s[j].isalpha(): out.append("*")
                    i=j
                else:
                    out.append(ch)
                    if ch.isdigit() and i+1<len(s) and s[i+1].isalpha(): out.append("*")
                    i+=1
            return "".join(out)
        try:
            for i in range(len(gens)):
                m=re.search(rf"@@COF_BEGIN_{i}\s*(.*?)\s*@@COF_END_{i}",combined,re.S)
                if not m: raise ValueError(f"missing cofactor {i}")
                txt="".join(line.strip() for line in m.group(1).splitlines())
                try:
                    cof.append(sp.sympify(txt.replace("^","**"),locals=loc))
                except Exception:
                    if not short_ok: raise
                    cof.append(sp.sympify(_expand_short(txt),locals=loc))
        except Exception as ex:
            last=f"{method} output parse failed: {ex}"; continue
        if sp.expand(sum(c*g for c,g in zip(cof,gens))-1)!=0:
            last=f"{method} returned cofactors that failed producer-side exact expansion"; continue
        return {"status":"CERTIFICATE-FOUND","lift_method":method,"lift_wall_seconds":wall,"cofactors":cof}
    return {"status":"NOT-YET-CERTIFICATED","reason":last or "Singular lift failed"}


def safe_name(kill_id:str)->str:
    return re.sub(r"[^A-Za-z0-9_.-]+","__",kill_id)+".json"


def preflight()->str|None:
    try:
        cp=subprocess.run(WSL_SINGULAR,input='ring R=0,(x),dp; "@@OK"; quit;\n',text=True,
                          encoding="utf-8",stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=20)
    except Exception as ex: return f"Singular preflight launch failed: {type(ex).__name__}: {ex}"
    text=(cp.stdout or "")+(cp.stderr or "")
    if cp.returncode or "@@OK" not in text: return f"Singular preflight failed (exit {cp.returncode}): {text.replace(chr(0),'').strip()[-500:]}"
    return None


def main()->int:
    ap=argparse.ArgumentParser(); ap.add_argument("--manifest",default="kill_manifest.json")
    ap.add_argument("--all",action="store_true",help="attempt every manifest entry")
    ap.add_argument("--id",action="append",default=[]); ap.add_argument("--timeout",type=float)
    ns=ap.parse_args(); mp=(ROOT/ns.manifest).resolve(); manifest=json.loads(mp.read_text(encoding="utf-8-sig"))
    entries=manifest["entries"]
    if ns.id: entries=[e for e in entries if e["id"] in set(ns.id)]
    elif not ns.all: ap.error("pass --all or one or more --id")
    outdir=(ROOT/manifest.get("output_dir","kill_certificates")); outdir.mkdir(parents=True,exist_ok=True)
    timeout=float(ns.timeout or manifest.get("timeout_seconds",300)); unavailable=preflight(); statuses=[]
    for n,e in enumerate(entries,1):
        print(f"[{n}/{len(entries)}] {e['id']}",flush=True); base={"schema":"d2-kill-certificate-v1","kill_id":e["id"],
              "category":e["category"],"manifest_recipe":e["recipe"],"timeout_seconds":timeout}
        started=time.monotonic()
        if unavailable:
            result={**base,"status":"NOT-YET-CERTIFICATED","reason":unavailable,"attempt_wall_seconds":round(time.monotonic()-started,3)}
        else:
            try:
                raw,material=resolve(e); variables=ring_vars(raw); gens=[primitive(g,variables) for g in raw if g!=0]
                lr=run_lift(gens,variables,timeout); result={**base,"status":lr["status"],"variable_order":[v.name for v in variables],
                    "generator_normalization":"primitive integer associate with positive leading coefficient",
                    "generating_recipe":material,"generators":[poly_json(g,variables) for g in gens],
                    "generator_sha256":hashlib.sha256(json.dumps([poly_json(g,variables) for g in gens],sort_keys=True,separators=(",",":")).encode()).hexdigest(),
                    "lift_method":lr.get("lift_method"),"lift_wall_seconds":lr.get("lift_wall_seconds"),
                    "attempt_wall_seconds":round(time.monotonic()-started,3)}
                if lr["status"]=="CERTIFICATE-FOUND": result["cofactors"]=[poly_json(c,variables) for c in lr["cofactors"]]
                else: result["reason"]=lr.get("reason","lift failed")
            except Exception as ex:
                result={**base,"status":"NOT-YET-CERTIFICATED","reason":f"recipe/build error: {type(ex).__name__}: {ex}",
                        "attempt_wall_seconds":round(time.monotonic()-started,3)}
        path=outdir/safe_name(e["id"]); path.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
        statuses.append({k:result.get(k) for k in ("kill_id","category","status","reason","lift_method","lift_wall_seconds","attempt_wall_seconds")})
        print(f"  {result['status']}: {result.get('reason','lift '+str(result.get('lift_wall_seconds'))+'s')}",flush=True)
    log={"schema":"d2-kill-status-log-v1","manifest":mp.name,"timeout_seconds":timeout,"entries":statuses}
    logpath=ROOT/manifest.get("status_log","kill_certificates/status_log.json"); logpath.parent.mkdir(parents=True,exist_ok=True)
    logpath.write_text(json.dumps(log,indent=2,sort_keys=True)+"\n")
    from collections import Counter
    print("PRODUCTION CENSUS",dict(Counter(x["status"] for x in statuses)),flush=True)
    return 0


if __name__=="__main__": raise SystemExit(main())




