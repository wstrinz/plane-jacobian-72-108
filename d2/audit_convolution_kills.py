#!/usr/bin/env python3
"""Spec-only independent audit of the 60 candidate convolution kills.

Reads only f31_graded.txt and the JSON artifact. It independently parses h_f,
checks homogeneity before and after d0=(d2**2+sigma)/4, reconstructs every
degree-exact ansatz, and walks exact coefficients of the master identity.

Degree-exactness makes gamma and all stated leading coefficients NONZERO by
definition. We write e as gamma*(y+1)**a times a monic generic polynomial.
Thus forcing a leading coefficient to zero is a state kill.

Only the requested moves are used: unit*X**k; unit*(A*X+B)**2 with A a unit;
and exact incompatibility of residual equations in nonzero parameters. A unit
is a nonzero rational times a Laurent monomial in declared-nonzero symbols.
Other stalls are conservative UNDECIDED-BY-AUDIT results. Each record runs in
a fresh child interpreter under a hard parent-side timeout.
"""
from __future__ import annotations
import argparse, json, math, os, re, subprocess, sys, time
from dataclasses import dataclass
from pathlib import Path
import sympy as sp

HERE=Path(__file__).resolve().parent
GRADED_FILE="f31_graded.txt"
DATA_FILE="batch_convolution_sub2.json"
DEFAULT_TIMEOUT=120.0
INITIAL_DEPTH=24

class AuditError(RuntimeError): pass
def need(ok,msg):
    if not ok: raise AuditError(msg)

@dataclass(frozen=True)
class Graded:
    terms: dict
    counts: tuple

@dataclass
class Ansatz:
    reverse: dict
    degrees: dict
    symbols: set
    nonzero: set
    leading: set
    gamma: sp.Symbol
    description: str

def parse_graded():
    """Own parser and original/sigma-coordinate homogeneity checks."""
    text=(HERE/GRADED_FILE).read_text(encoding="utf-8-sig")
    need("f31 = sum_{f=0}^7 Phi^f * dm1^(21-3f) * h_f" in text,
         "master-identity header missing")
    rows=re.findall(r"^h_(\d+) \(weight (\d+), dm1-power (\d+)\) = (.+)$",text,re.M)
    need(len(rows)==8,f"expected 8 h_f rows, found {len(rows)}")
    d2,d1,d0,e,sigma=sp.symbols("d2 d1 d0 dm1 sigma")
    loc={"d2":d2,"d1":d1,"d0":d0,"dm1":e}; result={}; counts={}
    for fs,ws,ps,source in rows:
        f,w,p=int(fs),int(ws),int(ps); expected=20-2*f
        need(w==expected and p==21-3*f,f"h_{f}: bad metadata")
        need(f not in result,f"duplicate h_{f}")
        try: original=sp.Poly(sp.sympify(source,locals=loc),d2,d1,d0,e,domain=sp.QQ)
        except Exception as exc: raise AuditError(f"h_{f}: parse failed: {exc}") from exc
        need(bool(original.terms()),f"h_{f} is zero")
        for mon,c in original.terms():
            actual=sum(a*b for a,b in zip(mon,(2,3,4,5)))
            need(c!=0 and actual==expected,f"h_{f}: monomial {mon} has weight {actual}")
        rewritten=sp.Poly(sp.expand(original.as_expr().subs(d0,(d2**2+sigma)/4)),
                          d2,d1,sigma,e,domain=sp.QQ)
        packed=[]
        for mon,c in rewritten.terms():
            actual=sum(a*b for a,b in zip(mon,(2,3,4,5)))
            need(c!=0 and actual==expected,f"sigma h_{f}: monomial {mon} has weight {actual}")
            packed.append((tuple(map(int,mon)),sp.Rational(c)))
        result[f]=tuple(packed); counts[f]=len(packed)
    need(set(result)==set(range(8)),f"h indices are {sorted(result)}")
    polys={f:sp.Poly.from_dict(dict(result[f]),(d2,d1,sigma,e),domain=sp.QQ) for f in (5,6,7)}
    h7=8192*d1**2
    h6=14336*d1**2*d2+8192*d1*e-3072*sigma**2
    h5=(-12288*d1**2*d2**2+32256*d1**2*sigma+18432*d1*d2*e
        -9216*d2*sigma**2+2048*e**2)
    need(sp.expand(polys[7].as_expr()-h7)==0,"h_7 spot check failed")
    need(sp.expand(polys[6].as_expr()-h6)==0,"h_6 spot check failed")
    need(sp.expand(polys[5].as_expr()-h5)==0,"h_5 spot check failed")
    return Graded(result,tuple(counts[f] for f in range(8)))

def load_records():
    data=json.loads((HERE/DATA_FILE).read_text(encoding="utf-8-sig"))
    need(isinstance(data,dict),"JSON root is not an object")
    rows=data.get("kills_pending_audit")
    need(isinstance(rows,list),"kills_pending_audit is not a list")
    need(len(rows)==60,f"expected 60 records, found {len(rows)}")
    need(data.get("kill_count")==60,"kill_count is not 60")
    return data,rows

def state_degree(row,key,zero_key):
    value,zero=row.get(key),bool(row.get(zero_key))
    if value in ("-inf",None):
        need(zero,f"{key}=-inf but {zero_key}=false"); return None
    need(not zero,f"{key} finite but {zero_key}=true")
    value=int(value); need(value>=0,f"{key} is negative"); return value

def conv(left,right,depth):
    size=min(depth+1,len(left)+len(right)-1); out=[sp.S.Zero]*size
    for i,a in enumerate(left):
        if i>=size or a==0: continue
        for j in range(min(len(right),size-i)):
            if right[j]!=0: out[i+j]+=a*right[j]
    return out

def reverse_generic(prefix,degree):
    symbols=list(sp.symbols(f"{prefix}0:{degree+1}"))
    symbols[-1]=sp.Symbol(f"{prefix}{degree}",nonzero=True)
    return list(reversed(symbols)),set(symbols),symbols[-1]

def reverse_e(a,m,gamma,depth):
    shifted=[sp.Integer(math.comb(a,k)) for k in range(a+1)]
    lower=list(sp.symbols(f"z0:{m}")) if m else []
    return [gamma*x for x in conv(shifted,[sp.S.One,*reversed(lower)],depth)],set(lower)

def reconstruct(row,depth):
    branch=row.get("branch"); need(branch in ("T1","T2"),"unknown branch")
    need(bool(row.get("d1_zero"))==(branch=="T2"),"branch/d1_zero mismatch")
    degrees={"d2":state_degree(row,"deg_d2","d2_zero"),
             "d1":state_degree(row,"deg_d1","d1_zero"),
             "sigma":state_degree(row,"deg_sigma","sigma_zero"),
             "e":int(row["deg_e"])}
    a=int(row["a_t"]); m=degrees["e"]-a
    need(m>=0 and row.get("m")==m,"bad generic-e degree")
    reverse={}; symbols=set(); nonzero=set(); leading=set()
    for name,prefix in (("d2","a"),("d1","b"),("sigma","s")):
        degree=degrees[name]
        if degree is None: reverse[name]=None
        else:
            co,syms,lead=reverse_generic(prefix,degree); reverse[name]=co
            symbols|=syms; nonzero.add(lead); leading.add(lead)
    gamma=sp.Symbol("gamma",nonzero=True)
    reverse["e"],lower=reverse_e(a,m,gamma,depth)
    symbols|=lower|{gamma}; nonzero.add(gamma); leading.add(gamma)
    def show(name): return f"{name}=0" if degrees[name] is None else f"deg({name})={degrees[name]} exact"
    desc=", ".join((f"branch={branch}",show("d2"),show("d1"),show("sigma"),
                    f"e=gamma*(y+1)^{a}*monic-generic(deg={m}), deg(e)={degrees['e']}"))
    return Ansatz(reverse,degrees,symbols,nonzero,leading,gamma,desc)

def power(base,exponent,depth):
    need(exponent>=0,"negative polynomial power")
    result=[sp.S.One]; factor=base; n=exponent
    while n:
        if n&1: result=conv(result,factor,depth)
        n//=2
        if n: factor=conv(factor,factor,depth)
    return result

def phi_reverse(depth):
    y=sp.Symbol("y"); q=2048*y**4-512*y**3+320*y**2-240*y+195
    p=sp.Poly(sp.Rational(-1,6630)*(y+1)**30*q,y,domain=sp.QQ)
    need(p.degree()==34,"Phi degree is not 34")
    return [p.nth(34-k) for k in range(min(depth,34)+1)]

def master_top(graded,ansatz,depth):
    """Exact reverse sparse convolution retaining the requested top slice."""
    phi=phi_reverse(depth); names=("d2","d1","sigma","e"); cache={}
    def getpow(name,exponent):
        key=(name,exponent)
        if key not in cache:
            base=phi if name=="phi" else ansatz.reverse[name]
            need(base is not None,f"positive power of zero {name}")
            cache[key]=power(base,exponent,depth)
        return cache[key]
    terms=[]; top=-1
    for f in range(8):
        outer=21-3*f
        for exponents,c in graded.terms[f]:
            if any(exponents[i] and ansatz.degrees[name] is None for i,name in enumerate(names)):
                continue
            degree=34*f+outer*ansatz.degrees["e"]
            degree+=sum(exponents[i]*ansatz.degrees[name] for i,name in enumerate(names)
                        if exponents[i])
            top=max(top,degree); terms.append((f,exponents,c,degree))
    need(top>=0 and terms,"all master terms vanished")
    out=[sp.S.Zero]*(depth+1)
    for f,exponents,c,degree in terms:
        offset=top-degree
        if offset>depth: continue
        current=[c]
        if f: current=conv(current,getpow("phi",f),depth-offset)
        for i,name in enumerate(names):
            exponent=exponents[i]+(21-3*f if name=="e" else 0)
            if exponent: current=conv(current,getpow(name,exponent),depth-offset)
        for k,value in enumerate(current):
            if offset+k<=depth: out[offset+k]+=value
    return top,out

def is_unit(expression,nonzero):
    expression=sp.factor(expression); coefficient,rest=expression.as_coeff_Mul()
    if not coefficient.is_Rational or coefficient==0: return False
    return all(base==1 or (base in nonzero and exponent.is_Integer)
               for base,exponent in rest.as_powers_dict().items())

def substituted(expression,substitutions):
    for symbol,value in substitutions: expression=expression.subs(symbol,value)
    return sp.factor(sp.cancel(expression))

def factored_numerator(expression):
    numerator,denominator=sp.fraction(sp.cancel(expression))
    coefficient,factors=sp.factor_list(numerator)
    return sp.sympify(coefficient),denominator,[(sp.sympify(a),int(b)) for a,b in factors]

def forced_move(expression,active,nonzero,gamma):
    """Recognize only unit*X**k and unit*(linear-in-X)**2."""
    coefficient,denominator,factors=factored_numerator(expression)
    if not is_unit(coefficient/denominator,nonzero): return None
    # Try a monomial force before regarding a declared leading symbol as a unit.
    for target in sorted(active-{gamma},key=str):
        rest=coefficient; exponent=0
        for factor,multiplicity in factors:
            if factor==target: exponent+=multiplicity
            else: rest*=factor**multiplicity
        if exponent and is_unit(rest/denominator,nonzero-{target}):
            return target,sp.S.Zero,f"unit*{target}**{exponent}"
    for position,(factor,multiplicity) in enumerate(factors):
        if multiplicity!=2: continue
        others=coefficient
        for other_position,(other,n) in enumerate(factors):
            if other_position!=position: others*=other**n
        if not is_unit(others/denominator,nonzero): continue
        for target in sorted(active-{gamma},key=str):
            try: polynomial=sp.Poly(factor,target)
            except sp.PolynomialError: continue
            if polynomial.degree()!=1: continue
            slope,intercept=polynomial.nth(1),polynomial.nth(0)
            if not is_unit(slope,nonzero-{target}): continue
            value=sp.factor(sp.cancel(-intercept/slope))
            return target,value,f"unit*({sp.factor(factor)})**2"
    return None

def strip_nonzero_monomial(expression,nonzero):
    numerator,_=sp.fraction(sp.cancel(expression))
    symbols=sorted(numerator.free_symbols,key=str)
    if not symbols: return sp.factor(numerator)
    poly=sp.Poly(numerator,*symbols,domain=sp.QQ)
    minima=[min(mon[i] for mon,_ in poly.terms()) for i in range(len(symbols))]
    divisor=sp.S.One
    for symbol,exponent in zip(symbols,minima):
        if symbol in nonzero and exponent: divisor*=symbol**exponent
    return sp.factor(sp.cancel(numerator/divisor))

def residual_contradiction(residuals,nonzero):
    normalized=[strip_nonzero_monomial(x,nonzero) for x in residuals]
    for item in normalized:
        if is_unit(item,nonzero):
            return True,f"nonzero Laurent monomial residual {item}"
    involved=set().union(*(x.free_symbols for x in normalized)) if normalized else set()
    need(involved<=nonzero,"residual contains an unrestricted unknown")
    if len(involved)==1 and len(normalized)>=2:
        symbol=next(iter(involved)); gcd=sp.Poly(normalized[0],symbol,domain=sp.QQ)
        for item in normalized[1:]: gcd=sp.gcd(gcd,sp.Poly(item,symbol,domain=sp.QQ))
        gcd_expr=sp.factor(gcd.as_expr())
        if not strip_nonzero_monomial(gcd_expr,{symbol}).free_symbols:
            return True,f"univariate residual gcd is {gcd_expr} (no nonzero root)"
    # Sound multivariate fallback: saturate by the product of nonzero symbols.
    if 1<len(involved)<=3 and 2<=len(normalized)<=8:
        variables=sorted(involved,key=str)
        maxdeg=max(sp.Poly(x,*variables).total_degree() for x in normalized)
        if maxdeg<=20:
            inverse=sp.Symbol("torus_inverse")
            equations=[*normalized,inverse*sp.prod(variables)-1]
            basis=sp.groebner(equations,inverse,*variables,order="grevlex")
            if basis.contains(sp.S.One):
                return True,"saturated nonzero-parameter ideal is the unit ideal"
    return False,""

def terminal_solution(residuals,nonzero,active):
    """Return only an explicit, exactly verified nonzero-parameter solution."""
    if active-nonzero: return False,None
    normalized=[strip_nonzero_monomial(x,nonzero) for x in residuals]
    involved=set().union(*(x.free_symbols for x in normalized)) if normalized else set()
    if not normalized: return True,"all coefficients vanish; set declared parameters to 1"
    if len(involved)==1:
        symbol=next(iter(involved)); gcd=sp.Poly(normalized[0],symbol,domain=sp.QQ)
        for item in normalized[1:]: gcd=sp.gcd(gcd,sp.Poly(item,symbol,domain=sp.QQ))
        for root in sp.roots(gcd.as_expr(),symbol):
            if root!=0 and all(sp.simplify(item.subs(symbol,root))==0 for item in normalized):
                return True,f"verified common nonzero solution {symbol}={root}"
    return False,None

def identifier(index,row):
    def degree(key,zero): return "zero" if row.get(zero) else str(row.get(key))
    return (f"K{index+1:03d}:a{row.get('a_t')}:{row.get('branch')}:"
            f"d2={degree('deg_d2','d2_zero')}:d1={degree('deg_d1','d1_zero')}:"
            f"sigma={degree('deg_sigma','sigma_zero')}:e={row.get('deg_e')}")

def audit_record(index,row,graded):
    started=time.perf_counter(); name=identifier(index,row)
    expected=row.get("final_verdict")
    need(expected in ("CONTRADICTION","STATE_KILLED_BY_DEGREE_DROP"),f"{name}: bad claim")
    depth=INITIAL_DEPTH
    while True:
        ansatz=reconstruct(row,depth); top,coefficients=master_top(graded,ansatz,depth)
        substitutions=[]; residuals=[]; trace=[]; terminal=None
        for offset,raw in enumerate(coefficients):
            degree=top-offset; expression=substituted(raw,substitutions)
            if expression==0:
                if len(trace)<40: trace.append({"degree":degree,"move":"IDENTITY"})
                continue
            used={x for x,_ in substitutions}; active=ansatz.symbols-used
            move=forced_move(expression,active,ansatz.nonzero,ansatz.gamma)
            if move is not None:
                target,value,reason=move
                if target in ansatz.leading and value==0:
                    trace.append({"degree":degree,"move":"STATE_KILL","symbol":str(target),
                                  "residual":str(expression),"reason":reason})
                    terminal={"mechanism":"STATE_KILL","degree":degree,
                              "detail":f"declared-nonzero leading coefficient {target} forced to zero"}
                    break
                if target in ansatz.leading and not is_unit(value,ansatz.nonzero-{target}):
                    terminal={"mechanism":"UNDECIDED","degree":degree,
                              "detail":f"unresolved nonzero side condition: {target}={value}",
                              "residual":str(expression)}
                    trace.append({"degree":degree,"move":"STALL","residual":str(expression)})
                    break
                substitutions.append((target,value))
                trace.append({"degree":degree,"move":"FORCED",
                              "substitution":[str(target),str(value)],"reason":reason})
                continue
            if expression.free_symbols<=ansatz.nonzero:
                normalized=strip_nonzero_monomial(expression,ansatz.nonzero)
                residuals.append(normalized)
                trace.append({"degree":degree,"move":"PARAMETER_CONSTRAINT",
                              "residual":str(normalized)})
                contradiction,reason=residual_contradiction(residuals,ansatz.nonzero)
                if contradiction:
                    trace.append({"degree":degree,"move":"CONTRADICTION","reason":reason})
                    terminal={"mechanism":"CONTRADICTION","degree":degree,"detail":reason}
                    break
                continue
            terminal={"mechanism":"UNDECIDED","degree":degree,
                      "detail":"coefficient is not decidable by the three permitted inference steps",
                      "residual":str(expression)}
            trace.append({"degree":degree,"move":"STALL","residual":str(expression)})
            break
        if terminal is not None: break
        if depth>=top:
            active=ansatz.symbols-{x for x,_ in substitutions}
            solvable,witness=terminal_solution(residuals,ansatz.nonzero,active)
            terminal={"mechanism":"SOLVABLE" if solvable else "UNDECIDED","degree":0,
                      "detail":witness or "full walk ended without proof either way"}
            break
        depth=min(top,depth*2)
    mechanism=terminal["mechanism"]
    if mechanism=="SOLVABLE": classification="DISAGREEMENT"
    elif ((mechanism=="CONTRADICTION" and expected=="CONTRADICTION") or
          (mechanism=="STATE_KILL" and expected=="STATE_KILLED_BY_DEGREE_DROP")):
        classification="CONFIRMED"
    else:
        classification="UNDECIDED-BY-AUDIT"
        if mechanism in ("CONTRADICTION","STATE_KILL"):
            terminal["detail"]=(f"independent mechanism {mechanism} does not match "
                                f"recorded verdict {expected}; "+terminal["detail"])
    return {"identifier":name,"classification":classification,"expected":expected,
            "mechanism":mechanism,"degree":terminal.get("degree"),
            "detail":terminal.get("detail"),"residual":terminal.get("residual"),
            "ansatz":ansatz.description,
            "substitutions":[[str(a),str(b)] for a,b in substitutions],
            "parameter_constraints":[str(x) for x in residuals],"trace":trace,
            "runtime_seconds":time.perf_counter()-started}

def worker(index):
    try:
        graded=parse_graded(); _,rows=load_records()
        need(0<=index<len(rows),f"worker index {index} out of range")
        result=audit_record(index,rows[index],graded)
    except Exception as exc:
        result={"identifier":f"K{index+1:03d}","classification":"UNDECIDED-BY-AUDIT",
                "mechanism":"WORKER_ERROR","detail":f"{type(exc).__name__}: {exc}",
                "trace":[],"runtime_seconds":0.0}
    print(json.dumps(result,sort_keys=True))
    return 0

def run_child(index,timeout):
    command=[sys.executable,str(Path(__file__).resolve()),"--worker-index",str(index)]
    environment=dict(os.environ); environment["PYTHONDONTWRITEBYTECODE"]="1"
    started=time.perf_counter()
    try:
        done=subprocess.run(command,cwd=HERE,env=environment,capture_output=True,
                            text=True,timeout=timeout,check=False)
    except subprocess.TimeoutExpired:
        return {"identifier":f"K{index+1:03d}","classification":"UNDECIDED-BY-AUDIT",
                "mechanism":"TIMEOUT","detail":f"child exceeded {timeout:g}s and was terminated",
                "trace":[],"runtime_seconds":time.perf_counter()-started}
    if done.returncode:
        return {"identifier":f"K{index+1:03d}","classification":"UNDECIDED-BY-AUDIT",
                "mechanism":"WORKER_ERROR",
                "detail":f"child exit {done.returncode}; stderr={done.stderr.strip()}",
                "trace":[],"runtime_seconds":time.perf_counter()-started}
    try: return json.loads(done.stdout)
    except json.JSONDecodeError as exc:
        return {"identifier":f"K{index+1:03d}","classification":"UNDECIDED-BY-AUDIT",
                "mechanism":"WORKER_ERROR",
                "detail":f"invalid child JSON: {exc}; stdout={done.stdout[-1000:]!r}",
                "trace":[],"runtime_seconds":time.perf_counter()-started}

def print_record(result):
    print(f"{result['identifier']} {result['classification']} "
          f"mechanism={result.get('mechanism')} degree={result.get('degree')} "
          f"seconds={result.get('runtime_seconds',0.0):.3f}")
    print(f"  ansatz: {result.get('ansatz','(worker did not reconstruct ansatz)')}")
    for step in result.get("trace",[]):
        move=step.get("move")
        if move=="IDENTITY": continue
        fields=", ".join(f"{k}={v}" for k,v in step.items() if k!="move")
        print(f"  {move}: {fields}")
    if result.get("classification")!="CONFIRMED":
        print(f"  DETAIL: {result.get('detail')}")
        if result.get("residual"): print(f"  RESIDUAL: {result['residual']}")
        if result.get("substitutions"): print(f"  SUBSTITUTIONS: {result['substitutions']}")
        if result.get("parameter_constraints"):
            print(f"  PARAMETER_CONSTRAINTS: {result['parameter_constraints']}")

def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quiet",action="store_true",help="suppress per-record traces")
    parser.add_argument("--timeout",type=float,default=DEFAULT_TIMEOUT,
                        help=f"hard child timeout per record (default {DEFAULT_TIMEOUT:g}s)")
    parser.add_argument("--worker-index",type=int,help=argparse.SUPPRESS)
    args=parser.parse_args(argv)
    if args.worker_index is not None: return worker(args.worker_index)
    if args.timeout<=0:
        print("AUDIT SETUP FAILED: --timeout must be positive",file=sys.stderr); return 2
    started=time.perf_counter()
    try: graded=parse_graded(); _,rows=load_records()
    except Exception as exc:
        print(f"AUDIT SETUP FAILED: {type(exc).__name__}: {exc}",file=sys.stderr); return 2
    print("parser_and_homogeneity: PASS "
          f"(h_0..h_7; sigma term counts={graded.counts}; h_5/h_6/h_7 spot checks)")
    print("ansatz_semantics: degree-exact leading coefficients declared NONZERO")
    print(f"process_isolation: PASS (hard timeout {args.timeout:g}s per record)")
    results=[]
    for index,row in enumerate(rows):
        result=run_child(index,args.timeout)
        if result.get("identifier")==f"K{index+1:03d}":
            result["identifier"]=identifier(index,row)
        results.append(result)
        if not args.quiet: print_record(result)
    labels=("CONFIRMED","UNDECIDED-BY-AUDIT","DISAGREEMENT")
    counts={x:sum(r.get("classification")==x for r in results) for x in labels}
    undecided=[r["identifier"] for r in results if r.get("classification")==labels[1]]
    disagreements=[r["identifier"] for r in results if r.get("classification")==labels[2]]
    elapsed=time.perf_counter()-started
    print(f"FINAL CENSUS: CONFIRMED={counts[labels[0]]} "
          f"UNDECIDED-BY-AUDIT={counts[labels[1]]} DISAGREEMENT={counts[labels[2]]} "
          f"TOTAL={len(results)}")
    print("UNDECIDED IDENTIFIERS: "+(", ".join(undecided) if undecided else "none"))
    print("DISAGREEMENT IDENTIFIERS: "+(", ".join(disagreements) if disagreements else "none"))
    print(f"TOTAL RUNTIME SECONDS: {elapsed:.3f}")
    return 0 if counts["CONFIRMED"]==60 else 1

if __name__=="__main__":
    raise SystemExit(main())
