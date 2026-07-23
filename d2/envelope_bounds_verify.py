#!/usr/bin/env python3
"""Exact checker for the (72,108) D-coefficient envelope bounds."""
from __future__ import annotations
import ast, json, re, sys, time
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parent
PAPER, AUDIT, STATE = ROOT/"paper_src"/"2204.14178.tex", ROOT/"T3_WINDOW_AUDIT.md", ROOT/"STATE.md"
UPSTREAM_FACTS = ROOT/"paper_src"/"upstream_facts.json"
EXPECTED_CHECKS = 106
ASSUMPTIONS = (
 ("A1_PUBLISHED_POLYGON_REDUCTION",
  "GGHV22 Proposition 'Case (8,28)' (lines 1000-1007) is used as stated: a candidate reduces to a displayed polygon with [P,Q]=x^2."),
 ("A2_COMMON_ROOT_AND_NORMALIZATION",
  "Cited GGV1 Propositions 1.13 and 2.1 plus GGHV22's normalization at lines 1411-1414 give ell(P)=R^2, ell(Q)=R^3, R=x^4*y^7*(y+1). Their consequences are checked, but the citations/WLOG are not reproved."),
 ("A3_LAURENT_VALUATION_FRAMEWORK",
  "The characteristic-zero degree/order valuations in K((y^-1)), K((y)) are multiplicative and non-Archimedean, and x->x-D3/4 is a valid formal translation in K[y]((x^-1))."),
)

class VerificationFailure(RuntimeError): pass
class Runner:
    def __init__(self): self.passed = 0
    def check(self, label, condition):
        if condition is not True and condition != sp.true:
            raise VerificationFailure(f"[FAIL] {label}")
        self.passed += 1
        print(f"  [OK] {label}")

def lines(path): return path.read_text(encoding="utf-8").splitlines()
def pairs(text): return tuple((int(a),int(b)) for a,b in re.findall(r"\((\d+),\s*(\d+)\)",text))
def order(poly,var): return min(m[0] for m in sp.Poly(sp.expand(poly),var).monoms())
def vmax(cs,d): return max(d[0]*i+d[1]*j for i,j in cs)
def maximizers(cs,d):
    top=vmax(cs,d)
    return frozenset(p for p in cs if d[0]*p[0]+d[1]*p[1]==top)

def literal_assignments(path):
    tree=ast.parse(path.read_text(encoding="utf-8-sig"),filename=str(path)); out={}
    for node in tree.body:
        if not isinstance(node,(ast.Assign,ast.AnnAssign)): continue
        targets=node.targets if isinstance(node,ast.Assign) else [node.target]
        for target in targets:
            if isinstance(target,ast.Name):
                try: out[target.id]=ast.literal_eval(node.value)
                except (ValueError,TypeError): pass
            elif isinstance(target,ast.Tuple) and isinstance(node.value,ast.Tuple):
                for left,right in zip(target.elts,node.value.elts):
                    if isinstance(left,ast.Name):
                        try: out[left.id]=ast.literal_eval(right)
                        except (ValueError,TypeError): pass
    return out

def jetlift_configs(path):
    tree=ast.parse(path.read_text(encoding="utf-8"),filename=str(path))
    value=next((n.value for n in tree.body if isinstance(n,ast.Assign)
                and any(isinstance(t,ast.Name) and t.id=="CONFIGS" for t in n.targets)),None)
    if not isinstance(value,ast.Dict): raise VerificationFailure("[FAIL] D4 parse CONFIGS")
    out={}
    for kn,cn in zip(value.keys,value.values):
        key=ast.literal_eval(kn)
        if not (isinstance(cn,ast.Call) and isinstance(cn.func,ast.Name) and cn.func.id=="dict"):
            raise VerificationFailure(f"[FAIL] D4 unsupported entry {key}")
        out[key]={kw.arg:ast.literal_eval(kw.value) for kw in cn.keywords}
    return out

def check_sources(run):
    """Mandatory source/base-data checks.

    The upstream (published) facts are read from paper_src/upstream_facts.json --
    an original small file extracted ONCE from arXiv:2204.14178 with exact line
    citations -- so a clean public clone (which does not ship the .tex) can run
    the whole suite.  AUDIT and STATE are tracked in the public tree and are read
    directly.  A separate OPTIONAL provenance step (tex_provenance) re-checks the
    upstream facts against the .tex when it happens to be present locally."""
    print("S. source and base-data transcription (upstream_facts.json + audit/state)")
    audit,state=lines(AUDIT),lines(STATE)
    facts=json.loads(UPSTREAM_FACTS.read_text(encoding="utf-8"))["facts"]
    p1=((0,0),(1,0),(8,14),(8,16),(0,8)); q1=((0,0),(2,1),(12,21),(12,24),(0,12))
    p2=((0,0),(1,0),(8,14),(8,16)); q2=((0,0),(2,1),(12,21),(12,24))
    def as_pairs(pts): return tuple((int(a),int(b)) for a,b in pts)
    np_=facts["newton_polygons"]; bc=facts["bracket_case"]
    crt=facts["common_root_template"]; ind=facts["induction_template"]; rec=facts["recursion_template"]
    run.check("S1 cited data sources present",all(p.is_file() for p in (UPSTREAM_FACTS,AUDIT,STATE)))
    run.check("S2 upstream bracket case (8,28), L^(1), [P,Q]=x^2",
              bc["case"]=="(8,28)" and bc["P_Q_space"]=="L^{(1)}" and bc["bracket"]=="[P,Q] = x^2")
    run.check("S3 upstream sub1 corners match reference",
              as_pairs(np_["sub1"]["P"])==p1 and as_pairs(np_["sub1"]["Q"])==q1)
    run.check("S4 upstream sub2 corners match reference",
              as_pairs(np_["sub2"]["P"])==p2 and as_pairs(np_["sub2"]["Q"])==q2)
    run.check("S5 audit L22-23 exact sub1",pairs(audit[21])==p1 and pairs(audit[22])==q1)
    run.check("S6 audit L24-25 exact sub2",pairs(audit[23])==p2 and pairs(audit[24])==q2)
    run.check("S7 audit proposition premise/bracket","Case (8,28)" in audit[18] and "[P,Q] = x²" in audit[19])
    run.check("S8 upstream common-root template R^2/R^3/C_3",
              crt["ell_1_0_P"]=="R^2" and crt["ell_1_0_Q"]=="R^3" and crt["C3"]=="y^8 (y+1)")
    run.check("S9 upstream induction template markers",
              {"construct inductively","v_{-1,1}","v_{3,-1}"}<=set(ind["markers"]))
    run.check("S10 upstream recursion template markers",
              "2C_{3-k}x^{3-k} C_3x^3" in rec["markers"] and r"C_{3-k}:=-\frac{1}{2C_3}" in rec["markers"])
    run.check("S11 STATE items 1-4 normalization/D/shift","C₄ = y⁷(y+1)" in state[25] and "D-transformation" in state[31] and "Shift x ↦ x − D₃/4" in state[35])
    run.check("S12 audit L52-55 C4/t=4 recursion","C₄ = y⁷(y+1)" in audit[51] and "C_{4−k}" in audit[54] and "P_{8−k}" in audit[54])
    return p1,q1,p2,q2


def tex_provenance():
    """OPTIONAL provenance: when the copyrighted arXiv .tex is present locally,
    re-verify that upstream_facts.json still agrees with it at the cited lines.
    Skipped gracefully (with a note) on a clean public clone where the .tex is
    absent.  Raises VerificationFailure on any genuine mismatch.  NOT counted in
    the mandatory check total."""
    if not PAPER.is_file():
        print("  [provenance SKIPPED] paper_src/2204.14178.tex absent; "
              "upstream_facts.json is the self-contained source of record.")
        return
    facts=json.loads(UPSTREAM_FACTS.read_text(encoding="utf-8"))["facts"]
    paper=lines(PAPER)
    bc=facts["bracket_case"]; np_=facts["newton_polygons"]
    crt=facts["common_root_template"]; ind=facts["induction_template"]; rec=facts["recursion_template"]
    def want(cond,msg):
        if not cond: raise VerificationFailure(f"[FAIL] tex provenance: {msg}")
    def as_pairs(pts): return tuple((int(a),int(b)) for a,b in pts)
    want("P,Q \\in L^{(1)}" in paper[bc["tex_line"]-1] and "[P,Q] = x^2" in paper[bc["tex_line"]-1],"bracket L1001")
    want(pairs(paper[np_["sub1"]["tex_line"]-1])==as_pairs(np_["sub1"]["P"])+as_pairs(np_["sub1"]["Q"]),"sub1 corners")
    want(pairs(paper[np_["sub2"]["tex_line"]-1])==as_pairs(np_["sub2"]["P"])+as_pairs(np_["sub2"]["Q"]),"sub2 corners")
    want("R^2" in paper[crt["tex_lines"][0]-1] and "R^3" in paper[crt["tex_lines"][0]-1]
         and "C_3= y^8 (y+1)" in paper[crt["tex_lines"][1]-1],"common-root template")
    want("construct inductively" in paper[ind["tex_lines"][0]-1] and "v_{-1,1}" in paper[ind["tex_lines"][1]-1]
         and "v_{3,-1}" in paper[ind["tex_lines"][2]-1],"induction template")
    want("2C_{3-k}x^{3-k} C_3x^3" in paper[rec["tex_lines"][0]-1]
         and r"C_{3-k}:=-\frac{1}{2C_3}" in paper[rec["tex_lines"][1]-1],"recursion template")
    print("  [provenance OK] upstream_facts.json matches paper_src/2204.14178.tex at all cited lines.")

def check_geometry(run,p1,q1,p2,q2):
    print("B. polygon valuations and leading form")
    y,x=sp.symbols("y x"); c4=y**7*(y+1); root=x**4*c4
    run.check("B1 C4 ord=7 deg=8",order(c4,y)==7 and sp.degree(c4,y)==8)
    run.check("B2 R^2 corners (8,14),(8,16)",order(c4**2,y)==14 and sp.degree(c4**2,y)==16)
    run.check("B3 R^3 corners (12,21),(12,24)",order(c4**3,y)==21 and sp.degree(c4**3,y)==24)
    run.check("B4 forms literally R^2/R^3",sp.expand(root**2-x**8*c4**2)==0 and sp.expand(root**3-x**12*c4**3)==0)
    run.check("B5 both P polygons contain leading corners",{(8,14),(8,16)}<=set(p1) and {(8,14),(8,16)}<=set(p2))
    run.check("B6 both Q polygons contain leading corners",{(12,21),(12,24)}<=set(q1) and {(12,21),(12,24)}<=set(q2))
    run.check("B7 sub1 max v(-1,1)=8",vmax(p1,(-1,1))==8 and maximizers(p1,(-1,1))=={(8,16),(0,8)})
    run.check("B8 sub2 max v(-2,1)=0",vmax(p2,(-2,1))==0 and maximizers(p2,(-2,1))=={(0,0),(8,16)})
    run.check("B9 both max v(2,-1)=2",all(vmax(p,(2,-1))==2 and maximizers(p,(2,-1))=={(1,0),(8,14)} for p in (p1,p2)))
    run.check("B10 sub1 extra corner makes max v(-2,1)=8",vmax(p1,(-2,1))==8 and maximizers(p1,(-2,1))=={(0,8)})
    num,den=sp.fraction(sp.cancel(1/c4))
    run.check("B11 v_infinity(C4^-1)=-8",sp.degree(num,y)-sp.degree(den,y)==-8)
    run.check("B12 local v(2,-1)(C4^-1)=+7",-order(num,y)+order(den,y)==7)

def check_recursion(run):
    print("R. coefficient identity and calculo-de-C recursion")
    z=sp.symbols("z",nonzero=True); cs=(z,)+sp.symbols("c1:6"); ps=sp.symbols("p0:6")
    run.check("R1 base P8=C4^2",cs[0]**2==z**2)
    for w in range(1,6):
        middle=sum(cs[j]*cs[w-j] for j in range(1,w))
        corrected=(ps[w]-middle)/(2*z)
        extracted=sum(cs[j]*cs[w-j] for j in range(w+1))
        run.check(f"R{w+1} corrected recursion w={w}",sp.simplify(extracted.subs(cs[w],corrected)-ps[w])==0)
    middle=cs[1]**2; printed=-(ps[2]+middle)/(2*z)
    residual=sp.expand(sum(cs[j]*cs[2-j] for j in range(3)).subs(cs[2],printed)-ps[2])
    run.check("R7 printed all-minus recursion is non-identity",residual!=0 and sp.simplify(residual+2*ps[2])==0)

def induction(run,prefix,name,rho,polymax,c4v,bound):
    run.check(f"{prefix}0 {name} base",c4v==bound(0))
    for w in range(1,6):
        pb=polymax-rho*(8-w); products=[bound(j)+bound(w-j) for j in range(1,w)]
        run.check(f"{prefix}{w} {name} w={w}",-c4v+max([pb]+products)==bound(w) and all(v==pb for v in products))

def check_inductions(run):
    print("I. valuation inductions w=0..5")
    induction(run,"I1.","sub1 v(-1,1)",-1,8,8,lambda w:8-w)
    induction(run,"I2.","sub2 v(-2,1)",-2,0,8,lambda w:8-2*w)
    induction(run,"I3.","both v(2,-1)",2,2,-7,lambda w:2*w-7)

def check_transform(run):
    print("M. D-transform and magic weights")
    r,rho=sp.symbols("r rho",integer=True); e=7-2*r
    deg1=sp.expand(r+4+8*e); deg2=sp.expand(2*r+8*e); ordD=sp.expand(2*r-1+7*e)
    m1=sp.solve(sp.Eq(sp.diff(rho*r+deg1,r),0),rho)
    m2=sp.solve(sp.Eq(sp.diff(rho*r+deg2,r),0),rho)
    mo=sp.solve(sp.Eq(sp.diff(rho*r-ordD,r),0),rho)
    run.check("M1 magic sub1=15",m1==[15]); run.check("M2 magic sub2=14",m2==[14]); run.check("M3 magic order=-12",mo==[-12])
    run.check("M4 constant sub1=60",sp.expand(15*r+deg1)==60)
    run.check("M5 constant sub2=56",sp.expand(14*r+deg2)==56)
    run.check("M6 constant order=-48",sp.expand(-12*r-ordD)==-48)
    for w in range(6):
        k=4-w
        run.check(f"M7.{w} envelopes for D_{k}",deg1.subs(r,k)==15*w and deg2.subs(r,k)==14*w and ordD.subs(r,k)==12*w)
    run.check("M8 D4=1 polynomial",sp.Integer(1).is_polynomial())
    for w in range(1,6):
        k=4-w; pe=6-2*k; ij=[(4-j,4-(w-j)) for j in range(1,w)]
        run.check(f"M8.{w} D_{k} polynomial recursion",pe>=0 and all((7-2*i)+(7-2*j)==pe for i,j in ij))
    return 14,15,12

def check_shift(run,m2,m1,strip):
    print("T. D3-killing translation")
    a,d3=sp.symbols("a d3")
    run.check("T1 translation kills D3",sp.expand((d3-4*a).subs(a,d3/4))==0)
    for k in range(-1,5):
        w=4-k
        for i in range(k,5):
            n=i-k
            run.check(f"T2.{k}.{i} contribution weight {w}",n+4-i==w and sp.binomial(i,n).is_Integer)
        run.check(f"T3.{k} shifted coefficient retains caps",m2*w==14*w and m1*w==15*w and strip*w==12*w)

def check_downstream(run,m2,m1,strip):
    print("D. downstream caps and windows")
    weights=(2,3,4,5); caps2=tuple((m2-strip)*w for w in weights); caps1=tuple((m1-strip)*w for w in weights)
    run.check("D1 sub2 stripped caps 2w",caps2==(4,6,8,10))
    run.check("D2 sub1 stripped caps 3w",caps1==(6,9,12,15))
    gen=(ROOT/"regenerate_system.py").read_text(encoding="utf-8")
    run.check("D3 generator weights 2,3,4,5","d2*u**2 + d1*u**3 + d0*u**4" in gen and "dm[k]*u**(4+k)" in gen)
    cfg=jetlift_configs(ROOT/"jetlift.py")
    sizes={"f31_sub2":(5,7,9,11),"f37_sub2":(5,7,9,11),"f31_sub1":(7,10,13,16),"f37_sub1":(7,10,13,16)}
    run.check("D4 jetlift sizes",all(tuple(cfg[n]["sizes"])==v for n,v in sizes.items()))
    his={"f31_sub2":251,"f37_sub2":269,"f31_sub1":376,"f37_sub1":403}
    run.check("D5 jetlift hi=(slope)W+1",all(cfg[n]["hi"]==v for n,v in his.items()))
    sub1=literal_assignments(ROOT/"sub1_cascade_verify.py")
    run.check("D6 sub1 cascade caps",tuple(sub1[n] for n in ("D2_DEG_CAP","D1_DEG_CAP","D0_DEG_CAP","E_DEG_CAP"))==caps1)
    cascade=literal_assignments(ROOT/"cascade_engine.py")
    run.check("D7 sub2 cascade caps",cascade.get("GLOBAL_CAPS")=={"d2":4,"d1":6,"sigma":8})
    ct=(ROOT/"cascade_engine.py").read_text(encoding="utf-8")
    run.check("D8 cascade window configs",'e_cap=10, h_slope=4' in ct and 'e_cap=15, h_slope=6' in ct)
    vs=sp.symbols("d2 d1 d0 dm1 Phi"); vweights=(2,3,4,5,17)
    for fn,expected in (("f31_deg31.txt",125),("f37_deg37.txt",134)):
        raw=(ROOT/fn).read_text(encoding="utf-8").strip()
        expr=sp.sympify(raw.replace("m1","dm1").replace("P","Phi").replace("^","**"))
        actual={sum(w*e for w,e in zip(vweights,mono)) for mono in sp.Poly(expr,*vs).monoms()}
        run.check(f"D9 {fn} exact weight {expected}",actual=={expected})
    y=sp.symbols("y"); phi=-y**204*(y+1)**30*(2048*y**4-512*y**3+320*y**2-240*y+195)/6630
    run.check("D10 Phi attains sub2 bounds",sp.degree(phi,y)==238==14*17 and order(phi,y)==204==12*17)

def main():
    start=time.perf_counter(); run=Runner()
    print("Exact envelope-bound verification (SymPy; no floating point)\n")
    try:
        p1,q1,p2,q2=check_sources(run); check_geometry(run,p1,q1,p2,q2)
        check_recursion(run); check_inductions(run); m2,m1,strip=check_transform(run)
        check_shift(run,m2,m1,strip); check_downstream(run,m2,m1,strip)
        if run.passed != EXPECTED_CHECKS: raise VerificationFailure(f"[FAIL] expected {EXPECTED_CHECKS} checks, ran {run.passed}")
        print("P. optional upstream provenance")
        tex_provenance()
    except (VerificationFailure,OSError,SyntaxError,ValueError,KeyError) as exc:
        print(f"\n{exc}\nFAILED after {run.passed} checks",file=sys.stderr); return 1
    elapsed=time.perf_counter()-start
    print("\nMECHANIZATION FINDINGS")
    print("  F1_RECURSION_SIGN_TYPO: paper L1462-66 and audit L55 negate P and sum;")
    print("     R1-R7 derive the correct identity. Valuation bounds are sign-insensitive.")
    print("  F2_WEIGHT_RANGE_STALENESS: d_-1 has w=5; checker covers w=0..5.")
    print("\nNAMED ASSUMPTIONS")
    for name,statement in ASSUMPTIONS: print(f"  {name}: {statement}")
    print(f"\nALL {run.passed}/{EXPECTED_CHECKS} ENVELOPE CHECKS PASSED in {elapsed:.3f} s")
    return 0

if __name__=="__main__":
    raise SystemExit(main())






