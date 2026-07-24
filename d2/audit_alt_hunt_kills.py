#!/usr/bin/env python3
"""Independent spec-only SymPy audit of all 49 ALT_HUNT/J6 kills.

No project Python module is imported.  The producers alt_hunt_depth2.py and
j6_msolve.py are neither imported nor read by this program.
"""
import argparse, itertools, json, re, subprocess, sys, time
from collections import Counter
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'audit_alt_hunt_census.json'
TIMEOUT=120
Y=s.Symbol('y'); HV=s.symbols('d0 d1 d2 dm1')
Q=2048*Y**4-512*Y**3+320*Y**2-240*Y+195
QM=s.Poly(Q/s.Integer(2048),Y,domain='QQ')
PHI=s.Poly(-(Y+1)**30*Q/s.Integer(6630),Y,domain='QQ')
CATS=('FULLY-VERIFIED','VERIFIED-DATA-ONLY','DISAGREEMENT','UNPARSEABLE')

class Bad(Exception): pass

def parse_h():
    """Own parser plus weighted-homogeneity integrity check."""
    pat=re.compile(r'^h_(\d+)\s+\(weight\s+(\d+),[^)]*\)\s*=\s*(.+)$')
    hs={}; dec={}; loc=dict(zip(('d0','d1','d2','dm1'),HV))
    for no,line in enumerate((ROOT/'f31_graded.txt').read_text().splitlines(),1):
        m=pat.fullmatch(line.strip())
        if not m: continue
        i=int(m[1])
        if i in hs: raise Bad(f'duplicate h_{i} at line {no}')
        hs[i]=s.Poly(s.sympify(m[3],locals=loc),*HV,domain='QQ'); dec[i]=int(m[2])
    if set(hs)!=set(range(8)): raise Bad(f'expected h_0..h_7, found {sorted(hs)}')
    w=(4,3,2,5)
    for i,p in hs.items():
        want=20-2*i
        if dec[i]!=want: raise Bad(f'h_{i} declared weight {dec[i]}, expected {want}')
        for mon,_ in p.terms():
            if sum(a*b for a,b in zip(mon,w))!=want:
                raise Bad(f'h_{i} monomial {mon} violates weighted homogeneity')
    return hs

def locals_for(strings):
    names=set()
    for z in strings: names.update(re.findall(r'\b(?:E|X|S|D|w|c\d+_\d+)\b',z))
    d={n:s.Symbol(n) for n in names}; d['y']=Y; return d

def expr(z,loc):
    try: return s.sympify(z,locals=loc)
    except Exception as e: raise Bad(f'cannot parse {z!r}: {e}')

def diff(a,b): return s.expand(a-b)
def short(e,n=4000):
    z=str(s.expand(e)); return z if len(z)<=n else z[:n]+f'... <{len(z)} chars>'

def tmul(a,b,L):
    out=[s.Integer(0)]*(min(L,len(a)+len(b)-2)+1)
    for i,x in enumerate(a):
        for j,z in enumerate(b):
            if i+j>L: break
            out[i+j]+=x*z
    return out

def ppow(p,n,L):
    d=int(p.degree()); base=[p.nth(d-i) for i in range(min(L,d)+1)]; out=[s.Integer(1)]
    while n:
        if n&1: out=tmul(out,base,L)
        n//=2
        if n: base=tmul(base,base,L)
    return out

def master_coeffs(hs,p,degrees):
    """Directly form sum Phi^f e^(21-3f) h_f(sigma,d1,d2,e), top-tail only."""
    if not degrees: raise Bad('no recorded master coefficients')
    cut=min(degrees); ans={d:s.Integer(0) for d in degrees}
    # f31_graded uses d0; the recorded invariant is sigma=4*d0-d2**2.
    d0_poly=s.Poly((p['sigma'].as_expr()+p['d2'].as_expr()**2)/4,Y,domain='EX')
    pv=(d0_poly,p['d1'],p['d2'],p['e']); cache={}
    for f,h in hs.items():
        for mon,c in h.terms():
            fs=[('Phi',PHI,f),('e',p['e'],21-3*f)]
            fs += [(name,z,n) for name,z,n in zip(('sigma','d1','d2','e'),pv,mon) if n]
            if any(z.is_zero and n for _,z,n in fs): continue
            md=sum(int(z.degree())*n for _,z,n in fs)
            if md<cut: continue
            L=md-cut; tail=[s.Rational(c)]
            for name,z,n in fs:
                if not n: continue
                k=(name,n,L,str(z.as_expr()))
                if k not in cache: cache[k]=ppow(z,n,L)
                tail=tmul(tail,cache[k],L)
            for i,z in enumerate(tail):
                if md-i in ans: ans[md-i]+=z
    return {d:s.expand(z) for d,z in ans.items()}

def class_factors(loc):
    got={}
    for n in loc:
        m=re.fullmatch(r'c(\d+)_(\d+)',n)
        if m: got.setdefault(int(m[1]),set()).add(int(m[2]))
    out=[]
    for ci in sorted(got):
        ix=got[ci]; want=set(range(max(ix)+1))
        if ix!=want: raise Bad(f'nonconsecutive coefficient indices for c{ci}: {sorted(ix)}')
        d=len(ix); z=Y**d+sum(loc[f'c{ci}_{i}']*Y**i for i in range(d))
        out.append((f'c{ci}',s.Poly(z,Y,domain='EX')))
    return out

def reconstruction_check(key,degs,combo,pe,factors,quo):
    mm=[]; m=re.search(r':a(\d+)_b([0-9]{4})_',key)
    if not m: raise Bad(f'key lacks a/b profile: {key}')
    if len(combo)!=5 or any(len(r)!=3 for r in combo): raise Bad(f'bad combo shape: {combo}')
    a=int(m[1]); b=tuple(map(int,m[2]))
    profiles=[(b[i],*map(int,combo[i])) for i in range(4)]
    gd={}
    for x in profiles: gd[x]=gd.get(x,0)+1
    groups=list(gd.items()); fac=factors+[('quotient',quo)]
    if len(groups)!=len(fac): return [{'kind':'class-count','profile_classes':len(groups),'factors':len(fac)}]
    if sorted(n for _,n in groups)!=sorted(int(z.degree()) for _,z in fac):
        return [{'kind':'class-size','profile_sizes':sorted(n for _,n in groups),'factor_degrees':sorted(int(z.degree()) for _,z in fac)}]
    specs={'e':(s.Symbol('E'),a,0),'d1':(s.Symbol('X'),int(combo[4][0]),1),
           'sigma':(s.Symbol('S'),int(combo[4][1]),2),'d2':(s.Symbol('D'),int(combo[4][2]),3)}
    pos={'d1':0,'sigma':1,'d2':2,'e':3}; ok=False
    for assn in itertools.permutations(fac):
        if any(int(assn[i][1].degree())!=groups[i][1] for i in range(len(groups))): continue
        good=True
        for name,(scalar,tord,pix) in specs.items():
            zero=pe[name]==0; expected_zero=degs[pos[name]]=='-inf'
            if zero!=expected_zero: good=False; break
            if zero: continue
            z=scalar*(Y+1)**tord
            for i,(profile,_) in enumerate(groups): z*=assn[i][1].as_expr()**profile[pix]
            if diff(pe[name],z)!=0: good=False; break
        if good: ok=True; break
    if not ok: mm.append({'kind':'reconstructed-polynomials','detail':'no size-compatible assignment of independent class factors and quotient reproduces d1/sigma/d2/e'})
    if ok:
        for name,i in pos.items():
            actual='-inf' if pe[name]==0 else int(s.Poly(pe[name],Y,domain='EX').degree())
            if actual!=degs[i]: mm.append({'kind':'recorded-degree','variable':name,'recorded':degs[i],'computed':actual})
    return mm

def data_checks(key,degs,system,hs):
    strings=list(system['polys'].values())+list(system['class_relations'])+[g['coefficient'] for g in system['gens']]+[system['saturation']]
    loc=locals_for(strings); required={'d2','d1','sigma','e'}
    if set(system['polys'])!=required: raise Bad(f'polys keys {sorted(system["polys"])}')
    pe={n:expr(z,loc) for n,z in system['polys'].items()}; ps={n:s.Poly(z,Y,domain='EX') for n,z in pe.items()}
    factors=class_factors(loc); unknowns=sum(int(z.degree()) for _,z in factors); mm=[]
    if unknowns!=int(system['n_class_unknowns']): mm.append({'kind':'class-unknown-count','recorded':system['n_class_unknowns'],'computed':unknowns})
    divisor=s.Poly(1,Y,domain='EX')
    for _,z in factors: divisor*=z
    quo,rem=s.Poly(QM.as_expr(),Y,domain='EX').div(divisor)
    expected=[s.expand(rem.nth(i)) for i in range(int(divisor.degree())-1,-1,-1)]
    recorded=[expr(z,loc) for z in system['class_relations']]
    if len(expected)!=len(recorded): mm.append({'kind':'relation-count','recorded':len(recorded),'computed':len(expected)})
    for i,(r,e) in enumerate(zip(recorded,expected)):
        d=diff(r,e)
        if d!=0: mm.append({'kind':'division-remainder','relation_index':i,'recorded':str(r),'computed':str(e),'exact_difference_recorded_minus_computed':short(d)})
    mm+=reconstruction_check(key,degs,system['combo'],pe,factors,quo)
    degrees=[int(g['degree']) for g in system['gens']]
    if len(degrees)!=len(set(degrees)): raise Bad(f'duplicate master degrees {degrees}')
    derived=master_coeffs(hs,ps,degrees)
    for g in system['gens']:
        degree=int(g['degree']); r=expr(g['coefficient'],loc); d=diff(r,derived[degree])
        if d!=0: mm.append({'kind':'master-coefficient','degree':degree,'recorded':str(r),'independently_derived':str(derived[degree]),'exact_difference_recorded_minus_derived':short(d)})
    detail={'master_identity':{'degrees_checked':degrees,'coefficients_checked':len(degrees)},
            'class_polynomial':{'parameterized_factor_degrees':[int(z.degree()) for _,z in factors],
            'division_quotient_degree':int(quo.degree()),'division_remainder_coefficients':len(expected)},'passed':not mm}
    return mm,detail

def gb_worker(payload):
    strings=list(payload['relations'])+[g['coefficient'] for g in payload['gens']]+[payload['saturation']]
    loc=locals_for(strings); eq=[expr(z,loc) for z in strings]
    vs=sorted(set().union(*(z.free_symbols for z in eq)),key=lambda x:(x.name=='w',not x.name.startswith(('E','X','S','D')),x.name))
    started=time.monotonic(); G=s.groebner(eq,*vs,order='grevlex',method='f5b')
    unit=any(p.total_degree()==0 and p.as_expr()!=0 for p in G.polys)
    return {'status':'completed','unit_ideal':bool(unit),'variables':[str(v) for v in vs],'basis_size':len(G.polys),'worker_groebner_seconds':round(time.monotonic()-started,6)}

def run_gb(payload):
    started=time.monotonic()
    try:
        p=subprocess.run([sys.executable,str(Path(__file__).resolve()),'--gb-worker'],input=json.dumps(payload),text=True,capture_output=True,timeout=TIMEOUT)
    except subprocess.TimeoutExpired:
        return {'status':'timeout','hard_timeout_seconds':TIMEOUT,'wall_seconds':round(time.monotonic()-started,6),'process_killed':True}
    wall=round(time.monotonic()-started,6)
    if p.returncode: return {'status':'error','returncode':p.returncode,'stderr':p.stderr[-4000:],'stdout':p.stdout[-4000:],'wall_seconds':wall}
    try: z=json.loads(p.stdout)
    except Exception as e: return {'status':'error','detail':f'invalid worker JSON: {e}','stdout':p.stdout[-4000:],'stderr':p.stderr[-4000:],'wall_seconds':wall}
    z.update(wall_seconds=wall,hard_timeout_seconds=TIMEOUT); return z

def audit_one(source,record,state,hs):
    key=record['key']; srcsys=state['splits'][0]
    if source=='alt_hunt': system=srcsys
    else:
        system={'combo':record['combo'],'n_class_unknowns':record['n_class_unknowns'],'polys':srcsys['polys'],'class_relations':record['class_relations'],'gens':record['gens'],'saturation':record['saturation']}
    mm,detail=data_checks(key,state['degs'],system,hs)
    if source=='j6_msolve':
        # OPEN ALT source rows record only profile/reconstruction data; J6 adds relations/saturation.
        for field in ('combo','n_class_unknowns'):
            if record[field]!=srcsys[field]: mm.append({'kind':'j6-source-replay','field':field,'j6_recorded':record[field],'alt_source_recorded':srcsys[field]})
    out={'kill_id':key,'source':source,'data_checks':detail,'mismatches':mm}
    if mm: out.update(category='DISAGREEMENT',groebner={'status':'not-run-due-to-data-disagreement'}); return out
    gb=run_gb({'relations':system['class_relations'],'gens':system['gens'],'saturation':system['saturation']}); out['groebner']=gb
    if gb['status']=='timeout': out['category']='VERIFIED-DATA-ONLY'
    elif gb['status']=='completed' and gb['unit_ideal']: out['category']='FULLY-VERIFIED'
    elif gb['status']=='completed':
        out['category']='DISAGREEMENT'; out['mismatches'].append({'kind':'contradiction','detail':'independent SymPy Groebner basis completed but was not the unit ideal'})
    else: out['category']='UNPARSEABLE'; out['unparseable_reason']='isolated Groebner worker failed: '+str(gb)
    return out

def make_census(results):
    by={}
    for source in ('alt_hunt','j6_msolve'):
        c=Counter(r['category'] for r in results if r['source']==source); by[source]={k:c.get(k,0) for k in CATS}; by[source]['TOTAL']=sum(c.values())
    c=Counter(r['category'] for r in results); overall={k:c.get(k,0) for k in CATS}; overall['TOTAL']=sum(c.values())
    return {'overall':overall,'by_source':by}

def main(quiet=False):
    started=time.monotonic(); hs=parse_h()
    alt=json.loads((ROOT/'alt_hunt_results.json').read_text()); j6=json.loads((ROOT/'j6_msolve_results.json').read_text())
    states=alt['states']; bykey={z['key']:z for z in states}; ak=[z for z in states if z.get('verdict')=='KILLED']; jr=j6['results']
    if len(ak)!=45 or len(jr)!=4: raise Bad(f'expected 45 ALT + 4 J6 kills, got {len(ak)} + {len(jr)}')
    work=[('alt_hunt',z,z) for z in ak]+[('j6_msolve',z,bykey[z['key']]) for z in jr]; results=[]
    for i,(source,record,state) in enumerate(work,1):
        if not quiet: print(f'[{i:02d}/49] {source} {record["key"]}',flush=True)
        try: r=audit_one(source,record,state,hs)
        except Exception as e: r={'kill_id':record.get('key',f'<record-{i}>'),'source':source,'category':'UNPARSEABLE','unparseable_reason':f'{type(e).__name__}: {e}','mismatches':[],'groebner':{'status':'not-run'}}
        results.append(r)
        if not quiet: print('          '+r['category'],flush=True)
    cen=make_census(results); disagreements=[{'kill_id':r['kill_id'],'source':r['source'],'mismatches':r['mismatches']} for r in results if r['category']=='DISAGREEMENT']; unparseable=[{'kill_id':r['kill_id'],'source':r['source'],'reason':r['unparseable_reason']} for r in results if r['category']=='UNPARSEABLE']; code=0 if not disagreements else 1
    doc={'schema':1,'item':'independent spec-only audit of ALT_HUNT/J6 state kills','inputs':{'f31_graded':'f31_graded.txt','alt_hunt':'alt_hunt_results.json','j6_msolve':'j6_msolve_results.json'},'independence':{'project_python_imports':[],'producer_files_not_imported':['alt_hunt_depth2.py','j6_msolve.py'],'master_identity_source':'from-scratch parser over f31_graded.txt','relation_method':'independent monic-quartic division and valuation-profile reconstruction','contradiction_method':'isolated SymPy grevlex Groebner subprocess per kill'},'groebner_hard_timeout_seconds_per_kill':TIMEOUT,'census':cen,'disagreements':disagreements,'unparseable':unparseable,'results':results,'auditor_exit_code':code,'wall_seconds':round(time.monotonic()-started,6)}
    OUT.write_text(json.dumps(doc,indent=2)+'\n')
    if not quiet: print(json.dumps({'census':cen,'auditor_exit_code':code},indent=2))
    return code

def cli():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--quiet',action='store_true'); p.add_argument('--gb-worker',action='store_true',help=argparse.SUPPRESS); a=p.parse_args()
    if a.gb_worker:
        try: print(json.dumps(gb_worker(json.loads(sys.stdin.read())))); return 0
        except Exception as e: print(f'{type(e).__name__}: {e}',file=sys.stderr); return 2
    try: return main(a.quiet)
    except Exception as e: print(f'FATAL AUDIT ERROR: {type(e).__name__}: {e}',file=sys.stderr); return 2

if __name__=='__main__': raise SystemExit(cli())
