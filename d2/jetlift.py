"""
jetlift.py — validated jet-lifting feasibility harness for the (72,108) endgame
================================================================================
Decides (numerically) whether the window system  f(d2,d1,d0,dm1,Phi) === 0  in
K[y] admits solutions, for f in {f31, f37} and Newton-polygon subcase in {2,1}.
See HANDOFF.md and STATE.md for full context. This file is the version that
passed the positive control (converges to FD floor when solutions exist).

Usage:
  python3 jetlift.py control f31_sub2          # must reach ~1e-6 (validates tool)
  python3 jetlift.py stats   f31_sub2 600      # N seconds of lifting runs
  python3 jetlift.py polish  f31_sub2          # central-diff polish of best saved runs

HARD-WON NUMERICAL RULES (violate these and you will reproduce our dead ends):
  * ALL tolerances relative to term scale (scales hit 1e24 at |e0|~7 via 25th powers)
  * per-slice normalization M_s from the absolute-value evaluator (higher slices
    carry ~1e6 multinomial factors; global normalization is wrong)
  * null bases of the gradient row via SVD, never QR of a projector
  * balanced gauge for Psi: divide by sqrt(max|coeff| * |Psi(0)|)
  * scipy 'lm' ignores its budget here; use the custom damped Gauss-Newton below
  * forward-diff floor ~1e-6, central-diff floor ~1e-10: only central-diff
    plateaus count as evidence of true positive minima
"""
import sys, time, pickle
import numpy as np
import sympy as sp
from sympy import symbols, Poly

# ---------------------------------------------------------------- configuration
# sizes = coefficient-window lengths for (d2,d1,d0,dm1); hi = number of G slices.
# Subcase 2 bounds are harness-verified (deg<=14w, ord>=12w). Subcase 1 bounds
# (deg<=15w from the (0,8)/(0,12) corners) were derived hastily — AUDIT (task T3)
# before trusting any subcase-1 conclusion.
CONFIGS = {
  'f31_sub2': dict(factor='f31_deg31.txt', W=125, sizes=[5,7,9,11],  hi=251, N=512),
  'f37_sub2': dict(factor='f37_deg37.txt', W=134, sizes=[5,7,9,11],  hi=269, N=512),
  'f31_sub1': dict(factor='f31_deg31.txt', W=125, sizes=[7,10,13,16], hi=376, N=1024),
  'f37_sub1': dict(factor='f37_deg37.txt', W=134, sizes=[7,10,13,16], hi=403, N=1024),
}

import os
_HERE = os.path.dirname(os.path.abspath(__file__))

def setup(cfgname):
    cfg = CONFIGS[cfgname]
    y = symbols('y'); V = symbols('d2 d1 d0 dm1 Phi')
    # resolve the factor file relative to this script, so the harness runs from
    # any working directory (needed once the D2 files live in a subdirectory).
    s = open(os.path.join(_HERE, cfg['factor'])).read().strip()
    s = s.replace('m1','dm1').replace('P','Phi').replace('^','**')
    fe = sp.sympify(s); f = Poly(fe, *V)
    E5 = np.array(f.monoms()); CF = np.array([float(c) for c in f.coeffs()])
    st = dict(cfg); st.update(V=V, fe=fe, E5=E5, CF=CF, ACF=np.abs(CF),
                              mx=E5.max(axis=0))
    st['grads'] = [sp.lambdify(V, sp.diff(fe, v), 'numpy') for v in V[:4]]
    st['ecoef_f'] = [sp.lambdify((V[0],V[1],V[2],V[4]), c, 'numpy')
                     for c in Poly(fe, V[3]).all_coeffs()]
    # Phi = f1 * C4^28 exactly; Psi = Phi/y^204, rescaled y -> y/4, balanced gauge
    Psi = sp.expand((-(y+1)**30*(2048*y**4-512*y**3+320*y**2-240*y+195)/6630)
                    .subs(y, y/4))
    co = np.array([float(Psi.coeff(y,k)) for k in range(35)])
    gauge = np.sqrt(np.abs(co).max()*abs(co[0]))
    N = st['N']
    PsiV = np.zeros(N,complex); PsiV[:35] = co/gauge
    st['Fpsi'] = np.fft.fft(PsiV); st['Psi0'] = float(PsiV[0].real)
    APsiV = np.zeros(N); APsiV[:35] = np.abs(co)/gauge
    st['FApsi'] = np.fft.fft(APsiV.astype(complex))
    sizes = st['sizes']
    st['avail'] = {k:[i for i,sz in enumerate(sizes) if k <= sz-1]
                   for k in range(1, max(sizes))}
    st['maxslice'] = max(sizes)-1          # last slice with new unknowns
    st['NDOF'] = sum(len(st['avail'][k])-1 for k in st['avail'])
    # JETLIFT_SEED lets independent stats chunks explore different bases
    # (the fixed default would make every fresh process repeat the same runs).
    st['seed'] = int(os.environ.get('JETLIFT_SEED', 2026))
    st['rng'] = np.random.default_rng(st['seed'])
    return st

def core(st, vecs, Fp, coefs):
    N, mx, E5 = st['N'], st['mx'], st['E5']
    Fs = [np.fft.fft(np.concatenate([v, np.zeros(N-len(v))])) for v in vecs]+[Fp]
    P = []
    for i in range(5):
        T = np.empty((mx[i]+1, N), complex); T[0] = 1
        for e in range(mx[i]): T[e+1] = T[e]*Fs[i]
        P.append(T)
    prod = coefs[:,None]*P[0][E5[:,0]]*P[1][E5[:,1]]*P[2][E5[:,2]]\
                        *P[3][E5[:,3]]*P[4][E5[:,4]]
    return np.fft.ifft(prod.sum(axis=0))[:st['hi']]

def Geval(st, vecs): return core(st, vecs, st['Fpsi'], st['CF'])
def Mslice(st, vecs):
    av = [np.abs(v).astype(complex) for v in vecs]
    return np.abs(core(st, av, st['FApsi'], st['ACF'].astype(complex))).real + 1.0

def lift(st, base, nullco, g):
    vecs = [np.zeros(sz, complex) for sz in st['sizes']]
    for i in range(4): vecs[i][0] = base[i]
    ptr = 0
    for k in range(1, st['maxslice']+1):
        av = st['avail'][k]; gv = g[av]
        Gs = Geval(st, vecs)[k]
        if not np.isfinite(Gs): return None
        part = -Gs*np.conj(gv)/(np.abs(gv)**2).sum()
        _,_,Vh = np.linalg.svd(gv.reshape(1,-1))
        nb = Vh.conj().T[:,1:]
        newj = part + nb @ nullco[ptr:ptr+len(av)-1]; ptr += len(av)-1
        for j,i in enumerate(av): vecs[i][k] = newj[j]
        if max(np.abs(v).max() for v in vecs) > 25: return None
    return vecs

def make_base(st):
    rng, ecf, Psi0 = st['rng'], st['ecoef_f'], st['Psi0']
    while True:
        a0,b0,c0 = rng.normal(0,0.6,3) + 1j*rng.normal(0,0.6,3)
        ec = np.array([complex(f(complex(a0),complex(b0),complex(c0),Psi0))
                       for f in ecf])
        if not np.all(np.isfinite(ec)): continue
        deg = len(ec)-1
        der = ec[:-1]*np.arange(deg, 0, -1)
        for e0 in [r for r in np.roots(ec) if 0.02 < abs(r) < 4.0]:
            for _ in range(3):
                dv = np.polyval(der, e0)
                if abs(dv) < 1e-14: break
                e0 -= np.polyval(ec, e0)/dv
            base = np.array([a0,b0,c0,e0], complex)
            g = np.array([gr(*base, Psi0) for gr in st['grads']], complex)
            Tm = float((st['ACF']*np.prod(np.abs(np.append(base,Psi0))**st['E5'],
                                          axis=1)).sum())
            if abs(np.polyval(ec, e0)) < 1e-10*Tm and np.linalg.norm(g) > 1e-8*Tm:
                return base, g

def resid(st, p, base, g, hi=None):
    hi = hi or st['hi']
    vecs = lift(st, base, p[0::2]+1j*p[1::2], g)
    if vecs is None: return None
    G = Geval(st, vecs); M = Mslice(st, vecs)
    if not np.all(np.isfinite(G)): return None
    lo = st['maxslice']+1
    r = G[lo:hi]/M[lo:hi]
    return np.concatenate([r.real, r.imag])

def gauss_newton(st, base, g, p=None, iters=25, h=1e-6, central=False, hi=None):
    rng = st['rng']
    if p is None:
        for _ in range(6):
            p = rng.normal(0, 0.2, 2*st['NDOF'])
            r = resid(st, p, base, g, hi)
            if r is not None: break
        else: return None, None, []
    else:
        r = resid(st, p, base, g, hi)
        if r is None: return None, p, []
    lam = 1e-4; traj = [np.linalg.norm(r)]
    for it in range(iters):
        J = np.empty((len(r), len(p))); bad = False
        for j in range(len(p)):
            pj = p.copy(); pj[j] += h
            rj = resid(st, pj, base, g, hi)
            if rj is None: bad = True; break
            if central:
                pk = p.copy(); pk[j] -= h
                rk = resid(st, pk, base, g, hi)
                if rk is None: bad = True; break
                J[:,j] = (rj-rk)/(2*h)
            else:
                J[:,j] = (rj-r)/h
        if bad: break
        ok = False
        for _ in range(9):
            dp = np.linalg.lstsq(np.vstack([J, np.sqrt(lam)*np.eye(len(p))]),
                                 np.concatenate([-r, np.zeros(len(p))]),
                                 rcond=None)[0]
            rn = resid(st, p+dp, base, g, hi)
            if rn is not None and np.linalg.norm(rn) < np.linalg.norm(r):
                p, r = p+dp, rn; lam = max(lam/3, 1e-10); ok = True; break
            lam *= 5
        traj.append(np.linalg.norm(r))
        if not ok: break
    return float(np.linalg.norm(r)), p, traj

def cmd_control(st, name):
    lo = st['maxslice']+1
    hi = lo + st['NDOF'] + 2 - 3        # fewer conditions than dof -> must solve
    print(f"[control] slices y^{lo}..y^{hi-1} ({hi-lo} conds vs {st['NDOF']+3} dof)")
    # A given base can diverge on every random start (lift rejected); those are
    # not failures of the tool, so skip them and keep sampling bases until we
    # have 4 completed optimizations (or run out of attempts).
    got = 0; attempts = 0
    while got < 4 and attempts < 200:
        attempts += 1
        base, g = make_base(st)
        v,_,_ = gauss_newton(st, base, g, iters=30, hi=hi)
        if v is None: continue
        got += 1
        print(f"  control {got}: {v:.3e}"+("  OK (<=1e-5)" if v<1e-5 else "  ** FAIL **"))
    if got < 4:
        print(f"  ** only {got}/4 bases produced a surviving lift in {attempts} tries **")

def cmd_stats(st, name, seconds):
    # per-seed output so parallel/chunked stats runs never clobber each other;
    # checkpoint after every completed run so a killed process loses nothing.
    out = os.path.join(_HERE, f'best_{name}_s{st["seed"]}.pkl')
    t0 = time.time(); res = []
    while time.time()-t0 < seconds:
        base, g = make_base(st)
        v, p, _ = gauss_newton(st, base, g)
        if v is not None:
            res.append((v, base, p))
            pickle.dump(sorted(res, key=lambda x: x[0])[:10], open(out,'wb'))
            if len(res) % 5 == 1:
                print(f"run {len(res):3d}: {v:.3e}  t={time.time()-t0:.0f}s", flush=True)
    res.sort(key=lambda x: x[0])
    print(f"\n{name} (seed {st['seed']}): {len(res)} runs | min {res[0][0]:.3e} | "
          f"median {res[len(res)//2][0]:.3e}")
    if res[0][0] < 1e-8:
        print("*** ALERT: converging run — possible SOLUTION. Verify exactly. ***")

def cmd_merge(st, name):
    # pickle is safe here: these files are produced by this same harness in
    # this workspace (session artifacts), never untrusted input.
    import glob
    res = []
    for fn in sorted(glob.glob(os.path.join(_HERE, f'best_{name}_s*.pkl'))):
        chunk = pickle.load(open(fn,'rb'))
        res += chunk
        print(f"  {os.path.basename(fn)}: {len(chunk)} kept, min {min(v for v,_,_ in chunk):.3e}")
    res.sort(key=lambda x: x[0])
    vals = [v for v,_,_ in res]
    print(f"{name}: merged {len(res)} | min {vals[0]:.3e} | median {vals[len(vals)//2]:.3e}")
    pickle.dump(res[:10], open(os.path.join(_HERE, f'best_{name}.pkl'),'wb'))

def cmd_polish(st, name):
    best = pickle.load(open(f'best_{name}.pkl','rb'))
    for k,(v0, base, p) in enumerate(best[:4]):
        g = np.array([gr(*base, st['Psi0']) for gr in st['grads']], complex)
        v, p2, traj = gauss_newton(st, base, g, p=p, iters=60, h=1e-5, central=True)
        print(f"cand {k}: {v0:.3e} -> {v:.3e}  tail {[f'{t:.1e}' for t in traj[-5:]]}")
        pickle.dump((base, p2, v), open(f'polished_{name}_{k}.pkl','wb'))

if __name__ == '__main__':
    np.seterr(all='ignore')
    mode, name = sys.argv[1], sys.argv[2]
    st = setup(name)
    if mode == 'control': cmd_control(st, name)
    elif mode == 'stats': cmd_stats(st, name, int(sys.argv[3]) if len(sys.argv)>3 else 600)
    elif mode == 'merge': cmd_merge(st, name)
    elif mode == 'polish': cmd_polish(st, name)
