"""
incremental_lift.py — EXACT incremental jet-lift for the T5 certificate.

Unlike exact_jetlift.py (which forms the whole substituted polynomial and chokes
on the degree-25 anchor), this computes slice by slice from truncated series,
never forming the full product — the discipline that makes jetlift.py fast,
now done exactly over F_p with a fixed generic base and symbolic free params.

Pipeline (f31, sub2, stripped vars d~_k of degree 2w_k):
  * fix base (a0,b0,c0,e0) on the bottom slice over F_p;
  * lift slices j=1..10: each is LINEAR in the new coeffs; pivot on one
    (gradient component), the rest become free params (18 total);
  * consistency slices j=11..M: polynomials in the 18 free params; add them to
    an ideal one at a time and Groebner-test for 1 (inconsistency) => the fixed
    base admits no completion => obstruction height M found.

Truncated arithmetic mod y^{M+1}; coefficients are Poly over GF(p) in the free
params.  Usage: python3 incremental_lift.py [p] [Mmax]
"""
import sys, time, random
import sympy as sp
from sympy import symbols, GF, Poly

P = int(sys.argv[1]) if len(sys.argv) > 1 else 32003
MMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 30
DOM = GF(P)
y = symbols('y')

# --- f31 as a list of (coeff, (a,b,c,e,f)) exponent tuples on (d2,d1,d0,m1,P)
def load_f31():
    d2,d1,d0,m1,Pi = symbols('d2 d1 d0 m1 P')
    f = sp.Poly(sp.sympify(open('f31_deg31.txt').read().replace('^','**')),
                d2,d1,d0,m1,Pi)
    return [(int(c) % P, m) for c, m in zip(f.coeffs(), f.monoms())]

F31 = load_f31()

def phit_series(M):
    # Phi~ = -(y+1)^30 (2048y^4-...+195)/6630 truncated to degree M (coeffs in F_p)
    q = sp.Poly(sp.expand(-(y+1)**30*(2048*y**4-512*y**3+320*y**2-240*y+195))
                * sp.invert(6630, P), y, modulus=P)
    co = [0]*(M+1)
    for c, (e,) in zip(q.coeffs(), q.monoms()):
        if e <= M: co[e] = int(c) % P
    return co

def cmul(u, v, M):
    """truncated product of coeff-lists (entries: python ints OR Poly), mod y^{M+1}."""
    r = [0]*(M+1)
    for i, ui in enumerate(u):
        if ui == 0 or i > M: continue
        for j, vj in enumerate(v):
            if vj == 0 or i+j > M: continue
            r[i+j] = (r[i+j] + ui*vj)
    return [ (t % P if isinstance(t, int) else t) for t in r ]

def cpow(u, n, M):
    r = [1] + [0]*M
    base = u[:]
    while n:
        if n & 1: r = cmul(r, base, M)
        n >>= 1
        if n: base = cmul(base, base, M)
    return r

ALLSYMS = []   # filled in main(); used to reduce coeffs mod P

def redp(e):
    if isinstance(e, int): return e % P
    if e == 0 or e.is_number: return int(e) % P
    return sp.Poly(e, *ALLSYMS, domain=DOM).as_expr()

def slice_j(dt, phit, j):
    """coeff of y^j in f31(d~,Phi~), computed mod y^{j+1} (truncate at j)."""
    m = j
    tot = 0
    for c, (a,b,cc,e,f) in F31:
        term = [c] + [0]*m
        for u, p in zip(dt, (a,b,cc,e)):
            if p: term = cmul(term, cpow(u, p, m), m)
        if f: term = cmul(term, cpow(phit, f, m), m)
        tot = tot + term[j]
    return redp(tot)

def main():
    rng = random.Random(2026)
    t0 = time.time()
    # free-param symbols for the non-pivot new coeffs; we allocate all then pick
    # base over F_p
    d2,d1,d0,m1,Pi = symbols('d2 d1 d0 m1 P')
    f31e = sp.sympify(open('f31_deg31.txt').read().replace('^','**'))
    phi0 = (-sp.invert(34, P)) % P
    base = None
    for _ in range(60):
        a0,b0,c0 = (rng.randrange(1,P) for _ in range(3))
        poly = sp.Poly(f31e.subs({d2:a0,d1:b0,d0:c0,Pi:phi0}), m1, modulus=P)
        rts = sp.ground_roots(poly)
        for r in rts:
            if int(r) % P: base = (a0%P,b0%P,c0%P,int(r)%P); break
        if base: break
    if not base:
        print("no base found"); return
    a0,b0,c0,e0 = base
    print(f"[incremental lift] p={P} base={base}", flush=True)

    M = MMAX
    # unknown coeffs are fresh symbols, base index 0 fixed
    asy = [symbols(f'a{i}') for i in range(1,5)]
    bsy = [symbols(f'b{i}') for i in range(1,7)]
    csy = [symbols(f'c{i}') for i in range(1,9)]
    esy = [symbols(f'e{i}') for i in range(1,11)]
    global ALLSYMS; ALLSYMS = asy+bsy+csy+esy
    phit = phit_series(M)
    d2c = [a0] + asy + [0]*(M-4)
    d1c = [b0] + bsy + [0]*(M-6)
    d0c = [c0] + csy + [0]*(M-8)
    emc = [e0] + esy + [0]*(M-10)
    dt = [d2c, d1c, d0c, emc]
    newarr = {'a':(d2c,4),'b':(d1c,6),'c':(d0c,8),'e':(emc,10)}

    free = []
    def subst_all(s, sol):
        for arr in dt:
            for t in range(len(arr)):
                v = arr[t]
                if not isinstance(v,int) and v.has(s):
                    arr[t] = redp(v.subs(s, sol))
    # LIFT slices 1..10: each linear in the new coeffs; pivot on one
    for j in range(1, 11):
        sj = slice_j(dt, phit, j)                    # expr, linear in new coeffs
        avail = [(nm+str(j), nm) for nm,(arr,mx) in newarr.items() if j <= mx]
        piv = None
        for name, nm in avail:
            s = symbols(name)
            co = sj.coeff(s, 1) if not isinstance(sj,int) else 0
            if co != 0 and co.is_number:
                piv = (name, nm, s, int(co) % P); break
        if piv is None:
            print(f"  slice {j}: no numeric pivot — stop"); break
        name, nm, s, alpha = piv
        beta = sj - sj.coeff(s,1)*s                  # sj = alpha*s + beta
        sol = redp(-beta * pow(alpha, -1, P))
        newarr[nm][0][j] = sol
        subst_all(s, sol)
        for name2, nm2 in avail:
            s2 = symbols(name2)
            if name2 != name and s2 not in free: free.append(s2)
        print(f"  slice {j}: pivot {name}; free={len(free)}  t={time.time()-t0:.0f}s", flush=True)

    print(f"residual free dof = {len(free)}", flush=True)
    # CONSISTENCY slices 11..M: incremental Groebner inconsistency test
    G = []
    for j in range(11, M+1):
        sj = slice_j(dt, phit, j)
        if isinstance(sj, int) or sj.is_number:
            if int(sj) % P != 0:
                print(f"  consistency slice {j}: nonzero CONSTANT -> INCONSISTENT"); return
            print(f"  consistency slice {j}: trivial 0   t={time.time()-t0:.0f}s", flush=True); continue
        G.append(sj)
        gb = sp.groebner(G, *free, order='grevlex', modulus=P)
        one = (1 in [g for g in gb.exprs]) or (len(gb.exprs)==1 and gb.exprs[0]==1)
        print(f"  consistency slice {j}: gens {len(G)}, GB {len(gb.exprs)}"
              f"{'  *** INCONSISTENT ***' if one else ''}  t={time.time()-t0:.0f}s", flush=True)
        if one:
            print(f"\nOBSTRUCTION: no completion over this base; first inconsistency "
                  f"at consistency slice y^{j}.")
            return
    print(f"\nno inconsistency through slice {M} (need higher M or special base).")

if __name__ == '__main__':
    main()
