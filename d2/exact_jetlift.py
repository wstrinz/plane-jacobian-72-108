"""
exact_jetlift.py — exact (finite-field) jet-lift pilot for the T5 certificate.

Mirrors jetlift.py's numeric lift, but EXACT over F_p, to turn "positive floor"
into "provably no solution". For a fixed generic base (trailing coeffs
a0,b0,c0,e0 on the bottom-slice variety), the remaining window coefficients are
kept SYMBOLIC; we impose the y-slice equations of f31(d~(y),Phi~(y)) == 0 and
ask Groebner whether they are inconsistent (1 in ideal) — i.e. no lift completes
over that base — and at which slice height the inconsistency first appears.

Stripped variables (ord bound y^{12w} divided out): d~_k has degree 2w_k (sub2):
  d~2 deg4 (a0..a4), d~1 deg6 (b0..b6), d~0 deg8 (c0..c8), d~-1 deg10 (e0..e10).
  Phi~ = -(y+1)^30 (2048y^4-512y^3+320y^2-240y+195)/6630  (deg 34).
Identity f31(d~,Phi~) == 0 has y-degree <= 250; slice j = coeff of y^j.

Usage:  python3 exact_jetlift.py [p] [Mmax]
Reports: base found, residual dof, and the minimal slice height M at which the
fixed-base free-coefficient system becomes inconsistent (the obstruction slice),
or "no inconsistency up to Mmax".
"""
import sys, subprocess, random
import sympy as sp

P = int(sys.argv[1]) if len(sys.argv) > 1 else 32003
MMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 60
FACTOR = 'f31_deg31.txt'
SIZES = dict(a=5, b=7, c=9, e=11)          # coeff counts (degrees 4,6,8,10)

def phit_coeffs(p):
    y = sp.symbols('y')
    Phit = sp.Poly(sp.expand(-(y+1)**30*(2048*y**4-512*y**3+320*y**2-240*y+195))
                   * sp.invert(6630, p), y, modulus=p)
    return [int(c) % p for c in Phit.all_coeffs()[::-1]]   # ascending

def find_base(p, f31_str, rng):
    """Pick random a0,b0,c0 in F_p; solve the bottom slice f31(a0,b0,c0,e0,-1/34)=0
    for a root e0 in F_p. Returns (a0,b0,c0,e0) or None."""
    d2,d1,d0,m1,Pi = sp.symbols('d2 d1 d0 m1 P')
    f31 = sp.sympify(f31_str.replace('^','**'))
    phi0 = (-sp.invert(34, p)) % p
    for _ in range(40):
        a0,b0,c0 = (rng.randrange(1,p) for _ in range(3))
        poly = sp.Poly(f31.subs({d2:a0,d1:b0,d0:c0,Pi:phi0}), m1, modulus=p)
        roots = sp.ground_roots(poly)      # roots in F_p with multiplicity
        for r,_ in roots.items():
            e0 = int(r) % p
            if e0 != 0:
                return (a0%p, b0%p, c0%p, e0)
    return None

def build_singular(p, base, phit, f31_str, Mmax):
    a0,b0,c0,e0 = base
    # symbolic coeffs: a1..a4,b1..b6,c1..c8,e1..e10 (base index 0 fixed numeric)
    avars = ([f'a{i}' for i in range(1,5)] + [f'b{i}' for i in range(1,7)]
             + [f'c{i}' for i in range(1,9)] + [f'e{i}' for i in range(1,11)])
    ndof = len(avars)                       # 4+6+8+10 = 28 symbolic coeffs
    d2s = f'{a0}+' + '+'.join(f'a{i}*y^{i}' for i in range(1,5))
    d1s = f'{b0}+' + '+'.join(f'b{i}*y^{i}' for i in range(1,7))
    d0s = f'{c0}+' + '+'.join(f'c{i}*y^{i}' for i in range(1,9))
    ems = f'{e0}+' + '+'.join(f'e{i}*y^{i}' for i in range(1,11))
    phit = '+'.join(f'{c}*y^{i}' for i,c in enumerate(phit) if c)
    sc = f"""ring R={p},(d2,d1,d0,m1,P,{','.join(avars)},y),dp;
poly f31={f31_str.strip()};
poly F=subst(f31,d2,{d2s},d1,{d1s},d0,{d0s},m1,{ems},P,{phit});
matrix C=coeffs(F,y);
// slice j = coeff of y^j is C[j+1,1]. slice 0 is 0 by base construction.
ideal I;
int j; for(j=2;j<={Mmax+1};j++){{ I[j-1]=C[j,1]; }}
"residual symbolic dof:",{ndof};
"num slice equations imposed:",{Mmax};
option(redSB);
ideal G=groebner(I);
"GB size:",size(G);
if(G[1]==1){{ "INCONSISTENT: 1 in ideal => no window solution over this base"; }}
else{{ "consistent so far; GB dim:",dim(G),"  deg:",vdim(G); }}
quit;
"""
    return sc, ndof

def main():
    rng = random.Random(2026)
    f31_str = open(FACTOR).read()
    base = find_base(P, f31_str, rng)
    if base is None:
        print("no base found over F_p; try another prime"); return
    print(f"[exact jet-lift] p={P}, base (a0,b0,c0,e0)={base}")
    sc, ndof = build_singular(P, base, phit_coeffs(P), f31_str, MMAX)
    open('/tmp/exact_jl.sing','w').write(sc)
    print(f"residual symbolic dof = {ndof}; running Singular (Mmax={MMAX})...", flush=True)
    out = subprocess.run(['Singular','-q','/tmp/exact_jl.sing'],
                         capture_output=True, text=True, timeout=100000)
    print(out.stdout.strip())
    if out.stderr.strip():
        print("STDERR:", out.stderr.strip()[:500])

if __name__ == '__main__':
    main()
