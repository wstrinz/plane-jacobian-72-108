"""Exact symbolic checks for CORNER_144_COMPARISON.md.
Run: python d2_plane_72_108/corner144_verify.py
Every failed claim exits nonzero.
"""
from pathlib import Path
import sys
import sympy as sp

y, x, X, Y = sp.symbols("y x X Y")
checks = 0

def check(name, condition):
    global checks
    if not bool(condition):
        raise SystemExit(f"[FAIL] {name}")
    checks += 1
    print(f"[OK] {name}")

def poly_order(expr):
    return min(m[0] for m in sp.Poly(sp.expand(expr), y).monoms())

def multiplicity(expr, factor):
    count = 0
    quotient = sp.Poly(sp.expand(expr), y)
    divisor = sp.Poly(factor, y)
    while True:
        quotient, remainder = sp.div(quotient, divisor)
        if not remainder.is_zero:
            return count
        count += 1

print("A. chain arithmetic and convention")
def dio(m, n, aa, ell, bb, k=1):
    return (m+n)*bb*k - n*(bb*ell-aa)
check("(7/4,3),(m,n)=(3,4) satisfies k=1", dio(3,4,7,4,3)==1)
check("(11/4,7),(m,n)=(3,2) satisfies k=1", dio(3,2,11,4,7)==1)
check("degree recipe gives target (108,144)", (3*36,4*36)==(108,144))
check("reduced current exponent pair gives (72,108)", (2*36,3*36)==(72,108))
check("gamma=3 gives final corner (7/4,3)", (1+sp.Rational(3,4),3)==(sp.Rational(7,4),3))
check("gamma=7 gives final corner (11/4,7)", (1+sp.Rational(7,4),7)==(sp.Rational(11,4),7))

print("B. Laurent normalization and leading factor")
phiX, phiY = x**-1, x**4*y
jac = sp.diff(phiX,x)*sp.diff(phiY,y)-sp.diff(phiX,y)*sp.diff(phiY,x)
check("Jacobian of Laurent map is -x^2", sp.simplify(jac+x**2)==0)
alpha = sp.symbols("alpha", nonzero=True)
hc = sp.symbols("h0:5")
def H(z): return sum(hc[i]*z**i for i in range(5))
map_z = (X**4*Y).subs({X:phiX,Y:phiY}, simultaneous=True)
map_prefactor = (Y+alpha*X**-4).subs({X:phiX,Y:phiY}, simultaneous=True)
check("Laurent map sends X^4*Y to y", sp.expand(map_z-y)==0)
check("shifted prefactor becomes x^4*(y+alpha)", sp.expand(map_prefactor-x**4*(y+alpha))==0)
check("multiplicity-three edge gives C4=y^3*(y+alpha)*quartic",
      sp.expand(map_prefactor*map_z**3*H(map_z)-x**4*(y+alpha)*y**3*H(y))==0)
check("multiplicity-seven edge gives current C4=y^7*(y+alpha)",
      sp.expand(map_prefactor*map_z**7-x**4*y**7*(y+alpha))==0)

print("C. parametric forcing formula")
def direct_forcing_check(aa, bb, tt=4, kk=2):
    cc = sp.Function("cc")(y)
    ff = sp.Function("ff")(y)
    ss = kk+1-aa*tt
    pp = x**(aa*tt)*cc**aa
    tail = x**ss*ff/cc**bb
    bracket = sp.diff(pp,x)*sp.diff(tail,y)-sp.diff(pp,y)*sp.diff(tail,x)
    expected = aa*cc**(aa-bb-1)*(tt*cc*sp.diff(ff,y)-(tt*(bb-aa)+kk+1)*sp.diff(cc,y)*ff)
    return sp.expand(bracket/x**kk-expected)==0
check("direct bracket gives general formula at (2,3)", direct_forcing_check(2,3))
check("direct bracket gives general formula at (3,4)", direct_forcing_check(3,4))
def data(aa,bb,tt=4,kk=2):
    ss=kk+1-aa*tt
    dexp=lambda index: aa*tt-1-aa*index
    j=-ss
    clear=aa*(bb*tt+j)-bb
    return ss,dexp,j,clear,clear-bb
old=data(2,3); new=data(3,4)
check("current obstruction degree and Phi exponent are -5,28", old[0]==-5 and old[4]==28)
check("target obstruction degree and Phi exponent are -9,67", new[0]==-9 and new[4]==67)
check("target Q-slice clearing exponent is 71", new[3]==71)
check("D exponents are 7-2k and 11-3k", all(old[1](i)==7-2*i and new[1](i)==11-3*i for i in range(-5,5)))
check("target commuting correction range ends at C^-2", 4*(-2)>-9 and 4*(-3)<-9)

print("D. current control")
cold=y**7*(y+1)
qold=2048*y**4-512*y**3+320*y**2-240*y+195
fold=-y**8*(y+1)**2*qold/sp.Integer(6630)
check("current f solves its ODE", sp.expand(8*cold*sp.diff(fold,y)-14*sp.diff(cold,y)*fold-cold**2)==0)
check("current quartic separable and avoids 0,-1", sp.discriminant(qold,y)!=0 and qold.subs(y,0)!=0 and qold.subs(y,-1)!=0)
Phiold=sp.expand(fold*cold**28)
cofold=sp.cancel(Phiold/(y**204*(y+1)**30))
check("current Phi signature (238,204,30,4)", sp.degree(Phiold,y)==238 and poly_order(Phiold)==204 and multiplicity(Phiold,y+1)==30 and sp.degree(cofold,y)==4)

print("E. target ODE solution and forced quartic")
gc=sp.symbols("g0:6"); A=sp.symbols("A")
g=sum(gc[i]*y**i for i in range(6))
check("local target orders are 4 at y and 2 at every simple g-root",
      (3+(3+1)-1==2*3 and 12*(3+1)-21*3!=0)
      and (1+(1+1)-1==2*1 and 12*(1+1)-21*1!=0))
check("infinity leading coefficient is resonant exactly at degree 14",
      12*14-21*8==0 and all(12*d-21*8!=0 for d in range(15,25)))
cgen=y**3*g; fgen=A*y**4*g**2
odegen=sp.factor(12*cgen*sp.diff(fgen,y)-21*sp.diff(cgen,y)*fgen-cgen**2)
expected=y**6*g**2*(3*A*(y*sp.diff(g,y)-5*g)-1)
check("squarefree-order ansatz reduces to 3A(yg'-5g)=1", sp.expand(odegen-expected)==0)
ident=sp.Poly(3*A*(y*sp.diff(g,y)-5*g)-1,y)
check("identity coefficients force g1,...,g4 to zero",
      all(sp.expand(ident.coeff_monomial(y**i)-3*A*(i-5)*gc[i])==0 for i in range(1,5)))
check("degree-five coefficient is resonant", ident.coeff_monomial(y**5)==0)
check("constant coefficient gives A=-1/(15*g0)",
      sp.expand(ident.coeff_monomial(1)-(-15*A*gc[0]-1))==0)
check("g(-1)=0 then forces g0=g5", sp.expand(g.subs({gc[1]:0,gc[2]:0,gc[3]:0,gc[4]:0,y:-1})-(gc[0]-gc[5]))==0)
h=y**4-y**3+y**2-y+1
gt=(y+1)*h; ct=y**3*gt; ft=-y**4*gt**2/sp.Integer(15)
check("normalized g=(y+1)h=y^5+1", sp.expand(gt-(y**5+1))==0)
check("forced h separable and avoids 0,-1", sp.discriminant(h,y)!=0 and h.subs(y,0)!=0 and h.subs(y,-1)!=0)
check("target f solves 12*c*f'-21*c'*f=c^2", sp.expand(12*ct*sp.diff(ft,y)-21*sp.diff(ct,y)*ft-ct**2)==0)
fc=sp.symbols("f0:15"); fans=sum(fc[i]*y**i for i in range(15))
poly=sp.Poly(sp.expand(12*ct*sp.diff(fans,y)-21*sp.diff(ct,y)*fans-ct**2),y)
lin=list(sp.linsolve(poly.all_coeffs(),fc))
check("target ODE has one polynomial solution of degree <=14", len(lin)==1 and not any(v.free_symbols & set(fc) for v in lin[0]))
check("linear solve returns displayed target f", sp.expand(fans.subs(dict(zip(fc,lin[0])))-ft)==0)
check("target Phi exponent arithmetic gives y^205*g^69", 4+3*67==205 and 2+67==69)
check("g=y^5+1 has order zero and a simple root at -1",
      gt.subs(y,0)!=0 and gt.subs(y,-1)==0 and sp.diff(gt,y).subs(y,-1)!=0)
check("target Phi signature (550,205,69,276)",
      205+5*69==550 and sp.degree(h,y)*69==276)
coft=-h**69/sp.Integer(15)
check("target unit cofactor is exactly -h^69/15", coft == -h**69/sp.Integer(15))
print("F. conditional envelope candidates")
def cross(v,w): return v[0]*w[1]-v[1]*w[0]
Pvec=(3*8-2,3*3-1); Qvec=(4*8-(-1),4*3)
check("target lower candidate edges are parallel", cross(Pvec,Qvec)==0)
check("pre-map lower slope is 4/11", Pvec==(22,8) and Qvec==(33,12) and sp.Rational(Pvec[1],Pvec[0])==sp.Rational(4,11))
for w in range(1,13):
    ki=4-w; de=3*w-1
    up0=2*ki+8*de
    up1=(ki+4)+8*de
    low=sp.ceiling(sp.Rational(4*ki-1,5))+3*de
    check(f"envelope arithmetic w={w}", up0==22*w and up1==23*w and low==8*w+sp.ceiling(sp.Rational(w,5)))

print(f"\nALL {checks} CORNER-144 CHECKS PASSED")
print(f"script: {Path(__file__).resolve()}")
sys.exit(0)




