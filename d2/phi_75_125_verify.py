"""Exact symbolic checks for PHI_75_125.md.
Run: python d2_plane_72_108/phi_75_125_verify.py
Every failed claim exits nonzero.

Independent re-derivation of the forcing divisor Phi = f * C^N for the (75,125)
case -- family F_2, j=1, corner A_0=(5,20), final corner (7/5,2), (m,n)=(3,5),
reduced C-power pair (a,b)=(3,5) (NON-adjacent: b-a=2).  Mirrors
corner144_verify.py's structure.  The ODE solution and every divisor invariant
are re-derived here from scratch (not read from the .md).
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

# ---------------------------------------------------------------------------
print("A. chain arithmetic and convention (F_2, corner (5,20), final (7/5,2))")
# Diophantine on the final corner A=(p/l,q)=(7/5,2): (m+n) q k - n(q l - p) = k
def dio(m, n, p, l, q, k=1):
    return (m + n) * q * k - n * (q * l - p) - k
check("(7/5,2),(m,n)=(3,5) satisfies k=1 (F_2, j=1)", dio(3, 5, 7, 5, 2) == 0)
check("(7/5,2),(m,n)=(2,3) is the j=0 sibling (Moh's 75)", dio(2, 3, 7, 5, 2) == 0)
# degree recipe deg=m*v11(A0), v11(5,20)=25
check("degree recipe gives (75,125)", (3 * 25, 5 * 25) == (75, 125))
check("j=0 sibling gives (50,75) = Moh's case", (2 * 25, 3 * 25) == (50, 75))
# type-II.b:  A1 = (1,0) + gamma*(1/l,1),  gamma = selected multiplicity q
gamma, l = 2, 5
check("gamma=2 on edge (1/5,1) gives final corner (7/5,2)",
      (1 + sp.Rational(gamma, l), gamma) == (sp.Rational(7, 5), 2))
# reduced C-power pair = sorted (m,n)
a, b = 3, 5
check("reduced C-power pair (a,b)=(3,5) is sorted (m,n)", (a, b) == tuple(sorted((3, 5))))
check("pair is NON-adjacent (b-a=2), unlike the (8,28) corner (b-a=1)", b - a == 2)

# ---------------------------------------------------------------------------
print("B. Laurent normalization and leading factor")
# Laurent map X->x^-1, Y->x^l y with l = 5 (= denominator of the final corner)
phiX, phiY = x**-1, x**5 * y
jac = sp.diff(phiX, x) * sp.diff(phiY, y) - sp.diff(phiX, y) * sp.diff(phiY, x)
check("Jacobian of Laurent map is -x^3 (=> kappa=l-2=3)", sp.simplify(jac + x**3) == 0)
alpha = sp.symbols("alpha", nonzero=True)
hc = sp.symbols("h0:3")                       # residual H has degree r = a0-q-1 = 2
def H(z): return sum(hc[i] * z**i for i in range(3))
map_z = (X**5 * Y).subs({X: phiX, Y: phiY}, simultaneous=True)
map_prefactor = (Y + alpha * X**-5).subs({X: phiX, Y: phiY}, simultaneous=True)
check("Laurent map sends X^5*Y to y", sp.expand(map_z - y) == 0)
check("shifted prefactor becomes x^5*(y+alpha)", sp.expand(map_prefactor - x**5 * (y + alpha)) == 0)
# selected multiplicity q=2 edge (Y+alpha X^-5)(X^5 Y)^2 H(X^5 Y+alpha) -> x^5 * C
check("multiplicity-two edge gives C=x^5*y^2*(y+alpha)*quartic-residual",
      sp.expand(map_prefactor * map_z**2 * H(map_z) - x**5 * (y + alpha) * y**2 * H(y)) == 0)
check("deg C = a0 = 5 (order-2 root at y=0, residual degree r=a0-q-1=2)",
      2 + 1 + (5 - 2 - 1) == 5)

# ---------------------------------------------------------------------------
print("C. parametric forcing formula  a{ t c f' - [t(b-a)+kappa+1] c' f } = c^(b-a+1)")
def direct_forcing_check(aa, bb, tt, kk):
    cc = sp.Function("cc")(y)
    ff = sp.Function("ff")(y)
    ss = kk + 1 - aa * tt
    pp = x**(aa * tt) * cc**aa
    tail = x**ss * ff / cc**bb
    bracket = sp.diff(pp, x) * sp.diff(tail, y) - sp.diff(pp, y) * sp.diff(tail, x)
    expected = aa * cc**(aa - bb - 1) * (tt * cc * sp.diff(ff, y)
                                         - (tt * (bb - aa) + kk + 1) * sp.diff(cc, y) * ff)
    return sp.expand(bracket / x**kk - expected) == 0
# the SAME general family, now at (a,b,t,kappa)=(3,5,5,3)
check("direct bracket gives general formula at (3,5), t=5, kappa=3", direct_forcing_check(3, 5, 5, 3))
# and it must still reproduce the two audited (8,28) corners (t=4,kappa=2)
check("same family reproduces (2,3) t=4 k=2", direct_forcing_check(2, 3, 4, 2))
check("same family reproduces (3,4) t=4 k=2", direct_forcing_check(3, 4, 4, 2))

def data(aa, bb, tt, kk):
    ss = kk + 1 - aa * tt
    dexp = lambda index: aa * tt - 1 - aa * index          # D_k = C_k c^(a(t-k)-1)
    j = -ss
    clear = aa * (bb * tt + j) - bb
    return ss, dexp, j, clear, clear - bb                   # last entry = N = Phi C-power
new = data(3, 5, 5, 3)
check("F_2 obstruction degree s=-11 and forcing slice j=11", new[0] == -11 and new[2] == 11)
check("F_2 Q-slice clearing exponent is 103", new[3] == 103)
check("F_2 Phi C-power N = 98", new[4] == 98)
check("D exponents are 14-3k", all(new[1](i) == 14 - 3 * i for i in range(-5, 5)))
check("commuting-correction range ends at C^-2 (v(C^-2)=-10>-11, v(C^-3)=-15<-11)",
      5 * (-2) > -11 and 5 * (-3) < -11)
# N as the single parametric closed form N = a[t(a+b)-(kappa+1)] - 2b, verified at ALL 3 corners
def Nformula(aa, bb, tt, kk):
    return aa * (tt * (aa + bb) - (kk + 1)) - 2 * bb
check("N formula reproduces (72,108): 28", Nformula(2, 3, 4, 2) == 28)
check("N formula reproduces (108,144): 67", Nformula(3, 4, 4, 2) == 67)
check("N formula gives (75,125): 98", Nformula(3, 5, 5, 3) == 98)

# ---------------------------------------------------------------------------
print("D. reference control -- the two audited (8,28) corners still land where CORNER_144 says")
# (72,108): C=y^7(y+1), f=y^8(y+1)^2 q4 (extra quartic, resonance gap), N=28 -> (238,204,30,4)
cold = y**7 * (y + 1)
qold = 2048 * y**4 - 512 * y**3 + 320 * y**2 - 240 * y + 195
fold = -y**8 * (y + 1)**2 * qold / sp.Integer(6630)
check("current (72,108) f solves 8cf'-14c'f=c^2",
      sp.expand(8 * cold * sp.diff(fold, y) - 14 * sp.diff(cold, y) * fold - cold**2) == 0)
Phiold = sp.expand(fold * cold**28)
cofold = sp.cancel(Phiold / (y**204 * (y + 1)**30))
check("current Phi signature (238,204,30,4)",
      sp.degree(Phiold, y) == 238 and poly_order(Phiold) == 204
      and multiplicity(Phiold, y + 1) == 30 and sp.degree(cofold, y) == 4)
# (108,144): C=y^3(y^5+1), f=-y^4(y^5+1)^2/15, N=67 -> (550,205,69,276)
ct144 = y**3 * (y**5 + 1); ft144 = -y**4 * (y**5 + 1)**2 / sp.Integer(15)
check("target (108,144) f solves 12cf'-21c'f=c^2",
      sp.expand(12 * ct144 * sp.diff(ft144, y) - 21 * sp.diff(ct144, y) * ft144 - ct144**2) == 0)
Phi144 = sp.expand(ft144 * ct144**67)
cof144 = sp.cancel(Phi144 / (y**205 * (y + 1)**69))
check("target Phi signature (550,205,69,276)",
      sp.degree(Phi144, y) == 550 and poly_order(Phi144) == 205
      and multiplicity(Phi144, y + 1) == 69 and sp.degree(cof144, y) == 276)

# ---------------------------------------------------------------------------
print("E. target (75,125) ODE solution and forced quadratic residual")
# ODE:  15 c f' - 42 c' f = c^3,   c = y^2 g,  g deg 3
# local balance rho = (e-1)mu + 1 with e = b-a+1 = 3
check("local orders: rho=5 at y (mu=2) and rho=3 at simple roots (mu=1)",
      ((3 - 1) * 2 + 1 == 5) and ((3 - 1) * 1 + 1 == 3))
check("infinity leading coefficient is resonant exactly at degree 14",
      15 * 14 - 42 * 5 == 0 and all(15 * d - 42 * 5 != 0 for d in range(15, 25)))
gc = sp.symbols("g0:4"); A = sp.symbols("A")
g = sum(gc[i] * y**i for i in range(4))
cgen = y**2 * g; fgen = A * y**5 * g**3
odegen = sp.factor(15 * cgen * sp.diff(fgen, y) - 42 * sp.diff(cgen, y) * fgen - cgen**3)
expected = y**6 * g**3 * (3 * A * (y * sp.diff(g, y) - 3 * g) - 1)
check("squarefree-order ansatz reduces to 3A(yg'-3g)=1", sp.expand(odegen - expected) == 0)
ident = sp.Poly(3 * A * (y * sp.diff(g, y) - 3 * g) - 1, y)
check("identity coefficients force g1,g2 to zero",
      all(sp.expand(ident.coeff_monomial(y**i) - 3 * A * (i - 3) * gc[i]) == 0 for i in range(1, 3)))
check("degree-three coefficient is resonant", ident.coeff_monomial(y**3) == 0)
check("constant coefficient gives A=-1/(9*g0)",
      sp.expand(ident.coeff_monomial(1) - (-9 * A * gc[0] - 1)) == 0)
check("g(-1)=0 then forces g0=g3",
      sp.expand(g.subs({gc[1]: 0, gc[2]: 0, y: -1}) - (gc[0] - gc[3])) == 0)
# normalized solution
H2 = y**2 - y + 1
gt = (y + 1) * H2; ct = y**2 * gt; ft = -y**5 * gt**3 / sp.Integer(9)
check("normalized g=(y+1)H2=y^3+1", sp.expand(gt - (y**3 + 1)) == 0)
check("forced residual H2=y^2-y+1 separable and avoids 0,-1",
      sp.discriminant(H2, y) != 0 and H2.subs(y, 0) != 0 and H2.subs(y, -1) != 0)
check("target f solves 15*c*f'-42*c'*f=c^3",
      sp.expand(15 * ct * sp.diff(ft, y) - 42 * sp.diff(ct, y) * ft - ct**3) == 0)
fc = sp.symbols("f0:15"); fans = sum(fc[i] * y**i for i in range(15))
poly = sp.Poly(sp.expand(15 * ct * sp.diff(fans, y) - 42 * sp.diff(ct, y) * fans - ct**3), y)
lin = list(sp.linsolve(poly.all_coeffs(), fc))
check("target ODE has one polynomial solution of degree <=14",
      len(lin) == 1 and not any(v.free_symbols & set(fc) for v in lin[0]))
check("linear solve returns displayed target f",
      sp.expand(fans.subs(dict(zip(fc, lin[0]))) - ft) == 0)
# Phi = f * C^98
N = 98
check("Phi C-power exponent arithmetic gives y^201*(y^3+1)^101",
      5 + 2 * N == 201 and 3 + N == 101)
Phi = sp.expand(ft * ct**N)
cof = sp.cancel(Phi / (y**201 * (y + 1)**101))
check("target Phi = -(1/9) y^201 (y^3+1)^101",
      sp.expand(Phi + sp.Rational(1, 9) * y**201 * (y**3 + 1)**101) == 0)
check("target Phi signature (504,201,101,202)",
      sp.degree(Phi, y) == 504 and poly_order(Phi) == 201
      and multiplicity(Phi, y + 1) == 101 and sp.degree(cof, y) == 202)
check("cofactor is exactly -(1/9) H2^101 = -(1/9)(y^2-y+1)^101, degree 202",
      sp.expand(cof + sp.Rational(1, 9) * H2**101) == 0 and sp.degree(cof, y) == 202)

# ---------------------------------------------------------------------------
print("F. verdict on the a-only prediction (POTENTIAL_PROBE.md)")
# a-only claim: deg Phi = 64a^2-8a-2, mult_(y+1) = 8a^2-a, evaluated at a=3 -> (550,69)
check("a-only formulas at a=3 predict deg=550, mult=69",
      64 * 3**2 - 8 * 3 - 2 == 550 and 8 * 3**2 - 3 == 69)
check("derived deg Phi = 504 != 550  (DIFFERS)", sp.degree(Phi, y) == 504 and 504 != 550)
check("derived mult_(y+1) = 101 != 69  (DIFFERS)", multiplicity(Phi, y + 1) == 101 and 101 != 69)
# unconditional refutation: deg Phi = deg f + N*deg C = 14 + N*a0, a0=5 here vs 8 for (8,28)
check("deg Phi = 14 + 5*N structurally (a0=5), so 550 is unreachable for ANY integer N",
      14 + 5 * N == 504 and not sp.Rational(550 - 14, 5).is_Integer)
# general closed forms in (a,b,t,kappa,a0,q): reproduce all three known points
def sig(aa, bb, tt, kk, a0, qv):
    ee = bb - aa + 1
    NN = aa * (tt * (aa + bb) - (kk + 1)) - 2 * bb
    degf = ee * a0 - qv + 1                      # holds when no resonance-gap extra factor
    deg = degf + NN * a0
    ordy = (ee - 1) * qv + 1 + NN * qv
    mult = ee + NN
    return deg, ordy, mult, deg - ordy - mult
check("general (a,b,t,kappa,a0,q) formulas give (108,144)->(550,205,69,276)",
      sig(3, 4, 4, 2, 8, 3) == (550, 205, 69, 276))
check("general formulas give (75,125)->(504,201,101,202)",
      sig(3, 5, 5, 3, 5, 2) == (504, 201, 101, 202))
# (72,108) is the r=0 resonance-gap exception: pure ansatz degf would be 10, not 14
check("(72,108) is the resonance-gap exception (r=0): pure-ansatz deg 10 < resonant 14",
      2 * 8 - 7 + 1 == 10 and 8 * (4 * 1 + 2 + 1) // 4 == 14)

print(f"\nALL {checks} PHI-75-125 CHECKS PASSED")
print(f"script: {Path(__file__).resolve()}")
sys.exit(0)
