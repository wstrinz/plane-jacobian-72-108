"""t6_premises_verify.py — finite checks backing T6_PREMISES.md.

The two remaining outline-only premises of the (72,108) derivation are the
transcription, from t=3 to t=4, of the machinery GGHV22 (arXiv:2204.14178)
uses for its *twin* closed case (9,27):

  Premise 1  ell_{1,0}(P) = R^2, ell_{1,0}(Q) = R^3, R = x^4 C4, C4 = y^7(y+1)
             (GGHV22 lines 1411-1414, via GGV1 Props 1.13 + 2.1).
  Premise 2  the alpha-strip WLOG reducing Q to C^3 + lambda C^-1 + F, v(F) = -5
             (GGHV22 lines 1508-1546).

Both premises are structural inputs (published GGV1 propositions) plus finite
arithmetic. This script verifies EVERY finite piece:

  P1a  the valuation gap 2 < v(P)+v(Q)-1 that triggers GGV1 Prop 1.13
  P1b  the alignment exponents (m,n)=(2,3): n*v(P) = m*v(Q)
  P1c  corner arithmetic: R^2, R^3 reproduce the subcase-2 leading corners,
       and force C4 = y^7(a0+a1 y) with a0*a1 != 0
  P1d  primitivity: x^4 y^7 (y+1) is not S^d for any d>=2  (=> R primitive,
       needed for Prop 2.1 to give a single power alpha_k R^k)
  P1e  the linear change of variables normalizing y^7(a0+a1 y) to y^7(y+1)
  P2a  powers of C Poisson-commute with P=C^2  (so stripping alpha_k C^k
       leaves [.,P] unchanged, and [F,P]=[Q,P]=-x^2)
  P2b  v(F) = -5 from the Prop-1.13 equality  2 = v(F)+v(P)-1
  P2c  the strip range k in (-2,3) = {-1,0,1,2}, terminating because
       v(F)=-5 is not a multiple of v_{1,0}(R)=4
  P2d  the Remark shift: P~ = P + (2/3)alpha1 gives C~ = C + (1/3)alpha1 C^-1
       + ... and C~^3 = C^3 + alpha1 C + ..., absorbing the alpha1 C term
  P2e  bracket preservation: [P~, Q~] = [P,Q] = x^2

Run:  python3 t6_premises_verify.py     (ends with ALL ... PASSED)
"""
import math
import sympy as sp
from sympy import symbols, expand, Rational, Poly, gcd

x, y, e, t = symbols('x y e t')
a0, a1, alpha1 = symbols('a0 a1 alpha1')

ok = [0]
def check(name, cond):
    if not cond:
        raise SystemExit(f"  [FAIL] {name}")
    ok[0] += 1
    print(f"  [OK] {name}")

# =============================================================== Premise 1
print("Premise 1: ell(P)=R^2, ell(Q)=R^3, R=x^4 C4, C4=y^7(y+1)")

# --- our (1,0)-valuations, read off the subcase-2 corners (STATE.md item, 4.3):
#   N(P) top-corners (8,14),(8,16);  N(Q) top-corners (12,21),(12,24).
#   v_{1,0} = rho*i + sigma*j with (rho,sigma)=(1,0) is just the x-exponent.
vP  = 8      # v_{1,0}(P): x-exponent of the leading corner (8,16)/(8,14)
vQ  = 12     # v_{1,0}(Q): x-exponent of the leading corner (12,24)/(12,21)
vC  = 4      # v_{1,0}(C) = v_{1,0}(R),  R = x^4 C4
vXX = 2      # v_{1,0}([P,Q]) = v_{1,0}(x^2)

# P1a: the strict gap that makes GGV1 Prop 1.13 force [ell P, ell Q] = 0.
#   Prop 1.13:  v([P,Q]) <= v(P)+v(Q)-(rho+sigma),  equality  <=>  [ellP,ellQ]!=0.
#   (rho,sigma)=(1,0) => rho+sigma = 1.
check("P1a: 2 = v([P,Q]) < v(P)+v(Q)-1 = 19  (=> [ellP,ellQ]=0 by Prop 1.13)",
      vXX < vP + vQ - 1 and vP + vQ - 1 == 19 and vXX == 2)

# P1b: alignment exponents.  Prop 2.1 with tau=v(P), mu=v(Q): n*tau = m*mu,
#   gcd(m,n)=1.  Here 3*8 = 2*12, so (m,n) = (2,3): ell(P)=R^2, ell(Q)=R^3.
m, n = 2, 3
check("P1b: (m,n)=(2,3) solves n*v(P)=m*v(Q)  (3*8 = 2*12 = 24), gcd=1",
      n*vP == m*vQ and math.gcd(m, n) == 1)

# P1c: corner arithmetic.  R = x^4 C4, C4 = y^7 (a0 + a1 y) (general).
#   R^2 = x^8 C4^2, R^3 = x^12 C4^3; read off the y-support corners.
C4 = y**7*(a0 + a1*y)
def yspan(poly_in_y):
    P = Poly(sp.expand(poly_in_y), y)
    ms = [mm[0] for mm in P.monoms()]
    return min(ms), max(ms)
lo2, hi2 = yspan(C4**2)          # y-degrees of C4^2
lo3, hi3 = yspan(C4**3)          # y-degrees of C4^3
check("P1c: R^2 = x^8 C4^2 has y-corners (8,14),(8,16)  [matches N(P)]",
      (8, lo2) == (8, 14) and (8, hi2) == (8, 16))
check("P1c: R^3 = x^12 C4^3 has y-corners (12,21),(12,24)  [matches N(Q)]",
      (12, lo3) == (12, 21) and (12, hi3) == (12, 24))
# both distinct corners present  <=>  a0 != 0 (gives ord 14, not 16) and
# a1 != 0 (gives deg 16, not 14).  Degenerate checks:
lo2_a0, hi2_a0 = yspan((y**7*a1*y)**2)      # a0 -> 0
lo2_a1, hi2_a1 = yspan((y**7*a0)**2)        # a1 -> 0
check("P1c: a0=0 collapses (8,14) [ord 16 not 14]; a1=0 collapses (8,16) "
      "[deg 14 not 16]  => a0*a1 != 0 forced by the two corners",
      lo2_a0 == 16 and hi2_a1 == 14)

# P1d: primitivity of x^4 y^7 (y+1).  It is S^d (d>=2) only if d divides the
#   multiplicity of every irreducible factor: x(4), y(7), (y+1)(1).
g = math.gcd(math.gcd(4, 7), 1)
check("P1d: gcd of exponent multiplicities (4,7,1) = 1  => x^4 C4 is not a "
      "proper power  (needed for Prop 2.1's single-power conclusion)", g == 1)

# P1e: linear change of variables y = c*y~ with c = a0/a1 sends
#   y^7 (a0 + a1 y)  to  (const) * y~^7 (y~ + 1).
c = a0/a1
ytil = symbols('ytil')
C4_sub = sp.expand(C4.subs(y, c*ytil))
target = ytil**7*(ytil + 1)
ratio = sp.simplify(C4_sub/target)
check("P1e: y=(a0/a1) y~ turns y^7(a0+a1 y) into const * y~^7(y~+1)",
      ratio.free_symbols.isdisjoint({ytil}) and ratio != 0)

# =============================================================== Premise 2
print("Premise 2: alpha-strip WLOG, Q = C^3 + lambda C^-1 + F, v(F) = -5")

# P2a: powers of C Poisson-commute; in particular each stripped term
#   alpha_k C^k (k=2,1,0,-1) commutes with P = C^2, so [F,P] = [Q,P].
#   Verify [C^a, C^2] = 0 for a in {3,1,0,-1} with a generic C(x,y).
Cgen = x*y**2 + x**2 - 3*y + 7          # generic, non-symmetric
def br(f, h):
    return sp.diff(f, x)*sp.diff(h, y) - sp.diff(f, y)*sp.diff(h, x)
allzero = all(sp.simplify(br(Cgen**a, Cgen**2)) == 0 for a in (3, 1, 0, -1))
check("P2a: [C^k, C^2] = 0 for k in {3,1,0,-1}  => [F,P]=[Q,P]=-x^2", allzero)

# P2b: v(F) from Prop 1.13 equality.  [F,P] = [Q,P] = -x^2 (P2a), so
#   v([F,P]) = v(x^2) = 2; [ellF,ellP] != 0 by construction (strip terminates),
#   so equality gives 2 = v(F) + v(P) - 1, i.e. v(F) = -5.
vF = vXX - vP + 1
check("P2b: v(F) = v([F,P]) - v(P) + 1 = 2 - 8 + 1 = -5", vF == -5)

# P2c: strip range.  Strippable leading forms are the pure powers alpha_k R^k
#   with v_{1,0} = 4k, in the open range -4 <= 4k < v(C^3)=12, i.e. k in
#   {-1,0,1,2} = (-2,3).  The descent halts at F because v(F)=-5 is NOT a
#   multiple of v_{1,0}(R)=4 (so ell(F) is not any alpha_k R^k).
krange = [k for k in range(-3, 4) if -4 <= vC*k < vQ]
check("P2c: k in (-2,3) = {-1,0,1,2}: powers -4,0,4,8 all in [-4,12)",
      krange == [-1, 0, 1, 2] and [vC*k for k in krange] == [-4, 0, 4, 8])
check("P2c: v(F) = -5 is not a multiple of v_{1,0}(R)=4  => strip terminates",
      (-5) % 4 != 0)

# P2d: the Remark shift.  P~ = P + (2/3) alpha1  (alpha1 a constant), and
#   C~ := sqrt(P~) = sqrt(C^2 + (2/3)alpha1).  Expand in u = 1/C (C large):
#   C~ = C + (1/3) alpha1 C^-1 + ...,  and  C~^3 = C^3 + alpha1 C + ....
tshift = Rational(2, 3)*alpha1
# C = 1/e ; sqrt(C^2 + tshift) = (1/e) sqrt(1 + tshift e^2)
Ctil = sp.series((1/e)*sp.sqrt(1 + tshift*e**2), e, 0, 6).removeO()
coeff_Cm1 = Ctil.coeff(e, 1)                       # coefficient of C^-1 = e^1
check("P2d: C~ = C + (1/3)alpha1 C^-1 + ...  (P~ = P + (2/3)alpha1)",
      sp.simplify(coeff_Cm1 - Rational(1, 3)*alpha1) == 0)
Ctil3 = sp.series(Ctil**3, e, 0, 4).removeO()
coeff_C1 = Ctil3.coeff(e, -1)                      # coefficient of C^+1 = e^-1
coeff_C3 = Ctil3.coeff(e, -3)                      # coefficient of C^3  = e^-3
check("P2d: C~^3 = C^3 + alpha1 C + ...  (the alpha1 C term is absorbed)",
      sp.simplify(coeff_C1 - alpha1) == 0 and sp.simplify(coeff_C3 - 1) == 0)

# P2e: bracket preservation.  P~ = P + (2/3)alpha1, Q~ = Q - alpha2 P - alpha0
#   (alpha0, alpha2 constants).  [P~,Q~] = [P,Q] since brackets with constants
#   vanish and [P, alpha2 P] = 0.  Verify with generic P, Q.
alpha0, alpha2 = symbols('alpha0 alpha2')
Pg = x**2*y + 5*x - y**3
Qg = x*y - 2*x**3 + y
Ptil = Pg + Rational(2, 3)*alpha1
Qtil = Qg - alpha2*Pg - alpha0
check("P2e: [P~,Q~] = [P,Q]  (constants and alpha2 P drop out)",
      sp.expand(br(Ptil, Qtil) - br(Pg, Qg)) == 0)

print(f"\nALL {ok[0]} PREMISE CHECKS PASSED")
