#!/usr/bin/env python3
"""phi_75_125.py  (NEW, uncommitted, read-only over all audited artifacts)

Derive the forcing divisor Phi = f * C^N for the (75,125) case, following EXACTLY
the corner-144 template (CORNER_144_COMPARISON.md + corner144_verify.py) but for a
NON-adjacent reduced C-power pair.

Case data (family F_2, j=1 -- the rung above Moh's 75 case):
  corner        A_0 = (5,20)         (GGV5 line 1679; v11 = 25)
  final corner  A_1 = (7/5, 2)       (GGV5 line 1679: p/l=7/5, q=2, k=1)
  (m,n)         = (3,5)              => degrees (75,125)
  reduced pair  (a,b) = (3,5)         sorted (m,n); NON-adjacent, b-a = 2

The independent PASS/FAIL checker is phi_75_125_verify.py (re-derives the ODE
solution and every divisor invariant from scratch).  Everything below is exact.
"""
import sympy as sp

y = sp.symbols("y")

# ---------------------------------------------------------------------------
# 1. Corner geometry -> template parameters
# ---------------------------------------------------------------------------
a0, l, q = 5, 5, 2          # a0 = deg C ; l = Laurent denominator ; q = selected mult
t = l                       # ell(C) = x^t c   (Laurent map Y -> x^l y)
kappa = l - 2               # [P,Q] = x^kappa ; Jacobian of (x^-1, x^l y) is -x^(l-2)
a, b = 3, 5                 # P = C^a , Q = C^b + (commuting C-powers) + F
e = b - a + 1               # RHS exponent c^e (=3 here; =2 on the adjacent (8,28) corner)
r = a0 - q - 1              # residual degree of the leading form
s = kappa + 1 - a * t       # v(F)   (obstruction valuation)
N = a * (t * (a + b) - (kappa + 1)) - 2 * b   # Phi C-power  = clear - b

print("=" * 74)
print("(75,125)  F_2 j=1   corner (5,20) -> (7/5,2)   (m,n)=(3,5)  (a,b)=(3,5)")
print("=" * 74)
print(f"  t = l = {t}   kappa = l-2 = {kappa}   e = b-a+1 = {e}   residual deg r = {r}")
print(f"  s = v(F) = kappa+1-a*t = {s}   N = a[t(a+b)-(kappa+1)]-2b = {N}")

# ---------------------------------------------------------------------------
# 2. The forcing ODE from the parametric family
#      a { t c f' - [t(b-a)+kappa+1] c' f } = c^(b-a+1)
# ---------------------------------------------------------------------------
coef = t * (b - a) + kappa + 1
print(f"\nForcing ODE (family F):  {a*t} c f' - {a*coef} c' f = c^{e}")
print(f"   i.e.   15 y^2(y^3+1) f' - 42 (...)' f = [y^2(y^3+1)]^3")

# ---------------------------------------------------------------------------
# 3. Solve the ODE.  Local orders (rho = (e-1)mu + 1) force
#      f = A y^5 g^3 ,   c = y^2 g ,   deg g = a0-q = 3 ,  g(-1)=0 , g(0)!=0
#    and the y=infinity leading coefficient is resonant at deg f = 14.
# ---------------------------------------------------------------------------
gc = sp.symbols("g0:4")
A = sp.symbols("A")
g = sum(gc[i] * y**i for i in range(4))
c = y**q * g
f = A * y**((e - 1) * q + 1) * g**e         # = A y^5 g^3
reduction = sp.factor(a * t * c * sp.diff(f, y) - a * coef * sp.diff(c, y) * f - c**e)
print(f"\nAnsatz f = A y^5 g^3 collapses the ODE to:   3A(y g' - 3 g) = 1")
print(f"   (sympy) residual/(-y^6 g^3) = {sp.simplify(-reduction/(y**6*g**3))}")

# coefficient extraction:  g1 = g2 = 0 (forced), g3 free (resonant),
# g(-1)=0 => g0 = g3, monic normalization g3 = 1, then A = -1/(9 g0) = -1/9.
g_sol = y**3 + 1
H2 = sp.factor(g_sol / (y + 1))             # residual quadratic
A_sol = sp.Rational(-1, 9)
c_sol = y**q * g_sol
f_sol = A_sol * y**5 * g_sol**3

assert sp.expand(a*t*c_sol*sp.diff(f_sol, y) - a*coef*sp.diff(c_sol, y)*f_sol - c_sol**e) == 0
print(f"\n  g = y^3 + 1 = (y+1)(y^2 - y + 1)      H2 = {H2}")
print(f"  C = y^2 (y^3+1)                       deg C = {sp.degree(c_sol,y)} (= a0)")
print(f"  f = -(1/9) y^5 (y^3+1)^3              deg f = {sp.degree(f_sol,y)}")

# ---------------------------------------------------------------------------
# 4. Phi = f * C^N  and its divisor signature (deg, ord_y, mult_{y+1}, cofactor)
# ---------------------------------------------------------------------------
Phi = sp.expand(f_sol * c_sol**N)

def order(expr):
    return min(m[0] for m in sp.Poly(expr, y).monoms())

def mult(expr, fac):
    cnt, qq, d = 0, sp.Poly(expr, y), sp.Poly(fac, y)
    while True:
        qq, rem = sp.div(qq, d)
        if not rem.is_zero:
            return cnt
        cnt += 1

deg = sp.degree(Phi, y)
ordy = order(Phi)
m1 = mult(Phi, y + 1)
cof = deg - ordy - m1
print("\n" + "-" * 74)
print(f"  Phi = f * C^{N} = -(1/9) y^{ordy} (y^3+1)^{m1}")
print(f"      = -(1/9) y^{ordy} (y+1)^{m1} (y^2-y+1)^{m1}")
print(f"  SIGNATURE (deg, ord_y, mult_(y+1), cofactor) = "
      f"({deg}, {ordy}, {m1}, {cof})")
print("-" * 74)

# ---------------------------------------------------------------------------
# 5. Verdict on the POTENTIAL_PROBE a-ONLY prediction
# ---------------------------------------------------------------------------
pred_deg = 64 * a**2 - 8 * a - 2      # 550
pred_mult = 8 * a**2 - a             # 69
print(f"\na-ONLY prediction (POTENTIAL_PROBE.md, a=3):  deg Phi = {pred_deg},"
      f"  mult_(y+1) = {pred_mult}")
print(f"derived (75,125)                          :  deg Phi = {deg},"
      f"  mult_(y+1) = {m1}")
print(f"VERDICT: DIFFERS."
      f"  deg {deg} != {pred_deg},  mult {m1} != {pred_mult}.")
print(f"  Unconditional: deg Phi = deg f + N*deg C = 14 + {N}*{a0} = {14+N*a0};")
print(f"  reaching 550 needs N=(550-14)/{a0}={sp.Rational(550-14,a0)} (non-integer) -- "
      f"impossible for a0={a0}.")
print(f"  The corner data (a0: 8->5, q=2, e=3) drives the difference, NOT a alone.")

if __name__ == "__main__":
    pass
