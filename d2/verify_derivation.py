"""verify_derivation.py — T6 audit: independent verification of the system derivation.

Verifies, by direct symbolic computation, every mechanically-checkable step
between "P,Q polynomial, [P,Q] = x^2, Prop-4.3 polygons" and the 12-equation
system regenerate_system.py feeds to Singular. Companion to T3_WINDOW_AUDIT.md
(window bounds) and T6_SELECTION_AUDIT.md (prose + remaining outline premises).

Setting: C in K[y,C4^-1]((x^-1)) with C^2 = P, leading term x^4*C4,
C4 = y^7(y+1); Q = C^3 + lambda*C^-1 + F with v_{1,0}(F) = -5. Writing
C = x^4*(C4 + c3 u + c2 u^2 + ...) with u = x^-1 and c_k the coefficient of
x^k, the polynomial d-variables are d_k := c_k * C4^(7-2k) (d4 = 1), and the
d3-killing shift x -> x - D3/4 sets c3 = 0.

Checks:
  A. the forcing ODE for f1 := C4^3 F_{-5} falls out of the commutator route
     [P,Q^2] = 2Q x^2 (paper's template, lines 1555-1596 of 2204.14178 src);
     its unique polynomial solution is STATE.md's f1; Phi = f1 C4^28 has the
     exact stats jetlift assumes, and jetlift's hard-coded Psi = Phi/y^204.
  B. lambda-isolation: (C^-1)_{-4} = C4^-1 (a unit) and (C^-1)_{-5} = -c3/C4^2,
     which vanishes after the shift. So lambda enters ONLY the dropped j=4
     slice, and the j=5 slice is (C^3)_{-5} + F_{-5} = 0 exactly.
  C. D_k := c_k C4^(7-2k) satisfies D_k = 1/2 P_{k+4} C4^(6-2k) - 1/2 sum D_i D_j
     (over i+j = k+4, i,j < 4): all C4 exponents cancel, so D_k in K[y] by
     induction on k downward (mirror of the paper's D_k-polynomiality prop).
  D. slice bridge: (C^2)_{-k} * C4^(14+2k) and (C^3)_{-j} * C4^(21+2j), under
     c_k = d_k C4^(2k-7), equal regenerate_system.py's D2(k), D3(j) EXACTLY;
     the F-term clearing at j=5 is F_{-5} C4^31 = f1 C4^28 = Phi.
  E. selection soundness: dm12 appears only in the dropped equations D2(8) and
     D3(4); slices beyond the used range each introduce a fresh unknown
     (dm14, dm15, ...) linearly, so truncation drops only always-satisfiable
     definitions. Dropping equations weakens a necessary condition => sound
     for the infeasibility direction.

Run:  python3 verify_derivation.py     (~1 min; must end with ALL ... PASSED)
"""
import sympy as sp
from sympy import symbols, Function, expand, together, cancel, simplify, Poly, Rational

y, x, u, lam = symbols('y x u lambda')
C4 = y**7*(y + 1)          # explicit, for section A
C4s = symbols('C4s')       # formal, for exponent bookkeeping in B-D

ok = [0]
def check(name, cond):
    if not cond:
        raise SystemExit(f"  [FAIL] {name}")
    ok[0] += 1
    print(f"  [OK] {name}")

# ---------------------------------------------------------------- A. the ODE
print("A. forcing ODE via the commutator route")
f1 = Function('f1')(y)
br = lambda g, h: sp.diff(g, x)*sp.diff(h, y) - sp.diff(g, y)*sp.diff(h, x)
# Q^2 - P^3 - 2*lambda*P = 2C^3F + (lower order); leading forms:
#   ell(P) = x^8 C4^2,  ell(Q) = x^12 C4^3,  ell(2C^3F) = 2x^7 C4^3 F_{-5} = 2x^7 f1.
# [P, Q^2 - P^3 - 2 lambda P] = [P,Q^2] = 2Q[P,Q] = 2Q x^2, leading 2 x^14 C4^3.
lhs = br(x**8*C4**2, 2*x**7*f1)
rhs = 2*x**2*x**12*C4**3
ode = cancel((lhs - rhs)/(2*x**14*C4*y**6))
state_ode = 8*y*(y+1)*sp.diff(f1, y) - 14*(8*y+7)*f1 - y**8*(y+1)**2
check("bracket identity reduces exactly to STATE.md's ODE",
      expand(ode - expand(state_ode)) == 0)

a = symbols('a0:16')
ansatz = sum(a[i]*y**i for i in range(16))
eqs = Poly(expand(8*y*(y+1)*sp.diff(ansatz, y) - 14*(8*y+7)*ansatz
                  - y**8*(y+1)**2), y).all_coeffs()
sol = sp.solve(eqs, a, dict=True)
check("ODE has a unique polynomial solution (deg <= 15 ansatz)", len(sol) == 1)
f1_state = -y**8*(y+1)**2*(2048*y**4 - 512*y**3 + 320*y**2 - 240*y + 195)/6630
check("unique solution == STATE.md's f1",
      expand(ansatz.subs(sol[0]) - f1_state) == 0)
quartic = 2048*y**4 - 512*y**3 + 320*y**2 - 240*y + 195
check("quartic separable; y,(y+1) do not divide it",
      sp.discriminant(quartic, y) != 0
      and quartic.subs(y, 0) != 0 and quartic.subs(y, -1) != 0)

Phi = expand(f1_state*C4**28)
Pp = Poly(Phi, y)
check("Phi = f1*C4^28: deg 238, ord 204, (y+1)-mult 30, coeffs -1/34 & -1024/3315",
      Pp.degree() == 238
      and min(m[0] for m in Pp.monoms()) == 204
      and simplify(Phi/(y+1)**30).subs(y, -1) != 0
      and Pp.coeff_monomial(y**204) == Rational(-1, 34)
      and Pp.coeff_monomial(y**238) == Rational(-1024, 3315))
psi_jetlift = -(y+1)**30*(2048*y**4-512*y**3+320*y**2-240*y+195)/6630
check("jetlift's hard-coded Psi == Phi/y^204", expand(Phi/y**204 - psi_jetlift) == 0)

# ------------------------------------------------- B. lambda-isolation slices
print("B. lambda-isolation: slices of C^-1")
c = {k: symbols(f'c{k}') for k in range(-13, 4)}
NT = 8
unit = C4s + sum(c[3-i]*u**(i+1) for i in range(NT-1))
inv = sp.series(1/unit, u, 0, NT).removeO()
inv4, inv5 = inv.coeff(u, 0), inv.coeff(u, 1)
check("(C^-1)_{-4} = C4^-1: unit => the dropped j=4 slice just determines lambda",
      simplify(inv4 - 1/C4s) == 0)
check("(C^-1)_{-5} = -c3/C4^2: zero after the d3-killing shift => j=5 has no lambda",
      simplify(inv5 + c[3]/C4s**2) == 0)
check("cleared j=4 lambda-term: lambda*(C^-1)_{-4}*C4^29 = lambda*C4^28",
      simplify(inv4*C4s**29 - C4s**28) == 0)

# --------------------------------------------- C. D_k in K[y]: the recursion
print("C. D_k := c_k C4^(7-2k) polynomial recursion (exponent cancellation)")
for k in range(3, -14, -1):
    pairs = [(i, k+4-i) for i in range(k+1, 4) if k+4-i <= 3 and k+4-i >= i]
    def cget(i): return c[i] if i in c else 0
    Pk4_rest = sum((2 if i != j else 1)*cget(i)*cget(j) for i, j in pairs)
    Pk4 = 2*C4s*cget(k) + Pk4_rest                       # P_{k+4} = (C^2)_{k+4}
    DD = sum((2 if i != j else 1)*(cget(i)*C4s**(7-2*i))*(cget(j)*C4s**(7-2*j))
             for i, j in pairs)
    rec = Rational(1,2)*Pk4*C4s**(6-2*k) - Rational(1,2)*DD
    check(f"D_({k:+d}) = 1/2 P_({k+4:+d}) C4^{6-2*k} - 1/2 sum D_i D_j",
          expand(rec - cget(k)*C4s**(7-2*k)) == 0)

# ------------------------------------------------------- D. the slice bridge
print("D. slice bridge: cleared C-side slices == regenerate_system.py's D2/D3")
d2s, d1s, d0s = symbols('d2 d1 d0')
dm = {k: symbols(f'dm{k}') for k in range(1, 14)}
S = 1 + d2s*u**2 + d1s*u**3 + d0s*u**4 + sum(dm[k]*u**(4+k) for k in range(1, 14))
S2 = Poly(expand(S*S), u); S3 = Poly(expand(S2.as_expr()*S), u)
D2 = lambda k: S2.coeff_monomial(u**(8+k))
D3 = lambda j: S3.coeff_monomial(u**(12+j))

cc = {4: C4s, 3: 0, 2: c[2], 1: c[1], 0: c[0]}
cc.update({-k: c[-k] for k in range(1, 14)})
U = sum(cc[4-i]*u**i for i in range(0, 18))              # C = x^4 * U
U2 = Poly(expand(U*U), u); U3 = Poly(expand(U2.as_expr()*U), u)
# d_k = c_k C4^(7-2k)  <=>  c_k = d_k C4^(2k-7)   (c4 = C4 <-> d4 = 1)
subs_cd = {c[2]: d2s*C4s**(-3), c[1]: d1s*C4s**(-5), c[0]: d0s*C4s**(-7)}
subs_cd.update({c[-k]: dm[k]*C4s**(-2*k-7) for k in range(1, 14)})

for k in range(1, 10):
    lhs = expand(cancel(together(U2.coeff_monomial(u**(8+k)).subs(subs_cd))
                        * C4s**(14+2*k)))
    check(f"(C^2)_(-{k}) * C4^{14+2*k} == D2({k})", expand(lhs - D2(k)) == 0)
for j in [1, 2, 3, 4, 5]:
    lhs = expand(cancel(together(U3.coeff_monomial(u**(12+j)).subs(subs_cd))
                        * C4s**(21+2*j)))
    check(f"(C^3)_(-{j}) * C4^{21+2*j} == D3({j})", expand(lhs - D3(j)) == 0)
check("F-term clearing at j=5: F_{-5} C4^31 = f1 C4^28 = Phi   (f1 = C4^3 F_{-5})",
      expand((f1_state/C4**3)*C4**31 - f1_state*C4**28) == 0)

# ------------------------------------------------- E. selection soundness
print("E. equation-selection soundness")
used = [D2(k) for k in [1,2,3,4,5,6,7,9]] + [D3(1), D3(2), D3(3), D3(5)]
check("dm12 absent from all 12 used equations",
      not any(e.has(dm[12]) for e in used))
check("dm12 present (linearly) in the dropped D2(8) and D3(4)",
      D2(8).has(dm[12]) and D3(4).has(dm[12])
      and D2(8).coeff(dm[12]) == 2 and D3(4).coeff(dm[12]) == 3)
dmx = {k: symbols(f'dm{k}') for k in range(14, 18)}
Sx = S + sum(dmx[k]*u**(4+k) for k in range(14, 18))
S2x = Poly(expand(Sx*Sx), u); S3x = Poly(expand(S2x.as_expr()*Sx), u)
for k in [10, 11, 12]:
    check(f"(D^2)_(-{k}) introduces fresh dm{k+4} linearly (coeff 2)",
          S2x.coeff_monomial(u**(8+k)).coeff(dmx[k+4]) == 2)
for j in [6, 7]:
    check(f"(D^3)_(-{j}) introduces fresh dm{j+8} linearly (coeff 3)",
          S3x.coeff_monomial(u**(12+j)).coeff(dmx[j+8]) == 3)

print(f"\nALL {ok[0]} DERIVATION CHECKS PASSED")
