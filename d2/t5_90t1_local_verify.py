#!/usr/bin/env python3
"""Exact local certificate for the nonconstant-E a_t=9 T1 cell."""

from __future__ import annotations

import sympy as sp

import t5_90t1_verify as base


# First link the local formulas to the repository source decomposition.
H = base.load_h()
d0, d1_source, d2_source, e_source = base.d0, base.d1, base.d2, base.e
sigma_source = 4*d0-d2_source**2
h4_source = -16*(
    756*sigma_source**3 + 324*sigma_source**2*d2_source**2
    + 1476*sigma_source*d1_source**2*d2_source
    - 2160*sigma_source*d1_source*e_source + 13797*d1_source**4
    + 1952*d1_source**2*d2_source**3
    + 192*d1_source*d2_source**2*e_source - 352*d2_source*e_source**2
)
assert sp.expand(H[4]-h4_source) == 0


x, T = sp.symbols("x T", nonzero=True)
c, gamma, delta = sp.symbols("c gamma delta", nonzero=True)
Q0 = sp.symbols("Q0", nonzero=True)
Q1, Q2, Q3, Q4 = sp.symbols("Q1 Q2 Q3 Q4")
W1, W2 = sp.symbols("W1 W2")
D0, D1, D2 = sp.symbols("D0 D1 D2")
S0, S1, S2 = sp.symbols("S0 S1 S2")

t = T+x
q = Q0+Q1*x+Q2*x**2+Q3*x**3+Q4*x**4
W0 = -gamma**4*T**6/(2*c*delta*Q0)
W = W0+W1*x+W2*x**2
D = D0+D1*x+D2*x**2
S = S0+S1*x+S2*x**2
eta = -8192*c**7*delta**2/gamma**3

# Level 6 defines G6=x^2*K after q^6 is removed.
h6_over_x5 = (
    -3072*x*S**2 + 14336*delta**2*x**3*W**2*D
    + 8192*delta*gamma*t**9*W
)
K = sp.expand((t**3*q*eta*W**2-c**6*h6_over_x5)/gamma**3)

# Level 5 has G5=N5/(gamma^3*x). Its existence plus level 4 force
# ord_x(N5)>=3.
h5_over_x2 = (
    -9216*x**4*D*S**2 + 32256*delta**2*x**9*W**2*S
    - 12288*delta**2*x**6*W**2*D**2
    + 18432*delta*gamma*x**3*W*D*t**9 + 2048*gamma**2*t**18
)
N5 = sp.expand(t**3*q*K-c**5*h5_over_x2)
R = sp.expand(gamma**4*t**6+2*c*delta*q*W)


def coefficient(expr: sp.Expr, degree: int) -> sp.Expr:
    return sp.Poly(expr, x).coeff_monomial(x**degree)


N51 = sp.factor(coefficient(N5, 1))
assert sp.factor(N51-3072*Q0*S0**2*T**3*c**6/gamma**3) == 0
R1 = sp.factor(coefficient(R, 1))
N52 = sp.factor(coefficient(N5, 2).subs(S0, 0))
assert sp.factor(N52+2048*T**6*c**5*R1**2/gamma**6) == 0

# R_1=0 fixes W1. Build level 4: G4=N4/(gamma^3*x^3).
W1_forced = sp.factor(gamma**4*T**5*(T*Q1-6*Q0)/(2*c*delta*Q0**2))
conditions = {S0: 0, W1: W1_forced}
sigma_local = x**3*S
d1_local = delta*x**4*W
e_local = gamma*x*t**9
h4 = -16*(
    756*sigma_local**3 + 324*sigma_local**2*D**2
    + 1476*sigma_local*d1_local**2*D - 2160*sigma_local*d1_local*e_local
    + 13797*d1_local**4 + 1952*d1_local**2*D**3
    + 192*d1_local*D**2*e_local - 352*D*e_local**2
)
A = sp.expand(t**3*q)


def n4_coefficient(degree: int) -> sp.Expr:
    value = sum(
        coefficient(A, i)*coefficient(N5, degree-i+1)/gamma**3
        for i in range(degree+1)
    ) - c**4*coefficient(h4, degree)
    return sp.factor(value.subs(conditions))


N42 = n4_coefficient(2)
assert sp.factor(N42-3072*Q0**2*S1**2*T**6*c**6/gamma**6) == 0
conditions[S1] = 0
R2 = sp.factor(coefficient(R, 2).subs(conditions))
N43 = n4_coefficient(3)
assert sp.factor(N43+2048*Q0*T**9*c**5*R2**2/gamma**9) == 0

# R_2=0 fixes the last coefficient W2. Level 3 forces ord_x(G4)>=2,
# hence the x^4 coefficient of N4 also vanishes.
W2_forced = sp.factor(
    -(15*gamma**4*T**4/(2*c*delta)+Q1*W1_forced+Q2*W0)/Q0
)
conditions[W2] = W2_forced
N44 = n4_coefficient(4)
assert sp.factor(N44-3072*Q0**2*S2**2*T**6*c**6/gamma**6) == 0

print("a_t=9 geometric q-coprime T1, nonconstant E: LOCAL REDUCTION PASS")
print("  level 5/4 alternation forces s^6 | sigma and s^3 | R")
print("  R=gamma^4*t^6+2*c*delta*q*W; all coefficients of quadratic W are fixed")
