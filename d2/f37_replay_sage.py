#!/usr/bin/env sage
"""Independent Sage replay of the f37 pre-resultant elimination theorem.

This script deliberately does not read generators.json, factor_*.txt,
f31_graded.txt, a pickle, a certificate, or generated CAS input. It rebuilds
the four generators by formal-series convolution and separately rebuilds f31
through the classical resultant chain.
"""
from sage.all import PolynomialRing, QQ


def report(label, ok):
    print(("PASS" if ok else "FAIL") + ": " + label, flush=True)
    return bool(ok)


# ---------------------------------------------------------------------------
# 1. Build the pre-resultant generators from the mathematical definition.
# ---------------------------------------------------------------------------
source_names = ["d2", "d1", "d0"] + [f"dm{k}" for k in range(1, 14)] + ["Phi"]
A = PolynomialRing(QQ, names=source_names, order="degrevlex")
(
    d2A,
    d1A,
    d0A,
    dm1A,
    dm2A,
    dm3A,
    dm4A,
    dm5A,
    dm6A,
    dm7A,
    dm8A,
    dm9A,
    dm10A,
    dm11A,
    dm12A,
    dm13A,
    PhiA,
) = A.gens()
dmA = {
    1: dm1A,
    2: dm2A,
    3: dm3A,
    4: dm4A,
    5: dm5A,
    6: dm6A,
    7: dm7A,
    8: dm8A,
    9: dm9A,
    10: dm10A,
    11: dm11A,
    12: dm12A,
    13: dm13A,
}


def s_coeff(n):
    """Coefficient of u^n in S = 1+d2*u^2+d1*u^3+d0*u^4+... ."""
    if n == 0:
        return A.one()
    if n == 2:
        return d2A
    if n == 3:
        return d1A
    if n == 4:
        return d0A
    if 5 <= n <= 17:
        return dmA[n - 4]
    return A.zero()


def coeff_square(n):
    return A(sum(s_coeff(i) * s_coeff(n - i) for i in range(n + 1)))


def coeff_cube(n):
    return A(
        sum(
            s_coeff(i) * s_coeff(j) * s_coeff(n - i - j)
            for i in range(n + 1)
            for j in range(n - i + 1)
        )
    )


# Published t=3 checkpoint: [u^7]S3^2 has this form. It independently checks
# the indexing convention before the t=4 construction is used.
def t3_coeff(n):
    if n == 0:
        return A.one()
    if n == 2:
        return d1A
    if n == 3:
        return d0A
    if 4 <= n <= 13:
        return dmA[n - 3]
    return A.zero()


t3_u7 = A(sum(t3_coeff(i) * t3_coeff(7 - i) for i in range(8)))
t3_expected = 2 * d0A * dm1A + 2 * d1A * dm2A + 2 * dm4A
t3_ok = t3_u7 == t3_expected

# D2(k)=[u^(8+k)]S^2 is linear in the displayed fresh variable. Solve the
# eight equations successively, exactly as the mathematical linear phase says.
substitutions = {}
linear_steps = (
    (1, dm5A),
    (2, dm6A),
    (3, dm7A),
    (4, dm8A),
    (5, dm9A),
    (6, dm10A),
    (7, dm11A),
    (9, dm13A),
)
linear_ok = True
for k, target in linear_steps:
    equation = A(coeff_square(8 + k).subs(substitutions))
    coefficient = equation.derivative(target)
    linear_ok = linear_ok and coefficient != 0
    linear_ok = linear_ok and coefficient.derivative(target) == 0
    linear_ok = linear_ok and equation.derivative(target, 2) == 0
    substitutions[target] = A(-equation.subs({target: 0}) / coefficient)

G1A, G2A, G3A, G5bodyA = [
    A(coeff_cube(12 + j).subs(substitutions)) for j in (1, 2, 3, 5)
]

half3 = QQ(3) / QQ(2)
expected = (
    half3 * d1A * dm1A**2
    + 3 * d2A * dm1A * dm2A
    + 3 * dm1A * dm4A
    + 3 * dm2A * dm3A,
    -half3 * d0A * dm1A**2
    + half3 * d2A * dm2A**2
    + 3 * dm2A * dm4A
    + half3 * dm3A**2,
    -3 * d0A * dm1A * dm2A
    - half3 * d1A * dm2A**2
    - QQ(1) / QQ(2) * dm1A**3
    + 3 * dm3A * dm4A,
    -3 * d0A * dm1A * dm4A
    - 3 * d0A * dm2A * dm3A
    - 3 * d1A * dm2A * dm4A
    - half3 * d1A * dm3A**2
    - 3 * d2A * dm3A * dm4A
    - half3 * dm1A**2 * dm3A
    - half3 * dm1A * dm2A**2,
)
construction_ok = t3_ok and linear_ok and all(
    got == want for got, want in zip((G1A, G2A, G3A, G5bodyA), expected)
)
report("formal-series construction of G1,G2,G3,G5body", construction_ok)

# Move only the retained variables into a lex ring with the three variables to
# eliminate first. The explicit homomorphism sends discarded source variables
# to zero; they have already disappeared from the four generated polynomials.
R = PolynomialRing(
    QQ,
    names=("dm2", "dm3", "dm4", "d2", "d1", "d0", "dm1", "Phi"),
    order="lex",
)
dm2, dm3, dm4, d2, d1, d0, dm1, Phi = R.gens()
to_R = A.hom(
    [
        d2,
        d1,
        d0,
        dm1,
        dm2,
        dm3,
        dm4,
        R.zero(),
        R.zero(),
        R.zero(),
        R.zero(),
        R.zero(),
        R.zero(),
        R.zero(),
        R.zero(),
        R.zero(),
        Phi,
    ],
    R,
)
G1, G2, G3, G5body = map(to_R, (G1A, G2A, G3A, G5bodyA))
G5 = G5body + Phi
I = R.ideal([G1, G2, G3, G5])

# ---------------------------------------------------------------------------
# 2. Independently regenerate f31 by the historical resultant route.
# ---------------------------------------------------------------------------
# These combinations are the denominator-cleared substitution solving G1 for
# dm4. They eliminate dm4 without importing sol4 or any serialized H-system.
H2 = dm1 * G2 - dm2 * G1
H3 = dm1 * G3 - dm3 * G1
H5 = dm1 * G5 + (d0 * dm1 + d1 * dm2 + d2 * dm3) * G1
assert all(h.degree(dm4) == 0 for h in (H2, H3, H5))

print("INFO: computing the two dm3 resultants", flush=True)
RA = H2.resultant(H3, dm3)
RB = H2.resultant(H5, dm3)
Ah_candidates = [q for q, exponent in RA.factor() if q.degree(dm2) > 0]
Bh_candidates = [q for q, exponent in RB.factor() if q.degree(dm2) > 0]
resultant_inputs_ok = len(Ah_candidates) == 1 and len(Bh_candidates) == 1
report("unique dm2-bearing factors in the two intermediate resultants", resultant_inputs_ok)
if not resultant_inputs_ok:
    raise SystemExit(1)
Ah = Ah_candidates[0]
Bh = Bh_candidates[0]

print("INFO: computing and factoring the final dm2 resultant", flush=True)
master = Ah.resultant(Bh, dm2)
f31_candidates = [
    q
    for q, exponent in master.factor()
    if q.total_degree() == 31 and q.degree(Phi) > 0
]
f31_unique = len(f31_candidates) == 1
report("unique Phi-bearing total-degree-31 resultant factor", f31_unique)
if not f31_unique:
    raise SystemExit(1)
f31 = f31_candidates[0]

# ---------------------------------------------------------------------------
# 3. Compute the true elimination ideal and replay the theorem.
# ---------------------------------------------------------------------------
print("INFO: computing I intersect QQ[d2,d1,d0,dm1,Phi]", flush=True)
E = I.elimination_ideal([dm2, dm3, dm4])
E_basis = list(E.groebner_basis())
principal_ok = len(E_basis) == 1 and E_basis[0] != 0
e_generator = E_basis[0] if principal_ok else R.zero()
match_ok = principal_ok and R.ideal([e_generator]) == R.ideal([f31])
membership_ok = f31.reduce(I.groebner_basis()) == 0

ok_a = report(
    "(a) elimination ideal in (d2,d1,d0,dm1,Phi) is principal",
    principal_ok,
)
ok_b = report(
    "(b) its generator equals independently reconstructed f31 up to a QQ unit",
    match_ok,
)
ok_c = report(
    "(c) independently reconstructed f31 lies in the pre-resultant ideal",
    membership_ok,
)

print(
    "INFO: elimination generators={}, f31 total degree={}, f31 terms={}".format(
        len(E_basis), f31.total_degree(), len(f31.dict())
    ),
    flush=True,
)
if not (construction_ok and resultant_inputs_ok and f31_unique and ok_a and ok_b and ok_c):
    raise SystemExit(1)
print("PASS: f37 theorem replay complete over QQ", flush=True)
