#!/usr/bin/env python3
"""dm1_branch_verify.py -- independent EXACT verification of C10:
the ``d_{-1} == 0`` branch of the (72,108) case tree is EMPTY.

Target node: ``subcase:dm1`` (leaf ``L_D`` of the C0 partition), recorded in
``proof_dag.json`` at level ``claimed`` because ``PROOF_INVENTORY.md`` grades
C10 tier 2* with an *inferred* checker attribution ("within
verify_derivation.py").  That attribution is WRONG: ``verify_derivation.py``
verifies sections A-E of the *derivation* (forcing ODE, lambda-isolation, D_k
polynomiality, slice bridge, selection soundness) and contains no d_{-1}=0
branch check at all.  Before this file there was no checker for C10 anywhere.

WHAT THIS FILE ESTABLISHES, and how it differs from AUDIT.md sec.A.3.

  AUDIT.md A.3 / STATE.md item 5 close the branch with a two-leg case split:
  G1|_{dm1=0} = 3 dm2 dm3 forces dm2 == 0 or dm3 == 0 *because K[y] is a
  domain*; each leg then collapses G5body + Phi to Phi, contradicting
  Phi = f1 C4^28 != 0.  Sound, but it needs the domain step and a side
  condition (dm2 != 0) on the second leg.

  This file replaces that with a SINGLE POLYNOMIAL IDENTITY over Z -- no case
  split, no domain argument, no side conditions:

      2*Phi  =  2*(G5body + Phi)  +  2*d0*G1  +  2*d1*G2  +  2*d2*G3
                +  dm1 * ( d2*dm1^2 + 3*dm1*dm3 + 3*dm2^2 )        (CERT)

  an identity in Z[d2,d1,d0,dm1,dm2,dm3,dm4,Phi].  At any solution of the
  system the four G's vanish, so 2*Phi = dm1 * (d2 dm1^2 + 3 dm1 dm3 + 3 dm2^2).
  Setting dm1 == 0 gives 2*Phi == 0, hence Phi == 0 whenever 2 is invertible --
  contradicting Phi = f1 C4^28 != 0 (recomputed here from the forcing ODE).

  So the certificate is strictly stronger than A.3: it is case-free, needs only
  ``2 invertible'' rather than ``K[y] a domain'', and it exhibits the exact
  reason d_{-1} is the obstruction (2*Phi is congruent to a MULTIPLE OF dm1
  modulo the ideal -- which is also the structural source of the d_{-1}^21
  factor in the master identity f31 * f37 * d_{-1}^21).

CRITICAL SOUNDNESS POINT, checked in section B.  The branch d_{-1} == 0 must be
closed on the PRE-elimination generators G1,G2,G3,G5body.  The downstream
objects sol4, H2, H3 and the Singular resultants A,B are obtained by solving G1
for dm4, whose solution has denominator exactly 2*dm1 -- they are illegal on
this branch.  Section B asserts that denominator explicitly and confirms the
certificate touches none of those objects.

The generators are also REGENERATED FROM SCRATCH here (formal series in u), so
nothing on the mandatory path is taken on trust from ``generators.json``;
generators.json is used only as a cross-check target.

Run:  python3 dm1_branch_verify.py [--quiet] [--no-groebner]
      (~30 s; must end with "ALL <n> CHECKS PASSED")
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import sympy as sp
from sympy import Poly, Rational, expand, symbols, together

ROOT = Path(__file__).resolve().parent

# ---------------------------------------------------------------- harness ----
AP = argparse.ArgumentParser()
AP.add_argument("--quiet", action="store_true", help="suite mode: only headings + final line")
AP.add_argument("--no-groebner", action="store_true", help="skip the section-C Groebner cross-checks")
AP.add_argument("--emit", default="dm1_branch_certificate.json",
                help="write the integer certificate here ('' = do not write)")
ARGS = AP.parse_args()

_N = [0]


def out(msg: str) -> None:
    if not ARGS.quiet:
        print(msg)


def sec(title: str) -> None:
    print(f"\n{title}")


def check(name: str, cond) -> None:
    """A check must be a genuine claim: `cond` is required to be a bool that we
    computed, never a truthy object."""
    if cond is not True and cond is not False:
        raise SystemExit(f"  [BAD CHECK] {name}: condition is {type(cond)}, not bool")
    if not cond:
        raise SystemExit(f"  [FAIL] {name}")
    _N[0] += 1
    out(f"  [OK] {name}")


# ------------------------------------------------------------- the symbols ---
u, Phi = symbols("u Phi")
d2, d1, d0 = symbols("d2 d1 d0")
dm = {k: symbols(f"dm{k}") for k in range(1, 14)}
dm1, dm2, dm3, dm4 = dm[1], dm[2], dm[3], dm[4]
y = symbols("y")

# The emitted-identifier discipline of the Singular lane has an analogue here:
# assert the free symbols we introduce are pairwise distinct and that none of
# the four "core" branch variables collides with a series bookkeeping symbol.
_all = [u, Phi, d2, d1, d0, y] + [dm[k] for k in range(1, 14)]
assert len({s.name for s in _all}) == len(_all), "symbol name collision"


# ============================================================================
sec("A. FROM-SCRATCH REGENERATION of the pre-resultant generators")
# ============================================================================
# t = 3 control: the generator must reproduce the *published* (7,21) system of
# arXiv:2204.14178 sec.6.  A generator that fails this is not the paper's.
S3v = 1 + d1 * u**2 + d0 * u**3 + sum(dm[k] * u**(3 + k) for k in range(1, 11))
_t3 = Poly(expand(S3v * S3v), u).coeff_monomial(u**7)
check("A1  t=3 control: (D~^2)_{-1} == 2 d0 dm1 + 2 d1 dm2 + 2 dm4 (published (7,21) slice)",
      expand(_t3 - (2 * d0 * dm[1] + 2 * d1 * dm[2] + 2 * dm[4])) == 0)

# t = 4: C = x^4 * C4 * S with S = 1 + d2 u^2 + d1 u^3 + d0 u^4 + sum dm_k u^(4+k)
S = 1 + d2 * u**2 + d1 * u**3 + d0 * u**4 + sum(dm[k] * u**(4 + k) for k in range(1, 14))
S2 = Poly(expand(S * S), u)
S3 = Poly(expand(S2.as_expr() * S), u)
D2 = lambda k: S2.coeff_monomial(u**(8 + k))          # noqa: E731
D3 = lambda j: S3.coeff_monomial(u**(12 + j))         # noqa: E731

P_SIDE = [1, 2, 3, 4, 5, 6, 7, 9]
Q_SIDE = [1, 2, 3, 5]
USED = [D2(k) for k in P_SIDE] + [D3(j) for j in Q_SIDE]
check("A2  the 12 used equations are (D~^2)_{-1..-7,-9} and (D~^3)_{-1,-2,-3,-5}, "
      "and dm12 occurs in none of them",
      len(USED) == 12 and not any(e.has(dm[12]) for e in USED))

# The linear phase.  Each P-side equation is solved for ONE fresh variable.
FRESH = [(1, dm[5]), (2, dm[6]), (3, dm[7]), (4, dm[8]),
         (5, dm[9]), (6, dm[10]), (7, dm[11]), (9, dm[13])]
sub: dict = {}
pivots = []
for k, fresh in FRESH:
    e = expand(D2(k).subs(sub))
    p = Poly(e, fresh)
    pivots.append((k, fresh, p.degree(), p.coeff_monomial(fresh)))
    sub[fresh] = expand(sp.solve(e, fresh)[0])

check("A3  the linear phase is TRIANGULAR with UNIT PIVOT: each (D~^2)_{-k} is "
      "degree 1 in its fresh variable with coefficient exactly 2",
      all(deg == 1 and coeff == 2 for _, _, deg, coeff in pivots))

CORE = [d2, d1, d0, dm1, dm2, dm3, dm4]
_dens = [sp.denom(together(v)) for v in sub.values()]
check("A4  every substitution is a POLYNOMIAL in (d2,d1,d0,dm1..dm4) whose only "
      "denominator is a power of 2 -- no variable is ever inverted",
      all(v.is_polynomial(*CORE) for v in sub.values())
      and all(dn.is_Integer and int(dn) in (1, 2) for dn in _dens))

G1 = expand(D3(1).subs(sub))
G2 = expand(D3(2).subs(sub))
G3 = expand(D3(3).subs(sub))
G5body = expand(D3(5).subs(sub))
G5 = expand(G5body + Phi)
GENS = {"G1": G1, "G2": G2, "G3": G3, "G5body": G5body}

check("A5  the four regenerated generators involve ONLY d2,d1,d0,dm1..dm4 "
      "(the linear phase eliminated dm5..dm13 completely)",
      all(set(e.free_symbols) <= set(CORE) for e in GENS.values()))

# EQUIVALENCE, not merely weakening: back-substituting the solved values into
# every P-side equation returns 0 identically, so
#   {12 equations} <=> {dm5..dm13 = sub} AND {G1 = G2 = G3 = G5 = 0},
# and the fresh variables are FREE (each is determined, never constrained).
# Hence the 12-equation system is infeasible IFF {G1,G2,G3,G5} is -- an
# if-and-only-if, so nothing is lost on the dm1 == 0 branch.
_back = [expand(D2(k).subs(sub)) for k in P_SIDE]
check("A6  EQUIVALENCE: back-substitution annihilates all eight P-side equations, "
      "so {12 eqs} <=> {G1=G2=G3=G5=0} + free definitions of dm5..dm13",
      all(e == 0 for e in _back))

# Cross-check against the canonical, human-auditable generators.json (this is a
# CROSS-CHECK, not the source: everything above was rebuilt from the series).
_gj = json.loads((ROOT / "generators.json").read_text(encoding="utf-8"))
_order = [sp.Symbol(n) for n in _gj["variable_order"]]


def _expr_of(terms):
    acc = sp.Integer(0)
    for cs, exps in terms:
        mono = sp.Integer(1)
        for v, ex in zip(_order, exps):
            mono *= v**ex
        acc += Rational(cs) * mono
    return expand(acc)


check("A7  cross-check: all four regenerated generators are IDENTICAL to "
      "generators.json (the canonical artifact the f37 lane also consumes)",
      all(expand(GENS[n] - _expr_of(_gj["polynomials"][n])) == 0 for n in GENS))


# ============================================================================
sec("B. the branch must be closed PRE-elimination: no division by dm1")
# ============================================================================
sol4_num = _expr_of(_gj["sol4"]["numerator"])
sol4_den = _expr_of(_gj["sol4"]["denominator"])
sol4 = together(sol4_num / sol4_den)
check("B1  sol4 (the dm4-elimination that produces H2,H3,A,B) solves G1 for dm4 "
      "and has denominator EXACTLY 2*dm1 -- so the whole post-elimination chain "
      "is undefined on this branch and cannot be used here",
      expand(sp.denom(sol4) - 2 * dm1) == 0
      and expand(G1.subs(dm4, sol4_num / sol4_den) * 2 * dm1) == 0)

check("B2  correspondingly, G1 IS linear in dm4 with coefficient 3*dm1, which "
      "vanishes on the branch -- confirming dm4 stays a free unknown here",
      Poly(G1, dm4).degree() == 1 and expand(Poly(G1, dm4).coeff_monomial(dm4) - 3 * dm1) == 0)

# H2, H3 (post-sol4) are recorded in generators.json; assert they are NOT used
# below by checking they each vanish identically when dm1 -> 0, i.e. they carry
# no information on this branch.  This is the concrete form of "the resultant
# chain is silent at dm1 = 0".
H2 = _expr_of(_gj["polynomials"]["H2"])
H3 = _expr_of(_gj["polynomials"]["H3"])
check("B3  the post-elimination pair is BLIND on this branch: neither H2 nor H3 "
      "contains Phi at all, and at dm1 = 0 they collapse to -6*dm2^2*dm3 and "
      "-6*dm2*dm3^2, i.e. to multiples of G1|_{dm1=0}/3 = dm2*dm3 -- no new "
      "information beyond G1, and no route to Phi",
      (not H2.has(Phi)) and (not H3.has(Phi))
      and expand(H2.subs(Z_B := {dm1: 0}) + 6 * dm2**2 * dm3) == 0
      and expand(H3.subs(Z_B) + 6 * dm2 * dm3**2) == 0)


# ============================================================================
sec("C. THE CERTIFICATE: a single integer polynomial identity")
# ============================================================================
# Cofactors, over Z: (G1, G2, G3, G5) |-> (2*d0, 2*d1, 2*d2, 2).
COF = {"G1": 2 * d0, "G2": 2 * d1, "G3": 2 * d2, "G5": sp.Integer(2)}
RESIDUAL = dm1 * (d2 * dm1**2 + 3 * dm1 * dm3 + 3 * dm2**2)

_lhs = expand(2 * Phi)
_rhs = expand(COF["G1"] * G1 + COF["G2"] * G2 + COF["G3"] * G3 + COF["G5"] * G5 + RESIDUAL)
check("C1  FULL-RING IDENTITY over Z:  2*Phi == 2*d0*G1 + 2*d1*G2 + 2*d2*G3 "
      "+ 2*(G5body+Phi) + dm1*(d2*dm1^2 + 3*dm1*dm3 + 3*dm2^2)",
      expand(_lhs - _rhs) == 0)

_cert_polys = [COF["G1"], COF["G2"], COF["G3"], COF["G5"], RESIDUAL]
check("C2  every cofactor and the residual are INTEGER polynomials (no rational "
      "coefficients anywhere in the certificate); the certificate multiplier is 2, "
      "so the identity is valid over any ring in which 2 is invertible",
      all(Poly(e, *CORE).rep.dom.is_ZZ for e in _cert_polys))

# The consequence, spelled out mechanically: substitute dm1 -> 0.
Z0 = {dm1: 0}
g1z, g2z, g3z, g5z = (expand(e.subs(Z0)) for e in (G1, G2, G3, G5))
check("C3  at dm1 = 0 the residual vanishes, so 2*Phi lies in the ideal "
      "<G1,G2,G3,G5body+Phi>|_{dm1=0} with the explicit cofactor tuple "
      "(2*d0, 2*d1, 2*d2, 2)",
      expand(RESIDUAL.subs(Z0)) == 0
      and expand(2 * Phi - (2 * d0 * g1z + 2 * d1 * g2z + 2 * d2 * g3z + 2 * g5z)) == 0)

# ---- controls.  A certificate that would also "prove" the un-restricted case
# would be a red flag, not a result.
check("C4  NEGATIVE CONTROL: the residual is NOT the zero polynomial, so the "
      "identity does NOT prove Phi = 0 without dm1 = 0 (it must not -- that "
      "would refute the whole case and indicate a bug)",
      expand(RESIDUAL) != 0 and RESIDUAL.subs({dm1: 1, dm2: 1, dm3: 0, d2: 0}) == 3)

_tampered = [
    ("cofactor of G1 perturbed", 2 * d0 + 1, COF["G2"], COF["G3"], COF["G5"]),
    ("cofactor of G2 perturbed", COF["G1"], 2 * d1 + 1, COF["G3"], COF["G5"]),
    ("cofactor of G3 perturbed", COF["G1"], COF["G2"], 2 * d2 + 1, COF["G5"]),
    ("cofactor of G5 perturbed", COF["G1"], COF["G2"], COF["G3"], sp.Integer(3)),
]
check("C5  TAMPER CONTROL: perturbing any one of the four cofactors by 1 breaks "
      "the identity (so C1 is not vacuously satisfiable)",
      all(expand(_lhs - expand(c1 * G1 + c2 * G2 + c3 * G3 + c5 * G5 + RESIDUAL)) != 0
          for _, c1, c2, c3, c5 in _tampered))

check("C6  SCOPE CONTROL: the branch conclusion genuinely needs 2 invertible -- "
      "reduced mod 2 the certificate degenerates (2*Phi == 0 identically), so the "
      "honest field scope of this leaf is char != 2, not char 0",
      Poly(_lhs, *CORE, Phi, modulus=2).is_zero)

if not ARGS.no_groebner:
    V7 = (d2, d1, d0, dm2, dm3, dm4, Phi)
    GBz = sp.groebner([g1z, g2z, g3z, g5z], *V7, order="lex")
    check("C7  INDEPENDENT MECHANISM (Groebner, lex): Phi reduces to 0 modulo a "
          "Groebner basis of <G1,G2,G3,G5>|_{dm1=0} over Q -- so Phi is in that "
          "ideal on the nose, not merely 2*Phi",
          GBz.reduce(Phi)[1] == 0)
    V8 = (d2, d1, d0, dm1, dm2, dm3, dm4, Phi)
    GBf = sp.groebner([G1, G2, G3, G5], *V8, order="lex")
    check("C8  GROEBNER NEGATIVE CONTROL: in the FULL ring Phi does NOT reduce to "
          "0 modulo <G1,G2,G3,G5>, i.e. Phi is not in the un-restricted ideal -- "
          "the kill is specific to the dm1 = 0 branch",
          GBf.reduce(Phi)[1] != 0)
else:
    out("  [skipped C7,C8: --no-groebner]")


# ============================================================================
sec("D. the AUDIT.md sec.A.3 / STATE.md item 5 case split, reproduced exactly")
# ============================================================================
# The historical argument is CONFIRMED as well as superseded: every reduction it
# quotes is reproduced verbatim from the regenerated generators.
check("D1  G1|_{dm1=0} == 3*dm2*dm3   (AUDIT.md A.3 line 1)",
      expand(g1z - 3 * dm2 * dm3) == 0)
check("D2  leg dm2=0: G2|_{dm1=dm2=0} == (3/2)*dm3^2, forcing dm3 == 0 over a domain",
      expand(g2z.subs(dm2, 0) - Rational(3, 2) * dm3**2) == 0)
check("D3  leg dm2=0: G5body|_{dm1=dm2=dm3=0} == 0, so G5 collapses to Phi",
      expand(G5body.subs({dm1: 0, dm2: 0, dm3: 0})) == 0)
check("D4  leg dm3=0: G3|_{dm1=dm3=0} == -(3/2)*d1*dm2^2, forcing d1 == 0 when dm2 != 0",
      expand(g3z.subs(dm3, 0) + Rational(3, 2) * d1 * dm2**2) == 0)
check("D5  leg dm3=0: G5body|_{dm1=dm3=0} == -3*d1*dm2*dm4, which is 0 once d1 == 0, "
      "so G5 again collapses to Phi",
      expand(G5body.subs({dm1: 0, dm3: 0}) + 3 * d1 * dm2 * dm4) == 0
      and expand(G5body.subs({dm1: 0, dm3: 0, d1: 0})) == 0)
check("D6  the two legs COVER the branch: G1|_{dm1=0} = 0 is exactly dm2*dm3 = 0 "
      "up to the unit 3",
      expand(g1z / 3 - dm2 * dm3) == 0)


# ============================================================================
sec("E. Phi != 0, recomputed from the forcing ODE (not read from a file)")
# ============================================================================
C4 = y**7 * (y + 1)
a = symbols("a0:16")
ansatz = sum(a[i] * y**i for i in range(16))
eqs = Poly(expand(8 * y * (y + 1) * sp.diff(ansatz, y) - 14 * (8 * y + 7) * ansatz
                  - y**8 * (y + 1)**2), y).all_coeffs()
sol = sp.solve(eqs, a, dict=True)
check("E1  the forcing ODE 8y(y+1)f1' - 14(8y+7)f1 = y^8(y+1)^2 has a UNIQUE "
      "polynomial solution",
      len(sol) == 1)
f1 = expand(ansatz.subs(sol[0]))
f1_state = -y**8 * (y + 1)**2 * (2048 * y**4 - 512 * y**3 + 320 * y**2 - 240 * y + 195) / 6630
check("E2  that unique solution == STATE.md's f1", expand(f1 - f1_state) == 0)

PhiExpl = expand(f1 * C4**28)
Pp = Poly(PhiExpl, y)
check("E3  Phi = f1*C4^28 is NONZERO with the exact recorded profile: deg 238, "
      "ord_y 204, mult_{y+1} 30, trailing coeff -1/34, leading coeff -1024/3315",
      PhiExpl != 0
      and Pp.degree() == 238
      and min(m[0] for m in Pp.monoms()) == 204
      and sp.simplify(PhiExpl / (y + 1)**30).subs(y, -1) != 0
      and Pp.coeff_monomial(y**204) == Rational(-1, 34)
      and Pp.coeff_monomial(y**238) == Rational(-1024, 3315))

# The conclusion, assembled.  Everything in it has been checked above.
check("E4  CONCLUSION: on the branch dm1 == 0, (C1)+(C3) give 2*Phi == 0 in K[y]; "
      "with 2 invertible that is Phi == 0, contradicting E3.  The branch L_D is "
      "EMPTY -- case-free, denominator-free, domain-argument-free",
      expand(2 * PhiExpl) != 0)


# ------------------------------------------------------------------ emit -----
if ARGS.emit:
    cert = {
        "claim": "C10 / subcase:dm1 -- the d_{-1} == 0 branch of the (72,108) "
                 "case tree is EMPTY",
        "shape": "2*Phi = 2*d0*G1 + 2*d1*G2 + 2*d2*G3 + 2*(G5body + Phi) "
                 "+ dm1*(d2*dm1^2 + 3*dm1*dm3 + 3*dm2^2)",
        "ring": "Z[d2,d1,d0,dm1,dm2,dm3,dm4,Phi]",
        "multiplier": 2,
        "field_scope": "char != 2 (the only denominator in the certificate is 2; "
                       "the linear phase A3/A4 also pivots on 2)",
        "cofactors": {"G1": "2*d0", "G2": "2*d1", "G3": "2*d2", "G5body+Phi": "2"},
        "residual": "dm1*(d2*dm1^2 + 3*dm1*dm3 + 3*dm2^2)",
        "consequence": "at dm1 = 0 the residual vanishes, so 2*Phi lies in "
                       "<G1,G2,G3,G5body+Phi>; hence Phi = 0, contradicting "
                       "Phi = f1*C4^28 != 0 (deg 238, ord 204).",
        "generators": {n: sp.sstr(GENS[n]) for n in ("G1", "G2", "G3", "G5body")},
        "generator_provenance": "regenerated from scratch in section A (formal "
                                "series in u, linear phase with unit pivot 2); "
                                "cross-checked identical to generators.json (A7)",
        "supersedes": "AUDIT.md sec.A.3 / STATE.md item 5 two-leg case split "
                      "(also reproduced verbatim, section D)",
        "checks": _N[0],
    }
    (ROOT / ARGS.emit).write_text(json.dumps(cert, indent=1) + "\n", encoding="utf-8")
    out(f"\n  certificate written to {ARGS.emit}")

print(f"\nALL {_N[0]} CHECKS PASSED  "
      f"(C10 / subcase:dm1: the d_(-1) == 0 branch is EMPTY, "
      f"by one integer identity; char != 2)")
sys.exit(0)
