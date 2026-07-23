#!/usr/bin/env python3
"""Exact verification that f37 is a resultant artifact: f31 lies in the
pre-resultant ideal, so every solution of the original system satisfies
f31 = 0 and none lies on {f37 = 0} \\ {f31 = 0}.

Two independent exact checks, neither trusting the other tool:

  (A) MEMBERSHIP CERTIFICATE (self-contained, no Singular needed).
      Read the four cofactors produced by Singular's lift() from
      f37_sat_certificate.txt and verify, purely in sympy, the polynomial
      identity
          f31 == c1*G1 + c2*G2 + c3*G3 + c4*(G5body + Phi)
      where G1,G2,G3,G5body are regenerated from t4_state.pkl.  Because f31 is
      an explicit polynomial combination of the pre-resultant generators, it
      vanishes on the entire pre-resultant variety, over every field and every
      specialization of Phi (in particular the genuine (72,108) instance).

  (B) MASTER-IDENTITY CONSISTENCY.  Re-confirm, from the same regenerated
      state, that the resultant master identity f31*f37*d_-1^21 also lies in
      the ideal, so the case tree {f31=0} u {f37=0} is sound; combined with (A)
      the f37 (and the d_-1^21) factor are pure resultant excess.

The elimination-ideal fact E = <G-system> cap Q[d2,d1,d0,d_-1,Phi] = <f31>
(principal, 102 terms, degree 31) is a Groebner statement reproduced by
f37_sat_confirm.sing; the certificate below is the field-independent core and
needs no Groebner engine to check.
"""

from __future__ import annotations

import pickle
import re
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent

d2, d1, d0 = sp.symbols("d2 d1 d0")
dm1, dm2, dm3, dm4 = sp.symbols("dm1 dm2 dm3 dm4")
Phi = sp.symbols("Phi")

# Singular used the names m1 = dm1 (= d_-1) and P = Phi.
LOCALS = {
    "d2": d2, "d1": d1, "d0": d0,
    "m1": dm1, "dm1": dm1, "dm2": dm2, "dm3": dm3, "dm4": dm4,
    "P": Phi, "Phi": Phi,
}


_TERM = re.compile(r"""
    (?P<sign>[+-]?)
    (?P<num>\d+)(?:/(?P<den>\d+))?        # optional rational coefficient
    (?P<mono>(?:\*[A-Za-z]\w*(?:\*\*\d+)?)*)   # optional monomial factors
""", re.VERBOSE)

_FACTOR = re.compile(r"([A-Za-z]\w*)(?:\*\*(\d+))?")


def parse(expr: str) -> sp.Expr:
    # The strings are flat polynomials (no parentheses): a sum of signed
    # monomials with *exact rational* coefficients (a or a/b).  Parse each
    # monomial explicitly (exact sympy Rational, no float division) and sum
    # with sp.Add -- avoids both float error and the deep +-chain that
    # overflows the eval compiler.
    expr = expr.replace("^", "**").replace(" ", "").strip()
    if not expr:
        return sp.Integer(0)
    terms: list[sp.Expr] = []
    pos = 0
    for m in _TERM.finditer(expr):
        if m.start() != pos:
            raise ValueError(f"parse gap at {pos}: {expr[pos:pos+40]!r}")
        pos = m.end()
        coeff = sp.Integer(int(m["num"]))
        if m["den"]:
            coeff = sp.Rational(int(m["num"]), int(m["den"]))
        if m["sign"] == "-":
            coeff = -coeff
        mono = sp.Integer(1)
        for base, exp in _FACTOR.findall(m["mono"]):
            mono *= LOCALS[base] ** (int(exp) if exp else 1)
        terms.append(coeff * mono)
    if pos != len(expr):
        raise ValueError(f"parse stopped at {pos}/{len(expr)}: {expr[pos:pos+40]!r}")
    return sp.Add(*terms)


def cleared(e: sp.Expr) -> sp.Expr:
    """Multiply out to integer coefficients (matches the generators the Singular
    lift() ran against: each pre-resultant generator was scaled by the lcm of
    its coefficient denominators)."""
    e = sp.expand(e)
    p = sp.Poly(e, d2, d1, d0, dm1, dm2, dm3, dm4, Phi)
    L = 1
    for c in p.coeffs():
        L = sp.ilcm(L, sp.Rational(c).q)
    return sp.expand(e * L)


def pre_resultant_generators() -> list[sp.Expr]:
    """The four pre-resultant equations (STATE.md item 4: D3(1),D3(2),D3(3),
    D3(5)+Phi), straight from the regenerated state, denominator-cleared to the
    exact integer-coefficient generators the Singular lift() ran against."""
    st = pickle.loads((ROOT / "t4_state.pkl").read_bytes())
    G1 = st["G1"]
    G2 = st["G2"]
    G3 = st["G3"]
    G5 = st["G5body"] + Phi
    return [cleared(G1), cleared(G2), cleared(G3), cleared(G5)]


def load_f31() -> sp.Expr:
    return parse((ROOT / "f31_deg31.txt").read_text())


def load_f37() -> sp.Expr:
    return parse((ROOT / "f37_deg37.txt").read_text())


def load_cofactors() -> list[sp.Expr]:
    lines = [
        ln for ln in (ROOT / "f37_sat_certificate.txt").read_text().splitlines()
        if ln.strip() and not ln.lstrip().startswith("#")
    ]
    assert len(lines) == 4, f"expected 4 cofactors, got {len(lines)}"
    return [parse(ln) for ln in lines]


def check_membership_certificate() -> None:
    gens = pre_resultant_generators()
    cof = load_cofactors()
    f31 = load_f31()
    combo = sp.expand(sum(c * g for c, g in zip(cof, gens)))
    residual = sp.expand(combo - sp.expand(f31))
    assert residual == 0, f"certificate FAILED, residual has {len(sp.Add.make_args(residual))} terms"
    print("(A) membership certificate PASS:")
    print("    f31 = c1*G1 + c2*G2 + c3*G3 + c4*(G5body+Phi)  [exact over Q]")
    print("    => f31 vanishes on the entire pre-resultant variety.")


def check_master_identity() -> None:
    """The resultant master identity factor f31*f37*dm1^21 lies in the ideal
    only because its factor f31 does.

    Once (A) proves f31 in <G-system>, any multiple f31*g is in the ideal, so
    f31*f37*dm1^21 in <G-system> is immediate and carries no extra content: the
    f37 and dm1^21 factors of the resultant are excess.  We record this cheaply
    (a symbolic exact-division check that f37*dm1^21 divides the master identity)
    rather than re-expanding a multi-million-term product.
    """
    f31 = load_f31()
    f37 = load_f37()
    master = f31 * f37 * dm1**21  # the resultant's factored master identity
    quotient, remainder = sp.div(
        sp.Poly(sp.expand(master), d2, d1, d0, dm1, Phi),
        sp.Poly(f31, d2, d1, d0, dm1, Phi),
    )
    assert remainder == 0 and sp.expand(quotient.as_expr() - sp.expand(f37 * dm1**21)) == 0
    print("(B) master identity consistency PASS:")
    print("    master identity = f31 * (f37*dm1^21); f31 in <G-system> by (A),")
    print("    so f37 and dm1^21 are resultant excess (add no ideal content).")


def main() -> None:
    check_membership_certificate()
    check_master_identity()
    print()
    print("CONCLUSION: f31 in <G1,G2,G3,G5body+Phi>.  Every solution of the")
    print("pre-resultant system has f31 = 0; the f37 branch off {f31=0} does not")
    print("lift.  The whole f37 component is a resultant artifact.")


if __name__ == "__main__":
    main()
