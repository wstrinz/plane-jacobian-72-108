#!/usr/bin/env python3
"""SERIALIZER ROUND-TRIP GUARD -- regression check for SERIALIZER_BUG.md.

Feeds rational-coefficient polynomials through BOTH production serializers

    blowup_diagnosis.sing_poly_intcoeff   (char 0 -> integer-coefficient string)
    modular_triage.poly_to_singular_modp  (char p -> F_p-coefficient string)

and verifies the emitted STRING re-parses to exactly the intended polynomial.
Nothing here trusts the serializer's own arithmetic: the expected value is
recomputed independently from sp.Poly coefficient dictionaries.

Three parts, all of which must pass for exit 0:

  A  live serializers round-trip on the corpus.
  B  a verbatim re-implementation of the PRE-FIX serializers is fed through the
     SAME checker and must be REPORTED AS FAILING.  This is what proves the
     guard has teeth -- if part B ever passes, the checker has stopped
     detecting the original defect and is worthless.
  C  with sp.together monkeypatched to the identity (which is precisely how
     sp.cancel behaved under sympy 1.14.0 + python-flint 0.9.0, and the reason
     the denominator-clearing step became a no-op), the live serializers must
     RAISE -- never silently truncate.  This checks the fail-loud guard is
     actually wired in.

Usage:  python serializer_roundtrip_verify.py [--quiet]
Exit 0 = all three parts as expected.  Exit 1 = regression.

NEW file.  READ-ONLY on every existing module and artifact.
"""
from __future__ import annotations

import argparse
import math
import sys
from fractions import Fraction

import sympy as sp

import blowup_diagnosis as bd
import modular_triage as mt

P = 10007          # the triage prime used throughout the repo
BIGP = 100019


# ==========================================================================
#  corpus
# ==========================================================================
def corpus() -> list[tuple[str, sp.Expr, list[sp.Symbol]]]:
    c0_0, c0_1, E, X, r, w = sp.symbols('c0_0 c0_1 E X r w')
    out = [
        # the class quartic from J6 -- the polynomial that exposed the bug.
        # Pre-fix this serialized to the bare monomial c0_0^4.
        ('J6 class quartic',
         c0_0**4 + c0_0**3/4 + sp.Rational(5, 32)*c0_0**2
         + sp.Rational(15, 128)*c0_0 + sp.Rational(195, 2048), [c0_0]),
        # the second documented J6 case: huge, unequal denominators
        ('J6 E-relation',
         sp.Rational(3981312, 221)*E**21
         - sp.Rational(2305843009213693952, 400329564123571875)*E**8, [E]),
        # multivariate with mixed denominators
        ('J6 class relation (2 vars)',
         c0_0*c0_1 + c0_0/8 - c0_1**3/2 - c0_1**2/8
         - sp.Rational(5, 64)*c0_1 - sp.Rational(15, 256), [c0_0, c0_1]),
        # negative rationals: int() truncates TOWARD ZERO, so sign matters
        ('negative rationals',
         -sp.Rational(7, 3)*E**2 + sp.Rational(1, 3)*E - sp.Rational(2, 5),
         [E]),
        # coefficients with |c| < 1 : pre-fix these vanished entirely
        ('sub-unit coefficients',
         E**3 + E**2/7 + E/11 + sp.Rational(1, 13), [E]),
        # a genuine denominator-1 case: the bug is INERT here, must still pass
        ('integral (bug inert)',
         2048*c0_0**4 + 512*c0_0**3 + 320*c0_0**2 + 240*c0_0 + 195, [c0_0]),
        # Rabinowitsch saturation generator: integral, appears in every system
        ('saturation generator', E*X*w - 1, [E, X, w]),
        # zero
        ('zero', sp.Integer(0), [E]),
        # number-field variable kept as a ring variable
        ('marked-root polynomial',
         r**4/2048 - r**3/4 + sp.Rational(5, 32)*r**2 - sp.Rational(15, 128)*r,
         [r]),
        # denominator sharing a factor with nothing; large prime-ish denom
        ('coprime denominators',
         X**2/9 + X/49 + sp.Rational(1, 121), [X]),
    ]
    return out


# ==========================================================================
#  independent expectations (never call the serializers)
# ==========================================================================
def coeffs_of(expr, gens) -> dict[tuple, Fraction]:
    expr = sp.expand(sp.sympify(expr))
    if expr == 0:
        return {}
    p = sp.Poly(expr, *gens, domain=sp.QQ)
    return {m: Fraction(int(c.p), int(c.q)) for m, c in p.terms() if c != 0}


def common_den(cs: dict[tuple, Fraction]) -> int:
    d = 1
    for c in cs.values():
        d = d * c.denominator // math.gcd(d, c.denominator)
    return d


def parse_singular(s: str, gens) -> sp.Expr:
    """Re-parse an emitted Singular/msolve polynomial string."""
    loc = {g.name: g for g in gens}
    return sp.expand(sp.sympify(s.replace('^', '**'), locals=loc))


# ==========================================================================
#  the two round-trip contracts
# ==========================================================================
def check_char0(fn, expr, gens) -> str | None:
    """fn(expr,gens) must emit an INTEGER-coefficient string equal to D*expr
    for a single positive rational-clearing constant D.  Returns None on
    success, else a failure description."""
    cs = coeffs_of(expr, gens)
    s = fn(expr, gens)
    got = parse_singular(s, gens)
    if not cs:
        return None if got == 0 else f'zero mis-serialized as {s!r}'
    # every emitted coefficient must be an integer
    gc = coeffs_of(got, gens)
    if any(c.denominator != 1 for c in gc.values()):
        return f'emitted non-integer coefficients: {s[:70]!r}'
    if not gc:
        return f'nonzero polynomial serialized to zero: {s[:70]!r}'
    # got must be a nonzero RATIONAL multiple of expr, with the same support
    if set(gc) != set(cs):
        lost = sorted(set(cs) - set(gc))
        return (f'support changed: {len(cs)} terms in, {len(gc)} out, '
                f'{len(lost)} monomial(s) dropped -- {s[:60]!r}')
    m0 = next(iter(cs))
    ratio = gc[m0] / cs[m0]
    if ratio == 0:
        return f'serialized to zero multiple: {s[:70]!r}'
    for m in cs:
        if gc[m] != ratio * cs[m]:
            return (f'coefficient of {m} is {gc[m]}, expected {ratio*cs[m]} '
                    f'(ratio {ratio}) -- not a scalar multiple of the input')
    if ratio.denominator != 1:
        return f'clearing constant {ratio} is not an integer'
    return None


def check_modp(fn, expr, gens, p) -> str | None:
    """fn(expr,gens,p) must emit a string whose F_p coefficient vector is a
    nonzero scalar multiple of the true reduction of expr mod p."""
    cs = coeffs_of(expr, gens)
    s = fn(expr, gens, p)
    got = parse_singular(s, gens)
    if not cs:
        return None if got == 0 else f'zero mis-serialized as {s!r}'
    # true reduction: num * inverse(den) mod p, dropping terms that vanish
    want = {}
    for m, c in cs.items():
        if c.denominator % p == 0:
            return None          # reduction undefined; not this check's job
        v = (c.numerator % p) * pow(c.denominator % p, -1, p) % p
        if v:
            want[m] = v
    gc = {m: int(c) % p for m, c in coeffs_of(got, gens).items()}
    gc = {m: v for m, v in gc.items() if v}
    if not want:
        return None if not gc else f'expected 0 mod {p}, got {s[:60]!r}'
    if set(gc) != set(want):
        return (f'mod-{p} support changed: {len(want)} terms expected, '
                f'{len(gc)} emitted -- {s[:60]!r}')
    m0 = next(iter(want))
    scale = gc[m0] * pow(want[m0], -1, p) % p
    if scale == 0:
        return f'serialized to zero mod {p}: {s[:60]!r}'
    for m in want:
        if (gc[m] - scale * want[m]) % p:
            return (f'mod-{p} coefficient of {m} is {gc[m]}, expected '
                    f'{scale*want[m]%p} -- not a scalar multiple')
    return None


# ==========================================================================
#  PART B: verbatim pre-fix re-implementations (must be caught)
# ==========================================================================
def prefix_intcoeff(expr, gens):
    """blowup_diagnosis.sing_poly_intcoeff AS IT WAS before the fix."""
    expr = sp.sympify(expr)
    if expr == 0:
        return '0'
    num, den = sp.fraction(sp.cancel(expr))       # <-- no-op in this env
    poly = sp.Poly(sp.expand(num), *gens)
    terms = []
    for monom, coeff in poly.terms():
        c = int(coeff)                            # <-- silent truncation
        if c == 0:
            continue
        factors = [str(c)] if (c != 1 or all(e == 0 for e in monom)) else []
        for g, e in zip(gens, monom):
            if e == 1:
                factors.append(g.name)
            elif e > 1:
                factors.append(f'{g.name}^{e}')
        terms.append('*'.join(factors) if factors else str(c))
    return ('+'.join(terms)).replace('+-', '-') or '0'


def prefix_modp(expr, gens, p):
    """modular_triage.poly_to_singular_modp AS IT WAS before the fix."""
    expr = sp.sympify(expr)
    if expr == 0:
        return '0'
    num, den = sp.fraction(sp.cancel(expr))       # <-- no-op in this env
    poly = sp.Poly(sp.expand(num), *gens)
    den_poly = sp.Poly(sp.expand(den), *gens)
    den_val = int(den_poly.coeff_monomial(1))
    den_inv = pow(den_val % p, -1, p)
    terms = []
    for monom, coeff in poly.terms():
        c = int(coeff)                            # <-- silent truncation
        if c % p == 0:
            continue
        cm = (c * den_inv) % p
        factors = [str(cm)] if (cm != 1 or all(e == 0 for e in monom)) else []
        for g, e in zip(gens, monom):
            if e == 1:
                factors.append(g.name)
            elif e > 1:
                factors.append(f'{g.name}^{e}')
        terms.append('*'.join(factors) if factors else str(cm))
    return ('+'.join(terms)).replace('+-', '-') or '0'


# ==========================================================================
#  driver
# ==========================================================================
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--quiet', action='store_true')
    args = ap.parse_args()
    say = (lambda *a: None) if args.quiet else (lambda *a: print(*a))

    cases = corpus()
    failures: list[str] = []

    # ---------------- PART A : live serializers must round-trip ------------
    say('== A. live serializers ==')
    for name, expr, gens in cases:
        for label, err in (
                ('char0    ', check_char0(bd.sing_poly_intcoeff, expr, gens)),
                (f'mod {P}', check_modp(mt.poly_to_singular_modp, expr, gens, P)),
                (f'mod {BIGP}', check_modp(mt.poly_to_singular_modp, expr, gens, BIGP)),
        ):
            if err:
                failures.append(f'A/{name}/{label.strip()}: {err}')
                say(f'  FAIL {name:28s} {label.strip():10s} {err}')
            else:
                say(f'  ok   {name:28s} {label.strip()}')

    # ---------------- PART B : the guard must catch the pre-fix code -------
    say('== B. pre-fix re-implementation must be DETECTED ==')
    caught_q = caught_p = 0
    for name, expr, gens in cases:
        cs = coeffs_of(expr, gens)
        if common_den(cs) == 1:
            continue                       # bug is provably inert here
        try:
            e_q = check_char0(prefix_intcoeff, expr, gens)
        except Exception as ex:
            e_q = f'raised {type(ex).__name__}'
        try:
            e_p = check_modp(prefix_modp, expr, gens, P)
        except Exception as ex:
            e_p = f'raised {type(ex).__name__}'
        caught_q += bool(e_q)
        caught_p += bool(e_p)
        say(f'  {name:28s} char0={"CAUGHT" if e_q else "MISSED"} '
            f'modp={"CAUGHT" if e_p else "MISSED"}')
        if not e_q:
            failures.append(f'B/{name}/char0: guard did NOT catch pre-fix output')
        if not e_p:
            failures.append(f'B/{name}/modp: guard did NOT catch pre-fix output')
    if caught_q == 0 or caught_p == 0:
        failures.append('B: no rational corpus case exercised the pre-fix path')
    say(f'  caught {caught_q} char-0 and {caught_p} mod-p pre-fix corruptions')

    # ---------------- PART C : fail-loud guard must be live ----------------
    say('== C. with sp.together neutralised the live code must RAISE ==')
    name, expr, gens = cases[0]
    for mod, fn, call in ((bd, 'sing_poly_intcoeff',
                           lambda: bd.sing_poly_intcoeff(expr, gens)),
                          (mt, 'poly_to_singular_modp',
                           lambda: mt.poly_to_singular_modp(expr, gens, P))):
        orig = sp.together
        sp.together = lambda e, *a, **k: e        # reproduce the 1.14 no-op
        try:
            out = call()
            failures.append(f'C/{mod.__name__}.{fn}: silently returned {out[:50]!r} '
                            f'instead of raising')
            say(f'  FAIL {mod.__name__}.{fn} returned {out[:40]!r}')
        except Exception as ex:
            say(f'  ok   {mod.__name__}.{fn} raised '
                f'{type(ex).__name__}: {str(ex)[:60]}')
        finally:
            sp.together = orig

    # ---------------- verdict ---------------------------------------------
    if failures:
        print(f'SERIALIZER ROUND-TRIP: FAIL ({len(failures)} problem(s))')
        for f in failures:
            print(f'  - {f}')
        return 1
    if not args.quiet:
        print(f'SERIALIZER ROUND-TRIP: PASS '
              f'({len(cases)} polynomials x 3 encodings, pre-fix behaviour '
              f'detected, fail-loud guard live)')
    return 0


if __name__ == '__main__':
    sys.exit(main())
