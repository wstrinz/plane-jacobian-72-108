#!/usr/bin/env python3
'''Automated exact convolution descent for the f31 master identity.

Requested coefficients of ``sum(Phi**f*e**(21-3*f)*h_f for f in range(8))``
are extracted from source-linked ``h_f`` formulas by sparse convolution; the
full degree-242 polynomial is never expanded.  The executable gate proves the
constant-E T5_90_T1 ansatz in ``T5_90_T1.md`` section 3 impossible: degrees
238 through 227 force the listed coefficients and degree 226 is nonzero,
including at the fixed ``c=-1/6630``.

For a general ansatz, ``FORCED`` is only a reduction and ``UNRESOLVED`` is the
first exact equation this driver cannot soundly force; neither is a kill.  The
gate does not settle nonconstant-E T1 or the open T2 components.
'''

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Mapping, Sequence

import sympy as sp
import t5_90t1_verify as base

def _require(_cond, _msg):
    """Proof-critical check: fails loudly and exits nonzero, unaffected by python -O."""
    if not _cond:
        import sys as _sys
        print("FAIL: " + str(_msg))
        _sys.exit(1)


y = base.y
SparsePolynomial = dict[int, sp.Expr]


def from_expr(expr: sp.Expr, variable: sp.Symbol = y) -> SparsePolynomial:
    '''Convert an exact univariate expression to degree: coefficient.'''
    poly = sp.Poly(sp.expand(sp.sympify(expr)), variable)
    return {monomial[0]: coefficient for monomial, coefficient in poly.terms()}


def add(left: SparsePolynomial, right: SparsePolynomial) -> SparsePolynomial:
    '''Add sparse polynomials without expanding coefficient expressions.'''
    out = dict(left)
    for degree, coefficient in right.items():
        out[degree] = out.get(degree, sp.Integer(0)) + coefficient
    return {degree: coefficient for degree, coefficient in out.items()
            if coefficient != 0}


def multiply(left: SparsePolynomial, right: SparsePolynomial) -> SparsePolynomial:
    '''Convolve two sparse polynomials exactly.'''
    out: SparsePolynomial = {}
    for i, left_coefficient in left.items():
        for j, right_coefficient in right.items():
            degree = i + j
            out[degree] = (out.get(degree, sp.Integer(0))
                           + left_coefficient*right_coefficient)
    return out


def power(poly: SparsePolynomial, exponent: int) -> SparsePolynomial:
    '''Raise a sparse polynomial to a nonnegative integer power.'''
    if exponent < 0:
        raise ValueError('exponent must be nonnegative')
    result, factor, remaining = {0: sp.Integer(1)}, poly, exponent
    while remaining:
        if remaining & 1:
            result = multiply(result, factor)
        remaining >>= 1
        if remaining:
            factor = multiply(factor, factor)
    return result


def evaluate_h(
    expr: sp.Expr,
    substitutions: Mapping[sp.Symbol, SparsePolynomial],
    power_cache: dict[tuple[sp.Symbol, int], SparsePolynomial] | None = None,
) -> SparsePolynomial:
    '''Evaluate one source h_f by sparse convolution.'''
    cache = power_cache if power_cache is not None else {}
    variables = (base.d0, base.d1, base.d2, base.e)
    out: SparsePolynomial = {}
    for monomial, coefficient in sp.Poly(expr, *variables).terms():
        term = {0: coefficient}
        for symbol, exponent in zip(variables, monomial):
            key = (symbol, exponent)
            if key not in cache:
                cache[key] = power(substitutions[symbol], exponent)
            term = multiply(term, cache[key])
        out = add(out, term)
    return out


@dataclass(frozen=True)
class Ansatz:
    d2: sp.Expr
    d1: sp.Expr
    d0: sp.Expr
    e: sp.Expr
    sigma: sp.Expr | None
    substitutions: Mapping[sp.Symbol, SparsePolynomial]
    unknowns: tuple[sp.Symbol, ...]
    parameters: frozenset[sp.Symbol]


def build_ansatz(
    *, d2: sp.Expr | None = None, d1: sp.Expr | None = None,
    e: sp.Expr | None = None, d0: sp.Expr | None = None,
    sigma: sp.Expr | None = None, degrees: Mapping[str, int] | None = None,
    prefixes: Mapping[str, str] | None = None, zero: Iterable[str] = (),
    parameters: Iterable[sp.Symbol] = (),
    unknowns: Iterable[sp.Symbol] | None = None,
) -> Ansatz:
    '''Build an exact ansatz from expressions and/or prescribed degrees.

    Omitted polynomials named in degrees receive generic coefficients. Exactly
    one of d0 and sigma is required; the latter sets d0=(d2**2+sigma)/4.
    zero accepts d2, d1, d0, sigma, and e, notably all T2 zero flags. Symbols
    in parameters are constants rather than forceable unknowns.
    '''
    degree_map = dict(degrees or {})
    prefix_map = {'d2': 'a', 'd1': 'b', 'd0': 'r', 'sigma': 's', 'e': 'g'}
    prefix_map.update(prefixes or {})
    zero_names = set(zero)
    if zero_names - set(prefix_map):
        raise ValueError(f'unknown zero flags: {sorted(zero_names-set(prefix_map))}')
    supplied = {'d2': d2, 'd1': d1, 'e': e, 'd0': d0, 'sigma': sigma}
    generated: list[sp.Symbol] = []

    def resolve(name: str) -> sp.Expr | None:
        expression = supplied[name]
        if name in zero_names:
            if expression is not None and sp.sympify(expression) != 0:
                raise ValueError(f'{name} is supplied nonzero and zero-flagged')
            if name in degree_map:
                raise ValueError(f'{name} has a degree and a zero flag')
            return sp.Integer(0)
        if expression is not None:
            if name in degree_map:
                raise ValueError(f'{name} has an expression and a degree')
            return sp.sympify(expression)
        if name not in degree_map:
            return None
        degree = degree_map[name]
        if degree < 0:
            raise ValueError(f'degree of {name} must be nonnegative')
        coefficients = tuple(sp.symbols(f'{prefix_map[name]}0:{degree+1}'))
        generated.extend(coefficients)
        return sum(coefficient*y**index
                   for index, coefficient in enumerate(coefficients))

    d2_expr, d1_expr, e_expr = resolve('d2'), resolve('d1'), resolve('e')
    missing = [name for name, value in
               (('d2', d2_expr), ('d1', d1_expr), ('e', e_expr))
               if value is None]
    if missing:
        raise ValueError(f'missing expressions or degrees for {missing}')
    d0_expr, sigma_expr = resolve('d0'), resolve('sigma')
    if (d0_expr is None) == (sigma_expr is None):
        raise ValueError('specify exactly one of d0 and sigma')
    _require(d2_expr is not None and d1_expr is not None and e_expr is not None, "d2_expr is not None and d1_expr is not None and e_expr is not None")
    if sigma_expr is not None:
        d0_expr = sp.expand((d2_expr**2 + sigma_expr)/4)
    _require(d0_expr is not None, "d0_expr is not None")

    parameter_set = frozenset(parameters)
    if unknowns is None:
        inferred = set(generated)
        for expression in (d2_expr, d1_expr, d0_expr, e_expr):
            inferred.update(expression.free_symbols)
        inferred.difference_update(parameter_set | {y})
        unknown_tuple = tuple(sorted(inferred, key=sp.default_sort_key))
    else:
        unknown_tuple = tuple(unknowns)
        if set(unknown_tuple) & parameter_set:
            raise ValueError('a symbol cannot be both unknown and parameter')
    substitutions = {
        base.d0: from_expr(d0_expr), base.d1: from_expr(d1_expr),
        base.d2: from_expr(d2_expr), base.e: from_expr(e_expr),
    }
    return Ansatz(d2_expr, d1_expr, d0_expr, e_expr, sigma_expr,
                  substitutions, unknown_tuple, parameter_set)


@dataclass(frozen=True)
class DescentStep:
    degree: int
    verdict: str
    factored: sp.Expr
    substitution: tuple[sp.Symbol, sp.Expr] | None = None
    equation: sp.Equality | None = None


@dataclass(frozen=True)
class DescentResult:
    verdict: str
    steps: tuple[DescentStep, ...]
    substitutions: Mapping[sp.Symbol, sp.Expr] = field(default_factory=dict)
    stopping_degree: int | None = None


class ConvolutionDescent:
    '''Cached exact coefficient engine for all eight master terms.'''
    def __init__(self, ansatz: Ansatz, *, c: sp.Expr | None = None,
                 phi: sp.Expr | None = None,
                 h: Mapping[int, sp.Expr] | None = None) -> None:
        if phi is None:
            if c is None:
                raise ValueError('supply c or phi')
            phi = sp.sympify(c)*(y+1)**30*base.q
        elif c is not None:
            raise ValueError('supply c or phi, not both')
        self.ansatz, self.h = ansatz, dict(base.load_h() if h is None else h)
        if sorted(self.h) != list(range(8)):
            raise ValueError(f'expected h_0,...,h_7; found {sorted(self.h)}')
        self.phi, self.e_poly = from_expr(phi), ansatz.substitutions[base.e]
        self._source_powers: dict[tuple[sp.Symbol, int], SparsePolynomial] = {}
        self._h: dict[int, SparsePolynomial] = {}
        self._phi_powers: dict[int, SparsePolynomial] = {}
        self._e_powers: dict[int, SparsePolynomial] = {}
        self._coefficients: dict[int, sp.Expr] = {}

    @staticmethod
    def _cached_power(cache: dict[int, SparsePolynomial],
                      poly: SparsePolynomial, exponent: int) -> SparsePolynomial:
        if exponent not in cache:
            cache[exponent] = power(poly, exponent)
        return cache[exponent]

    def term_coefficient(self, f: int, target: int) -> sp.Expr:
        '''Exact y**target coefficient of master term f.'''
        if f not in range(8):
            raise ValueError('f must be in 0,...,7')
        if f not in self._h:
            self._h[f] = evaluate_h(self.h[f], self.ansatz.substitutions,
                                     self._source_powers)
        phi_power = self._cached_power(self._phi_powers, self.phi, f)
        e_power = self._cached_power(self._e_powers, self.e_poly, 21-3*f)
        total = sp.Integer(0)
        for i, phi_coefficient in phi_power.items():
            for j, e_coefficient in e_power.items():
                h_degree = target-i-j
                if h_degree in self._h[f]:
                    total += (phi_coefficient*e_coefficient
                              * self._h[f][h_degree])
        return total

    def master_coefficient(self, target: int) -> sp.Expr:
        '''Exact coefficient, always summing every f=0,...,7.'''
        if target not in self._coefficients:
            # Lower f terms are never truncated; f=4 first enters at T1 y^239.
            self._coefficients[target] = sum(
                (self.term_coefficient(f, target) for f in range(8)),
                sp.Integer(0))
        return self._coefficients[target]

    @staticmethod
    def _forced_square(residual: sp.Expr, active: Sequence[sp.Symbol]
                       ) -> tuple[sp.Symbol, sp.Expr] | None:
        if len(active) != 1:
            return None
        unknown = active[0]
        try:
            polynomial = sp.Poly(residual, unknown)
        except sp.PolynomialError:
            return None
        if polynomial.degree() != 2 or polynomial.nth(2).is_zero is not False:
            return None
        a, b, constant = polynomial.nth(2), polynomial.nth(1), polynomial.nth(0)
        if sp.factor(b**2-4*a*constant) != 0:
            return None
        value = sp.factor(-b/(2*a))
        return ((unknown, value)
                if sp.factor(residual.subs(unknown, value)) == 0 else None)

    def descend(self, start_degree: int, floor: int, *,
                initial_substitutions: Mapping[sp.Symbol, sp.Expr] | None = None
                ) -> DescentResult:
        '''Walk downward, forcing only unique exact square equations.'''
        if start_degree < floor:
            raise ValueError('start_degree must be at least floor')
        substitutions: dict[sp.Symbol, sp.Expr] = dict(initial_substitutions or {})
        steps: list[DescentStep] = []
        if not set(substitutions).issubset(self.ansatz.unknowns):
            raise ValueError('initial substitutions include a non-unknown')
        for target in range(start_degree, floor-1, -1):
            residual = sp.factor(self.master_coefficient(target).subs(substitutions))
            if residual == 0:
                steps.append(DescentStep(target, 'IDENTITY', sp.Integer(0)))
                continue
            active = tuple(u for u in self.ansatz.unknowns
                           if u not in substitutions and residual.has(u))
            forced = self._forced_square(residual, active)
            if forced:
                substitutions[forced[0]] = forced[1]
                steps.append(DescentStep(target, 'FORCED', residual, forced))
                continue
            equation = sp.Eq(residual, 0, evaluate=False)
            if not active and (residual.is_zero is False
                               or sp.ask(sp.Q.nonzero(residual)) is True):
                steps.append(DescentStep(target, 'CONTRADICTION', residual,
                                         equation=equation))
                return DescentResult('CONTRADICTION', tuple(steps),
                                     dict(substitutions), target)
            steps.append(DescentStep(target, 'UNRESOLVED', residual,
                                     equation=equation))
            return DescentResult('UNRESOLVED', tuple(steps),
                                 dict(substitutions), target)
        return DescentResult('FORCED', tuple(steps), dict(substitutions), floor)


def master_coefficient(ansatz: Ansatz, target: int, *, c: sp.Expr | None = None,
                       phi: sp.Expr | None = None,
                       h: Mapping[int, sp.Expr] | None = None) -> sp.Expr:
    '''Functional wrapper for one all-level coefficient.'''
    return ConvolutionDescent(ansatz, c=c, phi=phi, h=h).master_coefficient(target)


def descend(ansatz: Ansatz, start_degree: int, floor: int, *,
            c: sp.Expr | None = None, phi: sp.Expr | None = None,
            h: Mapping[int, sp.Expr] | None = None,
            initial_substitutions: Mapping[sp.Symbol, sp.Expr] | None = None
            ) -> DescentResult:
    '''Functional wrapper for a complete descent.'''
    return ConvolutionDescent(ansatz, c=c, phi=phi, h=h).descend(
        start_degree, floor, initial_substitutions=initial_substitutions)


def _print_result(result: DescentResult) -> None:
    print('T5_90_T1 constant-E automated convolution descent')
    for step in result.steps:
        print(f'  degree {step.degree}: {step.verdict}')
        if step.verdict != 'IDENTITY':
            print(f'    coefficient = {step.factored}')
        if step.substitution:
            print(f'    substitute {step.substitution[0]} = {step.substitution[1]}')
        if step.verdict == 'UNRESOLVED':
            print(f'    residual equation: {step.equation}')


def _self_test() -> None:
    c, gamma, v0 = sp.symbols('c gamma v0', nonzero=True)
    a0, a1, a2, a3, a4 = sp.symbols('a0:5')
    s0, s1, s2, s3 = sp.symbols('s0:4')
    ansatz = build_ansatz(
        d2=a4*y**4+a3*y**3+a2*y**2+a1*y+a0,
        sigma=s3*y**3+s2*y**2+s1*y+s0,
        d1=-gamma**4/(4096*c)*(y**2+sp.Rational(25, 4)*y+v0),
        e=gamma*(y+1)**9, parameters=(c, gamma))
    result = ConvolutionDescent(ansatz, c=c).descend(238, 226)
    _print_result(result)
    expected = {
        v0: sp.Rational(525, 32), s3: 0,
        a4: -sp.Integer(95200)*c/gamma**3, s2: 0,
        a3: sp.Integer(255850)*c/gamma**3, s1: 0,
        a2: -sp.Integer(513451)*c/gamma**3, s0: 0,
        a1: -sp.Rational(10656467, 8)*c/gamma**3,
        a0: sp.Rational(132899897, 8)*c/gamma**3,
    }
    expected_kinds = {
        238: 'FORCED', 237: 'FORCED', 236: 'FORCED', 235: 'FORCED',
        234: 'FORCED', 233: 'FORCED', 232: 'FORCED', 231: 'FORCED',
        230: 'FORCED', 229: 'IDENTITY', 228: 'FORCED', 227: 'IDENTITY',
        226: 'CONTRADICTION',
    }
    final = sp.Integer(29570349989420274657771126784)*c**5*gamma**8
    failures = []
    if (result.verdict, result.stopping_degree) != ('CONTRADICTION', 226):
        failures.append(f'verdict/degree={result.verdict}/{result.stopping_degree}')
    if {step.degree: step.verdict for step in result.steps} != expected_kinds:
        failures.append('step verdicts differ')
    if result.substitutions != expected:
        failures.append(f'substitutions={result.substitutions}')
    if not result.steps or sp.factor(result.steps[-1].factored-final) != 0:
        actual = result.steps[-1].factored if result.steps else None
        failures.append(f'degree-226={actual}')
    if failures:
        print('FAIL: T5_90_T1 constant-E convolution descent gate')
        for failure in failures:
            print(f'  {failure}')
        raise AssertionError('; '.join(failures))
    _require(final.subs(c, -sp.Rational(1, 6630)).is_zero is False, "final.subs(c, -sp.Rational(1, 6630)).is_zero is False")
    print('PASS: T5_90_T1 constant-E CONTRADICTION at degree 226')


if __name__ == '__main__':
    _self_test()
