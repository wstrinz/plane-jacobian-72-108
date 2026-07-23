#!/usr/bin/env python3
'''Exact q-support elimination for the a9 b1000 T2 priority cell.

The seven states are those denoted (R9) in ``T5_T2_COLUMN.md``:

    (f,z,g;D,Sigma) = (0,z,2*z;10,2+z),  z = 0,...,6.

For one chosen root r of the fixed irreducible quartic q, this module imposes

    e     = gamma*(y+1)**9*(y-r),
    sigma = (y-r)**2*G(y),       deg(G) = z,
    d1    = 0,
    deg(d2) <= 4.

The quotient Q[r]/(q(r)) is represented exactly as a polynomial quotient:
``r`` remains a Groebner-ring variable, ``q(r)`` is a generator of every
ideal, and master coefficients are replaced by their remainders modulo q.
The latter replacement is exact because q is retained in the ideal.

This is deliberately a client of :mod:`convolution_elim`. In particular it
uses that landed module's high-coefficient engine, normalization convention,
unit-ideal test, and killable Groebner subprocess. A small accumulation loop
is needed because the landed public ``eliminate`` API cannot seed q(r)=0 and
only accepts symbols (not G(r)) as Rabinowitsch nonzero constraints. Unlike
the landed public API's single total-time budget, each basis update here gets
the requested hard 120-second subprocess budget.

CONTRADICTION is printed only after the returned basis has been checked to be
exactly {1}. Such a result is still labelled a candidate new kill pending an
independent audit. This module intentionally does not test any chain-level
F**3*G = 3072*c**6*Z**2 relation.
'''

from __future__ import annotations

from dataclasses import dataclass
import time
from typing import Iterable, Mapping, Sequence

import sympy as sp

import convolution_descent as cd
import convolution_elim as landed


y = cd.y
r = sp.Symbol('r')
gamma = sp.Symbol('gamma')
Q_R = 2048*r**4 - 512*r**3 + 320*r**2 - 240*r + 195

START_DEGREE = 251
COEFFICIENT_LIMIT = 12
GROEBNER_TIMEOUT = 120.0
ORDER = 'grevlex'
PLAIN_R9_Z0_BASIS_SIZE = 34


@dataclass(frozen=True)
class QSupportState:
    z: int
    ansatz: cd.Ansatz
    g_coefficients: tuple[sp.Symbol, ...]
    G: sp.Expr
    saturation_factors: tuple[sp.Expr, ...]
    saturation_symbols: tuple[sp.Symbol, ...]
    saturation_equations: tuple[sp.Expr, ...]


@dataclass(frozen=True)
class QSupportStep:
    degree: int
    action: str
    basis_size: int
    substitution: tuple[sp.Symbol, sp.Expr] | None = None


@dataclass(frozen=True)
class QSupportResult:
    z: int
    verdict: str
    steps: tuple[QSupportStep, ...]
    substitutions: Mapping[sp.Symbol, sp.Expr]
    groebner_basis: tuple[sp.Expr, ...]
    residual_equations: tuple[sp.Expr, ...]
    generators: tuple[sp.Symbol, ...]
    consumed: int
    requested: int
    start_degree: int
    first_unconsumed_degree: int | None
    reason: str
    elapsed_seconds: float
    timed_out: bool
    unit_certificate_checked: bool


@dataclass(frozen=True)
class ValuationCheck:
    z: int
    e_value: sp.Expr
    e_first_unit_part: sp.Expr
    sigma_value: sp.Expr
    sigma_first_value: sp.Expr
    sigma_second_unit_part: sp.Expr
    r_plus_one_is_unit: bool
    gamma_saturated: bool
    G_at_r_saturated: bool
    passed: bool


@dataclass(frozen=True)
class RefinementCheck:
    plain_basis_size: int
    qsupport_basis_size: int
    qsupport_verdict: str
    coefficient_specialization_matches: bool
    qsupport_generators_verified: bool
    support_relation_is_new: bool
    genuine_refinement: bool
    timed_out: bool
    elapsed_seconds: float


def quotient_reduce(expression: sp.Expr) -> sp.Expr:
    '''Return the canonical degree-<4 representative modulo q(r).'''
    expression = sp.cancel(sp.sympify(expression))
    if expression == 0 or not expression.has(r):
        return expression
    numerator, denominator = sp.fraction(expression)
    if denominator.has(r):
        raise ValueError('unexpected r-dependent denominator')
    remainder = sp.rem(sp.Poly(sp.expand(numerator), r), sp.Poly(Q_R, r))
    return sp.cancel(remainder.as_expr()/denominator)


def _unique_expressions(expressions: Iterable[sp.Expr]) -> tuple[sp.Expr, ...]:
    unique: list[sp.Expr] = []
    for expression in expressions:
        expression = quotient_reduce(expression)
        if not any(sp.expand(expression-old) == 0 for old in unique):
            unique.append(expression)
    return tuple(unique)


def build_qsupport_ansatz(z: int) -> QSupportState:
    '''Build one R9 state and all exact-order/nonzero constraints.'''
    if z not in range(7):
        raise ValueError('z must be in 0,...,6')
    g_coefficients = tuple(sp.symbols(f'g0:{z+1}'))
    G = sum(coefficient*y**degree
            for degree, coefficient in enumerate(g_coefficients))
    sigma = sp.expand((y-r)**2*G)
    e = gamma*(y+1)**9*(y-r)
    ansatz = cd.build_ansatz(
        degrees={'d2': 4}, d1=sp.Integer(0), e=e, sigma=sigma,
        prefixes={'d2': 'a'})

    # G(r) is essential: a nonzero leading coefficient does not by itself
    # exclude q | G for z >= 4, hence would only prove v_r(sigma) >= 2.
    factors = _unique_expressions(
        (gamma, g_coefficients[-1], G.subs(y, r)))
    rab_symbols = tuple(sp.Symbol(f'_qs_rab_{index}')
                        for index in range(len(factors)))
    rab_equations = tuple(
        quotient_reduce(unit*factor-1)
        for unit, factor in zip(rab_symbols, factors))
    return QSupportState(z, ansatz, g_coefficients, G, factors,
                         rab_symbols, rab_equations)


def _fixed_generators(state: QSupportState,
                      substitutions: Mapping[sp.Symbol, sp.Expr]
                      ) -> tuple[sp.Symbol, ...]:
    active_unknowns = set(state.ansatz.unknowns)-set(substitutions)
    active_unknowns.update(state.saturation_symbols)
    return tuple(sorted(active_unknowns, key=sp.default_sort_key))


def _current_equations(
    state: QSupportState,
    coefficient_equations: Sequence[sp.Expr],
    substitutions: Mapping[sp.Symbol, sp.Expr],
) -> tuple[sp.Expr, ...]:
    equations: list[sp.Expr] = [sp.cancel(Q_R.subs(substitutions))]
    for equation in (*state.saturation_equations, *coefficient_equations):
        reduced = quotient_reduce(sp.sympify(equation).subs(substitutions))
        if reduced != 0:
            equations.append(reduced)
    return tuple(equations)


def _timed_basis(
    equations: Sequence[sp.Expr], generators: Sequence[sp.Symbol],
    *, timeout: float = GROEBNER_TIMEOUT,
) -> sp.GroebnerBasis:
    '''Normalize and invoke the landed killable Groebner subprocess.'''
    generators = tuple(generators)
    if not generators:
        dummy = sp.Symbol('_qs_unit_dummy')
        return sp.groebner([1 if equations else 0], dummy,
                           order=ORDER, domain=sp.QQ)
    normalized = tuple(
        landed._normalize(equation, generators, sp.QQ)
        for equation in equations if equation != 0)
    return landed._timed_groebner(
        normalized, generators, (), ORDER, timeout)


def run_qsupport_state(
    z: int, *, start_degree: int = START_DEGREE,
    target_count: int = COEFFICIENT_LIMIT,
    groebner_timeout: float = GROEBNER_TIMEOUT,
) -> QSupportResult:
    '''Consume at most target_count master coefficients for one q-state.'''
    if target_count <= 0:
        raise ValueError('target_count must be positive')
    if groebner_timeout <= 0:
        raise ValueError('groebner_timeout must be positive')
    started = time.monotonic()
    state = build_qsupport_ansatz(z)
    engine = landed.HighCoefficientEngine(
        state.ansatz, start_degree=start_degree,
        target_count=target_count, c=landed.DEFAULT_C)
    substitutions: dict[sp.Symbol, sp.Expr] = {}
    coefficient_equations: list[sp.Expr] = []
    steps: list[QSupportStep] = []
    generators = _fixed_generators(state, substitutions)
    equations = _current_equations(state, coefficient_equations,
                                   substitutions)

    def finish(verdict: str, reason: str, *, timed_out: bool = False,
               certificate_checked: bool = False) -> QSupportResult:
        basis_expressions = tuple(poly.as_expr() for poly in basis.polys)
        consumed = len(steps)
        first_unconsumed = (None if verdict == 'CONTRADICTION'
                            else start_degree-consumed)
        return QSupportResult(
            z, verdict, tuple(steps), dict(substitutions), basis_expressions,
            tuple(equations), tuple(generators), consumed, target_count,
            start_degree, first_unconsumed, reason,
            time.monotonic()-started, timed_out, certificate_checked)

    try:
        basis = _timed_basis(equations, generators,
                             timeout=groebner_timeout)
    except landed._BudgetExpired:
        # A placeholder proper basis is used only to serialize the cutoff;
        # no CONTRADICTION can be emitted on this path.
        basis = sp.groebner([0], *generators, order=ORDER, domain=sp.QQ)
        return finish('REDUCED', 'TIMEOUT in initial q/saturation basis',
                      timed_out=True)
    if landed._is_one(basis):
        checked = (len(basis.polys) == 1
                   and basis.polys[0].as_expr() == 1)
        if not checked:
            raise AssertionError('unit-ideal predicate without basis {1}')
        return finish('CONTRADICTION',
                      'unit ideal in initial q/saturation constraints',
                      certificate_checked=True)

    for degree in range(start_degree, start_degree-target_count, -1):
        coefficient = quotient_reduce(
            engine.master_coefficient(degree).subs(substitutions))
        reduced = sp.cancel(basis.reduce(coefficient)[1])
        if reduced == 0:
            steps.append(QSupportStep(
                degree, 'IDENTITY_MOD_IDEAL', len(basis.polys)))
            continue

        # Preserve r as the explicit quotient-ring generator. Other unique
        # perfect-square equations use the same landed forcing test.
        active = tuple(
            symbol for symbol in state.ansatz.unknowns
            if symbol != r and symbol not in substitutions
            and reduced.has(symbol))
        forced = cd.ConvolutionDescent._forced_square(reduced, active)
        old_substitutions = dict(substitutions)
        old_equation_count = len(coefficient_equations)
        old_basis, old_generators, old_equations = (
            basis, generators, equations)
        if forced is None:
            coefficient_equations.append(coefficient)
            action = 'ADDED_TO_IDEAL'
        else:
            substitutions[forced[0]] = quotient_reduce(forced[1])
            action = 'FORCED'

        generators = _fixed_generators(state, substitutions)
        equations = _current_equations(state, coefficient_equations,
                                       substitutions)
        try:
            basis = _timed_basis(equations, generators,
                                 timeout=groebner_timeout)
        except landed._BudgetExpired:
            substitutions.clear()
            substitutions.update(old_substitutions)
            del coefficient_equations[old_equation_count:]
            basis, generators, equations = (
                old_basis, old_generators, old_equations)
            return finish(
                'REDUCED',
                f'TIMEOUT while testing unconsumed degree {degree}',
                timed_out=True)

        steps.append(QSupportStep(
            degree, action, len(basis.polys), forced))
        if landed._is_one(basis):
            # Verification gate: a timeout/cutoff can never reach this arm,
            # and the actual returned certificate must literally be {1}.
            checked = (len(basis.polys) == 1
                       and basis.polys[0].as_expr() == 1)
            if not checked:
                raise AssertionError('unit-ideal predicate without basis {1}')
            return finish(
                'CONTRADICTION',
                f'1 in saturated q-support ideal after degree {degree}',
                certificate_checked=True)

    return finish('REDUCED', 'COEFFICIENT_LIMIT')


def check_valuations(z: int) -> ValuationCheck:
    '''Verify the exact r-adic orders by quotient-ring derivatives.'''
    state = build_qsupport_ansatz(z)
    e = state.ansatz.e
    sigma = sp.sympify(state.ansatz.sigma)
    e_value = quotient_reduce(e.subs(y, r))
    e_first = quotient_reduce(sp.diff(e, y).subs(y, r))
    sigma_value = quotient_reduce(sigma.subs(y, r))
    sigma_first = quotient_reduce(sp.diff(sigma, y).subs(y, r))
    sigma_second_unit = quotient_reduce(
        sp.diff(sigma, y, 2).subs(y, r)/2)

    inverse = sp.invert(sp.Poly(r+1, r), sp.Poly(Q_R, r)).as_expr()
    r_plus_one_is_unit = quotient_reduce((r+1)*inverse) == 1
    gamma_saturated = any(
        sp.expand(factor-gamma) == 0 for factor in state.saturation_factors)
    G_at_r = quotient_reduce(state.G.subs(y, r))
    G_at_r_saturated = any(
        sp.expand(factor-G_at_r) == 0
        for factor in state.saturation_factors)
    passed = (
        e_value == 0 and e_first != 0 and r_plus_one_is_unit
        and gamma_saturated and sigma_value == 0 and sigma_first == 0
        and sigma_second_unit == G_at_r and sigma_second_unit != 0
        and G_at_r_saturated)
    return ValuationCheck(
        z, e_value, e_first, sigma_value, sigma_first,
        sigma_second_unit, r_plus_one_is_unit, gamma_saturated,
        G_at_r_saturated, passed)


def check_plain_refinement(
    *, groebner_timeout: float = GROEBNER_TIMEOUT,
) -> RefinementCheck:
    '''Compare the q-supported z=0 window with the landed 34-poly residual.'''
    started = time.monotonic()
    plain = landed.run_r9(
        0, target_count=2, time_budget=groebner_timeout)
    if (plain.verdict != 'REDUCED' or plain.consumed != 2
            or len(plain.groebner_basis) != PLAIN_R9_Z0_BASIS_SIZE):
        raise AssertionError(
            'landed plain R9 z=0 residual no longer has 34 polynomials')

    plain_ansatz, _ = landed.build_r9_ansatz(0)
    plain_engine = landed.HighCoefficientEngine(
        plain_ansatz, start_degree=250, target_count=1,
        c=landed.DEFAULT_C)
    state = build_qsupport_ansatz(0)
    q_engine = landed.HighCoefficientEngine(
        state.ansatz, start_degree=250, target_count=1,
        c=landed.DEFAULT_C)

    e0, e1 = sp.symbols('e0 e1')
    s0, s1, s2 = sp.symbols('s0 s1 s2')
    g0 = state.g_coefficients[0]
    specialization = {
        e1: gamma, e0: -gamma*r,
        s2: g0, s1: -2*g0*r, s0: g0*r**2,
    }
    specialized_coefficient = quotient_reduce(
        plain_engine.master_coefficient(250).subs(specialization))
    qsupport_coefficient = quotient_reduce(
        q_engine.master_coefficient(250))
    matches = quotient_reduce(
        specialized_coefficient-qsupport_coefficient) == 0

    q_result = run_qsupport_state(
        0, start_degree=251, target_count=2,
        groebner_timeout=groebner_timeout)
    if q_result.timed_out:
        return RefinementCheck(
            PLAIN_R9_Z0_BASIS_SIZE, 0, 'REDUCED', matches,
            False, False, False, True, time.monotonic()-started)

    q_basis = sp.groebner(
        q_result.groebner_basis, *q_result.generators,
        order=ORDER, domain=sp.QQ)
    generators_verified = all(
        sp.cancel(q_basis.reduce(equation)[1]) == 0
        for equation in (*state.saturation_equations,
                          specialized_coefficient, Q_R))

    # Work over QQ(r) only for this non-membership check. The already-landed
    # 34 polynomials are a Groebner basis, so reconstruction is cheap. A
    # nonzero remainder for e0+e1*r proves that root support is a genuinely
    # new restriction, rather than merely a rewritten plain equation.
    plain_over_r = sp.groebner(
        plain.groebner_basis, *plain.generators, order=ORDER,
        domain=sp.QQ.frac_field(r))
    support_remainder = sp.cancel(
        plain_over_r.reduce(e0+e1*r)[1])
    support_is_new = support_remainder != 0
    genuine = matches and generators_verified and support_is_new
    return RefinementCheck(
        PLAIN_R9_Z0_BASIS_SIZE, len(q_result.groebner_basis),
        q_result.verdict, matches, generators_verified, support_is_new,
        genuine, False,
        time.monotonic()-started)


def _print_state(result: QSupportResult) -> None:
    if result.verdict == 'CONTRADICTION':
        label = ('CONTRADICTION (candidate NEW KILL; pending audit, '
                 'not yet verified independently)')
        first = 'n/a'
    else:
        label = 'REDUCED'
        first = str(result.first_unconsumed_degree)
    status = 'TIMEOUT' if result.timed_out else 'within budget'
    print(
        f'z={result.z}: {label}; basis={len(result.groebner_basis)}; '
        f'consumed={result.consumed}/{result.requested}; '
        f'first_unconsumed={first}; status={status}; '
        f'certificate={{1}}={result.unit_certificate_checked}; '
        f'reason={result.reason}; time={result.elapsed_seconds:.3f}s',
        flush=True)
    for step in result.steps:
        forced = ('' if step.substitution is None else
                  f'; {step.substitution[0]}={step.substitution[1]}')
        print(f'  y^{step.degree}: {step.action}; '
              f'basis={step.basis_size}{forced}', flush=True)


def main() -> None:
    overall_started = time.monotonic()
    print('q-support validation gate', flush=True)
    valuations = tuple(check_valuations(z) for z in range(7))
    for check in valuations:
        valuation_label = 'PASS' if check.passed else 'FAIL'
        print(
            f'  z={check.z}: v_r(e)=1 and v_r(sigma)=2: '
            f'{valuation_label}; '
            f'e(r)={check.e_value}; sigma(r)={check.sigma_value}; '
            f'sigma\'(r)={check.sigma_first_value}; '
            f'(sigma\'\'/2)(r)={check.sigma_second_unit_part}; '
            f'r+1 unit={check.r_plus_one_is_unit}; '
            f'gamma saturated={check.gamma_saturated}; '
            f'G(r) saturated={check.G_at_r_saturated}', flush=True)
    if not all(check.passed for check in valuations):
        raise AssertionError('r-adic valuation validation failed')

    refinement = check_plain_refinement()
    print(
        '  refinement vs landed plain R9 z=0: '
        f'plain_basis={refinement.plain_basis_size}; '
        f'qsupport_basis={refinement.qsupport_basis_size}; '
        f'qsupport_verdict={refinement.qsupport_verdict}; '
        f'coefficient_specialization_matches='
        f'{refinement.coefficient_specialization_matches}; '
        f'qsupport_generators_verified='
        f'{refinement.qsupport_generators_verified}; '
        f'support_relation_is_new={refinement.support_relation_is_new}; '
        f'genuine_refinement={refinement.genuine_refinement}; '
        f'timeout={refinement.timed_out}; '
        f'time={refinement.elapsed_seconds:.3f}s', flush=True)
    if not refinement.genuine_refinement:
        raise AssertionError('q-support refinement validation failed')

    print('main z-sweep (degrees 251 downward; cap 12)', flush=True)
    results: list[QSupportResult] = []
    for z in range(7):
        print(f'starting z={z}', flush=True)
        result = run_qsupport_state(z)
        results.append(result)
        _print_state(result)

    full_column = all(
        result.verdict == 'CONTRADICTION'
        and result.unit_certificate_checked and not result.timed_out
        for result in results)
    cutoffs = sum(result.timed_out for result in results)
    full_column_label = 'YES' if full_column else 'NO'
    print(
        'FULL-COLUMN CONTRADICTION ACROSS z=0..6: '
        f'{full_column_label}', flush=True)
    print(f'total_runtime={time.monotonic()-overall_started:.3f}s; '
          f'budget_cutoffs={cutoffs}', flush=True)


if __name__ == '__main__':
    main()
