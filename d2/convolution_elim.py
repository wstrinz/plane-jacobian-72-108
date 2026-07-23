#!/usr/bin/env python3
'''Exact coefficient-ideal elimination for the f31 convolution identity.

This module upgrades :mod:`convolution_descent` without changing it. It walks
master coefficients downward, keeps the landed forced-square reductions, and
otherwise adjoins each coefficient to an ideal. A Groebner basis over QQ (or
over the rational-function field in declared parameters) is recomputed after
every consumed coefficient. Thus CONTRADICTION means that the saturated ideal
contains 1; REDUCED means the last ideal is proper and the reported residual
system remains.

Leading coefficients are made units by Rabinowitsch equations u_i*x_i-1.
Declared parameters are units in the coefficient field. No shipped run uses a
scale gauge: after c=-1/6630 is fixed, the natural weighted scaling also
rescales c, so setting gamma or e1 to 1 is unjustified. Both are retained,
which costs a larger ring.

The optional gauge_fix API is strict. The caller must supply a nonempty
gauge_justification asserting a genuine G_m action preserving the identity
and all hypotheses, with every nonzero orbit meeting x=1. Only then does
substitution x=1 lose no nonzero solutions. Cosmetic homogeneity is not
enough.

Arithmetic is exact. Groebner updates run in killable subprocesses so the
budget is hard at those bottlenecks; coefficient construction is checked
between exact operations. Every return reports the exact number of
coefficients consumed.
'''

from __future__ import annotations

from dataclasses import dataclass
import pickle
import subprocess
import sys
import time
from typing import Iterable, Mapping, Sequence

import sympy as sp

import convolution_descent as cd

y = cd.y
DEFAULT_C = -sp.Rational(1, 6630)


def _tmul(left: Mapping[int, sp.Expr], right: Mapping[int, sp.Expr],
          cutoff: int) -> dict[int, sp.Expr]:
    '''Multiply reverse/codegree polynomials through the requested cutoff.'''
    out: dict[int, sp.Expr] = {}
    for i, a in left.items():
        for j, b in right.items():
            if i+j <= cutoff:
                out[i+j] = out.get(i+j, sp.Integer(0)) + a*b
    return {i: a for i, a in out.items() if a != 0}


def _tpow(poly: Mapping[int, sp.Expr], exponent: int,
          cutoff: int) -> dict[int, sp.Expr]:
    if exponent < 0:
        raise ValueError('exponent must be nonnegative')
    result: dict[int, sp.Expr] = {0: sp.Integer(1)}
    factor, remaining = dict(poly), exponent
    while remaining:
        if remaining & 1:
            result = _tmul(result, factor, cutoff)
        remaining >>= 1
        if remaining:
            factor = _tmul(factor, factor, cutoff)
    return result


def _reverse(poly: Mapping[int, sp.Expr], maximum: int) -> dict[int, sp.Expr]:
    return {maximum-degree: coefficient for degree, coefficient in poly.items()
            if degree <= maximum}


class HighCoefficientEngine:
    '''Exact high-coefficient window using the landed convolution data.

    The landed engine constructs complete intermediate polynomials. For the
    degree-250 R9 systems that is needlessly expensive. This class uses its
    parsed h_f formulas and Ansatz substitutions, but reverses every factor
    and retains only requested codegrees. Every f=0,...,7 is still summed.
    '''

    def __init__(self, ansatz: cd.Ansatz, *, start_degree: int,
                 target_count: int, c: sp.Expr = DEFAULT_C,
                 h: Mapping[int, sp.Expr] | None = None) -> None:
        if target_count <= 0:
            raise ValueError('target_count must be positive')
        self.start_degree = start_degree
        self.floor = start_degree-target_count+1
        self.h = dict(cd.base.load_h() if h is None else h)
        if sorted(self.h) != list(range(8)):
            raise ValueError('expected h_0,...,h_7')
        phi = sp.sympify(c)*(y+1)**30*cd.base.q
        self._source = dict(ansatz.substitutions)
        self._source_max = {s: max(p) if p else 0
                            for s, p in self._source.items()}
        self._source_rev = {s: _reverse(p, self._source_max[s])
                            for s, p in self._source.items()}
        self._phi = cd.from_expr(phi)
        self._phi_max = max(self._phi)
        self._phi_rev = _reverse(self._phi, self._phi_max)
        self._e_max = self._source_max[cd.base.e]
        self._e_rev = self._source_rev[cd.base.e]
        self._terms: dict[int, tuple[int, dict[int, sp.Expr]]] = {}
        self._coefficients: dict[int, sp.Expr] = {}
        self._build_terms()

    def _build_h_reverse(self, f: int, cutoff: int
                         ) -> tuple[int, dict[int, sp.Expr]]:
        variables = (cd.base.d0, cd.base.d1, cd.base.d2, cd.base.e)
        terms = sp.Poly(self.h[f], *variables).terms()
        maxima = [sum(power*self._source_max[symbol]
                      for symbol, power in zip(variables, monomial))
                  for monomial, _ in terms]
        h_max = max(maxima)
        out: dict[int, sp.Expr] = {}
        cache: dict[tuple[sp.Symbol, int], dict[int, sp.Expr]] = {}
        for ((monomial, scalar), term_max) in zip(terms, maxima):
            shift = h_max-term_max
            if shift > cutoff:
                continue
            term: dict[int, sp.Expr] = {shift: scalar}
            for symbol, exponent in zip(variables, monomial):
                key = (symbol, exponent)
                if key not in cache:
                    cache[key] = _tpow(self._source_rev[symbol], exponent,
                                       cutoff)
                term = _tmul(term, cache[key], cutoff)
            for codegree, coefficient in term.items():
                out[codegree] = out.get(codegree, sp.Integer(0)) + coefficient
        return h_max, {i: a for i, a in out.items() if a != 0}

    def _build_terms(self) -> None:
        variables = (cd.base.d0, cd.base.d1, cd.base.d2, cd.base.e)
        for f in range(8):
            e_exponent = 21-3*f
            monomials = sp.Poly(self.h[f], *variables).monoms()
            h_max = max(sum(power*self._source_max[symbol]
                            for symbol, power in zip(variables, monomial))
                        for monomial in monomials)
            term_max = f*self._phi_max + e_exponent*self._e_max + h_max
            cutoff = term_max-self.floor
            if cutoff < 0:
                self._terms[f] = (term_max, {})
                continue
            _, h_rev = self._build_h_reverse(f, cutoff)
            phi_power = _tpow(self._phi_rev, f, cutoff)
            e_power = _tpow(self._e_rev, e_exponent, cutoff)
            reverse_term = _tmul(_tmul(phi_power, e_power, cutoff),
                                  h_rev, cutoff)
            self._terms[f] = (term_max, reverse_term)

    def master_coefficient(self, target: int) -> sp.Expr:
        if not self.floor <= target <= self.start_degree:
            raise ValueError('target lies outside the constructed window')
        if target not in self._coefficients:
            total = sp.Integer(0)
            for term_max, reverse_term in self._terms.values():
                total += reverse_term.get(term_max-target, sp.Integer(0))
            self._coefficients[target] = total
        return self._coefficients[target]


@dataclass(frozen=True)
class EliminationStep:
    degree: int
    action: str
    coefficient: sp.Expr
    reduced: sp.Expr
    substitution: tuple[sp.Symbol, sp.Expr] | None = None
    basis_size: int = 0


@dataclass(frozen=True)
class EliminationResult:
    verdict: str
    steps: tuple[EliminationStep, ...]
    substitutions: Mapping[sp.Symbol, sp.Expr]
    groebner_basis: tuple[sp.Expr, ...]
    residual_equations: tuple[sp.Expr, ...]
    generators: tuple[sp.Symbol, ...]
    consumed: int
    requested: int
    start_degree: int
    last_consumed_degree: int | None
    reason: str
    elapsed_seconds: float
    zero_dimensional: bool | None = None


def _is_one(basis: sp.GroebnerBasis) -> bool:
    return len(basis.polys) == 1 and basis.polys[0].as_expr() == 1


def _normalize(expr: sp.Expr, generators: Sequence[sp.Symbol], domain) -> sp.Expr:
    expr = sp.cancel(expr)
    if expr == 0:
        return sp.Integer(0)
    return sp.Poly(expr, *generators, domain=domain).monic().as_expr()


class _BudgetExpired(RuntimeError):
    pass


def _groebner_worker(payload) -> tuple:
    '''Subprocess worker: a stuck Buchberger computation is terminable.'''
    equations, generators, parameters, order = payload
    try:
        domain = sp.QQ if not parameters else sp.QQ.frac_field(*parameters)
        basis = sp.groebner(equations, *generators,
                            order=order, domain=domain)
        return ('ok', tuple(poly.as_expr() for poly in basis.polys))
    except BaseException as error:  # propagated with type/message to parent
        return ('error', type(error).__name__, str(error))


def _worker_main() -> None:
    payload = pickle.loads(sys.stdin.buffer.read())
    sys.stdout.buffer.write(pickle.dumps(_groebner_worker(payload)))


def _timed_groebner(equations: tuple[sp.Expr, ...],
                    generators: tuple[sp.Symbol, ...],
                    parameters: tuple[sp.Symbol, ...], order: str,
                    timeout: float) -> sp.GroebnerBasis:
    if timeout <= 0:
        raise _BudgetExpired
    expires = time.monotonic()+timeout
    process = subprocess.Popen(
        [sys.executable, __file__, '--groebner-worker'],
        stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    data = pickle.dumps((equations, generators, parameters, order))
    try:
        remaining = expires-time.monotonic()
        if remaining <= 0:
            process.kill()
            process.communicate()
            raise _BudgetExpired
        stdout, stderr = process.communicate(data, timeout=remaining)
    except subprocess.TimeoutExpired:
        process.kill()
        process.communicate()
        raise _BudgetExpired
    if process.returncode:
        raise RuntimeError(
            f'Groebner worker exited {process.returncode}: '
            f'{stderr!r}')
    payload = pickle.loads(stdout)
    if payload[0] == 'error':
        raise RuntimeError(f'Groebner worker {payload[1]}: {payload[2]}')
    domain = sp.QQ if not parameters else sp.QQ.frac_field(*parameters)
    # The returned polynomials already form a Groebner basis, so this
    # reconstruction is cheap and gives the parent the normal-form API.
    return sp.groebner(payload[1], *generators, order=order, domain=domain)


def eliminate(
    ansatz: cd.Ansatz, *, start_degree: int, target_count: int,
    nonzero: Iterable[sp.Symbol] = (), c: sp.Expr = DEFAULT_C,
    time_budget: float = 60.0,
    initial_substitutions: Mapping[sp.Symbol, sp.Expr] | None = None,
    gauge_fix: sp.Symbol | None = None,
    gauge_justification: str | None = None,
    order: str = 'grevlex',
) -> EliminationResult:
    '''Accumulate a saturated coefficient ideal over an exact domain.

    gauge_fix=x substitutes x=1 only when the caller supplies a nonempty
    gauge_justification certifying the G_m orbit argument in the module
    docstring. Otherwise every listed nonzero unknown gets its own
    Rabinowitsch variable. Parameters are already units in the fraction field.
    '''
    if target_count <= 0:
        raise ValueError('target_count must be positive')
    if time_budget <= 0:
        raise ValueError('time_budget must be positive')
    c = sp.sympify(c)
    parameters = tuple(sorted(ansatz.parameters, key=sp.default_sort_key))
    if c.free_symbols-set(parameters):
        raise ValueError('c may involve only declared parameters')
    substitutions = dict(initial_substitutions or {})
    if not set(substitutions).issubset(ansatz.unknowns):
        raise ValueError('initial substitutions include a non-unknown')
    nonzero_tuple = tuple(dict.fromkeys(nonzero))
    allowed = set(ansatz.unknowns) | set(parameters)
    if not set(nonzero_tuple).issubset(allowed):
        raise ValueError('nonzero constraints must be unknowns or parameters')
    if gauge_fix is not None:
        if gauge_fix not in nonzero_tuple:
            raise ValueError('gauge-fixed symbol must be declared nonzero')
        if gauge_fix not in ansatz.unknowns:
            raise ValueError('only an unknown can be gauge-fixed')
        if not gauge_justification or not gauge_justification.strip():
            raise ValueError('gauge_fix requires a G_m gauge justification')
        substitutions[gauge_fix] = sp.Integer(1)

    domain = sp.QQ if not parameters else sp.QQ.frac_field(*parameters)
    started = time.monotonic()
    deadline = started+time_budget
    coefficient_equations: list[sp.Expr] = []
    rab_symbols: dict[sp.Symbol, sp.Symbol] = {}
    occupied = set(ansatz.unknowns) | set(parameters)
    for index, symbol in enumerate(nonzero_tuple):
        if symbol in parameters or symbol == gauge_fix:
            continue
        name = f'_rab_{index}'
        while sp.Symbol(name) in occupied:
            name = '_'+name
        rab_symbols[symbol] = sp.Symbol(name)
        occupied.add(rab_symbols[symbol])

    def current_equations() -> list[sp.Expr]:
        equations = [sp.sympify(eq).subs(substitutions)
                     for eq in coefficient_equations]
        equations.extend((u*symbol-1).subs(substitutions)
                         for symbol, u in rab_symbols.items())
        return [sp.cancel(eq) for eq in equations if sp.cancel(eq) != 0]

    def make_basis(
        extra_symbols: Iterable[sp.Symbol] = (),
    ) -> tuple[sp.GroebnerBasis, tuple[sp.Symbol, ...], tuple[sp.Expr, ...]]:
        current = current_equations()
        free = (set().union(*(eq.free_symbols for eq in current))
                if current else set())
        free.update(extra_symbols)
        generators = tuple(sorted(free-set(parameters),
                                  key=sp.default_sort_key))
        if not generators:
            dummy = sp.Symbol('_unit_dummy')
            basis = sp.groebner([1 if current else 0], dummy,
                                order=order, domain=domain)
            return basis, (), tuple(current)
        normalized = tuple(_normalize(eq, generators, domain)
                           for eq in current)
        basis = _timed_groebner(
            normalized, generators, parameters, order,
            max(0.0, deadline-time.monotonic()))
        return basis, generators, normalized

    basis, generators, equations = make_basis()
    engine = HighCoefficientEngine(ansatz, start_degree=start_degree,
                                   target_count=target_count, c=c)
    steps: list[EliminationStep] = []

    def finish(verdict: str, reason: str) -> EliminationResult:
        basis_exprs = tuple(poly.as_expr() for poly in basis.polys)
        return EliminationResult(
            verdict, tuple(steps), dict(substitutions), basis_exprs,
            tuple(equations), tuple(generators), len(steps), target_count,
            start_degree, (steps[-1].degree if steps else None), reason,
            time.monotonic()-started, None)

    if _is_one(basis):
        return finish('CONTRADICTION',
                      'initial saturated ideal is the unit ideal')

    for target in range(start_degree, start_degree-target_count, -1):
        if time.monotonic()-started >= time_budget:
            return finish('REDUCED', 'TIME_BUDGET before next coefficient')
        coefficient = sp.cancel(
            engine.master_coefficient(target).subs(substitutions))
        if time.monotonic()-started >= time_budget:
            return finish(
                'REDUCED',
                f'TIME_BUDGET while computing unconsumed degree {target}')
        new_symbols = coefficient.free_symbols-set(parameters)-set(generators)
        if new_symbols:
            generators = tuple(sorted(set(generators) | new_symbols,
                                      key=sp.default_sort_key))
            basis = sp.groebner(
                tuple(poly.as_expr() for poly in basis.polys),
                *generators, order=order, domain=domain)
        reduced = (sp.cancel(basis.reduce(coefficient)[1])
                   if generators else coefficient)
        if reduced == 0:
            steps.append(EliminationStep(
                target, 'IDENTITY_MOD_IDEAL', coefficient, sp.Integer(0),
                basis_size=len(basis.polys)))
            continue

        active = tuple(symbol for symbol in ansatz.unknowns
                       if symbol not in substitutions and reduced.has(symbol))
        forced = cd.ConvolutionDescent._forced_square(reduced, active)
        old_basis, old_generators, old_equations = basis, generators, equations
        old_substitutions = dict(substitutions)
        old_equation_count = len(coefficient_equations)
        if forced is not None:
            substitutions[forced[0]] = forced[1]
            action = 'FORCED'
        else:
            coefficient_equations.append(coefficient)
            action = 'ADDED_TO_IDEAL'
        try:
            basis, generators, equations = make_basis()
        except _BudgetExpired:
            substitutions.clear()
            substitutions.update(old_substitutions)
            del coefficient_equations[old_equation_count:]
            basis, generators, equations = (
                old_basis, old_generators, old_equations)
            return finish(
                'REDUCED',
                f'TIME_BUDGET while testing unconsumed degree {target}')
        steps.append(EliminationStep(target, action, coefficient, reduced,
                                     forced, len(basis.polys)))
        if _is_one(basis):
            return finish('CONTRADICTION',
                          f'1 in saturated ideal after degree {target}')
    return finish('REDUCED', 'TARGET_LIMIT')


def _constant_t1_ansatz(
    *, sigma_degree: int | None = 3, prepared: bool = True,
) -> tuple[cd.Ansatz, dict[str, sp.Symbol]]:
    gamma = sp.Symbol('gamma', nonzero=True)
    symbols: dict[str, sp.Symbol] = {'gamma': gamma}
    if prepared:
        v0 = sp.Symbol('v0')
        symbols['v0'] = v0
        ansatz = cd.build_ansatz(
            degrees={'d2': 4, 'sigma': 3},
            d1=-gamma**4/(4096*DEFAULT_C)
               *(y**2+sp.Rational(25, 4)*y+v0),
            e=gamma*(y+1)**9, parameters=(gamma,))
        return ansatz, symbols
    if sigma_degree is None:
        raise ValueError('sigma_degree is required for an unprepared ansatz')
    ansatz = cd.build_ansatz(
        degrees={'d2': 4, 'd1': 2, 'sigma': sigma_degree},
        e=gamma*(y+1)**9,
        prefixes={'d2': 'a', 'd1': 'b', 'sigma': 's'})
    symbols['d1_lc'] = sp.Symbol('b2')
    symbols['sigma_lc'] = sp.Symbol(f's{sigma_degree}')
    return ansatz, symbols


def run_plain_gate() -> cd.DescentResult:
    '''Run the landed constant-E forced-chain gate at exact fixed c.'''
    ansatz, _ = _constant_t1_ansatz(prepared=True)
    result = cd.ConvolutionDescent(ansatz, c=DEFAULT_C).descend(238, 226)
    if (result.verdict, result.stopping_degree) != ('CONTRADICTION', 226):
        raise AssertionError('landed constant-E gate did not contradict at 226')
    return result


def run_tied_gate(sigma_degree: int, *, time_budget: float = 60.0
                  ) -> EliminationResult:
    '''Run an unprepared d=2 tied cell from its master coefficients only.'''
    if sigma_degree not in (4, 5):
        raise ValueError('the new-power gate covers sigma degrees 4 and 5')
    ansatz, names = _constant_t1_ansatz(sigma_degree=sigma_degree,
                                       prepared=False)
    count = 2 if sigma_degree == 5 else 4
    return eliminate(ansatz, start_degree=242, target_count=count,
                     nonzero=(names['d1_lc'], names['sigma_lc'],
                              names['gamma']),
                     c=DEFAULT_C, time_budget=time_budget)


def build_r9_ansatz(z: int) -> tuple[cd.Ansatz, tuple[sp.Symbol, sp.Symbol]]:
    '''Build the T5_T2_COLUMN.md (R9) pattern-B degree state.'''
    if z not in range(7):
        raise ValueError('z must be in 0,...,6')
    e0, e1 = sp.symbols('e0 e1')
    ansatz = cd.build_ansatz(
        degrees={'d2': 4, 'sigma': 2+z}, d1=sp.Integer(0),
        e=(y+1)**9*(e1*y+e0),
        prefixes={'d2': 'a', 'sigma': 's'})
    return ansatz, (e1, sp.Symbol(f's{2+z}'))


def run_r9(z: int, *, target_count: int = 10,
           time_budget: float = 90.0) -> EliminationResult:
    ansatz, leading = build_r9_ansatz(z)
    return eliminate(ansatz, start_degree=250, target_count=target_count,
                     nonzero=leading, c=DEFAULT_C,
                     time_budget=time_budget)


def _brief(label: str, result: EliminationResult) -> None:
    print(f'{label}: {result.verdict}; consumed {result.consumed}/'
          f'{result.requested}; last={result.last_consumed_degree}; '
          f'basis={len(result.groebner_basis)}; '
          f'zero_dimensional={result.zero_dimensional}; '
          f'reason={result.reason}; time={result.elapsed_seconds:.3f}s')
    for step in result.steps:
        print(f'  y^{step.degree}: {step.action}; basis={step.basis_size}')
        if step.substitution:
            print(f'    {step.substitution[0]} = {step.substitution[1]}')
    if result.verdict == 'REDUCED':
        print('  Groebner basis:')
        for polynomial in result.groebner_basis:
            print(f'    {polynomial}')


def main() -> None:
    print('convolution_elim exact gates; c = -1/6630')
    plain = run_plain_gate()
    print(f'plain gate: PASS; {plain.verdict} at degree '
          f'{plain.stopping_degree}; forced={len(plain.substitutions)}')
    for sigma_degree in (5, 4):
        result = run_tied_gate(sigma_degree)
        _brief(f'tied deg(sigma)={sigma_degree}', result)
        if result.verdict != 'CONTRADICTION':
            raise AssertionError(f'tied sigma degree {sigma_degree} not killed')
    if '--gates-only' in sys.argv:
        print('convolution_elim gates: PASS (R9 exploration skipped)')
        return
    for z in range(7):
        _brief(f'R9 z={z}', run_r9(z))


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '--groebner-worker':
        _worker_main()
    else:
        main()
