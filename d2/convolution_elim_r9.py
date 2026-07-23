#!/usr/bin/env python3
'''Focused elimination experiments for the q-supported R9 z=0 state.

This client reuses ``build_qsupport_ansatz``, ``quotient_reduce``, and the
landed ``HighCoefficientEngine``. It attacks the landed degree-250 ideal for
e=gamma*(y+1)^9*(y-r), sigma=g0*(y-r)^2, q(r)=0, gamma*g0 != 0.
Groebner and resultant bottlenecks run in killable subprocesses on Windows.

IMPORTANT: strategy 3's specializations are always labelled
``HEURISTIC GUIDANCE ONLY -- NOT A PROOF``.
'''

from __future__ import annotations

from dataclasses import dataclass
import argparse
import pickle
import subprocess
import sys
import time
from typing import Iterable, Mapping, Sequence

import sympy as sp

import convolution_elim as landed
import convolution_elim_qsupport as qs

r, gamma, Q_R = qs.r, qs.gamma, qs.Q_R
START_DEGREE, LAST_DEGREE = 251, 244
STRATEGY_BUDGET = 300.0
SHALLOW_GROEBNER_BUDGET = 900.0
BASELINE_BUDGET = 120.0


class BudgetExpired(RuntimeError):
    pass


@dataclass(frozen=True)
class Attempt:
    label: str
    status: str
    elapsed_seconds: float
    detail: str
    basis: tuple[sp.Expr, ...] = ()
    generators: tuple[sp.Symbol, ...] = ()


@dataclass(frozen=True)
class StrategyReport:
    number: int
    name: str
    outcome: str
    reason: str
    elapsed_seconds: float
    attempts: tuple[Attempt, ...]
    best_basis: tuple[sp.Expr, ...]
    best_generators: tuple[sp.Symbol, ...]
    proof_eligible: bool = True


def _basis_exprs(basis: sp.GroebnerBasis) -> tuple[sp.Expr, ...]:
    return tuple(poly.as_expr() for poly in basis.polys)


def _normalize(expr: sp.Expr, generators: Sequence[sp.Symbol]) -> sp.Expr:
    expr = sp.cancel(expr)
    return sp.Integer(0) if expr == 0 else landed._normalize(
        expr, tuple(generators), sp.QQ)


def _task_worker_main() -> None:
    kind, payload = pickle.loads(sys.stdin.buffer.read())
    try:
        if kind != 'resultant':
            raise ValueError(f'unknown worker task {kind!r}')
        numerator = sp.fraction(sp.cancel(payload))[0]
        answer = sp.cancel(sp.resultant(Q_R, numerator, r))
        result = ('ok', answer)
    except BaseException as error:
        result = ('error', type(error).__name__, str(error))
    sys.stdout.buffer.write(pickle.dumps(result))


def _timed_resultant(equation: sp.Expr, timeout: float) -> sp.Expr:
    if timeout <= 0:
        raise BudgetExpired
    process = subprocess.Popen(
        [sys.executable, '-B', __file__, '--task-worker'], stdin=subprocess.PIPE,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    try:
        stdout, stderr = process.communicate(
            pickle.dumps(('resultant', equation)), timeout=timeout)
    except subprocess.TimeoutExpired:
        process.kill()
        process.communicate()
        raise BudgetExpired
    if process.returncode:
        raise RuntimeError(
            f'resultant worker exited {process.returncode}: {stderr!r}')
    result = pickle.loads(stdout)
    if result[0] == 'error':
        raise RuntimeError(f'resultant worker {result[1]}: {result[2]}')
    return result[1]


def _timed_basis(equations: Iterable[sp.Expr],
                 generators: Sequence[sp.Symbol], *, order: str,
                 timeout: float) -> sp.GroebnerBasis:
    generators = tuple(generators)
    normalized = tuple(_normalize(eq, generators) for eq in equations
                       if sp.cancel(eq) != 0)
    try:
        return landed._timed_groebner(
            normalized, generators, (), order, timeout)
    except landed._BudgetExpired as error:
        raise BudgetExpired from error


def _is_one(basis: sp.GroebnerBasis) -> bool:
    return len(basis.polys) == 1 and basis.polys[0].as_expr() == 1


def _remaining(deadline: float) -> float:
    return max(0.0, deadline-time.monotonic())


def resultant_self_test() -> bool:
    '''Check Res_r(q(r),(u+1)(r-t))=(u+1)^4 q(t).'''
    started = time.monotonic()
    t, u = sp.symbols('t u')
    toy_relation = (u+1)*(r-t)
    actual = _timed_resultant(toy_relation, 30.0)
    expected = (u+1)**4*Q_R.subs(r, t)
    factor_check = sp.factor(actual-expected) == 0
    substitution_check = sp.cancel(
        expected/(u+1)**4-Q_R.subs(r, t)) == 0
    passed = factor_check and substitution_check
    label = 'PASS' if passed else 'FAIL'
    print(f'SELF-TEST: {label}', flush=True)
    print(f'  toy_relation={toy_relation}', flush=True)
    print(f'  resultant_factored={sp.factor(actual)}', flush=True)
    print(f'  expected_factored={sp.factor(expected)}', flush=True)
    print(f'  factor_check={factor_check}; '
          f'substitute_r=t_check={substitution_check}; '
          f'time={time.monotonic()-started:.3f}s', flush=True)
    return passed


class R9Problem:
    '''Shared exact ansatz, coefficient window, and saturation data.'''

    def __init__(self) -> None:
        self.state = qs.build_qsupport_ansatz(0)
        self.g0 = self.state.g_coefficients[0]
        self.a = tuple(sp.Symbol(f'a{i}') for i in range(5))
        if not set(self.a).issubset(set(self.state.ansatz.unknowns)):
            raise AssertionError('unexpected d2 coefficient names')
        self.original_generators = tuple(sorted(
            set(self.state.ansatz.unknowns)
            | set(self.state.saturation_symbols), key=sp.default_sort_key))
        self.original_fixed = (Q_R, *self.state.saturation_equations)
        self.rfree_saturation = self.state.saturation_equations
        self.rfree_generators = tuple(
            symbol for symbol in self.original_generators if symbol != r)
        self.engine = landed.HighCoefficientEngine(
            self.state.ansatz, start_degree=START_DEGREE,
            target_count=START_DEGREE-LAST_DEGREE+1,
            c=landed.DEFAULT_C)
        self.coefficients = {
            degree: qs.quotient_reduce(self.engine.master_coefficient(degree))
            for degree in range(START_DEGREE, LAST_DEGREE-1, -1)}
        self.resultant_cache: dict[int, sp.Expr] = {}
        self.baseline_basis: tuple[sp.Expr, ...] = ()

    def original_equations(
        self, through_degree: int,
        substitutions: Mapping[sp.Symbol, sp.Expr] | None = None,
    ) -> tuple[sp.Expr, ...]:
        substitutions = {} if substitutions is None else substitutions
        equations = [sp.cancel(Q_R.subs(substitutions))]
        for equation in self.state.saturation_equations:
            value = qs.quotient_reduce(equation.subs(substitutions))
            if value != 0:
                equations.append(value)
        for degree in range(250, through_degree-1, -1):
            value = qs.quotient_reduce(
                self.coefficients[degree].subs(substitutions))
            if value != 0:
                equations.append(value)
        return tuple(equations)

    def resultant(self, degree: int, timeout: float) -> sp.Expr:
        if degree not in self.resultant_cache:
            raw = _timed_resultant(self.coefficients[degree], timeout)
            self.resultant_cache[degree] = _normalize(
                raw, self.rfree_generators)
        return self.resultant_cache[degree]


def print_basis(label: str, basis: Sequence[sp.Expr],
                generators: Sequence[sp.Symbol]) -> None:
    print(f'{label}: size={len(basis)}; '
          f'generators={tuple(map(str, generators))}', flush=True)
    if not basis:
        print('    <no completed basis>', flush=True)
    for index, polynomial in enumerate(basis, 1):
        print(f'    B{index} = {polynomial}', flush=True)


def reconstruct_baseline(problem: R9Problem,
                         timeout: float) -> sp.GroebnerBasis:
    started = time.monotonic()
    basis = _timed_basis(
        problem.original_equations(250), problem.original_generators,
        order='grevlex', timeout=timeout)
    expressions = _basis_exprs(basis)
    print('LANDED DEGREE-250 BASELINE: '
          f'basis={len(expressions)}; expected=5; '
          f'match={len(expressions) == 5}; '
          f'time={time.monotonic()-started:.3f}s', flush=True)
    print_basis('  baseline basis', expressions,
                problem.original_generators)
    if len(expressions) != 5:
        raise AssertionError('degree-250 q-supported basis is no longer size 5')
    problem.baseline_basis = expressions
    return basis


def _audit_contradiction(problem: R9Problem, degrees: Sequence[int],
                         timeout: float) -> tuple[bool, str]:
    '''Require the original, not resultant-only, generators to yield {1}.'''
    if not degrees:
        return False, 'no master coefficients were present to audit'
    check_degrees = tuple(degrees[:2])
    originals_match = all(
        qs.quotient_reduce(problem.engine.master_coefficient(degree))
        == problem.coefficients[degree] for degree in check_degrees)
    if not originals_match:
        return False, 'original master-coefficient substitution check failed'
    try:
        direct = _timed_basis(
            problem.original_equations(min(degrees)),
            problem.original_generators, order='grevlex', timeout=timeout)
    except BudgetExpired:
        return False, 'direct original-generator audit timed out'
    if not _is_one(direct):
        return False, 'direct original-generator audit did not return {1}'
    return True, (f'original coefficients {check_degrees} re-fetched exactly; '
                  'direct original-generator basis is literally {1}')


def strategy_1(problem: R9Problem, budget: float) -> StrategyReport:
    started, attempts = time.monotonic(), []
    deadline = started+budget
    equations = list(problem.rfree_saturation)
    best = _timed_basis(equations, problem.rfree_generators,
                        order='grevlex', timeout=_remaining(deadline))
    used_degrees: list[int] = []
    reason, outcome = '', 'BUDGET EXHAUSTED'
    for degree in range(250, LAST_DEGREE-1, -1):
        attempt_started = time.monotonic()
        try:
            consequence = problem.resultant(degree, _remaining(deadline))
            if consequence == 0:
                attempts.append(Attempt(
                    f'degree {degree}', 'IDENTITY',
                    time.monotonic()-attempt_started,
                    'r-resultant is zero; no equation added',
                    _basis_exprs(best), problem.rfree_generators))
                continue
            trial = _timed_basis(
                (*equations, consequence), problem.rfree_generators,
                order='grevlex', timeout=_remaining(deadline))
        except BudgetExpired:
            reason = f'wall budget expired during degree {degree}'
            attempts.append(Attempt(
                f'degree {degree}', 'TIMEOUT',
                time.monotonic()-attempt_started, reason,
                _basis_exprs(best), problem.rfree_generators))
            break
        equations.append(consequence)
        used_degrees.append(degree)
        best = trial
        status = 'CONTRADICTION' if _is_one(trial) else 'PROPER'
        term_count = len(sp.Poly(
            consequence, *problem.rfree_generators).terms())
        attempts.append(Attempt(
            f'degree {degree}', status,
            time.monotonic()-attempt_started,
            f'resultant terms={term_count}; basis={len(trial.polys)}',
            _basis_exprs(trial), problem.rfree_generators))
        print(f'  strategy 1 degree {degree}: {status}; '
              f'basis={len(trial.polys)}; '
              f'elapsed={time.monotonic()-started:.3f}s', flush=True)
        if _is_one(trial):
            audit_ok, detail = _audit_contradiction(
                problem, used_degrees, max(BASELINE_BUDGET, budget))
            if audit_ok:
                outcome = 'CONTRADICTION'
                reason = 'PENDING AUDIT -- basis is {1}; '+detail
            else:
                reason = 'candidate {1} rejected by verification gate: '+detail
            break
    else:
        reason = (f'collected through degree {LAST_DEGREE} without a '
                  'contradiction or exhibited exact solution')
    return StrategyReport(
        1, 'ELIMINATE r FIRST', outcome, reason,
        time.monotonic()-started, tuple(attempts), _basis_exprs(best),
        problem.rfree_generators)


def strategy_2(problem: R9Problem, budget: float) -> StrategyReport:
    started, attempts = time.monotonic(), []
    u0, u1 = problem.state.saturation_symbols
    a0, a1, a2, a3, a4 = problem.a
    orderings = (
        (problem.g0, gamma, a0, a1, a2, a3, a4, r, u0, u1),
        (problem.g0, gamma, a4, a3, a2, a1, a0, r, u0, u1),
        (problem.g0, gamma, a0, a2, a4, a1, a3, r, u0, u1),
    )
    equations = problem.original_equations(249)
    best_basis = problem.baseline_basis
    best_generators = problem.original_generators
    completed_ordering = False
    outcome = 'BUDGET EXHAUSTED'
    reason = 'all lex ordering experiments exhausted their attack budget'
    slice_budget = budget/len(orderings)
    for index, generators in enumerate(orderings, 1):
        attempt_started = time.monotonic()
        label = f'lex ordering {index}: {tuple(map(str, generators))}'
        try:
            basis = _timed_basis(
                equations, generators, order='lex', timeout=slice_budget)
        except BudgetExpired:
            attempts.append(Attempt(
                label, 'TIMEOUT', time.monotonic()-attempt_started,
                f'hard timeout={slice_budget:.3f}s'))
            print(f'  strategy 2 ordering {index}: TIMEOUT; '
                  f'time={time.monotonic()-attempt_started:.3f}s', flush=True)
            continue
        expressions = _basis_exprs(basis)
        if not completed_ordering or len(expressions) < len(best_basis):
            best_basis, best_generators = expressions, generators
            completed_ordering = True
        status = 'CONTRADICTION' if _is_one(basis) else 'PROPER'
        attempts.append(Attempt(label, status,
                                time.monotonic()-attempt_started,
                                f'basis={len(expressions)}', expressions,
                                generators))
        print(f'  strategy 2 ordering {index}: {status}; '
              f'basis={len(expressions)}; '
              f'time={time.monotonic()-attempt_started:.3f}s', flush=True)
        if _is_one(basis):
            audit_ok, detail = _audit_contradiction(
                problem, [250, 249], max(BASELINE_BUDGET, budget))
            if audit_ok:
                outcome = 'CONTRADICTION'
                reason = 'PENDING AUDIT -- lex basis is {1}; '+detail
            else:
                reason = 'candidate {1} rejected by verification gate: '+detail
            best_basis, best_generators = expressions, generators
            break
    return StrategyReport(
        2, 'WEIGHT-ORDERED elimination', outcome, reason,
        time.monotonic()-started, tuple(attempts), best_basis,
        best_generators)


def strategy_3(problem: R9Problem, budget: float) -> StrategyReport:
    started, attempts = time.monotonic(), []
    # HEURISTIC GUIDANCE ONLY -- NOT A PROOF. Neither gamma=1 nor g0=1 is
    # justified by a scale action at fixed c.
    specializations = (
        ('gamma=1', {gamma: sp.Integer(1)}),
        ('gamma=1, g0=1',
         {gamma: sp.Integer(1), problem.g0: sp.Integer(1)}),
    )
    slice_budget = budget/len(specializations)
    best_basis: tuple[sp.Expr, ...] = ()
    best_generators: tuple[sp.Symbol, ...] = ()
    saw_one = False
    print('  HEURISTIC GUIDANCE ONLY -- NOT A PROOF', flush=True)
    for label, substitutions in specializations:
        attempt_started = time.monotonic()
        generators = tuple(symbol for symbol in problem.original_generators
                           if symbol not in substitutions)
        equations = problem.original_equations(249, substitutions)
        try:
            basis = _timed_basis(equations, generators, order='grevlex',
                                 timeout=slice_budget)
        except BudgetExpired:
            status = 'TIMEOUT -- HEURISTIC GUIDANCE ONLY -- NOT A PROOF'
            attempts.append(Attempt(label, status,
                                    time.monotonic()-attempt_started,
                                    f'hard timeout={slice_budget:.3f}s'))
            print(f'  {label}: {status}; '
                  f'time={time.monotonic()-attempt_started:.3f}s', flush=True)
            continue
        expressions = _basis_exprs(basis)
        if not best_basis or len(expressions) < len(best_basis):
            best_basis, best_generators = expressions, generators
        if _is_one(basis):
            saw_one = True
            status = ('SPECIALIZED CONTRADICTION -- HEURISTIC GUIDANCE ONLY '
                      '-- NOT A PROOF')
        else:
            status = 'PROPER -- HEURISTIC GUIDANCE ONLY -- NOT A PROOF'
        attempts.append(Attempt(label, status,
                                time.monotonic()-attempt_started,
                                f'basis={len(expressions)}', expressions,
                                generators))
        print(f'  {label}: {status}; basis={len(expressions)}; '
              f'time={time.monotonic()-attempt_started:.3f}s', flush=True)
    reason = (('a specialized fiber has basis {1}; this does not imply the '
               'general state is contradictory') if saw_one else
              ('no specialized contradiction completed; no exact general '
               'solution was exhibited'))
    return StrategyReport(
        3, 'SUCCESSIVE SPECIALIZATION SANITY', 'BUDGET EXHAUSTED',
        reason, time.monotonic()-started, tuple(attempts), best_basis,
        best_generators, proof_eligible=False)


def strategy_4(problem: R9Problem, resultant_budget: float,
               groebner_budget: float) -> StrategyReport:
    started, attempts = time.monotonic(), []
    equations = list(problem.rfree_saturation)
    used_degrees: list[int] = []
    # No intermediate Groebner reduction occurs in this collection loop.
    for degree in range(250, LAST_DEGREE-1, -1):
        attempt_started = time.monotonic()
        try:
            consequence = problem.resultant(degree, resultant_budget)
        except BudgetExpired:
            reason = f'resultant collection timed out at degree {degree}'
            attempts.append(Attempt(f'collect degree {degree}', 'TIMEOUT',
                                    time.monotonic()-attempt_started, reason))
            return StrategyReport(
                4, 'MORE COEFFICIENTS SHALLOW', 'BUDGET EXHAUSTED', reason,
                time.monotonic()-started, tuple(attempts),
                problem.baseline_basis, problem.original_generators)
        if consequence != 0:
            equations.append(consequence)
            used_degrees.append(degree)
        attempts.append(Attempt(
            f'collect degree {degree}', 'COLLECTED (no Groebner)',
            time.monotonic()-attempt_started, f'zero={consequence == 0}'))
        print(f'  strategy 4 collected degree {degree}; '
              f'time={time.monotonic()-attempt_started:.3f}s', flush=True)
    attempt_started = time.monotonic()
    try:
        basis = _timed_basis(equations, problem.rfree_generators,
                             order='grevlex', timeout=groebner_budget)
    except BudgetExpired:
        reason = f'final Groebner call hit its {groebner_budget:.3f}s budget'
        attempts.append(Attempt('final degree 250..244 basis', 'TIMEOUT',
                                time.monotonic()-attempt_started, reason))
        return StrategyReport(
            4, 'MORE COEFFICIENTS SHALLOW', 'BUDGET EXHAUSTED', reason,
            time.monotonic()-started, tuple(attempts),
            problem.baseline_basis, problem.original_generators)
    expressions = _basis_exprs(basis)
    status = 'CONTRADICTION' if _is_one(basis) else 'PROPER'
    attempts.append(Attempt('final degree 250..244 basis', status,
                            time.monotonic()-attempt_started,
                            f'basis={len(expressions)}', expressions,
                            problem.rfree_generators))
    if _is_one(basis):
        audit_ok, detail = _audit_contradiction(
            problem, used_degrees, max(BASELINE_BUDGET, groebner_budget))
        if audit_ok:
            outcome, reason = 'CONTRADICTION', (
                'PENDING AUDIT -- final basis is {1}; '+detail)
        else:
            outcome, reason = 'BUDGET EXHAUSTED', (
                'candidate {1} rejected by verification gate: '+detail)
    else:
        outcome, reason = 'BUDGET EXHAUSTED', (
            'proper basis completed, but no exact solution point was '
            'exhibited and the coefficient window is exhausted')
    return StrategyReport(
        4, 'MORE COEFFICIENTS SHALLOW', outcome, reason,
        time.monotonic()-started, tuple(attempts), expressions,
        problem.rfree_generators)


def print_report(report: StrategyReport) -> None:
    outcome = (report.outcome+' -- PENDING AUDIT'
               if report.outcome == 'CONTRADICTION'
               and report.proof_eligible else report.outcome)
    print(f'STRATEGY {report.number} -- {report.name}: {outcome}; '
          f'time={report.elapsed_seconds:.3f}s', flush=True)
    print(f'  reason={report.reason}', flush=True)
    print(f'  proof_eligible={report.proof_eligible}', flush=True)
    for attempt in report.attempts:
        print(f'  attempt={attempt.label}; status={attempt.status}; '
              f'time={attempt.elapsed_seconds:.3f}s; '
              f'detail={attempt.detail}', flush=True)
    print_basis('  best/smallest completed Groebner basis',
                report.best_basis, report.best_generators)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--strategy-budget', type=float,
                        default=STRATEGY_BUDGET)
    parser.add_argument('--shallow-budget', type=float,
                        default=SHALLOW_GROEBNER_BUDGET)
    parser.add_argument('--baseline-budget', type=float,
                        default=BASELINE_BUDGET)
    parser.add_argument('--resultant-budget', type=float,
                        default=STRATEGY_BUDGET)
    args = parser.parse_args(argv)
    for name, value in vars(args).items():
        if value <= 0:
            parser.error(f'{name} must be positive')
    overall_started = time.monotonic()
    print('R9 z=0 q-support focused elimination', flush=True)
    print('budgets: '
          f'strategies1-3={args.strategy_budget:.3f}s; '
          f'strategy4_final={args.shallow_budget:.3f}s; '
          f'per_resultant={args.resultant_budget:.3f}s', flush=True)

    # Mandatory gate: real master coefficients are not constructed first.
    try:
        self_test_ok = resultant_self_test()
    except BaseException as error:
        print(f'SELF-TEST: FAIL; {type(error).__name__}: {error}', flush=True)
        print('STOPPED before all real strategies.', flush=True)
        return 2
    if not self_test_ok:
        print('STOPPED before all real strategies.', flush=True)
        return 2

    setup_started = time.monotonic()
    problem = R9Problem()
    nonzero_ok = (problem.state.saturation_factors == (gamma, problem.g0)
                  and len(problem.state.saturation_equations) == 2)
    print('REAL SETUP: '
          f'unknowns={tuple(map(str, problem.state.ansatz.unknowns))}; '
          f'nonzero_factors={problem.state.saturation_factors}; '
          f'saturation_check={nonzero_ok}; '
          f'time={time.monotonic()-setup_started:.3f}s', flush=True)
    if not nonzero_ok:
        raise AssertionError('unexpected z=0 saturation structure')
    reconstruct_baseline(problem, args.baseline_budget)

    reports: list[StrategyReport] = []
    runners = (
        ('STRATEGY 1 START -- ELIMINATE r FIRST',
         lambda: strategy_1(problem, args.strategy_budget)),
        ('STRATEGY 2 START -- WEIGHT-ORDERED elimination',
         lambda: strategy_2(problem, args.strategy_budget)),
        ('STRATEGY 3 START -- SUCCESSIVE SPECIALIZATION SANITY -- '
         'HEURISTIC GUIDANCE ONLY -- NOT A PROOF',
         lambda: strategy_3(problem, args.strategy_budget)),
        ('STRATEGY 4 START -- MORE COEFFICIENTS SHALLOW',
         lambda: strategy_4(problem, args.resultant_budget,
                            args.shallow_budget)),
    )
    for heading, runner in runners:
        print(heading, flush=True)
        attempt_started = time.monotonic()
        try:
            report = runner()
        except BaseException as error:
            number = len(reports)+1
            report = StrategyReport(
                number, heading.split('--', 1)[-1].strip(),
                'BUDGET EXHAUSTED',
                f'FAILED honestly with {type(error).__name__}: {error}',
                time.monotonic()-attempt_started, (), (), ())
        reports.append(report)
        print_report(report)

    proof_kills = [x for x in reports if x.outcome == 'CONTRADICTION'
                   and x.proof_eligible]
    exact_solutions = [x for x in reports if x.outcome == 'SOLUTION'
                       and x.proof_eligible]
    if proof_kills:
        verdict = ('CONTRADICTION -- PENDING AUDIT: z=0 is a candidate NEW '
                   'KILL with a literal {1} basis and original-generator gate')
    elif exact_solutions:
        verdict = 'SOLUTION: z=0 survives at the printed exact point'
    else:
        verdict = ('UNRESOLVED: z=0 is neither killed nor supplied with an '
                   'exact solution by these bounded attacks')
    print('OVERALL VERDICT: '+verdict, flush=True)
    print('RECOMMENDATION: audit any PENDING AUDIT certificate independently; '
          'otherwise move the r-free degree-250..244 system to a faster CAS '
          '(Singular/Magma), preserving both nonzero saturations.', flush=True)
    print(f'TOTAL TIME: {time.monotonic()-overall_started:.3f}s', flush=True)
    return 0


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '--task-worker':
        _task_worker_main()
    else:
        raise SystemExit(main())
