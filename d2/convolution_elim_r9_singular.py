'''Phase F3: saturated R9 residue ideals over Q(alpha) via Singular/WSL.

Equations come from the landed q-support ansatz and high-coefficient engine.
Singular receives program text on stdin because /mnt/c is unusable in WSL.
The C08 bridge gate must pass before any residue run is trusted.
'''

from __future__ import annotations

import argparse
from dataclasses import dataclass
import re
import subprocess
import sys
import time
from typing import Iterable, Sequence

import sympy as sp

import convolution_descent as cd
import convolution_elim as landed
import convolution_elim_qsupport as qsupport


START_DEGREE = 251
DEFAULT_COEFFICIENT_COUNT = 6
DEFAULT_TIMEOUT_SECONDS = 20 * 60.0
MINPOLY_SINGULAR = '2048*a^4-512*a^3+320*a^2-240*a+195'
WSL_COMMAND = (
    'wsl.exe', '-d', 'Ubuntu', '--', 'bash', '-lc',
    'cd $HOME && Singular -q',
)


@dataclass(frozen=True)
class SingularRun:
    label: str
    program: str
    status: str
    returncode: int | None
    stdout: str
    stderr: str
    elapsed_seconds: float
    saturation_ran: bool
    unit_ideal: bool | None
    dimension: int | None
    vector_space_dimension: int | None
    basis_text: str | None


@dataclass(frozen=True)
class GeneratedState:
    z: int
    degrees: tuple[int, ...]
    equations: tuple[sp.Expr, ...]
    singular_program: str
    saturation_factors: tuple[sp.Expr, ...]
    generation_seconds: float
    coefficient_crosscheck_passed: bool


def _decode_timeout_stream(value: str | bytes | None) -> str:
    if value is None:
        return ''
    if isinstance(value, bytes):
        return value.decode('utf-8', errors='replace')
    return value.replace('\x00', '')


def _clean_completed_stream(value: str) -> str:
    # WSL service failures are emitted as UTF-16-like text even when the
    # child process was requested in UTF-8 text mode.
    return value.replace('\x00', '')


def _marker_int(output: str, marker: str) -> int | None:
    match = re.search(re.escape(marker) + r'\s*\r?\n\s*(-?\d+)', output)
    return None if match is None else int(match.group(1))


def _marker_block(output: str, begin: str, end: str) -> str | None:
    match = re.search(
        re.escape(begin) + r'\s*\r?\n(.*?)\r?\n' + re.escape(end),
        output, flags=re.DOTALL)
    return None if match is None else match.group(1).strip()


def run_singular(label: str, program: str, timeout: float) -> SingularRun:
    '''Pipe one program to WSL Singular, retaining timeout partial output.'''
    started = time.monotonic()
    try:
        completed = subprocess.run(
            WSL_COMMAND, input=program, text=True, encoding='utf-8',
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            timeout=timeout, check=False)
        stdout = _clean_completed_stream(completed.stdout)
        stderr = _clean_completed_stream(completed.stderr)
        returncode: int | None = completed.returncode
        status = 'completed' if completed.returncode == 0 else 'process_error'
    except subprocess.TimeoutExpired as exc:
        stdout = _decode_timeout_stream(exc.stdout)
        stderr = _decode_timeout_stream(exc.stderr)
        returncode = None
        status = 'timeout'
    elapsed = time.monotonic() - started
    combined = stdout + '\n' + stderr
    sat_ran = '@@SATURATION_RAN' in combined
    unit_value = _marker_int(combined, '@@UNIT_IDEAL')
    dimension = _marker_int(combined, '@@DIMENSION')
    vdim = _marker_int(combined, '@@VDIM')
    basis = _marker_block(combined, '@@BASIS_BEGIN', '@@BASIS_END')
    return SingularRun(
        label, program, status, returncode, stdout, stderr, elapsed, sat_ran,
        None if unit_value is None else bool(unit_value), dimension, vdim,
        basis)


def c08_program() -> str:
    template = '''LIB @DQ@elim.lib@DQ@;
ring R = (0,a),(X,D,E),dp;
minpoly = 2048*a^4-512*a^3+320*a^2-240*a+195;
ideal I = 6*X^2*D^2-9*X*D*E-E^2;
poly nonzero_product = X*D*E;
list SatData = sat(I,nonzero_product);
ideal Isat = SatData[1];
int saturation_exponent = SatData[2];
@DQ@@@SATURATION_RAN@DQ@;
saturation_exponent;
ideal G = std(Isat);
int is_unit = (reduce(1,G)==0);
@DQ@@@UNIT_IDEAL@DQ@;
is_unit;
@DQ@@@DIMENSION@DQ@;
dim(G);
@DQ@@@BASIS_BEGIN@DQ@;
G;
@DQ@@@BASIS_END@DQ@;
quit;
'''
    return template.replace('@DQ@', chr(34))


def validate_c08_extension_point() -> None:
    '''Check the advertised torus point over a quadratic extension formally.'''
    t = sp.Symbol('t')
    relation = t**2 + 9*t - 6
    value = 6 - 9*t - t**2
    assert sp.rem(sp.Poly(value, t), sp.Poly(relation, t)).is_zero
    assert relation.subs(t, 0) != 0


def validate_c08_run(run: SingularRun) -> tuple[bool, str]:
    validate_c08_extension_point()
    if run.status != 'completed':
        return False, f'bridge process status is {run.status}'
    if not run.saturation_ran:
        return False, 'saturation marker was not emitted'
    if run.unit_ideal is not False:
        return False, f'expected a proper ideal, got unit={run.unit_ideal}'
    if run.dimension is None or run.dimension < 1:
        return False, f'expected positive-dimensional proper ideal, dim={run.dimension}'
    return True, ('proper saturated ideal; (X,D,E)=(1,1,t), '
                  't^2+9*t-6=0 is a torus point over an extension of Q(alpha)')


def _singular_expression(expression: sp.Expr) -> str:
    '''Translate a quotient-reduced SymPy polynomial to Singular syntax.'''
    alpha = sp.Symbol('a')
    gm = sp.Symbol('gm')
    mapped = sp.sympify(expression).subs({qsupport.r: alpha,
                                         qsupport.gamma: gm})
    text = sp.sstr(mapped).replace('**', '^')
    if 'r' in {symbol.name for symbol in mapped.free_symbols}:
        raise AssertionError('the quotient generator r was not mapped to a')
    return text


def _state_variables(state: qsupport.QSupportState) -> tuple[str, ...]:
    names = [symbol.name for symbol in state.g_coefficients]
    names.append('gm')
    names.extend(('a4', 'a3', 'a2', 'a1', 'a0'))
    if len(names) != len(set(names)):
        raise AssertionError('duplicate Singular variable name')
    return tuple(names)


def _expected_saturation_factors(
    state: qsupport.QSupportState,
) -> tuple[sp.Expr, ...]:
    G_at_alpha = qsupport.quotient_reduce(state.G.subs(cd.y, qsupport.r))
    return qsupport._unique_expressions(
        (qsupport.gamma, state.g_coefficients[-1], G_at_alpha))


def _audit_ansatz(state: qsupport.QSupportState) -> None:
    expected_e = qsupport.gamma*(cd.y+1)**9*(cd.y-qsupport.r)
    expected_G = sum(coefficient*cd.y**degree for degree, coefficient in
                     enumerate(state.g_coefficients))
    expected_sigma = (cd.y-qsupport.r)**2*expected_G
    assert sp.expand(state.ansatz.e-expected_e) == 0
    assert sp.expand(state.G-expected_G) == 0
    assert sp.expand(state.ansatz.sigma-expected_sigma) == 0
    assert sp.sympify(state.ansatz.d1) == 0
    assert tuple(state.saturation_factors) == _expected_saturation_factors(
        state)


def _independent_coefficient_check(
    state: qsupport.QSupportState,
    fast_engine: landed.HighCoefficientEngine,
    degrees: Iterable[int],
) -> None:
    '''Compare required coefficients with the full all-f convolution engine.'''
    direct = cd.ConvolutionDescent(state.ansatz, c=landed.DEFAULT_C)
    for degree in degrees:
        fast = fast_engine.master_coefficient(degree)
        slow = direct.master_coefficient(degree)
        if sp.expand(fast-slow) != 0:
            raise AssertionError(
                f'high-coefficient mismatch at master degree {degree}')


def build_state_program(z: int, coefficient_count: int) -> GeneratedState:
    if z not in (0, 1, 2):
        raise ValueError('this Phase F3 driver runs z=0,1,2 only')
    if coefficient_count < 5:
        raise ValueError('at least degrees 251..247 must be collected')
    started = time.monotonic()
    state = qsupport.build_qsupport_ansatz(z)
    _audit_ansatz(state)
    degrees = tuple(range(START_DEGREE,
                          START_DEGREE-coefficient_count, -1))
    engine = landed.HighCoefficientEngine(
        state.ansatz, start_degree=START_DEGREE,
        target_count=coefficient_count, c=landed.DEFAULT_C)
    raw_equations = tuple(engine.master_coefficient(d) for d in degrees)

    # This independent all-f engine check catches reversed-codegree,
    # sign, or indexing errors in every mandatory coefficient.
    _independent_coefficient_check(state, engine, degrees[:5])
    equations = tuple(qsupport.quotient_reduce(eq) for eq in raw_equations)
    for raw, reduced in zip(raw_equations, equations):
        assert qsupport.quotient_reduce(raw-reduced) == 0

    variables = _state_variables(state)
    variable_text = ','.join(variables)
    coefficient_lines = [
        f'poly coeff_{degree} = {_singular_expression(equation)};'
        for degree, equation in zip(degrees, equations)
    ]
    ideal_members = ','.join(f'coeff_{degree}' for degree in degrees)
    saturation_factors = tuple(state.saturation_factors)
    nonzero_product = '*'.join(
        f'({_singular_expression(factor)})' for factor in saturation_factors)
    program = '\n'.join((
        'LIB @DQ@elim.lib@DQ@;',
        'LIB @DQ@solve.lib@DQ@;',
        f'ring R = (0,a),({variable_text}),dp;',
        f'minpoly = {MINPOLY_SINGULAR};',
        *coefficient_lines,
        f'ideal I = {ideal_members};',
        f'poly nonzero_product = {nonzero_product};',
        'list SatData = sat(I,nonzero_product);',
        'ideal Isat = SatData[1];',
        'int saturation_exponent = SatData[2];',
        '@DQ@@@SATURATION_RAN@DQ@;',
        'saturation_exponent;',
        'ideal G = std(Isat);',
        'int is_unit = (reduce(1,G)==0);',
        '@DQ@@@UNIT_IDEAL@DQ@;',
        'is_unit;',
        '@DQ@@@DIMENSION@DQ@;',
        'int ideal_dimension = dim(G);',
        'ideal_dimension;',
        '@DQ@@@BASIS_BEGIN@DQ@;',
        'G;',
        '@DQ@@@BASIS_END@DQ@;',
        'if (is_unit==0 && ideal_dimension==0)',
        '{',
        '  @DQ@@@VDIM@DQ@;',
        '  vdim(G);',
        '  @DQ@@@TRIANGL_BEGIN@DQ@;',
        '  list triangular_decomposition = triangL(G);',
        '  triangular_decomposition;',
        '  @DQ@@@TRIANGL_END@DQ@;',
        '}',
        'quit;',
        '',
    )).replace('@DQ@', chr(34))
    if 'minpoly = ' + MINPOLY_SINGULAR not in program:
        raise AssertionError('number-field minpoly missing from program')
    return GeneratedState(
        z, degrees, equations, program, saturation_factors,
        time.monotonic()-started, True)


def _print_program(label: str, program: str) -> None:
    print(f'===== EXACT SINGULAR PROGRAM {label} BEGIN =====', flush=True)
    print(program, end='', flush=True)
    print(f'===== EXACT SINGULAR PROGRAM {label} END =====', flush=True)


def _print_run(run: SingularRun) -> None:
    if run.label.startswith('z=') and run.unit_ideal is True:
        print('NEW KILL ' + chr(8212) + ' PENDING AUDIT', flush=True)
    print(
        f'{run.label}: status={run.status}; returncode={run.returncode}; '
        f'saturation_ran={run.saturation_ran}; unit_ideal={run.unit_ideal}; '
        f'dim={run.dimension}; vdim={run.vector_space_dimension}; '
        f'wall={run.elapsed_seconds:.3f}s', flush=True)
    if run.basis_text is not None:
        print(f'{run.label} basis:\n{run.basis_text}', flush=True)
    if run.stdout:
        print(f'{run.label} stdout:\n{run.stdout}', flush=True)
    if run.stderr:
        print(f'{run.label} stderr:\n{run.stderr}', flush=True)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--coefficient-count', type=int, default=DEFAULT_COEFFICIENT_COUNT,
        help='coefficients from degree 251 downward (minimum 5)')
    parser.add_argument(
        '--timeout-seconds', type=float, default=DEFAULT_TIMEOUT_SECONDS,
        help='bounded wall time for each Singular subprocess')
    parser.add_argument(
        '--programs-only', action='store_true',
        help='generate/audit and print programs without invoking WSL')
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    if args.coefficient_count < 5:
        raise SystemExit('--coefficient-count must be at least 5')
    if args.timeout_seconds <= 0:
        raise SystemExit('--timeout-seconds must be positive')

    c08 = c08_program()
    _print_program('C08', c08)
    if args.programs_only:
        for z in (0, 1, 2):
            generated = build_state_program(z, args.coefficient_count)
            print(
                f'z={z}: generated degrees {generated.degrees}; '
                f'crosscheck={generated.coefficient_crosscheck_passed}; '
                f'generation={generated.generation_seconds:.3f}s', flush=True)
            _print_program(f'z={z}', generated.singular_program)
        return 0

    bridge = run_singular('C08', c08, args.timeout_seconds)
    _print_run(bridge)
    gate_ok, gate_reason = validate_c08_run(bridge)
    gate_label = 'PASS' if gate_ok else 'FAIL'
    print(f'C08 validation gate: {gate_label}; '
          f'{gate_reason}', flush=True)
    if not gate_ok:
        print('Residue runs suppressed because the bridge gate did not pass.',
              flush=True)
        return 2

    generated0 = build_state_program(0, args.coefficient_count)
    print(f'z=0 coefficient audit: degrees={generated0.degrees}; '
          f'crosscheck=PASS; generation={generated0.generation_seconds:.3f}s',
          flush=True)
    _print_program('z=0', generated0.singular_program)
    z0 = run_singular('z=0', generated0.singular_program,
                      args.timeout_seconds)
    _print_run(z0)
    if z0.status == 'timeout':
        print('z=1 and z=2 suppressed because z=0 hit its wall timeout.',
              flush=True)
        return 0
    if z0.status != 'completed' or not z0.saturation_ran:
        print('z=1 and z=2 suppressed because z=0 did not complete cleanly '
              'with a saturation marker.', flush=True)
        return 2

    for z in (1, 2):
        generated = build_state_program(z, args.coefficient_count)
        print(f'z={z} coefficient audit: degrees={generated.degrees}; '
              f'crosscheck=PASS; generation={generated.generation_seconds:.3f}s',
              flush=True)
        _print_program(f'z={z}', generated.singular_program)
        result = run_singular(
            f'z={z}', generated.singular_program, args.timeout_seconds)
        _print_run(result)
    return 0


if __name__ == '__main__':
    sys.exit(main())
