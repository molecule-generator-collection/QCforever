"""Read optimization logs and write calculation records without changing results.

This diagnostic layer does not change legacy LAQA parsing or convergence rules.
Classification uses explicit log messages; unrecognized failures stay unknown.
"""
import json
from pathlib import Path
import math
import re


class PM6CalculationError(RuntimeError):
    pass


def read_xtb_trace(folder):
    """Keep partial trajectories even if a native process failed or timed out."""
    path = folder/'xtb.out'
    text = path.read_text(errors='replace') if path.is_file() else ''
    cycles = []
    trajectory = folder/'xtbopt.log'
    if trajectory.is_file():
        for line in trajectory.read_text().splitlines():
            match = re.search(r'energy:\s*(\S+)\s+gnorm:\s*(\S+)', line)
            if match:
                try:
                    cycles.append({'record': len(cycles)+1, 'energy_hartree': parse_finite_number(match[1]),
                                   'gradient_norm_hartree_per_bohr': parse_finite_number(match[2])})
                except ValueError:
                    continue  # A truncated last record is not invented or filled.
    steps = [int(x) for x in re.findall(r'CYCLE\s+(\d+)', text)]
    trace = {'native_converged': 'GEOMETRY OPTIMIZATION CONVERGED' in text,
             'optimizer_cycles': max(steps) if steps else None, 'trajectory_records': cycles,
             'record_axis': 'native_trajectory_record_not_assumed_optimizer_cycle'}
    write_json(folder/'trace.json', trace)
    return text, trace


def read_pm6_trace(log):
    """Return Gaussian convergence and SCF records, not assumed optimizer cycles."""
    converged = 'Stationary point found' in log and 'Normal termination' in log
    energies = [parse_finite_number(x) for x in re.findall(r'SCF Done:.*?=\s*([-+0-9.EeDd]+)', log)]
    return {'native_converged': converged, 'energy_evaluations_hartree': energies,
            'record_axis': 'SCF_evaluation'}

def diagnose_pm6_failure(log, original_error=None):
    rules = (
        ('invalid_interatomic_distances', ('small interatomic distances', 'problem with the distance matrix')),
        ('scf_not_converged', ('convergence failure', 'scf has not converged')),
        ('optimization_step_limit', ('number of steps exceeded',)),
    )
    lower = log.lower()
    reason = next((name for name, messages in rules if any(m in lower for m in messages)), None)
    has_energy = bool(re.search(r'SCF Done:.*?=\s*[-+0-9.]', log))
    if reason is None:
        reason = ('missing_native_log' if not log else 'gaussian_error_termination'
                  if 'error termination' in lower else 'missing_scf_energy'
                  if not has_energy else 'output_parse_failure'
                  if original_error is not None else 'optimization_not_converged')
    return {'reason': reason, 'native_log': 'Gau_molecule.log',
            'log_tail': log.splitlines()[-20:],
            'adapter_error': (f'{type(original_error).__name__}: {original_error}'
                              if original_error is not None else None)}


def record_xtb_call(path, jobtype, requested_cycles, wall_seconds, error=None):
    rows = []
    log = Path('xtbopt.log')
    if jobtype == 'opt' and log.is_file():
        for line in log.read_text().splitlines():
            match = re.match(r'\s*energy:\s*(\S+)\s+gnorm:\s*(\S+)', line)
            if match:
                rows.append({'cycle': len(rows)+1, 'energy_hartree': float(match[1]),
                             'gradient_norm_hartree_per_bohr': float(match[2])})
    out = Path('result.out').read_text(errors='replace') if Path('result.out').is_file() else ''
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(f'Refusing to overwrite xTB trace: {path}')
    write_json(path, {'jobtype': jobtype, 'requested_cycles': requested_cycles,
                'actual_cycles': len(rows), 'wall_seconds': wall_seconds,
                'native_converged': 'GEOMETRY OPTIMIZATION CONVERGED' in out,
                'error': error, 'cycles': rows})


def parse_finite_number(value):
    number = float(value.replace('D', 'E'))
    if not math.isfinite(number):
        raise ValueError('Nonfinite energy or gradient')
    return number


def write_json(path, value):
    """Write a calculation record; reject nonfinite JSON values."""
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')
