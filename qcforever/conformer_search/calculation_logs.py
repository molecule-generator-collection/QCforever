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


def read_pm6_block(log, atom_count, *, include_geometry=False):
    """Read scalar PM6 force/energy records without extra electronic calculations.

    A final SCF Done line can precede an SCF error. Preserve it for cost
    accounting; native convergence and continuation eligibility are separate.
    Missing forces stay None, never zero. Forces are mean atom-vector norms.
    Optional geometry is for post-run storage only; native execution is unchanged.
    """
    records = []
    geometry = None
    lines = log.splitlines()
    for i, line in enumerate(lines):
        if include_geometry and re.search(r'(Input|Standard|Z-Matrix) orientation:', line):
            geometry = None
            atoms, positions = [], []
            for row in lines[i + 5:i + 5 + atom_count]:
                columns = row.split()
                if len(columns) != 6 or not columns[0].isdigit():
                    break
                atoms.append(int(columns[1]))
                positions.append([parse_finite_number(value) for value in columns[3:]])
            if len(atoms) == atom_count:
                geometry = dict(atomic_numbers=atoms, coordinates_angstrom=positions)
        match = re.search(r'SCF Done:.*?E\(([^)]+)\)\s*=\s*(\S+).*?after\s+(\d+)', line)
        if match:
            try:
                energy = parse_finite_number(match[2])
            except ValueError:
                energy = None  # Still charge the failed/nonfinite SCF evaluation.
            records.append(dict(method=match[1], energy_hartree=energy,
                                scf_cycles=int(match[3]), force_mean_hartree_per_bohr=None))
            if include_geometry:
                records[-1].update(geometry=geometry, force_vectors_hartree_per_bohr=None)
        if 'Forces (Hartrees/Bohr)' in line and records:
            magnitudes = []
            vectors = []
            for row in lines[i + 3:]:
                columns = row.split()
                if len(columns) != 5 or not columns[0].isdigit():
                    break
                try:
                    vector = [parse_finite_number(v) for v in columns[2:]]
                except ValueError:
                    magnitudes = []
                    break
                magnitudes.append(math.sqrt(sum(v * v for v in vector)))
                vectors.append(vector)
            if len(magnitudes) == atom_count:
                records[-1]['force_mean_hartree_per_bohr'] = sum(magnitudes) / atom_count
                if include_geometry:
                    records[-1]['force_vectors_hartree_per_bohr'] = vectors
    native_limits = [int(value) for value in re.findall(r'Step number\s+\d+\s+out of a maximum of\s+(\d+)', log)]
    return dict(records=records,
                native_converged='Stationary point found' in log and 'Normal termination' in log,
                cycle_limit='Number of steps exceeded' in log,
                scf_failed='Convergence failure' in log or 'SCF has not converged' in log,
                step_numbers=[int(x) for x in re.findall(r'Step number\s+(\d+)', log)],
                native_maxcycles=native_limits[-1] if native_limits else None)

def parse_finite_number(value):
    number = float(value.replace('D', 'E'))
    if not math.isfinite(number):
        raise ValueError('Nonfinite energy or gradient')
    return number


def write_json(path, value):
    """Write a calculation record; reject nonfinite JSON values."""
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')
