"""Coordinate-free records of native xTB calls, captured before cleanup."""
from pathlib import Path
import re

from .pipeline import save


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
    save(path, {'jobtype': jobtype, 'requested_cycles': requested_cycles,
                'actual_cycles': len(rows), 'wall_seconds': wall_seconds,
                'native_converged': 'GEOMETRY OPTIMIZATION CONVERGED' in out,
                'error': error, 'cycles': rows})
