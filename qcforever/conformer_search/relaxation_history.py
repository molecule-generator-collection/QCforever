"""Keep readable relaxation histories, then remove finished native work files.

Call only after all native workers have stopped. Each candidate retains
trajectory.jsonl (units, step/energy/force records), trajectory.xyz (structures),
and status.json. Search-level settings, selection, cost and final SDF files are
untouched. A failed export preserves that candidate's original work files.
"""
import json
from pathlib import Path
import re
import shutil

from rdkit import Chem

from .calculation_logs import parse_finite_number, read_pm6_block, write_json
from .structure_file_io import write_verified_text


def compact_relaxation_outputs(electronic_directory, backend):
    """Compact this run's candidate directories without changing chemistry."""
    root = Path(electronic_directory)
    if root.is_symlink() or not root.is_dir():
        return
    if root.name != 'electronic' or backend not in ('pm6', 'xtb'):
        raise ValueError('Expected an electronic directory and pm6/xtb backend')
    report = dict(backend=backend, candidates=[], errors=[])
    # A live xTB process must never lose its output or restart state.
    if (root / 'cleanup_errors.json').exists():
        print('WARNING: Native process cleanup failed; relaxation files were preserved')
        return
    for folder in sorted(root.iterdir()):
        if not re.fullmatch(r'candidate_\d+', folder.name) or folder.is_symlink() or not folder.is_dir():
            continue
        try:
            result = _compact_candidate(folder, backend)
            report['candidates'].append(dict(candidate=folder.name, **result))
        except Exception as exc:
            report['errors'].append(dict(candidate=folder.name, error=f'{type(exc).__name__}: {exc}'))
    write_json(root / 'history_cleanup.json', report)
    if report['errors']:
        print(f'WARNING: Some relaxation files were preserved; see {root}/history_cleanup.json')


def _compact_candidate(folder, backend):
    history = folder / 'trajectory.jsonl'
    structures = folder / 'trajectory.xyz'
    # Do not run the conversion twice or silently replace an unrelated export.
    if structures.exists():
        raise FileExistsError(f'History already exists: {structures}')
    molecule = next((mol for mol in Chem.SDMolSupplier(
        str(folder / 'initial.sdf'), removeHs=False, sanitize=False) if mol is not None), None)
    if molecule is None:
        raise ValueError('Missing readable initial structure')
    atoms = [atom.GetAtomicNum() for atom in molecule.GetAtoms()]
    if backend == 'pm6':
        rows, frames, truncated = _read_pm6_history(folder, atoms)
        optimizer_rows = []
        step_axis = 'SCF evaluation order across restart calls'
        force_definition = 'mean of Cartesian atom force-vector norms'
        coordinate_precision = 'Gaussian printed orientation, normally 0.000001 angstrom'
    else:
        rows, frames, truncated = _read_xtb_history(folder, atoms)
        # Stdout ANC gradients and XYZ Cartesian gradients are distinct data.
        # Retain both without assuming that their row numbers are identical.
        optimizer_rows = [json.loads(line) for line in history.read_text().splitlines()] if history.exists() else []
        step_axis = 'saved native XYZ frame; not assumed to equal optimizer cycle'
        force_definition = 'Cartesian gradient norm from native XYZ header'
        coordinate_precision = 'native xtbopt.log coordinates'
    metadata = dict(record_type='metadata', schema_version=1, backend=backend,
        candidate_id=molecule.GetProp('candidate_id') if molecule.HasProp('candidate_id') else folder.name,
        atomic_numbers=atoms, step_axis=step_axis, energy_unit='hartree',
        force_unit='hartree/bohr', coordinate_unit='angstrom',
        force_definition=force_definition, coordinate_precision=coordinate_precision,
        optimizer_stdout_definition='Original xTB cycle/ANC-gradient/evaluation records, separate from Cartesian XYZ frames',
        structure_file=structures.name, omitted_incomplete_xyz_tail=truncated,
        records=len(rows), structures=len(frames))
    data = [metadata, *rows, *[dict(row, record_type='optimizer_stdout') for row in optimizer_rows]]
    json_text = ''.join(json.dumps(row, allow_nan=False) + '\n' for row in data)
    symbols = [Chem.GetPeriodicTable().GetElementSymbol(number) for number in atoms]
    xyz_lines = []
    for row, coordinates in zip(rows, frames):
        xyz_lines.append(f'{len(atoms)}\nstep={row["step"]} energy_hartree={row["energy_hartree"]}\n')
        for symbol, xyz in zip(symbols, coordinates):
            xyz_lines.append(symbol + ' ' + ' '.join(f'{value:.14f}' for value in xyz) + '\n')
    write_verified_text(structures, ''.join(xyz_lines))
    write_verified_text(history, json_text)
    for path in list(folder.iterdir()):
        if path.name in ('trajectory.jsonl', 'trajectory.xyz', 'status.json'):
            continue
        if path.is_dir() and not path.is_symlink():
            shutil.rmtree(path)
        else:
            path.unlink()
    return dict(records=len(rows))


def _read_pm6_history(folder, atoms):
    rows, frames = [], []
    calls = sorted(path for path in folder.iterdir()
                   if re.fullmatch(r'call_\d+', path.name) and path.is_dir() and not path.is_symlink())
    for call in calls:
        log = call / 'input.log'
        if not log.is_file():
            raise ValueError(f'Missing native log: {log}')
        parsed = read_pm6_block(log.read_text(errors='replace'), len(atoms), include_geometry=True)
        trace = call / 'trace.json'
        if trace.exists() and len(json.loads(trace.read_text())['records']) != len(parsed['records']):
            raise ValueError(f'SCF record count changed: {log}')
        for record in parsed['records']:
            geometry = record['geometry']
            if geometry is None or geometry['atomic_numbers'] != atoms:
                raise ValueError(f'Missing geometry or changed atom order: {log}')
            frames.append(geometry['coordinates_angstrom'])
            rows.append(dict(record_type='step', step=len(rows) + 1, call=int(call.name[5:]),
                energy_hartree=record['energy_hartree'], scf_cycles=record['scf_cycles'],
                force_mean_hartree_per_bohr=record['force_mean_hartree_per_bohr'],
                force_vectors_hartree_per_bohr=record['force_vectors_hartree_per_bohr']))
    return rows, frames, False


def _read_xtb_history(folder, atoms):
    rows, frames = [], []
    path = folder / 'xtbopt.log'
    if not path.exists():
        # No optimization data should be discarded merely because XYZ is absent.
        if (folder / 'xtb.out').exists():
            raise ValueError('xTB output exists without a readable trajectory')
        return rows, frames, False
    with path.open() as stream:
        while True:
            line = stream.readline()
            if not line:
                return rows, frames, False
            if not line.endswith('\n'):
                return rows, frames, True
            if int(line.strip()) != len(atoms):
                raise ValueError('xTB trajectory atom count changed')
            header = stream.readline()
            coordinates = [stream.readline() for _ in atoms]
            if not header.endswith('\n') or any(not row.endswith('\n') for row in coordinates):
                return rows, frames, True
            match = re.search(r'energy:\s*(\S+)\s+gnorm:\s*(\S+)', header)
            if match is None:
                raise ValueError('Unrecognized xTB trajectory header')
            positions, numbers = [], []
            for line in coordinates:
                symbol, x, y, z = line.split()
                numbers.append(Chem.GetPeriodicTable().GetAtomicNumber(symbol))
                positions.append([parse_finite_number(value) for value in (x, y, z)])
            if numbers != atoms:
                raise ValueError('xTB trajectory atom order changed')
            frames.append(positions)
            rows.append(dict(record_type='step', step=len(rows) + 1,
                energy_hartree=parse_finite_number(match[1]),
                gradient_norm_hartree_per_bohr=parse_finite_number(match[2])))
