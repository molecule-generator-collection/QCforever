"""Keep readable relaxation histories, then remove finished native work files.

Call only after all native workers have stopped. Each candidate retains
trajectory.jsonl (units, step/energy/force records), trajectory.xyz (structures),
and status.json. Failed calculations also retain failure.log with native log
tails. With keep_logs=True, native.log instead contains complete text logs from
all calls. After successful export, archive_relaxation_details bundles histories
and execution records, leaving structures, audit and settings directly readable.
A failed export preserves the candidate's original work files.
"""
import json
from pathlib import Path
import re
import shutil
import tarfile

from rdkit import Chem

from .calculation_logs import parse_finite_number, read_pm6_block, write_json
from .structure_file_io import write_text_atomically


def compact_relaxation_outputs(electronic_directory, backend, *, keep_logs=False):
    """Compact this run's candidate directories without changing chemistry."""
    root = Path(electronic_directory)
    if root.is_symlink() or not root.is_dir():
        return
    if root.name != 'electronic' or backend not in ('pm6', 'xtb'):
        raise ValueError('Expected an electronic directory and pm6/xtb backend')
    report = dict(backend=backend, keep_logs=keep_logs, candidates=[], errors=[])
    # A live xTB process must never lose its output or restart state.
    if (root / 'cleanup_errors.json').exists():
        print('WARNING: Native process cleanup failed; relaxation files were preserved')
        return
    for folder in sorted(root.iterdir()):
        if not re.fullmatch(r'candidate_\d+', folder.name) or folder.is_symlink() or not folder.is_dir():
            continue
        try:
            result = _compact_candidate(folder, backend, keep_logs)
            report['candidates'].append(dict(candidate=folder.name, **result))
        except Exception as exc:
            report['errors'].append(dict(candidate=folder.name, error=f'{type(exc).__name__}: {exc}'))
    write_json(root / 'history_cleanup.json', report)
    if report['errors']:
        print(f'WARNING: Some relaxation files were preserved; see {root}/history_cleanup.json')
    return report


def archive_relaxation_details(electronic_directory):
    """Archive closed work records, then remove originals after successful writing.

    Call only after successful history export and after native workers stop.
    Relative paths inside details.tar.gz are rooted at electronic/. Existing
    archives are not overwritten; compression failure leaves sources intact.
    """
    root = Path(electronic_directory)
    if root.is_symlink() or not root.is_dir() or root.name != 'electronic':
        raise ValueError('Expected an electronic directory')
    archive = root / 'details.tar.gz'
    temporary = root / 'details.tar.gz.tmp'
    if archive.exists() or archive.is_symlink() or temporary.exists() or temporary.is_symlink():
        raise FileExistsError('Relaxation archive already exists')
    retained = {'all_converged.sdf', 'audit.json', 'settings.json'}
    entries = sorted(path for path in root.iterdir() if path.name not in retained)
    if not entries:
        return
    try:
        with tarfile.open(temporary, 'x:gz', dereference=False) as output:
            for entry in entries:
                output.add(entry, arcname=entry.name)
        temporary.rename(archive)
    finally:
        temporary.unlink(missing_ok=True)
    for entry in entries:
        if entry.is_dir() and not entry.is_symlink():
            shutil.rmtree(entry)
        else:
            entry.unlink()


def _compact_candidate(folder, backend, keep_logs):
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
        if rows and rows[0]['record_type'] == 'initial_structure':
            coordinate_precision = 'native input.xyz coordinates; optimization never started'
    metadata = dict(record_type='metadata', schema_version=1, backend=backend,
        candidate_id=molecule.GetProp('candidate_id') if molecule.HasProp('candidate_id') else folder.name,
        atomic_numbers=atoms, step_axis=step_axis, energy_unit='hartree',
        force_unit='hartree/bohr', coordinate_unit='angstrom',
        force_definition=force_definition, coordinate_precision=coordinate_precision,
        optimizer_stdout_definition='Original xTB cycle/ANC-gradient/evaluation records, separate from Cartesian XYZ frames',
        structure_file=structures.name, omitted_incomplete_xyz_tail=truncated,
        records=len(rows), structures=len(frames))
    if rows and rows[0]['record_type'] == 'initial_structure':
        metadata['failure_reason'] = rows[0]['failure_reason']
        metadata['optimization_started'] = False
    data = [metadata, *rows, *[dict(row, record_type='optimizer_stdout') for row in optimizer_rows]]
    json_text = ''.join(json.dumps(row, allow_nan=False) + '\n' for row in data)
    symbols = [Chem.GetPeriodicTable().GetElementSymbol(number) for number in atoms]
    xyz_lines = []
    for row, coordinates in zip(rows, frames):
        xyz_lines.append(f'{len(atoms)}\nstep={row["step"]} energy_hartree={row["energy_hartree"]}\n')
        for symbol, xyz in zip(symbols, coordinates):
            xyz_lines.append(symbol + ' ' + ' '.join(f'{value:.14f}' for value in xyz) + '\n')
    write_text_atomically(structures, ''.join(xyz_lines))
    write_text_atomically(history, json_text)
    _save_native_logs(folder, backend, keep_logs)
    for path in list(folder.iterdir()):
        if path.name in ('trajectory.jsonl', 'trajectory.xyz', 'status.json', 'failure.log', 'native.log'):
            continue
        if path.is_dir() and not path.is_symlink():
            shutil.rmtree(path)
        else:
            path.unlink()
    return dict(records=len(rows))


def _save_native_logs(folder, backend, keep_logs):
    """Save complete logs on request, otherwise only bounded failure diagnostics."""
    status_path = folder / 'status.json'
    status = json.loads(status_path.read_text()) if status_path.exists() else {}
    if not keep_logs and status.get('state') != 'failed':
        return
    logs = [folder / 'xtb.out']
    if backend == 'pm6':
        calls = sorted(path for path in folder.iterdir()
                       if re.fullmatch(r'call_\d+', path.name)
                       and path.is_dir() and not path.is_symlink())
        logs = [call / name for call in (calls if keep_logs else calls[-1:])
                for name in ('input.log', 'launcher.out', 'formchk.out')]
    output = folder / ('native.log' if keep_logs else 'failure.log')
    # Stream full logs to avoid holding large Gaussian outputs in memory.
    # Any write error propagates before native work files are removed.
    with output.open('wb') as saved:
        message = status.get('error') or status.get('state', 'Calculation failed')
        saved.write((str(message) + '\n').encode())
        for path in logs:
            if not path.is_file() or path.is_symlink():
                continue
            suffix = '' if keep_logs else ' (last 64 KiB)'
            saved.write(f'\n--- {path.relative_to(folder)}{suffix} ---\n'.encode())
            with path.open('rb') as stream:
                if not keep_logs:
                    stream.seek(max(0, path.stat().st_size - 65536))
                shutil.copyfileobj(stream, saved)
            saved.write(b'\n')


def _read_pm6_history(folder, atoms):
    rows, frames = [], []
    calls = sorted(path for path in folder.iterdir()
                   if re.fullmatch(r'call_\d+', path.name) and path.is_dir() and not path.is_symlink())
    for call in calls:
        log = call / 'input.log'
        if not log.is_file():
            raise ValueError(f'Missing native log: {log}')
        parsed = read_pm6_block(log.read_text(errors='replace'), len(atoms), include_geometry=True)
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
            return _read_xtb_initial_scc_failure(folder, atoms)
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


def _read_xtb_initial_scc_failure(folder, atoms):
    """Keep the input for a known pre-optimization failure, never invent a step.

    Missing trajectories after optimization started, unknown native failures,
    or inconsistent status still raise and preserve the original work files.
    """
    text = (folder / 'xtb.out').read_text(errors='replace')
    status = json.loads((folder / 'status.json').read_text())
    history = folder / 'trajectory.jsonl'
    if not (status.get('state') == 'failed'
            and (folder / '.sccnotconverged').is_file()
            and 'scf: Self consistent charge iterator did not converge' in text
            and not re.search(r'\bCYCLE\s+\d+', text)
            and not (history.exists() and history.stat().st_size)):
        raise ValueError('xTB output exists without a readable trajectory')
    initial = Chem.MolFromXYZFile(str(folder / 'input.xyz'))
    if initial is None or [atom.GetAtomicNum() for atom in initial.GetAtoms()] != atoms:
        raise ValueError('Missing or inconsistent xTB input geometry')
    row = dict(record_type='initial_structure', step=0, energy_hartree=None,
        gradient_norm_hartree_per_bohr=None,
        failure_reason='scc_not_converged_before_optimization')
    return [row], [initial.GetConformer().GetPositions().tolist()], False
