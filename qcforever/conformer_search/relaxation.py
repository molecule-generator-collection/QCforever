"""Continuous native relaxation, isolated per candidate; no LAQA scheduling.

xTB results keep the input graph and atom order while replacing coordinates.
Physical checks and coordinate-derived stereochemistry are audited afterwards.
"""
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor
import multiprocessing

from rdkit import Chem
from qcforever.util import job_timeout
from .pipeline import save, write_sdf
from .validation import choose_best


def _bind_relaxation_worker(core_groups):
    """Give each native worker a disjoint subset of scheduler-assigned CPUs."""
    if core_groups is not None:
        os.sched_setaffinity(0, core_groups.get())


def _relax_one(task, adapter=None):
    """Process-isolated candidate; PM6 cwd/environment never cross workers."""
    payload, directory, charge, multiplicity, cores, memory, settings, method, deadline = task
    mol = Chem.Mol(payload)
    directory = Path(directory)
    row = {'candidate_id': mol.GetProp('candidate_id'), 'prepared_index': int(directory.name.split('_')[-1])}
    if hasattr(os, 'sched_getaffinity'):
        row['cpu_affinity'] = sorted(os.sched_getaffinity(0))
    started = time.monotonic()
    try:
        timeout = deadline-time.monotonic() if deadline is not None else None
        if timeout is not None and timeout <= 0:
            raise job_timeout.QCforeverTimeoutError('xTB/PM6 overall deadline exceeded')
        write_sdf(directory/'initial.sdf', [mol])
        with job_timeout.overall_timeout(timeout):
            out, energy = (adapter or (xtb_relax if method == 'xtb' else pm6_relax))(
                mol, directory, charge, multiplicity, cores, memory, settings)
        if not math.isfinite(energy):
            raise ValueError('Nonfinite final energy')
        out.SetDoubleProp('energy_hartree', energy)
        write_sdf(directory/'optimized.sdf', [out])
        row.update(state='converged', energy_hartree=energy)
        result = out.ToBinary(Chem.PropertyPickleOptions.AllProps)
    except job_timeout.QCforeverTimeoutError as exc:
        result, energy = None, None
        row.update(state='timeout', error=str(exc))
    except Exception as exc:
        result, energy = None, None
        row.update(state='failed', error=f'{type(exc).__name__}: {exc}')
    row['wall_seconds'] = time.monotonic()-started
    save(directory/'status.json', row)
    return result, energy, row


def _finite(value):
    number = float(value.replace('D', 'E'))
    if not math.isfinite(number):
        raise ValueError('Nonfinite energy or gradient')
    return number


def xtb_trace(folder):
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
                    cycles.append({'record': len(cycles)+1, 'energy_hartree': _finite(match[1]),
                                   'gradient_norm_hartree_per_bohr': _finite(match[2])})
                except ValueError:
                    continue  # A truncated last record is not invented or filled.
    steps = [int(x) for x in re.findall(r'CYCLE\s+(\d+)', text)]
    trace = {'native_converged': 'GEOMETRY OPTIMIZATION CONVERGED' in text,
             'optimizer_cycles': max(steps) if steps else None, 'trajectory_records': cycles,
             'record_axis': 'native_trajectory_record_not_assumed_optimizer_cycle'}
    save(folder/'trace.json', trace)
    return text, trace


def xtb_relax(mol, folder, charge, multiplicity, cores, memory, settings):
    executable = shutil.which(settings.get('xtb_executable', 'xtb'))
    if executable is None:
        raise FileNotFoundError('xTB executable not found')
    Chem.MolToXYZFile(mol, str(folder/'input.xyz'))
    argv = [executable, 'input.xyz', '--gfn', '2', '--chrg', str(charge),
            '--uhf', str(multiplicity-1), '--opt', settings.get('xtb_opt_level', 'normal'),
            '--cycles', str(settings.get('maximum_cycles', 1000)), '--parallel', str(cores)]
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS=str(cores), MKL_NUM_THREADS=str(cores), OMP_STACKSIZE='256M')
    save(folder/'command.json', {'argv': argv, 'charge': charge, 'multiplicity': multiplicity})
    try:
        with (folder/'xtb.out').open('w') as out:
            job_timeout.run(argv, cwd=folder, env=env, stdout=out, stderr=subprocess.STDOUT, check=True)
    finally:
        text, trace = xtb_trace(folder)
    if not trace['native_converged']:
        raise RuntimeError('xTB native optimization did not converge')
    energies = re.findall(r'TOTAL ENERGY\s+([-+0-9.EeDd]+)', text)
    if not energies:
        raise ValueError('xTB final energy missing')
    final = Chem.MolFromXYZFile(str(folder/'xtbopt.xyz'))
    if final is None or [a.GetAtomicNum() for a in final.GetAtoms()] != [a.GetAtomicNum() for a in mol.GetAtoms()]:
        raise ValueError('xTB final atom order/composition changed')
    result = Chem.Mol(mol)
    result.RemoveAllConformers()
    result.AddConformer(Chem.Conformer(final.GetConformer()), assignId=True)
    return result, _finite(energies[-1])


def pm6_relax(mol, folder, charge, multiplicity, cores, memory, settings):
    """Use the existing Gaussian adapter, retaining its files in this folder."""
    from .bridge import working_directory
    from qcforever.laqa_fafoom.pyg16 import g16Object
    from .pm6_errors import PM6CalculationError, failure_diagnostic
    binary = shutil.which('g16')
    if binary is None:
        raise FileNotFoundError('Gaussian16 executable not found for PM6')
    previous = {key: os.environ.get(key) for key in ('GAUSS_EXEDIR', 'GAUSS_SCRDIR')}
    try:
        with working_directory(folder):
            obj = g16Object(Chem.MolToMolBlock(mol), str(Path(binary).parent), str(folder),
                            cores, memory or '1GB', 'opt', charge, multiplicity, 'pm6',
                            settings.get('maximum_cycles', 1000))
            obj.generate_input()
            try:
                obj.run_g16()
            except job_timeout.QCforeverTimeoutError:
                raise
            except Exception as exc:
                path = Path('Gau_molecule.log')
                diagnostic = failure_diagnostic(path.read_text(errors='replace') if path.is_file() else '', exc)
                save(Path('failure.json'), diagnostic)
                raise PM6CalculationError(
                    f"PM6: {diagnostic['reason']}; see {folder/'failure.json'} and Gau_molecule.log") from exc
            finally:
                path = Path('Gau_molecule.log')
                log = path.read_text(errors='replace') if path.is_file() else ''
                converged = 'Stationary point found' in log and 'Normal termination' in log
                energies = [_finite(x) for x in re.findall(r'SCF Done:.*?=\s*([-+0-9.EeDd]+)', log)]
                save(Path('trace.json'), {'native_converged': converged,
                     'energy_evaluations_hartree': energies, 'record_axis': 'SCF_evaluation'})
            if not converged:
                diagnostic = failure_diagnostic(log)
                save(Path('failure.json'), diagnostic)
                raise PM6CalculationError(f"PM6: {diagnostic['reason']}; see {folder/'failure.json'}")
            out = Chem.MolFromMolBlock(obj.get_sdf_string_opt(), removeHs=False)
            if out is None:
                raise ValueError('Unreadable PM6 final structure')
            if [a.GetAtomicNum() for a in out.GetAtoms()] != [a.GetAtomicNum() for a in mol.GetAtoms()]:
                raise ValueError('PM6 final atom order/composition changed')
            for key in mol.GetPropNames():
                out.SetProp(key, mol.GetProp(key))
            return out, obj.get_energy('hartree')
    finally:
        for key, value in previous.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value


def relax_candidates(prepared, reference, config, charge, multiplicity, method, cores, memory,
                     *, adapter=None):
    if config.relaxation.get('implementation') != 'continuous':
        raise NotImplementedError('Only all-candidate continuous relaxation is currently implemented')
    if method not in ('xtb', 'pm6'):
        raise ValueError('optconf method must be xtb or pm6')
    folder = prepared.status_path.parent/'electronic'
    folder.mkdir(exist_ok=False)
    records, energies, rows = [], [], []
    per_call = config.relaxation.get('cores_per_calculation', 4)
    if cores < per_call:
        raise ValueError('Total allocated cores must be >= relaxation.cores_per_calculation')
    workers = min(len(prepared.candidates), cores // per_call)
    started = time.monotonic()
    tasks = []
    remaining = job_timeout.remaining_time()
    deadline = time.monotonic()+remaining if remaining is not None else None
    for i, mol in enumerate(prepared.candidates):
        directory = folder/f'candidate_{i:05d}'
        directory.mkdir()
        tasks.append((mol.ToBinary(Chem.PropertyPickleOptions.AllProps), str(directory),
                      charge, multiplicity, per_call, memory, config.relaxation, method,
                      deadline))
    # Test-injected local adapters are intentionally serial (often closures).
    # Real native calculations use isolated processes, not threads, for PM6.
    context = multiprocessing.get_context('spawn')
    core_groups = None
    if workers > 1 and adapter is None and hasattr(os, 'sched_getaffinity'):
        available = sorted(os.sched_getaffinity(0))
        if len(available) < workers*per_call:
            raise ValueError('Scheduler CPU affinity is smaller than requested parallel allocation')
        core_groups = context.Queue()
        for worker in range(workers):
            core_groups.put(available[worker*per_call:(worker+1)*per_call])
    pool = ProcessPoolExecutor(max_workers=workers, mp_context=context,
        initializer=_bind_relaxation_worker, initargs=(core_groups,)) if workers > 1 and adapter is None else None
    try:
        outputs = pool.map(_relax_one, tasks) if pool else (_relax_one(t, adapter) for t in tasks)
        for payload, energy, row in outputs:
            records.append(Chem.Mol(payload) if payload is not None else None)
            energies.append(energy)
            rows.append(row)
            save(folder/'progress.json', rows)
            if row['state'] == 'timeout':
                raise job_timeout.QCforeverTimeoutError(row['error'])
    finally:
        if pool:
            pool.shutdown()
    audit = choose_best(records, energies, reference, config.validation)
    counts = {'input_candidates': len(prepared.candidates), 'attempted_candidates': len(rows),
              'converged_candidates': sum(r['state'] == 'converged' for r in rows),
              'failed_candidates': sum(r['state'] == 'failed' for r in rows)}
    audit.update(implementation='continuous', backend=method, energy_unit='hartree',
                 candidate_runs=rows, **counts, wall_seconds=time.monotonic()-started,
                 sum_candidate_wall_seconds=sum(row['wall_seconds'] for row in rows),
                 allocated_cores=cores, cores_per_calculation=per_call,
                 parallel_workers=workers if adapter is None else 1)
    primary = audit['primary_best_index']
    selected = primary if primary is not None else audit['best_any_index']
    if selected is None:
        available = [i for i, mol in enumerate(records) if mol is not None]
        selected = min(available, key=lambda i: energies[i]) if available else None
    audit['selected_index'] = selected
    audit['selected_structure_warning'] = (None if primary is not None else
        'stereo_mismatch' if audit['best_any_index'] is not None else 'geometry_invalid' if selected is not None else 'no_converged_structure')
    for mol, check in zip(records, audit['candidates']):
        if mol is not None:
            mol.SetProp('geometry_check_failure', check['geometry_failure'] or '')
            mol.SetBoolProp('geometry_valid', check['geometry_failure'] is None)
            for key in ('tetrahedral_stereo_match', 'ez_stereo_match', 'joint_stereo_match'):
                if key in check:
                    mol.SetBoolProp(key, check[key])
    save(folder/'audit.json', audit)
    write_sdf(folder/'all_converged.sdf', records)
    if selected is None:
        raise RuntimeError('No converged xTB/PM6 structure; see electronic/audit.json')
    records[selected].SetProp('structure_check_warning', audit['selected_structure_warning'] or '')
    write_sdf(Path.cwd()/'optimized_structures.sdf', [records[selected]])
    return audit
