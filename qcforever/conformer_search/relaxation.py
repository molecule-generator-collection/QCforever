"""Continuous native relaxation, isolated per candidate; no LAQA scheduling.

xTB results keep the input graph and atom order while replacing coordinates.
Physical checks and coordinate-derived stereochemistry are audited afterwards.
"""
import math
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
from pathlib import Path
import re
import shutil
import subprocess
import time

from rdkit import Chem
from qcforever.util import job_timeout
from .pipeline import save, write_sdf
from .validation import choose_best, copy_with_coordinates, clear_structure_audit, STEREO_MATCH_PROPERTIES


def native_thread_environment(cores):
    """Native relaxation uses the per-candidate allocation, not model threads.

    Set only on the native subprocess (xTB), or temporarily around the legacy
    Gaussian adapter (PM6), restoring the caller environment afterwards.
    """
    return {key: str(cores) for key in (
        'OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'MKL_NUM_THREADS',
        'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS')}


def _initialize_native_worker(cpu_sets):
    """Disjoint Linux affinity also keeps Gaussian's scheduler wrapper bounded."""
    if cpu_sets is not None:
        os.sched_setaffinity(0, cpu_sets.get())


def _relax_one(task, adapter=None):
    """One candidate in its own directory and, when parallel, its own process."""
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
        clear_structure_audit(out)
        out.SetProp('stereo_check_stage', 'after_relaxation')
        out.SetProp('stereo_check_status', 'pending')
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
    thread_env = native_thread_environment(cores)
    env.update(thread_env, OMP_STACKSIZE='256M')
    save(folder/'command.json', {'argv': argv, 'charge': charge, 'multiplicity': multiplicity,
                                'thread_environment': thread_env})
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
    return copy_with_coordinates(mol, final), _finite(energies[-1])


def pm6_relax(mol, folder, charge, multiplicity, cores, memory, settings):
    """Use the existing Gaussian adapter, retaining its files in this folder."""
    from .conformer_search import working_directory
    from qcforever.laqa_fafoom.pyg16 import g16Object
    from .pm6_errors import PM6CalculationError, failure_diagnostic
    binary = shutil.which('g16')
    if binary is None:
        raise FileNotFoundError('Gaussian16 executable not found for PM6')
    thread_env = native_thread_environment(cores)
    previous = {key: os.environ.get(key) for key in (*thread_env, 'GAUSS_EXEDIR', 'GAUSS_SCRDIR')}
    try:
        os.environ.update(thread_env)
        save(folder/'environment.json', thread_env)
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
            # The Gaussian adapter edits coordinates in an SDF template. Read
            # positions without chemical sanitization, then retain the input
            # graph/electronic annotations exactly, as for the xTB backend.
            out = Chem.MolFromMolBlock(obj.get_sdf_string_opt(), removeHs=False, sanitize=False)
            if out is None:
                raise ValueError('Unreadable PM6 final structure')
            return copy_with_coordinates(mol, out), obj.get_energy('hartree')
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
    if isinstance(cores, bool) or not isinstance(cores, int) or cores < 1:
        raise ValueError('nproc must be a positive integer')
    folder = prepared.status_path.parent/'electronic'
    folder.mkdir(exist_ok=False)
    records, energies, rows = [], [], []
    # Separate processes are essential: the Gaussian adapter changes cwd/env.
    # Never exceed nproc, including allocations smaller than four cores.
    per_call = min(4, cores)
    workers = min(max(1, cores // per_call), len(prepared.candidates))
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
    completed = {}

    def collect(result):
        payload, energy, row = result
        completed[row['prepared_index']] = (payload, energy, row)
        save(folder/'progress.json', [completed[i][2] for i in sorted(completed)])
        if row['state'] == 'timeout':
            raise job_timeout.QCforeverTimeoutError(row['error'])

    if workers <= 1:
        for task in tasks:
            collect(_relax_one(task, adapter))
    else:
        # Spawn avoids inheriting RDKit/BLAS threads and Gaussian global state.
        context = multiprocessing.get_context('spawn')
        cpu_sets = None
        if hasattr(os, 'sched_getaffinity') and hasattr(os, 'sched_setaffinity'):
            available = sorted(os.sched_getaffinity(0))
            if len(available) < workers * per_call:
                raise ValueError('Native workers exceed scheduler CPU affinity allocation')
            cpu_sets = context.Queue()
            for i in range(workers):
                cpu_sets.put(available[i * per_call:(i + 1) * per_call])
        with ProcessPoolExecutor(max_workers=workers,
                                 mp_context=context, initializer=_initialize_native_worker,
                                 initargs=(cpu_sets,)) as pool:
            futures = [pool.submit(_relax_one, task, adapter) for task in tasks]
            try:
                for future in as_completed(futures):
                    collect(future.result())
            except BaseException:
                for future in futures:
                    future.cancel()
                raise
    # Completion order must not change candidate IDs or best-index semantics.
    for i in sorted(completed):
        payload, energy, row = completed[i]
        records.append(Chem.Mol(payload) if payload is not None else None)
        energies.append(energy)
        rows.append(row)
    audit = choose_best(records, energies, reference, config.validation)
    counts = {'input_candidates': len(prepared.candidates), 'attempted_candidates': len(rows),
              'converged_candidates': sum(r['state'] == 'converged' for r in rows),
              'failed_candidates': sum(r['state'] == 'failed' for r in rows)}
    audit.update(implementation='continuous', backend=method, energy_unit='hartree',
                 candidate_runs=rows, **counts, wall_seconds=time.monotonic()-started,
                 sum_candidate_wall_seconds=sum(row['wall_seconds'] for row in rows),
                 allocated_cores=cores, cores_per_calculation=per_call,
                 parallel_workers=workers, scheduling='parallel_candidates_4_cores')
    primary = audit['primary_best_index']
    selected = primary if primary is not None else audit['best_any_index']
    if selected is None:
        available = [i for i, mol in enumerate(records) if mol is not None]
        selected = min(available, key=lambda i: energies[i]) if available else None
    audit['selected_index'] = selected
    audit['selected_structure_warning'] = (None if primary is not None else
        'stereo_mismatch' if audit['best_any_index'] is not None else 'geometry_invalid' if selected is not None else 'no_converged_structure')
    for i, (mol, check) in enumerate(zip(records, audit['candidates'])):
        if mol is not None:
            clear_structure_audit(mol)
            mol.SetProp('geometry_check_failure', check['geometry_failure'] or '')
            mol.SetBoolProp('geometry_valid', check['geometry_failure'] is None)
            mol.SetProp('stereo_check_stage', 'after_relaxation')
            mol.SetProp('stereo_check_status', check['stereo_check_status'])
            for key in STEREO_MATCH_PROPERTIES:
                if key in check:
                    mol.SetBoolProp(key, check[key])
            write_sdf(folder/f'candidate_{i:05d}'/'optimized.sdf', [mol])
    save(folder/'audit.json', audit)
    write_sdf(folder/'all_converged.sdf', records)
    if selected is None:
        raise RuntimeError('No converged xTB/PM6 structure; see electronic/audit.json')
    records[selected].SetProp('structure_check_warning', audit['selected_structure_warning'] or '')
    write_sdf(Path.cwd()/'optimized_structures.sdf', [records[selected]])
    return audit
