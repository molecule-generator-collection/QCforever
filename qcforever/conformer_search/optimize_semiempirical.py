"""Optimize candidates with xTB or PM6, then check and save the results.

Each candidate uses one core. Separate processes isolate Gaussian's working
directory and environment; parallel completion never changes candidate order.
"""
from contextlib import contextmanager
from dataclasses import dataclass
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
from .structure_file_io import write_sdf
from .calculation_logs import (
    write_json, parse_finite_number, read_xtb_trace, read_pm6_trace,
    PM6CalculationError, diagnose_pm6_failure,
)
from .check_structures import (
    evaluate_optimized_candidates, select_optimized_candidate,
    copy_with_coordinates, clear_structure_audit, STEREO_MATCH_PROPERTIES,
)


@dataclass(frozen=True)
class CandidateTask:
    """Serializable inputs for one candidate; no state is shared between processes."""
    molecule: bytes
    directory: str
    charge: int
    multiplicity: int
    cores: int
    memory: str
    settings: dict
    method: str
    deadline: float | None


def optimize_candidates(prepared, reference, config, charge, multiplicity, method, cores, memory,
                        *, adapter=None):
    """Optimize all candidates, audit final structures, and save the selected one."""
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
    # One core per native calculation; use nproc to parallelize candidates.
    # Model-generation threads are configured independently.
    cores_per_candidate = 1
    workers = min(cores, len(prepared.candidates))
    started = time.monotonic()
    tasks = []
    remaining = job_timeout.remaining_time()
    deadline = time.monotonic()+remaining if remaining is not None else None
    for i, mol in enumerate(prepared.candidates):
        directory = folder/f'candidate_{i:05d}'
        directory.mkdir()
        tasks.append(CandidateTask(
            molecule=mol.ToBinary(Chem.PropertyPickleOptions.AllProps),
            directory=str(directory), charge=charge, multiplicity=multiplicity,
            cores=cores_per_candidate, memory=memory, settings=config.relaxation,
            method=method, deadline=deadline))
    completed = _execute_tasks(tasks, workers, cores_per_candidate, folder, adapter)
    # Completion order must not change candidate IDs or best-index semantics.
    for i in sorted(completed):
        payload, energy, candidate_result = completed[i]
        records.append(Chem.Mol(payload) if payload is not None else None)
        energies.append(energy)
        rows.append(candidate_result)
    audit = evaluate_optimized_candidates(records, energies, reference, config.validation)
    counts = {'input_candidates': len(prepared.candidates), 'attempted_candidates': len(rows),
              'converged_candidates': sum(r['state'] == 'converged' for r in rows),
              'failed_candidates': sum(r['state'] == 'failed' for r in rows)}
    audit.update(implementation='continuous', backend=method, energy_unit='hartree',
                 candidate_runs=rows, **counts, wall_seconds=time.monotonic()-started,
                 sum_candidate_wall_seconds=sum(candidate_result['wall_seconds'] for candidate_result in rows),
                 allocated_cores=cores, cores_per_calculation=cores_per_candidate,
                 parallel_workers=workers, scheduling='parallel_candidates_1_core')
    select_optimized_candidate(records, energies, audit)
    _save_optimization_results(folder, records, audit)
    return audit


def _execute_tasks(tasks, workers, cores_per_candidate, folder, adapter):
    """Run independent candidates and write progress in original candidate order."""
    completed = {}

    def collect(result):
        payload, energy, candidate_result = result
        completed[candidate_result['prepared_index']] = (payload, energy, candidate_result)
        write_json(folder/'progress.json', [completed[i][2] for i in sorted(completed)])
        if candidate_result['state'] == 'timeout':
            raise job_timeout.QCforeverTimeoutError(candidate_result['error'])

    if workers <= 1:
        for task in tasks:
            collect(_optimize_one_candidate(task, adapter))
    else:
        # Spawn avoids inheriting RDKit/BLAS threads and Gaussian global state.
        context = multiprocessing.get_context('spawn')
        cpu_sets = None
        if hasattr(os, 'sched_getaffinity') and hasattr(os, 'sched_setaffinity'):
            available = sorted(os.sched_getaffinity(0))
            if len(available) < workers * cores_per_candidate:
                raise ValueError('Native workers exceed scheduler CPU affinity allocation')
            cpu_sets = context.Queue()
            for i in range(workers):
                cpu_sets.put(available[i * cores_per_candidate:(i + 1) * cores_per_candidate])
        with ProcessPoolExecutor(max_workers=workers,
                                 mp_context=context, initializer=_initialize_native_worker,
                                 initargs=(cpu_sets,)) as pool:
            futures = [pool.submit(_optimize_one_candidate, task, adapter) for task in tasks]
            try:
                for future in as_completed(futures):
                    collect(future.result())
            except BaseException:
                for future in futures:
                    future.cancel()
                raise
    return completed


def _optimize_one_candidate(task, adapter=None):
    """One candidate in its own directory and, when parallel, its own process."""
    mol = Chem.Mol(task.molecule)
    directory = Path(task.directory)
    candidate_result = {'candidate_id': mol.GetProp('candidate_id'), 'prepared_index': int(directory.name.split('_')[-1])}
    if hasattr(os, 'sched_getaffinity'):
        candidate_result['cpu_affinity'] = sorted(os.sched_getaffinity(0))
    started = time.monotonic()
    try:
        timeout = task.deadline-time.monotonic() if task.deadline is not None else None
        if timeout is not None and timeout <= 0:
            raise job_timeout.QCforeverTimeoutError('xTB/PM6 overall deadline exceeded')
        write_sdf(directory/'initial.sdf', [mol])
        with job_timeout.overall_timeout(timeout):
            optimize = adapter or (optimize_with_xtb if task.method == 'xtb' else optimize_with_pm6)
            out, energy = optimize(mol, directory, task.charge, task.multiplicity,
                                   task.cores, task.memory, task.settings)
        if not math.isfinite(energy):
            raise ValueError('Nonfinite final energy')
        clear_structure_audit(out)
        out.SetProp('stereo_check_stage', 'after_relaxation')
        out.SetProp('stereo_check_status', 'pending')
        out.SetDoubleProp('energy_hartree', energy)
        write_sdf(directory/'optimized.sdf', [out])
        candidate_result.update(state='converged', energy_hartree=energy)
        result = out.ToBinary(Chem.PropertyPickleOptions.AllProps)
    except job_timeout.QCforeverTimeoutError as exc:
        result, energy = None, None
        candidate_result.update(state='timeout', error=str(exc))
    except Exception as exc:
        result, energy = None, None
        candidate_result.update(state='failed', error=f'{type(exc).__name__}: {exc}')
    candidate_result['wall_seconds'] = time.monotonic()-started
    write_json(directory/'status.json', candidate_result)
    return result, energy, candidate_result


def optimize_with_xtb(mol, folder, charge, multiplicity, cores, memory, settings):
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
    write_json(folder/'command.json', {'argv': argv, 'charge': charge, 'multiplicity': multiplicity,
                                'thread_environment': thread_env})
    try:
        with (folder/'xtb.out').open('w') as out:
            job_timeout.run(argv, cwd=folder, env=env, stdout=out, stderr=subprocess.STDOUT, check=True)
    finally:
        text, trace = read_xtb_trace(folder)
    if not trace['native_converged']:
        raise RuntimeError('xTB native optimization did not converge')
    energies = re.findall(r'TOTAL ENERGY\s+([-+0-9.EeDd]+)', text)
    if not energies:
        raise ValueError('xTB final energy missing')
    final = Chem.MolFromXYZFile(str(folder/'xtbopt.xyz'))
    if final is None or [a.GetAtomicNum() for a in final.GetAtoms()] != [a.GetAtomicNum() for a in mol.GetAtoms()]:
        raise ValueError('xTB final atom order/composition changed')
    return copy_with_coordinates(mol, final), parse_finite_number(energies[-1])


def optimize_with_pm6(mol, folder, charge, multiplicity, cores, memory, settings):
    """Use the existing Gaussian adapter, retaining its files in this folder."""
    from qcforever.laqa_fafoom.pyg16 import g16Object
    binary = shutil.which('g16')
    if binary is None:
        raise FileNotFoundError('Gaussian16 executable not found for PM6')
    thread_env = native_thread_environment(cores)
    previous = {key: os.environ.get(key) for key in (*thread_env, 'GAUSS_EXEDIR', 'GAUSS_SCRDIR')}
    try:
        os.environ.update(thread_env)
        write_json(folder/'environment.json', thread_env)
        with working_directory(folder):
            gaussian_job = g16Object(Chem.MolToMolBlock(mol), str(Path(binary).parent), str(folder),
                            cores, memory or '1GB', 'opt', charge, multiplicity, 'pm6',
                            settings.get('maximum_cycles', 1000))
            gaussian_job.generate_input()
            try:
                gaussian_job.run_g16()
            except job_timeout.QCforeverTimeoutError:
                raise
            except Exception as exc:
                path = Path('Gau_molecule.log')
                diagnostic = diagnose_pm6_failure(path.read_text(errors='replace') if path.is_file() else '', exc)
                write_json(Path('failure.json'), diagnostic)
                raise PM6CalculationError(
                    f"PM6: {diagnostic['reason']}; see {folder/'failure.json'} and Gau_molecule.log") from exc
            finally:
                path = Path('Gau_molecule.log')
                log = path.read_text(errors='replace') if path.is_file() else ''
                trace = read_pm6_trace(log)
                converged = trace['native_converged']
                write_json(Path('trace.json'), trace)
            if not converged:
                diagnostic = diagnose_pm6_failure(log)
                write_json(Path('failure.json'), diagnostic)
                raise PM6CalculationError(f"PM6: {diagnostic['reason']}; see {folder/'failure.json'}")
            # The Gaussian adapter edits coordinates in an SDF template. Read
            # positions without chemical sanitization, then retain the input
            # graph/electronic annotations exactly, as for the xTB backend.
            out = Chem.MolFromMolBlock(gaussian_job.get_sdf_string_opt(), removeHs=False, sanitize=False)
            if out is None:
                raise ValueError('Unreadable PM6 final structure')
            return copy_with_coordinates(mol, out), gaussian_job.get_energy('hartree')
    finally:
        for key, value in previous.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value


def _save_optimization_results(folder, records, audit):
    """Save every converged structure and return the selected structure through SDF."""
    selected = audit['selected_index']
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
    write_json(folder/'audit.json', audit)
    write_sdf(folder/'all_converged.sdf', records)
    if selected is None:
        raise RuntimeError('No converged xTB/PM6 structure; see electronic/audit.json')
    records[selected].SetProp('structure_check_warning', audit['selected_structure_warning'] or '')
    write_sdf(Path.cwd()/'optimized_structures.sdf', [records[selected]])


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


@contextmanager
def working_directory(path):
    previous = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)
