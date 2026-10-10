"""Optimize conformers with xTB or PM6, then check and save the results.

Each candidate uses one core. Separate processes isolate Gaussian's working
directory and environment; parallel completion never changes candidate order.
"""
from dataclasses import asdict, dataclass, replace
import json
import math
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
from pathlib import Path
import time

from rdkit import Chem
from qcforever.util import job_timeout
from qcforever.util.check_resource import native_thread_environment
from .run_xtb import optimize_with_xtb, XTBCandidate, XTB_FORCE_KIND, pool_memory_mb, cleanup_on_signals
from .run_pm6 import optimize_with_pm6, PM6Candidate
from .settings import resolve_relaxation
from .graybox_relaxation import GrayboxSearch
from .run_candidate_blocks import CandidateBlockExecutor, resolve_block_workers
from .structure_file_io import write_sdf
from .calculation_logs import write_json
from .check_structures import (
    evaluate_optimized_candidates, select_optimized_candidate,
    clear_structure_audit, STEREO_MATCH_PROPERTIES,
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
    """Optimize the requested pool continuously or adaptively, then audit/select."""
    if method not in ('xtb', 'pm6'):
        raise ValueError('optconf method must be xtb or pm6')
    config = replace(config, relaxation=resolve_relaxation(config.relaxation, method))
    if isinstance(cores, bool) or not isinstance(cores, int) or cores < 1:
        raise ValueError('nproc must be a positive integer')
    if config.relaxation['implementation'] == 'graybox':
        if adapter is not None:
            raise ValueError('Continuous adapters cannot be used for graybox blocks')
        return optimize_graybox_candidates(prepared, reference, config, charge, multiplicity, cores, memory, method=method)
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


def optimize_graybox_candidates(prepared, reference, config, charge, multiplicity, cores, memory,
                               *, method='pm6', candidate_factory=None):
    """Prepare native candidates, run the shared scheduler, then audit/select.

    nproc bounds concurrent one-core candidates. All policies use the same batch
    executor and result accounting; structure eligibility remains unchanged.
    """
    folder = prepared.status_path.parent / 'electronic'
    folder.mkdir(exist_ok=False)
    started = time.monotonic()
    settings = config.relaxation
    workers = resolve_block_workers(cores, len(prepared.candidates), settings.get('parallel_candidates', 'auto'))
    if candidate_factory is None:
        candidate_factory = PM6Candidate if method == 'pm6' else XTBCandidate
    force_kind = 'mean_atom_force' if method == 'pm6' else XTB_FORCE_KIND
    peak_memory_mb = 0.
    runners = {}
    for i, mol in enumerate(prepared.candidates):
        key = f'candidate_{i:05d}'
        directory = folder / key
        directory.mkdir()
        write_sdf(directory / 'initial.sdf', [mol])
        runners[key] = candidate_factory(mol, directory, charge, multiplicity, memory, settings)
    write_json(folder / 'settings.json', dict(settings, backend=method, allocated_cores=cores,
               parallel_workers=workers, cores_per_calculation=1, force_kind=force_kind))

    def save_event(event):
        # Small JSONL events remain readable even after an interrupted search.
        with (folder / 'events.jsonl').open('a') as stream:
            stream.write(json.dumps(event, allow_nan=False) + '\n')
        write_json(folder / 'progress.json', dict(event, initial_candidates=len(runners)))

    def check_pool_memory():
        nonlocal peak_memory_mb
        if method == 'xtb':
            peak_memory_mb = max(peak_memory_mb, pool_memory_mb(runners.values()))
            if peak_memory_mb > settings.get('xtb_pool_memory_mb', 8192):
                raise MemoryError('xTB retained-process pool RSS exceeded xtb_pool_memory_mb; no candidates were dropped or restarted')

    atoms = {key: prepared.candidates[i].GetNumAtoms() for i, key in enumerate(runners)}
    executor = CandidateBlockExecutor(runners, workers, check_pool_memory)
    search = GrayboxSearch(atoms, None, settings, save_event, force_kind=force_kind,
                           workers=workers, advance_batch=executor.advance_batch)
    with cleanup_on_signals():
        try:
            with executor:
                search.run()
        finally:
            # Release resident native processes even after quota, timeout or a
            # persistence error. PM6 has no resident subprocess between blocks.
            try:
                if method == 'xtb':
                    cleanup_errors = {}
                    for key, runner in runners.items():
                        try:
                            runner.close()
                        except Exception as exc:
                            cleanup_errors[key] = f'{type(exc).__name__}: {exc}'
                    if cleanup_errors:
                        write_json(folder/'cleanup_errors.json', cleanup_errors)
                        raise RuntimeError('Some xTB children could not be cleaned up; see cleanup_errors.json')
            finally:
                write_json(folder / 'search_state.json', dict(
                    stop_reason=search.stop_reason or 'interrupted', cost=search.cost,
                    required_converged=search.required, initial_candidates=search.initial_count,
                    peak_pool_rss_mb=peak_memory_mb if method == 'xtb' else None,
                    states={key: asdict(state) for key, state in search.states.items()}))
    records, energies, rows = [], [], []
    for i, (key, state) in enumerate(search.states.items()):
        mol = runners[key].final_molecule if state.status == 'converged' else None
        if mol is not None:
            clear_structure_audit(mol)
            mol.SetDoubleProp('energy_hartree', state.energy)
        records.append(mol)
        energies.append(state.energy if mol is not None else None)
        row = dict(candidate_id=prepared.candidates[i].GetProp('candidate_id'), prepared_index=i,
                   state=state.status, energy_hartree=state.energy, error=state.error,
                   evaluations=state.evaluations, calls=state.calls, wall_seconds=state.wall_seconds)
        rows.append(row)
        write_json(folder / key / 'status.json', row)
    audit = evaluate_optimized_candidates(records, energies, reference, config.validation)
    audit.update(implementation='graybox', algorithm=settings['algorithm'], backend=method,
        score_definition=force_kind, algorithm_label='laqa_norm' if method == 'xtb' and settings['algorithm']=='laqa' else settings['algorithm'],
        peak_pool_rss_mb=peak_memory_mb if method == 'xtb' else None,
        algorithm_variant=('sequential_budget_extension_with_reentry' if workers == 1
                           else 'batched_budget_extension_with_reentry'), energy_unit='hartree',
        convergence_fraction=settings['convergence_fraction'], required_converged=search.required,
        stop_reason=search.stop_reason, quota_reached=search.done(), cost=search.cost,
        input_candidates=len(rows), attempted_candidates=sum(row['calls'] > 0 for row in rows),
        converged_candidates=sum(row['state'] == 'converged' for row in rows),
        failed_candidates=sum(row['state'] == 'failed' for row in rows),
        limit_candidates=sum(row['state'] == 'limit' for row in rows),
        unfinished_candidates=sum(row['state'] in ('pending', 'paused') for row in rows),
        candidate_runs=rows, wall_seconds=time.monotonic() - started,
        sum_candidate_wall_seconds=search.cost['wall_seconds'], allocated_cores=cores,
        cores_per_calculation=1, parallel_workers=workers,
        scheduling='sequential_graybox_1_core' if workers == 1 else 'parallel_graybox_1_core')
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


def _initialize_native_worker(cpu_sets):
    """Disjoint Linux affinity also keeps Gaussian's scheduler wrapper bounded."""
    if cpu_sets is not None:
        os.sched_setaffinity(0, cpu_sets.get())
