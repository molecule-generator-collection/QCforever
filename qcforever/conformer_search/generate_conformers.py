"""Prepare conformers: generate in profile order, check candidates, then apply MM.

Read prepare_candidates first. Each generator gets a fixed attempt cap;
accepted candidates survive fallback to the next generator.
"""
from concurrent.futures import ThreadPoolExecutor
from contextvars import copy_context
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

from .settings import SearchConfig
from .calculation_logs import write_json
from .structure_file_io import write_sdf
from .optimize_mm import optimize_etkdg_candidates
from .check_structures import IncrementalCandidateFilter, filter_candidates


class GenerationError(RuntimeError):
    pass


class MissingGeneratorDependency(GenerationError):
    pass


@dataclass
class PreparationResult:
    state: str
    candidates: list
    initial_sdf: Path | None
    status_path: Path


def prepare_candidates(molecule, output_directory, config, *, allocated_cores=1,
                       generators=None, mm_optimizer=None):
    """Generate, validate, and MM-prepare a bounded candidate pool.

    Request the unfilled quota first, then worker-sized batches, at most 2N/stage.

    Requested attempts count even if coordinates are not returned. All stages
    use the same validator; accepted candidates survive subsequent stages.
    """
    if not isinstance(config, SearchConfig):
        config = SearchConfig.resolve(override=config)
    reference = Chem.MolFromSmiles(molecule) if isinstance(molecule, str) else Chem.Mol(molecule)
    if reference is None:
        raise ValueError('Invalid molecule')
    reference = Chem.AddHs(Chem.RemoveHs(reference), addCoords=True)
    budget = config.budget.resolve(Chem.RemoveHs(reference))
    maximum, cap = budget['maximum_candidates'], budget['raw_attempt_cap_per_generator']
    threads = config.effective_threads(allocated_cores)
    workers = config.parallelism(allocated_cores)
    root = Path(output_directory).resolve()
    root.mkdir(parents=True, exist_ok=False)
    write_json(root/'resolved_config.json', config.to_mapping())
    write_sdf(root/'reference.sdf', [reference])
    started = time.monotonic()
    status_path = root/'status.json'
    status = {'state': 'running', 'profile': config.profile, 'budget': budget,
              'reference_smiles': Chem.MolToSmiles(Chem.RemoveHs(reference)),
              'rdkit_version': rdBase.rdkitVersion, 'stages': [],
              'allocated_cores': allocated_cores,
              'requested_threads_per_worker': config.threads,
              'workers': workers, 'threads_per_worker': threads,
              'candidate_retention': 'merge', 'stereo_generation_policy': 'require_input_specified_stereo',
              'selected_mm': config.mm_method}
    combined = []
    generation_filter = IncrementalCandidateFilter(reference, maximum, config.validation)
    try:
        for stage_index, (name, _) in enumerate(config.stages()):
            if len(combined) >= maximum:
                break
            combined = _generate_with_method(
                name, stage_index, reference, combined, generation_filter,
                root=root, status=status, config=config, generators=generators)
        write_sdf(root/'generated_candidates.sdf', combined)
        if not combined:
            status['state'] = 'generation_failed'
            return PreparationResult(status['state'], [], None, status_path)
        optimized, mm_attempted, mm_details = optimize_etkdg_candidates(
            combined, config, optimizer=mm_optimizer)
        status.update(mm_details)
        write_sdf(root/'mm_candidates.sdf', optimized)
        validation_started = time.monotonic()
        # Follow the execution route, not a structural-change detector. Even
        # skipped or unchanged MM outputs take the full post-MM validation path.
        if mm_attempted:
            filtered, audit = filter_candidates(optimized, reference, maximum, config.validation)
            status['post_mm_validation_mode'] = 'full_recheck'
        else:
            filtered, audit = generation_filter.snapshot()
            status['post_mm_validation_mode'] = 'reused_generation_no_mm'
        status['post_mm_validation_wall_seconds'] = time.monotonic()-validation_started
        status.update(post_mm_audit=audit,
                      final_candidate_count=len(filtered))
        if not filtered:
            status['state'] = 'no_valid_candidates_after_mm'
            return PreparationResult(status['state'], [], None, status_path)
        initial = root/'initial_structures.sdf'
        write_sdf(initial, filtered)
        status.update(state='ready_for_electronic_optimization', initial_sdf=str(initial),
                      initial_sdf_sha256=hashlib.sha256(initial.read_bytes()).hexdigest())
        return PreparationResult(status['state'], filtered, initial, status_path)
    except Exception as exc:
        status.update(state='failed', error=f'{type(exc).__name__}: {exc}')
        raise
    finally:
        status['total_wall_seconds'] = time.monotonic()-started
        write_json(status_path, status)


def _generate_with_method(name, stage_index, reference, combined, generation_filter,
                          *, root, status, config, generators):
    """Fill the remaining quota with one generator, without discarding prior candidates."""
    maximum = status['budget']['maximum_candidates']
    cap = status['budget']['raw_attempt_cap_per_generator']
    threads, workers = status['threads_per_worker'], status['workers']
    status_path = root/'status.json'
    directory = root/f'{stage_index:02d}_{name}'
    directory.mkdir()
    stage_started = time.monotonic()
    stage_status = {'generator': name, 'state': 'running', 'requested_raw': 0,
           'returned_raw': 0, 'maximum_raw_candidates': cap, 'batches': []}
    status['stages'].append(stage_status)
    options = dict(config.generators.get(name, {}))
    stage_workers = workers
    session = None
    if not (generators or {}).get(name) and options.get('persistent'):
        from .model_workers.model_process import ModelSession
        session = ModelSession(directory, options, threads)
    try:
        while len(combined) < maximum and stage_status['requested_raw'] < cap:
            batch_index = len(stage_status['batches'])
            remaining = maximum - len(combined)
            count = remaining if batch_index == 0 else min(stage_workers, remaining)
            count = min(count, cap-stage_status['requested_raw'])
            folder = directory/f'batch_{batch_index:03d}'
            folder.mkdir()
            seed = (config.seed+stage_index*1000003+batch_index*1009) % 2**31
            batch = {'seed': seed, 'requested_raw': count, 'state': 'running'}
            stage_status['batches'].append(batch)
            stage_status['requested_raw'] += count
            batch_started = time.monotonic()
            try:
                adapter = (generators or {}).get(name)
                if batch_index == 0 and name in ('ditmc', 'torsional_diffusion'):
                    requested_device = options.get('device', config.device)
                    device = ({'requested': requested_device, 'effective': 'cpu',
                               'reason': 'injected_test_adapter'} if adapter else
                              resolve_model_device(options, requested_device, directory, threads))
                    options['device'] = device['effective']
                    stage_workers = 1 if device['effective'] == 'gpu' else workers
                    stage_status.update(device=device, workers=stage_workers)
                raw = _generate_raw_batch(
                    name, reference, count, seed, folder, options, threads, stage_workers,
                    session, adapter, config.generators.get(name, {}), batch_index)
                combined, audit = generation_filter.extend(raw)
                stage_status['returned_raw'] += len(raw)
                batch.update(state='completed', returned_raw=len(raw), merge_audit=audit,
                             accumulated_valid=len(combined))
                write_sdf(folder/'accepted_pool.sdf', combined)
            except MissingGeneratorDependency as exc:
                batch.update(state='unavailable_dependency', error=str(exc))
                stage_status['state'] = batch['state']
                break
            except GenerationError as exc:
                batch.update(state='generation_failed', error=str(exc))
                stage_status['state'] = batch['state']
                break
            finally:
                batch['wall_seconds'] = time.monotonic()-batch_started
                write_json(folder/'status.json', batch)
                write_json(status_path, status)
        if stage_status['state'] == 'running':
            stage_status['state'] = 'target_reached' if len(combined) >= maximum else 'trial_cap_reached'
        stage_status.update(route_candidate_count=len(combined), wall_seconds=time.monotonic()-stage_started)
        write_json(directory/'status.json', stage_status)

        return combined
    finally:
        if session is not None:
            session.close()


def _generate_raw_batch(name, reference, count, seed, folder, options, threads,
                        workers, session, adapter, adapter_options, batch_index):
    """Generate and label raw structures before the common chemical checks."""
    if adapter:
        raw = list(adapter(reference, count, seed, folder,
                           adapter_options, threads))
    else:
        if session is not None:
            raw = session.generate(reference, count, seed, folder, workers)
        else:
            raw = generate_batch(name, reference, count, seed, folder,
                                 options, threads, workers)
    write_sdf(folder/'raw_candidates.sdf', raw)
    if len(raw) > count:
        raise GenerationError('Returned candidates exceed requested batch')
    for i, mol in enumerate(raw):
        if mol is not None:
            mol.SetProp('candidate_id', f'{name}:b{batch_index:03d}:r{i:05d}')
            mol.SetProp('generator', name)
            mol.SetIntProp('generator_seed', seed)
    return raw


def generate_batch(name, reference, count, seed, directory, options, threads, workers):
    """Parallel model subprocesses; deterministic merge independent of completion order.

    Threads only supervise isolated command processes. No model is loaded in the
    QCforever interpreter. Native ETKDG uses its own bounded RDKit thread pool.
    """
    if name == 'etkdgv3':
        return generate(name, reference, count, seed, directory, options, threads*workers)
    if not options.get('command'):
        raise MissingGeneratorDependency(f'{name}: configure a raw-SDF command')
    slots = min(workers, count)
    requests = []
    for i in range(slots):
        folder = directory/f'worker_{i:02d}'
        folder.mkdir()
        size = count//slots + (i < count % slots)
        requests.append((reference, size, (seed+i) % 2**31, folder, options, threads))
    with ThreadPoolExecutor(max_workers=slots) as pool:
        futures = [pool.submit(copy_context().run, command_generator, *args) for args in requests]
        raw = []
        failures = []
        for future, args in zip(futures, requests):
            try:
                records = future.result()
                if len(records) > args[1]:
                    raise GenerationError('Worker returned more records than requested')
                raw.extend(records)
            except GenerationError as exc:
                failures.append(str(exc))
                (args[3]/'error.json').write_text(json.dumps({'error': str(exc)}))
        if len(failures) == slots:
            raise GenerationError(f'All generation workers failed: {failures}')
        return raw

def generate(name, reference, attempt_cap, seed, directory, options, threads):
    if options.get('command'):
        return command_generator(reference, attempt_cap, seed, directory, options, threads)
    if name == 'etkdgv3':
        return etkdgv3(reference, attempt_cap, seed, directory, options, threads)
    if name == 'etflow':
        return etflow(reference, attempt_cap, seed, directory, options, threads)
    if name in ('ditmc', 'torsional_diffusion'):
        raise MissingGeneratorDependency(f'{name} requires an optional environment and raw-SDF command')
    raise ValueError(f'Unknown generator: {name}')


def etkdgv3(reference, attempt_cap, seed, directory, options, threads):
    mol = Chem.AddHs(Chem.RemoveHs(reference))
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.numThreads = threads
    params.pruneRmsThresh = -1.0
    params.enforceChirality = True
    AllChem.EmbedMultipleConfs(mol, numConfs=attempt_cap, params=params)
    return split_conformers(mol)


def command_generator(reference, attempt_cap, seed, directory, options, threads):
    """Optional model environment writes raw SDF; it must not filter or optimize.

    Command tokens may contain {request}, {output}, {input}. No shell expansion.
    The caller must supply a generator command; unconfigured models are reported
    as unavailable rather than being replaced behind the caller's back.
    """
    command = options.get('command')
    if not isinstance(command, list) or not command or not all(isinstance(x, str) for x in command):
        raise MissingGeneratorDependency('An explicit list-valued generator command is required')
    request = directory/'request.json'
    output = directory/'external_raw.sdf'
    inp = directory/'reference.sdf'
    with Chem.SDWriter(str(inp)) as writer:
        writer.write(reference)
    request.write_text(json.dumps({'smiles': Chem.MolToSmiles(Chem.RemoveHs(reference), isomericSmiles=True),
                                   'maximum_raw_candidates': attempt_cap, 'seed': seed,
                                   'threads': threads, 'device': options.get('device', 'cpu'), 'optimization': 'none',
                                   'options': {k:v for k,v in options.items() if k not in ('command', 'timeout_seconds')}}, indent=2))
    values = {'request': str(request), 'output': str(output), 'input': str(inp)}
    argv = [token.format(**values) for token in command]
    try:
        with (directory/'generator.out').open('w') as out:
            from qcforever.util import job_timeout
            env = os.environ.copy()
            for key in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
                env[key] = str(threads)
            job_timeout.run(argv, cwd=directory, env=env, stdout=out, stderr=subprocess.STDOUT,
                           check=True, timeout=options.get('timeout_seconds', 600))
    except (subprocess.SubprocessError, OSError) as exc:
        raise GenerationError(f'Generator command failed: {exc}') from exc
    if not output.is_file():
        raise GenerationError('Generator command did not produce the requested raw SDF')
    # A successful model command may legitimately generate zero candidates.
    # RDKit raises OSError for an empty SDF; return an empty pool explicitly so
    # the fixed attempt cap/fallback policy can handle it, not abort the case.
    if not output.read_text().strip():
        return []
    try:
        return list(Chem.SDMolSupplier(str(output), removeHs=False))
    except (OSError, ValueError) as exc:
        raise GenerationError(f'Unreadable model raw SDF: {exc}') from exc


def resolve_model_device(options, requested, directory, threads):
    """Probe the optional model Python, not QCforever's lightweight environment.

    Explicit cpu/gpu skips discovery. Auto requires the supplied worker's probe
    protocol; probe errors remain visible and are not called CPU successes.
    """
    if requested != 'auto':
        return {'requested': requested, 'effective': requested}
    command = options.get('command')
    if not isinstance(command, list) or not command:
        raise MissingGeneratorDependency('Device discovery requires a model command')
    folder = directory/'device_probe'
    folder.mkdir()
    request, output = folder/'request.json', folder/'device.json'
    request.write_text(json.dumps({'device': 'auto', 'threads': threads}))
    values = {'request': str(request), 'output': str(output), 'input': str(directory.parent/'reference.sdf')}
    argv = [token.format(**values) for token in command] + ['--probe-device']
    env = os.environ.copy()
    for key in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
        env[key] = str(threads)
    try:
        from qcforever.util import job_timeout
        with (folder/'probe.out').open('w') as log:
            job_timeout.run(argv, cwd=folder, env=env, stdout=log, stderr=subprocess.STDOUT,
                            check=True, timeout=options.get('timeout_seconds', 600))
        result = json.loads(output.read_text())
        if result.get('effective') not in ('cpu', 'gpu'):
            raise ValueError('Invalid effective device')
        return result
    except (subprocess.SubprocessError, OSError, ValueError) as exc:
        raise GenerationError(f'Model device discovery failed; see {folder}/probe.out: {exc}') from exc


def etflow(reference, attempt_cap, seed, directory, options, threads):
    try:
        from etflow import BaseFlow
        import torch
    except ImportError as exc:
        raise MissingGeneratorDependency('ET-Flow is not installed in this Python environment') from exc
    if not options.get('checkpoint_cache'):
        raise MissingGeneratorDependency('ET-Flow requires an explicit checkpoint_cache')
    device = options.get('device', 'cpu')
    if device == 'cpu':
        torch.set_num_threads(threads)
    try:
        model = BaseFlow.from_default(model=options.get('model', 'drugs-so3'),
                                      cache=str(Path(options['checkpoint_cache']).expanduser().resolve()))
        model = model.to(device).eval()
        smiles = Chem.MolToSmiles(Chem.RemoveHs(reference), isomericSmiles=True)
        generated = model.predict([smiles], max_batch_size=options.get('batch_size', 4),
                                  num_samples=attempt_cap, seed=seed, device=device, as_mol=True)
        value = generated[smiles]
        return [Chem.Mol(m) for m in value] if isinstance(value, list) else split_conformers(value)
    except (RuntimeError, ValueError, KeyError, OSError) as exc:
        raise GenerationError(f'ET-Flow generation failed: {exc}') from exc


def split_conformers(mol):
    records = []
    for conf in mol.GetConformers():
        record = Chem.Mol(mol)
        record.RemoveAllConformers()
        record.AddConformer(Chem.Conformer(conf), assignId=True)
        records.append(record)
    return records
