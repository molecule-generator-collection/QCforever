"""Bounded incremental generation, common validation, merging, then MM."""
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import time

from rdkit import Chem, rdBase
from .config import SearchConfig
from .generators import GenerationError, MissingGeneratorDependency, generate_batch, resolve_model_device
from .preoptimization import MissingMMDependency, UnsupportedParametersError, preoptimize
from .validation import IncrementalCandidateFilter, filter_candidates


def write_sdf(path, records):
    with Chem.SDWriter(str(path)) as writer:
        for mol in records:
            if mol is not None:
                writer.write(mol)


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


@dataclass
class PreparationResult:
    state: str
    candidates: list
    initial_sdf: Path | None
    status_path: Path


def prepare_candidates(molecule, output_directory, config, *, allocated_cores=1,
                       generators=None, mm_optimizer=None):
    """Request the unfilled quota first, then worker-sized batches, at most 2N/stage.

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
    workers = config.parallelism(allocated_cores)
    root = Path(output_directory).resolve()
    root.mkdir(parents=True, exist_ok=False)
    save(root/'resolved_config.json', config.to_mapping())
    write_sdf(root/'reference.sdf', [reference])
    started = time.monotonic()
    status_path = root/'status.json'
    status = {'state': 'running', 'profile': config.profile, 'budget': budget,
              'reference_smiles': Chem.MolToSmiles(Chem.RemoveHs(reference)),
              'rdkit_version': rdBase.rdkitVersion, 'stages': [],
              'workers': workers, 'threads_per_worker': config.threads,
              'candidate_retention': 'merge', 'stereo_generation_policy': 'require_input_specified_stereo',
              'selected_mm': config.mm_method}
    combined = []
    generation_filter = IncrementalCandidateFilter(reference, maximum, config.validation)
    sessions = []
    try:
        for stage_index, (name, _) in enumerate(config.stages()):
            if len(combined) >= maximum:
                break
            directory = root/f'{stage_index:02d}_{name}'
            directory.mkdir()
            tick = time.monotonic()
            row = {'generator': name, 'state': 'running', 'requested_raw': 0,
                   'returned_raw': 0, 'maximum_raw_candidates': cap, 'batches': []}
            status['stages'].append(row)
            options = dict(config.generators.get(name, {}))
            stage_workers = workers
            session = None
            if not (generators or {}).get(name) and options.get('persistent'):
                from .model_session import ModelSession
                session = ModelSession(directory, options, config.threads)
                sessions.append(session)
            while len(combined) < maximum and row['requested_raw'] < cap:
                bi = len(row['batches'])
                remaining = maximum - len(combined)
                count = remaining if bi == 0 else min(stage_workers, remaining)
                count = min(count, cap-row['requested_raw'])
                folder = directory/f'batch_{bi:03d}'
                folder.mkdir()
                seed = (config.seed+stage_index*1000003+bi*1009) % 2**31
                batch = {'seed': seed, 'requested_raw': count, 'state': 'running'}
                row['batches'].append(batch)
                row['requested_raw'] += count
                bt = time.monotonic()
                try:
                    adapter = (generators or {}).get(name)
                    if bi == 0 and name in ('ditmc', 'torsional_diffusion'):
                        requested_device = options.get('device', config.device)
                        device = ({'requested': requested_device, 'effective': 'cpu',
                                   'reason': 'injected_test_adapter'} if adapter else
                                  resolve_model_device(options, requested_device, directory, config.threads))
                        options['device'] = device['effective']
                        stage_workers = 1 if device['effective'] == 'gpu' else workers
                        row.update(device=device, workers=stage_workers)
                    if adapter:
                        raw = list(adapter(reference, count, seed, folder,
                                           config.generators.get(name, {}), config.threads))
                    else:
                        if session is not None:
                            raw = session.generate(reference, count, seed, folder, stage_workers)
                        else:
                            raw = generate_batch(name, reference, count, seed, folder,
                                                 options, config.threads, stage_workers)
                    write_sdf(folder/'raw_candidates.sdf', raw)
                    if len(raw) > count:
                        raise GenerationError('Returned candidates exceed requested batch')
                    for i, mol in enumerate(raw):
                        if mol is not None:
                            mol.SetProp('candidate_id', f'{name}:b{bi:03d}:r{i:05d}')
                            mol.SetProp('generator', name)
                            mol.SetIntProp('generator_seed', seed)
                    combined, audit = generation_filter.extend(raw)
                    row['returned_raw'] += len(raw)
                    batch.update(state='completed', returned_raw=len(raw), merge_audit=audit,
                                 accumulated_valid=len(combined))
                    write_sdf(folder/'accepted_pool.sdf', combined)
                except MissingGeneratorDependency as exc:
                    batch.update(state='unavailable_dependency', error=str(exc))
                    row['state'] = batch['state']
                    break
                except GenerationError as exc:
                    batch.update(state='generation_failed', error=str(exc))
                    row['state'] = batch['state']
                    break
                finally:
                    batch['wall_seconds'] = time.monotonic()-bt
                    save(folder/'status.json', batch)
                    save(status_path, status)
            if row['state'] == 'running':
                row['state'] = 'target_reached' if len(combined) >= maximum else 'trial_cap_reached'
            row.update(route_candidate_count=len(combined), wall_seconds=time.monotonic()-tick)
            if session is not None:
                session.close()
                sessions.remove(session)
            save(directory/'status.json', row)
        write_sdf(root/'generated_candidates.sdf', combined)
        if not combined:
            status['state'] = 'generation_failed'
            return PreparationResult(status['state'], [], None, status_path)
        tick = time.monotonic()
        # Only ETKDG candidates receive MM. Learned coordinates bypass MM even
        # when the pool contains both learned and ETKDG-generated candidates.
        optimized = [Chem.Mol(m) for m in combined]
        codes = [None]*len(combined)
        mm_indices = [i for i, mol in enumerate(combined) if mol.GetProp('generator') == 'etkdgv3']
        mm_runs = []
        for i in mm_indices:
            row = {'candidate_index': i, 'candidate_id': combined[i].GetProp('candidate_id')}
            try:
                subset, subset_codes = (mm_optimizer or preoptimize)(
                    [Chem.Mol(combined[i])], config.mm_method,
                    config.mm.get(config.mm_method, {}))
                if len(subset) != 1 or len(subset_codes) != 1:
                    raise ValueError('MM optimizer changed candidate count')
                optimized[i], codes[i] = subset[0], subset_codes[0]
                row.update(state='not_requested' if config.mm_method == 'none' else 'completed',
                           optimizer_status=codes[i])
            except (UnsupportedParametersError, MissingMMDependency, RuntimeError) as exc:
                row.update(state='skipped', reason=f'{type(exc).__name__}: {exc}')
            mm_runs.append(row)
        skipped = [r for r in mm_runs if r['state'] == 'skipped']
        status['mm_state'] = ('not_requested' if config.mm_method == 'none' else
                              'not_applicable_to_learned_generators' if not mm_indices else
                              'skipped' if len(skipped) == len(mm_indices) else
                              'partially_skipped' if skipped else 'completed')
        status['mm_candidate_runs'] = mm_runs
        if skipped:
            status['mm_skip_reason'] = '; '.join(dict.fromkeys(r['reason'] for r in skipped))
        status['mm_candidate_indices'] = mm_indices
        status['mm_scope'] = 'etkdgv3_only'
        status['mm_wall_seconds'] = time.monotonic()-tick
        if len(optimized) != len(codes) or len(optimized) != len(combined):
            raise ValueError('MM optimizer changed candidate count')
        write_sdf(root/'mm_candidates.sdf', optimized)
        filtered, audit = filter_candidates(optimized, reference, maximum, config.validation)
        status.update(mm_optimizer_status=codes, post_mm_audit=audit,
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
        for session in sessions:
            session.close()
        status['total_wall_seconds'] = time.monotonic()-started
        save(status_path, status)
