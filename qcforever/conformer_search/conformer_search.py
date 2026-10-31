"""QCforever entry point: prepare candidates, optimize them, and return a summary.

Start with run_conformer_search. The detailed generation and optimization records
remain under conformer_search; optimized_structures.sdf is the QCforever handoff.
"""
import json
import shutil
from contextlib import chdir
from pathlib import Path

from .settings import SearchConfig
from .generate_conformers import prepare_candidates
from .optimize_conformers import optimize_candidates
from .relaxation_history import compact_relaxation_outputs, archive_relaxation_details
from .structure_file_io import read_input_structure
from .calculation_logs import write_json


def run_conformer_search(infilename, charge, multiplicity, method, nproc, memory, config):
    """Called by optconf with resolved defaults and optional YAML overrides.

    This adapter uses process-local working directories, like legacy QCforever;
    run independent molecules in separate processes rather than Python threads.
    """
    if not isinstance(config, SearchConfig):
        config = SearchConfig.load(config) if isinstance(config, (str, Path)) else SearchConfig.from_mapping(config)
    if method not in ('xtb', 'pm6'):
        raise ValueError('optconf backend must be xtb or pm6')
    if config.profile == 'legacy':
        return _run_legacy_conformer_search(infilename, charge, multiplicity, method, nproc, memory)
    config = config.for_backend(method)
    reference = read_input_structure(infilename, charge)
    root = Path.cwd()/'conformer_search'
    prepared = prepare_candidates(reference, root, config, allocated_cores=nproc)
    if prepared.initial_sdf is None:
        raise RuntimeError(f'Conformer preparation stopped: {prepared.state}; see {prepared.status_path}')
    # Preparation created a fresh run directory. An existing run is rejected
    # before this point, so cleanup cannot remove an earlier run's checkpoints.
    try:
        audit = optimize_candidates(prepared, reference, config, charge, multiplicity, method, nproc, memory)
        return _save_search_summary(root, prepared, config, method, audit)
    finally:
        # Both execution modes have joined their workers before returning or
        # raising. Paused candidates no longer need their restart files here.
        try:
            report = compact_relaxation_outputs(root / 'electronic', method, keep_logs=config.keep_logs)
            if report is not None and not report['errors']:
                archive_relaxation_details(root / 'electronic')
        except Exception as exc:
            # Storage cleanup must not change a scientific result or hide its error.
            print(f'WARNING: Relaxation history cleanup did not finish: {exc}')


def _run_legacy_conformer_search(infilename, charge, multiplicity, method, nproc, memory):
    """Run the original FAFOOM generation and LAQA relaxation unchanged.

    Legacy initialization appends to its SDF, so require a fresh working
    directory. Only its newly produced result may become the DFT handoff.
    The Gaussian/GAMESS caller performs the usual read_selected_structure
    check; modern candidate filtering and convergence quotas do not apply.
    """
    from qcforever.laqa_fafoom.laqa_confopt_QCforever import LAQA_confopt_main

    source = Path(infilename).resolve()
    root = Path.cwd() / 'conformer_search'
    legacy_directory = root / 'legacy'
    legacy_directory.mkdir(parents=True, exist_ok=False)
    summary = {
        'state': 'succeeded', 'profile': 'legacy', 'backend': method,
        'relaxation_implementation': 'legacy',
        'relaxation_algorithm': 'legacy_laqa',
        'structure_check_state': 'not_run_legacy',
        'details_directory': 'conformer_search',
        'legacy_directory': 'conformer_search/legacy',
    }
    try:
        with chdir(legacy_directory):
            LAQA_confopt_main(str(source), charge, multiplicity, method, nproc, memory)
        optimized = legacy_directory / 'optimized_structures.sdf'
        if not optimized.is_file() or optimized.stat().st_size == 0:
            raise RuntimeError('Legacy FAFOOM/LAQA did not produce an optimized structure')
        shutil.copyfile(optimized, root.parent / 'optimized_structures.sdf')
    except Exception as exc:
        summary.update(state='failed', error=f'{type(exc).__name__}: {exc}')
        write_json(root / 'summary.json', summary)
        raise
    write_json(root / 'summary.json', summary)
    return summary


def _save_search_summary(root, prepared, config, method, audit):
    """Keep the existing QCforever result keys and selected-candidate semantics."""
    status = json.loads(prepared.status_path.read_text())
    index = audit['selected_index']
    search_warning = 'convergence_target_not_reached' if audit.get('quota_reached') is False else None
    state = 'succeeded'
    if search_warning:
        state = 'succeeded_with_search_warning'
    if audit['selected_structure_warning']:
        state = 'succeeded_with_structure_warning'
    summary = {'state': state,
               'search_warning': search_warning,
               'structure_check_warning': audit['selected_structure_warning'],
               'profile': config.profile, 'backend': method,
               'maximum_candidates': status['budget']['maximum_candidates'],
               'input_candidates': len(prepared.candidates),
               'relaxed_candidates': audit['attempted_candidates'],
               'attempted_candidates': audit['attempted_candidates'],
               'converged_candidates': audit['converged_candidates'],
               'failed_candidates': audit['failed_candidates'],
               'limit_candidates': audit.get('limit_candidates', 0),
               'unfinished_candidates': audit.get('unfinished_candidates', 0),
               'relaxation_implementation': audit['implementation'],
               'relaxation_algorithm': audit.get('algorithm'),
               'relaxation_algorithm_label': audit.get('algorithm_label', audit.get('algorithm')),
               'force_definition': audit.get('score_definition'),
               'relaxation_stop_reason': audit.get('stop_reason', 'all_candidates_attempted'),
               'convergence_fraction': audit.get('convergence_fraction'),
               'convergence_target_reached': audit.get('quota_reached'),
               'relaxation_cost': audit.get('cost'),
               'primary_valid_candidates': audit['primary_valid_candidates'],
               'selected_candidate_id': audit['candidate_runs'][index]['candidate_id'],
               'energy_hartree': audit['candidates'][index]['energy_hartree'],
               'mm_method': config.mm_method, 'mm_state': status['mm_state'],
               'mm_skip_reason': status.get('mm_skip_reason'),
               'generation_seconds': sum(s['wall_seconds'] for s in status['stages']),
               'mm_seconds': status['mm_wall_seconds'], 'relaxation_seconds': audit['wall_seconds'],
               'details_directory': 'conformer_search',
               'generation_stages': [{'generator': s['generator'], 'state': s['state'],
                                      'requested_raw': s['requested_raw'], 'returned_raw': s['returned_raw']}
                                     for s in status['stages']]}
    write_json(root/'summary.json', summary)
    if summary['structure_check_warning']:
        print('WARNING: Conformer search selected a structure with '
              f"{summary['structure_check_warning']}. It will still be used for the subsequent calculation. "
              f"Check {root/'electronic/audit.json'} before interpreting the result.")
    if summary['search_warning']:
        print('WARNING: Conformer search did not reach its convergence target. '
              f"A converged candidate was selected; see {root/'summary.json'}.")
    return summary


def read_selected_structure(summary):
    """QC handoff: distinguish a successful search from failed readback.

    Mutate the returned summary so the caller's optconf flag and detailed
    status agree. Never use an old optimized_structures.sdf after search failure.
    """
    from qcforever.util import read_mol_file
    if not summary or summary.get('state') == 'failed':
        raise RuntimeError((summary or {}).get('error', 'No successful conformation search'))
    try:
        result = read_mol_file.read_sdf('./optimized_structures.sdf')
    except Exception as exc:
        summary.update(search_state=summary['state'], state='failed',
                       failure_stage='selected_structure_readback',
                       error=f'{type(exc).__name__}: {exc}',
                       structure_handoff={'state': 'failed', 'path': 'optimized_structures.sdf'})
        if Path('conformer_search').is_dir():
            write_json(Path('conformer_search/summary.json'), summary)
        raise
    summary['structure_handoff'] = {'state': 'succeeded', 'path': 'optimized_structures.sdf'}
    if Path('conformer_search').is_dir():
        write_json(Path('conformer_search/summary.json'), summary)
    return result
