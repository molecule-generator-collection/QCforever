"""QCforever entry point: prepare candidates, optimize them, and return a summary.

Start with configured_confopt. The detailed generation and optimization records
remain under conformer_search; optimized_structures.sdf is the QCforever handoff.
"""
import json
from dataclasses import replace
from pathlib import Path

from .settings import SearchConfig, resolve_relaxation
from .generate_conformers import prepare_candidates
from .optimize_conformers import optimize_candidates
from .structure_file_io import read_input_structure
from .calculation_logs import write_json


def configured_confopt(infilename, charge, multiplicity, method, nproc, memory, config):
    """Called by optconf with resolved defaults and optional YAML overrides.

    This adapter uses process-local working directories, like legacy QCforever;
    run independent molecules in separate processes rather than Python threads.
    """
    if not isinstance(config, SearchConfig):
        config = SearchConfig.load(config) if isinstance(config, (str, Path)) else SearchConfig.from_mapping(config)
    if method not in ('xtb', 'pm6'):
        raise ValueError('optconf backend must be xtb or pm6')
    config = replace(config, relaxation=resolve_relaxation(config.relaxation, method))
    reference = read_input_structure(infilename, charge)
    root = Path.cwd()/'conformer_search'
    prepared = prepare_candidates(reference, root, config, allocated_cores=nproc)
    if prepared.initial_sdf is None:
        raise RuntimeError(f'Conformer preparation stopped: {prepared.state}; see {prepared.status_path}')
    audit = optimize_candidates(prepared, reference, config, charge, multiplicity, method, nproc, memory)
    return _save_search_summary(root, prepared, config, method, audit)


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
    return summary


def read_selected_structure(summary):
    """Gaussian handoff: distinguish a successful search from failed readback.

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
