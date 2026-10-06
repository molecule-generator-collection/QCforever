"""Connect configurable preparation/relaxation to the existing optconf handoff."""
from contextlib import contextmanager
import os
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdDetermineBonds

from .config import SearchConfig
from .pipeline import prepare_candidates, save, write_sdf


@contextmanager
def working_directory(path):
    previous = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


def configured_confopt(infilename, charge, multiplicity, method, nproc, memory, config):
    """Called by optconf with resolved defaults and optional YAML overrides.

    This adapter uses process-local working directories, like legacy QCforever;
    run independent molecules in separate processes rather than Python threads.
    """
    if not isinstance(config, SearchConfig):
        config = SearchConfig.load(config) if isinstance(config, (str, Path)) else SearchConfig.from_mapping(config)
    if method not in ('xtb', 'pm6'):
        raise ValueError('optconf backend must be xtb or pm6')
    path = Path(infilename).resolve()
    if path.suffix.lower() == '.sdf':
        records = [m for m in Chem.SDMolSupplier(str(path), removeHs=False) if m is not None]
        if not records:
            raise ValueError('No readable input molecule')
        reference = records[0]
    elif path.suffix.lower() == '.xyz':
        reference = Chem.MolFromXYZFile(str(path))
        if reference is None:
            raise ValueError('Unreadable XYZ input')
        rdDetermineBonds.DetermineBonds(reference, charge=charge)
    else:
        raise ValueError('Configured optconf accepts SDF or XYZ')
    root = Path.cwd()/'conformer_search'
    prepared = prepare_candidates(reference, root, config, allocated_cores=nproc)
    if prepared.initial_sdf is None:
        raise RuntimeError(f'Conformer preparation stopped: {prepared.state}; see {prepared.status_path}')
    from .relaxation import relax_candidates
    audit = relax_candidates(prepared, reference, config, charge, multiplicity, method, nproc, memory)
    import json
    status = json.loads(prepared.status_path.read_text())
    index = audit['selected_index']
    summary = {'state': 'succeeded_with_structure_warning' if audit['selected_structure_warning'] else 'succeeded',
               'structure_check_warning': audit['selected_structure_warning'],
               'profile': config.profile, 'backend': method,
               'maximum_candidates': status['budget']['maximum_candidates'],
               'relaxed_candidates': len(prepared.candidates),
               'attempted_candidates': audit['attempted_candidates'],
               'converged_candidates': audit['converged_candidates'],
               'failed_candidates': audit['failed_candidates'],
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
    save(root/'summary.json', summary)
    return summary
