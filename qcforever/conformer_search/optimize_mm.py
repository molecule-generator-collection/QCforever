"""MM parameters must cover the whole molecule; never silently switch MM."""
import time

from rdkit import Chem
from rdkit.Chem import AllChem
from .check_structures import copy_with_coordinates


class UnsupportedParametersError(RuntimeError):
    pass


class MissingMMDependency(RuntimeError):
    pass


def optimize_etkdg_candidates(candidates, config, *, optimizer=None):
    """Apply MM only to ETKDG candidates; return structures, status, and audit.

    Model-generated candidates bypass MM. Unsupported parameters skip the
    requested force field without silently substituting another one.
    """
    details = {}
    started = time.monotonic()
    # Only ETKDG candidates receive MM. Learned coordinates bypass MM even
    # when the pool contains both learned and ETKDG-generated candidates.
    optimized = [Chem.Mol(m) for m in candidates]
    codes = [None]*len(candidates)
    mm_indices = [i for i, mol in enumerate(candidates) if mol.GetProp('generator') == 'etkdgv3']
    mm_runs = []
    mm_attempted = False
    for i in mm_indices:
        candidate_result = {'candidate_index': i, 'candidate_id': candidates[i].GetProp('candidate_id')}
        if config.mm_method == 'none':
            candidate_result.update(state='not_requested', optimizer_status=None)
            mm_runs.append(candidate_result)
            continue
        mm_attempted = True
        try:
            subset, subset_codes = (optimizer or optimize_with_mm)(
                [Chem.Mol(candidates[i])], config.mm_method,
                config.mm.get(config.mm_method, {}))
            if len(subset) != 1 or len(subset_codes) != 1:
                raise ValueError('MM optimizer changed candidate count')
            optimized[i], codes[i] = subset[0], subset_codes[0]
            candidate_result.update(state='completed', optimizer_status=codes[i])
        except (UnsupportedParametersError, MissingMMDependency, RuntimeError) as exc:
            candidate_result.update(state='skipped', reason=f'{type(exc).__name__}: {exc}')
        mm_runs.append(candidate_result)
    skipped = [r for r in mm_runs if r['state'] == 'skipped']
    details['mm_state'] = ('not_requested' if config.mm_method == 'none' else
                          'not_applicable_to_learned_generators' if not mm_indices else
                          'skipped' if len(skipped) == len(mm_indices) else
                          'partially_skipped' if skipped else 'completed')
    details['mm_candidate_runs'] = mm_runs
    if skipped:
        details['mm_skip_reason'] = '; '.join(dict.fromkeys(r['reason'] for r in skipped))
    details['mm_candidate_indices'] = mm_indices
    details['mm_scope'] = 'etkdgv3_only'
    details['mm_wall_seconds'] = time.monotonic()-started
    if len(optimized) != len(codes) or len(optimized) != len(candidates):
        raise ValueError('MM optimizer changed candidate count')

    details['mm_optimizer_status'] = codes
    return optimized, mm_attempted, details


def optimize_with_mm(records, method, options):
    if method == 'none':
        return [Chem.Mol(m) for m in records], [None]*len(records)
    limit = options.get('maximum_iterations', 5000)
    if not isinstance(limit, int) or limit < 1:
        raise ValueError('MM maximum_iterations must be a positive integer')
    if method in ('uff', 'mmff94s'):
        has = AllChem.UFFHasAllMoleculeParams if method == 'uff' else AllChem.MMFFHasAllMoleculeParams
        result, statuses = [], []
        for source in records:
            # Even parameter-availability checks can alter RDKit's aromaticity
            # flags. Never pass the preserved source to a force-field API.
            mol = Chem.Mol(source)
            if not has(mol):
                raise UnsupportedParametersError(f'{method} does not cover all molecular parameters')
            if method == 'uff':
                status = AllChem.UFFOptimizeMolecule(mol, maxIters=limit)
            else:
                props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94s')
                ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=0)
                status = ff.Minimize(maxIts=limit)
            result.append(copy_with_coordinates(source, mol))
            statuses.append(int(status))
        return result, statuses
    raise ValueError(f'Unknown MM method: {method}')
