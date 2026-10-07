"""No-MM snapshots and unconditional post-MM checks have identical decisions."""
from dataclasses import replace
import json
from unittest.mock import patch

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms

from qcforever.conformer_search import pipeline, validation
from qcforever.conformer_search.config import SearchConfig
from qcforever.conformer_search.preoptimization import UnsupportedParametersError


def rotamers():
    mol = Chem.AddHs(Chem.MolFromSmiles('CCO'))
    assert AllChem.EmbedMolecule(mol, randomSeed=42) == 0
    atoms = mol.GetSubstructMatch(Chem.MolFromSmarts('[C]-[C]-[O]-[H]'))
    records = []
    for i, angle in enumerate((180, 60, -60)):
        m = Chem.Mol(mol)
        rdMolTransforms.SetDihedralDeg(m.GetConformer(), *atoms, angle)
        m.SetProp('candidate_id', f'candidate_{i}')
        records.append(m)
    return records


def assert_same(left, right):
    assert len(left) == len(right)
    for a, b in zip(left, right):
        assert a.GetPropsAsDict() == b.GetPropsAsDict()
        assert Chem.MolToSmiles(a) == Chem.MolToSmiles(b)
        np.testing.assert_array_equal(a.GetConformer().GetPositions(), b.GetConformer().GetPositions())


@pytest.mark.parametrize('threshold', [0, .1])
@pytest.mark.parametrize('maximum', [1, 2, 10])
def test_snapshot_matches_full_filter_without_rechecking(threshold, maximum):
    a, b, c = rotamers()
    cfg = replace(SearchConfig.resolve().validation, duplicate_rmsd_angstrom=threshold)
    cache = validation.IncrementalCandidateFilter(a, maximum, cfg)
    cache.extend([None, a, a])
    accepted, _ = cache.extend([b, c])
    expected, audit = validation.filter_candidates(accepted, a, maximum, cfg)
    with patch.object(validation, 'geometry_failure', side_effect=AssertionError('geometry')), \
         patch.object(validation, 'stereo_audit', side_effect=AssertionError('stereo')), \
         patch.object(validation, 'duplicate_rmsd', side_effect=AssertionError('RMSD')):
        actual, reused = cache.snapshot()
    assert audit == reused
    assert_same(actual, expected)
    # Callers cannot corrupt the private pool/audit by changing their copies.
    for mol in (a, accepted[0], actual[0]):
        mol.GetConformer().SetAtomPosition(0, (99, 99, 99))
        mol.SetBoolProp('joint_stereo_match', False)
    reused['candidates'][0]['reason'] = 'corrupted'
    again, again_audit = cache.snapshot()
    assert again_audit == audit
    assert_same(again, expected)


def test_empty_snapshot_matches_full_filter():
    ref = rotamers()[0]
    cfg = SearchConfig.resolve().validation
    cache = validation.IncrementalCandidateFilter(ref, 10, cfg)
    assert cache.snapshot() == validation.filter_candidates([], ref, 10, cfg)


@pytest.mark.parametrize('smiles', ['C[C@H](O)F', 'C/C=C/C'])
def test_snapshot_keeps_coordinate_based_stereo_audit(smiles):
    ref = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(ref, randomSeed=42) == 0
    cfg = SearchConfig.resolve().validation
    cache = validation.IncrementalCandidateFilter(ref, 10, cfg)
    accepted, _ = cache.extend([ref])
    assert len(accepted) == 1
    original, audit = validation.filter_candidates(accepted, ref, 10, cfg)
    result, snapshot_audit = cache.snapshot()
    assert snapshot_audit == audit
    assert_same(result, original)


@pytest.mark.parametrize('mode', ['learned_td', 'learned_ditmc', 'none', 'mixed_none',
    'skipped', 'mixed_skipped', 'mixed_unchanged', 'mm_unchanged', 'mixed_changed',
    'mm_duplicate', 'mm_collision', 'mm_nonfinite', 'mm_connectivity'])
@pytest.mark.parametrize('cores', [1, 3, 8])
def test_pipeline_routes_on_mm_calls_not_coordinate_equality(tmp_path, mode, cores):
    a, b, _ = rotamers()
    profile = 'high' if mode == 'learned_ditmc' else (
        'medium' if mode.startswith('mixed') or mode == 'learned_td' else 'low')
    cfg = SearchConfig.resolve(profile, {'budget': {'formula': 'fixed', 'fixed': 2},
        'threads': 4, 'workers': 8,
        'mm_method': 'none' if mode in ('none', 'mixed_none') else 'mmff94s'})
    def learned(*args):
        assert args[-1] == min(4, cores)
        return [Chem.Mol(a)] if mode.startswith('mixed') else [Chem.Mol(a), Chem.Mol(b)]
    def etkdg(*args):
        assert args[-1] == min(4, cores)
        return [Chem.Mol(b)] if mode.startswith('mixed') else [Chem.Mol(a), Chem.Mol(b)]
    calls = []
    def mm(records, *args):
        calls.append(records[0].GetProp('generator'))
        if mode in ('none', 'mixed_none') or mode.startswith('learned'):
            pytest.fail('No-MM route called the optimizer')
        if mode in ('skipped', 'mixed_skipped'):
            records[0].GetConformer().SetAtomPosition(0, (99, 99, 99))
            raise UnsupportedParametersError('test missing parameters')
        if mode in ('mixed_changed', 'mm_duplicate'):
            return [validation.copy_with_coordinates(records[0], a)], [0]
        if mode == 'mm_collision':
            records[0].GetConformer().SetAtomPosition(1, records[0].GetConformer().GetAtomPosition(0))
        if mode == 'mm_nonfinite':
            records[0].GetConformer().SetAtomPosition(0, (float('nan'), 0, 0))
        if mode == 'mm_connectivity':
            records[0].GetAtomWithIdx(2).SetFormalCharge(1)
        return records, [0]
    with patch.object(pipeline, 'filter_candidates', wraps=validation.filter_candidates) as full:
        result = pipeline.prepare_candidates(a, tmp_path/'run', cfg, allocated_cores=cores,
            generators={'ditmc': learned, 'torsional_diffusion': learned, 'etkdgv3': etkdg},
            mm_optimizer=mm)
    status = json.loads(result.status_path.read_text())
    assert status['threads_per_worker'] == min(4, cores)
    assert status['workers'] == cores // min(4, cores)
    attempted = not (mode.startswith('learned') or mode in ('none', 'mixed_none'))
    assert bool(calls) == attempted
    assert full.call_count == int(attempted)
    assert status['post_mm_validation_mode'] == ('full_recheck' if attempted else 'reused_generation_no_mm')
    assert status['post_mm_validation_wall_seconds'] >= 0
    if mode in ('mm_collision', 'mm_nonfinite', 'mm_connectivity'):
        assert result.state == 'no_valid_candidates_after_mm'
        reason = {'mm_collision': 'collision', 'mm_nonfinite': 'nonfinite_coordinates',
                  'mm_connectivity': 'connectivity'}[mode]
        assert status['post_mm_audit']['rejected'] == {reason: 2}
    elif mode in ('mixed_changed', 'mm_duplicate'):
        assert len(result.candidates) == 1
        assert status['post_mm_audit']['rejected'] == {'duplicate': 1}
    else:
        assert len(result.candidates) == 2


@pytest.mark.parametrize('smiles', ['C[C@H](O)F', 'C/C=C/C'])
def test_mm_stereo_change_is_still_rejected(tmp_path, smiles):
    ref = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(ref, randomSeed=42) == 0
    cfg = SearchConfig.resolve('low', {'budget': {'formula': 'fixed', 'fixed': 1},
        'threads': 1, 'validation': {'minimum_nonbonded_covalent_ratio': 0}})
    def mm(records, *args):
        mol = records[0]
        if '@' in smiles:
            coords = mol.GetConformer().GetPositions()
            coords[:, 0] *= -1
            for i, xyz in enumerate(coords):
                mol.GetConformer().SetAtomPosition(i, xyz)
        else:
            rdMolTransforms.SetDihedralDeg(mol.GetConformer(), 0, 1, 2, 3, 0)
        return [mol], [0]
    result = pipeline.prepare_candidates(ref, tmp_path/'run', cfg, allocated_cores=1,
        generators={'etkdgv3': lambda *a: [Chem.Mol(ref)]}, mm_optimizer=mm)
    status = json.loads(result.status_path.read_text())
    assert result.state == 'no_valid_candidates_after_mm'
    assert status['post_mm_validation_mode'] == 'full_recheck'
    assert status['post_mm_audit']['rejected'] == {'stereo_mismatch': 1}
