"""Incremental filtering must preserve the full-filter decisions and audit."""
from dataclasses import replace
from unittest.mock import patch

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms

from qcforever.conformer_search import validation
from qcforever.conformer_search.config import SearchConfig
from qcforever.conformer_search.pipeline import prepare_candidates


def rotamers():
    mol = Chem.AddHs(Chem.MolFromSmiles('CCO'))
    assert AllChem.EmbedMolecule(mol, randomSeed=42) == 0
    atoms = mol.GetSubstructMatch(Chem.MolFromSmarts('[C]-[C]-[O]-[H]'))
    result = []
    for i, angle in enumerate((180, 60, -60)):
        copy = Chem.Mol(mol)
        rdMolTransforms.SetDihedralDeg(copy.GetConformer(), *atoms, angle)
        copy.SetProp('candidate_id', str(i))
        result.append(copy)
    return result


def assert_same(left, right):
    assert len(left) == len(right)
    for a, b in zip(left, right):
        assert a.GetPropsAsDict() == b.GetPropsAsDict()
        assert Chem.MolToSmiles(a) == Chem.MolToSmiles(b)
        np.testing.assert_array_equal(a.GetConformer().GetPositions(), b.GetConformer().GetPositions())


@pytest.mark.parametrize('maximum', [1, 2, 10])
@pytest.mark.parametrize('threshold', [0, 0.1])
def test_batches_match_full_filter(maximum, threshold):
    a, b, c = rotamers()
    collision = Chem.Mol(a)
    collision.GetConformer().SetAtomPosition(1, collision.GetConformer().GetAtomPosition(0))
    cfg = replace(SearchConfig.resolve().validation, duplicate_rmsd_angstrom=threshold)
    incremental = validation.IncrementalCandidateFilter(a, maximum, cfg)
    full = []
    for batch in ([None, a, a], [], [b, b, collision], [c, a]):
        full, old_audit = validation.filter_candidates(full + batch, a, maximum, cfg)
        new, new_audit = incremental.extend(batch)
        assert old_audit == new_audit
        assert_same(full, new)


@pytest.mark.parametrize('smiles', ['C[C@H](O)F', 'C/C=C/C'])
def test_stereo_rejection_matches_full_filter(smiles):
    ref = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(ref, randomSeed=42) == 0
    wrong = Chem.Mol(ref)
    if '@' in smiles:
        coords = wrong.GetConformer().GetPositions()
        coords[:, 0] *= -1
        for i, xyz in enumerate(coords):
            wrong.GetConformer().SetAtomPosition(i, xyz)
    else:
        rdMolTransforms.SetDihedralDeg(wrong.GetConformer(), 0, 1, 2, 3, 0)
    cfg = replace(SearchConfig.resolve().validation, minimum_nonbonded_covalent_ratio=0)
    incremental = validation.IncrementalCandidateFilter(ref, 10, cfg)
    kept, _ = incremental.extend([ref])
    full, old = validation.filter_candidates(kept + [wrong], ref, 10, cfg)
    new, audit = incremental.extend([wrong])
    assert old == audit
    assert audit['rejected'] == {'stereo_mismatch': 1}
    assert_same(full, new)


def test_only_new_records_are_checked_and_copies_are_private():
    a, b, c = rotamers()
    cfg = SearchConfig.resolve().validation
    incremental = validation.IncrementalCandidateFilter(a, 10, cfg)
    accepted, audit = incremental.extend([a, b])
    # Neither returned objects nor the original input may corrupt the cache.
    for mol in (a, accepted[0]):
        mol.GetConformer().SetAtomPosition(0, (100, 100, 100))
    audit['candidates'][0]['reason'] = 'corrupted_by_caller'
    with patch.object(validation, 'geometry_failure', wraps=validation.geometry_failure) as geometry, \
         patch.object(validation, 'stereo_audit', wraps=validation.stereo_audit) as stereo, \
         patch.object(validation, 'duplicate_rmsd', wraps=validation.duplicate_rmsd) as rmsd:
        kept, new_audit = incremental.extend([c])
    assert len(kept) == 3
    assert geometry.call_count == stereo.call_count == 1
    assert rmsd.call_count == 2
    assert new_audit['candidates'][0]['reason'] is None


def test_pipeline_rechecks_changed_coordinates_after_mm(tmp_path):
    a, b, _ = rotamers()
    cfg = SearchConfig.resolve('low', {
        'budget': {'formula': 'fixed', 'fixed': 2}, 'mm_method': 'mmff94s'})

    def generate(*args):
        return [Chem.Mol(a), Chem.Mol(b)]

    def same_coordinates(records, *args):
        return [validation.copy_with_coordinates(records[0], a)], [0]

    result = prepare_candidates(a, tmp_path/'search', cfg, allocated_cores=4,
                                generators={'etkdgv3': generate}, mm_optimizer=same_coordinates)
    assert len(result.candidates) == 1
    import json
    status = json.loads(result.status_path.read_text())
    assert status['stages'][0]['batches'][0]['accumulated_valid'] == 2
    assert status['post_mm_audit']['rejected'] == {'duplicate': 1}
