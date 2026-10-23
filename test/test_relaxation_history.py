"""Protect numerical histories and unrelated files when removing native scratch."""
import json

import pytest
from rdkit import Chem

from qcforever.conformer_search.relaxation_history import compact_relaxation_outputs


def candidate_folder(tmp_path):
    folder = tmp_path / 'electronic/candidate_00000'
    folder.mkdir(parents=True)
    mol = Chem.MolFromSmiles('[He]')
    mol.AddConformer(Chem.Conformer(1))
    mol.SetProp('candidate_id', 'test:0')
    with Chem.SDWriter(str(folder / 'initial.sdf')) as writer:
        writer.write(mol)
    (folder / 'status.json').write_text('{"state":"converged"}')
    (tmp_path / 'DFT.chk').write_text('unrelated checkpoint')
    return folder


def pm6_log(energy, x):
    return f''' Input orientation:
 -----
 Center Atomic Atomic Coordinates (Angstroms)
 Number Number Type X Y Z
 -----
 1 2 0 {x} 0.0 0.0
 -----
 SCF Done: E(RPM6) = {energy} A.U. after 3 cycles
 Forces (Hartrees/Bohr)
 Center Atomic Forces
 Number Number X Y Z
 1 2 0.003 0.004 0.000
 -----
'''


def read_history(folder):
    return [json.loads(line) for line in (folder / 'trajectory.jsonl').read_text().splitlines()]


def test_pm6_history_preserves_every_call_and_removes_native_files(tmp_path):
    folder = candidate_folder(tmp_path)
    for i in range(2):
        call = folder / f'call_{i:04d}'
        (call / 'scratch').mkdir(parents=True)
        (call / 'input.log').write_text(pm6_log(-10 - i, i / 10))
        (call / 'state.chk').write_bytes(b'large checkpoint')
    compact_relaxation_outputs(folder.parent, 'pm6')
    header, a, b = read_history(folder)
    assert header['atomic_numbers'] == [2]
    assert [a['step'], b['step']] == [1, 2]
    assert [a['call'], b['call']] == [0, 1]
    assert [a['energy_hartree'], b['energy_hartree']] == [-10, -11]
    assert a['force_mean_hartree_per_bohr'] == pytest.approx(.005)
    assert b['force_vectors_hartree_per_bohr'] == [[.003, .004, 0.]]
    assert '0.10000000000000' in (folder / 'trajectory.xyz').read_text()
    assert {p.name for p in folder.iterdir()} == {'trajectory.jsonl', 'trajectory.xyz', 'status.json'}
    assert (tmp_path / 'DFT.chk').read_text() == 'unrelated checkpoint'


@pytest.mark.parametrize('truncated', [False, True])
def test_xtb_preserves_xyz_and_distinct_anc_records(tmp_path, truncated):
    folder = candidate_folder(tmp_path)
    xyz = '1\n energy: -10.0 gnorm: 0.005 xtb: 6.6.1\nHe 0.1 0.2 0.3\n'
    (folder / 'xtbopt.log').write_text(xyz + ('1\n energy: -11' if truncated else ''))
    original = dict(cycle=1, energy_hartree=-10., anc_gradient_norm=.004, evaluations=2)
    (folder / 'trajectory.jsonl').write_text(json.dumps(original) + '\n')
    (folder / 'xtb.out').write_text('native output')
    compact_relaxation_outputs(folder.parent, 'xtb')
    header, step, optimizer = read_history(folder)
    assert header['omitted_incomplete_xyz_tail'] is truncated
    assert step['gradient_norm_hartree_per_bohr'] == .005
    assert optimizer == dict(original, record_type='optimizer_stdout')
    assert not (folder / 'xtb.out').exists()


@pytest.mark.parametrize('problem', ['missing_geometry', 'wrong_atom_order'])
def test_unreadable_geometry_preserves_native_files(tmp_path, problem):
    folder = candidate_folder(tmp_path)
    call = folder / 'call_0000'
    call.mkdir()
    (call / 'state.chk').write_bytes(b'preserved after export failure')
    text = pm6_log(-10, 0)
    if problem == 'missing_geometry':
        text = text[text.index(' SCF Done:'):]
    else:
        text = text.replace('1 2 0 0', '1 1 0 0')
    (call / 'input.log').write_text(text)
    compact_relaxation_outputs(folder.parent, 'pm6')
    assert (call / 'input.log').read_text() == text
    assert (call / 'state.chk').read_bytes() == b'preserved after export failure'
    assert not (folder / 'trajectory.xyz').exists()
    assert json.loads((folder.parent / 'history_cleanup.json').read_text())['errors']


def test_candidate_symlinks_are_not_followed(tmp_path):
    folder = candidate_folder(tmp_path)
    other = tmp_path / 'other_job'
    other.mkdir()
    (other / 'state.chk').write_text('keep')
    (folder.parent / 'candidate_00001').symlink_to(other, target_is_directory=True)
    compact_relaxation_outputs(folder.parent, 'pm6')
    assert (other / 'state.chk').read_text() == 'keep'
