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


@pytest.mark.parametrize('backend', ['pm6', 'xtb'])
def test_failed_calculation_keeps_only_bounded_diagnostic_logs(tmp_path, backend):
    folder = candidate_folder(tmp_path)
    (folder / 'status.json').write_text(json.dumps(dict(state='failed', error='native failure')))
    text = 'large prefix\n' * 10000 + 'Error termination: diagnostic reason\n'
    if backend == 'pm6':
        for index in range(2):
            call = folder / f'call_{index:04d}'
            call.mkdir()
            (call / 'input.log').write_text(pm6_log(-10 - index, 0) + (text if index else 'old call\n'))
            (call / 'state.chk').write_bytes(b'large checkpoint')
        (call / 'launcher.out').write_text('launcher diagnostic\n')
    else:
        (folder / 'xtb.out').write_text(text)
        (folder / 'xtbopt.log').write_text('1\n energy: -10.0 gnorm: 0.005\nHe 0.0 0.0 0.0\n')
    compact_relaxation_outputs(folder.parent, backend)
    diagnostic = (folder / 'failure.log').read_text()
    assert 'native failure' in diagnostic and 'diagnostic reason' in diagnostic
    assert 'old call' not in diagnostic
    assert len(diagnostic.encode()) < 66000
    if backend == 'pm6':
        assert 'call_0001/input.log' in diagnostic and 'launcher diagnostic' in diagnostic
    assert {p.name for p in folder.iterdir()} == {'trajectory.jsonl', 'trajectory.xyz', 'status.json', 'failure.log'}


@pytest.mark.parametrize('backend', ['pm6', 'xtb'])
@pytest.mark.parametrize('state', ['converged', 'failed'])
def test_keep_logs_preserves_full_text_without_changing_history(tmp_path, backend, state):
    histories = []
    for keep in (False, True):
        folder = candidate_folder(tmp_path / str(keep))
        (folder / 'status.json').write_text(json.dumps(dict(state=state, error=None)))
        text = 'beginning of native log\n' + 'output\n' * 10000 + 'end of native log\n'
        if backend == 'pm6':
            for index in range(2):
                call = folder / f'call_{index:04d}'
                (call / 'scratch').mkdir(parents=True)
                (call / 'input.log').write_text(pm6_log(-10 - index, 0) + text)
                (call / 'state.chk').write_bytes(b'checkpoint')
                (call / 'state.fchk').write_bytes(b'formatted checkpoint')
            (call / 'launcher.out').write_text('launcher message\n')
            (call / 'formchk.out').write_text('formchk message\n')
        else:
            (folder / 'xtb.out').write_text(text)
            (folder / 'xtbopt.log').write_text('1\n energy: -10.0 gnorm: 0.005\nHe 0.0 0.0 0.0\n')
            (folder / 'xtbrestart').write_bytes(b'restart state')
        compact_relaxation_outputs(folder.parent, backend, keep_logs=keep)
        histories.append(((folder / 'trajectory.jsonl').read_bytes(),
                          (folder / 'trajectory.xyz').read_bytes()))
        if keep:
            saved = (folder / 'native.log').read_text()
            assert text in saved
            if backend == 'pm6':
                assert 'call_0000/input.log' in saved and 'call_0001/input.log' in saved
                assert 'launcher message' in saved and 'formchk message' in saved
            assert {p.name for p in folder.iterdir()} == {
                'trajectory.jsonl', 'trajectory.xyz', 'status.json', 'native.log'}
        else:
            assert not (folder / 'native.log').exists()
            assert (folder / 'failure.log').exists() is (state == 'failed')
    assert histories[0] == histories[1]


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


@pytest.mark.parametrize('problem', [None, 'started_optimization', 'unknown_failure', 'wrong_atom_order'])
def test_xtb_initial_scc_failure_retains_input_without_inventing_results(tmp_path, problem):
    folder = candidate_folder(tmp_path)
    status = dict(state='failed', energy_hartree=None, error='native_failure_or_missing_convergence')
    (folder / 'status.json').write_text(json.dumps(status))
    (folder / '.sccnotconverged').touch()
    text = 'scf: Self consistent charge iterator did not converge\nabnormal termination of xtb\n'
    if problem == 'started_optimization':
        text = 'CYCLE 1\n' + text
    elif problem == 'unknown_failure':
        text = 'unknown native failure\n'
    (folder / 'xtb.out').write_text(text)
    symbol = 'H' if problem == 'wrong_atom_order' else 'He'
    (folder / 'input.xyz').write_text(f'1\ninitial input\n{symbol} 0.123456 0.0 0.0\n')
    compact_relaxation_outputs(folder.parent, 'xtb')
    report = json.loads((folder.parent / 'history_cleanup.json').read_text())
    assert json.loads((folder / 'status.json').read_text()) == status
    if problem:
        assert report['errors']
        assert (folder / 'xtb.out').read_text() == text
        assert not (folder / 'trajectory.xyz').exists()
    else:
        metadata, row = read_history(folder)
        assert not report['errors']
        assert metadata['optimization_started'] is False
        assert metadata['failure_reason'] == 'scc_not_converged_before_optimization'
        assert row['record_type'] == 'initial_structure' and row['step'] == 0
        assert row['energy_hartree'] is None and row['gradient_norm_hartree_per_bohr'] is None
        assert '0.12345600000000' in (folder / 'trajectory.xyz').read_text()
        assert (folder / 'failure.log').read_text().endswith(text + '\n')
        assert {p.name for p in folder.iterdir()} == {'trajectory.jsonl', 'trajectory.xyz', 'status.json', 'failure.log'}
