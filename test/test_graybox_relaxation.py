"""Graybox policies, parallel resource ownership, and PM6 restart contracts.

Native calculations are replaced by saved blocks, synthetic logs, or sleeping
Python children. Option syntax is covered separately in test_conformer_search.
"""
from dataclasses import replace
import json
import os
import sys
import threading
import time
from types import SimpleNamespace

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from qcforever.conformer_search.settings import GRAYBOX_DEFAULTS, SearchConfig, resolve_relaxation
from qcforever.conformer_search.graybox_relaxation import BlockResult, GrayboxSearch
from qcforever.conformer_search.calculation_logs import read_pm6_block
from qcforever.conformer_search.run_pm6 import PM6Candidate
from qcforever.conformer_search.optimize_conformers import optimize_graybox_candidates
from qcforever.conformer_search.run_candidate_blocks import CandidateBlockExecutor, resolve_block_workers
from qcforever.util import job_timeout


def replay_blocks(blocks, workers=1, **settings):
    """Use the same saved trajectories for serial and reversed-arrival batches."""
    cursors = dict.fromkeys(blocks, 0)
    batches = []
    def advance_batch(keys):
        batches.append(list(keys))
        result = {}
        for key in reversed(keys):
            result[key] = replace(blocks[key][cursors[key]])
            cursors[key] += 1
        return result
    search = GrayboxSearch(dict.fromkeys(blocks, 10),
        lambda key: advance_batch([key])[key], dict(GRAYBOX_DEFAULTS, **settings),
        workers=workers, advance_batch=advance_batch if workers > 1 else None)
    return search.run(), batches


@pytest.mark.parametrize('workers', [1, 2])
def test_fixed_denominator_initialization_and_failure_cost(workers):
    blocks = {'a': [BlockResult(-1, .1, 1, status='failed', wall_seconds=2)],
              'b': [BlockResult(-2, .1, 1, status='converged', wall_seconds=3)],
              'c': [BlockResult(-3, .1, 1), BlockResult(-4, .01, 10, status='converged')],
              'd': [BlockResult(-2, .1, 1), BlockResult(-3, .01, 10, status='converged')]}
    search, _ = replay_blocks(blocks, workers, algorithm='sh', convergence_fraction=.5)
    assert search.required == 2 and search.done()
    assert search.cost['wall_seconds'] == 5
    assert search.events[-1]['converged_count'] == 1 + workers
    assert search.cost['evaluations'] == 4 + 10 * workers
    assert search.states['a'].status == 'failed'
    assert search.initial_count == 4
    full, _ = replay_blocks(blocks, workers, algorithm='sh', convergence_fraction=1)
    assert not full.done() and full.stop_reason == 'candidate_pool_exhausted'
    assert full.events[-1]['converged_count'] == 3


@pytest.mark.parametrize('workers', [1, 4])
def test_quota_charges_dispatched_batches_without_seeing_future(workers):
    blocks = {str(i): [BlockResult(-i, .1, 1, wall_seconds=2, status='converged')] for i in range(10)}
    search, batches = replay_blocks(blocks, workers, algorithm='sh', convergence_fraction=.2)
    completed = 2 if workers == 1 else 4
    assert search.required == 2 and search.done()
    assert search.cost['evaluations'] == completed
    assert search.cost['wall_seconds'] == 2 * completed
    assert search.states[str(completed)].status == 'pending'
    assert [event['candidate'] for event in search.events] == [str(i) for i in range(completed)]
    assert batches == ([['0'], ['1']] if workers == 1 else [['0', '1', '2', '3']])
    if workers > 1:
        assert all(event['batch_size'] == workers for event in search.events)
    blocks['9'] = [BlockResult(-1e9, .1, 1, status='converged')]
    changed, _ = replay_blocks(blocks, workers, algorithm='sh', convergence_fraction=.2)
    assert changed.events == search.events


@pytest.mark.parametrize('algorithm,workers', [('sh', 1), ('laqa', 2), ('sr', 2), ('sh', 2)])
def test_budget_extension_charges_whole_overshooting_blocks(algorithm, workers):
    blocks = {str(i): [BlockResult(-i, .1, 1)] + [BlockResult(-i, .1, 10) for _ in range(50)] for i in range(4)}
    search, _ = replay_blocks(blocks, workers, algorithm=algorithm, maximum_budget_rounds=3,
                             initial_steps_per_candidate=2, budget_increment_steps_per_candidate=2)
    assert search.stop_reason == 'global_budget_cap'
    assert 24 <= search.cost['evaluations'] < 24 + 10 * workers
    assert search.round == 3 and search.cost['calls'] == len(search.events)
    assert all(row.status == 'paused' for row in search.states.values())


def test_laqa_force_definition_and_missing_force():
    search, _ = replay_blocks({'a': [BlockResult(-10, None, 1)]}, algorithm='laqa', convergence_fraction=1)
    assert search.states['a'].status == 'failed'
    assert search.cost['evaluations'] == 1
    search = GrayboxSearch({'a': 10}, None, dict(GRAYBOX_DEFAULTS, algorithm='laqa'))
    state = search.states['a']
    state.energy, state.forces = -10, [2]
    assert search.rank('a', 'laqa') == (-3, 'a')


@pytest.mark.parametrize('algorithm', ['sr', 'sh'])
def test_elimination_waits_for_equal_stage_targets(algorithm):
    blocks = {str(i): [BlockResult(-i, .1, 1) for _ in range(101)] for i in range(6)}
    search, batches = replay_blocks(blocks, 3, algorithm=algorithm, convergence_fraction=1,
                                   maximum_budget_rounds=1, initial_steps_per_candidate=8)
    stage = [event for event in search.events if event['phase'].endswith('selection_0')]
    counts = {key: sum(event['candidate'] == key for event in stage) for key in blocks}
    assert min(counts.values()) > 0 and len(set(counts.values())) == 1
    assert all(len(keys) <= 3 and len(keys) == len(set(keys)) for keys in batches)


@pytest.mark.parametrize('algorithm', ['laqa', 'sh', 'sr'])
def test_hundred_candidates_are_selected_globally_not_in_worker_groups(algorithm):
    blocks = {f'{i:03d}': [BlockResult(-i, .1, 1) for _ in range(1000)] for i in range(100)}
    search, batches = replay_blocks(blocks, 32, algorithm=algorithm, convergence_fraction=1,
                                   maximum_budget_rounds=1, initial_steps_per_candidate=50)
    assert [event['candidate'] for event in search.events[:100]] == list(blocks)
    assert all(len(batch) <= 32 and len(batch) == len(set(batch)) for batch in batches)
    if algorithm == 'laqa':
        assert batches[4] == [f'{i:03d}' for i in range(99, 67, -1)]
        return
    stages = sorted({int(event['phase'].split('_')[-1]) for event in search.events
                     if event['phase'].startswith(f'{algorithm.upper()}_selection_')})
    members = {stage: {event['candidate'] for event in search.events
                      if event['phase'] == f'{algorithm.upper()}_selection_{stage}'} for stage in stages[:2]}
    assert members[0] == set(blocks)
    next_stage = stages[1]  # SR rounding can skip a stage without new evaluations.
    survivors = 50 if algorithm == 'sh' else 100 - next_stage
    assert members[next_stage] == {f'{i:03d}' for i in range(100 - survivors, 100)}


def test_worker_limits_respect_cpu_allocation(monkeypatch):
    monkeypatch.setattr(os, 'sched_getaffinity', lambda _: {4, 5, 6, 7}, raising=False)
    assert resolve_block_workers(32, 100) == 4
    assert resolve_block_workers(32, 2) == 2
    assert resolve_block_workers(32, 100, 1) == 1
    assert resolve_block_workers(2, 100, 8) == 2


def test_executor_overlaps_without_reordering_or_sharing_cpu_slots(monkeypatch):
    monkeypatch.setattr(os, 'sched_getaffinity', lambda _: {4, 5, 6}, raising=False)
    barrier = threading.Barrier(3)
    observed = []
    class Candidate:
        def __init__(self, key):
            self.key = key
        def advance(self):
            observed.append((self.key, self.cpu_id))
            barrier.wait(timeout=3)
            return BlockResult(-int(self.key), .1, 1)
    runners = {str(i): Candidate(str(i)) for i in range(3)}
    with CandidateBlockExecutor(runners, 3) as pool:
        result = pool.advance_batch(['2', '0', '1'])
    assert list(result) == ['2', '0', '1']
    assert set(observed) == {('2', 4), ('0', 5), ('1', 6)}


def test_parallel_deadline_terminates_owned_children():
    children = []
    class Candidate:
        def advance(self):
            child = job_timeout.popen([sys.executable, '-c', 'import time; time.sleep(30)'])
            children.append(child)
            job_timeout.wait(child)
            return BlockResult(status='failed')
    started = time.monotonic()
    with pytest.raises(job_timeout.QCforeverTimeoutError):
        with job_timeout.overall_timeout(.5):
            with CandidateBlockExecutor({'a': Candidate(), 'b': Candidate()}, 2) as pool:
                pool.advance_batch(['a', 'b'])
    assert children and all(child.poll() is not None for child in children)
    assert time.monotonic() - started < 5
    assert not any(child in job_timeout._processes for child in children)


def test_executor_error_does_not_kill_unrelated_children():
    other = job_timeout.popen([sys.executable, '-c', 'import time; time.sleep(30)'])
    class Candidate:
        def advance(self):
            raise OSError('test worker error')
    try:
        with pytest.raises(OSError, match='worker error'):
            with CandidateBlockExecutor({'a': Candidate(), 'b': Candidate()}, 2) as pool:
                pool.advance_batch(['a', 'b'])
        assert other.poll() is None
    finally:
        job_timeout.terminate_process(other)
        job_timeout.unregister_process(other)


def molecule():
    mol = Chem.AddHs(Chem.MolFromSmiles('O'))
    AllChem.EmbedMolecule(mol, randomSeed=42)
    mol.SetProp('candidate_id', 'test:0')
    return mol


def gaussian_log(n=1, converged=False, atom_count=3):
    force = ' Forces (Hartrees/Bohr)\n Center Atomic Forces\n Number Number X Y Z\n'
    force += ''.join(f' {i+1} 1 0.010000 0.000000 0.000000\n' for i in range(atom_count)) + ' ----\n'
    record = ' SCF Done: E(RPM6) = -10.000000 A.U. after 3 cycles\n' + force
    ending = 'Stationary point found\nNormal termination\n' if converged else 'Number of steps exceeded\nError termination\n'
    return record * n + ending


def test_log_force_mean_and_scf_failure():
    text = gaussian_log(atom_count=2).replace('0.010000 0.000000 0.000000', '3.0 4.0 0.0', 1)
    parsed = read_pm6_block(text, 2)
    assert parsed['records'][0]['force_mean_hartree_per_bohr'] == pytest.approx(2.505)
    assert not parsed['native_converged'] and parsed['cycle_limit']
    assert read_pm6_block(text, 3)['records'][0]['force_mean_hartree_per_bohr'] is None


def test_nonfinite_log_still_counts_evaluation():
    text = gaussian_log().replace('-10.000000', 'NaN').replace('0.010000', 'NaN')
    parsed = read_pm6_block(text, 3)
    assert len(parsed['records']) == 1
    assert parsed['records'][0]['scf_cycles'] == 3
    assert parsed['records'][0]['energy_hartree'] is None
    assert parsed['records'][0]['force_mean_hartree_per_bohr'] is None


def test_call_timeout_retains_partial_cost_and_previous_checkpoint(tmp_path, monkeypatch):
    import subprocess
    monkeypatch.setattr('shutil.which', lambda name: '/fake/' + name)
    runner = PM6Candidate(molecule(), tmp_path, 0, 1, '1GB', dict(GRAYBOX_DEFAULTS, maximum_cycles=100))
    old = tmp_path / 'last_good.chk'
    old.write_text('good checkpoint')
    runner.checkpoint, runner.target = old, 1
    def timeout(command, cwd, **kwargs):
        (cwd / 'input.log').write_text(gaussian_log(2).replace('Number of steps exceeded', 'interrupted'))
        raise subprocess.TimeoutExpired(command, 1)
    monkeypatch.setattr('qcforever.conformer_search.run_pm6.job_timeout.run', timeout)
    block = runner.advance()
    assert block.status == 'failed' and block.evaluations == 2 and block.scf_cycles == 6
    assert 'TimeoutExpired' in block.error and block.wall_seconds > 0
    assert runner.checkpoint == old and old.read_text() == 'good checkpoint'


def test_unexpected_short_block_is_not_accepted(tmp_path, monkeypatch):
    monkeypatch.setattr('shutil.which', lambda name: '/fake/' + name)
    runner = PM6Candidate(molecule(), tmp_path, 0, 1, '1GB', dict(GRAYBOX_DEFAULTS, first_interval=10))
    def short_call(folder):
        (folder / 'input.log').write_text(gaussian_log(2))
    monkeypatch.setattr(runner, '_run_call', short_call)
    block = runner.advance()
    assert block.status == 'failed' and block.evaluations == 2
    assert 'Unexpected evaluation count' in block.error


@pytest.mark.parametrize('previous_target,extra,records,printed_maximum', [(0, 1000, 3, 3), (1, 10, 2, 3)])
def test_gaussian_internal_ceiling_is_limit_not_parse_failure(tmp_path, monkeypatch, previous_target, extra, records, printed_maximum):
    monkeypatch.setattr('shutil.which', lambda name: '/fake/' + name)
    mol = molecule()
    runner = PM6Candidate(mol, tmp_path, 0, 1, '1GB', dict(GRAYBOX_DEFAULTS,
                           first_interval=extra, subsequent_interval=extra, maximum_cycles=1000))
    runner.target = previous_target
    if previous_target:
        runner.checkpoint = tmp_path / 'previous.chk'
        runner.checkpoint.write_text('previous')
    def native_limit(folder):
        text = gaussian_log(records) + f'\nStep number {previous_target + records} out of a maximum of {printed_maximum}\n'
        (folder / 'input.log').write_text(text)
    monkeypatch.setattr(runner, '_run_call', native_limit)
    monkeypatch.setattr(runner, '_read_final_coordinates', lambda folder: mol.GetConformer().GetPositions())
    block = runner.advance()
    assert block.status == 'limit' and block.error == 'native_iteration_limit'
    assert block.evaluations == records
    assert runner.final_molecule is None


def test_gaussian_step_storage_error_remains_failure(tmp_path, monkeypatch):
    monkeypatch.setattr('shutil.which', lambda name: '/fake/' + name)
    runner = PM6Candidate(molecule(), tmp_path, 0, 1, '1GB', dict(GRAYBOX_DEFAULTS, maximum_cycles=1000))
    runner.target = 291
    runner.checkpoint = tmp_path / 'previous.chk'
    runner.checkpoint.write_text('previous')
    def storage_error(folder):
        text = gaussian_log(4).replace('Number of steps exceeded',
            'Step requested is out of bounds\nIStep = 295 but MinStp = 1 and MaxStp = 294\nRdWrOT: IStep out of bounds')
        (folder / 'input.log').write_text(text)
    monkeypatch.setattr(runner, '_run_call', storage_error)
    block = runner.advance()
    assert block.status == 'failed' and block.evaluations == 4
    assert runner.checkpoint == tmp_path / 'previous.chk'
    assert runner.final_molecule is None


def test_native_restart_input_and_precision(tmp_path, monkeypatch):
    monkeypatch.setattr('shutil.which', lambda name: '/fake/' + name)
    mol = molecule()
    commands = []
    def run(command, cwd, **kwargs):
        commands.append(command)
        if command[0].endswith('g16'):
            text = (cwd / 'input.com').read_text()
            assert '%nprocshared=1' in text and 'RFO,Tight' in text
            if len(commands) == 1:
                assert 'Restart' not in text and 'MaxCycles=1' in text
                assert f'{0} {1}\n' in text
            else:
                assert 'Restart,RFO,Tight,MaxCycles=11' in text
                assert (cwd / 'state.chk').read_text() == 'checkpoint'
            (cwd / 'state.chk').write_text('checkpoint')
            (cwd / 'input.log').write_text(gaussian_log(1 if len(commands) == 1 else 3, len(commands) != 1))
        else:
            xyz = ' '.join(str(v / .529177210903) for row in mol.GetConformer().GetPositions() for v in row)
            (cwd / 'state.fchk').write_text('Atomic numbers I N= 3\n8 1 1\nCurrent cartesian coordinates R N= 9\n' + xyz)
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr('qcforever.conformer_search.run_pm6.job_timeout.run', run)
    runner = PM6Candidate(mol, tmp_path, 0, 1, '1GB', dict(GRAYBOX_DEFAULTS, maximum_cycles=100))
    first, second = runner.advance(), runner.advance()
    assert first.status == 'paused' and second.status == 'converged'
    assert first.evaluations == 1 and second.evaluations == 3
    assert runner.final_molecule.GetConformer().GetPositions() == pytest.approx(mol.GetConformer().GetPositions())
    with pytest.raises(RuntimeError):
        runner.advance()


def test_graybox_handoff_only_converged(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    cfg = SearchConfig(relaxation=resolve_relaxation(
        {'convergence_fraction': .2, 'parallel_candidates': 1}, 'pm6'))
    mols = [molecule() for _ in range(10)]
    for i, mol in enumerate(mols):
        mol.SetProp('candidate_id', str(i))
    prepared = SimpleNamespace(status_path=tmp_path / 'status.json', candidates=mols)
    class FakeCandidate:
        def __init__(self, mol, *args):
            self.final_molecule = mol
        def advance(self):
            return BlockResult(-1, .1, 1, status='converged')
    audit = optimize_graybox_candidates(prepared, mols[0], cfg, 0, 1, 32, '1GB', candidate_factory=FakeCandidate)
    assert audit['attempted_candidates'] == 2 and audit['converged_candidates'] == 2
    assert audit['unfinished_candidates'] == 8 and audit['failed_candidates'] == 0
    assert audit['parallel_workers'] == 1 and audit['allocated_cores'] == 32
    assert (tmp_path / 'optimized_structures.sdf').exists()
    assert len(list(Chem.SDMolSupplier(str(tmp_path / 'electronic/all_converged.sdf'), removeHs=False))) == 2


def test_summary_warns_when_quota_is_unreachable(tmp_path, monkeypatch):
    from qcforever.conformer_search.conformer_search import _save_search_summary
    monkeypatch.chdir(tmp_path)
    cfg = SearchConfig(relaxation=resolve_relaxation({'convergence_fraction': 1}, 'pm6'))
    mols = [molecule(), molecule()]
    status = tmp_path / 'status.json'
    status.write_text(json.dumps(dict(budget={'maximum_candidates': 2}, mm_state='skipped',
                                     mm_wall_seconds=0, stages=[])))
    prepared = SimpleNamespace(status_path=status, candidates=mols)
    class PartlyFailed:
        def __init__(self, mol, folder, *args):
            self.final_molecule = mol
            self.failed = folder.name.endswith('00001')
        def advance(self):
            return BlockResult(-1, .1, 1, status='failed' if self.failed else 'converged')
    audit = optimize_graybox_candidates(prepared, mols[0], cfg, 0, 1, 1, '1GB', candidate_factory=PartlyFailed)
    summary = _save_search_summary(tmp_path, prepared, cfg, 'pm6', audit)
    assert summary['state'] == 'succeeded_with_search_warning'
    assert summary['search_warning'] == 'convergence_target_not_reached'
    assert summary['convergence_target_reached'] is False
    assert summary['converged_candidates'] == 1 and summary['failed_candidates'] == 1
