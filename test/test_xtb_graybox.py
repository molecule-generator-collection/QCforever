"""Test live POSIX ownership and scheduling; native chemistry is tested on GENKAI."""
import json
import math
import os
from pathlib import Path
import shutil
import signal
import time
from types import SimpleNamespace

import psutil
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from qcforever.conformer_search.settings import parse_conformer_options
from qcforever.conformer_search.run_xtb import XTBCandidate, XTB_FORCE_KIND, cleanup_on_signals
from qcforever.conformer_search.graybox_relaxation import GrayboxSearch, BlockResult
from qcforever.conformer_search.optimize_conformers import optimize_graybox_candidates
from qcforever.util import job_timeout

pytestmark = pytest.mark.skipif(os.name != 'posix', reason='POSIX retained process tests')


@pytest.fixture
def setup_native(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv('QCFOREVER_MODEL_REGISTRY', str(tmp_path/'models.json'))
    binary = tmp_path/'xtb'
    shutil.copyfile(Path(__file__).with_name('fake_xtb_process.py'), binary)
    binary.chmod(0o755)
    mol = Chem.AddHs(Chem.MolFromSmiles('O'))
    AllChem.EmbedMolecule(mol, randomSeed=42)
    mol.SetProp('candidate_id', 'test:0')
    cfg = parse_conformer_options('optconf=xtb', {'relaxation': {
        'xtb_executable': str(binary), 'first_interval': 1, 'subsequent_interval': 2}})
    return mol, cfg


def test_same_process_retained_and_native_final(tmp_path, setup_native):
    mol, cfg = setup_native
    runner = XTBCandidate(mol, tmp_path, 0, 1, '1GB', cfg.relaxation)
    try:
        first = runner.advance()
        pid = runner.process.pid
        assert first.status == 'paused' and runner.stopped
        assert psutil.Process(pid).status() == psutil.STATUS_STOPPED
        assert first.force_kind == XTB_FORCE_KIND
        assert first.force == pytest.approx(.1/math.sqrt(3))
        assert first.energy == -10.01 and first.evaluations == 2
        blocks = [first]
        while not runner.finished:
            blocks.append(runner.advance())
            assert runner.process.pid == pid
        assert blocks[-1].status == 'converged'
        assert sum(b.evaluations for b in blocks) == 7
        assert runner.final_molecule.GetNumAtoms() == 3
        assert runner.process not in job_timeout._processes
        assert runner.process.poll() == 0
    finally:
        runner.close()


def test_close_stopped_process_is_fast_and_idempotent(tmp_path, setup_native):
    mol, cfg = setup_native
    runner = XTBCandidate(mol, tmp_path, 0, 1, '1GB', cfg.relaxation)
    runner.advance()
    started = time.monotonic()
    runner.close()
    runner.close()
    assert time.monotonic()-started < 1.5
    assert runner.process.poll() is not None
    assert runner.process not in job_timeout._processes


def test_native_failure_not_converged(tmp_path, setup_native, monkeypatch):
    monkeypatch.setenv('QCFOREVER_FAKE_XTB_MODE', 'fail')
    mol, cfg = setup_native
    runner = XTBCandidate(mol, tmp_path, 0, 1, '1GB', cfg.relaxation)
    try:
        blocks = [runner.advance()]
        while not runner.finished:
            blocks.append(runner.advance())
        assert blocks[-1].status == 'failed'
        assert runner.final_molecule is None
        assert sum(b.evaluations for b in blocks) == 7
    finally:
        runner.close()


def test_close_handles_permission_error_only_after_exit(tmp_path, setup_native, monkeypatch):
    mol, cfg = setup_native
    runner = XTBCandidate(mol, tmp_path, 0, 1, '1GB', cfg.relaxation)
    runner.advance()
    killpg = os.killpg
    def disappearing_group(pid, sig):
        if sig == signal.SIGCONT:
            killpg(pid, signal.SIGKILL)
            runner.process.wait(timeout=2)
            raise PermissionError('simulated exited process group')
        killpg(pid, sig)
    monkeypatch.setattr(os, 'killpg', disappearing_group)
    runner.close()
    assert runner.process.poll() is not None
    assert runner.process not in job_timeout._processes


def test_close_does_not_hide_permission_error_on_live_child(tmp_path, setup_native, monkeypatch):
    mol, cfg = setup_native
    runner = XTBCandidate(mol, tmp_path, 0, 1, '1GB', cfg.relaxation)
    runner.advance()
    killpg = os.killpg
    def denied_group(pid, sig):
        raise PermissionError('simulated live process permission failure')
    try:
        monkeypatch.setattr(os, 'killpg', denied_group)
        with pytest.raises(PermissionError, match='live process'):
            runner.close()
        assert runner.process.poll() is None
    finally:
        monkeypatch.setattr(os, 'killpg', killpg)
        runner.close()


def test_active_timeout_and_overall_timeout_cleanup(tmp_path, setup_native, monkeypatch):
    monkeypatch.setenv('QCFOREVER_FAKE_XTB_MODE', 'hang')
    mol, cfg = setup_native
    settings = dict(cfg.relaxation, first_interval=10, xtb_active_timeout_seconds=1)
    runner = XTBCandidate(mol, tmp_path, 0, 1, '1GB', settings)
    block = runner.advance()
    assert block.status == 'failed' and block.evaluations == 2
    assert 'active-time' in block.error and runner.process.poll() is not None
    folder = tmp_path/'overall'
    folder.mkdir()
    runner = XTBCandidate(mol, folder, 0, 1, '1GB', settings)
    with pytest.raises(job_timeout.QCforeverTimeoutError):
        with job_timeout.overall_timeout(.3):
            runner.advance()
    assert runner.process.poll() is not None
    assert runner.process not in job_timeout._processes


def test_force_proxy_explicit_and_not_mixed(setup_native):
    mol, cfg = setup_native
    settings = dict(cfg.relaxation, algorithm='laqa')
    def advance(key):
        return BlockResult(-10, .1, 2, status='converged', force_kind=XTB_FORCE_KIND)
    good = GrayboxSearch({'a': 3}, advance, settings, force_kind=XTB_FORCE_KIND).run()
    assert good.done()
    bad = GrayboxSearch({'a': 3}, advance, settings).run()
    assert not bad.done() and bad.states['a'].status == 'failed'


@pytest.mark.parametrize('memory_failure', [False, True])
def test_workflow_cleans_entire_pool(tmp_path, setup_native, monkeypatch, memory_failure):
    mol, cfg = setup_native
    mols = [Chem.Mol(mol) for _ in range(5)]
    for i,m in enumerate(mols):
        m.SetProp('candidate_id', str(i))
    prepared = SimpleNamespace(status_path=tmp_path/'status.json', candidates=mols)
    runners = []
    def factory(*args):
        runner = XTBCandidate(*args)
        runners.append(runner)
        return runner
    if memory_failure:
        monkeypatch.setattr('qcforever.conformer_search.optimize_conformers.pool_memory_mb', lambda _: 9000)
        with pytest.raises(MemoryError):
            optimize_graybox_candidates(prepared, mol, cfg, 0, 1, 32, '1GB', method='xtb', candidate_factory=factory)
    else:
        audit = optimize_graybox_candidates(prepared, mol, cfg, 0, 1, 32, '1GB', method='xtb', candidate_factory=factory)
        assert audit['backend'] == 'xtb' and audit['quota_reached']
        assert audit['score_definition'] == XTB_FORCE_KIND and audit['parallel_workers'] == 5
    assert all(r.process is None or r.process.poll() is not None for r in runners)
    assert (tmp_path/'electronic/search_state.json').exists()


def test_signal_handler_restored():
    previous = signal.getsignal(signal.SIGTERM)
    with pytest.raises(InterruptedError):
        with cleanup_on_signals():
            os.kill(os.getpid(), signal.SIGTERM)
    assert signal.getsignal(signal.SIGTERM) == previous


def test_cleanup_error_does_not_skip_remaining_candidates(tmp_path, setup_native, monkeypatch):
    mol, cfg = setup_native
    prepared = SimpleNamespace(status_path=tmp_path/'status.json', candidates=[mol,mol,mol])
    closed = []
    class Candidate:
        def __init__(self, mol, folder, *args):
            self.name = folder.name
        def advance(self):
            return BlockResult(-10,.1,2,force_kind=XTB_FORCE_KIND)
        def close(self):
            closed.append(self.name)
            if self.name.endswith('00000'):
                raise OSError('simulated cleanup error')
    monkeypatch.setattr('qcforever.conformer_search.optimize_conformers.pool_memory_mb',lambda _:9000)
    with pytest.raises(RuntimeError,match='could not be cleaned up'):
        optimize_graybox_candidates(prepared,mol,cfg,0,1,1,'1GB',method='xtb',candidate_factory=Candidate)
    assert len(closed)==3
    assert (tmp_path/'electronic/cleanup_errors.json').exists()
