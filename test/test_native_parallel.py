"""Spawned worker contracts without installed Gaussian or xTB binaries."""
import json
import os
from pathlib import Path
import subprocess
import time
from unittest.mock import patch

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from qcforever.conformer_search.settings import SearchConfig
from qcforever.conformer_search.generate_conformers import prepare_candidates
from qcforever.conformer_search.optimize_conformers import (
    native_thread_environment, optimize_with_pm6, optimize_candidates, _initialize_native_worker,
)
from qcforever.conformer_search.run_pm6 import PM6Candidate


def test_native_worker_affinity(monkeypatch):
    from queue import Queue
    calls = []
    monkeypatch.setattr(os, 'sched_setaffinity', lambda pid, cpus: calls.append((pid, cpus)), raising=False)
    queue = Queue()
    queue.put([8])
    queue.put([14])
    _initialize_native_worker(queue)
    _initialize_native_worker(queue)
    _initialize_native_worker(None)
    assert calls == [(0, [8]), (0, [14])]


def _copies(reference, count, *args):
    return [Chem.Mol(reference) for _ in range(count)]


def _pm6_adapter(mol, folder, charge, multiplicity, cores, memory, settings):
    """Exercise continuous RFO input, parsing, and isolation in spawned workers."""
    before_cwd, before_env = Path.cwd(), dict(os.environ)
    index = int(folder.name.split('_')[-1])
    assert cores == 1 and settings['pm6_optimizer'] == 'rfo'

    def run(command, cwd, env, **kwargs):
        assert Path.cwd() == before_cwd
        text = (cwd / 'input.com').read_text()
        assert '%nprocshared=1' in text and '%mem=1GB' in text
        assert 'PM6 NoSymm SCF=(Tight,MaxCycle=512)' in text
        assert 'Opt=(RFO,Tight,MaxCycles=1000)' in text and 'Restart' not in text
        assert all(env[key] == value for key, value in native_thread_environment(1).items())
        assert env['GAUSS_SCRDIR'] == str(cwd / 'scratch')
        start = time.monotonic()
        time.sleep(0.4 + (index % 3) * 0.1)
        (folder / 'interval.json').write_text(json.dumps({
            'start': start, 'end': time.monotonic(), 'pid': os.getpid()}))
        if index == 1:
            log = 'SCF has not converged\nError termination\n'
        else:
            log = f'SCF Done: E(RPM6) = {-10.0-index} A.U. after 3 cycles\n'
            log += ' Forces (Hartrees/Bohr)\n Center Atomic Forces\n Number Number X Y Z\n'
            log += ''.join(f'{i+1} {atom.GetAtomicNum()} 0.0 0.0 0.0\n'
                           for i, atom in enumerate(mol.GetAtoms()))
            log += ' ----\nStationary point found\nNormal termination of Gaussian\n'
        (cwd / 'input.log').write_text(log)
        return subprocess.CompletedProcess(command, 1 if index == 1 else 0)

    try:
        with patch('shutil.which', side_effect=lambda name: '/fake/' + name), \
             patch('qcforever.conformer_search.run_pm6.job_timeout.run', side_effect=run), \
             patch.object(PM6Candidate, '_read_final_coordinates',
                          lambda candidate, call: candidate.mol.GetConformer().GetPositions()):
            return optimize_with_pm6(mol, folder, charge, multiplicity, cores, memory, settings)
    finally:
        assert Path.cwd() == before_cwd
        assert dict(os.environ) == before_env


@pytest.mark.parametrize('nproc', [1, 2])
def test_pm6_parallel_isolation_order_and_failure(tmp_path, monkeypatch, nproc):
    monkeypatch.chdir(tmp_path)
    ref = Chem.AddHs(Chem.MolFromSmiles('CCCC'))
    assert AllChem.EmbedMolecule(ref, randomSeed=42) == 0
    cfg = SearchConfig.resolve('low', {
        'relaxation': {'implementation': 'continuous'},
        'budget': {'formula': 'fixed', 'fixed': 10}, 'mm_method': 'none', 'threads': 1,
        'validation': {'duplicate_rmsd_angstrom': 0}})
    prepared = prepare_candidates(ref, tmp_path/'run', cfg, allocated_cores=nproc,
                                  generators={'etkdgv3': _copies})
    audit = optimize_candidates(prepared, ref, cfg, 0, 1, 'pm6', nproc, '1GB', adapter=_pm6_adapter)
    workers = min(nproc, 10)
    assert audit['parallel_workers'] == workers
    assert audit['cores_per_calculation'] == 1
    assert audit['failed_candidates'] == 1 and audit['converged_candidates'] == 9
    assert audit['selected_index'] == 9
    assert [r['prepared_index'] for r in audit['candidate_runs']] == list(range(10))
    progress = json.loads((tmp_path/'run/electronic/progress.json').read_text())
    assert progress == audit['candidate_runs']
    blocks = [json.loads(path.read_text()) for path in (tmp_path/'run/electronic').glob('*/call_*/result.json')]
    assert len(blocks) == 10
    assert sum(block['evaluations'] for block in blocks) == 9
    assert sum(block['scf_cycles'] for block in blocks) == 27
    intervals = [json.loads(p.read_text()) for p in (tmp_path/'run/electronic').glob('*/interval.json')]
    events = sorted([(r['start'], 1) for r in intervals] + [(r['end'], -1) for r in intervals])
    active = maximum = 0
    for _, delta in events:
        active += delta
        maximum = max(maximum, active)
    assert (1 if workers == 1 else 2) <= maximum <= workers
    assert len({r['pid'] for r in intervals}) <= workers
