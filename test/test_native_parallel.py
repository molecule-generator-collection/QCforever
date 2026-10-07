"""Spawned worker contracts without installed Gaussian or xTB binaries."""
import json
import os
from pathlib import Path
import time
from unittest.mock import patch

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from qcforever.conformer_search.config import SearchConfig
from qcforever.conformer_search.pipeline import prepare_candidates
from qcforever.conformer_search.relaxation import (
    native_thread_environment, pm6_relax, relax_candidates, _initialize_native_worker,
)


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
    """Exercise the real PM6 input/environment adapter inside spawned workers."""
    from qcforever.laqa_fafoom.pyg16 import g16Object
    before_cwd, before_env = Path.cwd(), dict(os.environ)
    index = int(folder.name.split('_')[-1])

    def run(obj):
        assert Path.cwd() == folder
        text = Path('Gau_molecule.com').read_text()
        assert cores == 1
        assert '%nprocshared=1' in text and '%mem=1GB' in text
        assert all(os.environ[k] == v for k, v in native_thread_environment(1).items())
        start = time.monotonic()
        time.sleep(0.4 + (index % 3) * 0.1)
        Path('interval.json').write_text(json.dumps({
            'start': start, 'end': time.monotonic(), 'pid': os.getpid()}))
        if index == 1:
            Path('Gau_molecule.log').write_text('SCF failed\n')
            raise RuntimeError('synthetic SCF failure')
        obj.energy = -10 - index
        obj.sdf_string_opt = obj.sdf_string
        Path('Gau_molecule.log').write_text(
            f'SCF Done: E(PM6) = {obj.energy}.0 A.U.\n'
            'Stationary point found\nNormal termination of Gaussian\n')

    try:
        with patch('shutil.which', return_value='/fake/g16'), patch.object(g16Object, 'run_g16', run):
            return pm6_relax(mol, folder, charge, multiplicity, cores, memory, settings)
    finally:
        assert Path.cwd() == before_cwd
        assert dict(os.environ) == before_env


@pytest.mark.parametrize('nproc', [1, 2])
def test_pm6_parallel_isolation_order_and_failure(tmp_path, monkeypatch, nproc):
    monkeypatch.chdir(tmp_path)
    ref = Chem.AddHs(Chem.MolFromSmiles('CCCC'))
    assert AllChem.EmbedMolecule(ref, randomSeed=42) == 0
    cfg = SearchConfig.resolve('low', {
        'budget': {'formula': 'fixed', 'fixed': 10}, 'mm_method': 'none', 'threads': 1,
        'validation': {'duplicate_rmsd_angstrom': 0}})
    prepared = prepare_candidates(ref, tmp_path/'run', cfg, allocated_cores=nproc,
                                  generators={'etkdgv3': _copies})
    audit = relax_candidates(prepared, ref, cfg, 0, 1, 'pm6', nproc, '1GB', adapter=_pm6_adapter)
    workers = min(nproc, 10)
    assert audit['parallel_workers'] == workers
    assert audit['cores_per_calculation'] == 1
    assert audit['failed_candidates'] == 1 and audit['converged_candidates'] == 9
    assert audit['selected_index'] == 9
    assert [r['prepared_index'] for r in audit['candidate_runs']] == list(range(10))
    progress = json.loads((tmp_path/'run/electronic/progress.json').read_text())
    assert progress == audit['candidate_runs']
    intervals = [json.loads(p.read_text()) for p in (tmp_path/'run/electronic').glob('*/interval.json')]
    events = sorted([(r['start'], 1) for r in intervals] + [(r['end'], -1) for r in intervals])
    active = maximum = 0
    for _, delta in events:
        active += delta
        maximum = max(maximum, active)
    assert (1 if workers == 1 else 2) <= maximum <= workers
    assert len({r['pid'] for r in intervals}) <= workers
