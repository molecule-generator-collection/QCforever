"""Check real subprocess lifetime/IPC without requiring learned model packages."""
import json
import sys
import time
import pytest
from rdkit import Chem
from qcforever.conformer_search.model_workers.model_process import ModelSession
from qcforever.conformer_search.generate_conformers import GenerationError


@pytest.fixture
def command(tmp_path):
    worker = tmp_path/'fake_worker.py'
    worker.write_text('''import argparse,json,os,pathlib,time
from rdkit import Chem
p=argparse.ArgumentParser()
p.add_argument('--request'); p.add_argument('--output'); p.add_argument('--session')
a=p.parse_args(); root=pathlib.Path(a.session); number=0
while not (root/'stop').exists():
    path=root/'next.json'
    if not path.exists(): time.sleep(.01); continue
    msg=json.loads(path.read_text()); path.unlink(); number+=1
    r=json.loads(pathlib.Path(msg['request']).read_text())
    if r['seed']==999: raise SystemExit(7)
    with Chem.SDWriter(msg['output']) as w:
        for i in range(r['maximum_raw_candidates']):
            m=Chem.MolFromSmiles(r['smiles']); m.SetIntProp('pid',os.getpid()); m.SetIntProp('number',number); w.write(m)
    dst=pathlib.Path(msg['response']); temp=dst.with_suffix('.tmp')
    temp.write_text(json.dumps({'state':'completed'})); temp.replace(dst)
''')
    return [sys.executable, str(worker), '--request', '{request}', '--output', '{output}']


def test_workers_reused_when_additional_batch_is_smaller(tmp_path, command, monkeypatch):
    monkeypatch.chdir(tmp_path)
    from pathlib import Path
    root = Path('stage'); root.mkdir()
    session = ModelSession(root, {'command': command, 'timeout_seconds': 10}, 1)
    ref = Chem.MolFromSmiles('CCO')
    def batch(n, size):
        folder=root/f'batch_{n}'; folder.mkdir()
        return session.generate(ref, size, 10+n, folder, 2)
    try:
        first = batch(0, 2)
        second = batch(1, 1)
        third = batch(2, 2)
        pids = [m.GetIntProp('pid') for m in first]
        assert len(set(pids)) == 2
        assert second[0].GetIntProp('pid') == pids[0]
        assert [m.GetIntProp('pid') for m in third] == pids
        assert [m.GetIntProp('number') for m in third] == [3, 2]
    finally:
        session.close()
    assert all(p.poll() is not None for p, _, _ in session.workers.values())


def test_worker_crash_not_silently_restarted(tmp_path, command):
    root=tmp_path/'stage'; root.mkdir()
    session=ModelSession(root, {'command': command, 'timeout_seconds': 10}, 1)
    folder=root/'batch'; folder.mkdir()
    try:
        with pytest.raises(GenerationError, match='exited'):
            session.generate(Chem.MolFromSmiles('CCO'), 1, 999, folder, 1)
        process = session.workers[0][0]
        assert process.returncode == 7
    finally:
        session.close()
