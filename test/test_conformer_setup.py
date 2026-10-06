"""Installer safety/routing tests without downloading ML packages or weights."""
import hashlib
import io
import json
import os
from pathlib import Path
import stat
import sys
from types import SimpleNamespace
import zipfile

import pytest

from qcforever_model_workers import install, registry, smoke


@pytest.fixture(autouse=True)
def isolated_registry(tmp_path, monkeypatch):
    monkeypatch.setenv('QCFOREVER_MODEL_REGISTRY', str(tmp_path/'config/models.json'))


def test_registration_and_explicit_yaml_precedence():
    from qcforever.conformer_search.config import SearchConfig
    registry.register_models({'ditmc': {'persistent': True, 'command': ['registered']}})
    registry.register_models({'torsional_diffusion': {'persistent': True, 'command': ['td']}})
    cfg = SearchConfig.resolve('high')
    assert cfg.generators['ditmc']['command'] == ['registered']
    assert cfg.generators['torsional_diffusion']['persistent']
    cfg = SearchConfig.resolve('high', {'generators': {'ditmc': {'command': ['explicit']}}})
    assert cfg.generators['ditmc']['command'] == ['explicit']
    assert cfg.generators['ditmc']['persistent']
    assert SearchConfig.resolve().profile == 'low'
    assert SearchConfig.resolve(override={'profile': 'high'}).generators['ditmc']['command'] == ['registered']


def test_missing_and_malformed_registry():
    assert registry.read_registry()['generators'] == {}
    registry.write_json(registry.registry_path(), {'schema_version': 9})
    with pytest.raises(ValueError, match='Invalid model registry'):
        registry.read_registry()


def test_download_validated_and_cached(tmp_path, monkeypatch):
    data = b'verified test content'
    item = {'name': 'test.zip', 'url': 'https://example.invalid/a',
            'size': len(data), 'checksum': 'sha256:'+hashlib.sha256(data).hexdigest()}
    calls = []
    def urlopen(*args, **kwargs):
        calls.append(args)
        return io.BytesIO(data)
    monkeypatch.setattr(install.urllib.request, 'urlopen', urlopen)
    target = install.download(item, tmp_path)
    assert target.read_bytes() == data
    assert install.download(item, tmp_path) == target
    assert len(calls) == 1
    target.write_bytes(b'corrupt')
    with pytest.raises(RuntimeError, match='Cached download failed'):
        install.download(item, tmp_path)


def test_bad_download_never_promoted(tmp_path, monkeypatch):
    monkeypatch.setattr(install.urllib.request, 'urlopen', lambda *a, **k: io.BytesIO(b'HTML'))
    spec = {'name': 'model.pt', 'url': 'https://example.invalid', 'size': 4, 'checksum': 'sha256:'+'0'*64}
    with pytest.raises(RuntimeError, match='checksum'):
        install.download(spec, tmp_path)
    assert not (tmp_path/'model.pt').exists()


@pytest.mark.parametrize('name', ['../escape', '/absolute', 'dir/../../escape', 'dir\\escape'])
def test_zip_traversal_refused(tmp_path, name):
    archive = tmp_path/'unsafe.zip'
    with zipfile.ZipFile(archive, 'w') as handle:
        handle.writestr(name, b'bad')
    with pytest.raises(ValueError, match='Unsafe archive'):
        install.extract_zip(archive, tmp_path/'out')


def test_zip_symlink_refused(tmp_path):
    archive = tmp_path/'link.zip'
    item = zipfile.ZipInfo('link')
    item.external_attr = (stat.S_IFLNK | 0o777) << 16
    with zipfile.ZipFile(archive, 'w') as handle:
        handle.writestr(item, '../../outside')
    with pytest.raises(ValueError, match='symlink'):
        install.extract_zip(archive, tmp_path/'out')


def test_zip_selects_only_model_files(tmp_path):
    archive = tmp_path/'source.zip'
    with zipfile.ZipFile(archive, 'w') as handle:
        handle.writestr('source/model/a.py', b'code')
        handle.writestr('source/data/large.pkl', b'unneeded')
    install.extract_zip(archive, tmp_path/'out', ['source/model/'])
    assert (tmp_path/'out/source/model/a.py').read_bytes() == b'code'
    assert not (tmp_path/'out/source/data').exists()


def test_dry_run_has_no_writes(tmp_path, monkeypatch):
    monkeypatch.setattr(install, 'prepare_assets', lambda *a: pytest.fail('download attempted'))
    root = tmp_path/'not-created'
    assert install.main(['--dry-run', '--device', 'cpu', '--directory', str(root)]) == 0
    assert not root.exists()
    assert not registry.registry_path().exists()


def test_unsupported_platform_fails_before_writes(tmp_path, monkeypatch):
    monkeypatch.setattr(install.sys, 'platform', 'darwin')
    root = tmp_path/'not-created'
    assert install.main(['--device', 'cpu', '--directory', str(root), '--python', sys.executable]) == 1
    assert not root.exists()


def test_refuses_unmanaged_directory(tmp_path):
    (tmp_path/'existing.txt').write_text('preserve me')
    with pytest.raises(RuntimeError, match='unmanaged'):
        with install.setup_lock(tmp_path):
            pass
    assert (tmp_path/'existing.txt').read_text() == 'preserve me'


def test_partial_failure_preserves_old_registration(tmp_path, monkeypatch):
    registry.register_models({'ditmc': {'command': ['old-working']}})
    monkeypatch.setattr(install, 'check_platform', lambda *a: None)
    def setup(model, args, device, python, output):
        (output/model).mkdir()
        if model == 'ditmc':
            raise RuntimeError('test model failure')
        registry.register_models({model: {'command': ['new-working']}})
    monkeypatch.setattr(install, 'setup_model', setup)
    assert install.main(['--directory', str(tmp_path/'managed'), '--device', 'cpu', '--python', sys.executable]) == 1
    models = registry.read_registry()['generators']
    assert models['ditmc']['command'] == ['old-working']
    assert models['torsional_diffusion']['command'] == ['new-working']


def test_masked_gpu_selects_cpu(monkeypatch):
    monkeypatch.setenv('CUDA_VISIBLE_DEVICES', '')
    monkeypatch.setattr(install.subprocess, 'run', lambda *a, **k: pytest.fail('must honor mask'))
    assert install.choose_device('auto') == 'cpu'
    assert install.choose_device('gpu') == 'gpu'


def test_smoke_real_ipc_with_fake_model(tmp_path):
    worker = tmp_path/'fake.py'
    worker.write_text('''import argparse,json,os,pathlib,time
from rdkit import Chem
from rdkit.Chem import AllChem
p=argparse.ArgumentParser();p.add_argument('--request');p.add_argument('--output');p.add_argument('--session')
a=p.parse_args();s=pathlib.Path(a.session);n=0
while not (s/'stop').exists():
 f=s/'next.json'
 if not f.exists():time.sleep(.01);continue
 m=json.loads(f.read_text());f.unlink();n+=1
 mol=Chem.AddHs(Chem.MolFromSmiles('CCO'));AllChem.EmbedMolecule(mol,randomSeed=n)
 with Chem.SDWriter(m['output']) as w:w.write(mol)
 folder=pathlib.Path(m['output']).parent
 (folder/'model_execution.json').write_text(json.dumps({'worker_pid':os.getpid(),'model_reused':n>1,'initialization_seconds':0,'device':{'effective':'cpu'}}))
 r=pathlib.Path(m['response']);t=r.with_suffix('.tmp');t.write_text('{}');t.replace(r)
''')
    options = {'command': [sys.executable, str(worker), '--request', '{request}', '--output', '{output}']}
    summary = smoke.check(options, tmp_path/'smoke', 'cpu', 1, 15)
    assert summary['state'] == 'passed'
    assert summary['rows'][0]['worker_pid'] == summary['rows'][1]['worker_pid']


def test_collapsed_coordinates_do_not_pass(tmp_path):
    from rdkit import Chem
    m = Chem.AddHs(Chem.MolFromSmiles('CCO'))
    m.AddConformer(Chem.Conformer(m.GetNumAtoms()))
    path = tmp_path/'bad.sdf'
    with Chem.SDWriter(str(path)) as writer:
        writer.write(m)
    with pytest.raises(RuntimeError, match='collisions'):
        smoke.verify_records(path)


@pytest.mark.parametrize('passed', [True, False])
def test_setup_registers_only_after_confirmed_smoke(tmp_path, monkeypatch, passed):
    registry.register_models({'ditmc': {'command': ['old']}})
    monkeypatch.setattr(install, 'prepare_assets', lambda *a: (tmp_path/'source', tmp_path/'weights'))
    monkeypatch.setattr(install, 'install_environment', lambda *a: Path(sys.executable))
    def run(argv, log):
        output = Path(argv[argv.index('--output')+1])
        registry.write_json(output/'summary.json', {'state': 'passed' if passed else 'failed'})
    monkeypatch.setattr(install, 'run', run)
    output = tmp_path/'logs'
    output.mkdir()
    args = SimpleNamespace(directory=tmp_path/'managed', threads=4, timeout=15)
    if passed:
        install.setup_model('ditmc', args, 'cpu', sys.executable, output)
        options = registry.read_registry()['generators']['ditmc']
        assert options['persistent']
        assert options['command'][0] == sys.executable
    else:
        with pytest.raises(RuntimeError, match='did not confirm success'):
            install.setup_model('ditmc', args, 'cpu', sys.executable, output)
        assert registry.read_registry()['generators']['ditmc']['command'] == ['old']
