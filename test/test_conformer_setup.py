"""Basic installer safety tests; no packages or weights are downloaded."""
import hashlib
import io
import subprocess
import stat
import sys
import zipfile

import pytest

from qcforever.conformer_search.model_workers import install_models as install, installed_models as registry


@pytest.fixture(autouse=True)
def isolated_registry(tmp_path, monkeypatch):
    monkeypatch.setenv('QCFOREVER_MODEL_REGISTRY', str(tmp_path/'config/models.json'))


def test_registration_and_explicit_yaml_precedence():
    from qcforever.conformer_search.settings import SearchConfig
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


def test_dry_run_has_no_writes(tmp_path, monkeypatch):
    monkeypatch.setattr(install, 'prepare_assets', lambda *a: pytest.fail('download attempted'))
    root = tmp_path/'not-created'
    assert install.main(['--dry-run', '--device', 'cpu', '--directory', str(root)]) == 0
    assert not root.exists()
    assert not registry.registry_path().exists()


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


def test_nested_worker_deploys_without_core_dependencies(tmp_path, monkeypatch):
    """Use the real package copy, but do not install/download ML dependencies."""
    site = tmp_path/'model-site'
    def run(argv, log, *, capture=False):
        if argv[1:3] == ['-m', 'venv']:
            binary = argv[3]/'bin/python'
            binary.parent.mkdir()
            binary.touch()
        if capture:
            return str(site) if 'sysconfig' in str(argv) else ''
    monkeypatch.setattr(install, 'run', run)
    executable = install.install_environment('ditmc', 'cpu', sys.executable, tmp_path, None)
    assert executable.is_file()
    assert (site/'qcforever_model_workers/model_sources.json').is_file()
    assert (site/'qcforever_model_workers/requirements/ditmc-cpu.txt').is_file()
    # -S prevents installed core/ML packages from masking an import dependency.
    script = ('import sys; sys.path.insert(0, sys.argv[1]); '
              'from qcforever_model_workers import run_model, check_installation, torsional_diffusion; '
              'assert callable(torsional_diffusion.initialize_model); '
              'assert callable(torsional_diffusion.generate_conformers); '
              'assert not any(m in sys.modules for m in ("qcforever", "numpy", "torch", "jax"))')
    subprocess.run([sys.executable, '-S', '-c', script, str(site)], cwd=tmp_path, check=True)
    options = install.worker_options('ditmc', executable, tmp_path, tmp_path, tmp_path)
    assert options['command'][1:3] == ['-m', 'qcforever_model_workers.run_model']


@pytest.mark.parametrize('name', ['ditmc', 'torsional_diffusion'])
def test_model_dispatch_after_module_rename(name, monkeypatch, tmp_path):
    """Check the real dispatch path without downloading model weights."""
    from types import SimpleNamespace
    from qcforever.conformer_search.model_workers import run_model, ditmc, torsional_diffusion
    model = object()
    monkeypatch.setattr(ditmc, 'DiTMC', lambda *a, **k: model)
    monkeypatch.setattr(torsional_diffusion, 'initialize_model', lambda *a: model)
    args = SimpleNamespace(model=name, cache=tmp_path, source=tmp_path, checkpoint=tmp_path)
    assert run_model.initialize_model(args, {'threads': 1}) is model
