"""One-command, user-local installation of optional conformer generators.

The current Python environment is never modified. Downloads are pinned and
verified before extraction/deserialization. Only real smoke-test successes are
registered; a failed update leaves a previous registration intact.
"""
import argparse
from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import platform
import shutil
import stat
import subprocess
import sys
import time
import urllib.request
import zipfile

from .registry import read_registry, register_models, registry_path, write_json

PACKAGE = Path(__file__).resolve().parent
SOURCES = json.loads((PACKAGE/'model_sources.json').read_text())


def digest(path, algorithm='sha256'):
    value = hashlib.new(algorithm)
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b''):
            value.update(chunk)
    return value.hexdigest()


def download(spec, directory):
    """Cache only verified complete files; interrupted .part files are retried."""
    directory.mkdir(parents=True, exist_ok=True)
    target = directory/spec['name']
    algorithm, expected = spec['checksum'].split(':', 1)
    def valid(path):
        return path.stat().st_size == spec['size'] and digest(path, algorithm) == expected
    if target.exists():
        if not valid(target):
            raise RuntimeError(f'Cached download failed checksum; move it aside and retry: {target}')
        return target
    print(f"Downloading {spec['name']} ({spec['size']/1e6:.1f} MB)", flush=True)
    partial = target.with_suffix(target.suffix+'.part')
    request = urllib.request.Request(spec['url'], headers={'User-Agent': 'QCforever-conformer-setup'})
    with urllib.request.urlopen(request, timeout=120) as response, partial.open('wb') as stream:
        total, reported = 0, 0
        while chunk := response.read(1024*1024):
            total += len(chunk)
            if total > spec['size']:
                raise RuntimeError(f'Download exceeded expected size: {spec["name"]}')
            stream.write(chunk)
            if total-reported >= 100*1024*1024:
                print(f'  {total/1e6:.0f} / {spec["size"]/1e6:.0f} MB', flush=True)
                reported = total
    if not valid(partial):
        raise RuntimeError(f'Download failed size/checksum validation: {partial}')
    partial.replace(target)
    return target


def extract_zip(archive, directory, prefixes=None):
    """Extract regular files only, with traversal/symlink/size safeguards."""
    directory.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive) as handle:
        total = 0
        for item in handle.infolist():
            name = PurePosixPath(item.filename)
            if name.is_absolute() or '..' in name.parts or '\\' in item.filename:
                raise ValueError(f'Unsafe archive member: {item.filename}')
            if stat.S_ISLNK(item.external_attr >> 16):
                raise ValueError(f'Archive symlink refused: {item.filename}')
            if prefixes and not any(item.filename.startswith(p) for p in prefixes):
                continue
            total += item.file_size
            if total > 8*1024**3:
                raise ValueError('Archive exceeds extraction size limit')
            target = directory.joinpath(*name.parts)
            if not target.resolve().is_relative_to(directory.resolve()):
                raise ValueError(f'Archive path escapes destination: {target}')
            if item.is_dir():
                target.mkdir(parents=True, exist_ok=True)
            else:
                target.parent.mkdir(parents=True, exist_ok=True)
                with handle.open(item) as src, target.open('wb') as dst:
                    shutil.copyfileobj(src, dst)


def prepare_assets(model, root):
    spec = SOURCES[model]
    fingerprint = hashlib.sha256(json.dumps(spec, sort_keys=True).encode()).hexdigest()[:12]
    directory = root/'assets'/f'{model}-{fingerprint}'
    marker = directory/'complete.json'
    if not marker.exists():
        for item in spec['archives']:
            archive = download(item, root/'downloads')
            extract_zip(archive, directory, item.get('prefixes'))
        for item in spec['files']:
            src = download(item, root/'downloads')
            dst = directory/item['destination']
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)
        write_json(marker, spec)
    return directory/spec['source_directory'], directory/spec['checkpoint_directory']


def run(argv, log, *, capture=False):
    print('  '+ ' '.join(map(str, argv)), flush=True)
    # Never inherit a caller PYTHONPATH/user-site into the dedicated environments.
    env = {k: v for k, v in os.environ.items() if k not in ('PYTHONPATH', 'PYTHONHOME')}
    env['PYTHONNOUSERSITE'] = '1'
    result = subprocess.run(list(map(str, argv)), env=env, check=True,
                            stdout=subprocess.PIPE if capture else log, stderr=log, text=True)
    return result.stdout if capture else None


def install_environment(model, device, python, root, log):
    recipe = PACKAGE/'requirements'/('ditmc-cpu.txt' if model == 'ditmc' else 'torsional-cpu.txt')
    # A changed adapter gets a new environment. Never overwrite a registered
    # working installation during an update that may subsequently fail.
    content = recipe.read_bytes()+b''.join(p.read_bytes() for p in sorted(PACKAGE.glob('*.py')))
    fingerprint = hashlib.sha256(content).hexdigest()[:12]
    folder = root/'envs'/f'{model}-{device}-{fingerprint}'
    marker = folder/'.qcforever-environment.json'
    executable = folder/'bin/python'
    if folder.exists() and not marker.exists():
        raise RuntimeError(f'Refusing to change an unmarked environment: {folder}')
    if (folder/'installed.json').exists() and executable.exists():
        print(f'  Reusing completed environment: {folder}', flush=True)
        return executable
    if not folder.exists():
        folder.mkdir(parents=True)
        write_json(marker, {'model': model, 'device': device, 'recipe_sha256': digest(recipe)})
    if not (folder/'bin/python').exists():
        run([python, '-m', 'venv', folder], log)
    pip = [executable, '-m', 'pip', '--isolated', 'install']
    if model == 'torsional_diffusion':
        backend = 'cu124' if device == 'gpu' else 'cpu'
        run(pip+['torch==2.6.0', '--index-url', f'https://download.pytorch.org/whl/{backend}'], log)
        run(pip+['torch-scatter==2.1.2', 'torch-cluster==1.6.3', '--only-binary=:all:',
                 '-f', f'https://data.pyg.org/whl/torch-2.6.0+{backend}.html'], log)
    extra = (['jax-cuda12-plugin[with_cuda]==0.5.1', 'jax-cuda12-pjrt==0.5.1']
             if model == 'ditmc' and device == 'gpu' else [])
    run(pip+['-r', recipe]+extra, log)
    # Copy our small, pure-Python worker package from this QCforever installation.
    # No second checkout, editable install, or Gaussian dependencies are needed.
    site = Path(run([executable, '-c', 'import sysconfig; print(sysconfig.get_path("purelib"))'], log, capture=True).strip())
    shutil.copytree(PACKAGE, site/'qcforever_model_workers', dirs_exist_ok=True,
                    ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
    run([executable, '-m', 'pip', 'check'], log)
    versions = run([executable, '-m', 'pip', 'freeze'], log, capture=True)
    (folder/'installed-packages.txt').write_text(versions)
    write_json(folder/'installed.json', {'fingerprint': fingerprint, 'device': device})
    return executable


def worker_options(model, python, source, checkpoint, cache):
    return {'persistent': True, 'timeout_seconds': 1800, 'command': [
        str(python), '-m', 'qcforever_model_workers.worker', '--model', model,
        '--source', str(source), '--checkpoint', str(checkpoint), '--cache', str(cache),
        '--request', '{request}', '--output', '{output}']}


def choose_device(requested):
    if requested != 'auto':
        return requested
    if os.environ.get('CUDA_VISIBLE_DEVICES') in ('', '-1'):
        return 'cpu'
    try:
        result = subprocess.run(['nvidia-smi', '-L'], capture_output=True, text=True, timeout=10)
        if result.returncode == 0 and 'GPU ' in result.stdout:
            return 'gpu'
    except (OSError, subprocess.TimeoutExpired):
        pass
    return 'cpu'


def check_platform(python, models):
    if sys.platform != 'linux' or platform.machine() != 'x86_64':
        raise RuntimeError('Automatic model setup currently supports Linux x86_64 only')
    version = subprocess.check_output([python, '-c', 'import sys; print("%d.%d" % sys.version_info[:2])'], text=True).strip()
    if version != '3.11':
        raise RuntimeError('Model environments require Python 3.11; supply --python /path/to/python3.11')
    subprocess.run([python, '-c', 'import venv, ensurepip'], check=True)
    if 'ditmc' in models and not (shutil.which('c++') and shutil.which('cc')):
        raise RuntimeError('DiTMC requires a C/C++ compiler (cc and c++); no system packages are installed by this command')
    if 'ditmc' in models:
        header = subprocess.check_output([python, '-c',
            'import pathlib,sysconfig; print(pathlib.Path(sysconfig.get_path("include"))/"Python.h")'], text=True).strip()
        if not Path(header).is_file():
            raise RuntimeError(f'DiTMC requires Python development headers: missing {header}')


@contextmanager
def setup_lock(root):
    import fcntl
    if root.exists() and not (root/'.qcforever-setup.json').exists() and any(root.iterdir()):
        raise RuntimeError(f'Refusing to use a nonempty, unmanaged setup directory: {root}')
    root.mkdir(parents=True, exist_ok=True)
    with (root/'setup.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(f'Another setup is using {root}') from exc
        write_json(root/'.qcforever-setup.json', {'schema_version': 1})
        yield


def setup_model(model, args, device, python, run_root):
    folder = run_root/model
    folder.mkdir()
    print(f'[{model}] Setting up {device}; log: {folder}/setup.log', flush=True)
    with (folder/'setup.log').open('w') as log:
        source, checkpoint = prepare_assets(model, args.directory)
        executable = install_environment(model, device, python, args.directory, log)
        cache = args.directory/'cache'/model/device
        cache.mkdir(parents=True, exist_ok=True)
        options = worker_options(model, executable, source, checkpoint, cache)
        write_json(folder/'options.json', options)
        print(f'[{model}] Testing actual ethanol generation and model reuse (1 + 1 candidates)', flush=True)
        run([executable, '-m', 'qcforever_model_workers.smoke', '--options', folder/'options.json',
             '--output', folder/'smoke', '--device', device, '--threads', str(args.threads),
             '--timeout', str(args.timeout)], log)
        summary = json.loads((folder/'smoke/summary.json').read_text())
        if summary.get('state') != 'passed':
            raise RuntimeError('Smoke test did not confirm success')
        register_models({model: options})
        write_json(folder/'result.json', {'state': 'passed', 'source': SOURCES[model],
                                        'options': options, 'smoke': summary})
    print(f'[{model}] PASSED and registered', flush=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    base = Path(os.environ.get('XDG_DATA_HOME', Path.home()/'.local/share'))
    parser.add_argument('--directory', type=Path, default=base/'qcforever/conformers')
    parser.add_argument('--models', nargs='+', choices=tuple(SOURCES), default=list(SOURCES))
    parser.add_argument('--device', choices=('auto', 'cpu', 'gpu'), default='auto',
                        help='auto selects a visible NVIDIA GPU, otherwise CPU; use gpu explicitly on a GPU node')
    parser.add_argument('--python', default=sys.executable if sys.version_info[:2] == (3, 11) else 'python3.11')
    parser.add_argument('--threads', type=int, default=4, help='CPU cores for the sequential smoke tests')
    parser.add_argument('--timeout', type=int, default=1800, help='seconds per smoke-test batch')
    parser.add_argument('--dry-run', action='store_true', help='show the plan without downloading or writing files')
    args = parser.parse_args(argv)
    args.directory = args.directory.expanduser().resolve()
    args.models = list(dict.fromkeys(args.models))
    if args.threads < 1 or args.timeout < 1:
        parser.error('threads and timeout must be positive')
    device = choose_device(args.device)
    print(f'Install: {args.directory}\nRegister: {registry_path()}\nModels: {", ".join(args.models)}\nDevice: {device}', flush=True)
    print('Downloads include official third-party code and weights under their own licenses. '
          'Existing Python environments are not changed. Use --help for setup options.', flush=True)
    if args.dry_run:
        print(json.dumps({name: SOURCES[name] for name in args.models}, indent=2))
        return 0
    try:
        python = shutil.which(args.python)
        if not python:
            raise RuntimeError(f'Python executable not found: {args.python}')
        check_platform(python, args.models)
        read_registry()  # Report malformed registration before expensive installation.
        with setup_lock(args.directory):
            run_root = args.directory/'logs'/f'{time.strftime("%Y%m%d-%H%M%S")}-{os.getpid()}'
            run_root.mkdir(parents=True)
            failures = []
            for model in args.models:
                try:
                    setup_model(model, args, device, python, run_root)
                except Exception as exc:
                    failures.append(model)
                    write_json(run_root/model/'result.json', {'state': 'failed', 'error': str(exc)})
                    print(f'[{model}] FAILED: {exc}\nSee {run_root/model}. Previous registration was not replaced.', file=sys.stderr, flush=True)
            if failures:
                print('Setup incomplete. Fix the reported issue and rerun the same command; verified downloads are reused.', file=sys.stderr)
                return 1
        print('Setup complete. Use optconf_medium / optconf_high; no per-job YAML is required.')
        return 0
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as exc:
        print(f'Setup failed: {exc}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
