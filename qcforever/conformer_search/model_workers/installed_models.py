"""User-local model registration shared by setup and conformer search.

Only explicit setup writes this file. Normal calculations never install or
download anything. A per-job YAML takes precedence over registered models.
"""
import json
import os
from pathlib import Path
import tempfile


def registry_path():
    override = os.environ.get('QCFOREVER_MODEL_REGISTRY')
    base = Path(os.environ.get('XDG_CONFIG_HOME', Path.home()/'.config'))
    return Path(override).expanduser().resolve() if override else base/'qcforever/conformer_models.json'


def read_registry():
    path = registry_path()
    if not path.exists():
        return {'schema_version': 1, 'generators': {}}
    value = json.loads(path.read_text())
    if (not isinstance(value, dict) or value.get('schema_version') != 1
            or not isinstance(value.get('generators'), dict)):
        raise ValueError(f'Invalid model registry: {path}')
    return value


def write_json(path, value):
    """Atomic replacement; never expose a partially written registration."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode='w', dir=path.parent, delete=False) as stream:
        temporary = Path(stream.name)
        json.dump(value, stream, indent=2)
        stream.write('\n')
    try:
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def register_models(generators):
    # Different setup roots may share this registry. Lock the read/merge/write.
    import fcntl
    path = registry_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.with_suffix('.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        value = read_registry()
        value['generators'].update(generators)
        write_json(path, value)
