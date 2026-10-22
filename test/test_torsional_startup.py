"""Check shared-reader/exclusive-writer behavior without ML dependencies."""
import os
import subprocess
import sys

import pytest

from qcforever.conformer_search.model_workers import torsional_diffusion as td

pytestmark = pytest.mark.skipif(os.name != 'posix', reason='TD uses POSIX file locks')


@pytest.mark.parametrize('warm', [False, True])
def test_concurrent_initialization_and_unlocked_model_load(tmp_path, warm):
    source=tmp_path/'source'
    source.mkdir()
    (source/'diffusion').mkdir()
    (source/'diffusion/__init__.py').write_text('')
    # Fake the upstream writer, while using real process-level file locks.
    (source/'diffusion/torus.py').write_text('''
from pathlib import Path
import time
if not Path('.p.npy').exists():
    with Path('writer_count').open('a') as stream:
        stream.write('writer\\n')
    Path('.p.npy').write_text('p')
    time.sleep(.1)
    Path('.score.npy').write_text('score')
assert Path('.score.npy').read_text() == 'score'
''')
    cache=tmp_path/'cache'
    cache.mkdir()
    if warm:
        (cache/'.p.npy').write_text('p')
        (cache/'.score.npy').write_text('score')
    script='''
import fcntl, importlib.util, sys
from pathlib import Path
from types import SimpleNamespace
spec=importlib.util.spec_from_file_location('td',sys.argv[1])
td=importlib.util.module_from_spec(spec)
spec.loader.exec_module(td)
args=SimpleNamespace(source=Path(sys.argv[2]),cache=Path(sys.argv[3]))
def load(args,request):
    # No table lock may remain held during model loading.
    if sys.argv[4] == 'check_unlocked':
        with (args.cache/'initialize.lock').open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    return 'loaded'
td._load_model=load
assert td.initialize_model(args,{'threads':1}) == 'loaded'
assert args.initialization_timings['cache_lock_wait_seconds'] >= 0
'''
    # Test unlocked model load separately: another cold-start worker may still
    # briefly own an exclusive lock before downgrading to a shared reader.
    command = [sys.executable, '-c', script, td.__file__, str(source), str(cache)]
    workers = [subprocess.Popen(command + ['parallel']) for _ in range(4)]
    statuses = [process.wait(timeout=10) for process in workers]
    assert statuses == [0] * 4
    count=cache/'writer_count'
    if warm:
        assert not count.exists()
    else:
        assert count.read_text() == 'writer\n'
    subprocess.run(command + ['check_unlocked'], check=True, timeout=10)


def test_partial_cache_is_not_silently_overwritten(tmp_path):
    (tmp_path/'.p.npy').write_bytes(b'partial')
    with pytest.raises(RuntimeError,match='Incomplete TD lookup tables'):
        td._initialize_torsion_tables(tmp_path)
    assert (tmp_path/'.p.npy').read_bytes() == b'partial'


def test_warm_cache_uses_shared_lock(tmp_path, monkeypatch):
    import fcntl
    for name in ['.p.npy','.score.npy']:
        (tmp_path/name).touch()
    operations=[]
    real_flock=fcntl.flock
    def record(fd, operation):
        operations.append(operation)
        return real_flock(fd,operation)
    monkeypatch.setattr(fcntl,'flock',record)
    monkeypatch.setattr(td.importlib,'import_module',lambda name: None)
    td._initialize_torsion_tables(tmp_path)
    assert operations and all(op==fcntl.LOCK_SH for op in operations)
