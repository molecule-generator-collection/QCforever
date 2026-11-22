"""Conservative, stage-level reuse shared by Gaussian, GAMESS and ORCA.

The input (including electronic state, geometry and numerical settings) and
external restart dependencies identify a calculation. Logs and required output
files are checked independently; mere file existence never establishes success.
"""

from contextvars import ContextVar
import fcntl
from functools import wraps
import hashlib
import json
import os
from pathlib import Path
import pickle
import re
import shutil
import time
import uuid


_active = ContextVar('qcforever_calculation_cache', default=None)
DIRECTORY = '.qcforever'
_KTLC_RECORD = 'ktlc_reuse'


def _hash(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def _json(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + '.' + uuid.uuid4().hex + '.tmp')
    try:
        temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def _pickle(path, value):
    temporary = path.with_name(path.name + '.' + uuid.uuid4().hex + '.tmp')
    try:
        with temporary.open('wb') as stream:
            pickle.dump(value, stream, protocol=pickle.HIGHEST_PROTOCOL)
            stream.flush()
            os.fsync(stream.fileno())
        temporary.replace(path)
        descriptor = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(descriptor)
        finally:
            os.close(descriptor)
    finally:
        temporary.unlink(missing_ok=True)


class _ResultUnpickler(pickle.Unpickler):
    """Read QCforever result data without instantiating arbitrary classes."""

    def find_class(self, module, name):
        original = module.replace('numpy._core.', 'numpy.core.', 1)
        if (original, name) in {
            ('numpy', 'dtype'), ('numpy', 'ndarray'),
            ('numpy.core.multiarray', '_reconstruct'), ('numpy.core.multiarray', 'scalar'),
            ('numpy.core.numeric', '_frombuffer'),
        }:
            return super().find_class(module, name)
        raise pickle.UnpicklingError(f'Unsupported result type {module}.{name}')


def _read_result(path):
    if not path.is_file() or path.is_symlink():
        return {}
    try:
        with path.open('rb') as stream:
            result = _ResultUnpickler(stream).load()
    except (OSError, ValueError, ImportError, AttributeError, EOFError, pickle.UnpicklingError):
        return {}
    return result if isinstance(result, dict) else {}


def canonical_input(backend, text):
    """Ignore allocation limits, retaining all physical/numerical controls."""
    if backend == 'gaussian':
        text = re.sub(r'^\s*%(?:mem|nproc|nprocshared|cpu)\s*=.*$', '', text, flags=re.I | re.M)
    elif backend == 'orca':
        text = re.sub(r'%pal\b.*?\bend\b', '', text, flags=re.I | re.S)
        text = re.sub(r'^\s*%maxcore\s+[^\n]+', '', text, flags=re.I | re.M)
    elif backend == 'gamess':
        text = re.sub(r'\b(?:MWORDS|MEMDDI|TIMLIM)\s*=\s*[^\s$]+', '', text, flags=re.I)
    return ' '.join(text.lower().split())


def dependencies(backend, text):
    """Files whose content is read by this input, not the requested outputs."""
    names = []
    if backend == 'gaussian' and re.search(r'\b(?:guess\s*=\s*(?:\([^)]*read|read)|geom\s*=\s*(?:\([^)]*check|\w*check))', text, re.I):
        old = re.findall(r'^\s*%oldchk\s*=\s*(\S+)', text, re.I | re.M)
        current = re.findall(r'^\s*%chk\s*=\s*(\S+)', text, re.I | re.M)
        names = old or current
        names = [name if Path(name).suffix else name + '.chk' for name in names]
    elif backend == 'orca':
        names += re.findall(r'%moinp\s+"([^"]+)"', text, re.I)
        names += re.findall(r'\bInHessName\s+"([^"]+)"', text, re.I)
        names += re.findall(r'^\s*\*\s*xyzfile\s+\S+\s+\S+\s+(\S+)', text, re.I | re.M)
    return {name: _hash(name) if Path(name).is_file() else None for name in names}


def successful(backend, text, log):
    """Check completion/convergence rather than trusting launcher exit status."""
    if backend == 'orca':
        if 'ORCA TERMINATED NORMALLY' not in log[-12000:]:
            return False
        if re.search(r'SCF NOT CONVERGED|SCF did not converge', log, re.I):
            return False
        if re.search(r'\b(?:Opt|LooseOpt|TightOpt|VeryTightOpt)\b', text, re.I):
            if 'THE OPTIMIZATION HAS CONVERGED' not in log:
                return False
        if not re.search(r'FINAL SINGLE POINT ENERGY\s+[-+\d.]', log):
            return False
        from qcforever.orca_run.parse_log import parse_log
        parsed = parse_log.__new__(parse_log)
        parsed.filename, parsed.text = '<reuse candidate>', log
        try:
            parsed.energy()
            if re.search(r'\bNumFreq\b|\bFreq\b', text, re.I):
                parsed.frequencies()
            if re.search(r'\bPolar\s+(?:true|1)\b', text, re.I):
                parsed.polar()
            if re.search(r'\bNMR\b', text, re.I):
                parsed.shielding()
            if re.search(r'%tddft\b', text, re.I):
                parsed.transitions()
        except (ValueError, IndexError):
            return False
        return True
    if backend == 'gaussian':
        if 'Normal termination of Gaussian' not in log[-4000:]:
            return False
        if re.search(r'Error termination|Convergence failure', log, re.I):
            return False
        density = re.search(r'\bDensity\s*=\s*(Current|SCF)\b', text, re.I)
        if density:
            root = re.search(r'\bTD\s*=?\s*\([^)]*\broot\s*=\s*(\d+)', text, re.I)
            from qcforever.gaussian_run.parse_log import parse_log
            try:
                if density[1].lower() == 'scf':
                    parse_log.Extract_PopulationDipole(log.splitlines(), 'SCF')
                elif root:
                    parse_log.Extract_PopulationDipole(log.splitlines(), 'CI', state=int(root[1]))
            except ValueError:
                return False
        return True
    if 'EXECUTION OF GAMESS TERMINATED NORMALLY' not in log:
        return False
    if re.search(r'SCF IS UNCONVERGED|VIBRATIONAL ANALYSIS IS NOT VALID|ABNORMALLY', log, re.I):
        return False
    if re.search(r'RUNTYP\s*=\s*OPTIMIZE', text, re.I):
        return 'EQUILIBRIUM GEOMETRY LOCATED' in log
    return True


def _program(backend, executable=None):
    path = shutil.which(str(executable or {'gaussian': 'g16', 'gamess': 'rungms', 'orca': 'orca'}[backend]))
    if not path:
        return None
    path = Path(path).resolve()
    try:
        stat = path.stat()
    except OSError:
        return None
    return {'path': str(path), 'size': stat.st_size, 'mtime_ns': stat.st_mtime_ns}


def _required(backend, name, text):
    result = [name + '.log']
    if backend == 'gamess':
        result.append(name + '.dat')
    elif backend == 'gaussian':
        for value in re.findall(r'^\s*%chk\s*=\s*(\S+)', text, re.I | re.M):
            result.append(value if Path(value).suffix else value + '.chk')
    return result


class CalculationCache:
    def __init__(self, backend, enabled=True, variant=None, *, ktlc_job=None):
        if type(enabled) is not bool:
            raise ValueError('reuse must be True or False')
        self.backend, self.enabled = backend, enabled
        self.variant = variant
        self.directory = None
        self.steps = []
        self.protected = set()
        self.lock = None
        self.lock_descriptor = None
        self.index = {'schema': 1, 'entries': {}}
        self.ktlc_job = ktlc_job
        self.ktlc = {}
        self.previous_ktlc = {}
        self.previous_result = {}
        self.completed = {}
        self.request_key = None

    def prepare(self, directory, source=None):
        directory = Path(directory).absolute()
        if directory.is_symlink():
            raise ValueError('A calculation directory must not be a symbolic link')
        directory = directory.resolve()
        directory.mkdir(exist_ok=True)
        state = directory / DIRECTORY
        if state.is_symlink():
            raise ValueError('The reuse state must not be a symbolic link')
        marker = state / 'owner.json'
        if state.exists() and (not marker.is_file() or json.loads(marker.read_text()) != {'application': 'QCforever', 'schema': 1}):
            raise ValueError(f'Unrecognized reuse state: {state}')
        state.mkdir(exist_ok=True)
        _json(marker, {'application': 'QCforever', 'schema': 1})
        lock = state / 'running.lock'
        self._lock_directory(lock)
        self.directory = directory
        self.state = state
        index = state / 'index.json'
        if index.exists():
            self.index = json.loads(index.read_text())
            if self.index.get('schema') != 1 or not isinstance(self.index.get('entries'), dict):
                raise ValueError(f'Unsupported or damaged reuse index: {index}')
        # Preserve an input located inside the calculation directory itself.
        if source:
            try:
                self.protected.add(Path(source).resolve().relative_to(directory).parts[0])
            except (ValueError, IndexError):
                pass
        if self.ktlc_job:
            self.protected.add(self.ktlc_job + '.pkl')
            self.protected.add('qcforever-history')
            if self.enabled:
                self.previous_result = _read_result(directory / (self.ktlc_job + '.pkl'))
                previous = self.previous_result.get(_KTLC_RECORD)
                if (self.previous_result.get('log') == 'normal' and isinstance(previous, dict)
                        and previous.get('schema') == 1 and previous.get('backend') == self.backend):
                    self.previous_ktlc = previous
                record = state / 'ktlc.json'
                if record.is_file():
                    saved = json.loads(record.read_text())
                    if saved.get('backend') == self.backend:
                        self.ktlc = saved
        previous = [p for p in directory.iterdir() if p.name != DIRECTORY and p.name not in self.protected]
        self.legacy = None
        if previous:
            history = state / 'history' / (time.strftime('%Y%m%dT%H%M%S') + '-' + uuid.uuid4().hex[:8])
            history.mkdir(parents=True)
            for path in previous:
                shutil.move(str(path), history / path.name)
            self.legacy = history
        return directory

    def _lock_directory(self, path):
        # An advisory lock is released by the OS even on scheduler cancellation.
        # Check for a legacy PID-only lock before taking over an older directory.
        descriptor = os.open(path, os.O_RDWR | os.O_CREAT | os.O_NOFOLLOW, 0o600)
        try:
            fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
            if not os.path.samestat(os.fstat(descriptor), path.stat()):
                raise BlockingIOError('Calculation directory changed while taking the lock')
            record = os.read(descriptor, 128)
            if record and not record.startswith(b'flock '):
                raise BlockingIOError('Legacy PID-only lock')
            os.lseek(descriptor, 0, os.SEEK_SET)
            os.write(descriptor, f'flock {os.getpid()}\n'.encode())
            os.ftruncate(descriptor, os.lseek(descriptor, 0, os.SEEK_CUR))
        except (BlockingIOError, FileNotFoundError) as exc:
            os.close(descriptor)
            raise RuntimeError(f'Calculation lock exists: {path}; check the recorded process for a legacy lock') from exc
        except BaseException:
            os.close(descriptor)
            raise
        self.lock_descriptor = descriptor
        self.lock = path

    def save_ktlc(self, section, data):
        self.ktlc.update(schema=1, backend=self.backend)
        self.ktlc[section] = data
        _json(self.state / 'ktlc.json', self.ktlc)

    def selected_conformer(self, source, charge, multiplicity, identity=None):
        if not self.ktlc_job or not self.enabled:
            return None
        source_key = [_hash(source), int(charge), int(multiplicity)]
        for results in (self.ktlc, self.previous_ktlc):
            saved = results.get('conformer', {})
            if (saved.get('source') != source_key
                    or (identity is not None and saved['identity'] != identity)):
                continue
            Path('optimized_structures.sdf').write_text(saved['structure'])
            for root in self._history_roots():
                selected = root / 'optimized_structures.sdf'
                details = root / 'conformer_search'
                if (selected.is_file() and details.is_dir() and not details.is_symlink()
                        and selected.read_text() == saved['structure']):
                    shutil.move(details, Path.cwd() / details.name)
                    break
            self.save_ktlc('conformer', saved)
            self.steps.append(dict(job='optconf', action='reused',
                                   reason='verified KTLC conformer selection', key=saved['identity']))
            return saved['summary'].copy()

    def _saved_stage(self, key, program, required):
        entry = self.previous_ktlc.get('stages', {}).get(key)
        if not entry or entry['program'] != program:
            return None
        logfile = required[0]
        if entry['logfile'] != logfile:
            return None
        for root in self._history_roots():
            log = self._saved_log(entry, root)
            if log is None:
                continue
            files = {}
            for name, digest in entry['files'].items():
                artifact = root / name
                if Path(name).name == name and artifact.is_file() and _hash(artifact) == digest:
                    files[name] = artifact
            if all(name in files for name in required[1:]):
                return log, files

    def _history_roots(self):
        if self.legacy is not None:
            yield self.legacy
        history = self.state / 'history'
        if history.is_dir():
            yield from (root for root in sorted(history.iterdir(), reverse=True)
                        if root != self.legacy and root.is_dir() and not root.is_symlink())

    def _saved_log(self, entry, root):
        name = entry['logfile']
        path = root / name
        if Path(name).name != name or not path.is_file():
            return None
        with path.open('rb') as stream:
            start, size = entry['log_range']
            stream.seek(start)
            log = stream.read(size)
        if hashlib.sha256(log).hexdigest() != entry['log_hash']:
            return None
        if successful(self.backend, entry['input'], log.decode(errors='replace')):
            return log

    def _gaussian_reseed(self, key, program, logfile, recover):
        entry = self.previous_ktlc.get('stages', {}).get(key)
        geometry = self.previous_ktlc.get('gaussian_geometries', {}).get(key)
        if geometry is None or not entry or entry['program'] != program or entry['logfile'] != logfile:
            return None
        for root in self._history_roots():
            log = self._saved_log(entry, root)
            if log is not None and b'Stationary point found' in log:
                if recover(geometry):
                    self.save_ktlc('gaussian_geometries', dict(
                        self.ktlc.get('gaussian_geometries', {}), **{key: geometry}))
                    return log
                return None

    def _finished_stages(self):
        records = {}
        last_logfile, combined_log = None, b''
        artifact_hashes = {}
        for key, (name, entry) in sorted(self.completed.items(), key=lambda stage: stage[1][0]):
            root = self.state / entry['snapshot']
            logfile = name + '.log'
            saved = root / logfile
            current = self.directory / logfile
            if not saved.is_file() or not current.is_file():
                continue
            log = saved.read_bytes()
            if logfile != last_logfile:
                last_logfile, combined_log = logfile, current.read_bytes()
            start = combined_log.find(log)
            if start < 0:
                continue
            files = {}
            for filename, digest in entry['files'].items():
                artifact = self.directory / filename
                if filename == logfile or not artifact.is_file():
                    continue
                if filename not in artifact_hashes:
                    artifact_hashes[filename] = _hash(artifact)
                if artifact_hashes[filename] == digest:
                    files[filename] = digest
            records[key] = dict(input=entry['input'], program=entry['program'],
                                logfile=logfile, log_range=(start, len(log)),
                                log_hash=entry['files'][logfile], files=files)
        return records

    def _entry_ok(self, entry, program, required):
        if program is not None and entry.get('program') != program:
            return False, 'executable changed or its provenance is unknown'
        files = entry.get('files', {})
        if not all(name in files for name in required):
            return False, 'required restart output is missing'
        root = self.state / entry['snapshot']
        if not root.resolve().is_relative_to(self.state.resolve()):
            return False, 'invalid snapshot path'
        for name, digest in files.items():
            path = root / name
            if Path(name).name != name or not path.is_file() or _hash(path) != digest:
                return False, 'cached output is missing or was modified'
        log = (root / required[0]).read_text(errors='replace')
        if not successful(self.backend, entry['input'], log):
            return False, 'cached log is incomplete or unconverged'
        return True, 'matching input, dependencies, executable and verified outputs'

    def _legacy_candidate(self, name, text, dependencies_now, required):
        # Old input files alone are not evidence of what actually ran. Recover
        # complete executed inputs from ORCA/GAMESS echoes. Gaussian logs do not
        # echo checkpoint dependencies, so those legacy jobs are recomputed.
        if self.legacy is None or dependencies_now:
            return None
        logfile = self.legacy / (name + '.log')
        if not logfile.is_file() or not all((self.legacy / path).is_file() for path in required):
            return None
        log = logfile.read_text(errors='replace')
        if self.backend == 'orca':
            echoed = '\n'.join(re.findall(r'^\|\s*\d+> ?(.*)$', log, re.M))
            echoed = echoed.split('****END OF INPUT****')[0]
        elif self.backend == 'gamess':
            echoed = '\n'.join(re.findall(r'^.*?INPUT CARD> ?(.*)$', log, re.M))
        else:
            return None
        if canonical_input(self.backend, echoed) != canonical_input(self.backend, text) or not successful(self.backend, text, log):
            return None
        return self.legacy

    def run(self, name, text, execute, executable=None, *, required_outputs=(), recover=None):
        if Path(name).name != name:
            raise ValueError('Calculation job names must not contain directory components')
        deps = dependencies(self.backend, text)
        payload = {'backend': self.backend, 'input': canonical_input(self.backend, text),
                   'dependencies': deps, 'variant': self.variant}
        key = hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()
        program = _program(self.backend, executable)
        required = [*_required(self.backend, name, text), *required_outputs]
        entry = self.index['entries'].get(key)
        reason = 'no verified calculation with matching conditions'
        if self.enabled and entry:
            ok, reason = self._entry_ok(entry, program, required)
            if ok:
                for filename in entry['files']:
                    shutil.copyfile(self.state / entry['snapshot'] / filename, filename)
                self.steps.append(dict(job=name, action='reused', reason=reason, key=key))
                self.completed[key] = (name, entry)
                return 'normal' if self.backend != 'gamess' else ''
        if self.enabled and self.previous_ktlc:
            saved = self._saved_stage(key, program, required)
            if saved is not None:
                log, files = saved
                Path(required[0]).write_bytes(log)
                for filename, source in files.items():
                    shutil.copyfile(source, filename)
                self.steps.append(dict(job=name, action='reused', reason='verified stage from the completed KTLC result', key=key))
                self._store(key, text, program, name, required)
                return 'normal' if self.backend != 'gamess' else ''
            if recover is not None:
                log = self._gaussian_reseed(key, program, required[0], recover)
                if log is not None and all(Path(path).is_file() for path in required[1:]):
                    Path(required[0]).write_bytes(log)
                    self.steps.append(dict(job=name, action='reused', key=key,
                                           reason='restored optimized geometry and rebuilt Gaussian SCF'))
                    self._store(key, text, program, name, required)
                    return 'normal'
        if self.enabled and entry is None:
            legacy = self._legacy_candidate(name, text, deps, required)
            if legacy is not None:
                for path in self._outputs(legacy, name, required):
                    shutil.copyfile(path, path.name)
                self.steps.append(dict(job=name, action='reused', reason='executed input recovered from the complete legacy log', key=key))
                self._store(key, text, program, name, required)
                return 'normal' if self.backend != 'gamess' else ''
        if not self.enabled:
            reason = 'reuse=False requested a fresh calculation'
        step = dict(job=name, action='executed', reason=reason, key=key)
        self.steps.append(step)
        Path(name + '.log').unlink(missing_ok=True)
        try:
            result = execute()
        except BaseException:
            step['action'] = 'failed'
            raise
        logfile = Path(name + '.log')
        if (result in ('', 'normal') and logfile.is_file()
                and successful(self.backend, text, logfile.read_text(errors='replace'))
                and all(Path(path).is_file() for path in required)
                and all(value is not None for value in deps.values())):
            self._store(key, text, program, name, required)
        else:
            step['action'] = 'failed'
            step['reason'] = 'completion, convergence or required outputs could not be verified'
        return result

    @staticmethod
    def _outputs(root, name, required):
        allowed = {'.log', '.chk', '.fchk', '.dat', '.gbw', '.xyz', '.hess', '.txt'}
        return [path for path in root.iterdir() if path.is_file() and not path.is_symlink()
                and (path.name in required or (path.name.startswith(name + '.') and path.suffix in allowed))]

    def _store(self, key, text, program, name, required):
        previous = self.index['entries'].get(key)
        snapshot = Path('cache') / (key + '-' + uuid.uuid4().hex[:8])
        root = self.state / snapshot
        root.mkdir(parents=True)
        files = {}
        for path in self._outputs(Path.cwd(), name, required):
            shutil.copyfile(path, root / path.name)
            files[path.name] = _hash(root / path.name)
        self.index['entries'][key] = dict(input=text, program=program, variant=self.variant,
                                          files=files, snapshot=str(snapshot), created_ns=time.time_ns())
        self.completed[key] = (name, self.index['entries'][key])
        _json(self.state / 'index.json', self.index)
        if self.ktlc_job and previous:
            old = self.state / previous['snapshot']
            if old.resolve().is_relative_to((self.state / 'cache').resolve()) and old.is_dir():
                shutil.rmtree(old)

    def _clear_cache(self):
        root = self.state / 'cache'
        if root.exists() or root.is_symlink():
            shutil.rmtree(root)
        if self.index['entries']:
            self.index['entries'].clear()
            _json(self.state / 'index.json', self.index)

    def _finish_ktlc(self, result, pickle_name):
        path = self.directory / pickle_name if pickle_name else None
        if not isinstance(result, dict):
            return
        last_steps = {step['job']: step['action'] for step in self.steps}
        if (result.get('log') != 'normal' or result.get('optconf') is False
                or 'failed' in last_steps.values()):
            if path is not None and not path.exists():
                _pickle(path, result)
            return
        if path is not None:
            result[_KTLC_RECORD] = dict(self.ktlc, schema=1, backend=self.backend,
                                        stages=self._finished_stages(),
                                        request_key=self.request_key if self.ktlc.get('tuning') else None)
            _pickle(path, result)
        try:
            self._preserve_history_files()
            obsolete = self.state.with_name(DIRECTORY + '-' + uuid.uuid4().hex)
            self.state.rename(obsolete)
            self.lock = None
            shutil.rmtree(obsolete)
        except OSError as exc:
            print(f'Could not remove KTLC recovery data in {self.directory}: {exc}', flush=True)

    def _preserve_history_files(self):
        """Keep other files a user may have placed in an existing job directory."""
        history = self.state / 'history'
        outputs = {'.com', '.gjf', '.inp', '.log', '.out', '.chk', '.fchk',
                   '.xyz', '.gbw', '.hess', '.dat', '.pkl'}
        jobs = {self.ktlc_job}
        for entry in self.index['entries'].values():
            jobs.update(Path(name).stem for name in entry.get('files', {}) if Path(name).suffix.lower() in outputs)
        for entry in self.previous_ktlc.get('stages', {}).values():
            jobs.add(Path(entry['logfile']).stem)
            jobs.update(Path(name).stem for name in entry['files'])
        jobs.update(job + '.readfc' for job in tuple(jobs))
        for folder in history.iterdir() if history.is_dir() else ():
            for path in folder.iterdir():
                if (path.name in {'conformer_search', 'optimized_structures.sdf'}
                        and not path.is_symlink()):
                    continue
                if (path.is_file() and not path.is_symlink() and path.suffix.lower() in outputs
                        and path.stem in jobs):
                    continue
                destination = self.directory / path.name
                if destination.exists() or destination.is_symlink():
                    destination = self.directory / 'qcforever-history' / folder.name / path.name
                    destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.move(path, destination)

    def finish(self, result=None, *, pickle_name=None, cleanup_cache=False):
        report = {'policy': 'reuse' if self.enabled else 'fresh', 'backend': self.backend,
                  'reused': sum(s['action'] == 'reused' for s in self.steps),
                  'executed': sum(s['action'] == 'executed' for s in self.steps),
                  'failed': sum(s['action'] == 'failed' for s in self.steps), 'steps': self.steps}
        try:
            if isinstance(result, dict):
                result['reuse_report'] = report
            if self.directory is not None:
                _json(self.state / 'reuse-report.json', report)
                if self.ktlc_job:
                    self._finish_ktlc(result, pickle_name)
                elif isinstance(result, dict) and pickle_name is not None:
                    _pickle(self.directory / pickle_name, result)
                if (not self.ktlc_job and cleanup_cache and isinstance(result, dict) and result.get('log') == 'normal'
                        and report['failed'] == 0):
                    try:
                        self._clear_cache()
                    except OSError as exc:
                        print(f'Could not remove calculation cache in {self.state}: {exc}', flush=True)
        finally:
            if self.lock is not None:
                self.lock.unlink(missing_ok=True)
            if self.lock_descriptor is not None:
                os.close(self.lock_descriptor)
        return report


def prepare_directory(directory, source=None):
    session = _active.get()
    if session is None:
        Path(directory).mkdir(exist_ok=False)
    else:
        session.prepare(directory, source)


def previous_output(path):
    """Find an output moved into this invocation's history, if it exists."""
    path = Path(path).absolute()
    if path.is_file():
        return path
    session = _active.get()
    if session is not None and session.legacy is not None and path.is_relative_to(session.directory):
        return session.legacy / path.relative_to(session.directory)
    return path


def workflow(backend):
    def decorate(method):
        @wraps(method)
        def run(owner, *args, **kwargs):
            cleanup_cache = getattr(owner, 'cleanup_cache', None)
            if cleanup_cache is None:
                cleanup_cache = os.environ.get('QCFOREVER_CLEANUP_CACHE') == '1'
            if type(cleanup_cache) is not bool:
                raise ValueError('cleanup_cache must be True, False or None')
            original_functional = getattr(owner, 'functional', '')
            original_parameters = getattr(owner, 'para_functional', None)
            functional = str(original_functional).strip().upper()
            ktlc_job = Path(owner.in_file).stem if functional.startswith('KTLC-') else None
            session = CalculationCache(backend, getattr(owner, 'reuse', True),
                                       getattr(owner, 'gamessversion', None) if backend == 'gamess' else None,
                                       ktlc_job=ktlc_job)
            if ktlc_job:
                from .ktlc_reuse import request_key, unchanged_result
                session.request_key = request_key(owner, session)
                saved = unchanged_result(session)
                if saved is not None:
                    return saved
            token = _active.set(session)
            result = None
            try:
                result = method(owner, *args, **kwargs)
                return result
            finally:
                try:
                    pickle_name = None
                    if isinstance(result, dict) and getattr(owner, 'pklsave', False) and session.directory:
                        pickle_name = Path(owner.in_file).stem + '.pkl'
                    session.finish(result, pickle_name=pickle_name, cleanup_cache=cleanup_cache)
                finally:
                    if ktlc_job:
                        owner.functional = original_functional
                        owner.para_functional = original_parameters
                    _active.reset(token)
        return run
    return decorate


def execute(backend, name, text, callback, executable=None, *, required_outputs=(), recover=None):
    session = _active.get()
    if session is None or session.backend != backend or session.directory != Path.cwd():
        return callback()
    return session.run(name, text, callback, executable, required_outputs=required_outputs, recover=recover)


def remember_gaussian_geometry(name, formatted_checkpoint):
    """Save full-precision coordinates independently of restart file retention."""
    session = _active.get()
    if session is None or not session.ktlc_job:
        return
    matches = [key for key, (job, entry) in session.completed.items()
               if job == name and re.search(r'^\s*Opt=\(', entry['input'], re.I | re.M)]
    if len(matches) != 1:
        return
    from .read_mol_file import read_fchk
    try:
        atoms, x, y, z, charge, mult = read_fchk(formatted_checkpoint)
    except (OSError, ValueError):
        return
    geometry = dict(atoms=list(atoms), coordinates=[[float(v) for v in row] for row in zip(x, y, z)],
                    charge=int(charge), multiplicity=int(mult))
    geometries = session.ktlc.get('gaussian_geometries', {}).copy()
    geometries[matches[0]] = geometry
    session.save_ktlc('gaussian_geometries', geometries)


def has_gaussian_geometry():
    session = _active.get()
    return bool(session and session.enabled and session.previous_ktlc.get('gaussian_geometries'))


def previous_hessian(backend, name, text, executable=None):
    """Find a hash-verified native Hessian with matching electronic conditions."""
    from .optimization_hessian import recipe, calculates_hessian, gamess_block
    session = _active.get()
    if session is not None and session.backend == backend and session.directory == Path.cwd():
        expected = recipe(backend, text)
        program = _program(backend, executable)
        for entry in sorted(session.index['entries'].values(),
                            key=lambda item: item.get('created_ns', 0), reverse=True):
            old = entry['input']
            artifact = name + ('.hess' if backend == 'orca' else '.dat')
            if (entry.get('variant') != session.variant
                    or artifact not in entry.get('files', {}) or not calculates_hessian(backend, old)
                    or recipe(backend, old) != expected):
                continue
            required = _required(backend, name, old)
            ok, _ = session._entry_ok(entry, program, required)
            if not ok:
                continue
            path = session.state / entry['snapshot'] / artifact
            if backend == 'gamess':
                try:
                    gamess_block(path.read_text())
                except ValueError:
                    continue
            return path
        for entry in session.previous_ktlc.get('stages', {}).values():
            old = entry['input']
            artifact = name + ('.hess' if backend == 'orca' else '.dat')
            digest = entry['files'].get(artifact)
            if (not digest or entry['program'] != program or not calculates_hessian(backend, old)
                    or recipe(backend, old) != expected):
                continue
            for root in session._history_roots():
                path = root / artifact
                if not path.is_file() or _hash(path) != digest or session._saved_log(entry, root) is None:
                    continue
                if backend == 'gamess':
                    try:
                        gamess_block(path.read_text())
                    except ValueError:
                        continue
                return path
    raise ValueError(f'{backend} ReadFC: no verified Hessian for {name} with matching atom order, '
                     'method, basis, charge, multiplicity, solvent and excited state. '
                     'Run this stage with FC=CalcFC or CalcAll first and retain its restart files.')


def preserve_path(directory, path):
    session = _active.get()
    return path.name == DIRECTORY or bool(session and session.directory == Path(directory).absolute()
                                          and path.name in session.protected)


def selected_ktlc_conformer(source, charge, multiplicity):
    session = _active.get()
    if session is not None:
        return session.selected_conformer(source, charge, multiplicity)


def ktlc_scratch_root():
    session = _active.get()
    return session.state if session and session.ktlc_job else Path.cwd()


def conformer_workflow(method):
    """Reuse a completed conformer search, including its retained histories."""
    @wraps(method)
    def run(infilename, charge, multiplicity, backend, nproc, memory, config):
        session = _active.get()
        if session is None or session.directory != Path.cwd():
            return method(infilename, charge, multiplicity, backend, nproc, memory, config)
        from qcforever.conformer_search.settings import SearchConfig
        if not isinstance(config, SearchConfig):
            config = SearchConfig.load(config) if isinstance(config, (str, Path)) else SearchConfig.from_mapping(config)
        identity = {'source': _hash(infilename), 'charge': int(charge), 'multiplicity': int(multiplicity),
                    'backend': backend, 'config': config.to_mapping(),
                    'program': _program('orca', config.relaxation.get('xtb_executable', 'xtb')) if backend == 'xtb' else _program('gaussian')}
        key = 'conformer-' + hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
        summary = session.selected_conformer(infilename, charge, multiplicity, key)
        if summary is not None:
            return summary
        entry = session.index['entries'].get(key)
        reason = 'no completed search with matching structure, state and settings'
        if session.enabled and entry:
            root = session.state / entry['snapshot']
            files = entry.get('files', {})
            valid = root.resolve().is_relative_to(session.state.resolve()) and all(
                not Path(name).is_absolute() and '..' not in Path(name).parts
                and (root / name).is_file() and _hash(root / name) == digest
                for name, digest in files.items())
            valid = valid and {'optimized_structures.sdf', 'conformer_search/summary.json'} <= files.keys()
            if valid:
                for name in files:
                    destination = Path(name)
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copyfile(root / name, destination)
                result = json.loads(Path('conformer_search/summary.json').read_text())
                if str(result.get('state', '')).startswith('succeeded'):
                    if session.ktlc_job:
                        _store_ktlc_conformer(session, infilename, charge, multiplicity, key, result)
                    session.steps.append(dict(job='optconf', action='reused', reason='verified completed conformer search', key=key))
                    return result
            reason = 'conformer search outputs are incomplete or were modified'
        step = dict(job='optconf', action='executed', reason=reason if session.enabled else 'reuse=False', key=key)
        session.steps.append(step)
        try:
            result = method(infilename, charge, multiplicity, backend, nproc, memory, config)
        except BaseException:
            step['action'] = 'failed'
            raise
        selected = Path('optimized_structures.sdf')
        if str(result.get('state', '')).startswith('succeeded') and selected.is_file():
            if session.ktlc_job:
                _store_ktlc_conformer(session, infilename, charge, multiplicity, key, result)
                return result
            snapshot = Path('cache') / (key + '-' + uuid.uuid4().hex[:8])
            root = session.state / snapshot
            files = {}
            for path in [selected, *Path('conformer_search').rglob('*')]:
                if path.is_file() and not path.is_symlink():
                    destination = root / path
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copyfile(path, destination)
                    files[str(path)] = _hash(destination)
            session.index['entries'][key] = dict(snapshot=str(snapshot), files=files)
            _json(session.state / 'index.json', session.index)
        else:
            step['action'] = 'failed'
        return result
    return run


def _store_ktlc_conformer(session, source, charge, multiplicity, key, summary):
    session.save_ktlc('conformer', dict(
        source=[_hash(source), int(charge), int(multiplicity)], identity=key,
        structure=Path('optimized_structures.sdf').read_text(), summary=summary,
    ))
