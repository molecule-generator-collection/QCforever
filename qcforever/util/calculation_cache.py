"""Conservative, stage-level reuse shared by Gaussian, GAMESS and ORCA.

The input (including electronic state, geometry and numerical settings) and
external restart dependencies identify a calculation. Logs and required output
files are checked independently; mere file existence never establishes success.
"""

from contextvars import ContextVar
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
        return not re.search(r'Error termination|Convergence failure', log, re.I)
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
    def __init__(self, backend, enabled=True, variant=None):
        if type(enabled) is not bool:
            raise ValueError('reuse must be True or False')
        self.backend, self.enabled = backend, enabled
        self.variant = variant
        self.directory = None
        self.steps = []
        self.protected = set()
        self.lock = None
        self.index = {'schema': 1, 'entries': {}}

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
        # Never modify a directory which may still have a running calculation.
        try:
            descriptor = os.open(lock, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
        except FileExistsError as exc:
            raise RuntimeError(f'Calculation lock exists: {lock}; check the recorded process before removing a stale lock') from exc
        with os.fdopen(descriptor, 'w') as stream:
            stream.write(str(os.getpid()))
        self.lock = lock
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
        previous = [p for p in directory.iterdir() if p.name != DIRECTORY and p.name not in self.protected]
        self.legacy = None
        if previous:
            history = state / 'history' / (time.strftime('%Y%m%dT%H%M%S') + '-' + uuid.uuid4().hex[:8])
            history.mkdir(parents=True)
            for path in previous:
                shutil.move(str(path), history / path.name)
            self.legacy = history
        return directory

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

    def run(self, name, text, execute, executable=None):
        if Path(name).name != name:
            raise ValueError('Calculation job names must not contain directory components')
        deps = dependencies(self.backend, text)
        payload = {'backend': self.backend, 'input': canonical_input(self.backend, text),
                   'dependencies': deps, 'variant': self.variant}
        key = hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()
        program = _program(self.backend, executable)
        required = _required(self.backend, name, text)
        entry = self.index['entries'].get(key)
        reason = 'no verified calculation with matching conditions'
        if self.enabled and entry:
            ok, reason = self._entry_ok(entry, program, required)
            if ok:
                for filename in entry['files']:
                    shutil.copyfile(self.state / entry['snapshot'] / filename, filename)
                self.steps.append(dict(job=name, action='reused', reason=reason, key=key))
                return 'normal' if self.backend != 'gamess' else ''
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
        snapshot = Path('cache') / (key + '-' + uuid.uuid4().hex[:8])
        root = self.state / snapshot
        root.mkdir(parents=True)
        files = {}
        for path in self._outputs(Path.cwd(), name, required):
            shutil.copyfile(path, root / path.name)
            files[path.name] = _hash(root / path.name)
        self.index['entries'][key] = dict(input=text, program=program, variant=self.variant,
                                          files=files, snapshot=str(snapshot), created_ns=time.time_ns())
        _json(self.state / 'index.json', self.index)

    def finish(self):
        report = {'policy': 'reuse' if self.enabled else 'fresh', 'backend': self.backend,
                  'reused': sum(s['action'] == 'reused' for s in self.steps),
                  'executed': sum(s['action'] == 'executed' for s in self.steps),
                  'failed': sum(s['action'] == 'failed' for s in self.steps), 'steps': self.steps}
        try:
            if self.directory is not None:
                _json(self.state / 'reuse-report.json', report)
        finally:
            if self.lock is not None:
                self.lock.unlink(missing_ok=True)
        return report


def prepare_directory(directory, source=None):
    session = _active.get()
    if session is None:
        Path(directory).mkdir(exist_ok=False)
    else:
        session.prepare(directory, source)


def workflow(backend):
    def decorate(method):
        @wraps(method)
        def run(owner, *args, **kwargs):
            session = CalculationCache(backend, getattr(owner, 'reuse', True),
                                       getattr(owner, 'gamessversion', None) if backend == 'gamess' else None)
            token = _active.set(session)
            result = None
            try:
                result = method(owner, *args, **kwargs)
                return result
            finally:
                try:
                    report = session.finish()
                    if isinstance(result, dict):
                        result['reuse_report'] = report
                        if getattr(owner, 'pklsave', False) and session.directory:
                            name = Path(owner.in_file).stem + '.pkl'
                            with (session.directory / name).open('wb') as stream:
                                pickle.dump(result, stream)
                finally:
                    _active.reset(token)
        return run
    return decorate


def execute(backend, name, text, callback, executable=None):
    session = _active.get()
    if session is None or session.backend != backend or session.directory != Path.cwd():
        return callback()
    return session.run(name, text, callback, executable)


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
    raise ValueError(f'{backend} ReadFC: no verified Hessian for {name} with matching atom order, '
                     'method, basis, charge, multiplicity, solvent and excited state. '
                     'Run this stage with FC=CalcFC or CalcAll first; keep its .qcforever directory.')


def preserve_path(directory, path):
    session = _active.get()
    return path.name == DIRECTORY or bool(session and session.directory == Path(directory).absolute()
                                          and path.name in session.protected)


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
