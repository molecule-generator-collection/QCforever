"""Run GFN2-xTB continuously or retain its optimizer in a stopped POSIX process.

No checkpoint restart is used for graybox execution. Each candidate starts
once, owns its process group and log files, and must be closed by its caller.
"""
from contextlib import contextmanager
from dataclasses import asdict
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import select
import signal
import threading
import time
import psutil
from rdkit import Chem
from qcforever.util import job_timeout
from qcforever.util.check_resource import native_thread_environment
from .calculation_logs import write_json, parse_finite_number, read_xtb_trace
from .check_structures import copy_with_coordinates
from .graybox_relaxation import BlockResult


NUMBER = r'[-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?'
XTB_FORCE_KIND = 'anc_norm_per_sqrt_atom'


class XTBCandidate:
    """Advance a single native process at observable cycle endpoints.

The first cycle may include two evaluations. Charge actual work, including
boundary overshoot, rather than treating cycles as energy evaluations.
Construction does not launch a process; untouched candidates cost nothing.
"""

    def __init__(self, mol, folder, charge, multiplicity, memory, settings):
        if os.name != 'posix':
            raise NotImplementedError('xTB graybox requires POSIX signals; use laqa=off on Windows')
        self.executable = shutil.which(settings.get('xtb_executable', 'xtb'))
        if self.executable is None:
            raise FileNotFoundError('xTB executable not found')
        self.mol = Chem.Mol(mol)
        self.folder = Path(folder).resolve()
        self.charge, self.multiplicity = charge, multiplicity
        self.settings = settings
        self.process = self.output = self.trajectory = None
        self.pending = b''
        self.eof = self.stopped = self.finished = self.native_converged = False
        self.started_cycle = self.completed_cycle = self.evaluations = self.scf_cycles = 0
        self.energy = self.gradient = self.final_energy = None
        self.final_molecule = None
        self.calls = 0
        self.active_wall = 0.
        self.scalar_rows = []
        self.cpu_id = None
        self.cancel_requested = lambda: False

    def advance(self):
        """Run one slice and return observed cost, status and scalar endpoint."""
        if self.finished:
            raise RuntimeError('Cannot advance a terminal xTB candidate')
        started = time.monotonic()
        previous_evaluations, previous_scf = self.evaluations, self.scf_cycles
        increment = self.settings['first_interval'] if self.calls == 0 else self.settings['subsequent_interval']
        target = self.completed_cycle + increment
        block = BlockResult(force_kind=XTB_FORCE_KIND, status='failed')
        interrupted = None
        try:
            if self.process is None:
                self._start()
            elif self.stopped:
                self._bind_retained_process()
                os.killpg(self.process.pid, signal.SIGCONT)
                self.stopped = False
            while not self.eof:
                if self.cancel_requested():
                    raise InterruptedError('xTB batch cancelled')
                remaining = job_timeout.remaining_time()
                if remaining is not None and remaining <= 0:
                    raise job_timeout.QCforeverTimeoutError('xTB graybox overall deadline exceeded')
                if self.active_wall + time.monotonic() - started > self.settings.get('xtb_active_timeout_seconds', 1800):
                    raise TimeoutError('xTB candidate active-time limit exceeded')
                readable, _, _ = select.select([self.process.stdout], [], [], .1)
                if readable:
                    self._drain()
                if self.completed_cycle >= target and not self.eof:
                    self._stop()
                    break
            if self.eof:
                self.process.wait(timeout=5)
            self._read_scalars()
            # Only already-written complete XYZ frames may affect selection.
            # Stdout is retained explicitly if that frame is not yet available.
            block.energy = self.energy
            energy_source = 'stdout'
            if self.completed_cycle and len(self.scalar_rows) >= self.completed_cycle:
                block.energy = self.scalar_rows[self.completed_cycle-1]['energy_hartree']
                energy_source = 'trajectory_header'
            block.force = None if self.gradient is None else self.gradient / math.sqrt(self.mol.GetNumAtoms())
            if self.eof:
                if self.process.returncode == 0 and self.native_converged and self.final_energy is not None:
                    final = Chem.MolFromXYZFile(str(self.folder/'xtbopt.xyz'))
                    if final is None or [a.GetAtomicNum() for a in final.GetAtoms()] != [a.GetAtomicNum() for a in self.mol.GetAtoms()]:
                        raise ValueError('xTB final atom order/composition changed')
                    block.energy, energy_source = self.final_energy, 'final_total_energy'
                    self.final_molecule = copy_with_coordinates(self.mol, final)
                    block.status = 'converged'
                elif self.completed_cycle >= self.settings.get('maximum_cycles', 1000):
                    block.status, block.error = 'limit', 'native_iteration_limit'
                else:
                    block.error = 'native_failure_or_missing_convergence'
            else:
                block.status = 'paused'
            if block.status in ('paused', 'converged') and (
                    block.energy is None or not math.isfinite(block.energy) or
                    block.force is None or not math.isfinite(block.force)):
                raise ValueError('Missing or nonfinite xTB endpoint')
        except BaseException as exc:
            block.status, block.error = 'failed', f'{type(exc).__name__}: {exc}'
            energy_source = 'incomplete'
            if isinstance(exc, (job_timeout.QCforeverTimeoutError, KeyboardInterrupt, SystemExit, InterruptedError)):
                interrupted = exc
        finally:
            self.finished = block.status != 'paused'
            if self.finished:
                self.close()
            block.evaluations = self.evaluations - previous_evaluations
            block.scf_cycles = self.scf_cycles - previous_scf
            block.wall_seconds = time.monotonic() - started
            self.active_wall += block.wall_seconds
            self.calls += 1
            if self.output is not None and not self.output.closed:
                self.output.flush()
            record = dict(asdict(block), completed_cycle=self.completed_cycle,
                target_cycle=target, cycle_overshoot=max(0, self.completed_cycle-target),
                energy_source=energy_source, native_process_starts=int(self.process is not None))
            with (self.folder/'blocks.jsonl').open('a') as stream:
                stream.write(json.dumps(record, allow_nan=False)+'\n')
        if interrupted is not None:
            raise interrupted
        return block

    def _start(self):
        # Match the existing continuous backend's XYZ serialization exactly;
        # changing coordinate precision must not confound scheduler comparisons.
        Chem.MolToXYZFile(self.mol, str(self.folder/'input.xyz'))
        argv = [self.executable, 'input.xyz', '--gfn', '2', '--chrg', str(self.charge),
            '--uhf', str(self.multiplicity-1), '--parallel', '1', '--opt', self.settings.get('xtb_opt_level', 'normal'),
            '--cycles', str(self.settings.get('maximum_cycles', 1000)), '--acc', str(self.settings.get('xtb_accuracy', 1.0))]
        if 'xtb_scf_iterations' in self.settings:
            argv += ['--iterations', str(self.settings['xtb_scf_iterations'])]
        # Pin only the child, never change the caller's CPU affinity.
        if hasattr(os, 'sched_getaffinity') and shutil.which('taskset'):
            cpu_id = self.cpu_id if self.cpu_id is not None else min(os.sched_getaffinity(0))
            if cpu_id not in os.sched_getaffinity(0):
                raise ValueError('xTB CPU slot is outside the scheduler allocation')
            argv = ['taskset', '-c', str(cpu_id), *argv]
        env = dict(os.environ, **native_thread_environment(1), OMP_STACKSIZE='256M', GFORTRAN_UNBUFFERED_ALL='y')
        write_json(self.folder/'command.json', dict(argv=argv, force_kind=XTB_FORCE_KIND,
            execution='retained_process', thread_environment=native_thread_environment(1)))
        self.output = (self.folder/'xtb.out').open('xb')
        self.process = job_timeout.popen(argv, cwd=self.folder, env=env, stdout=subprocess.PIPE,
                                        stderr=subprocess.STDOUT, start_new_session=True)
        os.set_blocking(self.process.stdout.fileno(), False)

    def _bind_retained_process(self):
        """Move a stopped candidate to its current batch slot, preserving state."""
        if self.cpu_id is None or not hasattr(os, 'sched_setaffinity'):
            return
        if self.cpu_id not in os.sched_getaffinity(0):
            raise ValueError('xTB CPU slot is outside the scheduler allocation')
        # Existing native threads keep their old affinity unless changed too.
        for thread in psutil.Process(self.process.pid).threads():
            try:
                os.sched_setaffinity(thread.id, {self.cpu_id})
            except ProcessLookupError:
                pass

    def _drain(self):
        while True:
            try:
                data = os.read(self.process.stdout.fileno(), 65536)
            except BlockingIOError:
                return
            if not data:
                self.eof = True
                if self.pending:
                    self._parse_line(self.pending.decode(errors='replace'))
                    self.pending = b''
                return
            self.output.write(data)
            lines = (self.pending + data).split(b'\n')
            self.pending = lines.pop()
            for line in lines:
                self._parse_line(line.decode(errors='replace'))

    def _parse_line(self, line):
        match = re.search(r'\bCYCLE\s+(\d+)', line)
        if match:
            self.started_cycle = int(match[1])
        if re.match(r'^\s*(?:molecular )?gradient\s+\.\.\.', line):
            self.evaluations += 1
        if re.match(r'\s*\d+\s+-?\d+\.\d+\s+[-+]?\d*\.?\d+[EeDd]', line):
            self.scf_cycles += 1
        match = re.search(r'\* total energy\s*:\s*(' + NUMBER + ')', line)
        if match:
            self.energy = parse_finite_number(match[1])
        match = re.match(r'\s*gradient norm\s*:\s*(' + NUMBER + ')', line)
        if match:
            self.gradient = parse_finite_number(match[1])
            self.completed_cycle = self.started_cycle
            with (self.folder/'trajectory.jsonl').open('a') as stream:
                stream.write(json.dumps(dict(cycle=self.completed_cycle, energy_hartree=self.energy,
                    anc_gradient_norm=self.gradient, evaluations=self.evaluations, scf_cycles=self.scf_cycles))+'\n')
        if 'GEOMETRY OPTIMIZATION CONVERGED AFTER' in line:
            self.native_converged = True
        match = re.search(r'TOTAL ENERGY\s+(' + NUMBER + ')', line)
        if match:
            self.final_energy = parse_finite_number(match[1])

    def _stop(self):
        try:
            os.killpg(self.process.pid, signal.SIGSTOP)
        except ProcessLookupError:
            self._drain()
            return
        deadline = time.monotonic() + 2
        native = psutil.Process(self.process.pid)
        while self.process.poll() is None:
            if native.status() == psutil.STATUS_STOPPED:
                self.stopped = True
                break
            if time.monotonic() > deadline:
                raise RuntimeError('xTB SIGSTOP confirmation timed out')
            time.sleep(.0005)
        self._drain()

    def _read_scalars(self):
        path = self.folder/'xtbopt.log'
        if self.trajectory is None and path.exists():
            self.trajectory = path.open()
        if self.trajectory is None:
            return
        while True:
            position = self.trajectory.tell()
            line = self.trajectory.readline()
            if not line:
                return
            if not line.endswith('\n'):
                self.trajectory.seek(position)
                return
            count = int(line.strip())
            header = self.trajectory.readline()
            coordinates = [self.trajectory.readline() for _ in range(count)]
            if not header.endswith('\n') or any(not row.endswith('\n') for row in coordinates):
                self.trajectory.seek(position)
                return
            if count != self.mol.GetNumAtoms():
                raise ValueError('xTB trajectory atom count changed')
            match = re.search(r'energy:\s*('+NUMBER+r')\s+gnorm:\s*('+NUMBER+')', header)
            if match is None:
                raise ValueError('Unrecognized xTB trajectory header')
            row = dict(cycle=len(self.scalar_rows)+1, energy_hartree=parse_finite_number(match[1]),
                       cartesian_gradient_norm=parse_finite_number(match[2]))
            self.scalar_rows.append(row)
            with (self.folder/'trajectory_scalars.jsonl').open('a') as stream:
                stream.write(json.dumps(row)+'\n')

    def close(self):
        """Terminate and reap only our process group, including stopped children."""
        if self.process is not None:
            if self.process.poll() is None:
                try:
                    os.killpg(self.process.pid, signal.SIGTERM)
                    os.killpg(self.process.pid, signal.SIGCONT)
                except ProcessLookupError:
                    pass
                except PermissionError as error:
                    # On macOS the group can disappear between TERM and CONT
                    # with EPERM rather than ESRCH. Ignore it only after reaping
                    # the child; a permission failure on a live child is fatal.
                    try:
                        self.process.wait(timeout=.1)
                    except subprocess.TimeoutExpired:
                        raise error
                try:
                    self.process.wait(timeout=2)
                except subprocess.TimeoutExpired:
                    os.killpg(self.process.pid, signal.SIGKILL)
                    self.process.wait(timeout=5)
            job_timeout.unregister_process(self.process)
            if not self.process.stdout.closed:
                try:
                    self._drain()
                except Exception:
                    pass  # Keep the original failure; raw output is retained.
                self.process.stdout.close()
        for stream in (self.output, self.trajectory):
            if stream is not None and not stream.closed:
                stream.close()
        self.stopped = False


def pool_memory_mb(candidates):
    """Conservative summed RSS (shared pages may be counted more than once)."""
    total = psutil.Process().memory_info().rss
    for candidate in candidates:
        if candidate.process is not None and candidate.process.poll() is None:
            try:
                total += psutil.Process(candidate.process.pid).memory_info().rss
            except psutil.NoSuchProcess:
                pass
    return total / 1024**2


def _interrupt_search(signum, frame):
    raise InterruptedError(f'xTB search interrupted by signal {signum}')


@contextmanager
def cleanup_on_signals():
    """Let the caller's finally block clean the pool after TERM/HUP, not SIGKILL."""
    previous = {}
    if threading.current_thread() is threading.main_thread():
        for signum in (signal.SIGTERM, signal.SIGHUP):
            previous[signum] = signal.signal(signum, _interrupt_search)
    try:
        yield
    finally:
        for signum, handler in previous.items():
            signal.signal(signum, handler)


def optimize_with_xtb(mol, folder, charge, multiplicity, cores, memory, settings):
    executable = shutil.which(settings.get('xtb_executable', 'xtb'))
    if executable is None:
        raise FileNotFoundError('xTB executable not found')
    Chem.MolToXYZFile(mol, str(folder/'input.xyz'))
    argv = [executable, 'input.xyz', '--gfn', '2', '--chrg', str(charge),
            '--uhf', str(multiplicity-1), '--opt', settings.get('xtb_opt_level', 'normal'),
            '--cycles', str(settings.get('maximum_cycles', 1000)), '--parallel', str(cores)]
    if 'xtb_accuracy' in settings:
        argv += ['--acc', str(settings['xtb_accuracy'])]
    if 'xtb_scf_iterations' in settings:
        argv += ['--iterations', str(settings['xtb_scf_iterations'])]
    env = os.environ.copy()
    thread_env = native_thread_environment(cores)
    env.update(thread_env, OMP_STACKSIZE='256M')
    write_json(folder/'command.json', {'argv': argv, 'charge': charge, 'multiplicity': multiplicity,
                                'thread_environment': thread_env})
    try:
        with (folder/'xtb.out').open('w') as out:
            job_timeout.run(argv, cwd=folder, env=env, stdout=out, stderr=subprocess.STDOUT, check=True)
    finally:
        text, trace = read_xtb_trace(folder)
    if not trace['native_converged']:
        raise RuntimeError('xTB native optimization did not converge')
    energies = re.findall(r'TOTAL ENERGY\s+([-+0-9.EeDd]+)', text)
    if not energies:
        raise ValueError('xTB final energy missing')
    final = Chem.MolFromXYZFile(str(folder/'xtbopt.xyz'))
    if final is None or [a.GetAtomicNum() for a in final.GetAtoms()] != [a.GetAtomicNum() for a in mol.GetAtoms()]:
        raise ValueError('xTB final atom order/composition changed')
    return copy_with_coordinates(mol, final), parse_finite_number(energies[-1])
