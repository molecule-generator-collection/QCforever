"""Run PM6 with Gaussian, keeping native execution separate from candidate selection."""
from contextlib import contextmanager
from dataclasses import asdict
import math
import os
from pathlib import Path
import shutil
import subprocess
import time
from rdkit import Chem
from qcforever.util import job_timeout
from qcforever.util.check_resource import native_thread_environment
from .calculation_logs import (write_json, read_pm6_trace, read_pm6_block,
                               parse_finite_number, diagnose_pm6_failure, PM6CalculationError)
from .check_structures import copy_with_coordinates
from .graybox_relaxation import BlockResult


class PM6Candidate:
    """One PM6-RFO optimizer, continued from Gaussian's binary checkpoint.

    MaxCycles is cumulative, not the number of new evaluations. Native calls
    use one core. Paths and environment belong to the subprocess, not the caller.
    Existing call directories are never overwritten or resumed implicitly.
    """

    def __init__(self, mol, folder, charge, multiplicity, memory, settings):
        self.mol = Chem.Mol(mol)
        self.folder = Path(folder).resolve()
        self.charge, self.multiplicity = charge, multiplicity
        self.memory = memory or '1GB'
        self.settings = settings
        self.target = 0
        self.calls = 0
        self.checkpoint = None
        self.final_molecule = None
        self.finished = False
        self.cpu_id = None
        self.cancel_requested = lambda: False
        self.g16 = shutil.which('g16')
        self.formchk = shutil.which('formchk')
        if not self.g16 or not self.formchk:
            raise FileNotFoundError('PM6-RFO requires Gaussian16 g16 and formchk on PATH')

    def advance(self):
        """Advance by one configured interval, retaining failed-call costs."""
        if self.finished:
            raise RuntimeError('Cannot advance a terminal PM6 candidate')
        first = self.checkpoint is None
        interval = self.settings['first_interval'] if first else self.settings['subsequent_interval']
        limit = self.settings.get('maximum_cycles', 1000)
        extra = min(interval, limit - self.target)
        if extra <= 0:
            self.finished = True
            return BlockResult(status='limit', error='native_iteration_limit')
        self.target += extra
        folder = self.folder / f'call_{self.calls:04d}'
        folder.mkdir(exist_ok=False)
        self.calls += 1
        started = time.monotonic()
        result = BlockResult(status='failed')
        parsed = dict(records=[])
        try:
            self._run_call(folder)
        except Exception as exc:
            result.error = f'{type(exc).__name__}: {exc}'
        # Parse even after timeout/nonzero return code, so spent work is retained.
        try:
            log = (folder / 'input.log').read_text(errors='replace') if (folder / 'input.log').exists() else ''
            parsed = read_pm6_block(log, self.mol.GetNumAtoms())
            records = parsed['records']
            result.evaluations = len(records)
            result.scf_cycles = sum(row['scf_cycles'] for row in records)
            if records:
                result.energy = records[-1]['energy_hartree']
                result.force = records[-1]['force_mean_hartree_per_bohr']
            if result.error is None:
                native_short_stop = self._check_call(parsed, extra)
                positions = self._read_final_coordinates(folder)
                if parsed['native_converged']:
                    self.final_molecule = Chem.Mol(self.mol)
                    for i, position in enumerate(positions):
                        self.final_molecule.GetConformer().SetAtomPosition(i, position)
                    result.status = 'converged'
                elif self.target >= limit or native_short_stop:
                    result.status, result.error = 'limit', 'native_iteration_limit'
                else:
                    result.status = 'paused'
                self.checkpoint = folder / 'state.chk'
        except Exception as exc:
            result.status, result.error = 'failed', f'{type(exc).__name__}: {exc}'
        result.wall_seconds = time.monotonic() - started
        self.finished = result.status != 'paused'
        write_json(folder / 'trace.json', parsed)
        write_json(folder / 'result.json', dict(asdict(result),
                   requested_additional_evaluations=extra, requested_cumulative_maxcycles=self.target))
        return result

    def _run_call(self, folder):
        if self.cancel_requested():
            raise InterruptedError('PM6 batch cancelled')
        scratch = folder / 'scratch'
        scratch.mkdir()
        if self.checkpoint is not None:
            shutil.copy2(self.checkpoint, folder / 'state.chk')
        restart = 'Restart,' if self.checkpoint is not None else ''
        route = f'PM6 NoSymm SCF=(Tight,MaxCycle=512) Opt=({restart}RFO,Tight,MaxCycles={self.target})'
        text = f'%chk=state.chk\n%mem={self.memory}\n%nprocshared=1\n#p {route}\n\n'
        if self.checkpoint is None:
            text += f'QCforever PM6-RFO\n\n{self.charge} {self.multiplicity}\n'
            xyz = self.mol.GetConformer().GetPositions()
            for atom, (x, y, z) in zip(self.mol.GetAtoms(), xyz):
                text += f'{atom.GetSymbol()} {x:.14f} {y:.14f} {z:.14f}\n'
            text += '\n'
        (folder / 'input.com').write_text(text)
        env = dict(os.environ, **native_thread_environment(1), GAUSS_SCRDIR=str(scratch))
        command = [self.g16]
        # GENKAI's Gaussian launcher takes its CPU binding from -c, separately
        # from %nprocshared. Never expand the scheduler's available CPU set.
        if hasattr(os, 'sched_getaffinity'):
            cpu_id = self.cpu_id if self.cpu_id is not None else min(os.sched_getaffinity(0))
            if cpu_id not in os.sched_getaffinity(0):
                raise ValueError('PM6 CPU slot is outside the scheduler allocation')
            command.append('-c=' + str(cpu_id))
        command.append('input.com')
        write_json(folder / 'command.json', dict(argv=command, route=route,
                   charge=self.charge, multiplicity=self.multiplicity,
                   thread_environment=native_thread_environment(1)))
        with (folder / 'launcher.out').open('w') as stream:
            process = job_timeout.run(command, cwd=folder, env=env, stdout=stream,
                stderr=subprocess.STDOUT, timeout=self.settings['timeout_seconds_per_call'])
        write_json(folder / 'termination.json', dict(returncode=process.returncode))
        # Gaussian step-limit exits may be nonzero. _check_call distinguishes
        # the expected optimization limit from SCF or other native failures.

    def _check_call(self, parsed, extra):
        records = parsed['records']
        if parsed['scf_failed'] or not (parsed['native_converged'] or parsed['cycle_limit']):
            raise PM6CalculationError('Unexpected Gaussian termination; see input.log')
        if not records or any(row['method'] not in ('RPM6', 'UPM6') for row in records):
            raise PM6CalculationError('Missing PM6 SCF records or changed electronic method')
        if not all(row['energy_hartree'] is not None and math.isfinite(row['energy_hartree']) for row in records):
            raise PM6CalculationError('Nonfinite PM6 energy')
        if records[-1]['force_mean_hartree_per_bohr'] is None:
            raise PM6CalculationError('Missing endpoint force vectors')
        # Gaussian can lower the requested MaxCycles internally. Only accept an
        # explicit native step-limit exit with consistent cumulative counters;
        # an arbitrary short block or step-storage error remains a failure.
        if len(records) > extra:
            raise PM6CalculationError('Unexpected evaluation count in restarted PM6 block')
        if not parsed['native_converged'] and len(records) < extra:
            expected_last_step = self.target - extra + len(records)
            if not parsed['cycle_limit'] or not parsed['step_numbers'] or parsed['step_numbers'][-1] != expected_last_step:
                raise PM6CalculationError('Unexpected evaluation count in restarted PM6 block')
            return True
        return False

    def _read_final_coordinates(self, folder):
        if self.cancel_requested():
            raise InterruptedError('PM6 batch cancelled')
        with (folder / 'formchk.out').open('w') as stream:
            job_timeout.run([self.formchk, 'state.chk', 'state.fchk'], cwd=folder,
                env=dict(os.environ, **native_thread_environment(1)), stdout=stream,
                stderr=subprocess.STDOUT, timeout=60, check=True)
        atomic, positions = read_checkpoint_coordinates(folder / 'state.fchk')
        if atomic != [atom.GetAtomicNum() for atom in self.mol.GetAtoms()]:
            raise PM6CalculationError('Checkpoint atom order/composition changed')
        return positions


def read_checkpoint_coordinates(path):
    """Read final coordinates from fchk, avoiding rounded log-table coordinates."""
    lines = Path(path).read_text().splitlines()
    arrays = {}
    for i, line in enumerate(lines):
        for key in ('Atomic numbers', 'Current cartesian coordinates'):
            if line.startswith(key) and 'N=' in line:
                count = int(line.split('N=')[1])
                values = []
                for row in lines[i + 1:]:
                    values.extend(row.split())
                    if len(values) >= count:
                        break
                arrays[key] = [parse_finite_number(value) for value in values[:count]]
    atomic = [int(value) for value in arrays['Atomic numbers']]
    coordinates = arrays['Current cartesian coordinates']
    if len(coordinates) != 3 * len(atomic):
        raise PM6CalculationError('Incomplete checkpoint coordinates')
    bohr_to_angstrom = 0.529177210903
    positions = [[v * bohr_to_angstrom for v in coordinates[i:i + 3]]
                 for i in range(0, len(coordinates), 3)]
    return atomic, positions


def optimize_with_pm6(mol, folder, charge, multiplicity, cores, memory, settings):
    """Use the existing Gaussian adapter, retaining its files in this folder."""
    if settings.get('pm6_optimizer') == 'rfo':
        # Matched continuous-RFO reference; not a graybox scheduler.
        native_settings = dict(settings, first_interval=settings.get('maximum_cycles', 1000),
                               subsequent_interval=10,
                               timeout_seconds_per_call=settings.get('timeout_seconds_per_call', 1800))
        candidate = PM6Candidate(mol, folder, charge, multiplicity, memory, native_settings)
        block = candidate.advance()
        if block.status != 'converged':
            raise PM6CalculationError(block.error or 'Continuous RFO did not converge')
        return candidate.final_molecule, block.energy
    from qcforever.laqa_fafoom.pyg16 import g16Object
    binary = shutil.which('g16')
    if binary is None:
        raise FileNotFoundError('Gaussian16 executable not found for PM6')
    thread_env = native_thread_environment(cores)
    previous = {key: os.environ.get(key) for key in (*thread_env, 'GAUSS_EXEDIR', 'GAUSS_SCRDIR')}
    try:
        os.environ.update(thread_env)
        write_json(folder/'environment.json', thread_env)
        with working_directory(folder):
            gaussian_job = g16Object(Chem.MolToMolBlock(mol), str(Path(binary).parent), str(folder),
                            cores, memory or '1GB', 'opt', charge, multiplicity, 'pm6',
                            settings.get('maximum_cycles', 1000))
            gaussian_job.generate_input()
            try:
                gaussian_job.run_g16()
            except job_timeout.QCforeverTimeoutError:
                raise
            except Exception as exc:
                path = Path('Gau_molecule.log')
                diagnostic = diagnose_pm6_failure(path.read_text(errors='replace') if path.is_file() else '', exc)
                write_json(Path('failure.json'), diagnostic)
                raise PM6CalculationError(
                    f"PM6: {diagnostic['reason']}; see {folder/'failure.json'} and Gau_molecule.log") from exc
            finally:
                path = Path('Gau_molecule.log')
                log = path.read_text(errors='replace') if path.is_file() else ''
                trace = read_pm6_trace(log)
                converged = trace['native_converged']
                write_json(Path('trace.json'), trace)
            if not converged:
                diagnostic = diagnose_pm6_failure(log)
                write_json(Path('failure.json'), diagnostic)
                raise PM6CalculationError(f"PM6: {diagnostic['reason']}; see {folder/'failure.json'}")
            # The Gaussian adapter edits coordinates in an SDF template. Read
            # positions without chemical sanitization, then retain the input
            # graph/electronic annotations exactly, as for the xTB backend.
            out = Chem.MolFromMolBlock(gaussian_job.get_sdf_string_opt(), removeHs=False, sanitize=False)
            if out is None:
                raise ValueError('Unreadable PM6 final structure')
            return copy_with_coordinates(mol, out), gaussian_job.get_energy('hartree')
    finally:
        for key, value in previous.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value


@contextmanager
def working_directory(path):
    previous = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)
