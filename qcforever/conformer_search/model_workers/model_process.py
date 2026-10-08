"""Persistent model subprocesses, scoped to one generator stage.

File-backed atomic IPC keeps model stdout in logs. No silent restart or extra
sampling occurs after errors. Loaded models live until stage completion.
"""
from concurrent.futures import ThreadPoolExecutor
from contextvars import copy_context
import json
import os
from pathlib import Path
import signal
import subprocess
import time
from rdkit import Chem
from qcforever.util import job_timeout


class ModelSession:
    def __init__(self, root, options, threads):
        self.root, self.options, self.threads = Path(root).resolve(), options, threads
        self.workers = {}

    def generate(self, reference, count, seed, directory, workers):
        from ..generate_conformers import GenerationError, MissingGeneratorDependency
        directory = Path(directory).resolve()
        command = self.options.get('command')
        if not isinstance(command, list) or not command or not all(isinstance(x, str) for x in command):
            raise MissingGeneratorDependency('Persistent model requires a list-valued generation command')
        slots = min(workers, count)
        raw, errors = [], []
        with ThreadPoolExecutor(max_workers=slots) as pool:
            futures = [pool.submit(copy_context().run, self._request_batch,
                                   i, reference, count, seed, directory, slots) for i in range(slots)]
            for future in futures:
                try:
                    raw.extend(future.result())
                except (GenerationError, OSError, ValueError) as exc:
                    errors.append(str(exc))
        if len(errors) == slots:
            raise GenerationError(f'All persistent workers failed: {errors}')
        (directory/'worker_errors.json').write_text(json.dumps(errors))
        return raw

    def _request_batch(self, index, reference, count, seed, directory, slots):
        """Start a worker if needed, then request a batch from the same process."""
        from ..generate_conformers import GenerationError
        folder = directory/f'worker_{index:02d}'
        folder.mkdir()
        size = count//slots + (index < count % slots)
        with Chem.SDWriter(str(folder/'reference.sdf')) as writer:
            writer.write(reference)
        data = {'smiles': Chem.MolToSmiles(Chem.RemoveHs(reference), isomericSmiles=True),
                'maximum_raw_candidates': size, 'seed': (seed+index) % 2**31,
                'threads': self.threads, 'device': self.options.get('device', 'cpu'), 'optimization': 'none'}
        request_file, output = folder/'request.json', folder/'external_raw.sdf'
        request_file.write_text(json.dumps(data))
        if index not in self.workers:
            self._start_worker(index, folder, request_file, output)
        process, session, _ = self.workers[index]
        response = folder/'response.json'
        temporary = session/'next.tmp'
        temporary.write_text(json.dumps({'request': str(request_file), 'output': str(output), 'response': str(response)}))
        temporary.replace(session/'next.json')
        limit = self.options.get('timeout_seconds', 600)
        remaining = job_timeout.remaining_time()
        if remaining is not None:
            limit = min(limit, remaining)
        deadline = time.monotonic()+limit
        while not response.is_file():
            if process.poll() is not None:
                raise GenerationError(f'Persistent worker exited ({process.returncode}); see {session}/worker.out')
            if time.monotonic() >= deadline:
                raise GenerationError(f'Persistent generation timed out; see {session}/worker.out')
            time.sleep(0.05)
        reply = json.loads(response.read_text())
        if reply.get('error'):
            raise GenerationError(reply['error'])
        if not output.is_file():
            raise GenerationError('Persistent worker produced no raw SDF')
        records = list(Chem.SDMolSupplier(str(output), removeHs=False)) if output.read_text().strip() else []
        if len(records) > size:
            raise GenerationError('Persistent worker exceeded requested count')
        return records

    def _start_worker(self, index, folder, request_file, output):
        """Load a model once; subsequent batches reuse this worker."""
        session = self.root/f'session_worker_{index:02d}'
        session.mkdir()
        values = {'request': str(request_file), 'output': str(output),
                  'input': str(folder/'reference.sdf')}
        argv = [token.format(**values) for token in self.options['command']]
        argv += ['--session', str(session)]
        env = os.environ.copy()
        for key in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
            env[key] = str(self.threads)
        log = (session/'worker.out').open('w')
        try:
            process = subprocess.Popen(argv, cwd=folder, env=env, stdout=log,
                stderr=subprocess.STDOUT, start_new_session=True)
        except Exception:
            log.close()
            raise
        self.workers[index] = (process, session, log)

    def close(self):
        for process, session, log in self.workers.values():
            try:
                (session/'stop').touch()
                try:
                    process.wait(timeout=2)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGTERM)
                    try:
                        process.wait(timeout=2)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGKILL)
                        process.wait()
            finally:
                log.close()
