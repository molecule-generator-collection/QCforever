"""Real generation/reuse check, run inside an isolated model environment.

Uses only the worker's dependencies, not QCforever's Gaussian/GAMESS imports.
No xTB/PM6 calculation is performed and no benchmark result is claimed.
"""
import argparse
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

from .installed_models import write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--options', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--device', choices=('cpu', 'gpu'), required=True)
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--timeout', type=int, default=1800)
    args = parser.parse_args()
    check(json.loads(args.options.read_text()), args.output, args.device, args.threads, args.timeout)


def check(options, root, device, threads, timeout):
    root = Path(root).resolve()
    root.mkdir(parents=True, exist_ok=False)
    session = root/'session'
    session.mkdir()
    process, rows = None, []
    with (root/'worker.log').open('w') as log:
        try:
            for number in range(2):
                folder = root/f'batch_{number}'
                folder.mkdir()
                request, output = folder/'request.json', folder/'raw.sdf'
                write_json(request, {'smiles': 'CCO', 'maximum_raw_candidates': 1,
                                    'seed': 20261006+number, 'threads': threads,
                                    'device': device, 'optimization': 'none'})
                if process is None:
                    argv = [s.format(request=str(request), output=str(output), input='')
                            for s in options['command']]
                    env = os.environ.copy()
                    for key in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
                        env[key] = str(threads)
                    process = subprocess.Popen(argv+['--session', str(session)], cwd=root,
                        env=env, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
                response = folder/'response.json'
                write_json(session/'next.json', {'request': str(request), 'output': str(output),
                                                 'response': str(response)})
                deadline = time.monotonic()+timeout
                while not response.exists():
                    if process.poll() is not None:
                        raise RuntimeError(f'Worker exited {process.returncode}; see {root}/worker.log')
                    if time.monotonic() >= deadline:
                        raise RuntimeError(f'Smoke test timed out; see {root}/worker.log')
                    time.sleep(0.1)
                reply = json.loads(response.read_text())
                if reply.get('error'):
                    raise RuntimeError(reply['error'])
                verify_records(output)
                rows.append(json.loads((folder/'model_execution.json').read_text()))
            if rows[0]['worker_pid'] != rows[1]['worker_pid'] or not rows[1]['model_reused']:
                raise RuntimeError('Additional batch did not reuse the model process')
            if rows[1]['initialization_seconds'] != 0:
                raise RuntimeError('Additional batch reinitialized the model')
            if any(row['device']['effective'] != device for row in rows):
                raise RuntimeError('Worker used a different device than the setup target')
            result = {'state': 'passed', 'molecule': 'CCO', 'requested_batches': [1, 1],
                      'device': device, 'rows': rows}
            write_json(root/'summary.json', result)
            return result
        finally:
            if process is not None:
                (session/'stop').touch()
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGTERM)
                    try:
                        process.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGKILL)
                        process.wait()


def verify_records(path):
    import numpy as np
    from rdkit import Chem
    records = list(Chem.SDMolSupplier(str(path), removeHs=False))
    if len(records) != 1 or records[0] is None:
        raise RuntimeError('Smoke test did not return one readable conformer')
    mol = records[0]
    if mol.GetNumConformers() != 1 or mol.GetNumAtoms() != 9:
        raise RuntimeError('Smoke test ethanol has wrong atom/conformer count')
    if Chem.MolToSmiles(Chem.RemoveHs(mol)) != 'CCO':
        raise RuntimeError('Smoke test changed ethanol connectivity')
    xyz = mol.GetConformer().GetPositions()
    if not np.isfinite(xyz).all():
        raise RuntimeError('Smoke test contains nonfinite coordinates')
    # Reject collapsed structures, not just an existing SDF file.
    distances = np.linalg.norm(xyz[:, None, :]-xyz[None, :, :], axis=-1)
    if (distances[np.triu_indices(len(xyz), 1)] < 0.25).any():
        raise RuntimeError('Smoke test contains atom collisions')




if __name__ == '__main__':
    main()
