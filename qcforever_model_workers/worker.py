"""Generation-only worker: persistent sessions reuse model and lookup tables.

No validation, MM or semi-empirical optimization. Official DiTMC correction
is retained; TD reuses the official sample_confs function with MM/energy off.
"""
import argparse
import json
import os
from pathlib import Path
import random
import sys
import time


def write_sdf(path, records):
    from rdkit import Chem
    with Chem.SDWriter(str(path)) as writer:
        for mol in records:
            writer.write(mol)


def initialize(args, request):
    folder = args.output.resolve().parent
    if hasattr(os, 'sched_getaffinity'):
        available = sorted(os.sched_getaffinity(0))
        index = int(folder.name.split('_')[-1]) if folder.name.startswith('worker_') else 0
        threads = request['threads']
        start = index*threads
        if len(available) < start+threads:
            raise ValueError('Insufficient scheduler-assigned CPUs for model worker')
        os.sched_setaffinity(0, available[start:start+threads])
    os.environ.setdefault('CUDA_VISIBLE_DEVICES', '')
    os.environ.setdefault('JAX_PLATFORMS', 'cpu')
    if args.model == 'ditmc':
        from .ditmc import DiTMC
        return DiTMC(args.cache, {'source': str(args.source),
            'checkpoint_root': str(args.checkpoint), 'steps': 50}, batch_size=1)
    from .torsional import initialize as initialize_torsional
    return initialize_torsional(args, request)


def generate(args, model, message, initialization_seconds, number):
    import numpy as np
    from rdkit import Chem
    request = json.loads(Path(message['request']).read_text())
    count, seed = request['maximum_raw_candidates'], request['seed']
    folder = Path(message['output']).parent
    started = time.perf_counter()
    if args.model == 'ditmc':
        raw, records = model.generate(request['smiles'], count, seed)
        write_sdf(folder/'pre_correction_raw.sdf', raw)
    else:
        import torch
        random.seed(seed)
        np.random.seed(seed % 2**32)
        torch.manual_seed(seed)
        records = model(request['smiles'], count, request['smiles']) or []
        records = [Chem.Mol(m) for m in records]
    seconds = time.perf_counter()-started
    if len(records) > count:
        raise ValueError('Model returned more than requested raw candidates')
    write_sdf(message['output'], records)
    (folder/'model_execution.json').write_text(json.dumps({
        'model': args.model, 'requested': count, 'returned': len(records), 'seed': seed,
        'worker_pid': os.getpid(), 'request_number': number,
        'model_initialization_count': 1, 'model_reused': number > 1,
        'initialization_seconds': initialization_seconds if number == 1 else 0.0,
        'generation_seconds': seconds, 'source': str(args.source), 'checkpoint': str(args.checkpoint),
        'filtering': False, 'force_field_optimization': False, 'semiempirical_relaxation': False,
        'cpu_affinity': sorted(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else None}, indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--model', choices=['ditmc', 'torsional_diffusion'], required=True)
    parser.add_argument('--request', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--checkpoint', type=Path, required=True)
    parser.add_argument('--legacy-workflow', type=Path, help=argparse.SUPPRESS)
    parser.add_argument('--cache', type=Path)
    parser.add_argument('--session', type=Path)
    args = parser.parse_args()
    # Resolve all caller paths BEFORE an upstream loader changes directory.
    for name in ('request', 'output', 'source', 'checkpoint', 'session', 'cache'):
        value = getattr(args, name)
        if value is not None:
            setattr(args, name, value.expanduser().resolve())
    if args.cache is None:
        args.cache = args.output.parent / 'model_cache'
    for name in ('source', 'checkpoint'):
        if not getattr(args, name).is_dir():
            parser.error(f'--{name} must be an existing directory')
    if args.legacy_workflow:
        print('--legacy-workflow is deprecated and ignored; using bundled adapter', flush=True)
    request = json.loads(args.request.read_text())
    started = time.perf_counter()
    model = initialize(args, request)
    elapsed = time.perf_counter()-started
    if args.session is None:
        generate(args, model, {'request': str(args.request), 'output': str(args.output)}, elapsed, 1)
        return
    number = 0
    while not (args.session/'stop').exists():
        path = args.session/'next.json'
        if not path.exists():
            time.sleep(0.05)
            continue
        message = json.loads(path.read_text())
        path.unlink()
        number += 1
        try:
            generate(args, model, message, elapsed, number)
            response = {'state': 'completed'}
        except Exception as exc:
            response = {'error': f'{type(exc).__name__}: {exc}'}
        target = Path(message['response'])
        temporary = target.with_suffix('.tmp')
        temporary.write_text(json.dumps(response))
        temporary.replace(target)


if __name__ == '__main__':
    main()
