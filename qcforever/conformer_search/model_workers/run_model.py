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


def main():
    """Load one model, then serve one request or a persistent request loop."""
    args = _parse_arguments()
    request = json.loads(args.request.read_text())
    configure_worker_cpus(args, request)
    args.device_info = select_device(args.model, request.get('device', 'cpu'))
    if args.probe_device:
        args.output.write_text(json.dumps(args.device_info, indent=2))
        return
    started = time.perf_counter()
    model = initialize_model(args, request)
    args.parameter_devices = parameter_devices(args.model, model)
    if args.device_info['effective'] == 'gpu' and (not args.parameter_devices or not all(
            'cuda' in d.lower() or 'gpu' in d.lower() for d in args.parameter_devices)):
        raise RuntimeError(f'Model parameters are not on GPU: {args.parameter_devices}')
    elapsed = time.perf_counter()-started
    if args.session is None:
        generate_batch(args, model, {'request': str(args.request), 'output': str(args.output)}, elapsed, 1)
        return
    _serve_requests(args, model, elapsed)


def _parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--model', choices=['ditmc', 'torsional_diffusion'], required=True)
    parser.add_argument('--request', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--checkpoint', type=Path, required=True)
    parser.add_argument('--legacy-workflow', type=Path, help=argparse.SUPPRESS)
    parser.add_argument('--cache', type=Path)
    parser.add_argument('--session', type=Path)
    parser.add_argument('--probe-device', action='store_true',
                        help='Write device availability JSON without loading model weights')
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

    return args


def initialize_model(args, request):
    """Load the selected model once, after CPU/GPU selection."""
    if args.model == 'ditmc':
        from .ditmc import DiTMC
        return DiTMC(args.cache, {'source': str(args.source),
            'checkpoint_root': str(args.checkpoint), 'steps': 50}, batch_size=1)
    from .torsional_diffusion import initialize_model as initialize_torsional
    return initialize_torsional(args, request)


def _serve_requests(args, model, initialization_seconds):
    """Reuse model parameters until the parent requests a clean shutdown."""
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
            generate_batch(args, model, message, initialization_seconds, number)
            response = {'state': 'completed'}
        except Exception as exc:
            response = {'error': f'{type(exc).__name__}: {exc}'}
        target = Path(message['response'])
        temporary = target.with_suffix('.tmp')
        temporary.write_text(json.dumps(response))
        temporary.replace(target)


def generate_batch(args, model, message, initialization_seconds, number):
    """Generate unfiltered candidates and record the model execution metadata."""
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
        from .torsional_diffusion import generate_conformers as generate_torsional
        random.seed(seed)
        np.random.seed(seed % 2**32)
        torch.manual_seed(seed)
        records = generate_torsional(model, request['smiles'], count, seed, request['threads']) or []
        records = [Chem.Mol(m) for m in records]
    seconds = time.perf_counter()-started
    if len(records) > count:
        raise ValueError('Model returned more than requested raw candidates')
    write_sdf(message['output'], records)
    (folder/'model_execution.json').write_text(json.dumps({
        'model': args.model, 'requested': count, 'returned': len(records), 'seed': seed,
        'rdkit_embedding_seed': seed if args.model == 'torsional_diffusion' else None,
        'worker_pid': os.getpid(), 'request_number': number,
        'model_initialization_count': 1, 'model_reused': number > 1,
        'initialization_seconds': initialization_seconds if number == 1 else 0.0,
        'generation_seconds': seconds, 'source': str(args.source), 'checkpoint': str(args.checkpoint),
        'device': args.device_info, 'parameter_devices': args.parameter_devices,
        'filtering': False, 'force_field_optimization': False, 'semiempirical_relaxation': False,
        'cpu_affinity': sorted(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else None}, indent=2))



def select_device(model_name, requested):
    """Resolve in the model environment; never unmask scheduler-hidden GPUs."""
    if requested not in ('auto', 'cpu', 'gpu'):
        raise ValueError('device must be auto, cpu, or gpu')
    masked = os.environ.get('CUDA_VISIBLE_DEVICES') in ('', '-1')
    if requested == 'gpu' and masked:
        raise RuntimeError('GPU requested but hidden by CUDA_VISIBLE_DEVICES')
    if requested == 'cpu' or (requested == 'auto' and masked):
        os.environ['CUDA_VISIBLE_DEVICES'] = ''
        os.environ['JAX_PLATFORMS'] = 'cpu'
        return {'requested': requested, 'effective': 'cpu', 'devices': []}
    os.environ.setdefault('XLA_PYTHON_CLIENT_PREALLOCATE', 'false')
    if model_name == 'ditmc':
        import jax
        try:
            devices = [str(d) for d in jax.devices() if d.platform == 'gpu']
        except RuntimeError as exc:
            # CUDA-enabled JAX can raise rather than return CPU devices on a
            # machine with no GPU. Do not hide driver/library incompatibilities.
            no_device = any(s in str(exc) for s in ('No visible GPU devices', 'CUDA_ERROR_NO_DEVICE'))
            if requested != 'auto' or not no_device:
                raise
            jax.config.update('jax_platforms', 'cpu')
            devices = [str(d) for d in jax.devices() if d.platform == 'gpu']
    else:
        import torch
        devices = ([torch.cuda.get_device_name(i) for i in range(torch.cuda.device_count())]
                   if torch.cuda.is_available() else [])
    if requested == 'gpu' and not devices:
        raise RuntimeError('GPU requested but unavailable in this model environment')
    effective = 'gpu' if devices else 'cpu'
    if effective == 'cpu':
        os.environ['CUDA_VISIBLE_DEVICES'] = ''
    return {'requested': requested, 'effective': effective, 'devices': devices}


def parameter_devices(model_name, model):
    """Record actual parameter placement, independently of availability probes."""
    if model_name == 'ditmc':
        import jax
        return sorted({str(d) for leaf in jax.tree.leaves(model.params) for d in leaf.devices()})
    return sorted({str(p.device) for p in model.__globals__['model'].parameters()})


def configure_worker_cpus(args, request):
    """Bind before framework import, so newly created threads inherit the mask."""
    folder = args.output.resolve().parent
    if hasattr(os, 'sched_getaffinity'):
        available = sorted(os.sched_getaffinity(0))
        index = int(folder.name.split('_')[-1]) if folder.name.startswith('worker_') else 0
        threads = request['threads']
        start = index*threads
        if len(available) < start+threads:
            raise ValueError('Insufficient scheduler-assigned CPUs for model worker')
        os.sched_setaffinity(0, available[start:start+threads])


def write_sdf(path, records):
    from rdkit import Chem
    with Chem.SDWriter(str(path)) as writer:
        for mol in records:
            writer.write(mol)



if __name__ == '__main__':
    main()
