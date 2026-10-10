"""Isolated adapter around the pinned upstream TD sample_confs function.

Only bootstrap the loader with an empty CSV; subsequent requests reuse weights.
MM/energy/likelihood work is explicitly disabled. Changes are process-local,
not edits to the upstream checkout. Cache/output files stay outside the checkout.
"""
import csv
import importlib
import os
import runpy
import sys
import time
import types


def initialize_model(args, request):
    """Prepare shared tables, then load each worker's model independently."""
    args.cache.mkdir(parents=True, exist_ok=True)
    os.chdir(args.cache)
    sys.path.insert(0, str(args.source.resolve()))
    args.initialization_timings = _initialize_torsion_tables(args.cache)
    started = time.perf_counter()
    model = _load_model(args, request)
    args.initialization_timings['model_load_seconds'] = time.perf_counter() - started
    return model


def _initialize_torsion_tables(cache):
    """Only the first writer needs exclusivity; existing tables allow readers.

    The pinned upstream diffusion.torus writes .p.npy followed by .score.npy
    on import. Keep its computations unchanged and use the same lock filename
    as older workers, so a new reader cannot read an old writer's partial file.
    """
    import fcntl
    tables = [cache/'.p.npy', cache/'.score.npy']
    started = time.perf_counter()
    with (cache/'initialize.lock').open('a') as lock:
        ready = all(path.is_file() for path in tables)
        fcntl.flock(lock, fcntl.LOCK_SH if ready else fcntl.LOCK_EX)
        # Another worker may have finished while this worker waited.
        if all(path.is_file() for path in tables):
            fcntl.flock(lock, fcntl.LOCK_SH)
        elif any(path.exists() for path in tables):
            raise RuntimeError(f'Incomplete TD lookup tables in {cache}; use a fresh cache directory')
        acquired = time.perf_counter()
        importlib.import_module('diffusion.torus')
    return dict(cache_lock_wait_seconds=acquired-started,
                table_initialization_seconds=time.perf_counter()-acquired)


def generate_conformers(sample, smiles, count, seed, threads):
    """Seed RDKit embedding for this request, including persistent extra batches.

    Python/NumPy/Torch RNGs are seeded by the worker. RDKit has a separate RNG;
    do not rely on those seeds or retain the first batch's seed in a closure.
    """
    globals_ = sample.__globals__
    def embed(mol, numConfs):
        globals_['AllChem'].EmbedMultipleConfs(mol, numConfs=numConfs,
                                             numThreads=threads, randomSeed=seed)
        return mol
    globals_['embed_func'] = embed
    return sample(smiles, count, smiles)


def _load_model(args, request):
    folder = args.output.parent
    import torch
    torch.set_num_threads(request['threads'])
    # Execute the official loader once with no molecules. The sample_confs
    # function retains its model/args/imported lookup tables for later calls.
    inp = folder/'td_empty.csv'
    with inp.open('w', newline='') as handle:
        csv.writer(handle).writerow(['smiles', 'num_conformers', 'corrected_smiles'])
    sys.modules['utils.xtb'] = types.ModuleType('utils.xtb')
    import diffusion.sampling as sampling
    sampling.populate_likelihood = lambda *a, **k: None
    sys.argv = ['generate_confs.py', '--test_csv', str(inp), '--inference_steps', '20',
                '--model_dir', str(args.checkpoint.resolve()), '--batch_size', '1', '--no_energy']
    namespace = runpy.run_path(str(args.source.resolve()/'generate_confs.py'), run_name='__main__')
    return namespace['sample_confs']
