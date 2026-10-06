"""Isolated adapter around the pinned upstream TD sample_confs function.

Only bootstrap the loader with an empty CSV; subsequent requests reuse weights.
MM/energy/likelihood work is explicitly disabled. Changes are process-local,
not edits to the upstream checkout. Cache/output files stay outside the checkout.
"""
import csv
import os
import runpy
import sys
import types


def initialize(args, request):
    # Upstream lookup tables write relative cache files during import. Serialize
    # first construction when multiple workers share a cache directory.
    import fcntl
    args.cache.mkdir(parents=True, exist_ok=True)
    with (args.cache/'initialize.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        os.chdir(args.cache)
        return _load(args, request)


def _load(args, request):
    folder = args.output.parent
    import torch
    torch.set_num_threads(request['threads'])
    # Execute the official loader once with no molecules. The sample_confs
    # function retains its model/args/imported lookup tables for later calls.
    inp = folder/'td_empty.csv'
    with inp.open('w', newline='') as handle:
        csv.writer(handle).writerow(['smiles', 'num_conformers', 'corrected_smiles'])
    sys.modules['utils.xtb'] = types.ModuleType('utils.xtb')
    sys.path.insert(0, str(args.source.resolve()))
    import diffusion.sampling as sampling
    sampling.populate_likelihood = lambda *a, **k: None
    sys.argv = ['generate_confs.py', '--test_csv', str(inp), '--inference_steps', '20',
                '--model_dir', str(args.checkpoint.resolve()), '--batch_size', '1', '--no_energy']
    namespace = runpy.run_path(str(args.source.resolve()/'generate_confs.py'), run_name='__main__')
    sample = namespace['sample_confs']
    globals_ = sample.__globals__
    def embed(mol, numConfs):
        globals_['AllChem'].EmbedMultipleConfs(mol, numConfs=numConfs, numThreads=request['threads'])
        return mol
    globals_['embed_func'] = embed
    return sample
