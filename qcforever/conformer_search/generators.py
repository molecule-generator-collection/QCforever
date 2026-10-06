"""Generation adapters with an explicit fixed raw-attempt cap."""
import json
import os
from concurrent.futures import ThreadPoolExecutor
from contextvars import copy_context
from pathlib import Path
import subprocess

from rdkit import Chem
from rdkit.Chem import AllChem


class GenerationError(RuntimeError):
    pass


class MissingGeneratorDependency(GenerationError):
    pass


def split_conformers(mol):
    records = []
    for conf in mol.GetConformers():
        record = Chem.Mol(mol)
        record.RemoveAllConformers()
        record.AddConformer(Chem.Conformer(conf), assignId=True)
        records.append(record)
    return records


def etkdgv3(reference, attempt_cap, seed, directory, options, threads):
    mol = Chem.AddHs(Chem.RemoveHs(reference))
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.numThreads = threads
    params.pruneRmsThresh = -1.0
    params.enforceChirality = True
    AllChem.EmbedMultipleConfs(mol, numConfs=attempt_cap, params=params)
    return split_conformers(mol)


def etflow(reference, attempt_cap, seed, directory, options, threads):
    try:
        from etflow import BaseFlow
        import torch
    except ImportError as exc:
        raise MissingGeneratorDependency('ET-Flow is not installed in this Python environment') from exc
    if not options.get('checkpoint_cache'):
        raise MissingGeneratorDependency('ET-Flow requires an explicit checkpoint_cache')
    device = options.get('device', 'cpu')
    if device == 'cpu':
        torch.set_num_threads(threads)
    try:
        model = BaseFlow.from_default(model=options.get('model', 'drugs-so3'),
                                      cache=str(Path(options['checkpoint_cache']).expanduser().resolve()))
        model = model.to(device).eval()
        smiles = Chem.MolToSmiles(Chem.RemoveHs(reference), isomericSmiles=True)
        generated = model.predict([smiles], max_batch_size=options.get('batch_size', 4),
                                  num_samples=attempt_cap, seed=seed, device=device, as_mol=True)
        value = generated[smiles]
        return [Chem.Mol(m) for m in value] if isinstance(value, list) else split_conformers(value)
    except (RuntimeError, ValueError, KeyError, OSError) as exc:
        raise GenerationError(f'ET-Flow generation failed: {exc}') from exc


def command_generator(reference, attempt_cap, seed, directory, options, threads):
    """Optional model environment writes raw SDF; it must not filter or optimize.

    Command tokens may contain {request}, {output}, {input}. No shell expansion.
    The caller must supply a generator command; unconfigured models are reported
    as unavailable rather than being replaced behind the caller's back.
    """
    command = options.get('command')
    if not isinstance(command, list) or not command or not all(isinstance(x, str) for x in command):
        raise MissingGeneratorDependency('An explicit list-valued generator command is required')
    request = directory/'request.json'
    output = directory/'external_raw.sdf'
    inp = directory/'reference.sdf'
    with Chem.SDWriter(str(inp)) as writer:
        writer.write(reference)
    request.write_text(json.dumps({'smiles': Chem.MolToSmiles(Chem.RemoveHs(reference), isomericSmiles=True),
                                   'maximum_raw_candidates': attempt_cap, 'seed': seed,
                                   'threads': threads, 'optimization': 'none',
                                   'options': {k:v for k,v in options.items() if k not in ('command', 'timeout_seconds')}}, indent=2))
    values = {'request': str(request), 'output': str(output), 'input': str(inp)}
    argv = [token.format(**values) for token in command]
    try:
        with (directory/'generator.out').open('w') as out:
            from qcforever.util import job_timeout
            env = os.environ.copy()
            for key in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
                env[key] = str(threads)
            job_timeout.run(argv, cwd=directory, env=env, stdout=out, stderr=subprocess.STDOUT,
                           check=True, timeout=options.get('timeout_seconds', 600))
    except (subprocess.SubprocessError, OSError) as exc:
        raise GenerationError(f'Generator command failed: {exc}') from exc
    if not output.is_file():
        raise GenerationError('Generator command did not produce the requested raw SDF')
    # A successful model command may legitimately generate zero candidates.
    # RDKit raises OSError for an empty SDF; return an empty pool explicitly so
    # the fixed attempt cap/fallback policy can handle it, not abort the case.
    if not output.read_text().strip():
        return []
    try:
        return list(Chem.SDMolSupplier(str(output), removeHs=False))
    except (OSError, ValueError) as exc:
        raise GenerationError(f'Unreadable model raw SDF: {exc}') from exc


def generate(name, reference, attempt_cap, seed, directory, options, threads):
    if options.get('command'):
        return command_generator(reference, attempt_cap, seed, directory, options, threads)
    if name == 'etkdgv3':
        return etkdgv3(reference, attempt_cap, seed, directory, options, threads)
    if name == 'etflow':
        return etflow(reference, attempt_cap, seed, directory, options, threads)
    if name in ('ditmc', 'torsional_diffusion'):
        raise MissingGeneratorDependency(f'{name} requires an optional environment and raw-SDF command')
    raise ValueError(f'Unknown generator: {name}')


def generate_batch(name, reference, count, seed, directory, options, threads, workers):
    """Parallel model subprocesses; deterministic merge independent of completion order.

    Threads only supervise isolated command processes. No model is loaded in the
    QCforever interpreter. Native ETKDG uses its own bounded RDKit thread pool.
    """
    if name == 'etkdgv3':
        return generate(name, reference, count, seed, directory, options, threads*workers)
    if not options.get('command'):
        raise MissingGeneratorDependency(f'{name}: configure a raw-SDF command')
    slots = min(workers, count)
    requests = []
    for i in range(slots):
        folder = directory/f'worker_{i:02d}'
        folder.mkdir()
        size = count//slots + (i < count % slots)
        requests.append((reference, size, (seed+i) % 2**31, folder, options, threads))
    with ThreadPoolExecutor(max_workers=slots) as pool:
        futures = [pool.submit(copy_context().run, command_generator, *args) for args in requests]
        raw = []
        failures = []
        for future, args in zip(futures, requests):
            try:
                records = future.result()
                if len(records) > args[1]:
                    raise GenerationError('Worker returned more records than requested')
                raw.extend(records)
            except GenerationError as exc:
                failures.append(str(exc))
                (args[3]/'error.json').write_text(json.dumps({'error': str(exc)}))
        if len(failures) == slots:
            raise GenerationError(f'All generation workers failed: {failures}')
        return raw
