"""DiTMC Drugs aPE-B inference using the upstream model API.

Extracted from the previously validated comparison adapter. No benchmark project
imports, trained weights, or upstream source code are bundled here.
"""
from pathlib import Path
import sys
import time
import numpy as np
from rdkit import Chem


def initialize_ditmc_extension(root, source):
    """Serialize upstream Cython compilation before concurrent model workers.

    The extension code and compiler options are upstream defaults. A campaign
    cache avoids cross-job races in the shared home .pyxbld directory.
    """
    import fcntl
    import importlib
    import pyximport
    cache = Path(root) / 'initialization'
    cache.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(Path(source).resolve()))
    with (cache/'ditmc_algos.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        pyximport.install(build_dir=str(cache/'pyxbld'), setup_args={'include_dirs':np.get_include()})
        module = importlib.import_module('dit_mc.algos')
        fcntl.flock(lock, fcntl.LOCK_UN)
    return module.__file__


class DiTMC:
    def __init__(self, root, spec, batch_size):
        self.source = (root / spec['source']).resolve()
        initialize_ditmc_extension(root, self.source)
        sys.path.insert(0, str(self.source))
        import hydra
        from omegaconf import OmegaConf
        from dit_mc.training.checkpoint import create_checkpoint_manager_from_workdir
        import jax
        from orbax import checkpoint as ocp
        configs = list((root / spec['checkpoint_root']).glob('**/drugs/apeB/.hydra/config.yaml'))
        if len(configs) != 1:
            raise ValueError(f'Expected one Drugs aPE-B config, found {len(configs)}')
        cfg = OmegaConf.load(configs[0])
        self.process = hydra.utils.instantiate(cfg.generative_process)
        workdir = configs[0].parent.parent
        folders = [d for d in ('last_checkpoint', 'checkpoints') if (workdir / d).is_dir()]
        if len(folders) != 1:
            raise ValueError(f'Ambiguous checkpoint folders: {folders}')
        # Restore weights onto the active device, not the training-time cuda:0.
        # No checkpoint files or tensor values are modified.
        manager = create_checkpoint_manager_from_workdir(str(workdir), ckpt_dir_name=folders[0], create_ckpt_dir=False)
        try:
            step = manager.latest_step()
            metadata = manager.item_metadata(step)['params']
            sharding = jax.sharding.SingleDeviceSharding(jax.devices()[0])
            abstract = jax.tree.map(lambda x: jax.ShapeDtypeStruct(x.shape, x.dtype, sharding=sharding), metadata)
            restored = manager.restore(step, args=ocp.args.Composite(params=ocp.args.StandardRestore(abstract)))
            self.params = restored['params']
        finally:
            manager.close()
        self.batch_size, self.steps = batch_size, spec['steps']
        # Official cutoff=inf is valid YAML but is not a finite JSON number.
        self.metadata = dict(config_yaml=OmegaConf.to_yaml(cfg, resolve=True), step=step,
                             checkpoint_directory=str(workdir))

    def generate(self, smiles, count, seed, start_index=0):
        import jax
        import jraph
        import tensorflow as tf
        from dit_mc.prepare_dataset import get_graph_row_from_mol, get_prior_graph_row_from_mol
        from dit_mc.data_loader.utils import create_graph_tuples
        from dit_mc.generate_confs import get_chiral_centers_from_mol, switch_parity_of_pos, set_rdmol_positions
        # Use official graph features/prior construction, without reference coordinates.
        sys.path.insert(0, str(self.source / 'tf_datasets' / 'geom'))
        from preprocessing import rows_to_sample, cutoff_graph_to_bond_graph
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        mol.AddConformer(Chem.Conformer(mol.GetNumAtoms()))  # zero placeholder, not an initial conformer
        graph, _ = get_graph_row_from_mol(mol, smiles, 'drugs', {})
        prior, _, _ = get_prior_graph_row_from_mol(mol, smiles, {}, {})
        sample = rows_to_sample(graph, prior, cutoff_graph_to_bond_graph(graph), smiles, smiles, -1, False)
        tuples = create_graph_tuples({k: tf.convert_to_tensor(v) for k, v in sample.items()},
                                      cutoff=float('inf'), split='train', to_numpy=True)
        # jraph statistics require an explicit padding graph, as in upstream batching.
        def padded(g):
            bg = jraph.batch([g] * self.batch_size)
            return jraph.pad_with_graphs(bg, n_node=int(sum(bg.n_node)) + 1,
                                         n_edge=int(sum(bg.n_edge)) + 1, n_graph=self.batch_size + 1)
        latent, cond, prior = [padded(g) for g in tuples]
        _, nbr, tag = get_chiral_centers_from_mol(mol)
        raw, corrected = [], []
        self.last_batch_seconds = []
        key = jax.random.PRNGKey(seed)
        if start_index and self.batch_size != 1:
            raise ValueError('Random-access generation requires batch_size=1')
        for _ in range(start_index):
            key, _ = jax.random.split(key)
        for _ in range((count + self.batch_size - 1) // self.batch_size):
            started = time.perf_counter()
            key, sample_key = jax.random.split(key)
            prediction = self.process.sample(self.params, latent, prior, cond, sample_key, num_steps=self.steps)
            coords = np.asarray(prediction.nodes['positions'])[:self.batch_size * mol.GetNumAtoms()]
            self.last_batch_seconds.append(time.perf_counter() - started)
            for xyz in coords.reshape(self.batch_size, mol.GetNumAtoms(), 3):
                raw.append(set_rdmol_positions(mol, xyz))
                corrected.append(set_rdmol_positions(mol, np.asarray(switch_parity_of_pos(xyz, nbr, tag))))
                if len(raw) == count:
                    return raw, corrected
        return raw, corrected
