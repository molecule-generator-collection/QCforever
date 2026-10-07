"""Generation must respect nproc after the existing runner adjusts it."""
import json

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from qcforever.conformer_search.config import SearchConfig
from qcforever.conformer_search import pipeline


@pytest.mark.parametrize('cores,threads,workers', [
    (1, 1, 1), (2, 2, 1), (3, 3, 1), (4, 4, 1),
    (5, 4, 1), (8, 4, 2), (16, 4, 4), (64, 4, 8),
])
def test_default_generation_cpu_limit(cores, threads, workers):
    config = SearchConfig()
    assert config.effective_threads(cores) == threads
    assert config.parallelism(cores) == workers
    assert threads * workers <= cores
    assert config.threads == 4  # Do not mutate requested settings.


@pytest.mark.parametrize('cores', [0, -1, True, 1.5])
def test_invalid_cpu_limit(cores):
    with pytest.raises(ValueError, match='positive integer'):
        SearchConfig().parallelism(cores)


@pytest.mark.parametrize('engine', ['gaussian', 'gamess'])
def test_runner_reduced_nproc_can_generate(tmp_path, monkeypatch, engine):
    from qcforever.util import check_resource
    from qcforever.gaussian_run.GaussianRunPack import GaussianDFTRun
    from qcforever.gamess_run.GamessRunPack import GamessDFTRun
    # Reproduce the existing adjustment, not a change to that legacy policy.
    monkeypatch.setattr(check_resource.psutil, 'cpu_count', lambda: 4)
    monkeypatch.setattr(check_resource.psutil, 'cpu_percent', lambda **kwargs: 1.0)
    runner = GaussianDFTRun if engine == 'gaussian' else GamessDFTRun
    job = runner('B3LYP', '6-31G*', 4, 'optconf=xtb opt energy', 'molecule.sdf')
    assert job.nproc == 3
    result = pipeline.prepare_candidates('CCO', tmp_path/'run', SearchConfig(),
                                         allocated_cores=job.nproc)
    assert result.candidates
    assert json.loads(result.status_path.read_text())['threads_per_worker'] == 3


@pytest.mark.parametrize('profile', ['low', 'medium', 'high'])
@pytest.mark.parametrize('cores', [1, 2, 3])
def test_real_etkdg_with_reduced_nproc(tmp_path, profile, cores):
    # Medium/high deliberately exercise fallback here, not real model inference.
    config = SearchConfig(profile=profile, mm_method='none', device='cpu')
    result = pipeline.prepare_candidates('CCO', tmp_path/'run', config,
        allocated_cores=cores,
        generators={'ditmc': lambda *args: [], 'torsional_diffusion': lambda *args: []})
    assert result.candidates
    status = json.loads(result.status_path.read_text())
    assert status['threads_per_worker'] == cores
    assert status['requested_threads_per_worker'] == 4
    assert status['allocated_cores'] == cores
    assert status['workers'] == 1
    saved = json.loads((tmp_path/'run'/'resolved_config.json').read_text())
    assert saved['threads'] == 4


@pytest.mark.parametrize('name,profile', [('ditmc', 'high'), ('torsional_diffusion', 'medium')])
@pytest.mark.parametrize('persistent', [False, True])
@pytest.mark.parametrize('device', ['cpu', 'gpu'])
def test_model_paths_receive_effective_threads(tmp_path, monkeypatch, name, profile, persistent, device):
    from qcforever.conformer_search import model_session
    calls = []
    reference = Chem.AddHs(Chem.MolFromSmiles('CCO'))
    assert AllChem.EmbedMolecule(reference, randomSeed=42) == 0

    def probe(options, requested, directory, threads):
        assert threads == 3
        return {'requested': requested, 'effective': device}

    def generate(model_name, mol, count, seed, folder, options, threads, workers):
        assert model_name == name and threads == 3 and workers == 1
        calls.append(count)
        return [Chem.Mol(reference)]

    class Session:
        def __init__(self, directory, options, threads):
            assert threads == 3

        def generate(self, mol, count, seed, folder, workers):
            assert workers == 1
            calls.append(count)
            return [Chem.Mol(reference)]

        def close(self):
            pass

    monkeypatch.setattr(pipeline, 'resolve_model_device', probe)
    monkeypatch.setattr(pipeline, 'generate_batch', generate)
    monkeypatch.setattr(model_session, 'ModelSession', Session)
    config = SearchConfig.from_mapping({'profile': profile, 'mm_method': 'none',
        'budget': {'formula': 'fixed', 'fixed': 2},
        'validation': {'duplicate_rmsd_angstrom': 0},
        'generators': {name: {'persistent': persistent}}})
    result = pipeline.prepare_candidates(reference, tmp_path/'run', config, allocated_cores=3)
    assert len(result.candidates) == 2
    assert calls == [2, 1]  # First request and additional batch both work.
