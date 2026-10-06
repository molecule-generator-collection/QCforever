"""Regression tests for shared geometry bounds and readable native failures."""
import json

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from qcforever.conformer_search.config import SearchConfig
from qcforever.conformer_search.validation import geometry_failure
from qcforever.conformer_search.pm6_errors import failure_diagnostic, PM6CalculationError
from qcforever.conformer_search.relaxation import pm6_relax


def test_device_default_and_explicit_cpu_do_not_import_frameworks(monkeypatch):
    from qcforever_model_workers.worker import select_device
    import os
    assert SearchConfig.resolve().device == 'auto'
    with pytest.raises(ValueError, match='device'):
        SearchConfig.resolve(override={'device': 'mystery'})
    monkeypatch.setenv('CUDA_VISIBLE_DEVICES', '3')
    monkeypatch.setenv('JAX_PLATFORMS', 'cuda')
    assert select_device('ditmc', 'cpu')['effective'] == 'cpu'
    assert os.environ['CUDA_VISIBLE_DEVICES'] == ''
    assert os.environ['JAX_PLATFORMS'] == 'cpu'


@pytest.mark.parametrize('available', [False, True])
def test_torch_device_selection_preserves_scheduler_mask(monkeypatch, available):
    import os
    import sys
    from types import SimpleNamespace
    from qcforever_model_workers.worker import select_device
    cuda = SimpleNamespace(is_available=lambda: available, device_count=lambda: 1,
                           get_device_name=lambda i: 'test GPU')
    monkeypatch.setitem(sys.modules, 'torch', SimpleNamespace(cuda=cuda))
    monkeypatch.setenv('CUDA_VISIBLE_DEVICES', '3')
    result = select_device('torsional_diffusion', 'auto')
    assert result['effective'] == ('gpu' if available else 'cpu')
    assert os.environ['CUDA_VISIBLE_DEVICES'] == ('3' if available else '')
    if not available:
        with pytest.raises(RuntimeError, match='GPU requested'):
            select_device('torsional_diffusion', 'gpu')


@pytest.mark.parametrize('platform', ['cpu', 'gpu'])
def test_jax_device_selection_uses_model_framework(monkeypatch, platform):
    import sys
    from types import SimpleNamespace
    from qcforever_model_workers.worker import select_device
    monkeypatch.setitem(sys.modules, 'jax', SimpleNamespace(
        devices=lambda: [SimpleNamespace(platform=platform)]))
    monkeypatch.setenv('CUDA_VISIBLE_DEVICES', '2')
    assert select_device('ditmc', 'auto')['effective'] == platform


@pytest.mark.parametrize('model', ['ditmc', 'torsional_diffusion'])
def test_auto_masked_gpu_selects_cpu_without_importing_framework(monkeypatch, model):
    import sys
    import os
    from qcforever_model_workers.worker import select_device
    monkeypatch.setenv('CUDA_VISIBLE_DEVICES', '')
    monkeypatch.setitem(sys.modules, 'jax', None)
    monkeypatch.setitem(sys.modules, 'torch', None)
    assert select_device(model, 'auto')['effective'] == 'cpu'
    assert os.environ['JAX_PLATFORMS'] == 'cpu'


def test_jax_no_device_falls_back_but_driver_errors_are_not_hidden(monkeypatch):
    import sys
    from types import SimpleNamespace
    from qcforever_model_workers.worker import select_device
    monkeypatch.delenv('CUDA_VISIBLE_DEVICES', raising=False)
    updates = []
    def devices():
        if not updates:
            raise RuntimeError('No visible GPU devices')
        return [SimpleNamespace(platform='cpu')]
    monkeypatch.setitem(sys.modules, 'jax', SimpleNamespace(
        devices=devices, config=SimpleNamespace(update=lambda *a: updates.append(a))))
    assert select_device('ditmc', 'auto')['effective'] == 'cpu'
    assert updates == [('jax_platforms', 'cpu')]
    monkeypatch.delenv('CUDA_VISIBLE_DEVICES', raising=False)
    def driver_error():
        raise RuntimeError('incompatible CUDA driver')
    monkeypatch.setitem(sys.modules, 'jax', SimpleNamespace(devices=driver_error))
    with pytest.raises(RuntimeError, match='incompatible CUDA'):
        select_device('ditmc', 'auto')


def test_device_probe_failure_is_logged_as_stage_failure(tmp_path, monkeypatch):
    from qcforever.conformer_search import pipeline
    from qcforever.conformer_search.generators import GenerationError
    def fail(*args):
        raise GenerationError('probe failure')
    monkeypatch.setattr(pipeline, 'resolve_model_device', fail)
    cfg = SearchConfig.resolve('medium', {'threads': 1, 'mm_method': 'none',
        'budget': {'formula': 'fixed', 'fixed': 1}})
    result = pipeline.prepare_candidates(formaldehyde(), tmp_path/'run', cfg, allocated_cores=1)
    stage = json.loads(result.status_path.read_text())['stages'][0]
    assert stage['state'] == 'generation_failed'
    assert stage['batches'][0]['error'] == 'probe failure'
    assert result.candidates  # The explicit ETKDG fallback still runs.


@pytest.mark.parametrize('effective,expected_workers', [('gpu', 1), ('cpu', 2)])
def test_pipeline_auto_device_controls_worker_count(tmp_path, monkeypatch, effective, expected_workers):
    from qcforever.conformer_search import pipeline
    monkeypatch.setattr(pipeline, 'resolve_model_device', lambda *a: {
        'requested': 'auto', 'effective': effective})
    calls = []
    def generate(name, reference, count, seed, folder, options, threads, workers):
        calls.append((options['device'], workers))
        return [Chem.Mol(reference)]
    monkeypatch.setattr(pipeline, 'generate_batch', generate)
    cfg = SearchConfig.resolve('high', {'threads': 1, 'workers': 2, 'mm_method': 'none',
        'budget': {'formula': 'fixed', 'fixed': 1},
        'generators': {'ditmc': {'command': ['placeholder']}}})
    result = pipeline.prepare_candidates(formaldehyde(), tmp_path/'run', cfg, allocated_cores=2)
    assert calls == [(effective, expected_workers)]
    stage = json.loads(result.status_path.read_text())['stages'][0]
    assert stage['device']['effective'] == effective and stage['workers'] == expected_workers
    assert 'device' not in cfg.generators['ditmc']


def formaldehyde():
    mol = Chem.AddHs(Chem.MolFromSmiles('C=O'))
    assert AllChem.EmbedMolecule(mol, randomSeed=4) == 0
    return mol


def test_generation_parallelism_is_separate_from_native_nproc():
    cfg = SearchConfig.resolve()
    assert cfg.threads == 4
    assert [cfg.parallelism(n) for n in (4, 8, 16)] == [1, 2, 4]
    assert 'cores_per_calculation' not in cfg.relaxation
    with pytest.raises(ValueError, match='up to 4 cores'):
        SearchConfig.resolve(override={'relaxation': {'cores_per_calculation': 4}})


def test_legacy_entry_keeps_original_signature_and_route(tmp_path, monkeypatch):
    import inspect
    from qcforever.laqa_fafoom import laqa_confopt_QCforever as legacy
    from qcforever import laqa_fafoom
    monkeypatch.chdir(tmp_path)
    source = tmp_path/'molecule.sdf'
    with Chem.SDWriter(str(source)) as writer:
        writer.write(formaldehyde())
    assert list(inspect.signature(legacy.LAQA_confopt_main).parameters) == [
        'infilename', 'TotalCharge', 'SpinMulti', 'method', 'nproc', 'mem']
    calls = []
    monkeypatch.setattr(laqa_fafoom.initgeom, 'LAQA_initgeom', lambda *args: calls.append('legacy_generation'))
    monkeypatch.setattr(laqa_fafoom.laqa_optgeom, 'LAQA_optgeom', lambda *args: calls.append('legacy_laqa'))
    legacy.LAQA_confopt_main(str(source), 0, 1, 'xtb', 4, '1GB')
    assert calls == ['legacy_generation', 'legacy_laqa']


@pytest.mark.parametrize('engine', ['gaussian', 'gamess'])
def test_public_runners_call_conformer_search_directly(engine):
    import ast
    import inspect
    from qcforever.gaussian_run.GaussianRunPack import GaussianDFTRun
    from qcforever.gamess_run.GamessRunPack import GamessDFTRun
    cls = GaussianDFTRun if engine == 'gaussian' else GamessDFTRun
    module_tree = ast.parse(inspect.getsource(inspect.getmodule(cls)))
    tree = next(node for node in ast.walk(module_tree)
                if isinstance(node, ast.FunctionDef) and node.name == '_run_'+engine)
    calls = [node.func for node in ast.walk(tree) if isinstance(node, ast.Call)]
    imports = [node.module for node in ast.walk(tree) if isinstance(node, ast.ImportFrom)]
    assert 'qcforever.conformer_search.conformer_search' in imports
    assert any(isinstance(f, ast.Name) and f.id == 'configured_confopt' for f in calls)
    assert not any(isinstance(f, ast.Attribute) and f.attr == 'LAQA_confopt_main' for f in calls)


def test_compressed_bond_is_rejected_without_altering_normal_structure():
    mol = formaldehyde()
    cfg = SearchConfig.resolve().validation
    assert geometry_failure(mol, mol, cfg) is None
    h = next(a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 1)
    conf = mol.GetConformer()
    carbon = conf.GetAtomPosition(0)
    direction = conf.GetAtomPosition(h) - carbon
    conf.SetAtomPosition(h, carbon + direction * (0.40567 / direction.Length()))
    assert geometry_failure(mol, mol, cfg) == 'bond_compressed'
    with pytest.raises(ValueError, match='Minimum bonded'):
        SearchConfig.resolve(override={'validation': {'minimum_bonded_covalent_ratio': 2}})


@pytest.mark.parametrize('log,reason', [
    ('Small interatomic distances encountered\nError termination via Lnk1e in l202.exe', 'invalid_interatomic_distances'),
    ('Convergence failure -- run terminated.', 'scf_not_converged'),
    ('Number of steps exceeded', 'optimization_step_limit'),
    ('', 'missing_native_log'),
    ('Unrecognized output', 'missing_scf_energy'),
])
def test_pm6_failure_classification(log, reason):
    assert failure_diagnostic(log)['reason'] == reason


def test_pm6_parser_error_preserves_native_reason_and_trace(tmp_path, monkeypatch):
    from qcforever.laqa_fafoom.pyg16 import g16Object
    from pathlib import Path
    monkeypatch.setattr('shutil.which', lambda _: '/fake/g16')
    def failed_native(self):
        Path('Gau_molecule.log').write_text('Small interatomic distances encountered\nError termination\n')
        raise UnboundLocalError('energy was not assigned')
    monkeypatch.setattr(g16Object, 'run_g16', failed_native)
    with pytest.raises(PM6CalculationError, match='invalid_interatomic_distances'):
        pm6_relax(formaldehyde(), tmp_path, 0, 1, 4, '1GB', {})
    diagnostic = json.loads((tmp_path/'failure.json').read_text())
    assert 'UnboundLocalError' in diagnostic['adapter_error']
    assert (tmp_path/'Gau_molecule.log').exists() and (tmp_path/'trace.json').exists()


def test_optional_worker_cli_needs_no_core_or_ml_imports():
    import subprocess
    import sys
    from pathlib import Path
    root = Path(__file__).resolve().parents[1]
    code = ('import sys; import qcforever_model_workers.worker; '
            'assert all(x not in sys.modules for x in ("qcforever", "torch", "jax", "rdkit"))')
    # -S disables site-packages, so accidental eager model/core imports fail.
    subprocess.run([sys.executable, '-S', '-c', code], cwd=root, check=True)
    result = subprocess.run([sys.executable, '-S', '-m', 'qcforever_model_workers.worker', '--help'],
                            cwd=root, check=True, capture_output=True, text=True)
    assert '--checkpoint' in result.stdout and '--legacy-workflow' not in result.stdout


@pytest.mark.parametrize('sample', ['pyrroleradical', 'naphthaleneanionradical'])
@pytest.mark.parametrize('method', ['none', 'uff', 'mmff94s'])
def test_mm_preserves_charge_radicals_and_input_through_sdf(sample, method):
    from pathlib import Path
    from qcforever.conformer_search.preoptimization import preoptimize
    from qcforever.conformer_search.validation import graph_key
    path = Path(__file__).resolve().parents[1]/'Samples'/f'{sample}.sdf'
    source = Chem.SDMolSupplier(str(path), removeHs=False)[0]
    before = Chem.MolToMolBlock(source)
    result, statuses = preoptimize([source], method, {})
    reread = Chem.MolFromMolBlock(Chem.MolToMolBlock(result[0]), removeHs=False)
    assert Chem.MolToMolBlock(source) == before
    assert graph_key(reread) == graph_key(source)
    assert [a.GetNumRadicalElectrons() for a in reread.GetAtoms()] == [a.GetNumRadicalElectrons() for a in source.GetAtoms()]
    assert Chem.GetFormalCharge(reread) == Chem.GetFormalCharge(source)


def test_mm_runtime_failure_does_not_discard_other_candidates(tmp_path):
    from qcforever.conformer_search.pipeline import prepare_candidates
    def copies(ref, count, *args):
        return [Chem.Mol(ref) for _ in range(count)]
    calls = []
    def mm(records, method, options):
        assert len(records) == 1
        mol = records[0]
        calls.append(mol.GetProp('candidate_id'))
        if len(calls) == 2:
            mol.GetConformer().SetAtomPosition(0, (99, 99, 99))
            raise RuntimeError('one candidate failed')
        mol.SetBoolProp('mm_test_completed', True)
        return [mol], [0]
    cfg = SearchConfig.resolve('low', {'threads': 1, 'budget': {'formula': 'fixed', 'fixed': 3},
                                     'validation': {'duplicate_rmsd_angstrom': 0}})
    result = prepare_candidates(formaldehyde(), tmp_path/'run', cfg,
                                generators={'etkdgv3': copies}, mm_optimizer=mm)
    status = json.loads(result.status_path.read_text())
    assert status['mm_state'] == 'partially_skipped'
    assert [r['state'] for r in status['mm_candidate_runs']] == ['completed', 'skipped', 'completed']
    assert len(result.candidates) == 3
    assert [m.HasProp('mm_test_completed') for m in result.candidates] == [1, 0, 1]


def test_pm6_takes_coordinates_not_reinterpreted_radical_metadata(tmp_path, monkeypatch):
    from pathlib import Path
    from qcforever.laqa_fafoom.pyg16 import g16Object
    from qcforever.conformer_search.validation import graph_key
    source = Chem.SDMolSupplier(str(Path(__file__).resolve().parents[1]/'Samples/pyrroleradical.sdf'), removeHs=False)[0]
    monkeypatch.setattr('shutil.which', lambda _: '/fake/g16')
    def native(self):
        Path('Gau_molecule.log').write_text('SCF Done: E(UPM6) = -1.0 A.U.\nStationary point found\nNormal termination\n')
        self.energy = -1.0
        mol = Chem.MolFromMolBlock(self.sdf_string, removeHs=False)
        AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94s')
        mol.GetConformer().SetAtomPosition(0, (1, 2, 3))
        self.sdf_string_opt = Chem.MolToMolBlock(mol)
    monkeypatch.setattr(g16Object, 'run_g16', native)
    result, energy = pm6_relax(source, tmp_path, 0, 2, 4, '1GB', {})
    assert energy == -1.0
    assert tuple(result.GetConformer().GetAtomPosition(0)) == (1, 2, 3)
    reread = Chem.MolFromMolBlock(Chem.MolToMolBlock(result), removeHs=False)
    assert graph_key(reread) == graph_key(source)
    assert sum(a.GetNumRadicalElectrons() for a in reread.GetAtoms()) == 1
    assert '\n0 2\n' in (tmp_path/'Gau_molecule.com').read_text()


def test_final_invalid_geometry_has_no_inherited_stereo_pass(tmp_path, monkeypatch):
    from qcforever.conformer_search.pipeline import prepare_candidates
    from qcforever.conformer_search.relaxation import relax_candidates
    from qcforever.conformer_search.validation import STEREO_MATCH_PROPERTIES
    monkeypatch.chdir(tmp_path)
    ref = formaldehyde()
    cfg = SearchConfig.resolve('low', {'threads': 1, 'mm_method': 'none',
                                     'budget': {'formula': 'fixed', 'fixed': 1}})
    prepared = prepare_candidates(ref, tmp_path/'run', cfg)
    assert prepared.candidates[0].GetBoolProp('joint_stereo_match')
    def bad(mol, *args):
        mol.GetConformer().SetAtomPosition(1, mol.GetConformer().GetAtomPosition(0))
        return mol, -1.0
    audit = relax_candidates(prepared, ref, cfg, 0, 1, 'xtb', 1, '1GB', adapter=bad)
    assert audit['selected_structure_warning'] == 'geometry_invalid'
    for path in [tmp_path/'optimized_structures.sdf', tmp_path/'run/electronic/all_converged.sdf',
                 tmp_path/'run/electronic/candidate_00000/optimized.sdf']:
        mol = Chem.SDMolSupplier(str(path), removeHs=False)[0]
        assert not mol.GetBoolProp('geometry_valid')
        assert mol.GetProp('stereo_check_status') == 'not_evaluated_geometry_invalid'
        assert all(not mol.HasProp(key) for key in STEREO_MATCH_PROPERTIES)


def test_failed_gaussian_readback_updates_memory_and_saved_summary(tmp_path, monkeypatch):
    from qcforever.conformer_search.conformer_search import read_selected_structure
    monkeypatch.chdir(tmp_path)
    (tmp_path/'conformer_search').mkdir()
    summary = {'state': 'succeeded', 'energy_hartree': -1}
    with pytest.raises(Exception):
        read_selected_structure(summary)
    assert summary['state'] == 'failed' and summary['search_state'] == 'succeeded'
    assert summary['failure_stage'] == 'selected_structure_readback'
    assert json.loads((tmp_path/'conformer_search/summary.json').read_text()) == summary


def test_failed_search_never_reads_stale_selected_structure(monkeypatch):
    from qcforever.conformer_search.conformer_search import read_selected_structure
    from qcforever.util import read_mol_file
    monkeypatch.setattr(read_mol_file, 'read_sdf', lambda *a: pytest.fail('stale file read'))
    with pytest.raises(RuntimeError, match='generation failed'):
        read_selected_structure({'state': 'failed', 'error': 'generation failed'})


def test_td_embedding_seed_changes_per_request_and_repeats_exactly():
    import numpy as np
    from qcforever_model_workers.torsional import generate
    namespace = {'Chem': Chem, 'AllChem': AllChem}
    exec('def sample(smiles, count, corrected):\n'
         '    m = Chem.AddHs(Chem.MolFromSmiles(smiles))\n'
         '    return embed_func(m, count)\n', namespace)
    sample = namespace['sample']
    first = generate(sample, 'CCCCCC', 2, 43, 1)
    generate(sample, 'CCCCCC', 1, 789, 1)  # another persistent batch
    again = generate(sample, 'CCCCCC', 2, 43, 1)
    other = generate(sample, 'CCCCCC', 2, 44, 1)
    assert np.array_equal(first.GetConformer().GetPositions(), again.GetConformer().GetPositions())
    assert not np.array_equal(first.GetConformer().GetPositions(), other.GetConformer().GetPositions())
