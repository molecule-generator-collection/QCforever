"""Basic conformer-search contracts; no model weights or native xTB/Gaussian needed."""
import json
import sys

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms

from qcforever.conformer_search import pipeline
from qcforever.conformer_search.config import SearchConfig
from qcforever.conformer_search.options import resolve_options, calculation_tokens
from qcforever.conformer_search.pipeline import prepare_candidates
from qcforever.conformer_search.preoptimization import UnsupportedParametersError
from qcforever.conformer_search.relaxation import relax_candidates
from qcforever.conformer_search.validation import (
    choose_best, filter_candidates, duplicate_rmsd, xh3_hydrogen_indices,
    rmsd_comparison_molecule, IncrementalCandidateFilter,
)


@pytest.fixture(autouse=True)
def isolated_registry(tmp_path, monkeypatch):
    monkeypatch.setenv('QCFOREVER_MODEL_REGISTRY', str(tmp_path/'models.json'))


def molecule():
    mol = Chem.AddHs(Chem.MolFromSmiles('CCCC'))
    assert AllChem.EmbedMolecule(mol, randomSeed=42) == 0
    return mol


def settings(profile='high', n=5, **extra):
    return SearchConfig.resolve(profile, {'budget': {'formula': 'fixed', 'fixed': n},
        'device': 'cpu',  # Synthetic command fixtures intentionally have no GPU probe API.
        'mm_method': 'none', 'workers': 2, 'threads': 1,
        'validation': {'duplicate_rmsd_angstrom': 0}, **extra})


def copies(reference, count, *args):
    return [Chem.Mol(reference) for _ in range(count)]


def test_defaults_and_options():
    cfg = resolve_options('optconf=xtb opt energy uv')
    assert cfg.profile == 'low' and cfg.mm_method == 'mmff94s'
    assert cfg.threads == 4 and cfg.parallelism(8) == 2 and cfg.parallelism(4) == 1
    assert cfg.parallelism(3) == 1 and cfg.effective_threads(3) == 3
    assert cfg.budget.resolve(Chem.MolFromSmiles('CC'))['maximum_candidates'] == 10
    assert cfg.budget.resolve(Chem.MolFromSmiles('C1CCCCC1'))['maximum_candidates'] == 15
    assert resolve_options('optconf optconf_high uv').profile == 'high'
    assert calculation_tokens('optconf=xtb optconf_medium opt uv') == ['optconf=xtb', 'opt', 'uv']
    assert resolve_options('opt energy uv') is None
    with pytest.raises(ValueError):
        resolve_options('optconf=xtb optconf_high optconf_medium')
    with pytest.raises(ValueError):
        resolve_options('optconf_high uv')


@pytest.mark.parametrize('level,names', [
    ('low', ['etkdgv3']),
    ('medium', ['torsional_diffusion', 'etkdgv3']),
    ('high', ['ditmc', 'torsional_diffusion', 'etkdgv3']),
])
def test_canonical_level_names(level, names):
    options = f'optconf=xtb optconf_{level} opt energy uv'
    cfg = resolve_options(options)
    assert cfg.profile == level
    assert [name for name, _ in cfg.stages()] == names
    assert calculation_tokens(options) == ['optconf=xtb', 'opt', 'energy', 'uv']
    assert SearchConfig.resolve(override={'profile': level}).profile == level


def test_partial_yaml_and_cpu_cap(tmp_path):
    path = tmp_path/'conformer.yaml'
    path.write_text('mm_method: uff\nworkers: 8\nthreads: 2\n')
    cfg = resolve_options('optconf=xtb optconf_medium', str(path))
    assert cfg.profile == 'medium' and cfg.budget.base == 10 and cfg.mm_method == 'uff'
    assert cfg.parallelism(8) == 4
    assert cfg.parallelism(1) == 1
    assert cfg.effective_threads(1) == 1


def test_incremental_stops_at_n(tmp_path):
    requested = []
    def adapter(reference, count, *args):
        requested.append(count)
        return copies(reference, min(2, count))
    out = prepare_candidates(molecule(), tmp_path/'run', settings(), allocated_cores=2,
                             generators={'ditmc': adapter})
    assert requested == [5, 2, 1]
    assert len(out.candidates) == 5
    status = json.loads(out.status_path.read_text())
    assert len(status['stages']) == 1
    assert len({m.GetProp('candidate_id') for m in out.candidates}) == 5


def test_merge_and_fixed_attempt_caps(tmp_path):
    calls = {'ditmc': [], 'torsional_diffusion': []}
    def d(reference, count, *args):
        calls['ditmc'].append(count)
        return copies(reference, 1) if len(calls['ditmc']) == 1 else []
    def t(reference, count, *args):
        calls['torsional_diffusion'].append(count)
        return copies(reference, count)
    out = prepare_candidates(molecule(), tmp_path/'run', settings(), allocated_cores=2,
        generators={'ditmc': d, 'torsional_diffusion': t})
    assert sum(calls['ditmc']) == 10
    assert calls['torsional_diffusion'] == [4]
    assert out.candidates[0].GetProp('generator') == 'ditmc'
    assert len(out.candidates) == 5


def test_unavailable_models_and_mm_skip(tmp_path):
    def skip(*args):
        raise UnsupportedParametersError('unsupported test molecule')
    out = prepare_candidates(molecule(), tmp_path/'run', settings(mm_method='mmff94s'),
        allocated_cores=2, generators={'etkdgv3': copies}, mm_optimizer=skip)
    state = json.loads(out.status_path.read_text())
    assert [x['state'] for x in state['stages'][:2]] == ['unavailable_dependency']*2
    assert state['mm_state'] == 'skipped' and len(out.candidates) == 5


def test_underfilled_and_zero(tmp_path):
    def one(reference, count, seed, directory, *args):
        return copies(reference, 1) if 'batch_000' in str(directory) else []
    cfg = settings('low')
    out = prepare_candidates(molecule(), tmp_path/'one', cfg, generators={'etkdgv3': one})
    assert len(out.candidates) == 1 and out.initial_sdf is not None
    failed = prepare_candidates(molecule(), tmp_path/'zero', cfg, generators={'etkdgv3': lambda *a: []})
    assert failed.state == 'generation_failed' and failed.initial_sdf is None


@pytest.mark.parametrize('method', ['mmff94s', 'uff', 'none'])
def test_native_etkdg_mm(tmp_path, method):
    out = prepare_candidates('CCCC', tmp_path/'run', settings('low', n=3, mm_method=method, threads=4),
                             allocated_cores=3)
    assert out.initial_sdf is not None
    assert 1 <= len(out.candidates) <= 3
    state = json.loads(out.status_path.read_text())
    assert state['threads_per_worker'] == 3 and state['workers'] == 1
    mode = 'reused_generation_no_mm' if method == 'none' else 'full_recheck'
    assert state['post_mm_validation_mode'] == mode


def test_all_atom_duplicates_preserve_oh_rotamers_and_coordinates():
    mol = Chem.AddHs(Chem.MolFromSmiles('CCO'))
    assert AllChem.EmbedMolecule(mol, randomSeed=42) == 0
    atoms = mol.GetSubstructMatch(Chem.MolFromSmarts('[C]-[C]-[O]-[H]'))
    anti, gauche = Chem.Mol(mol), Chem.Mol(mol)
    rdMolTransforms.SetDihedralDeg(anti.GetConformer(), *atoms, 180)
    rdMolTransforms.SetDihedralDeg(gauche.GetConformer(), *atoms, 60)
    cfg = SearchConfig.resolve().validation
    before = gauche.GetConformer().GetPositions().copy()
    assert duplicate_rmsd(anti, gauche, cfg) > cfg.duplicate_rmsd_angstrom
    assert (before == gauche.GetConformer().GetPositions()).all()
    accepted, audit = filter_candidates([anti, gauche, Chem.Mol(anti)], mol, 10, cfg)
    assert len(accepted) == 2
    assert audit['rejected'] == {'duplicate': 1}
    assert audit['duplicate_rmsd_atoms'] == 'all_explicit_atoms_except_XH3_hydrogens'


@pytest.mark.parametrize('smiles,expected', [
    ('CCO', 3), ('C[NH3+]', 6), ('C[NH2+]C', 6),
    ('C', 0), ('O', 0), ('N', 3), ('[NH4+]', 0),
])
def test_xh3_hydrogen_exclusion_rule(smiles, expected):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert len(xh3_hydrogen_indices(mol)) == expected
    reduced = rmsd_comparison_molecule(mol)
    assert reduced.GetNumAtoms() == mol.GetNumAtoms() - expected
    assert reduced.GetNumHeavyAtoms() == mol.GetNumHeavyAtoms()
    assert Chem.MolToSmiles(Chem.RemoveHs(reduced)) == Chem.MolToSmiles(Chem.RemoveHs(mol))


def test_wrong_stereo_is_rejected_at_generation_and_not_primary_best():
    ref = Chem.AddHs(Chem.MolFromSmiles('C[C@H](O)F'))
    AllChem.EmbedMolecule(ref, randomSeed=21)
    wrong = Chem.Mol(ref)
    conf = wrong.GetConformer()
    for i in range(wrong.GetNumAtoms()):
        point = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, (-point.x, point.y, point.z))
    cfg = settings()
    accepted, generated = filter_candidates([ref, wrong], ref, 2, cfg.validation)
    assert len(accepted) == 1 and generated['rejected'] == {'stereo_mismatch': 1}
    assert generated['candidates'][1]['tetrahedral_stereo_match'] is False
    audit = choose_best([ref, wrong], [-10.0, -11.0], ref, cfg.validation)
    assert audit['best_any_index'] == 1 and audit['primary_best_index'] == 0


def test_ez_mismatch_rejected_at_generation():
    ref = Chem.AddHs(Chem.MolFromSmiles('F/C=C/F'))
    assert AllChem.EmbedMolecule(ref, randomSeed=42) == 0
    wrong = Chem.Mol(ref)
    rdMolTransforms.SetDihedralDeg(wrong.GetConformer(), 0, 1, 2, 3, 0)
    kept, audit = filter_candidates([wrong], ref, 2, SearchConfig.resolve().validation)
    assert not kept and audit['rejected'] == {'stereo_mismatch': 1}
    assert not audit['candidates'][0]['ez_stereo_match']


def test_all_candidates_attempted_even_if_one_fails(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    cfg = settings('low', n=5)
    ref = molecule()
    prepared = prepare_candidates(ref, tmp_path/'run', cfg, generators={'etkdgv3': copies})
    calls = []
    def backend(mol, directory, *args):
        i = int(directory.name.split('_')[-1])
        calls.append(i)
        if i == 1:
            raise UnboundLocalError('simulate old parser on incomplete native output')
        return Chem.Mol(mol), -10-i
    audit = relax_candidates(prepared, ref, cfg, 0, 1, 'xtb', 1, '1GB', adapter=backend)
    assert calls == [0, 1, 2, 3, 4]
    assert audit['attempted_candidates'] == 5
    assert audit['converged_candidates'] == 4 and audit['failed_candidates'] == 1
    assert audit['primary_valid_candidates'] == 4 and audit['primary_best_index'] == 4


def test_learned_candidates_never_receive_mm_even_in_mixed_pool(tmp_path):
    calls = []
    def optimizer(records, method, options):
        calls.extend(m.GetProp('generator') for m in records)
        return [Chem.Mol(m) for m in records], [0]*len(records)
    cfg = settings('high', n=3, mm_method='mmff94s')
    def one_learned(ref, count, seed, folder, *a):
        return [Chem.Mol(ref)] if folder.name == 'batch_000' else []
    out = prepare_candidates(molecule(), tmp_path/'mixed', cfg,
        generators={'ditmc': one_learned,
                    'torsional_diffusion': lambda *a: [], 'etkdgv3': copies},
        mm_optimizer=optimizer)
    # Restrict the learned stage explicitly to its first batch so ETKDG fills.
    # A separate learned-only pool must never call the optimizer at all.
    only = prepare_candidates(molecule(), tmp_path/'learned', cfg,
        generators={'ditmc': copies}, mm_optimizer=lambda *a: pytest.fail('MM called on DiTMC'))
    assert json.loads(only.status_path.read_text())['mm_state'] == 'not_applicable_to_learned_generators'
    assert calls == ['etkdgv3', 'etkdgv3']
    assert [m.GetProp('generator') for m in out.candidates] == ['ditmc', 'etkdgv3', 'etkdgv3']


def test_final_stereo_mismatch_returns_structure_with_warning(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    ref = Chem.AddHs(Chem.MolFromSmiles('C[C@H](O)F'))
    assert AllChem.EmbedMolecule(ref, randomSeed=42) == 0
    cfg = settings('low', n=1)
    prepared = prepare_candidates(ref, tmp_path/'run', cfg, generators={'etkdgv3': copies})
    def backend(mol, *a):
        wrong = Chem.Mol(mol)
        conf = wrong.GetConformer()
        for i, (x, y, z) in enumerate(conf.GetPositions()):
            conf.SetAtomPosition(i, (-x, y, z))
        return wrong, -1
    audit = relax_candidates(prepared, ref, cfg, 0, 1, 'xtb', 1, '1GB', adapter=backend)
    assert audit['primary_best_index'] is None and audit['selected_index'] == 0
    assert audit['selected_structure_warning'] == 'stereo_mismatch'
    mol = Chem.SDMolSupplier(str(tmp_path/'optimized_structures.sdf'), removeHs=False)[0]
    assert mol.GetProp('structure_check_warning') == 'stereo_mismatch'
    assert not mol.GetBoolProp('tetrahedral_stereo_match')


@pytest.mark.parametrize('nproc', [1, 2])
def test_native_relaxation_uses_single_core_workers(tmp_path, monkeypatch, nproc):
    monkeypatch.chdir(tmp_path)
    binary = tmp_path/'fake_xtb_parallel'
    binary.write_text(f'#!{sys.executable}\n' +
        'import json,os,shutil,time\nfrom pathlib import Path\n'
        'start=time.time()\ntime.sleep(0.5)\n'
        'shutil.copyfile("input.xyz","xtbopt.xyz")\n'
        'Path("interval.json").write_text(json.dumps({"start":start,"end":time.time(),"pid":os.getpid(),"omp":os.environ["OMP_NUM_THREADS"],"blas":os.environ["OPENBLAS_NUM_THREADS"]}))\n'
        'print("CYCLE 1\\nGEOMETRY OPTIMIZATION CONVERGED\\nTOTAL ENERGY -1.0 Eh")\n')
    binary.chmod(0o755)
    cfg = settings('low', n=2, relaxation={'xtb_executable': str(binary)})
    ref = molecule()
    prepared = prepare_candidates(ref, tmp_path/'run', cfg, generators={'etkdgv3': copies})
    audit = relax_candidates(prepared, ref, cfg, 0, 1, 'xtb', nproc, '1GB')
    per_call = 1
    workers = min(2, nproc)
    assert audit['parallel_workers'] == workers and audit['cores_per_calculation'] == per_call
    assert audit['scheduling'] == 'parallel_candidates_1_core'
    folder = tmp_path/'run/electronic'
    intervals = [json.loads((folder/f'candidate_{i:05d}/interval.json').read_text()) for i in range(2)]
    assert intervals[0]['pid'] != intervals[1]['pid']
    if workers == 1:
        assert intervals[0]['end'] <= intervals[1]['start']
    else:
        assert max(r['start'] for r in intervals) < min(r['end'] for r in intervals)
    assert all(r['omp'] == r['blas'] == str(per_call) for r in intervals)
    for i in range(2):
        argv = json.loads((folder/f'candidate_{i:05d}/command.json').read_text())['argv']
        assert argv[argv.index('--parallel')+1] == str(per_call)


def test_no_mm_snapshot_keeps_results_without_calling_optimizer(tmp_path):
    ref, cfg = molecule(), settings('low', n=2, threads=4)
    prepared = prepare_candidates(ref, tmp_path/'run', cfg, allocated_cores=3,
        generators={'etkdgv3': copies},
        mm_optimizer=lambda *a: pytest.fail('MM called despite mm_method=none'))
    assert len(prepared.candidates) == 2
    assert json.loads(prepared.status_path.read_text())['post_mm_validation_mode'] == 'reused_generation_no_mm'
    cache = IncrementalCandidateFilter(ref, 2, cfg.validation)
    accepted, _ = cache.extend(copies(ref, 2))
    expected, audit = filter_candidates(accepted, ref, 2, cfg.validation)
    actual, snapshot = cache.snapshot()
    assert snapshot == audit
    for left, right in zip(actual, expected):
        assert left.GetPropsAsDict() == right.GetPropsAsDict()
        assert (left.GetConformer().GetPositions() == right.GetConformer().GetPositions()).all()


def test_mm_outputs_are_rechecked(tmp_path):
    def collision(records, *args):
        mol = records[0]
        mol.GetConformer().SetAtomPosition(1, mol.GetConformer().GetAtomPosition(0))
        return records, [0]
    result = prepare_candidates(molecule(), tmp_path/'run', settings('low', n=1, mm_method='uff'),
        generators={'etkdgv3': copies}, mm_optimizer=collision)
    assert result.state == 'no_valid_candidates_after_mm'
    assert json.loads(result.status_path.read_text())['post_mm_audit']['rejected'] == {'collision': 1}


@pytest.mark.parametrize('device,workers', [('cpu', 2), ('gpu', 1)])
def test_auto_device_controls_model_worker_count(tmp_path, monkeypatch, device, workers):
    monkeypatch.setattr(pipeline, 'resolve_model_device', lambda *a: {'requested': 'auto', 'effective': device})
    def generate(name, ref, count, seed, folder, options, threads, actual_workers):
        assert options['device'] == device and actual_workers == workers and threads == 4
        return copies(ref, count)
    monkeypatch.setattr(pipeline, 'generate_batch', generate)
    cfg = settings('high', n=2, device='auto', threads=4)
    result = prepare_candidates(molecule(), tmp_path/'run', cfg, allocated_cores=8)
    assert len(result.candidates) == 2
