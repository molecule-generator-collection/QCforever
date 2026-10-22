"""Basic conformer-search contracts; no model weights or native xTB/Gaussian needed."""
import json
import sys
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms

from qcforever.conformer_search import generate_conformers as pipeline
from qcforever.conformer_search.settings import SearchConfig
from qcforever.conformer_search.settings import parse_conformer_options, calculation_tokens
from qcforever.conformer_search.generate_conformers import prepare_candidates
from qcforever.conformer_search.optimize_mm import UnsupportedParametersError
from qcforever.conformer_search.optimize_conformers import optimize_candidates
from qcforever.conformer_search.check_structures import (
    evaluate_optimized_candidates, filter_candidates, duplicate_rmsd, xh3_hydrogen_indices,
    rmsd_comparison_molecule, IncrementalCandidateFilter,
)


@pytest.fixture(autouse=True)
def isolated_registry(tmp_path, monkeypatch):
    monkeypatch.setenv('QCFOREVER_MODEL_REGISTRY', str(tmp_path/'models.json'))


def make_butane_conformer():
    mol = Chem.AddHs(Chem.MolFromSmiles('CCCC'))
    assert AllChem.EmbedMolecule(mol, randomSeed=42) == 0
    return mol


def make_test_config(profile='high', n=5, **extra):
    return SearchConfig.resolve(profile, {'budget': {'formula': 'fixed', 'fixed': n},
        'device': 'cpu',  # Synthetic command fixtures intentionally have no GPU probe API.
        'mm_method': 'none', 'workers': 2, 'threads': 1,
        'validation': {'duplicate_rmsd_angstrom': 0}, **extra})


def copy_reference_conformers(reference, count, *args):
    return [Chem.Mol(reference) for _ in range(count)]


def test_default_budget_and_cpu_limits():
    cfg = SearchConfig.resolve()
    assert cfg.profile == 'low' and cfg.mm_method == 'mmff94s'
    assert cfg.threads == 4 and cfg.parallelism(8) == 2 and cfg.parallelism(4) == 1
    assert cfg.parallelism(3) == 1 and cfg.effective_threads(3) == 3
    assert cfg.budget.resolve(Chem.MolFromSmiles('CC'))['maximum_candidates'] == 10
    assert cfg.budget.resolve(Chem.MolFromSmiles('C1CCCCC1'))['maximum_candidates'] == 15


@pytest.mark.parametrize('option,profile,backend_token,generators', [
    ('optconf', 'low', 'optconf', ['etkdgv3']),
    ('optconf laqa', 'low', 'optconf', ['etkdgv3']),
    ('optconf=xtb optconf_medium laqa=50', 'medium', 'optconf=xtb', ['torsional_diffusion', 'etkdgv3']),
    ('optconf=(conf_gen=high, conf_relax=xtb)', 'high', 'optconf=xtb', ['ditmc', 'torsional_diffusion', 'etkdgv3']),
    ('OPTCONF = (CONF_GEN=MEDIUM, CONF_RELAX=PM6, LAQA=50)', 'medium', 'optconf=pm6', ['torsional_diffusion', 'etkdgv3']),
])
def test_search_options_use_expected_defaults(option, profile, backend_token, generators):
    cfg = parse_conformer_options(option)
    assert cfg.profile == profile
    assert [name for name, _ in cfg.stages()] == generators
    assert cfg.relaxation['implementation'] == 'graybox'
    assert cfg.relaxation['convergence_fraction'] == .5
    assert SearchConfig.resolve(profile).relaxation['convergence_fraction'] == .5
    assert calculation_tokens(option + ' energy freq=IR.dat,Raman.dat') == [
        backend_token, 'energy', 'freq=IR.dat,Raman.dat']


def test_explicit_convergence_settings_remain_available():
    for percentage in (20, 80, 100):
        cfg = parse_conformer_options(f'optconf=(laqa={percentage})')
        assert cfg.relaxation['convergence_fraction'] == percentage / 100
    cfg = parse_conformer_options('optconf', {'relaxation': {'convergence_fraction': .3}})
    assert cfg.relaxation['convergence_fraction'] == .3
    cfg = parse_conformer_options('optconf=(laqa=80)', {
        'relaxation': {'algorithm': 'sr', 'convergence_fraction': .3}})
    assert cfg.relaxation['algorithm'] == 'sr' and cfg.relaxation['convergence_fraction'] == .8
    for backend in ('pm6', 'xtb'):
        cfg = parse_conformer_options(f'optconf=(conf_relax={backend}, laqa=off)')
        assert cfg.relaxation['implementation'] == 'continuous'
        if backend == 'pm6':
            assert cfg.relaxation['pm6_optimizer'] == 'rfo'


@pytest.mark.parametrize('option', [
    'optconf_high', 'optconf optconf_high optconf_medium',
    'optconf=(conf_gen=high, conf_gen=low)', 'optconf=(unknown=1)',
    'optconf = (conf_gen=high', 'optconf conf_gen=high',
    'optconf laqa=0', 'optconf laqa=101', 'optconf laqa=nan',
    'optconf laqa=foo', 'optconf laqa laqa=80',
])
def test_invalid_search_options_are_rejected(option):
    with pytest.raises(ValueError):
        parse_conformer_options(option)


@pytest.mark.parametrize('key,value', [
    ('parallel_candidates', 0), ('parallel_candidates', True),
    ('parallel_candidates', 1.5), ('convergence_fraction', True),
])
def test_invalid_relaxation_settings_are_rejected(key, value):
    with pytest.raises(ValueError, match=key):
        SearchConfig.resolve(override={'relaxation': {key: value}})


def test_legacy_is_explicit_and_does_not_accept_modern_overrides():
    assert parse_conformer_options('energy opt') is None
    cfg = parse_conformer_options('optconf=(conf_gen=legacy, conf_relax=xtb)')
    assert cfg.profile == 'legacy'
    with pytest.raises(ValueError, match='old generation and relaxation'):
        parse_conformer_options('optconf=(conf_gen=legacy)', {'seed': 42})


@pytest.mark.parametrize('state', ['success', 'missing_result', 'failed_with_stale_result'])
def test_gamess_handoff_keeps_elements_with_their_coordinates(tmp_path, monkeypatch, state):
    from qcforever.gamess_run.GamessRunPack import GamessDFTRun
    from qcforever.conformer_search import conformer_search as workflow

    monkeypatch.chdir(tmp_path)
    original = Chem.AddHs(Chem.MolFromSmiles('OCC'))
    assert AllChem.EmbedMolecule(original, randomSeed=22) == 0
    selected = Chem.RenumberAtoms(original, [2, 1, 0, *range(3, 9)])
    for index, position in enumerate(selected.GetConformer().GetPositions()):
        selected.GetConformer().SetAtomPosition(index, position + [4, 2, -3])
    with Chem.SDWriter('molecule.sdf') as writer:
        writer.write(original)

    def search(*args, **kwargs):
        Path('conformer_search').mkdir()
        if state != 'missing_result':
            with Chem.SDWriter('optimized_structures.sdf') as writer:
                writer.write(selected)
        return {'state': 'failed' if state == 'failed_with_stale_result' else 'succeeded'}

    monkeypatch.setattr(workflow, 'run_conformer_search', search)
    monkeypatch.setattr('qcforever.util.check_resource.respec_cores', lambda *a: 1)
    monkeypatch.setattr('qcforever.util.check_resource.respec_memory', lambda **kwargs: '1GB')
    monkeypatch.setattr('qcforever.gamess_run.Exe_Gamess.exe_Gamess', lambda *a: '')
    monkeypatch.setattr(GamessDFTRun, 'Extract_values', lambda *a: {})
    result = GamessDFTRun('B3LYP', 'STO-3G', 1,
        'energy optconf=(conf_gen=medium, conf_relax=xtb)', 'molecule.sdf').run_gamess()
    lines = (tmp_path/'molecule/molecule.inp').read_text().splitlines()
    start = lines.index('C1') + 1
    records = [line.split() for line in lines[start:start + original.GetNumAtoms()]]
    expected = selected if state == 'success' else original
    assert [row[0] for row in records] == [atom.GetSymbol() for atom in expected.GetAtoms()]
    np.testing.assert_allclose([[float(x) for x in row[2:5]] for row in records],
                               expected.GetConformer().GetPositions(), atol=6e-5)
    assert result['optconf'] is (state == 'success')
    if state == 'missing_result':
        summary = json.loads((tmp_path/'molecule/conformer_search/summary.json').read_text())
        assert summary['failure_stage'] == 'selected_structure_readback'
        assert summary['state'] == result['conformer_search']['state'] == 'failed'


def test_pm6_cleanup_preserves_other_jobs_and_analysis_files(tmp_path):
    from qcforever.conformer_search.run_pm6 import remove_pm6_checkpoints

    electronic = tmp_path/'electronic'
    call = electronic/'candidate_00000/call_0000'
    call.mkdir(parents=True)
    for name in ('state.chk', 'state.fchk', 'input.log', 'trajectory.json', 'final.sdf'):
        (call/name).write_text('retained analysis or removable restart')
    other_job = tmp_path/'DFT.chk'
    other_job.write_text('unrelated checkpoint')
    (electronic/'linked.chk').symlink_to(other_job)
    remove_pm6_checkpoints(electronic)
    assert not (call/'state.chk').exists() and not (call/'state.fchk').exists()
    assert all((call/name).is_file() for name in ('input.log', 'trajectory.json', 'final.sdf'))
    assert other_job.read_text() == 'unrelated checkpoint'
    assert (electronic/'linked.chk').is_symlink()
    report = json.loads((electronic/'checkpoint_cleanup.json').read_text())
    assert report['removed_files'] == 2 and report['removed_bytes'] > 0 and not report['errors']


def test_yaml_overrides_respect_allocated_cores(tmp_path):
    path = tmp_path/'conformer.yaml'
    path.write_text('mm_method: uff\nworkers: 8\nthreads: 2\n')
    cfg = parse_conformer_options('optconf=xtb optconf_medium', str(path))
    assert cfg.profile == 'medium' and cfg.budget.base == 10 and cfg.mm_method == 'uff'
    assert cfg.parallelism(8) == 4
    assert cfg.parallelism(1) == 1
    assert cfg.effective_threads(1) == 1


def test_generation_stops_at_candidate_limit(tmp_path):
    requested = []
    def adapter(reference, count, *args):
        requested.append(count)
        return copy_reference_conformers(reference, min(2, count))
    out = prepare_candidates(make_butane_conformer(), tmp_path/'run', make_test_config(), allocated_cores=2,
                             generators={'ditmc': adapter})
    assert requested == [5, 2, 1]
    assert len(out.candidates) == 5
    status = json.loads(out.status_path.read_text())
    assert len(status['stages']) == 1
    assert len({m.GetProp('candidate_id') for m in out.candidates}) == 5


def test_fallback_keeps_candidates_and_respects_attempt_limit(tmp_path):
    calls = {'ditmc': [], 'torsional_diffusion': []}
    def generate_one_ditmc_candidate(reference, count, *args):
        calls['ditmc'].append(count)
        return copy_reference_conformers(reference, 1) if len(calls['ditmc']) == 1 else []
    def generate_torsional_candidates(reference, count, *args):
        calls['torsional_diffusion'].append(count)
        return copy_reference_conformers(reference, count)
    out = prepare_candidates(make_butane_conformer(), tmp_path/'run', make_test_config(), allocated_cores=2,
        generators={'ditmc': generate_one_ditmc_candidate, 'torsional_diffusion': generate_torsional_candidates})
    assert sum(calls['ditmc']) == 10
    assert calls['torsional_diffusion'] == [4]
    assert out.candidates[0].GetProp('generator') == 'ditmc'
    assert len(out.candidates) == 5


def test_generation_continues_when_models_or_mm_are_unavailable(tmp_path):
    def skip(*args):
        raise UnsupportedParametersError('unsupported test molecule')
    out = prepare_candidates(make_butane_conformer(), tmp_path/'run', make_test_config(mm_method='mmff94s'),
        allocated_cores=2, generators={'etkdgv3': copy_reference_conformers}, mm_optimizer=skip)
    state = json.loads(out.status_path.read_text())
    assert [x['state'] for x in state['stages'][:2]] == ['unavailable_dependency']*2
    assert state['mm_state'] == 'skipped' and len(out.candidates) == 5


def test_generation_accepts_nonempty_candidate_sets(tmp_path):
    def one(reference, count, seed, directory, *args):
        return copy_reference_conformers(reference, 1) if 'batch_000' in str(directory) else []
    cfg = make_test_config('low')
    out = prepare_candidates(make_butane_conformer(), tmp_path/'one', cfg, generators={'etkdgv3': one})
    assert len(out.candidates) == 1 and out.initial_sdf is not None
    failed = prepare_candidates(make_butane_conformer(), tmp_path/'zero', cfg, generators={'etkdgv3': lambda *a: []})
    assert failed.state == 'generation_failed' and failed.initial_sdf is None


@pytest.mark.parametrize('method', ['mmff94s', 'uff', 'none'])
def test_etkdg_generation_with_mm_settings(tmp_path, method):
    out = prepare_candidates('CCCC', tmp_path/'run', make_test_config('low', n=3, mm_method=method, threads=4),
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
    ('CCO', 3), ('C[NH3+]', 6), ('C[NH2+]C', 6), ('C', 0),
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
    cfg = make_test_config()
    accepted, generated = filter_candidates([ref, wrong], ref, 2, cfg.validation)
    assert len(accepted) == 1 and generated['rejected'] == {'stereo_mismatch': 1}
    assert generated['candidates'][1]['tetrahedral_stereo_match'] is False
    audit = evaluate_optimized_candidates([ref, wrong], [-10.0, -11.0], ref, cfg.validation)
    assert audit['best_any_index'] == 1 and audit['primary_best_index'] == 0


def test_ez_mismatch_rejected_at_generation():
    ref = Chem.AddHs(Chem.MolFromSmiles('F/C=C/F'))
    assert AllChem.EmbedMolecule(ref, randomSeed=42) == 0
    wrong = Chem.Mol(ref)
    rdMolTransforms.SetDihedralDeg(wrong.GetConformer(), 0, 1, 2, 3, 0)
    kept, audit = filter_candidates([wrong], ref, 2, SearchConfig.resolve().validation)
    assert not kept and audit['rejected'] == {'stereo_mismatch': 1}
    assert not audit['candidates'][0]['ez_stereo_match']


def test_learned_candidates_never_receive_mm_even_in_mixed_pool(tmp_path):
    calls = []
    def optimizer(records, method, options):
        calls.extend(m.GetProp('generator') for m in records)
        return [Chem.Mol(m) for m in records], [0]*len(records)
    cfg = make_test_config('high', n=3, mm_method='mmff94s')
    def one_learned(ref, count, seed, folder, *a):
        return [Chem.Mol(ref)] if folder.name == 'batch_000' else []
    out = prepare_candidates(make_butane_conformer(), tmp_path/'mixed', cfg,
        generators={'ditmc': one_learned,
                    'torsional_diffusion': lambda *a: [], 'etkdgv3': copy_reference_conformers},
        mm_optimizer=optimizer)
    # Restrict the learned stage explicitly to its first batch so ETKDG fills.
    # A separate learned-only pool must never call the optimizer at all.
    only = prepare_candidates(make_butane_conformer(), tmp_path/'learned', cfg,
        generators={'ditmc': copy_reference_conformers}, mm_optimizer=lambda *a: pytest.fail('MM called on DiTMC'))
    assert json.loads(only.status_path.read_text())['mm_state'] == 'not_applicable_to_learned_generators'
    assert calls == ['etkdgv3', 'etkdgv3']
    assert [m.GetProp('generator') for m in out.candidates] == ['ditmc', 'etkdgv3', 'etkdgv3']


def test_final_stereo_mismatch_returns_structure_with_warning(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    ref = Chem.AddHs(Chem.MolFromSmiles('C[C@H](O)F'))
    assert AllChem.EmbedMolecule(ref, randomSeed=42) == 0
    cfg = make_test_config('low', n=1, relaxation={'implementation': 'continuous'})
    prepared = prepare_candidates(ref, tmp_path/'run', cfg, generators={'etkdgv3': copy_reference_conformers})
    def backend(mol, *a):
        wrong = Chem.Mol(mol)
        conf = wrong.GetConformer()
        for i, (x, y, z) in enumerate(conf.GetPositions()):
            conf.SetAtomPosition(i, (-x, y, z))
        return wrong, -1
    audit = optimize_candidates(prepared, ref, cfg, 0, 1, 'xtb', 1, '1GB', adapter=backend)
    assert audit['primary_best_index'] is None and audit['selected_index'] == 0
    assert audit['selected_structure_warning'] == 'stereo_mismatch'
    mol = Chem.SDMolSupplier(str(tmp_path/'optimized_structures.sdf'), removeHs=False)[0]
    assert mol.GetProp('structure_check_warning') == 'stereo_mismatch'
    assert not mol.GetBoolProp('tetrahedral_stereo_match')


def test_xtb_parallel_relaxation_uses_single_core_workers(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    binary = tmp_path/'fake_xtb_parallel'
    binary.write_text(f'#!{sys.executable}\n' +
        'import json,os,shutil,time\nfrom pathlib import Path\n'
        'start=time.time()\ntime.sleep(0.5)\n'
        'shutil.copyfile("input.xyz","xtbopt.xyz")\n'
        'Path("interval.json").write_text(json.dumps({"start":start,"end":time.time(),"pid":os.getpid(),"omp":os.environ["OMP_NUM_THREADS"],"blas":os.environ["OPENBLAS_NUM_THREADS"]}))\n'
        'print("CYCLE 1\\nGEOMETRY OPTIMIZATION CONVERGED\\nTOTAL ENERGY -1.0 Eh")\n')
    binary.chmod(0o755)
    cfg = make_test_config('low', n=2, relaxation={'implementation': 'continuous', 'xtb_executable': str(binary)})
    ref = make_butane_conformer()
    prepared = prepare_candidates(ref, tmp_path/'run', cfg, generators={'etkdgv3': copy_reference_conformers})
    audit = optimize_candidates(prepared, ref, cfg, 0, 1, 'xtb', 2, '1GB')
    assert audit['parallel_workers'] == 2 and audit['cores_per_calculation'] == 1
    assert audit['scheduling'] == 'parallel_candidates_1_core'
    folder = tmp_path/'run/electronic'
    intervals = [json.loads((folder/f'candidate_{i:05d}/interval.json').read_text()) for i in range(2)]
    assert intervals[0]['pid'] != intervals[1]['pid']
    assert max(r['start'] for r in intervals) < min(r['end'] for r in intervals)
    assert all(r['omp'] == r['blas'] == '1' for r in intervals)
    for i in range(2):
        argv = json.loads((folder/f'candidate_{i:05d}/command.json').read_text())['argv']
        assert argv[argv.index('--parallel')+1] == '1'


def test_no_mm_snapshot_keeps_results_without_calling_optimizer(tmp_path):
    ref, cfg = make_butane_conformer(), make_test_config('low', n=2, threads=4)
    prepared = prepare_candidates(ref, tmp_path/'run', cfg, allocated_cores=3,
        generators={'etkdgv3': copy_reference_conformers},
        mm_optimizer=lambda *a: pytest.fail('MM called despite mm_method=none'))
    assert len(prepared.candidates) == 2
    assert json.loads(prepared.status_path.read_text())['post_mm_validation_mode'] == 'reused_generation_no_mm'
    cache = IncrementalCandidateFilter(ref, 2, cfg.validation)
    accepted, _ = cache.extend(copy_reference_conformers(ref, 2))
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
    result = prepare_candidates(make_butane_conformer(), tmp_path/'run', make_test_config('low', n=1, mm_method='uff'),
        generators={'etkdgv3': copy_reference_conformers}, mm_optimizer=collision)
    assert result.state == 'no_valid_candidates_after_mm'
    assert json.loads(result.status_path.read_text())['post_mm_audit']['rejected'] == {'collision': 1}


@pytest.mark.parametrize('device,workers', [('cpu', 2), ('gpu', 1)])
def test_auto_device_controls_model_worker_count(tmp_path, monkeypatch, device, workers):
    monkeypatch.setattr(pipeline, 'resolve_model_device', lambda *a: {'requested': 'auto', 'effective': device})
    def generate(name, ref, count, seed, folder, options, threads, actual_workers):
        assert options['device'] == device and actual_workers == workers and threads == 4
        return copy_reference_conformers(ref, count)
    monkeypatch.setattr(pipeline, 'generate_batch', generate)
    cfg = make_test_config('high', n=2, device='auto', threads=4)
    result = prepare_candidates(make_butane_conformer(), tmp_path/'run', cfg, allocated_cores=8)
    assert len(result.candidates) == 2
