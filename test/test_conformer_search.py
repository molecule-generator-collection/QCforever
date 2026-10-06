"""CPU-only contracts: no learned model or external QC executable required."""
import json
import sys
from pathlib import Path

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms

from qcforever.conformer_search.config import SearchConfig
from qcforever.conformer_search.options import resolve_options, calculation_tokens
from qcforever.conformer_search.pipeline import prepare_candidates
from qcforever.conformer_search.generators import generate_batch, command_generator, MissingGeneratorDependency
from qcforever.conformer_search.preoptimization import UnsupportedParametersError
from qcforever.conformer_search.relaxation import relax_candidates, xtb_relax
from qcforever.conformer_search.validation import (
    choose_best, filter_candidates, duplicate_rmsd,
    xh3_hydrogen_indices, rmsd_comparison_molecule,
)
from qcforever.util import job_cleanup, job_timeout


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


def test_all_atom_duplicates_handle_reordered_equivalent_hydrogens():
    mol = Chem.AddHs(Chem.MolFromSmiles('CCO'))
    assert AllChem.EmbedMolecule(mol, randomSeed=42) == 0
    order = list(range(mol.GetNumAtoms()))
    # CH2 hydrogens are retained; their equivalent labels must still be matched.
    hydrogens = [a.GetIdx() for a in mol.GetAtomWithIdx(1).GetNeighbors() if a.GetAtomicNum() == 1]
    order[hydrogens[0]], order[hydrogens[1]] = order[hydrogens[1]], order[hydrogens[0]]
    reordered = Chem.RenumberAtoms(mol, order)
    cfg = SearchConfig.resolve().validation
    assert duplicate_rmsd(mol, reordered, cfg) < 1e-6
    accepted, _ = filter_candidates([mol, reordered], mol, 10, cfg)
    assert len(accepted) == 1


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


@pytest.mark.parametrize('smiles', ['CCO', 'CC[NH3+]'])
def test_xh3_rotation_is_ignored_without_changing_full_structure(smiles):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=42) == 0
    moved = Chem.Mol(mol)
    # A rigid turn of just the terminal XH3 hydrogen tripod.
    pattern = Chem.MolFromSmarts('[!#1]-[!#1]-[!#1]-[#1]')
    excluded = xh3_hydrogen_indices(moved)
    matches = moved.GetSubstructMatches(pattern)
    match = next(m for m in matches if m[-1] in excluded)
    angle = rdMolTransforms.GetDihedralDeg(moved.GetConformer(), *match)
    rdMolTransforms.SetDihedralDeg(moved.GetConformer(), *match, angle + 40)
    before = moved.GetConformer().GetPositions().copy()
    cfg = SearchConfig.resolve().validation
    assert duplicate_rmsd(mol, moved, cfg) < 1e-6
    assert (before == moved.GetConformer().GetPositions()).all()
    assert moved.GetNumAtoms() == mol.GetNumAtoms()


@pytest.mark.parametrize('value', [0, -1, 1.5, True])
def test_duplicate_mapping_cap_validation(value):
    with pytest.raises(ValueError, match='duplicate_max_matches'):
        SearchConfig.resolve(override={'validation': {'duplicate_max_matches': value}})


def test_defaults_and_options():
    cfg = resolve_options('optconf=xtb opt energy uv')
    assert cfg.profile == 'low' and cfg.mm_method == 'mmff94s'
    assert cfg.threads == 4 and cfg.parallelism(8) == 2 and cfg.parallelism(4) == 1
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


@pytest.mark.parametrize('old_name', ['light', 'middle'])
def test_old_level_names_are_rejected_not_aliased(tmp_path, old_name):
    with pytest.raises(ValueError, match='Use optconf_low, optconf_medium, or optconf_high'):
        resolve_options(f'optconf=xtb optconf_{old_name} opt energy')
    path = tmp_path/'conformer.yaml'
    path.write_text(f'profile: {old_name}\n')
    with pytest.raises(ValueError, match='profile must be low, medium, or high'):
        resolve_options('optconf=xtb', path)


def test_partial_yaml_and_cpu_cap(tmp_path):
    path = tmp_path/'conformer.yaml'
    path.write_text('mm_method: uff\nworkers: 8\nthreads: 2\n')
    cfg = resolve_options('optconf=xtb optconf_medium', str(path))
    assert cfg.profile == 'medium' and cfg.budget.base == 10 and cfg.mm_method == 'uff'
    assert cfg.parallelism(8) == 4
    with pytest.raises(ValueError):
        cfg.parallelism(1)


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


@pytest.mark.parametrize('profile', ['high', 'medium'])
def test_next_stage_first_batch_requests_only_remaining_quota(tmp_path, profile):
    cfg = settings(profile, n=20)
    names = [name for name, _ in cfg.stages()]
    calls = {name: [] for name in names}
    def first(reference, count, *args):
        calls[names[0]].append(count)
        return copies(reference, 18) if len(calls[names[0]]) == 1 else []
    def second(reference, count, *args):
        calls[names[1]].append(count)
        return copies(reference, count)
    out = prepare_candidates(molecule(), tmp_path/'run', cfg, allocated_cores=2,
        generators={names[0]: first, names[1]: second})
    assert calls[names[0]][0] == 20
    assert sum(calls[names[0]]) == 40  # The per-generator cap remains 2N.
    assert calls[names[1]] == [2]
    assert len(out.candidates) == 20
    assert [m.GetProp('generator') for m in out.candidates] == [names[0]]*18 + [names[1]]*2
    assert len(json.loads(out.status_path.read_text())['stages']) == 2


def test_quota_is_carried_across_all_three_stages(tmp_path):
    cfg = settings(n=5)
    calls = {name: [] for name, _ in cfg.stages()}
    def adapter(name, keep_first):
        def generate(reference, count, *args):
            calls[name].append(count)
            return copies(reference, keep_first) if len(calls[name]) == 1 else []
        return generate
    out = prepare_candidates(molecule(), tmp_path/'run', cfg, allocated_cores=2,
        generators={'ditmc': adapter('ditmc', 2),
                    'torsional_diffusion': adapter('torsional_diffusion', 1),
                    'etkdgv3': adapter('etkdgv3', 2)})
    assert [counts[0] for counts in calls.values()] == [5, 3, 2]
    assert sum(calls['ditmc']) == sum(calls['torsional_diffusion']) == 10
    assert calls['etkdgv3'] == [2]
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
    out = prepare_candidates('CCCC', tmp_path/'run', settings('low', n=3, mm_method=method),
                             allocated_cores=2)
    assert out.initial_sdf is not None
    assert 1 <= len(out.candidates) <= 3


def test_empty_model_sdf_is_a_zero_candidate_pool(tmp_path):
    command = [sys.executable, '-c',
               'from pathlib import Path; import sys; Path(sys.argv[1]).touch()', '{output}']
    assert command_generator(molecule(), 2, 42, tmp_path, {'command': command}, 4) == []


def test_parallel_external_commands(tmp_path):
    directory = tmp_path/'batch'
    directory.mkdir()
    code = ('import json,sys,os; from rdkit import Chem; '
            'r=json.load(open(sys.argv[1])); m=next(iter(Chem.SDMolSupplier(sys.argv[2],removeHs=False))); '
            'w=Chem.SDWriter(sys.argv[3]); '
            '[w.write(m) for _ in range(r["maximum_raw_candidates"])]; w.close(); '
            'assert os.environ["OMP_NUM_THREADS"]=="1"')
    command = [sys.executable, '-c', code, '{request}', '{input}', '{output}']
    raw = generate_batch('ditmc', molecule(), 5, 12, directory, {'command': command}, 1, 2)
    assert len(raw) == 5
    counts = [json.loads(p.read_text())['maximum_raw_candidates'] for p in sorted(directory.glob('*/request.json'))]
    assert counts == [3, 2]


def test_cleanup_preserves_search_only_when_requested(tmp_path):
    saved = tmp_path/'conformer_search'
    saved.mkdir()
    (saved/'trace.json').write_text('{}')
    (tmp_path/'optimized_structures.sdf').write_text('test')
    (tmp_path/'junk').mkdir()
    job_cleanup.cleanup_gaussian(tmp_path, preserve_conformers=True)
    assert saved.is_dir() and (tmp_path/'optimized_structures.sdf').is_file()
    assert not (tmp_path/'junk').exists()
    job_cleanup.cleanup_gaussian(tmp_path)
    assert not saved.exists()


def test_continuous_relaxation_audit_and_selection(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    cfg = settings('low', n=3)
    ref = molecule()
    prepared = prepare_candidates(ref, tmp_path/'run', cfg, generators={'etkdgv3': copies})
    def relax(mol, directory, *args):
        index = int(directory.name.split('_')[-1])
        if index == 1:
            raise RuntimeError('synthetic QC failure')
        return Chem.Mol(mol), -10-index
    audit = relax_candidates(prepared, ref, cfg, 0, 1, 'xtb', 2, '1GB', adapter=relax)
    assert audit['primary_best_index'] == 2
    best = next(iter(Chem.SDMolSupplier('optimized_structures.sdf', removeHs=False)))
    assert best.GetProp('candidate_id') == prepared.candidates[2].GetProp('candidate_id')
    assert len(audit['candidate_runs']) == 3


def test_xtb_trace_not_assumed_to_equal_cycles(tmp_path, monkeypatch):
    monkeypatch.setattr('shutil.which', lambda name: '/fake/xtb')
    mol = molecule()
    def run(argv, **kwargs):
        folder = kwargs['cwd']
        Chem.MolToXYZFile(mol, str(folder/'xtbopt.xyz'))
        (folder/'xtbopt.log').write_text('energy: -1.0 gnorm: 0.1\nenergy: -2.0 gnorm: 0.01\n')
        kwargs['stdout'].write('CYCLE 1\nGEOMETRY OPTIMIZATION CONVERGED\nTOTAL ENERGY -2.0 Eh\n')
    monkeypatch.setattr(job_timeout, 'run', run)
    _, energy = xtb_relax(mol, tmp_path, 0, 1, 2, '', {})
    trace = json.loads((tmp_path/'trace.json').read_text())
    assert energy == -2 and trace['optimizer_cycles'] == 1
    assert len(trace['trajectory_records']) == 2


def test_generator_local_timeout():
    import subprocess
    with pytest.raises(subprocess.TimeoutExpired):
        job_timeout.run([sys.executable, '-c', 'import time; time.sleep(2)'], timeout=0.05)


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


@pytest.mark.parametrize('engine', ['gaussian', 'gamess'])
def test_public_runner_resolves_settings_before_chdir(tmp_path, monkeypatch, engine):
    from qcforever.gaussian_run.GaussianRunPack import GaussianDFTRun
    from qcforever.gamess_run.GamessRunPack import GamessDFTRun
    cls = GaussianDFTRun if engine == 'gaussian' else GamessDFTRun
    job = cls.__new__(cls)
    job.value = 'optconf=xtb optconf_medium opt energy uv'
    job.timejob = None
    job.conformer_config = 'custom.yaml'
    monkeypatch.chdir(tmp_path)
    (tmp_path/'custom.yaml').write_text('mm_method: uff\n')
    def worker():
        monkeypatch.chdir(tmp_path.parent)
        return {'profile': job._conformer_settings.profile, 'mm': job._conformer_settings.mm_method}
    setattr(job, '_run_'+engine, worker)
    assert getattr(job, 'run_'+engine)() == {'profile': 'medium', 'mm': 'uff'}
    assert Path.cwd() == tmp_path


def test_laqa_option_is_not_silently_ignored():
    with pytest.raises(NotImplementedError, match='reserved'):
        resolve_options('optconf=xtb laqa opt uv')


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
    audit = relax_candidates(prepared, ref, cfg, 0, 1, 'xtb', 2, '1GB', adapter=backend)
    assert calls == [0, 1, 2, 3, 4]
    assert audit['attempted_candidates'] == 5
    assert audit['converged_candidates'] == 4 and audit['failed_candidates'] == 1
    assert audit['primary_valid_candidates'] == 4 and audit['primary_best_index'] == 4


@pytest.mark.parametrize('nproc', [1, 2, 4])
def test_pm6_existing_input_adapter_and_all_candidate_loop(tmp_path, monkeypatch, nproc):
    from qcforever.laqa_fafoom.pyg16 import g16Object
    import os
    monkeypatch.chdir(tmp_path)
    cfg = settings('low', n=3)
    ref = molecule()
    prepared = prepare_candidates(ref, tmp_path/'run', cfg, generators={'etkdgv3': copies})
    monkeypatch.setattr('shutil.which', lambda name: '/fake/g16')
    from qcforever.conformer_search.relaxation import native_thread_environment
    # An HPC script may still export 4 for model generation. Native PM6 must
    # receive nproc, and the caller's values must survive both success/failure.
    for key in native_thread_environment(4):
        monkeypatch.setenv(key, '4')
    previous = {k: os.environ.get(k) for k in (*native_thread_environment(4), 'GAUSS_EXEDIR', 'GAUSS_SCRDIR')}
    calls = []
    def native_run(self):
        i = len(calls)
        calls.append(i)
        inp = Path('Gau_molecule.com').read_text()
        assert all(os.environ[k] == v for k, v in native_thread_environment(nproc).items())
        assert f'%nprocshared={nproc}' in inp and '%mem=1GB' in inp
        assert 'pm6 opt=(maxcycle=1000)' in inp and '\n0 1\n' in inp
        log = f'SCF Done: E(PM6) = {-10-i}.0 A.U.\n'
        if i != 1:
            log += 'Stationary point found\nNormal termination of Gaussian\n'
        Path('Gau_molecule.log').write_text(log)
        if i == 1:
            raise ValueError('synthetic Gaussian failure')
        self.energy = -10-i
        self.sdf_string_opt = self.sdf_string
    monkeypatch.setattr(g16Object, 'run_g16', native_run)
    audit = relax_candidates(prepared, ref, cfg, 0, 1, 'pm6', nproc, '1GB')
    assert audit['parallel_workers'] == 1 and audit['cores_per_calculation'] == nproc
    assert calls == [0, 1, 2] and audit['primary_best_index'] == 2
    assert audit['failed_candidates'] == 1
    for i in range(3):
        folder = tmp_path/'run'/'electronic'/f'candidate_{i:05d}'
        assert (folder/'Gau_molecule.com').is_file() and (folder/'trace.json').is_file()
    assert {k: os.environ.get(k) for k in previous} == previous


def test_failed_xtb_keeps_partial_trace(tmp_path, monkeypatch):
    import subprocess
    monkeypatch.setattr('shutil.which', lambda name: '/fake/xtb')
    def native_run(argv, **kwargs):
        (kwargs['cwd']/'xtbopt.log').write_text('energy: -1.0 gnorm: 0.1\n')
        kwargs['stdout'].write('CYCLE 1\nfailed\n')
        raise subprocess.CalledProcessError(1, argv)
    monkeypatch.setattr(job_timeout, 'run', native_run)
    with pytest.raises(subprocess.CalledProcessError):
        xtb_relax(molecule(), tmp_path, 0, 1, 2, '', {})
    trace = json.loads((tmp_path/'trace.json').read_text())
    assert not trace['native_converged'] and len(trace['trajectory_records']) == 1


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


def test_ez_mismatch_rejected_at_generation():
    ref = Chem.AddHs(Chem.MolFromSmiles('F/C=C/F'))
    assert AllChem.EmbedMolecule(ref, randomSeed=42) == 0
    wrong = Chem.Mol(ref)
    rdMolTransforms.SetDihedralDeg(wrong.GetConformer(), 0, 1, 2, 3, 0)
    kept, audit = filter_candidates([wrong], ref, 2, SearchConfig.resolve().validation)
    assert not kept and audit['rejected'] == {'stereo_mismatch': 1}
    assert not audit['candidates'][0]['ez_stereo_match']


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


@pytest.mark.parametrize('nproc', [1, 4, 6, 8, 16, 32])
def test_native_relaxation_uses_four_core_workers(tmp_path, monkeypatch, nproc):
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
    per_call = min(4, nproc)
    workers = min(2, nproc // per_call)
    assert audit['parallel_workers'] == workers and audit['cores_per_calculation'] == per_call
    assert audit['scheduling'] == 'parallel_candidates_4_cores'
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


def test_no_primary_still_attempts_every_candidate(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    cfg = settings('low', n=3)
    ref = molecule()
    prepared = prepare_candidates(ref, tmp_path/'run', cfg, generators={'etkdgv3': copies})
    def fails(*args):
        raise RuntimeError('synthetic failure')
    with pytest.raises(RuntimeError, match='No converged'):
        relax_candidates(prepared, ref, cfg, 0, 1, 'xtb', 2, '', adapter=fails)
    audit = json.loads((tmp_path/'run'/'electronic'/'audit.json').read_text())
    assert audit['attempted_candidates'] == 3 and audit['failed_candidates'] == 3


def test_global_timeout_is_not_a_candidate_failure(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    cfg = settings('low', n=3)
    ref = molecule()
    prepared = prepare_candidates(ref, tmp_path/'run', cfg, generators={'etkdgv3': copies})
    def expired(*args):
        raise job_timeout.QCforeverTimeoutError('synthetic deadline')
    with pytest.raises(job_timeout.QCforeverTimeoutError):
        relax_candidates(prepared, ref, cfg, 0, 1, 'xtb', 2, '', adapter=expired)
    progress = json.loads((tmp_path/'run'/'electronic'/'progress.json').read_text())
    assert len(progress) == 1 and progress[0]['state'] == 'timeout'


def test_xtb_subprocess_and_new_optconf_entry_end_to_end(tmp_path, monkeypatch):
    from qcforever.conformer_search.conformer_search import configured_confopt
    from qcforever.conformer_search.pipeline import write_sdf
    monkeypatch.chdir(tmp_path)
    binary = tmp_path/'fake_xtb'
    binary.write_text('#!'+sys.executable+'\n'+
        'import pathlib,sys\n'
        'p=pathlib.Path.cwd(); args=sys.argv[1:]\n'
        'assert args[args.index("--gfn")+1]=="2"\n'
        'assert args[args.index("--chrg")+1]=="0"\n'
        'assert args[args.index("--uhf")+1]=="0"\n'
        'assert args[args.index("--parallel")+1]=="2"\n'
        'i=int(p.name.split("_")[-1]); energy=-10-i\n'
        '(p/"xtbopt.xyz").write_text((p/"input.xyz").read_text())\n'
        '(p/"xtbopt.log").write_text(f"energy: {energy} gnorm: 0.0001\\n")\n'
        'print("CYCLE 1\\nGEOMETRY OPTIMIZATION CONVERGED")\n'
        'print(f"TOTAL ENERGY {energy} Eh")\n')
    binary.chmod(0o755)
    inp = tmp_path/'molecule.sdf'
    write_sdf(inp, [molecule()])
    cfg = settings('low', n=3, relaxation={'xtb_executable': str(binary)})
    summary = configured_confopt(str(inp), 0, 1, 'xtb', 2, '1GB', config=cfg)
    assert summary['attempted_candidates'] == summary['converged_candidates'] == 3
    assert summary['energy_hartree'] == -12
    assert summary['backend'] == 'xtb' and summary['failed_candidates'] == 0
    root = tmp_path/'conformer_search'
    assert len(list((root/'electronic').glob('candidate_*/trace.json'))) == 3
    best = next(iter(Chem.SDMolSupplier('optimized_structures.sdf', removeHs=False)))
    assert best.GetProp('candidate_id') == summary['selected_candidate_id']
    assert json.loads((root/'summary.json').read_text()) == summary
