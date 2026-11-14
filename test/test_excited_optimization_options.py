"""Keep ground- and excited-state optimization controls independent."""

from pathlib import Path
import re

import pytest

from qcforever.conformer_search.settings import calculation_tokens, parse_conformer_options
from qcforever.util.calculation_options import parse_excited_optimization, excited_options
from qcforever.util import check_resource
from qcforever.gaussian_run import GaussianRunPack, Exe_Gaussian, chk2fchk
from qcforever.gamess_run import GamessRunPack
from qcforever.orca_run import OrcaRunPack


@pytest.mark.parametrize('token,state,fc', [
    ('fluor', 1, ''), ('fluor=3', 3, ''), ('fluor=-1', -1, ''),
    ('fluor=(FC=calcFC, state=2)', 2, 'CalcFC'),
    ('FLUOR=( STATE = 3 , fc = CALCALL )', 3, 'CalcAll'),
    ('fluor=(state=2)', 2, ''), ('fluor=()', 1, ''),
    ('fluor=(FC=none)', 1, ''), ('fluor=(FC=default)', 1, ''),
    ('fluor=(FC=RecalcFC=5, state=2)', 2, 'RecalcFC=5'),
    ('tadf=(FC=ReadFC,state=2)', 2, 'ReadFC'),
])
def test_excited_group_parser(token,state,fc):
    assert parse_excited_optimization(token)=={'state':state,'fc':fc}


@pytest.mark.parametrize('token', [
    'fluor=', 'fluor=0', 'fluor=-2', 'fluor=1.5', 'fluor=(state=1 2)',
    'fluor=(FC=)', 'fluor=(state=)', 'fluor=(state=1,state=2)',
    'fluor=(FC=CalcFC,fc=CalcAll)', 'fluor=(method=CalcFC)',
    'fluor=(FC=bad)', 'fluor=(FC=RecalcFC=0)', 'fluor=(FC=CalcFC,)', 'tadf=-1',
])
def test_invalid_state_settings_are_rejected(token):
    with pytest.raises(ValueError):
        parse_excited_optimization(token)


@pytest.mark.parametrize('text', ['fluor=(state=1', 'fluor=(FC=CalcFC,state=1))',
                                 'fluor=(state=1)energy'])
def test_malformed_groups_are_rejected(text):
    with pytest.raises(ValueError):
        calculation_tokens(text)


def test_property_and_conformer_groups_coexist_without_altering_paths():
    text='optconf=(conf_gen=low, conf_relax=xtb) opt=(CalcAll, Tight) fluor=(FC=CalcFC, state=2) uv=/tmp/a(b).dat freq=IR.dat,Raman.dat'
    tokens=calculation_tokens(text)
    assert tokens==['optconf=xtb','opt=(CalcAll, Tight)','fluor=(FC=CalcFC, state=2)',
                    'uv=/tmp/a(b).dat','freq=IR.dat,Raman.dat']
    assert parse_conformer_options(text).profile=='low'


def test_duplicate_excited_groups_are_rejected():
    with pytest.raises(ValueError,match='fluor only once'):
        excited_options(['fluor=1','fluor=(state=2)'])


@pytest.fixture
def gaussian_workflow(tmp_path,monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(check_resource,'respec_cores',lambda value:value)
    monkeypatch.setattr(check_resource,'respec_memory',lambda **kwargs:'1GB')
    Path('water.xyz').write_text('3\n0 1\nO 0 0 0\nH 0 0 .96\nH .92 0 -.24\n')
    calls=[]
    def launch(name,*args):
        calls.append((name,Path(name+'.com').read_text()))
        Path(name+'.log').write_text('Normal termination of Gaussian 16\n')
        return 'normal'
    monkeypatch.setattr(Exe_Gaussian,'exe_Gaussian',launch)
    monkeypatch.setattr(chk2fchk,'Get_chklist',lambda *args:None)
    def extract(self,name,options,*args):
        if 'fluor' in options:
            return {'MinEtarget':-74.0,'fluor':[[400],[.1]]}
        if 'tadf' in options:
            return {'T_Min':-74.1}
        return {'Energy':[-75.0,0.0], 'uv':[[250,200,150],[.1,.3,.05]],
                'state_index':[[1,2,3],[]]}
    monkeypatch.setattr(GaussianRunPack.GaussianDFTRun,'Extract_values',extract)
    def run(options):
        calls.clear()
        job=GaussianRunPack.GaussianDFTRun('B3LYP','STO-3G',1,options,'water.xyz')
        result=job.run_gaussian()
        assert result['log']=='normal',result
        return calls.copy()
    return run


@pytest.mark.parametrize('excited', ['fluor', 'fluor=2', 'tadf'])
def test_ground_calcall_never_leaks_to_excited_states(gaussian_workflow,excited):
    calls=gaussian_workflow('opt=calcall energy '+excited)
    assert 'opt=(calcall,' in calls[0][1].lower()
    assert all('calcall' not in text.lower() for _,text in calls[1:])


def test_force_constants_and_roots_are_independent_for_ground_singlet_and_triplet(gaussian_workflow):
    calls=gaussian_workflow('opt=calcall energy fluor=(FC=CalcFC, state=2) tadf=(state=1, FC=ReadFC)')
    assert len(calls)==3
    ground,singlet,triplet=calls
    assert 'opt=(calcall,' in ground[1].lower()
    assert '_State2_1_0' in singlet[0] and 'Opt=(CalcFC,' in singlet[1]
    assert '_State1_3_0' in triplet[0] and 'Opt=(ReadFC,' in triplet[1]
    assert 'root=2' in singlet[1] and 'root=1' in triplet[1]
    assert 'calcall' not in (singlet[1]+triplet[1]).lower()
    assert 'CalcFC' not in ground[1] and 'ReadFC' not in ground[1]


def test_tadf_group_only_applies_fc_to_triplet_and_keeps_implicit_s1_default(gaussian_workflow):
    calls=gaussian_workflow('opt=calcall energy tadf=(FC=CalcFC,state=2)')
    assert '_State1_1_0' in calls[1][0]
    assert 'CalcFC' not in calls[1][1] and 'calcall' not in calls[1][1].lower()
    assert '_State2_3_0' in calls[2][0] and 'Opt=(CalcFC,' in calls[2][1]


def test_fc_without_ground_opt_and_out_of_order_tadf_fluor(gaussian_workflow):
    calls=gaussian_workflow('tadf=(FC=CalcAll,state=1) energy fluor=(state=3,FC=CalcFC)')
    assert not re.search(r'\bOpt=',calls[0][1],re.I)
    assert '_State3_1_0' in calls[1][0] and 'Opt=(CalcFC,' in calls[1][1]
    assert '_State1_3_0' in calls[2][0] and 'Opt=(CalcAll,' in calls[2][1]


def test_nac_is_attached_to_fluor_even_after_tadf(gaussian_workflow):
    calls=gaussian_workflow('opt=calcall fluor=(FC=CalcFC,state=2) tadf nac')
    assert ',NAC)' in calls[1][1]
    assert ',NAC)' not in calls[0][1]+calls[2][1]


@pytest.fixture
def gaussian_dipole(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(check_resource, 'respec_cores', lambda value: value)
    monkeypatch.setattr(check_resource, 'respec_memory', lambda **kwargs: '1GB')
    Path('water.xyz').write_text('3\n0 1\nO 0 0 0\nH 0 0 .96\nH .92 0 -.24\n')
    calls = []
    invalid_density = set()

    def launch(name, *args):
        route = Path(name + '.com').read_text()
        calls.append((name, route))
        log = 'Population analysis using the SCF Density.\n'
        log += 'Dipole moment (field-independent basis, Debye):\n X=1.0 Y=0.0 Z=0.0 Tot=1.0\n'
        if 'Opt=' in route:
            log += ' Optimization completed.\n'
        if 'Density=Current' in route and name not in invalid_density:
            root = int(re.search(r'root=(\d+)', route)[1])
            log += (f' Excited State   {root}: Singlet-A\n'
                    ' This state for optimization and/or second-order correction.\n'
                    'Population analysis using the CI density.\n'
                    'Dipole moment (field-independent basis, Debye):\n X=2.0 Y=0.0 Z=0.0 Tot=2.0\n'
                    'Population analysis using the SCF Density.\n'
                    'Dipole moment (field-independent basis, Debye):\n X=1.0 Y=0.0 Z=0.0 Tot=1.0\n')
        Path(name + '.log').write_text(log + 'Normal termination of Gaussian 16\n')
        Path(name + '.chk').write_text(f'checkpoint for {name}')
        return 'normal'

    def format_checkpoints(remove, *, keep=None):
        for path in Path.cwd().glob('*.chk'):
            path.with_suffix('.fchk').write_text(path.read_text())
            if remove and path.name != keep:
                path.unlink()

    def from_fchk(source, destination):
        Path(destination).write_bytes(Path(source).read_bytes())

    def restore_checkpoint(name):
        from_fchk(name + '.fchk', name + '.chk')
        return 0, 1

    monkeypatch.setattr(Exe_Gaussian, '_exe_Gaussian', launch)
    monkeypatch.setattr(chk2fchk, 'Get_chklist', format_checkpoints)
    monkeypatch.setattr(GaussianRunPack.gaussian_run.fchk2chk, 'convert', from_fchk)
    monkeypatch.setattr(GaussianRunPack.gaussian_run.fchk2chk, 'Get_fchk', restore_checkpoint)
    monkeypatch.setattr(GaussianRunPack.GaussianDFTRun, 'Extract_values', lambda *args:
                        {'state_index': [[1, 2], []], 'fluor': [[400], [.1]], 'NAC': [.4, .1]})

    def run(options, solvent=None, filename='water.xyz'):
        job = GaussianRunPack.GaussianDFTRun('B3LYP', 'STO-3G', 1, options, filename,
                                             solvent=solvent, pklsave=True)
        return job.run_gaussian()

    return run, calls, invalid_density


def test_excited_dipole_at_saved_input_retries_missing_relaxed_density(gaussian_dipole):
    run, calls, invalid_density = gaussian_dipole
    option = 'dipole=(state=1, geom=input)'
    invalid_density.add('water_Dipole_input_S1')
    first = run(option)
    assert first['log_1'] == 'error'
    assert 'CI density' in first['error_1']
    assert 'excited_dipole' not in first
    assert first['reuse_report']['failed'] == 1
    invalid_density.clear()
    result = run(option)
    assert result['log'] == 'normal'
    assert result['reuse_report']['reused'] == 1
    assert result['excited_dipole']['ground'] == [1, 0, 0, 1]
    assert result['excited_dipole']['excited'] == [2, 0, 0, 2]
    assert result['excited_dipole']['solvent_response'] == 'vacuum'
    assert len(calls) == 3
    assert all('H       0.92000' in route and 'Opt=' not in route for _, route in calls)
    assert 'TD(NStates=5,Singlet,root=1)' in calls[1][1]


def test_excited_dipole_reuses_supplied_formatted_checkpoint(gaussian_dipole):
    import pickle

    run, calls, _ = gaussian_dipole
    checkpoint = Path('water.fchk')
    checkpoint.write_text('Charge     I   1\nMultiplicity      I   1\n')
    result = run('dipole=(state=1)', filename=checkpoint)
    assert result['log'] == 'normal', result
    assert len(calls) == 2
    for _, route in calls:
        assert '%Oldchk=water_Dipole_input_seed.chk' in route
        assert 'Geom=CheckPoint Guess=Read' in route
        assert re.search(r'^\s*1[ \t]+1[ \t]*$', route, re.M)
    assert pickle.loads(Path('water/water.pkl').read_bytes())['excited_dipole'] == result['excited_dipole']
    assert checkpoint.is_file() and not Path('water/water_Dipole_input_seed.chk').exists()


def test_excited_dipole_adds_independent_equilibrium_stages_after_fluor(gaussian_dipole):
    run, calls, _ = gaussian_dipole
    base = run('fluor=2 nac', solvent='Water')
    before = len(calls)
    result = run('dipole=(state=2, geom=fluor) fluor=2 nac', solvent='Water')
    assert result['log'] == 'normal', result
    assert result['fluor'] == base['fluor'] and result['NAC'] == base['NAC']
    assert result['reuse_report']['reused'] == base['reuse_report']['executed']
    assert all('Density=Current' not in route and 'EqSolv' not in route for _, route in calls[:before])
    ground, excited = [route for _, route in calls[before:]]
    assert 'Density=SCF' in ground
    assert 'Density=Current' in excited
    assert 'TD(NStates=6,Singlet,root=2,EqSolv)' in excited
    for route in (ground, excited):
        assert '%Oldchk=water_State2_1_0.chk' in route
        assert 'Geom=CheckPoint Guess=Read' in route
        assert 'SCRF=(PCM, solvent=Water)' in route
        assert 'Opt=' not in route
    assert result['excited_dipole']['geometry'] == 'fluor'
    assert result['excited_dipole']['solvent_response'] == 'equilibrium'


@pytest.mark.parametrize('options', [
    'dipole=(state=0, geom=input)', 'dipole=(state=1, geom=opt)',
    'dipole=(state=1, geom=input) dipole=(state=2, geom=input)',
    'dipole=(state=1, geom=fluor)', 'fluor=2 dipole=(state=1, geom=fluor)',
])
def test_invalid_excited_dipole_options_fail_before_running(gaussian_dipole, options):
    run, calls, _ = gaussian_dipole
    with pytest.raises(ValueError):
        run(options)
    assert not calls


@pytest.mark.parametrize('option', ['fluor=-1', 'fluor=(state=-1)',
                                  'fluor=(FC=CalcFC,state=-1)'])
def test_legacy_auto_state_selection_is_preserved(gaussian_workflow, option):
    calls=gaussian_workflow('opt=calcall ' + option)
    assert '_State2_1_0' in calls[1][0]
    assert 'calcall' not in calls[1][1].lower()
    assert ('Opt=(CalcFC,' in calls[1][1]) == ('FC=CalcFC' in option)


def test_changing_only_excited_fc_keeps_ground_recipe_identical(gaussian_workflow):
    first=gaussian_workflow('opt=calcall energy fluor=(FC=CalcFC,state=1)')
    second=gaussian_workflow('opt=calcall energy fluor=(FC=CalcAll,state=1)')
    assert first[0][1]==second[0][1]
    assert first[1][1]!=second[1][1]


def test_equal_signs_in_fc_values_and_ground_options_are_preserved(gaussian_workflow):
    calls=gaussian_workflow('opt=RecalcFC=5 fluor=(state=2,FC=RecalcFC=3)')
    assert 'Opt=(RecalcFC=5,' in calls[0][1]
    assert 'Opt=(RecalcFC=3,' in calls[1][1]


@pytest.mark.parametrize('backend',['orca','gamess'])
def test_other_backends_accept_portable_fc_and_reject_gaussian_rcfc(backend):
    assert parse_excited_optimization('fluor=(state=2)',backend)=={'state':2,'fc':''}
    assert parse_excited_optimization('fluor=(FC=CalcFC,state=2)',backend)=={'state':2,'fc':'CalcFC'}
    with pytest.raises(ValueError,match='RCFC'):
        parse_excited_optimization('fluor=(FC=RCFC,state=2)',backend)


def test_orca_real_parser_accepts_grouped_state(tmp_path,monkeypatch):
    monkeypatch.setattr(check_resource,'respec_cores',lambda value:value)
    job=OrcaRunPack.OrcaDFTRun('B3LYP','STO-3G',1,'opt=tight fluor=(state=2)',tmp_path/'water.xyz')
    assert job._options()['fluor']==2


def test_gamess_workflow_parses_legacy_and_grouped_root(tmp_path,monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(check_resource,'respec_cores',lambda value:value)
    monkeypatch.setattr(check_resource,'respec_memory',lambda **kwargs:'1GB')
    Path('water.xyz').write_text('3\n0 1\nO 0 0 0\nH 0 0 .96\nH .92 0 -.24\n')
    calls=[]
    def stop(run_type,charge,mult,filename,*args,**kwargs):
        calls.append((run_type,kwargs.get('target')))
        if kwargs.get('target') is not None:
            raise RuntimeError('input inspected')
    for option in ('fluor=2','fluor=(state=2)'):
        job=GamessRunPack.GamessDFTRun('B3LYP','STO-3G',1,option,'water.xyz')
        monkeypatch.setattr(job,'make_input',stop)
        monkeypatch.setattr(GamessRunPack.gamess_run.Exe_Gamess,'exe_Gamess',lambda *args:'')
        monkeypatch.setattr(GamessRunPack.job_cleanup,'cleanup_gamess',lambda *args,**kwargs:None)
        with pytest.raises(RuntimeError,match='input inspected'):
            job.run_gamess()
        assert calls[-1]==('OPTIMIZE',[2,1])
