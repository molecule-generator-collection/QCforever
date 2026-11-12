"""Reuse contracts across all engines, without installed Gaussian or GAMESS."""

from contextlib import contextmanager
import json
from pathlib import Path
import runpy
from types import SimpleNamespace

import pytest

from qcforever.util import calculation_cache as cache
from qcforever.gaussian_run import Exe_Gaussian, GaussianRunPack
from qcforever.gamess_run import Exe_Gamess, GamessRunPack
from qcforever.orca_run import OrcaRunPack, Exe_Orca


ORCA_LOG = Path(__file__).parent / 'fixtures/orca/water.txt'
TD_LOG = Path(__file__).parent / 'fixtures/orca/tddft.txt'


@contextmanager
def session(directory, backend, enabled=True):
    value = cache.CalculationCache(backend, enabled)
    token = cache._active.set(value)
    original = Path.cwd()
    try:
        value.prepare(directory)
        __import__('os').chdir(directory)
        yield value
    finally:
        __import__('os').chdir(original)
        value.finish()
        cache._active.reset(token)


def example(backend):
    if backend == 'gaussian':
        return '%chk=mol\n%mem=1GB\n#p b3lyp/sto-3g\n\ntitle\n\n0 1\nH 0 0 0\nH 0 0 .74\n\n'
    if backend == 'gamess':
        return ' $CONTRL RUNTYP=ENERGY DFTTYP=B3LYP ICHARG=0 MULT=1 $END\n $SYSTEM MWORDS=125 TIMLIM=60 $END\n $DATA\nmol\nC1\nH 1 0 0 0\nH 1 0 0 .74\n $END\n'
    return '! B3LYP STO-3G\n%maxcore 750\n%pal nprocs 1 end\n* xyz 0 1\nH 0 0 0\nH 0 0 .74\n*\n'


def complete(backend, name='mol'):
    if backend == 'gaussian':
        Path(name + '.log').write_text(' Initial command:\nSCF Done: E(RB3LYP) = -1.0\nNormal termination of Gaussian 16\n')
        Path(name + '.chk').write_text('verified orbitals')
        return 'normal'
    if backend == 'gamess':
        Path(name + '.log').write_text('EXECUTION OF GAMESS TERMINATED NORMALLY\n')
        Path(name + '.dat').write_text(' $DATA\nmol\nC1\nH 1.0 0 0 0\nH 1.0 0 0 .74\n $END\n $VEC\n 1 1 .1 .2\n 2 1 .2 .1\n $END\n')
        return ''
    Path(name + '.log').write_text(ORCA_LOG.read_text())
    return 'normal'


@pytest.mark.parametrize('backend', ['gaussian', 'gamess', 'orca'])
def test_matching_stage_is_reused_and_old_files_are_preserved(tmp_path, backend):
    folder = tmp_path / 'job'
    calls = []
    def execute():
        calls.append(1)
        return complete(backend)
    with session(folder, backend) as first:
        cache.execute(backend, 'mol', example(backend), execute)
        Path('notes.txt').write_text('user data')
    with session(folder, backend) as second:
        cache.execute(backend, 'mol', example(backend), execute)
        assert second.steps[0]['action'] == 'reused'
        assert list(Path('.qcforever/history').glob('*/notes.txt'))[0].read_text() == 'user data'
    assert len(calls) == 1
    assert json.loads((folder/'.qcforever/reuse-report.json').read_text())['reused'] == 1


@pytest.mark.parametrize('backend', ['gaussian', 'gamess', 'orca'])
@pytest.mark.parametrize('change', ['method', 'geometry', 'state', 'solvent', 'convergence'])
def test_changed_physical_or_numerical_conditions_execute_again(tmp_path, backend, change):
    original = example(backend)
    changed = {
        'method': original.replace('B3LYP', 'PBE0').replace('b3lyp', 'pbe0'),
        'geometry': original.replace('.74', '.75'),
        'state': original.replace('0 1', '1 2').replace('ICHARG=0 MULT=1', 'ICHARG=1 MULT=2'),
        'solvent': original + '\nsolvent=water\n',
        'convergence': original + '\nSCF=tight\n',
    }[change]
    calls = []
    def run():
        calls.append(1)
        return complete(backend)
    with session(tmp_path/'job', backend):
        cache.execute(backend, 'mol', original, run)
    with session(tmp_path/'job', backend) as second:
        cache.execute(backend, 'mol', changed, run)
        assert second.steps[0]['action'] == 'executed'
    assert len(calls) == 2


@pytest.mark.parametrize('backend', ['gaussian', 'gamess', 'orca'])
def test_allocation_changes_do_not_recompute(tmp_path, backend):
    text = example(backend)
    calls = []
    def run():
        calls.append(1)
        return complete(backend)
    with session(tmp_path/'job', backend):
        cache.execute(backend, 'mol', text, run)
    changed = text.replace('%mem=1GB', '%mem=4GB').replace('MWORDS=125 TIMLIM=60', 'MWORDS=500 TIMLIM=1')
    changed = changed.replace('%maxcore 750', '%maxcore 1000').replace('nprocs 1', 'nprocs 2')
    with session(tmp_path/'job', backend):
        cache.execute(backend, 'mol', changed, run)
    assert len(calls) == 1


@pytest.mark.parametrize('backend', ['gaussian', 'gamess', 'orca'])
def test_incomplete_log_is_never_cached(tmp_path, backend):
    calls = []
    def run():
        calls.append(1)
        complete(backend)
        Path('mol.log').write_text('SCF failure; interrupted\n')
        return '' if backend == 'gamess' else 'normal'
    for _ in range(2):
        with session(tmp_path/'job', backend) as current:
            cache.execute(backend, 'mol', example(backend), run)
            assert current.steps[0]['action'] == 'failed'
    assert len(calls) == 2


@pytest.mark.parametrize('backend,suffix', [('gaussian','.chk'), ('gamess','.dat'), ('orca','.log')])
def test_missing_or_modified_cached_output_is_recomputed(tmp_path, backend, suffix):
    calls = []
    def run():
        calls.append(1)
        return complete(backend)
    with session(tmp_path/'job', backend) as first:
        cache.execute(backend, 'mol', example(backend), run)
        entry = next(iter(first.index['entries'].values()))
        (first.state / entry['snapshot'] / ('mol'+suffix)).write_text('modified output')
    with session(tmp_path/'job', backend) as second:
        cache.execute(backend, 'mol', example(backend), run)
        assert second.steps[0]['action'] == 'executed'
    assert len(calls) == 2


def test_gaussian_checkpoint_dependency_invalidates_downstream(tmp_path):
    folder=tmp_path/'job'
    text='%chk=mol\n#p b3lyp/sto-3g geom=allcheck guess=read\n\n'
    calls=[]
    def run():
        calls.append(1)
        return complete('gaussian')
    for content in ('first geometry', 'changed geometry'):
        with session(folder,'gaussian'):
            Path('mol.chk').write_text(content)
            cache.execute('gaussian','mol',text,run)
    assert len(calls)==2


def test_orca_recomputes_a_cached_step_when_its_orbitals_are_required(tmp_path):
    folder = tmp_path / 'job'
    calls = []

    def run():
        calls.append(1)
        if len(calls) == 2:
            Path('mol.gbw').write_bytes(b'converged orbitals')
        return complete('orca')

    with session(folder, 'orca'):
        cache.execute('orca', 'mol', example('orca'), run)
    for action in ('executed', 'reused'):
        with session(folder, 'orca') as current:
            cache.execute('orca', 'mol', example('orca'), run, required_outputs=('mol.gbw',))
            assert current.steps[0]['action'] == action
            assert Path('mol.gbw').read_bytes() == b'converged orbitals'
    assert len(calls) == 2


def test_reuse_false_forces_fresh_work(tmp_path):
    calls=[]
    def run():
        calls.append(1)
        return complete('orca')
    for enabled in (True,False):
        with session(tmp_path/'job','orca',enabled):
            cache.execute('orca','mol',example('orca'),run)
    assert len(calls)==2


def test_gaussian_link1_completed_prefix_and_failure_recovery(tmp_path,monkeypatch):
    original=Path.cwd()
    monkeypatch.chdir(tmp_path)
    calls=[]
    def fake(name,*args):
        text=Path(name+'.com').read_text()
        calls.append(text)
        return complete('gaussian',name)
    monkeypatch.setattr(Exe_Gaussian,'_exe_Gaussian',fake)
    base=example('gaussian')
    freq='%chk=mol\n#p b3lyp/sto-3g freq=raman geom=allcheck guess=read\n\n'
    for text in (base,base+'--Link1--\n'+freq,base+'--Link1--\n'+freq):
        with session(tmp_path/'job','gaussian'):
            Path('mol.com').write_text(text)
            assert Exe_Gaussian.exe_Gaussian('mol',60)=='normal'
            assert Path('mol.com').read_text()==text
    assert len(calls)==2
    assert len(list((tmp_path/'job/.qcforever/history').iterdir()))==2


@pytest.mark.parametrize('solvent', ['0', 'water'])
@pytest.mark.parametrize('separator', ['--Link1--\n', ' \t--lInK1-- \t\n'])
def test_gaussian_link1_preserves_required_blank_lines_and_reuses_stages(
        tmp_path, monkeypatch, solvent, separator):
    monkeypatch.setattr(GaussianRunPack.check_resource, 'respec_cores', lambda count: count)
    job = GaussianRunPack.GaussianDFTRun('B3LYP', 'STO-3G', 1,
                                       'opt symm freq', 'water.xyz', solvent=solvent)
    job.mem = '2GB'
    calls = []
    def execute(name, *args):
        text = Path(name + '.com').read_text()
        # A coordinate line followed by EOF is rejected by Gaussian l101 as
        # "End of file in ZSymb." The blank section terminator must survive.
        assert text.endswith('\n\n'), repr(text[-100:])
        calls.append(text)
        return complete('gaussian', name)
    monkeypatch.setattr(Exe_Gaussian, '_exe_Gaussian', execute)
    for attempt in range(2):
        with session(tmp_path / 'job', 'gaussian') as current:
            job.chain_job('mol', True, {'opt': True, 'symm': True, 'freq': True},
                          0, 1, 0, 'xyz', element=['O', 'H', 'H'],
                          atomX=[0, 0, .92], atomY=[0, 0, 0], atomZ=[0, .96, -.24])
            original = Path('mol.com').read_text()
            expected = original.split('--Link1--\n')
            submitted = original.replace('--Link1--\n', separator)
            Path('mol.com').write_text(submitted)
            assert Exe_Gaussian.exe_Gaussian('mol', 60) == 'normal'
            assert Path('mol.com').read_text() == submitted
            assert calls == expected
            assert Path('mol.log').read_text().count('Normal termination') == 3
            assert {step['action'] for step in current.steps} == ({'executed'} if attempt == 0 else {'reused'})


@pytest.mark.parametrize('restart', [True, False])
@pytest.mark.parametrize('status', ['normal', 'abnormal', 'timeout', 'oscillation'])
def test_gaussian_checkpoint_retention_and_job_status(tmp_path, monkeypatch, status, restart):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(GaussianRunPack.check_resource, 'respec_cores', lambda count: count)
    Path('water.xyz').write_text('3\n0 1\nO 0 0 0\nH 0 0 .96\nH .92 0 -.24\n')
    Path('water').mkdir()
    Path('water/previous.FCHK').write_text('previous checkpoint')
    def execute(name, *args):
        if status == 'normal':
            complete('gaussian', name)
            Path(name + '.fchk').write_text('formatted checkpoint')
            return status
        Path(name + '.log').write_text('End of file in ZSymb.\nError termination via Lnk1e\n')
        Path(name + '.chk').write_text('partial checkpoint')
        Path(name + '.fchk').write_text('partial formatted checkpoint')
        return status
    monkeypatch.setattr(Exe_Gaussian, '_exe_Gaussian', execute)
    if status == 'normal':
        monkeypatch.setattr(GaussianRunPack.gaussian_run.chk2fchk, 'Get_chklist', lambda *args: None)
        monkeypatch.setattr(GaussianRunPack.gaussian_run.Get_MolCoordinate_fchk, 'Get_fchklist2xyz',
                            lambda: Path('optimized.xyz').write_text('geometry'))
        monkeypatch.setattr(GaussianRunPack.GaussianDFTRun, 'Extract_values', lambda *args: {'symm': 'C1'})
    else:
        monkeypatch.setattr(GaussianRunPack.gaussian_run.chk2fchk, 'Get_chklist',
                            lambda *args: pytest.fail('Do not convert a failed checkpoint'))
        monkeypatch.setattr(GaussianRunPack.gaussian_run.Get_MolCoordinate_fchk, 'Get_fchklist2xyz',
                            lambda *args: pytest.fail('Do not export failed geometry'))
        monkeypatch.setattr(GaussianRunPack.GaussianDFTRun, 'Extract_values',
                            lambda *args: pytest.fail('Do not parse properties from a failed calculation'))
    job = GaussianRunPack.GaussianDFTRun('B3LYP', 'STO-3G', 1,
                                       'opt symm homolumo energy', 'water.xyz', restart=restart, pklsave=True)
    job.mem = '2GB'
    result = job.run_gaussian()
    assert result['log_0'] == status
    if status == 'normal':
        assert result['log'] == 'normal' and result['symm'] == 'C1'
        assert Path('water/optimized.xyz').exists() == (not restart)
        assert result['reuse_report']['executed'] == 2
    else:
        assert result['log'] != 'normal'
        assert status in result['error_0'] and result['logfile'] in result['error_0']
        assert 'ZSymb' in Path(result['logfile']).read_text()
        assert result['reuse_report']['failed'] == 1
    checkpoints = [p for p in Path('water').rglob('*') if p.suffix.lower() in {'.chk', '.fchk'}]
    if restart:
        assert {p.name for p in checkpoints} >= {'water_State0_1_0.chk', 'water_State0_1_0.fchk', 'previous.FCHK'}
    else:
        assert not checkpoints
    assert list(Path('water').glob('*.log')) and Path('water/water.pkl').is_file()


@pytest.mark.parametrize('result,reason', [
    ({'log': 'Gaussian computation failed!', 'error_0': 'End of file in ZSymb.',
      'logfile': '/work/ch2o/ch2o_State0_1_0.log'}, 'ZSymb'),
    ({'log': 'normal'}, 'requested symm result'),
    ({'log': 'normal', 'symm': 'C2V'}, 'requested triplet results: T_Min'),
])
def test_manual_gaussian_comparison_reports_missing_result_before_keyerror(
        monkeypatch, result, reason):
    monkeypatch.setattr(GaussianRunPack, 'GaussianDFTRun',
                        lambda *args, **kwargs: SimpleNamespace(run_gaussian=lambda: result))
    with pytest.raises(RuntimeError, match=reason):
        runpy.run_path(str(Path(__file__).with_name('testcal.py')))


def test_legacy_orca_requires_executed_input_echo_and_success(tmp_path):
    folder=tmp_path/'job'; folder.mkdir()
    text=example('orca')
    echo='\n'.join(f'| {i+1}> {line}' for i,line in enumerate(text.splitlines()))
    (folder/'mol.log').write_text(echo+'\n'+ORCA_LOG.read_text())
    with session(folder,'orca') as current:
        cache.execute('orca','mol',text,lambda: pytest.fail('verified legacy job should be reused'))
        assert current.steps[0]['action']=='reused'


def test_gaussian_failed_link_is_retried_without_repeating_completed_link(tmp_path,monkeypatch):
    calls=[]
    failed=[False]
    def fake(name,*args):
        text=Path(name+'.com').read_text()
        calls.append('freq' if 'freq=raman' in text else 'base')
        if 'freq=raman' in text and not failed[0]:
            failed[0]=True
            Path(name+'.log').write_text(' Initial command:\nError termination\n')
            return 'abnormal'
        return complete('gaussian',name)
    monkeypatch.setattr(Exe_Gaussian,'_exe_Gaussian',fake)
    text=example('gaussian')+'--Link1--\n%chk=mol\n#p b3lyp/sto-3g freq=raman geom=allcheck guess=read\n\n'
    for expected in ('abnormal','normal'):
        with session(tmp_path/'job','gaussian'):
            Path('mol.com').write_text(text)
            assert Exe_Gaussian.exe_Gaussian('mol',60)==expected
    assert calls==['base','freq','freq']
    assert 'Error termination' not in (tmp_path/'job/mol.log').read_text()


@pytest.mark.parametrize('backend',['gaussian','gamess'])
def test_real_workflow_reuses_base_and_adds_only_new_stage(tmp_path,monkeypatch,backend):
    monkeypatch.chdir(tmp_path)
    Path('mol.xyz').write_text('2\n0 1\nH 0 0 0\nH 0 0 .74\n')
    monkeypatch.setattr(OrcaRunPack.check_resource,'respec_cores',lambda count:count)
    monkeypatch.setattr(OrcaRunPack.check_resource,'respec_memory',lambda **kwargs:'1GB')
    calls=[]
    if backend=='gaussian':
        klass,method=GaussianRunPack.GaussianDFTRun,'run_gaussian'
        monkeypatch.setattr(GaussianRunPack.gaussian_run.chk2fchk,'Get_chklist',lambda *args:None)
        def fake(name,*args):
            calls.append(Path(name+'.com').read_text())
            return complete(backend,name)
        monkeypatch.setattr(Exe_Gaussian,'_exe_Gaussian',fake)
        options=('energy','energy freq','energy freq')
    else:
        klass,method=GamessRunPack.GamessDFTRun,'run_gamess'
        def fake(name,*args):
            calls.append(Path(name+'.inp').read_text())
            return complete(backend,name)
        monkeypatch.setattr(Exe_Gamess,'_exe_Gamess',fake)
        options=('energy','energy uv','energy uv dipole','energy uv dipole')
    monkeypatch.setattr(klass,'Extract_values',lambda *args:{'energy':[-1,0]})
    for option in options:
        job=klass('B3LYP','STO-3G',1,option,'mol.xyz')
        result=getattr(job,method)()
        assert result['log']=='normal',result
    assert len(calls)==(2 if backend=='gaussian' else 3)
    assert result['reuse_report']['reused']==(2 if backend=='gaussian' else 3)


def test_orca_new_property_and_failed_stage_do_not_repeat_optimization(tmp_path,monkeypatch):
    monkeypatch.chdir(tmp_path)
    Path('water.xyz').write_text('3\n0 1\nO 0 0 0\nH 0 0 .96\nH .92 0 -.24\n')
    monkeypatch.setattr(OrcaRunPack.check_resource,'respec_cores',lambda n:n)
    calls=[]
    failed=[False]
    def fake(name,*args):
        calls.append(name)
        if name.endswith('_TD') and not failed[0]:
            failed[0]=True
            Path(name+'.log').write_text('interrupted')
            raise RuntimeError('interrupted')
        text=ORCA_LOG.read_text()
        if name.endswith('_TD'): text+='\n'+TD_LOG.read_text()
        Path(name+'.log').write_text(text)
        return 'normal'
    monkeypatch.setattr(Exe_Orca,'exe_Orca',fake)
    for option in ('opt energy','opt energy uv','opt energy uv polar','opt energy uv polar'):
        job=OrcaRunPack.OrcaDFTRun('B3LYP','STO-3G',1,option,'water.xyz')
        job.mem='1GB'
        result=job.run_orca()
    assert result['log']=='normal',result
    assert calls.count('water')==1
    assert calls.count('water_TD')==2
    assert calls.count('water_polar')==1
    assert len(calls)==4
    assert result['reuse_report']['reused']==3


def test_active_lock_blocks_concurrent_run_without_touching_outputs(tmp_path):
    folder=tmp_path/'job'
    first=cache.CalculationCache('orca'); first.prepare(folder)
    (folder/'important.log').write_text('active output')
    second=cache.CalculationCache('orca')
    try:
        with pytest.raises(RuntimeError,match='lock exists'):
            second.prepare(folder)
        assert (folder/'important.log').read_text()=='active output'
    finally:
        first.finish()


def test_completed_conformer_search_is_reused_but_settings_changes_are_not(tmp_path,monkeypatch):
    from qcforever.conformer_search.settings import SearchConfig
    monkeypatch.chdir(tmp_path)
    Path('input.sdf').write_text('source structure')
    source=tmp_path/'input.sdf'
    calls=[]
    @cache.conformer_workflow
    def search(*args):
        calls.append(1)
        Path('optimized_structures.sdf').write_text('selected structure')
        Path('conformer_search').mkdir()
        Path('conformer_search/summary.json').write_text('{"state":"succeeded"}')
        Path('conformer_search/trajectory.xyz').write_text('trajectory')
        return {'state':'succeeded'}
    for config in (SearchConfig(),SearchConfig(),SearchConfig(seed=123)):
        with session(tmp_path/'job','orca'):
            result=search(source,0,1,'xtb',1,'1GB',config)
            assert result['state']=='succeeded'
            assert Path('conformer_search/trajectory.xyz').read_text()=='trajectory'
    assert len(calls)==2


def test_new_executable_forces_recomputation(tmp_path,monkeypatch):
    runtime=tmp_path/'orca'
    runtime.write_text('first program')
    runtime.chmod(0o755)
    calls=[]
    def run():
        calls.append(1)
        return complete('orca')
    for version in ('first program','a changed program'):
        runtime.write_text(version)
        with session(tmp_path/'job','orca'):
            cache.execute('orca','mol',example('orca'),run,str(runtime))
    assert len(calls)==2
