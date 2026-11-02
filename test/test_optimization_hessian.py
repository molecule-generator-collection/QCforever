"""Portable Hessian controls and verified ReadFC dependencies."""

from contextlib import contextmanager
from pathlib import Path
import os

import pytest

from qcforever.orca_run.OrcaRunPack import OrcaDFTRun
from qcforever.gamess_run.GamessRunPack import GamessDFTRun
from qcforever.util import calculation_cache as cache, check_resource
from qcforever.gamess_run import Exe_Gamess
from qcforever.util.calculation_options import ground_optimization


@pytest.fixture
def jobs(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(check_resource, 'respec_cores', lambda value: value)
    monkeypatch.setattr(check_resource, 'respec_memory', lambda **kwargs: '1GB')
    result = {backend: cls('B3LYP', 'STO-3G', 1, 'energy', tmp_path / 'water.xyz')
              for backend, cls in [('orca', OrcaDFTRun), ('gamess', GamessDFTRun)]}
    result['orca'].mem = '1GB'
    return result


def make(job, backend, fc, *, root=None, charge=0, atom='O', distance=.96):
    atoms = [atom, 'H', 'H']
    xyz = [[0, 0, 0], [0, 0, distance], [.92, 0, -.24]]
    if backend == 'orca':
        job.make_input('water', atoms, xyz, charge, 1, optimize='normal',
                       td=root is not None, root=root, fc=fc)
    else:
        job.make_input('OPTIMIZE', charge, 1, 'water.inp', atoms,
                       *map(list, zip(*xyz)), TDDFT=root is not None,
                       target=[root or 1, 1], solvent=job.solvent, fc=fc)
    return Path('water.inp').read_text()


@pytest.mark.parametrize('backend', ['orca', 'gamess'])
@pytest.mark.parametrize('fc,interval', [('CalcFC', 0), ('CalcAll', 1), ('RecalcFC=4', 4)])
@pytest.mark.parametrize('root', [None, 2])
def test_native_hessian_frequency_and_state(jobs, backend, fc, interval, root):
    text = make(jobs[backend], backend, fc, root=root)
    if backend == 'orca':
        assert 'Calc_Hess true' in text and 'NumHess true' in text
        assert ('Recalc_Hess' not in text if not interval else f'Recalc_Hess {interval}' in text)
        assert (f'IRoot {root}' in text) == (root is not None)
    else:
        assert f'HESS=CALC IHREP={interval}' in text
        assert '$FORCE METHOD=SEMINUM NVIB=2 $END' in text
        assert ('TDDFT=EXCITE' in text) == (root is not None)
        if root:
            assert f'IROOT={root}' in text and 'SCFTYP=RHF' in text


@contextmanager
def session(directory, backend):
    previous = Path.cwd()
    value = cache.CalculationCache(backend)
    token = cache._active.set(value)
    try:
        value.prepare(directory)
        os.chdir(directory)
        yield value
    finally:
        os.chdir(previous)
        value.finish()
        cache._active.reset(token)


def calculate(backend, text):
    def launch():
        if backend == 'orca':
            Path('water.log').write_text('ORCA TERMINATED NORMALLY\nTHE OPTIMIZATION HAS CONVERGED\n'
                                        'FINAL SINGLE POINT ENERGY -75.0\n' +
                                        (Path(__file__).parent / 'fixtures/orca/tddft.txt').read_text())
            Path('water.hess').write_text('$hessian\n9\n0 1 2\n0 .1 .0 .0\n')
        else:
            Path('water.log').write_text('EXECUTION OF GAMESS TERMINATED NORMALLY\n'
                                        'EQUILIBRIUM GEOMETRY LOCATED\n')
            Path('water.dat').write_text(' $HESS\nENERGY IS -75.0\n 1 1  .1 .0 .0\n $END\n')
        return 'normal'
    cache.execute(backend, 'water', text, launch)


@pytest.mark.parametrize('backend', ['orca', 'gamess'])
@pytest.mark.parametrize('root', [None, 1])
def test_readfc_restores_verified_native_hessian_and_hashes_dependency(jobs, tmp_path, backend, root):
    job = jobs[backend]
    with session(tmp_path / backend, backend):
        text = make(job, backend, 'CalcFC', root=root)
        calculate(backend, text)
    with session(tmp_path / backend, backend):
        text = make(job, backend, 'ReadFC', distance=.97, root=root)
        if backend == 'orca':
            assert 'InHess Read' in text and 'Calc_Hess' not in text
            assert cache.dependencies(backend, text)['water.readfc.hess'] is not None
            before = cache.dependencies(backend, text)
            Path('water.readfc.hess').write_text('modified')
            assert cache.dependencies(backend, text) != before
        else:
            assert 'HESS=READ IHREP=0' in text
            assert '$HESS\nENERGY IS -75.0' in text and '$FORCE' not in text


@pytest.mark.parametrize('backend', ['orca', 'gamess'])
@pytest.mark.parametrize('change', ['root', 'ground_to_excited', 'charge', 'atom', 'method', 'basis', 'solvent', 'corrupt'])
def test_readfc_rejects_incompatible_or_modified_hessian(jobs, tmp_path, backend, change):
    job = jobs[backend]
    with session(tmp_path / backend, backend) as state:
        text = make(job, backend, 'CalcFC', root=1 if change == 'root' else None)
        calculate(backend, text)
        if change == 'corrupt':
            entry = next(iter(state.index['entries'].values()))
            suffix = '.hess' if backend == 'orca' else '.dat'
            (state.state / entry['snapshot'] / ('water' + suffix)).write_text('corrupted')
    with session(tmp_path / backend, backend):
        kwargs = {'root': 2} if change in ('root', 'ground_to_excited') else {'charge': 1} if change == 'charge' else {'atom': 'N'} if change == 'atom' else {}
        if change in ('method', 'basis', 'solvent'):
            setattr(job, {'method': 'functional'}.get(change, change),
                    {'method': 'PBE', 'basis': '6-31G', 'solvent': 'water'}[change])
        with pytest.raises(ValueError, match='no verified Hessian'):
            make(job, backend, 'ReadFC', **kwargs)


@pytest.mark.parametrize('backend', ['orca', 'gamess'])
def test_missing_hessian_does_not_silently_fall_back(jobs, backend):
    with pytest.raises(ValueError, match='ReadFC'):
        make(jobs[backend], backend, 'ReadFC')


def test_orca_ground_and_excited_fc_are_independent(jobs):
    job = jobs['orca']
    job.value = 'opt=(tight,FC=CalcAll) fluor=(FC=CalcFC,state=2)'
    options = job._options()
    assert options['opt'] == 'tight' and options['opt_fc'] == 'CalcAll'
    assert options['fluor'] == 2 and options['fluor_fc'] == 'CalcFC'


def test_gamess_unsupported_tddft_spin_is_rejected(jobs):
    with pytest.raises(ValueError, match='singlet'):
        jobs['gamess'].make_input('OPTIMIZE', 0, 1, 'water.inp', ['H'], [0], [0], [0],
                                  TDDFT=True, target=[1, 3], fc='CalcFC')


@pytest.mark.parametrize('option', ['(FC=none,FC=CalcFC)', '(tight,normal)', 'RCFC', 'RecalcFC=0'])
def test_invalid_ground_controls_are_rejected(option):
    with pytest.raises(ValueError):
        ground_optimization(option, 'orca')


def test_gamess_version_change_prevents_hessian_reuse(jobs, tmp_path):
    with session(tmp_path / 'gamess', 'gamess') as state:
        state.variant = '01'
        calculate('gamess', make(jobs['gamess'], 'gamess', 'CalcFC'))
    with session(tmp_path / 'gamess', 'gamess') as state:
        state.variant = '02'
        with pytest.raises(ValueError, match='no verified Hessian'):
            make(jobs['gamess'], 'gamess', 'ReadFC')


def test_gamess_workflow_keeps_fc_separate_and_reads_excited_hessian(jobs, monkeypatch):
    calls = []
    def launch(name, *args):
        text = Path(name + '.inp').read_text()
        calls.append((name, text))
        Path(name + '.log').write_text('EXECUTION OF GAMESS TERMINATED NORMALLY\n'
                                     'EQUILIBRIUM GEOMETRY LOCATED\n')
        data = (' $DATA\nwater\nC1\nO 8.0 0 0 0\nH 1.0 0 0 .96\nH 1.0 .92 0 -.24\n $END\n'
                ' $VEC\n 1 1 .1 .2\n 2 1 .2 .1\n $END\n')
        if 'HESS=CALC' in text:
            data += ' $HESS\nENERGY IS -75.0\n 1 1 .1 .0 .0\n $END\n'
        Path(name + '.dat').write_text(data)
        return ''
    monkeypatch.setattr(Exe_Gamess, '_exe_Gamess', launch)
    monkeypatch.setattr(GamessDFTRun, 'Extract_values', lambda *args: {'energy': [-75, 0]})
    Path('water.xyz').write_text('3\n0 1\nO 0 0 0\nH 0 0 .96\nH .92 0 -.24\n')
    job = jobs['gamess']
    job.value = 'opt=CalcAll energy fluor=(FC=CalcFC,state=2)'
    job.run_gamess()
    inputs = dict(calls)
    assert 'HESS=CALC IHREP=1' in inputs['water']
    assert 'HESS=' not in inputs['water_TD']
    assert 'HESS=CALC IHREP=0' in inputs['water_TDopt']
    assert 'IROOT=2' in inputs['water_TDopt']
    calls.clear()
    job.value = 'opt=CalcAll energy fluor=(FC=ReadFC,state=2)'
    result = job.run_gamess()
    assert result['reuse_report']['reused'] == 2
    assert len(calls) == 1 and calls[0][0] == 'water_TDopt'
    assert 'HESS=READ IHREP=0' in calls[0][1] and '$HESS' in calls[0][1]
