from pathlib import Path
from types import SimpleNamespace

import pytest

from qcforever.util import check_resource
from qcforever.util.calculation_options import requires_hessian


def test_molecular_basis_information_sto3g_water():
    information = check_resource.molecular_basis_information(
        ['O', 'H', 'H'], 'STO-3G', total_charge=0
    )

    assert information == {
        'n_basis': 7,
        'n_electrons': 10,
        'ecp_electrons': 0,
    }


def test_molecular_basis_information_lanl2dz_uses_ecp_electrons():
    information = check_resource.molecular_basis_information(
        ['I', 'H'], 'LANL2DZ', total_charge=0
    )

    assert information['n_electrons'] == 8
    assert information['ecp_electrons'] == 46


def test_recommend_memory_uses_safe_available_limit():
    memory = check_resource.recommend_memory(
        atoms=['O', 'H', 'H'],
        basis='STO-3G',
        options=['freq'],
        available_memory_bytes=768 * check_resource.MIB,
    )

    assert memory == '512MB'


def test_safe_available_limit_is_rounded_down(monkeypatch):
    monkeypatch.setattr(
        check_resource,
        'estimate_memory_bytes',
        lambda **kwargs: (10 * check_resource.GIB, None),
    )

    memory = check_resource.recommend_memory(
        atoms=['C'],
        basis='STO-3G',
        available_memory_bytes=3 * check_resource.GIB,
    )

    assert memory == '2GB'


def test_recommend_memory_has_fallback_for_unknown_basis():
    memory = check_resource.recommend_memory(
        atoms=['C'],
        basis='not-a-real-basis',
        available_memory_bytes=8 * check_resource.GIB,
    )

    assert memory == '1GB'


def test_respec_memory_preserves_an_explicit_user_value(monkeypatch):
    def fail_if_called(**kwargs):
        raise AssertionError('The estimator must not run for explicit memory.')

    monkeypatch.setattr(check_resource, 'recommend_memory', fail_if_called)

    memory = check_resource.respec_memory(
        spec_memory='4GB',
        atoms=['O', 'H', 'H'],
        basis='STO-3G',
        calculate_hessian=True,
    )

    assert memory == '4GB'


@pytest.fixture
def large_basis(monkeypatch):
    # Keep the estimate above the 1GB floor and whole-GB rounding boundaries.
    monkeypatch.setattr(check_resource, 'molecular_basis_information',
                        lambda *args: {'n_basis': 1000, 'n_electrons': 100, 'ecp_electrons': 0})
    monkeypatch.setattr(check_resource.psutil, 'virtual_memory',
                        lambda: SimpleNamespace(available=64 * check_resource.GIB))


@pytest.mark.parametrize('options', [[], ['fluor'], ['tadf']])
def test_optimization_hessian_matches_freq_memory_and_is_not_counted_twice(large_basis, options):
    def estimate(extra=(), hessian=False):
        return check_resource.estimate_memory_bytes(
            ['C'], 'STO-3G', options=options + list(extra), calculate_hessian=hessian)[0]
    baseline = estimate()
    frequency = estimate(['freq'])
    assert frequency == estimate(hessian=True) == estimate(['freq'], hessian=True)
    assert frequency == pytest.approx(1.5 * baseline, abs=2)


def test_automatic_hessian_memory_still_respects_available_ram(large_basis):
    memory = check_resource.respec_memory(
        '', ['C'], 'STO-3G', calculate_hessian=True,
        available_memory_bytes=2 * check_resource.GIB)
    assert memory == '1GB'


@pytest.mark.parametrize('controls,expected', [
    (['CalcFC'], True), (['CALCALL,Tight'], True), (['RecalcFC=5'], True),
    (['ReadFC', 'CalcFC'], True), (['', 'ReadFC', 'RCFC'], False),
    (['none', 'default', 'NoCalcFC'], False),
])
def test_hessian_computation_controls(controls, expected):
    assert requires_hessian(*controls) is expected


@pytest.mark.parametrize('backend', ['gaussian', 'gamess', 'orca'])
@pytest.mark.parametrize('option,hessian', [
    ('opt=CalcFC', True), ('opt=(CalcAll)', True), ('opt=RecalcFC=3', True),
    ('fluor=(FC=CalcFC,state=1)', True), ('fluor=(FC=CalcAll,state=1)', True),
    ('opt=CalcAll fluor=(FC=CalcFC,state=1) freq', True),
    ('fluor=(FC=ReadFC,state=1)', False), ('opt fluor', False),
    ('tadf=(FC=CalcAll,state=1)', True),
])
def test_backends_pass_optimization_hessian_to_automatic_memory(
        tmp_path, monkeypatch, large_basis, backend, option, hessian):
    if option.startswith('tadf') and backend != 'gaussian':
        pytest.skip('TADF is supported by Gaussian only')
    from qcforever.gaussian_run.GaussianRunPack import GaussianDFTRun
    from qcforever.gamess_run.GamessRunPack import GamessDFTRun
    from qcforever.orca_run.OrcaRunPack import OrcaDFTRun

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(check_resource, 'respec_cores', lambda value: value)
    source = Path('water.xyz')
    source.write_text('3\n0 1\nO 0 0 0\nH 0 0 .96\nH .92 0 -.24\n')
    cls = {'gaussian': GaussianDFTRun, 'gamess': GamessDFTRun, 'orca': OrcaDFTRun}[backend]
    job = cls('B3LYP', 'STO-3G', 1, option, source)
    original = check_resource.respec_memory
    observed = {}

    class MemoryInspected(Exception):
        pass

    def inspect(*args, **kwargs):
        observed['hessian'] = kwargs['calculate_hessian']
        observed['memory'] = original(*args, **kwargs)
        raise MemoryInspected

    monkeypatch.setattr(check_resource, 'respec_memory', inspect)
    with pytest.raises(MemoryInspected):
        if backend == 'orca':
            job._run_orca(job._options(), {})
        else:
            getattr(job, '_run_' + backend)()
    assert observed['hessian'] is hessian
    properties = ['fluor'] if 'fluor' in option else ['tadf'] if 'tadf' in option else []
    if 'freq' in option:
        properties.append('freq')
    expected = check_resource.recommend_memory(
        ['C'], 'STO-3G', options=properties, calculate_hessian=hessian)
    assert observed['memory'] == expected
