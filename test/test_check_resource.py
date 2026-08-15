from qcforever.util import check_resource


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
    )

    assert memory == '4GB'
