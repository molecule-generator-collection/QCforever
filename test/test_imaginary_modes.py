import numpy as np
import pytest

from qcforever.gamess_run.parse_log import parse_log as GamessLog
from qcforever.gaussian_run.Get_FreqPro import Extract_vibvec
from qcforever.util.imaginary_modes import StationaryPoint, displace_imaginary_mode, follow_imaginary_modes
from qcforever.util.job_timeout import QCforeverTimeoutError


def test_follow_modes_checks_both_signs_and_selects_a_minimum():
    coordinates = np.zeros((2, 3))
    mode = [[[1, 2]], [[0, 2]], [[0, 1]]]
    start = StationaryPoint('ground', ['H', 'H'], coordinates, -1, [-10], mode, {})
    attempts = {}

    def evaluate(label, atoms, xyz):
        attempts[label] = xyz
        if label.endswith('_plus'):
            return StationaryPoint(label, atoms, xyz, -3, [-5], mode, {})
        return StationaryPoint(label, atoms, xyz, -2, [5], mode, {})

    point = follow_imaginary_modes(start, evaluate)
    plus, minus = attempts['1_plus'], attempts['1_minus']
    np.testing.assert_allclose((plus + minus) / 2, coordinates)
    assert np.linalg.norm(plus - coordinates, axis=1).max() == pytest.approx(0.05)
    assert point.name == '1_minus' and point.minimum_frequency == 5


def test_follow_modes_can_recover_from_one_failed_branch():
    start = StationaryPoint('ground', ['H'], np.zeros((1, 3)), -1, [-5], [[[1]], [[0]], [[0]]], {})
    calls = []

    def evaluate(label, atoms, xyz):
        calls.append(label)
        if label.endswith('_plus'):
            raise RuntimeError('SCF failed')
        return StationaryPoint(label, atoms, xyz, -2, [10], [], {})

    assert follow_imaginary_modes(start, evaluate).name == '1_minus'
    assert calls == ['1_plus', '1_minus']

    def timeout(*args):
        raise QCforeverTimeoutError('deadline')

    with pytest.raises(QCforeverTimeoutError):
        follow_imaginary_modes(start, timeout)


def test_gaussian_single_mode_block():
    lines = [' Atom AN X Y Z', ' 1 8 0.0 0.1 0.2', ' 2 1 0.3 0.4 0.5', '']
    np.testing.assert_array_equal(Extract_vibvec(lines), [[[0.0, 0.3]], [[0.1, 0.4]], [[0.2, 0.5]]])


def test_gamess_imaginary_mode_precedes_rigid_modes(tmp_path):
    path = tmp_path / 'modes.log'
    block = [' MODES 2 TO 6 ARE TAKEN AS ROTATIONS AND TRANSLATIONS.\n']
    for frequencies in ('52.4 I 0.0 0.0 0.0 0.0', '0.0'):
        width = len(GamessLog._frequencies(frequencies))
        block.append(f' FREQUENCY: {frequencies}\n')
        for atom in range(1, 3):
            for i, axis in enumerate('XYZ', 1):
                numbers = ' '.join(str(float(atom * i * j)) for j in range(1, width + 1))
                block.append(f' {atom} H {axis} {numbers}\n' if axis == 'X' else f' {axis} {numbers}\n')
    path.write_text(' ATOM   ATOMIC   COORDINATES (BOHR)\n'
                    '        CHARGE   X   Y   Z\n'
                    ' H  1.0  0.0  0.0  0.0\n'
                    ' H  1.0  0.0  0.0  1.4\n' + ''.join(block))
    parser = GamessLog(path)
    atoms, coordinates = parser.hessian_coordinates(2)
    assert atoms == ['H', 'H'] and coordinates[1, 2] == pytest.approx(0.740848)
    frequencies, modes = parser.normal_modes(parser.lines, 2)
    assert frequencies == [-52.4] and modes == [[[1., 2.]], [[2., 4.]], [[3., 6.]]]
    assert parser.getFreq(parser.lines)[0] == frequencies
    assert displace_imaginary_mode(np.zeros((2, 3)), frequencies, modes) is not None


def test_gaussian_seed_pickle_uses_its_cartesian_frame_and_validates_input(tmp_path):
    import pickle
    from qcforever.gaussian_run.frequency_seed import FREQUENCY_KEYS, load

    atoms = ['H', 'F']
    start = [[0., 0., 0.], [1., 0., 0.]]
    optimized = [[0., 0., 0.], [0., 0., 0.95]]
    path = tmp_path / 'saved.pkl'
    context = dict(backend='gaussian', route='#uhf/sto-3g test Freq=(Raman)',
                   atoms=atoms, coordinates_angstrom=optimized, charge=0, multiplicity=1,
                   energy=[-98., 0.])
    result = dict.fromkeys(FREQUENCY_KEYS, 0.)
    result.update(frequency_context=context, freq=[1500.], freqmode=[[[0., 0.]], [[0., 0.]], [[1., -1.]]],
                  stable_input_geometry=start, stable_0=True)
    path.write_bytes(pickle.dumps(result))
    options = dict(atoms=atoms, charge=0, mult=1, method='uhf/sto-3g',
                   solvent_route='', solvent_read='', parameter_route='')
    with pytest.raises(ValueError, match='geometry differs'):
        load(path, 'molecule', coordinates=start, **options)
    saved = load(path, 'molecule', coordinates=start, allow_completed_stable=True, **options)
    np.testing.assert_array_equal(saved.coordinates, optimized)
    assert saved.minimum_frequency == 1500.
    with pytest.raises(ValueError, match='solvent differs'):
        load(path, 'molecule', coordinates=optimized, **dict(options, solvent_route='SCRF=(PCM,solvent=Water)'))
