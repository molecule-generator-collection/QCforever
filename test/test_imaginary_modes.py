import numpy as np
import pytest

from qcforever.gamess_run.parse_log import parse_log as GamessLog
from qcforever.gamess_run.GamessRunPack import GamessDFTRun
from qcforever.gaussian_run.Get_FreqPro import Extract_vibvec
from qcforever.util.imaginary_modes import StationaryPoint, displace_imaginary_mode, follow_imaginary_modes
from qcforever.util.job_timeout import QCforeverTimeoutError


def test_follow_modes_checks_both_signs_and_selects_a_minimum():
    coordinates = np.array([[0., 0., 0.], [1., 0., 0.]])
    mode = [[[-1, 1]], [[0, 0]], [[0, 0]]]
    start = StationaryPoint('ground', ['H', 'H'], coordinates, -1, [-10], mode, {})
    attempts = {}
    energies = {
        '1_sp_reference': -1.,
        '1_sp_1_plus': -0.9, '1_sp_1_minus': -0.95,
        '1_sp_2_plus': -1.2, '1_sp_2_minus': -1.01,
        '1_sp_3_plus': -1.1, '1_sp_3_minus': -1.02,
    }

    def single_point(label, center, xyz, reference):
        assert center is start
        assert reference == (None if label.endswith('_reference') else '1_sp_reference')
        attempts[label] = xyz
        return energies[label]

    def evaluate(label, atoms, xyz, guess):
        np.testing.assert_array_equal(xyz, attempts[guess])
        assert guess == ('1_sp_2_plus' if label.endswith('_plus') else '1_sp_3_minus')
        attempts[label] = xyz
        if label.endswith('_plus'):
            return StationaryPoint(label, atoms, xyz, -2, [5], [], {})
        return StationaryPoint(label, atoms, xyz, -3, [5], [], {})

    point = follow_imaginary_modes(start, single_point, evaluate)
    assert list(attempts) == [*energies, '1_plus', '1_minus']
    np.testing.assert_array_equal(attempts['1_sp_reference'], coordinates)
    np.testing.assert_array_equal(attempts['1_plus'], attempts['1_sp_2_plus'])
    np.testing.assert_array_equal(attempts['1_minus'], attempts['1_sp_3_minus'])
    assert np.linalg.norm(attempts['1_plus'] - coordinates, axis=1).max() == pytest.approx(0.003)
    assert np.linalg.norm(attempts['1_minus'] - coordinates, axis=1).max() == pytest.approx(0.006)
    assert point.name == '1_minus' and point.minimum_frequency == 5


def test_follow_modes_widens_when_optimization_returns_to_a_previous_geometry():
    xyz = np.array([[0., 0., 0.], [1., 0., 0.]])
    mode = [[[-1, 1]], [[0, 0]], [[0, 0]]]
    start = StationaryPoint('ground', ['H', 'H'], xyz, -1, [-10], mode, {})
    next_xyz = xyz.copy()
    next_xyz[1, 0] += 0.1
    attempts = {}

    def single_point(label, point, coords, reference):
        center = next_xyz if label.startswith('3_') else xyz
        expected_reference = '3_sp_reference' if label.startswith('3_') else '1_sp_reference'
        assert reference == (None if label.endswith('_reference') else expected_reference)
        np.testing.assert_array_equal(point.coordinates, center)
        return -1. - np.linalg.norm(coords - center)

    def evaluate(label, atoms, coords, guess):
        attempts[label] = coords
        if label.startswith('1_'):
            return StationaryPoint(label, atoms, xyz + 1., -1, [-10], mode, {})
        if label.startswith('2_'):
            return StationaryPoint(label, atoms, next_xyz, -2, [-5], mode, {})
        return StationaryPoint(label, atoms, coords, -3, [5], [], {})

    assert follow_imaginary_modes(start, single_point, evaluate).minimum_frequency == 5
    for i, center, expected in ((1, xyz, 0.006), (2, xyz, 0.012), (3, next_xyz, 0.006)):
        plus, minus = attempts[f'{i}_plus'], attempts[f'{i}_minus']
        np.testing.assert_allclose((plus + minus) / 2, center, atol=1e-12)
        assert np.linalg.norm(plus - center, axis=1).max() == pytest.approx(expected)


def test_follow_modes_caps_displacement_for_a_localized_mode():
    xyz = np.column_stack((np.arange(36), np.zeros((36, 2))))
    mode = np.zeros((3, 1, len(xyz)))
    mode[0, 0, 0] = 1.
    start = StationaryPoint('ground', ['H'] * len(xyz), xyz, -1, [-10], mode, {})
    attempts = []
    samples = []

    def single_point(label, center, coords, reference):
        distance = np.linalg.norm(coords - xyz, axis=1).max()
        samples.append(distance)
        return -1. - distance

    def evaluate(label, atoms, coords, guess):
        attempts.append(np.linalg.norm(coords - xyz, axis=1).max())
        return start

    assert follow_imaginary_modes(start, single_point, evaluate) is start
    assert samples == pytest.approx([0., 0.0075, 0.0075, 0.015, 0.015, 0.03, 0.03, 0.05, 0.05])
    assert attempts == pytest.approx([0.03, 0.03, 0.05, 0.05])


def test_follow_modes_needs_a_resolved_energy_drop_in_both_directions():
    start = StationaryPoint('ground', ['H'], np.zeros((1, 3)), -1, [-5], [[[1]], [[0]], [[0]]], {})

    def single_point(label, center, coords, reference):
        x = coords[0, 0]
        return -1. - x if x > 0 else -1. - 1e-9 * abs(x)

    def evaluate(*args):
        pytest.fail('Optimization started before confirming a downhill seed in both directions')

    assert follow_imaginary_modes(start, single_point, evaluate) is start


def test_follow_modes_can_recover_from_one_failed_branch():
    start = StationaryPoint('ground', ['H'], np.zeros((1, 3)), -1, [-5], [[[1]], [[0]], [[0]]], {})
    calls = []

    def single_point(label, center, xyz, reference):
        return -1. - abs(xyz[0, 0])

    def evaluate(label, atoms, xyz, guess):
        calls.append(label)
        if label.endswith('_plus'):
            raise RuntimeError('SCF failed')
        return StationaryPoint(label, atoms, xyz, -2, [10], [], {})

    assert follow_imaginary_modes(start, single_point, evaluate).name == '1_minus'
    assert calls == ['1_plus', '1_minus']

    def timeout(*args):
        raise QCforeverTimeoutError('deadline')

    for scan, optimize in ((single_point, timeout), (timeout, evaluate)):
        with pytest.raises(QCforeverTimeoutError):
            follow_imaginary_modes(start, scan, optimize)


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


def test_gamess_displaced_input_reads_vec_with_explicit_coordinates(tmp_path):
    source = tmp_path / 'reference.dat'
    orbitals = ' 1  1 0.12345678E+00\n 1  1 0.23456789E+00\n'
    source.write_text(' $VEC\n' + orbitals + ' $END\n')
    path = tmp_path / 'displaced.inp'
    job = GamessDFTRun('B3LYP', 'STO-3G', 1, 'stable', 'test.xyz')
    job.make_input('ENERGY', 0, 2, path, ['H'], [0.01], [0.02], [0.03], datfile=source)
    deck = path.read_text()
    assert 'GUESS=MOREAD NORB=1' in deck
    assert ' $VEC\n' + orbitals + ' $END\n' in deck
    assert 'H   1.0   0.0100000000   0.0200000000   0.0300000000' in deck


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
                  stable_input_geometry=start, stable_0=True, gaussian_log='subdir/saved.log')
    log = tmp_path / result['gaussian_log']
    log.parent.mkdir()
    log.touch()
    path.write_bytes(pickle.dumps(result))
    options = dict(atoms=atoms, charge=0, mult=1, method='uhf/sto-3g',
                   solvent_route='', solvent_read='', parameter_route='')
    with pytest.raises(ValueError, match='geometry differs'):
        load(path, 'molecule', coordinates=start, **options)
    saved = load(path, 'molecule', coordinates=start, allow_completed_stable=True, **options)
    np.testing.assert_array_equal(saved.coordinates, optimized)
    assert saved.minimum_frequency == 1500.
    assert saved.properties['stable_source'] == str(log.resolve())
    with pytest.raises(ValueError, match='solvent differs'):
        load(path, 'molecule', coordinates=optimized, **dict(options, solvent_route='SCRF=(PCM,solvent=Water)'))
