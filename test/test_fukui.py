"""Exercise Fukui workflows using synthetic Gaussian and GAMESS output."""

from pathlib import Path

import pytest

from qcforever.gamess_run import Exe_Gamess, GamessRunPack
from qcforever.gaussian_run import Exe_Gaussian, GaussianRunPack, chk2fchk
from qcforever.util import check_resource


ATOMS = ('O', 'H', 'H')
STATES = (
    (0, 1, -75.0, (-0.4, 0.2, 0.2)),
    (1, 2, -74.6, (0.3, 0.35, 0.35)),
    (-1, 2, -75.1, (-0.6, -0.2, -0.2)),
)


@pytest.fixture
def water(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(check_resource, 'respec_cores', lambda n: n)
    monkeypatch.setattr(check_resource, 'respec_memory', lambda **kwargs: '1GB')
    path = tmp_path / 'water.xyz'
    path.write_text('3\n0 1\nO 0 0 0\nH 0 0 .96\nH .92 0 -.24\n')
    return path


def check_fukui(result):
    assert result['log'] == 'normal', result
    assert result['vip'][0] == pytest.approx(0.4 * 27.211)
    assert result['vea'][0] == pytest.approx(0.1 * 27.211)
    assert result['cden'][1] == [-0.4, 0.2, 0.2]
    assert result['fukui']['atoms'] == list(ATOMS)
    assert result['fukui']['f_plus'] == pytest.approx([0.2, 0.4, 0.4])
    assert result['fukui']['f_minus'] == pytest.approx([0.7, 0.15, 0.15])


def gaussian_step(charge, mult, energy, charges):
    spin = ' and spin densities' if mult == 2 else ''
    rows = ''.join(f' {i} {atom} {q}{" 0.0" if spin else ""}\n'
                   for i, (atom, q) in enumerate(zip(ATOMS, charges), start=1))
    return (
        ' Initial command:\n #ub3lyp/sto-3g test\n ----\n'
        f' Charge = {charge} Multiplicity = {mult}\n'
        ' Standard orientation:\n'
        ' ----------------------------------------\n'
        ' Center Atomic Atomic Coordinates (Angstroms)\n'
        ' Number Number Type X Y Z\n'
        ' ----------------------------------------\n'
        ' 1 8 0 0.0 0.0 0.0\n 2 1 0 0.0 0.0 0.96\n 3 1 0 0.92 0.0 -0.24\n'
        ' ----------------------------------------\n'
        f' SCF Done:  E(UB3LYP) = {energy} A.U.\n'
        f' Mulliken charges{spin}:\n  1\n{rows}'
        f' Sum of Mulliken charges = {charge}\n'
        ' Normal termination of Gaussian 16\n'
    )


def test_gaussian_fukui_enables_vertical_links(water, monkeypatch):
    def launch(name, *_):
        steps = Path(name + '.com').read_text().split('--Link1--')
        assert len(steps) == 3
        for step, suffix, charge in zip(steps[1:], ('IP', 'EA'), ('1', '-1')):
            assert f'%chk={name}_{suffix}\n' in step
            assert f'%Oldchk={name}\n' in step
            assert [charge, '2'] in [line.split() for line in step.splitlines()]
        Path(name + '.log').write_text(''.join(gaussian_step(*state) for state in STATES))
        return 'normal'

    monkeypatch.setattr(Exe_Gaussian, 'exe_Gaussian', launch)
    monkeypatch.setattr(chk2fchk, 'Get_chklist', lambda *_: None)
    job = GaussianRunPack.GaussianDFTRun('B3LYP', 'STO-3G', 1, 'cden=(fukui)', water)
    check_fukui(job.run_gaussian())


def test_gamess_fukui_enables_vertical_jobs(water, monkeypatch):
    calls = []

    def launch(name, *_):
        charge, mult, energy, charges = STATES[len(calls)]
        calls.append(name)
        input_text = Path(name + '.inp').read_text()
        assert 'RUNTYP=ENERGY' in input_text
        assert f'ICHARG={charge} MULT={mult}' in input_text
        rows = ''.join(f' {i} {atom} 0.0 {q} 0.0 {q}\n'
                       for i, (atom, q) in enumerate(zip(ATOMS, charges), start=1))
        Path(name + '.log').write_text(
            f' FINAL UHF ENERGY IS {energy} AFTER 10 ITERATIONS\n'
            ' TOTAL MULLIKEN AND LOWDIN ATOMIC POPULATIONS\n'
            '           MULLIKEN        LOWDIN\n'
            f'{rows}\n EXECUTION OF GAMESS TERMINATED NORMALLY\n'
        )
        Path(name + '.dat').write_text(
            ' $DATA\nwater\nC1\nO 8 0 0 0\nH 1 0 0 .96\nH 1 .92 0 -.24\n $END\n'
            ' $VEC\n 1 1 0.0\n 2 1 0.0\n $END\n'
        )
        return ''

    monkeypatch.setattr(Exe_Gamess, 'exe_Gamess', launch)
    job = GamessRunPack.GamessDFTRun('B3LYP', 'STO-3G', 1, 'cden = ( Fukui )', water)
    check_fukui(job.run_gamess())
    assert calls == ['water', 'water_VIP', 'water_VEA']
