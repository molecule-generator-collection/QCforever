"""Geometry/state handoff from Gaussian text checkpoints to other engines."""

from pathlib import Path

import numpy as np
import pytest

from qcforever.util import read_mol_file, check_resource
from qcforever.orca_run import Exe_Orca, OrcaRunPack
from qcforever.gamess_run import Exe_Gamess, GamessRunPack


BOHR = 0.529177210903
WATER_LOG = Path(__file__).parent / 'fixtures/orca/water.txt'


def fchk_text(atoms=(8, 1, 1), charge=0, multiplicity=1, nuclear=None):
    values = [0.0, 0.0, 0.0, 0.0, 0.0, 0.96 / BOHR, 0.92 / BOHR, 0.0, -0.24 / BOHR]
    if len(atoms) != 3:
        values = [0.0, 0.0, 0.0] * len(atoms)
    def field(name, kind, value):
        return f'{name:40s}   {kind}   {value}\n'
    text = 'Molecular structure from Gaussian\nSP RHF 6-31G(d)\n'
    text += field('Number of atoms', 'I', len(atoms))
    text += field('Charge', 'I', charge) + field('Multiplicity', 'I', multiplicity)
    text += field('Atomic numbers', 'I', f'N= {len(atoms)}')
    text += ' '.join(str(value) for value in atoms) + '\n'
    if nuclear is not None:
        text += field('Nuclear charges', 'R', f'N= {len(nuclear)}')
        text += ' '.join(f'{value:.8E}' for value in nuclear) + '\n'
    text += field('Current cartesian coordinates', 'R', f'N= {len(values)}')
    for i in range(0, len(values), 5):
        # Fixed-width Fortran fields can put a minus sign immediately after
        # the preceding exponent. Also exercise Gaussian's D exponents.
        text += ''.join(f'{value:16.8E}'.replace('E', 'D') for value in values[i:i+5]) + '\n'
    text += field('SCF Energy', 'R', '-999.0')
    text += field('Alpha MO coefficients', 'R', 'N= 999999999')
    return text


def test_reader_converts_bohr_and_d_exponents_without_reading_orbitals(tmp_path):
    path = tmp_path / 'water.fchk'
    original = fchk_text()
    path.write_text(original)
    atoms, x, y, z, charge, mult = read_mol_file.read_fchk(path)
    assert atoms == ['O', 'H', 'H']
    assert np.column_stack((x, y, z)) == pytest.approx(np.array([[0, 0, 0], [0, 0, 0.96], [0.92, 0, -0.24]]))
    assert (charge, mult) == (0, 1)
    assert path.read_text() == original


def test_ecp_nuclear_charges_do_not_change_element_identity(tmp_path):
    path = tmp_path / 'hbr.fchk'
    path.write_text(fchk_text(atoms=(35, 1), nuclear=(7, 1)))
    assert read_mol_file.read_fchk(path)[0] == ['Br', 'H']


@pytest.mark.parametrize('problem', ['missing_atomic_numbers', 'missing_charge', 'short_array',
    'coordinate_count', 'invalid_atomic_number', 'ghost', 'zero_multiplicity',
    'nonfinite', 'wrong_type', 'ghost_nuclear', 'nuclear_count'])
def test_invalid_structure_is_rejected(tmp_path, problem):
    text = fchk_text()
    if problem == 'missing_atomic_numbers':
        text = text.replace('Atomic numbers', 'Unused numbers')
    elif problem == 'missing_charge':
        text = text.replace('Charge', 'Unused')
    elif problem == 'short_array':
        text = text[:text.index('Current cartesian coordinates')]
        text += 'Current cartesian coordinates              R   N= 9\n0 0 0\n'
    elif problem == 'coordinate_count':
        text = text.replace('R   N= 9', 'R   N= 10')
    elif problem == 'invalid_atomic_number':
        text = fchk_text(atoms=(119, 1, 1))
    elif problem == 'ghost':
        text = fchk_text(atoms=(0, 1, 1))
    elif problem == 'zero_multiplicity':
        text = fchk_text(multiplicity=0)
    elif problem == 'nonfinite':
        text = text.replace('0.00000000D+00', '1.00000000D+999', 1)
    elif problem == 'wrong_type':
        text = text.replace('Atomic numbers                             I', 'Atomic numbers                             R')
    elif problem == 'ghost_nuclear':
        text = fchk_text(nuclear=(0, 1, 1))
    elif problem == 'nuclear_count':
        text = fchk_text(nuclear=(8, 1))
    path = tmp_path / 'invalid.fchk'
    path.write_text(text)
    with pytest.raises(ValueError, match='fchk'):
        read_mol_file.read_fchk(path)


@pytest.fixture(params=['gamess', 'orca'])
def engine(request, monkeypatch):
    monkeypatch.setattr(check_resource, 'respec_cores', lambda n: n)
    calls = []
    if request.param == 'orca':
        klass, method = OrcaRunPack.OrcaDFTRun, 'run_orca'
        def launch(name, *args):
            calls.append(Path(name + '.inp').read_text())
            Path(name + '.log').write_text(WATER_LOG.read_text())
            return 'normal'
        monkeypatch.setattr(Exe_Orca, 'exe_Orca', launch)
    else:
        klass, method = GamessRunPack.GamessDFTRun, 'run_gamess'
        def launch(name, *args):
            calls.append(Path(name + '.inp').read_text())
            Path(name + '.log').write_text('EXECUTION OF GAMESS TERMINATED NORMALLY\n')
            Path(name + '.dat').write_text('restart data')
            return ''
        monkeypatch.setattr(Exe_Gamess, '_exe_Gamess', launch)
        monkeypatch.setattr(klass, 'Extract_values', lambda *args: {'energy': [-75, 0]})
        monkeypatch.setattr(GamessRunPack.job_cleanup, 'cleanup_gamess', lambda *args, **kwargs: None)
    return request.param, klass, method, calls


@pytest.mark.parametrize('suffix', ['.fchk', '.FCHK', '.fch'])
@pytest.mark.parametrize('override', [False, True])
def test_fchk_workflow_preserves_state_geometry_and_destination_level(tmp_path, monkeypatch, engine, suffix, override):
    backend, klass, method, calls = engine
    inputs, runs = tmp_path / 'inputs', tmp_path / 'runs'
    inputs.mkdir(); runs.mkdir()
    path = inputs / ('water' + suffix)
    original = fchk_text(charge=1, multiplicity=2)
    path.write_text(original)
    monkeypatch.chdir(runs)
    job = klass('B3LYP', 'STO-3G', 1, 'energy', path)
    job.mem = '1GB'
    if override:
        job.SpecTotalCharge, job.SpecSpinMulti = 0, 1
    result = getattr(job, method)()
    assert result['log'] == 'normal', result
    text = calls[0]
    if backend == 'orca':
        assert f'* xyz {0 if override else 1} {1 if override else 2}' in text
        assert 'LibXC(B3LYP/G) STO-3G' in text
    else:
        assert f'ICHARG={0 if override else 1} MULT={1 if override else 2}' in text
        assert 'DFTTYP=B3LYP' in text
        assert 'GBASIS=STO NGAUSS=3' in text
        assert 'GUESS=HUCKEL' in text
    assert '0.960000' in text
    assert path.read_text() == original
    assert not list(runs.rglob('*.chk'))
    assert not list(inputs.glob('*.xyz'))
    assert result['energy'][0] != -999.0


@pytest.mark.parametrize('profile', ['low', 'legacy'])
def test_fchk_optconf_uses_xyz_handoff_and_resolved_electronic_state(tmp_path, monkeypatch, engine, profile):
    from qcforever.conformer_search import conformer_search
    backend, klass, method, calls = engine
    path = tmp_path / 'water.fchk'
    path.write_text(fchk_text())
    monkeypatch.chdir(tmp_path)
    seen = []
    def search(filename, charge, mult, relax, nproc, memory, config):
        assert str(filename).endswith('.xyz')
        values = read_mol_file.read_xyz(filename)
        assert values[0] == ['O', 'H', 'H']
        assert values[4:] == (1, 2)
        seen.append((charge, mult, relax, config.profile))
        return {'state': 'succeeded'}
    monkeypatch.setattr(conformer_search, 'run_conformer_search', search)
    monkeypatch.setattr(conformer_search, 'read_selected_structure', lambda _: (
        ['O', 'H', 'H'], [0, 0, .92], [0, 0, 0], [0, .96, -.24], 0, 1, [], []))
    job = klass('B3LYP', 'STO-3G', 1,
                f'optconf=(conf_gen={profile}, conf_relax=xtb) energy', path)
    job.mem = '1GB'
    job.SpecTotalCharge, job.SpecSpinMulti = 1, 2
    result = getattr(job, method)()
    assert result['log'] == 'normal', result
    assert result['optconf'] is True
    assert seen == [(1, 2, 'xtb', profile)]
    assert ('* xyz 1 2' if backend == 'orca' else 'ICHARG=1 MULT=2') in calls[0]


def test_source_gaussian_energy_and_method_do_not_replace_destination_result_or_force_recompute(tmp_path, monkeypatch, engine):
    backend, klass, method, calls = engine
    path = tmp_path / 'water.fchk'
    original = fchk_text()
    path.write_text(original)
    monkeypatch.chdir(tmp_path)
    for text in (original, original.replace('RHF 6-31G(d)', 'MP2 cc-pVTZ').replace('-999.0', '-12345.0')):
        path.write_text(text)
        job = klass('B3LYP', 'STO-3G', 1, 'energy', path)
        job.mem = '1GB'
        result = getattr(job, method)()
        assert result['log'] == 'normal', result
        assert result['energy'][0] not in (-999.0, -12345.0)
    assert len(calls) == 1
    assert result['reuse_report']['reused'] == 1
