"""ORCA output excerpts and workflow tests; no ORCA installation is required."""

import os
from pathlib import Path
import pickle
import shutil
import sys

import numpy as np
import pytest

from qcforever.orca_run import Exe_Orca, OrcaRunPack
from qcforever.orca_run.parse_log import parse_log
from qcforever.util import job_cleanup, job_timeout


FIXTURES = Path(__file__).parent / 'fixtures' / 'orca'


@pytest.fixture
def job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(OrcaRunPack.check_resource, 'respec_cores', lambda value: value)
    source = tmp_path / 'water.xyz'
    source.write_text('3\n0 1\nO 0 0 0\nH 0 0 .96\nH .92 0 -.24\n')
    result = OrcaRunPack.OrcaDFTRun('B3LYP', 'STO-3G', 1, 'energy', source)
    result.mem = '1GB'
    return result


@pytest.fixture
def launcher(monkeypatch):
    calls = []

    def launch(name, *_):
        calls.append((name, Path(f'{name}.inp').read_text()))
        fixture = 'cation.txt' if name.endswith(('_VIP', '_AIP', '_VEA', '_AEA')) else 'water.txt'
        text = (FIXTURES / fixture).read_text()
        if name.endswith(('_TD', '_TDopt', '_Emission')):
            text += '\n' + (FIXTURES / 'tddft.txt').read_text()
        if name.endswith(('_AIP', '_AEA')):
            text += '\nTHE OPTIMIZATION HAS CONVERGED\n'
        Path(f'{name}.log').write_text(text)
        Path(f'{name}.gbw').write_text('checkpoint')
        Path(f'{name}.tmp').write_text('temporary')
        return 'normal'

    monkeypatch.setattr(Exe_Orca, 'exe_Orca', launch)
    return calls


def test_closed_and_open_shell_properties():
    water = parse_log(FIXTURES / 'water.txt')
    water.check(optimize=True)
    assert water.energy() == [-75.323553799694, 0]
    assert water.orbitals()['homolumo'] == pytest.approx(12.2642806266)
    assert water.dipole()[-1] == pytest.approx(1.592941334)
    assert np.linalg.norm(water.dipole()[:3]) == pytest.approx(water.dipole()[-1])
    assert water.polar()['polar_iso'] == pytest.approx(0.4313131565)
    assert water.quadrupole()['quadrupole_iso'] == pytest.approx(-5.281643535)
    ion = parse_log(FIXTURES / 'cation.txt')
    assert ion.energy()[1] == pytest.approx(0.002195)
    mo = ion.orbitals()
    assert len(mo['Alpha_MO'][0]) == 5
    assert len(mo['Beta_MO'][0]) == 4
    assert mo['homolumo'] == pytest.approx([17.2191747368, 6.6145709801])
    assert sum(ion.charges()[1]) == pytest.approx(1.0)


def test_vibrations_and_thermochemistry():
    result = parse_log(FIXTURES / 'water.txt').frequencies()
    assert len(result['freq']) == len(result['IR']) == len(result['Raman']) == 9
    assert result['IR'][-3:] == [2.73, 69.41, 53.48]
    assert result['Raman'][-3:] == pytest.approx([9.732553, 38.934224, 19.478787])
    assert np.shape(result['freqmode']) == (3, 9, 3)
    assert result['Ezp'] == pytest.approx(-75.30200196)
    assert result['E_free'] == pytest.approx(-75.31985227)
    assert result['Si'] == pytest.approx(45.5191352552)


def test_orca6_absorption_table():
    rows = parse_log(FIXTURES / 'tddft.txt').transitions()
    assert rows == [(1, 121.7, 0.001524730), (2, 105.1, 0.0), (3, 88.3, 0.069604279)]


def test_final_block_is_selected_and_incomplete_quadrupole_is_unavailable(tmp_path):
    path = tmp_path / 'multiple.log'
    path.write_text((FIXTURES / 'water.txt').read_text() + '\nFINAL SINGLE POINT ENERGY -80.5\n'
                    'QUADRUPOLE MOMENT\nMethod : SCF\n')
    parsed = parse_log(path)
    assert parsed.energy()[0] == -80.5
    assert parsed.quadrupole()['quadrupole'] is None
    assert parsed.quadrupole()['quadrupole_info']['status'] == 'unavailable'


@pytest.mark.parametrize('memory,cores,expected', [('4GB', 4, 750), ('2048MB', 2, 768), ('1GiB', 1, 805)])
def test_total_memory_is_converted_to_per_process_budget(memory, cores, expected):
    assert OrcaRunPack.memory_per_core(memory, cores) == expected


@pytest.mark.parametrize('value', ['0GB', '-4GB', 'NaN', '2', '2GB\n! HF'])
def test_invalid_memory_is_rejected(value):
    with pytest.raises(ValueError):
        OrcaRunPack.memory_per_core(value, 1)


@pytest.mark.parametrize('options', ['KTLC-BLYP-BO', 'optconf=pm6',
    'optconf=(conf_relax=pm6)', 'stable', 'opt=bad', 'energy=anything', 'uv=', 'nac'])
def test_unsupported_requests_do_not_launch(job, launcher, options):
    if options.startswith('KTLC'):
        job.functional = options
    else:
        job.value = options
    result = job.run_orca()
    assert result['log'] == 'error'
    assert launcher == []
    assert not Path('water').exists()


@pytest.mark.parametrize('option', ['optconf energy', 'optconf=xtb energy',
                                  'optconf=(conf_gen=low, laqa=off) energy'])
def test_conformer_default_is_xtb(job, option):
    job.value = option
    options = job._options()
    assert options['optconf'] == 'xtb'
    assert job._search_config.profile == 'low'


def test_complete_workflow_uses_optimized_geometry_and_keeps_electronic_state(job, launcher):
    job.value = 'opt energy homolumo dipole cden polar vip vea aip aea uv fluor=1 freq nmr'
    job.nmr_reference = {'H': 31.9}
    job.pklsave = True
    result = job.run_orca()
    assert result['log'] == 'normal', result
    assert result['vip'][0] == pytest.approx(9.5970772297)
    assert result['vea'][0] == pytest.approx(-result['vip'][0])
    assert result['aip'][0] == pytest.approx(result['vip'][0])
    assert result['aea'][0] == pytest.approx(result['vea'][0])
    assert result['fluor_state']['root'] == 1
    assert result['nmr'][1][1:] == pytest.approx([-0.199, -0.2])
    assert result['nmr_shielding'][1][0] == pytest.approx(289.983)
    assert result['charge'] == 0 and result['multiplicity'] == 1
    ion_inputs = dict(launcher)
    assert '* xyz 1 2' in ion_inputs['water_AIP']
    assert '* xyz -1 2' in ion_inputs['water_AEA']
    assert '-0.052317' in ion_inputs['water_VIP']
    assert 'NumFreq' in ion_inputs['water_freq']
    assert 'LibXC(B3LYP/G)' in ion_inputs['water']
    assert ion_inputs['water_nmr'].index('* xyz') < ion_inputs['water_nmr'].index('%eprnmr')
    assert Path('water/water.gbw').exists()
    assert Path('water/water.xyz').exists()
    assert not Path('water/water.tmp').exists()
    with Path('water/water.pkl').open('rb') as stream:
        assert pickle.load(stream) == result


def test_conformer_search_uses_xtb_and_preserves_charge_and_spin(job, launcher, monkeypatch):
    from qcforever.conformer_search import conformer_search
    seen = []

    def search(filename, charge, mult, method, *_):
        seen.append((filename, charge, mult, method))
        return {'state': 'succeeded'}

    monkeypatch.setattr(conformer_search, 'run_conformer_search', search)
    monkeypatch.setattr(conformer_search, 'read_selected_structure', lambda _: (
        ['O', 'H', 'H'], [0, 0, 0.92], [0, 0, 0], [0, 0.96, -0.24], 0, 1, [], []))
    job.value = 'optconf energy'
    job.SpecTotalCharge = 1
    job.SpecSpinMulti = 2
    result = job.run_orca()
    assert result['log'] == 'normal'
    assert seen[0][1:] == (1, 2, 'xtb')
    assert '* xyz 1 2' in launcher[0][1]
    assert result['optconf'] is True


def test_search_failure_is_reported_without_silent_fallback(job, launcher, monkeypatch):
    from qcforever.conformer_search import conformer_search

    def fail(*_):
        raise RuntimeError('xTB not installed')

    monkeypatch.setattr(conformer_search, 'run_conformer_search', fail)
    job.value = 'optconf energy'
    result = job.run_orca()
    assert result['log'] == 'error'
    assert result['optconf'] is False
    assert result['conformer_search']['state'] == 'failed'
    assert launcher == []


def test_existing_work_is_preserved(job, launcher):
    Path('water').mkdir()
    marker = Path('water/important.txt')
    marker.write_text('keep')
    result = job.run_orca()
    assert result['log'] == 'normal'
    history = list(Path('water/.qcforever/history').glob('*/important.txt'))
    assert len(history) == 1 and history[0].read_text() == 'keep'
    assert len(launcher) == 1


def test_absolute_input_and_current_directory_are_preserved(job, launcher, tmp_path):
    elsewhere = tmp_path / 'elsewhere'
    elsewhere.mkdir()
    os.chdir(elsewhere)
    assert job.run_orca()['log'] == 'normal'
    assert Path.cwd() == elsewhere
    assert (tmp_path / 'water.xyz').is_file()


def test_restart_false_drops_gbw_but_retains_geometry(job, launcher):
    job.restart = False
    assert job.run_orca()['log'] == 'normal'
    assert Path('water/water.xyz').exists()
    assert not Path('water/water.gbw').exists()


def test_inconsistent_electronic_state_is_rejected(job, launcher):
    job.SpecSpinMulti = 2
    result = job.run_orca()
    assert result['log'] == 'error'
    assert 'electron count' in result['error']
    assert launcher == []


def test_reference_spectrum_paths_are_resolved_before_chdir(job, launcher):
    Path('ir.dat').write_text('2016.09 2.73\n3609.91 69.41\n3834.45 53.48\n')
    Path('raman.dat').write_text('2016.09 9.732553\n3609.91 38.934224\n3834.45 19.478787\n')
    Path('nmr.dat').write_text('-.199 1\n-.200 1\n')
    Path('uv.dat').write_text('121.7 .001524730\n105.1 0\n88.3 .069604279\n')
    job.value = 'uv=uv.dat freq=ir.dat,raman.dat,nmr.dat'
    job.nmr_reference = {'H': 31.9}
    result = job.run_orca()
    assert result['log'] == 'normal', result
    assert result['Similarity/Dissimilarity'][0] == pytest.approx(1)
    assert result['Wasserstein_distance'][0] == pytest.approx(0)
    for kind in ('IR', 'Raman', 'NMR'):
        assert result['Spectrum_similarity'][kind]['similarity'] == pytest.approx(1)


def test_automatic_tms_reference_uses_same_calculation_settings(job, launcher, monkeypatch):
    original_launch = Exe_Orca.exe_Orca

    def launch(name, *args):
        original_launch(name, *args)
        if name.endswith('_TMS'):
            Path(f'{name}.log').write_text((FIXTURES / 'water.txt').read_text() +
                '\nCHEMICAL SHIELDING SUMMARY (ppm)\n'
                '0 Si 330.0 1.0\n1 C 185.0 1.0\n2 H 31.9 1.0\n')

    monkeypatch.setattr(Exe_Orca, 'exe_Orca', launch)
    job.value = 'nmr'
    result = job.run_orca()
    assert result['log'] == 'normal', result
    assert result['nmr_reference']['H'] == 31.9
    tms_input = dict(launcher)['water_TMS']
    assert 'LibXC(B3LYP/G) STO-3G' in tms_input
    assert '* xyz 0 1' in tms_input
    assert 'TightOpt' in tms_input


def test_normal_exit_without_energy_is_not_success(job, tmp_path):
    job.executable = executable(tmp_path, 'print("****ORCA TERMINATED NORMALLY****")\n')
    result = job.run_orca()
    assert result['log'] == 'error'
    assert 'FINAL SINGLE POINT ENERGY' in result['error']


def executable(tmp_path, body):
    program = tmp_path / 'fake-orca'
    program.write_text(f'#!{sys.executable}\n' + body)
    program.chmod(0o755)
    return str(program)


@pytest.mark.parametrize('exit_code', [0, 9])
def test_launcher_detects_nonzero_exit_and_missing_success(job, tmp_path, exit_code):
    job.executable = executable(tmp_path, f'print("SCF failed", flush=True)\nraise SystemExit({exit_code})\n')
    result = job.run_orca()
    assert result['log'] == 'error'
    assert result['returncode'] == exit_code
    assert 'SCF failed' in result['error']
    assert Path(result['logfile']).read_text().strip() == 'SCF failed'


@pytest.mark.parametrize('limit', ['timejob', 'timexe'])
def test_timeout_retains_partial_log_and_restores_directory(job, tmp_path, limit):
    job.executable = executable(tmp_path, 'import time\nprint("started", flush=True)\ntime.sleep(10)\n')
    setattr(job, limit, 0.5)
    original = Path.cwd()
    result = job.run_orca()
    assert result['log'] == 'timeout'
    assert Path.cwd() == original
    assert Path(result['logfile']).exists()
    assert Path('water/water.inp').exists()


def test_cleanup_keeps_conformer_histories(tmp_path):
    (tmp_path / 'conformer_search').mkdir()
    (tmp_path / 'conformer_search' / 'history.json').write_text('{}')
    (tmp_path / 'optimized_structures.sdf').touch()
    (tmp_path / 'run.hess').touch()
    (tmp_path / 'run.scratch').touch()
    job_cleanup.cleanup_orca(tmp_path, preserve_conformers=True)
    assert (tmp_path / 'conformer_search' / 'history.json').exists()
    assert (tmp_path / 'optimized_structures.sdf').exists()
    assert (tmp_path / 'run.hess').exists()
    assert not (tmp_path / 'run.scratch').exists()


@pytest.mark.skipif(not os.environ.get('QCFOREVER_TEST_ORCA') or not shutil.which('orca'),
                    reason='Set QCFOREVER_TEST_ORCA=1 to run installed ORCA')
def test_real_orca_energy_uv_and_ion(job):
    job.value = 'energy homolumo dipole cden polar nmr uv vip'
    job.nmr_reference = {'H': 31.9}
    job.nstates = 3
    job.timejob = 120
    result = job.run_orca()
    assert result['log'] == 'normal', result
    assert result['energy'][0] < -70
    assert 5 < result['vip'][0] < 20
    assert len(result['uv'][0]) == 3
    assert len(result['cden'][0]) == 3
    assert result['polar_iso'] > 0
    assert len(result['nmr'][0]) == 3
