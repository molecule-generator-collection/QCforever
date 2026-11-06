"""Orbital handoff: dimensions, angular functions and engine input contracts."""

from pathlib import Path
import os
import shutil
import subprocess

import numpy as np
import pytest

from qcforever.util import fchk_orbitals, check_resource, calculation_cache
from qcforever.gamess_run import GamessRunPack, Exe_Gamess
from qcforever.orca_run import OrcaRunPack, Exe_Orca
from qcforever.gaussian_run import GaussianRunPack, fchk2chk


def checkpoint(kind=-2, unrestricted=False):
    """Carbon, with one S, one SP and one D shell, in Gaussian AO order."""
    nbf = 5 + (5 if kind == -2 else (6 if kind == 2 else 7))
    fields = {
        'Number of atoms': 1, 'Charge': 0, 'Multiplicity': 1,
        'Number of alpha electrons': 3, 'Number of beta electrons': 3,
        'Atomic numbers': [6], 'Nuclear charges': [6.],
        'Current cartesian coordinates': [0., 0., 0.],
        'Number of basis functions': nbf, 'Number of independent functions': nbf,
        'Shell types': [0, -1, kind], 'Number of primitives per shell': [2, 1, 1],
        'Shell to atom map': [1, 1, 1],
        'Primitive exponents': [2., .5, .3, .7],
        'Contraction coefficients': [.5, .7, 1., 1.],
        'P(S=P) Contraction coefficients': [0., 0., 1., 0.],
        'Alpha Orbital Energies': list(np.linspace(-1, 1, nbf)),
        'Alpha MO coefficients': list(np.eye(nbf).T.flatten()),
    }
    if unrestricted:
        fields['Beta Orbital Energies'] = list(np.linspace(-.9, 1.1, nbf))
        fields['Beta MO coefficients'] = list(np.fliplr(np.eye(nbf)).T.flatten())
    return fields


def write_fchk(path, fields):
    text = 'Native handoff test\nSP UHF Custom\n'
    for label, values in fields.items():
        array = isinstance(values, list)
        first = values[0] if array else values
        kind = 'I' if isinstance(first, int) else 'R'
        text += f'{label:40s}   {kind}   ' + (f'N= {len(values)}' if array else str(values)) + '\n'
        if array:
            for i in range(0, len(values), 5):
                text += ''.join(f'{v:16.8E}'.replace('E', 'D') if kind == 'R' else f'{v:12d}' for v in values[i:i+5]) + '\n'
    path.write_text(text)
    return text


@pytest.mark.parametrize('kind', [-2, 2])
@pytest.mark.parametrize('unrestricted', [False, True])
def test_parser_keeps_mo_major_order_sp_and_spin(tmp_path, kind, unrestricted):
    path = tmp_path/'carbon.fchk'
    original = write_fchk(path, checkpoint(kind, unrestricted))
    orbitals = fchk_orbitals.read(path)
    assert [s.kind for s in orbitals.shells] == [0, 0, 1, kind]
    assert [s.offset for s in orbitals.shells] == [0, 1, 2, 5]
    np.testing.assert_array_equal(orbitals.alpha, np.eye(orbitals.nmo))
    if unrestricted:
        np.testing.assert_array_equal(orbitals.beta, np.fliplr(np.eye(orbitals.nmo)))
    else:
        assert orbitals.beta is None
    assert path.read_text() == original


def test_5d_expansion_preserves_angular_wavefunction_and_overlap():
    rng = np.random.default_rng(8)
    x, y, z = rng.normal(size=(3, 43))
    # Normalized angular polynomials, removing the common radial factor.
    cart = np.array([x*x, y*y, z*z, np.sqrt(3)*x*y, np.sqrt(3)*x*z, np.sqrt(3)*y*z]).T
    pure = np.array([z*z-(x*x+y*y)/2, np.sqrt(3)*x*z, np.sqrt(3)*y*z,
                     np.sqrt(3)*(x*x-y*y)/2, np.sqrt(3)*x*y]).T
    random_mos = rng.normal(size=(5, 13))
    np.testing.assert_allclose(cart @ fchk_orbitals.D_SPH_TO_CART @ random_mos,
                               pure @ random_mos, atol=1e-13)
    overlap = np.eye(6)
    overlap[:3, :3] = 1/3
    np.fill_diagonal(overlap, 1)
    t = fchk_orbitals.D_SPH_TO_CART
    np.testing.assert_allclose(t.T @ overlap @ t, np.eye(5), atol=1e-15)


def test_nonsymmetric_matrices_are_read_mo_major_and_written_in_mkl_blocks(tmp_path):
    f = checkpoint(unrestricted=True)
    n = f['Number of basis functions']
    alpha = np.arange(n*n, dtype=float).reshape(n,n)/100
    beta = alpha + .5
    f['Alpha MO coefficients'] = list(alpha.T.flatten())
    f['Beta MO coefficients'] = list(beta.T.flatten())
    source = tmp_path/'carbon.fchk'
    write_fchk(source, f)
    orbitals = fchk_orbitals.read(source)
    np.testing.assert_allclose(orbitals.alpha, alpha)
    np.testing.assert_allclose(orbitals.beta, beta)
    target = tmp_path/'carbon.mkl'
    orbitals.write_mkl(target)
    text = target.read_text()
    assert '\n 5 D 1.0\n' in text  # MKL shell-size field has width two.
    for spin, expected in [('ALPHA', alpha), ('BETA', beta)]:
        rows = text.split('$COEFF_'+spin+'\n')[1].split('$END')[0].strip().splitlines()
        blocks = []
        for start in range(0, len(rows), n+2):
            assert all(label == 'a1g' for label in rows[start].split())
            blocks.append(np.array([[float(v) for v in row.split()] for row in rows[start+2:start+2+n]]))
        np.testing.assert_allclose(np.hstack(blocks), expected)


@pytest.mark.parametrize('problem', ['missing_beta_energy', 'short_mo', 'bad_nbf', 'bad_primitive',
    'missing_sp', 'nonfinite', 'ecp', 'bad_electrons', 'f_shell', 'wrong_center', 'duplicate'])
def test_invalid_or_unsupported_orbitals_fail(tmp_path, problem):
    f = checkpoint(unrestricted=True)
    if problem == 'missing_beta_energy': del f['Beta Orbital Energies']
    elif problem == 'short_mo': f['Alpha MO coefficients'].pop()
    elif problem == 'bad_nbf': f['Number of basis functions'] += 1
    elif problem == 'bad_primitive': f['Primitive exponents'][0] = -1.
    elif problem == 'missing_sp': del f['P(S=P) Contraction coefficients']
    elif problem == 'nonfinite': f['Alpha MO coefficients'][0] = 1e300
    elif problem == 'ecp': f['Nuclear charges'] = [4.]
    elif problem == 'bad_electrons': f['Number of alpha electrons'] = 2
    elif problem == 'f_shell': f = checkpoint(-3)
    elif problem == 'wrong_center': f['Shell to atom map'] = [1, 2, 1]
    path = tmp_path/'bad.fchk'
    text = write_fchk(path, f)
    if problem == 'duplicate': path.write_text(text + 'Number of basis functions                  I   10\n')
    if problem == 'nonfinite': path.write_text(text.replace('1.00000000D+300', '1.00000000D+999'))
    with pytest.raises(ValueError, match='fchk'):
        fchk_orbitals.read(path)


def vec_arrays(text, rows, nmo):
    section = text.split('$VEC')[1].split('$END')[0].strip('\n').splitlines()
    arrays, values = [], []
    for row in section:
        values += [float(row[i:i+15]) for i in range(5, len(row), 15)]
        if len(values) == rows*nmo:
            arrays.append(np.array(values).reshape(nmo, rows).T)
            values = []
    assert not values
    return arrays


@pytest.mark.parametrize('kind', [-2, 2])
def test_gamess_native_guess_and_basis_survive_followup(tmp_path, monkeypatch, kind):
    monkeypatch.setattr(check_resource, 'respec_cores', lambda n: n)
    monkeypatch.chdir(tmp_path)
    source = tmp_path/'carbon.fch'
    write_fchk(source, checkpoint(kind, True))
    job = GamessRunPack.GamessDFTRun('B3LYP', 'FCHK', 1, 'energy', source)
    job.mem = '1GB'
    job._fchk_orbitals = fchk_orbitals.for_job(job, 'gamess', {})
    job.make_input('ENERGY', 0, 1, 'carbon.inp', ['C'], [0], [0], [0])
    text = Path('carbon.inp').read_text()
    assert f'ISPHER={1 if kind == -2 else -1}' in text
    assert '$BASIS' not in text
    assert 'GUESS=MOREAD' in text and 'HUCKEL' not in text
    alpha, beta = vec_arrays(text, 11, job._fchk_orbitals.nmo)
    assert alpha.shape == beta.shape == (11, job._fchk_orbitals.nmo)
    assert alpha[0, 0] == 1 and beta[0, -1] == 1
    if kind == -2:
        np.testing.assert_allclose(alpha[5:, 5:], fchk_orbitals.D_SPH_TO_CART, atol=1e-8)
    else:
        np.testing.assert_allclose(alpha, np.eye(11))
    Path('carbon.dat').write_text(text)
    job.make_input('ENERGY', 1, 2, 'carbon_ion.inp', datfile='carbon.dat')
    followup = Path('carbon_ion.inp').read_text()
    assert f'NORB={job._fchk_orbitals.nmo}' in followup
    assert ' S 2\n' in followup and ' D 1\n' in followup
    assert 'ICHARG=1 MULT=2' in followup


def test_gamess_workflow_cache_tracks_imported_coefficients(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(check_resource, 'respec_cores', lambda n: n)
    source = tmp_path/'carbon.fchk'
    f = checkpoint()
    original = write_fchk(source, f)
    calls = []
    def launch(name, *args):
        text = Path(name+'.inp').read_text()
        calls.append(text)
        assert 'ISPHER=1' in text and 'GUESS=MOREAD' in text and '$VEC' in text
        Path(name+'.log').write_text('EXECUTION OF GAMESS TERMINATED NORMALLY\n')
        Path(name+'.dat').write_text(text)
        return ''
    monkeypatch.setattr(Exe_Gamess, '_exe_Gamess', launch)
    monkeypatch.setattr(GamessRunPack.GamessDFTRun, 'Extract_values', lambda *args: {'energy': [-37.,0.]})
    monkeypatch.setattr(GamessRunPack.job_cleanup, 'cleanup_gamess', lambda *args, **kwargs: None)
    job = GamessRunPack.GamessDFTRun('B3LYP', 'FCHK', 1, 'energy', source)
    job.mem = '1GB'
    assert job.run_gamess()['log'] == 'normal'
    assert job.run_gamess()['reuse_report']['reused'] == 1
    assert len(calls) == 1 and source.read_text() == original
    f['Alpha MO coefficients'][0] = .8
    write_fchk(source, f)
    result = job.run_gamess()
    assert result['log'] == 'normal' and result['reuse_report']['executed'] == 1
    assert len(calls) == 2
    assert not source.with_suffix('.chk').exists()


def test_orca_workflow_converts_own_mkl_and_reads_gbw(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(check_resource, 'respec_cores', lambda n: n)
    source = tmp_path/'carbon.fchk'
    original = write_fchk(source, checkpoint(unrestricted=True))
    calls = []
    monkeypatch.setattr(fchk_orbitals.shutil, 'which', lambda name: '/bin/orca_2mkl' if name == 'orca_2mkl' else None)
    def convert(args, **kwargs):
        calls.append(args)
        assert args == ['/bin/orca_2mkl', 'carbon_fchk_guess', '-gbw']
        mkl = Path('carbon_fchk_guess.mkl').read_text()
        assert '$COEFF_ALPHA' in mkl and '$COEFF_BETA' in mkl
        assert '5 D 1.0' in mkl and mkl.count('a1') == 20
        Path('carbon_fchk_guess.gbw').write_bytes(b'converted orbitals')
        return subprocess.CompletedProcess(args, 0)
    monkeypatch.setattr(fchk_orbitals.job_timeout, 'run', convert)
    def launch(name, *args):
        text = Path(name+'.inp').read_text()
        assert 'MOREAD' in text and '%moinp "carbon_fchk_guess.gbw"' in text
        assert text.count(' NewGTO') == 1 and ' FCHK ' not in text
        assert ' S 2\n' in text and ' D 1\n' in text
        Path(name+'.log').write_text((Path(__file__).parent/'fixtures/orca/water.txt').read_text())
        return 'normal'
    monkeypatch.setattr(Exe_Orca, 'exe_Orca', launch)
    job = OrcaRunPack.OrcaDFTRun('B3LYP', 'FCHK', 1, 'energy', source)
    job.mem = '1GB'
    result = job.run_orca()
    assert result['log'] == 'normal', result
    assert len(calls) == 1
    assert source.read_text() == original
    assert not source.with_suffix('.chk').exists()


def test_orca_does_not_fall_back_when_imported_guess_is_rejected(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(check_resource, 'respec_cores', lambda n: n)
    source = tmp_path/'carbon.fchk'
    write_fchk(source, checkpoint())
    job = OrcaRunPack.OrcaDFTRun('B3LYP', 'FCHK', 1, 'energy', source)
    job.mem = '1GB'
    job._fchk_orbitals = fchk_orbitals.read(source)
    Path('fchk_guess.gbw').write_bytes(b'guess')
    attempts = []
    def reject(*args):
        attempts.append(args)
        Path('carbon.log').write_text('error termination in GUESS\n')
        raise Exe_Orca.OrcaExecutionError('carbon', 1, 'carbon.log', 'error termination in GUESS')
    monkeypatch.setattr(calculation_cache, 'execute', reject)
    with pytest.raises(Exe_Orca.OrcaExecutionError):
        job._calculate('carbon', ['C'], np.zeros((1,3)), 0, 1, moread='fchk_guess.gbw', require_moread=True)
    assert len(attempts) == 1


@pytest.mark.parametrize('backend', ['gamess', 'orca'])
@pytest.mark.parametrize('problem', ['override', 'optconf', 'angular', 'input'])
def test_import_policy_rejects_incompatible_requests(tmp_path, monkeypatch, backend, problem):
    monkeypatch.setattr(check_resource, 'respec_cores', lambda n: n)
    source = tmp_path/'carbon.fchk'
    write_fchk(source, checkpoint())
    klass = GamessRunPack.GamessDFTRun if backend == 'gamess' else OrcaRunPack.OrcaDFTRun
    job = klass('B3LYP', 'FCHK', 1, 'energy', source)
    options = {}
    if problem == 'override': job.SpecTotalCharge = 1
    elif problem == 'optconf': options = {'optconf': True}
    elif problem == 'angular': job.fchk_angular = 'cartesian'
    elif problem == 'input': job.in_file = str(source.with_suffix('.xyz'))
    with pytest.raises(ValueError):
        fchk_orbitals.for_job(job, backend, options)


def test_orca_rejects_6d_with_actionable_error(tmp_path):
    source = tmp_path/'carbon.fchk'
    write_fchk(source, checkpoint(2))
    with pytest.raises(ValueError, match='ORCA requires 5D'):
        fchk_orbitals.read(source).validate_target('orca')


@pytest.mark.parametrize('suffix', ['.fchk', '.FCHK', '.fch'])
def test_gaussian_full_path_conversion_preserves_source_and_uses_guess(tmp_path, monkeypatch, suffix):
    monkeypatch.setattr(check_resource, 'respec_cores', lambda n: n)
    source = tmp_path/'inputs'/('carbon'+suffix)
    source.parent.mkdir()
    original = write_fchk(source, checkpoint())
    run_dir = tmp_path/'runs'
    run_dir.mkdir()
    monkeypatch.chdir(run_dir)
    seen = []
    def convert(args, **kwargs):
        assert args[0] == 'unfchk' and Path(args[1]) == source
        Path(args[2]).write_bytes(b'checkpoint orbitals')
    monkeypatch.setattr(fchk2chk.job_timeout, 'run', convert)
    def execute(name, *args):
        text = Path(name+'.com').read_text()
        seen.append(text)
        assert 'ChkBasis' in text and 'Guess=Read' in text
        assert Path(name+'.chk').read_bytes() == b'checkpoint orbitals'
        Path(name+'.log').write_text('Normal termination of Gaussian\n')
        return 'normal'
    monkeypatch.setattr(GaussianRunPack.gaussian_run.Exe_Gaussian, 'exe_Gaussian', execute)
    monkeypatch.setattr(GaussianRunPack.gaussian_run.chk2fchk, 'Get_chklist', lambda *args: None)
    monkeypatch.setattr(GaussianRunPack.GaussianDFTRun, 'Extract_values', lambda *args: {'energy': [-37., 0.]})
    job = GaussianRunPack.GaussianDFTRun('B3LYP', 'FCHK', 1, 'energy', source)
    job.mem = '1GB'
    assert job.run_gaussian()['log'] == 'normal'
    assert len(seen) == 1
    assert source.read_text() == original
    assert not source.with_suffix('.chk').exists()


def test_failed_unfchk_never_reuses_stale_chk(tmp_path, monkeypatch):
    target = tmp_path/'guess.chk'
    target.write_bytes(b'stale')
    def fail(args, **kwargs): raise subprocess.CalledProcessError(1, args)
    monkeypatch.setattr(fchk2chk.job_timeout, 'run', fail)
    with pytest.raises(subprocess.CalledProcessError):
        fchk2chk.convert('source.fchk', target)
    assert not target.exists()


@pytest.mark.parametrize('unrestricted', [False, True])
@pytest.mark.skipif(not os.environ.get('QCFOREVER_TEST_ORCA') or not shutil.which('orca'),
                    reason='Set QCFOREVER_TEST_ORCA=1 to run installed ORCA')
def test_installed_orca_imported_5d_guess_matches_reference(tmp_path, monkeypatch, unrestricted):
    """A normalized He S+5D basis has an exactly known initial density."""
    from qcforever.util import job_timeout
    from qcforever.orca_run.parse_log import parse_log
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(check_resource, 'respec_cores', lambda n: n)
    f = checkpoint()
    f.update({'Atomic numbers': [2], 'Nuclear charges': [2.],
              'Number of alpha electrons': 1, 'Number of beta electrons': 1,
              'Shell types': [0, -2], 'Number of primitives per shell': [1,1],
              'Shell to atom map': [1,1], 'Primitive exponents': [.5,.7],
              'Contraction coefficients': [1.,1.],
              'P(S=P) Contraction coefficients': [0.,0.],
              'Number of basis functions': 6, 'Number of independent functions': 6,
              'Alpha Orbital Energies': list(np.linspace(-1,1,6)),
              'Alpha MO coefficients': list(np.eye(6).T.flatten())})
    if unrestricted:
        f['Beta MO coefficients'] = f['Alpha MO coefficients'][:]
        f['Beta Orbital Energies'] = f['Alpha Orbital Energies'][:]
    source = tmp_path/'helium.fchk'
    write_fchk(source, f)
    orbitals = fchk_orbitals.read(source)
    job = OrcaRunPack.OrcaDFTRun('HF', 'FCHK', 1, 'energy', source)
    job.mem = '1GB'
    job._fchk_orbitals = orbitals
    with job_timeout.overall_timeout(45):
        guess = orbitals.orca_gbw(job.executable)
        energies = []
        for name, moread in [('imported', guess), ('reference', None)]:
            job.make_input(name, orbitals.atoms, orbitals.xyz, 0, 1, moread=moread)
            assert ('UKS' in Path(name+'.inp').read_text()) == unrestricted
            Exe_Orca.exe_Orca(name, 30, job.executable)
            energies.append(parse_log(name+'.log').energy()[0])
    assert energies[0] == pytest.approx(energies[1], abs=1e-8)


@pytest.mark.skipif(not os.environ.get('QCFOREVER_TEST_ORCA') or not shutil.which('orca'),
                    reason='Set QCFOREVER_TEST_ORCA=1 to run installed ORCA')
def test_installed_orca_water_density_matches_source_without_scf_updates(tmp_path, monkeypatch):
    from qcforever.util import job_timeout
    from qcforever.orca_run.parse_log import parse_log
    source = Path(__file__).parent/'fixtures/fchk/water_rhf_6-31gstar_5d.fchk'
    orbitals = fchk_orbitals.read(source)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(check_resource, 'respec_cores', lambda n: n)
    job = OrcaRunPack.OrcaDFTRun('HF', 'FCHK', 1, 'energy', source)
    job.mem = '1GB'
    job._fchk_orbitals = orbitals
    with job_timeout.overall_timeout(45):
        guess = orbitals.orca_gbw(job.executable)
        job.make_input('density', orbitals.atoms, orbitals.xyz, 0, 1, moread=guess)
        path = Path('density.inp')
        path.write_text(path.read_text().replace('TightSCF', 'VeryTightSCF NoIter'))
        Exe_Orca.exe_Orca('density', 30, job.executable)
    assert parse_log('density.log').energy()[0] == pytest.approx(-76.00915487179412, abs=1e-7, rel=0)
