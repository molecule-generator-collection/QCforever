"""Physical invariants, origin shifts, and backend output/input contracts."""

import numpy as np
import pytest

from qcforever.gaussian_run import GaussianRunPack
from qcforever.gamess_run import GamessRunPack
from qcforever.util.quadrupole import (
    BOHR_TO_ANGSTROM,
    E_ANGSTROM_TO_DEBYE,
    quadrupole_properties,
    tensor_from_components,
)


@pytest.fixture
def extract_gaussian_dipole(tmp_path):
    def extract(lines):
        path = tmp_path / "dipole.log"
        path.write_text("\n".join(lines))
        parser = GaussianRunPack.gaussian_run.parse_log.parse_log(path)
        return parser.Extract_DipoleMoment(lines)
    return extract


@pytest.fixture
def extract_gaussian_quadrupole(tmp_path, monkeypatch):
    # Bypass only the existing Link1 classifier, not tensor extraction or math.
    parser = GaussianRunPack.gaussian_run.parse_log.parse_log
    monkeypatch.setattr(parser, "Check_task", lambda self: ({"ts": 0}, [self.lines]))
    runner = GaussianRunPack.GaussianDFTRun("B3LYP", "STO-3G", 1, "dipole", "job.xyz")

    def extract(lines, fchk_path=None):
        (tmp_path / "job.log").write_text("\n".join(lines))
        if fchk_path is not None:
            (tmp_path / "job.fchk").write_text(fchk_path.read_text())
        return runner.Extract_values(str(tmp_path / "job"), {"dipole": True}, [], [])
    return extract


@pytest.fixture
def extract_gamess_quadrupole(tmp_path):
    runner = GamessRunPack.GamessDFTRun("B3LYP", "STO-3G", 1, "dipole", "job.xyz")

    def extract(lines):
        if not any(line.strip() == "ELECTROSTATIC MOMENTS" for line in lines):
            lines = gamess_block(quadrupole=False)[:-1] + lines + [" ...... END OF PROPERTY EVALUATION ......"]
        (tmp_path / "job.log").write_text("\n".join(lines) + "\n")
        return runner.Extract_values(str(tmp_path / "job"), {"dipole": True})
    return extract


def point_moments(coordinates, charges):
    dipole = E_ANGSTROM_TO_DEBYE * np.einsum("i,ij->j", charges, coordinates)
    quad = E_ANGSTROM_TO_DEBYE * np.einsum("i,ij,ik->jk", charges, coordinates, coordinates)
    return dipole, quad


def gaussian_block(coordinates, charges, masses=(1.0, 2.0, 12.0), orientation="Standard"):
    dipole, tensor = point_moments(coordinates, charges)
    lines = [" Charge = %g Multiplicity = 1" % sum(charges),
             " AtmWgt= " + " ".join(map(str, masses)),
             f" {orientation} orientation:", " ---------",
             " Center Atomic Atomic Coordinates (Angstroms)",
             " Number Number Type X Y Z", " ---------"]
    for i, xyz in enumerate(coordinates):
        lines.append(f" {i+1} 1 0 " + " ".join(f"{v:.12e}" for v in xyz))
    lines += [" ---------", " Charge= %g electrons" % sum(charges),
              " Dipole moment (field-independent basis, Debye):",
              " X= %.12e Y= %.12e Z= %.12e Tot= %.12e" % (*dipole, np.linalg.norm(dipole)),
              " Quadrupole moment (field-independent basis, Debye-Ang):",
              " XX= %.12e YY= %.12e ZZ= %.12e" % tuple(np.diag(tensor)),
              " XY= %.12e XZ= %.12e YZ= %.12e" % (tensor[0, 1], tensor[0, 2], tensor[1, 2]),
              " Traceless Quadrupole moment (field-independent basis, Debye-Ang):",
              " XX= 999 YY= 999 ZZ= 999", " XY= 999 XZ= 999 YZ= 999"]
    return lines


def test_gaussian_parser_returns_native_tensor_and_origin_data(tmp_path):
    xyz = np.eye(3) + [2, -1, 3]
    charges = [1, -1, 1]
    lines = gaussian_block(xyz, charges)
    path = tmp_path / "job.log"
    path.write_text("\n".join(lines))
    parser = GaussianRunPack.gaussian_run.parse_log.parse_log(path)
    moments = parser.Extract_QuadrupoleMoment(parser.lines)
    dipole, native_tensor = point_moments(xyz, charges)
    np.testing.assert_allclose(moments["primitive"], native_tensor)
    np.testing.assert_allclose(moments["dipole"], dipole)
    np.testing.assert_allclose(moments["coordinates"], xyz)
    np.testing.assert_allclose(moments["masses"], [1, 2, 12])
    assert moments["charge"] == 1
    assert "quadrupole_iso" not in moments
    assert "quadrupole_anisotropy" not in moments


@pytest.mark.parametrize("charges", [[0.5, -0.2, -0.3], [0.7, 0.2, 0.1], [1, 0, -2]])
def test_gaussian_com_origin_is_translation_invariant(charges, extract_gaussian_quadrupole):
    xyz = np.array([[0.2, 1.0, -0.4], [2.0, -0.5, 0.3], [-1.0, 0.5, 1.2]])
    masses = [1, 2, 12]
    center = np.average(xyz, axis=0, weights=masses)
    _, expected = point_moments(xyz - center, charges)
    original = extract_gaussian_quadrupole(gaussian_block(xyz, charges, masses))
    moved = extract_gaussian_quadrupole(gaussian_block(xyz + [3, -2, 1], charges, masses))
    np.testing.assert_allclose(tensor_from_components(original["quadrupole"]), expected, atol=1e-10)
    np.testing.assert_allclose(tensor_from_components(moved["quadrupole"]), expected, atol=1e-10)
    assert original["quadrupole_iso"] == pytest.approx(np.trace(expected) / 3)
    assert moved["quadrupole_anisotropy"] == pytest.approx(original["quadrupole_anisotropy"])
    assert original["quadrupole_info"]["origin"] == "center_of_mass"
    assert original["quadrupole_info"]["origin_coordinates"] == pytest.approx(center)


def test_scalar_invariants_and_offdiagonal_weight():
    tensor = np.array([[-10, 0.3, 0], [0.3, -12, -0.2], [0, -0.2, -8]])
    rotation, _ = np.linalg.qr(np.array([[1., 2, 3], [3, 1, 2], [2, 3, 1]]))
    expected_delta = np.sqrt(((2**2 + (-4)**2 + 2**2) + 6 * (0.3**2 + 0.2**2)) / 2)
    original = quadrupole_properties(primitive=tensor)
    rotated = quadrupole_properties(primitive=rotation @ tensor @ rotation.T)
    assert original["quadrupole_iso"] == -10
    assert original["quadrupole_anisotropy"] == pytest.approx(expected_delta)
    assert rotated["quadrupole_iso"] == pytest.approx(original["quadrupole_iso"])
    assert rotated["quadrupole_anisotropy"] == pytest.approx(expected_delta)
    assert np.trace(tensor_from_components(original["quadrupole_traceless"])) == pytest.approx(0)
    isotropic = quadrupole_properties(primitive=np.eye(3) * -5)
    assert isotropic["quadrupole_anisotropy"] == 0


def test_gaussian_uses_last_block_and_preserves_native_dipole(extract_gaussian_dipole, extract_gaussian_quadrupole):
    xyz = np.array([[1, 0, 0], [0, 2, 0], [0, 0, 1.]])
    first = gaussian_block(xyz, [1, -1, 0])
    final = gaussian_block(xyz + 2, [0, 0, 1], orientation="Input")
    # Fortran exponents and an unrelated X= line must not break dipole extraction.
    lines = [line.replace("e+", "D+").replace("e-", "D-") for line in first + final]
    lines.append(" X= 999 Y= 999 Z= 999 Tot= 999")
    expected_dipole, _ = point_moments(xyz + 2, [0, 0, 1])
    assert extract_gaussian_dipole(lines) == pytest.approx([*expected_dipole, np.linalg.norm(expected_dipole)])
    actual = extract_gaussian_quadrupole(lines)
    expected = extract_gaussian_quadrupole(final)
    assert actual["quadrupole"] == expected["quadrupole"]


@pytest.mark.parametrize("missing", ["quadrupole", "component", "geometry", "mass"])
def test_gaussian_missing_data_is_explicit_and_does_not_reuse_earlier_values(missing, extract_gaussian_quadrupole):
    xyz = np.eye(3)
    earlier = gaussian_block(xyz, [1, -1, 0])
    final = gaussian_block(xyz * 2, [0, 1, -1])
    if missing == "quadrupole":
        final = final[:final.index(" Quadrupole moment (field-independent basis, Debye-Ang):")]
    elif missing == "component":
        final = [line for line in final if not line.startswith(" XY=")]
    elif missing == "geometry":
        earlier = []
        final = final[-9:]
    else:
        earlier = []
        final = [line for line in final if "AtmWgt" not in line]
    result = extract_gaussian_quadrupole(earlier + final)
    assert result["quadrupole"] is None
    assert result["quadrupole_anisotropy"] is None
    assert result["quadrupole_info"]["status"] == "unavailable"
    assert result["quadrupole_info"]["reason"]


def test_gaussian_fchk_mass_fallback_preserves_isotopes(tmp_path, extract_gaussian_quadrupole):
    xyz = np.eye(3)
    lines = [line for line in gaussian_block(xyz, [1, -1, 0]) if "AtmWgt" not in line]
    fchk = tmp_path / "job.fchk"
    fchk.write_text("Real atomic weights                        R   N=           3\n 1.0 2.0 12.0\n")
    actual = extract_gaussian_quadrupole(lines, fchk)
    expected = extract_gaussian_quadrupole(gaussian_block(xyz, [1, -1, 0]))
    assert actual["quadrupole"] == expected["quadrupole"]


def test_gaussian_does_not_reuse_geometry_after_malformed_orientation(extract_gaussian_quadrupole):
    lines = gaussian_block(np.eye(3), [1, -1, 0])
    final = gaussian_block(np.eye(3) * 2, [1, -1, 0])
    final = [" 2 1 0 ***** 0.0 0.0" if line.startswith(" 2 1 0") else line
             for line in final]
    result = extract_gaussian_quadrupole(lines + final)
    assert result["quadrupole"] is None
    assert "orientation" in result["quadrupole_info"]["reason"]


def gamess_block(scale=1, label="ELECTROSTATIC MOMENTS", quadrupole=True):
    lines = [" ---------------------", " " + label, " ---------------------", "",
             " POINT 1 X Y Z (BOHR) CHARGE", " 0.000000 0.000000 0.123456 0.00 (A.U.)",
             " DX DY DZ /D/ (DEBYE)", " 1.0 2.0 2.0 3.0"]
    if quadrupole:
        lines += [" QXX QYY QZZ (BUCKINGHAMS)",
                  f" {3*scale:.6f} {-scale:.6f} {-2*scale:.6f}",
                  " QXY QXZ QYZ (BUCKINGHAMS)",
                  f" {0.3*scale:.6f} {0.6*scale:.6f} {-0.9*scale:.6f}"]
    return lines + [" ...... END OF PROPERTY EVALUATION ......"]


def test_gamess_parser_returns_native_buckingham_tensor(tmp_path):
    lines = gamess_block()[:-1] + [" PRINCIPAL QUADRUPOLE MOMENTS",
                                   " QXX QYY QZZ (BUCKINGHAMS)", " 999 999 999",
                                   " ...... END OF PROPERTY EVALUATION ......"]
    path = tmp_path / "job.log"
    path.write_text("\n".join(lines) + "\n")
    parser = GamessRunPack.gamess_run.parse_log.parse_log(path)
    moments = parser.getQuadrupoleMoment(parser.getBlock("ELECTROSTATIC MOMENTS")[-1])
    expected = np.array([[3, 0.3, 0.6], [0.3, -1, -0.9], [0.6, -0.9, -2]])
    np.testing.assert_allclose(moments["tensor"], expected)
    assert "quadrupole_iso" not in moments
    assert "quadrupole_anisotropy" not in moments


def test_gamess_buckingham_normalization_and_missing_primitive_trace(extract_gamess_quadrupole):
    block = gamess_block() + [" PRINCIPAL QUADRUPOLE MOMENTS",
                              " QXX QYY QZZ (BUCKINGHAMS)", " 999 999 999"]
    result = extract_gamess_quadrupole(block)
    expected = np.array([[2, 0.2, 0.4], [0.2, -2/3, -0.6], [0.4, -0.6, -4/3]])
    np.testing.assert_allclose(tensor_from_components(result["quadrupole_traceless"]), expected)
    assert result["quadrupole_anisotropy"] == pytest.approx(np.sqrt(1.5 * np.sum(expected**2)))
    assert result["quadrupole"] is None
    assert result["quadrupole_iso"] is None
    assert result["quadrupole_info"]["origin_coordinates"] == pytest.approx([0, 0, 0.123456 * BOHR_TO_ANGSTROM])


def test_gamess_single_header_column_order_and_fortran_exponents(extract_gamess_quadrupole):
    block = [" QXX QXY QXZ QYY QYZ QZZ (BUCKINGHAMS)", "",
             " 3.0D+00 3.0D-01 6.0D-01 -1.0D+00 -9.0D-01 -2.0D+00"]
    result = extract_gamess_quadrupole(block)
    expected = extract_gamess_quadrupole(gamess_block())
    assert result["quadrupole_traceless"] == expected["quadrupole_traceless"]


@pytest.mark.parametrize("solvent,label", [("0", "ELECTROSTATIC MOMENTS"), ("WATER", "SYSTEM ELECTROSTATIC MOMENTS")])
def test_gamess_extraction_selects_last_matching_block(tmp_path, monkeypatch, solvent, label):
    monkeypatch.chdir(tmp_path)
    lines = gamess_block() + gamess_block(2, label)
    (tmp_path / "job.log").write_text("\n".join(lines) + "\n")
    runner = GamessRunPack.GamessDFTRun("B3LYP", "STO-3G", 1, "dipole", "job.xyz", solvent=solvent)
    result = runner.Extract_values("job", {"dipole": True})
    assert result["dipole"] == 3.0
    assert result["quadrupole_traceless"]["XX"] == 4.0
    assert runner.Extract_values("job", {}) == {}
    # A final dipole-only block must not acquire an earlier step's quadrupole.
    (tmp_path / "job.log").write_text("\n".join(lines + gamess_block(3, label, quadrupole=False)) + "\n")
    result = runner.Extract_values("job", {"dipole": True})
    assert result["dipole"] == 3.0
    assert result["quadrupole_anisotropy"] is None


def test_gaussian_extraction_is_additive_and_scoped_to_requested_state(tmp_path, monkeypatch, extract_gaussian_dipole):
    monkeypatch.chdir(tmp_path)
    lines = gaussian_block(np.eye(3), [1, -1, 0])
    # Mock only the pre-existing Link1 classifier, so the test needs no executable.
    parser = GaussianRunPack.gaussian_run.parse_log.parse_log
    monkeypatch.setattr(parser, "Check_task", lambda self: ({"ts": 0}, [lines]))
    (tmp_path / "job.log").write_text("unused")
    runner = GaussianRunPack.GaussianDFTRun("B3LYP", "STO-3G", 1, "dipole", "job.xyz")
    result = runner.Extract_values("job", {"dipole": True}, [], [])
    assert result["dipole"] == pytest.approx(extract_gaussian_dipole(lines))
    assert result["quadrupole_info"]["status"] == "ok"
    assert runner.Extract_values("job", {}, [], []) == {}


@pytest.mark.parametrize("option,expected", [("energy dipole", True), ("DiPoLe", True),
                                            ("energy", False), ("uv=/tmp/dipole.dat", False)])
@pytest.mark.parametrize("run_type", ["ENERGY", "OPTIMIZE"])
def test_gamess_input_requests_quadrupole_only_for_dipole(tmp_path, option, expected, run_type):
    runner = GamessRunPack.GamessDFTRun("B3LYP", "STO-3G", 1, option, "job.xyz")
    path = tmp_path / "job.inp"
    runner.make_input(run_type, 0, 1, path, ["H", "H"], [0, 0], [0, 0], [0, 0.7])
    content = path.read_text()
    assert (" $ELMOM IEMOM=2 WHERE=COMASS $END" in content) is expected
    assert content.count("$ELMOM") == int(expected)


def test_gamess_dat_restart_with_pcm_requests_quadrupole(tmp_path):
    dat = tmp_path / "job.dat"
    dat.write_text(" $DATA\njob\nC1\nH 1.0 0.0 0.0 0.0\nH 1.0 0.0 0.0 0.7\n $END\n"
                   " $VEC\n 1 1 1.0 0.0\n 2 1 0.0 1.0\n $END\n")
    runner = GamessRunPack.GamessDFTRun("B3LYP", "STO-3G", 1, "dipole uv", "job.xyz", solvent="WATER")
    path = tmp_path / "excited.inp"
    runner.make_input("ENERGY", 0, 1, path, TDDFT=True, datfile=dat, solvent="WATER")
    content = path.read_text()
    assert content.count("$ELMOM") == 1
    assert " $ELMOM IEMOM=2 WHERE=COMASS $END" in content
    assert "$PCM" in content
    assert "$VEC" in content
