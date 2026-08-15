from pathlib import Path

from qcforever.util import job_cleanup


def _write(path, text="data"):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def test_gaussian_cleanup_keeps_only_restart_and_requested_pickle(tmp_path):
    job_directory = tmp_path / "molecule"
    _write(job_directory / "state.chk")
    _write(job_directory / "state.fchk")
    _write(job_directory / "state.log")
    _write(job_directory / "state.com")
    _write(job_directory / "molecule.pkl")
    _write(job_directory / "conformer" / "xtbopt.xyz")

    job_cleanup.cleanup_gaussian(job_directory, preserve_pickle=True)

    assert {path.name for path in job_directory.iterdir()} == {
        "state.chk",
        "state.fchk",
        "molecule.pkl",
    }


def test_gamess_cleanup_collects_dat_and_cleans_both_scratch_dirs(tmp_path):
    job_directory = tmp_path / "molecule"
    scr = tmp_path / "scr"
    userscr = tmp_path / "userscr"
    _write(job_directory / "molecule.inp")
    _write(job_directory / "molecule.log")
    _write(job_directory / "molecule.pkl")
    _write(scr / "molecule.F05")
    _write(userscr / "molecule.dat", "restart")
    _write(userscr / "molecule_TD.dat", "td restart")
    _write(userscr / "molecule2.dat", "another job")
    _write(userscr / "unrelated.dat", "keep")

    job_cleanup.cleanup_gamess(
        job_directory,
        "molecule",
        preserve_pickle=True,
        scratch_directories=[scr, userscr],
    )

    assert {path.name for path in job_directory.iterdir()} == {
        "molecule.dat",
        "molecule_TD.dat",
        "molecule.pkl",
    }
    assert not (scr / "molecule.F05").exists()
    assert (userscr / "molecule2.dat").exists()
    assert (userscr / "unrelated.dat").exists()


def test_pickle_is_removed_unless_explicitly_requested(tmp_path):
    job_directory = tmp_path / "molecule"
    _write(job_directory / "state.chk")
    _write(job_directory / "molecule.pkl")

    job_cleanup.cleanup_gaussian(job_directory, preserve_pickle=False)

    assert [path.name for path in job_directory.iterdir()] == ["state.chk"]


def test_gaussian_xyz_is_kept_only_when_selected_for_geometry_restart(tmp_path):
    job_directory = tmp_path / "molecule"
    _write(job_directory / "state.xyz")
    _write(job_directory / "state.log")

    job_cleanup.cleanup_gaussian(job_directory, preserve_xyz=True)

    assert [path.name for path in job_directory.iterdir()] == ["state.xyz"]


def test_gamess_scratch_directories_use_final_assignment(tmp_path, monkeypatch):
    monkeypatch.setenv("QCFOREVER_TEST_SCR", str(tmp_path / "expanded"))
    rungms = tmp_path / "rungms"
    rungms.write_text(
        "set SCR=/old/scr\n"
        "set SCR=$QCFOREVER_TEST_SCR\n"
        "set USERSCR = '/user/scr'\n"
    )

    directories = job_cleanup._read_gamess_scratch_directories(rungms)

    assert directories == [tmp_path / "expanded", Path("/user/scr")]
