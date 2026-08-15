"""Cleanup helpers that retain only files needed to restart QCforever jobs."""

from __future__ import annotations

import os
import re
import shutil
from pathlib import Path


def _remove_path(path):
    if path.is_symlink() or path.is_file():
        path.unlink(missing_ok=True)
    elif path.is_dir():
        shutil.rmtree(path)


def _clean_job_directory(job_directory, restart_suffixes, preserve_pickle=False):
    job_directory = Path(job_directory)
    if not job_directory.is_dir():
        return

    preserved_suffixes = {suffix.lower() for suffix in restart_suffixes}
    if preserve_pickle:
        preserved_suffixes.add(".pkl")

    for path in job_directory.iterdir():
        if path.is_file() and path.suffix.lower() in preserved_suffixes:
            continue
        _remove_path(path)


def cleanup_gaussian(job_directory, preserve_pickle=False, preserve_xyz=False):
    """Keep Gaussian restart files and remove all other generated files."""
    restart_suffixes = {".chk", ".fchk"}
    if preserve_xyz:
        # ``restart=False`` intentionally retains geometry rather than the
        # electronic checkpoint; an xyz can restart QCforever geometrically.
        restart_suffixes.add(".xyz")
    _clean_job_directory(
        job_directory,
        restart_suffixes=restart_suffixes,
        preserve_pickle=preserve_pickle,
    )


def _read_gamess_scratch_directories(rungms_path):
    """Read SCR and USERSCR assignments from a rungms script."""
    directories = []
    try:
        lines = Path(rungms_path).read_text(errors="replace").splitlines()
    except (OSError, TypeError):
        return directories

    for variable in ("SCR", "USERSCR"):
        directory = None
        pattern = re.compile(
            rf"^\s*set\s+{variable}\s*=\s*['\"]?([^'\"\s]+)",
            re.IGNORECASE,
        )
        for line in lines:
            match = pattern.match(line)
            if match:
                directory = Path(os.path.expandvars(match.group(1))).expanduser()
        if directory is not None and directory not in directories:
            directories.append(directory)
    return directories


def cleanup_gamess(
    job_directory,
    jobname,
    preserve_pickle=False,
    scratch_directories=None,
):
    """Collect GAMESS .dat files locally, then remove other job artifacts."""
    job_directory = Path(job_directory)
    job_directory.mkdir(parents=True, exist_ok=True)

    if scratch_directories is None:
        rungms_path = shutil.which("rungms")
        scratch_directories = _read_gamess_scratch_directories(rungms_path)

    # GAMESS writes restart data to USERSCR (and some installations use SCR).
    # Move it before deleting scratch files so a timed-out job can be resumed.
    scratch_matches = []
    for scratch_directory in scratch_directories:
        scratch_directory = Path(scratch_directory)
        if not scratch_directory.is_dir():
            continue
        scratch_matches.extend(
            path for path in scratch_directory.iterdir()
            if path.name.startswith((f"{jobname}.", f"{jobname}_"))
        )

    for path in scratch_matches:
        if path.is_file() and path.suffix.lower() == ".dat":
            destination = job_directory / path.name
            if path.resolve() == destination.resolve():
                continue
            if destination.exists():
                _remove_path(destination)
            shutil.move(str(path), str(destination))

    for path in scratch_matches:
        if (
            path.suffix.lower() == ".dat"
            and path.parent.resolve() == job_directory.resolve()
        ):
            continue
        if path.exists() or path.is_symlink():
            _remove_path(path)

    _clean_job_directory(
        job_directory,
        restart_suffixes={".dat"},
        preserve_pickle=preserve_pickle,
    )
