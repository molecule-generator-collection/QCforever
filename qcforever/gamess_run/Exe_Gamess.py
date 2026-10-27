import os
import sys
import shutil
import subprocess
from collections import deque
from pathlib import Path

from qcforever.util import job_timeout
from qcforever.util.gamess_scratch import resolve_scratch_directories


class GamessExecutionError(RuntimeError):
    """A failed launcher or unavailable required GAMESS restart output."""

    def __init__(self, jobname, returncode, logfile, reason):
        self.jobname = jobname
        self.returncode = returncode
        self.logfile = Path(logfile).resolve()
        try:
            with self.logfile.open(errors="replace") as log:
                self.log_tail = "".join(deque(log, maxlen=30)).strip()
        except OSError:
            self.log_tail = ""
        message = (
            f"GAMESS job {jobname!r} failed: {reason}. "
            f"rungms exit code: {returncode}. Log: {self.logfile}"
        )
        if "limit: stacksize: Can't set limit" in self.log_tail:
            message += (
                "\nThe launcher could not set the stack size. Check the job's "
                "soft/hard stack limits and the limit stacksize setting in "
                "the site launcher; QCforever cannot raise the system hard limit."
            )
        if self.log_tail:
            message += "\nLast GAMESS log lines:\n" + self.log_tail
        super().__init__(message)


def read_rungms(gmspath):
    return tuple(str(path) for path in resolve_scratch_directories(gmspath))


def exe_Gamess(jobname, gamessversion, nproc):
    job_state = ""
    cwd = os.getcwd()
    print(cwd)
    gmspath = shutil.which('rungms')
    if gmspath is None:
        raise FileNotFoundError("GAMESS launcher 'rungms' was not found on PATH")
    scr, userscr = read_rungms(gmspath)
    for directory in dict.fromkeys((Path(scr), Path(userscr))):
        directory.mkdir(parents=True, exist_ok=True)
        for path in directory.iterdir():
            if path.name.startswith((jobname + '.', jobname + '_')) and path.is_file():
                print(path)
                path.unlink()
    infile = jobname + '.inp'
    outfile = jobname +'.log'
    datfile = jobname +'.dat'

    environment = os.environ.copy()
    environment.update(SCR=scr, USERSCR=userscr)
    # Write directly to disk so both the reported scratch paths and partial
    # calculation output survive a timeout or other interrupted execution.
    with open(outfile, 'wb') as log:
        result = job_timeout.run(
            [gmspath, infile, gamessversion, str(nproc)],
            stdout=log,
            stderr=subprocess.STDOUT,
            env=environment,
        )

    if result.returncode != 0:
        raise GamessExecutionError(
            jobname, result.returncode, outfile, "rungms exited unsuccessfully"
        )

    _, actual_userscr = resolve_scratch_directories(gmspath, log_path=outfile)
    path_datfile = actual_userscr / datfile
    if not path_datfile.is_file():
        raise GamessExecutionError(
            jobname, result.returncode, outfile,
            f"required restart file was not found at {path_datfile}; "
            "GAMESS may have stopped before producing it or written it elsewhere",
        )
    destination = os.path.join(cwd, datfile)
    if Path(path_datfile).resolve() != Path(destination).resolve():
        try:
            shutil.move(path_datfile, destination)
        except OSError as exc:
            raise GamessExecutionError(
                jobname, result.returncode, outfile,
                f"could not collect {path_datfile} into {destination}: {exc}",
            ) from exc

    return job_state 

if __name__ == '__main__':
    usage = f'Usage; {sys.argv[0]}'

    try:
        jobname = sys.argv[1]
    except:
        print (usage); sys.exit()

    print(exe_Gamess(jobname, '00', 2))

