"""Launch the ORCA driver directly; ORCA starts its own MPI workers."""

import os
import shutil
import subprocess
from collections import deque
from pathlib import Path

from qcforever.util import check_resource, job_timeout


class OrcaExecutionError(RuntimeError):
    def __init__(self, jobname, returncode, logfile, reason):
        self.jobname = jobname
        self.returncode = returncode
        self.logfile = str(Path(logfile).resolve())
        with open(logfile, errors='replace') as stream:
            tail = ''.join(deque(stream, maxlen=30)).strip()
        super().__init__(f'ORCA job {jobname!r} failed: {reason}. '
                         f'Exit code: {returncode}. Log: {self.logfile}\n{tail}')


def exe_Orca(jobname, timexe=None, executable='orca'):
    """Save partial output on failure/timeout and respect the overall deadline."""
    program = shutil.which(str(executable))
    if program is None:
        raise FileNotFoundError(f'ORCA executable {executable!r} was not found')
    # The absolute installation path lets the driver locate sibling programs.
    program = str(Path(program).resolve())
    logfile = Path(f'{jobname}.log')
    environment = os.environ.copy()
    environment.update(check_resource.native_thread_environment(1))
    with logfile.open('w') as stream:
        try:
            result = job_timeout.run(
                [program, f'{jobname}.inp'], stdout=stream,
                stderr=subprocess.STDOUT, env=environment,
                timeout=timexe if timexe and timexe > 0 else None,
            )
        except subprocess.TimeoutExpired as exc:
            raise job_timeout.QCforeverTimeoutError(
                f'ORCA calculation {jobname!r} exceeded timexe ({timexe:g} s); '
                f'log: {logfile.resolve()}') from exc
    with logfile.open(errors='replace') as stream:
        tail = ''.join(deque(stream, maxlen=100))
    if result.returncode or 'ORCA TERMINATED NORMALLY' not in tail:
        raise OrcaExecutionError(jobname, result.returncode, logfile,
                                 'unsuccessful exit or missing normal termination')
    return 'normal'
