import os
import time

import pytest

from qcforever.util import job_timeout
from qcforever.gamess_run.GamessRunPack import GamessDFTRun
from qcforever.gaussian_run.GaussianRunPack import GaussianDFTRun


def test_overall_timeout_interrupts_python_work():
    if not hasattr(job_timeout.signal, "SIGALRM"):
        pytest.skip("SIGALRM is required to interrupt Python work")

    started = time.monotonic()
    with pytest.raises(job_timeout.QCforeverTimeoutError):
        with job_timeout.overall_timeout(0.1):
            time.sleep(10)

    assert time.monotonic() - started < 2


def test_timeout_terminates_subprocess_group():
    if os.name != "posix":
        pytest.skip("process-group assertion is POSIX-specific")

    started = time.monotonic()
    with pytest.raises(job_timeout.QCforeverTimeoutError):
        with job_timeout.overall_timeout(0.1):
            process = job_timeout.popen(["/bin/sh", "-c", "sleep 10"])
            job_timeout.wait(process)

    assert process.poll() is not None
    assert time.monotonic() - started < 2


def test_nonpositive_limit_is_unlimited():
    with job_timeout.overall_timeout(0):
        time.sleep(0.01)


def test_timeout_cannot_be_permanently_swallowed_by_inner_handler():
    if not hasattr(job_timeout.signal, "SIGALRM"):
        pytest.skip("SIGALRM is required to interrupt Python work")

    with pytest.raises(job_timeout.QCforeverTimeoutError):
        with job_timeout.overall_timeout(0.05):
            try:
                time.sleep(1)
            except job_timeout.QCforeverTimeoutError:
                pass
            time.sleep(1)


@pytest.mark.parametrize(
    ("run_class", "method_name", "worker_name"),
    [
        (GaussianDFTRun, "run_gaussian", "_run_gaussian"),
        (GamessDFTRun, "run_gamess", "_run_gamess"),
    ],
)
def test_top_level_run_returns_timeout_result(run_class, method_name, worker_name):
    calculation = run_class.__new__(run_class)
    calculation.timejob = 0.05
    setattr(calculation, worker_name, lambda: time.sleep(1))

    result = getattr(calculation, method_name)()

    assert result["log"] == "timeout"
    assert "overall wall-clock time limit" in result["error"]
