"""Wall-clock timeout support for a complete QCforever calculation."""

from __future__ import annotations

import contextvars
import os
import signal
import subprocess
import threading
import time
from contextlib import contextmanager


class QCforeverTimeoutError(TimeoutError):
    """Raised when the wall-clock limit for a QCforever job is reached."""


_deadline = contextvars.ContextVar("qcforever_deadline", default=None)
_job_id = contextvars.ContextVar("qcforever_job_id", default=None)
_processes: dict[subprocess.Popen, object | None] = {}
_process_lock = threading.Lock()


def remaining_time():
    """Return seconds remaining in the current job, or ``None`` if unlimited."""
    deadline = _deadline.get()
    if deadline is None:
        return None
    return max(0.0, deadline - time.monotonic())


def terminate_process(process, grace_period=2.0):
    """Terminate a process and every descendant in its process group."""
    if process.poll() is not None:
        return
    try:
        if os.name == "posix":
            os.killpg(os.getpgid(process.pid), signal.SIGTERM)
        else:
            process.terminate()
        process.wait(timeout=grace_period)
    except (ProcessLookupError, subprocess.TimeoutExpired):
        if process.poll() is None:
            try:
                if os.name == "posix":
                    os.killpg(os.getpgid(process.pid), signal.SIGKILL)
                else:
                    process.kill()
                process.wait(timeout=grace_period)
            except ProcessLookupError:
                pass
            except subprocess.TimeoutExpired:
                # The OS will reap it asynchronously; SIGKILL cannot be ignored.
                pass


def terminate_children(job_id=None):
    """Terminate registered external processes belonging to one job."""
    with _process_lock:
        processes = tuple(
            process for process, owner in _processes.items()
            if job_id is None or owner is job_id
        )
    for process in processes:
        terminate_process(process)
    with _process_lock:
        for process in processes:
            _processes.pop(process, None)


def popen(*args, **kwargs):
    """Start and register a subprocess in its own process group."""
    if os.name == "posix":
        kwargs.setdefault("start_new_session", True)
    process = subprocess.Popen(*args, **kwargs)
    with _process_lock:
        _processes[process] = _job_id.get()
    return process


def wait(process):
    """Wait for a registered process while respecting the overall deadline."""
    timeout = remaining_time()
    try:
        if timeout is not None and timeout <= 0:
            raise subprocess.TimeoutExpired(process.args, timeout)
        return process.wait(timeout=timeout)
    except subprocess.TimeoutExpired as exc:
        terminate_process(process)
        raise QCforeverTimeoutError(
            "QCforever overall wall-clock time limit exceeded"
        ) from exc
    finally:
        with _process_lock:
            _processes.pop(process, None)


def run(*args, **kwargs):
    """Equivalent to ``subprocess.run`` with deadline and tree termination."""
    input_data = kwargs.pop("input", None)
    if input_data is not None:
        kwargs["stdin"] = subprocess.PIPE
    check = kwargs.pop("check", False)
    process = popen(*args, **kwargs)
    try:
        stdout, stderr = process.communicate(
            input=input_data,
            timeout=remaining_time(),
        )
    except subprocess.TimeoutExpired as exc:
        terminate_process(process)
        process.communicate()
        raise QCforeverTimeoutError(
            "QCforever overall wall-clock time limit exceeded"
        ) from exc
    finally:
        with _process_lock:
            _processes.pop(process, None)
    result = subprocess.CompletedProcess(process.args, process.returncode, stdout, stderr)
    if check:
        result.check_returncode()
    return result


@contextmanager
def overall_timeout(seconds):
    """Apply a wall-clock deadline and terminate registered child processes."""
    if seconds is None:
        yield
        return

    seconds = float(seconds)
    if seconds <= 0:
        yield
        return

    deadline_token = _deadline.set(time.monotonic() + seconds)
    current_job_id = object()
    job_token = _job_id.set(current_job_id)
    main_thread = threading.current_thread() is threading.main_thread()
    previous_handler = None
    timer = None

    def expire(*_):
        if main_thread and hasattr(signal, "SIGALRM"):
            signal.setitimer(signal.ITIMER_REAL, 0)
        terminate_children(current_job_id)
        if main_thread and hasattr(signal, "SIGALRM"):
            # Arm another one-shot alarm in case legacy code catches this one.
            signal.setitimer(signal.ITIMER_REAL, 0.1)
        raise QCforeverTimeoutError(
            f"QCforever overall wall-clock time limit ({seconds:g} s) exceeded"
        )

    if main_thread and hasattr(signal, "SIGALRM"):
        previous_handler = signal.getsignal(signal.SIGALRM)
        signal.signal(signal.SIGALRM, expire)
        signal.setitimer(signal.ITIMER_REAL, seconds)
    else:
        # In a worker thread Python cannot asynchronously raise an exception,
        # but terminating the active external calculation still unblocks it.
        timer = threading.Timer(seconds, terminate_children, args=(current_job_id,))
        timer.daemon = True
        timer.start()

    try:
        yield
        if remaining_time() == 0:
            expire()
    finally:
        if timer is not None:
            timer.cancel()
        if previous_handler is not None:
            signal.setitimer(signal.ITIMER_REAL, 0)
            signal.signal(signal.SIGALRM, previous_handler)
        terminate_children(current_job_id)
        _job_id.reset(job_token)
        _deadline.reset(deadline_token)
