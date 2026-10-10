"""Execute independent one-core optimization blocks with bounded concurrency.

The scheduler supplies a batch of candidate IDs. Each candidate owns its native
state and files; this pool supplies CPU slots, deadline propagation and cleanup.
Threads only supervise external calculations: native adapters must not change
the Python process's working directory or environment. Results are returned in
dispatch order, independent of native completion order.
"""
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
from contextvars import copy_context
import os
import threading

from qcforever.util import job_timeout


def resolve_block_workers(allocated_cores, candidate_count, requested='auto'):
    """Never exceed nproc, the candidate pool, or the scheduler CPU allocation."""
    limits = [allocated_cores, candidate_count]
    if requested != 'auto':
        limits.append(requested)
    if hasattr(os, 'sched_getaffinity'):
        limits.append(len(os.sched_getaffinity(0)))
    return max(1, min(limits))


class CandidateBlockExecutor:
    """A reusable batch executor shared by LAQA, SR and SH, for either backend."""

    def __init__(self, runners, workers, check_resources=None):
        self.runners = runners
        self.workers = workers
        self.check_resources = check_resources
        self.cancelled = threading.Event()
        self.pool = None
        self.owner_scope = None
        self.owner = None
        self.futures = []
        available = sorted(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else []
        self.cpu_slots = available[:workers] if available else [None] * workers

    def __enter__(self):
        self.owner_scope = job_timeout.child_process_scope()
        self.owner = self.owner_scope.__enter__()
        if self.workers > 1:
            self.pool = ThreadPoolExecutor(max_workers=self.workers, thread_name_prefix='conformer-block')
        return self

    def advance_batch(self, candidates):
        """Run one block per candidate, then return every measured result."""
        if len(candidates) > self.workers or len(candidates) != len(set(candidates)):
            raise ValueError('Invalid native candidate batch')
        self._check_resources()
        if self.pool is None:
            return {key: self._advance_one(key, self.cpu_slots[0]) for key in candidates}
        futures = {}
        for slot, key in enumerate(candidates):
            # ContextVars do not propagate into ThreadPoolExecutor by default.
            # Each task gets its own copy, sharing the same deadline/child owner.
            context = copy_context()
            futures[key] = self.pool.submit(context.run, self._advance_one, key, self.cpu_slots[slot])
            self.futures.append(futures[key])
        try:
            pending = set(futures.values())
            while pending:
                completed, pending = wait(pending, timeout=.1, return_when=FIRST_COMPLETED)
                self._check_resources()
                for future in completed:
                    future.result()  # Abort promptly on infrastructure errors.
            results = {key: future.result() for key, future in futures.items()}
            self.futures.clear()
            return results
        except BaseException:
            self._cancel_active_blocks()
            raise

    def _cancel_active_blocks(self):
        """Also catch a child launched concurrently with the cancellation signal."""
        self.cancelled.set()
        for future in self.futures:
            future.cancel()
        while any(not future.done() for future in self.futures):
            job_timeout.terminate_children(self.owner)
            wait(self.futures, timeout=.1)
        job_timeout.terminate_children(self.owner)

    def _advance_one(self, key, cpu_id):
        if self.cancelled.is_set():
            raise InterruptedError('Candidate batch cancelled')
        remaining = job_timeout.remaining_time()
        if remaining is not None and remaining <= 0:
            raise job_timeout.QCforeverTimeoutError('Graybox overall deadline exceeded')
        runner = self.runners[key]
        runner.cpu_id = cpu_id
        runner.cancel_requested = self.cancelled.is_set
        return runner.advance()

    def _check_resources(self):
        remaining = job_timeout.remaining_time()
        if remaining is not None and remaining <= 0:
            raise job_timeout.QCforeverTimeoutError('Graybox overall deadline exceeded')
        if self.check_resources is not None:
            self.check_resources()

    def __exit__(self, exc_type, exc, traceback):
        # Join supervisors before closing their files or native candidate state.
        try:
            if exc_type is not None:
                self._cancel_active_blocks()
            if self.pool is not None:
                self.pool.shutdown(wait=True, cancel_futures=True)
        finally:
            self.owner_scope.__exit__(exc_type, exc, traceback)
