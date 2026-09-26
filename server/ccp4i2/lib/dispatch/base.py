"""The two run-target interfaces. Structural (Protocol): nothing to subclass.

A deployment's target is any class with these methods; it is registered by
dotted path in ``CCP4I2_RUN_TARGETS`` and never imported until asked for, so
its platform SDK is loaded only where that platform is configured.
"""
from typing import Any, Dict, Protocol, runtime_checkable


@runtime_checkable
class JobTarget(Protocol):
    """Axis A: the whole CCP4i2 job runs somewhere.

    ``run_job`` starts (or queues) the job and returns the dict shape every
    caller of ``run_job_context_aware`` already handles::

        {"success": True,  "data": job,   "status": 200}
        {"success": False, "error": "...", "status": 500}

    It must set ``job.status`` to ``QUEUED`` once the job is handed over, and
    must not raise: every failure is a result with ``success: False``.
    ``synchronous=True`` asks it to block until the job finishes; a target
    that cannot (a queue) proceeds asynchronously and says so in its log.
    """

    def run_job(self, job: Any, *, synchronous: bool = False) -> Dict[str, Any]: ...


@runtime_checkable
class ProgramTarget(Protocol):
    """Axis B: one heavy external program runs elsewhere; the job stays here.

    The CCP4i2 job sits in ``RUNNING_REMOTELY`` holding the handle; a
    reconcile step asks ``poll`` and maps the answer onto job status.

    - ``submit`` takes the staged input tree, the program's argv, the output
      directory and the ``sizing_hint`` dict, and returns an opaque string
      handle. The handle is recorded target-tagged, never bare.
    - ``poll`` returns one of ``"queued"``, ``"running"``, ``"succeeded"``,
      ``"failed"``, ``"cancelled"``, ``"unknown"``.
    - ``cancel`` asks the platform to stop the run; it may be a no-op.
    - ``logs`` returns the path, under the job directory, of the program's
      stderr, present by the time ``poll`` reports a terminal state. A target
      never classifies a failure: CCP4i2's own catalogue reads that file.
    """

    def submit(self, tree: Any, argv: list, out_dir: Any, sizing_hint: dict) -> str: ...
    def poll(self, handle: str) -> str: ...
    def cancel(self, handle: str) -> None: ...
    def logs(self, handle: str) -> Any: ...
