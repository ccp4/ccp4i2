"""A job whose program runs on a program run target: the record and the reconcile.

Generic (decision 18 of docs/pandda-campaign-design.md): any task that hands
one heavy program to a registered program target (``ccp4i2.lib.dispatch``)
works the same way.

1. At submit, the plugin writes ``dispatch.json`` in its job directory --
   target name, target-tagged handle, when -- and returns
   ``CPluginScript.DISPATCHED``; the runner leaves the job in
   ``RUNNING_REMOTELY``.
2. :func:`reconcile` is idempotent and caller-agnostic: it reads the record,
   asks the target (``poll``), and maps the answer onto the job. Not
   terminal: nothing changes. Terminal: it records the state and where the
   program's stderr is (``logs``), then starts the job again through the
   deployment's job target. The plugin, seeing a terminal record, skips the
   submit and completes from the output the run left -- so postProcess,
   gleaning and the report run in a job process exactly as for a local run.
   Who calls the reconcile is a deployment choice: the job page, a
   management command, the fan-out, a cron.
3. A target never classifies a failure. It says where stderr is; the plugin's
   own catalogue reads it.
"""
import json
import logging
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path

logger = logging.getLogger(f"ccp4i2:{__name__}")

RECORD_NAME = "dispatch.json"
TERMINAL_STATES = frozenset({"succeeded", "failed", "cancelled"})
KNOWN_STATES = TERMINAL_STATES | {"submitted", "queued", "running", "unknown"}


def _now():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def record_path(job_dir) -> Path:
    return Path(job_dir) / RECORD_NAME


def read_record(job_dir):
    """The dispatch record, or None when this job never dispatched."""
    path = record_path(job_dir)
    if not path.is_file():
        return None
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError) as err:
        logger.warning("dispatch record %s unreadable: %s", path, err)
        return None


def write_record(job_dir, **fields) -> dict:
    """Merge ``fields`` into the record and write it atomically."""
    record = read_record(job_dir) or {}
    record.update(fields)
    path = record_path(job_dir)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=".dispatch-", suffix=".json", dir=str(path.parent))
    with os.fdopen(fd, "w") as fh:
        json.dump(record, fh, indent=2, sort_keys=True)
    os.replace(tmp, path)
    return record


def new_record(job_dir, target: str, handle: str, **extra) -> dict:
    """The record a plugin writes at submit."""
    return write_record(job_dir, target=target, handle=str(handle), state="submitted",
                        submitted_at=_now(), **extra)


def is_harvest(job_dir) -> bool:
    """True when the run has finished elsewhere and this job process is here
    to complete from it (the plugin's startProcess asks this)."""
    record = read_record(job_dir)
    return bool(record) and record.get("state") in TERMINAL_STATES


def _default_run(job):
    from .context_run import run_job_context_aware
    return run_job_context_aware(job, force_dispatch=True)


def reconcile(job, *, run=None) -> dict:
    """Ask the target how the job's program is doing and act on the answer.

    Returns a dict with ``action`` (``none`` / ``harvest_started`` /
    ``harvest_failed`` / ``error``), the ``state`` the target reported, and a
    ``reason`` in words. Never raises for a target's sake: a target that
    cannot be loaded or polled is an ``error`` result naming why.
    """
    from ccp4i2.db import models
    from ccp4i2.lib.dispatch import (
        RunTargetError, UnknownRunTarget, get_target, runs_programs,
    )

    job_dir = Path(job.directory)
    record = read_record(job_dir)
    if not record:
        return {"action": "none", "reason": "the job has no dispatch record", "state": None}
    if job.status != models.Job.Status.RUNNING_REMOTELY:
        return {"action": "none", "state": record.get("state"),
                "reason": f"the job is {job.get_status_display().lower()}, not running remotely"}
    if record.get("harvest_started_at"):
        return {"action": "none", "state": record.get("state"),
                "reason": f"harvest already started at {record['harvest_started_at']}"}

    try:
        target = get_target(record["target"])
    except (KeyError, UnknownRunTarget, RunTargetError) as err:
        return {"action": "error", "state": record.get("state"), "reason": str(err)}
    if not runs_programs(target):
        return {"action": "error", "state": record.get("state"),
                "reason": f"run target '{record['target']}' does not run programs"}

    try:
        state = str(target.poll(record["handle"])).lower()
    except Exception as err:  # noqa: BLE001 -- the target's failure is the result
        return {"action": "error", "state": record.get("state"),
                "reason": f"poll failed: {type(err).__name__}: {err}"}
    if state not in KNOWN_STATES:
        logger.warning("run target %s reported an unknown state %r; treating as unknown",
                       record["target"], state)
        state = "unknown"
    record = write_record(job_dir, state=state, polled_at=_now())
    if state not in TERMINAL_STATES:
        return {"action": "none", "state": state, "reason": f"the run is {state}"}

    stderr = None
    try:
        stderr = target.logs(record["handle"])
    except Exception as err:  # noqa: BLE001 -- a missing log is not a reason to stop
        logger.warning("run target %s could not say where stderr is: %s", record["target"], err)
    record = write_record(job_dir, stderr=str(stderr) if stderr else None, harvest_started_at=_now())

    # Start the job again; its plugin completes from the record.
    job.status = models.Job.Status.PENDING
    job.save(update_fields=["status"])
    result = (run or _default_run)(job)
    if not result.get("success"):
        write_record(job_dir, harvest_started_at=None, harvest_error=result.get("error"))
        job.status = models.Job.Status.RUNNING_REMOTELY
        job.save(update_fields=["status"])
        return {"action": "harvest_failed", "state": state,
                "reason": result.get("error", "the job could not be started")}
    return {"action": "harvest_started", "state": state,
            "reason": f"the run {state}; the job is completing from it"}


def reconcile_all(*, run=None):
    """Every job running remotely, reconciled; returns [(job, result)]."""
    from ccp4i2.db import models
    out = []
    for job in models.Job.objects.filter(status=models.Job.Status.RUNNING_REMOTELY).order_by("id"):
        out.append((job, reconcile(job, run=run)))
    return out
