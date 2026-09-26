"""
Context-dependent job execution: resolve the run target, then run.

Where a job runs is a *run target* (ccp4i2.lib.dispatch): CCP4i2 ships
``local`` (a subprocess) and a deployment registers its own in settings
(``CCP4I2_RUN_TARGETS`` / ``CCP4I2_JOB_TARGET``). Nothing here names a
platform. With no settings at all, jobs run locally.

Example Usage:
    from ccp4i2.lib.utils.jobs.context_run import run_job_context_aware
    result = run_job_context_aware(job)
    if result["success"]:
        return Response(result["data"])      # started / queued
    else:
        return Response({"error": result["error"]}, status=result["status"])
"""

import os
import logging
import functools
import shutil

from ccp4i2.lib.dispatch import (
    UnknownRunTarget, RunTargetError, get_target, job_target_name, runs_jobs,
)
# Kept as a name here for existing importers; the implementation lives with
# the other run targets.
from ccp4i2.lib.dispatch.local import run_job_local  # noqa: F401

logger = logging.getLogger(__name__)


@functools.lru_cache(maxsize=1)
def ccp4_available():
    """True if this server has a usable CCP4 installation (can run CCP4 binaries).

    A CCP4-free deployment (e.g. the slim Django API container) returns False, so
    callers can decide not to attempt local execution of tasks that need a binary.
    Cached: the environment is fixed for the life of the server process.
    """
    ccp4 = os.environ.get("CCP4")
    if not ccp4 or not os.path.isdir(ccp4):
        return False
    # Probe a core CCP4 program ('cad' is present in any real install).
    if shutil.which("cad"):
        return True
    return os.path.exists(os.path.join(ccp4, "bin", "cad"))


def task_requires_ccp4(task_name):
    """True unless the task is declared ccp4_free in the registry.

    Conservative: an unknown or unflagged task is assumed to need CCP4.
    """
    try:
        from ccp4i2.core.tasks import TASKS
        task = TASKS.get(task_name)
    except Exception:
        task = None
    return not bool(task and getattr(task, "ccp4_free", False))


def can_run_local(task_name):
    """Whether `task_name` can be executed locally in this environment.

    True if either the server has CCP4, or the task needs no CCP4 (ccp4_free).
    This is the server-side feasibility decision behind the run_local endpoint:
    the client may *request* local execution, but the server reports whether it
    is actually possible here.
    """
    return ccp4_available() or not task_requires_ccp4(task_name)


def program_checks_are_authoritative():
    """True when "binary not found here" really means the job will fail.

    A pre-run program-availability check is only trustworthy when this process
    is the one that will spawn the job, and it has a CCP4 installation to look
    in. Three deployments, three answers:

    * desktop / Electron / i2run — local mode with CCP4 mounted. The check is
      authoritative, so a missing binary should *block* submission: the job is
      certain to fail, and failing now with "shelxe was not found; set its
      location in Preferences" beats failing later with a silent empty result.
    * a remote job target (e.g. a queue to a worker) — the job runs on a
      filesystem we cannot see. A "not found" here says nothing about it.
    * the slim CCP4-free API server — there is nothing to look in at all.

    In the last two the absence is our ignorance, not a defect in the user's
    setup, so the check stays advisory and must not block Confirm.
    """
    return job_target_name() == "local" and ccp4_available()


def run_job_context_aware(job, force_local=False, synchronous=False,
                          force_dispatch=False):
    """
    Run a job on this deployment's job target.

    Resolves the target by name (``CCP4I2_JOB_TARGET``, default ``local``)
    through the run-target registry and hands the job to it. An interactive
    task opens its session instead, unless ``force_dispatch``.

    Args:
        job: Job model instance
        force_local (bool): run on the ``local`` target whatever the deployment's
            job target is (the run_local endpoint; the caller has checked
            feasibility with ``can_run_local``).
        synchronous (bool): block until the job completes. Honoured by the local
            target; a target that queues proceeds asynchronously and logs so.
        force_dispatch (bool): dispatch even an interactive task (used when its
            session is finished); otherwise Run opens the session instead.

    Returns:
        dict: ``{"success": True, "data": job, "status": 200}`` or
        ``{"success": False, "error": str, "status": int}``. Never raises.
    """
    # An interactive task (the recorded Moorhen session) has no process to
    # dispatch on Run: Run opens the session and the job is dispatched when
    # the session is finished (force_dispatch=True from finish_session).
    if not force_dispatch:
        from .interactive import SessionError, is_interactive_job, open_session

        if is_interactive_job(job):
            try:
                return {"success": True, "data": open_session(job)}
            except SessionError as err:
                return {"success": False, "error": str(err), "status": err.status}

    name = "local" if force_local else job_target_name()
    try:
        target = get_target(name)
    except (UnknownRunTarget, RunTargetError) as err:
        logger.error("Job %s not run: %s", job.id, err)
        return {"success": False, "error": str(err), "status": 500}
    if not runs_jobs(target):
        msg = (f"run target '{name}' does not run jobs (no run_job); "
               "set CCP4I2_JOB_TARGET to one that does")
        logger.error("Job %s not run: %s", job.id, msg)
        return {"success": False, "error": msg, "status": 500}

    logger.info("Executing job %s (uuid=%s, task=%s) on target '%s' (synchronous=%s%s)",
                job.id, job.uuid, job.task_name, name, synchronous,
                ", forced local" if force_local else "")
    return target.run_job(job, synchronous=synchronous)
