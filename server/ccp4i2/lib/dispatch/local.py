"""The one run target CCP4i2 ships: a job as a local subprocess.

Moved here from ``lib/utils/jobs/context_run.py`` unchanged in behaviour; that
module still exports :func:`run_job_local` for its existing callers.
"""
import logging
import os
import pathlib
import shutil
import subprocess
import sys

logger = logging.getLogger(f"ccp4i2:{__name__}")

# ccp4i2/lib/dispatch/local.py -> dispatch -> lib -> ccp4i2 -> scripts/
SCRIPTS_DIR = pathlib.Path(__file__).resolve().parent.parent.parent / "scripts"


def _find_python_interpreter() -> tuple:
    """The interpreter that runs the job.

    1. ccp4-python on PATH (after sourcing ccp4.setup-sh): the CCP4 environment.
    2. sys.executable: the interpreter running this process (has Django).
    """
    ccp4_python = shutil.which("ccp4-python")
    if ccp4_python:
        logger.info("Found ccp4-python on PATH: %s", ccp4_python)
        return ccp4_python, "ccp4-python"
    if sys.executable:
        logger.info("Using current interpreter: %s", sys.executable)
        return sys.executable, "sys.executable"
    return None, None


def run_job_local(job, synchronous=False):
    """Run the job as a detached (or, synchronously, a blocking) subprocess.

    Returns the result dict of :class:`ccp4i2.lib.dispatch.base.JobTarget`;
    never raises. Asynchronous runs go through the crash-safe wrapper script
    so a C-extension crash marks the job FAILED instead of leaving it RUNNING.
    """
    logger.info("Running job %s in LOCAL mode via subprocess (synchronous=%s)",
                job.id, synchronous)
    try:
        python_interpreter, interpreter_name = _find_python_interpreter()
        if python_interpreter is None:
            error_msg = ("No suitable Python interpreter found. "
                         "Either source ccp4.setup-sh to get ccp4-python on PATH, "
                         "or create a virtual environment at .venv")
            logger.error(error_msg)
            return {"success": False, "error": error_msg, "status": 500}

        env = os.environ.copy()

        from ccp4i2.db import models
        job.status = models.Job.Status.QUEUED
        job.save()

        if synchronous:
            logger.info("Running job %s (%s) synchronously using %s",
                        job.id, job.uuid, interpreter_name)
            result = subprocess.run(
                [python_interpreter, "-m", "django", "run_job", "-ju", str(job.uuid)],
                env=env, capture_output=True, text=True)
            job.refresh_from_db()
            if result.returncode != 0:
                logger.warning("Job %s completed with non-zero exit code %d: %s",
                               job.id, result.returncode, result.stderr)
            logger.info("Job %s (%s) completed synchronously with status %s",
                        job.id, job.uuid, job.status)
            return {"success": True, "data": job, "status": 200}

        if sys.platform == "win32":
            subprocess.Popen(
                [str(SCRIPTS_DIR / "run_job_safe.cmd"), python_interpreter, str(job.uuid)],
                creationflags=subprocess.CREATE_NEW_PROCESS_GROUP, env=env)
        else:
            subprocess.Popen(
                ["/bin/bash", str(SCRIPTS_DIR / "run_job_safe.sh"),
                 python_interpreter, str(job.uuid)],
                start_new_session=True, env=env)
        logger.info("Started job %s (%s) via crash-safe wrapper using %s",
                    job.id, job.uuid, interpreter_name)
        return {"success": True, "data": job, "status": 200}

    except Exception as error:  # noqa: BLE001 -- contract: never raise
        logger.exception("Failed to start job via subprocess", exc_info=error)
        return {"success": False, "error": f"Subprocess error: {str(error)}", "status": 500}


class LocalTarget:
    """Registered as ``local``; the default everywhere."""

    name = "local"

    def run_job(self, job, *, synchronous=False):
        # Through the module attribute, so a test can patch run_job_local.
        return run_job_local(job, synchronous=synchronous)
