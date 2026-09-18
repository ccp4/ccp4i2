"""Server-side half of the bridge handshake (Python 3, wrapper context).

One function, shared by every Coot task wrapper: publish connection
details and the job identity into ``os.environ`` so the module loaded
into Coot (which receives nothing else) can find and authenticate to
the API and identify its job.

Why mutating ``os.environ`` works: the plugin's subprocess environment
is a copy of this process's environment
(CPluginScript._prepareProcessExecution), and this job-runner process
inherited the server's environment - including, in Electron mode,
CCP4I2_LOCAL_SESSION_TOKEN and UVICORN_PORT.
"""

import os
from pathlib import Path

#: Directory holding this package - and the static Coot startup stubs.
COOTBRIDGE_DIR = os.path.dirname(os.path.abspath(__file__))


def export_handshake(plugin, work_dir: Path, drop_dir: Path) -> None:
    """Set the CCP4I2_* contract variables for a Coot task's subprocess.

    ``plugin`` is the CPluginScript instance (source of the db identity
    set by setDbData: _dbJobId/_dbProjectId are uuids).
    """
    env = os.environ

    # Where the in-Coot modules live: the py3 stub imports the package,
    # the py2 (0.9) stub imp.load_source's these files by path.
    env["CCP4I2_COOTBRIDGE_DIR"] = COOTBRIDGE_DIR

    if not env.get("CCP4I2_API_URL"):
        port = env.get("UVICORN_PORT")
        if port:
            env["CCP4I2_API_URL"] = f"http://127.0.0.1:{port}"
    token = env.get("CCP4I2_LOCAL_SESSION_TOKEN")
    if token and not env.get("CCP4I2_ACCESS_TOKEN"):
        env["CCP4I2_ACCESS_TOKEN"] = token

    job_uuid = getattr(plugin, "_dbJobId", None)
    if job_uuid:
        env["CCP4I2_JOB_UUID"] = str(job_uuid)
    if getattr(plugin, "_dbProjectId", None):
        env["CCP4I2_PROJECT_UUID"] = str(plugin._dbProjectId)
    if getattr(plugin, "_dbProjectName", None):
        env["CCP4I2_PROJECT_NAME"] = str(plugin._dbProjectName)
    env["CCP4I2_JOB_DIRECTORY"] = str(work_dir)
    env["CCP4I2_DROP_DIR"] = str(drop_dir)

    project_dir = _project_directory(work_dir)
    if project_dir:
        env["CCP4I2_PROJECT_DIRECTORY"] = str(project_dir)

    # The REST job endpoints key on integer ids; resolve them from the
    # uuids while we still have Django. Guarded: i2run and tests may run
    # without a database row for this job.
    if job_uuid:
        try:
            from ccp4i2.db import models

            job = models.Job.objects.get(uuid=job_uuid)
            env["CCP4I2_JOB_ID"] = str(job.id)
            env["CCP4I2_PROJECT_ID"] = str(job.project.id)
        except Exception:
            pass


def _project_directory(work_dir: Path):
    """The project root is the parent of the CCP4_JOBS path component."""
    parts = work_dir.parts
    if "CCP4_JOBS" in parts:
        return Path(*parts[: parts.index("CCP4_JOBS")])
    return None
