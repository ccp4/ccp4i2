"""For plugin methods that change a job's parameters (reached through the
object_method endpoint, which neither checks a job's status nor saves): the
job's row, a refusal unless it is still editable, and the save."""
import xml.etree.ElementTree as ET
from pathlib import Path

from ccp4i2.db import models

EDITABLE = (models.Job.Status.UNKNOWN, models.Job.Status.PENDING)


def job_row(plugin):
    """The plugin's job: from its database context, else from the jobId its
    parameters' header records; None if it is not a job in the database."""
    job_id = plugin.get_db_job_id() if hasattr(plugin, "get_db_job_id") else None
    if not job_id:
        job_id = getattr(plugin, "_dbJobId", None)
    if not job_id:
        for name in ("input_params.xml", "params.xml"):
            path = Path(str(plugin.workDirectory)) / name
            if path.is_file():
                try:
                    job_id = ET.parse(path).getroot().findtext("ccp4i2_header/jobId")
                except ET.ParseError:
                    job_id = None
                if job_id:
                    break
    if not job_id:
        return None
    return models.Job.objects.filter(uuid=job_id).first()


def editable_job(plugin):
    """(job, None) when the plugin's job may be changed, else (None, why)."""
    job = job_row(plugin)
    if job is None:
        return None, "this job is not in the database"
    if job.status not in EDITABLE:
        return None, (f"job {job.number} is {job.get_status_display().lower()}; "
                      "only a pending job can be changed")
    return job, None


def save(plugin, job):
    from ccp4i2.lib.utils.parameters.save_params import save_params_for_job
    save_params_for_job(plugin, job)
