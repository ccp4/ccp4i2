from ccp4x.db import models
import logging
import uuid
from pytz import timezone
from ccp4i2.core import CCP4TaskManager

logger = logging.getLogger(f"ccp4x:{__name__}")


def job_directory(jobId=None, projectName=None, jobNumber=None, create=False, projectId=None, projectDirectory=None ):
    the_job = None
    if jobId is not None:
        the_job = models.Job.objects.get(uuid=uuid.UUID(jobId))
    elif projectName is not None and jobNumber is not None:
        the_job = models.Job.objects.get(project__name=projectName, number=jobNumber)
    elif projectId is not None and jobNumber is not None:
        the_job = models.Job.objects.get(project__uuid=projectId, number=jobNumber)
    if the_job is None:
        logger.error(
            "No job found with jobId %s, projectName %s, jobNumber %s",
            jobId,
            projectName,
            jobNumber,
        )
        return None
    if the_job.directory is None:
        logger.error(
            "Job %s has no directory set, cannot return job directory",
            the_job.uuid,
        )
        return None
    logger.debug("Job directory is %s %s", the_job.directory, the_job.project)
    if create:
        import os

        os.makedirs(the_job.directory, exist_ok=True)
    return str(the_job.directory)
