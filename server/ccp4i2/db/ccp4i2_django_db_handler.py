import datetime
import logging
import sys
import uuid

from django.utils import timezone

from ccp4i2.core.CCP4PluginScript import CPluginScript

from . import models
from .ccp4i2_django_dbapi import CCP4i2DjangoDbApi
from .ccp4i2_static_data import (
    JOB_STATUS_FAILED,
    JOB_STATUS_FINISHED,
    JOB_STATUS_INTERRUPTED,
    JOB_STATUS_TO_DELETE,
    JOB_STATUS_UNSATISFACTORY,
)


logger = logging.getLogger(f"ccp4x:{__name__}")
logger.setLevel(logging.WARNING)


def plugin_status_to_job_status(finishStatus):
    status = JOB_STATUS_FAILED
    if isinstance(finishStatus, dict):
        finishStatus = finishStatus.get("finishStatus")
    if finishStatus == CPluginScript.SUCCEEDED:
        status = JOB_STATUS_FINISHED
    elif finishStatus == CPluginScript.FAILED:
        status = JOB_STATUS_FAILED
    elif finishStatus == CPluginScript.INTERRUPTED:
        status = JOB_STATUS_INTERRUPTED
    elif finishStatus == CPluginScript.MARK_TO_DELETE:
        status = JOB_STATUS_TO_DELETE
    elif finishStatus == CPluginScript.UNSATISFACTORY:
        status = JOB_STATUS_UNSATISFACTORY
    return status


class CCP4i2DjangoDbHandler:
    """
    Handler class for interacting with the CCP4i2 Django database.

    Methods
    -------
    __init__():
        Initializes the database handler with a CCP4i2DjangoDbApi instance.

    createJob(pluginName, jobTitle=None, parentJobId=None, jobNumber=None):
        Creates a new job in the database.

        Parameters
        ----------
        pluginName : str
            The name of the plugin associated with the job.
        jobTitle : str, optional
            The title of the job (default is None).
        parentJobId : str, optional
            The ID of the parent job (default is None).
        jobNumber : int, optional
            The job number (default is None).

        Returns
        -------
        Job
            The created job object.

        Raises
        ------
        Exception
            If there is an error during job creation.

    updateJobStatus(jobId=None, status=None, finishStatus=None, container=None, dbOutputData=None):
        Updates the status of an existing job in the database.

        Parameters
        ----------
        jobId : str
            The ID of the job to update.
        status : str, optional
            The new status of the job (default is None).
        finishStatus : str, optional
            The finish status of the job (default is None).
        container : object, optional
            The container associated with the job (default is None).
        dbOutputData : object, optional
            Additional data output from the database (default is None).

        Returns
        -------
        int
            A status code indicating success (CPluginScript.SUCCEEDED).

        Raises
        ------
        AssertionError
            If jobId is not a string.
        Exception
            If there is an error during job status update.
    """

    def __init__(self):
        self.db = CCP4i2DjangoDbApi()

    def createJob(self, pluginName, jobTitle=None, parentJobId=None, jobNumber=None):
        logger.info(
            "Creating job %s %s %s %s", pluginName, jobTitle, parentJobId, jobNumber
        )
        try:
            # Deferred import to avoid circular dependency
            from ..lib.utils.jobs.create import create_job

            return create_job(
                parentJobId=uuid.UUID(parentJobId) if parentJobId else None,
                taskName=pluginName,
                jobNumber=jobNumber,
                saveParams=False,
                title=jobTitle,
            )
        except Exception as err:
            logger.error("Failed in createJob %s", err)
            raise (err)

    def createSubJob(self, taskName, parentJobId, jobNumber, jobTitle=None):
        """
        Create a database job entry for a sub-job.

        This is called by CPluginScript.makePluginObject() when creating nested jobs
        in a database-backed context.

        Args:
            taskName: Name of the plugin/task for the sub-job
            parentJobId: UUID of the parent job (as string)
            jobNumber: Job number for this sub-job (e.g., "1", "2")
            jobTitle: Optional title for the sub-job

        Returns:
            str: UUID of the created job

        Raises:
            Exception: If job creation fails
        """
        logger.debug(
            f"[createSubJob] Creating sub-job: task={taskName}, parent={parentJobId}, number={jobNumber}"
        )

        try:
            # Create the job using existing createJob method
            job = self.createJob(
                pluginName=taskName,
                jobTitle=jobTitle or f"{taskName} sub-job {jobNumber}",
                parentJobId=parentJobId,
                jobNumber=jobNumber
            )

            # Extract job UUID from returned job object
            if hasattr(job, 'jobId'):
                job_id = str(job.jobId)
            elif hasattr(job, 'job_id'):
                job_id = str(job.job_id)
            elif hasattr(job, 'uuid'):
                job_id = str(job.uuid)
            elif isinstance(job, dict):
                job_id = str(job.get('jobId') or job.get('job_id') or job.get('uuid'))
            else:
                # Fallback: assume job is the UUID itself
                job_id = str(job)

            logger.debug(f"[createSubJob] ✅ Created sub-job with ID: {job_id}")
            return job_id

        except Exception as e:
            logger.error(f"[createSubJob] Failed to create sub-job: {e}")
            raise

    def updateJobStatus(
        self,
        jobId=None,
        status=None,
        finishStatus=None,
        container=None,
        dbOutputData=None,
    ):
        try:
            assert isinstance(jobId, str)
        except AssertionError as err:
            logger.exception("In updateJobStatus %s" % (type(jobId),), exc_info=err)
            jobId = str(jobId)

        logger.debug(
            "In update JobStatus %s %s %s %s %s",
            jobId,
            status,
            finishStatus,
            container,
            dbOutputData,
        )
        sys.stdout.flush()
        if dbOutputData is not None:
            logger.error("dbOutputData is not None %s", dbOutputData)

        job_uuid = uuid.UUID(jobId)

        try:
            if status is None and finishStatus is not None:
                status = plugin_status_to_job_status(finishStatus)
            try:
                the_job = models.Job.objects.get(uuid=job_uuid)
                # Only update status if explicitly provided (not None)
                if status is not None:
                    the_job.status = status
                    the_job.save()
                if status is not None and models.Job.Status(status).label == "Finished":
                    the_job.finish_time = timezone.now()
                    the_job.save()
                    self.db.gleanJobFiles(container=container, jobId=jobId)
                if the_job.parent is None and models.Job.Status(
                    the_job.status
                ).label in [
                    "Finished",
                    "Interrupted",
                    "Failed",
                    "Unsatisfactory",
                ]:
                    pass  # backupProjectDb(the_job.projectid)
            except Exception as err:
                logger.exception(
                    "Failed in updateJobStatus %s" % (the_job,),
                    exc_info=err,
                )
        except Exception as err:
            logger.error("Issue in reportStatus %s %s", err, the_job, exc_info=err)
        return CPluginScript.SUCCEEDED

    def getProjectDirectory(self, projectId):
        """
        Get the directory path for a project.

        Args:
            projectId: Project UUID (string, with or without hyphens)

        Returns:
            str: Absolute path to the project directory, or None if not found
        """
        try:
            # Normalize projectId to UUID
            if isinstance(projectId, str):
                # Remove hyphens if present, then add them back in standard format
                clean_id = projectId.replace('-', '')
                if len(clean_id) == 32:
                    projectId = f"{clean_id[0:8]}-{clean_id[8:12]}-{clean_id[12:16]}-{clean_id[16:20]}-{clean_id[20:]}"

            project_uuid = uuid.UUID(str(projectId))
            project = models.Project.objects.get(uuid=project_uuid)
            return str(project.directory)
        except Exception as err:
            logger.debug(f"Could not get project directory for {projectId}: {err}")
            return None
