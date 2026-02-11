import logging
from pathlib import Path
import datetime
from pytz import timezone
from ccp4i2.core import CCP4TaskManager
from ccp4i2.core.CCP4Container import CContainer
from ccp4i2.core.task_manager.metadata import TITLES
from ccp4i2.db import models
from ccp4i2.lib.utils.parameters.save_params import save_params_for_job
from ccp4i2.lib.utils.files.patch_paths import patch_output_file_paths
from ccp4i2.lib.response import Result

logger = logging.getLogger(f"ccp4i2:{__name__}")


def clone_job(jobId: str = None) -> Result[models.Job]:
    """
    Clone an existing job by creating a new job with the same parameters.

    Args:
        jobId (str, optional): The UUID of the job to be cloned. Defaults to None.

    Returns:
        Result[models.Job]: Result containing the newly created job instance or error

    Example:
        >>> result = clone_job("550e8400-e29b-41d4-a716-446655440000")
        >>> if result.success:
        ...     new_job = result.data
        ...     print(f"Cloned job {new_job.number}")
        ... else:
        ...     print(f"Clone failed: {result.error}")
    """
    try:
        # Get the job to clone
        try:
            old_job = models.Job.objects.get(uuid=jobId)
        except models.Job.DoesNotExist:
            return Result.fail(
                f"Job not found: {jobId}",
                details={"job_id": jobId}
            )

        the_project = old_job.project
        task_name = old_job.task_name

        # Calculate next job number
        try:
            last_job_number = max(
                int(job.number)
                for job in models.Job.objects.filter(project=the_project).filter(
                    parent__isnull=True
                )
            )
        except ValueError:
            last_job_number = 0

        next_job_number = str(last_job_number + 1)
        new_job_dir = Path(the_project.directory) / "CCP4_JOBS" / f"job_{next_job_number}"

        # Create plugin instance
        task_manager = CCP4TaskManager.CTaskManager()
        plugin_class = task_manager.get_plugin_class(task_name)
        new_job_dir.mkdir(exist_ok=True, parents=True)
        the_job_plugin: CContainer = plugin_class(workDirectory=str(new_job_dir))

        # Load cloned parameters using ParamsXmlHandler
        # This properly overlays parameters from input_params.xml onto the .def.xml defaults
        params_file = old_job.directory / "input_params.xml"

        if not params_file.exists():
            return Result.fail(
                f"Source job parameters not found: {params_file}",
                details={"params_file": str(params_file), "source_job_id": str(jobId)}
            )

        # Use the plugin's loadDataFromXml which uses ParamsXmlHandler for proper overlay
        error = the_job_plugin.loadDataFromXml(str(params_file))
        if error and hasattr(error, 'hasError') and error.hasError():
            return Result.fail(
                f"Failed to load parameters: {error}",
                details={"params_file": str(params_file), "error": str(error)}
            )

        # Create new job record
        # Use title from old job, or get from task manager, or fall back to task name
        title = old_job.title or TITLES.get(task_name, task_name)
        new_job = models.Job(
            number=str(next_job_number),
            finish_time=datetime.datetime.fromtimestamp(0, tz=timezone("UTC")),
            status=models.Job.Status.PENDING,
            evaluation=models.Job.Evaluation.UNKNOWN,
            title=title,
            project=the_project,
            task_name=task_name,
            parent=None,
        )
        new_job.save()

        # Update file paths and save parameters
        patch_output_file_paths(the_job_plugin, new_job)
        save_params_for_job(the_job_plugin, new_job)

        # Update project metadata
        the_project.last_access = datetime.datetime.now(tz=timezone("UTC"))
        the_project.last_job_number = new_job.number
        the_project.save()

        logger.info(
            "Cloned job %s to new job %s in project %s",
            jobId, new_job.uuid, the_project.name
        )

        return Result.ok(new_job)

    except Exception as err:
        logger.exception("Failed to clone job %s", jobId, exc_info=err)
        return Result.fail(
            f"Failed to clone job: {str(err)}",
            details={"job_id": jobId, "error_type": type(err).__name__}
        )
