from pathlib import Path
import logging
import uuid

from core import CCP4TaskManager
from core.CCP4Container import CContainer

from ccp4x.db import models
from ..containers.remove_defaults import remove_container_default_values
from ..parameters.save_params import save_params_for_job
from ..files.patch_paths import patch_output_file_paths
from ..plugins.get_plugin import get_job_plugin


logger = logging.getLogger(f"ccp4x:{__name__}")


def create_job(
    projectId: str = None,
    projectName: str = None,
    parentJobId: str = None,
    taskName: str = None,
    jobNumber: str = None,
    jobId: str = None,
    saveParams: bool = True,
    title: str = None,
):
    """
    Create a new job within a project.

    Args:
        projectId (str, optional): The UUID of the project. Defaults to None.
        projectName (str, optional): The name of the project. Defaults to None.
        parentJobId (str, optional): The UUID of the parent job. Defaults to None.
        taskName (str, optional): The name of the task to be executed. Defaults to None.
        jobNumber (str, optional): The job number. Defaults to None.
        jobId (str, optional): The UUID of the job. Defaults to None.
        saveParams (bool, optional): Whether to save parameters for the job. Defaults to True.
        title (str, optional): The title of the job. Defaults to None.

    Returns:
        str: The UUID of the newly created job.

    Raises:
        models.Job.DoesNotExist: If the parent job or project does not exist.
    """
    logger.debug("%s, %s", projectName, projectId)
    if parentJobId is not None and projectId is None:
        parent_job = models.Job.objects.get(uuid=parentJobId)
        the_project = parent_job.project
    elif projectId is None and projectName is not None:
        parent_job = None
        the_project = models.Project.objects.get(name=projectName)
        projectId = the_project.uuid
    else:
        parent_job = None
        the_project = models.Project.objects.get(uuid=projectId)

    if jobNumber is None:
        project_jobs = models.Job.objects.filter(project__uuid=projectId).filter(
            parent__isnull=True
        )
        if len(project_jobs) == 0:
            last_job_number = 0
        else:
            last_job_number = sorted([int(a.number) for a in project_jobs])[-1]
        last_job_number = str(last_job_number)
    else:
        job_number_elements = jobNumber.split(".")
        job_number_elements[-1] = str(int(job_number_elements[-1]) - 1)
        last_job_number = ".".join(job_number_elements)

    job_number_elements = last_job_number.split(".")
    job_number_elements[-1] = str(int(job_number_elements[-1]) + 1)
    next_job_number = ".".join(job_number_elements)

    new_job_dir = Path(the_project.directory).joinpath(
        *(["CCP4_JOBS"] + [f"job_{j_no}" for j_no in next_job_number.split(".")])
    )

    if jobId is None:
        new_job_id = uuid.uuid4()
    elif "-" not in jobId:
        new_job_id = uuid.UUID(jobId)
    else:
        new_job_id = jobId

    task_manager = CCP4TaskManager.CTaskManager()
    plugin_class = task_manager.get_plugin_class(taskName)
    if saveParams:
        new_job_dir.mkdir(exist_ok=True, parents=True)
    new_job_plugin = plugin_class(workDirectory=str(new_job_dir))

    if title is None:
        title = task_manager.getTitle(taskName)
        # Fallback to taskName if getTitle returns None
        if title is None:
            title = taskName
    arg_dict = dict(
        uuid=new_job_id,
        number=str(next_job_number),
        status=models.Job.Status.PENDING,
        evaluation=models.Job.Evaluation.UNKNOWN,
        title=title,
        project=the_project,
        task_name=taskName,
        parent=parent_job,
    )
    logger.info("arg_dict %s", arg_dict)
    new_job = models.Job(**arg_dict)
    if saveParams:
        remove_container_default_values(new_job_plugin.container)
        patch_output_file_paths(new_job_plugin, new_job)
        save_params_for_job(new_job_plugin, new_job, exclude_unset=True)
    new_job.save()

    return str(new_job.uuid)
