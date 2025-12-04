import logging

from core import CCP4TaskManager
from core import CCP4Container
from ....db import models

logger = logging.getLogger(f"ccp4x:{__name__}")


def get_job_container(the_job: models.Job):
    """
    Retrieves and loads a job container for the given job.

    This function looks up the definition file for the specified job task,
    creates a container, and loads its contents from the definition file.
    It then attempts to load additional data from either 'params.xml' or
    'input_params.xml' located in the job's directory.

    Args:
        the_job (Job): The job object containing task information and directory paths.

    Returns:
        CCP4Container.CContainer: The loaded job container.
    """
    defFile = CCP4TaskManager.CTaskManager().locate_def_xml(
        task_name=the_job.task_name, version=None
    )
    # print 'CProjectDirToDb.globJobs defFile',defFile
    container = CCP4Container.CContainer()
    container.loadContentsFromXml(defFile)

    params_path = the_job.directory / "params.xml"
    fallback_params_path = the_job.directory / "input_params.xml"
    if the_job.status in [models.Job.Status.UNKNOWN, models.Job.Status.PENDING]:
        params_path = the_job.directory / "input_params.xml"
        fallback_params_path = the_job.directory / "params.xml"
    if (params_path).exists():
        container.loadDataFromXml(str(params_path))
    else:
        container.loadDataFromXml(str(fallback_params_path))
    return container
