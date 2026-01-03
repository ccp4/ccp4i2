import logging

from ccp4i2.core import CCP4TaskManager

from ccp4i2.db.models import Job
from ccp4i2.db.async_db_handler import AsyncDatabaseHandler

logger = logging.getLogger(f"ccp4i2:{__name__}")


def get_job_plugin(the_job: Job, parent=None, dbHandler=None):
    """
    Retrieves and initializes a job plugin instance based on the provided job.

    Args:
        the_job (Job): The job object containing details such as task name and directory.
        parent (optional): The parent object, if any. Defaults to None.
        dbHandler (optional): The database handler for CCP4i2. If not provided, an
            AsyncDatabaseHandler will be created automatically using the job's project UUID.
            This ensures CDataFile objects can resolve paths via dbFileId.

    Returns:
        pluginInstance: An instance of the plugin class initialized with the job's parameters.
        None: If an error occurs during plugin initialization.

    Raises:
        Exception: If no parameter definition file (params.xml or input_params.xml) is found in the job directory.
    """

    taskManager = CCP4TaskManager.CTaskManager()

    # Create a default database handler if one wasn't provided
    # This ensures CDataFile objects can resolve paths via dbFileId
    if dbHandler is None:
        dbHandler = AsyncDatabaseHandler(project_uuid=str(the_job.project.uuid))

    pluginClass = taskManager.get_plugin_class(the_job.task_name)
    try:
        pluginInstance = pluginClass(
            workDirectory=str(the_job.directory), parent=parent, dbHandler=dbHandler
        )
    except Exception as err:
        logger.exception("Error in get_job_plugin was", exc_info=err)
        return None

    # Determine which params file to load based on job status
    # For UNKNOWN/PENDING jobs, prefer input_params.xml (user control stage)
    # For other jobs, prefer params.xml (plugin lifecycle output)
    params_path = the_job.directory / "params.xml"
    fallback_params_path = the_job.directory / "input_params.xml"
    if the_job.status in [Job.Status.UNKNOWN, Job.Status.PENDING]:
        params_path = the_job.directory / "input_params.xml"
        fallback_params_path = the_job.directory / "params.xml"

    params_file = params_path
    if not params_file.exists():
        params_file = fallback_params_path

    # Load params if file exists, otherwise return fresh plugin from .def.xml
    # This is critical for new jobs - loading from fresh .def.xml gives proper
    # CData wrappers. If we tried to load from an empty XML, we'd get plain types.
    if params_file.exists():
        # Use CPluginScript.loadDataFromXml() which handles ParamsXmlHandler format
        # (with <ccp4i2> wrapper and header) as well as legacy CContainer format
        error = pluginInstance.loadDataFromXml(str(params_file))
        if error and hasattr(error, 'hasError') and error.hasError():
            logger.error("Failed to load params from %s: %s", params_file, error)
            raise Exception(f"Failed to load params from {params_file}: {error}")
        logger.info(f"Loaded parameters from {params_file}")
    else:
        logger.info(f"No params file found - using fresh plugin from .def.xml")

    return pluginInstance
