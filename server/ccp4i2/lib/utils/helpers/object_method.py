from typing import List
from ccp4i2.db import models
from ..plugins.get_plugin import get_job_plugin


def object_method(
    the_job: models.Job,
    object_path: str,
    method_name: str,
    args: List[str] = None,
    kwargs: dict = None,
):
    """
    Execute a method on a CData object within a job's container.

    Args:
        the_job: The Job model instance
        object_path: Dot-separated path to the object (e.g., "task.inputData.ASU_CONTENT")
        method_name: Name of the method to call on the object
        args: Positional arguments to pass to the method
        kwargs: Keyword arguments to pass to the method

    Returns:
        The return value of the called method
    """
    if args is None:
        args = []
    if kwargs is None:
        kwargs = {}

    the_job_plugin = get_job_plugin(the_job)

    # get_job_plugin loads the parameters but does not tell the plugin which
    # job it is. A method that saves them back (fillDatasetsFromCampaign does)
    # then writes a header with no jobId, and every later call that resolves
    # the job from that header answers "this job is not in the database".
    # Stamp the identity the same way save_params does.
    if getattr(the_job_plugin, "_dbJobId", None) is None:
        the_job_plugin._dbJobId = str(the_job.uuid)
    if getattr(the_job_plugin, "_dbJobNumber", None) is None:
        the_job_plugin._dbJobNumber = the_job.number
    if getattr(the_job_plugin, "_dbProjectId", None) is None:
        the_job_plugin._dbProjectId = str(the_job.project.uuid)
    if getattr(the_job_plugin, "_dbProjectName", None) is None:
        the_job_plugin._dbProjectName = the_job.project.name

    # If the path is just the task name (no dots beyond the first segment),
    # call the method on the plugin itself rather than a container element.
    path_segments = object_path.split(".")
    if len(path_segments) <= 1 and hasattr(the_job_plugin, method_name):
        result = getattr(the_job_plugin, method_name)(*args, **kwargs)
        return result

    base_element = the_job_plugin.container.find_by_path(object_path, skip_first=True)
    # Call the method with provided args/kwargs
    result = getattr(base_element, method_name)(*args, **kwargs)
    return result
