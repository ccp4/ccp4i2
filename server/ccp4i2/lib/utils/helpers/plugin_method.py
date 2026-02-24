from typing import List

from ccp4i2.db import models
from ..plugins.get_plugin import get_job_plugin


def plugin_method(
    the_job: models.Job,
    method_name: str,
    args: List = None,
    kwargs: dict = None,
):
    """
    Execute a method on a job's CPluginScript instance.

    Args:
        the_job: The Job model instance
        method_name: Name of the method to call on the plugin
        args: Positional arguments to pass to the method
        kwargs: Keyword arguments to pass to the method

    Returns:
        The return value of the called method (must be JSON-serializable)
    """
    if args is None:
        args = []
    if kwargs is None:
        kwargs = {}

    the_job_plugin = get_job_plugin(the_job)
    if the_job_plugin is None:
        raise RuntimeError(
            f"Could not load plugin for job {the_job.id} (task: {the_job.task_name})"
        )

    method = getattr(the_job_plugin, method_name, None)
    if method is None:
        raise AttributeError(
            f"Plugin '{the_job.task_name}' has no method '{method_name}'"
        )

    result = method(*args, **kwargs)
    return result
