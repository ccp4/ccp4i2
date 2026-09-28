import logging

from ....db import models
from ..plugins.get_plugin import get_job_plugin

logger = logging.getLogger(f"ccp4i2:{__name__}")


def get_job_container(the_job: models.Job):
    """
    Return the parameter container for a job, populated from its params file.

    Delegates to :func:`get_job_plugin`, which is the loader the rest of the
    server uses. The obvious-looking alternative --- build a bare
    ``CContainer`` and call ``loadContentsFromXml`` on the task's ``.def.xml``
    --- does not work and does not say so: ``CContainer.loadContentsFromXml``
    routes to ``setEtree(root, ignore_missing=True)``, a ``.def.xml`` root
    (``<ccp4i2><ccp4i2_header/><ccp4i2_body/></ccp4i2>``) is not
    container-shaped, and ``ignore_missing`` swallows the mismatch. The result
    is a container with **zero children** and no error, which is what made the
    i2run-command button render ``<task> --project_name <proj>`` and nothing
    else.

    Args:
        the_job (Job): the job whose parameters to load.

    Returns:
        CContainer: the job's container, or None if the plugin could not load.
    """
    plugin = get_job_plugin(the_job)
    if plugin is None:
        logger.error("No plugin for job %s (task %s)", the_job.id, the_job.task_name)
        return None

    container = plugin.container

    # Keep the plugin alive for as long as the container is. A parent owns its
    # children here: HierarchicalObject.__del__ calls destroy(), which
    # recursively destroys children and deletes them from __dict__. So letting
    # the plugin be collected empties the container the caller is still
    # holding --- `get_job_plugin(job).container` has 0 children where
    # `p = get_job_plugin(job); p.container` has 4.
    container._i2run_owner = plugin

    return container
