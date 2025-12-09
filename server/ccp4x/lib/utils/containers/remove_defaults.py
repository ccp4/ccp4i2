import logging
from typing import Set

# Note that these seem to have to be imported from "core" rather than "ccp4i2.core" for isinstance to work
# MN
from ccp4i2.core import CCP4Container

logger = logging.getLogger(f"ccp4x:{__name__}")


def remove_container_default_values(
    container: CCP4Container.CContainer,
    preserve_containers: Set[str] = None
):
    """
    Remove unset default values from a container hierarchy.

    This function recursively removes objects that have only default values
    (not explicitly set by user) from the container. This is used when creating
    new jobs to clean up unused default parameters before saving.

    IMPORTANT: inputData children are preserved because they need to remain
    available for context-based population (set_input_by_context_job).
    Only outputData and controlParameters have their defaults cleaned up.

    Args:
        container: The CContainer to clean
        preserve_containers: Set of container names to skip (don't remove children from).
                           Defaults to {'inputData'} to preserve input file structure.
    """
    if preserve_containers is None:
        preserve_containers = {'inputData'}

    data_list = container.children()
    for dobj in data_list:
        # dobj = getattr(the_job_plugin.container.outputData, object_name)
        if isinstance(dobj, CCP4Container.CContainer):
            obj_name = dobj.objectName() if hasattr(dobj, "objectName") else None
            # Skip containers in the preserve list
            if obj_name in preserve_containers:
                logger.debug(f"Preserving container {obj_name} (not removing defaults)")
                continue
            if obj_name not in ["temporary"]:
                remove_container_default_values(dobj, preserve_containers)
        else:
            is_set = dobj.isSet(allowUndefined=False, allowDefault=False, allSet=True)
            if (
                not is_set
                and hasattr(dobj, "objectName")
                and dobj.objectName() not in "header"
            ):
                try:
                    container.deleteObject(dobj.objectName())
                except Exception as err:
                    logger.exception(
                        "Error deleting default values %s",
                        container.objectName(),
                        exc_info=err,
                    )
