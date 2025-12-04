import logging

# Note that these seem to have to be imported from "core" rather than "ccp4i2.core" for isinstance to work
# MN
from core import CCP4Container

logger = logging.getLogger(f"ccp4x:{__name__}")


def remove_container_default_values(container: CCP4Container.CContainer):
    # Unset the output file data, which should be recalculated for the new plugin, I guess
    data_list = container.children()
    for dobj in data_list:
        # dobj = getattr(the_job_plugin.container.outputData, object_name)
        if isinstance(dobj, CCP4Container.CContainer):
            if hasattr(dobj, "objectName") and dobj.objectName() not in [
                "temporary",
            ]:
                remove_container_default_values(dobj)
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
