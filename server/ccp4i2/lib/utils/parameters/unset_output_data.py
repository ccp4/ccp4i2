import logging
from ccp4i2.core import CCP4PluginScript

# Note that these seem to have to be imported from "core" rather than "ccp4i2.core" for isinstance to work
# MN
from ccp4i2.core import CCP4File
from ..containers.find_objects import find_objects

logger = logging.getLogger(f"ccp4i2:{__name__}")


def unset_output_data(the_job_plugin: CCP4PluginScript.CPluginScript):
    """
    Unsets the output data for a given job plugin.
    This function finds all objects of type `CCP4File.CDataFile` within the output data
    of the provided job plugin and calls the `unSet` method on each of them.
    Args:
        the_job_plugin (CCP4PluginScript.CPluginScript): The job plugin containing the output data.
    Returns:
        None
    """

    outputs = find_objects(
        the_job_plugin.container.outputData, lambda a: isinstance(a, CCP4File.CDataFile)
    )
    for output in outputs:
        output.unSet()
