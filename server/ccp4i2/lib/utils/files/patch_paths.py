import os
import logging
from django.utils.text import slugify
from ccp4i2.core.base_object.cdata_file import CDataFile
from ccp4i2.core.CCP4ModelData import CPdbDataFile
from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core import CCP4File
from ccp4i2.core import CCP4Data
from ccp4i2.core import CCP4ModelData
from ccp4i2.db import models

logger = logging.getLogger(f"ccp4i2:{__name__}")


def patch_output_file_paths(
    job_plugin: CPluginScript, job: models.Job, instance_string: str = ""
):
    container = job_plugin.container
    dataList = container.outputData.dataOrder()
    for objectName in dataList:
        dobj = container.outputData.find(objectName)
        if isinstance(dobj, (CDataFile, CCP4File.CDataFile)):
            handle_cdatafile(job, dobj)
        elif isinstance(dobj, CCP4Data.COutputFileList):
            # Empty the current output file list array
            while len(dobj) > 0:
                dobj.remove(dobj[-1])
            for iItem in range(dobj.qualifiers()["listMaxLength"]):
                new_instance = dobj.makeItem()
                dobj.append(new_instance)
                handle_cdatafile(
                    job,
                    new_instance,
                    instance_string="_" + str(iItem),
                )


def handle_cdatafile(job: models.Job, dobj, instance_string=""):
    partPath = os.path.join(
        f"{job.number}_{job.project.name}_{dobj.objectName()}_{job.task_name}"
    )
    partPath = str(slugify(str(partPath)))

    # Get file extensions list - check if method exists
    if hasattr(dobj, 'fileExtensions') and callable(getattr(dobj, 'fileExtensions')):
        extensions = dobj.fileExtensions()
    else:
        # Fallback: try to get from qualifiers if method doesn't exist
        extensions = dobj.get_qualifier('fileExtensions', [])

    # For PDB files with contentFlag=2 (mmCIF format), use second extension if available
    if (
        isinstance(dobj, (CPdbDataFile, CCP4ModelData.CPdbDataFile))
        and hasattr(dobj, 'contentFlag')
        and dobj.contentFlag.isSet()
        and int(dobj.contentFlag) == 2
        and len(extensions) > 1
    ):
        extension = extensions[1]
    else:
        # Use first extension (or empty string if no extensions defined)
        extension = extensions[0] if extensions else ""

    fullPath = os.path.join(
        str(job.directory),
        partPath + instance_string + "." + extension,
    )
    # MN 2020-07-09
    # Set the full path for the output file
    # NOTE: setFullPath() signature changed - checkDb parameter removed
    # The method now handles path setting directly without needing checkDb flag
    dobj.setFullPath(fullPath)
