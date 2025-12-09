import logging
from pathlib import Path
from django.utils.text import slugify
from ccp4i2.core import CCP4ErrorHandling
from ccp4i2.core import CCP4Container
from ccp4i2.core import CCP4File
from ccp4i2.core import CCP4ModelData
from ccp4x.db import models

logger = logging.getLogger(f"ccp4x:{__name__}")


def set_output_file_names(
    container: CCP4Container.CContainer = None,
    projectId: str = None,
    jobNumber: str = None,
    force: bool = True,
    fancy: bool = False,
):
    """
    Sets the output file names for a given job in a CCP4 project.
    Args:
        container (CCP4Container.CContainer, optional): The container holding the output data. Defaults to None.
        projectId (str, optional): The UUID of the project. Defaults to None.
        jobNumber (str, optional): The job number in the project. Defaults to None.
        force (bool, optional): If True, forces the setting of output paths even if they are already set. Defaults to True.
    Returns:
        CCP4ErrorHandling.CErrorReport: An error report object containing any errors encountered during the process.
    """
    myErrorReport = CCP4ErrorHandling.CErrorReport()
    relPath = Path("CCP4_JOBS").joinpath(
        *[f"job_{numberElement}" for numberElement in jobNumber.split(".")]
    )
    the_job = models.Job.objects.get(project__uuid=projectId, number=jobNumber)
    jobName = (
        f"{jobNumber}_{slugify(the_job.project.name)}_{the_job.task_name}_"
        if fancy
        else ""
    )
    dataList = container.outputData.dataOrder()
    logger.info(f"set_output_file_names: Found {len(dataList)} output objects: {dataList}")
    for objectName in dataList:
        try:
            dobj = container.outputData.find(objectName)
            logger.info(f"  Processing {objectName}: type={type(dobj).__name__}, isinstance(CDataFile)={isinstance(dobj, CCP4File.CDataFile)}")
            if isinstance(dobj, CCP4File.CDataFile) and (force or not dobj.isSet()):
                logger.info(f"    Calling setOutputPath for {objectName}")
                dobj.setOutputPath(
                    jobName=jobName, projectId=projectId, relPath=str(relPath)
                )
                logger.info(f"    setOutputPath completed for {objectName}")
            if isinstance(dobj, CCP4ModelData.CPdbDataFile):
                logger.debug(f"    [DEBUG] CPdbDataFile detected: {objectName}")
                logger.debug(f"    [DEBUG] Current baseName: {dobj.baseName}")
                oldBaseName = str(Path(str(dobj.baseName)).stem)
                logger.debug(f"    [DEBUG] oldBaseName (stem): {oldBaseName}")
                # Use fileExtensions() method which returns extensions based on contentFlag
                # This indirects through the class's extension array
                if hasattr(dobj, 'fileExtensions'):
                    extensions = dobj.fileExtensions()
                    logger.debug(f"    [DEBUG] fileExtensions() returned: {extensions}")
                    if extensions and len(extensions) > 0:
                        primary_ext = extensions[0]
                        logger.debug(f"    [DEBUG] primary_ext: {primary_ext}")
                        dobj.baseName.set(f"{oldBaseName}.{primary_ext}")
                        logger.info(f"    Set baseName to {oldBaseName}.{primary_ext} based on fileExtensions()")
                else:
                    # Fallback: default to .pdb
                    logger.debug(f"    [DEBUG] No fileExtensions method, using fallback .pdb")
                    dobj.baseName.set(f"{oldBaseName}.pdb")

        except Exception as err:
            logger.exception(
                "Exception in setOutputFileNames for %s",
                dobj.objectPath(),
                exc_info=err,
            )
    return myErrorReport
