import logging
import pathlib
import getpass
from xml.etree import ElementTree as ET
from ccp4i2.core import CCP4File
from ccp4i2.core import CCP4Utils
from ccp4i2.core import CCP4PluginScript
from ccp4i2.core import CCP4Container
from ccp4i2.db import models

logger = logging.getLogger(f"ccp4i2:{__name__}")


def save_params_for_job(
    the_job_plugin: CCP4PluginScript.CPluginScript,
    the_job: models.Job,
    mode="JOB_INPUT",
    exclude_unset=False,
):
    """
    Save parameters for a given job to an XML file.

    This function generates a file name based on the job plugin and mode,
    relocates the file path to the job's directory, and saves the job's
    parameters in XML format. If a file with the same name already exists,
    it creates a backup before saving the new file.

    Args:
        the_job_plugin (CCP4PluginScript.CPluginScript): The job plugin script instance.
        the_job (models.Job): The job instance containing job details.
        mode (str, optional): The mode for generating the file name. Defaults to "JOB_INPUT".
        exclude_unset (bool, optional): Flag to exclude unset parameters. Defaults to False.

    Returns:
        None
    """
    fileName = the_job_plugin.makeFileName(mode)
    # Rework to the directory of "the_job"
    relocated_file_path = the_job.directory / pathlib.Path(fileName).name

    #print(f"[DEBUG save_params_for_job] mode={mode}, fileName={fileName}, relocated_file_path={relocated_file_path.name}")

    if relocated_file_path.exists():
        CCP4Utils.backupFile(str(relocated_file_path), delete=False)

    # IMPORTANT: Convert pathlib.Path to absolute string path
    # Otherwise CI2XmlDataFile.getFullPath() returns just the filename
    f = CCP4File.CI2XmlDataFile(fullPath=str(relocated_file_path.absolute()))

    # Debug: check what setFullPath did
    pass  # DEBUG: print(f"[DEBUG save_params_for_job] After init: baseName.value={f.baseName.value}, relPath.value={getattr(f.relPath, 'value', 'N/A')}")
    pass  # DEBUG: print(f"[DEBUG save_params_for_job] After init: getFullPath()={f.getFullPath()}")

    # CI2XmlDataFile already has a header instance created during __init__
    # Just populate it with job information
    f.header.function.set("PARAMS")
    f.header.projectName.set(the_job.project.name)
    f.header.projectId.set(str(the_job.project.uuid))
    f.header.jobNumber.set(the_job.number)
    f.header.jobId.set(str(the_job.uuid))
    f.header.setCurrent()
    f.header.pluginName.set(the_job.task_name)
    f.header.userId.set(getpass.getuser())

    # Use the exclude_unset parameter passed to this function
    body_etree = the_job_plugin.container.getEtree(excludeUnset=exclude_unset)
    ET.indent(body_etree, space="  ")
    #print(f"[DEBUG save_params_for_job] Generated body_etree with excludeUnset={exclude_unset} {ET.tostring(body_etree, encoding='unicode')}")
    f.saveFile(bodyEtree=body_etree)
