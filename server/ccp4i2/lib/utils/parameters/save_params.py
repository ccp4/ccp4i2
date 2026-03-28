# Copyright (C) 2025-2026 University of York
# Copyright (C) 2025-2026 Newcastle University
#
# This file is part of CCP4i2.
#
# CCP4i2 is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License version 3,
# modified in accordance with the provisions of the license to address
# the requirements of UK law.
#
# See https://www.ccp4.ac.uk/ccp4license.php for details.
import logging
import pathlib

from ccp4i2.core import CCP4Utils
from ccp4i2.core import CCP4PluginScript
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

    Ensures the plugin has DB context so the header is fully populated,
    then delegates to saveDataToXml() — the single serialization path.

    Args:
        the_job_plugin: The job plugin script instance.
        the_job: The job instance containing job details.
        mode: The mode for generating the file name. Defaults to "JOB_INPUT".
        exclude_unset: Flag to exclude unset parameters. Defaults to False.
    """
    fileName = the_job_plugin.makeFileName(mode)
    # Rework to the directory of "the_job"
    relocated_file_path = the_job.directory / pathlib.Path(fileName).name

    if relocated_file_path.exists():
        CCP4Utils.backupFile(str(relocated_file_path), delete=False)

    # Ensure the plugin has DB context for a complete header
    if getattr(the_job_plugin, '_dbJobId', None) is None:
        the_job_plugin._dbJobId = str(the_job.uuid)
    if getattr(the_job_plugin, '_dbJobNumber', None) is None:
        the_job_plugin._dbJobNumber = the_job.number
    if getattr(the_job_plugin, '_dbProjectId', None) is None:
        the_job_plugin._dbProjectId = str(the_job.project.uuid)
    if getattr(the_job_plugin, '_dbProjectName', None) is None:
        the_job_plugin._dbProjectName = the_job.project.name

    the_job_plugin.saveDataToXml(
        str(relocated_file_path.absolute()), exclude_unset=exclude_unset
    )
