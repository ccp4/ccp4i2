import datetime
import logging
import os
import pathlib
import shutil
import uuid
import re

from django.utils import timezone

from ccp4i2.core import CCP4Container
from ccp4i2.core.base_object.cdata_file import CDataFile
from ccp4i2.core import CCP4File
from ccp4i2.core import CCP4PluginScript
from ccp4i2.core.base_object.fundamental_types import CList
from ccp4i2.core.CCP4File import CI2XmlDataFile as CI2XmlDataFile

from ccp4x.db import models
from ..parameters.save_params import save_params_for_job
from ..containers.find_objects import find_objects


logger = logging.getLogger(f"ccp4x:{__name__}")


def extract_from_first_bracketed(path: str) -> str:
    parts = path.split(".")
    for i, part in enumerate(parts):
        if re.search(r"\[.*\]", part):
            return ".".join(parts[i:])
    return parts[-1]  # fallback to the last part if no bracketed part found


def _process_input(
    theJob: models.Job,
    plugin: CCP4PluginScript.CPluginScript,
    input: CCP4File.CDataFile,
):
    theFile = None
    if input.dbFileId is not None and len(str(input.dbFileId)) != 0:
        theFile = models.File.objects.get(uuid=uuid.UUID(str(input.dbFileId)))
    else:
        if input.baseName is not None and len(str(input.baseName).strip()) != 0:
            sourceFilePath = pathlib.Path(str(input.relPath)) / str(input.baseName)
            if not sourceFilePath.exists():
                sourceFilePath = (
                    pathlib.Path(theJob.project.directory)
                    / str(input.relPath)
                    / str(input.baseName)
                )
            # Load file
            input.loadFile()
            input.setContentFlag()
            destFilePath = (
                pathlib.Path(theJob.project.directory)
                / "CCP4_IMPORTED_FILES"
                / sourceFilePath.name
            )
            while destFilePath.exists():
                fileRoot, fileExt = os.path.splitext(destFilePath.name)
                destFilePath = destFilePath.parent / "{}_1{}".format(fileRoot, fileExt)
            logger.debug("src %s, UniqueDestFilePath %s", sourceFilePath, destFilePath)
            shutil.copyfile(sourceFilePath, destFilePath)
            # Now have to change the plugin to reflect the new location

            try:

                file_mime_type = input.qualifiers("mimeTypeName")
                if len(file_mime_type) == 0:
                    logger.error(
                        "Class %s Does not have an associated mimeTypeName....ASK FOR DEVELOPER FIX",
                        str(input.objectName()),
                    )
                    if isinstance(input, CCP4File.CI2XmlDataFile) or isinstance(
                        input, CI2XmlDataFile
                    ):
                        file_mime_type = "application/xml"
                file_type_object = models.FileType.objects.get(name=file_mime_type)

                # print('What I know about import is', str(valueDictForObject(input)))
                if (
                    not hasattr(input, "annotation")
                    or input.annotation is None
                    or len(str(input.annotation).strip()) == 0
                ):
                    annotation = f"Imported from file {sourceFilePath.name}"
                else:
                    annotation = str(input.annotation)
                job_param_name = extract_from_first_bracketed(input.objectPath())
                createDict = {
                    "name": str(destFilePath.name),
                    "annotation": annotation,
                    "type": file_type_object,
                    "job": theJob,
                    "job_param_name": job_param_name,
                    "directory": 2,
                }
                # print(createDict)
                theFile = models.File(**createDict)
                theFile.save()

                input.dbFileId.set(theFile.uuid)
                input.project.set(str(theJob.project.uuid))
                input.relPath.set("CCP4_IMPORTED_FILES")
                input.baseName.set(destFilePath.name)
                input.loadFile()
                input.setContentFlag(reset=True)

                createDict = {
                    "file": theFile,
                    "name": str(sourceFilePath),
                    "time": timezone.now(),
                    "last_modified": timezone.now(),
                    "checksum": input.checksum(),
                }
                # print(createDict)
                newImportfile = models.FileImport(**createDict)
                newImportfile.save()
                # for key in createDict:
                #    setattr(newImportfile, key, createDict[key])
                #    newImportfile.save()
            except ValueError as err:
                theFile = None
                logger.error(
                    f"Encountered issue - {err} importing file {input.baseName}"
                )

    if theFile is not None:
        theRole = models.FileUse.Role.IN
        job_param_name = extract_from_first_bracketed(input.objectPath())
        createDict = {
            "file": theFile,
            "job": theJob,
            "role": theRole,
            "job_param_name": job_param_name,
        }
        already_exists = models.FileUse.objects.filter(**createDict).exists()
        if not already_exists:
            fileUse = models.FileUse(**createDict)
            fileUse.save()
        else:
            logger.warning(
                "FileUse already exists for file %s job %s role %s job_param_name %s",
                theFile,
                theJob,
                theRole,
                job_param_name,
            )


def import_files(theJob, plugin):
    inputs = find_objects(
        plugin.container.inputData,
        lambda a: isinstance(
            a,
            (
                CCP4File.CDataFile,
                CDataFile,
            ),
        ),
        True,
    )
    logger.debug("In import_files %s", len(inputs))
    for the_input in inputs:
        _process_input(theJob, plugin, the_input)

    # removeDefaults(plugin.container)
    save_params_for_job(plugin, theJob, mode="JOB_INPUT")

    return plugin
