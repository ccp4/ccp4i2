import logging
import uuid
from typing import Optional, Any

from ccp4i2.core.CCP4Container import CContainer
from ccp4i2.core.base_object.cdata_file import CDataFile

from ccp4x.db import models
from ..plugins.get_plugin import get_job_plugin
from .save_params import save_params_for_job
from ..files.get_by_context import get_file_by_job_context

logger = logging.getLogger(f"ccp4x:{__name__}")


def set_input_by_context_job(
    job_id: Optional[str] = None,
    context_job_id: Optional[str] = None,
    plugin: Optional[Any] = None,
    save_params: bool = True,
):
    """
    Populate job inputs from a context job's outputs.

    Finds all input files with `fromPreviousJob=True` qualifier and attempts
    to match them with output files from the context job (or its ancestors).

    Args:
        job_id: UUID string of the job to populate
        context_job_id: UUID string of the context job to draw inputs from
        plugin: Optional plugin instance to use. If not provided, creates a new one.
        save_params: Whether to save params after modification (default True).
                    Set to False if caller will handle saving.
    """
    assert job_id is not None

    logger.info(
        "In set_input_by_context_job for job_id %s context_job %s",
        job_id,
        context_job_id,
    )

    job_uuid = uuid.UUID(job_id)

    the_job = models.Job.objects.get(uuid=job_uuid)

    # Use provided plugin or create a new one
    if plugin is not None:
        the_job_plugin = plugin
    else:
        the_job_plugin = get_job_plugin(the_job)

    the_container: CContainer = the_job_plugin.container
    input_data: CContainer = the_container.inputData

    # Use modern find_all_files() to get all file objects in inputData
    # Then filter for those with fromPreviousJob=True qualifier
    all_input_files = input_data.find_all_files()
    dobj_list = [
        f for f in all_input_files
        if f.qualifiers("fromPreviousJob")
    ]

    logger.debug(
        "Found %d input files with fromPreviousJob=True",
        len(dobj_list)
    )

    dobj: CDataFile
    for dobj in dobj_list:
        if context_job_id is None:
            dobj.unSet()
            continue

        # Parse qualifiers - they may be comma-separated strings representing lists
        sub_type = dobj.qualifiers("requiredSubType")
        content_flag = dobj.qualifiers("requiredContentFlag")

        # Convert comma-separated strings to lists of integers
        if isinstance(sub_type, str) and "," in sub_type:
            sub_type = [int(x.strip()) for x in sub_type.split(",")]
        elif sub_type is not None and sub_type != "":
            try:
                sub_type = int(sub_type)
            except (ValueError, TypeError):
                sub_type = None

        if isinstance(content_flag, str) and "," in content_flag:
            content_flag = [int(x.strip()) for x in content_flag.split(",")]
        elif content_flag is not None and content_flag != "":
            try:
                content_flag = int(content_flag)
            except (ValueError, TypeError):
                content_flag = None

        file_id_list = get_file_by_job_context(
            contextJobId=context_job_id,
            fileType=dobj.qualifiers("mimeTypeName"),
            subType=sub_type,
            contentFlag=content_flag,
            projectId=str(the_job.project.uuid),
        )

        if len(file_id_list) > 0:
            the_file = models.File.objects.get(uuid=file_id_list[0])
            dobj.set(
                {
                    "baseName": str(the_file.name),
                    "relPath": str(the_file.rel_path),
                    "project": str(the_job.project.uuid).replace("-", ""),
                    "annotation": str(the_file.annotation),
                    "dbFileId": str(the_file.uuid).replace("-", ""),
                    "contentFlag": the_file.content,  # DB uses 'content', CDataFile uses 'contentFlag'
                    "subType": the_file.sub_type,
                }
            )
            dobj.loadFile()
            # Auto-detect content flag from file (base CDataFile.setContentFlag
            # handles the case where contentFlag is None or not applicable)
            dobj.setContentFlag()

    # Only save params if requested (caller may handle saving themselves)
    if save_params:
        save_params_for_job(the_job_plugin, the_job=the_job)
