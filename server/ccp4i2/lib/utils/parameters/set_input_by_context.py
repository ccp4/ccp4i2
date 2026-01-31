import logging
import uuid
from typing import Optional, Any, Union, List

from ccp4i2.core.CCP4Container import CContainer
from ccp4i2.core.base_object.cdata_file import CDataFile

from ccp4i2.db import models
from ..plugins.get_plugin import get_job_plugin
from .save_params import save_params_for_job
from ..files.get_by_context import get_file_by_job_context

logger = logging.getLogger(f"ccp4i2:{__name__}")


def _normalize_int_qualifier(value) -> Union[int, List[int], None]:
    """
    Normalize a qualifier value that should be an int or list of ints.

    Handles values that may come as:
    - list of ints (already parsed by def_xml_handler)
    - comma-separated string (legacy)
    - single int or string representation of int
    - None or empty string

    Returns:
        int, list of ints, or None
    """
    if value is None or value == "":
        return None

    # Already a list - return as-is (from def_xml_handler parsing)
    if isinstance(value, list):
        return value

    # Comma-separated string
    if isinstance(value, str) and "," in value:
        return [int(x.strip()) for x in value.split(",")]

    # Single value - convert to int
    try:
        return int(value)
    except (ValueError, TypeError):
        return None


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

        # Normalize qualifiers to int or list of ints
        # These may come as: list (from def_xml_handler), comma-separated string, single int/string, or None
        sub_type = _normalize_int_qualifier(sub_type)
        content_flag = _normalize_int_qualifier(content_flag)

        logger.info(
            "Looking for file to populate %s: mimeType=%s, subType=%s, contentFlag=%s",
            dobj.objectPath(), dobj.qualifiers("mimeTypeName"), sub_type, content_flag
        )

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
