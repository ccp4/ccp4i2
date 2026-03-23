import logging
import uuid
from typing import Optional, Any, Union, List

from ccp4i2.core.CCP4Container import CContainer
from ccp4i2.core.base_object.cdata_file import CDataFile
from ccp4i2.core.base_object.fundamental_types import CList

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


def _populate_file_from_context(dobj: CDataFile, context_job_id: str, the_job):
    """Populate a single CDataFile from context job outputs."""
    sub_type = dobj.qualifiers("requiredSubType")
    content_flag = dobj.qualifiers("requiredContentFlag")
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

    logger.info(
        "File search for %s returned %d result(s)",
        dobj._name, len(file_id_list)
    )
    if len(file_id_list) > 0:
        _set_file_from_db(dobj, file_id_list[0], the_job)


def _set_file_from_db(dobj: CDataFile, file_uuid: str, the_job):
    """Set a CDataFile's value from a database File record."""
    the_file = models.File.objects.get(uuid=file_uuid)
    dobj.set(
        {
            "baseName": str(the_file.name),
            "relPath": str(the_file.rel_path),
            "project": str(the_job.project.uuid).replace("-", ""),
            "annotation": str(the_file.annotation),
            "dbFileId": str(the_file.uuid).replace("-", ""),
            "contentFlag": the_file.content,
            "subType": the_file.sub_type,
        }
    )
    dobj.loadFile()
    dobj.setContentFlag()


def _populate_file_lists_from_context(input_data, context_job_id: str, the_job):
    """Find empty CList objects with file subItems and populate from context.

    When a CList like DICT_LIST starts empty, find_all_files() discovers nothing
    inside it. This function finds such lists, queries for matching files from
    the context job, and creates new list items for each match.
    """
    for child in input_data.children():
        if not isinstance(child, CList):
            continue
        # Already populated (items found by find_all_files were handled above)
        if len(child) > 0:
            continue

        sub_item_def = child.get_qualifier('subItem')
        if not sub_item_def or not isinstance(sub_item_def, dict):
            continue

        item_class = sub_item_def.get('class')
        item_qualifiers = sub_item_def.get('qualifiers', {})
        if not item_class or not issubclass(item_class, CDataFile):
            continue

        # Check fromPreviousJob on either the subItem or the CList itself
        has_from_previous = (
            item_qualifiers.get('fromPreviousJob')
            or child.get_qualifier('fromPreviousJob')
        )
        if not has_from_previous:
            continue

        # Get mimeTypeName from subItem qualifiers, or fall back to the
        # class-level QUALIFIERS (e.g. CPdbDataFile has 'chemical/x-pdb')
        mime_type = item_qualifiers.get('mimeTypeName')
        if not mime_type:
            class_quals = getattr(item_class, 'QUALIFIERS', {})
            mime_type = class_quals.get('mimeTypeName')
        if not mime_type:
            continue

        sub_type = _normalize_int_qualifier(
            item_qualifiers.get('requiredSubType')
        )
        content_flag = _normalize_int_qualifier(
            item_qualifiers.get('requiredContentFlag')
        )

        logger.info(
            "Looking for files to populate list %s: mimeType=%s, subType=%s, contentFlag=%s",
            child.objectName(), mime_type, sub_type, content_flag
        )

        file_id_list = get_file_by_job_context(
            contextJobId=context_job_id,
            fileType=mime_type,
            subType=sub_type,
            contentFlag=content_flag,
            projectId=str(the_job.project.uuid),
        )

        logger.info(
            "File search for list %s returned %d result(s)",
            child.objectName(), len(file_id_list)
        )

        for file_uuid in file_id_list:
            new_item = child.makeItem()
            child.append(new_item)
            _set_file_from_db(new_item, file_uuid, the_job)


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
    logger.info(
        "find_all_files() returned %d files: %s",
        len(all_input_files),
        [f"{f._name}({type(f).__name__})" for f in all_input_files]
    )
    dobj_list = [
        f for f in all_input_files
        if f.qualifiers("fromPreviousJob")
    ]

    logger.info(
        "Found %d input files with fromPreviousJob=True: %s",
        len(dobj_list),
        [f._name for f in dobj_list]
    )

    dobj: CDataFile
    for dobj in dobj_list:
        if context_job_id is None:
            dobj.unSet()
            continue

        _populate_file_from_context(dobj, context_job_id, the_job)

    # Handle CList objects whose subItem is a CDataFile with fromPreviousJob=True.
    # When a CList (e.g. DICT_LIST) starts empty, find_all_files() finds nothing
    # inside it because there are no items to traverse. We need to find such lists,
    # query for matching files from the context job, and populate them.
    if context_job_id is not None:
        _populate_file_lists_from_context(input_data, context_job_id, the_job)

    # Only save params if requested (caller may handle saving themselves)
    if save_params:
        save_params_for_job(the_job_plugin, the_job=the_job)
