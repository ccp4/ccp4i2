import logging
import uuid
from typing import List
from xml.etree import ElementTree as ET

from core import CCP4Container
from core.CCP4Container import CContainer
from core import CCP4File
from core import CCP4Data
from core.base_object.fundamental_types import CList
from core.base_object.cdata_file import CDataFile

from ccp4x.db import models
from ..plugins.get_plugin import get_job_plugin
from .save_params import save_params_for_job
from ..files.get_by_context import get_file_by_job_context
from ..containers.find_objects import find_objects
from typing import Optional

logger = logging.getLogger(f"ccp4x:{__name__}")


def set_input_by_context_job(
    job_id: Optional[str] = None, context_job_id: Optional[str] = None
):

    assert job_id is not None

    logger.info(
        "In set_input_by_context_job for job_id %s context_job %s",
        job_id,
        context_job_id,
    )

    job_uuid = uuid.UUID(job_id)

    the_job = models.Job.objects.get(uuid=job_uuid)
    the_job_plugin = get_job_plugin(the_job)
    the_container: CContainer = the_job_plugin.container
    input_data: CContainer = the_container.inputData

    def cdata_func(item):
        if isinstance(item, CCP4File.CDataFile):
            return item.qualifiers("fromPreviousJob")
        return False

    dobj_list = find_objects(
        input_data,
        cdata_func,
        multiple=True,
    )

    # Now find input data that is a list of files, add an item of such lists to the list of dObjs for which we need to set the input
    def clist_func(item):
        if isinstance(item, CCP4Data.CList):
            return item.qualifiers("fromPreviousJob")
        return False

    list_list = find_objects(
        input_data,
        clist_func,
        multiple=True,
    )
    a_list: CList
    for a_list in list_list:
        try:
            a_list.unSet()
            # We may need to make a new item of the list, as the list is empty
            if len(a_list) == 0:
                a_list_item = a_list.makeItem()
            else:
                a_list_item = a_list[0]
            # We need to check if the item is a file, as it may be a list of other things
            if isinstance(a_list_item, CCP4File.CDataFile):
                a_list.append(a_list_item)
                dobj_list.append(a_list_item)
        except Exception as err:
            logger.exception(
                "Exception in set_input_by_context_job for %s",
                a_list.objectPath(),
                exc_info=err,
            )

    dobj: CDataFile
    for dobj in dobj_list:
        if context_job_id is None:
            dobj.unSet()
            continue
        file_id_list = get_file_by_job_context(
            contextJobId=context_job_id,
            fileType=dobj.qualifiers("mimeTypeName"),
            subType=dobj.qualifiers("requiredSubType"),
            contentFlag=dobj.qualifiers("requiredContentFlag"),
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
                    "content": the_file.content,
                    "subType": the_file.sub_type,
                }
            )
            dobj.loadFile()
            dobj.setContentFlag(reset=True)

    save_params_for_job(the_job_plugin, the_job=the_job)
