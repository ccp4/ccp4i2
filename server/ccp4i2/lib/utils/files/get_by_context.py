import logging
import uuid
from ccp4i2.db import models

logger = logging.getLogger(f"ccp4i2:{__name__}")


def get_file_by_job_context(
    contextJobId: str = None,
    fileType: str = None,
    subType: int = None,
    contentFlag: int = None,
    projectId: str = None,
) -> list:

    assert contextJobId is not None
    assert fileType is not None

    context_job_uuid = uuid.UUID(contextJobId)
    project_uuid = None
    if projectId is not None and len(projectId.strip()) > 0:
        project_uuid = uuid.UUID(projectId)

    context_job_in = models.Job.objects.get(uuid=context_job_uuid)
    project_jobs = models.Job.objects.filter(
        parent__isnull=True, project__uuid=context_job_in.project.uuid
    )
    search_jobs_list = [
        project_job
        for project_job in project_jobs
        if int(project_job.number) <= int(context_job_in.number)
    ]
    search_jobs_list = sorted(search_jobs_list, key=lambda a: int(a.number))

    while len(search_jobs_list) > 0:
        context_job = search_jobs_list.pop()
        logger.info("Looking for context in %s", context_job.number)

        output_file_ids, import_file_ids = _get_job_files(
            context_job, fileType, subType, contentFlag, project_uuid
        )

        if len(output_file_ids) > 0:
            if context_job.task_name == "coot_rebuild":
                return output_file_ids[::-1]
            else:
                return output_file_ids
        elif len(import_file_ids) > 0:
            return import_file_ids

        input_file_ids = _get_file_uses(
            context_job, fileType, subType, contentFlag, project_uuid
        )

        if len(input_file_ids) > 0:
            if context_job.task_name == "coot_rebuild":
                return input_file_ids[::-1]
            else:
                return input_file_ids

    return []


def _get_job_files(
    context_job: models.Job,
    fileType: str,
    subType: int,
    contentFlag: int,
    project_uuid: uuid.UUID,
):
    filter_dict = {"type__name": fileType}
    if isinstance(subType, list) and 0 not in subType and subType is not None:
        filter_dict["sub_type__in"] = subType
    elif not isinstance(subType, list) and subType != 0 and subType is not None:
        filter_dict["sub_type"] = subType
    if contentFlag is not None and contentFlag is not NotImplemented:
        if not isinstance(contentFlag, list):
            filter_dict["content"] = contentFlag
        else:
            filter_dict["content__in"] = contentFlag
    if project_uuid is not None:
        filter_dict["job__project__uuid"] = project_uuid
    logger.debug("Using filter_dict %s", filter_dict)
    file_qs = models.File.objects.filter(**filter_dict)

    job_file_qs = file_qs.filter(job=context_job)

    # Filter out files that don't exist on disk
    # This prevents autopopulation from setting references to deleted/missing files
    existing_job_files = []
    for jobFile in job_file_qs:
        if jobFile.path and jobFile.path.exists():
            existing_job_files.append(jobFile)
        else:
            logger.warning(
                "Skipping file %s (%s) from job %s - file does not exist on disk at %s",
                jobFile.name, jobFile.uuid, context_job.number, jobFile.path
            )

    job_file_id_list = [str(jobFile.uuid) for jobFile in existing_job_files]

    file_imports = models.FileImport.objects.filter(file__uuid__in=job_file_id_list)
    import_file_ids = [str(importFile.file.uuid) for importFile in file_imports]
    output_file_ids = [
        str(jobFileId)
        for jobFileId in job_file_id_list
        if jobFileId not in import_file_ids
    ]

    return output_file_ids, import_file_ids


def _get_file_uses(
    context_job: models.Job,
    fileType: str,
    subType: int,
    contentFlag: int,
    project_uuid: uuid.UUID,
):
    filter_dict = {"file__type__name": fileType}
    if isinstance(subType, list) and 0 not in subType and subType is not None:
        filter_dict["file__sub_type__in"] = subType
    elif not isinstance(subType, list) and subType != 0 and subType is not None:
        filter_dict["file__sub_type"] = subType
    if contentFlag is not None and contentFlag is not NotImplemented:
        if not isinstance(contentFlag, list):
            filter_dict["file__content"] = contentFlag
        else:
            filter_dict["file__content__in"] = contentFlag
    if project_uuid is not None:
        filter_dict["job__project__uuid"] = project_uuid

    fileuse_qs = models.FileUse.objects.filter(**filter_dict)
    jobfileuse_qs = fileuse_qs.filter(job=context_job)

    # Helper to check if file exists on disk
    def file_exists(file_use):
        file_obj = file_use.file
        if file_obj.path and file_obj.path.exists():
            return True
        logger.warning(
            "Skipping file use %s (%s) from job %s - file does not exist on disk at %s",
            file_obj.name, file_obj.uuid, context_job.number, file_obj.path
        )
        return False

    output_file_qs = jobfileuse_qs.filter(role=0)
    output_id_list = [
        str(outputFileUse.file.uuid)
        for outputFileUse in output_file_qs
        if file_exists(outputFileUse)
    ]

    inputfile_qs = jobfileuse_qs.filter(role=1)
    input_id_list = [
        str(inputFile.file.uuid)
        for inputFile in inputfile_qs
        if file_exists(inputFile)
    ]

    return output_id_list + input_id_list
