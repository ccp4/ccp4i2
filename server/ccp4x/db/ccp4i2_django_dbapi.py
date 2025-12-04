from datetime import datetime
import logging
import pathlib
import sys
import traceback
import uuid

from ..lib.utils.files.get_by_context import get_file_by_job_context
from . import models
from .ccp4i2_static_data import FILETYPELIST
from ..lib.utils.jobs.directory import job_directory

logger = logging.getLogger(f"ccp4x:{__name__}")


def get_file_type_id_from_mimetype(mimetype):
    """Convert mimetype back to numeric file type ID using FILETYPELIST."""
    # Handle case where mimetype is a FileType model object
    if hasattr(mimetype, "name"):
        mimetype_str = mimetype.name
    else:
        mimetype_str = str(mimetype)

    for file_type_id, file_mimetype, file_type_name in FILETYPELIST:
        if file_mimetype == mimetype_str:
            return file_type_id
    # Return a default if not found
    return 0


def get_file_type_name_from_mimetype(mimetype):
    """Convert mimetype to readable name using FILETYPELIST."""
    # Handle case where mimetype is a FileType model object
    if hasattr(mimetype, "name"):
        mimetype_str = mimetype.name
    else:
        mimetype_str = str(mimetype)

    for file_type_id, file_mimetype, file_type_name in FILETYPELIST:
        if file_mimetype == mimetype_str:
            return file_type_name
    # Return empty string if not found
    return ""


project_field_old_to_new = {
    "followfromjobid": "follow_from_job",
    "i1projectdirectory": "i1_project_directory",
    "i1projectname": "i1_project_name",
    "lastaccess": "last_access",
    "lastjobnumber": "last_job_number",
    # "parentprojectid": "parentprojectid",
    "projectcreated": "creation_time",
    "projectdirectory": "directory",
    "projectid": "uuid",
    "projectname": "name",
    # "userid": "userid",
}
project_field_new_to_old = {
    item[1]: item[0] for item in project_field_old_to_new.items()
}
project_field_new_to_old["follow_from_job_id"] = "followfromjobid"

job_field_old_to_new = {
    "jobid": "uuid",
    "jobnumber": "number",
    "parentjobid": "parent__uuid",
    "projectid": "project__uuid",
    "taskname": "task_name",
    "title": "title",
    "status": "status",
    "comments": "comments",
    "creationtime": "creation_time",
    "finishtime": "finish_time",
    "projectname": "project__name",
    "projectdirectory": "project__directory",
    "jobtitle": "title",  # Alias for title
}
job_field_new_to_old = {item[1]: item[0] for item in job_field_old_to_new.items()}

# Fields that need special handling (not directly in the Django model)
# These are computed/derived fields that getJobInfo handles specially
JOB_SPECIAL_FIELDS = {
    "runtime",        # Compute from creation_time and finish_time
    "descendentjobs", # Compute from child jobs
    "taskversion",    # Not stored, return placeholder
}

file_field_old_to_new = {
    "filecontent": "content",
    "subtype": "sub_type",
    "fileid": "uuid",
    "jobid": "job__uuid",
    "jobparamname": "job_param_name",
    "pathflag": "directory",
    "filename": "name",
    "annotation": "annotation",
    "filetypeid": "type",
    "filetype": "type__name",  # Get the type name via FK relation
    "jobnumber": "job__number",  # Job number via FK
    "projectname": "job__project__name",  # Project name via FK chain
    "projectid": "job__project__uuid",  # Project UUID via FK chain
}
file_field_new_to_old = {item[1]: item[0] for item in file_field_old_to_new.items()}

# File fields that need special handling (computed from other fields)
FILE_SPECIAL_FIELDS = {
    "relpath",  # Computed from directory flag and job number
}

# Directory flag constants (matching legacy CCP4DbApi)
PATH_FLAG_JOB_DIR = 0
PATH_FLAG_IMPORT_DIR = 1


class CCP4i2DjangoDbApi(object):

    class FakeSignal:
        def emit(self, *arg, **kwarg):
            logger.info("CCP4i2DjangoDbApi been asked to emit %s, %s", arg, kwarg)

    def __init__(self):
        self.projectReset = CCP4i2DjangoDbApi.FakeSignal()
        super().__init__()

    def __getattribute__(self, __name):
        logger.debug("CCP4i2DjangoDbApi being interrogated for %s", __name)
        return super().__getattribute__(__name)

    def getFileByJobContext(
        self,
        contextJobId: str = None,
        fileType: str = None,
        subType: int = None,
        contentFlag: int = None,
        projectId: str = None,
    ) -> list:
        assert contextJobId is not None
        assert fileType is not None
        return get_file_by_job_context(
            contextJobId, fileType, subType, contentFlag, projectId
        )

    def getTaskNameLookup(self, projectId=None, jobId=None, extras=False):
        # Fixme....this should produce a lookup of subtasks sfor use in CCP4i2 purgeJob
        try:
            logger.warning(
                "In unimplemented routine getTaskNameLookup %s, %s, %s",
                projectId,
                jobId,
                extras,
            )
        except Exception as err:
            logger.exception(
                "Err in unimplemented routine getTaskNameLookup", exc_info=err
            )
        return {}

    def getProjectInfo(
        self, projectId=None, projectName=None, mode="all", checkPermission=True
    ):
        """
        Retrieve project information based on project ID or project name.

        Args:
            projectId (str, optional): The unique identifier of the project. Defaults to None.
            projectName (str, optional): The name of the project. Defaults to None.
            mode (str, optional): The mode of information retrieval. Defaults to "all".
            checkPermission (bool, optional): Flag to check permissions. Defaults to True.

        Returns:
            dict or any: The project information. If only one field is requested, returns the value of that field.
                         Returns None if an error occurs.

        Raises:
            Exception: Logs any exceptions that occur during the retrieval process.
        """
        if projectId is not None and "-" not in projectId:
            projectId = uuid.UUID(projectId)
        try:
            the_qs = self._get_project_queryset(projectId, projectName)
            arg = self._get_mode_arguments(mode)
            unpatched_values = the_qs.values(*arg)
            values = self._get_values_from_queryset(
                unpatched_values, project_field_new_to_old
            )
            result = list(values)[0]
            if len(arg) == 1:
                return result[project_field_new_to_old[arg[0]]]
            return result
        except Exception as err:
            logger.exception("Err in getProjectInfo", exc_info=err)
        return None

    def _get_project_queryset(self, projectId, projectName):
        if projectId is None:
            return models.Project.objects.filter(name=projectName)
        else:
            return models.Project.objects.filter(uuid=projectId)

    def _get_mode_arguments(self, mode):
        if isinstance(mode, list):
            return [item.lower() for item in mode]
        elif mode.lower() == "all":
            return [key for key in project_field_new_to_old.keys()]
        else:
            return [project_field_old_to_new[mode]]

    def _get_values_from_queryset(self, unpatched_values, substitution_dict):
        values = []
        for unPatchedValue in unpatched_values:
            value = {}
            for key in unPatchedValue:
                if key in substitution_dict:
                    value[substitution_dict[key]] = self._to_simple_types(
                        unPatchedValue[key]
                    )
            values.append(value)
        return values

    def _to_simple_types(self, value):
        if isinstance(value, uuid.UUID):
            return str(value)
        elif isinstance(value, datetime):
            return value.timestamp()
        else:
            return value

    def deleteFilesOnJobNumberAndParamName(self, projectId=None, jobNumberParamList=[]):
        try:
            for jobNumberParam in jobNumberParamList:
                # print('Requested to delete', jobNumberParam)
                try:
                    file = models.File.objects.get(
                        job__project__uuid=projectId,
                        job__number=jobNumberParam[0],
                        job_param_name=jobNumberParam[1],
                    )
                    file.delete()
                except models.File.DoesNotExist:
                    logger.warning(
                        "Err in deleteFilesOnJobNumberAndParamName %s %s %s",
                        projectId,
                        jobNumberParam[0],
                        jobNumberParam[1],
                    )
        except Exception as err:
            logger.exception(
                "Err in deleteFilesOnJobNumberAndParamName",
                exc_info=err,
                stack_info=True,
            )
        return None

    def getFileInfo(self, fileId=None, mode="all", returnType=None):
        assert fileId is not None
        the_file_qs = models.File.objects.filter(uuid=uuid.UUID(fileId))
        the_file = list(the_file_qs)[0]

        if isinstance(mode, list):
            requested_fields = [item.lower() for item in mode]
        elif mode.lower() == "all":
            requested_fields = [key for key in file_field_old_to_new.keys()]
        else:
            requested_fields = [mode.lower()]

        # Separate regular fields from special computed fields
        special_fields_requested = [f for f in requested_fields if f in FILE_SPECIAL_FIELDS]
        regular_fields = [f for f in requested_fields if f not in FILE_SPECIAL_FIELDS]

        # Will need corrected for some cases
        replacements = file_field_old_to_new

        def patch(label):
            return replacements.get(label, label)

        arg = list(map(patch, regular_fields))

        if arg:
            unpatched_values = the_file_qs.values(*arg)
            listOfDicts = []
            for unPatchedValue in unpatched_values:
                # outer loop over jobs matching jobId
                value = {}
                for key in unPatchedValue:
                    # inner loop over parameters
                    if key.endswith("_id"):
                        value[key[:-3]] = unPatchedValue[key]
                    else:
                        value[key] = unPatchedValue[key]
                listOfDicts.append(value)
            result = listOfDicts[0]
        else:
            result = {}

        # Add computed/special fields
        for field in special_fields_requested:
            if field == "relpath":
                # Compute relative path from directory flag and job number
                path_flag = the_file.directory
                if path_flag == PATH_FLAG_JOB_DIR:
                    # Build path like CCP4_JOBS/job_1/job_2 from job number "1.2"
                    job_number = the_file.job.number
                    job_path_parts = [f"job_{n}" for n in job_number.split(".")]
                    result["relpath"] = str(pathlib.Path("CCP4_JOBS") / pathlib.Path(*job_path_parts))
                elif path_flag == PATH_FLAG_IMPORT_DIR:
                    result["relpath"] = "CCP4_IMPORTED_FILES"
                else:
                    result["relpath"] = ""

        if len(requested_fields) == 1 and not special_fields_requested and returnType != dict:
            return result[arg[0]]
        elif returnType == list:
            return [item[1] for item in result.items()]
        return result

    def getJobInfo(
        self, jobId=None, mode="all", projectName=None, jobNumber=None, returnType=None
    ):
        try:
            logger.info(f"jobId is {jobId}")
            if jobId is None:
                the_job_qs = models.Job.objects.filter(
                    project__name=projectName, number=jobNumber
                )
            else:
                if "-" not in str(jobId):
                    jobId = uuid.UUID(jobId)
                the_job_qs = models.Job.objects.filter(uuid=jobId)
            assert (
                len(list(the_job_qs)) == 1
            ), f"Expected 1 job, got {len(list(the_job_qs))} for jobId {jobId}"

            the_job = list(the_job_qs)[0]

            if isinstance(mode, list):
                requested_fields = [item.lower() for item in mode]
            elif mode.lower() == "all":
                requested_fields = [item.lower() for item in job_field_old_to_new.keys()]
            else:
                requested_fields = [mode.lower()]

            # Separate regular fields from special computed fields
            special_fields_requested = [f for f in requested_fields if f in JOB_SPECIAL_FIELDS]
            regular_fields = [f for f in requested_fields if f not in JOB_SPECIAL_FIELDS]

            def patch(label):
                return job_field_old_to_new.get(label, label)

            arg = list(map(patch, regular_fields))

            # Query the database for regular fields
            if arg:
                unpatched_values = the_job_qs.values(*arg)
                values = self._get_values_from_queryset(
                    unpatched_values, job_field_new_to_old
                )
                result = list(values)[0]
            else:
                result = {}

            # Handle single field request
            if len(requested_fields) == 1 and not special_fields_requested:
                return result[job_field_new_to_old.get(arg[0], requested_fields[0])]

            # Add computed/special fields
            for field in special_fields_requested:
                if field == "runtime":
                    # Compute runtime from creation_time and finish_time
                    if the_job.finish_time and the_job.creation_time:
                        result["runtime"] = (the_job.finish_time - the_job.creation_time).total_seconds()
                    else:
                        result["runtime"] = None
                elif field == "descendentjobs":
                    # Get descendant jobs as list of (job_id, [child_ids]) tuples
                    result["descendentjobs"] = self._get_descendent_jobs(the_job)
                elif field == "taskversion":
                    # Task version not stored in Django model, return placeholder
                    result["taskversion"] = "1.0"

            result["fileroot"] = str(the_job.directory)

            jobFiles = models.File.objects.filter(job=the_job)
            result["filenames"] = {}
            for jobFile in jobFiles:
                result["filenames"][jobFile.job_param_name] = str(jobFile.path)
            return result
        except AssertionError as err:
            logger.exception("Assertion Error in getJobInfo", exc_info=err)
        except Exception as err:
            logger.exception("Err in getJobInfo", exc_info=err)
        return None

    def _get_descendent_jobs(self, job):
        """
        Get descendent jobs for pipeline report generation.

        Returns a list of tuples where each tuple represents a child job.
        Format: [(child_job_id, [grandchild_ids]), ...]

        This format allows iteration like:
            for descendentJob in jobInfo["descendentjobs"]:
                child_job_id = descendentJob[0]
        """
        result = []
        children = list(job.children.all())
        for child in children:
            grandchildren = list(child.children.all())
            grandchild_ids = [str(gc.uuid) for gc in grandchildren]
            result.append((str(child.uuid), grandchild_ids))
            # Recursively add grandchildren
            result.extend(self._get_descendent_jobs(child))
        return result

    def jobDirectory(self, jobId=None, projectName=None, jobNumber=None, create=False, projectId=None, projectDirectory=None):
        logger.debug("in CCP4i2DjangoDbApi %s, %s, %s", jobId, projectName, jobNumber)
        # Accept projectId and projectDirectory for compatibility but pass to job_directory utility
        return job_directory(jobId=jobId, projectName=projectName, jobNumber=jobNumber, create=create,
                           projectId=projectId, projectDirectory=projectDirectory)

    def getJobFiles(
        self, jobId=None, role=0, mode="fileId", searchFileUses=True, fileTypes=[]
    ):
        """
        Retrieve files associated with a job with flexible filtering and output options.

        Args:
            jobId: The unique identifier of the job whose files you want to retrieve
            role: Specifies which files to get:
                FILE_ROLE_OUT (0): Output files produced by the job
                FILE_ROLE_IN (1): Input files used by the job
            mode: What information to return about the files:
                'fileId': Just the file IDs
                'fullPath': Complete file paths
                'fileName': Just the file names
                'annotation': File annotations/descriptions
                'fileType' or 'fileTypeId': File type information
                'all': Complete file information
            searchFileUses: Whether to also search the FileUses table (default: True)
            fileTypes: Optional filter to only return specific file types

        Returns:
            List of file information based on mode parameter
        """
        print("In CCP4i2DjangoDbApi getJobFiles")
        try:
            if jobId is None:
                logger.warning("getJobFiles called with no jobId")
                return []

            # Handle both UUID formats (with and without dashes)
            if "-" not in str(jobId):
                jobId = uuid.UUID(jobId)

            # Convert numeric file type IDs to mimetypes if needed
            converted_file_types = []
            if fileTypes:
                for file_type in fileTypes:
                    # Check if it's a numeric ID that needs conversion
                    if isinstance(file_type, (int, str)) and str(file_type).isdigit():
                        file_type_id = int(file_type)
                        # Find the mimetype for this file type ID
                        mimetype_found = None
                        for ft_id, ft_mimetype, ft_name in FILETYPELIST:
                            if ft_id == file_type_id:
                                mimetype_found = ft_mimetype
                                break
                        if mimetype_found is not None:
                            converted_file_types.append(mimetype_found)
                        else:
                            logger.warning(f"Unknown file type ID: {file_type_id}")
                    else:
                        # Already a mimetype string, use as-is
                        converted_file_types.append(file_type)
                fileTypes = converted_file_types

            # Get the job
            job = models.Job.objects.get(uuid=jobId)

            # Phase 1: Direct File Search
            files_from_direct = []
            if role == 0:  # FILE_ROLE_OUT - output files
                direct_files = models.File.objects.filter(job=job).exclude(
                    fileimport__isnull=False  # Exclude files that have import records
                )

                # Apply file type filter if specified
                if fileTypes:
                    direct_files = direct_files.filter(type__in=fileTypes)

                files_from_direct = list(direct_files)
            # Phase 2: FileUses Search (if enabled)
            files_from_uses = []
            if searchFileUses:
                # Search FileUses table for file usage relationships
                file_uses = models.FileUse.objects.filter(job=job, role=role)

                # Apply file type filter if specified
                if fileTypes:
                    file_uses = file_uses.filter(file__type__in=fileTypes)

                # Get the actual files from file uses
                files_from_uses = [file_use.file for file_use in file_uses]

            # Combine and deduplicate files
            all_files = files_from_direct + files_from_uses
            unique_files = []
            seen_uuids = set()
            for file_obj in all_files:
                if file_obj.uuid not in seen_uuids:
                    unique_files.append(file_obj)
                    seen_uuids.add(file_obj.uuid)
            print(f"Total unique files found: {len(unique_files)} {mode}")

            # Format output based on mode parameter
            if mode == "fileId":
                return [str(file_obj.uuid) for file_obj in unique_files]
            elif mode == "fullPath":
                return [str(file_obj.path) for file_obj in unique_files]
            elif mode == "fileName":
                return [file_obj.objectName() for file_obj in unique_files]
            elif mode == "jobparamname":
                return [file_obj.job_param_name for file_obj in unique_files]
            elif mode == "annotation":
                return [file_obj.annotation or "" for file_obj in unique_files]
            elif mode in ["fileType", "fileTypeId"]:
                return [
                    get_file_type_id_from_mimetype(file_obj.type)
                    for file_obj in unique_files
                ]
            elif mode == "all":
                # Return complete file information
                result = []
                for file_obj in unique_files:
                    file_info = {
                        "fileid": str(file_obj.uuid),
                        "filename": file_obj.objectName(),
                        "fullpath": str(file_obj.path),
                        "filetype": get_file_type_id_from_mimetype(file_obj.type),
                        "filetypename": get_file_type_name_from_mimetype(file_obj.type),
                        "subtype": file_obj.sub_type,
                        "annotation": file_obj.annotation or "",
                        "jobparamname": file_obj.job_param_name or "",
                        "content": file_obj.content,
                        # "directory": file_obj.directory,
                    }
                    result.append(file_info)
                return result
            else:
                # Default to returning file objects for unknown modes
                return unique_files

        except models.Job.DoesNotExist:
            logger.error("Job not found for jobId %s", jobId)
            return []
        except Exception as err:
            logger.exception("Error in getJobFiles for jobId %s", jobId, exc_info=err)
            return []

    def getFullPath(self, fileId=None):
        import clipper

        try:
            if fileId is None:
                logger.warning("getFullPath called with no fileId")
                return None

            # Handle both UUID formats (with and without dashes)
            if "-" not in str(fileId):
                fileId = uuid.UUID(fileId)

            # Get the file
            file = models.File.objects.get(uuid=fileId)
            print(f"File found: {file}, path: {str(file.path)}")
            return str(file.path)

        except models.File.DoesNotExist:
            logger.error("File not found for fileId %s", fileId)
            return None
        except Exception as err:
            logger.exception("Error in getFullPath for fileId %s", fileId, exc_info=err)
            return None

    def getChildJobs(
        self, jobId: str, descendents: bool = False, details: bool = False
    ):
        """
        Retrieve child jobs for a given parent job.

        Args:
            jobId: The UUID of the parent job
            descendents: If True, recursively get all descendant jobs
            details: If True, return detailed info; if False, return just UUIDs

        Returns:
            If details=False: List of child job UUIDs (without dashes)
            If details=True: List of tuples (job_number, uuid_without_dashes, task_name)
        """
        try:
            # Handle both UUID formats (with and without dashes)
            if "-" not in str(jobId):
                jobId = uuid.UUID(jobId)
            else:
                jobId = uuid.UUID(jobId)

            # Get the parent job
            parent_job = models.Job.objects.get(uuid=jobId)

            def get_children_recursive(job, include_descendants=False):
                """Helper function to get children, optionally recursive."""
                # Get direct children
                children = list(models.Job.objects.filter(parent=job))
                result = []

                for child in children:
                    if details:
                        # Return tuple with (number, uuid_without_dashes, task_name)
                        uuid_str = str(child.uuid).replace("-", "")
                        result.append((child.number, uuid_str, child.task_name))
                    else:
                        # Return just UUID without dashes
                        uuid_str = str(child.uuid).replace("-", "")
                        result.append(uuid_str)

                    # If descendents is True, recursively get children of this child
                    if include_descendants:
                        descendant_results = get_children_recursive(
                            child, include_descendants
                        )
                        result.extend(descendant_results)

                return result

            # Get the results
            return get_children_recursive(parent_job, descendents)

        except models.Job.DoesNotExist:
            logger.error("Parent job not found for jobId %s", jobId)
            return []
        except Exception as err:
            logger.exception("Error in getChildJobs for jobId %s", jobId, exc_info=err)
            return []
