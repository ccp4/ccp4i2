"""
Resolve fileUse references to file metadata.

This module provides functionality to resolve fileUse strings (references to
files from previous jobs) into file metadata dictionaries that can be used
to set file parameters.

FileUse Syntax:
    task_name[jobIndex].jobParamName[paramIndex]  - Full format
    [jobIndex].jobParamName[paramIndex]           - Without task name
    task_name[jobIndex].jobParamName              - Without param index
    [jobIndex].jobParamName                       - Simplest format

Examples:
    [-1].XYZOUT[0]           - First XYZOUT from most recent job
    prosmart_refmac[-1].XYZOUT  - XYZOUT from most recent prosmart_refmac job
    refmac[-2].HKLOUT[0]     - HKLOUT from second-to-last refmac job

The jobIndex can be negative (counting from end) or positive (counting from start).
"""

import re
import logging
from typing import Optional
from pathlib import Path

from ....db import models
from ...response import Result

logger = logging.getLogger(f"ccp4i2:{__name__}")


# Regex patterns for parsing fileUse strings
# Order matters - try more specific patterns first
FILEUSE_PATTERNS = [
    # task_name[jobIndex].jobParamName[paramIndex]
    re.compile(
        r"^(?P<task_name>\w+)\[(?P<jobIndex>-?\d+)\]\.(?P<jobParamName>\w+)\[(?P<paramIndex>\d+)\]$"
    ),
    # [jobIndex].jobParamName[paramIndex] (no task name)
    re.compile(
        r"^\[(?P<jobIndex>-?\d+)\]\.(?P<jobParamName>\w+)\[(?P<paramIndex>\d+)\]$"
    ),
    # task_name[jobIndex].jobParamName (no param index)
    re.compile(
        r"^(?P<task_name>\w+)\[(?P<jobIndex>-?\d+)\]\.(?P<jobParamName>\w+)$"
    ),
    # [jobIndex].jobParamName (simplest - no task name, no param index)
    re.compile(
        r"^\[(?P<jobIndex>-?\d+)\]\.(?P<jobParamName>\w+)$"
    ),
]


def is_fileuse_pattern(value: str) -> bool:
    """
    Check if a string matches the fileUse pattern.

    This is used for auto-detection of fileUse syntax without the explicit
    'fileUse=' prefix.

    Args:
        value: String to check

    Returns:
        True if the string matches a fileUse pattern
    """
    if not value or not isinstance(value, str):
        return False

    # Quick pre-check: must contain [ and ] and .
    if '[' not in value or ']' not in value or '.' not in value:
        return False

    # Try each pattern
    for pattern in FILEUSE_PATTERNS:
        if pattern.match(value):
            return True

    return False


def parse_fileuse(fileuse: str) -> dict:
    """
    Parse a fileUse string into its components.

    Args:
        fileuse: FileUse string to parse

    Returns:
        Dict with keys: task_name, jobIndex, jobParamName, paramIndex

    Raises:
        ValueError: If the string doesn't match any fileUse pattern
    """
    # Strip fileUse= prefix if present
    if fileuse.startswith("fileUse="):
        fileuse = fileuse[8:]

    # Strip quotes if present
    if fileuse.startswith('"') and fileuse.endswith('"'):
        fileuse = fileuse[1:-1]
    if fileuse.startswith("'") and fileuse.endswith("'"):
        fileuse = fileuse[1:-1]

    for pattern in FILEUSE_PATTERNS:
        match = pattern.match(fileuse)
        if match:
            result = {
                "task_name": match.groupdict().get("task_name"),
                "jobIndex": int(match.group("jobIndex")),
                "jobParamName": match.group("jobParamName"),
                "paramIndex": int(match.groupdict().get("paramIndex", 0)),
            }
            logger.debug(f"Parsed fileUse '{fileuse}' -> {result}")
            return result

    raise ValueError(
        f"Invalid fileUse syntax: '{fileuse}'. "
        f"Expected format like '[-1].XYZOUT[0]' or 'task_name[-1].PARAM'"
    )


def resolve_fileuse(
    project: models.Project,
    fileuse: str,
) -> Result[dict]:
    """
    Resolve a fileUse string to file metadata.

    Args:
        project: Project to search within
        fileuse: FileUse string (e.g., "[-1].XYZOUT[0]")

    Returns:
        Result containing file metadata dict on success:
            {
                "project": str (UUID without hyphens),
                "baseName": str (filename),
                "dbFileId": str (file UUID without hyphens),
                "relPath": str (relative path),
                "fullPath": str (full filesystem path),
            }
    """
    try:
        parsed = parse_fileuse(fileuse)
    except ValueError as e:
        return Result.failure(str(e))

    task_name = parsed["task_name"]
    job_index = parsed["jobIndex"]
    job_param_name = parsed["jobParamName"]
    param_index = parsed["paramIndex"]

    # Build job query
    jobs_query = models.Job.objects.filter(project=project).select_related('project')

    if task_name:
        jobs_query = jobs_query.filter(task_name=task_name)

    # Order by creation time/number to support indexing
    jobs_query = jobs_query.order_by('number')

    job_count = jobs_query.count()

    if job_count == 0:
        if task_name:
            return Result.failure(
                f"No jobs with task_name='{task_name}' found in project '{project.name}'"
            )
        return Result.failure(f"No jobs found in project '{project.name}'")

    # Handle negative indexing
    if job_index < 0:
        actual_index = job_count + job_index
    else:
        actual_index = job_index

    if actual_index < 0 or actual_index >= job_count:
        return Result.failure(
            f"Job index {job_index} out of range. "
            f"Project has {job_count} {'matching ' if task_name else ''}jobs (indices 0 to {job_count-1}, or -1 to -{job_count})"
        )

    try:
        the_job = jobs_query[actual_index]
    except IndexError:
        return Result.failure(f"Could not retrieve job at index {job_index}")

    logger.info(f"Resolved job index {job_index} to job #{the_job.number} ({the_job.task_name})")

    # Try output files first
    output_files = models.File.objects.filter(
        job=the_job,
        job_param_name=job_param_name
    ).select_related('job__project')

    if output_files.exists():
        file_count = output_files.count()
        if param_index >= file_count:
            return Result.failure(
                f"Parameter index {param_index} out of range. "
                f"Job #{the_job.number} has {file_count} files with param_name='{job_param_name}'"
            )
        the_file = output_files[param_index]
    else:
        # Try input files via FileUse
        input_file_uses = models.FileUse.objects.filter(
            job=the_job,
            job_param_name=job_param_name
        ).select_related('file', 'file__job__project')

        if not input_file_uses.exists():
            # List available param names to help user
            available_params = set()
            for f in models.File.objects.filter(job=the_job).values_list('job_param_name', flat=True):
                if f:
                    available_params.add(f)
            for fu in models.FileUse.objects.filter(job=the_job).values_list('job_param_name', flat=True):
                if fu:
                    available_params.add(fu)

            available_str = ", ".join(sorted(available_params)) if available_params else "(none)"
            return Result.failure(
                f"No file with param_name='{job_param_name}' found in job #{the_job.number}. "
                f"Available param names: {available_str}"
            )

        file_count = input_file_uses.count()
        if param_index >= file_count:
            return Result.failure(
                f"Parameter index {param_index} out of range. "
                f"Job #{the_job.number} has {file_count} input files with param_name='{job_param_name}'"
            )

        the_file = input_file_uses[param_index].file

    # Build file metadata dictionary
    file_dict = {
        "project": str(the_file.job.project.uuid).replace("-", ""),
        "baseName": the_file.name,
        "dbFileId": str(the_file.uuid).replace("-", ""),
    }

    # Determine relPath based on file directory type
    if the_file.directory == models.File.Directory.IMPORT_DIR:
        file_dict["relPath"] = "CCP4_IMPORTED_FILES"
        full_path = Path(the_file.job.project.directory) / "CCP4_IMPORTED_FILES" / the_file.name
    else:
        # For job files, build path from job directory
        job_dir = Path(the_file.job.directory)
        project_dir = Path(the_file.job.project.directory)

        # Get relative path from project
        try:
            rel_path = job_dir.relative_to(project_dir)
            file_dict["relPath"] = str(rel_path)
        except ValueError:
            # Fallback if not relative
            if "CCP4_JOBS" in job_dir.parts:
                parts = job_dir.parts
                ccp4_jobs_index = parts.index("CCP4_JOBS")
                file_dict["relPath"] = str(Path(*parts[ccp4_jobs_index:]))
            else:
                file_dict["relPath"] = str(job_dir.name)

        full_path = job_dir / the_file.name

    file_dict["fullPath"] = str(full_path)

    logger.info(f"Resolved fileUse '{fileuse}' to file: {file_dict}")

    return Result.success(file_dict)


def resolve_fileuse_by_project_id(
    project_id: str,
    fileuse: str,
) -> Result[dict]:
    """
    Resolve a fileUse string given a project ID (UUID or name).

    Args:
        project_id: Project UUID or name
        fileuse: FileUse string

    Returns:
        Result containing file metadata dict on success
    """
    # Try to find project by UUID first
    try:
        # Handle UUID with or without hyphens
        clean_uuid = project_id.replace("-", "")
        if len(clean_uuid) == 32:
            # Looks like a UUID
            project = models.Project.objects.get(uuid__icontains=clean_uuid[:8])
        else:
            raise models.Project.DoesNotExist()
    except models.Project.DoesNotExist:
        # Try by name
        try:
            project = models.Project.objects.get(name=project_id)
        except models.Project.DoesNotExist:
            # Try case-insensitive
            try:
                project = models.Project.objects.get(name__iexact=project_id)
            except models.Project.DoesNotExist:
                return Result.failure(f"Project not found: '{project_id}'")

    return resolve_fileuse(project, fileuse)
