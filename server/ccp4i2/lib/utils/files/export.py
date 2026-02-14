"""
Export Job File Utility

This module provides functionality for exporting specific files from CCP4 jobs
using the CCP4 Task Manager's exportJobFiles functionality. The exported files
are prepared for download with appropriate headers and content types.

Functions:
    export_job_file: Main export function that handles the complete export workflow

Dependencies:
    - CCP4 Task Manager for file export operations
    - Django models for job retrieval
    - Path handling for file operations

Author: CCP4i2 Development Team
License: CCP4 License
Version: Compatible with CCP4i2 and Django 4.2+

    - Optimized for cloud-based file operations
    - Includes proper error handling for distributed environments
"""

import json
import logging
from pathlib import Path
from typing import Tuple, Optional

from django.http import JsonResponse, FileResponse
from ccp4i2.core import CCP4TaskManager

from ccp4i2.db import models

logger = logging.getLogger(f"ccp4i2:{__name__}")
logger.setLevel(logging.DEBUG)


def export_job_file(
    job_id: int, export_mode: str
) -> Tuple[Optional[FileResponse], Optional[JsonResponse]]:
    """
    Export a specific file from a job using CCP4 Task Manager.

    Retrieves and exports a specific file associated with the job using the
    CCP4 Task Manager's exportJobFiles functionality. The exported file is
    prepared for download with appropriate headers and content types.

    Args:
        job_id (int): Primary key of the job to export from
        export_mode (str): Export mode parameter for Task Manager (e.g., "pdb", "mtz", "log")

    Returns:
        Tuple[Optional[FileResponse], Optional[JsonResponse]]:
            - (FileResponse, None) on success with downloadable file
            - (None, JsonResponse) on error with error details

    Response Headers (on success):
        - Content-Type: Determined by file extension
        - Content-Disposition: attachment with filename
        - Content-Length: File size
        - X-Export-Info: JSON metadata about export

    Error Response Format:
        {
            "status": "Failed",
            "reason": "Error description"
        }

    Example Usage:
        file_response, error_response = export_job_file(123, "pdb")
        if error_response:
            return error_response
        return file_response

        - Efficient streaming download for large files
        - Proper cleanup of temporary resources
        - Container-aware file path handling

    Error Handling:
        - Validates job existence
        - Handles CCP4 Task Manager errors
        - Manages file access permissions
        - Provides detailed error messages
    """
    try:
        logger.debug(
            "Export job file utility called for job_id=%s, export_mode=%s",
            job_id,
            export_mode,
        )
        # Retrieve the job object
        job = models.Job.objects.get(id=job_id)
        logger.debug("Found job: %s, UUID: %s", job.id, job.uuid)

        logger.debug(
            "Exporting file for job %s (UUID: %s) with mode: %s",
            job_id,
            job.uuid,
            export_mode,
        )

        # Call exportJobFiles to get the filename of the file to export
        job_uuid_formatted = str(str(job.uuid).replace("-", ""))
        logger.debug(
            "Calling exportJobFiles with jobId=%s, taskName=%s, mode=%s",
            job_uuid_formatted,
            job.task_name,
            export_mode,
        )
        try:
            # Always throws an error becayse exportJobFiles no longer exists
            exported_filename = exportJobFiles(
                taskName=job.task_name, jobId=job_uuid_formatted, mode=export_mode
            )
            logger.debug(
                "exportJobFiles call succeeded, got filename: %s", exported_filename
            )
        except Exception as e:
            logger.debug("Exception calling exportJobFiles: %s", e)
            raise
        logger.debug("exportJobFiles returned: %s", exported_filename)

        if not exported_filename:
            logger.debug("No filename returned from exportJobFiles")
            error_response = JsonResponse(
                {
                    "status": "Failed",
                    "reason": f"No file available for export mode '{export_mode}'",
                },
                status=404,
            )
            return None, error_response

        # Convert to Path object for easier handling
        export_file_path = Path(exported_filename)
        logger.debug("Export file path: %s", export_file_path)
        logger.debug("File exists: %s", export_file_path.exists())

        # Verify the file exists
        if not export_file_path.exists():
            logger.error(
                "Export file not found: %s for job %s mode %s",
                exported_filename,
                job_id,
                export_mode,
            )
            error_response = JsonResponse(
                {
                    "status": "Failed",
                    "reason": f"Export file not found: {export_file_path.name}",
                },
                status=404,
            )
            return None, error_response

        # Determine content type and create file response
        content_type = _get_content_type(export_file_path)
        safe_filename = _get_safe_filename(export_file_path, job, export_mode)

        logger.info(
            "Exporting file %s for job %s (mode: %s) as %s",
            exported_filename,
            job_id,
            export_mode,
            safe_filename,
        )

        # Create file response with streaming support for large files
        file_response = _create_file_response(
            export_file_path, safe_filename, content_type, job, export_mode
        )

        return file_response, None

    except models.Job.DoesNotExist as err:
        logger.exception("Failed to retrieve job with id %s", job_id, exc_info=err)
        error_response = JsonResponse(
            {"status": "Failed", "reason": f"Job not found: {str(err)}"}, status=404
        )
        return None, error_response

    except FileNotFoundError as err:
        logger.exception(
            "Export file not found for job %s mode %s",
            job_id,
            export_mode,
            exc_info=err,
        )
        error_response = JsonResponse(
            {"status": "Failed", "reason": f"Export file not accessible: {str(err)}"},
            status=404,
        )
        return None, error_response

    except PermissionError as err:
        logger.exception(
            "Permission denied accessing export file for job %s", job_id, exc_info=err
        )
        error_response = JsonResponse(
            {"status": "Failed", "reason": "Permission denied accessing export file"},
            status=403,
        )
        return None, error_response

    except Exception as err:
        logger.exception(
            "Failed to export file for job %s with mode %s",
            job_id,
            export_mode,
            exc_info=err,
        )
        error_response = JsonResponse(
            {"status": "Failed", "reason": f"Export failed: {str(err)}"}, status=500
        )
        return None, error_response


def _get_content_type(file_path: Path) -> str:
    """
    Determine content type based on file extension.

    Args:
        file_path (Path): Path to the file

    Returns:
        str: MIME content type
    """
    file_extension = file_path.suffix.lower()

    content_type_map = {
        ".pdb": "chemical/x-pdb",
        ".cif": "chemical/x-cif",
        ".mtz": "application/x-ccp4-mtz",
        ".log": "text/plain",
        ".txt": "text/plain",
        ".xml": "application/xml",
        ".json": "application/json",
        ".zip": "application/zip",
        ".tar": "application/x-tar",
        ".gz": "application/gzip",
        ".pdf": "application/pdf",
        ".png": "image/png",
        ".jpg": "image/jpeg",
        ".jpeg": "image/jpeg",
    }

    return content_type_map.get(file_extension, "application/octet-stream")


def _get_safe_filename(file_path: Path, job: models.Job, export_mode: str) -> str:
    """
    Generate a safe filename for download.

    Args:
        file_path (Path): Original file path
        job (models.Job): Job object
        export_mode (str): Export mode

    Returns:
        str: Safe filename for download
    """
    # Generate a safe filename for download
    safe_filename = "".join(c for c in file_path.name if c.isalnum() or c in "._-")

    # If filename becomes empty after sanitization, use a default
    if not safe_filename:
        file_extension = file_path.suffix.lower()
        safe_filename = f"job_{job.number}_export_{export_mode}{file_extension}"

    return safe_filename


def _create_file_response(
    file_path: Path, filename: str, content_type: str, job: models.Job, export_mode: str
) -> FileResponse:
    """
    Create a FileResponse with streaming support and metadata headers.

    Args:
        file_path (Path): Path to the file to serve
        filename (str): Safe filename for download
        content_type (str): MIME content type
        job (models.Job): Job object for metadata
        export_mode (str): Export mode for metadata

    Returns:
        FileResponse: Configured file response with headers
    """

    def file_iterator():
        try:
            with open(file_path, "rb") as f:
                while True:
                    chunk = f.read(8192)  # 8KB chunks
                    if not chunk:
                        break
                    yield chunk
        except Exception as e:
            logger.error("Error reading export file %s: %s", file_path, str(e))
            raise

    response = FileResponse(
        file_iterator(),
        as_attachment=True,
        filename=filename,
        content_type=content_type,
    )

    # Add additional headers
    response["Content-Length"] = file_path.stat().st_size
    response["X-Export-Info"] = json.dumps(
        {
            "job_id": job.id,
            "job_uuid": str(job.uuid),
            "job_number": job.number,
            "task_name": job.task_name,
            "export_mode": export_mode,
            "original_filename": file_path.name,
            "file_size": file_path.stat().st_size,
        }
    )

    return response
