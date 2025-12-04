"""
Job export utilities.

Provides functions to export jobs as ZIP archives.
"""

import logging
from pathlib import Path
from ccp4x.db import models
from ccp4x.lib.response import Result
from ccp4x.db.export_project import export_project_to_zip

logger = logging.getLogger(f"ccp4x:{__name__}")


def export_job(job: models.Job, output_path: Path) -> Result[Path]:
    """
    Export a single job as a ZIP archive.

    Creates a downloadable ZIP archive containing the job data,
    associated files, and dependency information.

    Args:
        job: Job model instance to export
        output_path: Path where ZIP file should be created

    Returns:
        Result containing path to created ZIP file

    Example:
        >>> result = export_job(job, Path("/tmp/job_export.zip"))
        >>> if result.success:
        ...     print(f"Exported to: {result.data}")
    """
    try:
        project = job.project
        job_selection = {str(job.number)}

        # Use existing export functionality
        result_path = export_project_to_zip(
            project=project,
            output_path=output_path,
            job_selection=job_selection
        )

        logger.info(
            "Exported job %s (number: %s) from project %s to %s",
            job.uuid, job.number, project.name, result_path
        )

        return Result.ok(result_path)

    except Exception as err:
        logger.exception("Failed to export job %s", job.uuid, exc_info=err)
        return Result.fail(
            f"Failed to export job: {str(err)}",
            details={
                "job_id": str(job.uuid),
                "job_number": job.number,
                "output_path": str(output_path),
                "error_type": type(err).__name__
            }
        )
