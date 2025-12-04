"""
Parameter setting utilities using CPluginScript + dbHandler architecture.

Provides functions to set job parameters with proper type handling and persistence.
Uses CPluginScript to ensure proper file handling, database synchronization,
and validation.
"""

import logging
import json
from typing import Union, Any, Dict
from pathlib import Path

from ccp4x.db import models
from ccp4x.lib.response import Result
from ccp4x.lib.utils.plugins.plugin_context import get_plugin_with_context

logger = logging.getLogger(__name__)


def normalize_object_path(object_path: str) -> str:
    """
    Normalize object paths from frontend to match backend container structure.

    The frontend JSON encoder includes the full hierarchy path which includes
    `.container.` (e.g., `prosmart_refmac.container.inputData.XYZIN`), but the
    backend container structure doesn't have that extra level.

    This function strips the `.container.` segment if present after the task name.

    Args:
        object_path: Path like "prosmart_refmac.container.inputData.XYZIN"

    Returns:
        Normalized path like "prosmart_refmac.inputData.XYZIN"
    """
    # Split into parts
    parts = object_path.split('.')

    # If second element is 'container', remove it
    # e.g., ['prosmart_refmac', 'container', 'inputData', 'XYZIN']
    #    -> ['prosmart_refmac', 'inputData', 'XYZIN']
    if len(parts) >= 2 and parts[1] == 'container':
        parts = [parts[0]] + parts[2:]

    return '.'.join(parts)


def set_parameter(
    job: models.Job,
    object_path: str,
    value: Union[str, int, float, bool, dict, None]
) -> Result[Dict[str, Any]]:
    """
    Set a parameter value using CPluginScript + dbHandler architecture.

    **IMPORTANT**: This function only works on jobs in PENDING/UNKNOWN status.
    For finished jobs, use the clone API first to create an editable copy.

    This function uses CPluginScript to ensure:
    - Proper file handling (CDataFile.setFullPath() with DB awareness)
    - Database synchronization (dbHandler.updateJobStatus())
    - Correct object hierarchy
    - Validation support

    Args:
        job: Job model instance (must be PENDING or UNKNOWN status)
        object_path: Path to the parameter
                    Examples:
                    - "inputData.XYZIN" - Sets entire file object
                    - "inputData.XYZIN.baseName" - Sets just filename
                    - "container.NCYCLES" - Sets control parameter
        value: New value for the parameter

    Returns:
        Result[Dict] with parameter info or error

    Example:
        >>> # Set parameter on pending job
        >>> result = set_parameter(job, "inputData.XYZIN", "/path/to/file.pdb")
        >>> if result.success:
        ...     print(f"File path: {result.data['file_path']}")
        >>>
        >>> # For finished jobs, clone first
        >>> clone_result = clone_job(finished_job)
        >>> new_job = clone_result.data
        >>> result = set_parameter(new_job, "container.NCYCLES", 10)
    """
    # Validate job status - only allow parameter setting on pending jobs
    if job.status not in [models.Job.Status.UNKNOWN, models.Job.Status.PENDING]:
        return Result.fail(
            f"Cannot modify parameters on job with status '{job.status}'. "
            f"Use clone API to create an editable copy first.",
            details={
                "job_id": str(job.uuid),
                "job_status": job.status,
                "allowed_statuses": [models.Job.Status.UNKNOWN, models.Job.Status.PENDING]
            }
        )

    # Get plugin with database context
    plugin_result = get_plugin_with_context(job)
    if not plugin_result.success:
        return Result.fail(
            f"Failed to load plugin: {plugin_result.error}",
            details=plugin_result.error_details
        )

    plugin = plugin_result.data

    try:
        # Normalize path to strip .container. segment if present from frontend
        normalized_path = normalize_object_path(object_path)

        # Set parameter through plugin's container using modern context-aware method
        # This ensures proper file handling, validation, hierarchy, and database sync
        logger.debug(
            "Setting parameter %s (normalized: %s) = %s on job %s (task: %s)",
            object_path, normalized_path, value, job.uuid, job.task_name
        )

        # Use modern CContainer.set_parameter() which auto-detects CPluginScript parent
        # and enables database synchronization when appropriate
        obj = plugin.container.set_parameter(normalized_path, value, skip_first=True)

        # Save parameters to input_params.xml (user control stage)
        # Use CPluginScript.saveDataToXml which uses ParamsXmlHandler for proper filtering
        # This is different from params.xml which is written at plugin lifecycle stages:
        # - After checkInputData() - clears unpopulated/non-existent file inputs
        # - After checkOutputData() - includes candidate output file names
        # - After processOutputFiles() - weeds out non-existent output files
        input_params_file = job.directory / "input_params.xml"
        logger.debug("Saving parameters to %s", input_params_file)
        error = plugin.saveDataToXml(str(input_params_file))
        if error and hasattr(error, 'hasError') and error.hasError():
            logger.error("Failed to save parameters to %s: %s", input_params_file, error)
        else:
            logger.debug("Successfully saved parameters to %s", input_params_file)

        # Note: Database sync already handled by container.set_parameter() if in CPluginScript context
        # No need to manually call dbHandler.updateJobStatus() here

        # obj now contains the CData object that was set

        # Build result data
        result_data = {
            "path": object_path,
            "value": value,
            "object_type": type(obj).__name__ if obj else "Unknown",
        }

        # Add file-specific info if it's a CDataFile
        if obj and hasattr(obj, 'getFullPath'):
            full_path = obj.getFullPath()
            if full_path:
                result_data["file_path"] = full_path

            # Get dbFileId if available
            if hasattr(obj, 'dbFileId'):
                db_file_attr = getattr(obj, 'dbFileId')
                if hasattr(db_file_attr, 'value') and db_file_attr.value:
                    result_data["db_file_id"] = str(db_file_attr.value)

            # Get baseName if available
            if hasattr(obj, 'baseName'):
                base_name_attr = getattr(obj, 'baseName')
                if hasattr(base_name_attr, 'value'):
                    result_data["base_name"] = str(base_name_attr.value)
                else:
                    result_data["base_name"] = str(base_name_attr)

        logger.info(
            "Successfully set parameter %s (normalized: %s) on job %s: %s",
            object_path, normalized_path, job.uuid, result_data.get('file_path', value)
        )
        return Result.ok(result_data)

    except AttributeError as e:
        logger.error(
            "Parameter path '%s' (normalized: '%s') not found on job %s: %s",
            object_path, normalized_path, job.uuid, str(e)
        )
        return Result.fail(
            f"Parameter path '{normalized_path}' not found",
            details={
                "job_id": str(job.uuid),
                "task_name": job.task_name,
                "object_path": object_path,
                "normalized_path": normalized_path,
                "error": str(e)
            }
        )

    except Exception as e:
        logger.exception(
            "Failed to set parameter %s (normalized: %s) on job %s",
            object_path, normalized_path, job.uuid
        )
        return Result.fail(
            f"Error setting parameter: {str(e)}",
            details={
                "job_id": str(job.uuid),
                "task_name": job.task_name,
                "object_path": object_path,
                "normalized_path": normalized_path,
                "value": value,
                "error_type": type(e).__name__
            }
        )
