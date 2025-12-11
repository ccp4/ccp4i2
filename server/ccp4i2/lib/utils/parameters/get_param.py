"""
Parameter retrieval utilities using CPluginScript architecture.

Provides functions to get job parameters with proper type handling.
Uses CPluginScript to ensure proper object hierarchy and access.
"""

import logging
from typing import Any, Dict

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


def get_parameter(
    job: models.Job,
    object_path: str
) -> Result[Dict[str, Any]]:
    """
    Get a parameter value using CPluginScript architecture.

    This function uses CPluginScript to ensure proper object hierarchy
    and access to parameters.

    Args:
        job: Job model instance
        object_path: Path to the parameter
                    Examples:
                    - "inputData.XYZIN" - Gets file object info
                    - "inputData.XYZIN.baseName" - Gets just filename
                    - "container.NCYCLES" - Gets control parameter

    Returns:
        Result[Dict] with parameter info including:
        - path: The object path requested
        - value: The parameter value
        - object_type: Type name of the object
        - file_path: (for CDataFile) Full file path
        - db_file_id: (for CDataFile) Database file UUID
        - base_name: (for CDataFile) File base name

    Example:
        >>> result = get_parameter(job, "inputData.XYZIN")
        >>> if result.success:
        ...     print(f"File path: {result.data['file_path']}")
        ...     print(f"Value: {result.data['value']}")
        >>>
        >>> result = get_parameter(job, "container.NCYCLES")
        >>> if result.success:
        ...     print(f"NCYCLES = {result.data['value']}")
    """
    # Get plugin with database context
    plugin_result = get_plugin_with_context(job)
    if not plugin_result.success:
        return Result.fail(
            f"Failed to load plugin: {plugin_result.error}",
            details=plugin_result.error_details
        )

    plugin = plugin_result.data

    try:
        # Navigate to the object
        # Normalize path to strip .container. segment if present from frontend
        normalized_path = normalize_object_path(object_path)

        logger.debug(
            "Getting parameter %s (normalized: %s) from job %s (task: %s)",
            object_path, normalized_path, job.uuid, job.task_name
        )

        obj = plugin.container
        parts = normalized_path.split('.')

        # Skip first part if it matches task name (legacy path format)
        if parts and parts[0] == job.task_name:
            parts = parts[1:]

        for part in parts:
            if hasattr(obj, part):
                obj = getattr(obj, part)
            else:
                return Result.fail(
                    f"Parameter path '{object_path}' not found - missing attribute '{part}'",
                    details={
                        "job_id": str(job.uuid),
                        "task_name": job.task_name,
                        "object_path": object_path,
                        "missing_part": part
                    }
                )

        # Build result data
        result_data = {
            "path": object_path,
            "object_type": type(obj).__name__,
        }

        # Get value (handle different CData types)
        if hasattr(obj, 'value'):
            result_data["value"] = obj.value
        elif hasattr(obj, 'get'):
            result_data["value"] = obj.get()
        else:
            result_data["value"] = str(obj)

        # Add file-specific info if it's a CDataFile
        if hasattr(obj, 'getFullPath'):
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

        logger.debug(
            "Successfully retrieved parameter %s from job %s: %s",
            object_path, job.uuid, result_data.get('value')
        )
        return Result.ok(result_data)

    except AttributeError as e:
        logger.error(
            "Parameter path '%s' not found on job %s: %s",
            object_path, job.uuid, str(e)
        )
        return Result.fail(
            f"Parameter path '{object_path}' not found",
            details={
                "job_id": str(job.uuid),
                "task_name": job.task_name,
                "object_path": object_path,
                "error": str(e)
            }
        )

    except Exception as e:
        logger.exception(
            "Failed to get parameter %s from job %s",
            object_path, job.uuid
        )
        return Result.fail(
            f"Error getting parameter: {str(e)}",
            details={
                "job_id": str(job.uuid),
                "task_name": job.task_name,
                "object_path": object_path,
                "error_type": type(e).__name__
            }
        )
