"""
Job parameter validation utilities using CPluginScript architecture.

Provides validation of job parameters against task definitions,
returning detailed error reports. Uses CPluginScript.validity_as_xml()
for comprehensive validation including container and input file checks.
"""

import logging
from xml.etree import ElementTree as ET
from ccp4x.db import models
from ccp4x.lib.response import Result
from ccp4x.lib.utils.plugins.plugin_context import get_plugin_with_context

logger = logging.getLogger(__name__)


def validate_job(job: models.Job) -> Result[ET.Element]:
    """
    Validate job parameters using CPluginScript architecture.

    Uses CPluginScript.validity_as_xml() which provides comprehensive
    validation including:
    - Container validation (all child objects' validity() methods)
    - Input file existence checks (mustExist files)

    Args:
        job: Job model instance to validate

    Returns:
        Result containing XML error report element with errors/warnings

    Example:
        >>> result = validate_job(job)
        >>> if result.success:
        ...     error_tree = result.data
        ...     errors = error_tree.findall('.//errorReport')
        ...     print(f"Found {len(errors)} validation errors")
        ...
        ...     # Check severity
        ...     for error in errors:
        ...         severity = error.findtext('severity')
        ...         description = error.findtext('description')
        ...         print(f"{severity}: {description}")
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
        logger.debug(
            "Validating job %s (task: %s) using plugin.validity_as_xml()",
            job.uuid, job.task_name
        )

        # Use CPluginScript.validity_as_xml() for comprehensive validation
        # This includes container validation AND input file checks
        error_tree = plugin.validity_as_xml()

        # Count errors for logging
        error_reports = error_tree.findall('.//errorReport')
        logger.info(
            "Validated job %s (task: %s): found %d validation reports",
            job.uuid, job.task_name, len(error_reports)
        )

        return Result.ok(error_tree)

    except AttributeError as err:
        # This can happen if container.validity() or validate_container has issues
        logger.exception(
            "Validation error for job %s - possible API mismatch: %s",
            job.uuid, str(err)
        )
        return Result.fail(
            f"Validation failed - API error: {str(err)}",
            details={
                "job_id": str(job.uuid),
                "task_name": job.task_name,
                "error_type": "AttributeError",
                "error_message": str(err)
            }
        )

    except FileNotFoundError as err:
        logger.error("Job params file not found for job %s: %s", job.uuid, str(err))
        return Result.fail(
            f"Job parameters not found: {str(err)}",
            details={
                "job_id": str(job.uuid),
                "task_name": job.task_name,
                "error_type": "FileNotFoundError"
            }
        )

    except Exception as err:
        logger.exception("Failed to validate job %s", job.uuid)
        return Result.fail(
            f"Validation failed: {str(err)}",
            details={
                "job_id": str(job.uuid),
                "task_name": job.task_name,
                "error_type": type(err).__name__,
                "error_message": str(err)
            }
        )
