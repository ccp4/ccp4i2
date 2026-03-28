# Copyright (C) 2025-2026 University of York
# Copyright (C) 2025-2026 Newcastle University
#
# This file is part of CCP4i2.
#
# CCP4i2 is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License version 3,
# modified in accordance with the provisions of the license to address
# the requirements of UK law.
#
# See https://www.ccp4.ac.uk/ccp4license.php for details.
"""
Unified utility for getting CPluginScript instances with database context.

This is the CANONICAL way to get a plugin for any job operation.
All other utilities should use this rather than creating containers directly.
"""
import logging
from typing import Optional
from pathlib import Path

from ccp4i2.db.models import Job
from ccp4i2.db.async_db_handler import AsyncDatabaseHandler
from .get_plugin import get_job_plugin
from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.lib.response import Result

logger = logging.getLogger(__name__)


def get_plugin_with_context(
    job: Job,
    params_file: Optional[Path] = None,
    create_db_handler: bool = True
) -> Result[CPluginScript]:
    """
    Get CPluginScript instance with full database context.

    This is the CANONICAL way to get a plugin for any operation.
    It ensures:
    - Plugin has dbHandler attached for DB synchronization
    - Plugin knows its job UUID, job number, project context
    - Container is loaded from params.xml or input_params.xml
    - File operations have proper DB awareness

    Args:
        job: Django Job model instance
        params_file: Optional specific params file to load (otherwise auto-detects)
        create_db_handler: Whether to attach dbHandler (default True)

    Returns:
        Result[CPluginScript] with plugin instance or error

    Example:
        # Get plugin for parameter operations
        result = get_plugin_with_context(job)
        if result.success:
            plugin = result.data

            # Set parameter through plugin's container
            from ccp4i2.lib.utils.parameters.set_parameter import set_parameter
            set_parameter(job, "inputData.XYZIN", "/path/to/file.pdb")

            # Save parameters
            plugin.container.saveDataToXml(str(job.directory / "params.xml"))

            # Update database
            if plugin._dbHandler:
                plugin._dbHandler.updateJobStatus(
                    jobId=str(job.uuid),
                    container=plugin.container
                )
    """
    try:
        # Create dbHandler for database synchronization
        dbHandler = None
        if create_db_handler:
            logger.debug(
                "Creating dbHandler for job %s (project: %s)",
                job.uuid, job.project.uuid
            )
            dbHandler = AsyncDatabaseHandler(project_uuid=job.project.uuid)

        # Get plugin instance (automatically loads params.xml or input_params.xml)
        logger.debug(
            "Loading plugin for task '%s' (job: %s)",
            job.task_name, job.uuid
        )

        plugin = get_job_plugin(
            the_job=job,
            parent=None,
            dbHandler=dbHandler
        )

        if plugin is None:
            return Result.fail(
                f"Failed to load plugin for task '{job.task_name}'",
                details={
                    "job_id": str(job.uuid),
                    "task_name": job.task_name,
                    "job_directory": str(job.directory)
                }
            )

        # Set database context on plugin
        # This allows file operations to know they're in a DB-aware environment
        plugin.set_db_job_id(str(job.uuid))
        plugin.set_db_job_number(job.number)
        plugin._dbProjectId = str(job.project.uuid)
        plugin._dbProjectName = job.project.name

        logger.info(
            "Loaded plugin '%s' with DB context: job=%s (#%s), project=%s",
            job.task_name, job.uuid, job.number, job.project.name
        )

        return Result.ok(plugin)

    except FileNotFoundError as e:
        logger.error(
            "Params file not found for job %s: %s",
            job.uuid, str(e)
        )
        return Result.fail(
            f"Params file not found: {str(e)}",
            details={
                "job_id": str(job.uuid),
                "task_name": job.task_name,
                "directory": str(job.directory),
                "error_type": "FileNotFoundError"
            }
        )

    except Exception as e:
        logger.exception(
            "Failed to get plugin with context for job %s (task: %s)",
            job.uuid, job.task_name
        )
        return Result.fail(
            f"Error loading plugin: {str(e)}",
            details={
                "job_id": str(job.uuid),
                "task_name": job.task_name,
                "error_type": type(e).__name__,
                "error_message": str(e)
            }
        )
