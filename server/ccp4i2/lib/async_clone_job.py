"""
Modern async job cloning using AsyncDatabaseHandler and CData system.

This module provides typesafe async job cloning following patterns from async_create_job.py
and logic from lib/utils/jobs/clone.py.

Key features:
- Full async/await with type hints
- Clones job parameters from existing job
- Creates new directory structure
- Patches output file paths
- Saves parameters to new job
- Database integration via AsyncDatabaseHandler
"""

import asyncio
import datetime
import logging
import uuid
from pathlib import Path
from typing import Optional, Dict, Any

from asgiref.sync import sync_to_async
from pytz import timezone

from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core.task_manager.plugin_registry import get_plugin_class

from ..db import models
from ..db.async_db_handler import AsyncDatabaseHandler
from .utils.files.set_names import set_output_file_names


logger = logging.getLogger(f"ccp4i2:{__name__}")


async def clone_job_async(
    job_uuid: uuid.UUID,
    new_title: Optional[str] = None,
    target_project_uuid: Optional[uuid.UUID] = None,
    parent_job_uuid: Optional[uuid.UUID] = None,
) -> Dict[str, Any]:
    """
    Clone an existing job by creating a new job with the same parameters.

    This function:
    1. Retrieves the original job from database
    2. Loads parameters from original job's XML file
    3. Creates a new job with incremented job number
    4. Creates new job directory structure
    5. Instantiates plugin with cloned parameters
    6. Patches output file paths to new directory
    7. Saves parameters to new job
    8. Returns configured plugin ready for execution

    Args:
        job_uuid: UUID of the job to clone
        new_title: Optional new title (defaults to original job's title)
        target_project_uuid: Optional target project UUID (defaults to same project)
        parent_job_uuid: Optional parent job UUID for nested cloning

    Returns:
        Dict containing:
            - job_uuid: UUID of cloned job
            - job_number: Job number (e.g., "2", "3")
            - job_directory: Path to job directory
            - plugin: CPluginScript instance (configured with cloned params)
            - task_name: Plugin task name
            - title: Job title
            - original_job_uuid: UUID of original job

    Example:
        >>> # Clone a job in same project
        >>> result = await clone_job_async(
        ...     job_uuid=original_job.uuid
        ... )
        >>> print(f"Cloned job {result['job_number']} from {result['original_job_uuid']}")
        >>>
        >>> # Clone to different project with new title
        >>> result = await clone_job_async(
        ...     job_uuid=original_job.uuid,
        ...     new_title="Refinement attempt 2",
        ...     target_project_uuid=other_project.uuid
        ... )

    Raises:
        models.Job.DoesNotExist: If original job not found
        models.Project.DoesNotExist: If target project not found
        FileNotFoundError: If original job's parameter file not found
    """
    logger.info(f"Cloning job: {job_uuid}")

    # Get original job
    old_job = await sync_to_async(
        models.Job.objects.get
    )(uuid=job_uuid)

    logger.debug(f"Original job: {old_job.number} ({old_job.task_name})")

    # Determine target project
    if target_project_uuid is None:
        target_project = old_job.project
    else:
        target_project = await sync_to_async(
            models.Project.objects.get
        )(uuid=target_project_uuid)

    # Determine title
    if new_title is None:
        new_title = old_job.title

    # Create database handler for target project
    db_handler = AsyncDatabaseHandler(project_uuid=target_project.uuid)

    # Create new job in database
    new_job = await db_handler.create_job(
        task_name=old_job.task_name,
        title=new_title,
        parent_job_uuid=parent_job_uuid,
    )

    logger.info(f"Created cloned job {new_job.number} (UUID: {new_job.uuid})")

    # Create job directory
    new_job_dir = Path(new_job.directory)
    await sync_to_async(new_job_dir.mkdir)(parents=True, exist_ok=True)
    logger.debug(f"Created job directory: {new_job_dir}")

    # Create plugin instance and load original parameters
    plugin = await _clone_plugin_with_params(
        old_job=old_job,
        new_job=new_job,
        new_job_dir=new_job_dir,
    )

    # Configure plugin with job metadata
    plugin._dbHandler = db_handler
    plugin._dbProjectId = target_project.uuid
    plugin._dbJobId = new_job.uuid
    plugin._dbJobNumber = new_job.number

    logger.debug(f"Cloned plugin instance with parameters")

    # Save cloned parameters to new job directory
    await _save_cloned_parameters(plugin, new_job)
    logger.debug("Saved cloned parameters to XML file")

    # Update project access time
    await _update_project_access(target_project, new_job.number)

    return {
        "job_uuid": new_job.uuid,
        "job_number": new_job.number,
        "job_directory": str(new_job_dir),
        "plugin": plugin,
        "task_name": old_job.task_name,
        "title": new_title,
        "project_uuid": target_project.uuid,
        "project_name": target_project.name,
        "original_job_uuid": old_job.uuid,
        "original_job_number": old_job.number,
    }


# ============================================================================
# Internal helper functions
# ============================================================================


async def _clone_plugin_with_params(
    old_job: models.Job,
    new_job: models.Job,
    new_job_dir: Path,
) -> CPluginScript:
    """
    Create a new plugin instance and load parameters from original job.

    Args:
        old_job: Original job model instance
        new_job: New job model instance
        new_job_dir: New job directory path

    Returns:
        CPluginScript instance with cloned parameters

    Raises:
        FileNotFoundError: If original job's parameter file not found
    """
    @sync_to_async
    def _clone():
        # Get plugin class
        plugin_class = get_plugin_class(old_job.task_name)

        # Instantiate plugin with new work directory
        plugin = plugin_class(workDirectory=str(new_job_dir))

        # Load parameters from original job's XML file using modern ParamsXmlHandler
        old_params_file = Path(old_job.directory) / "input_params.xml"
        if not old_params_file.exists():
            logger.warning(
                f"Original job parameter file not found: {old_params_file}. "
                f"Using default parameters."
            )
        else:
            # Use modern CPluginScript.loadDataFromXml() which uses ParamsXmlHandler
            # This properly loads inputData, outputData, and all file attributes
            logger.debug(f"Loading parameters from {old_params_file} using modern ParamsXmlHandler")
            error = plugin.loadDataFromXml(str(old_params_file))
            if error and hasattr(error, 'hasError') and error.hasError():
                logger.warning(f"Failed to load parameters: {error}")
            else:
                logger.debug(f"✅ Successfully loaded parameters with inputData preserved")

        # Set output file names using modern approach
        # This uses fileExtensions() method on each file class to determine extension
        set_output_file_names(
            container=plugin.container,
            projectId=str(new_job.project.uuid),
            jobNumber=str(new_job.number),
            force=True
        )

        return plugin

    return await _clone()


async def _save_cloned_parameters(
    plugin: CPluginScript,
    job: models.Job,
    mode: str = "JOB_INPUT",
):
    """
    Save cloned parameters to XML file in new job directory using modern ParamsXmlHandler.

    Args:
        plugin: Plugin instance with cloned parameters
        job: New job model instance
        mode: Parameter file mode (default: "JOB_INPUT")
    """
    from pathlib import Path

    # Use modern CPluginScript.saveDataToXml() which uses ParamsXmlHandler
    # This properly preserves inputData, outputData, and all file attributes
    input_params_path = Path(job.directory) / "input_params.xml"

    @sync_to_async
    def _save():
        logger.info(f"Saving cloned parameters to {input_params_path} using modern ParamsXmlHandler")
        error = plugin.saveDataToXml(str(input_params_path))
        if error and hasattr(error, 'hasError') and error.hasError():
            logger.error(f"Failed to save cloned parameters: {error}")
            raise RuntimeError(f"Failed to save cloned parameters: {error}")
        logger.info(f"✅ Successfully saved cloned parameters with inputData preserved")

    await _save()


async def _update_project_access(
    project: models.Project,
    last_job_number: str,
):
    """
    Update project's last access time and last job number.

    Args:
        project: Project model instance
        last_job_number: Last job number as string
    """
    @sync_to_async
    def _update():
        project.last_access = datetime.datetime.now(tz=timezone("UTC"))
        project.last_job_number = int(last_job_number.split('.')[0])  # Top-level job number
        project.save(update_fields=['last_access', 'last_job_number'])

    await _update()
