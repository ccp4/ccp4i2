"""
Modern async job creation using AsyncDatabaseHandler and CData system.

This module provides typesafe async job creation following patterns from async_run_job.py
and logic from lib/utils/jobs/create.py.

Key features:
- Full async/await with type hints
- Automatic directory structure creation
- Plugin instantiation with parameters
- Parameter file saving (XML)
- Output file path patching
- Database integration via AsyncDatabaseHandler
"""

import asyncio
import datetime
import logging
import uuid
from pathlib import Path
from typing import Optional, Dict, Any, Union

from asgiref.sync import sync_to_async

from core import CCP4TaskManager
from core.CCP4PluginScript import CPluginScript

from ..db import models
from ..db.async_db_handler import AsyncDatabaseHandler
from .utils.containers.remove_defaults import remove_container_default_values
from .utils.parameters.save_params import save_params_for_job
from .utils.files.set_names import set_output_file_names


logger = logging.getLogger(f"ccp4x:{__name__}")


async def create_job_async(
    project_uuid: uuid.UUID,
    task_name: str,
    title: Optional[str] = None,
    parent_job_uuid: Optional[uuid.UUID] = None,
    job_number: Optional[str] = None,
    job_id: Optional[uuid.UUID] = None,
    save_params: bool = True,
    input_params: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Create a new job with full directory structure and parameter files.

    This is the main entry point for async job creation. It:
    1. Creates the job record in the database
    2. Creates the job directory structure
    3. Instantiates the plugin with correct work directory
    4. Patches output file paths to job directory
    5. Optionally saves parameters to XML file
    6. Sets input parameters if provided

    Args:
        project_uuid: UUID of the project
        task_name: Name of the plugin/task (e.g., "ctruncate", "refmac")
        title: Human-readable job title (defaults to plugin title)
        parent_job_uuid: Optional parent job UUID for nested jobs
        job_number: Optional explicit job number (e.g., "1", "1.2")
        job_id: Optional explicit job UUID (auto-generated if not provided)
        save_params: Whether to save parameter file to job directory (default: True)
        input_params: Optional dict of input parameters to set on plugin

    Returns:
        Dict containing:
            - job_uuid: UUID of created job
            - job_number: Job number (e.g., "1", "1.2")
            - job_directory: Path to job directory
            - plugin: CPluginScript instance (configured and ready)
            - task_name: Plugin task name
            - title: Job title

    Example:
        >>> result = await create_job_async(
        ...     project_uuid=project.uuid,
        ...     task_name="ctruncate",
        ...     title="Convert intensities to F",
        ...     input_params={
        ...         "HKLIN": "/path/to/input.mtz",
        ...         "COLANO": "/*/*/[IMEAN,SIGIMEAN]",
        ...     }
        ... )
        >>> print(f"Created job {result['job_number']} at {result['job_directory']}")
        >>> plugin = result['plugin']
        >>> await plugin.execute()

    Raises:
        models.Project.DoesNotExist: If project not found
        models.Job.DoesNotExist: If parent job not found
        Exception: If plugin class not found or instantiation fails
    """
    logger.info(f"Creating job: task={task_name}, project={project_uuid}")

    # Create database handler
    db_handler = AsyncDatabaseHandler(project_uuid=project_uuid)

    # Get project
    project = await db_handler.get_project()
    logger.debug(f"Project: {project.name} ({project.uuid})")

    # Get task manager for plugin info
    task_manager = CCP4TaskManager.CTaskManager()

    # Determine title if not provided
    if title is None:
        title = task_manager.getTitle(task_name)

    # Create job in database
    job = await db_handler.create_job(
        task_name=task_name,
        title=title,
        parent_job_uuid=parent_job_uuid,
        job_number=job_number,
    )

    logger.info(f"Created job {job.number} (UUID: {job.uuid})")

    # Update job UUID if explicitly provided
    if job_id is not None:
        await _update_job_uuid(job, job_id)
        job.uuid = job_id

    # Create job directory structure
    job_dir = Path(job.directory)
    await sync_to_async(job_dir.mkdir)(parents=True, exist_ok=True)
    logger.debug(f"Created job directory: {job_dir}")

    # Create plugin instance
    plugin = await _create_plugin_instance(
        task_name=task_name,
        job_dir=job_dir,
        db_handler=db_handler,
        job=job,
    )

    # Configure plugin with job metadata
    plugin._dbHandler = db_handler
    plugin._dbProjectId = project.uuid
    plugin._dbJobId = job.uuid
    plugin._dbJobNumber = job.number

    logger.debug(f"Created plugin instance: {plugin.__class__.__name__}")

    # Set input parameters if provided
    if input_params:
        await _set_input_parameters(plugin, input_params)

    # Save parameters to file if requested
    if save_params:
        await _save_job_parameters(plugin, job)
        logger.debug("Saved job parameters to XML file")

    return {
        "job_uuid": job.uuid,
        "job_number": job.number,
        "job_directory": str(job_dir),
        "plugin": plugin,
        "task_name": task_name,
        "title": title,
        "project_uuid": project.uuid,
        "project_name": project.name,
    }


async def create_job_with_params_async(
    project_uuid: uuid.UUID,
    task_name: str,
    input_params: Dict[str, Any],
    title: Optional[str] = None,
    parent_job_uuid: Optional[uuid.UUID] = None,
    save_params: bool = True,
) -> Dict[str, Any]:
    """
    Convenience function to create a job and set input parameters in one call.

    This is a simpler interface for common use cases where you just want to
    create a job with specific input parameters.

    Args:
        project_uuid: UUID of the project
        task_name: Name of the plugin/task
        input_params: Dict of input parameter names to values
        title: Optional job title
        parent_job_uuid: Optional parent job UUID
        save_params: Whether to save parameters (default: True)

    Returns:
        Same dict as create_job_async

    Example:
        >>> result = await create_job_with_params_async(
        ...     project_uuid=project.uuid,
        ...     task_name="ctruncate",
        ...     input_params={
        ...         "HKLIN": "/data/input.mtz",
        ...         "COLANO": "/*/*/[IMEAN,SIGIMEAN]",
        ...     }
        ... )
    """
    return await create_job_async(
        project_uuid=project_uuid,
        task_name=task_name,
        title=title,
        parent_job_uuid=parent_job_uuid,
        save_params=save_params,
        input_params=input_params,
    )


async def create_subjob_async(
    parent_job_uuid: uuid.UUID,
    task_name: str,
    title: Optional[str] = None,
    input_params: Optional[Dict[str, Any]] = None,
    save_params: bool = True,
) -> Dict[str, Any]:
    """
    Create a subjob (child job) under an existing parent job.

    This is a convenience wrapper that automatically:
    - Determines project from parent job
    - Sets parent relationship
    - Creates nested job number (e.g., "1.1", "1.2")

    Args:
        parent_job_uuid: UUID of the parent job
        task_name: Plugin task name
        title: Optional job title
        input_params: Optional input parameters
        save_params: Whether to save parameters

    Returns:
        Same dict as create_job_async

    Example:
        >>> # Create parent job
        >>> parent = await create_job_async(project.uuid, "copycell", "Pipeline")
        >>>
        >>> # Create child jobs
        >>> child1 = await create_subjob_async(
        ...     parent['job_uuid'],
        ...     "mtzdump",
        ...     "Extract cell parameters"
        ... )
        >>> child2 = await create_subjob_async(
        ...     parent['job_uuid'],
        ...     "pdbset",
        ...     "Apply cell to PDB"
        ... )
    """
    # Get parent job to determine project
    parent_job = await sync_to_async(
        models.Job.objects.get
    )(uuid=parent_job_uuid)

    return await create_job_async(
        project_uuid=parent_job.project.uuid,
        task_name=task_name,
        title=title,
        parent_job_uuid=parent_job_uuid,
        save_params=save_params,
        input_params=input_params,
    )


# ============================================================================
# Internal helper functions
# ============================================================================


async def _update_job_uuid(job: models.Job, new_uuid: uuid.UUID):
    """Update job UUID in database."""
    @sync_to_async
    def _update():
        job.uuid = new_uuid
        job.save(update_fields=['uuid'])

    await _update()


async def _create_plugin_instance(
    task_name: str,
    job_dir: Path,
    db_handler: AsyncDatabaseHandler,
    job: models.Job,
) -> CPluginScript:
    """
    Create and configure a plugin instance.

    Args:
        task_name: Plugin task name
        job_dir: Job directory path
        db_handler: Database handler instance
        job: Job model instance

    Returns:
        Configured CPluginScript instance
    """
    @sync_to_async
    def _create():
        # Get plugin class
        task_manager = CCP4TaskManager.CTaskManager()
        plugin_class = task_manager.get_plugin_class(task_name)

        # Instantiate plugin with work directory
        plugin = plugin_class(workDirectory=str(job_dir))

        # Remove default values from container
        remove_container_default_values(plugin.container)

        # Set output file names using modern approach
        # This uses fileExtensions() method on each file class to determine extension
        set_output_file_names(
            container=plugin.container,
            projectId=str(job.project.uuid),
            jobNumber=str(job.number),
            force=True
        )

        return plugin

    return await _create()


async def _set_input_parameters(
    plugin: CPluginScript,
    input_params: Dict[str, Any]
):
    """
    Set input parameters on plugin instance.

    Handles both simple attribute assignment and complex objects
    with .set() methods (like CProgramColumnGroup).

    Args:
        plugin: Plugin instance
        input_params: Dict of parameter names to values
    """
    @sync_to_async
    def _set_params():
        for param_name, param_value in input_params.items():
            if not hasattr(plugin.container.inputData, param_name):
                logger.warning(
                    f"Plugin {plugin.__class__.__name__} has no input parameter: {param_name}"
                )
                continue

            param_obj = getattr(plugin.container.inputData, param_name)

            # Handle objects with .set() method (like CProgramColumnGroup)
            if hasattr(param_obj, 'set') and callable(param_obj.set):
                if isinstance(param_value, dict):
                    param_obj.set(param_value)
                else:
                    param_obj.set(param_value)
                logger.debug(f"Set {param_name} via .set() method")

            # Handle file paths
            elif hasattr(param_obj, 'setFullPath'):
                param_obj.setFullPath(str(param_value))
                logger.debug(f"Set {param_name} = {param_value}")

            # Handle direct assignment
            else:
                if hasattr(param_obj, 'value'):
                    param_obj.value = param_value
                else:
                    setattr(plugin.container.inputData, param_name, param_value)
                logger.debug(f"Set {param_name} = {param_value}")

    await _set_params()


async def _save_job_parameters(
    plugin: CPluginScript,
    job: models.Job,
    mode: str = "JOB_INPUT",
    exclude_unset: bool = False,
):
    """
    Save plugin parameters to XML file in job directory.

    Args:
        plugin: Plugin instance
        job: Job model instance
        mode: Parameter file mode (default: "JOB_INPUT")
        exclude_unset: Whether to exclude unset parameters
    """
    @sync_to_async
    def _save():
        save_params_for_job(
            the_job_plugin=plugin,
            the_job=job,
            mode=mode,
            exclude_unset=exclude_unset,
        )

    await _save()


# ============================================================================
# Batch job creation
# ============================================================================


async def create_jobs_batch_async(
    project_uuid: uuid.UUID,
    job_specs: list[Dict[str, Any]],
) -> list[Dict[str, Any]]:
    """
    Create multiple jobs in batch.

    This is more efficient than creating jobs one by one when you need
    to create many jobs at once.

    Args:
        project_uuid: UUID of the project
        job_specs: List of dicts, each containing:
            - task_name (required)
            - title (optional)
            - input_params (optional)
            - parent_job_uuid (optional)
            - save_params (optional, default: True)

    Returns:
        List of job result dicts (same as create_job_async)

    Example:
        >>> jobs = await create_jobs_batch_async(
        ...     project_uuid=project.uuid,
        ...     job_specs=[
        ...         {
        ...             "task_name": "ctruncate",
        ...             "title": "Truncate data",
        ...             "input_params": {"HKLIN": "input.mtz"},
        ...         },
        ...         {
        ...             "task_name": "refmac",
        ...             "title": "Refinement",
        ...             "input_params": {"HKLIN": "truncated.mtz"},
        ...         },
        ...     ]
        ... )
    """
    results = []

    for spec in job_specs:
        result = await create_job_async(
            project_uuid=project_uuid,
            task_name=spec['task_name'],
            title=spec.get('title'),
            parent_job_uuid=spec.get('parent_job_uuid'),
            save_params=spec.get('save_params', True),
            input_params=spec.get('input_params'),
        )
        results.append(result)

    return results


# ============================================================================
# Example usage
# ============================================================================


async def example_usage():
    """Example demonstrating typical usage patterns."""

    # Example 1: Simple job creation
    result = await create_job_async(
        project_uuid=uuid.UUID("..."),
        task_name="ctruncate",
        title="Convert intensities to amplitudes",
    )
    print(f"Created job {result['job_number']} at {result['job_directory']}")

    # Example 2: Job with input parameters
    result = await create_job_with_params_async(
        project_uuid=uuid.UUID("..."),
        task_name="ctruncate",
        input_params={
            "HKLIN": "/path/to/input.mtz",
            "COLANO": "/*/*/[IMEAN,SIGIMEAN]",
        }
    )

    # Example 3: Create parent and child jobs
    parent = await create_job_async(
        project_uuid=uuid.UUID("..."),
        task_name="copycell",
        title="Data processing pipeline",
    )

    child1 = await create_subjob_async(
        parent_job_uuid=parent['job_uuid'],
        task_name="mtzdump",
        title="Extract cell parameters",
    )

    child2 = await create_subjob_async(
        parent_job_uuid=parent['job_uuid'],
        task_name="pdbset",
        title="Apply cell to PDB",
    )

    # Example 4: Batch job creation
    jobs = await create_jobs_batch_async(
        project_uuid=uuid.UUID("..."),
        job_specs=[
            {"task_name": "ctruncate", "title": "Step 1"},
            {"task_name": "refmac", "title": "Step 2"},
            {"task_name": "coot", "title": "Step 3"},
        ]
    )


if __name__ == "__main__":
    # Run examples
    asyncio.run(example_usage())
