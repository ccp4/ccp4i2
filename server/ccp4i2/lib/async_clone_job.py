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

from ccp4i2.core import CCP4TaskManager
from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core.CCP4Container import CContainer

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


async def clone_job_to_project_async(
    job_uuid: uuid.UUID,
    target_project_uuid: uuid.UUID,
    new_title: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Clone a job from one project to another project.

    This is a convenience wrapper for cross-project job cloning.

    Args:
        job_uuid: UUID of the job to clone
        target_project_uuid: UUID of the target project
        new_title: Optional new title for the cloned job

    Returns:
        Same dict as clone_job_async

    Example:
        >>> result = await clone_job_to_project_async(
        ...     job_uuid=job_from_project_a.uuid,
        ...     target_project_uuid=project_b.uuid,
        ...     new_title="Ported from Project A"
        ... )
    """
    return await clone_job_async(
        job_uuid=job_uuid,
        new_title=new_title,
        target_project_uuid=target_project_uuid,
    )


async def clone_job_as_subjob_async(
    job_uuid: uuid.UUID,
    parent_job_uuid: uuid.UUID,
    new_title: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Clone a job as a subjob of another job.

    This creates a child job with cloned parameters.

    Args:
        job_uuid: UUID of the job to clone
        parent_job_uuid: UUID of the parent job
        new_title: Optional new title

    Returns:
        Same dict as clone_job_async

    Example:
        >>> # Create a pipeline parent
        >>> parent = await create_job_async(project.uuid, "copycell", "Pipeline")
        >>>
        >>> # Clone existing job as subjob
        >>> child = await clone_job_as_subjob_async(
        ...     job_uuid=existing_job.uuid,
        ...     parent_job_uuid=parent['job_uuid'],
        ...     new_title="Step 1 (cloned)"
        ... )
    """
    return await clone_job_async(
        job_uuid=job_uuid,
        new_title=new_title,
        parent_job_uuid=parent_job_uuid,
    )


async def clone_jobs_batch_async(
    job_uuids: list[uuid.UUID],
    target_project_uuid: Optional[uuid.UUID] = None,
    new_titles: Optional[list[str]] = None,
) -> list[Dict[str, Any]]:
    """
    Clone multiple jobs in batch.

    Args:
        job_uuids: List of job UUIDs to clone
        target_project_uuid: Optional target project (defaults to original projects)
        new_titles: Optional list of new titles (same length as job_uuids)

    Returns:
        List of job result dicts (same as clone_job_async)

    Example:
        >>> jobs = await clone_jobs_batch_async(
        ...     job_uuids=[job1.uuid, job2.uuid, job3.uuid],
        ...     target_project_uuid=new_project.uuid
        ... )
    """
    results = []

    for i, job_uuid in enumerate(job_uuids):
        new_title = None
        if new_titles and i < len(new_titles):
            new_title = new_titles[i]

        result = await clone_job_async(
            job_uuid=job_uuid,
            new_title=new_title,
            target_project_uuid=target_project_uuid,
        )
        results.append(result)
        logger.info(f"Cloned {i+1}/{len(job_uuids)}: {result['job_number']}")

    return results


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
        task_manager = CCP4TaskManager.CTaskManager()
        plugin_class = task_manager.get_plugin_class(old_job.task_name)

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


# ============================================================================
# Advanced cloning patterns
# ============================================================================


async def clone_job_with_modifications_async(
    job_uuid: uuid.UUID,
    param_modifications: Dict[str, Any],
    new_title: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Clone a job and modify specific parameters before saving.

    This is useful when you want to clone a job but change some inputs.

    Args:
        job_uuid: UUID of the job to clone
        param_modifications: Dict of parameter names to new values
        new_title: Optional new title

    Returns:
        Same dict as clone_job_async

    Example:
        >>> # Clone a job but use different input file
        >>> result = await clone_job_with_modifications_async(
        ...     job_uuid=original.uuid,
        ...     param_modifications={
        ...         "HKLIN": "/path/to/different/input.mtz",
        ...         "COLANO": "/*/*/[I,SIGI]",
        ...     },
        ...     new_title="Re-run with different data"
        ... )
    """
    # Clone the job
    result = await clone_job_async(
        job_uuid=job_uuid,
        new_title=new_title,
    )

    # Get the plugin from result
    plugin = result['plugin']

    # Apply modifications
    @sync_to_async
    def _modify_params():
        for param_name, param_value in param_modifications.items():
            if not hasattr(plugin.container.inputData, param_name):
                logger.warning(
                    f"Plugin has no input parameter: {param_name}, skipping"
                )
                continue

            param_obj = getattr(plugin.container.inputData, param_name)

            # Handle objects with .set() method
            if hasattr(param_obj, 'set') and callable(param_obj.set):
                if isinstance(param_value, dict):
                    param_obj.set(param_value)
                else:
                    param_obj.set(param_value)

            # Handle file paths
            elif hasattr(param_obj, 'setFullPath'):
                param_obj.setFullPath(str(param_value))

            # Handle direct assignment
            else:
                if hasattr(param_obj, 'value'):
                    param_obj.value = param_value
                else:
                    setattr(plugin.container.inputData, param_name, param_value)

            logger.debug(f"Modified parameter {param_name} = {param_value}")

    await _modify_params()

    # Re-save parameters with modifications
    new_job = await sync_to_async(
        models.Job.objects.get
    )(uuid=result['job_uuid'])

    await _save_cloned_parameters(plugin, new_job)
    logger.info("Saved modified parameters")

    return result


async def clone_job_chain_async(
    job_uuids: list[uuid.UUID],
    target_project_uuid: uuid.UUID,
    parent_job_uuid: Optional[uuid.UUID] = None,
) -> list[Dict[str, Any]]:
    """
    Clone a chain of jobs as subjobs of a parent.

    This creates a pipeline by cloning multiple jobs as children of a parent job.

    Args:
        job_uuids: List of job UUIDs to clone in order
        target_project_uuid: Target project UUID
        parent_job_uuid: Optional existing parent job UUID (creates new if not provided)

    Returns:
        List of cloned job result dicts

    Example:
        >>> # Clone a series of jobs as a pipeline
        >>> pipeline = await clone_job_chain_async(
        ...     job_uuids=[ctruncate_job.uuid, refmac_job.uuid, coot_job.uuid],
        ...     target_project_uuid=project.uuid
        ... )
    """
    # Create parent job if not provided
    if parent_job_uuid is None:
        from .async_create_job import create_job_async
        parent_result = await create_job_async(
            project_uuid=target_project_uuid,
            task_name="pipeline",
            title="Cloned pipeline",
        )
        parent_job_uuid = parent_result['job_uuid']
        logger.info(f"Created parent pipeline job: {parent_result['job_number']}")

    # Clone each job as a subjob
    results = []
    for i, job_uuid in enumerate(job_uuids):
        result = await clone_job_as_subjob_async(
            job_uuid=job_uuid,
            parent_job_uuid=parent_job_uuid,
            new_title=f"Step {i+1} (cloned)",
        )
        results.append(result)
        logger.info(f"Cloned step {i+1}/{len(job_uuids)}: {result['job_number']}")

    return results


# ============================================================================
# Example usage
# ============================================================================


async def example_usage():
    """Example demonstrating typical cloning patterns."""

    # Example 1: Simple job cloning in same project
    result = await clone_job_async(
        job_uuid=uuid.UUID("...")
    )
    print(f"Cloned job {result['original_job_number']} → {result['job_number']}")

    # Example 2: Clone to different project with new title
    result = await clone_job_to_project_async(
        job_uuid=uuid.UUID("..."),
        target_project_uuid=uuid.UUID("..."),
        new_title="Ported refinement job"
    )

    # Example 3: Clone as subjob
    parent_uuid = uuid.UUID("...")
    child = await clone_job_as_subjob_async(
        job_uuid=uuid.UUID("..."),
        parent_job_uuid=parent_uuid,
        new_title="Cloned as step 1"
    )

    # Example 4: Clone with modifications
    result = await clone_job_with_modifications_async(
        job_uuid=uuid.UUID("..."),
        param_modifications={
            "HKLIN": "/path/to/new/input.mtz",
            "COLANO": "/*/*/[I,SIGI]",
        },
        new_title="Re-run with different data"
    )

    # Example 5: Batch clone
    results = await clone_jobs_batch_async(
        job_uuids=[
            uuid.UUID("..."),
            uuid.UUID("..."),
            uuid.UUID("..."),
        ],
        target_project_uuid=uuid.UUID("...")
    )

    # Example 6: Clone job chain as pipeline
    pipeline = await clone_job_chain_async(
        job_uuids=[
            uuid.UUID("..."),  # ctruncate
            uuid.UUID("..."),  # refmac
            uuid.UUID("..."),  # coot
        ],
        target_project_uuid=uuid.UUID("...")
    )


if __name__ == "__main__":
    # Run examples
    asyncio.run(example_usage())
