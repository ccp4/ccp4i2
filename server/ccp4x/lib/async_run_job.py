"""
Modern async job runner using new CData system and AsyncDatabaseHandler.

This module replaces the legacy run_job.py with:
- Async/await throughout (no Qt event loop)
- Modern CData introspection
- Automatic database tracking via signals
- Clean context manager patterns
"""

import asyncio
import datetime
import logging
import uuid
import xml.etree.ElementTree as ET
from contextlib import asynccontextmanager
from pathlib import Path
from typing import Optional

from asgiref.sync import sync_to_async
from django.utils import timezone

logger = logging.getLogger(f"ccp4x:{__name__}")


async def run_job_async(job_uuid: uuid.UUID, project_uuid: Optional[uuid.UUID] = None):
    """
    Modern async job execution with automatic database tracking.

    This function:
    1. Retrieves job from database
    2. Creates plugin instance
    3. Imports input files (async)
    4. Executes plugin (async)
    5. Gleans output files and KPIs (automatic via track_job)
    6. Updates job status (automatic via track_job)

    Args:
        job_uuid: UUID of the job to run
        project_uuid: Optional project UUID (will be fetched from job if not provided)

    Returns:
        Plugin execution result

    Example:
        >>> result = await run_job_async(job_uuid)
        >>> print(f"Job completed with status: {result}")
    """
    from ..db import models
    from ..db.async_db_handler import AsyncDatabaseHandler
    from .async_import_files import import_input_files_async
    from .utils.plugins.get_plugin import get_job_plugin

    # Get job from database with related project
    job = await sync_to_async(models.Job.objects.select_related('project').get)(uuid=job_uuid)
    logger.info(f"Running job {job.number} ({job.task_name}) in {job.directory}")

    # Determine project UUID (project is prefetched via select_related)
    if project_uuid is None:
        project_uuid = job.project.uuid

    # Create database handler
    db_handler = AsyncDatabaseHandler(project_uuid=project_uuid)

    # Create or retrieve plugin instance
    plugin = await create_plugin_for_job(job, db_handler)

    try:
        # Update job status to RUNNING
        await db_handler.update_job_status(job.uuid, models.Job.Status.RUNNING)

        # Import input files (async)
        logger.info("Importing input files...")
        files_imported = await import_input_files_async(job, plugin, db_handler)
        logger.info(f"Imported {files_imported} input files")

        # Set output file names (guaranteed step, even if plugin overrides process())
        # This ensures output files have project/relPath/baseName set
        from .utils.files.set_names import set_output_file_names

        logger.info("Setting output file names...")
        await sync_to_async(set_output_file_names)(
            container=plugin.container,
            projectId=str(job.project.uuid),
            jobNumber=job.number,
            force=True
        )
        logger.info("Output file names set")

        # Save params.xml with all file attributes before execution
        # This ensures nested jobs can load parent parameters
        logger.info("Saving params.xml before execution...")

        # DEBUG: Check what's in inputData before saving
        if hasattr(plugin.container, 'inputData'):
            input_files = []
            for attr_name in dir(plugin.container.inputData):
                if not attr_name.startswith('_'):
                    try:
                        attr = getattr(plugin.container.inputData, attr_name)
                        if hasattr(attr, 'baseName'):
                            basename = str(attr.baseName) if hasattr(attr, 'baseName') else ''
                            input_files.append(f"{attr_name}={basename}")
                    except:
                        pass
            logger.debug(f"[DEBUG async_run_job] inputData before save: {input_files}")
            logger.debug(f"[DEBUG] inputData before save: {input_files}")

        await sync_to_async(plugin.saveDataToXml)(str(job.directory / "params.xml"))
        logger.info("Saved params.xml")

        # Execute plugin with automatic tracking
        logger.info("Executing plugin...")
        async with db_handler.track_job(plugin):
            # This automatically:
            # - Updates status via signals
            # - Gleans files on completion
            # - Gleans KPIs on completion

            # Call process() in a thread since it's not async
            # process() handles: checkInputData, makeCommandAndScript, startProcess
            # (checkOutputData may or may not be called if plugin overrides process)
            result = await sync_to_async(plugin.process)()

        # Write diagnostic.xml with error report (always, for debugging)
        await write_diagnostic_xml(plugin, job.directory)

        logger.info(f"Job {job.number} completed successfully")
        return result

    except Exception as e:
        logger.exception(f"Job {job.number} failed: {e}")

        # Write diagnostic.xml before updating status (critical for debugging failures)
        await write_diagnostic_xml(plugin, job.directory)

        # Update status to FAILED
        await db_handler.update_job_status(job.uuid, models.Job.Status.FAILED)

        raise


async def write_diagnostic_xml(plugin, job_directory):
    """
    Write diagnostic.xml with error report from the plugin.

    This file contains debugging information collected during job execution,
    including any errors, warnings, and status messages from the plugin.

    Args:
        plugin: CPluginScript instance with errorReport attribute
        job_directory: Path to job directory where diagnostic.xml will be written
    """
    try:
        diagnostic_path = Path(job_directory) / "diagnostic.xml"
        error_report = plugin.errorReport.getEtree()
        ET.indent(error_report, space="\t", level=0)
        with open(diagnostic_path, "wb") as f:
            f.write(ET.tostring(error_report, encoding="utf-8"))
        logger.info(f"Wrote diagnostic.xml to {diagnostic_path}")
    except Exception as err:
        logger.warning(f"Failed to write diagnostic.xml: {err}")


async def create_plugin_for_job(job, db_handler):
    """
    Create a plugin instance for a job.

    Args:
        job: Django Job model instance
        db_handler: AsyncDatabaseHandler instance

    Returns:
        CPluginScript instance
    """
    from core import CCP4TaskManager
    from .utils.plugins.get_plugin import get_job_plugin

    # Get plugin class
    task_manager = CCP4TaskManager.CTaskManager()
    plugin_class = task_manager.get_plugin_class(job.task_name)

    # Create plugin instance
    plugin = plugin_class(
        workDirectory=str(job.directory),
        name=job.title,
    )

    # Force synchronous execution for i2run top-level jobs
    # This overrides ASYNCHRONOUS=True class variable to ensure subprocess completion
    # before returning control to tests. Pipelines can still use async for sub-jobs.
    plugin.doAsync = False
    logger.info(f"Set plugin.doAsync=False to force synchronous execution (was ASYNCHRONOUS={plugin.ASYNCHRONOUS})")

    # Set database context using the proper API
    plugin.setDbData(
        handler=db_handler,
        projectId=str(job.project.uuid),
        jobNumber=job.number,
        jobId=str(job.uuid)
    )
    logger.info(f"Set plugin database context: jobId={plugin.get_db_job_id()}, projectId={plugin._dbProjectId}, jobNumber={plugin._dbJobNumber}")

    # Load parameters if they exist
    params_file = job.directory / "input_params.xml"
    if params_file.exists():
        await load_plugin_params(plugin, params_file)
        logger.info(f"After loading params: _dbJobId={plugin._dbJobId}")

    return plugin


async def load_plugin_params(plugin, params_file: Path):
    """
    Load plugin parameters from XML file using modern ParamsXmlHandler.

    Args:
        plugin: CPluginScript instance
        params_file: Path to parameters XML file
    """
    try:
        # Use modern CPluginScript.loadDataFromXml() which uses ParamsXmlHandler
        # This properly loads inputData, outputData, and all file attributes
        logger.info(f"Loading parameters from {params_file} using modern ParamsXmlHandler")
        error = await sync_to_async(plugin.loadDataFromXml)(str(params_file))
        if error and hasattr(error, 'hasError') and error.hasError():
            logger.error(f"Failed to load parameters: {error}")
        else:
            logger.info(f"✅ Successfully loaded parameters with inputData preserved")
    except Exception as e:
        logger.warning(f"Could not load parameters: {e}")
        import traceback
        traceback.print_exc()


@asynccontextmanager
async def job_execution_context(job):
    """
    Context manager for job execution with proper cleanup.

    Handles:
    - Output redirection to job directory
    - Process ID tracking
    - Timing
    - Final status enforcement

    Args:
        job: Django Job model instance

    Example:
        >>> async with job_execution_context(job):
        ...     await execute_plugin(plugin)
    """
    from ..db import models

    # Record start time
    start_time = timezone.now()

    # Ensure job directory exists
    job_dir = Path(job.directory)
    await sync_to_async(job_dir.mkdir)(parents=True, exist_ok=True)

    try:
        yield job

    finally:
        # Reload job to get latest status
        @sync_to_async
        def reload_job():
            return models.Job.objects.get(id=job.id)

        updated_job = await reload_job()

        # Ensure job has terminal status
        @sync_to_async
        def finalize_status():
            if updated_job.status not in [
                models.Job.Status.FINISHED,
                models.Job.Status.FAILED,
                models.Job.Status.UNSATISFACTORY,
                models.Job.Status.INTERRUPTED,
            ]:
                updated_job.status = models.Job.Status.FAILED
                updated_job.save()

            # Clear process ID
            updated_job.process_id = None
            updated_job.save()

        await finalize_status()

        # Log execution time
        elapsed = timezone.now() - start_time
        logger.info(f"Job {job.number} execution time: {elapsed}")


async def run_pipeline_async(
    project_uuid: uuid.UUID,
    task_name: str,
    title: str,
    input_data: dict,
    parent_job_uuid: Optional[uuid.UUID] = None,
) -> uuid.UUID:
    """
    High-level function to create and run a job from scratch.

    This is useful for creating jobs programmatically without going through
    the Django management command.

    Args:
        project_uuid: UUID of the project
        task_name: Name of the plugin (e.g., "ctruncate")
        title: Human-readable job title
        input_data: Dictionary of input parameters
        parent_job_uuid: Optional parent job UUID for nested execution

    Returns:
        UUID of the created job

    Example:
        >>> job_uuid = await run_pipeline_async(
        ...     project_uuid=project.uuid,
        ...     task_name="ctruncate",
        ...     title="Convert intensities to amplitudes",
        ...     input_data={
        ...         "HKLIN": "/path/to/input.mtz",
        ...         "COLANO": "/*/*/[IMEAN,SIGIMEAN]",
        ...     }
        ... )
    """
    from ..db.async_db_handler import AsyncDatabaseHandler

    # Create database handler
    db_handler = AsyncDatabaseHandler(project_uuid=project_uuid)

    # Create job in database
    job = await db_handler.create_job(
        task_name=task_name,
        title=title,
        parent_job_uuid=parent_job_uuid,
    )

    logger.info(f"Created job {job.number} ({task_name})")

    # Create plugin instance
    plugin = await create_plugin_for_job(job, db_handler)

    # Set input parameters
    for param_name, param_value in input_data.items():
        if hasattr(plugin.inputData, param_name):
            param_obj = getattr(plugin.inputData, param_name)
            if hasattr(param_obj, 'set'):
                param_obj.set(param_value)
            else:
                setattr(plugin.inputData, param_name, param_value)
            logger.info(f"Set {param_name} = {param_value}")
        else:
            logger.warning(f"Plugin has no input parameter: {param_name}")

    # Run the job
    await run_job_async(job.uuid, project_uuid)

    return job.uuid


async def run_nested_jobs_async(
    project_uuid: uuid.UUID,
    parent_task_name: str,
    parent_title: str,
    child_specs: list,
) -> uuid.UUID:
    """
    Run a pipeline with nested child jobs.

    Args:
        project_uuid: UUID of the project
        parent_task_name: Name of parent plugin
        parent_title: Title of parent job
        child_specs: List of dicts with 'task_name', 'title', 'input_data'

    Returns:
        UUID of parent job

    Example:
        >>> parent_uuid = await run_nested_jobs_async(
        ...     project_uuid=project.uuid,
        ...     parent_task_name="copycell",
        ...     parent_title="Data processing pipeline",
        ...     child_specs=[
        ...         {
        ...             "task_name": "ctruncate",
        ...             "title": "Truncate",
        ...             "input_data": {"HKLIN": "input.mtz"},
        ...         },
        ...         {
        ...             "task_name": "refmac",
        ...             "title": "Refine",
        ...             "input_data": {"HKLIN": "truncated.mtz"},
        ...         },
        ...     ]
        ... )
    """
    from ..db.async_db_handler import AsyncDatabaseHandler

    # Create parent job
    db_handler = AsyncDatabaseHandler(project_uuid=project_uuid)
    parent_job = await db_handler.create_job(
        task_name=parent_task_name,
        title=parent_title,
    )

    logger.info(f"Created parent job {parent_job.number}")

    # Run child jobs
    for child_spec in child_specs:
        child_uuid = await run_pipeline_async(
            project_uuid=project_uuid,
            task_name=child_spec['task_name'],
            title=child_spec['title'],
            input_data=child_spec['input_data'],
            parent_job_uuid=parent_job.uuid,
        )
        logger.info(f"Completed child job {child_uuid}")

    # Mark parent as finished
    await db_handler.update_job_status(
        parent_job.uuid,
        finish_time=timezone.now(),
    )

    return parent_job.uuid


# Example usage
async def example_usage():
    """Example of how to use the modern async job runner."""

    # Example 1: Run a single job by UUID
    job_uuid = uuid.UUID("...")  # Get from database
    result = await run_job_async(job_uuid)

    # Example 2: Create and run a new job
    project_uuid = uuid.UUID("...")
    job_uuid = await run_pipeline_async(
        project_uuid=project_uuid,
        task_name="ctruncate",
        title="Convert to amplitudes",
        input_data={
            "HKLIN": "/data/input.mtz",
            "COLANO": "/*/*/[IMEAN,SIGIMEAN]",
        }
    )

    # Example 3: Run nested jobs
    parent_uuid = await run_nested_jobs_async(
        project_uuid=project_uuid,
        parent_task_name="copycell",
        parent_title="Processing pipeline",
        child_specs=[
            {
                "task_name": "ctruncate",
                "title": "Step 1: Truncate",
                "input_data": {"HKLIN": "/data/input.mtz"},
            },
            {
                "task_name": "refmac",
                "title": "Step 2: Refine",
                "input_data": {"HKLIN": "/data/truncated.mtz"},
            },
        ]
    )


if __name__ == "__main__":
    # Run example
    asyncio.run(example_usage())
