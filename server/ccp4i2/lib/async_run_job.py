"""
Modern async job runner using new CData system and AsyncDatabaseHandler.

This module replaces the legacy run_job.py with:
- Async/await throughout (no Qt event loop)
- Modern CData introspection
- Automatic database tracking via signals
- Clean context manager patterns
"""

import logging
import uuid
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Optional

from asgiref.sync import sync_to_async

from ccp4i2.core.tasks import get_plugin_class

logger = logging.getLogger(f"ccp4i2:{__name__}")


class _ProcessFailedError(Exception):
    """Raised when process() returns FAILED (not an unexpected exception)."""
    pass


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

        # Run validity() to allow plugins to adjust qualifiers
        # This is critical for pipelines like servalcat_pipe which set allowUndefined=True
        # on embedded wrappers (e.g., metalCoordWrapper.inputData.XYZIN)
        logger.info("Running validity check...")
        validity_error = await sync_to_async(plugin.validity)()
        if validity_error and hasattr(validity_error, 'maxSeverity'):
            from ccp4i2.core.base_object.error_reporting import Severity
            if validity_error.maxSeverity() >= Severity.WARNING:
                logger.warning(f"Validity check has warnings/errors (continuing): {validity_error.report()}")

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
            try:
                result = await sync_to_async(plugin.process)()

                # Check if process() returned FAILED
                # This catches failures in processInputFiles, makeCommandAndScript, etc.
                # that don't throw exceptions but return FAILED status
                if result == plugin.FAILED:
                    logger.error(f"Job {job.number} - process() returned FAILED")
                    # Error should already be in errorReport, so write diagnostic now
                    # and re-raise directly (don't fall through to the except block
                    # which would overwrite the real errors with a generic code 999)
                    await write_diagnostic_xml(plugin, job.directory)
                    await db_handler.update_job_status(job.uuid, models.Job.Status.FAILED)
                    raise _ProcessFailedError(f"Job failed during process() - see diagnostic.xml for details")

            except _ProcessFailedError:
                # process() returned FAILED — real errors are already in errorReport.
                # Re-raise without adding a generic code 999 entry.
                raise

            except Exception as proc_exc:
                # Capture unexpected Python exceptions into the error report
                # This ensures they appear in diagnostic.xml, not just cplusplus_stdout.txt
                import traceback
                tb_str = ''.join(traceback.format_exception(type(proc_exc), proc_exc, proc_exc.__traceback__))

                # Add to plugin's error report
                plugin.errorReport.append(
                    klass=plugin.__class__.__name__,
                    code=999,  # Generic Python exception code
                    details=f'Python exception during process(): {type(proc_exc).__name__}: {str(proc_exc)}\n\n{tb_str}',
                    name='process',
                    severity=4  # ERROR level
                )

                logger.exception(f"Job {job.number} failed with Python exception: {proc_exc}")

                # Re-raise to trigger outer exception handler
                raise

        # For async tasks, call postProcessCheck to examine exit codes
        # This must be done AFTER process() completes and BEFORE writing diagnostic.xml
        # so that exit code errors are captured in the error report
        logger.info(f"Calling postProcessCheck for job {job.number}")

        status = plugin.SUCCEEDED
        try:
            # Get processId for legacy plugins that override postProcessCheck(processId)
            # Pass it for backward compatibility with legacy plugin signatures
            process_id = getattr(plugin, '_runningProcessId', None)

            # postProcessCheck always returns tuple (status, exitStatus, exitCode)
            status, exit_status, exit_code = await sync_to_async(plugin.postProcessCheck)(process_id)

            logger.info(f"postProcessCheck returned status: {status}, exitStatus: {exit_status}, exitCode: {exit_code} "
                       f"(SUCCEEDED={plugin.SUCCEEDED}, FAILED={plugin.FAILED})")

        except Exception as post_exc:
            # If postProcessCheck itself fails, add to error report
            logger.exception(f"postProcessCheck failed for job {job.number}: {post_exc}")
            plugin.errorReport.append(
                klass=plugin.__class__.__name__,
                code=994,
                details=f"postProcessCheck failed: {str(post_exc)}",
                name="postProcessCheck",
                severity=4  # ERROR
            )
            status = plugin.FAILED

        # Write diagnostic.xml with error report (always, for debugging)
        # This includes any errors found by postProcessCheck or errors during postProcessCheck itself
        await write_diagnostic_xml(plugin, job.directory)

        # Check if postProcessCheck found errors
        if status == plugin.FAILED:
            logger.error(f"Job {job.number} failed - postProcessCheck returned FAILED")
            await db_handler.update_job_status(job.uuid, models.Job.Status.FAILED)
            raise Exception(f"Job failed - see diagnostic.xml for details")

        # Explicitly update job status to FINISHED
        # This is belt-and-braces: track_job should have done this, but some legacy
        # pipelines don't set plugin._status properly, causing track_job to skip the update
        await db_handler.update_job_status(job.uuid, models.Job.Status.FINISHED)

        logger.info(f"Job {job.number} completed successfully")
        return result

    except Exception as e:
        logger.exception(f"Job {job.number} failed: {e}")

        # Write diagnostic.xml before updating status (critical for debugging failures)
        # This will now include the exception details added above
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

    # Get plugin class
    plugin_class = get_plugin_class(job.task_name)

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
