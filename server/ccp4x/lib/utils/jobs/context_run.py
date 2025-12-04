"""
Context-Dependent Job Execution Module

Provides environment-aware job execution that adapts to deployment context:
- Local Mode: Executes jobs via subprocess (laptop/development)
- Azure Mode: Queues jobs via Azure Service Bus (container apps)

The execution mode is determined automatically from environment variables,
keeping Azure-specific dependencies isolated and only loading when needed.

Environment Variables:
    EXECUTION_MODE: Explicit mode ('local' or 'azure')
    SERVICE_BUS_CONNECTION_STRING: Azure connection (implies azure mode)
    SERVICE_BUS_QUEUE_NAME: Azure queue name (default: 'job-queue')
    CCP4: Path to CCP4 installation (required for local mode)

Example Usage:
    from ccp4x.lib.context_dependent_run import run_job_context_aware

    result = run_job_context_aware(job)
    if result["success"]:
        # Job started/queued successfully
        return Response(result["data"])
    else:
        # Handle error
        return Response({"error": result["error"]}, status=result["status"])
"""

import os
import json
import logging
import subprocess
import pathlib
import platform

logger = logging.getLogger(__name__)


def _lazy_import_azure_servicebus():
    """
    Lazy import Azure Service Bus dependencies.

    Only imports Azure libraries when running in Azure mode, keeping
    local environments free of Azure dependencies.

    Returns:
        tuple: (ServiceBusClient, ServiceBusMessage) classes

    Raises:
        ImportError: If Azure libraries not installed when needed
    """
    try:
        from azure.servicebus import ServiceBusClient, ServiceBusMessage

        return ServiceBusClient, ServiceBusMessage
    except ImportError:
        raise ImportError(
            "Azure Service Bus libraries not installed. "
            "Required for Azure execution mode. "
            "Install with: pip install azure-servicebus azure-identity"
        )


def get_execution_mode():
    """
    Determine execution mode from environment variables.

    Detection Priority:
    1. EXECUTION_MODE env var (explicit: 'local' or 'azure')
    2. Presence of SERVICE_BUS_CONNECTION_STRING (implicit azure)
    3. Default to 'local'

    Returns:
        str: 'local' or 'azure'

    Example:
        >>> os.environ["EXECUTION_MODE"] = "azure"
        >>> get_execution_mode()
        'azure'

        >>> os.environ["SERVICE_BUS_CONNECTION_STRING"] = "Endpoint=..."
        >>> get_execution_mode()
        'azure'
    """
    # Explicit mode setting takes precedence
    explicit_mode = os.getenv("EXECUTION_MODE", "").lower()
    if explicit_mode in ["local", "azure"]:
        logger.info("Using explicit execution mode: %s", explicit_mode)
        return explicit_mode

    # Implicit detection based on Azure configuration
    if os.getenv("SERVICE_BUS_CONNECTION_STRING"):
        logger.info("Detected Azure Service Bus config, using azure mode")
        return "azure"

    # Default to local mode
    logger.info("No Azure config detected, defaulting to local mode")
    return "local"


def run_job_azure(job):
    """
    Execute job via Azure Service Bus queue.

    Sends a message to the Azure Service Bus queue containing job details,
    allowing asynchronous processing by worker container apps.

    Args:
        job: Job model instance with attributes:
            - id: Job primary key
            - uuid: Job UUID
            - task_name: Name of the task to execute
            - project: Related project with uuid attribute

    Returns:
        dict: Result dictionary with keys:
            - success (bool): True if job queued successfully
            - data (dict): Serialized job data (if success)
            - error (str): Error message (if failure)
            - status (int): HTTP status code

    Message Format:
        {
            "action": "run_job",
            "job_uuid": "550e8400-e29b-41d4-a716-446655440000",
            "job_id": 123,
            "task_name": "refmac5",
            "project_uuid": "project-uuid-here"
        }

    Raises:
        No exceptions - all errors returned in result dict
    """
    logger.info("Running job %s in AZURE mode via Service Bus", job.id)

    try:
        # Lazy import Azure dependencies
        ServiceBusClient, ServiceBusMessage = _lazy_import_azure_servicebus()

        # Prepare message payload
        message_body = {
            "action": "run_job",
            "job_uuid": str(job.uuid),
            "job_id": job.id,
            "task_name": job.task_name,
            "project_uuid": str(job.project.uuid),
        }

        # Get Service Bus configuration
        connection_string = os.getenv("SERVICE_BUS_CONNECTION_STRING")
        queue_name = os.getenv("SERVICE_BUS_QUEUE_NAME", "job-queue")

        if not connection_string:
            error_msg = "Azure Service Bus connection string not configured"
            logger.error(error_msg)
            return {
                "success": False,
                "error": "Service Bus configuration missing",
                "status": 500,
            }

        # Send message to Service Bus
        with ServiceBusClient.from_connection_string(connection_string) as client:
            with client.get_queue_sender(queue_name) as sender:
                message = ServiceBusMessage(json.dumps(message_body))
                sender.send_messages(message)

        # Update job status to QUEUED
        # Import here to avoid circular imports
        from ccp4x.db import models

        job.status = models.Job.Status.QUEUED
        job.save()

        logger.info("Queued job %s (%s) via Azure Service Bus", job.id, job.uuid)

        return {
            "success": True,
            "data": job,
            "status": 200,
        }

    except ImportError as import_error:
        logger.exception("Azure libraries not available", exc_info=import_error)
        return {
            "success": False,
            "error": str(import_error),
            "status": 500,
        }
    except Exception as error:
        logger.exception("Failed to queue job via Service Bus", exc_info=error)
        return {
            "success": False,
            "error": f"Service Bus error: {str(error)}",
            "status": 500,
        }


def _find_python_interpreter(project_root: pathlib.Path) -> tuple:
    """
    Find the appropriate Python interpreter for job execution.

    Priority:
    1. ccp4-python (if available on PATH after sourcing ccp4.setup-sh)
    2. .venv/bin/python (development virtual environment)
    3. .venv.py311/bin/python (alternative naming convention)

    Args:
        project_root: Path to the project root directory

    Returns:
        tuple: (interpreter_path: str, interpreter_name: str) or (None, None) if not found
    """
    import shutil

    # Priority 1: ccp4-python on PATH (preferred for CCP4 environment)
    ccp4_python = shutil.which("ccp4-python")
    if ccp4_python:
        logger.info("Found ccp4-python on PATH: %s", ccp4_python)
        return ccp4_python, "ccp4-python"

    # Priority 2: Project virtual environment
    venv_python = project_root / ".venv" / "bin" / "python"
    if venv_python.exists():
        logger.info("Found .venv Python: %s", venv_python)
        return str(venv_python), ".venv/bin/python"

    # Priority 3: Alternative venv naming
    venv_python = project_root / ".venv.py311" / "bin" / "python"
    if venv_python.exists():
        logger.info("Found .venv.py311 Python: %s", venv_python)
        return str(venv_python), ".venv.py311/bin/python"

    return None, None


def run_job_local(job):
    """
    Execute job via local subprocess.

    Starts the job in a detached subprocess using the most appropriate
    Python interpreter:
    1. ccp4-python (preferred - includes CCP4 environment and site-packages)
    2. Project virtual environment (fallback for development)

    Args:
        job: Job model instance with attributes:
            - id: Job primary key
            - uuid: Job UUID

    Returns:
        dict: Result dictionary with keys:
            - success (bool): True if job started successfully
            - data (dict): Serialized job data (if success)
            - error (str): Error message (if failure)
            - status (int): HTTP status code

    Environment Requirements:
        - ccp4-python on PATH (after sourcing ccp4.setup-sh), OR
        - Project virtual environment with Django dependencies
        - CCP4 environment variables

    Raises:
        No exceptions - all errors returned in result dict
    """
    logger.info("Running job %s in LOCAL mode via subprocess", job.id)

    try:
        # Path: context_run.py -> jobs -> utils -> lib -> ccp4x -> server -> manage.py
        server_dir = pathlib.Path(__file__).parent.parent.parent.parent.parent
        manage_py = str(server_dir / "manage.py")
        project_root = server_dir.parent

        # Find appropriate Python interpreter
        python_interpreter, interpreter_name = _find_python_interpreter(project_root)

        if python_interpreter is None:
            error_msg = (
                "No suitable Python interpreter found. "
                "Either source ccp4.setup-sh to get ccp4-python on PATH, "
                "or create a virtual environment at .venv"
            )
            logger.error(error_msg)
            return {
                "success": False,
                "error": error_msg,
                "status": 500,
            }

        # Inherit current environment (includes CCP4 vars, PYTHONPATH, etc.)
        env = os.environ.copy()

        # Start job in detached process
        subprocess.Popen(
            [
                python_interpreter,
                manage_py,
                "run_job",
                "-ju",
                str(job.uuid),
            ],
            start_new_session=True,
            env=env,
        )

        logger.info(
            "Started job %s (%s) via subprocess using %s",
            job.id, job.uuid, interpreter_name
        )

        return {
            "success": True,
            "data": job,
            "status": 200,
        }

    except Exception as error:
        logger.exception("Failed to start job via subprocess", exc_info=error)
        return {
            "success": False,
            "error": f"Subprocess error: {str(error)}",
            "status": 500,
        }


def run_job_context_aware(job, force_local=False):
    """
    Execute job using environment-appropriate backend.

    Automatically detects execution context and routes to appropriate handler:
    - Azure Mode: Queues job via Azure Service Bus
    - Local Mode: Executes job via subprocess

    This is the main entry point for context-aware job execution.

    Args:
        job: Job model instance
        force_local (bool): If True, forces local execution regardless of environment

    Returns:
        dict: Result dictionary with keys:
            - success (bool): True if job started/queued successfully
            - data (dict): Job instance (if success)
            - error (str): Error message (if failure)
            - status (int): HTTP status code

    Example:
        from ccp4x.lib.context_dependent_run import run_job_context_aware

        # Normal context-aware execution
        result = run_job_context_aware(job)

        # Force local execution
        result = run_job_context_aware(job, force_local=True)

        if result["success"]:
            serializer = JobSerializer(result["data"])
            return Response(serializer.data)
        else:
            return Response(
                {"error": result["error"]},
                status=result["status"]
            )
    """
    if force_local:
        execution_mode = "local"
        logger.info(
            "Forcing local execution for job %s (uuid=%s, task=%s) via force_local=True",
            job.id,
            job.uuid,
            job.task_name,
        )
    else:
        execution_mode = get_execution_mode()

    logger.info(
        "Executing job %s (uuid=%s) in %s mode",
        job.id,
        job.uuid,
        execution_mode.upper(),
    )

    if execution_mode == "azure":
        return run_job_azure(job)
    else:
        return run_job_local(job)
