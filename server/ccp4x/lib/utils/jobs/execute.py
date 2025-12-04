"""
Job execution utilities (context-aware: local or remote).

Provides environment-aware job execution (local or remote).
"""

import logging
from ccp4x.db import models
from ccp4x.lib.response import Result
from .context_run import run_job_context_aware

logger = logging.getLogger(f"ccp4x:{__name__}")


def execute_job(job: models.Job, force_local: bool = False) -> Result[models.Job]:
    """
    Execute a job using environment-appropriate backend.

    Automatically adapts to deployment context:
    - Local Mode: Executes job via subprocess
    - Azure Mode: Queues job via Service Bus

    Args:
        job: Job model instance to execute
        force_local: If True, force local execution even in Azure

    Returns:
        Result containing updated job instance

    Example:
        >>> result = execute_job(job)
        >>> if result.success:
        ...     updated_job = result.data
        ...     print(f"Job {updated_job.number} started")
    """
    try:
        # Use existing context-dependent runner
        result_dict = run_job_context_aware(job, force_local=force_local)

        if result_dict["success"]:
            logger.info(
                "Started job %s (task: %s, mode: %s)",
                job.uuid, job.task_name,
                "local" if force_local else "auto"
            )
            return Result.ok(result_dict["data"])
        else:
            return Result.fail(
                result_dict["error"],
                details={"job_id": str(job.uuid), "status_code": result_dict.get("status")}
            )

    except Exception as err:
        logger.exception("Failed to run job %s", job.uuid, exc_info=err)
        return Result.fail(
            f"Failed to start job: {str(err)}",
            details={"job_id": str(job.uuid), "error_type": type(err).__name__}
        )
