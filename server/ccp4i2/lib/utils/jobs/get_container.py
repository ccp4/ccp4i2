"""
Utility to get job container from database.

Loads job parameters and returns the CContainer for manipulation.
"""

import logging
from ccp4i2.core.CCP4Container import CContainer
from ccp4i2.db import models
from ..containers.get_container import get_job_container as _legacy_get_job_container

logger = logging.getLogger(f"ccp4i2:{__name__}")


def get_job_container(job: models.Job) -> CContainer:
    """
    Get CContainer for a job.

    Loads the job's plugin and returns its container with parameters loaded.

    Args:
        job: Job model instance

    Returns:
        CContainer instance with job parameters

    Raises:
        FileNotFoundError: If params file not found
        Exception: If loading fails
    """
    # Delegate to existing implementation
    return _legacy_get_job_container(job)
