# Copyright (C) 2025 University of York
#
# This file is part of CCP4i2.
#
# CCP4i2 is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License version 3,
# modified in accordance with the provisions of the license to address
# the requirements of UK law.
#
# See https://www.ccp4.ac.uk/ccp4license.php for details.
import uuid
import logging
from .create_job import create_job
from .set_input_by_context_job import set_input_by_context_job
from ccp4i2.db import models

logger = logging.getLogger(f"ccp4i2:{__name__}")


def create_task(the_project: models.Project, arg: any):
    """
    Creates a new task within the given project.

    Args:
        the_project (models.Project): The project in which the task will be created.
        arg (dict): A dictionary containing the task details. Must include 'task_name'.

    Returns:
        models.Job: The newly created job object.

    Raises:
        KeyError: If 'task_name' is not found in arg.
        models.Job.DoesNotExist: If the job with the created UUID does not exist.
    """

    task_name = arg["task_name"]
    new_job_uuid = create_job(projectId=the_project.uuid, taskName=task_name)
    context_job_uuid = arg.get("context_job_uuid", None)
    new_job = models.Job.objects.get(uuid=new_job_uuid)

    context_job = None
    if context_job_uuid is not None:
        try:
            context_job = models.Job.objects.get(uuid=context_job_uuid)
        except models.Job.DoesNotExist:
            logger.warning("Context job %s does not exist", str(context_job_uuid))
            context_job = None
    else:
        logger.debug("No context job provided, looking for the last one")
        context_job = (
            models.Job.objects.filter(parent__isnull=True, project=the_project)
            .exclude(id=new_job.id)
            .last()
        )
        logger.debug("Context job is %s", str(context_job))

    if context_job is not None:
        logger.debug("Context job selected is %s", context_job.number)
        set_input_by_context_job(str(new_job_uuid), str(context_job.uuid))

    return new_job
