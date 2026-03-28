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
from typing import List
from ccp4i2.db import models
from ..plugins.get_plugin import get_job_plugin


def object_method(
    the_job: models.Job,
    object_path: str,
    method_name: str,
    args: List[str] = None,
    kwargs: dict = None,
):
    """
    Execute a method on a CData object within a job's container.

    Args:
        the_job: The Job model instance
        object_path: Dot-separated path to the object (e.g., "task.inputData.ASU_CONTENT")
        method_name: Name of the method to call on the object
        args: Positional arguments to pass to the method
        kwargs: Keyword arguments to pass to the method

    Returns:
        The return value of the called method
    """
    if args is None:
        args = []
    if kwargs is None:
        kwargs = {}

    the_job_plugin = get_job_plugin(the_job)
    base_element = the_job_plugin.container.find_by_path(object_path, skip_first=True)
    # Call the method with provided args/kwargs
    result = getattr(base_element, method_name)(*args, **kwargs)
    return result
