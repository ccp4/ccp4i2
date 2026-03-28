# Copyright (C) 2025-2026 University of York
# Copyright (C) 2025-2026 Newcastle University
#
# This file is part of CCP4i2.
#
# CCP4i2 is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License version 3,
# modified in accordance with the provisions of the license to address
# the requirements of UK law.
#
# See https://www.ccp4.ac.uk/ccp4license.php for details.
from ccp4i2.db import models
from ccp4i2.core.tasks import get_task_title
from ..plugins.get_plugin import get_job_plugin


def get_what_next(job: models.Job):
    the_job_plugin = get_job_plugin(job)
    if the_job_plugin is None:
        return {"Status": "Failed", "result": []}

    if hasattr(the_job_plugin, "WHATNEXT"):
        result = {
            "Status": "Success",
            "result": [
                {
                    "taskName": taskName,
                    "shortTitle": taskName,
                    "title": get_task_title(taskName),
                }
                for taskName in the_job_plugin.WHATNEXT
                if isinstance(taskName, str)
            ],
        }
    else:
        result = {"Status": "Failed", "result": []}
    return result
