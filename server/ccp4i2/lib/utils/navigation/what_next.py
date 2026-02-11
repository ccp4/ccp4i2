from ccp4i2.db import models
from ccp4i2.core.task_manager.metadata import TITLES
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
                    "title": TITLES.get(taskName),
                }
                for taskName in the_job_plugin.WHATNEXT
            ],
        }
    else:
        result = {"Status": "Failed", "result": []}
    return result
