from ccp4i2.db import models
from ccp4i2.core.tasks import TASKS
from ..plugins.get_plugin import get_job_plugin


def get_what_next(job: models.Job):
    the_job_plugin = get_job_plugin(job)
    if the_job_plugin is None:
        return {"Status": "Failed", "result": []}

    if not hasattr(the_job_plugin, "WHATNEXT"):
        return {"Status": "Failed", "result": []}

    return {
        "Status": "Success",
        "result": [
            _describe(taskName)
            for taskName in the_job_plugin.WHATNEXT
            if isinstance(taskName, str)
        ],
    }


def _describe(taskName: str) -> dict:
    """Human-readable labels for a suggested next task.

    Both labels are offered so the caller can pick one that fits: `title` is
    the full descriptive name, `shortTitle` the compact one used where space is
    tight (the pinned "What next" toolbar). Either may fall back to the raw
    task name for a task with no registry entry.
    """
    task = TASKS.get(taskName)
    short_title = (task.shortTitle if task else None) or taskName
    return {
        "taskName": taskName,
        "shortTitle": short_title,
        "title": (task.title if task else None) or short_title,
    }
