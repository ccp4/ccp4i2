import logging

from ccp4i2.db import models
from ccp4i2.core.tasks import TASKS, get_plugin_class
from ..plugins.get_plugin import get_job_plugin

logger = logging.getLogger(__name__)


def get_what_next(job: models.Job):
    the_job_plugin = get_job_plugin(job)
    if the_job_plugin is None:
        return {"Status": "Failed", "result": []}

    if not hasattr(the_job_plugin, "WHATNEXT"):
        return {"Status": "Failed", "result": []}

    names = [name for name in the_job_plugin.WHATNEXT if isinstance(name, str)]
    return {
        "Status": "Success",
        "result": [_describe(name) for name in names if _can_run_here(name)],
    }


def _can_run_here(taskName: str) -> bool:
    """False only for a task whose required program is certainly missing: one
    the pre-run check would refuse to run. Coot 0.9 (coot_rebuild) is not in
    CCP4 10 and Coot 1 (coot1) is not in CCP4 9, and a list offering both
    should offer the one that is there (#838).

    Where this process cannot see the host that runs jobs -- a remote job
    target, the slim CCP4-free server -- "not found" is ignorance, so every
    suggestion stays, as the pre-run check stays advisory there.
    """
    try:
        from ccp4i2.lib.utils.jobs.context_run import program_checks_are_authoritative
        if not program_checks_are_authoritative():
            return True
        from ccp4i2.config.program_discovery import resolve_program
        required, _ = get_plugin_class(taskName).programsByNeed()
        return all(resolve_program(name) is not None for name in required)
    except Exception:
        logger.debug("Could not check the programs of %s; offering it", taskName, exc_info=True)
        return True


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
