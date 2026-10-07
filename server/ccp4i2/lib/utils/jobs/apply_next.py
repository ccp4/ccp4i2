"""Act on one of a judgement's next steps: make the job it describes.

A next step is either a follow-on job (a task, with inputs drawn from this
job's outputs) or a rerun (this task again, as a clone with some parameters
changed). The job is made and its inputs set, and left pending for the user
to check and run: the steps come from draft judgements, so a person decides.

References to this job's own task are pinned to this job first
(``crank2[-1].XYZOUT`` -> ``[5].XYZOUT``): the app shows an older job's
judgement and means that job, not the latest of its task.
"""
import logging

from asgiref.sync import async_to_sync

from ccp4i2.db import models

logger = logging.getLogger(f"ccp4i2:{__name__}")

FINISHED = (models.Job.Status.FINISHED, models.Job.Status.FAILED,
            models.Job.Status.UNSATISFACTORY, models.Job.Status.INTERRUPTED)


class ApplyNextError(ValueError):
    pass


def next_steps(job):
    """This job's next steps, from its (kept) verdict, references pinned."""
    from ccp4i2.agent.judgement import judge, judge_finished, pin_references
    verdict = (judge_finished if job.status in FINISHED else judge)(job.task_name, job.directory)
    return pin_references(verdict.get("next") or [], job.task_name, job.number)


def _literal(value):
    """A judgement's input value as the parameter takes it: "True"/"False"
    are flags (a CBoolean set from the string "False" would read as true)."""
    if isinstance(value, str) and value in ("True", "False"):
        return value == "True"
    return value


def _set(new_job, name, value):
    """Set one input, finding its section; {name, value, ok, error}."""
    from ccp4i2.lib.utils.files.resolve_fileuse import is_fileuse_pattern, resolve_fileuse
    from ccp4i2.lib.utils.parameters.set_param import SECTIONS, set_parameter

    shown = value
    if isinstance(value, str) and is_fileuse_pattern(value):
        try:
            resolved = resolve_fileuse(new_job.project, value)
        except Exception as err:  # noqa: BLE001 - reported per input, others still set
            return {"name": name, "value": shown, "ok": False, "error": f"{shown}: {err}"}
        if not resolved.success:
            return {"name": name, "value": shown, "ok": False, "error": f"{shown}: {resolved.error}"}
        value = dict(resolved.data)
        value.pop("fullPath", None)
    else:
        value = _literal(value)
    paths = [name] if "." in name else [f"{section}.{name}" for section in SECTIONS]
    last_error = None
    for path in paths:
        result = set_parameter(new_job, path, value)
        if result.success:
            return {"name": name, "value": shown, "ok": True}
        last_error = result.error
        if "no parameter at" not in str(result.error):
            break  # the parameter exists; the value was refused
    return {"name": name, "value": shown, "ok": False, "error": last_error}


def apply_next(job, index):
    """Make the job for next step ``index`` of ``job``'s judgement.

    Returns {"job": the new job, "rerun": bool, "inputs": [per-input
    outcomes], "advice": the step's advice}.
    """
    steps = next_steps(job)
    if not 0 <= index < len(steps):
        raise ApplyNextError(f"no next step {index}: this judgement offers {len(steps)}")
    step = steps[index]
    if not step.get("task"):
        raise ApplyNextError("this step is advice only: there is no job to make")
    if step.get("rerun"):
        from ccp4i2.lib.utils.jobs.clone import clone_job
        result = clone_job(job.uuid)
        if not result.success:
            raise ApplyNextError(f"could not clone job {job.number}: {result.error}")
        new_job = result.data
    else:
        from ccp4i2.lib.async_create_job import create_job_async
        created = async_to_sync(create_job_async)(
            project_uuid=job.project.uuid, task_name=step["task"], save_params=True,
            context_job_uuid=job.uuid, auto_context=True)
        new_job = models.Job.objects.get(uuid=created["job_uuid"])
    inputs = [_set(new_job, name, value) for name, value in (step.get("inputs") or {}).items()]
    return {"job": new_job, "rerun": bool(step.get("rerun")), "inputs": inputs,
            "advice": " ".join(str(step.get("advice") or "").split())}
