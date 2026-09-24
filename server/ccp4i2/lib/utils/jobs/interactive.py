"""Interactive sessions: jobs whose "program" is a window in the app.

The recorded Moorhen task is the first. Its session is a database row
(``JobInteractiveSession``), not a process:

* ``open_session`` -- what Run does for an interactive task instead of
  dispatching: create the row, set the job RUNNING, return. The app opens
  the session route; no process exists.
* ``heartbeat`` / ``drop_file`` -- what the window does while open. Saved
  files land in the job's drop directory under the ``output<N>`` contract
  the Coot wrappers share, with a ``.meta.json`` sidecar for the annotation.
* ``finish_session`` -- the only thing that ends a session (the user, from
  the window or the job menu). If anything was saved the job is dispatched
  through the ordinary runner, whose plugin harvests the drop directory and
  gleans; if nothing was, the job is marked for deletion without dispatch.
  A job that i2run dispatched first is already owned by a waiting plugin
  (``dispatched``), so finishing only marks the row.

No timer ends a session on the desktop: a laptop asleep with a window open
must find its session still there on waking. See
docs/moorhen-task-design.md.
"""

import json
import logging
from datetime import timedelta
from pathlib import Path

from django.utils import timezone

from ccp4i2.core.tasks import is_interactive
from ccp4i2.cootbridge import api_client
from ccp4i2.db import models

logger = logging.getLogger(f"ccp4i2:{__name__}")

DROP_DIR_NAME = "MOORHEN_FILE_DROP"
#: A window that has sent a heartbeat this recently counts as attached.
ATTACHED_WINDOW_SECONDS = 30
_COORDINATE_SUFFIXES = {".pdb": "pdb", ".ent": "pdb", ".cif": "cif", ".mmcif": "cif"}


class SessionError(Exception):
    """A request the session cannot honour, with the HTTP status it earns."""

    def __init__(self, status, message):
        super().__init__(message)
        self.status = status


def is_interactive_job(job):
    return is_interactive(job.task_name)


def drop_dir_for(job) -> Path:
    return Path(job.directory) / DROP_DIR_NAME


def _session(job):
    try:
        return job.interactive_session
    except models.JobInteractiveSession.DoesNotExist:
        return None


# ---------------------------------------------------------------------------
# Opening
# ---------------------------------------------------------------------------

def open_session(job):
    """Run, for an interactive task: open the session instead of dispatching.

    Idempotent for a job whose session is already open (reopening the
    window). Refuses a job that is queued or running under a real process.
    """
    session = _session(job)
    if session is not None and not session.finished \
            and job.status == models.Job.Status.RUNNING:
        return job
    if job.status in (models.Job.Status.QUEUED,
                      models.Job.Status.RUNNING,
                      models.Job.Status.RUNNING_REMOTELY):
        raise SessionError(409, f"Job is already {job.get_status_display().lower()}")
    drop_dir_for(job).mkdir(parents=True, exist_ok=True)
    models.JobInteractiveSession.objects.update_or_create(
        job=job,
        defaults={
            "requested_at": timezone.now(),
            "last_heartbeat": None,
            "dispatched": False,
            "finished": False,
            "finished_at": None,
        },
    )
    job.status = models.Job.Status.RUNNING
    job.process_id = None
    job.save()
    logger.info("Opened interactive session for job %s (%s)", job.id, job.task_name)
    return job


def ensure_session_for_runner(job):
    """Called by the waiting plugin: the runner owns this job now.

    A job dispatched first (i2run) has no session yet; create one so a
    window can attach. Either way, record that a plugin is waiting, so
    finishing must not dispatch again.
    """
    drop_dir_for(job).mkdir(parents=True, exist_ok=True)
    session, _created = models.JobInteractiveSession.objects.get_or_create(job=job)
    if not session.dispatched:
        session.dispatched = True
        session.save(update_fields=["dispatched"])
    return session


def require_open_session(job):
    session = _session(job)
    if session is None:
        raise SessionError(409, "Job has no interactive session")
    if session.finished:
        raise SessionError(409, "Interactive session is finished")
    if job.status != models.Job.Status.RUNNING:
        raise SessionError(409, f"Job is not running (status: {job.get_status_display()})")
    return session


# ---------------------------------------------------------------------------
# While open
# ---------------------------------------------------------------------------

def heartbeat(job):
    session = require_open_session(job)
    session.last_heartbeat = timezone.now()
    session.save(update_fields=["last_heartbeat"])
    return session_state(job)


def _extension_for(kind, filename, head):
    if kind == "dictionary":
        return "cif"
    suffix = Path(filename or "").suffix.lower()
    if suffix in _COORDINATE_SUFFIXES:
        return _COORDINATE_SUFFIXES[suffix]
    return "cif" if head.lstrip().startswith(b"data_") else "pdb"


def drop_file(job, source, filename="", kind="model", annotation=""):
    """Save a file into the session's drop directory.

    ``source`` is a Django UploadedFile (has ``chunks()``) or a path. The
    file becomes ``output<N>.<pdb|cif>`` -- the contract every GUI task's
    harvest reads -- with ``output<N>.meta.json`` beside it carrying the
    kind and annotation for the harvest to apply.
    """
    require_open_session(job)
    if kind not in ("model", "dictionary"):
        raise SessionError(400, f"Unknown kind {kind!r}; expected model or dictionary")
    drop_dir = drop_dir_for(job)
    drop_dir.mkdir(parents=True, exist_ok=True)

    if hasattr(source, "chunks"):
        chunks = list(source.chunks())
        filename = filename or getattr(source, "name", "")
    else:
        with open(source, "rb") as handle:
            chunks = [handle.read()]
        filename = filename or Path(source).name
    head = chunks[0][:64] if chunks else b""

    extension = _extension_for(kind, filename, head)
    number = api_client.next_output_number(str(drop_dir))
    target = Path(api_client.output_path(str(drop_dir), number, extension))
    with open(target, "wb") as handle:
        for chunk in chunks:
            handle.write(chunk)
    meta = {
        "kind": kind,
        "annotation": annotation or "",
        "original_name": filename,
        "saved_at": timezone.now().isoformat(),
    }
    with open(target.with_suffix(".meta.json"), "w") as handle:
        json.dump(meta, handle)
    logger.info("Session drop for job %s: %s (%s)", job.id, target.name, kind)
    return {"number": number, "name": target.name, **meta}


def session_outputs(job):
    """What the window has saved so far, in save order."""
    from ccp4i2.cootbridge.harvest import read_drop_metadata

    outputs = []
    for number, path in api_client.harvestable_outputs(str(drop_dir_for(job))):
        meta = read_drop_metadata(path)
        outputs.append({
            "number": number,
            "name": Path(path).name,
            "kind": meta.get("kind", "model"),
            "annotation": meta.get("annotation", ""),
        })
    return outputs


# ---------------------------------------------------------------------------
# Ending
# ---------------------------------------------------------------------------

def finish_session(job, finished=True):
    """End the session, or record that a window detached.

    Returns the session state plus a ``disposition``:

    * ``dispatched`` -- saved files exist; the job was handed to the runner.
    * ``harvesting`` -- saved files exist; a waiting plugin owns the job
      and will harvest now that the row is marked.
    * ``deleted`` -- nothing was saved; the job is marked for deletion.
    * ``kept_open`` -- ``finished=False`` (window closed) with saved files;
      the session stays open for reconnect or Finish from the job menu.
    * ``already_finished``.
    """
    session = _session(job)
    if session is None:
        raise SessionError(409, "Job has no interactive session")
    if session.finished:
        return {**session_state(job), "disposition": "already_finished"}

    outputs = session_outputs(job)
    if not finished and outputs:
        return {**session_state(job), "disposition": "kept_open"}

    session.finished = True
    session.finished_at = timezone.now()
    session.save(update_fields=["finished", "finished_at"])

    if not outputs:
        job.status = models.Job.Status.TO_DELETE
        job.finish_time = timezone.now()
        job.save()
        logger.info("Interactive session of job %s finished with nothing saved", job.id)
        return {**session_state(job), "disposition": "deleted"}

    if session.dispatched:
        return {**session_state(job), "disposition": "harvesting"}

    session.dispatched = True
    session.save(update_fields=["dispatched"])
    from .context_run import run_job_context_aware

    result = run_job_context_aware(job, force_dispatch=True)
    if not result.get("success"):
        raise SessionError(result.get("status", 500), result.get("error", "dispatch failed"))
    job.refresh_from_db()
    return {**session_state(job), "disposition": "dispatched"}


def cancel_session(job):
    """Cancel: the session is over and nothing is harvested."""
    session = _session(job)
    if session is not None and not session.finished:
        session.finished = True
        session.finished_at = timezone.now()
        session.save(update_fields=["finished", "finished_at"])


# ---------------------------------------------------------------------------
# State and the load plan
# ---------------------------------------------------------------------------

def session_state(job):
    session = _session(job)
    attached = False
    if session is not None and session.last_heartbeat is not None:
        attached = (timezone.now() - session.last_heartbeat
                    <= timedelta(seconds=ATTACHED_WINDOW_SECONDS))
    return {
        "job_id": job.id,
        "task_name": job.task_name,
        "status": job.status,
        "session": None if session is None else {
            "requested_at": session.requested_at,
            "last_heartbeat": session.last_heartbeat,
            "attached": attached,
            "dispatched": session.dispatched,
            "finished": session.finished,
            "finished_at": session.finished_at,
        },
        "load_plan": load_plan_for_job(job),
        "outputs": session_outputs(job),
    }


def _params_xml_text(job):
    for name in ("params.xml", "input_params.xml"):
        candidate = Path(job.directory) / name
        if candidate.is_file():
            return candidate.read_text()
    return None


def load_plan_for_job(job):
    """What the window loads at startup, in load order.

    Derived from the job's input parameters by the same walk the Coot
    bridge uses (``api_client.parse_input_params``: dictionaries first,
    then coordinates, then maps), resolved to File rows so the window can
    fetch each by id. An input with no File row yet (a path not imported
    until dispatch) is reported with ``file_id`` null.
    """
    xml_text = _params_xml_text(job)
    if not xml_text:
        return []
    plan = []
    for entry in api_client.parse_input_params(xml_text):
        file_row = None
        db_file_id = (entry.get("db_file_id") or "").strip()
        if db_file_id:
            try:
                file_row = models.File.objects.select_related("type").get(uuid=db_file_id)
            except (models.File.DoesNotExist, ValueError):
                file_row = None
        sub_type = None
        if file_row is not None and file_row.sub_type is not None:
            sub_type = file_row.sub_type
        elif entry.get("sub_type"):
            try:
                sub_type = int(entry["sub_type"])
            except ValueError:
                sub_type = None
        plan.append({
            "kind": entry["kind"],
            "param": entry["param"],
            "file_id": file_row.id if file_row is not None else None,
            "file_uuid": str(file_row.uuid) if file_row is not None else (db_file_id or None),
            "name": file_row.name if file_row is not None else entry["base_name"],
            "type": file_row.type_id if file_row is not None else None,
            "sub_type": sub_type,
            "label": entry.get("annotation")
                     or (file_row.annotation if file_row is not None else "")
                     or entry["base_name"],
        })
    return plan
