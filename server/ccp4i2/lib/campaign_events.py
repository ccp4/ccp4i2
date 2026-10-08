"""``CampaignEvent``: writing the projection, and reading the site matrix from it.

The rules (which event is at which site, when a dataset's frame cannot be
trusted, which job holds the current model) are pure and live in
:mod:`ccp4i2.lib.campaign_matrix`. This module is the database side:

* :func:`record_receipt` projects one ``pandda_events`` receipt's events into
  rows. It is called from a ``post_save`` receiver when a receipt reaches a
  status whose outputs are recorded (``db/signals.py``), and by the
  ``backfill_campaign_events`` command. Re-recording replaces the receipt's
  rows, so either may run any number of times.
* :func:`site_matrix` builds the per-dataset cells the campaign overview
  (``member_projects``) serves, in a fixed number of queries and with no
  per-dataset file reads.

Nothing here writes a ``SiteEvaluation``: an event is a machine observation,
a verdict is a person's (docs/pandda-campaign-design.md, section 9.1).
"""
import logging
import uuid as uuid_module
from functools import lru_cache
from pathlib import Path
from typing import Dict, Iterable, List, Optional

from django.db import transaction
from django.db.models import Q

from ..db import models
from .campaign_matrix import (
    DIMPLE_TASKS,
    REFINEMENT_TASKS,
    choose_current_model,
    frame_mismatch,
    match_events_to_sites,
    read_model_cell,
)
from .utils.jobs.pandda_site_index import RECEIPT_TASK, read_receipt, receipt_body

logger = logging.getLogger(f"ccp4i2:{__name__}")

#: The statuses at which a receipt's outputs have been recorded: the same set
#: the gleaner publishes for (``async_db_handler.GLEANING_JOB_STATUSES``). A
#: short receipt is UNSATISFACTORY and its events are still real (section 7.2).
RECORDED_STATUSES = frozenset(
    {models.Job.Status.FINISHED, models.Job.Status.UNSATISFACTORY}
)

#: Coordinate file types, in the order campaign_scene prefers them.
_COORD_TYPES = ("chemical/x-pdb", "chemical/x-cif", "chemical/x-mmcif")


# ---------------------------------------------------------------------------
# Writing
# ---------------------------------------------------------------------------

def _parse_uuid(text: Optional[str]):
    if not text:
        return None
    try:
        return uuid_module.UUID(str(text))
    except (TypeError, ValueError):
        return None


def resolve_campaign(job: models.Job, run_job_uuid) -> Optional[models.ProjectGroup]:
    """The campaign a receipt was made for, or None when that is ambiguous.

    A receipt records its run job, not its campaign. The run job lives in the
    parent project of the campaign it ran over, so the campaign that has that
    project as parent and this receipt's project as a member is the answer.
    Failing that (a run produced elsewhere has no run job), a project that is
    a member of exactly one campaign settles it.
    """
    member_of = list(
        models.ProjectGroupMembership.objects.filter(
            project_id=job.project_id,
            type=models.ProjectGroupMembership.MembershipType.MEMBER,
        ).values_list("group_id", flat=True)
    )
    if not member_of:
        return None
    if run_job_uuid is not None:
        run_project = (
            models.Job.objects.filter(uuid=run_job_uuid)
            .values_list("project_id", flat=True).first()
        )
        if run_project is not None:
            groups = list(
                models.ProjectGroupMembership.objects.filter(
                    project_id=run_project,
                    type=models.ProjectGroupMembership.MembershipType.PARENT,
                    group_id__in=member_of,
                ).values_list("group_id", flat=True)
            )
            if len(groups) == 1:
                return models.ProjectGroup.objects.filter(pk=groups[0]).first()
    if len(member_of) == 1:
        return models.ProjectGroup.objects.filter(pk=member_of[0]).first()
    return None


def receipt_cell(job: models.Job) -> Optional[List[float]]:
    """The dataset's unit cell, from the receipt's apo model (``XYZIN_APO``).

    The apo model is the ``-pandda-input.pdb`` the run started from, in the
    dataset's own frame -- the frame the event centroids are stated in. Read
    once, when the receipt is recorded, so no request ever opens it.
    """
    candidates = []
    for row in (models.File.objects.filter(job=job, job_param_name="XYZIN_APO")
                .select_related("job__project").order_by("-id")):
        candidates.append(row.path)
    body = receipt_body(job)
    if body is not None:
        name = (body.findtext("outputData/XYZIN_APO/baseName") or "").strip()
        if name:
            candidates.append(Path(job.directory) / name)
    for path in candidates:
        if path is not None and Path(path).is_file():
            cell = read_model_cell(path)
            if cell is not None:
                return cell
    return None


def _rows(job: models.Job, record: dict, group, cell) -> List[models.CampaignEvent]:
    run_uuid = _parse_uuid(record.get("run_job_uuid"))
    cell = cell or [None] * 6
    rows = []
    seen = set()
    for event in record.get("events", []):
        idx = event.get("event_idx")
        if idx is None:
            logger.warning("Receipt job %s: an event without an EVENT_IDX was not "
                           "recorded", job.number)
            continue
        key = (event.get("site_idx"), idx)
        if key in seen:
            logger.warning("Receipt job %s: event %s at site %s appears twice; "
                           "the first is recorded", job.number, idx, key[0])
            continue
        seen.add(key)
        centroid = event.get("centroid") or [None, None, None]
        rows.append(models.CampaignEvent(
            receipt=job,
            project_id=job.project_id,
            group=group,
            run_job_uuid=run_uuid,
            dtag=(record.get("dtag") or "")[:255],
            event_idx=idx,
            site_idx=event.get("site_idx"),
            centroid_x=centroid[0], centroid_y=centroid[1], centroid_z=centroid[2],
            hit_probability=event.get("hit_probability"),
            has_pose=bool(event.get("has_pose")),
            has_map=bool(event.get("has_map")),
            cell_a=cell[0], cell_b=cell[1], cell_c=cell[2],
            cell_alpha=cell[3], cell_beta=cell[4], cell_gamma=cell[5],
        ))
    return rows


def forget_receipt(job: models.Job) -> int:
    """Drop a receipt's rows. Returns how many there were."""
    deleted, _ = models.CampaignEvent.objects.filter(receipt=job).delete()
    return deleted


def record_receipt(job: models.Job) -> int:
    """Project one receipt's events into ``CampaignEvent``, replacing any rows
    it already had. Returns the number of rows written.

    A receipt not in a recorded status (rerunning, failed, marked for
    deletion) keeps no rows: its ``params.xml`` is not, or no longer, a
    finished record. One savepoint around the lot, so a failure here never
    poisons a transaction the caller is in -- the status update that fired
    the signal, say.
    """
    if job.task_name != RECEIPT_TASK:
        return 0
    with transaction.atomic():
        if job.status not in RECORDED_STATUSES:
            forget_receipt(job)
            return 0
        record = read_receipt(job)
        if record is None:
            forget_receipt(job)
            return 0
        run_uuid = _parse_uuid(record.get("run_job_uuid"))
        group = resolve_campaign(job, run_uuid)
        cell = receipt_cell(job) if record.get("events") else None
        rows = _rows(job, record, group, cell)
        forget_receipt(job)
        models.CampaignEvent.objects.bulk_create(rows)
    return len(rows)


def backfill(group: Optional[models.ProjectGroup] = None) -> Dict[str, int]:
    """(Re)record every receipt, or every receipt of one campaign's members.

    Idempotent: each receipt's rows are replaced. Returns counts of receipts
    read, rows written, and receipts that could not be read.
    """
    jobs = (models.Job.objects.filter(task_name=RECEIPT_TASK)
            .exclude(status=models.Job.Status.TO_DELETE)
            .select_related("project").order_by("id"))
    if group is not None:
        members = group.memberships.filter(
            type=models.ProjectGroupMembership.MembershipType.MEMBER
        ).values_list("project_id", flat=True)
        jobs = jobs.filter(project_id__in=list(members))
    counts = {"receipts": 0, "events": 0, "failed": 0}
    for job in jobs:
        counts["receipts"] += 1
        try:
            counts["events"] += record_receipt(job)
        except Exception:      # noqa: BLE001 - one bad receipt must not stop the rest
            counts["failed"] += 1
            logger.exception("Could not record events of receipt job %s in %s",
                             job.number, job.project.name)
    return counts


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------

def current_model_job(project) -> Optional[models.Job]:
    """The job holding this project's current model (``choose_current_model``)."""
    candidates = models.Job.objects.filter(
        project=project, status=models.Job.Status.FINISHED,
    ).filter(
        Q(task_name__in=REFINEMENT_TASKS, parent__isnull=True)
        | Q(task_name__in=DIMPLE_TASKS)
    )
    return choose_current_model(candidates)


def latest_receipt_ids(jobs: Iterable) -> Dict[int, int]:
    """``{project_id: receipt job id}``: each project's newest recorded receipt.

    "Newest" is the highest job id, which is the most recent fan-out into
    that project -- normally the most recent run.
    """
    latest: Dict[int, int] = {}
    for job in jobs:
        if job.task_name != RECEIPT_TASK or job.status not in RECORDED_STATUSES:
            continue
        if job.id > latest.get(job.project_id, 0):
            latest[job.project_id] = job.id
    return latest


@lru_cache(maxsize=64)
def _cached_cell(path: str, mtime_ns: int, size: int):
    cell = read_model_cell(path)
    return tuple(cell) if cell else None


def file_cell(path) -> Optional[List[float]]:
    """``read_model_cell``, remembered while the file is unchanged."""
    try:
        stat = Path(path).stat()
    except OSError:
        return None
    cell = _cached_cell(str(path), stat.st_mtime_ns, stat.st_size)
    return list(cell) if cell else None


def parent_cell(group) -> Optional[List[float]]:
    """The unit cell of the campaign's reference model (its parent's).

    The same file campaign_scene draws as the reference: the parent's XYZOUT
    PDB, else its XYZOUT, else any coordinate file, newest first. One query,
    and the header read is cached on the file's mtime.
    """
    parent_id = (
        group.memberships.filter(
            type=models.ProjectGroupMembership.MembershipType.PARENT)
        .values_list("project_id", flat=True).first()
    )
    if parent_id is None:
        return None
    files = list(
        models.File.objects.filter(job__project_id=parent_id,
                                   type__name__in=_COORD_TYPES)
        .select_related("job__project", "type").order_by("-id")
    )

    def tier(f):
        if f.job_param_name == "XYZOUT" and f.type.name == "chemical/x-pdb":
            return 0
        return 1 if f.job_param_name == "XYZOUT" else 2

    for f in sorted(files, key=lambda f: (tier(f), -f.id)):
        path = f.path
        if path is not None and path.exists():
            return file_cell(path)
    return None


def _job_ref(job) -> Optional[dict]:
    if job is None:
        return None
    return {"id": job.id, "uuid": str(job.uuid), "number": job.number,
            "task_name": job.task_name}


def site_matrix(group, sites: List[models.CampaignSite],
                jobs_by_project: Dict[int, list],
                verdicts: Dict[tuple, str]) -> Dict[int, dict]:
    """Per member project: ``site_cells``, ``frame_mismatch`` and
    ``current_model_job``, as ``member_projects`` serves them.

    ``jobs_by_project`` holds every job of every member (the caller has
    already fetched them for the job list); ``verdicts`` maps
    ``(project_id, site_id)`` to a verdict. Adds two queries for the whole
    page -- the events of each member's latest receipt, and the parent's
    reference file -- and reads at most one file header, cached.
    """
    latest = latest_receipt_ids(j for jobs in jobs_by_project.values() for j in jobs)
    events_by_receipt: Dict[int, List[models.CampaignEvent]] = {}
    if latest:
        for event in models.CampaignEvent.objects.filter(
                receipt_id__in=list(latest.values())):
            events_by_receipt.setdefault(event.receipt_id, []).append(event)

    reference_cell = parent_cell(group) if events_by_receipt else None

    site_dicts = [{"uuid": str(s.uuid), "origin": s.origin, "radius": s.radius}
                  for s in sites]
    out: Dict[int, dict] = {}
    for project_id, jobs in jobs_by_project.items():
        events = events_by_receipt.get(latest.get(project_id), [])
        cell = next((e.cell for e in events if e.cell), None)
        mismatch = frame_mismatch(cell, reference_cell) if events else None
        matched = match_events_to_sites(
            ({"event_idx": e.event_idx, "centroid": e.centroid,
              "hit_probability": e.hit_probability, "has_pose": e.has_pose}
             for e in events),
            site_dicts,
            mismatch=mismatch,
        )
        out[project_id] = {
            "site_cells": {
                str(site.uuid): {
                    "event": matched[str(site.uuid)],
                    "verdict": verdicts.get((project_id, site.id)),
                }
                for site in sites
            },
            "frame_mismatch": mismatch,
            "current_model_job": _job_ref(choose_current_model(jobs)),
        }
    return out
