"""
Remember built campaign scenes while nothing they were built from changes.

A campaign scene (``campaign_scene``) is assembled on every request: several
queries per member, every member's coordinates parsed and fitted onto the
reference, every dictionary read. On a campaign of a few hundred datasets
that is a minute of work for a page that is opened many times a day and
changes a few times a week -- when a refinement finishes, a verdict is
recorded, a site is moved.

So the scene is kept in Django's cache under a key that *is* a digest of
everything the build reads from the database: the group, its memberships,
the jobs and files of its projects, its sites and their verdicts. Nothing
is invalidated by hand; a change to any of those rows changes the digest,
and the next request builds afresh under a new key while the stale entry
ages out. The digest is a handful of aggregate queries, which is cheap next
to the build, and it cannot go stale in a way a signal receiver that
somebody forgot to wire could.

What the digest does not see is a coordinate file rewritten in place under
the same ``File`` row. Jobs write new files under new rows, so in normal use
that does not happen; the ``refresh`` escape hatch (``?refresh=1``) is there
for when it does.

The backend is whatever ``settings.CACHES`` names. Production uses a
file-based cache so an entry survives gunicorn recycling its workers; the
test settings use the dummy backend so one test cannot be served another's
scene, and the cache tests opt in to a real backend explicitly.
"""
import hashlib
import json
import logging
import time
from typing import Callable, Optional

from django.core.cache import cache
from django.db.models import Count, Max

from ..db import models
from . import campaign_matrix

logger = logging.getLogger(f"ccp4i2:{__name__}")

# Bump when the scene format or the build changes in a way that makes an
# entry built by older code wrong for newer readers.
SCENE_CACHE_VERSION = 1
# A day: long enough to carry a campaign through a working day of page
# opens, short enough that an entry orphaned by a change the digest does
# not see (see the module docstring) is not kept for ever.
SCENE_CACHE_TIMEOUT = 24 * 60 * 60

# The job kinds the build selects a member's current model from, and the
# dictionary job it falls back to for a disk-only acedrg CIF. Their status
# and evaluation are listed row by row (a few thousand small tuples at the
# most), because a change of status on an *old* job -- a re-evaluation, a
# job marked failed -- moves no aggregate.
_MODEL_TASKS = tuple(campaign_matrix.REFINEMENT_TASKS) + tuple(campaign_matrix.DIMPLE_TASKS)


def scene_fingerprint(group) -> str:
    """A digest of every database row a scene for ``group`` is built from."""
    memberships = list(
        group.memberships.order_by("project_id").values_list(
            "project_id", "type", "project__name", "project__uuid"
        )
    )
    project_ids = [m[0] for m in memberships]

    jobs = models.Job.objects.filter(project_id__in=project_ids).aggregate(
        count=Count("id"), max_id=Max("id"), finished=Max("finish_time")
    )
    model_jobs = list(
        models.Job.objects.filter(
            project_id__in=project_ids, task_name__in=_MODEL_TASKS
        ).order_by("id").values_list("id", "status", "evaluation", "parent_id")
    )
    files = models.File.objects.filter(job__project_id__in=project_ids).aggregate(
        count=Count("id"), max_id=Max("id")
    )
    sites = list(
        group.site_set.order_by("id").values_list(
            "id", "name", "origin_x", "origin_y", "origin_z",
            "radius", "order", "zoom", "quat",
        )
    )
    verdicts = list(
        models.SiteEvaluation.objects.filter(site__group=group)
        .order_by("id")
        .values_list("id", "project_id", "site_id", "verdict", "evaluator")
    )

    payload = {
        "version": SCENE_CACHE_VERSION,
        "group": [group.id, group.name, group.type],
        "memberships": memberships,
        "jobs": jobs,
        "model_jobs": model_jobs,
        "files": files,
        "sites": sites,
        "verdicts": verdicts,
    }
    text = json.dumps(payload, sort_keys=True, default=str)
    return hashlib.sha1(text.encode("utf-8")).hexdigest()


def _params_key(params: Optional[dict]) -> str:
    if not params:
        return "-"
    return hashlib.sha1(
        json.dumps(params, sort_keys=True, default=str).encode("utf-8")
    ).hexdigest()[:12]


def cached_scene(group, kind: str, build: Callable[[], dict],
                 params: Optional[dict] = None, refresh: bool = False) -> dict:
    """``build()``, or the result of an earlier identical call.

    ``kind`` names the scene (``summary``, ``site``...) and ``params`` the
    request options that change what it draws, so that two scenes of one
    campaign never share an entry. ``refresh`` builds regardless and
    replaces whatever was cached.
    """
    started = time.monotonic()
    fingerprint = scene_fingerprint(group)
    key = (
        f"ccp4i2:campaign-scene:v{SCENE_CACHE_VERSION}:{kind}:"
        f"{group.id}:{fingerprint}:{_params_key(params)}"
    )

    if not refresh:
        try:
            hit = cache.get(key)
        except Exception as exc:  # noqa: BLE001 - a broken cache must not break the page
            logger.warning("Scene cache read failed (%s); building", exc)
            hit = None
        if hit is not None:
            logger.info(
                "Scene cache hit: %s scene of group %s (%.0f ms)",
                kind, group.id, 1000 * (time.monotonic() - started),
            )
            return hit

    result = build()
    try:
        cache.set(key, result, SCENE_CACHE_TIMEOUT)
    except Exception as exc:  # noqa: BLE001 - as above
        logger.warning("Scene cache write failed (%s); scene served unsaved", exc)
    logger.info(
        "Scene cache miss: built %s scene of group %s in %.0f ms",
        kind, group.id, 1000 * (time.monotonic() - started),
    )
    return result
