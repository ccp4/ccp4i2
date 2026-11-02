"""Where a PanDDA run's events are looked at on the campaign page.

A run that belongs to a campaign (its project is the campaign's parent)
has, for each dataset, a member project and that member's current model.
The report's "open in Moorhen" buttons open the campaign viewer on that
model, with the event scene applied over it. The plugin records these in
program.xml as it writes it; a report drawn from an older record (a run
before this existed, or the plugin not being the process that wrote it)
looks them up itself. Both come here.

Reads the database; never raises. A run outside a campaign gives an empty
answer, and the report then opens events on the project page.
"""
import logging
from typing import Dict, Optional, Tuple

logger = logging.getLogger(f"ccp4i2:{__name__}")


def campaign_views(job, manifest) -> Dict:
    """``{"campaign": <ProjectGroup id or None>, "members": {dtag: (project id,
    current-model job id or None)}}`` for the run ``job`` (a database row)
    described by ``manifest`` (the run's manifest.json, parsed)."""
    views: Dict = {"campaign": None, "members": {}}
    try:
        from ccp4i2.db import models
        from ccp4i2.lib.campaign_events import current_model_job
        from ccp4i2.lib.utils.jobs.pandda_fanin import _parent_campaigns
        if job is None or not manifest:
            return views
        campaigns = _parent_campaigns(job.project)
        if not campaigns:
            return views
        views["campaign"] = campaigns[0].id
        for entry in manifest.get("datasets", []):
            uuid = entry.get("project_uuid")
            if not uuid:
                continue
            project = models.Project.objects.filter(uuid=uuid).first()
            if project is None:
                continue
            current = current_model_job(project)
            views["members"][entry["xtal"]] = (project.id, current.id if current else None)
    except Exception as e:      # noqa: BLE001 - a report is never worth a job
        logger.debug("campaign views not read: %s", e)
    return views


def campaign_views_for_job_uuid(job_uuid: Optional[str]) -> Dict:
    """``campaign_views`` for a job known only by uuid, its manifest read from
    its directory: what a report drawn later has to go on."""
    import json
    from pathlib import Path
    try:
        from ccp4i2.db import models
        job = models.Job.objects.filter(uuid=job_uuid).first() if job_uuid else None
        if job is None:
            return {"campaign": None, "members": {}}
        manifest_path = Path(job.directory) / "manifest.json"
        manifest = json.loads(manifest_path.read_text()) if manifest_path.is_file() else None
        return campaign_views(job, manifest)
    except Exception as e:      # noqa: BLE001
        logger.debug("campaign views not read for %s: %s", job_uuid, e)
        return {"campaign": None, "members": {}}


def member_view(views: Dict, dtag: str) -> Tuple[Optional[int], Optional[int]]:
    """``(campaign id, member job id)`` to open ``dtag``'s event on, or
    ``(None, None)`` when it opens on the project page."""
    member = views.get("members", {}).get(dtag) if views.get("campaign") is not None else None
    if member and member[1] is not None:
        return views["campaign"], member[1]
    return None, None
