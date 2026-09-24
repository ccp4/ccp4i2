"""Fan-in: fill a pandda_campaign job's DATASETS from the campaign its
project is the parent of.

Design note: docs/pandda-campaign-design.md, section 3.1. The orchestrator
never asks the database which datasets a campaign has; campaign-awareness
lives here, in job construction. The list is prefilled from the same
``collect_pandda_datasets`` the export endpoint uses, the user sees and edits
it before running, and a dataset the default rule misses is added by hand.
"""
import logging
from pathlib import Path
from typing import Dict, List

from ccp4i2.db import models
from ccp4i2.lib.response import Result

logger = logging.getLogger(f"ccp4i2:{__name__}")

TASK = "pandda_campaign"
EDITABLE = (models.Job.Status.UNKNOWN, models.Job.Status.PENDING)


def _parent_campaigns(project):
    return [m.group for m in models.ProjectGroupMembership.objects.filter(
        project=project, type=models.ProjectGroupMembership.MembershipType.PARENT
    ).select_related("group")]


def _listed_labels(job) -> set:
    """Labels already in the job's DATASETS list, read from its parameters."""
    import xml.etree.ElementTree as ET
    for name in ("input_params.xml", "params.xml"):
        path = Path(job.directory) / name
        if not path.is_file():
            continue
        try:
            root = ET.parse(path).getroot()
        except ET.ParseError:
            continue
        return {(n.findtext("DTAG") or "").strip()
                for n in root.findall("ccp4i2_body/inputData/DATASETS/CPanddaDataset")} - {""}
    return set()


def campaign_candidates(job) -> Dict:
    """What a fill would add, and why anything is left out. Read-only.

    Returns ``campaigns`` (the campaigns this project is parent of),
    ``candidates`` (datasets with a finished dimple, flagged ``listed`` when
    already in the job), ``skipped`` (members with a reason) and ``reason``
    (why nothing can be filled at all, or None).
    """
    from ccp4i2.wrappers.pandda_campaign.script.pandda_staging import campaign_dataset_specs

    result = {"campaigns": [], "candidates": [], "skipped": [], "reason": None}
    if job.task_name != TASK:
        result["reason"] = f"job {job.number} is a {job.task_name} job, not {TASK}"
        return result
    if job.status not in EDITABLE:
        result["reason"] = f"job {job.number} is {job.get_status_display().lower()}; only a pending job can be filled"
        return result
    campaigns = _parent_campaigns(job.project)
    if not campaigns:
        result["reason"] = f"{job.project.name} is not the parent project of any campaign"
        return result

    listed = _listed_labels(job)
    for group in campaigns:
        skipped: list = []
        specs = campaign_dataset_specs(group, skipped=skipped)
        result["campaigns"].append({
            "name": group.name, "id": group.id,
            # uuid arrives with the campaign-persistence change (#608)
            "uuid": str(group.uuid) if getattr(group, "uuid", None) else None,
            "members": group.memberships.filter(
                type=models.ProjectGroupMembership.MembershipType.MEMBER).count(),
        })
        for spec in specs:
            result["candidates"].append({
                "label": spec.label,
                "project_uuid": spec.project_uuid,
                "source_job_uuid": spec.source_job_uuid,
                "xyzin": str(spec.xyzin),
                "hklin": str(spec.hklin),
                "dict": str(spec.dictionary) if spec.dictionary else None,
                "listed": spec.label in listed,
            })
        result["skipped"].extend({"project": name, "reason": why} for name, why in skipped)
    if not result["candidates"]:
        result["reason"] = "no member of the campaign has a finished dimple job yet"
    return result


def fill_datasets(plugin, job) -> Result[Dict]:
    """Append every candidate not already listed to ``plugin``'s DATASETS
    and save the job's parameters. Idempotent by label. ``plugin`` is the
    job's loaded plugin (the generic object_method endpoint hands its own),
    so the container written is the one the interface is showing."""
    preview = campaign_candidates(job)
    to_add = [c for c in preview["candidates"] if not c["listed"]]
    if preview["reason"] and not preview["candidates"]:
        return Result.fail(preview["reason"], details=preview)
    if not to_add:
        return Result.ok({**preview, "added": []})

    datasets = plugin.container.inputData.DATASETS
    added: List[str] = []
    for candidate in to_add:
        item = datasets.makeItem()
        item.DTAG.set(candidate["label"])
        item.XYZIN.setFullPath(candidate["xyzin"])
        item.HKLIN.setFullPath(candidate["hklin"])
        if candidate["dict"]:
            item.DICT.setFullPath(candidate["dict"])
        if candidate["project_uuid"]:
            item.PROJECT_UUID.set(candidate["project_uuid"])
        if candidate["source_job_uuid"]:
            item.SOURCE_JOB_UUID.set(candidate["source_job_uuid"])
        datasets.append(item)
        added.append(candidate["label"])

    params_file = Path(plugin.workDirectory) / "input_params.xml"
    error = plugin.saveDataToXml(str(params_file))
    if error and hasattr(error, "hasError") and error.hasError():
        return Result.fail(f"could not save {params_file}: {error}")
    logger.info("pandda_campaign job %s: filled %d dataset(s) from campaign", job.number, len(added))
    for candidate in preview["candidates"]:
        if candidate["label"] in added:
            candidate["listed"] = True
    return Result.ok({**preview, "added": added})


def fill_datasets_from_campaign(job) -> Result[Dict]:
    """``fill_datasets`` for a caller that has only the job row."""
    from ccp4i2.lib.utils.plugins.plugin_context import get_plugin_with_context

    plugin_result = get_plugin_with_context(job)
    if not plugin_result.success:
        return Result.fail(f"could not load job {job.number}: {plugin_result.error}")
    return fill_datasets(plugin_result.data, job)
