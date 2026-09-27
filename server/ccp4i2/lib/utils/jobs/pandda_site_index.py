"""A PanDDA run's own sites, and which receipts sit at each.

Two site axes exist and must not be confused.

``CampaignSite`` is the curated one: a durable place with a uuid that
``SiteEvaluation`` rows hang off, which ``campaign_scene.build_site_scene``
draws as the exemplar plus every dataset a person has called a *hit* there.
Membership is the verdict and nothing else. That is the right axis for
annotation and the wrong one for looking at what a run has just found,
because at that moment nobody has voted and every site is empty of verdicts.

This module is the other axis: the sites PanDDA itself found, numbered by the
run that found them. It exists because fan-out is dtag-keyed and site-blind
-- ``plan_fanout`` walks the staging manifest's ``datasets`` and nothing else
-- so once a run has been scattered into per-dataset receipts, the one
grouping worth navigating by survives only as ``SITE_IDX`` on individual
events inside individual receipts, and nothing joined them back up. That is
why a panoptic view across the ligands at one site was hard to get at: the
scatter dissolved the grouping.

**Nothing is stored.** The index is computed from what already exists: each
receipt's own parameters say which run and dataset it is for, where the tree
is, and what its events are, and the tree's sites table says where each site
is. Computing beats recording here, because a recorded index goes stale
against a re-run or a second fan-out, and because this way runs that were
fanned out before any of this existed are navigable too.

PanDDA renumbers sites every run -- site 3 of one run is not site 3 of the
next, and its numbering is a function of which datasets went in -- so a run
site is not a ``CampaignSite`` and cannot be made one by matching numbers.
The link between the axes is *adoption*: a person promotes a run site to a
campaign site. That is the invariant ``CampaignSite`` already states, that an
automated process may propose sites and record observations but may not
delete or move a site an evaluation refers to.
"""
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Optional

from ccp4i2.db import models

RECEIPT_TASK = "pandda_events"

#: Event fields read straight off a receipt, with the parser for each. The
#: receipt is authoritative for these: it recorded them when it ran, from the
#: tree as it was then, and it survives the tree being moved or deleted.
_EVENT_FLOATS = ("BDC", "SCORE", "BUILD_SCORE", "RSCC", "HIT_PROBABILITY",
                 "OPTIMAL_CONTOUR", "DISPLAY_CONTOUR")
_EVENT_INTS = ("EVENT_IDX", "SITE_IDX")


# ---------------------------------------------------------------------------
# Reading one receipt
# ---------------------------------------------------------------------------

def receipt_body(job) -> Optional[ET.Element]:
    """The ``ccp4i2_body`` of a job's parameter file, or None.

    ``params.xml`` is what the job ran with and is rewritten as it runs;
    ``input_params.xml`` is what it was created with. A receipt that has run
    has both and they agree on the inputs, but only ``params.xml`` carries
    the outputs, so it is tried first.
    """
    for name in ("params.xml", "input_params.xml"):
        path = Path(job.directory) / name
        if not path.is_file():
            continue
        try:
            root = ET.parse(path).getroot()
        except ET.ParseError:
            continue
        body = root.find("ccp4i2_body")
        if body is not None:
            return body
    return None


def _number(text, cast):
    if text is None:
        return None
    text = text.strip()
    if not text:
        return None
    try:
        return cast(text)
    except (TypeError, ValueError):
        return None


def _is_set(element: Optional[ET.Element]) -> bool:
    """Whether a file-valued parameter names a file.

    An unset ``CDataFile`` is written as an empty element or one whose
    ``baseName`` is empty, so the element existing proves nothing.
    """
    if element is None:
        return False
    return bool((element.findtext("baseName") or "").strip())


def _centroid(element: Optional[ET.Element]) -> Optional[List[float]]:
    """A ``CXyz`` as ``[x, y, z]``, or None when it is not all three."""
    if element is None:
        return None
    values = [_number(element.findtext(axis), float) for axis in ("x", "y", "z")]
    if any(v is None for v in values):
        return None
    return values


def read_receipt(job) -> Optional[dict]:
    """What one receipt job says about itself, or None if it cannot be read.

    ``events`` are in list order, and each carries its ``position`` in that
    list. The position, not the event number, is what a scene reference needs
    (``EVENTS[2].POSE``); the two differ whenever PanDDA's event numbering
    does not start at zero or skips, so both are kept and neither is derived
    from the other.
    """
    body = receipt_body(job)
    if body is None:
        return None
    inputs = body.find("inputData")
    if inputs is None:
        return None
    text = lambda parent, tag: (parent.findtext(tag) or "").strip()

    events: List[dict] = []
    outputs = body.find("outputData")
    container = outputs.find("EVENTS") if outputs is not None else None
    for position, item in enumerate(list(container) if container is not None else []):
        event = {"position": position}
        for tag in _EVENT_INTS:
            event[tag.lower()] = _number(item.findtext(tag), int)
        for tag in _EVENT_FLOATS:
            event[tag.lower()] = _number(item.findtext(tag), float)
        event["ligand_id"] = text(item, "LIGAND_ID") or None
        event["centroid"] = _centroid(item.find("CENTROID"))
        event["has_map"] = _is_set(item.find("EVENT_MAP"))
        event["has_pose"] = _is_set(item.find("POSE"))
        event["has_scene"] = _is_set(item.find("SCENE"))
        events.append(event)

    return {
        "run_job_uuid": text(inputs, "RUN_JOB_UUID") or None,
        "dtag": text(inputs, "DTAG"),
        "tree": text(inputs, "PANDDA_OUT_DIR") or None,
        "has_apo": _is_set(outputs.find("XYZIN_APO")) if outputs is not None else False,
        "has_dictionary": _is_set(outputs.find("DICT")) if outputs is not None else False,
        "events": events,
    }


# ---------------------------------------------------------------------------
# Reading every receipt of a campaign
# ---------------------------------------------------------------------------

def _member_projects(group) -> Dict[int, models.Project]:
    memberships = group.memberships.filter(
        type=models.ProjectGroupMembership.MembershipType.MEMBER
    ).select_related("project")
    return {m.project_id: m.project for m in memberships}


def _uuid_key(value: Optional[str]) -> str:
    """A uuid as a comparable key: a uuid reaches parameters in both the
    hyphenated and the bare spelling, and neither is canonical there."""
    return (value or "").replace("-", "").lower()


def read_receipts(group, run_job_uuid: Optional[str] = None,
                  tree: Optional[str] = None) -> List[dict]:
    """Every readable receipt in this campaign, newest job last.

    Only current members are read. A receipt outlives the membership it was
    created under, so a dataset dropped from the campaign would otherwise go
    on appearing in the campaign's site views -- the same rule
    ``build_site_scene`` applies to verdicts, for the same reason.

    ``run_job_uuid`` selects one run. A run fanned out from a tree produced
    elsewhere has no run job and so no uuid; ``tree`` selects that one, which
    is the same key fan-out falls back to for idempotence
    (``existing_receipt``). Passing neither reads every run.
    """
    projects = _member_projects(group)
    if not projects:
        return []
    wanted = _uuid_key(run_job_uuid) if run_job_uuid else None

    out: List[dict] = []
    jobs = (models.Job.objects
            .filter(project_id__in=projects, task_name=RECEIPT_TASK)
            .exclude(status=models.Job.Status.TO_DELETE)
            .order_by("id"))
    for job in jobs:
        record = read_receipt(job)
        if record is None or not record["dtag"]:
            continue
        own = _uuid_key(record["run_job_uuid"])
        if wanted is not None and own != wanted:
            continue
        if wanted is None and tree is not None and (record["tree"] or "") != tree:
            continue
        if wanted is None and tree is not None and own:
            # A tree-keyed selection means "the run with no job of its own";
            # a receipt naming a run job belongs to that run, not this one.
            continue
        record["project"] = projects[job.project_id]
        record["job"] = job
        out.append(record)
    return out


def list_runs(group) -> List[dict]:
    """The PanDDA runs this campaign holds receipts from, newest first.

    Runs are enumerated from the receipts rather than from campaign job
    records, because a receipt is the thing that proves a run reached this
    campaign: fan-out may have been driven by the management command over a
    tree produced on another machine, in which case no job here ran it.
    """
    runs: Dict[str, dict] = {}
    for record in read_receipts(group):
        key = record["run_job_uuid"] or f"tree:{record['tree'] or ''}"
        run = runs.setdefault(key, {
            "run_job_uuid": record["run_job_uuid"],
            "tree": record["tree"],
            "n_receipts": 0,
            "n_events": 0,
            "latest_receipt_id": 0,
            "title": None,
            "project": None,
            "number": None,
        })
        run["n_receipts"] += 1
        run["n_events"] += len(record["events"])
        run["latest_receipt_id"] = max(run["latest_receipt_id"], record["job"].id)

    for run in runs.values():
        uuid = run["run_job_uuid"]
        job = models.Job.objects.filter(uuid=uuid).select_related("project").first() if uuid else None
        if job is not None:
            run["title"] = job.title
            run["number"] = job.number
            run["project"] = job.project.name
    return sorted(runs.values(), key=lambda r: r["latest_receipt_id"], reverse=True)


# ---------------------------------------------------------------------------
# The index itself
# ---------------------------------------------------------------------------

def group_events_by_site(receipts: List[dict], centroids: Dict[int, List[float]]) -> dict:
    """Group every event of every receipt by its site number.

    Pure: ``receipts`` are ``read_receipts`` records with ``project`` and
    ``job`` reduced to what a caller wants in the output, and ``centroids``
    maps a site number to the run's own centroid for it. Split out from
    :func:`build_site_index` so the grouping can be tested without a
    database.

    An event PanDDA gave no site number is not silently dropped: it goes to
    ``unsited``, because a run whose sites table was never written puts every
    event there and the caller must be able to say so rather than show an
    empty campaign.
    """
    sites: Dict[int, dict] = {}
    unsited: List[dict] = []

    for record in receipts:
        for event in record["events"]:
            member = dict(event)
            member["dtag"] = record["dtag"]
            member["project"] = record["project"]
            member["receipt"] = record["receipt"]
            idx = event.get("site_idx")
            if idx is None:
                unsited.append(member)
                continue
            site = sites.setdefault(idx, {
                "site_idx": idx,
                "table_centroid": centroids.get(idx),
                "members": [],
            })
            site["members"].append(member)

    ordered = []
    for idx in sorted(sites):
        site = sites[idx]
        # Best first: the reason to open a site is its strongest event, and
        # the panel is read from the top.
        site["members"].sort(
            key=lambda m: (m.get("score") is None, -(m.get("score") or 0.0), m["dtag"]))
        site.update(_site_rollup(site["members"]))
        site["centroid"] = _site_centroid(site)
        ordered.append(site)

    # A site the run listed but no receipt reached is still a site of the
    # run; saying "0 events" is information, and hiding it is not.
    for idx, centroid in sorted(centroids.items()):
        if idx not in sites:
            empty = {"site_idx": idx, "table_centroid": centroid, "members": []}
            empty.update(_site_rollup([]))
            empty["centroid"] = _site_centroid(empty)
            ordered.append(empty)
    ordered.sort(key=lambda s: s["site_idx"])

    return {"sites": ordered, "unsited": unsited}


def _site_centroid(site: dict) -> Optional[List[float]]:
    """Where the site is, preferring the events to the run's own sites table.

    Design note section 9: derive a site centroid from its member events and
    not from the ``pandda_analyse_sites.csv`` column, which is frequently
    ``(0, 0, 0)``. The table value is kept as ``table_centroid`` so a caller
    can see the disagreement, but it is not what anything navigates by.

    The mean is over centroids each stated in its own dataset's frame, which
    is only meaningful because a campaign's members are near-isomorphous --
    good enough to say which pocket a site is, never good enough to fit on.
    A scene computes its own centre from the fits it made
    (``campaign_scene.build_run_site_scene``) and does not use this.
    """
    points = [m["centroid"] for m in site["members"] if m.get("centroid")]
    if points:
        return [sum(axis) / len(points) for axis in zip(*points)]
    table = site.get("table_centroid")
    if table and any(abs(v) > 1e-6 for v in table):
        return list(table)
    return None


def _site_rollup(members: List[dict]) -> dict:
    scores = [m["score"] for m in members if m.get("score") is not None]
    probabilities = [m["hit_probability"] for m in members
                     if m.get("hit_probability") is not None]
    return {
        "n_events": len(members),
        "n_datasets": len({m["dtag"] for m in members}),
        "n_poses": sum(1 for m in members if m.get("has_pose")),
        "best_score": max(scores) if scores else None,
        "best_hit_probability": max(probabilities) if probabilities else None,
    }


def _run_centroids(tree: Optional[str]) -> Dict[int, List[float]]:
    """The run's own site centroids, keyed by site number; ``{}`` if the tree
    has gone. The tree is consulted for this and nothing else, so the index
    degrades to "sites without positions" rather than failing."""
    if not tree or not Path(tree).is_dir():
        return {}
    try:
        from ccp4i2.wrappers.pandda_campaign.script.pandda_run_summary import read_sites
    except ImportError:
        # The run summary reads YAML for other purposes; the sites table is
        # CSV. Losing centroids is a worse index, not a broken one.
        return {}
    out: Dict[int, List[float]] = {}
    for site in read_sites(tree):
        centroid = site.get("centroid")
        out[site["site_idx"]] = list(centroid) if centroid else None
    return out


def build_site_index(group, run_job_uuid: Optional[str] = None) -> dict:
    """The site index of one PanDDA run over one campaign.

    Without ``run_job_uuid`` the most recent run is indexed, "most recent"
    being the run whose newest receipt is newest -- the receipts are what
    this campaign has, and a run nobody fanned out is not navigable.
    """
    runs = list_runs(group)
    if run_job_uuid is None:
        if not runs:
            return {"run": None, "runs": [], "tree": None, "trees_disagree": None,
                    "sites": [], "unsited": [], "stats": _stats([], [])}
        chosen = runs[0]
    else:
        wanted = _uuid_key(run_job_uuid)
        chosen = next((r for r in runs if _uuid_key(r["run_job_uuid"]) == wanted), None)

    # Select by uuid where the run has one, by tree where it does not. Asking
    # for "the latest run" and getting every run's events merged under one set
    # of site numbers would be a silent wrong answer, and PanDDA's numbering
    # makes it a plausible-looking one.
    if chosen is None:
        receipts = []
    elif chosen["run_job_uuid"]:
        receipts = read_receipts(group, run_job_uuid=chosen["run_job_uuid"])
    else:
        receipts = read_receipts(group, tree=chosen["tree"])

    trees = {r["tree"] for r in receipts if r["tree"]}
    tree = (chosen or {}).get("tree") or (sorted(trees)[0] if trees else None)

    reduced = [{
        "dtag": r["dtag"],
        "events": r["events"],
        "project": {"id": r["project"].id, "uuid": str(r["project"].uuid),
                    "name": r["project"].name},
        "receipt": {"uuid": str(r["job"].uuid), "number": r["job"].number,
                    "id": r["job"].id, "status": r["job"].status},
    } for r in receipts]

    grouped = group_events_by_site(reduced, _run_centroids(tree))
    return {
        "run": chosen,
        "runs": runs,
        "tree": tree,
        "trees_disagree": sorted(trees) if len(trees) > 1 else None,
        "sites": grouped["sites"],
        "unsited": grouped["unsited"],
        "stats": _stats(reduced, grouped["sites"]),
    }


def _stats(receipts: List[dict], sites: List[dict]) -> dict:
    return {
        "n_receipts": len(receipts),
        "n_sites": len(sites),
        "n_events": sum(len(r["events"]) for r in receipts),
        "n_poses": sum(1 for r in receipts for e in r["events"] if e.get("has_pose")),
    }
