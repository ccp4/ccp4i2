"""E2E tests for the PanDDA run-site axis.

    GET /api/ccp4i2/projectgroups/{id}/pandda-sites/
    GET /api/ccp4i2/projectgroups/{id}/pandda-sites/{site_idx}/scene/

This is the axis PanDDA itself found, as against the curated ``CampaignSite``
axis that ``test_site_scene_api`` covers. The distinction the tests exist to
hold: **membership here needs no verdict**. The same synthetic campaign is
used, with no ``SiteEvaluation`` rows anywhere, and a scene of two ligands at
one site still comes back -- which is the whole reason for the second axis,
because a run has found its sites before anybody has judged anything at them.

Each member gets a ``pandda_events`` receipt whose parameter file carries its
events and whose ``XYZIN_APO`` is a real file, because the fit that brings
every dataset into one frame is computed from the apo model and not from the
three dozen atoms of the pose it draws.

Uses the api/ conftest, which auto-applies django_db(transaction=True) and
sets AllowAny on the viewsets -- do NOT add @pytest.mark.django_db here.
"""

import uuid

from rest_framework.test import APIClient

from ccp4i2.db import models

from .test_summary_scene_api import LIG_SHIFT, _build_campaign, _write_pdb

RUN_UUID = "11111111-2222-3333-4444-555555555555"

# The site, in the parent's frame: inside the helix _write_pdb writes, so a
# pocket really is found there (the same point test_site_scene_api uses).
SITE_POSITION = (5.0, 5.0, 5.0)


def _receipt(project, dtag, events, run_uuid=RUN_UUID, shift=(0.0, 0.0, 0.0),
             number="2", apo=True):
    """A finished receipt: its parameter file, and an apo model on disk."""
    job = models.Job.objects.create(
        uuid=uuid.uuid4(), project=project, number=number,
        title=f"PanDDA events for {dtag}", task_name="pandda_events",
        status=models.Job.Status.FINISHED,
    )
    job.directory.mkdir(parents=True, exist_ok=True)

    rows = []
    for event in events:
        x, y, z = event["centroid"]
        rows.append(
            "<CPanddaEvent>"
            f"<EVENT_IDX>{event['idx']}</EVENT_IDX>"
            f"<SITE_IDX>{event['site']}</SITE_IDX>"
            f"<SCORE>{event['score']}</SCORE>"
            f"<HIT_PROBABILITY>{event.get('probability', 0.5)}</HIT_PROBABILITY>"
            f"<CENTROID><x>{x}</x><y>{y}</y><z>{z}</z></CENTROID>"
            + ("<POSE><baseName>pose.pdb</baseName></POSE>"
               if event.get("pose", True) else "<POSE></POSE>")
            + "<EVENT_MAP><baseName>event.ccp4</baseName></EVENT_MAP>"
            "</CPanddaEvent>"
        )
    apo_element = "<baseName>apo.pdb</baseName>" if apo else ""
    (job.directory / "params.xml").write_text(
        "<ccp4i2><ccp4i2_body id='pandda_events'>"
        "<inputData>"
        "<PANDDA_OUT_DIR>/runs/pandda2_out</PANDDA_OUT_DIR>"
        f"<DTAG>{dtag}</DTAG>"
        f"<RUN_JOB_UUID>{run_uuid}</RUN_JOB_UUID>"
        "</inputData>"
        f"<outputData><XYZIN_APO>{apo_element}</XYZIN_APO>"
        f"<EVENTS>{''.join(rows)}</EVENTS></outputData>"
        "</ccp4i2_body></ccp4i2>"
    )

    if apo:
        _write_pdb(job.directory / "apo.pdb", shift=shift)
        pdb_type, _ = models.FileType.objects.get_or_create(name="chemical/x-pdb")
        models.File.objects.create(
            uuid=uuid.uuid4(), name="apo.pdb",
            directory=models.File.Directory.JOB_DIR, type=pdb_type,
            job=job, job_param_name="XYZIN_APO",
        )
    return job


def _project(name):
    return models.Project.objects.get(name=name)


def _shifted(point, shift):
    return tuple(p + s for p, s in zip(point, shift))


def _campaign_with_a_run(root):
    """The synthetic campaign, plus a run that found two sites over it.

    frag_lig's model is deposited in its own frame (``LIG_SHIFT``), so its
    event centroid is stated in that frame too -- which is what the fit has
    to undo before the two ligands can be said to be at the same site.
    """
    root.mkdir(parents=True, exist_ok=True)
    group = _build_campaign(root)
    _receipt(_project("frag_drg"), "x0001", [
        {"idx": 1, "site": 1, "score": 0.91, "centroid": SITE_POSITION},
        {"idx": 2, "site": 2, "score": 0.30, "centroid": (20.0, 20.0, 20.0)},
    ])
    _receipt(_project("frag_lig"), "x0002", [
        {"idx": 1, "site": 1, "score": 0.44,
         "centroid": _shifted(SITE_POSITION, LIG_SHIFT)},
    ], shift=LIG_SHIFT)
    # An event PanDDA found no pose for: in the index, out of the scene.
    _receipt(_project("frag_apo"), "x0003", [
        {"idx": 1, "site": 1, "score": 0.10, "centroid": SITE_POSITION,
         "pose": False},
    ])
    return group


def _index(group, query=""):
    return APIClient().get(
        f"/api/ccp4i2/projectgroups/{group.id}/pandda-sites/{query}")


def _scene(group, site_idx, query=""):
    return APIClient().get(
        f"/api/ccp4i2/projectgroups/{group.id}/pandda-sites/{site_idx}/scene/{query}")


# --------------------------------------------------------------------------
# The index
# --------------------------------------------------------------------------

def test_index_groups_the_runs_events_by_its_own_sites(
        bypass_api_permissions, test_project_path):
    group = _campaign_with_a_run(test_project_path)

    response = _index(group)
    assert response.status_code == 200, response.content
    data = response.json()

    assert len(data["runs"]) == 1
    assert data["runs"][0]["run_job_uuid"] == RUN_UUID
    assert data["runs"][0]["n_receipts"] == 3
    assert data["stats"]["n_events"] == 4
    assert data["stats"]["n_poses"] == 3

    sites = {s["site_idx"]: s for s in data["sites"]}
    assert sorted(sites) == [1, 2]
    assert sites[1]["n_events"] == 3 and sites[1]["n_datasets"] == 3
    assert sites[1]["n_poses"] == 2          # frag_apo's event has no build
    assert sites[1]["best_score"] == 0.91
    # Best first, so the panel reads top-down in order of what is worth a look.
    assert [m["dtag"] for m in sites[1]["members"]] == ["x0001", "x0002", "x0003"]
    assert sites[2]["n_events"] == 1

    member = sites[1]["members"][0]
    assert member["project"]["name"] == "frag_drg"
    assert member["receipt"]["number"] == "2"
    # The position, not the event number: it is what a scene reference needs.
    assert member["position"] == 0


def test_index_needs_no_verdicts(bypass_api_permissions, test_project_path):
    """The point of the second axis. The curated site scene draws hits and
    this campaign has no SiteEvaluation rows at all."""
    group = _campaign_with_a_run(test_project_path)
    assert models.SiteEvaluation.objects.count() == 0
    assert _index(group).json()["sites"][0]["n_events"] == 3


def test_index_of_a_campaign_with_no_receipts_is_empty_not_an_error(
        bypass_api_permissions, test_project_path):
    test_project_path.mkdir(parents=True, exist_ok=True)
    group = _build_campaign(test_project_path)
    data = _index(group).json()
    assert data["run"] is None and data["runs"] == [] and data["sites"] == []


def test_a_dropped_member_leaves_the_index(bypass_api_permissions, test_project_path):
    """A receipt outlives the membership it was created under, so the index
    must read membership and not merely the receipt's project."""
    group = _campaign_with_a_run(test_project_path)
    models.ProjectGroupMembership.objects.filter(
        group=group, project=_project("frag_lig")).delete()

    site = _index(group).json()["sites"][0]
    assert [m["dtag"] for m in site["members"]] == ["x0001", "x0003"]


def test_a_second_run_does_not_merge_into_the_first(
        bypass_api_permissions, test_project_path):
    """PanDDA renumbers sites every run, so merging two runs under one set of
    site numbers would be a wrong answer that looks right."""
    group = _campaign_with_a_run(test_project_path)
    other = "99999999-8888-7777-6666-555555555555"
    _receipt(_project("frag_drg"), "x0001",
             [{"idx": 1, "site": 1, "score": 0.99, "centroid": SITE_POSITION}],
             run_uuid=other, number="3")

    data = _index(group).json()
    assert len(data["runs"]) == 2
    # Default is the run whose newest receipt is newest: the one just added.
    assert data["run"]["run_job_uuid"] == other
    assert data["stats"]["n_events"] == 1

    asked = _index(group, f"?run={RUN_UUID}").json()
    assert asked["run"]["run_job_uuid"] == RUN_UUID
    assert asked["stats"]["n_events"] == 4


# --------------------------------------------------------------------------
# The scene
# --------------------------------------------------------------------------

def test_site_scene_draws_every_pose_at_the_site_in_one_frame(
        bypass_api_permissions, test_project_path):
    group = _campaign_with_a_run(test_project_path)

    response = _scene(group, 1)
    assert response.status_code == 200, response.content
    data = response.json()
    scene, stats = data["scene"], data["stats"]

    assert stats["members_claimed"] == 3
    assert stats["members_drawn"] == 2       # the pose-less event is skipped
    assert [s["reason"] for s in stats["skipped"]] == ["no autobuilt pose for this event"]
    assert stats["parent_present"] is True

    # Poses are referenced by job number and parameter, the way a receipt's
    # own scenes reference their outputs: it survives a project move.
    poses = [f for f in scene["files"] if f.get("param", "").endswith(".POSE")]
    assert {f["param"] for f in poses} == {"EVENTS[0].POSE"}
    assert all(f["job"] == "2" and f["kind"] == "coordinates" for f in poses)
    assert len({f["projectId"] for f in poses}) == 2

    # Fitted onto the exemplar, and the fit is what puts them at one site.
    assert len(scene["superpose"]) == 2
    assert all(entry["method"] == "matrix" for entry in scene["superpose"])

    # The site's position is learnt from the members, in the parent's frame.
    centre = stats["centre"]
    assert centre is not None
    assert all(abs(c - p) < 1.0 for c, p in zip(centre, SITE_POSITION)), centre
    # Moorhen's view origin is the negation of the point to centre on.
    assert scene["view"]["origin"] == [-centre[0], -centre[1], -centre[2]]
    assert stats["pocket_residues"] > 0


def test_site_scene_carries_no_maps(bypass_api_permissions, test_project_path):
    """A pose can be moved into the exemplar's frame by a matrix and an event
    map cannot, so twenty datasets' maps would be nineteen maps in the wrong
    place. The map belongs in the receipt's own event scene."""
    group = _campaign_with_a_run(test_project_path)
    scene = _scene(group, 1).json()["scene"]
    assert scene.get("maps", []) == []
    assert not [f for f in scene["files"] if f.get("kind") == "map"]


def test_site_scene_without_superposition_draws_frames_as_they_are(
        bypass_api_permissions, test_project_path):
    group = _campaign_with_a_run(test_project_path)
    data = _scene(group, 1, "?superpose=none").json()
    assert "superpose" not in data["scene"]
    assert data["stats"]["members_drawn"] == 2


def test_a_site_the_run_does_not_have_is_404(
        bypass_api_permissions, test_project_path):
    group = _campaign_with_a_run(test_project_path)
    assert _scene(group, 97).status_code == 404
