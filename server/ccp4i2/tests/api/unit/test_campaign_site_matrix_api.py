"""The campaign's dataset x site matrix, end to end.

    GET   /api/ccp4i2/projectgroups/{id}/member_projects/   site_cells, frame_mismatch,
                                                            current_model_job
    PATCH /api/ccp4i2/projectgroups/{id}/sites/{site}/      radius

and the ``CampaignEvent`` projection under it: recorded when a receipt
reaches a recorded status, replaced when re-recorded, dropped when the
receipt leaves that status, and rebuilt by ``backfill_campaign_events``.

Built on the synthetic campaign and receipts of test_pandda_site_index_api:
three members with a refmac model each, and one run whose receipts put events
at the parent-frame point SITE_POSITION (frag_lig's 2.3 A off it, in its own
frame) and one more of frag_drg's at (20, 20, 20).

Uses the api/ conftest, which auto-applies django_db(transaction=True) and
sets AllowAny on the viewsets -- do NOT add @pytest.mark.django_db here.
"""
import io
import uuid

from django.core.management import call_command
from django.db import connection
from django.test.utils import CaptureQueriesContext
from rest_framework.test import APIClient

from ccp4i2.db import models
from ccp4i2.lib import campaign_events

from .test_pandda_site_index_api import (
    RUN_UUID,
    SITE_POSITION,
    _campaign_with_a_run,
    _project,
    _receipt,
)
from .test_summary_scene_api import _make_refined_project


def _member_rows(group):
    response = APIClient().get(
        f"/api/ccp4i2/projectgroups/{group.id}/member_projects/")
    assert response.status_code == 200, response.content
    return {row["name"]: row for row in response.json()}


def _site(group, name, position, order, **kw):
    """A site saved as the viewer saves one: ``origin`` is Moorhen's view
    origin, the NEGATED position."""
    x, y, z = position
    return models.CampaignSite.objects.create(
        group=group, name=name, origin_x=-x, origin_y=-y, origin_z=-z,
        order=order, **kw)


def _sites(group):
    """S1 on the run's site 1, S2 tight round frag_drg's second event, S3 far
    from every event."""
    s1 = _site(group, "S1", SITE_POSITION, 0)
    s2 = _site(group, "S2", (21.0, 20.0, 20.0), 1, radius=3.0)
    s3 = _site(group, "S3", (40.0, 40.0, 40.0), 2)
    return s1, s2, s3


def _receipt_job(project_name):
    return models.Job.objects.get(project=_project(project_name),
                                  task_name="pandda_events")


# --------------------------------------------------------------------------
# Population
# --------------------------------------------------------------------------

def test_a_receipt_reaching_finished_records_its_events(
        bypass_api_permissions, test_project_path):
    group = _campaign_with_a_run(test_project_path)
    receipt = _receipt_job("frag_drg")
    # The fixture creates the job FINISHED before writing its params.xml, so
    # the creation signal found nothing to record -- as for any job whose
    # record is not yet on disk.
    assert not models.CampaignEvent.objects.filter(receipt=receipt).exists()

    receipt.status = models.Job.Status.RUNNING
    receipt.save()
    receipt.status = models.Job.Status.FINISHED
    receipt.save()

    rows = {e.event_idx: e for e in models.CampaignEvent.objects.filter(receipt=receipt)}
    assert sorted(rows) == [1, 2]
    first = rows[1]
    assert first.project == _project("frag_drg")
    assert first.group == group                # the project's only campaign
    assert first.run_job_uuid == uuid.UUID(RUN_UUID)
    assert first.dtag == "x0001"
    assert first.site_idx == 1
    assert first.centroid == list(SITE_POSITION)
    assert first.hit_probability == 0.5
    assert first.has_pose and first.has_map
    # The dataset's cell, from the receipt's apo model.
    assert first.cell == [60.0, 60.0, 60.0, 90.0, 90.0, 90.0]


def test_an_unsatisfactory_receipt_is_recorded_too(
        bypass_api_permissions, test_project_path):
    """A short receipt is UNSATISFACTORY and its events are still real."""
    _campaign_with_a_run(test_project_path)
    receipt = _receipt_job("frag_lig")
    receipt.status = models.Job.Status.UNSATISFACTORY
    receipt.save()
    assert models.CampaignEvent.objects.filter(receipt=receipt).count() == 1


def test_re_recording_replaces_rows_and_leaving_the_status_drops_them(
        bypass_api_permissions, test_project_path):
    _campaign_with_a_run(test_project_path)
    receipt = _receipt_job("frag_drg")
    assert campaign_events.record_receipt(receipt) == 2
    assert campaign_events.record_receipt(receipt) == 2
    assert models.CampaignEvent.objects.filter(receipt=receipt).count() == 2

    # The receipt rerun with a different result: its rows are the new result.
    params = receipt.directory / "params.xml"
    text = params.read_text()
    cut = text.index("<CPanddaEvent>", text.index("<CPanddaEvent>") + 1)
    params.write_text(text[:cut] + text[text.index("</EVENTS>"):])
    assert campaign_events.record_receipt(receipt) == 1

    receipt.status = models.Job.Status.TO_DELETE
    receipt.save()
    assert not models.CampaignEvent.objects.filter(receipt=receipt).exists()


def test_backfill_command_is_idempotent_and_scoped_to_a_campaign(
        bypass_api_permissions, test_project_path):
    group = _campaign_with_a_run(test_project_path)
    # A receipt outside the campaign: a project that is nobody's member.
    stray = models.Project.objects.create(
        name="stray", directory=str(test_project_path / "stray"))
    _receipt(stray, "x0099", [{"idx": 1, "site": 1, "score": 0.5,
                               "centroid": SITE_POSITION}])

    for _ in range(2):
        out = io.StringIO()
        call_command("backfill_campaign_events", "--group", str(group.uuid), stdout=out)
        assert "[OK] 3 receipt(s)" in out.getvalue()
        assert models.CampaignEvent.objects.count() == 4
    assert not models.CampaignEvent.objects.filter(project=stray).exists()

    call_command("backfill_campaign_events", stdout=io.StringIO())
    assert models.CampaignEvent.objects.count() == 5


# --------------------------------------------------------------------------
# The payload
# --------------------------------------------------------------------------

def test_member_projects_carries_a_cell_for_every_site(
        bypass_api_permissions, test_project_path):
    group = _campaign_with_a_run(test_project_path)
    s1, s2, s3 = _sites(group)
    models.SiteEvaluation.objects.create(
        project=_project("frag_drg"), site=s1, verdict="hit")
    # A verdict where PanDDA found nothing: shown as an outlined box.
    models.SiteEvaluation.objects.create(
        project=_project("frag_apo"), site=s3, verdict="empty")
    campaign_events.backfill(group)

    rows = _member_rows(group)
    drg = rows["frag_drg"]
    assert set(drg["site_cells"]) == {str(s1.uuid), str(s2.uuid), str(s3.uuid)}
    assert drg["site_cells"][str(s1.uuid)] == {
        "event": {"event_idx": 1, "hit_probability": 0.5, "distance": 0.0,
                  "has_pose": True,
                  # Which receipt: a click opens that receipt's evidence.
                  "receipt_job_id": _receipt_job("frag_drg").id},
        "verdict": "hit",
    }
    # S2's radius is 3 A and the event is 1 A from its origin.
    assert drg["site_cells"][str(s2.uuid)]["event"]["event_idx"] == 2
    assert drg["site_cells"][str(s2.uuid)]["event"]["distance"] == 1.0
    assert drg["site_cells"][str(s3.uuid)] == {"event": None, "verdict": None}
    assert drg["frame_mismatch"] is None

    # Compared directly: frag_lig's centroid is in its own frame, 2.29 A off.
    lig = rows["frag_lig"]["site_cells"][str(s1.uuid)]
    assert lig["event"]["distance"] == 2.29 and lig["verdict"] is None

    apo = rows["frag_apo"]["site_cells"]
    assert apo[str(s1.uuid)]["event"]["has_pose"] is False
    assert apo[str(s3.uuid)] == {"event": None, "verdict": "empty"}

    # The model a click opens: each member's refmac job.
    model = drg["current_model_job"]
    assert model["task_name"] == "refmac" and model["number"] == "1"
    job = models.Job.objects.get(id=model["id"])
    assert model["uuid"] == str(job.uuid)

    # The older fields are unchanged.
    assert drg["site_evaluations"] == [
        {"site_id": s1.id, "site_name": "S1", "verdict": "hit"}]
    assert drg["sites_evaluated"] == 1 and drg["sites_total"] == 3


def test_a_dataset_with_no_receipt_has_empty_cells(
        bypass_api_permissions, test_project_path):
    group = _campaign_with_a_run(test_project_path)
    s1, _, _ = _sites(group)
    # Nothing recorded yet: the cells are there, the events are not.
    row = _member_rows(group)["frag_drg"]
    assert row["site_cells"][str(s1.uuid)] == {"event": None, "verdict": None}
    assert row["frame_mismatch"] is None


def test_a_newer_run_supersedes_the_older_one(
        bypass_api_permissions, test_project_path):
    group = _campaign_with_a_run(test_project_path)
    s1, s2, _ = _sites(group)
    newer = _receipt(_project("frag_drg"), "x0001",
             [{"idx": 7, "site": 3, "score": 0.9, "probability": 0.95,
               "centroid": (5.5, 5.0, 5.0)}],
             run_uuid="99999999-8888-7777-6666-555555555555", number="3")
    campaign_events.backfill(group)

    cells = _member_rows(group)["frag_drg"]["site_cells"]
    assert cells[str(s1.uuid)]["event"] == {
        "event_idx": 7, "hit_probability": 0.95, "distance": 0.5, "has_pose": True,
        "receipt_job_id": newer.id}
    assert cells[str(s2.uuid)]["event"] is None    # the old run's event 2


def test_a_dataset_off_the_parent_cell_matches_nothing_and_says_why(
        bypass_api_permissions, test_project_path):
    group = _campaign_with_a_run(test_project_path)
    s1, _, _ = _sites(group)
    apo = _receipt_job("frag_lig").directory / "apo.pdb"
    apo.write_text(apo.read_text().replace(
        "CRYST1   60.000   60.000   60.000", "CRYST1   64.000   60.000   60.000"))
    campaign_events.backfill(group)

    rows = _member_rows(group)
    assert "cell a" in rows["frag_lig"]["frame_mismatch"]
    assert rows["frag_lig"]["site_cells"][str(s1.uuid)]["event"] is None
    assert rows["frag_drg"]["frame_mismatch"] is None


def test_current_model_prefers_a_refinement_and_a_finished_pipeline(
        bypass_api_permissions, test_project_path):
    group = _campaign_with_a_run(test_project_path)
    project = _project("frag_drg")

    def job(number, task, status=models.Job.Status.FINISHED, parent=None):
        return models.Job.objects.create(
            uuid=uuid.uuid4(), project=project, number=number, title=task,
            task_name=task, status=status, parent=parent)

    job("4", "i2Dimple")                       # newer, but only a DIMPLE fit
    assert _member_rows(group)["frag_drg"]["current_model_job"]["number"] == "1"

    pipeline = job("5", "SubstituteLigand")
    job("5.1", "servalcat_pipe", parent=pipeline)
    assert _member_rows(group)["frag_drg"]["current_model_job"] == {
        "id": pipeline.id, "uuid": str(pipeline.uuid), "number": "5",
        "task_name": "SubstituteLigand"}

    other = _project("frag_apo")
    models.Job.objects.filter(project=other).update(status=models.Job.Status.FAILED)
    assert _member_rows(group)["frag_apo"]["current_model_job"] is None


def test_member_projects_costs_the_same_queries_for_more_members(
        bypass_api_permissions, test_project_path):
    group = _campaign_with_a_run(test_project_path)
    _sites(group)
    campaign_events.backfill(group)
    models.SiteEvaluation.objects.create(
        project=_project("frag_drg"), site=group.site_set.first(), verdict="hit")
    client = APIClient()
    url = f"/api/ccp4i2/projectgroups/{group.id}/member_projects/"

    with CaptureQueriesContext(connection) as before:
        assert client.get(url).status_code == 200

    pdb_type, _ = models.FileType.objects.get_or_create(name="chemical/x-pdb")
    dict_type, _ = models.FileType.objects.get_or_create(
        name="application/refmac-dictionary")
    key, _ = models.JobValueKey.objects.get_or_create(
        name="RFree", defaults={"description": "R-free"})
    for i in range(3):
        project = _make_refined_project(
            test_project_path, f"extra_{i}", pdb_type, dict_type)
        models.ProjectGroupMembership.objects.create(
            group=group, project=project,
            type=models.ProjectGroupMembership.MembershipType.MEMBER)
        models.JobFloatValue.objects.create(
            job=project.jobs.first(), key=key, value=0.25)
        _receipt(project, f"x01{i}", [{"idx": 1, "site": 1, "score": 0.5,
                                       "centroid": SITE_POSITION}])
        models.SiteEvaluation.objects.create(
            project=project, site=group.site_set.first(), verdict="unclear")
    campaign_events.backfill(group)

    with CaptureQueriesContext(connection) as after:
        response = client.get(url)
    assert response.status_code == 200
    assert len(response.json()) == 6
    # group, memberships, tags, sites, verdicts, jobs, float and char KPI
    # values, events, parent membership, parent files -- and the session's
    # own bookkeeping. Twelve on the day this was written.
    assert len(after.captured_queries) <= 12, [
        q["sql"][:90] for q in after.captured_queries]
    assert len(after.captured_queries) == len(before.captured_queries), (
        [q["sql"] for q in after.captured_queries])


# --------------------------------------------------------------------------
# The site radius
# --------------------------------------------------------------------------

def test_site_radius_defaults_is_served_and_is_editable(
        bypass_api_permissions, test_project_path):
    group = _campaign_with_a_run(test_project_path)
    client = APIClient()
    base = f"/api/ccp4i2/projectgroups/{group.id}/sites/"

    created = client.post(base, {"name": "P", "origin": [1, 2, 3]}, format="json")
    assert created.status_code == 201 and created.json()["radius"] == 8.0
    site_id = created.json()["id"]

    wide = client.post(base, {"name": "Q", "origin": [1, 2, 3], "radius": 12},
                       format="json")
    assert wide.json()["radius"] == 12.0

    patched = client.patch(f"{base}{site_id}/", {"radius": 4.5}, format="json")
    assert patched.status_code == 200 and patched.json()["radius"] == 4.5
    assert models.CampaignSite.objects.get(id=site_id).radius == 4.5

    for bad in (0, -1, "wide", None):
        response = client.patch(f"{base}{site_id}/", {"radius": bad}, format="json")
        assert response.status_code == 400, bad
    assert models.CampaignSite.objects.get(id=site_id).radius == 4.5

    listed = {s["name"]: s for s in client.get(base).json()}
    assert listed["P"]["radius"] == 4.5 and listed["Q"]["radius"] == 12.0


# --------------------------------------------------------------------------
# Self-healing: a receipt that expects events but has no rows is recorded by
# the page that needs it (a campaign from before the table, a project
# imported before its files were on disk).
# --------------------------------------------------------------------------

def _expect_events(receipt, count):
    key, _ = models.JobValueKey.objects.get_or_create(
        name=campaign_events.EXPECTED_EVENTS_KPI,
        defaults={"description": "Events expected"})
    models.JobFloatValue.objects.update_or_create(
        job=receipt, key=key, defaults={"value": count})


def test_the_page_records_a_receipt_that_expects_events_and_has_none(
        bypass_api_permissions, test_project_path):
    group = _campaign_with_a_run(test_project_path)
    s1, _, _ = _sites(group)
    drg, lig = _receipt_job("frag_drg"), _receipt_job("frag_lig")
    _expect_events(drg, 2)
    _expect_events(lig, 0)   # claims none: left alone, though it has events
    models.CampaignEvent.objects.all().delete()

    rows = _member_rows(group)

    assert models.CampaignEvent.objects.filter(receipt=drg).count() == 2
    assert rows["frag_drg"]["site_cells"][str(s1.uuid)]["event"]["event_idx"] == 1
    assert not models.CampaignEvent.objects.filter(receipt=lig).exists()
    assert rows["frag_lig"]["site_cells"][str(s1.uuid)]["event"] is None

    # Healed once: the next load finds the rows and records nothing again.
    with CaptureQueriesContext(connection) as again:
        _member_rows(group)
    assert not any("INSERT" in q["sql"] for q in again.captured_queries)


def test_a_site_is_matched_at_its_position_not_its_saved_view_origin(
        bypass_api_permissions, test_project_path):
    """Moorhen's view origin is the negated centre. Matching against the
    stored value measured every event against the site's reflection through
    the molecule origin, and a real campaign showed no events at all."""
    group = _campaign_with_a_run(test_project_path)
    at_site = _site(group, "At the events", SITE_POSITION, 0)
    # A site whose stored origin equals the events' position is at their
    # reflection, so nothing may match it.
    reflected = models.CampaignSite.objects.create(
        group=group, name="Reflection", origin_x=SITE_POSITION[0],
        origin_y=SITE_POSITION[1], origin_z=SITE_POSITION[2], order=1)
    campaign_events.backfill(group)

    drg = _member_rows(group)["frag_drg"]["site_cells"]
    assert drg[str(at_site.uuid)]["event"]["event_idx"] == 1
    assert drg[str(reflected.uuid)]["event"] is None
    assert at_site.position == [float(c) for c in SITE_POSITION]
