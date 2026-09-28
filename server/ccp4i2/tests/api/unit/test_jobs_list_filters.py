"""``GET /api/ccp4i2/jobs/`` filters: ``project`` (already contracted) and the
additive ``task_name``, ``status`` and ``parent`` filters. Without the latter
a caller asking for "the pandda_campaign jobs" received the whole table."""
import pytest
from rest_framework.test import APIClient

from ccp4i2.db import models


@pytest.fixture
def client(bypass_api_permissions):
    return APIClient()


@pytest.fixture
def jobs(test_project_path):
    test_project_path.mkdir(exist_ok=True)
    p1 = models.Project.objects.create(name="filters_p1", directory=str(test_project_path / "p1"))
    p2 = models.Project.objects.create(name="filters_p2", directory=str(test_project_path / "p2"))
    top = models.Job.objects.create(project=p1, number="1", task_name="SubstituteLigand",
                                    title="t", status=models.Job.Status.FINISHED)
    sub = models.Job.objects.create(project=p1, number="1.3", task_name="i2Dimple", title="d",
                                    status=models.Job.Status.FINISHED, parent=top)
    pending = models.Job.objects.create(project=p1, number="2", task_name="pandda_campaign",
                                        title="p", status=models.Job.Status.PENDING)
    other = models.Job.objects.create(project=p2, number="1", task_name="pandda_campaign",
                                      title="p", status=models.Job.Status.FINISHED)
    return {"p1": p1, "p2": p2, "top": top, "sub": sub, "pending": pending, "other": other}


def _ids(resp):
    assert resp.status_code == 200, resp.content
    body = resp.json()
    rows = body.get("results", body) if isinstance(body, dict) else body
    return {row["id"] for row in rows}


def test_task_name_filters_across_projects(client, jobs):
    got = _ids(client.get("/api/ccp4i2/jobs/", {"task_name": "pandda_campaign"}))
    assert got == {jobs["pending"].id, jobs["other"].id}


def test_filters_are_additive(client, jobs):
    got = _ids(client.get("/api/ccp4i2/jobs/", {"task_name": "pandda_campaign", "project": jobs["p1"].id}))
    assert got == {jobs["pending"].id}
    got = _ids(client.get("/api/ccp4i2/jobs/", {"task_name": "pandda_campaign",
                                                 "status": models.Job.Status.FINISHED}))
    assert got == {jobs["other"].id}


def test_parent_null_selects_top_level_jobs(client, jobs):
    got = _ids(client.get("/api/ccp4i2/jobs/", {"project": jobs["p1"].id, "parent": "null"}))
    assert got == {jobs["top"].id, jobs["pending"].id}
    got = _ids(client.get("/api/ccp4i2/jobs/", {"parent": jobs["top"].id}))
    assert got == {jobs["sub"].id}


def test_unparseable_values_match_nothing_rather_than_everything(client, jobs):
    assert _ids(client.get("/api/ccp4i2/jobs/", {"status": "finished"})) == set()
    assert _ids(client.get("/api/ccp4i2/jobs/", {"parent": "abc"})) == set()
