"""API tests for the bibliography endpoints.

Builds minimal Project/Job rows directly via the ORM (no external test data),
then exercises:
  - GET /api/ccp4i2/jobs/{id}/bibliography/
  - GET /api/ccp4i2/projects/{id}/bibliography/
verifying subjob burrowing (child task refs surface) and the response envelope.
"""

import pytest
from rest_framework.test import APIClient

from ccp4i2.db import models


@pytest.fixture
def client(bypass_api_permissions):
    return APIClient()


@pytest.fixture
def project(test_project_path):
    test_project_path.mkdir(exist_ok=True)
    return models.Project.objects.create(
        name="biblio_test", directory=str(test_project_path)
    )


def _make_job(project, number, task_name, parent=None):
    return models.Job.objects.create(
        project=project,
        number=number,
        title=f"{task_name} job",
        task_name=task_name,
        status=models.Job.Status.FINISHED,
        parent=parent,
    )


def test_job_bibliography_includes_subjob_refs(client, project):
    # A prosmart_refmac pipeline with a refmac child: the pipeline itself is
    # non-citable, but subjob burrowing must surface refmac's references.
    parent = _make_job(project, "1", "prosmart_refmac")
    _make_job(project, "1.1", "refmac", parent=parent)

    resp = client.get(f"/api/ccp4i2/jobs/{parent.id}/bibliography/")
    assert resp.status_code == 200
    body = resp.json()
    assert body["success"] is True
    refs = body["data"]["result"]
    titles = " ".join((r.get("title") or "") for r in refs)
    assert "REFMAC5" in titles, "refmac subjob references should surface"
    # Every ref carries the expected shape.
    assert all({"title", "authors", "source", "link", "cited_by"} <= set(r) for r in refs)


def test_job_bibliography_monolithic_pipeline(client, project):
    # crank2 drives programs internally (no subjobs) -> cites-map must supply them.
    job = _make_job(project, "1", "crank2")
    resp = client.get(f"/api/ccp4i2/jobs/{job.id}/bibliography/")
    assert resp.status_code == 200
    refs = resp.json()["data"]["result"]
    cited = {r["cited_by"] for r in refs}
    assert "buccaneer" in cited and "shelxd" in cited


def test_project_bibliography_union(client, project):
    _make_job(project, "1", "refmac")
    _make_job(project, "2", "pointless")
    resp = client.get(f"/api/ccp4i2/projects/{project.id}/bibliography/")
    assert resp.status_code == 200
    refs = resp.json()["data"]["result"]
    cited = {r["cited_by"] for r in refs}
    assert "refmac" in cited and "pointless" in cited


def test_job_bibliography_missing_job_404(client):
    resp = client.get("/api/ccp4i2/jobs/999999/bibliography/")
    assert resp.status_code == 404
