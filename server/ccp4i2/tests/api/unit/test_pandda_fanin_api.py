"""The fill-from-campaign methods over the generic object_method endpoint:
no PanDDA-specific route, the plugin's own methods do the work.

  - POST /api/ccp4i2/jobs/{id}/object_method/  {object_path: "pandda_campaign",
        method_name: "campaignCandidates" | "fillDatasetsFromCampaign"}
"""
import pytest
from rest_framework.test import APIClient

from ccp4i2.db import models
from ccp4i2.lib.utils.jobs.create import create_job


@pytest.fixture
def client(bypass_api_permissions):
    return APIClient()


@pytest.fixture
def project(test_project_path):
    test_project_path.mkdir(exist_ok=True)
    return models.Project.objects.create(name="fanin_api", directory=str(test_project_path))


def _call(client, job, method):
    return client.post(f"/api/ccp4i2/jobs/{job.id}/object_method/",
                       {"object_path": "pandda_campaign", "method_name": method, "args": [], "kwargs": {}},
                       format="json")


def test_a_project_with_no_campaign_says_so(client, project):
    job = models.Job.objects.get(uuid=create_job(projectId=str(project.uuid), taskName="pandda_campaign"))
    resp = _call(client, job, "campaignCandidates")
    assert resp.status_code == 200, resp.content
    result = resp.json()["data"]["result"]
    assert result["candidates"] == []
    assert "not the parent" in result["reason"]

    resp = _call(client, job, "fillDatasetsFromCampaign")
    assert resp.status_code == 200, resp.content
    result = resp.json()["data"]["result"]
    assert result["success"] is False and "not the parent" in result["error"]
