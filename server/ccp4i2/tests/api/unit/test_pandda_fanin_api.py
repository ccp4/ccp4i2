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


@pytest.fixture
def campaign(test_project_path):
    """A campaign whose one member has a finished dimple job with its outputs
    on disk but unregistered, as every job run before COMPLETE_MTZ tracking has."""
    test_project_path.mkdir(exist_ok=True)

    def make_project(name):
        directory = test_project_path / name
        (directory / "CCP4_JOBS").mkdir(parents=True, exist_ok=True)
        return models.Project.objects.create(name=name, directory=str(directory))

    parent, member = make_project("fanin_ref"), make_project("fanin_x1")
    group = models.ProjectGroup.objects.create(name="fanin_camp", type="fragment_set")
    models.ProjectGroupMembership.objects.create(
        group=group, project=parent, type=models.ProjectGroupMembership.MembershipType.PARENT)
    models.ProjectGroupMembership.objects.create(
        group=group, project=member, type=models.ProjectGroupMembership.MembershipType.MEMBER)
    dimple = models.Job.objects.create(project=member, number="1.3", task_name="i2Dimple",
                                       title="dimple", status=models.Job.Status.FINISHED)
    dimple.directory.mkdir(parents=True)
    (dimple.directory / "final.pdb").write_text("END\n")
    (dimple.directory / "final.mtz").write_bytes(b"MTZ ")
    job = models.Job.objects.get(uuid=create_job(projectId=str(parent.uuid), taskName="pandda_campaign"))
    return {"job": job, "member": member}


def _header(job):
    import xml.etree.ElementTree as ET
    root = ET.parse(job.directory / "input_params.xml").getroot()
    return {child.tag: child.text for child in root.find("ccp4i2_header")}


def test_a_fill_through_object_method_keeps_the_job_identity(client, campaign):
    """The fill saves the parameters back. Built by the generic endpoint, the
    plugin did not know which job it was, so the saved header lost its jobId
    and every later campaignCandidates call answered "not in the database":
    the interface's fill button went dead after its first use."""
    job = campaign["job"]
    assert _header(job).get("jobId") == str(job.uuid)

    resp = _call(client, job, "fillDatasetsFromCampaign")
    assert resp.status_code == 200, resp.content
    assert resp.json()["data"]["result"]["data"]["added"] == ["fanin_x1"]
    assert _header(job).get("jobId") == str(job.uuid)

    resp = _call(client, job, "campaignCandidates")
    result = resp.json()["data"]["result"]
    assert result["reason"] is None, result
    assert [c["listed"] for c in result["candidates"]] == [True]
