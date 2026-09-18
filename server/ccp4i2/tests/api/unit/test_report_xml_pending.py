"""GET jobs/{id}/report_xml/ on a job that is running without a program.xml
returns the calm pending report; a finished job without one still gets the
failure panel."""

import json

import pytest
from rest_framework.test import APIClient

from ccp4i2.db import models

API = "/api/ccp4i2"


@pytest.fixture
def client(bypass_api_permissions):
    return APIClient()


@pytest.fixture
def job(client, tmp_path):
    directory = tmp_path / "pending"
    directory.mkdir()
    project = models.Project.objects.create(name="pending", directory=str(directory))
    response = client.post(f"{API}/projects/{project.id}/create_task/",
                           data=json.dumps({"task_name": "molrep_map"}), content_type="application/json")
    assert response.status_code == 200, response.content
    return models.Job.objects.get(id=response.json()["data"]["new_job"]["id"])


def report_for(client, job, status):
    job.status = status
    job.save()
    response = client.get(f"{API}/jobs/{job.id}/report_xml/")
    assert response.status_code == 200, response.content
    return response.json()["data"]["xml"]


def test_running_job_without_xml_gets_the_pending_report(client, job):
    xml = report_for(client, job, models.Job.Status.RUNNING)
    assert 'reportPending="true"' in xml
    assert "PROGRAM_XML_NOT_FOUND" not in xml
    assert "errorReportList" not in xml
    assert f"Job {job.number} is running" in xml


def test_queued_job_says_queued(client, job):
    xml = report_for(client, job, models.Job.Status.QUEUED)
    assert "is queued" in xml and 'reportPending="true"' in xml


def test_finished_job_without_xml_still_gets_the_failure_panel(client, job):
    xml = report_for(client, job, models.Job.Status.FAILED)
    assert "PROGRAM_XML_NOT_FOUND" in xml
    assert 'reportFailed="true"' in xml
