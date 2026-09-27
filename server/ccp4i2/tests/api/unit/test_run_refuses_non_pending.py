"""POST jobs/{id}/run/ on a job that is not pending is a 409 with the reason,
never a silent no-op (a re-run of a failed job did nothing on Materia,
2026-09-27; the right move is to clone it, and the response says so)."""
import pytest
from rest_framework.test import APIClient

from ccp4i2.db import models

API = "/api/ccp4i2"


@pytest.fixture
def client(bypass_api_permissions):
    return APIClient()


@pytest.fixture
def project(tmp_path):
    directory = tmp_path / "runguard"
    (directory / "CCP4_JOBS").mkdir(parents=True)
    return models.Project.objects.create(name="runguard", directory=str(directory))


def _job(client, project):
    r = client.post(f"{API}/projects/{project.id}/create_task/",
                    data={"task_name": "ProvideAsuContents", "title": "guard"}, format="json")
    assert r.status_code == 200, r.content
    return models.Job.objects.get(id=r.json()["data"]["new_job"]["id"])


@pytest.mark.parametrize("status, advice", [
    (models.Job.Status.FAILED, "clone it"),
    (models.Job.Status.FINISHED, "clone it"),
    (models.Job.Status.RUNNING, "wait for it"),
    (models.Job.Status.RUNNING_REMOTELY, "reconcile"),
])
def test_run_refuses_a_job_that_is_not_pending(client, project, status, advice):
    job = _job(client, project)
    job.status = status
    job.save()
    for endpoint in ("run", "run_local"):
        r = client.post(f"{API}/jobs/{job.id}/{endpoint}/")
        assert r.status_code == 409, (endpoint, r.status_code, r.content)
        body = r.json()
        assert body["error"] == "job_not_runnable" and advice in body["reason"], body
    job.refresh_from_db()
    assert job.status == status                 # untouched


def test_run_starts_a_pending_job(client, project, monkeypatch):
    started = []
    monkeypatch.setattr("ccp4i2.lib.dispatch.local.run_job_local",
                        lambda job, synchronous=False: started.append(job.id) or
                        {"success": True, "data": job, "status": 200})
    job = _job(client, project)
    r = client.post(f"{API}/jobs/{job.id}/run/")
    assert r.status_code == 200, r.content
    assert started == [job.id]
