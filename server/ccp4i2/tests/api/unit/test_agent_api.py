"""The endpoints an agent works from (docs/agentic-knowledge.md): a job's
parameters one line each, a task's judgement, and the verdict on a job."""
import pytest
import yaml
from rest_framework.test import APIClient

from ccp4i2.agent import judgement
from ccp4i2.db import models

API = "/api/ccp4i2"


@pytest.fixture
def client(bypass_api_permissions):
    return APIClient()


@pytest.fixture
def project(tmp_path):
    directory = tmp_path / "agentapi"
    (directory / "CCP4_JOBS").mkdir(parents=True)
    return models.Project.objects.create(name="agentapi", directory=str(directory))


def _job(client, project, task="ProvideAsuContents"):
    r = client.post(f"{API}/projects/{project.id}/create_task/",
                    data={"task_name": task}, format="json")
    assert r.status_code == 200, r.content
    return models.Job.objects.get(id=r.json()["data"]["new_job"]["id"])


def test_parameters_are_one_line_each(client, project):
    job = _job(client, project)
    r = client.get(f"{API}/jobs/{job.id}/parameters/")
    assert r.status_code == 200, r.content
    params = r.json()["data"]["parameters"]
    by_path = {p["path"]: p for p in params}
    asu = by_path["inputData.ASU_CONTENT"]
    assert asu["class"] == "CAsuContentSeqList" and asu["set"] is False
    assert all(p["path"].split(".")[0] in ("inputData", "controlParameters", "outputData")
               for p in params)

    only_input = client.get(f"{API}/jobs/{job.id}/parameters/?section=inputData").json()
    assert {p["path"].split(".")[0] for p in only_input["data"]["parameters"]} == {"inputData"}
    found = client.get(f"{API}/jobs/{job.id}/parameters/?query=asu_content").json()
    assert "inputData.ASU_CONTENT" in [p["path"] for p in found["data"]["parameters"]]


def test_a_task_without_judgement_says_so(client, project):
    r = client.get(f"{API}/agent/tasks/freerflag/")  # a task nobody has written up
    assert r.status_code == 200
    assert r.json()["data"]["judgement"] is None
    assert client.get(f"{API}/agent/tasks/no_such_task/").status_code == 404

    job = _job(client, project, task="freerflag")
    verdict = client.get(f"{API}/jobs/{job.id}/judgement/").json()["data"]
    assert verdict["outcome"] is None and "No judgement" in verdict["note"]


def test_judgement_is_read_from_the_job(client, project, tmp_path, monkeypatch):
    job = _job(client, project)
    (job.directory / "program.xml").write_text("<R><Final><RFree>0.24</RFree></Final></R>")
    key, _ = models.JobValueKey.objects.get_or_create(name="RFactor", defaults={"description": "R"})
    models.JobFloatValue.objects.create(job=job, key=key, value=0.2)
    path = tmp_path / "ProvideAsuContents.agent.yaml"
    path.write_text(yaml.safe_dump({
        "task": "ProvideAsuContents", "status": "draft",
        "results": {"RFREE": {"xpath": ".//Final/RFree"},
                    "R": {"file": "kpi", "kpi": "RFactor"}},
        "verdict": [{"when": "RFREE < 0.3 and R < RFREE", "outcome": "good", "basis": "test"},
                    {"when": True, "outcome": "poor"}],
        "next": [{"when": 'outcome == "good"', "task": "servalcat_pipe"}],
    }))
    monkeypatch.setattr(judgement, "judgement_path", lambda task: path)

    described = client.get(f"{API}/agent/tasks/ProvideAsuContents/").json()["data"]
    assert described["judgement"]["results"]["RFREE"]["xpath"] == ".//Final/RFree"

    verdict = client.get(f"{API}/jobs/{job.id}/judgement/").json()["data"]
    assert verdict["results"] == {"RFREE": 0.24, "R": 0.2}
    assert verdict["outcome"] == "good"
    assert verdict["next"] == [{"when": 'outcome == "good"', "task": "servalcat_pipe"}]
    assert verdict["note"] == judgement.DRAFT_NOTE
