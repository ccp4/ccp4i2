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


def test_task_lookup_says_which_tasks_have_a_judgement(client):
    # An agent choosing between a pipeline and its stages sees which tasks
    # have been written up (Haiku twice started SAD phasing with ShelxCD)
    lookup = client.get(f"{API}/task_lookup/").json()
    assert lookup["shelx"]["hasJudgement"] is True
    assert lookup["ShelxCD"]["hasJudgement"] is True
    assert lookup["editbfac"]["hasJudgement"] is False


def test_a_value_of_the_wrong_shape_is_refused_when_set(client, project):
    # An agent's ASU_CONTENT sequence given as {"text": ...} was stored as the
    # dict's repr and failed only at validation
    import json as _json
    job = _job(client, project)
    r = client.post(f"{API}/jobs/{job.id}/set_parameter/", data=_json.dumps({
        "object_path": "ProvideAsuContents.inputData.ASU_CONTENT",
        "value": [{"name": "A", "sequence": {"text": "MKV"}, "nCopies": 1}]}),
        content_type="application/json")
    assert r.status_code >= 400, r.content
    assert "ASU_CONTENT[0].sequence takes a single value" in r.content.decode()


def test_a_report_is_rendered_when_the_judgement_reads_it(client, project, monkeypatch):
    # MrBUMP's quick mode writes no program.xml; its judgement reads the
    # report, which existed only once someone had opened it, so a job an
    # agent ran itself judged "unjudged"
    from pathlib import Path
    from ccp4i2.lib.utils.jobs import reports
    rendered = []

    def fake_report(job, regenerate=False):
        rendered.append(job.id)
        (Path(job.directory) / "report_xml.xml").write_text("<report><title>MrBUMP</title></report>")

    monkeypatch.setattr(reports, "get_job_report_xml", fake_report)
    job = _job(client, project, task="mrbump_basic")
    Path(job.directory).mkdir(parents=True, exist_ok=True)
    job.status = models.Job.Status.PENDING
    job.save()
    client.get(f"{API}/jobs/{job.id}/judgement/")
    assert rendered == []  # not before it has finished
    job.status = models.Job.Status.FINISHED
    job.save()
    client.get(f"{API}/jobs/{job.id}/judgement/")
    assert rendered == [job.id]
    client.get(f"{API}/jobs/{job.id}/judgement/")
    assert rendered == [job.id]  # once: then it is there
