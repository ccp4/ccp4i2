"""The recorded Moorhen session over the REST API, with no browser.

Run on an interactive task opens a session and dispatches nothing; the
window heartbeats and drops files; Finish dispatches the job (or marks it
for deletion if nothing was saved); the runner then harvests and gleans.
"""

import json
from pathlib import Path

import pytest
from asgiref.sync import async_to_sync
from rest_framework.test import APIClient

import ccp4i2
from ccp4i2.db import models
from ccp4i2.lib.utils.jobs import interactive

API = "/api/ccp4i2"
GAMMA = Path(ccp4i2.__file__).parent / "demo_data" / "gamma"
MODEL = GAMMA / "gamma_model.pdb"


@pytest.fixture
def client(bypass_api_permissions):
    return APIClient()


@pytest.fixture
def project(tmp_path):
    directory = tmp_path / "moorhen"
    directory.mkdir()
    return models.Project.objects.create(name="moorhen", directory=str(directory))


@pytest.fixture
def no_dispatch(monkeypatch):
    """Record dispatches instead of spawning a runner process."""
    calls = []

    def fake_run_job_local(job, synchronous=False):
        calls.append(job.id)
        return {"success": True, "data": job}

    monkeypatch.setattr("ccp4i2.lib.dispatch.local.run_job_local", fake_run_job_local)
    return calls


def _create_moorhen_job(client, project):
    response = client.post(f"{API}/projects/{project.id}/create_task/",
                           data=json.dumps({"task_name": "moorhen"}),
                           content_type="application/json")
    assert response.status_code == 200, response.content
    job = models.Job.objects.get(id=response.json()["data"]["new_job"]["id"])
    with open(MODEL, "rb") as handle:
        response = client.post(
            f"{API}/jobs/{job.id}/upload_file_param/",
            {"file": handle, "object_path": "moorhen.inputData.XYZIN_LIST[0]"},
            format="multipart")
    assert response.status_code == 200, response.content
    return job


def _run(client, job):
    response = client.post(f"{API}/jobs/{job.id}/run/")
    assert response.status_code == 200, response.content
    job.refresh_from_db()
    return job


def test_run_opens_a_session_and_dispatches_nothing(client, project, no_dispatch):
    job = _create_moorhen_job(client, project)
    job = _run(client, job)

    assert job.status == models.Job.Status.RUNNING
    assert job.process_id is None
    assert no_dispatch == []
    session = job.interactive_session
    assert session.finished is False and session.dispatched is False
    assert (job.directory / interactive.DROP_DIR_NAME).is_dir()

    state = client.get(f"{API}/jobs/{job.id}/interactive_session/").json()["data"]
    assert state["session"]["attached"] is False
    plan = state["load_plan"]
    assert [entry["kind"] for entry in plan] == ["coordinates"]
    assert plan[0]["file_id"] is not None
    assert plan[0]["type"] == "chemical/x-pdb"
    assert state["outputs"] == []

    # Running again while the session is open just reopens it.
    assert client.post(f"{API}/jobs/{job.id}/run/").status_code == 200
    assert no_dispatch == []


def test_session_endpoints_refuse_a_job_without_an_open_session(client, project, no_dispatch):
    job = _create_moorhen_job(client, project)
    with open(MODEL, "rb") as handle:
        response = client.post(f"{API}/jobs/{job.id}/interactive_drop/",
                               {"file": handle}, format="multipart")
    assert response.status_code == 409
    assert client.post(f"{API}/jobs/{job.id}/interactive_heartbeat/").status_code == 409
    assert client.post(f"{API}/jobs/{job.id}/interactive_finish/").status_code == 409


def test_session_endpoints_refuse_a_non_interactive_task(client, project):
    response = client.post(f"{API}/projects/{project.id}/create_task/",
                           data=json.dumps({"task_name": "coot1"}),
                           content_type="application/json")
    job_id = response.json()["data"]["new_job"]["id"]
    assert client.get(f"{API}/jobs/{job_id}/interactive_session/").status_code == 400


def test_heartbeat_marks_the_window_attached(client, project, no_dispatch):
    job = _run(client, _create_moorhen_job(client, project))
    state = client.post(f"{API}/jobs/{job.id}/interactive_heartbeat/").json()["data"]
    assert state["session"]["attached"] is True


def test_drop_writes_the_output_contract_with_a_sidecar(client, project, no_dispatch):
    job = _run(client, _create_moorhen_job(client, project))
    with open(MODEL, "rb") as handle:
        response = client.post(
            f"{API}/jobs/{job.id}/interactive_drop/",
            {"file": handle, "kind": "model", "annotation": "fitted the ligand"},
            format="multipart")
    assert response.status_code == 200, response.content
    dropped = response.json()["data"]
    assert dropped["name"] == f"output{dropped['number']}.pdb"  # the Coot contract

    drop_dir = job.directory / interactive.DROP_DIR_NAME
    assert (drop_dir / dropped["name"]).read_bytes() == MODEL.read_bytes()
    meta = json.loads((drop_dir / f"output{dropped['number']}.meta.json").read_text())
    assert meta["annotation"] == "fitted the ligand"

    outputs = client.get(f"{API}/jobs/{job.id}/interactive_session/").json()["data"]["outputs"]
    assert [o["annotation"] for o in outputs] == ["fitted the ligand"]


def test_finish_with_nothing_saved_marks_the_job_for_deletion(client, project, no_dispatch):
    job = _run(client, _create_moorhen_job(client, project))
    response = client.post(f"{API}/jobs/{job.id}/interactive_finish/")
    assert response.status_code == 200, response.content
    assert response.json()["data"]["disposition"] == "deleted"
    job.refresh_from_db()
    assert job.status == models.Job.Status.TO_DELETE
    assert job.interactive_session.finished is True
    assert no_dispatch == []


def test_closing_a_window_with_saved_work_keeps_the_session_open(client, project, no_dispatch):
    job = _run(client, _create_moorhen_job(client, project))
    with open(MODEL, "rb") as handle:
        client.post(f"{API}/jobs/{job.id}/interactive_drop/", {"file": handle}, format="multipart")
    response = client.post(f"{API}/jobs/{job.id}/interactive_finish/",
                           data=json.dumps({"finished": False}),
                           content_type="application/json")
    assert response.json()["data"]["disposition"] == "kept_open"
    job.refresh_from_db()
    assert job.status == models.Job.Status.RUNNING
    assert job.interactive_session.finished is False
    assert no_dispatch == []


def test_finish_with_saved_work_dispatches_the_job_once(client, project, no_dispatch):
    job = _run(client, _create_moorhen_job(client, project))
    with open(MODEL, "rb") as handle:
        client.post(f"{API}/jobs/{job.id}/interactive_drop/", {"file": handle}, format="multipart")
    response = client.post(f"{API}/jobs/{job.id}/interactive_finish/")
    assert response.status_code == 200, response.content
    assert response.json()["data"]["disposition"] == "dispatched"
    assert no_dispatch == [job.id]
    job.refresh_from_db()
    assert job.interactive_session.finished is True
    assert job.interactive_session.dispatched is True
    # Finishing again is a no-op, not a second dispatch.
    assert client.post(f"{API}/jobs/{job.id}/interactive_finish/").json()["data"]["disposition"] == "already_finished"
    assert no_dispatch == [job.id]


def test_cancel_ends_the_session_without_harvesting(client, project, no_dispatch):
    job = _run(client, _create_moorhen_job(client, project))
    response = client.post(f"{API}/jobs/{job.id}/cancel/")
    assert response.status_code == 200, response.content
    job.refresh_from_db()
    assert job.status == models.Job.Status.INTERRUPTED
    assert job.interactive_session.finished is True


def test_dispatched_job_harvests_and_gleans_the_saved_model(client, project, no_dispatch):
    """The runner, run in-process after Finish, sees a finished session,
    returns at once from startProcess, harvests the drop directory and
    gleans XYZOUT with the annotation the window gave."""
    from ccp4i2.lib.async_run_job import run_job_async

    job = _run(client, _create_moorhen_job(client, project))
    with open(MODEL, "rb") as handle:
        client.post(f"{API}/jobs/{job.id}/interactive_drop/",
                    {"file": handle, "annotation": "built by hand"}, format="multipart")
    assert client.post(f"{API}/jobs/{job.id}/interactive_finish/").json()["data"]["disposition"] == "dispatched"

    async_to_sync(run_job_async)(job.uuid)

    job.refresh_from_db()
    assert job.status == models.Job.Status.FINISHED, job.get_status_display()
    assert (job.directory / "XYZOUT_0.pdb").exists()
    outputs = models.File.objects.filter(job=job, directory=models.File.Directory.JOB_DIR)
    assert [f.annotation for f in outputs] == ["built by hand"]
    assert outputs[0].job_param_name.startswith("XYZOUT")
