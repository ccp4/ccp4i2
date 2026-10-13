"""Upload refusals say what to do next.

Agents (gpt-oss:20b, 2026-10-13) stopped at three refusals that named
nothing they could act on: an upload to a list item rather than its file,
unmerged data offered to a merged-data input (gemmi's ensure_asu message),
and a raw import refused for its kind, which read as "this file can never
be used here".
"""
from pathlib import Path

import pytest
from django.core.files.uploadedfile import SimpleUploadedFile
from rest_framework.test import APIClient

import ccp4i2
from ccp4i2.db import models

API = "/api/ccp4i2"
DEMO = Path(ccp4i2.__file__).parent / "demo_data"


@pytest.fixture
def client(bypass_api_permissions):
    return APIClient()


@pytest.fixture
def project(tmp_path):
    directory = tmp_path / "uploaderrors"
    (directory / "CCP4_JOBS").mkdir(parents=True)
    return models.Project.objects.create(name="uploaderrors", directory=str(directory))


def _job(client, project, task):
    r = client.post(f"{API}/projects/{project.id}/create_task/", data={"task_name": task}, format="json")
    assert r.status_code == 200, r.content
    return models.Job.objects.get(id=r.json()["data"]["new_job"]["id"])


def _upload(client, job, path, source):
    return client.post(f"{API}/jobs/{job.id}/upload_file_param/",
                       {"file": SimpleUploadedFile(source.name, source.read_bytes()), "objectPath": path},
                       format="multipart")


def test_an_item_is_not_a_file_and_its_file_is_named(client, project):
    job = _job(client, project, "aimless_pipe")
    r = _upload(client, job, "inputData.UNMERGEDFILES[0]", DEMO / "gamma" / "gamma_native.mtz")
    assert r.status_code >= 400
    assert "UNMERGEDFILES[0].file" in r.json()["error"]


def test_unmerged_data_for_a_merged_input_points_to_aimless(client, project):
    pytest.importorskip("gemmi")
    job = _job(client, project, "molrep_mr")
    r = _upload(client, job, "inputData.F_SIGF", DEMO / "gamma" / "gamma_native.mtz")
    assert r.status_code >= 400
    assert "unmerged" in r.json()["error"] and "aimless_pipe" in r.json()["error"]
