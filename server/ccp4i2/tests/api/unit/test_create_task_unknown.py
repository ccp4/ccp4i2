"""create_task refuses a task name the registry does not know, before it
writes anything: no job row, no job directory."""
from pathlib import Path

import pytest
from rest_framework.test import APIClient

from ccp4i2.db import models

API = "/api/ccp4i2"


@pytest.fixture
def client(bypass_api_permissions):
    return APIClient()


@pytest.fixture
def project(tmp_path):
    directory = tmp_path / "unknowntask"
    (directory / "CCP4_JOBS").mkdir(parents=True)
    return models.Project.objects.create(name="unknowntask", directory=str(directory))


def test_unknown_task_is_refused_and_leaves_nothing(client, project):
    r = client.post(f"{API}/projects/{project.id}/create_task/",
                    data={"task_name": "phaser_mr"}, format="json")
    assert r.status_code == 400, r.content
    error = r.json()["error"]
    assert "phaser_mr" in error and "did you mean" in error
    assert "phaser_mr_auto_phil" in error
    assert not models.Job.objects.filter(project=project).exists()
    assert list((Path(project.directory) / "CCP4_JOBS").iterdir()) == []


def test_known_task_is_still_created(client, project):
    r = client.post(f"{API}/projects/{project.id}/create_task/",
                    data={"task_name": "ProvideAsuContents"}, format="json")
    assert r.status_code == 200, r.content
    assert models.Job.objects.filter(project=project).count() == 1
