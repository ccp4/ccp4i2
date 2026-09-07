"""Cloning a job of a superseded task makes a job of its successor, with
the old job's front page adopted. The project is built from nothing; the
PHIL successor needs libtbx, so this skips CCP4-free."""
import json
from pathlib import Path

import pytest
from rest_framework.test import APIClient

pytest.importorskip("libtbx.phil", reason="needs libtbx (CCP4/cctbx)")

import ccp4i2
from ccp4i2.core.tasks import get_plugin_class
from ccp4i2.db import models
from ccp4i2.lib.utils.parameters.save_params import save_params_for_job

API_PREFIX = "/api/ccp4i2"
GAMMA = Path(ccp4i2.__file__).parent / "demo_data" / "gamma"


@pytest.fixture
def client(bypass_api_permissions):
    return APIClient()


@pytest.fixture
def project(tmp_path):
    directory = tmp_path / "successor"
    directory.mkdir()
    return models.Project.objects.create(name="successor", directory=str(directory))


def _create_job(client, project, task):
    response = client.post(f"{API_PREFIX}/projects/{project.id}/create_task/",
                           data=json.dumps({"task_name": task}), content_type="application/json")
    assert response.status_code == 200, response.content
    return models.Job.objects.get(id=response.json()["data"]["new_job"]["id"])


def test_clone_of_a_classic_phaser_job_is_a_phil_job(client, project):
    job = _create_job(client, project, "phaser_simple")
    plugin = get_plugin_class("phaser_simple")(workDirectory=str(job.directory))
    inp = plugin.container.inputData
    inp.XYZIN.setFullPath(str(GAMMA / "gamma_model.pdb"))
    inp.NCOPIES.set(2)
    inp.RESOLUTION_HIGH.set(2.5)
    save_params_for_job(plugin, job)

    response = client.post(f"{API_PREFIX}/jobs/{job.id}/clone/")
    assert response.status_code == 201, response.content
    clone = response.json()
    assert clone["task_name"] == "phaser_simple_phil"
    assert clone["status"] == models.Job.Status.PENDING
    new_job = models.Job.objects.get(id=clone["id"])
    assert new_job.title != job.title

    successor = get_plugin_class("phaser_simple_phil")(workDirectory=str(new_job.directory))
    successor.loadDataFromXml(str(Path(new_job.directory) / "input_params.xml"))
    assert str(successor.container.inputData.XYZIN.fullPath) == str(GAMMA / "gamma_model.pdb")
    assert int(successor.container.inputData.NCOPIES) == 2
    assert float(successor.find_phil("phaser.keywords.resolution.high")) == 2.5


def test_clone_of_a_current_task_stays_itself(client, project):
    job = _create_job(client, project, "phaser_simple_phil")
    response = client.post(f"{API_PREFIX}/jobs/{job.id}/clone/")
    assert response.status_code == 201, response.content
    assert response.json()["task_name"] == "phaser_simple_phil"


def test_task_lookup_says_who_supersedes(client):
    lookup = client.get(f"{API_PREFIX}/task_lookup/").json()
    assert lookup["phaser_simple"]["supersededBy"] == "phaser_simple_phil"
    assert lookup["phaser_simple_phil"]["supersededBy"] is None
