"""Creating a Phaser step in the context of the step before it carries the
search models across, through the same population that fills file inputs.
The PHIL tasks need libtbx, so this skips CCP4-free."""
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
    directory = tmp_path / "inherit"
    directory.mkdir()
    return models.Project.objects.create(name="inherit", directory=str(directory))


def _create_job(client, project, task, **extra):
    response = client.post(f"{API_PREFIX}/projects/{project.id}/create_task/",
                           data=json.dumps({"task_name": task, **extra}), content_type="application/json")
    assert response.status_code == 200, response.content
    return models.Job.objects.get(id=response.json()["data"]["new_job"]["id"])


def _loaded(job):
    plugin = get_plugin_class(job.task_name)(workDirectory=str(job.directory))
    plugin.loadDataFromXml(str(Path(job.directory) / "input_params.xml"))
    return plugin


def test_ftf_created_from_frf_has_its_models(client, project):
    frf = _create_job(client, project, "phaser_mr_frf_phil", auto_context=False)
    plugin = get_plugin_class("phaser_mr_frf_phil")(workDirectory=str(frf.directory))
    inp = plugin.container.inputData
    inp.ENSEMBLES.append(inp.ENSEMBLES.makeItem())
    inp.ENSEMBLES[-1].label.set("beta")
    inp.ENSEMBLES[-1].number.set(1)
    item = inp.ENSEMBLES[-1].pdbItemList.makeItem()
    inp.ENSEMBLES[-1].pdbItemList.append(item)
    item.structure.setFullPath(str(GAMMA / "gamma_model.pdb"))
    item.identity_to_target.set(0.9)
    inp.COMP_BY.set("SOLVENT")
    inp.SOLVENT_FRACTION.set(0.42)
    save_params_for_job(plugin, frf)

    ftf = _create_job(client, project, "phaser_mr_ftf_phil", context_job_uuid=str(frf.uuid))
    got = _loaded(ftf).container.inputData
    assert [str(e.label) for e in got.ENSEMBLES] == ["beta"]
    assert str(got.ENSEMBLES[0].pdbItemList[0].structure.fullPath) == str(GAMMA / "gamma_model.pdb")
    assert str(got.COMP_BY) == "SOLVENT" and float(got.SOLVENT_FRACTION) == 0.42


def test_changing_the_context_job_later_inherits_too(client, project):
    frf = _create_job(client, project, "phaser_mr_frf_phil", auto_context=False)
    plugin = get_plugin_class("phaser_mr_frf_phil")(workDirectory=str(frf.directory))
    inp = plugin.container.inputData
    inp.ENSEMBLES.append(inp.ENSEMBLES.makeItem())
    inp.ENSEMBLES[-1].label.set("blip")
    save_params_for_job(plugin, frf)
    pak = _create_job(client, project, "phaser_mr_pak_phil", auto_context=False)
    assert len(_loaded(pak).container.inputData.ENSEMBLES) == 0
    response = client.post(f"{API_PREFIX}/jobs/{pak.id}/set_context_job/",
                           data=json.dumps({"context_job_uuid": str(frf.uuid)}), content_type="application/json")
    assert response.status_code in (200, 201), response.content
    assert [str(e.label) for e in _loaded(pak).container.inputData.ENSEMBLES] == ["blip"]
