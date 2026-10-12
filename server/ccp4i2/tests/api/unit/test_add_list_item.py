"""add_list_item: a list is built an item at a time, as the interface's
"add" button builds it, down to a list inside an item (a Phaser ensemble's
models)."""
import pytest
from rest_framework.test import APIClient

from ccp4i2.db import models

API = "/api/ccp4i2"


@pytest.fixture
def client(bypass_api_permissions):
    return APIClient()


@pytest.fixture
def project(tmp_path):
    directory = tmp_path / "addlistitem"
    (directory / "CCP4_JOBS").mkdir(parents=True)
    return models.Project.objects.create(name="addlistitem", directory=str(directory))


def _job(client, project, task):
    r = client.post(f"{API}/projects/{project.id}/create_task/", data={"task_name": task}, format="json")
    assert r.status_code == 200, r.content
    return models.Job.objects.get(id=r.json()["data"]["new_job"]["id"])


def _add(client, job, path, values=None):
    return client.post(f"{API}/jobs/{job.id}/add_list_item/",
                       data={"object_path": path, **({"values": values} if values else {})}, format="json")


def _ensembles(client, job):
    r = client.get(f"{API}/jobs/{job.id}/parameters/", {"query": "ENSEMBLES"})
    return next(p for p in r.json()["data"]["parameters"] if p["path"] == "inputData.ENSEMBLES")


def test_an_ensemble_and_its_model_are_added_an_item_at_a_time(client, project):
    job = _job(client, project, "phaser_simple")
    before = len(_ensembles(client, job).get("value") or [])
    r = _add(client, job, "inputData.ENSEMBLES", {"label": "beta"})
    assert r.status_code == 200, r.content
    data = r.json()["data"]
    assert data["item"] == f"inputData.ENSEMBLES[{before}]" and data["length"] == before + 1
    assert f"inputData.ENSEMBLES[{before}].pdbItemList" in data["lists"]
    assert f"inputData.ENSEMBLES[{before}].label" in data["fields"]

    r = _add(client, job, f"inputData.ENSEMBLES[{before}].pdbItemList")
    assert r.status_code == 200, r.content
    model = r.json()["data"]
    assert model["item"].startswith(f"inputData.ENSEMBLES[{before}].pdbItemList[")
    assert f"{model['item']}.structure" in model["files"]

    # Saved: the job, read afresh, has the ensemble with its label
    value = _ensembles(client, job)["value"]
    assert value[before]["label"] == "beta"


def test_a_refused_value_adds_nothing(client, project):
    job = _job(client, project, "ProvideAsuContents")
    r = _add(client, job, "inputData.ASU_CONTENT", {"nonsense_field": 1})
    assert r.status_code == 400
    assert "Item not added" in r.json()["error"] and "nCopies" in r.json()["error"]
    r = client.get(f"{API}/jobs/{job.id}/parameters/", {"query": "ASU_CONTENT"})
    asu = next(p for p in r.json()["data"]["parameters"] if p["path"] == "inputData.ASU_CONTENT")
    assert not asu.get("value")


def test_only_a_list_takes_an_item(client, project):
    job = _job(client, project, "ProvideAsuContents")
    r = _add(client, job, "inputData.ASUCONTENTIN")
    assert r.status_code == 400 and "not a list" in r.json()["error"]
    r = _add(client, job, "inputData.NO_SUCH_THING")
    assert r.status_code == 400 and "is nothing" in r.json()["error"]


def test_a_finished_job_is_not_changed(client, project):
    job = _job(client, project, "ProvideAsuContents")
    job.status = models.Job.Status.FINISHED
    job.save()
    r = _add(client, job, "inputData.ASU_CONTENT")
    assert r.status_code == 400 and "clone" in r.json()["error"]


def test_an_item_with_nothing_set_is_kept_for_the_next_call(client, project):
    # Each call reloads the job from input_params.xml: an item the save
    # dropped for being unset would be gone by the next one.
    job = _job(client, project, "phaser_simple")
    before = len(_ensembles(client, job).get("value") or [])
    r = _add(client, job, "inputData.ENSEMBLES")
    assert r.status_code == 200, r.content
    r = _add(client, job, f"inputData.ENSEMBLES[{before}].pdbItemList")
    assert r.status_code == 200, r.content
    r = _add(client, job, "inputData.ENSEMBLES")
    assert r.status_code == 200, r.content
    assert r.json()["data"]["index"] == before + 1


def test_an_unset_item_of_a_simple_list_is_kept(client, project):
    job = _job(client, project, "ProvideAsuContents")
    first = _add(client, job, "inputData.ASU_CONTENT").json()["data"]
    second = _add(client, job, "inputData.ASU_CONTENT").json()["data"]
    assert second["index"] == first["index"] + 1


def test_a_model_is_uploaded_into_the_item_just_added(client, project):
    # The skipped test_api.py::test_pdb_item_file_upload could not upload to
    # ENSEMBLES[0].pdbItemList[0] because the list was empty: added first,
    # the path exists, and the file stays through the next reload.
    from pathlib import Path
    from django.core.files.uploadedfile import SimpleUploadedFile
    import ccp4i2

    job = _job(client, project, "phaser_simple")
    ensemble = _add(client, job, "inputData.ENSEMBLES", {"label": "beta"}).json()["data"]
    model = _add(client, job, f"{ensemble['item']}.pdbItemList").json()["data"]
    pdb = Path(ccp4i2.__file__).parent / "demo_data" / "beta_blip" / "beta.pdb"
    r = client.post(f"{API}/jobs/{job.id}/upload_file_param/",
                    {"file": SimpleUploadedFile("beta.pdb", pdb.read_bytes()),
                     "objectPath": f"{model['item']}.structure"}, format="multipart")
    assert r.status_code == 200, r.content
    _add(client, job, "inputData.ENSEMBLES", {"label": "blip"})  # one more reload and save
    value = _ensembles(client, job)["value"]
    index = ensemble["index"]
    assert value[index]["label"] == "beta"
    assert value[index]["pdbItemList"][0]["structure"]["file"] == "beta.pdb"
    assert value[index + 1]["label"] == "blip"
