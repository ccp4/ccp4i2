"""Fetching an EMDB map into a project over the API, with EBI replaced by
fixtures and a locally generated map."""

import json
import pathlib

import pytest
from rest_framework.test import APIClient

gemmi = pytest.importorskip("gemmi", reason="needs gemmi to write a test map")

from ccp4i2.db import models
from ccp4i2.lib.utils.files import repository_fetch as repo

API = "/api/ccp4i2"
FIXTURES = pathlib.Path(__file__).parents[2] / "unit" / "lib" / "fixtures" / "emdb"


def entry(name):
    return json.loads((FIXTURES / f"{name}.json").read_text())


@pytest.fixture
def client(bypass_api_permissions):
    return APIClient()


@pytest.fixture
def project(tmp_path):
    directory = tmp_path / "emdb"
    directory.mkdir()
    return models.Project.objects.create(name="emdb", directory=str(directory))


@pytest.fixture
def offline_emdb(monkeypatch):
    """EBI replaced: the entry JSON comes from fixtures and the 'download'
    writes a small real CCP4 map (so type detection sees a map)."""
    downloads = []

    def fetch_entry(entry_id):
        try:
            return entry(entry_id)
        except FileNotFoundError:
            raise repo.RepositoryError(404, f"{entry_id} is not an EMDB entry")

    def download(url, target, gunzip):
        downloads.append((url, gunzip))
        grid = gemmi.FloatGrid(8, 8, 8)
        grid.set_unit_cell(gemmi.UnitCell(16, 16, 16, 90, 90, 90))
        grid.spacegroup = gemmi.find_spacegroup_by_name("P1")
        ccp4 = gemmi.Ccp4Map()
        ccp4.grid = grid
        ccp4.update_ccp4_header()
        ccp4.write_ccp4_map(str(target))
        return target

    monkeypatch.setattr(repo, "fetch_emdb_entry", fetch_entry)
    monkeypatch.setattr(repo, "download_to", download)
    return downloads


def _import_map_job(client, project):
    response = client.post(f"{API}/projects/{project.id}/create_task/",
                           data=json.dumps({"task_name": "ImportMap"}), content_type="application/json")
    assert response.status_code == 200, response.content
    return models.Job.objects.get(id=response.json()["data"]["new_job"]["id"])


def test_listing_an_entry(client, offline_emdb):
    response = client.get(f"{API}/repositories/emdb/11638/")
    assert response.status_code == 200, response.content
    data = response.json()["data"]
    assert data["entry"] == "EMD-11638"
    assert [f["label"] for f in data["files"]] == ["Main map", "Half map 1", "Half map 2", "Mask 1"]
    assert client.get(f"{API}/repositories/emdb/EMD-99999999/").status_code == 400  # not an id shape
    assert client.get(f"{API}/repositories/pdbe/1cbs/").status_code == 400


def test_fetch_half_map_into_a_parameter(client, project, offline_emdb):
    job = _import_map_job(client, project)
    response = client.post(f"{API}/jobs/{job.id}/fetch_repository_file/", data=json.dumps({
        "object_path": "ImportMap.inputData.MAPIN",
        "repository": "emdb", "entry": "EMD-11638",
        "file": "emd_11638_half_map_1.map.gz", "sub_type": 5,
    }), content_type="application/json")
    assert response.status_code == 200, response.content
    data = response.json()["data"]
    assert data["source_url"].endswith("/EMD-11638/other/emd_11638_half_map_1.map.gz")
    assert data["kind"] == "half_map"
    assert offline_emdb == [(data["source_url"], True)]  # gunzipped on the way

    the_file = models.File.objects.get(job=job, job_param_name="MAPIN")
    assert the_file.sub_type == 5
    assert the_file.type_id == "application/CCP4-map"
    assert the_file.annotation == "EMD-11638 half map 1, 0.53 A/px, 256^3, 1.22 A"
    assert the_file.name.endswith(".map") and not the_file.name.endswith(".gz")
    assert the_file.path.exists()
    fi = models.FileImport.objects.get(file=the_file)
    assert fi.description.startswith("Fetched from https://ftp.ebi.ac.uk/")
    # the scratch directory is gone
    assert not list((pathlib.Path(project.directory) / "CCP4_IMPORTED_FILES").glob(".fetch-*"))
    assert data["updated_item"]["_value"]["baseName"]["_value"] == the_file.name


def test_fetch_refuses_wrong_subtype_unlisted_file_and_unknown_repository(client, project, offline_emdb):
    job = _import_map_job(client, project)
    post = lambda body: client.post(f"{API}/jobs/{job.id}/fetch_repository_file/",
                                    data=json.dumps({"object_path": "ImportMap.inputData.MAPIN",
                                                     "repository": "emdb", "entry": "EMD-11638", **body}),
                                    content_type="application/json")
    r = post({"file": "emd_11638_half_map_1.map.gz", "sub_type": 1})
    assert r.status_code == 400 and "half map" in r.json()["error"]
    r = post({"file": "emd_11638_half_map_3.map.gz"})
    assert r.status_code == 400 and "does not list" in r.json()["error"]
    r = post({"file": "emd_11638.map.gz", "repository": "pdbe"})
    assert r.status_code == 400
    assert offline_emdb == []  # nothing was downloaded for any refusal


def test_fetch_main_map_takes_its_own_subtype(client, project, offline_emdb):
    job = _import_map_job(client, project)
    response = client.post(f"{API}/jobs/{job.id}/fetch_repository_file/", data=json.dumps({
        "object_path": "ImportMap.inputData.MAPIN", "entry": "30210", "file": "emd_30210.map.gz",
    }), content_type="application/json")
    assert response.status_code == 200, response.content
    assert models.File.objects.get(job=job, job_param_name="MAPIN").sub_type == 1
