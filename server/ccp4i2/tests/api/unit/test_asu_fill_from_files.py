"""ProvideAsuContents fills its table from a file set by any route.

The page used to do it in the browser after setting SEQIN, XYZIN or
ASUCONTENTIN, so a job set up through the API had the file and an empty
table; agents gave up there (2026-10-13). The server does it now, the way the
page did: a sequence already in the table gains copies, not a second row.
"""
from pathlib import Path

import pytest
from django.core.files.uploadedfile import SimpleUploadedFile
from rest_framework.test import APIClient

import ccp4i2
from ccp4i2.db import models
from ccp4i2.lib.utils.sequences.asu_entries import (
    chains_from_digest, deduplicate_chains, entries_from_digest, merge_asu_entries)

API = "/api/ccp4i2"
DEMO = Path(ccp4i2.__file__).parent / "demo_data"


# -- the rules, as the page had them ----------------------------------------

def test_identical_chains_become_one_row_with_copies():
    digest = {"composition": {"peptides": ["A", "B", "C"]},
              "sequences": {"A": "MKV", "B": "MKV", "C": "GGS"}}
    entries = deduplicate_chains(chains_from_digest(digest, "m.pdb"))
    assert [(e["name"], e["nCopies"]) for e in entries] == [("Chain_A_B", 2), ("Chain_C", 1)]


def test_a_sequence_already_present_gains_copies():
    existing = [{"name": "a", "sequence": "MKV", "polymerType": "PROTEIN", "nCopies": 1}]
    merged = merge_asu_entries(existing, [{"name": "b", "sequence": "M K V", "polymerType": "protein",
                                           "nCopies": 2}, {"name": "c", "sequence": "GGS",
                                                           "polymerType": "PROTEIN", "nCopies": 1}])
    assert [(e["name"], e["nCopies"]) for e in merged] == [("a", 3), ("c", 1)]


def test_a_sequence_digest_and_an_asu_digest():
    entries, chains = entries_from_digest({"moleculeType": "PROTEIN", "sequence": "MK V",
                                           "name": "gamma"}, "gamma.pir")
    assert chains == 1 and entries[0]["sequence"] == "MKV" and entries[0]["polymerType"] == "PROTEIN"
    entries, chains = entries_from_digest({"sequences": [{"name": "x", "sequence": "MKV",
                                                          "polymerType": "PROTEIN", "nCopies": 4}]}, "f")
    assert chains == 1 and entries[0]["nCopies"] == 4


# -- every route through the API --------------------------------------------

@pytest.fixture
def client(bypass_api_permissions):
    return APIClient()


@pytest.fixture
def job(client, tmp_path):
    directory = tmp_path / "asufill"
    (directory / "CCP4_JOBS").mkdir(parents=True)
    project = models.Project.objects.create(name="asufill", directory=str(directory))
    r = client.post(f"{API}/projects/{project.id}/create_task/",
                    data={"task_name": "ProvideAsuContents"}, format="json")
    assert r.status_code == 200, r.content
    return models.Job.objects.get(id=r.json()["data"]["new_job"]["id"])


def _upload(client, job, path, source):
    r = client.post(f"{API}/jobs/{job.id}/upload_file_param/",
                    {"file": SimpleUploadedFile(source.name, source.read_bytes()), "objectPath": path},
                    format="multipart")
    assert r.status_code == 200, r.content
    return r.json()["data"]


def _table(client, job):
    r = client.get(f"{API}/jobs/{job.id}/parameters/", {"query": "ASU_CONTENT"})
    entry = next(p for p in r.json()["data"]["parameters"] if p["path"] == "inputData.ASU_CONTENT")
    return entry.get("value") or []


def test_a_sequence_file_fills_the_table_and_again_adds_copies(client, job):
    data = _upload(client, job, "inputData.SEQIN", DEMO / "gamma" / "gamma.pir")
    assert "1 new row" in data["note"]
    table = _table(client, job)
    assert len(table) == 1 and len(table[0]["sequence"]) > 50 and table[0]["nCopies"] == 1
    data = _upload(client, job, "inputData.SEQIN", DEMO / "gamma" / "gamma.pir")
    assert "copies" in data["note"]
    assert [e["nCopies"] for e in _table(client, job)] == [2]


def test_a_model_and_an_asu_file_fill_it_too(client, job):
    data = _upload(client, job, "inputData.XYZIN", DEMO / "beta_blip" / "blip.pdb")
    assert "new row" in data["note"]
    data = _upload(client, job, "inputData.ASUCONTENTIN", DEMO / "gamma" / "gamma.asu.xml")
    assert "new row" in data["note"]
    assert len(_table(client, job)) == 2


def test_a_row_with_no_sequence_gives_way(client, job):
    r = client.post(f"{API}/jobs/{job.id}/add_list_item/",
                    data={"object_path": "inputData.ASU_CONTENT"}, format="json")
    assert r.status_code == 200, r.content
    _upload(client, job, "inputData.SEQIN", DEMO / "gamma" / "gamma.pir")
    table = _table(client, job)
    assert len(table) == 1 and table[0]["sequence"]


def test_a_project_file_picked_fills_it(client, job):
    # The page's pulldown sets a project file through set_parameter
    _upload(client, job, "inputData.SEQIN", DEMO / "gamma" / "gamma.pir")
    r = client.get(f"{API}/projects/{job.project.id}/resolve_fileuse/",
                   {"fileuse": f"[{job.number}].SEQIN"})
    value = r.json()["data"]
    value.pop("fullPath", None)
    r = client.post(f"{API}/projects/{job.project.id}/create_task/",
                    data={"task_name": "ProvideAsuContents"}, format="json")
    other = models.Job.objects.get(id=r.json()["data"]["new_job"]["id"])
    r = client.post(f"{API}/jobs/{other.id}/set_parameter/",
                    data={"object_path": "inputData.SEQIN", "value": value}, format="json")
    assert r.status_code == 200, r.content
    assert "new row" in r.json()["data"]["note"]
    assert len(_table(client, other)) == 1
