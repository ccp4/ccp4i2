"""ProvideSequence takes its sequence from a file given through the API, as
the interface does in the browser: no text has to be typed to pass
validation (and with neither, it is still required)."""
from pathlib import Path

import pytest
from django.core.files.uploadedfile import SimpleUploadedFile
from rest_framework.test import APIClient

import ccp4i2
from ccp4i2.db import models
from ccp4i2.lib.utils.plugins.plugin_context import get_plugin_with_context

API = "/api/ccp4i2"
BETA = Path(ccp4i2.__file__).parent / "demo_data" / "beta_blip" / "beta.seq"


@pytest.fixture
def client(bypass_api_permissions):
    return APIClient()


@pytest.fixture
def job(client, tmp_path):
    directory = tmp_path / "provideseq"
    (directory / "CCP4_JOBS").mkdir(parents=True)
    project = models.Project.objects.create(name="provideseq", directory=str(directory))
    r = client.post(f"{API}/projects/{project.id}/create_task/",
                    data={"task_name": "ProvideSequence"}, format="json")
    assert r.status_code == 200, r.content
    return models.Job.objects.get(id=r.json()["data"]["new_job"]["id"])


def _errors(client, job):
    r = client.get(f"{API}/jobs/{job.id}/validation/")
    assert r.status_code == 200, r.content
    return r.json()["data"]["xml"]


def _upload(client, job):
    r = client.post(f"{API}/jobs/{job.id}/upload_file_param/",
                    {"file": SimpleUploadedFile("beta.seq", BETA.read_bytes()),
                     "objectPath": "inputData.SEQIN"}, format="multipart")
    assert r.status_code == 200, r.content


def test_a_sequence_file_needs_no_text(client, job):
    assert "SEQUENCETEXT" in _errors(client, job)  # nothing given: required
    _upload(client, job)
    assert "SEQUENCETEXT" not in _errors(client, job)


def test_the_text_is_read_from_the_file_at_run_time(client, job):
    _upload(client, job)
    job.refresh_from_db()
    plugin = get_plugin_with_context(job).data
    text = plugin._text_from_files()
    assert text == BETA.read_text()


def test_whitespace_is_still_no_sequence(client, job):
    r = client.post(f"{API}/jobs/{job.id}/set_parameter/",
                    data={"object_path": "controlParameters.SEQUENCETEXT", "value": "\n"}, format="json")
    assert r.status_code == 200, r.content
    assert "SEQUENCETEXT" in _errors(client, job)


def _resolved(client, job, reference):
    r = client.get(f"{API}/projects/{job.project.id}/resolve_fileuse/", {"fileuse": reference})
    assert r.status_code == 200, r.content
    value = r.json()["data"]
    value.pop("fullPath", None)
    return value


def _new(client, project, task):
    r = client.post(f"{API}/projects/{project.id}/create_task/", data={"task_name": task}, format="json")
    assert r.status_code == 200, r.content
    return models.Job.objects.get(id=r.json()["data"]["new_job"]["id"])


def test_a_project_file_of_the_wrong_kind_is_refused(client, job):
    # A sequence file offered to an ASU-contents input, as an agent did to
    # mrbump_basic's ASUIN; the interface's file menu would not list it
    _upload(client, job)
    seq = _resolved(client, job, f"[{job.number}].SEQIN")
    asu_job = _new(client, job.project, "ProvideAsuContents")
    r = client.post(f"{API}/jobs/{asu_job.id}/set_parameter/",
                    data={"object_path": "inputData.ASUCONTENTIN", "value": seq}, format="json")
    assert r.status_code == 400, r.content
    error = r.json()["error"]
    assert "application/CCP4-asu-content" in error and "beta.seq" in error and "application/CCP4-seq" in error

    # The same file where a sequence file goes is taken
    other = _new(client, job.project, "ProvideSequence")
    r = client.post(f"{API}/jobs/{other.id}/set_parameter/",
                    data={"object_path": "inputData.SEQIN", "value": seq}, format="json")
    assert r.status_code == 200, r.content


def test_an_asu_entry_is_filled_from_its_source_file(client, job):
    # In the interface, choosing an entry's sequence file fills the entry in
    # the browser; an agent's upload left the sequence empty, and validation
    # then failed on it
    asu_job = _new(client, job.project, "ProvideAsuContents")
    r = client.post(f"{API}/jobs/{asu_job.id}/add_list_item/",
                    data={"object_path": "inputData.ASU_CONTENT", "values": {"nCopies": 2}}, format="json")
    assert r.status_code == 200, r.content
    item = r.json()["data"]["item"]
    r = client.post(f"{API}/jobs/{asu_job.id}/upload_file_param/",
                    {"file": SimpleUploadedFile("beta.seq", BETA.read_bytes()),
                     "objectPath": f"{item}.source"}, format="multipart")
    assert r.status_code == 200, r.content
    r = client.get(f"{API}/jobs/{asu_job.id}/parameters/", {"query": "ASU_CONTENT"})
    asu = next(p for p in r.json()["data"]["parameters"] if p["path"] == "inputData.ASU_CONTENT")
    entry = asu["value"][int(item.rsplit("[", 1)[1].rstrip("]"))]
    expected = "".join(BETA.read_text().splitlines()[1:]).replace(" ", "")
    assert entry["sequence"] == expected
    assert entry["name"] and entry["polymerType"] == "PROTEIN" and entry["nCopies"] == 2
    assert "sequence" not in _errors(client, asu_job).split(item)[-1][:200]


def test_a_kind_of_the_file_an_input_takes_is_taken():
    # An observed-data MTZ is an MTZ: ProvideAsuContents' HKLIN (a plain
    # MTZ) refused the import's OBSOUT on the exact-type rule alone
    from ccp4i2.core.cdata_registry import cdata_classes
    from ccp4i2.lib.utils.parameters.set_parameter import _is_a
    mtz = cdata_classes()["CMtzDataFile"]()
    assert _is_a("application/CCP4-mtz-observed", mtz)
    assert not _is_a("application/CCP4-seq", mtz)
    seq = cdata_classes()["CSeqDataFile"]()
    assert not _is_a("application/CCP4-mtz-observed", seq)
