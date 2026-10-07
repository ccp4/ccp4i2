"""Importing a project zip that clashes with an existing project.

The import runs detached, and the importer cannot create a project whose name or
directory another project holds, so without a check before dispatch the import
failed in a subprocess while the upload reported success.
"""

import uuid

import pytest
from django.conf import settings
from django.core.files.uploadedfile import SimpleUploadedFile
from rest_framework.test import APIClient

from ccp4i2.db import models

IMPORT_URL = "/api/ccp4i2/projects/import_project/"


@pytest.fixture
def client(bypass_api_permissions, settings, tmp_path):
    settings.MEDIA_ROOT = tmp_path
    return APIClient()


@pytest.fixture
def existing(db):
    return models.Project.objects.create(
        name="lysozyme", directory=str(settings.CCP4I2_PROJECTS_DIR / "lysozyme")
    )


@pytest.fixture
def dispatched(monkeypatch):
    calls = []
    monkeypatch.setattr(
        "ccp4i2.api.ProjectViewSet.call_command",
        lambda name, *a, **k: calls.append((name, a, k)),
    )
    return calls


def _archive_holds(monkeypatch, **summary):
    summary = {"project_uuid": uuid.uuid4().hex, "jobs": 3, **summary}
    monkeypatch.setattr(
        "ccp4i2.db.import_i2xml.inspect_ccp4_project_zip", lambda zip_path: summary
    )


def _upload(client, *names):
    files = [SimpleUploadedFile(name, b"PK\x03\x04") for name in names]
    return client.post(IMPORT_URL, {"files": files}, format="multipart")


def test_a_taken_name_is_rejected(client, existing, dispatched, monkeypatch):
    _archive_holds(
        monkeypatch, project_name="lysozyme", recorded_directory="/elsewhere/lyso2"
    )

    response = _upload(client, "lysozyme.ccp4_project.zip")

    assert response.status_code == 409
    assert "A project named 'lysozyme' already exists" in response.json()["error"]
    assert dispatched == []


def test_a_taken_directory_is_rejected(client, existing, dispatched, monkeypatch):
    _archive_holds(
        monkeypatch, project_name="lysozyme_b", recorded_directory="/elsewhere/lysozyme"
    )

    response = _upload(client, "lysozyme_b.ccp4_project.zip")

    assert response.status_code == 409
    assert "already uses the directory" in response.json()["error"]
    assert dispatched == []


def test_reimporting_the_same_project_is_allowed(
    client, existing, dispatched, monkeypatch
):
    _archive_holds(
        monkeypatch,
        project_name="lysozyme",
        project_uuid=existing.uuid.hex,
        recorded_directory="/elsewhere/lysozyme",
    )

    response = _upload(client, "lysozyme.ccp4_project.zip")

    assert response.status_code == 200, response.content
    assert len(dispatched) == 1


def test_one_clash_dispatches_none_of_the_request(
    client, existing, dispatched, monkeypatch
):
    summaries = iter(
        [
            {"project_name": "thaumatin", "recorded_directory": "/x/thaumatin"},
            {"project_name": "lysozyme", "recorded_directory": "/x/lyso2"},
        ]
    )
    monkeypatch.setattr(
        "ccp4i2.db.import_i2xml.inspect_ccp4_project_zip",
        lambda zip_path: {"project_uuid": uuid.uuid4().hex, "jobs": 1, **next(summaries)},
    )

    response = _upload(client, "thaumatin.zip", "lysozyme.zip")

    assert response.status_code == 409
    assert dispatched == []
