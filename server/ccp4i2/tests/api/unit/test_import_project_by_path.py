"""Importing a project *zip* on a served deployment, via a staged handle.

A served deployment imports a large zip by an owner-bound staged handle
(``staged_upload=<uuid>``), never a client-named path. This checks the endpoint
resolves the handle to the staged zip, dispatches the (stubbed) detached import
against it, and that a client-named ``local_path`` is no longer honoured in cloud
mode. The real import runs detached, so inspection and the command are stubbed:
this is about *which zip the endpoint acts on*.
"""

import pytest
from rest_framework.test import APIClient

from ccp4i2.db import models

IMPORT_URL = "/api/ccp4i2/projects/import_project/"


@pytest.fixture
def client(bypass_api_permissions):
    return APIClient()


@pytest.fixture(autouse=True)
def _stub_import(monkeypatch):
    calls = []
    monkeypatch.setattr(
        "ccp4i2.db.import_i2xml.inspect_ccp4_project_zip",
        lambda zip_path: {"project_name": "staged_proj", "jobs": []},
    )
    monkeypatch.setattr(
        "ccp4i2.api.ProjectViewSet.call_command",
        lambda name, *a, **k: calls.append((name, a, k)),
    )
    return calls


def _staging(tmp_path, monkeypatch):
    staging = tmp_path / "staging"
    staging.mkdir()
    monkeypatch.setenv("CCP4I2_IMPORT_STAGING_DIR", str(staging))
    return staging


def _ready_handle(staging, owner="None", name="project.ccp4_project.zip"):
    """A ready StagedUpload row + its assembled file (owner is str(user.pk),
    which is 'None' for the anonymous test client)."""
    row = models.StagedUpload.objects.create(
        owner=owner, filename=name, size_bytes=13,
        state=models.StagedUpload.State.READY,
    )
    d = staging / str(row.uuid)
    d.mkdir()
    zip_path = d / name
    zip_path.write_bytes(b"PK\x03\x04staged")
    return row, zip_path


def test_staged_handle_imports_the_staged_zip(client, tmp_path, monkeypatch, _stub_import):
    staging = _staging(tmp_path, monkeypatch)
    row, zip_path = _ready_handle(staging)

    resp = client.post(IMPORT_URL, data={"staged_upload": str(row.uuid)})

    assert resp.status_code == 200, resp.content
    # Imported straight from the staged path -- no copy into MEDIA_ROOT.
    assert _stub_import == [
        ("import_ccp4_project_zip", (str(zip_path), "--detach"), {})
    ]
    # The handle is retired so it can't be reused (file kept for the detached importer).
    row.refresh_from_db()
    assert row.state == models.StagedUpload.State.CONSUMED


def test_foreign_handle_is_404(client, tmp_path, monkeypatch, _stub_import):
    staging = _staging(tmp_path, monkeypatch)
    row, _ = _ready_handle(staging, owner="someone-else")

    resp = client.post(IMPORT_URL, data={"staged_upload": str(row.uuid)})

    assert resp.status_code == 404
    assert _stub_import == []


def test_cloud_local_path_is_not_honoured(client, tmp_path, monkeypatch, _stub_import):
    # In cloud mode a client-named local_path is dead -- even inside the staging
    # dir. It falls through to the (absent) body -> 400.
    staging = _staging(tmp_path, monkeypatch)
    inside = staging / "sneaky.ccp4_project.zip"
    inside.write_bytes(b"PK\x03\x04x")

    resp = client.post(IMPORT_URL, data={"local_path": str(inside)})

    assert resp.status_code == 400
    assert _stub_import == []


def test_no_mode_ignores_local_path(client, tmp_path, monkeypatch, _stub_import):
    monkeypatch.delenv("CCP4I2_IMPORT_STAGING_DIR", raising=False)
    monkeypatch.delenv("CCP4I2_LOCAL_SESSION_TOKEN", raising=False)
    real = tmp_path / "real.ccp4_project.zip"
    real.write_bytes(b"PK\x03\x04real")

    resp = client.post(IMPORT_URL, data={"local_path": str(real)})

    assert resp.status_code == 400
    assert _stub_import == []
