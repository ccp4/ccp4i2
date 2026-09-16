"""Importing a project *zip* by path (#512 follow-up).

The zip-import endpoint reuses the upload-by-path gate: a staged/desktop
deployment may hand the server a path to a zip already on its disk instead of
pushing the bytes through the ingress cap. The same gate that protects
upload_file_param protects this route -- a path outside the staging dir is
refused, so a shared deployment cannot be steered at an arbitrary server file.

The real import runs detached, so both the archive inspection and the
management command are stubbed here: this test is about *which zip the endpoint
acts on*, not about unpacking one.
"""

import pytest
from rest_framework.test import APIClient

IMPORT_URL = "/api/ccp4i2/projects/import_project/"


@pytest.fixture
def client(bypass_api_permissions):
    return APIClient()


@pytest.fixture(autouse=True)
def _stub_import(monkeypatch):
    """Stub the archive inspection and the detached import command."""
    calls = []

    def fake_inspect(zip_path):
        return {"project_name": "staged_proj", "jobs": []}

    def fake_call_command(name, *args, **kwargs):
        calls.append((name, args, kwargs))

    monkeypatch.setattr(
        "ccp4i2.db.import_i2xml.inspect_ccp4_project_zip", fake_inspect
    )
    monkeypatch.setattr(
        "ccp4i2.api.ProjectViewSet.call_command", fake_call_command
    )
    return calls


def _staging(tmp_path, monkeypatch, name="project.ccp4_project.zip"):
    staging = tmp_path / "staging"
    staging.mkdir()
    zip_path = staging / name
    zip_path.write_bytes(b"PK\x03\x04staged")
    monkeypatch.setenv("CCP4I2_IMPORT_STAGING_DIR", str(staging))
    return staging, zip_path


def test_staged_zip_imports_by_path_without_copying(
    client, tmp_path, monkeypatch, _stub_import
):
    _, zip_path = _staging(tmp_path, monkeypatch)

    resp = client.post(IMPORT_URL, data={"local_path": str(zip_path)})

    assert resp.status_code == 200, resp.content
    # The import command ran against the staged path itself -- had the endpoint
    # copied it into the upload store first, the argument would be a
    # secure_storage path under MEDIA_ROOT, not the staged zip.
    assert _stub_import == [
        ("import_ccp4_project_zip", (str(zip_path.resolve()), "--detach"), {})
    ]


def test_path_outside_staging_is_refused_and_falls_back(
    client, tmp_path, monkeypatch, _stub_import
):
    # A sibling of the staging dir: the gate returns None, so local_path is
    # ignored and the endpoint looks for a body upload -- of which there is none.
    _staging(tmp_path, monkeypatch)
    outside = tmp_path / "secret.ccp4_project.zip"
    outside.write_bytes(b"PK\x03\x04secret")

    resp = client.post(IMPORT_URL, data={"local_path": str(outside)})

    assert resp.status_code == 400
    assert _stub_import == []  # nothing imported


def test_no_mode_ignores_local_path(client, tmp_path, monkeypatch, _stub_import):
    # With neither env signal set (the web default) a local_path is never
    # trusted, even for a real file -> fall through to the (absent) body.
    monkeypatch.delenv("CCP4I2_IMPORT_STAGING_DIR", raising=False)
    monkeypatch.delenv("CCP4I2_LOCAL_SESSION_TOKEN", raising=False)
    real = tmp_path / "real.ccp4_project.zip"
    real.write_bytes(b"PK\x03\x04real")

    resp = client.post(IMPORT_URL, data={"local_path": str(real)})

    assert resp.status_code == 400
    assert _stub_import == []
