"""Importing a project directory must not touch the original.

Restoring from a directory re-roots the snapshot onto wherever the directory is
and, when repairing paths, rewrites the absolute paths inside the project's own
files. That is right for RECOVERY (your projects, your database) and wrong for
IMPORT (a folder someone handed you, or one your existing CCP4i2 is still
using). So import copies first and adopts the copy.
"""

import json
from pathlib import Path

import pytest
from django.test import override_settings
from rest_framework.test import APIClient

from ccp4i2.db.restore_project import copy_project_directory


@pytest.fixture
def client(bypass_api_permissions, monkeypatch):
    """A desktop-shaped client.

    The restore endpoint refuses a server-side path with 409 unless this is a
    desktop launch — so without this the endpoint tests below never reach the
    code they are about, and pass for the wrong reason. They did exactly that
    until the refusal test noticed a 409 where it wanted a 403.
    """
    monkeypatch.setenv("CCP4I2_LOCAL_SESSION_TOKEN", "desktop")
    return APIClient()


def _post(client, payload):
    return client.post(
        "/api/ccp4i2/projects/restore/",
        data=json.dumps(payload),
        content_type="application/json",
    )


def _fake_project(root: Path, name: str) -> Path:
    directory = root / name
    (directory / "CCP4_JOBS" / "job_1").mkdir(parents=True)
    (directory / "CCP4_JOBS" / "job_1" / "params.xml").write_text(
        f"<params><path>{directory}/CCP4_JOBS/job_1</path></params>"
    )
    return directory


# ---------------------------------------------------------------------------
# the copy helper itself
# ---------------------------------------------------------------------------

def test_copy_leaves_the_original_byte_for_byte(tmp_path):
    source = _fake_project(tmp_path / "elsewhere", "gamma")
    original = (source / "CCP4_JOBS" / "job_1" / "params.xml").read_text()

    copy = copy_project_directory(source, tmp_path / "store")

    assert copy == tmp_path / "store" / "gamma"
    assert (copy / "CCP4_JOBS" / "job_1" / "params.xml").is_file()
    assert (source / "CCP4_JOBS" / "job_1" / "params.xml").read_text() == original


def test_copy_refuses_an_occupied_destination(tmp_path):
    """Merging two same-named projects would interleave job directories, and
    the damage would surface much later as missing jobs."""
    source = _fake_project(tmp_path / "elsewhere", "gamma")
    (tmp_path / "store" / "gamma").mkdir(parents=True)

    with pytest.raises(FileExistsError):
        copy_project_directory(source, tmp_path / "store")


def test_copy_refuses_a_destination_inside_the_source(tmp_path):
    source = _fake_project(tmp_path / "elsewhere", "gamma")
    with pytest.raises(ValueError):
        copy_project_directory(source, source / "inner")


def test_copy_reports_a_missing_source(tmp_path):
    with pytest.raises(FileNotFoundError):
        copy_project_directory(tmp_path / "nope", tmp_path / "store")


# ---------------------------------------------------------------------------
# the endpoint's policy
# ---------------------------------------------------------------------------

def test_in_place_import_is_refused_by_default(client, tmp_path):
    source = _fake_project(tmp_path / "elsewhere", "gamma")
    resp = _post(client, {"source": "directory", "path": str(source), "copy": "false"})
    assert resp.status_code == 403
    assert "in place" in resp.json().get("error", resp.json().get("message", "")).lower()


@override_settings(CCP4I2_ALLOW_INPLACE_MIGRATION=True)
def test_the_flag_permits_in_place_import(client, tmp_path):
    """Kept for whoever is making in-place safe; not reachable from the UI."""
    source = _fake_project(tmp_path / "elsewhere", "gamma")
    resp = _post(client, {"source": "directory", "path": str(source), "copy": "false"})
    assert resp.status_code != 403


def test_a_dry_run_copies_nothing(client, tmp_path, settings):
    source = _fake_project(tmp_path / "elsewhere", "gamma")
    _post(client, {"source": "directory", "path": str(source), "dry_run": "true"})
    assert not (Path(settings.CCP4I2_PROJECTS_DIR) / "gamma").exists()
