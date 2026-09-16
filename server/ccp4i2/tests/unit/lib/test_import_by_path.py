"""The import-by-path gate (#512): which paths a deployment may import by copy.

Security-critical: the cloud/staged mode must import ONLY files inside the
configured staging directory, never an arbitrary server path -- otherwise a
shared deployment becomes an arbitrary-file-read. Desktop mode trusts the local
user. With neither signal set (the web default), no path is trusted.

CCP4-free: pure env + pathlib logic (imports gemmi transitively via upload_param).
"""

import pathlib

import pytest

from ccp4i2.lib.utils.files.upload_param import resolve_importable_path

LOCAL = "CCP4I2_LOCAL_SESSION_TOKEN"
STAGING = "CCP4I2_IMPORT_STAGING_DIR"


@pytest.fixture(autouse=True)
def _clear_env(monkeypatch):
    monkeypatch.delenv(LOCAL, raising=False)
    monkeypatch.delenv(STAGING, raising=False)


def _make(tmp_path, name="map.mrc"):
    p = tmp_path / name
    p.write_bytes(b"data")
    return p


def test_no_mode_never_imports_by_path(tmp_path):
    # The web default: even a real path is not trusted -> fall back to the body.
    assert resolve_importable_path(str(_make(tmp_path))) is None


def test_empty_path_is_none(monkeypatch):
    monkeypatch.setenv(LOCAL, "tok")
    assert resolve_importable_path("") is None
    assert resolve_importable_path(None) is None


def test_desktop_local_allows_any_readable_file(tmp_path, monkeypatch):
    monkeypatch.setenv(LOCAL, "tok")
    f = _make(tmp_path)
    assert resolve_importable_path(str(f)) == f


def test_desktop_local_rejects_missing_file(tmp_path, monkeypatch):
    monkeypatch.setenv(LOCAL, "tok")
    assert resolve_importable_path(str(tmp_path / "nope.mrc")) is None


def test_served_never_honours_local_path(tmp_path, monkeypatch):
    # A served deployment (staging dir set, no desktop token) does NOT trust a
    # client-named local_path -- even one inside the staging directory. Cloud
    # imports go through an owner-bound staged handle instead (StagedUpload), so
    # the "any file in the staging dir is importable by anyone who can name it"
    # hole is closed. resolve_importable_path returns None for every path here.
    staging = tmp_path / "staging"
    staging.mkdir()
    inside = _make(staging)
    outside = _make(tmp_path, "outside.mrc")
    monkeypatch.setenv(STAGING, str(staging))
    assert resolve_importable_path(str(inside)) is None
    assert resolve_importable_path(str(outside)) is None
    assert resolve_importable_path(str(staging / ".." / "outside.mrc")) is None
