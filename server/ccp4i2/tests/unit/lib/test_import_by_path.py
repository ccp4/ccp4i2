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


def test_staged_allows_inside_staging_dir(tmp_path, monkeypatch):
    staging = tmp_path / "staging"
    staging.mkdir()
    f = _make(staging)
    monkeypatch.setenv(STAGING, str(staging))
    assert resolve_importable_path(str(f)) == f.resolve()


def test_staged_refuses_outside_staging_dir(tmp_path, monkeypatch):
    staging = tmp_path / "staging"
    staging.mkdir()
    outside = _make(tmp_path, "outside.mrc")   # sibling of staging, not inside
    monkeypatch.setenv(STAGING, str(staging))
    assert resolve_importable_path(str(outside)) is None


def test_staged_refuses_symlink_escape(tmp_path, monkeypatch):
    # A symlink inside staging pointing outside must not smuggle an arbitrary
    # read: resolve() follows it, and the is_relative_to check then fails.
    staging = tmp_path / "staging"
    staging.mkdir()
    secret = _make(tmp_path, "secret.mrc")
    link = staging / "innocent.mrc"
    try:
        link.symlink_to(secret)
    except OSError:
        pytest.skip("symlinks not supported here")
    monkeypatch.setenv(STAGING, str(staging))
    assert resolve_importable_path(str(link)) is None


def test_staged_refuses_traversal(tmp_path, monkeypatch):
    staging = tmp_path / "staging"
    staging.mkdir()
    _make(tmp_path, "secret.mrc")
    monkeypatch.setenv(STAGING, str(staging))
    assert resolve_importable_path(str(staging / ".." / "secret.mrc")) is None
