"""version_info advertises the staging transport only when it is enabled."""

import pytest
from rest_framework.test import APIClient

VERSION_URL = "/api/ccp4i2/version/"


@pytest.fixture
def client():
    return APIClient()


def test_capability_present_when_staging_dir_set(client, tmp_path, monkeypatch, settings):
    d = tmp_path / "staging"
    d.mkdir()
    monkeypatch.setenv("CCP4I2_IMPORT_STAGING_DIR", str(d))
    settings.CCP4I2_IMPORT_STAGING_CHUNK_BYTES = 16 * 1024 * 1024
    settings.CCP4I2_IMPORT_STAGING_MAX_BYTES = 2 * 1024 * 1024 * 1024
    settings.CCP4I2_IMPORT_STAGING_THRESHOLD_BYTES = 32 * 1024 * 1024

    body = client.get(VERSION_URL).json()
    assert body["import_staging"] == {
        "chunk_bytes": 16 * 1024 * 1024,
        "max_bytes": 2 * 1024 * 1024 * 1024,
        "threshold_bytes": 32 * 1024 * 1024,
    }


def test_capability_absent_when_not_configured(client, monkeypatch):
    monkeypatch.delenv("CCP4I2_IMPORT_STAGING_DIR", raising=False)
    body = client.get(VERSION_URL).json()
    assert "import_staging" not in body
