"""The staging transport endpoints (begin / chunk / finish / status) over HTTP.

Exercises the thin viewset: routing, the raw-bytes chunk parser, the status-code
contract, and that a finished handle is ready to import. The lifecycle logic
itself is covered at the library level in tests/db/test_staged_upload.py.
"""

import hashlib

import pytest
from rest_framework.permissions import AllowAny
from rest_framework.test import APIClient

from ccp4i2.api.StagedUploadViewSet import StagedUploadViewSet

BASE = "/api/ccp4i2/staged-uploads/"


@pytest.fixture
def client(monkeypatch, tmp_path, settings):
    monkeypatch.setattr(StagedUploadViewSet, "permission_classes", [AllowAny])
    d = tmp_path / "staging"
    d.mkdir()
    monkeypatch.setenv("CCP4I2_IMPORT_STAGING_DIR", str(d))
    settings.CCP4I2_IMPORT_STAGING_CHUNK_BYTES = 8
    settings.CCP4I2_IMPORT_STAGING_MAX_BYTES = 1000
    return APIClient()


def _put_chunk(client, uid, index, data):
    return client.put(f"{BASE}{uid}/chunks/{index}/", data=data,
                      content_type="application/octet-stream")


def test_begin_chunk_finish_then_ready(client):
    data = b"a-staged-cryo-map!!!"        # 20 bytes -> 3 chunks (8,8,4)
    h = hashlib.sha256(data).hexdigest()

    begin = client.post(BASE, data={"filename": "m.mrc", "size_bytes": len(data),
                                    "sha256": h})
    assert begin.status_code == 201, begin.content
    uid = begin.json()["upload_id"]
    assert begin.json()["chunk_bytes"] == 8

    for i, off in enumerate(range(0, len(data), 8)):
        assert _put_chunk(client, uid, i, data[off:off + 8]).status_code == 204

    status = client.get(f"{BASE}{uid}/").json()
    assert status["present_chunks"] == [0, 1, 2]
    assert status["state"] == "staging"

    fin = client.post(f"{BASE}{uid}/finish/")
    assert fin.status_code == 200, fin.content
    assert fin.json()["state"] == "ready"


def test_over_max_is_413_at_begin(client):
    resp = client.post(BASE, data={"filename": "big.mrc", "size_bytes": 99999})
    assert resp.status_code == 413


def test_unknown_handle_is_404(client):
    resp = client.get(f"{BASE}00000000-0000-0000-0000-000000000000/")
    assert resp.status_code == 404


def test_finish_with_missing_chunk_is_409(client):
    begin = client.post(BASE, data={"filename": "m.mrc", "size_bytes": 20})
    uid = begin.json()["upload_id"]
    _put_chunk(client, uid, 0, b"eightby.")
    _put_chunk(client, uid, 2, b"xyz")        # index 1 missing -> non-contiguous
    resp = client.post(f"{BASE}{uid}/finish/")
    assert resp.status_code == 409
