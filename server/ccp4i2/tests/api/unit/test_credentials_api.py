"""API tests for the credential endpoints.

The behaviour that matters here is not the happy path — it is that secrets go
*in* and never come *out*, and that a shared deployment cannot be talked into
accepting one.
"""

import json

import pytest
from rest_framework.test import APIClient

from ccp4i2.config import credentials as cred

VALUES = {"token_id": "tokenid1234", "token_secret": "supersecret"}


@pytest.fixture
def client(bypass_api_permissions):
    return APIClient()


@pytest.fixture
def desktop(tmp_path, monkeypatch):
    """An isolated desktop-shaped environment with no real keychain."""
    monkeypatch.setenv("CCP4I2_HOME", str(tmp_path))
    monkeypatch.setenv("CCP4I2_LOCAL_SESSION_TOKEN", "desktop")
    monkeypatch.setattr(cred, "_keyring", lambda: None)
    for field in cred.PDB_REDO.fields:
        monkeypatch.delenv(field.env, raising=False)
        monkeypatch.delenv(
            cred._session_env_key(cred._storage_key("pdb_redo"), field.name),
            raising=False,
        )
    cred._validation_cache.clear()
    return tmp_path


def _post(client, url, payload):
    return client.post(url, data=json.dumps(payload), content_type="application/json")


def test_list_describes_without_disclosing(client, desktop):
    _post(client, "/api/ccp4i2/config/credentials/pdb_redo/set/", {"values": VALUES})

    resp = client.get("/api/ccp4i2/config/credentials/")
    assert resp.status_code == 200
    body = resp.json()["data"]
    assert body["editable"] is True

    # The whole response, serialised, must not contain the secret — this is the
    # invariant the endpoint exists to hold.
    blob = json.dumps(body)
    assert VALUES["token_secret"] not in blob
    assert VALUES["token_id"] not in blob

    descriptor = {c["name"]: c for c in body["credentials"]}["pdb_redo"]
    assert descriptor["isSet"] is True
    assert descriptor["hint"] == "1234"
    assert [f["name"] for f in descriptor["fields"]] == ["token_id", "token_secret"]


def test_set_and_clear_roundtrip(client, desktop):
    resp = _post(
        client, "/api/ccp4i2/config/credentials/pdb_redo/set/", {"values": VALUES}
    )
    assert resp.status_code == 200
    data = resp.json()["data"]
    assert data["storedIn"] == "file"  # keychain stubbed out in this fixture
    assert data["credential"]["isSet"] is True
    assert cred.get_credential("pdb_redo") == VALUES

    resp = _post(client, "/api/ccp4i2/config/credentials/pdb_redo/clear/", {})
    assert resp.status_code == 200
    assert resp.json()["data"]["credential"]["isSet"] is False
    assert cred.get_credential("pdb_redo") is None


def test_set_rejects_unknown_credential(client, desktop):
    resp = _post(
        client, "/api/ccp4i2/config/credentials/nope/set/", {"values": VALUES}
    )
    assert resp.status_code == 404


def test_set_rejects_empty_values(client, desktop):
    resp = _post(client, "/api/ccp4i2/config/credentials/pdb_redo/set/", {"values": {}})
    assert resp.status_code == 400


def test_set_rejects_unknown_persistence(client, desktop):
    resp = _post(
        client,
        "/api/ccp4i2/config/credentials/pdb_redo/set/",
        {"values": VALUES, "persistence": "forever"},
    )
    assert resp.status_code == 400


def test_session_persistence_leaves_no_file(client, desktop):
    resp = _post(
        client,
        "/api/ccp4i2/config/credentials/pdb_redo/set/",
        {"values": VALUES, "persistence": "session"},
    )
    assert resp.json()["data"]["storedIn"] == "session"
    assert not (desktop / cred.CREDENTIALS_FILENAME).exists()
    assert cred.get_credential("pdb_redo") == VALUES


def test_validate_stored_credential(client, desktop, monkeypatch):
    import dataclasses

    monkeypatch.setitem(
        cred.CREDENTIALS,
        "pdb_redo",
        dataclasses.replace(
            cred.PDB_REDO, validator=lambda values: (True, "Token accepted.")
        ),
    )
    _post(client, "/api/ccp4i2/config/credentials/pdb_redo/set/", {"values": VALUES})

    resp = _post(client, "/api/ccp4i2/config/credentials/pdb_redo/validate/", {})
    assert resp.status_code == 200
    assert resp.json()["data"] == {"ok": True, "message": "Token accepted."}

    # The outcome is cached for display, so the panel can say when it last
    # checked without re-probing on every render.
    resp = client.get("/api/ccp4i2/config/credentials/")
    descriptor = {
        c["name"]: c for c in resp.json()["data"]["credentials"]
    }["pdb_redo"]
    assert descriptor["lastValidationOk"] is True
    assert descriptor["lastValidated"] is not None


def test_validate_reports_failure_without_echoing_the_secret(
    client, desktop, monkeypatch
):
    import dataclasses

    monkeypatch.setitem(
        cred.CREDENTIALS,
        "pdb_redo",
        dataclasses.replace(
            cred.PDB_REDO, validator=lambda values: (False, "PDB-REDO rejected this token.")
        ),
    )
    resp = _post(
        client,
        "/api/ccp4i2/config/credentials/pdb_redo/validate/",
        {"values": VALUES},
    )
    body = resp.json()["data"]
    assert body["ok"] is False
    assert VALUES["token_secret"] not in json.dumps(body)


# ---------------------------------------------------------------------------
# The gate
# ---------------------------------------------------------------------------


@pytest.fixture
def cloud(tmp_path, monkeypatch):
    """A shared/server-shaped deployment: no desktop token, no opt-in."""
    monkeypatch.setenv("CCP4I2_HOME", str(tmp_path))
    monkeypatch.delenv("CCP4I2_LOCAL_SESSION_TOKEN", raising=False)
    monkeypatch.delenv("CCP4I2_ALLOW_CREDENTIAL_WRITE", raising=False)
    monkeypatch.setattr(cred, "_keyring", lambda: None)
    return tmp_path


def test_cloud_cannot_set_a_credential(client, cloud):
    """A loopback REMOTE_ADDR (which the test client sends, and which a reverse
    proxy also produces) must NOT be enough to accept a secret."""
    resp = _post(
        client, "/api/ccp4i2/config/credentials/pdb_redo/set/", {"values": VALUES}
    )
    assert resp.status_code == 409
    assert resp.json()["success"] is False
    assert cred.get_credential("pdb_redo") is None


def test_cloud_cannot_clear_a_credential(client, cloud):
    resp = _post(client, "/api/ccp4i2/config/credentials/pdb_redo/clear/", {})
    assert resp.status_code == 409


def test_cloud_cannot_test_arbitrary_values(client, cloud):
    """Testing supplied values means accepting a secret over the wire."""
    resp = _post(
        client,
        "/api/ccp4i2/config/credentials/pdb_redo/validate/",
        {"values": VALUES},
    )
    assert resp.status_code == 409


def test_cloud_list_is_read_only(client, cloud):
    resp = client.get("/api/ccp4i2/config/credentials/")
    assert resp.status_code == 200
    assert resp.json()["data"]["editable"] is False


def test_explicit_opt_in_re_enables_writes(client, cloud, monkeypatch):
    """A trusted single-user server deployment can opt back in."""
    monkeypatch.setenv("CCP4I2_ALLOW_CREDENTIAL_WRITE", "1")
    resp = _post(
        client, "/api/ccp4i2/config/credentials/pdb_redo/set/", {"values": VALUES}
    )
    assert resp.status_code == 200
    assert cred.get_credential("pdb_redo") == VALUES
