"""Unit tests for the credential store (``ccp4i2.config.credentials``).

No database, no network, no CCP4 binaries. The live-service validator is
stubbed; what is tested here is the store's contract: resolution precedence,
the three persistence modes, and — most importantly — that no secret material
escapes through a descriptor.
"""

import json
import os
import stat

import pytest

from ccp4i2.config import credentials as cred
from ccp4i2.config import preferences as prefs


@pytest.fixture
def home(tmp_path, monkeypatch):
    """An isolated CCP4i2 home, with no keyring and no env overrides."""
    monkeypatch.setenv("CCP4I2_HOME", str(tmp_path))
    # Force the file backend so these tests never touch the developer's real
    # keychain (and so they behave identically on CI).
    monkeypatch.setattr(cred, "_keyring", lambda: None)
    for field in cred.PDB_REDO.fields:
        monkeypatch.delenv(field.env, raising=False)
        monkeypatch.delenv(
            cred._session_env_key(cred._storage_key("pdb_redo"), field.name),
            raising=False,
        )
    cred._validation_cache.clear()
    return tmp_path


VALUES = {"token_id": "tokenid1234", "token_secret": "supersecret"}


# ---------------------------------------------------------------------------
# Resolution
# ---------------------------------------------------------------------------


def test_unset_credential_resolves_to_none(home):
    """The bug this store exists to fix: no token must mean None, not "None"."""
    assert cred.get_credential("pdb_redo") is None
    assert cred.resolve_credential("pdb_redo").source == "unset"


def test_partial_credential_is_not_usable(home):
    """Half a token pair is worse than none — it produces an opaque 401."""
    cred.set_credential("pdb_redo", {"token_id": "only-the-id"})
    assert cred.get_credential("pdb_redo") is None
    resolved = cred.resolve_credential("pdb_redo")
    assert resolved.values["token_id"] == "only-the-id"
    assert resolved.complete is False


def test_roundtrip_through_file_backend(home):
    assert cred.set_credential("pdb_redo", VALUES) == "file"
    assert cred.get_credential("pdb_redo") == VALUES
    assert cred.resolve_credential("pdb_redo").source == "file"


def test_file_backend_is_owner_only(home):
    cred.set_credential("pdb_redo", VALUES)
    path = home / cred.CREDENTIALS_FILENAME
    assert stat.S_IMODE(path.stat().st_mode) == 0o600


def test_environment_wins_over_stored(home, monkeypatch):
    """Cloud invariance: env-supplied values must beat anything on disk."""
    cred.set_credential("pdb_redo", VALUES)
    monkeypatch.setenv("PDB_REDO_TOKEN_ID", "envid")
    monkeypatch.setenv("PDB_REDO_TOKEN_SECRET", "envsecret")
    resolved = cred.resolve_credential("pdb_redo")
    assert resolved.source == "environment"
    assert resolved.values["token_id"] == "envid"


def test_session_wins_over_stored_but_not_environment(home, monkeypatch):
    cred.set_credential("pdb_redo", VALUES)
    cred.set_credential(
        "pdb_redo",
        {"token_id": "sessid", "token_secret": "sesssecret"},
        persistence=cred.PERSIST_SESSION,
    )
    assert cred.resolve_credential("pdb_redo").source == "session"
    monkeypatch.setenv("PDB_REDO_TOKEN_ID", "envid")
    monkeypatch.setenv("PDB_REDO_TOKEN_SECRET", "envsecret")
    assert cred.resolve_credential("pdb_redo").source == "environment"


def test_session_persistence_never_touches_disk(home):
    cred.set_credential("pdb_redo", VALUES, persistence=cred.PERSIST_SESSION)
    assert cred.get_credential("pdb_redo") == VALUES
    assert not (home / cred.CREDENTIALS_FILENAME).exists()


def test_persistence_none_stores_nothing(home):
    assert cred.set_credential("pdb_redo", VALUES, persistence="none") == "none"
    assert cred.get_credential("pdb_redo") is None


def test_legacy_preferences_are_read_then_migrated(home):
    """A hand-set token in preferences.json keeps working, then moves out.

    Leaving a plaintext copy behind in the file users are asked to share would
    defeat the point of the store.
    """
    prefs.save_preferences(
        {
            "userPreferences": {
                "PDB_REDO_TOKEN_ID": "legacyid",
                "PDB_REDO_TOKEN_SECRET": "legacysecret",
            }
        }
    )
    resolved = cred.resolve_credential("pdb_redo")
    assert resolved.source == "preferences"
    assert resolved.values["token_id"] == "legacyid"

    cred.set_credential("pdb_redo", VALUES)
    assert prefs.load_preferences().get("userPreferences") == {}


def test_clear_removes_every_writable_layer(home):
    cred.set_credential("pdb_redo", VALUES)
    cred.set_credential("pdb_redo", VALUES, persistence=cred.PERSIST_SESSION)
    cred.clear_credential("pdb_redo")
    assert cred.get_credential("pdb_redo") is None


def test_scopes_are_independent(home):
    """Project-scoped tokens must not collide with the global one."""
    cred.set_credential("pdb_redo", VALUES)
    cred.set_credential(
        "pdb_redo",
        {"token_id": "projid", "token_secret": "projsecret"},
        scope=cred.SCOPE_PROJECT,
        scope_id="abc-123",
    )
    assert cred.get_credential("pdb_redo")["token_id"] == "tokenid1234"
    project = cred.get_credential(
        "pdb_redo", scope=cred.SCOPE_PROJECT, scope_id="abc-123"
    )
    assert project["token_id"] == "projid"


def test_unknown_credential(home):
    assert cred.get_spec("nope") is None
    assert cred.get_credential("nope") is None
    with pytest.raises(KeyError):
        cred.set_credential("nope", VALUES)


# ---------------------------------------------------------------------------
# Descriptors — the load-bearing invariant
# ---------------------------------------------------------------------------


def test_descriptor_never_contains_secret_material(home):
    cred.set_credential("pdb_redo", VALUES)
    blob = json.dumps(cred.describe_credential("pdb_redo"))
    assert VALUES["token_secret"] not in blob
    # Not even the non-secret id is echoed in full — only a 4-character hint.
    assert VALUES["token_id"] not in blob


def test_descriptor_reports_state(home):
    descriptor = cred.describe_credential("pdb_redo")
    assert descriptor["isSet"] is False
    assert descriptor["secure"] is False  # keyring stubbed out
    assert [f["name"] for f in descriptor["fields"]] == ["token_id", "token_secret"]
    assert descriptor["canValidate"] is True

    cred.set_credential("pdb_redo", VALUES)
    descriptor = cred.describe_credential("pdb_redo")
    assert descriptor["isSet"] is True
    assert descriptor["source"] == "file"
    assert descriptor["hint"] == "1234"


def test_env_supplied_credential_is_not_editable(home, monkeypatch):
    """In cloud the token belongs to the deployment, not the user."""
    monkeypatch.setenv("PDB_REDO_TOKEN_ID", "envid")
    monkeypatch.setenv("PDB_REDO_TOKEN_SECRET", "envsecret")
    assert cred.describe_credential("pdb_redo", editable=True)["editable"] is False


def test_storage_label_is_honest_about_the_fallback(home):
    """If we are not getting OS protection the user must be told so."""
    label = cred.storage_label()
    assert "no secret store" in label
    assert cred.CREDENTIALS_FILENAME in label


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def test_validate_without_a_credential(home):
    ok, message = cred.validate_credential("pdb_redo")
    assert ok is False
    assert "No credential" in message


def _with_validator(monkeypatch, validator):
    """Swap in a stub validator (CredentialSpec is frozen, so replace the spec)."""
    import dataclasses

    spec = dataclasses.replace(cred.PDB_REDO, validator=validator)
    monkeypatch.setitem(cred.CREDENTIALS, "pdb_redo", spec)


def test_validate_uses_the_spec_validator(home, monkeypatch):
    seen = {}

    def fake(values):
        seen.update(values)
        return True, "fine"

    _with_validator(monkeypatch, fake)
    cred.set_credential("pdb_redo", VALUES)
    ok, message = cred.validate_credential("pdb_redo")
    assert (ok, message) == (True, "fine")
    assert seen == VALUES


def test_validator_exceptions_do_not_escape(home, monkeypatch):
    def boom(values):
        raise RuntimeError("network on fire")

    _with_validator(monkeypatch, boom)
    cred.set_credential("pdb_redo", VALUES)
    ok, message = cred.validate_credential("pdb_redo")
    assert ok is False
    assert "network on fire" not in message  # class name only, never the detail


def test_noted_validation_appears_in_the_descriptor(home):
    cred.note_validation("pdb_redo", True, "Token accepted by PDB-REDO.")
    descriptor = cred.describe_credential("pdb_redo")
    assert descriptor["lastValidationOk"] is True
    assert descriptor["lastValidated"] is not None


def test_pdb_redo_validator_maps_status_codes(home, monkeypatch):
    """The validator must distinguish a bad token from an unreachable service.

    A user who is told "token rejected" when the network is down will go hunting
    for a typo that is not there.
    """
    import requests

    class Response:
        def __init__(self, status_code, payload=None):
            self.status_code = status_code
            self._payload = payload

        def json(self):
            if self._payload is None:
                raise ValueError("no json")
            return self._payload

    monkeypatch.setattr(
        requests, "get", lambda *a, **k: Response(200, [{"id": 1}, {"id": 2}])
    )
    ok, message = cred._validate_pdb_redo(VALUES)
    assert ok is True and "2 run(s)" in message

    monkeypatch.setattr(requests, "get", lambda *a, **k: Response(401))
    ok, message = cred._validate_pdb_redo(VALUES)
    assert ok is False and "rejected" in message

    def unreachable(*a, **k):
        raise requests.ConnectionError("dns")

    monkeypatch.setattr(requests, "get", unreachable)
    ok, message = cred._validate_pdb_redo(VALUES)
    assert ok is False and "Could not reach" in message

    ok, message = cred._validate_pdb_redo({"token_id": "x"})
    assert ok is False and "required" in message


def test_pdb_redo_validator_never_echoes_the_secret(home, monkeypatch):
    import requests

    class Response:
        status_code = 500

        def json(self):
            raise ValueError

    monkeypatch.setattr(requests, "get", lambda *a, **k: Response())
    _, message = cred._validate_pdb_redo(VALUES)
    assert VALUES["token_secret"] not in message
    assert VALUES["token_id"] not in message
