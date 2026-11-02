"""A validator may hand back a replacement credential (token rotation).

OneDep's API key is a refresh token that its first use consumes, so the probe
that says "this key is good" also makes it the wrong key. The store has to
keep the replacement, in the layer the original came from, or the Test button
locks the user out. No network: the OneDep validator is stubbed.
"""

import pytest

from ccp4i2.config import credentials as cred


@pytest.fixture
def home(tmp_path, monkeypatch):
    monkeypatch.setenv("CCP4I2_HOME", str(tmp_path))
    monkeypatch.setattr(cred, "_keyring", lambda: None)
    for spec in (cred.ONEDEP, cred.PDB_REDO):
        for field in spec.fields:
            monkeypatch.delenv(field.env, raising=False)
            monkeypatch.delenv(
                cred._session_env_key(cred._storage_key(spec.name), field.name),
                raising=False,
            )
    cred._validation_cache.clear()
    return tmp_path


def _rotating(monkeypatch, new_key="rotated-key", ok=True):
    """Stub the OneDep validator: accept, and report a rotated key."""
    calls = []

    def validator(values):
        calls.append(dict(values))
        return ok, "Key accepted by OneDep.", ({"api_key": new_key} if ok else None)

    spec = cred.ONEDEP
    monkeypatch.setitem(cred.CREDENTIALS, "onedep", cred.CredentialSpec(
        **{**spec.__dict__, "validator": validator}))
    return calls


def test_onedep_is_registered_with_one_secret_field_and_an_optional_site():
    spec = cred.get_spec("onedep")
    assert spec is not None
    assert [f.name for f in spec.fields] == ["api_key", "site_host"]
    assert spec.fields[0].secret is True and spec.fields[0].required is True
    assert spec.fields[0].env == "ONEDEP_API_KEY"
    assert spec.fields[1].secret is False and spec.fields[1].required is False


def test_the_credential_is_complete_without_the_optional_site(home):
    cred.set_credential("onedep", {"api_key": "k"})
    assert cred.get_credential("onedep") == {"api_key": "k"}
    descriptor = cred.describe_credential("onedep")
    assert descriptor["isSet"] is True
    assert [f["required"] for f in descriptor["fields"]] == [True, False]


def test_rotation_carries_the_site_with_the_key(home, monkeypatch):
    cred.set_credential("onedep", {"api_key": "first"})
    spec = cred.ONEDEP
    monkeypatch.setitem(cred.CREDENTIALS, "onedep", cred.CredentialSpec(
        **{**spec.__dict__, "validator": lambda values: (
            True, "ok", {"api_key": "second", "site_host": "https://deposit-pdbe.wwpdb.org/deposition"})}))
    ok, _ = cred.validate_credential("onedep")
    assert ok
    assert cred.get_credential("onedep") == {
        "api_key": "second", "site_host": "https://deposit-pdbe.wwpdb.org/deposition"}
    assert cred.describe_credential("onedep")["hint"] == "tion"


def test_the_deployment_key_is_not_probed(home, monkeypatch):
    """A probe would spend a key nothing here can replace."""
    monkeypatch.setenv("ONEDEP_API_KEY", "deployment-key")
    called = []
    from ccp4i2.wrappers.onedep_submit.script import onedep_client
    monkeypatch.setattr(onedep_client, "check_key",
                        lambda *a, **k: called.append(1) or (True, "probed", None))
    ok, message = cred.validate_credential("onedep")
    assert ok and not called
    assert "not probed" in message


def test_validation_persists_the_rotated_key_in_the_file_store(home, monkeypatch):
    cred.set_credential("onedep", {"api_key": "first-key"})
    calls = _rotating(monkeypatch, "second-key")
    ok, message = cred.validate_credential("onedep")
    assert ok and "accepted" in message
    assert calls == [{"api_key": "first-key"}]
    assert cred.get_credential("onedep") == {"api_key": "second-key"}
    assert cred.resolve_credential("onedep").source == "file"


def test_validation_keeps_a_session_credential_in_the_session(home, monkeypatch):
    cred.set_credential("onedep", {"api_key": "first-key"}, persistence=cred.PERSIST_SESSION)
    _rotating(monkeypatch, "second-key")
    ok, _ = cred.validate_credential("onedep")
    assert ok
    resolved = cred.resolve_credential("onedep")
    assert resolved.values == {"api_key": "second-key"}
    assert resolved.source == "session"
    # nothing reached the file store
    assert cred._load_file_store() == {}


def test_environment_credential_is_not_rewritten_and_the_user_is_told(home, monkeypatch):
    monkeypatch.setenv("ONEDEP_API_KEY", "deployment-key")
    _rotating(monkeypatch, "second-key")
    ok, message = cred.validate_credential("onedep")
    assert ok
    assert "could not be stored" in message
    assert cred.get_credential("onedep") == {"api_key": "deployment-key"}


def test_a_rejected_key_is_not_replaced(home, monkeypatch):
    cred.set_credential("onedep", {"api_key": "first-key"})
    _rotating(monkeypatch, ok=False)
    ok, _ = cred.validate_credential("onedep")
    assert not ok
    assert cred.get_credential("onedep") == {"api_key": "first-key"}


def test_two_value_validators_still_work(home, monkeypatch):
    """pdb_redo's validator returns (ok, message); the contract extension is additive."""
    cred.set_credential("pdb_redo", {"token_id": "id", "token_secret": "s"})
    spec = cred.PDB_REDO
    monkeypatch.setitem(cred.CREDENTIALS, "pdb_redo", cred.CredentialSpec(
        **{**spec.__dict__, "validator": lambda values: (True, "fine")}))
    assert cred.validate_credential("pdb_redo") == (True, "fine")


def test_descriptor_never_echoes_the_key(home):
    cred.set_credential("onedep", {"api_key": "do-not-show"})
    descriptor = cred.describe_credential("onedep")
    assert "do-not-show" not in repr(descriptor)
    assert descriptor["isSet"] is True
    assert descriptor["hint"] == ""
