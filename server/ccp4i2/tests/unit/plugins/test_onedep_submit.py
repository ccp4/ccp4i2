"""The Submit-to-OneDep task, driven against a fake OneDep (no network).

Covered: creating a deposition and waiting for processing; replacing the
files on an existing deposition; the three ways waiting ends (processed,
rejected, lost contact) and the backoff that gets there; key checking with
rotation; and the task's own validity rules. ``onedep_lib`` is required (the
deposition extra); without it the module is skipped.
"""

import os
from datetime import datetime

import pytest

onedep_lib = pytest.importorskip("onedep_lib")

from onedep_lib.apis.deposit.enums import Status  # noqa: E402
from onedep_lib.apis.deposit.models import DepositedFile, DepositStatus, WwPDBDeposition  # noqa: E402
from onedep_lib.enums import FileType  # noqa: E402

from ccp4i2.config import credentials as cred  # noqa: E402
from ccp4i2.core.tasks import get_plugin_class  # noqa: E402
from ccp4i2.wrappers.onedep_submit.script import onedep_client  # noqa: E402

DEP_ID = "D_1000000001"
SITE = "https://deposit-pdbe.wwpdb.org/deposition/api/v1/depositions/%s/view" % DEP_ID


def _status(state, details):
    return DepositStatus(status=state, action="submit", step="upload", details=details,
                         date=datetime(2026, 11, 1, 12, 0, 0))


class FakeClient:
    """Implements onedep_lib's ApiClient protocol from a scripted status list."""

    def __init__(self, statuses, fail_upload=None):
        self.statuses = list(statuses)
        self.calls = []
        self.files = []
        self.next_file_id = 100
        self.fail_upload = fail_upload

    @staticmethod
    def _deposition(email):
        return WwPDBDeposition(dep_id=DEP_ID, email=email, pdb_id=None, emdb_id=None,
                               bmrb_id=None, title="", hold_exp_date=None,
                               created=datetime(2026, 11, 1), last_login=datetime(2026, 11, 1),
                               site="PDBE", status=Status.DEP, site_url=SITE)

    def create_deposition(self, email, users, country, experiments, password=""):
        self.calls.append(("create", email, list(users), country.value,
                           [e.exp_type.value for e in experiments]))
        return self._deposition(email)

    def get_deposition(self, dep_id):
        self.calls.append(("get", dep_id))
        return self._deposition("x@y.z")

    def upload_file(self, dep_id, file_path, file_type, overwrite=False, **kw):
        if self.fail_upload:
            raise self.fail_upload
        if overwrite:
            self.files = [f for f in self.files if f.file_type != file_type]
        self.calls.append(("upload", dep_id, os.path.basename(file_path), file_type.value, overwrite))
        self.next_file_id += 1
        f = DepositedFile(file_id=self.next_file_id, name=os.path.basename(file_path),
                          file_type=file_type, created=datetime(2026, 11, 1))
        self.files.append(f)
        return f

    def get_files(self, dep_id):
        return list(self.files)

    def remove_file(self, dep_id, file_id):
        self.files = [f for f in self.files if f.file_id != file_id]
        return True

    def process(self, dep_id):
        self.calls.append(("process", dep_id))
        return _status("started", "Processing started")

    def get_status(self, dep_id):
        self.calls.append(("status", dep_id))
        if len(self.statuses) > 1:
            return self.statuses.pop(0)
        return self.statuses[0]

    def get_all_depositions(self):
        return []

    def update_metadata(self, *a, **k):
        raise AssertionError("not used")


class FakeCheckRunner:
    def __getattr__(self, name):
        def check(*args, **kwargs):
            class Report:
                issues = []
            return Report()
        return check


@pytest.fixture
def files(tmp_path):
    model = tmp_path / "model.cif"
    sf = tmp_path / "sf.cif"
    model.write_text("data_x\n_entry.id x\n")
    sf.write_text("data_x\n_entry.id x\n")
    return str(model), str(sf)


@pytest.fixture
def config(tmp_path):
    return onedep_client.make_config("the-key", tmp_path / "state")


def _submit(client, files, config, tmp_path, **kw):
    log = []
    result = onedep_client.submit(files[0], files[1], "who@example.org",
                                  ["0000-0002-3595-9807"], "United Kingdom", config,
                                  log=log.append, sleep=lambda s: None,
                                  _api_client=client, _check_runner=FakeCheckRunner(),
                                  _base_dir=tmp_path / "sessions", **kw)
    return result, log


# ---------------------------------------------------------------------------
# submit()
# ---------------------------------------------------------------------------


def test_creates_a_deposition_uploads_both_files_and_waits(files, config, tmp_path):
    client = FakeClient([_status("running", "Processing model.cif"),
                         _status("finished", "Processing done.")])
    result, log = _submit(client, files, config, tmp_path)
    assert result["dep_id"] == DEP_ID
    assert result["site_url"] == SITE
    assert result["status"] == "finished"
    kinds = [c[0] for c in client.calls]
    assert kinds[:3] == ["create", "upload", "upload"]
    assert "process" in kinds
    assert client.calls[0][2] == ["0000-0002-3595-9807"]
    assert client.calls[0][3] == "United Kingdom"
    assert client.calls[0][4] == ["xray"]
    assert sorted(f["type"] for f in result["files"]) == sorted([
        FileType.MMCIF_COORD.value, FileType.CRYSTAL_STRUC_FACTORS.value])
    assert any("Processing done" in line for line in log)
    # the local session lives where we said, not in the job directory
    assert (tmp_path / "sessions").exists()


def test_a_rejected_upload_fails_with_the_servers_words_and_the_split_category_hint(
        files, config, tmp_path):
    client = FakeClient([_status("error", "Bad model file")])
    with pytest.raises(onedep_client.OneDepProcessingError) as err:
        _submit(client, files, config, tmp_path)
    assert err.value.dep_id == DEP_ID
    assert "Bad model file" in str(err.value)
    assert "two places" in str(err.value)


def test_other_rejections_carry_no_hint(files, config, tmp_path):
    client = FakeClient([_status("error", "Reflection file unreadable")])
    with pytest.raises(onedep_client.OneDepProcessingError) as err:
        _submit(client, files, config, tmp_path)
    assert "Reflection file unreadable" in str(err.value)
    assert "two places" not in str(err.value)


def test_replacing_files_on_an_existing_deposition_overwrites_and_reprocesses(
        files, config, tmp_path):
    client = FakeClient([_status("finished", "Processing done.")])
    client.files = [DepositedFile(file_id=7, name="old_model.cif", file_type=FileType.MMCIF_COORD,
                                  created=datetime(2026, 10, 31))]
    result, _ = _submit(client, files, config, tmp_path, replace_dep_id=DEP_ID)
    assert result["dep_id"] == DEP_ID
    assert ("create",) not in [c[:1] for c in client.calls]
    uploads = [c for c in client.calls if c[0] == "upload"]
    assert all(c[4] is True for c in uploads)
    assert ("process", DEP_ID) in client.calls
    # the old model is gone, the two new files are there
    assert sorted(f.name for f in client.files) == ["model.cif", "sf.cif"]


def test_auth_failure_is_an_auth_error_not_a_crash(files, config, tmp_path):
    client = FakeClient([_status("finished", "x")],
                        fail_upload=onedep_lib.exceptions.AuthError("expired"))
    with pytest.raises(onedep_client.OneDepAuthError):
        _submit(client, files, config, tmp_path, replace_dep_id=DEP_ID)


# ---------------------------------------------------------------------------
# wait_for_processing()
# ---------------------------------------------------------------------------


def test_polling_backs_off_caps_and_resets_on_change():
    states = [_status("running", "a")] * 12 + [_status("running", "b")] * 3 + [_status("finished", "done")]
    slept = []
    status = onedep_client.wait_for_processing(
        lambda: states.pop(0), sleep=slept.append, log=lambda s: None, clock=lambda: 0.0)
    assert status.status == "finished"
    assert max(slept) == onedep_client.POLL_MAX_SECONDS
    assert slept[0] == onedep_client.POLL_INITIAL_SECONDS
    # the change from "a" to "b" reset the interval
    assert slept[12] == onedep_client.POLL_INITIAL_SECONDS


def test_a_transient_failure_does_not_end_the_wait():
    items = [RuntimeError("blip"), _status("running", "a"), RuntimeError("blip"),
             _status("finished", "done")]

    def get():
        item = items.pop(0)
        if isinstance(item, Exception):
            raise item
        return item

    status = onedep_client.wait_for_processing(get, sleep=lambda s: None,
                                               log=lambda s: None, clock=lambda: 0.0)
    assert status.status == "finished"


def test_consecutive_failures_give_up_as_unreachable():
    def get():
        raise RuntimeError("down")

    with pytest.raises(onedep_client.OneDepUnreachable):
        onedep_client.wait_for_processing(get, sleep=lambda s: None, log=lambda s: None,
                                          clock=lambda: 0.0)


# ---------------------------------------------------------------------------
# check_key()
# ---------------------------------------------------------------------------


def _rotate_to(config, refresh, site=None):
    """What the library does on a successful call: rotate, maybe move site, persist."""
    from onedep_lib.auths.token import TokenStore
    store = TokenStore(config)
    if site:
        config.hostname = site
    store._store_tokens_for_key(onedep_lib.config._hostname_to_fqdn_key(config.hostname),
                               "access", refresh)


def test_check_key_reports_the_rotated_key_and_its_site(tmp_path):
    def fake_check(config):
        _rotate_to(config, "rotated-refresh", "https://deposit-pdbe.wwpdb.org/deposition")
        return True

    ok, message, new = onedep_client.check_key("the-key", base=tmp_path, _check=fake_check)
    assert ok
    assert new == {"api_key": "rotated-refresh",
                   "site_host": "https://deposit-pdbe.wwpdb.org/deposition"}
    assert "the-key" not in message and "rotated-refresh" not in message
    # the token landed where we said, owner-readable only
    toml = tmp_path / "config.toml"
    assert toml.exists()
    import stat
    assert stat.S_IMODE(toml.stat().st_mode) == 0o600


def test_check_key_probes_in_a_throwaway_directory_by_default(tmp_path, monkeypatch):
    probe_dir = tmp_path / "probe"
    probe_dir.mkdir()
    monkeypatch.setattr(onedep_client, "ephemeral_dir", lambda: probe_dir)
    seen = {}

    def fake_check(config):
        seen["path"] = str(config.config_path)
        _rotate_to(config, "rotated-refresh")
        return True

    ok, _, new = onedep_client.check_key("the-key", _check=fake_check)
    assert ok and new["api_key"] == "rotated-refresh"
    assert seen["path"].startswith(str(probe_dir))
    assert not probe_dir.exists()           # nothing persists after the probe


def test_check_key_rejection_mentions_single_use(tmp_path):
    ok, message, new = onedep_client.check_key("the-key", base=tmp_path,
                                               _check=lambda config: False)
    assert not ok and new is None
    assert "generate a new one" in message
    assert "the-key" not in message


def test_check_key_without_a_key():
    ok, message, new = onedep_client.check_key("")
    assert not ok and new is None


def test_rotated_credential_ignores_environment_tokens(tmp_path, monkeypatch):
    """Only the file's token is a rotation; a deployment's env token is not."""
    monkeypatch.setenv("ONEDEP_REFRESH_TOKEN", "env-token")
    config = onedep_client.make_config("the-key", tmp_path)
    assert onedep_client.rotated_credential(config, "the-key") is None
    _rotate_to(config, "file-token")
    assert onedep_client.rotated_credential(config, "the-key")["api_key"] == "file-token"


def test_make_config_targets_the_stored_site(tmp_path):
    config = onedep_client.make_config("k", tmp_path,
                                       hostname="https://deposit-pdbe.wwpdb.org/deposition")
    assert config.hostname == "https://deposit-pdbe.wwpdb.org/deposition"
    assert onedep_client.make_config("k", tmp_path).hostname.startswith("https://deposit.wwpdb.org")


# ---------------------------------------------------------------------------
# Error shapes from the API
# ---------------------------------------------------------------------------


def test_a_deposit_error_answer_ends_the_wait():
    from onedep_lib.apis.deposit.models import DepositError
    answers = [_status("running", "a"), DepositError(code="E42", message="session expired")]
    with pytest.raises(onedep_client.OneDepProcessingError) as err:
        onedep_client.wait_for_processing(lambda: answers.pop(0), sleep=lambda s: None,
                                          log=lambda s: None, clock=lambda: 0.0)
    assert "E42" in str(err.value) and "session expired" in str(err.value)


def test_replace_fails_when_process_itself_answers_with_an_error(files, config, tmp_path):
    from onedep_lib.apis.deposit.models import DepositError

    class RefusingProcess(FakeClient):
        def process(self, dep_id):
            self.calls.append(("process", dep_id))
            return DepositError(code="E1", message="cannot process")

    # a stale "finished" from the previous processing must not be read as success
    client = RefusingProcess([_status("finished", "Processing done.")])
    with pytest.raises(onedep_client.OneDepProcessingError) as err:
        _submit(client, files, config, tmp_path, replace_dep_id=DEP_ID)
    assert "cannot process" in str(err.value)
    assert ("status", DEP_ID) not in client.calls


def test_replace_of_a_missing_or_foreign_deposition_is_explained(files, config, tmp_path):
    class NotFound(FakeClient):
        def get_deposition(self, dep_id):
            raise onedep_lib.exceptions.ApiError("not found", 404)

    client = NotFound([_status("finished", "x")])
    with pytest.raises(onedep_client.OneDepRequestError) as err:
        _submit(client, files, config, tmp_path, replace_dep_id=DEP_ID)
    assert "404" in str(err.value) and DEP_ID in str(err.value)
    assert "ORCID" in str(err.value)


# ---------------------------------------------------------------------------
# The task's validity rules
# ---------------------------------------------------------------------------


@pytest.fixture
def home(tmp_path, monkeypatch):
    monkeypatch.setenv("CCP4I2_HOME", str(tmp_path))
    monkeypatch.setattr(cred, "_keyring", lambda: None)
    for field in cred.ONEDEP.fields:
        monkeypatch.delenv(field.env, raising=False)
        monkeypatch.delenv(cred._session_env_key(cred._storage_key("onedep"), field.name),
                           raising=False)
    return tmp_path


def _plugin(tmp_path, files):
    (tmp_path / "job").mkdir(exist_ok=True)
    plugin = get_plugin_class("onedep_submit")(workDirectory=str(tmp_path / "job"), name="a")
    inp = plugin.container.inputData
    inp.MODEL_CIF.setFullPath(files[0])
    inp.SF_CIF.setFullPath(files[1])
    return plugin


def _codes(plugin):
    return sorted(r["code"] for r in plugin.validity()._reports)


def test_validity_wants_a_key_and_sane_depositor_details(home, files, tmp_path):
    plugin = _plugin(tmp_path, files)
    cpar = plugin.container.controlParameters
    cpar.EMAIL.set("not-an-email")
    assert 200 in _codes(plugin)     # no key
    assert 202 in _codes(plugin)     # bad e-mail, no ORCID


def test_validity_is_clean_with_key_email_and_orcid(home, files, tmp_path):
    cred.set_credential("onedep", {"api_key": "k"})
    plugin = _plugin(tmp_path, files)
    cpar = plugin.container.controlParameters
    cpar.EMAIL.set("who@example.org")
    cpar.ORCIDS.set(["0000-0002-3595-9807"])
    assert 200 not in _codes(plugin)
    assert 202 not in _codes(plugin)


def test_validity_rejects_a_malformed_orcid_and_deposition_id(home, files, tmp_path):
    cred.set_credential("onedep", {"api_key": "k"})
    plugin = _plugin(tmp_path, files)
    cpar = plugin.container.controlParameters
    cpar.EMAIL.set("who@example.org")
    cpar.ORCIDS.set(["0000-0002-3595-980"])
    cpar.REPLACE_DEPOSITION_ID.set("1292162556")
    reports = [r for r in plugin.validity()._reports if r["code"] == 202]
    assert reports
    text = reports[0]["details"]
    assert "ORCID" in text and "D_" in text


def test_run_time_validity_never_probes_the_key(home, files, tmp_path, monkeypatch):
    """A probe rotates the key; the only authenticated call per job is the submission."""
    cred.set_credential("onedep", {"api_key": "k"})

    def forbidden(*a, **k):
        raise AssertionError("runTimeValidity probed the credential")

    monkeypatch.setattr(cred, "validate_credential", forbidden)
    monkeypatch.setattr(onedep_client, "check_key", forbidden)
    plugin = _plugin(tmp_path, files)
    cpar = plugin.container.controlParameters
    cpar.EMAIL.set("who@example.org")
    cpar.ORCIDS.set(["0000-0002-3595-9807"])
    codes = [r["code"] for r in plugin.runTimeValidity()._reports]
    assert 201 not in codes


def test_run_time_validity_rejects_an_unknown_country(home, files, tmp_path):
    cred.set_credential("onedep", {"api_key": "k"})
    plugin = _plugin(tmp_path, files)
    cpar = plugin.container.controlParameters
    cpar.EMAIL.set("who@example.org")
    cpar.ORCIDS.set(["0000-0002-3595-9807"])
    cpar.COUNTRY.set("Narnia")
    assert 202 in [r["code"] for r in plugin.runTimeValidity()._reports]


def _run_with_fake_submit(plugin, monkeypatch, outcome):
    """Drive startProcess with a submit() that rotates the key, then does `outcome`."""
    def fake_submit(model, sf, email, orcids, country, config, **kw):
        _rotate_to(config, "rotated-by-run", "https://deposit-pdbe.wwpdb.org/deposition")
        if isinstance(outcome, Exception):
            raise outcome
        return outcome

    monkeypatch.setattr(onedep_client, "submit", fake_submit)
    cpar = plugin.container.controlParameters
    cpar.EMAIL.set("who@example.org")
    cpar.ORCIDS.set(["0000-0002-3595-9807"])
    return plugin.startProcess()


def test_a_failed_run_still_keeps_the_rotated_key(home, files, tmp_path, monkeypatch, capsys):
    """The ephemeral state dir must not be removed before the rotated key is read back."""
    cred.set_credential("onedep", {"api_key": "first"})
    monkeypatch.delenv("CCP4I2_LOCAL_SESSION_TOKEN", raising=False)   # hosted: ephemeral dir
    plugin = _plugin(tmp_path, files)
    rc = _run_with_fake_submit(plugin, monkeypatch,
                               onedep_client.OneDepProcessingError("rejected", dep_id=DEP_ID))
    assert rc == plugin.FAILED
    assert cred.get_credential("onedep") == {"api_key": "rotated-by-run",
                                             "site_host": "https://deposit-pdbe.wwpdb.org/deposition"}
    assert plugin._ephemeral is None       # and the state dir is gone afterwards
    assert "new key has been stored" in capsys.readouterr().out


def test_a_successful_run_keeps_the_rotated_key_too(home, files, tmp_path, monkeypatch):
    cred.set_credential("onedep", {"api_key": "first"})
    plugin = _plugin(tmp_path, files)
    rc = _run_with_fake_submit(plugin, monkeypatch, {
        "dep_id": DEP_ID, "site_url": SITE, "session_id": "s", "files": [],
        "status": "finished", "details": "Processing done.", "replaced": False})
    assert rc == plugin.SUCCEEDED
    assert cred.get_credential("onedep")["api_key"] == "rotated-by-run"


def test_an_environment_key_cannot_be_rotated_and_the_log_says_so(home, files, tmp_path,
                                                                     monkeypatch, capsys):
    monkeypatch.setenv("ONEDEP_API_KEY", "deployment-key")
    plugin = _plugin(tmp_path, files)
    rc = _run_with_fake_submit(plugin, monkeypatch, {
        "dep_id": DEP_ID, "site_url": SITE, "session_id": "s", "files": [],
        "status": "finished", "details": "done", "replaced": False})
    assert rc == plugin.SUCCEEDED
    assert cred.get_credential("onedep")["api_key"] == "deployment-key"
    assert "WARNING" in capsys.readouterr().out


def test_the_report_escapes_what_the_server_said():
    from ccp4i2.wrappers.onedep_submit.script.onedep_submit_report import safe
    assert safe("<b>Bad</b> & worse") == "&lt;b&gt;Bad&lt;/b&gt; &amp; worse"
    assert safe(None) == ""


def test_the_key_never_appears_in_the_params_xml(home, files, tmp_path):
    """Rule zero: the credential is not a task parameter."""
    cred.set_credential("onedep", {"api_key": "secret-key-value"})
    plugin = _plugin(tmp_path, files)
    xml = plugin.container.__str__() if hasattr(plugin.container, "__str__") else ""
    from lxml import etree
    root = etree.Element("probe")
    try:
        plugin.container.saveDataToXml(fileName=str(tmp_path / "params.xml"))
        text = (tmp_path / "params.xml").read_text()
    except Exception:
        text = xml
    assert "secret-key-value" not in text
    assert "api_key" not in text.lower()
