"""Thin, testable layer over ``onedep_lib`` for the Submit-to-OneDep task.

Everything that talks to wwPDB lives here, behind functions that take the
API client as an argument, so the wrapper and the credential validator can be
exercised with a fake client and no network. ``onedep_lib`` itself is imported
lazily: the task plugin must import on a server with neither CCP4 nor the
optional ``ccp4i2[deposition]`` extra installed.

Two facts about OneDep shape this module:

* The DepUI "Generate Key" value is a long-lived *refresh* token bound to the
  depositor's ORCID. Every use rotates it: the old key dies and a new pair
  comes back. So a key can never be validated "read-only", and whoever holds
  the key must be told the replacement. :func:`check_key` returns it.
* The server's quick upload check answers ``Bad model file`` within seconds,
  with no per-file detail, when an mmCIF category is split into two places in
  the block (a trap gemmi's ``set_pair`` lays by appending at the end).
  :func:`explain_failure` says so, because the message alone sends people
  hunting through the atoms.
"""

import logging
import os
import stat
import tempfile
import time
from pathlib import Path

logger = logging.getLogger(f"ccp4i2.{__name__}")

INSTALL_HINT = (
    "The onedep_lib package is not installed. Install the deposition extra "
    "(pip install 'ccp4i2[deposition]' or pip install onedep_lib) in the "
    "Python that runs CCP4i2 jobs."
)

# Polling, per WEB_API_TASK_GUIDE: responsive at first, then back off.
POLL_INITIAL_SECONDS = 5.0
POLL_MAX_SECONDS = 60.0
POLL_BACKOFF = 1.5
POLL_MAX_CONSECUTIVE_FAILURES = 5
POLL_HEARTBEAT_SECONDS = 600.0
POLL_TIMEOUT_SECONDS = 6 * 3600.0

TERMINAL_STATES = ("finished", "error", "failed")

ORCID_RE = r"^\d{4}-\d{4}-\d{4}-\d{3}[\dX]$"


class OneDepError(Exception):
    """Base class; ``dep_id`` carries the remote deposition when there is one."""

    def __init__(self, message, dep_id=None, site_url=None):
        super().__init__(message)
        self.dep_id = dep_id
        self.site_url = site_url


class OneDepAuthError(OneDepError):
    """The key was rejected. Retrying will not help; the user must act."""


class OneDepUnreachable(OneDepError):
    """We lost contact. Nothing is lost remotely; the deposition persists."""


class OneDepProcessingError(OneDepError):
    """OneDep processed the upload and rejected it."""


class OneDepMissing(OneDepError):
    """``onedep_lib`` is not importable."""


def lib():
    """Import and return the ``onedep_lib`` package, or raise OneDepMissing."""
    try:
        import onedep_lib  # noqa: F401
        from onedep_lib import config, dsp, enums  # noqa: F401
    except ImportError as err:
        raise OneDepMissing(INSTALL_HINT) from err
    return onedep_lib


def library_available():
    try:
        lib()
    except OneDepMissing:
        return False
    return True


def valid_orcid(value):
    import re

    return bool(re.match(ORCID_RE, (value or "").strip()))


def valid_email(value):
    import re

    return bool(re.match(r"^[^@\s]+@[^@\s]+\.[^@\s]+$", (value or "").strip()))


def country_names():
    """Country labels OneDep accepts, in the library's order."""
    enums = lib().enums
    return [c.value for c in enums.Country]


def country_from_name(name):
    enums = lib().enums
    for c in enums.Country:
        if c.value == name:
            return c
    raise ValueError(f"OneDep does not list a country called {name!r}")


# ---------------------------------------------------------------------------
# Where the library's own state goes
# ---------------------------------------------------------------------------


def private_dir(base=None):
    """A 0700 directory for onedep_lib's config and session files.

    Never the job directory: that is exported in project zips. On a desktop
    this is under the per-user CCP4i2 home; a hosted server passes a
    temporary ``base`` that it removes when the job ends.
    """
    if base is None:
        from ccp4i2.config.preferences import ccp4i2_home

        base = ccp4i2_home()
    path = Path(base) / "onedep"
    path.mkdir(parents=True, exist_ok=True)
    try:
        os.chmod(str(path), stat.S_IRWXU)
    except OSError:
        logger.debug("could not chmod %s", path, exc_info=True)
    return path


def ephemeral_dir():
    """A private directory that exists only for one job (hosted deployments)."""
    return Path(tempfile.mkdtemp(prefix="ccp4i2-onedep-"))


def make_config(api_key, base):
    """A DepositConfig that keeps every token and session under ``base``.

    ``api_key`` is used as the refresh token; the library writes the rotated
    pair into ``base/config.toml`` (0600, its own doing) and session state into
    ``base/sessions``.
    """
    config_mod = lib().config
    base = Path(base)
    return config_mod.DepositConfig.load(
        refresh_token=api_key,
        config_path=base / "config.toml",
        session_dir=base / "sessions",
        schema_cache_dir=base / "schemas",
    )


def rotated_key(config, original):
    """The refresh token now on disk for ``config``'s site, if it changed.

    After any authenticated call the library has replaced the key it was given
    (on disk and on ``config`` itself, so compare with ``original``, the key
    the caller started from); this reads the replacement back so the caller
    can store it. None when nothing is on disk or it is still ``original``.
    """
    config_mod = lib().config
    try:
        reloaded = config_mod.DepositConfig.load(config_path=config.config_path,
                                                 hostname=config.hostname)
    except Exception:
        logger.debug("could not reload onedep config", exc_info=True)
        return None
    token = reloaded.refresh_token
    if not token or token == original:
        return None
    return token


# ---------------------------------------------------------------------------
# Credential check
# ---------------------------------------------------------------------------


def check_key(api_key, base=None, _check=None):
    """Probe a DepUI key against OneDep.

    Returns ``(ok, message, new_key)``. ``new_key`` is the rotated refresh
    token when the probe succeeded and the library rotated it, else None.
    The message never contains the key.
    """
    api_key = (api_key or "").strip()
    if not api_key:
        return False, "A OneDep API key is required.", None
    try:
        onedep = lib()
    except OneDepMissing as err:
        return False, str(err), None

    base = Path(base) if base else private_dir()
    config = make_config(api_key, base)
    check = _check or onedep.dsp.check_auth_key
    try:
        ok = bool(check(config))
    except Exception as err:  # the library wraps most things; be safe
        return False, f"Could not check the key ({err.__class__.__name__}).", None
    if not ok:
        return False, (
            "OneDep rejected this key. Keys are single-use until rotated: if "
            "this one was used elsewhere, generate a new one in the deposition "
            "interface."
        ), None
    return True, "Key accepted by OneDep.", rotated_key(config, api_key)


# ---------------------------------------------------------------------------
# Submission
# ---------------------------------------------------------------------------


def wait_for_processing(get_status, sleep=time.sleep, log=print, clock=time.monotonic,
                        initial=POLL_INITIAL_SECONDS, maximum=POLL_MAX_SECONDS,
                        backoff=POLL_BACKOFF, timeout=POLL_TIMEOUT_SECONDS,
                        max_failures=POLL_MAX_CONSECUTIVE_FAILURES):
    """Poll ``get_status()`` until OneDep reports a terminal state.

    Returns the final status object. Prints on change (and a heartbeat every
    ten minutes), backs off between polls, resets the interval when the status
    changes, survives transient failures, and gives up only after
    ``max_failures`` consecutive ones or ``timeout`` seconds.
    """
    delay = initial
    failures = 0
    last_key = None
    started = clock()
    last_heartbeat = started
    status = None
    while True:
        try:
            status = get_status()
            failures = 0
        except Exception as err:
            failures += 1
            log("WARNING: status poll failed (%s), attempt %d of %d"
                % (err.__class__.__name__, failures, max_failures))
            if failures >= max_failures:
                raise OneDepUnreachable(
                    "Lost contact with OneDep after %d consecutive failed polls."
                    % failures) from err
            sleep(delay)
            delay = min(delay * backoff, maximum)
            continue

        state = getattr(status, "status", None)
        key = (state, getattr(status, "step", None), getattr(status, "details", None))
        now = clock()
        if key != last_key:
            log("OneDep: %s - %s" % (state, getattr(status, "details", "")))
            last_key = key
            delay = initial
            last_heartbeat = now
        elif now - last_heartbeat >= POLL_HEARTBEAT_SECONDS:
            log("OneDep: still %s (%d min)" % (state, int((now - started) // 60)))
            last_heartbeat = now
        if state in TERMINAL_STATES:
            return status
        if now - started > timeout:
            raise OneDepUnreachable(
                "Gave up waiting for OneDep processing after %d minutes; the "
                "deposition is still there." % int(timeout // 60))
        sleep(delay)
        delay = min(delay * backoff, maximum)


def explain_failure(status):
    """A user-facing sentence for a terminal ``error`` status."""
    details = getattr(status, "details", "") or ""
    text = "OneDep rejected the upload: %s." % (details.rstrip(".") or "no details given")
    if "bad model file" in details.lower():
        text += (
            " OneDep answers 'Bad model file' when the model mmCIF is not "
            "parseable as one block of whole categories; the usual cause is a "
            "category whose items appear in two places (for example _refine "
            "items added after _atom_site). Rewrite the file so every category "
            "appears exactly once, and check it with 'gemmi validate'."
        )
    return text


def _file_types(enums):
    return {
        "model": enums.FileType.MMCIF_COORD,
        "sf": enums.FileType.CRYSTAL_STRUC_FACTORS,
    }


def submit(model_path, sf_path, email, orcids, country_name, config,
           experiment_type="xray", replace_dep_id=None, log=print,
           sleep=time.sleep, _api_client=None, _check_runner=None, _base_dir=None):
    """Create (or refresh) a OneDep deposition and wait for processing.

    Returns a dict: ``dep_id``, ``site_url``, ``session_id``, ``files``
    (list of {name, type, remote_id}), ``status``, ``details``.
    Raises OneDepAuthError / OneDepUnreachable / OneDepProcessingError.
    """
    onedep = lib()
    enums = onedep.enums
    exceptions = onedep.exceptions
    types = _file_types(enums)
    country = country_from_name(country_name)
    exp_type = enums.ExperimentType(experiment_type)

    def auth_wrapped(fn, *args, **kwargs):
        try:
            return fn(*args, **kwargs)
        except exceptions.AuthError as err:
            raise OneDepAuthError(
                "OneDep rejected the API key (%s). Set a fresh key and re-run."
                % err, dep_id=replace_dep_id) from err
        except exceptions.ApiUnreachableError as err:
            raise OneDepUnreachable(
                "Could not reach OneDep (%s)." % err, dep_id=replace_dep_id) from err

    result = {"dep_id": None, "site_url": None, "session_id": None, "files": [],
              "status": None, "details": ""}

    if replace_dep_id:
        # Replace the files on an existing deposition: upload with overwrite
        # (the library removes the server's file of the same type first), then
        # ask OneDep to process. No local session is needed for this path.
        if _api_client is not None:
            client = _api_client
        else:
            from onedep_lib.apis.deposit.client import HttpApiClient
            from onedep_lib.auths.token import TokenStore

            client = HttpApiClient(config, auth_provider=TokenStore(config))
        dep_id = replace_dep_id
        result["dep_id"] = dep_id
        dep = auth_wrapped(client.get_deposition, dep_id)
        result["site_url"] = getattr(dep, "site_url", None)
        for label, path in (("model", model_path), ("sf", sf_path)):
            log("Uploading %s as %s (replacing)" % (os.path.basename(path), types[label].value))
            uploaded = auth_wrapped(client.upload_file, dep_id, path, types[label], overwrite=True)
            result["files"].append({"name": os.path.basename(path), "type": types[label].value,
                                    "remote_id": getattr(uploaded, "file_id", None)})
        auth_wrapped(client.process, dep_id)
        get_status = lambda: auth_wrapped(client.get_status, dep_id)  # noqa: E731
    else:
        dep = onedep.dsp.deposit_init(
            email=email, users=list(orcids), country=country, config=config,
            experiment_type=exp_type, coordinates=True,
            _base_dir=_base_dir, _api_client=_api_client, _check_runner=_check_runner)
        result["session_id"] = dep.session_id
        dep.add_file(model_path, types["model"])
        dep.add_file(sf_path, types["sf"])
        report = dep.check_required_files()
        for issue in getattr(report, "issues", []) or []:
            log("check: %s %s" % (getattr(issue, "severity", ""), getattr(issue, "message", "")))
        log("Creating the deposition and uploading %s and %s"
            % (os.path.basename(model_path), os.path.basename(sf_path)))
        dep_id = auth_wrapped(dep.deposit)
        result["dep_id"] = dep_id
        result["site_url"] = dep.site_url
        client = dep._api_client
        try:
            for f in client.get_files(dep_id):
                result["files"].append({"name": f.name, "type": f.file_type.value,
                                        "remote_id": f.file_id})
        except Exception:
            logger.debug("could not list remote files", exc_info=True)
            for label, path in (("model", model_path), ("sf", sf_path)):
                result["files"].append({"name": os.path.basename(path),
                                        "type": types[label].value, "remote_id": None})
        get_status = lambda: auth_wrapped(dep.get_status)  # noqa: E731

    log("OneDep deposition %s" % dep_id)
    try:
        status = wait_for_processing(get_status, sleep=sleep, log=log)
    except OneDepUnreachable as err:
        err.dep_id = dep_id
        err.site_url = result["site_url"]
        raise
    result["status"] = getattr(status, "status", None)
    result["details"] = getattr(status, "details", "") or ""
    if result["status"] != "finished":
        raise OneDepProcessingError(explain_failure(status), dep_id=dep_id,
                                    site_url=result["site_url"])
    return result
