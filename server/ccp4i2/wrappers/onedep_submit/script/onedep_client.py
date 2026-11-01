"""Thin, testable layer over ``onedep_lib`` for the Submit-to-OneDep task.

Everything that talks to wwPDB lives here, behind functions that take the
API client as an argument, so the wrapper and the credential validator can be
exercised with a fake client and no network. ``onedep_lib`` itself is imported
lazily: the task plugin must import on a server with neither CCP4 nor the
optional ``ccp4i2[deposition]`` extra installed.

Three facts about OneDep shape this module:

* The DepUI "Generate Key" value is a long-lived *refresh* token bound to the
  depositor's ORCID. Every authenticated call rotates it: the old key dies and
  a new pair comes back. So a key can never be probed "read-only", two
  processes sharing one key invalidate each other, and whoever holds the key
  must be told the replacement. :func:`check_key` and
  :func:`rotated_credential` return it; the task makes exactly one
  authenticated call per job, the submission itself.
* The central host (deposit.wwpdb.org) redirects to a regional site
  (deposit-pdbe, deposit-rcsb, ...) and the library exchanges the key for one
  scoped to that site, recording the site on the config. The site is kept
  with the key so the next run talks to the right host from the start.
  Whether the central host would accept a site-scoped refresh token is
  unverified, which is the reason not to rely on it.
* The server's quick upload check answers ``Bad model file`` within seconds,
  with no per-file detail, when an mmCIF category is split into two places in
  the block (a trap gemmi's ``set_pair`` lays by appending at the end).
  :func:`explain_failure` says so, because the message alone sends people
  hunting through the atoms.
"""

import logging
import os
import shutil
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
    """OneDep processed the upload and rejected it, or answered with an error."""


class OneDepRequestError(OneDepError):
    """OneDep refused a request outright (HTTP 4xx): wrong id, not ours, gone."""


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
    """A private directory that exists only for one job or one probe."""
    return Path(tempfile.mkdtemp(prefix="ccp4i2-onedep-"))


def make_config(api_key, base, hostname=None):
    """A DepositConfig that keeps every token and session under ``base``.

    ``api_key`` is used as the refresh token; the library writes the rotated
    pair into ``base/config.toml`` and session state into ``base/sessions``.
    ``hostname`` is the deposition site the key was last rotated against;
    without it the central host is used and the library follows its redirect.
    """
    config_mod = lib().config
    base = Path(base)
    overrides = dict(
        refresh_token=api_key,
        config_path=base / "config.toml",
        session_dir=base / "sessions",
        schema_cache_dir=base / "schemas",
    )
    if hostname:
        overrides["hostname"] = hostname
    return config_mod.DepositConfig.load(**overrides)


def secure_config(config):
    """Make the library's config file owner-readable only.

    onedep_lib writes config.toml with the default umask (0644); only the 0700
    directory protects it. Tighten the file itself whenever we have written it.
    """
    path = Path(config.config_path)
    if path.exists():
        try:
            os.chmod(str(path), stat.S_IRUSR | stat.S_IWUSR)
        except OSError:
            logger.debug("could not chmod %s", path, exc_info=True)


def _read_toml(path):
    try:
        import tomllib
    except ImportError:  # Python < 3.11: onedep_lib depends on tomli there
        import tomli as tomllib
    try:
        with open(path, "rb") as handle:
            return tomllib.load(handle)
    except (FileNotFoundError, NotADirectoryError):
        return {}


def rotated_credential(config, original):
    """The key OneDep rotated to, and the site it is scoped to, if it changed.

    After any authenticated call the library has replaced the key it was given
    (on disk, and on ``config`` itself, so compare with ``original``, the key
    the caller started from) and may have moved ``config.hostname`` to the
    regional site. The on-disk TOML is read directly rather than through
    ``DepositConfig.load``, which would merge ``ONEDEP_REFRESH_TOKEN`` from the
    environment and report a deployment's token as the rotation.

    Returns ``{"api_key": ..., "site_host": ...}`` or None.
    """
    config_mod = lib().config
    host = config.hostname
    key = config_mod._hostname_to_fqdn_key(host)
    if not key:
        return None
    entry = _read_toml(config.config_path).get("auths", {}).get(key)
    token = entry.get("refresh_token") if isinstance(entry, dict) else None
    if not token or not isinstance(token, str) or token == original:
        return None
    return {"api_key": token, "site_host": host}


# ---------------------------------------------------------------------------
# Credential check
# ---------------------------------------------------------------------------


def check_key(api_key, base=None, _check=None, hostname=None):
    """Probe a DepUI key against OneDep.

    Returns ``(ok, message, new_values)``. ``new_values`` is the rotated
    credential (key and site) when the probe succeeded and the library
    rotated it, else None. The probe runs in a throwaway directory removed
    afterwards: the rotated key goes back to the credential store, nothing
    needs to persist here. The message never contains the key.
    """
    api_key = (api_key or "").strip()
    if not api_key:
        return False, "A OneDep API key is required.", None
    try:
        onedep = lib()
    except OneDepMissing as err:
        return False, str(err), None

    own_base = base is None
    base = Path(base) if base else ephemeral_dir()
    try:
        config = make_config(api_key, base, hostname=hostname)
        check = _check or onedep.dsp.check_auth_key
        try:
            ok = bool(check(config))
        except Exception as err:  # the library wraps most things; be safe
            return False, f"Could not check the key ({err.__class__.__name__}).", None
        secure_config(config)
        if not ok:
            return False, (
                "OneDep rejected this key. Keys are single-use until rotated: if "
                "this one was used elsewhere, generate a new one in the deposition "
                "interface."
            ), None
        return True, "Key accepted by OneDep.", rotated_credential(config, api_key)
    finally:
        if own_base:
            shutil.rmtree(str(base), ignore_errors=True)


# ---------------------------------------------------------------------------
# Submission
# ---------------------------------------------------------------------------


def _is_error_object(status):
    """A DepositError (code/message) rather than a DepositStatus (status)."""
    return not hasattr(status, "status") and hasattr(status, "message")


def _error_text(status):
    code = getattr(status, "code", "")
    message = getattr(status, "message", "") or ""
    return "OneDep error %s: %s" % (code, message) if code else "OneDep error: %s" % message


def wait_for_processing(get_status, sleep=time.sleep, log=print, clock=time.monotonic,
                        initial=POLL_INITIAL_SECONDS, maximum=POLL_MAX_SECONDS,
                        backoff=POLL_BACKOFF, timeout=POLL_TIMEOUT_SECONDS,
                        max_failures=POLL_MAX_CONSECUTIVE_FAILURES):
    """Poll ``get_status()`` until OneDep reports a terminal state.

    Returns the final status object. Prints on change (and a heartbeat every
    ten minutes), backs off between polls, resets the interval when the status
    changes, survives transient failures, and gives up only after
    ``max_failures`` consecutive ones or ``timeout`` seconds. A DepositError
    answer (the API's error shape, with no status) is terminal and raises
    OneDepProcessingError.
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

        if _is_error_object(status):
            raise OneDepProcessingError(_error_text(status))
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
    (list of {name, type, remote_id}), ``status``, ``details``, ``replaced``.
    Raises OneDepAuthError / OneDepUnreachable / OneDepRequestError /
    OneDepProcessingError.
    """
    onedep = lib()
    enums = onedep.enums
    exceptions = onedep.exceptions
    types = _file_types(enums)
    country = country_from_name(country_name)
    exp_type = enums.ExperimentType(experiment_type)

    def guarded(fn, *args, **kwargs):
        try:
            return fn(*args, **kwargs)
        except exceptions.AuthError as err:
            raise OneDepAuthError(
                "OneDep rejected the API key (%s). Set a fresh key and re-run."
                % err, dep_id=replace_dep_id) from err
        except exceptions.ApiUnreachableError as err:
            raise OneDepUnreachable(
                "Could not reach OneDep (%s)." % err, dep_id=replace_dep_id) from err
        except exceptions.ApiError as err:
            code = getattr(err, "status_code", None)
            if code in (403, 404) and replace_dep_id:
                raise OneDepRequestError(
                    "OneDep answered HTTP %s for deposition %s: it does not exist, "
                    "or it is not open to the ORCID this key belongs to."
                    % (code, replace_dep_id), dep_id=replace_dep_id) from err
            if code in (401, 403):
                raise OneDepAuthError(
                    "OneDep refused the request (HTTP %s): %s" % (code, err),
                    dep_id=replace_dep_id) from err
            raise OneDepRequestError(
                "OneDep refused the request (HTTP %s): %s" % (code, err),
                dep_id=replace_dep_id) from err

    result = {"dep_id": None, "site_url": None, "session_id": None, "files": [],
              "status": None, "details": "", "replaced": bool(replace_dep_id)}

    if replace_dep_id:
        # Replace the files on an existing deposition. overwrite=True makes the
        # library delete the server's file of the same type before uploading,
        # so a failed upload leaves the deposition without that file until the
        # task is re-run. No local session is needed for this path.
        if _api_client is not None:
            client = _api_client
        else:
            from onedep_lib.apis.deposit.client import HttpApiClient
            from onedep_lib.auths.token import TokenStore

            client = HttpApiClient(config, auth_provider=TokenStore(config))
        dep_id = replace_dep_id
        result["dep_id"] = dep_id
        dep = guarded(client.get_deposition, dep_id)
        result["site_url"] = getattr(dep, "site_url", None)
        log("Replacing the files on deposition %s; each file is removed from the "
            "server before its replacement is uploaded" % dep_id)
        for label, path in (("model", model_path), ("sf", sf_path)):
            log("Uploading %s as %s (replacing)" % (os.path.basename(path), types[label].value))
            uploaded = guarded(client.upload_file, dep_id, path, types[label], overwrite=True)
            result["files"].append({"name": os.path.basename(path), "type": types[label].value,
                                    "remote_id": getattr(uploaded, "file_id", None)})
        started = guarded(client.process, dep_id)
        if _is_error_object(started):
            # Do not go on to poll: get_status would hand back the previous
            # processing's "finished" and the failure would pass as success.
            raise OneDepProcessingError(_error_text(started), dep_id=dep_id,
                                        site_url=result["site_url"])
        get_status = lambda: guarded(client.get_status, dep_id)  # noqa: E731
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
        dep_id = guarded(dep.deposit)
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
        get_status = lambda: guarded(dep.get_status)  # noqa: E731

    log("OneDep deposition %s" % dep_id)
    try:
        status = wait_for_processing(get_status, sleep=sleep, log=log)
    except OneDepError as err:
        err.dep_id = dep_id
        err.site_url = result["site_url"]
        raise
    result["status"] = getattr(status, "status", None)
    result["details"] = getattr(status, "details", "") or ""
    if result["status"] != "finished":
        raise OneDepProcessingError(explain_failure(status), dep_id=dep_id,
                                    site_url=result["site_url"])
    return result
