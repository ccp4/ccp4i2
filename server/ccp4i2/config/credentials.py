"""
Credential store — where CCP4i2 keeps API tokens, passwords and passphrases.

This is the secret-bearing sibling of :mod:`ccp4i2.config.preferences`. The two
are deliberately separate stores: ``preferences.json`` holds configuration that
is safe to read, copy and paste into a support thread, and must stay that way.

Design (see ``docs/CREDENTIALS_DESIGN.md``) in three sentences:

* **A credential is never a task parameter.** Task parameters are persisted to
  the job container, ``input_params.xml``, the project database and every
  exported project zip. Secrets live here instead, orthogonal to projects.
* **Secrets are write-only from the outside.** Callers inside the job process
  read them; the REST API only ever exposes descriptors and status.
* **Three axes, one store**: *kind* (what fields it has), *scope* (global /
  project / target) and *persistence* (keychain / session / none). ssh and
  OAuth reuse this rather than growing their own systems.

Resolution precedence for every field::

    environment variable
      > session (this process's environ)
        > keychain (or the 0600-file fallback)
          > preferences.json userPreferences   [legacy, read-only]
            > None

Env-first preserves cloud invariance: containers set everything from Key Vault
and neither the keychain nor any file is consulted.

Storage backends, in order of preference:

``keyring``
    The OS secret store (macOS Keychain, Windows Credential Manager, libsecret /
    KWallet). One JSON blob per credential, so a credential costs one keychain
    access rather than one per field.
``file``
    Fallback for headless Linux / containers / CI where no secret service
    exists: mode-0600 JSON at ``~/.ccp4i2-django/credentials.json``. Reported to
    the UI as ``secure=False``. Deliberately *not* encrypted with a local key —
    a key stored beside its own ciphertext is theatre, and claiming encryption
    we are not providing is worse than being honest.
``session``
    This process's ``os.environ``. Detached job subprocesses inherit it, so an
    ad-hoc secret (an ssh password typed at run time) is reachable where it is
    needed, never touches disk, and dies with the app.

This module is pure stdlib plus an optional ``keyring`` import — it must not
import Django, because it is reachable from plugin code that runs CCP4-side.
"""

import json
import logging
import os
import stat
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, Iterable, Optional, Tuple

from .preferences import ccp4i2_home, load_preferences

logger = logging.getLogger(f"ccp4i2.{__name__}")

# The service name credentials are filed under in the OS keychain.
KEYRING_SERVICE = "CCP4i2"

# Name of the 0600 fallback file, kept *separate* from preferences.json so that
# preferences.json remains safe to share.
CREDENTIALS_FILENAME = "credentials.json"


# ---------------------------------------------------------------------------
# Descriptors — what credentials exist, and what shape they have
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class CredentialField:
    """One field of a credential.

    Args:
        name:     key within the credential (e.g. ``token_secret``).
        label:    human label for the input.
        secret:   True if the value must never be displayed or returned. A
                  non-secret field (a token *id*, a username) may be echoed
                  back to the UI; a secret one never is.
        env:      environment variable that overrides this field. Also the
                  ``userPreferences`` key consulted for legacy read-through, so
                  a user who hand-set ``PDB_REDO_TOKEN_ID`` keeps working.
        help:     one-line hint shown under the input.
    """

    name: str
    label: str
    secret: bool
    env: str
    help: str = ""


@dataclass(frozen=True)
class CredentialSpec:
    """Declaration of one credential the application knows how to use.

    Args:
        name:        stable identifier used in URLs and storage keys.
        label:       human name of the service.
        description: one or two sentences shown in the dialog.
        fields:      the fields, in the order they should be presented.
        signup_url:  where the user obtains the credential. Opened externally.
        validator:   optional callable taking the resolved ``{field: value}``
                     dict and returning ``(ok, message)``. Runs server-side.
        kind:        one of ``token_pair``, ``password``, ``key_file``, ``oauth``
                     — currently informational, used by the UI to pick affordances.
        paste_split: True if the dialog should accept a single paste containing
                     every field's value and split it on whitespace/punctuation.
                     Two long opaque strings copied separately is exactly where
                     users fail.
    """

    name: str
    label: str
    description: str
    fields: Tuple[CredentialField, ...]
    signup_url: str = ""
    validator: Optional[Callable[[Dict[str, str]], Tuple[bool, str]]] = None
    kind: str = "token_pair"
    paste_split: bool = False
    # Fields whose *value* may be echoed back as a hint (never secret ones).
    hint_field: str = ""


def _validate_pdb_redo(values: Dict[str, str]) -> Tuple[bool, str]:
    """Probe a PDB-REDO token pair against the live service.

    PDB-REDO publishes no token-management endpoint, but ``GET /api/run`` — the
    caller's own job list — is a signed request that answers exactly the
    question we care about: will this token be accepted? 200 means yes, 401
    means revoked or mistyped.

    Deliberately imports the wrapper's auth helper lazily so this module stays
    importable without the wrapper tree.
    """
    import requests

    from ..wrappers.pdb_redo_api.script.PDBRedoAPIAuth import PDBRedoAPIAuth

    token_id = (values.get("token_id") or "").strip()
    token_secret = (values.get("token_secret") or "").strip()
    if not token_id or not token_secret:
        return False, "Both the token ID and the token secret are required."

    auth = PDBRedoAPIAuth(token_id, token_secret)
    try:
        response = requests.get(
            "https://pdb-redo.eu/api/run", auth=auth, timeout=20
        )
    except requests.RequestException as err:
        # Network failure is NOT a bad token — say so, so the user does not go
        # hunting for a typo that is not there.
        return False, f"Could not reach pdb-redo.eu ({err.__class__.__name__})."

    if response.status_code == 200:
        try:
            runs = response.json()
            count = len(runs) if isinstance(runs, list) else None
        except ValueError:
            count = None
        if count is None:
            return True, "Token accepted by PDB-REDO."
        return True, f"Token accepted by PDB-REDO ({count} run(s) on this account)."

    if response.status_code in (401, 403):
        return False, (
            "PDB-REDO rejected this token. It may have been revoked, or the "
            "token ID and secret may not match."
        )
    # Never log or echo the token; the status code is the whole story.
    return False, f"PDB-REDO returned HTTP {response.status_code}."


PDB_REDO = CredentialSpec(
    name="pdb_redo",
    label="PDB-REDO web service",
    description=(
        "PDB-REDO re-refines and rebuilds your structure on their servers. "
        "Submitting a job needs an API token, which is free. Your token secret "
        "never leaves this machine — it is used to sign requests locally."
    ),
    fields=(
        CredentialField(
            name="token_id",
            label="Token ID",
            secret=False,
            env="PDB_REDO_TOKEN_ID",
            help="The public half of the token pair.",
        ),
        CredentialField(
            name="token_secret",
            label="Token secret",
            secret=True,
            env="PDB_REDO_TOKEN_SECRET",
            help="Used only to sign requests; never transmitted.",
        ),
    ),
    signup_url="https://pdb-redo.eu/token",
    validator=_validate_pdb_redo,
    kind="token_pair",
    paste_split=True,
    hint_field="token_id",
)


# The registry. Adding a credential = adding a CredentialSpec here; the API and
# the whole frontend are descriptor-driven and need no change.
CREDENTIALS: Dict[str, CredentialSpec] = {
    PDB_REDO.name: PDB_REDO,
}


def get_spec(name: str) -> Optional[CredentialSpec]:
    """Look up a credential declaration by name, or None if unknown."""
    return CREDENTIALS.get(name)


# ---------------------------------------------------------------------------
# Storage keys — (name, scope, scope_id)
# ---------------------------------------------------------------------------

SCOPE_GLOBAL = "global"
SCOPE_PROJECT = "project"
SCOPE_TARGET = "target"

PERSIST_KEYCHAIN = "keychain"
PERSIST_SESSION = "session"
PERSIST_NONE = "none"


def _storage_key(name: str, scope: str = SCOPE_GLOBAL, scope_id: str = "") -> str:
    """Build the account key a credential is filed under.

    Scoping is in the key from day one so that project-scoped tokens (a token
    per CCP4i2 project, as PDB-REDO suggested, keeping their job list legible)
    become a UI affordance later rather than a storage migration.
    """
    return f"{name}:{scope}:{scope_id or '-'}"


def _session_env_key(storage_key: str, field_name: str) -> str:
    """Environment variable name used for a session-persisted field."""
    safe = storage_key.replace(":", "__").replace("-", "_").upper()
    return f"CCP4I2_CRED__{safe}__{field_name.upper()}"


# ---------------------------------------------------------------------------
# Backends
# ---------------------------------------------------------------------------


def _keyring():
    """Return a usable ``keyring`` module, or None.

    None covers three distinct situations that all mean the same thing to us:
    the package is not installed, no backend could be selected, or the selected
    backend is one of keyring's null/fail placeholders.
    """
    try:
        import keyring
        from keyring.backends import fail as keyring_fail
    except Exception:
        return None
    try:
        backend = keyring.get_keyring()
    except Exception:
        return None
    if backend is None or isinstance(backend, keyring_fail.Keyring):
        return None
    # keyring.backends.null.Keyring silently discards writes — worse than no
    # backend at all, because it looks like success.
    if type(backend).__module__.endswith("backends.null"):
        return None
    return keyring


def keyring_available() -> bool:
    """True if an OS secret store is usable on this machine."""
    return _keyring() is not None


def _credentials_file() -> Path:
    return ccp4i2_home() / CREDENTIALS_FILENAME


def _load_file_store() -> dict:
    """Read the 0600 fallback file; ``{}`` if absent or unreadable."""
    try:
        with open(_credentials_file(), encoding="utf-8") as handle:
            data = json.load(handle)
    except (FileNotFoundError, NotADirectoryError, json.JSONDecodeError, OSError):
        return {}
    return data if isinstance(data, dict) else {}


def _save_file_store(store: dict) -> None:
    """Write the fallback file with owner-only permissions.

    The file is created with mode 0600 *before* anything is written to it, so
    there is no window in which a secret sits on disk world-readable.
    """
    home = ccp4i2_home()
    home.mkdir(parents=True, exist_ok=True)
    path = _credentials_file()
    flags = os.O_WRONLY | os.O_CREAT | os.O_TRUNC
    fd = os.open(str(path), flags, 0o600)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(store, handle, indent=2, sort_keys=True)
            handle.write("\n")
    finally:
        try:
            os.chmod(str(path), stat.S_IRUSR | stat.S_IWUSR)
        except OSError:
            logger.debug("could not chmod %s", path, exc_info=True)


# ---------------------------------------------------------------------------
# Read / write
# ---------------------------------------------------------------------------


@dataclass
class Resolved:
    """A resolved credential: its values plus where they came from."""

    values: Dict[str, str] = field(default_factory=dict)
    source: str = "unset"  # environment | session | keychain | file | preferences
    complete: bool = False


def _read_stored(storage_key: str) -> Tuple[Dict[str, str], str]:
    """Read the persisted blob for a storage key from keychain, then file."""
    kr = _keyring()
    if kr is not None:
        try:
            blob = kr.get_password(KEYRING_SERVICE, storage_key)
        except Exception:
            logger.debug("keyring read failed for %s", storage_key, exc_info=True)
            blob = None
        if blob:
            try:
                values = json.loads(blob)
                if isinstance(values, dict):
                    return {k: str(v) for k, v in values.items()}, "keychain"
            except ValueError:
                logger.debug("corrupt keyring blob for %s", storage_key)
    values = _load_file_store().get(storage_key)
    if isinstance(values, dict):
        return {k: str(v) for k, v in values.items()}, "file"
    return {}, "unset"


def resolve_credential(
    name: str, scope: str = SCOPE_GLOBAL, scope_id: str = ""
) -> Resolved:
    """Resolve every field of a credential, honouring the precedence chain.

    Precedence is applied *per field*, so a deployment can override just the
    secret from Key Vault while the id comes from the keychain.

    Args:
        name:     credential name (a key of :data:`CREDENTIALS`).
        scope:    ``global``, ``project`` or ``target``.
        scope_id: project id / host name, when scope is not global.

    Returns:
        Resolved: ``values`` maps field name to value for every field that
        resolved to something; ``complete`` is True only when *all* fields did.
        ``source`` names the highest-precedence layer that contributed.

    Example:
        >>> resolve_credential("pdb_redo").complete
        False
    """
    spec = get_spec(name)
    if spec is None:
        return Resolved()

    key = _storage_key(name, scope, scope_id)
    stored, stored_source = _read_stored(key)
    legacy = load_preferences().get("userPreferences", {})
    legacy = legacy if isinstance(legacy, dict) else {}

    values: Dict[str, str] = {}
    sources = set()
    for f in spec.fields:
        env_val = os.environ.get(f.env)
        if env_val:
            values[f.name] = env_val
            sources.add("environment")
            continue
        session_val = os.environ.get(_session_env_key(key, f.name))
        if session_val:
            values[f.name] = session_val
            sources.add("session")
            continue
        if stored.get(f.name):
            values[f.name] = stored[f.name]
            sources.add(stored_source)
            continue
        legacy_val = legacy.get(f.env)
        if legacy_val:
            values[f.name] = str(legacy_val)
            sources.add("preferences")

    # Report the highest-precedence contributing layer.
    for candidate in ("environment", "session", "keychain", "file", "preferences"):
        if candidate in sources:
            source = candidate
            break
    else:
        source = "unset"

    return Resolved(
        values=values,
        source=source,
        complete=all(values.get(f.name) for f in spec.fields),
    )


def get_credential(
    name: str, scope: str = SCOPE_GLOBAL, scope_id: str = ""
) -> Optional[Dict[str, str]]:
    """Return a *complete* credential's values, or None.

    This is the accessor plugin code should use: it returns None rather than a
    partially-filled dict, so a task can never accidentally sign a request with
    a missing half of a token pair (the ``str(None)`` -> ``"None"`` failure mode
    this store was written to remove).

    Example:
        >>> creds = get_credential("pdb_redo")
        >>> token_id = creds["token_id"] if creds else None
    """
    resolved = resolve_credential(name, scope, scope_id)
    return dict(resolved.values) if resolved.complete else None


def set_credential(
    name: str,
    values: Dict[str, str],
    scope: str = SCOPE_GLOBAL,
    scope_id: str = "",
    persistence: str = PERSIST_KEYCHAIN,
) -> str:
    """Store a credential. Returns the layer it was written to.

    Args:
        name:        credential name.
        values:      ``{field_name: value}``. Fields absent from the spec are
                     dropped; empty values are dropped.
        scope:       ``global``, ``project`` or ``target``.
        scope_id:    project id / host name, when scope is not global.
        persistence: ``keychain`` (durable), ``session`` (this process and the
                     job subprocesses it spawns, until the app quits) or
                     ``none`` (a no-op; the caller passes the value straight to
                     the transport).

    Returns:
        str: ``keychain``, ``file``, ``session`` or ``none``.

    Raises:
        KeyError: if ``name`` is not a registered credential.
    """
    spec = get_spec(name)
    if spec is None:
        raise KeyError(name)

    allowed = {f.name for f in spec.fields}
    clean = {
        k: str(v).strip()
        for k, v in values.items()
        if k in allowed and v is not None and str(v).strip()
    }
    key = _storage_key(name, scope, scope_id)

    if persistence == PERSIST_NONE:
        return PERSIST_NONE

    if persistence == PERSIST_SESSION:
        # In this process's environment, which detached job subprocesses
        # inherit. Never touches disk; dies with the process.
        for field_name, value in clean.items():
            os.environ[_session_env_key(key, field_name)] = value
        return PERSIST_SESSION

    kr = _keyring()
    if kr is not None:
        try:
            kr.set_password(KEYRING_SERVICE, key, json.dumps(clean))
            # A value that had previously landed in the 0600 file (or in
            # legacy preferences.json) would otherwise shadow nothing but
            # confuse `source` reporting — drop it now that the keychain has it.
            _forget_file(key)
            _forget_legacy(spec)
            return "keychain"
        except Exception:
            logger.warning(
                "OS keychain write failed for %s; falling back to file store",
                name,
                exc_info=True,
            )

    store = _load_file_store()
    store[key] = clean
    _save_file_store(store)
    _forget_legacy(spec)
    return "file"


def _forget_file(key: str) -> None:
    store = _load_file_store()
    if key in store:
        del store[key]
        _save_file_store(store)


def _forget_legacy(spec: CredentialSpec) -> None:
    """Remove hand-set values from ``preferences.json``.

    Called after a successful write to a real store: leaving a plaintext copy
    behind in a file users are asked to share would defeat the whole exercise.
    Best-effort — never raises into the caller.
    """
    try:
        from .preferences import save_preferences

        prefs = load_preferences()
        bag = prefs.get("userPreferences")
        if not isinstance(bag, dict):
            return
        removed = [f.env for f in spec.fields if f.env in bag]
        if not removed:
            return
        for env_name in removed:
            bag.pop(env_name, None)
        prefs["userPreferences"] = bag
        save_preferences(prefs)
        logger.info(
            "migrated %s out of preferences.json into the credential store",
            spec.name,
        )
    except Exception:
        logger.debug("legacy credential cleanup skipped", exc_info=True)


def clear_credential(
    name: str, scope: str = SCOPE_GLOBAL, scope_id: str = ""
) -> None:
    """Remove a credential from every writable layer.

    Environment-variable values are left alone — they are not ours to remove,
    and in cloud they are the point.
    """
    spec = get_spec(name)
    key = _storage_key(name, scope, scope_id)

    kr = _keyring()
    if kr is not None:
        try:
            kr.delete_password(KEYRING_SERVICE, key)
        except Exception:
            logger.debug("keyring delete skipped for %s", key, exc_info=True)

    _forget_file(key)

    for field_name in [f.name for f in spec.fields] if spec else []:
        os.environ.pop(_session_env_key(key, field_name), None)

    if spec is not None:
        _forget_legacy(spec)


def validate_credential(
    name: str,
    values: Optional[Dict[str, str]] = None,
    scope: str = SCOPE_GLOBAL,
    scope_id: str = "",
) -> Tuple[bool, str]:
    """Probe a credential against its live service.

    Args:
        name:   credential name.
        values: values to test. If omitted, the stored/resolved credential is
                tested — which is how the *Test* button re-checks something the
                user cannot see.

    Returns:
        (ok, message): ``message`` is user-facing and must never contain secret
        material.
    """
    spec = get_spec(name)
    if spec is None:
        return False, f"Unknown credential '{name}'."
    if values is None:
        resolved = resolve_credential(name, scope, scope_id)
        if not resolved.complete:
            return False, "No credential is configured."
        values = resolved.values
    if spec.validator is None:
        return True, "No validation available for this credential."
    try:
        return spec.validator(values)
    except Exception as err:  # never let a probe break the caller
        logger.debug("credential validation raised", exc_info=True)
        return False, f"Validation failed ({err.__class__.__name__})."


# ---------------------------------------------------------------------------
# Descriptors for the API / UI
# ---------------------------------------------------------------------------

# Cached outcome of the last validation, per storage key. Purely informational:
# it lets the UI say "checked 2 minutes ago" without re-probing on every render.
_validation_cache: Dict[str, dict] = {}


def note_validation(
    name: str, ok: bool, message: str, scope: str = SCOPE_GLOBAL, scope_id: str = ""
) -> None:
    """Record the outcome of a validation probe for display."""
    _validation_cache[_storage_key(name, scope, scope_id)] = {
        "ok": bool(ok),
        "message": message,
        "at": time.time(),
    }


def describe_credential(
    name: str, scope: str = SCOPE_GLOBAL, scope_id: str = "", editable: bool = True
) -> Optional[dict]:
    """Build the JSON descriptor the API returns for one credential.

    Contains everything the UI needs to render a form and a status line, and
    **no secret material** — not even for a loopback client. The only value
    echoed back is a hint derived from a non-secret field (the last four
    characters of a token id), enough for a user to tell two tokens apart.
    """
    spec = get_spec(name)
    if spec is None:
        return None

    resolved = resolve_credential(name, scope, scope_id)
    hint = ""
    if spec.hint_field:
        hint_value = resolved.values.get(spec.hint_field, "")
        if hint_value:
            hint = hint_value[-4:] if len(hint_value) > 4 else hint_value

    cached = _validation_cache.get(_storage_key(name, scope, scope_id))

    return {
        "name": spec.name,
        "label": spec.label,
        "description": spec.description,
        "kind": spec.kind,
        "signupUrl": spec.signup_url,
        "pasteSplit": spec.paste_split,
        "canValidate": spec.validator is not None,
        "fields": [
            {
                "name": f.name,
                "label": f.label,
                "secret": f.secret,
                "help": f.help,
                "isSet": bool(resolved.values.get(f.name)),
            }
            for f in spec.fields
        ],
        "isSet": resolved.complete,
        "source": resolved.source,
        "secure": keyring_available(),
        "storageLabel": storage_label(),
        # Env-supplied credentials are not ours to edit — in cloud they are the
        # deployment's configuration, not the user's.
        "editable": bool(editable) and resolved.source != "environment",
        "hint": hint,
        "lastValidated": cached["at"] if cached else None,
        "lastValidationOk": cached["ok"] if cached else None,
        "lastValidationMessage": cached["message"] if cached else "",
    }


def storage_label() -> str:
    """A plain-language statement of where a secret will actually be kept.

    Shown verbatim in the dialog. Being specific (and honest about the fallback)
    is the point: users deserve to know whether the platform is really
    protecting the value.
    """
    import sys

    if keyring_available():
        if sys.platform == "darwin":
            return "your macOS Keychain"
        if sys.platform.startswith("win"):
            return "the Windows Credential Manager"
        return "your desktop secret store (libsecret/KWallet)"
    return (
        "a file readable only by you (%s) - this system provides no secret store"
        % _credentials_file()
    )


def describe_all(editable: bool = True) -> Iterable[dict]:
    """Descriptors for every registered credential."""
    return [describe_credential(name, editable=editable) for name in CREDENTIALS]
