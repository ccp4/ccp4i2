"""
REST surface for the credential store — the single choke point through which
API tokens, passwords and passphrases are set, tested and cleared.

Mirrors the ``config/program-preferences/`` endpoints in :mod:`.views`,
including the ``editable`` flag that keeps cloud deployments honest, with two
extra rules that the secret-bearing case demands:

* **Write-only secrets.** ``set/`` takes values in; *nothing* here ever returns
  secret material, not even to a loopback client. A user proves a credential is
  correct by pressing *Test*, never by reading it back.
* **Desktop gate.** Writes are refused except on a desktop launch. Accepting a
  secret over a shared Django would route one user's credentials for a
  third-party service through a multi-tenant host — the same boundary
  ``docs/REMOTE_JOB_EXECUTION_PLAN.md`` draws for ssh.

See ``docs/CREDENTIALS_DESIGN.md``.
"""

import logging
import os

from django.http import JsonResponse
from rest_framework.decorators import api_view

from ..config import credentials as cred
from ..config.preferences import is_desktop

logger = logging.getLogger(__name__)

# Opt-in for a trusted single-user server deployment (someone running the
# server for themselves on a machine they control) that nonetheless wants the
# credential UI. Mirrors the CCP4I2_ALLOW_INTERACTIVE_SSH opt-in in the
# remote-execution plan: off by default, and a deliberate act to turn on.
_WRITE_OPT_IN = "CCP4I2_ALLOW_CREDENTIAL_WRITE"


def _can_write(request) -> bool:
    """True when this deployment may accept a secret over its own API.

    The signal is ``CCP4I2_LOCAL_SESSION_TOKEN`` — set only by the desktop
    launcher, and the same switch ``settings.py`` uses to select local-session
    auth and ``set_program_preferences`` uses to gate preference writes.

    Peer address is deliberately *not* consulted. Behind a reverse proxy or a
    sidecar, ``REMOTE_ADDR`` is the proxy's address and can perfectly well be
    ``127.0.0.1`` for a request that arrived from the far side of the internet;
    treating that as "this machine" would open credential writes on exactly the
    multi-tenant deployments this gate exists to protect.
    """
    if is_desktop():
        return True
    return os.environ.get(_WRITE_OPT_IN, "").lower() in ("true", "1", "yes")


def _forbidden_write():
    return JsonResponse(
        {
            "success": False,
            "error": (
                "Credentials can only be set from the local application. In a "
                "server deployment they are supplied as environment variables "
                "or per-user secrets."
            ),
        },
        status=409,
    )


@api_view(["GET"])
def list_credentials(request):
    """Describe every registered credential and its current status.

    GET /api/ccp4i2/config/credentials/

    Response::

        {"success": true, "data": {"credentials": [ <descriptor>, ... ],
                                   "editable": <bool>, "secure": <bool>,
                                   "storageLabel": "your macOS Keychain"}}

    A descriptor carries the field *shape*, whether each field is set, where the
    value came from, and the last validation outcome — never a secret. The only
    value echoed back is a short hint from a non-secret field, so a user can
    tell two tokens apart.
    """
    editable = _can_write(request)
    return JsonResponse(
        {
            "success": True,
            "data": {
                "credentials": list(cred.describe_all(editable=editable)),
                "editable": editable,
                "secure": cred.keyring_available(),
                "storageLabel": cred.storage_label(),
            },
        }
    )


@api_view(["POST", "PATCH"])
def set_credential_view(request, name):
    """Store a credential.

    POST /api/ccp4i2/config/credentials/{name}/set/
        {"values": {"token_id": "...", "token_secret": "..."},
         "persistence": "keychain" | "session" | "none",
         "scope": "global", "scopeId": ""}

    Returns the refreshed descriptor and the layer actually written to
    (``keychain``, ``file``, ``session`` or ``none``) — the UI reports this
    verbatim so a fallback to the 0600 file is never silent.
    """
    if not _can_write(request):
        return _forbidden_write()
    if cred.get_spec(name) is None:
        return JsonResponse(
            {"success": False, "error": f"Unknown credential '{name}'."}, status=404
        )

    payload = request.data if isinstance(request.data, dict) else {}
    values = payload.get("values")
    if not isinstance(values, dict) or not values:
        return JsonResponse(
            {"success": False, "error": "No values supplied."}, status=400
        )

    persistence = payload.get("persistence") or cred.PERSIST_KEYCHAIN
    if persistence not in (
        cred.PERSIST_KEYCHAIN,
        cred.PERSIST_SESSION,
        cred.PERSIST_NONE,
    ):
        return JsonResponse(
            {"success": False, "error": f"Unknown persistence '{persistence}'."},
            status=400,
        )

    scope = payload.get("scope") or cred.SCOPE_GLOBAL
    scope_id = str(payload.get("scopeId") or "")

    written = cred.set_credential(
        name, values, scope=scope, scope_id=scope_id, persistence=persistence
    )
    # Deliberately not logged with values, and not at INFO with the name alone
    # being enough to be useful.
    logger.info("credential '%s' stored (%s)", name, written)

    return JsonResponse(
        {
            "success": True,
            "data": {
                "storedIn": written,
                "credential": cred.describe_credential(
                    name, scope=scope, scope_id=scope_id, editable=True
                ),
            },
        }
    )


@api_view(["POST"])
def validate_credential_view(request, name):
    """Probe a credential against its live service.

    POST /api/ccp4i2/config/credentials/{name}/validate/
        {"values": {...}}   # optional: test these instead of the stored ones

    With no ``values`` the *stored* credential is tested, which is how the UI
    re-checks something the user cannot see. With ``values`` the user can test
    before saving — the point of Test-before-Save is that a user cannot
    otherwise distinguish a mistyped secret from a service outage.
    """
    if cred.get_spec(name) is None:
        return JsonResponse(
            {"success": False, "error": f"Unknown credential '{name}'."}, status=404
        )

    payload = request.data if isinstance(request.data, dict) else {}
    values = payload.get("values")
    if values is not None and not isinstance(values, dict):
        return JsonResponse(
            {"success": False, "error": "values must be an object."}, status=400
        )
    if values is not None and not _can_write(request):
        # Testing arbitrary values means accepting a secret over the wire.
        return _forbidden_write()

    scope = payload.get("scope") or cred.SCOPE_GLOBAL
    scope_id = str(payload.get("scopeId") or "")

    ok, message = cred.validate_credential(
        name, values=values or None, scope=scope, scope_id=scope_id
    )
    # Only cache the outcome when it describes the *stored* credential.
    if values is None:
        cred.note_validation(name, ok, message, scope=scope, scope_id=scope_id)

    return JsonResponse({"success": True, "data": {"ok": ok, "message": message}})


@api_view(["POST", "DELETE"])
def clear_credential_view(request, name):
    """Remove a credential from every writable layer.

    POST|DELETE /api/ccp4i2/config/credentials/{name}/clear/

    Environment-supplied values are left alone: they belong to the deployment,
    not to the user, and in cloud they are the whole configuration mechanism.
    """
    if not _can_write(request):
        return _forbidden_write()
    if cred.get_spec(name) is None:
        return JsonResponse(
            {"success": False, "error": f"Unknown credential '{name}'."}, status=404
        )

    payload = request.data if isinstance(request.data, dict) else {}
    scope = payload.get("scope") or cred.SCOPE_GLOBAL
    scope_id = str(payload.get("scopeId") or "")

    cred.clear_credential(name, scope=scope, scope_id=scope_id)
    logger.info("credential '%s' cleared", name)
    return JsonResponse(
        {
            "success": True,
            "data": {
                "credential": cred.describe_credential(
                    name, scope=scope, scope_id=scope_id, editable=True
                )
            },
        }
    )
