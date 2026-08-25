"""Scoped, short-lived read grants for browser-issued file requests.

Every CCP4i2 API path is authenticated by the middleware in this package,
which expects an ``Authorization: Bearer`` header. Some requests can never
carry one: when a task's HTML report is opened in a browser tab (ProSMART,
PDB-REDO, MrParse, xia2), the page's own ``<img>``, ``<link>``, ``fetch()``
and relative-link requests are issued by the browser itself, with no
opportunity to attach a header. Those requests currently 401.

A *file grant* is a capability, not an identity: a signed, expiring token
that authorises **read-only** access to one directory subtree of the
path-based file-serving endpoint, on behalf of the user who minted it. It
is deliberately not a credential — it carries no bearer token, so putting
it in a URL or a path-scoped cookie does not expose the caller's identity
token (an Azure AD access token would also blow the 4KB cookie limit once
group claims are included).

Signing key
-----------
Grants are signed with the per-launch desktop session secret when one is
present, falling back to ``SECRET_KEY``. On the desktop that matters:
``SECRET_KEY`` there is a well-known hardcoded development default, so
signing with it would let any local process mint grants and read project
files. Binding to ``CCP4I2_LOCAL_SESSION_TOKEN`` — random per launch —
keeps grants as unguessable as the session itself, and expires every
outstanding grant when the app restarts.
"""

import os
import posixpath
from typing import Optional

from django.conf import settings
from django.core import signing

# Signature namespace — keeps grants from being confused with any other
# signed value produced from the same key.
GRANT_SALT = "ccp4i2.file-grant"

# How long a grant stays valid. Long enough to browse a report, short
# enough that a leaked URL is not a lasting capability.
DEFAULT_TTL_SECONDS = 3600

# Where a grant may arrive. The header is what the Next proxy forwards;
# the query parameter seeds the initial navigation; the cookie (set by the
# proxy, scoped to the report's directory) carries it on the page's own
# subsequent relative requests.
GRANT_HEADER = "HTTP_X_CCP4I2_FILE_GRANT"
GRANT_QUERY_PARAM = "file_grant"
GRANT_COOKIE = "ccp4i2_file_grant"

# Grants authorise reads and nothing else.
SAFE_METHODS = ("GET", "HEAD")

# ...and only ever the path-based file-serving endpoint. A grant can never
# be replayed against the REST API proper, whatever prefix it claims.
REQUIRED_PATH_SEGMENT = "/files_by_path/"


def grant_ttl() -> int:
    """Grant lifetime in seconds (``CCP4I2_FILE_GRANT_TTL`` overrides)."""
    return int(getattr(settings, "CCP4I2_FILE_GRANT_TTL", DEFAULT_TTL_SECONDS))


def _signing_key() -> str:
    return os.environ.get("CCP4I2_LOCAL_SESSION_TOKEN") or settings.SECRET_KEY


def normalise_prefix(path: str) -> str:
    """Canonicalise a URL path to a directory prefix with a trailing slash.

    Collapses ``.``/``..`` segments so a grant cannot be minted for — or
    checked against — a prefix that traverses out of its subtree.
    """
    if not path.startswith("/"):
        path = "/" + path
    return posixpath.normpath(path).rstrip("/") + "/"


def mint_grant(*, user_pk, path_prefix: str) -> str:
    """Sign a grant authorising reads under ``path_prefix``.

    ``path_prefix`` is a URL path on the Django side (the caller should
    build it with ``reverse()``), not a filesystem path.
    """
    prefix = normalise_prefix(path_prefix)
    if REQUIRED_PATH_SEGMENT not in prefix:
        raise ValueError(
            f"file grants may only be scoped to {REQUIRED_PATH_SEGMENT} paths"
        )
    return signing.dumps(
        {"u": str(user_pk), "s": prefix},
        key=_signing_key(),
        salt=GRANT_SALT,
    )


def extract_grant(request) -> Optional[str]:
    """Pull a grant off a request: header, then query parameter, then cookie."""
    return (
        request.META.get(GRANT_HEADER)
        or request.GET.get(GRANT_QUERY_PARAM)
        or request.COOKIES.get(GRANT_COOKIE)
    )


def grant_user_pk(token: str, *, path: str, method: str) -> Optional[str]:
    """Validate a grant against a request, returning the minting user's pk.

    Returns ``None`` — never raises — if the grant is unsigned, tampered
    with, expired, presented for an unsafe method, or scoped to a subtree
    that does not contain ``path``. Callers treat ``None`` as "no grant"
    and fall through to normal authentication.
    """
    if method not in SAFE_METHODS:
        return None
    try:
        payload = signing.loads(
            token,
            key=_signing_key(),
            salt=GRANT_SALT,
            max_age=grant_ttl(),
        )
    except signing.BadSignature:
        return None
    if not isinstance(payload, dict):
        return None
    prefix = payload.get("s")
    user_pk = payload.get("u")
    if not prefix or not user_pk or REQUIRED_PATH_SEGMENT not in prefix:
        return None
    if not _within(path, prefix):
        return None
    return user_pk


def _within(path: str, prefix: str) -> bool:
    """True if ``path`` lies inside the ``prefix`` subtree.

    Both sides are normalised first, so ``/a/b/../../etc`` cannot smuggle
    its way past a ``/a/`` prefix.
    """
    normalised_path = normalise_prefix(path)
    normalised_prefix = normalise_prefix(prefix)
    return normalised_path.startswith(normalised_prefix)
