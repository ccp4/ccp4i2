"""Tests for ``LocalSessionAuthMiddleware``.

Covers the desktop-side per-launch-token contract: activation gating
on env-var presence, missing/non-matching token → 401, matching token
→ get-or-create desktop user with the canonical email, fallback to
``DEFAULT_DESKTOP_EMAIL`` when ``CCP4I2_LOCAL_USER_EMAIL`` is unset.

The middleware reads ``CCP4I2_LOCAL_SESSION_TOKEN`` at ``__init__``
time and caches it on ``self.expected_token`` (the env var doesn't
change at runtime; it's set once by the Electron main process when
spawning Django). Tests construct a fresh middleware instance after
each env mutation.
"""

import pytest
from django.http import HttpResponse
from django.test import RequestFactory

from ccp4i2_api.middleware.base import REQUEST_FLAG_ATTR
from ccp4i2_api.middleware.local_session import (
    DEFAULT_DESKTOP_EMAIL,
    LocalSessionAuthMiddleware,
)


@pytest.fixture
def captured_request():
    captured: dict = {}

    def get_response(request):
        captured["request"] = request
        captured["user"] = getattr(request, "user", None)
        captured["flag"] = getattr(request, REQUEST_FLAG_ATTR, False)
        return HttpResponse("ok", status=200)

    return get_response, captured


def _build(monkeypatch, token=None, email=None):
    """Set env to a known state and return a fresh middleware bound
    to the captured-request handler. Returns (middleware, captured).
    """
    if token is None:
        monkeypatch.delenv("CCP4I2_LOCAL_SESSION_TOKEN", raising=False)
    else:
        monkeypatch.setenv("CCP4I2_LOCAL_SESSION_TOKEN", token)
    if email is None:
        monkeypatch.delenv("CCP4I2_LOCAL_USER_EMAIL", raising=False)
    else:
        monkeypatch.setenv("CCP4I2_LOCAL_USER_EMAIL", email)


# ---------------------------------------------------------------------------
# is_active() — activation gating
# ---------------------------------------------------------------------------


def test_is_active_when_token_env_var_set(monkeypatch):
    _build(monkeypatch, token="some-secret")
    mw = LocalSessionAuthMiddleware(lambda r: HttpResponse("ok"))
    assert mw.is_active() is True


def test_is_inactive_when_token_env_var_unset(monkeypatch):
    _build(monkeypatch, token=None)
    mw = LocalSessionAuthMiddleware(lambda r: HttpResponse("ok"))
    assert mw.is_active() is False


def test_is_active_with_empty_string_token(monkeypatch):
    """Edge case: env var set but empty. The middleware activates
    (env var is "present") but every request will fail token comparison
    against the empty string, which is the safe outcome — at most you
    fail closed."""
    monkeypatch.setenv("CCP4I2_LOCAL_SESSION_TOKEN", "")
    mw = LocalSessionAuthMiddleware(lambda r: HttpResponse("ok"))
    assert mw.is_active() is True
    assert mw.expected_token == ""


# ---------------------------------------------------------------------------
# __call__ when inactive → no-op pass-through
# ---------------------------------------------------------------------------


def test_inactive_passes_through_with_no_flag_set(monkeypatch, captured_request):
    _build(monkeypatch, token=None)
    handler, captured = captured_request
    mw = LocalSessionAuthMiddleware(handler)

    response = mw(RequestFactory().get("/api/things"))

    assert response.status_code == 200
    assert captured["flag"] is False


# ---------------------------------------------------------------------------
# authenticate() — failure paths
# ---------------------------------------------------------------------------


def test_no_authorization_header_returns_401(monkeypatch):
    _build(monkeypatch, token="server-secret")
    mw = LocalSessionAuthMiddleware(lambda r: HttpResponse("ok"))

    response = mw(RequestFactory().get("/api/things"))

    assert response.status_code == 401


def test_authorization_header_without_bearer_prefix_returns_401(monkeypatch):
    _build(monkeypatch, token="server-secret")
    mw = LocalSessionAuthMiddleware(lambda r: HttpResponse("ok"))

    response = mw(
        RequestFactory().get(
            "/api/things",
            HTTP_AUTHORIZATION="server-secret",  # no "Bearer " prefix
        )
    )

    assert response.status_code == 401


def test_non_matching_token_returns_401(monkeypatch):
    _build(monkeypatch, token="server-secret")
    mw = LocalSessionAuthMiddleware(lambda r: HttpResponse("ok"))

    response = mw(
        RequestFactory().get(
            "/api/things",
            HTTP_AUTHORIZATION="Bearer different-secret",
        )
    )

    assert response.status_code == 401


# ---------------------------------------------------------------------------
# authenticate() — success path
# ---------------------------------------------------------------------------


@pytest.mark.django_db
def test_matching_token_authenticates_and_sets_flag(
    monkeypatch, captured_request
):
    _build(
        monkeypatch,
        token="exact-match",
        email="martin@ccp4i2.invalid",
    )
    handler, captured = captured_request
    mw = LocalSessionAuthMiddleware(handler)

    response = mw(
        RequestFactory().get(
            "/api/things",
            HTTP_AUTHORIZATION="Bearer exact-match",
        )
    )

    assert response.status_code == 200
    assert captured["flag"] is True
    assert captured["user"].email == "martin@ccp4i2.invalid"


# ---------------------------------------------------------------------------
# Desktop-user identity
# ---------------------------------------------------------------------------


@pytest.mark.django_db
def test_desktop_user_uses_local_user_email_env_var(
    monkeypatch, captured_request
):
    _build(
        monkeypatch,
        token="t",
        email="alice@ccp4i2.invalid",
    )
    handler, captured = captured_request
    mw = LocalSessionAuthMiddleware(handler)

    mw(RequestFactory().get("/api/x", HTTP_AUTHORIZATION="Bearer t"))

    assert captured["user"].email == "alice@ccp4i2.invalid"
    assert captured["user"].username == "alice"


@pytest.mark.django_db
def test_desktop_user_falls_back_to_default_email_when_env_unset(
    monkeypatch, captured_request
):
    """When CCP4I2_LOCAL_USER_EMAIL is unset (e.g., Electron's
    sanitisation collapsed to empty for some reason, or a deployment
    didn't set it), the middleware uses the RFC-6761-reserved
    .invalid TLD fallback. Always well-formed under Django's
    EmailValidator on every platform."""
    _build(monkeypatch, token="t", email=None)
    handler, captured = captured_request
    mw = LocalSessionAuthMiddleware(handler)

    mw(RequestFactory().get("/api/x", HTTP_AUTHORIZATION="Bearer t"))

    assert captured["user"].email == DEFAULT_DESKTOP_EMAIL
    # username defaults to local-part of the fallback email.
    assert captured["user"].username == "desktop"


@pytest.mark.django_db
def test_desktop_user_is_superuser(monkeypatch, captured_request):
    """Desktop deployments are single-user; the OS user running the
    Electron app is the only user the local Django ever sees, and is
    granted superuser to match the trust boundary (you're literally
    running this app on your own machine)."""
    _build(monkeypatch, token="t")
    handler, captured = captured_request
    mw = LocalSessionAuthMiddleware(handler)

    mw(RequestFactory().get("/api/x", HTTP_AUTHORIZATION="Bearer t"))

    user = captured["user"]
    assert user.is_staff is True
    assert user.is_superuser is True


@pytest.mark.django_db
def test_subsequent_requests_reuse_existing_user(monkeypatch, captured_request):
    """get_or_create — first request creates the user, second returns
    the same one. Verifies the migration's worth: not a new row per
    request, just stable identity bound to the email."""
    _build(monkeypatch, token="t", email="bob@ccp4i2.invalid")
    handler, captured = captured_request
    mw = LocalSessionAuthMiddleware(handler)

    mw(RequestFactory().get("/api/a", HTTP_AUTHORIZATION="Bearer t"))
    user_a = captured["user"]

    captured.clear()
    mw(RequestFactory().get("/api/b", HTTP_AUTHORIZATION="Bearer t"))
    user_b = captured["user"]

    assert user_a.pk == user_b.pk
    assert user_a.email == user_b.email == "bob@ccp4i2.invalid"


# ---------------------------------------------------------------------------
# Constant-time comparison — sanity check that wrong tokens of the
# right length still fail.
# ---------------------------------------------------------------------------


def test_same_length_wrong_token_returns_401(monkeypatch):
    """Sanity: hmac.compare_digest doesn't accidentally accept tokens
    of equal length but different content."""
    _build(monkeypatch, token="aaaaaa")
    mw = LocalSessionAuthMiddleware(lambda r: HttpResponse("ok"))

    response = mw(
        RequestFactory().get("/api/x", HTTP_AUTHORIZATION="Bearer bbbbbb")
    )

    assert response.status_code == 401
