"""Tests for ``AzureADAuthMiddleware``.

Covers the full request-lifecycle contract: activation gating, exempt-
path bypass, missing/invalid bearer → 401, misconfigured server → 401,
group authorization (pass / deny / claims-overage) → 200/403, user
get-or-create from JWT claims, X-User-Email header fallback when the
JWT lacks email claims.

JWT validation itself is not tested here — the middleware delegates
to ``AzureADTokenValidator`` which involves Azure AD's JWKS endpoint
(network) and PyJWT internals. Those are mocked at the
``get_validator`` boundary so the tests stay fast and offline.
"""

from unittest.mock import MagicMock, patch

import pytest
from django.http import HttpResponse
from django.test import RequestFactory

import ccp4i2_api.middleware.azure_ad as azure_ad_module
from ccp4i2_api.middleware.azure_ad import AzureADAuthMiddleware
from ccp4i2_api.middleware.base import REQUEST_FLAG_ATTR


@pytest.fixture(autouse=True)
def reset_validator_singleton():
    """Reset ``_validator`` between tests so env-var changes are seen
    by ``get_validator``'s cache."""
    azure_ad_module._validator = None
    yield
    azure_ad_module._validator = None


@pytest.fixture
def captured_request():
    captured: dict = {}

    def get_response(request):
        captured["request"] = request
        captured["user"] = getattr(request, "user", None)
        captured["flag"] = getattr(request, REQUEST_FLAG_ATTR, False)
        return HttpResponse("ok", status=200)

    return get_response, captured


def _claims(**overrides):
    """Default valid-token claims; tests override fields they care about."""
    base = {
        "sub": "abcd1234efgh5678ijkl9012mnop3456",
        "email": "alice@example.com",
        "given_name": "Alice",
        "family_name": "Smith",
    }
    base.update(overrides)
    return base


# ---------------------------------------------------------------------------
# is_active() — activation gating
# ---------------------------------------------------------------------------


def test_is_active_when_require_auth_true(monkeypatch):
    monkeypatch.setenv("CCP4I2_REQUIRE_AUTH", "true")
    mw = AzureADAuthMiddleware(lambda r: HttpResponse("ok"))
    assert mw.is_active() is True


def test_is_inactive_when_require_auth_unset(monkeypatch):
    monkeypatch.delenv("CCP4I2_REQUIRE_AUTH", raising=False)
    mw = AzureADAuthMiddleware(lambda r: HttpResponse("ok"))
    assert mw.is_active() is False


def test_is_inactive_when_require_auth_explicitly_false(monkeypatch):
    monkeypatch.setenv("CCP4I2_REQUIRE_AUTH", "false")
    mw = AzureADAuthMiddleware(lambda r: HttpResponse("ok"))
    assert mw.is_active() is False


# ---------------------------------------------------------------------------
# __call__ inactive → no-op pass-through (no flag, no user touch)
# ---------------------------------------------------------------------------


def test_inactive_passes_through_with_no_flag_set(monkeypatch, captured_request):
    monkeypatch.delenv("CCP4I2_REQUIRE_AUTH", raising=False)
    handler, captured = captured_request
    mw = AzureADAuthMiddleware(handler)

    response = mw(RequestFactory().get("/api/things"))

    assert response.status_code == 200
    assert captured["flag"] is False


# ---------------------------------------------------------------------------
# Exempt paths — active middleware skips auth for health/version checks
# ---------------------------------------------------------------------------


def test_active_exempt_path_skips_auth(monkeypatch, captured_request):
    monkeypatch.setenv("CCP4I2_REQUIRE_AUTH", "true")
    handler, captured = captured_request
    mw = AzureADAuthMiddleware(handler)

    response = mw(RequestFactory().get("/health"))

    assert response.status_code == 200
    assert captured["flag"] is False  # exempt = no flag


def test_active_exempt_path_with_trailing_slash(monkeypatch, captured_request):
    monkeypatch.setenv("CCP4I2_REQUIRE_AUTH", "true")
    handler, captured = captured_request
    mw = AzureADAuthMiddleware(handler)

    response = mw(RequestFactory().get("/api/health/"))

    assert response.status_code == 200
    assert captured["flag"] is False


# ---------------------------------------------------------------------------
# authenticate() — failure paths
# ---------------------------------------------------------------------------


def test_no_authorization_header_returns_401(monkeypatch):
    monkeypatch.setenv("CCP4I2_REQUIRE_AUTH", "true")
    mw = AzureADAuthMiddleware(lambda r: HttpResponse("ok"))

    response = mw(RequestFactory().get("/api/things"))

    assert response.status_code == 401


def test_invalid_token_returns_401(monkeypatch):
    monkeypatch.setenv("CCP4I2_REQUIRE_AUTH", "true")
    monkeypatch.setenv("AZURE_AD_TENANT_ID", "tenant-id")
    monkeypatch.setenv("AZURE_AD_CLIENT_ID", "client-id")

    with patch.object(azure_ad_module, "get_validator") as mock_get_v:
        mock_v = MagicMock()
        mock_v.validate_token.return_value = (False, None, "Token has expired")
        mock_get_v.return_value = mock_v

        mw = AzureADAuthMiddleware(lambda r: HttpResponse("ok"))
        request = RequestFactory().get(
            "/api/things",
            HTTP_AUTHORIZATION="Bearer expired-token",
        )
        response = mw(request)

    assert response.status_code == 401


def test_misconfigured_validator_returns_401(monkeypatch):
    """Validator unavailable (tenant/client env vars missing). Caller
    treats this the same as an unauthenticated request — 401, not 500
    — because the user-visible state is identical: 'we can't auth you'.
    The operator-facing log line distinguishes."""
    monkeypatch.setenv("CCP4I2_REQUIRE_AUTH", "true")
    # AZURE_AD_TENANT_ID and AZURE_AD_CLIENT_ID intentionally not set.

    with patch.object(azure_ad_module, "get_validator") as mock_get_v:
        mock_get_v.return_value = None

        mw = AzureADAuthMiddleware(lambda r: HttpResponse("ok"))
        request = RequestFactory().get(
            "/api/things",
            HTTP_AUTHORIZATION="Bearer some-token",
        )
        response = mw(request)

    assert response.status_code == 401


# ---------------------------------------------------------------------------
# Group authorization — passes / denies / overage 403
# ---------------------------------------------------------------------------


@pytest.mark.django_db
def test_group_authorization_passes(monkeypatch, captured_request):
    monkeypatch.setenv("CCP4I2_REQUIRE_AUTH", "true")
    monkeypatch.setenv("AZURE_AD_TENANT_ID", "tenant-id")
    monkeypatch.setenv("AZURE_AD_CLIENT_ID", "client-id")
    monkeypatch.setenv("ALLOWED_AZURE_AD_GROUPS", "allowed-group-id-1,other")

    handler, captured = captured_request

    with patch.object(azure_ad_module, "get_validator") as mock_get_v:
        mock_v = MagicMock()
        mock_v.validate_token.return_value = (
            True,
            _claims(groups=["allowed-group-id-1"]),
            None,
        )
        mock_get_v.return_value = mock_v

        mw = AzureADAuthMiddleware(handler)
        request = RequestFactory().get(
            "/api/things",
            HTTP_AUTHORIZATION="Bearer good-token",
        )
        response = mw(request)

    assert response.status_code == 200
    assert captured["flag"] is True
    assert captured["user"].email == "alice@example.com"


def test_group_authorization_denies_when_user_not_in_allowed(monkeypatch):
    monkeypatch.setenv("CCP4I2_REQUIRE_AUTH", "true")
    monkeypatch.setenv("AZURE_AD_TENANT_ID", "tenant-id")
    monkeypatch.setenv("AZURE_AD_CLIENT_ID", "client-id")
    monkeypatch.setenv("ALLOWED_AZURE_AD_GROUPS", "allowed-group-id-1")

    with patch.object(azure_ad_module, "get_validator") as mock_get_v:
        mock_v = MagicMock()
        mock_v.validate_token.return_value = (
            True,
            _claims(groups=["unrelated-group-id"]),
            None,
        )
        mock_get_v.return_value = mock_v

        mw = AzureADAuthMiddleware(lambda r: HttpResponse("ok"))
        request = RequestFactory().get(
            "/api/things",
            HTTP_AUTHORIZATION="Bearer wrong-group-token",
        )
        response = mw(request)

    assert response.status_code == 403


def test_group_claims_overage_returns_403(monkeypatch):
    """Azure AD substitutes _claim_names/_claim_sources for the full
    groups list when a user is in >200 groups. Without the full list
    we can't prove membership; surface a 403 with the canonical message
    pointing at the dedicated-app-access-group remediation."""
    monkeypatch.setenv("CCP4I2_REQUIRE_AUTH", "true")
    monkeypatch.setenv("AZURE_AD_TENANT_ID", "tenant-id")
    monkeypatch.setenv("AZURE_AD_CLIENT_ID", "client-id")
    monkeypatch.setenv("ALLOWED_AZURE_AD_GROUPS", "allowed-group-id-1")

    with patch.object(azure_ad_module, "get_validator") as mock_get_v:
        mock_v = MagicMock()
        mock_v.validate_token.return_value = (
            True,
            _claims(
                _claim_names={"groups": "src1"},
                _claim_sources={"src1": {"endpoint": "https://graph"}},
            ),
            None,
        )
        mock_get_v.return_value = mock_v

        mw = AzureADAuthMiddleware(lambda r: HttpResponse("ok"))
        request = RequestFactory().get(
            "/api/things",
            HTTP_AUTHORIZATION="Bearer overage-token",
        )
        response = mw(request)

    assert response.status_code == 403


@pytest.mark.django_db
def test_no_groups_required_user_is_authenticated(monkeypatch, captured_request):
    """When ALLOWED_AZURE_AD_GROUPS is unset the group check is
    skipped — any valid-JWT user is accepted."""
    monkeypatch.setenv("CCP4I2_REQUIRE_AUTH", "true")
    monkeypatch.setenv("AZURE_AD_TENANT_ID", "tenant-id")
    monkeypatch.setenv("AZURE_AD_CLIENT_ID", "client-id")
    monkeypatch.delenv("ALLOWED_AZURE_AD_GROUPS", raising=False)

    handler, captured = captured_request

    with patch.object(azure_ad_module, "get_validator") as mock_get_v:
        mock_v = MagicMock()
        mock_v.validate_token.return_value = (True, _claims(), None)
        mock_get_v.return_value = mock_v

        mw = AzureADAuthMiddleware(handler)
        request = RequestFactory().get(
            "/api/things",
            HTTP_AUTHORIZATION="Bearer good-token",
        )
        response = mw(request)

    assert response.status_code == 200
    assert captured["flag"] is True


# ---------------------------------------------------------------------------
# User get-or-create from claims
# ---------------------------------------------------------------------------


@pytest.mark.django_db
def test_user_created_with_email_and_name_from_claims(
    monkeypatch, captured_request
):
    monkeypatch.setenv("CCP4I2_REQUIRE_AUTH", "true")
    monkeypatch.setenv("AZURE_AD_TENANT_ID", "tenant-id")
    monkeypatch.setenv("AZURE_AD_CLIENT_ID", "client-id")

    handler, captured = captured_request

    with patch.object(azure_ad_module, "get_validator") as mock_get_v:
        mock_v = MagicMock()
        mock_v.validate_token.return_value = (
            True,
            _claims(
                sub="abcd1234efgh5678ijkl9012mnop3456",
                email="ALICE@example.com",  # uppercase to verify lowercasing
                given_name="Alice",
                family_name="Smith",
            ),
            None,
        )
        mock_get_v.return_value = mock_v

        mw = AzureADAuthMiddleware(handler)
        request = RequestFactory().get(
            "/api/things",
            HTTP_AUTHORIZATION="Bearer good-token",
        )
        mw(request)

    user = captured["user"]
    assert user.email == "alice@example.com"  # normalised
    assert user.username.startswith("aad_abcd1234")
    assert user.first_name == "Alice"
    assert user.last_name == "Smith"


@pytest.mark.django_db
def test_x_user_email_header_fallback_when_claim_missing(
    monkeypatch, captured_request
):
    """When the JWT has no email/preferred_username/upn/unique_name
    claim, the middleware accepts X-User-Email as a display fallback.
    Safe because the JWT is already validated and 'sub' is the
    cryptographic primary key."""
    monkeypatch.setenv("CCP4I2_REQUIRE_AUTH", "true")
    monkeypatch.setenv("AZURE_AD_TENANT_ID", "tenant-id")
    monkeypatch.setenv("AZURE_AD_CLIENT_ID", "client-id")

    handler, captured = captured_request

    claims = {"sub": "user-sub-no-email"}  # no email-bearing claims

    with patch.object(azure_ad_module, "get_validator") as mock_get_v:
        mock_v = MagicMock()
        mock_v.validate_token.return_value = (True, claims, None)
        mock_get_v.return_value = mock_v

        mw = AzureADAuthMiddleware(handler)
        request = RequestFactory().get(
            "/api/things",
            HTTP_AUTHORIZATION="Bearer good-token",
            HTTP_X_USER_EMAIL="header-supplied@example.com",
        )
        mw(request)

    assert captured["user"].email == "header-supplied@example.com"


@pytest.mark.django_db
def test_user_email_falls_back_to_placeholder_when_no_source_available(
    monkeypatch, captured_request
):
    """Last-resort: when the JWT has no email-bearing claims AND no
    X-User-Email header is provided, the middleware synthesises an
    azuread.local placeholder from the verified ``sub``."""
    monkeypatch.setenv("CCP4I2_REQUIRE_AUTH", "true")
    monkeypatch.setenv("AZURE_AD_TENANT_ID", "tenant-id")
    monkeypatch.setenv("AZURE_AD_CLIENT_ID", "client-id")

    handler, captured = captured_request

    claims = {"sub": "user-sub-no-email"}

    with patch.object(azure_ad_module, "get_validator") as mock_get_v:
        mock_v = MagicMock()
        mock_v.validate_token.return_value = (True, claims, None)
        mock_get_v.return_value = mock_v

        mw = AzureADAuthMiddleware(handler)
        request = RequestFactory().get(
            "/api/things",
            HTTP_AUTHORIZATION="Bearer good-token",
        )
        mw(request)

    assert captured["user"].email == "user_user-sub@azuread.local"
