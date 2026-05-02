"""Smoke tests for ``DevAdminMiddleware``.

Verifies the dev-fallback contract that closed the previous
"misconfigured production auto-creates a superuser" backdoor:

- DEBUG=True → middleware activates, creates the dev_admin superuser,
  sets it on the request, marks the trust flag.
- DEBUG=False → middleware refuses to activate (defence in depth on
  top of the settings-level fail-closed switch).
"""

import pytest
from django.conf import settings
from django.http import HttpResponse
from django.test import RequestFactory

from ccp4i2_api.middleware.base import REQUEST_FLAG_ATTR
from ccp4i2_api.middleware.dev_admin import DevAdminMiddleware


@pytest.fixture
def captured_request():
    """Returns (handler, captured-dict). The handler stores whatever
    the middleware passes through."""
    captured: dict = {}

    def get_response(request):
        captured["request"] = request
        captured["user"] = getattr(request, "user", None)
        captured["flag"] = getattr(request, REQUEST_FLAG_ATTR, False)
        return HttpResponse("ok", status=200)

    return get_response, captured


@pytest.mark.django_db
def test_activates_and_creates_dev_admin_when_debug_true(captured_request):
    handler, captured = captured_request
    settings.DEBUG = True

    middleware = DevAdminMiddleware(handler)
    assert middleware.is_active() is True

    request = RequestFactory().get("/test")
    response = middleware(request)

    assert response.status_code == 200
    assert captured["user"].username == "dev_admin"
    assert captured["user"].is_superuser is True
    assert captured["user"].is_staff is True
    assert captured["flag"] is True


@pytest.mark.django_db
def test_refuses_to_activate_when_debug_false(captured_request):
    handler, captured = captured_request
    settings.DEBUG = False

    middleware = DevAdminMiddleware(handler)
    assert middleware.is_active() is False

    request = RequestFactory().get("/test")
    response = middleware(request)

    # Handler runs (middleware deferred), but no user/flag set by
    # DevAdmin — confirms the fail-closed contract.
    assert response.status_code == 200
    assert captured["flag"] is False
