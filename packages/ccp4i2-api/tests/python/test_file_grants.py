"""Tests for scoped file grants.

A grant is a capability: it authorises read-only access to one directory
subtree of the path-based file-serving endpoint, for a bounded time, on
behalf of the user who minted it. These tests pin the properties that make
that safe — scope containment (including against ``..`` traversal), method
restriction, expiry, signature binding — plus the middleware integration
that lets a grant stand in for a bearer token on browser-issued reads.
"""

import pytest
from django.http import HttpResponse
from django.test import RequestFactory, override_settings

from ccp4i2_api.file_grants import (
    GRANT_COOKIE,
    extract_grant,
    grant_user_pk,
    mint_grant,
    normalise_prefix,
)
from ccp4i2_api.middleware.base import REQUEST_FLAG_ATTR
from ccp4i2_api.middleware.local_session import LocalSessionAuthMiddleware

PREFIX = "/api/ccp4i2/projects/12/files_by_path/CCP4_JOBS/job_31/0000000001/"
PAGE = PREFIX + "test-page.html"


def _grant(user_pk="7", prefix=PREFIX):
    return mint_grant(user_pk=user_pk, path_prefix=prefix)


# ---------------------------------------------------------------------------
# normalise_prefix()
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "raw,expected",
    [
        ("/a/b", "/a/b/"),
        ("/a/b/", "/a/b/"),
        ("a/b", "/a/b/"),
        ("/a/b/../c", "/a/c/"),
        ("/a/./b", "/a/b/"),
    ],
)
def test_normalise_prefix_canonicalises(raw, expected):
    assert normalise_prefix(raw) == expected


# ---------------------------------------------------------------------------
# Scope
# ---------------------------------------------------------------------------


def test_grant_authorises_a_file_in_its_directory():
    assert grant_user_pk(_grant(), path=PAGE, method="GET") == "7"


def test_grant_authorises_a_deeper_subdirectory():
    assert (
        grant_user_pk(_grant(), path=PREFIX + "css/w3.css", method="GET") == "7"
    )


def test_grant_authorises_the_directory_itself():
    assert grant_user_pk(_grant(), path=PREFIX, method="GET") == "7"


def test_grant_refuses_a_sibling_directory():
    other = "/api/ccp4i2/projects/12/files_by_path/CCP4_JOBS/job_30/secrets.pdb"
    assert grant_user_pk(_grant(), path=other, method="GET") is None


def test_grant_refuses_another_project():
    other = "/api/ccp4i2/projects/99/files_by_path/CCP4_JOBS/job_31/0000000001/x"
    assert grant_user_pk(_grant(), path=other, method="GET") is None


def test_grant_refuses_dot_dot_traversal_out_of_its_subtree():
    escape = PREFIX + "../../job_30/secrets.pdb"
    assert grant_user_pk(_grant(), path=escape, method="GET") is None


def test_grant_refuses_a_prefix_match_that_is_not_a_path_boundary():
    # /…/0000000001-other/ must not be admitted by a /…/0000000001/ grant.
    sneaky = PREFIX.rstrip("/") + "-other/x.json"
    assert grant_user_pk(_grant(), path=sneaky, method="GET") is None


def test_mint_refuses_a_prefix_outside_the_file_serving_endpoint():
    with pytest.raises(ValueError):
        mint_grant(user_pk="7", path_prefix="/api/ccp4i2/projects/12/")


# ---------------------------------------------------------------------------
# Method, expiry, signature
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("method", ["GET", "HEAD"])
def test_grant_allows_safe_methods(method):
    assert grant_user_pk(_grant(), path=PAGE, method=method) == "7"


@pytest.mark.parametrize("method", ["POST", "PUT", "PATCH", "DELETE"])
def test_grant_refuses_unsafe_methods(method):
    assert grant_user_pk(_grant(), path=PAGE, method=method) is None


def test_grant_expires():
    token = _grant()
    with override_settings(CCP4I2_FILE_GRANT_TTL=-1):
        assert grant_user_pk(token, path=PAGE, method="GET") is None


def test_tampered_grant_is_rejected():
    token = _grant()
    assert grant_user_pk(token[:-1] + "x", path=PAGE, method="GET") is None


def test_garbage_is_rejected_rather_than_raising():
    assert grant_user_pk("not-a-token", path=PAGE, method="GET") is None


def test_grant_is_bound_to_the_launch_secret(monkeypatch):
    monkeypatch.setenv("CCP4I2_LOCAL_SESSION_TOKEN", "launch-one")
    token = _grant()
    assert grant_user_pk(token, path=PAGE, method="GET") == "7"
    # A new launch mints a new secret; grants from the old one die with it.
    monkeypatch.setenv("CCP4I2_LOCAL_SESSION_TOKEN", "launch-two")
    assert grant_user_pk(token, path=PAGE, method="GET") is None


# ---------------------------------------------------------------------------
# extract_grant()
# ---------------------------------------------------------------------------


def test_extract_grant_prefers_header_then_query_then_cookie():
    factory = RequestFactory()
    request = factory.get(
        PAGE, {"file_grant": "from-query"}, HTTP_X_CCP4I2_FILE_GRANT="from-header"
    )
    request.COOKIES[GRANT_COOKIE] = "from-cookie"
    assert extract_grant(request) == "from-header"

    request = factory.get(PAGE, {"file_grant": "from-query"})
    request.COOKIES[GRANT_COOKIE] = "from-cookie"
    assert extract_grant(request) == "from-query"

    request = factory.get(PAGE)
    request.COOKIES[GRANT_COOKIE] = "from-cookie"
    assert extract_grant(request) == "from-cookie"


def test_extract_grant_returns_none_when_absent():
    assert extract_grant(RequestFactory().get(PAGE)) is None


# ---------------------------------------------------------------------------
# Middleware integration
# ---------------------------------------------------------------------------


@pytest.fixture
def captured():
    state: dict = {}

    def get_response(request):
        state["user"] = getattr(request, "user", None)
        state["flag"] = getattr(request, REQUEST_FLAG_ATTR, False)
        return HttpResponse("ok", status=200)

    return get_response, state


@pytest.mark.django_db
def test_middleware_accepts_a_grant_cookie_without_a_bearer_token(
    monkeypatch, captured
):
    from django.contrib.auth import get_user_model

    monkeypatch.setenv("CCP4I2_LOCAL_SESSION_TOKEN", "the-launch-secret")
    user = get_user_model().objects.create(username="reader", email="r@x.invalid")

    get_response, state = captured
    middleware = LocalSessionAuthMiddleware(get_response)
    request = RequestFactory().get(PAGE)
    request.COOKIES[GRANT_COOKIE] = mint_grant(
        user_pk=user.pk, path_prefix=PREFIX
    )

    response = middleware(request)

    assert response.status_code == 200
    assert state["user"] == user
    assert state["flag"] is True


@pytest.mark.django_db
def test_middleware_rejects_a_grant_for_a_path_outside_its_scope(
    monkeypatch, captured
):
    from django.contrib.auth import get_user_model

    monkeypatch.setenv("CCP4I2_LOCAL_SESSION_TOKEN", "the-launch-secret")
    user = get_user_model().objects.create(username="reader2", email="r2@x.invalid")

    get_response, _ = captured
    middleware = LocalSessionAuthMiddleware(get_response)
    request = RequestFactory().get(
        "/api/ccp4i2/projects/12/files_by_path/CCP4_JOBS/job_30/secrets.pdb"
    )
    request.COOKIES[GRANT_COOKIE] = mint_grant(
        user_pk=user.pk, path_prefix=PREFIX
    )

    response = middleware(request)

    # Falls through to bearer-token auth, which this request has none of.
    assert response.status_code == 401


@pytest.mark.django_db
def test_a_stale_grant_does_not_shadow_a_valid_bearer_token(monkeypatch, captured):
    monkeypatch.setenv("CCP4I2_LOCAL_SESSION_TOKEN", "the-launch-secret")

    get_response, state = captured
    middleware = LocalSessionAuthMiddleware(get_response)
    request = RequestFactory().get(
        PAGE, HTTP_AUTHORIZATION="Bearer the-launch-secret"
    )
    request.COOKIES[GRANT_COOKIE] = "expired-or-garbage"

    response = middleware(request)

    assert response.status_code == 200
    assert state["flag"] is True


@pytest.mark.django_db
def test_credentials_win_over_a_perfectly_good_grant(monkeypatch, captured):
    """A request that carries credentials is the identity they name.

    Following someone else's grant URL while signed in must not silently
    turn you into them for that subtree.
    """
    from django.contrib.auth import get_user_model

    monkeypatch.setenv("CCP4I2_LOCAL_SESSION_TOKEN", "the-launch-secret")
    minter = get_user_model().objects.create(username="minter", email="m@x.invalid")

    get_response, state = captured
    middleware = LocalSessionAuthMiddleware(get_response)
    request = RequestFactory().get(
        PAGE, HTTP_AUTHORIZATION="Bearer the-launch-secret"
    )
    request.COOKIES[GRANT_COOKIE] = mint_grant(
        user_pk=minter.pk, path_prefix=PREFIX
    )

    response = middleware(request)

    assert response.status_code == 200
    # The desktop user the bearer token authenticates, not the grant's minter.
    assert state["user"] != minter
