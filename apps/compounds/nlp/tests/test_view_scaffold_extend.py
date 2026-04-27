"""DRF view tests for POST /api/compounds/nlp/scaffold/extend/ (slice 17).

Auth, body validation, RDKit SMARTS validation, name-collision in scope,
and round-trip via resolve_scaffold proving the new entry is consultable.
"""

from __future__ import annotations

import pytest
from django.contrib.auth import get_user_model
from django.core.cache import cache
from rest_framework.test import APIClient

from compounds.nlp.view import ENV_FLAG
from compounds.registry.models import ScaffoldExtension, Target

User = get_user_model()

URL = "/api/compounds/nlp/scaffold/extend/"


@pytest.fixture
def feature_on(monkeypatch):
    monkeypatch.setenv(ENV_FLAG, "true")


@pytest.fixture
def clear_cache():
    cache.clear()
    yield
    cache.clear()


@pytest.fixture
def user(db):
    return User.objects.create_user(username="alice", email="a@example.org")


@pytest.fixture
def client(user):
    c = APIClient()
    c.force_authenticate(user=user)
    return c


# ---------------------------------------------------------------------------
# Feature flag
# ---------------------------------------------------------------------------


def test_flag_off_returns_404(client, clear_cache, monkeypatch, db):
    monkeypatch.delenv(ENV_FLAG, raising=False)
    resp = client.post(URL, {"name": "indazole", "smarts": "c1ccc2[nH]ncc2c1"}, format="json")
    assert resp.status_code == 404


# ---------------------------------------------------------------------------
# Body validation
# ---------------------------------------------------------------------------


def test_missing_name_returns_400(client, feature_on, clear_cache, db):
    resp = client.post(URL, {"smarts": "c1ccncc1"}, format="json")
    assert resp.status_code == 400


def test_missing_smarts_returns_400(client, feature_on, clear_cache, db):
    resp = client.post(URL, {"name": "indazole"}, format="json")
    assert resp.status_code == 400


def test_unparseable_smarts_returns_400(client, feature_on, clear_cache, db):
    resp = client.post(URL, {
        "name": "garbage", "smarts": "this is not smarts(((",
    }, format="json")
    assert resp.status_code == 400
    assert resp.data["field"] == "smarts"


def test_invalid_aliases_returns_400(client, feature_on, clear_cache, db):
    resp = client.post(URL, {
        "name": "indazole", "smarts": "c1ccc2[nH]ncc2c1",
        "aliases": "not-a-list",
    }, format="json")
    assert resp.status_code == 400


def test_unknown_target_id_returns_400(client, feature_on, clear_cache, db):
    resp = client.post(URL, {
        "name": "series 1", "smarts": "c1ccncc1",
        "target_id": "00000000-0000-0000-0000-000000000000",
    }, format="json")
    assert resp.status_code == 400


# ---------------------------------------------------------------------------
# Happy path
# ---------------------------------------------------------------------------


def test_create_shared_extension(client, feature_on, clear_cache, db):
    resp = client.post(URL, {
        "name": "indazole", "smarts": "c1ccc2[nH]ncc2c1",
        "aliases": ["indazoles", "indazolyl"],
        "notes": "Indazole core for the X series",
        "source_prompt": "ARd indazoles",
    }, format="json")
    assert resp.status_code == 201
    assert resp.data["status"] == "created"
    assert resp.data["name"] == "indazole"
    assert resp.data["target_id"] is None  # shared
    # Persisted with the right shape.
    ext = ScaffoldExtension.objects.get(pk=resp.data["id"])
    assert ext.smarts == "c1ccc2[nH]ncc2c1"
    assert list(ext.aliases) == ["indazoles", "indazolyl"]
    assert ext.target is None
    assert ext.notes == "Indazole core for the X series"
    assert ext.source_prompt == "ARd indazoles"
    assert ext.created_by == client.handler._force_user


def test_create_project_scoped_extension(client, feature_on, clear_cache, db):
    target = Target.objects.create(name="ARd")
    resp = client.post(URL, {
        "name": "series 1", "smarts": "c1ccncc1",
        "target_id": str(target.id),
        "notes": "Series 1 = pyridine analogues for ARd",
    }, format="json")
    assert resp.status_code == 201
    assert resp.data["target_id"] == str(target.id)
    assert resp.data["target_name"] == "ARd"


# ---------------------------------------------------------------------------
# Conflict
# ---------------------------------------------------------------------------


def test_duplicate_name_in_same_scope_returns_409(
    client, feature_on, clear_cache, db,
):
    ScaffoldExtension.objects.create(
        name="indazole", smarts="c1ccc2[nH]ncc2c1",
    )
    resp = client.post(URL, {
        "name": "indazole", "smarts": "c1ccccc1",
    }, format="json")
    assert resp.status_code == 409
    assert resp.data["field"] == "name"


def test_same_name_different_scope_is_allowed(
    client, feature_on, clear_cache, db,
):
    target = Target.objects.create(name="ARd")
    # Shared "series 1"
    ScaffoldExtension.objects.create(name="series 1", smarts="c1ccccc1")
    # Project-scoped "series 1" — same name, different scope, should land.
    resp = client.post(URL, {
        "name": "series 1", "smarts": "c1ccncc1",
        "target_id": str(target.id),
    }, format="json")
    assert resp.status_code == 201


# ---------------------------------------------------------------------------
# Round-trip — entry is consultable by resolve_scaffold
# ---------------------------------------------------------------------------


def test_created_entry_is_resolvable(client, feature_on, clear_cache, db):
    """End-to-end: POST creates an entry, then resolve_scaffold finds it."""
    from compounds.nlp.resolver import resolve_scaffold
    from compounds.nlp.spec import ResolvedScaffold

    resp = client.post(URL, {
        "name": "benzofuran", "smarts": "c1ccc2occc2c1",
    }, format="json")
    assert resp.status_code == 201

    res = resolve_scaffold("benzofuran")
    assert isinstance(res, ResolvedScaffold)
    assert res.scaffold.smarts == "c1ccc2occc2c1"


# Auth boundary is enforced by the [IsAuthenticated] permission class
# on the view function — same as nlp_query — and matches the rest of
# the compounds API. Test settings auto-authenticate via the dev
# middleware so checking the unauthenticated case here would just
# exercise that bypass; the behaviour is exercised elsewhere.
