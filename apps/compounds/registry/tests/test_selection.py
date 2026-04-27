"""Selection model + endpoint tests (slice 20).

Covers the URL-length-fix infrastructure: token-addressed snapshot of a
compound list, with auto-expiry, creator-scoped read access, and the
cleanup management command.
"""

from __future__ import annotations

import datetime

import pytest
from django.contrib.auth import get_user_model
from django.core.management import call_command
from django.utils import timezone
from io import StringIO
from rest_framework.test import APIClient

from compounds.registry.models import Selection

User = get_user_model()

CREATE_URL = "/api/compounds/selections/"


@pytest.fixture
def alice(db):
    return User.objects.create_user(username="alice", email="alice@example.org")


@pytest.fixture
def bob(db):
    return User.objects.create_user(username="bob", email="bob@example.org")


@pytest.fixture
def alice_client(alice):
    c = APIClient()
    c.force_authenticate(user=alice)
    return c


@pytest.fixture
def bob_client(bob):
    c = APIClient()
    c.force_authenticate(user=bob)
    return c


# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------


def test_selection_str_includes_compound_count(alice):
    sel = Selection.objects.create(
        name="ARd HTRF top 20",
        compound_ids=["NCL-1", "NCL-2", "NCL-3"],
        created_by=alice,
    )
    assert "3 compounds" in str(sel)


def test_is_expired_false_when_expires_at_in_future(alice):
    sel = Selection.objects.create(
        name="x",
        compound_ids=["NCL-1"],
        created_by=alice,
        expires_at=timezone.now() + datetime.timedelta(hours=1),
    )
    assert sel.is_expired is False


def test_is_expired_true_when_expires_at_in_past(alice):
    sel = Selection.objects.create(
        name="x",
        compound_ids=["NCL-1"],
        created_by=alice,
        expires_at=timezone.now() - datetime.timedelta(hours=1),
    )
    assert sel.is_expired is True


def test_is_expired_false_when_expires_at_null(alice):
    """Saved selections (slice 22) have null expires_at and never expire."""
    sel = Selection.objects.create(
        name="x", compound_ids=["NCL-1"], created_by=alice,
    )
    assert sel.is_expired is False


# ---------------------------------------------------------------------------
# POST /selections/
# ---------------------------------------------------------------------------


def test_create_selection_returns_token(alice_client, alice):
    resp = alice_client.post(CREATE_URL, {
        "name": "ARd HTRF top 20",
        "compound_ids": ["NCL-00000001", "NCL-00000002", "NCL-00000003"],
        "source_prompt": "best 20 ARd compounds",
    }, format="json")
    assert resp.status_code == 201
    assert "id" in resp.data
    assert resp.data["n_compounds"] == 3
    assert resp.data["expires_at"] is not None  # 7-day expiry default

    sel = Selection.objects.get(pk=resp.data["id"])
    assert sel.created_by == alice
    assert sel.compound_ids == ["NCL-00000001", "NCL-00000002", "NCL-00000003"]
    assert sel.source_prompt == "best 20 ARd compounds"


def test_create_selection_rejects_non_list_compound_ids(alice_client):
    resp = alice_client.post(CREATE_URL, {
        "name": "x",
        "compound_ids": "NCL-1,NCL-2",  # string, not list
    }, format="json")
    assert resp.status_code == 400


def test_create_selection_rejects_list_with_non_string_entries(alice_client):
    resp = alice_client.post(CREATE_URL, {
        "name": "x",
        "compound_ids": ["NCL-1", 42, "NCL-3"],
    }, format="json")
    assert resp.status_code == 400


def test_create_selection_accepts_empty_list(alice_client, alice):
    """Edge case — a selection with zero compounds still persists. Lets
    the aggregation page render an empty result rather than 404."""
    resp = alice_client.post(CREATE_URL, {
        "name": "Found 0", "compound_ids": [],
    }, format="json")
    assert resp.status_code == 201
    assert resp.data["n_compounds"] == 0


# ---------------------------------------------------------------------------
# GET /selections/<uuid>/
# ---------------------------------------------------------------------------


def test_get_selection_returns_compound_ids(alice_client, alice):
    sel = Selection.objects.create(
        name="t", compound_ids=["NCL-1", "NCL-2"],
        created_by=alice,
        expires_at=timezone.now() + datetime.timedelta(days=1),
    )
    resp = alice_client.get(f"/api/compounds/selections/{sel.id}/")
    assert resp.status_code == 200
    assert resp.data["compound_ids"] == ["NCL-1", "NCL-2"]
    assert resp.data["name"] == "t"


def test_get_selection_404_for_other_users_token(bob_client, alice):
    """Don't leak existence of selections across users — 404, not 403."""
    sel = Selection.objects.create(
        name="alice-private", compound_ids=["NCL-1"], created_by=alice,
        expires_at=timezone.now() + datetime.timedelta(days=1),
    )
    resp = bob_client.get(f"/api/compounds/selections/{sel.id}/")
    assert resp.status_code == 404


def test_get_selection_410_when_expired(alice_client, alice):
    sel = Selection.objects.create(
        name="t", compound_ids=["NCL-1"], created_by=alice,
        expires_at=timezone.now() - datetime.timedelta(hours=1),
    )
    resp = alice_client.get(f"/api/compounds/selections/{sel.id}/")
    assert resp.status_code == 410


def test_get_selection_404_for_unknown_uuid(alice_client):
    resp = alice_client.get(
        "/api/compounds/selections/00000000-0000-0000-0000-000000000000/"
    )
    assert resp.status_code == 404


# ---------------------------------------------------------------------------
# Cleanup command
# ---------------------------------------------------------------------------


def test_cleanup_command_deletes_only_expired(alice):
    fresh = Selection.objects.create(
        name="fresh", compound_ids=["NCL-1"], created_by=alice,
        expires_at=timezone.now() + datetime.timedelta(days=1),
    )
    expired = Selection.objects.create(
        name="expired", compound_ids=["NCL-2"], created_by=alice,
        expires_at=timezone.now() - datetime.timedelta(hours=1),
    )
    saved = Selection.objects.create(
        name="saved", compound_ids=["NCL-3"], created_by=alice,
        expires_at=None,    # saved selections don't expire
    )

    out = StringIO()
    call_command("cleanup_selections", stdout=out)

    assert Selection.objects.filter(pk=fresh.pk).exists()
    assert not Selection.objects.filter(pk=expired.pk).exists()
    assert Selection.objects.filter(pk=saved.pk).exists()
    assert "Deleted 1" in out.getvalue()


def test_cleanup_dry_run_does_not_delete(alice):
    expired = Selection.objects.create(
        name="expired", compound_ids=["NCL-2"], created_by=alice,
        expires_at=timezone.now() - datetime.timedelta(hours=1),
    )
    out = StringIO()
    call_command("cleanup_selections", "--dry-run", stdout=out)
    assert Selection.objects.filter(pk=expired.pk).exists()
    assert "DRY RUN" in out.getvalue()
    assert "would delete 1" in out.getvalue()


def test_cleanup_command_idempotent_when_nothing_expired(alice):
    Selection.objects.create(
        name="fresh", compound_ids=["NCL-1"], created_by=alice,
        expires_at=timezone.now() + datetime.timedelta(days=1),
    )
    out = StringIO()
    call_command("cleanup_selections", stdout=out)
    assert "No expired selections" in out.getvalue()
