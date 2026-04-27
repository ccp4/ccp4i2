"""Selection model + endpoint tests.

Covers the token-addressed compound-list snapshot: auto-expiry,
creator-scoped read access, the cleanup management command, and the
saved-selection / session-list CRUD surface.
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
    """Saved selections have null expires_at and never expire."""
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


def test_cleanup_command_skips_saved_even_with_past_expires_at(alice):
    """A saved selection with a stale expires_at must NOT be deleted —
    ``is_saved=True`` always wins over expiry. (Defensive; the toggle-off
    path normally re-applies expiry, but the cleanup must not depend on
    it.)"""
    saved_with_stale = Selection.objects.create(
        name="saved-but-stale", compound_ids=["NCL-1"], created_by=alice,
        is_saved=True,
        expires_at=timezone.now() - datetime.timedelta(days=10),
    )
    call_command("cleanup_selections", stdout=StringIO())
    assert Selection.objects.filter(pk=saved_with_stale.pk).exists()


# ---------------------------------------------------------------------------
# GET /selections/  — session list
# ---------------------------------------------------------------------------


LIST_URL = "/api/compounds/selections/"


def test_list_selections_returns_only_creators_rows(alice_client, alice, bob):
    Selection.objects.create(
        name="alice-1", compound_ids=["NCL-1"], created_by=alice,
        expires_at=timezone.now() + datetime.timedelta(days=1),
    )
    Selection.objects.create(
        name="bob-1", compound_ids=["NCL-2"], created_by=bob,
        expires_at=timezone.now() + datetime.timedelta(days=1),
    )
    resp = alice_client.get(LIST_URL)
    assert resp.status_code == 200
    names = [s["name"] for s in resp.data["selections"]]
    assert names == ["alice-1"]


def test_list_selections_ordered_most_recent_first(alice_client, alice):
    older = Selection.objects.create(
        name="older", compound_ids=["NCL-1"], created_by=alice,
        expires_at=timezone.now() + datetime.timedelta(days=1),
    )
    Selection.objects.filter(pk=older.pk).update(
        created_at=timezone.now() - datetime.timedelta(hours=1),
    )
    Selection.objects.create(
        name="newer", compound_ids=["NCL-2"], created_by=alice,
        expires_at=timezone.now() + datetime.timedelta(days=1),
    )
    resp = alice_client.get(LIST_URL)
    names = [s["name"] for s in resp.data["selections"]]
    assert names == ["newer", "older"]


def test_list_selections_omits_compound_ids(alice_client, alice):
    """The list view's payload shouldn't carry every selection's full
    compound list — that's the detail endpoint's job."""
    Selection.objects.create(
        name="x", compound_ids=["NCL-1", "NCL-2", "NCL-3"], created_by=alice,
        expires_at=timezone.now() + datetime.timedelta(days=1),
    )
    resp = alice_client.get(LIST_URL)
    body = resp.data["selections"][0]
    assert "compound_ids" not in body
    assert body["n_compounds"] == 3


def test_list_selections_hides_expired_unsaved_rows(alice_client, alice):
    Selection.objects.create(
        name="fresh", compound_ids=["NCL-1"], created_by=alice,
        expires_at=timezone.now() + datetime.timedelta(days=1),
    )
    Selection.objects.create(
        name="expired", compound_ids=["NCL-2"], created_by=alice,
        expires_at=timezone.now() - datetime.timedelta(hours=1),
    )
    resp = alice_client.get(LIST_URL)
    names = [s["name"] for s in resp.data["selections"]]
    assert names == ["fresh"]


def test_list_selections_includes_expired_saved_rows(alice_client, alice):
    """Saved selections never expire — even if a row carries a stale
    expires_at, ``is_saved=True`` keeps it in the list."""
    Selection.objects.create(
        name="saved-stale", compound_ids=["NCL-1"], created_by=alice,
        is_saved=True,
        expires_at=timezone.now() - datetime.timedelta(days=10),
    )
    resp = alice_client.get(LIST_URL)
    names = [s["name"] for s in resp.data["selections"]]
    assert names == ["saved-stale"]


def test_list_selections_filters_is_saved_true(alice_client, alice):
    Selection.objects.create(
        name="saved", compound_ids=["NCL-1"], created_by=alice,
        is_saved=True, expires_at=None,
    )
    Selection.objects.create(
        name="ephemeral", compound_ids=["NCL-2"], created_by=alice,
        is_saved=False,
        expires_at=timezone.now() + datetime.timedelta(days=1),
    )
    resp = alice_client.get(LIST_URL + "?is_saved=true")
    names = [s["name"] for s in resp.data["selections"]]
    assert names == ["saved"]


def test_list_selections_filters_is_saved_false(alice_client, alice):
    Selection.objects.create(
        name="saved", compound_ids=["NCL-1"], created_by=alice,
        is_saved=True, expires_at=None,
    )
    Selection.objects.create(
        name="ephemeral", compound_ids=["NCL-2"], created_by=alice,
        is_saved=False,
        expires_at=timezone.now() + datetime.timedelta(days=1),
    )
    resp = alice_client.get(LIST_URL + "?is_saved=false")
    names = [s["name"] for s in resp.data["selections"]]
    assert names == ["ephemeral"]


# ---------------------------------------------------------------------------
# PATCH /selections/<uuid>/  — save / rename / scope
# ---------------------------------------------------------------------------


def test_patch_save_clears_expires_at(alice_client, alice):
    sel = Selection.objects.create(
        name="x", compound_ids=["NCL-1"], created_by=alice,
        expires_at=timezone.now() + datetime.timedelta(days=1),
    )
    resp = alice_client.patch(
        f"/api/compounds/selections/{sel.id}/",
        {"is_saved": True}, format="json",
    )
    assert resp.status_code == 200
    sel.refresh_from_db()
    assert sel.is_saved is True
    assert sel.expires_at is None
    assert resp.data["is_saved"] is True
    assert resp.data["expires_at"] is None


def test_patch_unsave_reapplies_expires_at(alice_client, alice):
    sel = Selection.objects.create(
        name="x", compound_ids=["NCL-1"], created_by=alice,
        is_saved=True, expires_at=None,
    )
    resp = alice_client.patch(
        f"/api/compounds/selections/{sel.id}/",
        {"is_saved": False}, format="json",
    )
    assert resp.status_code == 200
    sel.refresh_from_db()
    assert sel.is_saved is False
    assert sel.expires_at is not None
    assert sel.expires_at > timezone.now()


def test_patch_rename(alice_client, alice):
    sel = Selection.objects.create(
        name="old name", compound_ids=["NCL-1"], created_by=alice,
        expires_at=timezone.now() + datetime.timedelta(days=1),
    )
    resp = alice_client.patch(
        f"/api/compounds/selections/{sel.id}/",
        {"name": "new name"}, format="json",
    )
    assert resp.status_code == 200
    sel.refresh_from_db()
    assert sel.name == "new name"


def test_patch_rejects_empty_name(alice_client, alice):
    sel = Selection.objects.create(
        name="x", compound_ids=["NCL-1"], created_by=alice,
        expires_at=timezone.now() + datetime.timedelta(days=1),
    )
    resp = alice_client.patch(
        f"/api/compounds/selections/{sel.id}/",
        {"name": "   "}, format="json",
    )
    assert resp.status_code == 400


def test_patch_rejects_non_bool_is_saved(alice_client, alice):
    sel = Selection.objects.create(
        name="x", compound_ids=["NCL-1"], created_by=alice,
        expires_at=timezone.now() + datetime.timedelta(days=1),
    )
    resp = alice_client.patch(
        f"/api/compounds/selections/{sel.id}/",
        {"is_saved": "yes"}, format="json",
    )
    assert resp.status_code == 400


def test_patch_set_target(alice_client, alice):
    from compounds.registry.models import Target
    target = Target.objects.create(name="ARd")
    sel = Selection.objects.create(
        name="x", compound_ids=["NCL-1"], created_by=alice,
        expires_at=timezone.now() + datetime.timedelta(days=1),
    )
    resp = alice_client.patch(
        f"/api/compounds/selections/{sel.id}/",
        {"target_id": str(target.id)}, format="json",
    )
    assert resp.status_code == 200
    sel.refresh_from_db()
    assert sel.target_id == target.id
    assert resp.data["target_name"] == "ARd"


def test_patch_clear_target(alice_client, alice):
    from compounds.registry.models import Target
    target = Target.objects.create(name="ARd")
    sel = Selection.objects.create(
        name="x", compound_ids=["NCL-1"], created_by=alice,
        target=target,
        expires_at=timezone.now() + datetime.timedelta(days=1),
    )
    resp = alice_client.patch(
        f"/api/compounds/selections/{sel.id}/",
        {"target_id": None}, format="json",
    )
    assert resp.status_code == 200
    sel.refresh_from_db()
    assert sel.target_id is None


def test_patch_unknown_target_id_400(alice_client, alice):
    sel = Selection.objects.create(
        name="x", compound_ids=["NCL-1"], created_by=alice,
        expires_at=timezone.now() + datetime.timedelta(days=1),
    )
    resp = alice_client.patch(
        f"/api/compounds/selections/{sel.id}/",
        {"target_id": "00000000-0000-0000-0000-000000000000"}, format="json",
    )
    assert resp.status_code == 400


def test_patch_other_users_selection_404(bob_client, alice):
    sel = Selection.objects.create(
        name="alice-private", compound_ids=["NCL-1"], created_by=alice,
        expires_at=timezone.now() + datetime.timedelta(days=1),
    )
    resp = bob_client.patch(
        f"/api/compounds/selections/{sel.id}/",
        {"is_saved": True}, format="json",
    )
    assert resp.status_code == 404


# ---------------------------------------------------------------------------
# DELETE /selections/<uuid>/
# ---------------------------------------------------------------------------


def test_delete_selection(alice_client, alice):
    sel = Selection.objects.create(
        name="x", compound_ids=["NCL-1"], created_by=alice,
        expires_at=timezone.now() + datetime.timedelta(days=1),
    )
    resp = alice_client.delete(f"/api/compounds/selections/{sel.id}/")
    assert resp.status_code == 204
    assert not Selection.objects.filter(pk=sel.pk).exists()


def test_delete_other_users_selection_404(bob_client, alice):
    sel = Selection.objects.create(
        name="alice-private", compound_ids=["NCL-1"], created_by=alice,
        expires_at=timezone.now() + datetime.timedelta(days=1),
    )
    resp = bob_client.delete(f"/api/compounds/selections/{sel.id}/")
    assert resp.status_code == 404
    assert Selection.objects.filter(pk=sel.pk).exists()


# ---------------------------------------------------------------------------
# POST /selections/  — saved-on-create + target
# ---------------------------------------------------------------------------


def test_create_selection_with_is_saved_true_has_no_expiry(alice_client):
    resp = alice_client.post(CREATE_URL, {
        "name": "saved on create",
        "compound_ids": ["NCL-1"],
        "is_saved": True,
    }, format="json")
    assert resp.status_code == 201
    assert resp.data["is_saved"] is True
    assert resp.data["expires_at"] is None


def test_create_selection_with_target(alice_client, alice):
    from compounds.registry.models import Target
    target = Target.objects.create(name="ARd")
    resp = alice_client.post(CREATE_URL, {
        "name": "x", "compound_ids": ["NCL-1"],
        "target_id": str(target.id),
    }, format="json")
    assert resp.status_code == 201
    assert resp.data["target_name"] == "ARd"


def test_create_selection_unknown_target_400(alice_client):
    resp = alice_client.post(CREATE_URL, {
        "name": "x", "compound_ids": ["NCL-1"],
        "target_id": "00000000-0000-0000-0000-000000000000",
    }, format="json")
    assert resp.status_code == 400
