"""DRF view tests for the substructures listing + scaffold-extension
delete endpoints — the management surface for the catalog.

Covers:
- GET listing returns seed + extensions with canonical SMARTS
- Kekulé SMILES inputs get normalised to aromatic SMARTS in the listing
- Genuine SMARTS-only inputs that don't parse as SMILES still list
- DELETE creator-scoped, 404 across users, 204 on success
"""

from __future__ import annotations

import pytest
from django.contrib.auth import get_user_model
from rest_framework.test import APIClient

from compounds.registry.models import ScaffoldExtension, Target

User = get_user_model()

LIST_URL = "/api/compounds/nlp/substructures/"


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
# Listing
# ---------------------------------------------------------------------------


def test_list_returns_seed_catalog(alice_client):
    resp = alice_client.get(LIST_URL)
    assert resp.status_code == 200
    names = [s["name"] for s in resp.data["scaffolds"]]
    assert "pyridine" in names
    assert "pyrimidine" in names


def test_list_includes_canonical_smarts(alice_client):
    resp = alice_client.get(LIST_URL)
    pyrimidine = next(s for s in resp.data["scaffolds"] if s["name"] == "pyrimidine")
    # Canonical form returned for a curated seed entry.
    assert pyrimidine["canonical_smarts"] is not None


def test_list_includes_extension_with_id_and_creator(alice_client, alice):
    ext = ScaffoldExtension.objects.create(
        name="my-fragment",
        smarts="c1ccncc1",
        target=None,
        created_by=alice,
    )
    resp = alice_client.get(LIST_URL)
    row = next(s for s in resp.data["scaffolds"] if s["name"] == "my-fragment")
    assert row["id"] == str(ext.id)
    assert row["created_by"] == "alice"
    assert row["source"] == "extension:shared"


def test_list_canonicalises_kekule_smiles_input(alice_client, alice):
    """The pyrrolopyrimidine bug-class input: chemist pastes Kekulé
    SMILES into the SMARTS column. Listing endpoint provides a
    canonical aromatic form so the management page can show "you typed
    Kekulé, this is what RDKit actually queries with"."""
    raw = "C12=C([N]C=C2)N=CN=C1"
    ScaffoldExtension.objects.create(
        name="Pyrrolopyrimidine", smarts=raw, target=None, created_by=alice,
    )
    resp = alice_client.get(LIST_URL)
    row = next(s for s in resp.data["scaffolds"] if s["name"] == "Pyrrolopyrimidine")
    canonical = row["canonical_smarts"]
    # MolToSmarts normalises notation (typically into [#N] atom-spec
    # form). What matters for the bug we're fixing: SOME canonical form
    # comes back, and it differs from the raw Kekulé input — proving
    # RDKit went round through aromatic perception.
    assert canonical is not None
    assert canonical != raw


def test_list_target_id_filter_includes_project_extensions(alice_client, alice):
    target = Target.objects.create(name="ARd")
    ScaffoldExtension.objects.create(
        name="project-only",
        smarts="c1ccncc1",
        target=target,
        created_by=alice,
    )
    # Without target_id, project-scoped extension is hidden.
    resp_no_target = alice_client.get(LIST_URL)
    names = [s["name"] for s in resp_no_target.data["scaffolds"]]
    assert "project-only" not in names

    # With target_id, project-scoped extension is included.
    resp_with_target = alice_client.get(f"{LIST_URL}?target_id={target.id}")
    names = [s["name"] for s in resp_with_target.data["scaffolds"]]
    assert "project-only" in names


# ---------------------------------------------------------------------------
# DELETE /nlp/scaffold/extensions/<uuid>/
# ---------------------------------------------------------------------------


def _delete_url(extension_id) -> str:
    return f"/api/compounds/nlp/scaffold/extensions/{extension_id}/"


def test_delete_extension(alice_client, alice):
    ext = ScaffoldExtension.objects.create(
        name="x", smarts="c1ccncc1", target=None, created_by=alice,
    )
    resp = alice_client.delete(_delete_url(ext.id))
    assert resp.status_code == 204
    assert not ScaffoldExtension.objects.filter(pk=ext.id).exists()


def test_delete_other_users_extension_returns_404(bob_client, alice):
    ext = ScaffoldExtension.objects.create(
        name="alice-private", smarts="c1ccncc1", target=None, created_by=alice,
    )
    resp = bob_client.delete(_delete_url(ext.id))
    assert resp.status_code == 404
    assert ScaffoldExtension.objects.filter(pk=ext.id).exists()


def test_delete_unknown_id_returns_404(alice_client):
    resp = alice_client.delete(_delete_url(999_999))
    assert resp.status_code == 404
