"""Poly-source registrant resolver tests (slice 16).

The resolver matches a typed phrase against BOTH the User pool (chemists
who clicked Register) and the Supplier pool (vendors / personal-synthesis
sources). De-dups user-linked Suppliers against their Users.
"""

from __future__ import annotations

import pytest
from django.contrib.auth import get_user_model

from compounds.nlp.resolver import resolve_registrant
from compounds.nlp.spec import (
    FIELD_REGISTERED_BY,
    ResolvedSupplier,
    ResolvedUser,
    SupplierCandidate,
    UserCandidate,
    UserClarify,
    UserMiss,
)
from compounds.registry.models import Compound, Supplier, Target

User = get_user_model()


@pytest.fixture
def registrant_world(db):
    """A mix of users and suppliers with realistic shapes:

    - Alice Jones — User who registers compounds
    - Bob Smith — User who registers compounds; ALSO a personal Supplier
      linked to Bob (de-dup case)
    - Enamine — pure vendor Supplier with compounds, no User link
    - Sigma — pure vendor Supplier with no compounds (filtered out)
    - Charlie Smith — User who registers compounds (collides with Bob
      Smith on surname — keeps the Clarify path honest)
    """
    t = Target.objects.create(name="ARd")

    alice = User.objects.create_user(
        username="alice.jones", email="alice.jones@ncl.ac.uk",
        first_name="Alice", last_name="Jones",
    )
    bob = User.objects.create_user(
        username="bob.smith", email="bob.smith@ncl.ac.uk",
        first_name="Bob", last_name="Smith",
    )
    charlie = User.objects.create_user(
        username="charlie.smith", email="charlie.smith@ncl.ac.uk",
        first_name="Charlie", last_name="Smith",
    )
    bob_supplier = Supplier.objects.create(name="Bob Smith", user=bob)
    enamine = Supplier.objects.create(name="Enamine", initials="ENM")
    Supplier.objects.create(name="Sigma")  # no compounds → filtered out

    Compound.objects.create(target=t, smiles="CCO", registered_by=alice)
    Compound.objects.create(target=t, smiles="CCN", registered_by=bob)
    Compound.objects.create(target=t, smiles="CCC", registered_by=charlie)
    Compound.objects.create(target=t, smiles="CCO", supplier=enamine)
    # A compound supplied by Bob's personal-supplier — should de-dup
    # against Bob the User.
    Compound.objects.create(target=t, smiles="CCN", supplier=bob_supplier)

    return {
        "alice": alice, "bob": bob, "charlie": charlie,
        "bob_supplier": bob_supplier, "enamine": enamine,
    }


# ---------------------------------------------------------------------------
# User-only resolution (existing behaviour preserved)
# ---------------------------------------------------------------------------


def test_resolves_unique_user_full_name(registrant_world):
    res = resolve_registrant("Alice Jones", field=FIELD_REGISTERED_BY)
    assert isinstance(res, ResolvedUser)
    assert res.user.pk == registrant_world["alice"].pk


def test_resolves_unique_user_username(registrant_world):
    res = resolve_registrant("alice.jones", field=FIELD_REGISTERED_BY)
    assert isinstance(res, ResolvedUser)
    assert res.user.pk == registrant_world["alice"].pk


# ---------------------------------------------------------------------------
# Supplier resolution
# ---------------------------------------------------------------------------


def test_resolves_pure_vendor_supplier_by_name(registrant_world):
    res = resolve_registrant("Enamine", field=FIELD_REGISTERED_BY)
    assert isinstance(res, ResolvedSupplier)
    assert res.supplier.pk == registrant_world["enamine"].pk


def test_resolves_supplier_by_initials(registrant_world):
    res = resolve_registrant("ENM", field=FIELD_REGISTERED_BY)
    assert isinstance(res, ResolvedSupplier)
    assert res.supplier.pk == registrant_world["enamine"].pk


def test_supplier_with_no_compounds_is_filtered_out(registrant_world):
    """Sigma exists but has no compounds → not in the pool → Miss."""
    res = resolve_registrant("Sigma", field=FIELD_REGISTERED_BY)
    assert isinstance(res, UserMiss)


# ---------------------------------------------------------------------------
# De-duplication of user-linked suppliers
# ---------------------------------------------------------------------------


def test_user_linked_supplier_dedupes_against_user(registrant_world):
    """Bob Smith is both a User AND a personal-Supplier linked to that
    same User. Typing 'Bob Smith' must surface only the User candidate
    — chemist sees one chip, not two for the same person."""
    res = resolve_registrant("Bob Smith", field=FIELD_REGISTERED_BY)
    # With Charlie Smith also matching 'Smith' partially, this is a
    # Clarify across {Bob (user), Charlie (user)}. Bob's supplier
    # should NOT appear in the chips.
    if isinstance(res, ResolvedUser):
        assert res.user.pk == registrant_world["bob"].pk
    else:
        assert isinstance(res, UserClarify)
        kinds = [getattr(c, "kind", "user") for c in res.candidates]
        # No supplier chip for Bob (de-duped against Bob the user).
        supplier_ids = {
            c.id for c, k in zip(res.candidates, kinds) if k == "supplier"
        }
        assert str(registrant_world["bob_supplier"].pk) not in supplier_ids


# ---------------------------------------------------------------------------
# Mixed clarify across pools
# ---------------------------------------------------------------------------


def test_clarify_mixes_user_and_supplier_when_both_match(db):
    """Real-world collision: a User named 'Martin' and an unrelated
    Supplier named 'Martin' (e.g. a small vendor). Typing 'Martin'
    surfaces both as separate chips, with the supplier flagged via
    `kind='supplier'` so the frontend pins to the right field."""
    t = Target.objects.create(name="ARd")
    martin_user = User.objects.create_user(
        username="martin.noble", email="martin@ncl.ac.uk",
        first_name="Martin", last_name="Noble",
    )
    martin_supplier = Supplier.objects.create(name="Martin")
    Compound.objects.create(target=t, smiles="CCO", registered_by=martin_user)
    Compound.objects.create(target=t, smiles="CCC", supplier=martin_supplier)

    res = resolve_registrant("Martin", field=FIELD_REGISTERED_BY)
    assert isinstance(res, UserClarify)
    kinds = {getattr(c, "kind", "user") for c in res.candidates}
    assert "user" in kinds
    assert "supplier" in kinds


# ---------------------------------------------------------------------------
# Pinning continuation
# ---------------------------------------------------------------------------


def test_pinned_user_id_resolves_to_user(registrant_world):
    res = resolve_registrant(
        "Martin",   # arbitrary — pin overrides the typed phrase
        field=FIELD_REGISTERED_BY,
        pinned_user_id=str(registrant_world["alice"].pk),
    )
    assert isinstance(res, ResolvedUser)
    assert res.user.pk == registrant_world["alice"].pk


def test_pinned_supplier_id_resolves_to_supplier(registrant_world):
    res = resolve_registrant(
        "Whoever",
        field=FIELD_REGISTERED_BY,
        pinned_supplier_id=str(registrant_world["enamine"].pk),
    )
    assert isinstance(res, ResolvedSupplier)
    assert res.supplier.pk == registrant_world["enamine"].pk


# ---------------------------------------------------------------------------
# Miss path
# ---------------------------------------------------------------------------


def test_unknown_typed_phrase_misses_with_suggestions(registrant_world):
    res = resolve_registrant("Zachary", field=FIELD_REGISTERED_BY)
    assert isinstance(res, UserMiss)
    assert len(res.suggestions) > 0


def test_miss_suggestions_include_both_users_and_suppliers(registrant_world):
    """The Miss suggestion list should include candidates from both
    pools so the chemist can scan recent options."""
    res = resolve_registrant("Zachary", field=FIELD_REGISTERED_BY)
    assert isinstance(res, UserMiss)
    kinds = {getattr(c, "kind", "user") for c in res.suggestions}
    assert "user" in kinds
    assert "supplier" in kinds
