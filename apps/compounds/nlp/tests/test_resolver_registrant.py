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


# ---------------------------------------------------------------------------
# Slice 18: same-name User+Supplier merger (Hannah Stewart case)
# ---------------------------------------------------------------------------


@pytest.fixture
def hannah_world(db):
    """Hannah Stewart is a User who has registered three compounds AND
    has a separate (unlinked) Supplier row "Hannah Stewart" supplying
    two compounds. Pre-slice-18 this surfaced as two chips for the
    same person; slice 18 merges them into one ``UnionCandidate``."""
    from compounds.registry.models import Supplier
    t = Target.objects.create(name="ARd")
    hannah_user = User.objects.create_user(
        username="hannah.stewart", email="hannah.stewart@ncl.ac.uk",
        first_name="Hannah", last_name="Stewart",
    )
    # Unlinked Supplier with the same name — the slice-18 case.
    hannah_sup = Supplier.objects.create(name="Hannah Stewart")
    # 3 registered + 2 supplied.
    for _ in range(3):
        Compound.objects.create(target=t, smiles="CCO", registered_by=hannah_user)
    for _ in range(2):
        Compound.objects.create(target=t, smiles="CCN", supplier=hannah_sup)
    return {"user": hannah_user, "supplier": hannah_sup}


def test_same_name_user_and_unlinked_supplier_merge_into_resolved_union(hannah_world):
    """Single-resolution case: only Hannah matches, and she matches in
    BOTH pools. Should resolve to ``ResolvedUnion``, not Clarify."""
    from compounds.nlp.spec import ResolvedUnion
    res = resolve_registrant("Hannah Stewart", field=FIELD_REGISTERED_BY)
    assert isinstance(res, ResolvedUnion)
    assert res.user.pk == hannah_world["user"].pk
    assert res.supplier.pk == hannah_world["supplier"].pk


def test_union_resolution_persists_compound_count(hannah_world):
    """When the resolved-union case appears in a Clarify (because of
    other ambiguities), the chip carries the registered + supplied
    counts so the chemist sees the merger explicitly."""
    from compounds.nlp.spec import (
        ResolvedUnion, UnionCandidate, UserClarify,
    )
    # Add a second "Hannah" to force a Clarify rather than direct resolution.
    User.objects.create_user(
        username="hannah.jones", email="hannah.jones@ncl.ac.uk",
        first_name="Hannah", last_name="Jones",
    )
    from compounds.registry.models import Compound
    Compound.objects.create(
        target=Target.objects.first(),
        smiles="CCO",
        registered_by=User.objects.get(username="hannah.jones"),
    )
    res = resolve_registrant("Hannah", field=FIELD_REGISTERED_BY)
    assert isinstance(res, UserClarify)
    union_chips = [
        c for c in res.candidates
        if isinstance(c, UnionCandidate)
    ]
    # Hannah Stewart is the only same-name pair → one union chip.
    assert len(union_chips) == 1
    chip = union_chips[0]
    assert chip.n_compounds_user == 3
    assert chip.n_compounds_supplier == 2


def test_unlinked_supplier_with_different_name_does_not_merge(db):
    """If the User name and the unlinked Supplier name don't match
    after normalisation, no merger — they're surfaced as two separate
    chips."""
    from compounds.registry.models import Supplier
    from compounds.nlp.spec import (
        ResolvedUnion, UnionCandidate, UserClarify,
    )
    t = Target.objects.create(name="ARd")
    user = User.objects.create_user(
        username="hannah.stewart", email="hannah.stewart@ncl.ac.uk",
        first_name="Hannah", last_name="Stewart",
    )
    sup = Supplier.objects.create(name="Stewart Lab")  # NOT "Hannah Stewart"
    Compound.objects.create(target=t, smiles="CCO", registered_by=user)
    Compound.objects.create(target=t, smiles="CCN", supplier=sup)

    res = resolve_registrant("Stewart", field=FIELD_REGISTERED_BY)
    # Both match "stewart" (last name on user, first word on supplier).
    # No merge because full normalised names differ. Two chips.
    assert isinstance(res, UserClarify)
    union_chips = [c for c in res.candidates if isinstance(c, UnionCandidate)]
    assert len(union_chips) == 0


def test_user_linked_supplier_does_not_merge_via_slice_18(db):
    """Slice 16 already de-duped user-linked Suppliers via the user_id
    FK. Slice 18 should NOT also try to merge them — the linked
    case is handled upstream and the Supplier never reaches the
    same-name merger."""
    from compounds.registry.models import Supplier
    from compounds.nlp.spec import ResolvedUser
    t = Target.objects.create(name="ARd")
    bob = User.objects.create_user(
        username="bob.smith", email="bob.smith@ncl.ac.uk",
        first_name="Bob", last_name="Smith",
    )
    bob_sup = Supplier.objects.create(name="Bob Smith", user=bob)
    Compound.objects.create(target=t, smiles="CCO", registered_by=bob)
    Compound.objects.create(target=t, smiles="CCN", supplier=bob_sup)
    res = resolve_registrant("Bob Smith", field=FIELD_REGISTERED_BY)
    # Single User entity (Supplier is de-duped via user_id), not a Union.
    assert isinstance(res, ResolvedUser)
    assert res.user.pk == bob.pk


def test_pinned_union_round_trip(hannah_world):
    """Frontend continuation: clarify chip pinned both ids → resolver
    returns ResolvedUnion."""
    from compounds.nlp.spec import ResolvedUnion
    res = resolve_registrant(
        "ignored",
        field=FIELD_REGISTERED_BY,
        pinned_user_id=str(hannah_world["user"].pk),
        pinned_supplier_id=str(hannah_world["supplier"].pk),
    )
    assert isinstance(res, ResolvedUnion)
    assert res.user.pk == hannah_world["user"].pk
    assert res.supplier.pk == hannah_world["supplier"].pk
