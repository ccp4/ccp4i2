"""Selection-name resolver tests.

Maps user-typed phrases to UUIDs of the chemist's saved Selection rows
for the scatter view's colour-by overlay. Per-user scope (don't leak
across users), three-tier match (exact / substring / fuzzy), silently
drops phrases that don't match.
"""

from __future__ import annotations

import datetime

import pytest
from django.contrib.auth import get_user_model
from django.utils import timezone

from compounds.nlp.resolver import resolve_selection_names
from compounds.registry.models import Selection

User = get_user_model()


@pytest.fixture
def alice(db):
    return User.objects.create_user(username="alice", email="alice@example.org")


@pytest.fixture
def bob(db):
    return User.objects.create_user(username="bob", email="bob@example.org")


@pytest.fixture
def alice_selections(alice):
    a = Selection.objects.create(
        name="CDK4 scorecard top 20", compound_ids=["NCL-1"],
        created_by=alice, is_saved=True, expires_at=None,
    )
    b = Selection.objects.create(
        name="Mike's recent screen", compound_ids=["NCL-2"],
        created_by=alice, is_saved=True, expires_at=None,
    )
    c = Selection.objects.create(
        name="ARd analogues", compound_ids=["NCL-3"],
        created_by=alice,
        expires_at=timezone.now() + datetime.timedelta(days=1),
    )
    return a, b, c


# ---------------------------------------------------------------------------
# Empty / boundary inputs
# ---------------------------------------------------------------------------


def test_no_phrases_returns_empty(alice):
    assert resolve_selection_names([], alice) == []


def test_empty_strings_skipped(alice):
    assert resolve_selection_names(["", "   "], alice) == []


def test_anonymous_user_returns_empty(alice_selections):
    """Unauthenticated user can't have saved selections — no leak."""
    assert resolve_selection_names(["CDK4 scorecard top 20"], None) == []


def test_no_matching_selections_returns_empty(alice):
    assert resolve_selection_names(["doesn't exist"], alice) == []


# ---------------------------------------------------------------------------
# Exact match
# ---------------------------------------------------------------------------


def test_exact_match(alice, alice_selections):
    a, _b, _c = alice_selections
    assert resolve_selection_names(["CDK4 scorecard top 20"], alice) == [str(a.id)]


def test_exact_match_case_insensitive(alice, alice_selections):
    a, _b, _c = alice_selections
    assert resolve_selection_names(["cdk4 scorecard top 20"], alice) == [str(a.id)]


def test_exact_match_punctuation_insensitive(alice, alice_selections):
    """Match normalisation strips non-alphanumerics — apostrophes don't block."""
    _a, b, _c = alice_selections
    assert resolve_selection_names(["mikes recent screen"], alice) == [str(b.id)]


# ---------------------------------------------------------------------------
# Substring match
# ---------------------------------------------------------------------------


def test_substring_match_query_in_name(alice, alice_selections):
    """User's typed phrase is a substring of the selection name."""
    _a, _b, c = alice_selections
    assert resolve_selection_names(["analogues"], alice) == [str(c.id)]


def test_substring_match_name_in_query(alice):
    """Selection name is a substring of the user's typed phrase."""
    sel = Selection.objects.create(
        name="ARd", compound_ids=["NCL-1"],
        created_by=alice, is_saved=True, expires_at=None,
    )
    assert resolve_selection_names(["the ARd hits"], alice) == [str(sel.id)]


# ---------------------------------------------------------------------------
# Fuzzy match (Levenshtein-1)
# ---------------------------------------------------------------------------


def test_fuzzy_match_one_typo(alice):
    """A single-character typo still resolves."""
    sel = Selection.objects.create(
        name="scorecard", compound_ids=["NCL-1"],
        created_by=alice, is_saved=True, expires_at=None,
    )
    assert resolve_selection_names(["scorcard"], alice) == [str(sel.id)]


def test_fuzzy_skipped_below_min_length(alice):
    """Phrases below 4 chars (post-normalisation) skip the fuzzy tier
    to avoid promiscuous matching of short tokens."""
    Selection.objects.create(
        name="abc", compound_ids=["NCL-1"],
        created_by=alice, is_saved=True, expires_at=None,
    )
    # "axc" is one edit from "abc" but only 3 chars — below threshold.
    assert resolve_selection_names(["axc"], alice) == []


# ---------------------------------------------------------------------------
# Per-user scope
# ---------------------------------------------------------------------------


def test_does_not_match_other_users_selections(alice, bob, alice_selections):
    """Bob can't see Alice's selection, even by exact name."""
    assert resolve_selection_names(["CDK4 scorecard top 20"], bob) == []


# ---------------------------------------------------------------------------
# Multi-phrase + dedup
# ---------------------------------------------------------------------------


def test_multiple_phrases_match_in_order(alice, alice_selections):
    a, b, _c = alice_selections
    assert resolve_selection_names(
        ["CDK4 scorecard top 20", "Mike's recent screen"], alice,
    ) == [str(a.id), str(b.id)]


def test_dedup_when_phrases_match_same_selection(alice, alice_selections):
    a, _b, _c = alice_selections
    # Two phrases that both resolve to "a" — only one id in output.
    assert resolve_selection_names(
        ["CDK4 scorecard top 20", "scorecard"], alice,
    ) == [str(a.id)]


# ---------------------------------------------------------------------------
# Expiry
# ---------------------------------------------------------------------------


def test_excludes_expired_unsaved_selections(alice):
    Selection.objects.create(
        name="expired thing", compound_ids=["NCL-1"],
        created_by=alice, is_saved=False,
        expires_at=timezone.now() - datetime.timedelta(hours=1),
    )
    assert resolve_selection_names(["expired thing"], alice) == []


def test_includes_saved_with_stale_expires_at(alice):
    """A saved row with a stale expires_at is still resolvable —
    is_saved trumps expiry."""
    sel = Selection.objects.create(
        name="saved oddity", compound_ids=["NCL-1"],
        created_by=alice, is_saved=True,
        expires_at=timezone.now() - datetime.timedelta(days=10),
    )
    assert resolve_selection_names(["saved oddity"], alice) == [str(sel.id)]
