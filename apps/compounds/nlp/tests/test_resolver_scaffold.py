"""Scaffold resolver tests (§19.3 slice 13).

Deterministic name-to-SMARTS mapping from ``compounds.nlp.substructures``.
Three-tier match: exact → substring → fuzzy miss with suggestions.
Pinned-id bypasses the fuzzy path.
"""

from __future__ import annotations

from compounds.nlp.resolver import resolve_scaffold
from compounds.nlp.spec import (
    FIELD_SCAFFOLD_HINT,
    ResolvedScaffold,
    ScaffoldClarify,
    ScaffoldMiss,
)


# ---------------------------------------------------------------------------
# Exact match (tier 1)
# ---------------------------------------------------------------------------


def test_exact_canonical_name_resolves():
    res = resolve_scaffold("pyrimidine")
    assert isinstance(res, ResolvedScaffold)
    assert res.scaffold.name == "pyrimidine"
    assert res.matched_via == "name"


def test_exact_alias_resolves():
    # "pyridinyl" is an alias of "pyridine".
    res = resolve_scaffold("pyridinyl")
    assert isinstance(res, ResolvedScaffold)
    assert res.scaffold.name == "pyridine"
    assert res.matched_via == "alias"


def test_plural_alias_resolves():
    # "pyrimidines" → pyrimidine.
    res = resolve_scaffold("pyrimidines")
    assert isinstance(res, ResolvedScaffold)
    assert res.scaffold.name == "pyrimidine"


def test_normalisation_ignores_case_and_punctuation():
    # Mixed case + hyphens + whitespace should all normalise away.
    res = resolve_scaffold(" Pyrimidine ")
    assert isinstance(res, ResolvedScaffold)
    assert res.scaffold.name == "pyrimidine"


# ---------------------------------------------------------------------------
# Substring match (tier 2)
# ---------------------------------------------------------------------------


def test_substring_single_hit_resolves():
    # "pyrimidi" is a substring of exactly "pyrimidine" / "pyrimidines" /
    # "pyrimidinyl" — all the same scaffold — so should resolve cleanly.
    res = resolve_scaffold("pyrimidi")
    assert isinstance(res, ResolvedScaffold)
    assert res.scaffold.name == "pyrimidine"


def test_substring_multiple_hits_clarifies():
    # "amid" is a substring of both "amide" and "sulfonamide" — ambiguous.
    res = resolve_scaffold("amid")
    assert isinstance(res, ScaffoldClarify)
    names = {c.name for c in res.candidates}
    assert "amide" in names
    assert "sulfonamide" in names


def test_clarify_tags_scaffold_index_and_field():
    res = resolve_scaffold("amid", scaffold_index=2)
    assert isinstance(res, ScaffoldClarify)
    assert res.scaffold_index == 2
    assert res.field == FIELD_SCAFFOLD_HINT


# ---------------------------------------------------------------------------
# Miss (tier 3)
# ---------------------------------------------------------------------------


def test_unknown_name_misses_with_suggestions():
    res = resolve_scaffold("zzzzznothing")
    assert isinstance(res, ScaffoldMiss)
    assert len(res.suggestions) > 0


def test_empty_query_misses():
    res = resolve_scaffold("")
    assert isinstance(res, ScaffoldMiss)


def test_miss_tags_scaffold_index_and_field():
    res = resolve_scaffold("zzzzznothing", scaffold_index=1)
    assert isinstance(res, ScaffoldMiss)
    assert res.scaffold_index == 1
    assert res.field == FIELD_SCAFFOLD_HINT


# ---------------------------------------------------------------------------
# Pinned id — clarify continuation
# ---------------------------------------------------------------------------


def test_pinned_id_bypasses_fuzzy_match():
    # Even though "amid" is ambiguous, supplying a pinned id resolves
    # directly to that scaffold.
    res = resolve_scaffold("amid", pinned_id="sulfonamide")
    assert isinstance(res, ResolvedScaffold)
    assert res.scaffold.name == "sulfonamide"
    assert res.matched_via == "pinned_id"


def test_pinned_id_unknown_misses():
    res = resolve_scaffold("pyridine", pinned_id="not-a-real-name")
    assert isinstance(res, ScaffoldMiss)
