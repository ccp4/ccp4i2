"""ScaffoldExtension priority + override tests (slice 17).

resolve_scaffold consults a unified pool: project-scoped extensions
first, then shared extensions, then the Python seed. This module pins
the priority semantics — same-name overrides, project-vs-shared,
substring/miss surfacing of DB rows alongside seed entries.
"""

from __future__ import annotations

import pytest

from compounds.nlp.resolver import resolve_scaffold
from compounds.nlp.spec import (
    FIELD_SCAFFOLD_HINT,
    ResolvedScaffold,
    ScaffoldClarify,
    ScaffoldMiss,
)
from compounds.registry.models import ScaffoldExtension, Target


pytestmark = pytest.mark.django_db


@pytest.fixture
def target_a(db):
    return Target.objects.create(name="ARd")


@pytest.fixture
def target_b(db):
    return Target.objects.create(name="CDK4")


# ---------------------------------------------------------------------------
# Shared catalog
# ---------------------------------------------------------------------------


def test_shared_extension_resolves_for_any_project(db, target_a):
    ScaffoldExtension.objects.create(
        name="indazole", smarts="c1ccc2[nH]ncc2c1", target=None,
    )
    res = resolve_scaffold("indazole", target=target_a)
    assert isinstance(res, ResolvedScaffold)
    assert res.scaffold.name == "indazole"
    assert res.scaffold.smarts == "c1ccc2[nH]ncc2c1"
    # Also resolves with no target (shared works everywhere).
    res2 = resolve_scaffold("indazole", target=None)
    assert isinstance(res2, ResolvedScaffold)


def test_shared_extension_aliases_resolve(db):
    ScaffoldExtension.objects.create(
        name="indazole", smarts="c1ccc2[nH]ncc2c1",
        aliases=["indazoles", "indazolyl"],
    )
    res = resolve_scaffold("indazoles")
    assert isinstance(res, ResolvedScaffold)
    assert res.scaffold.name == "indazole"


# ---------------------------------------------------------------------------
# Project-scoped catalog
# ---------------------------------------------------------------------------


def test_project_scoped_extension_only_resolves_under_its_target(
    db, target_a, target_b,
):
    ScaffoldExtension.objects.create(
        name="series 1", smarts="c1ccncc1", target=target_a,
    )
    # Resolves under ARd...
    res = resolve_scaffold("series 1", target=target_a)
    assert isinstance(res, ResolvedScaffold)
    # ...but not under CDK4 (different project).
    res_b = resolve_scaffold("series 1", target=target_b)
    assert isinstance(res_b, ScaffoldMiss)
    # ...and not when called with no target (shared+seed only).
    res_none = resolve_scaffold("series 1", target=None)
    assert isinstance(res_none, ScaffoldMiss)


def test_same_name_can_exist_in_two_projects(db, target_a, target_b):
    """'Series 1' in ARd and 'Series 1' in CDK4 are different fragments."""
    ScaffoldExtension.objects.create(
        name="series 1", smarts="c1ccncc1", target=target_a,
    )
    ScaffoldExtension.objects.create(
        name="series 1", smarts="c1ccccc1", target=target_b,
    )
    a = resolve_scaffold("series 1", target=target_a)
    b = resolve_scaffold("series 1", target=target_b)
    assert isinstance(a, ResolvedScaffold)
    assert isinstance(b, ResolvedScaffold)
    assert a.scaffold.smarts != b.scaffold.smarts


# ---------------------------------------------------------------------------
# Priority: project > shared > seed
# ---------------------------------------------------------------------------


def test_project_extension_overrides_shared(db, target_a):
    """When a name exists in both project and shared, the project entry
    wins — project-specific nomenclature beats global vocabulary."""
    ScaffoldExtension.objects.create(
        name="custom-frag", smarts="c1ccccc1", target=None,
    )
    ScaffoldExtension.objects.create(
        name="custom-frag", smarts="c1ccncc1", target=target_a,
    )
    res = resolve_scaffold("custom-frag", target=target_a)
    assert isinstance(res, ResolvedScaffold)
    # Project SMARTS, not shared.
    assert res.scaffold.smarts == "c1ccncc1"


def test_shared_extension_overrides_seed(db):
    """A shared extension named 'pyrimidine' would shadow the seed
    entry. (Curators take note — the seed is the canonical default,
    but extensions deliberately let you override locally.)"""
    ScaffoldExtension.objects.create(
        name="pyrimidine", smarts="c1cnnnn1",  # nonsense SMARTS for the test
    )
    res = resolve_scaffold("pyrimidine")
    assert isinstance(res, ResolvedScaffold)
    assert res.scaffold.smarts == "c1cnnnn1"   # extension wins, not seed's c1cncnc1


def test_seed_resolves_when_no_extension_shadows(db):
    """No extension with that name → fall through to the Python seed."""
    res = resolve_scaffold("pyrimidine")
    assert isinstance(res, ResolvedScaffold)
    assert res.scaffold.smarts == "c1cncnc1"   # seed value


# ---------------------------------------------------------------------------
# Substring + miss with extensions in the pool
# ---------------------------------------------------------------------------


def test_substring_match_finds_extension(db):
    ScaffoldExtension.objects.create(
        name="benzofuran", smarts="c1ccc2occc2c1",
    )
    # Substring "benzo" hits "benzofuran" only (no seed entry contains that).
    res = resolve_scaffold("benzofuran")
    assert isinstance(res, ResolvedScaffold)
    assert res.scaffold.name == "benzofuran"


def test_miss_surfaces_extension_in_suggestions(db):
    """Extensions participate in the miss-rank alongside seed entries.
    A query that's clearly closer to the extension than to any seed
    should put the extension at the top of suggestions."""
    ScaffoldExtension.objects.create(name="indazoline", smarts="C1CC2NCCc2N1")
    # 'indazoo' is a near-miss for 'indazoline' (starts the same) and
    # not for any seed entry.
    res = resolve_scaffold("indazoo")
    assert isinstance(res, ScaffoldMiss)
    names = {s.name for s in res.suggestions}
    assert "indazoline" in names


# ---------------------------------------------------------------------------
# Pinning round-trip with extensions
# ---------------------------------------------------------------------------


def test_pinned_extension_pk_resolves(db):
    ext = ScaffoldExtension.objects.create(
        name="indazole", smarts="c1ccc2[nH]ncc2c1",
    )
    res = resolve_scaffold("anything", pinned_id=str(ext.pk))
    assert isinstance(res, ResolvedScaffold)
    assert res.scaffold.name == "indazole"


def test_pinned_seed_canonical_name_still_resolves(db):
    """Pinning by a non-integer canonical name continues to hit the seed."""
    res = resolve_scaffold("anything", pinned_id="pyrimidine")
    assert isinstance(res, ResolvedScaffold)
    assert res.scaffold.name == "pyrimidine"
