"""Compound-reference resolver tests.

Deterministic parsing via ``compounds.formatting.extract_reg_number`` —
no clarify path, only Resolved / Miss.
"""

from __future__ import annotations

import pytest

from compounds.nlp.resolver import resolve_compound_ref
from compounds.nlp.spec import (
    FIELD_COMPOUND_REF,
    CompoundMiss,
    ResolvedCompound,
)
from compounds.registry.models import Compound, Target


@pytest.fixture
def compound_world(db):
    t = Target.objects.create(name="ARd")
    c = Compound.objects.create(target=t, smiles="CCO", reg_number=26007)
    return {"t": t, "c": c}


def test_full_id_resolves(compound_world):
    res = resolve_compound_ref("NCL-00026007")
    assert isinstance(res, ResolvedCompound)
    assert res.compound.pk == compound_world["c"].pk


def test_compact_id_resolves(compound_world):
    res = resolve_compound_ref("NCL26007")
    assert isinstance(res, ResolvedCompound)
    assert res.compound.pk == compound_world["c"].pk


def test_lowercase_prefix_resolves(compound_world):
    res = resolve_compound_ref("ncl-00026007")
    assert isinstance(res, ResolvedCompound)
    assert res.compound.pk == compound_world["c"].pk


def test_bare_number_resolves(compound_world):
    res = resolve_compound_ref("26007")
    assert isinstance(res, ResolvedCompound)
    assert res.compound.pk == compound_world["c"].pk


def test_unparseable_string_misses(compound_world):
    res = resolve_compound_ref("the lead compound")
    assert isinstance(res, CompoundMiss)
    assert res.field == FIELD_COMPOUND_REF


def test_nonexistent_id_misses(compound_world):
    res = resolve_compound_ref("NCL-00099999")
    assert isinstance(res, CompoundMiss)
    assert res.field == FIELD_COMPOUND_REF


def test_empty_string_misses(compound_world):
    res = resolve_compound_ref("")
    assert isinstance(res, CompoundMiss)


def test_miss_carries_ref_index(compound_world):
    res = resolve_compound_ref("not-a-compound", ref_index=2)
    assert isinstance(res, CompoundMiss)
    assert res.ref_index == 2
