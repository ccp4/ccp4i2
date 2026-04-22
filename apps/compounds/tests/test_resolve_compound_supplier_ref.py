"""
Tests for resolve_compound's supplier_ref matching, including the
dash/space-tolerant fallback used for Pharmaron-style refs where users
often drop separators (e.g. "PH-NTU-00-00-000-0" typed as "PHNTU0000000").
"""

from __future__ import annotations

import pytest

from compounds.registry.models import Compound, Target
from compounds.utils import resolve_compound, resolve_compound_batch


@pytest.fixture
def pharmaron_compound(db):
    target = Target.objects.create(name='Pharmaron target')
    return Compound.objects.create(
        target=target,
        smiles='CCO',
        supplier_ref='PH-NTU-00-00-000-0',
    )


@pytest.mark.django_db
def test_exact_supplier_ref_matches(pharmaron_compound):
    assert resolve_compound('PH-NTU-00-00-000-0') == pharmaron_compound


@pytest.mark.django_db
def test_supplier_ref_case_insensitive(pharmaron_compound):
    assert resolve_compound('ph-ntu-00-00-000-0') == pharmaron_compound


@pytest.mark.parametrize('typed', [
    'PHNTU00000000',          # all dashes dropped (5 letters + 8 zeros)
    'phntu00000000',          # lowercase, no dashes
    'PH NTU 00 00 000 0',     # spaces instead of dashes
    'PH_NTU_00_00_000_0',     # underscores
    'PH/NTU/00/00/000/0',     # slashes
    'PH-NTU 00-00-000-0',     # mixed separators
])
@pytest.mark.django_db
def test_supplier_ref_tolerates_missing_or_substituted_separators(pharmaron_compound, typed):
    assert resolve_compound(typed) == pharmaron_compound


@pytest.mark.django_db
def test_unrelated_string_does_not_match(pharmaron_compound):
    assert resolve_compound('totally-unrelated-id') is None


@pytest.mark.django_db
def test_batch_resolver_applies_tolerant_supplier_ref_match(pharmaron_compound):
    results = resolve_compound_batch([
        'PH-NTU-00-00-000-0',   # exact
        'PHNTU00000000',        # dashes dropped
        'ph ntu 00 00 000 0',   # spaces, lowercase
        'unknown',              # should be None
    ])
    assert results['PH-NTU-00-00-000-0'] == pharmaron_compound
    assert results['PHNTU00000000'] == pharmaron_compound
    assert results['ph ntu 00 00 000 0'] == pharmaron_compound
    assert results['unknown'] is None
