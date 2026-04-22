"""
Tests for TargetSerializer's gene_symbols write field and the genes
nested read field. See TARGET_MODEL_PROPOSAL.md §7.
"""

from __future__ import annotations

import pytest

from compounds.registry.models import Gene, Target
from compounds.registry.serializers import TargetSerializer


@pytest.mark.django_db
def test_create_target_with_gene_symbols_creates_pending_genes():
    serializer = TargetSerializer(data={
        'name': 'EGFR programme',
        'gene_symbols': ['EGFR', 'ERBB2'],
    })
    assert serializer.is_valid(), serializer.errors
    target = serializer.save()

    symbols = set(target.genes.values_list('symbol', flat=True))
    assert symbols == {'EGFR', 'ERBB2'}

    egfr = Gene.objects.get(symbol='EGFR')
    assert egfr.hydration_status == 'pending'
    assert egfr.aliases == []


@pytest.mark.django_db
def test_create_target_with_gene_symbols_reuses_existing_genes():
    pre = Gene.objects.create(
        symbol='EGFR',
        hydration_status='ok',
        name='epidermal growth factor receptor',
        aliases=['HER1'],
    )
    serializer = TargetSerializer(data={
        'name': 'AR degraders',
        'gene_symbols': ['EGFR'],
    })
    assert serializer.is_valid(), serializer.errors
    target = serializer.save()

    assert target.genes.count() == 1
    assert target.genes.first().pk == pre.pk
    # Existing row not overwritten
    pre.refresh_from_db()
    assert pre.name == 'epidermal growth factor receptor'


@pytest.mark.django_db
def test_create_target_normalises_symbols_to_upper():
    serializer = TargetSerializer(data={
        'name': 'normalisation target',
        'gene_symbols': ['egfr', ' akt1 ', ''],  # mixed case, whitespace, empty
    })
    assert serializer.is_valid(), serializer.errors
    target = serializer.save()

    symbols = set(target.genes.values_list('symbol', flat=True))
    assert symbols == {'EGFR', 'AKT1'}  # empty filtered, others normalised


@pytest.mark.django_db
def test_update_target_gene_symbols_replaces_set():
    gene_a = Gene.objects.create(symbol='A', hydration_status='ok')
    target = Target.objects.create(name='Replaceable')
    target.genes.add(gene_a)

    serializer = TargetSerializer(
        instance=target,
        data={'gene_symbols': ['B', 'C']},
        partial=True,
    )
    assert serializer.is_valid(), serializer.errors
    serializer.save()

    symbols = set(target.genes.values_list('symbol', flat=True))
    assert symbols == {'B', 'C'}  # A dropped, B and C added


@pytest.mark.django_db
def test_update_target_omits_gene_symbols_leaves_genes_alone():
    gene_a = Gene.objects.create(symbol='A', hydration_status='ok')
    target = Target.objects.create(name='Sticky')
    target.genes.add(gene_a)

    serializer = TargetSerializer(
        instance=target,
        data={'name': 'Sticky renamed'},
        partial=True,
    )
    assert serializer.is_valid(), serializer.errors
    serializer.save()

    assert target.genes.count() == 1
    assert target.genes.first().symbol == 'A'


@pytest.mark.django_db
def test_update_target_empty_list_clears_genes():
    gene = Gene.objects.create(symbol='X', hydration_status='ok')
    target = Target.objects.create(name='Clearable')
    target.genes.add(gene)

    serializer = TargetSerializer(
        instance=target,
        data={'gene_symbols': []},
        partial=True,
    )
    assert serializer.is_valid(), serializer.errors
    serializer.save()

    assert target.genes.count() == 0


@pytest.mark.django_db
def test_read_target_includes_nested_genes():
    gene = Gene.objects.create(
        symbol='EGFR',
        hydration_status='ok',
        hgnc_id='HGNC:3236',
        name='epidermal growth factor receptor',
        aliases=['ERBB', 'HER1'],
    )
    target = Target.objects.create(name='with-gene')
    target.genes.add(gene)

    data = TargetSerializer(target).data
    assert len(data['genes']) == 1
    assert data['genes'][0]['symbol'] == 'EGFR'
    assert data['genes'][0]['hgnc_id'] == 'HGNC:3236'
    assert data['genes'][0]['aliases'] == ['ERBB', 'HER1']
    # gene_symbols is write_only so shouldn't appear in output
    assert 'gene_symbols' not in data
