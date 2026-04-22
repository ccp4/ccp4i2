"""
Tests for the hydrate_genes management command.

Uses pytest-django's `django_db` fixture; mocks HGNC by monkey-patching
`compounds.registry.hgnc.fetch_by_symbol`.
"""

from __future__ import annotations

from io import StringIO
from unittest import mock

import pytest
from django.core.management import call_command

from compounds.registry import hgnc
from compounds.registry.models import Gene


@pytest.fixture
def egfr_record():
    return hgnc.GeneRecord(
        symbol="EGFR",
        hgnc_id="HGNC:3236",
        name="epidermal growth factor receptor",
        aliases=["ERBB", "ERBB1", "HER1"],
        uniprot_ids=["P00533"],
        ensembl_gene_id="ENSG00000146648",
    )


@pytest.fixture
def patch_hgnc_fetch(egfr_record):
    """Patch fetch_by_symbol to return egfr_record for EGFR, None otherwise."""

    def _fetch(symbol, timeout=hgnc.DEFAULT_TIMEOUT):
        if symbol.upper() == "EGFR":
            return egfr_record
        return None

    with mock.patch.object(hgnc, "fetch_by_symbol", side_effect=_fetch) as m:
        yield m


@pytest.mark.django_db
def test_hydrate_pending_gene_fills_fields(patch_hgnc_fetch):
    gene = Gene.objects.create(symbol="EGFR", hydration_status="pending")
    out = StringIO()

    call_command("hydrate_genes", stdout=out, stderr=StringIO())

    gene.refresh_from_db()
    assert gene.hgnc_id == "HGNC:3236"
    assert gene.name == "epidermal growth factor receptor"
    assert gene.aliases == ["ERBB", "ERBB1", "HER1"]
    assert gene.uniprot_ids == ["P00533"]
    assert gene.ensembl_gene_id == "ENSG00000146648"
    assert gene.hydration_status == "ok"
    assert gene.hydrated_at is not None
    assert gene.hydration_source == "hgnc-api"


@pytest.mark.django_db
def test_hydrate_unknown_symbol_marks_failed(patch_hgnc_fetch):
    gene = Gene.objects.create(symbol="NOPE", hydration_status="pending")
    call_command("hydrate_genes", stdout=StringIO(), stderr=StringIO())

    gene.refresh_from_db()
    assert gene.hydration_status == "failed"
    assert gene.hydrated_at is not None


@pytest.mark.django_db
def test_hydrate_dry_run_does_not_persist(patch_hgnc_fetch):
    gene = Gene.objects.create(symbol="EGFR", hydration_status="pending")
    call_command("hydrate_genes", "--dry-run", stdout=StringIO(), stderr=StringIO())

    gene.refresh_from_db()
    assert gene.hydration_status == "pending"
    assert gene.hgnc_id == ""
    assert gene.aliases == []


@pytest.mark.django_db
def test_hydrate_only_symbol_filters(patch_hgnc_fetch):
    Gene.objects.create(symbol="EGFR", hydration_status="pending")
    other = Gene.objects.create(symbol="AKT1", hydration_status="pending")

    call_command(
        "hydrate_genes", "--symbol", "egfr", stdout=StringIO(), stderr=StringIO()
    )

    other.refresh_from_db()
    # AKT1 was pending but --symbol scoped work to EGFR only
    assert other.hydration_status == "pending"


@pytest.mark.django_db
def test_hydrate_hgnc_error_marks_failed(egfr_record):
    gene = Gene.objects.create(symbol="EGFR", hydration_status="pending")

    def _raise(symbol, timeout=hgnc.DEFAULT_TIMEOUT):
        raise hgnc.HGNCError("simulated network failure")

    with mock.patch.object(hgnc, "fetch_by_symbol", side_effect=_raise):
        call_command("hydrate_genes", stdout=StringIO(), stderr=StringIO())

    gene.refresh_from_db()
    assert gene.hydration_status == "failed"


@pytest.mark.django_db
def test_hydrate_skips_fresh_ok_genes_by_default(patch_hgnc_fetch):
    """A Gene already 'ok' with a recent hydrated_at should be skipped."""
    from django.utils import timezone

    gene = Gene.objects.create(
        symbol="EGFR",
        hydration_status="ok",
        hydrated_at=timezone.now(),
        hgnc_id="HGNC:3236",
        name="existing name",
    )

    call_command("hydrate_genes", stdout=StringIO(), stderr=StringIO())

    # fetch_by_symbol should not have been called for this gene
    assert patch_hgnc_fetch.call_count == 0
    gene.refresh_from_db()
    assert gene.name == "existing name"  # untouched
