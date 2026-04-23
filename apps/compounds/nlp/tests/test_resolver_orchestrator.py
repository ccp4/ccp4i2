"""Two-target orchestrator tests.

Covers §6.5 (scope derivation from which target fields are populated) and
Q21 (sequential clarification — registration field first, then assay field).
"""

from __future__ import annotations

import pytest

from compounds.nlp.resolver import resolve_targets
from compounds.nlp.spec import (
    FIELD_ASSAY_TARGET,
    FIELD_REGISTRATION_TARGET,
    SCOPE_ASSAY_ONLY,
    SCOPE_BOTH_SAME,
    SCOPE_CROSS,
    SCOPE_REG_ONLY,
    CompoundSelector,
    ResolvedTargets,
    ScopeError,
    TargetClarify,
    TargetMiss,
)
from compounds.registry.models import Gene, Target


@pytest.fixture
def two_target_world(db):
    egfr_gene = Gene.objects.create(symbol="EGFR", hydration_status="ok")
    akt_gene = Gene.objects.create(symbol="AKT1", hydration_status="ok")
    myc_gene = Gene.objects.create(symbol="MYC", hydration_status="ok")

    ar = Target.objects.create(name="AR degraders")  # name-only, no genes
    egfr = Target.objects.create(name="EGFR programme")
    egfr.genes.add(egfr_gene)
    akt = Target.objects.create(name="AKT programme")
    akt.genes.add(akt_gene)
    myc_aur = Target.objects.create(name="Myc-Aur")
    myc_aur.genes.add(myc_gene)
    myc_rna = Target.objects.create(name="Myc RNA")
    myc_rna.genes.add(myc_gene)

    return {
        "ar": ar,
        "egfr": egfr,
        "akt": akt,
        "myc_aur": myc_aur,
        "myc_rna": myc_rna,
    }


# ---------------------------------------------------------------------------
# Scope-kind derivation (§6.5 table, rows 1–4)
# ---------------------------------------------------------------------------


def test_both_same_when_fields_equal(two_target_world):
    spec = CompoundSelector(
        registration_target_as_typed="AR degraders",
        assay_target_as_typed="AR degraders",
    )
    out = resolve_targets(spec)
    assert isinstance(out, ResolvedTargets)
    assert out.registration.pk == two_target_world["ar"].pk
    assert out.assay.pk == two_target_world["ar"].pk
    assert out.scope_kind == SCOPE_BOTH_SAME


def test_reg_only_when_assay_empty(two_target_world):
    spec = CompoundSelector(registration_target_as_typed="EGFR")
    out = resolve_targets(spec)
    assert isinstance(out, ResolvedTargets)
    assert out.registration.pk == two_target_world["egfr"].pk
    assert out.assay is None
    assert out.scope_kind == SCOPE_REG_ONLY


def test_assay_only_when_reg_empty(two_target_world):
    spec = CompoundSelector(assay_target_as_typed="AKT1")
    out = resolve_targets(spec)
    assert isinstance(out, ResolvedTargets)
    assert out.registration is None
    assert out.assay.pk == two_target_world["akt"].pk
    assert out.scope_kind == SCOPE_ASSAY_ONLY


def test_cross_when_fields_differ(two_target_world):
    spec = CompoundSelector(
        registration_target_as_typed="EGFR",
        assay_target_as_typed="AKT1",
    )
    out = resolve_targets(spec)
    assert isinstance(out, ResolvedTargets)
    assert out.registration.pk == two_target_world["egfr"].pk
    assert out.assay.pk == two_target_world["akt"].pk
    assert out.scope_kind == SCOPE_CROSS


def test_both_fields_empty_is_error(two_target_world):
    spec = CompoundSelector()
    out = resolve_targets(spec)
    assert isinstance(out, ScopeError)


def test_both_fields_whitespace_is_error(two_target_world):
    spec = CompoundSelector(
        registration_target_as_typed="   ",
        assay_target_as_typed="",
    )
    out = resolve_targets(spec)
    assert isinstance(out, ScopeError)


# ---------------------------------------------------------------------------
# Q21 — sequential clarification: reg first, then assay
# ---------------------------------------------------------------------------


def test_registration_ambiguity_short_circuits_before_assay(two_target_world):
    # MYC resolves to two targets (Myc-Aur and Myc RNA). If reg is
    # ambiguous, we return a TargetClarify tagged with the registration
    # field and stop — we don't even try to resolve assay yet.
    spec = CompoundSelector(
        registration_target_as_typed="MYC",
        assay_target_as_typed="also_ambiguous_or_whatever",  # never reached
    )
    out = resolve_targets(spec)
    assert isinstance(out, TargetClarify)
    assert out.field == FIELD_REGISTRATION_TARGET
    assert len(out.candidates) == 2


def test_assay_ambiguity_surfaces_only_when_reg_is_clean(two_target_world):
    spec = CompoundSelector(
        registration_target_as_typed="EGFR",
        assay_target_as_typed="MYC",
    )
    out = resolve_targets(spec)
    assert isinstance(out, TargetClarify)
    assert out.field == FIELD_ASSAY_TARGET
    assert len(out.candidates) == 2


def test_registration_miss_short_circuits_before_assay(two_target_world):
    spec = CompoundSelector(
        registration_target_as_typed="nonsense_target_name",
        assay_target_as_typed="EGFR",  # would resolve fine, but never tried
    )
    out = resolve_targets(spec)
    assert isinstance(out, TargetMiss)
    assert out.field == FIELD_REGISTRATION_TARGET


def test_assay_miss_surfaces_only_when_reg_is_clean(two_target_world):
    spec = CompoundSelector(
        registration_target_as_typed="EGFR",
        assay_target_as_typed="nonsense_target_name",
    )
    out = resolve_targets(spec)
    assert isinstance(out, TargetMiss)
    assert out.field == FIELD_ASSAY_TARGET


# ---------------------------------------------------------------------------
# Field-name tagging is reliable for UI dispatch
# ---------------------------------------------------------------------------


def test_clarify_field_tag_survives_resolution_path(two_target_world):
    spec = CompoundSelector(assay_target_as_typed="MYC")
    out = resolve_targets(spec)
    assert isinstance(out, TargetClarify)
    assert out.field == FIELD_ASSAY_TARGET
    # Candidates still populated; the field tag is additive.
    assert len(out.candidates) == 2
