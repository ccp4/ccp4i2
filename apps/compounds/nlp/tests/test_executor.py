"""Executor tests — composes resolvers + evaluator end-to-end (§7).

Covers spec validation, pass-through of clarify/miss/scope errors,
phys-chem column expansion (§6.6), row classification into matched /
excluded-footer / filtered-silent buckets, and scope-sentence generation
(§10 echo-back) across all four scope kinds.
"""

from __future__ import annotations

import pytest

from compounds.assays.models import AnalysisResult, Assay, DataSeries, Protocol
from compounds.nlp.executor import execute
from compounds.nlp.spec import (
    COL_PRESET_LIPINSKI,
    COL_PRESET_PHYS_CHEM,
    EXCLUDE_UNIT_UNKNOWN,
    FILTER_KPI_MISMATCH,
    FILTER_THRESHOLD_NOT_MET,
    FIELD_PROTOCOL_HINT,
    FIELD_REGISTRATION_TARGET,
    LIPINSKI_COLUMNS,
    PHYS_CHEM_COLUMNS,
    ProtocolClarify,
    ProtocolMiss,
    QuerySpec,
    ScopeError,
    SpecError,
    TablePayload,
    TargetClarify,
    TargetMiss,
    Threshold,
)
from compounds.registry.models import Compound, Gene, MolecularProperties, Target


# ---------------------------------------------------------------------------
# Fixture helpers
# ---------------------------------------------------------------------------


def _mk_compound(target, smiles, props=None) -> Compound:
    # Compound.save() auto-calculates MolecularProperties from SMILES via a
    # signal — we override its values or delete the row depending on what the
    # test needs.
    c = Compound.objects.create(target=target, smiles=smiles)
    if props is not None:
        MolecularProperties.objects.filter(compound=c).update(**props)
    else:
        MolecularProperties.objects.filter(compound=c).delete()
    return c


def _mk_result(assay, compound, *, status, results, row=0) -> AnalysisResult:
    ar = AnalysisResult.objects.create(status=status, results=results)
    DataSeries.objects.create(
        assay=assay,
        compound=compound,
        row=row,
        start_column=0,
        end_column=10,
        analysis=ar,
    )
    return ar


@pytest.fixture
def ar_world(db):
    """A compact AR programme catalog with one canonical HTRF protocol."""
    ar = Target.objects.create(name="AR degraders")
    akt = Target.objects.create(name="AKT programme")

    c1 = _mk_compound(
        ar,
        "CCO",
        {
            "molecular_weight": 400.5,
            "heavy_atom_count": 28,
            "hbd": 2,
            "hba": 5,
            "clogp": 3.2,
            "tpsa": 75.1,
            "rotatable_bonds": 6,
            "fraction_sp3": 0.35,
        },
    )
    c2 = _mk_compound(ar, "CCN", None)  # no MolecularProperties — null handling

    protocol = Protocol.objects.create(name="AR binding HTRF")
    assay = Assay.objects.create(protocol=protocol, target=ar)

    _mk_result(
        assay, c1,
        status="valid",
        results={"KPI": "IC50", "IC50": 5.0, "kpi_unit": "nM"},
        row=0,
    )
    _mk_result(
        assay, c2,
        status="valid",
        results={"KPI": "IC50", "IC50": 200.0, "kpi_unit": "nM"},
        row=1,
    )
    # Pre-filtered out at the DB level (status != valid).
    _mk_result(
        assay, c1,
        status="invalid",
        results={"KPI": "IC50", "IC50": 99.0, "kpi_unit": "nM"},
        row=2,
    )
    # Filtered silently — KPI mismatch.
    _mk_result(
        assay, c1,
        status="valid",
        results={"KPI": "EC50", "EC50": 12.0, "kpi_unit": "nM"},
        row=3,
    )
    # Excluded (footer) — missing unit under a real-unit threshold.
    _mk_result(
        assay, c1,
        status="valid",
        results={"KPI": "IC50", "IC50": 7.0, "kpi_unit": ""},
        row=4,
    )

    return {
        "ar": ar,
        "akt": akt,
        "c1": c1,
        "c2": c2,
        "protocol": protocol,
        "assay": assay,
    }


# ---------------------------------------------------------------------------
# Spec-level validation (no DB needed)
# ---------------------------------------------------------------------------


def test_missing_metric_returns_spec_error():
    res = execute(QuerySpec(
        registration_target_as_typed="AR degraders",
        protocol_hint="HTRF",
    ))
    assert isinstance(res, SpecError)
    assert res.field == "metric"


def test_missing_protocol_hint_returns_spec_error():
    res = execute(QuerySpec(
        registration_target_as_typed="AR degraders",
        metric="IC50",
    ))
    assert isinstance(res, SpecError)
    assert res.field == "protocol_hint"


# ---------------------------------------------------------------------------
# Pass-through of resolver responses
# ---------------------------------------------------------------------------


def test_scope_error_passthrough(db):
    res = execute(QuerySpec(metric="IC50", protocol_hint="HTRF"))
    assert isinstance(res, ScopeError)


def test_target_miss_passthrough(ar_world):
    res = execute(QuerySpec(
        registration_target_as_typed="zzzzznotarealtarget",
        metric="IC50",
        protocol_hint="HTRF",
    ))
    assert isinstance(res, TargetMiss)
    assert res.field == FIELD_REGISTRATION_TARGET


def test_protocol_miss_passthrough(ar_world):
    res = execute(QuerySpec(
        registration_target_as_typed="AR degraders",
        metric="IC50",
        protocol_hint="AlphaLISA",
    ))
    assert isinstance(res, ProtocolMiss)
    assert res.field == FIELD_PROTOCOL_HINT


def test_target_clarify_passthrough(db):
    # Two targets sharing a gene symbol — typing the shared symbol must hit
    # both pools and produce Clarify (matches the "typing MYC hits Myc-Aur
    # and Myc RNA" example from §6.1).
    myc = Gene.objects.create(symbol="MYC")
    t1 = Target.objects.create(name="Myc-Aur")
    t2 = Target.objects.create(name="Myc regulation")
    t1.genes.add(myc)
    t2.genes.add(myc)
    res = execute(QuerySpec(
        registration_target_as_typed="MYC",
        metric="IC50",
        protocol_hint="HTRF",
    ))
    assert isinstance(res, TargetClarify)
    assert res.field == FIELD_REGISTRATION_TARGET


def test_protocol_clarify_passthrough(db):
    ar = Target.objects.create(name="AR degraders")
    c = Compound.objects.create(target=ar, smiles="CCO")
    for name in ("AR binding HTRF", "AR degradation HTRF"):
        p = Protocol.objects.create(name=name)
        a = Assay.objects.create(protocol=p, target=ar)
        DataSeries.objects.create(
            assay=a, compound=c, row=0, start_column=0, end_column=10,
        )
    res = execute(QuerySpec(
        registration_target_as_typed="AR degraders",
        metric="IC50",
        protocol_hint="HTRF",
    ))
    assert isinstance(res, ProtocolClarify)
    assert {c.name for c in res.candidates} == {"AR binding HTRF", "AR degradation HTRF"}


# ---------------------------------------------------------------------------
# Happy path — row classification + properties pluck
# ---------------------------------------------------------------------------


def test_happy_path_single_compound_below_threshold(ar_world):
    res = execute(QuerySpec(
        registration_target_as_typed="AR degraders",
        metric="IC50",
        protocol_hint="HTRF",
        threshold=Threshold(op="<", value=10.0, unit="nM"),
    ))
    assert isinstance(res, TablePayload)
    assert [r.formatted_id for r in res.rows] == [ar_world["c1"].formatted_id]
    (matched,) = res.rows
    assert matched.value == 5.0
    assert matched.value_unit == "nM"
    assert matched.value_in_query_unit == 5.0
    # Phys-chem columns populated from MolecularProperties.
    assert matched.properties["molecular_weight"] == pytest.approx(400.5)
    assert matched.properties["clogp"] == pytest.approx(3.2)


def test_filtered_and_excluded_buckets(ar_world):
    res = execute(QuerySpec(
        registration_target_as_typed="AR degraders",
        metric="IC50",
        protocol_hint="HTRF",
        threshold=Threshold(op="<", value=10.0, unit="nM"),
    ))
    assert isinstance(res, TablePayload)
    # c1 passes, c2 fails threshold (200 >= 10) → filtered_silent[THRESHOLD_NOT_MET] == 1
    assert res.filtered_silent.get(FILTER_THRESHOLD_NOT_MET) == 1
    # One row with KPI='EC50' → filtered_silent[KPI_MISMATCH] == 1
    assert res.filtered_silent.get(FILTER_KPI_MISMATCH) == 1
    # One row with kpi_unit='' → footer_excluded[UNIT_UNKNOWN] == 1
    assert res.footer_excluded.get(EXCLUDE_UNIT_UNKNOWN) == 1


def test_null_molecular_properties_emits_none_cells(ar_world):
    # Drop the threshold so c2 passes through even at IC50=200 nM.
    res = execute(QuerySpec(
        registration_target_as_typed="AR degraders",
        metric="IC50",
        protocol_hint="HTRF",
    ))
    assert isinstance(res, TablePayload)
    c2_row = next(r for r in res.rows if r.formatted_id == ar_world["c2"].formatted_id)
    assert all(v is None for v in c2_row.properties.values())


# ---------------------------------------------------------------------------
# Column expansion (§6.6)
# ---------------------------------------------------------------------------


def test_default_columns_phys_chem(ar_world):
    res = execute(QuerySpec(
        registration_target_as_typed="AR degraders",
        metric="IC50",
        protocol_hint="HTRF",
    ))
    assert isinstance(res, TablePayload)
    assert res.property_columns == list(PHYS_CHEM_COLUMNS)


def test_lipinski_preset(ar_world):
    res = execute(QuerySpec(
        registration_target_as_typed="AR degraders",
        metric="IC50",
        protocol_hint="HTRF",
        columns=COL_PRESET_LIPINSKI,
    ))
    assert isinstance(res, TablePayload)
    assert res.property_columns == list(LIPINSKI_COLUMNS)


def test_columns_explicit_wins(ar_world):
    res = execute(QuerySpec(
        registration_target_as_typed="AR degraders",
        metric="IC50",
        protocol_hint="HTRF",
        columns=COL_PRESET_PHYS_CHEM,
        columns_explicit=["clogp", "molecular_weight"],
    ))
    assert isinstance(res, TablePayload)
    assert res.property_columns == ["clogp", "molecular_weight"]


def test_columns_explicit_silently_drops_unknowns(ar_world):
    res = execute(QuerySpec(
        registration_target_as_typed="AR degraders",
        metric="IC50",
        protocol_hint="HTRF",
        columns_explicit=["clogp", "nonsense_property"],
    ))
    assert isinstance(res, TablePayload)
    assert res.property_columns == ["clogp"]


def test_columns_explicit_all_unknown_falls_back_to_phys_chem(ar_world):
    res = execute(QuerySpec(
        registration_target_as_typed="AR degraders",
        metric="IC50",
        protocol_hint="HTRF",
        columns_explicit=["nonsense_a", "nonsense_b"],
    ))
    assert isinstance(res, TablePayload)
    assert res.property_columns == list(PHYS_CHEM_COLUMNS)


# ---------------------------------------------------------------------------
# Scope sentence (§10)
# ---------------------------------------------------------------------------


def test_scope_sentence_both_same_with_threshold(ar_world):
    res = execute(QuerySpec(
        registration_target_as_typed="AR degraders",
        assay_target_as_typed="AR degraders",
        metric="IC50",
        protocol_hint="HTRF",
        threshold=Threshold(op="<", value=10.0, unit="nM"),
    ))
    assert isinstance(res, TablePayload)
    assert res.scope_sentence == (
        "Showing compounds registered to AR degraders and tested in "
        "AR binding HTRF against AR degraders, with IC50 < 10.0 nM"
    )


def test_scope_sentence_reg_only_no_threshold(ar_world):
    res = execute(QuerySpec(
        registration_target_as_typed="AR degraders",
        metric="IC50",
        protocol_hint="HTRF",
    ))
    assert isinstance(res, TablePayload)
    assert res.scope_sentence == (
        "Showing compounds registered to AR degraders, tested in AR binding HTRF"
        ", IC50 values"
    )


def test_scope_sentence_assay_only(ar_world):
    res = execute(QuerySpec(
        assay_target_as_typed="AR degraders",
        metric="IC50",
        protocol_hint="HTRF",
    ))
    assert isinstance(res, TablePayload)
    assert res.scope_sentence == (
        "Showing compounds tested in AR binding HTRF against AR degraders"
        ", IC50 values"
    )


def test_scope_sentence_cross_project(db):
    ar = Target.objects.create(name="AR degraders")
    akt = Target.objects.create(name="AKT programme")
    c = Compound.objects.create(target=ar, smiles="CCO")
    protocol = Protocol.objects.create(name="AKT selectivity HTRF")
    assay = Assay.objects.create(protocol=protocol, target=akt)
    _mk_result(
        assay, c,
        status="valid",
        results={"KPI": "IC50", "IC50": 150.0, "kpi_unit": "nM"},
    )
    res = execute(QuerySpec(
        registration_target_as_typed="AR degraders",
        assay_target_as_typed="AKT programme",
        metric="IC50",
        protocol_hint="HTRF",
    ))
    assert isinstance(res, TablePayload)
    assert res.scope_sentence == (
        "Showing compounds registered to AR degraders, tested in "
        "AKT selectivity HTRF against AKT programme, IC50 values"
    )
