"""Executor tests (post-pivot).

Composes resolvers + evaluator end-to-end and verifies the compound-
intersection semantics introduced by the pivot: filters AND together,
each filter contributes its protocol name to the redirect URL, scope
sentence reflects all filters.
"""

from __future__ import annotations

import pytest

from compounds.assays.models import AnalysisResult, Assay, DataSeries, Protocol
from compounds.nlp.executor import execute
from compounds.nlp.spec import (
    CompoundSelector,
    CompoundSelection,
    FIELD_PROTOCOL_HINT,
    FIELD_REGISTRATION_TARGET,
    MeasurementFilter,
    ProtocolClarify,
    ProtocolMiss,
    ScopeError,
    SpecError,
    TargetClarify,
    TargetMiss,
    Threshold,
)
from compounds.registry.models import Compound, Gene, Target


# ---------------------------------------------------------------------------
# Fixture helpers
# ---------------------------------------------------------------------------


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
    """AR programme with one HTRF protocol + three compounds with measurements."""
    ar = Target.objects.create(name="AR degraders")
    c1 = Compound.objects.create(target=ar, smiles="CCO")
    c2 = Compound.objects.create(target=ar, smiles="CCN")
    c3 = Compound.objects.create(target=ar, smiles="CCC")

    protocol = Protocol.objects.create(name="AR binding HTRF")
    assay = Assay.objects.create(protocol=protocol, target=ar)

    _mk_result(assay, c1, status="valid",
               results={"KPI": "IC50", "IC50": 5.0, "kpi_unit": "nM"}, row=0)
    _mk_result(assay, c2, status="valid",
               results={"KPI": "IC50", "IC50": 200.0, "kpi_unit": "nM"}, row=1)
    _mk_result(assay, c3, status="valid",
               results={"KPI": "IC50", "IC50": 50.0, "kpi_unit": "nM"}, row=2)

    return {"ar": ar, "c1": c1, "c2": c2, "c3": c3, "protocol": protocol}


@pytest.fixture
def selectivity_world(db):
    """Target with WT + TM protocols on the same 3 compounds — tests AND across filters."""
    t = Target.objects.create(name="mEGFR")
    c1 = Compound.objects.create(target=t, smiles="CCO")
    c2 = Compound.objects.create(target=t, smiles="CCN")
    c3 = Compound.objects.create(target=t, smiles="CCC")

    p_wt = Protocol.objects.create(name="mEGFR TR-FRET WT")
    p_tm = Protocol.objects.create(name="mEGFR TR-FRET TM")
    assay_wt = Assay.objects.create(protocol=p_wt, target=t)
    assay_tm = Assay.objects.create(protocol=p_tm, target=t)

    # c1: potent on WT (5 nM), off on TM (5000 nM) — ideal selectivity profile
    # c2: potent on both (5 nM WT, 50 nM TM) — no selectivity window
    # c3: weak on WT (1000 nM), off on TM — fails WT filter
    _mk_result(assay_wt, c1, status="valid",
               results={"KPI": "IC50", "IC50": 5.0, "kpi_unit": "nM"}, row=0)
    _mk_result(assay_tm, c1, status="valid",
               results={"KPI": "IC50", "IC50": 5000.0, "kpi_unit": "nM"}, row=0)
    _mk_result(assay_wt, c2, status="valid",
               results={"KPI": "IC50", "IC50": 5.0, "kpi_unit": "nM"}, row=1)
    _mk_result(assay_tm, c2, status="valid",
               results={"KPI": "IC50", "IC50": 50.0, "kpi_unit": "nM"}, row=1)
    _mk_result(assay_wt, c3, status="valid",
               results={"KPI": "IC50", "IC50": 1000.0, "kpi_unit": "nM"}, row=2)
    _mk_result(assay_tm, c3, status="valid",
               results={"KPI": "IC50", "IC50": 5000.0, "kpi_unit": "nM"}, row=2)

    return {"t": t, "c1": c1, "c2": c2, "c3": c3, "p_wt": p_wt, "p_tm": p_tm}


# ---------------------------------------------------------------------------
# Spec-level validation (no DB state needed beyond the db marker)
# ---------------------------------------------------------------------------


def test_missing_target_returns_scope_error(db):
    res = execute(CompoundSelector())
    assert isinstance(res, ScopeError)


def test_threshold_without_metric_returns_spec_error(db):
    res = execute(CompoundSelector(
        registration_target_as_typed="ARd",
        measurement_filters=[
            MeasurementFilter(
                protocol_hint="HTRF",
                metric=None,
                threshold=Threshold(op="<", value=10.0, unit="nM"),
            ),
        ],
    ))
    assert isinstance(res, SpecError)
    assert res.field == "measurement_filters[0].metric"


# ---------------------------------------------------------------------------
# Pass-through of resolver responses
# ---------------------------------------------------------------------------


def test_target_miss_passthrough(ar_world):
    res = execute(CompoundSelector(
        registration_target_as_typed="zzzzznotarealtarget",
    ))
    assert isinstance(res, TargetMiss)
    assert res.field == FIELD_REGISTRATION_TARGET


def test_target_clarify_passthrough(db):
    gene = Gene.objects.create(symbol="MYC")
    t1 = Target.objects.create(name="Myc-Aur")
    t2 = Target.objects.create(name="Myc regulation")
    t1.genes.add(gene)
    t2.genes.add(gene)
    res = execute(CompoundSelector(registration_target_as_typed="MYC"))
    assert isinstance(res, TargetClarify)


def test_protocol_miss_passthrough(ar_world):
    res = execute(CompoundSelector(
        registration_target_as_typed="AR degraders",
        measurement_filters=[
            MeasurementFilter(protocol_hint="AlphaLISA"),
        ],
    ))
    assert isinstance(res, ProtocolMiss)
    assert res.field == FIELD_PROTOCOL_HINT
    assert res.filter_index == 0


def test_protocol_clarify_tags_filter_index(db):
    """When a later filter has an ambiguous protocol, the clarify result
    carries the right filter index so the UI pins the correct filter."""
    ar = Target.objects.create(name="AR degraders")
    c = Compound.objects.create(target=ar, smiles="CCO")
    for name in ("AR binding HTRF", "AR degradation HTRF"):
        p = Protocol.objects.create(name=name)
        a = Assay.objects.create(protocol=p, target=ar)
        DataSeries.objects.create(
            assay=a, compound=c, row=0, start_column=0, end_column=10,
        )
    # First filter has an unambiguous protocol; second has the ambiguous HTRF.
    # (Not a realistic prompt but tests index tagging.)
    solo = Protocol.objects.create(name="SolitaryBind")
    solo_a = Assay.objects.create(protocol=solo, target=ar)
    DataSeries.objects.create(
        assay=solo_a, compound=c, row=0, start_column=0, end_column=10,
    )
    res = execute(CompoundSelector(
        registration_target_as_typed="AR degraders",
        measurement_filters=[
            MeasurementFilter(protocol_hint="SolitaryBind"),
            MeasurementFilter(protocol_hint="HTRF"),
        ],
    ))
    assert isinstance(res, ProtocolClarify)
    assert res.filter_index == 1


# ---------------------------------------------------------------------------
# Happy path — single filter
# ---------------------------------------------------------------------------


def test_single_filter_selects_passing_compounds(ar_world):
    res = execute(CompoundSelector(
        registration_target_as_typed="AR degraders",
        measurement_filters=[
            MeasurementFilter(
                protocol_hint="HTRF",
                metric="IC50",
                threshold=Threshold(op="<", value=100.0, unit="nM"),
            ),
        ],
    ))
    assert isinstance(res, CompoundSelection)
    # c1 (5 nM) and c3 (50 nM) pass; c2 (200 nM) does not.
    expected = {ar_world["c1"].formatted_id, ar_world["c3"].formatted_id}
    assert set(res.compound_formatted_ids) == expected
    assert res.n_matched == 2
    assert res.target_names == ["AR degraders"]
    assert res.protocol_names == ["AR binding HTRF"]
    assert "AR binding HTRF" in res.scope_sentence
    assert "IC50 < 100.0 nM" in res.scope_sentence


def test_no_filters_returns_all_compounds_in_target(ar_world):
    res = execute(CompoundSelector(registration_target_as_typed="AR degraders"))
    assert isinstance(res, CompoundSelection)
    assert res.n_matched == 3
    assert res.protocol_names == []  # no filters → no protocols in URL
    # Order must be stable via reg_number ascending.
    assert res.compound_formatted_ids == sorted(res.compound_formatted_ids)


def test_zero_matches_returns_empty_selection(ar_world):
    res = execute(CompoundSelector(
        registration_target_as_typed="AR degraders",
        measurement_filters=[
            MeasurementFilter(
                protocol_hint="HTRF",
                metric="IC50",
                threshold=Threshold(op="<", value=0.001, unit="nM"),
            ),
        ],
    ))
    assert isinstance(res, CompoundSelection)
    assert res.n_matched == 0
    assert res.compound_formatted_ids == []
    # Still produces a scope sentence and protocol name for the UI.
    assert res.protocol_names == ["AR binding HTRF"]


# ---------------------------------------------------------------------------
# Cross-protocol selectivity — the pivot's headline feature
# ---------------------------------------------------------------------------


def test_cross_protocol_selectivity_intersects(selectivity_world):
    """c1 passes both filters (potent WT + inactive TM); c2 fails TM filter; c3 fails WT."""
    res = execute(CompoundSelector(
        registration_target_as_typed="mEGFR",
        measurement_filters=[
            MeasurementFilter(
                protocol_hint="WT",
                metric="IC50",
                threshold=Threshold(op="<", value=100.0, unit="nM"),
            ),
            MeasurementFilter(
                protocol_hint="TM",
                metric="IC50",
                threshold=Threshold(op=">", value=1000.0, unit="nM"),
            ),
        ],
    ))
    assert isinstance(res, CompoundSelection)
    assert res.compound_formatted_ids == [selectivity_world["c1"].formatted_id]
    assert res.n_matched == 1
    # Both protocols surface in the URL list.
    assert set(res.protocol_names) == {
        "mEGFR TR-FRET WT", "mEGFR TR-FRET TM",
    }
    # Scope sentence names both filters joined by AND.
    assert " AND " in res.scope_sentence
    assert "WT" in res.scope_sentence and "TM" in res.scope_sentence


def test_selectivity_empty_set_when_no_compound_passes_all(selectivity_world):
    res = execute(CompoundSelector(
        registration_target_as_typed="mEGFR",
        measurement_filters=[
            # Impossible for any compound: WT < 0.1 nM AND TM > 0.01 nM
            MeasurementFilter(
                protocol_hint="WT", metric="IC50",
                threshold=Threshold(op="<", value=0.1, unit="nM"),
            ),
            MeasurementFilter(
                protocol_hint="TM", metric="IC50",
                threshold=Threshold(op=">", value=0.01, unit="nM"),
            ),
        ],
    ))
    assert isinstance(res, CompoundSelection)
    assert res.n_matched == 0


# ---------------------------------------------------------------------------
# Filter with no threshold — "tested in X" / "has ANY measurement"
# ---------------------------------------------------------------------------


def test_filter_without_threshold_selects_any_valid_measurement(ar_world):
    res = execute(CompoundSelector(
        registration_target_as_typed="AR degraders",
        measurement_filters=[
            MeasurementFilter(protocol_hint="HTRF", metric="IC50"),
        ],
    ))
    assert isinstance(res, CompoundSelection)
    # All three compounds have a valid HTRF IC50, so all three survive.
    assert res.n_matched == 3


def test_filter_without_protocol_matches_any_protocol_in_scope(selectivity_world):
    """No protocol_hint → compound needs any valid IC50 in the scope;
    each of c1/c2/c3 has both WT and TM measurements, so all pass."""
    res = execute(CompoundSelector(
        registration_target_as_typed="mEGFR",
        measurement_filters=[
            MeasurementFilter(metric="IC50"),
        ],
    ))
    assert isinstance(res, CompoundSelection)
    assert res.n_matched == 3


# ---------------------------------------------------------------------------
# Scope sentence across scope_kinds
# ---------------------------------------------------------------------------


def test_scope_sentence_reg_only_no_filters(ar_world):
    res = execute(CompoundSelector(registration_target_as_typed="AR degraders"))
    assert isinstance(res, CompoundSelection)
    assert res.scope_sentence == "Showing compounds registered to AR degraders"


def test_scope_sentence_both_same_with_filter(ar_world):
    res = execute(CompoundSelector(
        registration_target_as_typed="AR degraders",
        assay_target_as_typed="AR degraders",
        measurement_filters=[
            MeasurementFilter(
                protocol_hint="HTRF", metric="IC50",
                threshold=Threshold(op="<", value=100.0, unit="nM"),
            ),
        ],
    ))
    assert isinstance(res, CompoundSelection)
    assert res.scope_sentence == (
        "Showing compounds registered to AR degraders and tested against "
        "AR degraders where AR binding HTRF IC50 < 100.0 nM"
    )


def test_scope_sentence_assay_only(db):
    ar = Target.objects.create(name="AR degraders")
    akt = Target.objects.create(name="AKT")
    # A compound registered to ARd tested in an AKT protocol.
    c = Compound.objects.create(target=ar, smiles="CCO")
    p = Protocol.objects.create(name="AKT HTRF")
    a = Assay.objects.create(protocol=p, target=akt)
    _mk_result(a, c, status="valid",
               results={"KPI": "IC50", "IC50": 50.0, "kpi_unit": "nM"})

    res = execute(CompoundSelector(assay_target_as_typed="AKT"))
    assert isinstance(res, CompoundSelection)
    assert "tested against AKT" in res.scope_sentence
