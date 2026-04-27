"""Executor tests (post-pivot).

Composes resolvers + evaluator end-to-end and verifies the compound-
intersection semantics introduced by the pivot: filters AND together,
each filter contributes its protocol name to the redirect URL, scope
sentence reflects all filters.
"""

from __future__ import annotations

import pytest

from datetime import datetime, timezone
from django.contrib.auth import get_user_model
from compounds.assays.models import AnalysisResult, Assay, DataSeries, Protocol
from compounds.nlp.executor import execute
from compounds.nlp.spec import (
    CompoundMiss,
    CompoundSelector,
    CompoundSelection,
    DateRange,
    FIELD_ASSAYED_BY,
    FIELD_COMPOUND_REF,
    FIELD_METRIC,
    FIELD_PROTOCOL_HINT,
    FIELD_REGISTERED_BY,
    FIELD_REGISTRATION_TARGET,
    FIELD_SCAFFOLD_HINT,
    MeasurementFilter,
    MetricMiss,
    ProtocolClarify,
    ProtocolMiss,
    ScaffoldClarify,
    ScaffoldMiss,
    ScopeError,
    SpecError,
    TargetClarify,
    TargetMiss,
    Threshold,
    UserClarify,
    UserMiss,
)
from compounds.registry.models import Compound, Gene, Target

User = get_user_model()


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


def test_fully_empty_selector_returns_scope_error(db):
    """A selector with no target AND no other narrowing predicate would
    trivially list every compound — guard with a helpful ScopeError."""
    res = execute(CompoundSelector())
    assert isinstance(res, ScopeError)
    assert "narrow by" in res.message.lower() or "narrow" in res.message.lower()


def test_target_less_selector_with_scaffold_runs(db):
    t = Target.objects.create(name="AR degraders")
    c_pyr = Compound.objects.create(target=t, smiles="c1ccncc1")
    Compound.objects.create(target=t, smiles="CCO")
    res = execute(CompoundSelector(scaffold_hints=["pyridine"]))
    assert isinstance(res, CompoundSelection)
    assert res.compound_formatted_ids == [c_pyr.formatted_id]
    assert res.target_names == []


def test_target_less_selector_with_filter_runs(ar_world):
    """Cross-target HTRF IC50 threshold with no project restriction —
    all AR compounds under 100 nM survive."""
    res = execute(CompoundSelector(
        measurement_filters=[
            MeasurementFilter(
                protocol_hint="HTRF", metric="IC50",
                threshold=Threshold(op="<", value=100.0, unit="nM"),
            ),
        ],
    ))
    assert isinstance(res, CompoundSelection)
    assert res.n_matched == 2  # c1 (5 nM) and c3 (50 nM)


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


# ---------------------------------------------------------------------------
# Date filters — registered_date_range on CompoundSelector, assay_date_range
# on MeasurementFilter
# ---------------------------------------------------------------------------


@pytest.fixture
def dated_world(db):
    """Compounds registered across three calendar years + HTRF measurements
    across two quarters. Lets us exercise both date filters independently."""
    t = Target.objects.create(name="CDK4")
    # Compounds — override registered_at post-save (auto_now_add means we
    # can't set it on create(); update() bypasses that).
    c_2024 = Compound.objects.create(target=t, smiles="CCO")
    c_2025 = Compound.objects.create(target=t, smiles="CCN")
    c_2026 = Compound.objects.create(target=t, smiles="CCC")
    Compound.objects.filter(pk=c_2024.pk).update(
        registered_at=datetime(2024, 6, 15, tzinfo=timezone.utc),
    )
    Compound.objects.filter(pk=c_2025.pk).update(
        registered_at=datetime(2025, 8, 20, tzinfo=timezone.utc),
    )
    Compound.objects.filter(pk=c_2026.pk).update(
        registered_at=datetime(2026, 2, 10, tzinfo=timezone.utc),
    )
    # Refresh so in-memory copies are current (some tests read them).
    c_2024.refresh_from_db()
    c_2025.refresh_from_db()
    c_2026.refresh_from_db()

    p = Protocol.objects.create(name="CDK4 HTRF")
    # Two assay runs on different dates, each with measurements for some
    # compounds. c_2024 measured in Q1 2026, c_2025 in Q2 2026.
    assay_q1 = Assay.objects.create(protocol=p, target=t)
    Assay.objects.filter(pk=assay_q1.pk).update(
        created_at=datetime(2026, 2, 15, tzinfo=timezone.utc),
    )
    assay_q2 = Assay.objects.create(protocol=p, target=t)
    Assay.objects.filter(pk=assay_q2.pk).update(
        created_at=datetime(2026, 5, 20, tzinfo=timezone.utc),
    )
    _mk_result(assay_q1, c_2024, status="valid",
               results={"KPI": "IC50", "IC50": 5.0, "kpi_unit": "nM"})
    _mk_result(assay_q1, c_2025, status="valid",
               results={"KPI": "IC50", "IC50": 10.0, "kpi_unit": "nM"})
    _mk_result(assay_q2, c_2026, status="valid",
               results={"KPI": "IC50", "IC50": 15.0, "kpi_unit": "nM"})

    return {
        "t": t, "c_2024": c_2024, "c_2025": c_2025, "c_2026": c_2026,
        "assay_q1": assay_q1, "assay_q2": assay_q2,
    }


def test_registered_date_range_narrows_base_scope(dated_world):
    res = execute(CompoundSelector(
        registration_target_as_typed="CDK4",
        registered_date_range=DateRange(after="2025-01-01", before="2026-01-01"),
    ))
    assert isinstance(res, CompoundSelection)
    assert res.compound_formatted_ids == [dated_world["c_2025"].formatted_id]


def test_registered_date_after_only(dated_world):
    res = execute(CompoundSelector(
        registration_target_as_typed="CDK4",
        registered_date_range=DateRange(after="2025-01-01", before=None),
    ))
    assert isinstance(res, CompoundSelection)
    assert set(res.compound_formatted_ids) == {
        dated_world["c_2025"].formatted_id,
        dated_world["c_2026"].formatted_id,
    }


def test_registered_date_before_only(dated_world):
    res = execute(CompoundSelector(
        registration_target_as_typed="CDK4",
        registered_date_range=DateRange(after=None, before="2025-01-01"),
    ))
    assert isinstance(res, CompoundSelection)
    assert res.compound_formatted_ids == [dated_world["c_2024"].formatted_id]


def test_assay_date_range_on_filter(dated_world):
    """Only compounds whose HTRF measurement falls in Q1 2026 survive."""
    res = execute(CompoundSelector(
        registration_target_as_typed="CDK4",
        measurement_filters=[
            MeasurementFilter(
                protocol_hint="HTRF",
                metric="IC50",
                threshold=Threshold(op="<", value=100, unit="nM"),
                assay_date_range=DateRange(after="2026-01-01", before="2026-04-01"),
            ),
        ],
    ))
    assert isinstance(res, CompoundSelection)
    # c_2024 and c_2025 were both measured in Q1 2026. c_2026's
    # measurement was in Q2 so it's excluded despite passing the threshold.
    assert set(res.compound_formatted_ids) == {
        dated_world["c_2024"].formatted_id,
        dated_world["c_2025"].formatted_id,
    }


def test_combined_registration_and_assay_dates(dated_world):
    """Registration in 2025 (narrows to c_2025) AND measurement in Q1 2026
    (c_2025's measurement IS in Q1 2026) → exactly c_2025."""
    res = execute(CompoundSelector(
        registration_target_as_typed="CDK4",
        registered_date_range=DateRange(after="2025-01-01", before="2026-01-01"),
        measurement_filters=[
            MeasurementFilter(
                protocol_hint="HTRF",
                metric="IC50",
                threshold=Threshold(op="<", value=100, unit="nM"),
                assay_date_range=DateRange(after="2026-01-01", before="2026-04-01"),
            ),
        ],
    ))
    assert isinstance(res, CompoundSelection)
    assert res.compound_formatted_ids == [dated_world["c_2025"].formatted_id]


def test_scope_sentence_includes_date_phrases(dated_world):
    res = execute(CompoundSelector(
        registration_target_as_typed="CDK4",
        registered_date_range=DateRange(after="2025-01-01", before="2026-01-01"),
        measurement_filters=[
            MeasurementFilter(
                protocol_hint="HTRF",
                metric="IC50",
                threshold=Threshold(op="<", value=100, unit="nM"),
                assay_date_range=DateRange(after="2026-01-01", before="2026-04-01"),
            ),
        ],
    ))
    assert isinstance(res, CompoundSelection)
    assert "registered between 2025-01-01 and 2026-01-01" in res.scope_sentence
    assert "measured between 2026-01-01 and 2026-04-01" in res.scope_sentence


# ---------------------------------------------------------------------------
# User filters — registered_by on CompoundSelector, assayed_by on
# MeasurementFilter
# ---------------------------------------------------------------------------


@pytest.fixture
def user_world(db):
    """Two chemists (Alice, Bob) each registering and assaying compounds.
    Enables us to exercise both filters and their clarify paths."""
    alice = User.objects.create_user(
        username="alice.jones", email="alice.jones@ncl.ac.uk",
        first_name="Alice", last_name="Jones",
    )
    bob = User.objects.create_user(
        username="bob.smith", email="bob.smith@ncl.ac.uk",
        first_name="Bob", last_name="Smith",
    )
    # Second "Alice" — different last name. Clarify trigger for "Alice".
    alice2 = User.objects.create_user(
        username="alice.brown", email="alice.brown@ncl.ac.uk",
        first_name="Alice", last_name="Brown",
    )

    t = Target.objects.create(name="CDK4")
    c_alice = Compound.objects.create(target=t, smiles="CCO", registered_by=alice)
    c_bob = Compound.objects.create(target=t, smiles="CCN", registered_by=bob)
    c_alice2 = Compound.objects.create(target=t, smiles="CCC", registered_by=alice2)

    protocol = Protocol.objects.create(name="CDK4 HTRF")
    assay_by_alice = Assay.objects.create(protocol=protocol, target=t, created_by=alice)
    assay_by_bob = Assay.objects.create(protocol=protocol, target=t, created_by=bob)
    # alice2 also runs one assay — needed so "Alice" on the assayed_by
    # path is ambiguous across the assay_creators_qs scope.
    assay_by_alice2 = Assay.objects.create(protocol=protocol, target=t, created_by=alice2)
    # Alice's compound measured by Bob (cross-hatched to test we don't
    # conflate registrant with assayer).
    _mk_result(assay_by_bob, c_alice,
               status="valid", results={"KPI": "IC50", "IC50": 5.0, "kpi_unit": "nM"})
    # Bob's compound measured by Alice.
    _mk_result(assay_by_alice, c_bob,
               status="valid", results={"KPI": "IC50", "IC50": 10.0, "kpi_unit": "nM"})
    # alice2's compound measured by herself.
    _mk_result(assay_by_alice2, c_alice2,
               status="valid", results={"KPI": "IC50", "IC50": 15.0, "kpi_unit": "nM"})

    return {
        "alice": alice, "bob": bob, "alice2": alice2,
        "c_alice": c_alice, "c_bob": c_bob, "c_alice2": c_alice2,
        "protocol": protocol,
    }


def test_registered_by_unique_full_name_resolves(user_world):
    res = execute(CompoundSelector(
        registration_target_as_typed="CDK4",
        registered_by_as_typed="Alice Jones",
    ))
    assert isinstance(res, CompoundSelection)
    assert res.compound_formatted_ids == [user_world["c_alice"].formatted_id]


def test_registered_by_unique_username_resolves(user_world):
    res = execute(CompoundSelector(
        registration_target_as_typed="CDK4",
        registered_by_as_typed="bob.smith",
    ))
    assert isinstance(res, CompoundSelection)
    assert res.compound_formatted_ids == [user_world["c_bob"].formatted_id]


def test_registered_by_unique_email_resolves(user_world):
    res = execute(CompoundSelector(
        registration_target_as_typed="CDK4",
        registered_by_as_typed="bob.smith@ncl.ac.uk",
    ))
    assert isinstance(res, CompoundSelection)
    assert res.compound_formatted_ids == [user_world["c_bob"].formatted_id]


def test_registered_by_ambiguous_first_name_clarifies(user_world):
    res = execute(CompoundSelector(
        registration_target_as_typed="CDK4",
        registered_by_as_typed="Alice",
    ))
    assert isinstance(res, UserClarify)
    assert res.field == FIELD_REGISTERED_BY
    displays = {c.display for c in res.candidates}
    assert displays == {"Alice Jones", "Alice Brown"}


def test_registered_by_miss_surfaces_suggestions(user_world):
    res = execute(CompoundSelector(
        registration_target_as_typed="CDK4",
        registered_by_as_typed="Zachary",
    ))
    assert isinstance(res, UserMiss)
    assert res.field == FIELD_REGISTERED_BY
    assert len(res.suggestions) > 0


def test_registered_by_pinned_id_bypasses_resolver(user_world):
    res = execute(CompoundSelector(
        registration_target_as_typed="CDK4",
        registered_by_as_typed="Alice",
        registered_by_id=str(user_world["alice"].pk),
    ))
    assert isinstance(res, CompoundSelection)
    assert res.compound_formatted_ids == [user_world["c_alice"].formatted_id]


def test_assayed_by_filters_per_measurement(user_world):
    """c_alice was measured BY Bob; filtering assayed_by=Alice Jones
    should select only compounds Alice herself assayed."""
    res = execute(CompoundSelector(
        registration_target_as_typed="CDK4",
        measurement_filters=[
            MeasurementFilter(
                protocol_hint="HTRF",
                metric="IC50",
                assayed_by_as_typed="Alice Jones",
            ),
        ],
    ))
    assert isinstance(res, CompoundSelection)
    # Alice measured only c_bob — so c_bob survives.
    assert res.compound_formatted_ids == [user_world["c_bob"].formatted_id]


def test_assayed_by_ambiguous_clarifies_with_filter_index(user_world):
    res = execute(CompoundSelector(
        registration_target_as_typed="CDK4",
        measurement_filters=[
            MeasurementFilter(
                protocol_hint="HTRF",
                metric="IC50",
                assayed_by_as_typed="Alice",
            ),
        ],
    ))
    assert isinstance(res, UserClarify)
    assert res.field == FIELD_ASSAYED_BY
    assert res.filter_index == 0


def test_registered_and_assayed_by_combined(user_world):
    """Registered by Alice Jones AND assayed by Bob Smith → only
    c_alice (registered by Alice, measured by Bob)."""
    res = execute(CompoundSelector(
        registration_target_as_typed="CDK4",
        registered_by_as_typed="Alice Jones",
        measurement_filters=[
            MeasurementFilter(
                protocol_hint="HTRF",
                metric="IC50",
                assayed_by_as_typed="Bob Smith",
            ),
        ],
    ))
    assert isinstance(res, CompoundSelection)
    assert res.compound_formatted_ids == [user_world["c_alice"].formatted_id]


# ---------------------------------------------------------------------------
# Assay-selection path (slice 12)
# ---------------------------------------------------------------------------


def test_assay_query_target_and_date_filters(dated_world):
    from compounds.nlp.executor import execute_assay_query
    from compounds.nlp.spec import AssaySelector, AssaySelection

    res = execute_assay_query(AssaySelector(
        target_as_typed="CDK4",
        date_range=DateRange(after="2026-01-01", before="2026-04-01"),
    ))
    assert isinstance(res, AssaySelection)
    # dated_world fixture has assay_q1 (Feb) and assay_q2 (May). Only q1
    # falls in the Q1 2026 window.
    assert res.n_matched == 1
    assert res.assay_ids == [str(dated_world["assay_q1"].pk)]
    assert res.target_names == ["CDK4"]
    assert "carried out between 2026-01-01 and 2026-04-01" in res.scope_sentence


def test_assay_query_by_creator(user_world):
    from compounds.nlp.executor import execute_assay_query
    from compounds.nlp.spec import AssaySelector, AssaySelection

    res = execute_assay_query(AssaySelector(
        target_as_typed="CDK4",
        created_by_as_typed="Alice Jones",
    ))
    assert isinstance(res, AssaySelection)
    # Only one assay created by Alice Jones in user_world.
    assert res.n_matched == 1
    assert "conducted by Alice Jones" in res.scope_sentence


def test_assay_query_nickname_prefix_matches(user_world):
    """'Jess' prefix of 'Jessica' — existing partial-tier resolution
    should handle it. Here the fixture has Alice/Bob so we test with 'Al'
    (prefix of Alice) to exercise the prefix-match path against an
    ambiguous-first-name clarify."""
    from compounds.nlp.executor import execute_assay_query
    from compounds.nlp.spec import AssaySelector

    res = execute_assay_query(AssaySelector(
        target_as_typed="CDK4",
        created_by_as_typed="Al",
    ))
    # 'Al' matches Alice Jones (partial tier) AND Alice Brown (same).
    # Both are in assay_creators_qs since both created assays in
    # user_world. So: UserClarify.
    assert isinstance(res, UserClarify)
    displays = {c.display for c in res.candidates}
    assert "Alice Jones" in displays


def test_assay_query_target_clarify_passthrough(db):
    """Ambiguous target on an assay query surfaces via TargetClarify,
    tagged with FIELD_ASSAY_TARGET so the UI knows which field to pin."""
    from compounds.nlp.executor import execute_assay_query
    from compounds.nlp.spec import AssaySelector

    gene = Gene.objects.create(symbol="MYC")
    t1 = Target.objects.create(name="Myc-Aur")
    t2 = Target.objects.create(name="Myc regulation")
    t1.genes.add(gene)
    t2.genes.add(gene)

    res = execute_assay_query(AssaySelector(target_as_typed="MYC"))
    assert isinstance(res, TargetClarify)
    assert res.field == "assay_target_as_typed"


def test_scope_sentence_includes_user_phrases(user_world):
    res = execute(CompoundSelector(
        registration_target_as_typed="CDK4",
        registered_by_as_typed="Alice Jones",
        measurement_filters=[
            MeasurementFilter(
                protocol_hint="HTRF",
                metric="IC50",
                assayed_by_as_typed="Bob Smith",
            ),
        ],
    ))
    assert isinstance(res, CompoundSelection)
    assert "registered by Alice Jones" in res.scope_sentence
    assert "assayed by Bob Smith" in res.scope_sentence


# ---------------------------------------------------------------------------
# Scaffold / substructure filters (slice 13)
# ---------------------------------------------------------------------------


@pytest.fixture
def scaffold_world(db):
    """Compounds spanning a pyridine, a pyrimidine, an amide, and a plain
    aliphatic. Backend-side RDKit HasSubstructMatch is the authority so we
    use real (aromatic) SMILES and rely on the catalog SMARTS.
    """
    t = Target.objects.create(name="CDK4")
    c_pyridine = Compound.objects.create(target=t, smiles="c1ccncc1")     # pyridine
    c_pyrimidine = Compound.objects.create(target=t, smiles="c1cncnc1")   # pyrimidine
    c_amide = Compound.objects.create(target=t, smiles="CC(=O)N")         # acetamide
    c_plain = Compound.objects.create(target=t, smiles="CCO")             # ethanol

    return {
        "t": t,
        "c_pyridine": c_pyridine,
        "c_pyrimidine": c_pyrimidine,
        "c_amide": c_amide,
        "c_plain": c_plain,
    }


def test_scaffold_hint_filters_by_substructure(scaffold_world):
    res = execute(CompoundSelector(
        registration_target_as_typed="CDK4",
        scaffold_hints=["pyridine"],
    ))
    assert isinstance(res, CompoundSelection)
    assert res.compound_formatted_ids == [
        scaffold_world["c_pyridine"].formatted_id,
    ]


def test_scaffold_hint_unresolved_clarifies(scaffold_world):
    # "amid" hits both amide and sulfonamide → clarify.
    res = execute(CompoundSelector(
        registration_target_as_typed="CDK4",
        scaffold_hints=["amid"],
    ))
    assert isinstance(res, ScaffoldClarify)
    assert res.scaffold_index == 0
    assert res.field == FIELD_SCAFFOLD_HINT
    names = {c.name for c in res.candidates}
    assert "amide" in names and "sulfonamide" in names


def test_scaffold_hint_miss_surfaces_suggestions(scaffold_world):
    res = execute(CompoundSelector(
        registration_target_as_typed="CDK4",
        scaffold_hints=["zzzznothing"],
    ))
    assert isinstance(res, ScaffoldMiss)
    assert res.field == FIELD_SCAFFOLD_HINT
    assert len(res.suggestions) > 0


def test_scaffold_hint_pinned_id_bypasses_resolver(scaffold_world):
    # Ambiguous hint "amid", but pinned to "amide" — resolves directly.
    res = execute(CompoundSelector(
        registration_target_as_typed="CDK4",
        scaffold_hints=["amid"],
        scaffold_ids=["amide"],
    ))
    assert isinstance(res, CompoundSelection)
    assert res.compound_formatted_ids == [
        scaffold_world["c_amide"].formatted_id,
    ]


def test_multiple_scaffolds_and_intersect(scaffold_world):
    """Require BOTH pyridine AND amide — none of the fixture compounds
    contain both, so the selection is empty. Exercises the AND semantics."""
    res = execute(CompoundSelector(
        registration_target_as_typed="CDK4",
        scaffold_hints=["pyridine", "amide"],
    ))
    assert isinstance(res, CompoundSelection)
    assert res.n_matched == 0


def test_scope_sentence_includes_scaffold_phrase(scaffold_world):
    res = execute(CompoundSelector(
        registration_target_as_typed="CDK4",
        scaffold_hints=["pyridine"],
    ))
    assert isinstance(res, CompoundSelection)
    assert "containing pyridine" in res.scope_sentence


def test_scope_sentence_multiple_scaffolds_joined_with_and(scaffold_world):
    res = execute(CompoundSelector(
        registration_target_as_typed="CDK4",
        scaffold_hints=["pyridine", "amide"],
    ))
    assert isinstance(res, CompoundSelection)
    assert "containing pyridine AND amide" in res.scope_sentence


def test_scaffold_filter_combines_with_measurement_filter(db):
    """Scaffold narrowing composes with measurement filters — compounds
    must contain the substructure AND pass the measurement threshold."""
    t = Target.objects.create(name="CDK4")
    c_pyr = Compound.objects.create(target=t, smiles="c1ccncc1")
    c_plain = Compound.objects.create(target=t, smiles="CCO")
    p = Protocol.objects.create(name="CDK4 HTRF")
    a = Assay.objects.create(protocol=p, target=t)
    _mk_result(a, c_pyr, status="valid",
               results={"KPI": "IC50", "IC50": 5.0, "kpi_unit": "nM"}, row=0)
    _mk_result(a, c_plain, status="valid",
               results={"KPI": "IC50", "IC50": 5.0, "kpi_unit": "nM"}, row=1)

    res = execute(CompoundSelector(
        registration_target_as_typed="CDK4",
        scaffold_hints=["pyridine"],
        measurement_filters=[
            MeasurementFilter(
                protocol_hint="HTRF", metric="IC50",
                threshold=Threshold(op="<", value=100.0, unit="nM"),
            ),
        ],
    ))
    assert isinstance(res, CompoundSelection)
    # Both pass the threshold, only c_pyr contains pyridine.
    assert res.compound_formatted_ids == [c_pyr.formatted_id]


def test_assay_query_with_scaffold(db):
    """Assay selector with a scaffold hint narrows to assays whose compounds
    contain the substructure."""
    from compounds.nlp.executor import execute_assay_query
    from compounds.nlp.spec import AssaySelector, AssaySelection

    t = Target.objects.create(name="CDK4")
    c_pyr = Compound.objects.create(target=t, smiles="c1ccncc1")
    c_plain = Compound.objects.create(target=t, smiles="CCO")
    p = Protocol.objects.create(name="CDK4 HTRF")
    a_pyr = Assay.objects.create(protocol=p, target=t)
    a_plain = Assay.objects.create(protocol=p, target=t)
    DataSeries.objects.create(assay=a_pyr, compound=c_pyr,
                              row=0, start_column=0, end_column=10)
    DataSeries.objects.create(assay=a_plain, compound=c_plain,
                              row=0, start_column=0, end_column=10)

    res = execute_assay_query(AssaySelector(
        target_as_typed="CDK4",
        scaffold_hints=["pyridine"],
    ))
    assert isinstance(res, AssaySelection)
    # Only a_pyr ran on a pyridine-containing compound.
    assert res.assay_ids == [str(a_pyr.pk)]
    assert "compounds containing pyridine" in res.scope_sentence


# ---------------------------------------------------------------------------
# Compound-ID pinning (slice 14) — UNION semantics, not narrowing
# ---------------------------------------------------------------------------


@pytest.fixture
def pin_world(db):
    """Two unrelated programmes with one compound each, plus an HTRF assay
    on the ARd compound. Lets us test that pinning a CDK4 compound makes
    it visible alongside a query scoped to ARd."""
    ar = Target.objects.create(name="ARd")
    cdk4 = Target.objects.create(name="CDK4")
    c_ar = Compound.objects.create(target=ar, smiles="CCO", reg_number=10001)
    c_cdk4 = Compound.objects.create(target=cdk4, smiles="CCN", reg_number=20002)
    p = Protocol.objects.create(name="ARd binding HTRF")
    a = Assay.objects.create(protocol=p, target=ar)
    _mk_result(a, c_ar, status="valid",
               results={"KPI": "IC50", "IC50": 5.0, "kpi_unit": "nM"})
    return {"ar": ar, "cdk4": cdk4, "c_ar": c_ar, "c_cdk4": c_cdk4}


def test_pinned_only_runs_with_no_other_predicate(pin_world):
    res = execute(CompoundSelector(
        compound_refs_as_typed=["NCL-00010001"],
    ))
    assert isinstance(res, CompoundSelection)
    assert res.compound_formatted_ids == [pin_world["c_ar"].formatted_id]
    assert "plus NCL-00010001 (pinned)" in res.scope_sentence


def test_pinned_compound_appears_alongside_target_query(pin_world):
    """Pinning a CDK4 compound while scoped to ARd surfaces both: ARd's
    compound (from the scope) plus the CDK4 compound (from the pin)."""
    res = execute(CompoundSelector(
        registration_target_as_typed="ARd",
        compound_refs_as_typed=["NCL-00020002"],  # the CDK4 compound
    ))
    assert isinstance(res, CompoundSelection)
    formatted = set(res.compound_formatted_ids)
    assert pin_world["c_ar"].formatted_id in formatted
    assert pin_world["c_cdk4"].formatted_id in formatted


def test_pinned_compound_survives_strict_filter(pin_world):
    """Even when the measurement filter excludes the pinned compound (it
    has no HTRF data), it still appears in the final selection."""
    res = execute(CompoundSelector(
        registration_target_as_typed="ARd",
        measurement_filters=[
            MeasurementFilter(
                protocol_hint="HTRF", metric="IC50",
                threshold=Threshold(op="<", value=10.0, unit="nM"),
            ),
        ],
        compound_refs_as_typed=["NCL-00020002"],
    ))
    assert isinstance(res, CompoundSelection)
    formatted = set(res.compound_formatted_ids)
    # c_ar passes the filter; c_cdk4 has no HTRF data but is pinned.
    assert pin_world["c_ar"].formatted_id in formatted
    assert pin_world["c_cdk4"].formatted_id in formatted


def test_pinned_bare_number_resolves(pin_world):
    res = execute(CompoundSelector(compound_refs_as_typed=["10001"]))
    assert isinstance(res, CompoundSelection)
    assert res.compound_formatted_ids == [pin_world["c_ar"].formatted_id]


def test_pinned_unknown_id_returns_compound_miss(pin_world):
    res = execute(CompoundSelector(compound_refs_as_typed=["NCL-99999999"]))
    assert isinstance(res, CompoundMiss)
    assert res.field == FIELD_COMPOUND_REF


def test_pinned_unparseable_ref_returns_compound_miss(pin_world):
    res = execute(CompoundSelector(compound_refs_as_typed=["the lead compound"]))
    assert isinstance(res, CompoundMiss)


def test_first_unresolvable_pin_short_circuits(pin_world):
    """Mixed list: one valid pin and one bogus. The miss surfaces with
    its index so the user sees which entry went wrong."""
    res = execute(CompoundSelector(
        compound_refs_as_typed=["NCL-00010001", "garbage"],
    ))
    assert isinstance(res, CompoundMiss)
    assert res.ref_index == 1


def test_multiple_pins_all_appear(pin_world):
    res = execute(CompoundSelector(
        compound_refs_as_typed=["NCL-00010001", "NCL-00020002"],
    ))
    assert isinstance(res, CompoundSelection)
    assert set(res.compound_formatted_ids) == {
        pin_world["c_ar"].formatted_id,
        pin_world["c_cdk4"].formatted_id,
    }


def test_scope_sentence_lists_multiple_pins(pin_world):
    res = execute(CompoundSelector(
        compound_refs_as_typed=["NCL-00010001", "NCL-00020002"],
    ))
    assert isinstance(res, CompoundSelection)
    assert "plus NCL-00010001, NCL-00020002 (pinned)" in res.scope_sentence


# ---------------------------------------------------------------------------
# Metric mismatch surfacing (slice 15)
# ---------------------------------------------------------------------------


@pytest.fixture
def kpi_world(db):
    """A protocol whose rows record a KPI other than the user's typed
    one — exposes the slice-15 MetricMiss path. ARd binding HTRF
    cellular records pIC50 / Hill (not IC50)."""
    t = Target.objects.create(name="ARd")
    c1 = Compound.objects.create(target=t, smiles="CCO")
    c2 = Compound.objects.create(target=t, smiles="CCN")
    p = Protocol.objects.create(name="ARd binding HTRF cellular")
    a = Assay.objects.create(protocol=p, target=t)
    _mk_result(a, c1, status="valid",
               results={"KPI": "pIC50", "pIC50": 7.5, "kpi_unit": "unitless"}, row=0)
    _mk_result(a, c2, status="valid",
               results={"KPI": "Hill", "Hill": 1.1, "kpi_unit": "unitless"}, row=1)
    return {"t": t, "c1": c1, "c2": c2, "protocol": p}


def test_metric_mismatch_returns_metric_miss_with_available_kpis(kpi_world):
    """User typed 'IC50' but rows record pIC50 / Hill. Pre-slice-15
    this silently returned 'Found 0'; now it returns a MetricMiss
    with the available KPIs so the response answers itself."""
    res = execute(CompoundSelector(
        registration_target_as_typed="ARd",
        measurement_filters=[
            MeasurementFilter(
                protocol_hint="HTRF",
                metric="IC50",
                threshold=Threshold(op="<", value=50.0, unit=None),
            ),
        ],
    ))
    assert isinstance(res, MetricMiss)
    assert res.field == FIELD_METRIC
    assert res.query == "IC50"
    assert res.filter_index == 0
    assert set(res.available_metrics) == {"pIC50", "Hill"}


def test_metric_mismatch_filter_index_tags_correctly(kpi_world):
    """Two filters; first has a valid KPI, second misses. The Miss
    carries filter_index=1 so the UI knows which filter went wrong."""
    res = execute(CompoundSelector(
        registration_target_as_typed="ARd",
        measurement_filters=[
            MeasurementFilter(protocol_hint="HTRF", metric="pIC50"),
            MeasurementFilter(protocol_hint="HTRF", metric="IC50"),
        ],
    ))
    assert isinstance(res, MetricMiss)
    assert res.filter_index == 1


def test_lenient_metric_match_finds_rows_with_punctuation(db):
    """User types 'IC50', rows are stored as 'ic-50'. The lenient
    match in the evaluator + the metric-existence check both treat
    them as the same KPI."""
    t = Target.objects.create(name="ARd")
    c = Compound.objects.create(target=t, smiles="CCO")
    p = Protocol.objects.create(name="ARd HTRF")
    a = Assay.objects.create(protocol=p, target=t)
    _mk_result(a, c, status="valid",
               results={"KPI": "ic-50", "ic-50": 5.0, "kpi_unit": "nM"})
    res = execute(CompoundSelector(
        registration_target_as_typed="ARd",
        measurement_filters=[
            MeasurementFilter(
                protocol_hint="HTRF", metric="IC50",
                threshold=Threshold(op="<", value=10.0, unit="nM"),
            ),
        ],
    ))
    assert isinstance(res, CompoundSelection)
    assert res.compound_formatted_ids == [c.formatted_id]


def test_metric_miss_when_no_rows_in_scope(db):
    """No rows in scope at all → available_metrics is empty.
    Frontend renders 'no measurements in scope under any KPI'."""
    Target.objects.create(name="ARd")
    res = execute(CompoundSelector(
        registration_target_as_typed="ARd",
        measurement_filters=[
            MeasurementFilter(protocol_hint=None, metric="IC50"),
        ],
    ))
    assert isinstance(res, MetricMiss)
    assert res.available_metrics == []


# ---------------------------------------------------------------------------
# Poly-source registered_by (slice 16) — Users AND Suppliers
# ---------------------------------------------------------------------------


@pytest.fixture
def registrant_executor_world(db):
    """Two compounds: one registered by Alice (user), one supplied by
    Enamine (vendor). Lets us prove "registered by Alice" returns the
    Alice-registered one and "from Enamine" returns the Enamine-supplied
    one."""
    from compounds.registry.models import Supplier
    t = Target.objects.create(name="ARd")
    alice = User.objects.create_user(
        username="alice.jones", email="alice.jones@ncl.ac.uk",
        first_name="Alice", last_name="Jones",
    )
    bob = User.objects.create_user(
        username="bob.smith", email="bob.smith@ncl.ac.uk",
        first_name="Bob", last_name="Smith",
    )
    bob_supplier = Supplier.objects.create(name="Bob Smith", user=bob)
    enamine = Supplier.objects.create(name="Enamine", initials="ENM")
    c_alice = Compound.objects.create(target=t, smiles="CCO", registered_by=alice)
    c_enamine = Compound.objects.create(target=t, smiles="CCN", supplier=enamine)
    # Bob has both a registered compound AND a personal-supplier-supplied
    # one — both should surface under "registered by Bob" thanks to the
    # OR'd filter.
    c_bob_reg = Compound.objects.create(target=t, smiles="CCC", registered_by=bob)
    c_bob_sup = Compound.objects.create(target=t, smiles="CCCC", supplier=bob_supplier)
    return {
        "alice": alice, "bob": bob, "enamine": enamine, "bob_supplier": bob_supplier,
        "c_alice": c_alice, "c_enamine": c_enamine,
        "c_bob_reg": c_bob_reg, "c_bob_sup": c_bob_sup,
    }


def test_registered_by_user_filters_user_compounds(registrant_executor_world):
    res = execute(CompoundSelector(
        registration_target_as_typed="ARd",
        registered_by_as_typed="Alice Jones",
    ))
    assert isinstance(res, CompoundSelection)
    assert res.compound_formatted_ids == [
        registrant_executor_world["c_alice"].formatted_id,
    ]


def test_registered_by_supplier_filters_supplied_compounds(registrant_executor_world):
    res = execute(CompoundSelector(
        registration_target_as_typed="ARd",
        registered_by_as_typed="Enamine",
    ))
    assert isinstance(res, CompoundSelection)
    assert res.compound_formatted_ids == [
        registrant_executor_world["c_enamine"].formatted_id,
    ]


def test_registered_by_user_with_personal_supplier_returns_both(registrant_executor_world):
    """Bob has registered one compound AND supplied one via his own
    personal-supplier. The Q-filter union surfaces BOTH under 'made
    by Bob' — that's what a chemist actually means."""
    res = execute(CompoundSelector(
        registration_target_as_typed="ARd",
        registered_by_as_typed="Bob Smith",
    ))
    assert isinstance(res, CompoundSelection)
    formatted = set(res.compound_formatted_ids)
    assert registrant_executor_world["c_bob_reg"].formatted_id in formatted
    assert registrant_executor_world["c_bob_sup"].formatted_id in formatted


def test_scope_sentence_renders_supplier_with_from_phrasing(registrant_executor_world):
    res = execute(CompoundSelector(
        registration_target_as_typed="ARd",
        registered_by_as_typed="Enamine",
    ))
    assert isinstance(res, CompoundSelection)
    assert "from Enamine" in res.scope_sentence


def test_pinned_supplier_id_bypasses_resolver(registrant_executor_world):
    """The frontend's clarify continuation pins via
    registered_by_supplier_id; the executor honours it directly."""
    res = execute(CompoundSelector(
        registration_target_as_typed="ARd",
        registered_by_as_typed="Enamine",
        registered_by_supplier_id=str(registrant_executor_world["enamine"].pk),
    ))
    assert isinstance(res, CompoundSelection)
    assert res.compound_formatted_ids == [
        registrant_executor_world["c_enamine"].formatted_id,
    ]


def test_unknown_registrant_returns_user_miss(registrant_executor_world):
    res = execute(CompoundSelector(
        registration_target_as_typed="ARd",
        registered_by_as_typed="ZachariahNonExistent",
    ))
    assert isinstance(res, UserMiss)


# ---------------------------------------------------------------------------
# Slice 18: same-name User+Supplier merger end-to-end (Hannah Stewart)
# ---------------------------------------------------------------------------


@pytest.fixture
def hannah_executor_world(db):
    """Hannah Stewart registered three compounds AND has an unlinked
    Supplier 'Hannah Stewart' supplying two more. Slice-18 merges
    them into a Union; the executor's filter ORs all three paths so
    the chemist sees ALL FIVE compounds."""
    from compounds.registry.models import Supplier
    t = Target.objects.create(name="ARd")
    hannah_user = User.objects.create_user(
        username="hannah.stewart", email="hannah.stewart@ncl.ac.uk",
        first_name="Hannah", last_name="Stewart",
    )
    hannah_sup = Supplier.objects.create(name="Hannah Stewart")
    registered = [
        Compound.objects.create(target=t, smiles="CCO", registered_by=hannah_user),
        Compound.objects.create(target=t, smiles="CCN", registered_by=hannah_user),
        Compound.objects.create(target=t, smiles="CCC", registered_by=hannah_user),
    ]
    supplied = [
        Compound.objects.create(target=t, smiles="CCCC", supplier=hannah_sup),
        Compound.objects.create(target=t, smiles="CCCN", supplier=hannah_sup),
    ]
    return {
        "user": hannah_user, "supplier": hannah_sup,
        "registered": registered, "supplied": supplied,
    }


def test_union_filter_returns_both_registered_and_supplied(hannah_executor_world):
    """The Hannah Stewart bug: pre-slice-18, "made by Hannah Stewart"
    asked the chemist to pick a role; post-slice-18 it auto-resolves
    to a Union and returns ALL FIVE compounds."""
    res = execute(CompoundSelector(
        registration_target_as_typed="ARd",
        registered_by_as_typed="Hannah Stewart",
    ))
    assert isinstance(res, CompoundSelection)
    expected = {c.formatted_id for c in (
        hannah_executor_world["registered"] + hannah_executor_world["supplied"]
    )}
    assert set(res.compound_formatted_ids) == expected


def test_union_scope_sentence_says_registered_or_supplied(hannah_executor_world):
    res = execute(CompoundSelector(
        registration_target_as_typed="ARd",
        registered_by_as_typed="Hannah Stewart",
    ))
    assert isinstance(res, CompoundSelection)
    assert "registered or supplied by Hannah Stewart" in res.scope_sentence


def test_union_pinned_continuation_returns_same_set(hannah_executor_world):
    """Frontend's clarify continuation sends both pin ids — the
    executor honours that as a Union without re-running the resolver."""
    res = execute(CompoundSelector(
        registration_target_as_typed="ARd",
        registered_by_as_typed="ignored",
        registered_by_id=str(hannah_executor_world["user"].pk),
        registered_by_supplier_id=str(hannah_executor_world["supplier"].pk),
    ))
    assert isinstance(res, CompoundSelection)
    expected = {c.formatted_id for c in (
        hannah_executor_world["registered"] + hannah_executor_world["supplied"]
    )}
    assert set(res.compound_formatted_ids) == expected


# ---------------------------------------------------------------------------
# Ranking — slice 19
# ---------------------------------------------------------------------------


@pytest.fixture
def scorecard_world(db):
    """A target with a configured scorecard (single protocol axis,
    log-scale, lower-is-better) and four compounds with a spread of
    HTRF IC50 measurements. Ranking should put the most potent first."""
    t = Target.objects.create(
        name="ARd",
        scorecard_config={
            "axes": [
                {
                    "kind": "protocol",
                    "label": "HTRF IC50",
                    "protocol_id": None,        # set below once we have the protocol
                    "target_value": 10,
                    "poor_value": 10000,
                    "threshold_scale": "log",
                },
            ],
        },
    )
    p = Protocol.objects.create(
        name="ARd HTRF",
        target_value=10, poor_value=10000, threshold_scale="log",
    )
    # Patch the scorecard config to point at the actual protocol pk.
    t.scorecard_config["axes"][0]["protocol_id"] = str(p.id)
    t.save()
    a = Assay.objects.create(protocol=p, target=t)

    # Four compounds at different potencies. With log normalisation
    # between target=10 and poor=10000, the t-scores should be
    # ordered: 10 nM (t=1.0) > 100 nM (t=2/3) > 1000 nM (t=1/3) > 5000 nM (t≈0.1).
    compounds = []
    for i, ic50 in enumerate([5000, 100, 1000, 10]):
        c = Compound.objects.create(target=t, smiles="CCO" + "C" * i)
        _mk_result(a, c, status="valid",
                   results={"KPI": "IC50", "IC50": ic50, "kpi_unit": "nM"},
                   row=i)
        compounds.append((c, ic50))

    return {"target": t, "protocol": p, "assay": a, "compounds": compounds}


def test_ranking_orders_compounds_by_scorecard_score(scorecard_world):
    """The most potent compound should come first; the weakest last."""
    res = execute(CompoundSelector(
        registration_target_as_typed="ARd",
        rank_by="scorecard",
    ))
    assert isinstance(res, CompoundSelection)
    # Compounds sorted by their potency: 10 nM > 100 nM > 1000 nM > 5000 nM.
    by_potency = sorted(scorecard_world["compounds"], key=lambda pair: pair[1])
    expected_order = [c.formatted_id for c, _ic50 in by_potency]
    assert res.compound_formatted_ids == expected_order


def test_ranking_caps_to_top_n(scorecard_world):
    res = execute(CompoundSelector(
        registration_target_as_typed="ARd",
        rank_by="scorecard",
        rank_top_n=2,
    ))
    assert isinstance(res, CompoundSelection)
    assert res.n_matched == 2
    by_potency = sorted(scorecard_world["compounds"], key=lambda pair: pair[1])
    expected_top_2 = [c.formatted_id for c, _ic50 in by_potency[:2]]
    assert res.compound_formatted_ids == expected_top_2


def test_ranking_default_top_n_is_twenty(scorecard_world):
    """Four compounds in scope, default cap → all four return (cap is
    20 by default, so it doesn't trim a small set)."""
    res = execute(CompoundSelector(
        registration_target_as_typed="ARd",
        rank_by="scorecard",
    ))
    assert isinstance(res, CompoundSelection)
    assert res.n_matched == 4


def test_ranking_composes_with_threshold_filter(scorecard_world):
    """Filter narrows first; rank applies only to survivors."""
    res = execute(CompoundSelector(
        registration_target_as_typed="ARd",
        measurement_filters=[
            MeasurementFilter(
                protocol_hint="HTRF", metric="IC50",
                threshold=Threshold(op="<", value=500, unit="nM"),
            ),
        ],
        rank_by="scorecard",
    ))
    assert isinstance(res, CompoundSelection)
    # Only the 10 nM and 100 nM compounds pass the filter.
    assert res.n_matched == 2
    # Ranked: 10 nM first.
    by_potency = sorted(scorecard_world["compounds"], key=lambda pair: pair[1])
    expected = [c.formatted_id for c, ic50 in by_potency if ic50 < 500]
    assert res.compound_formatted_ids == expected


def test_ranking_scope_sentence_includes_top_n(scorecard_world):
    res = execute(CompoundSelector(
        registration_target_as_typed="ARd",
        rank_by="scorecard",
        rank_top_n=2,
    ))
    assert isinstance(res, CompoundSelection)
    assert res.scope_sentence.lower().startswith("top 2 by scorecard")


def test_ranking_without_scorecard_returns_spec_error(db):
    """Target has no scorecard configured → SpecError, not a silent fallback."""
    Target.objects.create(name="UnscoredTarget")  # scorecard_config defaults to None
    res = execute(CompoundSelector(
        registration_target_as_typed="UnscoredTarget",
        rank_by="scorecard",
    ))
    assert isinstance(res, SpecError)
    assert "scorecard" in res.message.lower()


def test_unknown_rank_by_returns_spec_error(scorecard_world):
    res = execute(CompoundSelector(
        registration_target_as_typed="ARd",
        rank_by="popularity_contest",
    ))
    assert isinstance(res, SpecError)


def test_similarity_narrows_to_neighbours_of_anchor(db):
    """End-to-end: 'compounds similar to NCL-X' returns only compounds
    above the threshold — the anchor itself plus close analogues."""
    from compounds.nlp.spec import CompoundSelection
    t = Target.objects.create(name="ARd")
    # Aspirin (anchor), salicylic acid (close analogue, ~0.45 Tanimoto),
    # caffeine (unrelated, ~0.09 Tanimoto).
    anchor = Compound.objects.create(
        target=t, smiles="CC(=O)Oc1ccccc1C(=O)O", reg_number=10001,
    )
    analogue = Compound.objects.create(
        target=t, smiles="Oc1ccccc1C(=O)O", reg_number=10002,
    )
    distant = Compound.objects.create(
        target=t, smiles="Cn1cnc2c1c(=O)n(C)c(=O)n2C", reg_number=10003,
    )
    res = execute(CompoundSelector(
        registration_target_as_typed="ARd",
        similar_to_as_typed=["NCL-00010001"],
        similar_threshold=0.4,
    ))
    assert isinstance(res, CompoundSelection)
    formatted = set(res.compound_formatted_ids)
    assert anchor.formatted_id in formatted
    assert analogue.formatted_id in formatted
    assert distant.formatted_id not in formatted


def test_similarity_uses_default_threshold_when_unspecified(db):
    """No threshold → backend default of 0.7. Self-similarity is 1.0
    so the anchor lands; an unrelated compound is excluded."""
    from compounds.nlp.spec import CompoundSelection
    t = Target.objects.create(name="ARd")
    anchor = Compound.objects.create(
        target=t, smiles="CC(=O)Oc1ccccc1C(=O)O", reg_number=10001,
    )
    Compound.objects.create(
        target=t, smiles="Cn1cnc2c1c(=O)n(C)c(=O)n2C", reg_number=10002,
    )
    res = execute(CompoundSelector(
        registration_target_as_typed="ARd",
        similar_to_as_typed=["NCL-00010001"],
    ))
    assert isinstance(res, CompoundSelection)
    assert res.compound_formatted_ids == [anchor.formatted_id]


def test_similarity_unknown_anchor_returns_compound_miss(db):
    """An anchor that doesn't resolve surfaces as CompoundMiss with the
    ref_index, same as compound_refs_as_typed."""
    Target.objects.create(name="ARd")
    res = execute(CompoundSelector(
        registration_target_as_typed="ARd",
        similar_to_as_typed=["NCL-99999999"],
    ))
    assert isinstance(res, CompoundMiss)
    assert res.field == FIELD_COMPOUND_REF


def test_similarity_composes_with_other_predicates(db):
    """Similar-to is a compound-level narrowing predicate that composes
    cleanly with measurement filters: filter to similar compounds AND
    those passing a threshold."""
    from compounds.nlp.spec import CompoundSelection
    t = Target.objects.create(name="ARd")
    # Aspirin anchor (NCL-10001), close analogues (salicylic, methyl
    # salicylate), unrelated potent (caffeine).
    anchor = Compound.objects.create(
        target=t, smiles="CC(=O)Oc1ccccc1C(=O)O", reg_number=10001,
    )
    analogue_potent = Compound.objects.create(
        target=t, smiles="Oc1ccccc1C(=O)O", reg_number=10002,
    )
    analogue_weak = Compound.objects.create(
        target=t, smiles="COC(=O)c1ccccc1O", reg_number=10003,
    )
    distant_potent = Compound.objects.create(
        target=t, smiles="Cn1cnc2c1c(=O)n(C)c(=O)n2C", reg_number=10004,
    )
    p = Protocol.objects.create(name="ARd HTRF")
    a = Assay.objects.create(protocol=p, target=t)
    _mk_result(a, anchor, status="valid",
               results={"KPI": "IC50", "IC50": 5.0, "kpi_unit": "nM"}, row=0)
    _mk_result(a, analogue_potent, status="valid",
               results={"KPI": "IC50", "IC50": 8.0, "kpi_unit": "nM"}, row=1)
    _mk_result(a, analogue_weak, status="valid",
               results={"KPI": "IC50", "IC50": 500.0, "kpi_unit": "nM"}, row=2)
    _mk_result(a, distant_potent, status="valid",
               results={"KPI": "IC50", "IC50": 1.0, "kpi_unit": "nM"}, row=3)
    res = execute(CompoundSelector(
        registration_target_as_typed="ARd",
        similar_to_as_typed=["NCL-00010001"],
        similar_threshold=0.4,    # aspirin/salicylic ~0.45, others below
        measurement_filters=[
            MeasurementFilter(
                protocol_hint="HTRF", metric="IC50",
                threshold=Threshold(op="<", value=100, unit="nM"),
            ),
        ],
    ))
    assert isinstance(res, CompoundSelection)
    formatted = set(res.compound_formatted_ids)
    # Similar AND under threshold.
    assert anchor.formatted_id in formatted
    assert analogue_potent.formatted_id in formatted
    # Below similarity threshold (~0.35).
    assert analogue_weak.formatted_id not in formatted
    # Passes the threshold but not similar.
    assert distant_potent.formatted_id not in formatted


def test_similarity_scope_sentence_includes_anchors(db):
    from compounds.nlp.spec import CompoundSelection
    t = Target.objects.create(name="ARd")
    Compound.objects.create(target=t, smiles="c1ccncc1", reg_number=10001)
    res = execute(CompoundSelector(
        registration_target_as_typed="ARd",
        similar_to_as_typed=["NCL-00010001"],
        similar_threshold=0.5,
    ))
    assert isinstance(res, CompoundSelection)
    assert "similar to NCL-00010001" in res.scope_sentence
    assert "Tanimoto" in res.scope_sentence


def test_ranking_skips_compounds_with_no_signal(db):
    """A compound with no measurements scores None; it's dropped from
    the ranking rather than appearing at the bottom."""
    t = Target.objects.create(
        name="ARd",
        scorecard_config={
            "axes": [{
                "kind": "protocol", "label": "HTRF",
                "protocol_id": None,
                "target_value": 10, "poor_value": 10000,
                "threshold_scale": "log",
            }],
        },
    )
    p = Protocol.objects.create(
        name="ARd HTRF",
        target_value=10, poor_value=10000, threshold_scale="log",
    )
    t.scorecard_config["axes"][0]["protocol_id"] = str(p.id)
    t.save()
    a = Assay.objects.create(protocol=p, target=t)
    # One compound with data, one without.
    c_scored = Compound.objects.create(target=t, smiles="CCO")
    c_silent = Compound.objects.create(target=t, smiles="CCN")
    _mk_result(a, c_scored, status="valid",
               results={"KPI": "IC50", "IC50": 50.0, "kpi_unit": "nM"})

    res = execute(CompoundSelector(
        registration_target_as_typed="ARd",
        rank_by="scorecard",
    ))
    assert isinstance(res, CompoundSelection)
    # Only the scored compound makes it into the ranking.
    assert res.compound_formatted_ids == [c_scored.formatted_id]
    assert c_silent.formatted_id not in res.compound_formatted_ids
