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
    CompoundSelector,
    CompoundSelection,
    DateRange,
    FIELD_ASSAYED_BY,
    FIELD_PROTOCOL_HINT,
    FIELD_REGISTERED_BY,
    FIELD_REGISTRATION_TARGET,
    MeasurementFilter,
    ProtocolClarify,
    ProtocolMiss,
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
