"""DRF view tests for POST /api/compounds/nlp/query/ (post-pivot).

Exercises the HTTP surface: feature flag, body validation, prompt path
through mocked parse_prompt → execute, continuation path via selector,
pinning-ID round trip, error HTTP mappings, daily cap.
"""

from __future__ import annotations

from unittest.mock import patch

import pytest
from django.contrib.auth import get_user_model
from django.core.cache import cache
from rest_framework.test import APIClient

from compounds.assays.models import AnalysisResult, Assay, DataSeries, Protocol
from compounds.nlp.spec import (
    CompoundSelector,
    MeasurementFilter,
    NotAQuery,
    ParseError,
    Threshold,
)
from compounds.nlp.view import DAILY_CAP, ENV_FLAG
from compounds.registry.models import Compound, Gene, Target

User = get_user_model()

URL = "/api/compounds/nlp/query/"


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def feature_on(monkeypatch):
    monkeypatch.setenv(ENV_FLAG, "true")


@pytest.fixture
def clear_cache():
    cache.clear()
    yield
    cache.clear()


@pytest.fixture
def user(db):
    return User.objects.create_user(username="alice", email="alice@example.org")


@pytest.fixture
def client(user):
    c = APIClient()
    c.force_authenticate(user=user)
    return c


@pytest.fixture
def world(db):
    ar = Target.objects.create(name="AR degraders")
    compound = Compound.objects.create(target=ar, smiles="CCO")
    protocol = Protocol.objects.create(name="AR binding HTRF")
    assay = Assay.objects.create(protocol=protocol, target=ar)
    ar_row = AnalysisResult.objects.create(
        status="valid",
        results={"KPI": "IC50", "IC50": 5.0, "kpi_unit": "nM"},
    )
    DataSeries.objects.create(
        assay=assay, compound=compound, row=0, start_column=0, end_column=10,
        analysis=ar_row,
    )
    return {"ar": ar, "compound": compound, "protocol": protocol}


def _selector_body(**overrides):
    base = {
        "registration_target_as_typed": "AR degraders",
        "assay_target_as_typed": None,
        "measurement_filters": [],
    }
    base.update(overrides)
    return base


# ---------------------------------------------------------------------------
# Feature flag
# ---------------------------------------------------------------------------


def test_flag_off_returns_404(client, clear_cache, monkeypatch, db):
    monkeypatch.delenv(ENV_FLAG, raising=False)
    resp = client.post(URL, {"prompt": "anything"}, format="json")
    assert resp.status_code == 404
    assert resp.data["kind"] == "disabled"


# ---------------------------------------------------------------------------
# Body validation
# ---------------------------------------------------------------------------


def test_missing_body_returns_400(client, feature_on, clear_cache, db):
    resp = client.post(URL, {}, format="json")
    assert resp.status_code == 400
    assert resp.data["kind"] == "bad_request"


def test_both_prompt_and_selector_returns_400(client, feature_on, clear_cache, db):
    resp = client.post(URL, {"prompt": "x", "selector": _selector_body()}, format="json")
    assert resp.status_code == 400


def test_non_string_prompt_returns_400(client, feature_on, clear_cache, db):
    resp = client.post(URL, {"prompt": 42}, format="json")
    assert resp.status_code == 400


def test_malformed_selector_threshold_returns_400(client, feature_on, clear_cache, db):
    body = _selector_body(measurement_filters=[
        {"protocol_hint": "HTRF", "metric": "IC50",
         "threshold": {"op": "<"}},  # missing value
    ])
    resp = client.post(URL, {"selector": body}, format="json")
    assert resp.status_code == 400


# ---------------------------------------------------------------------------
# Happy path — prompt → parse → execute → redirect
# ---------------------------------------------------------------------------


def test_prompt_happy_path_returns_redirect_url(client, feature_on, clear_cache, world):
    parsed = CompoundSelector(
        registration_target_as_typed="AR degraders",
        measurement_filters=[
            MeasurementFilter(
                protocol_hint="HTRF",
                metric="IC50",
                threshold=Threshold(op="<", value=10.0, unit="nM"),
            ),
        ],
    )
    with patch("compounds.nlp.view.parse_prompt", return_value=parsed):
        resp = client.post(URL, {
            "prompt": "ARd compounds with HTRF IC50 < 10 nM",
        }, format="json")

    assert resp.status_code == 200
    assert resp.data["status"] == "selection"
    assert resp.data["n_matched"] == 1
    assert resp.data["compound_formatted_ids"] == [world["compound"].formatted_id]
    assert resp.data["target_names"] == ["AR degraders"]
    assert resp.data["protocol_names"] == ["AR binding HTRF"]

    redirect_url = resp.data["redirect_url"]
    assert redirect_url.startswith("/assays/aggregate?")
    assert "targets=AR%20degraders" in redirect_url or "targets=AR degraders" in redirect_url
    assert "format=cards" in redirect_url
    assert f"compound={world['compound'].formatted_id}" in redirect_url
    assert "protocols=AR%20binding%20HTRF" in redirect_url or "protocols=AR binding HTRF" in redirect_url


def test_selection_with_no_filters_omits_protocols_param(client, feature_on, clear_cache, world):
    parsed = CompoundSelector(registration_target_as_typed="AR degraders")
    with patch("compounds.nlp.view.parse_prompt", return_value=parsed):
        resp = client.post(URL, {"prompt": "show all ARd compounds"}, format="json")
    assert resp.status_code == 200
    assert resp.data["status"] == "selection"
    assert resp.data["protocol_names"] == []
    assert "protocols=" not in resp.data["redirect_url"]


def test_empty_selection_still_200(client, feature_on, clear_cache, world):
    """Zero matching compounds is a valid result — the UI shows 'Found 0'."""
    parsed = CompoundSelector(
        registration_target_as_typed="AR degraders",
        measurement_filters=[
            MeasurementFilter(
                protocol_hint="HTRF",
                metric="IC50",
                threshold=Threshold(op="<", value=0.001, unit="nM"),
            ),
        ],
    )
    with patch("compounds.nlp.view.parse_prompt", return_value=parsed):
        resp = client.post(URL, {"prompt": "impossible"}, format="json")
    assert resp.status_code == 200
    assert resp.data["status"] == "selection"
    assert resp.data["n_matched"] == 0
    # redirect_url still contains the protocol (user can see what was filtered)
    # but the compound list is empty.
    assert "compound=" not in resp.data["redirect_url"]


# ---------------------------------------------------------------------------
# Continuation path — selector instead of prompt, no LLM call
# ---------------------------------------------------------------------------


def test_selector_continuation_skips_llm(client, feature_on, clear_cache, world):
    body = _selector_body(measurement_filters=[
        {"protocol_hint": "HTRF", "metric": "IC50",
         "threshold": {"op": "<", "value": 10, "unit": "nM"}},
    ])
    with patch("compounds.nlp.view.parse_prompt") as mock_parse:
        resp = client.post(URL, {"selector": body}, format="json")
    mock_parse.assert_not_called()
    assert resp.status_code == 200
    assert resp.data["status"] == "selection"


def test_selector_with_pinning_ids_bypasses_resolution(client, feature_on, clear_cache, db):
    """Ambiguous target via shared gene → pin with registration_target_id."""
    gene = Gene.objects.create(symbol="MYC")
    t1 = Target.objects.create(name="Myc-Aur")
    t2 = Target.objects.create(name="Myc regulation")
    t1.genes.add(gene)
    t2.genes.add(gene)
    compound = Compound.objects.create(target=t1, smiles="CCO")
    protocol = Protocol.objects.create(name="Myc HTRF")
    assay = Assay.objects.create(protocol=protocol, target=t1)
    ar = AnalysisResult.objects.create(
        status="valid", results={"KPI": "IC50", "IC50": 5.0, "kpi_unit": "nM"},
    )
    DataSeries.objects.create(
        assay=assay, compound=compound, row=0, start_column=0, end_column=10,
        analysis=ar,
    )
    body = _selector_body(
        registration_target_as_typed="MYC",
        registration_target_id=str(t1.id),
        measurement_filters=[{"protocol_hint": "Myc HTRF", "metric": "IC50", "threshold": None}],
    )
    resp = client.post(URL, {"selector": body}, format="json")
    assert resp.status_code == 200
    assert resp.data["status"] == "selection"
    assert resp.data["n_matched"] == 1


# ---------------------------------------------------------------------------
# Clarify round-trip (partial_selector included, pick + re-POST skips LLM)
# ---------------------------------------------------------------------------


def test_clarify_response_includes_partial_selector(client, feature_on, clear_cache, db):
    gene = Gene.objects.create(symbol="MYC")
    t1 = Target.objects.create(name="Myc-Aur")
    t2 = Target.objects.create(name="Myc regulation")
    t1.genes.add(gene)
    t2.genes.add(gene)
    parsed = CompoundSelector(
        registration_target_as_typed="MYC",
        measurement_filters=[
            MeasurementFilter(protocol_hint="HTRF", metric="IC50"),
        ],
    )
    with patch("compounds.nlp.view.parse_prompt", return_value=parsed):
        resp = client.post(URL, {"prompt": "MYC compounds"}, format="json")
    assert resp.data["status"] == "clarify"
    assert "partial_selector" in resp.data
    assert resp.data["partial_selector"]["registration_target_as_typed"] == "MYC"
    assert len(resp.data["partial_selector"]["measurement_filters"]) == 1


def test_clarify_response_round_trips_into_continuation(client, feature_on, clear_cache, db):
    """Pick a target chip → build pinned selector → re-POST → execute, no LLM."""
    gene = Gene.objects.create(symbol="MYC")
    t1 = Target.objects.create(name="Myc-Aur")
    t2 = Target.objects.create(name="Myc regulation")
    t1.genes.add(gene)
    t2.genes.add(gene)
    compound = Compound.objects.create(target=t1, smiles="CCO")
    protocol = Protocol.objects.create(name="MYC HTRF")
    assay = Assay.objects.create(protocol=protocol, target=t1)
    ar = AnalysisResult.objects.create(
        status="valid", results={"KPI": "IC50", "IC50": 5.0, "kpi_unit": "nM"},
    )
    DataSeries.objects.create(
        assay=assay, compound=compound, row=0, start_column=0, end_column=10,
        analysis=ar,
    )

    parsed = CompoundSelector(
        registration_target_as_typed="MYC",
        measurement_filters=[
            MeasurementFilter(protocol_hint="HTRF", metric="IC50"),
        ],
    )
    with patch("compounds.nlp.view.parse_prompt", return_value=parsed):
        resp1 = client.post(URL, {"prompt": "MYC compounds"}, format="json")
    assert resp1.data["status"] == "clarify"

    partial = dict(resp1.data["partial_selector"])
    partial["registration_target_id"] = str(t1.id)
    with patch("compounds.nlp.view.parse_prompt") as mock_parse:
        resp2 = client.post(URL, {"selector": partial}, format="json")
    mock_parse.assert_not_called()
    assert resp2.data["status"] == "selection"
    assert resp2.data["n_matched"] == 1


# ---------------------------------------------------------------------------
# Error mapping
# ---------------------------------------------------------------------------


def test_target_miss_returns_200_with_miss(client, feature_on, clear_cache, world):
    parsed = CompoundSelector(registration_target_as_typed="zzzznotareal")
    with patch("compounds.nlp.view.parse_prompt", return_value=parsed):
        resp = client.post(URL, {"prompt": "whatever"}, format="json")
    assert resp.status_code == 200
    assert resp.data["status"] == "miss"


def test_fully_empty_selector_returns_400_scope_error(client, feature_on, clear_cache, db):
    """A selector with literally no narrowing predicate still errors,
    with a user-friendly "narrow by something" message."""
    parsed = CompoundSelector()
    with patch("compounds.nlp.view.parse_prompt", return_value=parsed):
        resp = client.post(URL, {"prompt": "whatever"}, format="json")
    assert resp.status_code == 400
    assert resp.data["kind"] == "scope"
    assert "narrow" in resp.data["message"].lower()


def test_target_less_selector_with_filter_runs(client, feature_on, clear_cache, world):
    """A target-less selector with a measurement filter runs unscoped —
    the ScopeError only fires for truly-empty selectors."""
    parsed = CompoundSelector(measurement_filters=[
        MeasurementFilter(
            protocol_hint="HTRF", metric="IC50",
            threshold=Threshold(op="<", value=10.0, unit="nM"),
        ),
    ])
    with patch("compounds.nlp.view.parse_prompt", return_value=parsed):
        resp = client.post(URL, {"prompt": "anything with HTRF IC50 < 10 nM"}, format="json")
    assert resp.status_code == 200
    assert resp.data["status"] == "selection"


def test_not_a_query_returns_200(client, feature_on, clear_cache, db):
    with patch(
        "compounds.nlp.view.parse_prompt",
        return_value=NotAQuery(reason="Asked for narrative, not a selection."),
    ):
        resp = client.post(URL, {"prompt": "summarise the SAR"}, format="json")
    assert resp.status_code == 200
    assert resp.data["status"] == "not_a_query"


def test_parse_error_returns_502(client, feature_on, clear_cache, db):
    with patch(
        "compounds.nlp.view.parse_prompt",
        return_value=ParseError(message="LLM returned junk", raw="{no"),
    ):
        resp = client.post(URL, {"prompt": "anything"}, format="json")
    assert resp.status_code == 502
    assert resp.data["kind"] == "parse"


# ---------------------------------------------------------------------------
# Daily cap
# ---------------------------------------------------------------------------


def test_daily_cap_blocks_after_limit(client, feature_on, clear_cache, world):
    parsed = CompoundSelector(registration_target_as_typed="AR degraders")
    from compounds.nlp.view import _cap_key
    cache.set(_cap_key(client.handler._force_user.pk), DAILY_CAP, timeout=3600)
    with patch("compounds.nlp.view.parse_prompt", return_value=parsed):
        resp = client.post(URL, {"prompt": "anything"}, format="json")
    assert resp.status_code == 429
    assert resp.data["kind"] == "rate_limited"


def test_separate_users_have_separate_caps(feature_on, clear_cache, world, db):
    alice = User.objects.create_user(username="alice2", email="a2@example.org")
    bob = User.objects.create_user(username="bob2", email="b2@example.org")
    from compounds.nlp.view import _cap_key
    cache.set(_cap_key(alice.pk), DAILY_CAP, timeout=3600)

    parsed = CompoundSelector(registration_target_as_typed="AR degraders")
    alice_c = APIClient(); alice_c.force_authenticate(user=alice)
    bob_c = APIClient(); bob_c.force_authenticate(user=bob)
    with patch("compounds.nlp.view.parse_prompt", return_value=parsed):
        assert alice_c.post(URL, {"prompt": "x"}, format="json").status_code == 429
        assert bob_c.post(URL, {"prompt": "x"}, format="json").status_code == 200
