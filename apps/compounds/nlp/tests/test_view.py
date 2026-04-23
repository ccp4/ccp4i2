"""DRF view tests for POST /api/compounds/nlp/query/ (slice 6).

Exercises the HTTP surface end-to-end: feature-flag gate, auth, the prompt
path (parse_prompt → execute) with the LLM mocked, the spec path
(clarify continuation) with pinning-ID honour, error mapping onto HTTP
status codes, and the daily-call cap.
"""

from __future__ import annotations

import os
from unittest.mock import patch

import pytest
from django.contrib.auth import get_user_model
from django.core.cache import cache
from rest_framework.test import APIClient

from compounds.assays.models import AnalysisResult, Assay, DataSeries, Protocol
from compounds.nlp.spec import (
    NotAQuery,
    ParseError,
    QuerySpec,
    Threshold,
)
from compounds.nlp.view import DAILY_CAP, ENV_FLAG
from compounds.registry.models import Compound, MolecularProperties, Target

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
    """Minimal fixture: one AR target + compound + HTRF protocol + valid
    measurement. Enough to exercise the happy path end-to-end."""
    ar = Target.objects.create(name="AR degraders")
    compound = Compound.objects.create(target=ar, smiles="CCO")
    MolecularProperties.objects.filter(compound=compound).update(
        molecular_weight=400.5, clogp=3.2, hbd=2, hba=5,
        heavy_atom_count=28, tpsa=75.1, rotatable_bonds=6, fraction_sp3=0.35,
    )
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
    return {
        "ar": ar,
        "compound": compound,
        "protocol": protocol,
    }


def _spec_body(**overrides):
    """Build a valid spec dict with sensible defaults, overridable per-test."""
    base = {
        "registration_target_as_typed": "AR degraders",
        "assay_target_as_typed": "AR degraders",
        "protocol_hint": "HTRF",
        "metric": "IC50",
        "threshold": None,
        "columns": "phys_chem",
        "columns_explicit": None,
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


def test_flag_explicit_false_returns_404(client, clear_cache, monkeypatch, db):
    monkeypatch.setenv(ENV_FLAG, "false")
    resp = client.post(URL, {"prompt": "anything"}, format="json")
    assert resp.status_code == 404


# ---------------------------------------------------------------------------
# Body validation
# ---------------------------------------------------------------------------


def test_missing_body_returns_400(client, feature_on, clear_cache, db):
    resp = client.post(URL, {}, format="json")
    assert resp.status_code == 400
    assert resp.data["kind"] == "bad_request"


def test_both_prompt_and_spec_returns_400(client, feature_on, clear_cache, db):
    resp = client.post(URL, {"prompt": "x", "spec": _spec_body()}, format="json")
    assert resp.status_code == 400


def test_non_string_prompt_returns_400(client, feature_on, clear_cache, db):
    resp = client.post(URL, {"prompt": 42}, format="json")
    assert resp.status_code == 400


def test_malformed_spec_threshold_returns_400(client, feature_on, clear_cache, db):
    resp = client.post(URL, {
        "spec": _spec_body(threshold={"op": "<"}),  # missing value
    }, format="json")
    assert resp.status_code == 400


# ---------------------------------------------------------------------------
# Happy path — prompt → parse → execute
# ---------------------------------------------------------------------------


def test_prompt_happy_path(client, feature_on, clear_cache, world):
    parsed = QuerySpec(
        registration_target_as_typed="AR degraders",
        assay_target_as_typed="AR degraders",
        protocol_hint="HTRF",
        metric="IC50",
        threshold=Threshold(op="<", value=10.0, unit="nM"),
    )
    with patch("compounds.nlp.view.parse_prompt", return_value=parsed):
        resp = client.post(URL, {
            "prompt": "phys chem for ARd compounds with HTRF IC50 < 10 nM",
        }, format="json")

    assert resp.status_code == 200
    assert resp.data["status"] == "table"
    assert len(resp.data["rows"]) == 1
    row = resp.data["rows"][0]
    assert row["formatted_id"] == world["compound"].formatted_id
    assert row["value"] == 5.0
    assert row["value_unit"] == "nM"
    assert row["properties"]["molecular_weight"] == pytest.approx(400.5)
    assert "scope_sentence" in resp.data
    assert "AR degraders" in resp.data["scope_sentence"]


# ---------------------------------------------------------------------------
# Happy path — spec (continuation) → execute, no LLM
# ---------------------------------------------------------------------------


def test_spec_continuation_skips_llm(client, feature_on, clear_cache, world):
    with patch("compounds.nlp.view.parse_prompt") as mock_parse:
        resp = client.post(URL, {"spec": _spec_body()}, format="json")

    assert resp.status_code == 200
    assert resp.data["status"] == "table"
    mock_parse.assert_not_called()


def test_spec_with_pinning_ids_bypasses_resolution(client, feature_on, clear_cache, db):
    # Build two same-gene targets so the typed string "MYC" would clarify —
    # but the pinned id picks one unambiguously.
    from compounds.registry.models import Gene
    gene = Gene.objects.create(symbol="MYC")
    t1 = Target.objects.create(name="Myc-Aur")
    t2 = Target.objects.create(name="Myc regulation")
    t1.genes.add(gene)
    t2.genes.add(gene)
    compound = Compound.objects.create(target=t1, smiles="CCO")
    protocol = Protocol.objects.create(name="Myc HTRF")
    assay = Assay.objects.create(protocol=protocol, target=t1)
    ar = AnalysisResult.objects.create(
        status="valid",
        results={"KPI": "IC50", "IC50": 5.0, "kpi_unit": "nM"},
    )
    DataSeries.objects.create(
        assay=assay, compound=compound, row=0, start_column=0, end_column=10,
        analysis=ar,
    )

    resp = client.post(URL, {
        "spec": _spec_body(
            registration_target_as_typed="MYC",
            assay_target_as_typed="MYC",
            registration_target_id=str(t1.id),
            assay_target_id=str(t1.id),
            protocol_hint="HTRF",
        ),
    }, format="json")

    assert resp.status_code == 200
    assert resp.data["status"] == "table"
    assert len(resp.data["rows"]) == 1


# ---------------------------------------------------------------------------
# Error mapping
# ---------------------------------------------------------------------------


def test_target_clarify_returns_clarify_body(client, feature_on, clear_cache, db):
    from compounds.registry.models import Gene
    gene = Gene.objects.create(symbol="MYC")
    t1 = Target.objects.create(name="Myc-Aur")
    t2 = Target.objects.create(name="Myc regulation")
    t1.genes.add(gene)
    t2.genes.add(gene)

    parsed = QuerySpec(
        registration_target_as_typed="MYC",
        assay_target_as_typed="MYC",
        protocol_hint="HTRF",
        metric="IC50",
    )
    with patch("compounds.nlp.view.parse_prompt", return_value=parsed):
        resp = client.post(URL, {"prompt": "MYC compounds"}, format="json")

    assert resp.status_code == 200
    assert resp.data["status"] == "clarify"
    assert resp.data["field"] == "registration_target_as_typed"
    candidate_names = {c["name"] for c in resp.data["candidates"]}
    assert candidate_names == {"Myc-Aur", "Myc regulation"}
    # partial_spec lets the client re-POST with a pinned ID without re-calling
    # the LLM (§18.3 slice-6 decision — decision 5 forbids LLM on continuation).
    assert "partial_spec" in resp.data
    assert resp.data["partial_spec"]["registration_target_as_typed"] == "MYC"
    assert resp.data["partial_spec"]["protocol_hint"] == "HTRF"
    assert resp.data["partial_spec"]["metric"] == "IC50"


def test_clarify_response_round_trips_into_continuation(client, feature_on, clear_cache, db):
    """Pick a chip → build a pinned spec from partial_spec → re-POST → execute."""
    from compounds.registry.models import Gene
    gene = Gene.objects.create(symbol="MYC")
    t1 = Target.objects.create(name="Myc-Aur")
    t2 = Target.objects.create(name="Myc regulation")
    t1.genes.add(gene)
    t2.genes.add(gene)
    # Wire t1 up with a compound + protocol + measurement so the continuation
    # actually resolves to a table.
    compound = Compound.objects.create(target=t1, smiles="CCO")
    MolecularProperties.objects.filter(compound=compound).update(
        molecular_weight=400.0, clogp=3.0, hbd=1, hba=3,
        heavy_atom_count=28, tpsa=60.0, rotatable_bonds=5, fraction_sp3=0.3,
    )
    protocol = Protocol.objects.create(name="MYC HTRF")
    assay = Assay.objects.create(protocol=protocol, target=t1)
    ar = AnalysisResult.objects.create(
        status="valid", results={"KPI": "IC50", "IC50": 5.0, "kpi_unit": "nM"},
    )
    DataSeries.objects.create(
        assay=assay, compound=compound, row=0, start_column=0, end_column=10,
        analysis=ar,
    )

    # First round: LLM-driven, returns Clarify.
    parsed = QuerySpec(
        registration_target_as_typed="MYC",
        assay_target_as_typed="MYC",
        protocol_hint="HTRF",
        metric="IC50",
    )
    with patch("compounds.nlp.view.parse_prompt", return_value=parsed):
        resp1 = client.post(URL, {"prompt": "MYC compounds"}, format="json")
    assert resp1.data["status"] == "clarify"

    # Second round: client pins t1 and re-POSTs — must skip the LLM.
    partial = dict(resp1.data["partial_spec"])
    partial["registration_target_id"] = str(t1.id)
    partial["assay_target_id"] = str(t1.id)
    with patch("compounds.nlp.view.parse_prompt") as mock_parse:
        resp2 = client.post(URL, {"spec": partial}, format="json")
    mock_parse.assert_not_called()
    assert resp2.data["status"] == "table"
    assert len(resp2.data["rows"]) == 1


def test_target_miss_returns_miss_body(client, feature_on, clear_cache, world):
    parsed = QuerySpec(
        registration_target_as_typed="zzzznotarealtarget",
        protocol_hint="HTRF",
        metric="IC50",
    )
    with patch("compounds.nlp.view.parse_prompt", return_value=parsed):
        resp = client.post(URL, {"prompt": "whatever"}, format="json")
    assert resp.status_code == 200
    assert resp.data["status"] == "miss"


def test_scope_error_returns_400(client, feature_on, clear_cache, db):
    parsed = QuerySpec(protocol_hint="HTRF", metric="IC50")
    with patch("compounds.nlp.view.parse_prompt", return_value=parsed):
        resp = client.post(URL, {"prompt": "whatever"}, format="json")
    assert resp.status_code == 400
    assert resp.data["kind"] == "scope"


def test_spec_error_missing_metric_returns_400(client, feature_on, clear_cache, world):
    parsed = QuerySpec(
        registration_target_as_typed="AR degraders",
        protocol_hint="HTRF",
        # metric deliberately omitted
    )
    with patch("compounds.nlp.view.parse_prompt", return_value=parsed):
        resp = client.post(URL, {"prompt": "whatever"}, format="json")
    assert resp.status_code == 400
    assert resp.data["kind"] == "spec"
    assert resp.data["field"] == "metric"


def test_not_a_query_returns_200_with_reason(client, feature_on, clear_cache, db):
    with patch(
        "compounds.nlp.view.parse_prompt",
        return_value=NotAQuery(reason="Asked for a written summary, not a table."),
    ):
        resp = client.post(URL, {"prompt": "summarise the SAR"}, format="json")
    assert resp.status_code == 200
    assert resp.data["status"] == "not_a_query"
    assert "summary" in resp.data["reason"]


def test_parse_error_returns_502(client, feature_on, clear_cache, db):
    with patch(
        "compounds.nlp.view.parse_prompt",
        return_value=ParseError(message="LLM returned junk", raw="{not json"),
    ):
        resp = client.post(URL, {"prompt": "anything"}, format="json")
    assert resp.status_code == 502
    assert resp.data["kind"] == "parse"


# ---------------------------------------------------------------------------
# Daily cap
# ---------------------------------------------------------------------------


def test_daily_cap_blocks_after_limit(client, feature_on, clear_cache, world):
    parsed = QuerySpec(
        registration_target_as_typed="AR degraders",
        protocol_hint="HTRF",
        metric="IC50",
    )
    # Pre-seed the cache to the cap so the next request tips over.
    from compounds.nlp.view import _cap_key
    cache.set(_cap_key(client.handler._force_user.pk), DAILY_CAP, timeout=3600)

    with patch("compounds.nlp.view.parse_prompt", return_value=parsed):
        resp = client.post(URL, {"prompt": "whatever"}, format="json")
    assert resp.status_code == 429
    assert resp.data["kind"] == "rate_limited"


def test_separate_users_have_separate_caps(feature_on, clear_cache, world, db):
    # Alice at the cap; Bob should still be allowed.
    alice = User.objects.create_user(username="alice2", email="a2@example.org")
    bob = User.objects.create_user(username="bob2", email="b2@example.org")
    from compounds.nlp.view import _cap_key
    cache.set(_cap_key(alice.pk), DAILY_CAP, timeout=3600)

    parsed = QuerySpec(
        registration_target_as_typed="AR degraders",
        protocol_hint="HTRF",
        metric="IC50",
    )

    alice_c = APIClient(); alice_c.force_authenticate(user=alice)
    bob_c = APIClient(); bob_c.force_authenticate(user=bob)

    with patch("compounds.nlp.view.parse_prompt", return_value=parsed):
        assert alice_c.post(URL, {"prompt": "x"}, format="json").status_code == 429
        assert bob_c.post(URL, {"prompt": "x"}, format="json").status_code == 200
