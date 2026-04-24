"""LLM parse tests.

The real Azure OpenAI client is never instantiated — tests patch
`compounds.nlp.azure_client.get_azure_openai_client` (and `get_model_name`)
to inject a MagicMock that plays back a hand-written JSON payload. This
verifies our glue code around the SDK without touching the network and
without requiring `openai` to be pip-installed locally.
"""

from __future__ import annotations

import json
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest

from compounds.nlp.llm import (
    PROMPT_SCHEMA,
    SYSTEM_PROMPT,
    _to_parse_result,
    parse_prompt,
)
from compounds.nlp.spec import (
    CompoundSelector,
    MeasurementFilter,
    NotAQuery,
    ParseError,
    Threshold,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _empty_payload() -> dict:
    """Valid schema-conforming payload with every field null / empty
    (strict-mode requires every property present even if null)."""
    return {
        "not_a_query": None,
        "reason": None,
        "registration_target_as_typed": None,
        "assay_target_as_typed": None,
        "measurement_filters": [],
        "assay_selector": None,
    }


def _mock_client(payload: dict | str) -> MagicMock:
    content = payload if isinstance(payload, str) else json.dumps(payload)
    response = SimpleNamespace(
        choices=[SimpleNamespace(message=SimpleNamespace(content=content))]
    )
    client = MagicMock()
    client.chat.completions.create.return_value = response
    return client


def _patched_parse(prompt: str, payload: dict | str):
    client = _mock_client(payload)
    with patch("compounds.nlp.llm.get_azure_openai_client", return_value=client), \
         patch("compounds.nlp.llm.get_model_name", return_value="gpt-4o"):
        result = parse_prompt(prompt)
    return result, client


# ---------------------------------------------------------------------------
# _to_parse_result — pure translator, no mocking
# ---------------------------------------------------------------------------


def test_single_filter_with_threshold():
    payload = _empty_payload()
    payload.update({
        "registration_target_as_typed": "ARd",
        "measurement_filters": [
            {"protocol_hint": "HTRF", "metric": "IC50",
             "threshold": {"op": "<", "value": 10, "unit": "uM"}},
        ],
    })
    result = _to_parse_result(payload)
    assert isinstance(result, CompoundSelector)
    assert result.registration_target_as_typed == "ARd"
    assert result.assay_target_as_typed is None
    assert len(result.measurement_filters) == 1
    flt = result.measurement_filters[0]
    assert flt.protocol_hint == "HTRF"
    assert flt.metric == "IC50"
    assert flt.threshold == Threshold(op="<", value=10.0, unit="uM")


def test_cross_protocol_selectivity():
    """The headline v2 capability — multiple filters AND'd together."""
    payload = _empty_payload()
    payload.update({
        "registration_target_as_typed": "mEGFR",
        "measurement_filters": [
            {"protocol_hint": "mEGFR TR-FRET WT", "metric": "IC50",
             "threshold": {"op": "<", "value": 10, "unit": "uM"}},
            {"protocol_hint": "TM", "metric": "IC50",
             "threshold": {"op": ">", "value": 1, "unit": "uM"}},
        ],
    })
    result = _to_parse_result(payload)
    assert isinstance(result, CompoundSelector)
    assert len(result.measurement_filters) == 2
    assert result.measurement_filters[0].threshold.op == "<"
    assert result.measurement_filters[1].threshold.op == ">"


def test_no_filters_just_target():
    payload = _empty_payload()
    payload["registration_target_as_typed"] = "ARd"
    result = _to_parse_result(payload)
    assert isinstance(result, CompoundSelector)
    assert result.measurement_filters == []


def test_filter_without_threshold():
    payload = _empty_payload()
    payload.update({
        "registration_target_as_typed": "ARd",
        "measurement_filters": [
            {"protocol_hint": "HTRF", "metric": None, "threshold": None},
        ],
    })
    result = _to_parse_result(payload)
    assert isinstance(result, CompoundSelector)
    assert result.measurement_filters[0].threshold is None


def test_cross_project_both_targets():
    payload = _empty_payload()
    payload.update({
        "registration_target_as_typed": "ARd",
        "assay_target_as_typed": "AKT",
        "measurement_filters": [
            {"protocol_hint": None, "metric": "EC50",
             "threshold": {"op": "<", "value": 1, "unit": "uM"}},
        ],
    })
    result = _to_parse_result(payload)
    assert isinstance(result, CompoundSelector)
    assert result.registration_target_as_typed == "ARd"
    assert result.assay_target_as_typed == "AKT"


def test_unitless_threshold():
    payload = _empty_payload()
    payload.update({
        "registration_target_as_typed": "Caco-2",
        "measurement_filters": [
            {"protocol_hint": "Caco-2 permeability", "metric": "efflux_ratio",
             "threshold": {"op": ">", "value": 2.0, "unit": None}},
        ],
    })
    result = _to_parse_result(payload)
    assert isinstance(result, CompoundSelector)
    assert result.measurement_filters[0].threshold.unit is None


def test_not_a_query():
    payload = _empty_payload()
    payload["not_a_query"] = True
    payload["reason"] = "Asking for a written SAR summary, not a selection."
    result = _to_parse_result(payload)
    assert isinstance(result, NotAQuery)
    assert "SAR summary" in result.reason


def test_not_a_query_blank_reason_defaults():
    payload = _empty_payload()
    payload["not_a_query"] = True
    payload["reason"] = ""
    result = _to_parse_result(payload)
    assert isinstance(result, NotAQuery)
    assert result.reason == "Not a compound-selection query."


def test_malformed_threshold_returns_parse_error():
    payload = _empty_payload()
    payload.update({
        "registration_target_as_typed": "ARd",
        "measurement_filters": [
            {"protocol_hint": "HTRF", "metric": "IC50",
             "threshold": {"op": "<", "value": None, "unit": "nM"}},
        ],
    })
    result = _to_parse_result(payload)
    assert isinstance(result, ParseError)
    assert "op or value" in result.message


def test_measurement_filters_not_a_list_returns_parse_error():
    payload = _empty_payload()
    payload["measurement_filters"] = "oops"
    result = _to_parse_result(payload)
    assert isinstance(result, ParseError)


# ---------------------------------------------------------------------------
# parse_prompt — full path with mocked client
# ---------------------------------------------------------------------------


def test_parse_prompt_happy_path():
    payload = _empty_payload()
    payload.update({
        "registration_target_as_typed": "ARd",
        "measurement_filters": [
            {"protocol_hint": "HTRF", "metric": "IC50",
             "threshold": {"op": "<", "value": 10, "unit": "uM"}},
        ],
    })
    result, client = _patched_parse(
        "ARd compounds with HTRF IC50 < 10 uM",
        payload,
    )
    assert isinstance(result, CompoundSelector)

    client.chat.completions.create.assert_called_once()
    kwargs = client.chat.completions.create.call_args.kwargs
    assert kwargs["model"] == "gpt-4o"
    assert kwargs["temperature"] == 0
    assert kwargs["response_format"]["type"] == "json_schema"
    assert kwargs["response_format"]["json_schema"]["schema"] is PROMPT_SCHEMA
    assert kwargs["response_format"]["json_schema"]["strict"] is True
    messages = kwargs["messages"]
    assert messages[0]["role"] == "system"
    assert messages[0]["content"] == SYSTEM_PROMPT
    assert messages[1]["role"] == "user"
    assert "ARd" in messages[1]["content"]


def test_parse_prompt_not_a_query():
    payload = _empty_payload()
    payload["not_a_query"] = True
    payload["reason"] = "Request is for a written summary, not a selection."
    result, _ = _patched_parse("Summarise the SAR of my ARd series", payload)
    assert isinstance(result, NotAQuery)
    assert "summary" in result.reason.lower()


def test_parse_prompt_malformed_json_returns_parse_error():
    result, _ = _patched_parse("anything", "{not valid json")
    assert isinstance(result, ParseError)
    assert "not valid JSON" in result.message
    assert result.raw == "{not valid json"


def test_parse_prompt_json_not_object_returns_parse_error():
    result, _ = _patched_parse("anything", "[1, 2, 3]")
    assert isinstance(result, ParseError)
    assert "not a JSON object" in result.message


def test_parse_prompt_missing_content_returns_parse_error():
    empty_response = SimpleNamespace(choices=[])
    client = MagicMock()
    client.chat.completions.create.return_value = empty_response
    with patch("compounds.nlp.llm.get_azure_openai_client", return_value=client), \
         patch("compounds.nlp.llm.get_model_name", return_value="gpt-4o"):
        result = parse_prompt("something")
    assert isinstance(result, ParseError)
    assert "no content" in result.message


def test_parse_prompt_empty_input_returns_parse_error():
    with patch("compounds.nlp.llm.get_azure_openai_client") as factory:
        result = parse_prompt("   ")
    assert isinstance(result, ParseError)
    factory.assert_not_called()


# ---------------------------------------------------------------------------
# Schema-shape guards
# ---------------------------------------------------------------------------


def test_schema_declares_every_top_level_property_as_required():
    assert set(PROMPT_SCHEMA["required"]) == set(PROMPT_SCHEMA["properties"].keys())


def test_assay_selector_subschema_strictly_shaped():
    assay = PROMPT_SCHEMA["properties"]["assay_selector"]
    # Nullable object — LLM populates when the user asks about assays,
    # leaves null otherwise.
    assert assay["type"] == ["object", "null"]
    assert assay["additionalProperties"] is False
    assert set(assay["required"]) == {
        "target_as_typed", "protocol_hint", "date_range", "created_by_as_typed",
    }


def test_assay_selector_dispatch():
    """When assay_selector is populated, _to_parse_result returns an
    AssaySelector (not a CompoundSelector)."""
    from compounds.nlp.spec import AssaySelector
    payload = _empty_payload()
    payload["assay_selector"] = {
        "target_as_typed": "CDK4",
        "protocol_hint": None,
        "date_range": {"after": "2026-04-17", "before": "2026-04-24"},
        "created_by_as_typed": None,
    }
    result = _to_parse_result(payload)
    assert isinstance(result, AssaySelector)
    assert result.target_as_typed == "CDK4"
    assert result.date_range is not None
    assert result.date_range.after == "2026-04-17"


def test_compound_selector_still_dispatches_when_assay_selector_null():
    """A top-level compound-selector prompt must still produce a
    CompoundSelector even now that assay_selector is a sibling field."""
    payload = _empty_payload()
    payload["registration_target_as_typed"] = "CDK4"
    result = _to_parse_result(payload)
    assert isinstance(result, CompoundSelector)
    assert result.registration_target_as_typed == "CDK4"


def test_measurement_filters_items_strictly_shaped():
    filters = PROMPT_SCHEMA["properties"]["measurement_filters"]
    assert filters["type"] == "array"
    item = filters["items"]
    assert item["additionalProperties"] is False
    assert set(item["required"]) == {
        "protocol_hint", "metric", "threshold",
        "assay_date_range", "assayed_by_as_typed",
    }


def test_threshold_subschema_additional_properties_disabled():
    threshold = PROMPT_SCHEMA["properties"]["measurement_filters"]["items"]["properties"]["threshold"]
    assert threshold["additionalProperties"] is False
    assert set(threshold["required"]) == {"op", "value", "unit"}


def test_today_is_prepended_to_user_message():
    """The LLM needs today's date to resolve relative phrasings like
    "last 30 days". Injected in the user message (not the system prompt)
    so the system prompt stays static and cache-friendly."""
    from compounds.nlp.llm import _wrap_with_today

    wrapped = _wrap_with_today("CDK4 compounds registered last month", today="2026-04-24")
    assert wrapped.startswith("[Today: 2026-04-24]\n")
    assert "CDK4 compounds registered last month" in wrapped


def test_parse_prompt_wraps_with_today_in_user_message():
    """The user message the SDK sees must carry the [Today: ...] prefix."""
    payload = _empty_payload()
    payload["registration_target_as_typed"] = "CDK4"
    _, client = _patched_parse("CDK4 compounds", payload)
    messages = client.chat.completions.create.call_args.kwargs["messages"]
    user_content = messages[1]["content"]
    assert user_content.startswith("[Today: ")
    assert "CDK4 compounds" in user_content


def test_system_prompt_mentions_key_discipline_rules():
    assert "verbatim" in SYSTEM_PROMPT
    assert "not_a_query" in SYSTEM_PROMPT
    assert "measurement_filters" in SYSTEM_PROMPT
    # Emphasise that the LLM's job ends at selection.
    assert "selection" in SYSTEM_PROMPT.lower()
    # Filler-noun discipline (§19.6)
    assert "filler" in SYSTEM_PROMPT or "hits" in SYSTEM_PROMPT
