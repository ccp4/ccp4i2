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
from compounds.nlp.spec import NotAQuery, ParseError, QuerySpec, Threshold


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _empty_payload() -> dict:
    """Valid schema-conforming payload with every field null (for
    strict-mode shapes). Individual tests override the fields they care
    about."""
    return {
        "not_a_query": None,
        "reason": None,
        "registration_target_as_typed": None,
        "assay_target_as_typed": None,
        "protocol_hint": None,
        "metric": None,
        "threshold": None,
        "columns": None,
        "columns_explicit": None,
    }


def _mock_client(payload: dict | str) -> MagicMock:
    """Build a MagicMock that returns ``payload`` as the LLM's JSON content."""
    content = payload if isinstance(payload, str) else json.dumps(payload)
    response = SimpleNamespace(
        choices=[SimpleNamespace(message=SimpleNamespace(content=content))]
    )
    client = MagicMock()
    client.chat.completions.create.return_value = response
    return client


def _patched_parse(prompt: str, payload: dict | str):
    """Run parse_prompt with the factory patched to return a mock client."""
    client = _mock_client(payload)
    with patch("compounds.nlp.llm.get_azure_openai_client", return_value=client), \
         patch("compounds.nlp.llm.get_model_name", return_value="gpt-4o"):
        result = parse_prompt(prompt)
    return result, client


# ---------------------------------------------------------------------------
# _to_parse_result — pure translator, no mocking needed
# ---------------------------------------------------------------------------


def test_to_parse_result_single_target_with_threshold():
    payload = _empty_payload()
    payload.update({
        "registration_target_as_typed": "ARd",
        "assay_target_as_typed": "ARd",
        "protocol_hint": "HTRF",
        "metric": "IC50",
        "threshold": {"op": "<", "value": 10, "unit": "uM"},
        "columns": "phys_chem",
    })
    result = _to_parse_result(payload)
    assert isinstance(result, QuerySpec)
    assert result.registration_target_as_typed == "ARd"
    assert result.assay_target_as_typed == "ARd"
    assert result.protocol_hint == "HTRF"
    assert result.metric == "IC50"
    assert result.threshold == Threshold(op="<", value=10.0, unit="uM")
    assert result.columns == "phys_chem"


def test_to_parse_result_assay_only_no_threshold():
    payload = _empty_payload()
    payload.update({
        "assay_target_as_typed": "AKT",
        "protocol_hint": "HTRF",
        "metric": "IC50",
    })
    result = _to_parse_result(payload)
    assert isinstance(result, QuerySpec)
    assert result.registration_target_as_typed is None
    assert result.assay_target_as_typed == "AKT"
    assert result.threshold is None


def test_to_parse_result_cross_project():
    payload = _empty_payload()
    payload.update({
        "registration_target_as_typed": "ARd",
        "assay_target_as_typed": "AKT",
        "protocol_hint": "HTRF",
        "metric": "IC50",
        "threshold": {"op": "<=", "value": 5.0, "unit": "nM"},
    })
    result = _to_parse_result(payload)
    assert isinstance(result, QuerySpec)
    assert result.registration_target_as_typed == "ARd"
    assert result.assay_target_as_typed == "AKT"
    assert result.threshold.op == "<="
    assert result.threshold.value == 5.0
    assert result.threshold.unit == "nM"


def test_to_parse_result_unitless_threshold():
    payload = _empty_payload()
    payload.update({
        "registration_target_as_typed": "Caco-2",
        "assay_target_as_typed": "Caco-2",
        "protocol_hint": "Caco-2 efflux",
        "metric": "efflux_ratio",
        "threshold": {"op": ">", "value": 2.0, "unit": None},
    })
    result = _to_parse_result(payload)
    assert isinstance(result, QuerySpec)
    assert result.threshold.unit is None


def test_to_parse_result_not_a_query():
    payload = _empty_payload()
    payload["not_a_query"] = True
    payload["reason"] = "The user is asking for a written SAR summary, not a table."
    result = _to_parse_result(payload)
    assert isinstance(result, NotAQuery)
    assert "SAR summary" in result.reason


def test_to_parse_result_not_a_query_with_blank_reason_defaults():
    payload = _empty_payload()
    payload["not_a_query"] = True
    payload["reason"] = ""
    result = _to_parse_result(payload)
    assert isinstance(result, NotAQuery)
    assert result.reason == "Not a tabular query."


def test_to_parse_result_malformed_threshold_returns_parse_error():
    payload = _empty_payload()
    payload["threshold"] = {"op": "<", "value": None, "unit": "nM"}
    result = _to_parse_result(payload)
    assert isinstance(result, ParseError)
    assert "op or value" in result.message


def test_to_parse_result_columns_explicit_list_preserved():
    payload = _empty_payload()
    payload.update({
        "registration_target_as_typed": "ARd",
        "assay_target_as_typed": "ARd",
        "protocol_hint": "HTRF",
        "metric": "IC50",
        "columns_explicit": ["molecular_weight", "clogp"],
    })
    result = _to_parse_result(payload)
    assert isinstance(result, QuerySpec)
    assert result.columns_explicit == ["molecular_weight", "clogp"]


# ---------------------------------------------------------------------------
# parse_prompt — full path with a mocked client
# ---------------------------------------------------------------------------


def test_parse_prompt_happy_path():
    payload = _empty_payload()
    payload.update({
        "registration_target_as_typed": "ARd",
        "assay_target_as_typed": "ARd",
        "protocol_hint": "HTRF",
        "metric": "IC50",
        "threshold": {"op": "<", "value": 10, "unit": "uM"},
        "columns": "phys_chem",
    })
    result, client = _patched_parse(
        "phys chem for ARd compounds with HTRF IC50 < 10 uM",
        payload,
    )
    assert isinstance(result, QuerySpec)

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
    payload["reason"] = "Request is for a written summary, not a table."
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
    # Empty / whitespace prompts should short-circuit before hitting the SDK.
    with patch("compounds.nlp.llm.get_azure_openai_client") as factory:
        result = parse_prompt("   ")
    assert isinstance(result, ParseError)
    factory.assert_not_called()


# ---------------------------------------------------------------------------
# Schema shape sanity — guards against accidental drift
# ---------------------------------------------------------------------------


def test_schema_declares_every_spec_field_as_required():
    # Strict mode demands every property appears in `required`, even nullable
    # ones. If we add a field to the schema, this test catches the missing
    # `required` entry before Azure does.
    assert set(PROMPT_SCHEMA["required"]) == set(PROMPT_SCHEMA["properties"].keys())


def test_threshold_subschema_additional_properties_disabled():
    threshold = PROMPT_SCHEMA["properties"]["threshold"]
    assert threshold["additionalProperties"] is False
    assert set(threshold["required"]) == {"op", "value", "unit"}


def test_system_prompt_mentions_key_discipline_rules():
    # Spot-check that the prompt still carries the load-bearing constraints.
    # If the prompt is rewritten these assertions will need updating; that's
    # intentional — a bare rewrite without checking these rules is a regression.
    assert "verbatim" in SYSTEM_PROMPT
    assert "not_a_query" in SYSTEM_PROMPT
    assert "IC50" in SYSTEM_PROMPT
