"""Prompt → CompoundSelector via Azure OpenAI structured output (§8).

The LLM is a parser, not a query engine (§3):
- It sees a short static system prompt and the user's one-line query.
- It emits exactly one JSON object matching ``PROMPT_SCHEMA`` — either the
  CompoundSelector fields (with unspecified fields left null/empty) or
  ``{"not_a_query": true, "reason": "..."}`` for requests that aren't
  compound-selection queries.
- It never sees catalogs, compound data, assay results, or IDs.

**Pivoted 2026-04-23**: the LLM no longer decides display columns / KPIs
for output — it only describes **which compounds** to pick. Display is
owned downstream by the aggregation page. The schema collapses
accordingly: no columns / columns_explicit / output_metric; instead a
list of measurement_filters, each with its own protocol_hint / metric /
threshold, ANDed together.

For tests, patch ``compounds.nlp.azure_client.get_azure_openai_client``
to return a MagicMock whose ``chat.completions.create`` returns a shaped
response. See ``test_llm.py`` for the pattern.
"""

from __future__ import annotations

import json
from typing import Any, List

from .azure_client import get_azure_openai_client, get_model_name
from .spec import (
    CompoundSelector,
    MeasurementFilter,
    NotAQuery,
    ParseError,
    PromptParseResult,
    Threshold,
)


SYSTEM_PROMPT = """\
You translate a user's English request about registered compounds into a
CompoundSelector JSON object. You do not answer the question; you do not
fetch data; you only parse intent into structured fields.

The user is looking for a SET OF COMPOUNDS to view. Your job is to describe:
- WHICH compounds (by target registration, and optionally assay exposure)
- WHAT FILTERS those compounds must satisfy (measurement thresholds)

The downstream system intersects the filters to compute the matching
compound list, then redirects the user to a project-card view that
automatically renders spider plots from configured data relationships.
You DO NOT decide what data columns to show, what metric to tabulate on
the output axis, or what display format to use. You only describe the
selection.

Rules:
- Emit strings verbatim as the user typed them for
  registration_target_as_typed, assay_target_as_typed, protocol_hint,
  and metric. Do not canonicalise, translate case, or guess IDs.
- For single-target prompts like "mEGFR compounds ...", set
  registration_target_as_typed = "mEGFR" and leave assay_target_as_typed
  null — unless the user distinguishes WHERE compounds are tested from
  WHO designed them (e.g. "ARd compounds tested against KRAS" sets
  registration to "ARd" and assay to "KRAS").
- measurement_filters is a list. Multiple filters are ANDed — the
  compound must satisfy ALL of them. This is how cross-protocol
  selectivity is expressed:
    "mEGFR compounds with WT IC50 < 10 uM AND TM IC50 > 1 uM"
    → two filters, one per protocol, both on metric "IC50".
- Each filter has optional protocol_hint, metric, threshold. Null fields:
  - protocol_hint null → any protocol in scope is eligible
  - metric null → any numeric KPI counts
  - threshold null → just requires a valid measurement exists
- threshold is {op, value, unit}. op is one of <, <=, >, >=, =, !=.
  Echo the unit verbatim as typed ("uM" not "µM").
- Conversational filler nouns — "hits", "compounds", "molecules",
  "series", "analogues" — carry no structural meaning unless the user
  names an explicit cutoff. Do not synthesise a threshold from filler
  words alone.
- If the request asks for rows of specific column values, data columns,
  phys-chem properties, or display formats — ignore those preferences.
  Those live on the aggregation page the user lands on; your job ends
  at the compound selection.
- If the request is not a compound-selection query (summarise the SAR,
  register this compound, delete something, ask a meta-question), emit
  {"not_a_query": true, "reason": "<short explanation>"}.

Never emit compound IDs, numeric result values, counts, or catalog
entries. Only the selector shape.
"""


# JSON Schema for the response. Every property listed in `required` per
# Azure OpenAI strict-mode rules. `additionalProperties: false` throughout.
PROMPT_SCHEMA: dict = {
    "type": "object",
    "additionalProperties": False,
    "properties": {
        "not_a_query": {"type": ["boolean", "null"]},
        "reason": {"type": ["string", "null"]},
        "registration_target_as_typed": {"type": ["string", "null"]},
        "assay_target_as_typed": {"type": ["string", "null"]},
        "measurement_filters": {
            "type": "array",
            "items": {
                "type": "object",
                "additionalProperties": False,
                "properties": {
                    "protocol_hint": {"type": ["string", "null"]},
                    "metric": {"type": ["string", "null"]},
                    "threshold": {
                        "type": ["object", "null"],
                        "additionalProperties": False,
                        "properties": {
                            "op": {"type": "string"},
                            "value": {"type": "number"},
                            "unit": {"type": ["string", "null"]},
                        },
                        "required": ["op", "value", "unit"],
                    },
                },
                "required": ["protocol_hint", "metric", "threshold"],
            },
        },
    },
    "required": [
        "not_a_query",
        "reason",
        "registration_target_as_typed",
        "assay_target_as_typed",
        "measurement_filters",
    ],
}


def _threshold_from_dict(data: Any) -> Threshold:
    return Threshold(
        op=data["op"],
        value=float(data["value"]),
        unit=data.get("unit"),
    )


def _filter_from_dict(data: Any) -> MeasurementFilter:
    if not isinstance(data, dict):
        raise ValueError(f"filter not a dict: {data!r}")
    threshold_raw = data.get("threshold")
    if threshold_raw is None:
        threshold = None
    else:
        if not isinstance(threshold_raw, dict):
            raise ValueError(f"threshold not a dict: {threshold_raw!r}")
        if threshold_raw.get("op") is None or threshold_raw.get("value") is None:
            raise ValueError("threshold missing op or value")
        threshold = _threshold_from_dict(threshold_raw)
    return MeasurementFilter(
        protocol_hint=data.get("protocol_hint"),
        metric=data.get("metric"),
        threshold=threshold,
    )


def _to_parse_result(data: dict) -> PromptParseResult:
    """Map a validated LLM JSON object to a CompoundSelector / NotAQuery /
    ParseError."""
    if data.get("not_a_query"):
        reason = data.get("reason") or "Not a compound-selection query."
        return NotAQuery(reason=reason)

    raw_filters = data.get("measurement_filters") or []
    if not isinstance(raw_filters, list):
        return ParseError(
            message="measurement_filters must be a list",
            raw=json.dumps(data),
        )

    try:
        filters = [_filter_from_dict(f) for f in raw_filters]
    except ValueError as e:
        return ParseError(
            message=f"Malformed measurement_filter: {e}",
            raw=json.dumps(raw_filters),
        )

    return CompoundSelector(
        registration_target_as_typed=data.get("registration_target_as_typed"),
        assay_target_as_typed=data.get("assay_target_as_typed"),
        measurement_filters=filters,
    )


def parse_prompt(prompt: str) -> PromptParseResult:
    """Send ``prompt`` to Azure OpenAI in structured-output mode and parse
    the response. Returns a CompoundSelector, NotAQuery, or ParseError."""
    if not prompt or not prompt.strip():
        return ParseError(message="Empty prompt.")

    client = get_azure_openai_client()
    response = client.chat.completions.create(
        model=get_model_name(),
        response_format={
            "type": "json_schema",
            "json_schema": {
                "name": "CompoundSelector",
                "schema": PROMPT_SCHEMA,
                "strict": True,
            },
        },
        messages=[
            {"role": "system", "content": SYSTEM_PROMPT},
            {"role": "user", "content": prompt},
        ],
        temperature=0,
    )
    content = _extract_content(response)
    if content is None:
        return ParseError(message="LLM response contained no content.")

    try:
        data = json.loads(content)
    except json.JSONDecodeError as e:
        return ParseError(message=f"LLM output was not valid JSON: {e}", raw=content)

    if not isinstance(data, dict):
        return ParseError(message="LLM output was not a JSON object.", raw=content)

    return _to_parse_result(data)


def _extract_content(response: Any) -> str | None:
    """Pull the message content out of an OpenAI SDK response, tolerant of
    both dict-shaped and object-shaped responses so mocks can be either."""
    try:
        choice = response.choices[0]
    except (AttributeError, IndexError, TypeError):
        return None
    message = getattr(choice, "message", None)
    if message is None:
        return None
    return getattr(message, "content", None)
