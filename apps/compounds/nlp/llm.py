"""Prompt → QuerySpec via Azure OpenAI structured output (§8).

The LLM is a parser, not a query engine (§3):
- It sees a short static system prompt and the user's one-line query.
- It emits exactly one JSON object matching ``PROMPT_SCHEMA`` — either the
  QuerySpec fields (with any unspecified field left null) or
  ``{"not_a_query": true, "reason": "..."}`` for requests that aren't
  tabular queries at all.
- It never sees catalogs, compound data, assay results, or IDs.

Everything downstream — target/protocol resolution, threshold conversion,
row evaluation, ORM walk — is deterministic and LLM-free.

For tests, patch ``compounds.nlp.azure_client.get_azure_openai_client`` to
return a MagicMock whose ``chat.completions.create`` returns a shaped
response. See ``test_llm.py`` for the pattern.
"""

from __future__ import annotations

import json
from typing import Any

from .azure_client import get_azure_openai_client, get_model_name
from .spec import NotAQuery, ParseError, PromptParseResult, QuerySpec, Threshold


SYSTEM_PROMPT = """\
You translate a user's English request about registered compounds and assay
results into a QuerySpec JSON object. You do not answer the question; you do
not fetch data; you only parse intent into structured fields.

Rules:
- Emit strings verbatim as the user typed them for `registration_target_as_typed`,
  `assay_target_as_typed`, and `protocol_hint`. Do not canonicalise, translate,
  or guess IDs.
- For single-target prompts like "compounds in the ARd project ...", set both
  target fields to the same typed string. For cross-project prompts like
  "ARd compounds tested against AKT", set registration to "ARd" and assay to
  "AKT". If the user names only an assay context ("compounds tested against X"),
  leave `registration_target_as_typed` null.
- `metric` is the KPI name: IC50, EC50, Ki, Kd, pIC50, % inhibition, DC50,
  CLint, etc. Emit verbatim as typed.
- `threshold` is `{op, value, unit}`. `op` is one of `<`, `<=`, `>`, `>=`,
  `=`, `!=`. Echo the unit verbatim as typed (e.g. "uM" not "µM", "nM",
  "%", etc.). If no threshold is in the prompt, leave `threshold` null.
- `columns` defaults to `"phys_chem"` (all eight physicochemical properties)
  unless the user asks for Lipinski-only (`"lipinski"`) or names fields
  explicitly in `columns_explicit`.
- Conversational filler nouns like "hits", "compounds", "molecules", "series"
  carry no structural meaning unless the user also names a cutoff. Do not
  synthesise a threshold from filler words alone.
- If the request is not a tabular query of compounds/assay results, emit
  `{"not_a_query": true, "reason": "<short explanation>"}`. Examples:
  "summarise the SAR", "register this compound", "delete X" — all not queries.

Never emit compound IDs, numeric result values, counts, or catalog entries.
Only the spec.
"""


# JSON Schema for the response. Designed for Azure OpenAI strict-mode
# structured output: every property is listed in `required` and nullable
# values use union-with-null `type` lists. `additionalProperties: false`.
PROMPT_SCHEMA: dict = {
    "type": "object",
    "additionalProperties": False,
    "properties": {
        "not_a_query": {"type": ["boolean", "null"]},
        "reason": {"type": ["string", "null"]},
        "registration_target_as_typed": {"type": ["string", "null"]},
        "assay_target_as_typed": {"type": ["string", "null"]},
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
        "columns": {"type": ["string", "null"]},
        "columns_explicit": {
            "type": ["array", "null"],
            "items": {"type": "string"},
        },
    },
    "required": [
        "not_a_query",
        "reason",
        "registration_target_as_typed",
        "assay_target_as_typed",
        "protocol_hint",
        "metric",
        "threshold",
        "columns",
        "columns_explicit",
    ],
}


def _to_parse_result(data: dict) -> PromptParseResult:
    """Map a validated LLM JSON object to a QuerySpec / NotAQuery / ParseError."""
    if data.get("not_a_query"):
        reason = data.get("reason") or "Not a tabular query."
        return NotAQuery(reason=reason)

    threshold_raw = data.get("threshold")
    threshold = None
    if threshold_raw is not None:
        op = threshold_raw.get("op")
        value = threshold_raw.get("value")
        if op is None or value is None:
            return ParseError(
                message="Threshold object missing op or value.",
                raw=json.dumps(threshold_raw),
            )
        threshold = Threshold(
            op=op,
            value=float(value),
            unit=threshold_raw.get("unit"),
        )

    return QuerySpec(
        registration_target_as_typed=data.get("registration_target_as_typed"),
        assay_target_as_typed=data.get("assay_target_as_typed"),
        protocol_hint=data.get("protocol_hint"),
        metric=data.get("metric"),
        threshold=threshold,
        columns=data.get("columns"),
        columns_explicit=data.get("columns_explicit"),
    )


def parse_prompt(prompt: str) -> PromptParseResult:
    """Send ``prompt`` to Azure OpenAI in structured-output mode and parse the
    response. Returns a QuerySpec, NotAQuery, or ParseError."""
    if not prompt or not prompt.strip():
        return ParseError(message="Empty prompt.")

    client = get_azure_openai_client()
    response = client.chat.completions.create(
        model=get_model_name(),
        response_format={
            "type": "json_schema",
            "json_schema": {
                "name": "QuerySpec",
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
