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

import datetime
import json
from typing import Any, List, Optional

from .azure_client import get_azure_openai_client, get_model_name
from .spec import (
    AssaySelector,
    CompoundSelector,
    DateRange,
    MeasurementFilter,
    NotAQuery,
    ParseError,
    PromptParseResult,
    Threshold,
)
from .substructures import all_scaffold_names


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
- If the user names NO target at all ("all pyrimidines", "compounds
  registered by Alice in 2025", "compounds with HTRF IC50 < 10 nM")
  then leave both target fields null. The backend runs the selection
  across the whole registry, narrowed only by the other predicates
  you emit. Do NOT invent a target to satisfy a schema nicety; null is
  the correct answer here.
- measurement_filters is a list. Multiple filters are ANDed — the
  compound must satisfy ALL of them. This is how cross-protocol
  selectivity is expressed:
    "mEGFR compounds with WT IC50 < 10 uM AND TM IC50 > 1 uM"
    → two filters, one per protocol, both on metric "IC50".
- Each filter has optional protocol_hint, metric, threshold, assay_date_range.
  Null fields:
  - protocol_hint null → any protocol in scope is eligible
  - metric null → any numeric KPI counts
  - threshold null → just requires a valid measurement exists
  - assay_date_range null → no date constraint on the assay
- threshold is {op, value, unit}. op is one of <, <=, >, >=, =, !=.
  Echo the unit verbatim as typed ("uM" not "µM").
- Substructures / scaffolds / functional groups: when the user names
  a chemical substructure — "compounds containing pyrimidine", "with a
  piperidine and a sulfonamide", "tetrazole-containing ARd compounds"
  — emit the names in `scaffold_hints` (a list of strings, ANDed).
  Emit the name the user typed, verbatim (singular or plural as
  they wrote it); the backend resolves the name against its curated
  SMARTS catalog. You MUST NOT emit SMARTS yourself, guess a ring
  structure, or synthesise a substructure when the user didn't name
  one. If the user names multiple substructures ("a pyrimidine and
  a piperidine"), each goes in `scaffold_hints` as its own entry.
  The same applies inside `assay_selector` — "HTRF assays on
  pyrimidines" sets `assay_selector.scaffold_hints = ["pyrimidines"]`.
- People AND suppliers: `registered_by_as_typed` (on the selector)
  filters on WHO registered / made / synthesised / SUPPLIED the
  compound. The backend resolves it against both internal user
  accounts AND the chemical-supplier catalogue, so chemist names
  ("Alice Jones") and vendor names ("Enamine", "Sigma", "ChemBridge")
  are equally valid. `assayed_by_as_typed` (per filter) is User-only
  and filters on WHO ran the assay — use for "tested by", "assayed
  by", "run by". Emit the name / vendor / username / email VERBATIM
  as typed; the backend resolves. Examples:
    "compounds registered by Alice"            → registered_by: "Alice"
    "compounds made by Alice Jones"            → registered_by: "Alice Jones"
    "compounds from Enamine"                   → registered_by: "Enamine"
    "compounds supplied by Sigma"              → registered_by: "Sigma"
    "ChemBridge compounds"                     → registered_by: "ChemBridge" (when ChemBridge is a known supplier)
    "IC50 in HTRF tested by alice.jones"       → filter.assayed_by: "alice.jones"
    "assays run by alice.jones@ncl.ac.uk"      → filter.assayed_by: "alice.jones@ncl.ac.uk"
  If the user names themselves via "I", "me", "my" without naming a
  person, leave the field null — the backend cannot resolve first-person.
- Dates: both `registered_date_range` (on the selector) and
  `assay_date_range` (on each filter) are {after, before} half-open
  ISO-date ranges — `after` is the inclusive lower bound, `before`
  is the EXCLUSIVE upper bound. Map calendar units exactly:
    "in 2025"         → {after: "2025-01-01", before: "2026-01-01"}
    "in Q1 2026"      → {after: "2026-01-01", before: "2026-04-01"}
    "in March 2026"   → {after: "2026-03-01", before: "2026-04-01"}
    "since 2025-03-15"→ {after: "2025-03-15", before: null}
    "before 2024"     → {after: null, before: "2024-01-01"}
  Relative phrasings — "last 30 days", "this week", "recently" — resolve
  against the date given on the `[Today: YYYY-MM-DD]` line at the START
  of each user message. For example, with Today: 2026-04-24:
    "in the last 30 days" → {after: "2026-03-25", before: null}
    "this quarter"        → {after: "2026-04-01", before: "2026-07-01"}
    "yesterday"           → {after: "2026-04-23", before: "2026-04-24"}
  Emit null dates rather than guessing when the phrasing is ambiguous
  (e.g. "recently" without a time anchor).
- Ranking ("best X compounds" / "top 20 X compounds"): when the user
  asks for the *best* compounds, set `rank_by = "scorecard"` and
  `rank_top_n` to the number they named (or null when they didn't —
  the backend defaults to 20). The scorecard is the project-defined
  merit metric, so "best" only makes sense when the prompt names a
  target. Examples:
    "best CDK4 compounds"                  → rank_by: "scorecard", rank_top_n: null
    "top 10 ARd compounds"                 → rank_by: "scorecard", rank_top_n: 10
    "show me the top 50 EGFR compounds"    → rank_by: "scorecard", rank_top_n: 50
    "best ARd compounds with HTRF IC50 < 10 nM"
      → rank_by: "scorecard", rank_top_n: null
        registration_target_as_typed: "ARd"
        measurement_filters: [{ protocol_hint: "HTRF", metric: "IC50",
                                threshold: { op: "<", value: 10, unit: "nM" }}]
  When the user asks for ranking by something OTHER than scorecard
  ("best by potency" / "highest IC50") — leave rank_by null in v1;
  the chemist can sort manually on the aggregation page.
- Compound-ID pinning: when the user names specific compounds by ID —
  "compound 26007", "NCL-00026007", "ncl26007", "the lead compound
  NCL-00026007", "26007 and 26012" — emit each typed reference into
  `compound_refs_as_typed` as a single string. These are ADDITIVE
  (pinned alongside whatever the other predicates select), so they
  combine with filters via UNION rather than AND. Examples:
    "compound 26007 and compounds made since January"
      → compound_refs_as_typed: ["26007"]
        registered_date_range: {after: "<computed>", before: null}
    "show NCL-00026007 alongside ARd pyrimidines"
      → compound_refs_as_typed: ["NCL-00026007"]
        registration_target_as_typed: "ARd"
        scaffold_hints: ["pyrimidines"]
    "just compound 12345"
      → compound_refs_as_typed: ["12345"]
        (no other predicates)
  Emit the typed reference verbatim — full prefix or bare integer is
  fine. The backend's parser handles every variant. Do NOT emit a
  reference for vague phrases like "my compounds" or "the lead series"
  unless the user actually says an ID.
- Conversational filler nouns — "hits", "compounds", "molecules",
  "series", "analogues" — carry no structural meaning unless the user
  names an explicit cutoff. Do not synthesise a threshold from filler
  words alone.
- If the request asks for rows of specific column values, data columns,
  phys-chem properties, or display formats — ignore those preferences.
  Those live on the aggregation page the user lands on; your job ends
  at the compound selection.
- The two query families: by default the user is asking about
  COMPOUNDS ("mEGFR compounds with HTRF IC50 < 10 uM") — populate the
  top-level compound-selector fields and leave `assay_selector` null.
  If the user is instead asking about a list of ASSAY EXPERIMENTS
  ("CDK4 assays carried out last week", "HTRF assays conducted by
  Alice", "assays against BRD4 in Q1 2026") — populate the nested
  `assay_selector` object and leave ALL the top-level compound-selector
  fields null (registration_target_as_typed, assay_target_as_typed,
  registered_date_range, registered_by_as_typed, measurement_filters=[]).
  The grammatical subject of the prompt is the tell: "compounds" →
  compound selector; "assays" / "experiments" / "runs" → assay selector.
- Inside `assay_selector`, the fields mirror the compound filters but
  single-valued (one target, one protocol, one creator, one date
  range). "tested/conducted/carried out by X" → created_by_as_typed.
- If the request is not any kind of selection query (summarise the
  SAR, register this compound, delete something, ask a meta-question),
  emit {"not_a_query": true, "reason": "<short explanation>"}.

Never emit compound IDs, numeric result values, counts, or catalog
entries. Only the selector shape.
"""


# JSON Schema for the response. Every property listed in `required` per
# Azure OpenAI strict-mode rules. `additionalProperties: false` throughout.
# Reusable date-range subschema shared by `registered_date_range` (on the
# selector) and `assay_date_range` (on each measurement filter). Dates are
# "YYYY-MM-DD" strings; both ends nullable (half-open ranges).
_DATE_RANGE_SUBSCHEMA: dict = {
    "type": ["object", "null"],
    "additionalProperties": False,
    "properties": {
        "after": {"type": ["string", "null"]},
        "before": {"type": ["string", "null"]},
    },
    "required": ["after", "before"],
}


_SCAFFOLD_HINTS_SUBSCHEMA: dict = {
    "type": "array",
    "items": {"type": "string"},
}


_COMPOUND_REFS_SUBSCHEMA: dict = {
    "type": "array",
    "items": {"type": "string"},
}


_ASSAY_SELECTOR_SUBSCHEMA: dict = {
    "type": ["object", "null"],
    "additionalProperties": False,
    "properties": {
        "target_as_typed": {"type": ["string", "null"]},
        "protocol_hint": {"type": ["string", "null"]},
        "date_range": _DATE_RANGE_SUBSCHEMA,
        "created_by_as_typed": {"type": ["string", "null"]},
        "scaffold_hints": _SCAFFOLD_HINTS_SUBSCHEMA,
    },
    "required": [
        "target_as_typed", "protocol_hint", "date_range",
        "created_by_as_typed", "scaffold_hints",
    ],
}


PROMPT_SCHEMA: dict = {
    "type": "object",
    "additionalProperties": False,
    "properties": {
        "not_a_query": {"type": ["boolean", "null"]},
        "reason": {"type": ["string", "null"]},
        "registration_target_as_typed": {"type": ["string", "null"]},
        "assay_target_as_typed": {"type": ["string", "null"]},
        "registered_date_range": _DATE_RANGE_SUBSCHEMA,
        "registered_by_as_typed": {"type": ["string", "null"]},
        "scaffold_hints": _SCAFFOLD_HINTS_SUBSCHEMA,
        "compound_refs_as_typed": _COMPOUND_REFS_SUBSCHEMA,
        "rank_by": {"type": ["string", "null"]},
        "rank_top_n": {"type": ["integer", "null"]},
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
                    "assay_date_range": _DATE_RANGE_SUBSCHEMA,
                    "assayed_by_as_typed": {"type": ["string", "null"]},
                },
                "required": [
                    "protocol_hint", "metric", "threshold",
                    "assay_date_range", "assayed_by_as_typed",
                ],
            },
        },
        # Populated INSTEAD OF the top-level compound-selector fields
        # when the user is asking for a list of assays rather than
        # compounds. The LLM leaves whichever isn't applicable null.
        "assay_selector": _ASSAY_SELECTOR_SUBSCHEMA,
    },
    "required": [
        "not_a_query",
        "reason",
        "registration_target_as_typed",
        "assay_target_as_typed",
        "registered_date_range",
        "registered_by_as_typed",
        "scaffold_hints",
        "compound_refs_as_typed",
        "rank_by",
        "rank_top_n",
        "measurement_filters",
        "assay_selector",
    ],
}


def _threshold_from_dict(data: Any) -> Threshold:
    return Threshold(
        op=data["op"],
        value=float(data["value"]),
        unit=data.get("unit"),
    )


def _date_range_from_dict(data: Any) -> Optional[DateRange]:
    """Build a DateRange from an LLM-emitted `{after, before}` object.
    Returns None when the range is absent or both ends are null (treat
    as "no date filter")."""
    if data is None:
        return None
    if not isinstance(data, dict):
        raise ValueError(f"date range not a dict: {data!r}")
    after = data.get("after")
    before = data.get("before")
    if after is None and before is None:
        return None
    return DateRange(after=after, before=before)


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
        assay_date_range=_date_range_from_dict(data.get("assay_date_range")),
        assayed_by_as_typed=data.get("assayed_by_as_typed"),
    )


def _assay_selector_from_dict(data: Any) -> AssaySelector:
    if not isinstance(data, dict):
        raise ValueError(f"assay_selector not a dict: {data!r}")
    return AssaySelector(
        target_as_typed=data.get("target_as_typed"),
        protocol_hint=data.get("protocol_hint"),
        date_range=_date_range_from_dict(data.get("date_range")),
        created_by_as_typed=data.get("created_by_as_typed"),
        scaffold_hints=list(data.get("scaffold_hints") or []),
    )


def _to_parse_result(data: dict) -> PromptParseResult:
    """Map a validated LLM JSON object to a CompoundSelector /
    AssaySelector / NotAQuery / ParseError.

    Dispatch rule: if ``assay_selector`` is a non-null object, treat
    the query as an assay-selection query and build an AssaySelector.
    Otherwise build a CompoundSelector from the top-level fields."""
    if data.get("not_a_query"):
        reason = data.get("reason") or "Not a compound-selection query."
        return NotAQuery(reason=reason)

    assay_selector_raw = data.get("assay_selector")
    if assay_selector_raw is not None:
        try:
            return _assay_selector_from_dict(assay_selector_raw)
        except ValueError as e:
            return ParseError(
                message=f"Malformed assay_selector: {e}",
                raw=json.dumps(assay_selector_raw),
            )

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
        registered_date_range=_date_range_from_dict(data.get("registered_date_range")),
        registered_by_as_typed=data.get("registered_by_as_typed"),
        scaffold_hints=list(data.get("scaffold_hints") or []),
        compound_refs_as_typed=list(data.get("compound_refs_as_typed") or []),
        rank_by=data.get("rank_by"),
        rank_top_n=data.get("rank_top_n"),
        measurement_filters=filters,
    )


def _wrap_with_today(prompt: str, today: Optional[str] = None) -> str:
    """Prepend a `[Today: YYYY-MM-DD]` line to the user message so the
    LLM can resolve relative date phrasings against a known reference.

    Today is prepended to the USER message, not the system prompt —
    that way the system prompt stays static and cache-friendly (decision
    16's intent) while still giving the model a moving reference point."""
    if today is None:
        today = datetime.date.today().isoformat()
    return f"[Today: {today}]\n{prompt}"


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
            {"role": "user", "content": _wrap_with_today(prompt)},
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
