"""Golden-set evaluation (§12).

YAML at ``golden/prompts.yaml`` is the curated specification of what the
LLM must emit for a range of prompts. Two kinds of test live here:

* **Offline** (runs every pytest pass, ~40ms): structural validation of
  the golden file + round-trip each entry's expected shape through
  ``_to_parse_result`` to verify it parses back to the expected
  CompoundSelector / NotAQuery. Protects against YAML typos and drift.
  Does NOT touch the network.
* **Online** (skipped unless ``RUN_NLP_ONLINE_EVAL=true``): calls
  Azure OpenAI via ``parse_prompt`` for each entry and measures match
  rate against the golden set. §12 bars: ≥90% selector match, 100%
  not_a_query detection, zero parse errors.

Run the online eval (from a server container with openai + MI access):

    az containerapp exec --name ccp4i2-bicep-server \\
        --resource-group ccp4i2-bicep-rg-uksouth --command /bin/bash

    # inside the shell
    RUN_NLP_ONLINE_EVAL=true /usr/src/app/entrypoint.sh \\
      ccp4-python -m pytest /usr/src/app/compounds/nlp/tests/test_golden.py -v -s
"""

from __future__ import annotations

import math
import os
from pathlib import Path
from typing import List, Optional

import pytest
import yaml

from compounds.nlp.llm import _to_parse_result
from compounds.nlp.spec import (
    AssaySelector,
    CompoundSelector,
    DateRange,
    MeasurementFilter,
    NotAQuery,
    ParseError,
    Threshold,
)


GOLDEN_PATH = Path(__file__).parent / "golden" / "prompts.yaml"

VALID_ENTRY_TYPES = {
    "selector", "selector_with_clarify", "assay_selector", "not_a_query",
}
VALID_THRESHOLD_OPS = {"<", "<=", ">", ">=", "=", "==", "!="}

RUN_ONLINE = os.environ.get("RUN_NLP_ONLINE_EVAL", "").lower() in ("true", "1", "yes")

SPEC_MATCH_BAR = 0.90
NOT_A_QUERY_MATCH_BAR = 1.0


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------


def _load_golden() -> List[dict]:
    with GOLDEN_PATH.open() as f:
        data = yaml.safe_load(f)
    assert isinstance(data, list), f"{GOLDEN_PATH.name} must be a YAML list at the top level"
    return data


GOLDEN = _load_golden()


# ---------------------------------------------------------------------------
# Helpers — translate an expected entry to the LLM-output shape, diff selectors
# ---------------------------------------------------------------------------


def _expected_to_llm_output(expected: dict) -> dict:
    """Build the JSON object the LLM would emit for an entry's expected shape.
    Strict-mode schema requires every property listed even if null."""
    t = expected.get("type", "selector")
    base = {
        "not_a_query": None,
        "reason": None,
        "registration_target_as_typed": None,
        "assay_target_as_typed": None,
        "registered_date_range": None,
        "registered_by_as_typed": None,
        "scaffold_hints": [],
        "compound_refs_as_typed": [],
        "measurement_filters": [],
        "assay_selector": None,
    }
    if t == "not_a_query":
        base["not_a_query"] = True
        base["reason"] = expected.get("reason") or "Not a compound-selection query."
        return base

    if t == "assay_selector":
        raw = expected.get("assay_selector") or {}
        base["assay_selector"] = {
            "target_as_typed": raw.get("target_as_typed"),
            "protocol_hint": raw.get("protocol_hint"),
            "date_range": _date_range_payload(raw.get("date_range")),
            "created_by_as_typed": raw.get("created_by_as_typed"),
            "scaffold_hints": list(raw.get("scaffold_hints") or []),
        }
        return base

    if "registration_target_as_typed" in expected:
        base["registration_target_as_typed"] = expected["registration_target_as_typed"]
    if "assay_target_as_typed" in expected:
        base["assay_target_as_typed"] = expected["assay_target_as_typed"]
    if "registered_date_range" in expected:
        base["registered_date_range"] = _date_range_payload(expected["registered_date_range"])
    if "registered_by_as_typed" in expected:
        base["registered_by_as_typed"] = expected["registered_by_as_typed"]
    if "scaffold_hints" in expected:
        base["scaffold_hints"] = list(expected["scaffold_hints"])
    if "compound_refs_as_typed" in expected:
        base["compound_refs_as_typed"] = list(expected["compound_refs_as_typed"])

    raw_filters = expected.get("measurement_filters") or []
    filters: List[dict] = []
    for flt in raw_filters:
        threshold = flt.get("threshold")
        threshold_out = None if threshold is None else {
            "op": threshold["op"],
            "value": threshold["value"],
            "unit": threshold.get("unit"),
        }
        filters.append({
            "protocol_hint": flt.get("protocol_hint"),
            "metric": flt.get("metric"),
            "threshold": threshold_out,
            "assay_date_range": _date_range_payload(flt.get("assay_date_range")),
            "assayed_by_as_typed": flt.get("assayed_by_as_typed"),
        })
    base["measurement_filters"] = filters
    return base


def _date_range_payload(dr: Optional[dict]) -> Optional[dict]:
    """Convert a golden-YAML date-range dict to the strict-mode LLM-output
    shape. Returns None when the range itself is null/absent; otherwise
    returns `{after, before}` with both keys (null if unset)."""
    if dr is None:
        return None
    return {"after": dr.get("after"), "before": dr.get("before")}


def _date_range_diff(
    actual: Optional[DateRange],
    expected: Optional[dict],
    path: str,
) -> List[str]:
    """Diff a golden-YAML date-range block against an actual DateRange.
    Treat both-null / both-absent as equivalent."""
    a_after = actual.after if actual is not None else None
    a_before = actual.before if actual is not None else None
    e_after = expected.get("after") if expected else None
    e_before = expected.get("before") if expected else None
    diffs: List[str] = []
    if a_after != e_after:
        diffs.append(f"{path}.after: expected {e_after!r}, got {a_after!r}")
    if a_before != e_before:
        diffs.append(f"{path}.before: expected {e_before!r}, got {a_before!r}")
    return diffs


def _filter_diff(
    actual: MeasurementFilter,
    expected: dict,
    idx: int,
) -> List[str]:
    diffs: List[str] = []
    for field in ("protocol_hint", "metric", "assayed_by_as_typed"):
        a = getattr(actual, field)
        e = expected.get(field)
        if a != e:
            diffs.append(f"filter[{idx}].{field}: expected {e!r}, got {a!r}")
    a_t = actual.threshold
    e_t = expected.get("threshold")
    if (a_t is None) != (e_t is None):
        diffs.append(f"filter[{idx}].threshold: expected {e_t!r}, got {a_t!r}")
    elif a_t is not None and e_t is not None:
        if a_t.op != e_t["op"]:
            diffs.append(f"filter[{idx}].threshold.op: expected {e_t['op']!r}, got {a_t.op!r}")
        if not math.isclose(a_t.value, float(e_t["value"])):
            diffs.append(f"filter[{idx}].threshold.value: expected {e_t['value']}, got {a_t.value}")
        if a_t.unit != e_t.get("unit"):
            diffs.append(f"filter[{idx}].threshold.unit: expected {e_t.get('unit')!r}, got {a_t.unit!r}")
    diffs.extend(_date_range_diff(
        actual.assay_date_range,
        expected.get("assay_date_range"),
        f"filter[{idx}].assay_date_range",
    ))
    return diffs


def _selector_diff(actual: CompoundSelector, expected: dict) -> Optional[str]:
    diffs: List[str] = []
    for field in (
        "registration_target_as_typed",
        "assay_target_as_typed",
        "registered_by_as_typed",
    ):
        a = getattr(actual, field)
        e = expected.get(field)
        if a != e:
            diffs.append(f"{field}: expected {e!r}, got {a!r}")

    diffs.extend(_date_range_diff(
        actual.registered_date_range,
        expected.get("registered_date_range"),
        "registered_date_range",
    ))

    expected_scaffolds = list(expected.get("scaffold_hints") or [])
    actual_scaffolds = list(actual.scaffold_hints or [])
    if expected_scaffolds != actual_scaffolds:
        diffs.append(
            f"scaffold_hints: expected {expected_scaffolds!r}, got {actual_scaffolds!r}"
        )

    expected_pins = list(expected.get("compound_refs_as_typed") or [])
    actual_pins = list(actual.compound_refs_as_typed or [])
    if expected_pins != actual_pins:
        diffs.append(
            f"compound_refs_as_typed: expected {expected_pins!r}, got {actual_pins!r}"
        )

    expected_filters = expected.get("measurement_filters") or []
    if len(actual.measurement_filters) != len(expected_filters):
        diffs.append(
            f"measurement_filters length: expected {len(expected_filters)}, "
            f"got {len(actual.measurement_filters)}"
        )
    else:
        for idx, (a_flt, e_flt) in enumerate(zip(actual.measurement_filters, expected_filters)):
            diffs.extend(_filter_diff(a_flt, e_flt, idx))

    return "; ".join(diffs) if diffs else None


# ---------------------------------------------------------------------------
# Offline structural validation
# ---------------------------------------------------------------------------


def test_golden_file_is_nonempty():
    assert len(GOLDEN) >= 20, (
        f"Golden set has {len(GOLDEN)} entries — §12 calls for ~30+. "
        "Add more prompts to golden/prompts.yaml."
    )


def test_every_entry_has_required_top_level_keys():
    for i, entry in enumerate(GOLDEN):
        assert isinstance(entry, dict), f"entry #{i} is not a dict: {entry!r}"
        for key in ("name", "prompt", "expected"):
            assert key in entry, f"entry #{i} missing '{key}': {entry!r}"


def test_entry_names_are_unique():
    names = [e["name"] for e in GOLDEN]
    duplicates = {n for n in names if names.count(n) > 1}
    assert not duplicates, f"duplicate golden entry names: {duplicates}"


def test_every_expected_type_is_valid():
    for entry in GOLDEN:
        t = entry["expected"].get("type", "selector")
        assert t in VALID_ENTRY_TYPES, (
            f"entry {entry['name']!r}: unknown type {t!r} — must be one of {VALID_ENTRY_TYPES}"
        )


def test_every_threshold_uses_a_valid_op():
    for entry in GOLDEN:
        for flt in entry["expected"].get("measurement_filters") or []:
            threshold = flt.get("threshold")
            if threshold is None:
                continue
            op = threshold.get("op")
            assert op in VALID_THRESHOLD_OPS, (
                f"entry {entry['name']!r}: threshold.op={op!r} not in {VALID_THRESHOLD_OPS}"
            )


def test_every_selector_entry_has_some_narrowing_predicate(db):
    """A selector entry must either name a target or include at least
    one other narrowing predicate (scaffold / user / date / measurement
    filter). A fully-empty selector would be a golden-set typo, not a
    legitimate chemist prompt."""
    for entry in GOLDEN:
        t = entry["expected"].get("type", "selector")
        if t not in ("selector", "selector_with_clarify"):
            continue
        reg = entry["expected"].get("registration_target_as_typed")
        assay = entry["expected"].get("assay_target_as_typed")
        if reg or assay:
            continue
        has_other_narrowing = any([
            entry["expected"].get("scaffold_hints"),
            entry["expected"].get("registered_by_as_typed"),
            entry["expected"].get("registered_date_range"),
            entry["expected"].get("measurement_filters"),
            entry["expected"].get("compound_refs_as_typed"),
        ])
        assert has_other_narrowing, (
            f"entry {entry['name']!r}: target-less selector must have at "
            "least one other narrowing predicate (scaffold / user / date / filter)"
        )


def test_threshold_with_value_has_numeric_value():
    for entry in GOLDEN:
        for flt in entry["expected"].get("measurement_filters") or []:
            threshold = flt.get("threshold")
            if threshold is None:
                continue
            assert isinstance(threshold.get("value"), (int, float)), (
                f"entry {entry['name']!r}: threshold.value must be numeric"
            )


def test_every_selector_entry_round_trips_via_to_parse_result():
    """A correct LLM output built from an expected selector must parse back
    to a CompoundSelector matching the entry."""
    for entry in GOLDEN:
        t = entry["expected"].get("type", "selector")
        if t not in ("selector", "selector_with_clarify"):
            continue
        llm_output = _expected_to_llm_output(entry["expected"])
        result = _to_parse_result(llm_output)
        assert isinstance(result, CompoundSelector), (
            f"entry {entry['name']!r}: round-trip produced {type(result).__name__}, not CompoundSelector"
        )
        diff = _selector_diff(result, entry["expected"])
        assert diff is None, f"entry {entry['name']!r}: round-trip diff: {diff}"


def test_every_assay_selector_entry_round_trips():
    """Assay-selection entries must parse back to an AssaySelector (not
    a CompoundSelector) via the `assay_selector != null` dispatch."""
    for entry in GOLDEN:
        if entry["expected"].get("type") != "assay_selector":
            continue
        llm_output = _expected_to_llm_output(entry["expected"])
        result = _to_parse_result(llm_output)
        assert isinstance(result, AssaySelector), (
            f"entry {entry['name']!r}: produced {type(result).__name__}, not AssaySelector"
        )
        expected_as = entry["expected"].get("assay_selector") or {}
        for field in ("target_as_typed", "protocol_hint", "created_by_as_typed"):
            assert getattr(result, field) == expected_as.get(field), (
                f"entry {entry['name']!r}: {field} diff "
                f"— expected {expected_as.get(field)!r}, got {getattr(result, field)!r}"
            )
        expected_date = expected_as.get("date_range")
        actual_date = result.date_range
        if expected_date is None:
            assert actual_date is None, f"entry {entry['name']!r}: date_range should be null"
        else:
            assert actual_date is not None, f"entry {entry['name']!r}: date_range missing"
            assert actual_date.after == expected_date.get("after"), (
                f"entry {entry['name']!r}: date_range.after diff"
            )
            assert actual_date.before == expected_date.get("before"), (
                f"entry {entry['name']!r}: date_range.before diff"
            )
        expected_scaffolds = list(expected_as.get("scaffold_hints") or [])
        actual_scaffolds = list(result.scaffold_hints or [])
        assert expected_scaffolds == actual_scaffolds, (
            f"entry {entry['name']!r}: scaffold_hints diff "
            f"— expected {expected_scaffolds!r}, got {actual_scaffolds!r}"
        )


def test_every_not_a_query_entry_round_trips():
    for entry in GOLDEN:
        if entry["expected"].get("type") != "not_a_query":
            continue
        llm_output = _expected_to_llm_output(entry["expected"])
        result = _to_parse_result(llm_output)
        assert isinstance(result, NotAQuery), (
            f"entry {entry['name']!r}: round-trip produced {type(result).__name__}, not NotAQuery"
        )


# ---------------------------------------------------------------------------
# Online eval — gated on RUN_NLP_ONLINE_EVAL + AZURE_OPENAI_ENDPOINT
# ---------------------------------------------------------------------------


@pytest.mark.skipif(
    not RUN_ONLINE,
    reason="online eval disabled; set RUN_NLP_ONLINE_EVAL=true to run against Azure",
)
def test_online_eval():
    from compounds.nlp.llm import parse_prompt

    if not os.environ.get("AZURE_OPENAI_ENDPOINT"):
        pytest.skip("AZURE_OPENAI_ENDPOINT not set")

    selector_entries = [
        e for e in GOLDEN
        if e["expected"].get("type", "selector") in ("selector", "selector_with_clarify")
    ]
    naq_entries = [e for e in GOLDEN if e["expected"].get("type") == "not_a_query"]

    sel_matches = 0
    sel_mismatches: List[tuple] = []
    for entry in selector_entries:
        result = parse_prompt(entry["prompt"])
        if isinstance(result, CompoundSelector):
            diff = _selector_diff(result, entry["expected"])
            if diff is None:
                sel_matches += 1
            else:
                sel_mismatches.append((entry["name"], entry["prompt"], diff))
        else:
            sel_mismatches.append(
                (entry["name"], entry["prompt"],
                 f"expected CompoundSelector, got {type(result).__name__}: {result!r}")
            )

    naq_matches = 0
    naq_mismatches: List[tuple] = []
    for entry in naq_entries:
        result = parse_prompt(entry["prompt"])
        if isinstance(result, NotAQuery):
            naq_matches += 1
        else:
            naq_mismatches.append(
                (entry["name"], entry["prompt"],
                 f"expected NotAQuery, got {type(result).__name__}: {result!r}")
            )

    sel_rate = sel_matches / len(selector_entries) if selector_entries else 1.0
    naq_rate = naq_matches / len(naq_entries) if naq_entries else 1.0

    print()
    print("=" * 72)
    print(f"Selector match rate:     {sel_matches}/{len(selector_entries)} = {sel_rate:.1%}   (bar ≥ {SPEC_MATCH_BAR:.0%})")
    print(f"Not-a-query match rate:  {naq_matches}/{len(naq_entries)} = {naq_rate:.1%}   (bar = {NOT_A_QUERY_MATCH_BAR:.0%})")
    print("=" * 72)
    if sel_mismatches:
        print("\nSelector mismatches:")
        for name, prompt, diff in sel_mismatches:
            print(f"  • {name}: {prompt!r}")
            print(f"      {diff}")
    if naq_mismatches:
        print("\nNot-a-query mismatches (LLM gave us a selector where it shouldn't have):")
        for name, prompt, msg in naq_mismatches:
            print(f"  • {name}: {prompt!r}")
            print(f"      {msg}")

    assert sel_rate >= SPEC_MATCH_BAR, (
        f"Selector match rate {sel_rate:.1%} below bar {SPEC_MATCH_BAR:.0%}"
    )
    assert naq_rate >= NOT_A_QUERY_MATCH_BAR, (
        f"Not-a-query match rate {naq_rate:.1%} below bar {NOT_A_QUERY_MATCH_BAR:.0%}"
    )


@pytest.mark.skipif(
    not RUN_ONLINE,
    reason="online eval disabled; set RUN_NLP_ONLINE_EVAL=true to run against Azure",
)
def test_online_parse_error_rate():
    from compounds.nlp.llm import parse_prompt

    if not os.environ.get("AZURE_OPENAI_ENDPOINT"):
        pytest.skip("AZURE_OPENAI_ENDPOINT not set")

    parse_errors: List[tuple] = []
    for entry in GOLDEN:
        result = parse_prompt(entry["prompt"])
        if isinstance(result, ParseError):
            parse_errors.append((entry["name"], entry["prompt"], result.message))

    if parse_errors:
        print("\nParseError cases (glue-code / LLM-output failures):")
        for name, prompt, msg in parse_errors:
            print(f"  • {name}: {prompt!r} → {msg}")
    assert not parse_errors, (
        f"{len(parse_errors)} prompts produced ParseError — expected 0 (schema/glue issue)"
    )
