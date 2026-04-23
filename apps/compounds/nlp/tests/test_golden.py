"""Golden-set evaluation (slice 8, §12).

The YAML at ``golden/prompts.yaml`` is a curated specification of what the LLM
must emit for a range of prompts. Two kinds of test live here:

* **Offline** (run in every pytest pass): structural validation of the golden
  file + a round-trip check that each entry's expected shape is reachable
  from ``_to_parse_result`` given a correct LLM output. Protects against YAML
  typos and drift as the set grows. Does NOT touch the network.
* **Online** (skipped unless ``RUN_NLP_ONLINE_EVAL=true``): actually calls
  Azure OpenAI via ``parse_prompt`` for each entry and measures agreement
  against the golden specs. Success bar per §12:
    - ≥90% exact-match rate on ``type: spec`` / ``type: spec_with_ambiguity``
    - 100% of ``type: not_a_query`` entries correctly rejected

Run the online eval manually:

    source ../../ccp4-20251105/bin/ccp4.setup-sh
    AZURE_OPENAI_ENDPOINT=https://ddu-openai.openai.azure.com/ \\
    RUN_NLP_ONLINE_EVAL=true \\
    PYTHONPATH="$PWD:$PWD/../apps" DJANGO_SETTINGS_MODULE=compounds.settings \\
      ccp4-python -m pytest ../apps/compounds/nlp/tests/test_golden.py::test_online_eval -v -s

Requires ``openai`` installed locally plus either managed-identity via
``az login`` or ``AZURE_OPENAI_API_KEY`` set as the dev-only fallback.
"""

from __future__ import annotations

import math
import os
from pathlib import Path
from typing import List, Optional

import pytest
import yaml

from compounds.nlp.llm import _to_parse_result
from compounds.nlp.spec import NotAQuery, ParseError, QuerySpec, Threshold


GOLDEN_PATH = Path(__file__).parent / "golden" / "prompts.yaml"

VALID_ENTRY_TYPES = {"spec", "spec_with_ambiguity", "not_a_query"}
VALID_THRESHOLD_OPS = {"<", "<=", ">", ">=", "=", "==", "!="}
VALID_COLUMN_PRESETS = {"phys_chem", "lipinski", None}

# Online eval gate — default off so CI / local runs don't burn LLM calls.
RUN_ONLINE = os.environ.get("RUN_NLP_ONLINE_EVAL", "").lower() in ("true", "1", "yes")

# Thresholds from §12.
SPEC_MATCH_BAR = 0.90           # at least 90% exact spec matches
NOT_A_QUERY_MATCH_BAR = 1.0     # 100% — any "parses as spec" is a false positive


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
# Helpers — translate an expected entry to the shape _to_parse_result consumes,
# and diff two specs for human-readable mismatch messages.
# ---------------------------------------------------------------------------


def _expected_to_llm_output(expected: dict) -> dict:
    """Build the JSON object the LLM would emit for an entry's expected shape.
    Strict-mode schema requires every property listed even if null."""
    t = expected.get("type", "spec")
    base = {
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
    if t == "not_a_query":
        base["not_a_query"] = True
        base["reason"] = expected.get("reason") or "Not a tabular query."
        return base
    for key in (
        "registration_target_as_typed", "assay_target_as_typed",
        "protocol_hint", "metric", "columns", "columns_explicit",
    ):
        if key in expected:
            base[key] = expected[key]
    threshold = expected.get("threshold")
    if threshold is not None:
        base["threshold"] = {
            "op": threshold["op"],
            "value": threshold["value"],
            "unit": threshold.get("unit"),
        }
    return base


def _spec_diff(actual: QuerySpec, expected: dict) -> Optional[str]:
    """Return None if actual matches expected; otherwise a ``"; "``-joined diff."""
    diffs: List[str] = []
    for field in (
        "registration_target_as_typed", "assay_target_as_typed",
        "protocol_hint", "metric", "columns", "columns_explicit",
    ):
        a = getattr(actual, field)
        e = expected.get(field)
        if a != e:
            diffs.append(f"{field}: expected {e!r}, got {a!r}")

    a_t = actual.threshold
    e_t = expected.get("threshold")
    if (a_t is None) != (e_t is None):
        diffs.append(f"threshold: expected {e_t!r}, got {a_t!r}")
    elif a_t is not None and e_t is not None:
        if a_t.op != e_t["op"]:
            diffs.append(f"threshold.op: expected {e_t['op']!r}, got {a_t.op!r}")
        if not math.isclose(a_t.value, float(e_t["value"])):
            diffs.append(f"threshold.value: expected {e_t['value']}, got {a_t.value}")
        if a_t.unit != e_t.get("unit"):
            diffs.append(f"threshold.unit: expected {e_t.get('unit')!r}, got {a_t.unit!r}")

    return "; ".join(diffs) if diffs else None


# ---------------------------------------------------------------------------
# Offline structural validation
# ---------------------------------------------------------------------------


def test_golden_file_is_nonempty():
    assert len(GOLDEN) >= 20, (
        f"Golden set has {len(GOLDEN)} entries — §12 calls for ~30+ to be meaningful. "
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
        t = entry["expected"].get("type", "spec")
        assert t in VALID_ENTRY_TYPES, (
            f"entry {entry['name']!r}: unknown type {t!r} — must be one of {VALID_ENTRY_TYPES}"
        )


def test_every_threshold_uses_a_valid_op():
    for entry in GOLDEN:
        threshold = entry["expected"].get("threshold")
        if threshold is None:
            continue
        op = threshold.get("op")
        assert op in VALID_THRESHOLD_OPS, (
            f"entry {entry['name']!r}: threshold.op={op!r} not in {VALID_THRESHOLD_OPS}"
        )


def test_every_columns_value_is_a_valid_preset_or_null():
    for entry in GOLDEN:
        columns = entry["expected"].get("columns")
        assert columns in VALID_COLUMN_PRESETS, (
            f"entry {entry['name']!r}: columns={columns!r} not in {VALID_COLUMN_PRESETS}"
        )


def test_every_spec_entry_has_a_metric():
    for entry in GOLDEN:
        t = entry["expected"].get("type", "spec")
        if t not in ("spec", "spec_with_ambiguity"):
            continue
        # Every spec entry needs a metric so the executor can evaluate rows.
        assert entry["expected"].get("metric"), (
            f"entry {entry['name']!r}: spec entries must name a metric"
        )


def test_every_spec_entry_names_at_least_one_target():
    for entry in GOLDEN:
        t = entry["expected"].get("type", "spec")
        if t not in ("spec", "spec_with_ambiguity"):
            continue
        reg = entry["expected"].get("registration_target_as_typed")
        assay = entry["expected"].get("assay_target_as_typed")
        assert reg or assay, (
            f"entry {entry['name']!r}: at least one target field must be non-null (§6.5)"
        )


def test_columns_explicit_and_columns_are_mutually_consistent():
    # columns_explicit overrides columns at executor time; it's fine for both
    # to be set (the LLM may default `columns` to phys_chem even when naming
    # explicit fields), but not all-null — a spec with neither means the
    # executor's phys_chem default kicks in, which is also fine. This test
    # just sanity-checks the shapes aren't malformed.
    for entry in GOLDEN:
        t = entry["expected"].get("type", "spec")
        if t not in ("spec", "spec_with_ambiguity"):
            continue
        columns_explicit = entry["expected"].get("columns_explicit")
        if columns_explicit is not None:
            assert isinstance(columns_explicit, list), (
                f"entry {entry['name']!r}: columns_explicit must be a list or null"
            )
            assert all(isinstance(c, str) for c in columns_explicit), (
                f"entry {entry['name']!r}: columns_explicit entries must be strings"
            )


def test_every_spec_entry_round_trips_via_to_parse_result():
    """A correct LLM output built from an expected spec must parse back to that spec."""
    for entry in GOLDEN:
        t = entry["expected"].get("type", "spec")
        if t not in ("spec", "spec_with_ambiguity"):
            continue
        llm_output = _expected_to_llm_output(entry["expected"])
        result = _to_parse_result(llm_output)
        assert isinstance(result, QuerySpec), (
            f"entry {entry['name']!r}: round-trip produced {type(result).__name__}, "
            f"not QuerySpec. Raw: {result!r}"
        )
        diff = _spec_diff(result, entry["expected"])
        assert diff is None, f"entry {entry['name']!r}: round-trip diff: {diff}"


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
# Online eval — skipped unless RUN_NLP_ONLINE_EVAL=true and the endpoint env
# is populated.
# ---------------------------------------------------------------------------


@pytest.mark.skipif(
    not RUN_ONLINE,
    reason="online eval disabled; set RUN_NLP_ONLINE_EVAL=true to run against Azure",
)
def test_online_eval():
    """Iterate the golden set, call the real LLM, measure agreement."""
    # Imported here so the module stays importable when `openai` isn't installed.
    from compounds.nlp.llm import parse_prompt

    if not os.environ.get("AZURE_OPENAI_ENDPOINT"):
        pytest.skip("AZURE_OPENAI_ENDPOINT not set — can't reach Azure OpenAI")

    spec_entries = [e for e in GOLDEN if e["expected"].get("type", "spec") in ("spec", "spec_with_ambiguity")]
    naq_entries = [e for e in GOLDEN if e["expected"].get("type") == "not_a_query"]

    spec_matches = 0
    spec_mismatches: List[tuple] = []
    for entry in spec_entries:
        result = parse_prompt(entry["prompt"])
        if isinstance(result, QuerySpec):
            diff = _spec_diff(result, entry["expected"])
            if diff is None:
                spec_matches += 1
            else:
                spec_mismatches.append((entry["name"], entry["prompt"], diff))
        else:
            spec_mismatches.append(
                (entry["name"], entry["prompt"], f"expected QuerySpec, got {type(result).__name__}: {result!r}")
            )

    naq_matches = 0
    naq_mismatches: List[tuple] = []
    for entry in naq_entries:
        result = parse_prompt(entry["prompt"])
        if isinstance(result, NotAQuery):
            naq_matches += 1
        else:
            naq_mismatches.append(
                (entry["name"], entry["prompt"], f"expected NotAQuery, got {type(result).__name__}: {result!r}")
            )

    # Report — this runs with -s so it lands on stdout rather than being swallowed.
    print()
    print("=" * 72)
    print(f"Spec match rate:         {spec_matches}/{len(spec_entries)} = {spec_matches / len(spec_entries):.1%}   (bar ≥ {SPEC_MATCH_BAR:.0%})")
    print(f"Not-a-query match rate:  {naq_matches}/{len(naq_entries)} = {naq_matches / len(naq_entries):.1%}   (bar = {NOT_A_QUERY_MATCH_BAR:.0%})")
    print("=" * 72)
    if spec_mismatches:
        print("\nSpec mismatches:")
        for name, prompt, diff in spec_mismatches:
            print(f"  • {name}: {prompt!r}")
            print(f"      {diff}")
    if naq_mismatches:
        print("\nNot-a-query mismatches (LLM gave us a spec where it shouldn't have):")
        for name, prompt, msg in naq_mismatches:
            print(f"  • {name}: {prompt!r}")
            print(f"      {msg}")

    # Bar assertions last so the report prints even on failure.
    spec_rate = spec_matches / len(spec_entries) if spec_entries else 1.0
    naq_rate = naq_matches / len(naq_entries) if naq_entries else 1.0
    assert spec_rate >= SPEC_MATCH_BAR, (
        f"Spec match rate {spec_rate:.1%} below bar {SPEC_MATCH_BAR:.0%}"
    )
    assert naq_rate >= NOT_A_QUERY_MATCH_BAR, (
        f"Not-a-query match rate {naq_rate:.1%} below bar {NOT_A_QUERY_MATCH_BAR:.0%}"
    )


@pytest.mark.skipif(
    not RUN_ONLINE,
    reason="online eval disabled; set RUN_NLP_ONLINE_EVAL=true to run against Azure",
)
def test_online_parse_error_rate():
    """Catch systematic LLM / glue failures separately from spec mismatches.

    A ParseError means our glue got junk from the upstream — distinct from the
    LLM emitting a spec that simply disagrees with the golden file. We tolerate
    none of these.
    """
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
