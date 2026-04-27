"""Pure-Python scorecard evaluation tests (slice 19).

The ranker mirrors the frontend's ``scorecard.ts`` exactly so a
compound's NLP-ranked position lines up with its position on the
project's spider/radar in the aggregation page UI. These tests pin
that behaviour at the per-axis-value and per-axis-normalisation
layers without a DB.
"""

from __future__ import annotations

import math

from compounds.nlp.scorecard_rank import (
    _axis_value,
    _normalise_axis,
    _row_mean_t,
    scorecard_data_needs,
)


# ---------------------------------------------------------------------------
# Axis value extraction
# ---------------------------------------------------------------------------


def test_protocol_axis_extracts_geomean():
    axis = {"kind": "protocol", "protocol_id": "p1"}
    row = {"protocols": {"p1": {"geomean": 12.5}}}
    assert _axis_value(axis, row) == 12.5


def test_protocol_axis_returns_none_for_missing_protocol():
    axis = {"kind": "protocol", "protocol_id": "p1"}
    row = {"protocols": {}}
    assert _axis_value(axis, row) is None


def test_protocol_axis_coerces_string_decimals():
    """Django decimal fields arrive as strings; the ranker should coerce."""
    axis = {"kind": "protocol", "protocol_id": "p1"}
    row = {"protocols": {"p1": {"geomean": "10.5"}}}
    assert _axis_value(axis, row) == 10.5


def test_ratio_axis_divides_geomeans():
    axis = {"kind": "ratio", "numerator_id": "p1", "denominator_id": "p2"}
    row = {"protocols": {"p1": {"geomean": 100.0}, "p2": {"geomean": 4.0}}}
    assert _axis_value(axis, row) == 25.0


def test_ratio_axis_returns_none_when_denominator_zero():
    axis = {"kind": "ratio", "numerator_id": "p1", "denominator_id": "p2"}
    row = {"protocols": {"p1": {"geomean": 100.0}, "p2": {"geomean": 0.0}}}
    assert _axis_value(axis, row) is None


def test_worst_of_axis_uses_max_when_target_lower_than_poor():
    """Lower-is-better axis: the 'worst' value is the max."""
    axis = {
        "kind": "worst_of",
        "protocol_ids": ["p1", "p2"],
        "target_value": 10, "poor_value": 1000,
    }
    row = {"protocols": {"p1": {"geomean": 50}, "p2": {"geomean": 500}}}
    assert _axis_value(axis, row) == 500


def test_worst_of_axis_uses_min_when_target_higher_than_poor():
    """Higher-is-better axis: the 'worst' value is the min."""
    axis = {
        "kind": "worst_of",
        "protocol_ids": ["p1", "p2"],
        "target_value": 1000, "poor_value": 10,
    }
    row = {"protocols": {"p1": {"geomean": 50}, "p2": {"geomean": 500}}}
    assert _axis_value(axis, row) == 50


def test_lipinski_axis_counts_passes():
    axis = {"kind": "lipinski"}
    row = {"properties": {
        "molecular_weight": 350,    # ≤500 ✓
        "clogp": 3.5,               # ≤5 ✓
        "hbd": 2,                   # ≤5 ✓
        "hba": 12,                  # >10 ✗
    }}
    assert _axis_value(axis, row) == 3.0


def test_lipinski_axis_returns_zero_when_no_properties():
    axis = {"kind": "lipinski"}
    row = {"properties": {}}
    assert _axis_value(axis, row) == 0.0


def test_unknown_axis_kind_returns_none():
    assert _axis_value({"kind": "unknown"}, {"protocols": {}}) is None


# ---------------------------------------------------------------------------
# Normalisation
# ---------------------------------------------------------------------------


def test_log_normalisation_at_target_returns_one():
    axis = {"target_value": 10, "poor_value": 1000, "threshold_scale": "log"}
    assert _normalise_axis(axis, 10) == 1.0


def test_log_normalisation_at_poor_returns_zero():
    axis = {"target_value": 10, "poor_value": 1000, "threshold_scale": "log"}
    assert _normalise_axis(axis, 1000) == 0.0


def test_log_normalisation_clamps_above_target():
    """A value better than target still maxes at 1 — no super-credit
    for being SUPER good."""
    axis = {"target_value": 10, "poor_value": 1000, "threshold_scale": "log"}
    assert _normalise_axis(axis, 1) == 1.0


def test_log_normalisation_clamps_below_poor():
    axis = {"target_value": 10, "poor_value": 1000, "threshold_scale": "log"}
    assert _normalise_axis(axis, 100000) == 0.0


def test_log_normalisation_geometric_midpoint_is_half():
    """At the geometric midpoint of target/poor, t == 0.5."""
    axis = {"target_value": 10, "poor_value": 1000, "threshold_scale": "log"}
    midpoint = math.sqrt(10 * 1000)
    assert math.isclose(_normalise_axis(axis, midpoint), 0.5, abs_tol=1e-9)


def test_linear_normalisation_arithmetic_midpoint_is_half():
    axis = {
        "target_value": 0, "poor_value": 100,
        "threshold_scale": "linear", "kind": "lipinski",
    }
    assert _normalise_axis(axis, 50) == 0.5


def test_normalisation_returns_none_for_equal_target_and_poor():
    axis = {"target_value": 10, "poor_value": 10}
    assert _normalise_axis(axis, 5) is None


# ---------------------------------------------------------------------------
# Row → mean t
# ---------------------------------------------------------------------------


def test_row_mean_t_averages_non_null_axes():
    config = {"axes": [
        {"kind": "protocol", "protocol_id": "p1",
         "target_value": 10, "poor_value": 1000, "threshold_scale": "log"},
        {"kind": "protocol", "protocol_id": "p2",
         "target_value": 10, "poor_value": 1000, "threshold_scale": "log"},
    ]}
    row = {"protocols": {
        "p1": {"geomean": 10},     # t = 1.0
        "p2": {"geomean": 1000},   # t = 0.0
    }}
    assert _row_mean_t(config, row) == 0.5


def test_row_mean_t_skips_axes_with_null_values():
    config = {"axes": [
        {"kind": "protocol", "protocol_id": "p1",
         "target_value": 10, "poor_value": 1000, "threshold_scale": "log"},
        {"kind": "protocol", "protocol_id": "p2",
         "target_value": 10, "poor_value": 1000, "threshold_scale": "log"},
    ]}
    row = {"protocols": {"p1": {"geomean": 10}}}   # p2 missing
    # Only p1 contributes — mean over one value is that value.
    assert _row_mean_t(config, row) == 1.0


def test_row_mean_t_returns_none_when_no_axis_evaluates():
    config = {"axes": [
        {"kind": "protocol", "protocol_id": "p1",
         "target_value": 10, "poor_value": 1000, "threshold_scale": "log"},
    ]}
    row = {"protocols": {}}
    assert _row_mean_t(config, row) is None


# ---------------------------------------------------------------------------
# scorecard_data_needs
# ---------------------------------------------------------------------------


def test_data_needs_collects_protocol_ids_across_axis_kinds():
    config = {"axes": [
        {"kind": "protocol", "protocol_id": "p1",
         "target_value": 1, "poor_value": 10},
        {"kind": "ratio", "numerator_id": "p2", "denominator_id": "p3",
         "target_value": 1, "poor_value": 10},
        {"kind": "worst_of", "protocol_ids": ["p4", "p5"],
         "target_value": 1, "poor_value": 10},
    ]}
    needs = scorecard_data_needs(config)
    assert set(needs["protocol_ids"]) == {"p1", "p2", "p3", "p4", "p5"}


def test_data_needs_includes_lipinski_properties():
    config = {"axes": [{"kind": "lipinski",
                        "target_value": 4, "poor_value": 0}]}
    needs = scorecard_data_needs(config)
    assert needs["properties"] == [
        "molecular_weight", "clogp", "hbd", "hba",
    ]


def test_data_needs_no_lipinski_means_no_properties():
    config = {"axes": [{"kind": "protocol", "protocol_id": "p1",
                        "target_value": 1, "poor_value": 10}]}
    needs = scorecard_data_needs(config)
    assert needs["properties"] == []


def test_data_needs_empty_config_returns_empty():
    needs = scorecard_data_needs({})
    assert needs == {"protocol_ids": [], "properties": []}
