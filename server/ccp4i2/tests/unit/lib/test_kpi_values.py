"""ccp4i2.lib.kpi_values -- the gate on KPI values JSON cannot carry.

These are the properties the three layers (write gate, read gate, render
backstop) all rest on, so they are pinned here once, without a database.
"""

import json

import pytest

from ccp4i2.lib.kpi_values import (
    drop_unstorable,
    is_storable_kpi_value,
    kpi_map,
)


class _Row:
    """Stands in for a JobFloatValue/JobCharValue row.

    `key_id` rather than `key.name` because JobValueKey.name is that model's
    primary key, so the FK column already holds the name.
    """

    def __init__(self, key_id, value):
        self.key_id = key_id
        self.value = value


class TestIsStorableKpiValue:
    @pytest.mark.parametrize(
        "value", [0.0, 1.5, -1.5, 1e300, 0, 42, "0.21", "", None, True],
    )
    def test_accepts_everything_json_can_carry(self, value):
        assert is_storable_kpi_value(value)

    @pytest.mark.parametrize(
        "value", [float("nan"), float("inf"), float("-inf")],
    )
    def test_rejects_floats_json_cannot_spell(self, value):
        assert not is_storable_kpi_value(value)

    def test_rejects_nan_despite_it_comparing_unequal_to_everything(self):
        # The trap in backfill_kpis: `val != 0.0` is True for NaN, so a
        # "skip the zeroes" test lets NaN straight through.
        nan = float("nan")
        assert nan != 0.0
        assert not is_storable_kpi_value(nan)


class TestKpiMap:
    def test_builds_the_map_from_key_id(self):
        rows = [_Row("RFactor", 0.21), _Row("RFree", 0.25)]
        assert kpi_map(rows) == {"RFactor": 0.21, "RFree": 0.25}

    def test_omits_rather_than_nulls_a_non_finite_value(self):
        # Omission is load-bearing: the client tests `!== undefined` before
        # calling .toFixed(), so a null would render a confident "0.000".
        rows = [_Row("RFactor", 0.21), _Row("RFree", float("nan"))]
        result = kpi_map(rows)
        assert result == {"RFactor": 0.21}
        assert "RFree" not in result

    def test_char_rows_pass_through(self):
        assert kpi_map([_Row("spaceGroup", "P 21 21 21")]) == {
            "spaceGroup": "P 21 21 21"
        }

    def test_result_survives_strict_json(self):
        rows = [_Row("RFactor", 0.21), _Row("RFree", float("inf"))]
        # allow_nan=False is what DRF's renderer does under STRICT_JSON.
        assert json.dumps(kpi_map(rows), allow_nan=False)


class TestDropUnstorable:
    def test_keeps_the_good_and_drops_the_bad(self):
        values = {"a": 1.0, "b": float("nan"), "c": "text", "d": float("-inf")}
        assert drop_unstorable(values, context="job 3") == {"a": 1.0, "c": "text"}
