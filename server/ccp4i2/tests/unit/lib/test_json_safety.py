"""ccp4i2.lib.json_safety -- the shared scrub for floats JSON cannot spell.

This is the one implementation behind `digest.json_safe`, the KPI read gate and
`SafeJSONRenderer`, so its properties are pinned here. CCP4-free: stdlib only.
"""

import json
import math

import pytest

from ccp4i2.lib.json_safety import is_finite_number, replace_non_finite


class TestIsFiniteNumber:
    @pytest.mark.parametrize("value", [0.0, -1.5, 1e300, 0, 42, "nan", None, True])
    def test_true_for_anything_json_can_carry(self, value):
        assert is_finite_number(value)

    @pytest.mark.parametrize("value", [float("nan"), float("inf"), float("-inf")])
    def test_false_for_the_three_floats_json_cannot(self, value):
        assert not is_finite_number(value)


class TestReplaceNonFinite:
    def test_nulls_a_non_finite_float(self):
        assert replace_non_finite(float("nan")) is None
        assert replace_non_finite(float("inf")) is None
        assert replace_non_finite(float("-inf")) is None

    def test_leaves_ordinary_values_alone(self):
        assert replace_non_finite(0.21) == 0.21
        assert replace_non_finite("nan") == "nan"
        assert replace_non_finite(None) is None

    def test_walks_nested_structures(self):
        data = {
            "jobs": [
                {"kpis": {"RFree": float("nan"), "RFactor": 0.21}},
                {"kpis": {"RFree": 0.25}},
            ],
            "total": 2,
        }
        assert replace_non_finite(data) == {
            "jobs": [
                {"kpis": {"RFree": None, "RFactor": 0.21}},
                {"kpis": {"RFree": 0.25}},
            ],
            "total": 2,
        }

    def test_output_survives_strict_json(self):
        # allow_nan=False is what DRF's renderer does under STRICT_JSON.
        deep = {"a": [1.0, float("nan"), {"b": (float("inf"), 2.0)}]}
        assert json.loads(json.dumps(replace_non_finite(deep), allow_nan=False)) == {
            "a": [1.0, None, {"b": [None, 2.0]}]
        }

    def test_does_not_mutate_its_input(self):
        data = {"x": float("nan")}
        replace_non_finite(data)
        assert math.isnan(data["x"])


def test_digest_json_safe_is_this_function():
    """`digest.json_safe` kept its name but delegates here; one behaviour."""
    gemmi = pytest.importorskip("gemmi", reason="digest imports gemmi")
    from ccp4i2.lib.utils.files.digest import json_safe

    raw = {"a": float("nan"), "b": [1.5, float("inf")], "c": "x"}
    assert json_safe(raw) == replace_non_finite(raw) == {
        "a": None, "b": [1.5, None], "c": "x"
    }
