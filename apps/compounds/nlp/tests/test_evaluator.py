"""Row-level evaluator tests (§6.3 metric, §6.4 units, Q25 tri-state).

These tests deliberately avoid the ORM — evaluate_row() only reads .status
and .results, so a SimpleNamespace stands in for AnalysisResult. Keeps the
suite fast and makes the logic easy to reason about.
"""

from __future__ import annotations

import math
from types import SimpleNamespace

import pytest

from compounds.nlp.evaluator import (
    QUERY_UNIT_NONE,
    QUERY_UNIT_UNKNOWN,
    ROW_UNIT_UNITLESS,
    ROW_UNIT_UNKNOWN,
    classify_query_unit,
    classify_row_unit,
    convert,
    evaluate_row,
)
from compounds.nlp.spec import (
    EXCLUDE_QUERY_MISSING_UNIT,
    EXCLUDE_UNIT_INCOMPATIBLE,
    EXCLUDE_UNIT_TYPE_MISMATCH,
    EXCLUDE_UNIT_UNKNOWN,
    FILTER_INVALID_STATUS,
    FILTER_KPI_MISMATCH,
    FILTER_THRESHOLD_NOT_MET,
    FILTER_VALUE_NOT_NUMERIC,
    RowExcluded,
    RowFiltered,
    RowMatched,
    Threshold,
)


def row(*, status="valid", kpi="IC50", value=None, kpi_unit="nM", **extra):
    """Build a lightweight AnalysisResult stand-in."""
    results = {"KPI": kpi, **extra}
    if value is not None:
        results[kpi] = value
    if kpi_unit is not ...:  # allow absence by passing ...
        if kpi_unit is None:
            results["kpi_unit"] = None
        elif kpi_unit != "__absent__":
            results["kpi_unit"] = kpi_unit
    return SimpleNamespace(status=status, results=results)


# ===========================================================================
# convert()
# ===========================================================================


@pytest.mark.parametrize(
    "value,from_u,to_u,expected",
    [
        (1.0, "uM", "nM", 1000.0),
        (1000.0, "nM", "uM", 1.0),
        (0.005, "uM", "mM", 5e-6),
        (1.0, "M", "nM", 1e9),
        (60.0, "s", "min", 1.0),
        (1.0, "h", "min", 60.0),
        (1.0, "h", "s", 3600.0),
        (1.0, "cm/s", "1e-6 cm/s", 1e6),
        (50.0, "%", "%", 50.0),      # identity
        (10.0, "nM", "nM", 10.0),    # identity
    ],
)
def test_convert_within_family(value, from_u, to_u, expected):
    result = convert(value, from_u, to_u)
    assert result is not None
    assert math.isclose(result, expected, rel_tol=1e-12, abs_tol=1e-20)


@pytest.mark.parametrize(
    "from_u,to_u",
    [
        ("nM", "min"),         # conc vs time
        ("%", "nM"),           # percent vs conc
        ("uL/min/mg", "mL/min/kg"),  # different clearance families
        ("cm/s", "nM"),        # permeability vs conc
        ("unknown", "nM"),     # unknown unit
        ("nM", "unknown"),
    ],
)
def test_convert_cross_family_returns_none(from_u, to_u):
    assert convert(1.0, from_u, to_u) is None


# ===========================================================================
# classify_row_unit() — the tri-state
# ===========================================================================


@pytest.mark.parametrize(
    "kpi_unit,expected",
    [
        ("nM", "nM"),
        ("uM", "uM"),
        ("µM", "uM"),          # unicode mu normalises to 'u'
        ("μM", "uM"),          # greek mu also
        ("  nM  ", "nM"),      # whitespace-tolerant
        ("unitless", ROW_UNIT_UNITLESS),
        ("UNITLESS", ROW_UNIT_UNITLESS),
        (None, ROW_UNIT_UNITLESS),      # backward compat — None means unitless
        ("", ROW_UNIT_UNKNOWN),         # empty string = absent kpi_unit
        ("  ", ROW_UNIT_UNKNOWN),
        ("bananas", ROW_UNIT_UNKNOWN),  # unrecognised
        (5, ROW_UNIT_UNKNOWN),          # non-string
    ],
)
def test_classify_row_unit(kpi_unit, expected):
    assert classify_row_unit(kpi_unit) == expected


# ===========================================================================
# classify_query_unit()
# ===========================================================================


@pytest.mark.parametrize(
    "unit,expected",
    [
        ("nM", "nM"),
        ("µM", "uM"),
        (None, QUERY_UNIT_NONE),
        ("", QUERY_UNIT_NONE),
        ("  ", QUERY_UNIT_NONE),
        ("unitless", QUERY_UNIT_NONE),   # users don't type "unitless" — treat as absent
        ("bananas", QUERY_UNIT_UNKNOWN),
    ],
)
def test_classify_query_unit(unit, expected):
    assert classify_query_unit(unit) == expected


# ===========================================================================
# §6.3 — metric + status checks (filter reasons, not footer-noted)
# ===========================================================================


def test_valid_row_matching_metric_with_no_threshold_matches():
    out = evaluate_row(row(value=5.0), metric="IC50", threshold=None)
    assert isinstance(out, RowMatched)
    assert out.value == 5.0
    assert out.row_unit == "nM"
    assert out.value_in_query_unit == 5.0


def test_invalid_status_filters():
    out = evaluate_row(row(status="invalid", value=5.0), metric="IC50", threshold=None)
    assert isinstance(out, RowFiltered)
    assert out.reason == FILTER_INVALID_STATUS


def test_unassigned_status_filters():
    out = evaluate_row(row(status="unassigned", value=5.0), metric="IC50", threshold=None)
    assert isinstance(out, RowFiltered)
    assert out.reason == FILTER_INVALID_STATUS


def test_kpi_mismatch_filters():
    # Row targets EC50 but user asked for IC50 — even if the row happens to
    # carry an IC50 key, it wasn't the fit's primary KPI.
    r = SimpleNamespace(
        status="valid",
        results={"KPI": "EC50", "EC50": 12.0, "IC50": 5.0, "kpi_unit": "nM"},
    )
    out = evaluate_row(r, metric="IC50", threshold=None)
    assert isinstance(out, RowFiltered)
    assert out.reason == FILTER_KPI_MISMATCH


def test_value_not_numeric_filters():
    r = SimpleNamespace(status="valid", results={"KPI": "IC50", "IC50": "N/A", "kpi_unit": "nM"})
    out = evaluate_row(r, metric="IC50", threshold=None)
    assert isinstance(out, RowFiltered)
    assert out.reason == FILTER_VALUE_NOT_NUMERIC


def test_value_is_bool_filters():
    r = SimpleNamespace(status="valid", results={"KPI": "IC50", "IC50": True, "kpi_unit": "nM"})
    out = evaluate_row(r, metric="IC50", threshold=None)
    assert isinstance(out, RowFiltered)
    assert out.reason == FILTER_VALUE_NOT_NUMERIC


def test_value_is_nan_filters():
    r = SimpleNamespace(
        status="valid",
        results={"KPI": "IC50", "IC50": float("nan"), "kpi_unit": "nM"},
    )
    out = evaluate_row(r, metric="IC50", threshold=None)
    assert isinstance(out, RowFiltered)
    assert out.reason == FILTER_VALUE_NOT_NUMERIC


def test_value_missing_filters():
    r = SimpleNamespace(status="valid", results={"KPI": "IC50", "kpi_unit": "nM"})
    out = evaluate_row(r, metric="IC50", threshold=None)
    assert isinstance(out, RowFiltered)
    assert out.reason == FILTER_VALUE_NOT_NUMERIC


# ===========================================================================
# §6.4 — threshold comparison + unit conversion (happy paths)
# ===========================================================================


def test_same_unit_threshold_pass():
    out = evaluate_row(
        row(value=5.0, kpi_unit="nM"),
        metric="IC50",
        threshold=Threshold(op="<", value=10.0, unit="nM"),
    )
    assert isinstance(out, RowMatched)
    assert out.value == 5.0
    assert out.row_unit == "nM"
    assert out.value_in_query_unit == 5.0


def test_cross_unit_threshold_converts_and_passes():
    # Row IC50 = 5 nM; query IC50 < 10 uM → 5 nM = 0.005 uM, passes.
    out = evaluate_row(
        row(value=5.0, kpi_unit="nM"),
        metric="IC50",
        threshold=Threshold(op="<", value=10.0, unit="uM"),
    )
    assert isinstance(out, RowMatched)
    assert out.value == 5.0                       # native-unit value preserved
    assert out.row_unit == "nM"
    assert math.isclose(out.value_in_query_unit, 0.005, rel_tol=1e-12)


def test_cross_unit_threshold_converts_and_fails():
    # Row IC50 = 10 uM; query IC50 < 10 nM → 10 uM = 10000 nM, fails.
    out = evaluate_row(
        row(value=10.0, kpi_unit="uM"),
        metric="IC50",
        threshold=Threshold(op="<", value=10.0, unit="nM"),
    )
    assert isinstance(out, RowFiltered)
    assert out.reason == FILTER_THRESHOLD_NOT_MET


def test_unicode_micro_sign_converts_cleanly():
    # Row stores kpi_unit as unicode 'µM'; query uses ASCII 'uM'. After
    # normalize_unit both collapse to 'uM' — no EXCLUDE result.
    out = evaluate_row(
        row(value=1.0, kpi_unit="µM"),
        metric="IC50",
        threshold=Threshold(op="<", value=10.0, unit="uM"),
    )
    assert isinstance(out, RowMatched)
    assert out.row_unit == "uM"


@pytest.mark.parametrize(
    "op,value,threshold_val,expected",
    [
        ("<", 5.0, 10.0, RowMatched),
        ("<", 10.0, 10.0, RowFiltered),        # strict
        ("<=", 10.0, 10.0, RowMatched),
        (">", 15.0, 10.0, RowMatched),
        (">", 5.0, 10.0, RowFiltered),
        (">=", 10.0, 10.0, RowMatched),
        ("=", 10.0, 10.0, RowMatched),
        ("==", 10.0, 10.0, RowMatched),
        ("!=", 10.0, 10.0, RowFiltered),
        ("!=", 5.0, 10.0, RowMatched),
    ],
)
def test_comparison_operators(op, value, threshold_val, expected):
    out = evaluate_row(
        row(value=value, kpi_unit="nM"),
        metric="IC50",
        threshold=Threshold(op=op, value=threshold_val, unit="nM"),
    )
    assert isinstance(out, expected)


# ===========================================================================
# §6.4 — unit tri-state (unitless / unknown / real)
# ===========================================================================


def test_unitless_row_matches_unitless_query():
    # Caco-2 efflux ratio stored as kpi_unit='unitless'; query has no unit.
    out = evaluate_row(
        row(kpi="Efflux", value=0.8, kpi_unit="unitless"),
        metric="Efflux",
        threshold=Threshold(op="<", value=2.0, unit=None),
    )
    assert isinstance(out, RowMatched)
    assert out.row_unit is None
    assert out.value_in_query_unit == 0.8


def test_unitless_row_with_backward_compat_none_matches_unitless_query():
    out = evaluate_row(
        row(kpi="Efflux", value=0.8, kpi_unit=None),
        metric="Efflux",
        threshold=Threshold(op="<", value=2.0, unit=None),
    )
    assert isinstance(out, RowMatched)
    assert out.row_unit is None


def test_unitless_row_vs_real_unit_query_is_type_mismatch():
    out = evaluate_row(
        row(kpi="Efflux", value=0.8, kpi_unit="unitless"),
        metric="Efflux",
        threshold=Threshold(op="<", value=2.0, unit="nM"),
    )
    assert isinstance(out, RowExcluded)
    assert out.reason == EXCLUDE_UNIT_TYPE_MISMATCH


def test_real_unit_row_vs_unitless_query_is_query_missing_unit():
    # Query omits unit but row has a real unit → ambiguous, footer-excluded.
    out = evaluate_row(
        row(value=5.0, kpi_unit="nM"),
        metric="IC50",
        threshold=Threshold(op="<", value=10.0, unit=None),
    )
    assert isinstance(out, RowExcluded)
    assert out.reason == EXCLUDE_QUERY_MISSING_UNIT


def test_unknown_row_unit_is_excluded():
    # results['kpi_unit'] is empty string — the "missing" state from §6.4.
    r = SimpleNamespace(
        status="valid",
        results={"KPI": "IC50", "IC50": 5.0, "kpi_unit": ""},
    )
    out = evaluate_row(
        r,
        metric="IC50",
        threshold=Threshold(op="<", value=10.0, unit="nM"),
    )
    assert isinstance(out, RowExcluded)
    assert out.reason == EXCLUDE_UNIT_UNKNOWN


def test_absent_kpi_unit_key_is_excluded():
    # No kpi_unit key at all → also the "unknown" state.
    r = SimpleNamespace(status="valid", results={"KPI": "IC50", "IC50": 5.0})
    out = evaluate_row(
        r,
        metric="IC50",
        threshold=Threshold(op="<", value=10.0, unit="nM"),
    )
    # Note: an absent kpi_unit key is equivalent to is_unitless(None)==True
    # per kpi_utils' backward-compat rule. So this row classifies as UNITLESS
    # under a threshold query, and unitless × real-unit query = type mismatch.
    assert isinstance(out, RowExcluded)
    assert out.reason == EXCLUDE_UNIT_TYPE_MISMATCH


def test_unrecognised_row_unit_is_excluded():
    out = evaluate_row(
        row(value=5.0, kpi_unit="bananas"),
        metric="IC50",
        threshold=Threshold(op="<", value=10.0, unit="nM"),
    )
    assert isinstance(out, RowExcluded)
    assert out.reason == EXCLUDE_UNIT_UNKNOWN


def test_unrecognised_query_unit_is_excluded():
    out = evaluate_row(
        row(value=5.0, kpi_unit="nM"),
        metric="IC50",
        threshold=Threshold(op="<", value=10.0, unit="bananas"),
    )
    assert isinstance(out, RowExcluded)
    assert out.reason == EXCLUDE_UNIT_UNKNOWN


def test_cross_family_units_are_incompatible():
    out = evaluate_row(
        row(kpi="t_half", value=30.0, kpi_unit="min"),
        metric="t_half",
        threshold=Threshold(op="<", value=10.0, unit="nM"),
    )
    assert isinstance(out, RowExcluded)
    assert out.reason == EXCLUDE_UNIT_INCOMPATIBLE


def test_different_clearance_families_are_incompatible():
    out = evaluate_row(
        row(kpi="CLint", value=10.0, kpi_unit="uL/min/mg"),
        metric="CLint",
        threshold=Threshold(op="<", value=5.0, unit="mL/min/kg"),
    )
    assert isinstance(out, RowExcluded)
    assert out.reason == EXCLUDE_UNIT_INCOMPATIBLE


# ===========================================================================
# No-threshold path — any valid, metric-matching row passes regardless of unit
# ===========================================================================


def test_no_threshold_unknown_unit_still_matches():
    # Without a threshold there's no comparison to do, so even an unknown
    # unit shouldn't footer-exclude; payload just reports row_unit=None.
    r = SimpleNamespace(
        status="valid",
        results={"KPI": "IC50", "IC50": 5.0, "kpi_unit": ""},
    )
    out = evaluate_row(r, metric="IC50", threshold=None)
    assert isinstance(out, RowMatched)
    assert out.row_unit is None


def test_no_threshold_integer_value_is_coerced_to_float():
    out = evaluate_row(row(value=5, kpi_unit="nM"), metric="IC50", threshold=None)
    assert isinstance(out, RowMatched)
    assert isinstance(out.value, float)
    assert out.value == 5.0


# ===========================================================================
# Operator validation
# ===========================================================================


def test_unsupported_operator_raises():
    with pytest.raises(ValueError):
        evaluate_row(
            row(value=5.0, kpi_unit="nM"),
            metric="IC50",
            threshold=Threshold(op="LIKE", value=10.0, unit="nM"),
        )
