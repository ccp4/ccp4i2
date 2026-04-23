"""Row-level evaluation for threshold queries.

Implements §6.3 (metric match: results['KPI']==metric, results[metric] numeric,
status='valid') and §6.4 (kpi_unit tri-state per Q25: known real unit /
"unitless" sentinel / absent-or-empty unknown). Threshold conversion is
per-row via the unit family tables below.

Pure Python over an AnalysisResult-shaped object. No ORM query here — the
future executor will walk the scoped queryset and call `evaluate_row()` on
each row to classify it into matched / filtered / excluded buckets.
"""

from __future__ import annotations

import math
from typing import Dict, Optional, Tuple

from compounds.assays.kpi_utils import UNITLESS, VALID_UNITS, is_unitless, normalize_unit

from .spec import (
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
    RowOutcome,
    Threshold,
)


# Row unit-classification sentinels (distinct from the stored kpi_unit values).
ROW_UNIT_UNITLESS = "unitless"
ROW_UNIT_UNKNOWN = "unknown"

# Query unit classifications.
QUERY_UNIT_NONE = "none"
QUERY_UNIT_UNKNOWN = "unknown"


# Unit families and the multiplier to the family's canonical unit. Conversion
# from unit A to unit B within the same family is `value * (A_factor / B_factor)`.
# Cross-family conversion is undefined (returns None).
_UNIT_FAMILIES: Dict[str, Tuple[str, float]] = {
    # Concentration (canonical: M)
    "pM": ("concentration", 1e-12),
    "nM": ("concentration", 1e-9),
    "uM": ("concentration", 1e-6),
    "mM": ("concentration", 1e-3),
    "M": ("concentration", 1.0),
    # Time (canonical: seconds)
    "s": ("time", 1.0),
    "min": ("time", 60.0),
    "h": ("time", 3600.0),
    # Percentage (no conversion to/from anything else)
    "%": ("percent", 1.0),
    # Permeability (canonical: cm/s)
    "1e-6 cm/s": ("permeability", 1e-6),
    "cm/s": ("permeability", 1.0),
    # Clearance families each stand alone — the denominators don't compose.
    "uL/min/mg": ("clearance_microsomal", 1.0),
    "mL/min/kg": ("clearance_hepatic", 1.0),
}


def convert(value: float, from_unit: str, to_unit: str) -> Optional[float]:
    """Convert a value between units in the same family. None if incompatible."""
    if from_unit == to_unit:
        return value
    a = _UNIT_FAMILIES.get(from_unit)
    b = _UNIT_FAMILIES.get(to_unit)
    if a is None or b is None or a[0] != b[0]:
        return None
    return value * (a[1] / b[1])


def classify_row_unit(kpi_unit) -> str:
    """Return a normalized real-unit string, ROW_UNIT_UNITLESS, or ROW_UNIT_UNKNOWN.

    Implements the Q25 tri-state:
    - `"unitless"` (or `None` for backward compat) → ROW_UNIT_UNITLESS
    - known real unit (after normalize_unit) → that normalized string
    - absent key, empty string, or unrecognised unit → ROW_UNIT_UNKNOWN
    """
    if is_unitless(kpi_unit):
        return ROW_UNIT_UNITLESS
    if not isinstance(kpi_unit, str) or not kpi_unit.strip():
        return ROW_UNIT_UNKNOWN
    normalized = normalize_unit(kpi_unit.strip())
    if normalized in VALID_UNITS and normalized != UNITLESS:
        return normalized
    return ROW_UNIT_UNKNOWN


def classify_query_unit(unit) -> str:
    """Return a normalized real-unit string, QUERY_UNIT_NONE, or QUERY_UNIT_UNKNOWN.

    A query unit of `"unitless"` is treated as QUERY_UNIT_NONE — users don't
    type "unitless" in natural language; absence means dimensionless-or-raw.
    """
    if unit is None:
        return QUERY_UNIT_NONE
    if not isinstance(unit, str):
        return QUERY_UNIT_UNKNOWN
    stripped = unit.strip()
    if not stripped:
        return QUERY_UNIT_NONE
    normalized = normalize_unit(stripped)
    if normalized == UNITLESS:
        return QUERY_UNIT_NONE
    if normalized in VALID_UNITS:
        return normalized
    return QUERY_UNIT_UNKNOWN


def _compare(value: float, op: str, threshold: float) -> bool:
    if op == "<":
        return value < threshold
    if op == "<=":
        return value <= threshold
    if op == ">":
        return value > threshold
    if op == ">=":
        return value >= threshold
    if op == "=" or op == "==":
        return value == threshold
    if op == "!=":
        return value != threshold
    raise ValueError(f"Unsupported comparison operator: {op!r}")


def _is_numeric(value) -> bool:
    # bool is a subclass of int — filter it out so True/False can't sneak past.
    if isinstance(value, bool):
        return False
    if isinstance(value, (int, float)):
        return not (isinstance(value, float) and math.isnan(value))
    return False


def evaluate_row(row, metric: str, threshold: Optional[Threshold]) -> RowOutcome:
    """Classify one AnalysisResult row against a metric + optional threshold.

    `row` is duck-typed — anything with `.status` and `.results` works, so
    tests can pass lightweight stand-ins when that's simpler than creating
    real AnalysisResult rows.
    """
    if row.status != "valid":
        return RowFiltered(reason=FILTER_INVALID_STATUS)

    results = row.results or {}
    if results.get("KPI") != metric:
        return RowFiltered(reason=FILTER_KPI_MISMATCH)

    value = results.get(metric)
    if not _is_numeric(value):
        return RowFiltered(reason=FILTER_VALUE_NOT_NUMERIC)
    value = float(value)

    row_unit_class = classify_row_unit(results.get("kpi_unit"))

    # No threshold → the caller just wants to know the row matches the metric.
    if threshold is None:
        row_unit_for_payload = None if row_unit_class == ROW_UNIT_UNITLESS else (
            None if row_unit_class == ROW_UNIT_UNKNOWN else row_unit_class
        )
        return RowMatched(
            value=value,
            row_unit=row_unit_for_payload,
            value_in_query_unit=value,
        )

    query_unit_class = classify_query_unit(threshold.unit)

    # Row-unit-unknown is always excluded under a threshold query — we
    # cannot safely compare numeric magnitudes of unknown scale.
    if row_unit_class == ROW_UNIT_UNKNOWN:
        return RowExcluded(reason=EXCLUDE_UNIT_UNKNOWN)

    if query_unit_class == QUERY_UNIT_UNKNOWN:
        # The LLM emitted something we don't recognise. Treat symmetrically
        # to unknown row unit — excluded at the row level.
        return RowExcluded(reason=EXCLUDE_UNIT_UNKNOWN)

    # Row is unitless (dimensionless KPI).
    if row_unit_class == ROW_UNIT_UNITLESS:
        if query_unit_class != QUERY_UNIT_NONE:
            return RowExcluded(reason=EXCLUDE_UNIT_TYPE_MISMATCH)
        # Query unit absent → compare raw.
        if _compare(value, threshold.op, threshold.value):
            return RowMatched(value=value, row_unit=None, value_in_query_unit=value)
        return RowFiltered(reason=FILTER_THRESHOLD_NOT_MET)

    # Row has a known real unit.
    if query_unit_class == QUERY_UNIT_NONE:
        # Query lacks a unit but row has one — ambiguous; footer it so the
        # user can correct by re-typing with a unit.
        return RowExcluded(reason=EXCLUDE_QUERY_MISSING_UNIT)

    # Both sides have real units. Try to convert the threshold into the row's
    # unit. Cross-family (e.g. nM ↔ min) fails.
    converted_threshold = convert(threshold.value, query_unit_class, row_unit_class)
    if converted_threshold is None:
        return RowExcluded(reason=EXCLUDE_UNIT_INCOMPATIBLE)

    # Same-family → also express the row value in the query's unit for the
    # payload, so the caller can render a consistent column regardless of
    # which per-row unit was stored.
    value_in_query = convert(value, row_unit_class, query_unit_class)
    assert value_in_query is not None  # same family → cannot fail

    if _compare(value, threshold.op, converted_threshold):
        return RowMatched(
            value=value,
            row_unit=row_unit_class,
            value_in_query_unit=value_in_query,
        )
    return RowFiltered(reason=FILTER_THRESHOLD_NOT_MET)
