"""Server-side scorecard score computation for "best X compounds" ranking.

Mirrors the logic in ``apps/compounds/frontend/lib/compounds/scorecard.ts``
exactly — same per-axis value extraction, same target/poor normalisation,
same log-scale default for protocol/ratio/worst_of axes — so a compound
ranked "best" by the NLP executor lines up with the same compound's
position on the project's spider/radar in the aggregation page UI.

Used by the NLP executor when ``CompoundSelector.rank_by == "scorecard"``.
The catalogue authority is the Target's ``scorecard_config`` JSON; this
module is a pure consumer.

Returns a dict ``{compound_pk: mean_t_score}`` where ``mean_t_score`` is
the arithmetic mean of non-null per-axis t-values (each ∈ [0, 1]) for
that compound. Compounds whose every axis is null are dropped — there's
no signal to rank them on.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional


# ---------------------------------------------------------------------------
# Per-axis value extraction
# ---------------------------------------------------------------------------


def _as_finite(value: Any) -> Optional[float]:
    """Coerce number-ish inputs (Decimal strings, None, NaN) to a finite
    float or None — same defensive coercion as the frontend's
    ``asFiniteNumber``."""
    if value is None:
        return None
    try:
        n = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(n):
        return None
    return n


def _axis_value_protocol(axis: dict, row: dict) -> Optional[float]:
    protocols = row.get("protocols") or {}
    pdata = protocols.get(axis.get("protocol_id"))
    if not pdata:
        return None
    return _as_finite(pdata.get("geomean"))


def _axis_value_ratio(axis: dict, row: dict) -> Optional[float]:
    protocols = row.get("protocols") or {}
    num = _as_finite((protocols.get(axis.get("numerator_id")) or {}).get("geomean"))
    den = _as_finite((protocols.get(axis.get("denominator_id")) or {}).get("geomean"))
    if num is None or den is None or den == 0:
        return None
    return num / den


def _axis_value_worst_of(axis: dict, row: dict) -> Optional[float]:
    protocols = row.get("protocols") or {}
    ids = axis.get("protocol_ids") or []
    if not isinstance(ids, list) or not ids:
        return None
    values: List[float] = []
    for pid in ids:
        v = _as_finite((protocols.get(pid) or {}).get("geomean"))
        if v is not None:
            values.append(v)
    if not values:
        return None
    target = _as_finite(axis.get("target_value"))
    poor = _as_finite(axis.get("poor_value"))
    if target is None or poor is None:
        return None
    # When target < poor (lower-is-better), the worst observed value is
    # the maximum; otherwise the minimum. Same convention as the
    # frontend's worst_of branch.
    return max(values) if target < poor else min(values)


def _axis_value_lipinski(axis: dict, row: dict) -> Optional[float]:
    """Count of passed Ro5 rules (0–4). Pre-defined thresholds match
    the frontend's hard-coded Lipinski check."""
    props = row.get("properties") or {}
    mw = _as_finite(props.get("molecular_weight"))
    clogp = _as_finite(props.get("clogp"))
    hbd = _as_finite(props.get("hbd"))
    hba = _as_finite(props.get("hba"))
    pass_count = 0
    if mw is not None and mw <= 500:
        pass_count += 1
    if clogp is not None and clogp <= 5:
        pass_count += 1
    if hbd is not None and hbd <= 5:
        pass_count += 1
    if hba is not None and hba <= 10:
        pass_count += 1
    return float(pass_count)


_AXIS_DISPATCH = {
    "protocol": _axis_value_protocol,
    "ratio": _axis_value_ratio,
    "worst_of": _axis_value_worst_of,
    "lipinski": _axis_value_lipinski,
}


def _axis_value(axis: dict, row: dict) -> Optional[float]:
    fn = _AXIS_DISPATCH.get(axis.get("kind"))
    if fn is None:
        return None
    try:
        return fn(axis, row)
    except Exception:
        # Fail-soft on unexpected config shapes — same discipline as the
        # frontend: one bad axis doesn't crash the whole ranking.
        return None


# ---------------------------------------------------------------------------
# Per-axis normalisation (raw value → 0–1 t-score)
# ---------------------------------------------------------------------------


def _normalise_axis(axis: dict, value: float) -> Optional[float]:
    """Linear or log normalisation between (poor_value → 0) and
    (target_value → 1), clamped to [0, 1]. Mirrors the frontend's
    ``normaliseAxis``."""
    target = _as_finite(axis.get("target_value"))
    poor = _as_finite(axis.get("poor_value"))
    if target is None or poor is None or target == poor:
        return None

    scale = axis.get("threshold_scale") or (
        "linear" if axis.get("kind") == "lipinski" else "log"
    )
    can_log = (
        scale == "log" and value > 0 and target > 0 and poor > 0
    )
    v = math.log10(value) if can_log else value
    t = math.log10(target) if can_log else target
    p = math.log10(poor) if can_log else poor
    raw = (v - p) / (t - p)
    return max(0.0, min(1.0, raw))


# ---------------------------------------------------------------------------
# Per-compound mean t-score
# ---------------------------------------------------------------------------


def _row_mean_t(config: dict, row: dict) -> Optional[float]:
    """Arithmetic mean of non-null per-axis t-values for one compound's
    aggregated row. None when no axis produces a t-value (no signal to
    rank on — drop the compound from the ranking)."""
    axes = (config or {}).get("axes") or []
    ts: List[float] = []
    for axis in axes:
        v = _axis_value(axis, row)
        if v is None:
            continue
        t = _normalise_axis(axis, v)
        if t is not None:
            ts.append(t)
    if not ts:
        return None
    return sum(ts) / len(ts)


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------


def scorecard_data_needs(config: dict) -> Dict[str, Any]:
    """Return the protocol_ids and molecular-property names a scorecard
    config depends on — needed to drive the aggregation request that
    feeds the ranker. Mirrors frontend ``scorecardDataNeeds``."""
    protocol_ids: List[str] = []
    seen: set = set()
    needs_lipinski = False
    for axis in (config or {}).get("axes") or []:
        kind = axis.get("kind")
        if kind == "protocol":
            pid = axis.get("protocol_id")
            if pid and pid not in seen:
                seen.add(pid)
                protocol_ids.append(pid)
        elif kind == "ratio":
            for pid in (axis.get("numerator_id"), axis.get("denominator_id")):
                if pid and pid not in seen:
                    seen.add(pid)
                    protocol_ids.append(pid)
        elif kind == "worst_of":
            for pid in axis.get("protocol_ids") or []:
                if pid and pid not in seen:
                    seen.add(pid)
                    protocol_ids.append(pid)
        elif kind == "lipinski":
            needs_lipinski = True

    properties: List[str] = []
    if needs_lipinski:
        properties = ["molecular_weight", "clogp", "hbd", "hba"]
    return {"protocol_ids": protocol_ids, "properties": properties}


def rank_compounds_by_scorecard(
    target,
    compound_pks: set,
) -> Dict[Any, float]:
    """Compute mean-t scorecard scores for the given compounds against
    the target's scorecard config.

    Returns ``{compound_pk: mean_t}`` for compounds that produced a
    score — compounds where every axis is null are simply absent from
    the dict (the caller filters them out of the ranked top-N).

    Raises ``ValueError`` if the target has no scorecard config — the
    caller (executor) maps this to a SpecError that the chemist sees as
    a structured error response.
    """
    config = target.scorecard_config
    if not config or not config.get("axes"):
        raise ValueError(
            f"Target {target.name!r} has no scorecard configured; "
            "cannot rank by scorecard."
        )

    needs = scorecard_data_needs(config)

    # Lazy import: the aggregation module pulls in the assays graph;
    # we don't want to load it until ranking is actually requested.
    from compounds.assays.aggregation import aggregate_compact_from_compounds
    from compounds.registry.models import Compound

    if not compound_pks:
        return {}

    # The aggregator iterates ``compound.filtered_data_series`` — set up
    # via a Prefetch with ``to_attr``. Restrict to valid analyses on the
    # protocols the scorecard actually needs; everything else is noise.
    from django.db.models import Prefetch
    from compounds.assays.models import DataSeries
    valid_ds = DataSeries.objects.filter(
        analysis__status="valid",
        assay__protocol_id__in=needs["protocol_ids"] or [],
    ).select_related("assay__protocol", "analysis")
    compound_qs = (
        Compound.objects
        .filter(pk__in=compound_pks)
        .prefetch_related(
            Prefetch(
                "assay_results", queryset=valid_ds,
                to_attr="filtered_data_series",
            ),
        )
    )

    payload = aggregate_compact_from_compounds(
        compound_qs,
        protocol_ids=needs["protocol_ids"],
        aggregations=["geomean"],
        include_properties=needs["properties"],
        include_identifiers=False,
    )

    scores: Dict[Any, float] = {}
    for row in payload.get("data") or []:
        # The compact row's "compound_id" is the pk; double-check vs the
        # aggregator's actual key.
        pk = row.get("compound_id") or row.get("compound", {}).get("id")
        if pk is None:
            continue
        score = _row_mean_t(config, row)
        if score is None:
            continue
        scores[pk] = score
    return scores
