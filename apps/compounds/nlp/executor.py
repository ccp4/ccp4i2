"""Compound-selection execution (§7 executor after the pivot).

Runs a fully-formed CompoundSelector against the database and produces
a list of compound IDs that satisfy **all** measurement filters (AND).

Composes earlier-slice modules:
    resolve_targets  →  resolve_protocol (once per filter) → scoped
    AnalysisResult walk through evaluate_row → intersect compound sets

Returns one of:
- CompoundSelection       — compound IDs + resolved target + filter-used
                            protocols + scope sentence
- TargetClarify / TargetMiss / ScopeError   — passed through from resolve_targets
- ProtocolClarify / ProtocolMiss            — passed through from resolve_protocol,
                                               tagged with the filter index so
                                               the UI knows which filter to pin
- SpecError               — filter with threshold but no metric, etc.

The aggregation page (`/assays/aggregate`) is the display surface — the
executor emits only enough to construct the redirect URL there. Display
format, columns, aggregation mode etc. all live downstream.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Set

from django.db.models import QuerySet

from compounds.assays.models import AnalysisResult
from compounds.registry.models import Compound, Target

from .evaluator import evaluate_row
from .resolver import resolve_protocol, resolve_targets
from .spec import (
    SCOPE_ASSAY_ONLY,
    SCOPE_BOTH_SAME,
    SCOPE_CROSS,
    SCOPE_REG_ONLY,
    CompoundSelection,
    CompoundSelector,
    ExecutionResult,
    MeasurementFilter,
    ProtocolClarify,
    ProtocolMiss,
    ResolvedProtocol,
    ResolvedTargets,
    RowMatched,
    ScopeError,
    SpecError,
    TargetClarify,
    TargetMiss,
    Threshold,
)


def execute(selector: CompoundSelector) -> ExecutionResult:
    """Run a fully-formed CompoundSelector and return a CompoundSelection
    or a clarify/miss/error response to surface back to the user."""
    # Validate filters before we touch the DB — a threshold without a
    # metric is meaningless and caught here rather than silently mis-
    # filtering downstream.
    for idx, flt in enumerate(selector.measurement_filters):
        if flt.threshold is not None and not flt.metric:
            return SpecError(
                field=f"measurement_filters[{idx}].metric",
                message="A threshold requires a metric to compare against.",
            )

    rt = resolve_targets(selector)
    if not isinstance(rt, ResolvedTargets):
        return rt

    # Resolve each filter's protocol hint (if any) in order. First
    # ambiguity short-circuits — the UI clarifies one filter at a time.
    resolved_protocols: List[Optional[ResolvedProtocol]] = []
    for idx, flt in enumerate(selector.measurement_filters):
        if flt.protocol_hint is None and flt.protocol_id is None:
            # No protocol constraint — filter matches any protocol in scope.
            resolved_protocols.append(None)
            continue
        rp = resolve_protocol(
            flt.protocol_hint or "", rt,
            pinned_id=flt.protocol_id,
            filter_index=idx,
        )
        if not isinstance(rp, ResolvedProtocol):
            return rp  # ProtocolClarify / ProtocolMiss already tagged with filter_index
        resolved_protocols.append(rp)

    # Apply each filter and intersect the surviving compound sets.
    base_compound_ids = _base_scope_compound_ids(rt)
    surviving: Optional[Set[str]] = None

    for flt, rp in zip(selector.measurement_filters, resolved_protocols):
        matching = _filter_matching_compounds(rt, rp, flt)
        surviving = matching if surviving is None else surviving & matching

    if surviving is None:
        # No filters — every compound in the base scope is in the selection.
        compound_ids = base_compound_ids
    else:
        compound_ids = surviving & base_compound_ids

    # Hydrate compound IDs into formatted identifiers, ordered by
    # reg_number for stable URLs + predictable UI.
    formatted = list(
        Compound.objects
        .filter(pk__in=compound_ids)
        .order_by("reg_number")
        .values_list("reg_number", flat=True)
    )
    formatted_ids = [_format_compound_id(rn) for rn in formatted]

    target_names = _target_names(rt)
    protocol_names = [rp.protocol.name for rp in resolved_protocols if rp is not None]

    return CompoundSelection(
        compound_formatted_ids=formatted_ids,
        target_names=target_names,
        protocol_names=protocol_names,
        n_matched=len(formatted_ids),
        scope_sentence=_scope_sentence(rt, selector, resolved_protocols),
    )


# ---------------------------------------------------------------------------
# Scope construction
# ---------------------------------------------------------------------------


def _base_scope_compound_ids(rt: ResolvedTargets) -> Set[str]:
    """Compound IDs that pass the target scope (before any filters).

    Each scope_kind maps to a different base query:
    - SCOPE_REG_ONLY / SCOPE_BOTH_SAME / SCOPE_CROSS: compounds registered
      to the registration target.
    - SCOPE_ASSAY_ONLY: compounds that have at least one DataSeries under
      an Assay with the given target.
    """
    if rt.registration is not None:
        qs = Compound.objects.filter(target=rt.registration)
    else:
        # assay-only scope
        qs = Compound.objects.filter(assay_results__assay__target=rt.assay).distinct()

    if rt.scope_kind == SCOPE_CROSS and rt.assay is not None:
        # Registration target filter above; also require the compound to
        # have at least one DataSeries against the assay target so we don't
        # surface "my ARd compounds" with no actual AKT data.
        qs = qs.filter(assay_results__assay__target=rt.assay).distinct()

    return set(qs.values_list("pk", flat=True))


def _filter_matching_compounds(
    rt: ResolvedTargets,
    rp: Optional[ResolvedProtocol],
    flt: MeasurementFilter,
) -> Set[str]:
    """Compound IDs whose AnalysisResult rows satisfy this filter.

    A row satisfies the filter iff evaluate_row returns RowMatched for it.
    The row itself is per-measurement; a compound is in the set if *any*
    of its rows pass (compound-level OR within a filter).
    """
    qs = _scoped_analysis_results(rt)
    if rp is not None:
        qs = qs.filter(data_series__assay__protocol=rp.protocol)

    matching: Set[str] = set()
    for ar in qs.iterator():
        outcome = evaluate_row(ar, flt.metric, flt.threshold)
        if isinstance(outcome, RowMatched):
            ds = getattr(ar, "data_series", None)
            if ds is None or ds.compound_id is None:
                continue
            matching.add(ds.compound_id)
    return matching


def _scoped_analysis_results(rt: ResolvedTargets) -> QuerySet:
    """AnalysisResult queryset filtered to the target scope (registration
    and/or assay target). Status is pre-filtered to 'valid' — invalid rows
    carry no selection value."""
    qs = AnalysisResult.objects.filter(status="valid")
    if rt.registration is not None:
        qs = qs.filter(data_series__compound__target=rt.registration)
    if rt.assay is not None:
        qs = qs.filter(data_series__assay__target=rt.assay)
    return qs.select_related("data_series__compound", "data_series__assay__protocol")


# ---------------------------------------------------------------------------
# Presentation helpers
# ---------------------------------------------------------------------------


def _format_compound_id(reg_number: int) -> str:
    """Mirror the logic in registry/models.Compound.formatted_id —
    uses the configured COMPOUND_ID_PREFIX and zero-padding.

    Imported lazily via the Compound class since the values_list avoids
    the instance hit."""
    from django.conf import settings
    prefix = getattr(settings, "COMPOUND_ID_PREFIX", "NCL")
    digits = getattr(settings, "COMPOUND_ID_DIGITS", 8)
    return f"{prefix}-{reg_number:0{digits}d}"


def _target_names(rt: ResolvedTargets) -> List[str]:
    """Names to pass into the aggregation page's `targets=` URL param.

    In §6.5 terms:
    - SCOPE_BOTH_SAME / SCOPE_REG_ONLY / SCOPE_ASSAY_ONLY → single target
    - SCOPE_CROSS → both targets (UI can pick whichever makes sense)
    """
    if rt.scope_kind == SCOPE_CROSS:
        return [t.name for t in (rt.registration, rt.assay) if t is not None]
    target = rt.registration or rt.assay
    return [target.name] if target is not None else []


def _scope_sentence(
    rt: ResolvedTargets,
    selector: CompoundSelector,
    resolved_protocols: List[Optional[ResolvedProtocol]],
) -> str:
    """Human echo-back (§10). Two clauses: the target scope, and the
    AND-joined filter list (if any)."""
    reg = rt.registration.name if rt.registration is not None else ""
    assay = rt.assay.name if rt.assay is not None else ""
    if rt.scope_kind == SCOPE_BOTH_SAME:
        base = f"Showing compounds registered to {reg} and tested against {reg}"
    elif rt.scope_kind == SCOPE_CROSS:
        base = f"Showing compounds registered to {reg}, tested against {assay}"
    elif rt.scope_kind == SCOPE_REG_ONLY:
        base = f"Showing compounds registered to {reg}"
    elif rt.scope_kind == SCOPE_ASSAY_ONLY:
        base = f"Showing compounds tested against {assay}"
    else:
        base = "Showing compounds"

    filter_clauses: List[str] = []
    for flt, rp in zip(selector.measurement_filters, resolved_protocols):
        clause = _filter_clause(flt, rp)
        if clause:
            filter_clauses.append(clause)
    if filter_clauses:
        base += " where " + " AND ".join(filter_clauses)
    return base


def _filter_clause(
    flt: MeasurementFilter,
    rp: Optional[ResolvedProtocol],
) -> str:
    """Render one MeasurementFilter as human prose for the scope sentence."""
    protocol_name = rp.protocol.name if rp is not None else None
    parts: List[str] = []
    if protocol_name:
        parts.append(protocol_name)
    if flt.metric:
        parts.append(flt.metric)
    if flt.threshold is not None:
        t: Threshold = flt.threshold
        unit = f" {t.unit}" if t.unit else ""
        parts.append(f"{t.op} {t.value}{unit}")
    return " ".join(parts)
