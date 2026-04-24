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
from .resolver import (
    assay_creators_qs,
    compound_registrars_qs,
    resolve_protocol,
    resolve_targets,
    resolve_user,
)
from .spec import (
    FIELD_ASSAYED_BY,
    FIELD_REGISTERED_BY,
    SCOPE_ASSAY_ONLY,
    SCOPE_BOTH_SAME,
    SCOPE_CROSS,
    SCOPE_REG_ONLY,
    CompoundSelection,
    CompoundSelector,
    DateRange,
    ExecutionResult,
    MeasurementFilter,
    ProtocolClarify,
    ProtocolMiss,
    ResolvedProtocol,
    ResolvedTargets,
    ResolvedUser,
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

    # Selector-level user resolution (registered-by). Sequential clarify —
    # any unresolved user short-circuits the whole selector.
    resolved_registered_by = None
    if selector.registered_by_as_typed or selector.registered_by_id:
        ru = resolve_user(
            selector.registered_by_as_typed or "",
            compound_registrars_qs(),
            field=FIELD_REGISTERED_BY,
            pinned_id=selector.registered_by_id,
        )
        if not isinstance(ru, ResolvedUser):
            return ru
        resolved_registered_by = ru

    # Resolve each filter's protocol hint + assayed-by user (if any) in
    # order. First ambiguity short-circuits — the UI clarifies one
    # filter at a time.
    resolved_protocols: List[Optional[ResolvedProtocol]] = []
    resolved_assayed_bys: List[Optional[ResolvedUser]] = []
    for idx, flt in enumerate(selector.measurement_filters):
        if flt.protocol_hint is None and flt.protocol_id is None:
            resolved_protocols.append(None)
        else:
            rp = resolve_protocol(
                flt.protocol_hint or "", rt,
                pinned_id=flt.protocol_id,
                filter_index=idx,
            )
            if not isinstance(rp, ResolvedProtocol):
                return rp
            resolved_protocols.append(rp)

        if flt.assayed_by_as_typed or flt.assayed_by_id:
            ru = resolve_user(
                flt.assayed_by_as_typed or "",
                assay_creators_qs(),
                field=FIELD_ASSAYED_BY,
                filter_index=idx,
                pinned_id=flt.assayed_by_id,
            )
            if not isinstance(ru, ResolvedUser):
                return ru
            resolved_assayed_bys.append(ru)
        else:
            resolved_assayed_bys.append(None)

    # Apply each filter and intersect the surviving compound sets.
    base_compound_ids = _base_scope_compound_ids(rt, selector, resolved_registered_by)
    surviving: Optional[Set[str]] = None

    for flt, rp, ru in zip(
        selector.measurement_filters, resolved_protocols, resolved_assayed_bys,
    ):
        matching = _filter_matching_compounds(rt, rp, flt, ru)
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
        scope_sentence=_scope_sentence(
            rt, selector, resolved_protocols,
            resolved_registered_by=resolved_registered_by,
            resolved_assayed_bys=resolved_assayed_bys,
        ),
    )


# ---------------------------------------------------------------------------
# Scope construction
# ---------------------------------------------------------------------------


def _apply_date_range(qs, field: str, dr: Optional[DateRange]):
    """Apply a half-open [after, before) date range to ``qs`` against the
    date/datetime column named ``field``. Returns the queryset unchanged
    when dr is None or both bounds are null."""
    if dr is None:
        return qs
    if dr.after:
        qs = qs.filter(**{f"{field}__date__gte": dr.after})
    if dr.before:
        qs = qs.filter(**{f"{field}__date__lt": dr.before})
    return qs


def _base_scope_compound_ids(
    rt: ResolvedTargets,
    selector: CompoundSelector,
    resolved_registered_by: Optional[ResolvedUser],
) -> Set[str]:
    """Compound IDs that pass the target scope + registered-date filter
    + registered-by filter.

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

    qs = _apply_date_range(qs, "registered_at", selector.registered_date_range)

    if resolved_registered_by is not None:
        qs = qs.filter(registered_by=resolved_registered_by.user)

    return set(qs.values_list("pk", flat=True))


def _filter_matching_compounds(
    rt: ResolvedTargets,
    rp: Optional[ResolvedProtocol],
    flt: MeasurementFilter,
    ru: Optional[ResolvedUser],
) -> Set[str]:
    """Compound IDs whose AnalysisResult rows satisfy this filter.

    A row satisfies the filter iff evaluate_row returns RowMatched for it.
    The row itself is per-measurement; a compound is in the set if *any*
    of its rows pass (compound-level OR within a filter).

    Optional per-filter constraints:
    - ``assay_date_range`` — filters rows by Assay.created_at.
    - ``ru`` — resolved assayed-by user; filters on Assay.created_by.
    """
    qs = _scoped_analysis_results(rt)
    if rp is not None:
        qs = qs.filter(data_series__assay__protocol=rp.protocol)
    qs = _apply_date_range(qs, "data_series__assay__created_at", flt.assay_date_range)
    if ru is not None:
        qs = qs.filter(data_series__assay__created_by=ru.user)

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
    *,
    resolved_registered_by: Optional[ResolvedUser] = None,
    resolved_assayed_bys: Optional[List[Optional[ResolvedUser]]] = None,
) -> str:
    """Human echo-back (§10). Target scope, optional registered-date and
    registered-by phrases, and the AND-joined filter list (if any)."""
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

    parens: List[str] = []
    reg_date_phrase = _date_range_phrase(selector.registered_date_range, "registered")
    if reg_date_phrase:
        parens.append(reg_date_phrase)
    if resolved_registered_by is not None:
        from .resolver import _user_display  # avoid circular at module load
        parens.append(f"registered by {_user_display(resolved_registered_by.user)}")
    if parens:
        base += " (" + ", ".join(parens) + ")"

    if resolved_assayed_bys is None:
        resolved_assayed_bys = [None] * len(selector.measurement_filters)

    filter_clauses: List[str] = []
    for flt, rp, ru in zip(
        selector.measurement_filters, resolved_protocols, resolved_assayed_bys,
    ):
        clause = _filter_clause(flt, rp, ru)
        if clause:
            filter_clauses.append(clause)
    if filter_clauses:
        base += " where " + " AND ".join(filter_clauses)
    return base


def _filter_clause(
    flt: MeasurementFilter,
    rp: Optional[ResolvedProtocol],
    ru: Optional[ResolvedUser] = None,
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
    date_phrase = _date_range_phrase(flt.assay_date_range, "measured")
    if date_phrase:
        parts.append(date_phrase)
    if ru is not None:
        from .resolver import _user_display
        parts.append(f"assayed by {_user_display(ru.user)}")
    return " ".join(parts)


def _date_range_phrase(dr: Optional[DateRange], verb: str) -> str:
    """Human rendering of a DateRange — "<verb> after X", "<verb> before Y",
    or "<verb> between X and Y". Empty string when the range is absent or
    has neither bound set."""
    if dr is None:
        return ""
    if dr.after and dr.before:
        return f"{verb} between {dr.after} and {dr.before}"
    if dr.after:
        return f"{verb} after {dr.after}"
    if dr.before:
        return f"{verb} before {dr.before}"
    return ""
