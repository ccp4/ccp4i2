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

from typing import Any, Dict, List, Optional, Set, Tuple, Union

from django.db.models import QuerySet

from compounds.assays.models import AnalysisResult
from compounds.formatting import format_compound_id
from compounds.registry.models import Compound, Target

from .evaluator import evaluate_row
from .resolver import (
    assay_creators_qs,
    compound_registrars_qs,
    resolve_compound_ref,
    resolve_protocol,
    resolve_registrant,
    resolve_scaffold,
    resolve_target,
    resolve_targets,
    resolve_user,
)
from .spec import (
    DEFAULT_RANK_TOP_N,
    FIELD_ASSAY_TARGET,
    FIELD_ASSAYED_BY,
    FIELD_METRIC,
    FIELD_RANK_BY,
    FIELD_REGISTERED_BY,
    RANK_BY_SCORECARD,
    SCOPE_ASSAY_ONLY,
    SCOPE_BOTH_SAME,
    SCOPE_CROSS,
    SCOPE_REG_ONLY,
    AssaySelection,
    AssaySelector,
    CompoundSelection,
    CompoundSelector,
    DateRange,
    ExecutionResult,
    MeasurementFilter,
    MetricMiss,
    ProtocolClarify,
    ProtocolMiss,
    ResolvedCompound,
    ResolvedProtocol,
    ResolvedScaffold,
    ResolvedSupplier,
    ResolvedTarget,
    ResolvedTargets,
    ResolvedUnion,
    ResolvedUser,
    RowMatched,
    ScopeError,
    SpecError,
    TargetClarify,
    TargetMiss,
    Threshold,
)

from compounds.assays.models import Assay


def _selector_has_non_pin_narrowing(selector: CompoundSelector) -> bool:
    """True iff the selector has any narrowing predicate other than
    ``compound_refs_as_typed``. Used to decide whether an unscoped
    selector should default its base set to "every compound in the
    registry" (so a non-pin predicate can narrow it) or to the empty
    set (so the pin set is the entire selection).

    Slice 21: ``similar_to_as_typed`` counts as a narrowing predicate
    — *"compounds similar to NCL-26007"* without a target should run
    against every compound in the registry and narrow via Tanimoto.
    """
    return bool(
        selector.registration_target_as_typed
        or selector.registration_target_id
        or selector.assay_target_as_typed
        or selector.assay_target_id
        or selector.scaffold_hints
        or selector.scaffold_ids
        or selector.registered_by_as_typed
        or selector.registered_by_id
        or selector.registered_date_range is not None
        or selector.measurement_filters
        or selector.similar_to_as_typed
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

    # A fully-empty selector would hand back every compound in the
    # registry with no narrowing — almost certainly not what the user
    # asked for. Guard with a friendly error that actually tells them
    # what they need to add. A target-less selector is fine AS LONG AS
    # some other predicate narrows (scaffold / user / date / filter).
    if not _selector_has_non_pin_narrowing(selector) and not selector.compound_refs_as_typed:
        return ScopeError(
            message=(
                "Your query doesn't name anything to narrow by — no target "
                "(like CDK4 / ARd / mEGFR), no substructure, no person, no "
                "date, and no measurement threshold. Try adding at least "
                "one. Examples: \"CDK4 compounds with HTRF IC50 < 100 nM\", "
                "\"compounds registered by Suzannah\", "
                "\"all pyrimidines\", or \"compounds with HTRF IC50 < 10 nM\"."
            )
        )

    rt = resolve_targets(selector)
    if not isinstance(rt, ResolvedTargets):
        return rt

    # Selector-level registrant resolution — poly-source (User OR
    # Supplier, slice 16). The chemist's "made by" / "registered by" /
    # "from" phrasing matches against both pools; the resolver de-dupes
    # user-linked Suppliers against their Users.
    resolved_registered_by = None
    if (
        selector.registered_by_as_typed
        or selector.registered_by_id
        or selector.registered_by_supplier_id
    ):
        rr = resolve_registrant(
            selector.registered_by_as_typed or "",
            field=FIELD_REGISTERED_BY,
            pinned_user_id=selector.registered_by_id,
            pinned_supplier_id=selector.registered_by_supplier_id,
        )
        if not isinstance(rr, (ResolvedUser, ResolvedSupplier, ResolvedUnion)):
            return rr
        resolved_registered_by = rr

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

    # Scaffold filter — if any scaffold_hints are named, resolve each in
    # order and RDKit-match against the base compound scope. All
    # scaffolds ANDed. First unresolved scaffold short-circuits the
    # whole selector.
    resolved_scaffolds: List[ResolvedScaffold] = []
    hints = list(selector.scaffold_hints or [])
    pinned = list(selector.scaffold_ids or [])
    # Align pinned list to hint positions — pad with empty strings.
    while len(pinned) < len(hints):
        pinned.append("")
    # Slice 17: project-scoped extensions take precedence over shared
    # extensions and seed. The registration target (when set) drives
    # "project context"; assay-only scope falls through to shared+seed
    # because per-project nomenclature is naturally tied to the
    # registration programme, not the assay readout target.
    scaffold_target = rt.registration if rt.registration is not None else None
    for idx, hint in enumerate(hints):
        rs = resolve_scaffold(
            hint, target=scaffold_target, scaffold_index=idx,
            pinned_id=pinned[idx] if pinned[idx] else None,
        )
        if not isinstance(rs, ResolvedScaffold):
            return rs
        resolved_scaffolds.append(rs)
    # Apply scaffold narrowing BEFORE measurement filters — compound-
    # level predicate, cheap to intersect first.
    if resolved_scaffolds:
        for rs in resolved_scaffolds:
            base_compound_ids = (
                base_compound_ids & _match_compounds_by_scaffold(base_compound_ids, rs)
            )

    # Slice 21: Tanimoto similarity narrowing. Anchors are typed
    # compound IDs (same resolver as compound_refs_as_typed); the
    # narrowing is "≥ threshold similar to ANY anchor" (UNION across
    # anchors). Applied after scaffold but before measurement filters
    # because compound-level predicates compose cheaply at this layer.
    # First unresolved anchor short-circuits to CompoundMiss.
    similar_refs = list(selector.similar_to_as_typed or [])
    resolved_similar_to: List[ResolvedCompound] = []
    if similar_refs:
        for idx, ref in enumerate(similar_refs):
            rc = resolve_compound_ref(ref, ref_index=idx)
            if not isinstance(rc, ResolvedCompound):
                return rc
            resolved_similar_to.append(rc)
        from .similarity import neighbours_within_threshold, DEFAULT_SIMILAR_THRESHOLD
        threshold = (
            selector.similar_threshold
            if selector.similar_threshold is not None
            else DEFAULT_SIMILAR_THRESHOLD
        )
        base_compound_ids = base_compound_ids & neighbours_within_threshold(
            [rc.compound for rc in resolved_similar_to], base_compound_ids, threshold,
        )

    # Pinned compound refs — UNION semantics, not narrowing. The user
    # named specific compounds by ID and they should appear in the
    # final selection regardless of whether they pass the other
    # predicates. A first unresolved ref short-circuits to a
    # CompoundMiss so the user sees which ID went wrong.
    resolved_pins: List[ResolvedCompound] = []
    for idx, ref in enumerate(selector.compound_refs_as_typed or []):
        rc = resolve_compound_ref(ref, ref_index=idx)
        if not isinstance(rc, ResolvedCompound):
            return rc
        resolved_pins.append(rc)

    surviving: Optional[Set[str]] = None

    for idx, (flt, rp, ru) in enumerate(zip(
        selector.measurement_filters, resolved_protocols, resolved_assayed_bys,
    )):
        # Slice 15: surface metric mismatches as a self-diagnostic Miss
        # rather than letting the strict-match filter silently drop every
        # row. If the user typed "IC50" but the rows in scope record
        # "pIC50" / "EC50" / etc., the response answers itself with the
        # available KPIs.
        if flt.metric:
            metric_miss = _check_metric_in_scope(rt, rp, flt, ru, idx)
            if metric_miss is not None:
                return metric_miss
        matching = _filter_matching_compounds(rt, rp, flt, ru)
        surviving = matching if surviving is None else surviving & matching

    if surviving is None:
        # No filters — every compound in the base scope is in the selection.
        compound_ids = base_compound_ids
    else:
        compound_ids = surviving & base_compound_ids
    if resolved_pins:
        # UNION the pinned set on AFTER all narrowing — the lead compound
        # appears in the final selection even when it doesn't pass the
        # measurement filter or sit in the target scope.
        compound_ids |= {rc.compound.pk for rc in resolved_pins}

    # Slice 19: ranking. When the chemist says "best X compounds" the
    # LLM sets rank_by="scorecard"; we rank the surviving set by the
    # target's scorecard score and take top N. Ranking happens AFTER
    # all narrowing so the chemist sees the best of what matched.
    rank_order: Optional[List[Any]] = None  # ordered list of compound pks (top first)
    if selector.rank_by:
        rank_result = _apply_ranking(rt, selector, compound_ids)
        if isinstance(rank_result, SpecError):
            return rank_result
        compound_ids, rank_order = rank_result

    # Hydrate compound IDs into formatted identifiers. When ranked, use
    # the rank order; otherwise fall back to reg_number ASC for stable
    # URLs.
    if rank_order is not None:
        formatted = list(
            Compound.objects
            .filter(pk__in=compound_ids)
            .in_bulk(rank_order, field_name="pk")  # preserves dict order
            .values()
        )
        formatted_ids = [format_compound_id(c.reg_number) for c in formatted]
        # in_bulk returns a dict keyed by pk; iterate rank_order to keep
        # the ranking order in the URL.
        compound_lookup = {
            c.pk: format_compound_id(c.reg_number) for c in formatted
        }
        formatted_ids = [compound_lookup[pk] for pk in rank_order if pk in compound_lookup]
    else:
        formatted = list(
            Compound.objects
            .filter(pk__in=compound_ids)
            .order_by("reg_number")
            .values_list("reg_number", flat=True)
        )
        formatted_ids = [format_compound_id(rn) for rn in formatted]

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
            resolved_scaffolds=resolved_scaffolds,
            resolved_pins=resolved_pins,
            resolved_similar_to=resolved_similar_to,
            similar_threshold=selector.similar_threshold,
        ),
        view_format=selector.view_format,
        # categorisation_selection_ids is resolved in view.py against
        # the request user's saved selections (per-user scope).
    )


# ---------------------------------------------------------------------------
# Scope construction
# ---------------------------------------------------------------------------


def _match_compounds_by_scaffold(
    compound_ids: Set[str],
    resolved_scaffold: ResolvedScaffold,
) -> Set[str]:
    """Return the subset of ``compound_ids`` whose SMILES contains the
    scaffold's SMARTS substructure (RDKit ``HasSubstructMatch``).

    Imported lazily so the substructures module and this function only
    pull in RDKit when a prompt actually names a scaffold — the rest of
    the executor stays RDKit-free.
    """
    from rdkit import Chem
    scaffold = resolved_scaffold.scaffold
    patt = Chem.MolFromSmarts(scaffold.smarts)
    if patt is None:
        # Malformed SMARTS in the catalog — should be caught at commit-
        # time, but safe-fail here by excluding everything.
        return set()
    matching: Set[str] = set()
    rows = (
        Compound.objects
        .filter(pk__in=compound_ids)
        .exclude(smiles__isnull=True).exclude(smiles="")
        .values_list("pk", "smiles")
    )
    for pk, smiles in rows:
        mol = Chem.MolFromSmiles(smiles)
        if mol and mol.HasSubstructMatch(patt):
            matching.add(pk)
    return matching


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
    resolved_registered_by: Optional[object],   # ResolvedUser | ResolvedSupplier | None
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
    elif rt.assay is not None:
        # assay-only scope
        qs = Compound.objects.filter(assay_results__assay__target=rt.assay).distinct()
    elif _selector_has_non_pin_narrowing(selector):
        # Unscoped, but some non-pin predicate (scaffold / user / date /
        # filter) will narrow the result — start from every compound.
        qs = Compound.objects.all()
    else:
        # Pin-only selector: no target and no other narrowing predicate.
        # Skip the "all compounds" base so the final selection is just
        # the pinned set, not the entire registry.
        return set()

    if rt.scope_kind == SCOPE_CROSS and rt.assay is not None:
        # Registration target filter above; also require the compound to
        # have at least one DataSeries against the assay target so we don't
        # surface "my ARd compounds" with no actual AKT data.
        qs = qs.filter(assay_results__assay__target=rt.assay).distinct()

    qs = _apply_date_range(qs, "registered_at", selector.registered_date_range)

    if resolved_registered_by is not None:
        from django.db.models import Q
        if isinstance(resolved_registered_by, ResolvedUser):
            user = resolved_registered_by.user
            # Match compounds where EITHER the registered_by FK is this
            # user OR the supplier is this user's user-linked Supplier.
            # Chemists conflate "registered by" with "supplied by self"
            # routinely; exposing both via one filter mirrors that.
            qs = qs.filter(Q(registered_by=user) | Q(supplier__user=user))
        elif isinstance(resolved_registered_by, ResolvedUnion):
            # Slice 18: same-name (User, unlinked-Supplier) pair —
            # union of all three paths. Hannah-Stewart-the-User's three
            # registered compounds + Hannah-Stewart-the-Supplier's two
            # supplied compounds, all surfaced together.
            user = resolved_registered_by.user
            sup = resolved_registered_by.supplier
            qs = qs.filter(
                Q(registered_by=user) | Q(supplier__user=user) | Q(supplier=sup)
            )
        else:
            # ResolvedSupplier — filter on the supplier FK directly.
            qs = qs.filter(supplier=resolved_registered_by.supplier)

    return set(qs.values_list("pk", flat=True))


def _filter_scope_qs(
    rt: ResolvedTargets,
    rp: Optional[ResolvedProtocol],
    flt: MeasurementFilter,
    ru: Optional[ResolvedUser],
) -> QuerySet:
    """The AnalysisResult queryset narrowed to the filter's scope —
    target + protocol + assay date + assayed-by. Shared by metric
    pre-check and matching-compound iteration."""
    qs = _scoped_analysis_results(rt)
    if rp is not None:
        qs = qs.filter(data_series__assay__protocol=rp.protocol)
    qs = _apply_date_range(qs, "data_series__assay__created_at", flt.assay_date_range)
    if ru is not None:
        qs = qs.filter(data_series__assay__created_by=ru.user)
    return qs


def _check_metric_in_scope(
    rt: ResolvedTargets,
    rp: Optional[ResolvedProtocol],
    flt: MeasurementFilter,
    ru: Optional[ResolvedUser],
    filter_index: int,
) -> Optional[MetricMiss]:
    """If no row in the filter's scope has a KPI matching ``flt.metric``
    (lenient: case + non-alphanumeric punctuation insensitive), return
    a MetricMiss listing every distinct KPI that IS in scope. Otherwise
    return None and let the matching loop run.

    Cheap: one pass over the scope queryset, projecting just the JSON
    KPI key. The filter scope is usually small (target-narrowed, often
    protocol-narrowed too)."""
    from compounds.utils import normalize_ref as normalize  # avoid module-level cycle
    metric = flt.metric or ""
    if not metric:
        return None
    target_norm = normalize(metric)
    seen: Set[str] = set()
    qs = _filter_scope_qs(rt, rp, flt, ru)
    for results in qs.values_list("results", flat=True):
        kpi = (results or {}).get("KPI") if isinstance(results, dict) else None
        if not kpi:
            continue
        kpi_str = str(kpi)
        if normalize(kpi_str) == target_norm:
            return None  # found it — the regular evaluator will handle it
        seen.add(kpi_str)
    return MetricMiss(
        query=metric,
        available_metrics=sorted(seen),
        filter_index=filter_index,
        field=FIELD_METRIC,
    )


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
    qs = _filter_scope_qs(rt, rp, flt, ru)

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
# Ranking (slice 19) — "best X compounds" / "top N X compounds"
# ---------------------------------------------------------------------------


def _apply_ranking(
    rt: ResolvedTargets,
    selector: CompoundSelector,
    compound_ids: Set,
) -> Union[Tuple[Set, List], SpecError]:
    """Rank the surviving compound set by the configured ``rank_by``
    metric and trim to ``rank_top_n``. Returns ``(top_pks, ordered_pks)``
    on success — ``top_pks`` is the trimmed set for downstream consumers
    that don't care about order; ``ordered_pks`` is the rank-ordered
    list for the redirect URL.

    Returns a ``SpecError`` when:
    - ``rank_by`` is unrecognised
    - ``rank_by == "scorecard"`` but the resolved target has no
      ``scorecard_config`` (the chemist meant "best", but the project
      hasn't defined what "best" means yet — surface a clear error
      rather than silently fall back to some other metric).
    """
    rank_by = (selector.rank_by or "").strip()
    if rank_by != RANK_BY_SCORECARD:
        return SpecError(
            field=FIELD_RANK_BY,
            message=(
                f"Unknown rank_by value {rank_by!r}; only "
                f"{RANK_BY_SCORECARD!r} is supported."
            ),
        )

    # Pick which target's scorecard to rank by. The registration target
    # is the natural anchor (project context); fall back to the assay
    # target if only that's set.
    rank_target = rt.registration or rt.assay
    if rank_target is None:
        return SpecError(
            field=FIELD_RANK_BY,
            message=(
                'Cannot rank by scorecard without a target — name a '
                'project (e.g. "best EGFR compounds").'
            ),
        )

    if not rank_target.scorecard_config or not (
        rank_target.scorecard_config or {}
    ).get("axes"):
        return SpecError(
            field=FIELD_RANK_BY,
            message=(
                f'Target {rank_target.name!r} has no scorecard '
                'configured. Configure a scorecard for this project, '
                'or drop the "best" phrasing and use an explicit '
                'threshold (e.g. "HTRF IC50 < 10 nM").'
            ),
        )

    from .scorecard_rank import rank_compounds_by_scorecard
    try:
        scores = rank_compounds_by_scorecard(rank_target, compound_ids)
    except ValueError as e:
        # Defensive — should be caught by the no-scorecard check above.
        return SpecError(field=FIELD_RANK_BY, message=str(e))

    # Sort: highest score first, then lower reg_number for stable ties.
    # The aggregate_compact_from_compounds output keys by stringified
    # UUIDs while compound_ids may contain UUID objects — normalise.
    scored_pks = [
        (pk, score) for pk, score in scores.items()
    ]
    scored_pks.sort(key=lambda row: (-row[1], str(row[0])))

    top_n = selector.rank_top_n if selector.rank_top_n else DEFAULT_RANK_TOP_N
    top = scored_pks[:top_n]

    # Map back to the original pk type — compound_ids may contain UUID
    # objects but the aggregator returns strings. Build a set of stringified
    # pks for membership check, and use the original UUID where possible.
    str_to_orig: Dict[str, Any] = {str(pk): pk for pk in compound_ids}
    ordered_pks: List = []
    for pk_str, _score in top:
        orig = str_to_orig.get(str(pk_str))
        if orig is not None:
            ordered_pks.append(orig)

    return set(ordered_pks), ordered_pks


# ---------------------------------------------------------------------------
# Presentation helpers
# ---------------------------------------------------------------------------


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
    resolved_registered_by: Optional[object] = None,   # ResolvedUser | ResolvedSupplier | None
    resolved_assayed_bys: Optional[List[Optional[ResolvedUser]]] = None,
    resolved_scaffolds: Optional[List[ResolvedScaffold]] = None,
    resolved_pins: Optional[List[ResolvedCompound]] = None,
    resolved_similar_to: Optional[List[ResolvedCompound]] = None,
    similar_threshold: Optional[float] = None,
) -> str:
    """Human echo-back. Target scope, optional registered-date and
    registered-by phrases, scaffold narrowing, the AND-joined filter
    list, similarity narrowing, and any pinned-compound UNION clause."""
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
        if isinstance(resolved_registered_by, ResolvedUser):
            from .resolver import _user_display  # avoid circular at module load
            parens.append(
                f"registered by {_user_display(resolved_registered_by.user)}"
            )
        elif isinstance(resolved_registered_by, ResolvedUnion):
            # Slice 18: same-name union — both roles ORed in the filter,
            # so the prose reads "registered or supplied by".
            from .resolver import _user_display
            parens.append(
                f"registered or supplied by "
                f"{_user_display(resolved_registered_by.user)}"
            )
        else:
            # ResolvedSupplier — render the supplier name verbatim. "from"
            # reads better than "registered by" for a vendor.
            parens.append(f"from {resolved_registered_by.supplier.name}")
    if parens:
        base += " (" + ", ".join(parens) + ")"

    if resolved_scaffolds:
        names = [rs.scaffold.name for rs in resolved_scaffolds]
        if len(names) == 1:
            base += f" containing {names[0]}"
        else:
            base += " containing " + " AND ".join(names)

    if resolved_similar_to:
        ids = [
            format_compound_id(rc.compound.reg_number) for rc in resolved_similar_to
        ]
        from .similarity import DEFAULT_SIMILAR_THRESHOLD
        t = similar_threshold if similar_threshold is not None else DEFAULT_SIMILAR_THRESHOLD
        anchors = ids[0] if len(ids) == 1 else " or ".join(ids)
        base += f" similar to {anchors} (Tanimoto ≥ {t:g})"

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

    if resolved_pins:
        pinned_ids = [format_compound_id(rc.compound.reg_number) for rc in resolved_pins]
        base += f", plus {', '.join(pinned_ids)} (pinned)"

    # Slice 19: ranking phrase prefixed onto the base sentence so the
    # chemist sees "top 20" before reading the scope predicates.
    if selector.rank_by == RANK_BY_SCORECARD:
        n = selector.rank_top_n or DEFAULT_RANK_TOP_N
        base = f"Top {n} by scorecard, " + base[0].lower() + base[1:]
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


# ---------------------------------------------------------------------------
# Assay-selection path — introduced in the §19.7 "assay list" follow-up.
# Parallel to the compound-selection path above but simpler: a single
# queryset against Assay, filtered by target / protocol / created_by /
# date, returning Assay IDs for the /assays?ids=... redirect.
# ---------------------------------------------------------------------------


def execute_assay_query(selector: AssaySelector) -> ExecutionResult:
    """Run an AssaySelector and return an AssaySelection, or a
    clarify/miss/error to surface back to the user.

    Unlike the compound path, an AssaySelector has no list of filters
    and no multi-predicate intersection — it's a straight queryset
    against Assay with optional target/protocol/creator/date
    predicates. The result is a flat list of assay IDs for the
    aggregation-page's sibling, `/assays?ids=...`.
    """
    # Target resolution. The assay query has a single target field
    # (not the two-target dance of CompoundSelector).
    resolved_target: Optional[ResolvedTarget] = None
    if selector.target_as_typed or selector.target_id:
        rt_single = resolve_target(
            selector.target_as_typed or "",
            pinned_id=selector.target_id,
        )
        if isinstance(rt_single, TargetClarify):
            # Tag the field so the UI routes the picker correctly.
            from dataclasses import replace
            return replace(rt_single, field=FIELD_ASSAY_TARGET)
        if isinstance(rt_single, TargetMiss):
            from dataclasses import replace
            return replace(rt_single, field=FIELD_ASSAY_TARGET)
        resolved_target = rt_single

    # Protocol resolution — scope-bound to the resolved target (if any).
    # Reuse resolve_protocol with a ResolvedTargets shaped as assay-only
    # so _scope_assays narrows protocols to ones that have actually been
    # run against this target.
    resolved_protocol: Optional[ResolvedProtocol] = None
    if selector.protocol_hint or selector.protocol_id:
        rt_proxy = ResolvedTargets(
            registration=None,
            assay=resolved_target.target if resolved_target is not None else None,
            scope_kind=SCOPE_ASSAY_ONLY if resolved_target is not None else SCOPE_REG_ONLY,
        )
        rp = resolve_protocol(
            selector.protocol_hint or "",
            rt_proxy,
            pinned_id=selector.protocol_id,
        )
        if not isinstance(rp, ResolvedProtocol):
            return rp
        resolved_protocol = rp

    # User resolution — scope-bound to users who have actually created
    # any assay (not the whole Django auth table).
    resolved_creator: Optional[ResolvedUser] = None
    if selector.created_by_as_typed or selector.created_by_id:
        ru = resolve_user(
            selector.created_by_as_typed or "",
            assay_creators_qs(),
            field=FIELD_ASSAYED_BY,
            pinned_id=selector.created_by_id,
        )
        if not isinstance(ru, ResolvedUser):
            return ru
        resolved_creator = ru

    # Scaffold filter — resolve each hint in order, RDKit-match each
    # against the set of compounds that have been assayed in this scope,
    # then narrow the assay queryset to assays whose DataSeries hit a
    # matching compound.
    resolved_scaffolds: List[ResolvedScaffold] = []
    hints = list(selector.scaffold_hints or [])
    pinned = list(selector.scaffold_ids or [])
    while len(pinned) < len(hints):
        pinned.append("")
    # Slice 17: assay-side scope is the assay target, since assay
    # queries don't have a registration-target concept distinct from
    # the assay one. Project-scoped extensions resolve against it.
    scaffold_target = resolved_target.target if resolved_target is not None else None
    for idx, hint in enumerate(hints):
        rs = resolve_scaffold(
            hint, target=scaffold_target, scaffold_index=idx,
            pinned_id=pinned[idx] if pinned[idx] else None,
        )
        if not isinstance(rs, ResolvedScaffold):
            return rs
        resolved_scaffolds.append(rs)

    # Build the assay queryset.
    qs = Assay.objects.all()
    if resolved_target is not None:
        qs = qs.filter(target=resolved_target.target)
    if resolved_protocol is not None:
        qs = qs.filter(protocol=resolved_protocol.protocol)
    if resolved_creator is not None:
        qs = qs.filter(created_by=resolved_creator.user)
    qs = _apply_date_range(qs, "created_at", selector.date_range)

    if resolved_scaffolds:
        # Compounds in scope of the current assay queryset.
        candidate_compound_ids = set(
            Compound.objects
            .filter(assay_results__assay__in=qs)
            .values_list("pk", flat=True)
            .distinct()
        )
        for rs in resolved_scaffolds:
            candidate_compound_ids = (
                candidate_compound_ids & _match_compounds_by_scaffold(candidate_compound_ids, rs)
            )
        qs = qs.filter(data_series__compound__in=candidate_compound_ids).distinct()

    assay_ids = [
        str(pk) for pk in qs.order_by("-created_at").values_list("pk", flat=True)
    ]

    target_names = [resolved_target.target.name] if resolved_target is not None else []
    protocol_names = (
        [resolved_protocol.protocol.name] if resolved_protocol is not None else []
    )

    return AssaySelection(
        assay_ids=assay_ids,
        n_matched=len(assay_ids),
        target_names=target_names,
        protocol_names=protocol_names,
        scope_sentence=_assay_scope_sentence(
            selector, resolved_target, resolved_protocol, resolved_creator,
            resolved_scaffolds=resolved_scaffolds,
        ),
    )


def _assay_scope_sentence(
    selector: AssaySelector,
    target: Optional[ResolvedTarget],
    protocol: Optional[ResolvedProtocol],
    creator: Optional[ResolvedUser],
    *,
    resolved_scaffolds: Optional[List[ResolvedScaffold]] = None,
) -> str:
    parts: List[str] = ["Showing assays"]
    if protocol is not None:
        parts.append(f"using {protocol.protocol.name}")
    if target is not None:
        parts.append(f"against {target.target.name}")
    if resolved_scaffolds:
        names = [rs.scaffold.name for rs in resolved_scaffolds]
        parts.append("on compounds containing " + " AND ".join(names))
    date_phrase = _date_range_phrase(selector.date_range, "carried out")
    if date_phrase:
        parts.append(date_phrase)
    if creator is not None:
        from .resolver import _user_display
        parts.append(f"conducted by {_user_display(creator.user)}")
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
