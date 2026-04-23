"""End-to-end execution of a QuerySpec against the compounds database.

Implements §7 (executor) + §6.6 (phys-chem column expansion) + §10 (scope
echo-back sentence) of NLP_QUERY_PROPOSAL.md.

Composes the resolvers + row evaluator from earlier slices:
    resolve_targets  →  resolve_protocol  →  scoped AnalysisResult walk
                                              → evaluate_row per row

Returns one of:
- TablePayload            — rows + property columns + footer counts + scope sentence
- TargetClarify / TargetMiss / ScopeError    — passed through from resolve_targets
- ProtocolClarify / ProtocolMiss             — passed through from resolve_protocol
- SpecError               — spec-level validation failure (missing metric etc.)

Still no LLM. Slice 5 will build the QuerySpec from an LLM call; slice 6
will host the DRF view that drives this function.
"""

from __future__ import annotations

from typing import Dict, List, Optional

from django.db.models import QuerySet

from compounds.assays.models import AnalysisResult
from compounds.registry.models import Compound, Target

from .evaluator import evaluate_row
from .resolver import resolve_protocol, resolve_targets
from .spec import (
    COL_PRESET_LIPINSKI,
    COL_PRESET_PHYS_CHEM,
    LIPINSKI_COLUMNS,
    PHYS_CHEM_COLUMNS,
    SCOPE_ASSAY_ONLY,
    SCOPE_BOTH_SAME,
    SCOPE_CROSS,
    SCOPE_REG_ONLY,
    ExecutionResult,
    ProtocolClarify,
    ProtocolMiss,
    QuerySpec,
    ResolvedProtocol,
    ResolvedTargets,
    RowExcluded,
    RowFiltered,
    RowMatched,
    ScopeError,
    SpecError,
    TablePayload,
    TableRow,
    TargetClarify,
    TargetMiss,
)


_PHYS_CHEM_SET = frozenset(PHYS_CHEM_COLUMNS)


def _resolve_columns(spec: QuerySpec) -> List[str]:
    """Expand the spec's column choice into an ordered list of MolecularProperties
    field names (§6.6). Unknown explicit names are silently dropped; an empty
    outcome falls through to the phys-chem default so the payload isn't blank."""
    if spec.columns_explicit:
        kept = [c for c in spec.columns_explicit if c in _PHYS_CHEM_SET]
        if kept:
            return kept
    if spec.columns == COL_PRESET_LIPINSKI:
        return list(LIPINSKI_COLUMNS)
    return list(PHYS_CHEM_COLUMNS)


def _scoped_results(rp: ResolvedProtocol, rt: ResolvedTargets) -> QuerySet:
    """AnalysisResult queryset filtered to the resolved protocol + target scope.

    Pre-filters `status='valid'` — the evaluator would have classified invalid
    rows as FILTER_INVALID_STATUS, but those are silent anyway and pre-filtering
    avoids hauling them back from the DB.
    """
    qs = AnalysisResult.objects.filter(
        status="valid",
        data_series__assay__protocol=rp.protocol,
    )
    if rt.registration is not None:
        qs = qs.filter(data_series__compound__target=rt.registration)
    if rt.assay is not None:
        qs = qs.filter(data_series__assay__target=rt.assay)
    return qs.select_related(
        "data_series__compound__molecular_properties",
    )


def _pluck_properties(compound: Compound, columns: List[str]) -> Dict[str, Optional[float]]:
    props = getattr(compound, "molecular_properties", None)
    if props is None:
        return {c: None for c in columns}
    return {c: getattr(props, c, None) for c in columns}


def _threshold_clause(spec: QuerySpec) -> str:
    if spec.metric and spec.threshold is not None:
        t = spec.threshold
        unit = f" {t.unit}" if t.unit else ""
        return f", with {spec.metric} {t.op} {t.value}{unit}"
    if spec.metric:
        return f", {spec.metric} values"
    return ""


def _target_name(target: Optional[Target]) -> str:
    return target.name if target is not None else ""


def _scope_sentence(rt: ResolvedTargets, rp: ResolvedProtocol, spec: QuerySpec) -> str:
    protocol = rp.protocol.name
    reg = _target_name(rt.registration)
    assay = _target_name(rt.assay)
    if rt.scope_kind == SCOPE_BOTH_SAME:
        base = f"Showing compounds registered to {reg} and tested in {protocol} against {reg}"
    elif rt.scope_kind == SCOPE_CROSS:
        base = f"Showing compounds registered to {reg}, tested in {protocol} against {assay}"
    elif rt.scope_kind == SCOPE_REG_ONLY:
        base = f"Showing compounds registered to {reg}, tested in {protocol}"
    elif rt.scope_kind == SCOPE_ASSAY_ONLY:
        base = f"Showing compounds tested in {protocol} against {assay}"
    else:
        base = f"Showing compounds tested in {protocol}"
    return base + _threshold_clause(spec)


def execute(spec: QuerySpec) -> ExecutionResult:
    """Run a fully-formed QuerySpec and return the table payload or a
    clarify/miss/error response to surface back to the user."""
    if not spec.metric:
        return SpecError(field="metric", message="metric is required")
    if not spec.protocol_hint and not spec.protocol_id:
        return SpecError(field="protocol_hint", message="protocol_hint is required")

    rt = resolve_targets(spec)
    if not isinstance(rt, ResolvedTargets):
        return rt

    rp = resolve_protocol(spec.protocol_hint or "", rt, pinned_id=spec.protocol_id)
    if not isinstance(rp, ResolvedProtocol):
        return rp

    columns = _resolve_columns(spec)
    qs = _scoped_results(rp, rt)

    rows: List[TableRow] = []
    footer_excluded: Dict[str, int] = {}
    filtered_silent: Dict[str, int] = {}

    for ar in qs.iterator():
        outcome = evaluate_row(ar, spec.metric, spec.threshold)
        if isinstance(outcome, RowMatched):
            ds = getattr(ar, "data_series", None)
            compound = ds.compound if ds is not None else None
            if compound is None:
                continue
            rows.append(
                TableRow(
                    compound_id=str(compound.id),
                    formatted_id=compound.formatted_id,
                    smiles=compound.smiles,
                    value=outcome.value,
                    value_unit=outcome.row_unit,
                    value_in_query_unit=outcome.value_in_query_unit,
                    properties=_pluck_properties(compound, columns),
                )
            )
        elif isinstance(outcome, RowExcluded):
            footer_excluded[outcome.reason] = footer_excluded.get(outcome.reason, 0) + 1
        elif isinstance(outcome, RowFiltered):
            filtered_silent[outcome.reason] = filtered_silent.get(outcome.reason, 0) + 1

    return TablePayload(
        rows=rows,
        property_columns=columns,
        footer_excluded=footer_excluded,
        filtered_silent=filtered_silent,
        scope_sentence=_scope_sentence(rt, rp, spec),
    )
