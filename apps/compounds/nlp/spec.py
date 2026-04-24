"""CompoundSelector + resolver + executor result types.

**Pivot (2026-04-23)**: this subsystem now produces a *compound selection*
rather than a rendered table. The LLM describes *which compounds* to view;
the executor intersects filter sets to produce a list of compound IDs; the
view returns a redirect URL to the existing `/assays/aggregate` page where
the user picks their display format (cards, compact, pivot, …). The
aggregation page already owns everything the previous `TablePayload` tried
to reproduce (phys-chem columns, unit display, geomean aggregations,
spider plots via project cards).

See NLP_QUERY_PROPOSAL.md §5–§7 for the data model and §10 for the UI.

Shape mirrors to:
- Backend: `apps/compounds/nlp/spec.py`     (this file)
- Frontend: `apps/compounds/frontend/lib/compounds/nlp-api.ts`
Keep the two in sync.
"""

from __future__ import annotations

import datetime
from dataclasses import dataclass, field
from typing import List, Optional, Union

from compounds.assays.models import Protocol
from compounds.registry.models import Target


# ---------------------------------------------------------------------------
# String constants — scope_kind, field-name tags, matched_via labels.
# Plain strings keep JSON round-trip trivial and avoid enum-in-dict dance.
# ---------------------------------------------------------------------------

SCOPE_BOTH_SAME = "both_same"
SCOPE_REG_ONLY = "reg_only"
SCOPE_ASSAY_ONLY = "assay_only"
SCOPE_CROSS = "cross"

FIELD_REGISTRATION_TARGET = "registration_target_as_typed"
FIELD_ASSAY_TARGET = "assay_target_as_typed"
FIELD_PROTOCOL_HINT = "protocol_hint"

# Filter / exclusion reasons from the row evaluator. Preserved from the pre-
# pivot shape because the evaluator still classifies rows during selection;
# only FILTER_* and EXCLUDE_* counts are no longer surfaced to the client
# (the selection payload just carries the surviving compound IDs).
FILTER_INVALID_STATUS = "invalid_status"
FILTER_KPI_MISMATCH = "kpi_mismatch"
FILTER_VALUE_NOT_NUMERIC = "value_not_numeric"
FILTER_THRESHOLD_NOT_MET = "threshold_not_met"

EXCLUDE_UNIT_UNKNOWN = "unit_unknown"
EXCLUDE_UNIT_TYPE_MISMATCH = "unit_type_mismatch"
EXCLUDE_UNIT_INCOMPATIBLE = "unit_incompatible"
EXCLUDE_QUERY_MISSING_UNIT = "query_missing_unit"


# ---------------------------------------------------------------------------
# Input — Threshold + MeasurementFilter + CompoundSelector
# ---------------------------------------------------------------------------


@dataclass
class Threshold:
    op: str                          # one of <, <=, >, >=, =, ==, !=
    value: float
    unit: Optional[str] = None       # None → unitless / inherit row unit


@dataclass
class DateRange:
    """ISO-date half-open filter range: [after, before).

    Either end may be null (unbounded on that side). Both null is
    treated as no filter. The half-open convention makes calendar-unit
    phrasings exact — "in 2025" is {after: 2025-01-01, before: 2026-01-01}
    with no fence-post ambiguity.

    Dates are strings in "YYYY-MM-DD" form. The LLM interprets relative
    phrasings ("last 30 days", "since March") against the `[Today: ...]`
    line prepended to the user message.
    """

    after: Optional[str] = None      # inclusive lower bound
    before: Optional[str] = None     # exclusive upper bound


@dataclass
class MeasurementFilter:
    """A single condition compounds must satisfy to be included in the
    selection. Multiple filters on a `CompoundSelector` are ANDed — the
    resulting compound set is the intersection of compounds passing each.

    Null-field semantics:
    - `protocol_hint = None` → filter doesn't scope to a specific protocol
      (any protocol under the selector's target scope is eligible).
    - `metric = None` → filter doesn't scope to a specific metric (any
      numeric KPI counts).
    - `threshold = None` → filter requires only that *some* valid
      measurement exists matching the protocol/metric constraints; any
      numeric value passes. Useful for "compounds tested in X" prompts.
    """

    protocol_hint: Optional[str] = None
    metric: Optional[str] = None
    threshold: Optional[Threshold] = None
    # Date range on the Assay (Assay.created_at). Limits the compound set
    # to those whose measurements under this filter's protocol fall in
    # the window — "HTRF IC50 < 100 nM measured in Q1 2026".
    assay_date_range: Optional[DateRange] = None
    # Pinned protocol id on clarify continuation — LLM never emits this.
    protocol_id: Optional[str] = None


@dataclass
class CompoundSelector:
    """Root of the v2 spec — describes *which compounds* to surface.

    Target fields follow the same §6.5 semantics as the pre-pivot spec:
    - both populated: registered to X AND tested against Y (default X==Y)
    - only registration: all compounds in that programme, whatever their assays
    - only assay: any compound with measurements against that target

    Pinning fields are populated by the view on clarify continuation (§9);
    the LLM never emits them.
    """

    registration_target_as_typed: Optional[str] = None
    assay_target_as_typed: Optional[str] = None
    measurement_filters: List[MeasurementFilter] = field(default_factory=list)
    # Date range on Compound.registered_at — filters the base compound
    # scope BEFORE any measurement filters intersect.
    registered_date_range: Optional[DateRange] = None

    # Pinnings from clarify continuation — the view injects these, LLM doesn't.
    registration_target_id: Optional[str] = None
    assay_target_id: Optional[str] = None


# ---------------------------------------------------------------------------
# Target resolution (single field)
# ---------------------------------------------------------------------------


@dataclass
class TargetCandidate:
    """Thin JSON-serialisable view of a Target for clarify / miss payloads."""

    id: str
    name: str
    gene_symbols: List[str] = field(default_factory=list)


@dataclass
class ResolvedTarget:
    target: Target
    matched_via: str
    query: str


@dataclass
class TargetClarify:
    query: str
    candidates: List[TargetCandidate]
    field: str = ""  # FIELD_REGISTRATION_TARGET / FIELD_ASSAY_TARGET — orchestrator sets


@dataclass
class TargetMiss:
    query: str
    suggestions: List[TargetCandidate]
    field: str = ""


TargetResolution = Union[ResolvedTarget, TargetClarify, TargetMiss]


# ---------------------------------------------------------------------------
# Two-target orchestrator results
# ---------------------------------------------------------------------------


@dataclass
class ResolvedTargets:
    registration: Optional[Target]
    assay: Optional[Target]
    scope_kind: str  # one of SCOPE_* above


@dataclass
class ScopeError:
    """At-least-one-target rule violated — selector had no target and no
    registration_* / assay_* pinning id."""

    message: str


TargetsResolution = Union[ResolvedTargets, TargetClarify, TargetMiss, ScopeError]


# ---------------------------------------------------------------------------
# Protocol resolution (per MeasurementFilter)
# ---------------------------------------------------------------------------


@dataclass
class ProtocolCandidate:
    """Per §9 clarify payload: chips show n_runs / last_run / n_compounds."""

    id: str
    name: str
    n_runs: int
    last_run: Optional[datetime.datetime]
    n_compounds: int


@dataclass
class ResolvedProtocol:
    protocol: Protocol
    query: str


@dataclass
class ProtocolClarify:
    query: str
    candidates: List[ProtocolCandidate]
    # The index of the MeasurementFilter this clarify applies to, so the
    # UI knows which filter to pin. The overall request-body field path is
    # `measurement_filters[<filter_index>].protocol_id`.
    filter_index: int = 0
    field: str = FIELD_PROTOCOL_HINT


@dataclass
class ProtocolMiss:
    query: str
    suggestions: List[ProtocolCandidate]
    filter_index: int = 0
    field: str = FIELD_PROTOCOL_HINT


ProtocolResolution = Union[ResolvedProtocol, ProtocolClarify, ProtocolMiss]


# ---------------------------------------------------------------------------
# Row-level evaluation (metric + unit tri-state). Evaluator module consumes
# these — public because tests over the evaluator also do.
# ---------------------------------------------------------------------------


@dataclass
class RowMatched:
    value: float
    row_unit: Optional[str]
    value_in_query_unit: float


@dataclass
class RowFiltered:
    reason: str  # one of FILTER_*


@dataclass
class RowExcluded:
    reason: str  # one of EXCLUDE_*


RowOutcome = Union[RowMatched, RowFiltered, RowExcluded]


# ---------------------------------------------------------------------------
# Executor output — CompoundSelection
# ---------------------------------------------------------------------------


@dataclass
class CompoundSelection:
    """Result of intersecting all MeasurementFilters against the target-
    scoped compound set. View layer converts this into a redirect URL to
    `/assays/aggregate`."""

    compound_formatted_ids: List[str]       # ordered by reg_number for stable URLs
    target_names: List[str]                  # what goes into `targets=` (usually 1)
    protocol_names: List[str]                # what goes into `protocols=` (0+ — only filter-used protocols)
    n_matched: int                           # == len(compound_formatted_ids)
    scope_sentence: str                      # human echo-back


@dataclass
class SpecError:
    """Executor-level spec-validation failure (e.g. filter with threshold
    but no metric). Distinct from ScopeError which is specifically about
    the at-least-one-target §6.5 rule."""

    field: str
    message: str


ExecutionResult = Union[
    CompoundSelection,
    TargetClarify,
    TargetMiss,
    ScopeError,
    ProtocolClarify,
    ProtocolMiss,
    SpecError,
]


# ---------------------------------------------------------------------------
# LLM prompt-parse result
# ---------------------------------------------------------------------------


@dataclass
class NotAQuery:
    """LLM-emitted signal that the user's prompt isn't a compound selection."""

    reason: str


@dataclass
class ParseError:
    """LLM output failed downstream validation (malformed JSON, schema
    violation, etc.). Distinct from NotAQuery — this is our glue-code
    failure, not the user's."""

    message: str
    raw: Optional[str] = None


PromptParseResult = Union[CompoundSelector, NotAQuery, ParseError]
