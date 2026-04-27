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
# Neither target field is set — the selector is narrowed only by its
# other predicates (scaffold, registered_by, registered_date_range,
# measurement filters). Runs across the whole registry. Prompts like
# *"all pyrimidines"*, *"compounds registered by Suzannah"*, or *"anything
# with HTRF IC50 < 10 nM"* land here. The executor still requires at
# least one non-target predicate so we don't silently dump the whole
# compound list for a fully-empty selector.
SCOPE_UNSCOPED = "unscoped"

FIELD_REGISTRATION_TARGET = "registration_target_as_typed"
FIELD_ASSAY_TARGET = "assay_target_as_typed"
FIELD_PROTOCOL_HINT = "protocol_hint"
FIELD_REGISTERED_BY = "registered_by_as_typed"
FIELD_ASSAYED_BY = "assayed_by_as_typed"
FIELD_SCAFFOLD_HINT = "scaffold_hint"
FIELD_COMPOUND_REF = "compound_ref"
FIELD_METRIC = "metric"

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
    # Name of the person who ran the assay — filters on Assay.created_by.
    # Resolved via the user resolver; clarifies when ambiguous.
    assayed_by_as_typed: Optional[str] = None
    # Pinned IDs on clarify continuation — LLM never emits these.
    protocol_id: Optional[str] = None
    assayed_by_id: Optional[str] = None


@dataclass
class AssaySelector:
    """Root of the assay-selection query family — "CDK4 assays carried
    out last week", "Assays conducted by Alice against ARd", etc.

    Distinct from CompoundSelector because the output is a list of
    **Assay** records (per §19.7 follow-up), not compounds. The LLM
    emits either an AssaySelector OR a CompoundSelector, not both.
    """

    target_as_typed: Optional[str] = None
    protocol_hint: Optional[str] = None
    date_range: Optional[DateRange] = None      # Assay.created_at
    created_by_as_typed: Optional[str] = None   # "conducted by", "run by", "assayed by"
    # Scaffold filter — "HTRF assays on pyrimidines": narrow the output
    # to assays that ran on compounds containing ALL named scaffolds.
    # Entries are names typed by the user; resolved via the curated
    # substructures catalog and applied with RDKit HasSubstructMatch.
    scaffold_hints: List[str] = field(default_factory=list)

    # Pinnings from clarify continuation — the view injects these.
    target_id: Optional[str] = None
    protocol_id: Optional[str] = None
    created_by_id: Optional[str] = None
    scaffold_ids: List[str] = field(default_factory=list)    # pinned canonical names


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
    # Name of the person OR supplier who registered / made / synthesised /
    # supplied the compound — filters on Compound.registered_by (User FK)
    # or Compound.supplier (Supplier FK, possibly user-linked). Slice 16
    # broadened the resolver to be poly-source: a single typed phrase
    # ("Martin", "Enamine", "Alice Jones") is matched against both pools,
    # de-duped when a Supplier links to a User of the same name.
    registered_by_as_typed: Optional[str] = None
    # Substructure / scaffold filter — "ARd compounds containing pyrimidine".
    # Multiple hints are ANDed (each must match). Resolved via the curated
    # substructures catalog + RDKit HasSubstructMatch at query time.
    scaffold_hints: List[str] = field(default_factory=list)
    # Compound-ID pins — *additive* (UNION), not a filter. The user names
    # specific compounds by ID ("NCL-00026007", "compound 26007", "the
    # lead compound NCL26007") and they appear in the final selection
    # alongside whatever the other predicates select. Useful for
    # "compound 26007 and compounds made since January" — the lead
    # compound stays available for comparison alongside the recent set.
    # Strings are emitted verbatim by the LLM and parsed by the
    # resolver via ``compounds.formatting.extract_reg_number``, which
    # handles every prefix variant (NCL-00026007 / NCL26007 / NCL 26007
    # / 26007 / NCL000-26007).
    compound_refs_as_typed: List[str] = field(default_factory=list)

    # Pinnings from clarify continuation — the view injects these, LLM doesn't.
    registration_target_id: Optional[str] = None
    assay_target_id: Optional[str] = None
    registered_by_id: Optional[str] = None         # User pk pinned by clarify
    # Slice 16: when the clarify chip pinned a Supplier (rather than a
    # User), the supplier id round-trips here. Mutually exclusive with
    # registered_by_id at any one time, but the schema allows either.
    registered_by_supplier_id: Optional[str] = None
    scaffold_ids: List[str] = field(default_factory=list)    # pinned canonical names


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
# User resolution (compound registrant / assay creator)
# ---------------------------------------------------------------------------


@dataclass
class UserCandidate:
    """Thin JSON-serialisable view of a User for clarify / miss payloads.

    ``display`` is the best human-readable rendering we can offer in a
    chip — falls through first_name+last_name → display_name → email →
    username. ``email`` is always shown as the secondary disambiguator
    in the picker since names collide (Alice J, Alice J).

    ``kind`` is a discriminator (slice 16) — pickers under
    ``registered_by_as_typed`` may also surface ``SupplierCandidate``s,
    so the frontend dispatches on this field to pick the right
    pin-field on continuation (``registered_by_id`` vs
    ``registered_by_supplier_id``).
    """

    id: str
    display: str
    email: Optional[str] = None
    n_compounds: int = 0          # how many compounds this user registered (hint for registered-by picker)
    n_assays: int = 0             # how many assays this user created (hint for assayed-by picker)
    kind: str = "user"


@dataclass
class ResolvedUser:
    # Type-any because importing django.contrib.auth.User here would
    # pull Django models into this spec module — keep it agnostic.
    user: object
    matched_via: str
    query: str


@dataclass
class UserClarify:
    query: str
    # Slice 16/18: registered_by clarifies may mix User, Supplier, and
    # Union candidates. The kind field on each candidate disambiguates
    # — frontend pins to registered_by_id (user), registered_by_supplier_id
    # (supplier), or BOTH (union) accordingly.
    candidates: List[Union["UserCandidate", "SupplierCandidate", "UnionCandidate"]]
    field: str = ""                # FIELD_REGISTERED_BY / FIELD_ASSAYED_BY — set by caller
    filter_index: int = 0          # only meaningful for FIELD_ASSAYED_BY


@dataclass
class UserMiss:
    query: str
    suggestions: List[Union["UserCandidate", "SupplierCandidate", "UnionCandidate"]]
    field: str = ""
    filter_index: int = 0


UserResolution = Union[
    ResolvedUser, "ResolvedSupplier", "ResolvedUnion", UserClarify, UserMiss,
]


# ---------------------------------------------------------------------------
# Supplier (slice 16) — sibling of UserCandidate / ResolvedUser. The
# registered_by resolver is poly-source: a typed phrase like "Enamine"
# resolves to a Supplier, while "Alice Jones" resolves to a User. A
# single Clarify can carry candidates of either kind so picker chips
# render uniformly. Note: when a Supplier is user-linked (Supplier.user
# == User), the resolver de-dupes — that's one entity from the
# chemist's perspective.
# ---------------------------------------------------------------------------


@dataclass
class SupplierCandidate:
    """Thin JSON-serialisable view of a Supplier for clarify / miss
    payloads. ``kind`` is the discriminator the frontend uses to pick
    the right pin field on continuation."""

    id: str
    name: str
    kind: str = "supplier"        # discriminator vs UserCandidate (kind="user")
    initials: Optional[str] = None
    n_compounds: int = 0
    is_user_linked: bool = False  # de-duped against User candidate when True


@dataclass
class ResolvedSupplier:
    # Type-any to avoid importing the Supplier model into spec.
    supplier: object
    matched_via: str
    query: str


@dataclass
class UnionCandidate:
    """Slice 18: surfaces when a typed phrase matches BOTH a User and a
    same-named unlinked Supplier — e.g. "Hannah Stewart" who registered
    three compounds AND has a personal-Supplier row supplying two more.
    The chemist's intent is the *union* of both sets, not picking one
    role. The frontend renders this as a single chip ("3 registered +
    2 supplied") and the executor's Q-filter ORs all three paths."""

    id: str                          # composite, opaque to the wire
    display: str
    user_id: str
    supplier_id: str
    n_compounds_user: int = 0        # how many compounds this user registered
    n_compounds_supplier: int = 0    # how many compounds this supplier supplied
    kind: str = "union"              # discriminator vs UserCandidate / SupplierCandidate


@dataclass
class ResolvedUnion:
    """A merged User+Supplier resolution where the same-named entities
    are treated as one entity for filtering purposes. The executor
    ORs Q(registered_by=user) | Q(supplier__user=user) | Q(supplier=supplier)."""

    user: object
    supplier: object
    matched_via: str
    query: str


SupplierResolution = Union[ResolvedSupplier, UserClarify, UserMiss]


# ---------------------------------------------------------------------------
# Scaffold / substructure resolution (slice 13)
# ---------------------------------------------------------------------------


@dataclass
class ScaffoldCandidate:
    """Thin view of a scaffold catalog entry for clarify / miss payloads.
    The ``id`` is the canonical name — stable across the catalog and
    directly usable by RDKit-match-side lookups."""

    id: str                # == canonical name
    name: str              # canonical (same as id; separate field for
                           # consistency with the other candidate types)
    aliases: List[str] = field(default_factory=list)
    smarts: str = ""       # exposed so a future SVG preview can render


@dataclass
class ResolvedScaffold:
    # Type-any to keep spec.py free of the substructures module import
    # cycle (spec is the common dependency).
    scaffold: object
    matched_via: str
    query: str


@dataclass
class ScaffoldClarify:
    query: str
    candidates: List[ScaffoldCandidate]
    scaffold_index: int = 0        # which entry of scaffold_hints this clarifies
    field: str = ""                # FIELD_SCAFFOLD_HINT


@dataclass
class ScaffoldMiss:
    query: str
    suggestions: List[ScaffoldCandidate]
    scaffold_index: int = 0
    field: str = ""


ScaffoldResolution = Union[ResolvedScaffold, ScaffoldClarify, ScaffoldMiss]


# ---------------------------------------------------------------------------
# Compound-ID pinning resolution — no clarify because compound IDs are
# deterministic. A typed reference either resolves to a single Compound
# or doesn't exist; the user revises their prompt on a miss.
# ---------------------------------------------------------------------------


@dataclass
class ResolvedCompound:
    # Type-any on the model to keep spec.py free of a registry import cycle.
    compound: object
    query: str               # the original typed reference, for round-trip


@dataclass
class CompoundMiss:
    """Typed reference didn't parse as a compound ID, or referred to a
    compound that doesn't exist. The view surfaces this as a generic
    miss; clients should advise the user to check the ID and retry."""

    query: str
    ref_index: int = 0       # which entry of compound_refs_as_typed missed
    field: str = ""          # FIELD_COMPOUND_REF


CompoundResolution = Union[ResolvedCompound, CompoundMiss]


# ---------------------------------------------------------------------------
# Metric (KPI) miss (slice 15) — diagnostic surfacing of empty selections
# caused by a typed metric that's not present in any row's stored KPI.
# No clarify dataclass: the executor either finds the metric in scope and
# proceeds, or returns a Miss listing what KPIs ARE there.
# ---------------------------------------------------------------------------


@dataclass
class MetricMiss:
    """No row in the filter's scope has a KPI matching the user's typed
    ``metric`` (lenient match: case + punctuation insensitive). The
    ``available_metrics`` field carries the KPIs that ARE recorded in
    scope, so the response answers itself ("did you mean pIC50?")."""

    query: str
    available_metrics: List[str]
    filter_index: int = 0
    field: str = ""


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
class AssaySelection:
    """Result of filtering Assays by target/protocol/date/creator. View
    layer converts this into a redirect URL to `/assays?ids=...`."""

    assay_ids: List[str]                     # Assay UUIDs (for the /assays?ids= URL)
    n_matched: int
    target_names: List[str]                  # for scope sentence (usually 1)
    protocol_names: List[str]                # for scope sentence (0 or 1)
    scope_sentence: str


@dataclass
class SpecError:
    """Executor-level spec-validation failure (e.g. filter with threshold
    but no metric). Distinct from ScopeError which is specifically about
    the at-least-one-target §6.5 rule."""

    field: str
    message: str


ExecutionResult = Union[
    CompoundSelection,
    AssaySelection,
    TargetClarify,
    TargetMiss,
    ScopeError,
    ProtocolClarify,
    ProtocolMiss,
    UserClarify,
    UserMiss,
    ScaffoldClarify,
    ScaffoldMiss,
    CompoundMiss,
    MetricMiss,
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


PromptParseResult = Union[CompoundSelector, AssaySelector, NotAQuery, ParseError]
