"""QuerySpec and resolver result dataclasses.

The shape mirrors §5 of NLP_QUERY_PROPOSAL.md. Slices 1–2 populate the
target fields and protocol hint; the rest are typed here so callers can
start building specs by hand and so future slices add behaviour without
reshaping the schema.
"""

from __future__ import annotations

import datetime
from dataclasses import dataclass, field
from typing import List, Optional, Union

from compounds.assays.models import Protocol
from compounds.registry.models import Target


# Scope-kind string constants (derived from which target fields are populated;
# see §6.5). Strings keep serialisation trivial — no enum-in-json dance.
SCOPE_BOTH_SAME = "both_same"
SCOPE_REG_ONLY = "reg_only"
SCOPE_ASSAY_ONLY = "assay_only"
SCOPE_CROSS = "cross"


# Field-name constants — used by the orchestrator to tag a Clarify/Miss so
# the UI knows which picker to show (see §9 `field` key in the JSON).
FIELD_REGISTRATION_TARGET = "registration_target_as_typed"
FIELD_ASSAY_TARGET = "assay_target_as_typed"
FIELD_PROTOCOL_HINT = "protocol_hint"


@dataclass
class Threshold:
    op: str
    value: float
    unit: Optional[str] = None


@dataclass
class QuerySpec:
    registration_target_as_typed: Optional[str] = None
    assay_target_as_typed: Optional[str] = None
    protocol_hint: Optional[str] = None
    metric: Optional[str] = None
    threshold: Optional[Threshold] = None
    columns: Optional[str] = None
    columns_explicit: Optional[List[str]] = None


# ---------------------------------------------------------------------------
# Target-field resolution (single field)
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
    field: str = ""  # set by the orchestrator to FIELD_REGISTRATION_TARGET / FIELD_ASSAY_TARGET


@dataclass
class TargetMiss:
    query: str
    suggestions: List[TargetCandidate]
    field: str = ""


TargetResolution = Union[ResolvedTarget, TargetClarify, TargetMiss]


# ---------------------------------------------------------------------------
# Two-target orchestrator results (§6.5 scope derivation + Q21 sequential)
# ---------------------------------------------------------------------------


@dataclass
class ResolvedTargets:
    registration: Optional[Target]
    assay: Optional[Target]
    scope_kind: str  # one of SCOPE_* above


@dataclass
class ScopeError:
    """At-least-one-target rule violated (§6.5 last row)."""

    message: str


TargetsResolution = Union[ResolvedTargets, TargetClarify, TargetMiss, ScopeError]


# ---------------------------------------------------------------------------
# Protocol resolution (§6.2)
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
    field: str = FIELD_PROTOCOL_HINT


@dataclass
class ProtocolMiss:
    query: str
    suggestions: List[ProtocolCandidate]
    field: str = FIELD_PROTOCOL_HINT


ProtocolResolution = Union[ResolvedProtocol, ProtocolClarify, ProtocolMiss]


# ---------------------------------------------------------------------------
# Row-level evaluation (§6.3 metric, §6.4 units tri-state)
# ---------------------------------------------------------------------------

# Filter reasons — a row didn't match the filter for expected reasons (the
# fit targeted a different KPI, the status isn't valid, or the value simply
# doesn't satisfy the threshold). These are *silent* — not footer-noted.
FILTER_INVALID_STATUS = "invalid_status"
FILTER_KPI_MISMATCH = "kpi_mismatch"
FILTER_VALUE_NOT_NUMERIC = "value_not_numeric"
FILTER_THRESHOLD_NOT_MET = "threshold_not_met"

# Exclusion reasons — a row couldn't be evaluated due to data-quality
# issues. These are the footer-count categories surfaced in §6.4.
EXCLUDE_UNIT_UNKNOWN = "unit_unknown"              # row's kpi_unit absent / empty / unrecognised
EXCLUDE_UNIT_TYPE_MISMATCH = "unit_type_mismatch"  # row unitless vs query has real unit, or vice versa
EXCLUDE_UNIT_INCOMPATIBLE = "unit_incompatible"    # different unit family (e.g. nM vs min)
EXCLUDE_QUERY_MISSING_UNIT = "query_missing_unit"  # query has no unit but row's unit is real


@dataclass
class RowMatched:
    """Row passes the filter. Values expressed in both row-native and query-native units."""

    value: float
    row_unit: Optional[str]          # normalized; None if row is genuinely unitless
    value_in_query_unit: float       # == value when no unit conversion happened


@dataclass
class RowFiltered:
    reason: str  # one of FILTER_*


@dataclass
class RowExcluded:
    reason: str  # one of EXCLUDE_*


RowOutcome = Union[RowMatched, RowFiltered, RowExcluded]


# ---------------------------------------------------------------------------
# Executor payload (§7 execution, §6.6 columns, §10 echo-back)
# ---------------------------------------------------------------------------

# Column-preset tokens the LLM may emit for QuerySpec.columns (§6.6).
COL_PRESET_PHYS_CHEM = "phys_chem"
COL_PRESET_LIPINSKI = "lipinski"

# Canonical ordered field sets on MolecularProperties. Order is the render
# order in the payload; the frontend is free to override.
PHYS_CHEM_COLUMNS: tuple = (
    "molecular_weight",
    "heavy_atom_count",
    "hbd",
    "hba",
    "clogp",
    "tpsa",
    "rotatable_bonds",
    "fraction_sp3",
)
LIPINSKI_COLUMNS: tuple = ("molecular_weight", "hbd", "hba", "clogp")


@dataclass
class TableRow:
    """One matched measurement. Per-measurement, not per-compound — the
    frontend may aggregate (geomean etc.) if desired."""

    compound_id: str                               # UUID string
    formatted_id: str                              # e.g. "NCL-00026042"
    smiles: Optional[str]
    value: float                                   # KPI value in the row's own unit
    value_unit: Optional[str]                      # None if the row is genuinely unitless
    value_in_query_unit: float                     # threshold-comparable — == value when no conversion happened
    properties: dict                               # {phys_chem_key: float | None}


@dataclass
class TablePayload:
    rows: List[TableRow]
    property_columns: List[str]                    # ordered phys-chem keys shown on each row
    footer_excluded: dict                          # {EXCLUDE_* reason: count} — surfaced in UI footer (§6.4)
    filtered_silent: dict                          # {FILTER_* reason: count} — debug/telemetry, not UI-surfaced
    scope_sentence: str                            # human echo-back (§10)


@dataclass
class SpecError:
    """Executor-level spec validation failure (missing metric, unknown column, …)."""

    field: str
    message: str


ExecutionResult = Union[
    TablePayload,
    TargetClarify,
    TargetMiss,
    ScopeError,
    ProtocolClarify,
    ProtocolMiss,
    SpecError,
]


# ---------------------------------------------------------------------------
# LLM prompt-parse result (§8 — the LLM emits either a QuerySpec or this)
# ---------------------------------------------------------------------------


@dataclass
class NotAQuery:
    """LLM-emitted signal that the user's prompt isn't a tabular query."""

    reason: str


@dataclass
class ParseError:
    """LLM output failed downstream validation (malformed JSON, schema violation,
    both branches empty, …). Distinct from NotAQuery — this is *our* fault, not
    the user's. Slice 6 view will map to a 502-ish response; the caller is
    expected to retry or report rather than show a user-friendly message."""

    message: str
    raw: Optional[str] = None


PromptParseResult = Union[QuerySpec, NotAQuery, ParseError]
