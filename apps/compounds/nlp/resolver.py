"""Deterministic target/protocol resolution for the NLP query subsystem.

Implements §6.1 (target), §6.2 (protocol), and §6.5 (two-target scope
derivation) of NLP_QUERY_PROPOSAL.md.

Target (§6.1):
- Matching pool per Target: normalize(name) ∪ {normalize(g.symbol)} ∪
  {normalize(alias)} ∪ {normalize(g.name)} for each linked Gene.
- Exact-match pass; on miss, allow Levenshtein distance 1 **only if the
  normalized query is ≥4 characters**.
- Name-only fallback for targets without Gene links.
- Multiple hits → TargetClarify; no hits → TargetMiss with top-5 suggestions.

Protocol (§6.2):
- Scope = protocols that have been run on compounds registered to the
  resolved registration target (and/or against the resolved assay target).
- Tokenize name + query (lowercase, strip punctuation); require all query
  tokens to appear in the protocol name. Rank ties by fewer extra tokens
  then recency of last run. Ambiguity is the common case.

Scope (§6.5):
- Derive `scope_kind` from which target fields are populated.
- Q21: sequential clarification — reg first, then assay.

Pure Python. No LLM. No network. Unit-testable with pytest-django fixtures.
"""

from __future__ import annotations

import re
from dataclasses import replace
from difflib import SequenceMatcher
from typing import Dict, Iterable, List, Optional, Tuple

from django.db.models import Count, Max, QuerySet

from compounds.assays.models import Assay, DataSeries, Protocol
from compounds.registry.models import Gene, Target
from compounds.utils import normalize_ref as normalize

from .spec import (
    FIELD_ASSAY_TARGET,
    FIELD_PROTOCOL_HINT,
    FIELD_REGISTRATION_TARGET,
    SCOPE_ASSAY_ONLY,
    SCOPE_BOTH_SAME,
    SCOPE_CROSS,
    SCOPE_REG_ONLY,
    ProtocolCandidate,
    ProtocolClarify,
    ProtocolMiss,
    ProtocolResolution,
    QuerySpec,
    ResolvedProtocol,
    ResolvedTarget,
    ResolvedTargets,
    ScopeError,
    TargetCandidate,
    TargetClarify,
    TargetMiss,
    TargetResolution,
    TargetsResolution,
)


_FUZZY_MIN_QUERY_LEN = 4
_MISS_SUGGESTION_COUNT = 5


def _edit_distance_le_1(a: str, b: str) -> bool:
    """True iff Levenshtein distance between a and b is at most 1."""
    if a == b:
        return True
    la, lb = len(a), len(b)
    if abs(la - lb) > 1:
        return False
    if la == lb:
        return sum(x != y for x, y in zip(a, b)) == 1
    if la < lb:
        a, b = b, a
        la, lb = lb, la
    # la == lb + 1: b must be obtainable by deleting one char from a
    for i in range(la):
        if a[:i] + a[i + 1 :] == b:
            return True
    return False


def _target_candidate(target: Target, genes: Iterable[Gene]) -> TargetCandidate:
    return TargetCandidate(
        id=str(target.id),
        name=target.name,
        gene_symbols=sorted(g.symbol for g in genes if g.symbol),
    )


def _build_pool(target: Target, genes: List[Gene]) -> Dict[str, str]:
    """Return normalized-entry → matched_via-tag for a single target.

    Name-only fallback: targets with no genes get just {normalize(name)}.
    """
    pool: Dict[str, str] = {}
    name_norm = normalize(target.name)
    if name_norm:
        pool[name_norm] = "name"
    for g in genes:
        if g.symbol:
            pool.setdefault(normalize(g.symbol), "gene_symbol")
        for alias in g.aliases or []:
            alias_norm = normalize(alias)
            if alias_norm:
                pool.setdefault(alias_norm, "gene_alias")
        if g.name:
            name_norm = normalize(g.name)
            if name_norm:
                pool.setdefault(name_norm, "gene_name")
    pool.pop("", None)
    return pool


def _iter_target_pools() -> List[Tuple[Target, List[Gene], Dict[str, str]]]:
    """Load every Target with prefetched genes and build its matching pool."""
    out: List[Tuple[Target, List[Gene], Dict[str, str]]] = []
    for t in Target.objects.prefetch_related("genes").all():
        genes = list(t.genes.all())
        out.append((t, genes, _build_pool(t, genes)))
    return out


def resolve_target(as_typed: str, pinned_id: Optional[str] = None) -> TargetResolution:
    """Resolve a user-typed target string against the Gene-hydrated pool.

    If ``pinned_id`` is supplied (clarify continuation), look up the Target
    directly and skip fuzzy matching — the user has already disambiguated.
    """
    if pinned_id:
        target = Target.objects.filter(pk=pinned_id).first()
        if target is None:
            return TargetMiss(query=as_typed or "", suggestions=[])
        return ResolvedTarget(target=target, matched_via="pinned_id", query=as_typed or "")

    norm_query = normalize(as_typed or "")
    if not norm_query:
        return TargetMiss(query=as_typed or "", suggestions=[])

    target_pools = _iter_target_pools()

    exact: List[Tuple[Target, List[Gene], str]] = []
    for target, genes, pool in target_pools:
        via = pool.get(norm_query)
        if via is not None:
            exact.append((target, genes, via))

    if len(exact) == 1:
        target, genes, via = exact[0]
        return ResolvedTarget(target=target, matched_via=via, query=as_typed)
    if len(exact) > 1:
        return TargetClarify(
            query=as_typed,
            candidates=[_target_candidate(t, genes) for t, genes, _ in exact],
        )

    if len(norm_query) >= _FUZZY_MIN_QUERY_LEN:
        fuzzy: List[Tuple[Target, List[Gene], str]] = []
        for target, genes, pool in target_pools:
            for entry, via in pool.items():
                if entry != norm_query and _edit_distance_le_1(norm_query, entry):
                    fuzzy.append((target, genes, f"{via}_fuzzy"))
                    break
        if len(fuzzy) == 1:
            target, genes, via = fuzzy[0]
            return ResolvedTarget(target=target, matched_via=via, query=as_typed)
        if len(fuzzy) > 1:
            return TargetClarify(
                query=as_typed,
                candidates=[_target_candidate(t, genes) for t, genes, _ in fuzzy],
            )

    scored: List[Tuple[float, Target, List[Gene]]] = []
    for target, genes, pool in target_pools:
        if not pool:
            continue
        best = max(
            SequenceMatcher(None, norm_query, entry).ratio() for entry in pool
        )
        scored.append((best, target, genes))
    scored.sort(key=lambda row: row[0], reverse=True)
    suggestions = [
        _target_candidate(t, genes)
        for _, t, genes in scored[:_MISS_SUGGESTION_COUNT]
    ]
    return TargetMiss(query=as_typed, suggestions=suggestions)


# ---------------------------------------------------------------------------
# Two-target orchestrator (§6.5 + Q21)
# ---------------------------------------------------------------------------


def _derive_scope_kind(reg: Optional[Target], assay: Optional[Target]) -> str:
    if reg is not None and assay is not None:
        return SCOPE_BOTH_SAME if reg.pk == assay.pk else SCOPE_CROSS
    if reg is not None:
        return SCOPE_REG_ONLY
    return SCOPE_ASSAY_ONLY


def resolve_targets(spec: QuerySpec) -> TargetsResolution:
    """Resolve both target fields, with sequential clarification per Q21.

    Returns the first unresolved field's clarify/miss if any, otherwise a
    fully-populated ResolvedTargets with a derived scope_kind. Returns
    ScopeError if neither target field is set (§6.5 last row).
    """
    reg_typed = (spec.registration_target_as_typed or "").strip()
    assay_typed = (spec.assay_target_as_typed or "").strip()

    if (
        not reg_typed and not assay_typed
        and not spec.registration_target_id and not spec.assay_target_id
    ):
        return ScopeError(
            message="At least one of registration_target_as_typed or "
            "assay_target_as_typed must be set."
        )

    reg: Optional[Target] = None
    if reg_typed or spec.registration_target_id:
        r = resolve_target(reg_typed, pinned_id=spec.registration_target_id)
        if isinstance(r, (TargetClarify, TargetMiss)):
            return replace(r, field=FIELD_REGISTRATION_TARGET)
        reg = r.target

    assay: Optional[Target] = None
    if assay_typed or spec.assay_target_id:
        r = resolve_target(assay_typed, pinned_id=spec.assay_target_id)
        if isinstance(r, (TargetClarify, TargetMiss)):
            return replace(r, field=FIELD_ASSAY_TARGET)
        assay = r.target

    return ResolvedTargets(
        registration=reg,
        assay=assay,
        scope_kind=_derive_scope_kind(reg, assay),
    )


# ---------------------------------------------------------------------------
# Protocol resolution (§6.2)
# ---------------------------------------------------------------------------


_TOKEN_BOUNDARY = re.compile(r"[^a-z0-9]+")


def tokenize(s: str) -> set:
    """Lowercase and split on runs of non-alphanumerics.

    Exported for testing; the protocol match rule uses set-containment
    over the result.
    """
    if not s:
        return set()
    return {tok for tok in _TOKEN_BOUNDARY.split(s.lower()) if tok}


def _scope_assays(rt: ResolvedTargets) -> QuerySet:
    """Build the Assay queryset matching the derived scope (§6.5 table)."""
    qs = Assay.objects.all()
    if rt.assay is not None:
        qs = qs.filter(target=rt.assay)
    if rt.registration is not None:
        qs = qs.filter(data_series__compound__target=rt.registration)
    return qs.distinct()


def _protocol_candidate(protocol: Protocol, scoped_assays: QuerySet) -> ProtocolCandidate:
    assays_of = scoped_assays.filter(protocol=protocol)
    stats = assays_of.aggregate(
        n_runs=Count("id", distinct=True),
        last_run=Max("created_at"),
    )
    n_compounds = (
        DataSeries.objects
        .filter(assay__in=assays_of, compound__isnull=False)
        .values("compound")
        .distinct()
        .count()
    )
    return ProtocolCandidate(
        id=str(protocol.id),
        name=protocol.name,
        n_runs=stats["n_runs"] or 0,
        last_run=stats["last_run"],
        n_compounds=n_compounds,
    )


def resolve_protocol(
    hint: str,
    resolved_targets: ResolvedTargets,
    pinned_id: Optional[str] = None,
) -> ProtocolResolution:
    """Resolve a protocol hint against the scope-filtered protocol pool.

    If ``pinned_id`` is supplied (clarify continuation), look up the Protocol
    directly — we still verify it falls within the resolved scope so a stale
    UUID from a different target's clarify response can't slip through.
    """
    if pinned_id:
        scoped_assays = _scope_assays(resolved_targets)
        if not scoped_assays.filter(protocol_id=pinned_id).exists():
            return ProtocolMiss(query=hint or "", suggestions=[])
        protocol = Protocol.objects.filter(pk=pinned_id).first()
        if protocol is None:
            return ProtocolMiss(query=hint or "", suggestions=[])
        return ResolvedProtocol(protocol=protocol, query=hint or "")

    query_tokens = tokenize(hint or "")
    if not query_tokens:
        return ProtocolMiss(query=hint or "", suggestions=[])

    scoped_assays = _scope_assays(resolved_targets)
    protocol_ids = list(
        scoped_assays.values_list("protocol_id", flat=True).distinct()
    )
    if not protocol_ids:
        return ProtocolMiss(query=hint, suggestions=[])

    protocols = list(Protocol.objects.filter(id__in=protocol_ids))

    # Pass 1: protocols whose tokens ⊇ query tokens.
    matches: List[Tuple[Protocol, int]] = []
    for p in protocols:
        p_tokens = tokenize(p.name)
        if query_tokens.issubset(p_tokens):
            extra = len(p_tokens - query_tokens)
            matches.append((p, extra))

    if len(matches) == 1:
        return ResolvedProtocol(protocol=matches[0][0], query=hint)

    if len(matches) > 1:
        # Rank by (fewer extra tokens, then recency of last run desc).
        candidates_with_extra = [
            (extra, _protocol_candidate(p, scoped_assays))
            for p, extra in matches
        ]
        candidates_with_extra.sort(
            key=lambda row: (
                row[0],
                -(row[1].last_run.timestamp() if row[1].last_run else 0),
            )
        )
        return ProtocolClarify(
            query=hint,
            candidates=[c for _, c in candidates_with_extra],
        )

    # Miss — rank in-scope protocols by token overlap, take top 5.
    scored: List[Tuple[int, Protocol]] = []
    for p in protocols:
        overlap = len(query_tokens & tokenize(p.name))
        if overlap > 0:
            scored.append((overlap, p))
    scored.sort(key=lambda row: row[0], reverse=True)
    suggestions = [
        _protocol_candidate(p, scoped_assays)
        for _, p in scored[:_MISS_SUGGESTION_COUNT]
    ]
    return ProtocolMiss(query=hint, suggestions=suggestions)
