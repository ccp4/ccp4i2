"""Deterministic target + protocol resolution for the compound-selector path.

Implements §6.1 (target), §6.2 (protocol), and §6.5 (two-target scope
derivation) of NLP_QUERY_PROPOSAL.md.

**Pivoted 2026-04-23** from QuerySpec (single protocol) to CompoundSelector
(list of MeasurementFilters, each with its own protocol_hint). The target
resolver is unchanged; the protocol resolver gains a `filter_index`
parameter so clarify/miss responses can tell the UI which filter to pin.

Target (§6.1):
- Matching pool per Target: normalize(name) ∪ {normalize(g.symbol)} ∪
  {normalize(alias)} ∪ {normalize(g.name)} for each linked Gene.
- Exact pass, then distance-1 fuzzy (query length ≥ 4 only).
- Name-only fallback for targets without Gene links.

Protocol (§6.2):
- Scope: protocols run on compounds registered to the resolved registration
  target and/or against the resolved assay target.
- Token match: all query tokens must appear in the protocol name (set
  containment). Ranked by (fewer extras, then recency of last run).
- Pinning: when `pinned_id` is set, skip fuzzy match and look up directly;
  still scope-validate so a stale UUID from another target can't slip in.

Pure Python. No LLM. No network. Unit-testable with pytest-django fixtures.
"""

from __future__ import annotations

import re
from dataclasses import replace
from difflib import SequenceMatcher
from typing import Dict, Iterable, List, Optional, Tuple

from django.contrib.auth import get_user_model
from django.db.models import Count, Max, Q, QuerySet

from compounds.assays.models import Assay, DataSeries, Protocol
from compounds.registry.models import Gene, Target
from compounds.utils import normalize_ref as normalize

from .spec import (
    FIELD_ASSAY_TARGET,
    FIELD_ASSAYED_BY,
    FIELD_PROTOCOL_HINT,
    FIELD_REGISTERED_BY,
    FIELD_REGISTRATION_TARGET,
    FIELD_SCAFFOLD_HINT,
    SCOPE_ASSAY_ONLY,
    SCOPE_BOTH_SAME,
    SCOPE_CROSS,
    SCOPE_REG_ONLY,
    SCOPE_UNSCOPED,
    CompoundSelector,
    ProtocolCandidate,
    ProtocolClarify,
    ProtocolMiss,
    ProtocolResolution,
    ResolvedProtocol,
    ResolvedTarget,
    ResolvedTargets,
    ResolvedUser,
    ScopeError,
    TargetCandidate,
    TargetClarify,
    TargetMiss,
    TargetResolution,
    TargetsResolution,
    UserCandidate,
    UserClarify,
    UserMiss,
    UserResolution,
    ScaffoldCandidate,
    ScaffoldClarify,
    ScaffoldMiss,
    ScaffoldResolution,
    ResolvedScaffold,
)
from .substructures import (
    SCAFFOLDS,
    Scaffold,
    lookup_scaffold_by_id,
    lookup_scaffold_by_typed,
    lookup_scaffolds_by_substring,
)

User = get_user_model()


_FUZZY_MIN_QUERY_LEN = 4
_MISS_SUGGESTION_COUNT = 5
# Cap on how many chips we put in a Clarify picker — the UI wraps chips so
# long lists work, but beyond ~10 the picker becomes a scrolling exercise
# rather than a decision aid.
_CLARIFY_MAX_CANDIDATES = 10


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
    if assay is not None:
        return SCOPE_ASSAY_ONLY
    return SCOPE_UNSCOPED


def resolve_targets(selector: CompoundSelector) -> TargetsResolution:
    """Resolve both target fields on the selector, with sequential clarify
    per Q21. Returns the first unresolved field's clarify/miss if any,
    otherwise a ResolvedTargets with a derived scope_kind.

    Either-empty is permitted — the selector is simply unscoped w.r.t.
    target and relies on its other predicates (scaffold / user / date /
    measurement filter) to narrow. The executor enforces that at least
    one narrowing predicate exists when the scope is unscoped.
    """
    reg_typed = (selector.registration_target_as_typed or "").strip()
    assay_typed = (selector.assay_target_as_typed or "").strip()

    reg: Optional[Target] = None
    if reg_typed or selector.registration_target_id:
        r = resolve_target(reg_typed, pinned_id=selector.registration_target_id)
        if isinstance(r, (TargetClarify, TargetMiss)):
            return replace(r, field=FIELD_REGISTRATION_TARGET)
        reg = r.target

    assay: Optional[Target] = None
    if assay_typed or selector.assay_target_id:
        r = resolve_target(assay_typed, pinned_id=selector.assay_target_id)
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
    filter_index: int = 0,
) -> ProtocolResolution:
    """Resolve a protocol hint against the scope-filtered protocol pool.

    ``filter_index`` is tagged onto any Clarify/Miss result so the UI
    knows *which* measurement_filter in the selector triggered the picker.

    ``pinned_id`` short-circuits to a direct lookup (clarify continuation);
    still scope-validates.
    """
    if pinned_id:
        scoped_assays = _scope_assays(resolved_targets)
        if not scoped_assays.filter(protocol_id=pinned_id).exists():
            return ProtocolMiss(
                query=hint or "", suggestions=[], filter_index=filter_index,
            )
        protocol = Protocol.objects.filter(pk=pinned_id).first()
        if protocol is None:
            return ProtocolMiss(
                query=hint or "", suggestions=[], filter_index=filter_index,
            )
        return ResolvedProtocol(protocol=protocol, query=hint or "")

    query_tokens = tokenize(hint or "")
    if not query_tokens:
        # Empty hint is OK at the selector level (filter doesn't scope to a
        # protocol — "compounds with ANY valid IC50 measurement"). But for a
        # direct resolve_protocol call with an empty hint, treat as miss.
        return ProtocolMiss(
            query=hint or "", suggestions=[], filter_index=filter_index,
        )

    scoped_assays = _scope_assays(resolved_targets)
    protocol_ids = list(
        scoped_assays.values_list("protocol_id", flat=True).distinct()
    )
    if not protocol_ids:
        return ProtocolMiss(
            query=hint, suggestions=[], filter_index=filter_index,
        )

    protocols = list(Protocol.objects.filter(id__in=protocol_ids))

    # Per-protocol last_run for recency ranking — one query for the whole
    # scope rather than one per candidate.
    recency: Dict = {
        row["protocol_id"]: row["last_run"]
        for row in scoped_assays
        .values("protocol_id")
        .annotate(last_run=Max("created_at"))
    }

    def _recency_ts(p: Protocol) -> float:
        dt = recency.get(p.id)
        return dt.timestamp() if dt else 0.0

    norm_query = normalize(hint or "")

    # Bucket every in-scope protocol by match tier. A protocol lands in the
    # strongest tier it qualifies for:
    #   strict     — protocol tokens ⊇ query tokens (original §6.2 rule)
    #   substring  — normalised query is a substring of the normalised name
    #   overlap    — at least one shared token (weak but still relevant)
    # Anything else misses entirely.
    strict: List[Tuple[Protocol, int]] = []        # (p, extra_token_count)
    substring: List[Protocol] = []
    overlap: List[Tuple[Protocol, int]] = []       # (p, overlap_size)
    for p in protocols:
        p_tokens = tokenize(p.name)
        if query_tokens.issubset(p_tokens):
            strict.append((p, len(p_tokens - query_tokens)))
            continue
        p_norm = normalize(p.name)
        if norm_query and norm_query in p_norm:
            substring.append(p)
            continue
        shared = query_tokens & p_tokens
        if shared:
            overlap.append((p, len(shared)))

    # Single strict hit resolves without clarification (preserve §6.2 rule).
    if len(strict) == 1:
        return ResolvedProtocol(protocol=strict[0][0], query=hint)

    # Build the tiered clarify candidate list. Ordering: strict (ranked by
    # fewer extras then recency), then substring (recency), then overlap
    # (more shared tokens then recency). Cap at 10 so the picker stays
    # legible for very broad queries.
    tiered: List[Protocol] = []
    strict.sort(key=lambda row: (row[1], -_recency_ts(row[0])))
    tiered.extend(p for p, _ in strict)
    tiered.extend(sorted(substring, key=lambda p: -_recency_ts(p)))
    overlap.sort(key=lambda row: (-row[1], -_recency_ts(row[0])))
    tiered.extend(p for p, _ in overlap)
    tiered = tiered[: _CLARIFY_MAX_CANDIDATES]

    if tiered:
        return ProtocolClarify(
            query=hint,
            candidates=[_protocol_candidate(p, scoped_assays) for p in tiered],
            filter_index=filter_index,
        )

    # Truly nothing matched in any tier. Miss, but surface the most-recently-
    # used protocols in scope as suggestions — better than a dead end, and
    # lets the user eyeball what's actually available for the target.
    suggestions = sorted(protocols, key=lambda p: -_recency_ts(p))[
        : _MISS_SUGGESTION_COUNT
    ]
    return ProtocolMiss(
        query=hint,
        suggestions=[_protocol_candidate(p, scoped_assays) for p in suggestions],
        filter_index=filter_index,
    )


# ---------------------------------------------------------------------------
# User resolution (compound registrant / assay creator)
# ---------------------------------------------------------------------------


def _user_display(user) -> str:
    """Best human-readable rendering for a User — first_name+last_name,
    falling through to display_name / email / username. Keeps picker
    chips meaningful even when accounts are sparsely populated."""
    first = (user.first_name or "").strip()
    last = (user.last_name or "").strip()
    if first and last:
        return f"{first} {last}"
    if first:
        return first
    if last:
        return last
    profile = getattr(user, "profile", None)
    display_name = getattr(profile, "display_name", None) if profile else None
    if display_name:
        return display_name
    return user.email or user.username


def _user_candidate(user, *, n_compounds: int = 0, n_assays: int = 0) -> UserCandidate:
    return UserCandidate(
        id=str(user.pk),
        display=_user_display(user),
        email=user.email or None,
        n_compounds=n_compounds,
        n_assays=n_assays,
    )


def _user_pool(user) -> set:
    """Normalised strings a user's name/email/username should match."""
    pool: set = set()
    for field in ("first_name", "last_name", "username", "email"):
        value = getattr(user, field, None) or ""
        norm = normalize(value)
        if norm:
            pool.add(norm)
    # Full name "First Last" normalises to "firstlast" — users often type
    # names in this combined form.
    first = (user.first_name or "").strip()
    last = (user.last_name or "").strip()
    if first and last:
        combined = normalize(f"{first} {last}")
        if combined:
            pool.add(combined)
    # Email local-part too: "alice.jones@ncl.ac.uk" → also match "alice.jones"
    email = (user.email or "").strip()
    if "@" in email:
        local = email.split("@", 1)[0]
        norm = normalize(local)
        if norm:
            pool.add(norm)
    # UserProfile.display_name if present.
    profile = getattr(user, "profile", None)
    display_name = getattr(profile, "display_name", None) if profile else None
    if display_name:
        norm = normalize(display_name)
        if norm:
            pool.add(norm)
    pool.discard("")
    return pool


def resolve_user(
    as_typed: str,
    scope_queryset: QuerySet,
    *,
    field: str,
    filter_index: int = 0,
    pinned_id: Optional[str] = None,
) -> UserResolution:
    """Resolve a user name/email/username against a pre-filtered User
    queryset (typically users who have registered any compound, or users
    who have created any assay).

    ``field`` and ``filter_index`` are tagged onto Clarify/Miss results
    so the UI knows which place on the selector to pin.
    """
    if pinned_id:
        user = scope_queryset.filter(pk=pinned_id).first()
        if user is None:
            return UserMiss(
                query=as_typed or "", suggestions=[],
                field=field, filter_index=filter_index,
            )
        return ResolvedUser(user=user, matched_via="pinned_id", query=as_typed or "")

    norm_query = normalize(as_typed or "")
    if not norm_query:
        return UserMiss(
            query=as_typed or "", suggestions=[],
            field=field, filter_index=filter_index,
        )

    users = list(
        scope_queryset.select_related("profile")
        .annotate(
            _n_compounds=Count("registered_compounds", distinct=True),
            _n_assays=Count("created_assays", distinct=True),
        )
    )

    exact = []
    partial = []  # substring of normalised query in any pool entry
    for user in users:
        pool = _user_pool(user)
        if norm_query in pool:
            exact.append(user)
            continue
        # Partial: user's query is contained in some pool entry (e.g. "alice"
        # in "alicejones"). Cheap and catches common prefix typing.
        if any(norm_query in entry for entry in pool):
            partial.append(user)

    def _to_candidate(user) -> UserCandidate:
        return _user_candidate(
            user,
            n_compounds=getattr(user, "_n_compounds", 0) or 0,
            n_assays=getattr(user, "_n_assays", 0) or 0,
        )

    if len(exact) == 1 and not partial:
        return ResolvedUser(user=exact[0], matched_via="exact", query=as_typed)
    if exact or partial:
        # Exact matches first; within each tier, more-prolific users first
        # (likely the intended "Alice" is the one with many compounds/assays).
        def _rank(user):
            if field == FIELD_ASSAYED_BY:
                return -(getattr(user, "_n_assays", 0) or 0)
            return -(getattr(user, "_n_compounds", 0) or 0)

        exact.sort(key=_rank)
        partial.sort(key=_rank)
        tiered = exact + partial
        tiered = tiered[: _CLARIFY_MAX_CANDIDATES]
        return UserClarify(
            query=as_typed,
            candidates=[_to_candidate(u) for u in tiered],
            field=field,
            filter_index=filter_index,
        )

    # Miss — surface the top few most-prolific users in scope as a hint.
    users.sort(
        key=lambda u: (
            -(getattr(u, "_n_assays" if field == FIELD_ASSAYED_BY else "_n_compounds", 0) or 0)
        )
    )
    suggestions = [_to_candidate(u) for u in users[:_MISS_SUGGESTION_COUNT]]
    return UserMiss(
        query=as_typed, suggestions=suggestions,
        field=field, filter_index=filter_index,
    )


# ---------------------------------------------------------------------------
# Scaffold / substructure resolution (slice 13)
# ---------------------------------------------------------------------------


def _scaffold_candidate(scaffold: Scaffold) -> ScaffoldCandidate:
    return ScaffoldCandidate(
        id=scaffold.name,
        name=scaffold.name,
        aliases=list(scaffold.aliases),
        smarts=scaffold.smarts,
    )


def resolve_scaffold(
    as_typed: str,
    *,
    scaffold_index: int = 0,
    pinned_id: Optional[str] = None,
) -> ScaffoldResolution:
    """Resolve a typed substructure / scaffold / functional-group name
    against the curated catalog. Three-tier match:

    1. Exact normalised match against name or alias → Resolved.
    2. Substring match on normalised pool entries → Clarify (multiple)
       or Resolved (single) — catches "pyrimidi" → "pyrimidine" and
       also surfaces the "pyrimidine vs pyrimidinone" ambiguity when
       the user types a substring that hits both.
    3. Miss with the closest matches by SequenceMatcher ratio — better
       than a blank dead end.
    """
    if pinned_id:
        scaffold = lookup_scaffold_by_id(pinned_id)
        if scaffold is None:
            return ScaffoldMiss(
                query=as_typed or "", suggestions=[],
                scaffold_index=scaffold_index, field=FIELD_SCAFFOLD_HINT,
            )
        return ResolvedScaffold(scaffold=scaffold, matched_via="pinned_id", query=as_typed or "")

    norm_query = normalize(as_typed or "")
    if not norm_query:
        return ScaffoldMiss(
            query=as_typed or "", suggestions=[],
            scaffold_index=scaffold_index, field=FIELD_SCAFFOLD_HINT,
        )

    # Tier 1: exact normalised match.
    exact = lookup_scaffold_by_typed(as_typed)
    if exact is not None:
        scaffold, matched_via = exact
        return ResolvedScaffold(scaffold=scaffold, matched_via=matched_via, query=as_typed)

    # Tier 2: substring on normalised pool entries.
    substring_hits = list(lookup_scaffolds_by_substring(as_typed))
    if len(substring_hits) == 1:
        scaffold, matched_via = substring_hits[0]
        return ResolvedScaffold(scaffold=scaffold, matched_via=matched_via, query=as_typed)
    if len(substring_hits) > 1:
        return ScaffoldClarify(
            query=as_typed,
            candidates=[_scaffold_candidate(s) for s, _ in substring_hits],
            scaffold_index=scaffold_index,
            field=FIELD_SCAFFOLD_HINT,
        )

    # Miss — rank whole catalog by SequenceMatcher ratio.
    scored: List[Tuple[float, Scaffold]] = []
    for s in SCAFFOLDS:
        # Use the best ratio across name + aliases.
        best = SequenceMatcher(None, norm_query, normalize(s.name)).ratio()
        for alias in s.aliases:
            best = max(best, SequenceMatcher(None, norm_query, normalize(alias)).ratio())
        scored.append((best, s))
    scored.sort(key=lambda row: row[0], reverse=True)
    suggestions = [_scaffold_candidate(s) for _, s in scored[:_MISS_SUGGESTION_COUNT]]
    return ScaffoldMiss(
        query=as_typed, suggestions=suggestions,
        scaffold_index=scaffold_index, field=FIELD_SCAFFOLD_HINT,
    )


def compound_registrars_qs() -> QuerySet:
    """All users who have registered at least one compound — the eligible
    pool for `registered_by` resolution."""
    return User.objects.filter(registered_compounds__isnull=False).distinct()


def assay_creators_qs() -> QuerySet:
    """All users who have created at least one assay — the eligible pool
    for per-filter `assayed_by` resolution."""
    return User.objects.filter(created_assays__isnull=False).distinct()
