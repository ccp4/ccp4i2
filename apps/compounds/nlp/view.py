"""DRF endpoint: POST /api/compounds/nlp/query/ (§7 view, §9 clarify, §10 redirect).

After the pivot (2026-04-23), the endpoint no longer renders its own
table. Successful selections return a ``redirect_url`` pointing at the
existing ``/assays/aggregate`` page with the matching compound IDs
pre-populated in the query string. The aggregation page owns the display
(cards, compact, pivot, …) and the chemist lands there ready to eyeball.

Two request shapes, unchanged in spirit:

- Fresh:        {"prompt": "<english>"}           → parse_prompt → execute
- Continuation: {"selector": {<CompoundSelector>}} → skip LLM; execute with
                                                     the client-supplied
                                                     pinning IDs (clarify)

Response shapes discriminated by ``status``:
- selection    → redirect_url + n_matched + compound_formatted_ids + target/protocol names + scope_sentence
- clarify      → target or protocol picker; partial_selector round-trips back
- miss         → suggestions
- not_a_query  → LLM's reason
- error        → structured error payloads (disabled / bad_request / scope /
                 spec / rate_limited / parse)

Feature-gated via ``COMPOUNDS_NLP_ENABLED``. Daily call cap per user via
Django's default cache (§11, decision 17).
"""

from __future__ import annotations

import dataclasses
import datetime
import logging
import os
from typing import Any, List, Optional, Tuple
from urllib.parse import quote

from django.core.cache import cache
from rest_framework import status as http_status
from rest_framework.decorators import api_view, permission_classes
from rest_framework.permissions import IsAuthenticated
from rest_framework.request import Request
from rest_framework.response import Response

from .executor import execute, execute_assay_query
from .llm import parse_prompt
from .spec import (
    AssaySelection,
    AssaySelector,
    CompoundSelection,
    CompoundSelector,
    DateRange,
    MeasurementFilter,
    NotAQuery,
    ParseError,
    ProtocolClarify,
    ProtocolMiss,
    ScaffoldClarify,
    CompoundMiss,
    MetricMiss,
    ScaffoldMiss,
    ScopeError,
    SpecError,
    TargetClarify,
    TargetMiss,
    Threshold,
    UserClarify,
    UserMiss,
)

logger = logging.getLogger(__name__)

ENV_FLAG = "COMPOUNDS_NLP_ENABLED"

# §11 runaway-loop cap (decision 17). Per-process LocMemCache — acceptable
# because the cap is anti-runaway, not cost-protection. See §18.3 slice-6.
DAILY_CAP = 500
_CAP_TTL_SECONDS = 26 * 60 * 60

# Default destination for successful selections. Per §10 pivot note: the
# chemist workflow is "find compounds → view spider-plot cards".
REDIRECT_BASE_PATH = "/assays/aggregate"
REDIRECT_DEFAULT_FORMAT = "cards"

# Assay-selection landing page (sibling of the aggregation page).
ASSAYS_REDIRECT_BASE_PATH = "/assays"


# ---------------------------------------------------------------------------
# Rate-limit helpers
# ---------------------------------------------------------------------------


def _cap_key(user_pk: Any) -> str:
    today = datetime.date.today().isoformat()
    return f"nlp_cap:{user_pk}:{today}"


def _check_and_increment_cap(user_pk: Any) -> bool:
    key = _cap_key(user_pk)
    count = cache.get(key, 0)
    if count >= DAILY_CAP:
        return False
    cache.set(key, count + 1, timeout=_CAP_TTL_SECONDS)
    return True


# ---------------------------------------------------------------------------
# Selector deserialization (continuation path)
# ---------------------------------------------------------------------------


def _threshold_from_dict(data: Any) -> Optional[Threshold]:
    if data is None:
        return None
    if not isinstance(data, dict):
        raise ValueError("threshold must be an object")
    op = data.get("op")
    value = data.get("value")
    if not isinstance(op, str) or not isinstance(value, (int, float)):
        raise ValueError("threshold.op must be a string and threshold.value a number")
    return Threshold(op=op, value=float(value), unit=data.get("unit"))


def _date_range_from_dict(data: Any) -> Optional[DateRange]:
    if data is None:
        return None
    if not isinstance(data, dict):
        raise ValueError("date range must be an object")
    after = data.get("after")
    before = data.get("before")
    if after is None and before is None:
        return None
    if after is not None and not isinstance(after, str):
        raise ValueError("date range 'after' must be a string")
    if before is not None and not isinstance(before, str):
        raise ValueError("date range 'before' must be a string")
    return DateRange(after=after, before=before)


def _filter_from_dict(data: Any) -> MeasurementFilter:
    if not isinstance(data, dict):
        raise ValueError("measurement_filter must be an object")
    return MeasurementFilter(
        protocol_hint=data.get("protocol_hint"),
        metric=data.get("metric"),
        threshold=_threshold_from_dict(data.get("threshold")),
        assay_date_range=_date_range_from_dict(data.get("assay_date_range")),
        assayed_by_as_typed=data.get("assayed_by_as_typed"),
        protocol_id=data.get("protocol_id"),
        assayed_by_id=data.get("assayed_by_id"),
    )


def _selector_from_dict(data: Any) -> CompoundSelector:
    if not isinstance(data, dict):
        raise ValueError("selector must be an object")
    raw_filters = data.get("measurement_filters") or []
    if not isinstance(raw_filters, list):
        raise ValueError("measurement_filters must be a list")
    return CompoundSelector(
        registration_target_as_typed=data.get("registration_target_as_typed"),
        assay_target_as_typed=data.get("assay_target_as_typed"),
        measurement_filters=[_filter_from_dict(f) for f in raw_filters],
        registered_date_range=_date_range_from_dict(data.get("registered_date_range")),
        registered_by_as_typed=data.get("registered_by_as_typed"),
        scaffold_hints=list(data.get("scaffold_hints") or []),
        compound_refs_as_typed=list(data.get("compound_refs_as_typed") or []),
        registration_target_id=data.get("registration_target_id"),
        assay_target_id=data.get("assay_target_id"),
        registered_by_id=data.get("registered_by_id"),
        registered_by_supplier_id=data.get("registered_by_supplier_id"),
        scaffold_ids=list(data.get("scaffold_ids") or []),
    )


# ---------------------------------------------------------------------------
# Redirect URL construction
# ---------------------------------------------------------------------------


def _build_assay_redirect_url(selection: AssaySelection) -> str:
    """Construct the /assays?ids=... URL from an AssaySelection. Comma-
    separated UUID list is the minimum the /assays page needs to render
    exactly the filtered set."""
    if not selection.assay_ids:
        # No matches — send the user to the plain /assays listing rather
        # than /assays?ids=<empty> (which would show nothing).
        return ASSAYS_REDIRECT_BASE_PATH
    return (
        f"{ASSAYS_REDIRECT_BASE_PATH}?ids="
        + quote(",".join(selection.assay_ids), safe=",")
    )


def _selector_for_echo_from_dict(data: dict):
    """Client continuation body carries a selector; pick which type based
    on whether it's shaped like an AssaySelector or CompoundSelector."""
    if not isinstance(data, dict):
        raise ValueError("selector must be an object")
    # Heuristic: an assay-shaped selector has `target_as_typed` (single
    # field) rather than the two-target CompoundSelector fields. Use
    # the presence of `measurement_filters` as the compound-shape tell.
    if "measurement_filters" in data:
        return _selector_from_dict(data)
    return _assay_selector_from_dict(data)


def _assay_selector_from_dict(data: Any) -> AssaySelector:
    if not isinstance(data, dict):
        raise ValueError("selector must be an object")
    return AssaySelector(
        target_as_typed=data.get("target_as_typed"),
        protocol_hint=data.get("protocol_hint"),
        date_range=_date_range_from_dict(data.get("date_range")),
        created_by_as_typed=data.get("created_by_as_typed"),
        scaffold_hints=list(data.get("scaffold_hints") or []),
        target_id=data.get("target_id"),
        protocol_id=data.get("protocol_id"),
        created_by_id=data.get("created_by_id"),
        scaffold_ids=list(data.get("scaffold_ids") or []),
    )


def _build_redirect_url(selection: CompoundSelection) -> str:
    """Construct the /assays/aggregate URL from a selection.

    Only emits query params that differ from the aggregation page's
    defaults (per the URL contract documented alongside that page).
    `format=cards` IS emitted even though it's a non-default override
    from `compact`, because the NLP pivot's default destination is the
    project-card / spider-plot view (§10).
    """
    params: List[Tuple[str, str]] = []

    if selection.target_names:
        params.append(("targets", ",".join(selection.target_names)))
    if selection.protocol_names:
        params.append(("protocols", ",".join(selection.protocol_names)))
    if selection.compound_formatted_ids:
        params.append(("compound", ",".join(selection.compound_formatted_ids)))
    params.append(("format", REDIRECT_DEFAULT_FORMAT))

    encoded = "&".join(f"{k}={quote(v, safe=',')}" for k, v in params)
    return f"{REDIRECT_BASE_PATH}?{encoded}"


# ---------------------------------------------------------------------------
# Result serialization
# ---------------------------------------------------------------------------


def _serialize(
    result: Any,
    selector: Optional[Any] = None,
) -> Tuple[dict, int]:
    """Map an execution / parse result to (response_body, http_status).

    ``selector`` is the CompoundSelector OR AssaySelector that produced
    this result — echoed back as ``partial_selector`` on clarify responses
    so the client can re-POST with a pinning ID without re-calling the
    LLM (§9, decision 5).
    """
    if isinstance(result, CompoundSelection):
        body = {
            "status": "selection",
            "redirect_url": _build_redirect_url(result),
            **dataclasses.asdict(result),
        }
        return body, http_status.HTTP_200_OK

    if isinstance(result, AssaySelection):
        body = {
            "status": "assay_selection",
            "redirect_url": _build_assay_redirect_url(result),
            **dataclasses.asdict(result),
        }
        return body, http_status.HTTP_200_OK

    if isinstance(result, (TargetClarify, ProtocolClarify, UserClarify, ScaffoldClarify)):
        body = {"status": "clarify", **dataclasses.asdict(result)}
        if selector is not None:
            body["partial_selector"] = dataclasses.asdict(selector)
        return body, http_status.HTTP_200_OK

    if isinstance(
        result,
        (TargetMiss, ProtocolMiss, UserMiss, ScaffoldMiss, CompoundMiss, MetricMiss),
    ):
        body = {"status": "miss", **dataclasses.asdict(result)}
        if isinstance(result, CompoundMiss):
            # Compound IDs are deterministic — no fuzzy hits to surface.
            # Add an empty list so the frontend's miss renderer is uniform.
            body["suggestions"] = []
        if isinstance(result, MetricMiss):
            # MetricMiss carries `available_metrics` (the KPIs that ARE
            # in scope) rather than fuzzy candidate suggestions; provide
            # an empty `suggestions` so the frontend's uniform miss
            # renderer doesn't trip.
            body.setdefault("suggestions", [])
        return body, http_status.HTTP_200_OK

    if isinstance(result, NotAQuery):
        return {
            "status": "not_a_query",
            "reason": result.reason,
        }, http_status.HTTP_200_OK

    if isinstance(result, ScopeError):
        return {
            "status": "error",
            "kind": "scope",
            "message": result.message,
        }, http_status.HTTP_400_BAD_REQUEST

    if isinstance(result, SpecError):
        return {
            "status": "error",
            "kind": "spec",
            "field": result.field,
            "message": result.message,
        }, http_status.HTTP_400_BAD_REQUEST

    if isinstance(result, ParseError):
        return {
            "status": "error",
            "kind": "parse",
            "message": result.message,
        }, http_status.HTTP_502_BAD_GATEWAY

    logger.error("nlp_query: unknown result type %r", type(result))
    return {
        "status": "error",
        "kind": "internal",
        "message": "Unknown result type.",
    }, http_status.HTTP_500_INTERNAL_SERVER_ERROR


# ---------------------------------------------------------------------------
# View
# ---------------------------------------------------------------------------


def _feature_enabled() -> bool:
    return os.environ.get(ENV_FLAG, "").lower() in ("true", "1", "yes")


@api_view(["POST"])
@permission_classes([IsAuthenticated])
def nlp_query(request: Request) -> Response:
    """Translate an English prompt (or a pinned CompoundSelector) into a
    redirect URL for the aggregation page."""
    if not _feature_enabled():
        return Response(
            {"status": "error", "kind": "disabled",
             "message": "NLP query feature is not enabled on this instance."},
            status=http_status.HTTP_404_NOT_FOUND,
        )

    if not _check_and_increment_cap(request.user.pk):
        return Response(
            {"status": "error", "kind": "rate_limited",
             "message": f"Daily limit of {DAILY_CAP} NLP queries reached."},
            status=http_status.HTTP_429_TOO_MANY_REQUESTS,
        )

    body = request.data or {}
    prompt = body.get("prompt")
    selector_dict = body.get("selector")

    if prompt is not None and selector_dict is not None:
        return Response(
            {"status": "error", "kind": "bad_request",
             "message": "Provide either 'prompt' or 'selector', not both."},
            status=http_status.HTTP_400_BAD_REQUEST,
        )

    selector_for_echo: Optional[Any] = None
    if prompt is not None:
        if not isinstance(prompt, str):
            return Response(
                {"status": "error", "kind": "bad_request",
                 "message": "'prompt' must be a string."},
                status=http_status.HTTP_400_BAD_REQUEST,
            )
        parse_result = parse_prompt(prompt)
        if isinstance(parse_result, CompoundSelector):
            selector_for_echo = parse_result
            result: Any = execute(parse_result)
        elif isinstance(parse_result, AssaySelector):
            selector_for_echo = parse_result
            result = execute_assay_query(parse_result)
        else:
            result = parse_result  # NotAQuery or ParseError
    elif selector_dict is not None:
        try:
            selector_for_echo = _selector_for_echo_from_dict(selector_dict)
        except ValueError as e:
            return Response(
                {"status": "error", "kind": "bad_request",
                 "message": f"Invalid selector: {e}"},
                status=http_status.HTTP_400_BAD_REQUEST,
            )
        if isinstance(selector_for_echo, AssaySelector):
            result = execute_assay_query(selector_for_echo)
        else:
            result = execute(selector_for_echo)
    else:
        return Response(
            {"status": "error", "kind": "bad_request",
             "message": "Body must contain 'prompt' or 'selector'."},
            status=http_status.HTTP_400_BAD_REQUEST,
        )

    body_out, status_code = _serialize(result, selector=selector_for_echo)
    return Response(body_out, status=status_code)


# ---------------------------------------------------------------------------
# Scaffold-extension endpoint (slice 17) — runtime extension of the
# substructure catalog. Lets chemists add fragments the seed Python
# module doesn't cover, scoped to a project (Target) or shared.
# ---------------------------------------------------------------------------


@api_view(["POST"])
@permission_classes([IsAuthenticated])
def nlp_scaffold_extend(request: Request) -> Response:
    """Create a ScaffoldExtension row.

    Body:
        name (str, required): canonical name chemists will type.
        smarts (str, required): SMARTS pattern.
        aliases (list[str], optional): plurals / abbreviations.
        target_id (str|null, optional): UUID of the project Target. Null
            (or omitted) = shared catalog (visible to every project).
        notes (str, optional): human gloss.
        source_prompt (str, optional): the original NLP prompt that
            triggered the addition, for audit.

    Same feature flag as /nlp/query (so a deployment with NLP disabled
    can't grow the catalog through this surface either). Authenticated
    users only. Validates SMARTS via RDKit before persisting.
    """
    if not _feature_enabled():
        return Response(
            {"status": "error", "kind": "disabled",
             "message": "NLP query feature is not enabled on this instance."},
            status=http_status.HTTP_404_NOT_FOUND,
        )

    body = request.data or {}
    name = (body.get("name") or "").strip()
    smarts = (body.get("smarts") or "").strip()
    if not name or not smarts:
        return Response(
            {"status": "error", "kind": "bad_request",
             "message": "Both 'name' and 'smarts' are required."},
            status=http_status.HTTP_400_BAD_REQUEST,
        )

    aliases = body.get("aliases") or []
    if not isinstance(aliases, list) or not all(isinstance(a, str) for a in aliases):
        return Response(
            {"status": "error", "kind": "bad_request",
             "message": "'aliases' must be a list of strings."},
            status=http_status.HTTP_400_BAD_REQUEST,
        )

    # Validate SMARTS server-side via RDKit before persisting — a
    # parse failure here is much cheaper than a runtime crash later
    # when the resolver tries to use it.
    try:
        from rdkit import Chem
    except ImportError:
        return Response(
            {"status": "error", "kind": "internal",
             "message": "RDKit is not available on the server."},
            status=http_status.HTTP_500_INTERNAL_SERVER_ERROR,
        )
    if Chem.MolFromSmarts(smarts) is None:
        return Response(
            {"status": "error", "kind": "bad_request", "field": "smarts",
             "message": f"SMARTS did not parse: {smarts!r}"},
            status=http_status.HTTP_400_BAD_REQUEST,
        )

    # Resolve target (if provided).
    target_id = body.get("target_id")
    target_obj = None
    if target_id:
        from compounds.registry.models import Target
        target_obj = Target.objects.filter(pk=target_id).first()
        if target_obj is None:
            return Response(
                {"status": "error", "kind": "bad_request", "field": "target_id",
                 "message": f"Target {target_id!r} not found."},
                status=http_status.HTTP_400_BAD_REQUEST,
            )

    from compounds.registry.models import ScaffoldExtension
    # Reject same-name collision in the same scope (the unique_together
    # would do it via IntegrityError, but a clean 400 is friendlier).
    if ScaffoldExtension.objects.filter(name=name, target=target_obj).exists():
        scope = target_obj.name if target_obj else "shared"
        return Response(
            {"status": "error", "kind": "conflict", "field": "name",
             "message": f"A scaffold extension named {name!r} already exists in scope {scope!r}."},
            status=http_status.HTTP_409_CONFLICT,
        )

    ext = ScaffoldExtension.objects.create(
        name=name,
        smarts=smarts,
        aliases=aliases,
        target=target_obj,
        notes=(body.get("notes") or "").strip(),
        created_by=request.user if request.user.is_authenticated else None,
        source_prompt=(body.get("source_prompt") or "").strip(),
        llm_generated=bool(body.get("llm_generated", False)),
    )
    return Response(
        {
            "status": "created",
            "id": str(ext.pk),
            "name": ext.name,
            "smarts": ext.smarts,
            "aliases": list(ext.aliases or []),
            "target_id": str(target_obj.id) if target_obj else None,
            "target_name": target_obj.name if target_obj else None,
            "notes": ext.notes,
        },
        status=http_status.HTTP_201_CREATED,
    )
