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

from .executor import execute
from .llm import parse_prompt
from .spec import (
    CompoundSelection,
    CompoundSelector,
    MeasurementFilter,
    NotAQuery,
    ParseError,
    ProtocolClarify,
    ProtocolMiss,
    ScopeError,
    SpecError,
    TargetClarify,
    TargetMiss,
    Threshold,
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


def _filter_from_dict(data: Any) -> MeasurementFilter:
    if not isinstance(data, dict):
        raise ValueError("measurement_filter must be an object")
    return MeasurementFilter(
        protocol_hint=data.get("protocol_hint"),
        metric=data.get("metric"),
        threshold=_threshold_from_dict(data.get("threshold")),
        protocol_id=data.get("protocol_id"),
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
        registration_target_id=data.get("registration_target_id"),
        assay_target_id=data.get("assay_target_id"),
    )


# ---------------------------------------------------------------------------
# Redirect URL construction
# ---------------------------------------------------------------------------


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
    selector: Optional[CompoundSelector] = None,
) -> Tuple[dict, int]:
    """Map an execution / parse result to (response_body, http_status).

    ``selector`` is the CompoundSelector that produced this result — echoed
    back as ``partial_selector`` on clarify responses so the client can
    re-POST with a pinning ID without re-calling the LLM (§9, decision 5).
    """
    if isinstance(result, CompoundSelection):
        body = {
            "status": "selection",
            "redirect_url": _build_redirect_url(result),
            **dataclasses.asdict(result),
        }
        return body, http_status.HTTP_200_OK

    if isinstance(result, (TargetClarify, ProtocolClarify)):
        body = {"status": "clarify", **dataclasses.asdict(result)}
        if selector is not None:
            body["partial_selector"] = dataclasses.asdict(selector)
        return body, http_status.HTTP_200_OK

    if isinstance(result, (TargetMiss, ProtocolMiss)):
        body = {"status": "miss", **dataclasses.asdict(result)}
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

    selector_for_echo: Optional[CompoundSelector] = None
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
        else:
            result = parse_result  # NotAQuery or ParseError
    elif selector_dict is not None:
        try:
            selector_for_echo = _selector_from_dict(selector_dict)
        except ValueError as e:
            return Response(
                {"status": "error", "kind": "bad_request",
                 "message": f"Invalid selector: {e}"},
                status=http_status.HTTP_400_BAD_REQUEST,
            )
        result = execute(selector_for_echo)
    else:
        return Response(
            {"status": "error", "kind": "bad_request",
             "message": "Body must contain 'prompt' or 'selector'."},
            status=http_status.HTTP_400_BAD_REQUEST,
        )

    body_out, status_code = _serialize(result, selector=selector_for_echo)
    return Response(body_out, status=status_code)
