"""DRF endpoint: POST /api/compounds/nlp/query/ (§7 view, §9 clarify loop).

Composes slices 4 (executor) and 5 (LLM parse) into a single HTTP surface:

- Fresh query: body = {"prompt": "<english>"} — runs parse_prompt → execute.
- Continuation: body = {"spec": <QuerySpec dict>} — skips the LLM entirely
  (decision 5) and runs execute() on the spec the client round-tripped,
  optionally with `registration_target_id` / `assay_target_id` /
  `protocol_id` filled in from a clarify picker.

Feature-gated via `COMPOUNDS_NLP_ENABLED` env var. Daily-call cap per user
via Django's default cache (§11, decision 17 — best-effort runaway guard,
not a cost gate; per-process LocMemCache is fine for that framing).
"""

from __future__ import annotations

import dataclasses
import datetime
import logging
import os
from typing import Any, Optional, Tuple

from django.core.cache import cache
from rest_framework import status as http_status
from rest_framework.decorators import api_view, permission_classes
from rest_framework.permissions import IsAuthenticated
from rest_framework.request import Request
from rest_framework.response import Response

from .executor import execute
from .llm import parse_prompt
from .spec import (
    NotAQuery,
    ParseError,
    ProtocolClarify,
    ProtocolMiss,
    QuerySpec,
    ScopeError,
    SpecError,
    TablePayload,
    TargetClarify,
    TargetMiss,
    Threshold,
)

logger = logging.getLogger(__name__)

ENV_FLAG = "COMPOUNDS_NLP_ENABLED"

# §11 runaway-loop cap: decision 17 calls this out as "catch retry storms,
# not cost protection". Per-process counter via LocMemCache is acceptable at
# that framing — a user who happens to hit N workers could reach N × 500,
# which is still orders of magnitude above real human workload.
DAILY_CAP = 500
_CAP_TTL_SECONDS = 26 * 60 * 60  # just over a day to survive DST edges


# ---------------------------------------------------------------------------
# Rate-limit helpers
# ---------------------------------------------------------------------------


def _cap_key(user_pk: Any) -> str:
    today = datetime.date.today().isoformat()
    return f"nlp_cap:{user_pk}:{today}"


def _check_and_increment_cap(user_pk: Any) -> bool:
    """True if the call is allowed (counter incremented); False at cap."""
    key = _cap_key(user_pk)
    count = cache.get(key, 0)
    if count >= DAILY_CAP:
        return False
    # get→set race is tolerable for a runaway-loop guard; see module docstring.
    cache.set(key, count + 1, timeout=_CAP_TTL_SECONDS)
    return True


# ---------------------------------------------------------------------------
# Spec deserialization (continuation path)
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


def _spec_from_dict(data: Any) -> QuerySpec:
    if not isinstance(data, dict):
        raise ValueError("spec must be an object")
    return QuerySpec(
        registration_target_as_typed=data.get("registration_target_as_typed"),
        assay_target_as_typed=data.get("assay_target_as_typed"),
        protocol_hint=data.get("protocol_hint"),
        metric=data.get("metric"),
        threshold=_threshold_from_dict(data.get("threshold")),
        columns=data.get("columns"),
        columns_explicit=data.get("columns_explicit"),
        registration_target_id=data.get("registration_target_id"),
        assay_target_id=data.get("assay_target_id"),
        protocol_id=data.get("protocol_id"),
    )


# ---------------------------------------------------------------------------
# Result serialization
# ---------------------------------------------------------------------------


def _serialize(result: Any) -> Tuple[dict, int]:
    """Map an execution / parse result to (response_body, http_status)."""
    if isinstance(result, TablePayload):
        body = {"status": "table", **dataclasses.asdict(result)}
        return body, http_status.HTTP_200_OK

    if isinstance(result, (TargetClarify, ProtocolClarify)):
        # `field` is on the dataclass already — orchestrator sets it for
        # TargetClarify, ProtocolClarify has FIELD_PROTOCOL_HINT by default.
        body = {"status": "clarify", **dataclasses.asdict(result)}
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

    # Should not reach here; the union is closed. Log loudly and fail safe.
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
    """Translate an English prompt (or a pinned QuerySpec) into a table."""
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
    spec_dict = body.get("spec")

    if prompt is not None and spec_dict is not None:
        return Response(
            {"status": "error", "kind": "bad_request",
             "message": "Provide either 'prompt' or 'spec', not both."},
            status=http_status.HTTP_400_BAD_REQUEST,
        )

    if prompt is not None:
        if not isinstance(prompt, str):
            return Response(
                {"status": "error", "kind": "bad_request",
                 "message": "'prompt' must be a string."},
                status=http_status.HTTP_400_BAD_REQUEST,
            )
        parse_result = parse_prompt(prompt)
        if isinstance(parse_result, QuerySpec):
            result: Any = execute(parse_result)
        else:
            result = parse_result  # NotAQuery or ParseError
    elif spec_dict is not None:
        try:
            spec = _spec_from_dict(spec_dict)
        except ValueError as e:
            return Response(
                {"status": "error", "kind": "bad_request",
                 "message": f"Invalid spec: {e}"},
                status=http_status.HTTP_400_BAD_REQUEST,
            )
        result = execute(spec)
    else:
        return Response(
            {"status": "error", "kind": "bad_request",
             "message": "Body must contain 'prompt' or 'spec'."},
            status=http_status.HTTP_400_BAD_REQUEST,
        )

    body_out, status_code = _serialize(result)
    return Response(body_out, status=status_code)
