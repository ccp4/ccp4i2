"""Service contract guards for CCP4i2's REST API.

These tests enforce the v0 surface documented in
``docs/CCP4I2_SERVICE_CONTRACT.md`` and typed in
``packages/ccp4i2-auth/src/contracts/ccp4i2.ts``. They are
deliberately *narrow* — they assert that documented fields exist on
the relevant DRF serializers — so future refactors that silently
rename or remove a documented field fail this test rather than
silently breaking external consumers (Materia, third-party
integrators).

Tests in this file should never assert on undocumented fields. The
contract is the *promise*; everything else is internal and may
change. If you find yourself adding a contract test for a field
that's not in the docs, the right move is to either (a) add the
field to the contract document and TS types first, or (b) drop the
guard.

Implementation note: these tests inspect the serializer classes
directly via ``serializer.fields`` rather than driving HTTP requests
through the API client, so they don't need the heavyweight
project-zip-import test fixtures used by the rest of ``tests/api/``.
Fast, hermetic, no DB.
"""

import json

from ccp4i2.api import serializers


# Fields each contract type documents. Source of truth:
# packages/ccp4i2-auth/src/contracts/ccp4i2.ts.
PROJECT_LIST_ITEM_FIELDS = {
    "id",
    "uuid",
    "name",
    "creation_time",
    "last_access",
    "tags",
}

PROJECT_DETAIL_FIELDS = {
    "id",
    "uuid",
    "name",
    "description",
    "directory",
    "creation_time",
    "creation_user",
    "creation_host",
    "last_access",
    "last_job_number",
    "follow_from_job",
    "i1_project_name",
    "i1_project_directory",
    "tags",
}

JOB_FIELDS = {
    "id",
    "uuid",
    "project",
    "parent",
    "number",
    "title",
    "status",
    "evaluation",
    "comments",
    "creation_time",
    "finish_time",
    "task_name",
    "process_id",
    "float_values",
    "char_values",
}

FILE_FIELDS = {
    "id",
    "uuid",
    "name",
    "directory",
    "type",
    "sub_type",
    "content",
    "annotation",
    "job",
    "job_param_name",
    "path",
}


def _serializer_fields(serializer_class) -> set:
    """Return the bound field names of a DRF serializer class.
    Doesn't require an instance, an HTTP request, or a database.
    """
    return set(serializer_class().fields.keys())


# ---------------------------------------------------------------------------
# Serializer-shape guards
# ---------------------------------------------------------------------------


def test_project_list_item_serializer_publishes_documented_fields():
    """Failure means a refactor has changed the shape of
    ``GET /api/ccp4i2/projects/`` in a way that may break external
    consumers. Fix the regression OR update the contract document AND
    bump the package's major version.
    """
    fields = _serializer_fields(serializers.ProjectListSerializer)
    missing = PROJECT_LIST_ITEM_FIELDS - fields
    assert not missing, (
        f"ProjectListItem fields missing from ProjectListSerializer: "
        f"{missing}. See docs/CCP4I2_SERVICE_CONTRACT.md."
    )


def test_project_detail_serializer_publishes_documented_fields():
    fields = _serializer_fields(serializers.ProjectSerializer)
    missing = PROJECT_DETAIL_FIELDS - fields
    assert not missing, (
        f"Project fields missing from ProjectSerializer: {missing}. "
        f"See docs/CCP4I2_SERVICE_CONTRACT.md."
    )


def test_job_serializer_publishes_documented_fields():
    fields = _serializer_fields(serializers.JobSerializer)
    missing = JOB_FIELDS - fields
    assert not missing, (
        f"Job fields missing from JobSerializer: {missing}. "
        f"See docs/CCP4I2_SERVICE_CONTRACT.md."
    )


def test_file_serializer_publishes_documented_fields():
    fields = _serializer_fields(serializers.FileSerializer)
    missing = FILE_FIELDS - fields
    assert not missing, (
        f"File fields missing from FileSerializer: {missing}. "
        f"See docs/CCP4I2_SERVICE_CONTRACT.md."
    )


# ---------------------------------------------------------------------------
# Numeric enum stability — pattern-matched by external consumers
# ---------------------------------------------------------------------------


def test_job_status_numeric_values_are_stable():
    """Materia and third-party consumers pattern-match on these
    integer values (e.g., ``status === 3`` for Running). Changing
    them is a major-version-bump-required breaking change.
    """
    from ccp4i2.db.models import Job

    expected_status_codes = {
        "UNKNOWN": 0,
        "PENDING": 1,
        "QUEUED": 2,
        "RUNNING": 3,
        "INTERRUPTED": 4,
        "FAILED": 5,
        "FINISHED": 6,
        "RUNNING_REMOTELY": 7,
        "FILE_HOLDER": 8,
        "TO_DELETE": 9,
        "UNSATISFACTORY": 10,
    }
    for label, value in expected_status_codes.items():
        assert getattr(Job.Status, label).value == value, (
            f"Job.Status.{label} value changed from {value} to "
            f"{getattr(Job.Status, label).value}. Breaking change."
        )


def test_job_evaluation_numeric_values_are_stable():
    from ccp4i2.db.models import Job

    expected = {"UNKNOWN": 0, "BEST": 1, "GOOD": 2, "REJECTED": 3}
    for label, value in expected.items():
        assert getattr(Job.Evaluation, label).value == value, (
            f"Job.Evaluation.{label} value changed; breaking change."
        )


def test_file_directory_numeric_values_are_stable():
    from ccp4i2.db.models import File

    expected = {"JOB_DIR": 1, "IMPORT_DIR": 2}
    for label, value in expected.items():
        assert getattr(File.Directory, label).value == value, (
            f"File.Directory.{label} value changed; breaking change."
        )


# ---------------------------------------------------------------------------
# AuthErrorBody (any 401/403 response across the API)
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Action endpoints (existence + URL-name stability)
# ---------------------------------------------------------------------------


def test_report_xml_action_is_registered_on_jobviewset():
    """``GET /jobs/{id}/report_xml/`` is contracted; the action must
    remain a registered @action on JobViewSet (renaming or removing it
    is a breaking change). The body shape is XML-flavoured and not
    typed in TS, so this is presence-only.
    """
    from ccp4i2.api.JobViewSet import JobViewSet

    method = getattr(JobViewSet, "report_xml", None)
    assert method is not None, "JobViewSet.report_xml action removed"
    bind = getattr(method, "mapping", None) or getattr(method, "bind_to_methods", None)
    # DRF stores the http-method binding on @action-decorated methods;
    # we don't assert the exact attribute name (varies across DRF versions),
    # only that *something* is bound — the action exists, not just any callable.
    assert bind, "JobViewSet.report_xml is not a registered DRF @action"


def test_resolve_fileuse_parser_returns_documented_shape():
    """The fileUse DSL is contracted; the parser's return shape (the
    four keys ``task_name`` / ``jobIndex`` / ``jobParamName`` /
    ``paramIndex``) is what the resolver consumes and what consumers
    pattern-match on. Pure-function, no DB.
    """
    from ccp4i2.lib.utils.files.resolve_fileuse import parse_fileuse

    expected_keys = {"task_name", "jobIndex", "jobParamName", "paramIndex"}

    # Each of the four contracted syntactic forms.
    cases = [
        # task_name[jobIndex].jobParamName[paramIndex]
        ("refmac[-2].HKLOUT[0]", {"task_name": "refmac", "jobIndex": -2, "jobParamName": "HKLOUT", "paramIndex": 0}),
        # [jobIndex].jobParamName[paramIndex]
        ("[-1].XYZOUT[0]", {"task_name": None, "jobIndex": -1, "jobParamName": "XYZOUT", "paramIndex": 0}),
        # task_name[jobIndex].jobParamName
        ("prosmart_refmac[-1].XYZOUT", {"task_name": "prosmart_refmac", "jobIndex": -1, "jobParamName": "XYZOUT", "paramIndex": 0}),
        # [jobIndex].jobParamName
        ("[-1].HKLOUT", {"task_name": None, "jobIndex": -1, "jobParamName": "HKLOUT", "paramIndex": 0}),
    ]

    for fileuse, expected in cases:
        parsed = parse_fileuse(fileuse)
        assert set(parsed.keys()) == expected_keys, (
            f"parse_fileuse({fileuse!r}) returned keys {set(parsed.keys())}; "
            f"contract requires {expected_keys}"
        )
        for key, value in expected.items():
            assert parsed[key] == value, (
                f"parse_fileuse({fileuse!r}) -> {key}={parsed[key]!r}; expected {value!r}"
            )


def test_auth_error_response_shape_is_stable():
    """The ``{success: false, error: string}`` shape is referenced
    across consumers (api-fetch's AUTH_ERROR_EVENT, the Materia
    binding). Encoded into ``BaseAuthMiddleware._error_response`` —
    this guard fails if the shape drifts.
    """
    from ccp4i2_auth.middleware.base import BaseAuthMiddleware

    response = BaseAuthMiddleware._error_response("nope", status=401)
    body = json.loads(response.content.decode())

    assert body["success"] is False
    assert isinstance(body["error"], str)
    assert response.status_code == 401

    # Same shape with 403.
    response = BaseAuthMiddleware._error_response(
        "denied", status=403, prefix="Access denied"
    )
    body = json.loads(response.content.decode())
    assert body["success"] is False
    assert isinstance(body["error"], str)
    assert response.status_code == 403
