"""A KPI that is NaN or infinite must not take an endpoint down with it.

JSON has no NaN and no infinity, and DRF renders strictly (``STRICT_JSON``), so
one such value used to raise during rendering -- after the view had returned,
where its own try/except could not see it -- and the endpoint answered 500.
Three layers now stop that, and each is tested here on its own terms:

* the write gate, which keeps such a value out of the database;
* the read gate, which omits (never nulls) anything unservable from a KPI map;
* `SafeJSONRenderer`, which nulls any non-finite float anywhere in any payload.

The seeded values are infinities, not NaN, because SQLite refuses NaN in a
NOT NULL REAL column outright -- which is itself half of the original bug, and
is covered by the write-gate tests rather than by seeding.
"""

import json

import pytest
from rest_framework.test import APIClient

from ccp4i2.db.models import (
    Job, JobCharValue, JobFloatValue, JobValueKey, Project,
    ProjectGroup, ProjectGroupMembership,
)

API_PREFIX = "/api/ccp4i2"


@pytest.fixture
def project_with_a_bad_kpi(tmp_path):
    """One job carrying a good KPI, a bad one, and a string one."""
    directory = tmp_path / "toxd"
    directory.mkdir(parents=True, exist_ok=True)
    project = Project.objects.create(name="toxd", directory=str(directory))
    job = Job.objects.create(
        project=project, number="1", task_name="servalcat_pipe",
        title="Refinement", status=Job.Status.FINISHED,
    )
    for name in ("RFactor", "RFree", "spaceGroup"):
        JobValueKey.objects.get_or_create(name=name, defaults={"description": name})

    JobFloatValue.objects.create(
        job=job, key=JobValueKey.objects.get(name="RFactor"), value=0.21
    )
    JobFloatValue.objects.create(
        job=job, key=JobValueKey.objects.get(name="RFree"), value=float("inf")
    )
    JobCharValue.objects.create(
        job=job, key=JobValueKey.objects.get(name="spaceGroup"), value="P 21 21 21"
    )
    return {"project": project, "job": job}


class TestEndpointsThatEmbedKPIs:
    """Every reader of the KPI tables, not just the one that was reported."""

    @pytest.fixture(autouse=True)
    def setup(self, bypass_api_permissions):
        self.client = APIClient()

    def test_job_tree_serves_the_project(self, project_with_a_bad_kpi):
        project = project_with_a_bad_kpi["project"]
        response = self.client.get(f"{API_PREFIX}/projects/{project.id}/job_tree/")
        assert response.status_code == 200

        kpis = response.json()["job_tree"][0]["kpis"]
        assert kpis["float_values"] == {"RFactor": 0.21}
        assert kpis["char_values"] == {"spaceGroup": "P 21 21 21"}

    def test_job_tree_omits_the_bad_kpi_rather_than_nulling_it(
        self, project_with_a_bad_kpi
    ):
        # Load-bearing: the client tests `!== undefined` before .toFixed(), so
        # a null would render a confident "0.000" instead of a dash.
        project = project_with_a_bad_kpi["project"]
        response = self.client.get(f"{API_PREFIX}/projects/{project.id}/job_tree/")

        floats = response.json()["job_tree"][0]["kpis"]["float_values"]
        assert "RFree" not in floats

    def test_job_list_serves_the_job(self, project_with_a_bad_kpi):
        # JobSerializer embeds the same KPIs on every job list and detail.
        job = project_with_a_bad_kpi["job"]
        response = self.client.get(f"{API_PREFIX}/jobs/{job.id}/")
        assert response.status_code == 200

        payload = response.json()
        assert payload["float_values"] == {"RFactor": 0.21}
        assert "RFree" not in payload["float_values"]

    def test_campaign_member_projects_serves_the_row(self, project_with_a_bad_kpi):
        project = project_with_a_bad_kpi["project"]
        group = ProjectGroup.objects.create(
            name="Campaign", type=ProjectGroup.GroupType.FRAGMENT_SET
        )
        ProjectGroupMembership.objects.create(
            group=group, project=project,
            type=ProjectGroupMembership.MembershipType.MEMBER,
        )

        response = self.client.get(
            f"{API_PREFIX}/projectgroups/{group.id}/member_projects/"
        )
        assert response.status_code == 200

        payload = response.json()
        rows = payload["results"] if isinstance(payload, dict) else payload
        row = rows[0] if isinstance(rows, list) else rows
        kpis = row["kpis"] if "kpis" in row else row["results"][0]["kpis"]
        assert kpis["RFactor"] == 0.21
        assert "RFree" not in kpis

    def test_raw_float_values_endpoint_serves_null_not_a_500(
        self, project_with_a_bad_kpi
    ):
        # This endpoint dumps rows rather than building a KPI map, and its
        # serializer is also used for writes on the import path, so it is left
        # alone and the renderer backstop is what saves it. Null is right here:
        # the row exists, its value cannot be expressed.
        project = project_with_a_bad_kpi["project"]
        response = self.client.get(
            f"{API_PREFIX}/projects/{project.id}/job_float_values/"
        )
        assert response.status_code == 200

        values = {row["key"]: row["value"] for row in response.json()}
        assert values["RFactor"] == 0.21
        assert values["RFree"] is None


class TestSafeJSONRenderer:
    """The backstop, tested directly: it must not need a bad KPI to exist."""

    def test_renders_an_ordinary_payload_unchanged(self):
        from ccp4i2.api.renderers import SafeJSONRenderer

        data = {"a": 1, "b": [1.5, "x"], "c": None}
        assert json.loads(SafeJSONRenderer().render(data)) == data

    def test_nulls_a_nan_instead_of_raising(self):
        from ccp4i2.api.renderers import SafeJSONRenderer

        rendered = SafeJSONRenderer().render({"RFree": float("nan"), "n": 3})
        assert json.loads(rendered) == {"RFree": None, "n": 3}

    def test_emits_no_bare_nan_token(self):
        # The failure mode of "just set STRICT_JSON=False": valid to Python,
        # rejected by JSON.parse in the browser.
        from ccp4i2.api.renderers import SafeJSONRenderer

        rendered = SafeJSONRenderer().render({"x": float("-inf")})
        assert b"Infinity" not in rendered and b"NaN" not in rendered

    def test_a_value_error_from_elsewhere_still_propagates(self):
        from ccp4i2.api.renderers import SafeJSONRenderer

        class Unserialisable:
            pass

        with pytest.raises((TypeError, ValueError)):
            SafeJSONRenderer().render({"x": Unserialisable()})


class TestTheWriteGate:
    """Nothing should be able to put such a row there in the first place."""

    def test_sqlite_still_cannot_hold_a_nan(self, project_with_a_bad_kpi):
        # Not a rule we impose: the column is NOT NULL and SQLite maps NaN to
        # NULL. Pinned because it is why the same defect looks like a silently
        # truncated KPI list on the desktop and a 500 on the deployed servers.
        from django.db import IntegrityError, connection

        if connection.vendor != "sqlite":
            pytest.skip("SQLite-specific storage behaviour")

        job = project_with_a_bad_kpi["job"]
        key, _ = JobValueKey.objects.get_or_create(
            name="Rnan", defaults={"description": "Rnan"}
        )
        with pytest.raises(IntegrityError):
            JobFloatValue.objects.create(job=job, key=key, value=float("nan"))

    def test_gleaning_skips_a_non_finite_kpi(self, project_with_a_bad_kpi):
        from asgiref.sync import async_to_sync

        from ccp4i2.db.async_db_handler import AsyncDatabaseHandler

        job = project_with_a_bad_kpi["job"]
        handler = AsyncDatabaseHandler(project_uuid=job.project.uuid)

        async_to_sync(handler.register_job_float_value)(
            job_uuid=job.uuid, key="Rbad", value=float("inf")
        )
        assert not JobFloatValue.objects.filter(job=job, key_id="Rbad").exists()

    def test_gleaning_still_registers_a_good_kpi(self, project_with_a_bad_kpi):
        from asgiref.sync import async_to_sync

        from ccp4i2.db.async_db_handler import AsyncDatabaseHandler

        job = project_with_a_bad_kpi["job"]
        handler = AsyncDatabaseHandler(project_uuid=job.project.uuid)

        async_to_sync(handler.register_job_float_value)(
            job_uuid=job.uuid, key="Rgood", value=0.19
        )
        assert JobFloatValue.objects.get(job=job, key_id="Rgood").value == 0.19


class TestPruneCommand:
    """The rows already in deployed databases, which no gate can retro-fix."""

    def test_reports_without_changing_anything(self, project_with_a_bad_kpi):
        from io import StringIO

        from django.core.management import call_command

        out = StringIO()
        call_command("prune_nonfinite_kpis", stdout=out)

        assert "RFree" in out.getvalue()
        assert JobFloatValue.objects.filter(key_id="RFree").exists()

    def test_apply_deletes_only_the_bad_row(self, project_with_a_bad_kpi):
        from io import StringIO

        from django.core.management import call_command

        call_command("prune_nonfinite_kpis", "--apply", stdout=StringIO())

        assert not JobFloatValue.objects.filter(key_id="RFree").exists()
        assert JobFloatValue.objects.get(key_id="RFactor").value == 0.21
