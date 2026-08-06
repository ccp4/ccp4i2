"""API-level tests for the legacy-migration admin endpoints.

Asserts the wire contract (docs/LEGACY_MIGRATION_WIRE_CONTRACT.md): the
validate/import responses carry structure + plan, and errors use the
standardised {error: <machine_string>, reason: ...} envelope.

pytest-style (house convention for API unit tests). The admin views are
function-based with @permission_classes([IsAdminUser]); the autouse
``bypass_api_permissions`` fixture only patches viewsets, so real admin gating
still applies here and can be exercised with force_authenticate.
"""

import tempfile
import uuid
from pathlib import Path

import pytest
from django.contrib.auth import get_user_model
from rest_framework.test import APIRequestFactory, force_authenticate

from ccp4i2.api.admin_views import import_sqlite, validate_sqlite
from ccp4i2.db.models import Project
from ccp4i2.tests.db.test_import_sqlite import make_legacy_db, make_project_tree


@pytest.mark.django_db
class TestAdminMigrationApi:
    @pytest.fixture(autouse=True)
    def setup(self):
        User = get_user_model()
        self.admin = User.objects.create_superuser(
            username="admin", email="a@b.c", password="x")
        self.factory = APIRequestFactory()
        self.tmp = Path(tempfile.mkdtemp())
        self.legacy_root = self.tmp / "legacy"
        self.legacy_root.mkdir()
        self.dest_root = self.tmp / "CCP4X_PROJECTS"
        self.db_path = self.tmp / "database.sqlite"

    def _post(self, view, data, user="admin"):
        req = self.factory.post("/", data, format="multipart")
        if user:
            force_authenticate(req, user=self.admin)
        return view(req)

    def _make_db(self, name="gamma"):
        d = make_project_tree(self.legacy_root, name)
        pid = uuid.uuid4().hex
        make_legacy_db(self.db_path,
                       [{"id": pid, "name": name, "directory": d}],
                       [{"id": uuid.uuid4().hex, "number": "1", "project_id": pid}])
        return d

    def test_validate_returns_structure_and_plan(self):
        self._make_db()
        resp = self._post(validate_sqlite, {
            "db_path": str(self.db_path),
            "copy_files": "true",
            "dest_root": str(self.dest_root),
        })
        assert resp.status_code == 200
        for key in ("summary", "structure", "plan"):
            assert key in resp.data
        assert "blocking_issues" in resp.data["summary"]
        assert "plan_summary" in resp.data["summary"]
        assert resp.data["plan"][0]["mode"] == "copy"

    def test_bad_input_uses_standard_envelope(self):
        resp = self._post(validate_sqlite, {})
        assert resp.status_code == 400
        assert resp.data["error"] == "bad_input"
        assert "reason" in resp.data

    def test_db_not_found_envelope(self):
        resp = self._post(validate_sqlite, {"db_path": str(self.tmp / "nope.sqlite")})
        assert resp.status_code == 404
        assert resp.data["error"] == "db_not_found"

    def test_import_copy_mode_creates_project(self):
        self._make_db("delta")
        resp = self._post(import_sqlite, {
            "db_path": str(self.db_path),
            "copy_files": "true",
            "dest_root": str(self.dest_root),
        })
        assert resp.status_code == 200
        assert "plan" in resp.data
        assert Project.objects.filter(name="delta").exists()
        assert (self.dest_root / "delta" / "CCP4_JOBS").is_dir()

    def test_requires_admin(self):
        self._make_db()
        resp = self._post(validate_sqlite, {"db_path": str(self.db_path)}, user=None)
        assert resp.status_code in (401, 403)
