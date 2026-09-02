"""Adopting a legacy install in place is off for the alpha; copying is not.

Agreed with Paul Bond and Stuart McNicholas, 2026-09-02. The distinction is the
MODE, not the endpoint: copying a legacy project to a new root writes somewhere
else and cannot damage what the user is still working with, while adopting in
place makes the old directories live under the new app while the old app may
still be using them. Read-only validation stays available either way —
inspecting a database cannot corrupt it.
"""

import json

import pytest
from django.test import override_settings
from rest_framework.test import APIClient


@pytest.fixture
def client(bypass_api_permissions):
    """An admin client: these endpoints sit behind IsAdminUser, which would
    otherwise 403 first and make this suite pass for the wrong reason."""
    from django.contrib.auth import get_user_model

    user_model = get_user_model()
    admin = user_model(username="migration-test-admin", is_staff=True,
                       is_superuser=True)
    admin.set_unusable_password()
    admin.save()

    api = APIClient()
    api.force_authenticate(user=admin)
    return api


def _post(client, url, payload=None):
    return client.post(
        url, data=json.dumps(payload or {}), content_type="application/json"
    )


@pytest.mark.parametrize(
    "url",
    [
        "/api/ccp4i2/admin/import-sqlite/",
        "/api/ccp4i2/admin/import-legacy/",
    ],
)
def test_in_place_mode_is_refused_by_default(client, url):
    resp = _post(client, url, {"db_path": "/nowhere/database.sqlite"})
    assert resp.status_code == 403
    body = resp.json()
    assert body["error"] == "in_place_migration_disabled"
    # The message has to say what to do instead, or it is just a wall.
    assert "import" in body["message"].lower()


@pytest.mark.parametrize(
    "url",
    [
        "/api/ccp4i2/admin/import-sqlite/",
        "/api/ccp4i2/admin/import-legacy/",
    ],
)
@override_settings(CCP4I2_ALLOW_INPLACE_MIGRATION=True)
def test_the_flag_re_enables_in_place(client, url):
    """Off by default, not deleted: the importer and its ~50 consistency checks
    stay intact for whoever makes this safe."""
    resp = _post(client, url, {"db_path": "/nowhere/database.sqlite"})
    assert resp.status_code != 403 or resp.json().get("error") != (
        "in_place_migration_disabled"
    )


def test_read_only_validation_still_works(client):
    """Gating this too would be cargo-culted caution: it cannot write."""
    resp = _post(client, "/api/ccp4i2/admin/validate-sqlite/",
                 {"db_path": "/nowhere/database.sqlite"})
    assert resp.status_code != 403


def test_copy_mode_is_still_allowed(client):
    """The mode the email explicitly kept: read the old installation, write the
    copies somewhere new. Gating this too would have blocked the one migration
    route testers ARE meant to use — an over-block the existing
    test_admin_migration_api suite caught.
    """
    resp = _post(client, "/api/ccp4i2/admin/import-sqlite/", {
        "db_path": "/nowhere/database.sqlite",
        "copy_files": "true",
        "dest_root": "/tmp/wherever",
    })
    assert resp.status_code != 403 or resp.json().get("error") != (
        "in_place_migration_disabled"
    )
