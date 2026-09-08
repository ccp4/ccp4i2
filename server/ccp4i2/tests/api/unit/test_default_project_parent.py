"""The New Project dialog must be told where a project will actually go.

The server resolves an unspecified directory to the configured projects
directory — always, regardless of where any previous project landed. The
dialog cannot compute that itself without a second resolver for the same
question, so it asks. It can also change or reset that default (desktop
only) through the paired "set/" endpoint.
"""

import json

import pytest
from rest_framework.test import APIClient

from ccp4i2.db.models import Project


@pytest.fixture
def client(bypass_api_permissions):
    return APIClient()


URL = "/api/ccp4i2/config/default-project-parent/"
SET_URL = "/api/ccp4i2/config/default-project-parent/set/"


def test_proposes_the_configured_projects_directory(client, settings):
    resp = client.get(URL)
    assert resp.status_code == 200
    data = resp.json()["data"]
    assert data["directory"] == str(settings.CCP4I2_PROJECTS_DIR)
    assert data["configured"] == str(settings.CCP4I2_PROJECTS_DIR)


def test_ignores_where_the_most_recent_project_landed(client, tmp_path):
    """A one-off project made somewhere unusual must not nudge the default
    proposed to the next New Project dialog — that is what stranded the
    configured default in a hidden home folder nobody could find their way
    back to."""
    elsewhere = tmp_path / "somewhere_else"
    (elsewhere / "proj").mkdir(parents=True)
    Project.objects.create(name="proj", directory=str(elsewhere / "proj"))

    data = client.get(URL).json()["data"]
    assert data["directory"] != str(elsewhere)


def test_configured_root_is_always_reported(client, settings):
    assert client.get(URL).json()["data"]["configured"] == str(
        settings.CCP4I2_PROJECTS_DIR
    )


def test_editable_only_on_desktop(client, monkeypatch):
    monkeypatch.setenv("CCP4I2_LOCAL_SESSION_TOKEN", "desktop")
    assert client.get(URL).json()["data"]["editable"] is True

    monkeypatch.delenv("CCP4I2_LOCAL_SESSION_TOKEN", raising=False)
    assert client.get(URL).json()["data"]["editable"] is False


def test_set_and_reset_roundtrip_on_desktop(client, monkeypatch, tmp_path):
    """The New Project page's "make this the default" checkbox and the
    Preferences page's "Reset to default" both go through this endpoint, and
    the very next GET — from any window — must reflect it immediately, with
    no server restart needed."""
    monkeypatch.setenv("CCP4I2_HOME", str(tmp_path))
    monkeypatch.setenv("CCP4I2_LOCAL_SESSION_TOKEN", "desktop")
    monkeypatch.delenv("CCP4I2_PROJECTS_DIR", raising=False)

    chosen = tmp_path / "elsewhere"
    chosen.mkdir()

    resp = client.patch(
        SET_URL,
        data=json.dumps({"directory": str(chosen)}),
        content_type="application/json",
    )
    assert resp.status_code == 200
    assert resp.json()["data"]["directory"] == str(chosen)
    assert client.get(URL).json()["data"]["directory"] == str(chosen)

    # Resetting (no "directory") clears the preference; the built-in default
    # is something other than the directory just chosen.
    resp = client.patch(SET_URL, data=json.dumps({}), content_type="application/json")
    assert resp.status_code == 200
    reset_to = resp.json()["data"]["directory"]
    assert reset_to != str(chosen)
    assert client.get(URL).json()["data"]["directory"] == reset_to


def test_set_rejected_in_cloud(client, monkeypatch, tmp_path):
    monkeypatch.setenv("CCP4I2_HOME", str(tmp_path))
    monkeypatch.delenv("CCP4I2_LOCAL_SESSION_TOKEN", raising=False)
    monkeypatch.delenv("CCP4I2_PROJECTS_DIR", raising=False)

    resp = client.patch(
        SET_URL,
        data=json.dumps({"directory": str(tmp_path / "x")}),
        content_type="application/json",
    )
    assert resp.status_code == 409
    assert resp.json()["success"] is False
