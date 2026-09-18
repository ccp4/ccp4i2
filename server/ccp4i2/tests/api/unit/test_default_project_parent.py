"""The New Project dialog must be told where a project will actually go.

The server resolves an unspecified directory to the configured projects
directory — always, regardless of where any previous project landed. The
dialog cannot compute that itself without a second resolver for the same
question, so it asks. It can also change or reset that default (desktop
only) through the paired "set/" endpoint.
"""

import json
import os

import pytest
from rest_framework.test import APIClient

from ccp4i2.db.models import Project


@pytest.fixture
def client(bypass_api_permissions):
    return APIClient()


@pytest.fixture
def desktop(monkeypatch, tmp_path):
    """A desktop session: its own CCP4i2 home, and — as the Electron launcher
    really does it — CCP4I2_PROJECTS_DIR already in the server's environment,
    holding whatever preferences.json said at launch."""
    monkeypatch.setenv("CCP4I2_HOME", str(tmp_path / "home"))
    monkeypatch.setenv("CCP4I2_LOCAL_SESSION_TOKEN", "desktop")
    launched_with = tmp_path / "home" / "projects"
    launched_with.mkdir(parents=True)
    monkeypatch.setenv("CCP4I2_PROJECTS_DIR", str(launched_with))
    return launched_with


URL = "/api/ccp4i2/config/default-project-parent/"
SET_URL = "/api/ccp4i2/config/default-project-parent/set/"


def _set(client, **payload):
    return client.patch(
        SET_URL, data=json.dumps(payload), content_type="application/json"
    )


def test_proposes_the_configured_projects_directory(client, settings):
    resp = client.get(URL)
    assert resp.status_code == 200
    data = resp.json()["data"]
    assert data["directory"] == str(settings.CCP4I2_PROJECTS_DIR)


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


def test_reports_what_a_reset_would_restore(client, desktop):
    """Preferences can only say "this is already the default", and disable its
    Reset button, if it is told what the default is."""
    data = client.get(URL).json()["data"]
    assert data["default"] == str(desktop)
    assert data["directory"] == data["default"]


def test_editable_only_on_desktop(client, monkeypatch):
    monkeypatch.setenv("CCP4I2_LOCAL_SESSION_TOKEN", "desktop")
    assert client.get(URL).json()["data"]["editable"] is True

    monkeypatch.delenv("CCP4I2_LOCAL_SESSION_TOKEN", raising=False)
    assert client.get(URL).json()["data"]["editable"] is False


def test_set_and_reset_roundtrip_on_desktop(client, desktop, tmp_path):
    """The New Project page's "make this the default" checkbox and the
    Preferences page's "Reset to default" both go through this endpoint, and
    the very next GET — from any window — must reflect it immediately, with
    no server restart needed."""
    chosen = tmp_path / "elsewhere"
    chosen.mkdir()

    resp = _set(client, directory=str(chosen))
    assert resp.status_code == 200
    assert resp.json()["data"]["directory"] == str(chosen)
    assert client.get(URL).json()["data"]["directory"] == str(chosen)

    resp = _set(client)
    assert resp.status_code == 200
    assert resp.json()["data"]["directory"] == str(desktop)
    assert client.get(URL).json()["data"]["directory"] == str(desktop)


def test_the_launch_environment_does_not_outrank_the_file(client, desktop, tmp_path):
    """The Electron launcher passes CCP4I2_PROJECTS_DIR to the server it
    spawns, so on the desktop that variable is a copy of preferences.json
    taken at launch, not an instruction that outranks it.

    While it did outrank it, every change here was written to the file and
    then immediately read back over: the Preferences panel's Change and Reset
    buttons, and New Project's "make this the default", all appeared to do
    nothing at all for the life of the app.
    """
    chosen = tmp_path / "chosen"
    chosen.mkdir()
    _set(client, directory=str(chosen))

    assert os.environ["CCP4I2_PROJECTS_DIR"] == str(desktop)
    assert client.get(URL).json()["data"]["directory"] == str(chosen)


def test_a_new_project_goes_to_the_new_default(client, desktop, tmp_path):
    """The point of the setting: the *serializer* must land the project where
    the dialog said it would, without the server being restarted first."""
    from ccp4i2.api.serializers import default_project_parent

    chosen = tmp_path / "chosen"
    chosen.mkdir()
    _set(client, directory=str(chosen))

    assert default_project_parent() == chosen


def test_set_creates_the_directory(client, desktop, tmp_path):
    """Whatever is chosen is about to have a project written into it, and the
    New Project form asserts its parent exists."""
    chosen = tmp_path / "not" / "there" / "yet"
    assert _set(client, directory=str(chosen)).status_code == 200
    assert chosen.is_dir()


def test_set_rejects_a_relative_path(client, desktop):
    resp = _set(client, directory="projects")
    assert resp.status_code == 400
    assert resp.json()["success"] is False


def test_set_rejects_a_path_that_cannot_be_created(client, desktop, tmp_path):
    blocker = tmp_path / "a-file-not-a-directory"
    blocker.write_text("", encoding="utf-8")
    resp = _set(client, directory=str(blocker / "projects"))
    assert resp.status_code == 400
    assert resp.json()["success"] is False


def test_set_rejected_in_cloud(client, monkeypatch, tmp_path):
    monkeypatch.setenv("CCP4I2_HOME", str(tmp_path))
    monkeypatch.delenv("CCP4I2_LOCAL_SESSION_TOKEN", raising=False)
    monkeypatch.delenv("CCP4I2_PROJECTS_DIR", raising=False)

    resp = _set(client, directory=str(tmp_path / "x"))
    assert resp.status_code == 409
    assert resp.json()["success"] is False
