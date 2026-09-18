"""The New Project dialog must be told where a project will actually go.

The server resolves an unspecified directory to the configured projects
directory — always, regardless of where any previous project landed. The
dialog cannot compute that itself without a second resolver for the same
question, so it asks. The same endpoint changes and resets that default
(desktop only).
"""

import json
import os

import pytest
from rest_framework.test import APIClient

from ccp4i2.api.serializers import default_project_parent
from ccp4i2.db.models import Project

URL = "/api/ccp4i2/config/default-project-parent/"


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


def _set(client, **payload):
    return client.patch(
        URL, data=json.dumps(payload), content_type="application/json"
    )


def test_proposes_the_configured_projects_directory(client, settings):
    resp = client.get(URL)
    assert resp.status_code == 200
    assert resp.json()["data"]["directory"] == str(settings.CCP4I2_PROJECTS_DIR)


def test_ignores_where_the_most_recent_project_landed(client, tmp_path):
    """A one-off project made somewhere unusual must not nudge the default
    proposed to the next New Project dialog — that is what stranded the
    configured default in a hidden home folder nobody could find their way
    back to."""
    elsewhere = tmp_path / "somewhere_else"
    (elsewhere / "proj").mkdir(parents=True)
    Project.objects.create(name="proj", directory=str(elsewhere / "proj"))

    assert client.get(URL).json()["data"]["directory"] != str(elsewhere)


def test_editable_only_on_desktop(client, monkeypatch):
    monkeypatch.setenv("CCP4I2_LOCAL_SESSION_TOKEN", "desktop")
    assert client.get(URL).json()["data"]["editable"] is True

    monkeypatch.delenv("CCP4I2_LOCAL_SESSION_TOKEN", raising=False)
    assert client.get(URL).json()["data"]["editable"] is False


def test_set_and_reset_roundtrip_on_desktop(client, desktop, tmp_path):
    """New Project's "make this the default" and Preferences' Reset both come
    here, and the very next GET — from any window — must reflect it with no
    restart. The serializer must land a project there too, or the dialog and
    the project it creates disagree.

    The launch environment must not outrank the file: while it did, the
    directory written here was read straight back over and both controls
    appeared to do nothing at all.
    """
    chosen = tmp_path / "chosen"

    body = _set(client, directory=str(chosen)).json()["data"]
    assert body["directory"] == str(chosen)
    assert body["default"] == str(desktop)
    assert os.environ["CCP4I2_PROJECTS_DIR"] == str(desktop)
    assert client.get(URL).json()["data"]["directory"] == str(chosen)
    assert default_project_parent() == chosen
    # Created on the way: the New Project form asserts its parent exists.
    assert chosen.is_dir()

    assert _set(client).json()["data"]["directory"] == str(desktop)
    assert client.get(URL).json()["data"]["directory"] == str(desktop)


@pytest.mark.parametrize("bad", ["projects", "{tmp}/a-file/projects"])
def test_set_rejects_a_directory_it_cannot_use(client, desktop, tmp_path, bad):
    (tmp_path / "a-file").write_text("", encoding="utf-8")
    resp = _set(client, directory=bad.format(tmp=tmp_path))
    assert resp.status_code == 400
    assert resp.json()["success"] is False


def test_set_rejected_in_cloud(client, monkeypatch, tmp_path):
    monkeypatch.setenv("CCP4I2_HOME", str(tmp_path))
    monkeypatch.delenv("CCP4I2_LOCAL_SESSION_TOKEN", raising=False)
    monkeypatch.delenv("CCP4I2_PROJECTS_DIR", raising=False)

    resp = _set(client, directory=str(tmp_path / "x"))
    assert resp.status_code == 409
    assert resp.json()["success"] is False
