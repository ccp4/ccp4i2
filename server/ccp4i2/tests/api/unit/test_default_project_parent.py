"""The New Project dialog must be told where a project will actually go.

The server resolves an unspecified directory to the parent of the most recently
created project, falling back to the configured projects directory. The dialog
cannot compute that, and reimplementing the rule client-side would give two
answers to one question — so it asks.
"""

import pytest
from rest_framework.test import APIClient

from ccp4i2.db.models import Project


@pytest.fixture
def client(bypass_api_permissions):
    return APIClient()


URL = "/api/ccp4i2/config/default-project-parent/"


def test_falls_back_to_the_configured_store_when_there_are_no_projects(client, settings):
    resp = client.get(URL)
    assert resp.status_code == 200
    assert resp.json()["data"]["directory"] == str(settings.CCP4I2_PROJECTS_DIR)


def test_offers_the_parent_of_the_most_recent_project(client, tmp_path):
    elsewhere = tmp_path / "somewhere_else"
    (elsewhere / "proj").mkdir(parents=True)
    Project.objects.create(name="proj", directory=str(elsewhere / "proj"))

    assert client.get(URL).json()["data"]["directory"] == str(elsewhere)


def test_skips_a_project_whose_directory_has_gone(client, tmp_path, settings):
    """An unplugged external disk must not send the next project somewhere
    unwritable."""
    Project.objects.create(name="vanished", directory=str(tmp_path / "gone" / "p"))

    assert client.get(URL).json()["data"]["directory"] == str(
        settings.CCP4I2_PROJECTS_DIR
    )
