"""The job panel's "i2run command" button renders a configured job.

Every one of these assertions failed before the fix, and every failure was
silent: the endpoint either returned a command consisting of nothing but the
task name (the container it walked had zero children, because
``CContainer.loadContentsFromXml`` on a ``.def.xml`` routes to
``setEtree(..., ignore_missing=True)`` and swallows the mismatch) or raised
``AttributeError: 'CContainer' object has no attribute 'CONTENTS'``. The
existing coverage only asserted that a ``command`` key was present, which both
of those satisfy.

Built from nothing, so it needs no external project zips.
"""
import json
from pathlib import Path

import pytest
from rest_framework.test import APIClient

import ccp4i2
from ccp4i2.db import models
from ccp4i2.lib.utils.containers.get_container import get_job_container

API_PREFIX = "/api/ccp4i2"
GAMMA = Path(ccp4i2.__file__).parent / "demo_data" / "gamma"

# freerflag is ccp4_free (gemmi-native), so this tier runs without CCP4.
TASK = "freerflag"


@pytest.fixture
def client(bypass_api_permissions):
    return APIClient()


@pytest.fixture
def project(tmp_path):
    directory = tmp_path / "i2run_command"
    directory.mkdir()
    return models.Project.objects.create(
        name="i2run_command", directory=str(directory)
    )


@pytest.fixture
def job(client, project):
    response = client.post(
        f"{API_PREFIX}/projects/{project.id}/create_task/",
        data=json.dumps({"task_name": TASK}),
        content_type="application/json",
    )
    assert response.status_code == 200, response.content
    return models.Job.objects.get(id=response.json()["data"]["new_job"]["id"])


def _set(client, job, path, value):
    response = client.post(
        f"{API_PREFIX}/jobs/{job.id}/set_parameter/",
        data=json.dumps({"object_path": path, "value": value}),
        content_type="application/json",
    )
    assert response.status_code == 200, response.content


def _i2run(client, job):
    response = client.get(f"{API_PREFIX}/jobs/{job.id}/i2run_command/")
    assert response.status_code == 200, response.content
    body = response.json()
    assert body.get("success"), body
    return body["data"]


def test_the_container_a_job_renders_from_is_populated(job):
    """The root cause: a container with no children renders no parameters."""
    container = get_job_container(job)
    assert container is not None
    names = {child.objectName() for child in container.children()}
    assert {"inputData", "controlParameters"} <= names, names


def test_the_container_survives_its_plugin_being_collected(job):
    """A parent owns its children here: ``HierarchicalObject.__del__`` calls
    ``destroy()``, which empties ``__dict__``. So a container handed out
    without a reference to the plugin that owns it loses every child the
    moment that plugin is collected --- silently, and only sometimes."""
    import gc

    container = get_job_container(job)
    before = len(container.children())
    assert before > 0
    gc.collect()
    assert len(container.children()) == before


def test_a_configured_job_renders_its_parameters(client, job):
    _set(client, job, f"{TASK}.controlParameters.FRAC", 0.07)

    data = _i2run(client, job)

    # The task and project are the least of it; the parameters are the point.
    assert data["command"].startswith(f"{TASK} --project_name i2run_command")
    assert "--FRAC" in data["command"], data["command"]
    assert '"0.07"' in data["command"], data["command"]


def test_the_rendered_line_is_runnable_as_given(client, job):
    """``command`` is the argument list; ``command_line`` is what a person can
    paste. The dialog shows the latter, so it has to name an entry point that
    exists: not the ``i2run`` console script (shadowed by the legacy Qt
    ``$CCP4/bin/i2run``) and not ``manage.py`` (absent from a packaged tree)."""
    data = _i2run(client, job)

    assert data["command_line"].startswith("ccp4-python -m ccp4i2.cli.i2run ")
    assert data["command_line"].endswith(data["command"])
    # None in a packaged tree, the checkout's server/ directory in a dev one.
    if data["working_directory"] is not None:
        assert (Path(data["working_directory"]) / "manage.py").is_file()


def test_an_input_file_is_rendered(client, job):
    """A file parameter renders as the key=value pairs i2run parses back."""
    _set(
        client,
        job,
        f"{TASK}.inputData.F_SIGF",
        {"fullPath": str(GAMMA / "merged_intensities_Xe.mtz")},
    )

    command = _i2run(client, job)["command"]

    assert "--F_SIGF" in command, command
    assert "merged_intensities_Xe.mtz" in command, command


def test_a_registered_file_renders_as_its_id_alone(client, job, gamma_mtz):
    """dbFileId identifies the row and the row carries the rest, so anything
    else alongside it is redundant -- and makes an inconsistent command
    representable (edit baseName, leave dbFileId, and they disagree).

    contentFlag especially must NOT be restated: it is set by introspection
    when the file is registered, so a hand-editable copy of it in the command
    can contradict the file it names."""
    _set(
        client,
        job,
        f"{TASK}.inputData.F_SIGF",
        {
            "dbFileId": "00477e8c18224a099779746d384aee8e",
            "baseName": "merged_intensities_Xe.mtz",
            "relPath": "CCP4_IMPORTED_FILES",
            "contentFlag": 1,
        },
    )

    command = _i2run(client, job)["command"]

    assert "dbFileId=00477e8c18224a099779746d384aee8e" in command, command
    for redundant in ("baseName=", "relPath=", "project=", "contentFlag="):
        assert redundant not in command, f"{redundant} still rendered: {command}"


def test_the_environment_the_server_runs_with_is_reported(client, job, monkeypatch):
    """The GUI's server can be running on a projects directory a fresh terminal
    knows nothing about. Rendering the command without saying so hands someone
    a line that addresses a DIFFERENT database and looks like it worked."""
    monkeypatch.setenv("CCP4I2_PROJECTS_DIR", "/data/My Projects")
    monkeypatch.setenv("CCP4I2_JOB_TARGET", "local")

    data = _i2run(client, job)

    assert data["environment"]["CCP4I2_PROJECTS_DIR"] == "/data/My Projects"
    assert data["environment"]["CCP4I2_JOB_TARGET"] == "local"
    assert data["platform"]


def test_only_variables_that_are_actually_set_are_reported(client, job, monkeypatch):
    """Reporting a default as though it were a setting is noise, and goes stale
    the day the default changes."""
    monkeypatch.delenv("CCP4I2_HOME", raising=False)
    monkeypatch.delenv("CCP4I2_DB_FILE", raising=False)

    environment = _i2run(client, job)["environment"]

    assert "CCP4I2_HOME" not in environment
    assert "CCP4I2_DB_FILE" not in environment


def test_the_ccp4_setup_script_is_reported_when_there_is_one(client, job, monkeypatch, tmp_path):
    ccp4 = tmp_path / "ccp4"
    (ccp4 / "bin").mkdir(parents=True)
    setup = ccp4 / "bin" / "ccp4.setup-sh"
    setup.write_text("# stub\n")
    monkeypatch.setenv("CCP4", str(ccp4))

    assert _i2run(client, job)["ccp4_setup"] == str(setup)


def test_no_setup_script_is_claimed_when_there_is_none(client, job, monkeypatch, tmp_path):
    """Windows has no ccp4.setup-sh -- the desktop app builds that environment
    itself -- so there is nothing to tell someone to source."""
    monkeypatch.setenv("CCP4", str(tmp_path / "nothing here"))

    assert _i2run(client, job)["ccp4_setup"] is None
