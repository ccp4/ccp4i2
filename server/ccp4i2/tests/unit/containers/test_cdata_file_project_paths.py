"""
How CDataFile stores a path when the job has database context.

A path is only project-relative if it is genuinely inside *this* project.
Deciding that by looking for "CCP4_IMPORTED_FILES" anywhere in the string
re-roots any file chosen from another CCP4i2 project: the real location is
discarded and a path is rebuilt under the current project, where nothing
exists. Diffraction images meet this constantly, because they are referenced
in place and never copied into the project.

These tests drive the database-aware branch of setFullPath(), which only runs
when the plugin carries a _dbProjectId — the reason the bug survived earlier
testing.
"""

import os
import uuid

import pytest

from ccp4i2.core.base_object.cdata_file import CDataFile
from ccp4i2.core.tasks import get_plugin_class


@pytest.fixture(name="project_layout")
def project_layout_fixture(tmp_path):
    """A project, a *different* project, and somewhere with no markers."""
    project = tmp_path / "gammabysad"
    other = tmp_path / "ccp4x_clean" / "xia2xds_test_0"
    job_dir = project / "CCP4_JOBS" / "job_43"
    job_dir.mkdir(parents=True)

    paths = {}
    for name, base in (("inside", project), ("outside", other)):
        sweep = base / "CCP4_IMPORTED_FILES" / "th_8_2"
        sweep.mkdir(parents=True)
        paths[name] = sweep / "th_8_2_0001.pdb"
        paths[name].write_text("")

    unmarked = tmp_path / "elsewhere" / "images"
    unmarked.mkdir(parents=True)
    paths["unmarked"] = unmarked / "th_8_2_0001.pdb"
    paths["unmarked"].write_text("")

    paths["project"] = project
    paths["job_dir"] = job_dir
    return paths


@pytest.fixture(name="file_in_db_context")
def file_in_db_context_fixture(project_layout, monkeypatch):
    """A CDataFile belonging to a job in a project, as at parameter-set time."""

    class FakeDbHandler:
        def getProjectDirectory(self, _project_id):
            return str(project_layout["project"])

    monkeypatch.setattr(
        CDataFile, "_get_db_handler", lambda self: FakeDbHandler()
    )

    # The parent chain holds the plugin weakly, so the test must keep a strong
    # reference: without one the plugin is collected, _find_plugin_parent()
    # returns None and the database-aware branch never runs — which is exactly
    # how this bug stayed invisible.
    plugins = []

    def make():
        plugin = get_plugin_class("chainsaw")(
            workDirectory=str(project_layout["job_dir"])
        )
        plugin._dbProjectId = str(uuid.uuid4())
        plugins.append(plugin)
        return plugin.container.inputData.XYZIN

    return make


def test_file_from_another_project_keeps_its_absolute_path(
    project_layout, file_in_db_context
):
    """The job 43 regression: a path is not ours just because it says
    CCP4_IMPORTED_FILES."""
    xyzin = file_in_db_context()
    xyzin.set(str(project_layout["outside"]))

    assert str(xyzin.relPath) == ""
    assert str(xyzin.baseName) == str(project_layout["outside"])
    assert os.path.realpath(xyzin.getFullPath()) == os.path.realpath(
        project_layout["outside"]
    )


def test_file_inside_project_keeps_its_sub_directory(
    project_layout, file_in_db_context
):
    """Images sit in a per-sweep sub-directory, which must survive."""
    xyzin = file_in_db_context()
    xyzin.set(str(project_layout["inside"]))

    assert str(xyzin.relPath) == os.path.join("CCP4_IMPORTED_FILES", "th_8_2")
    assert str(xyzin.baseName) == "th_8_2_0001.pdb"
    assert os.path.realpath(xyzin.getFullPath()) == os.path.realpath(
        project_layout["inside"]
    )


def test_file_with_no_project_markers_keeps_its_absolute_path(
    project_layout, file_in_db_context
):
    xyzin = file_in_db_context()
    xyzin.set(str(project_layout["unmarked"]))

    assert str(xyzin.baseName) == str(project_layout["unmarked"])
    assert os.path.realpath(xyzin.getFullPath()) == os.path.realpath(
        project_layout["unmarked"]
    )
