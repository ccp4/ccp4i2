"""Deleting a project removes its records; its directory only on request.

Qt-era CCP4i2 asked whether to delete the files on disk when a project was
deleted, and defaulted to leaving them. DELETE /projects/<id>/ does the same:
the directory is left alone unless ``?delete_files=true``, and even then only a
directory that is recognisably a CCP4i2 project is removed.
"""

import pytest
from rest_framework.test import APIClient

from ccp4i2.db import models
from ccp4i2.db.delete_project_directory import project_directory_removal_problem


@pytest.fixture
def client(bypass_api_permissions):
    return APIClient()


def make_project(tmp_path, name="proj"):
    directory = tmp_path / name
    (directory / "CCP4_JOBS" / "job_1").mkdir(parents=True)
    (directory / "CCP4_JOBS" / "job_1" / "log.txt").write_text("hello")
    (directory / "DATABASE.db.xml").write_text("<db/>")
    project = models.Project.objects.create(name=name, directory=str(directory))
    job = models.Job.objects.create(
        project=project, number="1", title="a job", task_name="refmac"
    )
    return project, job, directory


def test_by_default_the_directory_is_left_on_disk(client, tmp_path):
    project, job, directory = make_project(tmp_path)

    resp = client.delete(f"/api/ccp4i2/projects/{project.id}/")

    assert resp.status_code == 200
    data = resp.json()["data"]
    assert data["deleted"] is True
    assert data["files_deleted"] is False
    assert not models.Project.objects.filter(id=project.id).exists()
    assert not models.Job.objects.filter(id=job.id).exists()
    assert (directory / "CCP4_JOBS" / "job_1" / "log.txt").read_text() == "hello"
    assert (directory / "DATABASE.db.xml").exists()


def test_delete_files_removes_the_directory(client, tmp_path):
    project, job, directory = make_project(tmp_path)

    resp = client.delete(f"/api/ccp4i2/projects/{project.id}/?delete_files=true")

    assert resp.status_code == 200
    data = resp.json()["data"]
    assert data["deleted"] is True
    assert data["files_deleted"] is True
    assert not directory.exists()
    assert not models.Project.objects.filter(id=project.id).exists()


def test_a_directory_that_is_not_a_project_is_never_removed(client, tmp_path):
    """Project.directory is a text column; a stale or mistaken path must not be
    able to take an unrelated directory with it."""
    directory = tmp_path / "just_some_folder"
    directory.mkdir()
    (directory / "notes.txt").write_text("keep me")
    project = models.Project.objects.create(name="odd", directory=str(directory))
    # Creating the row drops a DATABASE.db.xml snapshot into the directory,
    # which is exactly the marker the guard looks for; take it away again to
    # stand for a row whose directory was never a project.
    (directory / "DATABASE.db.xml").unlink()
    assert sorted(p.name for p in directory.iterdir()) == ["notes.txt"]

    resp = client.delete(f"/api/ccp4i2/projects/{project.id}/?delete_files=true")

    assert resp.status_code == 200
    data = resp.json()["data"]
    assert data["deleted"] is True
    assert data["files_deleted"] is False
    assert "does not look like a CCP4i2 project" in data["files_kept_because"]
    assert (directory / "notes.txt").read_text() == "keep me"
    assert not models.Project.objects.filter(id=project.id).exists()


def test_a_project_whose_directory_has_gone_still_deletes(client, tmp_path):
    project = models.Project.objects.create(
        name="vanished", directory=str(tmp_path / "gone")
    )

    resp = client.delete(f"/api/ccp4i2/projects/{project.id}/?delete_files=true")

    assert resp.status_code == 200
    assert resp.json()["data"]["files_kept_because"] == "directory does not exist"
    assert not models.Project.objects.filter(id=project.id).exists()


def test_xdata_rows_do_not_block_deletion(client, tmp_path):
    """XData -> Job is RESTRICT, so it has to be cleared before the cascade."""
    project, job, directory = make_project(tmp_path)
    models.XData.objects.create(data_class="CSomething", xml="<x/>", job=job)

    resp = client.delete(f"/api/ccp4i2/projects/{project.id}/")

    assert resp.status_code == 200
    assert not models.Job.objects.filter(id=job.id).exists()
    assert not models.XData.objects.filter(job_id=job.id).exists()


def test_deleting_a_missing_project_is_a_404(client):
    resp = client.delete("/api/ccp4i2/projects/999999/")
    assert resp.status_code == 404


@pytest.mark.parametrize(
    "make_path, expected",
    [
        (lambda tmp: tmp / "relative_is_impossible_here", "does not exist"),
        (lambda tmp: tmp.anchor, "filesystem root"),
    ],
)
def test_removal_guard_names_its_reason(tmp_path, make_path, expected):
    from pathlib import Path

    problem = project_directory_removal_problem(Path(make_path(tmp_path)))
    assert problem is not None and expected in problem


def test_removal_guard_refuses_the_projects_directory(settings):
    from pathlib import Path

    problem = project_directory_removal_problem(Path(settings.CCP4I2_PROJECTS_DIR))
    assert problem == "path is the projects directory itself"
