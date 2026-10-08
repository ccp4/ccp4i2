"""Deleting a job keeps the files it imported, unless asked not to (#587).

Qt-era CCP4i2 asked "delete imported files too?" when a job that had imported
files was deleted. Keeping them is the default here. A kept file stays in the
project and stays usable, so a job that merely used it is not a dependent and
survives; the importing job stays as a FILE_HOLDER, because a File row cannot
outlive its job (File.job cascades, and an imported file's path is resolved
through its job's project). Deleting the imported files takes every job that
used them with them, which is what deleting did before, unconditionally.
"""

import pytest
from rest_framework.test import APIClient

from ccp4i2.db import models

API = "/api/ccp4i2"


@pytest.fixture
def client(bypass_api_permissions):
    return APIClient()


class World:
    """Job 1 imported X (and made Y). Job 2 used X. Job 3 used Y.

    Job 4 is unrelated and imported nothing.
    """

    def __init__(self, tmp_path):
        self.directory = tmp_path / "proj"
        (self.directory / "CCP4_IMPORTED_FILES").mkdir(parents=True)
        self.project = models.Project.objects.create(
            name="delete_imports", directory=str(self.directory)
        )
        self.ftype, _ = models.FileType.objects.get_or_create(
            name="chemical/x-pdb", defaults={"description": "pdb"}
        )
        self.importer = self.job("1", "import job")
        self.x_user = self.job("2", "uses the import")
        self.y_user = self.job("3", "uses the output")
        self.bystander = self.job("4", "imports nothing")

        self.x = self.imported_file(self.importer, "X.pdb")
        self.y = self.output_file(self.importer, "Y.pdb")
        models.FileUse.objects.create(
            file=self.x, job=self.importer, role=models.FileUse.Role.IN,
            job_param_name="XYZIN",
        )
        models.FileUse.objects.create(
            file=self.x, job=self.x_user, role=models.FileUse.Role.IN,
            job_param_name="XYZIN",
        )
        models.FileUse.objects.create(
            file=self.y, job=self.y_user, role=models.FileUse.Role.IN,
            job_param_name="XYZIN",
        )

    def job(self, number, title, parent=None):
        job = models.Job.objects.create(
            project=self.project, number=number, title=title,
            task_name="prosmart_refmac", parent=parent,
            status=models.Job.Status.FINISHED,
            evaluation=models.Job.Evaluation.GOOD,
        )
        job.directory.mkdir(parents=True)
        (job.directory / "params.xml").write_text("<params/>")
        (job.directory / "log.txt").write_text("log")
        (job.directory / "sub").mkdir()
        (job.directory / "sub" / "scratch").write_text("scratch")
        return job

    def imported_file(self, job, name):
        f = models.File.objects.create(
            job=job, name=name, directory=models.File.Directory.IMPORT_DIR,
            type=self.ftype, job_param_name="XYZIN",
        )
        models.FileImport.objects.create(file=f, name=f"/elsewhere/{name}", checksum="")
        f.path.write_text("imported")
        return f

    def output_file(self, job, name):
        f = models.File.objects.create(
            job=job, name=name, directory=models.File.Directory.JOB_DIR,
            type=self.ftype, job_param_name="XYZOUT",
        )
        models.FileUse.objects.create(
            file=f, job=job, role=models.FileUse.Role.OUT, job_param_name="XYZOUT"
        )
        f.path.write_text("made")
        return f


def exists(job):
    return models.Job.objects.filter(id=job.id).exists()


@pytest.fixture
def world(tmp_path):
    return World(tmp_path)


def test_by_default_the_imported_files_are_kept(client, world):
    x_path, y_path = world.x.path, world.y.path

    resp = client.delete(f"{API}/jobs/{world.importer.id}/")

    assert resp.status_code == 200
    importer = models.Job.objects.get(id=world.importer.id)
    assert importer.status == models.Job.Status.FILE_HOLDER
    assert importer.evaluation == models.Job.Evaluation.UNKNOWN
    # The import is still there, on disk and in the database...
    assert x_path.read_text() == "imported"
    assert models.File.objects.filter(id=world.x.id, job=importer).exists()
    assert models.FileImport.objects.filter(file_id=world.x.id).exists()
    # ...and the job that used it survives, still recorded as using it.
    assert exists(world.x_user)
    assert models.FileUse.objects.filter(file=world.x, job=world.x_user).exists()
    # Everything else the job made goes, with the job that used it.
    assert not y_path.exists()
    assert not models.File.objects.filter(id=world.y.id).exists()
    assert not exists(world.y_user)
    # The holder keeps its use of what it holds (its report lists inputs
    # through those) and nothing else.
    assert list(
        models.FileUse.objects.filter(job=importer).values_list("file_id", "role")
    ) == [(world.x.id, models.FileUse.Role.IN)]
    # The holder keeps its parameters and nothing else of its directory.
    assert sorted(p.name for p in importer.directory.iterdir()) == ["params.xml"]
    assert exists(world.bystander)


def test_a_file_holder_has_a_report_not_a_failure(client, world):
    from ccp4i2.lib.utils.reporting.i2_report import (
        generate_job_report,
        report_is_failure,
    )

    client.delete(f"{API}/jobs/{world.importer.id}/")
    holder = models.Job.objects.get(id=world.importer.id)

    assert not report_is_failure(generate_job_report(holder))


def test_delete_imported_files_takes_them_and_their_users(client, world):
    x_path = world.x.path

    resp = client.delete(
        f"{API}/jobs/{world.importer.id}/?delete_imported_files=true"
    )

    assert resp.status_code == 200
    assert not exists(world.importer)
    assert not exists(world.x_user)
    assert not exists(world.y_user)
    assert not x_path.exists()
    assert not models.File.objects.filter(id=world.x.id).exists()
    assert exists(world.bystander)


def test_deleting_a_file_holder_deletes_its_files(client, world):
    x_path = world.x.path
    client.delete(f"{API}/jobs/{world.importer.id}/")
    assert exists(world.importer)

    resp = client.delete(f"{API}/jobs/{world.importer.id}/")

    assert resp.status_code == 200
    assert not exists(world.importer)
    assert not x_path.exists()
    assert not exists(world.x_user)


def test_a_job_that_imported_nothing_is_simply_deleted(client, world):
    resp = client.delete(f"{API}/jobs/{world.bystander.id}/")

    assert resp.status_code == 200
    assert not exists(world.bystander)
    assert not world.bystander.directory.exists()


def test_preview_describes_both_choices(client, world):
    resp = client.post(
        f"{API}/jobs/delete_preview/", {"job_ids": [world.importer.id]},
        format="json",
    )

    assert resp.status_code == 200
    data = resp.json()["data"]
    assert [j["id"] for j in data["selected_jobs"]] == [world.importer.id]
    files = data["imported_files"]
    assert [f["name"] for f in files] == ["X.pdb"]
    assert files[0]["source_name"] == "X.pdb"
    assert files[0]["job"]["number"] == "1"
    assert [u["number"] for u in files[0]["used_by"]] == ["2"]
    keep = data["keep_imported_files"]
    assert [j["number"] for j in keep["additional_dependents"]] == ["3"]
    assert keep["total_to_delete"] == 2
    delete = data["delete_imported_files"]
    assert [j["number"] for j in delete["additional_dependents"]] == ["2", "3"]
    assert delete["total_to_delete"] == 3


def test_preview_of_a_job_without_imports_offers_no_choice(client, world):
    resp = client.post(
        f"{API}/jobs/delete_preview/", {"job_ids": [world.bystander.id]},
        format="json",
    )

    data = resp.json()["data"]
    assert data["imported_files"] == []
    assert data["keep_imported_files"]["total_to_delete"] == 1


def test_dependent_jobs_follows_the_choice(client, world):
    kept = client.get(f"{API}/jobs/{world.importer.id}/dependent_jobs/").json()
    gone = client.get(
        f"{API}/jobs/{world.importer.id}/dependent_jobs/?delete_imported_files=true"
    ).json()

    assert sorted(j["number"] for j in kept) == ["3"]
    assert sorted(j["number"] for j in gone) == ["2", "3"]


@pytest.mark.parametrize("delete_imported_files", [False, True])
def test_bulk_delete(client, world, delete_imported_files):
    x_path = world.x.path

    resp = client.post(
        f"{API}/jobs/bulk_delete/",
        {
            "job_ids": [world.importer.id, world.bystander.id],
            "delete_imported_files": delete_imported_files,
        },
        format="json",
    )

    assert resp.status_code == 200
    assert not exists(world.bystander)
    assert not exists(world.y_user)
    if delete_imported_files:
        assert not exists(world.importer)
        assert not exists(world.x_user)
        assert not x_path.exists()
    else:
        assert (
            models.Job.objects.get(id=world.importer.id).status
            == models.Job.Status.FILE_HOLDER
        )
        assert exists(world.x_user)
        assert x_path.exists()


def test_a_sub_jobs_imports_are_handed_to_its_top_level_job(client, world):
    pipeline = world.job("5", "pipeline")
    sub = world.job("5.1", "sub-job", parent=pipeline)
    z = world.imported_file(sub, "Z.pdb")

    resp = client.delete(f"{API}/jobs/{pipeline.id}/")

    assert resp.status_code == 200
    assert not exists(sub)
    holder = models.Job.objects.get(id=pipeline.id)
    assert holder.status == models.Job.Status.FILE_HOLDER
    assert models.File.objects.get(id=z.id).job_id == pipeline.id
    assert z.path.exists()
