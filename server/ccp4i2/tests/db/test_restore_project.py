"""Tests for rebuilding the database from project directory snapshots.

The invariant a restore has to hold is not "the import ran without error" but
"the database now agrees with the disk". A restore that renumbered a job, or
rooted a project in the wrong place, would import perfectly cleanly and leave
every job directory reference pointing at nothing. So these tests check the
disk, not just the rows.
"""

import os
from pathlib import Path
from shutil import rmtree
from unittest.mock import patch

from django.conf import settings
from django.test import TestCase, override_settings

from ...db import project_snapshot
from ...db.models import File, FileType, Job, Project
from ...db.project_snapshot import SNAPSHOT_NAME, write_snapshot
from ...db.restore_project import (
    discover_restorable,
    inspect,
    restore_all,
    restore_from_directory,
    restore_from_registry,
    verify_restored,
)

PROJECTS_DIR = Path(__file__).parent.parent / "CCP4I2_RESTORE_TEST_DIR"
HOME_DIR = PROJECTS_DIR / "home"


def build_project_on_disk(name: str, job_numbers=("1", "1.1", "2")) -> Project:
    """A project with real job directories, so disk and database can be compared."""
    directory = Path(settings.CCP4I2_PROJECTS_DIR) / name
    (directory / "CCP4_JOBS").mkdir(parents=True, exist_ok=True)
    (directory / "CCP4_IMPORTED_FILES").mkdir(exist_ok=True)

    project = Project.objects.create(name=name, directory=str(directory))
    project.description = f"{name} description"
    project.save()

    file_type, _ = FileType.objects.get_or_create(
        name="chemical/x-pdb", defaults={"description": "Model coordinates"}
    )
    by_number = {}
    for number in job_numbers:
        parent = by_number.get(number.rsplit(".", 1)[0]) if "." in number else None
        job = Job.objects.create(
            project=project,
            number=number,
            parent=parent,
            task_name="prosmart_refmac",
            title=f"Job {number}",
            comments=f"Notes on job {number}",
            evaluation=Job.Evaluation.BEST,
            status=Job.Status.FINISHED,
        )
        by_number[number] = job
        directory.joinpath("CCP4_JOBS", *[f"job_{p}" for p in number.split(".")]).mkdir(
            parents=True, exist_ok=True
        )
        the_file = File.objects.create(
            name=f"XYZOUT_{number.replace('.', '_')}.pdb",
            directory=File.Directory.JOB_DIR,
            type=file_type,
            job=job,
            job_param_name="XYZOUT",
        )
        the_file.annotation = f"Model from job {number}"
        the_file.save()

    write_snapshot(project)
    return project


class RestoreTestBase(TestCase):
    def setUp(self):
        PROJECTS_DIR.mkdir(parents=True, exist_ok=True)
        HOME_DIR.mkdir(parents=True, exist_ok=True)
        patcher = self.settings(CCP4I2_PROJECTS_DIR=PROJECTS_DIR)
        patcher.enable()
        self.addCleanup(patcher.disable)

        # Restoring refreshes the registry, so every test here can reach it
        # whether or not it cares about it. Point it somewhere disposable, or
        # the suite rewrites the user's real ~/.ccp4i2-django.
        registry = patch.object(
            project_snapshot,
            "registry_path",
            return_value=HOME_DIR / project_snapshot.REGISTRY_NAME,
        )
        registry.start()
        self.addCleanup(registry.stop)
        self.addCleanup(rmtree, PROJECTS_DIR, ignore_errors=True)

    def lose_the_database(self):
        """Delete every row, the way a corrupt database file loses them.

        Under ``suspended()`` deliberately. Deleting a project through the ORM
        rewrites the registry, which would destroy the very thing the restore
        depends on -- but a database file that has been corrupted or removed
        never gets to run those signals, so a test that let them run would be
        testing a situation that cannot arise.
        """
        with project_snapshot.suspended():
            File.objects.all().delete()
            Job.objects.all().delete()
            Project.objects.all().delete()


@override_settings(CCP4I2_PROJECTS_DIR=PROJECTS_DIR)
class InspectTest(RestoreTestBase):
    """Reading a directory without touching anything."""

    def test_inspect_reports_what_the_snapshot_holds(self):
        project = build_project_on_disk("Alpha")
        report = inspect(Path(project.directory))

        self.assertEqual(report.project_name, "Alpha")
        self.assertEqual(report.jobs, 3)
        self.assertEqual(report.files, 3)
        self.assertFalse(report.relocated)
        self.assertIsNone(report.skipped_reason)

    def test_a_directory_with_no_snapshot_is_reported_not_raised(self):
        empty = PROJECTS_DIR / "Empty"
        empty.mkdir()
        report = inspect(empty)
        self.assertIn(SNAPSHOT_NAME, report.skipped_reason)

    def test_an_unreadable_snapshot_is_reported_not_raised(self):
        broken = PROJECTS_DIR / "Broken"
        broken.mkdir()
        (broken / SNAPSHOT_NAME).write_text("<not valid xml")
        report = inspect(broken)
        self.assertIn("not readable", report.skipped_reason)

    def test_inspect_notices_the_project_has_moved(self):
        project = build_project_on_disk("Alpha")
        moved = PROJECTS_DIR / "Moved"
        Path(project.directory).rename(moved)

        report = inspect(moved)
        self.assertTrue(report.relocated)
        self.assertEqual(report.recorded_directory, str(project.directory))


@override_settings(CCP4I2_PROJECTS_DIR=PROJECTS_DIR)
class RestoreOneProjectTest(RestoreTestBase):
    def setUp(self):
        super().setUp()
        self.project = build_project_on_disk("Alpha")
        self.directory = Path(self.project.directory)
        self.uuid = self.project.uuid

    def test_a_lost_database_is_rebuilt_from_the_directory(self):
        self.lose_the_database()

        report = restore_from_directory(self.directory)

        self.assertTrue(report.restored)
        project = Project.objects.get(uuid=self.uuid)
        self.assertEqual(project.name, "Alpha")
        self.assertEqual(Job.objects.filter(project=project).count(), 3)

    def test_the_database_ends_up_agreeing_with_the_disk(self):
        """The post-condition. A renumbered job would import cleanly and leave
        the database pointing at a directory that does not exist."""
        self.lose_the_database()
        restore_from_directory(self.directory)

        project = Project.objects.get(uuid=self.uuid)
        self.assertEqual(verify_restored(project), [])

    def test_job_numbers_come_back_exactly(self):
        self.lose_the_database()
        restore_from_directory(self.directory)

        project = Project.objects.get(uuid=self.uuid)
        numbers = set(
            Job.objects.filter(project=project).values_list("number", flat=True)
        )
        self.assertEqual(numbers, {"1", "1.1", "2"})

    def test_what_the_user_wrote_comes_back(self):
        self.lose_the_database()
        restore_from_directory(self.directory)

        project = Project.objects.get(uuid=self.uuid)
        job = Job.objects.get(project=project, number="1")
        self.assertEqual(job.title, "Job 1")
        self.assertEqual(job.comments, "Notes on job 1")
        self.assertEqual(job.evaluation, Job.Evaluation.BEST)
        self.assertEqual(
            File.objects.get(job=job).annotation, "Model from job 1"
        )

    def test_an_existing_project_is_left_alone(self):
        """Overwriting a live project with an older snapshot would lose work."""
        report = restore_from_directory(self.directory)

        self.assertFalse(report.restored)
        self.assertIn("already in the database", report.skipped_reason)

    def test_replace_rebuilds_an_existing_project(self):
        Job.objects.filter(project=self.project, number="2").delete()
        self.assertEqual(Job.objects.filter(project=self.project).count(), 2)

        report = restore_from_directory(self.directory, replace=True)

        self.assertTrue(report.restored)
        project = Project.objects.get(uuid=self.uuid)
        self.assertEqual(Job.objects.filter(project=project).count(), 3)
        self.assertEqual(verify_restored(project), [])

    def test_a_dry_run_changes_nothing(self):
        self.lose_the_database()

        report = restore_from_directory(self.directory, dry_run=True)

        self.assertFalse(report.restored)
        self.assertTrue(report.dry_run)
        self.assertEqual(report.jobs, 3)
        self.assertEqual(Project.objects.count(), 0)

    def test_a_name_held_by_a_different_project_is_refused(self):
        self.lose_the_database()
        Project.objects.create(name="Alpha", directory=str(PROJECTS_DIR / "Other"))

        report = restore_from_directory(self.directory)

        self.assertFalse(report.restored)
        self.assertIn("already called", report.skipped_reason)

    def test_missing_job_directories_are_reported(self):
        rmtree(self.directory / "CCP4_JOBS" / "job_2")
        self.lose_the_database()

        report = restore_from_directory(self.directory)

        self.assertTrue(report.restored)
        self.assertEqual(report.missing_job_directories, ["2"])
        self.assertTrue(report.warnings)

    def test_restoring_writes_a_fresh_snapshot_and_registry_entry(self):
        self.lose_the_database()
        restore_from_directory(self.directory)

        self.assertTrue((self.directory / SNAPSHOT_NAME).is_file())
        directories = {
            entry["directory"] for entry in project_snapshot.read_registry()
        }
        self.assertIn(str(self.directory), directories)


@override_settings(CCP4I2_PROJECTS_DIR=PROJECTS_DIR)
class RestoreRelocatedProjectTest(RestoreTestBase):
    """The directory is authoritative, even when the snapshot disagrees."""

    def setUp(self):
        super().setUp()
        self.project = build_project_on_disk("Alpha")
        self.original = Path(self.project.directory)
        self.uuid = self.project.uuid
        # A path baked into a job's files, as plugins do.
        self.params = self.original / "CCP4_JOBS" / "job_1" / "params.xml"
        self.params.write_text(
            f"<params><KILLFILEPATH>{self.original}/CCP4_JOBS/job_1/INTERRUPT"
            "</KILLFILEPATH></params>",
            encoding="utf-8",
        )
        write_snapshot(self.project)

        self.moved = PROJECTS_DIR / "SomewhereElse"
        self.original.rename(self.moved)
        self.lose_the_database()

    def test_the_directory_wins_over_the_recorded_path(self):
        restore_from_directory(self.moved)

        project = Project.objects.get(uuid=self.uuid)
        self.assertEqual(project.directory, str(self.moved))
        self.assertEqual(verify_restored(project), [])

    def test_paths_inside_the_files_are_repaired(self):
        report = restore_from_directory(self.moved)

        self.assertTrue(report.relocated)
        self.assertGreater(report.paths_rewritten, 0)
        moved_params = self.moved / "CCP4_JOBS" / "job_1" / "params.xml"
        self.assertIn(str(self.moved), moved_params.read_text())
        self.assertNotIn(str(self.original), moved_params.read_text())

    def test_path_repair_can_be_declined(self):
        report = restore_from_directory(self.moved, repair_paths=False)

        self.assertEqual(report.paths_rewritten, 0)
        moved_params = self.moved / "CCP4_JOBS" / "job_1" / "params.xml"
        self.assertIn(str(self.original), moved_params.read_text())


@override_settings(CCP4I2_PROJECTS_DIR=PROJECTS_DIR)
class RestoreEverythingTest(RestoreTestBase):
    def setUp(self):
        super().setUp()
        self.names = ["Alpha", "Beta", "Gamma"]
        self.projects = [build_project_on_disk(name) for name in self.names]
        project_snapshot.update_registry()

    def test_discover_finds_project_directories(self):
        found = discover_restorable(PROJECTS_DIR)
        self.assertEqual(
            {path.name for path in found}, {"Alpha", "Beta", "Gamma"}
        )

    def test_discover_ignores_directories_without_a_snapshot(self):
        (PROJECTS_DIR / "NotAProject").mkdir()
        found = discover_restorable(PROJECTS_DIR)
        self.assertNotIn("NotAProject", {path.name for path in found})

    def test_a_whole_database_is_rebuilt_by_scanning(self):
        directories = discover_restorable(PROJECTS_DIR)
        self.lose_the_database()

        result = restore_all(directories)

        self.assertEqual(len(result["restored"]), 3)
        self.assertEqual(set(Project.objects.values_list("name", flat=True)),
                         set(self.names))
        for project in Project.objects.all():
            self.assertEqual(verify_restored(project), [])

    def test_a_whole_database_is_rebuilt_from_the_registry(self):
        self.lose_the_database()

        result = restore_from_registry()

        self.assertEqual(result["registry_entries"], 3)
        self.assertEqual(len(result["restored"]), 3)
        self.assertEqual(Project.objects.count(), 3)

    def test_one_bad_project_does_not_abandon_the_rest(self):
        (Path(self.projects[1].directory) / SNAPSHOT_NAME).write_text("<broken")
        self.lose_the_database()

        result = restore_all(discover_restorable(PROJECTS_DIR))

        self.assertEqual(len(result["restored"]), 2)
        self.assertEqual(len(result["skipped"]), 1)
        self.assertEqual(Project.objects.count(), 2)

    def test_rerunning_skips_what_is_already_there(self):
        self.lose_the_database()
        restore_all(discover_restorable(PROJECTS_DIR))

        result = restore_all(discover_restorable(PROJECTS_DIR))

        self.assertEqual(len(result["restored"]), 0)
        self.assertEqual(len(result["skipped"]), 3)
        self.assertEqual(Project.objects.count(), 3)


# Reading a project out of a path means reading the *server's* filesystem, so
# the routes that take one are desktop-only. is_desktop() keys on the token the
# desktop launcher installs; a test exercising those routes has to say it is the
# desktop, and RefusedOffTheDesktopTest below covers the other side.
DESKTOP_ENV = {"CCP4I2_LOCAL_SESSION_TOKEN": "test-session-token"}


@override_settings(CCP4I2_PROJECTS_DIR=PROJECTS_DIR)
@patch.dict(os.environ, DESKTOP_ENV)
class RestoreEndpointTest(RestoreTestBase):
    """The REST surface the desktop app drives."""

    def setUp(self):
        super().setUp()
        from django.contrib.auth import get_user_model
        from rest_framework.test import APIClient

        self.project = build_project_on_disk("Alpha")
        project_snapshot.update_registry()

        self.client = APIClient()
        self.client.force_authenticate(
            get_user_model().objects.create(username="tester")
        )

    def test_restorable_reads_the_registry(self):
        response = self.client.get("/api/ccp4i2/projects/restorable/")

        self.assertEqual(response.status_code, 200, response.data)
        data = response.data["data"]
        self.assertEqual(data["source"]["kind"], "registry")
        self.assertEqual(data["restorable"], 1)
        self.assertEqual(data["candidates"][0]["project_name"], "Alpha")

    def test_restorable_can_scan_a_folder_instead(self):
        response = self.client.get(
            "/api/ccp4i2/projects/restorable/", {"scan": str(PROJECTS_DIR)}
        )

        self.assertEqual(response.status_code, 200, response.data)
        self.assertEqual(response.data["data"]["source"]["kind"], "scan")
        self.assertEqual(response.data["data"]["restorable"], 1)

    def test_restorable_changes_nothing(self):
        self.lose_the_database()
        self.client.get("/api/ccp4i2/projects/restorable/")
        self.assertEqual(Project.objects.count(), 0)

    def test_restore_rebuilds_from_the_registry(self):
        self.lose_the_database()

        response = self.client.post(
            "/api/ccp4i2/projects/restore/", {}, format="json"
        )

        self.assertEqual(response.status_code, 200, response.data)
        self.assertEqual(len(response.data["data"]["restored"]), 1)
        self.assertEqual(Project.objects.count(), 1)
        self.assertEqual(verify_restored(Project.objects.get()), [])

    def test_restore_dry_run_changes_nothing(self):
        self.lose_the_database()

        response = self.client.post(
            "/api/ccp4i2/projects/restore/", {"dry_run": True}, format="json"
        )

        self.assertEqual(response.status_code, 200, response.data)
        self.assertEqual(Project.objects.count(), 0)

    def test_restore_one_directory(self):
        directory = str(self.project.directory)
        self.lose_the_database()

        response = self.client.post(
            "/api/ccp4i2/projects/restore/",
            {"source": "directory", "path": directory},
            format="json",
        )

        self.assertEqual(response.status_code, 200, response.data)
        self.assertEqual(Project.objects.count(), 1)

    def test_a_path_is_required_when_scanning(self):
        response = self.client.post(
            "/api/ccp4i2/projects/restore/", {"source": "scan"}, format="json"
        )
        self.assertEqual(response.status_code, 400)

    def test_an_unknown_source_is_refused(self):
        response = self.client.post(
            "/api/ccp4i2/projects/restore/", {"source": "telepathy"}, format="json"
        )
        self.assertEqual(response.status_code, 400)


@override_settings(CCP4I2_PROJECTS_DIR=PROJECTS_DIR)
@patch.dict(os.environ, {"CCP4I2_LOCAL_SESSION_TOKEN": ""}, clear=False)
class RefusedOffTheDesktopTest(RestoreTestBase):
    """A served deployment must not read projects off the server's disk.

    Registry-driven recovery still has to work everywhere -- it reads only the
    directories this installation already recorded -- but a caller-supplied
    path is a different thing: in a browser the folder the user picked is not
    the folder the server would open, and every other logged-in user shares
    that filesystem.
    """

    def setUp(self):
        super().setUp()
        from django.contrib.auth import get_user_model
        from rest_framework.test import APIClient

        os.environ.pop("CCP4I2_LOCAL_SESSION_TOKEN", None)
        self.project = build_project_on_disk("Alpha")
        project_snapshot.update_registry()

        self.client = APIClient()
        self.client.force_authenticate(
            get_user_model().objects.create(username="tester")
        )

    def test_scanning_a_path_is_refused(self):
        response = self.client.get(
            "/api/ccp4i2/projects/restorable/", {"scan": str(PROJECTS_DIR)}
        )
        self.assertEqual(response.status_code, 409, response.data)

    def test_restoring_a_directory_is_refused(self):
        response = self.client.post(
            "/api/ccp4i2/projects/restore/",
            {"source": "directory", "path": str(self.project.directory)},
            format="json",
        )
        self.assertEqual(response.status_code, 409, response.data)
        self.assertEqual(Project.objects.count(), 1)

    def test_the_registry_route_still_works(self):
        response = self.client.get("/api/ccp4i2/projects/restorable/")
        self.assertEqual(response.status_code, 200, response.data)
        self.assertEqual(response.data["data"]["source"]["kind"], "registry")
