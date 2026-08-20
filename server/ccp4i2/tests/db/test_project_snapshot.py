"""Tests for the on-disk project snapshot and the project directory registry.

Two things matter here and they pull against each other: the snapshot must be
written whenever the user changes something that exists nowhere else, and it
must not be written once per row when a job finishes and drops thirty files
into the database. Most of these tests therefore count writes, not just check
content.
"""

from pathlib import Path
from shutil import rmtree
from unittest.mock import patch
from xml.etree import ElementTree as ET

from django.conf import settings
from django.db import transaction
from django.test import TestCase, TransactionTestCase, override_settings

from ...db import project_snapshot
from ...db.import_i2xml import import_i2xml_from_file
from ...db.models import File, FileType, Job, Project, ProjectTag
from ...db.project_snapshot import (
    GENERATOR_ELEMENT,
    GENERATOR_MARK,
    PRESERVED_NAME,
    SNAPSHOT_NAME,
    has_snapshot,
    read_registry,
    snapshot_all_projects,
    snapshot_path,
    snapshot_status,
    suspended,
    write_snapshot,
)

PROJECTS_DIR = Path(__file__).parent.parent / "CCP4I2_SNAPSHOT_TEST_DIR"
HOME_DIR = PROJECTS_DIR / "home"


class IsolatedHomeMixin:
    """Keep the registry out of the user's real ~/.ccp4i2-django.

    ``update_registry`` runs from the on-commit callback whether or not a test
    cares about the registry, so pointing it somewhere disposable has to be
    done for every test here, not only the ones that assert on it.
    """

    def setUp(self):
        PROJECTS_DIR.mkdir(parents=True, exist_ok=True)
        HOME_DIR.mkdir(parents=True, exist_ok=True)
        self._home_patcher = patch.object(
            project_snapshot,
            "registry_path",
            return_value=HOME_DIR / project_snapshot.REGISTRY_NAME,
        )
        self._home_patcher.start()
        self.addCleanup(self._home_patcher.stop)
        self.addCleanup(rmtree, PROJECTS_DIR, ignore_errors=True)
        super().setUp()


def make_project(name: str = "Snap") -> Project:
    directory = Path(settings.CCP4I2_PROJECTS_DIR) / name
    directory.mkdir(parents=True, exist_ok=True)
    return Project.objects.create(name=name, directory=str(directory))


def make_file(job: Job, name: str = "XYZOUT.pdb") -> File:
    file_type, _ = FileType.objects.get_or_create(
        name="chemical/x-pdb", defaults={"description": "Model coordinates"}
    )
    return File.objects.create(
        name=name,
        directory=File.Directory.JOB_DIR,
        type=file_type,
        job=job,
        job_param_name="XYZOUT",
    )


@override_settings(CCP4I2_PROJECTS_DIR=PROJECTS_DIR)
class SnapshotContentTest(IsolatedHomeMixin, TestCase):
    """What ends up in the file."""

    def setUp(self):
        super().setUp()
        self.project = make_project()

    def test_snapshot_is_written_into_the_project_directory(self):
        path = write_snapshot(self.project)
        self.assertEqual(path, snapshot_path(self.project))
        self.assertEqual(path.name, SNAPSHOT_NAME)
        self.assertTrue(path.is_file())

    def test_snapshot_captures_what_the_database_alone_holds(self):
        job = Job.objects.create(
            project=self.project,
            number="1",
            task_name="prosmart_refmac",
            title="Refinement against native",
            comments="Looks better than job 3",
            evaluation=Job.Evaluation.BEST,
            status=Job.Status.FINISHED,
        )
        the_file = make_file(job)
        the_file.annotation = "Final model, deposit this one"
        the_file.save()

        root = ET.parse(write_snapshot(self.project)).getroot()

        job_elem = root.find(".//jobTable/job")
        self.assertEqual(job_elem.get("title"), "Refinement against native")
        self.assertEqual(job_elem.get("evaluation"), str(Job.Evaluation.BEST))
        file_elem = root.find(".//fileTable/file")
        self.assertEqual(
            file_elem.get("annotation"), "Final model, deposit this one"
        )

    def test_snapshot_is_written_atomically(self):
        """A reader must never see a half-written snapshot."""
        write_snapshot(self.project)
        original = snapshot_path(self.project).read_bytes()

        with patch("os.replace", side_effect=OSError("disk full")):
            write_snapshot(self.project)

        self.assertEqual(snapshot_path(self.project).read_bytes(), original)
        leftovers = list(Path(self.project.directory).glob(f".{SNAPSHOT_NAME}-*"))
        self.assertEqual(leftovers, [])

    def test_a_missing_project_directory_is_not_an_error(self):
        rmtree(self.project.directory)
        self.assertIsNone(write_snapshot(self.project))

    def test_a_failure_to_write_never_propagates(self):
        """A job must not be reported as failed because a snapshot could not
        be written -- a safety net that breaks the thing it protects is worse
        than no safety net."""
        with patch(
            "ccp4i2.db.project_snapshot.generate_project_xml_tree",
            side_effect=RuntimeError("boom"),
        ):
            self.assertIsNone(write_snapshot(self.project))

    def test_snapshot_all_projects(self):
        make_project("Second")
        result = snapshot_all_projects()
        self.assertEqual(result["written"], 2)
        self.assertEqual(result["skipped"], 0)


@override_settings(CCP4I2_PROJECTS_DIR=PROJECTS_DIR)
class RetrofitTest(IsolatedHomeMixin, TestCase):
    """Projects that grew up without a snapshot have to be able to get one.

    Nothing will ever touch a finished, dormant project again, so signals alone
    would leave exactly the projects most worth protecting unprotected.
    """

    def setUp(self):
        super().setUp()
        self.first = make_project("Alpha")
        self.second = make_project("Beta")

    def test_a_project_starts_unprotected(self):
        self.assertFalse(has_snapshot(self.first))

    def test_retrofit_protects_every_project(self):
        result = snapshot_all_projects(missing_only=True)

        self.assertEqual(result["written"], 2)
        self.assertTrue(has_snapshot(self.first))
        self.assertTrue(has_snapshot(self.second))

    def test_retrofit_is_cheap_to_repeat(self):
        snapshot_all_projects(missing_only=True)
        result = snapshot_all_projects(missing_only=True)

        self.assertEqual(result["written"], 0)
        self.assertEqual(result["already_current"], 2)

    def test_all_rewrites_even_current_snapshots(self):
        snapshot_all_projects(missing_only=True)
        result = snapshot_all_projects(missing_only=False)
        self.assertEqual(result["written"], 2)

    def test_status_reports_what_is_protected(self):
        write_snapshot(self.first)
        by_name = {row["name"]: row for row in snapshot_status()}

        self.assertTrue(by_name["Alpha"]["has_snapshot"])
        self.assertFalse(by_name["Beta"]["has_snapshot"])
        self.assertIsNotNone(by_name["Alpha"]["written"])

    def test_status_notices_a_missing_directory(self):
        rmtree(self.second.directory)
        by_name = {row["name"]: row for row in snapshot_status()}
        self.assertFalse(by_name["Beta"]["directory_exists"])

    def test_a_legacy_snapshot_is_preserved_before_being_replaced(self):
        """A Qt-era DATABASE.db.xml may describe jobs this database never
        imported. Overwriting it blind would destroy the only record of them."""
        legacy = snapshot_path(self.first)
        legacy.write_text("<database><legacy>irreplaceable</legacy></database>")

        write_snapshot(self.first)

        preserved = Path(self.first.directory) / PRESERVED_NAME
        self.assertTrue(preserved.is_file())
        self.assertIn("irreplaceable", preserved.read_text())
        self.assertIn(GENERATOR_ELEMENT, legacy.read_text())

    def test_a_legacy_snapshot_is_only_preserved_once(self):
        """The backup must never be overwritten by one of our own snapshots."""
        legacy = snapshot_path(self.first)
        legacy.write_text("<database><legacy>the original</legacy></database>")

        write_snapshot(self.first)
        write_snapshot(self.first)
        write_snapshot(self.first)

        preserved = Path(self.first.directory) / PRESERVED_NAME
        self.assertIn("the original", preserved.read_text())

    def test_our_own_snapshot_is_not_treated_as_foreign(self):
        write_snapshot(self.first)
        write_snapshot(self.first)
        self.assertFalse((Path(self.first.directory) / PRESERVED_NAME).exists())

    def test_a_legacy_snapshot_mentioning_the_marker_is_still_foreign(self):
        """Found on the first real project this was tried on.

        Its legacy DATABASE.db.xml referenced a directory called
        ``ccp4i2-djangodbapi``, so a substring search for the marker matched
        and the file was taken for one of ours -- meaning it would have been
        overwritten without being preserved.
        """
        legacy = snapshot_path(self.first)
        legacy.write_text(
            "<database><file filename='/Users/me/ccp4i2-djangodbapi/x.pdb'/>"
            "</database>"
        )
        self.assertFalse(has_snapshot(self.first))

        write_snapshot(self.first)

        preserved = Path(self.first.directory) / PRESERVED_NAME
        self.assertTrue(preserved.is_file())
        self.assertIn("ccp4i2-djangodbapi", preserved.read_text())

    def test_status_flags_a_foreign_snapshot(self):
        snapshot_path(self.second).write_text("<database/>")
        by_name = {row["name"]: row for row in snapshot_status()}

        self.assertTrue(by_name["Beta"]["foreign_snapshot"])
        self.assertFalse(by_name["Beta"]["has_snapshot"])


@override_settings(CCP4I2_PROJECTS_DIR=PROJECTS_DIR)
class SnapshotTriggerTest(IsolatedHomeMixin, TransactionTestCase):
    """When it gets written.

    TransactionTestCase, not TestCase: on_commit callbacks never run inside the
    transaction TestCase wraps each test in, so a snapshot scheduled for commit
    would silently never happen.
    """

    def setUp(self):
        super().setUp()
        self.writes = []
        self.patcher = patch.object(
            project_snapshot,
            "write_snapshot",
            side_effect=lambda project: self.writes.append(project.pk),
        )
        self.patcher.start()
        self.project = make_project()
        self.writes.clear()

    def tearDown(self):
        self.patcher.stop()

    def test_creating_a_top_level_job_writes_one_snapshot(self):
        Job.objects.create(
            project=self.project, number="1", task_name="aimless"
        )
        self.assertEqual(len(self.writes), 1)

    def test_creating_a_subjob_writes_nothing(self):
        parent = Job.objects.create(
            project=self.project, number="1", task_name="phaser_pipeline"
        )
        self.writes.clear()

        Job.objects.create(
            project=self.project, number="1.1", task_name="phaser_MR_AUTO",
            parent=parent,
        )
        self.assertEqual(self.writes, [])

    def test_job_finishing_writes_a_snapshot(self):
        job = Job.objects.create(
            project=self.project, number="1", task_name="aimless"
        )
        self.writes.clear()

        job.status = Job.Status.FINISHED
        job.save()
        self.assertEqual(len(self.writes), 1)

    def test_a_subjob_changing_status_writes_nothing(self):
        parent = Job.objects.create(
            project=self.project, number="1", task_name="phaser_pipeline"
        )
        subjob = Job.objects.create(
            project=self.project, number="1.1", task_name="phaser_MR_AUTO",
            parent=parent,
        )
        self.writes.clear()

        subjob.status = Job.Status.FINISHED
        subjob.save()
        self.assertEqual(self.writes, [])

    def test_a_comment_on_a_subjob_is_still_captured(self):
        """Status churn is noise; what the user wrote is not."""
        parent = Job.objects.create(
            project=self.project, number="1", task_name="phaser_pipeline"
        )
        subjob = Job.objects.create(
            project=self.project, number="1.1", task_name="phaser_MR_AUTO",
            parent=parent,
        )
        self.writes.clear()

        subjob.comments = "This is the run that worked"
        subjob.save()
        self.assertEqual(len(self.writes), 1)

    def test_saving_a_job_unchanged_writes_nothing(self):
        job = Job.objects.create(
            project=self.project, number="1", task_name="aimless"
        )
        self.writes.clear()

        job.save()
        self.assertEqual(self.writes, [])

    def test_job_evaluation_writes_a_snapshot(self):
        job = Job.objects.create(
            project=self.project, number="1", task_name="aimless"
        )
        self.writes.clear()

        job.evaluation = Job.Evaluation.REJECTED
        job.save()
        self.assertEqual(len(self.writes), 1)

    def test_creating_files_writes_nothing(self):
        """The gleaner drops many files as a job finishes; the job's own
        transition covers them."""
        job = Job.objects.create(
            project=self.project, number="1", task_name="aimless"
        )
        self.writes.clear()

        for i in range(20):
            make_file(job, f"out{i}.pdb")
        self.assertEqual(self.writes, [])

    def test_annotating_a_file_writes_a_snapshot(self):
        job = Job.objects.create(
            project=self.project, number="1", task_name="aimless"
        )
        the_file = make_file(job)
        self.writes.clear()

        the_file.annotation = "Best map so far"
        the_file.save()
        self.assertEqual(len(self.writes), 1)

    def test_project_description_writes_a_snapshot(self):
        self.project.description = "MDM2 fragment campaign"
        self.project.save()
        self.assertEqual(len(self.writes), 1)

    def test_tagging_a_project_writes_a_snapshot(self):
        tag = ProjectTag.objects.create(text="fragment")
        self.writes.clear()

        tag.projects.add(self.project)
        self.assertEqual(len(self.writes), 1)

    def test_many_changes_in_one_transaction_write_once(self):
        """The point of coalescing: one request, one snapshot."""
        job = Job.objects.create(
            project=self.project, number="1", task_name="aimless"
        )
        files = [make_file(job, f"out{i}.pdb") for i in range(30)]
        self.writes.clear()

        with transaction.atomic():
            job.status = Job.Status.FINISHED
            job.save()
            job.comments = "Done"
            job.save()
            for index, the_file in enumerate(files):
                the_file.annotation = f"Output {index}"
                the_file.save()
            self.project.description = "Updated"
            self.project.save()

        self.assertEqual(len(self.writes), 1)

    def test_suspended_writes_nothing(self):
        with suspended():
            job = Job.objects.create(
                project=self.project, number="1", task_name="aimless"
            )
            job.status = Job.Status.FINISHED
            job.save()
        self.assertEqual(self.writes, [])

    def test_a_rolled_back_transaction_does_not_block_the_next_write(self):
        """De-duplication state must not outlive the transaction that set it.

        If it does, a project touched inside a transaction that then rolls back
        is skipped forever after -- a safety net with a hole nothing reports.
        """
        try:
            with transaction.atomic():
                self.project.description = "Doomed"
                self.project.save()
                raise RuntimeError("rolled back")
        except RuntimeError:
            pass
        self.assertEqual(self.writes, [])

        self.project.description = "This one sticks"
        self.project.save()
        self.assertEqual(len(self.writes), 1)

    def test_writing_resumes_after_suspension(self):
        with suspended():
            Job.objects.create(
                project=self.project, number="1", task_name="aimless"
            )
        Job.objects.create(project=self.project, number="2", task_name="aimless")
        self.assertEqual(len(self.writes), 1)


@override_settings(CCP4I2_PROJECTS_DIR=PROJECTS_DIR)
class RegistryTest(IsolatedHomeMixin, TransactionTestCase):
    """The list of where the projects are, kept outside the database."""

    def setUp(self):
        super().setUp()
        self.home = HOME_DIR

    def test_registry_lists_every_project(self):
        first = make_project("Alpha")
        second = make_project("Beta")

        entries = read_registry()
        directories = {entry["directory"] for entry in entries}
        self.assertEqual(
            directories, {first.directory, second.directory}
        )
        self.assertEqual({e["name"] for e in entries}, {"Alpha", "Beta"})

    def test_registry_records_the_uuid_so_a_rebuild_can_match_projects(self):
        project = make_project("Alpha")
        entry = read_registry()[0]
        self.assertEqual(entry["uuid"], str(project.uuid))

    def test_deleting_a_project_removes_it_from_the_registry(self):
        make_project("Alpha")
        second = make_project("Beta")
        second.delete()

        self.assertEqual([e["name"] for e in read_registry()], ["Alpha"])

    def test_missing_registry_reads_as_empty(self):
        self.assertEqual(read_registry(), [])

    def test_a_corrupt_registry_reads_as_empty_rather_than_raising(self):
        (self.home / project_snapshot.REGISTRY_NAME).write_text("{not json")
        self.assertEqual(read_registry(), [])


@override_settings(CCP4I2_PROJECTS_DIR=PROJECTS_DIR)
class RecoveryRoundTripTest(IsolatedHomeMixin, TestCase):
    """The point of the whole exercise: a snapshot has to be able to restore.

    Writing a file nobody can read back would be worse than writing nothing,
    because it would look like protection. So this loses the database outright
    and rebuilds from what was left on disk.

    A TestCase rather than a TransactionTestCase: the import resolves file
    types against the static rows a data migration installs, and a
    TransactionTestCase truncates those along with everything else. The cost is
    that on-commit callbacks do not fire, so the registry is written explicitly
    below; ``RegistryTest`` covers the signal-driven path.
    """

    def setUp(self):
        super().setUp()
        self.project = make_project("Recoverable")
        self.project.description = "MDM2 fragment campaign, autumn"
        self.project.save()

        self.job = Job.objects.create(
            project=self.project,
            number="1",
            task_name="prosmart_refmac",
            title="Refinement against native data",
            comments="Better than job 3; use this for deposition",
            evaluation=Job.Evaluation.BEST,
            status=Job.Status.FINISHED,
        )
        self.rejected = Job.objects.create(
            project=self.project,
            number="2",
            task_name="prosmart_refmac",
            title="Refinement with TLS",
            comments="TLS made it worse",
            evaluation=Job.Evaluation.REJECTED,
            status=Job.Status.FINISHED,
        )
        self.file = make_file(self.job)
        self.file.annotation = "Final model, deposit this one"
        self.file.sub_type = 1
        self.file.save()

        self.snapshot = write_snapshot(self.project)

    def _lose_the_database(self):
        """Everything the user typed lives here and nowhere else."""
        File.objects.all().delete()
        Job.objects.all().delete()
        ProjectTag.objects.all().delete()
        Project.objects.all().delete()

    def test_the_database_really_is_the_only_copy(self):
        """Guard the premise: if this data were recoverable anyway, the whole
        mechanism would be unnecessary."""
        self._lose_the_database()
        self.assertEqual(Project.objects.count(), 0)
        self.assertEqual(Job.objects.count(), 0)

    def test_a_lost_database_is_rebuilt_from_the_snapshot(self):
        self._lose_the_database()

        with suspended():
            import_i2xml_from_file(self.snapshot)

        project = Project.objects.get(name="Recoverable")
        self.assertEqual(project.directory, str(self.project.directory))
        self.assertEqual(Job.objects.filter(project=project).count(), 2)

    def test_what_the_user_wrote_survives(self):
        self._lose_the_database()
        with suspended():
            import_i2xml_from_file(self.snapshot)

        project = Project.objects.get(name="Recoverable")
        jobs = {job.number: job for job in Job.objects.filter(project=project)}

        self.assertEqual(jobs["1"].title, "Refinement against native data")
        self.assertEqual(jobs["1"].evaluation, Job.Evaluation.BEST)
        self.assertEqual(jobs["2"].title, "Refinement with TLS")
        self.assertEqual(jobs["2"].evaluation, Job.Evaluation.REJECTED)

        restored_file = File.objects.get(name="XYZOUT.pdb")
        self.assertEqual(restored_file.annotation, "Final model, deposit this one")
        self.assertEqual(restored_file.sub_type, 1)

    def test_identity_survives_so_a_rebuild_can_be_matched_up(self):
        """UUIDs, not primary keys: the file references inside every params.xml
        in the project are keyed on these."""
        original_project_uuid = self.project.uuid
        original_job_uuid = self.job.uuid
        original_file_uuid = self.file.uuid

        self._lose_the_database()
        with suspended():
            import_i2xml_from_file(self.snapshot)

        self.assertTrue(Project.objects.filter(uuid=original_project_uuid).exists())
        self.assertTrue(Job.objects.filter(uuid=original_job_uuid).exists())
        self.assertTrue(File.objects.filter(uuid=original_file_uuid).exists())

    def test_the_registry_says_where_to_look(self):
        """With the database gone, nothing else knows the project exists."""
        project_snapshot.update_registry()
        directories = {entry["directory"] for entry in read_registry()}
        self.assertIn(str(self.project.directory), directories)
