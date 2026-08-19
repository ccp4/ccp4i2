"""Tests for moving a project on disk and rebasing its absolute paths."""

from pathlib import Path
from shutil import rmtree
from tempfile import mkdtemp

from django.conf import settings
from django.test import TestCase, override_settings

from ...db.models import Job, Project
from ...db.move_project import (
    JOURNAL_NAME,
    MoveProjectError,
    apply_rebase,
    directory_is_missing,
    dry_run_move,
    find_stale_roots,
    looks_like_a_project,
    move_project,
    path_variants,
    plan_rebase,
    plan_root_rebase,
    rebase_projects_root,
    relocate_project,
    repair_project_paths,
    substitution_map,
)


def make_project_tree(root: Path, project_root_text: str) -> None:
    """Build a miniature project containing every class of path leakage."""
    job = root / "CCP4_JOBS" / "job_1"
    job.mkdir(parents=True)
    (root / "CCP4_IMPORTED_FILES").mkdir()

    # Portable: relPath + baseName. Must survive untouched.
    (job / "params.xml").write_text(
        "<params>\n"
        "  <CPdbDataFile>\n"
        "    <baseName>XYZOUT.pdb</baseName>\n"
        "    <relPath>CCP4_JOBS/job_1</relPath>\n"
        "  </CPdbDataFile>\n"
        f"  <KILLFILEPATH>{project_root_text}/CCP4_JOBS/job_1/INTERRUPT</KILLFILEPATH>\n"
        "</params>\n",
        encoding="utf-8",
    )
    # Program-authored XML.
    (job / "program.xml").write_text(
        f"<AIMLESS><MTZmergedfilename>{project_root_text}/CCP4_JOBS/job_1/HKLOUT.mtz"
        "</MTZmergedfilename></AIMLESS>\n",
        encoding="utf-8",
    )
    # Generated script re-read by a viewer.
    (job / "script.py").write_text(
        f"filePath = r'{project_root_text}/CCP4_JOBS/job_1/XYZOUT.pdb'\n",
        encoding="utf-8",
    )
    (job / "scene_1.scene.xml").write_text(
        f"<scene><filename>{project_root_text}/CCP4_JOBS/job_1/XYZOUT.pdb</filename></scene>\n",
        encoding="utf-8",
    )
    # Provenance: must be left exactly as it was.
    (job / "log.txt").write_text(
        f"Opening {project_root_text}/CCP4_JOBS/job_1/HKLIN.mtz\n", encoding="utf-8"
    )
    (job / "stdout.txt").write_text(
        f"wrote {project_root_text}/CCP4_JOBS/job_1/HKLOUT.mtz\n", encoding="utf-8"
    )
    # Regenerable caches: must be deleted.
    (job / "report_xml.xml").write_text(
        f"<report>{project_root_text}/CCP4_JOBS/job_1/XYZOUT.pdb</report>\n",
        encoding="utf-8",
    )
    (job / "report.html").write_text(
        f"<a href='{project_root_text}/x.pdb'>x</a>\n", encoding="utf-8"
    )
    # Binary: must never be opened for rewriting.
    (job / "HKLOUT.mtz").write_bytes(
        b"MTZ \x00\x00\x00" + project_root_text.encode() + b"\x00padding"
    )


class RebaseHelpersTest(TestCase):
    """The pure path-substitution helpers, no database involved."""

    def test_path_variants_cover_windows_spellings(self):
        variants = path_variants(r"C:\Users\me\Proj")
        self.assertEqual(variants["windows"], r"C:\Users\me\Proj")
        self.assertEqual(variants["posix"], "C:/Users/me/Proj")
        self.assertEqual(variants["windows_escaped"], r"C:\\Users\\me\\Proj")

    def test_trailing_separator_ignored(self):
        self.assertEqual(path_variants("/a/b/"), path_variants("/a/b"))

    def test_roots_of_different_shapes_pair_by_form_not_position(self):
        """A POSIX source and a Windows destination must not cross-wire."""
        pairs = dict(substitution_map("/mnt/data/Proj", r"D:\Projects\Proj"))
        self.assertEqual(pairs["/mnt/data/Proj"], "D:/Projects/Proj")
        self.assertEqual(pairs[r"\mnt\data\Proj"], r"D:\Projects\Proj")

    def test_substitutions_apply_longest_first(self):
        pairs = substitution_map(r"C:\old\proj", r"D:\new\proj")
        lengths = [len(old) for old, _ in pairs]
        self.assertEqual(lengths, sorted(lengths, reverse=True))

    def test_doubled_backslash_form_survives_substitution(self):
        old, new = r"C:\old\proj", r"D:\new\proj"
        text = r"path = 'C:\\old\\proj\\CCP4_JOBS'"
        for needle, replacement in substitution_map(old, new):
            text = text.replace(needle, replacement)
        self.assertEqual(text, r"path = 'D:\\new\\proj\\CCP4_JOBS'")


class PlanRebaseTest(TestCase):
    """Classification of a real-shaped project tree."""

    def setUp(self):
        self.root = Path(mkdtemp())
        self.old_root = "/old/place/MyProject"
        make_project_tree(self.root, self.old_root)

    def tearDown(self):
        rmtree(self.root, ignore_errors=True)

    def test_plan_selects_functional_files_only(self):
        plan = plan_rebase(self.root, self.old_root, "/new/place/MyProject")
        rewritten = {candidate.path.name for candidate in plan.rewrites}
        self.assertEqual(
            rewritten,
            {"params.xml", "program.xml", "script.py", "scene_1.scene.xml"},
        )

    def test_provenance_is_skipped(self):
        plan = plan_rebase(self.root, self.old_root, "/new/place/MyProject")
        skipped = {path.name for path in plan.skipped_provenance}
        self.assertIn("log.txt", skipped)
        self.assertIn("stdout.txt", skipped)

    def test_caches_are_listed_for_deletion(self):
        plan = plan_rebase(self.root, self.old_root, "/new/place/MyProject")
        caches = {path.name for path in plan.caches}
        self.assertEqual(caches, {"report_xml.xml", "report.html"})

    def test_binary_is_skipped(self):
        plan = plan_rebase(self.root, self.old_root, "/new/place/MyProject")
        self.assertIn("HKLOUT.mtz", {path.name for path in plan.skipped_binary})

    def test_apply_rewrites_and_deletes(self):
        new_root = "/new/place/MyProject"
        plan = plan_rebase(self.root, self.old_root, new_root)
        rewritten = apply_rebase(plan)

        self.assertEqual(len(rewritten), 4)
        job = self.root / "CCP4_JOBS" / "job_1"

        params = (job / "params.xml").read_text()
        self.assertIn(f"{new_root}/CCP4_JOBS/job_1/INTERRUPT", params)
        self.assertNotIn(self.old_root, params)
        # The portable half is untouched.
        self.assertIn("<relPath>CCP4_JOBS/job_1</relPath>", params)

        self.assertIn(new_root, (job / "script.py").read_text())
        self.assertIn(new_root, (job / "scene_1.scene.xml").read_text())

        # Provenance preserved verbatim.
        self.assertIn(self.old_root, (job / "log.txt").read_text())
        self.assertIn(self.old_root, (job / "stdout.txt").read_text())

        # Caches gone, binary untouched.
        self.assertFalse((job / "report_xml.xml").exists())
        self.assertFalse((job / "report.html").exists())
        self.assertIn(
            self.old_root.encode(), (job / "HKLOUT.mtz").read_bytes()
        )


class StaleRootTest(TestCase):
    """Detection of roots left behind by an earlier move or import."""

    def setUp(self):
        self.root = Path(mkdtemp())
        self.old_root = "/Volumes/OldDrive/MyProject"
        make_project_tree(self.root, self.old_root)

    def tearDown(self):
        rmtree(self.root, ignore_errors=True)

    def test_finds_the_root_the_files_point_at(self):
        found = find_stale_roots(self.root, str(self.root))
        self.assertIn(self.old_root, found)
        self.assertGreater(found[self.old_root], 0)

    def test_current_root_is_not_reported(self):
        make_project_tree(self.root / "nested", str(self.root / "nested"))
        found = find_stale_roots(self.root / "nested", str(self.root / "nested"))
        self.assertNotIn(str(self.root / "nested"), found)

    def test_windows_root_is_recovered(self):
        root = Path(mkdtemp())
        try:
            (root / "CCP4_JOBS").mkdir()
            (root / "CCP4_JOBS" / "params.xml").write_text(
                r"<p>C:\Data\Projects\Foo\CCP4_JOBS\job_1\XYZOUT.pdb</p>",
                encoding="utf-8",
            )
            found = find_stale_roots(root, str(root))
            self.assertIn(r"C:\Data\Projects\Foo", found)
        finally:
            rmtree(root, ignore_errors=True)

    def test_nothing_reported_when_everything_is_relative(self):
        root = Path(mkdtemp())
        try:
            (root / "CCP4_JOBS").mkdir()
            (root / "CCP4_JOBS" / "params.xml").write_text(
                "<p><relPath>CCP4_JOBS/job_1</relPath></p>", encoding="utf-8"
            )
            self.assertEqual(find_stale_roots(root, str(root)), {})
        finally:
            rmtree(root, ignore_errors=True)


@override_settings(
    CCP4I2_PROJECTS_DIR=Path(__file__).parent.parent / "CCP4I2_MOVE_TEST_DIR"
)
class MoveProjectTest(TestCase):
    def setUp(self):
        self.projects_dir = Path(settings.CCP4I2_PROJECTS_DIR)
        self.projects_dir.mkdir(parents=True, exist_ok=True)
        self.source = self.projects_dir / "MyProject"
        self.source.mkdir()
        make_project_tree(self.source, str(self.source))
        self.project = Project.objects.create(
            name="MyProject", directory=str(self.source)
        )

    def tearDown(self):
        rmtree(self.projects_dir, ignore_errors=True)

    def test_dry_run_changes_nothing(self):
        destination = self.projects_dir / "Moved"
        summary = dry_run_move(self.project, destination)

        self.assertEqual(summary["destination"], str(destination))
        self.assertEqual(len(summary["rewrites"]), 4)
        self.assertTrue(self.source.is_dir())
        self.assertFalse(destination.exists())
        self.project.refresh_from_db()
        self.assertEqual(self.project.directory, str(self.source))

    def test_move_relocates_and_rebases(self):
        destination = self.projects_dir / "Moved"
        summary = move_project(self.project, destination)

        self.assertFalse(self.source.exists())
        self.assertTrue(destination.is_dir())
        self.assertEqual(summary["rewritten"], 4)

        self.project.refresh_from_db()
        self.assertEqual(self.project.directory, str(destination))

        params = (destination / "CCP4_JOBS" / "job_1" / "params.xml").read_text()
        self.assertIn(str(destination), params)
        self.assertNotIn(str(self.source), params)

        # No journal left behind on success.
        self.assertFalse((destination / JOURNAL_NAME).exists())

    def test_move_refuses_existing_destination(self):
        destination = self.projects_dir / "Occupied"
        destination.mkdir()
        with self.assertRaises(MoveProjectError):
            move_project(self.project, destination)
        self.assertTrue(self.source.is_dir())

    def test_move_refuses_destination_inside_project(self):
        with self.assertRaises(MoveProjectError):
            move_project(self.project, self.source / "CCP4_JOBS" / "inner")

    def test_move_refuses_while_a_job_is_running(self):
        for number, status in (
            ("1", Job.Status.RUNNING),
            ("2", Job.Status.QUEUED),
            ("3", Job.Status.RUNNING_REMOTELY),
        ):
            with self.subTest(status=status):
                job = Job.objects.create(
                    project=self.project,
                    number=number,
                    status=status,
                    task_name="aimless",
                )
                with self.assertRaises(MoveProjectError) as caught:
                    move_project(self.project, self.projects_dir / "Moved")
                self.assertIn("running", str(caught.exception).lower())
                self.assertTrue(self.source.is_dir())
                job.delete()

    def test_unrun_tasks_do_not_block_a_move(self):
        """PENDING is a task being edited, not one that is doing anything."""
        for number in ("12", "13", "14", "15", "17"):
            Job.objects.create(
                project=self.project,
                number=number,
                status=Job.Status.PENDING,
                task_name="aimless",
            )
        destination = self.projects_dir / "Moved"

        move_project(self.project, destination)

        self.assertTrue(destination.is_dir())
        self.project.refresh_from_db()
        self.assertEqual(self.project.directory, str(destination))

    def test_finished_and_failed_jobs_do_not_block_a_move(self):
        for number, status in (
            ("1", Job.Status.FINISHED),
            ("2", Job.Status.FAILED),
            ("3", Job.Status.INTERRUPTED),
            ("4", Job.Status.UNKNOWN),
        ):
            Job.objects.create(
                project=self.project,
                number=number,
                status=status,
                task_name="aimless",
            )

        move_project(self.project, self.projects_dir / "Moved")

        self.assertTrue((self.projects_dir / "Moved").is_dir())

    def test_repair_rewrites_without_moving(self):
        """The renamed-mount case: directory intact, its old path stale."""
        stale_root = "/Volumes/OldDriveName/MyProject"
        params = self.source / "CCP4_JOBS" / "job_1" / "params.xml"
        params.write_text(
            f"<params><KILLFILEPATH>{stale_root}/INTERRUPT</KILLFILEPATH></params>\n",
            encoding="utf-8",
        )

        summary = repair_project_paths(self.project, stale_root)

        self.assertEqual(summary["rewritten"], 1)
        self.assertIn(str(self.source), params.read_text())
        self.assertNotIn(stale_root, params.read_text())
        self.project.refresh_from_db()
        self.assertEqual(self.project.directory, str(self.source))

    def test_move_reports_roots_it_could_not_know_about(self):
        """A project moved twice still points at its first home; say so."""
        params = self.source / "CCP4_JOBS" / "job_1" / "params.xml"
        params.write_text(
            "<params><f>/Volumes/FirstHome/MyProject/CCP4_JOBS/job_1/X.pdb</f></params>",
            encoding="utf-8",
        )
        summary = move_project(self.project, self.projects_dir / "Moved")
        self.assertIn("/Volumes/FirstHome/MyProject", summary["stale_roots"])

    def test_directory_cannot_be_changed_through_the_serializer(self):
        from ...api.serializers import ProjectSerializer

        serializer = ProjectSerializer(
            self.project,
            data={"directory": str(self.projects_dir / "Elsewhere")},
            partial=True,
        )
        self.assertFalse(serializer.is_valid())
        self.assertIn("directory", serializer.errors)

    def test_serializer_still_accepts_an_unchanged_directory(self):
        from ...api.serializers import ProjectSerializer

        serializer = ProjectSerializer(
            self.project,
            data={"directory": self.project.directory, "description": "hello"},
            partial=True,
        )
        self.assertTrue(serializer.is_valid(), serializer.errors)

    def test_endpoints_plan_then_move(self):
        """The REST surface, end to end: preview, move, then repair."""
        from django.contrib.auth import get_user_model
        from rest_framework.test import APIClient

        client = APIClient()
        client.force_authenticate(get_user_model().objects.create(username="tester"))
        destination = self.projects_dir / "Moved"

        preview = client.post(
            f"/api/ccp4i2/projects/{self.project.id}/move/",
            {"directory": str(destination), "dry_run": True},
            format="json",
        )
        self.assertEqual(preview.status_code, 200, preview.data)
        self.assertEqual(len(preview.data["data"]["rewrites"]), 4)
        self.assertTrue(self.source.is_dir())

        moved = client.post(
            f"/api/ccp4i2/projects/{self.project.id}/move/",
            {"directory": str(destination)},
            format="json",
        )
        self.assertEqual(moved.status_code, 200, moved.data)
        self.assertEqual(moved.data["data"]["rewritten"], 4)
        self.assertTrue(destination.is_dir())
        self.assertFalse(self.source.exists())

    def test_move_endpoint_reports_a_refusal_as_a_400(self):
        from django.contrib.auth import get_user_model
        from rest_framework.test import APIClient

        client = APIClient()
        client.force_authenticate(get_user_model().objects.create(username="tester"))
        occupied = self.projects_dir / "Occupied"
        occupied.mkdir()

        response = client.post(
            f"/api/ccp4i2/projects/{self.project.id}/move/",
            {"directory": str(occupied)},
            format="json",
        )
        self.assertEqual(response.status_code, 400)
        self.assertTrue(self.source.is_dir())

    def test_stale_roots_endpoint(self):
        from django.contrib.auth import get_user_model
        from rest_framework.test import APIClient

        client = APIClient()
        client.force_authenticate(get_user_model().objects.create(username="tester"))
        params = self.source / "CCP4_JOBS" / "job_1" / "params.xml"
        params.write_text(
            "<params><f>/Volumes/Gone/MyProject/CCP4_JOBS/job_1/X.pdb</f></params>",
            encoding="utf-8",
        )

        response = client.get(
            f"/api/ccp4i2/projects/{self.project.id}/stale_roots/"
        )
        self.assertEqual(response.status_code, 200, response.data)
        self.assertIn("/Volumes/Gone/MyProject", response.data["data"]["stale_roots"])

    def test_repair_dry_run_changes_nothing(self):
        stale_root = "/Volumes/OldDriveName/MyProject"
        params = self.source / "CCP4_JOBS" / "job_1" / "params.xml"
        params.write_text(
            f"<params><KILLFILEPATH>{stale_root}/INTERRUPT</KILLFILEPATH></params>\n",
            encoding="utf-8",
        )

        summary = repair_project_paths(self.project, stale_root, dry_run=True)

        self.assertEqual(len(summary["rewrites"]), 1)
        self.assertEqual(summary["rewritten"], 0)
        self.assertIn(stale_root, params.read_text())


@override_settings(
    CCP4I2_PROJECTS_DIR=Path(__file__).parent.parent / "CCP4I2_RELOCATE_TEST_DIR"
)
class RelocateProjectTest(TestCase):
    """The database is intact; the directory is not where it says it is.

    A renamed drive or a projects folder moved outside CCP4i2. Nothing needs
    moving - the record needs correcting, and the paths inside rebasing.
    """

    def setUp(self):
        self.projects_dir = Path(settings.CCP4I2_PROJECTS_DIR)
        self.projects_dir.mkdir(parents=True, exist_ok=True)

        # Where the project really is.
        self.actual = self.projects_dir / "NewDriveName" / "MyProject"
        self.actual.mkdir(parents=True)
        self.stale = "/Volumes/OldDriveName/MyProject"
        make_project_tree(self.actual, self.stale)

        # What the database still believes.
        self.project = Project.objects.create(
            name="MyProject", directory=self.stale
        )

    def tearDown(self):
        rmtree(self.projects_dir, ignore_errors=True)

    def test_missing_directory_is_detected(self):
        self.assertTrue(directory_is_missing(self.project))

    def test_relocate_repoints_and_rebases(self):
        summary = relocate_project(self.project, self.actual)

        self.project.refresh_from_db()
        self.assertEqual(self.project.directory, str(self.actual))
        self.assertTrue(summary["directory_was_missing"])
        self.assertEqual(summary["rewritten"], 4)

        params = (self.actual / "CCP4_JOBS" / "job_1" / "params.xml").read_text()
        self.assertIn(str(self.actual), params)
        self.assertNotIn(self.stale, params)
        self.assertEqual(summary["stale_roots"], {})

    def test_relocate_moves_nothing(self):
        relocate_project(self.project, self.actual)
        self.assertTrue((self.actual / "CCP4_JOBS" / "job_1").is_dir())

    def test_dry_run_changes_nothing(self):
        summary = relocate_project(self.project, self.actual, dry_run=True)

        self.assertEqual(len(summary["rewrites"]), 4)
        self.assertEqual(summary["rewritten"], 0)
        self.project.refresh_from_db()
        self.assertEqual(self.project.directory, self.stale)

    def test_refuses_a_folder_that_is_not_a_project(self):
        elsewhere = self.projects_dir / "Downloads"
        elsewhere.mkdir()
        with self.assertRaises(MoveProjectError) as caught:
            relocate_project(self.project, elsewhere)
        self.assertIn("does not look like", str(caught.exception))

    def test_refuses_a_directory_another_project_claims(self):
        Project.objects.create(name="Other", directory=str(self.actual))
        with self.assertRaises(MoveProjectError):
            relocate_project(self.project, self.actual)

    def test_looks_like_a_project(self):
        self.assertTrue(looks_like_a_project(self.actual))
        self.assertFalse(looks_like_a_project(self.projects_dir))

    def test_repair_points_at_relocate_when_the_directory_is_gone(self):
        """The message must not send the user somewhere they cannot go."""
        with self.assertRaises(MoveProjectError) as caught:
            repair_project_paths(self.project, "/somewhere/else")
        self.assertIn("re-point", str(caught.exception).lower())


@override_settings(
    CCP4I2_PROJECTS_DIR=Path(__file__).parent.parent / "CCP4I2_ROOT_REBASE_TEST_DIR"
)
class RebaseProjectsRootTest(TestCase):
    """One drive rename, every project stale at once."""

    def setUp(self):
        self.base = Path(settings.CCP4I2_PROJECTS_DIR)
        self.base.mkdir(parents=True, exist_ok=True)
        self.old_root = "/Volumes/OldName/CCP4I2_PROJECTS"
        self.new_root = self.base / "NewName" / "CCP4I2_PROJECTS"
        self.new_root.mkdir(parents=True)

        self.projects = []
        for name in ("Alpha", "Beta", "Gamma"):
            directory = self.new_root / name
            directory.mkdir()
            make_project_tree(directory, f"{self.old_root}/{name}")
            self.projects.append(
                Project.objects.create(
                    name=name, directory=f"{self.old_root}/{name}"
                )
            )
        # One project that was never on the renamed drive.
        self.untouched_dir = self.base / "Local" / "Delta"
        self.untouched_dir.mkdir(parents=True)
        make_project_tree(self.untouched_dir, str(self.untouched_dir))
        self.untouched = Project.objects.create(
            name="Delta", directory=str(self.untouched_dir)
        )

    def tearDown(self):
        rmtree(self.base, ignore_errors=True)

    def test_plan_maps_projects_onto_the_new_root(self):
        plan = plan_root_rebase(self.old_root, str(self.new_root))

        self.assertEqual({e["name"] for e in plan["matched"]}, {"Alpha", "Beta", "Gamma"})
        self.assertEqual(plan["missing"], [])
        self.assertEqual({e["name"] for e in plan["unaffected"]}, {"Delta"})

    def test_plan_reports_projects_it_cannot_find(self):
        rmtree(self.new_root / "Beta")
        plan = plan_root_rebase(self.old_root, str(self.new_root))

        self.assertEqual({e["name"] for e in plan["matched"]}, {"Alpha", "Gamma"})
        self.assertEqual([e["name"] for e in plan["missing"]], ["Beta"])
        self.assertEqual(plan["missing"][0]["reason"], "not found")

    def test_dry_run_changes_nothing(self):
        summary = rebase_projects_root(
            self.old_root, str(self.new_root), dry_run=True
        )
        self.assertEqual(summary["relocated"], [])
        for project in self.projects:
            project.refresh_from_db()
            self.assertTrue(project.directory.startswith(self.old_root))

    def test_all_three_projects_are_repointed(self):
        summary = rebase_projects_root(self.old_root, str(self.new_root))

        self.assertEqual(len(summary["relocated"]), 3)
        self.assertEqual(summary["failed"], [])
        for project in self.projects:
            project.refresh_from_db()
            self.assertEqual(
                project.directory, str(self.new_root / project.name)
            )
            self.assertFalse(directory_is_missing(project))

    def test_unrelated_project_is_left_alone(self):
        rebase_projects_root(self.old_root, str(self.new_root))
        self.untouched.refresh_from_db()
        self.assertEqual(self.untouched.directory, str(self.untouched_dir))

    def test_a_missing_project_does_not_abandon_the_others(self):
        rmtree(self.new_root / "Beta")
        summary = rebase_projects_root(self.old_root, str(self.new_root))

        self.assertEqual(len(summary["relocated"]), 2)
        self.assertEqual([e["name"] for e in summary["missing"]], ["Beta"])
        beta = Project.objects.get(name="Beta")
        self.assertTrue(directory_is_missing(beta))

    def test_paths_inside_the_files_are_rebased_too(self):
        rebase_projects_root(self.old_root, str(self.new_root))
        params = (
            self.new_root / "Alpha" / "CCP4_JOBS" / "job_1" / "params.xml"
        ).read_text()
        self.assertIn(str(self.new_root / "Alpha"), params)
        self.assertNotIn(self.old_root, params)

    def test_endpoints(self):
        from django.contrib.auth import get_user_model
        from rest_framework.test import APIClient

        client = APIClient()
        client.force_authenticate(get_user_model().objects.create(username="tester"))

        broken = client.get("/api/ccp4i2/projects/missing_directories/")
        self.assertEqual(broken.status_code, 200, broken.data)
        self.assertEqual(broken.data["data"]["count"], 3)

        preview = client.post(
            "/api/ccp4i2/projects/rebase_root/",
            {
                "old_root": self.old_root,
                "new_root": str(self.new_root),
                "dry_run": True,
            },
            format="json",
        )
        self.assertEqual(preview.status_code, 200, preview.data)
        self.assertEqual(len(preview.data["data"]["matched"]), 3)

        applied = client.post(
            "/api/ccp4i2/projects/rebase_root/",
            {
                "old_root": self.old_root,
                "new_root": str(self.new_root),
                "update_preference": False,
            },
            format="json",
        )
        self.assertEqual(applied.status_code, 200, applied.data)
        self.assertEqual(len(applied.data["data"]["relocated"]), 3)

        after = client.get("/api/ccp4i2/projects/missing_directories/")
        self.assertEqual(after.data["data"]["count"], 0)
