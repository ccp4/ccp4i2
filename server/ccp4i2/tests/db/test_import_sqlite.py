"""Tests for legacy CCP4i2 SQLite migration: structural analysis + plan + copy.

These are fast, need no CCP4 binaries and no external data — they synthesise a
minimal legacy sqlite database and a temp project tree on the fly.
"""

import sqlite3
import tempfile
import uuid
from pathlib import Path

from django.conf import settings
from django.test import TestCase, override_settings

from ...db.import_sqlite import (
    SQLiteImporter,
    SQLiteValidator,
    StructuralIssuesError,
)
from ...db import legacy_structure as ls
from ...db.models import Project


# ---------------------------------------------------------------------------
# Fixture helpers
# ---------------------------------------------------------------------------

LEGACY_SCHEMA = """
CREATE TABLE Projects (
    ProjectID TEXT PRIMARY KEY, ProjectName TEXT, ProjectCreated REAL,
    UserID TEXT, ParentProjectID TEXT, ProjectDirectory TEXT,
    LastJobNumber INTEGER, FollowFromJobID TEXT, I1ProjectName TEXT,
    I1ProjectDirectory TEXT, LastAccess REAL
);
CREATE TABLE Jobs (
    JobID TEXT PRIMARY KEY, JobNumber TEXT, ProjectID TEXT, ParentJobID TEXT,
    PreceedingJobID TEXT, Status INTEGER, TaskName TEXT, TaskVersion TEXT,
    JobTitle TEXT, Evaluation INTEGER, UserId TEXT, UserAgent TEXT,
    CreationTime REAL, FinishTime REAL, processId INTEGER
);
CREATE TABLE Files (
    FileID TEXT PRIMARY KEY, Filename TEXT, Annotation TEXT, FiletypeID INTEGER,
    FileSubType INTEGER, FileContent INTEGER, FilePath TEXT, JobID TEXT,
    JobParamName TEXT, PathFlag INTEGER
);
CREATE TABLE FileUses (FileID TEXT, JobID TEXT, RoleID INTEGER, JobParamName TEXT);
CREATE TABLE ImportFiles (
    ImportId TEXT PRIMARY KEY, FileID TEXT, SourceFilename TEXT,
    SourceFileID TEXT, ExportFileID TEXT, Annotation TEXT, CreationTime REAL,
    LastModifiedTime REAL, Checksum TEXT, ImportNumber INTEGER, Reference TEXT
);
CREATE TABLE ExportFiles (
    ExportId TEXT PRIMARY KEY, FileID TEXT, ExportFilename TEXT,
    Annotation TEXT, CreationTime REAL
);
CREATE TABLE XData (XDataID TEXT, JobID TEXT, XDataClass TEXT, XDataXml TEXT);
CREATE TABLE JobKeyValues (JobID TEXT, KeyTypeID INTEGER, Value REAL);
CREATE TABLE JobKeyCharValues (JobID TEXT, KeyTypeID INTEGER, Value TEXT);
CREATE TABLE Tags (TagID INTEGER PRIMARY KEY, ParentTagID INTEGER, Text TEXT);
CREATE TABLE ProjectTags (TagID INTEGER, ProjectID TEXT);
CREATE TABLE Comments (CommentID TEXT, JobID TEXT, Comment TEXT);
CREATE TABLE ProjectComments (
    ProjectCommentID TEXT PRIMARY KEY, ProjectID TEXT, UserID TEXT,
    TimeOfComment REAL, Comment TEXT
);
CREATE TABLE ServerJobs (
    JobId TEXT, ServerProcessId INTEGER, Machine TEXT, Username TEXT,
    Mechanism TEXT, RemotePath TEXT, CustomCodeFile TEXT, Validate TEXT,
    KeyFilename TEXT, ServerGroup TEXT
);
CREATE TABLE FileTypes (FileTypeID INTEGER PRIMARY KEY, FileTypeName TEXT, FileTypeDescription TEXT);
CREATE TABLE KeyTypes (KeyTypeID INTEGER PRIMARY KEY, KeyTypeName TEXT, KeyTypeDescription TEXT);
CREATE TABLE FileRoles (RoleID INTEGER PRIMARY KEY, RoleText TEXT);
CREATE TABLE ProjectExports (ProjectExportId TEXT, ProjectID TEXT, ProjectExportTime REAL);
CREATE TABLE ProjectImports (ProjectImportId TEXT, ProjectID TEXT, ProjectImportTime REAL);
"""


def _hex():
    return uuid.uuid4().hex


def make_legacy_db(db_path, projects, jobs=None, files=None):
    """Create a minimal legacy sqlite DB.

    projects: list of dicts {id, name, directory, i1_directory?, parent?}
    jobs: optional list of {id, number, project_id, task}
    files: optional list of {filename, job_id, pathflag?} — DB rows only; the
        test decides whether the file exists on disk.
    """
    conn = sqlite3.connect(str(db_path))
    conn.executescript(LEGACY_SCHEMA)
    for p in projects:
        conn.execute(
            "INSERT INTO Projects (ProjectID, ProjectName, ProjectCreated, UserID, "
            "ParentProjectID, ProjectDirectory, LastJobNumber, I1ProjectName, "
            "I1ProjectDirectory, LastAccess) VALUES (?,?,?,?,?,?,?,?,?,?)",
            (p["id"], p["name"], 1700000000.0, "tester", p.get("parent"),
             p["directory"], 1, "", p.get("i1_directory", ""), 1700000000.0),
        )
    for j in (jobs or []):
        conn.execute(
            "INSERT INTO Jobs (JobID, JobNumber, ProjectID, Status, TaskName, "
            "CreationTime) VALUES (?,?,?,?,?,?)",
            (j["id"], j["number"], j["project_id"], 6, j.get("task", "refmac"),
             1700000000.0),
        )
    for f in (files or []):
        conn.execute(
            "INSERT INTO Files (FileID, Filename, JobID, PathFlag) VALUES (?,?,?,?)",
            (_hex(), f["filename"], f["job_id"], f.get("pathflag", 1)),
        )
    conn.commit()
    conn.close()


def make_project_tree(base, name, with_jobs=True):
    d = Path(base) / name
    if with_jobs:
        jd = d / "CCP4_JOBS" / "job_1"
        jd.mkdir(parents=True, exist_ok=True)
        (jd / "XYZOUT.pdb").write_text("ATOM\n" * 10)
    else:
        d.mkdir(parents=True, exist_ok=True)
    (d / "CCP4_IMPORTED_FILES").mkdir(parents=True, exist_ok=True)
    return str(d)


# ---------------------------------------------------------------------------
# Pure structural-analysis tests (no Django DB)
# ---------------------------------------------------------------------------

class StructureAnalysisTests(TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.dest_root = self.tmp / "CCP4X_PROJECTS"

    def test_nested_project_flagged_and_planned_as_copy(self):
        a = make_project_tree(self.tmp, "alpha")
        # beta nested inside alpha's jobs dir
        b = str(Path(a) / "CCP4_JOBS" / "job_1" / "beta")
        (Path(b) / "CCP4_JOBS").mkdir(parents=True, exist_ok=True)
        projects = [
            {"id": "a" * 32, "name": "alpha", "directory": a, "i1_directory": ""},
            {"id": "b" * 32, "name": "beta", "directory": b, "i1_directory": ""},
        ]
        issues, plan = ls.analyse_structure(
            projects, copy_files=False, dest_root=self.dest_root)

        # beta is copied out even though copy_files=False (mixed mode)
        beta = next(e for e in plan if e["name"] == "beta")
        alpha = next(e for e in plan if e["name"] == "alpha")
        self.assertEqual(beta["mode"], "copy")
        self.assertEqual(beta["reason"], "nested")
        self.assertEqual(alpha["mode"], "in_place")
        # nesting is auto-resolved -> not blocking-unresolved
        self.assertEqual(ls.blocking_unresolved(issues), 0)
        self.assertTrue(any(i["type"] == ls.NESTED_PROJECT and i["resolution"]
                            for i in issues))

    def test_dest_collision_deduped_by_rename(self):
        a1 = make_project_tree(self.tmp / "one", "proj")
        a2 = make_project_tree(self.tmp / "two", "proj")  # same basename
        projects = [
            {"id": "1" * 32, "name": "proj", "directory": a1, "i1_directory": ""},
            {"id": "2" * 32, "name": "proj", "directory": a2, "i1_directory": ""},
        ]
        _issues, plan = ls.analyse_structure(
            projects, copy_files=True, dest_root=self.dest_root)
        dests = [e["dest"] for e in plan]
        self.assertEqual(len(set(dests)), 2, "dests must be distinct")
        self.assertTrue(any(e["renamed_to"] for e in plan))

    def test_not_project_root_warning(self):
        d = make_project_tree(self.tmp, "noroot", with_jobs=False)  # no CCP4_JOBS
        projects = [{"id": "n" * 32, "name": "noroot", "directory": d,
                     "i1_directory": ""}]
        issues, _plan = ls.analyse_structure(
            projects, copy_files=False, dest_root=self.dest_root)
        self.assertTrue(any(i["type"] == ls.NOT_PROJECT_ROOT for i in issues))

    def test_nested_excludes(self):
        a = make_project_tree(self.tmp, "alpha")
        b = str(Path(a) / "CCP4_JOBS" / "job_1" / "beta")
        (Path(b) / "CCP4_JOBS").mkdir(parents=True, exist_ok=True)
        projects = [
            {"id": "a" * 32, "name": "alpha", "directory": a, "i1_directory": ""},
            {"id": "b" * 32, "name": "beta", "directory": b, "i1_directory": ""},
        ]
        _issues, plan = ls.analyse_structure(
            projects, copy_files=True, dest_root=self.dest_root)
        resolved = {p["id"]: ls._resolve(p["directory"]) for p in projects}
        excludes = ls.nested_excludes_for(plan, resolved)
        self.assertTrue(any("beta" in str(d) for d in excludes["a" * 32]))


# ---------------------------------------------------------------------------
# Full validator + importer tests (Django DB)
# ---------------------------------------------------------------------------

@override_settings()
class ImporterCopyTests(TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.legacy_root = self.tmp / "legacy"
        self.legacy_root.mkdir()
        self.dest_root = self.tmp / "CCP4X_PROJECTS"
        self.db_path = self.tmp / "database.sqlite"

    def _validator(self, **kw):
        return SQLiteValidator(db_path=self.db_path, dest_root=self.dest_root, **kw)

    def _importer(self, **kw):
        kw.setdefault("dest_root", self.dest_root)
        kw.setdefault("continue_on_error", True)
        return SQLiteImporter(db_path=self.db_path, **kw)

    def test_validate_reports_plan_and_summary(self):
        d = make_project_tree(self.legacy_root, "gamma")
        pid = _hex()
        make_legacy_db(self.db_path,
                       [{"id": pid, "name": "gamma", "directory": d}],
                       [{"id": _hex(), "number": "1", "project_id": pid}])
        report = self._validator(copy_files=True).run()
        self.assertIn("structure", report)
        self.assertIn("plan", report)
        self.assertEqual(len(report["plan"]), 1)
        self.assertEqual(report["plan"][0]["mode"], "copy")
        self.assertEqual(report["summary"]["plan_summary"]["copy"], 1)

    def test_missing_files_split_by_preservation_contract(self):
        # A project with a top-level job (1) and a sub-job (1.1). Give each a
        # Files row whose file is NOT on disk. The top-level miss is an
        # in-contract violation (gates ok); the sub-job miss is out of contract.
        d = make_project_tree(self.legacy_root, "kappa")
        # sub-job dir must exist so only the FILE is missing, not the dir
        (Path(d) / "CCP4_JOBS" / "job_1" / "job_1").mkdir(parents=True, exist_ok=True)
        pid = _hex()
        top_job, sub_job = _hex(), _hex()
        make_legacy_db(
            self.db_path,
            [{"id": pid, "name": "kappa", "directory": d}],
            [{"id": top_job, "number": "1", "project_id": pid},
             {"id": sub_job, "number": "1.1", "project_id": pid}],
            files=[
                {"filename": "TOPLEVEL_RESULT.pdb", "job_id": top_job},
                {"filename": "sub_intermediate.mtz", "job_id": sub_job},
            ],
        )
        s = self._validator().run()["summary"]
        self.assertEqual(s["top_level_files_missing"], 1)
        self.assertEqual(s["sub_job_files_missing"], 1)
        # In-contract violation makes the run not-ok; sub-job miss alone would not.
        self.assertFalse(s["ok"])
        # But it is not a *blocking* structural issue — migration can proceed.
        self.assertEqual(s["blocking_issues"], 0)

    def test_copy_mode_copies_tree_and_repoints(self):
        d = make_project_tree(self.legacy_root, "delta")
        pid = _hex()
        make_legacy_db(self.db_path,
                       [{"id": pid, "name": "delta", "directory": d}],
                       [{"id": _hex(), "number": "1", "project_id": pid}])

        self._importer(copy_files=True).run()

        proj = Project.objects.get(name="delta")
        expected_dest = self.dest_root / "delta"
        self.assertEqual(Path(proj.directory), expected_dest)
        self.assertTrue((expected_dest / "CCP4_JOBS" / "job_1" / "XYZOUT.pdb").is_file())

    def test_dry_run_copies_nothing(self):
        d = make_project_tree(self.legacy_root, "epsilon")
        pid = _hex()
        make_legacy_db(self.db_path,
                       [{"id": pid, "name": "epsilon", "directory": d}])
        result = self._importer(copy_files=True, dry_run=True).run()
        self.assertFalse((self.dest_root / "epsilon").exists())
        self.assertEqual(Project.objects.filter(name="epsilon").count(), 0)
        # stats still report the intended copy
        self.assertEqual(result["stats"].get("projects_copied", 0), 1)

    def test_nested_mixed_mode_copies_only_inner(self):
        a = make_project_tree(self.legacy_root, "outer")
        b = str(Path(a) / "CCP4_JOBS" / "job_1" / "inner")
        (Path(b) / "CCP4_JOBS" / "job_1").mkdir(parents=True, exist_ok=True)
        (Path(b) / "CCP4_JOBS" / "job_1" / "IN.pdb").write_text("x")
        (Path(b) / "CCP4_IMPORTED_FILES").mkdir(parents=True, exist_ok=True)
        outer_id, inner_id = _hex(), _hex()
        make_legacy_db(self.db_path, [
            {"id": outer_id, "name": "outer", "directory": a},
            {"id": inner_id, "name": "inner", "directory": b},
        ])
        # schlurp overall (copy_files=False); inner must still be copied out
        self._importer(copy_files=False).run()

        outer = Project.objects.get(name="outer")
        inner = Project.objects.get(name="inner")
        # in place: directory is the (canonicalised) legacy dir
        self.assertEqual(Path(outer.directory).resolve(), Path(a).resolve())
        # copied out to dest_root even though the run was schlurp overall
        self.assertEqual(Path(inner.directory).resolve(),
                         (self.dest_root / "inner").resolve())
        self.assertTrue((self.dest_root / "inner" / "CCP4_JOBS" / "job_1" / "IN.pdb").is_file())

    def test_blocking_issue_refused_without_ack(self):
        # An unreadable source in copy mode is blocking-unresolved.
        d = make_project_tree(self.legacy_root, "locked")
        import os
        os.chmod(Path(d) / "CCP4_JOBS", 0o000)
        try:
            pid = _hex()
            make_legacy_db(self.db_path,
                           [{"id": pid, "name": "locked", "directory": d}])
            with self.assertRaises(StructuralIssuesError):
                self._importer(copy_files=True).run()
            # with ack it proceeds (copy of the unreadable subdir is skipped/errs
            # but the run completes)
            self._importer(copy_files=True,
                           allow_structural_issues=True).run()
        finally:
            os.chmod(Path(d) / "CCP4_JOBS", 0o755)
