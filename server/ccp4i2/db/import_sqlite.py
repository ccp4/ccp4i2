# Copyright (C) 2026 Newcastle University
#
# This file is part of CCP4i2.
#
# CCP4i2 is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License version 3,
# modified in accordance with the provisions of the license to address
# the requirements of UK law.
#
# See https://www.ccp4.ac.uk/ccp4license.php for details.
"""
Import and validate legacy CCP4i2 data from a SQLite database file.

This module reads the legacy CCP4i2 SQLite database (e.g. ~/.CCP4I2/db/database.sqlite)
as ground truth and can either validate or import the data into the Django schema.

Classes:
    SQLiteValidator - Validate the legacy database against the filesystem
    SQLiteImporter  - Import all data into Django models

Both are usable from management commands and API endpoints.
"""

import logging
import sqlite3
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from uuid import UUID

from django.db import transaction
from django.utils import timezone

from .models import (
    File,
    FileExport,
    FileImport,
    FileType,
    FileUse,
    Job,
    JobCharValue,
    JobFloatValue,
    JobValueKey,
    Project,
    ProjectExport,
    ProjectGroup,
    ProjectGroupMembership,
    ProjectImport,
    ProjectTag,
    ServerJob,
    XData,
)

logger = logging.getLogger(f"ccp4i2:{__name__}")


def convert_uuid(hex_str):
    """Convert 32-char UUID hex string to UUID object."""
    if not hex_str:
        return None
    try:
        return UUID(hex_str)
    except (ValueError, TypeError):
        return None


def convert_timestamp(epoch_float):
    """Convert epoch float to timezone-aware datetime."""
    if epoch_float is None:
        return None
    try:
        dt = datetime.fromtimestamp(float(epoch_float))
        return timezone.make_aware(dt, timezone.get_default_timezone())
    except (ValueError, TypeError, OSError):
        return None


def dict_factory(cursor, row):
    """sqlite3 row factory that returns dicts with lowercase keys."""
    return {col[0].lower(): row[i] for i, col in enumerate(cursor.description)}


class SQLiteValidator:
    """Validate a legacy CCP4i2 SQLite database against the filesystem.

    Checks that project directories, job directories, output files, and
    imported files all exist on disk. Also checks for referential integrity
    (orphan jobs, files pointing to missing jobs, etc.) and data quality
    (empty task names, null timestamps, invalid UUIDs).

    This is a read-only operation — it never writes to Django or the SQLite.
    """

    def __init__(self, db_path, remap_dirs=None, verbose=False, log_fn=None):
        self.db_path = Path(db_path)
        self.remap_dirs = remap_dirs
        self.verbose = verbose
        self.log_fn = log_fn or (lambda msg: logger.info(msg))

    def remap_directory(self, directory):
        if not directory or not self.remap_dirs:
            return directory
        from_path, to_path = self.remap_dirs
        if directory.startswith(from_path):
            return directory.replace(from_path, to_path, 1)
        return directory

    def run(self):
        """Execute all validations. Returns a structured report dict."""
        if not self.db_path.exists():
            raise FileNotFoundError(f"Database not found: {self.db_path}")

        conn = sqlite3.connect(str(self.db_path))
        conn.row_factory = dict_factory
        cur = conn.cursor()

        report = {
            "source": str(self.db_path),
            "counts": self._table_counts(cur),
            "projects": self._validate_projects(cur),
            "jobs": self._validate_jobs(cur),
            "files": self._validate_files(cur),
            "imported_files": self._validate_imported_files(cur),
            "integrity": self._validate_integrity(cur),
            "data_quality": self._validate_data_quality(cur),
        }

        conn.close()

        # Compute a top-level summary
        report["summary"] = self._summarise(report)
        return report

    # ------------------------------------------------------------------
    # Table counts
    # ------------------------------------------------------------------

    def _table_counts(self, cur):
        counts = {}
        for table in [
            "Projects", "Jobs", "Files", "FileUses", "ImportFiles",
            "ExportFiles", "XData", "JobKeyValues", "JobKeyCharValues",
            "Tags", "ProjectTags", "Comments", "ProjectComments",
            "ServerJobs", "FileTypes", "KeyTypes", "FileRoles",
            "ProjectExports", "ProjectImports",
        ]:
            try:
                cur.execute(f"SELECT COUNT(*) as cnt FROM [{table}]")
                counts[table] = cur.fetchone()["cnt"]
            except sqlite3.OperationalError:
                counts[table] = None  # table missing
        return counts

    # ------------------------------------------------------------------
    # Project validation
    # ------------------------------------------------------------------

    def _validate_projects(self, cur):
        cur.execute(
            "SELECT ProjectID, ProjectName, ProjectDirectory FROM Projects"
        )
        results = {"total": 0, "dir_exists": 0, "dir_missing": [], "dir_empty": []}
        for row in cur.fetchall():
            results["total"] += 1
            directory = self.remap_directory(row["projectdirectory"] or "")
            if not directory:
                results["dir_empty"].append(row["projectname"])
                continue
            p = Path(directory)
            if p.is_dir():
                results["dir_exists"] += 1
            else:
                results["dir_missing"].append({
                    "project": row["projectname"],
                    "directory": directory,
                })

        self._log_section("Projects", results["total"],
                          results["dir_exists"], results["dir_missing"])
        return results

    # ------------------------------------------------------------------
    # Job validation
    # ------------------------------------------------------------------

    def _validate_jobs(self, cur):
        # Build project directory lookup
        cur.execute("SELECT ProjectID, ProjectDirectory FROM Projects")
        project_dirs = {}
        for row in cur.fetchall():
            d = self.remap_directory(row["projectdirectory"] or "")
            project_dirs[row["projectid"]] = d

        cur.execute(
            "SELECT JobID, JobNumber, ProjectID, ParentJobID, Status, TaskName "
            "FROM Jobs"
        )
        results = {
            "total": 0, "dir_exists": 0,
            "dir_missing": [], "dir_missing_count": 0,
            "project_missing": 0,
        }
        for row in cur.fetchall():
            results["total"] += 1
            proj_dir = project_dirs.get(row["projectid"])
            if not proj_dir:
                results["project_missing"] += 1
                continue

            # Build job directory: project_dir/CCP4_JOBS/job_N/job_M/...
            number = row["jobnumber"]
            path_elements = [f"job_{e}" for e in number.split(".")]
            job_dir = Path(proj_dir) / "CCP4_JOBS" / Path(*path_elements)

            if job_dir.is_dir():
                results["dir_exists"] += 1
            else:
                results["dir_missing_count"] += 1
                if self.verbose or len(results["dir_missing"]) < 50:
                    results["dir_missing"].append({
                        "job_number": number,
                        "task_name": row["taskname"],
                        "expected_dir": str(job_dir),
                    })

        self._log_section("Jobs", results["total"],
                          results["dir_exists"], results["dir_missing"],
                          missing_count=results["dir_missing_count"])
        return results

    # ------------------------------------------------------------------
    # File validation
    # ------------------------------------------------------------------

    def _validate_files(self, cur):
        # Build project dir + job number lookups
        cur.execute("SELECT ProjectID, ProjectDirectory FROM Projects")
        project_dirs = {}
        for row in cur.fetchall():
            project_dirs[row["projectid"]] = self.remap_directory(
                row["projectdirectory"] or ""
            )

        cur.execute("SELECT JobID, JobNumber, ProjectID FROM Jobs")
        job_info = {}
        for row in cur.fetchall():
            job_info[row["jobid"]] = {
                "number": row["jobnumber"],
                "project_id": row["projectid"],
            }

        cur.execute(
            "SELECT FileID, Filename, JobID, PathFlag FROM Files"
        )
        results = {
            "total": 0, "exists": 0,
            "missing": [], "missing_count": 0,
            "orphan_job": 0,
        }
        for row in cur.fetchall():
            results["total"] += 1
            job_id = row["jobid"]
            ji = job_info.get(job_id)
            if not ji:
                results["orphan_job"] += 1
                continue

            proj_dir = project_dirs.get(ji["project_id"], "")
            if not proj_dir:
                continue

            pathflag = row["pathflag"]
            filename = row["filename"] or ""
            if pathflag == 2:
                # Import directory
                file_path = Path(proj_dir) / "CCP4_IMPORTED_FILES" / filename
            else:
                # Job directory
                path_elements = [f"job_{e}" for e in ji["number"].split(".")]
                file_path = Path(proj_dir) / "CCP4_JOBS" / Path(*path_elements) / filename

            if file_path.is_file():
                results["exists"] += 1
            else:
                results["missing_count"] += 1
                if self.verbose or len(results["missing"]) < 50:
                    results["missing"].append({
                        "filename": filename,
                        "expected_path": str(file_path),
                    })

        self._log_section("Files", results["total"],
                          results["exists"], results["missing"],
                          missing_count=results["missing_count"])
        return results

    # ------------------------------------------------------------------
    # Imported file validation (source files)
    # ------------------------------------------------------------------

    def _validate_imported_files(self, cur):
        cur.execute(
            "SELECT ImportId, FileID, SourceFilename FROM ImportFiles"
        )
        results = {
            "total": 0, "source_exists": 0,
            "source_missing": [], "source_missing_count": 0,
        }
        for row in cur.fetchall():
            results["total"] += 1
            src = row["sourcefilename"] or ""
            if not src:
                continue
            src = self.remap_directory(src)
            if Path(src).is_file():
                results["source_exists"] += 1
            else:
                results["source_missing_count"] += 1
                if self.verbose or len(results["source_missing"]) < 50:
                    results["source_missing"].append(src)

        self._log_section("ImportFiles (source)", results["total"],
                          results["source_exists"], results["source_missing"],
                          missing_count=results["source_missing_count"])
        return results

    # ------------------------------------------------------------------
    # Referential integrity
    # ------------------------------------------------------------------

    def _validate_integrity(self, cur):
        issues = []

        # Jobs referencing non-existent projects
        cur.execute(
            "SELECT COUNT(*) as cnt FROM Jobs j "
            "LEFT JOIN Projects p ON j.ProjectID = p.ProjectID "
            "WHERE p.ProjectID IS NULL"
        )
        orphan_jobs = cur.fetchone()["cnt"]
        if orphan_jobs:
            issues.append(f"{orphan_jobs} jobs reference non-existent projects")

        # Jobs referencing non-existent parent jobs
        cur.execute(
            "SELECT COUNT(*) as cnt FROM Jobs j "
            "LEFT JOIN Jobs p ON j.ParentJobID = p.JobID "
            "WHERE j.ParentJobID IS NOT NULL AND p.JobID IS NULL"
        )
        orphan_parents = cur.fetchone()["cnt"]
        if orphan_parents:
            issues.append(f"{orphan_parents} jobs reference non-existent parent jobs")

        # Files referencing non-existent jobs
        cur.execute(
            "SELECT COUNT(*) as cnt FROM Files f "
            "LEFT JOIN Jobs j ON f.JobID = j.JobID "
            "WHERE f.JobID IS NOT NULL AND j.JobID IS NULL"
        )
        orphan_files = cur.fetchone()["cnt"]
        if orphan_files:
            issues.append(f"{orphan_files} files reference non-existent jobs")

        # FileUses referencing non-existent files or jobs
        cur.execute(
            "SELECT COUNT(*) as cnt FROM FileUses fu "
            "LEFT JOIN Files f ON fu.FileID = f.FileID "
            "WHERE f.FileID IS NULL"
        )
        orphan_fileuses_file = cur.fetchone()["cnt"]
        if orphan_fileuses_file:
            issues.append(f"{orphan_fileuses_file} file uses reference non-existent files")

        cur.execute(
            "SELECT COUNT(*) as cnt FROM FileUses fu "
            "LEFT JOIN Jobs j ON fu.JobID = j.JobID "
            "WHERE j.JobID IS NULL"
        )
        orphan_fileuses_job = cur.fetchone()["cnt"]
        if orphan_fileuses_job:
            issues.append(f"{orphan_fileuses_job} file uses reference non-existent jobs")

        # ImportFiles referencing non-existent files
        cur.execute(
            "SELECT COUNT(*) as cnt FROM ImportFiles i "
            "LEFT JOIN Files f ON i.FileID = f.FileID "
            "WHERE f.FileID IS NULL"
        )
        orphan_imports = cur.fetchone()["cnt"]
        if orphan_imports:
            issues.append(f"{orphan_imports} import records reference non-existent files")

        # Duplicate job numbers within a project
        cur.execute(
            "SELECT ProjectID, JobNumber, COUNT(*) as cnt "
            "FROM Jobs GROUP BY ProjectID, JobNumber HAVING cnt > 1"
        )
        dupes = cur.fetchall()
        if dupes:
            issues.append(f"{len(dupes)} duplicate job numbers within projects")

        self.log_fn(f"  Integrity: {'OK' if not issues else f'{len(issues)} issues'}")
        for issue in issues:
            self.log_fn(f"    - {issue}")

        return {"ok": len(issues) == 0, "issues": issues}

    # ------------------------------------------------------------------
    # Data quality
    # ------------------------------------------------------------------

    def _validate_data_quality(self, cur):
        issues = []

        # Invalid UUIDs (not 32 hex chars)
        for table, col in [
            ("Projects", "ProjectID"), ("Jobs", "JobID"), ("Files", "FileID"),
        ]:
            cur.execute(
                f"SELECT COUNT(*) as cnt FROM [{table}] "
                f"WHERE length([{col}]) != 32"
            )
            bad = cur.fetchone()["cnt"]
            if bad:
                issues.append(f"{bad} rows in {table} have non-standard UUID length")

        # Jobs with empty task names
        cur.execute(
            "SELECT COUNT(*) as cnt FROM Jobs "
            "WHERE TaskName IS NULL OR TaskName = ''"
        )
        empty_tasks = cur.fetchone()["cnt"]
        if empty_tasks:
            issues.append(f"{empty_tasks} jobs have empty task names")

        # Projects with NULL creation time
        cur.execute(
            "SELECT COUNT(*) as cnt FROM Projects WHERE ProjectCreated IS NULL"
        )
        null_ts = cur.fetchone()["cnt"]
        if null_ts:
            issues.append(f"{null_ts} projects have NULL creation time")

        # Jobs with NULL creation time
        cur.execute(
            "SELECT COUNT(*) as cnt FROM Jobs WHERE CreationTime IS NULL"
        )
        null_job_ts = cur.fetchone()["cnt"]
        if null_job_ts:
            issues.append(f"{null_job_ts} jobs have NULL creation time")

        # Jobs with unknown status values
        cur.execute(
            "SELECT COUNT(*) as cnt FROM Jobs "
            "WHERE Status NOT IN (0,1,2,3,4,5,6,7,8,9,10)"
        )
        bad_status = cur.fetchone()["cnt"]
        if bad_status:
            issues.append(f"{bad_status} jobs have unrecognised status values")

        # Files with unknown file type IDs (not in FileTypes table)
        cur.execute(
            "SELECT COUNT(*) as cnt FROM Files f "
            "LEFT JOIN FileTypes ft ON f.FiletypeID = ft.FileTypeID "
            "WHERE ft.FileTypeID IS NULL"
        )
        unknown_ft = cur.fetchone()["cnt"]
        if unknown_ft:
            issues.append(f"{unknown_ft} files reference unknown file types")

        # Projects with empty directories
        cur.execute(
            "SELECT COUNT(*) as cnt FROM Projects "
            "WHERE ProjectDirectory IS NULL OR ProjectDirectory = ''"
        )
        empty_dirs = cur.fetchone()["cnt"]
        if empty_dirs:
            issues.append(f"{empty_dirs} projects have empty directories")

        self.log_fn(f"  Data quality: {'OK' if not issues else f'{len(issues)} issues'}")
        for issue in issues:
            self.log_fn(f"    - {issue}")

        return {"ok": len(issues) == 0, "issues": issues}

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _log_section(self, label, total, found, missing_list,
                     missing_count=None):
        actual_missing = missing_count if missing_count is not None else len(missing_list)
        self.log_fn(
            f"  {label}: {total} total, {found} found on disk, "
            f"{actual_missing} missing"
        )
        if self.verbose:
            for item in missing_list:
                if isinstance(item, dict):
                    detail = item.get("expected_path") or item.get("expected_dir") or item.get("directory") or str(item)
                    self.log_fn(f"    MISSING: {detail}")
                else:
                    self.log_fn(f"    MISSING: {item}")

    def _summarise(self, report):
        projects = report["projects"]
        jobs = report["jobs"]
        files = report["files"]
        imports = report["imported_files"]
        integrity = report["integrity"]
        quality = report["data_quality"]

        all_ok = (
            not projects["dir_missing"]
            and not projects["dir_empty"]
            and jobs["dir_missing_count"] == 0
            and files["missing_count"] == 0
            and integrity["ok"]
            and quality["ok"]
        )

        return {
            "ok": all_ok,
            "projects_on_disk": f"{projects['dir_exists']}/{projects['total']}",
            "jobs_on_disk": f"{jobs['dir_exists']}/{jobs['total']}",
            "files_on_disk": f"{files['exists']}/{files['total']}",
            "import_sources_on_disk": f"{imports['source_exists']}/{imports['total']}",
            "integrity_issues": len(integrity["issues"]),
            "data_quality_issues": len(quality["issues"]),
        }


class SQLiteImporter:
    """Import legacy CCP4i2 SQLite database into Django models."""

    def __init__(self, db_path, remap_dirs=None, dry_run=False,
                 continue_on_error=False, verbose=False, log_fn=None):
        self.db_path = Path(db_path)
        self.remap_dirs = remap_dirs  # (from_path, to_path) or None
        self.dry_run = dry_run
        self.continue_on_error = continue_on_error
        self.verbose = verbose
        self.log_fn = log_fn or (lambda msg: logger.info(msg))

        # Mapping caches: legacy ID -> Django object
        self.project_map = {}
        self.job_map = {}
        self.file_map = {}
        self.filetype_map = {}
        self.keytype_map = {}
        self.tag_map = {}
        self.filerole_map = {}

        # Parent relationships deferred for group creation
        self.parent_relationships = []  # [(child_legacy_id, parent_legacy_id)]

        self.stats = defaultdict(int)
        self.errors = []

    def log(self, message):
        if self.verbose:
            self.log_fn(message)

    def log_count(self, entity, count, extra=""):
        extra_msg = f" ({extra})" if extra else ""
        self.log_fn(f"  {entity}: {count}{extra_msg}")

    def remap_directory(self, directory):
        if not directory or not self.remap_dirs:
            return directory
        from_path, to_path = self.remap_dirs
        if directory.startswith(from_path):
            return directory.replace(from_path, to_path, 1)
        return directory

    def get_unique_project_name(self, name):
        if not Project.objects.filter(name=name).exists():
            return name
        for suffix in ["_legacy"] + [f"_legacy_{i}" for i in range(2, 1000)]:
            candidate = f"{name}{suffix}"
            if not Project.objects.filter(name=candidate).exists():
                self.stats["projects_renamed"] += 1
                return candidate
        raise RuntimeError(f"Unable to find unique name for project: {name}")

    def _record_error(self, phase, record_id, error):
        msg = f"{phase} {record_id}: {error}"
        self.errors.append(msg)
        logger.error(msg)
        if not self.continue_on_error:
            raise

    def run(self):
        """Execute the full import. Returns a summary dict."""
        if not self.db_path.exists():
            raise FileNotFoundError(f"Database not found: {self.db_path}")

        conn = sqlite3.connect(str(self.db_path))
        conn.row_factory = dict_factory

        try:
            if self.dry_run:
                self.log_fn("[DRY RUN - no changes will be committed]")

            with transaction.atomic():
                self._import_all(conn)
                if self.dry_run:
                    transaction.set_rollback(True)
        finally:
            conn.close()

        return self.summary()

    def _import_all(self, conn):
        """Run all import phases in order."""
        cur = conn.cursor()

        self.log_fn("Phase 1: Lookup tables")
        self._import_filetypes(cur)
        self._import_keytypes(cur)
        self._import_fileroles(cur)
        self._import_tags(cur)

        self.log_fn("Phase 2: Projects")
        self._import_projects(cur)
        self._import_projecttags(cur)

        self.log_fn("Phase 3: Project groups (parent relationships)")
        self._create_project_groups()

        self.log_fn("Phase 4: Jobs")
        self._import_jobs(cur)
        self._import_serverjobs(cur)
        self._import_jobkeyvalues(cur)
        self._import_jobkeycharvalues(cur)
        self._import_xdata(cur)

        self.log_fn("Phase 5: Files")
        self._import_files(cur)
        self._import_importfiles(cur)
        self._import_exportfiles(cur)
        self._import_fileuses(cur)

        self.log_fn("Phase 6: Project history")
        self._import_projectexports(cur)
        self._import_projectimports(cur)

        self.log_fn("Phase 7: Comments")
        self._import_project_comments(cur)

    # ------------------------------------------------------------------
    # Lookup tables
    # ------------------------------------------------------------------

    def _import_filetypes(self, cur):
        cur.execute("SELECT FileTypeID, FileTypeName, FileTypeDescription FROM FileTypes")
        for row in cur.fetchall():
            ft, created = FileType.objects.get_or_create(
                name=row["filetypename"],
                defaults={"description": row["filetypedescription"] or ""},
            )
            self.filetype_map[row["filetypeid"]] = ft
            if created:
                self.stats["filetypes"] += 1
        self.log_count("FileTypes", self.stats["filetypes"])

    def _import_keytypes(self, cur):
        cur.execute("SELECT KeyTypeID, KeyTypeName, KeyTypeDescription FROM KeyTypes")
        for row in cur.fetchall():
            kt, created = JobValueKey.objects.get_or_create(
                name=row["keytypename"],
                defaults={"description": row["keytypedescription"] or ""},
            )
            self.keytype_map[row["keytypeid"]] = kt
            if created:
                self.stats["keytypes"] += 1
        self.log_count("JobValueKeys", self.stats["keytypes"])

    def _import_fileroles(self, cur):
        cur.execute("SELECT RoleID, RoleText FROM FileRoles")
        for row in cur.fetchall():
            roletext = (row["roletext"] or "").upper()
            if roletext == "IN":
                self.filerole_map[row["roleid"]] = FileUse.Role.IN
            else:
                self.filerole_map[row["roleid"]] = FileUse.Role.OUT
        self.log_count("FileRoles mapped", len(self.filerole_map))

    def _import_tags(self, cur):
        cur.execute("SELECT TagID, ParentTagID, Text FROM Tags")
        rows = cur.fetchall()

        # First pass: create without parents
        tag_parent_map = {}
        for row in rows:
            tag, created = ProjectTag.objects.get_or_create(
                text=row["text"], parent=None,
            )
            self.tag_map[row["tagid"]] = tag
            if row["parenttagid"]:
                tag_parent_map[tag.id] = row["parenttagid"]
            if created:
                self.stats["tags"] += 1

        # Second pass: wire up parents
        for tag_id, legacy_parent_id in tag_parent_map.items():
            if legacy_parent_id in self.tag_map:
                tag = ProjectTag.objects.get(id=tag_id)
                tag.parent = self.tag_map[legacy_parent_id]
                tag.save()

        self.log_count("Tags", self.stats["tags"])

    # ------------------------------------------------------------------
    # Projects
    # ------------------------------------------------------------------

    def _import_projects(self, cur):
        cur.execute(
            "SELECT ProjectID, ProjectName, ProjectCreated, UserID, "
            "ParentProjectID, ProjectDirectory, LastJobNumber, "
            "FollowFromJobID, I1ProjectName, I1ProjectDirectory, LastAccess "
            "FROM Projects"
        )
        for row in cur.fetchall():
            try:
                pk = row["projectid"]
                directory = self.remap_directory(row["projectdirectory"] or "")

                # Skip if project already exists by directory
                if directory:
                    existing = Project.objects.filter(directory=directory).first()
                    if existing:
                        self.project_map[pk] = existing
                        self.stats["projects_skipped"] += 1
                        self.log(f"  Project exists (skipped): {existing.name}")
                        if row["parentprojectid"]:
                            self.parent_relationships.append((pk, row["parentprojectid"]))
                        continue

                name = self.get_unique_project_name(row["projectname"])
                project = Project.objects.create(
                    uuid=convert_uuid(pk) or None,
                    name=name,
                    description="",
                    directory=directory,
                    creation_time=convert_timestamp(row["projectcreated"]) or timezone.now(),
                    creation_user=row["userid"] or "legacy_import",
                    creation_host="legacy_import",
                    last_access=convert_timestamp(row["lastaccess"]) or timezone.now(),
                    last_job_number=int(row["lastjobnumber"] or 0),
                    i1_project_name=row["i1projectname"] or "",
                    i1_project_directory=row["i1projectdirectory"] or "",
                )
                self.project_map[pk] = project
                self.stats["projects"] += 1
                self.log(f"  Project: {name}")

                if row["parentprojectid"]:
                    self.parent_relationships.append((pk, row["parentprojectid"]))

            except Exception as e:
                if self.continue_on_error:
                    self._record_error("Project", pk, e)
                else:
                    raise

        parts = []
        if self.stats["projects_renamed"]:
            parts.append(f"{self.stats['projects_renamed']} renamed")
        if self.stats["projects_skipped"]:
            parts.append(f"{self.stats['projects_skipped']} skipped")
        self.log_count("Projects", self.stats["projects"], ", ".join(parts))

    def _import_projecttags(self, cur):
        cur.execute("SELECT TagID, ProjectID FROM ProjectTags")
        for row in cur.fetchall():
            pid = row["projectid"]
            tid = row["tagid"]
            if pid in self.project_map and tid in self.tag_map:
                self.project_map[pid].tags.add(self.tag_map[tid])
                self.stats["projecttags"] += 1
        self.log_count("Project-Tag associations", self.stats["projecttags"])

    def _create_project_groups(self):
        parent_children = defaultdict(list)
        for child_id, parent_id in self.parent_relationships:
            if child_id in self.project_map and parent_id in self.project_map:
                parent_children[parent_id].append(child_id)

        for parent_id, child_ids in parent_children.items():
            parent_project = self.project_map[parent_id]
            group_name = f"Legacy_{parent_project.name}_group"
            base_name = group_name
            counter = 1
            while ProjectGroup.objects.filter(name=group_name).exists():
                group_name = f"{base_name}_{counter}"
                counter += 1

            group = ProjectGroup.objects.create(
                name=group_name,
                type=ProjectGroup.GroupType.FRAGMENT_SET,
            )
            self.stats["groups"] += 1

            ProjectGroupMembership.objects.create(
                group=group, project=parent_project,
                type=ProjectGroupMembership.MembershipType.PARENT,
            )
            self.stats["memberships"] += 1

            for child_id in child_ids:
                ProjectGroupMembership.objects.create(
                    group=group, project=self.project_map[child_id],
                    type=ProjectGroupMembership.MembershipType.MEMBER,
                )
                self.stats["memberships"] += 1

        self.log_count("Groups", self.stats["groups"])
        self.log_count("Memberships", self.stats["memberships"])

    # ------------------------------------------------------------------
    # Jobs
    # ------------------------------------------------------------------

    def _import_jobs(self, cur):
        # Fetch sorted by job number so parents are created before children
        cur.execute(
            "SELECT JobID, JobNumber, CreationTime, FinishTime, Status, "
            "Evaluation, UserId, UserAgent, JobTitle, ProjectID, TaskName, "
            "TaskVersion, ParentJobID, PreceedingJobID, processId "
            "FROM Jobs ORDER BY length(JobNumber), JobNumber"
        )
        job_parent_map = {}

        for row in cur.fetchall():
            try:
                pk = row["jobid"]
                project_id = row["projectid"]

                if project_id not in self.project_map:
                    self.log(f"  Skipping job {pk}: project {project_id} not found")
                    continue

                project = self.project_map[project_id]
                job_number = row["jobnumber"]

                existing = Job.objects.filter(project=project, number=job_number).first()
                if existing:
                    self.job_map[pk] = existing
                    self.stats["jobs_skipped"] += 1
                    continue

                job = Job.objects.create(
                    uuid=convert_uuid(pk) or None,
                    project=project,
                    number=job_number,
                    title=row["jobtitle"] or "",
                    status=row["status"] or 0,
                    evaluation=row["evaluation"] or 0,
                    comments="",
                    creation_time=convert_timestamp(row["creationtime"]) or timezone.now(),
                    finish_time=convert_timestamp(row["finishtime"]),
                    task_name=row["taskname"] or "unknown",
                    process_id=row["processid"],
                )
                self.job_map[pk] = job
                self.stats["jobs"] += 1

                if row["parentjobid"]:
                    job_parent_map[job.id] = row["parentjobid"]

            except Exception as e:
                if self.continue_on_error:
                    self._record_error("Job", pk, e)
                else:
                    raise

        # Second pass: set parent references
        for job_id, legacy_parent_id in job_parent_map.items():
            if legacy_parent_id in self.job_map:
                Job.objects.filter(id=job_id).update(parent=self.job_map[legacy_parent_id])

        extra = f"{self.stats['jobs_skipped']} skipped" if self.stats.get("jobs_skipped") else ""
        self.log_count("Jobs", self.stats["jobs"], extra)

    def _import_serverjobs(self, cur):
        cur.execute(
            "SELECT JobId, ServerProcessId, Machine, Username, Mechanism, "
            "RemotePath, CustomCodeFile, Validate, KeyFilename, ServerGroup "
            "FROM ServerJobs"
        )
        for row in cur.fetchall():
            try:
                job_id = row["jobid"]
                if job_id not in self.job_map:
                    continue
                ServerJob.objects.create(
                    job=self.job_map[job_id],
                    server_process_id=row["serverprocessid"],
                    machine=row["machine"] or "",
                    username=row["username"] or "",
                    mechanism=row["mechanism"] or "",
                    remote_path=row["remotepath"] or "",
                    custom_code_file=row["customcodefile"] or "",
                    validate=row["validate"] or "",
                    key_file_name=row["keyfilename"] or "",
                    server_group=row["servergroup"] or "",
                )
                self.stats["serverjobs"] += 1
            except Exception as e:
                if self.continue_on_error:
                    self._record_error("ServerJob", job_id, e)
                else:
                    raise
        self.log_count("ServerJobs", self.stats["serverjobs"])

    def _import_jobkeyvalues(self, cur):
        cur.execute("SELECT JobID, KeyTypeID, Value FROM JobKeyValues")
        for row in cur.fetchall():
            try:
                job_id = row["jobid"]
                keytype_id = row["keytypeid"]
                if job_id not in self.job_map or keytype_id not in self.keytype_map:
                    continue
                job = self.job_map[job_id]
                key = self.keytype_map[keytype_id]
                if JobFloatValue.objects.filter(job=job, key=key).exists():
                    self.stats["jobfloatvalues_skipped"] += 1
                    continue
                JobFloatValue.objects.create(
                    job=job, key=key,
                    value=float(row["value"]) if row["value"] is not None else 0.0,
                )
                self.stats["jobfloatvalues"] += 1
            except Exception as e:
                if self.continue_on_error:
                    pass
                else:
                    raise
        extra = f"{self.stats['jobfloatvalues_skipped']} skipped" if self.stats.get("jobfloatvalues_skipped") else ""
        self.log_count("JobFloatValues", self.stats["jobfloatvalues"], extra)

    def _import_jobkeycharvalues(self, cur):
        cur.execute("SELECT JobID, KeyTypeID, Value FROM JobKeyCharValues")
        for row in cur.fetchall():
            try:
                job_id = row["jobid"]
                keytype_id = row["keytypeid"]
                if job_id not in self.job_map or keytype_id not in self.keytype_map:
                    continue
                job = self.job_map[job_id]
                key = self.keytype_map[keytype_id]
                if JobCharValue.objects.filter(job=job, key=key).exists():
                    self.stats["jobcharvalues_skipped"] += 1
                    continue
                JobCharValue.objects.create(
                    job=job, key=key, value=row["value"] or "",
                )
                self.stats["jobcharvalues"] += 1
            except Exception as e:
                if self.continue_on_error:
                    pass
                else:
                    raise
        extra = f"{self.stats['jobcharvalues_skipped']} skipped" if self.stats.get("jobcharvalues_skipped") else ""
        self.log_count("JobCharValues", self.stats["jobcharvalues"], extra)

    def _import_xdata(self, cur):
        cur.execute("SELECT XDataID, XDataClass, XDataXml, JobID FROM XData")
        for row in cur.fetchall():
            try:
                job = self.job_map.get(row["jobid"]) if row["jobid"] else None
                XData.objects.create(
                    data_class=row["xdataclass"] or "",
                    xml=row["xdataxml"] or "",
                    job=job,
                )
                self.stats["xdata"] += 1
            except Exception as e:
                if self.continue_on_error:
                    pass
                else:
                    raise
        self.log_count("XData", self.stats["xdata"])

    # ------------------------------------------------------------------
    # Files
    # ------------------------------------------------------------------

    def _import_files(self, cur):
        cur.execute(
            "SELECT FileID, Filename, Annotation, FiletypeID, FileSubType, "
            "FileContent, FilePath, JobID, JobParamName, PathFlag "
            "FROM Files"
        )
        for row in cur.fetchall():
            try:
                pk = row["fileid"]
                job_id = row["jobid"]
                if job_id not in self.job_map:
                    self.log(f"  Skipping file {pk}: job {job_id} not found")
                    continue

                file_uuid = convert_uuid(pk)
                if file_uuid and File.objects.filter(uuid=file_uuid).exists():
                    self.file_map[pk] = File.objects.get(uuid=file_uuid)
                    self.stats["files_skipped"] += 1
                    continue

                filetype_id = row["filetypeid"]
                if filetype_id not in self.filetype_map:
                    ft, _ = FileType.objects.get_or_create(
                        name=f"unknown_{filetype_id}",
                        defaults={"description": "Unknown file type from legacy import"},
                    )
                    self.filetype_map[filetype_id] = ft

                pathflag = row["pathflag"]
                directory = File.Directory.IMPORT_DIR if pathflag == 2 else File.Directory.JOB_DIR

                file_obj = File.objects.create(
                    uuid=file_uuid,
                    name=row["filename"] or "",
                    directory=directory,
                    type=self.filetype_map[filetype_id],
                    sub_type=row["filesubtype"],
                    content=row["filecontent"],
                    annotation=row["annotation"] or "",
                    job=self.job_map[job_id],
                    job_param_name=row["jobparamname"] or "",
                )
                self.file_map[pk] = file_obj
                self.stats["files"] += 1

            except Exception as e:
                if self.continue_on_error:
                    self._record_error("File", pk, e)
                else:
                    raise

        extra = f"{self.stats['files_skipped']} skipped" if self.stats.get("files_skipped") else ""
        self.log_count("Files", self.stats["files"], extra)

    def _import_importfiles(self, cur):
        cur.execute(
            "SELECT ImportId, FileID, SourceFilename, SourceFileID, "
            "ExportFileID, Annotation, CreationTime, LastModifiedTime, "
            "Checksum, ImportNumber, Reference "
            "FROM ImportFiles"
        )
        for row in cur.fetchall():
            try:
                file_id = row["fileid"]
                if file_id not in self.file_map:
                    continue
                file_obj = self.file_map[file_id]
                if FileImport.objects.filter(file=file_obj).exists():
                    self.stats["fileimports_skipped"] += 1
                    continue
                FileImport.objects.create(
                    file=file_obj,
                    time=convert_timestamp(row["creationtime"]) or timezone.now(),
                    name=row["sourcefilename"] or "",
                    checksum=row["checksum"] or "",
                    last_modified=convert_timestamp(row["lastmodifiedtime"]),
                )
                self.stats["fileimports"] += 1
            except Exception as e:
                if self.continue_on_error:
                    pass
                else:
                    raise
        extra = f"{self.stats['fileimports_skipped']} skipped" if self.stats.get("fileimports_skipped") else ""
        self.log_count("FileImports", self.stats["fileimports"], extra)

    def _import_exportfiles(self, cur):
        cur.execute(
            "SELECT ExportId, FileID, ExportFilename, Annotation, CreationTime "
            "FROM ExportFiles"
        )
        for row in cur.fetchall():
            try:
                file_id = row["fileid"]
                if file_id not in self.file_map:
                    continue
                FileExport.objects.create(
                    file=self.file_map[file_id],
                    time=convert_timestamp(row["creationtime"]) or timezone.now(),
                    name=row["exportfilename"] or "",
                )
                self.stats["fileexports"] += 1
            except Exception as e:
                if self.continue_on_error:
                    pass
                else:
                    raise
        self.log_count("FileExports", self.stats["fileexports"])

    def _import_fileuses(self, cur):
        cur.execute("SELECT FileID, JobID, RoleID, JobParamName FROM FileUses")
        for row in cur.fetchall():
            try:
                file_id = row["fileid"]
                job_id = row["jobid"]
                if file_id not in self.file_map or job_id not in self.job_map:
                    continue
                file_obj = self.file_map[file_id]
                job = self.job_map[job_id]
                role = self.filerole_map.get(row["roleid"], FileUse.Role.OUT)
                jpn = row["jobparamname"] or ""
                if FileUse.objects.filter(file=file_obj, job=job, role=role, job_param_name=jpn).exists():
                    self.stats["fileuses_skipped"] += 1
                    continue
                FileUse.objects.create(
                    file=file_obj, job=job, role=role, job_param_name=jpn,
                )
                self.stats["fileuses"] += 1
            except Exception as e:
                if self.continue_on_error:
                    pass
                else:
                    raise
        extra = f"{self.stats['fileuses_skipped']} skipped" if self.stats.get("fileuses_skipped") else ""
        self.log_count("FileUses", self.stats["fileuses"], extra)

    # ------------------------------------------------------------------
    # Project history
    # ------------------------------------------------------------------

    def _import_projectexports(self, cur):
        cur.execute(
            "SELECT ProjectExportId, ProjectID, ProjectExportTime "
            "FROM ProjectExports"
        )
        for row in cur.fetchall():
            try:
                pid = row["projectid"]
                if pid not in self.project_map:
                    continue
                ProjectExport.objects.create(
                    project=self.project_map[pid],
                    time=convert_timestamp(row["projectexporttime"]) or timezone.now(),
                )
                self.stats["projectexports"] += 1
            except Exception as e:
                if self.continue_on_error:
                    pass
                else:
                    raise
        self.log_count("ProjectExports", self.stats["projectexports"])

    def _import_projectimports(self, cur):
        cur.execute(
            "SELECT ProjectImportId, ProjectID, ProjectImportTime "
            "FROM ProjectImports"
        )
        for row in cur.fetchall():
            try:
                pid = row["projectid"]
                if pid not in self.project_map:
                    continue
                ProjectImport.objects.create(
                    project=self.project_map[pid],
                    time=convert_timestamp(row["projectimporttime"]) or timezone.now(),
                )
                self.stats["projectimports"] += 1
            except Exception as e:
                if self.continue_on_error:
                    pass
                else:
                    raise
        self.log_count("ProjectImports", self.stats["projectimports"])

    # ------------------------------------------------------------------
    # Comments (legacy table -> Job.comments or project description)
    # ------------------------------------------------------------------

    def _import_project_comments(self, cur):
        """Import ProjectComments into Project.description (concatenated)."""
        cur.execute(
            "SELECT ProjectCommentID, ProjectID, UserID, TimeOfComment, Comment "
            "FROM ProjectComments ORDER BY TimeOfComment"
        )
        comments_by_project = defaultdict(list)
        for row in cur.fetchall():
            pid = row["projectid"]
            if pid in self.project_map and row["comment"]:
                comments_by_project[pid].append(row["comment"])

        for pid, comments in comments_by_project.items():
            project = self.project_map[pid]
            combined = "\n".join(c for c in comments if c)
            if combined and not project.description:
                project.description = combined
                project.save(update_fields=["description"])
                self.stats["project_comments"] += 1

        self.log_count("Project comments", self.stats["project_comments"])

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------

    def summary(self):
        return {
            "dry_run": self.dry_run,
            "stats": dict(self.stats),
            "errors": self.errors,
            "source": str(self.db_path),
        }
