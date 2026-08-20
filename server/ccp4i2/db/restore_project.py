"""Rebuild database rows from the snapshots left in project directories.

There are two ways a project can arrive from outside, and they pull in opposite
directions:

**Import** brings a project in *alongside* an existing world. Another copy may
already be here, and its job numbers may collide with jobs that already exist.
Renumbering is then correct and necessary -- the incoming project has to be made
to fit, and ``import_ccp4_project_zip`` renames the job directories as it
extracts them so that disk and database still agree afterwards.

**Restore** is the opposite. The directory is the truth and the database is the
thing being repaired. Nothing may be renumbered, because the job directories are
already on disk under the names they have and nothing is going to move them.
Job numbers and UUIDs must come back exactly as they were, or every ``relPath``
in every ``params.xml``, and every job directory the database points at, is
wrong.

So a clash is an error to report here, not something to paper over, and the
post-condition is explicit: after a restore the database must agree with the
disk. :func:`verify_restored` checks that, rather than trusting that it held.

This module handles restore. Import already exists in :mod:`import_i2xml`.
"""

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional
from xml.etree import ElementTree as ET

from django.db import transaction

from . import project_snapshot
from .import_i2xml import import_i2xml
from .models import File, Job, Project
from .project_snapshot import SNAPSHOT_NAME, read_registry

logger = logging.getLogger(f"ccp4i2:{__name__}")


class RestoreError(Exception):
    """Raised when a project cannot be restored from its directory."""


@dataclass
class RestoreReport:
    """What a restore found, did, and could not do."""

    directory: str
    project_name: Optional[str] = None
    project_uuid: Optional[str] = None
    recorded_directory: Optional[str] = None
    jobs: int = 0
    files: int = 0
    restored: bool = False
    dry_run: bool = False
    skipped_reason: Optional[str] = None
    relocated: bool = False
    paths_rewritten: int = 0
    missing_job_directories: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)

    def as_dict(self) -> Dict:
        return {
            "directory": self.directory,
            "project_name": self.project_name,
            "project_uuid": self.project_uuid,
            "recorded_directory": self.recorded_directory,
            "jobs": self.jobs,
            "files": self.files,
            "restored": self.restored,
            "dry_run": self.dry_run,
            "skipped_reason": self.skipped_reason,
            "relocated": self.relocated,
            "paths_rewritten": self.paths_rewritten,
            "missing_job_directories": self.missing_job_directories,
            "warnings": self.warnings,
        }


# ---------------------------------------------------------------------------
# Reading what a directory offers
# ---------------------------------------------------------------------------


def snapshot_in(directory: Path) -> Optional[Path]:
    """The snapshot in ``directory``, if it has a readable one."""
    path = Path(directory) / SNAPSHOT_NAME
    return path if path.is_file() else None


def inspect(directory: Path) -> RestoreReport:
    """Read a project directory's snapshot and report what it holds.

    Purely read-only, so it can back a preview.
    """
    directory = Path(directory).expanduser()
    report = RestoreReport(directory=str(directory))

    path = snapshot_in(directory)
    if path is None:
        report.skipped_reason = f"no {SNAPSHOT_NAME} in this directory"
        return report

    try:
        root = ET.parse(path).getroot()
    except ET.ParseError as err:
        report.skipped_reason = f"{SNAPSHOT_NAME} is not readable: {err}"
        return report

    node = root.find("ccp4i2_body/projectTable/project")
    if node is None:
        report.skipped_reason = f"{SNAPSHOT_NAME} names no project"
        return report

    report.project_uuid = node.attrib.get("projectid")
    report.project_name = node.attrib.get("projectname")
    report.recorded_directory = node.attrib.get("projectdirectory")
    report.jobs = len(root.findall("ccp4i2_body/jobTable/job"))
    report.files = len(root.findall("ccp4i2_body/fileTable/file"))
    report.relocated = report.recorded_directory != str(directory)
    return report


# ---------------------------------------------------------------------------
# Restoring one project
# ---------------------------------------------------------------------------


def restore_from_directory(
    directory: Path,
    replace: bool = False,
    dry_run: bool = False,
    repair_paths: bool = True,
) -> RestoreReport:
    """Rebuild one project's rows from the snapshot in its directory.

    Args:
        directory: the project directory. Authoritative -- if the snapshot
            records a different location, the directory wins and the snapshot
            is re-rooted onto it.
        replace: rebuild a project that is already in the database, discarding
            the rows currently held for it. Without this an existing project
            is left alone, since overwriting a live project with an older
            snapshot would lose work rather than recover it.
        dry_run: report what would happen and change nothing.
        repair_paths: when the directory has moved since the snapshot was
            written, rewrite the absolute paths inside the project's files.
    """
    directory = Path(directory).expanduser().resolve()
    report = inspect(directory)
    report.dry_run = dry_run

    if report.skipped_reason is not None:
        return report

    existing = Project.objects.filter(uuid=report.project_uuid).first()
    by_name = Project.objects.filter(name=report.project_name).first()

    if existing is not None and not replace:
        report.skipped_reason = (
            f"project '{existing.name}' is already in the database; "
            "use replace to rebuild it from this directory"
        )
        return report
    if existing is None and by_name is not None:
        report.skipped_reason = (
            f"a different project is already called '{report.project_name}'"
        )
        return report

    clash = _job_number_clash(report.project_uuid, directory)
    if clash:
        report.skipped_reason = (
            f"job numbers {', '.join(clash[:5])} are already used by another "
            "project in this database; restoring would renumber them and the "
            "job directories on disk would no longer match"
        )
        return report

    if dry_run:
        report.missing_job_directories = _missing_job_directories_from_snapshot(
            directory
        )
        return report

    with project_snapshot.suspended():
        with transaction.atomic():
            if existing is not None:
                # Rows first, so the rebuilt numbering starts from a clean sheet.
                Job.objects.filter(project=existing).delete()
                existing.delete()
            _import_rerooted(directory)

    project = Project.objects.get(uuid=report.project_uuid)
    report.restored = True
    report.jobs = Job.objects.filter(project=project).count()
    report.files = File.objects.filter(job__project=project).count()

    if repair_paths and report.relocated and report.recorded_directory:
        report.paths_rewritten = _repair_after_relocation(
            project, report.recorded_directory
        )

    report.missing_job_directories = verify_restored(project)
    if report.missing_job_directories:
        report.warnings.append(
            f"{len(report.missing_job_directories)} job director"
            f"{'y is' if len(report.missing_job_directories) == 1 else 'ies are'} "
            "recorded but not present on disk"
        )

    # The project is in the database again; make sure its snapshot and the
    # registry reflect where it actually is now.
    project_snapshot.write_snapshot(project)
    project_snapshot.update_registry()

    logger.info(
        "Restored project '%s' from %s: %d job(s), %d file(s)",
        project.name,
        directory,
        report.jobs,
        report.files,
    )
    return report


def _import_rerooted(directory: Path) -> None:
    """Import the snapshot with ``directory`` substituted as the project root.

    ``import_i2xml``'s own ``relocate_path`` appends the recorded directory's
    basename to whatever it is given, which assumes the project folder kept its
    name. A restore knows exactly where the project is, so the path is set
    outright instead of being reconstructed.
    """
    root = ET.parse(directory / SNAPSHOT_NAME).getroot()
    for node in root.findall("ccp4i2_body/projectTable/project"):
        node.attrib["projectdirectory"] = str(directory)
    import_i2xml(root, relocate_path=None)


def _job_number_clash(project_uuid: str, directory: Path) -> List[str]:
    """Top-level job numbers this snapshot needs that another project holds.

    Only a clash *within the same project* matters, since job numbers are
    unique per project -- so this can only bite when the snapshot's project
    UUID differs from the row already occupying those numbers.
    """
    existing = Project.objects.filter(uuid=project_uuid).first()
    if existing is None:
        return []
    root = ET.parse(directory / SNAPSHOT_NAME).getroot()
    wanted = {
        node.attrib["jobnumber"]
        for node in root.findall("ccp4i2_body/jobTable/job")
        if "." not in node.attrib["jobnumber"]
    }
    held = set(
        Job.objects.filter(project=existing)
        .exclude(
            uuid__in=[
                node.attrib["jobid"]
                for node in root.findall("ccp4i2_body/jobTable/job")
            ]
        )
        .values_list("number", flat=True)
    )
    return sorted(wanted & held)


def _missing_job_directories_from_snapshot(directory: Path) -> List[str]:
    """Job numbers the snapshot records with no matching directory on disk."""
    root = ET.parse(directory / SNAPSHOT_NAME).getroot()
    missing = []
    for node in root.findall("ccp4i2_body/jobTable/job"):
        number = node.attrib["jobnumber"]
        parts = [f"job_{part}" for part in number.split(".")]
        if not (directory / "CCP4_JOBS").joinpath(*parts).is_dir():
            missing.append(number)
    return sorted(missing, key=lambda n: [int(p) for p in n.split(".")])


def verify_restored(project: Project) -> List[str]:
    """Job numbers whose directory the database claims but disk does not have.

    The post-condition of a restore is that the database agrees with the disk.
    Checking it directly catches renumbering, a mis-rooted directory and a
    snapshot that was written for a different tree, none of which would show up
    as an error during the import itself.
    """
    missing = []
    for job in Job.objects.filter(project=project):
        if not job.directory.is_dir():
            missing.append(job.number)
    return sorted(missing, key=lambda n: [int(p) for p in n.split(".")])


def _repair_after_relocation(project: Project, old_root: str) -> int:
    """Rewrite absolute paths inside a project restored from a new location."""
    from .move_project import apply_rebase, plan_rebase

    new_root = str(Path(project.directory))
    plan = plan_rebase(Path(new_root), old_root, new_root)
    return len(apply_rebase(plan))


# ---------------------------------------------------------------------------
# Restoring everything
# ---------------------------------------------------------------------------


def discover_restorable(root: Path) -> List[Path]:
    """Project directories under ``root`` that carry a snapshot.

    For when the registry is gone too, and all that is left is a folder full of
    project directories.
    """
    root = Path(root).expanduser()
    if not root.is_dir():
        return []
    found = []
    if snapshot_in(root) is not None:
        found.append(root)
    for child in sorted(root.iterdir()):
        if child.is_dir() and snapshot_in(child) is not None:
            found.append(child)
    return found


def restore_all(
    directories: List[Path],
    replace: bool = False,
    dry_run: bool = False,
) -> Dict:
    """Restore each directory in turn, reporting on all of them.

    One project failing does not abandon the rest: a partial recovery is worth
    having, and the ones that failed can be looked at individually.
    """
    reports = []
    for directory in directories:
        try:
            reports.append(
                restore_from_directory(directory, replace=replace, dry_run=dry_run)
            )
        except Exception as err:
            logger.exception("Could not restore %s", directory)
            failed = RestoreReport(directory=str(directory))
            failed.skipped_reason = str(err)
            reports.append(failed)

    return {
        "dry_run": dry_run,
        "restored": [r.as_dict() for r in reports if r.restored],
        "skipped": [r.as_dict() for r in reports if not r.restored],
        "total": len(reports),
    }


def restore_from_registry(replace: bool = False, dry_run: bool = False) -> Dict:
    """Restore every project the directory registry knows about.

    The registry is kept outside the database precisely so that it survives the
    database, and this is what it is for.
    """
    entries = read_registry()
    directories = [Path(entry["directory"]) for entry in entries if entry.get("directory")]
    result = restore_all(directories, replace=replace, dry_run=dry_run)
    result["registry"] = str(project_snapshot.registry_path())
    result["registry_entries"] = len(entries)
    return result
