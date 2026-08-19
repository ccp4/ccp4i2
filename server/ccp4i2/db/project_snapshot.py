"""Keep an on-disk record of each project, so the database can be rebuilt.

Everything CCP4i2 knows about a project's *files* is recoverable from the
project directory: job directories, ``params.xml``, the files themselves. What
is not recoverable is everything the user typed. Annotations, job titles and
comments, evaluations of good and rejected, tags, the project description --
none of that has an on-disk artefact. It lives only in ``db.sqlite3``, and if
that is lost or corrupted it is gone.

The Qt-era CCP4i2 wrote a per-project ``DATABASE.db.xml`` as jobs started and
finished, and kept a list of known project directories outside the database, so
a lost database could be reconstructed. The Django app has had neither. This
module restores both.

The snapshot is the same format the project export already produces and the
project import already reads -- the interchange format that project zips carry
and that legacy project directories on disk already contain. Reusing it means
one recovery path rather than two, and it keys on UUIDs rather than on primary
keys or on the current shape of the models, which is what a file written to be
read back after something went wrong needs.

Writing is driven by signals rather than by calls at each endpoint, so it
catches every write path -- REST, i2run, management commands, the job runner --
including ones that do not exist yet.

**On coalescing.** A single request can touch many rows: saving a task writes a
dozen files, a finishing job writes its outputs. Snapshots are therefore
registered with ``transaction.on_commit`` and de-duplicated per transaction, so
one request produces one write however much it changed. That is coalescing by
transaction, not by wall-clock: there is deliberately no timer anywhere here.
The desktop app runs uvicorn with two workers and the CLI runs in processes that
exit immediately, so there is no long-lived process for a timer to live in, and
an in-process one would be either duplicated or never fire.
"""

import json
import logging
import os
import shutil
import tempfile
import threading
from contextlib import contextmanager
from pathlib import Path
from typing import Dict, List, Optional
from xml.etree import ElementTree as ET

from django.db import transaction

from .export_project import generate_project_xml_tree
from .models import Project

logger = logging.getLogger(f"ccp4i2:{__name__}")


#: Written at the root of each project directory. The name the export and the
#: import already agree on, and the one legacy project directories carry.
SNAPSHOT_NAME = "DATABASE.db.xml"

#: Where a snapshot we did not write is set aside before being replaced.
PRESERVED_NAME = "DATABASE.db.xml.pre-ccp4i2-django"

#: Stamped into the header, so a later run can tell our snapshots from one a
#: Qt-era CCP4i2 left in the project directory.
GENERATOR_MARK = "ccp4i2-django"

#: The element as it appears in the file. Recognition matches this rather than
#: the marker alone, which occurs in ordinary directory paths.
GENERATOR_ELEMENT = f"<generator>{GENERATOR_MARK}</generator>"

#: The header is the first thing in the file; no need to read further.
GENERATOR_SEARCH_BYTES = 8192

#: A list of every project directory known to this installation, kept outside
#: the database so that a lost database can be told where to look.
REGISTRY_NAME = "project_directories.json"

#: Set to disable snapshot writing entirely -- for bulk imports, migrations and
#: tests that would otherwise snapshot after every row.
DISABLE_ENV = "CCP4I2_DISABLE_PROJECT_SNAPSHOT"


_state = threading.local()


@contextmanager
def suspended():
    """Suspend snapshot writing for the duration of the block.

    For bulk work -- importing a project zip, migrating a legacy database --
    where the caller will snapshot once at the end rather than after every row.
    """
    previous = getattr(_state, "suspended", False)
    _state.suspended = True
    try:
        yield
    finally:
        _state.suspended = previous


def is_enabled() -> bool:
    if getattr(_state, "suspended", False):
        return False
    return os.environ.get(DISABLE_ENV, "").lower() not in ("1", "true", "yes")


# ---------------------------------------------------------------------------
# Writing
# ---------------------------------------------------------------------------


def snapshot_path(project: Project) -> Path:
    return Path(project.directory) / SNAPSHOT_NAME


def write_snapshot(project: Project) -> Optional[Path]:
    """Write ``project``'s snapshot, atomically. Returns the path, or None.

    Never raises. A snapshot is a safety net: failing to write one must not
    fail the operation that triggered it, and a job that has just finished
    should not be reported as failed because its directory is read-only.
    """
    directory = Path(project.directory)
    if not directory.is_dir():
        logger.debug(
            "Not snapshotting project '%s': %s does not exist", project.name, directory
        )
        return None

    _preserve_foreign_snapshot(directory)

    try:
        root = generate_project_xml_tree(project)
        header = root.find("ccp4i2_header")
        if header is not None:
            ET.SubElement(header, "generator").text = GENERATOR_MARK
        ET.indent(root, space="  ")
        payload = ET.tostring(root, encoding="utf-8", xml_declaration=True)
        path = _atomic_write_bytes(directory, SNAPSHOT_NAME, payload)
        logger.debug("Wrote project snapshot %s (%d bytes)", path, len(payload))
        return path
    except Exception as err:
        logger.warning(
            "Could not write snapshot for project '%s': %s", project.name, err
        )
        return None


def _preserve_foreign_snapshot(directory: Path) -> None:
    """Keep a one-time copy of a snapshot this installation did not write.

    A project directory carried over from Qt-era CCP4i2 already has a
    ``DATABASE.db.xml``, and it may describe jobs or annotations that were
    never imported into this database. Overwriting it would destroy the only
    remaining record of them. So the first time we are about to write over a
    snapshot we did not author, the original is set aside.

    One-time by construction: once the backup exists it is never touched again,
    so this cannot creep over the real backup with our own output.
    """
    existing = directory / SNAPSHOT_NAME
    backup = directory / PRESERVED_NAME
    if not existing.is_file() or backup.exists():
        return
    try:
        if _was_written_here(existing):
            return
        shutil.copy2(existing, backup)
        logger.info("Preserved the pre-existing %s as %s", existing, backup.name)
    except OSError as err:
        logger.warning("Could not preserve %s: %s", existing, err)


def _was_written_here(path: Path) -> bool:
    """True if this module wrote ``path`` (it stamps a generator element).

    Matches the whole element, not the marker on its own. A bare substring
    search finds the marker inside any path that happens to contain it -- and
    real project directories do: the first project this was tried on referenced
    ``/Users/.../ccp4i2-djangodbapi/demo_data``, and was consequently mistaken
    for one of ours and would have been overwritten without being preserved.
    """
    try:
        with open(path, "rb") as stream:
            head = stream.read(GENERATOR_SEARCH_BYTES)
    except OSError:
        return False
    return GENERATOR_ELEMENT.encode("utf-8") in head


def _atomic_write_bytes(directory: Path, name: str, payload: bytes) -> Path:
    """Write via a temp file and os.replace, so a reader never sees a partial.

    The desktop app runs two uvicorn workers, so two processes can write the
    same snapshot at once. Both read the same database, so their content agrees
    and last-writer-wins is harmless -- but a half-written file would not be.
    """
    target = directory / name
    handle, temp_name = tempfile.mkstemp(dir=str(directory), prefix=f".{name}-")
    try:
        with os.fdopen(handle, "wb") as stream:
            stream.write(payload)
        os.replace(temp_name, target)
    finally:
        if os.path.exists(temp_name):
            os.unlink(temp_name)
    return target


def _already_scheduled(project_pk) -> bool:
    """True when this transaction has a snapshot pending for ``project_pk``.

    The pending set is Django's own ``run_on_commit`` list rather than one kept
    alongside it. That matters: Django clears that list on rollback as well as
    on commit, so the de-duplication state lives exactly as long as the
    transaction does. A separate set would keep an entry for a transaction that
    rolled back, and the next change to that project would be silently skipped
    -- a safety net with a hole in it that nothing would ever report.

    In autocommit there is no pending list, so every call schedules its own
    write, which is correct: each save is its own transaction.
    """
    for entry in transaction.get_connection().run_on_commit:
        callback = entry[1]
        if getattr(callback, "snapshot_project_pk", None) == project_pk:
            return True
    return False


def schedule_snapshot(project: Optional[Project]) -> None:
    """Arrange for ``project`` to be snapshotted when this transaction commits.

    Repeated calls within one transaction collapse to a single write, so a
    request that touches thirty files still writes once.
    """
    if project is None or not is_enabled():
        return

    project_pk = project.pk
    if _already_scheduled(project_pk):
        return

    def run():
        # Re-read: the instance that triggered this may be stale by commit time,
        # and may have been deleted altogether.
        current = Project.objects.filter(pk=project_pk).first()
        if current is not None:
            write_snapshot(current)
        update_registry()

    run.snapshot_project_pk = project_pk
    transaction.on_commit(run)


# ---------------------------------------------------------------------------
# The registry of known project directories
# ---------------------------------------------------------------------------


def registry_path() -> Path:
    from ..config.preferences import ccp4i2_home

    return ccp4i2_home() / REGISTRY_NAME


def update_registry() -> Optional[Path]:
    """Rewrite the list of known project directories. Never raises.

    A snapshot inside a project directory is only useful if something can find
    the project directory in the first place. With the database gone, this file
    is the only thing that knows where the projects are.
    """
    if not is_enabled():
        return None
    try:
        entries = [
            {
                "uuid": str(project.uuid),
                "name": project.name,
                "directory": str(project.directory),
            }
            for project in Project.objects.all().order_by("name")
        ]
        home = registry_path().parent
        home.mkdir(parents=True, exist_ok=True)
        payload = json.dumps({"projects": entries}, indent=2).encode("utf-8")
        return _atomic_write_bytes(home, REGISTRY_NAME, payload)
    except Exception as err:
        logger.warning("Could not update the project directory registry: %s", err)
        return None


def read_registry() -> List[Dict[str, str]]:
    """The known project directories, or an empty list if there is no registry."""
    try:
        with open(registry_path(), encoding="utf-8") as stream:
            data = json.load(stream)
    except (OSError, json.JSONDecodeError):
        return []
    projects = data.get("projects") if isinstance(data, dict) else None
    return projects if isinstance(projects, list) else []


# ---------------------------------------------------------------------------
# Rebuilding everything at once
# ---------------------------------------------------------------------------


def has_snapshot(project: Project) -> bool:
    """True when this installation has written a snapshot for ``project``.

    A Qt-era ``DATABASE.db.xml`` sitting in the directory does not count: it
    describes whatever that installation knew, which is not necessarily what
    this database holds.
    """
    path = snapshot_path(project)
    return path.is_file() and _was_written_here(path)


def snapshot_status() -> List[Dict]:
    """Per project, whether it is protected and when it was last written."""
    status = []
    for project in Project.objects.all().order_by("name"):
        path = snapshot_path(project)
        exists = path.is_file()
        status.append(
            {
                "id": project.pk,
                "name": project.name,
                "directory": str(project.directory),
                "directory_exists": Path(project.directory).is_dir(),
                "has_snapshot": exists and _was_written_here(path),
                "foreign_snapshot": exists and not _was_written_here(path),
                "written": (
                    int(path.stat().st_mtime) if exists else None
                ),
            }
        )
    return status


def snapshot_all_projects(missing_only: bool = False) -> Dict[str, int]:
    """Snapshot every project, and refresh the registry.

    Projects that predate this mechanism never trigger a signal -- a finished,
    dormant project is exactly the kind you most want protected and exactly the
    kind that nothing will touch again -- so there has to be a way to write
    them all out on demand.

    Args:
        missing_only: skip projects this installation has already snapshotted.
            Cheap enough to run on every start-up; use False to refresh all.
    """
    written, skipped, already = 0, 0, 0
    for project in Project.objects.all().order_by("name"):
        if missing_only and has_snapshot(project):
            already += 1
            continue
        if write_snapshot(project) is not None:
            written += 1
        else:
            skipped += 1
    update_registry()
    logger.info(
        "Snapshotted %d project(s); %d already current, %d could not be written",
        written,
        already,
        skipped,
    )
    return {"written": written, "skipped": skipped, "already_current": already}
